#include "MAPL_Generic.h"

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  SU_GridCompMod --- SU Grid Component Class
!
! !INTERFACE:
!

   module  SU_GridCompMod

! !USES:

   use ESMF
   use MAPL

   use Chem_Mod              ! Chemistry Base Class
   use Chem_StateMod         ! Chemistry State
   use Chem_ConstMod, only: grav, von_karman, cpd, &   ! Constants !
                            undefval => undef, airMolWght => airmw
   use Chem_UtilMod          ! I/O
   use Chem_MieMod           ! Aerosol LU Tables, calculator
   use m_inpak90             ! Resource file management
   use m_die, only: die
   USE Henrys_law_ConstantsMod, ONLY: get_HenrysLawCts

   use m_StrTemplate
   use SulfateChemDriverMod
   use ConvectionMod         ! Offline convective mixing/scavenging
   use Chem_SettlingMod      ! Gravitiational Settling

   implicit none

! !PUBLIC TYPES:
!
   PRIVATE
   PUBLIC  SU_GridComp       ! The SU object 
   PUBLIC  SU_GridComp1      ! Single instance SU object

!
! !PUBLIC MEMBER FUNCTIONS:
!

   PUBLIC  SU_GridCompSetServices
   PUBLIC  SU_GridCompInitialize
   PUBLIC  SU_GridCompRun1
   PUBLIC  SU_GridCompRun2
   PUBLIC  SU_GridCompFinalize

!
! !DESCRIPTION:
!
!  This module implements the (pre-ESMF) SU Grid Component. 
!
! !REVISION HISTORY:
!
!  16Sep2003 da Silva  First crack.
!  18May2006 da Silva  Removed ensure postive, now in GOCART_GridComp
!  25Aug2009 Nielsen   Connections, usage of GMI Combo OH, H2O2, and NO3
!
!EOP
!-------------------------------------------------------------------------

! Note that the dates associated with the input files are a real mess
! Chem_UtilMPread cares about the date!
! Arbitrarily I set so2biomass_src, so2anthro_l1_src, and so2anthro_l2_src to 1971
! (of course these are not really valid for 1971)
! DMSO is valid 2000
! OH, NO3, H2O2 files are valid 2001
! Go figure...this is what happens when I get inputs from other people
! who are not the primary sources (e.g., Mian and Bian instead of
! geoschem...I get what they've got).

  type SU_GridComp1
        character(len=255) :: name
        character(len=255) :: iname           ! instance name
        character(len=255) :: rcfilen         ! resource file name
        character(len=255) :: maskFileName
        character(len=255) :: regionsString   ! Comma-delimited string of regions
        real, pointer      :: regionMask(:,:) ! regional mask
        integer :: instance                   ! instance number
        logical :: run_alarm = .false.        ! run alarm

        type(Chem_Mie), pointer :: mie_tables => null()  ! aod LUTs
        real, pointer :: so2biomass_src_(:,:) ! before diurnal
        real, pointer :: so2biomass_src(:,:)
        real, pointer :: so2anthro_l1_src(:,:)  ! level 1
        real, pointer :: so2anthro_l2_src(:,:)  ! level 2
        real, pointer :: so2ship_src(:,:)
        real, pointer :: so4ship_src(:,:)
        real, pointer :: aircraft_fuel_src(:,:,:)
        real, pointer :: aviation_lto_src(:,:)  ! aviation - landing and takeoff
        real, pointer :: aviation_cds_src(:,:)  ! aviation - climbing and descent
        real, pointer :: aviation_crs_src(:,:)  ! aviation - cruise
        real          :: aviation_layers(4)     ! heights of the LTO, CDS and CRS layers
        real, pointer :: dmso_conc(:,:)
!       Special handling for volcanic emissions
        integer :: nvolc = 0
        real, pointer, dimension(:)    :: vLat   => null(), &
                                          vLon   => null(), &
                                          vSO2   => null(), &
                                          vElev  => null(), &
                                          vCloud => null()
        integer, pointer, dimension(:) :: vStart => null(), &
                                          vEnd   => null()
!       Note that the OH, NO3, and H2O2 are from a geoschem run
!       Ideally would be from a run of the fv chemistry package!
        real, pointer :: oh_conc(:,:,:)
        real, pointer :: no3_mr(:,:,:)
        real, pointer :: h2o2_mr(:,:,:)
!       OH and NO3 are scaled every timestep.  H2O2 is replaced every
!       3 hours with the monthly value.  Hence, we need to save a value
!       somewhere!  For now we save the instantaneous value here.
        real, pointer :: h2o2_int(:,:,:)
        real :: fSO4ant         ! Fraction of anthropogenic emissions are SO4
        real :: eAircraftFuel   ! Emission factor to go from fuel to SO2
        real :: fMassSulfur     ! gram molar weight of S
        real :: fMassSO2        ! gram molar weight of SO2
        real :: fMassSO4        ! gram molar weight of SO4
        real :: fMassDMS        ! gram molar weight of DMS
        real :: fMassMSA        ! gram molar weight of MSA
        integer :: nDMS
        integer :: nSO2
        integer :: nSO4
        integer :: nMSA
        integer :: nymd          ! Update the emissions?
        integer :: nymd_oxidants ! Update the oxidant files?
        logical :: using_GMI_OH
        logical :: using_GMI_NO3
        logical :: using_GMI_H2O2
        logical :: using_ACHEM_pSO2_OCS
        logical :: export_H2O2
        logical :: firstRun
        logical :: recycle_H2O2 = .false.
        character(len=255) :: volcano_srcfilen
!       parameters for sulfate gravitational settling
        integer :: rhFlag          !flag for sulfate growth parameterization
        real, pointer :: radius(:) !particle effective radius [um]
        real, pointer :: rhop(:)   ! SU class density [kg m-3]
! age of tracers
        integer :: myDOW = -1             ! my Day of the week: Sun=1, Mon=2,...,Sat=7
        logical :: doing_nei=.FALSE.      ! NEI08: National Emission Inventory (US+Canada)
        real    :: nei_lon(2), nei_lat(2) ! NEI bounding box; superseeds eocant1/2 inside
        character(len=255) :: nei_srcfilen(2) ! 1=bottom layer, 2=above bottom layer
        integer :: nei_hour = -1
        integer :: nei_year = 2010        ! Hardwire this for now

!       Workspace for any requested point emissions (handled in run)
!       ------------------------------------------------------------
        logical :: doing_point_emissions=.FALSE.  ! Providing pointwise emissions
        character(len=255) :: point_emissions_srcfilen   ! filename for pointwise emissions
        integer                         :: nPts = -1
        integer, pointer, dimension(:)  :: pstart => null(), pend => null()
        real, pointer, dimension(:)     :: pLat  => null(), &
                                           pLon  => null(), &
                                           pBase => null(), &
                                           pTop  => null(), &
                                           pEmis => null()
  end type SU_GridComp1

  type SU_GridComp
     integer                     :: n = 0                ! number of instances 
     type(Chem_Mie), pointer     :: mie_tables => null() ! aod LUTs
     type(SU_GridComp1), pointer :: gcs(:)     => null() ! instances
  end type SU_GridComp

  character(len=*), parameter :: rc_basename = 'SU_GridComp'

  real, parameter :: OCEAN=0.0, LAND = 1.0, SEA_ICE = 2.0
  real, parameter :: pi = 3.1415, rearth = 6.37e6
  real, parameter :: radTODeg = 57.2957795
  real, parameter :: rH2O2 = 34./airMolWght

CONTAINS

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  SU_GridCompSetServices ---  SetServices SU_GridComp
!
! !INTERFACE:
!

   subroutine SU_GridCompSetServices(  gc, chemReg, rc)
   type(ESMF_GridComp), intent(INOUT) :: GC
   type(Chem_Registry), intent(INOUT) :: chemReg
   integer,             intent(OUT  ) :: rc

   CHARACTER(LEN=255) :: name

   integer            :: status
   character(len=ESMF_MAXSTR) :: Iam

   integer :: n, i

   type(ESMF_Config) :: cfg

   Iam = "SU_GridCompSetServices"

!  Load resource file
!  ------------------
   cfg = ESMF_ConfigCreate(__RC__)
   call ESMF_ConfigLoadFile(cfg, trim(rc_basename)//'.rc', __RC__)

!  Parse resource file
!  -------------------
   n = ESMF_ConfigGetLen(cfg, label='SU_instances:', __RC__)

!  We have 4 tracers for each instance of SU
!  Chem_Registry provides the number (total)
!  of tracers to be run.  Therefore n*4 must
!  be >= to that number or else we don't have
!  enough instances requested.
!  --------------------------------------------------------
   if ( n*4 .lt. chemReg%n_SU ) then
        rc = 35
        return
   end if
   n = min(n,chemReg%n_SU/4 )

   call ESMF_ConfigFindLabel(cfg, 'SU_instances:', __RC__)

   do i = 1, n
      ! read the resource file name
      call ESMF_ConfigGetAttribute(cfg, name, __RC__)

      if (trim(name) == "full" ) then
       name = " "              ! blank instance name for full (1)
      else
       name = trim(name)       ! instance name for others
      END IF

      call SU_GridCompSetServices1_(gc, chemReg, name, __RC__)
   end do

   call MAPL_AddImportSpec(GC,           &
        SHORT_NAME = 'SU_regionMask',    &
        LONG_NAME  = 'source species'  , &
        UNITS      = '1',                &
        DIMS       = MAPL_DimsHorzOnly,  &
        VLOCATION  = MAPL_VLocationNone, &
        RESTART    = MAPL_RestartSkip,   &
        RC         = STATUS)
   VERIFY_(STATUS)

   call MAPL_AddImportSpec(GC,             &
        SHORT_NAME = 'pSO2_OCS',           &
        LONG_NAME  = 'source species'  ,   &
        UNITS      = '1',                  &
        DIMS       = MAPL_DimsHorzVert,    &
        VLOCATION  = MAPL_VLocationCenter, &
        RESTART    = MAPL_RestartSkip,     &
        RC         = STATUS)
   VERIFY_(STATUS)

   RETURN_(ESMF_SUCCESS)
   end subroutine SU_GridCompSetServices

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  SU_GridCompInitialize --- Initialize SU_GridComp
!
! !INTERFACE:
!

   subroutine SU_GridCompInitialize ( gcSU, w_c, impChem, expChem, &
                                      nymd, nhms, cdt, rc )

! !USES:

  implicit NONE

! !INPUT PARAMETERS:

   type(Chem_Bundle), intent(inout) :: w_c        ! Chemical tracer fields      
   integer, intent(in) :: nymd, nhms           ! time
   real, intent(in) :: cdt                     ! chemistry timestep (secs)

! !OUTPUT PARAMETERS:

   type(SU_GridComp), intent(inout) :: gcSU   ! Grid Component
   type(ESMF_State), intent(inout)  :: impChem  ! Import State
   type(ESMF_State), intent(inout)  :: expChem  ! Export State
   integer, intent(out) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 - 

! !DESCRIPTION: Initializes the SU Grid Component. It primarily sets
!               the import state for each active constituent package.
!
! !REVISION HISTORY:
!
!  18Sep2003 da Silva  First crack.
!
!EOP
!-------------------------------------------------------------------------

   character(len=*), parameter :: myname = 'SU_GridCompInitialize'
   CHARACTER(LEN=255) :: name
   
   integer :: i, ier, n, i_

!  Load resource file
!  ------------------
   call i90_loadf ( trim(rc_basename)//'.rc', ier )
   if ( ier .NE. 0 ) then
      rc = 10
      return
   end if

!  Parse resource file
!  -------------------
   CALL I90_label ( 'SU_instances:', ier )
   if ( ier .NE. 0 ) then
      rc = 20
      return
   end if

!  First determine how many instances we have
!  ------------------------------------------   
   n = 0
   do while ( ier .EQ. 0 )
      CALL I90_gtoken( name, ier )
      if(ier .eq. 0) n = n + 1
   end do
   if ( n .EQ. 0 ) then
      rc = 30
      return
   end if
   
!  We have 4 tracers for each instance of SU
!  Chem_Registry provides the number (total)
!  of tracers to be run.  Therefore n*4 must
!  be >= to that number or else we don't have
!  enough instances requested.
!  --------------------------------------------------------
   if ( n*4 .lt. w_c%reg%n_SU ) then
        rc = 35
        return
   end if
   n = min(n,w_c%reg%n_SU/4 )
   gcSU%n = n

!  Next allocate necessary memory
!  ------------------------------
   allocate ( gcSU%gcs(n), stat=ier )    
   if ( ier .NE. 0 ) then
      rc = 40
      return
   end if

!  Record name of each instance
!  ----------------------------
   CALL I90_label ( 'SU_instances:', ier )
   do i = 1, n
      CALL I90_gtoken( name, ier )
      if ( ier .NE. 0 ) then
         rc = 40
         return
      end if
                                            ! resource file name
      gcSU%gcs(i)%rcfilen = trim(rc_basename)//'---'//trim(name)//'.rc'
      gcSU%gcs(i)%instance = i              ! instance number

      if (trim(name) == "full") then
         gcSU%gcs(i)%iname = " "            ! blank instance name for full (1)
      else
         gcSU%gcs(i)%iname = trim(name)     ! instance name for others
      end if
      
   end do    

!  Next initialize each instance
!  -----------------------------
   do i = 1, gcSU%n
      IF(MAPL_AM_I_ROOT()) THEN
       PRINT *," "
       PRINT *,myname,": Initializing instance ",TRIM(gcSU%gcs(i)%iname)," [",gcSU%gcs(i)%instance,"]"
      END IF
      call SU_SingleInstance_ ( SU_GridCompInitialize1_, i, &
                                gcSU%gcs(i), w_c, impChem, expChem,  &
                                nymd, nhms, cdt, ier )
      if ( ier .NE. 0 ) then
         rc = 1000+ier
         return
      end if
      gcSU%gcs(i)%mie_tables => gcSU%mie_tables
   end do

!  Get Henrys Law cts for the parameterized convective wet removal
!  -----------------------------------------------------------
   do i = 1, gcSU%n
      !- DMS     
      i_ = w_c%reg%i_SU  + 4*(i - 1)
      CALL get_HenrysLawCts('DMS',w_c%reg%Hcts(1,i_),w_c%reg%Hcts(2,i_)&
                                 ,w_c%reg%Hcts(3,i_),w_c%reg%Hcts(4,i_))  
      !print*,"DMS=",w_c%reg%Hcts(1,i_),w_c%reg%Hcts(2,i_),w_c%reg%Hcts(3,i_),w_c%reg%Hcts(4,i_)

      !- SO2     
      i_ =  w_c%reg%i_SU  + 4*(i - 1) + 1
      CALL get_HenrysLawCts('SO2',w_c%reg%Hcts(1,i_),w_c%reg%Hcts(2,i_)&
				 ,w_c%reg%Hcts(3,i_),w_c%reg%Hcts(4,i_))  
      !print*,"SO2=",w_c%reg%Hcts(1,i_),w_c%reg%Hcts(2,i_),w_c%reg%Hcts(3,i_),w_c%reg%Hcts(4,i_)
      !call flush(6)
   ENDDO


!  All done
!  --------
   CALL I90_FullRelease( ier )
   IF( ier /= 0 ) THEN
    PRINT *,myname,": I90_FullRelease not successful."
    rc = 40
   END IF


 end subroutine SU_GridCompInitialize

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  SU_GridCompRun1 --- Run SU_GridComp
!
! !INTERFACE:
!

   subroutine SU_GridCompRun1 ( gcSU, w_c, impChem, expChem, &
                                      nymd, nhms, cdt, rc )

! !USES:

  IMPLICIT NONE

! !INPUT PARAMETERS:

   TYPE(Chem_Bundle), intent(in) :: w_c        ! Chemical tracer fields      
   INTEGER, INTENT(IN) :: nymd, nhms           ! time
   REAL,    INTENT(IN) :: cdt                  ! chemical timestep (secs)


! !OUTPUT PARAMETERS:

   TYPE(SU_GridComp), INTENT(INOUT) :: gcSU     ! Grid Component
   TYPE(ESMF_State), INTENT(INOUT)  :: impChem  ! Import State
   TYPE(ESMF_State), INTENT(INOUT)  :: expChem  ! Export State
   INTEGER, INTENT(OUT) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 - 

! !DESCRIPTION: Runs the CO Grid Component. Multiple instance
!               version.
!
! !REVISION HISTORY:
!
!  27Feb2008  da Silva  Introduced multiple instances
!
!EOP
!-------------------------------------------------------------------------

   integer :: i, ier

   do i = 1, gcSU%n
      call SU_SingleInstance_ ( SU_GridCompRun1_, i, &
                                gcSU%gcs(i), w_c, impChem, expChem, &
                                nymd, nhms, cdt, ier )
      if ( ier .NE. 0 ) then
         rc = i * 1000+ier
         return
      end if
   end do

 end subroutine SU_GridCompRun1


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  SU_GridCompRun2 --- Run SU_GridComp
!
! !INTERFACE:
!

   subroutine SU_GridCompRun2 ( gcSU, w_c, impChem, expChem, &
                                run_alarm, nymd, nhms, cdt, rc )

! !USES:

  IMPLICIT NONE

! !INPUT PARAMETERS:

   TYPE(Chem_Bundle), intent(in) :: w_c        ! Chemical tracer fields      
   LOGICAL, INTENT(IN) :: run_alarm            ! run alarm
   INTEGER, INTENT(IN) :: nymd, nhms           ! time
   REAL,    INTENT(IN) :: cdt                  ! chemical timestep (secs)


! !OUTPUT PARAMETERS:

   TYPE(SU_GridComp), INTENT(INOUT) :: gcSU     ! Grid Component
   TYPE(ESMF_State), INTENT(INOUT)  :: impChem  ! Import State
   TYPE(ESMF_State), INTENT(INOUT)  :: expChem  ! Export State
   INTEGER, INTENT(OUT) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 - 

! !DESCRIPTION: Runs the CO Grid Component. Multiple instance
!               version.
!
! !REVISION HISTORY:
!
!  27Feb2008  da Silva  Introduced multiple instances
!
!EOP
!-------------------------------------------------------------------------

   integer :: i, ier

   do i = 1, gcSU%n
      gcSU%gcs(i)%run_alarm = run_alarm

      call SU_SingleInstance_ ( SU_GridCompRun2_, i, &
                                gcSU%gcs(i), w_c, impChem, expChem, &
                                nymd, nhms, cdt, ier )
      if ( ier .NE. 0 ) then
         rc = i * 1000+ier
         return
      end if
   end do

 end subroutine SU_GridCompRun2



!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  SU_GridCompFinalize --- Initialize SU_GridComp
!
! !INTERFACE:
!

   subroutine SU_GridCompFinalize ( gcSU, w_c, impChem, expChem, &
                                      nymd, nhms, cdt, rc )

! !USES:

  IMPLICIT NONE

! !INPUT PARAMETERS:

   TYPE(Chem_Bundle), intent(in) :: w_c        ! Chemical tracer fields      
   INTEGER, INTENT(IN) :: nymd, nhms           ! time
   REAL,    INTENT(IN) :: cdt                  ! chemical timestep (secs)


! !OUTPUT PARAMETERS:

   TYPE(SU_GridComp), INTENT(INOUT) :: gcSU     ! Grid Component
   TYPE(ESMF_State), INTENT(INOUT)  :: impChem  ! Import State
   TYPE(ESMF_State), INTENT(INOUT)  :: expChem  ! Export State
   INTEGER, INTENT(OUT) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 - 

! !DESCRIPTION: Finalizes the SU Grid Component. Multiple instance
!               version.
!
! !REVISION HISTORY:
!
!  27Feb2008  da Silva  Introduced multiple instances
!
!EOP
!-------------------------------------------------------------------------

   integer i, ier

   do i = 1, gcSU%n
      call SU_SingleInstance_ ( SU_GridCompFinalize1_, i, &
                                gcSU%gcs(i), w_c, impChem, expChem, &
                                nymd, nhms, cdt, ier )
      if ( ier .NE. 0 ) then
         rc = i * 1000+ier
         return
      end if
   end do

   if (associated(gcSU%gcs)) deallocate ( gcSU%gcs, stat=ier )
   gcSU%n = -1

 end subroutine SU_GridCompFinalize


 subroutine SU_GridCompSetServices1_(  gc, chemReg, iname, rc)
 type(ESMF_GridComp), intent(INOUT) :: GC
 type(Chem_Registry), intent(INOUT) :: chemReg
 character(len=*),    intent(IN   ) :: iname
 integer,             intent(OUT  ) :: rc


 CHARACTER(LEN=255) :: name
 type(ESMF_Config)  :: cfg

 logical :: using_GMI_H2O2, using_GMI_OH, using_GMI_NO3
 logical :: using_ACHEM_pSO2_OCS, doing_NEI
 character(len=ESMF_MAXSTR) :: tmpchar

 integer :: Status
 character(len=ESMF_MAXSTR) :: Iam

 Iam ="SU_GridCOmpSetServices1_"

 if (trim(iname) == "") then
    name = trim(rc_basename)//'---full.rc'
 else
    name = trim(rc_basename)//'---'//trim(iname)//'.rc'
 end if
 
 cfg = ESMF_ConfigCreate(__RC__)
 call ESMF_ConfigLoadFile(cfg, name, __RC__)

 call ESMF_ConfigGetAttribute(cfg, using_GMI_H2O2, label='using_GMI_H2O2:', __RC__)
 call ESMF_ConfigGetAttribute(cfg, using_GMI_OH,   label='using_GMI_OH:',   __RC__)
 call ESMF_ConfigGetAttribute(cfg, using_GMI_NO3,  label='using_GMI_NO3:',  __RC__)

 call ESMF_ConfigGetAttribute(cfg, using_ACHEM_pSO2_OCS,  label='using_ACHEM_pSO2_OCS:',  __RC__)
 call ESMF_ConfigGetAttribute(cfg, tmpchar, label='nei_boundingbox:', RC=status)
 if (status == 0) then
     doing_NEI = .true.
 else
     doing_NEI = .false.
 end if


 if (doing_NEI) then
    call MAPL_AddImportSpec(GC,             &
         SHORT_NAME = 'SU_NEI_SRC1', &
         LONG_NAME  = 'source species'  ,   &
         UNITS      = '1',                  &
         DIMS       = MAPL_DimsHorzOnly,    &
         VLOCATION  = MAPL_VLocationNone,   &
         RESTART    = MAPL_RestartSkip,     &
         RC         = STATUS)
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,             &
         SHORT_NAME = 'SU_NEI_SRC2', &
         LONG_NAME  = 'source species'  ,   &
         UNITS      = '1',                  &
         DIMS       = MAPL_DimsHorzOnly,    &
         VLOCATION  = MAPL_VLocationNone,   &
         RESTART    = MAPL_RestartSkip,     &
         RC         = STATUS)
    VERIFY_(STATUS)
 end if
    
 if (.not. using_GMI_H2O2) then
    call MAPL_AddImportSpec(GC,               &
         SHORT_NAME = 'SU_H2O2'//trim(iname), &
         LONG_NAME  = 'source species'  ,     &
         UNITS      = '1',                    &
         DIMS       = MAPL_DimsHorzVert,      &
         VLOCATION  = MAPL_VLocationCenter,   &
         RESTART    = MAPL_RestartSkip,       &
         RC         = STATUS)
    VERIFY_(STATUS)
 end if

 if (.not. using_GMI_OH) then
    call MAPL_AddImportSpec(GC,             &
         SHORT_NAME = 'SU_OH'//trim(iname), &
         LONG_NAME  = 'source species'  ,   &
         UNITS      = '1',                  &
         DIMS       = MAPL_DimsHorzVert,    &
         VLOCATION  = MAPL_VLocationCenter, &
         RESTART    = MAPL_RestartSkip,     &
         RC         = STATUS)
         VERIFY_(STATUS)
 end if

 if (.not. using_GMI_NO3) then
    call MAPL_AddImportSpec(GC,              &
         SHORT_NAME = 'SU_NO3'//trim(iname), &
         LONG_NAME  = 'source species'  ,    &
         UNITS      = '1',                   &
         DIMS       = MAPL_DimsHorzVert,     &
         VLOCATION  = MAPL_VLocationCenter,  &
         RESTART    = MAPL_RestartSkip,      &
         RC         = STATUS)
    VERIFY_(STATUS)
 end if

 call MAPL_AddImportSpec(GC,            &
      SHORT_NAME = 'SU_BIOMASS'//iname, &
      LONG_NAME  = 'source species'  ,  &
      UNITS      = '1',                 &
      DIMS       = MAPL_DimsHorzOnly,   &
      VLOCATION  = MAPL_VLocationNone,  &
      RESTART    = MAPL_RestartSkip,    &
      RC         = STATUS)
 VERIFY_(STATUS)

 call MAPL_AddImportSpec(GC,             & 
      SHORT_NAME = 'SU_ANTHROL1'//iname, &
      LONG_NAME  = 'source species'  ,   &
      UNITS      = '1',                  &
      DIMS       = MAPL_DimsHorzOnly,    &
      VLOCATION  = MAPL_VLocationNone,   &
      RESTART    = MAPL_RestartSkip,     &
      RC         = STATUS)
 VERIFY_(STATUS)

 call MAPL_AddImportSpec(GC,             &
      SHORT_NAME = 'SU_ANTHROL2'//iname, &
      LONG_NAME  = 'source species'  ,   &
      UNITS      = '1',                  &
      DIMS       = MAPL_DimsHorzOnly,    &
      VLOCATION  = MAPL_VLocationNone,   &
      RESTART    = MAPL_RestartSkip,     &
      RC         = STATUS)
 VERIFY_(STATUS)

 call MAPL_AddImportSpec(GC,            &
      SHORT_NAME = 'SU_SHIPSO2'//iname, &
      LONG_NAME  = 'source species'  ,  &
      UNITS      = '1',                 & 
      DIMS       = MAPL_DimsHorzOnly,   &
      VLOCATION  = MAPL_VLocationNone,  &
      RESTART    = MAPL_RestartSkip,    &
      RC         = STATUS)
 VERIFY_(STATUS)

 call MAPL_AddImportSpec(GC,            &
      SHORT_NAME = 'SU_SHIPSO4'//iname, &
      LONG_NAME  = 'source species'  ,  &
      UNITS      = '1',                 &
      DIMS       = MAPL_DimsHorzOnly,   &
      VLOCATION  = MAPL_VLocationNone,  &
      RESTART    = MAPL_RestartSkip,    &
      RC         = STATUS)
 VERIFY_(STATUS)

 call MAPL_AddImportSpec(GC,           &
      SHORT_NAME = 'SU_DMSO'//iname,   &
      LONG_NAME  = 'source species'  , &
      UNITS      = '1',                &
      DIMS       = MAPL_DimsHorzOnly,  &
      VLOCATION  = MAPL_VLocationNone, &
      RESTART    = MAPL_RestartSkip,   &
      RC         = STATUS)
 VERIFY_(STATUS)

 call MAPL_AddImportSpec(GC,             &
      SHORT_NAME = 'SU_AIRCRAFT'//iname, &
      LONG_NAME  = 'source species'  ,   &
      UNITS      = '1',                  &
      DIMS       = MAPL_DimsHorzVert,    &
      VLOCATION  = MAPL_VLocationCenter, &
      RESTART    = MAPL_RestartSkip,     &
      RC         = STATUS)
 VERIFY_(STATUS)

 call MAPL_AddImportSpec(GC, &
      SHORT_NAME = 'SU_AVIATION_LTO'//trim(iname), &
      LONG_NAME  = 'su_aviation_lto' , &
      UNITS      = '1',                &
      DIMS       = MAPL_DimsHorzOnly,  &
      VLOCATION  = MAPL_VLocationNone, &
      RESTART    = MAPL_RestartSkip,   &
      RC         = STATUS)
 VERIFY_(STATUS)

 call MAPL_AddImportSpec(GC, &
      SHORT_NAME = 'SU_AVIATION_CDS'//trim(iname), &
      LONG_NAME  = 'su_aviation_cds' , &
      UNITS      = '1',                &
      DIMS       = MAPL_DimsHorzOnly,  &
      VLOCATION  = MAPL_VLocationNone, &
      RESTART    = MAPL_RestartSkip,   &
      RC         = STATUS)
 VERIFY_(STATUS)

 call MAPL_AddImportSpec(GC, &
      SHORT_NAME = 'SU_AVIATION_CRS'//trim(iname), &
      LONG_NAME  = 'su_aviation_crs' , &
      UNITS      = '1',                &
      DIMS       = MAPL_DimsHorzOnly,  &
      VLOCATION  = MAPL_VLocationNone, &
      RESTART    = MAPL_RestartSkip,   &
      RC         = STATUS)
 VERIFY_(STATUS)
 

 RETURN_(ESMF_SUCCESS)
 end subroutine SU_GridCompSetServices1_

!--------------------------------------------------------------------------

!                      Single Instance Methods

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  SU_GridCompInitialize --- Initialize SU_GridComp
!
! !INTERFACE:
!

   subroutine SU_GridCompInitialize1_ ( gcSU, w_c, impChem, expChem, &
                                        nymd, nhms, cdt, rc )

! !USES:

  implicit NONE

! !INPUT PARAMETERS:

   type(Chem_Bundle), intent(inout) :: w_c        ! Chemical tracer fields      
   integer, intent(in) :: nymd, nhms           ! time
   real, intent(in) :: cdt                     ! chemistry timestep (secs)

! !OUTPUT PARAMETERS:

   type(SU_GridComp1), intent(inout) :: gcSU   ! Grid Component
   type(ESMF_State), intent(inout)   :: impChem  ! Import State
   type(ESMF_State), intent(inout)   :: expChem  ! Export State
   integer, intent(out) ::  rc                   ! Error return code:
                                                 !  0 - all is well
                                                 !  1 - 

! !DESCRIPTION: Initializes the SU Grid Component. It primarily sets
!               the import state for each active constituent package.
!
! !REVISION HISTORY:
!
!  18Sep2003 da Silva  First crack.
!
!EOP
!-------------------------------------------------------------------------

   character(len=*), parameter :: myname = 'SU_GridCompInitialize1'


   character(len=255) :: rcfilen = 'SU_GridComp.rc'
   integer :: n 
   integer :: i1, i2, im, j1, j2, jm, km, nbins, n1, n2, nbins_rc
   integer, allocatable :: ier(:)
   real :: qmax, qmin
   CHARACTER(LEN=255) :: string
   logical :: NoRegionalConstraint 
   real    :: radius, rhop
   integer :: irhFlag

   rcfilen = gcSU%rcfilen
   gcSU%name = 'SU Constituent Package'
   gcSU%firstRun = .true.


!  Initialize local variables
!  --------------------------
   rc = 0
   i1 = w_c%grid%i1; i2 = w_c%grid%i2; im = w_c%grid%im
   j1 = w_c%grid%j1; j2 = w_c%grid%j2; jm = w_c%grid%jm
   km = w_c%grid%km
   nbins = w_c%reg%n_SU
   n1  = w_c%reg%i_SU
   n2  = w_c%reg%j_SU

!  Check on the number of bins
   if(nbins .ne. 4) then
    rc = 1
    return
   endif


   call init_()
   if ( rc /= 0 ) return

!  Set the bin assignments to the gcSU grid component
   gcSU%nDMS = 1
   gcSU%nSO2 = 2
   gcSU%nSO4 = 3
   gcSU%nMSA = 4


!                       -------------------
!                       Parse resource file
!                       -------------------

!  Load resource file
!  ------------------
   call i90_loadf ( rcfilen, ier(1) )
   if ( ier(1) .ne. 0 ) then
      call final_(10)
      return
   end if

   call i90_label ( 'number_su_classes:', ier(1) )
   nbins_rc = i90_gint ( ier(2) )
   if ( any(ier(1:2) /= 0) ) then
      call final_(20)
      return
   end if
   if ( nbins_rc /= nbins ) then
      call final_(25)
      return
   end if


!  SU sources files
!  ---------------------
   call i90_label ( 'volcano_srcfilen:', ier(1) )
   if ( ier(1) /= 0 ) then
      call final_(31)
      return
   else
      call i90_gtoken ( gcSU%volcano_srcfilen, ier(1) )
      if ( ier(1) /= 0 ) then
         call final_(41)
         return
      end if
   end if

!  Aircraft emissions
!  ------------------
   ier(:) = 0
   call i90_label  ( 'aviation_vertical_layers:', ier(1) )
   gcSU%aviation_layers(1) = i90_gfloat(ier(2))
   gcSU%aviation_layers(2) = i90_gfloat(ier(3))
   gcSU%aviation_layers(3) = i90_gfloat(ier(4))
   gcSU%aviation_layers(4) = i90_gfloat(ier(5))

   if ( any(ier(1:5) /= 0) ) then
         call final_(77)
         return
   end if

!  Handle Point-wise Emission Sources Specified in a Text File
!  -----------------------------------------------------------
   ier(:) = 0
   call i90_label  ( 'point_emissions_srcfilen:',   ier(1) )
   call i90_gtoken ( gcSU%point_emissions_srcfilen, ier(2) )
   if ( ier(1) /= 0 ) then
        gcSU%doing_point_emissions = .FALSE. ! if rc is missing, don't fuss
   else if ( any(ier(2:2) /= 0) ) then
         call final_(42) ! this means point emissions info is messed up, abort
         return
   else
         if ( (index(gcSU%point_emissions_srcfilen,'/dev/null')>0) ) then
               gcSU%doing_point_emissions = .FALSE. ! disable it if no file specified
         else
               gcSU%doing_point_emissions = .TRUE.  ! we are good to go
         end if
   end if

!  Handle NEI08 Emissions
!  ----------------------
   ier(:) = 0
   call i90_label  ( 'nei_bot_srcfilen:',  ier(1) )
   call i90_gtoken ( gcSU%nei_srcfilen(1), ier(2) )
   call i90_label  ( 'nei_top_srcfilen:',  ier(3) )
   call i90_gtoken ( gcSU%nei_srcfilen(2), ier(4) )
   call i90_label  ( 'nei_boundingbox:',   ier(5) )
   gcSU%nei_lon(1) = i90_gfloat(ier(6))
   gcSU%nei_lon(2) = i90_gfloat(ier(7))
   gcSU%nei_lat(1) = i90_gfloat(ier(8))
   gcSU%nei_lat(2) = i90_gfloat(ier(9))
   if ( ier(1) /= 0 ) then
        gcSU%doing_nei = .FALSE. ! if rc is missing, don't fuss
   else if ( any(ier(2:9) /= 0) ) then
         call final_(42) ! this means NEI info is messed up, abort
         return
   else
         if ( (index(gcSU%nei_srcfilen(1),'/dev/null')>0) .or. &
              (index(gcSU%nei_srcfilen(2),'/dev/null')>0) ) then 
               gcSU%doing_nei = .FALSE. ! disable it if no file specified
         else
               gcSU%doing_nei = .TRUE.  ! we are good to go
         end if
   end if

   if ( MAPL_AM_I_ROOT() ) then
    if ( gcSU%doing_nei ) then
      print *, 'SU_GridComp: using NEI08 Emissions over North America'
    else
      print *, 'SU_GridComp: skipping NEI08 Emissions over North America'
    end if
   end if

!                          -------

!  Day of the week to reset tracer to zero
!  ---------------------------------------
   call i90_label ( 'my_day_of_the_week:',ier(1))
   if ( ier(1) /= 0 ) then
        gcSU%myDOW = -1   ! by default never reset tracer to zero
   else
        gcSU%myDOW = i90_gint (ier(1))
        if ( ier(1) /= 0 ) then
           call final_(60)
           return
        end if
   end if

!                          -------
!  Fraction of anthropogenic emissions to SO4
!  ---------------
   call i90_label ( 'so4_anthropogenic_fraction:', ier(1) )
   gcSU%fSO4ant = i90_gfloat ( ier(2) )
   if ( any(ier(1:2) /= 0) ) then
      call final_(50)
      return
   end if


!  Aircraft Fuel Emission Factor
!  ---------------
   call i90_label ( 'aircraft_fuel_emission_factor:', ier(1) )
   gcSU%eAircraftFuel = i90_gfloat ( ier(2) )
   if ( any(ier(1:2) /= 0) ) then
      call final_(52)
      return
   end if


!  Scavenging Efficiency
!  To be used in convtran.F90, this parameter
!  is the scavenging efficiency of the tracer [km -1]
!  ---------------
   call i90_label ( 'fscav:', ier(1) )
   do n = 1, nbins
      w_c%reg%fscav(n1+n-1) = i90_gfloat ( ier(n+1) )
   end do
   if ( any(ier(1:nbins+1) /= 0) ) then
      call final_(53)
      return
   end if

!  Particle radius
!  ---------------
   call i90_label ( 'particle_radius:', ier(1) )
   do n = 1, nbins
      radius           = i90_gfloat ( ier(n+1) )
      gcSU%radius(n)   = radius  ! save radius in [um]
!      w_c%qa(n1+n-1)%r = radius * 1.e-6 !radius in [m]
   end do
   if ( any(ier(1:nbins+1) /= 0) ) then
      call final_(50)
      return
   end if

!  Particle affected by relative humidity?
!  ---------------
   call i90_label ( 'rhFlag:', ier(1) )
   irhFlag                    = i90_gint ( ier(2) )
   gcSU%rhFlag                = irhFlag
   w_c%qa(n1+n-1)%irhFlag     = irhFlag
   if ( any(ier(1:2) /= 0) ) then
      call final_(50)
      return
   end if

!  Particle density
!  To be used in droplet activation code
!  ---------------
   call i90_label ( 'particle_density:', ier(1) )
   do n = 1, nbins
      w_c%reg%rhop(n1+n-1)  = i90_gfloat ( ier(n+1) )
      gcSU%rhop(n)          = w_c%reg%rhop(n1+n-1)
   end do
   if ( any(ier(1:nbins+1) /= 0) ) then
      call final_(54)
      return
   end if
!                          -------

!  Number median radius
!  To be used in droplet activation code
!  ---------------
   call i90_label ( 'particle_radius_number:', ier(1) )
   do n = 1, nbins
      w_c%reg%rmed(n1+n-1)  = i90_gfloat ( ier(n+1) ) * 1e-6
   end do
   if ( any(ier(1:nbins+1) /= 0) ) then
      call final_(55)
      return
   end if
!                          -------

!  Sigma (lognormal mode width)
!  To be used in droplet activation code
!  ---------------
   call i90_label ( 'sigma:', ier(1) )
   do n = 1, nbins
      w_c%reg%sigma(n1+n-1)  = i90_gfloat ( ier(n+1) )
   end do
   if ( any(ier(1:nbins+1) /= 0) ) then
      call final_(56)
      return
   end if
!                          -------

!  Number to mass conversion factor
!  To be used in droplet activation code
!  ---------------
   call i90_label ( 'fnum:', ier(1) )
   do n = 1, nbins
      w_c%reg%fnum(n1+n-1)  = i90_gfloat ( ier(n+1) )
   end do
   if ( any(ier(1:nbins+1) /= 0) ) then
      call final_(57)
      return
   end if
!                          -------

!  Molecular weight
!  To be used in droplet activation code
!  ---------------
   call i90_label ( 'molecular_weight:', ier(1) )
   do n = 1, nbins
      w_c%reg%molwght(n1+n-1)  = i90_gfloat ( ier(n+1) )
   end do
   if ( any(ier(1:nbins+1) /= 0) ) then
      call final_(58)
      return
   end if
!                          -------


! Switches to select predicted OH H2O2 NO3 from the 
! GMI Combined Stratosphere Troposphere Chemical Mechanism
! --------------------------------------------------------
   gcSU%using_GMI_OH = .FALSE.
   CALL I90_Label("using_GMI_OH:",ier(1))
   IF(ier(1) /= 0) THEN
    CALL final_(81)
    RETURN
   ELSE
    CALL I90_GToken(string,ier(1))
    IF(ier(1) /= 0) THEN
     CALL final_(82)
     RETURN
    END IF
    IF(TRIM(string) == "yes" .AND. w_c%reg%doing_GMI) gcSU%using_GMI_OH = .TRUE.
   END IF

   gcSU%using_GMI_NO3 = .FALSE.
   CALL I90_Label("using_GMI_NO3:",ier(1))
   IF(ier(1) /= 0) THEN
    CALL final_(83)
    RETURN
   ELSE
    CALL I90_GToken(string,ier(1))
    IF(ier(1) /= 0) THEN
     CALL final_(84)
     RETURN
    END IF
    IF(TRIM(string) == "yes" .AND. w_c%reg%doing_GMI) gcSU%using_GMI_NO3 = .TRUE.
   END IF

   gcSU%using_GMI_H2O2 = .FALSE.
   CALL I90_Label("using_GMI_H2O2:",ier(1))
   IF(ier(1) /= 0) THEN
    CALL final_(85)
    RETURN
   ELSE
    CALL I90_GToken(string,ier(1))
    IF(ier(1) /= 0) THEN
     CALL final_(86)
     RETURN
    END IF
    IF(TRIM(string) == "yes" .AND. w_c%reg%doing_GMI) gcSU%using_GMI_H2O2 = .TRUE.
   END IF

   gcSU%export_H2O2 = .FALSE.
   CALL I90_Label("export_H2O2:",ier(1))
   IF(ier(1) /= 0) THEN
    CALL final_(87)
    RETURN
   ELSE
    CALL I90_GToken(string,ier(1))
    IF(ier(1) /= 0) THEN
     CALL final_(88)
     RETURN
    END IF
    IF(TRIM(string) == "yes" .AND. w_c%reg%doing_GMI) gcSU%export_H2O2 = .TRUE.
   END IF

! Switches to select import of production of SO2 from OCS
! provided in ACHEM mechanism
! --------------------------------------------------------
   gcSU%using_ACHEM_pSO2_OCS = .FALSE.
   CALL I90_Label("using_ACHEM_pSO2_OCS:",ier(1))
   IF(ier(1) /= 0) THEN
    CALL final_(89)
    RETURN
   ELSE
    CALL I90_GToken(string,ier(1))
    IF(ier(1) /= 0) THEN
     CALL final_(90)
     RETURN
    END IF
    IF(TRIM(string) == "yes" .AND. w_c%reg%doing_OCS) gcSU%using_ACHEM_pSO2_OCS = .TRUE.
   END IF

   IF(MAPL_AM_I_ROOT()) THEN
   PRINT *," "
   PRINT *,TRIM(myname)//":"
   PRINT *," Using GMI   OH: ",gcSU%using_GMI_OH
   PRINT *," Using GMI  NO3: ",gcSU%using_GMI_NO3
   PRINT *," Using GMI H2O2: ",gcSU%using_GMI_H2O2
   PRINT *," Using ACHEM pSO2_OCS: ", gcSU%using_ACHEM_pSO2_OCS
   PRINT *," Exporting updated H2O2 to GMI: ",gcSU%export_H2O2
   PRINT *," "
   END IF

!                          -------


!                          -------

!  Set the gram molecular weights of the species
!  ---------------------------------------------
   gcSU%fMassSulfur = 32.0
   gcSU%fMassSO2    = 64.0
   gcSU%fMassSO4    = 96.0
   gcSU%fMassDMS    = 62.0
   gcSU%fMassMSA    = 96.0

!  Initialize date for boundary conditions
!  ---------------------------------------
   gcSU%nymd = -1          ! nothing read yet
   gcSU%nymd_oxidants = -1

!  Grab the region string.
!  -----------------------
   call i90_label ( 'SU_regions_indices:', ier(1) )
   CALL I90_gtoken( gcSU%regionsString, ier(2) )
   IF( ANY(ier(1:2) < 0 ) ) THEN
    CALL final_(51)
    RETURN
   END IF

!  Is this instantiation a global case?
!  -----------------------------------
   IF(gcSU%regionsString(1:2) == "-1") THEN
    NoRegionalConstraint = .TRUE.
   ELSE
    SELECT CASE (ESMF_UtilStringLowerCase(gcSU%regionsString(1:2)))
     CASE ("gl") 
      NoRegionalConstraint = .TRUE.
     CASE ("al") 
      NoRegionalConstraint = .TRUE.
     CASE DEFAULT
      NoRegionalConstraint = .FALSE.
    END SELECT
   END IF

!  Set regionsString to "-1" for the global case
!  ---------------------------------------------
   IF(NoRegionalConstraint) gcSU%regionsString = "-1"

   IF(MAPL_AM_I_ROOT()) THEN
    IF(NoRegionalConstraint) THEN
     PRINT *,myname,": This instantiation has no regional constraints."
    ELSE
     PRINT *,myname,": This instantiation is regionally constrained."
     PRINT *,myname,": List of region numbers included: ",TRIM(gcSU%regionsString)
    END IF
   END IF

!  All done
!  --------
   call i90_release()
   deallocate(ier)

   return


CONTAINS

   subroutine init_()
   integer ios, nerr
   nerr = max ( 32, nbins+1 )
   allocate ( gcSU%so2biomass_src(i1:i2,j1:j2), gcSU%so2anthro_l1_src(i1:i2,j1:j2), &
              gcSU%so2biomass_src_(i1:i2,j1:j2), &
              gcSU%so2anthro_l2_src(i1:i2,j1:j2), gcSU%dmso_conc(i1:i2,j1:j2), &
              gcSU%so2ship_src(i1:i2,j1:j2), gcSU%so4ship_src(i1:i2,j1:j2), &
              gcSU%oh_conc(i1:i2,j1:j2,km), gcSU%no3_mr(i1:i2,j1:j2,km), &
              gcSU%h2o2_mr(i1:i2,j1:j2,km), gcSU%h2o2_int(i1:i2,j1:j2,km), &
              gcSU%aircraft_fuel_src(i1:i2,j1:j2,km), &
              gcSU%aviation_lto_src(i1:i2,j1:j2), &
              gcSU%aviation_cds_src(i1:i2,j1:j2), &
              gcSU%aviation_crs_src(i1:i2,j1:j2), &
              gcSU%regionMask(i1:i2,j1:j2), &
              gcSU%radius(nbins), gcSU%rhop(nbins), &
              ier(nerr), stat=ios )
   if ( ios /= 0 ) rc = 100
   end subroutine init_

   subroutine final_(ierr)
   integer :: ierr
   integer ios
   deallocate ( gcSU%so2biomass_src, gcSU%so2anthro_l1_src, gcSU%so2anthro_l2_src, &
                gcSU%so2biomass_src_, &
                gcSU%dmso_conc, gcSU%oh_conc, gcSU%no3_mr, &
                gcSU%so2ship_src, gcSU%so4ship_src, &
                gcSU%h2o2_mr, gcSU%h2o2_int, gcSU%aircraft_fuel_src, &
                gcSU%aviation_lto_src, gcSU%aviation_cds_src, gcSU%aviation_crs_src, &
                gcSU%regionMask, gcSU%radius, gcSU%rhop, &
                ier, stat=ios )

   call i90_release()
   rc = ierr
   end subroutine final_

   end subroutine SU_GridCompInitialize1_

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  SU_GridCompRun1_ --- The Chem Driver, run phase 1
!
! !INTERFACE:
!

   subroutine SU_GridCompRun1_ ( gcSU, w_c, impChem, expChem, &
                                 nymd, nhms, cdt, rc )

! !USES:

  implicit NONE

! !INPUT/OUTPUT PARAMETERS:

   type(SU_GridComp1), intent(inout) :: gcSU   ! Grid Component
   type(Chem_Bundle), intent(inout)  :: w_c    ! Chemical tracer fields   

! !INPUT PARAMETERS:

   type(ESMF_State), intent(inout) :: impChem  ! Import State
   integer, intent(in) :: nymd, nhms           ! time
   real, intent(in) :: cdt                     ! chemistry timestep (secs)

! !OUTPUT PARAMETERS:

   type(ESMF_State), intent(inout) :: expChem  ! Export State
   integer, intent(out) ::  rc                 ! Error return code:
                                               !  0 - all is well
                                               !  1 -
 
! !DESCRIPTION: This routine implements the so-called SU Driver. That 
!               is, adds chemical tendencies to each of the constituents,
!  Note: water wapor, the first constituent is not considered a chemical
!  constituents.
!
! !REVISION HISTORY:
!
!  18Sep2003 da Silva  First crack.
!
!EOP
!-------------------------------------------------------------------------

   character(len=*), parameter :: myname = 'SU_GridCompRun1_'
   character(len=*), parameter :: Iam = myname

   integer :: i1, i2, im, j1, j2, jm, nbins, n1, n2, km, n
   integer :: i, j, ijl, ijkl, ijk1l
   real :: qmax, qmin

!  Input fields from fvGCM
!  -----------------------
   real, pointer, dimension(:,:)   ::  pblh, oro, & 
                                       u10m, v10m, hsurf
   real, pointer, dimension(:,:,:) ::  tmpu, rhoa, hghte


   real, pointer                         :: var2d(:,:) => null()

   integer                               :: idow
   character(len=3)                      :: cdow


#define EXPORT        expChem
#define iNAME         TRIM(gcSU%iname)

#define ptrSUEM       SU_emis

#define ptrSO4EMAN    SU_SO4eman
#define ptrSO2EMAN    SU_SO2eman
#define ptrSO2EMBB    SU_SO2embb
#define ptrSO2EMVN    SU_SO2emvn
#define ptrSO2EMVE    SU_SO2emve

   integer :: STATUS

!  Indices for point emissions
   integer, pointer, dimension(:)  :: iPoint, jPoint
   real, dimension(w_c%grid%km)    :: point_column_emissions
   integer                         :: ios, ii

#include "SU_GetPointer___.h"

!  Initialize local variables
!  --------------------------
   rc = 0
   i1 = w_c%grid%i1; i2 = w_c%grid%i2; im = w_c%grid%im
   j1 = w_c%grid%j1; j2 = w_c%grid%j2; jm = w_c%grid%jm

   km    = w_c%grid%km
   nbins = w_c%reg%n_SU
   n1    = w_c%reg%i_SU
   n2    = w_c%reg%j_SU

   ijl   = ( i2 - i1 + 1 ) * ( j2 - j1 + 1 )
   ijkl  = ijl * km
   ijk1l = ijl * (km+1)

!  Reset tracer to zero at 0Z on specific day of week
!  --------------------------------------------------
   idow = Chem_UtilIdow(nymd)
   if ( (nhms==0) .and. (idow == gcSU%myDOW) ) then
        cdow = Chem_UtilCdow(nymd)
        do n = n1, n2
           w_c%qa(n)%data3d(i1:i2,j1:j2,1:km) = tiny(1.) ! avoid division by zero
        end do
        if ( MAPL_AM_I_ROOT() ) then
           print *, '<> SU '//cdow//' tracer being set to zero on ', nymd, nhms
        end if
   end if


   call MAPL_GetPointer(impChem, var2d, 'SU_regionMask', __RC__)
   gcSU%regionMask = var2d


   call SulfateUpdateEmissions (impChem, iNAME, i1, i2, im, j1, j2, jm, km, cdt, &
                                nymd, nhms, &
                                w_c%grid_esmf, w_c%grid%lon, w_c%grid%lat, &
                                gcSU%nymd, &
                                w_c%diurnal_bb, &
                                gcSU%so2biomass_src, gcSU%so2biomass_src_, &
                                gcSU%so2anthro_l1_src, &
                                gcSU%so2anthro_l2_src, &
                                gcSU%so2ship_src, &
                                gcSU%so4ship_src, &
                                gcSU%dmso_conc, &
                                gcSU%aircraft_fuel_src, &
                                gcSU%aviation_lto_src, &
                                gcSU%aviation_cds_src, &
                                gcSU%aviation_crs_src, &
                                gcSU%volcano_srcfilen, &
                                gcSU%nvolc, gcSU%vLat, gcSU%vLon, &
                                gcSU%vElev, gcSU%vCloud, gcSU%vSO2, &
                                gcSU%vStart, gcSU%vEnd, &
                                doing_NEI=gcSU%doing_NEI, & 
                                nei_hour=gcSU%nei_hour, &
                                nei_year=gcSU%nei_year, &
                                nei_srcfilen=gcSU%nei_srcfilen, &
                                nei_lon=gcSU%nei_lon, &
                                nei_lat=gcSU%nei_lat, &
                                lons=w_c%grid%lon, &
                                lats=w_c%grid%lat, &
                                maskString=trim(gcSU%regionsString), &
                                gridMask=gcSU%regionMask, &
                                rc=STATUS)
   VERIFY_(STATUS)

!  Read any pointwise emissions, if requested (hardcoded will go to sulfate)
!  -------------------------------------------------------------------------
   if(gcSU%doing_point_emissions) then
    call Chem_UtilPointEmissions( nymd, gcSU%point_emissions_srcfilen, &
                                  gcSU%nPts, gcSU%pLat, gcSU%pLon, &
                                  gcSU%pBase, gcSU%pTop, gcSU%pEmis, &
                                  gcSU%pStart, gcSU%pEnd )

!   In case pStart or pEnd were not specified in the file set to defaults
    where(gcSU%pStart < 0) gcSU%pStart = 000000
    where(gcSU%pEnd < 0)   gcSU%pEnd   = 240000
   endif

                           
!  Get 2D Imports
!  --------------
   call MAPL_GetPointer ( impChem, pblh,     'ZPBL',     __RC__ )
   call MAPL_GetPointer ( impChem, oro,      'LWI',      __RC__ )
   call MAPL_GetPointer ( impChem, u10m,     'U10M',     __RC__ )
   call MAPL_GetPointer ( impChem, v10m,     'V10M',     __RC__ )

!  Get 3D Imports
!  --------------
   call MAPL_GetPointer ( impChem, tmpu,     'T',        __RC__ )
   call MAPL_GetPointer ( impChem, rhoa,     'AIRDENS',  __RC__ )
   call MAPL_GetPointer ( impChem, hghte,    'ZLE',      __RC__ )


!  Unlike GEOS-4 hghte is defined for km+1
!  ---------------------------------------
   hsurf => hghte(i1:i2,j1:j2,km) ! in GEOS-5 hghte is in [0,km]

#ifdef DEBUG

   call pmaxmin('SU: oro        ', oro     , qmin, qmax, ijl,1, 1. )
   call pmaxmin('SU: u10m       ', u10m    , qmin, qmax, ijl,1, 1. )
   call pmaxmin('SU: v10m       ', v10m    , qmin, qmax, ijl,1, 1. )
   call pmaxmin('SU: pblh       ', pblh    , qmin, qmax, ijl,1, 1. )
   call pmaxmin('SU: hsurf      ', hsurf   , qmin, qmax, ijl,1, 1. )

   call pmaxmin('SU: tmpu       ', tmpu    , qmin, qmax, ijkl,1, 1. )
   call pmaxmin('SU: rhoa       ', rhoa    , qmin, qmax, ijkl,1, 1. )
   call pmaxmin('SU: hghte      ', hghte   , qmin, qmax, ijkl,1, 1. )

#endif

!  SU Source
!  -----------
   call SulfateDistributeEmissions ( i1, i2, j1, j2, km, nbins, cdt, nymd, nhms, &
                                     gcSU%fSO4ant, &
                                     gcSU%eAircraftFuel, &
                                     gcSU%so2anthro_l1_src, gcSU%so2anthro_l2_src, &
                                     gcSU%so2biomass_src, gcSU%dmso_conc, &
                                     gcSU%so2ship_src, gcSU%so4ship_src, &
                                     gcSU%aircraft_fuel_src, &
                                     gcSU%nvolc, gcSU%vLat, gcSU%vLon, &
                                     gcSU%vElev, gcSU%vCloud, gcSU%vSO2, &
                                     gcSU%vStart, gcSU%vEnd, &
                                     w_c%qa(n1+gcSU%nDMS-1)%data3d, &
                                     w_c%qa(n1+gcSU%nSO2-1)%data3d, &
                                     w_c%qa(n1+gcSU%nSO4-1)%data3d, &
                                     oro, u10m, v10m, hsurf, hghte, pblh, &
                                     tmpu, rhoa, w_c%delp, &
                                     w_c%grid%cell_area, &
                                     w_c%grid_esmf, &
                                     SU_emis, &
                                     SU_SO4eman, SU_SO2eman, SU_SO2embb, &
                                     SU_SO2emvn, SU_SO2emve, &
                                     rc, &
                                     aviation_layers=gcSU%aviation_layers,   &
                                     aviation_lto_src=gcSU%aviation_lto_src, &
                                     aviation_cds_src=gcSU%aviation_cds_src, &
                                     aviation_crs_src=gcSU%aviation_crs_src)

#ifdef DEBUG
   do n = n1, n2
      call pmaxmin('SU: q_emi', w_c%qa(n)%data3d(i1:i2,j1:j2,1:km), qmin, qmax, &
                    ijl, km, 1. )
   end do
#endif


!  Distribute pointwise sources if requested
!  -----------------------------------------
   POINTWISE_SOURCES: if( gcSU%doing_point_emissions .and. gcSU%nPts > 0) then

!    Get indices for point emissions
!    -------------------------------
     allocate(iPoint(gcSU%nPts), jPoint(gcSU%nPts), stat=ios)

     call MAPL_GetHorzIJIndex(gcSU%nPts, iPoint, jPoint, &
                              grid = w_c%grid_esmf,      &
                              lon  = gcSU%pLon/radToDeg, &
                              lat  = gcSU%pLat/radToDeg, &
                              rc   = rc)

     if ( rc /= 0 ) call die(myname,'cannot get indices for point emissions')

     do ii = 1, gcSU%nPts
      i = iPoint(ii)
      j = jPoint(ii)
      if( i<1 .OR. j<1 )              cycle    ! point emission not in this sub-domain
!      if( gcSU%regionMask(i,j) == 0 ) cycle    ! masked by region mask

!     Emissions not occurring in current time step
!     --------------------------------------------
      if(nhms < gcSU%pStart(ii) .or. nhms >= gcSU%pEnd(ii)) cycle

      call distribute_point_emissions(w_c%delp(i,j,:), rhoa(i,j,:), &
                                      gcSU%pBase(ii), gcSU%pTop(ii), gcSU%pEmis(ii), &
                                      point_column_emissions, km)
      w_c%qa(n1+gcSU%nSO4-1)%data3d(i,j,:) = w_c%qa(n1+gcSU%nSO4-1)%data3d(i,j,:) & 
         + cdt * grav / w_c%delp(i,j,:) &
               * point_column_emissions / w_c%grid%cell_area(i,j)

     enddo

     deallocate(iPoint, jPoint, stat=ios)

   endif POINTWISE_SOURCES


   return

CONTAINS

!  Abstracted from distribute_aviation_emissions, but called per column
   subroutine distribute_point_emissions(delp, rhoa, z_bot, z_top, emissions_point, &
                                         emissions, km)

    implicit none

    integer, intent(in) :: km

    real, dimension(:), intent(in) :: delp
    real, dimension(:), intent(in) :: rhoa
    real,               intent(in) :: emissions_point
    real, intent(in)                   :: z_bot
    real, intent(in)                   :: z_top
    real, dimension(:), intent(out):: emissions
    
!   local
    integer :: k
    integer :: k_bot, k_top
    real    :: z_
    real, dimension(km) :: z, dz, w_
    
!   find level height
    z = 0.0
    z_= 0.0 

    do k = km, 1, -1
       dz(k) = delp(k)/rhoa(k)/grav
       z_    = z_ + dz(k)
       z(k)  = z_
    end do

!   find the bottom level
    do k = km, 1, -1
       if (z(k) >= z_bot) then
           k_bot = k
           exit
       end if
    end do
            
!   find the top level
    do k = k_bot, 1, -1
       if (z(k) >= z_top) then
           k_top = k
           exit
       end if
    end do

!   find the weights
    w_ = 0

!   if (k_top > k_bot) then
!       need to bail - something went wrong here
!   end if

    if (k_bot .eq. k_top) then
        w_(k_bot) = z_top - z_bot
    else
     do k = k_bot, k_top, -1
        if ((k < k_bot) .and. (k > k_top)) then
             w_(k) = dz(k)
        else
             if (k == k_bot) then
                 w_(k) = (z(k) - z_bot)
             end if

             if (k == k_top) then
                 w_(k) = z_top - (z(k)-dz(k))
             end if
        end if
     end do
    end if
           
!   distribute emissions in the vertical 
    emissions(:) = (w_ / sum(w_)) * emissions_point

    end subroutine distribute_point_emissions


 end subroutine SU_GridCompRun1_



!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  SU_GridCompRun2_ --- The Chem Driver, run phase 2
!
! !INTERFACE:
!

   subroutine SU_GridCompRun2_ ( gcSU, w_c, impChem, expChem, &
                                 nymd, nhms, cdt, rc )

! !USES:

  implicit NONE

! !INPUT/OUTPUT PARAMETERS:

   type(SU_GridComp1), intent(inout) :: gcSU   ! Grid Component
   type(Chem_Bundle), intent(inout)  :: w_c    ! Chemical tracer fields   

! !INPUT PARAMETERS:

   type(ESMF_State), intent(inout) :: impChem    ! Import State
   integer, intent(in) :: nymd, nhms          ! time
   real, intent(in) :: cdt                    ! chemistry timestep (secs)

! !OUTPUT PARAMETERS:

   type(ESMF_State), intent(inout) :: expChem   ! Export State
   integer, intent(out) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 -
 
! !DESCRIPTION: This routine implements the so-called SU Driver. That 
!               is, adds chemical tendencies to each of the constituents,
!  Note: water wapor, the first constituent is not considered a chemical
!  constituents.
!
! !REVISION HISTORY:
!
!  18Sep2003 da Silva  First crack.
!
!EOP
!-------------------------------------------------------------------------

   character(len=*), parameter :: myname = 'SU_GridCompRun2_'
   character(len=*), parameter :: Iam = myname

   integer :: i1, i2, im, j1, j2, jm, nbins, n1, n2, km, n
   integer :: k, ijl, ijkl, ijk1l
   real :: qmax, qmin
   real, pointer :: SU_radius(:), SU_rhop(:)
   logical :: KIN

!  Input fields from fvGCM
!  -----------------------
   real, pointer, dimension(:,:)   ::  frlake, frocean, frseaice, &
                                       pblh, oro, shflux, ustar, precc, &
                                       precl, u10m, v10m, hsurf, z0h
   real, pointer, dimension(:,:,:) ::  tmpu, cloud, rhoa, u, v, hghte, ple
   real, pointer, dimension(:,:,:) ::  pfllsan, pfilsan

   REAL, POINTER, DIMENSION(:,:,:) ::  GMI_H2O2mr, GMI_OHmr, GMI_NO3mr, ACHEM_PSO2_OCS
   real, pointer, dimension(:,:,:) ::  xoh, xno3, xh2o2

!  Additional needs for GOCART convective diagnostic
   real, pointer, dimension(:,:,:)       ::  cmfmc, qlcn, qicn, dtrain
   real, pointer, dimension(:,:)         ::  area
   real*8, allocatable, dimension(:,:,:) ::  cmfmc_, qccu_, dtrain_, &
                                             airmass_, airmol_, vud_, &
                                             delz_, delp_, h2o2_, tmpu_, ple_
   real*8, allocatable                   ::  tc_(:,:,:,:), bcnv_(:,:,:)
   real*8, allocatable                   ::  area_(:,:), frlake_(:,:), &
                                             frocean_(:,:), frseaice_(:,:)
   integer*4                             ::  icdt

   real, pointer                         :: var2d(:,:) => null()


#define EXPORT        expChem
#define iNAME         TRIM(gcSU%iname)

#define ptrSUWT       SU_wet
#define ptrSUSV       SU_conv
#define ptrSUDP       SU_dep
#define ptrSUSD       SU_set

#define ptrSUPSO2     SU_PSO2
#define ptrSUPSO4     SU_PSO4
#define ptrSUPSO4G    SU_PSO4g
#define ptrSUPSO4AQ   SU_PSO4aq
#define ptrSUPSO4WT   SU_PSO4wet
#define ptrSUPMSA     SU_PMSA

#define ptrSO2SMASS   SU_SO2sfcmass
#define ptrSO2CMASS   SU_SO2colmass
#define ptrSO4SMASS   SU_SO4sfcmass
#define ptrSO4CMASS   SU_SO4colmass
#define ptrDMSSMASS   SU_DMSsfcmass
#define ptrDMSCMASS   SU_DMScolmass
#define ptrMSASMASS   SU_MSAsfcmass
#define ptrMSACMASS   SU_MSAcolmass
#define ptrSUCONC     SU_conc
#define ptrSUEXTCOEF  SU_extcoef
#define ptrSUSCACOEF  SU_scacoef
#define ptrSUANGSTR   SU_angstrom
#define ptrSUFLUXU    SU_fluxu
#define ptrSUFLUXV    SU_fluxv
#define ptrSO4MASS    SU_so4mass
#define ptrSUEXTTAU   SU_exttau
#define ptrSUSCATAU   SU_scatau 
#define ptrSO4SAREA   SU_sarea
#define ptrSO4SNUM    SU_snum

   integer :: STATUS

#include "SU_GetPointer___.h"

!  Initialize local variables
!  --------------------------
   rc = 0
   i1 = w_c%grid%i1; i2 = w_c%grid%i2; im = w_c%grid%im
   j1 = w_c%grid%j1; j2 = w_c%grid%j2; jm = w_c%grid%jm

   km    = w_c%grid%km
   nbins = w_c%reg%n_SU
   n1    = w_c%reg%i_SU
   n2    = w_c%reg%j_SU

   ijl   = ( i2 - i1 + 1 ) * ( j2 - j1 + 1 )
   ijkl  = ijl * km
   ijk1l = ijl * (km+1)


   allocate(xoh(i1:i2,j1:j2,km), xno3(i1:i2,j1:j2,km), xh2o2(i1:i2,j1:j2,km), __STAT__)

   call MAPL_GetPointer(impChem, var2d, 'SU_regionMask', __RC__)
   gcSU%regionMask = var2d

   if (gcSU%firstRun) then
      gcSU%h2o2_mr  = MAPL_UNDEF
      gcSU%h2o2_int = MAPL_UNDEF
      xh2o2         = MAPL_UNDEF
      gcSU%firstRun = .false.
   end if

   xoh   = 0.0
   xno3  = 0.0
   xh2o2 = gcSU%h2o2_int


!  Get 2D Imports
!  --------------
   call MAPL_GetPointer ( impChem, pblh,     'ZPBL',     __RC__ )
   call MAPL_GetPointer ( impChem, oro,      'LWI',      __RC__ )
   call MAPL_GetPointer ( impChem, shflux,   'SH',       __RC__ )
   call MAPL_GetPointer ( impChem, ustar,    'USTAR',    __RC__ )
   call MAPL_GetPointer ( impChem, precc,    'CN_PRCP',  __RC__ )
   call MAPL_GetPointer ( impChem, precl,    'NCN_PRCP', __RC__ )
   call MAPL_GetPointer ( impChem, u10m,     'U10M',     __RC__ )
   call MAPL_GetPointer ( impChem, v10m,     'V10M',     __RC__ )
   call MAPL_GetPointer ( impChem, z0h,      'Z0H',      __RC__ )
   call MAPL_GetPointer ( impChem, area,     'AREA',     __RC__ )
   call MAPL_GetPointer ( impChem, frocean,  'FROCEAN',  __RC__ )
   call MAPL_GetPointer ( impChem, frseaice, 'FRACI',    __RC__ )
   call MAPL_GetPointer ( impChem, frlake,   'FRLAKE',   __RC__ )

!  Get 3D Imports
!  --------------
   call MAPL_GetPointer ( impChem, tmpu,     'T',        __RC__ )
   call MAPL_GetPointer ( impChem, cloud,    'FCLD',     __RC__ )
   call MAPL_GetPointer ( impChem, rhoa,     'AIRDENS',  __RC__ )
   call MAPL_GetPointer ( impChem, u,        'U',        __RC__ )
   call MAPL_GetPointer ( impChem, v,        'V',        __RC__ )
   call MAPL_GetPointer ( impChem, hghte,    'ZLE',      __RC__ )
   call MAPL_GetPointer ( impChem, ple,      'PLE',      __RC__ )
   call MAPL_GetPointer ( impChem, qlcn,     'QLCN',     __RC__ )
   call MAPL_GetPointer ( impChem, qicn,     'QICN',     __RC__ )
   call MAPL_GetPointer ( impChem, cmfmc,    'CNV_MFC',  __RC__ )
   call MAPL_GetPointer ( impChem, dtrain,   'CNV_MFD',  __RC__ )
   call MAPL_GetPointer ( impChem, pfllsan,  'PFL_LSAN', __RC__ )
   call MAPL_GetPointer ( impChem, pfilsan,  'PFI_LSAN', __RC__ )

!  Oxidants from GMICHEM.  Get pointers first ...
!  ----------------------------------------------
   if (  gcSU%using_GMI_OH) call MAPL_GetPointer(impChem,   GMI_OHmr,   'OH', __RC__)
   if ( gcSU%using_GMI_NO3) call MAPL_GetPointer(impChem,  GMI_NO3mr,  'NO3', __RC__)
   if (gcSU%using_GMI_H2O2) call MAPL_GetPointer(impChem, GMI_H2O2mr, 'H2O2', __RC__)

!  Production of SO2 from OCS provided by ACHEM
!  --------------------------------------------
   if (gcSU%using_ACHEM_pSO2_OCS) then
       call MAPL_GetPointer(impChem, ACHEM_PSO2_OCS, 'pSO2_OCS', __RC__)
   end if 


!  Unlike GEOS-4 hghte is defined for km+1
!  ---------------------------------------
   hsurf => hghte(i1:i2,j1:j2,km) ! in GEOS-5 hghte is in [0,km]

#ifdef DEBUG

   call pmaxmin('SU: oro        ', oro     , qmin, qmax, ijl,1, 1. )
   call pmaxmin('SU: u10m       ', u10m    , qmin, qmax, ijl,1, 1. )
   call pmaxmin('SU: v10m       ', v10m    , qmin, qmax, ijl,1, 1. )
   call pmaxmin('SU: ustar      ', ustar   , qmin, qmax, ijl,1, 1. )
   call pmaxmin('SU: precc      ', precc   , qmin, qmax, ijl,1, 1. )
   call pmaxmin('SU: precl      ', precl   , qmin, qmax, ijl,1, 1. )
   call pmaxmin('SU: pblh       ', pblh    , qmin, qmax, ijl,1, 1. )
   call pmaxmin('SU: shfflux    ', shflux  , qmin, qmax, ijl,1, 1. )
   call pmaxmin('SU: hsurf      ', hsurf   , qmin, qmax, ijl,1, 1. )

   call pmaxmin('SU: tmpu       ', tmpu    , qmin, qmax, ijkl,1, 1. )
   call pmaxmin('SU: rhoa       ', rhoa    , qmin, qmax, ijkl,1, 1. )
   call pmaxmin('SU: u          ', u       , qmin, qmax, ijkl,1, 1. )
   call pmaxmin('SU: v          ', v       , qmin, qmax, ijkl,1, 1. )
   call pmaxmin('SU: hghte      ', hghte   , qmin, qmax, ijkl,1, 1. )
   call pmaxmin('SU: pfllsan    ', pfllsan , qmin, qmax, ijk1l,1, 1. )
   call pmaxmin('SU: pfilsan    ', pfilsan , qmin, qmax, ijk1l,1, 1. )

#endif

#ifdef DEBUG

   call pmaxmin('SU: h2o2', gcSU%h2o2_int(i1:i2,j1:j2,1:km), qmin, qmax, &
                 ijl, km, 1. )
   call pmaxmin('SU: oh', gcSU%oh_conc(i1:i2,j1:j2,1:km), qmin, qmax, &
                 ijl, km, 1. )
   call pmaxmin('SU: no3', gcSU%no3_mr(i1:i2,j1:j2,1:km), qmin, qmax, &
                 ijl, km, 1. )

#endif


   call SulfateUpdateOxidants ( impChem, iNAME, i1, i2, im, j1, j2, jm, km, cdt, &
                                gcSU%using_GMI_OH, gcSU%using_GMI_NO3, &
                                gcSU%using_GMI_H2O2, &
                                GMI_OHmr, GMI_NO3mr, GMI_H2O2mr, &
                                nymd, nhms, &
                                w_c%grid_esmf, w_c%grid%lon, w_c%grid%lat, &
                                rhoa, &
                                gcSU%nymd_oxidants, &
                                gcSU%oh_conc, gcSU%no3_mr, gcSU%h2o2_mr, &
                                xoh, xno3, xh2o2, gcSU%recycle_H2O2 )

#ifdef DEBUG
   CALL pmaxmin('SU: OH_conc', gcSU%oh_conc, qmin, qmax, ijl,km, 1. )
   CALL pmaxmin('SU: NO3_mr ', gcSU%no3_mr, qmin, qmax, ijl,km, 1. )
   CALL pmaxmin('SU: H2O2_mr', gcSU%h2o2_mr, qmin, qmax, ijl,km, 1. )
   CALL pmaxmin('SU: OH     ', xoh, qmin, qmax, ijl,km, 1. )
   CALL pmaxmin('SU: NO3    ', xno3, qmin, qmax, ijl,km, 1. )
   CALL pmaxmin('SU: H2O2   ', xh2o2, qmin, qmax, ijl,km, 1. )
#endif

!  Settling calculation
!  Sulfate particle radius [m] and density [kg m-3]
!  ---------------------------------------------
   allocate( SU_radius(nbins), SU_rhop(nbins) )
   SU_radius = 1.e-6*gcSU%radius
   SU_rhop   = gcSU%rhop

RUN_ALARM: if (gcSU%run_alarm) then

   call Chem_Settling ( i1, i2, j1, j2, km, n1, n2, nbins, gcSU%rhFlag, &
                          SU_radius, SU_rhop, cdt, w_c, tmpu, rhoa, hsurf,    &
                          hghte, SU_set, rc )

!  If doing the ACHEM provided pSO2 from OCS then add to SO2 here
!  --------------------------------------------------------------
   IF(gcSU%using_ACHEM_pSO2_OCS .and. associated(ACHEM_PSO2_OCS) ) THEN
      w_c%qa(n1+gcSU%nSO2-1)%data3d = &
                    w_c%qa(n1+gcSU%nSO2-1)%data3d + ACHEM_PSO2_OCS*cdt
   ENDIF

!  SU Chemistry Driver (dry deposition and chemistry)
!  -----------
   call SU_ChemDrv ( i1, i2, j1, j2, km, nbins, cdt, nymd, nhms, gcSU, w_c, &
                     ustar, u, v, shflux, oro, pblh, tmpu, cloud, rhoa, hghte, &
                     SU_dep, SU_PSO2, SU_PMSA, SU_pSO4, SU_PSO4g, SU_PSO4aq, & ! 2d diagnostics
                     pso2, pmsa, pso4, pso4g, pso4aq,  &                       ! 3d diagnostics
                     xoh, xno3, xh2o2, &                                       ! oxidants
                     rc)

!  Sulfate Large-scale Wet Removal
!  -------------------------------
   KIN = .TRUE.
   call SU_Wet_Removal ( i1, i2, j1, j2, km, nbins, cdt, 'sulfur', KIN, &
                         ple, rhoa, gcSU, w_c, &
                         precc, precl, pfllsan, pfilsan, &
                         tmpu, SU_wet, SU_pso4, SU_pso4wet, pso4, pso4wet, rc )

!  Sulfate Convective-scale Mixing and Wet Removal
!  -----------------------------------------------
   KIN = .TRUE.
   icdt = cdt
   allocate(cmfmc_(i1:i2,j1:j2,km+1), qccu_(i1:i2,j1:j2,km), &
            dtrain_(i1:i2,j1:j2,km), airmass_(i1:i2,j1:j2,km), &
            delz_(i1:i2,j1:j2,km), vud_(i1:i2,j1:j2,km), &
            tc_(i1:i2,j1:j2,km,n1:n2), delp_(i1:i2,j1:j2,km), &
            airmol_(i1:i2,j1:j2,km), tmpu_(i1:i2,j1:j2,km), &
            bcnv_(i1:i2,j1:j2,n1:n2), ple_(i1:i2,j1:j2,km+1), &
            area_(i1:i2,j1:j2), frlake_(i1:i2,j1:j2), &
            frocean_(i1:i2,j1:j2), frseaice_(i1:i2,j1:j2),&
            h2o2_(i1:i2,j1:j2,km), __STAT__ )

   bcnv_            = 0.0
   area_            = area
   frlake_          = frlake
   frocean_         = frocean
   frseaice_        = frseaice 
   do k = 1, km+1
    cmfmc_(:,:,k)   = cmfmc(:,:,km-k+1)
    ple_(:,:,k)     = ple(:,:,km-k+1)
   end do
   do k = 1, km
    dtrain_(:,:,k)  = dtrain(:,:,km-k+1)
    qccu_(:,:,k)    = qlcn(:,:,km-k+1) + qicn(:,:,km-k+1)
    delp_(:,:,k)    = w_c%delp(:,:,km-k+1)/100.
    airmass_(:,:,k) = w_c%delp(:,:,km-k+1)/grav*area_
    airmol_(:,:,k)  = airmass_(:,:,k)*1000./28.966
    delz_(:,:,k)    = w_c%delp(:,:,km-k+1)/grav/rhoa(:,:,km-k+1)
    tmpu_(:,:,k)    = tmpu(:,:,km-k+1)
   enddo
   do n = n1, n2
    do k = 1, km
     tc_(:,:,k,n)   = w_c%qa(n)%data3d(:,:,km-k+1)
    enddo
   enddo
!  H2O2 is in vmr and SU are mmr.  Convert H2O2 to mmr
   do k = 1, km
     h2o2_(:,:,k)                = gcSU%h2o2_int(:,:,km-k+1)*rH2O2
   enddo
   call set_vud(i1, i2, j1, j2, km, frlake_, frocean_, frseaice_, cmfmc_, qccu_, &
                airmass_, delz_, area_, vud_)
   call convection(i1, i2, j1, j2, km, n1, n2, icdt, 'sulfur', kin, &
                   tc_, cmfmc_, dtrain_, area_, delz_, delp_, vud_, &
                   airmass_, airmol_, tmpu_, ple_, &
                   bcnv_, h2o2_)

!  Return adjusted tracer to mixing ratio
   do n = n1, n2
    do k = 1, km
     w_c%qa(n)%data3d(:,:,km-k+1) = tc_(:,:,k,n)
    enddo
   enddo
!  Return adjusted h2o2
   do k = 1, km
     gcSU%h2o2_int(:,:,km-k+1) = h2o2_(:,:,k)/rH2O2
   enddo

!  Note GOCART returns bcnv_ as negative, recast for my diagnostic
   if(associated(SU_conv(1)%data2d)) SU_conv(1)%data2d = 0.0
   if(associated(SU_conv(2)%data2d)) SU_conv(2)%data2d = -bcnv_(:,:,n1+gcSU%nSO2-1)/area_/icdt
   if(associated(SU_conv(3)%data2d)) SU_conv(3)%data2d = -bcnv_(:,:,n1+gcSU%nSO4-1)/area_/icdt
   if(associated(SU_conv(4)%data2d)) SU_conv(4)%data2d = -bcnv_(:,:,n1+gcSU%nMSA-1)/area_/icdt

   deallocate(cmfmc_, qccu_, dtrain_, tc_, airmass_, &
              delz_, vud_, delp_, airmol_, tmpu_, bcnv_, ple_, &
              area_, frlake_, frocean_, frseaice_, h2o2_, __STAT__ )

! Update GMI Combo oxidants before exiting.
! Note: H2O2 is the only one modified as of this writing.
! -------------------------------------------------------
   IF(gcSU%using_GMI_H2O2 .AND. gcSU%export_H2O2) &
       GMI_H2O2mr(i1:i2,j1:j2,1:km) = gcSU%h2o2_int(i1:i2,j1:j2,1:km)

   end if RUN_ALARM


!  Compute the desired output diagnostics here
!  Ideally this will go where chemout is called in fvgcm.F since that
!  will reflect the distributions after transport, etc.
!  -----------
   call SU_Compute_Diags(i1, i2, j1, j2, km, nbins, gcSU, w_c, tmpu, rhoa, u, v, &
                         SU_DMSsfcmass, SU_DMScolmass, &
                         SU_MSAsfcmass, SU_MSAcolmass, &
                         SU_SO2sfcmass, SU_SO2colmass, &
                         SU_SO4sfcmass, SU_SO4colmass, &
                         SU_exttau, SU_scatau, SU_SO4mass,  &
                         SU_conc, SU_extcoef, SU_scacoef, &
                         SU_angstrom, SU_fluxu, SU_fluxv, &
                         SU_sarea, SU_snum, rc)


   deallocate(xoh, xno3, xh2o2, SU_radius, SU_rhop, stat=STATUS)
   VERIFY_(STATUS)

   RETURN

CONTAINS

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  SU_ChemDrv - Do SU cycle chemistry following GOCART
!
! !INTERFACE:
!

   subroutine SU_ChemDrv ( i1, i2, j1, j2, km, nbins, cdt, nymd, nhms, gcSU, &
                           w_c, ustar, u, v, shflux, oro, pblh, tmpu, &
                           cloud, rhoa, hghte, &
                           su_dep, &
                           su_pSO2, su_pMSA, su_pSO4, su_pSO4g, su_pSO4aq, &   ! 2d diagnostics
                           pSO2, pMSA, pSO4, pSO4g, pSO4aq,  &                 ! 3d diagnostics
                           xoh, xno3, xh2o2, &
                           rc)

! !USES:

  implicit NONE

! !INPUT PARAMETERS:

   integer, intent(in) :: i1, i2, j1, j2, km, nbins, nymd, nhms
   real, intent(in)    :: cdt
   type(SU_GridComp1), intent(inout)   :: gcSU       ! SU Grid Component
   real, pointer, dimension(:,:,:)     :: tmpu, cloud, rhoa, u, v, hghte
   real, pointer, dimension(:,:)       :: ustar, shflux, oro, pblh
   real, pointer, dimension(:,:,:)     :: xoh, xno3, xh2o2

! !OUTPUT PARAMETERS:

   type(Chem_Bundle), intent(inout) :: w_c            ! Chemical tracer fields
   type(Chem_Array), intent(inout)  :: su_dep(nbins)  ! Mass lost by deposition
                                                      ! to surface, kg/m2/s
!  chemical production terms d(mixing ratio) /s
   type(Chem_Array), intent(inout)  :: su_pSO2, su_pMSA, su_pSO4, su_pSO4g, su_pSO4aq
   type(Chem_Array), intent(inout)  :: pSO2, pMSA, pSO4, pSO4g, pSO4aq 

   integer, intent(out)             :: rc          ! Error return code:
                                                   !  0 - all is well
                                                   !  1 - 
   character(len=*), parameter :: myname = 'SU_ChemDrv'

! !DESCRIPTION: Updates the SU concentration due to chemistry
!  The SU grid component is currently established with 4 different
!  species (bins) following this convection:
!   1) DMS
!   2) SO2
!   3) SO4
!   4) MSA
!  Accordingly we have 4 chemical cycles to follow through, which are
!  sub-subroutines under this one.
!  The chemistry is a function of OH, NO3, and H2O2 concentrations
!  as well as DMS, SO2, SO4, MSA concentrations.  It is also a function
!  of solar zenith angle and temperature.  We pass in temperature.  SZA
!  will be a function of time of day and lat/lon.  For now we simply add
!  this to the grid component before calculating it.  I bet this is
!  somewhere else in the model.
!
! !REVISION HISTORY:
!
!  06Nov2003, Colarco
!
!EOP
!-------------------------------------------------------------------------

! !Local Variables
   integer :: ndystep, i, j, k, im
   real :: pSO2_DMS(i1:i2,j1:j2,1:km), pMSA_DMS(i1:i2,j1:j2,1:km), &
           pSO4g_SO2(i1:i2,j1:j2,1:km), pSO4aq_SO2(i1:i2,j1:j2,1:km)

!  Variables used in chemistry step
   real :: drydepf(i1:i2,j1:j2)
   real :: qmin, qmax
   integer :: ijl, ijkl, n1, STATUS
   integer :: nDMS, nSO2, nSO4, nMSA

   ijl  = ( i2 - i1 + 1 ) * ( j2 - j1 + 1 )
   ijkl = ijl * km
   n1 = w_c%reg%i_SU

   nDMS = gcSU%nDMS
   nSO2 = gcSU%nSO2
   nSO4 = gcSU%nSO4
   nMSA = gcSU%nMSA

!  Reset the production terms
   pSO2_DMS(i1:i2,j1:j2,1:km) = 0.
   pMSA_DMS(i1:i2,j1:j2,1:km) = 0.
   pSO4g_SO2(i1:i2,j1:j2,1:km) = 0.
   pSO4aq_SO2(i1:i2,j1:j2,1:km) = 0.
   if( associated(su_pSO2%data2d) )   su_pSO2%data2d(i1:i2,j1:j2) = 0.
   if( associated(su_pMSA%data2d) )   su_pMSA%data2d(i1:i2,j1:j2) = 0.
   if( associated(su_pSO4%data2d) )   su_pSO4%data2d(i1:i2,j1:j2) = 0.
   if( associated(su_pSO4g%data2d) )  su_pSO4g%data2d(i1:i2,j1:j2) = 0.
   if( associated(su_pSO4aq%data2d) ) su_pSO4aq%data2d(i1:i2,j1:j2) = 0.
   if( associated(pSO2%data3d) )      pSO2%data3d(i1:i2,j1:j2,1:km) = 0.
   if( associated(pMSA%data3d) )      pMSA%data3d(i1:i2,j1:j2,1:km) = 0.
   if( associated(pSO4%data3d) )      pSO4%data3d(i1:i2,j1:j2,1:km) = 0.
   if( associated(pSO4g%data3d) )     pSO4g%data3d(i1:i2,j1:j2,1:km) = 0.
   if( associated(pSO4aq%data3d) )    pSO4aq%data3d(i1:i2,j1:j2,1:km) = 0.


!  Now call the chemistry packages...
!  ----------------------------------
   call SulfateChemDriverGOCART ( i1, i2, j1, j2, km, n1, &
                                  nbins, cdt, nymd, nhms, &
                                  w_c%grid%lon, w_c%grid%lat, &
                                  w_c%qa(n1+nDMS-1)%data3d, &
                                  w_c%qa(n1+nSO2-1)%data3d, &
                                  w_c%qa(n1+nSO4-1)%data3d, &
                                  w_c%qa(n1+nMSA-1)%data3d, &
                                  xoh, xno3, xh2o2, &
                                  u, v, w_c%delp, tmpu, cloud, rhoa, hghte, &
                                  ustar, shflux, oro, pblh, z0h, &
                                  SU_dep, SU_PSO2, SU_PMSA, &
                                  SU_PSO4, SU_PSO4g, SU_PSO4aq, &     ! 2d diagnostics
                                  pso2, pmsa, pso4, pso4g, pso4aq,  & ! 3d diagnostics
                                  rc)


!  Save the h2o2 value after chemistry
   gcSU%h2o2_int = xh2o2

#ifdef DEBUG
   if(associated(su_pso2%data2d)) call pmaxmin('SU: su_pso2',su_pso2%data2d,qmin,qmax,ijl,1,1.)
   if(associated(su_pmsa%data2d)) call pmaxmin('SU: su_pmsa',su_pmsa%data2d,qmin,qmax,ijl,1,1.)
   if(associated(su_pso4g%data2d)) call pmaxmin('SU: su_pso4g',su_pso4g%data2d,qmin,qmax,ijl,1,1.)
   if(associated(su_pso4aq%data2d)) call pmaxmin('SU: su_pso4aq',su_pso4aq%data2d,qmin,qmax,ijl,1,1.)
   call pmaxmin('SU:  pSO4g_SO2',  pSO4g_SO2, qmin, qmax, ijl, km, 1. )
   call pmaxmin('SU: pSO4aq_SO2', pSO4aq_SO2, qmin, qmax, ijl, km, 1. )
   call pmaxmin('SU:   pSO2_DMS',   pSO2_DMS, qmin, qmax, ijl, km, 1. )
   call pmaxmin('SU:   pMSA_DMS',   pMSA_DMS, qmin, qmax, ijl, km, 1. )
#endif

   rc = 0

   end subroutine SU_ChemDrv

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  SU_Wet_Removal - Removal of dust by precipitation
!  NOTE: For the removal term, fluxout is the sum of the in-cloud
!        convective and large-scale washout and the total flux across
!        the surface due to below-cloud (rainout) convective and
!        large-scale precipitation reaching the surface.  The fluxout
!        is initialized to zero at the beginning and then at each i, j
!        grid point it is added to.
!        See Chin et al. 1996 for some of the logic of this.  SO4 and
!        MSA are scavenged "normally."  DMS is not scavenged at all.
!        SO2 is weakly soluble in water, but some fraction can be
!        removed because of rapid aqueous phase reaction with H2O2.
!        Accordingly, we compare the mixing ratios of H2O2 and SO2 and
!        only scavenge that fraction of SO2 that is less than the
!        H2O2 mixing ratio.  If any of the scavenged SO2 is released
!        by re-evaporation is emerges as SO4
!        
!
! !INTERFACE:
!

   subroutine SU_Wet_Removal ( i1, i2, j1, j2, km, nbins, cdt, aero_type, kin, &
                               ple, rhoa, gcSU, w_c, &
                               precc, precl, pfllsan, pfilsan, tmpu, &
                               fluxout, pSO4_colflux, pSO4wet_colflux, &
                               pso4, pso4wet, rc )

! !USES:

  implicit NONE

! !INPUT PARAMETERS:

   integer, intent(in) :: i1, i2, j1, j2, km, nbins
   real, intent(in)    :: cdt
   character(len=*)    :: aero_type
   logical, intent(in) :: KIN                 ! true for aerosol
   real, pointer, dimension(:,:)   :: precc   ! total convective precip, [mm day-1]
   real, pointer, dimension(:,:)   :: precl   ! total large-scale prec,  [mm day-1]
   real, pointer, dimension(:,:,:) :: pfllsan ! 
   real, pointer, dimension(:,:,:) :: pfilsan ! 
   real, pointer, dimension(:,:,:) :: tmpu    ! temperature, [K]
   real, pointer, dimension(:,:,:) :: rhoa    ! air density, [kg m-3]
   real, pointer, dimension(:,:,:) :: ple     ! level edge air pressure

! !OUTPUT PARAMETERS:

   type(SU_GridComp1), intent(inout) :: gcSU            ! SU Grid Component
   type(Chem_Bundle), intent(inout)  :: w_c             ! Chemical tracer fields
   type(Chem_Array), intent(inout)   :: fluxout(nbins)  ! Mass lost by wet dep
                                                        ! to surface, kg/m2/s
   type(Chem_Array), intent(inout)   :: pSO4_colflux    ! total chemical
                                                        ! production of SO4
                                                        ! from SO2 
                                                        ! (column integrated)
   type(Chem_Array), intent(inout)   :: pSO4wet_colflux ! aqueous chemical
                                                        ! production of SO4
                                                        ! from SO2 
                                                        ! (column integrated)
   type(Chem_Array), intent(inout)   :: pSO4            ! total chemical 
                                                        ! production of SO4 from SO2
   type(Chem_Array), intent(inout)   :: pSO4wet         ! aqueous chemical 
                                                        ! production of SO4 from SO2
   integer, intent(out)              :: rc              ! Error return code:
                                                        !  0 - all is well
                                                        !  1 - 
   character(len=*), parameter :: myname = 'SU_Wet_Removal'

! !DESCRIPTION: Updates the dust concentration in each vertical layer
!               due to wet removal
!
! !REVISION HISTORY:
!
!  17Nov2003, Colarco
!
!EOP
!-------------------------------------------------------------------------

! !Local Variables
   integer  ::  i, j, k, iit, n, LH, kk, ios
   integer  ::  n1, n2
   real :: pdog(i1:i2,j1:j2,km)      ! air mass factor dp/g [kg m-2]
   real*8 :: Td_ls, Td_cv              ! ls and cv timescales [s]
   real*8 :: pls, pcv, pac             ! ls, cv, tot precip [mm day-1]
   real*8 :: qls(km), qcv(km)          ! ls, cv portion of moisture tendency [kg m-3 s-1]
   real*8 :: qmx, qd, A                ! temporary variables on moisture
   real*8 :: F, B, BT                  ! temporary variables on cloud, freq.
   real*8, allocatable :: fd(:,:)      ! flux across layers [kg m-2]
   real*8, allocatable :: dpfli(:,:,:) ! 
   real*8, allocatable :: DC(:)        ! scavenge change in mass mixing ratio
   real :: c_h2o(i1:i2,j1:j2,km), cldliq(i1:i2,j1:j2,km), cldice(i1:i2,j1:j2,km)
   real, parameter :: kb = 1.3807e-23 ! Boltzmann constant [kg m2 s-1 K-1 mol-1]
   real, parameter :: m_air = 4.8096e-26 ! Mass of <avg> air molecule [kg]

!  Rain parameters (from where?)
   real, parameter :: B0_ls = 1.0e-4
   real, parameter :: F0_ls = 1.0
   real, parameter :: XL_ls = 5.0e-4
   real, parameter :: B0_cv = 1.5e-3
   real, parameter :: F0_cv = 0.3
   real, parameter :: XL_cv = 2.0e-3
   real, parameter :: one = 1.0, zero = 0.0

!  Integer locations of SO2, etc. species
   integer :: nDMS, nSO2, nSO4, nMSA

!  Conversion of SO2 mmr to SO2 vmr (since H2O2 is carried around like
!  a volume mixing ratio)
   real*8 :: fmr, SO2Soluble
   fMR = airMolWght / gcSU%fMassSO2

   rc=0

!  Initialize local variables
!  --------------------------
!  c_h2o, cldliq, and cldice are respectively intended to be the 
!  water mixing ratio (liquid or vapor?, in or out of cloud?)
!  cloud liquid water mixing ratio
!  cloud ice water mixing ratio
   c_h2o  = (10d0**(-2663.5d0/tmpu(:,:,:) + 12.537d0 ) ) /  &
                   (ple(:,:,0:km-1)+ple(:,:,1:km)) /2d0   
   cldliq = 0.d0
   where(tmpu >  248.) cldliq = 1.d-6 * ( ( tmpu - 248.d0) / 20.d0 )
   where(tmpu >= 268.) cldliq = 1.d-6
   cldice = 1.d-6 - cldliq

   do n = 1, nbins
    if( associated(fluxout(n)%data2d) ) fluxout(n)%data2d(i1:i2,j1:j2) = 0.0
   end do
   if( associated(pso4wet_colflux%data2d)) pso4wet_colflux%data2d(i1:i2,j1:j2) = 0.
   if( associated(pso4wet%data3d) ) pso4wet%data3d(i1:i2,j1:j2,1:km) = 0.  

   n1  = w_c%reg%i_SU
   n2  = w_c%reg%j_SU
   nDMS = gcSU%nDMS
   nSO2 = gcSU%nSO2
   nSO4 = gcSU%nSO4
   nMSA = gcSU%nMSA

!  Allocate the dynamic arrays
   allocate(fd(km,nbins),stat=ios)
   if(ios .ne. 0) stop
   allocate(dc(nbins),stat=ios)
   if(ios .ne. 0) stop
   allocate(dpfli(i1:i2, j1:j2, km),stat=ios)
   if(ios .ne. 0) stop

!  Duration of rain: ls = model timestep, cv = 1800 s (<= cdt)
   Td_ls = cdt
   Td_cv = 1800.

!  Accumulate the 3-dimensional arrays of rhoa and pdog
   pdog = w_c%delp/grav
   dpfli = pfllsan(:,:,1:km)-pfllsan(:,:,0:km-1)+pfilsan(:,:,1:km)-pfilsan(:,:,0:km-1)

!  Loop over spatial indices
   do j = j1, j2
    do i = i1, i2

!    Check for total precipitation amount
!    Assume no precip in column if precl+precc = 0
     pac = precl(i,j) + precc(i,j)
     if(pac .le. 0.) goto 100
     pls = precl(i,j)
     pcv = precc(i,j)

!    Initialize the precipitation fields
     qls(:)  = 0.
     qcv(:)  = 0.
     fd(:,:) = 0.
     Dc(:)   = 0.

!    Find the highest model layer experiencing rainout.  Assumes no
!    scavenging if T < 258 K
     LH = 0
     do k = 1, km
      if(dpfli(i,j,k) .gt. 0. .and. tmpu(i,j,k) .gt. 258.) then
       LH = k
       goto 15
      endif
     end do
 15  continue
     if(LH .lt. 1) goto 100

     do k = LH, km
      qls(k) = dpfli(i,j,k)/pdog(i,j,k)*rhoa(i,j,k)
     end do

!    Loop over vertical to do the scavenging!
     do k = LH, km

!-----------------------------------------------------------------------------
!   (1) LARGE-SCALE RAINOUT:             
!       Tracer loss by rainout = TC0 * F * exp(-B*dt)
!         where B = precipitation frequency,
!               F = fraction of grid box covered by precipitating clouds.
!       We assume that tracer scavenged by rain is falling down to the
!       next level, where a fraction could be re-evaporated to gas phase
!       if Qls is less then 0 in that level.
!-----------------------------------------------------------------------------
      if (qls(k) .gt. 0.) then
       F  = F0_ls / (1. + F0_ls*B0_ls*XL_ls/(qls(k)*cdt/Td_ls))
       B  = B0_ls/F0_ls +1./(F0_ls*XL_ls/qls(k)) 
       BT = B * Td_ls
       if (BT.gt.10.) BT = 10.               !< Avoid overflow >
!      What is the soluble amount of SO2?
       SO2Soluble = min(fmr*w_c%qa(n1+nSO2-1)%data3d(i,j,k),gcSU%h2o2_int(i,j,k)*one)/fmr
       if(SO2Soluble .lt. 0.) SO2Soluble = 0.

!      Adjust SU amounts
       DC(nDMS) = 0.
       DC(nSO2) = SO2Soluble * F * (1.-exp(-BT))
       DC(nSO4) = w_c%qa(n1+nSO4-1)%data3d(i,j,k) * F * (1.-exp(-BT))
       DC(nMSA) = w_c%qa(n1+nMSA-1)%data3d(i,j,k) * F * (1.-exp(-BT))

!      Adjust H2O2 concentration in cloudy portion of cell
       if(fmr*w_c%qa(n1+nSO2-1)%data3d(i,j,k) .gt. gcSU%h2o2_int(i,j,k)) then
!        gcSU%h2o2_int(i,j,k) = max(zero,(1.-F)*gcSU%h2o2_int(i,j,k))
! GOCART removes all
        gcSU%h2o2_int(i,j,k) = 0.
       else
        gcSU%h2o2_int(i,j,k) &
          = gcSU%h2o2_int(i,j,k) - F*fmr*w_c%qa(n1+nSO2-1)%data3d(i,j,k)
       endif

       do n = 1, nbins
        if (DC(n).lt.0.) DC(n) = 0.
        w_c%qa(n1+n-1)%data3d(i,j,k) = w_c%qa(n1+n-1)%data3d(i,j,k)-DC(n)
        if (w_c%qa(n1+n-1)%data3d(i,j,k) .lt. 1.0E-32) w_c%qa(n1+n-1)%data3d(i,j,k) = 1.0E-32
       end do

!      Flux down:  unit is kg m-2
!      Formulated in terms of production in the layer.  In the revaporation step
!      we consider possibly adding flux from above...
       do n = 1, nbins
        Fd(k,n) = DC(n) * pdog(i,j,k)
       end do

      end if                                    ! if Qls > 0  >>>

!-----------------------------------------------------------------------------
! * (2) LARGE-SCALE WASHOUT:
! *     Occurs when rain at this level is less than above.
!-----------------------------------------------------------------------------
      if(k .gt. LH .and. qls(k) .ge. 0.) then
       if(qls(k) .lt. qls(k-1)) then
!       Find a maximum F overhead until the level where Qls<0.
        Qmx   = 0.
        do kk = k-1,LH,-1
         if (Qls(kk).gt.0.) then
          Qmx = max(Qmx,Qls(kk))
         else
          goto 333
         end if
        end do

 333    continue
        F = F0_ls / (1. + F0_ls*B0_ls*XL_ls/(Qmx*cdt/Td_ls))
        if (F.lt.0.01) F = 0.01
!-----------------------------------------------------------------------------
!  The following is to convert Q(k) from kgH2O/m3/sec to mm/sec in order
!  to use the Harvard formula.  Convert back to mixing ratio by multiplying
!  by rhoa.  Multiply by pdog gives kg/m2/s of precip.  Divide by density
!  of water (=1000 kg/m3) gives m/s of precip and multiply by 1000 gives
!  units of mm/s (omit the multiply and divide by 1000).
!-----------------------------------------------------------------------------

        Qd = Qmx /rhoa(i,j,k)*pdog(i,j,k)
        if (Qd.ge.50.) then
         B = 0.
        else
         B = Qd * 0.1
        end if
        BT = B * cdt
        if (BT.gt.10.) BT = 10.

!       Adjust SO2 for H2O2 oxidation
        SO2Soluble = min(fmr*w_c%qa(n1+nSO2-1)%data3d(i,j,k),gcSU%h2o2_int(i,j,k)*one)/fmr
        if(SO2Soluble .lt. 0.) SO2Soluble = 0.

!       Adjust SU amounts
        DC(nDMS) = 0.
        DC(nSO2) = SO2Soluble * F * (1.-exp(-BT))
        DC(nSO4) = w_c%qa(n1+nSO4-1)%data3d(i,j,k) * F * (1.-exp(-BT))
        DC(nMSA) = w_c%qa(n1+nMSA-1)%data3d(i,j,k) * F * (1.-exp(-BT))

!       Adjust H2O2 concentration in cloudy portion of cell
        if(fmr*w_c%qa(n1+nSO2-1)%data3d(i,j,k) .gt. gcSU%h2o2_int(i,j,k)) then
         gcSU%h2o2_int(i,j,k) = max(zero,(one-F)*gcSU%h2o2_int(i,j,k))
!  GOCART removes all
!         gcSU%h2o2_int(i,j,k) = 0.
        else
         gcSU%h2o2_int(i,j,k) &
           = gcSU%h2o2_int(i,j,k) - F*fmr*w_c%qa(n1+nSO2-1)%data3d(i,j,k)
        endif
 
        do n = 1, nbins
         if (DC(n).lt.0.) DC(n) = 0.
         w_c%qa(n1+n-1)%data3d(i,j,k) = w_c%qa(n1+n-1)%data3d(i,j,k)-DC(n)
         if (w_c%qa(n1+n-1)%data3d(i,j,k) .lt. 1.0E-32) w_c%qa(n1+n-1)%data3d(i,j,k) = 1.0E-32
        end do
 
        do n = 1, nbins
         if( associated(fluxout(n)%data2d) ) then
          fluxout(n)%data2d(i,j) = fluxout(n)%data2d(i,j)+DC(n)*pdog(i,j,k)/cdt
         endif
        end do

       end if
      end if                                    ! if ls washout  >>>

!-----------------------------------------------------------------------------
!  (3) CONVECTIVE RAINOUT:
!      Tracer loss by rainout = dd0 * F * exp(-B*dt)
!        where B = precipitation frequency,
!              F = fraction of grid box covered by precipitating clouds.
!-----------------------------------------------------------------------------

      if (qcv(k) .gt. 0.) then
       F  = F0_cv / (1. + F0_cv*B0_cv*XL_cv/(Qcv(k)*cdt/Td_cv))
       B  = B0_cv
       BT = B * Td_cv
       if (BT.gt.10.) BT = 10.               !< Avoid overflow >

!      Adjust SO2 for H2O2 oxidation
       SO2Soluble = min(fmr*w_c%qa(n1+nSO2-1)%data3d(i,j,k),gcSU%h2o2_int(i,j,k)*one)/fmr
       if(SO2Soluble .lt. 0.) SO2Soluble = 0.

!      Adjust SU amounts
       DC(nDMS) = 0.
       DC(nSO2) = SO2Soluble * F * (1.-exp(-BT))
       DC(nSO4) = w_c%qa(n1+nSO4-1)%data3d(i,j,k) * F * (1.-exp(-BT))
       DC(nMSA) = w_c%qa(n1+nMSA-1)%data3d(i,j,k) * F * (1.-exp(-BT))
       DC(nSO4) = 0.
       DC(nMSA) = 0.

!      Adjust H2O2 concentration in cloudy portion of cell
       if(fmr*w_c%qa(n1+nSO2-1)%data3d(i,j,k) .gt. gcSU%h2o2_int(i,j,k)) then
        gcSU%h2o2_int(i,j,k) = max(zero,(one-F)*gcSU%h2o2_int(i,j,k))
       else
        gcSU%h2o2_int(i,j,k) &
          = gcSU%h2o2_int(i,j,k) - F*fmr*w_c%qa(n1+nSO2-1)%data3d(i,j,k)
       endif

       do n = 1, nbins
        if (DC(n).lt.0.) DC(n) = 0.
        w_c%qa(n1+n-1)%data3d(i,j,k) = w_c%qa(n1+n-1)%data3d(i,j,k)-DC(n)
        if (w_c%qa(n1+n-1)%data3d(i,j,k) .lt. 1.0E-32) w_c%qa(n1+n-1)%data3d(i,j,k) = 1.0E-32
       end do

!      Flux down:  unit is kg m-2
!      Formulated in terms of production in the layer.  In the revaporation step
!      we consider possibly adding flux from above...
       do n = 1, nbins
        Fd(k,n) = Fd(k,n) + DC(n)*pdog(i,j,k)
       end do

      end if                                  ! if Qcv > 0   >>>

!-----------------------------------------------------------------------------
!  (4) CONVECTIVE WASHOUT:
!      Occurs when rain at this level is less than above.
!-----------------------------------------------------------------------------

      if (k.gt.LH .and. Qcv(k).ge.0.) then
       if (Qcv(k).lt.Qcv(k-1)) then
!-----  Find a maximum F overhead until the level where Qls<0.
        Qmx   = 0.
        do kk = k-1, LH, -1
         if (Qcv(kk).gt.0.) then
          Qmx = max(Qmx,Qcv(kk))
         else
          goto 444
         end if
        end do

 444    continue
        F = F0_cv / (1. + F0_cv*B0_cv*XL_cv/(Qmx*cdt/Td_cv))
        if (F.lt.0.01) F = 0.01
!-----------------------------------------------------------------------------
!  The following is to convert Q(k) from kgH2O/m3/sec to mm/sec in order
!  to use the Harvard formula.  Convert back to mixing ratio by multiplying
!  by rhoa.  Multiply by pdog gives kg/m2/s of precip.  Divide by density
!  of water (=1000 kg/m3) gives m/s of precip and multiply by 1000 gives
!  units of mm/s (omit the multiply and divide by 1000).
!-----------------------------------------------------------------------------

        Qd = Qmx / rhoa(i,j,k)*pdog(i,j,k)
        if (Qd.ge.50.) then
         B = 0.
        else
         B = Qd * 0.1
        end if
        BT = B * cdt
        if (BT.gt.10.) BT = 10.

!       Adjust SO2 for H2O2 oxidation
        SO2Soluble = min(fmr*w_c%qa(n1+nSO2-1)%data3d(i,j,k),gcSU%h2o2_int(i,j,k)*one)/fmr
        if(SO2Soluble .lt. 0.) SO2Soluble = 0.

!       Adjust SU amounts
        DC(nDMS) = 0.
        DC(nSO2) = SO2Soluble * F * (1.-exp(-BT))
! Sulfate scavenged in moist
!        DC(nSO4) = w_c%qa(n1+nSO4-1)%data3d(i,j,k) * F * (1.-exp(-BT))
!        DC(nMSA) = w_c%qa(n1+nMSA-1)%data3d(i,j,k) * F * (1.-exp(-BT))
        DC(nSO4) = 0.
        DC(nMSA) = 0.

!       Adjust H2O2 concentration in cloudy portion of cell
        if(fmr*w_c%qa(n1+nSO2-1)%data3d(i,j,k) .gt. gcSU%h2o2_int(i,j,k)) then
         gcSU%h2o2_int(i,j,k) = max(zero,(one-F)*gcSU%h2o2_int(i,j,k))
        else
         gcSU%h2o2_int(i,j,k) &
           = gcSU%h2o2_int(i,j,k) - F*fmr*w_c%qa(n1+nSO2-1)%data3d(i,j,k)
        endif
 
        do n = 1, nbins
         if (DC(n).lt.0.) DC(n) = 0.
         w_c%qa(n1+n-1)%data3d(i,j,k) = w_c%qa(n1+n-1)%data3d(i,j,k)-DC(n)
         if (w_c%qa(n1+n-1)%data3d(i,j,k) .lt. 1.0E-32) w_c%qa(n1+n-1)%data3d(i,j,k) = 1.0E-32
        end do

        do n = 1, nbins
         if( associated(fluxout(n)%data2d) ) then
          fluxout(n)%data2d(i,j) = fluxout(n)%data2d(i,j)+DC(n)*pdog(i,j,k)/cdt
         endif
        end do

       end if
      end if                                    ! if cv washout  >>>

!-----------------------------------------------------------------------------
!  (5) RE-EVAPORATION.  Assume that SO2 is re-evaporated as SO4 since it
!      has been oxidized by H2O2 at the level above. 
!-----------------------------------------------------------------------------
!     Add in the flux from above, which will be subtracted if reevaporation occurs
      if(k .gt. LH) then
       do n = 1, nbins
        Fd(k,n) = Fd(k,n) + Fd(k-1,n)
       end do

!      Is there evaporation in the currect layer?
       if (dpfli(i,j,k) .lt. 0.) then
!       Fraction evaporated = H2O(k)evap / H2O(next condensation level).
        if (dpfli(i,j,k-1) .gt. 0.) then

          A =  abs(  dpfli(i,j,k) /  dpfli(i,j,k-1)  )
          if (A .gt. 1.) A = 1.

!         Adjust tracer in the level
!         For the SO2 tracer we do not allow re-evaporation.
!         We compute DC(nSO2) solely to add this to DC(nSO4) and to remove
!         from Fd(k,nSO2)
!         Instead, the SO2 gets re-evaporated to the SO4 bin because of
!         previous H2O2 oxidation

          DC(nDMS) = 0.
          DC(nSO2) = Fd(k-1,nSO2) / pdog(i,j,k) * A
          DC(nSO4) = Fd(k-1,nSO4) / pdog(i,j,k) * A
          DC(nMSA) = Fd(k-1,nMSA) / pdog(i,j,k) * A
          do n = 1, nbins
           if (DC(n).lt.0.) DC(n) = 0.
          end do

          w_c%qa(n1+nMSA-1)%data3d(i,j,k) = w_c%qa(n1+nMSA-1)%data3d(i,j,k) + DC(nMSA)

!         SO2 gets added to SO4, but remember to remove the SO2 from FD!
          w_c%qa(n1+nSO4-1)%data3d(i,j,k) =  w_c%qa(n1+nSO4-1)%data3d(i,j,k) + DC(nSO4) &
                                    + DC(nSO2)*gcSU%fMassSO4/gcSU%fMassSO2
          if( associated(pso4wet_colflux%data2d)) &
             pso4wet_colflux%data2d(i,j) = pso4wet_colflux%data2d(i,j) &
              + DC(nSO2)*gcSU%fMassSO4/gcSU%fMassSO2 / cdt * w_c%delp(i,j,k)/grav
          if( associated(pso4wet%data3d) ) &
             pso4wet%data3d(i,j,k) = DC(nSO2)*gcSU%fMassSO4/gcSU%fMassSO2 / cdt

          if( associated(pso4_colflux%data2d)) &
             pso4_colflux%data2d(i,j) = pso4_colflux%data2d(i,j) &
              + DC(nSO2)*gcSU%fMassSO4/gcSU%fMassSO2 / cdt * w_c%delp(i,j,k)/grav
          if( associated(pso4%data3d) ) &
             pso4%data3d(i,j,k) = pso4%data3d(i,j,k) + DC(nSO2)*gcSU%fMassSO4/gcSU%fMassSO2 / cdt

!         Adjust the flux out of the bottom of the layer--remove SO2 here!
          do n = 1, nbins
           w_c%qa(n1+n-1)%data3d(i,j,k) = &
             max(w_c%qa(n1+n-1)%data3d(i,j,k),tiny(1.0))
           Fd(k,n)  = Fd(k,n) - DC(n)*pdog(i,j,k)
          end do

        endif
       endif                                   ! if -moistq < 0
      endif
     end do  ! k

     do n = 1, nbins
      if( associated(fluxout(n)%data2d) ) then
       fluxout(n)%data2d(i,j) = fluxout(n)%data2d(i,j)+Fd(km,n)/cdt
      endif
     end do

 100 continue
    end do   ! i
   end do    ! j

   deallocate(fd,DC,dpfli,stat=ios)

   end subroutine SU_Wet_Removal


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  SU_Compute_Diags - Calculate dust 2D diagnostics
!
! !INTERFACE:
!

   subroutine SU_Compute_Diags ( i1, i2, j1, j2, km, nbins, gcSU, w_c, tmpu, rhoa, u, v, &
                                 dmssfcmass, dmscolmass, &
                                 msasfcmass, msacolmass, &
                                 so2sfcmass, so2colmass, &
                                 so4sfcmass, so4colmass, &
                                 exttau, scatau, so4mass, so4conc, extcoef, &
                                 scacoef, angstrom, fluxu, fluxv, sarea, snum, rc )

! !USES:

  implicit NONE

! !INPUT PARAMETERS:
   integer, intent(in) :: i1, i2, j1, j2, km, nbins
   type(SU_GridComp1), intent(inout):: gcSU     ! SU Grid Component
   type(Chem_Bundle), intent(in)   :: w_c
   real, pointer, dimension(:,:,:) :: tmpu      ! temperature [K]
   real, pointer, dimension(:,:,:) :: rhoa      ! air density [kg m-3]
   real, pointer, dimension(:,:,:) :: u         ! east-west wind [m s-1]
   real, pointer, dimension(:,:,:) :: v         ! north-south wind [m s-1]

! !OUTPUT PARAMETERS:
   type(Chem_Array), intent(inout)  :: dmssfcmass ! sfc mass concentration kg/m3
   type(Chem_Array), intent(inout)  :: dmscolmass ! col mass density kg/m2
   type(Chem_Array), intent(inout)  :: msasfcmass ! sfc mass concentration kg/m3
   type(Chem_Array), intent(inout)  :: msacolmass ! col mass density kg/m2
   type(Chem_Array), intent(inout)  :: so2sfcmass ! sfc mass concentration kg/m3
   type(Chem_Array), intent(inout)  :: so2colmass ! col mass density kg/m2
   type(Chem_Array), intent(inout)  :: so4sfcmass ! sfc mass concentration kg/m3
   type(Chem_Array), intent(inout)  :: so4colmass ! col mass density kg/m2
   type(Chem_Array), intent(inout)  :: exttau     ! ext. AOT at 550 nm
   type(Chem_Array), intent(inout)  :: scatau     ! sct. AOT at 550 nm
   type(Chem_Array), intent(inout)  :: so4mass    ! 3D sulfate mass mr
   type(Chem_Array), intent(inout)  :: so4conc    ! 3D mass concentration, kg/m3
   type(Chem_Array), intent(inout)  :: extcoef    ! 3D ext. coefficient, 1/m
   type(Chem_Array), intent(inout)  :: scacoef    ! 3D scat.coefficient, 1/m
   type(Chem_Array), intent(inout)  :: angstrom   ! 470-870 nm Angstrom parameter
   type(Chem_Array), intent(inout)  :: fluxu      ! Column mass flux in x direction
   type(Chem_Array), intent(inout)  :: fluxv      ! Column mass flux in y direction
   type(Chem_Array), intent(inout)  :: sarea      ! Sulfate surface area density [m2 m-3]
   type(Chem_Array), intent(inout)  :: snum       ! Sulfate number density [# m-2]
   integer, intent(out)             :: rc         ! Error return code:
                                                  !  0 - all is well
                                                  !  1 - 

! !DESCRIPTION: Calculates some simple 2d diagnostics from the SU fields
!  NOTE: For now this operates solely on the sulfate bin!!!!
!               Surface concentration (dry)
!               Column mass load (dry)
!               Extinction aot 550 (wet)
!               Scattering aot 550 (wet)
!               For the moment, this is hardwired.
!
! !REVISION HISTORY:
!
!  16APR2004, Colarco
!
!EOP
!-------------------------------------------------------------------------

! !Local Variables
   character(len=*), parameter :: myname = 'SU_Compute_Diags'
   integer :: i, j, k, n, n1, n2, nSO4, nSO2, nDMS, nMSA, ios, nch, idx
   real :: tau, ssa
   character(len=255) :: qname
   real, dimension(i1:i2,j1:j2) :: tau470, tau870
   real    :: ilam550, ilam470, ilam870
   logical :: do_angstrom
   real    :: rh, gf, rwet, rmed, sigma, svol


!  Initialize local variables
!  --------------------------
   n1  = w_c%reg%i_SU
   n2  = w_c%reg%j_SU
   nch   = gcSU%mie_tables%nch
   nSO4  = gcSU%nSO4
   nSO2  = gcSU%nSO2
   nDMS  = gcSU%nDMS
   nMSA  = gcSU%nMSA

!  Get the wavelength indices
!  --------------------------
!  Must provide ilam550 for AOT calculation
   ilam550 = 1.
   ilam470 = 0.
   ilam870 = 0.
   if(nch .gt. 1) then
    do i = 1, nch
     if ( gcSU%mie_tables%channels(i) .ge. 5.49e-7 .and. &
          gcSU%mie_tables%channels(i) .le. 5.51e-7) ilam550 = i
     if ( gcSU%mie_tables%channels(i) .ge. 4.69e-7 .and. &
          gcSU%mie_tables%channels(i) .le. 4.71e-7) ilam470 = i
     if ( gcSU%mie_tables%channels(i) .ge. 8.69e-7 .and. &
          gcSU%mie_tables%channels(i) .le. 8.71e-7) ilam870 = i
    enddo
   endif

!  Determine if going to do Angstrom parameter calculation
!  -------------------------------------------------------
   do_angstrom = .false.
!  If both 470 and 870 channels provided (and not the same) then
!  possibly will do Angstrom parameter calculation
   if(ilam470 .ne. 0. .and. &
      ilam870 .ne. 0. .and. &
      ilam470 .ne. ilam870) do_angstrom = .true.


!  Calculate the diagnostic variables if requested
!  -----------------------------------------------

!  Calculate the surface mass concentration
   if( associated(so4sfcmass%data2d) ) then
      so4sfcmass%data2d(i1:i2,j1:j2) = 0.
      so4sfcmass%data2d(i1:i2,j1:j2) &
       =   w_c%qa(nSO4+n1-1)%data3d(i1:i2,j1:j2,km)*rhoa(i1:i2,j1:j2,km)
   endif
   if( associated(so2sfcmass%data2d) ) then
      so2sfcmass%data2d(i1:i2,j1:j2) = 0.
      so2sfcmass%data2d(i1:i2,j1:j2) &
       =   w_c%qa(nSO2+n1-1)%data3d(i1:i2,j1:j2,km)*rhoa(i1:i2,j1:j2,km)
   endif
   if( associated(dmssfcmass%data2d) ) then
      dmssfcmass%data2d(i1:i2,j1:j2) = 0.
      dmssfcmass%data2d(i1:i2,j1:j2) &
       =   w_c%qa(nDMS+n1-1)%data3d(i1:i2,j1:j2,km)*rhoa(i1:i2,j1:j2,km)
   endif
   if( associated(msasfcmass%data2d) ) then
      msasfcmass%data2d(i1:i2,j1:j2) = 0.
      msasfcmass%data2d(i1:i2,j1:j2) &
       =   w_c%qa(nMSA+n1-1)%data3d(i1:i2,j1:j2,km)*rhoa(i1:i2,j1:j2,km)
   endif


!  Initialize the diagnostic variables
!  -----------------------------------

!  Calculate the column loading
   if( associated(so4colmass%data2d) ) then
      so4colmass%data2d(i1:i2,j1:j2) = 0.
      do k = 1, km
       so4colmass%data2d(i1:i2,j1:j2) &
        =   so4colmass%data2d(i1:i2,j1:j2) &
          + w_c%qa(nSO4+n1-1)%data3d(i1:i2,j1:j2,k)*w_c%delp(i1:i2,j1:j2,k)/grav
      enddo
   endif
   if( associated(so2colmass%data2d) ) then
      so2colmass%data2d(i1:i2,j1:j2) = 0.
      do k = 1, km
       so2colmass%data2d(i1:i2,j1:j2) &
        =   so2colmass%data2d(i1:i2,j1:j2) &
          + w_c%qa(nSO2+n1-1)%data3d(i1:i2,j1:j2,k)*w_c%delp(i1:i2,j1:j2,k)/grav
      enddo
   endif
   if( associated(dmscolmass%data2d) ) then
      dmscolmass%data2d(i1:i2,j1:j2) = 0.
      do k = 1, km
       dmscolmass%data2d(i1:i2,j1:j2) &
        =   dmscolmass%data2d(i1:i2,j1:j2) &
          + w_c%qa(nDMS+n1-1)%data3d(i1:i2,j1:j2,k)*w_c%delp(i1:i2,j1:j2,k)/grav
      enddo
   endif
   if( associated(msacolmass%data2d) ) then
      msacolmass%data2d(i1:i2,j1:j2) = 0.
      do k = 1, km
       msacolmass%data2d(i1:i2,j1:j2) &
        =   msacolmass%data2d(i1:i2,j1:j2) &
          + w_c%qa(nMSA+n1-1)%data3d(i1:i2,j1:j2,k)*w_c%delp(i1:i2,j1:j2,k)/grav
      enddo
   endif


!  Calculate the mass concentration of sulfate 
   if( associated(so4conc%data3d) ) then
      so4conc%data3d(i1:i2,j1:j2,1:km) = 0.
      so4conc%data3d(i1:i2,j1:j2,1:km) = w_c%qa(nSO4+n1-1)%data3d(i1:i2,j1:j2,1:km)*rhoa(i1:i2,j1:j2,1:km)
   endif

!  Mass mixing ratio of sulfate
   if( associated(so4mass%data3d) ) then
      so4mass%data3d(i1:i2,j1:j2,1:km) = 0.
      so4mass%data3d(i1:i2,j1:j2,1:km) = w_c%qa(nSO4+n1-1)%data3d(i1:i2,j1:j2,1:km)
   endif

!  Calculate the column mass flux in x direction
   if( associated(fluxu%data2d) ) then
      fluxu%data2d(i1:i2,j1:j2) = 0.
       do k = 1, km
        fluxu%data2d(i1:i2,j1:j2) &
         =   fluxu%data2d(i1:i2,j1:j2) &
           + w_c%qa(nSO4+n1-1)%data3d(i1:i2,j1:j2,k)*w_c%delp(i1:i2,j1:j2,k)/grav*u(i1:i2,j1:j2,k)
       end do
   endif   
   
!  Calculate the column mass flux in y direction
   if( associated(fluxv%data2d) ) then
      fluxv%data2d(i1:i2,j1:j2) = 0.
       do k = 1, km
        fluxv%data2d(i1:i2,j1:j2) &
         =   fluxv%data2d(i1:i2,j1:j2) &
           + w_c%qa(nSO4+n1-1)%data3d(i1:i2,j1:j2,k)*w_c%delp(i1:i2,j1:j2,k)/grav*v(i1:i2,j1:j2,k)
       end do
   endif      

!  Calculate the extinction and/or scattering AOD
   if( associated(exttau%data2d) .or. associated(scatau%data2d) ) then

      if( associated(exttau%data2d) ) then
       exttau%data2d(i1:i2,j1:j2) = 0.
      endif
      if( associated(scatau%data2d) ) then
       scatau%data2d(i1:i2,j1:j2) = 0.
      endif

      if( associated(extcoef%data3d) ) then
          extcoef%data3d(i1:i2,j1:j2,1:km) = 0.
      endif
      if( associated(scacoef%data3d) ) then
          scacoef%data3d(i1:i2,j1:j2,1:km) = 0.
      endif

!     Note the binning is different for SO4
      do n = nSO4, nSO4

!      Select the name for species
       qname = trim(w_c%reg%vname(w_c%reg%i_SU+n-1))
       idx = Chem_MieQueryIdx(gcSU%mie_tables,qname,rc)
       if(rc .ne. 0) call die(myname, 'cannot find proper Mie table index')

       do k = 1, km
        do j = j1, j2
         do i = i1, i2
          call Chem_MieQuery(gcSU%mie_tables, idx, ilam550, &
              w_c%qa(n1+n-1)%data3d(i,j,k)*w_c%delp(i,j,k)/grav, &
              w_c%rh(i,j,k), tau=tau, ssa=ssa)

!         Calculate the total ext. and scat. coefficients
          if( associated(extcoef%data3d) ) then
              extcoef%data3d(i,j,k) = extcoef%data3d(i,j,k) + &
                                      tau * (grav * rhoa(i,j,k) / w_c%delp(i,j,k))
          endif
          if( associated(scacoef%data3d) ) then
              scacoef%data3d(i,j,k) = scacoef%data3d(i,j,k) + &
                                      ssa * tau * (grav * rhoa(i,j,k) / w_c%delp(i,j,k))
          endif

!         Integrate in the vertical
          if( associated(exttau%data2d) ) then
           exttau%data2d(i,j) = exttau%data2d(i,j) + tau
          endif
          if( associated(scatau%data2d) ) then
           scatau%data2d(i,j) = scatau%data2d(i,j) + tau*ssa
          endif

         enddo
        enddo
       enddo

      enddo  ! nbins

   endif


!  Calculate the 470-870 Angstrom parameter
   if( associated(angstrom%data2d) .and. do_angstrom ) then

      angstrom%data2d(i1:i2,j1:j2) = 0.
!     Set tau to small number by default
      tau470(i1:i2,j1:j2) = tiny(1.0)
      tau870(i1:i2,j1:j2) = tiny(1.0)

      do n = nSO4, nSO4

!      Select the name for species
       qname = trim(w_c%reg%vname(n+n1-1))
       idx = Chem_MieQueryIdx(gcSU%mie_tables,qname,rc)
       if(rc .ne. 0) call die(myname, 'cannot find proper Mie table index')

       do k = 1, km
        do j = j1, j2
         do i = i1, i2

          call Chem_MieQuery(gcSU%mie_tables, idx, ilam470, &
              w_c%qa(n+n1-1)%data3d(i,j,k)*w_c%delp(i,j,k)/grav, &
              w_c%rh(i,j,k), tau=tau)
          tau470(i,j) = tau470(i,j) + tau

          call Chem_MieQuery(gcSU%mie_tables, idx, ilam870, &
              w_c%qa(n+n1-1)%data3d(i,j,k)*w_c%delp(i,j,k)/grav, &
              w_c%rh(i,j,k), tau=tau)
          tau870(i,j) = tau870(i,j) + tau

         enddo
        enddo
       enddo

      enddo  ! nbins
      angstrom%data2d(i1:i2,j1:j2) = &
        -log(tau470(i1:i2,j1:j2)/tau870(i1:i2,j1:j2)) / &
         log(470./870.)
   endif

!  Calculate the sulfate surface area density [m2 m-3], possibly for use in
!  StratChem or other component.  Assumption here is a specified effective
!  radius (gcSU%radius for sulfate) and standard deviation of lognormal
!  distribution.  Hydration is by grid box provided RH and is follows Petters
!  and Kreeidenweis (ACP2007)
   if(associated(sarea%data3d) .or. associated(snum%data3d)) then
        rmed   = w_c%reg%rmed(n1+nSO4-1)                    ! median radius, m
        if(rmed > 0.) then 
         sigma  = w_c%reg%sigma(n1+nSO4-1)                  ! width of lognormal distribution
         do k = 1, km
         do j = j1, j2
          do i = i1, i2
           rh = min(0.95,w_c%rh(i,j,k))
           gf = (1. + 1.19*rh/(1.-rh) )                   ! ratio of wet/dry volume, eq. 5
           rwet = rmed * gf**(1./3.)                      ! wet effective radius, m
!          Wet particle volume m3 m-3
           svol = w_c%qa(n1+nSO4-1)%data3d(i,j,k) * rhoa(i,j,k) / gcSU%rhop(nSO4) * gf
!          Integral of lognormal surface area m2 m-3
           if(associated(sarea%data3d)) sarea%data3d(i,j,k) = 3./rwet*svol*exp(-5./2.*alog(sigma)**2.)
!          Integral of lognormal number density # m-3
           if(associated(snum%data3d)) snum%data3d(i,j,k) = svol / (rwet**3) * exp(-9/2.*alog(sigma)**2.) * 3./4./pi
          enddo
         enddo
        enddo
       endif
   endif   

   rc = 0

   end subroutine SU_Compute_Diags

 end subroutine SU_GridCompRun2_




!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  SU_GridCompFinalize --- The Chem Driver 
!
! !INTERFACE:
!

   subroutine SU_GridCompFinalize1_ ( gcSU, w_c, impChem, expChem, &
                                      nymd, nhms, cdt, rc )

! !USES:

  implicit NONE

! !INPUT/OUTPUT PARAMETERS:

   type(SU_GridComp1), intent(inout) :: gcSU   ! Grid Component

! !INPUT PARAMETERS:

   type(Chem_Bundle), intent(in)  :: w_c      ! Chemical tracer fields   
   integer, intent(in) :: nymd, nhms          ! time
   real,    intent(in) :: cdt                 ! chemical timestep (secs)


! !OUTPUT PARAMETERS:

   type(ESMF_State), intent(inout) :: impChem   ! Import State
   type(ESMF_State), intent(inout) :: expChem   ! Import State
   integer, intent(out) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 -
 
! !DESCRIPTION: This routine finalizes this Grid Component.
!
! !REVISION HISTORY:
!
!  18Sep2003 da Silva  First crack.
!
!EOP
!-------------------------------------------------------------------------

   integer :: ios
   character(len=*), parameter :: myname = 'SU_GridCompFinalize1_'

!  If initialized volcanic emissions from daily tables, clean-up
   if (associated(gcSU%vLat))   deallocate(gcSU%vLat,   stat=ios)
   if (associated(gcSU%vLon))   deallocate(gcSU%vLon,   stat=ios)
   if (associated(gcSU%vSO2))   deallocate(gcSU%vSO2,   stat=ios)
   if (associated(gcSU%vElev))  deallocate(gcSU%vElev,  stat=ios)
   if (associated(gcSU%vCloud)) deallocate(gcSU%vCloud, stat=ios)
   if (associated(gcSU%vStart)) deallocate(gcSU%vStart, stat=ios)
   if (associated(gcSU%vEnd))   deallocate(gcSU%vEnd,   stat=ios)
   rc=0
   return

 end subroutine SU_GridCompFinalize1_

 end module SU_GridCompMod

!-----------------------------------------------------------------------

!                     Single Instance Wrapper

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  SU_SingleInstance_ --- Runs single instance of method
!
! !INTERFACE:
!
  subroutine SU_SingleInstance_ ( Method_, instance, &
                                  gcSU, w_c, impChem, expChem, &
                                  nymd, nhms, cdt, rc )

! !USES:

  Use SU_GridCompMod
  Use ESMF
  Use MAPL
  Use Chem_Mod 

  IMPLICIT NONE

! !INPUT PARAMETERS:

!  Input "function pointer"
!  -----------------------
   interface 
     subroutine Method_ (gc, w, imp, exp, ymd, hms, dt, rcode )
       Use SU_GridCompMod
       Use ESMF
       Use MAPL
       Use Chem_Mod 
       type(SU_GridComp1),  intent(inout)  :: gc
       type(Chem_Bundle),   intent(in)     :: w
       type(ESMF_State),    intent(inout)  :: imp
       type(ESMF_State),    intent(inout)  :: exp
       integer,             intent(in)     :: ymd, hms
       real,                intent(in)     :: dt
       integer,             intent(out)    :: rcode
     end subroutine Method_
   end interface

   integer, intent(in)           :: instance   ! instance number

   TYPE(Chem_Bundle), intent(inout) :: w_c     ! Chemical tracer fields      
   INTEGER, INTENT(IN) :: nymd, nhms           ! time
   REAL,    INTENT(IN) :: cdt                  ! chemical timestep (secs)


! !OUTPUT PARAMETERS:

   TYPE(SU_GridComp1), INTENT(INOUT) :: gcSU    ! Grid Component
   TYPE(ESMF_State), INTENT(INOUT)  :: impChem  ! Import State
   TYPE(ESMF_State), INTENT(INOUT)  :: expChem  ! Export State
   INTEGER, INTENT(OUT) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 - 

! !DESCRIPTION: Finalizes the CO Grid Component. Multiple instance
!               version.
!
! !REVISION HISTORY:
!
!  27Feb2008  da Silva  Introduced multiple instances
!
!EOP
!-------------------------------------------------------------------------

  integer n_SU, i_SU, j_SU
  character(len=255) :: dmsname, so2name, so4name, msaname

! Save overall CO indices
! -----------------------
  n_SU = w_c%reg%n_SU
  i_SU = w_c%reg%i_SU
  j_SU = w_c%reg%j_SU

! Save the name of the variables in this instance
! -----------------------------------------------
  dmsname = trim(w_c%reg%vname(i_SU + 4*(instance - 1)))
  so2name = trim(w_c%reg%vname(i_SU + 4*(instance - 1) + 1))
  so4name = trim(w_c%reg%vname(i_SU + 4*(instance - 1) + 2))
  msaname = trim(w_c%reg%vname(i_SU + 4*(instance - 1) + 3))
  
! Customize indices for this particular instance
! ----------------------------------------------
  w_c%reg%n_SU = 4
  w_c%reg%i_SU = i_SU + 4*(instance - 1)
  w_c%reg%j_SU = i_SU + 4*(instance - 1) + 3

! Update names to "full" version names iff so4name != 'SO4v'
  if(trim(so4name) .ne. 'SO4v') then
   w_c%reg%vname(i_SU + 4*(instance - 1))     = w_c%reg%vname(i_SU)
   w_c%reg%vname(i_SU + 4*(instance - 1) + 1) = w_c%reg%vname(i_SU + 1)
   w_c%reg%vname(i_SU + 4*(instance - 1) + 2) = w_c%reg%vname(i_SU + 2)
   w_c%reg%vname(i_SU + 4*(instance - 1) + 3) = w_c%reg%vname(i_SU + 3)
  endif

  
! Execute the instance method
! ---------------------------
  call Method_ ( gcSU, w_c, impChem, expChem, &
                 nymd, nhms, cdt, rc )

! Restore the overall SU indices
! ------------------------------
  w_c%reg%vname(i_SU + 4*(instance - 1))     = dmsname
  w_c%reg%vname(i_SU + 4*(instance - 1) + 1) = so2name
  w_c%reg%vname(i_SU + 4*(instance - 1) + 2) = so4name
  w_c%reg%vname(i_SU + 4*(instance - 1) + 3) = msaname
  w_c%reg%n_SU = n_SU
  w_c%reg%i_SU = i_SU
  w_c%reg%j_SU = j_SU

  end subroutine SU_SingleInstance_

!-----------------------------------------------------------------------
