#include "MAPL_Generic.h"

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  BRC_GridCompMod --- BRC Grid Component Class
!
! !INTERFACE:
!

   module  BRC_GridCompMod

! !USES:

   USE ESMF
   USE MAPL

   use Chem_Mod              ! Chemistry Base Class
   use Chem_StateMod         ! Chemistry State
   use Chem_ConstMod, only: grav, von_karman, cpd, &
                            undefval => undef         ! Constants !
   use Chem_UtilMod          ! I/O
   use Chem_MieMod           ! Aerosol LU Tables, calculator
   use m_inpak90             ! Resource file management
   use m_die, only: die
   use Chem_SettlingMod      ! Settling
   use DryDepositionMod      ! Dry Deposition
   use WetRemovalMod         ! Large-scale Wet Removal
   use ConvectionMod         ! Offline convective mixing/scavenging

   implicit none

! !PUBLIC TYPES:
!
   PRIVATE
   PUBLIC  BRC_GridComp       ! The BRC object 
   PUBLIC  BRC_GridComp1      ! Single instance BRC object 

!
! !PUBLIC MEMBER FUNCTIONS:
!

   PUBLIC  BRC_GridCompSetServices
   PUBLIC  BRC_GridCompInitialize
   PUBLIC  BRC_GridCompRun1
   PUBLIC  BRC_GridCompRun2
   PUBLIC  BRC_GridCompFinalize

!
! !DESCRIPTION:
!
!  This module implements the (pre-ESMF) BRC Grid Component. 
!
! !REVISION HISTORY:
!
!  16Sep2003 da Silva  First crack.
!
!EOP
!-------------------------------------------------------------------------

  type BRC_GridComp1
        character(len=255) :: name
        character(len=255) :: iname           ! instance name
        character(len=255) :: rcfilen         ! resource file name
        character(len=255) :: regionsString   ! Comma-delimited string of regions

        integer :: instance                   ! instance number

        logical :: run_alarm = .false.        ! run alarm

        type(Chem_Mie), pointer :: mie_tables  => null() ! aod LUTs
        real, pointer :: biofuel_src(:,:)
        real, pointer :: biomass_src(:,:)
        real, pointer :: biomass_src_(:,:)
        real, pointer :: ebrcantt1_src(:,:)  ! level 1
        real, pointer :: ebrcantt2_src(:,:)  ! level 2
        real, pointer :: terpene_src(:,:)  ! level 2
        real, pointer :: brc_ship_src(:,:)
        real, pointer :: psoa_biob_voc(:,:,:) ! production of SOA from anthropogenic VOC
        real, pointer :: aviation_lto_src(:,:)  ! aviation - landing and takeoff
        real, pointer :: aviation_cds_src(:,:)  ! aviation - climbing and descent
        real, pointer :: aviation_crs_src(:,:)  ! aviation - cruise
        real          :: aviation_layers(4)     ! heights of the LTO, CDS and CRS layers
        
        real :: ratPOM               ! Ratio of POM to BRC mass
        real :: fHydrophobic         ! Fraction of emissions hydrophobic
        real :: fTerpene             ! Fraction of terpene emissions -> aerosol
        integer :: myDOW = -1             ! my Day of the week: Sun=1, Mon=2,...,Sat=7
        logical :: doing_nei=.FALSE.      ! NEI08: National Emission Inventory (US+Canada)
        real    :: nei_lon(2), nei_lat(2) ! NEI bounding box; superseeds ebrcantt1/2 inside
!       Workspace for any requested point emissions
!       -------------------------------------------
        logical :: doing_point_emissions=.FALSE.  ! Providing pointwise emissions
        character(len=255) :: point_emissions_srcfilen   ! filename for pointwise emissions
        integer                         :: nPts = -1
        integer, pointer, dimension(:)  :: vstart => null(), vend => null()
        real, pointer, dimension(:)     :: vLat  => null(), &
                                           vLon  => null(), &
                                           vBase => null(), &
                                           vTop  => null(), &
                                           vEmis => null()
 end type BRC_GridComp1

  type BRC_GridComp
     integer                     :: n = 0                ! number of instances 
     type(Chem_Mie), pointer     :: mie_tables => null() ! aod LUTs
     type(BRC_GridComp1), pointer :: gcs(:)     => null() ! instances
  end type BRC_GridComp

  character(len=*), parameter :: rc_basename = 'BRC_GridComp'

  real, parameter :: OCEAN=0.0, LAND = 1.0, SEA_ICE = 2.0
  real, parameter :: radToDeg = 57.2957795

CONTAINS

   subroutine BRC_GridCompSetServices(  gc, chemReg, rc)
   type(ESMF_GridComp), intent(INOUT) :: GC
   type(Chem_Registry), intent(INOUT) :: chemReg
   integer,             intent(OUT  ) :: rc

   CHARACTER(LEN=255) :: name

   integer            :: status
   character(len=ESMF_MAXSTR) :: Iam

   integer :: n,i

   type(ESMF_Config) :: cfg

   Iam = "BRC_GridCompSetServices"

!  Load resource file
!  ------------------
   cfg = ESMF_ConfigCreate(rc=status)
   VERIFY_(STATUS)
   call ESMF_ConfigLoadFile(cfg,trim(rc_basename)//'.rc',rc=status)
   VERIFY_(STATUS)

!  Parse resource file
!  -------------------
   n = ESMF_ConfigGetLen(cfg,label='BRC_instances:',rc=status)
   VERIFY_(STATUS)


!  We have 2 tracers for each instance of BRC
!  We cannot have fewer instances than half the number of
!   BRC bins in the registry (it is OK to have less, though)
!  --------------------------------------------------------
   if ( n .LT. chemReg%n_BRC/2 ) then
        rc = 35
        return
   else if ( n .GT. chemReg%n_BRC/2 ) then
        if (MAPL_AM_I_ROOT()) &
        print *, trim(Iam)// &
                 ': fewer BRC bins than possible BRC instances: ',&
                 n, chemReg%n_BRC/2
   end if
   n = min(n,chemReg%n_BRC/2 )

!  Record name of each instance
!  ----------------------------
   call ESMF_ConfigFindLabel(cfg,'BRC_instances:',rc=status)
   VERIFY_(STATUS)
   do i = 1, n
      call ESMF_ConfigGetAttribute(cfg,name,rc=status)
      VERIFY_(STATUS)
                                            ! resource file name
      IF(TRIM(name) == "full" ) THEN
       name = " "              ! blank instance name for full (1)
      ELSE
       name = TRIM(name)       ! instance name for others
      END IF
      call BRC_GridCompSetServices1_(gc,chemReg,name,rc=status)
      VERIFY_(STATUS)
   end do

!  Set profiling timers
!  --------------------
   call MAPL_TimerAdd(GC, name = '-BRC_TOTAL',           __RC__)
   call MAPL_TimerAdd(GC, name = '-BRC_RUN',             __RC__)
   call MAPL_TimerAdd(GC, name = '-BRC_INITIALIZE',      __RC__)
   call MAPL_TimerAdd(GC, name = '-BRC_FINALIZE',        __RC__)

   call MAPL_TimerAdd(GC, name = '-BRC_RUN1',            __RC__)
   call MAPL_TimerAdd(GC, name = '--BRC_EMISSIONS',      __RC__)

   call MAPL_TimerAdd(GC, name = '-BRC_RUN2',            __RC__)
   call MAPL_TimerAdd(GC, name = '--BRC_SETTLING',       __RC__)
   call MAPL_TimerAdd(GC, name = '--BRC_DRY_DEPOSITION', __RC__)
   call MAPL_TimerAdd(GC, name = '--BRC_WET_LS',         __RC__)
   call MAPL_TimerAdd(GC, name = '--BRC_WET_CV',         __RC__)
   call MAPL_TimerAdd(GC, name = '--BRC_DIAGNOSTICS',    __RC__)

!  All done
!  --------

   RETURN_(ESMF_SUCCESS)
   end subroutine BRC_GridCompSetServices

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  BRC_GridCompInitialize --- Initialize BRC_GridComp
!
! !INTERFACE:
!

   subroutine BRC_GridCompInitialize ( gcBRC, w_c, impChem, expChem, ggState, &
                                       nymd, nhms, cdt, rc )

! !USES:

  implicit NONE

! !INPUT PARAMETERS:

   type(Chem_Bundle), intent(inout) :: w_c        ! Chemical tracer fields      
   integer, intent(in) :: nymd, nhms           ! time
   real, intent(in) :: cdt                     ! chemistry timestep (secs)

! !OUTPUT PARAMETERS:

   type(BRC_GridComp), intent(inout) :: gcBRC   ! Grid Component
   type(ESMF_State), intent(inout)  :: impChem  ! Import State
   type(ESMF_State), intent(inout)  :: expChem  ! Export State
   type(MAPL_MetaComp), intent(inout) :: ggState
   integer, intent(out) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 - 

! !DESCRIPTION: Initializes the BRC Grid Component. It primarily sets
!               the import state for each active constituent package.
!
! !REVISION HISTORY:
!
!  18Sep2003 da Silva  First crack.
!
!EOP
!-------------------------------------------------------------------------

   character(len=*), parameter :: myname = 'BRC_GridCompInitialize'
   CHARACTER(LEN=255) :: name
   
   integer :: i, ier, n

   call MAPL_TimerOn(ggState, '-BRC_TOTAL')
   call MAPL_TimerOn(ggState, '-BRC_INITIALIZE')


!  Load resource file
!  ------------------
   call i90_loadf ( trim(rc_basename)//'.rc', ier )
   if ( ier .NE. 0 ) then
      rc = 10
      return
   end if

!  Parse resource file
!  -------------------
   CALL I90_label ( 'BRC_instances:', ier )
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
   
!  We have 2 tracers for each instance of BRC
!  We cannot have fewer instances than half the number of
!   BRC bins in the registry (it is OK to have less, though)
!  --------------------------------------------------------
   if ( n .LT. w_c%reg%n_BRC/2 ) then
        rc = 35
        return
   else if ( n .GT. w_c%reg%n_BRC/2 ) then
        if (MAPL_AM_I_ROOT()) &
        print *, trim(myname)// &
                 ': fewer BRC bin sets than possible BRC instances'//&
                 ' (2 bins per instance): ',&
                 n, w_c%reg%n_BRC
   end if
   n = min(n,w_c%reg%n_BRC/2 )
   gcBRC%n = n

!  Next allocate necessary memory
!  ------------------------------
   allocate ( gcBRC%gcs(n), stat=ier )    
   if ( ier .NE. 0 ) then
      rc = 40
      return
   end if

!  Record name of each instance
!  ----------------------------
   CALL I90_label ( 'BRC_instances:', ier )
   do i = 1, n
      CALL I90_gtoken( name, ier )
      if ( ier .NE. 0 ) then
         rc = 40
         return
      end if
                                            ! resource file name
      gcBRC%gcs(i)%rcfilen = trim(rc_basename)//'---'//trim(name)//'.rc'
      gcBRC%gcs(i)%instance = i              ! instance number 
      IF(TRIM(name) == "full" ) THEN
       gcBRC%gcs(i)%iname = " "              ! blank instance name for full (1)
      ELSE
       gcBRC%gcs(i)%iname = TRIM(name)       ! instance name for others
      END IF
   end do    

!  Next initialize each instance
!  -----------------------------
   do i = 1, gcBRC%n
      IF(MAPL_AM_I_ROOT()) THEN
       PRINT *," "
       PRINT *,myname,": Initializing instance ",TRIM(gcBRC%gcs(i)%iname)," [",gcBRC%gcs(i)%instance,"]"
      END IF
      call BRC_SingleInstance_ ( BRC_GridCompInitialize1_, i, &
                                gcBRC%gcs(i), w_c, impChem, expChem, ggState, &
                                nymd, nhms, cdt, ier )
      if ( ier .NE. 0 ) then
         rc = 1000+ier
         return
      end if
      gcBRC%gcs(i)%mie_tables => gcBRC%mie_tables
   end do

!  All done
!  --------
   CALL I90_FullRelease( ier )
   IF( ier /= 0 ) THEN
    PRINT *,myname,": I90_FullRelease not successful."
    rc = 40
   END IF

   call MAPL_TimerOff(ggState, '-BRC_INITIALIZE')
   call MAPL_TimerOff(ggState, '-BRC_TOTAL')

 end subroutine BRC_GridCompInitialize

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  BRC_GridCompRun1 --- Run BRC_GridComp
!
! !INTERFACE:
!

   subroutine BRC_GridCompRun1 ( gcBRC, w_c, impChem, expChem, ggState, &
                                 nymd, nhms, cdt, rc )

! !USES:

  IMPLICIT NONE

! !INPUT PARAMETERS:

   TYPE(Chem_Bundle), intent(in) :: w_c        ! Chemical tracer fields      
   INTEGER, INTENT(IN) :: nymd, nhms           ! time
   REAL,    INTENT(IN) :: cdt                  ! chemical timestep (secs)


! !OUTPUT PARAMETERS:

   TYPE(BRC_GridComp), INTENT(INOUT) :: gcBRC     ! Grid Component
   TYPE(ESMF_State), INTENT(INOUT)  :: impChem  ! Import State
   TYPE(ESMF_State), INTENT(INOUT)  :: expChem  ! Export State
   TYPE(MAPL_MetaComp), INTENT(INOUT) :: ggState
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

   call MAPL_TimerOn(ggState, '-BRC_TOTAL')
   call MAPL_TimerOn(ggState, '-BRC_RUN')

   do i = 1, gcBRC%n
      call BRC_SingleInstance_ ( BRC_GridCompRun1_, i, &
                                gcBRC%gcs(i), w_c, impChem, expChem, ggState, &
                                nymd, nhms, cdt, ier )
      if ( ier .NE. 0 ) then
         rc = i * 1000+ier
         return
      end if
   end do

   call MAPL_TimerOff(ggState, '-BRC_RUN')
   call MAPL_TimerOff(ggState, '-BRC_TOTAL')

 end subroutine BRC_GridCompRun1


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  BRC_GridCompRun2 --- Run BRC_GridComp
!
! !INTERFACE:
!

   subroutine BRC_GridCompRun2 ( gcBRC, w_c, impChem, expChem, ggState, &
                                 run_alarm, nymd, nhms, cdt, rc )

! !USES:

  IMPLICIT NONE

! !INPUT PARAMETERS:

   TYPE(Chem_Bundle), intent(in) :: w_c        ! Chemical tracer fields      
   LOGICAL, INTENT(IN) :: run_alarm            ! run alarm
   INTEGER, INTENT(IN) :: nymd, nhms           ! time
   REAL,    INTENT(IN) :: cdt                  ! chemical timestep (secs)


! !OUTPUT PARAMETERS:

   TYPE(BRC_GridComp), INTENT(INOUT) :: gcBRC     ! Grid Component
   TYPE(ESMF_State), INTENT(INOUT)  :: impChem  ! Import State
   TYPE(ESMF_State), INTENT(INOUT)  :: expChem  ! Export State
   TYPE(MAPL_MetaComp), INTENT(INOUT) :: ggState
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

   call MAPL_TimerOn(ggState, '-BRC_TOTAL')
   call MAPL_TimerOn(ggState, '-BRC_RUN')

   do i = 1, gcBRC%n
      gcBRC%gcs(i)%run_alarm = run_alarm

      call BRC_SingleInstance_ ( BRC_GridCompRun2_, i, &
                                gcBRC%gcs(i), w_c, impChem, expChem, ggState, &
                                nymd, nhms, cdt, ier )
      if ( ier .NE. 0 ) then
         rc = i * 1000+ier
         return
      end if
   end do

   call MAPL_TimerOff(ggState, '-BRC_RUN')
   call MAPL_TimerOff(ggState, '-BRC_TOTAL')

 end subroutine BRC_GridCompRun2



!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  BRC_GridCompFinalize --- Initialize BRC_GridComp
!
! !INTERFACE:
!

   subroutine BRC_GridCompFinalize ( gcBRC, w_c, impChem, expChem, ggState, &
                                     nymd, nhms, cdt, rc )

! !USES:

  IMPLICIT NONE

! !INPUT PARAMETERS:

   TYPE(Chem_Bundle), intent(in) :: w_c        ! Chemical tracer fields      
   INTEGER, INTENT(IN) :: nymd, nhms           ! time
   REAL,    INTENT(IN) :: cdt                  ! chemical timestep (secs)


! !OUTPUT PARAMETERS:

   TYPE(BRC_GridComp), INTENT(INOUT) :: gcBRC     ! Grid Component
   TYPE(ESMF_State), INTENT(INOUT)  :: impChem  ! Import State
   TYPE(ESMF_State), INTENT(INOUT)  :: expChem  ! Export State
   TYPE(MAPL_MetaComp), INTENT(INOUT) :: ggState
   INTEGER, INTENT(OUT) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 - 

! !DESCRIPTION: Finalizes the BRC Grid Component. Multiple instance
!               version.
!
! !REVISION HISTORY:
!
!  27Feb2008  da Silva  Introduced multiple instances
!
!EOP
!-------------------------------------------------------------------------

   integer i, ier

   call MAPL_TimerOn(ggState, '-BRC_TOTAL')
   call MAPL_TimerOn(ggState, '-BRC_FINALIZE')

   do i = 1, gcBRC%n
      call BRC_SingleInstance_ ( BRC_GridCompFinalize1_, i, &
                                gcBRC%gcs(i), w_c, impChem, expChem, ggState, &
                                nymd, nhms, cdt, ier )
      if ( ier .NE. 0 ) then
         rc = i * 1000+ier
         return
      end if
   end do

   if (associated(gcBRC%gcs)) deallocate ( gcBRC%gcs, stat=ier )
   gcBRC%n = -1

   call MAPL_TimerOff(ggState, '-BRC_FINALIZE')
   call MAPL_TimerOff(ggState, '-BRC_TOTAL')

 end subroutine BRC_GridCompFinalize


 subroutine BRC_GridCompSetServices1_(  gc, chemReg, iname, rc)
   type(ESMF_GridComp), intent(INOUT) :: GC
   type(Chem_Registry), intent(INOUT) :: chemReg
   character(len=*),    intent(IN   ) :: iname
   integer,             intent(OUT  ) :: rc

   ! local
   logical:: doing_nei

   integer :: Status
   character(len=ESMF_MAXSTR) :: Iam

   Iam ="BRC_GridCompSetServices1_"
 	
   call MAPL_AddImportSpec(GC, &
      SHORT_NAME = 'BRC_BIOMASS'//trim(iname), &
      LONG_NAME  = 'source species'  , &
      UNITS      = '1',                &
      DIMS       = MAPL_DimsHorzOnly,  &
      VLOCATION  = MAPL_VLocationNone, &
      RESTART    = MAPL_RestartSkip,   &
      RC         = STATUS)
   VERIFY_(STATUS)

   call MAPL_AddImportSpec(GC, &
      SHORT_NAME = 'BRC_TERPENE'//trim(iname), &
      LONG_NAME  = 'source species'  , &
      UNITS      = '1',                &
      DIMS       = MAPL_DimsHorzOnly,  &
      VLOCATION  = MAPL_VLocationNone, &
      RESTART    = MAPL_RestartSkip,   &
      RC         = STATUS)
   VERIFY_(STATUS)

   call MAPL_AddImportSpec(GC, &
      SHORT_NAME = 'BRC_BIOFUEL'//trim(iname), &
      LONG_NAME  = 'source species'  , &
      UNITS      = '1',                &
      DIMS       = MAPL_DimsHorzOnly,  &
      VLOCATION  = MAPL_VLocationNone, &
      RESTART    = MAPL_RestartSkip,   &
      RC         = STATUS)
   VERIFY_(STATUS)

   call MAPL_AddImportSpec(GC, &
      SHORT_NAME = 'BRC_ANTEBRC1'//trim(iname), &
      LONG_NAME  = 'source species'  , &
      UNITS      = '1',                &
      DIMS       = MAPL_DimsHorzOnly,  &
      VLOCATION  = MAPL_VLocationNone, &
      RESTART    = MAPL_RestartSkip,   &
      RC         = STATUS)
   VERIFY_(STATUS)

   call MAPL_AddImportSpec(GC, &
      SHORT_NAME = 'BRC_ANTEBRC2'//trim(iname), &
      LONG_NAME  = 'source species'  , &
      UNITS      = '1',                &
      DIMS       = MAPL_DimsHorzOnly,  &
      VLOCATION  = MAPL_VLocationNone, &
      RESTART    = MAPL_RestartSkip,   &
      RC         = STATUS)
   VERIFY_(STATUS)

   call MAPL_AddImportSpec(GC, &
      SHORT_NAME = 'BRC_SHIP'//trim(iname), &
      LONG_NAME  = 'source species'  , &
      UNITS      = '1',                &
      DIMS       = MAPL_DimsHorzOnly,  &
      VLOCATION  = MAPL_VLocationNone, &
      RESTART    = MAPL_RestartSkip,   &
      RC         = STATUS)
   VERIFY_(STATUS)

   call MAPL_AddImportSpec(GC, &
      SHORT_NAME = 'BRC_AVIATION_LTO'//trim(iname), &
      LONG_NAME  = 'brc_aviation_lto' , &
      UNITS      = '1',                &
      DIMS       = MAPL_DimsHorzOnly,  &
      VLOCATION  = MAPL_VLocationNone, &
      RESTART    = MAPL_RestartSkip,   &
      RC         = STATUS)
   VERIFY_(STATUS)

   call MAPL_AddImportSpec(GC, &
      SHORT_NAME = 'BRC_AVIATION_CDS'//trim(iname), &
      LONG_NAME  = 'brc_aviation_cds' , &
      UNITS      = '1',                &
      DIMS       = MAPL_DimsHorzOnly,  &
      VLOCATION  = MAPL_VLocationNone, &
      RESTART    = MAPL_RestartSkip,   &
      RC         = STATUS)
   VERIFY_(STATUS)

   call MAPL_AddImportSpec(GC, &
      SHORT_NAME = 'BRC_AVIATION_CRS'//trim(iname), &
      LONG_NAME  = 'brc_aviation_crs' , &
      UNITS      = '1',                &
      DIMS       = MAPL_DimsHorzOnly,  &
      VLOCATION  = MAPL_VLocationNone, &
      RESTART    = MAPL_RestartSkip,   &
      RC         = STATUS)
   VERIFY_(STATUS)

   call MAPL_AddImportSpec(GC, &
     SHORT_NAME = 'pSOA_BIOB_VOC'//trim(iname), &
     LONG_NAME  = 'Production of SOA from Biomass Burning VOC' , &
     UNITS      = 'kg m-3 s-1',       &
     DIMS       = MAPL_DimsHorzVert,  &
     VLOCATION  = MAPL_VLocationCenter, &
     RESTART    = MAPL_RestartSkip,   &
     RC         = STATUS)
  VERIFY_(STATUS)
  

!  Parse the resource file to see if NEI imports are required
!  ----------------------------------------------------------
   call doing_nei_(trim(rc_basename), trim(iname), doing_nei, __RC__)

   NEI_EMISSIONS: if (doing_nei) then
   call MAPL_AddImportSpec(GC, &
       SHORT_NAME = 'BRC_NEI_BOT'//trim(iname), &
       LONG_NAME  = 'brc_nei_bot' , &
       UNITS      = '1',                &
       DIMS       = MAPL_DimsHorzOnly,  &
       VLOCATION  = MAPL_VLocationNone, &
       RESTART    = MAPL_RestartSkip,   &
       RC         = STATUS)
   VERIFY_(STATUS)

   call MAPL_AddImportSpec(GC, &
       SHORT_NAME = 'BRC_NEI_TOP'//trim(iname), &
       LONG_NAME  = 'brc_nei_top' , &
       UNITS      = '1',                &
       DIMS       = MAPL_DimsHorzOnly,  &
       VLOCATION  = MAPL_VLocationNone, &
       RESTART    = MAPL_RestartSkip,   &
       RC         = STATUS)
   VERIFY_(STATUS)
   end if NEI_EMISSIONS


  RETURN_(ESMF_SUCCESS)

 contains
   subroutine doing_nei_(rcbasen, iname, result, rc)

   character(len=*), intent(in) :: rcbasen
   character(len=*), intent(in) :: iname
   logical, intent(out)         :: result
   integer, intent(out)         :: rc

   ! local
   type(ESMF_Config)  :: cfg
   character(len=255) :: name
   logical            :: isPresent
   integer            :: status
   character(len=255) :: Iam

   Iam = 'BRC_GridCOmpSetServices1_::doing_nei_'

   if (iname == '') then
       name = 'full'
   else 
       name = trim(iname)
   end if

   name = trim(rcbasen)//'---'//trim(name)//'.rc'

   cfg = ESMF_ConfigCreate(__RC__)
   call ESMF_ConfigLoadFile(cfg, trim(name), __RC__)
   call ESMF_ConfigFindLabel(cfg, 'nei_boundingbox:', isPresent=isPresent, __RC__)

   if (isPresent) then
       result = .true.
   else 
       result = .false.
   end if

   RETURN_(ESMF_SUCCESS)
   end subroutine doing_nei_

 end subroutine BRC_GridCompSetServices1_

!--------------------------------------------------------------------------

!                      Single Instance Methods

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  BRC_GridCompInitialize --- Initialize BRC_GridComp
!
! !INTERFACE:
!

   subroutine BRC_GridCompInitialize1_ ( gcBRC, w_c, impChem, expChem, ggState, &
                                         nymd, nhms, cdt, rc )

! !USES:

  implicit NONE

! !INPUT PARAMETERS:

   type(Chem_Bundle), intent(inout) :: w_c     ! Chemical tracer fields      
   integer, intent(in) :: nymd, nhms           ! time
   real, intent(in) :: cdt                     ! chemistry timestep (secs)

! !OUTPUT PARAMETERS:

   type(BRC_GridComp1), intent(inout) :: gcBRC    ! Grid Component
   type(ESMF_State), intent(inout)  :: impChem  ! Import State
   type(ESMF_State), intent(inout)  :: expChem  ! Export State
   type(MAPL_MetaComp), intent(inout) :: ggState
   integer, intent(out) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 - 

! !DESCRIPTION: Initializes the BRC Grid Component. It primarily sets
!               the import state for each active constituent package.
!
! !REVISION HISTORY:
!
!  18Sep2003 da Silva  First crack.
!
!EOP
!-------------------------------------------------------------------------

   character(len=*), parameter :: myname = 'BRC_GridCompInitialize1'


   character(len=255) :: rcfilen
   integer :: n
   integer :: i1, i2, im, j1, j2, jm, nbins, n1, n2, nbins_rc, km
   integer, allocatable :: ier(:)
   real :: qmax, qmin
   LOGICAL :: NoRegionalConstraint 

   rcfilen = gcBRC%rcfilen
   gcBRC%name = 'BRC Constituent Package'

!  Initialize local variables
!  --------------------------
   rc = 0
   i1 = w_c%grid%i1; i2 = w_c%grid%i2; im = w_c%grid%im
   j1 = w_c%grid%j1; j2 = w_c%grid%j2; jm = w_c%grid%jm
   km = w_c%grid%km
   nbins = w_c%reg%n_BRC
   n1  = w_c%reg%i_BRC
   n2  = w_c%reg%j_BRC

   call init_()
   if ( rc /= 0 ) return


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

   call i90_label ( 'number_brc_classes:', ier(1) )
   nbins_rc = i90_gint ( ier(2) )
   if ( any(ier(1:2) /= 0) ) then
      call final_(20)
      return
   end if
   if ( nbins_rc /= nbins ) then
      call final_(25)
      return
   end if

!  Aircraft emissions
!  ------------------
   ier(:) = 0
   call i90_label  ( 'aviation_vertical_layers:', ier(1) )
   gcBRC%aviation_layers(1) = i90_gfloat(ier(2))
   gcBRC%aviation_layers(2) = i90_gfloat(ier(3))
   gcBRC%aviation_layers(3) = i90_gfloat(ier(4))
   gcBRC%aviation_layers(4) = i90_gfloat(ier(5))

   if ( any(ier(1:5) /= 0) ) then
         call final_(77)
         return
   end if

!  Handle Point-wise Emission Sources Specified in a Text File
!  -----------------------------------------------------------
   ier(:) = 0
   call i90_label  ( 'point_emissions_srcfilen:',   ier(1) )
   call i90_gtoken ( gcBRC%point_emissions_srcfilen, ier(2) )
   if ( ier(1) /= 0 ) then
        gcBRC%doing_point_emissions = .FALSE. ! if rc is missing, don't fuss
   else if ( any(ier(2:2) /= 0) ) then
         call final_(42) ! this means point emissions info is messed up, abort
         return
   else
         if ( (index(gcBRC%point_emissions_srcfilen,'/dev/null')>0) ) then
               gcBRC%doing_point_emissions = .FALSE. ! disable it if no file specified
         else
               gcBRC%doing_point_emissions = .TRUE.  ! we are good to go
         end if
   end if

!  Handle NEI08 Emissions
!  ----------------------
   ier(:) = 0
   call i90_label  ( 'nei_boundingbox:',   ier(1) )
   gcBRC%nei_lon(1) = i90_gfloat(ier(2))
   gcBRC%nei_lon(2) = i90_gfloat(ier(3))
   gcBRC%nei_lat(1) = i90_gfloat(ier(4))
   gcBRC%nei_lat(2) = i90_gfloat(ier(5))
   if ( ier(1) /= 0 ) then
        gcBRC%doing_nei = .FALSE. ! if rc is missing, don't fuss
   else if ( any(ier(2:5) /= 0) ) then
         call final_(42) ! this means NEI info is messed up, abort
         return
   else
! --------------------------------------------------------------------------   
!        if ( (index(gcBRC%nei_srcfilen(1),'/dev/null')>0) .or. &
!             (index(gcBRC%nei_srcfilen(2),'/dev/null')>0) ) then 
!              gcBRC%doing_nei = .FALSE. ! disable it if no file specified
!        else
!              gcBRC%doing_nei = .TRUE.  ! we are good to go
!        end if
! -------------------------------------------------------------------------- 
! TODO: Need to parse the ExtData file to replicate the above logic,
!       until then do not include the NOI datasets in the ExtData primary 
!       export tables
! --------------------------------------------------------------------------

         gcBRC%doing_nei = .TRUE.  ! we are good to go
   end if

   if ( MAPL_AM_I_ROOT() ) then
    if ( gcBRC%doing_nei ) then
      print *, 'BRC_GridComp: using NEI08 Emissions over North America'
    else
      print *, 'BRC_GridComp: skipping NEI08 Emissions over North America'
    end if
   end if

!                          -------

!  Day of the week to reset tracer to zero
!  ---------------------------------------
   call i90_label ( 'my_day_of_the_week:',ier(1))
   if ( ier(1) /= 0 ) then
        gcBRC%myDOW = -1   ! by default never reset tracer to zero
   else
        gcBRC%myDOW = i90_gint (ier(1))
        if ( ier(1) /= 0 ) then
           call final_(60)
           return
        end if
   end if

!                          -------


!  Ratio of POM to BRC mass
!  -----------------------
   call i90_label ( 'pom_brc_ratio:', ier(1) )
   gcBRC%ratPOM = i90_gfloat ( ier(2) )
   if ( any(ier(1:2) /= 0) ) then
      call final_(50)
      return
   end if


!  Hydrophilic fraction
!  ---------------
   call i90_label ( 'hydrophobic_fraction:', ier(1) )
   gcBRC%fHydrophobic = i90_gfloat ( ier(2) )
   if ( any(ier(1:2) /= 0) ) then
      call final_(50)
      return
   end if


!  Terpene Emission Factor
!  ---------------
   call i90_label ( 'terpene_emission_fraction:', ier(1) )
   gcBRC%fTerpene = i90_gfloat ( ier(2) )
   if ( any(ier(1:2) /= 0) ) then
      call final_(50)
      return
   end if


!                          -------


!  Scavenging Efficiency
!  To be used in convtran.F90, this parameter
!  is the scavenging efficiency of the tracer [km -1]
!  ---------------
   call i90_label ( 'fscav:', ier(1) )
   do n = 1, nbins
      w_c%reg%fscav(n1+n-1) = i90_gfloat ( ier(n+1) )
   end do
   if ( any(ier(1:nbins+1) /= 0) ) then
      call final_(50)
      return
   end if

!  Particle density
!  To be used in droplet activation code
!  ---------------
   call i90_label ( 'particle_density:', ier(1) )
   do n = 1, nbins
      w_c%reg%rhop(n1+n-1)  = i90_gfloat ( ier(n+1) )
   end do
   if ( any(ier(1:nbins+1) /= 0) ) then
      call final_(50)
      return
   end if

!  Number median radius
!  To be used in droplet activation code
!  ---------------
   call i90_label ( 'particle_radius_number:', ier(1) )
   do n = 1, nbins
      w_c%reg%rmed(n1+n-1)  = i90_gfloat ( ier(n+1) ) * 1e-6
   end do
   if ( any(ier(1:nbins+1) /= 0) ) then
      call final_(50)
      return
   end if

!  Sigma (lognormal mode width)
!  To be used in droplet activation code
!  ---------------
   call i90_label ( 'sigma:', ier(1) )
   do n = 1, nbins
      w_c%reg%sigma(n1+n-1)  = i90_gfloat ( ier(n+1) )
   end do
   if ( any(ier(1:nbins+1) /= 0) ) then
      call final_(50)
      return
   end if

!  Number to mass conversion factor
!  To be used in droplet activation code
!  ---------------
   call i90_label ( 'fnum:', ier(1) )
   do n = 1, nbins
      w_c%reg%fnum(n1+n-1)  = i90_gfloat ( ier(n+1) )
   end do
   if ( any(ier(1:nbins+1) /= 0) ) then
      call final_(50)
      return
   end if

!  Molecular weight
!  To be used in droplet activation code
!  ---------------
   call i90_label ( 'molecular_weight:', ier(1) )
   do n = 1, nbins
      w_c%reg%molwght(n1+n-1)  = i90_gfloat ( ier(n+1) )
   end do
   if ( any(ier(1:nbins+1) /= 0) ) then
      call final_(50)
      return
   end if

!                          -------

!  Grab the region string.
!  -----------------------
   ier(:)=0
   call i90_label ( 'BRC_regions_indices:', ier(1) )
   CALL I90_gtoken( gcBRC%regionsString, ier(2) )
   IF( ANY(ier(1:2) < 0 ) ) THEN
    CALL final_(51)
    RETURN
   END IF

!  Is this instantiation a global case?
!  -----------------------------------
   IF(gcBRC%regionsString(1:2) == "-1") THEN
    NoRegionalConstraint = .TRUE.
   ELSE
    SELECT CASE (ESMF_UtilStringLowerCase(gcBRC%regionsString(1:2)))
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
   IF(NoRegionalConstraint) gcBRC%regionsString = "-1"

   IF(MAPL_AM_I_ROOT()) THEN
    IF(NoRegionalConstraint) THEN
     PRINT *,myname,": This instantiation has no regional constraints."
    ELSE
     PRINT *,myname,": This instantiation is regionally constrained."
     PRINT *,myname,": List of region numbers included: ",TRIM(gcBRC%regionsString)
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
   allocate ( gcBRC%biomass_src(i1:i2,j1:j2), gcBRC%biofuel_src(i1:i2,j1:j2), &
              gcBRC%biomass_src_(i1:i2,j1:j2), &
              gcBRC%ebrcantt1_src(i1:i2,j1:j2), gcBRC%ebrcantt2_src(i1:i2,j1:j2), &
              gcBRC%terpene_src(i1:i2,j1:j2), gcBRC%brc_ship_src(i1:i2,j1:j2), &
              gcBRC%psoa_biob_voc(i1:i2,j1:j2,km), &
              gcBRC%aviation_lto_src(i1:i2,j1:j2), &
              gcBRC%aviation_cds_src(i1:i2,j1:j2), &
              gcBRC%aviation_crs_src(i1:i2,j1:j2), ier(nerr), stat=ios )
   if ( ios /= 0 ) rc = 100
   end subroutine init_

   subroutine final_(ierr)
   integer :: ierr
   integer ios
   deallocate ( gcBRC%biomass_src, gcBRC%biofuel_src, &
                gcBRC%biomass_src_, &
                gcBRC%ebrcantt1_src, gcBRC%ebrcantt2_src, &
                gcBRC%terpene_src, gcBRC%brc_ship_src, &
                gcBRC%psoa_biob_voc, &
                gcBRC%aviation_lto_src, &
                gcBRC%aviation_cds_src, &
                gcBRC%aviation_crs_src, ier, stat=ios )
   call i90_release()
   rc = ierr
   end subroutine final_

   end subroutine BRC_GridCompInitialize1_

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  BRC_GridCompRun1_ --- The Chem Driver, run phase 1
!
! !INTERFACE:
!

   subroutine BRC_GridCompRun1_ ( gcBRC, w_c, impChem, expChem, ggState, &
                                  nymd, nhms, cdt, rc )

! !USES:

  implicit NONE

! !INPUT/OUTPUT PARAMETERS:

   type(BRC_GridComp1), intent(inout) :: gcBRC   ! Grid Component
   type(Chem_Bundle), intent(inout) :: w_c      ! Chemical tracer fields
   type(MAPL_MetaComp), intent(inout) :: ggState

! !INPUT PARAMETERS:

   type(ESMF_State), intent(inout) :: impChem    ! Import State
   integer, intent(in) :: nymd, nhms          ! time
   real, intent(in) :: cdt                    ! chemistry timestep (secs)

! !OUTPUT PARAMETERS:

   type(ESMF_State), intent(inout) :: expChem   ! Export State
   integer, intent(out) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 -
 
! !DESCRIPTION: This routine implements the so-called BRC Driver. That 
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

   character(len=*), parameter :: myname = 'BRC_GridCompRun1_'
   character(len=*), parameter :: Iam = myname

   integer :: i1, i2, im, j1, j2, jm, nbins, n1, n2, km, n
   integer :: i, j, ijl, ijkl, ijk1l
   real :: qmax, qmin

!  Input fields from fvGCM
!  -----------------------
   real, pointer, dimension(:,:)   :: pblh
   real, pointer, dimension(:,:,:) :: tmpu, rhoa, ple, hghte

!  Workspace for NEI emissions
!  ---------------------------
   real, pointer, dimension(:,:)         ::  nei_src1, nei_src2

   integer          :: idow
   character(len=3) :: cdow

   real, pointer :: var2d(:,:) => null()


#define EXPORT        expChem
#define iNAME         TRIM(gcBRC%iname)

#define ptrBRCEM       BRC_emis

#define ptrBRCEMAN     BRC_emisAN
#define ptrBRCEMBB     BRC_emisBB
#define ptrBRCEMBF     BRC_emisBF
#define ptrBRCEMBG     BRC_emisBG

   integer :: STATUS

#include "BRC_GetPointer___.h"

   call MAPL_TimerOn(ggState, '-BRC_RUN1')

!  Initialize local variables
!  --------------------------
   rc = 0
   i1 = w_c%grid%i1; i2 = w_c%grid%i2; im = w_c%grid%im
   j1 = w_c%grid%j1; j2 = w_c%grid%j2; jm = w_c%grid%jm

   km    = w_c%grid%km
   nbins = w_c%reg%n_BRC
   n1    = w_c%reg%i_BRC
   n2    = w_c%reg%j_BRC

   ijl   = ( i2 - i1 + 1 ) * ( j2 - j1 + 1 )
   ijkl  = ijl * km
   ijk1l = ijl * (km+1)

! Reset tracer to zero at 0Z on specific day of week
! --------------------------------------------------
  idow = Chem_UtilIdow(nymd)
  if ( (nhms==0) .and. (idow == gcBRC%myDOW) ) then
        cdow = Chem_UtilCdow(nymd)
        do n = n1, n2
           w_c%qa(n)%data3d(i1:i2,j1:j2,1:km) = tiny(1.) ! avoid division by zero
        end do
        if ( MAPL_AM_I_ROOT() ) then
           print *, '<> BRC '//cdow//' tracer being set to zero on ', nymd, nhms
        end if
  end if

! Update emissions/production if necessary (daily)
! ------------------------------------------------
   call MAPL_TimerOn(ggState, '--BRC_EMISSIONS')

!   Biomass Burning -- select on known inventories
!   ----------------------------------------------

    call MAPL_GetPointer(impChem, var2d, 'BRC_BIOMASS'//iNAME, __RC__)
    gcBRC%biomass_src = var2d


!   Terpene, biofuel and anthropogenic emissions (inventories)
!   ----------------------------------------------------------
    call MAPL_GetPointer(impChem, var2d, 'BRC_TERPENE'//iNAME, __RC__)
    gcBRC%terpene_src = var2d

    call MAPL_GetPointer(impChem, var2d, 'BRC_BIOFUEL'//iNAME, __RC__)
    gcBRC%biofuel_src = var2d

    call MAPL_GetPointer(impChem, var2d, 'BRC_ANTEBRC1'//iNAME, __RC__)
    gcBRC%ebrcantt1_src = var2d

    call MAPL_GetPointer(impChem, var2d, 'BRC_ANTEBRC2'//iNAME, __RC__)
    gcBRC%ebrcantt2_src = var2d

!   Ship based BRC emissions
    call MAPL_GetPointer(impChem, var2d, 'BRC_SHIP'//iNAME, __RC__)
    gcBRC%brc_ship_src = var2d

!   Aircraft emissions during the three phases of flight
    call MAPL_GetPointer(impChem, var2d, 'BRC_AVIATION_LTO'//iNAME, __RC__)
    gcBRC%aviation_lto_src = var2d

    call MAPL_GetPointer(impChem, var2d, 'BRC_AVIATION_CDS'//iNAME, __RC__)
    gcBRC%aviation_cds_src = var2d

    call MAPL_GetPointer(impChem, var2d, 'BRC_AVIATION_CRS'//iNAME, __RC__)
    gcBRC%aviation_crs_src = var2d
    
    
!   As a safety check, where value is undefined set to 0
    do j = j1, j2
     do i = i1, i2
      if(1.01*gcBRC%biomass_src(i,j) .gt. undefval) gcBRC%biomass_src(i,j) = 0.
      if(1.01*gcBRC%terpene_src(i,j) .gt. undefval) gcBRC%terpene_src(i,j) = 0.
      if(1.01*gcBRC%biofuel_src(i,j) .gt. undefval) gcBRC%biofuel_src(i,j) = 0.
      if(1.01*gcBRC%ebrcantt1_src(i,j) .gt. undefval) gcBRC%ebrcantt1_src(i,j) = 0.
      if(1.01*gcBRC%ebrcantt2_src(i,j) .gt. undefval) gcBRC%ebrcantt2_src(i,j) = 0.
      if(1.01*gcBRC%brc_ship_src(i,j) .gt. undefval) gcBRC%brc_ship_src(i,j) = 0.
      if(1.01*gcBRC%aviation_lto_src(i,j) .gt. undefval) gcBRC%aviation_lto_src(i,j) = 0.
      if(1.01*gcBRC%aviation_cds_src(i,j) .gt. undefval) gcBRC%aviation_cds_src(i,j) = 0.
      if(1.01*gcBRC%aviation_crs_src(i,j) .gt. undefval) gcBRC%aviation_crs_src(i,j) = 0.
     enddo
    enddo


#ifdef DEBUG
    call pmaxmin('BRC: biomass', gcBRC%biomass_src, qmin, qmax, ijl,1, 1. )
    call pmaxmin('BRC: biofuel', gcBRC%biofuel_src, qmin, qmax, ijl,1, 1. )
    call pmaxmin('BRC: ebrcantt1', gcBRC%ebrcantt1_src, qmin, qmax, ijl,1,1.)
    call pmaxmin('BRC: ebrcantt2', gcBRC%ebrcantt2_src, qmin, qmax, ijl,1,1.)
    call pmaxmin('BRC: terpene', gcBRC%terpene_src, qmin, qmax, ijl,1, 1.)
    call pmaxmin('BRC: brc_ship', gcBRC%brc_ship_src, qmin, qmax, ijl,1, 1.)
    call pmaxmin('BRC: avi_lto', gcBRC%aviation_lto_src, qmin, qmax, ijl,1,1.)
    call pmaxmin('BRC: avi_cds', gcBRC%aviation_cds_src, qmin, qmax, ijl,1,1.)
    call pmaxmin('BRC: avi_crs', gcBRC%aviation_crs_src, qmin, qmax, ijl,1,1.)
#endif

!   Save this in case we need to apply diurnal cycle
!   ------------------------------------------------
   if ( w_c%diurnal_bb ) then
        gcBRC%biomass_src_(:,:) = gcBRC%biomass_src(:,:)
   end if

!  Apply diurnal cycle if so desired
!  ---------------------------------
   if ( w_c%diurnal_bb ) then
      call Chem_BiomassDiurnal ( gcBRC%biomass_src, gcBRC%biomass_src_,   &
                                 w_c%grid%lon(:,:)*radToDeg, &
                                 w_c%grid%lat(:,:)*radToDeg, nhms, cdt )      
   end if

!  Read any pointwise emissions, if requested
!  ------------------------------------------
   if(gcBRC%doing_point_emissions) then
    call Chem_UtilPointEmissions( nymd, gcBRC%point_emissions_srcfilen, &
                                  gcBRC%nPts, gcBRC%vLat, gcBRC%vLon, &
                                  gcBRC%vBase, gcBRC%vTop, gcBRC%vEmis, &
                                  gcBRC%vStart, gcBRC%vEnd )

!   In case vStart or vEnd were not specified in the file set to defaults
    where(gcBRC%vStart < 0) gcBRC%vStart = 000000
    where(gcBRC%vEnd < 0)   gcBRC%vEnd   = 240000
   endif


!  Apply NEI emissions over North America if so desired
!  ----------------------------------------------------
   if (gcBRC%doing_NEI) then

       allocate(nei_src1(i1:i2,j1:j2), nei_src2(i1:i2,j1:j2), __STAT__)

       call MAPL_GetPointer(impChem,var2d,'BRC_NEI_BOT'//iNAME, __RC__)
       nei_src1 = var2d

       call MAPL_GetPointer(impChem,var2d,'BRC_NEI_TOP'//iNAME, __RC__)
       nei_src2 = var2d

       where ( (w_c%grid%lon >= gcBRC%nei_lon(1)) .and. &
               (w_c%grid%lon <= gcBRC%nei_lon(2)) .and. &
               (w_c%grid%lat >= gcBRC%nei_lat(1)) .and. &
               (w_c%grid%lat <= gcBRC%nei_lat(2))   )

               gcBRC%ebrcantt1_src = nei_src1
               gcBRC%ebrcantt2_src = nei_src2
       end where

#ifdef DEBUG
            call pmaxmin('BRC: nei_bot', nei_src1, qmin, qmax, ijl,1, 1. )
            call pmaxmin('BRC: nei_top', nei_src2, qmin, qmax, ijl,1, 1. )
#endif

            deallocate(nei_src1, nei_src2)

   end if ! doing NEI

#ifdef DEBUG
   do n = n1, n2
      call pmaxmin('BRC: q_beg', w_c%qa(n)%data3d(i1:i2,j1:j2,1:km), qmin, qmax, &
                   ijl, km, 1. )
   end do
#endif


!  Get 2D Imports
!  --------------
   call MAPL_GetPointer ( impChem, pblh,     'ZPBL',     __RC__ )

!  Get 3D Imports
!  --------------
   call MAPL_GetPointer ( impChem, tmpu,     'T',        __RC__ )
   call MAPL_GetPointer ( impChem, rhoa,     'AIRDENS',  __RC__ )
   call MAPL_GetPointer ( impChem, ple,      'PLE',      __RC__ )
   call MAPL_GetPointer ( impChem, hghte,    'ZLE',      __RC__ )

  

#ifdef DEBUG

   call pmaxmin('BRC: pblh       ', pblh    , qmin, qmax, ijl,1, 1. )

   call pmaxmin('BRC: tmpu       ', tmpu    , qmin, qmax, ijkl,1, 1. )
   call pmaxmin('BRC: rhoa       ', rhoa    , qmin, qmax, ijkl,1, 1. )
   call pmaxmin('BRC: hghte      ', hghte   , qmin, qmax, ijkl,1, 1. )

#endif

!  BRC Source
!  -----------
   call BRC_Emission ( i1, i2, j1, j2, km, nbins, cdt, gcBRC, w_c, &
                      pblh, tmpu, rhoa, hghte, BRC_emis, &
                      BRC_emisAN, BRC_emisBB, BRC_emisBF, BRC_emisBG, rc )
#ifdef DEBUG
   do n = n1, n2
      call pmaxmin('BRC: q_emi', w_c%qa(n)%data3d(i1:i2,j1:j2,1:km), qmin, qmax, &
                    ijl, km, 1. )
   end do
#endif

   call MAPL_TimerOff(ggState, '--BRC_EMISSIONS')
   
   call MAPL_TimerOff(ggState, '-BRC_RUN1')

!  All done
!  --------
   return

CONTAINS

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  BRC_Emission - Adds Organic Carbon emission for one timestep
!             We have emissions from 6 sources, which are distributed
!             differently in the vertical
!             1) biomass burning - uniformly mixed in PBL
!             2) biofuel sources - emitted into lowest 100 m
!             3) anthropogenic l1 - emitted into lowest 100 m
!             4) anthropogenic l2 - emitted into 100 - 500 m levels
!             5) terpene          - emitted to surface (hydrophilic only)
!             6) point sources    - emitted in altitudes specified in input
!
! !INTERFACE:
!

   subroutine BRC_Emission ( i1, i2, j1, j2, km, nbins, cdt, gcBRC, w_c, &
                             pblh, tmpu, rhoa, hghte, BRC_emis, &
                             BRC_emisAN, BRC_emisBB, BRC_emisBF, BRC_emisBG, rc )

! !USES:

  implicit NONE

! !INPUT PARAMETERS:

   integer, intent(in) :: i1, i2, j1, j2, km, nbins
   real, intent(in)    :: cdt
   type(BRC_GridComp1), intent(in)    :: gcBRC       ! BRC Grid Component
   real, pointer, dimension(:,:)    :: pblh
   real, pointer, dimension(:,:,:)  :: tmpu
   real, pointer, dimension(:,:,:)  :: rhoa
   real, pointer, dimension(:,:,:)  :: hghte

! !OUTPUT PARAMETERS:

   type(Chem_Bundle), intent(inout) :: w_c         ! Chemical tracer fields
   type(Chem_Array), intent(inout)  :: BRC_emis(nbins) ! BRC emissions, kg/m2/s
   type(Chem_Array), intent(inout)  :: BRC_emisAN      ! BRC emissions, kg/m2/s
   type(Chem_Array), intent(inout)  :: BRC_emisBB      ! BRC emissions, kg/m2/s
   type(Chem_Array), intent(inout)  :: BRC_emisBF      ! BRC emissions, kg/m2/s
   type(Chem_Array), intent(inout)  :: BRC_emisBG      ! BRC emissions, kg/m2/s
   integer, intent(out)             :: rc          ! Error return code:
                                                   !  0 - all is well
                                                   !  1 - 
   character(len=*), parameter :: myname = 'BRC_Emission'

! !DESCRIPTION: Updates the BRC concentration with emissions every timestep
!
! !REVISION HISTORY:
!
!  06Nov2003, Colarco
!  Based on Ginoux
!
!EOP
!-------------------------------------------------------------------------

! !Local Variables
   integer  ::  i, j, k, m, n, ios, ijl, ii
   integer  ::  n1, n2
!  pressure at 100m, 500m, & PBLH
   real, dimension(i1:i2,j1:j2) :: p100, p500, pPBL  
   real, dimension(i1:i2,j1:j2) :: p0, z0, ps
   real :: p1, z1, dz, delz, delp, f100, f500, fPBL, fBot
   real :: qmax, qmin, eBiofuel, eBiomass, eTerpene, eAnthro

   real, dimension(i1:i2,j1:j2) :: factor, srcHydrophobic, srcHydrophilic
   real, dimension(i1:i2,j1:j2) :: srcBiofuel, srcBiomass, srcAnthro, srcBiogenic
   real                         :: srcTmp, zpbl, maxAll

   real, dimension(i1:i2,j1:j2,km) :: emis_aviation
   real, dimension(i1:i2,j1:j2,km) :: srcAviation
   real                            :: z_lto_bot, z_lto_top
   real                            :: z_cds_bot, z_cds_top
   real                            :: z_crs_bot, z_crs_top

   real, dimension(i1:i2,j1:j2)          :: f_bb_        ! scaling factor for BB emissions based on maximum allowed exttau
   real, dimension(i1:i2,j1:j2)          :: exttau_bb_   ! increment of exttau due to BB during the current time step
   real, allocatable, dimension(:,:,:,:) :: qa_bb_       ! increment of qa due to BB during the current time step (nbins,i1:i2,j1:j2:km)
   real                                  :: cutoff_bb_exttau
   integer                               :: nch, idx
   real                                  :: ilam550
   real                                  :: tau, ssa
   character(len=255)                    :: qname
   real, parameter                       :: max_bb_exttau = 30.0

!  Indices for point emissions
   integer, pointer, dimension(:)  :: iPoint, jPoint
   real, dimension(km)             :: point_column_emissions

!  Source function terms for SOA from Anthropogenic VOCs
   real :: srcSOAanthro = 0.0

!  Initialize local variables
!  --------------------------
   n1  = w_c%reg%i_BRC
   n2  = w_c%reg%j_BRC
   ijl = ( i2 - i1 + 1 ) * ( j2 - j1 + 1 )

!  Emission factors scaling from source files to desired mass quantity
   eBiomass = gcBRC%ratPOM
   eBiofuel = gcBRC%ratPOM 
   eTerpene = gcBRC%ratPOM * gcBRC%fTerpene
   eAnthro  = gcBRC%ratPOM

!  Zero diagnostic accumulators
   do n = 1, nbins
     if( associated(BRC_emis(n)%data2d) ) BRC_emis(n)%data2d = 0.0
   end do
     if(associated(BRC_emisAN%data2d) )   BRC_emisAN%data2d  = 0.0
     if(associated(BRC_emisBF%data2d) )   BRC_emisBF%data2d  = 0.0
     if(associated(BRC_emisBB%data2d) )   BRC_emisBB%data2d  = 0.0
     if(associated(BRC_emisBG%data2d) )   BRC_emisBG%data2d  = 0.0

!  Distribute aircraft emissions from LTO, CDS and CRS layers
!  ----------------------------------------------------------
   z_lto_bot = max(1e-3, gcBRC%aviation_layers(1))
   z_lto_top = max(2e-3, gcBRC%aviation_layers(2))

   z_cds_bot = max(2e-3, gcBRC%aviation_layers(2))
   z_cds_top = max(3e-3, gcBRC%aviation_layers(3))

   z_crs_bot = max(3e-3, gcBRC%aviation_layers(3))
   z_crs_top = max(4e-3, gcBRC%aviation_layers(4))

   emis_aviation = 0.0
   srcAviation   = 0.0

   call distribute_aviation_emissions(w_c%delp, rhoa, z_lto_bot, z_lto_top, gcBRC%aviation_lto_src, emis_aviation, i1, i2, j1, j2, km)
   srcAviation = srcAviation + emis_aviation

   call distribute_aviation_emissions(w_c%delp, rhoa, z_cds_bot, z_cds_top, gcBRC%aviation_cds_src, emis_aviation, i1, i2, j1, j2, km)
   srcAviation = srcAviation + emis_aviation

   call distribute_aviation_emissions(w_c%delp, rhoa, z_crs_bot, z_crs_top, gcBRC%aviation_crs_src, emis_aviation, i1, i2, j1, j2, km)
   srcAviation = srcAviation + emis_aviation

!  Determine surface pressure
!  AMS Note: pass this in
!  --------------------------
   ps = 0.0
   do k = 1, km
    ps(i1:i2,j1:j2) = ps(i1:i2,j1:j2) + w_c%delp(i1:i2,j1:j2,k)
   end do

!  Find the pressure of the 100m, 500m, and PBLH altitudes
!  AMS Note: this could be greatly simplified by using ze/zm and having a
!      generic routine from the bottom up with an early exit condition
!  -----------------------------------------------------------------------
   p0 = ps  
   z0(i1:i2,j1:j2) = 0.
   do k = km, 1, -1
    do j = j1, j2
     do i = i1, i2
      p1 = p0(i,j) - w_c%delp(i,j,k)
      dz = w_c%delp(i,j,k)/rhoa(i,j,k)/grav
      z1 = z0(i,j)+dz
      if(z0(i,j) .lt. 100 .and. z1 .ge. 100.) then
       delz = z1-100.
       delp = delz*rhoa(i,j,k)*grav
       p100(i,j) = p1+delp
      endif
      if(z0(i,j) .lt. 500 .and. z1 .ge. 500.) then
       delz = z1-500.
       delp = delz*rhoa(i,j,k)*grav
       p500(i,j) = p1+delp
      endif
      zpbl = max ( pblh(i,j), 100. )
      if(z0(i,j) .lt. zpbl .and. z1 .ge. zpbl ) then
       delz = z1-zpbl
       delp = delz*rhoa(i,j,k)*grav
       pPBL(i,j) = p1+delp
      endif
      p0(i,j) = p1
      z0(i,j) = z1
     end do
    end do
   end do

#if 0
   call pmaxmin ( 'BRC: p100   ', p100,  qmin, qmax, ijl, 1, 1. )
   call pmaxmin ( 'BRC: p500   ', p500,  qmin, qmax, ijl, 1, 1. )
   call pmaxmin ( 'BRC: pPBL   ', pPBLh, qmin, qmax, ijl, 1, 1. )
#endif


!   Limit biomass burning emissions
!   -------------------------------
    allocate(qa_bb_(nbins,i1:i2,j1:j2,km), __STAT__)
    qa_bb_ = 0.0

    p0 = ps
K_LOOP_BB: do k = km, 1, -1

!   First determine emissions for this layer
!   ----------------------------------------
    maxAll = 0.0
    do j = j1, j2
     do i = i1, i2

      p1 = p0(i,j) - w_c%delp(i,j,k)

!     Pressure @ PBL height
!     ---------------------
      fPBL = 0.
      if(p1 .ge. pPBL(i,j)) fPBL = w_c%delp(i,j,k)/(ps(i,j)-pPBL(i,j))
      if(p1 .lt. pPBL(i,j) .and. p0(i,j) .ge. pPBL(i,j)) fPBL = (p0(i,j)-pPBL(i,j))/(ps(i,j)-pPBL(i,j))

!     Sources by class in kg m-2 s-1
!     ------------------------------
      srcBiomass(i,j)  = fPBL * eBiomass * gcBRC%biomass_src(i,j)

      srcHydrophobic(i,j) =     gcBRC%fHydrophobic  * srcBiomass(i,j)
      srcHydrophilic(i,j) = (1.-gcBRC%fHydrophobic) * srcBiomass(i,j)

!     Update pressure of lower level
!     ------------------------------
      p0(i,j) = p1

     end do ! i
    end do  ! j

#if (0)
!   Determine global max/min
!   ------------------------
    call pmaxmin ( 'BRC: Phobic ', srcHydrophobic, qmin, qmax, ijl, 1, 0. )
    maxAll = abs(qmax) + abs(qmin)
    call pmaxmin ( 'BRC: Philic ', srcHydrophilic, qmin, qmax, ijl, 1, 0. )
    maxAll = max ( maxAll, abs(qmax) + abs(qmin) )

!   If emissions are zero at this level (globally), we are done
!   -----------------------------------------------------------
    if ( maxAll .eq. 0.0 ) exit K_LOOP_BB
#endif
    
!   Update concentrations at this layer
!   The "1" element is hydrophobic 
!   The "2" element is hydrophilic
!   -----------------------------------    
    factor = cdt * grav / w_c%delp(:,:,k)

    qa_bb_(1,:,:,k) = factor * srcHydrophobic
    qa_bb_(2,:,:,k) = factor * srcHydrophilic

   end do K_LOOP_BB


    nch   = gcBRC%mie_tables%nch

!   Get the wavelength indices
!   --------------------------
!   Must provide ilam550 for AOT calculation
    ilam550 = 1.
    if(nch .gt. 1) then
     do i = 1, nch
      if ( gcBRC%mie_tables%channels(i) .ge. 5.49e-7 .and. &
           gcBRC%mie_tables%channels(i) .le. 5.51e-7) ilam550 = i
     enddo
    endif

!  Calculate the extinction and/or scattering AOD

   exttau_bb_(i1:i2,j1:j2) = 0.0

   do n = 1, nbins

!     Select the name for species and the index
      qname = trim(w_c%reg%vname(n1+n-1))
      idx = Chem_MieQueryIdx(gcBRC%mie_tables,qname,rc)
      if(rc .ne. 0) call die(myname, 'cannot find proper Mie table index')

      do k = 1, km
       do j = j1, j2
        do i = i1, i2
         call Chem_MieQuery(gcBRC%mie_tables, idx, ilam550, &
              qa_bb_(n,i,j,k)*w_c%delp(i,j,k)/grav, &
              w_c%rh(i,j,k), tau=tau, ssa=ssa)

!        Integrate in the vertical
         exttau_bb_(i,j) = exttau_bb_(i,j) + tau

        enddo
       enddo
      enddo

   enddo  ! nbins


   f_bb_ = 1.0
   cutoff_bb_exttau = (cdt / (24 * 3600.0)) * max_bb_exttau

   do j = j1, j2
    do i = i1, i2
     if (exttau_bb_(i,j) > cutoff_bb_exttau) then
      f_bb_(i,j) = cutoff_bb_exttau / exttau_bb_(i,j)
     end if
    enddo
   enddo
 
   deallocate(qa_bb_, __STAT__)


!  Now update the tracer mixing ratios with the aerosol sources
!  ------------------------------------------------------------
   p0 = ps
K_LOOP: do k = km, 1, -1

!!!    print *, 'BRC_Emissions: getting emissions for layer ', k

!   First determine emissions for this layer
!   ----------------------------------------
    maxAll = 0.0
    do j = j1, j2
     do i = i1, i2

      p1 = p0(i,j) - w_c%delp(i,j,k)

!     Pressure @ 100m
!     ---------------
      f100 = 0.
      if(p1 .ge. p100(i,j)) f100 = w_c%delp(i,j,k)/(ps(i,j)-p100(i,j))
      if(p1 .lt. p100(i,j) .and. p0(i,j) .ge. p100(i,j)) &
       f100 = (p0(i,j)-p100(i,j))/(ps(i,j)-p100(i,j))

!     Pressure @ 500m
!     ---------------
      f500 = 0.
      if ( p0(i,j) .ge. p100(i,j) .and. p1 .lt. p100(i,j) .and. p1 .ge. p500(i,j)) &
       f500 = (p100(i,j)-p1)/(p100(i,j)-p500(i,j))
      if(p0(i,j) .lt. p100(i,j) .and. p1 .ge. p500(i,j)) &
       f500 = w_c%delp(i,j,k)/(p100(i,j)-p500(i,j))
      if(p0(i,j) .ge. p500(i,j) .and. p1 .lt. p500(i,j)) &
       f500 = (p0(i,j)-p500(i,j))/(p100(i,j)-p500(i,j))

!     Pressure @ PBL height
!     ---------------------
      fPBL = 0.
      if(p1 .ge. pPBL(i,j)) fPBL = w_c%delp(i,j,k)/(ps(i,j)-pPBL(i,j))
      if(p1 .lt. pPBL(i,j) .and. p0(i,j) .ge. pPBL(i,j)) &
       fPBL = (p0(i,j)-pPBL(i,j))/(ps(i,j)-pPBL(i,j))

!     Terpene is tree-top emission; only add in bottom layer
!     ------------------------------------------------------
      if ( k .eq. km ) then
           fBot = 1.0
      else
           fBot = 0.0
      end if

!     Sources by class in kg m-2 s-1
!     ------------------------------
      srcBiofuel(i,j)  = f100 * eBiofuel * gcBRC%biofuel_src(i,j)
      srcAnthro(i,j)   = f100 * eAnthro  * gcBRC%ebrcantt1_src(i,j) &
                       + f500 * eAnthro  * gcBRC%ebrcantt2_src(i,j) &
                       + f100 * eAnthro  * gcBRC%brc_ship_src(i,j) &
                       +        eAnthro  * srcAviation(i,j,k)
      srcBiomass(i,j)  = fPBL * eBiomass * gcBRC%biomass_src(i,j) * f_bb_(i,j)
      srcBiogenic(i,j) = fBot * eTerpene * gcBRC%terpene_src(i,j)

      srcTmp = srcBiofuel(i,j) + srcAnthro(i,j) + srcBiomass(i,j)

      srcHydrophobic(i,j) =     gcBRC%fHydrophobic  * srcTmp
      srcHydrophilic(i,j) = (1.-gcBRC%fHydrophobic) * srcTmp + srcBiogenic(i,j)

!     Update pressure of lower level
!     ------------------------------
      p0(i,j) = p1

     end do ! i
    end do  ! j

#if (0)
!   Determine global max/min
!   ------------------------
    call pmaxmin ( 'BRC: Phobic ', srcHydrophobic, qmin, qmax, ijl, 1, 0. )
    maxAll = abs(qmax) + abs(qmin)
    call pmaxmin ( 'BRC: Philic ', srcHydrophilic, qmin, qmax, ijl, 1, 0. )
    maxAll = max ( maxAll, abs(qmax) + abs(qmin) )

!   If emissions are zero at this level (globally), we are done
!   -----------------------------------------------------------
    if ( maxAll .eq. 0.0 ) exit K_LOOP
#endif

!   Update concentrations at this layer
!   The "1" element is hydrophobic 
!   The "2" element is hydrophilic
!   -----------------------------------    
    factor = cdt * grav / w_c%delp(:,:,k)

    w_c%qa(n1)%data3d(:,:,k) = w_c%qa(n1)%data3d(:,:,k) & 
                             + factor * srcHydrophobic 

    w_c%qa(n2)%data3d(:,:,k) = w_c%qa(n2)%data3d(:,:,k) & 
                             + factor * srcHydrophilic

!   Fill in diagnostics if requested
!   --------------------------------
    if ( associated(BRC_emis(1)%data2d)) &
                    BRC_emis(1)%data2d = BRC_emis(1)%data2d + srcHydrophobic

    if ( associated(BRC_emis(2)%data2d)) &
                    BRC_emis(2)%data2d = BRC_emis(2)%data2d + srcHydrophilic

    if ( associated(BRC_emisBF%data2d)) &
                    BRC_emisBF%data2d  = BRC_emisBF%data2d  + srcBiofuel

    if ( associated(BRC_emisBB%data2d)) &
                    BRC_emisBB%data2d  = BRC_emisBB%data2d  + srcBiomass

    if ( associated(BRC_emisAN%data2d)) &
                    BRC_emisAN%data2d  = BRC_emisAN%data2d  + srcAnthro

   if ( associated(BRC_emisBG%data2d)) &
                   BRC_emisBG%data2d   = BRC_emisBG%data2d  + srcBiogenic 

   end do K_LOOP

!  Distribute pointwise sources if requested
!  -----------------------------------------
   if( gcBRC%doing_point_emissions .and. gcBRC%nPts > 0) then

!    Get indices for point emissions
!    -------------------------------
     allocate(iPoint(gcBRC%nPts), jPoint(gcBRC%nPts), stat=ios)

     call MAPL_GetHorzIJIndex(gcBRC%nPts, iPoint, jPoint, &
                              grid = w_c%grid_esmf,      &
                              lon  = gcBRC%vLon/radToDeg, &
                              lat  = gcBRC%vLat/radToDeg, &
                              rc   = rc)

     if ( rc /= 0 ) call die(myname,'cannot get indices for point emissions')

     do ii = 1, gcBRC%nPts
      i = iPoint(ii)
      j = jPoint(ii)
      if( i<1 .OR. j<1 )              cycle    ! point emission not in this sub-domain
!      if( gcBRC%regionMask(i,j) == 0 ) cycle    ! masked by region mask
      
!     Emissions not occurring in current time step
!     --------------------------------------------
      if(nhms < gcBRC%vStart(ii) .or. nhms >= gcBRC%vEnd(ii)) cycle

      call Chem_UtilDistributePointEmissions(hghte(i,j,:), &
                                             gcBRC%vBase(ii), gcBRC%vTop(ii), gcBRC%vEmis(ii), &
                                             point_column_emissions, km)
      w_c%qa(n1)%data3d(i,j,:) = w_c%qa(n1)%data3d(i,j,:) & 
         + gcBRC%fHydrophobic * cdt * grav / w_c%delp(i,j,:) &
                             * point_column_emissions / w_c%grid%cell_area(i,j)
      w_c%qa(n2)%data3d(i,j,:) = w_c%qa(n2)%data3d(i,j,:) & 
         + (1-gcBRC%fHydrophobic) * cdt * grav / w_c%delp(i,j,:) &
                                 * point_column_emissions / w_c%grid%cell_area(i,j)

     enddo
     deallocate(iPoint, jPoint, stat=ios)
   endif


   rc = 0

   end subroutine BRC_Emission

   subroutine distribute_aviation_emissions(delp, rhoa, z_bot, z_top, emissions_layer, emissions, i1, i2, j1, j2, km)

    implicit none

    integer, intent(in) :: i1, i2, j1, j2, km

    real, dimension(:,:,:), intent(in) :: delp
    real, dimension(:,:,:), intent(in) :: rhoa
    real, dimension(:,:),   intent(in) :: emissions_layer
    real, intent(in)                   :: z_bot
    real, intent(in)                   :: z_top
    real, dimension(:,:,:), intent(out):: emissions
    
!   local
    integer :: i, j, k
    integer :: k_bot, k_top
    real    :: z_
    real, dimension(km) :: z, dz, w_
    
    do j = j1, j2
        do i = i1, i2
            ! find level height
            z = 0.0
            z_= 0.0 

            do k = km, 1, -1
                dz(k) = delp(i,j,k)/rhoa(i,j,k)/grav
                z_    = z_ + dz(k)
                z(k)  = z_
            end do

            ! find the bottom level
            do k = km, 1, -1
                if (z(k) >= z_bot) then
                    k_bot = k
                    exit
                end if
            end do
            
            ! find the top level
            do k = k_bot, 1, -1
                if (z(k) >= z_top) then
                    k_top = k
                    exit
                end if
            end do

            ! find the weights
            w_ = 0

!           if (k_top > k_bot) then
!               need to bail - something went wrong here
!           end if

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
           
            ! distribute emissions in the vertical 
            emissions(i,j,:) = (w_ / sum(w_)) * emissions_layer(i,j)
        end do 
    end do

    end subroutine distribute_aviation_emissions

 end subroutine BRC_GridCompRun1_



!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  BRC_GridCompRun2_ --- The Chem Driver, run phase 2 
!
! !INTERFACE:
!

   subroutine BRC_GridCompRun2_ ( gcBRC, w_c, impChem, expChem, ggState, &
                                  nymd, nhms, cdt, rc )

! !USES:

  implicit NONE

! !INPUT/OUTPUT PARAMETERS:

   type(BRC_GridComp1), intent(inout) :: gcBRC   ! Grid Component
   type(Chem_Bundle), intent(inout)  :: w_c    ! Chemical tracer fields
   type(MAPL_MetaComp), intent(inout) :: ggState

! !INPUT PARAMETERS:

   type(ESMF_State), intent(inout) :: impChem  ! Import State
   integer, intent(in) :: nymd, nhms           ! time
   real, intent(in) :: cdt                     ! chemistry timestep (secs)

! !OUTPUT PARAMETERS:

   type(ESMF_State), intent(inout) :: expChem  ! Export State
   integer, intent(out) :: rc                  ! Error return code:
                                               !  0 - all is well
                                               !  1 -
 
! !DESCRIPTION: This routine implements the so-called BRC Driver. That 
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

   character(len=*), parameter :: myname = 'BRC_GridCompRun2_'
   character(len=*), parameter :: Iam = myname

   integer :: i1, i2, im, j1, j2, jm, nbins, n1, n2, km, n, ios
   integer :: i, j, k, ijl, ijkl, ijk1l
   real :: qmax, qmin
   real :: qUpdate, delq
   real, pointer :: dqa(:,:), drydepositionfrequency(:,:)
   type(Chem_Array), pointer :: fluxout
   logical :: KIN

!  Input fields from fvGCM
!  -----------------------
   real, pointer, dimension(:,:)   :: frlake, frocean, frseaice, &
                                      oro, u10m, v10m, &
                                      ustar, precc, precl,                &
                                      pblh, shflux, z0h, hsurf
   real, pointer, dimension(:,:,:) :: tmpu, rhoa, u, v, hghte, ple
   real, pointer, dimension(:,:,:) :: pfllsan, pfilsan


!  Additional needs for GOCART convective diagnostic
   real, pointer, dimension(:,:,:)       ::  cmfmc, qlcn, qicn, dtrain
   real, pointer, dimension(:,:)         ::  area
   real*8, allocatable, dimension(:,:,:) ::  cmfmc_, qccu_, dtrain_, &
                                             airmass_, airmol_, vud_, &
                                             delz_, delp_, tmpu_, ple_
   real*8, allocatable                   ::  tc_(:,:,:,:), bcnv_(:,:,:)
   real*8, allocatable                   ::  area_(:,:), frlake_(:,:), &
                                             frocean_(:,:), frseaice_(:,:)
   integer*4                             ::  icdt

   real, pointer :: BRC_radius(:), BRC_rhop(:)
   integer       :: rhFlag
  
   real, pointer :: var3d(:,:,:) => null() 


#define EXPORT     expChem
#define iNAME      TRIM(gcBRC%iname)

#define ptrBRCWT       BRC_wet
#define ptrBRCSV       BRC_conv
#define ptrBRCEM       BRC_emis
#define ptrBRCDP       BRC_dep
#define ptrBRCSD       BRC_set

#define ptrBRCMASS     BRC_mass
#define ptrBRCEMAN     BRC_emisAN
#define ptrBRCEMBB     BRC_emisBB
#define ptrBRCEMBF     BRC_emisBF
#define ptrBRCEMBG     BRC_emisBG
#define ptrBRCPSOA     BRC_pSOA
#define ptrBRCHYPHIL   BRC_toHydrophilic
#define ptrBRCSMASS    BRC_sfcmass
#define ptrBRCCMASS    BRC_colmass
#define ptrBRCEXTTAU   BRC_exttau
#define ptrBRCSCATAU   BRC_scatau
#define ptrBRCCONC     BRC_conc
#define ptrBRCEXTCOEF  BRC_extcoef
#define ptrBRCSCACOEF  BRC_scacoef
#define ptrBRCANGSTR   BRC_angstrom
#define ptrBRCFLUXU    BRC_fluxu
#define ptrBRCFLUXV    BRC_fluxv

   integer :: STATUS

#include "BRC_GetPointer___.h"

   call MAPL_TimerOn(ggState, '-BRC_RUN2')

!  Initialize local variables
!  --------------------------
   rc = 0
   i1 = w_c%grid%i1; i2 = w_c%grid%i2; im = w_c%grid%im
   j1 = w_c%grid%j1; j2 = w_c%grid%j2; jm = w_c%grid%jm

   km    = w_c%grid%km
   nbins = w_c%reg%n_BRC
   n1    = w_c%reg%i_BRC
   n2    = w_c%reg%j_BRC

   ijl   = ( i2 - i1 + 1 ) * ( j2 - j1 + 1 )
   ijkl  = ijl * km
   ijk1l = ijl * (km+1)


#ifdef DEBUG
   do n = n1, n2
      call pmaxmin('BRC: q_beg', w_c%qa(n)%data3d(i1:i2,j1:j2,1:km), qmin, qmax, &
                   ijl, km, 1. )
   end do
#endif



!  Get 2D Imports
!  --------------
   call MAPL_GetPointer ( impChem, frlake,   'FRLAKE',   __RC__ )
   call MAPL_GetPointer ( impChem, oro,      'LWI',      __RC__ )
   call MAPL_GetPointer ( impChem, u10m,     'U10M',     __RC__ )
   call MAPL_GetPointer ( impChem, v10m,     'V10M',     __RC__ )
   call MAPL_GetPointer ( impChem, ustar,    'USTAR',    __RC__ )
   call MAPL_GetPointer ( impChem, precc,    'CN_PRCP',  __RC__ )
   call MAPL_GetPointer ( impChem, precl,    'NCN_PRCP', __RC__ )
   call MAPL_GetPointer ( impChem, pblh,     'ZPBL',     __RC__ )
   call MAPL_GetPointer ( impChem, shflux,   'SH',       __RC__ )
   call MAPL_GetPointer ( impChem, z0h,      'Z0H',      __RC__ )
   call MAPL_GetPointer ( impChem, area,     'AREA',     __RC__ )
   call MAPL_GetPointer ( impChem, frocean,  'FROCEAN',  __RC__ )
   call MAPL_GetPointer ( impChem, frseaice, 'FRACI',    __RC__ )

!  Get 3D Imports
!  --------------
   call MAPL_GetPointer ( impChem, tmpu,     'T',        __RC__ )
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

!  Unlike GEOS-4 hghte is defined for km+1
!  ---------------------------------------
   hsurf => hghte(i1:i2,j1:j2,km) ! Recall: GEOS-5 has edges with k in [0,km]
    
  

#ifdef DEBUG

   call pmaxmin('BRC: frlake     ', frlake  , qmin, qmax, ijl,1, 1. )
   call pmaxmin('BRC: frocean    ', frocean , qmin, qmax, ijl,1, 1. )
   call pmaxmin('BRC: frseaice   ', frseaice, qmin, qmax, ijl,1, 1. )
   call pmaxmin('BRC: area       ', area    , qmin, qmax, ijl,1, 1. )
   call pmaxmin('BRC: oro        ', oro     , qmin, qmax, ijl,1, 1. )
   call pmaxmin('BRC: u10m       ', u10m    , qmin, qmax, ijl,1, 1. )
   call pmaxmin('BRC: v10m       ', v10m    , qmin, qmax, ijl,1, 1. )
   call pmaxmin('BRC: ustar      ', ustar   , qmin, qmax, ijl,1, 1. )
   call pmaxmin('BRC: precc      ', precc   , qmin, qmax, ijl,1, 1. )
   call pmaxmin('BRC: precl      ', precl   , qmin, qmax, ijl,1, 1. )
   call pmaxmin('BRC: pblh       ', pblh    , qmin, qmax, ijl,1, 1. )
   call pmaxmin('BRC: shflux     ', shflux  , qmin, qmax, ijl,1, 1. )
   call pmaxmin('BRC: z0h        ', z0h     , qmin, qmax, ijl,1, 1. )
   call pmaxmin('BRC: hsurf      ', hsurf   , qmin, qmax, ijl,1, 1. )

   call pmaxmin('BRC: tmpu       ', tmpu    , qmin, qmax, ijkl,1, 1. )
   call pmaxmin('BRC: rhoa       ', rhoa    , qmin, qmax, ijkl,1, 1. )
   call pmaxmin('BRC: u          ', u       , qmin, qmax, ijkl,1, 1. )
   call pmaxmin('BRC: v          ', v       , qmin, qmax, ijkl,1, 1. )
   call pmaxmin('BRC: hghte      ', hghte   , qmin, qmax, ijkl,1, 1. )
   call pmaxmin('BRC: qlcn       ', qlcn    , qmin, qmax, ijkl,1, 1. )
   call pmaxmin('BRC: qicn       ', qicn    , qmin, qmax, ijkl,1, 1. )
   call pmaxmin('BRC: cmfmc      ', cmfmc   , qmin, qmax, ijkl,1, 1. )
   call pmaxmin('BRC: dtrain     ', dtrain  , qmin, qmax, ijkl,1, 1. )
   call pmaxmin('BRC: pfllsan    ', pfllsan , qmin, qmax, ijk1l,1, 1. )
   call pmaxmin('BRC: pfilsan    ', pfilsan , qmin, qmax, ijk1l,1, 1. )

#endif

RUN_ALARM: if (gcBRC%run_alarm) then

   allocate( fluxout )
   allocate( fluxout%data2d(i1:i2,j1:j2), dqa(i1:i2,j1:j2), &
             drydepositionfrequency(i1:i2,j1:j2), stat=STATUS)
   VERIFY_(STATUS)


!  SOA production from oxidation of biomass burning VOC
   call MAPL_GetPointer(impChem, var3d, 'pSOA_BIOB_VOC'//iNAME, __RC__)
   gcBRC%psoa_biob_voc = var3d

   where( 1.01 * gcBRC%psoa_biob_voc .gt. undefval) gcBRC%psoa_biob_voc = 0.0 


!  Add on SOA from Anthropogenic VOC oxidation
!  -------------------------------------------
   w_c%qa(n2)%data3d = w_c%qa(n2)%data3d + cdt*gcBRC%psoa_biob_voc/rhoa  ! hydrophilic

   if ( associated(BRC_pSOA%data2d)) &
       BRC_pSOA%data2d = sum(gcBRC%psoa_biob_voc*w_c%delp/rhoa/grav, 3)


!  Ad Hoc transfer of hydrophobic to hydrophilic aerosols
!  Following Chin's parameterization, the rate constant is
!  k = 4.63e-6 s-1 (.4 day-1; e-folding time = 2.5 days)
   if(associated(BRC_toHydrophilic%data2d)) &
     BRC_toHydrophilic%data2d(i1:i2,j1:j2) = 0.0

   do k = 1, km
    do j = j1, j2
     do i = i1, i2
      qUpdate = w_c%qa(n1)%data3d(i,j,k)*exp(-4.63e-6*cdt)
      qUpdate = max(qUpdate,1.e-32)
      delq = max(0.,w_c%qa(n1)%data3d(i,j,k)-qUpdate)
      w_c%qa(n1)%data3d(i,j,k) = qUpdate
      w_c%qa(n2)%data3d(i,j,k) = w_c%qa(n2)%data3d(i,j,k)+delq
      if(associated(BRC_toHydrophilic%data2d)) &
       BRC_toHydrophilic%data2d(i,j) = BRC_toHydrophilic%data2d(i,j) &
        + delq*w_c%delp(i,j,k)/grav/cdt
     end do
    end do
   end do

!  BRC Settling
!  -----------
   call MAPL_TimerOn(ggState, '--BRC_SETTLING')

   allocate( BRC_radius(nbins), BRC_rhop(nbins) )
   BRC_radius(:) = 0.35e-6  ! radius for settling [m]
   BRC_rhop(:)   = 1800.    ! density for setting [kg m-3]
   rhFlag       = 0        ! settle like dry particles
   call Chem_Settling ( i1, i2, j1, j2, km, n1, n2, nbins, rhFlag, &
                        BRC_radius, BRC_rhop, cdt, w_c, tmpu, rhoa, hsurf,    &
                        hghte, BRC_set, rc )
   deallocate( BRC_radius, BRC_rhop)

   call MAPL_TimerOff(ggState, '--BRC_SETTLING')

!  BRC Deposition
!  --------------
   call MAPL_TimerOn(ggState, '--BRC_DRY_DEPOSITION')

   drydepositionfrequency = 0.
   call DryDepositionGOCART( i1, i2, j1, j2, km, &
                             tmpu, rhoa, hghte, oro, ustar, &
                             pblh, shflux, z0h, drydepositionfrequency, rc )
    
   do n = 1, nbins
    dqa = 0.
    dqa = max(0.0, w_c%qa(n1+n-1)%data3d(:,:,km)*(1.-exp(-drydepositionfrequency*cdt)))
    w_c%qa(n1+n-1)%data3d(:,:,km) = &
            w_c%qa(n1+n-1)%data3d(:,:,km) - dqa
    if( associated(BRC_dep(n)%data2d) ) &
     BRC_dep(n)%data2d = dqa*w_c%delp(:,:,km)/grav/cdt
   end do

   call MAPL_TimerOff(ggState, '--BRC_DRY_DEPOSITION')

#ifdef DEBUG
   do n = n1, n2
      call pmaxmin('BRC: q_dry', w_c%qa(n)%data3d(i1:i2,j1:j2,1:km), qmin, qmax, &
                    ijl, km, 1. )
   end do
#endif


!  Organic Carbon Large-scale Wet Removal
!  --------------------------------------
   call MAPL_TimerOn(ggState, '--BRC_WET_LS')

!  Hydrophobic mode (first tracer) is not removed
   if(associated(BRC_wet(1)%data2d)) BRC_wet(1)%data2d = 0.
!  Hydrophilic mode (second tracer) is removed
   KIN = .TRUE.
   do n = nbins, nbins
    w_c%qa(n1+n-1)%fwet = 1.
    call WetRemovalGOCART(i1, i2, j1, j2, km, n1+n-1, n1+n-1, cdt, 'OC', KIN, &
                          w_c%qa, ple, tmpu, rhoa, pfllsan, pfilsan, &
                          precc, precl, fluxout, rc )
    if(associated(BRC_wet(n)%data2d)) BRC_wet(n)%data2d = fluxout%data2d
   end do

   call MAPL_TimerOff(ggState, '--BRC_WET_LS')

#ifdef DEBUG
   do n = n1, n2
      call pmaxmin('BRC: q_wet', w_c%qa(n)%data3d(i1:i2,j1:j2,1:km), qmin, qmax, &
                    ijl, km, 1. )
   end do
#endif

!  Organic Carbon Convective-scale Mixing and Wet Removal
!  ------------------------------------------------------
   call MAPL_TimerOn(ggState, '--BRC_WET_CV')

   KIN = .TRUE.
   icdt = cdt
   allocate(cmfmc_(i1:i2,j1:j2,km+1), qccu_(i1:i2,j1:j2,km), &
            dtrain_(i1:i2,j1:j2,km), airmass_(i1:i2,j1:j2,km), &
            delz_(i1:i2,j1:j2,km), vud_(i1:i2,j1:j2,km), &
            tc_(i1:i2,j1:j2,km,n1:n2), delp_(i1:i2,j1:j2,km), &
            airmol_(i1:i2,j1:j2,km), tmpu_(i1:i2,j1:j2,km), &
            bcnv_(i1:i2,j1:j2,n1:n2), ple_(i1:i2,j1:j2,km+1), &
            area_(i1:i2,j1:j2), frlake_(i1:i2,j1:j2), &
            frocean_(i1:i2,j1:j2), frseaice_(i1:i2,j1:j2), __STAT__ )

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
   call set_vud(i1, i2, j1, j2, km, frlake_, frocean_, frseaice_, cmfmc_, qccu_, &
                airmass_, delz_, area_, vud_)
   call convection(i1, i2, j1, j2, km, n1, n2, icdt, 'OC', kin, &
                   tc_, cmfmc_, dtrain_, area_, delz_, delp_, vud_, &
                   airmass_, airmol_, tmpu_, ple_, &
                   bcnv_)

!  Return adjusted tracer to mixing ratio
   do n = n1, n2
    do k = 1, km
     w_c%qa(n)%data3d(:,:,km-k+1) = tc_(:,:,k,n)
    enddo
   enddo

!  Note GOCART returns bcnv_ as negative, recast for my diagnostic
   if(associated(BRC_conv(1)%data2d)) BRC_conv(1)%data2d = 0.0
   if(associated(BRC_conv(2)%data2d)) BRC_conv(2)%data2d = -bcnv_(:,:,n2)/area_/icdt

!  Clean up
!  --------
   deallocate(cmfmc_, qccu_, dtrain_, tc_, airmass_, &
              delz_, vud_, delp_, airmol_, tmpu_, bcnv_, ple_, &
              area_, frlake_, frocean_, frseaice_, __STAT__ )

   deallocate(fluxout%data2d)
   deallocate(fluxout, dqa, drydepositionfrequency, stat=ios )

   call MAPL_TimerOff(ggState, '--BRC_WET_CV')

   end if RUN_ALARM


!  Compute the desired output diagnostics here
!  Ideally this will go where chemout is called in fvgcm.F since that
!  will reflect the distributions after transport, etc.
!  ------------------------------------------------------------------
   call MAPL_TimerOn(ggState, '--BRC_DIAGNOSTICS')

   call BRC_Compute_Diags(i1, i2, j1, j2, km, nbins, gcBRC, w_c, tmpu, rhoa, u, v, &
                         BRC_sfcmass, BRC_colmass, BRC_mass, BRC_exttau, &
                         BRC_scatau, BRC_conc, BRC_extcoef, BRC_scacoef, BRC_angstrom, &
                         BRC_fluxu, BRC_fluxv, rc)

   call MAPL_TimerOff(ggState, '--BRC_DIAGNOSTICS')

   call MAPL_TimerOff(ggState, '-BRC_RUN2')

   return

CONTAINS

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  BRC_Compute_Diags - Calculate dust 2D diagnostics
!
! !INTERFACE:
!

   subroutine BRC_Compute_Diags ( i1, i2, j1, j2, km, nbins, gcBRC, w_c, tmpu, rhoa, u, v, &
                                  sfcmass, colmass, mass, exttau, scatau, &
                                  conc, extcoef, scacoef, angstrom, fluxu, fluxv, rc )

! !USES:

  implicit NONE

! !INPUT PARAMETERS:
   integer, intent(in) :: i1, i2, j1, j2, km, nbins
   type(BRC_GridComp1), intent(inout):: gcBRC     ! BRC Grid Component
   type(Chem_Bundle), intent(in)   :: w_c      ! Chem Bundle
   real, pointer, dimension(:,:,:) :: tmpu     ! temperature [K]
   real, pointer, dimension(:,:,:) :: rhoa     ! air density [kg m-3]
   real, pointer, dimension(:,:,:) :: u        ! east-west wind [m s-1]
   real, pointer, dimension(:,:,:) :: v        ! north-south wind [m s-1]

! !OUTPUT PARAMETERS:
   type(Chem_Array), intent(inout)  :: sfcmass  ! sfc mass concentration kg/m3
   type(Chem_Array), intent(inout)  :: colmass  ! col mass density kg/m2
   type(Chem_Array), intent(inout)  :: mass     ! 3d mass mixing ratio kg/kg
   type(Chem_Array), intent(inout)  :: exttau   ! ext. AOT at 550 nm
   type(Chem_Array), intent(inout)  :: scatau   ! sct. AOT at 550 nm
   type(Chem_Array), intent(inout)  :: conc     ! 3d mass concentration, kg/m3
   type(Chem_Array), intent(inout)  :: extcoef  ! 3d ext. coefficient, 1/m
   type(Chem_Array), intent(inout)  :: scacoef  ! 3d scat.coefficient, 1/m
   type(Chem_Array), intent(inout)  :: angstrom ! 470-870 nm Angstrom parameter
   type(Chem_Array), intent(inout)  :: fluxu    ! Column mass flux in x direction
   type(Chem_Array), intent(inout)  :: fluxv    ! Column mass flux in y direction
   integer, intent(out)             :: rc       ! Error return code:
                                                !  0 - all is well
                                                !  1 - 

! !DESCRIPTION: Calculates some simple 2d diagnostics from the BRC fields
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
   character(len=*), parameter :: myname = 'BRC_Compute_Diags'
   integer :: i, j, k, n, n1, n2, ios, nch, idx
   real :: tau, ssa
   character(len=255) :: qname
   real, dimension(i1:i2,j1:j2) :: tau470, tau870
   real    :: ilam550, ilam470, ilam870
   logical :: do_angstrom


!  Initialize local variables
!  --------------------------
   n1  = w_c%reg%i_BRC
   n2  = w_c%reg%j_BRC
   nch   = gcBRC%mie_tables%nch

!  Get the wavelength indices
!  --------------------------
!  Must provide ilam550 for AOT calculation
   ilam550 = 1.
   ilam470 = 0.
   ilam870 = 0.
   if(nch .gt. 1) then
    do i = 1, nch
     if ( gcBRC%mie_tables%channels(i) .ge. 5.49e-7 .and. &
          gcBRC%mie_tables%channels(i) .le. 5.51e-7) ilam550 = i
     if ( gcBRC%mie_tables%channels(i) .ge. 4.69e-7 .and. &
          gcBRC%mie_tables%channels(i) .le. 4.71e-7) ilam470 = i
     if ( gcBRC%mie_tables%channels(i) .ge. 8.69e-7 .and. &
          gcBRC%mie_tables%channels(i) .le. 8.71e-7) ilam870 = i
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
   if( associated(sfcmass%data2d) ) then
      sfcmass%data2d(i1:i2,j1:j2) = 0.
      do n = 1, nbins
         sfcmass%data2d(i1:i2,j1:j2) &
              =   sfcmass%data2d(i1:i2,j1:j2) &
              + w_c%qa(n+n1-1)%data3d(i1:i2,j1:j2,km)*rhoa(i1:i2,j1:j2,km)
      end do
   endif

!  Calculate the dust column loading
   if( associated(colmass%data2d) ) then
      colmass%data2d(i1:i2,j1:j2) = 0.
      do n = 1, nbins
       do k = 1, km
        colmass%data2d(i1:i2,j1:j2) &
         =   colmass%data2d(i1:i2,j1:j2) &
           + w_c%qa(n+n1-1)%data3d(i1:i2,j1:j2,k)*w_c%delp(i1:i2,j1:j2,k)/grav
       end do
      end do
   endif

!  Calculate the total mass concentration
   if( associated(conc%data3d) ) then
      conc%data3d(i1:i2,j1:j2,1:km) = 0.
      do n = 1, nbins
       conc%data3d(i1:i2,j1:j2,1:km) &
         =   conc%data3d(i1:i2,j1:j2,1:km) &
           + w_c%qa(n+n1-1)%data3d(i1:i2,j1:j2,1:km)*rhoa(i1:i2,j1:j2,1:km)
      end do
   endif

!  Calculate the total mass mixing ratio
   if( associated(mass%data3d) ) then
      mass%data3d(i1:i2,j1:j2,1:km) = 0.
      do n = 1, nbins
       mass%data3d(i1:i2,j1:j2,1:km) &
         =   mass%data3d(i1:i2,j1:j2,1:km) &
           + w_c%qa(n+n1-1)%data3d(i1:i2,j1:j2,1:km)
      end do
   endif

!  Calculate the column mass flux in x direction
   if( associated(fluxu%data2d) ) then
      fluxu%data2d(i1:i2,j1:j2) = 0.
      do n = 1, nbins
       do k = 1, km
        fluxu%data2d(i1:i2,j1:j2) &
         =   fluxu%data2d(i1:i2,j1:j2) &
           + w_c%qa(n+n1-1)%data3d(i1:i2,j1:j2,k)*w_c%delp(i1:i2,j1:j2,k)/grav*u(i1:i2,j1:j2,k)
       end do
      end do
   endif   
   
!  Calculate the column mass flux in y direction
   if( associated(fluxv%data2d) ) then
      fluxv%data2d(i1:i2,j1:j2) = 0.
      do n = 1, nbins
       do k = 1, km
        fluxv%data2d(i1:i2,j1:j2) &
         =   fluxv%data2d(i1:i2,j1:j2) &
           + w_c%qa(n+n1-1)%data3d(i1:i2,j1:j2,k)*w_c%delp(i1:i2,j1:j2,k)/grav*v(i1:i2,j1:j2,k)
       end do
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

      if( associated(extcoef%data3d)) then 
       extcoef%data3d(i1:i2,j1:j2,1:km) = 0.
      endif
      if( associated(scacoef%data3d)) then
       scacoef%data3d(i1:i2,j1:j2,1:km) = 0.
      endif 

      do n = 1, nbins

!      Select the name for species and the index
       qname = trim(w_c%reg%vname(n1+n-1))
       idx = Chem_MieQueryIdx(gcBRC%mie_tables,qname,rc)
       if(rc .ne. 0) call die(myname, 'cannot find proper Mie table index')

       do k = 1, km
        do j = j1, j2
         do i = i1, i2
          call Chem_MieQuery(gcBRC%mie_tables, idx, ilam550, &
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

      do n = 1, nbins

!      Select the name for species
       qname = trim(w_c%reg%vname(n+n1-1))
       idx = Chem_MieQueryIdx(gcBRC%mie_tables,qname,rc)
       if(rc .ne. 0) call die(myname, 'cannot find proper Mie table index')

       do k = 1, km
        do j = j1, j2
         do i = i1, i2

          call Chem_MieQuery(gcBRC%mie_tables, idx, ilam470, &
              w_c%qa(n+n1-1)%data3d(i,j,k)*w_c%delp(i,j,k)/grav, &
              w_c%rh(i,j,k), tau=tau)
          tau470(i,j) = tau470(i,j) + tau

          call Chem_MieQuery(gcBRC%mie_tables, idx, ilam870, &
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


   rc = 0

   end subroutine BRC_Compute_Diags

 end subroutine BRC_GridCompRun2_



!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  BRC_GridCompFinalize --- The Chem Driver 
!
! !INTERFACE:
!

   subroutine BRC_GridCompFinalize1_ ( gcBRC, w_c, impChem, expChem, ggState, &
                                       nymd, nhms, cdt, rc )

! !USES:

  implicit NONE

! !INPUT/OUTPUT PARAMETERS:

   type(BRC_GridComp1), intent(inout) :: gcBRC   ! Grid Component

! !INPUT PARAMETERS:

   type(Chem_Bundle), intent(in)  :: w_c      ! Chemical tracer fields   
   integer, intent(in) :: nymd, nhms          ! time
   real,    intent(in) :: cdt                 ! chemical timestep (secs)


! !OUTPUT PARAMETERS:

   type(ESMF_State), intent(inout) :: impChem   ! Import State
   type(ESMF_State), intent(inout) :: expChem   ! Import State
   type(MAPL_MetaComp), intent(inout) :: ggState
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
   character(len=*), parameter :: myname = 'BRC_GridCompFinalize'

!  If initialized pointwise emissions from daily tables, clean-up
   if(associated(gcBRC%vLat))    deallocate(gcBRC%vLat, stat=ios)
   if(associated(gcBRC%vLon))    deallocate(gcBRC%vLon, stat=ios)
   if(associated(gcBRC%vEmis))   deallocate(gcBRC%vEmis, stat=ios)
   if(associated(gcBRC%vBase))   deallocate(gcBRC%vBase, stat=ios)
   if(associated(gcBRC%vTop))    deallocate(gcBRC%vTop, stat=ios)
   if(associated(gcBRC%vStart))  deallocate(gcBRC%vStart, stat=ios)
   if(associated(gcBRC%vEnd))    deallocate(gcBRC%vEnd, stat=ios)

   rc=0
   return

 end subroutine BRC_GridCompFinalize1_

 end module BRC_GridCompMod


!-----------------------------------------------------------------------

!                     Single Instance Wrapper

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  BRC_SingleInstance_ --- Runs single instance of method
!
! !INTERFACE:
!
  subroutine BRC_SingleInstance_ ( Method_, instance, &
                                   gcBRC, w_c, impChem, expChem, ggState, &
                                   nymd, nhms, cdt, rc )

! !USES:

  Use BRC_GridCompMod
  Use ESMF
  Use MAPL
  Use Chem_Mod 

  IMPLICIT NONE

! !INPUT PARAMETERS:

!  Input "function pointer"
!  -----------------------
   interface 
     subroutine Method_ (gc, w, imp, exp, state, ymd, hms, dt, rcode )
       Use BRC_GridCompMod
       Use ESMF
       Use MAPL
       Use Chem_Mod 
       type(BRC_GridComp1),  intent(inout)  :: gc
       type(Chem_Bundle),   intent(in)     :: w
       type(ESMF_State),    intent(inout)  :: imp
       type(ESMF_State),    intent(inout)  :: exp
       type(MAPL_MetaComp), intent(inout)  :: state
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

   TYPE(BRC_GridComp1), INTENT(INOUT) :: gcBRC    ! Grid Component
   TYPE(ESMF_State), INTENT(INOUT)  :: impChem  ! Import State
   TYPE(ESMF_State), INTENT(INOUT)  :: expChem  ! Export State
   TYPE(MAPL_MetaComp), intent(INOUT) :: ggState
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

  integer n_BRC, i_BRC, j_BRC
  character(len=255) :: i_qname, j_qname

! Save overall CO indices
! -----------------------
  n_BRC = w_c%reg%n_BRC
  i_BRC = w_c%reg%i_BRC
  j_BRC = w_c%reg%j_BRC

! Save the name of the variables in this instance
! -----------------------------------------------
  i_qname = trim(w_c%reg%vname(i_BRC + 2*(instance - 1)))
  j_qname = trim(w_c%reg%vname(i_BRC + 2*(instance - 1) + 1))
  
! Customize indices for this particular instance
! ----------------------------------------------
  w_c%reg%n_BRC = 2
  w_c%reg%i_BRC = i_BRC + 2*(instance - 1)
  w_c%reg%j_BRC = i_BRC + 2*(instance - 1) + 1
  w_c%reg%vname(i_BRC + 2*(instance - 1))     = w_c%reg%vname(i_BRC)
  w_c%reg%vname(i_BRC + 2*(instance - 1) + 1) = w_c%reg%vname(i_BRC + 1)
  
! Execute the instance method
! ---------------------------
  call Method_ ( gcBRC, w_c, impChem, expChem, ggState, &
                 nymd, nhms, cdt, rc )

! Restore the overall BRC indices
! ------------------------------
  w_c%reg%vname(i_BRC + 2*(instance - 1))     = i_qname
  w_c%reg%vname(i_BRC + 2*(instance - 1) + 1) = j_qname
  w_c%reg%n_BRC = n_BRC
  w_c%reg%i_BRC = i_BRC
  w_c%reg%j_BRC = j_BRC

  end subroutine BRC_SingleInstance_

!-----------------------------------------------------------------------
