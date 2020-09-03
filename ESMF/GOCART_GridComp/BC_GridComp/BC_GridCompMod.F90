#include "MAPL_Generic.h"

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  BC_GridCompMod --- BC Grid Component Class
!
! !INTERFACE:
!

   module  BC_GridCompMod

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
   PUBLIC  BC_GridComp       ! The BC object 
   PUBLIC  BC_GridComp1      ! Single instance BC object 

!
! !PUBLIC MEMBER FUNCTIONS:
!

   PUBLIC  BC_GridCompSetServices
   PUBLIC  BC_GridCompInitialize
   PUBLIC  BC_GridCompRun1
   PUBLIC  BC_GridCompRun2
   PUBLIC  BC_GridCompFinalize

!
! !DESCRIPTION:
!
!  This module implements the (pre-ESMF) BC Grid Component. 
!
! !REVISION HISTORY:
!
!  16Sep2003 da Silva  First crack.
!
!EOP
!-------------------------------------------------------------------------

  type BC_GridComp1
        character(len=255) :: name
        character(len=255) :: iname           ! instance name
        character(len=255) :: rcfilen         ! resource file name
        character(len=255) :: regionsString   ! Comma-delimited string of regions

        integer :: instance                   ! instance number

        logical :: run_alarm=.false.          ! run alarm

        type(Chem_Mie), pointer :: mie_tables => null()     ! aod LUTs
        real, pointer :: biofuel_src(:,:)
        real, pointer :: biomass_src_(:,:) ! before diurnal
        real, pointer :: biomass_src(:,:)
        real, pointer :: ebcant1_src(:,:)  ! level 1
        real, pointer :: ebcant2_src(:,:)  ! level 2
        real, pointer :: bc_ship_src(:,:)
        real, pointer :: aviation_lto_src(:,:)  ! aviation - landing and takeoff
        real, pointer :: aviation_cds_src(:,:)  ! aviation - climbing and descent
        real, pointer :: aviation_crs_src(:,:)  ! aviation - cruise
        real          :: aviation_layers(4)     ! heights of the LTO, CDS and CRS layers
        real :: fHydrophobic         ! Fraction of emissions hydrophobic
        integer :: myDOW = -1             ! my Day of the week: Sun=1, Mon=2,...,Sat=7
        logical :: doing_nei=.FALSE.      ! NEI08: National Emission Inventory (US+Canada)
        real    :: nei_lon(2), nei_lat(2) ! NEI bounding box; superseeds eocant1/2 inside
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
  end type BC_GridComp1

  type BC_GridComp
     integer                     :: n = 0                ! number of instances
     type(Chem_Mie), pointer     :: mie_tables => null() ! aod LUTs
     type(BC_GridComp1), pointer :: gcs(:)     => null() ! instances
  end type BC_GridComp

  character(len=*), parameter :: rc_basename = 'BC_GridComp'

  real, parameter :: OCEAN=0.0, LAND = 1.0, SEA_ICE = 2.0
  real, parameter :: radToDeg = 57.2957795

CONTAINS

   subroutine BC_GridCompSetServices(  gc, chemReg, rc)
   type(ESMF_GridComp), intent(INOUT) :: GC
   type(Chem_Registry), intent(INOUT) :: chemReg
   integer,             intent(OUT  ) :: rc

   CHARACTER(LEN=255) :: name

   integer            :: status
   character(len=ESMF_MAXSTR) :: Iam

   integer :: ier,n,i

   type(ESMF_Config) :: cfg

   Iam = "BC_GridCompSetServices"

!  Load resource file
!  ------------------
   cfg = ESMF_ConfigCreate(rc=status)
   VERIFY_(STATUS)
   call ESMF_ConfigLoadFile(cfg,trim(rc_basename)//'.rc',rc=status)
   VERIFY_(STATUS)

!  Parse resource file
!  -------------------
   n = ESMF_ConfigGetLen(cfg,label='BC_instances:',rc=status)
   VERIFY_(STATUS)


!  We have 2 tracers for each instance of BC
!  We cannot have fewer instances than half the number of
!   BC bins in the registry (it is OK to have less, though)
!  --------------------------------------------------------
   if ( n .LT. chemReg%n_BC/2 ) then
        rc = 35
        return
   else if ( n .GT. chemReg%n_BC/2 ) then
        if (MAPL_AM_I_ROOT()) &
        print *, trim(Iam)// &
                 ': fewer BC bins than possible BC instances: ',&
                 n, chemReg%n_BC/2
   end if
   n = min(n,chemReg%n_BC/2 )

!  Record name of each instance
!  ----------------------------
   call ESMF_ConfigFindLabel(cfg,'BC_instances:',rc=status)
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
      call BC_GridCompSetServices1_(gc,chemReg,name,rc=status)
      VERIFY_(STATUS)
   end do

   RETURN_(ESMF_SUCCESS)
   end subroutine BC_GridCompSetServices


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  BC_GridCompInitialize --- Initialize BC_GridComp
!
! !INTERFACE:
!

   subroutine BC_GridCompInitialize ( gcBC, w_c, impChem, expChem, &
                                      nymd, nhms, cdt, rc )

! !USES:

  implicit NONE

! !INPUT PARAMETERS:

   type(Chem_Bundle), intent(inout) :: w_c     ! Chemical tracer fields      
   integer, intent(in) :: nymd, nhms           ! time
   real, intent(in) :: cdt                     ! chemistry timestep (secs)

! !OUTPUT PARAMETERS:

   type(BC_GridComp), intent(inout) :: gcBC     ! Grid Component
   type(ESMF_State), intent(inout)  :: impChem  ! Import State
   type(ESMF_State), intent(inout)  :: expChem  ! Export State
   integer, intent(out) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 - 

! !DESCRIPTION: Initializes the BC Grid Component. It primarily sets
!               the import state for each active constituent package.
!
! !REVISION HISTORY:
!
!  18Sep2003 da Silva  First crack.
!
!EOP
!-------------------------------------------------------------------------

   character(len=*), parameter :: myname = 'BC_GridCompInitialize'
   CHARACTER(LEN=255) :: name
   
   integer :: i, ier, n

!  Load resource file
!  ------------------
   call i90_loadf ( trim(rc_basename)//'.rc', ier )
   if ( ier .NE. 0 ) then
      rc = 10
      return
   end if

!  Parse resource file
!  -------------------
   CALL I90_label ( 'BC_instances:', ier )
   if ( ier .NE. 0 ) then
      rc = 20
      return
   end if

!  First determine how many instances we have
!  ------------------------------------------   
   n = 0
   do while ( ier .EQ. 0 )
      CALL I90_gtoken( name, ier )
      if (ier .eq. 0) n = n + 1
   end do
   if ( n .EQ. 0 ) then
      rc = 30
      return
   end if
   
!  We have 2 tracers for each instance of BC
!  We cannot have fewer instances than half the number of
!   BC bins in the registry (it is OK to have less, though)
!  --------------------------------------------------------
   if ( n .LT. w_c%reg%n_BC/2 ) then
        rc = 35
        return
   else if ( n .GT. w_c%reg%n_BC/2 ) then
        if (MAPL_AM_I_ROOT()) &
        print *, trim(myname)// &
                 ': fewer BC bin sets than possible BC instances'//&
                 ' (2 bins per instance): ',&
                 n, w_c%reg%n_BC
   end if
   n = min(n,w_c%reg%n_BC/2 )
   gcBC%n = n

!  Next allocate necessary memory
!  ------------------------------
   allocate ( gcBC%gcs(n), stat=ier )    
   if ( ier .NE. 0 ) then
      rc = 40
      return
   end if

!  Record name of each instance
!  ----------------------------
   CALL I90_label ( 'BC_instances:', ier )
   do i = 1, n
      CALL I90_gtoken( name, ier )
      if ( ier .NE. 0 ) then
         rc = 40
         return
      end if
                                            ! resource file name
      gcBC%gcs(i)%rcfilen = trim(rc_basename)//'---'//trim(name)//'.rc'
      gcBC%gcs(i)%instance = i              ! instance number 
      IF(TRIM(name) == "full" ) THEN
       gcBC%gcs(i)%iname = " "              ! blank instance name for full (1)
      ELSE
       gcBC%gcs(i)%iname = TRIM(name)       ! instance name for others
      END IF
   end do    

!  Next initialize each instance
!  -----------------------------
   do i = 1, gcBC%n
      IF(MAPL_AM_I_ROOT()) THEN
       PRINT *," "
       PRINT *,myname,": Initializing instance ",TRIM(gcBC%gcs(i)%iname)," [",gcBC%gcs(i)%instance,"]"
      END IF
      call BC_SingleInstance_ ( BC_GridCompInitialize1_, i, &
                                gcBC%gcs(i), w_c, impChem, expChem,  &
                                nymd, nhms, cdt, ier )
      if ( ier .NE. 0 ) then
         rc = 1000+ier
         return
      end if
      gcBC%gcs(i)%mie_tables => gcBC%mie_tables
   end do

!  All done
!  --------
   CALL I90_FullRelease( ier )
   IF( ier /= 0 ) THEN
    PRINT *,myname,": I90_FullRelease not successful."
    rc = 40
   END IF


 end subroutine BC_GridCompInitialize

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  BC_GridCompRun1 --- Run BC_GridComp
!
! !INTERFACE:
!

   subroutine BC_GridCompRun1 ( gcBC, w_c, impChem, expChem, &
                                      nymd, nhms, cdt, rc )

! !USES:

  IMPLICIT NONE

! !INPUT PARAMETERS:

   TYPE(Chem_Bundle), intent(in) :: w_c        ! Chemical tracer fields      
   INTEGER, INTENT(IN) :: nymd, nhms           ! time
   REAL,    INTENT(IN) :: cdt                  ! chemical timestep (secs)


! !OUTPUT PARAMETERS:

   TYPE(BC_GridComp), INTENT(INOUT) :: gcBC     ! Grid Component
   TYPE(ESMF_State), INTENT(INOUT)  :: impChem  ! Import State
   TYPE(ESMF_State), INTENT(INOUT)  :: expChem  ! Export State
   INTEGER, INTENT(OUT) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 - 

! !DESCRIPTION: Runs the BC Grid Component. Multiple instance
!               version.
!
! !REVISION HISTORY:
!
!  27Feb2008  da Silva  Introduced multiple instances
!
!EOP
!-------------------------------------------------------------------------

   integer :: i, ier

   do i = 1, gcBC%n
      call BC_SingleInstance_ ( BC_GridCompRun1_, i, &
                                gcBC%gcs(i), w_c, impChem, expChem, &
                                nymd, nhms, cdt, ier )
      if ( ier .NE. 0 ) then
         rc = i * 1000+ier
         return
      end if
   end do

 end subroutine BC_GridCompRun1


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  BC_GridCompRun2 --- Run BC_GridComp
!
! !INTERFACE:
!

   subroutine BC_GridCompRun2 ( gcBC, w_c, impChem, expChem, &
                                      run_alarm, nymd, nhms, cdt, rc )

! !USES:

  IMPLICIT NONE

! !INPUT PARAMETERS:

   TYPE(Chem_Bundle), intent(in) :: w_c        ! Chemical tracer fields      
   LOGICAL, INTENT(IN) :: run_alarm            ! run alarm
   INTEGER, INTENT(IN) :: nymd, nhms           ! time
   REAL,    INTENT(IN) :: cdt                  ! chemical timestep (secs)


! !OUTPUT PARAMETERS:

   TYPE(BC_GridComp), INTENT(INOUT) :: gcBC     ! Grid Component
   TYPE(ESMF_State), INTENT(INOUT)  :: impChem  ! Import State
   TYPE(ESMF_State), INTENT(INOUT)  :: expChem  ! Export State
   INTEGER, INTENT(OUT) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 - 

! !DESCRIPTION: Runs the BC Grid Component. Multiple instance
!               version.
!
! !REVISION HISTORY:
!
!  27Feb2008  da Silva  Introduced multiple instances
!
!EOP
!-------------------------------------------------------------------------

   integer :: i, ier

   do i = 1, gcBC%n
      gcBC%gcs(i)%run_alarm = run_alarm

      call BC_SingleInstance_ ( BC_GridCompRun2_, i, &
                                gcBC%gcs(i), w_c, impChem, expChem, &
                                nymd, nhms, cdt, ier )
      if ( ier .NE. 0 ) then
         rc = i * 1000+ier
         return
      end if
   end do

 end subroutine BC_GridCompRun2



!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  BC_GridCompFinalize --- Initialize BC_GridComp
!
! !INTERFACE:
!

   subroutine BC_GridCompFinalize ( gcBC, w_c, impChem, expChem, &
                                      nymd, nhms, cdt, rc )

! !USES:

  IMPLICIT NONE

! !INPUT PARAMETERS:

   TYPE(Chem_Bundle), intent(in) :: w_c        ! Chemical tracer fields      
   INTEGER, INTENT(IN) :: nymd, nhms           ! time
   REAL,    INTENT(IN) :: cdt                  ! chemical timestep (secs)


! !OUTPUT PARAMETERS:

   TYPE(BC_GridComp), INTENT(INOUT) :: gcBC     ! Grid Component
   TYPE(ESMF_State), INTENT(INOUT)  :: impChem  ! Import State
   TYPE(ESMF_State), INTENT(INOUT)  :: expChem  ! Export State
   INTEGER, INTENT(OUT) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 - 

! !DESCRIPTION: Finalizes the BC Grid Component. Multiple instance
!               version.
!
! !REVISION HISTORY:
!
!  27Feb2008  da Silva  Introduced multiple instances
!
!EOP
!-------------------------------------------------------------------------

   integer :: i, ier

   do i = 1, gcBC%n
      call BC_SingleInstance_ ( BC_GridCompFinalize1_, i, &
                                gcBC%gcs(i), w_c, impChem, expChem, &
                                nymd, nhms, cdt, ier )
      if ( ier .NE. 0 ) then
         rc = i * 1000+ier
         return
      end if
   end do

   if (associated(gcBC%gcs)) deallocate ( gcBC%gcs, stat=ier )
   gcBC%n = -1

 end subroutine BC_GridCompFinalize


 subroutine BC_GridCompSetServices1_(  gc, chemReg, iname, rc)
   type(ESMF_GridComp), intent(INOUT) :: GC
   type(Chem_Registry), intent(INOUT) :: chemReg
   character(len=*),    intent(IN   ) :: iname
   integer,             intent(OUT  ) :: rc

   ! local
   logical:: doing_nei
 
   integer :: Status
   character(len=ESMF_MAXSTR) :: Iam

   Iam ="BC_GridCOmpSetServices1_"

   call MAPL_AddImportSpec(GC, &
       SHORT_NAME = 'BC_BIOMASS'//trim(iname), &
       LONG_NAME  = 'source species'  , &
       UNITS      = '1',                &
       DIMS       = MAPL_DimsHorzOnly,  &
       VLOCATION  = MAPL_VLocationNone, &
       RESTART    = MAPL_RestartSkip,   &
       RC         = STATUS)
   VERIFY_(STATUS)

   call MAPL_AddImportSpec(GC, &
       SHORT_NAME = 'BC_BIOFUEL'//trim(iname), &
        LONG_NAME  = 'source species'  , &
        UNITS      = '1',                &
        DIMS       = MAPL_DimsHorzOnly,  &
        VLOCATION  = MAPL_VLocationNone, &
        RESTART    = MAPL_RestartSkip,   &
        RC         = STATUS)
   VERIFY_(STATUS)

   call MAPL_AddImportSpec(GC, &
        SHORT_NAME = 'BC_ANTEBC1'//trim(iname), &
        LONG_NAME  = 'source species'  , &
        UNITS      = '1',                &
        DIMS       = MAPL_DimsHorzOnly,  &
        VLOCATION  = MAPL_VLocationNone, &
        RESTART    = MAPL_RestartSkip,   &
        RC         = STATUS)
   VERIFY_(STATUS)

   call MAPL_AddImportSpec(GC, &
        SHORT_NAME = 'BC_ANTEBC2'//trim(iname), &
        LONG_NAME  = 'source species'  , &
        UNITS      = '1',                &
        DIMS       = MAPL_DimsHorzOnly,  &
        VLOCATION  = MAPL_VLocationNone, &
        RESTART    = MAPL_RestartSkip,   &
        RC         = STATUS)
   VERIFY_(STATUS)

   call MAPL_AddImportSpec(GC, &
        SHORT_NAME = 'BC_SHIP'//trim(iname), &
        LONG_NAME  = 'source species'  , &
        UNITS      = '1',                &
        DIMS       = MAPL_DimsHorzOnly,  &
        VLOCATION  = MAPL_VLocationNone, &
        RESTART    = MAPL_RestartSkip,   &
        RC         = STATUS)
   VERIFY_(STATUS)

   call MAPL_AddImportSpec(GC, &
        SHORT_NAME = 'BC_AVIATION_LTO'//trim(iname), &
        LONG_NAME  = 'bc_aviation_lto' , &
        UNITS      = '1',                &
        DIMS       = MAPL_DimsHorzOnly,  &
        VLOCATION  = MAPL_VLocationNone, &
        RESTART    = MAPL_RestartSkip,   &
        RC         = STATUS)
   VERIFY_(STATUS)

   call MAPL_AddImportSpec(GC, &
        SHORT_NAME = 'BC_AVIATION_CDS'//trim(iname), &
        LONG_NAME  = 'bc_aviation_cds' , &
        UNITS      = '1',                &
        DIMS       = MAPL_DimsHorzOnly,  &
        VLOCATION  = MAPL_VLocationNone, &
        RESTART    = MAPL_RestartSkip,   &
        RC         = STATUS)
   VERIFY_(STATUS)

   call MAPL_AddImportSpec(GC, &
        SHORT_NAME = 'BC_AVIATION_CRS'//trim(iname), &
        LONG_NAME  = 'bc_aviation_crs' , &
        UNITS      = '1',                &
        DIMS       = MAPL_DimsHorzOnly,  &
        VLOCATION  = MAPL_VLocationNone, &
        RESTART    = MAPL_RestartSkip,   &
        RC         = STATUS)
   VERIFY_(STATUS)

!  Parse the resource file to see if NEI imports are required
!  ----------------------------------------------------------
   call doing_nei_(trim(rc_basename), trim(iname), doing_nei, __RC__)

   NEI_EMISSIONS: if (doing_nei) then
   call MAPL_AddImportSpec(GC, &
       SHORT_NAME = 'BC_NEI_BOT'//trim(iname), &
       LONG_NAME  = 'bc_nei_bot' , &
       UNITS      = '1',                &
       DIMS       = MAPL_DimsHorzOnly,  &
       VLOCATION  = MAPL_VLocationNone, &
       RESTART    = MAPL_RestartSkip,   &
       RC         = STATUS)
   VERIFY_(STATUS)

   call MAPL_AddImportSpec(GC, &
       SHORT_NAME = 'BC_NEI_TOP'//trim(iname), &
       LONG_NAME  = 'bc_nei_top' , &
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
   integer            :: status
   logical            :: isPresent
   character(len=255) :: Iam

   Iam = 'BC_GridCOmpSetServices1_::doing_nei_'
    
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

 end subroutine BC_GridCompSetServices1_


!                      Single Instance Methods

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  BC_GridCompInitialize --- Initialize BC_GridComp
!
! !INTERFACE:
!

   subroutine BC_GridCompInitialize1_ ( gcBC, w_c, impChem, expChem, &
                                      nymd, nhms, cdt, rc )

! !USES:

  implicit NONE

! !INPUT PARAMETERS:

   type(Chem_Bundle), intent(inout) :: w_c     ! Chemical tracer fields      
   integer, intent(in) :: nymd, nhms           ! time
   real, intent(in) :: cdt                     ! chemistry timestep (secs)

! !OUTPUT PARAMETERS:

   type(BC_GridComp1), intent(inout) :: gcBC    ! Grid Component
   type(ESMF_State), intent(inout)  :: impChem  ! Import State
   type(ESMF_State), intent(inout)  :: expChem  ! Export State
   integer, intent(out) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 - 

! !DESCRIPTION: Initializes the BC Grid Component. It primarily sets
!               the import state for each active constituent package.
!
! !REVISION HISTORY:
!
!  18Sep2003 da Silva  First crack.
!
!EOP
!-------------------------------------------------------------------------

   character(len=*), parameter :: myname = 'BC_GridCompInitialize1'


   character(len=255) :: rcfilen
   integer :: n
   integer :: i1, i2, im, j1, j2, jm, nbins, n1, n2, nbins_rc
   integer, allocatable :: ier(:)
   LOGICAL :: NoRegionalConstraint 



   rcfilen = gcBC%rcfilen
   gcBC%name = 'BC Constituent Package'

!  Initialize local variables
!  --------------------------
   rc = 0
   i1 = w_c%grid%i1; i2 = w_c%grid%i2; im = w_c%grid%im
   j1 = w_c%grid%j1; j2 = w_c%grid%j2; jm = w_c%grid%jm
   nbins = w_c%reg%n_BC
   n1  = w_c%reg%i_BC
   n2  = w_c%reg%j_BC

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

   call i90_label ( 'number_bc_classes:', ier(1) )
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
   gcBC%aviation_layers(1) = i90_gfloat(ier(2))
   gcBC%aviation_layers(2) = i90_gfloat(ier(3))
   gcBC%aviation_layers(3) = i90_gfloat(ier(4))
   gcBC%aviation_layers(4) = i90_gfloat(ier(5))

   if ( any(ier(1:5) /= 0) ) then
         call final_(77)
         return
   end if

!  Handle Point-wise Emission Sources Specified in a Text File
!  -----------------------------------------------------------
   ier(:) = 0
   call i90_label  ( 'point_emissions_srcfilen:',   ier(1) )
   call i90_gtoken ( gcBC%point_emissions_srcfilen, ier(2) )
   if ( ier(1) /= 0 ) then
        gcBC%doing_point_emissions = .FALSE. ! if rc is missing, don't fuss
   else if ( any(ier(2:2) /= 0) ) then
         call final_(42) ! this means point emissions info is messed up, abort
         return
   else
         if ( (index(gcBC%point_emissions_srcfilen,'/dev/null')>0) ) then
               gcBC%doing_point_emissions = .FALSE. ! disable it if no file specified
         else
               gcBC%doing_point_emissions = .TRUE.  ! we are good to go
         end if
   end if

!  Handle NEI08 Emissions
!  ----------------------
   ier(:) = 0
   call i90_label  ( 'nei_boundingbox:',   ier(1) )
   gcBC%nei_lon(1) = i90_gfloat(ier(2))
   gcBC%nei_lon(2) = i90_gfloat(ier(3))
   gcBC%nei_lat(1) = i90_gfloat(ier(4))
   gcBC%nei_lat(2) = i90_gfloat(ier(5))
   if ( ier(1) /= 0 ) then
        gcBC%doing_nei = .FALSE. ! if rc is missing, don't fuss
   else if ( any(ier(2:5) /= 0) ) then
         call final_(42) ! this means NEI info is messed up, abort
         return
   else
! --------------------------------------------------------------------------   
!         if ( (index(gcBC%nei_srcfilen(1),'/dev/null')>0) .or. &
!              (index(gcBC%nei_srcfilen(2),'/dev/null')>0) ) then 
!               gcBC%doing_nei = .FALSE. ! disable it if no file specified
!         else
!               gcBC%doing_nei = .TRUE.  ! we are good to go
!         end if
! -------------------------------------------------------------------------- 
! TODO: Need to parse the ExtData file to replicate the above logic,
!       until then do not include the NOI datasets in the ExtData primary 
!       export tables
! --------------------------------------------------------------------------

         gcBC%doing_nei = .TRUE.  ! we are good to go
   end if

   if ( MAPL_AM_I_ROOT() ) then
    if ( gcBC%doing_nei ) then
      print *, 'BC_GridComp: using NEI08 Emissions over North America'
    else
      print *, 'BC_GridComp: skipping NEI08 Emissions over North America'
    end if
   end if

!                          -------

!  Day of the week to reset tracer to zero
!  ---------------------------------------
   call i90_label ( 'my_day_of_the_week:',ier(1))
   if ( ier(1) /= 0 ) then
        gcBC%myDOW = -1   ! by default never reset tracer to zero
   else
        gcBC%myDOW = i90_gint (ier(1))
        if ( ier(1) /= 0 ) then
           call final_(60)
           return
        end if
   end if

!                          -------

!  Hydrophobic fraction
!  ---------------
   call i90_label ( 'hydrophobic_fraction:', ier(1) )
   gcBC%fHydrophobic = i90_gfloat ( ier(2) )
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
   call i90_label ( 'BC_regions_indices:', ier(1) )
   CALL I90_gtoken( gcBC%regionsString, ier(2) )
   IF( ANY(ier(1:2) < 0 ) ) THEN
    CALL final_(51)
    RETURN
   END IF

!  Is this instantiation a global case?
!  -----------------------------------
   IF(gcBC%regionsString(1:2) == "-1") THEN
    NoRegionalConstraint = .TRUE.
   ELSE
    SELECT CASE (ESMF_UtilStringLowerCase(gcBC%regionsString(1:2)))
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
   IF(NoRegionalConstraint) gcBC%regionsString = "-1"

   IF(MAPL_AM_I_ROOT()) THEN
    IF(NoRegionalConstraint) THEN
     PRINT *,myname,": This instantiation has no regional constraints."
    ELSE
     PRINT *,myname,": This instantiation is regionally constrained."
     PRINT *,myname,": List of region numbers included: ",TRIM(gcBC%regionsString)
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
   allocate ( gcBC%biomass_src(i1:i2,j1:j2), gcBC%biofuel_src(i1:i2,j1:j2), &
              gcBC%biomass_src_(i1:i2,j1:j2), &
              gcBC%ebcant1_src(i1:i2,j1:j2), gcBC%ebcant2_src(i1:i2,j1:j2), &
              gcBC%bc_ship_src(i1:i2,j1:j2), &
              gcBC%aviation_lto_src(i1:i2,j1:j2), &
              gcBC%aviation_cds_src(i1:i2,j1:j2), &
              gcBC%aviation_crs_src(i1:i2,j1:j2), ier(nerr), stat=ios )
   if ( ios /= 0 ) rc = 100
   end subroutine init_

   subroutine final_(ierr)
   integer :: ierr
   integer ios
   deallocate ( gcBC%biomass_src, gcBC%biofuel_src, &
                gcBC%biomass_src_, &
                gcBC%ebcant1_src, gcBC%ebcant2_src, &
                gcBC%bc_ship_src, &
                gcBC%aviation_lto_src, &
                gcBC%aviation_cds_src, &
                gcBC%aviation_crs_src, ier, stat=ios )
   call i90_release()
   rc = ierr
   end subroutine final_

   end subroutine BC_GridCompInitialize1_

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  BC_GridCompRun1_ --- The Chem Driver, run phase 1
!
! !INTERFACE:
!

   subroutine BC_GridCompRun1_ ( gcBC, w_c, impChem, expChem, &
                                 nymd, nhms, cdt, rc )

! !USES:

  implicit NONE

! !INPUT/OUTPUT PARAMETERS:

   type(BC_GridComp1), intent(inout) :: gcBC   ! Grid Component
   type(Chem_Bundle), intent(inout) :: w_c      ! Chemical tracer fields   

! !INPUT PARAMETERS:

   type(ESMF_State), intent(inout) :: impChem ! Import State
   integer, intent(in) :: nymd, nhms          ! time
   real, intent(in) :: cdt                    ! chemistry timestep (secs)

! !OUTPUT PARAMETERS:

   type(ESMF_State), intent(inout) :: expChem   ! Export State
   integer, intent(out) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 -
 
! !DESCRIPTION: This routine implements the so-called BC Driver. That 
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

   character(len=*), parameter :: myname = 'BC_GridCompRun1_'
   character(len=*), parameter :: Iam = myname

   integer :: i1, i2, im, j1, j2, jm, nbins, n1, n2, km, n
   integer :: i, j, ijl, ijkl, ijk1l
   real :: qmax, qmin

!  Input fields from fvGCM
!  -----------------------
   real, pointer, dimension(:,:)   :: pblh
   real, pointer, dimension(:,:,:) :: tmpu, rhoa, ple

!  Workspace for NEI emissions
!  ---------------------------
   real, pointer, dimension(:,:)         ::  nei_src1, nei_src2


   integer          :: idow
   character(len=3) :: cdow

   real, pointer :: var2D(:,:) => null()

#define EXPORT        expChem
#define iNAME         TRIM(gcBC%iname)

#define ptrBCEM       BC_emis

#define ptrBCEMAN     BC_emisAN
#define ptrBCEMBB     BC_emisBB
#define ptrBCEMBF     BC_emisBF

   
   integer :: STATUS

#include "BC_GetPointer___.h"


!  Initialize local variables
!  --------------------------
   rc = 0
   i1 = w_c%grid%i1; i2 = w_c%grid%i2; im = w_c%grid%im
   j1 = w_c%grid%j1; j2 = w_c%grid%j2; jm = w_c%grid%jm

   km = w_c%grid%km
   nbins = w_c%reg%n_BC
   n1  = w_c%reg%i_BC
   n2  = w_c%reg%j_BC

   ijl   = ( i2 - i1 + 1 ) * ( j2 - j1 + 1 )
   ijkl  = ijl * km
   ijk1l = ijl * (km+1)

! Reset tracer to zero at 0Z on specific day of week
! --------------------------------------------------
  idow = Chem_UtilIdow(nymd)
  if ( (nhms==0) .and. (idow == gcBC%myDOW) ) then
        cdow = Chem_UtilCdow(nymd)
        do n = n1, n2
           w_c%qa(n)%data3d(i1:i2,j1:j2,1:km) = tiny(1.) ! avoid division by zero
        end do
        if ( MAPL_AM_I_ROOT() ) then
           print *, '<> BC '//cdow//' tracer being set to zero on ', nymd, nhms
        end if
  end if

! Update emissions/production if necessary (daily)
!  -----------------------------------------------

!   Biomass Burning -- select on known inventories
!   ----------------------------------------------

    call MAPL_GetPointer(impChem,var2d,'BC_BIOMASS'//iNAME,rc=status)
    VERIFY_(STATUS)
    gcBC%biomass_src = var2d

!   Biofuel and anthropogenic emissions (inventories)
!   -------------------------------------------------
    call MAPL_GetPointer(impChem,var2d,'BC_BIOFUEL'//iNAME,rc=status)
    VERIFY_(STATUS)
    gcBC%biofuel_src = var2d

    call MAPL_GetPointer(impChem,var2d,'BC_ANTEBC1'//iNAME,rc=status)
    VERIFY_(STATUS)
    gcBC%ebcant1_src = var2d

    call MAPL_GetPointer(impChem,var2d,'BC_ANTEBC2'//iNAME,rc=status)
    VERIFY_(STATUS)
    gcBC%ebcant2_src = var2d

!   Ship based BC emissions
    call MAPL_GetPointer(impChem,var2d,'BC_SHIP'//iNAME,rc=status)
    VERIFY_(STATUS)
    gcBC%bc_ship_src = var2d

!   Aircraft emissions during LTO, CDS and CRS phases of flight
    call MAPL_GetPointer(impChem,var2d,'BC_AVIATION_LTO'//iNAME,rc=status)
    VERIFY_(STATUS)
    gcBC%aviation_lto_src = var2d

    call MAPL_GetPointer(impChem,var2d,'BC_AVIATION_CDS'//iNAME,rc=status)
    VERIFY_(STATUS)
    gcBC%aviation_cds_src = var2d

    call MAPL_GetPointer(impChem,var2d,'BC_AVIATION_CRS'//iNAME,rc=status)
    VERIFY_(STATUS)
    gcBC%aviation_crs_src = var2d

!   As a safety check, where value is undefined set to 0
    do j = j1, j2
     do i = i1, i2
      if(1.01*gcBC%biomass_src(i,j) .gt. undefval) gcBC%biomass_src(i,j) = 0.
      if(1.01*gcBC%biofuel_src(i,j) .gt. undefval) gcBC%biofuel_src(i,j) = 0.
      if(1.01*gcBC%ebcant1_src(i,j) .gt. undefval) gcBC%ebcant1_src(i,j) = 0.
      if(1.01*gcBC%ebcant2_src(i,j) .gt. undefval) gcBC%ebcant2_src(i,j) = 0.
      if(1.01*gcBC%bc_ship_src(i,j) .gt. undefval) gcBC%bc_ship_src(i,j) = 0.
      if(1.01*gcBC%aviation_lto_src(i,j) .gt. undefval) gcBC%aviation_lto_src(i,j) = 0.
      if(1.01*gcBC%aviation_cds_src(i,j) .gt. undefval) gcBC%aviation_cds_src(i,j) = 0.
      if(1.01*gcBC%aviation_crs_src(i,j) .gt. undefval) gcBC%aviation_crs_src(i,j) = 0.
     enddo
    enddo


#ifdef DEBUG
    call pmaxmin ( 'BC: biomass', gcBC%biomass_src, qmin, qmax, ijl,1, 1. )
    call pmaxmin ( 'BC: biofuel', gcBC%biofuel_src, qmin, qmax, ijl,1, 1. )
    call pmaxmin ( 'BC: ebcant1', gcBC%ebcant1_src, qmin, qmax, ijl,1,1.)
    call pmaxmin ( 'BC: ebcant2', gcBC%ebcant2_src, qmin, qmax, ijl,1,1.)
    call pmaxmin ( 'BC: bc_ship', gcBC%bc_ship_src, qmin, qmax, ijl,1,1.)
    call pmaxmin ( 'BC: avi_lto', gcBC%aviation_lto_src, qmin, qmax, ijl,1,1.)
    call pmaxmin ( 'BC: avi_cds', gcBC%aviation_cds_src, qmin, qmax, ijl,1,1.)
    call pmaxmin ( 'BC: avi_crs', gcBC%aviation_crs_src, qmin, qmax, ijl,1,1.)
#endif

!   Save this in case we need to apply diurnal cycle
!   ------------------------------------------------
   if ( w_c%diurnal_bb ) then
        gcBC%biomass_src_(:,:) = gcBC%biomass_src(:,:)
   end if

!  Apply diurnal cycle if so desired
!  ---------------------------------
   if ( w_c%diurnal_bb ) then
      call Chem_BiomassDiurnal ( gcBC%biomass_src, gcBC%biomass_src_,   &
                                 w_c%grid%lon(:,:)*radToDeg, &
                                 w_c%grid%lat(:,:)*radToDeg, nhms, cdt )      
   end if

!  Read any pointwise emissions, if requested
!  ------------------------------------------
   if(gcBC%doing_point_emissions) then
    call Chem_UtilPointEmissions( nymd, gcBC%point_emissions_srcfilen, &
                                  gcBC%nPts, gcBC%vLat, gcBC%vLon, &
                                  gcBC%vBase, gcBC%vTop, gcBC%vEmis, &
                                  gcBC%vStart, gcBC%vEnd )

!   In case vStart or vEnd were not specified in the file set to defaults
    where(gcBC%vStart < 0) gcBC%vStart = 000000
    where(gcBC%vEnd < 0)   gcBC%vEnd   = 240000
   endif


!  Apply NEI emissions over North America if so desired
!  ----------------------------------------------------
   if (gcBC%doing_NEI) then

       allocate(nei_src1(i1:i2,j1:j2), nei_src2(i1:i2,j1:j2), __STAT__)

       call MAPL_GetPointer(impChem,var2d,'BC_NEI_BOT'//iNAME, __RC__)
       nei_src1 = var2d

       call MAPL_GetPointer(impChem,var2d,'BC_NEI_TOP'//iNAME, __RC__)
       nei_src2 = var2d

       where ( (w_c%grid%lon.ge.gcBC%nei_lon(1)) .and. &
               (w_c%grid%lon.le.gcBC%nei_lon(2)) .and. &
               (w_c%grid%lat.ge.gcBC%nei_lat(1)) .and. &
               (w_c%grid%lat.le.gcBC%nei_lat(2))   )

               gcBC%ebcant1_src = nei_src1
               gcBC%ebcant2_src = nei_src2
       end where

#ifdef DEBUG
            call pmaxmin('BC: nei_bot', nei_src1, qmin, qmax, ijl,1, 1. )
            call pmaxmin('BC: nei_top', nei_src2, qmin, qmax, ijl,1, 1. )
#endif

            deallocate(nei_src1, nei_src2)

   end if ! doing NEI

#ifdef DEBUG
   do n = n1, n2
      call pmaxmin ( 'BC: q_beg', w_c%qa(n)%data3d(i1:i2,j1:j2,1:km), qmin, qmax, &
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



#ifdef DEBUG

   call pmaxmin('BC: pblh       ', pblh    , qmin, qmax, ijl,1, 1. )

   call pmaxmin('BC: tmpu       ', tmpu    , qmin, qmax, ijkl,1, 1. )
   call pmaxmin('BC: rhoa       ', rhoa    , qmin, qmax, ijkl,1, 1. )

#endif

!  BC Source
!  -----------
   call BC_Emission ( i1, i2, j1, j2, km, nbins, cdt, gcBC, w_c, &
                      pblh, tmpu, rhoa, BC_emis, &
                      BC_emisAN, BC_emisBB, BC_emisBF, rc )

#ifdef DEBUG
   do n = n1, n2
      call pmaxmin('BC: q_emi', w_c%qa(n)%data3d(i1:i2,j1:j2,1:km), &
                    qmin, qmax, ijl, km, 1. )
   end do
#endif

   return

CONTAINS

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  BC_Emission - Adds Black Carbon emission for one timestep
!             We have emissions from 5 sources, which are distributed
!             differently in the vertical
!             1) biomass burning - uniformly mixed in PBL
!             2) biofuel sources - emitted into lowest 100 m
!             3) anthropogenic l1 - emitted into lowest 100 m
!             4) anthropogenic l2 - emitted into 100 - 500 m levels
!             5) point sources    - emitted in altitudes specified in input
!
! !INTERFACE:
!

   subroutine BC_Emission ( i1, i2, j1, j2, km, nbins, cdt, gcBC, w_c, &
                            pblh, tmpu, rhoa, BC_emis, &
                            BC_emisAN, BC_emisBB, BC_emisBF, rc )

! !USES:

  implicit NONE

! !INPUT PARAMETERS:

   integer, intent(in) :: i1, i2, j1, j2, km, nbins
   real, intent(in)    :: cdt
   type(BC_GridComp1), intent(in)    :: gcBC       ! BC Grid Component
   real, pointer, dimension(:,:)    :: pblh
   real, pointer, dimension(:,:,:)  :: tmpu
   real, pointer, dimension(:,:,:)  :: rhoa

! !OUTPUT PARAMETERS:

   type(Chem_Bundle), intent(inout) :: w_c         ! Chemical tracer fields
   type(Chem_Array), intent(inout)  :: BC_emis(nbins) ! BC emissions, kg/m2/s
   type(Chem_Array), intent(inout)  :: BC_emisAN      ! BC emissions, kg/m2/s
   type(Chem_Array), intent(inout)  :: BC_emisBB      ! BC emissions, kg/m2/s
   type(Chem_Array), intent(inout)  :: BC_emisBF      ! BC emissions, kg/m2/s
   integer, intent(out)             :: rc          ! Error return code:
                                                   !  0 - all is well
                                                   !  1 - 
   character(len=*), parameter :: myname = 'BC_Emission'

! !DESCRIPTION: Updates the BC concentration with emissions every timestep
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
                                       ! pressure at 100m, 500m, & PBLH
   real, dimension(i1:i2,j1:j2) :: p100, p500, pPblh  
   real, dimension(i1:i2,j1:j2) :: p0, z0, ps
   real :: p1, z1, dz, delz, delp, f100, f500, fPblh
   real :: qmax, qmin 

   real, dimension(i1:i2,j1:j2) :: factor, srcHydrophobic, srcHydrophilic
   real, dimension(i1:i2,j1:j2) :: srcBiofuel, srcBiomass, srcAnthro
   real                         :: srcAll, zpbl, maxAll

   real, dimension(i1:i2,j1:j2,km) :: emis_aviation
   real, dimension(i1:i2,j1:j2,km) :: srcAviation
   real                            :: z_lto_bot, z_lto_top
   real                            :: z_cds_bot, z_cds_top
   real                            :: z_crs_bot, z_crs_top

!  Indices for point emissions
   integer, pointer, dimension(:)  :: iPoint, jPoint
   real, dimension(km)             :: point_column_emissions

!  Initialize local variables
!  --------------------------
   n1  = w_c%reg%i_BC
   n2  = w_c%reg%j_BC
   ijl = ( i2 - i1 + 1 ) * ( j2 - j1 + 1 )

!  Zero diagnostic accumulators
   do n = 1, nbins
     if( associated(BC_emis(n)%data2d) ) BC_emis(n)%data2d = 0.0
   end do
     if(associated(BC_emisAN%data2d) )   BC_emisAN%data2d  = 0.0
     if(associated(BC_emisBF%data2d) )   BC_emisBF%data2d  = 0.0
     if(associated(BC_emisBB%data2d) )   BC_emisBB%data2d  = 0.0

!  Distribute aircraft emissions from LTO, CDS and CRS layers
!  ----------------------------------------------------------
   z_lto_bot = max(1e-3, gcBC%aviation_layers(1))
   z_lto_top = max(2e-3, gcBC%aviation_layers(2))

   z_cds_bot = max(2e-3, gcBC%aviation_layers(2))
   z_cds_top = max(3e-3, gcBC%aviation_layers(3))

   z_crs_bot = max(3e-3, gcBC%aviation_layers(3))
   z_crs_top = max(4e-3, gcBC%aviation_layers(4))

   emis_aviation = 0.0
   srcAviation   = 0.0

   call distribute_aviation_emissions(w_c%delp, rhoa, z_lto_bot, z_lto_top, gcBC%aviation_lto_src, emis_aviation, i1, i2, j1, j2, km)
   srcAviation = srcAviation + emis_aviation

   call distribute_aviation_emissions(w_c%delp, rhoa, z_cds_bot, z_cds_top, gcBC%aviation_cds_src, emis_aviation, i1, i2, j1, j2, km)
   srcAviation = srcAviation + emis_aviation

   call distribute_aviation_emissions(w_c%delp, rhoa, z_crs_bot, z_crs_top, gcBC%aviation_crs_src, emis_aviation, i1, i2, j1, j2, km)
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
       pPblh(i,j) = p1+delp
      endif
      p0(i,j) = p1
      z0(i,j) = z1
     end do
    end do
   end do

#if 0
   call pmaxmin ( 'BC: p100   ', p100,  qmin, qmax, ijl, 1, 1. )
   call pmaxmin ( 'BC: p500   ', p500,  qmin, qmax, ijl, 1, 1. )
   call pmaxmin ( 'BC: pPBL   ', pPBLh, qmin, qmax, ijl, 1, 1. )
#endif

!  Now update the tracer mixing ratios with the aerosol sources
!  ------------------------------------------------------------
   p0 = ps
K_LOOP: do k = km, 1, -1

!!!    print *, 'BC_Emissions: getting emissions for layer ', k

!   First determine emissions for this layer
!   ----------------------------------------
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
      fPblh = 0.
      if(p1 .ge. pPblh(i,j)) fPblh = w_c%delp(i,j,k)/(ps(i,j)-pPblh(i,j))
      if(p1 .lt. pPblh(i,j) .and. p0(i,j) .ge. pPblh(i,j)) &
       fPblh = (p0(i,j)-pPblh(i,j))/(ps(i,j)-pPblh(i,j))

!     Sources by class in kg m-2 s-1
!     ------------------------------
      srcBiofuel(i,j) = f100 * gcBC%biofuel_src(i,j)
      srcAnthro(i,j)  = f100 * gcBC%ebcant1_src(i,j) &
                      + f500 * gcBC%ebcant2_src(i,j) &
                      + f100 * gcBC%bc_ship_src(i,j) &
                      +        srcAviation(i,j,k)

      srcBiomass(i,j) = fPblh*gcBC%biomass_src(i,j)

      srcAll = srcBiofuel(i,j) + srcAnthro(i,j) + srcBiomass(i,j)
      srcHydrophobic(i,j) =     gcBC%fHydrophobic  * srcAll
      srcHydrophilic(i,j) = (1.-gcBC%fHydrophobic) * srcAll

!     Update pressure of lower level
!     ------------------------------
      p0(i,j) = p1

     end do ! i
    end do  ! j

!   Determine global max/min
!   ------------------------
    call pmaxmin ( 'BC: Phobic ', srcHydrophobic, qmin, qmax, ijl, 1, 0. )
    maxAll = abs(qmax) + abs(qmin)
    call pmaxmin ( 'BC: Philic ', srcHydrophilic, qmin, qmax, ijl, 1, 0. )
    maxAll = max ( maxAll, abs(qmax) + abs(qmin) )

!   If emissions are zero at this level (globally), we are done
!   -----------------------------------------------------------
    if ( maxAll .eq. 0.0 ) exit K_LOOP

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
    if ( associated(BC_emis(1)%data2d)) &
                    BC_emis(1)%data2d = BC_emis(1)%data2d + srcHydrophobic

    if ( associated(BC_emis(2)%data2d)) &
                    BC_emis(2)%data2d = BC_emis(2)%data2d + srcHydrophilic

    if ( associated(BC_emisBF%data2d)) &
                    BC_emisBF%data2d  = BC_emisBF%data2d  + srcBiofuel

    if ( associated(BC_emisBB%data2d)) &
                    BC_emisBB%data2d  = BC_emisBB%data2d  + srcBiomass

    if ( associated(BC_emisAN%data2d)) &
                    BC_emisAN%data2d  = BC_emisAN%data2d  + srcAnthro

   end do K_LOOP

!  Distribute pointwise sources if requested
!  -----------------------------------------
   if( gcBC%doing_point_emissions .and. gcBC%nPts > 0) then

!    Get indices for point emissions
!    -------------------------------
     allocate(iPoint(gcBC%nPts), jPoint(gcBC%nPts), stat=ios)

     call MAPL_GetHorzIJIndex(gcBC%nPts, iPoint, jPoint, &
                              grid = w_c%grid_esmf,      &
                              lon  = gcBC%vLon/radToDeg, &
                              lat  = gcBC%vLat/radToDeg, &
                              rc   = rc)

     if ( rc /= 0 ) call die(myname,'cannot get indices for point emissions')

     do ii = 1, gcBC%nPts
      i = iPoint(ii)
      j = jPoint(ii)
      if( i<1 .OR. j<1 )              cycle    ! point emission not in this sub-domain
!      if( gcBC%regionMask(i,j) == 0 ) cycle    ! masked by region mask
      
!     Emissions not occurring in current time step
!     --------------------------------------------
      if(nhms < gcBC%vStart(ii) .or. nhms >= gcBC%vEnd(ii)) cycle

      call distribute_point_emissions(w_c%delp(i,j,:), rhoa(i,j,:), &
                                      gcBC%vBase(ii), gcBC%vTop(ii), gcBC%vEmis(ii), &
                                      point_column_emissions, km)
      w_c%qa(n1)%data3d(i,j,:) = w_c%qa(n1)%data3d(i,j,:) & 
         + gcBC%fHydrophobic * cdt * grav / w_c%delp(i,j,:) &
                             * point_column_emissions / w_c%grid%cell_area(i,j)
      w_c%qa(n2)%data3d(i,j,:) = w_c%qa(n2)%data3d(i,j,:) & 
         + (1-gcBC%fHydrophobic) * cdt * grav / w_c%delp(i,j,:) &
                                 * point_column_emissions / w_c%grid%cell_area(i,j)

     enddo
     deallocate(iPoint, jPoint, stat=ios)
   endif



   rc = 0

   end subroutine BC_Emission

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


!  Abstracted from distribute_aviation_emissions above, but called per column
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

 end subroutine BC_GridCompRun1_


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  BC_GridCompRun2_ --- The Chem Driver, run phase 2
!
! !INTERFACE:
!

   subroutine BC_GridCompRun2_ ( gcBC, w_c, impChem, expChem, &
                                 nymd, nhms, cdt, rc )

! !USES:

  implicit NONE

! !INPUT/OUTPUT PARAMETERS:

   type(BC_GridComp1), intent(inout) :: gcBC   ! Grid Component
   type(Chem_Bundle), intent(inout) :: w_c      ! Chemical tracer fields   

! !INPUT PARAMETERS:

   type(ESMF_State), intent(inout) :: impChem ! Import State
   integer, intent(in) :: nymd, nhms          ! time
   real, intent(in) :: cdt                    ! chemistry timestep (secs)

! !OUTPUT PARAMETERS:

   type(ESMF_State), intent(inout) :: expChem   ! Export State
   integer, intent(out) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 -
 
! !DESCRIPTION: This routine implements the so-called BC Driver. That 
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

   character(len=*), parameter :: myname = 'BC_GridCompRun2_'
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

   real, pointer    :: BC_radius(:), BC_rhop(:)
   integer          :: rhFlag


#define EXPORT        expChem
#define iNAME         TRIM(gcBC%iname)

#define ptrBCWT       BC_wet
#define ptrBCSV       BC_conv
#define ptrBCEM       BC_emis
#define ptrBCDP       BC_dep
#define ptrBCSD       BC_set

#define ptrBCMASS     BC_mass
#define ptrBCEMAN     BC_emisAN
#define ptrBCEMBB     BC_emisBB
#define ptrBCEMBF     BC_emisBF
#define ptrBCHYPHIL   BC_toHydrophilic
#define ptrBCSMASS    BC_sfcmass
#define ptrBCCMASS    BC_colmass
#define ptrBCEXTTAU   BC_exttau
#define ptrBCSCATAU   BC_scatau
#define ptrBCCONC     BC_conc
#define ptrBCEXTCOEF  BC_extcoef
#define ptrBCSCACOEF  BC_scacoef
#define ptrBCANGSTR   BC_angstrom
#define ptrBCFLUXU    BC_fluxu
#define ptrBCFLUXV    BC_fluxv


   
   integer :: STATUS

#include "BC_GetPointer___.h"


!  Initialize local variables
!  --------------------------
   rc = 0
   i1 = w_c%grid%i1; i2 = w_c%grid%i2; im = w_c%grid%im
   j1 = w_c%grid%j1; j2 = w_c%grid%j2; jm = w_c%grid%jm

   km = w_c%grid%km
   nbins = w_c%reg%n_BC
   n1  = w_c%reg%i_BC
   n2  = w_c%reg%j_BC

   ijl   = ( i2 - i1 + 1 ) * ( j2 - j1 + 1 )
   ijkl  = ijl * km
   ijk1l = ijl * (km+1)


#ifdef DEBUG
   do n = n1, n2
      call pmaxmin ( 'BC: q_beg', w_c%qa(n)%data3d(i1:i2,j1:j2,1:km), qmin, qmax, &
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

   call pmaxmin('BC: frlake     ', frlake  , qmin, qmax, ijl,1, 1. )
   call pmaxmin('BC: frocean    ', frocean , qmin, qmax, ijl,1, 1. )
   call pmaxmin('BC: frseaice   ', frseaice, qmin, qmax, ijl,1, 1. )
   call pmaxmin('BC: area       ', area    , qmin, qmax, ijl,1, 1. )
   call pmaxmin('BC: oro        ', oro     , qmin, qmax, ijl,1, 1. )
   call pmaxmin('BC: u10m       ', u10m    , qmin, qmax, ijl,1, 1. )
   call pmaxmin('BC: v10m       ', v10m    , qmin, qmax, ijl,1, 1. )
   call pmaxmin('BC: ustar      ', ustar   , qmin, qmax, ijl,1, 1. )
   call pmaxmin('BC: precc      ', precc   , qmin, qmax, ijl,1, 1. )
   call pmaxmin('BC: precl      ', precl   , qmin, qmax, ijl,1, 1. )
   call pmaxmin('BC: pblh       ', pblh    , qmin, qmax, ijl,1, 1. )
   call pmaxmin('BC: shflux     ', shflux  , qmin, qmax, ijl,1, 1. )
   call pmaxmin('BC: z0h        ', z0h     , qmin, qmax, ijl,1, 1. )
   call pmaxmin('BC: hsurf      ', hsurf   , qmin, qmax, ijl,1, 1. )

   call pmaxmin('BC: tmpu       ', tmpu    , qmin, qmax, ijkl,1, 1. )
   call pmaxmin('BC: rhoa       ', rhoa    , qmin, qmax, ijkl,1, 1. )
   call pmaxmin('BC: u          ', u       , qmin, qmax, ijkl,1, 1. )
   call pmaxmin('BC: v          ', v       , qmin, qmax, ijkl,1, 1. )
   call pmaxmin('BC: hghte      ', hghte   , qmin, qmax, ijkl,1, 1. )
   call pmaxmin('BC: qlcn       ', qlcn    , qmin, qmax, ijkl,1, 1. )
   call pmaxmin('BC: qicn       ', qicn    , qmin, qmax, ijkl,1, 1. )
   call pmaxmin('BC: cmfmc      ', cmfmc   , qmin, qmax, ijkl,1, 1. )
   call pmaxmin('BC: dtrain     ', dtrain  , qmin, qmax, ijkl,1, 1. )
   call pmaxmin('BC: pfllsan    ', pfllsan , qmin, qmax, ijk1l,1, 1. )
   call pmaxmin('BC: pfilsan    ', pfilsan , qmin, qmax, ijk1l,1, 1. )

#endif

RUN_ALARM: if (gcBC%run_alarm) then 

   allocate( fluxout )
   allocate( fluxout%data2d(i1:i2,j1:j2), dqa(i1:i2,j1:j2), &
             drydepositionfrequency(i1:i2,j1:j2), stat=STATUS)
   VERIFY_(STATUS)



!  Ad Hoc transfer of hydrophobic to hydrophilic aerosols
!  Following Chin's parameterization, the rate constant is
!  k = 4.63e-6 s-1 (.4 day-1; e-folding time = 2.5 days)
   if(associated(BC_toHydrophilic%data2d)) &
     BC_toHydrophilic%data2d(i1:i2,j1:j2) = 0.0

   do k = 1, km
    do j = j1, j2
     do i = i1, i2
      qUpdate = w_c%qa(n1)%data3d(i,j,k)*exp(-4.63e-6*cdt)
      qUpdate = max(qUpdate,1e-32)
      delq = max(0.,w_c%qa(n1)%data3d(i,j,k)-qUpdate)
      w_c%qa(n1)%data3d(i,j,k) = qUpdate
      w_c%qa(n2)%data3d(i,j,k) = w_c%qa(n2)%data3d(i,j,k)+delq
      if(associated(BC_toHydrophilic%data2d)) &
       BC_toHydrophilic%data2d(i,j) = BC_toHydrophilic%data2d(i,j) &
        + delq*w_c%delp(i,j,k)/grav/cdt
     end do
    end do
   end do

!  BC Settling
!  -----------
   allocate( BC_radius(nbins), BC_rhop(nbins) )
   BC_radius(:) = 0.35e-6  ! radius for settling [m]
   BC_rhop(:)   = 1800.    ! density for setting [kg m-3]
   rhFlag       = 0        ! settle like dry particles
   call Chem_Settling ( i1, i2, j1, j2, km, n1, n2, nbins, rhFlag, &
                        BC_radius, BC_rhop, cdt, w_c, tmpu, rhoa, hsurf,    &
                        hghte, BC_set, rc )
   deallocate( BC_radius, BC_rhop)

!  BC Deposition
!  -----------
   drydepositionfrequency = 0.
   call DryDepositionGOCART( i1, i2, j1, j2, km, &
                             tmpu, rhoa, hghte, oro, ustar, &
                             pblh, shflux, z0h, drydepositionfrequency, rc )
    
   do n = 1, nbins
    dqa = 0.
    dqa = max(0.0, w_c%qa(n1+n-1)%data3d(:,:,km)*(1.-exp(-drydepositionfrequency*cdt)))
    w_c%qa(n1+n-1)%data3d(:,:,km) = &
            w_c%qa(n1+n-1)%data3d(:,:,km) - dqa
    if( associated(BC_dep(n)%data2d) ) &
     BC_dep(n)%data2d = dqa*w_c%delp(:,:,km)/grav/cdt
   end do

#ifdef DEBUG
   do n = n1, n2
      call pmaxmin('BC: q_dry', w_c%qa(n)%data3d(i1:i2,j1:j2,1:km), qmin, qmax, &
                    ijl, km, 1. )
   end do
#endif


!  BC Large-scale Wet Removal
!  --------------------------
!  Hydrophobic mode (first tracer) is not removed
   if(associated(BC_wet(1)%data2d)) BC_wet(1)%data2d = 0.
!  Hydrophilic mode (second tracer) is removed
   KIN = .TRUE.
   do n = nbins, nbins
    w_c%qa(n1+n-1)%fwet = 1.
    call WetRemovalGOCART(i1, i2, j1, j2, km, n1+n-1, n1+n-1, cdt, 'BC', KIN, &
                          w_c%qa, ple, tmpu, rhoa, pfllsan, pfilsan, &
                          precc, precl, fluxout, rc )
    if(associated(BC_wet(n)%data2d)) BC_wet(n)%data2d = fluxout%data2d
   end do

#ifdef DEBUG
   do n = n1, n2
      call pmaxmin('BC: q_wet', w_c%qa(n)%data3d(i1:i2,j1:j2,1:km), qmin, qmax, &
                    ijl, km, 1. )
   end do
#endif

!  Black Carbon Convective-scale Mixing and Wet Removal
!  ----------------------------------------------------
   KIN = .TRUE.
   icdt = cdt
   allocate(cmfmc_(i1:i2,j1:j2,km+1), qccu_(i1:i2,j1:j2,km), &
            dtrain_(i1:i2,j1:j2,km), airmass_(i1:i2,j1:j2,km), &
            delz_(i1:i2,j1:j2,km), vud_(i1:i2,j1:j2,km), &
            tc_(i1:i2,j1:j2,km,n1:n2), delp_(i1:i2,j1:j2,km), ple_(i1:i2,j1:j2,km+1), &
            airmol_(i1:i2,j1:j2,km), tmpu_(i1:i2,j1:j2,km), bcnv_(i1:i2,j1:j2,n1:n2), &
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
   call convection(i1, i2, j1, j2, km, n1, n2, icdt, 'BC', kin, &
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
   if(associated(BC_conv(1)%data2d)) BC_conv(1)%data2d = 0.0
   if(associated(BC_conv(2)%data2d)) BC_conv(2)%data2d = -bcnv_(:,:,n2)/area_/icdt


!  Clean up
!  --------
   deallocate(cmfmc_, qccu_, dtrain_, tc_, airmass_, &
              delz_, vud_, delp_, airmol_, tmpu_, bcnv_, ple_, &
              area_, frlake_, frocean_, frseaice_, __STAT__ )

   deallocate(fluxout%data2d)
   deallocate(fluxout, dqa, drydepositionfrequency, stat=ios )

   end if RUN_ALARM

!  Compute the desired output diagnostics here
!  Ideally this will go where chemout is called in fvgcm.F since that
!  will reflect the distributions after transport, etc.
!  ------------------------------------------------------------------
   call BC_Compute_Diags(i1, i2, j1, j2, km, nbins, gcBC, w_c, tmpu, rhoa, u, v, &
                         BC_sfcmass, BC_colmass, BC_mass, BC_exttau, &
                         BC_scatau, BC_conc, BC_extcoef, BC_scacoef, BC_angstrom, &
                         BC_fluxu, BC_fluxv, rc)

   return

CONTAINS

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  BC_Compute_Diags - Calculate dust 2D diagnostics
!
! !INTERFACE:
!

   subroutine BC_Compute_Diags ( i1, i2, j1, j2, km, nbins, gcBC, w_c, tmpu, rhoa, u, v, &
                                 sfcmass, colmass, mass, exttau, scatau, &
                                 conc, extcoef, scacoef, angstrom, fluxu, fluxv, rc )

! !USES:

  implicit NONE

! !INPUT PARAMETERS:
   integer, intent(in) :: i1, i2, j1, j2, km, nbins
   type(BC_GridComp1), intent(inout):: gcBC     ! BC Grid Component
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

! !DESCRIPTION: Calculates some simple 2d diagnostics from the BC fields
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
   character(len=*), parameter :: myname = 'BC_Compute_Diags'
   integer :: i, j, k, n, n1, n2, ios, nch, idx
   real :: tau, ssa
   character(len=255) :: qname
   real, dimension(i1:i2,j1:j2) :: tau470, tau870
   real    :: ilam550, ilam470, ilam870
   logical :: do_angstrom


!  Initialize local variables
!  --------------------------
   n1  = w_c%reg%i_BC
   n2  = w_c%reg%j_BC
   nch   = gcBC%mie_tables%nch

!  Get the wavelength indices
!  --------------------------
!  Must provide ilam550 for AOT calculation
   ilam550 = 1.
   ilam470 = 0.
   ilam870 = 0.
   if(nch .gt. 1) then
    do i = 1, nch
     if ( gcBC%mie_tables%channels(i) .ge. 5.49e-7 .and. &
          gcBC%mie_tables%channels(i) .le. 5.51e-7) ilam550 = i
     if ( gcBC%mie_tables%channels(i) .ge. 4.69e-7 .and. &
          gcBC%mie_tables%channels(i) .le. 4.71e-7) ilam470 = i
     if ( gcBC%mie_tables%channels(i) .ge. 8.69e-7 .and. &
          gcBC%mie_tables%channels(i) .le. 8.71e-7) ilam870 = i
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

!      Select the name for species
       qname = trim(w_c%reg%vname(n+n1-1))
       idx = Chem_MieQueryIdx(gcBC%mie_tables,qname,rc)
       if(rc .ne. 0) call die(myname, 'cannot find proper Mie table index')

       do k = 1, km
        do j = j1, j2
         do i = i1, i2
          call Chem_MieQuery(gcBC%mie_tables, idx, ilam550, &
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
       idx = Chem_MieQueryIdx(gcBC%mie_tables,qname,rc)
       if(rc .ne. 0) call die(myname, 'cannot find proper Mie table index')

       do k = 1, km
        do j = j1, j2
         do i = i1, i2

          call Chem_MieQuery(gcBC%mie_tables, idx, ilam470, &
              w_c%qa(n+n1-1)%data3d(i,j,k)*w_c%delp(i,j,k)/grav, &
              w_c%rh(i,j,k), tau=tau)
          tau470(i,j) = tau470(i,j) + tau

          call Chem_MieQuery(gcBC%mie_tables, idx, ilam870, &
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

   end subroutine BC_Compute_Diags

 end subroutine BC_GridCompRun2_



!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  BC_GridCompFinalize --- The Chem Driver 
!
! !INTERFACE:
!

   subroutine BC_GridCompFinalize1_ ( gcBC, w_c, impChem, expChem, &
                                      nymd, nhms, cdt, rc )

! !USES:

  implicit NONE

! !INPUT/OUTPUT PARAMETERS:

   type(BC_GridComp1), intent(inout) :: gcBC   ! Grid Component

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
   character(len=*), parameter :: myname = 'BC_GridCompFinalize'

!  If initialized pointwise emissions from daily tables, clean-up
   if(associated(gcBC%vLat))    deallocate(gcBC%vLat, stat=ios)
   if(associated(gcBC%vLon))    deallocate(gcBC%vLon, stat=ios)
   if(associated(gcBC%vEmis))   deallocate(gcBC%vEmis, stat=ios)
   if(associated(gcBC%vBase))   deallocate(gcBC%vBase, stat=ios)
   if(associated(gcBC%vTop))    deallocate(gcBC%vTop, stat=ios)
   if(associated(gcBC%vStart))  deallocate(gcBC%vStart, stat=ios)
   if(associated(gcBC%vEnd))    deallocate(gcBC%vEnd, stat=ios)

   rc=0
   return

 end subroutine BC_GridCompFinalize1_

 end module BC_GridCompMod


!-----------------------------------------------------------------------

!                     Single Instance Wrapper

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  BC_SingleInstance_ --- Runs single instance of method
!
! !INTERFACE:
!
  subroutine BC_SingleInstance_ ( Method_, instance, &
                                  gcBC, w_c, impChem, expChem, &
                                  nymd, nhms, cdt, rc )

! !USES:

  Use BC_GridCompMod
  Use ESMF
  Use MAPL
  Use Chem_Mod 

  IMPLICIT NONE

! !INPUT PARAMETERS:

!  Input "function pointer"
!  -----------------------
   interface 
     subroutine Method_ (gc, w, imp, exp, ymd, hms, dt, rcode )
       Use BC_GridCompMod
       Use ESMF
       Use MAPL
       Use Chem_Mod 
       type(BC_GridComp1),  intent(inout)  :: gc
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

   TYPE(BC_GridComp1), INTENT(INOUT) :: gcBC    ! Grid Component
   TYPE(ESMF_State), INTENT(INOUT)  :: impChem  ! Import State
   TYPE(ESMF_State), INTENT(INOUT)  :: expChem  ! Export State
   INTEGER, INTENT(OUT) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 - 

! !DESCRIPTION: Finalizes the BC Grid Component. Multiple instance
!               version.
!
! !REVISION HISTORY:
!
!  27Feb2008  da Silva  Introduced multiple instances
!
!EOP
!-------------------------------------------------------------------------

  integer n_BC, i_BC, j_BC
  character(len=255) :: i_qname, j_qname

! Save overall BC indices
! -----------------------
  n_BC = w_c%reg%n_BC
  i_BC = w_c%reg%i_BC
  j_BC = w_c%reg%j_BC

! Save the name of the variables in this instance
! -----------------------------------------------
  i_qname = trim(w_c%reg%vname(i_BC + 2*(instance - 1)))
  j_qname = trim(w_c%reg%vname(i_BC + 2*(instance - 1) + 1))
  
! Customize indices for this particular instance
! ----------------------------------------------
  w_c%reg%n_BC = 2
  w_c%reg%i_BC = i_BC + 2*(instance - 1)
  w_c%reg%j_BC = i_BC + 2*(instance - 1) + 1
  w_c%reg%vname(i_BC + 2*(instance - 1))     = w_c%reg%vname(i_BC)
  w_c%reg%vname(i_BC + 2*(instance - 1) + 1) = w_c%reg%vname(i_BC + 1)
  
! Execute the instance method
! ---------------------------
  call Method_ ( gcBC, w_c, impChem, expChem, &
                 nymd, nhms, cdt, rc )

! Restore the overall BC indices
! ------------------------------
  w_c%reg%vname(i_BC + 2*(instance - 1))     = i_qname
  w_c%reg%vname(i_BC + 2*(instance - 1) + 1) = j_qname
  w_c%reg%n_BC = n_BC
  w_c%reg%i_BC = i_BC
  w_c%reg%j_BC = j_BC

  end subroutine BC_SingleInstance_

!-----------------------------------------------------------------------
