#include "MAPL_Generic.h"

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  NI_GridCompMod --- NI Grid Component Class
!
! !INTERFACE:
!

   module  NI_GridCompMod

! !USES:

   USE ESMF
   USE MAPL
   USE MAPL_ConstantsMod, only: MAPL_AIRMW, MAPL_AVOGAD, MAPL_PI

   use Chem_Mod              ! Chemistry Base Class
   use Chem_StateMod         ! Chemistry State
   use Chem_ConstMod, only: grav, von_karman, cpd, &
                            undefval => undef         ! Constants !
   use Chem_UtilMod          ! I/O
   use Chem_MieMod           ! Aerosol LU Tables, calculator
   use m_inpak90             ! Resource file management
   use m_die, only: die
   use NitrateChemDriverMod, only: RPMARES, sktrs_hno3, sktrs_sslt
   use Chem_SettlingMod      ! Settling
   use DryDepositionMod      ! Dry Deposition
   use WetRemovalMod         ! Large-scale Wet Removal
   use ConvectionMod         ! Offline convective mixing/scavenging
   USE Henrys_law_ConstantsMod, ONLY: get_HenrysLawCts

   implicit none

! !PUBLIC TYPES:
!
   PRIVATE
   PUBLIC  NI_GridComp       ! The NI object 
   PUBLIC  NI_GridComp1      ! Single instance NI object 

!
! !PUBLIC MEMBER FUNCTIONS:
!

   PUBLIC  NI_GridCompSetServices
   PUBLIC  NI_GridCompInitialize
   PUBLIC  NI_GridCompRun1
   PUBLIC  NI_GridCompRun2
   PUBLIC  NI_GridCompFinalize

!
! !DESCRIPTION:
!
!  This module implements the (pre-ESMF) NI Grid Component. 
!
! !REVISION HISTORY:
!
!  24Sep2014 - Colarco, first crack
!  Code here is based on existing GOCART components with Nitrate specific
!  functions provided by Huisheng Bian.  Turn on "doing_NI" in Chem_Registry.
!  As default, we are running 5 tracers under Nitrate (see specification of
!  "global"tracer indices below):
!   - NH3    - Gas phase ammonia
!   - NH4    - Aerosol phase ammonium ion
!   - NO3an1 - Aerosol nitrate (radius <= 0.5 um)
!   - NO3an2 - Aerosol nitrate (0.5 < radius < 4 um)
!   - NO3an3 - Aerosol nitrate (4 < radius < 10 um)
!  Needed inputs are mixing ratio of nitric acid (HNO3) and emissions of
!  ammonia (NH3).  At present we take input of HNO3 from offline file, but
!  an unanswered question is how important it is to couple to HNO3 from
!  chemistry module.  Proper coupling also requires running GOCART
!  sulfate, dust, and sea salt (other species neglected at present)
!  with compatible assumptions of particle sizes (e.g., typical 5 dust
!  size bins and 5 sea salt size bins).
!  Pathway is:
!   1) emit HN3
!   2) Call thermodynamic module RPMARES to compute updates to 
!      NH3, NH4, and NO3an1 (and potentially update HNO3, SO4, etc.)
!   3) Call heterogeneous reactions to update NO3anX species
!   4) Loss processes via dry deposition, sedimentation, and convective
!      and wet removal
!
!  Parameters
   real*8, parameter :: gamma_seasalt = 0.2
!EOP
!-------------------------------------------------------------------------

  type NI_GridComp1
        character(len=255) :: name
        character(len=255) :: iname           ! instance name
        character(len=255) :: rcfilen         ! resource file name
        character(len=255) :: regionsString   ! Comma-delimited string of regions
        real, pointer      :: regionMask(:,:) ! regional mask
        integer :: instance                   ! instance number
        logical :: run_alarm = .false.        ! run alarm
        type(Chem_Mie), pointer :: mie_tables => null() ! aod LUTs
        real, pointer :: radius(:) !particle effective radius [um]
        real, pointer :: rhop(:)   ! NI class density [kg m-3]
        !real, pointer :: hno3(:,:,;)
        real, pointer :: xhno3(:,:,:)
        logical       :: first
        logical       :: recycle_HNO3 = .false.
  end type NI_GridComp1

  type NI_GridComp
     integer                     :: n = 0                ! number of instances 
     type(Chem_Mie), pointer     :: mie_tables => null() ! aod LUTs
     type(NI_GridComp1), pointer :: gcs(:)     => null() ! instances
  end type NI_GridComp

! Tracer assignments
  integer, parameter :: globalnNH3    = 1
  integer, parameter :: globalnNH4a   = 2
  integer, parameter :: globalnNO3an1 = 3
  integer, parameter :: globalnNO3an2 = 4
  integer, parameter :: globalnNO3an3 = 5
  real, parameter :: OCEAN=0.0, LAND = 1.0, SEA_ICE = 2.0
  real, parameter :: radToDeg = 57.2957795
  real, parameter :: fMassHNO3 = 63., fMassNO3 = 62., fMassAir = 29.

CONTAINS

 subroutine NI_GridCompSetServices(  gc, chemReg, rc)
   type(ESMF_GridComp), intent(INOUT) :: GC
   type(Chem_Registry), intent(INOUT) :: chemReg
   integer,             intent(OUT  ) :: rc

   CHARACTER(LEN=255) :: rcbasen = 'NI_GridComp'
   CHARACTER(LEN=255) :: name

   integer            :: status
   character(len=ESMF_MAXSTR) :: Iam

   integer :: n,i

   type(ESMF_Config) :: cfg

   Iam = "NI_GridCompSetServices"

!  Load resource file
!  ------------------
   cfg = ESMF_ConfigCreate(rc=status)
   VERIFY_(STATUS)
   call ESMF_ConfigLoadFile(cfg,TRIM(rcbasen)//'.rc',rc=status)
   VERIFY_(STATUS)

!  Parse resource file
!  -------------------
   n = ESMF_ConfigGetLen(cfg,label='NI_instances:',rc=status)
   VERIFY_(STATUS)


!  We have 5 tracers for each instance of BC
!  We cannot have fewer instances than half the number of
!   BC bins in the registry (it is OK to have less, though)
!  --------------------------------------------------------
   if ( 5*n .LT. chemReg%n_NI ) then
        rc = 35
        return
   else if ( 5*n .GT. chemReg%n_NI ) then
        if (MAPL_AM_I_ROOT()) &
        print *, trim(Iam)// &
                 ': fewer NI bins than possible NI instances: ',&
                 n, chemReg%n_NI/5
   end if
   n = min(n,chemReg%n_NI/5)

!  Record name of each instance
!  ----------------------------
   call ESMF_ConfigFindLabel(cfg,'NI_instances:',rc=status)
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

      call NI_GridCompSetServices1_(gc,chemReg,name,rc=status)
      VERIFY_(STATUS)
   end do

   call MAPL_AddImportSpec(GC,           &
        SHORT_NAME = 'NI_regionMask',    &
        LONG_NAME  = 'source species'  , &
        UNITS      = '1',                &
        DIMS       = MAPL_DimsHorzOnly,  &
        VLOCATION  = MAPL_VLocationNone, &
        RESTART    = MAPL_RestartSkip,   &
        RC         = STATUS)
   VERIFY_(STATUS)

   RETURN_(ESMF_SUCCESS)
 end subroutine NI_GridCompSetServices


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  NI_GridCompInitialize --- Initialize NI_GridComp
!
! !INTERFACE:
!

   subroutine NI_GridCompInitialize ( gcNI, w_c, impChem, expChem, &
                                      nymd, nhms, cdt, rc )

! !USES:

  implicit NONE

! !INPUT PARAMETERS:

   type(Chem_Bundle), intent(inout) :: w_c     ! Chemical tracer fields      
   integer, intent(in) :: nymd, nhms           ! time
   real, intent(in) :: cdt                     ! chemistry timestep (secs)

! !OUTPUT PARAMETERS:

   type(NI_GridComp), intent(inout) :: gcNI     ! Grid Component
   type(ESMF_State), intent(inout)  :: impChem  ! Import State
   type(ESMF_State), intent(inout)  :: expChem  ! Export State
   integer, intent(out) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 - 

! !DESCRIPTION: Initializes the NI Grid Component. It primarily sets
!               the import state for each active constituent package.
!
! !REVISION HISTORY:
!
!  26Aug2014 Colarco  First crack.
!
!EOP
!-------------------------------------------------------------------------

   character(len=*), parameter :: myname = 'NI_GridCompInitialize'
   character(len=255) :: rcbasen = 'NI_GridComp'
   CHARACTER(LEN=255) :: name
   
   integer i, ier, n,i_

!  Load resource file
!  ------------------
   call i90_loadf ( trim(rcbasen)//'.rc', ier )
   if ( ier .NE. 0 ) then
      rc = 10
      return
   end if

!  Parse resource file
!  -------------------
   CALL I90_label ( 'NI_instances:', ier )
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
   
!  We have 5 tracers for each instance of NI
!  Chem_Registry provides the number (total)
!  of tracers to be run.  Therefore n*5 must
!  be >= to that number or else we don't have
!  enough instances requested.
!  --------------------------------------------------------
   if ( n*5 .lt. w_c%reg%n_NI ) then
        rc = 35
        return
   end if
   n = min(n,w_c%reg%n_NI/5 )
   gcNI%n = n

!  Next allocate necessary memory
!  ------------------------------
   allocate ( gcNI%gcs(n), stat=ier )    
   if ( ier .NE. 0 ) then
      rc = 40
      return
   end if

!  Record name of each instance
!  ----------------------------
   CALL I90_label ( 'NI_instances:', ier )
   do i = 1, n
      CALL I90_gtoken( name, ier )
      if ( ier .NE. 0 ) then
         rc = 40
         return
      end if
                                            ! resource file name
      gcNI%gcs(i)%rcfilen = trim(rcbasen)//'---'//trim(name)//'.rc'
      gcNI%gcs(i)%instance = i              ! instance number 
      IF(TRIM(name) == "full" ) THEN
       gcNI%gcs(i)%iname = " "              ! blank instance name for full (1)
      ELSE
       gcNI%gcs(i)%iname = TRIM(name)       ! instance name for others
      END IF
   end do    

!  Next initialize each instance
!  -----------------------------
   do i = 1, gcNI%n
      IF(MAPL_AM_I_ROOT()) THEN
       PRINT *," "
       PRINT *,myname,": Initializing instance ",TRIM(gcNI%gcs(i)%iname)," [",gcNI%gcs(i)%instance,"]"
      END IF
      call NI_SingleInstance_ ( NI_GridCompInitialize1_, i, &
                                gcNI%gcs(i), w_c, impChem, expChem,  &
                                nymd, nhms, cdt, ier )
      if ( ier .NE. 0 ) then
         rc = 1000+ier
         return
      end if
      gcNI%gcs(i)%mie_tables => gcNI%mie_tables
   end do

!  Get Henrys Law cts for the parameterized convective wet removal
!  -----------------------------------------------------------
   do i = 1, gcNI%n
      !- NH3    
      i_ = w_c%reg%i_NI  + 4*(i - 1)
      CALL get_HenrysLawCts('NH3',w_c%reg%Hcts(1,i_),w_c%reg%Hcts(2,i_)&
                                 ,w_c%reg%Hcts(3,i_),w_c%reg%Hcts(4,i_))  
      !IF(MAPL_AM_I_ROOT()) THEN
      !  print*,"NH3=",i,w_c%reg%Hcts(1,i_),w_c%reg%Hcts(2,i_),w_c%reg%Hcts(3,i_),w_c%reg%Hcts(4,i_)
      !  call FLUSH(6)
      !ENDIF
   end do

!  All done
!  --------
   CALL I90_FullRelease( ier )
   IF( ier /= 0 ) THEN
    PRINT *,myname,": I90_FullRelease not successful."
    rc = 40
   END IF


 end subroutine NI_GridCompInitialize

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  NI_GridCompRun1 --- Run NI_GridComp
!
! !INTERFACE:
!

   subroutine NI_GridCompRun1 ( gcNI, w_c, impChem, expChem, &
                                      nymd, nhms, cdt, rc )

! !USES:

  IMPLICIT NONE

! !INPUT PARAMETERS:

   TYPE(Chem_Bundle), intent(in) :: w_c        ! Chemical tracer fields      
   INTEGER, INTENT(IN) :: nymd, nhms           ! time
   REAL,    INTENT(IN) :: cdt                  ! chemical timestep (secs)


! !OUTPUT PARAMETERS:

   TYPE(NI_GridComp), INTENT(INOUT) :: gcNI     ! Grid Component
   TYPE(ESMF_State), INTENT(INOUT)  :: impChem  ! Import State
   TYPE(ESMF_State), INTENT(INOUT)  :: expChem  ! Export State
   INTEGER, INTENT(OUT) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 - 

! !DESCRIPTION: Runs the NI Grid Component. Multiple instance
!               version.
!
! !REVISION HISTORY:
!
!  27Feb2008  da Silva  Introduced multiple instances
!
!EOP
!-------------------------------------------------------------------------

   integer i, ier

   do i = 1, gcNI%n
      call NI_SingleInstance_ ( NI_GridCompRun1_, i, &
                                gcNI%gcs(i), w_c, impChem, expChem, &
                                nymd, nhms, cdt, ier )
      if ( ier .NE. 0 ) then
         rc = i * 1000+ier
         return
      end if
   end do

 end subroutine NI_GridCompRun1



!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  NI_GridCompRun2 --- Run NI_GridComp
!
! !INTERFACE:
!

   subroutine NI_GridCompRun2 ( gcNI, w_c, impChem, expChem, &
                                run_alarm, nymd, nhms, cdt, rc )

! !USES:

  IMPLICIT NONE

! !INPUT PARAMETERS:

   TYPE(Chem_Bundle), intent(in) :: w_c        ! Chemical tracer fields
   LOGICAL, INTENT(IN) :: run_alarm            ! run alarm
   INTEGER, INTENT(IN) :: nymd, nhms           ! time
   REAL,    INTENT(IN) :: cdt                  ! chemical timestep (secs)


! !OUTPUT PARAMETERS:

   TYPE(NI_GridComp), INTENT(INOUT) :: gcNI     ! Grid Component
   TYPE(ESMF_State), INTENT(INOUT)  :: impChem  ! Import State
   TYPE(ESMF_State), INTENT(INOUT)  :: expChem  ! Export State
   INTEGER, INTENT(OUT) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 - 

! !DESCRIPTION: Runs the NI Grid Component. Multiple instance
!               version.
!
! !REVISION HISTORY:
!
!  27Feb2008  da Silva  Introduced multiple instances
!
!EOP
!-------------------------------------------------------------------------

   integer i, ier

   do i = 1, gcNI%n
      gcNI%gcs(i)%run_alarm = run_alarm

      call NI_SingleInstance_ ( NI_GridCompRun2_, i, &
                                gcNI%gcs(i), w_c, impChem, expChem, &
                                nymd, nhms, cdt, ier )
      if ( ier .NE. 0 ) then
         rc = i * 1000+ier
         return
      end if
   end do

 end subroutine NI_GridCompRun2


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  NI_GridCompFinalize --- Initialize NI_GridComp
!
! !INTERFACE:
!

   subroutine NI_GridCompFinalize ( gcNI, w_c, impChem, expChem, &
                                      nymd, nhms, cdt, rc )

! !USES:

  IMPLICIT NONE

! !INPUT PARAMETERS:

   TYPE(Chem_Bundle), intent(in) :: w_c        ! Chemical tracer fields      
   INTEGER, INTENT(IN) :: nymd, nhms           ! time
   REAL,    INTENT(IN) :: cdt                  ! chemical timestep (secs)


! !OUTPUT PARAMETERS:

   TYPE(NI_GridComp), INTENT(INOUT) :: gcNI     ! Grid Component
   TYPE(ESMF_State), INTENT(INOUT)  :: impChem  ! Import State
   TYPE(ESMF_State), INTENT(INOUT)  :: expChem  ! Export State
   INTEGER, INTENT(OUT) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 - 

! !DESCRIPTION: Finalizes the NI Grid Component. Multiple instance
!               version.
!
! !REVISION HISTORY:
!
!  27Feb2008  da Silva  Introduced multiple instances
!
!EOP
!-------------------------------------------------------------------------

   integer i, ier

   do i = 1, gcNI%n
      call NI_SingleInstance_ ( NI_GridCompFinalize1_, i, &
                                gcNI%gcs(i), w_c, impChem, expChem, &
                                nymd, nhms, cdt, ier )
      if ( ier .NE. 0 ) then
         rc = i * 1000+ier
         return
      end if
   end do

   deallocate ( gcNI%gcs, stat=ier )    
   gcNI%n = -1

 end subroutine NI_GridCompFinalize


   subroutine NI_GridCompSetServices1_(  gc, chemReg, iname, rc)

! !USES:

  implicit NONE

! !INPUT PARAMETERS:
  
   type(ESMF_GridComp), intent(INOUT) :: GC
   type(Chem_Registry), intent(INOUT) :: chemReg
   character(len=*),    intent(IN   ) :: iname
   integer,             intent(OUT  ) :: rc

   integer :: Status
   character(len=ESMF_MAXSTR) :: Iam

   Iam ="NI_GridCOmpSetServices1_"

   call MAPL_AddImportSpec(GC, &
        SHORT_NAME = 'EMI_NH3_AG'//trim(iname), &
        LONG_NAME  = ''  , &
        UNITS      = 'kg m-2 s-1',       &
        DIMS       = MAPL_DimsHorzOnly,  &
        VLOCATION  = MAPL_VLocationNone, &
        RESTART    = MAPL_RestartSkip,   &
        RC         = STATUS)
   VERIFY_(STATUS)

   call MAPL_AddImportSpec(GC, &
        SHORT_NAME = 'EMI_NH3_BB'//trim(iname), &
        LONG_NAME  = ''  , &
        UNITS      = 'kg m-2 s-1',       &
        DIMS       = MAPL_DimsHorzOnly,  &
        VLOCATION  = MAPL_VLocationNone, &
        RESTART    = MAPL_RestartSkip,   &
        RC         = STATUS)
   VERIFY_(STATUS)

   call MAPL_AddImportSpec(GC, &
        SHORT_NAME = 'EMI_NH3_EN'//trim(iname), &
        LONG_NAME  = ''  , &
        UNITS      = 'kg m-2 s-1',       &
        DIMS       = MAPL_DimsHorzOnly,  &
        VLOCATION  = MAPL_VLocationNone, &
        RESTART    = MAPL_RestartSkip,   &
        RC         = STATUS)
   VERIFY_(STATUS)

   call MAPL_AddImportSpec(GC, &
        SHORT_NAME = 'EMI_NH3_IN'//trim(iname), &
        LONG_NAME  = ''  , &
        UNITS      = 'kg m-2 s-1',       &
        DIMS       = MAPL_DimsHorzOnly,  &
        VLOCATION  = MAPL_VLocationNone, &
        RESTART    = MAPL_RestartSkip,   &
        RC         = STATUS)
   VERIFY_(STATUS)
   
   call MAPL_AddImportSpec(GC, &
        SHORT_NAME = 'EMI_NH3_OC'//trim(iname), &
        LONG_NAME  = ''  , &
        UNITS      = 'kg m-2 s-1',       &
        DIMS       = MAPL_DimsHorzOnly,  &
        VLOCATION  = MAPL_VLocationNone, &
        RESTART    = MAPL_RestartSkip,   &
        RC         = STATUS)
   VERIFY_(STATUS)

   call MAPL_AddImportSpec(GC, &
        SHORT_NAME = 'EMI_NH3_RE'//trim(iname), &
        LONG_NAME  = ''  , &
        UNITS      = 'kg m-2 s-1',       &
        DIMS       = MAPL_DimsHorzOnly,  &
        VLOCATION  = MAPL_VLocationNone, &
        RESTART    = MAPL_RestartSkip,   &
        RC         = STATUS)
   VERIFY_(STATUS)

   call MAPL_AddImportSpec(GC, &
        SHORT_NAME = 'EMI_NH3_TR'//trim(iname), &
        LONG_NAME  = ''  , &
        UNITS      = 'kg m-2 s-1',       &
        DIMS       = MAPL_DimsHorzOnly,  &
        VLOCATION  = MAPL_VLocationNone, &
        RESTART    = MAPL_RestartSkip,   &
        RC         = STATUS)
   VERIFY_(STATUS)

   call MAPL_AddImportSpec(GC, &
        SHORT_NAME = 'NITRATE_HNO3'//trim(iname), &
        LONG_NAME  = ''  , &
        UNITS      = '',                   &
        DIMS       = MAPL_DimsHorzVert,    &
        VLOCATION  = MAPL_VLocationCenter, &
        RESTART    = MAPL_RestartSkip,     &
        RC         = STATUS)
   VERIFY_(STATUS)

   RETURN_(ESMF_SUCCESS)

 end subroutine NI_GridCompSetServices1_


!--------------------------------------------------------------------------

!                      Single Instance Methods

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  NI_GridCompInitialize --- Initialize NI_GridComp
!
! !INTERFACE:
!

   subroutine NI_GridCompInitialize1_ ( gcNI, w_c, impChem, expChem, &
                                      nymd, nhms, cdt, rc )

! !USES:

  implicit NONE

! !INPUT PARAMETERS:

   type(Chem_Bundle), intent(inout) :: w_c     ! Chemical tracer fields      
   integer, intent(in) :: nymd, nhms           ! time
   real, intent(in) :: cdt                     ! chemistry timestep (secs)

! !OUTPUT PARAMETERS:

   type(NI_GridComp1), intent(inout) :: gcNI    ! Grid Component
   type(ESMF_State), intent(inout)  :: impChem  ! Import State
   type(ESMF_State), intent(inout)  :: expChem  ! Export State
   integer, intent(out) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 - 

! !DESCRIPTION: Initializes the NI Grid Component. It primarily sets
!               the import state for each active constituent package.
!
! !REVISION HISTORY:
!
!  18Sep2003 da Silva  First crack.
!
!EOP
!-------------------------------------------------------------------------

   character(len=*), parameter :: myname = 'NI_GridCompInitialize1'


   character(len=255) :: rcfilen
   integer :: n
   integer :: i1, i2, im, j1, j2, jm, km, nbins, n1, n2
   integer, allocatable :: ier(:)
   real :: qmax, qmin
   LOGICAL :: NoRegionalConstraint 


   rcfilen = gcNI%rcfilen
   gcNI%name = 'NI Constituent Package'

!  Initialize local variables
!  --------------------------
   rc = 0
   i1 = w_c%grid%i1; i2 = w_c%grid%i2; im = w_c%grid%im
   j1 = w_c%grid%j1; j2 = w_c%grid%j2; jm = w_c%grid%jm
   km = w_c%grid%km
   nbins = w_c%reg%n_NI
   n1  = w_c%reg%i_NI
   n2  = w_c%reg%j_NI

   gcNI%first = .True.

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

!  Scavenging Efficiency
!  To be used in convtran.F90, this parameter
!  is the scavenging efficiency of the tracer [km -1]
!  ---------------
   call i90_label ( 'fscav:', ier(1) )
   do n = 1, nbins
      w_c%reg%fscav(n1+n-1) = i90_gfloat ( ier(n+1) )
   end do
   if ( any(ier(1:nbins+1) /= 0) ) then
      call final_(20)
      return
   end if
!                          -------

!  Particle radius [um]
!  To be used in settling code
!  ---------------
   call i90_label ( 'particle_radius:', ier(1) )
   do n = 1, nbins
      gcNI%radius(n)        = i90_gfloat ( ier(n+1) )
   end do
   if ( any(ier(1:nbins+1) /= 0) ) then
      call final_(30)
      return
   end if
!                          -------

!  Particle density
!  To be used in droplet activation code
!  ---------------
   call i90_label ( 'particle_density:', ier(1) )
   do n = 1, nbins
      w_c%reg%rhop(n1+n-1)  = i90_gfloat ( ier(n+1) )
      gcNI%rhop(n)          = w_c%reg%rhop(n1+n-1)
   end do
   if ( any(ier(1:nbins+1) /= 0) ) then
      call final_(40)
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
      call final_(50)
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
      call final_(60)
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
      call final_(70)
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
      call final_(80)
      return
   end if
!                          -------

!  Grab the region string.
!  -----------------------
   call i90_label ( 'NI_regions_indices:', ier(1) )
   CALL I90_gtoken( gcNI%regionsString, ier(2) )
   IF( ANY(ier(1:2) < 0 ) ) THEN
    CALL final_(90)
    RETURN
   END IF

!  Is this instantiation a global case?
!  -----------------------------------
   IF(gcNI%regionsString(1:2) == "-1") THEN
    NoRegionalConstraint = .TRUE.
   ELSE
    SELECT CASE (ESMF_UtilStringLowerCase(gcNI%regionsString(1:2)))
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
   IF(NoRegionalConstraint) gcNI%regionsString = "-1"

   IF(MAPL_AM_I_ROOT()) THEN
    IF(NoRegionalConstraint) THEN
     PRINT *,myname,": This instantiation has no regional constraints."
    ELSE
     PRINT *,myname,": This instantiation is regionally constrained."
     PRINT *,myname,": List of region numbers included: ",TRIM(gcNI%regionsString)
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
   allocate( gcNI%regionmask(i1:i2,j1:j2), ier(nerr), &
             gcNI%radius(nbins), gcNI%rhop(nbins), &
             gcNI%xhno3(i1:i2,j1:j2,km), stat=ios )
   if ( ios /= 0 ) rc = 100
   end subroutine init_

   subroutine final_(ierr)
   integer :: ierr
   integer ios
   deallocate( gcNI%regionmask, ier, &
               gcNI%radius, gcNI%rhop, &
               gcNI%xhno3, stat=ios )
   call i90_release()
   rc = ierr
   end subroutine final_

   end subroutine NI_GridCompInitialize1_

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  NI_GridCompRun1_ --- The Chem Driver, run phase 1 
!
! !INTERFACE:
!

   subroutine NI_GridCompRun1_ ( gcNI, w_c, impChem, expChem, &
                                 nymd, nhms, cdt, rc )

! !USES:

  implicit NONE

! !INPUT/OUTPUT PARAMETERS:

   type(NI_GridComp1), intent(inout) :: gcNI   ! Grid Component
   type(Chem_Bundle), intent(inout)  :: w_c    ! Chemical tracer fields   

! !INPUT PARAMETERS:

   type(ESMF_State), intent(inout) :: impChem  ! Import State
   integer, intent(in) :: nymd, nhms           ! time
   real, intent(in) :: cdt                     ! chemistry timestep (secs)

! !OUTPUT PARAMETERS:

   type(ESMF_State), intent(inout) :: expChem  ! Export State
   integer, intent(out) :: rc                  ! Error return code:
                                               !  0 - all is well
                                               !  1 -
 
! !DESCRIPTION: This routine implements the so-called NI Driver. That 
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

   character(len=*), parameter :: myname = 'NI_GridCompRun1_'
   character(len=*), parameter :: Iam = myname

   integer :: i1, i2, im, j1, j2, jm, nbins, n1, n2, km
   integer :: ijl, ijkl, ijk1l
   real :: qmax, qmin

!  Input fields from fvGCM
!  -----------------------
   real, pointer, dimension(:,:)   :: emi_nh3_ag, emi_nh3_en, emi_nh3_tr, &
                                      emi_nh3_oc, emi_nh3_in, emi_nh3_re, &
                                      emi_nh3_bb


   real, pointer    :: var2D(:,:) => null()

!  Tracer assignments (local)
   integer :: nNH3, nNH4a, nNO3an1, nNO3an2, nNO3an3


#define EXPORT        expChem
#define iNAME         TRIM(gcNI%iname)

#define ptrNH3EM      NH3_emis
#define ptrNIEM       NI_emis

   
   integer :: STATUS

#include "NI_GetPointer___.h"


!  Initialize local variables
!  --------------------------
   rc = 0
   i1 = w_c%grid%i1; i2 = w_c%grid%i2; im = w_c%grid%im
   j1 = w_c%grid%j1; j2 = w_c%grid%j2; jm = w_c%grid%jm

   km = w_c%grid%km
   nbins = w_c%reg%n_NI
   n1  = w_c%reg%i_NI
   n2  = w_c%reg%j_NI

   nNH3    = n1 + globalnNH3    - 1
   nNH4a   = n1 + globalnNH4a   - 1
   nNO3an1 = n1 + globalnNO3an1 - 1
   nNO3an2 = n1 + globalnNO3an2 - 1
   nNO3an3 = n1 + globalnNO3an3 - 1

   ijl  = ( i2 - i1 + 1 ) * ( j2 - j1 + 1 )
   ijkl = ijl * km
   ijk1l = ijl * (km+1)

   call MAPL_GetPointer(impChem, var2D, 'NI_regionMask', __RC__)
   gcNI%regionMask = var2D


#ifdef DEBUG
   do n = n1, n2
      call pmaxmin ( 'NI: q_beg', w_c%qa(n)%data3d(i1:i2,j1:j2,1:km), qmin, qmax, &
                   ijl, km, 1. )
   end do
#endif


!  Emissions
!  ---------
   call MAPL_GetPointer ( impChem, emi_nh3_bb, 'EMI_NH3_BB'//iNAME, __RC__ )
   call MAPL_GetPointer ( impChem, emi_nh3_ag, 'EMI_NH3_AG'//iNAME, __RC__ )
   call MAPL_GetPointer ( impChem, emi_nh3_en, 'EMI_NH3_EN'//iNAME, __RC__ )
   call MAPL_GetPointer ( impChem, emi_nh3_re, 'EMI_NH3_RE'//iNAME, __RC__ )
   call MAPL_GetPointer ( impChem, emi_nh3_tr, 'EMI_NH3_TR'//iNAME, __RC__ )
   call MAPL_GetPointer ( impChem, emi_nh3_in, 'EMI_NH3_IN'//iNAME, __RC__ )
   call MAPL_GetPointer ( impChem, emi_nh3_oc, 'EMI_NH3_OC'//iNAME, __RC__ )


!  NH3 Emissions
!  -------------
   if(associated(NH3_emis%data2d)) then 
    NH3_emis%data2d = 0.
    if(associated(emi_nh3_bb)) NH3_emis%data2d = NH3_emis%data2d + emi_nh3_bb
    if(associated(emi_nh3_ag)) NH3_emis%data2d = NH3_emis%data2d + emi_nh3_ag
    if(associated(emi_nh3_en)) NH3_emis%data2d = NH3_emis%data2d + emi_nh3_en
    if(associated(emi_nh3_tr)) NH3_emis%data2d = NH3_emis%data2d + emi_nh3_tr
    if(associated(emi_nh3_re)) NH3_emis%data2d = NH3_emis%data2d + emi_nh3_re
    if(associated(emi_nh3_in)) NH3_emis%data2d = NH3_emis%data2d + emi_nh3_in
    if(associated(emi_nh3_oc)) NH3_emis%data2d = NH3_emis%data2d + emi_nh3_oc
   endif

   if(associated(emi_nh3_bb)) &
    w_c%qa(nNH3)%data3d(:,:,km) = w_c%qa(nNH3)%data3d(:,:,km) &
                                + cdt * grav / w_c%delp(:,:,km) * emi_nh3_bb
   if(associated(emi_nh3_ag)) &
    w_c%qa(nNH3)%data3d(:,:,km) = w_c%qa(nNH3)%data3d(:,:,km) &
                                + cdt * grav / w_c%delp(:,:,km) * emi_nh3_ag
   if(associated(emi_nh3_en)) &
    w_c%qa(nNH3)%data3d(:,:,km) = w_c%qa(nNH3)%data3d(:,:,km) &
                                + cdt * grav / w_c%delp(:,:,km) * emi_nh3_en
   if(associated(emi_nh3_in)) &
    w_c%qa(nNH3)%data3d(:,:,km) = w_c%qa(nNH3)%data3d(:,:,km) &
                                + cdt * grav / w_c%delp(:,:,km) * emi_nh3_in
   if(associated(emi_nh3_re)) &
    w_c%qa(nNH3)%data3d(:,:,km) = w_c%qa(nNH3)%data3d(:,:,km) &
                                + cdt * grav / w_c%delp(:,:,km) * emi_nh3_re
   if(associated(emi_nh3_tr)) &
    w_c%qa(nNH3)%data3d(:,:,km) = w_c%qa(nNH3)%data3d(:,:,km) &
                                + cdt * grav / w_c%delp(:,:,km) * emi_nh3_tr
   if(associated(emi_nh3_oc)) &
    w_c%qa(nNH3)%data3d(:,:,km) = w_c%qa(nNH3)%data3d(:,:,km) &
                                + cdt * grav / w_c%delp(:,:,km) * emi_nh3_oc

   return

 end subroutine NI_GridCompRun1_



!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  NI_GridCompRun2 --- The Chem Driver, run phase 2
!
! !INTERFACE:
!

   subroutine NI_GridCompRun2_ ( gcNI, w_c, impChem, expChem, &
                                 nymd, nhms, cdt, rc )

! !USES:

  implicit NONE

! !INPUT/OUTPUT PARAMETERS:

   type(NI_GridComp1), intent(inout) :: gcNI   ! Grid Component
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
 
! !DESCRIPTION: This routine implements the so-called NI Driver. That 
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

   character(len=*), parameter :: myname = 'NI_GridCompRun2_'
   character(len=*), parameter :: Iam = myname

   integer :: i1, i2, im, j1, j2, jm, nbins, n1, n2, km, n, ios
   integer :: i, j, k, ijl, ijkl, ijk1l
   real :: qmax, qmin
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
   real, pointer, dimension(:,:,:) :: hno3

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

   real             :: NI_radius, NI_rhop
   integer          :: rhFlag

   real, pointer    :: var2D(:,:) => null()

!  Tracer assignments (local)
   integer :: nNH3, nNH4a, nNO3an1, nNO3an2, nNO3an3, nSO4, na

!  variables for call to RPMARES
   real   :: fmmr_to_conc
   real*8 :: SO4, GNO3, GNH3, RH, TEMP, ASO4, AHSO4, AH2O, ANO3, ANH4

!  variables for call to heterogeneous chemistry
   real*8 :: kan1, kan2, kan3, sad, ad, rad, deltahno3

   character(len=5), allocatable :: duname(:)
   character(len=5), allocatable :: ssname(:)

#define EXPORT        expChem
#define iNAME         TRIM(gcNI%iname)

#define ptrNH3EM      NH3_emis
#define ptrNH3WT      NH3_wet
#define ptrNH3SV      NH3_conv
#define ptrNH3DP      NH3_dep
#define ptrNH3MASS    NH3_mass
#define ptrNH4WT      NH4_wet
#define ptrNH4SV      NH4_conv
#define ptrNH4DP      NH4_dep
#define ptrNH4SD      NH4_set
#define ptrNH4MASS    NH4_mass
#define ptrNIPNO3AQ   NI_pno3aq
#define ptrNIPNH4AQ   NI_pnh4aq
#define ptrNIPNH3AQ   NI_pnh3aq
#define ptrNIWT       NI_wet
#define ptrNISV       NI_conv
#define ptrNIEM       NI_emis
#define ptrNIDP       NI_dep
#define ptrNISD       NI_set
#define ptrNIHT       NI_phet

#define ptrHNO3SMASS  HNO3_sfcmass
#define ptrHNO3CMASS  HNO3_colmass
#define ptrHNO3CONC   HNO3_conc
#define ptrNH3SMASS   NH3_sfcmass
#define ptrNH3CMASS   NH3_colmass
#define ptrNH3CONC    NH3_conc
#define ptrNH4SMASS   NH4_sfcmass
#define ptrNH4CMASS   NH4_colmass
#define ptrNH4CONC    NH4_conc
#define ptrNISMASS25  NI_sfcmass25
#define ptrNICMASS25  NI_colmass25
#define ptrNISMASS    NI_sfcmass
#define ptrNICMASS    NI_colmass
#define ptrNIEXTT25   NI_exttau25
#define ptrNISCAT25   NI_scatau25
#define ptrNIEXTTFM   NI_exttaufm
#define ptrNISCATFM   NI_scataufm
#define ptrNIEXTTAU   NI_exttau
#define ptrNISCATAU   NI_scatau
#define ptrNIMASS25   NI_mass25
#define ptrNICONC25   NI_conc25
#define ptrNIMASS     NI_mass
#define ptrNICONC     NI_conc
#define ptrNIEXTCOEF  NI_extcoef
#define ptrNISCACOEF  NI_scacoef
#define ptrNIANGSTR   NI_angstrom
#define ptrNIFLUXU    NI_fluxu
#define ptrNIFLUXV    NI_fluxv


   
   integer :: STATUS

#include "NI_GetPointer___.h"


!  Initialize local variables
!  --------------------------
   rc = 0
   i1 = w_c%grid%i1; i2 = w_c%grid%i2; im = w_c%grid%im
   j1 = w_c%grid%j1; j2 = w_c%grid%j2; jm = w_c%grid%jm

   km = w_c%grid%km
   nbins = w_c%reg%n_NI
   n1  = w_c%reg%i_NI
   n2  = w_c%reg%j_NI

   nNH3    = n1 + globalnNH3    - 1
   nNH4a   = n1 + globalnNH4a   - 1
   nNO3an1 = n1 + globalnNO3an1 - 1
   nNO3an2 = n1 + globalnNO3an2 - 1
   nNO3an3 = n1 + globalnNO3an3 - 1

   ijl  = ( i2 - i1 + 1 ) * ( j2 - j1 + 1 )
   ijkl = ijl * km
   ijk1l = ijl * (km+1)

   call MAPL_GetPointer(impChem, var2D, 'NI_regionMask', __RC__)
   gcNI%regionMask = var2D


#ifdef DEBUG
   do n = n1, n2
      call pmaxmin ( 'NI: q_beg', w_c%qa(n)%data3d(i1:i2,j1:j2,1:km), qmin, qmax, &
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

#ifdef DEBUG

   call pmaxmin('NI: frlake     ', frlake  , qmin, qmax, ijl,1, 1. )
   call pmaxmin('NI: frocean    ', frocean , qmin, qmax, ijl,1, 1. )
   call pmaxmin('NI: frseaice   ', frseaice, qmin, qmax, ijl,1, 1. )
   call pmaxmin('NI: area       ', area    , qmin, qmax, ijl,1, 1. )
   call pmaxmin('NI: oro        ', oro     , qmin, qmax, ijl,1, 1. )
   call pmaxmin('NI: u10m       ', u10m    , qmin, qmax, ijl,1, 1. )
   call pmaxmin('NI: v10m       ', v10m    , qmin, qmax, ijl,1, 1. )
   call pmaxmin('NI: ustar      ', ustar   , qmin, qmax, ijl,1, 1. )
   call pmaxmin('NI: precc      ', precc   , qmin, qmax, ijl,1, 1. )
   call pmaxmin('NI: precl      ', precl   , qmin, qmax, ijl,1, 1. )
   call pmaxmin('NI: pblh       ', pblh    , qmin, qmax, ijl,1, 1. )
   call pmaxmin('NI: shflux     ', shflux  , qmin, qmax, ijl,1, 1. )
   call pmaxmin('NI: z0h        ', z0h     , qmin, qmax, ijl,1, 1. )
   call pmaxmin('NI: hsurf      ', hsurf   , qmin, qmax, ijl,1, 1. )

   call pmaxmin('NI: tmpu       ', tmpu    , qmin, qmax, ijkl,1, 1. )
   call pmaxmin('NI: rhoa       ', rhoa    , qmin, qmax, ijkl,1, 1. )
   call pmaxmin('NI: u          ', u       , qmin, qmax, ijkl,1, 1. )
   call pmaxmin('NI: v          ', v       , qmin, qmax, ijkl,1, 1. )
   call pmaxmin('NI: hghte      ', hghte   , qmin, qmax, ijkl,1, 1. )
   call pmaxmin('NI: qlcn       ', qlcn    , qmin, qmax, ijkl,1, 1. )
   call pmaxmin('NI: qicn       ', qicn    , qmin, qmax, ijkl,1, 1. )
   call pmaxmin('NI: cmfmc      ', cmfmc   , qmin, qmax, ijkl,1, 1. )
   call pmaxmin('NI: dtrain     ', dtrain  , qmin, qmax, ijkl,1, 1. )
   call pmaxmin('NI: pfllsan    ', pfllsan , qmin, qmax, ijk1l,1, 1. )
   call pmaxmin('NI: pfilsan    ', pfilsan , qmin, qmax, ijk1l,1, 1. )

#endif



!  Nitric Acid
!  -----------
   call MAPL_GetPointer ( impChem, hno3, 'NITRATE_HNO3'//iNAME, __RC__ )

!  Save local copy of HNO3 for first pass through run method regardless
   if (gcNI%first) then
       gcNI%xhno3 = MAPL_UNDEF
       gcNI%first = .False.
   end if

   ! Recycle HNO3 every 3 hours
   if (gcNI%recycle_HNO3) then
       gcNI%xhno3 = hno3
       gcNI%recycle_HNO3 = .false.
   end if


RUN_ALARM: if (gcNI%run_alarm) then

   allocate( fluxout )
   allocate( fluxout%data2d(i1:i2,j1:j2), dqa(i1:i2,j1:j2), &
             drydepositionfrequency(i1:i2,j1:j2), stat=STATUS)
   VERIFY_(STATUS)


!  Unlike GEOS-4 hghte is defined for km+1
!  ---------------------------------------
   hsurf => hghte(i1:i2,j1:j2,km) ! Recall: GEOS-5 has edges with k in [0,km]
    

!  Nitrate Chemistry
!  -----------------
   if(associated(NI_pno3aq%data2d)) NI_pno3aq%data2d(:,:) = 0.
   if(associated(NI_pnh4aq%data2d)) NI_pnh4aq%data2d(:,:) = 0.
   if(associated(NI_pnh3aq%data2d)) NI_pnh3aq%data2d(:,:) = 0.


!  RPMARES - thermodynamic module
!  ------------------------------
!  Take as input GOCART provided SO4, model provided RH,
!  and HNO3, NH3, NH4, and fine-mode nitrate (NO3an1).
!  At present we update NH3, NH4, and NO3an1.
!  Check we are running GOCART sulfate

   nSO4 = -1
   if(w_c%reg%doing_SU) then
    do n = w_c%reg%i_SU, w_c%reg%j_SU
     if(trim(w_c%reg%vname(n)) .eq. 'SO4') nSO4 = n
    enddo
   endif

   do k = 1, km
    do j = j1, j2
     do i = i1, i2

!     Conversion of mass mixing ratio to concentration (ug m-3)
      fmmr_to_conc = 1.e9 * rhoa(i,j,k)

!     Unit conversion for input to thermodynamic module
!     Per grid box call to RPMARES thermodynamic module
!     We do not presently treat chemistry of sulfate completely,
!     hence we ignore terms for ASO4, AHSO4, AH2O, and we do
!     not update SO4 on output from RPMARES.
!     At present we are importing HNO3 from offline file, so we
!     do not update on return.
      SO4 = 1.d-32
      if(nSO4 > 0) SO4  = max(1.d-32,w_c%qa(nSO4)%data3d(i,j,k) * fmmr_to_conc)
      GNO3              = max(1.d-32,gcNI%xhno3(i,j,k) * fMassHNO3 / fMassAir * fmmr_to_conc)
      GNH3              = max(1.d-32,w_c%qa(nNH3)%data3d(i,j,k)  * fmmr_to_conc)
      RH                = w_c%rh(i,j,k)
      TEMP              = tmpu(i,j,k)
      ASO4              = 1.d-32
      AHSO4             = 1.d-32
      ANO3              = max(1.d-32,w_c%qa(nNO3an1)%data3d(i,j,k) * fmmr_to_conc)
      AH2O              = 1.d-32
      ANH4              = max(1.d-32,w_c%qa(nNH4a)%data3d(i,j,k) * fmmr_to_conc)

      call RPMARES (  SO4,  GNO3,  GNH3, RH,   TEMP, &
                      ASO4, AHSO4, ANO3, AH2O, ANH4 )
!if(mapl_am_i_root()) print*,'RPMARES GNH3 = ',GNH3
!if(mapl_am_i_root())print*,'RPMARES ANO3 = ',ANO3
!if(mapl_am_i_root())print*,'RPMARES ANH4 = ',ANH4

!     Diagnostic terms
      if(associated(NI_pno3aq%data2d)) &
       NI_pno3aq%data2d(i,j) = NI_pno3aq%data2d(i,j) &
        + (ANO3 / fmmr_to_conc - w_c%qa(nNO3an1)%data3d(i,j,k)) &
          * w_c%delp(i,j,k)/grav/cdt
      if(associated(NI_pnh4aq%data2d)) &
       NI_pnh4aq%data2d(i,j) = NI_pnh4aq%data2d(i,j) &
        + (ANH4 / fmmr_to_conc - w_c%qa(nNH4a)%data3d(i,j,k)) &
          * w_c%delp(i,j,k)/grav/cdt
      if(associated(NI_pnh3aq%data2d)) &
       NI_pnh3aq%data2d(i,j) = NI_pnh3aq%data2d(i,j) &
        + (GNH3 / fmmr_to_conc - w_c%qa(nNH3)%data3d(i,j,k)) &
          * w_c%delp(i,j,k)/grav/cdt

!     Unit conversion back on return from thermodynamic module
      w_c%qa(nNH3)%data3d(i,j,k)    = GNH3 / fmmr_to_conc
      w_c%qa(nNO3an1)%data3d(i,j,k) = ANO3 / fmmr_to_conc
      w_c%qa(nNH4a)%data3d(i,j,k)   = ANH4 / fmmr_to_conc
      gcNI%xhno3(i,j,k) = max(1.d-32, GNO3 * fMassAir / fMassHNO3 / fmmr_to_conc)

     enddo
    enddo
   enddo

   ! prepare the variable names for comparison
   if(w_c%reg%doing_DU) then
      allocate(duname(w_c%reg%i_DU:w_c%reg%j_DU))
      do n = w_c%reg%i_DU, w_c%reg%j_DU
        duname(n) = ESMF_UtilStringUpperCase(trim(w_c%reg%vname(n)))
      end do
   end if

   if(w_c%reg%doing_SS) then
      allocate(ssname(w_c%reg%i_SS:w_c%reg%j_SS))
      do n = w_c%reg%i_SS, w_c%reg%j_SS
        ssname(n) = ESMF_UtilStringUpperCase(trim(w_c%reg%vname(n)))
      end do
   end if

!  Heterogeneous chemistry
!  -----------------------
!  Heterogeneous chemistry wants to know about GOCART dust and sea
!  salt tracers.  This code is not at the moment generalized as it
!  seems very wedded to the traditional GOCART arrangement (5 dust,
!  5 sea salt) and the particulars of the nitrate aerosol arrangement.
   if(associated(NI_phet(1)%data2d)) NI_phet(1)%data2d = 0.
   if(associated(NI_phet(2)%data2d)) NI_phet(2)%data2d = 0.
   if(associated(NI_phet(3)%data2d)) NI_phet(3)%data2d = 0.
   do k = 1, km
    do j = j1, j2
     do i = i1, i2
      kan1 = 0.
      kan2 = 0.
      kan3 = 0.
      ad = 1.e-6*rhoa(i,j,k)*MAPL_AVOGAD/MAPL_AIRMW  ! air number density # cm-3
      temp = tmpu(i,j,k)
      rh = w_c%rh(i,j,k)
!     Dust
      if(w_c%reg%doing_DU) then
       do n = w_c%reg%i_DU, w_c%reg%j_DU
        sad = 0.01*4.*MAPL_PI*w_c%reg%rmed(n)**2.*w_c%reg%fnum(n) * &
              rhoa(i,j,k) * w_c%qa(n)%data3d(i,j,k)       ! surface area density cm2 cm-3
        rad = 100.*w_c%reg%rmed(n)                        ! radius cm

        if (sad > 0.) then
           if(duname(n) .eq. 'DU001') &
            kan1 = kan1 + sktrs_hno3(temp,rh,sad,ad,rad)
           if(duname(n) .eq. 'DU002') &
            kan2 = kan2 + sktrs_hno3(temp,rh,sad,ad,rad)
           if(duname(n) .eq. 'DU003') &
            kan2 = kan2 + sktrs_hno3(temp,rh,sad,ad,rad)
           if(duname(n) .eq. 'DU004') &
            kan3 = kan3 + sktrs_hno3(temp,rh,sad,ad,rad)
           if(duname(n) .eq. 'DU005') &
            kan3 = kan3 + sktrs_hno3(temp,rh,sad,ad,rad)
        end if

       enddo
      endif

!     Sea salt
      if(w_c%reg%doing_SS) then
       do n = w_c%reg%i_SS, w_c%reg%j_SS
        sad = 0.01*4.*MAPL_PI*w_c%reg%rmed(n)**2.*w_c%reg%fnum(n) * &
              rhoa(i,j,k) * w_c%qa(n)%data3d(i,j,k)       ! surface area density cm2 cm-3
        rad = 100.*w_c%reg%rmed(n)                        ! radius cm

        if (sad > 0.) then
           if(ssname(n) .eq. 'SS001') &
            kan1 = kan1 + sktrs_sslt(temp,rh,sad,ad,rad)
           if(ssname(n) .eq. 'SS002') &
            kan1 = kan1 + sktrs_sslt(temp,rh,sad,ad,rad)
           if(ssname(n) .eq. 'SS003') &
            kan2 = kan2 + sktrs_sslt(temp,rh,sad,ad,rad)
           if(ssname(n) .eq. 'SS004') &
            kan2 = kan2 + sktrs_sslt(temp,rh,sad,ad,rad)
           if(ssname(n) .eq. 'SS005') &
            kan3 = kan3 + sktrs_sslt(temp,rh,sad,ad,rad)
        end if

       enddo
      endif

!     Compute the nitric acid loss (but don't actually update)
      if( (kan1+kan2+kan3) > 0.) then
       deltahno3 = gcNI%xhno3(i,j,k) * fMassHNO3 / fMassAir * (1.-exp(-(kan1+kan2+kan3)*cdt))
       gcNI%xhno3(i,j,k) = gcNI%xhno3(i,j,k) - deltahno3 * fMassAir / fMassHNO3
       w_c%qa(nNO3an1)%data3d(i,j,k) = &
         w_c%qa(nNO3an1)%data3d(i,j,k) + kan1/(kan1+kan2+kan3)*deltahno3*fMassNO3/fMassHNO3
       w_c%qa(nNO3an2)%data3d(i,j,k) = &
         w_c%qa(nNO3an2)%data3d(i,j,k) + kan2/(kan1+kan2+kan3)*deltahno3*fMassNO3/fMassHNO3
       w_c%qa(nNO3an3)%data3d(i,j,k) = &
         w_c%qa(nNO3an3)%data3d(i,j,k) + kan3/(kan1+kan2+kan3)*deltahno3*fMassNO3/fMassHNO3
       if(associated(NI_phet(1)%data2d)) &
         NI_phet(1)%data2d(i,j) = NI_phet(1)%data2d(i,j) + kan1/(kan1+kan2+kan3)*deltahno3*w_c%delp(i,j,k)/grav/cdt
       if(associated(NI_phet(2)%data2d)) &
         NI_phet(2)%data2d(i,j) = NI_phet(2)%data2d(i,j) + kan2/(kan1+kan2+kan3)*deltahno3*w_c%delp(i,j,k)/grav/cdt
       if(associated(NI_phet(3)%data2d)) &
         NI_phet(3)%data2d(i,j) = NI_phet(3)%data2d(i,j) + kan3/(kan1+kan2+kan3)*deltahno3*w_c%delp(i,j,k)/grav/cdt
      endif
     enddo
    enddo
   enddo   

   if(w_c%reg%doing_DU) deallocate(duname)
   if(w_c%reg%doing_SS) deallocate(ssname)


!  Output diagnostic HNO3
!  ----------------------
!  Calculate the HNO3 mass concentration
   if( associated(HNO3_conc%data3d) ) then
      HNO3_conc%data3d = gcNI%xhno3 * fMassHNO3 / fMassAir*rhoa
   endif
!  Calculate the HNO3 surface mass concentration
   if( associated(HNO3_sfcmass%data2d) ) then
      HNO3_sfcmass%data2d(i1:i2,j1:j2) = gcNI%xhno3(i1:i2,j1:j2,km) * fMassHNO3 / fMassAir*rhoa(i1:i2,j1:j2,km)
   endif
!  Calculate the HNO3 column loading
   if( associated(HNO3_colmass%data2d) ) then
      HNO3_colmass%data2d(i1:i2,j1:j2) = 0.
      do k = 1, km
        HNO3_colmass%data2d(i1:i2,j1:j2) &
         =   HNO3_colmass%data2d(i1:i2,j1:j2) + gcNI%xhno3(i1:i2,j1:j2,k)*w_c%delp(i1:i2,j1:j2,k)/grav
      end do
   endif

!  NI Settling
!  -----------
!  Because different bins having different swelling coefficients I need to
!  handle the call to settling differently.

!  Ammonium - settles like ammonium sulfate (rhflag = 3)
   n = globalnNH4a
   rhflag = 3
   NI_radius = 1.e-6*gcNI%radius(n)   ! radius in [m]
   NI_rhop   = gcNI%rhop(n)
   call Chem_SettlingSimple ( i1, i2, j1, j2, km, nNH4a, rhFlag, &
                        NI_radius, NI_rhop, cdt, w_c, tmpu, rhoa, hsurf,    &
                        hghte, fluxout, rc )
   if(associated(NH4_set%data2d)) NH4_set%data2d = fluxout%data2d

!  Nitrate bin 1 - settles like ammonium sulfate (rhflag = 3)
   n = globalnNO3an1
   rhflag = 3
   NI_radius = 1.e-6*gcNI%radius(n)   ! radius in [m]
   NI_rhop   = gcNI%rhop(n)
   call Chem_SettlingSimple ( i1, i2, j1, j2, km, nNO3an1, rhFlag, &
                        NI_radius, NI_rhop, cdt, w_c, tmpu, rhoa, hsurf,    &
                        hghte, fluxout, rc )
   if(associated(NI_set(1)%data2d)) NI_set(1)%data2d = fluxout%data2d

!  Nitrate bin 2 - settles like sea salt (rhflag = 2)
   n = globalnNO3an2
   rhflag = 2
   NI_radius = 1.e-6*gcNI%radius(n)   ! radius in [m]
   NI_rhop   = gcNI%rhop(n)
   call Chem_SettlingSimple ( i1, i2, j1, j2, km, nNO3an2, rhFlag, &
                        NI_radius, NI_rhop, cdt, w_c, tmpu, rhoa, hsurf,    &
                        hghte, fluxout, rc )
   if(associated(NI_set(2)%data2d)) NI_set(2)%data2d = fluxout%data2d

!  Nitrate bin 3 - settles like dust (rhflag = 0)
   n = globalnNO3an3
   rhflag = 0
   NI_radius = 1.e-6*gcNI%radius(n)   ! radius in [m]
   NI_rhop   = gcNI%rhop(n)
   call Chem_SettlingSimple ( i1, i2, j1, j2, km, nNO3an3, rhFlag, &
                        NI_radius, NI_rhop, cdt, w_c, tmpu, rhoa, hsurf,    &
                        hghte, fluxout, rc )
   if(associated(NI_set(3)%data2d)) NI_set(3)%data2d = fluxout%data2d


!  NI Deposition
!  -----------
   drydepositionfrequency = 0.
   call DryDepositionGOCART( i1, i2, j1, j2, km, &
                             tmpu, rhoa, hghte, oro, ustar, &
                             pblh, shflux, z0h, drydepositionfrequency, rc )
    
   n = globalnNH3
   dqa = 0.
   where (abs(oro - OCEAN) < 0.5)
       dqa = max(0.0, w_c%qa(n1+n-1)%data3d(:,:,km)*(1.-exp(-10.0*drydepositionfrequency*cdt)))
   elsewhere
       dqa = max(0.0, w_c%qa(n1+n-1)%data3d(:,:,km)*(1.-exp( -3.0*drydepositionfrequency*cdt)))
   end where
   w_c%qa(n1+n-1)%data3d(:,:,km) = w_c%qa(n1+n-1)%data3d(:,:,km) - dqa
   if( associated(NH3_dep%data2d) ) NH3_dep%data2d = dqa*w_c%delp(:,:,km)/grav/cdt

   n = globalnNH4a
   dqa = 0.
   dqa = max(0.0, w_c%qa(n1+n-1)%data3d(:,:,km)*(1.-exp(-drydepositionfrequency*cdt)))
   w_c%qa(n1+n-1)%data3d(:,:,km) = w_c%qa(n1+n-1)%data3d(:,:,km) - dqa
   if( associated(NH4_dep%data2d) ) NH4_dep%data2d = dqa*w_c%delp(:,:,km)/grav/cdt

   do n = globalnNO3an1, globalnNO3an3
    dqa = 0.
    dqa = max(0.0, w_c%qa(n1+n-1)%data3d(:,:,km)*(1.-exp(-drydepositionfrequency*cdt)))
    w_c%qa(n1+n-1)%data3d(:,:,km) = w_c%qa(n1+n-1)%data3d(:,:,km) - dqa
    if( associated(NI_dep(n-2)%data2d) ) NI_dep(n-2)%data2d = dqa*w_c%delp(:,:,km)/grav/cdt
   end do

#ifdef DEBUG
   do n = n1, n2
      call pmaxmin('NI: q_dry', w_c%qa(n)%data3d(i1:i2,j1:j2,1:km), qmin, qmax, &
                    ijl, km, 1. )
   end do
#endif

!  NI Large-scale Wet Removal
!  --------------------------
   w_c%qa(nNH3)%fwet = 1.
   KIN = .FALSE.   ! treat ammonia as gas
   call WetRemovalGOCART(i1, i2, j1, j2, km, nNH3, nNH3, cdt, 'NH3', KIN, &
                         w_c%qa, ple, tmpu, rhoa, pfllsan, pfilsan, &
                         precc, precl, fluxout, rc )
   if(associated(NH3_wet%data2d)) NH3_wet%data2d = fluxout%data2d

   w_c%qa(nNH4a)%fwet = 1.
   KIN = .TRUE.
   call WetRemovalGOCART(i1, i2, j1, j2, km, nNH4a, nNH4a, cdt, 'NH4a', KIN, &
                         w_c%qa, ple, tmpu, rhoa, pfllsan, pfilsan, &
                         precc, precl, fluxout, rc )
   if(associated(NH4_wet%data2d)) NH4_wet%data2d = fluxout%data2d

   do n = nNO3an1, nNO3an3
    w_c%qa(n)%fwet = 1.
    if(n .eq. nNO3an3) w_c%qa(n)%fwet = 0.3  ! treat coarse mode like dust
    KIN = .TRUE.
    call WetRemovalGOCART(i1, i2, j1, j2, km, n, n, cdt, 'nitrate', KIN, &
                         w_c%qa, ple, tmpu, rhoa, pfllsan, pfilsan, &
                         precc, precl, fluxout, rc )
    na = n - n1 - 1
    if(associated(NI_wet(na)%data2d)) NI_wet(na)%data2d = fluxout%data2d
   end do

#ifdef DEBUG
   do n = n1, n2
      call pmaxmin('NI: q_wet', w_c%qa(n)%data3d(i1:i2,j1:j2,1:km), qmin, qmax, &
                    ijl, km, 1. )
   end do
#endif


!  Nitrate Convective-scale Mixing and Wet Removal
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
   call convection(i1, i2, j1, j2, km, n1, n2, icdt, 'nitrate', kin, &
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
   if(associated(NH3_conv%data2d)) NH3_conv%data2d = -bcnv_(:,:,nNH3)/area_/icdt
   if(associated(NH4_conv%data2d)) NH4_conv%data2d = -bcnv_(:,:,nNH4a)/area_/icdt
   if(associated(NI_conv(1)%data2d)) NI_conv(1)%data2d = -bcnv_(:,:,nNO3an1)/area_/icdt
   if(associated(NI_conv(2)%data2d)) NI_conv(2)%data2d = -bcnv_(:,:,nNO3an2)/area_/icdt
   if(associated(NI_conv(3)%data2d)) NI_conv(3)%data2d = -bcnv_(:,:,nNO3an3)/area_/icdt


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
   call NI_Compute_Diags(i1, i2, j1, j2, km, nbins, gcNI, w_c, tmpu, rhoa, u, v, &
                         NH3_sfcmass, NH3_colmass, NH3_mass, NH3_conc, &
                         NH4_sfcmass, NH4_colmass, NH4_mass, NH4_conc, &
                         NI_sfcmass, NI_colmass, NI_mass, NI_conc, &
                         NI_sfcmass25, NI_colmass25, NI_mass25, NI_conc25, &
                         NI_exttau,  NI_scatau, NI_extcoef, NI_scacoef,  NI_angstrom, &
                         NI_exttau25,  NI_scatau25, NI_exttauFM,  NI_scatauFM, &
                         NI_fluxu, NI_fluxv, rc)


   return

CONTAINS

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  NI_Compute_Diags - Calculate dust 2D diagnostics
!
! !INTERFACE:
!

   subroutine NI_Compute_Diags ( i1, i2, j1, j2, km, nbins, gcNI, w_c, tmpu, rhoa, u, v, &
                                 NH3sfcmass, NH3colmass, NH3mass, NH3conc, &
                                 NH4sfcmass, NH4colmass, NH4mass, NH4conc, &
                                 sfcmass, colmass, mass, conc, &
                                 sfcmass25, colmass25, mass25, conc25, &
                                 exttau, scatau, extcoef, scacoef, angstrom, &
                                 exttau25, scatau25, exttaufm, scataufm, &
                                 fluxu, fluxv, rc )

! !USES:

  implicit NONE

! !INPUT PARAMETERS:
   integer, intent(in) :: i1, i2, j1, j2, km, nbins
   type(NI_GridComp1), intent(inout):: gcNI     ! NI Grid Component
   type(Chem_Bundle), intent(in)   :: w_c      ! Chem Bundle
   real, pointer, dimension(:,:,:) :: tmpu     ! temperature [K]
   real, pointer, dimension(:,:,:) :: rhoa     ! air density [kg m-3]
   real, pointer, dimension(:,:,:) :: u        ! east-west wind [m s-1]
   real, pointer, dimension(:,:,:) :: v        ! north-south wind [m s-1]

! !OUTPUT PARAMETERS:
   type(Chem_Array), intent(inout)  :: NH3sfcmass  ! NH3 sfc mass concentration kg/m3
   type(Chem_Array), intent(inout)  :: NH3colmass  ! NH3 col mass density kg/m2
   type(Chem_Array), intent(inout)  :: NH3mass     ! NH3 3d mass mixing ratio kg/kg
   type(Chem_Array), intent(inout)  :: NH3conc     ! NH3 3d mass concentration, kg/m3
   type(Chem_Array), intent(inout)  :: NH4sfcmass  ! NH4 sfc mass concentration kg/m3
   type(Chem_Array), intent(inout)  :: NH4colmass  ! NH4 col mass density kg/m2
   type(Chem_Array), intent(inout)  :: NH4mass     ! NH4 3d mass mixing ratio kg/kg
   type(Chem_Array), intent(inout)  :: NH4conc     ! NH4 3d mass concentration, kg/m3
   type(Chem_Array), intent(inout)  :: sfcmass     ! nitrate sfc mass concentration kg/m3
   type(Chem_Array), intent(inout)  :: colmass     ! nitrate col mass density kg/m2
   type(Chem_Array), intent(inout)  :: conc        ! nitrate 3d mass concentration, kg/m3
   type(Chem_Array), intent(inout)  :: mass        ! 3d mass mixing ratio kg/kg
   type(Chem_Array), intent(inout)  :: sfcmass25   ! nitrate sfc mass concentration kg/m3 [pm2.5]
   type(Chem_Array), intent(inout)  :: colmass25   ! nitrate col mass density kg/m2 [pm2.5]
   type(Chem_Array), intent(inout)  :: conc25      ! nitrate 3d mass concentration, kg/m3 [pm2.5]
   type(Chem_Array), intent(inout)  :: mass25      ! 3d mass mixing ratio kg/kg [pm2.5]
   type(Chem_Array), intent(inout)  :: exttau25    ! ext. AOT at 550 nm [pm2.5]
   type(Chem_Array), intent(inout)  :: scatau25    ! sct. AOT at 550 nm [pm2.5]
   type(Chem_Array), intent(inout)  :: exttaufm    ! ext. AOT at 550 nm [pm1.0]
   type(Chem_Array), intent(inout)  :: scataufm    ! sct. AOT at 550 nm [pm1.0]
   type(Chem_Array), intent(inout)  :: exttau      ! ext. AOT at 550 nm
   type(Chem_Array), intent(inout)  :: scatau      ! sct. AOT at 550 nm
   type(Chem_Array), intent(inout)  :: extcoef     ! 3d ext. coefficient, 1/m
   type(Chem_Array), intent(inout)  :: scacoef     ! 3d scat.coefficient, 1/m
   type(Chem_Array), intent(inout)  :: angstrom    ! 470-870 nm Angstrom parameter
   type(Chem_Array), intent(inout)  :: fluxu       ! Column mass flux in x direction
   type(Chem_Array), intent(inout)  :: fluxv       ! Column mass flux in y direction
   integer, intent(out)             :: rc          ! Error return code:
                                                   !  0 - all is well
                                                   !  1 - 

! !DESCRIPTION: Calculates some simple 2d diagnostics from the NI fields
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
   character(len=*), parameter :: myname = 'NI_Compute_Diags'
   integer :: i, j, k, n, n1, n2, ios, nch, idx
   integer :: nNH3, nNH4a, nNO3an1, nNO3an2, nNO3an3
   real :: tau, ssa
!  For now we do not try to explicitly resolve the PM fractions per bin;
!  This could be implemented as in dust if we provide bin edges
!  For now we simply use the first size bin of nitrate as the PM1 and PM2.5
!  component
!   real :: fPMfm(nbins)  ! fraction of bin with particles diameter < 1.0 um
!   real :: fPM25(nbins)  ! fraction of bin with particles diameter < 2.5 um
   character(len=255) :: qname
   real, dimension(i1:i2,j1:j2) :: tau470, tau870
   real    :: ilam550, ilam470, ilam870
   logical :: do_angstrom


!  Initialize local variables
!  --------------------------
   n1  = w_c%reg%i_NI
   n2  = w_c%reg%j_NI

   nNH3    = n1 + globalnNH3    - 1
   nNH4a   = n1 + globalnNH4a   - 1
   nNO3an1 = n1 + globalnNO3an1 - 1
   nNO3an2 = n1 + globalnNO3an2 - 1
   nNO3an3 = n1 + globalnNO3an3 - 1

   nch   = gcNI%mie_tables%nch

!  Get the wavelength indices
!  --------------------------
!  Must provide ilam550 for AOT calculation
   ilam550 = 1.
   ilam470 = 0.
   ilam870 = 0.
   if(nch .gt. 1) then
    do i = 1, nch
     if ( gcNI%mie_tables%channels(i) .ge. 5.49e-7 .and. &
          gcNI%mie_tables%channels(i) .le. 5.51e-7) ilam550 = i
     if ( gcNI%mie_tables%channels(i) .ge. 4.69e-7 .and. &
          gcNI%mie_tables%channels(i) .le. 4.71e-7) ilam470 = i
     if ( gcNI%mie_tables%channels(i) .ge. 8.69e-7 .and. &
          gcNI%mie_tables%channels(i) .le. 8.71e-7) ilam870 = i
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

!  NH3 diagnostics
!  ---------------
!  Calculate the NH3 mass mixing ratio
   if( associated(NH3mass%data3d) ) then
      NH3mass%data3d = w_c%qa(nNH3)%data3d
   endif
!  Calculate the NH3 mass concentration
   if( associated(NH3conc%data3d) ) then
      NH3conc%data3d = w_c%qa(nNH3)%data3d*rhoa
   endif
!  Calculate the NH3 surface mass concentration
   if( associated(NH3sfcmass%data2d) ) then
      NH3sfcmass%data2d(i1:i2,j1:j2) = w_c%qa(nNH3)%data3d(i1:i2,j1:j2,km)*rhoa(i1:i2,j1:j2,km)
   endif
!  Calculate the NH3 column loading
   if( associated(NH3colmass%data2d) ) then
      NH3colmass%data2d(i1:i2,j1:j2) = 0.
      do k = 1, km
        NH3colmass%data2d(i1:i2,j1:j2) &
         =   NH3colmass%data2d(i1:i2,j1:j2) + w_c%qa(nNH3)%data3d(i1:i2,j1:j2,k)*w_c%delp(i1:i2,j1:j2,k)/grav
      end do
   endif

!  NH4 diagnostics
!  ---------------
!  Calculate the NH4 mass mixing ratio
   if( associated(NH4mass%data3d) ) then
      NH4mass%data3d = w_c%qa(nNH4a)%data3d
   endif
!  Calculate the NH4 mass concentration
   if( associated(NH4conc%data3d) ) then
      NH4conc%data3d = w_c%qa(nNH4a)%data3d*rhoa
   endif
!  Calculate the NH4 surface mass concentration
   if( associated(NH4sfcmass%data2d) ) then
      NH4sfcmass%data2d(i1:i2,j1:j2) = w_c%qa(nNH4a)%data3d(i1:i2,j1:j2,km)*rhoa(i1:i2,j1:j2,km)
   endif
!  Calculate the NH4 column loading
   if( associated(NH4colmass%data2d) ) then
      NH4colmass%data2d(i1:i2,j1:j2) = 0.
      do k = 1, km
        NH4colmass%data2d(i1:i2,j1:j2) &
         =   NH4colmass%data2d(i1:i2,j1:j2) + w_c%qa(nNH4a)%data3d(i1:i2,j1:j2,k)*w_c%delp(i1:i2,j1:j2,k)/grav
      end do
   endif

!  Nitrate mass diagnostics
!  -----------------------------------------------
!  Calculate the nitrate surface mass concentration
   if( associated(sfcmass%data2d) ) then
      sfcmass%data2d(i1:i2,j1:j2) = 0.
      do n = globalnNO3an1, globalnNO3an3
         sfcmass%data2d(i1:i2,j1:j2) &
              =   sfcmass%data2d(i1:i2,j1:j2) &
              + w_c%qa(n1+n-1)%data3d(i1:i2,j1:j2,km)*rhoa(i1:i2,j1:j2,km)
      end do
   endif

!  Calculate the nitrate column loading
   if( associated(colmass%data2d) ) then
      colmass%data2d(i1:i2,j1:j2) = 0.
      do n = globalnNO3an1, globalnNO3an3
       do k = 1, km
        colmass%data2d(i1:i2,j1:j2) &
         =   colmass%data2d(i1:i2,j1:j2) &
           + w_c%qa(n1+n-1)%data3d(i1:i2,j1:j2,k)*w_c%delp(i1:i2,j1:j2,k)/grav
       end do
      end do
   endif

!  Calculate the nitrate total mass concentration
   if( associated(conc%data3d) ) then
      conc%data3d(i1:i2,j1:j2,1:km) = 0.
      do n = globalnNO3an1, globalnNO3an3
       conc%data3d(i1:i2,j1:j2,1:km) &
         =   conc%data3d(i1:i2,j1:j2,1:km) &
           + w_c%qa(n1+n-1)%data3d(i1:i2,j1:j2,1:km)*rhoa(i1:i2,j1:j2,1:km)
      end do
   endif

!  Calculate the nitrate total mass mixing ratio
   if( associated(mass%data3d) ) then
      mass%data3d(i1:i2,j1:j2,1:km) = 0.
      do n = globalnNO3an1, globalnNO3an3
       mass%data3d(i1:i2,j1:j2,1:km) &
         =   mass%data3d(i1:i2,j1:j2,1:km) &
           + w_c%qa(n1+n-1)%data3d(i1:i2,j1:j2,1:km)
      end do
   endif

!  Calculate the PM2.5 diagnostics
   n = nNO3an1
   if( associated(sfcmass25%data2d) ) &
       sfcmass25%data2d(i1:i2,j1:j2) = w_c%qa(n)%data3d(i1:i2,j1:j2,km)*rhoa(i1:i2,j1:j2,km)
   if( associated(colmass25%data2d) ) then
      colmass25%data2d(i1:i2,j1:j2) = 0.
       do k = 1, km
        colmass25%data2d(i1:i2,j1:j2) &
         =   colmass25%data2d(i1:i2,j1:j2) &
           + w_c%qa(n)%data3d(i1:i2,j1:j2,k)*w_c%delp(i1:i2,j1:j2,k)/grav
       end do
   endif
   if( associated(conc25%data3d) ) &
       conc25%data3d(i1:i2,j1:j2,1:km) = w_c%qa(n)%data3d(i1:i2,j1:j2,1:km)*rhoa(i1:i2,j1:j2,1:km)
   if( associated(mass25%data3d) ) &
       mass25%data3d(i1:i2,j1:j2,1:km) = w_c%qa(n)%data3d(i1:i2,j1:j2,1:km)

!  Calculate the nitrate column mass flux in x direction
   if( associated(fluxu%data2d) ) then
      fluxu%data2d(i1:i2,j1:j2) = 0.
      do n = globalnNO3an1, globalnNO3an3
       do k = 1, km
        fluxu%data2d(i1:i2,j1:j2) &
         =   fluxu%data2d(i1:i2,j1:j2) &
           + w_c%qa(n1+n-1)%data3d(i1:i2,j1:j2,k)*w_c%delp(i1:i2,j1:j2,k)/grav*u(i1:i2,j1:j2,k)
       end do
      end do
   endif   
   
!  Calculate the nitrate column mass flux in y direction
   if( associated(fluxv%data2d) ) then
      fluxv%data2d(i1:i2,j1:j2) = 0.
      do n = globalnNO3an1, globalnNO3an3
       do k = 1, km
        fluxv%data2d(i1:i2,j1:j2) &
         =   fluxv%data2d(i1:i2,j1:j2) &
           + w_c%qa(n1+n-1)%data3d(i1:i2,j1:j2,k)*w_c%delp(i1:i2,j1:j2,k)/grav*v(i1:i2,j1:j2,k)
       end do
      end do
   endif      

!  Calculate the nitrate optical quantities
!  ----------------------------------------
   if( associated(exttau%data2d) .or. associated(scatau%data2d) ) then

      if( associated(exttau%data2d) )   exttau%data2d(i1:i2,j1:j2) = 0.
      if( associated(scatau%data2d) )   scatau%data2d(i1:i2,j1:j2) = 0.

      if( associated(exttau25%data2d) ) exttau25%data2d(i1:i2,j1:j2) = 0.
      if( associated(scatau25%data2d) ) scatau25%data2d(i1:i2,j1:j2) = 0.

      if( associated(exttaufm%data2d) ) exttaufm%data2d(i1:i2,j1:j2) = 0.
      if( associated(scataufm%data2d) ) scataufm%data2d(i1:i2,j1:j2) = 0.

      if( associated(extcoef%data3d))   extcoef%data3d(i1:i2,j1:j2,1:km) = 0.
      if( associated(scacoef%data3d))   scacoef%data3d(i1:i2,j1:j2,1:km) = 0.

      do n = globalnNO3an1, globalnNO3an3

!      Select the name for species
       qname = trim(w_c%reg%vname(n+n1-1))
       idx = Chem_MieQueryIdx(gcNI%mie_tables,qname,rc)
       if(rc .ne. 0) call die(myname, 'cannot find proper Mie table index')

       do k = 1, km
        do j = j1, j2
         do i = i1, i2
          call Chem_MieQuery(gcNI%mie_tables, idx, ilam550, &
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
          if( associated(exttau%data2d) ) exttau%data2d(i,j) = exttau%data2d(i,j) + tau
          if( associated(scatau%data2d) ) scatau%data2d(i,j) = scatau%data2d(i,j) + tau*ssa
          if( n .eq. globalnNO3an1) then
           if( associated(exttau25%data2d) ) exttau25%data2d(i,j) = exttau25%data2d(i,j) + tau
           if( associated(scatau25%data2d) ) scatau25%data2d(i,j) = scatau25%data2d(i,j) + tau*ssa
           if( associated(exttaufm%data2d) ) exttaufm%data2d(i,j) = exttaufm%data2d(i,j) + tau
           if( associated(scataufm%data2d) ) scataufm%data2d(i,j) = scataufm%data2d(i,j) + tau*ssa
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

      do n = globalnNO3an1, globalnNO3an3

!      Select the name for species
       qname = trim(w_c%reg%vname(n+n1-1))
       idx = Chem_MieQueryIdx(gcNI%mie_tables,qname,rc)
       if(rc .ne. 0) call die(myname, 'cannot find proper Mie table index')

       do k = 1, km
        do j = j1, j2
         do i = i1, i2

          call Chem_MieQuery(gcNI%mie_tables, idx, ilam470, &
              w_c%qa(n+n1-1)%data3d(i,j,k)*w_c%delp(i,j,k)/grav, &
              w_c%rh(i,j,k), tau=tau)
          tau470(i,j) = tau470(i,j) + tau

          call Chem_MieQuery(gcNI%mie_tables, idx, ilam870, &
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

   end subroutine NI_Compute_Diags

 end subroutine NI_GridCompRun2_



!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  NI_GridCompFinalize --- The Chem Driver 
!
! !INTERFACE:
!

   subroutine NI_GridCompFinalize1_ ( gcNI, w_c, impChem, expChem, &
                                      nymd, nhms, cdt, rc )

! !USES:

  implicit NONE

! !INPUT/OUTPUT PARAMETERS:

   type(NI_GridComp1), intent(inout) :: gcNI   ! Grid Component

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

   character(len=*), parameter :: myname = 'NI_GridCompFinalize'
   rc=0
   return

 end subroutine NI_GridCompFinalize1_

 end module NI_GridCompMod


!-----------------------------------------------------------------------

!                     Single Instance Wrapper

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  NI_SingleInstance_ --- Runs single instance of method
!
! !INTERFACE:
!
  subroutine NI_SingleInstance_ ( Method_, instance, &
                                  gcNI, w_c, impChem, expChem, &
                                  nymd, nhms, cdt, rc )

! !USES:

  Use NI_GridCompMod
  Use ESMF
  Use MAPL
  Use Chem_Mod 

  IMPLICIT NONE

! !INPUT PARAMETERS:

!  Input "function pointer"
!  -----------------------
   interface 
     subroutine Method_ (gc, w, imp, exp, ymd, hms, dt, rcode )
       Use NI_GridCompMod
       Use ESMF
       Use MAPL
       Use Chem_Mod 
       type(NI_GridComp1),  intent(inout)  :: gc
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

   TYPE(NI_GridComp1), INTENT(INOUT) :: gcNI    ! Grid Component
   TYPE(ESMF_State), INTENT(INOUT)  :: impChem  ! Import State
   TYPE(ESMF_State), INTENT(INOUT)  :: expChem  ! Export State
   INTEGER, INTENT(OUT) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 - 

! !DESCRIPTION: Finalizes the NI Grid Component. Multiple instance
!               version.
!
! !REVISION HISTORY:
!
!  27Feb2008  da Silva  Introduced multiple instances
!
!EOP
!-------------------------------------------------------------------------

  integer n_NI, i_NI, j_NI
  character(len=255) :: nh3_qname, nh4_qname, no3an1_qname, no3an2_qname, no3an3_qname

! Save overall NI indices
! -----------------------
  n_NI = w_c%reg%n_NI
  i_NI = w_c%reg%i_NI
  j_NI = w_c%reg%j_NI

! Save the name of the variables in this instance
! -----------------------------------------------
  nh3_qname    = trim(w_c%reg%vname(i_NI + 5*(instance - 1)))
  nh4_qname    = trim(w_c%reg%vname(i_NI + 5*(instance - 1)+1))
  no3an1_qname = trim(w_c%reg%vname(i_NI + 5*(instance - 1)+2))
  no3an2_qname = trim(w_c%reg%vname(i_NI + 5*(instance - 1)+3))
  no3an3_qname = trim(w_c%reg%vname(i_NI + 5*(instance - 1)+4))
  
! Customize indices for this particular instance
! ----------------------------------------------
  w_c%reg%n_NI = 5
  w_c%reg%i_NI = i_NI + 5*(instance - 1)
  w_c%reg%j_NI = i_NI + 5*(instance - 1) + 4
  w_c%reg%vname(i_NI + 5*(instance - 1))     = w_c%reg%vname(i_NI)
  w_c%reg%vname(i_NI + 5*(instance - 1)+1)   = w_c%reg%vname(i_NI+1)
  w_c%reg%vname(i_NI + 5*(instance - 1)+2)   = w_c%reg%vname(i_NI+2)
  w_c%reg%vname(i_NI + 5*(instance - 1)+3)   = w_c%reg%vname(i_NI+3)
  w_c%reg%vname(i_NI + 5*(instance - 1)+4)   = w_c%reg%vname(i_NI+4)
  
! Execute the instance method
! ---------------------------
  call Method_ ( gcNI, w_c, impChem, expChem, &
                 nymd, nhms, cdt, rc )

! Restore the overall NI indices
! ------------------------------
  w_c%reg%vname(i_NI + 5*(instance - 1))     = nh3_qname
  w_c%reg%vname(i_NI + 5*(instance - 1)+1)   = nh4_qname
  w_c%reg%vname(i_NI + 5*(instance - 1)+2)   = no3an1_qname
  w_c%reg%vname(i_NI + 5*(instance - 1)+3)   = no3an2_qname
  w_c%reg%vname(i_NI + 5*(instance - 1)+4)   = no3an3_qname
  w_c%reg%n_NI = n_NI
  w_c%reg%i_NI = i_NI
  w_c%reg%j_NI = j_NI

  end subroutine NI_SingleInstance_

!-----------------------------------------------------------------------
