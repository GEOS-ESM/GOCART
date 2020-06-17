#include "MAPL_Generic.h"

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  SS_GridCompMod --- SS Grid Component Class
!
! !INTERFACE:
!

   module  SS_GridCompMod

! !USES:

   USE ESMF
   USE MAPL

   use Chem_Mod              ! Chemistry Base Class
   use Chem_StateMod         ! Chemistry State
   use Chem_ConstMod, only: grav  ! Constants 
   use Chem_UtilMod          ! I/O
   use Chem_MieMod           ! Aerosol LU Tables, calculator
   use m_inpak90             ! Resource file management
   use m_die, only: die
   use SeasaltEmissionMod    ! Emissions
   use Chem_SettlingMod      ! Settling
   use DryDepositionMod      ! Dry deposition
   use WetRemovalMod         ! Large scale wet removal
   use ConvectionMod         ! Offline convective mixing/scavenging

   implicit none

! !PUBLIC TYPES:
!
   PRIVATE
   PUBLIC  SS_GridComp       ! The SS object
   PUBLIC  SS_GridComp1      ! Single instance SS object



!
! !PUBLIC MEMBER FUNCTIONS:
!

   PUBLIC  SS_GridCompSetServices
   PUBLIC  SS_GridCompInitialize
   PUBLIC  SS_GridCompRun1
   PUBLIC  SS_GridCompRun2
   PUBLIC  SS_GridCompFinalize

!
! !DESCRIPTION:
!
!  This module implements the (pre-ESMF) SS Grid Component. 
!
! !REVISION HISTORY:
!
!  16Sep2003 da Silva  First crack.
!
!EOP
!-------------------------------------------------------------------------

  type SS_GridComp1
        character(len=255) :: name
        character(len=255) :: iname           ! instance name
        character(len=255) :: rcfilen         ! resource file name

        integer :: instance                   ! instance number

        logical :: run_alarm = .false.        ! run alarm

        type(Chem_Mie), pointer :: mie_tables => null() ! aod LUTs
        integer       :: rhFlag                ! Choice of relative humidity parameterization for radius
        integer       :: sstemisFlag           ! Choice of SST correction to emissions: 0 - none; 1 - Jaegle et al. 2011; 2 - GEOS5
        logical       :: hoppelFlag            ! Apply the Hoppel correction to emissions (Fan and Toon, 2011)
        logical       :: weibullFlag           ! Apply the Weibull distribution to wind speed for emissions (Fan and Toon, 2011)
        integer       :: emission_scheme       ! Emission scheme to use (see SeasaltEmissionMod)
        real          :: emission_scale        ! Global tuning factor for emissions
        real, pointer :: radius(:) => null()   ! particle effective radius [um]
        real, pointer :: rLow(:) => null()     ! lower radius of particle bin [um]
        real, pointer :: rUp(:) => null()      ! upper radius of particle bin [um]
        real, pointer :: rhop(:) => null()     ! dry salt particle density [kg m-3]
        real, pointer :: deep_lakes_mask(:,:) => null() 
                                               ! mask used to supress emissions from lakes that are OCEAN tiles
  end type SS_GridComp1

 type SS_GridComp
     integer                     :: n = 0                ! number of instances
     type(Chem_Mie), pointer     :: mie_tables => null() ! aod LUTs
     type(SS_GridComp1), pointer :: gcs(:)     => null() ! instances
  end type SS_GridComp

  character(len=*), parameter :: rc_basename = 'SS_GridComp'



  real, parameter :: OCEAN=0.0, LAND = 1.0, SEA_ICE = 2.0

  real, parameter :: radToDeg = 57.2957795

CONTAINS

   subroutine SS_GridCompSetServices(  gc, chemReg, rc)
   type(ESMF_GridComp), intent(INOUT) :: GC
   type(Chem_Registry), intent(INOUT) :: chemReg
   integer,             intent(OUT  ) :: rc

   CHARACTER(LEN=255) :: name

   integer            :: status
   character(len=ESMF_MAXSTR) :: Iam

   integer :: n,i

   type(ESMF_Config) :: cfg

   Iam = 'SS_GridCompSetServices'

!  Load resource file
!  ------------------
   cfg = ESMF_ConfigCreate(rc=status)
   VERIFY_(STATUS)
   call ESMF_ConfigLoadFile(cfg,trim(rc_basename)//'.rc',rc=status)
   VERIFY_(STATUS)

!  Parse resource file
!  -------------------
   n = ESMF_ConfigGetLen(cfg,label='SS_instances:',rc=status)
   VERIFY_(STATUS)


!  We have 5 tracers for each instance of SS
!  We cannot have fewer instances than half the number of
!  SS bins in the registry (it is OK to have less, though)
!  --------------------------------------------------------
   if ( 5*n .LT. chemReg%n_SS ) then
        rc = 35
        return
   else if ( 5*n .GT. chemReg%n_SS ) then
        if (MAPL_AM_I_ROOT()) &
        print *, trim(Iam)// &
                 ': fewer SS bins than possible SS instances: ',&
                 n, chemReg%n_SS/5
   end if
   n = min(n,chemReg%n_SS/5 )

!  Record name of each instance
!  ----------------------------
   call ESMF_ConfigFindLabel(cfg,'SS_instances:',rc=status)
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
      call SS_GridCompSetServices1_(gc,chemReg,name,rc=status)
      VERIFY_(STATUS)
   end do

   RETURN_(ESMF_SUCCESS)
   end subroutine SS_GridCompSetServices


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  SS_GridCompInitialize --- Initialize SS_GridComp
!
! !INTERFACE:
!

   subroutine SS_GridCompInitialize ( gcSS, w_c, impChem, expChem, &
                                      nymd, nhms, cdt, rc )

! !USES:

   implicit NONE

! !INPUT PARAMETERS:

   type(Chem_Bundle), intent(inout) :: w_c     ! Chemical tracer fields      
   integer, intent(in) :: nymd, nhms           ! time
   real, intent(in) :: cdt                     ! chemistry timestep (secs)

! !OUTPUT PARAMETERS:

   type(SS_GridComp), intent(inout) :: gcSS     ! Grid Component
   type(ESMF_State), intent(inout)  :: impChem  ! Import State
   type(ESMF_State), intent(inout)  :: expChem  ! Export State
   integer, intent(out) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 - 

! !DESCRIPTION: Initializes the SS Grid Component. It primarily sets
!               the import state for each active constituent package.
!
! !REVISION HISTORY:
!
!  18Sep2003 da Silva  First crack.
!
!EOP
!-------------------------------------------------------------------------

   character(len=*), parameter :: myname = 'SS_GridCompInitialize'
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
   CALL I90_label ( 'SS_instances:', ier )
   if ( ier .NE. 0 ) then
      rc = 20
      return
   end if

!  First determine how many instances we have
!  ------------------------------------------   
   n = 0
   do while ( ier .EQ. 0 )
      CALL I90_gtoken( name, ier )
      n = n + 1
   end do
   if ( n .EQ. 0 ) then
      rc = 30
      return
   end if
   
!  We have 5 tracers for each instance of SS
!  We cannot have fewer instances than half the number of
!  SS bins in the registry (it is OK to have less, though)
!  --------------------------------------------------------
   if ( 5*n .LT. w_c%reg%n_SS ) then
        rc = 35
        return
   else if ( 5*n .GT. w_c%reg%n_SS ) then
        if (MAPL_AM_I_ROOT()) &
        print *, trim(myname)// &
                 ': fewer SS bin sets than possible SS instances'//&
                 ' (5 bins per instance): ',&
                 n, w_c%reg%n_SS
   end if
   n = min(n,w_c%reg%n_SS/5)
   gcSS%n = n

!  Next allocate necessary memory
!  ------------------------------
   allocate ( gcSS%gcs(n), stat=ier )
   if ( ier .NE. 0 ) then
      rc = 40
      return
   end if

!  Record name of each instance
!  ----------------------------
   CALL I90_label ( 'SS_instances:', ier )
   do i = 1, n
      CALL I90_gtoken( name, ier )
      if ( ier .NE. 0 ) then
         rc = 40
         return
      end if
                                            ! resource file name
      gcSS%gcs(i)%rcfilen = trim(rc_basename)//'---'//trim(name)//'.rc'
      gcSS%gcs(i)%instance = i              ! instance number 
      IF(TRIM(name) == "full" ) THEN
       gcSS%gcs(i)%iname = " "              ! blank instance name for full (1)
      ELSE
       gcSS%gcs(i)%iname = TRIM(name)       ! instance name for others
      END IF
   end do    

!  Next initialize each instance
!  -----------------------------
   do i = 1, gcSS%n
      IF(MAPL_AM_I_ROOT()) THEN
       PRINT *," "
       PRINT *,myname,": Initializing instance ",TRIM(gcSS%gcs(i)%iname)," [",gcSS%gcs(i)%instance,"]"
      END IF
      call SS_SingleInstance_ ( SS_GridCompInitialize1_, i, &
                                gcSS%gcs(i), w_c, impChem, expChem,  &
                                nymd, nhms, cdt, ier )
      if ( ier .NE. 0 ) then
         rc = 1000+ier
         return
      end if
      gcSS%gcs(i)%mie_tables => gcSS%mie_tables
   end do

!  All done
!  --------
   CALL I90_FullRelease( ier )
   IF( ier /= 0 ) THEN
    PRINT *,myname,": I90_FullRelease not successful."
    rc = 40
   END IF

   end subroutine SS_GridCompInitialize


   
   
!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  SS_GridCompRun1 --- Run SS_GridComp
!
! !INTERFACE:
!

   subroutine SS_GridCompRun1 ( gcSS, w_c, impChem, expChem, &
                                      nymd, nhms, cdt, rc )

! !USES:

   IMPLICIT NONE

! !INPUT PARAMETERS:

   TYPE(Chem_Bundle), intent(in) :: w_c        ! Chemical tracer fields      
   INTEGER, INTENT(IN) :: nymd, nhms           ! time
   REAL,    INTENT(IN) :: cdt                  ! chemical timestep (secs)


! !OUTPUT PARAMETERS:

   TYPE(SS_GridComp), INTENT(INOUT) :: gcSS     ! Grid Component
   TYPE(ESMF_State), INTENT(INOUT)  :: impChem  ! Import State
   TYPE(ESMF_State), INTENT(INOUT)  :: expChem  ! Export State
   INTEGER, INTENT(OUT) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 - 

! !DESCRIPTION: Runs the SS Grid Component. Multiple instance
!               version.
!
! !REVISION HISTORY:
!
!  27Feb2008  da Silva  Introduced multiple instances
!
!EOP
!-------------------------------------------------------------------------

   integer :: i, ier

   do i = 1, gcSS%n
      call SS_SingleInstance_ ( SS_GridCompRun1_, i, &
                                gcSS%gcs(i), w_c, impChem, expChem, &
                                nymd, nhms, cdt, ier )
      if ( ier .NE. 0 ) then
         rc = i * 1000+ier
         return
      end if
   end do

   end subroutine SS_GridCompRun1

 
   
   
!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  SS_GridCompRun2 --- Run SS_GridComp
!
! !INTERFACE:
!

   subroutine SS_GridCompRun2 ( gcSS, w_c, impChem, expChem, &
                                      run_alarm, nymd, nhms, cdt, rc )

! !USES:

   IMPLICIT NONE

! !INPUT PARAMETERS:

   TYPE(Chem_Bundle), intent(in) :: w_c        ! Chemical tracer fields      
   LOGICAL, INTENT(IN) :: run_alarm            ! run alarm
   INTEGER, INTENT(IN) :: nymd, nhms           ! time
   REAL,    INTENT(IN) :: cdt                  ! chemical timestep (secs)


! !OUTPUT PARAMETERS:

   TYPE(SS_GridComp), INTENT(INOUT) :: gcSS     ! Grid Component
   TYPE(ESMF_State), INTENT(INOUT)  :: impChem  ! Import State
   TYPE(ESMF_State), INTENT(INOUT)  :: expChem  ! Export State
   INTEGER, INTENT(OUT) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 - 

! !DESCRIPTION: Runs the SS Grid Component. Multiple instance
!               version.
!
! !REVISION HISTORY:
!
!  27Feb2008  da Silva  Introduced multiple instances
!
!EOP
!-------------------------------------------------------------------------

   integer :: i, ier

   do i = 1, gcSS%n
      gcSS%gcs(i)%run_alarm = run_alarm

      call SS_SingleInstance_ ( SS_GridCompRun2_, i, &
                                gcSS%gcs(i), w_c, impChem, expChem, &
                                nymd, nhms, cdt, ier )
      if ( ier .NE. 0 ) then
         rc = i * 1000+ier
         return
      end if
   end do

   end subroutine SS_GridCompRun2


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  SS_GridCompFinalize --- Finalize SS_GridComp
!
! !INTERFACE:
!

   subroutine SS_GridCompFinalize ( gcSS, w_c, impChem, expChem, &
                                      nymd, nhms, cdt, rc )

! !USES:

  IMPLICIT NONE

! !INPUT PARAMETERS:

   TYPE(Chem_Bundle), intent(in) :: w_c        ! Chemical tracer fields      
   INTEGER, INTENT(IN) :: nymd, nhms           ! time
   REAL,    INTENT(IN) :: cdt                  ! chemical timestep (secs)


! !OUTPUT PARAMETERS:

   TYPE(SS_GridComp), INTENT(INOUT) :: gcSS     ! Grid Component
   TYPE(ESMF_State), INTENT(INOUT)  :: impChem  ! Import State
   TYPE(ESMF_State), INTENT(INOUT)  :: expChem  ! Export State
   INTEGER, INTENT(OUT) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 - 

! !DESCRIPTION: Finalizes the SS Grid Component. Multiple instance
!               version.
!
! !REVISION HISTORY:
!
!  27Feb2008  da Silva  Introduced multiple instances
!
!EOP
!-------------------------------------------------------------------------

   integer :: i, ier

   do i = 1, gcSS%n
      call SS_SingleInstance_ ( SS_GridCompFinalize1_, i, &
                                gcSS%gcs(i), w_c, impChem, expChem, &
                                nymd, nhms, cdt, ier )
      if ( ier .NE. 0 ) then
         rc = i * 1000+ier
         return
      end if
   end do

   if (associated(gcSS%gcs)) deallocate ( gcSS%gcs, stat=ier )
   gcSS%n = -1

   end subroutine SS_GridCompFinalize


   subroutine SS_GridCompSetServices1_(  gc, chemReg, iname, rc)
   type(ESMF_GridComp), intent(INOUT) :: GC
   type(Chem_Registry), intent(INOUT) :: chemReg
   character(len=*),    intent(IN   ) :: iname
   integer,             intent(OUT  ) :: rc

   character(len=ESMF_MAXSTR) :: Iam

   Iam = "SS_GridCompSetServices1_"

   ! Import spec goes here... 

   RETURN_(ESMF_SUCCESS)

   end subroutine SS_GridCompSetServices1_


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  SS_GridCompInitialize --- Initialize SS_GridComp
!
! !INTERFACE:
!

   subroutine SS_GridCompInitialize1_ ( gcSS, w_c, impChem, expChem, &
                                      nymd, nhms, cdt, rc )

! !USES:

  implicit NONE

! !INPUT PARAMETERS:

   type(Chem_Bundle), intent(inout) :: w_c     ! Chemical tracer fields      
   integer, intent(in) :: nymd, nhms           ! time
   real,    intent(in) :: cdt                  ! chemistry timestep (secs)

! !OUTPUT PARAMETERS:

   type(SS_GridComp1), intent(inout) :: gcSS    ! Grid Component
   type(ESMF_State), intent(inout)   :: impChem ! Import State
   type(ESMF_State), intent(inout)   :: expChem ! Export State
   integer, intent(out) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 - 

! !DESCRIPTION: Initializes the SS Grid Component. It primarily sets
!               the import state for each active constituent package.
!
! !REVISION HISTORY:
!
!  18Sep2003 da Silva  First crack.
!
!EOP
!-------------------------------------------------------------------------

   character(len=*), parameter :: myname = 'SS_GridCompInitialize1_'


   character(len=255) :: rcfilen
   integer :: n
   integer, allocatable :: ier(:)
   integer :: i1, i2, im, j1, j2, jm, nbins, n1, n2, nbins_rc, i, j
   real :: qmin, qmax, dummylon
   real :: radius, rlow, rup, rhop, fscav, fnum, molwght, rnum
   integer :: irhFlag, isstemisFlag, ihoppelFlag, iweibullFlag, iemission_scheme

   integer, parameter :: nhres = 6   ! number of horizontal model resolutions: a,b,c,d,e
   real    :: emission_scale(nhres)  ! scale factor buffer



   rcfilen = trim(gcSS%rcfilen)
   gcSS%name = 'SS Constituent Package'

!  Initialize local variables
!  --------------------------
   rc = 0
   i1 = w_c%grid%i1; i2 = w_c%grid%i2; im = w_c%grid%im
   j1 = w_c%grid%j1; j2 = w_c%grid%j2; jm = w_c%grid%jm
   nbins = w_c%reg%n_SS
   n1  = w_c%reg%i_SS
   n2  = w_c%reg%j_SS

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

   call i90_label ( 'number_SS_bins:', ier(1) )
   nbins_rc = i90_gint ( ier(2) )
   if ( any(ier(1:2) /= 0) ) then
      call final_(20)
      return
   end if
   if ( nbins_rc /= nbins ) then
      call final_(25)
      return
   end if

!  Particle radius
!  ---------------
   call i90_label ( 'particle_radius:', ier(1) )
   do n = 1, nbins
      radius           = i90_gfloat ( ier(n+1) )
      gcSS%radius(n)   = radius
   end do
   if ( any(ier(1:nbins+1) /= 0) ) then
      call final_(50)
      return
   end if

!  Particle radius (lower bound)
!  ---------------
   call i90_label ( 'radius_lower:', ier(1) )
   do n = 1, nbins
      rlow                  = i90_gfloat ( ier(n+1) )
      gcSS%rlow(n)          = rlow
      w_c%reg%rlow(n1+n-1)  = rlow * 1.e-6
   end do
   if ( any(ier(1:nbins+1) /= 0) ) then
      call final_(50)
      return
   end if

!  Particle radius (upper bound)
!  ---------------
   call i90_label ( 'radius_upper:', ier(1) )
   do n = 1, nbins
      rup                 = i90_gfloat ( ier(n+1) )
      gcSS%rup(n)         = rup
      w_c%reg%rup(n1+n-1) = rup * 1.e-6
   end do
   if ( any(ier(1:nbins+1) /= 0) ) then
      call final_(50)
      return
   end if
!                          -------

!  Dry Particle Density
!  ---------------
   call i90_label ( 'SS_density:', ier(1) )
   do n = 1, nbins
      rhop                 = i90_gfloat ( ier(n+1) )
      gcSS%rhop(n)         = rhop
      w_c%reg%rhop(n1+n-1) = rhop
   end do
   if ( any(ier(1:nbins+1) /= 0) ) then
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
      fscav                   = i90_gfloat ( ier(n+1) )
      w_c%reg%fscav(n1+n-1)   = fscav
      w_c%qa(n1+n-1)%fscav    = fscav
   end do
   if ( any(ier(1:nbins+1) /= 0) ) then
      call final_(50)
      return
   end if
!                          -------

!  Number median radius
!  To be used in droplet activation code
!  ---------------
   call i90_label ( 'particle_radius_number:', ier(1) )
   do n = 1, nbins
      rnum                    = i90_gfloat ( ier(n+1) )
      w_c%reg%rmed(n1+n-1)    = rnum * 1e-6
   end do
   if ( any(ier(1:nbins+1) /= 0) ) then
      call final_(50)
      return
   end if
!                          -------

!  Number to mass conversion factor
!  To be used in droplet activation code
!  ---------------
   call i90_label ( 'fnum:', ier(1) )
   do n = 1, nbins
      fnum                    = i90_gfloat ( ier(n+1) )
      w_c%reg%fnum(n1+n-1)    = fnum
   end do
   if ( any(ier(1:nbins+1) /= 0) ) then
      call final_(50)
      return
   end if
!                          -------

!  Molecular weight
!  To be used in droplet activation code
!  ---------------
   call i90_label ( 'molecular_weight:', ier(1) )
   do n = 1, nbins
      molwght                 = i90_gfloat ( ier(n+1) )
      w_c%reg%molwght(n1+n-1) = molwght
   end do
   if ( any(ier(1:nbins+1) /= 0) ) then
      call final_(50)
      return
   end if
!                          -------

!  Particle affected by relative humidity?
!  ---------------
   call i90_label ( 'rhFlag:', ier(1) )
   irhFlag                    = i90_gint ( ier(2) )
   gcSS%rhFlag                = irhFlag
   if ( any(ier(1:2) /= 0) ) then
      call final_(50)
      return
   end if

!  Which Emissions Scheme to Use (see SeasaltEmissionMod)
!  ---------------
   call i90_label ( 'emission_scheme:', ier(1) )
   iemission_scheme           = i90_gint ( ier(2) )
   gcSS%emission_scheme       = iemission_scheme
   if ( any(ier(1:2) /= 0) ) then
      call final_(50)
      return
   end if

!  Emissions Efficiency
!  Scaling factor to multiply calculated
!  emissions by.  Applies to all size bins.
!  ---------------
   CALL I90_Label ( 'emission_scale:', ier(1) )
   do n = 1, nhres
      emission_scale(n) = i90_gfloat ( ier(n+1) )
   end do
   gcSS%emission_scale = Chem_UtilResVal(im, jm, emission_scale(:), ier(nhres + 2))
   if ( any(ier(1:nhres+2) /= 0) ) then
      call final_(50)
      return
   end if
!                          -------

!  SST correction to emission strength following Jaegle et al, 2011
!  ---------------
   call i90_label ( 'sstemisFlag:', ier(1) )
   isstemisFlag = i90_gint ( ier(2) )
   if ((isstemisFlag < 0) .or. (isstemisFlag > 2)) then
      gcSS%sstemisFlag = 0     ! unknown correction - fall back to no SST correction
   else
      gcSS%sstemisFlag = isstemisFlag
   end if
   if ( any(ier(1:2) /= 0) ) then
      call final_(50)
      return
   end if

!  Hoppel 2005 correction to emissions (Fan and Toon 2011)
!  ---------------
   call i90_label ( 'hoppelFlag:', ier(1) )
   ihoppelFlag = i90_gint ( ier(2) )
   if (ihoppelFlag /= 0) then
      gcSS%hoppelFlag = .True.
   else
      gcSS%hoppelFlag = .False.
   end if
   if ( any(ier(1:2) /= 0) ) then
      call final_(50)
      return
   end if

!  Weibull wind speed adjustment to emissions (Fan and Toon 2011)
!  ---------------
   call i90_label ( 'weibullFlag:', ier(1) )
   iweibullFlag = i90_gint ( ier(2) )
   if (iweibullFlag /= 0) then
      gcSS%weibullFlag = .True.
   else
      gcSS%weibullFlag = .False.
   end if
   if ( any(ier(1:2) /= 0) ) then
      call final_(50)
      return
   end if

!  Mask to prevent emissions from the Great Lakes and the Caspian Sea
!  ------------------------------------------------------------------
   gcSS%deep_lakes_mask = 1.0

   do j = j1, j2
     do i = i1, i2
                             dummylon = w_c%grid%lon(i,j)*radToDeg
        if( dummylon < 0.0 ) dummylon = dummylon + 360.0
       ! The Great Lakes: lon = [91W,75W], lat = [40.5N, 50N]
       if ((dummylon > 267.0) .and. &
           (dummylon < 285.0) .and. &
           (w_c%grid%lat(i,j)*radToDeg >  40.5) .and. &
           (w_c%grid%lat(i,j)*radToDeg <  50.0)) gcSS%deep_lakes_mask(i,j) = 0.0 

       ! The Caspian Sea: lon = [45.0, 56], lat = 35, 48]
       if ((dummylon >  45.0) .and. &
           (dummylon <  56.0) .and. &
           (w_c%grid%lat(i,j)*radToDeg >  35.0) .and. &
           (w_c%grid%lat(i,j)*radToDeg <  48.0)) gcSS%deep_lakes_mask(i,j) = 0.0
     end do
   end do



!  All done
!  --------
   call i90_release()
   deallocate(ier)

   return

CONTAINS

   subroutine init_()
   integer ios, nerr
   nerr = max ( 32, nbins+1 )
   allocate ( gcSS%radius(nbins), gcSS%rLow(nbins), gcSS%rUp(nbins), &
              gcSS%rhop(nbins), gcSS%deep_lakes_mask(i1:i2,j1:j2), ier(nerr), stat=ios )
   if ( ios /= 0 ) rc = 100
   end subroutine init_

   subroutine final_(ierr)
   integer :: ierr
   integer ios
   deallocate ( gcSS%radius, gcSS%rhop, gcSS%rLow, gcSS%rUp, gcSS%deep_lakes_mask, ier, stat=ios )
   call i90_release()
   rc = ierr
   end subroutine final_

   end subroutine SS_GridCompInitialize1_




!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  SS_GridCompRun1_ --- The Chem Driver, run phase 1
!
! !INTERFACE:
!

   subroutine SS_GridCompRun1_ ( gcSS, w_c, impChem, expChem, &
                               nymd, nhms, cdt, rc )

! !USES:

  implicit NONE

! !INPUT/OUTPUT PARAMETERS:

   type(SS_GridComp1), intent(inout) :: gcSS    ! Grid Component
   type(Chem_Bundle), intent(inout)  :: w_c     ! Chemical tracer fields   

! !INPUT PARAMETERS:

   type(ESMF_State), intent(inout) :: impChem   ! Import State
   integer, intent(in) :: nymd, nhms            ! time
   real,    intent(in) :: cdt                   ! chemistry timestep (secs)

! !OUTPUT PARAMETERS:

   type(ESMF_State), intent(inout) :: expChem   ! Export State
   integer, intent(out) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 -
 
! !DESCRIPTION: This routine implements the so-called SS Driver. That 
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

   character(len=*), parameter :: Iam = 'SS_GridCompRun1_'

   integer :: i1, i2, im, j1, j2, jm, nbins, n1, n2, km, n, ijl, ijkl, ijk1l, i, j
   real :: qmin, qmax
   real, pointer :: SS_radius(:), SS_rhop(:)
   real, pointer :: memissions(:,:), nemissions(:,:), w10m(:,:), dqa(:,:)

!  Input fields from fvGCM
!  -----------------------
   real, pointer, dimension(:,:)     :: frlake, frocean, frseaice, oro
   real, pointer, dimension(:,:)     :: u10m, v10m, ustar
   real, pointer, dimension(:,:)     :: dz
   real, pointer, dimension(:,:)     :: tskin
   real, pointer, dimension(:,:,:)   :: tmpu, rhoa
   real, pointer, dimension(:,:)     :: area

!  Modifications to source function
   real, allocatable, dimension(:,:) :: tskin_c
   real, allocatable, dimension(:,:) :: fsstemis
   real, allocatable, dimension(:,:) :: fgridefficiency
   real, allocatable, dimension(:,:) :: fhoppel, vsettle
   real                              :: radius_wet, rhop_wet, diff_coef
   double precision                  :: a, c, k, wt, x
   double precision, allocatable, dimension(:,:) :: gweibull, wm

#define EXPORT        expChem
#define iNAME         TRIM(gcSS%iname)

#define ptrSSEM       SS_emis

   integer :: STATUS

#include "SS_GetPointer___.h"

!  Initialize local variables
!  --------------------------
   rc = 0
   i1 = w_c%grid%i1; i2 = w_c%grid%i2; im = w_c%grid%im
   j1 = w_c%grid%j1; j2 = w_c%grid%j2; jm = w_c%grid%jm

   km = w_c%grid%km
   nbins = w_c%reg%n_SS
   n1  = w_c%reg%i_SS
   n2  = w_c%reg%j_SS

   ijl  = ( i2 - i1 + 1 ) * ( j2 - j1 + 1 )
   ijkl = ijl * km
   ijk1l = ijl * (km+1)


!  Seasalt particle radius [m] and density [kg m-3]
!  ---------------------------------------------
   allocate(SS_radius(nbins), SS_rhop(nbins), __STAT__)

   SS_radius = 1.e-6*gcSS%radius
   SS_rhop   = gcSS%rhop

   allocate(w10m(i1:i2,j1:j2), dqa(i1:i2,j1:j2), __STAT__)
   allocate(memissions(i1:i2,j1:j2), nemissions(i1:i2,j1:j2), __STAT__)


!  Get 2D Imports
!  --------------
   call MAPL_GetPointer ( impChem, frocean,  'FROCEAN',  __RC__ )
   call MAPL_GetPointer ( impChem, frseaice, 'FRACI',    __RC__ )
   call MAPL_GetPointer ( impChem, oro,      'LWI',      __RC__ )
   call MAPL_GetPointer ( impChem, u10m,     'U10M',     __RC__ )
   call MAPL_GetPointer ( impChem, v10m,     'V10M',     __RC__ )
   call MAPL_GetPointer ( impChem, ustar,    'USTAR',    __RC__ )
   call MAPL_GetPointer ( impChem, tskin,    'TS',       __RC__ )
   call MAPL_GetPointer ( impChem, dz,       'DZ',       __RC__ )
   call MAPL_GetPointer ( impChem, frlake,   'FRLAKE',   __RC__ )
   call MAPL_GetPointer ( impChem, area,     'AREA',     __RC__ )

!  Define 10-m wind speed
   w10m = sqrt(u10m*u10m + v10m*v10m)

!  Get 3D Imports
!  --------------
   call MAPL_GetPointer ( impChem, tmpu,     'T',        __RC__ )
   call MAPL_GetPointer ( impChem, rhoa,     'AIRDENS',  __RC__ )

    

#ifdef DEBUG

   call pmaxmin('SS: frocean    ', frocean,  qmin, qmax, ijl,1, 1. )
   call pmaxmin('SS: frseaice   ', frseaice, qmin, qmax, ijl,1, 1. )
   call pmaxmin('SS: oro        ', oro     , qmin, qmax, ijl,1, 1. )
   call pmaxmin('SS: u10m       ', u10m    , qmin, qmax, ijl,1, 1. )
   call pmaxmin('SS: v10m       ', v10m    , qmin, qmax, ijl,1, 1. )
   call pmaxmin('SS: ustar      ', ustar   , qmin, qmax, ijl,1, 1. )
   call pmaxmin('SS: tskin      ', tskin   , qmin, qmax, ijl,1, 1. )
   call pmaxmin('SS: rhoa       ', rhoa    , qmin, qmax, ijkl,1, 1. )

#endif

!  Seasalt Source (and modifications)
!  -----------
!  Grid box efficiency to emission (fraction of sea water)
   allocate(fgridefficiency(i1:i2,j1:j2), __STAT__ )
   fgridefficiency = min(max(0.,(frocean-frseaice)*gcSS%deep_lakes_mask),1.)


!  Apply SST correction to emissions
   allocate(fsstemis(i1:i2,j1:j2), __STAT__ )

   fsstemis = 1.0

   if (gcSS%sstemisFlag == 1) then          ! SST correction folowing Jaegle et al. 2011
       fsstemis = 0.0

       allocate( tskin_c(i1:i2,j1:j2), __STAT__ )
       tskin_c  = tskin - 273.15
       fsstemis = (0.3 + 0.1*tskin_c - 0.0076*tskin_c**2 + 0.00021*tskin_c**3)

       where(fsstemis < 0.0) fsstemis = 0.0

       deallocate( tskin_c, __STAT__ )
   else if (gcSS%sstemisFlag == 2) then     ! GEOS5 SST correction
       fsstemis = 0.0

       allocate( tskin_c(i1:i2,j1:j2), __STAT__ )
       tskin_c  = tskin - 273.15
    
       where(tskin_c < -0.1) tskin_c = -0.1    ! temperature range (0, 36) C 
       where(tskin_c > 36.0) tskin_c = 36.0    !

       fsstemis = (-1.107211 -0.010681*tskin_c -0.002276*tskin_c**2 + 60.288927*1.0/(40.0 - tskin_c))
       where(fsstemis < 0.0) fsstemis = 0.0
       where(fsstemis > 7.0) fsstemis = 7.0

       deallocate( tskin_c, __STAT__ )
   end if

!  Apply a Weibull distribution to emissions wind speeds
!  The Weibull distribution correction ends up being a multiplicative constant (g) times
!  our present source function (see Eq. 12 in Fan & Toon, 2011 and notes for 9/22/11).
!  This constant is derived from the incomplete and complete forms of the gamma
!  function, hence the utilities pasted below.  The Weibull function and shape
!  parameters (k, c) assumed are from Justus 1978.

   allocate(gweibull(i1:i2,j1:j2), wm(i1:i2,j1:j2), __STAT__ )
   gweibull = 1.0
   wm = sqrt(u10m**2 + v10m**2)   ! mean wind speed
   wt = 4.d0                      ! a threshold (Fan & Toon, 2011)

   if (gcSS%weibullFlag) then
       gweibull = 0.0

       do j = j1, j2
       do i = i1, i2
           if (wm(i,j) > 0.01) then
               k = 0.94d0 * sqrt(wm(i,j))         ! Weibull shape parameter
               c = wm(i,j) / gamma(1.d0 + 1.d0/k) ! Weibull shape parameter
               x = (wt / c) ** k
               a = 3.41d0 / k + 1.d0
               gweibull(i,j)  = (c / wm(i,j))**3.41d0 * igamma(a,x)
           end if
       end do ! i
       end do ! j
   endif

!  Loop over bins and do emission calculation
!  Possibly apply the Hoppel correction based on fall speed (Fan and Toon, 2011)
   allocate(fhoppel(i1:i2,j1:j2), vsettle(i1:i2,j1:j2), __STAT__ )
   fhoppel = 1.0
   do n = 1, nbins
    memissions = 0.
    nemissions = 0.
    dqa = 0.
    call SeasaltEmission( gcSS%rLow(n), gcSS%rUp(n), gcSS%emission_scheme, w10m, ustar, &
                          memissions, nemissions, rc )
!   For the Hoppel correction need to compute the wet radius and settling velocity
!   in the surface
    if(gcSS%hoppelFlag) then
     do j = j1, j2
      do i = i1, i2
       call wet_radius ( SS_radius(n), SS_rhop(n), w_c%rh(i,j,km), gcSS%rhFlag, &
                         radius_wet, rhop_wet )
       call Chem_CalcVsettle ( radius_wet, rhop_wet, rhoa(i,j,km), tmpu(i,j,km), &
                               diff_coef, vsettle(i,j) )
       fhoppel(i,j) = (10./dz(i,j)) ** (vsettle(i,j)/MAPL_KARMAN/ustar(i,j))
      end do
     end do
    endif

!   For the moment, do not apply these corrections to the emissions
    memissions = gcSS%emission_scale * fgridefficiency * fsstemis * fhoppel * gweibull * memissions

    dqa = memissions * cdt * grav / w_c%delp(:,:,km)

    w_c%qa(n1+n-1)%data3d(:,:,km) = w_c%qa(n1+n-1)%data3d(:,:,km) + dqa

    if (associated(SS_emis(n)%data2d)) then
        SS_emis(n)%data2d = memissions
    end if 
   end do

#ifdef DEBUG
   do n = n1, n2
      call pmaxmin('SS: q_emi', w_c%qa(n)%data3d(i1:i2,j1:j2,1:km), &
                                qmin, qmax, ijl, km, 1. )
   end do
#endif


!  Clean up
!  --------
   deallocate(SS_radius, SS_rhop, __STAT__)
   deallocate(memissions, nemissions, __STAT__)
   deallocate(w10m, dqa, __STAT__)
   deallocate(fsstemis, fgridefficiency, __STAT__)
   deallocate(fhoppel, vsettle, gweibull, wm, __STAT__)

   return

 end subroutine SS_GridCompRun1_




!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  SS_GridCompRun2 --- The Chem Driver 
!
! !INTERFACE:
!

   subroutine SS_GridCompRun2_ ( gcSS, w_c, impChem, expChem, &
                                 nymd, nhms, cdt, rc )

! !USES:

  implicit NONE

! !INPUT/OUTPUT PARAMETERS:

   type(SS_GridComp1), intent(inout) :: gcSS    ! Grid Component
   type(Chem_Bundle), intent(inout)  :: w_c     ! Chemical tracer fields   

! !INPUT PARAMETERS:

   type(ESMF_State), intent(inout) :: impChem   ! Import State
   integer, intent(in) :: nymd, nhms            ! time
   real,    intent(in) :: cdt                   ! chemistry timestep (secs)

! !OUTPUT PARAMETERS:

   type(ESMF_State), intent(inout) :: expChem   ! Export State
   integer, intent(out) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 -
 
! !DESCRIPTION: This routine implements the so-called SS Driver. That 
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

   character(len=*), parameter :: Iam = 'SS_GridCompRun2_'

   integer :: i1, i2, im, j1, j2, jm, nbins, n1, n2, km, n, ijl, ijkl, ijk1l, k
   real :: qmin, qmax
   real, pointer :: SS_radius(:), SS_rhop(:)
   real, pointer :: dqa(:,:), drydepositionfrequency(:,:)
   type(Chem_Array), pointer :: fluxout
   logical :: KIN

!  Input fields from fvGCM
!  -----------------------
   real, pointer, dimension(:,:)   ::  frlake, frocean, frseaice, &
                                       oro, u10m, v10m, &
                                       ustar, precc, precl, pblh, dz,      &
                                       shflux, z0h, hsurf, tskin
   real, pointer, dimension(:,:,:) ::  tmpu, rhoa, u, v, hghte, ple
   real, pointer, dimension(:,:,:) ::  pfllsan, pfilsan

!  Additional needs for GOCART convective diagnostic - Run2
   real, pointer, dimension(:,:,:)       ::  cmfmc, qlcn, qicn, dtrain
   real, pointer, dimension(:,:)         ::  area
   real*8, allocatable, dimension(:,:,:) ::  cmfmc_, qccu_, dtrain_, &
                                             airmass_, airmol_, vud_, &
                                             delz_, delp_, tmpu_, ple_
   real*8, allocatable                   ::  tc_(:,:,:,:), bcnv_(:,:,:)
   real*8, allocatable                   ::  area_(:,:), frlake_(:,:), &
                                             frocean_(:,:), frseaice_(:,:)
   integer*4                             ::  icdt


#define EXPORT        expChem
#define iNAME         TRIM(gcSS%iname)

#define ptrSSWT       SS_wet
#define ptrSSSV       SS_conv
#define ptrSSEM       SS_emis
#define ptrSSDP       SS_dep
#define ptrSSSD       SS_set

#define ptrSSMASS     SS_mass
#define ptrSSMASS25   SS_mass25
#define ptrSSSMASS    SS_sfcmass
#define ptrSSCMASS    SS_colmass
#define ptrSSEXTTAU   SS_exttau
#define ptrSSSCATAU   SS_scatau
#define ptrSSSMASS25  SS_sfcmass25
#define ptrSSCMASS25  SS_colmass25
#define ptrSSEXTT25   SS_exttau25
#define ptrSSSCAT25   SS_scatau25
#define ptrSSAERIDX   SS_aeridx
#define ptrSSCONC     SS_conc
#define ptrSSEXTCOEF  SS_extcoef
#define ptrSSSCACOEF  SS_scacoef
#define ptrSSEXTTFM   SS_exttaufm
#define ptrSSSCATFM   SS_scataufm
#define ptrSSANGSTR   SS_angstrom
#define ptrSSFLUXU    SS_fluxu
#define ptrSSFLUXV    SS_fluxv

   integer :: STATUS

#include "SS_GetPointer___.h"

!  Initialize local variables
!  --------------------------
   rc = 0
   i1 = w_c%grid%i1; i2 = w_c%grid%i2; im = w_c%grid%im
   j1 = w_c%grid%j1; j2 = w_c%grid%j2; jm = w_c%grid%jm

   km = w_c%grid%km
   nbins = w_c%reg%n_SS
   n1  = w_c%reg%i_SS
   n2  = w_c%reg%j_SS

   ijl  = ( i2 - i1 + 1 ) * ( j2 - j1 + 1 )
   ijkl = ijl * km
   ijk1l = ijl * (km+1)


!  Get 2D Imports
!  --------------
   call MAPL_GetPointer ( impChem, frocean,  'FROCEAN',  __RC__ )
   call MAPL_GetPointer ( impChem, frseaice, 'FRACI',    __RC__ )
   call MAPL_GetPointer ( impChem, oro,      'LWI',      __RC__ )
   call MAPL_GetPointer ( impChem, u10m,     'U10M',     __RC__ )
   call MAPL_GetPointer ( impChem, v10m,     'V10M',     __RC__ )
   call MAPL_GetPointer ( impChem, ustar,    'USTAR',    __RC__ )
   call MAPL_GetPointer ( impChem, precc,    'CN_PRCP',  __RC__ )
   call MAPL_GetPointer ( impChem, precl,    'NCN_PRCP', __RC__ )
   call MAPL_GetPointer ( impChem, tskin,    'TS',       __RC__ )
   call MAPL_GetPointer ( impChem, z0h,      'Z0H',      __RC__ )
   call MAPL_GetPointer ( impChem, pblh,     'ZPBL',     __RC__ )
   call MAPL_GetPointer ( impChem, shflux,   'SH',       __RC__ )
   call MAPL_GetPointer ( impChem, dz,       'DZ',       __RC__ )
   call MAPL_GetPointer ( impChem, frlake,   'FRLAKE',   __RC__ )
   call MAPL_GetPointer ( impChem, area,     'AREA',     __RC__ )


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
   hsurf => hghte(i1:i2,j1:j2,km) ! in GEOS-5 hghte is in [0,km]
    

#ifdef DEBUG

   call pmaxmin('SS: frocean    ', frocean,  qmin, qmax, ijl,1, 1. )
   call pmaxmin('SS: frseaice   ', frseaice, qmin, qmax, ijl,1, 1. )
   call pmaxmin('SS: oro        ', oro     , qmin, qmax, ijl,1, 1. )
   call pmaxmin('SS: u10m       ', u10m    , qmin, qmax, ijl,1, 1. )
   call pmaxmin('SS: v10m       ', v10m    , qmin, qmax, ijl,1, 1. )
   call pmaxmin('SS: ustar      ', ustar   , qmin, qmax, ijl,1, 1. )
   call pmaxmin('SS: precc      ', precc   , qmin, qmax, ijl,1, 1. )
   call pmaxmin('SS: precl      ', precl   , qmin, qmax, ijl,1, 1. )
   call pmaxmin('SS: pblh       ', pblh    , qmin, qmax, ijl,1, 1. )
   call pmaxmin('SS: shfflux    ', shflux  , qmin, qmax, ijl,1, 1. )
   call pmaxmin('SS: z0h        ', z0h     , qmin, qmax, ijl,1, 1. )
   call pmaxmin('SS: hsurf      ', hsurf   , qmin, qmax, ijl,1, 1. )
   call pmaxmin('SS: tskin      ', tskin   , qmin, qmax, ijl,1, 1. )
   call pmaxmin('SS: tmpu       ', tmpu    , qmin, qmax, ijkl,1, 1. )
   call pmaxmin('SS: rhoa       ', rhoa    , qmin, qmax, ijkl,1, 1. )
   call pmaxmin('SS: u          ', u       , qmin, qmax, ijkl,1, 1. )
   call pmaxmin('SS: v          ', v       , qmin, qmax, ijkl,1, 1. )
   call pmaxmin('SS: hghte      ', hghte   , qmin, qmax, ijkl,1, 1. )
   call pmaxmin('SS: pfllsan    ', pfllsan , qmin, qmax, ijk1l,1, 1. )
   call pmaxmin('SS: pfilsan    ', pfilsan , qmin, qmax, ijk1l,1, 1. )

#endif


RUN_ALARM: if (gcSS%run_alarm) then

!  Seasalt particle radius [m] and density [kg m-3]
!  ---------------------------------------------
   allocate( SS_radius(nbins), SS_rhop(nbins), __STAT__ )
   SS_radius = 1.e-6*gcSS%radius
   SS_rhop   = gcSS%rhop

   allocate( fluxout, __STAT__)
   allocate( fluxout%data2d(i1:i2,j1:j2), dqa(i1:i2,j1:j2), &
             drydepositionfrequency(i1:i2,j1:j2), __STAT__)


!  Seasalt Settling
!  ----------------
   call Chem_Settling ( i1, i2, j1, j2, km, n1, n2, nbins, gcSS%rhFlag, &
                        SS_radius, SS_rhop, cdt, w_c, tmpu, rhoa, hsurf,    &
                        hghte, SS_set, rc )

#ifdef DEBUG
   do n = n1, n2
      call pmaxmin('SS: q_set', w_c%qa(n)%data3d(i1:i2,j1:j2,1:km), qmin, qmax, &
                    ijl, km, 1. )
   end do
#endif

!  Seasalt Deposition
!  -----------
   drydepositionfrequency = 0.
   call DryDepositionGOCART( i1, i2, j1, j2, km, &
                             tmpu, rhoa, hghte, oro, ustar, &
                             pblh, shflux, z0h, drydepositionfrequency, rc )

   ! increase deposition velocity over land
   where (abs(oro - LAND) < 0.5) 
       drydepositionfrequency = 5.0 * drydepositionfrequency
   end where

   do n = 1, nbins
    dqa = 0.
    dqa = max(0.0, w_c%qa(n1+n-1)%data3d(:,:,km)*(1.-exp(-drydepositionfrequency*cdt)))
    w_c%qa(n1+n-1)%data3d(:,:,km) = &
            w_c%qa(n1+n-1)%data3d(:,:,km) - dqa
    if( associated(SS_dep(n)%data2d) ) &
     SS_dep(n)%data2d = dqa*w_c%delp(:,:,km)/grav/cdt
   end do

#ifdef DEBUG
   do n = n1, n2
      call pmaxmin('SS: q_dry', w_c%qa(n)%data3d(i1:i2,j1:j2,1:km), qmin, qmax, &
                    ijl, km, 1. )
   end do
#endif

!  Seasalt Large-scale Wet Removal
!  -------------------------------
   KIN = .TRUE.
   do n = 1, nbins
    w_c%qa(n1+n-1)%fwet = 1.
    call WetRemovalGOCART(i1, i2, j1, j2, km, n1+n-1, n1+n-1, cdt, 'sea_salt', KIN, &
                          w_c%qa, ple, tmpu, rhoa, pfllsan, pfilsan, &
                          precc, precl, fluxout, rc )
    if(associated(SS_wet(n)%data2d)) SS_wet(n)%data2d = fluxout%data2d
   end do

#ifdef DEBUG
   do n = n1, n2
      call pmaxmin('SS: q_wet', w_c%qa(n)%data3d(i1:i2,j1:j2,1:km), qmin, qmax, &
                    ijl, km, 1. )
   end do
#endif


!  Seasalt Convective-scale Mixing and Wet Removal
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
   call convection(i1, i2, j1, j2, km, n1, n2, icdt, 'sea_salt', kin, &
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
   do n = 1, nbins
    if(associated(SS_conv(n)%data2d)) SS_conv(n)%data2d = -bcnv_(:,:,n1+n-1)/area_/icdt
   end do

!  Clean up
!  --------
   
   deallocate(cmfmc_, qccu_, dtrain_, tc_, airmass_, &
              delz_, vud_, delp_, airmol_, tmpu_, bcnv_, ple_, &
              area_, frlake_, frocean_, frseaice_, __STAT__ )

   deallocate(fluxout%data2d, __STAT__)
   deallocate(fluxout, __STAT__)

   deallocate(dqa, drydepositionfrequency, __STAT__ ) 

   deallocate(SS_radius, SS_rhop, __STAT__)

   end if RUN_ALARM

!  Compute the desired output diagnostics here
!  Ideally this will go where chemout is called in fvgcm.F since that
!  will reflect the distributions after transport, etc.
!  ------------------------------------------------------------------
   call SS_Compute_Diags(i1, i2, j1, j2, km, nbins, gcSS, w_c, tmpu, rhoa, u, v, &
                         SS_sfcmass, SS_colmass, SS_mass, SS_exttau, SS_scatau, &
                         SS_sfcmass25, SS_colmass25, SS_mass25, SS_exttau25, SS_scatau25, &
                         SS_conc, SS_extcoef, SS_scacoef, SS_exttaufm, SS_scataufm, &
                         SS_angstrom, SS_fluxu, SS_fluxv, rc)

   return

CONTAINS
!##############################################################################
!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  SS_Compute_Diags - Calculate seasalt 2D diagnostics
!
! !INTERFACE:
!

   subroutine SS_Compute_Diags ( i1, i2, j1, j2, km, nbins, gcSS, w_c, tmpu, rhoa, u, v, &
                                 sfcmass, colmass, mass, exttau, scatau, &
                                 sfcmass25, colmass25, mass25, exttau25, scatau25, &
                                 conc, extcoef, scacoef, exttaufm, scataufm, &
                                 angstrom, fluxu, fluxv, rc )

! !USES:

  implicit NONE

! !INPUT PARAMETERS:
   integer, intent(in) :: i1, i2, j1, j2, km, nbins
   type(SS_GridComp1), intent(inout) :: gcSS   ! SS Grid Component
   type(Chem_Bundle), intent(in)   :: w_c      ! Chem Bundle
   real, pointer, dimension(:,:,:) :: tmpu     ! temperature [K]
   real, pointer, dimension(:,:,:) :: rhoa     ! air density [kg m-3]
   real, pointer, dimension(:,:,:) :: u        ! east-west wind [m s-1]
   real, pointer, dimension(:,:,:) :: v        ! north-south wind [m s-1]

! !OUTPUT PARAMETERS:
   type(Chem_Array), intent(inout)  :: sfcmass   ! sfc mass concentration kg/m3
   type(Chem_Array), intent(inout)  :: colmass   ! col mass density kg/m2
   type(Chem_Array), intent(inout)  :: mass      ! 3d mass mixing ratio kg/kg
   type(Chem_Array), intent(inout)  :: exttau    ! ext. AOT at 550 nm
   type(Chem_Array), intent(inout)  :: scatau    ! sct. AOT at 550 nm
   type(Chem_Array), intent(inout)  :: sfcmass25 ! sfc mass concentration kg/m3 (pm2.5)
   type(Chem_Array), intent(inout)  :: colmass25 ! col mass density kg/m2 (pm2.5)
   type(Chem_Array), intent(inout)  :: mass25    ! 3d mass mixing ratio kg/kg (pm2.5)
   type(Chem_Array), intent(inout)  :: exttau25  ! ext. AOT at 550 nm (pm2.5)
   type(Chem_Array), intent(inout)  :: scatau25  ! sct. AOT at 550 nm (pm2.5)
   type(Chem_Array), intent(inout)  :: conc      ! 3d mass concentration, kg/m3
   type(Chem_Array), intent(inout)  :: extcoef   ! 3d ext. coefficient, 1/m
   type(Chem_Array), intent(inout)  :: scacoef   ! 3d scat.coefficient, 1/m
   type(Chem_Array), intent(inout)  :: exttaufm  ! fine mode (sub-micron) ext. AOT at 550 nm
   type(Chem_Array), intent(inout)  :: scataufm  ! fine mode (sub-micron) sct. AOT at 550 nm
   type(Chem_Array), intent(inout)  :: angstrom  ! 470-870 nm Angstrom parameter
   type(Chem_Array), intent(inout)  :: fluxu     ! Column mass flux in x direction
   type(Chem_Array), intent(inout)  :: fluxv     ! Column mass flux in y direction
   integer, intent(out)             :: rc        ! Error return code:
                                                 !  0 - all is well
                                                 !  1 - 

! !DESCRIPTION: Calculates some simple 2d diagnostics from the SS fields
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
   character(len=*), parameter :: myname = 'SS_Compute_Diags'
   integer :: i, j, k, n, n1, n2, nch, idx
   real :: tau, ssa
   real :: fPMfm(nbins)  ! fraction of bin with particles diameter < 1.0 um
   real :: fPM25(nbins)  ! fraction of bin with particles diameter < 2.5 um
   character(len=255) :: qname
   real, dimension(i1:i2,j1:j2) :: tau470, tau870
   real    :: ilam550, ilam470, ilam870
   logical :: do_angstrom


!  Initialize local variables
!  --------------------------
   n1  = w_c%reg%i_SS
   n2  = w_c%reg%j_SS
   nch   = gcSS%mie_tables%nch

!  Get the wavelength indices
!  --------------------------
!  Must provide ilam550 for AOT calculation
   ilam550 = 1.
   ilam470 = 0.
   ilam870 = 0.
   if(nch .gt. 1) then
    do i = 1, nch
     if ( gcSS%mie_tables%channels(i) .ge. 5.49e-7 .and. &
          gcSS%mie_tables%channels(i) .le. 5.51e-7) ilam550 = i
     if ( gcSS%mie_tables%channels(i) .ge. 4.69e-7 .and. &
          gcSS%mie_tables%channels(i) .le. 4.71e-7) ilam470 = i
     if ( gcSS%mie_tables%channels(i) .ge. 8.69e-7 .and. &
          gcSS%mie_tables%channels(i) .le. 8.71e-7) ilam870 = i
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



!  Compute the fine mode (sub-micron) and PM2.5 bin-wise fractions
!  ------------------------------------
   call SS_Binwise_PM_Fractions(fPMfm, 0.50, gcSS%rlow, gcSS%rup, nbins)   ! 2*r < 1.0 um
   call SS_Binwise_PM_Fractions(fPM25, 1.25, gcSS%rlow, gcSS%rup, nbins)   ! 2*r < 2.5 um


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
   if( associated(sfcmass25%data2d) ) then
      sfcmass25%data2d(i1:i2,j1:j2) = 0.
      do n = 1, nbins
         sfcmass25%data2d(i1:i2,j1:j2) &
              =   sfcmass25%data2d(i1:i2,j1:j2) &
              + w_c%qa(n+n1-1)%data3d(i1:i2,j1:j2,km)*rhoa(i1:i2,j1:j2,km)*fPM25(n)
      end do
   endif

!  Calculate the seasalt column loading
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
   if( associated(colmass25%data2d)) then
      colmass25%data2d(i1:i2,j1:j2) = 0.
      do n = 1, nbins
       do k = 1, km
        colmass25%data2d(i1:i2,j1:j2) &
         =   colmass25%data2d(i1:i2,j1:j2) &
           + w_c%qa(n+n1-1)%data3d(i1:i2,j1:j2,k)*w_c%delp(i1:i2,j1:j2,k)/grav*fPM25(n)
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
   if( associated(mass25%data3d) ) then
      mass25%data3d(i1:i2,j1:j2,1:km) = 0.
      do n = 1, nbins
       mass25%data3d(i1:i2,j1:j2,1:km) &
         =   mass25%data3d(i1:i2,j1:j2,1:km) &
           + w_c%qa(n+n1-1)%data3d(i1:i2,j1:j2,1:km)*fPM25(n)
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

      if( associated(exttau%data2d)) exttau%data2d(i1:i2,j1:j2) = 0.
      if( associated(scatau%data2d)) scatau%data2d(i1:i2,j1:j2) = 0.

      if( associated(exttau25%data2d)) exttau25%data2d(i1:i2,j1:j2) = 0.
      if( associated(scatau25%data2d)) scatau25%data2d(i1:i2,j1:j2) = 0.

      if( associated(exttaufm%data2d)) exttaufm%data2d(i1:i2,j1:j2) = 0.
      if( associated(scataufm%data2d)) scataufm%data2d(i1:i2,j1:j2) = 0.

      if( associated(extcoef%data3d)) extcoef%data3d(i1:i2,j1:j2,1:km) = 0.
      if( associated(scacoef%data3d)) scacoef%data3d(i1:i2,j1:j2,1:km) = 0.

      do n = 1, nbins

!      Select the name for species
       qname = trim(w_c%reg%vname(n+n1-1))
       idx = Chem_MieQueryIdx(gcSS%mie_tables,qname,rc)
       if(rc .ne. 0) call die(myname, 'cannot find proper Mie table index')

       do k = 1, km
        do j = j1, j2
         do i = i1, i2
          call Chem_MieQuery(gcSS%mie_tables, idx, ilam550, &
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
          if( associated(exttaufm%data2d)) &
                         exttaufm%data2d(i,j) = exttaufm%data2d(i,j) + tau*fPMfm(n)          
          if( associated(exttau25%data2d)) &
                         exttau25%data2d(i,j) = exttau25%data2d(i,j) + tau*fPM25(n)

          if( associated(scatau%data2d) ) scatau%data2d(i,j) = scatau%data2d(i,j) + tau*ssa
          if( associated(scataufm%data2d) ) &
                         scataufm%data2d(i,j) = scataufm%data2d(i,j) + tau*ssa*fPMfm(n)
          if( associated(scatau25%data2d) ) &
                         scatau25%data2d(i,j) = scatau25%data2d(i,j) + tau*ssa*fPM25(n)

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
       idx = Chem_MieQueryIdx(gcSS%mie_tables,qname,rc)
       if(rc .ne. 0) call die(myname, 'cannot find proper Mie table index')

       do k = 1, km
        do j = j1, j2
         do i = i1, i2

          call Chem_MieQuery(gcSS%mie_tables, idx, ilam470, &
              w_c%qa(n+n1-1)%data3d(i,j,k)*w_c%delp(i,j,k)/grav, &
              w_c%rh(i,j,k), tau=tau)
          tau470(i,j) = tau470(i,j) + tau

          call Chem_MieQuery(gcSS%mie_tables, idx, ilam870, &
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

   end subroutine SS_Compute_Diags


!##############################################################################
!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  SS_Binwise_PM_Fractions - Calculate bin-wise PM fractions
!
! !INTERFACE:
!

   subroutine SS_Binwise_PM_Fractions(fPM, rPM, r_low, r_up, nbins)

! !USES:

  implicit NONE

! !INPUT/OUTPUT PARAMETERS:

  real, dimension(:), intent(inout) :: fPM     ! bin-wise PM fraction (r < rPM)

! !INPUT PARAMETERS:

   real,    intent(in)              :: rPM     ! PM radius
   integer, intent(in)              :: nbins   ! number of bins
   real, dimension(:), intent(in)   :: r_low   ! bin radii - low bounds
   real, dimension(:), intent(in)   :: r_up    ! bin radii - upper bounds

! !OUTPUT PARAMETERS:


! !DESCRIPTION: Calculates bin-wise PM fractions
!
! !REVISION HISTORY:
!
!EOP
!-------------------------------------------------------------------------

! !Local Variables

   integer :: n

   character(len=*), parameter :: myname = 'SS_Binwise_PM_Fractions'

   do n = 1, nbins
     if(r_up(n) < rPM) then
       fPM(n) = 1.0
     else
       if(r_low(n) < rPM) then
!        Assume dm/dlnr = constant, i.e., dm/dr ~ 1/r
         fPM(n) = log(rPM/r_low(n)) / log(r_up(n)/r_low(n))
       else
         fPM(n) = 0.0
       endif
     endif
   enddo

   end subroutine SS_Binwise_PM_Fractions

 end subroutine SS_GridCompRun2_

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  SS_GridCompFinalize --- The Chem Driver 
!
! !INTERFACE:
!

   subroutine SS_GridCompFinalize1_ ( gcSS, w_c, impChem, expChem, &
                                    nymd, nhms, cdt, rc )

! !USES:

  implicit NONE

! !INPUT/OUTPUT PARAMETERS:

   type(SS_GridComp1), intent(inout) :: gcSS  ! Grid Component

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

   character(len=*), parameter :: myname = 'SS_GridCompFinalize1_'

   if (associated(gcSS%radius))          deallocate(gcSS%radius, stat=rc)
   if (associated(gcSS%rhop))            deallocate(gcSS%rhop, stat=rc)
   if (associated(gcSS%rLow))            deallocate(gcSS%rLow, stat=rc)
   if (associated(gcSS%rUp))             deallocate(gcSS%rUp, stat=rc)
   if (associated(gcSS%deep_lakes_mask)) deallocate(gcSS%deep_lakes_mask, stat=rc)

   rc = 0

   return

 end subroutine SS_GridCompFinalize1_

 subroutine wet_radius (radius, rhop, rh, flag, &
                        radius_wet, rhop_wet)

! Compute the wet radius of sea salt particle
! Inputs:  radius - dry radius [m]
!          rhop   - dry density [kg m-3]
!          rh     - relative humidity [0-1]
!          flag   - 1 (Fitzgerald, 1975)
!                 - 2 (Gerber, 1985)
! Outputs: radius_wet - humidified radius [m]
!          rhop_wet   - wet density [kg m-3]

       real, intent(in)  :: radius, rhop, rh
       integer           :: flag                  ! control method of humidification
       real, intent(out) :: radius_wet, rhop_wet

!      Local
       real :: sat, rcm, rrat
       real, parameter ::  rhow = 1000.  ! Density of water [kg m-3]

!      The following parameters relate to the swelling of seasalt like particles
!      following Fitzgerald, Journal of Applied Meteorology, 1975.
       real, parameter :: epsilon = 1.   ! soluble fraction of deliqeuscing particle
       real, parameter :: alphaNaCl = 1.35
       real :: alpha, alpha1, alpharat, beta, theta, f1, f2

!      parameter from Gerber 1985 (units require radius in cm, see rcm)
       real, parameter :: c1=0.7674, c2=3.079, c3=2.573e-11, c4=-1.424

!      Default is to return radius as radius_wet, rhop as rhop_wet
       radius_wet = radius
       rhop_wet   = rhop

!      Make sure saturation ratio (RH) is sensible
       sat = max(rh,tiny(1.0)) ! to avoid zero FPE

!      Fitzgerald Scheme
       if(flag .eq. 1 .and. sat .ge. 0.80) then
!       parameterization blows up for RH > 0.995, so set that as max
!       rh needs to be scaled 0 - 1
        sat = min(0.995,sat)
!       Calculate the alpha and beta parameters for the wet particle
!       relative to amonium sulfate
        beta = exp( (0.00077*sat) / (1.009-sat) )
        if(sat .le. 0.97) then
         theta = 1.058
        else
         theta = 1.058 - (0.0155*(sat-0.97)) /(1.02-sat**1.4)
        endif
        alpha1 = 1.2*exp( (0.066*sat) / (theta-sat) )
        f1 = 10.2 - 23.7*sat + 14.5*sat**2.
        f2 = -6.7 + 15.5*sat - 9.2*sat**2.
        alpharat = 1. - f1*(1.-epsilon) - f2*(1.-epsilon**2.)
        alpha = alphaNaCl * (alpha1*alpharat)
!       radius_wet is the radius of the wet particle
        radius_wet = alpha * radius**beta
        rrat       = (radius/radius_wet)**3.
        rhop_wet   = rrat*rhop + (1.-rrat)*rhow
       elseif(flag .eq. 2) then   ! Gerber
        sat = min(0.995,sat)
        rcm = radius*100.
        radius_wet = 0.01 * (   c1*rcm**c2 / (c3*rcm**c4-alog10(sat)) &
                              + rcm**3.)**(1./3.)
        rrat       = (radius/radius_wet)**3.
        rhop_wet   = rrat*rhop + (1.-rrat)*rhow
       endif

 end subroutine wet_radius

!===============================================================================
!gamma and incomplete gamma functions from Tianyi Fan, but she cannot recall
!origin of the code.  I did verify against IDL for accuracy.  --prc
 double precision function gamma(X)

!----------------------------------------------------------------------- 
! Gamma function
!----------------------------------------------------------------------- 
 implicit none
double precision, intent(in)  :: X
! local variable
double precision G(26)
double precision M1, Z, M, R, GR
integer K
double precision, parameter :: PI = 4.d0 * atan(1.d0)

        IF (X.EQ.INT(X)) THEN
           IF (X.GT.0.0) THEN
              gamma=1.0d0
              M1=X-1
              DO K=2,M1
                 gamma=gamma*K
              end do
           ELSE
              gamma=1.0d+300
           ENDIF
        ELSE
           IF (ABS(X).GT.1.0) THEN
              Z=ABS(X)
              M=INT(Z)
              R=1.0
              DO  K=1,M
                 R=R*(Z-K)
              END DO
              Z=Z-M
           ELSE
              Z=X
           ENDIF

           DATA G/1.0D0,0.5772156649015329D0,                 &
               -0.6558780715202538D0, -0.420026350340952D-1,  &
               0.1665386113822915D0,-.421977345555443D-1,     &
               -.96219715278770D-2, .72189432466630D-2,       &
               -.11651675918591D-2, -.2152416741149D-3,       &
               .1280502823882D-3, -.201348547807D-4,          &
               -.12504934821D-5, .11330272320D-5,             &
               -.2056338417D-6, .61160950D-8,                 &
               .50020075D-8, -.11812746D-8,                   &
               .1043427D-9, .77823D-11,                       &
               -.36968D-11, .51D-12,                          &
               -.206D-13, -.54D-14, .14D-14, .1D-15/
           GR=G(26)
           DO K=25,1,-1
              GR=GR*Z+G(K)
           END DO
           gamma =1.0/(GR*Z)
           IF (ABS(X).GT.1.0) THEN
              gamma=gamma*R
              IF (X.LT.0.0D0) gamma =-PI/(X*gamma*SIN(PI*X))
           ENDIF
        ENDIF
        RETURN

 
 
 end function gamma
 
!===============================================================================
 DOUBLE PRECISION function igamma(A, X)
!----------------------------------------------------------------------- 
! incomplete Gamma function
!----------------------------------------------------------------------- 
 IMPLICIT NONE
 double precision, intent(in) :: 	A
 DOUBLE PRECISION, INTENT(IN) ::      X
! LOCAL VARIABLE
 DOUBLE PRECISION :: XAM, GIN,  S, R, T0
 INTEGER K
        XAM=-X+A*LOG(X)
        IF (XAM.GT.700.0.OR.A.GT.170.0) THEN
           WRITE(*,*)'IGAMMA: a and/or x too large, X = ', X
	   WRITE(*,*) 'A = ', A
           STOP
	   
        ENDIF

        IF (X.EQ.0.0) THEN           
           IGAMMA=GAMMA(A)
           
        ELSE IF (X.LE.1.0+A) THEN
           S=1.0/A
           R=S
           DO  K=1,60
              R=R*X/(A+K)
              S=S+R
              IF (ABS(R/S).LT.1.0e-15) EXIT
           END DO
           GIN=EXP(XAM)*S           
           IGAMMA=GAMMA(A)-GIN
        ELSE IF (X.GT.1.0+A) THEN
           T0=0.0
           DO K=60,1,-1
              T0=(K-A)/(1.0+K/(X+T0))
           end do

           IGAMMA=EXP(XAM)/(X+T0)

        ENDIF
 
 end function igamma

 end module SS_GridCompMod


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  SS_SingleInstance_ --- Runs single instance of method
!
! !INTERFACE:
!
  subroutine SS_SingleInstance_ ( Method_, instance, &
                                  gcSS, w_c, impChem, expChem, &
                                  nymd, nhms, cdt, rc )

! !USES:

  Use SS_GridCompMod
  Use ESMF
  Use MAPL
  Use Chem_Mod 

  IMPLICIT NONE

! !INPUT PARAMETERS:

!  Input "function pointer"
!  -----------------------
   interface 
     subroutine Method_ (gc, w, imp, exp, ymd, hms, dt, rcode )
       Use SS_GridCompMod
       Use ESMF
       Use MAPL
       Use Chem_Mod 
       type(SS_GridComp1),  intent(inout)  :: gc
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

   TYPE(SS_GridComp1), INTENT(INOUT) :: gcSS    ! Grid Component
   TYPE(ESMF_State), INTENT(INOUT)  :: impChem  ! Import State
   TYPE(ESMF_State), INTENT(INOUT)  :: expChem  ! Export State
   INTEGER, INTENT(OUT) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 - 

! !DESCRIPTION: Finalizes the SS Grid Component. Multiple instance
!               version.
!
! !REVISION HISTORY:
!
!  27Feb2008  da Silva  Introduced multiple instances
!
!EOP
!-------------------------------------------------------------------------

  integer :: n, n_SS, i_SS, j_SS
  integer :: status
  character(len=255), allocatable :: qname(:)
  character(len=ESMF_MAXSTR) :: Iam
  integer, parameter :: n_bins = 5

  Iam = 'SS_SingleInstance_'

! Save overall SS indices
! -----------------------
  n_SS = w_c%reg%n_SS
  i_SS = w_c%reg%i_SS
  j_SS = w_c%reg%j_SS

! Save the name of the variables in this instance
! -----------------------------------------------
  allocate(qname(n_bins), __STAT__)

  do n = 1, n_bins
      qname(n) = trim(w_c%reg%vname(i_SS + n_bins*(instance - 1) + n - 1))
  end do
  
! Customize indices for this particular instance
! ----------------------------------------------
  w_c%reg%n_SS = n_bins
  w_c%reg%i_SS = i_SS + n_bins*(instance - 1)
  w_c%reg%j_SS = i_SS + n_bins*(instance - 1) + (n_bins - 1)

  do n = 1, n_bins
      w_c%reg%vname(i_SS + n_bins*(instance - 1) + n - 1) = w_c%reg%vname(i_SS + n - 1)
  end do
  
! Execute the instance method
! ---------------------------
  call Method_ ( gcSS, w_c, impChem, expChem, &
                 nymd, nhms, cdt, rc )

! Restore the overall SS indices
! ------------------------------
  do n = 1, n_bins
      w_c%reg%vname(i_SS + n_bins*(instance - 1) + n - 1) = qname(n)
  end do

  w_c%reg%n_SS = n_SS
  w_c%reg%i_SS = i_SS
  w_c%reg%j_SS = j_SS

  deallocate(qname, __STAT__)

  end subroutine SS_SingleInstance_

!-----------------------------------------------------------------------

