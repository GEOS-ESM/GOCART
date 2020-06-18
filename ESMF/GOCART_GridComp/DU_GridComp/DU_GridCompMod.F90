#include "MAPL_Generic.h"

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  DU_GridCompMod --- DU Grid Component Class
!
! !INTERFACE:
!

   module  DU_GridCompMod

! !USES:

   USE ESMF
   USE MAPL

   use Chem_Mod              ! Chemistry Base Class
   use Chem_StateMod         ! Chemistry State
   use Chem_ConstMod, only: grav, undefval => undef  ! Constants
   use Chem_UtilMod          ! I/O
   use Chem_MieMod           ! Aerosol LU Tables, calculator
   use m_inpak90             ! Resource file management
   use m_die, only: die
   use DustEmissionMod       ! Emissions
   use Chem_SettlingMod      ! Settling
   use DryDepositionMod      ! Dry deposition
   use WetRemovalMod         ! Large-scale wet removal
   use ConvectionMod         ! Offline convective mixing/scavenging

   implicit none

! !PUBLIC TYPES:
!
   PRIVATE
   PUBLIC  DU_GridComp       ! The DU object
   PUBLIC  DU_GridComp1      ! Single instance DU object

!
! !PUBLIC MEMBER FUNCTIONS:
!

   PUBLIC  DU_GridCompSetServices
   PUBLIC  DU_GridCompInitialize
   PUBLIC  DU_GridCompRun1
   PUBLIC  DU_GridCompRun2
   PUBLIC  DU_GridCompFinalize

!
! !DESCRIPTION:
!
!  This module implements the (pre-ESMF) DU Grid Component. 
!
! !REVISION HISTORY:
!
!  16Sep2003 da Silva  First crack.
!  16Aug2005 da Silva  Passed ESMF grid to MPread().
!
!EOP
!-------------------------------------------------------------------------

  type DU_GridComp1
        character(len=255) :: name
        character(len=255) :: iname           ! instance name
        character(len=255) :: rcfilen         ! resource file name

        integer :: instance                   ! instance number

        logical :: run_alarm = .false.        ! run alarm
        
        type(Chem_Mie), pointer :: mie_tables => null() ! aod LUTs
        integer       :: rhFlag         ! choice of relative humidity parameterization for radius
        logical       :: maringFlag     ! settling velocity correction
        real, pointer :: src(:,:)       ! Ginoux dust sources
        real, pointer :: radius(:)      ! particle effective radius [um]
        real, pointer :: rlow(:)        ! particle effective radius lower bound [um]
        real, pointer :: rup(:)         ! particle effective radius upper bound [um]
        real, pointer :: sfrac(:)       ! fraction of total source
        real, pointer :: rhop(:)        ! soil class density [kg m-3]
        integer       :: nymd
        real          :: Ch_DU          ! dust emission tuning coefficient [kg s2 m-5].

!       Workspace for any requested point emissions (handled in run)
!       ------------------------------------------------------------
        logical :: doing_point_emissions=.FALSE.         ! providing pointwise emissions
        character(len=255) :: point_emissions_srcfilen   ! filename for pointwise emissions
        integer                         :: nPts = -1
        integer, pointer, dimension(:)  :: pstart => null(), pend => null()
        real, pointer, dimension(:)     :: pLat  => null(), &
                                           pLon  => null(), &
                                           pBase => null(), &
                                           pTop  => null(), &
                                           pEmis => null()
  end type DU_GridComp1

  type DU_GridComp
     integer                     :: n = 0                ! number of instances
     type(Chem_Mie), pointer     :: mie_tables => null() ! aod LUTs
     type(DU_GridComp1), pointer :: gcs(:)     => null() ! instances
  end type DU_GridComp

  character(len=*), parameter :: rc_basename = 'DU_GridComp'


  real, parameter :: OCEAN=0.0, LAND = 1.0, SEA_ICE = 2.0
  real, parameter :: radTODeg = 57.2957795

CONTAINS

   subroutine DU_GridCompSetServices(  gc, chemReg, rc)
   type(ESMF_GridComp), intent(INOUT) :: GC
   type(Chem_Registry), intent(INOUT) :: chemReg
   integer,             intent(OUT  ) :: rc

   CHARACTER(LEN=255) :: name

   integer            :: status
   character(len=ESMF_MAXSTR) :: Iam

   integer :: n, i

   type(ESMF_Config) :: cfg

   Iam = "DU_GridCompSetServices"

!  Load resource file
!  ------------------
   cfg = ESMF_ConfigCreate(rc=status)
   VERIFY_(STATUS)
   call ESMF_ConfigLoadFile(cfg,trim(rc_basename)//'.rc',rc=status)
   VERIFY_(STATUS)

!  Parse resource file
!  -------------------
   n = ESMF_ConfigGetLen(cfg,label='DU_instances:',rc=status)
   VERIFY_(STATUS)


!  We have 5 tracers for each instance of DU
!  We cannot have fewer instances than half the number of
!  DU bins in the registry (it is OK to have less, though)
!  --------------------------------------------------------
   if ( 5*n .LT. chemReg%n_DU ) then
        rc = 35
        return
   else if ( 5*n .GT. chemReg%n_DU ) then
        if (MAPL_AM_I_ROOT()) &
        print *, trim(Iam)// &
                 ': fewer DU bins than possible DU instances: ',&
                 n, chemReg%n_DU/5
   end if
   n = min(n,chemReg%n_DU/5 )

!  Record name of each instance
!  ----------------------------
   call ESMF_ConfigFindLabel(cfg,'DU_instances:',rc=status)
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
      call DU_GridCompSetServices1_(gc,chemReg,name,rc=status)
      VERIFY_(STATUS)
   end do

   RETURN_(ESMF_SUCCESS)
   end subroutine DU_GridCompSetServices


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  DU_GridCompInitialize --- Initialize DU_GridComp
!
! !INTERFACE:
!

   subroutine DU_GridCompInitialize ( gcDU, w_c, impChem, expChem, &
                                      nymd, nhms, cdt, rc )

! !USES:

   implicit NONE

! !INPUT PARAMETERS:

   type(Chem_Bundle), intent(inout) :: w_c     ! Chemical tracer fields      
   integer, intent(in) :: nymd, nhms           ! time
   real, intent(in) :: cdt                     ! chemistry timestep (secs)

! !OUTPUT PARAMETERS:

   type(DU_GridComp), intent(inout) :: gcDU     ! Grid Component
   type(ESMF_State), intent(inout)  :: impChem  ! Import State
   type(ESMF_State), intent(inout)  :: expChem  ! Export State
   integer, intent(out) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 - 

! !DESCRIPTION: Initializes the DU Grid Component. It primarily sets
!               the import state for each active constituent package.
!
! !REVISION HISTORY:
!
!  18Sep2003 da Silva  First crack.
!
!EOP
!-------------------------------------------------------------------------

   character(len=*), parameter :: myname = 'DU_GridCompInitialize'
   CHARACTER(LEN=255) :: name
   
   integer i, ier, n

!  Load resource file
!  ------------------
   call i90_loadf ( trim(rc_basename)//'.rc', ier )
   if ( ier .NE. 0 ) then
      rc = 10
      return
   end if

!  Parse resource file
!  -------------------
   CALL I90_label ( 'DU_instances:', ier )
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
   
!  We have 5 tracers for each instance of DU
!  We cannot have fewer instances than half the number of
!  DU bins in the registry (it is OK to have less, though)
!  --------------------------------------------------------
   if ( 5*n .LT. w_c%reg%n_DU ) then
        rc = 35
        return
   else if ( 5*n .GT. w_c%reg%n_DU ) then
        if (MAPL_AM_I_ROOT()) &
        print *, trim(myname)// &
                 ': fewer DU bin sets than possible DU instances'//&
                 ' (5 bins per instance): ',&
                 n, w_c%reg%n_DU
   end if
   n = min(n,w_c%reg%n_DU/5)
   gcDU%n = n

!  Next allocate necessary memory
!  ------------------------------
   allocate ( gcDU%gcs(n), stat=ier )
   if ( ier .NE. 0 ) then
      rc = 40
      return
   end if

!  Record name of each instance
!  ----------------------------
   CALL I90_label ( 'DU_instances:', ier )
   do i = 1, n
      CALL I90_gtoken( name, ier )
      if ( ier .NE. 0 ) then
         rc = 40
         return
      end if
                                            ! resource file name
      gcDU%gcs(i)%rcfilen = trim(rc_basename)//'---'//trim(name)//'.rc'
      gcDU%gcs(i)%instance = i              ! instance number 
      IF(TRIM(name) == "full" ) THEN
       gcDU%gcs(i)%iname = " "              ! blank instance name for full (1)
      ELSE
       gcDU%gcs(i)%iname = TRIM(name)       ! instance name for others
      END IF
   end do    

!  Next initialize each instance
!  -----------------------------
   do i = 1, gcDU%n
      IF(MAPL_AM_I_ROOT()) THEN
       PRINT *," "
       PRINT *,myname,": Initializing instance ",TRIM(gcDU%gcs(i)%iname)," [",gcDU%gcs(i)%instance,"]"
      END IF
      call DU_SingleInstance_ ( DU_GridCompInitialize1_, i, &
                                gcDU%gcs(i), w_c, impChem, expChem,  &
                                nymd, nhms, cdt, ier )
      if ( ier .NE. 0 ) then
         rc = 1000+ier
         return
      end if
      gcDU%gcs(i)%mie_tables => gcDU%mie_tables
   end do

!  All done
!  --------
   CALL I90_FullRelease( ier )
   IF( ier /= 0 ) THEN
    PRINT *,myname,": I90_FullRelease not successful."
    rc = 40
   END IF

   end subroutine DU_GridCompInitialize



!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  DU_GridCompRun1 --- Run DU_GridComp
!
! !INTERFACE:
!

   subroutine DU_GridCompRun1 ( gcDU, w_c, impChem, expChem, &
                                      nymd, nhms, cdt, rc )

! !USES:

   IMPLICIT NONE

! !INPUT PARAMETERS:

   TYPE(Chem_Bundle), intent(in) :: w_c        ! Chemical tracer fields      
   INTEGER, INTENT(IN) :: nymd, nhms           ! time
   REAL,    INTENT(IN) :: cdt                  ! chemical timestep (secs)


! !OUTPUT PARAMETERS:

   TYPE(DU_GridComp), INTENT(INOUT) :: gcDU     ! Grid Component
   TYPE(ESMF_State), INTENT(INOUT)  :: impChem  ! Import State
   TYPE(ESMF_State), INTENT(INOUT)  :: expChem  ! Export State
   INTEGER, INTENT(OUT) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 - 

! !DESCRIPTION: Runs the DU Grid Component. Multiple instance
!               version.
!
! !REVISION HISTORY:
!
!  27Feb2008  da Silva  Introduced multiple instances
!
!EOP
!-------------------------------------------------------------------------

   integer :: i, ier

   do i = 1, gcDU%n
      call DU_SingleInstance_ ( DU_GridCompRun1_, i, &
                                gcDU%gcs(i), w_c, impChem, expChem, &
                                nymd, nhms, cdt, ier )
      if ( ier .NE. 0 ) then
         rc = i * 1000+ier
         return
      end if
   end do

   end subroutine DU_GridCompRun1


   
   
!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  DU_GridCompRun2 --- Run DU_GridComp
!
! !INTERFACE:
!

   subroutine DU_GridCompRun2 ( gcDU, w_c, impChem, expChem, &
                                      run_alarm, nymd, nhms, cdt, rc )

! !USES:

   IMPLICIT NONE

! !INPUT PARAMETERS:

   TYPE(Chem_Bundle), intent(in) :: w_c        ! Chemical tracer fields      
   LOGICAL, INTENT(IN) :: run_alarm            ! run alarm
   INTEGER, INTENT(IN) :: nymd, nhms           ! time
   REAL,    INTENT(IN) :: cdt                  ! chemical timestep (secs)


! !OUTPUT PARAMETERS:

   TYPE(DU_GridComp), INTENT(INOUT) :: gcDU     ! Grid Component
   TYPE(ESMF_State), INTENT(INOUT)  :: impChem  ! Import State
   TYPE(ESMF_State), INTENT(INOUT)  :: expChem  ! Export State
   INTEGER, INTENT(OUT) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 - 

! !DESCRIPTION: Runs the DU Grid Component. Multiple instance
!               version.
!
! !REVISION HISTORY:
!
!  27Feb2008  da Silva  Introduced multiple instances
!
!EOP
!-------------------------------------------------------------------------

   integer :: i, ier

   do i = 1, gcDU%n
      gcDU%gcs(i)%run_alarm = run_alarm

      call DU_SingleInstance_ ( DU_GridCompRun2_, i, &
                                gcDU%gcs(i), w_c, impChem, expChem, &
                                nymd, nhms, cdt, ier )
      if ( ier .NE. 0 ) then
         rc = i * 1000+ier
         return
      end if
   end do

   end subroutine DU_GridCompRun2


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  DU_GridCompFinalize --- Finalize DU_GridComp
!
! !INTERFACE:
!

   subroutine DU_GridCompFinalize ( gcDU, w_c, impChem, expChem, &
                                      nymd, nhms, cdt, rc )

! !USES:

  IMPLICIT NONE

! !INPUT PARAMETERS:

   TYPE(Chem_Bundle), intent(in) :: w_c        ! Chemical tracer fields      
   INTEGER, INTENT(IN) :: nymd, nhms           ! time
   REAL,    INTENT(IN) :: cdt                  ! chemical timestep (secs)


! !OUTPUT PARAMETERS:

   TYPE(DU_GridComp), INTENT(INOUT) :: gcDU     ! Grid Component
   TYPE(ESMF_State), INTENT(INOUT)  :: impChem  ! Import State
   TYPE(ESMF_State), INTENT(INOUT)  :: expChem  ! Export State
   INTEGER, INTENT(OUT) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 - 

! !DESCRIPTION: Finalizes the DU Grid Component. Multiple instance
!               version.
!
! !REVISION HISTORY:
!
!  27Feb2008  da Silva  Introduced multiple instances
!
!EOP
!-------------------------------------------------------------------------

   integer :: i, ier

   do i = 1, gcDU%n
      call DU_SingleInstance_ ( DU_GridCompFinalize1_, i, &
                                gcDU%gcs(i), w_c, impChem, expChem, &
                                nymd, nhms, cdt, ier )
      if ( ier .NE. 0 ) then
         rc = i * 1000+ier
         return
      end if
   end do

   if (associated(gcDU%gcs)) deallocate ( gcDU%gcs, stat=ier )
   gcDU%n = -1

   end subroutine DU_GridCompFinalize


   subroutine DU_GridCompSetServices1_(  gc, chemReg, iname, rc)
   type(ESMF_GridComp), intent(INOUT) :: GC
   type(Chem_Registry), intent(INOUT) :: chemReg
   character(len=*),    intent(IN   ) :: iname
   integer,             intent(OUT  ) :: rc


   integer            :: status
   character(len=ESMF_MAXSTR) :: Iam

   Iam = "DU_GridCompSetServices1_"

   call MAPL_AddImportSpec(GC,              &
        SHORT_NAME = 'DU_SRC'//trim(iname), &
        LONG_NAME  = 'erod'  ,              &
        UNITS      = '1',                   &
        DIMS       = MAPL_DimsHorzOnly,     &
        VLOCATION  = MAPL_VLocationNone,    &
        RESTART    = MAPL_RestartSkip,      &
        RC         = STATUS)
   VERIFY_(STATUS)

  RETURN_(ESMF_SUCCESS)

   end subroutine DU_GridCompSetServices1_


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  DU_GridCompInitialize --- Initialize DU_GridComp
!
! !INTERFACE:
!

   subroutine DU_GridCompInitialize1_ ( gcDU, w_c, impChem, expChem, &
                                      nymd, nhms, cdt, rc )

! !USES:

  implicit NONE

! !INPUT PARAMETERS:

   type(Chem_Bundle), intent(inout) :: w_c     ! Chemical tracer fields      
   integer, intent(in) :: nymd, nhms           ! time
   real, intent(in) :: cdt                     ! chemistry timestep (secs)

! !OUTPUT PARAMETERS:

   type(DU_GridComp1), intent(inout) :: gcDU    ! Grid Component
   type(ESMF_State), intent(inout)   :: impChem ! Import State
   type(ESMF_State), intent(inout)   :: expChem ! Export State
   integer, intent(out) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 - 

! !DESCRIPTION: Initializes the DU Grid Component. It primarily sets
!               the import state for each active constituent package.
!
! !REVISION HISTORY:
!
!  18Sep2003 da Silva  First crack.
!
!EOP
!-------------------------------------------------------------------------

   character(len=*), parameter :: myname = 'DU_GridCompInitialize1_'


   character(len=255) :: rcfilen
   integer :: n
   integer :: i1, i2, im, j1, j2, jm, nbins, n1, n2, nbins_rc
   integer, allocatable :: ier(:)
   real :: qmax, qmin
   real :: radius, rlow, rup, rhop, fscav, fnum, molwght
   integer :: irhFlag
   integer :: imaringFlag

   integer, parameter :: nhres = 6   ! number of horizontal model resolutions: a,b,c,d,e
   real    :: Ch_DU(nhres)           ! emission tuning coefficient buffer


   rcfilen = trim(gcDU%rcfilen)
   gcDU%name = 'DU Constituent Package'

!  Initialize local variables
!  --------------------------
   rc = 0
   i1 = w_c%grid%i1; i2 = w_c%grid%i2; im = w_c%grid%im
   j1 = w_c%grid%j1; j2 = w_c%grid%j2; jm = w_c%grid%jm
   nbins = w_c%reg%n_DU
   n1  = w_c%reg%i_DU
   n2  = w_c%reg%j_DU

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

   call i90_label ( 'number_dust_bins:', ier(1) )
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
      radius               = i90_gfloat ( ier(n+1) )
      gcDU%radius(n)       = radius
      w_c%reg%rmed(n1+n-1) = radius * 1.e-6
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
      gcDU%rlow(n)          = rlow
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
      gcDU%rup(n)         = rup
      w_c%reg%rup(n1+n-1) = rup * 1.e-6
   end do
   if ( any(ier(1:nbins+1) /= 0) ) then
      call final_(50)
      return
   end if
!                          -------

!  Source fraction
!  ---------------
   call i90_label ( 'source_fraction:', ier(1) )
   do n = 1, nbins
      gcDU%sfrac(n) = i90_gfloat ( ier(n+1) )
   end do
   if ( any(ier(1:nbins+1) /= 0) ) then
      call final_(50)
      return
   end if
!                          -------

!  Soil Density
!  ---------------
   call i90_label ( 'soil_density:', ier(1) )
   do n = 1, nbins
      rhop                 = i90_gfloat ( ier(n+1) )
      gcDU%rhop(n)         = rhop
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
   gcDU%rhFlag                = irhFlag
   if ( any(ier(1:2) /= 0) ) then
      call final_(50)
      return
   end if
!                          -------

!  Dust emission tuning coefficient [kg s2 m-5]. NOT bin specific.
!  ---------------------------------------------------------------
   CALL I90_Label ( 'Ch_DU:', ier(1) )
   do n = 1, nhres
      Ch_DU(n) = i90_gfloat ( ier(n+1) )
   end do
   gcDU%Ch_DU = Chem_UtilResVal(im, jm, Ch_DU(:), ier(nhres + 2))
   gcDU%Ch_DU = gcDU%Ch_DU * 1.00E-09
   if ( any(ier(1:nhres+2) /= 0) ) then
      call final_(50)
      return
   end if

!  Settling velocity correction following Maring et al, 2003
!  ---------------
   call i90_label ( 'maringFlag:', ier(1) )
   imaringFlag = i90_gint ( ier(2) )
   if (imaringFlag /= 0) then
      gcDU%maringFlag = .True.
   else
      gcDU%maringFlag = .False.
   end if
   if ( any(ier(1:2) /= 0) ) then
      call final_(50)
      return
   end if

!  Handle Point-wise Emission Sources Specified in a Text File
!  -----------------------------------------------------------
   ier(:) = 0
   call i90_label  ( 'point_emissions_srcfilen:',   ier(1) )
   call i90_gtoken ( gcDU%point_emissions_srcfilen, ier(2) )
   if ( ier(1) /= 0 ) then
        gcDU%doing_point_emissions = .FALSE. ! if rc is missing, don't fuss
   else if ( any(ier(2:2) /= 0) ) then
         call final_(42) ! this means point emissions info is messed up, abort
         return
   else
         if ( (index(gcDU%point_emissions_srcfilen,'/dev/null')>0) ) then
               gcDU%doing_point_emissions = .FALSE. ! disable it if no file specified
         else
               gcDU%doing_point_emissions = .TRUE.  ! we are good to go
         end if
   end if


!                          -------
!  Initialize date for BCs
!  -----------------------
   gcDU%nymd = -1   ! nothing read yet

!  All done
!  --------
   call i90_release()
   deallocate(ier)

   return


CONTAINS

   subroutine init_()
   integer ios, nerr
   nerr = max ( 32, nbins+1 )
   allocate ( gcDU%radius(nbins), gcDU%src(i1:i2,j1:j2), &
              gcDU%rlow(nbins), gcDU%rup(nbins), &
              gcDU%sfrac(nbins), gcDU%rhop(nbins), ier(nerr), &
              stat=ios )
   if ( ios /= 0 ) rc = 100
   end subroutine init_

   subroutine final_(ierr)
   integer :: ierr
   integer ios
   deallocate ( gcDU%radius, gcDU%src, gcDU%sfrac, gcDU%rhop, &
                gcDU%rlow, gcDU%rup, ier, stat=ios )
   call i90_release()
   rc = ierr
   end subroutine final_

   end subroutine DU_GridCompInitialize1_



!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  DU_GridCompRun1_ --- The Chem Driver, run phase 1 
!
! !INTERFACE:
!

   subroutine DU_GridCompRun1_ ( gcDU, w_c, impChem, expChem, &
                                 nymd, nhms, cdt, rc )

! !USES:

  implicit NONE

! !INPUT/OUTPUT PARAMETERS:

   type(DU_GridComp1), intent(inout) :: gcDU    ! Grid Component
   type(Chem_Bundle), intent(inout)  :: w_c     ! Chemical tracer fields   

! !INPUT PARAMETERS:

   type(ESMF_State), intent(inout) :: impChem   ! Import State
   integer, intent(in) :: nymd, nhms            ! time
   real, intent(in) :: cdt                      ! chemistry timestep (secs)

! !OUTPUT PARAMETERS:

   type(ESMF_State), intent(inout) :: expChem   ! Export State
   integer, intent(out) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 -
 
! !DESCRIPTION: This routine implements the so-called DU Driver. That 
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

   character(len=*), parameter :: myname = 'DU_GridCompRun1_'
   character(len=*), parameter :: Iam = myname
   integer :: i1, i2, im, j1, j2, jm, nbins, n1, n2, km, n
   integer :: i, j, ijl, ijkl
   real :: qmax, qmin
   real, pointer :: DU_radius(:), DU_rhop(:)
   real, pointer :: emissions(:,:), dqa(:,:) 

!  Input fields from fvGCM
!  -----------------------
   real, pointer, dimension(:,:)   :: gwettop, oro, u10m, v10m, &
                                      ustar
   real, pointer, dimension(:,:,:) :: rhoa
   real, pointer, dimension(:,:)   :: frlake
   real, pointer, dimension(:,:)   :: du_src => null()

#define EXPORT        expChem
#define iNAME         TRIM(gcDU%iname)

#define ptrDUEM       DU_emis

   integer :: STATUS

!  Indices for point emissions
   integer, pointer, dimension(:)  :: iPoint, jPoint
   real, dimension(w_c%grid%km)    :: point_column_emissions
   integer                         :: ios, ii


#include "DU_GetPointer___.h"

!  Initialize local variables
!  --------------------------
   rc = 0
   i1 = w_c%grid%i1; i2 = w_c%grid%i2; im = w_c%grid%im
   j1 = w_c%grid%j1; j2 = w_c%grid%j2; jm = w_c%grid%jm

   km = w_c%grid%km
   nbins = w_c%reg%n_DU
   n1    = w_c%reg%i_DU
   n2    = w_c%reg%j_DU

   ijl  = ( i2 - i1 + 1 ) * ( j2 - j1 + 1 )
   ijkl = ijl * km

   if ( nbins /= NBIN_DUEM .OR. nbins /= NBIN_DUWT .OR. &
        nbins /= NBIN_DUDP .OR. nbins /= NBIN_DUSD ) then
      call die(myname,'inconsistent bins in resource file and registry')
   endif

! Update emissions/production if necessary (daily)
!  ------------------------------------------
   if(gcDU%nymd < 0) then

   call MAPL_GetPointer( impChem, du_src, 'DU_SRC'//iNAME, rc=status)
   VERIFY_(STATUS)
   gcDU%src = du_src

!   As a safety check, where du_src is undefined set to 0
!   -----------------------------------------------------
    do j = j1, j2
     do i = i1, i2
      if(1.01*gcDU%src(i,j) .gt. undefval) gcDU%src(i,j) = 0.
     enddo
    enddo

#ifdef DEBUG
    call pmaxmin('DU: src', gcDU%src, qmin, qmax, ijl, 1, 1. )
#endif

    gcDU%nymd = nymd

   endif

!  Read any pointwise emissions, if requested
!  ------------------------------------------
   if(gcDU%doing_point_emissions) then
    call Chem_UtilPointEmissions( nymd, gcDU%point_emissions_srcfilen, &
                                  gcDU%nPts, gcDU%pLat, gcDU%pLon, &
                                  gcDU%pBase, gcDU%pTop, gcDU%pEmis, &
                                  gcDU%pStart, gcDU%pEnd )

!   In case pStart or pEnd were not specified in the file set to defaults
    where(gcDU%pStart < 0) gcDU%pStart = 000000
    where(gcDU%pEnd < 0)   gcDU%pEnd   = 240000
   endif


!  Dust particle radius [m] and density [kg m-3]
!  ---------------------------------------------
   allocate( DU_radius(nbins), DU_rhop(nbins) )
   DU_radius = 1.e-6*gcDU%radius
   DU_rhop   = gcDU%rhop
   allocate( emissions(i1:i2,j1:j2), dqa(i1:i2,j1:j2), stat=STATUS)
   VERIFY_(STATUS)


#ifdef DEBUG
   do n = n1, n2
      call pmaxmin('DU: q_beg', w_c%qa(n)%data3d(i1:i2,j1:j2,1:km), qmin, qmax, &
                   ijl, km, 1. )
   end do
#endif

!  Get 2D Imports
!  --------------
   call MAPL_GetPointer ( impChem, frlake,   'FRLAKE',   __RC__ )
   call MAPL_GetPointer ( impChem, gwettop,  'WET1',     __RC__ )
   call MAPL_GetPointer ( impChem, oro,      'LWI',      __RC__ )
   call MAPL_GetPointer ( impChem, u10m,     'U10M',     __RC__ )
   call MAPL_GetPointer ( impChem, v10m,     'V10M',     __RC__ )
   call MAPL_GetPointer ( impChem, ustar,    'USTAR',    __RC__ )

!  Get 3D Imports
!  --------------
   call MAPL_GetPointer ( impChem, rhoa,     'AIRDENS',  __RC__ )

#ifdef DEBUG

   call pmaxmin('DU: frlake     ', frlake  , qmin, qmax, ijl,1, 1. )
   call pmaxmin('DU: gwtop      ', gwettop , qmin, qmax, ijl,1, 1. )
   call pmaxmin('DU: oro        ', oro     , qmin, qmax, ijl,1, 1. )
   call pmaxmin('DU: u10m       ', u10m    , qmin, qmax, ijl,1, 1. )
   call pmaxmin('DU: v10m       ', v10m    , qmin, qmax, ijl,1, 1. )
   call pmaxmin('DU: ustar      ', ustar   , qmin, qmax, ijl,1, 1. )

#endif

!  Dust Source
!  -----------
   do n = 1, nbins
       emissions = 0.0
       dqa = 0.0

       call DustEmissionGOCART( i1, i2, j1, j2, km, DU_radius(n), &
                                frlake, gwettop, oro, u10m, v10m, &
                                emissions, rc )

       dqa = gcDU%Ch_DU * gcDU%sfrac(n)*gcDU%src * emissions * cdt * grav / w_c%delp(:,:,km)

       w_c%qa(n1+n-1)%data3d(:,:,km) = w_c%qa(n1+n-1)%data3d(:,:,km) + dqa

       if (associated(DU_emis(n)%data2d)) then
           DU_emis(n)%data2d = gcDU%Ch_DU*gcDU%sfrac(n)*gcDU%src * emissions
       end if
   end do


#ifdef DEBUG
   do n = n1, n2
      call pmaxmin('DU: q_emi', w_c%qa(n)%data3d(i1:i2,j1:j2,1:km), qmin, qmax, &
                    ijl, km, 1. )
   end do
#endif


!  Distribute pointwise sources if requested
!  -----------------------------------------
   POINTWISE_SOURCES: if( gcDU%doing_point_emissions .and. gcDU%nPts > 0) then

!    Get indices for point emissions
!    -------------------------------
     allocate(iPoint(gcDU%nPts), jPoint(gcDU%nPts), stat=ios)

     call MAPL_GetHorzIJIndex(gcDU%nPts, iPoint, jPoint, &
                              grid = w_c%grid_esmf,      &
                              lon  = gcDU%pLon/radToDeg, &
                              lat  = gcDU%pLat/radToDeg, &
                              rc   = rc)

     if ( rc /= 0 ) call die(myname,'cannot get indices for point emissions')

     do ii = 1, gcDU%nPts
      i = iPoint(ii)
      j = jPoint(ii)
      if( i<1 .OR. j<1 )              cycle    ! point emission not in this sub-domain
!      if( gcDU%regionMask(i,j) == 0 ) cycle    ! masked by region mask

!     Emissions not occurring in current time step
!     --------------------------------------------
      if(nhms < gcDU%pStart(ii) .or. nhms >= gcDU%pEnd(ii)) cycle

      call distribute_point_emissions(w_c%delp(i,j,:), rhoa(i,j,:), &
                                      gcDU%pBase(ii), gcDU%pTop(ii), gcDU%pEmis(ii), &
                                      point_column_emissions, km)
      do n = 1, nbins

       w_c%qa(n1+n-1)%data3d(i,j,:) = w_c%qa(n1+n-1)%data3d(i,j,:) & 
          + cdt * grav / w_c%delp(i,j,:) &
                * gcDU%sfrac(n) * point_column_emissions / w_c%grid%cell_area(i,j)
      enddo

     enddo

     deallocate(iPoint, jPoint, stat=ios)

   endif POINTWISE_SOURCES


!  Clean up
!  --------
   deallocate ( DU_radius, DU_rhop, emissions, dqa, stat=STATUS )

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

 end subroutine DU_GridCompRun1_



!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  DU_GridCompRun2_ --- The Chem Driver, run phase 2 
!
! !INTERFACE:
!

   subroutine DU_GridCompRun2_ ( gcDU, w_c, impChem, expChem, &
                                 nymd, nhms, cdt, rc )

! !USES:

  implicit NONE

! !INPUT/OUTPUT PARAMETERS:

   type(DU_GridComp1), intent(inout) :: gcDU    ! Grid Component
   type(Chem_Bundle), intent(inout)  :: w_c     ! Chemical tracer fields   

! !INPUT PARAMETERS:

   type(ESMF_State), intent(inout) :: impChem   ! Import State
   integer, intent(in) :: nymd, nhms            ! time
   real, intent(in) :: cdt                      ! chemistry timestep (secs)

! !OUTPUT PARAMETERS:

   type(ESMF_State), intent(inout) :: expChem   ! Export State
   integer, intent(out) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 -
 
! !DESCRIPTION: This routine implements the so-called DU Driver. That 
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

   character(len=*), parameter :: myname = 'DU_GridCompRun2_'
   character(len=*), parameter :: Iam = myname
   integer :: i1, i2, im, j1, j2, jm, nbins, n1, n2, km, n
   integer :: k, ijl, ijkl
   real :: qmax, qmin
   real, pointer :: DU_radius(:), DU_rhop(:)
   real, pointer :: dqa(:,:), drydepositionfrequency(:,:)
   type(Chem_Array), pointer :: fluxout
   logical :: KIN


!  Input fields from fvGCM
!  -----------------------
   real, pointer, dimension(:,:)   ::  gwettop, oro, u10m, v10m, &
                                       ustar, precc, precl, pblh,          &
                                       shflux, z0h, hsurf, frocean, frseaice
   real, pointer, dimension(:,:,:) ::  tmpu, rhoa, u, v, hghte, ple
   real, pointer, dimension(:,:,:) ::  pfllsan, pfilsan

!  Additional needs for GOCART convective diagnostic
   real, pointer, dimension(:,:,:)     ::  cmfmc, qlcn, qicn, dtrain
   real, pointer, dimension(:,:)       ::  frlake, area
   real*8, allocatable, dimension(:,:,:) ::  cmfmc_, qccu_, dtrain_, &
                                             airmass_, airmol_, vud_, &
                                             delz_, delp_, tmpu_, ple_
   real*8, allocatable                   ::  tc_(:,:,:,:), bcnv_(:,:,:)
   real*8, allocatable                   ::  area_(:,:), frlake_(:,:), &
                                             frocean_(:,:), frseaice_(:,:)
   integer*4                             ::  icdt

#define EXPORT     expChem
#define iNAME      TRIM(gcDU%iname)

#define ptrDUWT       DU_wet
#define ptrDUSV       DU_conv
#define ptrDUEM       DU_emis
#define ptrDUDP       DU_dep
#define ptrDUSD       DU_set

#define    DUSMASS    DU_sfcmass
#define    DUCMASS    DU_colmass
#define    DUMASS     DU_mass
#define    DUEXTTAU   DU_exttau
#define    DUSCATAU   DU_scatau
#define    DUSMASS25  DU_sfcmass25
#define    DUCMASS25  DU_colmass25
#define    DUMASS25   DU_mass25
#define    DUEXTT25   DU_exttau25
#define    DUSCAT25   DU_scatau25
#define    DUAERIDX   DU_aeridx
#define    DUFLUXU    DU_fluxu
#define    DUFLUXV    DU_fluxv
#define    DUCONC     DU_conc
#define    DUEXTCOEF  DU_extcoef
#define    DUSCACOEF  DU_scacoef
#define    DUEXTTFM   DU_exttaufm
#define    DUSCATFM   DU_scataufm
#define    DUANGSTR   DU_angstrom

   integer :: STATUS

#include "DU_GetPointer___.h"

!  Initialize local variables
!  --------------------------
   rc = 0
   i1 = w_c%grid%i1; i2 = w_c%grid%i2; im = w_c%grid%im
   j1 = w_c%grid%j1; j2 = w_c%grid%j2; jm = w_c%grid%jm

   km = w_c%grid%km
   nbins = w_c%reg%n_DU
   n1    = w_c%reg%i_DU
   n2    = w_c%reg%j_DU

   ijl  = ( i2 - i1 + 1 ) * ( j2 - j1 + 1 )
   ijkl = ijl * km

   if ( nbins /= NBIN_DUEM .OR. nbins /= NBIN_DUWT .OR. &
        nbins /= NBIN_DUDP .OR. nbins /= NBIN_DUSD ) then
      call die(myname,'inconsistent bins in resource file and registry')
   endif


#ifdef DEBUG
   do n = n1, n2
      call pmaxmin('DU: q_beg', w_c%qa(n)%data3d(i1:i2,j1:j2,1:km), qmin, qmax, &
                   ijl, km, 1. )
   end do
#endif

!  Get 2D Imports
!  --------------
   call MAPL_GetPointer ( impChem, frlake,   'FRLAKE',   __RC__ )
   call MAPL_GetPointer ( impChem, gwettop,  'WET1',     __RC__ )
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

   call pmaxmin('DU: frlake     ', frlake  , qmin, qmax, ijl,1, 1. )
   call pmaxmin('DU: gwtop      ', gwettop , qmin, qmax, ijl,1, 1. )
   call pmaxmin('DU: oro        ', oro     , qmin, qmax, ijl,1, 1. )
   call pmaxmin('DU: u10m       ', u10m    , qmin, qmax, ijl,1, 1. )
   call pmaxmin('DU: v10m       ', v10m    , qmin, qmax, ijl,1, 1. )
   call pmaxmin('DU: ustar      ', ustar   , qmin, qmax, ijl,1, 1. )
   call pmaxmin('DU: precc      ', precc   , qmin, qmax, ijl,1, 1. )
   call pmaxmin('DU: precl      ', precl   , qmin, qmax, ijl,1, 1. )
   call pmaxmin('DU: pblh       ', pblh    , qmin, qmax, ijl,1, 1. )
   call pmaxmin('DU: shfflux    ', shflux  , qmin, qmax, ijl,1, 1. )
   call pmaxmin('DU: z0h        ', z0h     , qmin, qmax, ijl,1, 1. )
   call pmaxmin('DU: hsurf      ', hsurf   , qmin, qmax, ijl,1, 1. )

   call pmaxmin('DU: tmpu       ', tmpu    , qmin, qmax, ijl,km, 1. )
   call pmaxmin('DU: rhoa       ', rhoa    , qmin, qmax, ijl,km, 1. )
   call pmaxmin('DU: u          ', u       , qmin, qmax, ijl,km, 1. )
   call pmaxmin('DU: v          ', v       , qmin, qmax, ijl,km, 1. )
   call pmaxmin('DU: hghte      ', hghte   , qmin, qmax, ijl,km, 1. )
   call pmaxmin('DU: rh         ', w_c%rh  , qmin, qmax, ijl,km, 1. )
   call pmaxmin('DU: pfllsan    ', pfllsan , qmin, qmax, ijl,km+1, 1. )
   call pmaxmin('DU: pfilsan    ', pfilsan , qmin, qmax, ijl,km+1, 1. )

#endif


RUN_ALARM: if (gcDU%run_alarm) then

!  Dust particle radius [m] and density [kg m-3]
!  ---------------------------------------------
   allocate( DU_radius(nbins), DU_rhop(nbins) )
   DU_radius = 1.e-6*gcDU%radius
   DU_rhop   = gcDU%rhop

   allocate( fluxout )
   allocate( fluxout%data2d(i1:i2,j1:j2), dqa(i1:i2,j1:j2), &
             drydepositionfrequency(i1:i2,j1:j2), stat=STATUS)
   VERIFY_(STATUS)

!  Dust Settling
!  -----------
   call Chem_Settling ( i1, i2, j1, j2, km, n1, n2, nbins, gcDU%rhFlag, &
                        DU_radius, DU_rhop, cdt, w_c, tmpu, rhoa, hsurf,    &
                        hghte, DU_set, rc, correctionMaring=gcDU%maringFlag )

#ifdef DEBUG
   do n = n1, n2
      call pmaxmin('DU: q_set', w_c%qa(n)%data3d(i1:i2,j1:j2,1:km), qmin, qmax, &
                    ijl, km, 1. )
   end do
#endif

!  Dust Deposition
!  -----------
   do n = 1, nbins
    drydepositionfrequency = 0.
    call DryDepositionGOCART( i1, i2, j1, j2, km, &
                              tmpu, rhoa, hghte, oro, ustar, &
                              pblh, shflux, z0h, drydepositionfrequency, rc, &
                              DU_radius(n), DU_rhop(n), u10m, v10m, frlake, gwettop )
    
    dqa = 0.
    dqa = max(0.0, w_c%qa(n1+n-1)%data3d(:,:,km)*(1.-exp(-drydepositionfrequency*cdt)))
    w_c%qa(n1+n-1)%data3d(:,:,km) = &
            w_c%qa(n1+n-1)%data3d(:,:,km) - dqa
    if( associated(DU_dep(n)%data2d) ) &
     DU_dep(n)%data2d = dqa*w_c%delp(:,:,km)/grav/cdt
   end do

#ifdef DEBUG
   do n = n1, n2
      call pmaxmin('DU: q_dry', w_c%qa(n)%data3d(i1:i2,j1:j2,1:km), qmin, qmax, &
                    ijl, km, 1. )
   end do
#endif

!  Dust Large-scale Wet Removal
!  ----------------------------
   KIN = .TRUE. 
   do n = 1, nbins
   !w_c%qa(n1+n-1)%fwet = 1.0   ! GEOS-Chem
    w_c%qa(n1+n-1)%fwet = 0.3
    call WetRemovalGOCART(i1, i2, j1, j2, km, n1+n-1, n1+n-1, cdt, 'dust', KIN, &
                          w_c%qa, ple, tmpu, rhoa, pfllsan, pfilsan, &
                          precc, precl, fluxout, rc )
    if(associated(DU_wet(n)%data2d)) DU_wet(n)%data2d = fluxout%data2d
   end do

#ifdef DEBUG
   do n = n1, n2
      call pmaxmin('DU: q_wet', w_c%qa(n)%data3d(i1:i2,j1:j2,1:km), qmin, qmax, &
                    ijl, km, 1. )
   end do
#endif

!  Dust Convective-scale Mixing and Wet Removal
!  --------------------------------------------
   KIN = .TRUE.
   icdt = cdt

   allocate(cmfmc_(i1:i2,j1:j2,km+1), qccu_(i1:i2,j1:j2,km), &
            dtrain_(i1:i2,j1:j2,km), airmass_(i1:i2,j1:j2,km), &
            delz_(i1:i2,j1:j2,km), vud_(i1:i2,j1:j2,km), &
            tc_(i1:i2,j1:j2,km,n1:n2), delp_(i1:i2,j1:j2,km), &
            airmol_(i1:i2,j1:j2,km), tmpu_(i1:i2,j1:j2,km),&
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
   call convection(i1, i2, j1, j2, km, n1, n2, icdt, 'dust', kin, &
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
    if(associated(DU_conv(n)%data2d)) DU_conv(n)%data2d = -bcnv_(:,:,n1+n-1)/area_/icdt
   end do

!  Clean up
!  --------
   deallocate(cmfmc_, qccu_, dtrain_, tc_, airmass_, &
              delz_, vud_, delp_, airmol_, tmpu_, bcnv_, ple_, &
              area_, frlake_, frocean_, frseaice_, __STAT__ )

   deallocate ( fluxout%data2d )
   deallocate ( fluxout, DU_radius, DU_rhop, &
                dqa, drydepositionfrequency, stat=STATUS )

   end if RUN_ALARM

!  Compute the desired output diagnostics here
!  Ideally this will go where chemout is called in fvgcm.F since that
!  will reflect the distributions after transport, etc.
!  -----------
   call DU_Compute_Diags(i1, i2, j1, j2, km, nbins, gcDU, w_c, tmpu, rhoa,    &
                         u, v, DU_sfcmass,  DU_colmass, DU_mass, DU_exttau,   &
                         DU_scatau,   DU_sfcmass25, DU_colmass25, DU_mass25,  &
                         DU_exttau25, DU_scatau25,  DU_aeridx, DU_fluxu,      &
                         DU_fluxv, DU_conc, DU_extcoef, DU_scacoef,           &
                         DU_exttaufm, DU_scataufm, DU_angstrom, rc)

   return

CONTAINS

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  DU_Compute_Diags - Calculate dust 2D diagnostics
!
! !INTERFACE:
!

   subroutine DU_Compute_Diags ( i1, i2, j1, j2, km, nbins, gcDU, w_c, tmpu, rhoa, &
                                 u, v, sfcmass, colmass, mass, exttau, scatau,     &
                                 sfcmass25, colmass25, mass25, exttau25, scatau25, &
                                 aerindx, fluxu, fluxv, conc, extcoef, scacoef,    &
                                 exttaufm, scataufm, angstrom, rc )

! !USES:

  implicit NONE

! !INPUT PARAMETERS:
   integer, intent(in) :: i1, i2, j1, j2, km, nbins
   type(DU_GridComp1), intent(inout) :: gcDU   ! DU Grid Component
   type(Chem_Bundle), intent(in)   :: w_c      ! Chem Bundle
   real, pointer, dimension(:,:,:) :: tmpu     ! temperature [K]
   real, pointer, dimension(:,:,:) :: rhoa     ! air density [kg m-3]
   real, pointer, dimension(:,:,:) :: u        ! east-west wind [m s-1]
   real, pointer, dimension(:,:,:) :: v        ! north-south wind [m s-1]
   

! !OUTPUT PARAMETERS:
!  Total mass
   type(Chem_Array), intent(inout)  :: sfcmass   ! sfc mass concentration kg/m3
   type(Chem_Array), intent(inout)  :: colmass   ! col mass density kg/m2
   type(Chem_Array), intent(inout)  :: mass      ! 3d mass mixing ratio kg/kg
!  Total optical properties
   type(Chem_Array), intent(inout)  :: exttau    ! ext. AOT at 550 nm
   type(Chem_Array), intent(inout)  :: scatau    ! sct. AOT at 550 nm
   type(Chem_Array), intent(inout)  :: sfcmass25 ! sfc mass concentration kg/m3 (pm2.5)
   type(Chem_Array), intent(inout)  :: colmass25 ! col mass density kg/m2 (pm2.5)
   type(Chem_Array), intent(inout)  :: mass25    ! 3d mass mixing ratio kg/kg (pm2.5)
   type(Chem_Array), intent(inout)  :: exttau25  ! ext. AOT at 550 nm (pm2.5)
   type(Chem_Array), intent(inout)  :: scatau25  ! sct. AOT at 550 nm (pm2.5)
   type(Chem_Array), intent(inout)  :: aerindx   ! TOMS UV AI
   type(Chem_Array), intent(inout)  :: fluxu     ! Column mass flux in x direction
   type(Chem_Array), intent(inout)  :: fluxv     ! Column mass flux in y direction
   type(Chem_Array), intent(inout)  :: conc      ! 3d mass concentration, kg/m3
   type(Chem_Array), intent(inout)  :: extcoef   ! 3d ext. coefficient, 1/m
   type(Chem_Array), intent(inout)  :: scacoef   ! 3d scat.coefficient, 1/m
   type(Chem_Array), intent(inout)  :: exttaufm  ! fine mode (sub-micron) ext. AOT at 550 nm
   type(Chem_Array), intent(inout)  :: scataufm  ! fine mode (sub-micron) sct. AOT at 550 nm
   type(Chem_Array), intent(inout)  :: angstrom  ! 470-870 nm Angstrom parameter
   integer, intent(out)             :: rc        ! Error return code:
                                                 !  0 - all is well
                                                 !  1 - 

! !DESCRIPTION: Calculates some simple 2d diagnostics from the dust fields
!
! !REVISION HISTORY:
!
!  16APR2004, Colarco
!  11MAR2010, Nowottnick  
!
!EOP
!-------------------------------------------------------------------------

! !Local Variables
   character(len=*), parameter :: myname = 'DU_Compute_Diags'
   integer :: i, j, k, n, n1, n2, ios, nch, idx
   real :: ilam550, ilam470, ilam870
   real :: tau, ssa
   real :: fPMfm(nbins)  ! fraction of bin with particles diameter < 1.0 um
   real :: fPM25(nbins)  ! fraction of bin with particles diameter < 2.5 um
   character(len=255) :: qname
   logical :: do_angstrom
   real, dimension(i1:i2,j1:j2) :: tau470, tau870

!  Initialize local variables
!  --------------------------
   n1    = w_c%reg%i_DU
   n2    = w_c%reg%j_DU
   nch   = gcDU%mie_tables%nch

!  Get the wavelength indices
!  --------------------------
!  Must provide ilam550 for AOT calculation
   ilam550 = 1.
   ilam470 = 0.
   ilam870 = 0.
   if(nch .gt. 1) then
    do i = 1, nch
     if ( gcDU%mie_tables%channels(i) .ge. 5.49e-7 .and. &
          gcDU%mie_tables%channels(i) .le. 5.51e-7) ilam550 = i
     if ( gcDU%mie_tables%channels(i) .ge. 4.69e-7 .and. &
          gcDU%mie_tables%channels(i) .le. 4.71e-7) ilam470 = i
     if ( gcDU%mie_tables%channels(i) .ge. 8.69e-7 .and. &
          gcDU%mie_tables%channels(i) .le. 8.71e-7) ilam870 = i
    enddo
   endif

   do_angstrom = .false.
!  If both 470 and 870 channels provided (and not the same) then
!  possibly will do Angstrom parameter calculation
   if(ilam470 .ne. 0. .and. &
      ilam870 .ne. 0. .and. &
      ilam470 .ne. ilam870) do_angstrom = .true.

!  Compute the fine mode (sub-micron) and PM2.5 bin-wise fractions
!  ------------------------------------
   call DU_Binwise_PM_Fractions(fPMfm, 0.50, gcDU%rlow, gcDU%rup, nbins)   ! 2*r < 1.0 um
   call DU_Binwise_PM_Fractions(fPM25, 1.25, gcDU%rlow, gcDU%rup, nbins)   ! 2*r < 2.5 um


   if ( associated(aerindx%data2d) )  aerindx%data2d = 0.0  ! for now

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
       qname = trim(w_c%reg%vname(w_c%reg%i_DU+n-1))
       idx = Chem_MieQueryIdx(gcDU%mie_tables,qname,rc)
       if(rc .ne. 0) call die(myname, 'cannot find proper Mie table index')

       do k = 1, km
        do j = j1, j2
         do i = i1, i2
          call Chem_MieQuery(gcDU%mie_tables, idx, ilam550, &
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
       qname = trim(w_c%reg%vname(w_c%reg%i_DU+n-1))
       idx = Chem_MieQueryIdx(gcDU%mie_tables,qname,rc)
       if(rc .ne. 0) call die(myname, 'cannot find proper Mie table index')

       do k = 1, km
        do j = j1, j2
         do i = i1, i2

          call Chem_MieQuery(gcDU%mie_tables, idx, ilam470, &
              w_c%qa(n1+n-1)%data3d(i,j,k)*w_c%delp(i,j,k)/grav, &
              w_c%rh(i,j,k), tau=tau)
          tau470(i,j) = tau470(i,j) + tau

          call Chem_MieQuery(gcDU%mie_tables, idx, ilam870, &
              w_c%qa(n1+n-1)%data3d(i,j,k)*w_c%delp(i,j,k)/grav, &
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

   end subroutine DU_Compute_Diags


!##############################################################################
!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  DU_Binwise_PM_Fractions - Calculate bin-wise PM fractions
!
! !INTERFACE:
!

   subroutine DU_Binwise_PM_Fractions(fPM, rPM, r_low, r_up, nbins)

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
!EOP

! !Local Variables

   integer :: n

   character(len=*), parameter :: myname = 'DU_Binwise_PM_Fractions'

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

   end subroutine DU_Binwise_PM_Fractions

 end subroutine DU_GridCompRun2_

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  DU_GridCompFinalize --- The Chem Driver 
!
! !INTERFACE:
!

   subroutine DU_GridCompFinalize1_ ( gcDU, w_c, impChem, expChem, &
                                    nymd, nhms, cdt, rc )

! !USES:

  implicit NONE

! !INPUT/OUTPUT PARAMETERS:

   type(DU_GridComp1), intent(inout) :: gcDU  ! Grid Component

! !INPUT PARAMETERS:

   type(Chem_Bundle), intent(in)  :: w_c      ! Chemical tracer fields   
   integer, intent(in) :: nymd, nhms          ! time
   real, intent(in) :: cdt                    ! chemistry timestep (secs)

! !OUTPUT PARAMETERS:

   type(ESMF_State), intent(inout) :: impChem   ! Import State
   type(ESMF_State), intent(inout) :: expChem   ! Export State
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

    rc=0
!   integer :: ios

!   deallocate ( gcDU%radius, gcDU%src, stat=ios )
!   if ( ios /= 0 ) then
!      rc = 1
!      return
!   end if

   return

 end subroutine DU_GridCompFinalize1_

 end module DU_GridCompMod

!-----------------------------------------------------------------------

!                     Single Instance Wrapper

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  DU_SingleInstance_ --- Runs single instance of method
!
! !INTERFACE:
!
  subroutine DU_SingleInstance_ ( Method_, instance, &
                                  gcDU, w_c, impChem, expChem, &
                                  nymd, nhms, cdt, rc )

! !USES:

  Use DU_GridCompMod
  Use ESMF
  Use MAPL
  Use Chem_Mod 

  IMPLICIT NONE

! !INPUT PARAMETERS:

!  Input "function pointer"
!  -----------------------
   interface 
     subroutine Method_ (gc, w, imp, exp, ymd, hms, dt, rcode )
       Use DU_GridCompMod
       Use ESMF
       Use MAPL
       Use Chem_Mod 
       type(DU_GridComp1),  intent(inout)  :: gc
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

   TYPE(DU_GridComp1), INTENT(INOUT) :: gcDU    ! Grid Component
   TYPE(ESMF_State), INTENT(INOUT)  :: impChem  ! Import State
   TYPE(ESMF_State), INTENT(INOUT)  :: expChem  ! Export State
   INTEGER, INTENT(OUT) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 - 

! !DESCRIPTION: Finalizes the DU Grid Component. Multiple instance
!               version.
!
! !REVISION HISTORY:
!
!  27Feb2008  da Silva  Introduced multiple instances
!
!EOP
!-------------------------------------------------------------------------

  integer :: n, n_DU, i_DU, j_DU
  integer :: status
  character(len=255), allocatable :: qname(:)
  character(len=ESMF_MAXSTR) :: Iam
  integer, parameter :: n_bins = 5

  Iam = 'DU_SingleInstance_'

! Save overall DU indices
! -----------------------
  n_DU = w_c%reg%n_DU
  i_DU = w_c%reg%i_DU
  j_DU = w_c%reg%j_DU

! Save the name of the variables in this instance
! -----------------------------------------------
  allocate(qname(n_bins), __STAT__)

  do n = 1, n_bins
      qname(n) = trim(w_c%reg%vname(i_DU + n_bins*(instance - 1) + n - 1))
  end do
  
! Customize indices for this particular instance
! ----------------------------------------------
  w_c%reg%n_DU = n_bins
  w_c%reg%i_DU = i_DU + n_bins*(instance - 1)
  w_c%reg%j_DU = i_DU + n_bins*(instance - 1) + (n_bins - 1)

  do n = 1, n_bins
      w_c%reg%vname(i_DU + n_bins*(instance - 1) + n - 1) = w_c%reg%vname(i_DU + n - 1)
  end do
  
! Execute the instance method
! ---------------------------
  call Method_ ( gcDU, w_c, impChem, expChem, &
                 nymd, nhms, cdt, rc )

! Restore the overall DU indices
! ------------------------------
  do n = 1, n_bins
      w_c%reg%vname(i_DU + n_bins*(instance - 1) + n - 1) = qname(n)
  end do

  w_c%reg%n_DU = n_DU
  w_c%reg%i_DU = i_DU
  w_c%reg%j_DU = j_DU

  deallocate(qname, __STAT__)

  end subroutine DU_SingleInstance_

!-----------------------------------------------------------------------

