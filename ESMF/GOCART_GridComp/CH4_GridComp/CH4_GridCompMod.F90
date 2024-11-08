#include "MAPL_Generic.h"
!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  CH4_GridCompMod --- CH4 Grid Component Class
!
! !INTERFACE:
!

   MODULE  CH4_GridCompMod

! !USES:

   USE ESMF
   USE MAPL
   USE Chem_Mod              ! Chemistry Base Class
   USE Chem_StateMod         ! Chemistry State
   USE Chem_UtilMod          ! I/O
   USE m_inpak90             ! Resource file management
   USE Henrys_law_ConstantsMod, ONLY: get_HenrysLawCts
   USE m_chars, ONLY: lowercase, uppercase

   IMPLICIT NONE

! !PUBLIC TYPES:
!
   PRIVATE
   PUBLIC  CH4_GridComp       ! Multiple instance CH4 object 
   PUBLIC  CH4_GridComp1      ! Single instance CH4 object

!
! !PUBLIC MEMBER FUNCTIONS:
!

   PUBLIC  CH4_GridCompSetServices
   PUBLIC  CH4_GridCompInitialize
   PUBLIC  CH4_GridCompRun
   PUBLIC  CH4_GridCompFinalize

!
! !DESCRIPTION:
!
!  This module implements the CH4 Grid Component. 
!
! !REVISION HISTORY:
!
!  24 Jun 2010 Nielsen:  First crack.
!  25 Oct 2012 Nielsen:  Added photolysis.
!  06 Oct 2017 Weir:     Fixing whitespace (STOP USING TABS IN FORTRAN!).
!
!EOP
!-------------------------------------------------------------------------

  TYPE CH4_GridComp1

   CHARACTER(LEN=ESMF_MAXSTR) :: name            ! generic name of the package
   CHARACTER(LEN=ESMF_MAXSTR) :: iname           ! instance name
   CHARACTER(LEN=ESMF_MAXSTR) :: rcfilen         ! resource file name

   INTEGER :: instance                 ! Instantiation number

   REAL, POINTER ::      CH4(:,:,:)    ! CH4 mixing ratio mol/mol
   REAL, POINTER ::     OHnd(:,:,:)    ! OH number density (cm^{-3})
   REAL, POINTER ::     Clnd(:,:,:)    ! Cl number density (cm^{-3})
   REAL, POINTER ::    O1Dnd(:,:,:)    ! O(1D) number density (cm^{-3})

   REAL, POINTER :: eCH4_agw(:,:)      ! kgCH4/m2/s, Earth surface
   REAL, POINTER :: eCH4_ind(:,:)      ! kgCH4/m2/s, Earth surface
   REAL, POINTER :: eCH4_wetl(:,:)     ! kgCH4/m2/s, Earth surface
   REAL, POINTER :: eCH4_mnat(:,:)     ! kgCH4/m2/s, Earth surface
   REAL, POINTER :: eCH4_bb(:,:)       ! kgCH4/m2/s, PBL (before diurnal)
   REAL, POINTER :: eCH4_bb_(:,:)      ! kgCH4/m2/s, PBL
   REAL, POINTER :: eCH4_bf(:,:)       ! kgCH4/m2/s, Earth surface

   LOGICAL :: DebugIsOn     ! Run-time debug switch
   LOGICAL :: CH4FeedBack   ! Permit increments to CH4 from CH4 + hv => 2H2O + CO
   LOGICAL :: H2OFeedBack   ! Permit increments to   Q from CH4 + hv => 2H2O + CO
   CHARACTER(LEN=ESMF_MAXSTR) :: units_oh ! Units for OH 

   REAL :: szaCutoff        ! Largest solar zenith angle (degrees) allowed as daytime

  END TYPE CH4_GridComp1

  TYPE CH4_GridComp
     INTEGER                      ::  n        ! number of instances 
     TYPE(CH4_GridComp1), pointer ::  gcs(:)   ! instances
  END TYPE CH4_GridComp

  real, parameter :: radToDeg = 180./MAPL_PI
  real, parameter :: mwtCH4   = 16.0422

CONTAINS

   subroutine CH4_GridCompSetServices(  gc, chemReg, rc)
   type(ESMF_GridComp), intent(INOUT) :: GC
   type(Chem_Registry), intent(INOUT) :: chemReg
   integer,             intent(OUT  ) :: rc

   CHARACTER(LEN=255) :: rcbasen = 'CH4_GridComp'
   CHARACTER(LEN=255) :: name

   integer            :: status
   character(len=ESMF_MAXSTR) :: Iam

   integer :: ier,n,i

   type(ESMF_Config) :: cfg

   Iam = "CH4_GridCompSetServices"

!  Load resource file
!  ------------------
   cfg = ESMF_ConfigCreate(rc=status)
   VERIFY_(STATUS)
   call ESMF_ConfigLoadFile(cfg,TRIM(rcbasen)//'.rc',rc=status)
   VERIFY_(STATUS)

!  Parse resource file
!  -------------------
   n = ESMF_ConfigGetLen(cfg,label='CH4_instances:',rc=status)
   VERIFY_(STATUS)

!  We cannot have fewer instances than the number of
!  CH4 bins in the registry (it is OK to have less, though)
!  --------------------------------------------------------
   if ( n .LT. chemReg%n_CH4 ) then
        rc = 35
        return
   else if ( n .GT. chemReg%n_CH4 ) then
        if (MAPL_AM_I_ROOT()) &
        PRINT *, TRIM(Iam)//": Bins = ",chemReg%n_CH4," of ",n," expected."
   end if
   n = min(n, chemReg%n_CH4)

!  Record name of each instance
!  ----------------------------
   call ESMF_ConfigFindLabel(cfg,'CH4_instances:',rc=status)
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
      call CH4_GridCompSetServices1_(gc,chemReg,name,rc=status)
      VERIFY_(STATUS)
   end do

   RETURN_(ESMF_SUCCESS)

   end subroutine CH4_GridCompSetServices

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  CH4_GridCompInitialize --- Initialize CH4_GridComp
!
! !INTERFACE:
!

   SUBROUTINE CH4_GridCompInitialize ( gcCH4, w_c, impChem, expChem, nymd, nhms, cdt, rc )

! !USES:

  IMPLICIT NONE

! !INPUT PARAMETERS:

   TYPE(Chem_Bundle), intent(in) :: w_c        ! Chemical tracer fields      
   INTEGER, INTENT(IN) :: nymd, nhms           ! time
   REAL,    INTENT(IN) :: cdt                  ! chemical timestep (secs)


! !OUTPUT PARAMETERS:

   TYPE(CH4_GridComp), INTENT(INOUT) :: gcCH4   ! Grid Component
   TYPE(ESMF_State), INTENT(INOUT)  :: impChem  ! Import State
   TYPE(ESMF_State), INTENT(INOUT)  :: expChem  ! Export State
   INTEGER, INTENT(OUT) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 - 

! !DESCRIPTION: Initializes the CH4 Grid Component. Multiple instance
!               version.
!
! !REVISION HISTORY:
!
!  24 Jun 2010 Nielsen:  First crack.
!  25 Oct 2012 Nielsen:  Added photolysis.
!
!EOP
!-------------------------------------------------------------------------

   CHARACTER(LEN=*), PARAMETER :: Iam = 'CH4_GridCompInitialize'
   CHARACTER(LEN=ESMF_MAXSTR) :: rcbasen = 'CH4_GridComp'
   CHARACTER(LEN=ESMF_MAXSTR) :: name
   
   integer :: i, n, status
   real :: c1, c2, c3, c4

!  Load resource file
!  ------------------
   CALL I90_loadf ( TRIM(rcbasen)//'.rc', status )
   VERIFY_(status)

!  Parse resource file
!  -------------------
   CALL I90_label ( 'CH4_instances:', status )
   VERIFY_(status)

!  First determine how many instances we have
!  ------------------------------------------   
   n = 0
   status = 0

   DO WHILE ( status == 0 )
    CALL I90_gtoken( name, status )
    IF(status == 0) n = n + 1
   END DO
   IF ( n == 0 ) THEN
    status = 1
    VERIFY_(status)
   END IF

!  We cannot have fewer instances than the number of
!  CH4 bins in the registry (it is OK to have less, though)
!  --------------------------------------------------------
   IF ( n < w_c%reg%n_CH4 ) THEN
    status = 1
    VERIFY_(status)
   ELSE IF ( n >= w_c%reg%n_CH4 ) THEN
    IF(MAPL_AM_I_ROOT()) PRINT *, TRIM(Iam)//": Bins = ",w_c%reg%n_CH4," of ",n," expected."
   END IF
   n = min(n,w_c%reg%n_CH4 )
   gcCH4%n = n

!  Next allocate necessary memory
!  ------------------------------
   ALLOCATE ( gcCH4%gcs(n), STAT=status )    
   VERIFY_(status)

!  Record name of each instance
!  ----------------------------
   CALL I90_label ( 'CH4_instances:', status )
   VERIFY_(status)
   DO i = 1, n
      CALL I90_gtoken( name, status )
      VERIFY_(status)
                                            ! resource file name
      gcCH4%gcs(i)%rcfilen = trim(rcbasen)//'---'//trim(name)//'.rc'
      gcCH4%gcs(i)%instance = i              ! instance number 
      IF(TRIM(name) == "full" ) THEN
       gcCH4%gcs(i)%iname = " "              ! blank instance name for full (1)
      ELSE
       gcCH4%gcs(i)%iname = TRIM(name)       ! instance name for others
      END IF
   END DO    

!  Next initialize each instance
!  -----------------------------
   DO i = 1, gcCH4%n
      IF(MAPL_AM_I_ROOT()) THEN
       PRINT *," "
       PRINT *,TRIM(Iam)//": Initializing instance ",TRIM(gcCH4%gcs(i)%iname)," [",gcCH4%gcs(i)%instance,"]"
      END IF
      CALL CH4_SingleInstance_ ( CH4_GridCompInitialize1_, i, &
                                gcCH4%gcs(i), w_c, impChem, expChem,  &
                                nymd, nhms, cdt, status )
      VERIFY_(status)
   END DO

!  Get Henrys Law cts for the parameterized convective wet removal
!  ---------------------------------------------------------------
   CALL get_HenrysLawCts('CH4', c1, c2, c3, c4)
   w_c%reg%Hcts(1,w_c%reg%i_CH4:w_c%reg%j_CH4) = c1
   w_c%reg%Hcts(2,w_c%reg%i_CH4:w_c%reg%j_CH4) = c2
   w_c%reg%Hcts(3,w_c%reg%i_CH4:w_c%reg%j_CH4) = c3
   w_c%reg%Hcts(4,w_c%reg%i_CH4:w_c%reg%j_CH4) = c4

!  All done
!  --------
   CALL I90_FullRelease( status )
   VERIFY_(status)

 END SUBROUTINE CH4_GridCompInitialize

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  CH4_GridCompRun --- Run CH4_GridComp
!
! !INTERFACE:
!

   SUBROUTINE CH4_GridCompRun ( gcCH4, w_c, impChem, expChem, nymd, nhms, cdt, rc )

! !USES:

  IMPLICIT NONE

! !INPUT PARAMETERS:

   TYPE(Chem_Bundle), intent(in) :: w_c        ! Chemical tracer fields      
   INTEGER, INTENT(IN) :: nymd, nhms           ! time
   REAL,    INTENT(IN) :: cdt                  ! chemical timestep (secs)


! !OUTPUT PARAMETERS:

   TYPE(CH4_GridComp), INTENT(INOUT) :: gcCH4   ! Grid Component
   TYPE(ESMF_State), INTENT(INOUT)  :: impChem  ! Import State
   TYPE(ESMF_State), INTENT(INOUT)  :: expChem  ! Export State
   INTEGER, INTENT(OUT) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 - 

! !DESCRIPTION: Runs the CH4 Grid Component. Multiple instance
!               version.
!
! !REVISION HISTORY:
!
!  24 Jun 2010 Nielsen:  First crack.
!  25 Oct 2012 Nielsen:  Added photolysis.
!
!EOP
!-------------------------------------------------------------------------

   CHARACTER(LEN=ESMF_MAXSTR), PARAMETER :: Iam = 'CH4_GridCompRun'
   INTEGER :: i, status

   DO i = 1, gcCH4%n
      CALL CH4_SingleInstance_ ( CH4_GridCompRun1_, i, &
                                gcCH4%gcs(i), w_c, impChem, expChem, &
                                nymd, nhms, cdt, status )
      VERIFY_(status)
   END DO

 END SUBROUTINE CH4_GridCompRun


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  CH4_GridCompFinalize --- Initialize CH4_GridComp
!
! !INTERFACE:
!

   SUBROUTINE CH4_GridCompFinalize ( gcCH4, w_c, impChem, expChem, nymd, nhms, cdt, rc )

! !USES:

  IMPLICIT NONE

! !INPUT PARAMETERS:

   TYPE(Chem_Bundle), intent(in) :: w_c        ! Chemical tracer fields      
   INTEGER, INTENT(IN) :: nymd, nhms           ! time
   REAL,    INTENT(IN) :: cdt                  ! chemical timestep (secs)


! !OUTPUT PARAMETERS:

   TYPE(CH4_GridComp), INTENT(INOUT) :: gcCH4   ! Grid Component
   TYPE(ESMF_State), INTENT(INOUT)  :: impChem  ! Import State
   TYPE(ESMF_State), INTENT(INOUT)  :: expChem  ! Export State
   INTEGER, INTENT(OUT) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 - 

! !DESCRIPTION: Finalizes the CH4 Grid Component. Multiple instance
!               version.
!
! !REVISION HISTORY:
!
!  27Feb2008  da Silva  Introduced multiple instances
!
!EOP
!-------------------------------------------------------------------------

   CHARACTER(LEN=ESMF_MAXSTR), PARAMETER :: Iam = 'CH4_GridCompFinalize'
   INTEGER :: i, status

   DO i = 1, gcCH4%n
      CALL CH4_SingleInstance_ ( CH4_GridCompFinalize1_, i, &
                                gcCH4%gcs(i), w_c, impChem, expChem, &
                                nymd, nhms, cdt, status )
      VERIFY_(status)
   END DO

   DEALLOCATE ( gcCH4%gcs, stat=status )    
   gcCH4%n = -1

 END SUBROUTINE CH4_GridCompFinalize

 subroutine CH4_GridCompSetServices1_(  gc, chemReg, iname, rc)
 type(ESMF_GridComp), intent(INOUT) :: GC
 type(Chem_Registry), intent(INOUT) :: chemReg
 character(len=*),    intent(IN   ) :: iname
 integer,             intent(OUT  ) :: rc

 integer :: Status
 character(len=ESMF_MAXSTR) :: Iam

 Iam ="CH4_GridCOmpSetServices1_"

  call MAPL_AddImportSpec(GC,             &
       SHORT_NAME = 'CH4_industr'//iname, &
       LONG_NAME  = 'source species'  ,   &
       UNITS      = 'kg CH4 m-2 s-1',     &
       DIMS       = MAPL_DimsHorzOnly,    &
       VLOCATION  = MAPL_VLocationNone,   &
       RESTART    = MAPL_RestartSkip,     &
       RC         = STATUS)
  call MAPL_AddImportSpec(GC,             &
       SHORT_NAME = 'CH4_agwaste'//iname, &
       LONG_NAME  = 'source species'  ,   &
       UNITS      = 'kg CH4 m-2 s-1',     &
       DIMS       = MAPL_DimsHorzOnly,    &
       VLOCATION  = MAPL_VLocationNone,   &
       RESTART    = MAPL_RestartSkip,     &
       RC         = STATUS)
  call MAPL_AddImportSpec(GC,             &
       SHORT_NAME = 'CH4_bioburn'//iname, &
       LONG_NAME  = 'source species'  ,   &
       UNITS      = 'kg CH4 m-2 s-1',     &
       DIMS       = MAPL_DimsHorzOnly,    &
       VLOCATION  = MAPL_VLocationNone,   &
       RESTART    = MAPL_RestartSkip,     &
       RC         = STATUS)
  call MAPL_AddImportSpec(GC,             &
       SHORT_NAME = 'CH4_biofuel'//iname, &
       LONG_NAME  = 'source species'  ,   &
       UNITS      = 'kg CH4 m-2 s-1',     &
       DIMS       = MAPL_DimsHorzOnly,    &
       VLOCATION  = MAPL_VLocationNone,   &
       RESTART    = MAPL_RestartSkip,     &
       RC         = STATUS)
  call MAPL_AddImportSpec(GC,             &
       SHORT_NAME = 'CH4_minnatl'//iname, &
       LONG_NAME  = 'source species'  ,   &
       UNITS      = 'kg CH4 m-2 s-1',     &
       DIMS       = MAPL_DimsHorzOnly,    &
       VLOCATION  = MAPL_VLocationNone,   &
       RESTART    = MAPL_RestartSkip,     &
       RC         = STATUS)
  call MAPL_AddImportSpec(GC,             &
       SHORT_NAME = 'CH4_wetland'//iname, &
       LONG_NAME  = 'source species',     &
       UNITS      = 'kg CH4 m-2 s-1',     &
       DIMS       = MAPL_DimsHorzOnly,    &
       VLOCATION  = MAPL_VLocationNone,   &
       RESTART    = MAPL_RestartSkip,     &
       RC         = STATUS)
  call MAPL_AddImportSpec(GC,             &
       SHORT_NAME = 'CH4_oh'//iname,      &
       LONG_NAME  = 'source species',     &
       UNITS      = 'molecules cm-3',     &
       DIMS       = MAPL_DimsHorzVert,    &
       VLOCATION  = MAPL_VLocationCenter, &
       RESTART    = MAPL_RestartSkip,     &
       RC         = STATUS)
  call MAPL_AddImportSpec(GC,             &
       SHORT_NAME = 'CH4_cl'//iname,      &
       LONG_NAME  = 'source species',     &
       UNITS      = 'molecules cm-3',     &
       DIMS       = MAPL_DimsHorzVert,    &
       VLOCATION  = MAPL_VLocationCenter, &
       RESTART    = MAPL_RestartSkip,     &
       RC         = STATUS)
  call MAPL_AddImportSpec(GC,             &
       SHORT_NAME = 'CH4_o1d'//iname,     &
       LONG_NAME  = 'source species',     &
       UNITS      = 'molecules cm-3',     &
       DIMS       = MAPL_DimsHorzVert,    &
       VLOCATION  = MAPL_VLocationCenter, &
       RESTART    = MAPL_RestartSkip,     &
       RC         = STATUS)
  VERIFY_(STATUS)

 RETURN_(ESMF_SUCCESS)

 end subroutine CH4_GridCompSetServices1_

!--------------------------------------------------------------------------

!                      Single Instance Methods

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  CH4_GridCompInitialize --- Initialize CH4_GridComp
!
! !INTERFACE:
!

   SUBROUTINE CH4_GridCompInitialize1_ ( gcCH4, w_c, impChem, expChem, nymd, nhms, cdt, rc )

! !USES:

  IMPLICIT NONE

! !INPUT PARAMETERS:

   TYPE(Chem_Bundle), intent(in) :: w_c        ! Chemical tracer fields      
   INTEGER, INTENT(IN) :: nymd, nhms           ! time
   REAL,    INTENT(IN) :: cdt                  ! chemical timestep (secs)


! !OUTPUT PARAMETERS:

   TYPE(CH4_GridComp1), INTENT(INOUT) :: gcCH4  ! Grid Component
   TYPE(ESMF_State), INTENT(INOUT)  :: impChem  ! Import State
   TYPE(ESMF_State), INTENT(INOUT)  :: expChem  ! Export State
   INTEGER, INTENT(OUT) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 - 

! !DESCRIPTION: Initializes the CH4 Grid Component. It primarily sets
!               the import state for each active constituent package.
!
! !REVISION HISTORY:
!
!  18Sep2003 da Silva  First crack.
!  31May2005  Nielsen  Mods for 7 CH4 bins, 5 region masks
!  04Nov2005     Bian  CO tagged to 4 regions 
!                      (global, North America, South America, and Africa)
!                      for CR-AVE
!  25Oct2012  Nielsen  Added photolysis.
!
!EOP
!-------------------------------------------------------------------------

   CHARACTER(LEN=*), PARAMETER :: Iam = 'CH4_GridCompInitialize1_'

   CHARACTER(LEN=ESMF_MAXSTR) :: rcfilen

   INTEGER :: j, n, status
   INTEGER :: i1, i2, im, j1, j2, jm, km
   INTEGER :: nTimes, begTime, incSecs
   INTEGER :: nbeg, nend, nymd1, nhms1

   REAL :: limitN, limitS, radTODeg
   REAL, ALLOCATABLE :: var2D(:,:)
   real, pointer     :: ptr2d(:,:) => null()

   rcfilen = gcCH4%rcfilen
   gcCH4%name = 'GEOS-5/GOCART Parameterized CH4 Package'
   radTODeg = 57.2957795

!  Initialize local variables
!  --------------------------
   rc = 0
   i1 = w_c%grid%i1
   i2 = w_c%grid%i2
   im = w_c%grid%im

   j1 = w_c%grid%j1
   j2 = w_c%grid%j2
   jm = w_c%grid%jm

   km = w_c%grid%km

   nbeg  = w_c%reg%i_CH4
   nend  = w_c%reg%j_CH4

!  It requires 1 bin
!  -----------------
   if ( nbeg /= nend ) then
      IF(MAPL_AM_I_ROOT()) PRINT *,TRIM(Iam)//": Must have only 1 bin at the single instance level"
      status = 1
      VERIFY_(status)
   end if

!  Load resource file
!  ------------------
   CALL I90_loadf ( TRIM(rcfilen), status )
   VERIFY_(status)

!  Maximum allowed solar zenith angle for "daylight"
!  -------------------------------------------------
   CALL I90_label ( 'solar_ZA_cutoff:', status )
   VERIFY_(status)
   gcCH4%szaCutoff = I90_gfloat( status )
   VERIFY_(status)

!  Run-time debug switch
!  ---------------------
   CALL I90_label ( 'DEBUG:', status )
   VERIFY_(status)
   n = I90_gint ( status )
   VERIFY_(status)
   IF(n /= 0) THEN
    gcCH4%DebugIsOn = .TRUE.
   ELSE
    gcCH4%DebugIsOn = .FALSE.
   END IF

!  Methane photolysis feedback switch
!  ----------------------------------
   CALL I90_label ( 'CH4_Feedback:', status )
   VERIFY_(status)
   n = I90_gint ( status )
   VERIFY_(status)
   IF(n /= 0) THEN
    gcCH4%CH4FeedBack = .TRUE.
   ELSE
    gcCH4%CH4FeedBack = .FALSE.
   END IF

!  Water vapor feedback switch
!  ---------------------------
   CALL I90_label ( 'H2O_Feedback:', status )
   VERIFY_(status)
   n = I90_gint ( status )
   VERIFY_(status)
   IF(n /= 0) THEN
    gcCH4%H2OFeedBack = .TRUE.
   ELSE
    gcCH4%H2OFeedBack = .FALSE.
   END IF

!  Possibly oxidants are in a different unit
!  Allowable choices are: "mol/mol" or "mol mol-1"
!  or else behavior is as though input in 
!  "molecules cm-3"
!  Should this checking be done now in ExtData?
!  --------------------------------------------
   call i90_label('units_oh:', status)
   if (status /= 0) then
      gcCH4%units_oh = " "
   else
      call I90_gtoken(gcCH4%units_oh, status)
   end if
   VERIFY_(status)

   RETURN

 END SUBROUTINE CH4_GridCompInitialize1_

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  CH4_GridCompRun
!
! !INTERFACE:
!

   SUBROUTINE CH4_GridCompRun1_ ( gcCH4, w_c, impChem, expChem, nymd, nhms, cdt, rc )

! !USES:

  IMPLICIT NONE

! !INPUT/OUTPUT PARAMETERS:

   TYPE(CH4_GridComp1), INTENT(INOUT) :: gcCH4   ! Grid Component
   TYPE(Chem_Bundle), INTENT(INOUT) :: w_c       ! Chemical tracer fields   

! !INPUT PARAMETERS:

   TYPE(ESMF_State), INTENT(inout) :: impChem    ! Import State
   INTEGER, INTENT(IN) :: nymd, nhms             ! time
   REAL,    INTENT(IN) :: cdt                    ! chemical timestep (secs)


! !OUTPUT PARAMETERS:

   TYPE(ESMF_State), intent(inout) :: expChem    ! Export State
   INTEGER, INTENT(OUT) ::  rc                   ! Error return code:
                                                 !  0 - all is well
                                                 !  1 -
 
! !DESCRIPTION: This routine implements the CH4 Driver for GOCART.
!
! !REVISION HISTORY:
!
!  24 Jun 2010 Nielsen:  First crack.
!  25 Oct 2012 Nielsen:  Added photolysis.
!
!EOP
!-------------------------------------------------------------------------

   CHARACTER(LEN=*), PARAMETER :: Iam = 'CH4_GridCompRun1_'

!  Input fields from fvGCM
!  -----------------------
   REAL, POINTER, DIMENSION(:,:)   ::  cellArea => null()
   REAL, POINTER, DIMENSION(:,:)   ::  pblh     => null()
   REAL, POINTER, DIMENSION(:,:,:) ::  zle      => null()
   REAL, POINTER, DIMENSION(:,:,:) ::  T        => null()
   REAL, POINTER, DIMENSION(:,:,:) ::  Q        => null()
   REAL, POINTER, DIMENSION(:,:,:) ::  qtot     => null()
   REAL, POINTER, DIMENSION(:,:,:) ::  rhowet   => null()
   REAL, POINTER, DIMENSION(:,:,:) ::  ple      => null()

   INTEGER :: i1, i2, im, j1, j2, jm, km, idiag, iXj
   INTEGER :: i, j, k, kReverse, n, nbeg, nend
   INTEGER :: nymd1, nhms1, status

   REAL    :: qmin, qmax

   REAL, ALLOCATABLE :: cellDepth(:,:,:), cellVolume(:,:,:)
   REAL, ALLOCATABLE ::      rkoh(:,:,:),       rkcl(:,:,:),    rko1d(:,:,:)
   REAL, ALLOCATABLE ::         p(:,:,:),      ndwet(:,:,:)
   REAL, ALLOCATABLE ::    dCH4ox(:,:,:),      photJ(:,:,:), dCH4Phot(:,:,:)

   REAL, POINTER :: ptr2d(:,:)   => null()
   REAL, POINTER :: ptr3d(:,:,:) => null()

   REAL, ALLOCATABLE, DIMENSION(:,:) :: psdry, psch4

#define EXPORT   expChem
#define iNAME    TRIM(gcCH4%iname)

#define CH4EM    CH4_emis
#define CH4SC    CH4_surface
#define CH4CL    CH4_column
#define CH4DRY   CH4_dry
#define CH4SCDRY CH4_surfdry
#define CH4CLDRY CH4_coldry
#define CH4PD    CH4_prod
#define CH4LS    CH4_loss
#define CH4JL    CH4_phot
#define CH4QP    CH4_qprod

#include "CH4_GetPointer___.h"

!  Initialize local variables
!  --------------------------
   rc = 0
   i1 = w_c%grid%i1
   i2 = w_c%grid%i2
   im = w_c%grid%im

   j1 = w_c%grid%j1
   j2 = w_c%grid%j2
   jm = w_c%grid%jm
   km = w_c%grid%km

   iXj = ( i2 - i1 + 1 ) * ( j2 - j1 + 1 )

   nbeg  = w_c%reg%i_CH4
   nend  = w_c%reg%j_CH4

!  It requires 1 bin
!  -----------------
   if (nbeg /= nend) then
      if (MAPL_AM_I_ROOT()) then
         print *, trim(iNAME)//": Must have only 1 bin at the single instance level"
      endif
      status = 1
      VERIFY_(status)
   endif

!  Get imports
!  -----------
   CALL MAPL_GetPointer(impChem, PBLH,     'ZPBL',        RC=status)
   VERIFY_(status)
   CALL MAPL_GetPointer(impChem, ZLE,      'ZLE',         RC=status)
   VERIFY_(status)
   CALL MAPL_GetPointer(impChem, T,        'T',           RC=status)
   VERIFY_(status)
   CALL MAPL_GetPointer(impChem, Q,        'Q',           RC=status)
   VERIFY_(status)
   CALL MAPL_GetPointer(impChem, qtot,     'QTOT',        RC=status)
   VERIFY_(status)
   CALL MAPL_GetPointer(impChem, rhowet,   'AIRDENS',     RC=status)
   VERIFY_(status)
   CALL MAPL_GetPointer(impChem, cellArea, 'AREA',        RC=status)
   VERIFY_(status)
   CALL MAPL_GetPointer(impChem, ple,      'PLE',         RC=status)
   VERIFY_(status)

   if (gcCH4%DebugIsOn) then
      call pmaxmin('CH4:AREA', cellArea, qmin, qmax, iXj,  1,   1. )
      call pmaxmin('CH4:ZPBL',     pblh, qmin, qmax, iXj,  1,   1. )
      call pmaxmin('CH4:ZLE',       zle, qmin, qmax, iXj, km+1, 1. )
      call pmaxmin('CH4:T',           T, qmin, qmax, iXj, km,   1. )
      call pmaxmin('CH4:Q',           q, qmin, qmax, iXj, km,   1. )
      call pmaxmin('CH4:QTOT',     qtot, qmin, qmax, iXj, km,   1. )
      call pmaxmin('CH4:RHOWET', rhowet, qmin, qmax, iXj, km,   1. )
      call pmaxmin('CH4:PLE',       ple, qmin, qmax, iXj, km+1, 1. )
   endif

!  Allocate memory, etc
!  --------------------
   ALLOCATE ( gcCH4%eCH4_ind(i1:i2,j1:j2),   &
              gcCH4%eCH4_wetl(i1:i2,j1:j2),  &
              gcCH4%eCH4_agw(i1:i2,j1:j2),   &
              gcCH4%eCH4_bb(i1:i2,j1:j2),    &
              gcCH4%eCH4_bb_(i1:i2,j1:j2),   &
              gcCH4%eCH4_bf(i1:i2,j1:j2),    &
              gcCH4%eCH4_mnat(i1:i2,j1:j2),  &
              gcCH4%OHnd(i1:i2,j1:j2,km),    &
              gcCH4%Clnd(i1:i2,j1:j2,km),    &
              gcCH4%O1Dnd(i1:i2,j1:j2,km), STAT=status )
   VERIFY_(status)

!  Update CH4 emissions and OH, Cl, and O(1D) number densities
!  (in molec cm^-3) once each day
!  ------------------------------------------------------------

!  Industrial Sources
!  ------------------
   call MAPL_GetPointer(impChem, ptr2d, 'CH4_industr'//iNAME,rc=status)
   VERIFY_(STATUS)
   gcCH4%eCH4_ind = ptr2d

!  Ag/Waste Sources
!  ----------------
   call MAPL_GetPointer(impChem, ptr2d, 'CH4_agwaste'//iNAME,rc=status)
   VERIFY_(STATUS)
   gcCH4%eCH4_agw = ptr2d

!  Biomass Burning
!  ---------------
   call MAPL_GetPointer(impChem, ptr2d, 'CH4_bioburn'//iNAME,rc=status)
   VERIFY_(STATUS)
   gcCH4%eCH4_bb = ptr2d

!  Biofuel
!  -------
   call MAPL_GetPointer(impChem, ptr2d, 'CH4_biofuel'//iNAME,rc=status)
   VERIFY_(STATUS)
   gcCH4%eCH4_bf = ptr2d

!  Wetlands
!  --------
   call MAPL_GetPointer(impChem, ptr2d, 'CH4_wetland'//iNAME,rc=status)
   VERIFY_(STATUS)
   gcCH4%eCH4_wetl = ptr2d

!  Minor Natural (currently termites, soil absorption)
!  ---------------------------------------------------
   call MAPL_GetPointer(impChem, ptr2d, 'CH4_minnatl'//iNAME,rc=status)
   VERIFY_(STATUS)
   gcCH4%eCH4_mnat = ptr2d

!  Background OH, Cl, and O(1D) for loss term
!  ------------------------------------------
   call MAPL_GetPointer(impChem,ptr3d,'CH4_oh'//trim(iNAME),rc=status)
   VERIFY_(STATUS)
   gcCH4%OHnd=ptr3d

   call MAPL_GetPointer(impChem,ptr3d,'CH4_cl'//trim(iNAME),rc=status)
   VERIFY_(STATUS)
   gcCH4%Clnd=ptr3d

   call MAPL_GetPointer(impChem,ptr3d,'CH4_o1d'//trim(iNAME),rc=status)
   VERIFY_(STATUS)
   gcCH4%O1Dnd=ptr3d

!  Allocate temporary workspace
!  ----------------------------
   allocate( p(i1:i2,j1:j2,km), &
         ndwet(i1:i2,j1:j2,km), &
          rkoh(i1:i2,j1:j2,km), &
          rkcl(i1:i2,j1:j2,km), &
         rko1d(i1:i2,j1:j2,km), &
        dCH4ox(i1:i2,j1:j2,km), &
     cellDepth(i1:i2,j1:j2,km), &
    cellVolume(i1:i2,j1:j2,km), STAT=status)
   VERIFY_(status)

!  Layer mean pressures
!  --------------------
   do k=1,km
      p(i1:i2,j1:j2,k) = (ple(i1:i2,j1:j2,k-1)+ple(i1:i2,j1:j2,k))*0.50
   enddo

!  Wet-air number density
!  ----------------------
   ndwet(i1:i2,j1:j2,1:km) = rhowet(i1:i2,j1:j2,1:km)*MAPL_AVOGAD/MAPL_AIRMW
 
!  Handle mole fraction or number density units of oxidants
!  --------------------------------------------------------
   if (trim(gcCH4%units_oh) == 'mol/mol' .or. trim(gcCH4%units_oh) == 'mol mol-1') then
       gcCH4%OHnd(i1:i2,j1:j2,1:km)  = gcCH4%OHnd(i1:i2,j1:j2,1:km)  &
                                     * ndwet(i1:i2,j1:j2,1:km)
       gcCH4%Clnd(i1:i2,j1:j2,1:km)  = gcCH4%Clnd(i1:i2,j1:j2,1:km)  &
                                     * ndwet(i1:i2,j1:j2,1:km)
       gcCH4%O1Dnd(i1:i2,j1:j2,1:km) = gcCH4%O1Dnd(i1:i2,j1:j2,1:km) &
                                     * ndwet(i1:i2,j1:j2,1:km)
   else
!      Otherwise, assume units are molec cm^-3 and convert to molec m^-3
       gcCH4%OHnd(i1:i2,j1:j2,1:km)  =  gcCH4%OHnd(i1:i2,j1:j2,1:km)*1.00E+06
       gcCH4%Clnd(i1:i2,j1:j2,1:km)  =  gcCH4%Clnd(i1:i2,j1:j2,1:km)*1.00E+06
       gcCH4%O1Dnd(i1:i2,j1:j2,1:km) = gcCH4%O1Dnd(i1:i2,j1:j2,1:km)*1.00E+06
   end if

!  Cell depth and volume
!  ---------------------
   do k=1,km
      cellDepth(i1:i2,j1:j2,k)  = (ple(i1:i2,j1:j2,k)-ple(i1:i2,j1:j2,k-1)) &
                                / (rhowet(i1:i2,j1:j2,k)*MAPL_GRAV)
      cellVolume(i1:i2,j1:j2,k) = cellArea(i1:i2,j1:j2)*cellDepth(i1:i2,j1:j2,k)
   enddo

   if (gcCH4%DebugIsOn) then
      call pmaxmin('CH4: eCH4_ind',   gcCH4%eCH4_ind, qmin, qmax, iXj, 1,  1. )
      call pmaxmin('CH4: eCH4_agw',   gcCH4%eCH4_agw, qmin, qmax, iXj, 1,  1. )
      call pmaxmin('CH4: eCH4_bb',     gcCH4%eCH4_bb, qmin, qmax, iXj, 1,  1. )
      call pmaxmin('CH4: eCH4_bf',     gcCH4%eCH4_bf, qmin, qmax, iXj, 1,  1. )
      call pmaxmin('CH4: eCH4_wetl', gcCH4%eCH4_wetl, qmin, qmax, iXj, 1,  1. )
      call pmaxmin('CH4: eCH4_mnat', gcCH4%eCH4_mnat, qmin, qmax, iXj, 1,  1. )
      call pmaxmin('CH4: OH Conc',        gcCH4%OHnd, qmin, qmax, iXj, 1,  1. )
      call pmaxmin('CH4: Cl Conc',        gcCH4%Clnd, qmin, qmax, iXj, 1,  1. )
      call pmaxmin('CH4: O(1D) Conc',    gcCH4%O1Dnd, qmin, qmax, iXj, 1,  1. )
      call pmaxmin('CH4: Cell Depth',      cellDepth, qmin, qmax, iXj, km, 1. )
      call pmaxmin('CH4: Cell Vol',       cellVolume, qmin, qmax, iXj, km, 1. )
      call pmaxmin('CH4: Wet-air ND',          ndwet, qmin, qmax, iXj, km, 1. )
   endif

!  Loss rate [m^3 s^-1] for OH + CH4 => CO + products
!  --------------------------------------------------
   rkoh(i1:i2,j1:j2,1:km)  = 2.45E-12*1.00E-06*exp(-1775./T(i1:i2,j1:j2,1:km))

!  Loss rate [m^3 s^-1] for Cl + CH4 => CO + products
!  --------------------------------------------------
   rkcl(i1:i2,j1:j2,1:km)  = 7.10E-12*1.00E-06*exp(-1270./T(i1:i2,j1:j2,1:km))

!  Loss rate [m^3 s^-1] for O(1D) + CH4 => CO + products
!  -----------------------------------------------------
   rko1d(i1:i2,j1:j2,1:km) = 1.75E-10*1.00E-06

!  Change in CH4 mole fraction due to oxidation
!  --------------------------------------------
   dCH4ox(i1:i2,j1:j2,1:km) = w_c%qa(nbeg)%data3d(i1:i2,j1:j2,1:km)                         &
                              * (    rkoh(i1:i2,j1:j2,1:km)*gcCH4%OHnd(i1:i2,j1:j2,1:km)    &
                                  +  rkcl(i1:i2,j1:j2,1:km)*gcCH4%Clnd(i1:i2,j1:j2,1:km)    &
                                  + rko1d(i1:i2,j1:j2,1:km)*gcCH4%O1Dnd(i1:i2,j1:j2,1:km) )

!  Vertically integrated CH4 loss due to oxidation (only)
!  ------------------------------------------------------
   if (associated(CH4_loss)) then
      CH4_loss(i1:i2,j1:j2) = 0.
      do k = 1,km
         CH4_loss(i1:i2,j1:j2) = CH4_loss(i1:i2,j1:j2)                       &
                               +     dCH4ox(i1:i2,j1:j2,k)*mwtCH4/MAPL_AIRMW &
                                 * w_c%delp(i1:i2,j1:j2,k)/MAPL_GRAV
      enddo
   endif

!  CH4 production (none)
!  ---------------------
   if (associated(CH4_prod)) CH4_prod(i1:i2,j1:j2) = 0.

!  Calculate photolytic loss rates, J [s^-1], for CH4 + hv => 2H2O + CO
!  Notice that J and the losses are always computed. However, the setting 
!  of the feedback switch(es) determines if the increments are actually applied
!  ----------------------------------------------------------------------------
   ALLOCATE(photJ(i1:i2,j1:j2,1:km), dCH4Phot(i1:i2,j1:j2,1:km), STAT=status)
   VERIFY_(STATUS)

   photJ(i1:i2,j1:j2,1:km) = 0.

   CALL getJRates(status)
   VERIFY_(status)

!  Change in CH4 due to photolysis
!  -------------------------------
   dCH4Phot(i1:i2,j1:j2,1:km) = photJ(i1:i2,j1:j2,1:km) * w_c%qa(nbeg)%data3d(i1:i2,j1:j2,1:km)

   if (associated(CH4_phot)) THEN
      CH4_phot(i1:i2,j1:j2,1:km) = dCH4Phot(i1:i2,j1:j2,1:km)
   endif

!  Increment the CH4 mole fraction due to oxidation 
!  ------------------------------------------------
   w_c%qa(nbeg)%data3d(i1:i2,j1:j2,1:km) = w_c%qa(nbeg)%data3d(i1:i2,j1:j2,1:km) - cdt*dCH4ox(i1:i2,j1:j2,1:km)

!  Increment the CH4 mole fraction due to photolysis when the switch is on
!  -----------------------------------------------------------------------
   if (gcCH4%CH4FeedBack) then
      w_c%qa(nbeg)%data3d(i1:i2,j1:j2,1:km) = w_c%qa(nbeg)%data3d(i1:i2,j1:j2,1:km) - cdt*dCH4Phot(i1:i2,j1:j2,1:km)
   endif

!  If both feedback switches are on, increment water vapor by
!  adding two molecules of H2O for each CH4 molecule lost
!  ----------------------------------------------------------
   if (associated(CH4_qprod)) CH4_qprod(i1:i2,j1:j2,1:km) = 0.

   if (gcCH4%CH4FeedBack .and. gcCH4%H2OFeedBack) then
      Q(i1:i2,j1:j2,1:km) = Q(i1:i2,j1:j2,1:km) + 2.00*cdt*dCH4Phot(i1:i2,j1:j2,1:km)*MAPL_H2OMW/MAPL_AIRMW

!     Water vapor tendency [kg kg^-1 s^-1]
!     ------------------------------------
      if (associated(CH4_qprod)) then
         CH4_qprod(i1:i2,j1:j2,1:km) = 2.00*dCH4Phot(i1:i2,j1:j2,1:km)*MAPL_H2OMW/MAPL_AIRMW
      endif
   endif

!  Compute and add surface emissions
!  ---------------------------------
   call CH4_Emission(rc)

!  Surface concentration [ppbv]
!  ----------------------------
   if (associated(CH4_surface)) then
      CH4_surface(i1:i2,j1:j2) = w_c%qa(nbeg)%data3d(i1:i2,j1:j2,km)*1.00E+09
   endif

!  Column burden [kg m-2]
!  ----------------------
   if (associated(CH4_column)) then
      CH4_column(i1:i2,j1:j2) = 0.
      do k = 1,km
        CH4_column(i1:i2,j1:j2) = CH4_column(i1:i2,j1:j2)                              &
                                + w_c%qa(nbeg)%data3d(i1:i2,j1:j2,k)*mwtCH4/MAPL_AIRMW &
                                           * w_c%delp(i1:i2,j1:j2,k)/MAPL_GRAV
      enddo
   endif

!  Dry-air mole fraction [mol mol-1]
!  ---------------------------------
   if (associated(CH4_dry)) then
      CH4_dry(i1:i2,j1:j2,1:km) = w_c%qa(nbeg)%data3d(i1:i2,j1:j2,1:km) &
                                         / (1. - qtot(i1:i2,j1:j2,1:km))
   endif

!  Dry-air surface concentration [mol mol-1]
!  -----------------------------------------
   if (associated(CH4_surfdry)) then
      CH4_surfdry(i1:i2,j1:j2) = w_c%qa(nbeg)%data3d(i1:i2,j1:j2,km)*(1. - qtot(i1:i2,j1:j2,km))
   endif

!  Dry-air column average [mol mol-1]
!  ----------------------------------
   allocate(psdry(i1:i2,j1:j2), stat=status)    ! dry-air surface pressure
   allocate(psch4(i1:i2,j1:j2), stat=status)    ! ch4     surface pressure

   if (associated(CH4_coldry)) then
      psdry(i1:i2,j1:j2) = 0.
      psch4(i1:i2,j1:j2) = 0.
      do k = 1,km
         psdry(i1:i2,j1:j2) = psdry(i1:i2,j1:j2) + w_c%delp(i1:i2,j1:j2,k)*(1. - qtot(i1:i2,j1:j2,k))
         psch4(i1:i2,j1:j2) = psch4(i1:i2,j1:j2) + w_c%delp(i1:i2,j1:j2,k)*w_c%qa(nbeg)%data3d(i1:i2,j1:j2,k)
      enddo
      CH4_coldry(i1:i2,j1:j2) = psch4(i1:i2,j1:j2)/psdry(i1:i2,j1:j2)
   endif

   deallocate(psdry, psch4)

   IF(gcCH4%DebugIsOn) THEN
     IF(ASSOCIATED(CH4_emis))    CALL pmaxmin('CH4: emis',     CH4_emis,    qmin, qmax, iXj,  1, 1. )
     IF(ASSOCIATED(CH4_loss))    CALL pmaxmin('CH4: loss',     CH4_loss,    qmin, qmax, iXj,  1, 1. )
     IF(ASSOCIATED(CH4_prod))    CALL pmaxmin('CH4: prod',     CH4_prod,    qmin, qmax, iXj,  1, 1. )
     IF(ASSOCIATED(CH4_column))  CALL pmaxmin('CH4: column',   CH4_column,  qmin, qmax, iXj,  1, 1. )
     IF(ASSOCIATED(CH4_surface)) CALL pmaxmin('CH4: surface',  CH4_surface, qmin, qmax, iXj,  1, 1. )
     IF(ASSOCIATED(CH4_qprod))   CALL pmaxmin('CH4: qprod',    CH4_qprod,   qmin, qmax, iXj, km, 1. )
     IF(ASSOCIATED(CH4_phot))    CALL pmaxmin('CH4: dch4phot', dCH4phot,    qmin, qmax, iXj, km, 1. )
     IF(ASSOCIATED(CH4_dry))     CALL pmaxmin('CH4: dry',      CH4_dry,     qmin, qmax, iXj, km, 1. )
   END IF

!  Housekeeping
!  ------------
   DEALLOCATE(gcCH4%eCH4_ind,  gcCH4%eCH4_wetl,  gcCH4%eCH4_agw, &
              gcCH4%eCH4_bb,   gcCH4%eCH4_bb_,   gcCH4%eCH4_bf,  &
              gcCH4%eCH4_mnat, gcCH4%OHnd,       gcCH4%Clnd,     &
              gcCH4%O1Dnd,     STAT=status)
   VERIFY_(status)
   DEALLOCATE(p, ndwet, rkoh, rkcl, rko1d, dCH4ox, cellDepth, cellVolume, STAT=status)
   VERIFY_(status)
   DEALLOCATE(photJ, dCH4Phot, STAT=status)
   VERIFY_(status)

   RETURN

 CONTAINS
 
!-------------------------------------------------------------------------
  SUBROUTINE CH4_Emission ( rc )

  IMPLICIT NONE

! OUTPUT VARIABLES
  INTEGER, INTENT(OUT) :: rc ! Error return code

! LOCAL VARIABLES
  CHARACTER(LEN=*), PARAMETER :: myname = 'CH4_Emission'

  INTEGER ::  i, j, k, kt, minkPBL
  INTEGER, ALLOCATABLE :: index(:)

  REAL, ALLOCATABLE :: pblLayer(:,:), sfcFlux(:,:)
  REAL, ALLOCATABLE :: fPBL(:,:,:)

  rc = 0

  allocate(sfcFlux(i1:i2,j1:j2),      STAT=rc)      ! emissions
  allocate(   fPBL(i1:i2,j1:j2,1:km), STAT=rc)      ! partitioning of BB

! Apply biomass burning diurnal cycle if desired
! ----------------------------------------------
  if (w_c%diurnal_bb) then
     gcCH4%eCH4_bb_(:,:) = gcCH4%eCH4_bb(:,:)

     call Chem_BiomassDiurnal ( gcCH4%eCH4_bb, gcCH4%eCH4_bb_,   &
                                w_c%grid%lon(:,:)*radToDeg,      &
                                w_c%grid%lat(:,:)*radToDeg, nhms, cdt )      
  endif

! Find the layer that contains the PBL
! Layer thicknesses are ZLE(:,:,0:km)
! ------------------------------------
  allocate(index(0:km), STAT=rc)
  allocate(pblLayer(i1:i2,j1:j2), STAT=rc)
  do j = j1,j2
     do i = i1,i2
        index(0:km) = 0
        where(zle(i,j,0:km) - zle(i,j,km) > pblh(i,j)) index(0:km) = 1
        pblLayer(i,j) = sum(index)
     enddo
  enddo
  minkPBL = minval(pblLayer)

! Determine partitioning fraction based on layer thicknesses
! ----------------------------------------------------------
  fPBL(i1:i2,j1:j2,1:km) = 0.
  do j = j1,j2
     do i = i1,i2
        kt = pblLayer(i,j)
        do k = kt,km
           fPBL(i,j,k) = (zle(i,j, k-1) - zle(i,j,k))  &
                       / (zle(i,j,kt-1) - zle(i,j,km))
        enddo
     enddo
  enddo

  deallocate(index, pblLayer, STAT=rc)

! Establish range of layers on which to work
! ------------------------------------------
  kt = minkPBL

  Layer: do k = kt,km

!    Emissions: Weighted biomass burning
!    -----------------------------------
     sfcFlux(i1:i2,j1:j2) = gcCH4%eCH4_BB(i1:i2,j1:j2)*fPBL(i1:i2,j1:j2,k)

!    Add other emission components when in surface layer
!    ---------------------------------------------------
     if (k == km) sfcFlux(i1:i2,j1:j2) = sfcFlux(i1:i2,j1:j2)         &
                                       + gcCH4%eCH4_bf(i1:i2,j1:j2)   &
                                       + gcCH4%eCH4_ind(i1:i2,j1:j2)  &
                                       + gcCH4%eCH4_agw(i1:i2,j1:j2)  &
                                       + gcCH4%eCH4_mnat(i1:i2,j1:j2) &
                                       + gcCH4%eCH4_wetl(i1:i2,j1:j2)

!    Update CH4 at this level
!    ------------------------
     w_c%qa(nbeg)%data3d(i1:i2,j1:j2,k) = w_c%qa(nbeg)%data3d(i1:i2,j1:j2,k)           &
                                        + cdt * sfcFlux(i1:i2,j1:j2)*MAPL_AIRMW/mwtCH4 &
                                              / (w_c%delp(i1:i2,j1:j2,k)/MAPL_GRAV)
  enddo Layer

! Update Surface flux diagnostic for this bin
! -------------------------------------------
  if (associated(CH4_emis)) then
     CH4_emis(i1:i2,j1:j2) = gcCH4%eCH4_bb(i1:i2,j1:j2)   + gcCH4%eCH4_bf(i1:i2,j1:j2)   &
                           + gcCH4%eCH4_ind(i1:i2,j1:j2)  + gcCH4%eCH4_agw(i1:i2,j1:j2)  &
                           + gcCH4%eCH4_mnat(i1:i2,j1:j2) + gcCH4%eCH4_wetl(i1:i2,j1:j2)
  endif

  deallocate(fPBL, sfcFlux, STAT=rc)

  return

  end subroutine CH4_Emission

!-------------------------------------------------------------------------
! Borrowed from meso_phot.F of StratChem, where number densities are cgs [cm^{-3}]
  SUBROUTINE getJRates(rc)

  IMPLICIT NONE

  INTEGER, INTENT(out) :: rc

  REAL, ALLOCATABLE :: o2Column(:,:,:)
  REAL, ALLOCATABLE :: SZARad(:,:)
  REAL, ALLOCATABLE :: SZADeg(:,:)
  REAL, ALLOCATABLE :: sinSZA(:,:)
  REAL, ALLOCATABLE :: zgrz(:,:)
  REAL, ALLOCATABLE :: sfaca(:,:)
  REAL, ALLOCATABLE :: arg(:,:)

  REAL, PARAMETER :: wavel = 1215.7
  REAL, PARAMETER :: O2xs  = 1.000E-20
  REAL, PARAMETER :: CH4xs = 2.000E-17
  REAL, PARAMETER :: sflux = 4.006E+11

! Constants for Chapman function at high solar zenith angle
! ---------------------------------------------------------
  REAL, PARAMETER :: hbar = 6.79
  REAL, PARAMETER :: zbar = 30.0
  REAL, PARAMETER :: r0   = 6.371E+03
  REAL, PARAMETER :: zp   = 60.0

  REAL, PARAMETER :: d1 = 1.060693
  REAL, PARAMETER :: d2 = 0.55643831
  REAL, PARAMETER :: d3 = 1.0619896
  REAL, PARAMETER :: d4 = 1.7245609
  REAL, PARAMETER :: d5 = 0.56498823
  REAL, PARAMETER :: d6 = 0.06651874

  REAL, PARAMETER :: O2Abv80km = 7.072926E+19 ![cm^{-2}]
  REAL, PARAMETER :: O2VMR = 0.20946

  REAL :: b
  REAL :: f
  REAL :: r
  REAL :: s

  INTEGER :: status

  CHARACTER(LEN=ESMF_MAXSTR) :: Iam = "CH4::getJRates"

  rc = 0
  b = SQRT(0.50*r0/hbar)

! O2 overhead number density profile [cm^{-2}]
! --------------------------------------------
  ALLOCATE(O2Column(i1:i2,j1:j2,1:km), STAT=status)
  VERIFY_(status)

  f = O2VMR*5.00E-05
  O2Column(:,:,1) = O2Abv80km+cellDepth(:,:,1)*ndwet(:,:,1)*f

  DO k = 2,km
   O2Column(:,:,k) = O2Column(:,:,k-1)+(cellDepth(:,:,k-1)*ndwet(:,:,k-1)+ &
                                        cellDepth(:,:,  k)*ndwet(:,:,  k))*f
  END DO

  IF(gcCH4%DebugIsOn) THEN
   CALL pmaxmin('CH4: O2Column', O2Column, qmin, qmax, iXj, km,  1. )
  END IF

! Grab some memory
! ----------------
  ALLOCATE(SZARad(i1:i2,j1:j2), STAT=status)
  VERIFY_(status)
  ALLOCATE(SZADeg(i1:i2,j1:j2), STAT=status)
  VERIFY_(status)
  ALLOCATE(sinSZA(i1:i2,j1:j2), STAT=status)
  VERIFY_(status)

  WHERE(w_c%cosz(i1:i2,j1:j2) > 1.00)
    SZARad(i1:i2,j1:j2) = 0.00
  ELSEWHERE
    SZARad(i1:i2,j1:j2) = ACOS(w_c%cosz(i1:i2,j1:j2))
  ENDWHERE
  SZADeg(i1:i2,j1:j2) = SZARad(i1:i2,j1:j2)*radToDeg
  sinSZA(i1:i2,j1:j2) = SIN(SZARad(i1:i2,j1:j2))

  ALLOCATE(zgrz(i1:i2,j1:j2), STAT=status)
  VERIFY_(status)

  WHERE(SZADeg(i1:i2,j1:j2) <= 90.00)
    zgrz(i1:i2,j1:j2) = 1000.00
  ELSEWHERE
    zgrz(i1:i2,j1:j2) = sinSZA(i1:i2,j1:j2)*(zp+r0)-r0
  ENDWHERE

  IF(gcCH4%DebugIsOn) THEN
   CALL pmaxmin('CH4: zgrz', zgrz, qmin, qmax, iXj, 1,  1. )
   CALL pmaxmin('CH4: cosz', w_c%cosz, qmin, qmax, iXj, 1,  1. )
  END IF

  ALLOCATE(sfaca(i1:i2,j1:j2), STAT=status)
  VERIFY_(status)
  sfaca(i1:i2,j1:j2) = 0.00

! Chapman function calculation from ACDB 2-D model
! ------------------------------------------------
  DO j = j1,j2
   DO i = i1,i2

    Daytime: IF(SZADeg(i,j) < gcCH4%szaCutoff) THEN

     IF(SZADeg(i,j) < 70.00) THEN

      sfaca(i,j) = 1.00/w_c%cosz(i,j)

     ELSE IF(zgrz(i,j) > 0.00) THEN

      s = b*ABS(w_c%cosz(i,j))

      IF(s <= 8.00) THEN
       s = (d1+d2*s)/(d3+d4*s+s**2)
      ELSE
       s = d5/(d6+s)
      END IF

      r = b*SQRT(MAPL_PI)
      sfaca(i,j) = r*s

      IF(SZADeg(i,j) > 90.00) THEN
       sfaca(i,j) = 2.00*r*EXP((r0+zbar)*(1.00-sinSZA(i,j))/hbar)-sfaca(i,j)
      END IF

     END IF

    END IF Daytime

   END DO
  END DO

  IF(gcCH4%DebugIsOn) THEN
   CALL pmaxmin('CH4: sfaca', sfaca, qmin, qmax, iXj, 1,  1. )
  END IF

  ALLOCATE(arg(i1:i2,j1:j2), STAT=status)
  VERIFY_(status)

! At each layer, compute the rate constant, J [s^{-1}], if the sun is up
! ----------------------------------------------------------------------
  DO k = 1,km

   WHERE(SZADeg(i1:i2,j1:j2) < gcCH4%szaCutoff)
    arg(i1:i2,j1:j2) = O2Column(i1:i2,j1:j2,k)*O2xs*sfaca(i1:i2,j1:j2)
    photJ(i1:i2,j1:j2,k) = sflux*EXP(-arg(i1:i2,j1:j2))*CH4xs
   END WHERE

  END DO

  IF(gcCH4%DebugIsOn) THEN
   CALL pmaxmin('CH4: photJ', photJ, qmin, qmax, iXj, km,  1. )
  END IF

  DEALLOCATE(SZARad, STAT=status)
  VERIFY_(status)
  DEALLOCATE(SZADeg, STAT=status)
  VERIFY_(status)
  DEALLOCATE(sinSZA, STAT=status)
  VERIFY_(status)
  DEALLOCATE(zgrz, STAT=status)
  VERIFY_(status)
  DEALLOCATE(sfaca, STAT=status)
  VERIFY_(status)
  DEALLOCATE(arg, STAT=status)
  VERIFY_(status)
  DEALLOCATE(O2Column, STAT=status)
  VERIFY_(status)

  RETURN

  END SUBROUTINE getJRates

 END SUBROUTINE CH4_GridCompRun1_

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  CH4_GridCompFinalize --- The Chem Driver 
!
! !INTERFACE:
!

   SUBROUTINE CH4_GridCompFinalize1_ ( gcCH4, w_c, impChem, expChem, nymd, nhms, cdt, rc )

! !USES:

  IMPLICIT NONE

! !INPUT/OUTPUT PARAMETERS:

   TYPE(CH4_GridComp1), INTENT(INOUT) :: gcCH4   ! Grid Component

! !INPUT PARAMETERS:

   TYPE(Chem_Bundle), INTENT(IN)  :: w_c      ! Chemical tracer fields   
   INTEGER, INTENT(IN) :: nymd, nhms          ! time
   REAL,    INTENT(IN) :: cdt                 ! chemical timestep (secs)


! !OUTPUT PARAMETERS:

   TYPE(ESMF_State), INTENT(INOUT) :: impChem   ! Import State
   TYPE(ESMF_State), INTENT(INOUT) :: expChem   ! Import State
   INTEGER, INTENT(OUT) ::  rc                  ! Error return code:
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

   CHARACTER(LEN=*), PARAMETER :: Iam = 'CH4_GridCompFinalize1_'
   rc = 0

   RETURN

 END SUBROUTINE CH4_GridCompFinalize1_

 END MODULE CH4_GridCompMod

!-----------------------------------------------------------------------

!                     Single Instance Wrapper

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  CH4_SingleInstance_ --- Runs single instance of method
!
! !INTERFACE:
!
  SUBROUTINE CH4_SingleInstance_ ( Method_, instance, &
                                  gcCH4, w_c, impChem, expChem, &
                                  nymd, nhms, cdt, rc )

! !USES:

  use CH4_GridCompMod
  use ESMF
  use MAPL
  use Chem_Mod 

  IMPLICIT NONE

! !INPUT PARAMETERS:

!  Input "function pointer"
!  -----------------------
   interface 
     subroutine Method_ (gc, w, imp, exp, ymd, hms, dt, rcode )
       use CH4_GridCompMod
       use ESMF
       use MAPL
       use Chem_Mod 
       type(CH4_GridComp1), intent(inout)  :: gc
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

   TYPE(CH4_GridComp1), INTENT(INOUT) :: gcCH4  ! Grid Component
   TYPE(ESMF_State), INTENT(INOUT)  :: impChem  ! Import State
   TYPE(ESMF_State), INTENT(INOUT)  :: expChem  ! Export State
   INTEGER, INTENT(OUT) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 - 

! !DESCRIPTION: Finalizes the CH4 Grid Component. Multiple instance
!               version.
!
! !REVISION HISTORY:
!
!  27Feb2008  da Silva  Introduced multiple instances
!
!EOP
!-------------------------------------------------------------------------

  integer n_CH4, i_CH4, j_CH4

! Save overall CH4 indices
! -----------------------
  n_CH4 = w_c%reg%n_CH4
  i_CH4 = w_c%reg%i_CH4
  j_CH4 = w_c%reg%j_CH4
  
! Customize indices for this particular instance
! ----------------------------------------------
  w_c%reg%n_CH4 = 1
  w_c%reg%i_CH4 = i_CH4 + instance - 1
  w_c%reg%j_CH4 = i_CH4 + instance - 1
  
! Execute the instance method
! ---------------------------
  call Method_ ( gcCH4, w_c, impChem, expChem, &
                 nymd, nhms, cdt, rc )

! Restore the overall CH4 indices
! ------------------------------
  w_c%reg%n_CH4 = n_CH4
  w_c%reg%i_CH4 = i_CH4
  w_c%reg%j_CH4 = j_CH4

  END SUBROUTINE CH4_SingleInstance_

!-----------------------------------------------------------------------
