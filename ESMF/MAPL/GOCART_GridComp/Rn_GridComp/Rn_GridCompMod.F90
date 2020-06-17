#include "MAPL_Generic.h"

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1     !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  Rn_GridCompMod --- Rn Grid Component Class
!
! !INTERFACE:
!

   MODULE  Rn_GridCompMod

! !USES:

   USE ESMF
   USE MAPL

   USE Chem_Mod 	     ! Chemistry Base Class
   USE Chem_StateMod	     ! Chemistry State
   USE Chem_ConstMod, only: grav
   USE Chem_UtilMod	     ! I/O

   USE m_inpak90	     ! Resource file management
   USE m_die, ONLY: die

   IMPLICIT NONE

! !PUBLIC TYPES:
!
   PRIVATE
   PUBLIC  Rn_GridComp       ! Multiple instance Radon object 
   PUBLIC  Rn_GridComp1      ! Single instance Radon object

!
! !PUBLIC MEMBER FUNCTIONS:
!

   PUBLIC  Rn_GridCompSetServices
   PUBLIC  Rn_GridCompInitialize
   PUBLIC  Rn_GridCompRun
   PUBLIC  Rn_GridCompFinalize

!
! !DESCRIPTION:
!
!  This module implements the Rn Grid Component. 
!
! !REVISION HISTORY:
!
!  16Sep2003 da Silva  First crack.
!  01Aug2006 da Silva  Extensions for GEOS-5.
!  10Mar2008 da Silva  Multiple instances for ARCTAS.
!  12Apr2008 Nielsen   Configured for radon.
!
!EOP
!-------------------------------------------------------------------------

  TYPE Rn_GridComp1

        CHARACTER(LEN=ESMF_MAXSTR) :: name            ! generic name of the package
        CHARACTER(LEN=ESMF_MAXSTR) :: iname           ! instance name
        CHARACTER(LEN=ESMF_MAXSTR) :: rcfilen         ! resource file name
        CHARACTER(LEN=ESMF_MAXSTR) :: regionsString   ! Comma-delimited string of regions

        INTEGER :: instance                   ! instance number
        INTEGER :: BCnymd                     ! Date of last emissions update

        REAL :: halfLife                      ! Half-life
        CHARACTER(LEN=ESMF_MAXSTR) :: halfLifeUnit    ! Half-life unit: years, days, or seconds
        REAL :: decayConstant                 ! Decay constant, inverse seconds.
	REAL :: emission                      ! kg m^{-2} s^{-1}

        REAL, POINTER :: regionMask(:,:)      ! regional mask
        REAL, POINTER :: RnsfcFlux(:,:)       ! Rn surface flux kg m^-2 s^-1
        REAL, POINTER :: ScheryEmission(:,:)  ! Monthly mean emission mBq m^{-2} s^{-1}

        LOGICAL :: DebugIsOn                  ! Run-time debug switch

  END TYPE Rn_GridComp1

  TYPE Rn_GridComp
     integer                     ::  n        ! number of instances 
     TYPE(Rn_GridComp1), pointer ::  gcs(:)   ! instances
  END TYPE Rn_GridComp

CONTAINS

   subroutine RN_GridCompSetServices(  gc, chemReg, rc)
   type(ESMF_GridComp), intent(INOUT) :: GC
   type(Chem_Registry), intent(INOUT) :: chemReg
   integer,             intent(OUT  ) :: rc

   CHARACTER(LEN=255) :: rcbasen = 'Rn_GridComp'
   CHARACTER(LEN=255) :: name

   integer            :: status
   character(len=ESMF_MAXSTR) :: Iam

   integer :: ier,n,i
   CHARACTER(LEN=1) :: sOrP

   type(ESMF_Config) :: cfg

   Iam = "RN_GridCompSetServices"

!  Load resource file
!  ------------------
   cfg = ESMF_ConfigCreate(rc=status)
   VERIFY_(STATUS)
   call ESMF_ConfigLoadFile(cfg,TRIM(rcbasen)//'.rc',rc=status)
   VERIFY_(STATUS)

!  Parse resource file
!  -------------------
   n = ESMF_ConfigGetLen(cfg,label='Rn_instances:',rc=status)
   VERIFY_(STATUS)

!  We cannot have fewer instances than the number of
!  Rn bins in the registry (it is OK to have less, though)
!  --------------------------------------------------------
   IF( n < chemReg%n_Rn ) THEN
    status = 1
    VERIFY_(status)
   ELSE IF ( n >= chemReg%n_Rn ) THEN
    IF(MAPL_AM_I_ROOT()) THEN
     sOrP = " "
     IF(chemReg%n_Rn > 1) sOrP = "s"
     PRINT *, " "
     PRINT *, TRIM(Iam)//": Rn has ",chemReg%n_Rn," instantiation"//TRIM(sOrP)
    END IF
   END IF
   n = MIN(n,chemReg%n_Rn)

!  Record name of each instance
!  ----------------------------
   call ESMF_ConfigFindLabel(cfg,'Rn_instances:',rc=status)
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
      call RN_GridCompSetServices1_(gc,chemReg,name,rc=status)
      VERIFY_(STATUS)
   end do

   call MAPL_AddImportSpec(GC,           &
        SHORT_NAME = 'Rn_regionMask',    &
        LONG_NAME  = 'source species'  , &
        UNITS      = '1',                &
        DIMS       = MAPL_DimsHorzOnly,  &
        VLOCATION  = MAPL_VLocationNone, &
        RESTART    = MAPL_RestartSkip,   &
        RC         = STATUS)
   VERIFY_(STATUS)

   RETURN_(ESMF_SUCCESS)

   end subroutine RN_GridCompSetServices

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Rn_GridCompInitialize --- Initialize Rn_GridComp
!
! !INTERFACE:
!

   subroutine Rn_GridCompInitialize ( gcRn, w_c, impChem, expChem, &
                                      nymd, nhms, cdt, rc )

! !USES:

  IMPLICIT NONE

! !INPUT PARAMETERS:

   TYPE(Chem_Bundle), intent(in) :: w_c        ! Chemical tracer fields      
   INTEGER, INTENT(IN) :: nymd, nhms	       ! time
   REAL,    INTENT(IN) :: cdt		       ! chemical timestep (secs)


! !OUTPUT PARAMETERS:

   TYPE(Rn_GridComp), INTENT(INOUT) :: gcRn   ! Grid Component
   TYPE(ESMF_State), INTENT(INOUT)  :: impChem  ! Import State
   TYPE(ESMF_State), INTENT(INOUT)  :: expChem  ! Export State
   INTEGER, INTENT(OUT) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 - 

! !DESCRIPTION: Initializes the Rn Grid Component. Multiple instance
!               version.
!
! !REVISION HISTORY:
!
!  27Feb2008  da Silva  Introduced multiple instances
!  12Apr2008  Nielsen   Configured for radon
!
!EOP
!-------------------------------------------------------------------------

   CHARACTER(LEN=ESMF_MAXSTR), PARAMETER :: Iam = 'Rn_GridCompInitialize'
   CHARACTER(LEN=ESMF_MAXSTR) :: rcbasen = 'Rn_GridComp'
   CHARACTER(LEN=ESMF_MAXSTR) :: name
   
   INTEGER :: i, n, status
   CHARACTER(LEN=1) :: sOrP

!  Load resource file
!  ------------------
   CALL I90_loadf ( TRIM(rcbasen)//'.rc', status )
   VERIFY_(status)

!  Parse resource file
!  -------------------
   CALL I90_label ( 'Rn_instances:', status )
   VERIFY_(status)

!  First determine how many instances we have
!  ------------------------------------------   
   n = 0
   status = 0
   DO WHILE( status == 0 )
    CALL I90_gtoken( name, status )
    n = n + 1
   END DO
   IF( n == 0 ) THEN
    status = 1
    VERIFY_(status)
   END IF
   
!  We cannot have fewer instances than the number of
!  Rn bins in the registry (it is OK to have less, though)
!  --------------------------------------------------------
   IF( n < w_c%reg%n_Rn ) THEN
    status = 1
    VERIFY_(status)
   ELSE IF ( n >= w_c%reg%n_Rn ) THEN
    IF(MAPL_AM_I_ROOT()) THEN
     sOrP = " "
     IF(w_c%reg%n_Rn > 1) sOrP = "s"
     PRINT *, " " 
     PRINT *, TRIM(Iam)//": Rn has ",w_c%reg%n_Rn," instantiation"//TRIM(sOrP)
    END IF
   END IF
   n = MIN(n,w_c%reg%n_Rn)
   gcRn%n = n

!  Next allocate necessary memory
!  ------------------------------
   ALLOCATE ( gcRn%gcs(n), STAT=status )    
   VERIFY_(status)

!  Record name of each instance
!  ----------------------------
   CALL I90_label ( 'Rn_instances:', status )
   VERIFY_(status)
   DO i = 1, n
    CALL I90_gtoken( name, status )
    VERIFY_(status)
                                            ! resource file name
      gcRn%gcs(i)%rcfilen = TRIM(rcbasen)//'---'//TRIM(name)//'.rc'
      gcRn%gcs(i)%instance = i              ! instance number 
      IF(TRIM(name) == "full" ) THEN
       gcRn%gcs(i)%iname = " "              ! blank instance name for full (1)
      ELSE
       gcRn%gcs(i)%iname = TRIM(name)       ! instance name for others
      END IF
   END DO

!  Next initialize each instance
!  -----------------------------
   DO i = 1, gcRn%n
    IF(MAPL_AM_I_ROOT()) THEN
     PRINT *," "
     PRINT *, TRIM(Iam)//": Initializing instance ",TRIM(gcRn%gcs(i)%iname)," [",gcRn%gcs(i)%instance,"]"
    END IF
    CALL Rn_SingleInstance_ ( Rn_GridCompInitialize1_, i, &
                              gcRn%gcs(i), w_c, impChem, expChem,  &
                              nymd, nhms, cdt, status )
    VERIFY_(status)
   END DO

!  All done
!  --------
   CALL I90_FullRelease( status )
   VERIFY_(status)

 END SUBROUTINE Rn_GridCompInitialize

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Rn_GridCompRun --- Run Rn_GridComp
!
! !INTERFACE:
!

   SUBROUTINE Rn_GridCompRun ( gcRn, w_c, impChem, expChem, &
                                      nymd, nhms, cdt, rc )

! !USES:

  IMPLICIT NONE

! !INPUT PARAMETERS:

   TYPE(Chem_Bundle), intent(in) :: w_c        ! Chemical tracer fields      
   INTEGER, INTENT(IN) :: nymd, nhms	       ! time
   REAL,    INTENT(IN) :: cdt		       ! chemical timestep (secs)


! !OUTPUT PARAMETERS:

   TYPE(Rn_GridComp), INTENT(INOUT) :: gcRn   ! Grid Component
   TYPE(ESMF_State), INTENT(INOUT)  :: impChem  ! Import State
   TYPE(ESMF_State), INTENT(INOUT)  :: expChem  ! Export State
   INTEGER, INTENT(OUT) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 - 

! !DESCRIPTION: Runs the Rn Grid Component. Multiple instance
!               version.
!
! !REVISION HISTORY:
!
!  27Feb2008  da Silva  Introduced multiple instances
!  12Apr2008  Nielsen   Configured for radon
!
!EOP
!-------------------------------------------------------------------------

   CHARACTER(LEN=ESMF_MAXSTR), PARAMETER :: Iam = 'Rn_GridCompRun'
   INTEGER :: i, status

   DO i = 1, gcRn%n
    CALL Rn_SingleInstance_ ( Rn_GridCompRun1_, i, &
                              gcRn%gcs(i), w_c, impChem, expChem, &
                              nymd, nhms, cdt, status )
    VERIFY_(status)
   END DO

 END SUBROUTINE Rn_GridCompRun


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Rn_GridCompFinalize --- Initialize Rn_GridComp
!
! !INTERFACE:
!

   SUBROUTINE Rn_GridCompFinalize ( gcRn, w_c, impChem, expChem, &
                                      nymd, nhms, cdt, rc )

! !USES:

  IMPLICIT NONE

! !INPUT PARAMETERS:

   TYPE(Chem_Bundle), intent(in) :: w_c        ! Chemical tracer fields      
   INTEGER, INTENT(IN) :: nymd, nhms	       ! time
   REAL,    INTENT(IN) :: cdt		       ! chemical timestep (secs)


! !OUTPUT PARAMETERS:

   TYPE(Rn_GridComp), INTENT(INOUT) :: gcRn   ! Grid Component
   TYPE(ESMF_State), INTENT(INOUT)  :: impChem  ! Import State
   TYPE(ESMF_State), INTENT(INOUT)  :: expChem  ! Export State
   INTEGER, INTENT(OUT) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 - 

! !DESCRIPTION: Finalizes the Rn Grid Component. Multiple instance
!               version.
!
! !REVISION HISTORY:
!
!  27Feb2008  da Silva  Introduced multiple instances
!  12Apr2008  Nielsen   Configured for radon
!
!EOP
!-------------------------------------------------------------------------

   CHARACTER(LEN=ESMF_MAXSTR), PARAMETER :: Iam = 'Rn_GridCompFinalize'
   INTEGER :: i, status

   DO i = 1, gcRn%n
    CALL Rn_SingleInstance_ ( Rn_GridCompFinalize1_, i, &
                              gcRn%gcs(i), w_c, impChem, expChem, &
                              nymd, nhms, cdt, status )
    VERIFY_(status)
   END DO

   DEALLOCATE( gcRn%gcs, STAT=status )    
   gcRn%n = -1

 END SUBROUTINE Rn_GridCompFinalize

subroutine RN_GridCompSetServices1_(  gc, chemReg, iname, rc)
 type(ESMF_GridComp), intent(INOUT) :: GC
 type(Chem_Registry), intent(INOUT) :: chemReg
 character(len=*),    intent(IN   ) :: iname
 integer,             intent(OUT  ) :: rc

 integer :: Status
 character(len=ESMF_MAXSTR) :: Iam

 Iam ="RN_GridCOmpSetServices1_"

  call MAPL_AddImportSpec(GC,             &
       SHORT_NAME = 'Rn_EMISSION'//iname, &
       LONG_NAME  = 'source species'  ,   &
       UNITS      = '1',                  &
       DIMS       = MAPL_DimsHorzOnly,    &
       VLOCATION  = MAPL_VLocationNone,   &
       RESTART    = MAPL_RestartSkip,     &
       RC         = STATUS)
  VERIFY_(STATUS)

 RETURN_(ESMF_SUCCESS)

 end subroutine RN_GridCompSetServices1_

!--------------------------------------------------------------------------

!                      Single Instance Methods

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Rn_GridCompInitialize --- Initialize Rn_GridComp
!
! !INTERFACE:
!

   subroutine Rn_GridCompInitialize1_ ( gcRn, w_c, impChem, expChem, &
                                        nymd, nhms, cdt, rc )

! !USES:

  IMPLICIT NONE

! !INPUT PARAMETERS:

   TYPE(Chem_Bundle), intent(in) :: w_c        ! Chemical tracer fields      
   INTEGER, INTENT(IN) :: nymd, nhms	       ! time
   REAL,    INTENT(IN) :: cdt		       ! chemical timestep (secs)


! !OUTPUT PARAMETERS:

   TYPE(Rn_GridComp1), INTENT(INOUT) :: gcRn   ! Grid Component
   TYPE(ESMF_State), INTENT(INOUT)  :: impChem  ! Import State
   TYPE(ESMF_State), INTENT(INOUT)  :: expChem  ! Export State
   INTEGER, INTENT(OUT) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 - 

! !DESCRIPTION: Initializes the Rn Grid Component. It primarily sets
!               the import state for each active constituent package.
!
! !REVISION HISTORY:
!
!  18Sep2003 da Silva  First crack.
!  31May2005  Nielsen  Mods for 7 CO bins, 5 region masks
!  04Nov2005     Bian  CO tagged to 4 regions 
!                      (global, North America, South America, and Africa)
!                      for CR-AVE
!  12Apr2008  Nielsen  Configured for radon
!
!EOP
!-------------------------------------------------------------------------

   CHARACTER(LEN=ESMF_MAXSTR), PARAMETER :: Iam = 'Rn_GridCompInitialize1_'
   CHARACTER(LEN=ESMF_MAXSTR) :: rcfilen 

   INTEGER :: i1, i2, im, j, j1, j2, jm, km, n, status
   INTEGER :: nTimes, begTime, incSecs
   INTEGER :: nbeg, nend, nymd1, nhms1
   LOGICAL :: NoRegionalConstraint
   LOGICAL :: unitOK

   REAL :: conFac, limitN, limitS, log10Emission, radTODeg
   REAL, ALLOCATABLE :: var2d(:,:)

   rcfilen     = gcRn%rcfilen
   gcRn%name   = 'GEOS-5/GOCART Parameterized Radon Package'
   radTODeg    = 57.2957795
   gcRn%BCnymd = -1

!  Initialize local variables
!  --------------------------
   rc = 0
   status = 0

   i1 = w_c%grid%i1
   i2 = w_c%grid%i2
   im = w_c%grid%im

   j1 = w_c%grid%j1
   j2 = w_c%grid%j2
   jm = w_c%grid%jm

   km = w_c%grid%km

   nbeg  = w_c%reg%i_Rn
   nend  = w_c%reg%j_Rn

!  It requires 1 bin
!  -----------------
   IF( nbeg /= nend ) THEN
    IF(MAPL_AM_I_ROOT()) PRINT *, TRIM(Iam)//": Must have only 1 bin at the single instance level"
    status = 1
    VERIFY_(status)
   END IF

!  Load resource file
!  ------------------
   CALL I90_loadf ( TRIM(rcfilen), status )
   VERIFY_(status)

! Run-time debug switch
! ---------------------
   CALL I90_label ( 'DEBUG:', status )
   VERIFY_(status)
   n = I90_gint ( status )
   VERIFY_(status)
   IF(n /= 0) THEN
    gcRn%DebugIsOn = .TRUE.
   ELSE
    gcRn%DebugIsOn = .FALSE.
   END IF

!  Allocate
!  --------
   ALLOCATE( gcRn%RnsfcFlux(i1:i2,j1:j2), gcRn%regionMask(i1:i2,j1:j2), &
	     gcRn%ScheryEmission(i1:i2,j1:j2), STAT=status )
   VERIFY_(status)

!  Obtain half life.  Unit must be "years", "days", or "seconds".
!  --------------------------------------------------------------
   CALL I90_label ( 'HalfLife:', status )
   VERIFY_(status)
   gcRn%halfLife = I90_gfloat ( status )
   VERIFY_(status)
   CALL I90_label ( 'HalfLifeUnit:', status )
   VERIFY_(status)
   CALL I90_gtoken( gcRn%halfLifeUnit, status )
   VERIFY_(status)

!  Validate the specified half-life and units, and find
!  the constant needed to convert half-life to seconds.
!  ----------------------------------------------------
   unitOK = .FALSE.
   IF(TRIM(gcRn%halfLifeUnit) ==   "years") THEN
    unitOK = .TRUE.
    conFac = 86400.00*365.25
   END IF
   IF(TRIM(gcRn%halfLifeUnit) ==    "days") THEN
    unitOK = .TRUE.
    conFac = 86400.00
   END IF
   IF(TRIM(gcRn%halfLifeUnit) == "seconds") THEN
    unitOK = .TRUE.
    conFac = 1.00
   END IF
   IF( .NOT. unitOK ) THEN
    IF(MAPL_AM_I_ROOT()) PRINT *, TRIM(Iam)//": Invalid unit specified for radon half-life."
    status = 1
    VERIFY_(status)
   END IF
   IF(gcRn%halfLife <= 0.00) THEN
    IF(MAPL_AM_I_ROOT()) PRINT *, TRIM(Iam)//": Radon half-life must be greater than zero."
    VERIFY_(status)
   END IF

!  Compute the decay constant (inverse seconds) from the half-life:
!    ln(N/No) = ln(1/2) = -decayConstant * halfLife
!  ----------------------------------------------------------------
   gcRn%decayConstant = 0.693147/(gcRn%halfLife*conFac)

!  Grab the region string.
!  -----------------------
   CALL I90_label ( 'Rn_regions_indices:', status )
   VERIFY_(status)
   CALL I90_gtoken( gcRn%regionsString, status )
   VERIFY_(status)

!  Is this instantiation a global case?
!  -----------------------------------
   IF(gcRn%regionsString(1:2) == "-1") THEN
    NoRegionalConstraint = .TRUE.
   ELSE
    SELECT CASE (ESMF_UtilStringLowerCase(gcRn%regionsString(1:2)))
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
   IF(NoRegionalConstraint) gcRn%regionsString = "-1"

   IF(MAPL_AM_I_ROOT()) THEN
    IF(NoRegionalConstraint) THEN
     PRINT *, TRIM(Iam)//": This instantiation has no regional constraints."
    ELSE
     PRINT *, TRIM(Iam)//": This instantiation is regionally constrained."
     PRINT *, TRIM(Iam)//": List of region numbers included: ",TRIM(gcRn%regionsString)
    END IF
   END IF

!  Set the initial radon surface fluxes to zero
!  --------------------------------------------
   gcRn%RnsfcFlux(i1:i2,j1:j2) = 0.00

   RETURN

 END SUBROUTINE Rn_GridCompInitialize1_

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Rn_GridCompRun
!
! !INTERFACE:
!

   SUBROUTINE Rn_GridCompRun1_ ( gcRn, w_c, impChem, expChem, &
                                 nymd, nhms, cdt, rc )

! !USES:

  IMPLICIT NONE

! !INPUT/OUTPUT PARAMETERS:

   TYPE(Rn_GridComp1), INTENT(INOUT) :: gcRn   ! Grid Component
   TYPE(Chem_Bundle), INTENT(INOUT) :: w_c	! Chemical tracer fields   

! !INPUT PARAMETERS:

   TYPE(ESMF_State), INTENT(inout) :: impChem    ! Import State
   INTEGER, INTENT(IN) :: nymd, nhms	      ! time
   REAL,    INTENT(IN) :: cdt		      ! chemical timestep (secs)


! !OUTPUT PARAMETERS:

   TYPE(ESMF_State), intent(inout) :: expChem     ! Export State
   INTEGER, INTENT(OUT) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 -
 
! !DESCRIPTION: This routine implements the Rn driver.
!
! !REVISION HISTORY:
!
!  18Sep2003 da Silva  First crack.
!  12Apr2008  Nielsen  Configured for radon
!
!EOP
!-------------------------------------------------------------------------

   CHARACTER(LEN=ESMF_MAXSTR), PARAMETER :: Iam = 'Rn_GridCompRun1_'

!  Input fields from fvGCM
!  -----------------------
   REAL, POINTER, DIMENSION(:,:,:) :: T       => null()
   REAL, POINTER, DIMENSION(:,:,:) :: zle     => null()
   REAL, POINTER, DIMENSION(:,:)   :: soilT   => null()
   REAL, POINTER, DIMENSION(:,:)   :: fracIce => null()

   INTEGER :: i1, i2, im, j1, j2, jm, km, idiag, iXj
   INTEGER :: i, j, k, kReverse, n, nbeg, nend, nymd1
   INTEGER :: status

   REAL, PARAMETER :: nsuba=6.022E+26
   REAL, PARAMETER :: mwtAir=28.97
   REAL, PARAMETER :: mwtRn=222.00
   REAL, PARAMETER :: rstar=8.3143E+03
   REAL, PARAMETER :: rpstd=1.00E-05

   REAL    :: decadence, qmin, qmax, toND
   REAL, ALLOCATABLE :: F(:,:),nd(:,:,:),p(:,:,:),pe(:,:,:),dZ(:,:,:)
   INTEGER, ALLOCATABLE :: mask(:,:)
   real, pointer        :: ptr2d(:,:) => null()

#define EXPORT   expChem
#define iNAME    TRIM(gcRn%iname)

#define RnEM     Rn_emis
#define RnCL     Rn_column
#define RnSC     Rn_surface
#define RnLS     Rn_loss

#include "Rn_GetPointer___.h"

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

   nbeg  = w_c%reg%i_Rn
   nend  = w_c%reg%j_Rn

!  Get the region mask
!  -------------------
   call MAPL_GetPointer(impChem,ptr2d,'Rn_regionMask',rc=status)
   VERIFY_(STATUS)
   gcRn%regionMask=ptr2d

!  It requires 1 bin
!  -----------------
   IF ( nbeg /= nend ) THEN
    IF(MAPL_AM_I_ROOT()) PRINT *, TRIM(Iam)//": Must have only 1 bin at the single instance level"
    status = 1
    VERIFY_(status)
   END IF

   call MAPL_GetPointer(impChem, ptr2d, 'Rn_EMISSION'//iNAME,rc=status)
   VERIFY_(STATUS)
   gcRn%ScheryEmission = ptr2d

!  Conversion factor: mBq m^{-2} to atoms m^{-2} s^{-1}
!  ----------------------------------------------------
   toND = 0.001/gcRn%decayConstant

!  Allocate temporary workspace
!  ----------------------------
   ALLOCATE(pe(i1:i2,j1:j2,km+1), p(i1:i2,j1:j2,km), nd(i1:i2,j1:j2,km), &
            dZ(i1:i2,j1:j2,km), F(i1:i2,j1:j2), mask(i1:i2,j1:j2), STAT=status)
   VERIFY_(status)

!  Get imports
!  -----------
   CALL MAPL_GetPointer( impChem,	T,	'T', RC=status ) 
   VERIFY_(status)
   CALL MAPL_GetPointer( impChem,     zle,    'ZLE', RC=status ) 
   VERIFY_(status)
   CALL MAPL_GetPointer( impChem,   soilT, 'TSOIL1', RC=status ) 
   VERIFY_(status)
   CALL MAPL_GetPointer( impChem, fracIce,  'FRACI', RC=status ) 
   VERIFY_(status)

!  Layer thicknesses.  ZLE(:,:,0:km).
!  ----------------------------------
   DO k=1,km
    dZ(i1:i2,j1:j2,k) = zle(i1:i2,j1:j2,k-1)-zle(i1:i2,j1:j2,k)
   END DO

!  Layer interface pressures
!  -------------------------
   pe(i1:i2,j1:j2,1)=w_c%grid%ptop
   DO k=2,km+1
    pe(i1:i2,j1:j2,k)=pe(i1:i2,j1:j2,k-1)+w_c%delp(i1:i2,j1:j2,k-1)
   END DO

!  Layer mean pressures
!  --------------------
   DO k=1,km
    p(i1:i2,j1:j2,k)=(pe(i1:i2,j1:j2,k)+pe(i1:i2,j1:j2,k+1))*0.50
   END DO

!  Number density
!  --------------
   nd(i1:i2,j1:j2,1:km)= nsuba*p(i1:i2,j1:j2,1:km)/ &
                        (rstar*T(i1:i2,j1:j2,1:km))

!  Validate
!  --------
   IF(gcRn%DebugIsOn) THEN
    CALL pmaxmin('Rn: T     ',      T, qmin, qmax, iXj,   km, 1. )
    CALL pmaxmin('Rn: TSOIL1',  soilT, qmin, qmax, iXj,    1, 1. )
    CALL pmaxmin('Rn: FRACI ',fracIce, qmin, qmax, iXj,    1, 1. )
    CALL pmaxmin('Rn: ZLE   ',    zle, qmin, qmax, iXj, km+1, 1. )
    CALL pmaxmin('Rn: dZ    ',     dZ, qmin, qmax, iXj,   km, 1. )
    CALL pmaxmin('Rn: Edge p',     pe, qmin, qmax, iXj,   km, 1. )
    CALL pmaxmin('Rn: Mid p ',      p, qmin, qmax, iXj,   km, 1. )
    CALL pmaxmin('Rn: Numden',     nd, qmin, qmax, iXj,   km, 1. )
   END IF

!  Convert Radon from mole fraction to number density
!  --------------------------------------------------
   w_c%qa(nbeg)%data3d(i1:i2,j1:j2,1:km) = &
                  w_c%qa(nbeg)%data3d(i1:i2,j1:j2,1:km)*nd(i1:i2,j1:j2,1:km)

!  Clear surface flux diagnostic
!  -----------------------------
   IF(ASSOCIATED(Rn_emis)) Rn_emis(i1:i2,j1:j2) = 0.00

!  Fraction of emissions to accept from each box
!  ---------------------------------------------
   F(i1:i2,j1:j2) = 0.00

!  Find land boxes from regional mask file.
!  ----------------------------------------
   CALL setLandMask(status)
   VERIFY_(status)

!  For the global intantiation, include ocean emissions.
!  -----------------------------------------------------
   IF(gcRn%instance == 1) THEN
    WHERE(mask(i1:i2,j1:j2) == 0) F(i1:i2,j1:j2) = 1.00-fracIce(i1:i2,j1:j2)
   END IF

!  Assume frozen soil emits no radon
!  ---------------------------------
   WHERE(mask(i1:i2,j1:j2) == 1 .AND. soilT(i1:i2,j1:j2) < 273.00) mask(i1:i2,j1:j2) = 0

!  Account for bad-valued soil temperatures
!  ----------------------------------------
   WHERE(mask(i1:i2,j1:j2) == 1 .AND. soilT(i1:i2,j1:j2) > 500.00) mask(i1:i2,j1:j2) = 0

!  Finalize fraction from land boxes.
!  ----------------------------------
   WHERE(mask(i1:i2,j1:j2) == 1) F(i1:i2,j1:j2) = 1.00

!  Place emissions into the surface layer, adding the number of
!  atoms released in one time step to the surface layer number density.
!  --------------------------------------------------------------------
   w_c%qa(nbeg)%data3d(i1:i2,j1:j2,km)=w_c%qa(nbeg)%data3d(i1:i2,j1:j2,km)+F(i1:i2,j1:j2)* &
                                       toND*gcRn%ScheryEmission(i1:i2,j1:j2)*cdt/dZ(i1:i2,j1:j2,km)

!  Diagnostic emissions, kg m^{-2} s^{-1}.
!  ---------------------------------------
   IF(ASSOCIATED(Rn_emis)) THEN
    Rn_emis(i1:i2,j1:j2) = F(i1:i2,j1:j2)*toND*mwtRn*gcRn%ScheryEmission(i1:i2,j1:j2)/nsuba
   END IF

!  Diagnostic loss, vertically integrated, kg m^{-2} s^{-1}.  Compute
!  before applying radioactive decay to the three-dimensional radon field.
!  -----------------------------------------------------------------------
   n = gcRn%instance
   decadence = EXP(-gcRn%decayConstant*cdt)
   IF(ASSOCIATED(Rn_loss)) Rn_loss(i1:i2,j1:j2) = 0.00
   DO k = 1, km

    IF(ASSOCIATED(Rn_loss)) &
    Rn_loss(i1:i2,j1:j2) = w_c%qa(nbeg)%data3d(i1:i2,j1:j2,k)*(1.00-decadence)* &
   	                   mwtRn*dZ(i1:i2,j1:j2,k)/(nsuba*cdt)

!  Apply radioactive decay, q(f) = q(i)EXP(-c delta t), to number density.
!  -----------------------------------------------------------------------
    w_c%qa(nbeg)%data3d(i1:i2,j1:j2,k) = w_c%qa(nbeg)%data3d(i1:i2,j1:j2,k)*decadence

   END DO ! Next layer, k

!  Column burden in kg m^{-2}
!  --------------------------
   n = gcRn%instance 
   IF(ASSOCIATED(Rn_column)) then
    Rn_column(i1:i2,j1:j2) = 0.
    DO k = 1, km
     Rn_column(i1:i2,j1:j2) = Rn_column(i1:i2,j1:j2) +  w_c%qa(nbeg)%data3d(i1:i2,j1:j2,k)* &
   			      mwtRn*dZ(i1:i2,j1:j2,k)/nsuba
    END DO
   END IF

!  Return to mole fraction
!  -----------------------
   w_c%qa(nbeg)%data3d(i1:i2,j1:j2,1:km)=w_c%qa(nbeg)%data3d(i1:i2,j1:j2,1:km)/nd(i1:i2,j1:j2,1:km)

!  Surface concentration in mole fraction
!  --------------------------------------
   n = gcRn%instance 
   IF(ASSOCIATED(Rn_surface)) Rn_surface(i1:i2,j1:j2) = w_c%qa(nbeg)%data3d(i1:i2,j1:j2,km)

   IF(gcRn%DebugIsOn) THEN
    n = gcRn%instance 
    IF(ASSOCIATED(   Rn_emis)) &
     CALL pmaxmin(   'Rn_emis',	   Rn_emis(i1:i2,j1:j2), qmin, qmax, iXj, 1, 1. )
    IF(ASSOCIATED(   Rn_loss)) &
     CALL pmaxmin(   'Rn_loss',	   Rn_loss(i1:i2,j1:j2), qmin, qmax, iXj, 1, 1. )
    IF(ASSOCIATED( Rn_column)) &
     CALL pmaxmin( 'Rn_column',  Rn_column(i1:i2,j1:j2), qmin, qmax, iXj, 1, 1. )
    IF(ASSOCIATED(Rn_surface)) &
     CALL pmaxmin('Rn:surface', Rn_surface(i1:i2,j1:j2), qmin, qmax, iXj, 1, 1. )
   END IF

!  Housekeeping
!  ------------
   DEALLOCATE(F, mask, dZ, nd, p, pe, STAT=status)
   VERIFY_(status)

   RETURN

CONTAINS

  SUBROUTINE setLandMask(rc)
   IMPLICIT NONE
   INTEGER, INTENT(OUT) ::  rc
   INTEGER :: i, k, status
   INTEGER, ALLOCATABLE :: regionNumbers(:),flag(:)
   CHARACTER(LEN=ESMF_MAXSTR), PARAMETER :: Iam = 'Rn::setLandMask'

   rc = 0
   k = 32
   ALLOCATE(regionNumbers(k), flag(k), STAT=status)
   VERIFY_(status)

! Obtain region numbers from delimited list of integers
! -----------------------------------------------------
   regionNumbers(:) = 0
   CALL Chem_UtilExtractIntegers(gcRn%regionsString, k, regionNumbers, RC=status)
   VERIFY_(status)

! How many integers were found?
! -----------------------------
   flag(:) = 1
   WHERE(regionNumbers(:) == 0) flag(:) = 0
   k = SUM(flag)
   DEALLOCATE(flag, STAT=status)
   VERIFY_(status)

! Set local mask to 1 where gridMask matches each integer (within precision!).
! ----------------------------------------------------------------------------
   mask(i1:i2,j1:j2) = 0
   IF(regionNumbers(1) == -1) THEN
    WHERE(gcRn%regionMask(i1:i2,j1:j2) /= 0) mask(i1:i2,j1:j2) = 1
   ELSE
    DO i = 1,k
     WHERE(       regionNumbers(i)-0.01 <= gcRn%regionMask(i1:i2,j1:j2) .AND. &
           gcRn%regionMask(i1:i2,j1:j2) <= regionNumbers(i)+0.01) mask(i1:i2,j1:j2) = 1
    END DO
   END IF

   RETURN
  END SUBROUTINE setLandMask

 END SUBROUTINE Rn_GridCompRun1_

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Rn_GridCompFinalize --- The Chem Driver 
!
! !INTERFACE:
!

   SUBROUTINE Rn_GridCompFinalize1_ ( gcRn, w_c, impChem, expChem, &
                                      nymd, nhms, cdt, rc )

! !USES:

  IMPLICIT NONE

! !INPUT/OUTPUT PARAMETERS:

   TYPE(Rn_GridComp1), INTENT(INOUT) :: gcRn   ! Grid Component

! !INPUT PARAMETERS:

   TYPE(Chem_Bundle), INTENT(IN)  :: w_c      ! Chemical tracer fields   
   INTEGER, INTENT(IN) :: nymd, nhms	      ! time
   REAL,    INTENT(IN) :: cdt		      ! chemical timestep (secs)


! !OUTPUT PARAMETERS:

   TYPE(ESMF_State), INTENT(INOUT) :: impChem	! Import State
   TYPE(ESMF_State), INTENT(INOUT) :: expChem	! Import State
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

   CHARACTER(LEN=ESMF_MAXSTR), PARAMETER :: Iam = 'Rn_GridCompFinalize1_'
   INTEGER :: status
   rc = 0

   DEALLOCATE ( gcRn%RnsfcFlux,  gcRn%regionMask, gcRn%ScheryEmission, STAT=status )
   VERIFY_(status)

   RETURN

 END SUBROUTINE Rn_GridCompFinalize1_

 END MODULE Rn_GridCompMod

!-----------------------------------------------------------------------

!                     Single Instance Wrapper

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Rn_SingleInstance_ --- Runs single instance of method
!
! !INTERFACE:
!
  SUBROUTINE Rn_SingleInstance_ ( Method_, instance, &
                                  gcRn, w_c, impChem, expChem, &
                                  nymd, nhms, cdt, rc )

! !USES:

  USE Rn_GridCompMod
  USE ESMF
  USE MAPL
  USE Chem_Mod 

  IMPLICIT NONE

! !INPUT PARAMETERS:

!  Input "function pointer"
!  -----------------------
   INTERFACE 
     SUBROUTINE Method_ (gc, w, imp, exp, ymd, hms, dt, rcode )
       USE Rn_GridCompMod
       USE ESMF
       USE MAPL
       USE Chem_Mod 
       TYPE(Rn_GridComp1),  INTENT(INOUT)  :: gc
       TYPE(Chem_Bundle),   INTENT(IN)     :: w
       TYPE(ESMF_State),    INTENT(INOUT)  :: imp
       TYPE(ESMF_State),    INTENT(INOUT)  :: exp
       INTEGER,    	    INTENT(IN)     :: ymd, hms
       REAL,	   	    INTENT(IN)     :: dt	
       INTEGER,    	    INTENT(OUT)    :: rcode
     END SUBROUTINE Method_
   END INTERFACE

   INTEGER, INTENT(IN)           :: instance   ! instance number

   TYPE(Chem_Bundle), INTENT(INOUT) :: w_c     ! Chemical tracer fields      
   INTEGER, INTENT(IN) :: nymd, nhms	       ! time
   REAL,    INTENT(IN) :: cdt		       ! chemical timestep (secs)


! !OUTPUT PARAMETERS:

   TYPE(Rn_GridComp1), INTENT(INOUT) :: gcRn    ! Grid Component
   TYPE(ESMF_State), INTENT(INOUT)  :: impChem  ! Import State
   TYPE(ESMF_State), INTENT(INOUT)  :: expChem  ! Export State
   INTEGER, INTENT(OUT) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 - 

! !DESCRIPTION: Finalizes the Rn Grid Component. Multiple instance
!               version.
!
! !REVISION HISTORY:
!
!  27Feb2008  da Silva  Introduced multiple instances
!  12Apr2008  Nielsen   Configured for radon.
!
!EOP
!-------------------------------------------------------------------------

  INTEGER :: n_Rn, i_Rn, j_Rn

! Save overall Rn indices
! -----------------------
  n_Rn = w_c%reg%n_Rn
  i_Rn = w_c%reg%i_Rn
  j_Rn = w_c%reg%j_Rn
  
! Customize indices for this particular instance
! ----------------------------------------------
  w_c%reg%n_Rn = 1
  w_c%reg%i_Rn = i_Rn + instance - 1
  w_c%reg%j_Rn = i_Rn + instance - 1
  
! Execute the instance method
! ---------------------------
  CALL Method_ ( gcRn, w_c, impChem, expChem, nymd, nhms, cdt, rc )

! Restore the overall Rn indices
! ------------------------------
  w_c%reg%n_Rn = n_Rn
  w_c%reg%i_Rn = i_Rn
  w_c%reg%j_Rn = j_Rn

  END SUBROUTINE Rn_SingleInstance_

!-----------------------------------------------------------------------
