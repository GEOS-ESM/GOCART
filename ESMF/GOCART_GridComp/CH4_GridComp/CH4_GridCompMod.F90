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
   USE Chem_Mod 	     ! Chemistry Base Class
   USE Chem_StateMod	     ! Chemistry State
   USE Chem_UtilMod	     ! I/O
   USE m_inpak90	     ! Resource file management
   USE Henrys_law_ConstantsMod, ONLY: get_HenrysLawCts

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
!
!EOP
!-------------------------------------------------------------------------

  TYPE CH4_GridComp1

   CHARACTER(LEN=ESMF_MAXSTR) :: name		 ! generic name of the package
   CHARACTER(LEN=ESMF_MAXSTR) :: iname  	 ! instance name
   CHARACTER(LEN=ESMF_MAXSTR) :: rcfilen	 ! resource file name
   CHARACTER(LEN=ESMF_MAXSTR) :: regionsString   ! Comma-delimited string of regions
   CHARACTER(LEN=ESMF_MAXSTR) :: CH4Source	 ! Source name on emission file (CH4_ANIMLS, for example)

   INTEGER :: instance                 ! Instantiation number
   INTEGER :: nymd_oh
   INTEGER :: nymd_ch4
   INTEGER :: BCnymd                   ! Date of last emissions/prodction read

   REAL, POINTER :: regionMask(:,:)    ! regional mask
   REAL, POINTER ::	 CH4(:,:,:)    ! CH4 mixing ratio mol/mol
   REAL, POINTER ::	OHnd(:,:,:)    ! OH number density (cm^{-3})
   REAL, POINTER :: CH4sfcFlux(:,:)    ! CH4 surface flux kg m^-2 s^-1

   LOGICAL :: DebugIsOn     ! Run-time debug switch
   LOGICAL :: CH4FeedBack   ! Permit increments to CH4 from CH4 + hv => 2H2O + CO
   LOGICAL :: H2OFeedBack   ! Permit increments to   Q from CH4 + hv => 2H2O + CO

   REAL :: szaCutoff        ! Largest solar zenith angle (degrees) allowed as daytime
	
  END TYPE CH4_GridComp1

  TYPE CH4_GridComp
     INTEGER                      ::  n        ! number of instances 
     TYPE(CH4_GridComp1), pointer ::  gcs(:)   ! instances
  END TYPE CH4_GridComp

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
!   CH4 bins in the registry (it is OK to have less, though)
!  --------------------------------------------------------
   if ( n .LT. chemReg%n_CH4 ) then
        rc = 35
        return
   else if ( n .GT. chemReg%n_CH4 ) then
        if (MAPL_AM_I_ROOT()) &
        PRINT *, TRIM(Iam)//": Bins = ",chemReg%n_CH4," of ",n," expected."
   end if
   n = min(n,chemReg%n_CH4 )

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

   call MAPL_AddImportSpec(GC,           &
        SHORT_NAME = 'CH4_regionMask',   &
        LONG_NAME  = 'source species'  , &
        UNITS      = '1',                &
        DIMS       = MAPL_DimsHorzOnly,  &
        VLOCATION  = MAPL_VLocationNone, &
        RESTART    = MAPL_RestartSkip,   &
        RC         = STATUS)
   VERIFY_(STATUS)

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
   INTEGER, INTENT(IN) :: nymd, nhms	       ! time
   REAL,    INTENT(IN) :: cdt		       ! chemical timestep (secs)


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
   
   INTEGER :: i, n, status
   REAL :: c1,c2,c3,c4

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
!   CH4 bins in the registry (it is OK to have less, though)
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
!  -----------------------------------------------------------
   CALL get_HenrysLawCts('CH4',c1,c2,c3,c4)  
   w_c%reg%Hcts(1,w_c%reg%i_CH4 : w_c%reg%j_CH4)=c1
   w_c%reg%Hcts(2,w_c%reg%i_CH4 : w_c%reg%j_CH4)=c2
   w_c%reg%Hcts(3,w_c%reg%i_CH4 : w_c%reg%j_CH4)=c3
   w_c%reg%Hcts(4,w_c%reg%i_CH4 : w_c%reg%j_CH4)=c4

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
   INTEGER, INTENT(IN) :: nymd, nhms	       ! time
   REAL,    INTENT(IN) :: cdt		       ! chemical timestep (secs)


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
   INTEGER, INTENT(IN) :: nymd, nhms	       ! time
   REAL,    INTENT(IN) :: cdt		       ! chemical timestep (secs)


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

  call MAPL_AddImportSpec(GC, &
       SHORT_NAME = 'CH4_sfcFlux'//iname, &
       LONG_NAME  = 'source species'  , &
       UNITS      = '1', &
       DIMS       = MAPL_DimsHorzOnly, &
       VLOCATION  = MAPL_VLocationNone, &
       RC         = STATUS)
  VERIFY_(STATUS)
  call MAPL_AddImportSpec(GC, &
       SHORT_NAME = 'CH4_oh'//iname, &
       LONG_NAME  = 'source species'  , &
       UNITS      = '1', &
       DIMS       = MAPL_DimsHorzVert, &
       VLOCATION  = MAPL_VLocationCenter, &
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
   INTEGER, INTENT(IN) :: nymd, nhms	       ! time
   REAL,    INTENT(IN) :: cdt		       ! chemical timestep (secs)


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
   LOGICAL :: NoRegionalConstraint 

   REAL :: limitN, limitS, radTODeg
   REAL, ALLOCATABLE :: var2D(:,:)
   real, pointer     :: ptr2d(:,:) => null()

   rcfilen = gcCH4%rcfilen
   gcCH4%name = 'GEOS-5/GOCART Parameterized CH4 Package'
   gcCH4%BCnymd = -1
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

!  Allocate memory, etc
!  --------------------
   ALLOCATE ( gcCH4%CH4sfcFlux(i1:i2,j1:j2), &
              gcCH4%regionMask(i1:i2,j1:j2), &
              gcCH4%OHnd(i1:i2,j1:j2,km), STAT=status )
   VERIFY_(status)

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

   call MAPL_GetPointer(impChem,ptr2D,'CH4_regionMask',rc=status)
   VERIFY_(STATUS)

!  Grab the region string.
!  -----------------------
   CALL I90_label ( 'CH4_regions_indices:', status )
   VERIFY_(status)
   CALL I90_gtoken( gcCH4%regionsString, status )
   VERIFY_(status)

!  Is this instantiation a global case?
!  -----------------------------------
   IF(gcCH4%regionsString(1:2) == "-1") THEN
    NoRegionalConstraint = .TRUE.
   ELSE
    SELECT CASE (ESMF_UtilStringLowerCase(gcCH4%regionsString(1:2)))
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
   IF(NoRegionalConstraint) gcCH4%regionsString = "-1"

   IF(MAPL_AM_I_ROOT()) THEN
    IF(NoRegionalConstraint) THEN
     PRINT *,TRIM(Iam)//": This instantiation has no regional constraints."
    ELSE
     PRINT *,TRIM(Iam)//": This instantiation is regionally constrained."
     PRINT *,TRIM(Iam)//": List of region numbers included: ",TRIM(gcCH4%regionsString)
    END IF
   END IF

!  Set the initial CH4 surface fluxes to zero
!  ------------------------------------------
   gcCH4%CH4sfcFlux(i1:i2,j1:j2) = 0.00

!  Use instance name as key to CH4 emission source
!  -----------------------------------------------
   gcCH4%CH4Source = "CH4_"//TRIM(ESMF_UtilStringUpperCase(gcCH4%iname))

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
   REAL, POINTER, DIMENSION(:,:,:) ::  T        => null()
   REAL, POINTER, DIMENSION(:,:,:) ::  Q        => null()
   REAL, POINTER, DIMENSION(:,:,:) ::  rhoWet   => null()
   REAL, POINTER, DIMENSION(:,:,:) ::  rhoDry   => null()
   REAL, POINTER, DIMENSION(:,:,:) ::  ple      => null()

   INTEGER :: i1, i2, im, j1, j2, jm, km, idiag, iXj
   INTEGER :: i, j, k, kReverse, n, nbeg, nend
   INTEGER :: nymd1, nhms1, status

   REAL, PARAMETER :: mwtCH4 = 16.043

   REAL    :: qmin, qmax

   REAL, ALLOCATABLE :: cellDepth(:,:,:), cellVolume(:,:,:),     rkoh(:,:,:)
   REAL, ALLOCATABLE ::         p(:,:,:),      ndWet(:,:,:),    ndDry(:,:,:)
   REAL, ALLOCATABLE ::    dCH4nd(:,:,:),      photJ(:,:,:), dCH4Phot(:,:,:)

   CHARACTER(LEN=256) :: CH4Source
 
   REAL, POINTER :: ptr2d(:,:)   => null()
   REAL, POINTER :: ptr3d(:,:,:) => null()

#define EXPORT   expChem
#define iNAME    TRIM(gcCH4%iname)

#define CH4EM	CH4_emis
#define CH4CL	CH4_column
#define CH4SC	CH4_surface
#define CH4PD	CH4_prod
#define CH4LS	CH4_loss
#define CH4JL	CH4_phot
#define CH4QP	CH4_qprod
#define CH4DRY	CH4_dry

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
   IF ( nbeg /= nend ) THEN
    IF(MAPL_AM_I_ROOT()) PRINT *,TRIM(Iam)//": Must have only 1 bin at the single instance level"
    status = 1
    VERIFY_(status)
   END IF

!  Get imports
!  -----------
   CALL MAPL_GetPointer(impChem, T,	   'T',           RC=status)
   VERIFY_(status)
   CALL MAPL_GetPointer(impChem, Q,	   'Q',           RC=status)
   VERIFY_(status)
   CALL MAPL_GetPointer(impChem, rhoWet,   'AIRDENS',     RC=status)
   VERIFY_(status)
   CALL MAPL_GetPointer(impChem, rhoDry,   'AIRDENS_DRYP',RC=status)
   VERIFY_(status)
   CALL MAPL_GetPointer(impChem, cellArea, 'AREA',        RC=status)
   VERIFY_(status)
   CALL MAPL_GetPointer(impChem, ple,	   'PLE',         RC=status)
   VERIFY_(status)

   IF(gcCH4%DebugIsOn) THEN
    CALL pmaxmin('CH4:AREA', cellArea, qmin, qmax, iXj,  1,   1. )
    CALL pmaxmin('CH4:T',           T, qmin, qmax, iXj, km,   1. )
    CALL pmaxmin('CH4:Q',           q, qmin, qmax, iXj, km,   1. )
    CALL pmaxmin('CH4:rhoWet', rhoWet, qmin, qmax, iXj, km,   1. )
    CALL pmaxmin('CH4:rhoDry', rhoDry, qmin, qmax, iXj, km,   1. )
    CALL pmaxmin('CH4:ple',       ple, qmin, qmax, iXj, km+1, 1. )
   END IF

!  Update CH4 emissions and OH number density once each day.
!  The latter appears to be in molecules cm^-3.
!  ---------------------------------------------------------
    call MAPL_GetPointer(impChem,ptr2d,'CH4_sfcFlux'//trim(iNAME),rc=status)
    VERIFY_(STATUS)
    gcCH4%CH4sfcFlux=ptr2d

    call MAPL_GetPointer(impChem,ptr3d,'CH4_oh'//trim(iNAME),rc=status)
    VERIFY_(STATUS)
    gcCH4%OHnd=ptr3d

!  OH number density from molecules cm^-3 to molecules m^-3
!  --------------------------------------------------------
    gcCH4%OHnd(i1:i2,j1:j2,1:km) = gcCH4%OHnd(i1:i2,j1:j2,1:km)*1.00E+06

!  Allocate temporary workspace
!  ----------------------------
   ALLOCATE( p(i1:i2,j1:j2,km), &
         ndWet(i1:i2,j1:j2,km), &
         ndDry(i1:i2,j1:j2,km), &
          rkoh(i1:i2,j1:j2,km), &
        dCH4nd(i1:i2,j1:j2,km), &
    cellVolume(i1:i2,j1:j2,km), &
     cellDepth(i1:i2,j1:j2,km), STAT=status)
   VERIFY_(status)

!  Layer mean pressures
!  --------------------
   DO k=1,km
    p(i1:i2,j1:j2,k) = (ple(i1:i2,j1:j2,k-1)+ple(i1:i2,j1:j2,k))*0.50
   END DO
 
!  Moist air number density
!  ------------------------
   ndWet(i1:i2,j1:j2,1:km) = rhoWet(i1:i2,j1:j2,1:km)*MAPL_AVOGAD/MAPL_AIRMW
 
!  Dry air number density
!  ----------------------
   ndDry(i1:i2,j1:j2,1:km) = rhoDry(i1:i2,j1:j2,1:km)*MAPL_AVOGAD/MAPL_AIRMW

!  Cell depth and volume
!  ---------------------
   DO k=1,km
    cellDepth(i1:i2,j1:j2,k) = (ple(i1:i2,j1:j2,k)-ple(i1:i2,j1:j2,k-1))/(rhoWet(i1:i2,j1:j2,k)*MAPL_GRAV)
    cellVolume(i1:i2,j1:j2,k) = cellArea(i1:i2,j1:j2)*cellDepth(i1:i2,j1:j2,k)
   END DO

   IF(gcCH4%DebugIsOn) THEN
    CALL pmaxmin('CH4: SfcFlux', gcCH4%CH4sfcFlux, qmin, qmax, iXj, 1,  1. )
    CALL pmaxmin('CH4: OH Conc',       gcCH4%OHnd, qmin, qmax, iXj, 1,  1. )
    CALL pmaxmin('CH4: Cell Vol',      cellVolume, qmin, qmax, iXj, km, 1. )
    CALL pmaxmin('CH4: Cell Depth',     cellDepth, qmin, qmax, iXj, km, 1. )
    CALL pmaxmin('CH4: Wet air nd',         ndWet, qmin, qmax, iXj, km, 1. )
    CALL pmaxmin('CH4: Dry air nd',         ndDry, qmin, qmax, iXj, km, 1. )
   END IF

!  Convert methane from mole fraction to number density
!  ----------------------------------------------------
   w_c%qa(nbeg)%data3d(i1:i2,j1:j2,1:km) = &
                  w_c%qa(nbeg)%data3d(i1:i2,j1:j2,1:km)*ndWet(i1:i2,j1:j2,1:km)

!  Loss rate [m^{3} s^{-1}] for OH + CH4 => CH3 + H2O
!  --------------------------------------------------
   rkoh(i1:i2,j1:j2,1:km) = 2.45e-12*1.00E-06*exp(-1775./T(i1:i2,j1:j2,1:km))

!  Change in CH4 number density [m^{-3} s^{-1}] due to oxydation
!  -------------------------------------------------------------
   dCH4nd(i1:i2,j1:j2,1:km) = rkoh(i1:i2,j1:j2,1:km)*gcCH4%OHnd(i1:i2,j1:j2,1:km)*w_c%qa(nbeg)%data3d(i1:i2,j1:j2,1:km)

!  Increment the CH4 number density 
!  --------------------------------
   w_c%qa(nbeg)%data3d(i1:i2,j1:j2,1:km) = w_c%qa(nbeg)%data3d(i1:i2,j1:j2,1:km)-cdt*dCH4nd(i1:i2,j1:j2,1:km)

!  Calculate photolytic loss rates, J [s^{-1}], for CH4 + hv => 2H2O + CO.
!  Notice that J and the losses are always computed. However, the setting 
!  of the feedback switch(es) determines if the increments are actually applied.
!  -----------------------------------------------------------------------------
   ALLOCATE(photJ(i1:i2,j1:j2,1:km), STAT=status)
   VERIFY_(STATUS)
   photJ(:,:,:) = 0.00

   CALL getJRates(status)
   VERIFY_(status)

!  Change in CH4 number density [m^{-3} s^{-1}] due to photolysis
!  --------------------------------------------------------------
   ALLOCATE(dCH4Phot(i1:i2,j1:j2,1:km), STAT=status)
   VERIFY_(STATUS)

   dCH4Phot(i1:i2,j1:j2,1:km) = photJ(i1:i2,j1:j2,1:km)*w_c%qa(nbeg)%data3d(i1:i2,j1:j2,1:km)

   IF(ASSOCIATED(CH4_phot)) THEN
    CH4_phot(i1:i2,j1:j2,1:km) = dCH4Phot(i1:i2,j1:j2,1:km)
   END IF

   DEALLOCATE(photJ, STAT=status)
   VERIFY_(STATUS)

!  Increment the CH4 number density when the switch is on
!  ------------------------------------------------------
   IF(gcCH4%CH4FeedBack) &
    w_c%qa(nbeg)%data3d(i1:i2,j1:j2,1:km) = w_c%qa(nbeg)%data3d(i1:i2,j1:j2,1:km)-cdt*dCH4Phot(i1:i2,j1:j2,1:km)

!  Vertically integrated CH4 loss due to oxydation (only)
!  ------------------------------------------------------
   IF(ASSOCIATED(CH4_loss)) THEN
    CH4_loss(i1:i2,j1:j2) = 0.
    DO k = 1, km
     CH4_loss(i1:i2,j1:j2) = CH4_loss(i1:i2,j1:j2) &
   	 + w_c%qa(nbeg)%data3d(i1:i2,j1:j2,k)*rkoh(i1:i2,j1:j2,k) &
   	 * gcCH4%OHnd(i1:i2,j1:j2,k)/ndWet(i1:i2,j1:j2,k) &
   	 * mwtCH4/MAPL_AIRMW*w_c%delp(i1:i2,j1:j2,k)/MAPL_GRAV
    END DO
   END IF

!  CH4 production (None)
!  ---------------------
   IF(ASSOCIATED(CH4_prod)) CH4_prod(i1:i2,j1:j2) = 0.00

!  CH4 Emission: kg cell^{-1} s^{-1}.  Convert to m^{-3} s^{-1}, 
!  multiply by dt, and add to the number density.  Note: No need
!  for regionMask when using a B. Duncan GMI emission dataset.
!  -------------------------------------------------------------
   w_c%qa(nbeg)%data3d(i1:i2,j1:j2,km) = w_c%qa(nbeg)%data3d(i1:i2,j1:j2,km) + cdt*MAPL_AVOGAD*gcCH4%CH4sfcFlux(i1:i2,j1:j2)/ &
                                         (mwtCH4 * cellVolume(i1:i2,j1:j2,km) )

!  Column burden [kg m^{-2}]
!  -------------------------
   IF(ASSOCIATED(CH4_column)) THEN
    CH4_column(i1:i2,j1:j2) = 0.00
    DO k = 1, km
     CH4_column(i1:i2,j1:j2) = CH4_column(i1:i2,j1:j2) + w_c%qa(nbeg)%data3d(i1:i2,j1:j2,k)* &
                               cellDepth(i1:i2,j1:j2,k)*mwtCH4/MAPL_AVOGAD
    END DO
   END IF

!  Fill export state for CH4 mole fraction in dry air
!  --------------------------------------------------
   IF(ASSOCIATED(CH4_dry)) THEN
    CH4_dry(i1:i2,j1:j2,1:km) = w_c%qa(nbeg)%data3d(i1:i2,j1:j2,1:km)/ndDry(i1:i2,j1:j2,1:km)
   END IF

!  Return internal state CH4 to moist-air mole fraction
!  ----------------------------------------------------
   w_c%qa(nbeg)%data3d(i1:i2,j1:j2,1:km) = w_c%qa(nbeg)%data3d(i1:i2,j1:j2,1:km)/ndWet(i1:i2,j1:j2,1:km)

!  Surface concentration in ppmv
!  -----------------------------
   IF(ASSOCIATED(CH4_surface)) &
     CH4_surface(i1:i2,j1:j2) = w_c%qa(nbeg)%data3d(i1:i2,j1:j2,km)*1.00E+06

!  CH4 surface flux [kg m^{-2} s^{-1}] diagnostic
!  ----------------------------------------------
   IF(ASSOCIATED(CH4_emis)) &
        CH4_emis(i1:i2,j1:j2) = gcCH4%CH4sfcFlux(i1:i2,j1:j2)/cellArea(i1:i2,j1:j2)

!  Allow changes to water vapor when the switch is on
!  --------------------------------------------------
   ChangeH2O: IF(gcCH4%CH4FeedBack .AND. gcCH4%H2OFeedBack) THEN

!  Convert water vapor to number density
!  -------------------------------------
    Q(i1:i2,j1:j2,1:km) = Q(i1:i2,j1:j2,1:km)*rhoWet(i1:i2,j1:j2,1:km)*MAPL_AVOGAD/MAPL_H2OMW

!  Increment the water vapor number density by
!  adding two molecules for each CH4 molecule lost.
!  ------------------------------------------------
    Q(i1:i2,j1:j2,1:km) = Q(i1:i2,j1:j2,1:km)+2.00*cdt*dCH4Phot(i1:i2,j1:j2,1:km)

!  Convert water vapor back to mass mixing ratio
!  ---------------------------------------------
    Q(i1:i2,j1:j2,1:km) = Q(i1:i2,j1:j2,1:km)*MAPL_H2OMW/(rhoWet(i1:i2,j1:j2,1:km)*MAPL_AVOGAD)

   END IF ChangeH2O

!  Water vapor tendency [kg kg^{-1} s^{-1}]
!  ----------------------------------------
   IF(ASSOCIATED(CH4_qprod)) THEN
    CH4_qprod(i1:i2,j1:j2,1:km) = 2.00*dCH4Phot(i1:i2,j1:j2,1:km)*MAPL_H2OMW/(ndWet(i1:i2,j1:j2,1:km)*MAPL_AIRMW)
   END IF

   IF(gcCH4%DebugIsOn) THEN
     IF(ASSOCIATED(CH4_emis))    CALL pmaxmin(    'CH4: emis', CH4_emis,    qmin, qmax, iXj,  1, 1. )
     IF(ASSOCIATED(CH4_loss))    CALL pmaxmin(    'CH4: loss', CH4_loss,    qmin, qmax, iXj,  1, 1. )
     IF(ASSOCIATED(CH4_prod))    CALL pmaxmin(    'CH4: prod', CH4_prod,    qmin, qmax, iXj,  1, 1. )
     IF(ASSOCIATED(CH4_column))  CALL pmaxmin(  'CH4: column', CH4_column,  qmin, qmax, iXj,  1, 1. )
     IF(ASSOCIATED(CH4_surface)) CALL pmaxmin( 'CH4: surface', CH4_surface, qmin, qmax, iXj,  1, 1. )
     IF(ASSOCIATED(CH4_qprod))   CALL pmaxmin('CH4: H2O_prod', CH4_qprod,   qmin, qmax, iXj, km, 1. )
     IF(ASSOCIATED(CH4_phot))    CALL pmaxmin('CH4: dCH4phot', dCH4phot,    qmin, qmax, iXj, km, 1. )
     IF(ASSOCIATED(CH4_dry))     CALL pmaxmin(     'CH4: dry', CH4_dry,     qmin, qmax, iXj, km, 1. )
   END IF

!  Housekeeping
!  ------------
   DEALLOCATE(ndDry, ndWet, p, rkoh, dCH4nd, cellDepth, cellVolume, STAT=status)
   VERIFY_(status)
   DEALLOCATE(dCH4Phot, STAT=status)
   VERIFY_(status)

   RETURN

 CONTAINS
 
 SUBROUTINE getJRates(rc)

! Borrowed from meso_phot.F of StratChem, where number densities are cgs [cm^{-3}]
! --------------------------------------------------------------------------------
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
  REAL :: r, radToDeg
  REAL :: s

  INTEGER :: status

  CHARACTER(LEN=ESMF_MAXSTR) :: Iam = "CH4::getJRates"

  radToDeg  = 180.00/MAPL_PI
  rc = 0
  b = SQRT(0.50*r0/hbar)

! O2 overhead number density profile [cm^{-2}]
! --------------------------------------------
  ALLOCATE(O2Column(i1:i2,j1:j2,1:km), STAT=status)
  VERIFY_(status)

  f = O2VMR*5.00E-05
  O2Column(:,:,1) = O2Abv80km+cellDepth(:,:,1)*ndWet(:,:,1)*f

  DO k = 2,km
   O2Column(:,:,k) = O2Column(:,:,k-1)+(cellDepth(:,:,k-1)*ndWet(:,:,k-1)+ &
                                        cellDepth(:,:,  k)*ndWet(:,:,  k))*f
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

   CHARACTER(LEN=*), PARAMETER :: Iam = 'CH4_GridCompFinalize1_'
   INTEGER :: status
   rc = 0

   DEALLOCATE ( gcCH4%CH4sfcFlux, gcCH4%regionMask, gcCH4%OHnd, STAT=status )
   rc = status
   VERIFY_(status)

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

  Use CH4_GridCompMod
  Use ESMF
  Use MAPL
  Use Chem_Mod 

  IMPLICIT NONE

! !INPUT PARAMETERS:

!  Input "function pointer"
!  -----------------------
   interface 
     subroutine Method_ (gc, w, imp, exp, ymd, hms, dt, rcode )
       Use CH4_GridCompMod
       Use ESMF
       Use MAPL
       Use Chem_Mod 
       type(CH4_GridComp1),  intent(inout)  :: gc
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
   INTEGER, INTENT(IN) :: nymd, nhms	       ! time
   REAL,    INTENT(IN) :: cdt		       ! chemical timestep (secs)


! !OUTPUT PARAMETERS:

   TYPE(CH4_GridComp1), INTENT(INOUT) :: gcCH4    ! Grid Component
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
