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
!  14 Jul 2017 Manyin:   Merged ECCOH functionality into Icarus
!  11 Aug 2017 Manyin:   Revert to using TV from DYN
!
!EOP
!-------------------------------------------------------------------------

  TYPE CH4_GridComp1

   CHARACTER(LEN=ESMF_MAXSTR) :: name            ! generic name of the package
   CHARACTER(LEN=ESMF_MAXSTR) :: iname           ! instance name
   CHARACTER(LEN=ESMF_MAXSTR) :: rcfilen         ! resource file name
   CHARACTER(LEN=ESMF_MAXSTR) :: regionsString   ! Comma-delimited string of regions

   INTEGER :: instance                 ! Instantiation number

!  REAL, POINTER :: regionMask(:,:)    ! regional mask
   REAL, POINTER :: CH4sfcFlux(:,:)    ! CH4 surface flux [kg m^-2 s^-1]  (do not alloc)

   LOGICAL :: DebugIsOn     ! Run-time debug switch

   INTEGER :: C_isotope     ! 12, 13 or 0,  set based on the instance name
   REAL*8  :: mwtCH4        ! set based on isotopologue

!  REAL :: szaCutoff        ! Largest solar zenith angle (degrees) allowed as daytime

  END TYPE CH4_GridComp1

  TYPE CH4_GridComp
     INTEGER                      ::  n        ! number of instances 
     TYPE(CH4_GridComp1), pointer ::  gcs(:)   ! instances
  END TYPE CH4_GridComp

CONTAINS

!------------------------------------------------------------------------
!     NASA/GSFC, Atmospheric Chemistry and Dynamics Lab, Code 614       !
!------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  CH4_GridCompSetServices --- Set Services for CH4_GridComp
!
! !INTERFACE:
!
   subroutine CH4_GridCompSetServices( gc, chemReg, rc)

! !USES:

  IMPLICIT NONE

! !INPUT PARAMETERS:

! !OUTPUT PARAMETERS:

   type(ESMF_GridComp), intent(INOUT) :: gc
   type(Chem_Registry), intent(INOUT) :: chemReg
   integer,             intent(OUT  ) :: rc

! !DESCRIPTION: Set Services for the CH4 Grid Component.
!
! !REVISION HISTORY:
!
!
!EOP
!-------------------------------------------------------------------------


   CHARACTER(LEN=255) :: rcbasen = 'CH4_GridComp'
   CHARACTER(LEN=255) :: name

   integer            :: status
   character(len=ESMF_MAXSTR) :: Iam

   integer :: ier,n,i

   type(ESMF_Config) :: cfg

   Iam = "CH4_GridCompSetServices"

!  Load resource file
!  ------------------
   cfg = ESMF_ConfigCreate(__RC__)

   call ESMF_ConfigLoadFile(cfg,TRIM(rcbasen)//'.rc',__RC__)

!  Parse resource file
!  -------------------
   n = ESMF_ConfigGetLen(cfg,label='CH4_instances:',__RC__)

!  We cannot have fewer instances than the number of
!   CH4 bins in the registry
!  --------------------------------------------------------------------
   if ( n .LT. chemReg%n_CH4 ) then
        if (MAPL_AM_I_ROOT()) &
        PRINT *, TRIM(Iam)//": HALTING - Fewer instances (",n,") than Registry bins (",chemReg%n_CH4,")"
        rc = 35
        return
   else if ( n .GT. chemReg%n_CH4 ) then
        if (MAPL_AM_I_ROOT()) then
          PRINT *, TRIM(Iam)//": More instances (",n,") than Registry bins (",chemReg%n_CH4,")"
          PRINT *, TRIM(Iam)//": Only running the first ",chemReg%n_CH4," instances"
        end if
        n = chemReg%n_CH4
   end if

!  Record name of each instance
!  ----------------------------
   call ESMF_ConfigFindLabel(cfg,'CH4_instances:',__RC__)

   do i = 1, n

      call ESMF_ConfigGetAttribute(cfg,name,__RC__)

      IF(TRIM(name) == "full" ) THEN
       name = " "              ! blank instance name for full (1)
      END IF


      call CH4_GridCompSetServices1_(gc,chemReg,name,__RC__)

   end do

   
   call MAPL_AddImportSpec(GC,             &
        SHORT_NAME = 'CH4_Cl',             &
        LONG_NAME  = 'Cl for methane'  ,   &
        UNITS      = '1',                  &
        DIMS       = MAPL_DimsHorzVert,    &
        VLOCATION  = MAPL_VLocationCenter, &
        __RC__)

   call MAPL_AddImportSpec(GC,             &
        SHORT_NAME = 'CH4_O1D',            &
        LONG_NAME  = 'O1D for methane'  ,  &
        UNITS      = '1',                  &
        DIMS       = MAPL_DimsHorzVert,    &
        VLOCATION  = MAPL_VLocationCenter, &
        __RC__)

   call MAPL_AddImportSpec(GC,             &
        SHORT_NAME = 'CH4_oh',             &
        LONG_NAME  = 'OH for methane'  ,   &
        UNITS      = '1',                  &
        DIMS       = MAPL_DimsHorzVert,    &
        VLOCATION  = MAPL_VLocationCenter, &
        __RC__)

   call MAPL_AddImportSpec(GC,           &
        SHORT_NAME = 'CH4_regionMask',   &
        LONG_NAME  = 'mask'  ,           &
        UNITS      = '1',                &
        DIMS       = MAPL_DimsHorzOnly,  &
        VLOCATION  = MAPL_VLocationNone, &
        RESTART    = MAPL_RestartSkip,   &
        __RC__)

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

   TYPE(Chem_Bundle), intent(in) :: w_c         ! Chemical tracer fields      
   INTEGER, INTENT(IN) :: nymd, nhms            ! time
   REAL,    INTENT(IN) :: cdt                   ! chemical timestep (secs)


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
   
   INTEGER :: i, status

!  Load resource file
!  ------------------
   CALL I90_loadf ( TRIM(rcbasen)//'.rc', status )
   VERIFY_(status)

!  Parse resource file
!  -------------------
   CALL I90_label ( 'CH4_instances:', status )
   VERIFY_(status)

!  We already verified the instance count in the registry
!  ------------------------------------------------------
   gcCH4%n = w_c%reg%n_CH4

!  Next allocate necessary memory
!  ------------------------------
   ALLOCATE ( gcCH4%gcs(gcCH4%n), STAT=status )    
   VERIFY_(status)

!  Record name of each instance
!  ----------------------------
   CALL I90_label ( 'CH4_instances:', status )
   VERIFY_(status)
   DO i = 1, gcCH4%n
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

   TYPE(Chem_Bundle), intent(in) :: w_c         ! Chemical tracer fields      
   INTEGER, INTENT(IN) :: nymd, nhms            ! time
   REAL,    INTENT(IN) :: cdt                   ! chemical timestep (secs)


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
! !IROUTINE:  CH4_GridCompFinalize --- Finalize CH4_GridComp
!
! !INTERFACE:
!

   SUBROUTINE CH4_GridCompFinalize ( gcCH4, w_c, impChem, expChem, nymd, nhms, cdt, rc )

! !USES:

  IMPLICIT NONE

! !INPUT PARAMETERS:

   TYPE(Chem_Bundle), intent(in) :: w_c         ! Chemical tracer fields      
   INTEGER, INTENT(IN) :: nymd, nhms            ! time
   REAL,    INTENT(IN) :: cdt                   ! chemical timestep (secs)


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

!--------------------------------------------------------------------------

!                      Single Instance Methods

 subroutine CH4_GridCompSetServices1_(  gc, chemReg, iname, rc)
 type(ESMF_GridComp), intent(INOUT) :: GC
 type(Chem_Registry), intent(INOUT) :: chemReg
 character(len=*),    intent(IN   ) :: iname
 integer,             intent(OUT  ) :: rc

 integer :: Status
 character(len=ESMF_MAXSTR) :: Iam

 Iam ="CH4_GridCOmpSetServices1_"

  call MAPL_AddImportSpec(GC,             &
       SHORT_NAME = 'CH4_sfcFlux'//iname, &
       LONG_NAME  = 'source species'  ,   &
       UNITS      = '1',                  &
       DIMS       = MAPL_DimsHorzOnly,    &
       VLOCATION  = MAPL_VLocationNone,   &
       RC         = STATUS)
  VERIFY_(STATUS)

 RETURN_(ESMF_SUCCESS)

 end subroutine CH4_GridCompSetServices1_

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  CH4_GridCompInitialize1_ --- Initialize CH4_GridComp
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
   INTEGER :: nbeg, nend
   LOGICAL :: NoRegionalConstraint 
   INTEGER :: length

   rcfilen = gcCH4%rcfilen
   gcCH4%name = 'GEOS-5/GOCART Parameterized CH4 Package'

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
!  ALLOCATE ( gcCH4%regionMask(i1:i2,j1:j2), &
!             __STAT__ )

!  Load resource file
!  ------------------
   CALL I90_loadf ( TRIM(rcfilen), status )
   VERIFY_(status)

!  Maximum allowed solar zenith angle for "daylight"
!  -------------------------------------------------
!  CALL I90_label ( 'solar_ZA_cutoff:', status )
!  VERIFY_(status)
!  gcCH4%szaCutoff = I90_gfloat( status )
!  VERIFY_(status)

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

!  13C vs 12C isotopologue -sas
!  ---------------------------
   ! default values:
   gcCH4%C_isotope = 0
   gcCH4%mwtCH4 = 1.008*4 + 12.011

   ! conditionally overwrite
   length = LEN( TRIM(gcCH4%iname) )
   IF ( length > 3 ) THEN
     IF      ( gcCH4%iname(length-2:length) == "12C" ) THEN
       gcCH4%C_isotope = 12
       gcCH4%mwtCH4 = 1.008*4 + 12
     ELSE IF ( gcCH4%iname(length-2:length) == "13C" ) THEN
       gcCH4%C_isotope = 13
       gcCH4%mwtCH4 = 1.008*4 + 13
     END IF
   END IF

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
    SELECT CASE (lowercase(gcCH4%regionsString(1:2)))
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

!  All done -sas, release to save space
!  --------
   CALL I90_FullRelease( status )
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
   TYPE(Chem_Bundle),   INTENT(INOUT) :: w_c     ! Chemical tracer fields   

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
   REAL, POINTER, DIMENSION(:,:,:) ::  Q        => null()
   REAL, POINTER, DIMENSION(:,:,:) ::  T        => null()
   REAL, POINTER, DIMENSION(:,:,:) ::  rhoa     => null()
   REAL, POINTER, DIMENSION(:,:,:) ::  rhoDry   => null()
   REAL, POINTER, DIMENSION(:,:,:) ::  ple      => null()
   REAL, POINTER, DIMENSION(:,:,:) ::  zle      => null()
!  REAL, POINTER, DIMENSION(:,:,:) ::  CH4_rat  => null()      ! if we need a special export to RADIATION

   REAL, POINTER, DIMENSION(:,:,:) ::  CL_nd    => null()      !  CL number density (molec cm^-3)
   REAL, POINTER, DIMENSION(:,:,:) ::  O1D_nd   => null()      ! O1D number density (molec cm^-3)
   REAL, POINTER, DIMENSION(:,:,:) ::  OH_nd    => null()      !  OH number density (molec cm^-3)

   INTEGER :: i1, i2, im, j1, j2, jm, km, idiag, iXj
   INTEGER :: i, j, k, kReverse, n, nbeg, nend
   INTEGER :: status


   REAL    :: qmin, qmax
   REAL    :: a_OH, a_O1D, a_Cl    !sas, alpha for isotope rate differences

   REAL, ALLOCATABLE :: cellDepth(:,:,:), cellVolume(:,:,:),     rkoh(:,:,:)
   REAL, ALLOCATABLE ::         p(:,:,:),         nd(:,:,:)
   REAL, ALLOCATABLE ::     ndDry(:,:,:)
   REAL, ALLOCATABLE ::        TV(:,:,:)

   REAL, PARAMETER ::     epsilon = (MAPL_H2OMW/MAPL_AIRMW)


!  + EYegorova - updated Apr 27, 2011
   REAL, ALLOCATABLE ::      rkcl(:,:,:)
   REAL, ALLOCATABLE ::  rkohloss(:,:,:),   rkclloss(:,:,:), rko1dloss(:,:,:) 
!  - EYegorova - updated Apr 27, 2011

   REAL              ::  rko1d_scalar

   REAL, POINTER :: ptr2d(:,:)   => null()

#define EXPORT   expChem
#define iNAME    TRIM(gcCH4%iname)

#define CH4EM   CH4_emis
#define CH4CL   CH4_column
#define CH4SC   CH4_surface
#define CH4PD   CH4_prod
#define CH4LS   CH4_loss
#define CH4DRY  CH4_dry

#include "CH4_GetPointer___.h"

! WITH 'full', WE ALWAYS EXPORT 'CH4'
!#define RATSpro_INSTANCE 'tot'
!
!!  To export total CH4 to the Radiation code:
!!  ------------------------------------------
!   IF ( iNAME == RATSpro_INSTANCE ) THEN
!     call MAPL_GetPointer ( expChem, CH4_rat,  'CH4', RC=STATUS )
!     VERIFY_(STATUS)
!   END IF 

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
   CALL MAPL_GetPointer(impChem, Q,        'Q',            __RC__ )
   CALL MAPL_GetPointer(impChem, T,        'T',            __RC__ )
   CALL MAPL_GetPointer(impChem, rhoa,     'AIRDENS',      __RC__ )
   CALL MAPL_GetPointer(impChem, rhoDry,   'AIRDENS_DRYP', __RC__ )
   CALL MAPL_GetPointer(impChem, cellArea, 'AREA',         __RC__ )
   CALL MAPL_GetPointer(impChem, ple,      'PLE',          __RC__ )
   CALL MAPL_GetPointer(impChem, zle,      'ZLE',          __RC__ )

   IF(gcCH4%DebugIsOn) THEN
    CALL pmaxmin('CH4:AREA', cellArea, qmin, qmax, iXj,  1,   1. )
    CALL pmaxmin('CH4:T',           T, qmin, qmax, iXj, km,   1. )
    CALL pmaxmin('CH4:Q',           Q, qmin, qmax, iXj, km,   1. )
    CALL pmaxmin('CH4:rhoa',     rhoa, qmin, qmax, iXj, km,   1. )
    CALL pmaxmin('CH4:rhoDry', rhoDry, qmin, qmax, iXj, km,   1. )
    CALL pmaxmin('CH4:ple',       ple, qmin, qmax, iXj, km+1, 1. )
    CALL pmaxmin('CH4:zle',       zle, qmin, qmax, iXj, km+1, 1. )
   END IF


!   CH4 emissions [kg m^-2 s^-1]
!   ----------------------------
    call MAPL_GetPointer(impChem,ptr2d,'CH4_sfcFlux'//trim(iNAME),__RC__)
!   KEEP as a pointer:
    gcCH4%CH4sfcFlux => ptr2d

!   NOTE:
!   For Cl, O1D, OH   we do not need separate "instance specific" files

!   Cl [molec cm^-3]
!   ----------------
    call MAPL_GetPointer(impChem,  CL_nd, 'CH4_Cl',  __RC__)

!   O1D [molec cm^-3]
!   -----------------
    call MAPL_GetPointer(impChem, O1D_nd, 'CH4_O1D', __RC__)

!   OH
!   May be from ExtData or online OH
!   __

    call MAPL_GetPointer(impChem, OH_nd, 'CH4_oh', __RC__)




!  Allocate temporary workspace
!  ----------------------------
   ALLOCATE( p(i1:i2,j1:j2,km), &
            nd(i1:i2,j1:j2,km), &
            TV(i1:i2,j1:j2,km), &
         ndDry(i1:i2,j1:j2,km), &
          rkoh(i1:i2,j1:j2,km), &
          rkcl(i1:i2,j1:j2,km), &
      rkohloss(i1:i2,j1:j2,km), &
      rkclloss(i1:i2,j1:j2,km), &
     rko1dloss(i1:i2,j1:j2,km), &
    cellVolume(i1:i2,j1:j2,km), &
     cellDepth(i1:i2,j1:j2,km), STAT=status)
   VERIFY_(status)

!  Layer mean pressures (top-down)
!  -------------------------------
   DO k=1,km
    p(i1:i2,j1:j2,k) = (ple(i1:i2,j1:j2,k-1)+ple(i1:i2,j1:j2,k))*0.50
   END DO
 
   !  Virtual Temperature (K)
   TV = T * (1.0 + Q/MAPL_EPSILON)/(1.0 + Q)

!  Moist air number density  (molec/m3)
!  ------------------------------------
   nd = (MAPL_AVOGAD * p) / (MAPL_RUNIV * TV)

!    Dry air number density  (molec/m3)
!    Start with the mass of only the Dry Air in the gridbox...
!    (Note that this is not consistent with fields like T and P!)
!  ------------------------
   ndDry(:,:,1:km) = rhoDry(:,:,1:km)*MAPL_AVOGAD/MAPL_AIRMW



!  Cell depth and volume
!  ---------------------
   DO k=1,km
!   MANYIN - previous approach:
!   cellDepth(:,:,k) = (ple(:,:,k)-ple(:,:,k-1))/(rhoa(:,:,k)*MAPL_GRAV)
    cellDepth(:,:,k) = zle(:,:,k-1)-zle(:,:,k)
    cellVolume(:,:,k) = cellArea(:,:)*cellDepth(:,:,k)
   END DO

   IF(gcCH4%DebugIsOn) THEN !-sas
      CALL pmaxmin('surf con mr1',w_c%qa(nbeg)%data3d(i1:i2,j1:j2,km) , qmin, qmax, iXj, 1,  1. )
   ENDIF

!  Convert internal state CH4 from volume mixing ratio wrt moist air to number density (molec/m3)
!  ----------------------------------------------------------------------------------------------
   w_c%qa(nbeg)%data3d(:,:,:) = &
   w_c%qa(nbeg)%data3d(:,:,:) * nd(:,:,:)

   IF(gcCH4%DebugIsOn) THEN !-sas
      CALL pmaxmin('surf con nd1',w_c%qa(nbeg)%data3d(i1:i2,j1:j2,km) , qmin, qmax, iXj, 1,  1. )
   ENDIF

!  Debug
!  -----
   IF(gcCH4%DebugIsOn) THEN
    CALL pmaxmin('CH4: SfcFlux', gcCH4%CH4sfcFlux, qmin, qmax, iXj, 1,  1. )
    CALL pmaxmin('CH4: Cell Vol',      cellVolume, qmin, qmax, iXj, km, 1. )
    CALL pmaxmin('CH4: Cell Depth',     cellDepth, qmin, qmax, iXj, km, 1. )
    CALL pmaxmin('CH4: OH Conc',            OH_nd, qmin, qmax, iXj, km, 1. )
   END IF

   ! adjust rate constants for 13C- vs. 12C-CH4 -sas
   IF (gcCH4%C_isotope == 13) THEN
      a_OH  = 0.9946
      a_O1D = 0.987
      a_Cl  = 0.938
   ELSE
      ! Values for both 12C and global avg C
      a_OH  = 1.0
      a_O1D = 1.0
      a_Cl  = 1.0
   END IF

!  Loss rates [m^{-3} s^{-1}]

!  CH4 + OH => CH3 + H2O
!  --------------------------------------------------
   rkoh(:,:,:) = 2.45e-12*1.00E-06*exp(-1775./T(:,:,:)) * a_OH

!  CH4 + O1D =>  products
!  --------------------------------------------------
   rko1d_scalar = 1.75e-10*1.00E-06 * a_O1D   ! MANYIN  just a constant

!  CH4 + Cl => HCl + MO2
!  --------------------------------------------------
   rkcl(:,:,:) =  7.3e-12*1.00E-06*exp(-1280./T(:,:,:)) * a_Cl


!  Compute loss due to OH (molec/m3)
!  (convert OH number density from molecules cm^-3 to molecules m^-3)
   rkohloss = cdt*rkoh*(OH_nd*1.00E+06)

!  Compute loss due to O1D (molec/m3)
!  (Convert O1D_nd from molecules cm^-3 to molecules m^-3)
   rko1dloss = cdt*rko1d_scalar*(O1D_nd * 1.00E+06)

!  Compute loss due to Cl (molec/m3)
!  (Convert CL_nd from molecules cm^-3 to molecules m^-3)
   rkclloss = cdt*rkcl*(CL_nd * 1.00E+06)


!  jsw: modified the original version in the following to account for
!  CH4 concentration, to have units of molec/m3/s rather than molec/m3 to  
!  facilitate time averaging, and to contain all layers rather than just lowest.
   IF(ASSOCIATED(CH4_loss)) CH4_loss(:,:,:) = w_c%qa(nbeg)%data3d(:,:,:)*(rkohloss(:,:,:)+rko1dloss(:,:,:)+rkclloss(:,:,:))/cdt

   w_c%qa(nbeg)%data3d(:,:,:) = &
   w_c%qa(nbeg)%data3d(:,:,:) * (1.00-rkohloss(:,:,:)-rko1dloss(:,:,:)-rkclloss(:,:,:))


!  CH4 production (None)
!  ---------------------
   IF(ASSOCIATED(CH4_prod)) CH4_prod(:,:) = 0.00


!  CH4 Emission -  Flux: kg m^-2 s^-1.  Convert to molec m^-3 
!  -----------------------------------------------------------
   w_c%qa(nbeg)%data3d(:,:,km) = &
   w_c%qa(nbeg)%data3d(:,:,km) + (cdt*MAPL_AVOGAD* (gcCH4%CH4sfcFlux(:,:)/gcCH4%mwtCH4) /cellDepth(:,:,km))


!  Column burden [kg m^{-2}]
!  -------------------------
   IF(ASSOCIATED(CH4_column)) THEN
    CH4_column(:,:) = 0.00
    DO k = 1, km
     CH4_column(:,:) = CH4_column(:,:) + w_c%qa(nbeg)%data3d(:,:,k)* &
                               cellDepth(:,:,k)*gcCH4%mwtCH4/MAPL_AVOGAD
    END DO
   END IF

!  Fill export state for CH4 mole fraction in dry air
!  --------------------------------------------------
   IF(ASSOCIATED(CH4_dry)) THEN
    CH4_dry(:,:,:) = w_c%qa(nbeg)%data3d(:,:,:)/ndDry(:,:,:)   ! MANYIN - dry?
   END IF


!!  CH4 for the radiation code, stored in parts-per-part by volume (wrt dry air)
!!  ----------------------------------------------------------------------------
!   IF( iNAME==RATSpro_INSTANCE .AND. ASSOCIATED(CH4_rat)) THEN
!     IF(MAPL_AM_I_ROOT()) PRINT *, 'Provide GOCART::CH4 to Radiation (RATs provider)'
!     CH4_rat(:,:,:) = w_c%qa(nbeg)%data3d(:,:,:)/ndDry(:,:,:)
!   ENDIF

   IF(gcCH4%DebugIsOn) THEN !-sas
      CALL pmaxmin('surf con nd2',w_c%qa(nbeg)%data3d(i1:i2,j1:j2,km) , qmin, qmax, iXj, 1,  1. )
   ENDIF

!  Return internal state CH4 from number density (molec/m3) to volume mixing ratio wrt moist air
!  ---------------------------------------------------------------------------------------------
   w_c%qa(nbeg)%data3d(:,:,:) = &
   w_c%qa(nbeg)%data3d(:,:,:) / nd(:,:,:)

   IF(gcCH4%DebugIsOn) THEN !-sas
      CALL pmaxmin('surf con mr2',w_c%qa(nbeg)%data3d(i1:i2,j1:j2,km) , qmin, qmax, iXj, 1,  1. )
   ENDIF

!  Surface concentration in ppbv (wrt moist air)
!  ---------------------------------------------
   IF(ASSOCIATED(CH4_surface)) &
                 CH4_surface(:,:) = w_c%qa(nbeg)%data3d(:,:,km)*1.00E+09

!  CH4 surface flux diagnostic [kg m-2 s-1]
!  -----------------------------------------------------------------------
   IF(ASSOCIATED(CH4_emis)) &
                 CH4_emis(:,:) = gcCH4%CH4sfcFlux(:,:)

   IF(gcCH4%DebugIsOn) THEN
     IF(ASSOCIATED(CH4_emis))    CALL pmaxmin( 'CH4: emis',     CH4_emis,    qmin, qmax, iXj,  1, 1. )
     IF(ASSOCIATED(CH4_loss))    CALL pmaxmin( 'CH4: loss',     CH4_loss,    qmin, qmax, iXj,  1, 1. )
     IF(ASSOCIATED(CH4_prod))    CALL pmaxmin( 'CH4: prod',     CH4_prod,    qmin, qmax, iXj,  1, 1. )
     IF(ASSOCIATED(CH4_column))  CALL pmaxmin( 'CH4: column',   CH4_column,  qmin, qmax, iXj,  1, 1. )
     IF(ASSOCIATED(CH4_surface)) CALL pmaxmin( 'CH4: surface',  CH4_surface, qmin, qmax, iXj,  1, 1. )
     IF(ASSOCIATED(CH4_dry))     CALL pmaxmin( 'CH4: dry',      CH4_dry,     qmin, qmax, iXj, km, 1. )
   END IF

!  Housekeeping
!  ------------
   DEALLOCATE( p, TV, nd, ndDry, rkoh, rkcl, rkohloss, rkclloss, rko1dloss, cellVolume, cellDepth, __STAT__)


   RETURN

 END SUBROUTINE CH4_GridCompRun1_

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  CH4_GridCompFinalize1_ --- Finalize CH4_GridComp
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

!  DEALLOCATE ( gcCH4%regionMask, STAT=status )
!  rc = status
!  VERIFY_(status)

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
  Use MAPL_Mod
  Use Chem_Mod 

  IMPLICIT NONE

! !INPUT PARAMETERS:

!  Input "function pointer"
!  -----------------------
   interface 
     subroutine Method_ (gc, w, imp, exp, ymd, hms, dt, rcode )
       Use CH4_GridCompMod
       Use ESMF
       Use MAPL_Mod
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
   INTEGER, INTENT(IN) :: nymd, nhms           ! time
   REAL,    INTENT(IN) :: cdt                  ! chemical timestep (secs)


! !OUTPUT PARAMETERS:

   TYPE(CH4_GridComp1), INTENT(INOUT) :: gcCH4    ! Grid Component
   TYPE(ESMF_State),    INTENT(INOUT) :: impChem  ! Import State
   TYPE(ESMF_State),    INTENT(INOUT) :: expChem  ! Export State
   INTEGER,             INTENT(OUT)   :: rc       ! Error return code:
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
