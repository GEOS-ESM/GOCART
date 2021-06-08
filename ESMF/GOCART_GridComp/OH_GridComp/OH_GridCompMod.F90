#include "MAPL_Generic.h"

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  OH_GridCompMod --- OH Parameterization Grid Component Class
!
! !INTERFACE:
!

   MODULE  OH_GridCompMod

! !USES:

   USE ESMF
   USE MAPL

   USE Chem_Mod              ! Chemistry Base Class
   USE Chem_StateMod         ! Chemistry State
   USE Chem_UtilMod          ! I/O

   USE m_inpak90             ! Resource file management
   USE m_die, ONLY: die
   USE ohParameterizationMethod_mod
   USE solarZenithAngle_mod
   USE computeReflectivity_mod

   USE cblock_size_mod         ! "CMN_SIZE"       ! Size parameters
   USE cblock_CO_mod           ! "CMN_CO"         ! CO arrays
   USE cblock_CO_budget_mod    ! "CMN_CO_BUDGET"  ! FMOL_CO
   USE cblock_OH_mod           ! "CMN_OH"

   IMPLICIT NONE

! !PUBLIC TYPES:
!
   PRIVATE
   PUBLIC  OH_GridComp       ! OH object 

!
! !PUBLIC MEMBER FUNCTIONS:
!
   PUBLIC  OH_GridCompSetServices
   PUBLIC  OH_GridCompInitialize
   PUBLIC  OH_GridCompRun
   PUBLIC  OH_GridCompFinalize

!
!
! !DESCRIPTION:
!
!  This module implements the OH Parameterization Grid Component. 
!
! !REVISION HISTORY:
!
!  06September2011 Kouatchou  First crack.
!  11 Aug 2017 Manyin:   Revert to using TV from DYN
!
!EOP
!-------------------------------------------------------------------------

! MANYIN- Notes (things to change/test in future)
!   LPAUSE should be 36
!   o3up should be 0.0 at level km
!   i_OH+6 (DryAirMass_24h) never gets used in the current GETINFO module
!   Current year is always used in Acquire_GMIspecies, but
!           2001 is always used in Acquire_OHratios

  TYPE OH_GridComp

       CHARACTER(LEN=255) :: name            ! generic name of the package
       CHARACTER(LEN=255) :: iname           ! instance name
       CHARACTER(LEN=255) :: rcfilen         ! resource file name
       CHARACTER(LEN=255) :: OHFileName      ! OH concentration file name
!      CHARACTER(LEN=255) :: GMI_infile_name ! GMI file name
       CHARACTER(LEN=255) :: Qscl_infile_name ! Q scalar file name

       INTEGER :: min_time                       ! smallest time >= 00Z that the Run routine will see  (sec)

       LOGICAL :: DBG

       REAL*8, POINTER ::       Temp_24h(:,:,:)  ! 24h average temperature                               bottom-up
       REAL*8, POINTER ::      Press_24h(:,:,:)  ! 24h average pressure                                  bottom-up
                                                 ! MANYIN: changed from edge to midpoint pressure 6.9.17
     ! REAL*8, POINTER ::      Tropp_24h(:,:)    ! 24h average tropopause pressure
     ! REAL*8, POINTER ::  CloudFlux_24h(:,:,:)  ! 24h average cloud flux
       REAL*8, POINTER ::    SpecHum_24h(:,:,:)  ! 24h average specific humidity                         bottom-up
     ! REAL*8, POINTER :: DryAirMass_24h(:,:,:)  ! 24h average dry air mass                                        top-down
     ! REAL*8, POINTER ::   optDepth_24h(:,:,:)  ! 24h average optical depth

    ! REAL*8, POINTER ::        TroppAcc(:,:)    ! Accumulation of tropopause pressure
    ! REAL*8, POINTER ::    CloudFluxAcc(:,:,:)  ! Accumulation of cloud flux
    ! REAL*8, POINTER ::     optDepthAcc(:,:,:)  ! Accumulation of optical depth

      REAL*8, POINTER :: Ravga(:,:,:)            !                                                       bottom-up
      REAL*8, POINTER :: Ravgb(:,:,:)            !                                                       bottom-up

      INTEGER, POINTER :: LPAUSE(:,:)            ! tropopause level

      REAL*8, POINTER :: troppresIn(:,:)         ! tropopause pressure read from file
      REAL*8, POINTER :: qIn(:,:,:)              ! GMI water vapor pressure read from file        unused
      REAL*8, POINTER :: CH4_24h(:,:,:)          ! 24h average of CH4                                    bottom-up
      REAL*8, POINTER :: CO_24h (:,:,:)          ! 24h average of CO                                     bottom-up
      REAL, ALLOCATABLE :: rxnRateRatio(:,:,:)   ! ratio to adjust for updated  reaction rate constants  bottom-up


      type(t_OHparam) :: OHparam
  END TYPE OH_GridComp

  REAL,    PARAMETER :: mwtWat = MAPL_H2OMW
  REAL,    PARAMETER :: mwtCO  = 28.01
  REAL,    PARAMETER :: mwtCH4 = 16.043
  REAL,    PARAMETER :: Pa2hPa = 0.01
  REAL,    PARAMETER :: secPerDay = 86400.00
  INTEGER, PARAMETER :: negativeOne = -1

! GMI species names for the BBIJ array:
! MANYIN-- note: only NFIELDS2 (12) species are currently allocated
! MANYIN         see cblock_CO_mod
!
!   ALK4  = alkanes
!   ISOP  = Isoprene
!   ACET  = Acetone
!   PRPE  = Propene
!   C3H8  = Propane
!   C2H6  = Ethane
  CHARACTER(LEN=20) :: gmiSpeciesName3D(17) = (/'OH_NO      ', 'OH_NO2     ', 'OH_N2O5    ', 'OH_NO3     ', 'OH_HNO2    ',   &
                                                'OH_HNO4    ', 'OH_ALK4    ', 'OH_ISOP    ', 'OH_ACET    ', 'OH_PRPE    ',   &
                                                'OH_C3H8    ', 'OH_C2H6    ', 'OH_O3      ', 'OH_CO      ', 'OH_CH4     ',   &
                                                'OH_OH      ', 'OH_METWATER'/)
  CHARACTER(LEN=10) :: gmiSpeciesName2D( 1) = (/'OH_TRPPS'/)

  INTEGER, PARAMETER :: IC_TRPPS = 1


!-------------------------------------------------------------------------
CONTAINS
!------------------------------------------------------------------------
!     NASA/GSFC, Atmospheric Chemistry and Dynamics Lab, Code 614       !
!------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  OH_GridCompSetServices --- Set Services for OH_GridComp
!
! !INTERFACE:
!
   subroutine OH_GridCompSetServices( gc, chemReg, rc)

! !USES:

  IMPLICIT NONE

! !INPUT PARAMETERS:

! !OUTPUT PARAMETERS:

   type(ESMF_GridComp), intent(INOUT) :: gc
   type(Chem_Registry), intent(INOUT) :: chemReg
   integer,             intent(OUT  ) :: rc

! !DESCRIPTION: Set Services for the OH Grid Component.
!
! !REVISION HISTORY:
!
!
!EOP
!-------------------------------------------------------------------------


   CHARACTER(LEN=255) :: rcbasen = 'OH_GridComp'
   CHARACTER(LEN=255) :: name

   integer            :: status
   character(len=ESMF_MAXSTR) :: Iam

   integer :: ier,n,i

   type(ESMF_Config) :: cfg

   Iam = "OH_GridCompSetServices"

!!!
!!!  There are currently 7 instances in the Registry, but only the first is OH...
!!!


!!  Load resource file
!!  ------------------
!   cfg = ESMF_ConfigCreate(__RC__)
!
!   call ESMF_ConfigLoadFile(cfg,TRIM(rcbasen)//'.rc',__RC__)
!
!!  Parse resource file
!!  -------------------
!   n = ESMF_ConfigGetLen(cfg,label='OH_instances:',__RC__)

!  MANYIN - We enforce the number of "instances"
!  ---------------------------------------------
   if ( chemReg%n_OH /= 1 ) THEN
        IF (MAPL_AM_I_ROOT()) THEN
          PRINT *, TRIM(Iam)//": Registry must list 1 bin"
          PRINT *, TRIM(Iam)//": Registry bins = ",chemReg%n_OH
         END IF
        VERIFY_(35)
   end if

!!  Record name of each instance
!!  ----------------------------
!   call ESMF_ConfigFindLabel(cfg,'OH_instances:',__RC__)
!
!   do i = 1, n
!
!      call ESMF_ConfigGetAttribute(cfg,name,__RC__)
!
!!  We choose not to use the 'full' instance approach because
!!  this would mean that the full field, stored in mass mixing ratio,
!!  would be exported as 'CH4'.  But that name is needed for
!!  export to RADIATION, in volume mixing ratio.
!!     IF(TRIM(name) == "full" ) THEN
!!      name = " "              ! blank instance name for full (1)
!!     ELSE
!!      name = TRIM(name)       ! instance name for others
!!     END IF
!
!
!      call OH_GridCompSetServices1_(gc,chemReg,name,__RC__)
!
!   end do


!      Preserve between runs, in the Internal State
!      Do not transport
!      --------------------------------------------

       call MAPL_AddInternalSpec(GC,                        &
          SHORT_NAME = 'GOCART::CO_24h',                    &
          LONG_NAME  = '24h avg CO for OH',                 &
          UNITS      = 'mol/mol',                           &
          FRIENDLYTO = ' ',                                 &
          RESTART    = MAPL_RestartOptional,                &
          DEFAULT    = 0.0,                                 &
          DIMS       = MAPL_DimsHorzVert,                   &
          VLOCATION  = MAPL_VLocationCenter, __RC__)

       call MAPL_AddInternalSpec(GC,                        &
          SHORT_NAME = 'GOCART::CH4_24h',                   &
          LONG_NAME  = '24h avg CH4 for OH',                &
          UNITS      = 'mol/mol',                           &
          FRIENDLYTO = ' ',                                 &
          RESTART    = MAPL_RestartOptional,                &
          DEFAULT    = 0.0,                                 &
          DIMS       = MAPL_DimsHorzVert,                   &
          VLOCATION  = MAPL_VLocationCenter, __RC__)

       call MAPL_AddInternalSpec(GC,                        &
          SHORT_NAME = 'GOCART::Temp_24h',                  &
          LONG_NAME  = '24h avg Temperature for OH',        &
          UNITS      = 'K',                                 &
          FRIENDLYTO = ' ',                                 &
          RESTART    = MAPL_RestartOptional,                &
          DEFAULT    = 0.0,                                 &
          DIMS       = MAPL_DimsHorzVert,                   &
          VLOCATION  = MAPL_VLocationCenter, __RC__)

       call MAPL_AddInternalSpec(GC,                        &
          SHORT_NAME = 'GOCART::Press_24h',                 &
          LONG_NAME  = '24h avg Pressure for OH',           &
          UNITS      = 'Pa',                                &
          FRIENDLYTO = ' ',                                 &
          RESTART    = MAPL_RestartOptional,                &
          DEFAULT    = 0.0,                                 &
          DIMS       = MAPL_DimsHorzVert,                   &
          VLOCATION  = MAPL_VLocationCenter, __RC__)

       call MAPL_AddInternalSpec(GC,                        &
          SHORT_NAME = 'GOCART::SpecHum_24h',               &
          LONG_NAME  = '24h avg Specific Humidity for OH',  &
          UNITS      = '1',                                 &
          FRIENDLYTO = ' ',                                 &
          RESTART    = MAPL_RestartOptional,                &
          DEFAULT    = 0.0,                                 &
          DIMS       = MAPL_DimsHorzVert,                   &
          VLOCATION  = MAPL_VLocationCenter, __RC__)

! We accumulate fields over 24 hrs in the following:
! (stored in restart files to pass regression tests)

       call MAPL_AddInternalSpec(GC,                        &
          SHORT_NAME = 'GOCART::COqMassAcc',                &
          LONG_NAME  = 'CO Accumulation',                   &
          UNITS      = '1',                                 &
          FRIENDLYTO = ' ',                                 &
          RESTART    = MAPL_RestartOptional,                &
          DEFAULT    = 0.0,                                 &
          DIMS       = MAPL_DimsHorzVert,                   &
          VLOCATION  = MAPL_VLocationCenter, __RC__)

       call MAPL_AddInternalSpec(GC,                        &
          SHORT_NAME = 'GOCART::CH4qMassAcc',               &
          LONG_NAME  = 'CH4 Accumulation',                  &
          UNITS      = '1',                                 &
          FRIENDLYTO = ' ',                                 &
          RESTART    = MAPL_RestartOptional,                &
          DEFAULT    = 0.0,                                 &
          DIMS       = MAPL_DimsHorzVert,                   &
          VLOCATION  = MAPL_VLocationCenter, __RC__)

       call MAPL_AddInternalSpec(GC,                        &
          SHORT_NAME = 'GOCART::TempAcc',                   &
          LONG_NAME  = 'Temperature Accumulation',          &
          UNITS      = '1',                                 &
          FRIENDLYTO = ' ',                                 &
          RESTART    = MAPL_RestartOptional,                &
          DEFAULT    = 0.0,                                 &
          DIMS       = MAPL_DimsHorzVert,                   &
          VLOCATION  = MAPL_VLocationCenter, __RC__)

       call MAPL_AddInternalSpec(GC,                        &
          SHORT_NAME = 'GOCART::PressAcc',                  &
          LONG_NAME  = 'Pressure Accumulation',             &
          UNITS      = '1',                                 &
          FRIENDLYTO = ' ',                                 &
          RESTART    = MAPL_RestartOptional,                &
          DEFAULT    = 0.0,                                 &
          DIMS       = MAPL_DimsHorzVert,                   &
          VLOCATION  = MAPL_VLocationCenter, __RC__)

       call MAPL_AddInternalSpec(GC,                        &
          SHORT_NAME = 'GOCART::SpecHumAcc',                &
          LONG_NAME  = 'Specific Humidity Accumulation',    &
          UNITS      = '1',                                 &
          FRIENDLYTO = ' ',                                 &
          RESTART    = MAPL_RestartOptional,                &
          DEFAULT    = 0.0,                                 &
          DIMS       = MAPL_DimsHorzVert,                   &
          VLOCATION  = MAPL_VLocationCenter, __RC__)

       call MAPL_AddInternalSpec(GC,                        &
          SHORT_NAME = 'GOCART::DryAirMassAcc',             &
          LONG_NAME  = 'Dry Air Mass Accumulation',         &
          UNITS      = '1',                                 &
          FRIENDLYTO = ' ',                                 &
          RESTART    = MAPL_RestartOptional,                &
          DEFAULT    = 1.0,                                 &
          DIMS       = MAPL_DimsHorzVert,                   &
          VLOCATION  = MAPL_VLocationCenter, __RC__)


   DO i = 1, 17
   
     call MAPL_AddImportSpec(GC,                  &
          SHORT_NAME = TRIM(gmiSpeciesName3D(i)), &
          LONG_NAME  = TRIM(gmiSpeciesName3D(i)), &
          UNITS      = '1',                       &
          DIMS       = MAPL_DimsHorzVert,         &
          VLOCATION  = MAPL_VLocationCenter,      &
          __RC__)

   END DO

   DO i = 1, 1
   
     call MAPL_AddImportSpec(GC,                  &
          SHORT_NAME = TRIM(gmiSpeciesName2D(i)), &
          LONG_NAME  = TRIM(gmiSpeciesName2D(i)), &
          UNITS      = '1',                       &
          DIMS       = MAPL_DimsHorzOnly,         &
          VLOCATION  = MAPL_VLocationCenter,      &
          __RC__)

   END DO

   call MAPL_AddImportSpec(GC,                  &
        SHORT_NAME = 'ohRatio',                 &
        LONG_NAME  = 'OH_ratio',                &
        UNITS      = '1',                       &
        DIMS       = MAPL_DimsHorzVert,         &
        VLOCATION  = MAPL_VLocationCenter,      &
        __RC__)

   call MAPL_AddImportSpec(GC,                                         &
        SHORT_NAME = 'TAUCLI',                                         &
        LONG_NAME  = 'in_cloud_optical_thickness_for_ice_clouds',      &
        UNITS      = '1' ,                                             &
        DIMS       = MAPL_DimsHorzVert,                                &
        VLOCATION  = MAPL_VLocationCenter,                             &
        __RC__)

   call MAPL_AddImportSpec(GC,                                         &
        SHORT_NAME = 'TAUCLW',                                         &
        LONG_NAME  = 'in_cloud_optical_thickness_for_liquid_clouds',   &
        UNITS      = '1' ,                                             &
        DIMS       = MAPL_DimsHorzVert,                                &
        VLOCATION  = MAPL_VLocationCenter,                             &
        __RC__)


!!!
!!! EXPORTS
!!!

!   THIS COULD GO INTO   OH_ExportSpec___.h:

!    call MAPL_AddExportSpec(GC,  &
!       SHORT_NAME         = 'OHVMR',  &
!       LONG_NAME          = 'OH',  &
!       UNITS              = 'mol/mol', &
!       DIMS               = MAPL_DimsHorzVert,    &
!       VLOCATION          = MAPL_VLocationCenter,    &
!                                                      RC=STATUS  )
!    VERIFY_(STATUS)
     call MAPL_AddExportSpec(GC,  &
        SHORT_NAME         = 'O3upIn',  &
        LONG_NAME          = 'Overhead ozone',  &
        UNITS              = 'DU', &
        DIMS               = MAPL_DimsHorzVert,    &
        VLOCATION          = MAPL_VLocationCenter,    &
                                                       RC=STATUS  )
     VERIFY_(STATUS)


     call MAPL_AddExportSpec(GC,  &
        SHORT_NAME         = 'NOyIn',  &
        LONG_NAME          = 'Total Nitrogen Oxide',  &
        UNITS              = 'pptv', &
        DIMS               = MAPL_DimsHorzVert,    &
        VLOCATION          = MAPL_VLocationCenter,    &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,  &
        SHORT_NAME         = 'IsopreneIn',  &
        LONG_NAME          = 'Overhead ozone',  &
        UNITS              = 'pptv', &
        DIMS               = MAPL_DimsHorzVert,    &
        VLOCATION          = MAPL_VLocationCenter,    &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,  &
        SHORT_NAME         = 'ALK4In',  &
        LONG_NAME          = 'ALK4',  &
        UNITS              = 'pptv', &
        DIMS               = MAPL_DimsHorzVert,    &
        VLOCATION          = MAPL_VLocationCenter,    &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,  &
        SHORT_NAME         = 'AcetoneIn',  &
        LONG_NAME          = 'Acetone',  &
        UNITS              = 'pptv', &
        DIMS               = MAPL_DimsHorzVert,    &
        VLOCATION          = MAPL_VLocationCenter,    &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,  &
        SHORT_NAME         = 'PropaneIn',  &
        LONG_NAME          = 'Propane',  &
        UNITS              = 'pptv', &
        DIMS               = MAPL_DimsHorzVert,    &
        VLOCATION          = MAPL_VLocationCenter,    &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,  &
        SHORT_NAME         = 'PropeneIn',  &
        LONG_NAME          = 'Propene',  &
        UNITS              = 'pptv', &
        DIMS               = MAPL_DimsHorzVert,    &
        VLOCATION          = MAPL_VLocationCenter,    &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,  &
        SHORT_NAME         = 'EthaneIn',  &
        LONG_NAME          = 'Ethane',  &
        UNITS              = 'pptv', &
        DIMS               = MAPL_DimsHorzVert,    &
        VLOCATION          = MAPL_VLocationCenter,    &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,  &
        SHORT_NAME         = 'OzoneIn',  &
        LONG_NAME          = 'Ozone',  &
        UNITS              = 'ppbv', &
        DIMS               = MAPL_DimsHorzVert,    &
        VLOCATION          = MAPL_VLocationCenter,    &
                                                       RC=STATUS  )

! MANYIN - The name of this export MUST be changed:

     VERIFY_(STATUS)
     call MAPL_AddExportSpec(GC,  &
        SHORT_NAME         = 'varOut',  &
        LONG_NAME          = 'diagnostic variable',  &
        UNITS              = 'none', &
        DIMS               = MAPL_DimsHorzVert,    &
        VLOCATION          = MAPL_VLocationCenter,    &
                                                       RC=STATUS  )
     VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC,  &
        SHORT_NAME         = 'SpecHumIn',  &
        LONG_NAME          = 'Specific Humidity',  &
        UNITS              = '%', &
        DIMS               = MAPL_DimsHorzVert,    &
        VLOCATION          = MAPL_VLocationCenter,    &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,  &
        SHORT_NAME         = 'COIn',  &
        LONG_NAME          = 'CO',  &
        UNITS              = 'ppbv', &
        DIMS               = MAPL_DimsHorzVert,    &
        VLOCATION          = MAPL_VLocationCenter,    &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,  &
        SHORT_NAME         = 'CH4In',  &
        LONG_NAME          = 'CH4',  &
        UNITS              = 'ppbv', &
        DIMS               = MAPL_DimsHorzVert,    &
        VLOCATION          = MAPL_VLocationCenter,    &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,  &
        SHORT_NAME         = 'PressIn',  &
        LONG_NAME          = 'Pressure',  &
        UNITS              = 'hPa', &
        DIMS               = MAPL_DimsHorzVert,    &
        VLOCATION          = MAPL_VLocationCenter,    &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,  &
        SHORT_NAME         = 'TempIn',  &
        LONG_NAME          = 'Temperature',  &
        UNITS              = 'K', &
        DIMS               = MAPL_DimsHorzVert,    &
        VLOCATION          = MAPL_VLocationCenter,    &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,  &
        SHORT_NAME         = 'TRPPSIn',  &
        LONG_NAME          = 'Tropopause Level',  &
        UNITS              = 'hPa', &
        DIMS               = MAPL_DimsHorzOnly,    &
        VLOCATION          = MAPL_VLocationCenter,    &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,  &
        SHORT_NAME         = 'AlbedoIn',  &
        LONG_NAME          = 'Albedo',  &
        UNITS              = 'Fraction', &
        DIMS               = MAPL_DimsHorzOnly,    &
        VLOCATION          = MAPL_VLocationCenter,    &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,  &
        SHORT_NAME         = 'RavgaIn',  &
        LONG_NAME          = 'Reflectance Above a box',  &
        UNITS              = 'Fraction', &
        DIMS               = MAPL_DimsHorzVert,    &
        VLOCATION          = MAPL_VLocationCenter,    &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,  &
        SHORT_NAME         = 'RavgbIn',  &
        LONG_NAME          = 'Reflectance Below a box',  &
        UNITS              = 'Fraction', &
        DIMS               = MAPL_DimsHorzVert,    &
        VLOCATION          = MAPL_VLocationCenter,    &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

! MANYIN - These should be set up in
!#     include "OH_ExportSpec___.h"
! and that can be included here or in Aero_GridCompMod.F90

   RETURN_(ESMF_SUCCESS)

   end subroutine OH_GridCompSetServices

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  OH_GridCompInitialize --- Initialize OH_GridComp
!
! !INTERFACE:
!

      subroutine OH_GridCompInitialize ( this, w_c, gc, impChem, expChem, &
                                      nymd, nhms, cdt, rc )
! !USES:

     IMPLICIT NONE

! !INPUT PARAMETERS:
      INTEGER, INTENT(IN) :: nymd, nhms        ! time
      REAL,    INTENT(IN) :: cdt               ! chemical timestep (secs)
      TYPE(Chem_Bundle), intent(in) :: w_c     ! Chemical tracer fields      
!
! !OUTPUT PARAMETERS:

      TYPE(OH_GridComp),   INTENT(INOUT) :: this     ! Grid Component
      TYPE(ESMF_GridComp), INTENT(INOUT) :: gc       ! use the ESMF state handle to access internal fields
      TYPE(ESMF_State),    INTENT(INOUT) :: impChem  ! Import State
      TYPE(ESMF_State),    INTENT(INOUT) :: expChem  ! Export State
      INTEGER, INTENT(OUT) ::  rc                  ! Error return code:
                                                   !  0 - all is well
                                                   !  1 - 
!
! !DESCRIPTION: 
! Initializes the OH Parameterization Grid Component. 
!
! !REVISION HISTORY:
!
!  24June2010 Nielsen  First crack.
!
!EOP
!-------------------------------------------------------------------------

      CHARACTER(LEN=*), PARAMETER :: myname = 'OH_GridCompInitialize'
      CHARACTER(LEN=*), PARAMETER :: Iam    = myname
      CHARACTER(LEN=255) :: name
      integer            :: status
   
      INTEGER :: ios, j, i, ic, n
      INTEGER, ALLOCATABLE :: ier(:)

      INTEGER :: iDT, steps, time_of_day
      INTEGER, PARAMETER :: SecPerDay = 24*60*60

      INTEGER :: JDAY, curMonth
      INTEGER :: i1, i2, im, j1, j2, jm, km
      INTEGER :: nbeg, nend, begTime, nTimes, incSecs

      REAL, POINTER, DIMENSION(:,:,:) ::  ptr3d      => null()

      type(MAPL_MetaComp), pointer :: genState    ! MAPL generic state
      type(ESMF_State)             :: INTERNAL
      type(ESMF_Field)             :: field

      this%rcfilen = 'OH_GridComp.rc'
      this%name = 'GEOS-5/GOCART OH Parameterization Package'

      IF(MAPL_AM_I_ROOT()) PRINT *,myname,": Begin..."

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

      !------------------------------------------------
      ! i_OH+0: the only OH bin
      ! i_OH+1: 24 hour average of CO
      ! i_OH+2: 24 hour average of CH4
      ! i_OH+3: 24 hour average of temperature
      ! i_OH+4: 24 hour average of atmospheric pressure
      ! i_OH+5: 24 hour average of specific humidity
      ! i_OH+6: 24 hour average of dry air mass     -> was this ONLY USED IN AN OLD VERSION OF GETINFO ?  MANYIN
      ! j_OH = i_OH + 6
      !------------------------------------------------

      nbeg  = w_c%reg%i_OH
      nend  = w_c%reg%j_OH

!  It requires exactly 1 instance
!  -----------------
!     if ( (nend - nbeg) /= 6 ) then     MANYIN
      if ( nend /= nbeg ) then
         IF(MAPL_AM_I_ROOT()) PRINT *,myname,": Must use exactly 1 OH instance"
         VERIFY_(10)
      end if

    CALL MAPL_GetObjectFromGC( gc, genState, __RC__ )

    CALL MAPL_Get( genState, INTERNAL_ESMF_STATE=INTERNAL, __RC__ )

    call MAPL_GetPointer(INTERNAL, ptr3d, 'GOCART::CO_24h'      ,  __RC__ )
    call MAPL_GetPointer(INTERNAL, ptr3d, 'GOCART::CH4_24h'     ,  __RC__ )
    call MAPL_GetPointer(INTERNAL, ptr3d, 'GOCART::Temp_24h'    ,  __RC__ )
    call MAPL_GetPointer(INTERNAL, ptr3d, 'GOCART::Press_24h'   ,  __RC__ )
    call MAPL_GetPointer(INTERNAL, ptr3d, 'GOCART::SpecHum_24h' ,  __RC__ )

!!  Alternate versions:

!!    call MAPL_GetPointer(INTERNAL, NAME='GOCART::CO_24h'      , ptr=ptr3d,  __RC__ )
!!    call MAPL_GetPointer(INTERNAL, NAME='GOCART::CH4_24h'     , ptr=ptr3d,  __RC__ )
!!    call MAPL_GetPointer(INTERNAL, NAME='GOCART::Temp_24h'    , ptr=ptr3d,  __RC__ )
!!    call MAPL_GetPointer(INTERNAL, NAME='GOCART::Press_24h'   , ptr=ptr3d,  __RC__ )
!!    call MAPL_GetPointer(INTERNAL, NAME='GOCART::SpecHum_24h' , ptr=ptr3d,  __RC__ )
!!
!!! MANYIN:  probably do not need these calls;
!!
!!    call ESMF_StateGet(internal, 'GOCART::CO_24h', field, __RC__)
!!    call ESMF_AttributeSet(field, NAME='ScavengingFractionPerKm', VALUE=0.0, __RC__)
!!    call ESMF_StateGet(internal, 'GOCART::CH4_24h', field, __RC__)
!!    call ESMF_AttributeSet(field, NAME='ScavengingFractionPerKm', VALUE=0.0, __RC__)
!!    call ESMF_StateGet(internal, 'GOCART::Temp_24h', field, __RC__)
!!    call ESMF_AttributeSet(field, NAME='ScavengingFractionPerKm', VALUE=0.0, __RC__)
!!    call ESMF_StateGet(internal, 'GOCART::Press_24h', field, __RC__)
!!    call ESMF_AttributeSet(field, NAME='ScavengingFractionPerKm', VALUE=0.0, __RC__)
!!    call ESMF_StateGet(internal, 'GOCART::SpecHum_24h', field, __RC__)
!!    call ESMF_AttributeSet(field, NAME='ScavengingFractionPerKm', VALUE=0.0, __RC__)


      allocate(ier(199))
      ier(:) = 0

!  Compute min_time
!  ----------------
!     #steps to test HMS values over the course of 48 hours  (much more than needed!)
      iDT = INT(cdt)
      steps = (48*60*60) / iDT

      time_of_day = convertTimeToSeconds(nhms)

      this%min_time = time_of_day
      DO i=1,steps
        time_of_day = MOD(time_of_day+iDT,SecPerDay)
        this%min_time = MIN(this%min_time,time_of_day)
      END DO
      

!  Load resource file
!  ------------------
      CALL I90_loadf ( TRIM(this%rcfilen), ier(1) )
      if ( ier(1) .NE. 0 ) then
         VERIFY_(11)
      end if

!  Parse resource file
!  -------------------
! now paths are listed in ExtData
!     CALL I90_label ( 'GMI_infile_name:', ier(3) )
!     CALL I90_gtoken( this%GMI_infile_name, ier(4) )

! only needed if we are calling adjustH2OMER
!     CALL I90_label ( 'Qscale_infile_name:', ier(5) )
!     CALL I90_gtoken( this%Qscl_infile_name, ier(6) )

      IF ( ANY(ier(:) /= 0)) THEN
!        PRINT*, TRIM(myname),': Failed to parse GMI or Qscale infile name.'
         PRINT*, TRIM(myname),': Failed to parse GMI infile name.'
         VERIFY_(14)
      END IF
      ier(:)=0

      allocate(        this%CO_24h(i1:i2,j1:j2,1:km), STAT=ier(10))
      allocate(       this%CH4_24h(i1:i2,j1:j2,1:km), STAT=ier(11))
      allocate(      this%Temp_24h(i1:i2,j1:j2,1:km), STAT=ier(12))
      allocate(     this%Press_24h(i1:i2,j1:j2,1:km), STAT=ier(13))
      allocate(   this%SpecHum_24h(i1:i2,j1:j2,1:km), STAT=ier(14))
      allocate(  this%rxnRateRatio(i1:i2,j1:j2,1:km), STAT=ier(16))


      IF ( ANY(ier(:) /= 0)) THEN
         PRINT*, TRIM(myname),': Failed to allocate work space.'
         VERIFY_(12)
      END IF
      ier(:) = 0

!
!     Restore from INTERNAL state
!     Convert to REAL*8, bottom up
!

      call MAPL_GetPointer(INTERNAL, ptr3d, 'GOCART::CO_24h', __RC__ )
      this%CO_24h(:,:,1:km) = ptr3d(:,:,km:1:-1)

      call MAPL_GetPointer(INTERNAL, ptr3d, 'GOCART::CH4_24h', __RC__ )
      this%CH4_24h(:,:,1:km) = ptr3d(:,:,km:1:-1)

      call MAPL_GetPointer(INTERNAL, ptr3d, 'GOCART::Temp_24h', __RC__ )
      this%Temp_24h(:,:,1:km) = ptr3d(:,:,km:1:-1)

      call MAPL_GetPointer(INTERNAL, ptr3d, 'GOCART::Press_24h', __RC__ )
      this%Press_24h(:,:,1:km) = ptr3d(:,:,km:1:-1)

      call MAPL_GetPointer(INTERNAL, ptr3d, 'GOCART::SpecHum_24h', __RC__ )
      this%SpecHum_24h(:,:,1:km) = ptr3d(:,:,km:1:-1)

      this%rxnRateRatio   =       1.0d0


      allocate( this%Ravga(i1:i2, j1:j2, km), this%Ravgb(i1:i2, j1:j2, km), STAT=ier(26))
      allocate(this%LPAUSE(i1:i2, j1:j2), this%troppresIn(i1:i2, j1:j2), this%qIn(i1:i2,j1:j2,km), STAT=ier(27))

      IF ( ANY(ier(:) /= 0)) THEN
         PRINT*, TRIM(myname),': Failed to allocate space.'
         VERIFY_(15)
      END IF
      ier(:) = 0
   
      this%Ravga = 0.0d0
      this%Ravgb = 0.0d0

          this%LPAUSE(i1:i2, j1:j2)     = 0
      this%troppresIn(i1:i2, j1:j2)     = 0.0d0
             this%qIn(i1:i2, j1:j2, : ) = 0.0d0

      CALL initializeOHparam(this%OHparam, i2, j2, km, this%rcfilen)

      DEALLOCATE(ier)

      RETURN

 end subroutine OH_GridCompInitialize

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  OH_GridCompRun --- Run OH_GridComp
!
! !INTERFACE:
!

   subroutine OH_GridCompRun ( this, w_c, gc, impChem, expChem, &
                                      nymd, nhms, cdt, rc )

! !USES:

  IMPLICIT NONE

! !INPUT PARAMETERS:

   TYPE(Chem_Bundle), intent(in) :: w_c        ! Chemical tracer fields      
   INTEGER, INTENT(IN) :: nymd, nhms           ! time
   REAL,    INTENT(IN) :: cdt                  ! chemical timestep (secs)


! !OUTPUT PARAMETERS:

   TYPE(OH_GridComp),   INTENT(INOUT) :: this     ! Grid Component
   TYPE(ESMF_GridComp), INTENT(INOUT) :: gc       ! use the ESMF state handle to access internal fields
   TYPE(ESMF_State),    INTENT(INOUT) :: impChem  ! Import State
   TYPE(ESMF_State),    INTENT(INOUT) :: expChem  ! Export State
   INTEGER, INTENT(OUT) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 - 

! !DESCRIPTION: Runs the OH Parameterization Grid Component. Multiple instance
!               version.
!
! !REVISION HISTORY:
!
!  24June2010 Nielsen  First crack.
!
!EOP
!-------------------------------------------------------------------------

      INTEGER, ALLOCATABLE :: ier(:)
      CHARACTER(LEN=*), PARAMETER :: myname = 'OH_GridCompRun'
      CHARACTER(LEN=*), PARAMETER :: Iam = myname

!  Input fields from fvGCM
!  -----------------------

      REAL, POINTER, DIMENSION(:,:)   ::  tropp     => null()
!     REAL, POINTER, DIMENSION(:,:)   ::  albedo    => null()
      REAL, POINTER, DIMENSION(:,:)   ::  albedoVF  => null()
!     REAL, POINTER, DIMENSION(:,:)   ::  albedoVR  => null()
!     REAL, POINTER, DIMENSION(:,:)   ::  albedoNF  => null()
!     REAL, POINTER, DIMENSION(:,:)   ::  albedoNR  => null()      
      REAL, POINTER, DIMENSION(:,:)   ::  cellArea  => null()

!  All are  top-down:
      REAL, POINTER, DIMENSION(:,:,:) ::  T          => null()
      REAL, POINTER, DIMENSION(:,:,:) ::  TV         => null()
      REAL, POINTER, DIMENSION(:,:,:) ::  Q          => null()
      REAL, POINTER, DIMENSION(:,:,:) ::  exO3up     => null()
!     REAL, POINTER, DIMENSION(:,:,:) ::  exOHVMR    => null()
      REAL, POINTER, DIMENSION(:,:,:) ::  exNOy      => null()
      REAL, POINTER, DIMENSION(:,:,:) ::  exvarOut   => null()
      REAL, POINTER, DIMENSION(:,:,:) ::  exALK4     => null()
      REAL, POINTER, DIMENSION(:,:,:) ::  exIsoprene => null()
      REAL, POINTER, DIMENSION(:,:,:) ::  exAcetone  => null()
      REAL, POINTER, DIMENSION(:,:,:) ::  exPropane  => null()
      REAL, POINTER, DIMENSION(:,:,:) ::  exPropene  => null()
      REAL, POINTER, DIMENSION(:,:,:) ::  exEthane   => null()
      REAL, POINTER, DIMENSION(:,:,:) ::  exOzone    => null()
      REAL, POINTER, DIMENSION(:,:,:) ::  exSpecHum  => null()
      REAL, POINTER, DIMENSION(:,:,:) ::  exCO       => null()
      REAL, POINTER, DIMENSION(:,:,:) ::  exCH4      => null()
      REAL, POINTER, DIMENSION(:,:,:) ::  exPress    => null()
      REAL, POINTER, DIMENSION(:,:,:) ::  exTemp     => null()
      REAL, POINTER, DIMENSION(:,:)   ::  exTRPPS    => null()
      REAL, POINTER, DIMENSION(:,:)   ::  exAlbedo   => null()
      REAL, POINTER, DIMENSION(:,:,:) ::  exRavga    => null()
      REAL, POINTER, DIMENSION(:,:,:) ::  exRavgb    => null()
      REAL, POINTER, DIMENSION(:,:,:) ::  exRLAT     => null()
      REAL, POINTER, DIMENSION(:,:,:) ::  exDEC      => null()
      REAL, POINTER, DIMENSION(:,:,:) ::  ple        => null()
      REAL, POINTER, DIMENSION(:,:,:) ::  zle        => null()
      REAL, POINTER, DIMENSION(:,:,:) ::  rhoWet     => null()
      REAL, POINTER, DIMENSION(:,:,:) ::  rhoDry     => null()
      REAL, POINTER, DIMENSION(:,:,:) ::  taucli     => null()
      REAL, POINTER, DIMENSION(:,:,:) ::  tauclw     => null()
      REAL, POINTER, DIMENSION(:,:,:) ::  fcld       => null()
      REAL, POINTER, DIMENSION(:,:)   ::  ptr2d      => null()
      REAL, POINTER, DIMENSION(:,:,:) ::  ptr3d      => null()

! Pointers to INTERNAL fields  (all top-down)
!     24 hour averages
      REAL, POINTER, DIMENSION(:,:,:) ::  pInternal_CO_24h        => null()
      REAL, POINTER, DIMENSION(:,:,:) ::  pInternal_CH4_24h       => null()
      REAL, POINTER, DIMENSION(:,:,:) ::  pInternal_Temp_24h      => null()
      REAL, POINTER, DIMENSION(:,:,:) ::  pInternal_Press_24h     => null()
      REAL, POINTER, DIMENSION(:,:,:) ::  pInternal_SpecHum_24h   => null()
!     Accumulation fields
      REAL, POINTER, DIMENSION(:,:,:) ::  pInternal_COqMassAcc    => null()
      REAL, POINTER, DIMENSION(:,:,:) ::  pInternal_CH4qMassAcc   => null()
      REAL, POINTER, DIMENSION(:,:,:) ::  pInternal_TempAcc       => null()
      REAL, POINTER, DIMENSION(:,:,:) ::  pInternal_PressAcc      => null()
      REAL, POINTER, DIMENSION(:,:,:) ::  pInternal_SpecHumAcc    => null()
      REAL, POINTER, DIMENSION(:,:,:) ::  pInternal_DryAirMassAcc => null()


      REAL, ALLOCATABLE :: pe(:,:,:),p(:,:,:),ndDry(:,:,:),ndWet(:,:,:)   !   top-down
      REAL, ALLOCATABLE :: mw_moist_air(:,:,:)                  !   bottom-up

      REAL*8,  ALLOCATABLE :: solarZenithAngle    (:,:)
      REAL*8,  ALLOCATABLE :: o3up       (:,:,:)                !   bottom-up
      REAL*8,  ALLOCATABLE :: STT        (:,:,:,:)    !  ppbv       bottom-up
      REAL*8,  ALLOCATABLE :: Qppmv      (:,:,:)      !             bottom-up
      REAL*8,  ALLOCATABLE :: OHparamNum (:,:,:)      !             bottom-up
      REAL*8,  ALLOCATABLE :: CH4VAR     (:,:,:)      !             bottom-up
      REAL*8,  ALLOCATABLE :: RLON       (:,:)   ! Array of longitudes
      REAL*8,  ALLOCATABLE :: RLAT       (:,:)   ! Array of latitudes
      REAL*8,  ALLOCATABLE :: BBIJ       (:, :, :,:)  !             bottom-up
      REAL*8,  ALLOCATABLE :: gmiOH(:,:,:) !-sas                    bottom-up
      REAL*8,  ALLOCATABLE :: mass(:,:,:)             !  of dry air           top-down
      REAL*8,  ALLOCATABLE :: gridBoxThickness(:,:,:) !                       top-down
      REAL*8,  ALLOCATABLE :: PL(:,:,:)               !  never used          (top-down)
      REAL*8,  ALLOCATABLE :: optDepth(:,:,:)         !             bottom-up
      REAL*8,  ALLOCATABLE :: ALBD(:,:)
      REAL*8,  ALLOCATABLE :: tropPress(:,:)
      REAL,    ALLOCATABLE :: var3d(:,:,:)

      INTEGER :: i1, i2, im, j1, j2, jm, km, ios, iXj
      INTEGER :: i, j, k, kReverse, n, nbeg, nend, ic

      INTEGER :: JDAY, curMonth
      INTEGER :: ix, jx
      REAL*8  :: time_sec, offset_sec
      integer :: STATUS
      REAL    :: radToDeg, degToRad, pi
      REAL    :: qmin, qmax, h20factor

      REAL, PARAMETER ::     epsilon = (MAPL_H2OMW/MAPL_AIRMW)

      type(MAPL_MetaComp), pointer :: genState    ! MAPL generic state
      type(ESMF_State)             :: INTERNAL

      REAL :: massFactor = 1.0E15   ! divide MASS by this, to keep the accumulation from approaching UNDEF
      REAL :: avgFactor
      INTEGER, PARAMETER :: SecPerDay = 24*60*60


      avgFactor = cdt/secPerDay
     
      pi       = 4.00*ATAN(1.00)
      degToRad = pi/180.00
      radToDeg = 180.00/pi

      radToDeg = 57.2957795

      IF(MAPL_AM_I_ROOT()) PRINT *,myname,": begin...."

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

      allocate(ier(128))
      ier(:) = 0

      nbeg  = w_c%reg%i_OH
      nend  = w_c%reg%j_OH

      !  It requires a single bin
      !  ------------------------
      if ( nend /= nbeg ) then
         IF(MAPL_AM_I_ROOT()) PRINT *,myname,": Must have 1 OH bin"
         VERIFY_(20)
      end if

      !Get imports

      CALL MAPL_GetPointer(impChem,    rhoWet,      'AIRDENS',  rc=ier(40))   ! top-down
      CALL MAPL_GetPointer(impChem,    rhoDry, 'AIRDENS_DRYP',  rc=ier(41))   ! top-down
      CALL MAPL_GetPointer(impChem,         Q,            'Q',  rc=ier(42))   ! top-down
      CALL MAPL_GetPointer(impChem,         T,            'T',  rc=ier(43))   ! top-down
      CALL MAPL_GetPointer(impChem,        TV,           'TV',  rc=ier(44))   ! top-down
      CALL MAPL_GetPointer(impChem,       zle,          'ZLE',  rc=ier(45))   ! top-down
      CALL MAPL_GetPointer(impChem,       ple,          'PLE',  rc=ier(46))   ! top-down
      CALL MAPL_GetPointer(impChem,     tropp,        'TROPP',  rc=ier(47))
      CALL MAPL_GetPointer(impChem,    taucli,       'TAUCLI',  rc=ier(48))   ! top-down
      CALL MAPL_GetPointer(impChem,    tauclw,       'TAUCLW',  rc=ier(49))   ! top-down
      CALL MAPL_GetPointer(impChem,  cellArea,         'AREA',  rc=ier(50)) 
      CALL MAPL_GetPointer(impChem,      fcld,         'FCLD',  rc=ier(51))   ! top-down
      CALL MAPL_GetPointer(impChem,  albedoVF,        'ALBVF',  rc=ier(52))
!     CALL MAPL_GetPointer(impChem,  albedoVR,        'ALBVR',  rc=ier(53))
!     CALL MAPL_GetPointer(impChem,  albedoNF,        'ALBNF',  rc=ier(54))
!     CALL MAPL_GetPointer(impChem,  albedoNR,        'ALBNR',  rc=ier(55)) 

      if ( any(ier(:) /= 0) ) then
!       IF (MAPL_AM_I_ROOT()) THEN
          DO i=40,55
            print *, 'GetPtr error code ',i,' = ', ier(i)
          END DO
!       END IF
        PRINT *,TRIM(myname)//':'//TRIM(Iam),": Failed to obtain the import state"
        VERIFY_(10)
      end if
      ier(:) = 0

      CALL MAPL_GetObjectFromGC( gc, genState, __RC__ )
      CALL MAPL_Get( genState, INTERNAL_ESMF_STATE=INTERNAL, __RC__ )

      CALL MAPL_GetPointer(INTERNAL, pInternal_CO_24h,        'GOCART::CO_24h',        __RC__ )
      CALL MAPL_GetPointer(INTERNAL, pInternal_CH4_24h,       'GOCART::CH4_24h',       __RC__ )
      CALL MAPL_GetPointer(INTERNAL, pInternal_Temp_24h,      'GOCART::Temp_24h',      __RC__ )
      CALL MAPL_GetPointer(INTERNAL, pInternal_Press_24h,     'GOCART::Press_24h',     __RC__ )
      CALL MAPL_GetPointer(INTERNAL, pInternal_SpecHum_24h,   'GOCART::SpecHum_24h',   __RC__ )

      CALL MAPL_GetPointer(INTERNAL, pInternal_COqMassAcc,    'GOCART::COqMassAcc',    __RC__ )
      CALL MAPL_GetPointer(INTERNAL, pInternal_CH4qMassAcc,   'GOCART::CH4qMassAcc',   __RC__ )
      CALL MAPL_GetPointer(INTERNAL, pInternal_TempAcc,       'GOCART::TempAcc',       __RC__ )
      CALL MAPL_GetPointer(INTERNAL, pInternal_PressAcc,      'GOCART::PressAcc',      __RC__ )
      CALL MAPL_GetPointer(INTERNAL, pInternal_SpecHumAcc,    'GOCART::SpecHumAcc',    __RC__ )
      CALL MAPL_GetPointer(INTERNAL, pInternal_DryAirMassAcc, 'GOCART::DryAirMassAcc', __RC__ )

      ! i_OH+1: 24 hour average of CO
      ! i_OH+2: 24 hour average of CH4
      ! i_OH+3: 24 hour average of temperature
      ! i_OH+4: 24 hour average of atmospheric pressure
      ! i_OH+5: 24 hour average of specific humidity
      ! i_OH+6: 24 hour average of dry air mass     -> was this ONLY USED IN AN OLD VERSION OF GETINFO ?  MANYIN
      !                                                So for now, we do not store pInternal_DryAirMass_24h



      ! Reserve Some local work space
      !------------------------------
      allocate(            RLON(i1:i2,j1:j2),STAT=ier(11))
      allocate(            RLAT(i1:i2,j1:j2),STAT=ier(12))
      allocate(            ALBD(i1:i2,j1:j2),STAT=ier(13))
      allocate(       tropPress(i1:i2,j1:j2),STAT=ier(14))
      allocate(solarZenithAngle(i1:i2,j1:j2),STAT=ier(15))

      allocate(            o3up(i1:i2,j1:j2,1:km),STAT=ier(16))
      allocate(            mass(i1:i2,j1:j2,1:km),STAT=ier(17))
      allocate(        optDepth(i1:i2,j1:j2,1:km),STAT=ier(18))
      allocate(gridBoxThickness(i1:i2,j1:j2,1:km),STAT=ier(19))
      allocate(           var3d(i1:i2,j1:j2,1:km),STAT=ier(20))
      allocate(              PL(i1:i2,j1:j2,1:km         ),STAT=ier(22))
      allocate(              pe(i1:i2,j1:j2,1:km+1       ),STAT=ier(23))
      allocate(               p(i1:i2,j1:j2,1:km         ),STAT=ier(24))
      allocate(             STT(i1:i2,j1:j2,1:km,NNPAR   ),STAT=ier(25))
      allocate(            BBIJ(i1:i2,j1:j2,1:km,NFIELDS2),STAT=ier(26))
      allocate(          CH4VAR(i1:i2,j1:j2,1:km         ),STAT=ier(27))
      allocate(           Qppmv(i1:i2,j1:j2,1:km         ),STAT=ier(28))
      allocate(           ndDry(i1:i2,j1:j2,1:km         ),STAT=ier(29))
      allocate(           ndWet(i1:i2,j1:j2,1:km         ),STAT=ier(30))
      allocate(    mw_moist_air(i1:i2,j1:j2,1:km         ),STAT=ier(31))
      allocate(      OHparamNum(i1:i2,j1:j2,1:km         ),STAT=ier(32))
      allocate(           gmiOH(i1:i2,j1:j2,1:km         ),STAT=ier(33)) !sas

      IF ( ANY(ier(:) /= 0)) THEN
         PRINT*, TRIM(myname),': Failed to allocate work space.'
         VERIFY_(21)
      END IF
      ier(:) = 0

      ! Extract latitude and longitude - convert to degrees ! EY
      RLON(:,:) = w_c%grid%lon(:,:)*radToDeg 
      RLAT(:,:) = w_c%grid%lat(:,:)*radToDeg 

   !  Layer interface pressures (top-down)
   !  ------------------------------------
      pe(i1:i2,j1:j2,1)=w_c%grid%ptop

      DO k=2,km+1
         pe(i1:i2,j1:j2,k)=pe(i1:i2,j1:j2,k-1)+w_c%delp(i1:i2,j1:j2,k-1)
      END DO

   !  Layer mean pressures  (Pa)
   !  --------------------
      DO k=1,km
         p(i1:i2,j1:j2,k)=(pe(i1:i2,j1:j2,k)+pe(i1:i2,j1:j2,k+1))*0.50
      END DO

   !    Dry air number density  (molec/m3)
   !    (Note that this is not consistent with fields like T and P!)
   !  ------------------------
      ndDry(:,:,1:km) = rhoDry(:,:,1:km)*MAPL_AVOGAD/MAPL_AIRMW

   !  Moist air number density  (molec/m3)
   !  ------------------------
      ndWet = (MAPL_AVOGAD * p) / (MAPL_RUNIV * TV)

   !  Molecular Weight of Moist Air  (kg/mol)
   !  MANYIN - This works well in theory, but less so in practice;
   !           e.g. the value differs a bit from 28.97 in the stratosphere
   !           for now just use Dry Air MW
   !  -----------------------------------------------------------------------
!     mw_moist_air(:,:,1:km) = rhoWet(:,:,km:1:-1)*MAPL_AVOGAD/ndWet(:,:,km:1:-1)
      mw_moist_air(:,:,1:km) = MAPL_AIRMW

 
   !  Units for OH: w_c%qa(nbeg)%data3D(i1:i2,j1:j2,1:km) is always in terms of molec/cm3


   !  Cell depth (in meters)
   !  ---------------------
   !  (Another approach:)
   !  gridBoxThickness(:,:,:) = w_c%delp(:,:,:)/(rhoWet(:,:,:)*MAPL_GRAV)
      gridBoxThickness(:,:,1:km) = zle(:,:,0:km-1)-zle(:,:,1:km)

   ! Read monthly mean tropPress from GMI file
      CALL MAPL_GetPointer(impChem, ptr2d, gmiSpeciesName2D(IC_TRPPS), __RC__ )
      this%troppresIn(:,:) = ptr2d(:,:)

      ! tropPress(i1:i2,j1:j2)       =  tropp(i1:i2,j1:j2)*Pa2hPa

   ! Also read in GMI water vapor	
      ! Read in 3D field - monthly mean gmi 

      !ic = 17
      !CALL MAPL_GetPointer(impChem, ptr3d, gmiSpeciesName3D(ic), __RC__ )
      ! This proc reversed vertical levels into GEOS5 grid; reverse back to GMI levels	
      !this%qIn(:,:,1:km) = ptr3d(:,:,km:1:-1)


      ALBD(i1:i2,j1:j2) = albedoVF(i1:i2,j1:j2) !+albedoNF(i1:i2,j1:j2))/2.0

      ! get OH ratios for aerosol correction

          ! This proc retrieves arrays on GEOS-5 levels
          CALL MAPL_GetPointer(impChem, ptr3d, 'ohRatio', __RC__ )

          ! Reverse levels back to GMI (& OH param) convention 1(at surface) , 72 (at top of atmos) 
          this%OHparam%dustrat(:,:,1:km) = ptr3d(:,:,km:1:-1)


      IF (this%min_time == nhms ) THEN

          ! Store the accumulated averages into the INTERNAL fields
               pInternal_Temp_24h(:,:,:) =       pInternal_TempAcc(:,:,:)
              pInternal_Press_24h(:,:,:) =      pInternal_PressAcc(:,:,:)
            pInternal_SpecHum_24h(:,:,:) =    pInternal_SpecHumAcc(:,:,:)

          ! And install them in the REAL*8 working copies  (reverse in vertical)
                this%Temp_24h(:,:,1:km) =    pInternal_Temp_24h(:,:,km:1:-1)
               this%Press_24h(:,:,1:km) =   pInternal_Press_24h(:,:,km:1:-1)
             this%specHum_24h(:,:,1:km) = pInternal_SpecHum_24h(:,:,km:1:-1)

          ! Reset the accumulators
               pInternal_TempAcc(:,:,:) = 0.0
              pInternal_PressAcc(:,:,:) = 0.0
            pInternal_SpecHumAcc(:,:,:) = 0.0

        ! OHparamGrid Comp can run on its own, or together with CO and/or CH4
        ! If running CO Grid Comp, monitor 24 hr average CO 
        IF ( w_c%reg%doing_CO ) THEN
          ! Desired units for CO_24h - parts per part wrt Dry Air
          ! Store INTERNAL field of CO avg, in GEOS-5 vertical coords
          pInternal_CO_24h(:,:,:) = pInternal_COqMassAcc(:,:,:)/pInternal_DryAirMassAcc(:,:,:)

          ! And install it in our REAL*8 working copy  (reverse in vertical)
          this%CO_24h(:,:,1:km) = pInternal_CO_24h(:,:,km:1:-1)

          ! Reset the accumulator
          pInternal_COqMassAcc(:,:,:) = 0.0
        END IF

        ! If running CH4 Grid Comp, monitor 24 hr average CH4 
        IF ( w_c%reg%doing_CH4) THEN
          ! Desired units for CH4_24h - parts per part wrt Dry Air
          ! Store INTERNAL field of CH4 avg, in GEOS-5 vertical coords
          pInternal_CH4_24h(:,:,:) = pInternal_CH4qMassAcc(:,:,:)/pInternal_DryAirMassAcc(:,:,:)

          ! And install it in our REAL*8 working copy  (reverse in vertical)
          this%CH4_24h(:,:,1:km) = pInternal_CH4_24h(:,:,km:1:-1)

          ! Reset the accumulator
          pInternal_CH4qMassAcc(:,:,:) = 0.0
        END IF

        pInternal_DryAirMassAcc(:,:,:) = 0.0

      END IF

      ! Update the Accumulators

         pInternal_TempAcc(:,:,:) =    pInternal_TempAcc(:,:,:) + avgFactor * T(:,:,:)

        pInternal_PressAcc(:,:,:) =   pInternal_PressAcc(:,:,:) + avgFactor * 0.5*(ple(:,:,1:km)+ &
                                                                                   ple(:,:,0:km-1))*Pa2hPa

      pInternal_SpecHumAcc(:,:,:) = pInternal_SpecHumAcc(:,:,:) + avgFactor * Q(:,:,:)

      !--------------------------------------------------------------				     
      ! Sum up mass for the current day; divide by a factor to keep the number smaller; then remultiply at the end 
      ! These quantities are on GEOS-5 vertical grid
      !--------------------------------------------------------------				     
      DO k=1,km
         mass(:,:,k)=rhoDry(:,:,k)*cellArea(:,:)*(zle(:,:,k-1)-zle(:,:,k))  ! kg
      END DO

      pInternal_DryAirMassAcc(:,:,:) = &
      pInternal_DryAirMassAcc(:,:,:) + mass(:,:,:)/massFactor
    
      IF ( w_c%reg%doing_CO) THEN !sas -if not doing_CO, skip next line
         ! Manyin - convert CO units:  VMR wrt moist  to  VMR wrt dry, assuming MW_moistair == MW_dryair
         ! Manyin - Also (apparently) we do a mass-weighted average
         pInternal_COqMassAcc(:,:,:) = &
         pInternal_COqMassAcc(:,:,:) + &
             ( w_c%qa(w_c%reg%i_CO)%data3D(:,:,:)  / (1. - Q(:,:,:)) ) *  mass(:,:,:)/massFactor
      ENDIF

      IF ( w_c%reg%doing_CH4) THEN !sas -if not doing_CH4, skip next line
         ! Manyin - convert CH4 units:  VMR wrt moist  to  VMR wrt dry, assuming MW_moistair == MW_dryair
         ! Manyin - Also (apparently) we do a mass-weighted average
         pInternal_CH4qMassAcc(:,:,:) = &
         pInternal_CH4qMassAcc(:,:,:) + &
             ( w_c%qa(w_c%reg%i_CH4)%data3D(:,:,:)  / (1. - Q(:,:,:)) ) *  mass(:,:,:)/massFactor
      ENDIF


      optDepth(:,:,1:km)   = (taucli(:,:,km:1:-1) + &
                              tauclw(:,:,km:1:-1))*fcld(:,:,km:1:-1)**1.5



      !-------------------------------------------------
      ! Determine the tropopause level for each grid box
      !-------------------------------------------------

      this%LPAUSE = 0   

!     DO k = km,2,-1
!        DO jx = j1, j2
!           DO ix = i1, i2
!              IF ((tropPress(ix,jx) .LE. PL(ix,jx,k)) .AND. &
!                  (tropPress(ix,jx) .GT. PL(ix,jx,k-1))) THEN
!                 this%LPAUSE(ix,jx) = km - k + 1
!              END IF
!           END DO
!        END DO
!     END DO

      !CALL calcTroppLevel()
      ! Use a constant value of ~72 hPa for tropopause; the correct OH in the
      ! stratosphere within the CH4 and CO code
! MANYIN:  This should be changed from 37 to 36 (which is ~ 66.6 - 78.5 hPa)
      this%LPAUSE = 36

      call calcRxnRateFix()

      JDAY     = JulianDay(nymd)
      curMonth = MOD(nymd, 10000) / 100

      offset_sec = 0.0d0
      time_sec   = convertTimeToSeconds(nhms)

!     solarZenithAngle(i1:i2,j1:j2) = computeSolarZenithAngle (JDAY, time_sec, &
!                                     offset_sec, RLAT, RLON, i1, i2, j1, j2)
      solarZenithAngle(i1:i2,j1:j2) = acos(w_c%cosz(i1:i2,j1:j2))*radToDeg    

      CALL avgrefl(optDepth, solarZenithAngle, this%LPAUSE, this%Ravga, this%Ravgb, &
                            nhms, i2, j2, km)

      OHparamNum(:,:,:) = -999

      !===============================
      ! Do calculations every 24 hours
      !===============================
      IF (this%min_time == nhms ) THEN

         CALL Acquire_GMIspecies(rc)

        ! if not doing CO, use CO from input file (already in ppbv)    
        !else convert CO from mol/mol (GEOS-5) to ppbv in parameterization 
        if  ( .NOT. w_c%reg%doing_CO ) then
           STT(:,:,1:km,1) = BBIJ(:,:,:,9) 
        else
           STT(:,:,1:km,1) = this%CO_24h(:,:,1:km)*1.0d9   
        end if

        ! Same for CH4
        if  ( .NOT. w_c%reg%doing_CH4 ) then
            CH4VAR(:,:,:) = BBIJ(:,:,:,10) 
        else
            CH4VAR(:,:,:) = this%CH4_24h(:,:,:)*1.0d9    
        endif

        ! convert Spec Humidity (kg vapor/kg moist air) to ppmv 
        Qppmv(:,:,1:km) = this%SpecHum_24h(:,:,1:km)*mw_moist_air(:,:,1:km)/mwtWat*1.0d6 

        ! test with water vapor from GMI on GMI levels run (in ppmv)
        ! Qppmv(:,:,1:km) = this%qIn(:,:,1:km)*1.0d6 

        ! H2O adjustment to account for H2O bias in CCM in comparison with AIRS	
! yelshorbany, oct, 27, 2014, change the "print" statement below to OFF and CALL adjustH2OMER()
        IF (MAPL_AM_I_ROOT()) print *, 'H2O adjustment TURNED OFF'
!       IF (MAPL_AM_I_ROOT()) print *, 'H2O adjustment TURNED ON'
!	CALL adjustH20() !sas turned this back on
!       CALL adjustH2OMER() !sas version

        CALL runOHparam(this%OHparam, curMonth, ALBD, this%LPAUSE,    &
                 Qppmv, this%Press_24h, this%Temp_24h, &
                 this%Ravga, this%Ravgb, &
                 STT, BBIJ, o3up, &
                 CH4VAR, JDAY, RLON, RLAT,OHparamNum)
         ! runOH param returns OH in molec/cm-3 
         ! Reverse coordinates to pass into GEOS-5 environment
         IF (MAPL_AM_I_ROOT())  print *, ' After runOHparam: maxval BOH ', maxval(this%OHparam%BOH(:,:,:))

         ! Apply reaction rate ratio adjustment, still on GMI levels
         
         this%OHparam%BOH(:,:,:)= this%rxnRateRatio(:,:,:)*this%OHparam%BOH(:,:,:)

         ! Set BOH to 0 above the tropopause
         ! now set it to fullchem value instead of 0 -sas
         ! also do it above 300hPa, not just above tropopause -sas
         DO jx = j1, j2
            DO ix = i1, i2
       !        WHERE(this%Press_24h(ix,jx,:) < this%troppresIn(ix,jx)) this%OHparam%BOH(ix,jx,:)=0.
       !        WHERE(this%Press_24h(ix,jx,:) < this%troppresIn(ix,jx)) this%OHparam%BOH(ix,jx,:)= gmiOH(ix,jx,:) !-sas
                WHERE(this%Press_24h(ix,jx,:) < 300) this%OHparam%BOH(ix,jx,:)= gmiOH(ix,jx,:) !sas
            END DO
         END DO

         ! reverse levels into GEOS-5 format 
         ! MANYIN - OH is stored in data3D and used by CH4 and CO; it is not transported, and
         !          at the next change of day it will be overwritten
         w_c%qa(nbeg)%data3D(:,:,1:km) = this%OHparam%BOH(:,:,km:1:-1)  

         ! Ouput input variables for debugging
         call getExportStateVariables(rc)

         this%Ravga = 0.0d0
         this%Ravgb = 0.0d0

      ENDIF
  
      CALL pmaxmin('OH: OH Conc  ',       w_c%qa(nbeg)%data3D(i1:i2,j1:j2,1:km), qmin, qmax, iXj, km,  1. )

  !    CALL pmaxmin('fulchem OH: OH Conc  ',       gmiOH(i1:i1,j1:j1,1:km), qmin, qmax, iXj, km, 1. ) !sas


!  Housekeeping
!  ------------

      DEALLOCATE(ndDry, ndWet, o3up, STT, CH4VAR, Qppmv,RLON, RLAT, BBIJ, var3d, ALBD,PL, pe,p,OHparamNum)
      DEALLOCATE(mass, gridBoxThickness, solarZenithAngle, tropPress, optDepth, mw_moist_air)
      DEALLOCATE (gmiOH) !-sas

   
      DEALLOCATE(ier)

      RETURN

!---------------------------------------------------------------------------
   CONTAINS
!---------------------------------------------------------------------------
!BOP
      SUBROUTINE adjustH20()
! !DESCRIPTION:
! Adjust Qppmv based on CCM/T1R1 vs AIRS comparison ratios 
!
!EOP
!----------------------------------------------------------------------------
!BOC
        ! Adjust vapor based on bias wrt to AIRS data (2002-2009)
        !----------------------------------------------------------------------------
        ! ****** DJF months ********************************************************
        !----------------------------------------------------------------------------
        IF (curMonth .EQ. 1. .OR. curMonth .EQ. 2. .OR. curMonth .EQ. 12. ) THEN
            DO jx = j1, j2
               DO ix = i1, i2
                   DO k = 1,km
                       !----------------------------------------------------------------------------
                       ! NH  
                       !---------------------------------------------------------------------------- 
                       IF (RLAT(ix,jx) .GE. 0.) THEN
                           ! 800 hPa to 600 hPa level	 
                           IF (this%press_24h(ix,jx,k) .LE. 800 .AND. this%press_24h(ix,jx,k) .GT. 600) THEN
                                Qppmv(ix,jx,k) = (0.8)*this%SpecHum_24h(ix,jx,k)*mw_moist_air(ix,jx,k)/mwtWat*1.0d6 
                           ENDIF
                           ! 600 hPa to 400 hPa level	 
                           IF (this%press_24h(ix,jx,k) .LE. 600 .AND. this%press_24h(ix,jx,k) .GT. 400) THEN
                                Qppmv(ix,jx,k) = (0.6)*this%SpecHum_24h(ix,jx,k)*mw_moist_air(ix,jx,k)/mwtWat*1.0d6 
                           ENDIF
                           ! 400 hPa to 200 hPa level	 
                           IF (this%press_24h(ix,jx,k) .LE. 400 .AND. this%press_24h(ix,jx,k) .GT. 200) THEN
                                Qppmv(ix,jx,k) = (0.4)*this%SpecHum_24h(ix,jx,k)*mw_moist_air(ix,jx,k)/mwtWat*1.0d6 
                           ENDIF
                           ! above 200 hPa 	 
                           IF (this%press_24h(ix,jx,k) .LT. 200) THEN
                                Qppmv(ix,jx,k) = (1-0.45)*this%SpecHum_24h(ix,jx,k)*mw_moist_air(ix,jx,k)/mwtWat*1.0d6 
                           ENDIF
                       ENDIF
                       !----------------------------------------------------------------------------
                       ! SH
                       !----------------------------------------------------------------------------  
                       IF (RLAT(ix,jx) .LT. 0.) THEN
                           ! 800 hPa to 300 hPa level	 
                           IF (this%press_24h(ix,jx,k) .LE. 800 .AND. this%press_24h(ix,jx,k) .GT. 300) THEN
                                Qppmv(ix,jx,k) = (0.7)*this%SpecHum_24h(ix,jx,k)*mw_moist_air(ix,jx,k)/mwtWat*1.0d6 
                           ENDIF
                           ! above 300 hPa 	 
                           IF (this%press_24h(ix,jx,k) .LT. 300) THEN
                                Qppmv(ix,jx,k) = (1-0.45)*this%SpecHum_24h(ix,jx,k)*mw_moist_air(ix,jx,k)/mwtWat*1.0d6 
                           ENDIF
                       ENDIF       
                ! E141 run version 	
                IF (.FALSE.) THEN   ! NH >  60N	
                       IF (RLAT(ix,jx) .GT. 60.) THEN
                           ! 800 hPa to 500 hPa level	 
                           IF (this%press_24h(ix,jx,k) .LE. 800 .AND. this%press_24h(ix,jx,k) .GT. 500) THEN
                                Qppmv(ix,jx,k) = (1.1)*this%SpecHum_24h(ix,jx,k)*mw_moist_air(ix,jx,k)/mwtWat*1.0d6 
                           ENDIF
                           ! 500 hPa to 300 hPa level	 
                           IF (this%press_24h(ix,jx,k) .LE. 500 .AND. this%press_24h(ix,jx,k) .GT. 300) THEN
                                Qppmv(ix,jx,k) = (1.2)*this%SpecHum_24h(ix,jx,k)*mw_moist_air(ix,jx,k)/mwtWat*1.0d6 
                           ENDIF
                           ! above 300 hPa 	 
                           IF (this%press_24h(ix,jx,k) .LT. 300) THEN
                                Qppmv(ix,jx,k) = (0.90)*this%SpecHum_24h(ix,jx,k)*mw_moist_air(ix,jx,k)/mwtWat*1.0d6 
                           ENDIF
                       ENDIF
                       ! NH from 30N to 60N 
                       IF (RLAT(ix,jx) .GE. 30 .AND. RLAT(ix,jx) .LT. 60) THEN
                           ! 800 hPa to 700 hPa level	 
                           IF (this%press_24h(ix,jx,k) .LE. 800 .AND. this%press_24h(ix,jx,k) .GT. 700) THEN
                                Qppmv(ix,jx,k) = (0.85)*this%SpecHum_24h(ix,jx,k)*mw_moist_air(ix,jx,k)/mwtWat*1.0d6 
                           ENDIF
                           ! 700 hPa to 400 hPa level	 
                           IF (this%press_24h(ix,jx,k) .LE. 700 .AND. this%press_24h(ix,jx,k) .GT. 400) THEN
                                Qppmv(ix,jx,k) = (0.65)*this%SpecHum_24h(ix,jx,k)*mw_moist_air(ix,jx,k)/mwtWat*1.0d6 
                           ENDIF
                           ! 400 hPa to 300 hPa level	 
                           IF (this%press_24h(ix,jx,k) .LE. 400 .AND. this%press_24h(ix,jx,k) .GT. 300) THEN
                                Qppmv(ix,jx,k) = (0.85)*this%SpecHum_24h(ix,jx,k)*mw_moist_air(ix,jx,k)/mwtWat*1.0d6 
                           ENDIF
                           ! above 300 hPa 	 
                           IF (this%press_24h(ix,jx,k) .LT. 300) THEN
                                Qppmv(ix,jx,k) = (0.60)*this%SpecHum_24h(ix,jx,k)*mw_moist_air(ix,jx,k)/mwtWat*1.0d6 
                           ENDIF
                       ENDIF
                       ! from 30N to 30S 
                       IF (RLAT(ix,jx) .GE. -30 .AND. RLAT(ix,jx) .LT. 30) THEN
                           ! 800 hPa to 600 hPa level	 
                           IF (this%press_24h(ix,jx,k) .LE. 800 .AND. this%press_24h(ix,jx,k) .GT. 600) THEN
                                Qppmv(ix,jx,k) = (0.75)*this%SpecHum_24h(ix,jx,k)*mw_moist_air(ix,jx,k)/mwtWat*1.0d6 
                           ENDIF
                           ! 600 hPa to 500 hPa level	 
                           IF (this%press_24h(ix,jx,k) .LE. 600 .AND. this%press_24h(ix,jx,k) .GT. 500) THEN
                                Qppmv(ix,jx,k) = (0.6)*this%SpecHum_24h(ix,jx,k)*mw_moist_air(ix,jx,k)/mwtWat*1.0d6 
                           ENDIF
                           ! 500 hPa to 300 hPa level	 
                           IF (this%press_24h(ix,jx,k) .LE. 500 .AND. this%press_24h(ix,jx,k) .GT. 300) THEN
                                Qppmv(ix,jx,k) = (0.20)*this%SpecHum_24h(ix,jx,k)*mw_moist_air(ix,jx,k)/mwtWat*1.0d6 
                           ENDIF
                           ! above 300 hPa 	 
                           IF (this%press_24h(ix,jx,k) .LT. 300) THEN
                                Qppmv(ix,jx,k) = (0.10)*this%SpecHum_24h(ix,jx,k)*mw_moist_air(ix,jx,k)/mwtWat*1.0d6 
                           ENDIF
                       ENDIF
                ENDIF
                        ! >60S no correction ratio = 1
                   END DO
                END DO
            END DO  
        ENDIF
        !----------------------------------------------------------------------------
        ! ****** MAM months ********************************************************
        !----------------------------------------------------------------------------
        IF (curMonth .EQ. 3. .OR. curMonth .EQ. 4. .OR. curMonth .EQ. 5. ) THEN
            DO jx = j1, j2
               DO ix = i1, i2
                   DO k = 1,km
                       !----------------------------------------------------------------------------
                       ! NH  
                       !---------------------------------------------------------------------------- 
                       IF (RLAT(ix,jx) .GE. 0.) THEN
                           ! 800 hPa to 500 hPa level	 
                           IF (this%press_24h(ix,jx,k) .LE. 800 .AND. this%press_24h(ix,jx,k) .GT. 500) THEN
                                Qppmv(ix,jx,k) = (0.8)*this%SpecHum_24h(ix,jx,k)*mw_moist_air(ix,jx,k)/mwtWat*1.0d6 
                           ENDIF
                           ! 500 hPa to 300 hPa level	 
                           IF (this%press_24h(ix,jx,k) .LE. 500 .AND. this%press_24h(ix,jx,k) .GT. 300) THEN
                                Qppmv(ix,jx,k) = (0.6)*this%SpecHum_24h(ix,jx,k)*mw_moist_air(ix,jx,k)/mwtWat*1.0d6 
                           ENDIF
                           ! 300 hPa to 200 hPa level	 
                           IF (this%press_24h(ix,jx,k) .LE. 300 .AND. this%press_24h(ix,jx,k) .GT. 200) THEN
                                Qppmv(ix,jx,k) = (0.4)*this%SpecHum_24h(ix,jx,k)*mw_moist_air(ix,jx,k)/mwtWat*1.0d6 
                           ENDIF
                           ! above 200 hPa 	 
                           IF (this%press_24h(ix,jx,k) .LT. 200) THEN
                                Qppmv(ix,jx,k) = (0.2)*this%SpecHum_24h(ix,jx,k)*mw_moist_air(ix,jx,k)/mwtWat*1.0d6 
                           ENDIF
                       ENDIF
                       !----------------------------------------------------------------------------
                       ! SH
                       !----------------------------------------------------------------------------  
                       IF (RLAT(ix,jx) .LT. 0.) THEN
                           ! 800 hPa to 300 hPa level	 
                           IF (this%press_24h(ix,jx,k) .LE. 800 .AND. this%press_24h(ix,jx,k) .GT. 300) THEN
                                Qppmv(ix,jx,k) = (0.7)*this%SpecHum_24h(ix,jx,k)*mw_moist_air(ix,jx,k)/mwtWat*1.0d6 
                           ENDIF
                           ! above 300 hPa 	 
                           IF (this%press_24h(ix,jx,k) .LT. 300) THEN
                                Qppmv(ix,jx,k) = (1-0.45)*this%SpecHum_24h(ix,jx,k)*mw_moist_air(ix,jx,k)/mwtWat*1.0d6 
                           ENDIF
                       ENDIF
                   END DO
                END DO
            END DO  
        ENDIF

        !----------------------------------------------------------------------------
        ! ****** JJA months ********************************************************
        !----------------------------------------------------------------------------
        IF (curMonth .GE. 6 .AND. curMonth .LE. 8) THEN
            DO jx = j1, j2
               DO ix = i1, i2
                   DO k = 1,km
                       !----------------------------------------------------------------------------
                       ! NH from O to 60N
                       !----------------------------------------------------------------------------	
                       IF (RLAT(ix,jx) .GE. 0 .AND. RLAT(ix,jx) .LT. 60) THEN
                           ! 700 hPa to 600 hPa level	 
                           IF (this%press_24h(ix,jx,k) .LE. 700 .AND. this%press_24h(ix,jx,k) .GT. 600) THEN
                                Qppmv(ix,jx,k) = (1-0.30)*this%SpecHum_24h(ix,jx,k)*mw_moist_air(ix,jx,k)/mwtWat*1.0d6 
                           ENDIF
                           ! 600 hPa to 500 hPa level	 
                           IF (this%press_24h(ix,jx,k) .LE. 600 .AND. this%press_24h(ix,jx,k) .GT. 500) THEN
                                Qppmv(ix,jx,k) = (1-0.40)*this%SpecHum_24h(ix,jx,k)*mw_moist_air(ix,jx,k)/mwtWat*1.0d6 
                           ENDIF
                           ! above 500 hPa 	 
                           IF (this%press_24h(ix,jx,k) .LT. 500) THEN
                                Qppmv(ix,jx,k) = (1-0.60)*this%SpecHum_24h(ix,jx,k)*mw_moist_air(ix,jx,k)/mwtWat*1.0d6 
                           ENDIF
                           
                       ENDIF
                   END DO
                END DO
            END DO  
        ENDIF
        !----------------------------------------------------------------------------
        ! ****** SON months ********************************************************
        !----------------------------------------------------------------------------
        IF (curMonth .EQ. 9. .OR. curMonth .EQ. 10. .OR. curMonth .EQ. 11. ) THEN
            DO jx = j1, j2
               DO ix = i1, i2
                   DO k = 1,km
                       !----------------------------------------------------------------------------
                       ! NH  
                       !---------------------------------------------------------------------------- 
                       IF (RLAT(ix,jx) .GE. 0.) THEN
                           ! 800 hPa to 500 hPa level	 
                           IF (this%press_24h(ix,jx,k) .LE. 800 .AND. this%press_24h(ix,jx,k) .GT. 500) THEN
                                Qppmv(ix,jx,k) = (0.8)*this%SpecHum_24h(ix,jx,k)*mw_moist_air(ix,jx,k)/mwtWat*1.0d6 
                           ENDIF
                           ! 500 hPa to 300 hPa level	 
                           IF (this%press_24h(ix,jx,k) .LE. 500 .AND. this%press_24h(ix,jx,k) .GT. 300) THEN
                                Qppmv(ix,jx,k) = (0.6)*this%SpecHum_24h(ix,jx,k)*mw_moist_air(ix,jx,k)/mwtWat*1.0d6 
                           ENDIF
                           ! 300 hPa to 200 hPa level	 
                           IF (this%press_24h(ix,jx,k) .LE. 300 .AND. this%press_24h(ix,jx,k) .GT. 200) THEN
                                Qppmv(ix,jx,k) = (0.4)*this%SpecHum_24h(ix,jx,k)*mw_moist_air(ix,jx,k)/mwtWat*1.0d6 
                           ENDIF
                           ! above 200 hPa 	 
                           IF (this%press_24h(ix,jx,k) .LT. 200) THEN
                                Qppmv(ix,jx,k) = (0.2)*this%SpecHum_24h(ix,jx,k)*mw_moist_air(ix,jx,k)/mwtWat*1.0d6 
                           ENDIF
                       ENDIF
                       !----------------------------------------------------------------------------
                       ! SH
                       !----------------------------------------------------------------------------  
                       IF (RLAT(ix,jx) .LT. 0.) THEN
                           ! 800 hPa to 300 hPa level	 
                           IF (this%press_24h(ix,jx,k) .LE. 800 .AND. this%press_24h(ix,jx,k) .GT. 300) THEN
                                Qppmv(ix,jx,k) = (0.8)*this%SpecHum_24h(ix,jx,k)*mw_moist_air(ix,jx,k)/mwtWat*1.0d6 
                           ENDIF
                           ! above 300 hPa 	 
                           IF (this%press_24h(ix,jx,k) .LT. 300) THEN
                                Qppmv(ix,jx,k) = (1-0.45)*this%SpecHum_24h(ix,jx,k)*mw_moist_air(ix,jx,k)/mwtWat*1.0d6 
                           ENDIF
                       ENDIF
                   END DO
                END DO
            END DO  
        ENDIF

      
      END SUBROUTINE adjustH20
!EOC
!----------------------------------------------------------------------------

!!!
!!!  MANYIN - IF THIS ROUTINE IS TO BE USED, IT MUST BE REFACTORED TO MAKE
!!!           THE 7-LEVEL FILE INTO MULTIPLE SINGLE-LEVEL FILES
!!!
!!!!---------------------------------------------------------------------------
!!!!BOP
!!!      SUBROUTINE adjustH2OMER()
!!!! !DESCRIPTION:
!!!! Adjust Qppmv based on CCM/ACCMIP vs MERRA comparison ratios 
!!!!
!!!!EOP
!!!!----------------------------------------------------------------------------
!!!!BOC
!!!        ! Adjust vapor based on bias wrt to MERRA data (1997-2010)
!!!        !----------------------------------------------------------------------------
!!!
!!!        INTEGER :: nymd1, nhms1
!!!
!!!        REAL, ALLOCATABLE :: Qscale(:,:,:)
!!!        REAL, ALLOCATABLE :: Qtest(:,:,:)
!!!
!!!        !  Initialize local variables
!!!        !  --------------------------
!!! 
!!!        ! allocate local variables
!!!        allocate( Qscale(i1:i2,j1:j2,1:7),STAT=ier(38))
!!!        allocate( Qtest(i1:i2,j1:j2,1:km),STAT=ier(39))
!!!
!!!        nymd1 = 2000*10000 + MOD ( nymd, 10000 )  ! assumes 2000
!!!        nhms1 = 120000
!!!
!!!        var3d = 0.0d0
!!!
!!!! MANYIN:
!!!! Since this file is 7 levels, ExtData will not work
!!!! We could store as 7 individual 2D datasets, better for ExtData
!!!
!!!        ! read in Qratio: note that reader will flip to make lowest pressure 
!!!        ! 1st value of array
!!!        CALL Chem_UtilMPread ( this%Qscl_infile_name, 'Qratio', nymd1, nhms1, &
!!!                           i1, i2, 0, im, j1, j2, 0, jm, 7, &
!!!                           var3d=Qscale, cyclic=.true., &
!!!                           grid=w_c%grid_esmf )
!!! 
!!!        CALL pmaxmin('Qratio:  ', Qscale, qmin, qmax, iXj, 7,  1. )
!!!
!!!        Qtest = Qppmv
!!!        CALL pmaxmin('Qppmv orig:  ', Qtest(i1:i2,j1:j2,50), qmin, qmax, iXj, 1,  1. )
!!!
!!!        DO jx = j1, j2
!!!           DO ix = i1, i2
!!!              DO k = 1,km
!!!                 IF (this%press_24h(ix,jx,k) .LE. 300 .AND. this%press_24h(ix,jx,k) .GT. 200) THEN 
!!!                    Qppmv(ix,jx,k) = this%SpecHum_24h(ix,jx,k) * Qscale(ix,jx,1) * mw_moist_air(ix,jx,k)/mwtWat*1.0d6 
!!!                 ENDIF
!!!                 IF (this%press_24h(ix,jx,k) .LE. 400 .AND. this%press_24h(ix,jx,k) .GT. 300) THEN 
!!!                    Qppmv(ix,jx,k) = this%SpecHum_24h(ix,jx,k) * Qscale(ix,jx,2) * mw_moist_air(ix,jx,k)/mwtWat*1.0d6 
!!!                 ENDIF
!!!                 IF (this%press_24h(ix,jx,k) .LE. 500 .AND. this%press_24h(ix,jx,k) .GT. 400) THEN 
!!!                    Qppmv(ix,jx,k) = this%SpecHum_24h(ix,jx,k) * Qscale(ix,jx,3) * mw_moist_air(ix,jx,k)/mwtWat*1.0d6 
!!!                 ENDIF
!!!                 IF (this%press_24h(ix,jx,k) .LE. 600 .AND. this%press_24h(ix,jx,k) .GT. 500) THEN 
!!!                    Qppmv(ix,jx,k) = this%SpecHum_24h(ix,jx,k) * Qscale(ix,jx,4) * mw_moist_air(ix,jx,k)/mwtWat*1.0d6 
!!!                 ENDIF
!!!                 IF (this%press_24h(ix,jx,k) .LE. 700 .AND. this%press_24h(ix,jx,k) .GT. 600) THEN 
!!!                    Qppmv(ix,jx,k) = this%SpecHum_24h(ix,jx,k) * Qscale(ix,jx,5) * mw_moist_air(ix,jx,k)/mwtWat*1.0d6 
!!!                 ENDIF
!!!                 IF (this%press_24h(ix,jx,k) .LE. 800 .AND. this%press_24h(ix,jx,k) .GT. 700) THEN 
!!!                    Qppmv(ix,jx,k) = this%SpecHum_24h(ix,jx,k) * Qscale(ix,jx,6) * mw_moist_air(ix,jx,k)/mwtWat*1.0d6 
!!!                 ENDIF
!!!                 IF (this%press_24h(ix,jx,k) .LE. 900 .AND. this%press_24h(ix,jx,k) .GT. 800) THEN 
!!!                    Qppmv(ix,jx,k) = this%SpecHum_24h(ix,jx,k) * Qscale(ix,jx,7) * mw_moist_air(ix,jx,k)/mwtWat*1.0d6 
!!!                 ENDIF
!!!
!!!              END DO
!!!           END DO
!!!        END DO
!!!
!!!        Qtest = Qppmv
!!!        CALL pmaxmin('Qppmv scaled:  ', Qtest(i1:i2,j1:j2,50), qmin, qmax, iXj, 1,  1. )
!!!
!!!        DEALLOCATE( Qscale, Qtest)
!!!
!!!      END SUBROUTINE adjustH2OMER
!!!!EOC
!!!!----------------------------------------------------------------------------


!---------------------------------------------------------------------------
!BOP
      SUBROUTINE calcRxnRateFix()

      INTEGER, SAVE :: savedDay = negativeOne
      INTEGER       :: currDay, k1strat
      REAL,    ALLOCATABLE :: k1_new(:,:,:)           !   bottom-up
      REAL,    ALLOCATABLE :: k1_old(:,:,:)           !   bottom-up
      REAL,    ALLOCATABLE :: k2_new(:,:,:)           !   bottom-up
      REAL,    ALLOCATABLE :: k2_old(:,:,:)           !   bottom-up
      REAL,    ALLOCATABLE :: k3_new(:,:,:)           !   bottom-up
      REAL,    ALLOCATABLE :: k3_old(:,:,:)           !   bottom-up
      REAL,    ALLOCATABLE :: N2molec(:,:,:)          !   bottom-up
      REAL,    ALLOCATABLE :: O2molec(:,:,:)          !   bottom-up
      REAL,    ALLOCATABLE :: H2Omolec(:,:,:)         !   bottom-up
      REAL,    ALLOCATABLE :: ADENSITY(:,:,:)         !   bottom-up
      REAL,    ALLOCATABLE :: newarr(:,:,:)           !   bottom-up
      REAL,    ALLOCATABLE :: oldarr(:,:,:)           !   bottom-up
      REAL,    ALLOCATABLE :: oldarr2(:,:,:)          !   bottom-up

      INTEGER, ALLOCATABLE :: ier(:)
! !DESCRIPTION:
! Calculate BOH adjustment ratio based on new and old reaction rates
! Use monthly mean water vapor, and N2, O2 (based on AIR DENSITY)
!
!
!EOP
!----------------------------------------------------------------------------
!BOC
 CHARACTER(LEN=255) :: MyName
  INTEGER :: ic
  currDay = MOD(nymd,100)
  rc=0
  MyName="calcRateFix"
  allocate(ier(199))
  ier(:) = 0
  allocate(k1_new(i1:i2,j1:j2,1:km),STAT=ier(10))
  allocate(k2_new(i1:i2,j1:j2,1:km),STAT=ier(11))
  allocate(k3_new(i1:i2,j1:j2,1:km),STAT=ier(12))
  allocate(k1_old(i1:i2,j1:j2,1:km),STAT=ier(13))
  allocate(k2_old(i1:i2,j1:j2,1:km),STAT=ier(14))
  allocate(k3_old(i1:i2,j1:j2,1:km),STAT=ier(15))
  allocate(N2molec(i1:i2,j1:j2,1:km),STAT=ier(16))
  allocate(O2molec(i1:i2,j1:j2,1:km),STAT=ier(17))
  allocate(H2Omolec(i1:i2,j1:j2,1:km),STAT=ier(18))
  allocate(ADENSITY(i1:i2,j1:j2,1:km),STAT=ier(19))
  allocate(newarr(i1:i2,j1:j2,1:km),STAT=ier(20))
  allocate(oldarr(i1:i2,j1:j2,1:km),STAT=ier(21))
  allocate(oldarr2(i1:i2,j1:j2,1:km),STAT=ier(22))
      IF ( ANY(ier(:) /= 0)) THEN
         PRINT*, TRIM(myname),': Failed to allocate work space.'
         rc = 21
         RETURN
      END IF
      ier(:) = 0

! Execute updates at 0:00 UTC and/or on first pass
! ------------------------------------------------

   IF (currDay /= savedDay) THEN
      var3d = 0.0d0



      !  reverse coordinates into GMI levels and convert number density molec/m3 to molec/cm3
      ADENSITY(:,:,1:km) = ndWet(:,:,km:1:-1)/(1E6)
      ! convert Spec Humidity (on GMI levels) (kg H20/kg air) to VMR and them to molec/cm3
      H2Omolec(:,:,:) = this%SpecHum_24h(:,:,:)*mw_moist_air(:,:,:)/mwtWat*ADENSITY(:,:,:)
      N2molec(:,:,:) = ADENSITY*0.78
      O2molec(:,:,:) = ADENSITY*0.21

      ! Calculate the change due to change in reaction 
      !  TempIn is in Kelvin
      ! R1 : O1D + H2O -> 2OH 
      ! R2 : O1D + N2  -> O3P + N2
      ! R3 : O1D + O2 -> O + O2

      ! IF (	min (this%Temp_24h(:,:,:)-273) .LE. 0) print *, 'Error; Check your Temperature Array in Restart'

      !Temperature on GMI levels already surface 1-> 72 (surface --> toa)
      ! rate constants in  units of cm3/molec/s
      ! JPL11
      k1_new(:,:,:) = (1.63E-10)*EXP(60./(this%Temp_24h(:,:,:)))
      ! JPL2000
      k1_old(:,:,:) = (2.2E-10)*EXP(100/(this%Temp_24h(:,:,:)))

      ! JPL11
      k2_new(:,:,:) =(2.15E-11)*EXP(110/(this%Temp_24h(:,:,:)))
      !Atkins 1997
      k2_old(:,:,:) =  (2.8E-11)*EXP(107/(this%Temp_24h(:,:,:)))

      !JPL11
      k3_new(:,:,:) =(3.3E-11)*EXP(55/(this%Temp_24h(:,:,:)))
      !Atkins 1997
      k3_old(:,:,:) = (3.2E-11)*EXP(67/(this%Temp_24h(:,:,:)))

      newarr(:,:,:) = 2*(H2Omolec(:,:,:)*k1_new(:,:,:))/ &
         (H2Omolec(:,:,:)*(k1_new(:,:,:))+N2molec(:,:,:)*(k2_new(:,:,:))+O2molec(:,:,:)*k3_new(:,:,:))
      
      oldarr2(:,:,:) = (H2Omolec(:,:,:)*(k1_old(:,:,:))+N2molec(:,:,:)*(k2_old(:,:,:))+O2molec(:,:,:)*(k3_old(:,:,:)))
      
      oldarr(:,:,:) = 2*(H2Omolec(:,:,:)*k1_old(:,:,:))/oldarr2(:,:,:)
        
        
        ! Adjust for a O1D reaction rate changes
       DO jx = j1, j2
            DO ix = i1, i2
               DO k = 1,km
                    IF(oldarr2(ix,jx,k) == 0) THEN
                       oldarr(ix,jx,k) = newarr(ix,jx,k)
                    ELSE 
                       oldarr(ix,jx,k) = 2*(H2Omolec(ix,jx,k)*k1_old(ix,jx,k))/oldarr2(ix,jx,k)
                    END IF
                    IF(oldarr(ix,jx,k) == 0) THEN
                       this%rxnRateRatio(ix,jx,k) = 1.0
                    ELSE 
                       this%rxnRateRatio(ix,jx,k) = newarr(ix,jx,k)/oldarr(ix,jx,k)
                    END IF
            
               END DO
            END DO
         END DO

      CALL pmaxmin('Phot. Adjust Ratio', this%rxnRateRatio  , qmin, qmax, iXj,   km, 1. )
      
      DEALLOCATE(k1_new, k2_new, k3_new, k1_old, k2_old,k3_old,N2molec,O2molec, H2Omolec, ADENSITY)
      DEALLOCATE(newarr, oldarr)
      savedDay = currDay
      END IF
      DEALLOCATE(ier)

      END SUBROUTINE calcRxnRateFix
!EOC
!----------------------------------------------------------------------------
!---------------------------------------------------------------------------
!BOP
      SUBROUTINE calcTroppLevel()

      INTEGER, SAVE :: savedDay = negativeOne
      INTEGER       :: currDay, k1strat
      LOGICAL       :: tropo(im,jm,km)

!
! !DESCRIPTION:
! Determine the tropopause level of each grid box.
!
!EOP
!----------------------------------------------------------------------------
!BOC

!      currDay = nymd /1000000
       currDay = MOD(nymd, 100)
      IF (currDay /= savedDay) THEN

        this%LPAUSE = 0   
   
        tropo(:,:,:) = .FALSE.
         
        DO jx = j1, j2
               DO ix = i1, i2
                    WHERE(this%Press_24h(ix,jx,:) >= this%troppresIn(ix,jx)) tropo(ix,jx,:)=.TRUE.
               END DO
        END DO
        ! itrop pres index
        k1strat = km
        
           DO jx = j1, j2
               DO ix = i1, i2
                   DO k = km,1, -1
                    IF (tropo(ix,jx,k)) EXIT
                   END DO
                   !this%LPAUSE(ix,jx) = k
                END DO
            END DO   
        
        
! Start from the top of atmosphere
         DO k = km,2, -1       
            DO jx = j1, j2
               DO ix = i1, i2
                 IF ((this%troppresIn(ix,jx) .GE. this%Press_24h(ix,jx,k))) THEN 
                      this%LPAUSE(ix,jx) = k
                 END IF
               END DO
            END DO
         END DO
        
         savedDay = currDay
      END IF

      RETURN 

      END SUBROUTINE calcTroppLevel
!EOC
!BOP
!
! !ROUTINE:  Acquire_GMIspecies
!
! !INTERFACE:

  SUBROUTINE Acquire_GMIspecies(rc)

  INTEGER, INTENT(OUT) :: rc
!
! !DESCRIPTION:
! Reads the GMI species concentrations every 24 hours
! from an atmosphere with and without aerosols.
!
!EOP
!---------------------------------------------------------------------------

  CHARACTER(LEN=255) :: MyName
  INTEGER :: ic

  rc=0
  MyName="Acquire_GMIspecies"

! Execute updates at 0:00 UTC and/or on first pass
! ------------------------------------------------
   IF(nhms == 0) THEN

    IF( MAPL_AM_I_ROOT() ) THEN
     PRINT *, ' '
     PRINT *,'OH_GridCompRun: Obtaining GMI species concentration from monthly means ...'
    END IF


      !=========================================
      ! Reading species concentrations from file
      !=========================================
      ! BBIJ stores species concentration for:
      ! 1. NOt             2. ALK4          3. Isoprene
      ! 4. Acetone         5. Propene       6. Propane
      ! 7. Ethane          8. O3            9. CO (not used)
      !10. CH4 (not used)
      BBIJ = 0.0d0
      DO ic = 1, 15                             ! 14 if reading CO too, 15 if CH4

      ! This provides top-down field (GEOS-5 levels)
         CALL MAPL_GetPointer(impChem, ptr3d, gmiSpeciesName3D(ic), __RC__ )

! MANYIN - elim the loop
      ! Reverse vertical to get back to GMI and OH param level convention (1 at surface , 72 at top)
         DO k=1,km
            kReverse=km-k+1
            IF (ic .LE. 6) THEN
               BBIJ(:,:,k,1) = BBIJ(:,:,k,1) + ptr3d(:,:,kReverse)
               ! ADD in N2O5 the second time into the NOt calculation
               IF(ic .EQ. 3) THEN 
                 BBIJ(:,:,k,1) = BBIJ(:,:,k,1) + ptr3d(:,:,kReverse)
               ENDIF
            ELSE
               BBIJ(:,:,k,ic-5)=ptr3d(:,:,kReverse)
            END IF
         END DO
      END DO

      BBIJ(:,:,:, 1) = BBIJ(:,:,:, 1) * 1.D12     ! NOt  pptv
      BBIJ(:,:,:, 2) = BBIJ(:,:,:, 2) * 1.D12     ! ALK4 (alkanes pptv)
      BBIJ(:,:,:, 3) = BBIJ(:,:,:, 3) * 1.D12     ! ISOP (isoprene pptv)
      BBIJ(:,:,:, 4) = BBIJ(:,:,:, 4) * 1.D12     ! ACET (acetone pptv)
      BBIJ(:,:,:, 5) = BBIJ(:,:,:, 5) * 1.D12     ! PRPE (propene pptv)
      BBIJ(:,:,:, 6) = BBIJ(:,:,:, 6) * 1.D12     ! C3H8 (propane pptv)
      BBIJ(:,:,:, 7) = BBIJ(:,:,:, 7) * 1.D12     ! C2H6 (ethane pptv)
      BBIJ(:,:,:, 8) = BBIJ(:,:,:, 8) * 1.D9      ! O3 (ozone ppbv)
      BBIJ(:,:,:, 9) = BBIJ(:,:,:, 9) * 1.D9      ! CO (co ppbv)
      BBIJ(:,:,:,10) = BBIJ(:,:,:,10) * 1.D9      ! CH4 (ppbv)

      ! Also read in full-chem OH for use in the stratosphere -sas
      ic = 16
      CALL MAPL_GetPointer(impChem, ptr3d, gmiSpeciesName3D(ic), __RC__ )

      ! convert OH from mol/mol to molec/cm3
      ! (assume that GMI output is volume mixing ratio with respect to dry air)  MANYIN
      var3d(:,:,:) = ptr3d(:,:,:) * ndDry(:,:,:) * 1.E-6

      ! Unreverse to get back to GMI and OH param level convention (1 at surface , 72 at top)
      DO k=1,km
         kReverse=km-k+1
         gmiOH(:,:,k)=var3d(:,:,kReverse) !molec/cm3
      END DO

      IF (MAPL_AM_I_ROOT())  print *, ' maxval gmiOH ', maxval(gmiOH(:,:,:))

      ! Determine the overhead O3 column
      !---------------------------------
      ! convert (ppbv) to ozone (DU) for each box
      ! O3 * air density * box thickness * conversion to DU
      ! ppbv *  kg/m3 * m * 1 mol air/0.02897 kg * 6.023e23 molec/mol air * 1DU/(2.69e20 molec/m2)
! MANYIN: previous version did not vertically reverse gridBoxThickness
!   --- but I am testing with it, to reproduce Yasins results better
!     var3d(:,:,:) = BBIJ(:,:,:,8)*1.0d-9*rhoa(:,:,km:1:-1)*gridBoxThickness(:,:,:)/0.02897*6.023e23/ (2.69d20)
      ! (The rhoWet just cancels out with the same term used in computing gridBoxThickness)
      var3d(:,:,1:km) = BBIJ(:,:,:,8)*1.0d-9*rhoWet(:,:,km:1:-1)*gridBoxThickness(:,:,km:1:-1)/0.02897*6.023E23/(2.69E20)

! MANYIN  var3d and o3up are bottom-up;  this code is slightly wrong and very inefficient
      o3up(:,:,:) = 0.0d0
      o3up(:,:,km) = var3d(:,:,km)
      DO ic = 1, km-1
         DO k = ic+1, km
            o3up(:,:,ic) = o3up(:,:,ic) + var3d(:,:,k)
         END DO
      END DO
! MANYIN  instead do this, for 'overhead column including the current gridbox':
!     o3up(:,:,km) = var3d(:,:,km)
!     DO k = km-1,1,-1
!           o3up(:,:,k) = var3d(:,:,k) + o3up(:,:,k+1)
!     END DO
! MANYIN  or do this, for 'overhead column not including the current gridbox':
!     Note: the O3 above the top gridbox should not really be zero !
!     o3up(:,:,km) = 0.0d0
!     DO k = km-1,1,-1
!           o3up(:,:,k) = var3d(:,:,k+1) + o3up(:,:,k+1)
!     END DO

  END IF

  RETURN

  END SUBROUTINE Acquire_GMIspecies
!EOC
!------------------------------------------------------------------------------------

 !---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  getExportStateVariables
!
! !INTERFACE:
!
      SUBROUTINE getExportStateVariables(rc)

      IMPLICIT NONE

      INTEGER, INTENT(OUT) :: rc
!
! !DESCRIPTION:
!
!  Find pointers to import and export states
!
!EOP
!---------------------------------------------------------------------------

      CHARACTER(LEN=255) :: Iam
      REAL :: qmin,qmax

      rc=0
      Iam="getExportStateVariables"

! MANYIN - These had included  ALLOC=.TRUE.  - unneeded since we check association before using

      CALL MAPL_GetPointer(expChem, exO3up,     'O3upIn',     __RC__)
      CALL MAPL_GetPointer(expChem, exNOy,      'NOyIn',      __RC__)
      CALL MAPL_GetPointer(expChem, exALK4,     'ALK4In',     __RC__)
      CALL MAPL_GetPointer(expChem, exIsoprene, 'IsopreneIn', __RC__)
      CALL MAPL_GetPointer(expChem, exAcetone,  'AcetoneIn',  __RC__)
      CALL MAPL_GetPointer(expChem, exPropene,  'PropeneIn',  __RC__)
      CALL MAPL_GetPointer(expChem, exPropane,  'PropaneIn',  __RC__)
      CALL MAPL_GetPointer(expChem, exEthane,   'EthaneIn',   __RC__)
      CALL MAPL_GetPointer(expChem, exOzone,    'OzoneIn',    __RC__)
      CALL MAPL_GetPointer(expChem, exSpecHum,  'SpecHumIn',  __RC__)
      CALL MAPL_GetPointer(expChem, exCO,       'COIn',       __RC__)
      CALL MAPL_GetPointer(expChem, exCH4,      'CH4In',      __RC__)
      CALL MAPL_GetPointer(expChem, exPress,    'PressIn',    __RC__)
      CALL MAPL_GetPointer(expChem, exTemp,     'TempIn',     __RC__)
      CALL MAPL_GetPointer(expChem, exTRPPS,    'TRPPSIn',    __RC__)
      CALL MAPL_GetPointer(expChem, exAlbedo,   'AlbedoIn',   __RC__)
      CALL MAPL_GetPointer(expChem, exRavga,    'RavgaIn',    __RC__)
      CALL MAPL_GetPointer(expChem, exRavgb,    'RavgbIn',    __RC__)
      CALL MAPL_GetPointer(expChem, exvarOut,   'varOut',     __RC__)
      

      ! flip these back to GEOS5 (toa-down) convention -sas
      IF (ASSOCIATED(exO3up))         exO3up(:,:,:) = O3up(:,:,km:1:-1)
      IF (ASSOCIATED(exNOy))           exNOy(:,:,:) = BBIJ(:,:,km:1:-1,1)
      IF (ASSOCIATED(exALK4))         exALK4(:,:,:) = BBIJ(:,:,km:1:-1,2)
      IF (ASSOCIATED(exIsoprene)) exIsoprene(:,:,:) = BBIJ(:,:,km:1:-1,3)
      IF (ASSOCIATED(exAcetone))   exAcetone(:,:,:) = BBIJ(:,:,km:1:-1,4)
      IF (ASSOCIATED(exPropene))   exPropene(:,:,:) = BBIJ(:,:,km:1:-1,5)      
      IF (ASSOCIATED(exPropane))   exPropane(:,:,:) = BBIJ(:,:,km:1:-1,6)
      IF (ASSOCIATED(exEthane))     exEthane(:,:,:) = BBIJ(:,:,km:1:-1,7)
      IF (ASSOCIATED(exOzone))       exOzone(:,:,:) = BBIJ(:,:,km:1:-1,8)
      IF (ASSOCIATED(exSpecHum))   exSpecHum(:,:,:) = Qppmv(i1:i2,j1:j2,km:1:-1)
      IF (ASSOCIATED(exCO))             exCO(:,:,:) = STT(:,:,km:1:-1,1)
      IF (ASSOCIATED(exCH4))           exCH4(:,:,:) = CH4VAR(:,:,km:1:-1)
      IF (ASSOCIATED(exPress))       exPress(:,:,:) = this%Press_24h(i1:i2,j1:j2,km:1:-1)
      IF (ASSOCIATED(exTemp))         exTemp(:,:,:) = this%Temp_24h(i1:i2,j1:j2,km:1:-1)
      IF (ASSOCIATED(exTRPPS))       exTRPPS(:,:)   = this%LPAUSE(i1:i2,j1:j2) * 1.0
      IF (ASSOCIATED(exAlbedo))     exAlbedo(:,:)   = ALBD(i1:i2,j1:j2)
      IF (ASSOCIATED(exRavga))       exRavga(:,:,:) = this%Ravga(i1:i2,j1:j2,km:1:-1)
      IF (ASSOCIATED(exRavgb))       exRavgb(:,:,:) = this%Ravgb(i1:i2,j1:j2,km:1:-1)
      IF (ASSOCIATED(exvarOut))     exvarOut(:,:,:) = OHparamNum(:,:,km:1:-1)
      RETURN

      END SUBROUTINE getExportStateVariables     


 END SUBROUTINE OH_GridCompRun
!EOP
!---------------------------------------------------------------------------

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  OH_GridCompFinalize --- Finalize OH_GridComp
!
! !INTERFACE:
!

   subroutine OH_GridCompFinalize ( this, w_c, impChem, expChem, &
                                      nymd, nhms, cdt, rc )

! !USES:

  IMPLICIT NONE

! !INPUT PARAMETERS:

   TYPE(Chem_Bundle), intent(in) :: w_c        ! Chemical tracer fields      
   INTEGER, INTENT(IN) :: nymd, nhms           ! time
   REAL,    INTENT(IN) :: cdt                  ! chemical timestep (secs)


! !OUTPUT PARAMETERS:

   TYPE(OH_GridComp), INTENT(INOUT) :: this   ! Grid Component
   TYPE(ESMF_State), INTENT(INOUT)  :: impChem  ! Import State
   TYPE(ESMF_State), INTENT(INOUT)  :: expChem  ! Export State
   INTEGER, INTENT(OUT) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 - 

! !DESCRIPTION: Finalizes the OH Parameterization Grid Component. Multiple instance
!               version.
!
! !REVISION HISTORY:
!
!  27Feb2008  da Silva  Introduced multiple instances
!
!EOP
!-------------------------------------------------------------------------

   CHARACTER(LEN=*), PARAMETER :: myname = 'OH_GridCompFinalize'
   INTEGER :: ios

   rc = 0
   DEALLOCATE ( this%Temp_24h, this%Press_24h, this%troppresIn,this%qIn, this%rxnRateRatio, &
                this%SpecHum_24h, this%LPAUSE, this%CO_24h, this%CH4_24h,  &
                this%Ravga ,  this%Ravgb, &
                STAT=ios )
   IF ( ios /= 0 ) rc = 1

   CALL finalizeOHparam(this%OHparam)

 end subroutine OH_GridCompFinalize
!EOC
!--------------------------------------------------------------------------
!BOP
      INTEGER FUNCTION JulianDay(nymd)
         INTEGER :: nymd
!
! !DESCRIPTION:
! Determines the Julian day: number between 1 and 365 (or 366).
         INTEGER :: ny, mm, dd
         INTEGER :: m, ds
         INTEGER :: days(12)
!EOP
!--------------------------------------------------------------------------
!BOC
         data days /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/

         ny = nymd / 10000
         mm = mod(nymd, 10000) / 100
         dd = mod(nymd,   100)

         ds = dd

         if ( mm .ne. 1) then
            do m=1, mm-1
               if ( m.eq.2  .and. leap_year(ny) ) then
                  ds = ds + 29
               else
                  ds = ds + days(m)
               endif
            enddo
         endif

         JulianDay = ds
         
      END FUNCTION JulianDay
!EOC
!--------------------------------------------------------------------------
!BOP
      function leap_year(ny)
!
! Determine if year ny is a leap year
!
! Author: S.-J. Lin
      implicit none
      logical leap_year
      integer ny
      integer ny00
!EOP
!--------------------------------------------------------------------------
!BOC

!
! No leap years prior to 0000
!
      parameter ( ny00 = 0000 )   ! The threshold for starting leap-year

      if( ny >= ny00 ) then
         if( mod(ny,100) == 0. .and. mod(ny,400) == 0. ) then
             leap_year = .true.
         elseif( mod(ny,4) == 0. .and. mod(ny,100) /= 0.  ) then
             leap_year = .true.
         else
             leap_year = .false.
         endif
      else
          leap_year = .false.
      endif

      return
    end function leap_year

!EOC
!--------------------------------------------------------------------------
!BOP
      FUNCTION convertTimeToSeconds (nhms) result(this_)
!
! !INPUT PARAMETERS
      INTEGER, intent(in) :: nhms
!
! !RETURNED VALUE
      INTEGER :: this_
!
! !DESCRIPTION:
! Converts the current time into seconds.
!
! !DEFINED PARAMETERS:
      REAL*8, PARAMETER :: SECPHR = 3600.0d0
      REAL*8, PARAMETER :: SECPMN =   60.0d0
      REAL*8, PARAMETER :: SECPDY = 24.0d0*SECPHR
!EOP
!------------------------------------------------------------------------------
!BOC
      this_ = (nhms / 10000) * SECPHR +              &
              (Mod (nhms, 10000) / 100) * SECPMN +   &
               Mod (nhms, 100)

      RETURN

      END FUNCTION convertTimeToSeconds
!EOC
!------------------------------------------------------------------------------
END MODULE  OH_GridCompMod
