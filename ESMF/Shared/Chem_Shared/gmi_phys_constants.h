
!=============================================================================
!
! $Id$
!
! CODE DEVELOPER
!   John Tannahill, LLNL
!   jrt@llnl.gov
!
! FILE
!   gmi_phys_constants.h
!
! DESCRIPTION
!   This include file contains some physical constants.
!
!   Mar 30, 2017: Moved this file from GmiInclude/ to Chem_Shared/ for TR
!
!=============================================================================


!     -------------------
!     Physical constants.
!     -------------------

      real*8,  parameter ::  &
     &  ABS_ZERO    = -273.15d0,           & ! absolute zero (degC)
     &  AVG_SRFPRS  = 1000.0d0,            & ! average surface pressure (mb)
     &  AVOGAD      =    6.0221367d+23,    & ! Avogadro number (mole^-1)
     &  BOLTZMN_E   =    1.380662d-16,     & ! Boltzman constant
!    &                                     & ! (erg  *degK^-1*mole^-1)
     &  BOLTZMN_J   =    1.380662d-23,     & ! Boltzman constant
!    &                                     & ! (joule*degK^-1*mole^-1)
!c   &  GAS_CONST_E =    8.314409d+07,   ! gas constant (erg  *degK^-1*mole^-1)
     &  GAS_CONST_J =    8.314409d0,       & ! gas constant (joule*degK^-1*mole^-1)
     &  GMI_G       =    9.81d0,           & ! mean surface gravity accel.  (m/s^2)
     &  GMI_LN2     =    0.693147180559945d0,    & ! natural logarithm of 2.0
     &  GMI_PI      =    3.141592653589793d0,    & ! pi
     &  MASSATM     =    5.14d18,          & ! total mass of the atmosphere (kg)
     &  RADEAR      =    6.371d+06,        & ! radius of the earth (m)
     &  MXRN2       =    0.781d0,          & ! N2 atmospheric volume mixing ratio
     &  MXRO2       =    0.209d0,          & ! O2 atmospheric volume mixing ratio
     &  MXRH2       =    0.5d-06             ! H2 atmospheric volume mixing ratio


!     ---------------------------------------------
!     SAREA_EARTH : surface area of the earth (m^2)
!     ---------------------------------------------

      real*8,  parameter ::  &
     &  SAREA_EARTH = 4.0d0 * GMI_PI * RADEAR * RADEAR


!     --------------------------
!     Molecular weights (g/mol).
!     --------------------------

      real*8,  parameter ::  &
     &  MWTAIR = 28.96d0,  &
     &  MWTH2O = 18.02d0


!     ----------------
!     Other constants.
!     ----------------

      real*8,  parameter ::  &
     &  CGS2MB  =      0.001d0,        & ! cgs units to millibars (cm^2/dyne)
     &  MB2CGS  =   1000.0d0         ! millibars to cgs units (dyne/cm^2)

      real*8,  parameter ::  &
     &  DEGPRAD = 180.0d0 / GMI_PI,    & ! degrees     per radian
     &  RADPDEG =  GMI_PI / 180.0d0  ! radians     per degree

      real*8,  parameter ::  &
     &  KMPCM   =      0.00001d0,      & ! kilometers  per centimeter
     &  KMPM    =      0.001d0,        & ! kilometers  per meter
     &  MPCM    =      0.01d0,         & ! meters      per centimeter
     &  CMPM    =    100.0d0,          & ! centimeters per meter
     &  MPKM    =   1000.0d0,          & ! meters      per kilometer
     &  CMPKM   = 100000.0d0         ! centimeters per kilometer

      real*8,  parameter ::  &
     &  TGPKG   =      1.0d-09,        & ! teragrams   per kilogram
     &  KGPG    =      0.001d0,        & ! kilograms   per gram
     &  GPKG    =   1000.0d0,          & ! grams       per kilogram
     &  KGPTG   =      1.0d+09       ! kilograms   per teragram

      real*8,  parameter ::  &
     &  BPMB    =      0.001d0,        & ! bars        per millibar
     &  MBPPAS  =      0.01d0,         & ! millibars   per pascal
     &  PASPMB  =    100.0d0,          & ! pascals     per millibar
     &  MBPB    =   1000.0d0         ! millibars   per bar

      real*8,  parameter ::  &
     &  PPT_FAC =      1.0d-12,        & ! parts       per trillion factor
     &  PPB_FAC =      1.0d-09,        & ! parts       per billion  factor
     &  PPM_FAC =      1.0d-06       ! parts       per million  factor

