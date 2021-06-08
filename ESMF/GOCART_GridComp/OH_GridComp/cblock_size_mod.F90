   MODULE  cblock_size_mod

!
! !DESCRIPTION:
!
!  This module replaces CMN_SIZE, storing common blocks and parameters
!
! !REVISION HISTORY:
!
!  13October2015 Manyin  First crack.
!

   IMPLICIT NONE

      !=======================================================================
      ! CMN_SIZE -- size parameters for GEOS-CHEM arrays (bmy, 10/2/00)
      !=======================================================================

      ! C Preprocessor #define statements for conditional compilation
!#     include "define.h"

      !=======================================================================
      ! IM   - Longitudinal Dimension I=1,IM
      ! JM   - Latitudinal  Dimension J=1,JM
      ! I0   - Zero-PT start of Longitude grid (W. R. T. absolute global grid)
      ! J0   - Zero-PT start of Latitude  grid (W. R. T. absolute global grid)
      ! IMX  - Maximum Longitude Dimension
      ! JMX  - Maximum Latitude  Dimension
      !=======================================================================
!KNJR      INTEGER ::      IM, JM, I0, J0, IMX, JMX
!KNJR      COMMON /CMNSZ1/ IM, JM, I0, J0, IMX, JMX  ! needed for CTM

      !=======================================================================
      ! DISIZE = size (in degrees) of a longitude grid box
      ! DJSIZE = size (in degrees) of a latitude  grid box
      !=======================================================================
#if   defined( GRID4x5  ) 
      REAL*8, PARAMETER :: DISIZE = 5.0d0
      REAL*8, PARAMETER :: DJSIZE = 4.0d0

#elif defined( GRID2x25 )
      REAL*8, PARAMETER :: DISIZE = 2.5d0 
      REAL*8, PARAMETER :: DJSIZE = 2.0d0

#elif defined( GRID1x1 )
      REAL*8, PARAMETER :: DISIZE = 1.0d0 
      REAL*8, PARAMETER :: DJSIZE = 1.0d0

#endif

      !=================================================================
      ! GRID PARAMETERS
      !
      ! IGLOB  = global longitude dimension
      ! JGLOB  = global latitude dimension
      ! LGLOB  = max number of sigma levels 
      ! IIPAR  = window longitude dimension
      ! JJPAR  = window latitude dimension
      ! LLPAR  = window vertical dimension
      ! LLTROP = number of tropospheric levels 
      ! PTOP   = model top pressure (mb)
      !
      ! NOTES: 
      ! (1) GEOS-CHEM is usually set up for a global run, thus,
      !     IIPAR=IGLOB, JJPAR=JGLOB, LLPAR=LGLOB (bmy, 4/9/99)
      !
      ! (2) IIPAR, JJPAR, LLPAR may be smaller than or equal to 
      !     IGLOB, JGLOB, LGLOB (bmy, 4/12/99)    
      !
      ! (3) PTOP is now correct (0.01 hPa) for GEOS-2, GEOS-3 grids
      !     (bmy, 10/2/00)
      !=================================================================
!KNJR      INTEGER, PARAMETER :: IGLOB  = 1
!KNJR      INTEGER, PARAMETER :: JGLOB  = 1
!KNJR      INTEGER, PARAMETER :: LGLOB  = 2

!KNJR      INTEGER, PARAMETER :: IIPAR  = IGLOB
!KNJR      INTEGER, PARAMETER :: JJPAR  = JGLOB
!KNJR      INTEGER, PARAMETER :: LLPAR  = LGLOB

      INTEGER, PARAMETER :: LLTROP = 1

      REAL*8,  PARAMETER :: PTOP   = 10d0




      ! IGCMPAR, JGCMPAR, LGCMPAR - synonyms for IGLOB, JGLOB, LGLOB
      ! These are needed for backwards compatibility w/ CTM code (bmy, 4/9/99)
!KNJR      INTEGER, PARAMETER :: IGCMPAR = IGLOB
!KNJR      INTEGER, PARAMETER :: JGCMPAR = JGLOB
!KNJR      INTEGER, PARAMETER :: LGCMPAR = LGLOB

      !=================================================================
      ! TRACER & EMISSION SPECIES PARAMETERS
      !
      ! NNPAR   = max number of tracers
      ! NEMPARA = max number of anthropogenic emission species
      ! NEMPARB = max number of biogenic      emission species
      !
      ! NOTE: Need to define these for LGEOSCO, even if FULLCHEM
      !       and SMALLCHEM switches are turned off (bmy, 10/2/00)
      !=================================================================
#if   defined( SMALLCHEM )
!KNJR      INTEGER, PARAMETER :: NNPAR   = 6
      INTEGER, PARAMETER :: NEMPARA = 10
      INTEGER, PARAMETER :: NEMPARB = 1

#elif defined( FULLCHEM  )
!KNJR      INTEGER, PARAMETER :: NNPAR   = 24
      INTEGER, PARAMETER :: NEMPARA = 10
      INTEGER, PARAMETER :: NEMPARB = 1

#elif defined( LGEOSCO ) && !defined( FULLCHEM ) && !defined( SMALLCHEM )
!KNJR      INTEGER, PARAMETER :: NNPAR   = 10
      INTEGER, PARAMETER :: NEMPARA = 10
      INTEGER, PARAMETER :: NEMPARB = 1
      
#endif

      INTEGER, PARAMETER :: NNPAR   = 1

      !=================================================================
      ! OTHER PARAMETERS 
      !=================================================================

      ! NVEGTYPE - Maximum number of surface types: 74 olson
      ! NTYPE    - Maximum number of veg types in a 4x5 box
      ! NPOLY    - Number of coefficients for polynomial fits
      INTEGER, PARAMETER :: NVEGTYPE = 74
      INTEGER, PARAMETER :: NTYPE    = 15
      INTEGER, PARAMETER :: NPOLY    = 20

      ! NAIR     - Maximum number (km) for aircraft NOx emissions   
      ! LAIREMS  - Maximum number of layers for aircraft emissions AIREMIS
      INTEGER, PARAMETER :: NAIR    = 20
      INTEGER, PARAMETER :: LAIREMS = LLTROP

      ! NNSTA = max number of time series stations (in inptr.ctm)
      INTEGER, PARAMETER :: NNSTA = 800

      ! MAXIJ - Maximum number of 1st level grid boxes
!KNJR      INTEGER, PARAMETER :: MAXIJ = IIPAR * JJPAR

      ! MAXDEP - Maximum number of depositing species for drydep
      INTEGER, PARAMETER :: MAXDEP = 25

      ! Now define NUMDEP in DRYDEP_SETUP.F, but put it into
      ! a common block here (bmy, 4/9/99)
      INTEGER            :: NUMDEP
!MEM  COMMON /BMYDEP1/      NUMDEP

      ! NBIOMAX - Max number of species emitted in biomass burning
      INTEGER, PARAMETER :: NBIOMAX = 9  

      ! LLCONVM - Max number of layers for convection
!KNJR      INTEGER, PARAMETER :: LLCONVM = LGLOB-1       

      ! NOXLEVELS = Number of levels of anthro NOx emission 
      !             (e.g. surface and 100m)
      ! NOXEXTENT = Highest sigma level that receives anthro NOx emission 
      INTEGER, PARAMETER :: NOXLEVELS = 2
      INTEGER, PARAMETER :: NOXEXTENT = 2 

      ! MAXFAM -- Max number of families for prod and loss output
      INTEGER, PARAMETER :: MAXFAM = 8



!-------------------------------------------------------------------------
CONTAINS
!------------------------------------------------------------------------------
END MODULE  cblock_size_mod
