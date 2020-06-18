
!=============================================================================
!
! $Id$
!
! CODE DEVELOPER
!   Dan Bergmann, LLNL
!   dbergmann@llnl.gov
!
! FILE
!   gmi_drydep_data.h
!
! DESCRIPTION
!   This include file sets the data required for the dry deposition algorithm.
!   It includes land types, resistances, and polynomial coefficients from
!   files "drydep_table" and "drydep_coef", which were supplied by Daniel
!   Jacob at Harvard.  The original files follow these parameter statements.
!
!   NOTE THAT SURFACE DEPENDENT PARAMETERS WERE ADDED ON 2/11/02 TO CALCULATE
!   SURFACE RESISTANCE FOR AEROSOLS.
!
!   Mar 30, 2017: Moved this file from GmiInclude/ to Chem_Shared/ for TR
!
!=============================================================================


!     -------------------
!     Integer parameters.
!     -------------------

!     --------------------------------------------------------------------
!     IDEP   : deposition surface type for each Olson surface type
!     IRAC   : resistance that depends on canopy height and density (s^-1)
!     IRCLO  : resistance for leaves, twig, bark in lower canopy    (s^-1)
!     IRCLS  : resistance for leaves, twig, bark in lower canopy    (s^-1)
!     IRGSO  : ground    resistance (s^-1)
!     IRGSS  : ground    resistance (s^-1)
!     IRI    : internal  resistance (s^-1)
!     IRLU   : cuticular resistance (s^-1)
!     IWATER : id's for surface types that are water
!     IZO    : roughness height (m/10000)
!     NWATER : number of Olson's surface types that are water
!     --------------------------------------------------------------------

      integer, parameter :: IDEP(NVEGTYPE) = (/  &
     &  11, 10,  5,  1,  1,  1,  2,  1,  8,  1,  &
     &   1,  1,  1,  1,  1,  1,  5,  1,  1,  1,  &
     &   3,  3,  3,  3,  2,  2,  2,  3,  2,  2,  &
     &   4,  4,  2,  6,  1,  1,  9,  4,  4,  4,  &
     &   5,  5,  5,  5,  5,  9,  5,  5,  5,  5,  &
     &   8,  8,  5,  7,  6,  2,  2,  2,  2,  2,  &
     &   3,  3,  3,  5,  5, 11, 11, 11, 11,  8,  &
     &   1,  8,  9, 11 /)

      integer, parameter :: IRAC(NTYPE) = (/  &
     &     0, 2000, 2000,  200,  100,  &
     &  2000,    0,    0,  300,  100,  &
     &     0,    0,    0,    0,    0 /)

      integer, parameter :: IRCLO(NTYPE) = (/  &
     &  1000, 1000, 1000, 1000, 1000,  &
     &  9999, 9999, 9999, 1000, 9999,  &
     &  9999,    0,    0,    0,    0 /)

      integer, parameter :: IRCLS(NTYPE) = (/  &
     &  9999, 2000, 2000, 2000, 2000,  &
     &  9999, 9999, 9999, 2500, 9999,  &
     &  9999,    0,    0,    0,    0 /)

      integer, parameter :: IRGSO(NTYPE) = (/  &
     &  3500,  200,  200,  150,  200,  &
     &   200,  340,  400, 1000,  300,  &
     &  2000,    0,    0,    0,    0 /)

      integer, parameter :: IRGSS(NTYPE) = (/  &
     &   100,  500,  500,  150,  350,  &
     &   200,  340, 1000,    0,  400,  &
     &     0,    0,    0,    0,    0 /)

      integer, parameter :: IRI(NTYPE) = (/  &
     &  9999,  200,  400,  200,  200,  &
     &   200,  200, 9999,  200, 9999,  &
     &  9999,    0,    0,    0,    0 /)

      integer, parameter :: IRLU(NTYPE) = (/  &
     &  9999, 9000, 9000, 9000, 9000,  &
     &  1000, 4000, 9999, 9000, 9999,  &
     &  9999,    0,    0,    0,    0 /)

      integer, parameter :: IWATER(NVEGTYPE) = (/  &
     &  1, 66, 67, 68, 69, 74,  0,  0,  0,  0,  &
     &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
     &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
     &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
     &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
     &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
     &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
     &  0,  0,  0,  0 /)

      integer, parameter :: IZO(NVEGTYPE) = (/  &
     &    10, 25000,   100,  1000,  1000,  1000, 10000,  1000,   10,  1000,  &
     &  1000,  1000,  1000,  1000,  1000,  1000,  1000,     1, 1000,  1000,  &
     & 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 1000, 10000,  &
     &  1000,  1000,  2000, 10000,  1000,  1000,   100,  1000, 1000,  1000,  &
     &   100,   100,   100,   100,   100,   100,  1000,  1000, 1000,  1000,  &
     &    10,    10,   100,    50, 10000,  2000,  2000,  2000, 2000,  2000,  &
     & 10000, 10000, 10000,  2000,    50,   100,   100,   100,  100,    10,  &
     &     1,     1,   500,    10 /)

      integer, parameter :: NWATER = 6


!     --------------------------------------------------
!     From Zhang et al. [Atmos. Env., p549-560, 2001] =>
!
!     IWET : flag for wet surface = 1, otherwise = -1
!     --------------------------------------------------

      integer, parameter :: IWET(NTYPE) = (/  &
     &     1,   -1,   -1,   -1,   -1,  &
     &    -1,   -1,   -1,    1,   -1,  &
     &     1,    0,    0,    0,    0 /)


!     ----------------
!     Real parameters.
!     ----------------

!     ------------------------------------------
!     DRYCOEFF : polynomial fitting coefficients
!     ------------------------------------------

      real*8, parameter  :: DRYCOEFF(NPOLY) = (/  &
     &  -3.58d-01,  3.02d+00,  3.85d+00, -9.78d-02, -3.66d+00,  &
     &   1.20d+01,  2.52d-01, -7.80d+00,  2.26d-01,  2.74d-01,  &
     &   1.14d+00, -2.19d+00,  2.61d-01, -4.62d+00,  6.85d-01,  &
     &  -2.54d-01,  4.37d+00, -2.66d-01, -1.59d-01, -2.06d-01 /)


!     ------------------------------------------------------------------------
!     From Zhang et al. [Atmos. Env., p549-560, 2001] =>
!
!     A_VEG   : characteristic radius of collectors over vegetated surface (m)
!     ALFA_IM : value for calculating collection efficiency from impaction (m)
!     GAMA_BR : value for calculating collection efficiency from Brownian  (m)
!     ------------------------------------------------------------------------

      real*8, parameter  :: A_VEG(NTYPE) = (/  &
     &  -1.00d+00,  5.10d-03,  3.50d-03,  3.20d-03, 10.00d-03,  &
     &   5.00d-03, -1.00d+00, -1.00d+00, 10.00d-03, 10.00d-03,  &
     &  -1.00d+00,  0.00d+00,  0.00d+00,  0.00d+00,  0.00d+00 /)

      real*8, parameter  :: ALFA_IM(NTYPE) = (/  &
     &  50.00d+00,  0.95d+00,  0.80d+00,  1.20d+00,  1.30d+00,  &
     &   0.80d+00, 50.00d+00, 50.00d+00,  2.00d+00,  1.50d+00,  &
     & 100.00d+00,  0.00d+00,  0.00d+00,  0.00d+00,  0.00d+00 /)

      real*8, parameter  :: GAMA_BR(NTYPE) = (/  &
     &   0.54d+00,  0.56d+00,  0.57d+00,  0.54d+00,  0.54d+00,  &
     &   0.56d+00,  0.54d+00,  0.54d+00,  0.54d+00,  0.56d+00,  &
     &   0.50d+00,  0.00d+00,  0.00d+00,  0.00d+00,  0.00d+00 /)


! ****************Daniel J. Jacob 6/27/94 ****************************************
! *** Olson surface types, corresponding surface types for deposition,
! ***  and corresponding roughness height
! ********************************************************************************
! Olson ID# - deposition ID# - z0 in 1.E-4 m
!      1    11     1    Water
!      2    10 10000    Urban
!      3     5    50    Shrub
!      4     1  1000    Not used
!      5     1  1000    ibid.
!      6     1  1000    ibid.
!      7     2 10000    Tropical evergreen
!      8     1  1000    Not used
!      9     8     1    Desert
!     10     1  1000    Not used
!     11     1  1000    ibid.
!     12     1  1000    ibid.
!     13     1  1000    ibid.
!     14     1  1000    ibid.
!     15     1  1000    ibid.
!     16     1  1000    ibid.
!     17     5  1000    Scrub
!     18     1     1    Ice
!     19     1  1000    Not used
!     20     1  1000    ibid.
!     21     3 10000    Conifer
!     22     3 10000    Conifer
!     23     3 10000    Conifer
!     24     3 10000    Conifer/deciduous
!     25     2 10000    Deciduous/conifer
!     26     2 10000    Deciduous
!     27     2 10000    Deciduous
!     28     3 10000    Conifer
!     29     2  1000    Dwarf forest
!     30     2 10000    Tropical broadleaf
!     31     4  1000    Agric.
!     32     4  1000    ibid.
!     33     2  1000    Dec. woodland
!     34     6 10000    Tropical rainforest
!     35     1  1000    Not used
!     36     1  1000    Not used
!     37     9   100    Rice paddies
!     38     4  1000    Agric.
!     39     4  1000    Agric.
!     40     4  1000    Agric.
!     41     5   100    Shrub/grass
!     42     5   100    Shrub/grass
!     43     5   100    Shrub/grass
!     44     5   100    Shrub/grass
!     45     5   100    Shrub/grass
!     46     9   100    Wetland
!     47     5  1000    Scrub
!     48     5  1000    Scrub
!     49     5  1000    Scrub
!     50     5  1000    Scrub
!     51     8   100    Desert
!     52     8   100    Desert
!     53     5   100    Steppe
!     54     7    50    Tundra
!     55     6 10000    Rainforest
!     56     2  1000    Mixed wood/open
!     57     2  1000    Mixed wood/open
!     58     2  1000    Mixed wood/open
!     59     2  1000    Mixed wood/open
!     60     2  1000    Mixed wood/open
!     61     3 10000    Conifers
!     62     3 10000    Conifers
!     63     3 10000    Conifers
!     64     5  1000    Wooded tundra
!     65     5    50    Moor
!     66    11     1    Coastal
!     67    11     1    Coastal
!     68    11     1    Coastal
!     69    11     1    Coastal
!     70     8    10    Desert
!     71     1    10    Ice
!     72     8     1    Salt flats
!     73     9   500    Wetland
!     74    11     1    Water
! **
!   6  1 66 67 68 69 74                    #/Olson ID's of water surface types
! **  Resistances (s m-1) for each deposition surface type,
! **  and maximum deposition velocity Vsmax(1.d-2 cm s-1) for aerosol
!   Type  Ri  Rlu  Rac RgsS RgsO RclS RclO Vsmax
!     1 9999 9999    0  100 3500 9999 1000  100   Snow/Ice (Wesely, AE 1989)-
!                                                   listed first.
!     2  200 9000 2000  500  200 2000 1000  100   Deciduous  forest (Wesely)
!     3  400 9000 2000  500  200 2000 1000  100   Coniferous forest (Wesely)
!     4  200 9000  200  150  150 2000 1000  100   Agricultural land (Wesely)
!     5  200 9000  100  350  200 2000 1000  100   Shrub/grassland   (Wesely)
!     6  200 1000 2000  200  200 9999 9999  100   Amazon forest
!                                                   (Jacob & Wofsy, JGR 1990)
!     7  200 4000    0  340  340 9999 9999  100   Tundra
!                                                   (Jacob et al., JGR 1992)
!     8 9999 9999    0 1000  400 9999 9999   10   Desert  (Wesely)
!     9  200 9000  300    0 1000 2500 1000  100   Wetland (Wesely)
!    10 9999 9999  100  400  300 9999 9999  100   Urban   (Wesely)
!    11 9999 9999    0    0 2000 9999 9999   10   Water   (Wesely)
!

