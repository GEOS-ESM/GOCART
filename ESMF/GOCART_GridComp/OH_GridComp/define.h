!
!******************************************************************************
!  Include file "define.h" specifies C-preprocessor "switches" that are 
!  used to include or exclude certain sections of code.  
!  (bmy, bdf, 1/30/98, 10/3/00)
!
!  List of "Switches"
!  ===========================================================================
!  (1 ) GEOS_1      -> Enables code for GEOS-1 met fields & chemistry
!  (2 ) GEOS_STRAT  -> Enables code for GEOS-STRAT met fields & chemistry
!  (3 ) GEOS_2      -> Enables code for GEOS-2 met fields & chemistry
!  (4 ) GEOS_3      -> Enables code for GEOS-3 met fields & chemistry
!  (5 ) GRID1x1     -> Enables code for 1 x 1   grid
!  (6 ) GRID2x25    -> Enables code for 2 x 2.5 grid 
!  (7 ) GRID4x5     -> Enables code for 4 x 5   grid 
!  (8 ) FULLCHEM    -> Enables code for "Full" Chemistry (ISOP and NMHC)
!  (9 ) SMALLCHEM   -> Enables code for "Small" chemistry (No ISOP, no NMHC)
!  (10) LGEOSCO     -> Enables code for CO run w/ parameterized OH
!  (11) LFASTJ      -> Enables code for FAST-J photolysis
!  (12) LSLOWJ      -> Enables code for SLOW-J photolysis
!
!  NOTES:
!  (1) "define.h" is #include'd at the top of CMN_SIZE.  All subroutines
!      that normally reference CMN_SIZE will also reference "define.h".
!
!  (2) Only define the "switches" that are *absolutely* needed for a
!      given implementation, as the criteria for code inclusion/exclusion 
!      is the #if defined() statement.  Undefined "switches" are "off".
!
!  (3) To turn off a switch, comment that line of code out.
!
!  (4) As of 11/30/99, DO_MASSFLUX is obsolete, since the mass flux
!      arrays are now declared allocatable in "diag_mod.f".  
!
!  (5) Eliminate DO_MASSB switch -- ND63 diagnostic is now obsolete.
!      (bmy, 4/12/00)
!
!  (6) Add GEOS_3 and GRID1x1 switches for future use (bmy, 7/7/00)
!
!  (7) Make sure that one of FULLCHEM, SMALLCHEM, or LGEOSCO is turned on.
!      Also cosmetic changes. (bmy, 10/3/00)
!
!  Undefine all "switches" so that they cannot be accidentally reset  
!******************************************************************************
!
#undef GEOS_1
#undef GEOS_STRAT
#undef GEOS_2
#undef GEOS_3
#undef GRID2x25  
#undef GRID4x5
#undef GRID1x1   
#undef FULLCHEM  
#undef SMALLCHEM 
#undef LGEOSCO
#undef LFASTJ
#undef LSLOWJ
#undef SGI
#undef COMPAQ
#undef LINUX
#undef SPARC
#undef IBM_AIX
!
!******************************************************************************
!  Define the necessary "switches" for GEOS-CTM
!  Give each switch it's own name as a value, since this will prevent
!  the C-preprocessor overwriting the name everywhere in the code.
!******************************************************************************
!
#define GEOS_1      GEOS_1       
!#define GEOS_STRAT  GEOS_STRAT
!#define GEOS_2      GEOS_2
!#define GEOS_3      GEOS_3

!#define GRID1x1     GRID1x1
#define GRID2x25    GRID2x25
!#define GRID4x5     GRID4x5

!#define SMALLCHEM   SMALLCHEM
!#define FULLCHEM    FULLCHEM
#define LGEOSCO     LGEOSCO

!#define LFASTJ      LFASTJ
!#define LSLOWJ      LSLOWJ

!#define SGI         'SGI'
!#define COMPAQ      'COMPAQ'
#define LINUX       'LINUX'
!#define SPARC       'SPARC'
!#define IBM_AIX     'IBM_AIX'

!
!****************************************************************************
!  Force a compile-time error if switches GEOS_1, GEOS_STRAT, 
!  and GEOS_2, and GEOS_3 are all undefined. 
!****************************************************************************
!
#if !defined( GEOS_1 ) && !defined( GEOS_STRAT ) && !defined( GEOS_2 ) && !defined( GEOS_3 )
#error "ERROR: GEOS_1, GEOS_STRAT, GEOS_2, and GEOS_3"
#error "are ALL undefined in header file define.h"
#endif
!
!****************************************************************************
!  Force a compile-time error if switches GRID1x1, GRID2x25,
!  and GRID4x5 are all undefined. 
!****************************************************************************
!
#if !defined( GRID2x25 ) && !defined( GRID4x5 ) && !defined( GRID1x1 )
#error "ERROR: GRID2x25, GRID4x5, and GRID1x1"
#error "are ALL undefined in header file define.h"
#endif
!
!****************************************************************************
!  Force a compile-time error if switches FULLCHEM, SMALLCHEM,
!  and LGEOSCO are all undefined
!****************************************************************************
!
#if !defined( FULLCHEM ) && !defined( SMALLCHEM ) && !defined( LGEOSCO )
#error "ERROR: One of FULLCHEM, SMALLCHEM, LGEOSCO 
#error "needs to be defined in header file define.h"
#endif
