   MODULE  cblock_OH_mod

!
! !DESCRIPTION:
!
!  This module replaces CMN_OH, storing common blocks and parameters
!
! !REVISION HISTORY:
!
!  13October2015 Manyin  First crack.
!

   IMPLICIT NONE

!
!  Created by Bryan Duncan.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The user needs to supply: 
!
! 1)  Details of the model grid resolution.
!
! Array dimensions of troposphere:
!
!   MVRTBX  = # of boxes in vertical
!   MLATBX  = # of boxes from pole to pole (following one longitude band)
!   MLONBX  = # of boxes circling globe (following one latitude band)
!
!KNJR      INTEGER, PARAMETER :: MVRTBX = 2
!KNJR      INTEGER, PARAMETER :: MLATBX = 1
!KNJR      INTEGER, PARAMETER :: MLONBX = 1
!EY DEBUG
!
! 2)  Filenumbers to be used by the parameterization code for reading
!     input data.  Make sure these numbers do not conflict with numbers
!     within the user's code.
!
      INTEGER, PARAMETER :: NGENERICFILE = 65
!
! 3)  Error Check.  If you want error or warning messages printed
!     set ERRORON = 0, otherwise set it to 1. This flag will only
!     stop some error and warning messages from being printed that 
!     are found in SR CALC_OH.  All other messages, the user will
!     need to comment out manually in the code.
!
      INTEGER, PARAMETER :: ERRORON = 1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The user should not change anything below this point!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!***************************************************************
! Parameters used in SR READCOEFF, SR GETOH, & others.
!***************************************************************
      INTEGER, PARAMETER :: NPARAMA = 27
      INTEGER, PARAMETER :: NPARAMB = 27
      INTEGER, PARAMETER :: NPARAMC = 22
      INTEGER, PARAMETER :: NPARAMD = 22
      INTEGER, PARAMETER :: NPARAME = 13
      INTEGER, PARAMETER :: NPARAMF = 11
!
! NPARAM = # of parameterizations
      INTEGER, PARAMETER :: NPARAM = 3*NPARAMA+3*NPARAMB+3*NPARAMC+ &
                                     NPARAMD+NPARAME+NPARAMF

      INTEGER NDFUNCSA1,NDFUNCSA2,NDFUNCSA3,NDFUNCSB1, &
     &        NDFUNCSC1,NDFUNCSC2,NDFUNCSC3,NDFUNCSD1, &
     &        NDFUNCSE1,NDFUNCSF1,NDFUNCSB2,NDFUNCSB3

!MEM  COMMON /NDFUNCS/NDFUNCSA1,NDFUNCSA2,NDFUNCSA3,NDFUNCSB1, &
!MEM &        NDFUNCSC1,NDFUNCSC2,NDFUNCSC3,NDFUNCSD1,NDFUNCSE1, &
!MEM &        NDFUNCSF1,NDFUNCSB2,NDFUNCSB3

!***************************************************************
! Variables used in SR READ_COEFF.
!***************************************************************
! MXVAR = maximum # of independent variables upon
!    which the parameterized [OH] is dependent.
      INTEGER, PARAMETER :: MXDVD=21
      INTEGER, PARAMETER :: MXELM=(2**((MXDVD/3)+1)-1)
      INTEGER, PARAMETER :: MXVAR=20
      INTEGER, PARAMETER :: MXCOL=300
      INTEGER, PARAMETER :: MXCLM=MXVAR+2
      INTEGER, PARAMETER :: MXROW=MXELM+(MXVAR*2)
      INTEGER, PARAMETER :: MXDONE=45

      INTEGER NENDW(NPARAM,MXDONE),IROWST(NPARAM,MXELM,0:MXROW)
      INTEGER NDONE(NPARAM),IDENTOLD(NPARAM,MXDONE)
      REAL*8 RANGEM(NPARAM,MXVAR,2)
      INTEGER MOVAR(NPARAM)

!MEM  COMMON /STARTUP/ RANGEM
!MEM  COMMON /ISTARTUP/ MOVAR

      REAL*8 ELTODO(NPARAM,MXROW,MXCLM),COEFF(NPARAM,MXDONE,MXCOL)
      REAL*8 ACCUR(NPARAM)    !  MEM: this is only used during read routine

!MEM  COMMON /EQUAL/ ELTODO,COEFF,ACCUR
!MEM  COMMON /INEQUAL/ NENDW,IROWST
!MEM  COMMON /NDONES/ NDONE,IDENTOLD

!***************************************************************
! Variables used in SR GETINFO & SR GETOH.
!***************************************************************
! Input data provided by the user is stored in INDVARA-D.
!
!KNJR      REAL*8 INDVARA(MLONBX,MLATBX,MVRTBX,MXVAR), &
!KNJR     &       INDVARB(MLONBX,MLATBX,MVRTBX,MXVAR), &
!KNJR     &       INDVARC(MLONBX,MLATBX,MVRTBX,MXVAR), &
!KNJR     &       INDVARD(MLONBX,MLATBX,MVRTBX,MXVAR)
      REAL*8 INDVAR(MXVAR),PARAMOH
!KNJR      COMMON /COb6/ INDVARA,INDVARB,INDVARC,INDVARD
!MEM  COMMON /COb8/ INDVAR,PARAMOH
!
!***************************************************************
! The tropospheric domain is divided into subdomains based on
! season, latitude, pressure, NOx & isoprene concentrations.
!***************************************************************
!
! Variables defining parameterizations.
!
! NSEAS   = number of seasons.
! NLATS   = number of latitude bands.
! COCOUNT = selected parameterization out of NPARAM parameterizations. 
!
      INTEGER, PARAMETER :: NSEAS=4
      INTEGER, PARAMETER :: NLATS=7
      INTEGER COCOUNT
!MEM  COMMON /COUNTERSA/COCOUNT
!
! ALATS(NLATS) = northernmost boundary for latitude band.
!
      REAL*8 ALATS(NLATS)
      DATA ALATS /-60.,-40.,-30.,0.,30.,40,60./
!      DATA ALATS /-60.,-40.,-29.,0.,29.,40,60./
!
! PRESSES = pressure levels.
! CNOXS   = NOx levels.
!
      REAL*8 PRESSES(3),CNOXS(7)
      DATA PRESSES/800.,350.,100./
      DATA CNOXS/300.,500.,1000.,5000.,5000.,100.,40./ 
!
! Variables specific to box in question (see SR GETOH):
!
! OH_NOX    = NOx concentration
! OH_PRESS  = pressure
! OH_LAT    = latitude 
! OH_ISOP   = isoprene concentration
! OH_O3     = ozone
! OH_LON    = longitude
! OH_MONTH = month (January-December = 1-12)
! OH_SEASON = (Northern Hemispheric) season
!        1 --> winter (Dec, Jan, Feb)
!        2 --> spring (Mar, Apr, May)
!        3 --> summer (Jun, Jul, Aug)
!        4 --> autumn (Sep, Oct, Nov)
!
      REAL*8 OH_NOX,OH_PRESS,OH_LAT,OH_ISOP,OH_LON,OH_O3
      INTEGER OH_SEASON,OH_MONTH
!MEM  COMMON /OHa/ OH_NOX,OH_PRESS,OH_LAT,OH_ISOP,OH_LON,OH_O3
!MEM  COMMON /iOHb/ OH_SEASON
!
!***************************************************************
! Parameters & arrays used in SR INTERPOH & SR READAVGOH.
!***************************************************************
      INTEGER, PARAMETER :: NCMSALTS=7
      INTEGER, PARAMETER :: NCMSLATS=24
      REAL*8 CMSALTS(NCMSALTS),CMSLATS(NCMSLATS)
      REAL*8 avgOH(NSEAS,NCMSLATS,NCMSALTS)
      DATA CMSALTS/1000.,900.,800.,700.,500.,300.,200./
      DATA CMSLATS/89.,84.,76.,68.,60.,52.,44.,36.,28.,20.,12.,4., &
     & -4.,-12.,-20.,-28.,-36.,-44.,-52.,-60.,-68.,-76.,-84.,-89./
!MEM  COMMON /CMSOH/ avgOH
!
!***************************************************************
! Parameters & arrays used in SR AerosolOH & SR READ_Percentages.
!***************************************************************

      INTEGER, PARAMETER :: numStratLevels = 18
!KNJR      REAL*8 dustrat(MLONBX,MLATBX,numStratLevels,12)
!KNJR      COMMON /CMSAerosol/ dustrat
! Correction factors for tropical UT
      INTEGER, PARAMETER :: NLEV_UT=4

! MEM These are not used
!     REAL*8 CF_UT_NH(NLEV_UT,12)
!     REAL*8 CF_UT_SH(NLEV_UT,12)
!

      REAL*8 valNH(NLEV_UT*12),valSH(NLEV_UT*12)
      DATA valNH / 0.71,0.49,0.33,0.24, &
                   0.66,0.46,0.31,0.23, &
                   0.63,0.42,0.29,0.21, &
                   0.62,0.41,0.27,0.20, &
                   0.67,0.46,0.30,0.21, &
                   0.67,0.45,0.28,0.20, &
                   0.65,0.44,0.29,0.20, &
                   0.72,0.50,0.32,0.22, &
                   0.79,0.53,0.32,0.20, &
                   0.79,0.52,0.30,0.20, &
                   0.81,0.52,0.31,0.21, &
                   0.83,0.57,0.37,0.26 /

      DATA valSH / 0.67,0.43,0.26,0.19, &
                   0.67,0.42,0.26,0.19, &
                   0.83,0.53,0.33,0.24, &
                   0.80,0.52,0.32,0.23, &
                   0.69,0.46,0.30,0.22, &
                   0.62,0.42,0.29,0.21, &
                   0.56,0.39,0.27,0.20, &
                   0.58,0.40,0.27,0.19, &
                   0.58,0.38,0.24,0.16, &
                   0.65,0.41,0.26,0.18, &
                   0.67,0.44,0.28,0.20, &
                   0.72,0.45,0.28,0.21 /


!-------------------------------------------------------------------------
CONTAINS
!------------------------------------------------------------------------------
END MODULE  cblock_OH_mod
