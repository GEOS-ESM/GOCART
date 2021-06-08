   MODULE  cblock_CO_mod

!
! !DESCRIPTION:
!
!  This module replaces CMN_CO, storing common blocks and parameters
!
! !REVISION HISTORY:
!
!  13October2015 Manyin  First crack.
!

   IMPLICIT NONE

!
!  Created by bnd.
!
!MEM  INTEGER NSEASON_last2, &
!MEM &        NOH,NNCOUNT,NNEW,MNCOUNT,MNEW,&
!MEM &        NTALDT,LNCOUNT,LNEW,NCLIMATOLOGY,NCLIMATOLOGY2, &
!MEM &        NCLIMATOLOGY3,LMN_last

      INTEGER, PARAMETER :: NFIELDS2=12
      INTEGER, PARAMETER :: NJVREAD=65
      INTEGER, PARAMETER :: NO3READ=65

!KNJR      INTEGER Rcount(IIPAR,JJPAR,LLPAR)
!KNJR     REAL*8 BBIJ(IIPAR,JJPAR,LLPAR,NFIELDS2),
!KNJR    &       o3up(IIPAR,JJPAR,LLPAR,12),
!KNJR    &       toms_climO3(IIPAR,JJPAR,12)
!KNJR      REAL*8 KRATE(IIPAR,JJPAR,LLPAR),
!KNJR     &       SUMACETCO(IIPAR,JJPAR),
!KNJR     &       SUMISOPCO(IIPAR,JJPAR),SUMMONOCO(IIPAR,JJPAR),
!KNJR     &       EMACET(IIPAR,JJPAR,12)
!KNJR      REAL*8 GCO(IIPAR,JJPAR,LLPAR),GO3(IIPAR,JJPAR,LLPAR)
!KNJR      REAL*8 GISOP(IIPAR,JJPAR,LLPAR)
!KNJR      REAL*8 Tavg(IIPAR,JJPAR,LLPAR),Pavg(IIPAR,JJPAR,LLPAR),
!KNJR      REAL*8      Ravga(IIPAR,JJPAR,LLPAR),Ravgb(IIPAR,JJPAR,LLPAR)
!KNJR      REAL*8     Wavg(IIPAR,JJPAR,LLPAR)
!KNJR       REAL*8     Ravg(IIPAR,JJPAR,LLPAR)
!KNJR      REAL*8 STTTOGCO(IIPAR,JJPAR,LLPAR)
!KNJR     REAL*8 STTTOPPB(IIPAR,JJPAR,LLPAR)
!KNJR     REAL*8 bairdens(IIPAR,JJPAR,LLPAR),AVGNO
          REAL*8, PARAMETER :: AVGNO=6.022E23
!KNJR     COMMON /COb1/ o3up,BBIJ,bairdens,toms_climO3
!KNJR     COMMON /COb2/ NNCOUNT,NNEW,MNCOUNT,MNEW,NSEASON_last2,
!KNJR    &              LMN_last,NTALDT
!KNJR     COMMON /COb3/ Rcount
!KNJR     COMMON /COb4/ Ravg,Ravga,Ravgb,Wavg
!KNJR     COMMON /COb5/ SUMISOPCO,SUMMONOCO,SUMACETCO,
!KNJR    &              EMACET
!KNJR     COMMON /COb7/ LNCOUNT,LNEW

!KNJR      REAL*8 CO_prod(JJPAR,LLPAR),CO_loss(JJPAR,LLPAR)
!KNJR      COMMON /COb8/ CO_prod,CO_loss

      INTEGER NREFL,NREFLCOL
      PARAMETER (NREFL=226,NREFLCOL=5)
      REAL*8 OPTDEPTH(NREFL),RFLC(NREFL,NREFLCOL)
!MEM  COMMON /COb9/OPTDEPTH,RFLC

!      REAL*8 FASTJO3(31,18,12),FASTJT(41,18,12),
!     *    FASTJZ(IIPAR,JJPAR,LLPAR)
!      COMMON/OTHERCLIM/ FASTJO3,FASTJT,FASTJZ
      
! TOMS O3 column interannual variation
!KNJR     REAL*8 scaleO3(IIPAR,JJPAR)
!KNJR     COMMON/TOMSO3/scaleO3

!***************************************************************
! Parameters & arrays used in SR INTERPOH2 & SR READAVGOH2.
!***************************************************************
!KNJR      REAL*8 avgOH2(IIPAR,JJPAR,LLPAR)
!KNJR      COMMON /CMSOH2/ avgOH2

!KNJR      REAL*8 OH94(IIPAR,JJPAR,LLPAR)
!KNJR      COMMON /CMSOH3/ OH94


!-------------------------------------------------------------------------
CONTAINS
!------------------------------------------------------------------------------
END MODULE  cblock_CO_mod
