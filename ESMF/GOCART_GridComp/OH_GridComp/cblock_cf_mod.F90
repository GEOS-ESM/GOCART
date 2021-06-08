   MODULE  cblock_cf_mod

!
! !DESCRIPTION:
!
!  This module replaces CMN_CF, storing common blocks and parameters
!
! !REVISION HISTORY:
!
!  13October2015 Manyin  First crack.
!

   IMPLICIT NONE

!
!  CMN_CF Created by Bryan Duncan.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!***************************************************************
! Parameters & arrays used in SR READ_CORRECTION & SR CORRECTOH.
!***************************************************************
      INTEGER, PARAMETER :: NCMSALTS2 = 9
      INTEGER, PARAMETER :: NCMSLATS2 = 24
      INTEGER, PARAMETER :: NSEAS_CF  = 4
      INTEGER, PARAMETER :: NFILENUM  = 65

      REAL*8 CMSALTS2(NCMSALTS2),CMSLATS2(NCMSLATS2)
      REAL*8 correction(NSEAS_CF,NCMSLATS2,NCMSALTS2)

      DATA CMSALTS2/1000.,900.,800.,700.,500.,300.,200.,150.,100./

      DATA CMSLATS2/89.,84.,76.,68.,60.,52.,44.,36.,28.,20.,12.,4., &
     & -4.,-12.,-20.,-28.,-36.,-44.,-52.,-60.,-68.,-76.,-84.,-89./

!MEM  COMMON /CMSCFS/ correction


!-------------------------------------------------------------------------
CONTAINS
!------------------------------------------------------------------------------
END MODULE  cblock_cf_mod
