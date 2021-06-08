#include "MAPL_Generic.h"

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  computeParameters_mod ---
!
! !INTERFACE:
!

   MODULE  computeParameters_mod

! !USES:

   USE ESMF
   USE MAPL

   IMPLICIT NONE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PRIVATE
      PUBLIC  :: PARAMA, PARAMB, PARAMC
      PUBLIC  :: PARAMD, PARAME, PARAMF
                
!**********************************************************************
CONTAINS
!**********************************************************************
      SUBROUTINE PARAMA(BLAT,BSEASON,COCOUNT,NSEAS,NLATS,ALATS)
!**********************************************************************
      IMPLICIT NONE
      INTEGER K,M,J,BSEASON,COCOUNT,NSEAS,NLATS
      REAL*8 BLAT,ALATS(NLATS)
!
! Information used in Error Check.
!
      INTEGER NERROR,CSEASON1(18),CSEASON2(21),CSEASON3(21), &
     &        CSEASON4(21),NSTOP
      REAL*8 DLATS1(6),DLATS2(12),DLATS3(12), &
     & DLATS4(12),DLATS5(12),DLATS6(12),DLATS7(9),DLATS8(6) 
      DATA CSEASON1/1,3,7,11,15,19,28,30,34,38,42,46,55,57,61,65, &
     &              69,73/
      DATA CSEASON2/4,8,12,16,20,23,26,31,35,39,43,47,50,53,58,62, &
     &              66,70,74,77,80/
      DATA CSEASON3/5,9,13,17,21,24,27,32,36,40,44,48,51,54, &
     &              59,63,67,71,75,78,81/
      DATA CSEASON4/2,6,10,14,18,22,25,29,33,37,41,45,49,52,56, &
     &              60,64,48,72,76,79/
      DATA DLATS1/1,2,28,29,55,56/
      DATA DLATS2/3,30,57,4,31,58,5,32,50,6,33,60/
      DATA DLATS3/7,34,61,8,35,62,9,36,63,10,37,64/
      DATA DLATS4/11,38,65,12,39,66,13,40,67,14,41,68/
      DATA DLATS5/15,42,69,16,43,70,17,44,71,18,45,72/
      DATA DLATS6/19,46,73,20,47,74,21,48,75,22,49,76/
      DATA DLATS7/23,50,77,24,51,78,25,52,79/
      DATA DLATS8/27,54,81,26,53,80/
!**********************************************************************
! Created by Bryan Duncan.
!**********************************************************************
! SR PARAMA is called by SR CLEAN.  Based on season and latitude, this
! subroutine choses the appropriate identification number of the 
! parameterization required.  
!
!  ----------------------------------------------------------
!  | Latitude   |   DJF    |   MAM    |   JJA    |   SON    | 
!  ----------------------------------------------------------
!  |   > 80 N   |   AVG    |   AVG    |   AVG    |   AVG    |
!  |  60-80 N   |   AVG    | 26 43 80 | 27 54 81 |   AVG    |
!  |  40-60 N   |   AVG    | 23 50 77 | 24 51 78 | 25 52 79 |  
!  |  30-40 N   | 19 46 73 | 20 47 74 | 21 48 75 | 22 49 76 |
!  |   0-30 N   | 15 42 69 | 16 43 70 | 17 44 71 | 18 45 72 |  
!  |   0-30 S   | 11 38 65 | 12 39 66 | 13 40 67 | 14 41 68 |  
!  |  30-40 S   |  7 34 61 |  8 35 62 |  9 36 63 | 10 37 64 |  
!  |  40-60 S   |  3 30 57 |  4 31 58 |  5 32 59 |  6 33 60 |  
!  |  60-80 S   |  1 28 55 |   AVG    |   AVG    |  2 29 56 |  
!  |   > 80 S   |   AVG    |   AVG    |   AVG    |    AVG   |  
!  ----------------------------------------------------------
!
! The table above illustrates the spatial and temporal coverage of the
! parameterizations found in the "CLEAN" subdomain.  In the boxes labeled 
! "AVG", the OH is sufficiently low so that the climatological mean OH for
! that season and latitude is used (Spivakovsky et al., 1999).  (For these 
! boxes, the parameterized OH is not calculated.)  For all other boxes,
! the parameterized OH is calculated.  The three numbers in the box correspond
! to the identification numbers of the parameterization for that box.
! The first number is for the surface layer, the second for the middle 
! tropospheric layer and the third for the upper tropospheric layer.
! For instance, the parameterization the code would use to predict OH for
! April, 15N, middle troposphere has the identification number 43.
!**********************************************************************
! List of Variables & Arrays
!
! COCOUNT = identification number of parameterization. 
!
!**********************************************************************
! Loop over latitude bands and then seasons to find appropriate
! parameterization.
!
      DO K=1,NLATS
!
!***************
!
        IF(BLAT.GE.ALATS(K).AND.K.NE.NLATS) THEN
          GOTO 298
        ENDIF
!
           NSTOP=0
           NERROR=1
!
        IF(BLAT.LT.ALATS(1)) THEN
          NERROR=0
!	  IF (MAPL_AM_I_ROOT()) print *, 'PARAMA: polar night check SH'
          IF(BSEASON.EQ.2) NSTOP=1
          IF(BSEASON.EQ.3) NSTOP=1
        ENDIF 
        IF(BLAT.LT.ALATS(7).AND.BLAT.GE.ALATS(1)) THEN
         COCOUNT=(K-1)*NSEAS-2
!	 IF (MAPL_AM_I_ROOT()) print *, 'PARAMA: COCOUNT updated1 COCOunt, K ,nseas :', COCOUNT, K, NSEAS
        ENDIF 
        IF(BLAT.GE.ALATS(7)) THEN
!	  IF (MAPL_AM_I_ROOT()) print *, 'PARAMA: polar night check NH'
          IF(BSEASON.EQ.1) NSTOP=1
          IF(BSEASON.EQ.4) NSTOP=1
         COCOUNT=(K)*NSEAS-2-1
!	 IF (MAPL_AM_I_ROOT()) print *, 'PARAMA: COCOUNT updated2 COCOUNT , K, NSEAS:', COCOUNT
        ENDIF 
!
!**********************************************************************
!
! Error Check.
!
        IF(NERROR.EQ.1.AND.COCOUNT.EQ.0) THEN
          NSTOP=1
!	  IF (MAPL_AM_I_ROOT()) print *, 'PARAMA: NERROR = 1, COCOUNT = 0'
        ENDIF
!
! NSTOP is an error check to see if a parameterization which
!       does not exist has been called.  This error indicates
!       that there is a problem in SR SKIPS.

        IF(NSTOP.EQ.1) THEN
          PRINT*,'***********************************'
          PRINT*,'Stopped in SR PARAMA/SR CLEAN!'
          PRINT*,'A parameterization which does not exist has been'
          PRINT*,'called.'
          PRINT*,'Season=',BSEASON
          PRINT*,'Latitude=',BLAT
          PRINT*,'***********************************'
          STOP
        ENDIF
!
! End Error Check.
!
!**********************************************************************
! Loop over seasons.
!
        DO M=1,NSEAS
!
!***********
!
!	IF (MAPL_AM_I_ROOT()) PRINT*,'PARAMA: Before 299:  ' , COCOUNT
        IF(BLAT.LT.ALATS(7).AND.BLAT.GE.ALATS(6).AND.M.EQ.1) GOTO 299 
        IF(BLAT.LT.ALATS(1).AND.(M.EQ.2.OR.M.EQ.3)) GOTO 299 
        IF(BLAT.GE.ALATS(7).AND.M.EQ.1) GOTO 299
        IF(BLAT.GE.ALATS(7).AND.M.EQ.4) GOTO 299

          IF(BSEASON.GT.M) THEN
            COCOUNT=COCOUNT+1
!	    IF (MAPL_AM_I_ROOT()) print *, 'PARAMA: COCOUNT updated :', COCOUNT
            GOTO 299
          ENDIF
!
                 COCOUNT=COCOUNT+1
!	IF (MAPL_AM_I_ROOT()) PRINT*,'PARAMA: Updated COCOUNT   ' , COCOUNT
		 
!
                 GOTO 300
!
 299        CONTINUE
!***********
        ENDDO ! End Loop over Seasons.
!***********
 298        CONTINUE
!***************
      ENDDO   ! End Loop over Latitude Bands.
!***************
!
 300  CONTINUE
!
!**********************************************************************
!
! Error Check.
!
!	IF (MAPL_AM_I_ROOT()) PRINT*,'PARAMA: Before 301 GOTO: ' , COCOUNT
      IF(BSEASON.EQ.1) THEN
        DO J=1,18
         IF(COCOUNT.EQ.CSEASON1(J)) GOTO 301
        ENDDO
 !       IF (MAPL_AM_I_ROOT()) PRINT*,'Season is DJF:  Did not call a DJF parameterization!'
      ENDIF
      IF(BSEASON.EQ.2) THEN
        DO J=1,21
         IF(COCOUNT.EQ.CSEASON2(J)) GOTO 301
        ENDDO
  !      IF (MAPL_AM_I_ROOT()) PRINT*,'Season is MAM:  Did not call a MAM parameterization!'
      ENDIF
      IF(BSEASON.EQ.3) THEN
        DO J=1,21
         IF(COCOUNT.EQ.CSEASON3(J)) GOTO 301
        ENDDO
   !     IF (MAPL_AM_I_ROOT()) PRINT*,'Season is JJA:  Did not call a JJA parameterization!'
      ENDIF
      IF(BSEASON.EQ.4) THEN
        DO J=1,21
         IF(COCOUNT.EQ.CSEASON4(J)) GOTO 301
        ENDDO
    !    IF (MAPL_AM_I_ROOT()) PRINT*,'Season is SON:  Did not call a SON parameterization!'
      ENDIF

        PRINT*,'Stopped in SR PARAMA.'
        PRINT*,'COCOUNT = ',COCOUNT
        STOP
!
 301  CONTINUE
!
! EY DEBUG31 this where we exit the loop in the first 302 GOTO statement
!	IF (MAPL_AM_I_ROOT()) print *, 'PARAMA: BEFORE GOTO 302 statements ', COCOUNT
      IF(BLAT.LT.ALATS(1)) THEN
        DO J=1,9
         IF(COCOUNT.EQ.DLATS1(J)) GOTO 302 
        ENDDO
 !     IF (MAPL_AM_I_ROOT()) PRINT*,'Latitude is not lt 60S: Did not call 60-80S param!'
      ENDIF

      IF(BLAT.LT.ALATS(2).AND.BLAT.GE.ALATS(1)) THEN
        DO J=1,12
         IF(COCOUNT.EQ.DLATS2(J)) GOTO 302 
        ENDDO
  !    IF (MAPL_AM_I_ROOT()) PRINT*,'Latitude is not 40-60S: Did not call 40-60S param!'
      ENDIF

      IF(BLAT.LT.ALATS(3).AND.BLAT.GE.ALATS(2)) THEN
        DO J=1,12
         IF(COCOUNT.EQ.DLATS3(J)) GOTO 302 
        ENDDO
   !   IF (MAPL_AM_I_ROOT()) PRINT*,'Latitude is not 30-40S: Did not call 30-40S param!'
      ENDIF

      IF(BLAT.LT.ALATS(4).AND.BLAT.GE.ALATS(3)) THEN
        DO J=1,12
         IF(COCOUNT.EQ.DLATS4(J)) GOTO 302 
        ENDDO
    !  IF (MAPL_AM_I_ROOT()) PRINT*,'Latitude is not 0-30S: Did not call 0-30S param!'
      ENDIF

      IF(BLAT.LT.ALATS(5).AND.BLAT.GE.ALATS(4)) THEN
        DO J=1,12
         IF(COCOUNT.EQ.DLATS5(J)) GOTO 302 
        ENDDO
     ! IF (MAPL_AM_I_ROOT()) PRINT*,'Latitude is not 0-30N: Did not call 0-30N param!'
      ENDIF

      IF(BLAT.LT.ALATS(6).AND.BLAT.GE.ALATS(5)) THEN
        DO J=1,12
         IF(COCOUNT.EQ.DLATS6(J)) GOTO 302 
        ENDDO
      !IF (MAPL_AM_I_ROOT()) PRINT*,'Latitude is not 30-40N: Did not call 30-40N param!'
      ENDIF

      IF(BLAT.LT.ALATS(7).AND.BLAT.GE.ALATS(6)) THEN
        DO J=1,9
         IF(COCOUNT.EQ.DLATS7(J)) GOTO 302 
        ENDDO
      !IF (MAPL_AM_I_ROOT()) PRINT*,'Latitude is not 40-60N: Did not call 40-60N param!'
      ENDIF

      IF(BLAT.GE.ALATS(7)) THEN
        DO J=1,6
         IF(COCOUNT.EQ.DLATS8(J)) GOTO 302 
        ENDDO
      !IF (MAPL_AM_I_ROOT()) PRINT*,'Latitude is not gt 60N: Did not call 60-80N param!'
      ENDIF

        PRINT*,'Stopped in SR PARAMA.'
        PRINT*,'COCOUNT = ',COCOUNT
        STOP
!
! End Error Check.
!
!**********************************************************************
!

 302  RETURN

      END SUBROUTINE PARAMA
!**********************************************************************
      SUBROUTINE PARAMB(BLAT,BSEASON,COCOUNT,NSEAS,NLATS,ALATS)
!**********************************************************************
      IMPLICIT NONE
      INTEGER K,M,J,BSEASON,COCOUNT,NSEAS,NLATS
      REAL*8 BLAT,ALATS(NLATS)
!
! Information used in Error Check.
!
      INTEGER NERROR,CSEASON1(18),CSEASON2(21),CSEASON3(21), &
     &        CSEASON4(21),NSTOP
      REAL*8 DLATS1(6),DLATS2(12),DLATS3(12), &
     & DLATS4(12),DLATS5(12),DLATS6(12),DLATS7(9),DLATS8(6)
      DATA CSEASON1/1,3,7,11,15,19,28,30,34,38,42,46,55,57,61,65, &
     &              69,73/
      DATA CSEASON2/4,8,12,16,20,23,26,31,35,39,43,47,50,53,58,62, &
     &              66,70,74,77,80/
      DATA CSEASON3/5,9,13,17,21,24,27,32,36,40,44,48,51,54, &
     &              59,63,67,71,75,78,81/
      DATA CSEASON4/2,6,10,14,18,22,25,29,33,37,41,45,49,52,56, &
     &              60,64,48,72,76,79/
      DATA DLATS1/1,2,28,29,55,56/
      DATA DLATS2/3,30,57,4,31,58,5,32,50,6,33,60/
      DATA DLATS3/7,34,61,8,35,62,9,36,63,10,37,64/
      DATA DLATS4/11,38,65,12,39,66,13,40,67,14,41,68/
      DATA DLATS5/15,42,69,16,43,70,17,44,71,18,45,72/
      DATA DLATS6/19,46,73,20,47,74,21,48,75,22,49,76/
      DATA DLATS7/23,50,77,24,51,78,25,52,79/
      DATA DLATS8/27,54,81,26,53,80/
!
!**********************************************************************
! Created by Bryan Duncan.
!**********************************************************************
! SR PARAMB is called by SR REMOTE.  Based on season and latitude, this
! subroutine choses the appropriate identification number of the 
! parameterization required.  
!
!  ----------------------------------------------------------
!  | Latitude   |   DJF    |   MAM    |   JJA    |   SON    |
!  ----------------------------------------------------------
!  |   > 80 N   |   AVG    |   AVG    |   AVG    |   AVG    |
!  |  60-80 N   |   AVG    | 26 43 80 | 27 54 81 |   AVG    |
!  |  40-60 N   |   AVG    | 23 50 77 | 24 51 78 | 25 52 79 |
!  |  30-40 N   | 19 46 73 | 20 47 74 | 21 48 75 | 22 49 76 |
!  |   0-30 N   | 15 42 69 | 16 43 70 | 17 44 71 | 18 45 72 |
!  |   0-30 S   | 11 38 65 | 12 39 66 | 13 40 67 | 14 41 68 |
!  |  30-40 S   |  7 34 61 |  8 35 62 |  9 36 63 | 10 37 64 |
!  |  40-60 S   |  3 30 57 |  4 31 58 |  5 32 59 |  6 33 60 |
!  |  60-80 S   |  1 28 55 |   AVG    |   AVG    |  2 29 56 |
!  |   > 80 S   |   AVG    |   AVG    |   AVG    |    AVG   |
!  ----------------------------------------------------------
!
! The table above illustrates the spatial and temporal coverage of the
! parameterizations found in the "REMOTE" subdomain.  In the boxes labeled 
! "AVG", the OH is sufficiently low so that the climatological mean OH for
! that season and latitude is used (Spivakovsky et al., 1999).  (For these 
! boxes, the parameterized OH is not calculated.)  For all other boxes,
! the parameterized OH is calculated.  The three numbers in the box correspond
! to the identification numbers of the parameterization for that box.
! The first number is for the surface layer, the second for the middle 
! tropospheric layer and the third for the upper tropospheric layer.
! For instance, the parameterization the code would use to predict OH for
! April, 15N, has the identification number 16.
!**********************************************************************
! List of Variables & Arrays
!
! COCOUNT = identification number of parameterization. 
!
!**********************************************************************
! Loop over latitude bands and then seasons to find appropriate
! parameterization.
!
      DO K=1,NLATS
!
!***************
!
        IF(BLAT.GE.ALATS(K).AND.K.NE.NLATS) THEN
          GOTO 298
        ENDIF
!
           NSTOP=0
           NERROR=1
!
        IF(BLAT.LT.ALATS(1)) THEN
          NERROR=0
          IF(BSEASON.EQ.2) NSTOP=1
          IF(BSEASON.EQ.3) NSTOP=1
        ENDIF 
        IF(BLAT.LT.ALATS(7).AND.BLAT.GE.ALATS(1)) THEN
         COCOUNT=(K-1)*NSEAS-2
        ENDIF 
        IF(BLAT.GE.ALATS(7)) THEN
          IF(BSEASON.EQ.1) NSTOP=1
          IF(BSEASON.EQ.4) NSTOP=1
         COCOUNT=(K)*NSEAS-2-1
        ENDIF 
!
!**********************************************************************
!
! Error Check.
!
        IF(NERROR.EQ.1.AND.COCOUNT.EQ.0) THEN
          NSTOP=1
        ENDIF
!
! NSTOP is an error check to see if a parameterization which
!       does not exist has been called.  This error indicates
!       that there is a problem in SR SKIPS.

        IF(NSTOP.EQ.1) THEN
          PRINT*,'***********************************'
          PRINT*,'Stopped in SR PARAMB/SR REMOTE!'
          PRINT*,'A parameterization which does not exist has been'
          PRINT*,'called.'
          PRINT*,'Season=',BSEASON
          PRINT*,'Latitude=',BLAT
          PRINT*,'***********************************'
          STOP
        ENDIF
!
! End Error Check.
!
!**********************************************************************
! Loop over seasons.
!
        DO M=1,NSEAS
!
!***********
!
        IF(BLAT.LT.ALATS(7).AND.BLAT.GE.ALATS(6).AND.M.EQ.1) GOTO 299 
        IF(BLAT.LT.ALATS(1).AND.(M.EQ.2.OR.M.EQ.3)) GOTO 299 
        IF(BLAT.GE.ALATS(7).AND.M.EQ.1) GOTO 299
        IF(BLAT.GE.ALATS(7).AND.M.EQ.4) GOTO 299

          IF(BSEASON.GT.M) THEN
            COCOUNT=COCOUNT+1
            GOTO 299
          ENDIF
!
                 COCOUNT=COCOUNT+1
!
                 GOTO 300
!
 299        CONTINUE
!***********
        ENDDO ! End Loop over Seasons.
!***********
 298        CONTINUE
!***************
      ENDDO   ! End Loop over Latitude Bands.
!***************
!
 300  CONTINUE
!
!**********************************************************************
!
! Error Check.
!
      IF(BSEASON.EQ.1) THEN
        DO J=1,18
         IF(COCOUNT.EQ.CSEASON1(J)) GOTO 301
        ENDDO
        PRINT*,'Season is DJF:  Did not call a DJF parameterization!'
      ENDIF
      IF(BSEASON.EQ.2) THEN
        DO J=1,21
         IF(COCOUNT.EQ.CSEASON2(J)) GOTO 301
        ENDDO
        PRINT*,'Season is MAM:  Did not call a MAM parameterization!'
      ENDIF
      IF(BSEASON.EQ.3) THEN
        DO J=1,21
         IF(COCOUNT.EQ.CSEASON3(J)) GOTO 301
        ENDDO
        PRINT*,'Season is JJA:  Did not call a JJA parameterization!'
      ENDIF
      IF(BSEASON.EQ.4) THEN
        DO J=1,21
         IF(COCOUNT.EQ.CSEASON4(J)) GOTO 301
        ENDDO
        PRINT*,'Season is SON:  Did not call a SON parameterization!'
      ENDIF

        PRINT*,'Stopped in SR PARAMA.'
        PRINT*,'COCOUNT = ',COCOUNT
        STOP
!
 301  CONTINUE
!
      IF(BLAT.LT.ALATS(1)) THEN
        DO J=1,9
         IF(COCOUNT.EQ.DLATS1(J)) GOTO 302
        ENDDO
      PRINT*,'Latitude is not lt 60S: Did not call 60-80S param!'
      ENDIF

      IF(BLAT.LT.ALATS(2).AND.BLAT.GE.ALATS(1)) THEN
        DO J=1,12
         IF(COCOUNT.EQ.DLATS2(J)) GOTO 302
        ENDDO
      PRINT*,'Latitude is not 40-60S: Did not call 40-60S param!'
      ENDIF

      IF(BLAT.LT.ALATS(3).AND.BLAT.GE.ALATS(2)) THEN
        DO J=1,12
         IF(COCOUNT.EQ.DLATS3(J)) GOTO 302
        ENDDO
      PRINT*,'Latitude is not 30-40S: Did not call 30-40S param!'
      ENDIF

      IF(BLAT.LT.ALATS(4).AND.BLAT.GE.ALATS(3)) THEN
        DO J=1,12
         IF(COCOUNT.EQ.DLATS4(J)) GOTO 302
        ENDDO
      PRINT*,'Latitude is not 0-30S: Did not call 0-30S param!'
      ENDIF

      IF(BLAT.LT.ALATS(5).AND.BLAT.GE.ALATS(4)) THEN
        DO J=1,12
         IF(COCOUNT.EQ.DLATS5(J)) GOTO 302
        ENDDO
      PRINT*,'Latitude is not 0-30N: Did not call 0-30N param!'
      ENDIF

      IF(BLAT.LT.ALATS(6).AND.BLAT.GE.ALATS(5)) THEN
        DO J=1,12
         IF(COCOUNT.EQ.DLATS6(J)) GOTO 302
        ENDDO
      PRINT*,'Latitude is not 30-40N: Did not call 30-40N param!'
      ENDIF

      IF(BLAT.LT.ALATS(7).AND.BLAT.GE.ALATS(6)) THEN
        DO J=1,9
         IF(COCOUNT.EQ.DLATS7(J)) GOTO 302
        ENDDO
      PRINT*,'Latitude is not 40-60N: Did not call 40-60N param!'
      ENDIF

      IF(BLAT.GE.ALATS(7)) THEN
        DO J=1,6
         IF(COCOUNT.EQ.DLATS8(J)) GOTO 302
        ENDDO
      PRINT*,'Latitude is not gt 60N: Did not call 60-80N param!'
      ENDIF

        PRINT*,'Stopped in SR PARAMA.'
        PRINT*,'COCOUNT = ',COCOUNT
        STOP
!
! End Error Check.
!
!**********************************************************************
!
 302  RETURN
      END SUBROUTINE PARAMB
!**********************************************************************
      SUBROUTINE PARAMC(BLAT,BSEASON,COCOUNT,NSEAS,NLATS,ALATS)
!**********************************************************************
      IMPLICIT NONE
      INTEGER K,M,J,BSEASON,COCOUNT,NSEAS,NLATS
      REAL*8 BLAT,ALATS(NLATS)
!
! Information used in Error Check.
!
      INTEGER NSTOP,NERROR,CSEASON1(5),CSEASON2(5),CSEASON3(6), &
     &        CSEASON4(6)
      REAL*8 DLATS2(2),DLATS3(4),DLATS4(4),DLATS5(4),DLATS6(4), &
     &        DLATS7(3),DLATS8(1)
      DATA CSEASON1/1,3,7,11,15/
      DATA CSEASON2/4,8,12,16,19/
      DATA CSEASON3/5,9,13,17,20,22/
      DATA CSEASON4/2,6,10,14,18,21/
      DATA DLATS2/1,2/
      DATA DLATS3/3,4,5,6/
      DATA DLATS4/7,8,9,10/
      DATA DLATS5/11,12,13,14/
      DATA DLATS6/15,16,17,18/
      DATA DLATS7/19,20,21/
      DATA DLATS8/22/
!**********************************************************************
! Created by Bryan Duncan.
!**********************************************************************
! SR PARAMC is called by SR CLEANISOP.  Based on season and latitude, this
! subroutine choses the appropriate identification number of the
! parameterization required.
!
!  ----------------------------------------------------------
!  | Latitude   |   DJF    |   MAM    |   JJA    |   SON    |
!  ----------------------------------------------------------
!  |   > 80 N   |   AVG    |   AVG    |    *     |   AVG    |
!  |  60-80 N   |   AVG    |    *     | 22 44 66 |   AVG    |
!  |  40-60 N   |   AVG    | 19 41 63 | 20 42 64 | 21 43 65 |
!  |  30-40 N   | 15 37 59 | 16 38 60 | 17 39 61 | 18 40 62 |
!  |   0-30 N   | 11 33 55 | 12 34 56 | 13 35 57 | 14 36 58 |
!  |   0-30 S   |  7 29 51 |  8 30 52 |  9 31 53 | 10 32 54 |
!  |  30-40 S   |  3 25 47 |  4 26 48 |  5 27 49 |  6 28 50 |
!  |  40-60 S   |  1 23 45 |    *     |    *     |  2 24 46 |
!  |  60-80 S   |     *    |   AVG    |   AVG    |     *    |
!  |   > 80 S   |     *    |   AVG    |   AVG    |    AVG   |
!  ----------------------------------------------------------
!
! The table above illustrates the spatial and temporal coverage of the
! parameterizations found in the "CLEANISOP" subdomain.  In the boxes labeled
! "AVG", the OH is sufficiently low so that the climatological mean OH for
! that season and latitude is used (Spivakovsky et al., 1999).  (For these
! boxes, the parameterized OH is not calculated.)  For all other boxes, the
! parameterized OH is calculated. The 3 numbers in the box correspond
! to the identification numbers of the parameterization for that box.
! The first number is for the surface layer, the second for the middle
! tropospheric layer and the third for the upper tropospheric layer.
! For instance, the parameterization the code would use to predict OH for
! April, 15N, middle troposphere has the identification number 34.
! The boxes designated by an asterisk are not covered by this domain.
!**********************************************************************
! List of Variables & Arrays
!
! COCOUNT = identification number of parameterization.
!
!**********************************************************************
! Loop over latitudes.
!
      DO K=2,NLATS
!
!***************
!
        IF(BLAT.GE.ALATS(K).AND.K.NE.NLATS) THEN
          GOTO 298
        ENDIF
!
           NSTOP=0
           NERROR=1
!
        IF(BLAT.LT.ALATS(1)) NSTOP=1
        IF(BLAT.LT.ALATS(2).AND.BLAT.GE.ALATS(1)) THEN
           NERROR=0
          IF(BSEASON.EQ.2) NSTOP=1
          IF(BSEASON.EQ.3) NSTOP=1
        ENDIF
!
        IF(BLAT.LT.ALATS(3).AND.BLAT.GE.ALATS(2)) THEN
          COCOUNT=2
        ENDIF
!
        IF(BLAT.LT.ALATS(7).AND.BLAT.GE.ALATS(3)) THEN
          COCOUNT=(K-1)*NSEAS-6
          IF(BLAT.GE.ALATS(6)) THEN
!              COCOUNT=COCOUNT-1
              IF(BSEASON.EQ.1) NSTOP=1
          ENDIF
        ENDIF
        IF(BLAT.GE.ALATS(7)) THEN
          IF(BSEASON.NE.3) NSTOP=1
          COCOUNT=(K-1)*NSEAS-2-1
        ENDIF
!
! Error Check.
!
        IF(NERROR.EQ.1.AND.COCOUNT.EQ.0) THEN
          IF (MAPL_AM_I_ROOT()) PRINT*,'Latitude band not found!'
          NSTOP=1
        ENDIF
!
! NSTOP is an error check to see if a parameterization which
!       does not exist has been called.  This error indicates
!       that there is a problem in SR SKIPS.

        IF(NSTOP.EQ.1) THEN
          PRINT*,'STOPPED IN SR PARAMC'
          STOP
        ENDIF
!
!***********
! Loop over seasons.
!
        DO M=1,NSEAS
!
!***********
        IF(BLAT.LT.ALATS(2).AND.BLAT.GE.ALATS(1).AND. &
     &       (M.EQ.2.OR.M.EQ.3)) GOTO 299
        IF(BLAT.LT.ALATS(7).AND.BLAT.GE.ALATS(6).AND.M.EQ.1) GOTO 299
        IF(BLAT.GE.ALATS(7).AND.M.NE.3) GOTO 299
!
          IF(BSEASON.GT.M) THEN
            COCOUNT=COCOUNT+1
            GOTO 299
          ENDIF
!
                 COCOUNT=COCOUNT+1
!
                 GOTO 300
!
 299        CONTINUE
!***********
        ENDDO ! End Loop over Seasons.
!***********
 298     CONTINUE
!***************
      ENDDO   ! End Loop over Latitude Bands.
!***************
!
 300        CONTINUE
!
! Error Check.

      IF(BSEASON.EQ.1) THEN
        DO J=1,5
         IF(COCOUNT.EQ.CSEASON1(J)) GOTO 301
         IF(COCOUNT.EQ.CSEASON1(J)+22) GOTO 301
         IF(COCOUNT.EQ.CSEASON1(J)+44) GOTO 301
        ENDDO
        PRINT*,'Season is DJF:  Did not call a DJF param!'
      ENDIF
      IF(BSEASON.EQ.2) THEN
        DO J=1,5
         IF(COCOUNT.EQ.CSEASON2(J)) GOTO 301
         IF(COCOUNT.EQ.CSEASON2(J)+22) GOTO 301
         IF(COCOUNT.EQ.CSEASON2(J)+44) GOTO 301
        ENDDO
        PRINT*,'Season is MAM:  Did not call a MAM param!'
      ENDIF
      IF(BSEASON.EQ.3) THEN
        DO J=1,6
         IF(COCOUNT.EQ.CSEASON3(J)) GOTO 301
         IF(COCOUNT.EQ.CSEASON3(J)+22) GOTO 301
         IF(COCOUNT.EQ.CSEASON3(J)+44) GOTO 301
        ENDDO
        PRINT*,'Season is JJA:  Did not call a JJA param!'
      ENDIF
      IF(BSEASON.EQ.4) THEN
        DO J=1,6
         IF(COCOUNT.EQ.CSEASON4(J)) GOTO 301
         IF(COCOUNT.EQ.CSEASON4(J)+22) GOTO 301
         IF(COCOUNT.EQ.CSEASON4(J)+44) GOTO 301
        ENDDO
        PRINT*,'Season is SON:  Did not call a SON param!'
      ENDIF

        PRINT*,'Stopped in SR paramC.'
        PRINT*,'COCOUNT = ',COCOUNT
        STOP

 301  CONTINUE
!
      IF(BLAT.LT.ALATS(2).AND.BLAT.GE.ALATS(1)) THEN
        DO J=1,2
         IF(COCOUNT.EQ.DLATS2(J)) GOTO 302
        ENDDO
        PRINT*,'Latitude is not 40-60S: Did not call a 40-60S param!'
      ENDIF

      IF(BLAT.LT.ALATS(3).AND.BLAT.GE.ALATS(2)) THEN
        DO J=1,4
         IF(COCOUNT.EQ.DLATS3(J)) GOTO 302
        ENDDO
        PRINT*,'Latitude is not 30-40S: Did not call a 30-40S param!'
      ENDIF

      IF(BLAT.LT.ALATS(4).AND.BLAT.GE.ALATS(3)) THEN
        DO J=1,4
         IF(COCOUNT.EQ.DLATS4(J)) GOTO 302
        ENDDO
        PRINT*,'Latitude is not 0-30S: Did not call a 0-30S param!'
      ENDIF

      IF(BLAT.LT.ALATS(5).AND.BLAT.GE.ALATS(4)) THEN
        DO J=1,4
         IF(COCOUNT.EQ.DLATS5(J)) GOTO 302
        ENDDO
        PRINT*,'Latitude is not 0-30N: Did not call a 0-30N param!'
      ENDIF

      IF(BLAT.LT.ALATS(6).AND.BLAT.GE.ALATS(5)) THEN
        DO J=1,4
         IF(COCOUNT.EQ.DLATS6(J)) GOTO 302
        ENDDO
        PRINT*,'Latitude is not 30-40N: Did not call a 30-40N param!'
      ENDIF

      IF(BLAT.LT.ALATS(7).AND.BLAT.GE.ALATS(6)) THEN
        DO J=1,3
         IF(COCOUNT.EQ.DLATS7(J)) GOTO 302
        ENDDO
        PRINT*,'Latitude is not 40-60N: Did not call a 40-60N param!'
      ENDIF

      IF(BLAT.GE.ALATS(7)) THEN
        DO J=1,1
         IF(COCOUNT.EQ.DLATS8(J)) GOTO 302
        ENDDO
        PRINT*,'Latitude is not gt 60N: Did not call a 60-80N param!'
      ENDIF

        PRINT*,'Stopped in SR paramC.'
        PRINT*,'COCOUNT = ',COCOUNT
        STOP
!
 302  RETURN

      END SUBROUTINE PARAMC
!**********************************************************************
      SUBROUTINE PARAMD(BLAT,BSEASON,COCOUNT,NSEAS,NLATS,ALATS)
!**********************************************************************
      IMPLICIT NONE
      INTEGER K,M,J,BSEASON,COCOUNT
      INTEGER NSEAS,NLATS
      REAL*8 BLAT,ALATS(NLATS)
! 
! Information used in Error Check.
!
      REAL*8 DLATS2(2),DLATS3(4),DLATS4(4),DLATS5(4),DLATS6(4), &
     &       DLATS7(3),DLATS8(1)
      INTEGER NSTOP,CSEASON1(5),CSEASON2(5),CSEASON3(6),CSEASON4(6)
      DATA CSEASON1/1,3,7,11,15/
      DATA CSEASON2/4,8,12,16,19/
      DATA CSEASON3/5,9,13,17,20,22/
      DATA CSEASON4/2,6,10,14,18,21/
      DATA DLATS2/1,2/
      DATA DLATS3/3,4,5,6/
      DATA DLATS4/7,8,9,10/
      DATA DLATS5/11,12,13,14/
      DATA DLATS6/15,16,17,18/
      DATA DLATS7/19,20,21/
      DATA DLATS8/22/
!**********************************************************************
! Created by Bryan Duncan.
!**********************************************************************
! SR PARAMD is called by SR MODPOLLUTEDISOP.  Based on season and latitude,
! this subroutine choses the appropriate identification number of the
! parameterization required for the surface layer only.
!
!  ----------------------------------------------------------
!  | Latitude   |   DJF    |   MAM    |   JJA    |   SON    |
!  ----------------------------------------------------------
!  |   > 80 N   |   AVG    |   AVG    |    *     |   AVG    |
!  |  60-80 N   |   AVG    |    *     |    22    |   AVG    |
!  |  40-60 N   |   AVG    |    19    |    20    |    21    |
!  |  30-40 N   |     15   |    16    |    17    |    18    |
!  |   0-30 N   |     11   |    12    |    13    |    14    |
!  |   0-30 S   |     7    |    8     |    9     |    10    |
!  |  30-40 S   |     3    |    4     |    5     |    6     |
!  |  40-60 S   |     1    |    *     |    *     |    2     |
!  |  60-80 S   |     *    |   AVG    |   AVG    |    *     |
!  |   > 80 S   |     *    |   AVG    |   AVG    |   AVG    |
!  ----------------------------------------------------------
!
! The table above illustrates the spatial and temporal coverage of the
! parameterizations found in the "MODPOLLUTEDISOP" subdomain. In boxes labeled
! "AVG", the OH is sufficiently low so that the climatological mean OH for
! that season and latitude is used (Spivakovsky et al., 1999).  (For these
! boxes, the parameterized OH is not calculated.)  For all other boxes, the
! parameterized OH is calculated. The number in the box corresponds
! to the identification numbers of the parameterization for that box.
! For instance, the parameterization the code would use to predict OH for
! 15N in April has the identification number 12.
! The boxes designated by an asterisk are not covered by this domain.
!**********************************************************************
! List of Variables & Arrays
!
! COCOUNT = identification number of parameterization.
!
!**********************************************************************
! Loop over latitudes.
!

      DO K=2,NLATS
!
!***************
!
        IF(BLAT.GE.ALATS(K).AND.K.NE.NLATS) THEN
          GOTO 298
        ENDIF
!
           NSTOP=0
!
        IF(BLAT.LT.ALATS(2).AND.BLAT.GE.ALATS(1)) THEN
          IF(BSEASON.EQ.2) NSTOP=1
          IF(BSEASON.EQ.3) NSTOP=1
        ENDIF
        IF(BLAT.LT.ALATS(3).AND.BLAT.GE.ALATS(2)) THEN
          COCOUNT=2
        ENDIF
        IF(BLAT.LT.ALATS(7).AND.BLAT.GE.ALATS(3)) THEN
         COCOUNT=(K-1)*NSEAS-6
          IF(BLAT.GE.ALATS(6)) THEN
              COCOUNT=COCOUNT-1
              IF(BSEASON.EQ.1) NSTOP=1
          ENDIF
        ENDIF
        IF(BLAT.GE.ALATS(7)) THEN
          IF(BSEASON.NE.3) NSTOP=1
          COCOUNT=(K-1)*NSEAS-2-1
        ENDIF
!
! NSTOP is an error check to see if a parameterization which
!       does not exist has been called.  This error indicates
!       that there is a problem in SR SKIPS.

        IF(NSTOP.EQ.1) THEN
          PRINT*,'STOPPED IN SR PARAMD'
          STOP
        ENDIF
!
!***********
! Loop over seasons.
!
        DO M=1,NSEAS
!
!***********
!
        IF(BLAT.LT.ALATS(2).AND.BLAT.GE.ALATS(1).AND. &
     &       (M.EQ.2.OR.M.EQ.3)) GOTO 299
        IF(BLAT.GE.ALATS(7).AND.M.NE.3) GOTO 299
!
          IF(BSEASON.GT.M) THEN
            COCOUNT=COCOUNT+1
            GOTO 299
          ENDIF
!
                 COCOUNT=COCOUNT+1
!
                 GOTO 300

 299        CONTINUE
!***********
        ENDDO  ! End Loop over Seasons.
!***********
 298        CONTINUE
!***************
      ENDDO    ! End Loop over Latitude Bands.
!***************
!
 300  CONTINUE
!
!**********************************************************************
! Error Check.
!
      IF(BSEASON.EQ.1) THEN
        DO J=1,5
         IF(COCOUNT.EQ.CSEASON1(J)) GOTO 301
        ENDDO
        PRINT*,'Season is DJF:  Did not call a DJF param!'
      ENDIF
!
      IF(BSEASON.EQ.2) THEN
        DO J=1,5
         IF(COCOUNT.EQ.CSEASON2(J)) GOTO 301
        ENDDO
        PRINT*,'Season is MAM:  Did not call a MAM param!'
      ENDIF
!
      IF(BSEASON.EQ.3) THEN
        DO J=1,6
         IF(COCOUNT.EQ.CSEASON3(J)) GOTO 301
        ENDDO
        PRINT*,'Season is JJA:  Did not call a JJA param!'
      ENDIF
!
      IF(BSEASON.EQ.4) THEN
        DO J=1,6
         IF(COCOUNT.EQ.CSEASON4(J)) GOTO 301
        ENDDO
        PRINT*,'Season is SON:  Did not call a SON param!'
      ENDIF
!
        PRINT*,'Stopped in SR paramD.'
        PRINT*,'COCOUNT = ',COCOUNT
        STOP
!
 301  CONTINUE
!
      IF(BLAT.LT.ALATS(2).AND.BLAT.GE.ALATS(1)) THEN
        DO J=1,2
         IF(COCOUNT.EQ.DLATS2(J)) GOTO 302
        ENDDO
        PRINT*,'Latitude is not 40-60S: Did not call a 40-60S param!'
      ENDIF
!
      IF(BLAT.LT.ALATS(3).AND.BLAT.GE.ALATS(2)) THEN
        DO J=1,4
         IF(COCOUNT.EQ.DLATS3(J)) GOTO 302
        ENDDO
        PRINT*,'Latitude is not 30-40S: Did not call a 30-40S param!'
      ENDIF
!
      IF(BLAT.LT.ALATS(4).AND.BLAT.GE.ALATS(3)) THEN
        DO J=1,4
         IF(COCOUNT.EQ.DLATS4(J)) GOTO 302
        ENDDO
        PRINT*,'Latitude is not 0-30S: Did not call a 0-30S param!'
      ENDIF
!
      IF(BLAT.LT.ALATS(5).AND.BLAT.GE.ALATS(4)) THEN
        DO J=1,4
         IF(COCOUNT.EQ.DLATS5(J)) GOTO 302
        ENDDO
        PRINT*,'Latitude is not 0-30N: Did not call a 0-30N param!'
      ENDIF
!
      IF(BLAT.LT.ALATS(6).AND.BLAT.GE.ALATS(5)) THEN
        DO J=1,4
         IF(COCOUNT.EQ.DLATS6(J)) GOTO 302
        ENDDO
        PRINT*,'Latitude is not 30-40N: Did not call a 30-40N param!'
      ENDIF
!
      IF(BLAT.LT.ALATS(7).AND.BLAT.GE.ALATS(6)) THEN
        DO J=1,3
         IF(COCOUNT.EQ.DLATS7(J)) GOTO 302
        ENDDO
        PRINT*,'Latitude is not 40-60N: Did not call a 40-60N param!'
      ENDIF
!
      IF(BLAT.GE.ALATS(7)) THEN
        DO J=1,1
         IF(COCOUNT.EQ.DLATS8(J)) GOTO 302
        ENDDO
        PRINT*,'Latitude is not gt 60N: Did not call a 60-80N param!'
      ENDIF
!
        PRINT*,'Stopped in SR paramD.'
        PRINT*,'COCOUNT = ',COCOUNT
        STOP
!
!**********************************************************************
!
! EY DEBUG
 302  RETURN

      END SUBROUTINE PARAMD
!**********************************************************************
      SUBROUTINE PARAME(BLAT,BSEASON,COCOUNT,NSEAS,NLATS,ALATS)
!**********************************************************************
      IMPLICIT NONE
      INTEGER K,M,BSEASON,COCOUNT
      INTEGER NSEAS,NLATS,NSTOP
      REAL*8 BLAT,ALATS(NLATS)
!**********************************************************************
! Created by Bryan Duncan.
!**********************************************************************
! SR PARAME is called by SR HIPOLLUTEDISOP.  Based on season and latitude,
! this subroutine choses the appropriate identification number of the
! parameterization required for surface only.
!
!  ----------------------------------------------------------
!  | Latitude   |   DJF    |   MAM    |   JJA    |   SON    |
!  ----------------------------------------------------------
!  |   > 80 N   |   AVG    |   AVG    |    *     |   AVG    |
!  |  60-80 N   |   AVG    |    *     |    *     |   AVG    |
!  |  40-60 N   |   AVG    |    *     |   12     |   13     |
!  |  30-40 N   |     *    |    9     |   10     |   11     |
!  |   0-30 N   |     5    |    6     |    7     |    8     |
!  |   0-30 S   |     1    |    2     |    3     |    4     |
!  |  30-40 S   |     *    |    *     |    *     |    *     |
!  |  40-60 S   |     *    |    *     |    *     |    *     |
!  |  60-80 S   |     *    |   AVG    |   AVG    |    *     |
!  |   > 80 S   |     *    |   AVG    |   AVG    |   AVG    |
!  ----------------------------------------------------------
!
! The table above illustrates the spatial and temporal coverage of the
! parameterizations found in the "HIPOLLUTEDISOP" subdomain. In boxes labeled
! "AVG", the OH is sufficiently low so that the climatological mean OH for
! that season and latitude is used (Spivakovsky et al., 1999).  (For these
! boxes, the parameterized OH is not calculated.)  For all other boxes,
! the parameterized OH is calculated.  The number in the box corresponds
! to the identification number of the parameterization for that box.
! For instance, the parameterization the code would use to predict OH for
! 35N in April has the identification number 3.
! The boxes designated by an asterisk are not covered by this domain.
!**********************************************************************
! List of Variables & Arrays
!
! COCOUNT = identification number of parameterization.
!
!**********************************************************************
!
          NSTOP=0
!
        IF(BLAT.LE.ALATS(3)) NSTOP=1 
        IF(BLAT.GE.ALATS(7)) NSTOP=1
!
        IF(BLAT.LT.ALATS(7).AND.BLAT.GE.ALATS(6)) THEN
          IF(BSEASON.LE.2) NSTOP=1
          IF(BSEASON.EQ.3) COCOUNT=12
          IF(BSEASON.EQ.4) COCOUNT=13
        ENDIF
        IF(BLAT.LT.ALATS(6).AND.BLAT.GE.ALATS(5)) THEN
          IF(BSEASON.EQ.1) NSTOP=1
          IF(BSEASON.EQ.2) COCOUNT=9
          IF(BSEASON.EQ.3) COCOUNT=10
          IF(BSEASON.EQ.4) COCOUNT=11
        ENDIF
        IF(BLAT.LT.ALATS(5).AND.BLAT.GE.ALATS(4)) THEN
          IF(BSEASON.EQ.1) COCOUNT=5
          IF(BSEASON.EQ.2) COCOUNT=6
          IF(BSEASON.EQ.3) COCOUNT=7
          IF(BSEASON.EQ.4) COCOUNT=8
        ENDIF
        IF(BLAT.LT.ALATS(4).AND.BLAT.GE.ALATS(3)) THEN
          IF(BSEASON.EQ.1) COCOUNT=1
          IF(BSEASON.EQ.2) COCOUNT=2
          IF(BSEASON.EQ.3) COCOUNT=3
          IF(BSEASON.EQ.4) COCOUNT=4
        ENDIF
!
!**********************************************************************
! Error Check.
!
! NSTOP is an error check to see if a parameterization which
!       does not exist has been called.  This error indicates
!       that there is a problem in SR SKIPS.
!
        IF(NSTOP.EQ.1) THEN
          PRINT*,'STOPPED IN SR PARAME'
          STOP
        ENDIF
!
      IF(BSEASON.EQ.1)THEN
         IF(COCOUNT.NE.1.OR.COCOUNT.NE.5) GOTO 301
        PRINT*,'Season is DJF:  Did not call a DJF param!'
      ENDIF
      IF(BSEASON.EQ.2) THEN
         IF(COCOUNT.EQ.2.OR.COCOUNT.EQ.6.OR.COCOUNT.EQ.9) GOTO 301
        PRINT*,'Season is MAM:  Did not call a MAM param!'
      ENDIF
!
      IF(BSEASON.EQ.3) THEN
         IF(COCOUNT.EQ.3.OR.COCOUNT.EQ.7.OR.COCOUNT.EQ.10.OR.COCOUNT.EQ.12) &
            GOTO 301
        PRINT*,'Season is JJA:  Did not call a JJA param!'
      ENDIF
!
      IF(BSEASON.EQ.4) THEN
         IF(COCOUNT.EQ.4.OR.COCOUNT.EQ.8.OR.COCOUNT.EQ.11.OR.COCOUNT.EQ.13) &
            GOTO 301
        PRINT*,'Season is SON:  Did not call a SON param!'
      ENDIF
!
        PRINT*,'Stopped in SR paramE.'
        PRINT*,'COCOUNT = ',COCOUNT
        STOP
!
 301  CONTINUE
!
      IF(BLAT.LT.ALATS(7).AND.BLAT.GE.ALATS(6)) THEN
       IF(COCOUNT.GE.12) GOTO 302
      ENDIF
      IF(BLAT.LT.ALATS(6).AND.BLAT.GE.ALATS(5)) THEN
       IF(COCOUNT.GE.9.AND.COCOUNT.LE.11) GOTO 302
      ENDIF
      IF(BLAT.LT.ALATS(5).AND.BLAT.GE.ALATS(4)) THEN
       IF(COCOUNT.GE.5.AND.COCOUNT.LE.8) GOTO 302
      ENDIF
      IF(BLAT.LT.ALATS(4).AND.BLAT.GE.ALATS(3)) THEN
       IF(COCOUNT.LE.4) GOTO 302
      ENDIF
!
        PRINT*,'Stopped in SR paramE.'
        PRINT*,'COCOUNT = ',COCOUNT
        STOP
!
!**********************************************************************
!
 302  RETURN

      END SUBROUTINE PARAME
!**********************************************************************
      SUBROUTINE PARAMF(BLAT,BSEASON,COCOUNT,NSEAS,NLATS,ALATS)
!**********************************************************************
      IMPLICIT NONE
      INTEGER K,M,BSEASON,COCOUNT
      INTEGER NSEAS,NLATS,NSTOP
      REAL*8 BLAT,ALATS(NLATS)
!**********************************************************************
! Created by Bryan Duncan.
!**********************************************************************
! SR PARAMF is called by SR MODPOLLUTEDISOP.  Based on season and latitude,
! this subroutine choses the appropriate identification number of the
! parameterization required for middle troposphere only.
!
!  ----------------------------------------------------------
!  | Latitude   |   DJF    |   MAM    |   JJA    |   SON    |
!  ----------------------------------------------------------
!  |   > 80 N   |   AVG    |   AVG    |    *     |   AVG    |
!  |  60-80 N   |   AVG    |    *     |    *     |   AVG    |
!  |  40-60 N   |   AVG    |    9     |   10     |   11     |
!  |  30-40 N   |    5     |    6     |    7     |    8     |
!  |   0-30 N   |    3     |    4     |    *     |    *     |
!  |   0-30 S   |    *     |    *     |    1     |    2     |
!  |  30-40 S   |    *     |    *     |    *     |    *     |
!  |  40-60 S   |    *     |    *     |    *     |    *     |
!  |  60-80 S   |    *     |   AVG    |   AVG    |    *     |
!  |   > 80 S   |    *     |   AVG    |   AVG    |   AVG    |
!  ----------------------------------------------------------
!
! The table above illustrates the spatial and temporal coverage of the
! parameterizations found in the "MODPOLLUTEDISOP" subdomain. In boxes labeled
! "AVG", the OH is sufficiently low so that the climatological mean OH for
! that season and latitude is used (Spivakovsky et al., 1999).  (For these
! boxes, the parameterized OH is not calculated.)  For all other boxes,
! the parameterized OH is calculated.  The number in the box corresponds
! to the identification number of the parameterization for that box.
! For instance, the parameterization the code would use to predict OH for
! 35N in April has the identification number 2.
! The boxes designated by an asterisk are not covered by this domain.
!**********************************************************************
! List of Variables & Arrays
!
! COCOUNT = identification number of parameterization.
!
!**********************************************************************
!
          NSTOP=0
!
        IF(BLAT.LE.ALATS(3)) NSTOP=1 
        IF(BLAT.GE.ALATS(7)) NSTOP=1
!
        IF(BLAT.LE.ALATS(4)) THEN
          IF(BSEASON.EQ.3) COCOUNT=1
          IF(BSEASON.EQ.4) COCOUNT=2
        ENDIF
        IF(BLAT.LE.ALATS(5).AND.BLAT.GT.ALATS(4)) THEN
          IF(BSEASON.EQ.1) COCOUNT=3
          IF(BSEASON.EQ.2) COCOUNT=4
        ENDIF
        IF(BLAT.LE.ALATS(6).AND.BLAT.GT.ALATS(5)) THEN
          IF(BSEASON.EQ.1) COCOUNT=5
          IF(BSEASON.EQ.2) COCOUNT=6
          IF(BSEASON.EQ.3) COCOUNT=7
          IF(BSEASON.EQ.4) COCOUNT=8
        ENDIF
        IF(BLAT.GT.ALATS(6)) THEN
          IF(BSEASON.EQ.1) NSTOP=1
          IF(BSEASON.EQ.2) COCOUNT=9
          IF(BSEASON.EQ.3) COCOUNT=10
          IF(BSEASON.EQ.4) COCOUNT=11
        ENDIF
!
!**********************************************************************
! Error Check.
!
! NSTOP is an error check to see if a parameterization which
!       does not exist has been called.  This error indicates
!       that there is a problem in SR SKIPS.

        IF(NSTOP.EQ.1) THEN
          PRINT*,'STOPPED IN SR PARAMF'
          STOP
        ENDIF
!
      IF(BSEASON.EQ.1) THEN
         IF(COCOUNT.EQ.3.OR.COCOUNT.EQ.5) GOTO 301
        PRINT*,'Season is DJF:  Did not call a DJF param!'
      ENDIF
!
      IF(BSEASON.EQ.2) THEN
         IF(COCOUNT.EQ.4.OR.COCOUNT.EQ.6.OR.COCOUNT.EQ.9) GOTO 301
        PRINT*,'Season is MAM:  Did not call a MAM param!'
      ENDIF
!
      IF(BSEASON.EQ.3) THEN
         IF(COCOUNT.EQ.1.OR.COCOUNT.EQ.7.OR.COCOUNT.EQ.10) GOTO 301
        PRINT*,'Season is JJA:  Did not call a JJA param!'
      ENDIF
!
      IF(BSEASON.EQ.4) THEN
         IF(COCOUNT.EQ.2.OR.COCOUNT.EQ.8.OR.COCOUNT.EQ.11) GOTO 301
        PRINT*,'Season is SON:  Did not call a SON param!'
      ENDIF
!
        PRINT*,'Stopped in SR paramF - 1.'
        PRINT*,'COCOUNT = ',COCOUNT
        STOP
!
 301  CONTINUE
!
      IF(BLAT.LE.ALATS(4)) THEN
       IF(COCOUNT.LE.2) GOTO 302
      ENDIF
      IF(BLAT.LE.ALATS(5).AND.BLAT.GT.ALATS(4)) THEN
       IF(COCOUNT.EQ.3.OR.COCOUNT.EQ.4) GOTO 302
      ENDIF
      IF(BLAT.LE.ALATS(6).AND.BLAT.GT.ALATS(5)) THEN
       IF(COCOUNT.GE.5.AND.COCOUNT.LE.8) GOTO 302
      ENDIF
      IF(BLAT.GT.ALATS(6)) THEN
       IF(COCOUNT.GE.9) GOTO 302
      ENDIF

!
        PRINT*,'Stopped in SR paramF - 2.'
        PRINT*,'COCOUNT = ',COCOUNT
        STOP
!
!**********************************************************************
!
 302  RETURN

      END SUBROUTINE PARAMF

  END MODULE computeParameters_mod
