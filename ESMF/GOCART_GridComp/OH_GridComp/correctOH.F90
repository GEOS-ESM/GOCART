!**********************************************************************
      SUBROUTINE CORRECTOH(BOH,COLON,COLAT,COALT, iLon, iLat, iVert)
!**********************************************************************
      USE cblock_cf_mod    ! "CMN_CF"
      USE cblock_OH_mod    ! "CMN_OH"

      IMPLICIT NONE

      INTEGER :: iLon, iLat, iVert
      INTEGER N,I,J,K
      INTEGER COLON,COLAT,COALT
      REAL*8 BOH(iLon, iLat, iVert)
!**********************************************************************
! Created by Bryan Duncan.     
!**********************************************************************
! This SR corrects OH to remove the problem with high
! background aerosol (tau=0.1) in the OH parameterization.
! This file contains ratios of OH calculated assuming a uniform
! tau=0.01 over the globe and assuming tau=0.1. OH calculated by
! CMS, June, 2001. Ratio created by bnd, June, 2001.
!*************************************************************************
!
!      print*,'COLON,COLAT,COALT=',COLON,COLAT,COALT
!      print*,'BOH(COLON,COLAT,COALT)=',BOH(COLON,COLAT,COALT)
!      print*,'OH_SEASON=',OH_SEASON
!      print*,'OH_LAT=',OH_LAT
!      print*,'OH_PRESS=',OH_PRESS

      N=OH_SEASON
!
      DO J=1,NCMSLATS2
!
       DO K=NCMSALTS2,1,-1
!
!         print*,'a',j,k,OH_LAT,CMSLATS2(J)
         IF(OH_LAT.GE.CMSLATS2(J)) THEN         
!
!         print*,'b',j,k,OH_PRESS,CMSALTS2(K)
           IF(CMSALTS2(K).GE.OH_PRESS) THEN
!         print*,'c',j,k,n,correction(N,J,K)
            BOH(COLON,COLAT,COALT)=BOH(COLON,COLAT,COALT)* &
     &           correction(N,J,K)
              GOTO 2
           ENDIF
!
           IF(K.EQ.1) THEN
!         print*,'d',j,k,n,correction(N,J,K)
            BOH(COLON,COLAT,COALT)=BOH(COLON,COLAT,COALT)* &
     &           correction(N,J,K)
              GOTO 2
           ENDIF
!
         ENDIF
! EY treat polar points using the Spivakovsky formulation	    
! MeM  Changed -90 to CMSLATS2(NCMSLATS)
         IF ( OH_LAT .LE. CMSLATS2(NCMSLATS2)) THEN
 	     BOH(COLON,COLAT,COALT)=BOH(COLON,COLAT,COALT)* &
     &           correction(N,J,K)	    
	     GOTO 2
	 ENDIF
!
       ENDDO
!
      ENDDO
!
! Error Check.
      PRINT*,'STOPPED IN SR correctOH!'
      PRINT*,'Point lies nowhere!'
      STOP
!
 2    CONTINUE
!         print*,'e',BOH(COLON,COLAT,COALT)
      IF(BOH(COLON,COLAT,COALT).LT.0.) THEN
         print*,COLON,COLAT,COALT,BOH(COLON,COLAT,COALT)
         print*,'went 0!'
         stop
      ENDIF
!
      RETURN
      END SUBROUTINE CORRECTOH
