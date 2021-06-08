!**********************************************************************
      SUBROUTINE INTERPOH()
!**********************************************************************
      USE cblock_OH_mod    ! "CMN_OH"

      IMPLICIT NONE

      INTEGER I,J,K
!      REAL*8 :: avgOH(:,:,:)
!**********************************************************************
! Created by Bryan Duncan.     
!**********************************************************************
! In subdomains where OH is very low or negligible due to low
!  sunlight (e.g., high latitudes in winter), concentrations of OH are
!  set to climatological mean values as a function of latitude,
!  altitude and season. This SR picks the appropriate average OH
!  field described in Spivakovsky et al., "Three-dimensional climatological
!  distribution of tropospheric OH: update and evaluation", accepted
!  to JGR, 1999.  The fields are stored in array avgOH and read in
!  in SR readavgOH.
!
! avgOH = array containing the climatological OH values.
!
! NCMSALTS = number of altitude levels of climatology
!
! NCMSLATS = number of latitude bands of climatology
!
! NOTE:  CMSLATS go from 89 to -89; CMSALTS go from 1000 to 200
!
!**********************************************************************
! 
      I=OH_SEASON
!
      DO J=1,NCMSLATS
!
         DO K=NCMSALTS,1,-1
!
!bnd
            IF (J.GE.NCMSLATS-1.AND.OH_PRESS.GE.600.) THEN
               PARAMOH=avgOH(I,J,5)*1.D5
               GOTO 2
            ENDIF
!bnd
            IF (OH_LAT.GE.CMSLATS(J)) THEN         
               IF (CMSALTS(K).GE.OH_PRESS) THEN
                  PARAMOH=avgOH(I,J,K)*1.D5
                  GOTO 2
               ENDIF

               IF (K.EQ.1 ) THEN
                  PARAMOH=avgOH(I,J,K)*1.D5
                  GOTO 2
               ENDIF
            ENDIF
! EY treat polar points using the Spivakovsky formulation	    
! MeM  Changed -90 to CMSLATS(NCMSLATS)
	    IF ( OH_LAT .LE. CMSLATS(NCMSLATS)) THEN
	       PARAMOH=avgOH(I,J,K)*1.D5
	        GOTO 2
	    ENDIF
         ENDDO
!
      ENDDO
!
! Error Check.
      PRINT*,'STOPPED IN SR interpOH!'
      PRINT*,'Point lies nowhere!'
      STOP
!
 2    CONTINUE
!bnd      PRINT*,'Got it!'
!
      RETURN

      END SUBROUTINE INTERPOH
