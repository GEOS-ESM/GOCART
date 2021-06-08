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

   MODULE  identification_mod

! !USES:

   USE ESMF
   USE MAPL

   USE computeParameters_mod
   USE cblock_OH_mod    ! "CMN_OH"
   USE computeParameters_mod   ! paramA, paramB, ..., paramF

   IMPLICIT NONE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PRIVATE
      PUBLIC :: CLEAN, CLEANISOP
      PUBLIC :: REMOTE
      PUBLIC :: MODPOLLUTED, HIPOLLUTED

!**********************************************************************
CONTAINS
!**********************************************************************
      SUBROUTINE CLEAN(NDFUNCS)
!**********************************************************************
      IMPLICIT NONE

      INTEGER NDFUNCS,NALL,K
!**********************************************************************
! Select proper identification number of parameterization needed for 
!   box characterized by unpolluted conditions and no isoprene using
!   SR PARAMA. Store this value in COCOUNT. Then select the corresponding 
!   polynomials and store in NDFUNCS.
!**********************************************************************
!   Created by Bryan Duncan.
!**********************************************************************
! A) Clean
!***************************************************
!  A1) Surface
!***************************************************
!
      ! EY DEBUG - reinitialize COCOUNT 
      
      COCOUNT = 0


      IF(OH_PRESS.GE.PRESSES(1)) THEN
!
! Determine the identification number of the appropriate parameterization.
!

        CALL PARAMA(OH_LAT,OH_SEASON,COCOUNT,NSEAS,NLATS,ALATS)
!
! Error Check.
!
         IF(COCOUNT.EQ.0) THEN
           PRINT*,'COCOUNT=0!'
           PRINT*,'No parameterization was chosen!'
           PRINT*,'Stopped in SR CLEAN.'
           STOP 
         ENDIF
!
! Determine the polynomials associated with selected identification
!   number.
!
            NALL=3*NPARAMB
         IF(COCOUNT.GT.1) THEN
          DO K=NALL+1,NALL+COCOUNT-1
           NDFUNCS=NDFUNCS+NDONE(K)
          ENDDO
         ENDIF
!
            COCOUNT=COCOUNT+NALL
!
! Error Check.
!
            NALL=NALL+NPARAMA
         IF(COCOUNT.GT.NALL) THEN
           PRINT*,'COCOUNT=',COCOUNT
           PRINT*,'NALL   =',NALL
           PRINT*,'COCOUNT > NALL'
           PRINT*,'Stopped in SR CLEAN.'
           STOP 
         ENDIF
!
            GOTO 321
!
      ENDIF
!
!***************************************************
!  END A1
!***************************************************
! Sum of all NDONE functions in A).
!
            NDFUNCS=NDFUNCS+NDFUNCSA1
!
!***************************************************
!  A2) Middle Troposphere
!***************************************************
!
      IF(OH_PRESS.LT.PRESSES(1).AND.OH_PRESS.GE.PRESSES(2)) THEN
!
! Determine the identification number of the appropriate parameterization.
!

        CALL PARAMA(OH_LAT,OH_SEASON,COCOUNT,NSEAS,NLATS,ALATS)
!
! Error Check.
!
         IF(COCOUNT.EQ.0) THEN
           PRINT*,'COCOUNT=0!'
           PRINT*,'No parameterization was chosen!'
           PRINT*,'Stopped in SR CLEAN.'
           STOP 
         ENDIF
!
! Determine the polynomials associated with selected identification
!   number.
!
	

            NALL=3*NPARAMB+NPARAMA
         IF(COCOUNT.GT.1) THEN
          DO K=NALL+1,NALL+COCOUNT-1
           NDFUNCS=NDFUNCS+NDONE(K)
          ENDDO
         ENDIF
!
          COCOUNT=COCOUNT+NALL
!
! Error Check.
!
            NALL=NALL+NPARAMA
         IF(COCOUNT.GT.NALL) THEN
           PRINT*,'COCOUNT=',COCOUNT
           PRINT*,'NALL   =',NALL
           PRINT*,'COCOUNT > NALL in SR CLEAN'
           STOP
         ENDIF

!
            GOTO 321
!
      ENDIF
!
!***************************************************
!  END A2
!***************************************************
!
            NDFUNCS=NDFUNCS+NDFUNCSA2
!
!***************************************************
!  A3) Upper Troposphere
!***************************************************
!
      IF(OH_PRESS.LT.PRESSES(2)) THEN
!
! Determine the identification number of the appropriate parameterization.
!
        CALL PARAMA(OH_LAT,OH_SEASON,COCOUNT,NSEAS,NLATS,ALATS)
!
! Error Check.
!
         IF(COCOUNT.EQ.0) THEN
           PRINT*,'COCOUNT=0!'
           PRINT*,'No parameterization was chosen!'
           PRINT*,'Stopped in SR CLEAN.'
           STOP 
         ENDIF
!
! Determine the polynomials associated with selected identification
!   number.
!
            NALL=3*NPARAMB+2*NPARAMA
         IF(COCOUNT.GT.1) THEN
          DO K=NALL+1,NALL+COCOUNT-1
           NDFUNCS=NDFUNCS+NDONE(K)
          ENDDO
         ENDIF
!
          COCOUNT=COCOUNT+NALL
!
! Error Check.
!
           NALL=NALL+NPARAMA
         IF(COCOUNT.GT.NALL) THEN
           PRINT*,'COCOUNT=',COCOUNT
           PRINT*,'NALL   =',NALL
           PRINT*,'COCOUNT > NALL in SR CLEAN'
           STOP
         ENDIF
!
            GOTO 321
!
      ENDIF
!
!***************************************************
!  END A3
!***************************************************
! Error check.
!    
      PRINT*,'********************************'
      PRINT*,'Stopped in SR CLEAN!'
      PRINT*,'No parameterization was chosen!'
      PRINT*,'********************************'
      STOP
!
!***************************************************
 321  RETURN

      END SUBROUTINE CLEAN
!**********************************************************************
      SUBROUTINE CLEANISOP(NDFUNCS)
!**********************************************************************
      IMPLICIT NONE

      INTEGER NDFUNCS,NALL,K
!**********************************************************************
! Select proper identification number of parameterization needed for
!   box characterized by unpolluted conditions with isoprene using
!   SR PARAMC. Store this value in COCOUNT. Then select the corresponding
!   polynomials and store in NDFUNCS.
!**********************************************************************
!   Created by Bryan Duncan.
!**********************************************************************
! C) Clean
!***************************************************
!  C1) Surface
!***************************************************
!
      IF(OH_PRESS.GE.PRESSES(1)) THEN
!
! Determine the identification number of the appropriate parameterization.
!
        CALL PARAMC(OH_LAT,OH_SEASON,COCOUNT,NSEAS,NLATS,ALATS)
!
! Error Check.
!
         IF(COCOUNT.EQ.0) THEN
           PRINT*,'COCOUNT=0 in SR CLEANISOP'
           STOP
         ENDIF
!
! Determine the polynomials associated with selected identification
!   number.
!
            NALL=3*NPARAMA+3*NPARAMB
         IF(COCOUNT.GT.1) THEN
          DO K=NALL+1,NALL+COCOUNT-1
           NDFUNCS=NDFUNCS+NDONE(K)
          ENDDO
         ENDIF
!
          COCOUNT=COCOUNT+NALL
!
! Error Check.
!
            NALL=NALL+NPARAMC
         IF(COCOUNT.GT.NALL) THEN
           PRINT*,'COCOUNT=',COCOUNT
           PRINT*,'NALL   =',NALL
           PRINT*,'COCOUNT > NALL in SR CLEANISOP'
           STOP
         ENDIF
!
            GOTO 321
!
      ENDIF
!
!***************************************************
!  END C1
!***************************************************
! Sum of all NDONE functions in C).
!
            NDFUNCS=NDFUNCS+NDFUNCSC1
!
!***************************************************
!  C2) Middle Troposphere
!***************************************************
!
      IF(OH_PRESS.LT.PRESSES(1).AND.OH_PRESS.GE.PRESSES(2)) THEN
!
! Determine the identification number of the appropriate parameterization.
!
        CALL PARAMC(OH_LAT,OH_SEASON,COCOUNT,NSEAS,NLATS,ALATS)
!
! Error Check.
!
         IF(COCOUNT.EQ.0) THEN
           PRINT*,'COCOUNT=0 in SR CLEANISOP'
           STOP
         ENDIF
!
! Determine the polynomials associated with selected identification
!   number.
!
            NALL=3*NPARAMA+3*NPARAMB+NPARAMC
         IF(COCOUNT.GT.1) THEN
          DO K=NALL+1,NALL+COCOUNT-1
           NDFUNCS=NDFUNCS+NDONE(K)
          ENDDO
         ENDIF
!
          COCOUNT=COCOUNT+NALL
!
! Error Check.
!
            NALL=NALL+NPARAMC
         IF(COCOUNT.GT.NALL) THEN
           PRINT*,'COCOUNT=',COCOUNT
           PRINT*,'NALL   =',NALL
           PRINT*,'COCOUNT > NALL in SR CLEANISOP'
           STOP
         ENDIF
!
            GOTO 321
!
      ENDIF
!
!***************************************************
!  END C2
!***************************************************
!
            NDFUNCS=NDFUNCS+NDFUNCSC2
!
!***************************************************
!  C3) Upper Troposphere
!***************************************************
!
      IF(OH_PRESS.LT.PRESSES(2)) THEN
!
! Determine the identification number of the appropriate parameterization.
!
        CALL PARAMC(OH_LAT,OH_SEASON,COCOUNT,NSEAS,NLATS,ALATS)
!
! Error Check.
!
         IF(COCOUNT.EQ.0) THEN
           PRINT*,'COCOUNT=0 in SR CLEANISOP'
           STOP
         ENDIF
!
! Determine the polynomials associated with selected identification
!   number.
!
            NALL=3*NPARAMA+3*NPARAMB+2*NPARAMC
         IF(COCOUNT.GT.1) THEN
          DO K=NALL+1,NALL+COCOUNT-1
           NDFUNCS=NDFUNCS+NDONE(K)
          ENDDO
         ENDIF
!
          COCOUNT=COCOUNT+NALL
!
! Error Check.
!
            NALL=NALL+NPARAMC
         IF(COCOUNT.GT.NALL) THEN
           PRINT*,'COCOUNT=',COCOUNT
           PRINT*,'NALL   =',NALL
           PRINT*,'COCOUNT > NALL in SR CLEANISOP'
           STOP
         ENDIF
!
            GOTO 321
!
      ENDIF
!
!***************************************************
!  END C3
!***************************************************
!
            NDFUNCS=NDFUNCS+NDFUNCSC3
!
!***************************************************
! Error check.
!
      PRINT*,'********************************'
      PRINT*,'Stopped in SR CLEANISOP!'
      PRINT*,'No parameterization was chosen!'
      PRINT*,'********************************'
      STOP
!
!***************************************************
 321  RETURN

      END SUBROUTINE CLEANISOP
!**********************************************************************
      SUBROUTINE MODPOLLUTED(NDFUNCS)
!**********************************************************************
      IMPLICIT NONE

      INTEGER NDFUNCS,NALL,K
!**********************************************************************
! Select proper identification number of parameterization needed for
!   box characterized by unpolluted conditions with isoprene using
!   SR PARAMD-F. Store this value in COCOUNT. Then select the corresponding
!   polynomials and store in NDFUNCS.
!**********************************************************************
!   Created by Bryan Duncan.
!**********************************************************************
! D) Moderately Polluted
!***************************************************
!  D1) Surface
!***************************************************
!
      IF(OH_PRESS.GE.PRESSES(1)) THEN
!
! Determine the identification number of the appropriate parameterization.
!
        CALL PARAMD(OH_LAT,OH_SEASON,COCOUNT,NSEAS,NLATS,ALATS)
!
! Error Check.
!
         IF(COCOUNT.EQ.0) THEN
           PRINT*,'No parameterization was chosen!'
           PRINT*,'COCOUNT=0 in SR MODPOLLUTED'
           STOP
         ENDIF
!
! Determine the polynomials associated with selected identification
!   number.
!
            NALL=3*NPARAMA+3*NPARAMB+3*NPARAMC
         IF(COCOUNT.GT.1) THEN
          DO K=NALL+1,NALL+COCOUNT-1
           NDFUNCS=NDFUNCS+NDONE(K)
          ENDDO
         ENDIF
!
          COCOUNT=COCOUNT+NALL
!
! Error Check.
!
            NALL=NALL+NPARAMD
         IF(COCOUNT.GT.NALL) THEN
           PRINT*,'COCOUNT=',COCOUNT
           PRINT*,'NALL   =',NALL
           PRINT*,'COCOUNT > NALL in SR MODPOLLUTED'
           STOP
         ENDIF
!
            GOTO 321
!
      ENDIF
!
!***************************************************
!  END D1
!***************************************************
!  F1) Middle Troposphere
!***************************************************
!
      IF(OH_PRESS.LT.PRESSES(1).AND.OH_PRESS.GE.PRESSES(2)) THEN
!
! Determine the identification number of the appropriate parameterization.
!
        CALL PARAMF(OH_LAT,OH_SEASON,COCOUNT,NSEAS,NLATS,ALATS)
!
! Error Check.
!
         IF(COCOUNT.EQ.0) THEN
           PRINT*,'COCOUNT=0 in SR MODPOLLUTED'
           STOP
         ENDIF
!
! Determine the polynomials associated with selected identification
!   number.
!
            NALL=3*NPARAMA+3*NPARAMB+3*NPARAMC+NPARAMD+NPARAME
         IF(COCOUNT.GT.1) THEN
          DO K=NALL+1,NALL+COCOUNT-1
           NDFUNCS=NDFUNCS+NDONE(K)
          ENDDO
         ENDIF
!
          COCOUNT=COCOUNT+NALL
!
! Error Check.
!
            NALL=NALL+NPARAMF
         IF(COCOUNT.GT.NALL) THEN
           PRINT*,'COCOUNT=',COCOUNT
           PRINT*,'NALL   =',NALL
           PRINT*,'COCOUNT > NALL in SR MODPOLLUTED'
           STOP
         ENDIF
!
            GOTO 321
!
      ENDIF
!
!***************************************************
!  END F1
!***************************************************
! Error check.
!
      PRINT*,'********************************'
      PRINT*,'Stopped in SR MODPOLLUTED!'
      PRINT*,'No parameterization was chosen!'
      PRINT*,'********************************'
      STOP
!
!***************************************************
 321  RETURN

      END SUBROUTINE MODPOLLUTED

!**********************************************************************
      SUBROUTINE HIPOLLUTED(NDFUNCS)
!**********************************************************************
      IMPLICIT NONE

      INTEGER NDFUNCS,NALL,K
!**********************************************************************
! Select proper identification number of parameterization needed for
!   box characterized by highly polluted conditions with isoprene using
!   SR PARAME. Store this value in COCOUNT. Then select the corresponding
!   polynomials and store in NDFUNCS.
!**********************************************************************
!   Created by Bryan Duncan.
!**********************************************************************
! E) Highly Polluted
!***************************************************
!  E1) Surface
!***************************************************
!
      IF(OH_PRESS.GE.PRESSES(1)) THEN
!
! Determine the identification number of the appropriate parameterization.
!
        CALL PARAME(OH_LAT,OH_SEASON,COCOUNT,NSEAS,NLATS,ALATS)
!
! Error Check.
!
         IF(COCOUNT.EQ.0) THEN
           PRINT*,'COCOUNT=0 in SR HIPOLLUTEDISOP'
           STOP
         ENDIF
!
! Determine the polynomials associated with selected identification
!   number.
!
            NALL=3*NPARAMA+3*NPARAMB+3*NPARAMC+NPARAMD
         IF(COCOUNT.GT.1) THEN
          DO K=NALL+1,NALL+COCOUNT-1
           NDFUNCS=NDFUNCS+NDONE(K)
          ENDDO
         ENDIF
!
          COCOUNT=COCOUNT+NALL
!
! Error Check.
!
            NALL=NALL+NPARAME
         IF(COCOUNT.GT.NALL) THEN
           PRINT*,'COCOUNT=',COCOUNT
           PRINT*,'NALL   =',NALL
           PRINT*,'COCOUNT > NALL in SR HIPOLLUTEDISOP'
           STOP
         ENDIF
!
            GOTO 321
!
      ENDIF
!
!***************************************************
!  END E1
!***************************************************
! Sum of all NDONE functions in E).
!
            NDFUNCS=NDFUNCS+NDFUNCSE1
!
!***************************************************
! Error check.
!
      PRINT*,'********************************'
      PRINT*,'Stopped in SR HIPOLLUTEDISOP!'
      PRINT*,'No parameterization was chosen!'
      PRINT*,'********************************'
      STOP
!
!***************************************************
 321  RETURN

      END SUBROUTINE HIPOLLUTED

!**********************************************************************
      SUBROUTINE REMOTE(NDFUNCS)
!**********************************************************************
      IMPLICIT NONE

      INTEGER NDFUNCS,NALL,K
!**********************************************************************
! Select proper identification number of parameterization needed for 
!   box characterized by remote conditions and no isoprene using
!   SR PARAMB. Store this value in COCOUNT. Then select the corresponding 
!   polynomials and store in NDFUNCS.
!**********************************************************************
!   Created by Bryan Duncan.
!**********************************************************************
! B) Remote
!***************************************************
!  B1) Surface
!***************************************************
!
      IF(OH_PRESS.GE.PRESSES(1)) THEN
!
! Determine the identification number of the appropriate parameterization.
!
        CALL PARAMB(OH_LAT,OH_SEASON,COCOUNT,NSEAS,NLATS,ALATS)
!
! Error Check.
!
         IF(COCOUNT.EQ.0) THEN
           PRINT*,'COCOUNT=0!'
           PRINT*,'No parameterization was chosen!'
           PRINT*,'Stopped in SR REMOTE.'
           STOP 
         ENDIF
!
! Determine the polynomials associated with selected identification
!   number.
!
         IF(COCOUNT.GT.1) THEN
          DO K=1,COCOUNT-1
           NDFUNCS=NDFUNCS+NDONE(K)
          ENDDO
         ENDIF
!
! Error Check.
!
            NALL=NPARAMB
         IF(COCOUNT.GT.NALL) THEN
           PRINT*,'COCOUNT=',COCOUNT
           PRINT*,'NALL   =',NALL
           PRINT*,'COCOUNT > NALL in SR REMOTE'
           STOP
         ENDIF
!
            GOTO 321
!
      ENDIF
!
!***************************************************
!  END B1
!***************************************************
! Sum of all NDONE functions in B).
!
            NDFUNCS=NDFUNCS+NDFUNCSB1
!
!***************************************************
!  B2) Middle Troposphere
!***************************************************
!
     IF(OH_PRESS.LT.PRESSES(1).AND.OH_PRESS.GE.PRESSES(2)) THEN
!
! Determine the identification number of the appropriate parameterization.
!
        CALL PARAMB(OH_LAT,OH_SEASON,COCOUNT,NSEAS,NLATS,ALATS)
!
! Error Check.
!
         IF(COCOUNT.EQ.0) THEN
           PRINT*,'COCOUNT=0!'
           PRINT*,'No parameterization was chosen!'
           PRINT*,'Stopped in SR REMOTE.'
           STOP
         ENDIF
!
! Determine the polynomials associated with selected identification
!   number.
!
            NALL=NPARAMB
         IF(COCOUNT.GT.1) THEN
          DO K=NALL+1,NALL+COCOUNT-1
           NDFUNCS=NDFUNCS+NDONE(K)
          ENDDO
         ENDIF
!
          COCOUNT=COCOUNT+NALL
!
! Error Check.
!
            NALL=NALL+NPARAMB
         IF(COCOUNT.GT.NALL) THEN
           PRINT*,'COCOUNT=',COCOUNT
           PRINT*,'NALL   =',NALL
           PRINT*,'COCOUNT > NALL in SR REMOTE'
           STOP
         ENDIF
!
            GOTO 321
!
      ENDIF
!
!***************************************************
!  END B2
!***************************************************
! Sum of all NDONE functions in B).
!
            NDFUNCS=NDFUNCS+NDFUNCSB2
!
!***************************************************
!  B3) Upper Troposphere
!***************************************************
!
      IF(OH_PRESS.LT.PRESSES(2)) THEN
!
! Determine the identification number of the appropriate parameterization.
!
        CALL PARAMB(OH_LAT,OH_SEASON,COCOUNT,NSEAS,NLATS,ALATS)
!
! Error Check.
!
         IF(COCOUNT.EQ.0) THEN
           PRINT*,'COCOUNT=0!'
           PRINT*,'No parameterization was chosen!'
           PRINT*,'Stopped in SR REMOTE.'
           STOP
         ENDIF
!
! Determine the polynomials associated with selected identification
!   number.
!
            NALL=2*NPARAMB
         IF(COCOUNT.GT.1) THEN
          DO K=NALL+1,NALL+COCOUNT-1
           NDFUNCS=NDFUNCS+NDONE(K)
          ENDDO
         ENDIF
!
          COCOUNT=COCOUNT+NALL
!
! Error Check.
!
            NALL=NALL+NPARAMB
         IF(COCOUNT.GT.NALL) THEN
           PRINT*,'COCOUNT=',COCOUNT
           PRINT*,'NALL   =',NALL
           PRINT*,'COCOUNT > NALL in SR REMOTE'
           STOP
         ENDIF
!
            GOTO 321
!
      ENDIF
!
!***************************************************
!  END B3
!***************************************************
! Sum of all NDONE functions in B).
!
            NDFUNCS=NDFUNCS+NDFUNCSB3
!
!***************************************************
! Error check.
!    
      PRINT*,'********************************'
      PRINT*,'Stopped in SR REMOTE!'
      PRINT*,'No parameterization was chosen!'
      PRINT*,'********************************'
      STOP
!
!***************************************************
 321  RETURN
      END SUBROUTINE REMOTE

  END MODULE identification_mod
