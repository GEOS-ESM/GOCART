#include "MAPL_Generic.h"
! **********************************************************************
      SUBROUTINE CALC_OH(NDFUNCS,OOB)
! **********************************************************************

! USES: 
      USE MAPL

      USE cblock_OH_mod    ! "CMN_OH"

      IMPLICIT NONE

      INTEGER MM,M,KDONE,ID,JKEEP,NOROWS,KKRW,KROW,JVAR, &
     &       K,NDFUNCS,NOTDONE,OOB,OOA
      REAL*8 FNCVAL,SYMBOL,BMINSA,MIDWAY,C(MXCOL),X(MXVAR),ERRORMINMAX, &
     &       DIFFERR
      CHARACTER*18 ERRORCHECK(5),VARNAMES(19)
      DATA ERRORCHECK/'DECLINATION ANGLE','NOx','PRESSURE', &
     &                'LATITUDE','O3'/
      DATA VARNAMES/'DECLINATION ANGLE','O3 COLUMN','NOx','OZONE', &
     &   'CO','CH4','H2O','ACETONE','PRESSURE','PROPENE','ETHANE', &
     &   'PROPANE','ALK4','ALBEDO','LATITUDE','TEMPERATURE', &
     &   'CLOUD ABOVE','CLOUD BELOW','ISOPRENE'/
! **********************************************************************
! Created by Bryan Duncan.
! **********************************************************************
! Calculate the parameterized OH.
! **********************************************************************
       IF (0) THEN
           print*,'*************'
           print*,'inside SR calc_oh'
           print*,'NPARAM=',NPARAM
           print*,'COCOUNT=',COCOUNT
           print*,'NDONE=',NDONE(COCOUNT)
           print*,'*************'	
       ENDIF

       IF(COCOUNT.GT.NPARAM) THEN
            PRINT*,'*************'
            PRINT*,'COCOUNT=',COCOUNT
            PRINT*,'NPARAM=',NPARAM
            PRINT*,'SR CALC_OH: COCOUNT.GT.NPARAM!'
            PRINT*,'*************'
            STOP
       ENDIF
!
! NOTDONE = error check.
!         = 0 then done function not found.
!         = 1 then done function found.
!
       NOTDONE=0

       IF(COCOUNT.LT.1) THEN
           print*,'*************'
           print*,'STOPPED IN SR CALC_OH!'
           print*,'The box lies nowhere!' 
           print*,'NPARAM=',NPARAM
           print*,'COCOUNT=',COCOUNT
           print*,'NDONE=',NDONE(COCOUNT)
           print*,'*************'
           STOP
       ENDIF
!
! **********************************************************************
! 1) Independent variables are rescaled to fit range [-0.9, 0.9].  The 
!    parameterization was constructed with the independent variables rescaled.
! **********************************************************************
! Loop through MOVAR independent variables.
!
       DO M=1,MOVAR(COCOUNT)
!        print*,RANGEM(COCOUNT,M,1),INDVAR(M),RANGEM(COCOUNT,M,2)
!
        BMINSA=RANGEM(COCOUNT,M,2)-RANGEM(COCOUNT,M,1)
!
! If point falls out of prescribed ranges, reset to just in range.
!
        IF(INDVAR(M).LE.RANGEM(COCOUNT,M,1)) THEN
!
         IF(INDVAR(M).LT.RANGEM(COCOUNT,M,1)) THEN 
!
!************************
! Error Check.
!
! The following four independent variables should not be out of bounds
!  of the prescribed ranges.  An error here is a good indication that
!  the wrong parameterization has been chosen by the code.
!   M=1 : Declination Angle (i.e., Season)
!   M=3 : NOx
!   M=9 : Pressure (i.e., Altitude)
!   M=15: Latitude
!  
         IF(M.EQ.1.OR.M.EQ.15.OR.M.EQ.3.OR.M.EQ.9)THEN
           PRINT*,'*******************'
           PRINT*,'RANGE OUT OF BOUNDS FOR:'
           IF(M.EQ.1)  PRINT*,ERRORCHECK(1)
           IF(M.EQ.3)  PRINT*,ERRORCHECK(2)
           IF(M.EQ.9)  PRINT*,ERRORCHECK(3)
           IF(M.EQ.15) PRINT*,ERRORCHECK(4)
           IF(M.EQ.4)  PRINT*,ERRORCHECK(5)
           PRINT*,''
           PRINT*,'Minimum Range = ',RANGEM(COCOUNT,M,1) 
           PRINT*,'Maximum Range = ',RANGEM(COCOUNT,M,2) 
           PRINT*,'Variable Value= ',INDVAR(M)
           PRINT*,''
           PRINT*,'This variable should not be out of bounds'
           PRINT*,'of the prescribed ranges.  An error here is a good'
           PRINT*,'indication that the wrong parameterization has been'
           PRINT*,'chosen by the code.'
           PRINT*,''
           PRINT*,'STOPPED IN SR CALC_OH!'
           PRINT*,'*******************'
           STOP
         ENDIF
!
! End Error Check.
!************************
!
         ENDIF
!
             DIFFERR=RANGEM(COCOUNT,M,1)-INDVAR(M)
!
              OOA=M 
           IF(OOB.EQ.3.AND.M.EQ.18) OOA=19
!
             ERRORMINMAX=INDVAR(M)
             INDVAR(M)=RANGEM(COCOUNT,M,1)+1D-2*BMINSA
!
         IF(ERRORON.EQ.0.AND.DIFFERR.GT.1.D-8) THEN
           PRINT*,'*******************'
           PRINT*,'WARNING! RANGE OUT OF BOUNDS FOR:  ',VARNAMES(OOA)
           PRINT*,'Variable reset to value inside prescribed range.' 
           PRINT*,'Minimum Range        = ',RANGEM(COCOUNT,M,1) 
           PRINT*,'Maximum Range        = ',RANGEM(COCOUNT,M,2) 
           PRINT*,'Variable Value - NEW = ',INDVAR(M)
           PRINT*,'Variable Value - OLD = ',ERRORMINMAX
           PRINT*,'Latitude = ', INDVAR(15)
           PRINT*,'Longitude= ', INDVAR(20)
           PRINT*,'Pressure = ', INDVAR(9)
           PRINT*,'Dec. >   = ', INDVAR(1)
           PRINT*,'NOx      = ', INDVAR(3)

!      if(INDVAR(20).GT.-135.and.INDVAR(20).lt.-90.and.INDVAR(15).
!     *  lt.60.and.INDVAR(15).gt.40.) print*,'bryan'

           PRINT*,'*******************'
         ENDIF
!
        ENDIF
!
        IF(INDVAR(M).GE.RANGEM(COCOUNT,M,2)) THEN
!
!************************
!
         IF(INDVAR(M).GT.RANGEM(COCOUNT,M,2)) THEN
!
!************************
! Error Check.
!
         IF(M.EQ.1.OR.M.EQ.15.OR.M.EQ.3.OR.M.EQ.9)THEN
           PRINT*,'*******************'
           PRINT*,'RANGE OUT OF BOUNDS FOR:'
           IF(M.EQ.1)  PRINT*,ERRORCHECK(1)
           IF(M.EQ.3)  PRINT*,ERRORCHECK(2)
           IF(M.EQ.9)  PRINT*,ERRORCHECK(3)
           IF(M.EQ.15) PRINT*,ERRORCHECK(4)
           IF(M.EQ.4)  PRINT*,ERRORCHECK(5)
           PRINT*,''
           PRINT*,'Minimum Range = ',RANGEM(COCOUNT,M,1) 
           PRINT*,'Maximum Range = ',RANGEM(COCOUNT,M,2)
           PRINT*,'Variable Value= ',INDVAR(M)
           PRINT*,''
           PRINT*,'This variable should not be out of bounds'
           PRINT*,'of the prescribed ranges.  An error here is a good'
           PRINT*,'indication that the wrong parameterization has been'
           PRINT*,'chosen by the code.'
           PRINT*,''
           PRINT*,'STOPPED IN SR CALC_OH!'
           PRINT*,'*******************'
           STOP
         ENDIF
!************************
!
         ENDIF
!
             DIFFERR=INDVAR(M)-RANGEM(COCOUNT,M,2)
!
              OOA=M 
           IF(OOB.EQ.3.AND.M.EQ.18) OOA=19
!
             ERRORMINMAX=INDVAR(M)
             INDVAR(M)=RANGEM(COCOUNT,M,2)-1D-2*BMINSA
!
         IF(ERRORON.EQ.0.AND.DIFFERR.GT.1.D-8) THEN
           PRINT*,'*******************'
           PRINT*,'WARNING! RANGE OUT OF BOUNDS FOR:  ',VARNAMES(OOA)
           PRINT*,'Variable reset to value just prescribed range.' 
           PRINT*,'Minimum Range        = ',RANGEM(COCOUNT,M,1) 
           PRINT*,'Maximum Range        = ',RANGEM(COCOUNT,M,2) 
           PRINT*,'Variable Value - NEW = ',INDVAR(M)
           PRINT*,'Variable Value - OLD = ',ERRORMINMAX
           PRINT*,'Latitude = ', INDVAR(15)
           PRINT*,'Longitude= ', INDVAR(20)
           PRINT*,'Pressure = ', INDVAR(9)
           PRINT*,'Dec. >   = ', INDVAR(1)
           PRINT*,'NOx      = ', INDVAR(3)
!      if(INDVAR(20).GT.-135.and.INDVAR(20).lt.-90.and.INDVAR(15). &
!     *  lt.60.and.INDVAR(15).gt.40.) print*,'bryan'
           PRINT*,'*******************'
         ENDIF
!
        ENDIF
!
        MIDWAY=0.9
!
         X(M)=-MIDWAY + (((INDVAR(M)-RANGEM(COCOUNT,M,1))/BMINSA)*MIDWAY*2.)
!
      ENDDO
!
! ***********************************************************
! ***********************************************************
!
       DO 2 KDONE=1,NDONE(COCOUNT)
!
! ***********************************************************
! 2) See if point within range of done function.
! ***********************************************************
!
	IF (0) PRINT *, 'KDONE NDONE(COUNT) ' , KDONE, NDONE(COCOUNT)
          ID=IDENTOLD(COCOUNT,KDONE)
!
! *** Iterate number of inequalites (rows) this element.
!
          JKEEP=0
         NOROWS=IROWST(COCOUNT,ID,0)
!
         DO 30 KROW=1,NOROWS
            KKRW=IROWST(COCOUNT,ID,KROW)
!
! *** calculate plane equation values.
!
            FNCVAL=0.-ELTODO(COCOUNT,KKRW,MOVAR(COCOUNT)+2)

       DO 20 JVAR=1,MOVAR(COCOUNT)
           FNCVAL=FNCVAL+(ELTODO(COCOUNT,KKRW,JVAR)*INDVAR(JVAR))
   20  CONTINUE
!
! *** IF points above the "plane" are sought, verify that this point lies
! *** on or above plane.
!
            SYMBOL=ELTODO(COCOUNT,KKRW,MOVAR(COCOUNT)+1)
            IF (SYMBOL.GT.0..AND.FNCVAL.LT.0.) JKEEP=1
	    
! *** If points below the "plane" are sought, verify that this point lies
! *** below plane.
!
            IF (SYMBOL.LT.0..AND.FNCVAL.GE.0.) JKEEP=1
!            IF (SYMBOL.LT.0..AND.FNCVAL.GT.0.) JKEEP=1
!	    IF (MAPL_AM_I_ROOT() .AND. JKEEP.EQ.1) THEN 
!	       PRINT *, 'CALCOH: points not on plane NOTDONE : ', NOTDONE
!	    ENDIF   
            IF(JKEEP.EQ.1) GOTO 2

   30    CONTINUE
!
! ***********************************************************
! 3) Calculate parameterized OH.
! ***********************************************************
!
         IF(JKEEP.EQ.0) THEN 

          NOTDONE=1
!
         do 3 MM=1,NENDW(COCOUNT,KDONE)
	   C(MM)=0.
           C(MM)=COEFF(COCOUNT,KDONE,MM)
 3       CONTINUE
!
           NDFUNCS=NDFUNCS+KDONE
!

 
           CALL FUNC_IFS(NDFUNCS,X,C,PARAMOH)

!
! 
!          GOTO 2
          GOTO 33
!
         ENDIF
!
! ***********************************************************
 2      CONTINUE
! ***********************************************************
 33     CONTINUE
! Error Check
      IF(NOTDONE.EQ.0) THEN
        PRINT*,'SR CALC_OH:  Done function not found!!!!'
        STOP
      ENDIF
!
      RETURN

      END SUBROUTINE CALC_OH
