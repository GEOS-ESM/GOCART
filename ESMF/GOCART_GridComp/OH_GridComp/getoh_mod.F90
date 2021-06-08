#include "MAPL_Generic.h"

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  getoh_mod ---
!
! !INTERFACE:
!

   MODULE  getoh_mod

! !USES:

   USE ESMF
   USE MAPL
   USE identification_mod   ! clean, cleanisop, remote, hipolluted, modpolluted
   USE cblock_OH_mod        ! "CMN_OH"


   IMPLICIT NONE
!
! !PUBLIC MEMBER FUNCTIONS:
!
   PRIVATE
   PUBLIC  :: getoh 

!*****************************************************************************
   CONTAINS
!*****************************************************************************
      SUBROUTINE GETOH(COLON,COLAT,COALT,FIRSTDT,BOH,OHLMN,DOskip, &
                       INDVARA, INDVARB, INDVARC, INDVARD, &
                       iLon, iLat, iVert,OHparamNumOut)
!*****************************************************************************
      INTEGER :: iLon, iLat, iVert
      INTEGER :: COLON,COLAT,COALT,DOskip
      REAL*8 :: BOH(iLon, iLat, iVert)
      REAL*8 :: INDVARA(iLon, iLat, iVert,MXVAR)
      REAL*8 :: INDVARB(iLon, iLat, iVert,MXVAR)
      REAL*8 :: INDVARC(iLon, iLat, iVert,MXVAR)
      REAL*8 :: INDVARD(iLon, iLat, iVert,MXVAR)
!
      REAL*8,  intent(inOut) :: OHparamNumOut      ! number of    parametrization used

      INTEGER FIRSTDT
      INTEGER I,J,K,L,M
      INTEGER NDFUNCS,NALL,NLOOK,OHLMN
      INTEGER OOB
!
!*****************************************************************************
! Created by Bryan Duncan.
!*****************************************************************************
! This SR determines which parameterization is needed for an individual 
! box and calculates the parameterized [OH] (SR CALC_OH).
!
! The parameterizations are defined by latitude, pressure, [NOx], [ISOP]
! and season.  (The corresponding model variables are OH_LAT, OH_PRESS,
! OH_NOX, OH_ISOP and OH_SEASON, resp.)

      OH_LAT   = INDVARA(COLON,COLAT,COALT,15)
      OH_LON   = INDVARA(COLON,COLAT,COALT,20)
      OH_PRESS = INDVARA(COLON,COLAT,COALT,9)
      OH_NOX   = INDVARA(COLON,COLAT,COALT,3)
      OH_ISOP  = INDVARC(COLON,COLAT,COALT,18)
      OH_O3    = INDVARA(COLON,COLAT,COALT,4)
!
!*****************************************************************************
! List of Variables & Arrays
!
! BOH = array holding parameterized OH.
!
! Dimensions of troposphere.  User specified values in CMN_OH.
!   COALT = position of box in vertical in troposphere
!   COLAT = position of box from pole to pole (following one longitude band)
!   COLON = position of box circling globe (following one latitude band)
!
! NDFUNCS = number of polynomials.  NDFUNCS is greater than the number
!           of parameterizations (NPARAM). A parameterization of
!           a subdomain may be described by more than one polynomial
!           due to domain divisions.
!
      NDFUNCS=0
!
! FIRSTDT = 0 on first chemistry time step.
!         = 1 on any time chemistry step, but first. 
!
! PARAMOH = parameterized [OH].
!
      PARAMOH=100.
!
! COCOUNT = identification number of parameterization.
!
      COCOUNT=0
!
! OOB     = bookkeeping for error check in SR CALC_OH.
!
      OOB=0
!
! INDVAR  = independent variables needed to calculate parameterized [OH].
!
      INDVAR(1:MXVAR)=0.
!
! DOskip  = Counter to skip parameterization if necessary.
!         = 0 : use parameterization to get [OH].
!         = 1 : set to low value and skip parameterization.
!         = 2 : use climatological [OH] (see SR READAVGOH & SR INTERPOH).
!         = 3 : use [OH] = 1x10^5 if isoprene > 0.6 ppbv and skip parameterization
!         = 4 : use [OH] = 1x10^4 if isoprene > 3 ppbv and skip parameterization
!
! Variables specific to box in question:
!
! OH_NOX    = NOx concentration
! OH_PRESS  = pressure
! OH_LAT    = latitude
! OH_ISOP   = isoprene concentration
! OH_LON    = longitude
! OH_MONTH = month (January-December = 1-12)
! OH_SEASON = (Northern Hemispheric) season
!        1 --> winter (Dec, Jan, Feb)
!        2 --> spring (Mar, Apr, May)
!        3 --> summer (Jun, Jul, Aug)
!        4 --> autumn (Sep, Oct, Nov)
!
!*****************************************************************************
! 1)  Determine which parameterization the box needs.
!*****************************************************************************
! Error Check.
! If you'd like to see the specifics for the model box you're
! calculating OH, change value of NLOOK.
!
      NLOOK=0  ! Set to 1 for DEBUG 
      IF (NLOOK.EQ.1) THEN
	   IF (MAPL_AM_I_ROOT()) THEN
            print*,''
            print*,'OHparam (ppbv/pptv units)'
            print*,'OH_LAT   =',OH_LAT
            print*,'OH_LON   =',OH_LON
	    print*,'OH_DEC   =',INDVARA(COLON,COLAT,COALT,1)
            print*,'OH_PRESS =',OH_PRESS
            print*,'OH_NOX   =',OH_NOX
            print*,'OH_ISOP  =',OH_ISOP
            print*,'OH_O3    =',OH_O3
            print*,'OH_O3COL =',INDVARA(COLON,COLAT,COALT,2)
	    print*,'OH_CO =',INDVARA(COLON,COLAT,COALT,5)
	    print*,'OH_CH4 =',INDVARA(COLON,COLAT,COALT,6)
	    print*,'OH_H2O =',INDVARA(COLON,COLAT,COALT,7)
	    print*,'OH_ACET =',INDVARA(COLON,COLAT,COALT,8)
	    print*,'OH_PRPE =',INDVARA(COLON,COLAT,COALT,10)
	    print*,'OH_ETHA =',INDVARA(COLON,COLAT,COALT,11)	    
	    print*,'OH_PROP =',INDVARA(COLON,COLAT,COALT,12)
	    print*,'OH_ALK4 =',INDVARA(COLON,COLAT,COALT,13) 
	    print*,'OH_ALB =',INDVARA(COLON,COLAT,COALT,14)  
	    print*,'OH_TEMP =',INDVARA(COLON,COLAT,COALT,16)  	    
	    print*,'OH_CLOUDA =',INDVARA(COLON,COLAT,COALT,17)
	    print*,'OH_CLOUDB =',INDVARD(COLON,COLAT,COALT,18)	    
	      	    	
    	    print*,'OH_SEASON=',OH_SEASON
            print*,'COCOUNT =',COCOUNT
	 ENDIF
      ENDIF
!
!*****************************************************************************
! SR SKIP prevents calling OH parameterizations in regions where OH is very low,
! resets the values of isoprene and NOx to keep from calling nonexistent
! parameterizations, and performs other miscellaneous functions.
!

      CALL SKIP(DOskip,COLON,COLAT,COALT, &
                INDVARA, INDVARB, INDVARC, INDVARD, iLon, iLat, iVert)
     

!
! If DOskip = 1, set [OH] to low value (i.e., skip parameterization).
!
      IF (DOskip.EQ.1) GOTO 320
!
! In subdomains where OH is very low or negligible due to low
!  sunlight (e.g., high latitudes in winter), concentrations of OH are
!  set to climatological mean values as a function of latitude,
!  altitude and season. No parameterized OH is calculated. 
!
      IF (OH_O3.LT.1.) DOskip=2
      IF (DOskip.EQ.2) THEN
         CALL INTERPOH()
         GOTO 320
      ENDIF
!
!*************
! Error Check.
! Surface
!       IF(OH_PRESS.LT.PRESSES(1)) goto 320
!
! MT
!       IF(OH_PRESS.GE.PRESSES(1)) goto 320
!       IF(OH_PRESS.LT.PRESSES(2)) goto 320
!
! UT
!       IF(OH_PRESS.GE.PRESSES(2)) goto 320
!
!
!*****************************************************************************
! B) REMOTE:  Unpolluted Domain:  Low NOx, No Isoprene (<10 ppt)
!             Surface Layer, MT, & UT (NOx <50 ppt)
!*****************************************************************************
! Surface & MT
      IF (OH_ISOP.LE.10..AND.OH_O3.LE.45.AND. &
     &    OH_NOX.LE.CNOXS(7).AND.OH_PRESS.GE.PRESSES(2)) THEN
!     
         OHparamNumOut = 1. 
         CALL REMOTE(NDFUNCS)
!
!    Assign appropriate independent variables to be used in
!     OH calculation.
!
         IF (OH_PRESS.GE.PRESSES(1)) THEN
            INDVAR(1:MXVAR)=INDVARA(COLON,COLAT,COALT,1:MXVAR)
            OOB=1
         ELSE
            INDVAR(1:MXVAR)=INDVARB(COLON,COLAT,COALT,1:MXVAR)
            OOB=2
         ENDIF
!
         GOTO 321
!
      ENDIF
! UT
      IF (OH_ISOP.LE.10..AND.OH_O3.LE.100..AND. &
     &    OH_NOX.LE.CNOXS(6).AND.OH_PRESS.LT.PRESSES(2)) THEN
!
!        IF (MAPL_AM_I_ROOT()) print*, ' REMOTE PARAM '
	 OHparamNumOut = 1.
         CALL REMOTE(NDFUNCS)
!
!    Assign appropriate independent variables to be used in
!     OH calculation.
!
         INDVAR(1:MXVAR)=INDVARB(COLON,COLAT,COALT,1:MXVAR)
         OOB=2
!
         GOTO 321
!
      ENDIF
!
!*****************************************************************************
! END B) REMOTE
!*****************************************************************************
!
      NDFUNCS=NDFUNCS+NDFUNCSB1+NDFUNCSB2+NDFUNCSB3

!
!*****************************************************************************
! A) CLEAN:  Unpolluted Domain:  Low NOx, No Isoprene (<10 ppt) 
!            Surface Layer (NOx <300 ppt), Middle & Up Tropospheres
!*****************************************************************************
!
      IF (OH_ISOP.LE.10.) THEN
         IF((OH_NOX.LE.CNOXS(1).AND.OH_PRESS.GE.PRESSES(2)).OR.   &
     &        (OH_NOX.LE.CNOXS(2).AND.OH_PRESS.LT.PRESSES(2))) THEN
!     
!           IF (MAPL_AM_I_ROOT()) PRINT*, 'CLEAN PARAM'
	    OHparamNumOut = 2.
            CALL CLEAN(NDFUNCS)

!
!    Assign appropriate independent variables to be used in
!     OH calculation.
!
            IF (OH_PRESS.GE.PRESSES(1)) THEN
               INDVAR(1:MXVAR)=INDVARA(COLON,COLAT,COALT,1:MXVAR)
               OOB=1
            ELSE
               INDVAR(1:MXVAR)=INDVARB(COLON,COLAT,COALT,1:MXVAR)
               OOB=2
            ENDIF
!
            GOTO 321
!     
         ENDIF
!
      ENDIF
!
!*****************************************************************************
! END A) CLEAN
!*****************************************************************************
!
      NDFUNCS=NDFUNCS+NDFUNCSA1+NDFUNCSA2+NDFUNCSA3
!
!*****************************************************************************
! C) CLEAN with ISOPRENE:  Low NOx, Isoprene (>10 ppt) 
!            Surface Layer (NOx <300 ppt), Mid Troposphere 
!            (NOx <300 ppt), Up Troposphere (NOx <500ppt)
!*****************************************************************************
!
      IF (OH_ISOP.GT.10.) THEN
!
         IF ((OH_NOX.LE.CNOXS(1).AND.OH_PRESS.GE.PRESSES(2)).OR.   &
     &        (OH_NOX.LE.CNOXS(2).AND.OH_PRESS.LT.PRESSES(2))) THEN
!
            IF (OH_PRESS.GE.PRESSES(1)) THEN
               IF (DOskip.EQ.3) THEN
                  BOH(COLON,COLAT,COALT)=1.E5
                  GOTO 322
               ENDIF
!
               IF (DOskip.EQ.4) THEN
                  BOH(COLON,COLAT,COALT)=1.E4
                  GOTO 322
               ENDIF
            ENDIF
!
!           IF (MAPL_AM_I_ROOT()) PRINT*, 'CLEANISOP PARAM'
	    OHparamNumOut = 3.
            CALL CLEANISOP(NDFUNCS)
!
            IF (OH_PRESS.GE.PRESSES(1)) THEN
               INDVAR(1:MXVAR)=INDVARC(COLON,COLAT,COALT,1:MXVAR)
               OOB=3
            ELSE
               INDVAR(1:MXVAR)=INDVARD(COLON,COLAT,COALT,1:MXVAR)
               OOB=4
            ENDIF
!     
            GOTO 321
!     
         ENDIF
!     
      ENDIF
!
!*****************************************************************************
! END C) CLEAN with ISOPRENE
!*****************************************************************************
!
      NDFUNCS=NDFUNCS+NDFUNCSC1+NDFUNCSC2+NDFUNCSC3
!
!*****************************************************************************
! D) MIDNOx SUBDOMAIN with ISOPRENE: 
!            Mid NOx (>300 ppt <1000), Isoprene
!            Surface (see F for Middle Troposphere)
!*****************************************************************************
!
      IF (OH_NOX.GT.CNOXS(1).AND.   &
     &     OH_NOX.LE.CNOXS(3).AND.OH_PRESS.GE.PRESSES(1)) THEN
!
!        IF (MAPL_AM_I_ROOT()) PRINT*, 'MIDNOx with ISOP PARAM'
	 OHparamNumOut = 5.
         CALL MODPOLLUTED(NDFUNCS)
!
         INDVAR(1:MXVAR)=INDVARC(COLON,COLAT,COALT,1:MXVAR)
         OOB=3
!
         GOTO 321
!
      ENDIF
!     
!*****************************************************************************
! END D) MIDNOx SUBDOMAIN with ISOPRENE - Surface
!*****************************************************************************
!
      NDFUNCS=NDFUNCS+NDFUNCSD1
!
!*****************************************************************************
! E) HINOx SUBDOMAIN with ISOPRENE
!               Hi NOx (>1000), Isoprene
!*****************************************************************************
!
      IF (OH_NOX.GT.CNOXS(3).AND.OH_PRESS.GE.PRESSES(1)) THEN
!
!         IF (MAPL_AM_I_ROOT()) PRINT*, 'Hi NOX with ISOP PARAM'
	  OHparamNumOut = 6.
         CALL HIPOLLUTED(NDFUNCS)
!
         INDVAR(1:MXVAR)=INDVARC(COLON,COLAT,COALT,1:MXVAR)
         OOB=3
!
         GOTO 321
!     
      ENDIF
!
!*****************************************************************************
! END E) HINOx SUBDOMAIN with ISOPRENE
!*****************************************************************************
!
      NDFUNCS=NDFUNCS+NDFUNCSE1
!
!*****************************************************************************
! F) MIDNOx SUBDOMAIN with ISOPRENE
!            Mid NOx (>300 ppt <1000), Isoprene
!            Middle Troposphere 
!****************************************************************************
!     
      IF(OH_NOX.GT.CNOXS(1).AND.OH_NOX.LE.CNOXS(3).AND. &
     &   OH_PRESS.LT.PRESSES(1).AND.OH_PRESS.GE.PRESSES(2)) THEN
!     
!        IF (MAPL_AM_I_ROOT()) PRINT*, 'MIDNOx with ISOPRENE '
	 OHparamNumOut = 4.
         CALL MODPOLLUTED(NDFUNCS)
!     
         INDVAR(1:MXVAR)=INDVARD(COLON,COLAT,COALT,1:MXVAR)
         OOB=4
!     
         GOTO 321
!
      ENDIF
!
!*****************************************************************************
! END F) MIDNOx SUBDOMAIN with ISOPRENE - MT
!*****************************************************************************
!
      NDFUNCS=NDFUNCS+NDFUNCSF1
!
!*****************************************************************************
! Error Check.
!
      IF (MAPL_AM_I_ROOT()) THEN
         PRINT*,'BOX LIES NOWHERE!'
         PRINT*,'STOPPED IN SR GETOH!'
         PRINT*,'*********************'
         PRINT*,'OH_LAT=',OH_LAT
         PRINT*,'OH_PRESS=',OH_PRESS
         PRINT*,'OH_NOX=',OH_NOX
         PRINT*,'OH_ISOP=',OH_ISOP
         PRINT*,'OH_SEASON=',OH_SEASON
         PRINT*,'*********************'
      END IF
      STOP
!     
!*****************************************************************************
!*****************************************************************************
! 2)  Calculate the parameterized OH.
!*****************************************************************************
!*****************************************************************************
!
 321  CALL CALC_OH(NDFUNCS,OOB)

 320  BOH(COLON,COLAT,COALT)=PARAMOH

!bnd
      IF(DOskip.EQ.2) GOTO 322
!bnd
!
!***********************************************************
! SR correctOH: Correct aerosol problem from chem1d.
!***********************************************************
      IF(BOH(COLON,COLAT,COALT).GT.0.) THEN
  !sas       IF (MAPL_AM_I_ROOT()) PRINT *, 'GETOH_MOD: BOH before correctOH: ' , BOH(COLON,COLAT,COALT)
         CALL correctOH(BOH,COLON,COLAT,COALT, iLon, iLat, iVert)
  !sas    	 IF (MAPL_AM_I_ROOT()) PRINT *, 'BOH after correctOH: ' , BOH(COLON,COLAT,COALT)
      ENDIF

!bnd      GOTO 330

      
! Correct problem in the tropical UT
!  SH
       IF(OH_LAT.LE.ALATS(4).AND.OH_LAT.GT.-50.) THEN
         IF(OH_PRESS.LT.130) THEN
              BOH(COLON,COLAT,COALT)=BOH(COLON,COLAT,COALT)*valSH(OHLMN*4)
           GOTO 322
         ENDIF
         IF(OH_PRESS.LT.150) THEN
              BOH(COLON,COLAT,COALT)=BOH(COLON,COLAT,COALT)*valSH(OHLMN*4-1)
           GOTO 322
         ENDIF
         IF(OH_PRESS.LT.180) THEN
              BOH(COLON,COLAT,COALT)=BOH(COLON,COLAT,COALT)*valSH(OHLMN*4-2)
           GOTO 322
         ENDIF
         IF(OH_PRESS.LT.220) THEN
              BOH(COLON,COLAT,COALT)=BOH(COLON,COLAT,COALT)*valSH(OHLMN*4-3)
           GOTO 322
         ENDIF

       ENDIF
!  NH
       IF(OH_LAT.LE.50.AND.OH_LAT.GT.ALATS(4)) THEN
         IF(OH_PRESS.LT.130) THEN
              BOH(COLON,COLAT,COALT)=BOH(COLON,COLAT,COALT)*valNH(OHLMN*4)
           GOTO 322
         ENDIF
         IF(OH_PRESS.LT.150) THEN
              BOH(COLON,COLAT,COALT)=BOH(COLON,COLAT,COALT)*valNH(OHLMN*4-1)
           GOTO 322
         ENDIF
         IF(OH_PRESS.LT.180) THEN
              BOH(COLON,COLAT,COALT)=BOH(COLON,COLAT,COALT)*valNH(OHLMN*4-2)
           GOTO 322
         ENDIF
         IF(OH_PRESS.LT.220) THEN
              BOH(COLON,COLAT,COALT)=BOH(COLON,COLAT,COALT)*valNH(OHLMN*4-3)
           GOTO 322
         ENDIF
       ENDIF


! decrease OH by 20%
! this allows feedback of CO - amplification!!!!!!!
!       BOH(COLON,COLAT,COALT)=BOH(COLON,COLAT,COALT)*0.8
       
!***********************************************************
 322   CONTINUE
      !sas -comment out this part
!bnd
!!$       IF(OH_PRESS.LT.300.AND.BOH(COLON,COLAT,COALT).GT.1.5e6) THEN
!!$          BOH(COLON,COLAT,COALT)=BOH(COLON,COLAT,COALT)*0.33
!!$       ENDIF
!bnd
 330    CONTINUE        
!

       END SUBROUTINE GETOH
!**********************************************************************
!**********************************************************************
      SUBROUTINE SKIP(DOskip,COLON,COLAT,COALT, &
                      INDVARA, INDVARB, INDVARC, INDVARD, &
                      iLon, iLat, iVert)
!**********************************************************************
      INTEGER :: iLon, iLat, iVert
      REAL*8 :: INDVARA(iLon, iLat, iVert, MXVAR)
      REAL*8 :: INDVARB(iLon, iLat, iVert, MXVAR)
      REAL*8 :: INDVARC(iLon, iLat, iVert, MXVAR)
      REAL*8 :: INDVARD(iLon, iLat, iVert, MXVAR)
      INTEGER DOskip,COLON,COLAT,COALT,NODOMOD,NODOHI
!**********************************************************************
!   Created by Bryan Duncan.
!**********************************************************************
! SR SKIP performs the following functions:
! 1) Low OH :  Use if statements to keep calling OH parameterization
!              when OH is very low.  About 1/5 of the boxes
!              are skipped here.
! 2) RESET  :  Reset the values of isoprene and NOx to keep from calling
!              "nonexistent" parameterizations.
! 3) Other Miscellaneous Functions
!
!**********************************************************************
! List of Variables & Arrays
!
! DOskip = Counter to skip parameterization. 
!        = 0 : use parameterization to get [OH]
!        = 1 : set to low value and skip parameterization
!        = 2 : use average [OH] (see SR readavgOH)
!        = 3 : use [OH] = 1x10^5 if isoprene > 600 pptv and skip parameterization
!        = 4 : use [OH] = 1x10^4 if isoprene > 3000 pptv and skip parameterization
!
            DOskip=0
!
! Dimensions of troposphere.  User specified values in CMN_OH.
!   COALT = position of box in vertical in troposphere
!   COLAT = position of box from pole to pole (following one longitude band)
!   COLON = position of box circling globe (following one latitude band)
!
! PARAMOH = parameterized [OH].
!
! Variables specific to box in question:
!
! OH_NOX    = NOx concentration
! OH_PRESS  = pressure
! OH_LAT    = latitude
! OH_ISOP   = isoprene concentration
! OH_LON    = longitude
! OH_MONTH = month (January-December = 1-12)
! OH_SEASON = (Northern Hemispheric) season
!        1 --> winter (Dec, Jan, Feb)
!        2 --> spring (Mar, Apr, May)
!        3 --> summer (Jun, Jul, Aug)
!        4 --> autumn (Sep, Oct, Nov)
!
!***************************************************
! 1) Miscellaneous Functions.
!***************************************************
! PRESS > 100 mb:  Don't calculate OH over 100 mb.

       IF (OH_PRESS.LT.PRESSES(3)) THEN
!
! Error Check.
!       
       IF (0) THEN  
         IF (MAPL_AM_I_ROOT()) THEN
            PRINT*,'************************************************'
            PRINT*,'WARNING MESSAGE - SR SKIP!'
            PRINT*,'The pressure for this box is less than 100 mb.'
            PRINT*,'No parameterization is designed for pressures'
            PRINT*,'less than 100 mb. Box is skipped and [OH] set'
            PRINT*,'to low value. '
            PRINT*,'LON BOX=',COLON
            PRINT*,'LAT BOX=',COLAT
            PRINT*,'ALT BOX=',COALT
            PRINT*,'************************************************'
         END IF
 	END IF
!
         DOskip=1
!
         GOTO 320
!
       ENDIF
!
! NOX < 1 ppt:  Reset NOx to 1 ppt.
!
       IF (OH_NOX.LT.1.) THEN
         OH_NOX=1.1

         INDVARA(COLON,COLAT,COALT,3)=OH_NOX
         INDVARB(COLON,COLAT,COALT,3)=OH_NOX
         INDVARC(COLON,COLAT,COALT,3)=OH_NOX
         INDVARD(COLON,COLAT,COALT,3)=OH_NOX
       ENDIF
!
!***************************************************
! 2) LOW OH
!***************************************************
! SH (80-90S)
      IF (OH_LAT.LE.-80.) THEN
         COCOUNT=0
         DOskip=2
         GOTO 320
         GOTO 320
      ENDIF
! NH (80-90N)
      IF (OH_LAT.GE.80.) THEN
         COCOUNT=0
         DOskip=2
        GOTO 320
      ENDIF
! WIN NH (40-90N)
      IF (OH_LAT.GE.ALATS(6).AND.OH_SEASON.EQ.1) THEN
         COCOUNT=0
         DOskip=2
         GOTO 320
      ENDIF
! SUM & SPR SH (60-90S)
      IF (OH_LAT.LE.ALATS(2).AND.(OH_SEASON.EQ.3.OR.  &
     &   OH_SEASON.EQ.2)) THEN
         COCOUNT=0
         DOskip=2
         GOTO 320
      ENDIF
! AUT NH (60-90N)
      IF (OH_LAT.GE.ALATS(7).AND.OH_SEASON.EQ.4) THEN
         COCOUNT=0
         DOskip=2
         GOTO 320
      ENDIF
!
!***************************************************
! 3) RESET
!***************************************************
!***************************************************
! ISOPRENE
!***************************************************
! [NOx] < 50 ppt
!    [ISOP] > 600 ppt set [OH] = 1e10^5 and skip parameterization
!    [ISOP] > 3000 ppt set [OH] = 1e10^4 and skip parameterization
! Assumptions based on:
! Mickley, L. J., D. J. Jacob, and D. Rind,
! Uncertainty in preindustrial abundance of tropospheric ozone:
! Implications for radiative forcing calculations.
! In press Journal of Geophysical Research, 2000.
!
      IF ( OH_NOX.LT.50..AND.OH_PRESS.GE.PRESSES(1) ) THEN
       IF ( OH_ISOP.GE.600.) THEN
         DOskip=3
         GOTO 320
       ENDIF
       IF ( OH_ISOP.GE.3000.) THEN
         DOskip=4
         GOTO 320
       ENDIF
      ENDIF
!
! These statements ensure isoprene parameterizations not called
!    for domains with no isoprene parameterizations.
!    40-90N WIN, 60-90N SPR, 60-90N AUT, 60-90S All Seasons
!
      IF (OH_LAT.GE.ALATS(NLATS-1).AND.OH_SEASON.EQ.1) OH_ISOP=1.
      IF (OH_LAT.GE.ALATS(NLATS).AND.OH_SEASON.NE.3) OH_ISOP=1.
      IF (OH_LAT.LE.ALATS(2).AND.OH_SEASON.EQ.2) OH_ISOP=1.
      IF (OH_LAT.LE.ALATS(2).AND.OH_SEASON.EQ.3) OH_ISOP=1.
      IF (OH_LAT.LE.ALATS(1)) OH_ISOP=1.
!
!***************************************************
! These statements ensure the MODERATELY POLLUTED parameterizations
!    not called for domains with no parameterization.
!    Surface and Middle Troposphere.
!***************************************************
!
! For surface.
!
      IF ( OH_NOX.GT.CNOXS(1).AND.  &
     &  OH_NOX.LE.CNOXS(3).AND.OH_PRESS.GE.PRESSES(1)) THEN
!
         NODOMOD=0
!
       IF (OH_LAT.LE.ALATS(1)) NODOMOD=1
       IF (OH_LAT.LE.ALATS(2).AND.OH_LAT.GT.ALATS(1)) THEN
          IF (OH_SEASON.EQ.3.OR.OH_SEASON.EQ.2) NODOMOD=1
       ENDIF
       IF (OH_LAT.GT.ALATS(6).AND.OH_SEASON.EQ.1) NODOMOD=1
       IF (OH_LAT.GT.ALATS(7).AND.OH_SEASON.NE.3) NODOMOD=1
!
          IF (NODOMOD.NE.0) THEN
                OH_NOX=CNOXS(1)-10.
                INDVARA(COLON,COLAT,COALT,3)=OH_NOX
                INDVARB(COLON,COLAT,COALT,3)=OH_NOX
                INDVARC(COLON,COLAT,COALT,3)=OH_NOX
                INDVARD(COLON,COLAT,COALT,3)=OH_NOX
          ENDIF
!
      ENDIF
!
!******
!
! For Middle Troposphere.
!
      IF (OH_NOX.GT.CNOXS(1).AND.OH_PRESS.LT.PRESSES(1).AND. &
          OH_PRESS.GE.PRESSES(2)) THEN
!
         NODOMOD=0
!
       IF (OH_LAT.LE.ALATS(3)) NODOMOD=1
       IF (OH_LAT.LE.ALATS(4).AND.OH_LAT.GT.ALATS(3)) THEN
           IF (OH_SEASON.EQ.1) NODOMOD=1
           IF (OH_SEASON.EQ.2) NODOMOD=1
       ENDIF
       IF (OH_LAT.LE.ALATS(5).AND.OH_LAT.GT.ALATS(4)) THEN
           IF (OH_SEASON.EQ.3) NODOMOD=1
           IF (OH_SEASON.EQ.4) NODOMOD=1
       ENDIF
       IF (OH_LAT.GT.ALATS(6).AND.OH_SEASON.EQ.1) NODOMOD=1
       IF (OH_LAT.GT.ALATS(7)) NODOMOD=1
!
          IF (NODOMOD.NE.0) THEN
                OH_NOX=CNOXS(1)-10.
                INDVARA(COLON,COLAT,COALT,3)=OH_NOX
                INDVARB(COLON,COLAT,COALT,3)=OH_NOX
                INDVARC(COLON,COLAT,COALT,3)=OH_NOX
                INDVARD(COLON,COLAT,COALT,3)=OH_NOX
          ENDIF
!
! Reset boxes where NOx > 1000 in MT.
          IF (OH_NOX.GT.CNOXS(3)) THEN
                OH_NOX=CNOXS(3)-10.
                INDVARA(COLON,COLAT,COALT,3)=OH_NOX
                INDVARB(COLON,COLAT,COALT,3)=OH_NOX
                INDVARC(COLON,COLAT,COALT,3)=OH_NOX
                INDVARD(COLON,COLAT,COALT,3)=OH_NOX
          ENDIF
!
         GOTO 320
      ENDIF
!
!******
!
! For Upper Troposphere.
!
      IF (OH_NOX.GT.CNOXS(2).AND.OH_PRESS.LT.PRESSES(2)) THEN
          OH_NOX=CNOXS(2)-10.
          INDVARA(COLON,COLAT,COALT,3)=OH_NOX
          INDVARB(COLON,COLAT,COALT,3)=OH_NOX
          INDVARC(COLON,COLAT,COALT,3)=OH_NOX
          INDVARD(COLON,COLAT,COALT,3)=OH_NOX
         GOTO 320
      ENDIF
!
!***************************************************
! These statements ensure the E) HEAVILY POLLUTED parameterizations
!    not called for domains with no parameterization.
!    Very hi NOx associated with low OH, especially in wintertime.
!***************************************************
!
! Winter months:  high NOx=low OH.
!
      IF (OH_NOX.GT.CNOXS(3)) THEN
        IF (OH_LAT.GE.ALATS(7).AND.OH_SEASON.EQ.2) THEN
!bnd            PARAMOH=0.
            PARAMOH=1.e3
            DOskip=1
         GOTO 320
        ENDIF
        IF (OH_LAT.GE.ALATS(5).AND.OH_SEASON.EQ.1) THEN
!bnd            PARAMOH=0.
            PARAMOH=1.e3
            DOskip=1
         GOTO 320
        ENDIF
        IF (OH_LAT.LE.ALATS(3).AND.OH_SEASON.EQ.3) THEN
!bnd            PARAMOH=0.
            PARAMOH=1.e3
            DOskip=1
         GOTO 320
        ENDIF
      ENDIF
!
      IF (OH_NOX.GT.CNOXS(4)) THEN
!
! Assume hi NOx = low OH for areas without heavily polluted parameterizations.
!  
       IF (OH_LAT.LT.ALATS(3).OR.OH_LAT.GE.ALATS(7).OR.  &
     &    (OH_LAT.GE.ALATS(6).AND.  &
     &    OH_SEASON.LE.2).OR.(OH_LAT.GE.ALATS(5).AND.  &
     &    OH_SEASON.EQ.1)) THEN
         DOskip=1
         GOTO 320
!
       ELSE
!
! Reset NOx < 5 ppbv in regions with heavily polluted parameterizations.
!
          OH_NOX=CNOXS(4)-10.
          INDVARA(COLON,COLAT,COALT,3)=OH_NOX
          INDVARB(COLON,COLAT,COALT,3)=OH_NOX
          INDVARC(COLON,COLAT,COALT,3)=OH_NOX
          INDVARD(COLON,COLAT,COALT,3)=OH_NOX

!
       ENDIF
!
      ENDIF
!
! Make sure not to call nonexistent parameterization. 
!
      IF (OH_NOX.GT.CNOXS(3).AND.OH_NOX.LT.CNOXS(4)) THEN
!       
         NODOHI=0
!
       IF (OH_LAT.LE.ALATS(3))  NODOHI=1
       IF (OH_LAT.GT.ALATS(5).AND.OH_SEASON.EQ.1) NODOHI=1
       IF (OH_LAT.GT.ALATS(7))  NODOHI=1
       IF (OH_LAT.GT.ALATS(6).AND.OH_SEASON.EQ.2) NODOHI=1

       IF (NODOHI.NE.0) THEN
          OH_NOX=CNOXS(3)-10.
          INDVARA(COLON,COLAT,COALT,3)=OH_NOX
          INDVARB(COLON,COLAT,COALT,3)=OH_NOX
          INDVARC(COLON,COLAT,COALT,3)=OH_NOX
          INDVARD(COLON,COLAT,COALT,3)=OH_NOX
          GOTO 320
       ENDIF
!
      ENDIF
!
!***************************************************
!
 320  RETURN
      END SUBROUTINE SKIP

   END MODULE  getoh_mod
