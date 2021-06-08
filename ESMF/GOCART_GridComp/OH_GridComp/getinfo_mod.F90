#ifdef GEOS5
#include "MAPL_Generic.h"
#endif
!-------------------------------------------------------------------------
! NASA GSFC - SIVO Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: module ohParameterizationMethod_mod
!
! !INTERFACE:
!
      module getinfo_mod
!
! !USES:
!
#ifdef GEOS5
      USE ESMF
      USE MAPL
#endif

      USE cblock_size_mod         ! "CMN_SIZE"       ! Size parameters
      USE cblock_CO_mod           ! "CMN_CO"         ! CO arrays
      USE cblock_CO_budget_mod    ! "CMN_CO_BUDGET"  ! FMOL_CO
      USE cblock_OH_mod           ! "CMN_OH"

      implicit none

! ****************************************************************************
      CONTAINS
! ****************************************************************************
      SUBROUTINE GETINFO(OH_MONTH2,ALBD, &
                         INDVARA, INDVARB, INDVARC, INDVARD, &
                         Wavg, Pavg, Tavg, Ravga, Ravgb,     &
                         STT, BBIJ, o3up, &
                         METHVAR, RLON, RLAT, JDAY, &
                         MLONBX, MLATBX, MVRTBX)
! ****************************************************************************
! THIS SUBROUTINE IS TO BE MODIFIED BY THE USER!
! The user can interface it with a model by passing necessary
! information through the subroutine's argument list or
! through the user's own common blocks.
! ****************************************************************************
! ****************************************************************************
! In this subroutine, the user supplies the parameterization code
! with information it needs to calculate the 24-hour averaged [OH].
! ****************************************************************************
! Created by Bryan Duncan.
! ****************************************************************************
! The user is to put whatever common blocks and declaration 
! statements here which are necessary to supply this subroutine
! with the model input it requires.
!

      IMPLICIT NONE

      INTEGER, intent(in) :: MLONBX, MLATBX, MVRTBX
      INTEGER, intent(in) :: JDAY
      INTEGER, intent(in) :: OH_MONTH2
      REAL*8,  intent(in) :: Wavg(MLONBX, MLATBX, MVRTBX)  ! 24h avg spc. humidity
      REAL*8,  intent(in) :: Tavg(MLONBX, MLATBX, MVRTBX)  ! 24h avg temperature
      REAL*8,  intent(in) :: o3up       (MLONBX, MLATBX, MVRTBX)
      REAL*8,  intent(in) :: STT        (MLONBX, MLATBX, MVRTBX,NNPAR)
      REAL*8,  intent(in) :: RLON       (MLONBX, MLATBX)   ! Array of longitudes
      REAL*8,  intent(in) :: RLAT       (MLONBX, MLATBX)   ! Array of latitudes
      REAL*8,  intent(in) :: METHVAR    (MLONBX, MLATBX, MVRTBX)
      ! BBIJ stores species concentration for:
      ! 1. NOt             2. ALK4          3. Isoprene
      ! 4. Acetone         5. Propene       6. Propane
      ! 7. Ethane          8. O3            9. CO
      REAL*8,  intent(in) :: BBIJ       (MLONBX, MLATBX, MVRTBX,NFIELDS2)

      REAL*8,  intent(inOut) :: Ravga(MLONBX, MLATBX, MVRTBX) ! 24h avg cloud flux above box
      REAL*8,  intent(inOut) :: Ravgb(MLONBX, MLATBX, MVRTBX) ! 24h avg cloud flux below box
      REAL*8,  intent(inOut) :: Pavg(MLONBX, MLATBX, MVRTBX)  ! 24h avg pressure
      REAL*8,  intent(inOut) :: INDVARA(MLONBX, MLATBX, MVRTBX, MXVAR)
      REAL*8,  intent(inOut) :: INDVARB(MLONBX, MLATBX, MVRTBX, MXVAR)
      REAL*8,  intent(inOut) :: INDVARC(MLONBX, MLATBX, MVRTBX, MXVAR)
      REAL*8,  intent(inOut) :: INDVARD(MLONBX, MLATBX, MVRTBX, MXVAR)

      REAL*8   ALBD(MLONBX, MLATBX)

! !LOCAL VARIABLES:
      REAL*8   PI,DEC,A0,A1,A2,A3,B1,B2,B3,R,BLAT
      PARAMETER ( PI = 3.141592 )
      INTEGER  KK,NDOCO,NCHECKIT,NDOCH4
      REAL*8   AISOP,ddd,sumo3up
      REAL*8 dummyo3(MLONBX, MLATBX)

      REAL*8 :: STTTOPPB   (MLONBX, MLATBX, MVRTBX)

      ! Weight of air (taken from "comode.h") 
      REAL*8, PARAMETER :: WTAIR = 28.966d0
!
! ****************************************************************************
! This common block and these declaration statements are not to be
! modified by the user.
!
      INTEGER I,J,K,L,NSTOP
!
! ****************************************************************************
! The polynomials that are used to calculate the parameterized
! OH concentrations are a function of the following independent
! variables:
!
!   Independent Variables                   Units
!   -----------------------------------------------------------
!   Solar Declination Angle                 degrees
!   Total Ozone Column Above                DU
!   Total Nitrogen Oxides                   pptv
!   Ozone                                   ppbv
!   Carbon Monoxide                         ppbv
!   Methane                                 ppbv
!   Water Vapor                             ppmv
!   Acetone                                 pptv
!   Pressure                                mb
!   Propene                                 pptv
!   Ethane                                  pptv
!   Propane                                 pptv
!   ALK4                                    pptv
!   Surface Albedo                          unitless
!   Latitude                                degrees
!   Temperature                             K
!   Cloude Albedo - Above                   unitless
!   Cloude Albedo - Below                   unitless
!   Isoprene (daylight)                     pptv
!   -----------------------------------------------------------
!
! All these independent variables need to be 24-hour averages!!
!
! All the independent variables are REAL*8!
!
! ****************************************************************************
!
! The independent variables are stored in the arrays INDVARA-D:
!
! INDVARA is used in surface layer parameterizations without isoprene.
! INDVARB is used in middle and upper tropospheric layer parameterizations
!         without isoprene.
! INDVARC is used in surface layer parameterizations with isoprene.
! INDVARD is used in middle and upper tropospheric layer parameterizations
!         with isoprene.
!
! --------------------------------------------------------------------------
!   MXVAR    INDVARA        INDVARB        INDVARC         INDVARD
! --------------------------------------------------------------------------
!     1       DEC            DEC            DEC             DEC
!     2       O3COL          O3COL          O3COL           O3COL
!     3       NOt            NOt            NOt             NOt
!     4       O3             O3             O3              O3
!     5       CO             CO             CO              CO
!     6       CH4            CH4            CH4             CH4
!     7       H2O            H2O            H2O             H2O
!     8       ACET           ACET           ACET            ACET
!     9       PRESS          PRESS          PRESS           PRESS
!    10       PRPE           PRPE           PRPE            PRPE
!    11       ETHA           ETHA           ETHA            ETHA
!    12       PROP           PROP           PROP            PROP
!    13       ALK4           ALK4           ALK4            ALK4
!    14       ALB            ALB            ALB             ALB
!    15       LAT            LAT            LAT             LAT
!    16       TEMP           TEMP           TEMP            TEMP
!    17       CLOUDA         CLOUDA         CLOUDA          CLOUDA
!    18         *            CLOUDB         ISOP            CLOUDB
!    19         *              *              *             ISOP
! --------------------------------------------------------------------------
!   *No independent variable value placed in this position.
!
! ****************************************************************************
! Error Check.  If the user does not totally fill in the INDVARA-D arrays,
! then the code will stop to tell the user.
!
      IF(MAPL_AM_I_ROOT()) PRINT*,'ENTERING SR GETINFO'
     
      INDVARA(:,:,:,:)=-1000.D0
      INDVARB(:,:,:,:)=-1000.D0
      INDVARC(:,:,:,:)=-1000.D0
      INDVARD(:,:,:,:)=-1000.D0
!
      INDVARA(:,:,:,18:19)=0.D0
      INDVARB(:,:,:,19)=0.D0
      INDVARC(:,:,:,19)=0.D0
!
! ****************************************************************************

!
! (1) Solar declination angle (low precision formula)
!     JDAY = julian day
!
      A0 = 0.006918
      A1 = 0.399912
      A2 = 0.006758
      A3 = 0.002697
      B1 = 0.070257
      B2 = 0.000907
      B3 = 0.000148
 
      R  = 2.* PI * float(JDAY-1) / 365.
!
      DEC = A0 - A1*COS(  R) + B1*SIN(  R) &
     &         - A2*COS(2.*R) + B2*SIN(2.*R) &
     &         - A3*COS(3.*R) + B3*SIN(3.*R)
!
      DEC = DEC*180./PI
!      TESTDEC(:,:) = DEC
!
      IF (DEC.GT.24..OR.DEC.LT.-24.) THEN
         PRINT*,'Stopped in SR GETINFO!'
         PRINT*,'  Declination angle is too high!'
         PRINT*,'DEC = ',DEC
              CALL FLUSH( 6 )
         STOP
      ENDIF
!
      INDVARA(:,:,:,1)=DEC
!
! (2) O3 Column Above (DU)


!KNJR      INDVARA(:,:,:,2)=o3up(:,:,:,OH_MONTH2)
      INDVARA(:,:,:,2)=o3up(:,:,:)

! (3) NOt (pptv)
! Note:  if you change the units of NOx here, you need to update
!        the NOx in the isoprene yield calculation!
!G
!        NOt = NO + NO2 + 2N2O5 + NO3 + HNO2 + HNO4
!
      INDVARA(:,:,:,3)=BBIJ(:,:,:,1)

!
      DO I=1,MLONBX
         DO J=1,MLATBX
            DO L=1,5  
               IF (INDVARA(I,J,L,3).LT.1.D-5) THEN
                  PRINT*,'Stopped in SR getinfo.f!'
                  PRINT*,'NOx is mixing ratio for box:'
                  PRINT*,I,J,L,INDVARA(I,J,L,3)
                  CALL FLUSH( 6 )
!bnd               STOP
               ENDIF 
            ENDDO
         ENDDO
      ENDDO
!       
! (4) O3 (ppbv)
!

      INDVARA(:,:,:,4)=BBIJ(:,:,:,8)

      DO I=1,MLONBX
         DO J=1,MLATBX
            DO L=1,9           
               IF (INDVARA(I,J,L,4).LT.1.D-5) THEN
                  PRINT*,'Stopped in SR getinfo.f!'
                  PRINT*,'O3 is mixing ratio for box:'
                  PRINT*,I,J,L,INDVARA(I,J,L,4)
                  CALL FLUSH( 6 )
                  STOP
               ENDIF
            ENDDO
         ENDDO
      ENDDO

!
! (5) CO (ppbv)
!
      NDOCO = 0
      IF (NDOCO.EQ.0) THEN
!
! Variable CO.

         INDVARA(:,:,:,5)=STT(:,:,:,1)
!
      ENDIF

      DO I=1,MLONBX
         DO J=1,MLATBX
            DO L=1,9    
	         
               IF (INDVARA(I,J,L,5).LT.1.D-5) THEN
                  
		  IF(MAPL_AM_I_ROOT()) PRINT*,'Stopped in SR getinfo_mod.F90!'
                  IF(MAPL_AM_I_ROOT()) PRINT*,'CO is lt 1e-5 for box:'
                  IF(MAPL_AM_I_ROOT()) PRINT*,I,J,L,INDVARA(I,J,L,5)

                  CALL FLUSH( 6 )
                  !STOP
               ENDIF
            ENDDO
         ENDDO
      ENDDO

!
! (6) CH4 (ppbv)
! See SR CO_fromHCs to see where I got the [CH4] from.
!

!************
! Daily mean CH4 (ppbv) ppbv
      INDVARA(:,:,:,6) = METHVAR(:,:,:) 

! (7) H2O (ppmv)
!
      INDVARA(:,:,:,7)=Wavg(:,:,:)
      
!
! (8) Acetone (pptv) 
!
      
      INDVARA(:,:,:,8)=BBIJ(:,:,:,4)
      DO I=1,MLONBX
         DO J=1,MLATBX
            DO L=1,9          
               IF(INDVARA(I,J,L,8).LT.1.D-5) THEN
                  PRINT*,'Stopped in SR getinfo.f!'
                  PRINT*,'ACET is mixing ratio for box:'
                  PRINT*,I,J,L,INDVARA(I,J,L,8)
                  CALL FLUSH( 6 )
                  STOP
               ENDIF
            ENDDO
         ENDDO
      ENDDO

!
! (9) Pressure (mb) - 24 hr avg
!
      WHERE(Pavg(:,:,:).GT.1020.D0) Pavg(:,:,:)=1019.9D0

      DO L=1,MVRTBX
         INDVARA(:,:,L,9)=Pavg(:,:,L)  
      ENDDO
!
! (10) Propene (pptv)
!
      INDVARA(:,:,:,10)=BBIJ(:,:,:,5)
!
! (11) Ethane (pptv)
!
      INDVARA(:,:,:,11)=BBIJ(:,:,:,7)

      DO I=1,MLONBX
         DO J=1,MLATBX
            DO L=1,9          
               IF(INDVARA(I,J,L,11).LT.1.D-5) THEN
                  PRINT*,'Stopped in SR getinfo.f!'
                  PRINT*,'ETHA is mixing ratio for box:'
                  PRINT*,I,J,L,INDVARA(I,J,L,11)
                  CALL FLUSH( 6 )
                  STOP
               ENDIF
            ENDDO
         ENDDO
      ENDDO

!
! (12) Propane (pptv)
!
      INDVARA(:,:,:,12)=BBIJ(:,:,:,6)

      DO I=1,MLONBX
         DO J=1,MLATBX
            DO L=1,9         
               IF(INDVARA(I,J,L,12).LT.1.D-5) THEN
                  PRINT*,'Stopped in SR getinfo.f!'
                  PRINT*,'PROP is mixing ratio for box:'
                  PRINT*,I,J,L,INDVARA(I,J,L,12)
                  CALL FLUSH( 6 )
                  STOP
               ENDIF
            ENDDO
         ENDDO
      ENDDO

!
! (13) ALK4 (pptv)
!     
      INDVARA(:,:,:,13)=BBIJ(:,:,:,2)

      DO I=1,MLONBX
         DO J=1,MLATBX
            DO L=1,9
               IF(INDVARA(I,J,L,13).LT.1.D-5) THEN
                  PRINT*,'Stopped in SR getinfo.f!'
                  PRINT*,'ALK4 is mixing ratio for box:'
                  PRINT*,I,J,L,INDVARA(I,J,L,13)
                  CALL FLUSH( 6 )
!bnd              STOP
               ENDIF
            ENDDO
         ENDDO
      ENDDO
!
! (14) Albedo (fraction of 1)
! While visible albedo rarely goes below 0.06 (6%),
! the uv albedo does go below 0.06 over land occasionally.
! In general though, most of the earth's surface uv
! albedo is above this 0.06 minimum set for the params
! initially designed for use with visible albedo.
!
!vis       WHERE(AVGF(:,:).LT.0.06D0) AVGF(:,:)=0.061D0
      WHERE(ALBD(:,:).LT.0.06D0) ALBD(:,:)=0.061D0
      
!vis       WHERE(AVGF(:,:).GT.0.80D0) AVGF(:,:)=0.79D0
      WHERE(ALBD(:,:).GT.0.80D0) ALBD(:,:)=0.79D0

      DO L=1,MVRTBX
!vis       INDVARA(:,:,L,14)=AVGF(:,:)
         INDVARA(:,:,L,14)=ALBD(:,:)
      ENDDO
!
! (15) Latitude (deg = -90 to 90)
!
! EY DEBUG - INSERT LATITUDE ARRAY HERE  ; Originally - 2 degrees at the poles ; 4 degrees elsewehere (4x4.5 grid)
      !BLAT=35.D0  ! BLAT SET HERE
      DO J=1,MLATBX
       DO I=1,MLONBX
        !IF(J.NE.MLATBX.AND.J.NE.1) THEN
        !   BLAT=BLAT+4.D0
        !ELSE
        !   BLAT=BLAT+2.D0
        !ENDIF
        !IF(J.EQ.2) BLAT=BLAT-1.D0
        !IF(J.EQ.MLATBX) BLAT=BLAT+1.D0
         INDVARA(I,J,:,15)=RLAT(I,J) ! BLAT
       ENDDO
      ENDDO
      CONTINUE
      
!
! (16) Temperature (K)
!
      DO L=1,MVRTBX
         INDVARA(:,:,L,16)=Tavg(:,:,L)  
      ENDDO
!
! (17) Cloud flux fraction above box
!
      WHERE(Ravga(:,:,:).LT.0.D0) Ravga(:,:,:)=0.D0

      WHERE(Ravga(:,:,:).GT.0.6D0) Ravga(:,:,:)=0.59999D0

      INDVARA(:,:,:,17)=Ravga(:,:,:)
      
!*****************************************************************************
!
! The first 17 independent variables in INDVARA are the
! same for INDVARB, INDVARC and INDVARD.  However, they
! are different for the 18-19 independent variables.
!
      INDVARB(:,:,:,1:17)=INDVARA(:,:,:,1:17)
      INDVARC(:,:,:,1:17)=INDVARA(:,:,:,1:17)
      INDVARD(:,:,:,1:17)=INDVARA(:,:,:,1:17)
!
!*****************************************************************************
!
! (18) Cloud flux fraction below box
!
      WHERE(Ravgb(:,:,:).LT.0.D0) Ravgb(:,:,:)=0.D0

      WHERE(Ravgb(:,:,:).GT.0.6D0) Ravgb(:,:,:)=0.59999D0

      INDVARB(:,:,:,18)=Ravgb(:,:,:)
      INDVARD(:,:,:,18)=Ravgb(:,:,:)

!*****************************************************************************
! Resetting cloud flux fractions to 0 for test!
!      INDVARA(:,:,:,17)=0.0001
!      INDVARB(:,:,:,17)=0.0001
!      INDVARC(:,:,:,17)=0.0001
!      INDVARD(:,:,:,17)=0.0001
!      INDVARB(:,:,:,18)=0.0001
!      INDVARD(:,:,:,18)=0.0001
!*****************************************************************************
   

! (19) Isoprene (pptv) - daylight average?
!
      DO L=1,MVRTBX
         DO J=1,MLATBX
            DO I=1,MLONBX
!

               AISOP= BBIJ(I,J,L,3) ! sas changed back 2.D0 ! EY DEBUG BBIJ(I,J,L,3)
!
! are you reading in 24-Hr averages or daylight averages?
! see CO_fillfields.f

               IF(AISOP.LT.2.D0) AISOP= 2.D0
               INDVARC(I,J,L,18)=AISOP
               INDVARD(I,J,L,19)=AISOP
            ENDDO
         ENDDO
      ENDDO
!
! ****************************************************************************
! Error Check. Print out input variables.
!
       NCHECKIT=0
!
      IF (NCHECKIT.EQ.1) THEN       
!
         PRINT*,'DOING ERROR CHECK BY PRINTING fort.819!'
!
         open(65,file='fort.819',status='unknown')

         DO I=1,MLONBX
            DO J=1,MLATBX
               DO L=1,MVRTBX
                  DO KK=1,19
!c       if(INDVARD(I,J,L,3).le.350..and.
!c     *    INDVARD(I,J,L,19).gt.11.) then     
                     WRITE(65,*) INDVARD(I,J,L,KK)
!c       else
!c         ddd=0. 
!c        write(819,*) ddd
!c       endif
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         close(65)
         CALL FLUSH( 6 )
         STOP 
!
      ENDIF
!
! ***************************************************************************
! ***************************************************************************
! Other parameters needed.
! ***************************************************************************
! ***************************************************************************
!
! A) Longitude (deg = -180 - 180)
!   Longitude is not an independent variable used in the polynomials
!   to calculate the parameterized OH, but is used to locate 
!   regions of biomass burning.  There are separate parameterizations
!   for regions affected by biomass burning.  For convenience, however,
!   longitude is stored in the INDVAR* arrays.
!
! EY DEBUG SET LON HERE
      !BLON=-90.D0
      DO J=1,MLATBX
       DO I=1,MLONBX
         !BLON=BLON+5.D0

         INDVARA(I,J,:,20)=RLON(I,J)
         INDVARB(I,J,:,20)=RLON(I,J)
         INDVARC(I,J,:,20)=RLON(I,J)
         INDVARD(I,J,:,20)=RLON(I,J)

       ENDDO
      ENDDO
!
! B) Season
!     OH_MONTH2 = month (January-December = 1-12)
!     OH_SEASON = Northern Hemispheric season
!        1 --> winter (Dec, Jan, Feb)
!        2 --> spring (Mar, Apr, May)
!        3 --> summer (Jun, Jul, Aug)
!        4 --> autumn (Sep, Oct, Nov)
!
! From month, determine season.  The parameterizations are
! functions of season and not month.
!
      IF(OH_MONTH2.LE.2)  OH_SEASON=1
      IF(OH_MONTH2.EQ.12) OH_SEASON=1
      IF(OH_MONTH2.LE.5.AND.OH_MONTH2.GE.3)  OH_SEASON=2
      IF(OH_MONTH2.LE.8.AND.OH_MONTH2.GE.6)  OH_SEASON=3
      IF(OH_MONTH2.LE.11.AND.OH_MONTH2.GE.9) OH_SEASON=4
!     
!*****************************************************************************
! Error Check.
!     
      NSTOP=0
!
      DO I = 1,MXVAR
         DO J = 1,MVRTBX
            DO K = 1,MLATBX
               DO L = 1,MLONBX

!bnd         IF ( IEEE_IS_NAN( INDVARA(L,K,J,I) ) ) THEN
!bnd           PRINT*,'Stopped in SR getinfo.'
!bnd           PRINT*, 'INDVARA(:,:,:,I) is NaN!'
!bnd           PRINT*,I,L,K,J,INDVARA(L,K,J,I)
!bnd              CALL FLUSH( 6 )
!bnd           STOP
!bnd         ENDIF
!
                  IF(INDVARA(L,K,J,I).LT.-200.) THEN
                     PRINT*,'******************************************'
                     PRINT*,'ERROR IN SR GETINFO!'
                     PRINT*,'The array INDVARA not completely filled.'
                     PRINT*,'Check independent variable number = ',I
                     PRINT*,'vert/lat/lon = ',J,K,L
                     PRINT*,'INDVARA(L,K,J,I)=',INDVARA(L,K,J,I)
                     NSTOP=1
                  ENDIF
!
                  IF(INDVARB(L,K,J,I).LT.-200.) THEN
                     PRINT*,'******************************************'
                     PRINT*,'ERROR IN SR GETINFO!'
                     PRINT*,'The array INDVARB not completely filled.'
                     PRINT*,'Check independent variable number = ',I
                     PRINT*,'vert/lat/lon = ',J,K,L
                     PRINT*,'INDVARB(L,K,J,I)=',INDVARB(L,K,J,I)
                     NSTOP=1
                  ENDIF
!
                  IF(INDVARC(L,K,J,I).LT.-200.) THEN
                     PRINT*,'******************************************'
                     PRINT*,'ERROR IN SR GETINFO!'
                     PRINT*,'The array INDVARC not completely filled.'
                     PRINT*,'Check independent variable number = ',I
                     PRINT*,'vert/lat/lon = ',J,K,L
                     PRINT*,'INDVARC(L,K,J,I)=',INDVARC(L,K,J,I)
                     NSTOP=1
                  ENDIF
!
                  IF(INDVARD(L,K,J,I).LT.-200.) THEN
                     PRINT*,'******************************************'
                     PRINT*,'ERROR IN SR GETINFO!'
                     PRINT*,'The array INDVARD not completely filled.'
                     PRINT*,'Check independent variable number = ',I
                     PRINT*,'vert/lat/lon = ',J,K,L
                     PRINT*,'INDVARD(L,K,J,I)=',INDVARD(L,K,J,I)
                     NSTOP=1
                  ENDIF
!
                  IF(NSTOP.EQ.1) THEN
                     PRINT*,'L,K,J=',L,K,J
                     PRINT*,'******************************************'
                       CALL FLUSH( 6 )
                     STOP
                  ENDIF
!
               ENDDO
            ENDDO
         ENDDO
      ENDDO
!
      IF(OH_SEASON.LT.1.OR.OH_SEASON.GT.4) THEN
         PRINT*,'******************************************'
         PRINT*,'ERROR IN SR GETINFO!'
         PRINT*,'OH_SEASON does not have a value between 1 & 4.'
         PRINT*,'OH_SEASON = ',OH_SEASON
         PRINT*,'******************************************'
              CALL FLUSH( 6 )
         STOP
      ENDIF
!
!*****************************************************************************
!
      IF(MAPL_AM_I_ROOT()) PRINT*,'LEAVING SR GETINFO'

      RETURN
      END SUBROUTINE GETINFO
!
!*****************************************************************************
!
      end module getinfo_mod
