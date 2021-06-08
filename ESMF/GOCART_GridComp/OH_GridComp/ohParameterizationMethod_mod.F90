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
module ohParameterizationMethod_mod
!
! !USES:
!
#ifdef GEOS5
      USE ESMF
      USE MAPL
#endif
      USE m_inpak90 
      USE getoh_mod
      USE getinfo_mod
      USE inputFileReading_mod
      USE cblock_OH_mod    ! "CMN_OH"
!
      implicit none
!
!
! !PUBLIC MEMBER FUNCTIONS:
!
      private
      public  :: initializeOHparam
      public  :: runOHparam
      public  :: finalizeOHparam
!
! !PUBLIC DATA MEMBERS:
!
      public  :: t_OHparam

      type t_OHparam
         integer :: MVRTBX  ! # of boxes in vertical
         integer :: MLATBX  ! # of boxes following one longitude band
         integer :: MLONBX  ! # of boxes following one latitude band
              ! the vertical level of the tropopause.  Above this level,
              ! no [OH] is calculated.  The user can feed this SR
              ! a high value for LTPAUSE which effectively turns
              ! this option off (i.e., LTPAUSE > MVRTBX).  Note that
              ! this parameterization is only designed for pressures 
              ! greater than 100 mb, however.  If the
              ! [OH] = -9999 then the [OH] was not calculated.
              ! LTPAUSE assumes a dimension of longitude and latitude.
              ! The user may change this dimension to one- or zero-,
              ! however, the if statement surrounding the call to
              ! SR GETOH needs to be modified accordingly.
         integer, pointer :: LTPAUSE(:,:) => null()
              ! array holding tropospheric parameterized OH.
         real*8,  pointer :: BOH    (:,:,:)    => null()
         real*8,  pointer :: INDVARA(:,:,:,:)  => null()
         real*8,  pointer :: INDVARB(:,:,:,:)  => null()
         real*8,  pointer :: INDVARC(:,:,:,:)  => null()
         real*8,  pointer :: INDVARD(:,:,:,:)  => null()
         real*8,  pointer :: dustrat(:,:,:)  => null()
         character (len=255) :: srefl_inFileName
         character (len=255) :: avgOH_inFileName
         character (len=255) :: ohStart_inFileName
         character (len=255) :: correctRatio_inFileName
         character (len=255) :: allAerosolRatio_inFileName
         INTEGER :: FIRSTDT
      end type t_OHparam
! !DESCRIPTION:
! We present a parameterization for the tropospheric 
! concentration of the hydroxyl radical (OH) which can be used
! to overcome the costs of solving kinetic equations in chemical 
! tracer models.  This parameterization accurately represents OH
! predicted by a full chemical mechanism.  The 24-hour average
! concentration of OH is represented as a set of high-order 
! polynomials in variables such as temperature, latitude,
! declination and the concentrations of ozone, water vapor, 
! carbon monoxide, nitrogen oxides (as a family), and
! hydrocarbons.  Results include computer-written FORTRAN 
! functions for an efficient computation of the polynomials.
!
! !AUTHOR:
! Bryan Duncan
! David Portman, David_Portman@aer.com
! Clarissa Spivakovsky,  cms@io.harvard.edu
!
! !REVISION HISTORY:
!
!EOP
!-------------------------------------------------------------------------
  CONTAINS
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: initializeOHparam
!
! !INTERFACE:
!
      SUBROUTINE initializeOHparam(self, IM, JM, LM, rcfilen)
!
! !INPUT PARAMETERS:
      INTEGER, intent(in) :: IM, JM, LM
      character (len=*)   :: rcfilen
!
! !INPUT/OUTPUT PARAMETERS:
      type (t_OHparam),   intent(inOut) :: self
!
! !DESCRIPTION:
!
! !LOCAL VARIABLES:
      INTEGER :: unitNum
      INTEGER :: ier(50)
!EOP
!-----------------------------------------------------------------------------
!BOC
      IF(MAPL_AM_I_ROOT()) PRINT *,"initializeOHparam: Begin..."

      ! Set the dimensions
      self%MLONBX = IM
      self%MLATBX = JM
      self%MVRTBX = LM

      ! Allocate variables
      !allocate(self%LTPAUSE(self%MLONBX,self%MLATBX))
      allocate(self%BOH    (self%MLONBX,self%MLATBX,self%MVRTBX))
      allocate(self%dustrat(self%MLONBX,self%MLATBX,self%MVRTBX))
      allocate(self%INDVARA(self%MLONBX,self%MLATBX,self%MVRTBX,MXVAR))
      allocate(self%INDVARB(self%MLONBX,self%MLATBX,self%MVRTBX,MXVAR))
      allocate(self%INDVARC(self%MLONBX,self%MLATBX,self%MVRTBX,MXVAR))
      allocate(self%INDVARD(self%MLONBX,self%MLATBX,self%MVRTBX,MXVAR))

      !=====================
      ! Popopulate Variables
      !=====================

      self%BOH = 0.0d0

      !=======================
      ! Load the resource file
      !=======================

      ier(:) = 0

      CALL I90_loadf ( TRIM(rcfilen), ier(1) )
      IF ( ier(1) .NE. 0 ) THEN
         PRINT*,"Can not open the file: ", rcfilen
         CALL I90_release()
         RETURN
      END IF

      !========================
      ! Parse the resource file
      !========================

      CALL I90_label ( 'ohStart_inFileName:',           ier(2) )
      CALL I90_Gtoken( self%ohStart_inFileName,         ier(3) )   ! 'start_information'

      CALL I90_label ( 'avgOH_inFileName:',             ier(4) )
      CALL I90_Gtoken( self%avgOH_inFileName,           ier(5) )   ! 'avgOH'

      CALL I90_label ( 'correctRatio_inFileName:',      ier(6) )
      CALL I90_Gtoken( self%correctRatio_inFileName,    ier(7) )   ! 'CorrectRatio'

      CALL I90_label ( 'allAerosolRatio_inFileName:',   ier(8) )
      CALL I90_Gtoken( self%allAerosolRatio_inFileName, ier(9) )   ! 'oh.allaer.ratio'

      CALL I90_label ( 'srefl_inFileName:',             ier(10))
      CALL I90_Gtoken( self%srefl_inFileName,           ier(11))   ! 'Srefl.table'

      IF ( ANY( ier(:) /= 0 ) ) THEN
         PRINT*,"Can not locate OH Parameterization input files"
         CALL I90_release()
         RETURN
      END IF

      ! MANYIN added
      CALL I90_release()

      unitNum = 99

      !==========================================
      ! Read in reflectivity table (Srefl.table).
      !==========================================

      CALL read_srefl(unitNum, self%srefl_inFileName) 

      !=====================================================
      ! Read in start up information for OH parameterization.
      !=====================================================

      CALL READ_COEFF(unitNum, self%ohStart_inFileName) 

      !===============================
      ! Read in climatological OH data.
      !===============================

      CALL readavgOH (unitNum, self%avgOH_inFileName)

      !========================================================
      ! Read in ratio to correct aerosol optical depth problem.
      !========================================================

      CALL read_correction (unitNum, self%correctRatio_inFileName)

      self%FIRSTDT = 1

      RETURN

      END SUBROUTINE initializeOHparam
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: runOHparam
!
! !INTERFACE:
!
      SUBROUTINE runOHparam(self, curMonth, ALBD, LTPAUSE,      &
                         Wavg, Pavg, Tavg, Ravga, Ravgb,        &
                         STT, BBIJ, o3up, &
                         METHVAR, JDAY, RLON, RLAT, OHparamNum)
!
! !INPUT PARAMETERS:
      INTEGER, intent(in) :: curMonth
      INTEGER, intent(in) :: JDAY
      INTEGER, intent(in) :: LTPAUSE(:,:)
      REAL*8,  intent(in) :: Wavg       (:,:,:)  ! 24h avg spc. humidity
      REAL*8,  intent(in) :: Tavg       (:,:,:)  ! 24h avg temperature
      REAL*8,  intent(in) :: METHVAR    (:,:,:)  ! 24h avg methane concentration
      REAL*8,  intent(in) :: o3up       (:,:,:)
      REAL*8,  intent(in) :: STT        (:,:,:,:)
      REAL*8,  intent(in) :: RLON       (:,:)    ! Array of longitudes
      REAL*8,  intent(in) :: RLAT       (:,:)    ! Array of latitudes

      ! BBIJ stores species concentration for:
      ! 1. NOt             2. ALK4          3. Isoprene
      ! 4. Acetone         5. Propene       6. Propane
      ! 7. Ethane          8. O3            9. CO
      REAL*8,  intent(in) :: BBIJ       (:, :, :,:)

      REAL*8   ALBD(:, :)
!
! !INPUT/OUTPUT PARAMETERS:
      type (t_OHparam),   intent(inOut) :: self
      REAL*8,  intent(inOut) :: Pavg       (:,:,:)  ! 24h avg pressure
      REAL*8,  intent(inOut) :: Ravga      (:,:,:)  ! 24h avg reflectivities  above box
      REAL*8,  intent(inOut) :: Ravgb      (:,:,:)  ! 24h avg reflectivities  below box
      REAL*8,  intent(inOut) :: OHparamNum      (:,:,:) ! number of    parametrization used


!
! !DESCRIPTION:
! Driver of the code that calculates parameterized OH:
! \begin{enumerate}
! \item SR GETOH is called for each model grid box and the value
!    of the parameterized [OH] is placed in the array BOH.
!
! \item Parameterized values of OH that are less than zero are
!    set to a small number.
! \end{enumerate}
! The code is designed for use in a 3-D tropospheric model of transport 
! and chemistry. However, if the user wishes to use this code to
! calculate only one point, for instance, then specify MVRTBX,
! MLATBX & MLONBX (i.e., the dimensions specified by
! the user in CMN_OH) equal to 1, 1 & 1.  Then call SR OHPARAM
! for each point.
!
! !LOCAL VARIABLES:
      INTEGER I,J,L,NDODUST,DOskip
      REAL*8 CF, OHParamNumOut
      REAL*8 :: Rout(self%MLONBX, self%MLATBX, self%MVRTBX)
!
! !AUTHOR:
! Bryan Duncan, bryan.n.duncan@nasa.gov
! Jules Kouatchou, Jules.Kouatchou@nasa.gov
!
!EOP
!------------------------------------------------------------------------------
!BOC

      IF (MAPL_AM_I_ROOT()) THEN
         print*,'******************************'
         print*,'Beginning OH param calculation'
         print*,'******************************'
      END IF

      CALL GETINFO(curMonth, ALBD, self%INDVARA, self%INDVARB,   &
                   self%INDVARC, self%INDVARD, Wavg, Pavg, Tavg, &
                   Ravga, Ravgb, STT, BBIJ, o3up,                &
                   METHVAR, RLON, RLAT, JDAY,       &
                   self%MLONBX, self%MLATBX, self%MVRTBX)

      ! Set BOH to large number for error check.

      self%BOH(:,:,:)= 1.E10
      ! Loop over dimensions of troposphere one box at a time.
      OHparamNum(:,:,:) = -999
      DO L=1,self%MVRTBX
!       DO L=self%MVRTBX,1, -1
         DO J=1,self%MLATBX
            DO I=1,self%MLONBX
               IF ( L .LT. LTPAUSE(I,J) ) THEN 
                  DOskip=0
		  OHparamNumOut = -999
                  CALL GETOH(I, J, L, self%FIRSTDT, self%BOH, &
                             curMonth, DOskip, self%INDVARA,  &
                             self%INDVARB, self%INDVARC,      &
                             self%INDVARD, self%MLONBX,       &
                             self%MLATBX, self%MVRTBX, OHparamNumOut)


!		  IF (DOskip.GE.1) THEN 
!		     OHparamNum(I,J,L) = 0     
!		  ELSE 
		     OHparamNum(I,J,L) =  OHparamNumOut
!		  ENDIF
!bnd
                  IF(DOskip.GE.1) GOTO 300
!bnd
!
! SR AerosolOH: Adjust OH due to presence of absorbing/reflecting
!            & heterogenous chemistry on aerosol.
!***********************************************************
                  ! NDODUST=0 Don't reduce OH.
                  ! NDODUST=1 Do reduce OH.
!
                  NDODUST=1 !sas changed then unchanged - default is DODUST=1
                  IF (NDODUST.EQ.1) THEN
		     IF (L.LE.18) self%BOH(I,J,L) = self%BOH(I,J,L) * self%dustrat(I,J,L)
                  END IF

                  ! Error Check.

                  IF (self%BOH(I,J,L).GE.1.E9) THEN
                     IF (MAPL_AM_I_ROOT()) THEN
                        PRINT*,''
                        PRINT*,'OH not calculated for this box!'
                        PRINT*,'  BOH(',I,',',J,',',L,') = ',self%BOH(I,J,L)
                     END IF
                     IF (L.GE.16) THEN
                       self%BOH(I,J,L)=1000.
!		       PRINT *, 'runOH: OH > 1E9 and L > 16 reset OH to 1000'
                     ELSE
                        PRINT*,'Stopped in SR OHPARAM'
                        PRINT*,''
                        STOP
                     ENDIF
                  ENDIF

                  ! Where the parameterized OH is negative, 
                  ! set equal to small number.  
                  ! Sometimes for low OH, the parameterized OH 
                  ! is actually negative.

                   IF (self%BOH(I,J,L).LT.0.) THEN
		   	self%BOH(I,J,L)=1000.
!		   	PRINT *, 'runOH: OH < 0 reset OH to 1000'
		   ENDIF
	       ENDIF

	       
!
 300           CONTINUE
!
 
            ENDDO  !DO I=1,MLONBX
         ENDDO  !DO J=1,MLATBX
      ENDDO  !DO L=1,MVRTBX

      
      IF (MAPL_AM_I_ROOT()) THEN
         print*,'********************************************'
         print*,'Leaving OH param calculation in SR OHparam.'
         print*,'********************************************'
      END IF

      RETURN

      END SUBROUTINE runOHparam
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: finalizeOHparam
!
! !INTERFACE:
!
      SUBROUTINE finalizeOHparam(self)
!
! !INPUT/OUTPUT PARAMETERS:
      type (t_OHparam),   intent(inOut) :: self
!
! !DESCRIPTION:
!
!EOP
!------------------------------------------------------------------------------
!BOP
      deallocate(self%BOH    )
      deallocate(self%dustrat)
      deallocate(self%INDVARA)
      deallocate(self%INDVARB)
      deallocate(self%INDVARC)
      deallocate(self%INDVARD)

      return

      END SUBROUTINE finalizeOHparam
!EOC
!------------------------------------------------------------------------------


end module ohParameterizationMethod_mod
