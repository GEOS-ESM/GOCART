#include "MAPL_Generic.h"

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3
!
!  June 2017 - Manyin: Only read in Root process, then Broadcast to all
!
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  inputFileReading_mod ---
!
! !INTERFACE:
!

   MODULE  inputFileReading_mod

! !USES:

   USE ESMF
   USE MAPL


   USE cblock_OH_mod    ! "CMN_OH"
   USE cblock_cf_mod    ! "CMN_CF"
   USE cblock_CO_mod    ! "CMN_CO"

   IMPLICIT NONE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PRIVATE
      PUBLIC  :: readavgOH
      PUBLIC  :: READ_COEFF
      PUBLIC  :: read_correction
      PUBLIC  :: read_dustrat
      PUBLIC  :: read_srefl

  include "mpif.h"

!**********************************************************************
CONTAINS
!**********************************************************************
      SUBROUTINE readavgOH(unitNum, fileName)
!**********************************************************************
      IMPLICIT NONE
 
      INTEGER :: unitNum
      CHARACTER (len=*) :: fileName

      CHARACTER(LEN=*), PARAMETER :: Iam = 'OH:readavgOH'
      integer        :: status, rc
      type (ESMF_VM) :: VM
      integer        :: mpi_comm
      integer        :: rootID

      INTEGER I,J,K
!**********************************************************************
! Created by Bryan Duncan.
!**********************************************************************
! In subdomains where OH is very low or negligible due to low 
!  sunlight (e.g., high latitudes in winter), concentrations of OH are 
!  set to climatological mean values as a function of latitude,
!  altitude and season. This SR reads in the average OH
!  fields from Spivakovsky et al., "Three-dimensional climatological 
!  distribution of tropospheric OH: update and evaluation", accepted
!  to JGR, 1999.
!**********************************************************************
! List of Variables & Arrays
!
! avgOH = array containing the climatological OH values. 
!
! NCMSALTS = number of altitude levels of climatology.
!
! NCMSLATS = number of latitude bands of climatology.
!
! NSEAS    = number of seasons of climatlogy.
!
! NGENERICFILE   = filenumber specified by user in CMN_OH.
!**********************************************************************
! 
!      print*,'Reading in OH climatology of cms in readavgOH.'
!
      call ESMF_VMGetCurrent(VM, __RC__)
      call ESMF_VMGet(VM, mpiCommunicator=mpi_comm, __RC__)
      rootID = 0

      IF(MAPL_AM_I_ROOT()) THEN
!!!

      OPEN(UNIT=unitNum, FILE=TRIM(fileName),STATUS='OLD')
!
!      unitNum = NGENERICFILE
!      OPEN(UNIT=NGENERICFILE,FILE='avgOH',STATUS='OLD')
!
      READ(unitNum,2)
 2    FORMAT(//)
!
      DO I=1,NSEAS
!
         READ(unitNum,3)
 3       FORMAT(////)
!
         DO J=1,NCMSLATS
            READ(unitNum,4)(avgOH(I,J,K),K=1,NCMSALTS)
         ENDDO
!
      ENDDO
!
 4    FORMAT(5x,7(F6.2))
!

      CLOSE(unitNum)
!!!
      ENDIF

! broadcast to all processors
      CALL MPI_BCAST( avgOH, NSEAS*NCMSLATS*NCMSALTS, MPI_REAL8, rootID, mpi_comm, status )
      VERIFY_(status)

!
!**********************************************************************
      RETURN

      END SUBROUTINE readavgOH
! **********************************************************************
      SUBROUTINE READ_COEFF (unitNum, fileName)
! **********************************************************************
      IMPLICIT NONE

      INTEGER :: unitNum
      CHARACTER (len=*) :: fileName

      CHARACTER(LEN=*), PARAMETER :: Iam = 'OH:READ_COEFF'
      integer        :: status, rc
      type (ESMF_VM) :: VM
      integer        :: mpi_comm
      integer        :: rootID

      INTEGER JTERM,PARSUM,NFT
      INTEGER J,L,K,KK,LH,KR,IDONE,M
! **********************************************************************
! Created by Bryan Duncan.
! **********************************************************************
! This SR reads in information from the file 'start_information'. 
! The information is needed to calculate the OH using parameterizations
! and includes polynomial coefficients, ranges of independent variables,
! etc.  See the text of this code for descriptions of information. 
! *********************************************************************

      call ESMF_VMGetCurrent(VM, __RC__)
      call ESMF_VMGet(VM, mpiCommunicator=mpi_comm, __RC__)
      rootID = 0

      IF(MAPL_AM_I_ROOT()) THEN
!!!
!      print*,'Reading info for parameterization in SR read_coeff.'
!
! NGENERICFILE = filenumber specified by user in CMN_OH.
!
      OPEN(UNIT=unitNum, FILE=TRIM(fileName),STATUS='OLD')

!      unitNum = NGENERICFILE
!      OPEN(UNIT=NGENERICFILE,FILE='start_information',STATUS='OLD')
      REWIND(unitNum)
!
!*********************************************************************
! Loop over all NPARAM parameterizations.  There are over 200 individual
! parameterizations that are used to describe the tropospheric OH
! field.
!
! NPARAM  = total number of individual parameterizations.
!*********************************************************************
!
      DO M=1,NPARAM
!
!*********************************************************************
! The 24-hour average concentration of OH is represented as a set 
!  of high-order polynomials in variables such as temperature, latitude,
!  declination and the concentrations of ozone, water vapor, carbon
!  monoxide, nitrogen oxides (as a family), and hydrocarbons.  (See
!  SR GETINFO for a complete list of independent variables.)
!
! MOVAR = the total number of independent variables.
!         
      READ(unitNum,*) MOVAR(M)
!      IF (MAPL_AM_I_ROOT()) PRINT *, ' MOVAR:  : ', MOVAR(M)
!
! Loop over MOVAR independent variables reading in parameter ranges.
! These ranges for chemical species and physical parameters were estimated from
! observations and monthly average values predicted by a global model. The global
! model used for estimating ranges was the Harvard-GEOS model of transport and
! chemistry driven by assimilated meteorological data from the Goddard Earth
! Observing System Data Assimilation Office (GEOSDAO) (Bey et al. [1999]). 
!
! Bey, I., R. Yantosca and D. Jacob:  Export of Pollutants from Eastern Asia: A
! Simulation of the PEM-West (B) Aircraft Mission Using a 3-D Model Driven by
! Assimilated Meteorological Fields, Presentation at American Geophysical Union
! Spring Meeting, Boston, 1999.
!
! RANGEM = ranges for chemical species and physical parameters
!
        DO J=1,MOVAR(M)
!
          READ(unitNum,*)RANGEM(M,J,1),RANGEM(M,J,2)
!	  IF (MAPL_AM_I_ROOT()) PRINT *, ' RANGE : ',RANGEM(M,J,1),RANGEM(M,J,2)
!
        ENDDO
!
! IDENTOLD,ELTODO,IROWST,L = information of domain divisions. 
!
! NDONE = number of polynomials used in one parameterization. If there are
!         no domain divisions, NDONE=1.
!
      READ(unitNum,*) L
!      IF (MAPL_AM_I_ROOT()) PRINT *, ' L : ',L

      DO K=1,L 
        READ(unitNum,*) IROWST(M,K,0), &
     &       (IROWST(M,K,KR),KR=1,IROWST(M,K,0))
!      IF (MAPL_AM_I_ROOT()) PRINT *, 'IROWST: ', &
!         IROWST(M,K,0), (IROWST(M,K,KR),KR=1,IROWST(M,K,0))
      ENDDO
!
        READ(unitNum,*) KK
!	 IF (MAPL_AM_I_ROOT()) PRINT *, 'KK: ', KK
 
      DO LH=1,KK
        READ(unitNum,*) (ELTODO(M,LH,J),J=1,MOVAR(M)+2)
!	IF (MAPL_AM_I_ROOT()) PRINT *, 'ELTODO ', (ELTODO(M,LH,J),J=1,MOVAR(M)+2)
      ENDDO
!
        READ(unitNum,*)NDONE(M),(IDENTOLD(M,J),J=1,NDONE(M))
!	IF (MAPL_AM_I_ROOT()) PRINT *, 'NDONE: ' , NDONE(M),(IDENTOLD(M,J),J=1,NDONE(M))
!
! *********************************************
! Error Check.
       IF(NDONE(M).GT.MXDONE) THEN
         PRINT*,'SR READ_COEFF: NDONE > MXDONE'
         PRINT*,'NDONE(M)=',NDONE(M)
         PRINT*,'MXDONE  =',MXDONE
         PRINT*,'Increase MXDONE in CMN_OH to be greater than NDONE.'
         STOP
       ENDIF 
! *********************************************
! Read in number of terms and coefficients used in each polynomial.
!
! NENDW = number of terms in polynomial.
!
! COEFF = polynomial coefficients.
!
!	IF (MAPL_AM_I_ROOT()) PRINT *,'Print polymonial coeffs:'
       DO IDONE=1,NDONE(M)
        READ(unitNum,*) NFT,NENDW(M,IDONE)
!	IF (MAPL_AM_I_ROOT()) PRINT *, NFT,NENDW(M,IDONE)
          DO JTERM=1,NENDW(M,IDONE)
            READ(unitNum,*) COEFF(M,IDONE,JTERM)
!	    IF (MAPL_AM_I_ROOT()) PRINT *, COEFF(M,IDONE,JTERM)
          ENDDO
       ENDDO 
!
! ACCUR = RMS error associated with individual parameterization.
!
    
      READ(unitNum,*) ACCUR(M)
!     IF (MAPL_AM_I_ROOT()) PRINT *,' Accuracy: ',ACCUR(M)
!
!*********************************************************************
      ENDDO
!*********************************************************************
!
      CLOSE(unitNum)
!!!
      ENDIF

! broadcast to all processors
      CALL MPI_BCAST( MOVAR,    NPARAM,                 MPI_INTEGER,  rootID, mpi_comm, status )
      VERIFY_(status)
      CALL MPI_BCAST( RANGEM,   NPARAM*MXVAR*2,         MPI_REAL8,    rootID, mpi_comm, status )
      VERIFY_(status)
      CALL MPI_BCAST( IDENTOLD, NPARAM*MXDONE,          MPI_INTEGER,  rootID, mpi_comm, status )
      VERIFY_(status)
      CALL MPI_BCAST( ELTODO,   NPARAM*MXROW*MXCLM,     MPI_REAL8,    rootID, mpi_comm, status )
      VERIFY_(status)
      CALL MPI_BCAST( IROWST,   NPARAM*MXELM*(MXROW+1), MPI_INTEGER,  rootID, mpi_comm, status )
      VERIFY_(status)
      CALL MPI_BCAST( NDONE,    NPARAM,                 MPI_INTEGER,  rootID, mpi_comm, status )
      VERIFY_(status)
      CALL MPI_BCAST( NENDW,    NPARAM*MXDONE,          MPI_INTEGER,  rootID, mpi_comm, status )
      VERIFY_(status)
      CALL MPI_BCAST( COEFF,    NPARAM*MXDONE*MXCOL,    MPI_REAL8,    rootID, mpi_comm, status )
      VERIFY_(status)

!
!*********************************************************************
!
! Sum up the total number of individual parameterizations in PARSUM.  
! Since each parameterization may be described by more than one
! polynomial, the total number of polynomials are summed in the
! variables NDFUNCS*. NDFUNCS* are used only for bookkeeping. 
! See SR GETOH for a description of variables NPARAMA-F.
!
!*********************************************************************
        PARSUM=0
!*********************************************************************

        NDFUNCSB1=0
      DO J=1,NPARAMB
        NDFUNCSB1=NDFUNCSB1+NDONE(J)
      ENDDO

!*********************************************************************
        PARSUM=PARSUM+NPARAMB
!*********************************************************************

        NDFUNCSB2=0
      DO J=PARSUM+1,PARSUM+NPARAMB
        NDFUNCSB2=NDFUNCSB2+NDONE(J)
      ENDDO

!*********************************************************************
       PARSUM=PARSUM+NPARAMB
!*********************************************************************

        NDFUNCSB3=0
      DO J=PARSUM+1,PARSUM+NPARAMB
        NDFUNCSB3=NDFUNCSB3+NDONE(J)
      ENDDO

!*********************************************************************
       PARSUM=PARSUM+NPARAMB
!*********************************************************************
 
        NDFUNCSA1=0
      DO J=PARSUM+1,PARSUM+NPARAMA
        NDFUNCSA1=NDFUNCSA1+NDONE(J)       
      ENDDO

!*********************************************************************
       PARSUM=PARSUM+NPARAMA
!*********************************************************************

        NDFUNCSA2=0
      DO J=PARSUM+1,PARSUM+NPARAMA
        NDFUNCSA2=NDFUNCSA2+NDONE(J)       
      ENDDO

!*********************************************************************
       PARSUM=PARSUM+NPARAMA
!*********************************************************************

        NDFUNCSA3=0
      DO J=PARSUM+1,PARSUM+NPARAMA
        NDFUNCSA3=NDFUNCSA3+NDONE(J)       
      ENDDO

!*********************************************************************
       PARSUM=PARSUM+NPARAMA
!*********************************************************************

        NDFUNCSC1=0
      DO J=PARSUM+1,PARSUM+NPARAMC
        NDFUNCSC1=NDFUNCSC1+NDONE(J)       
      ENDDO

!*********************************************************************
       PARSUM=PARSUM+NPARAMC
!*********************************************************************

        NDFUNCSC2=0
      DO J=PARSUM+1,PARSUM+NPARAMC
        NDFUNCSC2=NDFUNCSC2+NDONE(J)       
      ENDDO

!*********************************************************************
       PARSUM=PARSUM+NPARAMC
!*********************************************************************

       NDFUNCSC3=0
      DO J=PARSUM+1,PARSUM+NPARAMC
        NDFUNCSC3=NDFUNCSC3+NDONE(J)
      ENDDO

!*********************************************************************
       PARSUM=PARSUM+NPARAMC
!*********************************************************************

        NDFUNCSD1=0
      DO J=PARSUM+1,PARSUM+NPARAMD
        NDFUNCSD1=NDFUNCSD1+NDONE(J)       
      ENDDO

!*********************************************************************
       PARSUM=PARSUM+NPARAMD
!*********************************************************************

        NDFUNCSE1=0
      DO J=PARSUM+1,PARSUM+NPARAME
        NDFUNCSE1=NDFUNCSE1+NDONE(J)       
      ENDDO

!*********************************************************************
       PARSUM=PARSUM+NPARAME
!*********************************************************************

        NDFUNCSF1=0
      DO J=PARSUM+1,PARSUM+NPARAMF
        NDFUNCSF1=NDFUNCSF1+NDONE(J)       
      ENDDO

!*********************************************************************
       PARSUM=PARSUM+NPARAMF
!*********************************************************************

! Error Check.
       IF(PARSUM.NE.NPARAM) THEN
          PRINT*,'SR READ_COEFF: PARSUM.NE.NPARAM!'
          PRINT*,'This means that NPARAM specified in CMN_OH does'
          PRINT*,'not equal NPARAM (i.e., PARSUM here) calculated'
          PRINT*,'in this SR.  Check to see why they differ and'
          PRINT*,'correct the problem.'
          STOP
       ENDIF
!
!*********************************************************************
 300  RETURN 
      END SUBROUTINE READ_COEFF
!**********************************************************************
      SUBROUTINE read_correction (unitNum, fileName)
!**********************************************************************
      IMPLICIT NONE

      INTEGER :: unitNum
      CHARACTER (len=*) :: fileName

      CHARACTER(LEN=*), PARAMETER :: Iam = 'OH:read_correction'
      integer        :: status, rc
      type (ESMF_VM) :: VM
      integer        :: mpi_comm
      integer        :: rootID

      INTEGER I,J,K
!**********************************************************************
! Created by Bryan Duncan.
!**********************************************************************
! This SR reads in correction factors to remove the problem with high
! background aerosol (tau=0.1) in the OH parameterization.
! This file contains ratios of OH calculated assuming a uniform
! tau=0.01 over the globe and assuming tau=0.1. OH calculated by CMS, June, 2001.
! Ratio created by bnd, June, 2001.
!*************************************************************************
! List of Variables & Arrays
!
! correction = correction factor 
!
! NCMSALTS2 = number of altitude levels.
!
! NCMSLATS2 = number of latitude bands.
!
! NSEAS_CF  = number of seasons.
!
! NFILENUM  = filenumber specified by user in CMN_OH.
!**********************************************************************
!      print*,'Reading in OH correction factors in SR read_correction.'
!
      call ESMF_VMGetCurrent(VM, __RC__)
      call ESMF_VMGet(VM, mpiCommunicator=mpi_comm, __RC__)
      rootID = 0

      IF(MAPL_AM_I_ROOT()) THEN
!!!
      OPEN(UNIT=unitNum, FILE=TRIM(fileName),STATUS='OLD')

!      unitNum = NFILENUM
!      OPEN(UNIT=NFILENUM,FILE='CorrectRatio',STATUS='OLD')
!
      correction(:,:,:)=1.
!
      READ(unitNum,2)
 2    FORMAT(//////////////)
!
      DO I=1,NSEAS_CF
!
         READ(unitNum,3)
 3       FORMAT()
!
         DO J=1,NCMSLATS2
            READ(unitNum,4)(correction(I,J,K),K=1,NCMSALTS2)
         ENDDO
!
      ENDDO
!
 4    FORMAT(4x,9(f5.3,2x)) 
!
      CLOSE(unitNum)
!
      DO I=1,NSEAS_CF
         DO J=1,NCMSLATS2
            DO K=1,NCMSALTS2
               IF (correction(I,J,K).LT.0..OR.correction(I,J,K).GT.2.) THEN
                  PRINT*,'Correction factor < 0!'
                  PRINT*,I,J,K,correction(I,J,K)
                  PRINT*,'Stopped in SR READ_CORRECTION.'
                  STOP
              ENDIF
            ENDDO
         ENDDO
      ENDDO
!!!
      ENDIF

! broadcast to all processors
      CALL MPI_BCAST( correction, NSEAS_CF*NCMSLATS2*NCMSALTS2, MPI_REAL8, rootID, mpi_comm, status )
      VERIFY_(status)
!
!      print*,'Leaving SR read_correction.'

!**********************************************************************
      RETURN

      END SUBROUTINE read_correction
!**********************************************************************
      SUBROUTINE read_dustrat (unitNum, fileName, dustrat, MLONBX, MLATBX)
!**********************************************************************
      IMPLICIT NONE

      INTEGER :: unitNum
      INTEGER :: MLONBX, MLATBX
      REAL*8  :: dustrat (MLONBX, MLATBX, 18, 12)
      CHARACTER (len=*) :: fileName

      CHARACTER(LEN=*), PARAMETER :: Iam = 'OH:read_dustrat'
      integer        :: status, rc
      type (ESMF_VM) :: VM
      integer        :: mpi_comm
      integer        :: rootID

      INTEGER I,J,L,K
!****************************************************
! Created by Bryan Duncan.
!****************************************************
! This SR reads in information needed for SR AerosolOH.
!****************************************************
! List of Variables & Arrays
!
!  = array holding the percent reduction in OH
!     from an atmosphere with and without aerosols
!
!****************************************************
!      print*,'Reading in info for OH reduction by aerosols'
!      print*,'in SR read_percentages.'
!
!bnd      OPEN(NGENERICFILE,FILE='oh.dust.ratio',STATUS='OLD') 

      call ESMF_VMGetCurrent(VM, __RC__)
      call ESMF_VMGet(VM, mpiCommunicator=mpi_comm, __RC__)
      rootID = 0

      IF(MAPL_AM_I_ROOT()) THEN
!!!
      OPEN(UNIT=unitNum, FILE=TRIM(fileName),STATUS='OLD')

      ! unitNum = NGENERICFILE
      ! OPEN(NGENERICFILE,FILE='oh.allaer.ratio',STATUS='OLD') 
!
      dustrat(:,:,:,:)=1.
!

      DO K=1,12
         DO I=1,MLONBX
            DO J=1,MLATBX
               DO L=1,numStratLevels
                  READ(unitNum,*) dustrat(I,J,L,K)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
!
      CLOSE(unitNum)
!!!
      ENDIF

! broadcast to all processors
      CALL MPI_BCAST( dustrat, MLONBX*MLATBX*numStratLevels*12, MPI_REAL8, rootID, mpi_comm, status )
      VERIFY_(status)

!      print*,'Leaving SR read_percentages.'
!       
!***************************************************************
!
      RETURN
      END SUBROUTINE read_dustrat
!**********************************************************************
      SUBROUTINE read_srefl (unitNum, fileName)
!**********************************************************************
      IMPLICIT NONE

      INTEGER :: unitNum
      CHARACTER (len=*) :: fileName

      CHARACTER(LEN=*), PARAMETER :: Iam = 'OH:read_srefl'
      integer        :: status, rc
      type (ESMF_VM) :: VM
      integer        :: mpi_comm
      integer        :: rootID

      INTEGER :: J, M
!****************************************************
! Created by Bryan Duncan.
!****************************************************
!*****************************************************************************
! Read in reflectivity table (Srefl.table)
! Table data used to convert user's optical depths to
!   flux fractions - see get_fluxfraction.f.
!*****************************************************************************

      call ESMF_VMGetCurrent(VM, __RC__)
      call ESMF_VMGet(VM, mpiCommunicator=mpi_comm, __RC__)
      rootID = 0

      IF(MAPL_AM_I_ROOT()) THEN
!!!
      OPEN(UNIT=unitNum, FILE=TRIM(fileName),STATUS='OLD')

      !OPEN (UNIT=unitNum,FILE='Srefl.table',STATUS='old')
      REWIND(unitNum)
!     
      READ(unitNum,9)
 9    FORMAT(////)

      DO M=1,NREFL
         READ(unitNum,100) OPTDEPTH(M),(RFLC(M,J),J=1,NREFLCOL)
      END DO
!     
 100     FORMAT(3X,F6.3,5(F8.3))

      CLOSE(unitNum)
!!!
      ENDIF

! broadcast to all processors
      CALL MPI_BCAST( OPTDEPTH, NREFL, MPI_REAL8, rootID, mpi_comm, status )
      VERIFY_(status)

      CALL MPI_BCAST( RFLC, NREFL*NREFLCOL, MPI_REAL8, rootID, mpi_comm, status )
      VERIFY_(status)
!       
!***************************************************************
!
      RETURN
      END SUBROUTINE read_srefl

  END MODULE inputFileReading_mod
