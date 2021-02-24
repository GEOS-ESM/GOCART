#include "MAPL_Generic.h"
!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  CFC_GridCompMod --- CFC Grid Component Class
!
! !INTERFACE:
!

   MODULE  CFC_GridCompMod

! !USES:

   USE ESMF
   USE MAPL
   USE Chem_Mod 	     ! Chemistry Base Class
   USE Chem_StateMod	     ! Chemistry State
   USE Chem_ConstMod, ONLY: grav
   USE Chem_UtilMod	     ! I/O
   USE m_inpak90	     ! Resource file management

   IMPLICIT NONE

! !PUBLIC TYPES:
!
   PRIVATE
   PUBLIC  CFC_GridComp       ! The CFC object 

   include "netcdf.inc"

!
! !PUBLIC MEMBER FUNCTIONS:
!

   PUBLIC  CFC_GridCompSetServices
   PUBLIC  CFC_GridCompInitialize
   PUBLIC  CFC_GridCompRun
   PUBLIC  CFC_GridCompFinalize

!
! !DESCRIPTION:
!
!  This module implements the CFC Grid Component. 
!
! !REVISION HISTORY:
!
!  16Sep2003 da Silva  First crack.
!  01Aug2006 da Silva  Extensions for GEOS-5.
!  12Feb2005  Nielsen  8 regions for INTEX-B 2006
!   1Jan2008  Nielsen  CFC-12 configuration for ARCTAS.
!   8Feb2008  Nielsen  Standard configuration call(s) from AeroChem.
!   1Nov2012  Nielsen  Accomodate cubed sphere for GEOS-5 Ganymed releases
!  27Jun2014  Nielsen  Added CFC-12 photorate to the export state
!
!EOP
!-------------------------------------------------------------------------

  TYPE CFC_GridComp

    CHARACTER(LEN=255) :: name

! For CFC-12 photolysis
! ---------------------
    INTEGER :: nlam
    INTEGER :: nsza
    INTEGER :: numo3
    INTEGER :: nx
    INTEGER :: nxdo
    INTEGER :: nts
    INTEGER :: photEquNumber

    REAL, POINTER :: sdat(:,:,:,:)
    REAL, POINTER :: sza_tab(:)
    REAL, POINTER :: o3_tab(:,:)
    REAL, POINTER :: xtab(:,:,:)

    REAL, POINTER :: CFCsfcFlux(:,:)	 ! CFC-12 surface flux kg m^-2 s^-1
    REAL, POINTER :: CFCloss(:,:,:,:)	 ! CFC loss due to photolysis m^-3 s^-1

    LOGICAL :: DebugIsOn

  END TYPE CFC_GridComp

CONTAINS

   subroutine CFC_GridCompSetServices(  gc, chemReg, rc)
   type(ESMF_GridComp), intent(INOUT) :: GC
   type(Chem_Registry), intent(INOUT) :: chemReg
   integer,             intent(OUT  ) :: rc

   CHARACTER(LEN=255) :: rcbasen = 'CFC_GridComp'
   CHARACTER(LEN=255) :: name

   integer            :: status
   character(len=ESMF_MAXSTR) :: Iam

   Iam = "CFC_GridCompSetServices"

   call MAPL_AddImportSpec(GC,           &
        SHORT_NAME = 'CFC12',            &
        LONG_NAME  = 'CFC 12 Emissions', &
        UNITS      = 'kg s-1 m-2',       &
        DIMS       = MAPL_DimsHorzOnly,  &
        VLOCATION  = MAPL_VLocationNone, &
        RESTART    = MAPL_RestartSkip,   &
        RC         = STATUS)
   VERIFY_(STATUS)

   RETURN_(ESMF_SUCCESS)
   end subroutine CFC_GridCompSetServices

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  CFC_GridCompInitialize --- Initialize CFC_GridComp
!
! !INTERFACE:
!

   SUBROUTINE CFC_GridCompInitialize( gcCFC, w_c, impChem, expChem, &
                                      nymd, nhms, cdt, rc )

! !USES:

  IMPLICIT NONE

! !INPUT PARAMETERS:

   TYPE(Chem_Bundle), intent(in) :: w_c        ! Chemical tracer fields      
   INTEGER, INTENT(IN) :: nymd, nhms	       ! time
   REAL,    INTENT(IN) :: cdt		       ! chemical timestep (secs)


! !OUTPUT PARAMETERS:

   TYPE(CFC_GridComp), INTENT(INOUT) :: gcCFC   ! Grid Component
   TYPE(ESMF_State), INTENT(INOUT)  :: impChem  ! Import State
   TYPE(ESMF_State), INTENT(INOUT)  :: expChem  ! Export State
   INTEGER, INTENT(OUT) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 - 

! !DESCRIPTION: Initializes the CFC Grid Component. It primarily sets
!               the import state for each active constituent package.
!
! !REVISION HISTORY:
!
!  18Sep2003 da Silva  First crack.
!  31May2005  Nielsen  Mods for 7 CO bins, 5 region masks
!  04Nov2005     Bian  CO tagged to 4 regions 
!                      (global, North America, South America, and Africa)
!                      for CR-AVE
!  12Feb2005  Nielsen  8 regions for INTEX-B 2006
!   1Jan2008  Nielsen  CFC-12 configuration for ARCTAS
!   1Nov2012  Nielsen  Accomodate cubed sphere for GEOS-5 Ganymed releases
!
!EOP
!-------------------------------------------------------------------------
#include "mpif.h"

   CHARACTER(LEN=*), PARAMETER :: Iam = 'CFC_GridCompInitialize'
   TYPE(ESMF_VM) :: vm

   CHARACTER(LEN=255) :: rcfilen = 'CFC_GridComp.rc'

   CHARACTER(LEN=255) :: fnPhoto, fileName

   REAL :: x
   REAL, ALLOCATABLE :: w(:)

   INTEGER :: i, i1, i2, im, j, j1, j2, jm, k, km, kr, n, nbins, status

   gcCFC%name = 'CFC-12 Chemistry for ARCTAS'

!  Initialize local variables
!  --------------------------
   rc = 0
   i1 = w_c%grid%i1
   i2 = w_c%grid%i2
   im = w_c%grid%im

   j1 = w_c%grid%j1
   j2 = w_c%grid%j2
   jm = w_c%grid%jm

   km = w_c%grid%km

   nbins = w_c%reg%n_CFC

! Grab the virtual machine
! ------------------------
   CALL ESMF_VMGetCurrent(vm, RC=status)
   VERIFY_(status)

! Load resource file
! ------------------
   CALL I90_loadf ( TRIM(rcfilen), status )
   VERIFY_(status)

   CALL I90_label ( 'photolysisFile:', status )
   VERIFY_(status)
   CALL I90_Gtoken( fnPhoto, status )
   VERIFY_(status)

   CALL I90_Label ( 'phot_Equation_number:', status )
   VERIFY_(status)
   gcCFC%photEquNumber = I90_Gint( status )
   VERIFY_(status)

! Run-time debug switch
! ---------------------
   CALL I90_label ( 'DEBUG:', status )
   VERIFY_(status)
   n = I90_gint ( status )
   VERIFY_(status)
   IF(n /= 0) THEN
    gcCFC%DebugIsOn = .TRUE.
   ELSE
    gcCFC%DebugIsOn = .FALSE.
   END IF

! Allocate space
! --------------
   ALLOCATE(gcCFC%CFCsfcFlux(i1:i2,j1:j2), STAT=status )
   VERIFY_(status)
   ALLOCATE(gcCFC%CFCloss(i1:i2,j1:j2,1:km,nbins), STAT=status )
   VERIFY_(status)

! Photolysis tables: Initialize from NetCDF file
! ----------------------------------------------
   fileName = TRIM(fnPhoto)
   CALL readPhotTables(fileName, status)
   VERIFY_(status)

! Reverse vertical ordering of the radiative
! source function and the overhead O3 reference
! ---------------------------------------------
   DO n = 1,gcCFC%nlam
    DO j = 1,gcCFC%numo3
     DO i = 1,gcCFC%nsza
      DO k = 1,km/2
       kr = km-k+1
       x = gcCFC%sdat(i,j,k,n)
       gcCFC%sdat(i,j,k,n) = gcCFC%sdat(i,j,kr,n)
       gcCFC%sdat(i,j,kr,n) = x
      END DO
     END DO
    END DO
   END DO

   ALLOCATE(w(gcCFC%numo3), STAT=status)
   VERIFY_(status)

   DO k = 1,km/2
    kr = km-k+1
    w(1:gcCFC%numo3) = gcCFC%o3_tab(1:gcCFC%numo3,k)
    gcCFC%o3_tab(1:gcCFC%numo3,k) = gcCFC%o3_tab(1:gcCFC%numo3,kr)
    gcCFC%o3_tab(1:gcCFC%numo3,kr) = w(1:gcCFC%numo3)
   END DO

   DEALLOCATE(w, STAT=status)
   VERIFY_(status)

   RETURN
  CONTAINS
!-------------------------------------------------------------------------
!NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS!
!-------------------------------------------------------------------------
!BOP
!
! !INTERFACE:

 SUBROUTINE readPhotTables(fileName, rc)

! !USES:

  IMPLICIT NONE

! !INPUT PARAMETERS:
!
  CHARACTER(LEN=*), INTENT(IN) :: fileName
!
! !OUTPUT PARAMETERS:
!
  INTEGER, INTENT(OUT) :: rc

! !DESCRIPTION:
!
! Read tables for photolysis in StratChem ... from a NetCDF file
!
! Restrictions:
!  ASSERT that the number of pressure layers in the dataset equals km.
!
! !REVISION HISTORY:
!  Nielsen     11 May 2012: First crack.
!
!EOP
!-----------------------------------------------------------------------

  CHARACTER(LEN=ESMF_MAXSTR) :: Iam = "CFC::readPhotTables"

  INTEGER :: comm, info, unit, status
  INTEGER :: dimid, i, n

  INTEGER :: length

  INTEGER, PARAMETER :: nD = 6
  CHARACTER(LEN=ESMF_MAXSTR) :: dimName(nD)= (/"nsza  ", &
             "numO3 ", "layers", "nlam  ", "nts   ", "nxdo  " /)

  INTEGER, PARAMETER :: nV = 4
  CHARACTER(LEN=ESMF_MAXSTR) :: varName(nV)= (/"sza  ", &
                           "O3TAB",  "SDAT ",  "XTAB " /)

  rc = 0

  CALL ESMF_VMGet(vm, MPICOMMUNICATOR=comm, rc=status)
  VERIFY_(status)

#undef H5_HAVE_PARALLEL
#ifdef H5_HAVE_PARALLEL

  CALL MPI_Info_create(info, status)
  CALL MPI_Info_set(info, "romio_cb_read", "automatic", status)
  VERIFY_(status)

#ifdef NETCDF_NEED_NF_MPIIO
  status = NF_OPEN_PAR(TRIM(fileName), IOR(NF_NOWRITE,NF_MPIIO), comm, info, unit)
#else
  status = NF_OPEN_PAR(TRIM(fileName), NF_NOWRITE, comm, info, unit)
#endif

#else

  IF(MAPL_AM_I_ROOT(vm)) THEN 
   status = NF_OPEN(TRIM(fileName), NF_NOWRITE, unit)

#endif

   IF(status /= NF_NOERR) THEN
    PRINT *,'Error opening file ',TRIM(fileName), status
    PRINT *, NF_STRERROR(status)
    VERIFY_(status)
   END IF

   DO i = 1,nD

    status = NF_INQ_DIMID(unit, TRIM(dimName(i)), dimid)
    IF(status /= NF_NOERR) THEN
     PRINT *,"Error inquiring dimension ID for ", TRIM(dimName(i)), status
     PRINT *, NF_STRERROR(status)
     VERIFY_(status)
    END IF

    status = NF_INQ_DIMLEN(unit, dimid, n)
    IF(status /= NF_NOERR) THEN
     PRINT *,"Error inquiring  dimension length for ", TRIM(dimName(i)), status
     PRINT *, NF_STRERROR(status)
    END IF

    SELECT CASE (i)
     CASE (1)
      gcCFC%nsza = n
     CASE (2)
      gcCFC%numO3 = n
     CASE (3)
      _ASSERT(n == km,'needs informative message')
     CASE (4)
      gcCFC%nlam = n
     CASE (5)
      gcCFC%nts = n
     CASE (6)
      gcCFC%nxdo = n
     CASE DEFAULT
    END SELECT

   END DO

#ifndef H5_HAVE_PARALLEL

  END IF ! MAPL_AM_I_ROOT

  CALL MAPL_CommsBcast(vm, gcCFC%nsza, 1, 0, RC=status)
  VERIFY_(status)
  CALL MAPL_CommsBcast(vm, gcCFC%numO3, 1, 0, RC=status)
  VERIFY_(status)
  CALL MAPL_CommsBcast(vm, gcCFC%nlam, 1, 0, RC=status)
  VERIFY_(status)
  CALL MAPL_CommsBcast(vm, gcCFC%nts, 1, 0, RC=status)
  VERIFY_(status)
  CALL MAPL_CommsBcast(vm, gcCFC%nxdo, 1, 0, RC=status)
  VERIFY_(status)

#endif

  ALLOCATE(gcCFC%sdat(gcCFC%nsza,gcCFC%numo3,km,gcCFC%nlam), STAT=status)
  VERIFY_(status)
  ALLOCATE(gcCFC%o3_tab(gcCFC%numo3,km), STAT=status)
  VERIFY_(status)
  ALLOCATE(gcCFC%xtab(gcCFC%nlam,gcCFC%nxdo,gcCFC%nts), STAT=status)
  VERIFY_(status)
  ALLOCATE(gcCFC%sza_tab(gcCFC%nsza), STAT=status)
  VERIFY_(status)

#ifndef H5_HAVE_PARALLEL

  IF(MAPL_AM_I_ROOT(vm)) THEN

#endif

   DO i = 1,nV

    status = NF_INQ_VARID(unit, TRIM(varName(i)), n)
    IF(status /= NF_NOERR) THEN
     PRINT *,"Error getting varid for ", TRIM(varName(i)), status
     PRINT *, NF_STRERROR(status)
     VERIFY_(status)
    END IF

    SELECT CASE (i)
     CASE (1)
      status = NF_GET_VAR_REAL(unit, n, gcCFC%sza_tab)
     CASE (2)
      status = NF_GET_VAR_REAL(unit, n, gcCFC%o3_tab)
     CASE (3)
      status = NF_GET_VAR_REAL(unit, n, gcCFC%sdat)
     CASE (4)
      status = NF_GET_VAR_REAL(unit, n, gcCFC%xtab)
     CASE DEFAULT
    END SELECT

    IF(status /= NF_NOERR) THEN
     PRINT *,"Error getting values for ", TRIM(varName(i)), status
     PRINT *, NF_STRERROR(status)
     VERIFY_(status)
    END IF

   END DO

#ifdef H5_HAVE_PARALLEL

   CALL MPI_Info_free(info, status)
   VERIFY_(status)

#else

   status = NF_CLOSE(unit)
   VERIFY_(status)

  END IF ! MAPL_AM_I_ROOT

  length = SIZE(gcCFC%sza_tab)
  CALL MPI_Bcast(gcCFC%sza_tab, length, MPI_REAL, 0, comm, status)
  VERIFY_(status)

  length = SIZE(gcCFC%o3_tab)
  CALL MPI_Bcast(gcCFC%o3_tab, length, MPI_REAL, 0, comm, status)
  VERIFY_(status)

  length = SIZE(gcCFC%sdat)
  CALL MPI_Bcast(gcCFC%sdat, length, MPI_REAL, 0, comm, status)
  VERIFY_(status)

  length = SIZE(gcCFC%xtab)
  CALL MPI_Bcast(gcCFC%xtab, length, MPI_REAL, 0, comm, status)
  VERIFY_(status)

#endif

  RETURN
 END SUBROUTINE readPhotTables

 END SUBROUTINE CFC_GridCompInitialize

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  CFC_GridCompRun --- The CFC Driver 
!
! !INTERFACE:
!

   SUBROUTINE CFC_GridCompRun( gcCFC, w_c, impChem, expChem, nymd, nhms, &
                               cdt, rc)

! !USES:

  IMPLICIT NONE

! !INPUT/OUTPUT PARAMETERS:

   TYPE(CFC_GridComp), INTENT(INOUT) :: gcCFC   ! Grid Component
   TYPE(Chem_Bundle), INTENT(INOUT) :: w_c	! Chemical tracer fields   

! !INPUT PARAMETERS:

   TYPE(ESMF_State), INTENT(INOUT) :: impChem   ! Import State
   INTEGER, INTENT(IN) :: nymd, nhms	        ! time
   REAL,    INTENT(IN) :: cdt		        ! chemical timestep (secs)

! !OUTPUT PARAMETERS:

   TYPE(ESMF_State), INTENT(INOUT) :: expChem   ! Export State
   INTEGER, INTENT(OUT) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 -
!EOP

   CHARACTER(LEN=*), PARAMETER :: Iam = 'CFC_GridCompRun'

   INTEGER ::  i1, i2, im, iXj, j1, j2, jm, km, status
   INTEGER ::  i, indt, j, k, m, n, nbeg, nbins, nend
   REAL :: o3c, qmin, qmax, r, rg, szan

!  Imports
!  -------
   REAL, POINTER, DIMENSION(:,:,:) ::  T  => null()
   REAL, POINTER, DIMENSION(:,:,:) ::  O3 => null()
   REAL, POINTER, DIMENSION(:,:)   ::  tropp => null()

!  Local Variables
!  ---------------
   REAL, PARAMETER :: badVal=2.00E+05
   REAL, PARAMETER :: grav=9.80
   REAL, PARAMETER :: mwtAir=28.97
   REAL, PARAMETER :: mwtCFC12=120.917
   REAL, PARAMETER :: Nsuba=6.022E+26
   REAL, PARAMETER :: rstar=8.3143E+03
   REAL, PARAMETER :: O3abv80km = 1.10E+15 !m^{-2}

   REAL, ALLOCATABLE :: emit2vmr(:,:)
   REAL, ALLOCATABLE :: tropPa(:,:)
   REAL, ALLOCATABLE :: pPa(:,:,:)
   REAL, ALLOCATABLE :: nd(:,:,:)
   REAL, ALLOCATABLE :: O3Col(:,:,:)
   REAL, ALLOCATABLE :: photoRate(:,:,:)
   REAL, ALLOCATABLE :: s(:,:,:,:)

! Disable the ACG'd CFC_GetPointer___.h for now. [Maybe fix it soon.]
! -------------------------------------------------------------------
#define EXPORT     expChem
#define ptrCFCEM   CFC_emis
#define ptrCFCPH   CFC_phot
#define ptrCFCLS   CFC_loss
#define ptrCFCCL   CFC_column

!JEN#include "CFC_GetPointer___.h"

!  Bin sizes
!  ---------
   integer, parameter		   :: NBIN_CFCEM = 1 ! CFC Emission
   integer, parameter		   :: NBIN_CFCPH = 1 ! CFC Photorate
   integer, parameter		   :: NBIN_CFCLS = 2 ! CFC Loss due to photolysis
   integer, parameter		   :: NBIN_CFCCL = 2 ! CFC Column

!  Bin-indexed Chem Arrays
!  -----------------------
   type(Chem_Array), target	   ::	 CFCEM(NBIN_CFCEM) ! Export: CFC Surface flux
   type(Chem_Array), pointer	   :: ptrCFCEM(:)	   ! Export: CFC Surface flux
   type(Chem_Array), target	   ::	 CFCPH(NBIN_CFCPH) ! Export: CFC Photorate
   type(Chem_Array), pointer	   :: ptrCFCPH(:)	   ! Export: CFC Photorate
   type(Chem_Array), target	   ::	 CFCLS(NBIN_CFCLS) ! Export: CFC Loss due to photolysis
   type(Chem_Array), pointer	   :: ptrCFCLS(:)	   ! Export: CFC Loss due to photolysis
   type(Chem_Array), target	   ::	 CFCCL(NBIN_CFCCL) ! Export: CFC Column
   type(Chem_Array), pointer	   :: ptrCFCCL(:)	   ! Export: CFC Column

!  Local array referencing the Import/Export states
!  ------------------------------------------------
   type(Chem_Array), target	   ::	 CFC12S ! Export: Stratospheric CFC-12 (CCl2F2)
   type(Chem_Array), pointer	   :: ptrCFC12S ! Export: Stratospheric CFC-12 (CCl2F2)
   type(Chem_Array), target	   ::	 CFC12T ! Export: Tropospheric CFC-12 (CCl2F2)
   type(Chem_Array), pointer	   :: ptrCFC12T ! Export: Tropospheric CFC-12 (CCl2F2)
   real, pointer :: ptr2d(:,:) => null()

!  Get pointers to data in state
!  -----------------------------
   call MAPL_GetPointer ( impChem,ptr2d,'CFC12', RC=STATUS)
   VERIFY_(STATUS) 
   gcCFC%CFCsfcFlux = ptr2d

   ptrCFC12S => CFC12S   ! Stratospheric CFC-12 (CCl2F2)
   call MAPL_GetPointer ( EXPORT, CFC12S%data3d,  'CFC12S', RC=STATUS )
   VERIFY_(STATUS)

   ptrCFC12T => CFC12T   ! Tropospheric CFC-12 (CCl2F2)
   call MAPL_GetPointer ( EXPORT, CFC12T%data3d,  'CFC12T', RC=STATUS )
   VERIFY_(STATUS)

   ptrCFCEM => CFCEM	 ! CFC-12 Surface flux
   call MAPL_GetPointer ( EXPORT, CFCEM(1)%data2d,  'CFC12EM', RC=STATUS )
   VERIFY_(STATUS)

   ptrCFCPH => CFCPH	 ! CFC-12 Photorate
   call MAPL_GetPointer ( EXPORT, CFCPH(1)%data3d,  'CFC12PH', RC=STATUS )
   VERIFY_(STATUS)

   ptrCFCLS => CFCLS	 ! CFC-12 Loss due to photolysis
   call MAPL_GetPointer ( EXPORT, CFCLS(1)%data3d,  'CFC12SLS', RC=STATUS )
   VERIFY_(STATUS)
   call MAPL_GetPointer ( EXPORT, CFCLS(2)%data3d,  'CFC12TLS', RC=STATUS )
   VERIFY_(STATUS)

   ptrCFCCL => CFCCL	 ! CFC-12 Column mass density
   call MAPL_GetPointer ( EXPORT, CFCCL(1)%data2d,  'CFC12SCL', RC=STATUS )
   VERIFY_(STATUS)
   call MAPL_GetPointer ( EXPORT, CFCCL(2)%data2d,  'CFC12TCL', RC=STATUS )
   VERIFY_(STATUS)

!  Initialize local variables
!  --------------------------
   rc = 0
   i1 = w_c%grid%i1
   i2 = w_c%grid%i2
   im = w_c%grid%im

   j1 = w_c%grid%j1
   j2 = w_c%grid%j2
   jm = w_c%grid%jm

   km = w_c%grid%km

   iXj = ( i2 - i1 + 1 ) * ( j2 - j1 + 1 )

   nbins = w_c%reg%n_CFC
   nbeg  = w_c%reg%i_CFC
   nend  = w_c%reg%j_CFC

!  Imports
!  -------
   CALL MAPL_GetPointer( impChem,     T,     'T', RC=status )
   VERIFY_(status) 
   CALL MAPL_GetPointer( impChem,    O3,    'O3', RC=status ) 
   VERIFY_(status) 
   CALL MAPL_GetPointer( impChem, tropp, 'TROPP', RC=status ) 
   VERIFY_(status) 

   IF(gcCFC%DebugIsOn) THEN
    CALL pmaxmin('CFC:     T',     T, qmin, qmax, iXj, km, 1. )
    CALL pmaxmin('CFC:    O3',    O3, qmin, qmax, iXj, km, 1. )
    CALL pmaxmin('CFC: TROPP', tropp, qmin, qmax, iXj,  1, 1. )
   END IF

!  Allocate temporary workspace
!  ----------------------------
   ALLOCATE(    emit2vmr(i1:i2,j1:j2), STAT=status)
   VERIFY_(status) 
   ALLOCATE(      tropPa(i1:i2,j1:j2), STAT=status)
   VERIFY_(status) 
   ALLOCATE(      pPa(i1:i2,j1:j2,km), STAT=status)
   VERIFY_(status) 
   ALLOCATE(       nd(i1:i2,j1:j2,km), STAT=status)
   VERIFY_(status) 
   ALLOCATE(    O3Col(i1:i2,j1:j2,km), STAT=status)
   VERIFY_(status) 
   ALLOCATE(photoRate(i1:i2,j1:j2,km), STAT=status)
   VERIFY_(status) 

!  Fix bad tropopause pressure values if they exist.
!  -------------------------------------------------
   CALL Chem_UtilTroppFixer(i2, j2, tropp, VERBOSE=.TRUE., NEWTROPP=tropPa, RC=status)
   VERIFY_(status)

!  Find the pressure at mid-layer
!  ------------------------------
   pPa(i1:i2,j1:j2,1) = w_c%grid%ptop + 0.50*w_c%delp(i1:i2,j1:j2,1)
   DO k = 2, km
    pPa(i1:i2,j1:j2,k) = pPa(i1:i2,j1:j2,k-1) + 0.50* &
                         (w_c%delp(i1:i2,j1:j2,k-1)+w_c%delp(i1:i2,j1:j2,k))
   END DO

!  Number density
!  --------------
   nd(i1:i2,j1:j2,1:km)= nsuba*pPa(i1:i2,j1:j2,1:km)/(rstar*T(i1:i2,j1:j2,1:km))

!  Compute the overlying ozone from mole fraction.  Result: m^{-2}
!  ---------------------------------------------------------------
   r = Nsuba*0.50/(mwtAir*grav)
   O3col(i1:i2,j1:j2,1) = O3abv80km + O3(i1:i2,j1:j2,1)*w_c%delp(i1:i2,j1:j2,1)*r
   DO k=2,km
    O3col(i1:i2,j1:j2,k) = O3col(i1:i2,j1:j2,k-1) + &
                             (O3(i1:i2,j1:j2,k-1) * w_c%delp(i1:i2,j1:j2,k-1) + &
                              O3(i1:i2,j1:j2,  k) * w_c%delp(i1:i2,j1:j2,  k))*r
   END DO

!  Enable the conversion from emission [kg CFC m^{-2} s^{-1}] 
!  to an incremental change in the mixing ratio [s^{-1}].
!  ----------------------------------------------------------
   emit2vmr(i1:i2,j1:j2) = mwtAir*grav/(mwtCFC12*w_c%delp(i1:i2,j1:j2,km))

!  Increment mixing ratio in surface layer of tropospheric CFC-12
!  --------------------------------------------------------------
   w_c%qa(nbeg+1)%data3d(i1:i2,j1:j2,km) = w_c%qa(nbeg+1)%data3d(i1:i2,j1:j2,km)+cdt* &
                                           gcCFC%CFCsfcFlux(i1:i2,j1:j2)*emit2vmr(i1:i2,j1:j2)

!  When tropospheric CFC-12 migrates to the stratosphere, reassign it
!  ------------------------------------------------------------------
   DO k = 1, km
    WHERE(pPa(i1:i2,j1:j2,k) < tropPa(i1:i2,j1:j2) .AND. &
          w_c%qa(nbeg+1)%data3d(i1:i2,j1:j2,k) > 0.00 )
     w_c%qa(nbeg)%data3d(i1:i2,j1:j2,k) = w_c%qa(nbeg)%data3d(i1:i2,j1:j2,k) + &
                                          w_c%qa(nbeg+1)%data3d(i1:i2,j1:j2,k)
     w_c%qa(nbeg+1)%data3d(i1:i2,j1:j2,k) = 0.00
    END WHERE
   END DO

!  Convert CFC-12 to number density
!  --------------------------------
   DO n=nbeg,nend
    w_c%qa(n)%data3d(i1:i2,j1:j2,1:km) = w_c%qa(n)%data3d(i1:i2,j1:j2,1:km)* &
                                         nd(i1:i2,j1:j2,1:km)
   END DO

   ALLOCATE(s(gcCFC%nlam,i1:i2,j1:j2,1:km), STAT=status)
   VERIFY_(status)

!  Photolysis:  Loop over horizontal domain
!  ----------------------------------------
   DO j=j1,j2
    DO i=i1,i2

!  Solar zenith angle (radians).  w_c%cosz has no negative values,
!  which are required for correct interpolation in the S-dat tables.
!  -----------------------------------------------------------------
     IF(w_c%cosz(i,j) <= 1.00E-06) THEN
      szan = ACOS(-0.50)
     ELSE if(w_c%cosz(i,j) < 1.0 ) THEN
      szan = ACOS(w_c%cosz(i,j))
     ELSE
      szan = 0.0
     END IF

     DO k=1,km
      o3c = O3Col(i,j,k)*1.00E-04 !to cm^{-2}

!  Interpolate radiative flux function values.
!  Call getS even when sun is below the horizon.
!  ---------------------------------------------
      CALL getS(k,km,szan,o3c,s(:,i,j,k))
      indt = T(i,j,k)-148.5
      indt = MAX(1,indt)
      indt = MIN(indt,200)

!  Rate constant is sum over wavelengths
!  -------------------------------------
      photoRate(i,j,k) = SUM(s(1:gcCFC%nlam,i,j,k)*gcCFC%xtab(1:gcCFC%nlam,gcCFC%photEquNumber,indt))

     END DO ! Layer
    END DO  ! Longitude
   END DO   ! Latitude

   DEALLOCATE(s, STAT=status)
   VERIFY_(status)
   m = 0

!  Apply photolysis
!  ----------------
   DO n=nbeg,nend
    m = m+1
    w_c%qa(n)%data3d(i1:i2,j1:j2,1:km) = w_c%qa(n)%data3d(i1:i2,j1:j2,1:km) - cdt * &
                                         w_c%qa(n)%data3d(i1:i2,j1:j2,1:km) * &
                                         photoRate(i1:i2,j1:j2,1:km)
    gcCFC%CFCloss(i1:i2,j1:j2,1:km,m) = w_c%qa(n)%data3d(i1:i2,j1:j2,1:km) * &
                                        photoRate(i1:i2,j1:j2,1:km)
   END DO
 
!  Return CFC-12 to mole fraction
!  ------------------------------
   DO n=nbeg,nend
    w_c%qa(n)%data3d(i1:i2,j1:j2,1:km) = w_c%qa(n)%data3d(i1:i2,j1:j2,1:km)/ &
                                         nd(i1:i2,j1:j2,1:km)
   END DO

!  Fill the export states.

!  CFC-12 Surface emission in kg m^{-2} s^{-1}
!  -------------------------------------------
   IF(ASSOCIATED(CFC_emis(1)%data2d)) &
    CFC_emis(1)%data2d(i1:i2,j1:j2) = gcCFC%CFCsfcFlux(i1:i2,j1:j2)

!  CFC-12 photorate s^{-1}
!  -----------------------
   IF(ASSOCIATED(CFC_phot(1)%data3d)) &
    CFC_phot(1)%data3d(i1:i2,j1:j2,1:km) = photoRate(i1:i2,j1:j2,1:km)

!  Loss due to photolysis: Currently m^{-3} s^(-1), and positive for loss.
!  -----------------------------------------------------------------------
   DO n = 1, nbins
    IF(ASSOCIATED(CFC_loss(n)%data3d)) &
     CFC_loss(n)%data3d(i1:i2,j1:j2,1:km) = gcCFC%CFCloss(i1:i2,j1:j2,1:km,n)
   END DO

!  Column burden in kg m(^-2)
!  --------------------------
   DO n = 1, nbins
     IF(ASSOCIATED(CFC_column(n)%data2d)) THEN
      CFC_column(n)%data2d(i1:i2,j1:j2) = 0.
      DO k = 1, km
       CFC_column(n)%data2d(i1:i2,j1:j2) &
        =   CFC_column(n)%data2d(i1:i2,j1:j2) &
          +   w_c%qa(nbeg+n-1)%data3d(i1:i2,j1:j2,k)*mwtCFC12/mwtAir &
            * w_c%delp(i1:i2,j1:j2,k)/grav
     END DO
    END IF
   END DO

!  Clean up
!  --------
   DEALLOCATE( emit2vmr, STAT=status)
   VERIFY_(status)
   DEALLOCATE(   tropPa, STAT=status)
   VERIFY_(status)
   DEALLOCATE(      pPa, STAT=status)
   VERIFY_(status)
   DEALLOCATE(       nd, STAT=status)
   VERIFY_(status)
   DEALLOCATE(    O3Col, STAT=status)
   VERIFY_(status)
   DEALLOCATE(photoRate, STAT=status)
   VERIFY_(status)

   RETURN

   CONTAINS

   SUBROUTINE getS(ik,levels,sza,o3column,s)
! --------------------------------------------------------------------------
! NAME:
!   interp_s
! PURPOSE:
!   Interpolate s values for each wavelength in table to specified O3
!   col and zenith angles
! CATEGORY:
! CALLING SEQUENCE:
!   Call interp_s(nlam,sza,o3column,s)
! INPUTS:
!   nlam     --  number of wavelength intervals used
!   sza      --  zenith angle
!   o3column --  overhead o3 column value 
! OPTIONAL INPUT PARAMETERS:
! OUTPUTS:
!   s -- array of s values (nlam) for each wavelength 
!	 at model p-level interpolated to o3column and sza values
! INTERNAL VARIABLES
!   sza_tab -- values of sza corresponding to sdat table values
!   o3_tab  -- array of overhead O3 values at each p-level (numo3s,np_ctm)
!		used to index sdat
!   sdat    -- input array of values of radiative source function 
!	       (nzens,numo3,np_ctm,nlam) gridded to ctm p layers
! COMMON BLOCKS:
! SIDE EFFECTS:
! PROCEDURE:
!   bi-linear interpolation, for sza>94 s=0, for o3 out of range use min/max
! RESTRICTIONS:
! REQUIRED ROUTINES:
! MODIFICATION HISTORY: 
!   Created 930825 - SR Kawa
!   Modified 960710 for 28 levels and to handle J(O2) separately
!   1Jan2008  Nielsen  CFC-12 configuration for ARCTAS.
!   1Nov2012  Nielsen  Accomodate cubed sphere for GEOS-5 Ganymed releases
! --------------------------------------------------------------------------

   IMPLICIT NONE

   INTEGER, INTENT(IN) :: ik,levels
   REAL, INTENT(IN) :: sza,o3column 
   REAL, INTENT(OUT) :: s(gcCFC%nlam)

   INTEGER :: ijj,ikk,ikkm,il,is
   REAL :: omt,omu,t,u
       
! For each input solar zenith angle, find the first element of
! tabled sza_tab values that is greater than it.  Use this
! table element and previous table element to determined
! interpolated value.
! ------------------------------------------------------------
   DO is=1,gcCFC%nsza
      ijj = is 
      if(gcCFC%sza_tab(is) > sza) EXIT
   END DO
      
! Location is dark, set s/jo2=0        
! -----------------------------
   IF(sza > gcCFC%sza_tab(gcCFC%nsza)) THEN
      s(1:gcCFC%nlam) = 0.
   ELSE
      t = (sza-gcCFC%sza_tab(ijj-1))/(gcCFC%sza_tab(ijj)-gcCFC%sza_tab(ijj-1))
      omt = 1.-t
         
! For each input overhead o3 column find the first element
! of tabled o3_tab values that is > than it.  Use this
! table element and previous table element to determine
! interpolated value
! --------------------------------------------------------
      DO is=1,gcCFC%numo3
  	 ikk = is 
  	 IF (gcCFC%o3_tab(is,ik) > o3column) EXIT
      END DO 

      ikkm = ikk-1

      IF(ikk > 1 .AND. o3column <= gcCFC%o3_tab(gcCFC%numo3,ik)) THEN 
  	 u = (o3column-gcCFC%o3_tab(ikkm,ik))/ &
  	     (gcCFC%o3_tab(ikk,ik)-gcCFC%o3_tab(ikkm,ik))
  	 omu = 1.-u

! Bilinear interpolation at ik for each wavelength
! ------------------------------------------------
    	 DO il=1,gcCFC%nlam	    
    	    s(il) = omt*omu*gcCFC%sdat(ijj-1,ikkm,ik,il) &
    		 +t*omu*gcCFC%sdat(ijj,ikkm,ik,il) &
    		 +t*u*gcCFC%sdat(ijj,ikk,ik,il) &
    		 +omt*u*gcCFC%sdat(ijj-1,ikk,ik,il)
    	 END DO
    
! Extrapolate before table
! ------------------------
      ELSE IF (ikk == 1) THEN
   	 DO il=1,gcCFC%nlam
   	    s(il) = omt*gcCFC%sdat(ijj-1,1,ik,il)+t*gcCFC%sdat(ijj,1,ik,il)
   	 END DO

! Extrapolate past table
! ----------------------
      ELSE 
  	 DO il=1,gcCFC%nlam
  	    s(il) = omt*gcCFC%sdat(ijj-1,gcCFC%numo3,ik,il)+ &
  		    t*gcCFC%sdat(ijj,gcCFC%numo3,ik,il)
  	 END DO 
      END IF  
   END IF
   
   RETURN
   END SUBROUTINE getS

   END SUBROUTINE CFC_GridCompRun

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  CFC_GridCompFinalize
!
! !INTERFACE:
!

   SUBROUTINE CFC_GridCompFinalize( gcCFC, w_c, impChem, expChem, &
                                    nymd, nhms, cdt, rc )

! !USES:

  IMPLICIT NONE

! !INPUT/OUTPUT PARAMETERS:

   TYPE(CFC_GridComp), INTENT(INOUT) :: gcCFC ! Grid Component

! !INPUT PARAMETERS:

   TYPE(Chem_Bundle), INTENT(IN)  :: w_c      ! Chemical tracer fields   
   INTEGER, INTENT(IN) :: nymd, nhms	      ! time
   REAL,    INTENT(IN) :: cdt		      ! chemical timestep (secs)


! !OUTPUT PARAMETERS:

   TYPE(ESMF_State), INTENT(INOUT) :: impChem	! Import State
   TYPE(ESMF_State), INTENT(INOUT) :: expChem	! Import State
   INTEGER, INTENT(OUT) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 -
 
! !DESCRIPTION: This routine finalizes this Grid Component.
!
! !REVISION HISTORY:
!
!  18Sep2003 da Silva  First crack.
!
!EOP
!-------------------------------------------------------------------------

   CHARACTER(LEN=*), PARAMETER :: Iam = 'CFC_GridCompFinalize'
   INTEGER :: status

   rc = 0

   DEALLOCATE(gcCFC%sdat, gcCFC%xtab, gcCFC%o3_tab, gcCFC%sza_tab, &
              gcCFC%CFCloss, gcCFC%CFCsfcFlux, STAT=status )
   VERIFY_(status)

   RETURN
   END SUBROUTINE CFC_GridCompFinalize

 END MODULE CFC_GridCompMod

