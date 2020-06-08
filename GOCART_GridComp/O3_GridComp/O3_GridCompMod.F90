#include "MAPL_Generic.h"
!-------------------------------------------------------------------------
!NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS!
!-------------------------------------------------------------------------
!BOP

! !MODULE:  O3_GridCompMod

! Grid Component class for parameterized Chemistry for ozone:

! !INTERFACE:
!

   MODULE  O3_GridCompMod

! !USES:

   USE ESMF
   USE MAPL

   USE Chem_Mod 	     ! Chemistry Base Class
   USE Chem_StateMod	     ! Chemistry State
   Use Chem_UtilMod, ONLY: pmaxmin             ! Utilities
   USE Chem_UtilMod, ONLY: Chem_UtilTroppFixer ! Fixes bad tropopause pressure values
   USE m_inpak90	     ! Resource file management  

   USE Henrys_law_ConstantsMod, ONLY: get_HenrysLawCts

   IMPLICIT NONE

! !PUBLIC TYPES:

   PRIVATE
#include "mpif.h"

   include "netcdf.inc"

   PUBLIC  O3_GridComp       ! The O3 object 

! !PUBLIC MEMBER FUNCTIONS:

   PUBLIC  O3_GridCompSetServices
   PUBLIC  O3_GridCompInitialize
   PUBLIC  O3_GridCompRun
   PUBLIC  O3_GridCompFinalize

!
! !DESCRIPTION:

!  This module implements a parameterized chemistry for ozone that includes
!  dry deposition based on the GMIchem handling.

! !REVISION HISTORY:

!       2000 Nielsen   Initial coding
!   4Mar2005 Nielsen   Implementation of parameterized ozone chemistry
!  31Jan2011 Nielsen   Add dry deposition and NetCDF reads from PCHEM

!EOP
!-------------------------------------------------------------------------

  TYPE O3_GridComp

  LOGICAL :: DebugIsOn    ! For echoing the state of the states

! PCHEM climatology particulars
! -----------------------------
  CHARACTER(LEN=ESMF_MAXSTR) :: PCHEMfileName   ! NetCDF file borrowed from PCHEM
  
  INTEGER :: NSPECIES = 7   ! Number of species.  Usually in order
                            !  OX, N2O, CFC-11, CFC-12, CH4, HCFC-22, and H2O.
  INTEGER :: climYears      ! Number of years
  INTEGER :: nlatsPCHEM     ! Number of latitudes in climatology
  INTEGER :: nlevsPCHEM     ! Number of layers in climatology
  INTEGER :: begClimYear    ! First year in PCHEM climatology
  INTEGER :: endClimYear    ! Last year
  INTEGER :: BCnymd
  INTEGER :: PCnymd

  REAL, POINTER, DIMENSION(:)	    :: lats => null() ! Latitudes
  REAL, POINTER, DIMENSION(:)	    :: levs => null() ! Layers
  REAL, POINTER, DIMENSION(:,:,:,:) :: mnpl => null() ! Production rates and loss frequencies, O3 only
  REAL, POINTER, DIMENSION(:,:,:)   :: mncv => null() ! Concentration (mole fraction), O3 only

! Dry deposition borrowed from GMIchem
! ------------------------------------
  INTEGER, ALLOCATABLE ::  ireg(:,:)
  INTEGER, ALLOCATABLE :: iland(:,:,:)
  INTEGER, ALLOCATABLE ::  iuse(:,:,:)
  REAL, ALLOCATABLE ::  xlai(:,:,:)

  CHARACTER(LEN=ESMF_MAXSTR) :: GMIvegFileName
  CHARACTER(LEN=ESMF_MAXSTR) :: GMIlaiFileName

! ----------------
! Ozone parameters
! ----------------

  REAL :: hstar = 0.01
  REAL :: oxidize = 1.00

! -------------------
! Integer parameters.
! -------------------

  INTEGER :: NPOLY    = 20 ! Number of coefficients for polynomial fits
  INTEGER :: NVEGTYPE = 74 ! Maximum number of surface types (Olson)
  INTEGER :: NTYPE    = 15 ! maximum number of vegetation types in a cell
  INTEGER :: NWATER   =  6

! ----------------
! Real parameters.
! ----------------

  REAL ::   KGPG = 0.001   ! Kilograms per gram

  LOGICAL :: firstRun

  END TYPE O3_GridComp

CONTAINS

   subroutine O3_GridCompSetServices(  gc, chemReg, rc)
   type(ESMF_GridComp), intent(INOUT) :: GC
   type(Chem_Registry), intent(INOUT) :: chemReg
   integer,             intent(OUT  ) :: rc

   CHARACTER(LEN=255) :: rcbasen = 'O3_GridComp'
   CHARACTER(LEN=255) :: name

   integer                    :: status
   character(len=ESMF_MAXSTR) :: Iam
   type(O3_GridComp)          :: gcO3
   integer                    :: ic
   CHARACTER(LEN= 3 ) :: vegID
   CHARACTER(LEN=255) :: vegName


   Iam = "O3_GridCompSetServices"

   call MAPL_AddImportSpec(GC,           &
        SHORT_NAME = 'O3_ireg',          &
        LONG_NAME  = 'O3_emissions'  ,   &
        UNITS      = 'kg s-1 m-2',       &
        DIMS       = MAPL_DimsHorzOnly,  &
        VLOCATION  = MAPL_VLocationNone, &
        RESTART    = MAPL_RestartSkip,   &
        RC         = STATUS)
   VERIFY_(STATUS)
   do ic = 1, gcO3%NTYPE

      WRITE(vegID,'(I3.3)') ic

      vegName = 'O3_iuseVegID'//vegID
      call MAPL_AddImportSpec(GC,          &
          SHORT_NAME = trim(vegName),      &
          LONG_NAME  = 'O3_emissions'  ,   &
          UNITS      = 'kg s-1 m-2',       &
          DIMS       = MAPL_DimsHorzOnly,  &
          VLOCATION  = MAPL_VLocationNone, &
          RESTART    = MAPL_RestartSkip,   &
          RC         = STATUS)
      VERIFY_(STATUS)
      vegName = 'O3_ilandVegID'//vegID
      call MAPL_AddImportSpec(GC,          &
          SHORT_NAME = trim(vegName),      &
          LONG_NAME  = 'O3_emissions'  ,   &
          UNITS      = 'kg s-1 m-2',       &
          DIMS       = MAPL_DimsHorzOnly,  &
          VLOCATION  = MAPL_VLocationNone, &
          RESTART    = MAPL_RestartSkip,   &
          RC         = STATUS)
      VERIFY_(STATUS)
      vegName = 'O3_laiVegID'//vegID
      call MAPL_AddImportSpec(GC,          &
          SHORT_NAME = trim(vegName),      &
          LONG_NAME  = 'O3_emissions'  ,   &
          UNITS      = 'kg s-1 m-2',       &
          DIMS       = MAPL_DimsHorzOnly,  &
          VLOCATION  = MAPL_VLocationNone, &
          RESTART    = MAPL_RestartSkip,   &
          RC         = STATUS)
      VERIFY_(STATUS)
   enddo

   RETURN_(ESMF_SUCCESS)

   end subroutine O3_GridCompSetServices

!-------------------------------------------------------------------------
!NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS!
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  O3_GridCompInitialize --- Initialize O3_GridComp
!
! !INTERFACE:
!

   SUBROUTINE O3_GridCompInitialize ( gcO3, w_c, impChem, expChem, &
                                      nymd, nhms, cdt, rc )

! !USES:

  IMPLICIT none

! !INPUT PARAMETERS:

   TYPE(Chem_Bundle), INTENT(IN) :: w_c     ! Chemical tracer fields, delp, +
   INTEGER, INTENT(IN) :: nymd, nhms	    ! time
   REAL,    INTENT(IN) :: cdt		    ! chemistry time step (secs)

! !OUTPUT PARAMETERS:

   TYPE(O3_GridComp), INTENT(INOUT) :: gcO3	! Grid Component
   TYPE(ESMF_State), INTENT(INOUT)  :: impChem  ! Import State
   TYPE(ESMF_State), INTENT(INOUT)  :: expChem  ! Export State
   INTEGER, INTENT(OUT) ::  rc                  ! Error return code:

! !DESCRIPTION: Initializes the O3 Grid Component.

! !REVISION HISTORY:

!  18Sep2003 da Silva  First crack.
!   4Mar2005 Nielsen   Implementation of parameterized ozone chemistry
!  31Jan2011 Nielsen   Add dry deposition and NetCDF reads from PCHEM

!EOP
!-------------------------------------------------------------------------

   CHARACTER(LEN=ESMF_MAXSTR), PARAMETER :: Iam = 'O3_GridCompInitialize'
   CHARACTER(LEN=ESMF_MAXSTR) :: rcFileName = 'O3_GridComp.rc'

   TYPE(ESMF_VM) :: vm

   INTEGER :: n, status
   INTEGER :: i, i1, i2, ic, im, j, j1, j2, jm, km
   REAL :: c1,c2,c3,c4
   
   REAL, ALLOCATABLE :: veg2D(:,:)
   CHARACTER(LEN= 3) :: vegID
   CHARACTER(LEN=15) :: vegName

   gcO3%BCnymd = -1
   gcO3%PCnymd = -1
   gcO3%firstRun = .true.

! Grab the virtual machine
! ------------------------
   CALL ESMF_VMGetCurrent(vm, RC=status)
   VERIFY_(status)

! Initialize local variables
! --------------------------
   rc = 0
   status = 0
   i1 = w_c%grid%i1
   i2 = w_c%grid%i2
   im = w_c%grid%im

   j1 = w_c%grid%j1
   j2 = w_c%grid%j2
   jm = w_c%grid%jm
   
   km = w_c%grid%km

! Load resource file
! ------------------
   CALL I90_loadf(TRIM(rcFileName), status)
   VERIFY_(status)

! Parse resource file
! -------------------   
   CALL I90_label("DEBUG:", status)
   VERIFY_(status)
   i = I90_gint(status)
   VERIFY_(status)
   IF(i == 0) THEN
    gcO3%DebugIsOn = .FALSE.
   ELSE
    gcO3%DebugIsOn = .TRUE.
   END IF

   CALL I90_label("PCHEMs_file_name:", status)
   VERIFY_(status)
   CALL I90_Gtoken(gcO3%PCHEMfileName, status)
   VERIFY_(status)

   CALL I90_label("pchem_clim_years:", status)
   VERIFY_(status)
   gcO3%climYears = I90_gint(status)
   VERIFY_(status)

! PCHEM: Perform the initialization for
! establishing the production rates and loss frequencies
! ------------------------------------------------------
   CALL setUpPandL(RC=status)
   VERIFY_(status)
   _ASSERT(gcO3%nlevsPCHEM == km,'needs informative message')

! GMIchem: Obtain static vegetation properties for dry deposition
! ---------------------------------------------------------------
   CALL I90_label("veg_file_name:", status)
   VERIFY_(status)
   CALL I90_Gtoken(gcO3%GMIvegFileName, status)
   VERIFY_(status)
   CALL I90_label("lai_file_name:", status)
   VERIFY_(status)
   CALL I90_Gtoken(gcO3%GMIlaiFileName, status)
   VERIFY_(status)

   ALLOCATE(gcO3%ireg (i1:i2, j1:j2), STAT=status)
   VERIFY_(status)
   ALLOCATE(gcO3%iland(i1:i2, j1:j2, gcO3%NTYPE), STAT=status)
   VERIFY_(status)
   ALLOCATE(gcO3%iuse (i1:i2, j1:j2, gcO3%NTYPE), STAT=status)
   VERIFY_(status)
   ALLOCATE(gcO3%xlai (i1:i2, j1:j2, gcO3%NTYPE), STAT=status)
   VERIFY_(status)
   gcO3%ireg  = 0
   gcO3%iland = 0
   gcO3%iuse  = 0

!  Get Henrys Law cts for parameterized convective wet removal
!  -----------------------------------------------------------
   CALL get_HenrysLawCts('O3',c1,c2,c3,c4)  
   w_c%reg%Hcts(1,w_c%reg%i_O3 : w_c%reg%j_O3)=c1
   w_c%reg%Hcts(2,w_c%reg%i_O3 : w_c%reg%j_O3)=c2
   w_c%reg%Hcts(3,w_c%reg%i_O3 : w_c%reg%j_O3)=c3
   w_c%reg%Hcts(4,w_c%reg%i_O3 : w_c%reg%j_O3)=c4

   RETURN

CONTAINS

!-------------------------------------------------------------------------
!NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS!
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  setUpPandL
!
! !INTERFACE:
!
   SUBROUTINE setUpPandL(rc)

! !USES:

  IMPLICIT NONE

! !INPUT PARAMETERS

! !OUTPUT PARAMETERS:

   INTEGER, OPTIONAL, INTENT(OUT) ::  rc  ! Error return code:

! !DESCRIPTION: Read PCHEM's NetCDF file and distribute. Code borrowed from
!               GEOS\_PChemGridComp.F90 with minor modifications. For use 
!               with one-year datasets ONLY!

! !REVISION HISTORY:

!EOP
!-------------------------------------------------------------------------

    INTEGER :: dimid, varid, nspecies, comm, info, climYears, unit
    CHARACTER(LEN=ESMF_MAXSTR), PARAMETER :: Iam = 'O3_GridCompInitialize::setUpPandL'

    rc = 0 

! Get the communicator from the virtual machine and open the PCHEM file
! ---------------------------------------------------------------------
    CALL ESMF_VMGet(vm, MPICOMMUNICATOR=comm, RC=status)
    VERIFY_(status)

#undef H5_HAVE_PARALLEL
#ifdef H5_HAVE_PARALLEL

    CALL MPI_Info_create(info, status)
    VERIFY_(status)
    CALL MPI_Info_set(info, "romio_cb_read", "automatic", status)
    VERIFY_(status)

#ifdef NETCDF_NEED_NF_MPIIO
    status = NF_OPEN_PAR(TRIM(gcO3%PCHEMfileName), IOR(NF_NOWRITE,NF_MPIIO), comm, info, unit)
#else
    status = NF_OPEN_PAR(TRIM(gcO3%PCHEMfileName), NF_NOWRITE, comm, info, unit)
#endif

#else

    IF(MAPL_AM_I_ROOT(vm) ) THEN
       STATUS = NF_OPEN(TRIM(gcO3%PCHEMfileName), NF_NOWRITE, unit)

#endif

    IF(status /= NF_NOERR) THEN
       PRINT *,"Error opening file ",TRIM(gcO3%PCHEMfileName), status
       PRINT *, NF_STRERROR(status)
       status = 1
       VERIFY_(status)
    END IF

! Obtain dimensions of the latitudes, layers, and species in the PCHEM file
! -------------------------------------------------------------------------
    status = NF_INQ_DIMID(unit, 'lat', dimid)
    IF(status /= NF_NOERR) THEN
       PRINT *,"Error getting dimid for lat", status
       PRINT *, NF_STRERROR(status)
       status = 1
       VERIFY_(status)
    END IF
    status = NF_INQ_DIMLEN(unit, dimid, gcO3%nlatsPCHEM)
    IF(status /= nf_noerr) then
       PRINT *,"Error getting dimlen for lat", status
       PRINT *, NF_STRERROR(status)
       status = 1
       VERIFY_(status)
    END IF
    status = NF_INQ_DIMID(unit, 'lev', dimid)
    IF(status /= nf_noerr) THEN
       PRINT *,"Error getting dimid for lev", status
       PRINT *, NF_STRERROR(status)
       status = 1
       VERIFY_(status)
    END IF
    status = NF_INQ_DIMLEN(unit, dimid, gcO3%nlevsPCHEM)
    IF(status /= nf_noerr) THEN
       PRINT *,"Error getting dimlen for lev", status
       PRINT *, NF_STRERROR(status)
       status = 1
       VERIFY_(status)
    END IF
    status = NF_GET_ATT_INT(unit, NF_GLOBAL, 'NSPECIES', nspecies)
    IF(status /= nf_noerr) THEN
       PRINT *,"Error getting nspecies", status
       PRINT *, NF_STRERROR(status)
       status = 1
       VERIFY_(status)
    END IF
    _ASSERT(gcO3%NSPECIES == nspecies,'needs informative message')

! Validate the length of the climatology
! --------------------------------------
    gcO3%begClimYear = 1
    gcO3%endClimYear = 1

#ifndef H5_HAVE_PARALLEL
    END IF ! MAPL_am_I_root

    CALL MAPL_CommsBcast(vm, gcO3%nlatsPCHEM, 1, 0, RC=status)
    VERIFY_(status)
    CALL MAPL_CommsBcast(vm, gcO3%nlevsPCHEM, 1, 0, RC=status)
    VERIFY_(status)
    CALL MAPL_CommsBcast(vm, gcO3%begClimYear, 1, 0, RC=status)
    VERIFY_(status)
    CALL MAPL_CommsBcast(vm, gcO3%endClimYear, 1, 0, RC=status)
    VERIFY_(status)

#endif

! Allocate and broadcast the latitudes and layers
! -----------------------------------------------
    ALLOCATE(gcO3%lats(gcO3%nlatsPCHEM), STAT=status)
    VERIFY_(status)
    ALLOCATE(gcO3%levs(gcO3%nlevsPCHEM), STAT=status)
    VERIFY_(status)

#ifndef H5_HAVE_PARALLEL
    IF ( MAPL_AM_I_ROOT(vm) ) THEN
#endif

    status = NF_INQ_VARID(unit, 'lat', varid)
    IF(status /= NF_NOERR) THEN
       PRINT *,"Error getting varid for lat", status
       PRINT *, NF_STRERROR(status)
       status = 1
       VERIFY_(status)
    END IF
    status = NF_GET_VAR_REAL(unit, varid, gcO3%lats)
    IF(status /= NF_NOERR) THEN
       PRINT *,'Error getting values for lat', status
       PRINT *, NF_STRERROR(status)
       status = 1
       VERIFY_(status)
    END IF
    status = NF_INQ_VARID(unit, 'lev', varid)
    IF(status /= NF_NOERR) THEN
       PRINT *,"Error getting varid for lev", status
       PRINT *, NF_STRERROR(status)
       status = 1
       VERIFY_(status)
    END IF
    status = NF_GET_VAR_REAL(unit, varid, gcO3%levs)
    IF(status /= NF_NOERR) THEN
       PRINT *,"Error getting values for lev", status
       PRINT *, NF_STRERROR(status)
       status = 1
       VERIFY_(status)
    END IF

#ifdef H5_HAVE_PARALLEL

    CALL MPI_Info_free(info, status)
    VERIFY_(status)

#else

    status = NF_CLOSE(unit)
    VERIFY_(status)

    END IF ! MAPL_am_I_root

    CALL MAPL_CommsBcast(vm, gcO3%lats, SIZE(gcO3%lats), 0, RC=status)
    VERIFY_(status)
    CALL MAPL_CommsBcast(vm, gcO3%levs, SIZE(gcO3%levs), 0, RC=status)
    VERIFY_(status)

#endif

! Allocate space for concentration and production rates and
! loss frequencies. Note that we will be working with ozone only.
!----------------------------------------------------------------
    ALLOCATE(gcO3%mncv(gcO3%nlatsPCHEM, gcO3%nlevsPCHEM, 2), STAT=status)
    VERIFY_(status)
    gcO3%mncv = Z'7FA00000'

    ALLOCATE(gcO3%mnpl(gcO3%nlatsPCHEM, gcO3%nlevsPCHEM, 2, 2), STAT=status)
    VERIFY_(status)
    gcO3%mnpl = Z'7FA00000'

    rc = 0 

   RETURN
   END SUBROUTINE setUpPandL

  END SUBROUTINE O3_GridCompInitialize

!-------------------------------------------------------------------------
!NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS!
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  O3_GridCompRun --- The O3 run method 
!
! !INTERFACE:
!

   SUBROUTINE O3_GridCompRun ( gcO3, w_c, impChem, expChem, &
                               nymd, nhms, cdt, rc )

! !USES:

  IMPLICIT none

! !INPUT/OUTPUT PARAMETERS:

   TYPE(O3_GridComp), INTENT(INOUT) :: gcO3   ! Grid Component
   TYPE(Chem_Bundle), INTENT(INOUT) :: w_c    ! Chemical tracer fields   

! !INPUT PARAMETERS:

   TYPE(ESMF_State), INTENT(INOUT) :: impChem ! Import State
   INTEGER, INTENT(IN) :: nymd, nhms	      ! time
   REAL,    INTENT(IN) :: cdt		      ! chemical timestep (secs)


! !OUTPUT PARAMETERS:

   TYPE(ESMF_State), INTENT(INOUT) :: expChem   ! Export State
   INTEGER, INTENT(OUT) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 -
 
! !DESCRIPTION: This routine implements a parameterized chemistry for
!               ozone. 
!
! !REVISION HISTORY:
!
!  18Sep2003 da Silva  First crack.
!   4Mar2005 Nielsen   Implementation of parameterized ozone chemistry
!  31Jan2012 Nielsen   Revisions for running dry deposition plagarized
!                      from GMIChem
!
!EOP
!-------------------------------------------------------------------------

   CHARACTER(LEN=ESMF_MAXSTR), PARAMETER :: Iam = 'O3_GridCompRun'
   TYPE(ESMF_VM) :: vm

! Local
! -----
   INTEGER :: BCymd, BChms
   INTEGER :: i1, i2, ic, im, i, j, j1, j2, jm, km, ixj
   INTEGER :: k, n
   INTEGER :: nCells
   INTEGER :: nFirstO3, numberO3s
   INTEGER :: status

   REAL :: qmin, qmax
   REAL, PARAMETER :: von_karman = 0.40      ! von Karman constant
   REAL, PARAMETER :: DOBSONS_PER_MOLE = MAPL_AVOGAD/2.69E+20

!  Input fields from fvGCM
!  -----------------------
   REAL, POINTER, DIMENSION(:,:)   ::  cldtt   => null()
   REAL, POINTER, DIMENSION(:,:)   ::  lwi     => null()
   REAL, POINTER, DIMENSION(:,:)   ::  pblh    => null()
   REAL, POINTER, DIMENSION(:,:)   ::  swndsrf => null()
   REAL, POINTER, DIMENSION(:,:)   ::  shFlux  => null()
   REAL, POINTER, DIMENSION(:,:)   ::  srfAirT => null()
   REAL, POINTER, DIMENSION(:,:)   ::  tropp   => null()
   REAL, POINTER, DIMENSION(:,:)   ::  ustar   => null()
   REAL, POINTER, DIMENSION(:,:)   ::  z0h     => null()
   REAL, POINTER, DIMENSION(:,:,:) ::  ple     => null()
   REAL, POINTER, DIMENSION(:,:,:) ::  T       => null()
   REAL, POINTER, DIMENSION(:,:,:) ::  rhoa    => null()
   REAL, POINTER, DIMENSION(:,:,:) ::  zle     => null()

! Quantities to be exported
! -------------------------
   REAL, POINTER, DIMENSION(:,:)   :: o3tot  => null()
   REAL, POINTER, DIMENSION(:,:)   :: o3ddp  => null()
   REAL, POINTER, DIMENSION(:,:)   :: o3ddv  => null()
   REAL, POINTER, DIMENSION(:,:,:) :: o3     => null()
   REAL, POINTER, DIMENSION(:,:,:) :: ox     => null()
   REAL, POINTER, DIMENSION(:,:,:) :: o3ppmv => null()
   REAL, POINTER, DIMENSION(:,:,:) :: o3tend => null()
   REAL, POINTER, DIMENSION(:,:,:) :: oxtend => null()

   REAL, ALLOCATABLE :: initialO3(:,:,:)
   REAL, ALLOCATABLE :: plPa(:,:,:)
   REAL, ALLOCATABLE :: cellDepth(:,:)
   REAL, ALLOCATABLE :: obk(:,:)
   REAL, ALLOCATABLE :: dvel(:,:)
   REAL, ALLOCATABLE :: dryDepFreq(:,:)
   REAL, ALLOCATABLE :: dO3(:,:)
   INTEGER, ALLOCATABLE :: oro(:,:)

   REAL, ALLOCATABLE :: lai2D(:,:)
   CHARACTER(LEN= 3) :: laiID
   CHARACTER(LEN= 3) :: vegID
   CHARACTER(LEN=15) :: laiName
   CHARACTER(LEN=255) :: vegName
   real, pointer     :: ptr2d(:,:) => null()

!     --------------------------------------------------------------------
!     IDEP   : deposition surface type for each Olson surface type
!     IRAC   : resistance that depends on canopy height and density (s^-1)
!     IRCLO  : resistance for leaves, twig, bark in lower canopy    (s^-1)
!     IRCLS  : resistance for leaves, twig, bark in lower canopy    (s^-1)
!     IRGSO  : ground    resistance (s^-1)
!     IRGSS  : ground    resistance (s^-1)
!     IRI    : internal  resistance (s^-1)
!     IRLU   : cuticular resistance (s^-1)
!     IWATER : id's for surface types that are water
!     IZO    : roughness height (m/10000)
!     NWATER : number of Olson's surface types that are water
!     --------------------------------------------------------------------

  INTEGER, PARAMETER :: NPOLY    = 20 ! Number of coefficients for polynomial fits
  INTEGER, PARAMETER :: NVEGTYPE = 74 ! Maximum number of surface types (Olson)
  INTEGER, PARAMETER :: NTYPE    = 15 ! maximum number of vegetation types in a cell

      INTEGER, PARAMETER :: IDEP(NVEGTYPE) = (/  &
      11, 10,  5,  1,  1,  1,  2,  1,  8,  1,  &
       1,  1,  1,  1,  1,  1,  5,  1,  1,  1,  &
       3,  3,  3,  3,  2,  2,  2,  3,  2,  2,  &
       4,  4,  2,  6,  1,  1,  9,  4,  4,  4,  &
       5,  5,  5,  5,  5,  9,  5,  5,  5,  5,  &
       8,  8,  5,  7,  6,  2,  2,  2,  2,  2,  &
       3,  3,  3,  5,  5, 11, 11, 11, 11,  8,  &
       1,  8,  9, 11 /)

      INTEGER, PARAMETER :: IRAC(NTYPE) = (/  &
     	  0, 2000, 2000,  200,  100,  &
       2000,	0,    0,  300,  100,  &
     	  0,	0,    0,    0,    0 /)

      INTEGER, PARAMETER :: IRCLO(NTYPE) = (/  &
       1000, 1000, 1000, 1000, 1000,  &
       9999, 9999, 9999, 1000, 9999,  &
       9999,	0,    0,    0,    0 /)

      INTEGER, PARAMETER :: IRCLS(NTYPE) = (/  &
       9999, 2000, 2000, 2000, 2000,  &
       9999, 9999, 9999, 2500, 9999,  &
       9999,	0,    0,    0,    0 /)

      INTEGER, PARAMETER :: IRGSO(NTYPE) = (/  &
       3500,  200,  200,  150,  200,  &
     	200,  340,  400, 1000,  300,  &
       2000,	0,    0,    0,    0 /)

      INTEGER, PARAMETER :: IRGSS(NTYPE) = (/  &
     	100,  500,  500,  150,  350,  &
     	200,  340, 1000,    0,  400,  &
     	  0,	0,    0,    0,    0 /)

      INTEGER, PARAMETER :: IRI(NTYPE) = (/  &
       9999,  200,  400,  200,  200,  &
     	200,  200, 9999,  200, 9999,  &
       9999,	0,    0,    0,    0 /)

      INTEGER, PARAMETER :: IRLU(NTYPE) = (/  &
       9999, 9000, 9000, 9000, 9000,  &
       1000, 4000, 9999, 9000, 9999,  &
       9999,	0,    0,    0,    0 /)

      INTEGER, PARAMETER :: IWATER(NVEGTYPE) = (/  &
       1, 66, 67, 68, 69, 74,  0,  0,  0,  0,  &
       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
       0,  0,  0,  0 /)

      INTEGER, PARAMETER :: IZO(NVEGTYPE) = (/  &
     	   1,10000,   50, 1000, 1000, 1000,10000, 1000,    1, 1000,  &
     	1000, 1000, 1000, 1000, 1000, 1000, 1000,    1, 1000, 1000,  &
       10000,10000,10000,10000,10000,10000,10000,10000, 1000,10000,  &
     	1000, 1000, 1000,10000, 1000, 1000,  100, 1000, 1000, 1000,  &
     	 100,  100,  100,  100,  100,  100, 1000, 1000, 1000, 1000,  &
     	 100,  100,  100,   50,10000, 1000, 1000, 1000, 1000, 1000,  &
       10000,10000,10000, 1000,   50,	 1,    1,    1,    1,	10,  &
     	  10,	 1,  500,    1 /)

!     ------------------------------------------
!     DRYCOEFF : polynomial fitting coefficients
!     ------------------------------------------

      REAL, PARAMETER :: DRYCOEFF(NPOLY) = (/  &
       -3.58E-01,  3.02E+00,  3.85E+00, -9.78E-02, -3.66E+00,  &
     	1.20E+01,  2.52E-01, -7.80E+00,  2.26E-01,  2.74E-01,  &
     	1.14E+00, -2.19E+00,  2.61E-01, -4.62E+00,  6.85E-01,  &
       -2.54E-01,  4.37E+00, -2.66E-01, -1.59E-01, -2.06E-01 /)

! Grab the virtual machine
! ------------------------
   CALL ESMF_VMGetCurrent(vm, RC=status)
   VERIFY_(status)

!  Grid specs from Chem_Bundle%grid
!  --------------------------------
   rc = 0
   i1 = w_c%grid%i1
   i2 = w_c%grid%i2
   im = w_c%grid%im
   
   j1 = w_c%grid%j1
   j2 = w_c%grid%j2
   jm = w_c%grid%jm
   
   km = w_c%grid%km

   ixj = (i2-i1+1)*(j2-j1+1)
   
!  Get pointers to imports
!  -----------------------
   CALL MAPL_GetPointer(impChem, cldtt,   'CLDTT',   RC=status)
   VERIFY_(status)
   CALL MAPL_GetPointer(impChem, lwi,     'LWI',     RC=status)
   VERIFY_(status)
   CALL MAPL_GetPointer(impChem, pblh,    'ZPBL',    RC=status)
   VERIFY_(status)
   CALL MAPL_GetPointer(impChem, shFlux,  'SH',      RC=status)
   VERIFY_(status)
   CALL MAPL_GetPointer(impChem, swndsrf, 'SWNDSRF', RC=status)
   VERIFY_(status)
   CALL MAPL_GetPointer(impChem, srfAirT, 'TA',      RC=status)
   VERIFY_(status)
   CALL MAPL_GetPointer(impChem, tropp,   'TROPP',   RC=status)
   VERIFY_(status)
   CALL MAPL_GetPointer(impChem, ustar,   'USTAR',   RC=status)
   VERIFY_(status)
   CALL MAPL_GetPointer(impChem, z0h,     'Z0H',     RC=status)
   VERIFY_(status)
   CALL MAPL_GetPointer(impChem, ple,     'PLE',     RC=status)
   VERIFY_(status)
   CALL MAPL_GetPointer(impChem, T,       'T',       RC=status)
   VERIFY_(status)
   CALL MAPL_GetPointer(impChem, rhoa,    'AIRDENS', RC=status)
   VERIFY_(status)
   CALL MAPL_GetPointer(impChem, zle,     'ZLE',     RC=status)
   VERIFY_(status)

   IF(gcO3%DebugIsOn) THEN
    CALL pmaxmin('O3: cldtt',   cldtt,   qmin, qmax, ixj, 1,    1.)
    CALL pmaxmin('O3: lwi',     lwi,     qmin, qmax, ixj, 1,    1.)
    CALL pmaxmin('O3: pblh',    pblh,    qmin, qmax, ixj, 1,    1.)
    CALL pmaxmin('O3: shFlux',  shFlux,  qmin, qmax, ixj, 1,    1.)
    CALL pmaxmin('O3: swndsrf', swndsrf, qmin, qmax, ixj, 1,    1.)
    CALL pmaxmin('O3: ustar',   ustar,   qmin, qmax, ixj, 1,    1.)
    CALL pmaxmin('O3: tropp',   tropp,   qmin, qmax, ixj, 1,    1.)
    CALL pmaxmin('O3: z0h',     z0h,     qmin, qmax, ixj, 1,    1.)
    CALL pmaxmin('O3: ple',     ple,     qmin, qmax, ixj, km+1, 1.)
    CALL pmaxmin('O3: T',       T,       qmin, qmax, ixj, km,   1.)
    CALL pmaxmin('O3: rhoa',    rhoa,    qmin, qmax, ixj, km,   1.)
    CALL pmaxmin('O3: zle',     zle,     qmin, qmax, ixj, km+1, 1.)
    CALL pmaxmin('O3: srfAirT', srfAirT, qmin, qmax, ixj,    1, 1.)
   END IF

!  Get pointers to exports
!  -----------------------
   CALL MAPL_GetPointer(expChem, o3,     'O3',      RC=status)
   VERIFY_(status)
   CALL MAPL_GetPointer(expChem, ox,     'OX',      RC=status)
   VERIFY_(status)
   CALL MAPL_GetPointer(expChem, o3tot,  'O3TOT',   RC=status)
   VERIFY_(status)
   CALL MAPL_GetPointer(expChem, o3ddp,  'O3DDP',   RC=status)
   VERIFY_(status)
   CALL MAPL_GetPointer(expChem, o3ddv,  'O3DDV',   RC=status)
   VERIFY_(status)
   CALL MAPL_GetPointer(expChem, o3ppmv, 'O3PPMV',  RC=status)
   VERIFY_(status)
   CALL MAPL_GetPointer(expChem, o3tend, 'DO3DT' ,  RC=status)
   VERIFY_(status)
   CALL MAPL_GetPointer(expChem, oxtend, 'OX_TEND', RC=status)
   VERIFY_(status)

   if (gcO3%firstRun) then
      gcO3%firstRun = .false.

      vegName = 'O3_ireg'
      call MAPL_GetPointer(impChem,ptr2d,vegName,rc=status)
      VERIFY_(STATUS)

      gcO3%ireg(:,:) = INT(ptr2D(:,:))

      DO ic = 1, gcO3%NTYPE

       WRITE(vegID,'(I3.3)') ic

       vegName = 'O3_iuseVegID'//vegID
       call MAPL_GetPointer(impChem,ptr2d,vegName,rc=status)
       VERIFY_(STATUS)
       gcO3%iuse(:,:,ic) = INT(ptr2D(:,:))

       vegName = 'O3_ilandVegID'//vegID
       call MAPL_GetPointer(impChem,ptr2d,vegName,rc=status)
       VERIFY_(STATUS)
       gcO3%iland(:,:,ic) = INT(ptr2D(:,:))

      END DO

!     Until the ireg = -1 issue is resolved, do the following bug fix
!     ---------------------------------------------------------------
      DO j = j1,j2
       DO i = i1,i2
        IF(gcO3%ireg(i,j) == -1) THEN
         gcO3%iuse(i,j,1:gcO3%NTYPE) = 0
         gcO3%iland(i,j,1) = 1000
         gcO3%iland(i,j,2:gcO3%NTYPE) = 0
        END IF
       END DO
      END DO
      WHERE(gcO3%ireg(:,:) == -1) gcO3%ireg(:,:) = 1

   end if

   DO ic = 1, gcO3%NTYPE

      WRITE(laiID,'(I3.3)') ic
      laiName = 'O3_laiVegID'//laiID
      call MAPL_GetPointer(impChem,ptr2d,laiName,rc=status)
      VERIFY_(STATUS)
      gcO3%xlai(:,:,ic) = INT(ptr2D(:,:))

   END DO

! Save current O3
! ---------------
   ALLOCATE(initialO3(i1:i2,j1:j2,km),STAT=status)
   VERIFY_(status)
   n = w_c%reg%i_O3
   initialO3(:,:,:) = w_c%qa(n)%data3d(:,:,:)

! Middle-layer pressures
! ----------------------
   ALLOCATE(plPa(i1:i2,j1:j2,km),STAT=status)
   VERIFY_(status)
   plPa = 0.50*(ple(:,:,0:km-1)+ple(:,:,1:km))
   IF(gcO3%DebugIsOn) THEN
    CALL pmaxmin('O3: plPa', plPa, qmin, qmax, ixj, km, 1.)
   END IF

! Repair bad tropopause pressures, if any exist
! ---------------------------------------------
   CALL Chem_UtilTroppFixer(i2, j2, tropp, VERBOSE=.TRUE., RC=status)
   VERIFY_(status)

! Perform parameterized production and loss chemistry
! ---------------------------------------------------
   CALL doProdLoss(status)
   VERIFY_(status)

! Grab some memory
! ----------------
   ALLOCATE(cellDepth(i1:i2,j1:j2), STAT=status)
   VERIFY_(status)
   ALLOCATE(obk(i1:i2,j1:j2), STAT=status)
   VERIFY_(status)
   ALLOCATE(oro(i1:i2,j1:j2), STAT=status)
   VERIFY_(status)
   ALLOCATE(dvel(i1:i2,j1:j2), STAT=status)
   VERIFY_(status)
   ALLOCATE(dryDepFreq(i1:i2,j1:j2), STAT=status)
   VERIFY_(status)
   ALLOCATE(dO3(i1:i2,j1:j2), STAT=status)
   VERIFY_(status)

! Thickness of the surface layer
! ------------------------------
   cellDepth(:,:) = zle(:,:,km-1) - zle(:,:,km)
   IF(gcO3%DebugIsOn) THEN
    CALL pmaxmin('O3: cellDepth', cellDepth, qmin, qmax, ixj, 1, 1.)
   END IF

! Calculate the Obukhov length scale, obk. Code is from Colarco.
! If obk < 0 the air is unstable, and if OBK > 0 the air is stable.
! For sensible heat flux == 0, obk is set to a large value, 1.00E+05.
! -------------------------------------------------------------------
   WHERE(ABS(shflux) > 1.00E-32)
    obk(:,:) = - rhoa(:,:,km) * MAPL_CP * T(:,:,km) * ustar**3. / (von_karman * MAPL_GRAV * shFlux)
   ELSEWHERE
    obk(:,:) = 1.00E+05
   END WHERE

! Add one to the land-water-ice mask
! ----------------------------------
   oro(:,:) = lwi(:,:)+1
   dvel(:,:) = 0.00
   
! What is the deposition velocity [m s^{-1}]?
! -------------------------------------------
   DO j = j1,j2
    nCells = i2-i1+1
    CALL DeposVelo(nCells, j, swndsrf(:,j), srfAirT(:,j), w_c%cosz(:,j), gcO3%oxidize, &
		   gcO3%hstar, MAPL_O3MW, ustar(:,j), cellDepth(:,j)*0.50, obk(:,j), &
		   cldtt(:,j), oro(:,j), rhoa(:,j,km), dvel(:,j), rc)
    IF(rc /= 0) THEN
     PRINT *,TRIM(Iam)//": ERROR in GOCART::O3 DeposVelo"
     status = rc
     VERIFY_(status)
    END IF

   END DO

   IF(gcO3%DebugIsOn) THEN
    CALL pmaxmin('O3: dvel', dvel, qmin, qmax, ixj, 1, 1.)
   END IF

! Set a minimum value
! -------------------
   WHERE(dvel(:,:) < 1.00E-04) dvel(:,:) = 1.00E-04

! Fill export state for dry deposition speed [m s^{-1}]
! -----------------------------------------------------
   IF(ASSOCIATED(o3ddv)) o3ddv(:,:) = dvel(:,:)

! Dry deposition frequency [s^{-1}] for the chemical removal term
! ---------------------------------------------------------------
   dryDepFreq(:,:) = dvel(:,:)/cellDepth(:,:)

! Concentration increment (mol/mol)
! ---------------------------------
   n = w_c%reg%i_O3
   dO3(:,:) = w_c%qa(n)%data3d(:,:,km)*(1.00-EXP(-dryDepFreq(:,:)*cdt))
   WHERE(dO3(:,:) < 0.00) dO3(:,:) = 0.00

! Update concentration internal state: GOCART::OX, volume mixing ratio
! --------------------------------------------------------------------
   w_c%qa(n)%data3d(:,:,km) =  w_c%qa(n)%data3d(:,:,km) - dO3(:,:)
   IF(gcO3%DebugIsOn) THEN
    CALL pmaxmin('O3: OX(ppmv)', w_c%qa(n)%data3d, qmin, qmax, ixj, km, 1.00E+06)
   END IF

! Fill export state for dry deposition [kg m^{-2} s^{-1}]
! -------------------------------------------------------
   IF(ASSOCIATED(o3ddp)) o3ddp(:,:) = dO3(:,:)*cellDepth(:,:)*rhoa(:,:,km)*MAPL_O3MW/(cdt*MAPL_AIRMW)

!  ------------------------------------------------------
!  Fill export states
!
!   Name    Units              Contents
!  ------- ------------------- --------------------------
!  O3      kg kg^{-1}          Ozone mass fraction
!  OX      mol mol^{-1}        Ozone volume mixing ratio*
!  O3PPMV  ppmv                OX*1.00E+06
!  OX_TEND mol mol^{-1} s^{-1} Tendency
!  O3_TEND kg kg^{-1} s^{-1}   Tendency
!
!  *OX is necessary so that CHEM sees an ANALYSIS_OX 
!   with the same name, irregardless of the PROVIDER.
!  ------------------------------------------------------
   n = w_c%reg%i_O3

   IF(ASSOCIATED(ox)) ox(i1:i2,j1:j2,1:km) = w_c%qa(n)%data3d(i1:i2,j1:j2,1:km)

   DO k=1,km

    WHERE(plPa(:,:,k) < 100.00 .AND. w_c%cosz > 0.00)
       dvel(:,:) = w_c%qa(n)%data3d(i1:i2,j1:j2,k)*EXP(-1.50*(LOG10(plPa(:,:,k))-2.00)**2)
    ELSEWHERE
       dvel(:,:) = w_c%qa(n)%data3d(i1:i2,j1:j2,k)
    END WHERE

    IF(ASSOCIATED(o3)) o3(i1:i2,j1:j2,k) = dvel(i1:i2,j1:j2)*MAPL_O3MW/MAPL_AIRMW
    IF(ASSOCIATED(o3ppmv)) o3ppmv(i1:i2,j1:j2,k) = dvel(i1:i2,j1:j2)*1.00E+06

   END DO

! Total ozone (Dobsons)
! ---------------------
   IF(ASSOCIATED(o3tot)) THEN
    o3tot(:,:) = 0.00
    DO k=1,km     
     o3tot(:,:) = o3tot(:,:)+w_c%qa(n)%data3d(i1:i2,j1:j2,k)*w_c%delp(:,:,k)*(DOBSONS_PER_MOLE/(MAPL_AIRMW*MAPL_GRAV))  
    END DO
   END IF

! Ozone tendency
! --------------
   n = w_c%reg%i_O3
   IF(ASSOCIATED(oxtend)) oxtend(:,:,:) = (w_c%qa(n)%data3d(:,:,:)-initialO3(:,:,:))/cdt
   IF(ASSOCIATED(o3tend)) o3tend(:,:,:) = (w_c%qa(n)%data3d(:,:,:)-initialO3(:,:,:))*MAPL_O3MW/(MAPL_AIRMW*cdt)

! Clean up
! --------
   DEALLOCATE(plPa, STAT=status)
   VERIFY_(status)
   DEALLOCATE(obk, STAT=status)
   VERIFY_(status)
   DEALLOCATE(dvel, STAT=status)
   VERIFY_(status)
   DEALLOCATE(oro, STAT=status)
   VERIFY_(status)
   DEALLOCATE(dryDepFreq, STAT=status)
   VERIFY_(status)
   DEALLOCATE(dO3, STAT=status)
   VERIFY_(status)
   DEALLOCATE(cellDepth, STAT=status)
   VERIFY_(status)
   DEALLOCATE(initialO3, STAT=status)
   VERIFY_(status)

   rc = 0

   RETURN

CONTAINS

!-------------------------------------------------------------------------
!NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS!
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: doProdLoss
!
!  Run the parameterized chemistry for Ox. Ozone is derived from Ox.
!
! !INTERFACE:
!
   SUBROUTINE doProdLoss(rc)

! !USES:

  IMPLICIT NONE

   INTEGER, OPTIONAL, INTENT(OUT) :: rc

! !DESCRIPTION:
!
!  This module implements a parameterized chemistry for ozone.  The NetCDF
!  file that contains the production rates and loss frequencies has coefficients 
!  for seven species, OX, N2O, CFC-11, CFC-12, CH4, HCFC-22, and H2O.\\
!  
!  Advection produces the "intermediate" constituent distribution
!  before this routine is called.\\
!  
!  USAGE NOTES:\\
!  
!  The resulting O3 mole fraction is the product of the Ox mole fraction
!  multiplied by the O3-to-Ox ratio, ro3ox. At pressures greater than
!  approximately 1 hPa, ro3ox = 1 everywhere.  At pressures less than 
!  approximately 0.1 hPa, Ox is mostly O3 at night and ro3ox = 1.  During 
!  the day, ro3ox in this region depends to first order on pressure.\\
!
!  Code is plagarized from GEOS\_PchemGridComp.F90
!
! !REVISION HISTORY:
!
!  31Jan2011 Nielsen
!
!EOP
!-------------------------------------------------------------------------

   CHARACTER(LEN=ESMF_MAXSTR), PARAMETER :: Iam = 'O3_GridCompRun::doProdLoss'

   INTEGER :: i, j, k
   INTEGER :: yy, mm, dd, indx1, indx2, daysThisMonth
   INTEGER :: dimid, varid, comm, info, start(3), cnt(3), unit
   INTEGER :: status

   INTEGER, ALLOCATABLE :: mask(:,:,:)

   REAL :: fac

   REAL, ALLOCATABLE :: Pi(:,:,:)
   REAL, ALLOCATABLE :: Li(:,:,:)
   REAL, ALLOCATABLE :: Pclim(:,:)
   REAL, ALLOCATABLE :: Lclim(:,:)
   REAL, ALLOCATABLE :: P1(:,:)
   REAL, ALLOCATABLE :: L1(:,:)

   INTEGER, PARAMETER :: monthLength(12) = (/ 31,28,31,30,31,30,31,31,30,31,30,31 /)

   rc = 0
   yy = nymd/10000
   mm = (nymd-yy*10000)/100
   dd = nymd-yy*10000-mm*100

   IF(dd < 16) THEN
    indx1 = mm-1
    indx2 = mm
   ELSE
    indx1 = mm
    indx2 = mm+1
   END IF

   IF(indx1 < 1) indx1 = 12
   IF(indx2 > 12) indx2 = 1

   IF(mm == 2) THEN
    daysThisMonth = 28
    IF(MOD(yy,  4) == 0) daysThisMonth = 29
    IF(MOD(yy,100) == 0) daysThisMonth = 28
    IF(MOD(yy,400) == 0) daysThisMonth = 29
   ELSE
    daysThisMonth = monthLength(mm)
   END IF

   IF(dd < 16) THEN
    fac = 0.50*(1.00+(dd-1.00)/15.00)
   ELSE
    fac = 0.50*(1.00+(daysThisMonth+1.00-dd)/(daysThisMonth-15.00))
   END IF

   ChangeOfDay: IF(gcO3%PCnymd /= nymd) THEN
    gcO3%PCnymd = nymd

    CALL ESMF_VMGet(vm, MPICOMMUNICATOR=comm, RC=status)
    VERIFY_(status)

#undef H5_HAVE_PARALLEL
#ifdef H5_HAVE_PARALLEL

    CALL MPI_Info_create(info, status)
    VERIFY_(status)
    CALL MPI_Info_set(info, "romio_cb_read", "automatic", status)
    VERIFY_(status)

#ifdef NETCDF_NEED_NF_MPIIO
    status = NF_OPEN_PAR(TRIM(gcO3%PCHEMfileName), IOR(NF_NOWRITE,NF_MPIIO), comm, info, unit)
#else
    status = NF_OPEN_PAR(TRIM(gcO3%PCHEMfileName), NF_NOWRITE, comm, info, unit)
#endif

#else
    IF(MAPL_AM_I_ROOT(vm) ) THEN
       status = NF_OPEN(TRIM(gcO3%PCHEMfileName), NF_NOWRITE, unit)
#endif
    IF(status /= NF_NOERR) THEN
       PRINT *,"Error opening file ",TRIM(gcO3%PCHEMfileName), status
       PRINT *, NF_STRERROR(status)
       status = 1
       VERIFY_(status)
    END IF

    start(1) = 1
    start(2) = 1
    cnt(1) = gcO3%nlatsPCHEM
    cnt(2) = gcO3%nlevsPCHEM
    cnt(3) = 1

    status = NF_INQ_VARID(unit, "OX_PROD", varid)
    IF(status /= NF_NOERR) THEN
       PRINT *,"Error getting varid for variable OX_PROD", status
       PRINT *, NF_STRERROR(status)
       status = 1
       VERIFY_(status)
    END IF
    start(3) = indx1
    status = NF_GET_VARA_REAL(unit, varid, start, cnt, gcO3%mnpl(:,:,1,1))
    IF(status /= NF_NOERR) THEN
       PRINT *,"Error reading lower bracket month for production ",status
       PRINT *, NF_STRERROR(status)
       status = 1
       VERIFY_(status)
    END IF
    start(3) = indx2
    status = NF_GET_VARA_REAL(unit, varid, start, cnt, gcO3%mnpl(:,:,1,2))
    IF(status /= NF_NOERR) THEN
       PRINT *,"Error reading upper bracket month for production ",status
       PRINT *, NF_STRERROR(status)
       status = 1
       VERIFY_(status)
    END IF

    status = NF_INQ_VARID(unit, "OX_LOSS", varid)
    IF(status /= NF_NOERR) THEN
       PRINT *,"Error getting varid for variable OX_LOSS", status
       PRINT *, NF_STRERROR(status)
       status = 1
       VERIFY_(status)
    END IF
    start(3) = indx1
    status = NF_GET_VARA_REAL(unit, varid, start, cnt, gcO3%mnpl(:,:,2,1))
    IF(status /= NF_NOERR) THEN
       PRINT *,"Error reading lower bracket month for loss ",status
       PRINT *, NF_STRERROR(status)
       status = 1
       VERIFY_(status)
    END IF
    start(3) = indx2
    status = NF_GET_VARA_REAL(unit, varid, start, cnt, gcO3%mnpl(:,:,2,2))
    IF(status /= NF_NOERR) THEN
       PRINT *,"Error reading upper bracket month for loss ",status
       PRINT *, NF_STRERROR(status)
       status = 1
       VERIFY_(status)
    END IF

#ifdef H5_HAVE_PARALLEL
    CALL MPI_Info_free(info, status)
    VERIFY_(status)
#else
    status = NF_CLOSE(unit)
    VERIFY_(status)

    END IF ! MAPL_am_I_root

    CALL MPI_Bcast(gcO3%mncv, SIZE(gcO3%mncv), MPI_REAL, 0, comm, status)
    VERIFY_(status)
    CALL MPI_Bcast(gcO3%mnpl, SIZE(gcO3%mnpl), MPI_REAL, 0, comm, status)
    VERIFY_(status)
#endif

   END IF ChangeOfDay
   
   ALLOCATE(Pclim(gcO3%nlatsPCHEM,gcO3%nlevsPCHEM),STAT=status)
   VERIFY_(status)
   ALLOCATE(Lclim(gcO3%nlatsPCHEM,gcO3%nlevsPCHEM),STAT=status)
   VERIFY_(status)

   ALLOCATE(P1(i1:i2,gcO3%nlevsPCHEM),STAT=status)
   VERIFY_(status)
   ALLOCATE(L1(i1:i2,gcO3%nlevsPCHEM),STAT=status)
   VERIFY_(status)
   
   ALLOCATE(Pi(i1:i2,j1:j2,km),STAT=status)
   VERIFY_(status)
   ALLOCATE(Li(i1:i2,j1:j2,km),STAT=status)
   VERIFY_(status)

   Pclim(:,:) = gcO3%mnpl(:,:,1,1)*fac + gcO3%mnpl(:,:,1,2)*(1.00-fac)
   Lclim(:,:) = gcO3%mnpl(:,:,2,1)*fac + gcO3%mnpl(:,:,2,2)*(1.00-fac)

   DO j=j1,j2
    DO k=1,gcO3%nlevsPCHEM
     CALL MAPL_Interp( P1(:,k), w_c%grid%lat(:,j), Pclim(:,k), gcO3%lats)
     CALL MAPL_Interp( L1(:,k), w_c%grid%lat(:,j), Lclim(:,k), gcO3%lats)
    END DO
    DO i=i1,i2
     CALL MAPL_Interp( Pi(i,j,:), plPa(i,j,:), P1(i,:), gcO3%levs)
     CALL MAPL_Interp( Li(i,j,:), plPa(i,j,:), L1(i,:), gcO3%levs)
    END DO
   END DO

! Turn off in the troposphere
! ---------------------------
   ALLOCATE(mask(i1:i2,j1:j2,1:km),STAT=status)
   VERIFY_(status)
   
   mask = 0

   DO k=1,km
    WHERE(plPa(:,:,k) <= tropp(:,:)) mask(:,:,k) = 1
    WHERE(tropp(:,:) == MAPL_UNDEF)  mask(:,:,k) = 0
   END DO
  
   n = w_c%reg%i_O3
   WHERE(mask == 1) w_c%qa(n)%data3d = (w_c%qa(n)%data3d + cdt*Pi)/(1.00 + cdt*Li)
   
   DEALLOCATE(mask, STAT=status)
   VERIFY_(status)
   DEALLOCATE(Pclim, STAT=status)
   VERIFY_(status)
   DEALLOCATE(Lclim, STAT=status)
   VERIFY_(status)
   DEALLOCATE(Pi, STAT=status)
   VERIFY_(status)
   DEALLOCATE(Li, STAT=status)
   VERIFY_(status)
   DEALLOCATE(P1, STAT=status)
   VERIFY_(status)
   DEALLOCATE(L1, STAT=status)
   VERIFY_(status)

  RETURN
  END SUBROUTINE doProdLoss

!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: DeposVelo
!

      !INTERFACE:
      SUBROUTINE DeposVelo(npts, ij, radiat, tempk, suncos, f0, hstar, xmw,  &
                           ustar, cz1, obk, cfrac, lsnow, rhoa, dvel, rc)

      IMPLICIT NONE

!     ----------------------
!     Argument declarations.
!     ----------------------

      !ARGUMENTS:
      INTEGER, INTENT(IN) :: npts, ij
      REAL, INTENT(IN)    :: radiat(npts)
      REAL, INTENT(IN)    :: tempk (npts)
      REAL, INTENT(IN)    :: suncos(npts)
      REAL, INTENT(IN)    :: f0
      REAL, INTENT(IN)    :: hstar
      REAL, INTENT(IN)    :: xmw
      REAL, INTENT(IN)    :: ustar (npts) ! Cannot be identically zero
      REAL, INTENT(IN)    :: cz1   (npts)
      REAL, INTENT(IN)    :: obk   (npts)
      REAL, INTENT(IN)    :: cfrac (npts)
      INTEGER, INTENT(IN) :: lsnow (npts)
      REAL, INTENT(IN)    :: rhoa  (npts)
      REAL, INTENT(OUT)   :: dvel  (npts)
      INTEGER, INTENT(OUT):: rc

!DESCRIPTION:
!   This routine computes the dry deposition velocities using a
!   resistance-in-series model.
!
!   Routine reads data which:
!   \begin{itemize}
!     \item converts land type id to deposition surface type id
!     \item gives roughness heights for each land type id
!     \item identifies water land type id's, for stability and z0 calculations
!     \item reads surface resistance data for each deposition surface type id
!   \end{itemize}
!
!   Changes from version 3.1 to version 3.2:
!   \begin{itemize}
!     \item  In unstable atmospheres with |zlmo| < zo, as can happen occasionally
!       under very low wind conditions with tall canopies, application of
!       Monin-Obukhov similarity yields negative values for ra.  This was a
!       problem in version 3.1.  In fact, Monin-Obukhov similarity does not
!       apply under such conditions, so we now set ra to zero and let the
!       boundary resistance rb define the overall aerodynamic resistance.
!       Since rb varies inversely with u*, it will impose a large aerodynamic
!       resistance under very low wind conditions.
!     \item The range of applicability of stability correction functions to
!       Monin-Obukhov similarity has been extended to -2.5 < z/zmo < 1.5,
!       based on Figure 2 of Businger et al. [1971].  The range used to be
!       -1 < z/zmo < 1 in version 3.1.
!   \end{itemize}
!
!  Literature cited:
!  \begin{enumerate}
!   \item Baldocchi, D.D., B.B. Hicks, and P. Camara, A canopy stomatal
!       resistance model for gaseous deposition to vegetated surfaces,
!       Atmos. Environ. 21, 91-101, 1987.
!   \item Brutsaert, W., Evaporation into the Atmosphere, Reidel, 1982.
!     Businger, J.A., et al., Flux-profile relationships in the atmospheric
!       surface layer, J. Atmos. Sci., 28, 181-189, 1971.
!   \item Dwight, H.B., Tables of integrals and other mathematical data,
!       MacMillan, 1957.
!   \item Guenther, A., and 15 others, A global model of natural volatile
!       organic compound emissions, J. Geophys. Res., 100, 8873-8892, 1995.
!   \item Hicks, B.B., and P.S. Liss, Transfer of SO2 and other reactive
!       gases across the air-sea interface, Tellus, 28, 348-354, 1976.
!   \item Jacob, D.J., and S.C. Wofsy, Budgets of reactive nitrogen,
!       hydrocarbons, and ozone over the Amazon forest during the wet season,
!       J.  Geophys. Res., 95, 16737-16754, 1990.
!   \item Jacob, D.J., and 9 others, Deposition of ozone to tundra,
!       J. Geophys. Res., 97, 16473-16479, 1992.
!   \item Levine, I.N., Physical Chemistry, 3rd ed., McGraw-Hill, New York, 1988.
!   \item Munger, J.W., and 8 others, Atmospheric deposition of reactive
!       nitrogen oxides and ozone in a temperate deciduous forest and a
!       sub-arctic woodland, J. Geophys. Res., in press, 1996.
!   \item Walcek, C.J., R.A. Brost, J.S. Chang, and M.L. Wesely, SO2, sulfate,
!       and HNO3 deposition velocities computed using regional landuse and
!       meteorological data, Atmos. Environ., 20, 949-964, 1986.
!   \item Wang, Y.H., paper in preparation, 1996.
!   \item Wesely, M.L, Improved parameterizations for surface resistance to
!       gaseous dry deposition in regional-scale numerical models,
!       Environmental Protection Agency Report EPA/600/3-88/025,
!       Research Triangle Park (NC), 1988.
!   \item Wesely, M.L., same title, Atmos. Environ., 23, 1293-1304, 1989.
!  \end{enumerate}
!
!EOP
!
! ARGUMENTS
!   npts    : (i2-i1+1) from AGCM
!   ij      : Current value of 2nd index for ireg: range j1 <= ij <= j2
!   radiat  : solar radiation (W*m^-2)
!   tempk   : surface air temperature (degK)
!   suncos  : cosine of solar zenith angle
!   f0      : reactivity factor for oxidation of biological substances
!   hstar   : Henry's Law constant
!   xmw     : molecular weights; used to calculate molecular diffusivities
!             (kg/mol)
!   airosol : 0 => gas-phase species, 1 => aerosol species (0 for O3)
!   ustar   : friction velocity (m/s)
!   cz1     : altitude at which deposition velocity is computed (m)
!   obk     : Monin-Obukhov length, set to 1.0E5 under neutral conditions (m). Cannot be zero.
!   cfrac   : fractional cloud cover
!   lsnow   : integer for snow and sea ice (1=>water, 2=>land, 3=>ice)
!   dvel    : deposition velocities (m/s)
!   ireg    : # of landtypes in grid square
!   iland   : land type id for elements ldt = 1, ireg;
!             could be from any source, mapped to deposition surface id
!   iuse    : fraction of gridbox area occupied by land type elements (mil^-1)
!   xlai    : leaf area index of land type elements
!   delh_298_over_r : temperature scaling term for Henry's law (GMI has 0.00 for O3)
!
!-----------------------------------------------------------------------------

!     -----------------------
!     Parameter declarations.
!     -----------------------

      REAL, PARAMETER :: XCKMAN = 0.40
      REAL, PARAMETER :: PRESS = 1.50E+05

!     ----------------------
!     Variable declarations.
!     ----------------------

      LOGICAL :: is_found
      LOGICAL :: lrgera(npts)  ! T -> stable atmosphere; a high aerodynamic
                               ! resistance (ra = 1.0d4 m/s) is imposed; else
                               ! ra is calculated

      INTEGER :: idep1
      INTEGER :: ijloop
      INTEGER :: iolson
      INTEGER :: iw
      INTEGER :: k
      INTEGER :: ldt
      INTEGER :: nwaterx

      REAL :: airden
      REAL :: c1
      REAL :: ckustr
      REAL :: corr1
      REAL :: cz       ! altitude where deposition velocity is computed (m)
      REAL :: dair
      REAL :: dummy1, dummy2, dummy3, dummy4
      REAL :: md
      REAL :: ra, rb
      REAL :: rdc
      REAL :: reyno
      REAL :: schno
      REAL :: stono
      REAL :: ebr, eim, ein, r1
      REAL :: rix
      REAL :: rt
      REAL :: tempc1   ! surface air temperatures in degC
      REAL :: tempk1   ! surface air temperatures in degK
      REAL :: xnu
      REAL :: z0obk

      REAL :: c1x   ! total resistance to deposition for
          	    ! each species (s/m)
      REAL :: vd  
      REAL :: vk  

      REAL :: henry_eff

      REAL :: rac (gcO3%NTYPE)
      REAL :: rclo(gcO3%NTYPE)
      REAL :: rcls(gcO3%NTYPE)
      REAL :: rgso(gcO3%NTYPE)
      REAL :: rgss(gcO3%NTYPE)
      REAL :: ri  (gcO3%NTYPE)
      REAL :: rlu (gcO3%NTYPE)

      REAL :: zo  (gcO3%NTYPE)  ! roughness height for specific surface types (m)

      REAL :: rsurfc(gcO3%NTYPE)  ! bulk surface resistance

!     ----------------
!     Begin execution.
!     ----------------

      rc = 0

      dvel(1:npts) = 0.00

!     --------------------------------------------------
!     In each grid cell ...
!     --------------------------------------------------

      IJLoopX: DO ijloop = 1, npts

        cz = cz1(ijloop)
	airden = rhoa(ijloop)

        tempk1 = tempk(ijloop)
        tempc1 = tempk1 - MAPL_TICE
	vd = 0.00
        rsurfc(1:gcO3%NTYPE) = 0.00

!       ---------------------------------------------------------------------
!       Calculate the kinematic viscosity xnu (m^2/s) of air as a function of
!       temperature.  The kinematic viscosity is used to calculate the
!       roughness heights over water surfaces and to diagnose whether such
!       surfaces are aerodynamically rough or smooth using a Reynolds number
!       criterion.  The expression for the temperature dependence of xnu is
!       from the Fortran code in Appendix II of Wesely [1988].
!       ---------------------------------------------------------------------

        c1  = tempk1 / MAPL_TICE
        xnu = 0.151 * (c1**1.77) * 1.00E-04

!       -----------------------------------------------------------------
!       Compute bulk surface resistance for gases.
!
!       Adjust external surface resistances for temperature;
!       from Wesely [1989], expression given in text on p. 1296.
!
!       There is no evidence that the resistance continues to increase at
!       temperatures below -18 C, so at colder temperatures, hold the
!       resistance fixed.
!       -----------------------------------------------------------------

        rt = 1000.00 * EXP (-tempc1 - 4.00)
        IF (tempc1 < -18.00) rt = 1.20E+09

!       ------------------------------------------------------------------
!       Get surface resistances - loop over land types ldt.
!
!       The land types within each grid square are defined using the Olson
!       land type database.  Each of the Olson land types is assigned a
!       corresponding "deposition land type" with characteristic values of
!       surface resistance components.  There are 74 Olson land-types, but
!       only 11 deposition land types (i.e., many of the Olson land types
!       share the same deposition characteristics).  Surface resistance
!       components for the "deposition land types" are from Wesely [1989],
!       except for tropical forests [Jacob and Wofsy, 1990] and for tundra
!       [Jacob et. al., 1992].  All surface resistance components are
!       normalized to a leaf area index of unity.
!
!       Olson land types, deposition land types, and surface resistance
!       components are read from the dry deposition file; check that file
!       for further details.
!       ------------------------------------------------------------------

        LDTLoop1: DO ldt = 1, gcO3%ireg(ijloop,ij)

          IF (gcO3%iuse(ijloop,ij,ldt) == 0) CYCLE LDTLoop1

          IF (lsnow(ijloop) == 3) THEN  ! snow or ice
            idep1  = 1
          ELSE
            iolson = gcO3%iland(ijloop,ij,ldt) + 1
            idep1  = IDEP(iolson)
          END IF

          CALL SurfaceResist(idep1, rac(ldt), rclo(ldt), rcls(ldt), rgso(ldt),  &
                             rgss(ldt), ri(ldt), rlu(ldt), rt, tempc1,  &
                             cfrac(ijloop), radiat(ijloop), suncos(ijloop),  &
                             gcO3%xlai(ijloop,ij,ldt), rix)

!         ----------------------------------------------------------------
!         Compute aerodynamic resistance to lower elements in lower part
!         of the canopy or structure, assuming level terrain; equation (5)
!         of Wesely [1989].
!         ----------------------------------------------------------------

          rdc = 100.00 * (1.00 + (1000.00 / (radiat(ijloop) + 10.00)))

              henry_eff = hstar

!             -----------------------------------------------------
!             Species-dependent corrections to resistances are from
!             equations (6)-(9) of Wesely [1989].
!             -----------------------------------------------------

              CALL CanopyResist(rdc, rix, airden, tempk1, f0, henry_eff, &
                                xmw, rac(ldt), rclo(ldt), rcls(ldt), rgso(ldt), &
                                rgss(ldt), rlu(ldt), rsurfc(ldt))

!           ----------------------------------------------------
!           Set min and max values for bulk surface resistances.
!           ----------------------------------------------------

            rsurfc(ldt) =  MAX (1.00, MIN (rsurfc(ldt), 9999.00))

        END DO LDTLoop1

        LDTLoop3: DO ldt = 1, gcO3%ireg(ijloop,ij)

!         -------------------------------------------------------
!         Loop through the different landuse types present in the
!         grid square.
!         -------------------------------------------------------

          IF (gcO3%iuse(ijloop,ij,ldt) == 0) CYCLE LDTLoop3

          iolson = gcO3%iland(ijloop,ij,ldt) + 1

!         -------------------------------------------------------------
!         Get roughness heights; they are specified constants for each
!         surface type, except over water where zo = f(u*).  The latter
!         dependence is from equation (6) of Hicks and Liss [1976].
!         -------------------------------------------------------------

          nwaterx = gcO3%NWATER
          is_found = .FALSE.

          IWLoop: DO iw = 1, nwaterx

            IF (iolson /= iwater(iw)) cycle IWLoop

            zo(ldt) = 1.40E-02 * ustar(ijloop) * ustar(ijloop) / 9.80 +  &
                      1.10E-01 * xnu / ustar(ijloop)

            is_found = .TRUE.

            EXIT IWLoop

          END DO IWLoop

          IF (.NOT. is_found) THEN
            zo(ldt) = izo(iolson)
            zo(ldt) = zo(ldt) * 1.00E-04
          END IF

!         ------------------------------------------------------------------
!         Get aerodynamic resistances ra and rb.
!
!         The aerodynamic resistance ra is integrated from altitude z0+d up
!         to the altitude z1, at which the dry deposition velocity is to be
!         referenced.  The integration corrects for stability using
!         Monin-Obukhov similarity formulas from Businger et. al. [1971],
!         which apply over the range -2.5 < z/zMO < 1.5  (see their
!         Figure 2).  Under very unstable conditions when z1 > -2.5 zMO, we
!         assume that there is no resistance to transfer in the convective
!         column between zmo and z1.  Under very stable conditions when
!         z1 > 1.5 zMO, we assume that vertical transfer in the column
!         between zmo and z1 is strongly suppressed so that the deposition
!         velocity at altitude z1 is very low.  Under these conditions, we
!         just specify a very large ra=1.0d4 s/m (lrgera = T).
!
!         The Reynolds number reyno diagnoses whether a surface is
!         aerodynamically rough (reyno > 10) or smooth.
!         NOTE: The criterion "reyno > 10" is now replaced by "reyno > 1".

!         Surface is rough in all cases except over water with low wind
!         speeds.  In the smooth case, vertical transport IN THE SUBLAYER
!         near the surface is limited by molecular diffusion and is
!         therefore very slow; we assign a large value of ra + rb to
!         account for this effect. (djj, hyl, bmy, 5/8/00)
!
!         In the aerodynamically rough case, the expression for ra is as
!         given in equation (5) of Jacob et. al. [1992]:
!           ra = (1/ku*)*int(from z0 to z1) (phi(x)/z)dz
!         where x = (z-d)/zmo, z is the height above ground, and d is the
!         displacement height, which is typically 70-80% of the canopy
!         height [Brutsaert, 1982].  We change the vertical coordinate so
!         that z=0 at the displacement height; that's OK since for all
!         practical applications, z1 >> d.  In this manner, we don't need to
!         assume any specific value for the displacement height.  Applying
!         the variable transformation z -> x = z/zmo, the equation above
!         becomes:
!           ra = (1/ku*)*int(from x0 to x1) (phi(x)/x)dx   with x=z/zmo
!         Here phi is a stability correction function originally formulated
!         by Businger et. al. [1971] and given in eqns 5a and 5b of Jacob
!         et. al. [1992].
!         For unstable conditions:
!           phi(x) = a/sqrt(1-bx)  where a=0.74, b = 9
!         The analytical solution to the integral is [Dwight, 1957, integral
!         192.11]:
!           int(dx/(x*sqrt(1-bx))) = log(abs((sqrt(1-bx)-1)/(sqrt(1-bx)+1)))
!         which yields the expression for ra used in the code for unstable
!         conditions.
!         For stable conditions,
!           phi(x) = a + bx        where a=0.74, b = 4.7
!         and the analytical solution to the integral is:
!           int((a/x)+b)dx = a*ln(x) + bx
!         which yields the expression of ra used in the code for stable
!         conditions.
!
!         The formulation of rb for gases is equation (12) of Walcek et. al. [1986].
!         --------------------------------------------------------------------------

          ckustr = XCKMAN * ustar(ijloop)

          reyno = ustar(ijloop) * zo(ldt) / xnu

          IF(obk(ijloop) == 0.00) THEN
	   PRINT *,"Obukhov length scale cannot be zero."
	   rc = 1
	   RETURN
	  END IF

          corr1 = cz / obk(ijloop)

          z0obk = zo(ldt) / obk(ijloop)

          lrgera(ijloop) = .FALSE.

          IF (corr1 >   1.50) lrgera(ijloop) = .TRUE.
          IF (corr1 <= -2.50) corr1 = -2.50

!         ------------------------------------------------------------------
!         Aerodynamically rough or smooth surface.
!
!         In the classic study by Nikuradse (1933) the transition from
!         smooth to rough was examined in pipe flow.  He introduced a
!         roughness Reynolds number rr = u* z0 / nu, and found the flow to
!         be smooth for rr < 0.13, and rough for rr > 2.5, with a transition
!         regime in between (E.B. Kraus and J.A. Businger, Atmosphere-Ocean
!         Interaction, second edition, p.144-145, 1994).  Similar statements
!         can be found in the books:
!           Evaporation into the Atmosphere, by Wilfried Brutsaert, p.59,89,
!             1982; or
!           Seinfeld & Pandis, p.858, 1998.
!         Here we assume a sudden transition point rr = 1 from smooth to
!         rough, following L. Merlivat (1978, The dependence of bulk
!         evaporation coefficients on air-water interfacial conditions as
!         determined by the isotopic method, J. Geophys. Res.,
!         Oceans & Atmos., 83, C6, 2977-2980).  Also refer to Brutsaert's
!         book, p.125.  we used to use the criterion "reyno > 10" for
!         aerodynamically rough surface and now change to "reyno > 1".
!         (hyl, 10/15/99)
!         ------------------------------------------------------------------

          RgeOne: IF (reyno >= 1.00) THEN

!           ------------------------------
!           Aerodynamically rough surface.
!           ------------------------------

            IF ((corr1 <= 0.00) .AND. (z0obk < -1.00)) THEN

!             ----------------------------------------
!             Unstable condition with Abs (zlmo) < Z0;
!             set ra to zero (version 3.2).
!             ----------------------------------------

              ra = 0.00

            ELSE IF ((corr1 <= 0.00) .and. (z0obk >= -1.00)) THEN

!             ---------------------------------------------------
!             Unstable conditions; compute ra as described above.
!             ---------------------------------------------------

              dummy1 = (1.00 - (9.00 * corr1))**0.50
              dummy2 = (1.00 - (9.00 * z0obk))**0.50

	      IF(ABS(dummy1) == 1.00 .OR. ABS(dummy2) .EQ. 1.00) THEN
	       ra = 0.00
	      ELSE
               dummy3 = ABS ((dummy1 - 1.00) / (dummy1 + 1.00))
               dummy4 = ABS ((dummy2 - 1.00) / (dummy2 + 1.00))
               ra = 0.740* (1.00 / ckustr) * LOG (dummy3 / dummy4)
	      END IF

            ELSE IF ((corr1 > 0.00) .AND. (.NOT. lrgera(ijloop))) THEN

!             -----------------------------------------
!             Moderately stable conditions (z/zmo < 1);
!             compute ra as described above
!             -----------------------------------------
              IF(ckustr == 0.00 .OR. corr1 == 0.00 .OR. z0obk == 0.00) THEN
	       rc = 2
               PRINT *,"Cannot calculate ra for moderately stable conditions"
	       RETURN
	      END IF

              ra = (1.00 / ckustr) * (0.740 * LOG (corr1 / z0obk) +  &
                   4.70  * (corr1 - z0obk))

            ELSE IF (lrgera(ijloop)) THEN

!             -----------------------
!             Very stable conditions.
!             -----------------------

              ra = 1.00E+04

            END IF

!             ------------------------------------------------------------
!             If ra is negative (as occasionally happened in version 3.1),
!             set it to zero
!             ------------------------------------------------------------

            ra = MAX(0.00,ra)

!             ------------------------------------
!             Get total resistance for deposition.
!             ------------------------------------

!               ------------------------------------------------------------
!               dair is the thermal diffusivity of air; value of
!               0.2*1.E-4 m^2/s, cited on p. 16,476 of Jacob et. al. [1992].
!               ------------------------------------------------------------

                dair = 0.20 * 1.00E-04
		CALL molDiff(tempk1, airden, xmw, md)

                IF(ckustr == 0.00 .OR. md == 0.00) THEN
	         rc = 3
		 PRINT *,"Zero denominator for rb"
	         RETURN
	        END IF

                rb = (2.00 / ckustr) * (dair / md)**0.667
                c1x = ra + rb + rsurfc(ldt)

          ELSE

!           -------------------------------
!           Aerodynamically smooth surface.
!           -------------------------------

!             --------------------------------------------------------------
!             Suppress drydep over smooth surfaces by setting ra to a
!             large value.  This prevents negative dry deposition
!             velocities when ustar is very small (djj, bmy, 5/8/00)
!             --------------------------------------------------------------

              ra = 1.00E+04
              c1x = ra + rsurfc(ldt)

          END IF RgeOne

!         ----------------------------------------------------------------
!         iuse is the fraction of the grid square occupied by surface ldt
!         in units of mil^-1 (iuse=500 => 50% of the grid square).  Add
!         the contribution of surface type ldt to the deposition velocity;
!         this is a loop over all surface types in the gridbox.
!         ----------------------------------------------------------------

          vd = vd + (0.001 * gcO3%iuse(ijloop,ij,ldt) / c1x)

        END DO LDTLoop3

        dvel(ijloop) = vd

      END DO IJLoopX

      RETURN
      END SUBROUTINE DeposVelo

!-----------------------------------------------------------------------------
!
! ROUTINE
!   CanopyResist
!
! DESCRIPTION
!   This routine calculates bulk surface resistance of the canopy from the
!   network of resistances in parallel and in series.
!
! ARGUMENTS
!   rdc     : tbd
!   rix     :
!   airden  :
!   tempk1  :
!   f01     :
!   hstar1  :
!   xmw1    :
!   rac1    :
!   rclo1   :
!   rcls1   :
!   rgso1   :
!   rgss1   :
!   rlu1    :
!   rsurfc1 :
!
!-----------------------------------------------------------------------------

      SUBROUTINE CanopyResist(rdc, rix, airden, tempk1, f01, hstar1, &
                              xmw1, rac1, rclo1, rcls1, rgso1, rgss1, &
                              rlu1, rsurfc1)

      IMPLICIT NONE

!     ----------------------
!     Argument declarations.
!     ----------------------

      REAL, INTENT(IN)  :: rdc
      REAL, INTENT(IN)  :: rix
      REAL, INTENT(IN)  :: airden
      REAL, INTENT(IN)  :: tempk1
      REAL, INTENT(IN)  :: f01
      REAL, INTENT(IN)  :: hstar1
      REAL, INTENT(IN)  :: xmw1
      REAL, INTENT(IN)  :: rac1
      REAL, INTENT(IN)  :: rclo1
      REAL, INTENT(IN)  :: rcls1
      REAL, INTENT(IN)  :: rgso1
      REAL, INTENT(IN)  :: rgss1
      REAL, INTENT(IN)  :: rlu1
      REAL, INTENT(OUT) :: rsurfc1

!     -----------------------
!     Parameter declarations.
!     -----------------------

      REAL, PARAMETER :: BIG1 = 1.00E+22 

!     ----------------------
!     Variable declarations.
!     ----------------------

      REAL :: dtmp1, dtmp2, dtmp3, dtmp4
      REAL :: mdw, mdx
      REAL :: rclx
      REAL :: rgsx
      REAL :: rixx
      REAL :: rluxx

!     ----------------
!     Begin execution.
!     ----------------

      CALL molDiff (tempk1, airden, MAPL_H2OMW, mdw)
      CALL molDiff (tempk1, airden,       xmw1, mdx)

      rixx = rix * (mdw/mdx) + 1.00/(hstar1/3000.00 + 100.00 * f01)

      IF (rlu1 < 9999.00) then
        rluxx = rlu1 / (hstar1*1.00E-05 + f01)
      ELSE
        rluxx = BIG1
      END IF

!     ---------------------------------------------------------------------
!     To prevent virtually zero resistance to species with huge HSTAR, such
!     as HNO3, a minimum value of rluxx needs to be set.  The rationality
!     of the existence of such a minimum is demonstrated by the observed
!     relationship between Vd(NOy-NOx) and ustar in Munger et. al. [1996];
!     Vd(HNO3) never exceeds 2 cm*s^-1 in observations.  The corresponding
!     minimum resistance is 50 s*m^-1.  This correction was introduced by
!     J.Y. Liang on 7/9/95.
!     ---------------------------------------------------------------------

      rluxx = MAX (rluxx, 50.00)

      rgsx = 1.00 / (hstar1*1.00E-05 / rgss1 + f01 / rgso1)
      rclx = 1.00 / (hstar1*1.00E-05 / rcls1 + f01 / rclo1)

!     -----------------------------------------------------------------------
!     Get the bulk surface resistance of the canopy, rsurfc, from the network
!     of resistances in parallel and in series (Fig. 1 of Wesely [1989]).
!     -----------------------------------------------------------------------

      dtmp1 = 1.00 / rixx
      dtmp2 = 1.00 / rluxx
      dtmp3 = 1.00 / (rac1 + rgsx)
      dtmp4 = 1.00 / (rdc + rclx)

      rsurfc1 = 1.00 / (dtmp1 + dtmp2 + dtmp3 + dtmp4)

      RETURN
      END SUBROUTINE CanopyResist

!-----------------------------------------------------------------------------
!
! ROUTINE
!   SurfaceResist
!
! DESCRIPTION
!   This routine tbd
!
! ARGUMENTS
!   idep1   : 
!   rac1    :
!   rclo1   :
!   rcls1   :
!   rgso1   :
!   rgss1   :
!   ri1     :
!   rlu1    :
!   rt      :
!   tempc1  :
!   cfrac1  :
!   radiat1 :
!   suncos1 :
!   xlai1   :
!   rix     :
!
!-----------------------------------------------------------------------------

      SUBROUTINE SurfaceResist(idep1, rac1, rclo1, rcls1, rgso1, rgss1, &
                               ri1, rlu1, rt, tempc1, cfrac1, radiat1, &
                               suncos1, xlai1, rix)

      IMPLICIT NONE

!     ----------------------
!     Argument declarations.
!     ----------------------

      INTEGER, INTENT(IN) :: idep1
      REAL, INTENT(OUT) :: rac1
      REAL, INTENT(OUT) :: rclo1
      REAL, INTENT(OUT) :: rcls1
      REAL, INTENT(OUT) :: rgso1
      REAL, INTENT(OUT) :: rgss1
      REAL, INTENT(OUT) :: ri1
      REAL, INTENT(OUT) :: rlu1
      REAL, INTENT(IN)  :: rt
      REAL, INTENT(IN)  :: tempc1
      REAL, INTENT(IN)  :: cfrac1
      REAL, INTENT(IN)  :: radiat1
      REAL, INTENT(IN)  :: suncos1
      REAL, INTENT(IN)  :: xlai1
      REAL, INTENT(OUT) :: rix

!     -----------------------
!     Parameter declarations.
!     -----------------------

      REAL, PARAMETER :: BIG1 = 1.00E+06
      REAL, PARAMETER :: BIG2 = 1.00E+22

!     ----------------------
!     Variable declarations.
!     ----------------------

      REAL :: gfaci, gfact, biofit
      REAL :: xdrycoeff(gcO3%NPOLY)

!     ----------------
!     Begin execution.
!     ----------------

!     ----------------------------------------------------------------------
!     Read the internal resistance ri (minimum stomatal resistance for water
!     vapor, per unit area of leaf) from the IRI array; a "9999" value means
!     no deposition to stomata, so we impose a very large value for ri.
!     ----------------------------------------------------------------------

      ri1 = iri(idep1)

      IF (ri1 >= 9999.00) ri1 = BIG2

!     ----------------------------------------------------------------------
!     Cuticular resistances IRLU read in from drydep table are per unit area
!     of leaf; divide them by the leaf area index to get a cuticular
!     resistance for the bulk canopy.  If IRLU is "9999", it means there are
!     no cuticular surfaces on which to deposit, so we impose a very large
!     value for rlu.
!     ----------------------------------------------------------------------

      IF ((irlu(idep1) >= 9999) .OR. (xlai1 <= 0.00)) THEN

        rlu1 = BIG1

      ELSE

        rlu1 = irlu(idep1)
        rlu1 = (rlu1 / xlai1) + rt

      END IF

!     ----------------------------------------------------------
!     The following are the remaining resistances for the Wesely
!     resistance-in-series model for a surface canopy
!     (see Atmos. Environ. paper, Fig.1).
!     ----------------------------------------------------------

      rac1  = MAX (irac(idep1), 1)
      if (rac1  >= 9999.00) rac1  = BIG2

      rgss1 = irgss(idep1)
      rgss1 = MAX (rgss1+rt, 1.00)

      rgso1 = irgso(idep1)
      rgso1 = MAX (rgso1+rt, 1.00)
      if (rgso1 >= 9999.00) rgso1 = BIG2

      rcls1 = ircls(idep1)
      rcls1 = rcls1 + rt
      if (rcls1 >= 9999.00) rcls1 = BIG2

      rclo1 = irclo(idep1)
      rclo1 = rclo1 + rt
      if (rclo1 >= 9999.00) rclo1 = BIG2

!     ----------------------------------------------------------------------
!     Adjust stomatal resistances for insolation and temperature =>
!
!     Temperature adjustment is from Wesely [1989], equation (3).
!
!     Light adjustment by SUBROUTINE lightCorr is described by Wang [1996].
!     It combines:
!       - local dependence of stomal resistance on the intensity I of light
!         impinging the leaf; this is expressed as a mutliplicative
!         factor I/(I+b) to the stomatal resistance where b = 50 W*m^-2
!         (equation (7) of Baldocchi et. al. [1987]);
!       - radiative transfer of direct and diffuse radiation in the
!         canopy using equations (12)-(16) from Guenther et. al. [1995];
!       - separate accounting of sunlit and shaded leaves using
!         equation (12) of Guenther et. al. [1995];
!       - partitioning of the radiation at the top of the canopy into direct
!         and diffuse components using a parameterization to results from
!         an atmospheric radiative transfer model [Wang, 1996].
!     The dependent variables of SUBROUTINE lightCorr are the leaf area
!     index (xylai), the cosine of zenith angle (suncos) and the fractional
!     cloud cover (cfrac).  The factor gfaci integrates the light
!     dependence over the canopy depth; sp even though ri is input per
!     unit area of leaf, it need not be scaled by lai to yield a bulk
!     canopy value because that is already done in the gfaci formulation.
!     ----------------------------------------------------------------------

      IF (ri1 >= 9999.00) THEN

        rix = ri1

      ELSE

        IF ((tempc1 > 0.00) .AND. (tempc1 < 40.00)) THEN
          gfact = 400.00 / tempc1 / (40.00 - tempc1)
        ELSE
          gfact = 100.00
        END IF

        IF ((radiat1 > 1.0E-05) .and. (xlai1 > 0.00)) then
          xdrycoeff(:) = DRYCOEFF(:)
          CALL lightCorr(cfrac1, suncos1, xlai1, xdrycoeff, biofit)
          gfaci = 1.00 /  biofit
        ELSE
          gfaci = 100.00
        END IF

        rix = ri1 * gfact * gfaci

      END IF

      RETURN
      END SUBROUTINE SurfaceResist

!-----------------------------------------------------------------------------
!
! ROUTINE
!   lightCorr
!
! DESCRIPTION
!   This routine calculates the light correction.
!
!   Light adjustment by SUBROUTINE lightCorr is described by Wang [1996].
!   It combines:
!     * local dependence of stomal resistance on the intensity I of light
!       impinging the leaf; this is expressed as a mutliplicative factor
!       I/(I+b) to the stomatal resistance where b = 50 (W*m^-2)
!       (equation (7) of Baldocchi et al. [1987])
!     * radiative transfer of direct and diffuse radiation in the canopy
!       using equations (12)-(16) from Guenther et al. [1995]
!     * separate accounting of sunlit and shaded leaves using equation (12)
!       of Guenther et al. [1995]
!     * partitioning of the radiation at the top of the canopy into direct
!       and diffuse components using a parameterization of results from an
!       atmospheric radiative transfer model [Wang, 1996].
!
! ARGUMENTS
!   cloud_frac1 : fractional cloud cover
!   suncos1     : cosine of the solar zenith angle
!   xlai1       : leaf area index of land type for current month
!   coeff       : factor that integrates the light dependence over the canopy
!                 depth; "sp" even though "ri" is input per unit area of leaf
!
! REVISION HISTORY
!   Unknown date Original from GMIchem
!   Jan 2012     Nielsen, revisions for GOCART O3
!
!-----------------------------------------------------------------------------

      SUBROUTINE lightCorr(cloud_frac1, suncos1, xlai1, coeff, f)

      IMPLICIT NONE

!     ----------------------
!     Argument declarations.
!     ----------------------

      REAL, INTENT(IN) :: cloud_frac1
      REAL, INTENT(IN) :: suncos1
      REAL, INTENT(IN) :: xlai1
      REAL, INTENT(IN) :: coeff(gcO3%NPOLY)
      REAL, INTENT(OUT) :: f

!     -----------------------
!     Parameter declarations.
!     -----------------------

      INTEGER, PARAMETER :: BIOFIT_DIM = 4   ! biofit dimension
      INTEGER, PARAMETER :: SPDIM = BIOFIT_DIM-1

!     ----------------------
!     Variable declarations.
!     ----------------------

      INTEGER :: ii, nn
      INTEGER :: k0, k1, k2, k3

      REAL  :: bigterm(gcO3%NPOLY)
      REAL  :: spterm (BIOFIT_DIM-1)
      REAL  :: term   (BIOFIT_DIM)

!     ----------------
!     Begin execution.
!     ----------------

      term(1) = 1.00
      term(2) = xlai1
      term(3) = suncos1
      term(4) = cloud_frac1

      DO ii = 1, SPDIM
        spterm(ii) = term(ii+1)
      END DO

      CALL ScaleBiofit(SPDIM, spterm)

      DO ii = 1, SPDIM
        term(ii+1) = spterm(ii)
      END DO

      k0 = 0

      DO k3 = 1, BIOFIT_DIM
        DO k2 = k3, BIOFIT_DIM
          DO k1 = k2, BIOFIT_DIM

            k0 = k0 + 1
            bigterm(k0) = term(k1) * term(k2) * term(k3)

     	  END DO
     	END DO
      END DO

      f = 0.00

      DO nn = 1, gcO3%NPOLY
        f = f + (coeff(nn) * bigterm(nn))
      END DO

      f = MAX (f, 0.10)

      RETURN
      END SUBROUTINE lightCorr

!-----------------------------------------------------------------------------
!
! ROUTINE
!   ScaleBiofit
!
! DESCRIPTION
!   This routine scales and constrains xlai, suncos, and cloud_frac; called
!   by lightCorr.
!
! ARGUMENTS
!   spdim  : spterm dimension
!   spterm : array of terms containing xlai, suncos, cloud_frac
!
! REVISION HISTORY
!   Unknown date Original from GMIchem
!   Jan 2012     Nielsen, revisions for GOCART O3
!
!-----------------------------------------------------------------------------

      SUBROUTINE ScaleBiofit(spdim, spterm)

      IMPLICIT NONE

!     ----------------------
!     Argument declarations.
!     ----------------------

      INTEGER, INTENT(IN) :: spdim
      REAL,INTENT(INOUT) :: spterm(spdim)

!     -----------------------
!     Parameter declarations.
!     -----------------------

      INTEGER, PARAMETER :: SPSIZE = 3
      INTEGER :: ND(SPSIZE) = (/ 55, 20, 11 /) ! Scaling factor for each variable
      REAL  ::  XHI(SPSIZE) = (/ 11.00, 1.00, 1.00 /) ! Maximum for each variable

!     ----------------------
!     Variable declarations.
!     ----------------------

      INTEGER :: ii

      REAL :: rnd
      REAL :: xlow ! Minimum for each variable

!     ----------------
!     Begin execution.
!     ----------------

      DO ii = 1, spdim

        spterm(ii) = MIN (spterm(ii), XHI(ii))

        IF (ii /= spdim) THEN

          rnd  = ND (ii)
          xlow = XHI(ii) / rnd

        ELSE

          xlow = 0.00

        END IF

        spterm(ii) = MAX (spterm(ii), xlow)
        spterm(ii) = spterm(ii) / XHI(ii)

      END DO

      RETURN
      END SUBROUTINE ScaleBiofit

!-----------------------------------------------------------------------------
!
! ROUTINE
!   molDiff
!
! DESCRIPTION
!   This routine calculates the molecular diffusivity of a gas in air
!   (m^2/s).
!
!   The molecular radius of air is given in a table on p. 479 of Levine
!   [1988]; the table also gives radii for some other molecules.  Rather
!   than using a specific molecular radius a generic value is used for all
!   molecules, which is good enough in terms of calculating the diffusivity
!   as long as the molecule is not too big.
!
! ARGUMENTS
!   tk     : temperature [K]
!   airden : air density [kg m^{-3}]
!   xm     : molecular weight of gas [kg kmol^{-1}]
!   radx   : hard-sphere molecular radius of the diffusing gas [m]
!
! REVISION HISTORY
!   Unknown date Original from GMIchem
!   Jan 2012     Nielsen, revisions for GOCART O3
!
!-----------------------------------------------------------------------------

      SUBROUTINE molDiff(tk, airden, xm, d)

      IMPLICIT NONE

!     ----------------------
!     Argument declarations.
!     ----------------------

      REAL, INTENT(IN) :: tk
      REAL, INTENT(IN) :: airden
      REAL, INTENT(IN) :: xm
      REAL, INTENT(OUT):: d

!     -----------------------
!     Parameter declarations.
!     -----------------------

      REAL, PARAMETER :: RADAIR = 1.2E-10 ! Hard-sphere molecular radii of air (m)
      REAL, PARAMETER :: RADX   = 1.5E-10 ! Hard-sphere molecular radius of diffusing gas
      REAL, PARAMETER :: PI = 3.14159

!     ----------------------
!     Variable declarations.
!     ----------------------

      REAL :: diam
      REAL :: frpath
      REAL :: speed
      REAL :: zz

!     ----------------
!     Begin execution.
!     ----------------

!     -------------------------------------------------------------------
!     Calculate the mean free path for gas X in air:  eq. 8.5 of Seinfeld
!     [1986]; diam is the collision diameter for gas X with air.
!     -------------------------------------------------------------------

      zz = xm / MAPL_AIRMW

      diam = radx + RADAIR

      frpath = 1.00 / (PI * SQRT (1.00 + zz ) * airden * (diam * diam))

!     -------------------------------------------------------------
!     Calculate average speed of gas X; eq. 15.47 of Levine [1988].
!     -------------------------------------------------------------

      speed = SQRT (8.00 * MAPL_RUNIV * tk / (PI * xm))

!     --------------------------------------------------------------------
!     Calculate diffusion coefficient of gas X in air; eq. 8.9 of Seinfeld
!     [1986].
!     --------------------------------------------------------------------

      d = (3.00 * PI / 32.00) * (1.00 + zz) * frpath * speed

      RETURN
      END SUBROUTINE molDiff

 END SUBROUTINE O3_GridCompRun

!-------------------------------------------------------------------------
!NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS!
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  O3_GridCompFinalize --- The Chem Driver 
!
! !INTERFACE:
!

   SUBROUTINE O3_GridCompFinalize ( gcO3, w_c, impChem, expChem, &
                                    nymd, nhms, cdt, rc )

! !USES:

  IMPLICIT none

! !INPUT/OUTPUT PARAMETERS:

   TYPE(O3_GridComp), INTENT(INOUT) :: gcO3   ! Grid Component

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
!   4Mar2005 Nielsen   Implementation of parameterized ozone chemistry
!  31Jan2012 Nielsen   Add dry deposition and NetCDF reads from PCHEM
!
!EOP
!-------------------------------------------------------------------------

   CHARACTER(LEN=*), PARAMETER :: Iam = 'O3_GridCompFinalize'
   INTEGER :: status
   
   status = 0
   rc = 0

   DEALLOCATE(gcO3%ireg, STAT=status)
   VERIFY_(status)
   DEALLOCATE(gcO3%iland, STAT=status)
   VERIFY_(status)
   DEALLOCATE(gcO3%iuse, STAT=status)
   VERIFY_(status)
   DEALLOCATE(gcO3%xlai, STAT=status)
   VERIFY_(status)
   DEALLOCATE(gcO3%lats, STAT=status)
   VERIFY_(status)
   DEALLOCATE(gcO3%levs, STAT=status)
   VERIFY_(status)
   DEALLOCATE(gcO3%mncv, STAT=status)
   VERIFY_(status)
   DEALLOCATE(gcO3%mnpl, STAT=status)
   VERIFY_(status)

   RETURN

 END SUBROUTINE O3_GridCompFinalize

 END MODULE O3_GridCompMod
