#include "MAPL_Generic.h"
!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 610.1, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  CO2_GridCompMod --- CO2 Grid Component Class
!
! !INTERFACE:
!

   module  CO2_GridCompMod

! !USES:

   USE ESMF
   USE MAPL

   use Chem_Mod              ! Chemistry Base Class
   use Chem_StateMod         ! Chemistry State
   use Chem_UtilMod
   use m_inpak90             ! Resource file management

   implicit none

! !PUBLIC TYPES:
!
   PRIVATE
   PUBLIC  CO2_GridComp       ! The CO2 object 

!
! !PUBLIIC MEMBER FUNCTIONS:
!

   PUBLIC  CO2_GridCompSetServices
   PUBLIC  CO2_GridCompInitialize
   PUBLIC  CO2_GridCompRun
   PUBLIC  CO2_GridCompFinalize

!
! !DESCRIPTION:
!
!  This module implements the (pre-ESMF) CO2 Grid Component. 
!
! !REVISION HISTORY:
!
!  16Sep2003 da Silva  First crack.
!  24Oct2005     Bian  tag CO2 to 4 regions 
!                      (total, north america, south america, africa)
!  19Dec2005 da Silva  Activated 3D diags for output
!  26Nov2010  Nielsen  Simplified PBL partitioning for biomass burning emissions   
!                      
!EOP
!-------------------------------------------------------------------------

  type CO2_GridComp
        character(len=255) :: name

        real               :: BBconFac           ! conversion factor of BB emissions to C
        real               :: FFconFac           ! conversion factor of FF emissions to C
        real               :: BioDrawDownFactor  ! Biosphere drawdown factor

        real,    pointer   ::    eCO2_FF(:,:)    ! kgC/m2/s, Earth surface
        real,    pointer   ::   eCO2_NEP(:,:)    ! kgC/m2/s, Earth surface
        real,    pointer   ::   eCO2_OCN(:,:)    ! kgC/m2/s, Earth surface
        real,    pointer   ::   eCO2_BB_(:,:)    ! kgC/m2/s, PBL (before diurnal)
        real,    pointer   ::    eCO2_BB(:,:)    ! kgC/m2/s, PBL
        real,    pointer   :: regionMask(:,:)    ! regional mask
        real,    pointer   ::        SST(:,:)    ! SST, C, for CMS ocean flux calc
        real,    pointer   ::        SSS(:,:)    ! Salinty, PSU, for CMS ocean flux calc
        real,    pointer   ::       pco2(:,:)    ! Ocean pCO2, for CMS ocean flux calc
        real,    pointer   ::       pice(:,:)    ! % grid covered by ice, for CMS ocean flux calc
        integer, pointer   ::  regionIndex(:)    ! desired regions from mask

        logical            :: DBG                ! Run-time debug switch
        logical            :: CMS_EMIS           ! Run-time switch to use CMS emissions
        logical            :: OCN_FLUX_CALC      ! Method of ocean flux calculation from input file name
  end type CO2_GridComp

! real, parameter :: radToDeg = 57.2957795
  real, parameter :: radToDeg = 180./MAPL_PI
  real, parameter :: mwtCO2   = 44.0098
  real, parameter :: mwtC     = 12.0110

CONTAINS

   subroutine CO2_GridCompSetServices(  gc, chemReg, rc)
   type(ESMF_GridComp), intent(INOUT) :: GC
   type(Chem_Registry), intent(INOUT) :: chemReg
   integer,             intent(OUT  ) :: rc

   CHARACTER(LEN=255) :: rcbasen = 'CO2_GridComp'
   CHARACTER(LEN=255) :: name

   integer            :: status
   character(len=ESMF_MAXSTR) :: Iam

   integer :: ier,doingCMS,ocnFlux

   type(ESMF_Config) :: cfg

   Iam = "CO2_GridCompSetServices"

!  Load resource file
!  ------------------
   cfg = ESMF_ConfigCreate(rc=status)
   VERIFY_(STATUS)
   call ESMF_ConfigLoadFile(cfg,TRIM(rcbasen)//'.rc',rc=status)
   VERIFY_(STATUS)

   call ESMF_ConfigGetAttribute(cfg, value=doingCMS, label='CMS_EMIS:', rc=status)
   VERIFY_(STATUS)
   call ESMF_ConfigGetAttribute(cfg, value=ocnFlux, label='OCN_FLUX_CALC:', default= 0, rc=status)
   VERIFY_(STATUS)

   call MAPL_AddImportSpec(GC,           &
        SHORT_NAME = 'CO2_regionMask',   &
        LONG_NAME  = 'source species'  , &
        UNITS      = '1',                &
        DIMS       = MAPL_DimsHorzOnly,  &
        VLOCATION  = MAPL_VLocationNone, &
        RESTART    = MAPL_RestartSkip,   &
        RC         = STATUS)
   VERIFY_(STATUS)

   call MAPL_AddImportSpec(GC,                            &
       SHORT_NAME = 'QTOT',                                &
       LONG_NAME  = 'mass_fraction_of_all_water',          &
       UNITS      = 'kg kg-1',                             &
       DIMS       = MAPL_DimsHorzVert,                     &
       VLOCATION  = MAPL_VLocationCenter,                  &
       RESTART    = MAPL_RestartSkip,     __RC__)

   if (doingCMS == 0) then 
      call MAPL_AddImportSpec(GC,           &
           SHORT_NAME = 'CO2_BIOMASS',      &
           LONG_NAME  = 'source species'  , &
           UNITS      = 'kg CO2 m-2 s-1',   &
           DIMS       = MAPL_DimsHorzOnly,  &
           VLOCATION  = MAPL_VLocationNone, &
           RESTART    = MAPL_RestartSkip,   &
           RC         = STATUS)
      VERIFY_(STATUS)
      call MAPL_AddImportSpec(GC,           &
           SHORT_NAME = 'CO2_FF',           &
           LONG_NAME  = 'source species'  , &
           UNITS      = 'kg C m-2 s-1',     &
           DIMS       = MAPL_DimsHorzOnly,  &
           VLOCATION  = MAPL_VLocationNone, &
           RESTART    = MAPL_RestartSkip,   &
           RC         = STATUS)
      VERIFY_(STATUS)
      call MAPL_AddImportSpec(GC,           &
           SHORT_NAME = 'CO2_NEP',          &
           LONG_NAME  = 'source species'  , &
           UNITS      = 'kg C m-2 s-1',     &
           DIMS       = MAPL_DimsHorzOnly,  &
           VLOCATION  = MAPL_VLocationNone, &
           RESTART    = MAPL_RestartSkip,   &
           RC         = STATUS)
      VERIFY_(STATUS)
      call MAPL_AddImportSpec(GC,           &
           SHORT_NAME = 'CO2_OCN',          &
           LONG_NAME  = 'source species'  , &
           UNITS      = 'kg C m-2 s-1',     &
           DIMS       = MAPL_DimsHorzOnly,  &
           VLOCATION  = MAPL_VLocationNone, &
           RESTART    = MAPL_RestartSkip,   &
           RC         = STATUS)
      VERIFY_(STATUS)

   else
      call MAPL_AddImportSpec(GC,           &
           SHORT_NAME = 'CO2_CMS_BIOMASS',  &
           LONG_NAME  = 'source species'  , &
           UNITS      = 'kg C m-2 s-1',     &
           DIMS       = MAPL_DimsHorzOnly,  &
           VLOCATION  = MAPL_VLocationNone, &
           RESTART    = MAPL_RestartSkip,   &
           RC         = STATUS)
      VERIFY_(STATUS)
      call MAPL_AddImportSpec(GC,           &
           SHORT_NAME = 'CO2_CMS_FF',       &
           LONG_NAME  = 'source species'  , &
           UNITS      = 'kg C m-2 s-1',     &
           DIMS       = MAPL_DimsHorzOnly,  &
           VLOCATION  = MAPL_VLocationNone, &
           RESTART    = MAPL_RestartSkip,   &
           RC         = STATUS)
      VERIFY_(STATUS)
      call MAPL_AddImportSpec(GC,           &
           SHORT_NAME = 'CO2_CMS_NEP',      &
           LONG_NAME  = 'source species'  , &
           UNITS      = 'kg C m-2 s-1',     &
           DIMS       = MAPL_DimsHorzOnly,  &
           VLOCATION  = MAPL_VLocationNone, &
           RESTART    = MAPL_RestartSkip,   &
           RC         = STATUS)
      VERIFY_(STATUS)

      if (ocnFlux /=0) then
         call MAPL_AddImportSpec(GC,           &
              SHORT_NAME = 'CO2_CMS_T',        &
              LONG_NAME  = 'source species'  , &
              UNITS      = '1',                &
              DIMS       = MAPL_DimsHorzOnly,  &
              VLOCATION  = MAPL_VLocationNone, &
              RESTART    = MAPL_RestartSkip,   &
              RC         = STATUS)
              VERIFY_(STATUS)
         call MAPL_AddImportSpec(GC,           &
              SHORT_NAME = 'CO2_CMS_S',        &
              LONG_NAME  = 'source species'  , &
              UNITS      = '1',                &
              DIMS       = MAPL_DimsHorzOnly,  &
              VLOCATION  = MAPL_VLocationNone, &
              RESTART    = MAPL_RestartSkip,   &
              RC         = STATUS)
         VERIFY_(STATUS)
         call MAPL_AddImportSpec(GC,           &
              SHORT_NAME = 'CO2_CMS_PICE',     &
              LONG_NAME  = 'source species'  , &
              UNITS      = '1',                &
              DIMS       = MAPL_DimsHorzOnly,  &
              VLOCATION  = MAPL_VLocationNone, &
              RESTART    = MAPL_RestartSkip,   &
              RC         = STATUS)
         VERIFY_(STATUS)
         call MAPL_AddImportSpec(GC,           &
              SHORT_NAME = 'CO2_CMS_PCO2',     &
              LONG_NAME  = 'source species'  , &
              UNITS      = '1',                &
              DIMS       = MAPL_DimsHorzOnly,  &
              VLOCATION  = MAPL_VLocationNone, &
              RESTART    = MAPL_RestartSkip,   &
              RC         = STATUS)
         VERIFY_(STATUS)
 
      else
         call MAPL_AddImportSpec(GC,           &
              SHORT_NAME = 'CO2_CMS_OCN',      &
              LONG_NAME  = 'source species'  , &
              UNITS      = 'kg C m-2 s-1',     &
              DIMS       = MAPL_DimsHorzOnly,  &
              VLOCATION  = MAPL_VLocationNone, &
              RESTART    = MAPL_RestartSkip,   &
              RC         = STATUS)
         VERIFY_(STATUS)
      end if
   end if

   RETURN_(ESMF_SUCCESS)

   end subroutine CO2_GridCompSetServices

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  CO2_GridCompInitialize --- Initialize CO2_GridComp
!
! !INTERFACE:
!

   subroutine CO2_GridCompInitialize ( gcCO2, w_c, impChem, expChem, &
                                      nymd, nhms, cdt, rc )

! !USES:

  implicit NONE

! !INPUT PARAMETERS:

   type(Chem_Bundle), intent(in) :: w_c        ! Chemical tracer fields      
   integer, intent(in) :: nymd, nhms           ! time
   real,    intent(in) :: cdt                  ! chemical timestep (secs)


! !OUTPUT PARAMETERS:

   type(CO2_GridComp), intent(inout) :: gcCO2   ! Grid Component
   type(ESMF_State), intent(inout)  :: impChem  ! Import State
   type(ESMF_State), intent(inout)  :: expChem  ! Export State
   integer, intent(out) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 - 

! !DESCRIPTION: Initializes the CO2 Grid Component. It primarily sets
!               the import state for each active constituent package.
!
! !REVISION HISTORY:
!
!  18Sep2003 da Silva  First crack.
!  24OCT2005     Bian  Mods for 5 tagged CO2  
!                      (total, fossil fuel, ecosystem, oceanic, and biomass)
!  25OCT2005     Bian  Mods for 5 regions
!  01SEP2015     Weir  Added dry-air mole fraction export

!
!EOP
!-------------------------------------------------------------------------

   character(len=*), parameter :: myname = 'CO2_GridCompInitialize'


   character(len=255) :: rcfilen = 'CO2_GridComp.rc'
   integer :: ios, n
   integer, allocatable :: ier(:)
   integer :: i1, i2, im, j1, j2, jm, km, ijl
   integer :: nbins, nbeg, nend, nbins_rc, nymd1, nhms1
   integer :: nTimes, begTime, incSecs
   real    :: qmin, qmax

   gcCO2%name = 'CO2 Constituent Package'

!  Initialize local variables
!  --------------------------
   rc = 0
   i1 = w_c%grid%i1; i2 = w_c%grid%i2; im = w_c%grid%im
   j1 = w_c%grid%j1; j2 = w_c%grid%j2; jm = w_c%grid%jm
   km = w_c%grid%km
   nbins = w_c%reg%n_CO2;  nbeg  = w_c%reg%i_CO2; nend  = w_c%reg%j_CO2

   ijl = ( i2 - i1 + 1 ) * ( j2 - j1 + 1 )

   call init_()
   if ( rc /= 0 ) return

   ier(:) =0

!                       -------------------
!                       Parse resource file
!                       -------------------

!  Load resource file
!  ------------------
   call i90_loadf ( TRIM(rcfilen), ier(1) )
   if ( ier(1) .ne. 0 ) then
      call final_(10)
      return
   end if
   ier(:)=0

   call i90_label ( 'number_CO2_bins:', ier(1) )
   nbins_rc = i90_gint ( ier(2) )
   if ( nbins_rc /= nbins ) then
      call final_(11)
      return
   end if

!  Run-time switch to use CMS emissions
!  ------------------------------------
   CALL I90_label ( 'CMS_EMIS:', ier(3) )
   n = I90_gint ( ier(4) )
   IF(n /= 0) THEN
    gcCO2%CMS_EMIS = .TRUE.
    if ( MAPL_AM_I_ROOT() ) then
       print *, 'CO2_GridComp: Reading  C M S  Emissions'
    end if
   ELSE
    gcCO2%CMS_EMIS = .FALSE.
    if ( MAPL_AM_I_ROOT() ) then
       print *, 'CO2_GridComp: Reading CLIMATOLOGICAL Emissions'
    end if
   END IF

   IF (gcCO2%CMS_EMIS) THEN
    CALL I90_label ( 'CMS_biomass_emission_factor:', ier(5) )
    gcCO2%BBconFac = i90_gfloat ( ier(6) )

    CALL I90_label ( 'CMS_fossilfuel_emissions_factor:', ier(7) )
    gcCO2%FFconFac = i90_gfloat ( ier(8) )
   ELSE
    CALL I90_label ( 'CO2_biomass_emission_factor:', ier(5) )
    gcCO2%BBconFac = i90_gfloat ( ier(6) )

    CALL I90_label ( 'CO2_fossilfuel_emissions_factor:', ier(7) )
    gcCO2%FFconFac = i90_gfloat ( ier(8) )
   END IF

! Biosphere drawdown factor, used for adjusting biosphere drawdown.
! Valid range: >0.  If < 1, reduces sink.  If = 1, neutral.  If > 1, enhances sink.
! ---------------------------------------------------------------------------------
   gcCO2%BioDrawDownFactor = -1.00
   CALL I90_label ( 'Biosphere_drawdown_factor:', ier(9) )
   gcCO2%BioDrawDownFactor = I90_gfloat ( ier(10) )
   IF(gcCO2%BioDrawDownFactor < 0.00) THEN
    IF(MAPL_AM_I_ROOT()) THEN
     PRINT *," "
     PRINT *,TRIM(myname)//": Invalid biosphere drawdown factor."
    END IF
    CALL final_(12)
    RETURN
   END IF

   IF(MAPL_AM_I_ROOT()) THEN
    PRINT *," "
    PRINT *,TRIM(myname)//": "
    PRINT *," Biomass emission factor:     ", gcCO2%BBconFac
    PRINT *," Fossil fuel emission factor: ", gcCO2%FFconFac
    PRINT *," Biosphere drawdown factor:   ", gcCO2%BioDrawDownFactor
   END IF

!  Run-time debug switch
!  ---------------------
   CALL I90_label ( 'DEBUG:', ier(11) )
   n = I90_gint ( ier(12) )
   IF(n /= 0) THEN
    gcCO2%DBG = .TRUE.
   ELSE
    gcCO2%DBG = .FALSE.
   END IF

!  Get the desired regions to run on
!  ---------------------------------
   call i90_label ( 'CO2_regions_indices:', ier(13) )
   do n = 1, nbins
      gcCO2%regionIndex(n) = i90_gint ( ier(13+n) )
   end do

   IF( ANY( ier(:) /= 0 ) ) THEN
    CALL final_(13)
    RETURN
   END IF
   ier(:)=0

   call I90_label ( 'OCN_FLUX_CALC:' , ier(1))
   if (ier(1) == 0) then
      n = I90_gint( ier(2) )
   else
      n =0
   end if
   if(n /= 0) then
     gcCO2%OCN_FLUX_CALC=.true.
   else
     gcCO2%OCN_FLUX_CALC=.false.
   end if

   DEALLOCATE(ier)

   return

CONTAINS

   subroutine init_()
   integer ios, nerr
   nerr = max ( 100, nbins+27 )
   allocate (   gcCO2%eCO2_FF(i1:i2,j1:j2), & 
               gcCO2%eCO2_NEP(i1:i2,j1:j2), & 
               gcCO2%eCO2_OCN(i1:i2,j1:j2), & 
                gcCO2%eCO2_BB(i1:i2,j1:j2), & 
               gcCO2%eCO2_BB_(i1:i2,j1:j2), & 
             gcCO2%regionMask(i1:i2,j1:j2), &
                  gcCO2%regionIndex(nbins), &
                    gcCO2%SST(i1:i2,j1:j2), & 
                    gcCO2%SSS(i1:i2,j1:j2), & 
                   gcCO2%pCO2(i1:i2,j1:j2), & 
                   gcCO2%pice(i1:i2,j1:j2), & 
                      ier(nerr), stat=ios )

   if ( ios /= 0 ) rc = 100
   end subroutine init_

   subroutine final_(ierr)
   integer :: ierr
   integer ios
   deallocate ( gcCO2%eCO2_FF, gcCO2%eCO2_NEP, gcCO2%eCO2_OCN, &
                gcCO2%eCO2_BB, gcCO2%regionMask, gcCO2%regionIndex, &
                gcCO2%eCO2_BB_, gcCO2%SST, gcCO2%SSS, gcCO2%pCO2, &
                gcCO2%pice, ier, stat=ios )
   call i90_release()
   rc = ierr
   end subroutine final_

   end subroutine CO2_GridCompInitialize

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  CO2_GridCompRun --- The Chem Driver 
!
! !INTERFACE:
!

   subroutine CO2_GridCompRun ( gcCO2, w_c, impChem, expChem, &
                               nymd, nhms, cdt, rc )

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(CO2_GridComp), intent(inout) :: gcCO2   ! Grid Component
   type(Chem_Bundle),  intent(inout) :: w_c     ! Chemical tracer fields   

! !INPUT PARAMETERS:

   type(ESMF_State), intent(inout) :: impChem   ! Import State
   integer, intent(in) :: nymd, nhms            ! time
   real,    intent(in) :: cdt                   ! chemical timestep (secs)


! !OUTPUT PARAMETERS:

   type(ESMF_State), intent(inout) :: expChem   ! Export State
   integer, intent(out) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 -
 
! !DESCRIPTION: This routine implements the so-called CO2 Driver. That 
!               is, adds chemical tendencies to each of the constituents,
!  Note: water wapor, the first constituent is not considered a chemical
!  constituents.
!
! !REVISION HISTORY:
!
!  18Sep2003 da Silva  First crack.
!  24OCT2005     Bian  Mods for 5 tagged CO2  
!                      (total, fossil fuel, ecosystem, oceanic, and biomass)
!  25OCT2005     Bian  Mods for 5 regions
!  29SEP2017     Weir  Lots of zero diff (hopefully) changes to make
!                      understandable
!  19SEP2019     Weir  More zero diff (hopefully) changes to make
!                      understandable

!
!EOP
!-------------------------------------------------------------------------

   character(len=*), parameter :: myname = 'CO2_GridCompRun'
   character(len=*), parameter :: Iam = myname

!  Input fields from fvGCM
!  -----------------------
   REAL, POINTER, DIMENSION(:,:)   ::  pblh  => null()
   REAL, POINTER, DIMENSION(:,:,:) ::  zle   => null()
   REAL, POINTER, DIMENSION(:,:)   ::  ps    => null()
   REAL, POINTER, DIMENSION(:,:)   ::  v10m  => null()
   REAL, POINTER, DIMENSION(:,:)   ::  u10m  => null()
   REAL, POINTER, DIMENSION(:,:,:) ::  T     => null()
   REAL, POINTER, DIMENSION(:,:,:) ::  qtot  => null()

   integer :: i1, i2, im, j1, j2, jm, km, idiag, ios, ijl
   integer :: i, j, k, n, nbins, nbeg, nend
   INTEGER :: nymd1, nhms1, ier(9)

   REAL    :: qmin, qmax, c2co2

   REAL, POINTER, DIMENSION(:,:)   :: ptr2d => null()
   REAL, POINTER, DIMENSION(:,:,:) :: ptr3d => null()

   REAL, ALLOCATABLE, DIMENSION(:,:) :: psdry, psco2

#define EXPORT     expChem

#define ptrCO2EM      CO2_emis
#define ptrCO2SC      CO2_surface
#define ptrCO2CL      CO2_column
#define ptrCO2DRY     CO2_dry
#define ptrCO2SCDRY   CO2_surfdry
#define ptrCO2CLDRY   CO2_coldry

   integer :: STATUS

#include "CO2_GetPointer___.h"

!  Initialize local variables
!  --------------------------
   rc = 0
   i1 = w_c%grid%i1; i2 = w_c%grid%i2; im = w_c%grid%im
   j1 = w_c%grid%j1; j2 = w_c%grid%j2; jm = w_c%grid%jm
   km = w_c%grid%km
   ijl = ( i2 - i1 + 1 ) * ( j2 - j1 + 1 )

!  Retreive mask
!  -------------
   call MAPL_GetPointer(impChem,ptr2d,'CO2_regionMask',rc=status)
   VERIFY_(STATUS)
   gcCO2%regionMask=ptr2d

   nbins = w_c%reg%n_CO2;  nbeg  = w_c%reg%i_CO2; nend  = w_c%reg%j_CO2

   c2co2 = mwtCO2 / mwtC

   if ( any((/NBIN_CO2CL,NBIN_CO2SC/)/=NBIN_CO2EM)) then
      call die(myname,'all emissions in registry must have same number of bins')
   endif
   if ( nbins > NBIN_CO2EM ) then
      call die(myname,'nbins in chem registry must be <= those in component registry')
   end if

!  Read climatological emissions, if specified
!  -------------------------------------------
   IF ( .NOT. gcCO2%CMS_EMIS ) THEN

!     Biomass burning
!     ---------------
      call MAPL_GetPointer(impChem, ptr2d, 'CO2_BIOMASS',rc=status)
      VERIFY_(STATUS)
      gcCO2%eCO2_BB = ptr2d

!     Fossil fuel emissions
!     ---------------------
      call MAPL_GetPointer(impChem, ptr2d, 'CO2_FF',rc=status)
      VERIFY_(STATUS)
      gcCO2%eCO2_FF = ptr2d

!     Biosphere flux
!     --------------
      call MAPL_GetPointer(impChem, ptr2d, 'CO2_NEP',rc=status)
      VERIFY_(STATUS)
      gcCO2%eCO2_NEP = ptr2d

!     Ocean flux
!     ----------
      call MAPL_GetPointer(impChem, ptr2d, 'CO2_OCN',rc=status)
      VERIFY_(STATUS)
      gcCO2%eCO2_OCN = ptr2d

!     Bian says that we need to adjust the uptake flux of CO2 in the
!     ecosystem database to reflect the emissions from biomass burning.
!     In principle this adds a factor that needs to be balanced on an
!     interannual basis.  For year 2000 TRMM (GFED v 1.2) emissions this
!     factor is 1.2448
!     bweir: WHA?!?
!     ------------------------------------------------------------------
      where (gcCO2%eCO2_NEP(i1:i2,j1:j2) .le. 0.0)
          gcCO2%eCO2_NEP(i1:i2,j1:j2) = gcCO2%eCO2_NEP(i1:i2,j1:j2)*gcCO2%BioDrawDownFactor
      end where

   ELSE ! TYPE OF EMISS IS CMS			     

      call MAPL_GetPointer(impChem, ptr2d, 'CO2_CMS_BIOMASS',rc=status)
      VERIFY_(STATUS)
      gcCO2%eCO2_BB = ptr2d
      
      call MAPL_GetPointer(impChem, ptr2d, 'CO2_CMS_FF',rc=status)
      VERIFY_(STATUS)
      gcCO2%eCO2_FF = ptr2d

      call MAPL_GetPointer(impChem, ptr2d, 'CO2_CMS_NEP',rc=status)
      VERIFY_(STATUS)
      gcCO2%eCO2_NEP = ptr2d

      if ( gcCO2%OCN_FLUX_CALC ) then
         call MAPL_GetPointer(impChem, ptr2d, 'CO2_CMS_T',rc=status)
         VERIFY_(STATUS)
         gcCO2%sst = ptr2d
         call MAPL_GetPointer(impChem, ptr2d, 'CO2_CMS_S',rc=status)
         VERIFY_(STATUS)
         gcCO2%sss = ptr2d
         call MAPL_GetPointer(impChem, ptr2d, 'CO2_CMS_PICE',rc=status)
         VERIFY_(STATUS)
         gcCO2%pice = ptr2d
         call MAPL_GetPointer(impChem, ptr2d, 'CO2_CMS_PCO2',rc=status)
         VERIFY_(STATUS)
         gcCO2%pco2 = ptr2d
      else
         call MAPL_GetPointer(impChem, ptr2d, 'CO2_CMS_OCN',rc=status)
         VERIFY_(STATUS)
         gcCO2%eCO2_OCN = ptr2d
      end if

   END IF ! type of emiss  

!  Apply conversion factors
!  ------------------------
   gcCO2%eCO2_BB(i1:i2,j1:j2) = gcCO2%eCO2_BB(i1:i2,j1:j2) * gcCO2%BBconFac
   gcCO2%eCO2_FF(i1:i2,j1:j2) = gcCO2%eCO2_FF(i1:i2,j1:j2) * gcCO2%FFconFac

!  Units for surface flux must be kgCO2 m^-2 s^-1
!  ----------------------------------------------
   gcCO2%eCO2_NEP(i1:i2,j1:j2) = gcCO2%eCO2_NEP(i1:i2,j1:j2) * c2co2
   gcCO2%eCO2_OCN(i1:i2,j1:j2) = gcCO2%eCO2_OCN(i1:i2,j1:j2) * c2co2
   gcCO2%eCO2_FF(i1:i2,j1:j2)  = gcCO2%eCO2_FF(i1:i2,j1:j2)  * c2co2
!  Mimicking original code: only convert biomass burning for CMS emissions case
   if ( gcCO2%CMS_EMIS ) then
      gcCO2%eCO2_BB(i1:i2,j1:j2) = gcCO2%eCO2_BB(i1:i2,j1:j2) * c2co2
   end if

!  Apply diurnal cycle if desired
!  ------------------------------
   if ( w_c%diurnal_bb ) then
      gcCO2%eCO2_BB_(:,:) = gcCO2%eCO2_BB(:,:)

      call Chem_BiomassDiurnal ( gcCO2%eCO2_BB, gcCO2%eCO2_BB_,   &
                                 w_c%grid%lon(:,:)*radToDeg,      &
                                 w_c%grid%lat(:,:)*radToDeg, nhms, cdt )      
   end if

!  Get imports
!  -----------
   call MAPL_GetPointer ( impChem, pblh, 'ZPBL',    rc=ier(1) ) ; VERIFY_(ier(1))
   call MAPL_GetPointer ( impChem, zle,  'ZLE',     rc=ier(4) ) ; VERIFY_(ier(4))
   call MAPL_GetPointer ( impChem, T,    'T',       rc=ier(2) ) ; VERIFY_(ier(2))
   call MAPL_GetPointer ( impChem, u10m, 'U10M',    rc=ier(5) ) ; VERIFY_(ier(5))
   call MAPL_GetPointer ( impChem, v10m, 'V10M',    rc=ier(6) ) ; VERIFY_(ier(6))
   call MAPL_GetPointer ( impChem, ps,   'PS',      rc=ier(7) ) ; VERIFY_(ier(7))
   call MAPL_GetPointer ( impChem, qtot, 'QTOT',    rc=ier(8) ) ; VERIFY_(ier(8))

   if (gcCO2%DBG) then
      call pmaxmin('CO2: e_ff',  gcCO2%eCO2_FF,  qmin, qmax, ijl, 1, 1. )
      call pmaxmin('CO2: e_nep', gcCO2%eCO2_NEP, qmin, qmax, ijl, 1, 1. )
      call pmaxmin('CO2: e_ocn', gcCO2%eCO2_OCN, qmin, qmax, ijl, 1, 1. )
      call pmaxmin('CO2: e_bb',  gcCO2%eCO2_BB,  qmin, qmax, ijl, 1, 1. )

      call pmaxmin('CO2: pblh', pblh, qmin, qmax, ijl, 1,    1. )
      call pmaxmin('CO2: ps',   ps,   qmin, qmax, ijl, 1,    1. )
      call pmaxmin('CO2: u10m', u10m, qmin, qmax, ijl, 1,    1. )
      call pmaxmin('CO2: v10m', v10m, qmin, qmax, ijl, 1,    1. )
      call pmaxmin('CO2: T',    T,    qmin, qmax, ijl, km,   1. )
      call pmaxmin('CO2: zle',  zle,  qmin, qmax, ijl, km+1, 1. )
      call pmaxmin('CO2: qtot', qtot, qmin, qmax, ijl, km,   1. )
   end if

!  Read and apply CO2 Emissions
!  ----------------------------
   call CO2_Emission(rc)

!  Fill the export states
!  ----------------------

!  Surface concentration [ppmv]
!  ----------------------------
   do n = 1,nbins
      if (associated(CO2_surface(n)%data2d)) then
         CO2_surface(n)%data2d(i1:i2,j1:j2) = w_c%qa(nbeg+n-1)%data3d(i1:i2,j1:j2,km)*1.00E6
      endif
   enddo

!  Column burden [kg m-2]
!  ----------------------
   do n = 1,nbins
      if (associated(CO2_column(n)%data2d)) then
         CO2_column(n)%data2d(i1:i2,j1:j2) = 0.
         do k = 1,km
           CO2_column(n)%data2d(i1:i2,j1:j2) = CO2_column(n)%data2d(i1:i2,j1:j2) &
                 + w_c%qa(nbeg+n-1)%data3d(i1:i2,j1:j2,k)*mwtCO2/MAPL_AIRMW      &
                            * w_c%delp(i1:i2,j1:j2,k)/MAPL_GRAV
         enddo
      endif
   enddo

!  Dry-air mole fraction [mol mol-1]
!  ---------------------------------
   do n = 1,nbins
      if (associated(CO2_dry(n)%data3d)) then
         CO2_dry(n)%data3d(i1:i2,j1:j2,1:km) = w_c%qa(nbeg+n-1)%data3d(i1:i2,j1:j2,1:km) &
                                             / (1. - qtot(i1:i2,j1:j2,1:km))
      endif
   enddo

!  Dry-air surface concentration [mol mol-1]
!  -----------------------------------------
   do n = 1,nbins
      if (associated(CO2_surfdry(n)%data2d)) then
         CO2_surfdry(n)%data2d(i1:i2,j1:j2) = w_c%qa(nbeg+n-1)%data3d(i1:i2,j1:j2,km)*(1. - qtot(i1:i2,j1:j2,km))
      endif
   enddo

!  Dry-air column average [mol mol-1]
!  ----------------------------------
   allocate(psdry(i1:i2,j1:j2), stat=ios)    ! dry-air surface pressure
   allocate(psco2(i1:i2,j1:j2), stat=ios)    ! co2     surface pressure

   do n = 1,nbins
      if (associated(CO2_coldry(n)%data2d)) then
         psdry(i1:i2,j1:j2) = 0.
         psco2(i1:i2,j1:j2) = 0.
         do k = 1,km
            psdry(i1:i2,j1:j2) = psdry(i1:i2,j1:j2) + w_c%delp(i1:i2,j1:j2,k)*(1. - qtot(i1:i2,j1:j2,k))
            psco2(i1:i2,j1:j2) = psco2(i1:i2,j1:j2) + w_c%delp(i1:i2,j1:j2,k)*w_c%qa(nbeg+n-1)%data3d(i1:i2,j1:j2,k)
         enddo
         CO2_coldry(n)%data2d(i1:i2,j1:j2) = psco2(i1:i2,j1:j2)/psdry(i1:i2,j1:j2)
      endif
   enddo

   deallocate(psdry, psco2)

!  Fill the export state with current mixing ratios
!  ------------------------------------------------
   do n = 1,nbins
      if (associated(CO2BIN(n)%data3d)) then
         CO2BIN(n)%data3d(i1:i2,j1:j2,1:km) = w_c%qa(nbeg+n-1)%data3d(i1:i2,j1:j2,1:km)
      endif
   enddo

   return

CONTAINS

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
! !IROUTINE:  CO2_Emission - Adds emissions for CO2 for one timestep
!             We have emissions from 4 sources, which are distributed
!             differently in the vertical
!             1) fossil fuel - emitted at surface
!             2) ecosystem   - fluxes at surface
!             3) oceanic     - fluxes at surface
!             4) biomass burning - uniformly mixed in PBL
!
! !DESCRIPTION: Updates the CO2 concentration with emissions every timestep
!
! !REVISION HISTORY:
!
!  24Oct2005, Bian
!  26Nov2010, Nielsen Simplified PBL partitioning for biomass burning emissions 
!  
! !INTERFACE:
!
!EOP
!-------------------------------------------------------------------------
     SUBROUTINE CO2_Emission ( rc )
!-------------------------------------------------------------------------

     IMPLICIT NONE

!    OUTPUT VARIABLES
     INTEGER, INTENT(OUT) :: rc ! Error return code

!    LOCAL VARIABLES
     CHARACTER(LEN=*), PARAMETER :: myname = 'CO2_Emission'

     INTEGER ::  i, j, k, kt, minkPBL
     INTEGER, ALLOCATABLE :: index(:)

     REAL, ALLOCATABLE :: pblLayer(:,:),sfcFlux(:,:),myMask(:,:)
     REAL, ALLOCATABLE :: fPBL(:,:,:)

     real :: scco2, scco2arg,wssq,rkwco2,tk,tk100,tk1002,ff
     real :: ffuatm,xco2,deltco2,wspd,flxmolm2

     rc = 0

     ALLOCATE(sfcFlux(i1:i2,j1:j2),STAT=ios)        ! emissions
     ALLOCATE(myMask(i1:i2,j1:j2),STAT=ios)         ! region mask
     ALLOCATE(fPBL(i1:i2,j1:j2,1:km),STAT=ios)      ! partitioning of BB

!    Find the layer that contains the PBL.
!    Layer thicknesses are ZLE(:,:,0:km).
!    -------------------------------------
     ALLOCATE(index(0:km),STAT=ios)
     ALLOCATE(pblLayer(i1:i2,j1:j2),STAT=ios)
     DO j=j1,j2
       DO i=i1,i2
        index(0:km)=0
        WHERE(zle(i,j,0:km)-zle(i,j,km) > pblh(i,j)) index(0:km)=1
        pblLayer(i,j)=SUM(index)
       END DO
     END DO
     minkPBL=MINVAL(pblLayer)

!    Determine partitioning fraction based on layer thicknesses
!    ----------------------------------------------------------
     fPBL(i1:i2,j1:j2,1:km)=0.00
     DO j=j1,j2
       DO i=i1,i2
        kt=pblLayer(i,j)
        DO k=kt,km
         fPBL(i,j,k)=(zle(i,j,k-1)-zle(i,j,k))/(zle(i,j,kt-1)-zle(i,j,km))
        END DO
       END DO
     END DO

!    Calculate NOBM fluxes using archived ocean fields.
!    --------------------------------------------------
     IF (gcCO2%OCN_FLUX_CALC ) then
        if ( gcCO2%sst(i,j) .lt. 1000.) THEN
         scco2 = 2073.1 - 125.62*gcCO2%sst(i,j) + 3.6276*gcCO2%sst(i,j)**2. - &
                 0.043219*gcCO2%sst(i,j)**3.
         scco2arg = (scco2/660.0)**(-0.5)
         wspd = (u10m(i,j)**2 + v10m(i,j)**2)**0.5
         wssq = wspd*wspd
         rkwco2 = wssq*scco2arg*0.337/(3.6E5) 
         tk = gcCO2%sst(i,j) + 273.15
         tk100 = tk*0.01
         tk1002 = tk100*tk100
         ff = exp(-162.8301 + 218.2968/tk100  +                       &
                 90.9241*log(tk100) - 1.47696*tk1002 +                &
                 gcCO2%sss(i,j) * (.025695 - .025225*tk100 +          &
                 0.0049867*tk1002))   
         ffuatm = ff*1.0E-6                    
         xco2 = w_c%qa(nbeg)%data3d(i,j,km)*ps(i,j)*1.0E4/1013.25
         deltco2 = (xco2-gcCO2%pco2(i,j))*ffuatm*1024.5   ! mol/m3
         flxmolm2 = rkwco2*deltco2                        ! units of mol/m2/s
         gcCO2%eCO2_OCN(i,j) = -flxmolm2*mwtC*1E-3*(100.-gcCO2%pice(i,j))/100. 
        ELSE  ! SST must be undef
         gcCO2%eCO2_OCN(i,j) = 0.   
        END IF 
     END IF

     DEALLOCATE(index,STAT=ios)
     DEALLOCATE(pblLayer,STAT=ios)

     CO2_Bin: DO n=1,nbins

!      Finalize the mask. Update CO2 globally if the region index is -1.
!      Otherwise update only where the mask's value is the region index.
!      -----------------------------------------------------------------
       IF (gcCO2%regionIndex(n) == -1) THEN
         myMask(i1:i2,j1:j2)=gcCO2%regionIndex(n)
       ELSE
         myMask(i1:i2,j1:j2)=gcCO2%regionMask(i1:i2,j1:j2)
       END IF

!      Establish range of layers on which to work
!      ------------------------------------------
       kt = minkPBL

       Layer: do k=kt,km

!        Emissions: Weighted biomass burning
!        -----------------------------------
         sfcFlux(i1:i2,j1:j2) = gcCO2%eCO2_BB(i1:i2,j1:j2)*fPBL(i1:i2,j1:j2,k)

!        Add Fossil fuel, net ecosystem production, and ocean source when in surface layer
!        ---------------------------------------------------------------------------------
         if (k == km) sfcFlux(i1:i2,j1:j2) = sfcFlux(i1:i2,j1:j2)        &
                                           + gcCO2%eCO2_FF(i1:i2,j1:j2)  &
                                           + gcCO2%eCO2_NEP(i1:i2,j1:j2) &
                                           + gcCO2%eCO2_OCN(i1:i2,j1:j2)

!        Update CO2 for this bin
!        -----------------------
         where (myMask(i1:i2,j1:j2) == gcCO2%regionIndex(n))
            w_c%qa(nbeg+n-1)%data3d(i1:i2,j1:j2,k) = w_c%qa(nbeg+n-1)%data3d(i1:i2,j1:j2,k) &
                  + cdt * sfcFlux(i1:i2,j1:j2)*MAPL_AIRMW/mwtCO2                            &
                        / (w_c%delp(i1:i2,j1:j2,k)/MAPL_GRAV)
         end where
       end do Layer

!      Update Surface flux diagnostic for this bin
!      -------------------------------------------
       if (ASSOCIATED(CO2_emis(n)%data2d)) then
         CO2_emis(n)%data2d(i1:i2,j1:j2) = 0.

         sfcFlux(i1:i2,j1:j2) = gcCO2%eCO2_FF(i1:i2,j1:j2)  + gcCO2%eCO2_NEP(i1:i2,j1:j2) &
                              + gcCO2%eCO2_OCN(i1:i2,j1:j2) + gcCO2%eCO2_BB(i1:i2,j1:j2)

         where (myMask(i1:i2,j1:j2) == gcCO2%regionIndex(n))
           CO2_emis(n)%data2d(i1:i2,j1:j2) = sfcFlux(i1:i2,j1:j2)
         end where
       end if
     end do CO2_Bin

     DEALLOCATE(fPBL,STAT=ios)
     DEALLOCATE(myMask,STAT=ios)
     DEALLOCATE(sfcFlux,STAT=ios)

     RETURN

     END SUBROUTINE CO2_Emission

   END SUBROUTINE CO2_GridCompRun

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  CO2_GridCompFinalize --- The Chem Driver 
!
! !INTERFACE:
!

   subroutine CO2_GridCompFinalize ( gcCO2, w_c, impChem, expChem, &
                                    nymd, nhms, cdt, rc )

! !USES:

  implicit NONE

! !INPUT/OUTPUT PARAMETERS:

   type(CO2_GridComp), intent(inout) :: gcCO2   ! Grid Component

! !INPUT PARAMETERS:

   type(Chem_Bundle), intent(in)  :: w_c      ! Chemical tracer fields   
   integer, intent(in) :: nymd, nhms          ! time
   real,    intent(in) :: cdt                 ! chemical timestep (secs)


! !OUTPUT PARAMETERS:

   type(ESMF_State), intent(inout) :: impChem   ! Import State
   type(ESMF_State), intent(inout) :: expChem   ! Import State
   integer, intent(out) ::  rc                  ! Error return code:
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

   character(len=*), parameter :: myname = 'CO2_GridCompFinalize'
   INTEGER :: ios

   DEALLOCATE ( gcCO2%eCO2_FF, gcCO2%eCO2_NEP, gcCO2%eCO2_OCN, &
                gcCO2%eCO2_BB, gcCO2%regionMask, gcCO2%SST, &
                gcCO2%SSS, gcCO2%pCO2, gcCO2%pice, STAT=ios )
   rc = 0
   IF ( ios /= 0 ) rc = 1

   return

 end subroutine CO2_GridCompFinalize

 end module CO2_GridCompMod

