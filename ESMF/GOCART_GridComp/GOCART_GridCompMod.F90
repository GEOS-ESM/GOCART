#include "MAPL_Generic.h"

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1     !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: GOCART_GridCompMod - The GOCART Aerosol Grid Component
!
! !INTERFACE:
!
   Module GOCART_GridCompMod
!
! !USES:
!
   use ESMF
   use MAPL

   use Chem_Mod              ! Chemistry Base Class
   use Chem_UtilMod, only: Chem_UtilNegFiller
   use Aero_GridCompMod      ! Parent Aerosol component with IRF methods but no SetServices()

   use ConvectionMod, only: Disable_Convection

   implicit none
   private

   type(Chem_Mie), dimension(2), save :: gocartMieTable
   integer, parameter :: instanceComputational = 1
   integer, parameter :: instanceData          = 2

   character(len=*), parameter :: H2O2_RECYCLE_ALARM = 'GOCART::RECYCLE_H2O2'
   character(len=*), parameter :: HNO3_RECYCLE_ALARM = 'GOCART::RECYCLE_HNO3'

!
! !PUBLIC MEMBER FUNCTIONS:

   public SetServices
!
! !DESCRIPTION: 
!
!   {\tt GOCART} is a gridded component from the GOCART model and includes 
!  dust, sea salt, sulfates, organic and black carbon. In addition, we
!  also include closely related components for CO and CO2 with relatively
!  simple parameterization of the chemical processes, but sharing
!  consistent emissions with the aerosols.
!
!  This code derives from the pre-ESMF Chem component from GEOS-4. This
!  GEOS-4 Chem "component" used ESMF like constructs (Chem component class, 
!  import/export states, etc) but no ESMF specific data types because of 
!  an odd incompatibility with the fvGCM code (the so-called 
!  {\tt oldworld} library. Unlike GEOS-4, the Stratospheric Chemistry
!  component is treated separately here.
!
! !REVISION HISTORY:
!
!  25feb2005  da Silva  First crack.
!  19jul2006  da Silva  First separate GOCART component.
!
!EOP
!-------------------------------------------------------------------------

  type GOCART_State
     private
     type(Chem_Registry), pointer :: chemReg => null()
     type(Aero_GridComp), pointer :: gcChem  => null()
     type(Chem_Bundle), pointer   :: w_c     => null()
     logical                      :: data_driven = .false.
   end type GOCART_State

  type GOCART_WRAP
     type (GOCART_State), pointer :: PTR => null()
  end type GOCART_WRAP

CONTAINS


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: SetServices --- Sets IRF services for GOCART Grid Component
!
! !INTERFACE:

   subroutine SetServices ( GC, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, optional                  :: RC  ! return code

! !DESCRIPTION: Sets Initialize, Run and Finalize services. 
!
! !REVISION HISTORY:
!
!  25feb2005  da Silva  First crack.
!
!EOP
!-------------------------------------------------------------------------


!   ErrLog Variables
!   ----------------
    character(len=ESMF_MAXSTR)      :: IAm = 'SetServices'
    integer                         :: STATUS
    character(len=ESMF_MAXSTR)      :: COMP_NAME

!   Local derived type aliases
!   --------------------------
    type (ESMF_Config)              :: CF
    type (GOCART_State), pointer    :: state   ! internal, that is
    type (GOCART_wrap)              :: wrap
    type(Chem_Registry), pointer    :: r

    integer                         :: n, nq
    integer                         :: DO_CO2CNNEE
    character(len=ESMF_MAXSTR)      :: FRIENDLIES
    character(len=ESMF_MAXSTR)      :: AEROFRIENDLY
    character(len=ESMF_MAXSTR)      :: providerName
    character(len=ESMF_MAXSTR)      :: short_name
    real                            :: DEFVAL
    real                            :: DEFVAL_CO2

    character(len=ESMF_MAXSTR)      :: field_name
    character(len=ESMF_MAXSTR)      :: chem_registry_file

!                              ------------

    
!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, CONFIG=CF, __RC__ )

    Iam = TRIM(COMP_NAME) // '::' // 'SetServices'

!   Wrap internal state for storing in GC; rename legacyState
!   -------------------------------------
    allocate ( state, __STAT__ )
    wrap%ptr => state

!   Is the component data driven
!   ----------------------------
    state%data_driven = IsDataDrivenGC_(GC, __RC__)

!   Start by loading the Chem Registry
!   ----------------------------------
    allocate ( state%chemReg )

    if (state%data_driven) then
        state%chemReg = Chem_RegistryCreate(STATUS, rcfile='GOCARTdata_AerRegistry.rc')
        VERIFY_(STATUS)
    else
       call ESMF_ConfigGetAttribute(cf, chem_registry_file, label = "Chem_Registry_File:", &
            default = "Chem_Registry.rc", rc = status)
       VERIFY_(status)
       state%chemReg = Chem_RegistryCreate(STATUS, rcfile=chem_registry_file)
        VERIFY_(STATUS)
    end if    

    r => state%chemReg   ! short hand
    

!                       ------------------------
!                       ESMF Functional Services
!                       ------------------------

!   Set the Initialize, Run, Finalize entry points
!   ----------------------------------------------
     if ( r%doing_GOCART ) then

        if(MAPL_AM_I_ROOT()) then
         print *, trim(Iam)//': ACTIVE'
         print *,' '
        end if

        call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_INITIALIZE,  Initialize_, __RC__ )
        call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,         Run1_,       __RC__ )
        call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,         Run2_,       __RC__ )

        call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_FINALIZE,    Finalize_,   __RC__ )
        
!       Store internal state in GC
!       --------------------------
        call ESMF_UserCompSetInternalState ( GC, 'GOCART_state', wrap, STATUS )
        VERIFY_(STATUS)

     else

        if (MAPL_AM_I_ROOT()) then
         print *, trim(Iam)//': NOT ACTIVE, defaulting to Generic No-op stubs'
        end if 

        call MAPL_GenericSetServices ( GC, __RC__ )
        RETURN_(ESMF_SUCCESS)

     endif


!                         ------------------
!                         GEOS Data Services
!                         ------------------

!  NOTE: For now, always define import state to avoid breaking connectivities.

!!BOS
!
! !IMPORT STATE:

!    Pressure at layer edges
!    -----------------------
    call MAPL_AddImportSpec(GC,                            &
       SHORT_NAME = 'PLE',                                 &
       LONG_NAME  = 'air_pressure',                        &
       UNITS      = 'Pa',                                  &
       DIMS       = MAPL_DimsHorzVert,                     &
       VLOCATION  = MAPL_VLocationEdge,                    &
       RESTART    = MAPL_RestartSkip,     __RC__) 

!   Pressure thickness
!   ------------------
    call MAPL_AddImportSpec(GC,                            &
       SHORT_NAME = 'DELP',                                &
       LONG_NAME  = 'pressure_thickness',                  &
       UNITS      = 'Pa',                                  &
       DIMS       = MAPL_DimsHorzVert,                     &
       VLOCATION  = MAPL_VLocationCenter,                  &
       RESTART    = MAPL_RestartSkip,     __RC__)

!   Height at the edges
!   -------------------
    call MAPL_AddImportSpec(GC,                            &
       SHORT_NAME = 'ZLE',                                 &
       LONG_NAME  = 'geopotential_height',                 &
       UNITS      = 'm',                                   &
       DIMS       = MAPL_DimsHorzVert,                     &
       VLOCATION  = MAPL_VLocationEdge,                    &
       RESTART    = MAPL_RestartSkip,     __RC__)

!    AIRDENS: moist air density
!    --------------------------
     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME = 'AIRDENS',                            &
        LONG_NAME  = 'moist_air_density',                  &
        UNITS      = 'kg/m^3',                             &
        DIMS       = MAPL_DimsHorzVert,                    &
        VLOCATION  = MAPL_VLocationCenter,                 &
        RESTART    = MAPL_RestartSkip,    __RC__)

!    AIRDENS_DRY: dry air density
!    ----------------------------
     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME = 'AIRDENS_DRYP',                       &
        LONG_NAME  = 'partial_dry_air_density',            &
        UNITS      = 'kg dry m-3 tot',                     &
        DIMS       = MAPL_DimsHorzVert,                    &
        VLOCATION  = MAPL_VLocationCenter,                 &
        RESTART    = MAPL_RestartSkip,     __RC__)

!    CLOUD
!    -----
     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME = 'FCLD'  ,                             &
        LONG_NAME  = 'Cloud fraction for radiation',       &
        UNITS      = '1',                                  &
        DIMS       = MAPL_DimsHorzVert,                    &
        VLOCATION  = MAPL_VLocationCenter,                 &
        RESTART    = MAPL_RestartSkip,     __RC__)


!   T
!   -
    call MAPL_AddImportSpec(GC,                            &
       SHORT_NAME = 'T',                                   &
       LONG_NAME  = 'air_temperature',                     &
       UNITS      = 'K',                                   &
       DIMS       = MAPL_DimsHorzVert,                     &
       VLOCATION  = MAPL_VLocationCenter,                  &
       RESTART    = MAPL_RestartSkip,     __RC__)

!   U
!   -
    call MAPL_AddImportSpec(GC,                            &
       SHORT_NAME = 'U',                                   &
       LONG_NAME  = 'eastward_wind',                       &
       UNITS      = 'm s-1',                               &
       DIMS       = MAPL_DimsHorzVert,                     &
       VLOCATION  = MAPL_VLocationCenter,                  &
       RESTART    = MAPL_RestartSkip,     __RC__)

!   V
!   -
    call MAPL_AddImportSpec(GC,                            &
       SHORT_NAME = 'V',                                   &
       LONG_NAME  = 'northward_wind',                      &
       UNITS      = 'm s-1',                               &
       DIMS       = MAPL_DimsHorzVert,                     &
       VLOCATION  = MAPL_VLocationCenter,                  &
       RESTART    = MAPL_RestartSkip,     __RC__)

!   Q
!   -
    call MAPL_AddImportSpec(GC,                            &
       SHORT_NAME = 'Q',                                   &
       LONG_NAME  = 'specific_humidity',                   &
       UNITS      = 'kg kg-1',                             &
       DIMS       = MAPL_DimsHorzVert,                     &
       VLOCATION  = MAPL_VLocationCenter,                  &
       RESTART    = MAPL_RestartSkip,     __RC__)  

!   QITOT + QITOT
!   -------------
    call MAPL_AddImportSpec(GC,                            &
       SHORT_NAME = 'QCTOT',                               &
       LONG_NAME  = 'mass_fraction_of_total_cloud_water',  &
       UNITS      = 'kg kg-1',                             &
       DIMS       = MAPL_DimsHorzVert,                     &
       VLOCATION  = MAPL_VLocationCenter, __RC__)

!   RH: is between 0 and 1
!   ----------------------
    call MAPL_AddImportSpec(GC,                            &
       SHORT_NAME = 'RH2',                                 &
       LONG_NAME  = 'Rel_Hum_after_moist',                 &
       UNITS      = '1',                                   &
       DIMS       = MAPL_DimsHorzVert,                     &
       VLOCATION  = MAPL_VLocationCenter,                  &
       RESTART    = MAPL_RestartSkip,     __RC__)

!   PFI_LAAN: this is nonconvective precipition
!   ----------------------------------------
     call MAPL_AddImportSpec(GC,                           &
         SHORT_NAME='PFI_LSAN',                            &
         LONG_NAME ='3D_flux_of_ice_nonconvective_precipitation',  &
         UNITS     ='kg/m2/s',                             &
         DIMS      = MAPL_DimsHorzVert,                    &
         VLOCATION = MAPL_VLocationEdge,                   &
         RESTART   = MAPL_RestartSkip,    __RC__)

!    PFL_LAAN: this is nonconvective precipition
!    ----------------------------------------
     call MAPL_AddImportSpec(GC,                           &
         SHORT_NAME='PFL_LSAN',                            &
         LONG_NAME ='3D_flux_of_liquid_nonconvective_precipitation',  &
         UNITS     ='kg/m2/s',                             &
         DIMS      = MAPL_DimsHorzVert,                    &
         VLOCATION = MAPL_VLocationEdge,                   &
         RESTART   = MAPL_RestartSkip,    __RC__)

!   Ozone from PCHEM for CFC-12 photolysis
!   --------------------------------------
  if (state%chemReg%doing_CFC) then
    call MAPL_AddImportSpec(GC,                            &
       SHORT_NAME = 'O3',                                  &
       LONG_NAME  = 'ozone_mass_mixing_ratio',             &
       UNITS      = 'kg/kg',                               &
       DIMS       = MAPL_DimsHorzVert,                     &
       VLOCATION  = MAPL_VLocationCenter,                  &
       RESTART    = MAPL_RestartSkip,     __RC__)
  end if

!    2-D Quantities
!    --------------    

!    TROPP - Connectivity from SDYN to PHYS is TROPP_BLENDED to TROPP
!    ----------------------------------------------------------------
     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME = 'TROPP',                              &
        LONG_NAME  = 'tropopause_pressure_based_on_blended_estimate', &
        UNITS      = 'Pa',                                 &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationNone,                   &
        RESTART    = MAPL_RestartSkip,    __RC__)

!    LWI
!    ---
     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME = 'LWI',                                &
        LONG_NAME  = 'land-ocean-ice_mask',                &
        UNITS      = '1',                                  &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationNone,                   &
      ! RESTART    = MAPL_RestartSkip,                     &
                                          __RC__)

!    PBL 
!    ---
     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME = 'ZPBL',                               &
        LONG_NAME  = 'Planetary boundary layer height',    &
        UNITS      = 'm',                                  &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationNone,                   &
        RESTART    = MAPL_RestartSkip,    __RC__)

!    FRACLAKE
!    --------
     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME = 'FRLAKE',                             &
        LONG_NAME  = 'fraction_of_lake',                   &
        UNITS      = '1',                                  &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationNone,                   &
        RESTART    = MAPL_RestartSkip,    __RC__)

!    FROCEAN
!    -------
     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME = 'FROCEAN',                            &
        LONG_NAME  = 'fraction_of_ocean',                  &
        UNITS      = '1',                                  &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationNone,                   &
        RESTART    = MAPL_RestartSkip,    __RC__)

!    FRACI
!    -----
     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME = 'FRACI',                              &
        LONG_NAME  = 'ice_covered_fraction_of_tile',       &
        UNITS      = '1',                                  &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationNone,                   &
      ! RESTART    = MAPL_RestartSkip,                     &
                                          __RC__)

!    GWETTOP
!    -------
     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME = 'WET1',                               &
        LONG_NAME  = 'surface_soil_wetness',               &
        UNITS      = '1',                                  &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationNone, __RC__)

!    LAI
!    ---
     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME = 'LAI',                                &
        LONG_NAME  = 'leaf_area_index',                    &
        UNITS      = '1',                                  &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationNone,                   &
        RESTART    = MAPL_RestartSkip,   __RC__)

!    This could be useful,  but it is not needed now
!    -----------------------------------------------
     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME = 'GRN',                                &
        LONG_NAME  = 'greeness_fraction',                  &
        UNITS      = '1',                                  &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationNone,                   &
        RESTART    = MAPL_RestartSkip,   __RC__)

!    PRECC: I hope this is defined over oceans
!    -----------------------------------------
     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME = 'CN_PRCP',                            &
        LONG_NAME  = 'Surface Conv. rain flux needed by land', &
        UNITS      = 'kg/m^2/s',                           &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationNone,                   &
        RESTART    = MAPL_RestartSkip,   __RC__)

!    PRECL: Non-convective precip, provided by Cinderella
!    ----------------------------------------------------
     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME = 'NCN_PRCP',                           &
        LONG_NAME  = 'Non-convective precipitation',       &
        UNITS      = 'kg/m^2/s',                           &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationNone,                   &
        RESTART    = MAPL_RestartSkip,   __RC__)

!    PS: from where???
!    -----------------
     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME = 'PS',                                 &
        LONG_NAME  = 'surface_pressure',                   &
        UNITS      = 'Pa',                                 &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationNone,                   &
        RESTART    = MAPL_RestartSkip,   __RC__)

!    SHFX (pos is up) - why not evap, Ri, ???
!    ----------------------------------------
     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME = 'SH',                                 &
        LONG_NAME  = 'sensible_heat_flux_from_turbulence', &
        UNITS      = 'W m-2',                              &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationNone,                   &
        RESTART    = MAPL_RestartSkip,   __RC__)

!    TA -- Surface Air Temperature
!    ----
     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME = 'TA',                                 &
        LONG_NAME  = 'surface_temperature_from_surface',   &
        UNITS      = 'K',                                  &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationNone,                   &
        RESTART    = MAPL_RestartSkip,   __RC__)

!    DZ -- Surface Mid-layer Height
!    ----
     call MAPL_AddImportSpec(GC,                           &
        LONG_NAME  = 'surface_layer_height',               &
        UNITS      = 'm',                                  &
        SHORT_NAME = 'DZ',                                 &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationNone,                   &
        RESTART    = MAPL_RestartSkip,   __RC__)

!    TSOIL1, from SURFACE
!    --------------------
     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME = 'TSOIL1',                             &
        LONG_NAME  = 'soil_temperatures_layer_1',          &
        UNITS      = 'K',                                  &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationNone,                   &
        RESTART    = MAPL_RestartSkip,   __RC__)

!    U10M
!    ----
     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME = 'U10M',                               &
        LONG_NAME  = '10-meter_eastward_wind',             &
        UNITS      = 'm s-1',                              &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationNone,                   &
        RESTART    = MAPL_RestartSkip,   __RC__)

!    V10M
!    ----
     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME = 'V10M',                               &
        LONG_NAME  = '10-meter_northward_wind',            &
        UNITS      = 'm s-1',                              &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationNone,                   &
        RESTART    = MAPL_RestartSkip,   __RC__)

!    U10N
!    ----
     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME = 'U10N',                               &
        LONG_NAME  = 'equivalent_neutral_10-meter_eastward_wind', &
        UNITS      = 'm s-1',                              &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationNone,                   &
        RESTART    = MAPL_RestartSkip,   __RC__)

!    V10N
!    ----
     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME = 'V10N',                               &
        LONG_NAME  = 'equivalent_neutral_10-meter_northward_wind', &
        UNITS      = 'm s-1',                              &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationNone,                   &
        RESTART    = MAPL_RestartSkip,   __RC__)

!    USTAR
!    -----
     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME = 'USTAR',                              &
        LONG_NAME  = 'surface_velocity_scale',             &
        UNITS      = 'm s-1',                              &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationNone,                   &
        RESTART    = MAPL_RestartSkip,   __RC__)

!    Z0H
!    ---
     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME = 'Z0H',                                &
        LONG_NAME  = 'surface_roughness_for_heat',         &
        UNITS      = 'm',                                  &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationNone,                   &
        RESTART    = MAPL_RestartSkip,   __RC__)

!    Cell area
!    ---------
     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME = 'AREA',                               &
        LONG_NAME  = 'agrid_cell_area',                    &
        UNITS      = 'm^2',                                &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationNone,                   &
        RESTART    = MAPL_RestartSkip,   __RC__)

     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME = 'TS',                                 &
        LONG_NAME  = 'surface skin temperature',           &
        UNITS      = 'K',                                  &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationNone,                   &
        RESTART    = MAPL_RestartSkip,   __RC__)

     call MAPL_AddImportSpec(GC,                           & 
         SHORT_NAME = 'CNV_MFD',                           & 
         LONG_NAME  = 'detraining_mass_flux',              &
         UNITS      = 'kg m-2 s-1',                        &    
         DIMS       = MAPL_DimsHorzVert,                   &  
         VLOCATION  = MAPL_VLocationCenter,                &
         RESTART    = MAPL_RestartSkip,  __RC__)

     call MAPL_AddImportSpec(GC,                           &                  
         SHORT_NAME = 'CNV_MFC',                           & 
         LONG_NAME  = 'cumulative_mass_flux',              &
         UNITS      = 'kg m-2 s-1',                        &
         DIMS       = MAPL_DimsHorzVert,                   &
         VLOCATION  = MAPL_VLocationEdge,                  &
         RESTART    = MAPL_RestartSkip,  __RC__)

     call MAPL_AddImportSpec(GC,                           &
         SHORT_NAME = 'QLCN',                              &
         LONG_NAME  = 'mass_fraction_of_convective_cloud_liquid_water', &
         UNITS      = 'kg kg-1',                           &
         DIMS       = MAPL_DimsHorzVert,                   &
         VLOCATION  = MAPL_VLocationCenter,                &
         RESTART    = MAPL_RestartSkip,  __RC__)

     call MAPL_AddImportSpec(GC,                           &
         SHORT_NAME = 'QICN',                              &
         LONG_NAME  = 'mass_fraction_of_convective_cloud_ice_water', &
         UNITS      = 'kg kg-1',                           &
         DIMS       = MAPL_DimsHorzVert,                   &
         VLOCATION  = MAPL_VLocationCenter,                &
         RESTART    = MAPL_RestartSkip,  __RC__)

     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME  = 'SWNDSRF',                           &
        LONG_NAME   = 'surface_net_downward_shortwave_flux', &
        UNITS       = 'W m^-2',                            &
        DIMS        = MAPL_DimsHorzOnly,                   &
        VLOCATION   = MAPL_VLocationNone,                  &
        RESTART     = MAPL_RestartSkip,  __RC__)
     
     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME = 'CLDTT',                              &
        LONG_NAME  = 'total_cloud_area_fraction',          &
        UNITS      = '1',                                  &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationNone,                   &
        RESTART    = MAPL_RestartSkip,   __RC__)


     call ESMF_ConfigGetAttribute(CF, DO_CO2CNNEE, label='USE_CNNEE:', default=0, __RC__)

     if (DO_CO2CNNEE == 1) then
        call MAPL_AddImportSpec(GC,                           &
             SHORT_NAME = 'CNNEE',                              &
             LONG_NAME  = 'CN_net_ecosystem_exchange',          &
             UNITS      = 'kg m-2 s-1',                         &
             DIMS       = MAPL_DimsHorzOnly,                    &
             VLOCATION  = MAPL_VLocationNone, __RC__)
     endif

     
if ( r%doing_GOCART ) then

! !INTERNAL STATE:

!
!  NOTES: 
!    1)  vtitle as it stands is as the CF definition of long name.
!        I may need to create a "standard name" in chemReg and pass
!        this to GEOS Generic
!    2)  Host model MUST provide convective transport as well
!

    nq = r%nq     ! total number of chemical tracers
    
! Get BOOTSTRAP Default Values for GOCART INTERNAL
! ------------------------------------------------
    CALL ESMF_ConfigGetAttribute(CF, DEFVAL_CO2, Default=380.0e-6, Label='DEFVAL_CO2:', __RC__)

! Is GOCART providing O3 to the ANALYSIS bundle?
! ----------------------------------------------
    CALL ESMF_ConfigGetAttribute(CF, providerName, Default="PCHEM", Label='ANALYSIS_OX_PROVIDER:', __RC__)

! r%doing_O3 must be TRUE if the ANALYSIS_OX_PROVIDER is GOCART.
! --------------------------------------------------------------
    IF (providerName == 'GOCART' .AND. .NOT. r%doing_O3) THEN
     IF(MAPL_AM_I_ROOT()) THEN
      PRINT *," "
      PRINT *,TRIM(Iam)//": Set doing_O3 to yes in Chem_Registry.rc if GOCART "
      PRINT *,"             is the ANALYSIS_OX_PROVIDER."
      PRINT *," "
     END IF
     STATUS = 1
     VERIFY_(STATUS)
    END IF
     
!   Loop over all constituents on registry
!   --------------------------------------
    do n = r%i_GOCART, r%j_GOCART 

       if (state%data_driven) then 
          FRIENDLIES = trim(COMP_NAME)

          call ESMF_ConfigGetAttribute(CF, AEROFRIENDLY, Label='AERO_FRIENDLIES:', default=FRIENDLIES, __RC__)
       else
          if (trim(r%vname(n)) == 'OX' .and. trim(providerName) == 'GOCART') then
             FRIENDLIES = 'ANALYSIS:DYNAMICS:TURBULENCE:MOIST'
          else
             FRIENDLIES = 'DYNAMICS:TURBULENCE:MOIST'
          end if

!         Set aerosol friendly attribute to MOIST as function of Convective Parameterization
!         ----------------------------------------------------------------------------------

          short_name = ESMF_UtilStringUpperCase(trim(r%vname(n)))

       end if ! data or computational GC

                                         DEFVAL = 0.0
       if ( short_name(1:3) .eq. 'CO2' ) DEFVAL = DEFVAL_CO2

!      Aerosol Tracers to be transported
!      ---------------------------------

       call MAPL_AddInternalSpec(GC,               &
          SHORT_NAME = trim(COMP_NAME)//'::'//trim(r%vname(n)), &
          LONG_NAME  = r%vtitle(n),                &
          UNITS      = r%vunits(n),                &     
          FRIENDLYTO = FRIENDLIES,                 &
          RESTART    = MAPL_RestartOptional,       &
          DEFAULT    = DEFVAL,                     &
          DIMS       = MAPL_DimsHorzVert,          &
          VLOCATION  = MAPL_VLocationCenter, __RC__)

    end do


!   Call Legacy Set Services
!   ----------------------
    if (.not.state%data_driven) then
       call Aero_GridCompSetServices ( gc, r, __RC__)
    end if


!!EOS

!   Set the Profiling timers
!   ------------------------
    call MAPL_TimerAdd ( GC, name = 'RUN',        __RC__ )
    call MAPL_TimerAdd ( GC, name = 'INITIALIZE', __RC__ )
    call MAPL_TimerAdd ( GC, name = 'FINALIZE',   __RC__ )
    call MAPL_TimerAdd ( GC, name = 'AERO1',   __RC__ )
    call MAPL_TimerAdd ( GC, name = 'AERO2',   __RC__ )

    call MAPL_TimerAdd (GC, name = 'SS', __RC__)
    call MAPL_TimerAdd (GC, name = 'O3', __RC__)
    call MAPL_TimerAdd (GC, name = 'DU', __RC__)
    call MAPL_TimerAdd (GC, name = 'BC', __RC__)
    call MAPL_TimerAdd (GC, name = 'OC', __RC__)
    call MAPL_TimerAdd (GC, name = 'SU', __RC__)
    call MAPL_TimerAdd (GC, name = 'CO', __RC__)
    call MAPL_TimerAdd (GC, name = 'CO2', __RC__)
    call MAPL_TimerAdd (GC, name = 'NI', __RC__)
    call MAPL_TimerAdd (GC, name = 'BRC', __RC__)
    call MAPL_TimerAdd (GC, name = 'CH4', __RC__)
    call MAPL_TimerAdd (GC, name = 'CFC', __RC__)

end if ! doing GOCART


!   Generic Set Services
!   --------------------
    call MAPL_GenericSetServices ( GC, __RC__ )

!   All done
!   --------

    RETURN_(ESMF_SUCCESS)
  
  end subroutine SetServices


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Initialize_ --- Initialize Aero_GridComp (ESMF)
!
! !INTERFACE:
!

   subroutine Initialize_ ( gc, impChem, expChem, clock, rc )

! !USES:

   implicit NONE

! !INPUT PARAMETERS:

   type(ESMF_Clock),  intent(inout) :: clock      ! The clock

! !OUTPUT PARAMETERS:

   type(ESMF_GridComp), intent(inout)  :: gc      ! Grid Component
   type(ESMF_State), intent(inout) :: impChem     ! Import State
   type(ESMF_State), intent(inout) :: expChem     ! Export State
   integer, intent(out) ::  rc                    ! Error return code:
                                                  !  0 - all is well
                                                  !  1 - 

! !DESCRIPTION: This is a simple ESMF wrapper.
!
! !REVISION HISTORY:
!
!  27Feb2005 da Silva  First crack.
!
!EOP
!-------------------------------------------------------------------------

!  ErrLog Variables
!  ----------------
   character(len=ESMF_MAXSTR)      :: IAm = 'Initialize_'
   integer                         :: STATUS
   character(len=ESMF_MAXSTR)      :: COMP_NAME

   type(Chem_Registry), pointer    :: chemReg
   type(Aero_GridComp), pointer    :: gcChem      ! Grid Component
   type(Chem_Bundle), pointer      :: w_c         ! Chemical tracer fields     
   integer                         :: nymd, nhms  ! time
   real                            :: cdt         ! chemistry timestep (secs)
   integer                         :: hdt         ! model     timestep (secs)
   integer                         :: rft

   type(ESMF_Grid)                 :: grid       
 
   integer                         :: i1=1, i2, ig=0, im  ! dist grid indices
   integer                         :: j1=1, j2, jg=0, jm  ! dist grid indices
   integer                         :: km, nq              ! dist grid indices
   integer                         :: n, dims(3), l

   type(ESMF_Config)               :: CF
   character(len=ESMF_MAXSTR)      :: diurnal_bb

   type(MAPL_MetaComp), pointer    :: ggState      ! GEOS Generic State
   type(GOCART_state),  pointer    :: myState       ! GOCART state
   type(ESMF_State)                :: internal
   type(ESMF_Field)                :: field
   type(ESMF_Field)                :: fld
   type(ESMF_FieldBundle)          :: bundle
   type(ESMF_State)                :: aero
   type(ESMF_FieldBundle)          :: aero_state_aerosols
   type(ESMF_State)                :: aero_aci
   type(ESMF_FieldBundle)          :: aero_aci_aerosols
   character(len=ESMF_MAXSTR)      :: fld_name
   integer                         :: n_aerosols
   integer                         :: n_modes
   integer, parameter              :: n_gocart_modes = 13
   character(len=ESMF_MAXSTR)      :: aero_aci_modes(n_gocart_modes)
   character(len=ESMF_MAXSTR)      :: short_name
   real                            :: f_aci_seasalt, maxclean, ccntuning
   character(LEN=ESMF_MAXSTR)      :: CLDMICRO
     
   type(MAPL_VarSpec), pointer     :: InternalSpec(:)
   integer                         :: instance

   real(ESMF_KIND_R4), pointer, dimension(:,:) :: LATS
   real(ESMF_KIND_R4), pointer, dimension(:,:) :: LONS
   real(ESMF_KIND_R4), pointer, dimension(:,:) :: CELL_AREA

   type(ESMF_Calendar)     :: calendar
   type(ESMF_Time)         :: currentTime
   type(ESMF_Alarm)        :: alarm_H2O2
   type(ESMF_Alarm)        :: alarm_HNO3
   type(ESMF_Time)         :: ringTime
   type(ESMF_TimeInterval) :: ringInterval
   integer                 :: year, month, day, hh, mm, ss

   type(ESMF_State)           :: providerState
   character(len=ESMF_MAXSTR) :: prefix
  
  real(ESMF_KIND_R4),  dimension(4) :: Vect_Hcts

!  Get my name and set-up traceback handle
!  ---------------------------------------
   call ESMF_GridCompGet( GC, NAME=COMP_NAME, CONFIG=CF, __RC__ )

   Iam = trim(COMP_NAME) // '::' // 'Initialize_'

   if (MAPL_AM_I_ROOT()) then
      print *, TRIM(Iam)//': Starting...'
      print *,' '
   end if

!  Get my internal MAPL_Generic state
!  -----------------------------------
   call MAPL_GetObjectFromGC ( GC, ggState, __RC__ )

   call MAPL_TimerOn(ggState, 'INITIALIZE')

!  Compute proper REFERENCE_TIME if not explicitly specified
!  ---------------------------------------------------------
   call MAPL_GetResource( ggState, hdt, Label='RUN_DT:',                                   __RC__ )
   call MAPL_GetResource( ggState, cdt, Label='GOCART_DT:',             default=real(hdt), __RC__ )
   call MAPL_GetResource( ggState, rft, Label='GOCART_REFERENCE_TIME:', default=-999,      __RC__ )

   if (rft == -999 .and. int(cdt) /= hdt) then
       hh  = int(  hdt/3600          )
       mm  = int( (hdt-hh*3600)/60   )
       ss  = int(  hdt-hh*3600-mm*60 )
       rft = hh*10000 + mm*100 + ss

       call MAPL_ConfigSetAttribute(cf, value=rft, Label="GOCART_REFERENCE_TIME:", __RC__ )
       if (MAPL_AM_I_ROOT()) write(*,"(21x,'Re-Setting GOCART_REFERENCE_TIME: ',i6.6)") rft
   endif


!  Initialize GEOS Generic
!  ------------------------
   call MAPL_GenericInitialize ( gc, impChem, expChem, clock, __RC__ )

   call MAPL_TimerOn(ggState, 'TOTAL')

!  Get pre-ESMF parameters from gc and clock
!  -----------------------------------------
   call extract_ ( gc, clock, chemReg, gcChem, w_c, nymd, nhms, cdt, STATUS, state=myState )
   VERIFY_(STATUS)


!  Create Chem Bundle
!  ------------------
   call ESMF_GridCompGet ( GC, grid=grid, __RC__ )

   call MAPL_GridGet ( grid, globalCellCountPerDim=DIMS, __RC__ )

   im = dims(1)
   jm = dims(2)
   nq = chemReg%nq

   call ESMF_GridGet(GRID, localDE = 0, &
                           staggerloc = ESMF_STAGGERLOC_CENTER, &
                           computationalCount = DIMS, __RC__)

!  Associate the Internal State fields with our legacy state 
!  ---------------------------------------------------------
   call MAPL_Get ( ggSTATE, INTERNALSPEC = InternalSpec, &
                            INTERNAL_ESMF_STATE = internal, &
                            LONS = LONS, &
                            LATS = LATS, __RC__ )

! A-Grid cell area
! ----------------
  if (myState%data_driven) then
      CELL_AREA => null()
  else
      call MAPL_GetPointer(impChem, NAME='AREA', ptr=CELL_AREA, __RC__)
  end if 

! Local sizes of three dimensions
!--------------------------------
   i2 = dims(1)
   j2 = dims(2)
   km = dims(3)

!  Initalize the legacy state but do not allocate memory for arrays
!  ----------------------------------------------------------------
   call Chem_BundleCreate_ ( chemReg, i1, i2, ig, im, j1, j2, jg, jm, km,  &
                             w_c, lon=LONS(:,:), lat=LATS(:,:), cell_area=CELL_AREA, &
                             skipAlloc=.true., __RC__ )

   w_c%grid_esmf = grid  ! Will need this for I/O later

!   Check whether to de-activate diurnal biomass burning (default is *on*)
!   ----------------------------------------------------------------------
    call ESMF_ConfigGetAttribute(CF, diurnal_bb, label='DIURNAL_BIOMASS_BURNING:', default='yes', __RC__)

    if ( diurnal_bb(1:3) .eq. 'yes' .or. &
         diurnal_bb(1:3) .eq. 'YES' .or. &
         diurnal_bb(1:3) .eq. 'Yes' ) then
         if (MAPL_AM_I_ROOT()) print *, trim(Iam)//': Diurnal Biomass Burning is ON'
         w_c%diurnal_bb = .true.
    else
         if (MAPL_AM_I_ROOT()) print *, trim(Iam)//': Diurnal Biomass Burning is OFF'
         w_c%diurnal_bb = .false.
    endif

!  Allocate these because they are not friendly tracers
!  ----------------------------------------------------
   allocate(w_c%delp(i1:i2,j1:j2,km), w_c%rh(i1:i2,j1:j2,km), __STAT__)

   _ASSERT( size(InternalSpec) == chemReg%n_GOCART, 'needs informative message' )

   do L = 1, size(InternalSpec)

      call MAPL_VarSpecGet(InternalSpec(L), SHORT_NAME=short_name, __RC__)

      N = chemReg%i_GOCART + L - 1
      call MAPL_GetPointer(internal, NAME=short_name, ptr=w_c%qa(N)%data3d, __RC__)

   end do

#ifdef PRINT_STATES

   if (MAPL_AM_I_ROOT()) then
       print *, trim(Iam)//': INTERNAL State during Initialize():' 
       call ESMF_StatePrint ( internal )
       print *, trim(Iam)//': IMPORT   State during Initialize():'
       call ESMF_StatePrint ( impChem  )
       print *, trim(Iam)//': EXPORT   State during Initialize():' 
       call ESMF_StatePrint ( expChem  )
    end if

#endif

!   Call Legacy Initialize
!   ----------------------
    call Aero_GridCompInitialize ( gcChem, w_c, gc, impChem, expChem, &
                                   nymd, nhms, cdt, myState%data_driven, STATUS )
    VERIFY_(STATUS)


!   Only at this point we have the scavenging coefficients filled,
!    so annotate the convection friendly internal state
!   Note: Move this to AddInternalSpec but first we need to have
!         the subcomponents as bonafide ESMF components
!-srf added Henrys law constants
!   --------------------------------------------------------------
    do n = ChemReg%i_GOCART, ChemReg%j_GOCART 
       call ESMF_StateGet(internal, trim(COMP_NAME)//'::'//trim(ChemReg%vname(n)), field, __RC__)
       call ESMF_AttributeSet(field, NAME='ScavengingFractionPerKm', VALUE=ChemReg%fscav(n), __RC__)
       Vect_Hcts(1:4)= ChemReg%hcts(1:4,n) 
  
       call ESMF_AttributeSet(field, 'SetofHenryLawCts', Vect_Hcts, __RC__)
    end do


    call MAPL_TimerOff(ggState, 'TOTAL')
    call MAPL_TimerOff(ggState, 'INITIALIZE')

    RETURN_(ESMF_SUCCESS)

CONTAINS

       subroutine AddFromExportToBundle_(STATE, BUNDLE, NAME, RC)
         type(ESMF_State)       :: STATE
         type(ESMF_FieldBundle) :: BUNDLE
         CHARACTER(LEN=*)       :: NAME
         integer, optional      :: RC
         type(ESMF_Field) :: FIELD
                   __Iam__('AddFromExportToBundle_')
         call ESMF_StateGet( STATE, NAME, FIELD, __RC__ )
         call MAPL_AllocateCoupling( FIELD, __RC__ )
         call MAPL_FieldBundleAdd ( BUNDLE, FIELD, __RC__ )
         RETURN_(ESMF_SUCCESS)
       end subroutine AddFromExportToBundle_

   end subroutine Initialize_


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Run1_ --- Runs Aero_GridComp (ESMF)
!
! !INTERFACE:
!

   subroutine Run1_ ( gc, impChem, expChem, clock, rc )

! !USES:

  implicit NONE

! !INPUT PARAMETERS:

   type(ESMF_Clock),  intent(inout) :: clock      ! The clock

! !OUTPUT PARAMETERS:

   type(ESMF_GridComp), intent(inout)  :: gc      ! Grid Component
   type(ESMF_State), intent(inout) :: impChem     ! Import State
   type(ESMF_State), intent(inout) :: expChem     ! Export State
   integer, intent(out) ::  rc                    ! Error return code:
                                                  !  0 - all is well
                                                  !  1 - 

! !DESCRIPTION: This is a simple ESMF wrapper.
!
! !REVISION HISTORY:
!
!  27Feb2005 da Silva  First crack.
!
!EOP
!-------------------------------------------------------------------------

!  ErrLog Variables
!  ----------------
   character(len=ESMF_MAXSTR)      :: IAm
   integer                         :: STATUS
   character(len=ESMF_MAXSTR)      :: COMP_NAME

   type(Chem_Registry), pointer    :: chemReg
   type(Aero_GridComp), pointer    :: gcChem      ! Grid Component
   type(Chem_Bundle), pointer      :: w_c         ! Chemical tracer fields     
   integer                         :: nymd, nhms  ! time
   real                            :: cdt         ! chemistry timestep (secs)
   real                            :: hdt         ! heartbeat time step (secs)
   real, pointer                   :: var(:,:,:)
   integer                         :: n

   type(ESMF_Config)               :: CF

   type(MAPL_MetaComp), pointer    :: ggState      ! GEOS Generic State
   type(ESMF_Alarm)                :: ALARM

   real(ESMF_KIND_R4), pointer, dimension(:,:) :: LATS
   real(ESMF_KIND_R4), pointer, dimension(:,:) :: LONS

   type (MAPL_SunOrbit)            :: ORBIT
   real, allocatable, target       :: ZTH(:,:)    ! can be R8
   real(ESMF_KIND_R4), allocatable :: r4ZTH(:,:)
   real(ESMF_KIND_R4), allocatable :: SLR(:,:)

   real, pointer                   :: rh2(:,:,:)
   integer                         :: in, jn

   type(GOCART_state), pointer     :: myState


!                               ---

!  Get my name and set-up traceback handle
!  ---------------------------------------
   call ESMF_GridCompGet( GC, NAME=COMP_NAME, CONFIG=CF, __RC__ )
   Iam = trim(COMP_NAME) // '::' // 'Run1_'

!  Get my internal MAPL_Generic state
!  -----------------------------------
   call MAPL_GetObjectFromGC ( GC, ggState, __RC__)

   call MAPL_TimerOn(ggState, 'TOTAL')
   call MAPL_TimerOn(ggState, 'RUN')

!  Get parameters from generic state.
!  ----------------------------------
   call MAPL_Get(ggState, LONS=LONS, LATS=LATS, ORBIT=ORBIT, RUNALARM=ALARM, __RC__)

!  Get heartbeat time step 
!  -----------------------
   call MAPL_GetResource(ggState, hdt, label='RUN_DT:', __RC__)
 
!  Get pre-ESMF parameters from gc and clock
!  -----------------------------------------
   call extract_ ( gc, clock, chemReg, gcChem, w_c, nymd, nhms, cdt, STATUS, state=myState )
   VERIFY_(STATUS)

!  Until all the gas phase species handle GOCART_DT correctly, we must run at the heartbeat
!  ----------------------------------------------------------------------------------------
!  Assume that DT is always an integral number of seconds
!  Add a fraction to both (and then truncate to int), to avoid cases like 900 /= 899.999999
   _ASSERT(abs(cdt-hdt) < 0.1, 'Implementation of GOCART_DT is problematic; set GOCART_DT = HEARTBEAT_DT')

   allocate(r4ZTH(SIZE(LATS,1), SIZE(LATS,2)), __STAT__)
   allocate(  ZTH(SIZE(LATS,1), SIZE(LATS,2)), __STAT__)
   allocate(  SLR(SIZE(LATS,1), SIZE(LATS,2)), __STAT__)

!  Update solar zenith angle
!  --------------------------
   call MAPL_SunGetInsolation(LONS, LATS, ORBIT, r4ZTH, SLR, CLOCK=CLOCK, __RC__)

!  Set pointers for sine/cosine zenith angle
!  -----------------------------------------

!  w_c%sinz => ...
   ZTH = r4ZTH
   w_c%cosz => zth

!  Fill in delp
!  ------------
   call MAPL_GetPointer ( impChem, var, 'DELP', __RC__ )
   w_c%delp = var(:,:,:)

!  Fill in RH
!  ----------
   call MAPL_GetPointer ( impChem, rh2, 'RH2', __RC__ )
   w_c%rh = rh2

!  Make sure tracers remain positive
!  ---------------------------------
   in = size(w_c%delp,1);   jn = size(w_c%delp,2)
   do n = ChemReg%i_GOCART, ChemReg%j_GOCART 
      call Chem_UtilNegFiller ( w_c%qa(n)%data3d, w_c%delp, in, jn, &
                                qmin=tiny(1.0) )
   end do

!  Call pre-ESMF version: runs at the heartbeat
!  --------------------------------------------
   call MAPL_TimerOn(ggState,'AERO1')
   call Aero_GridCompRun1 ( gcChem, w_c, gc, impChem, expChem, &
                            nymd, nhms, hdt, STATUS )
   VERIFY_(STATUS)
   call MAPL_TimerOff(ggState,'AERO1')

   deallocate(SLR,   __STAT__)
   deallocate(ZTH,   __STAT__)
   deallocate(r4ZTH, __STAT__)

   call MAPL_TimerOff(ggState, 'RUN')
   call MAPL_TimerOff(ggState, 'TOTAL')

   RETURN_(ESMF_SUCCESS)

   end subroutine Run1_


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Run_ --- Runs Aero_GridComp (ESMF)
!
! !INTERFACE:
!

   subroutine Run2_ ( gc, impChem, expChem, clock, rc )

! !USES:

  implicit NONE

! !INPUT PARAMETERS:

   type(ESMF_Clock),  intent(inout) :: clock      ! The clock

! !OUTPUT PARAMETERS:

   type(ESMF_GridComp), intent(inout)  :: gc      ! Grid Component
   type(ESMF_State), intent(inout) :: impChem     ! Import State
   type(ESMF_State), intent(inout) :: expChem     ! Export State
   integer, intent(out) ::  rc                    ! Error return code:
                                                  !  0 - all is well
                                                  !  1 - 

! !DESCRIPTION: This is a simple ESMF wrapper.
!
! !REVISION HISTORY:
!
!  27Feb2005 da Silva  First crack.
!
!EOP
!-------------------------------------------------------------------------

!  ErrLog Variables
!  ----------------
   character(len=ESMF_MAXSTR)      :: IAm
   integer                         :: STATUS
   character(len=ESMF_MAXSTR)      :: COMP_NAME

   type(Chem_Registry), pointer    :: chemReg
   type(Aero_GridComp), pointer    :: gcChem      ! Grid Component
   type(Chem_Bundle), pointer      :: w_c         ! Chemical tracer fields     
   integer                         :: nymd, nhms  ! time
   real                            :: hdt         ! heartbeat timestep (secs)
   real                            :: cdt         ! chemistry timestep (secs)
   real, pointer                   :: var(:,:,:)
   integer                         :: n

   type(ESMF_Config)               :: CF

   type(MAPL_MetaComp), pointer    :: ggState      ! GEOS Generic State
   type(ESMF_Alarm)                :: ALARM

   real(ESMF_KIND_R4), pointer, dimension(:,:) :: LATS
   real(ESMF_KIND_R4), pointer, dimension(:,:) :: LONS

   type (MAPL_SunOrbit)            :: ORBIT
   real, allocatable, target       :: ZTH(:,:)    ! can be R8
   real(ESMF_KIND_R4), allocatable :: r4ZTH(:,:)
   real(ESMF_KIND_R4), allocatable :: SLR(:,:)

   real, pointer                   :: rh2(:,:,:)
   integer                         :: in, jn

   type(ESMF_State)                :: internal
   type(GOCART_state), pointer     :: myState
   real, pointer, dimension(:,:,:) :: ptr3d_int
   real, pointer, dimension(:,:,:) :: ptr3d_imp
   
   logical                         :: run_alarm
   logical                         :: alarm_is_ringing


!  Diagnostics
   real, pointer, dimension(:,:)   :: totexttau, totscatau, &
                                      totextt25, totscat25, &
                                      totexttfm, totscatfm, &
                                      totangstr 
   real, pointer, dimension(:,:)   :: pm,        pm25,      &
                                      pm_rh35,   pm25_rh35, &
                                      pm_rh50,   pm25_rh50
   real, pointer, dimension(:,:)   :: duexttau, duscatau, &
                                      duextt25, duscat25, &
                                      duexttfm, duscatfm, &
                                      duangstr, dusmass,  &
                                      dusmass25
   real, pointer, dimension(:,:)   :: ssexttau, ssscatau, &
                                      ssextt25, ssscat25, &
                                      ssexttfm, ssscatfm, &
                                      ssangstr, sssmass,  &
                                      sssmass25
   real, pointer, dimension(:,:)   :: niexttau, niscatau, &
                                      niextt25, niscat25, &
                                      niexttfm, niscatfm, &
                                      niangstr, nismass,  &
                                      nismass25
   real, pointer, dimension(:,:)   :: nh4smass
   real, pointer, dimension(:,:)   :: suexttau, suscatau, &
                                      suangstr, susmass
   real, pointer, dimension(:,:)   :: suexttauvolc, suscatauvolc, &
                                      suangstrvolc, susmassvolc
   real, pointer, dimension(:,:)   :: bcexttau, bcscatau, &
                                      bcangstr, bcsmass
   real, pointer, dimension(:,:)   :: ocexttau, ocscatau, &
                                      ocangstr, ocsmass
   real, pointer, dimension(:,:)   :: brcexttau, brcscatau, &
                                      brcangstr, brcsmass
   real, pointer, dimension(:,:,:) :: rh2x, delpx
   real, pointer, dimension(:,:,:) :: pso4, pso4v, pso4t
   real, allocatable               :: tau1(:,:), tau2(:,:)
   real                            :: c1, c2, c3

!                               ---

!  Get my name and set-up traceback handle
!  ---------------------------------------
   call ESMF_GridCompGet( GC, NAME=COMP_NAME, CONFIG=CF, __RC__ )
   Iam = trim(COMP_NAME) // '::' // 'Run2_'

!  Get my internal MAPL_Generic state
!  -----------------------------------
   call MAPL_GetObjectFromGC ( GC, ggState, __RC__)

   call MAPL_TimerOn(ggState, 'TOTAL')
   call MAPL_TimerOn(ggState, 'RUN')

!  Get heartbeat time step 
!  -----------------------
   call MAPL_GetResource(ggState, hdt, label='RUN_DT:', __RC__)

!  Is time to recycle H2O2 and HNO3?
!  ---------------------------------
   call extract_ ( gc, clock, chemReg, gcChem, w_c, nymd, nhms, cdt, STATUS, state=myState )
   VERIFY_(STATUS)

!  Get parameters from generic state.
!  ----------------------------------
   call MAPL_Get(ggState, LONS=LONS, LATS=LATS, ORBIT=ORBIT, RUNALARM=ALARM, __RC__)

!  Get pre-ESMF parameters from gc and clock
!  -----------------------------------------
   call extract_ ( gc, clock, chemReg, gcChem, w_c, nymd, nhms, cdt, STATUS, state=myState )
   VERIFY_(STATUS)

   if (myState%data_driven) then
       
       call MAPL_Get ( ggState, INTERNAL_ESMF_STATE=internal, __RC__ )

       do n = chemReg%i_GOCART, chemReg%j_GOCART
           call MAPL_GetPointer ( internal, NAME=trim(COMP_NAME)//'::'//trim(chemReg%vname(n)), ptr=ptr3d_int, __RC__ )
           call MAPL_GetPointer ( impChem,  NAME='clim'//trim(chemReg%vname(n)), ptr=ptr3d_imp, __RC__ )
             
           ptr3d_int = ptr3d_imp
       end do

       call MAPL_TimerOff(ggState, 'RUN')
       call MAPL_TimerOff(ggState, 'TOTAL')

       RETURN_(ESMF_SUCCESS)
   end if 


   allocate(r4ZTH(SIZE(LATS,1), SIZE(LATS,2)), __STAT__)
   allocate(  ZTH(SIZE(LATS,1), SIZE(LATS,2)), __STAT__)
   allocate(  SLR(SIZE(LATS,1), SIZE(LATS,2)), __STAT__)

!  Update solar zenith angle
!  --------------------------
   call MAPL_SunGetInsolation(LONS, LATS, ORBIT, r4ZTH, SLR, CLOCK=CLOCK, __RC__)

!  Set pointers for sine/cosine zenith angle
!  -----------------------------------------

!  w_c%sinz => ...
   ZTH = r4ZTH
   w_c%cosz => zth

!  Fill in delp
!  ------------
   call MAPL_GetPointer ( impChem, var, 'DELP', __RC__ )
   w_c%delp = var(:,:,:)

!  Fill in RH
!  ----------
   call MAPL_GetPointer ( impChem, rh2, 'RH2', __RC__ )
   w_c%rh = rh2

!  Make sure tracers remain positive
!  ---------------------------------
   in = size(w_c%delp,1);   jn = size(w_c%delp,2)
   do n = ChemReg%i_GOCART, ChemReg%j_GOCART 
      call Chem_UtilNegFiller ( w_c%qa(n)%data3d, w_c%delp, in, jn, &
                                qmin=tiny(1.0) )
   end do

!  Call pre-ESMF version
!  ---------------------
   run_alarm = ESMF_AlarmIsRinging(ALARM, RC=STATUS)

   call MAPL_TimerOn(ggState,'AERO2')
   call Aero_GridCompRun2 ( gcChem, w_c, gc, impChem, expChem, &
                            run_alarm, nymd, nhms, cdt, STATUS )
   VERIFY_(STATUS)
   call MAPL_TimerOff(ggState,'AERO2')

   if (run_alarm) then
      call ESMF_AlarmRingerOff(ALARM, __RC__)
   end if

   deallocate(SLR,   __STAT__)
   deallocate(ZTH,   __STAT__)
   deallocate(r4ZTH, __STAT__)

   call MAPL_TimerOff(ggState, 'RUN')
   call MAPL_TimerOff(ggState, 'TOTAL')

   RETURN_(ESMF_SUCCESS)

   end subroutine Run2_


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Finalize_ --- Finalize Aero_GridComp (ESMF)
!
! !INTERFACE:
!

   subroutine Finalize_ ( gc, impChem, expChem, clock, rc )

! !USES:

  implicit NONE

! !INPUT PARAMETERS:

   type(ESMF_Clock),  intent(inout) :: clock      ! The clock

! !OUTPUT PARAMETERS:

   type(ESMF_GridComp), intent(inout)  :: gc      ! Grid Component
   type(ESMF_State), intent(inout) :: impChem     ! Import State
   type(ESMF_State), intent(inout) :: expChem     ! Export State
   integer, intent(out) ::  rc                    ! Error return code:
                                                  !  0 - all is well
                                                  !  1 - 

! !DESCRIPTION: This is a simple ESMF wrapper.
!
! !REVISION HISTORY:
!
!  27Feb2005 da Silva  First crack.
!
!EOP
!-------------------------------------------------------------------------


!  ErrLog Variables
!  ----------------
   character(len=ESMF_MAXSTR)      :: IAm = 'Finalize_'
   integer                         :: STATUS
   character(len=ESMF_MAXSTR)      :: COMP_NAME

   type(Chem_Registry), pointer    :: chemReg
   type(Aero_GridComp), pointer    :: gcChem      ! Grid Component
   type(Chem_Bundle), pointer      :: w_c         ! Chemical tracer fields     
   integer                         :: nymd, nhms  ! time
   real                            :: cdt         ! chemistry timestep (secs)
   type(MAPL_MetaComp), pointer    :: ggState     ! GEOS Generic State
   type(GOCART_state), pointer     :: state

!  Get my name and set-up traceback handle
!  ---------------------------------------
   call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
   VERIFY_(STATUS)
   Iam = trim(COMP_NAME) // '::' // 'Finalize_'

!  Get my internal MAPL_Generic state
!  -----------------------------------
   call MAPL_GetObjectFromGC ( GC, ggState, __RC__)

   call MAPL_TimerON(ggState, 'TOTAL')
   call MAPL_TimerON(ggState, 'FINALIZE')

!  Get pre-ESMF parameters from gc and clock
!  -----------------------------------------
   call extract_ ( gc, clock, chemReg, gcChem, w_c, nymd, nhms, cdt, STATUS, &
                   state = state )
   VERIFY_(STATUS)

!  Call pre-ESMF version
!  ---------------------
   call Aero_GridCompFinalize ( gcChem, w_c, impChem, expChem, &
                                nymd, nhms, cdt, STATUS )
   VERIFY_(STATUS)

!  Destroy Chem_Bundle
!  -------------------
   call Chem_BundleDestroy ( w_c, STATUS )
   VERIFY_(STATUS)

!  Destroy Chem_Registry
!  ---------------------
   call Chem_RegistryDestroy ( chemReg, STATUS ) 
   VERIFY_(STATUS)

!  Destroy Legacy state
!  --------------------
   deallocate ( state%chemReg, state%gcChem, state%w_c, stat = STATUS )
   VERIFY_(STATUS)

   call MAPL_TimerOff(ggState, 'FINALIZE')
   call MAPL_TimerOff(ggState, 'TOTAL')

!  Finalize GEOS Generic
!  ---------------------
!ALT: do not deallocate "foreign objects"
   call MAPL_GenericFinalize ( gc, impChem, expChem, clock, __RC__ )

   RETURN_(ESMF_SUCCESS)

   end subroutine Finalize_


!.......................................................................

    subroutine extract_ ( gc, clock, chemReg, gcChem, w_c, nymd, nhms, cdt, &
                          rc, state )

    type(ESMF_GridComp), intent(INout)  :: gc
    type(ESMF_Clock), intent(in)     :: clock
    type(Chem_Registry), pointer     :: chemReg
    type(Aero_GridComp), pointer     :: gcChem
    type(Chem_Bundle), pointer       :: w_c
    integer, intent(out)             :: nymd, nhms
    real, intent(out)                :: cdt
    integer, intent(out)             :: rc
    type(MAPL_MetaComp), pointer            :: ggState
    type(GOCART_state), pointer, optional   :: state


    type(GOCART_state), pointer    :: myState

!   ErrLog Variables
!   ----------------
    character(len=ESMF_MAXSTR)      :: IAm
    integer                         :: STATUS
    character(len=ESMF_MAXSTR)      :: COMP_NAME

    type(ESMF_Alarm)                :: ALARM
    type(ESMF_TimeInterval)         :: RingInterval

    type(ESMF_Time)      :: TIME
    type(ESMF_Config)    :: CF
    type(GOCART_Wrap)    :: wrap
    integer              :: IYR, IMM, IDD, IHR, IMN, ISC
    real(ESMF_KIND_R8)   :: dt_r8

!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // '::' // 'extract_'

    rc = 0

!  Get my internal MAPL_Generic state
!  -----------------------------------
   call MAPL_GetObjectFromGC ( GC, ggState, __RC__ )


!   Get my internal state
!   ---------------------
    call ESMF_UserCompGetInternalState(gc, 'GOCART_state', WRAP, STATUS)
    VERIFY_(STATUS)
    myState => wrap%ptr
    if ( present(state) ) then
         state => wrap%ptr
    end if

!   This is likely to be allocated during initialize only
!   -----------------------------------------------------
    if ( .not. associated(myState%chemReg) ) then
         allocate ( myState%chemReg, stat=STATUS )
         VERIFY_(STATUS)
    end if
    if ( .not. associated(myState%gcChem) ) then
         allocate ( myState%gcChem, stat=STATUS )
         VERIFY_(STATUS)
    end if
    if ( .not. associated(myState%w_c) ) then
         allocate ( myState%w_c, stat=STATUS )
         VERIFY_(STATUS)
    end if

    chemReg => myState%chemReg
    gcChem  => myState%gcChem
    w_c     => myState%w_c

!   Get the configuration
!   ---------------------
    call ESMF_GridCompGet ( GC, CONFIG = CF, RC=STATUS )
    VERIFY_(STATUS)

!   Get time step
!   -------------
    call MAPL_Get(ggState, RUNALARM=ALARM, __RC__ )
    call ESMF_AlarmGet(ALARM, ringInterval=RingInterval, __RC__)

    call ESMF_TimeIntervalGet(RingInterval, s_r8=dt_r8, __RC__)
    cdt = real(dt_r8)


!   Need code to extract nymd(20050205), nhms(120000) from clock
!   ------------------------------------------

    call ESMF_ClockGet(CLOCK,currTIME=TIME,rc=STATUS)
    VERIFY_(STATUS)

    call ESMF_TimeGet(TIME ,YY=IYR, MM=IMM, DD=IDD, H=IHR, M=IMN, S=ISC, rc=STATUS)
    VERIFY_(STATUS)

    call MAPL_PackTime(NYMD,IYR,IMM,IDD)
    call MAPL_PackTime(NHMS,IHR,IMN,ISC)

    RETURN_(ESMF_SUCCESS)

   end subroutine extract_


logical function isDataDrivenGC_(gc, rc)

   implicit none

   type(ESMF_GridComp), intent(INout) :: gc
   integer, intent(out)               :: rc
 
!  local 
   character(len=ESMF_MAXSTR)         :: IAm
   integer                            :: STATUS

   integer                            :: i
   character(len=ESMF_MAXSTR)         :: comp_name
   character(len=*), parameter        :: modifier = '.data'

   call ESMF_GridCompGet(gc, name=comp_name, __RC__)   
   i = index(trim(comp_name), trim(modifier), back=.true.)
 
   if (i > 0) then 
       ! lets be strict
       if (comp_name(i:) == modifier) then
           isDataDrivenGC_ = .true.
       else
           isDataDrivenGC_ = .false.
       end if
   else
       isDataDrivenGC_ = .false.
   end if

   RETURN_(ESMF_SUCCESS)

end function isDataDrivenGC_


end module GOCART_GridCompMod

