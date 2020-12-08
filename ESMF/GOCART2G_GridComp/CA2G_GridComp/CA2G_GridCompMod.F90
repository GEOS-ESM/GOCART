#include "MAPL_Generic.h"

!=============================================================================
!BOP

! !MODULE: CA2G_GridCompMod - GOCART Carbonaceous Aerosol gridded component 

! !INTERFACE:
module CA2G_GridCompMod

!  !USES:
   use ESMF
   use MAPL
   use Chem_MieTableMod2G
   use Chem_AeroGeneric
   use iso_c_binding, only: c_loc, c_f_pointer, c_ptr

   use GOCART2G_Process       ! GOCART2G process library
   use GA_GridCompMod
   use MAPL_StringTemplate, only: StrTemplate

   implicit none
   private

   integer, parameter :: instanceComputational = 1
   integer, parameter :: instanceData          = 2

! !PUBLIC MEMBER FUNCTIONS:
   public  SetServices


! !DESCRIPTION: This module implements GOCART2G's Carbonaceous Aerosol (CA) Gridded Component.

! !REVISION HISTORY:
! 15June2020  Sherman, da Silva, Darmenov, Clune -  First attempt at refactoring.

!EOP
!===========================================================================

!  !Carbonaceous aerosol state
      type, extends(GA_GridComp) :: CA2G_GridComp
       integer            :: myDOW = -1   ! my Day of the week: Sun=1, Mon=2,...,Sat=7
       real               :: ratPOM = 1.0  ! Ratio of POM to OC mass
       real               :: fTerpene = 0.0 ! Fraction of terpene emissions -> aerosol
       real               :: fHydrophobic ! Initially hydrophobic portion
       logical            :: diurnal_bb   ! diurnal biomass burning
       real               :: aviation_layers(4)  ! heights of the LTO, CDS and CRS layers
!      !Workspae for point emissions
       logical                :: doing_point_emissions = .false.
       character(len=255)     :: point_emissions_srcfilen   ! filename for pointwise emissions
       integer                         :: nPts = -1
       integer, allocatable, dimension(:)  :: pstart, pend
       real, allocatable, dimension(:)     :: pLat, &
                                              pLon, &
                                              pBase, &
                                              pTop, &
                                              pEmis

   end type CA2G_GridComp 

   type wrap_
      type (CA2G_GridComp), pointer     :: PTR => null()
   end type wrap_

contains

!============================================================================
!BOP

! !IROUTINE: SetServices 

! !INTERFACE:
  subroutine SetServices (GC, RC)

!   !ARGUMENTS:
    type (ESMF_GridComp), intent(INOUT)   :: GC  ! gridded component
    integer,              intent(  OUT)   :: RC  ! return code

!   !DESCRIPTION: This version uses MAPL_GenericSetServices, which sets
!     the Initialize and Finalize services to generic versions. It also
!     allocates our instance of a generic state and puts it in the 
!     gridded component (GC). Here we only set the two-stage run method
!     and declare the data services.

!   !REVISION HISTORY: 
!   june2019   E.Sherman, A.Da Silva, A.Darmenov, T.Clune  First attempt at refactoring

!EOP
!============================================================================

!   !Locals
    character (len=ESMF_MAXSTR)              :: COMP_NAME
    type (ESMF_Config)                       :: cfg
    type (wrap_)                             :: wrap
    type (CA2G_GridComp), pointer            :: self

    character (len=ESMF_MAXSTR)              :: field_name, GCsuffix

    integer                                  :: i, nbins
    real                                     :: DEFVAL
    logical                                  :: data_driven = .true.

    __Iam__('SetServices')

!****************************************************************************
!   Begin...

!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet (GC, NAME=COMP_NAME, __RC__)
    Iam = trim(COMP_NAME) //'::'// Iam

    if (comp_name(1:5) == 'CA.oc') then
       GCsuffix = 'OC'
    else if (comp_name(1:5) == 'CA.bc') then
       GCsuffix = 'BC'
    else if (comp_name(1:5) == 'CA.br') then
       GCsuffix = 'BR'
    end if

!   Wrap internal state for storing in GC
!   -------------------------------------
    allocate (self, __STAT__)
    wrap%ptr => self

!   Load resource file  
!   -------------------
    cfg = ESMF_ConfigCreate (__RC__)
    call ESMF_ConfigLoadFile (cfg, 'CA2G_GridComp_'//trim(COMP_NAME)//'.rc', rc=status)
    if (status /= 0) then
        if (mapl_am_i_root()) print*,'CA2G_GridComp_'//trim(COMP_NAME)//'.rc does not exist! &
                                      Loading CA2G_GridComp_CA.oc.rc instead'
        call ESMF_ConfigLoadFile (cfg, 'CA2G_GridComp_CA.oc.rc', __RC__)
    end if

!   process generic config items
    call self%GA_GridComp%load_from_config( cfg, __RC__)

    call ESMF_ConfigGetAttribute (cfg, self%nbins, label='nbins:', __RC__)
    nbins = self%nbins

!   Parse config file into private internal state
!   ----------------------------------------------
    call ESMF_ConfigGetAttribute (cfg, self%myDOW, label='my_day_of_week:', default=-1, __RC__)
    call ESMF_ConfigGetAttribute (cfg, self%fhydrophobic, label='hydrophobic_fraction:', __RC__)
    call ESMF_ConfigGetAttribute (cfg, self%ratPOM, label='pom_oc_ratio:', default=1.0, __RC__)
    call ESMF_ConfigGetAttribute (cfg, self%fTerpene, label='terpene_emission_fraction:', default=0.0, __RC__)

    call ESMF_ConfigFindLabel (cfg, 'aviation_vertical_layers:', __RC__)
    do i=1,size(self%aviation_layers)
       call ESMF_ConfigGetAttribute (cfg, self%aviation_layers(i), __RC__)
    end do

    call ESMF_ConfigGetAttribute (cfg, self%point_emissions_srcfilen, &
                                  label='point_emissions_srcfilen:', default='/dev/null', __RC__)
    if ( (index(self%point_emissions_srcfilen,'/dev/null')>0) ) then
       self%doing_point_emissions = .false. ! disable it if no file specified
    else
       self%doing_point_emissions = .true.  ! we are good to go
    end if

!   Is CA data driven?
!   ------------------
    call determine_data_driven (COMP_NAME, data_driven, __RC__)

!   Set entry points
!   ------------------------
    call MAPL_GridCompSetEntryPoint (GC, ESMF_Method_Initialize,  Initialize, __RC__)
    call MAPL_GridCompSetEntryPoint (GC, ESMF_Method_Run, Run, __RC__)
    if (data_driven .neqv. .true.) then
       call MAPL_GridCompSetEntryPoint (GC, ESMF_Method_Run, Run2, __RC__)
    end if

    DEFVAL = 0.0

!   IMPORT STATE
!   -------------
    if (data_driven) then

!      Pressure at layer edges
!      -----------------------
       call MAPL_AddImportSpec(GC,                            &
          SHORT_NAME = 'PLE',                                 &
          LONG_NAME  = 'air_pressure',                        &
          UNITS      = 'Pa',                                  &
          DIMS       = MAPL_DimsHorzVert,                     &
          VLOCATION  = MAPL_VLocationEdge,                    &
          RESTART    = MAPL_RestartSkip,     __RC__)

!      RH: is between 0 and 1
!      ----------------------
       call MAPL_AddImportSpec(GC,                            &
          SHORT_NAME = 'RH2',                                 &
          LONG_NAME  = 'Rel_Hum_after_moist',                 &
          UNITS      = '1',                                   &
          DIMS       = MAPL_DimsHorzVert,                     &
          VLOCATION  = MAPL_VLocationCenter,                  &
          RESTART    = MAPL_RestartSkip,     __RC__)

       call MAPL_AddInternalSpec(GC,                          &
          short_name = 'CAphobic'//trim(comp_name),           &
          long_name  = trim(GCsuffix)//' phobic Mixing Ratio',&
          units      = 'kg kg-1',                             &
          restart    = MAPL_RestartOptional,                  &
          dims       = MAPL_DimsHorzVert,                     &
          vlocation  = MAPL_VLocationCenter, __RC__) 

       call MAPL_AddInternalSpec(GC,                          &
          short_name = 'CAphilic'//trim(comp_name),           &
          long_name  = trim(GCsuffix)//' philic Mixing Ratio',&
          units      = 'kg kg-1',                             &
          restart    = MAPL_RestartOptional,                  &
          dims       = MAPL_DimsHorzVert,                     &
          vlocation  = MAPL_VLocationCenter, __RC__) 

       call MAPL_AddImportSpec(GC,                            &
          short_name = 'clim'//trim(GCsuffix)//'phobic',      &
          long_name  = trim(GCsuffix)//' phobic Mixing Ratio',&
          units      = 'kg kg-1',                             &
          restart    = MAPL_RestartOptional,                  &
          dims       = MAPL_DimsHorzVert,                     &
          vlocation  = MAPL_VLocationCenter, __RC__) 

       call MAPL_AddImportSpec(GC,                            &
          short_name = 'clim'//trim(GCsuffix)//'philic',      &
          long_name  = trim(GCsuffix)//' philic Mixing Ratio',&
          units      = 'kg kg-1',                             &
          restart    = MAPL_RestartOptional,                  &
          dims       = MAPL_DimsHorzVert,                     &
          vlocation  = MAPL_VLocationCenter, __RC__)

       do i = 1, self%nbins
         write (field_name, '(A, I0.3)') '', i
!        !dry deposition
          call MAPL_AddImportSpec(GC,                                       &
             short_name = 'clim'//trim(GCsuffix)//'DP'//trim(field_name),                      &
             long_name  = 'Organic Carbon Mixing Ratio (bin '//trim(field_name)//')',  &
             units      = 'kg kg-1',                                         &
             dims       = MAPL_DimsHorzOnly,                                 &
             vlocation  = MAPL_VLocationCenter,                              &
             restart    = MAPL_RestartSkip, __RC__)

!        !wet deposition    
          call MAPL_AddImportSpec(GC,                                       &
             short_name = 'clim'//trim(GCsuffix)//'WT'//trim(field_name),                     &
             long_name  = 'Organic Carbon Mixing Ratio (bin '//trim(field_name)//')', &
             units      = 'kg kg-1',                                        &
             dims       = MAPL_DimsHorzOnly,                                &
             vlocation  = MAPL_VLocationCenter,                             &
             restart    = MAPL_RestartSkip, __RC__)

!        !gravitational settling
          call MAPL_AddImportSpec(GC,                                       &
             short_name = 'clim'//trim(GCsuffix)//'SD'//trim(field_name),                     &
             long_name  = 'Organic Carbon Mixing Ratio (bin '//trim(field_name)//')', &
             units      = 'kg kg-1',                                        &
             dims       = MAPL_DimsHorzOnly,                                &
             vlocation  = MAPL_VLocationCenter,                             &
             restart    = MAPL_RestartSkip, __RC__)

!        !convective scavenging
          call MAPL_AddImportSpec(GC,                                       &
             short_name = 'clim'//trim(GCsuffix)//'SV'//trim(field_name),                     &
             long_name  = 'Organic Carbon Mixing Ratio (bin '//trim(field_name)//')', &
             units      = 'kg kg-1',                                        &
             dims       = MAPL_DimsHorzOnly,                                &
             vlocation  = MAPL_VLocationCenter,                             &
             restart    = MAPL_RestartSkip, __RC__)
!         end if
       end do
    end if ! (data_driven)


!   Computational Import and Export states
!   --------------------------------------
    if (.not. data_driven) then
#include "CA2G_Export___.h"
#include "CA2G_Import___.h"
#include "CA2G_Internal___.h"
    end if

!   This state holds fields needed by radiation
!   ---------------------------------------------
    call MAPL_AddExportSpec (GC,                             &
       short_name = trim(COMP_NAME)//'_AERO',                &
       long_name  = 'aerosols_from_'//trim(COMP_NAME),       &
       units      = 'kg kg-1',                               &
       dims       = MAPL_DimsHorzVert,                       &
       vlocation  = MAPL_VLocationCenter,                    &
       datatype   = MAPL_StateItem, __RC__)

!   This state is needed by MOIST - It will contain aerosols
!   ----------------------------------------------------------
    call MAPL_AddExportSpec (GC,                                                  &
       short_name = trim(COMP_NAME)//'_AERO_ACI',                                 &
       long_name  = 'aerosol_cloud_interaction_aerosols_from_'//trim(COMP_NAME),  &
       units      = 'kg kg-1',                                                    &
       dims       = MAPL_DimsHorzVert,                                            &
       vlocation  = MAPL_VLocationCenter,                                         &
       datatype   = MAPL_StateItem, __RC__)

!   This bundle is needed by surface for snow albedo modification
!   by aerosol settling and deposition
!   ~~~DEVELOPERS NOTE~~~ Change to StateItem when possible
!   ---------------------------------------------------------------
    call MAPL_AddExportSpec (GC,                                  &
       short_name = trim(COMP_NAME)//'_AERO_DP',                                    &
       long_name  = 'aerosol_deposition_from_'//trim(COMP_NAME),  &
       units      = 'kg m-2 s-1',                                 &
       dims       = MAPL_DimsHorzOnly,                            &
       datatype   = MAPL_BundleItem, __RC__)

!   Store internal state in GC
!   --------------------------
    call ESMF_UserCompSetInternalState ( GC, 'CA2G_GridComp', wrap, STATUS )
    VERIFY_(STATUS)

!   Set generic services
!   ----------------------------------
    call MAPL_GenericSetServices (GC, __RC__)

    RETURN_(ESMF_SUCCESS)

  end subroutine SetServices

!============================================================================
!BOP

! !IROUTINE: Initialize 

! !INTERFACE:
  subroutine Initialize (GC, import, export, clock, RC)

!   !ARGUMENTS:
    type (ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type (ESMF_State),    intent(inout) :: import ! Import state
    type (ESMF_State),    intent(inout) :: export ! Export state
    type (ESMF_Clock),    intent(inout) :: clock  ! The clock
    integer, optional,    intent(  out) :: RC     ! Error code

! !DESCRIPTION: This initializes CA's Grid Component. It primaryily fills 
!               GOCART's AERO states with its carbonaceous aerosol fields. 

! !REVISION HISTORY: 
! june2019   E.Sherman  First attempt at refactoring

!EOP
!============================================================================
!   !Locals
    character (len=ESMF_MAXSTR)          :: COMP_NAME
    type (MAPL_MetaComp), pointer        :: MAPL
    type (ESMF_Grid)                     :: grid
    type (ESMF_State)                    :: internal
    type (ESMF_State)                    :: aero, aero_aci
    type (ESMF_State)                    :: providerState
    type (ESMF_Config)                   :: cfg, cf
    type (ESMF_FieldBundle)              :: Bundle_DP
    type (wrap_)                         :: wrap
    type (CA2G_GridComp), pointer        :: self

    integer, allocatable                 :: mieTable_pointer(:)
    integer                              :: i, dims(3), km
    integer                              :: instance
    type (ESMF_Field)                    :: field, fld
    character (len=ESMF_MAXSTR)          :: prefix, GCsuffix, diurnal_bb, bin_index
    character (len=ESMF_MAXSTR),allocatable :: aerosol_names(:)
    real, pointer, dimension(:,:,:)      :: int_ptr
    real                                 :: CDT         ! chemistry timestep (secs)
    integer                              :: HDT         ! model     timestep (secs)
    real, pointer, dimension(:,:,:)      :: ple
    logical                              :: data_driven
    integer                              :: NUM_BANDS
    logical                              :: bands_are_present

    __Iam__('Initialize')

!****************************************************************************
!   Begin... 

!   Get the target components name and set-up traceback handle.
!   -----------------------------------------------------------
    call ESMF_GridCompGet (GC, grid=grid, name=COMP_NAME, config=cf, __RC__)
    Iam = trim(COMP_NAME) // '::' //trim(Iam)

    if (comp_name(1:5) == 'CA.oc') then
       GCsuffix = 'OC'
    else if (comp_name(1:5) == 'CA.bc') then
       GCsuffix = 'BC'
    else if (comp_name(1:5) == 'CA.br') then
       GCsuffix = 'BR'
    end if

!   Get my internal MAPL_Generic state
!   -----------------------------------
    call MAPL_GetObjectFromGC (GC, MAPL, __RC__)

!   Get my internal private state
!   -----------------------------
    call ESMF_UserCompGetInternalState(GC, 'CA2G_GridComp', wrap, STATUS)
    VERIFY_(STATUS)
    self => wrap%ptr

!   Get dimensions
!   ---------------
    call MAPL_GridGet (grid, globalCellCountPerDim=dims, __RC__ )
    km = dims(3)
    self%km = km

!   Get DTs
!   -------
    call MAPL_GetResource(mapl, HDT, Label='RUN_DT:', __RC__)
    call MAPL_GetResource(mapl, CDT, Label='GOCART_DT:', default=real(HDT), __RC__)
    self%CDT = CDT

!   Check whether to de-activate diurnal biomass burning (default is *on*)
!   ----------------------------------------------------------------------
    call ESMF_ConfigGetAttribute(cf, diurnal_bb, label='DIURNAL_BIOMASS_BURNING:', &
                                 default='YES', __RC__)
    diurnal_bb = ESMF_UtilStringUpperCase(diurnal_bb, __RC__)
    if (trim(diurnal_bb) == 'YES') then
       self%diurnal_bb = .true.
    else
       self%diurnal_bb = .false.
    end if

!   Load resource file  
!   -------------------
    cfg = ESMF_ConfigCreate (__RC__)
    call ESMF_ConfigLoadFile (cfg, 'CA2G_GridComp_'//trim(COMP_NAME)//'.rc', rc=status)
    if (status /= 0) then
        if (mapl_am_i_root()) print*,'ERROR: CA2G_GridComp_'//trim(COMP_NAME)//'.rc does not exist!' 
        return
    end if

!   Call Generic Initialize 
!   ------------------------
    call MAPL_GenericInitialize (GC, import, export, clock, __RC__)

!   Get parameters from generic state.
!   -----------------------------------
    call MAPL_Get (MAPL, INTERNAL_ESMF_STATE=internal, __RC__)

!   Is CA data driven?
!   ------------------
    call determine_data_driven (COMP_NAME, data_driven, __RC__)

!   If this is a data component, the data is provided in the import
!   state via ExtData instead of the actual GOCART children
!   ---------------------------------------------------------------- 
    if ( data_driven ) then
        providerState = import
        prefix = 'clim'
    else
        providerState = export
        prefix = ''
    end if

!   Fill AERO States with carbon fields
!   ------------------------------------
    call ESMF_StateGet (export, trim(COMP_NAME)//'_AERO'    , aero    , __RC__)
    call ESMF_StateGet (export, trim(COMP_NAME)//'_AERO_ACI', aero_aci, __RC__)
    call ESMF_StateGet (export, trim(COMP_NAME)//'_AERO_DP' , Bundle_DP, __RC__)

    call ESMF_StateGet (internal, 'CAphobic'//trim(comp_name), field, __RC__)
    call ESMF_AttributeSet (field, NAME='ScavengingFractionPerKm', VALUE=self%fscav(1), __RC__)
    fld = MAPL_FieldCreate (field, 'CAphobic'//trim(comp_name), __RC__)
    call MAPL_StateAdd (aero, fld, __RC__)
!    call MAPL_StateAdd (aero_aci, fld, __RC__)

    self%klid = 1 ! temporary value
!   Set internal CAphobic values to 0 where above klid
    call MAPL_GetPointer (internal, int_ptr, 'CAphobic'//trim(comp_name), __RC__)
    call setZeroKlid(self%km, self%klid, int_ptr)

    call ESMF_StateGet (internal, 'CAphilic'//trim(comp_name), field, __RC__)
    call ESMF_AttributeSet (field, NAME='ScavengingFractionPerKm', VALUE=self%fscav(2), __RC__)
    fld = MAPL_FieldCreate (field, 'CAphilic'//trim(comp_name), __RC__)
    call MAPL_StateAdd (aero, fld, __RC__)
    call MAPL_StateAdd (aero_aci, fld, __RC__)

    if (.not. data_driven) then
!      Set klid
       call MAPL_GetPointer(import, ple, 'PLE', __RC__)
       call findKlid (self%klid, self%plid, ple, __RC__)
!      Set internal CAphilic values to 0 where above klid
       call MAPL_GetPointer (internal, int_ptr, 'CAphilic'//trim(comp_name), __RC__)
       call setZeroKlid(self%km, self%klid, int_ptr)
    end if

    if (data_driven) then
       instance = instanceData

       do i = 1, self%nbins
          write (bin_index, '(A, I0.3)') '', i
!         Dry deposition
          call append_to_bundle(trim(GCsuffix)//'DP'//trim(bin_index), providerState, prefix, Bundle_DP, __RC__)

!         Wet deposition (Convective scavenging)
          call append_to_bundle(trim(GCsuffix)//'SV'//trim(bin_index), providerState, prefix, Bundle_DP, __RC__)

!         Wet deposition
          call append_to_bundle(trim(GCsuffix)//'WT'//trim(bin_index), providerState, prefix, Bundle_DP, __RC__)

!         Gravitational Settling
          call append_to_bundle(trim(GCsuffix)//'SD'//trim(bin_index), providerState, prefix, Bundle_DP, __RC__)
       end do
    else
       instance = instanceComputational

!      Dry deposition
       call append_to_bundle('CADP'//trim(comp_name), providerState, prefix, Bundle_DP, __RC__)

!      Wet deposition (Convective scavenging)
       call append_to_bundle('CASV'//trim(comp_name), providerState, prefix, Bundle_DP, __RC__)

!      Wet deposition
       call append_to_bundle('CAWT'//trim(comp_name), providerState, prefix, Bundle_DP, __RC__)

!      Gravitational Settling
       call append_to_bundle('CASD'//trim(comp_name), providerState, prefix, Bundle_DP, __RC__)
    end if

    self%instance = instance

!   Create Radiation Mie Table
!   --------------------------
    call MAPL_GetResource (MAPL, NUM_BANDS, 'NUM_BANDS:', __RC__)

!   Get file names for the optical tables
    call ESMF_ConfigGetAttribute (cfg, self%rad_MieTable(instance)%optics_file, &
                                  label="aerosol_radBands_optics_file:", __RC__ )

    allocate (self%rad_MieTable(instance)%channels(NUM_BANDS), __STAT__ )
    self%rad_MieTable(instance)%nch = NUM_BANDS

    call ESMF_ConfigFindLabel(cfg, label="BANDS:", isPresent=bands_are_present, __RC__)

    if (bands_are_present) then
       call ESMF_ConfigGetAttribute (cfg, self%rad_MieTable(instance)%channels, label= "BANDS:", &
                                    count=self%rad_MieTable(instance)%nch, __RC__)
    else
       do i = 1, NUM_BANDS
          self%rad_MieTable(instance)%channels(i) = i
       end do
    endif

    allocate (self%rad_MieTable(instance)%mie_aerosol, __STAT__)
    self%rad_MieTable(instance)%mie_aerosol = Chem_MieTableCreate (self%rad_MieTable(instance)%optics_file, __RC__)
    call Chem_MieTableRead (self%rad_MieTable(instance)%mie_aerosol, NUM_BANDS, self%rad_MieTable(instance)%channels, __RC__)

!   Create Diagnostics Mie Table
!   -----------------------------
!   Get file names for the optical tables
    call ESMF_ConfigGetAttribute (cfg, self%diag_MieTable(instance)%optics_file, &
                                  label="aerosol_monochromatic_optics_file:", __RC__ )
    call ESMF_ConfigGetAttribute (cfg, self%diag_MieTable(instance)%nch, label="n_channels:", __RC__)
    call ESMF_ConfigGetAttribute (cfg, self%diag_MieTable(instance)%nmom, label="n_moments:", default=0,  __RC__)
    allocate (self%diag_MieTable(instance)%channels(self%diag_MieTable(instance)%nch), __STAT__ )
    call ESMF_ConfigGetAttribute (cfg, self%diag_MieTable(instance)%channels, &
                                  label= "aerosol_monochromatic_optics_wavelength:", __RC__)

    allocate (self%diag_MieTable(instance)%mie_aerosol, __STAT__)
    self%diag_MieTable(instance)%mie_aerosol = Chem_MieTableCreate (self%diag_MieTable(instance)%optics_file, __RC__ )
    call Chem_MieTableRead (self%diag_MieTable(instance)%mie_aerosol, self%diag_MieTable(instance)%nch, &
                            self%diag_MieTable(instance)%channels, rc, nmom=self%diag_MieTable(instance)%nmom)

!   Finish creating AERO state
!   --------------------------
    ! Mie Table instance/index
    call ESMF_AttributeSet(aero, name='mie_table_instance', value=instance, __RC__)

    ! Add variables to CA instance's aero state. This is used in aerosol optics calculations
    call add_aero (aero, label='air_pressure_for_aerosol_optics',      label2='PLE', grid=grid, typekind=MAPL_R4, __RC__)
    call add_aero (aero, label='relative_humidity_for_aerosol_optics', label2='RH',  grid=grid, typekind=MAPL_R4,__RC__)
!   call ESMF_StateGet (import, 'PLE', field, __RC__)
!   call MAPL_StateAdd (aero, field, __RC__)
!   call ESMF_StateGet (import, 'RH2', field, __RC__)
!   call MAPL_StateAdd (aero, field, __RC__)

    call add_aero (aero, label='extinction_in_air_due_to_ambient_aerosol',    label2='EXT', grid=grid, typekind=MAPL_R8,__RC__)
    call add_aero (aero, label='single_scattering_albedo_of_ambient_aerosol', label2='SSA', grid=grid, typekind=MAPL_R8,__RC__)
    call add_aero (aero, label='asymmetry_parameter_of_ambient_aerosol',      label2='ASY', grid=grid, typekind=MAPL_R8,__RC__)

    call ESMF_AttributeSet(aero, name='band_for_aerosol_optics', value=0, __RC__)

    mieTable_pointer = transfer(c_loc(self), [1])
    call ESMF_AttributeSet(aero, name='mieTable_pointer', valueList=mieTable_pointer, itemCount=size(mieTable_pointer), __RC__)

    allocate(aerosol_names(self%nbins), __STAT__)
    aerosol_names(1) = 'CAphobic'//trim(comp_name)
    aerosol_names(2) = 'CAphilic'//trim(comp_name)
    call ESMF_AttributeSet(aero, name='internal_varaible_name', valueList=aerosol_names, &
                           itemCount=size(aerosol_names), __RC__)

    call ESMF_MethodAdd(AERO, label='aerosol_optics', userRoutine=aerosol_optics, __RC__)

    RETURN_(ESMF_SUCCESS)

  end subroutine Initialize

!============================================================================

!BOP
! !IROUTINE: Run 

! !INTERFACE:
  subroutine Run (GC, import, export, clock, rc)

!   !ARGUMENTS:
    type (ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type (ESMF_State),    intent(inout) :: import ! Import state
    type (ESMF_State),    intent(inout) :: export ! Export state
    type (ESMF_Clock),    intent(inout) :: clock  ! The clock
    integer, optional,    intent(  out) :: rc     ! Error code:

! !DESCRIPTION: Run method for the Sea Salt Grid Component. Determines whether to run
!               data or computational run method.

!EOP
!============================================================================
!   !Locals
    character (len=ESMF_MAXSTR)       :: COMP_NAME
    type (MAPL_MetaComp), pointer     :: MAPL
    type (ESMF_State)                 :: internal

    logical                           :: data_driven

    __Iam__('Run')

!*****************************************************************************
!   Begin... 

!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet (GC, NAME=COMP_NAME, __RC__)
    Iam = trim(COMP_NAME) //'::'// Iam

!   Get my internal MAPL_Generic state
!   -----------------------------------
    call MAPL_GetObjectFromGC (GC, MAPL, __RC__)
!   Get parameters from generic state.
!   -----------------------------------
    call MAPL_Get (MAPL, INTERNAL_ESMF_STATE=internal, __RC__)

!   Is SS data driven?
!   ------------------
    call determine_data_driven (COMP_NAME, data_driven, __RC__)

!   Update INTERNAL state variables with ExtData
!   ---------------------------------------------
    if (data_driven) then
       call Run_data (GC, import, export, internal, __RC__)
    else
       call Run1 (GC, import, export, clock, __RC__)
    end if

    RETURN_(ESMF_SUCCESS)

  end subroutine Run


!============================================================================
!BOP
! !IROUTINE: Run1 

! !INTERFACE:
  subroutine Run1 (GC, import, export, clock, RC)

!   !ARGUMENTS:
    type (ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type (ESMF_State),    intent(inout) :: import ! Import state
    type (ESMF_State),    intent(inout) :: export ! Export state
    type (ESMF_Clock),    intent(inout) :: clock  ! The clock
    integer, optional,    intent(  out) :: RC     ! Error code:

! !DESCRIPTION:  Computes emissions/sources for Sea Salt

!EOP
!============================================================================
!   !Locals
    character (len=ESMF_MAXSTR)       :: COMP_NAME
    type (MAPL_MetaComp), pointer     :: mapl
    type (ESMF_State)                 :: internal
    type (ESMF_Grid)                  :: grid
    type (wrap_)                      :: wrap
    type (CA2G_GridComp), pointer     :: self
    type(ESMF_Time)                   :: time

    character(len=3) :: cdow
    integer          :: idow
    integer          :: nymd, nhms, iyr, imm, idd, ihr, imn, isc
    real, pointer, dimension(:,:)     :: lats
    real, pointer, dimension(:,:)     :: lons
    real, dimension(:,:), allocatable :: biomass_src, terpene_src, biofuel_src, &
          eocant1_src, eocant2_src, oc_ship_src, aviation_lto_src, aviation_cds_src, &
          aviation_crs_src, biomass_src_
    real, dimension(:,:,:), allocatable :: emissions_point
    integer, pointer, dimension(:)    :: iPoint, jPoint
    character (len=ESMF_MAXSTR)  :: fname ! file name for point source emissions
    character(len=2)  :: GCsuffix

    real, pointer, dimension(:,:,:)  :: intPtr_phobic, intPtr_philic

#include "CA2G_DeclarePointer___.h"

   __Iam__('Run1')

!*****************************************************************************
!   Begin... 

!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet (GC, grid=grid, NAME=COMP_NAME, __RC__)
    Iam = trim(COMP_NAME) //'::'// Iam

!   Get my internal MAPL_Generic state
!   -----------------------------------
    call MAPL_GetObjectFromGC (GC, mapl, __RC__)

!   Get parameters from generic state.
!   -----------------------------------
    call MAPL_Get (mapl, INTERNAL_ESMF_STATE=internal, & 
                        LONS = LONS, &
                        LATS = LATS, __RC__ )

#include "CA2G_GetPointer___.h"

    if (comp_name(1:5) == 'CA.oc') then
       GCsuffix = 'OC'
    else if (comp_name(1:5) == 'CA.bc') then
       GCsuffix = 'BC'
    else if (comp_name(1:5) == 'CA.br') then
       GCsuffix = 'BR' 
    end if

    call MAPL_GetPointer (internal, intPtr_phobic, 'CAphobic'//trim(comp_name), __RC__)
    call MAPL_GetPointer (internal, intPtr_philic, 'CAphilic'//trim(comp_name), __RC__)

!   Get my private internal state
!   ------------------------------
    call ESMF_UserCompGetInternalState(GC, 'CA2G_GridComp', wrap, STATUS)
    VERIFY_(STATUS)
    self => wrap%ptr

    allocate(emissions_point, mold=delp,  __STAT__)
    emissions_point = 0.0

!   Extract nymd(yyyymmdd) from clock
!   ---------------------------------
    call ESMF_ClockGet (clock, currTime=time, __RC__)
    call ESMF_TimeGet (time ,YY=iyr, MM=imm, DD=idd, H=ihr, M=imn, S=isc, __RC__)
    call MAPL_PackTime (nymd, iyr, imm , idd)
    call MAPL_PackTime (nhms, ihr, imn, isc)

!   Reset tracer to zero at 0Z on specific day of week
!   --------------------------------------------------
    idow = Chem_UtilIdow(nymd)
    if ( (nhms==0) .and. (idow == self%myDOW) ) then
       cdow = Chem_UtilCdow(nymd)
       intPtr_phobic = tiny(1.) ! avoid division by zero
       intPtr_philic = tiny(1.) ! avoid division by zero
       if ( MAPL_AM_I_ROOT() ) then
          print *, '<> CA '//cdow//' tracer being set to zero on ', nymd, nhms
       end if
    end if

!   Implicit allocation with Fortran 2003
    if (trim(comp_name) == 'CA.oc') then
       biomass_src = OC_BIOMASS
       terpene_src = OC_TERPENE
       biofuel_src = OC_BIOFUEL
       eocant1_src = OC_ANTEOC1
       eocant2_src = OC_ANTEOC2
       oc_ship_src = OC_SHIP
       aviation_lto_src = OC_AVIATION_LTO
       aviation_cds_src = OC_AVIATION_CDS
       aviation_crs_src = OC_AVIATION_CRS
    else if (trim(comp_name) == 'CA.bc') then
       biomass_src = BC_BIOMASS
       biofuel_src = BC_BIOFUEL
       eocant1_src = BC_ANTEBC1
       eocant2_src = BC_ANTEBC2
       oc_ship_src = BC_SHIP
       aviation_lto_src = BC_AVIATION_LTO
       aviation_cds_src = BC_AVIATION_CDS
       aviation_crs_src = BC_AVIATION_CRS
       allocate(terpene_src, mold=BC_BIOMASS, __STAT__)
       terpene_src = 1.0
    else if (trim(comp_name) == 'CA.br') then
       biomass_src = BRC_BIOMASS
       terpene_src = BRC_TERPENE
       biofuel_src = BRC_BIOFUEL
       eocant1_src = BRC_ANTEBRC1
       eocant2_src = BRC_ANTEBRC2
       oc_ship_src = BRC_SHIP
       aviation_lto_src = BRC_AVIATION_LTO
       aviation_cds_src = BRC_AVIATION_CDS
       aviation_crs_src = BRC_AVIATION_CRS
    end if

!   As a safety check, where value is undefined set to 0
    where(1.01*biomass_src > MAPL_UNDEF) biomass_src = 0.
    where(1.01*terpene_src > MAPL_UNDEF) terpene_src = 0.
    where(1.01*biofuel_src > MAPL_UNDEF) biofuel_src = 0.
    where(1.01*eocant1_src > MAPL_UNDEF) eocant1_src = 0.
    where(1.01*eocant2_src > MAPL_UNDEF) eocant2_src = 0.
    where(1.01*oc_ship_src > MAPL_UNDEF) oc_ship_src = 0.
    where(1.01*aviation_lto_src > MAPL_UNDEF) aviation_lto_src = 0.
    where(1.01*aviation_cds_src > MAPL_UNDEF) aviation_cds_src = 0.
    where(1.01*aviation_crs_src > MAPL_UNDEF) aviation_crs_src = 0.

!   Save this in case we need to apply diurnal cycle
!   ------------------------------------------------
    if ( self%diurnal_bb ) then
       biomass_src_ = biomass_src
    end if

!   Apply diurnal cycle if so desired
!   ---------------------------------
    if ( self%diurnal_bb ) then
       call Chem_BiomassDiurnal (biomass_src, biomass_src_, &
                                 lons(:,:)*real(MAPL_RADIANS_TO_DEGREES), &
                                 lats(:,:)*real(MAPL_RADIANS_TO_DEGREES), &
                                 nhms, self%cdt)
    end if

    call CAEmission (self%diag_MieTable(self%instance), self%km, self%nbins, self%cdt, MAPL_GRAV, GCsuffix, self%ratPOM, &
                     self%fTerpene, aviation_lto_src, aviation_cds_src, aviation_crs_src, self%fHydrophobic, &
                     zpbl, t, airdens, rh2, intPtr_philic, intPtr_phobic, delp, self%aviation_layers, biomass_src, &
                     terpene_src, eocant1_src, eocant2_src, oc_ship_src, biofuel_src, &
                     CAEM, CAEMAN, CAEMBB, CAEMBF, CAEMBG, __RC__ )

!   Read any pointwise emissions, if requested
!   ------------------------------------------
    if(self%doing_point_emissions) then
       call StrTemplate(fname, self%point_emissions_srcfilen, xid='unknown', &
                         nymd=nymd, nhms=120000 )
       call ReadPointEmissions (nymd, fname, self%nPts, self%pLat, self%pLon, &
                                 self%pBase, self%pTop, self%pEmis, self%pStart, &
                                 self%pEnd, label='source')
    endif

!   Get indices for point emissions
!   -------------------------------
    if (self%nPts > 0) then
        allocate(iPoint(self%nPts), jPoint(self%nPts),  __STAT__)
        call MAPL_GetHorzIJIndex(self%nPts, iPoint, jPoint, &
                                 grid = grid,               &
                                 lon  = self%pLon/real(MAPL_RADIANS_TO_DEGREES), &
                                 lat  = self%pLat/real(MAPL_RADIANS_TO_DEGREES), &
                                 rc   = status)
            if ( status /= 0 ) then
                if (mapl_am_i_root()) print*, trim(Iam), ' - cannot get indices for point emissions'
                VERIFY_(status)
            end if

        call updatePointwiseEmissions (self%km, self%pBase, self%pTop, self%pEmis, self%nPts, &
                                       self%pStart, self%pEnd, zle, &
                                       area, iPoint, jPoint, nhms, emissions_point, __RC__)

       intPtr_phobic = intPtr_phobic + self%fHydrophobic * self%cdt * MAPL_GRAV / delp * emissions_point
       intPtr_philic = intPtr_philic + (1-self%fHydrophobic) * self%cdt * MAPL_GRAV / delp * emissions_point
    end if

    RETURN_(ESMF_SUCCESS)

  end subroutine Run1

!============================================================================
!BOP
! !IROUTINE: Run2 

! !INTERFACE:

  subroutine Run2 (GC, import, export, clock, RC)

    ! !ARGUMENTS:
    type (ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type (ESMF_State),    intent(inout) :: import ! Import state
    type (ESMF_State),    intent(inout) :: export ! Export state
    type (ESMF_Clock),    intent(inout) :: clock  ! The clock
    integer, optional,    intent(  out) :: RC     ! Error code:

! !DESCRIPTION: Run2 method for the CA Grid Component.

!EOP
!============================================================================
! Locals
    character (len=ESMF_MAXSTR)       :: COMP_NAME
    type (MAPL_MetaComp), pointer     :: MAPL
    type (ESMF_State)                 :: internal
    type (wrap_)                      :: wrap
    type (CA2G_GridComp), pointer     :: self
    type(MAPL_VarSpec), pointer       :: InternalSpec(:)

    integer                           :: n
    real, allocatable, dimension(:,:) :: drydepositionfrequency, dqa
    real                              :: fwet
    logical                           :: KIN
    real, allocatable, dimension(:,:,:)   :: pSOA_VOC
    real, pointer, dimension(:,:,:)       :: int_ptr
    real, allocatable, dimension(:,:,:,:) :: int_arr
    character(len=2)  :: GCsuffix
    character(len=ESMF_MAXSTR)      :: short_name
    real, pointer, dimension(:,:,:)  :: intPtr_phobic, intPtr_philic

    real, parameter ::  cpd    = 1004.16

#include "CA2G_DeclarePointer___.h"

    __Iam__('Run2')

!*****************************************************************************
!   Begin... 

!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet (GC, NAME=COMP_NAME, __RC__)
    Iam = trim(COMP_NAME) // '::' // Iam

!   Get my internal MAPL_Generic state
!   -----------------------------------
    call MAPL_GetObjectFromGC (GC, MAPL, __RC__)

!   Get parameters from generic state.
!   -----------------------------------
    call MAPL_Get (MAPL, INTERNAL_ESMF_STATE=internal, &
                         INTERNALSPEC = InternalSpec, __RC__)

    call MAPL_GetPointer (internal, intPtr_phobic, 'CAphobic'//trim(comp_name), __RC__)
    call MAPL_GetPointer (internal, intPtr_philic, 'CAphilic'//trim(comp_name), __RC__)

#include "CA2G_GetPointer___.h"

    if (comp_name(1:5) == 'CA.oc') then
       GCsuffix = 'OC'
    else if (comp_name(1:5) == 'CA.bc') then
       GCsuffix = 'BC'
    else if (comp_name(1:5) == 'CA.br') then
       GCsuffix = 'BR'
    end if

!   Get my private internal state
!   ------------------------------
    call ESMF_UserCompGetInternalState(GC, 'CA2G_GridComp', wrap, STATUS)
    VERIFY_(STATUS)
    self => wrap%ptr

!   Add on SOA from Anthropogenic VOC oxidation
!   -------------------------------------------
    if (trim(comp_name) == 'CA.oc') then
       pSOA_VOC = pSOA_ANTHRO_VOC
       where (1.01 * pSOA_VOC > MAPL_UNDEF) pSOA_VOC = 0.0

       intPtr_philic = intPtr_philic + self%cdt * pSOA_VOC/airdens
       if (associated(CAPSOA)) &
          CAPSOA = CAPSOA+sum(self%cdt*pSOA_VOC*delp/airdens/MAPL_GRAV, 3)
    end if

    if (trim(comp_name) == 'CA.br') then
       pSOA_VOC = pSOA_BIOB_VOC 
       where (1.01 * pSOA_VOC > MAPL_UNDEF) pSOA_VOC = 0.0

       intPtr_philic = intPtr_philic + self%cdt * pSOA_VOC/airdens
       if (associated(CAPSOA)) &
          CAPSOA = sum(self%cdt*pSOA_VOC*delp/airdens/MAPL_GRAV, 3)
    end if

!   Ad Hoc transfer of hydrophobic to hydrophilic aerosols
!   Following Chin's parameterization, the rate constant is
!   k = 4.63e-6 s-1 (.4 day-1; e-folding time = 2.5 days)
    call phobicTophilic (intPtr_phobic, intPtr_philic, CAHYPHIL, self%km, self%cdt, MAPL_GRAV, delp, __RC__)

!   CA Settling
!   -----------
    do n = 1, self%nbins
       call MAPL_VarSpecGet(InternalSpec(n), SHORT_NAME=short_name, __RC__)
       call MAPL_GetPointer(internal, NAME=short_name, ptr=int_ptr, __RC__)

       call Chem_Settling2Gorig (self%km, self%klid, self%rhFlag, n, int_ptr, MAPL_GRAV, delp, &
                                 self%radius(n)*1.e-6, self%rhop(n), self%cdt, t, airdens, &
                                 rh2, zle, CASD, __RC__)
    end do

!   CA Deposition
!   -----------
    allocate(dqa, mold=lwi, __STAT__)
    allocate(drydepositionfrequency, mold=lwi, __STAT__)

    drydepositionfrequency = 0.
    call DryDeposition(self%km, t, airdens, zle, lwi, ustar, zpbl, sh,&
                       MAPL_KARMAN, cpd, MAPL_GRAV, z0h, drydepositionfrequency, __RC__ )

    do n = 1, self%nbins
       call MAPL_VarSpecGet(InternalSpec(n), SHORT_NAME=short_name, __RC__)
       call MAPL_GetPointer(internal, NAME=short_name, ptr=int_ptr, __RC__)
       dqa = 0.
       dqa = max(0.0, int_ptr(:,:,self%km)*(1.-exp(-drydepositionfrequency*self%cdt)))
       int_ptr(:,:,self%km) = int_ptr(:,:,self%km) - dqa
       if (associated(CADP)) then
          CADP(:,:,n) = dqa * delp(:,:,self%km) / MAPL_GRAV / self%cdt
       end if
    end do

!   Large-scale Wet Removal
!   -------------------------------
!   Hydrophobic mode (first tracer) is not removed
    if (associated(CAWT)) CAWT(:,:,1)=0.0
    KIN = .true.
!   Hydrophilic mode (second tracer) is removed
    fwet = 1.
    call WetRemovalGOCART2G (self%km, self%klid, self%nbins, self%nbins, 2, self%cdt, GCsuffix, &
                             KIN, MAPL_GRAV, fwet, CAphilic, ple, t, airdens, &
                             pfl_lsan, pfi_lsan, cn_prcp, ncn_prcp, CAWT, rc)

!   Compute diagnostics
!   -------------------
    allocate(int_arr((ubound(intPtr_phobic,1)), (ubound(intPtr_phobic,2)), self%km, self%nbins), __STAT__)
    int_arr(:,:,:,1) = intPtr_phobic
    int_arr(:,:,:,2) = intPtr_philic

    call Aero_Compute_Diags (mie_table=self%diag_MieTable(self%instance), km=self%km, klid=self%klid, nbegin=1, nbins=2, &
                             channels=self%diag_MieTable(self%instance)%channels, aerosol=int_arr, grav=MAPL_GRAV, &
                             tmpu=t, rhoa=airdens, rh=rh2, u=u, v=v, delp=delp, sfcmass=CASMASS, colmass=CACMASS, &
                             mass=CAMASS, exttau=CAEXTTAU, scatau=CASCATAU, fluxu=CAFLUXU, fluxv=CAFLUXV, &
                             conc=CACONC, extcoef=CAEXTCOEF, scacoef=CASCACOEF, angstrom=CAANGSTR, aerindx=CAAERIDX,&
                             NO3nFlag=.false., __RC__)

if (trim(comp_name) == 'CA.oc') then
  if(mapl_am_i_root()) print*,'CA2G Run2 E sum(OCphobic) = ',sum(intPtr_phobic)
  if(mapl_am_i_root()) print*,'CA2G Run2 E sum(OCphilic) = ',sum(intPtr_philic)
!  if(mapl_am_i_root()) print*,'CA2G Run2 E array OCphobic = ',CAphobic
!  if(mapl_am_i_root()) print*,'CA2G Run2 E array OCphilic = ',CAphilic
!if(mapl_am_i_root()) print*,'CA2G OC CASMASS = ',sum(CASMASS)
!if(mapl_am_i_root()) print*,'CA2G OC CAMASS = ',sum(CAMASS)
!if(mapl_am_i_root()) print*,'CA2G OC CAEXTTAU = ',sum(CAEXTTAU)
!if(mapl_am_i_root()) print*,'CA2G OC CASCATAU = ',sum(CASCATAU)
else if (trim(comp_name) == 'CA.bc') then
  if(mapl_am_i_root()) print*,'CA2G Run2 E sum(BCphobic) = ',sum(intPtr_phobic)
  if(mapl_am_i_root()) print*,'CA2G Run2 E sum(BCphilic) = ',sum(intPtr_philic)
!  if(mapl_am_i_root()) print*,'CA2G Run2 E array BCphobic = ',CAphobic
!  if(mapl_am_i_root()) print*,'CA2G Run2 E array BCphilic = ',CAphilic
!if(mapl_am_i_root()) print*,'CA2G BC CASMASS = ',sum(CASMASS)
!if(mapl_am_i_root()) print*,'CA2G BC CAMASS = ',sum(CAMASS)
!if(mapl_am_i_root()) print*,'CA2G BC CAEXTTAU = ',sum(CAEXTTAU)
!if(mapl_am_i_root()) print*,'CA2G BC CASCATAU = ',sum(CASCATAU)
else if (trim(comp_name) == 'CA.br') then
  if(mapl_am_i_root()) print*,'CA2G Run2 E sum(BRphobic) = ',sum(intPtr_phobic)
  if(mapl_am_i_root()) print*,'CA2G Run2 E sum(BRphilic) = ',sum(intPtr_philic)
end if


    RETURN_(ESMF_SUCCESS)

  end subroutine Run2

!============================================================================
!BOP
! !IROUTINE: Run_data -- ExtData Sea Salt Grid Component

! !INTERFACE:

  subroutine Run_data (GC, IMPORT, EXPORT, INTERNAL, RC)

    ! !ARGUMENTS:

    type (ESMF_GridComp), intent(inout) :: GC       ! Gridded component 
    type (ESMF_State),    intent(inout) :: IMPORT   ! Import state
    type (ESMF_State),    intent(inout) :: EXPORT   ! Export state
    type (ESMF_State),    intent(inout) :: INTERNAL ! Interal state
    integer, optional,    intent(  out) :: RC       ! Error code:

! !DESCRIPTION: Updates pointers in Internal state with fields from ExtData. 

!EOP
!============================================================================
! Locals
    character (len=ESMF_MAXSTR)        :: COMP_NAME
    type (wrap_)                       :: wrap
    type (CA2G_GridComp), pointer      :: self

    character (len=ESMF_MAXSTR)        :: GCsuffix

    real, pointer, dimension(:,:,:)  :: ptr3d_int_phobic, ptr3d_int_philic
    real, pointer, dimension(:,:,:)  :: ptr3d_imp

    __Iam__('Run_data')

!*****************************************************************************
! Begin... 

! Get my name and set-up traceback handle
! ---------------------------------------
    call ESMF_GridCompGet (GC, NAME=COMP_NAME, __RC__)
    Iam = trim(COMP_NAME) //'::'//Iam

    if (comp_name(1:5) == 'CA.oc') then
       GCsuffix = 'OC'
    else if (comp_name(1:5) == 'CA.bc') then
       GCsuffix = 'BC'
    else if (comp_name(1:5) == 'CA.br') then
       GCsuffix = 'BR'
    end if

!   Get my private internal state
!   ------------------------------
    call ESMF_UserCompGetInternalState(GC, 'CA2G_GridComp', wrap, STATUS)
    VERIFY_(STATUS)
    self => wrap%ptr

!   Update interal data pointers with ExtData
!   -----------------------------------------
    call MAPL_GetPointer (internal, NAME='CAphobic'//trim(comp_name), ptr=ptr3d_int_phobic, __RC__)
    call MAPL_GetPointer (internal, NAME='CAphilic'//trim(comp_name), ptr=ptr3d_int_philic, __RC__)

    call MAPL_GetPointer (import, NAME='clim'//trim(GCsuffix)//'phobic', ptr=ptr3d_imp, __RC__)
    ptr3d_int_phobic = ptr3d_imp

    call MAPL_GetPointer (import, NAME='clim'//trim(GCsuffix)//'philic', ptr=ptr3d_imp, __RC__)
    ptr3d_int_philic = ptr3d_imp

    RETURN_(ESMF_SUCCESS)

  end subroutine Run_data

!-------------------------------------------------------------------------------------
  subroutine aerosol_optics(state, rc)

    implicit none

!   !ARGUMENTS:
    type (ESMF_State)                                :: state
    integer,            intent(out)                  :: rc

!   !Local
    integer, parameter                               :: DP=kind(1.0d0)
    real, dimension(:,:,:), pointer                  :: ple, rh
    real(kind=DP), dimension(:,:,:), pointer         :: var
    real, dimension(:,:,:), pointer                  :: q
    real, dimension(:,:,:,:), pointer                :: q_4d
    integer, allocatable                             :: opaque_self(:)
    type(C_PTR)                                      :: address
    type(CA2G_GridComp), pointer                     :: self

    character (len=ESMF_MAXSTR)                      :: fld_name
    type(ESMF_Field)                                 :: fld
    character (len=ESMF_MAXSTR),allocatable          :: aerosol_names(:)

    real(kind=DP), dimension(:,:,:), allocatable     :: ext_s, ssa_s, asy_s  ! (lon:,lat:,lev:)
    real                                             :: x
    integer                                          :: instance
    integer                                          :: n, nbins
    integer                                          :: i1, j1, i2, j2, km
    integer                                          :: band, offset
    integer, parameter                               :: n_bands = 1

    integer :: i, j, k

    __Iam__('CA2G::aerosol_optics')

!   Begin... 

!   Mie Table instance/index
!   ------------------------
    call ESMF_AttributeGet(state, name='mie_table_instance', value=instance, __RC__)

!   Get aerosol names
!   -----------------
    call ESMF_AttributeGet (state, name='internal_varaible_name', itemCount=nbins, __RC__)
    allocate (aerosol_names(nbins), __STAT__)
    call ESMF_AttributeGet (state, name='internal_varaible_name', valueList=aerosol_names, __RC__)

!   Radiation band
!   --------------
    band = 0
    call ESMF_AttributeGet(state, name='band_for_aerosol_optics', value=band, __RC__)
    offset = band - n_bands

!   Pressure at layer edges 
!   ------------------------
    call ESMF_AttributeGet(state, name='air_pressure_for_aerosol_optics', value=fld_name, __RC__)
    call MAPL_GetPointer(state, ple, trim(fld_name), __RC__)

!    call MAPL_GetPointer (state, ple, 'PLE', __RC__)

    i1 = lbound(ple, 1); i2 = ubound(ple, 1)
    j1 = lbound(ple, 2); j2 = ubound(ple, 2)
                         km = ubound(ple, 3)

!   Relative humidity
!   -----------------
    call ESMF_AttributeGet(state, name='relative_humidity_for_aerosol_optics', value=fld_name, __RC__)
    call MAPL_GetPointer(state, rh, trim(fld_name), __RC__)

!    call MAPL_GetPointer (state, rh, 'RH2', __RC__)

    allocate(ext_s(i1:i2, j1:j2, km), &
             ssa_s(i1:i2, j1:j2, km), &
             asy_s(i1:i2, j1:j2, km), __STAT__)

    allocate(q_4d(i1:i2, j1:j2, km, nbins), __STAT__)

    do n = 1, nbins
       call ESMF_StateGet (state, trim(aerosol_names(n)), field=fld, __RC__)
       call ESMF_FieldGet (fld, farrayPtr=q, __RC__)

        do k = 1, km
           do j = j1, j2
              do i = i1, i2
                 x = ((ple(i,j,k) - ple(i,j,k-1))*0.01)*(100./MAPL_GRAV)
                 q_4d(i,j,k,n) = x * q(i,j,k)
              end do
           end do
        end do
    end do

    call ESMF_AttributeGet(state, name='mieTable_pointer', itemCount=n, __RC__)
    allocate (opaque_self(n), __STAT__)
    call ESMF_AttributeGet(state, name='mieTable_pointer', valueList=opaque_self, __RC__)

    address = transfer(opaque_self, address)
    call c_f_pointer(address, self)

    call mie_ (self%rad_MieTable(instance), nbins, n_bands, offset, q_4d, rh, ext_s, ssa_s, asy_s, __RC__)

    call ESMF_AttributeGet(state, name='extinction_in_air_due_to_ambient_aerosol', value=fld_name, __RC__)
    if (fld_name /= '') then
        call MAPL_GetPointer(state, var, trim(fld_name), __RC__)
        var = ext_s(:,:,:)
    end if

    call ESMF_AttributeGet(state, name='single_scattering_albedo_of_ambient_aerosol', value=fld_name, __RC__)
    if (fld_name /= '') then
        call MAPL_GetPointer(state, var, trim(fld_name), __RC__)
        var = ssa_s(:,:,:)
    end if

   call ESMF_AttributeGet(state, name='asymmetry_parameter_of_ambient_aerosol', value=fld_name, __RC__)
    if (fld_name /= '') then
        call MAPL_GetPointer(state, var, trim(fld_name), __RC__)
        var = asy_s(:,:,:)
    end if

    deallocate(ext_s, ssa_s, asy_s, __STAT__)
    deallocate(q_4d, __STAT__)

    RETURN_(ESMF_SUCCESS)

  contains

    subroutine mie_(mie_table, nbins, nb, offset, q, rh, bext_s, bssa_s, basym_s, rc)

    implicit none

    type(Chem_Mie),                intent(inout) :: mie_table        ! mie table
    integer,                       intent(in   ) :: nbins            ! number of bins
    integer,                       intent(in )   :: nb               ! number of bands
    integer,                       intent(in )   :: offset           ! bands offset 
    real,                          intent(in )   :: q(:,:,:,:)       ! aerosol mass mixing ratio, kg kg-1
    real,                          intent(in )   :: rh(:,:,:)        ! relative humidity
    real(kind=8), intent(  out) :: bext_s (size(ext_s,1),size(ext_s,2),size(ext_s,3))
    real(kind=8), intent(  out) :: bssa_s (size(ext_s,1),size(ext_s,2),size(ext_s,3))
    real(kind=8), intent(  out) :: basym_s(size(ext_s,1),size(ext_s,2),size(ext_s,3))
    integer,           intent(out)  :: rc
    ! local
    integer                           :: l
    real                              :: bext (size(ext_s,1),size(ext_s,2),size(ext_s,3))  ! extinction
    real                              :: bssa (size(ext_s,1),size(ext_s,2),size(ext_s,3))  ! SSA
    real                              :: gasym(size(ext_s,1),size(ext_s,2),size(ext_s,3))  ! asymmetry parameter

    __Iam__('CA2G::aerosol_optics::mie_')

     bext_s  = 0.0d0
     bssa_s  = 0.0d0
     basym_s = 0.0d0

    do l = 1, nbins
       call Chem_MieQuery(mie_table, l, real(offset+1.), q(:,:,:,l), rh, bext, gasym=gasym, ssa=bssa)

       bext_s  = bext_s  +             bext     ! extinction
       bssa_s  = bssa_s  +       (bssa*bext)    ! scattering extinction
       basym_s = basym_s + gasym*(bssa*bext)    ! asymetry parameter multiplied by scatering extiction 

    end do

    RETURN_(ESMF_SUCCESS)

    end subroutine mie_

  end subroutine aerosol_optics

end module CA2G_GridCompMod


