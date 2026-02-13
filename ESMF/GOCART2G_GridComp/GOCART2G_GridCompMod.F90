#include "MAPL.h"

!BOP
!MODULE: GOCART2G_GridCompMod - The GOCART 2nd Generation Aerosol Grid Component
!INTERFACE:
module GOCART2G_GridCompMod

   !USES:
   use ESMF
   use mapl_ErrorHandling, only: MAPL_Verify, MAPL_Return, MAPL_Assert
   use MAPL_Constants, only: MAPL_GRAV, MAPL_PI

   use mapl3g_generic, only: MAPL_GridCompSetEntryPoint, MAPL_GridCompGet, MAPL_GridCompAddSpec
   use mapl3g_generic, only: MAPL_GridCompAddChild, MAPL_GridCompGetChildName, MAPL_GridCompRunChild
   use mapl3g_generic, only: MAPL_GridCompAddConnectivity
   use mapl3g_generic, only: MAPL_GridCompGetResource, MAPL_GridCompReexport
   use mapl3g_generic, only: MAPL_STATEITEM_STATE, MAPL_STATEITEM_FIELDBUNDLE
   use mapl3g_generic, only: MAPL_UserCompGetInternalState, MAPL_UserCompSetInternalState
   use mapl3g_RestartModes, only: MAPL_RESTART_SKIP
   use mapl3g_VerticalStaggerLoc, only: VERTICAL_STAGGER_NONE, VERTICAL_STAGGER_CENTER, VERTICAL_STAGGER_EDGE
   use mapl3g_FieldBundle_API, only: MAPL_FieldBundleAdd, MAPL_FieldBundleGet
   use mapl3g_State_API, only: MAPL_StateGetPointer
   use mapl3g_Geom_API, only: MAPL_GridGet
   use mapl3g_UngriddedDim, only: UngriddedDim
   use gftl2_StringVector, only: StringVector

   use Chem_AeroGeneric

   ! Establish the Childen's SetServices
   use DU2G_GridCompMod, only: DU2G_SetServices  => SetServices
   use SS2G_GridCompMod, only: SS2G_SetServices  => SetServices
   ! use SU2G_GridCompMod, only: SU2G_SetServices  => SetServices
   ! use CA2G_GridCompMod, only: CA2G_SetServices  => SetServices
   ! use NI2G_GridCompMod, only: NI2G_SetServices  => SetServices

   implicit none
   private

   !PUBLIC MEMBER FUNCTIONS:
   public  SetServices

   ! Private State
   type :: Instance
      character(:), allocatable :: name
      logical :: is_active
   end type Instance

   type Constituent
      character(:), allocatable :: name
      type(Instance), allocatable :: instances(:)
      integer :: n_active
   end type Constituent

   type GOCART_State
      private
      type(Constituent) :: DU
      type(Constituent) :: SS
      type(Constituent) :: SU
      type(Constituent) :: CA
      type(Constituent) :: NI
      real, allocatable :: wavelengths_profile(:) ! wavelengths for profile aop [nm]
      real, allocatable :: wavelengths_vertint(:) ! wavelengths for vertically integrated aop [nm]
   end type GOCART_State

   character(*), parameter :: PRIVATE_STATE = "GOCART_STATE"

   !DESCRIPTION:
   !
   !   {\tt GOCART} is a gridded component from the GOCART model and includes
   !  dust, sea salt, sulfates, nitrate, organic and black carbon.

   !REVISION HISTORY:
   !  25feb2005  da Silva   First crack.
   !  19jul2006  da Silva   First separate GOCART component.
   !  14Oct2019  E.Sherman, A.Darmenov, A. da Silva, T. Clune  First attempt at refactoring.

   !EOP

contains

   !BOP
   !IROUTINE: SetServices -- Sets ESMF services for this component
   !INTERFACE:
   subroutine SetServices(gc, rc)

      !ARGUMENTS:
      type (ESMF_GridComp), intent(inout) :: gc ! gridded component
      integer, optional, intent(out) :: rc ! return code

      !DESCRIPTION: This version uses MAPL_GenericSetServices, which sets
      !   the Initialize and Finalize services to generic versions. It also
      !   allocates our instance of a generic state and puts it in the
      !   gridded component (gc). Here we only set the two-stage run method and
      !   declare the data services.

      !REVISION HISTORY:
      !  14oct2019  Sherman, da Silva, Darmenov, Clune - First attempt at refactoring for ESMF compatibility

      !EOP
      type(ESMF_HConfig) :: hconfig
      type(GOCART_State), pointer :: self
      integer, allocatable :: wavelengths_diagmie(:)
      ! logical :: use_threads
      type(Instance), allocatable :: child
      character(len=:), allocatable :: child_items
      type(UngriddedDim) :: ungrd_wavelengths_profile, ungrd_wavelengths_vertint
      integer :: iter, status

      ! Wrap gridcomp's private state and store in gc
      _SET_NAMED_PRIVATE_STATE(gc, GOCART_state, PRIVATE_STATE)

      ! Set the Initialize, Run, Finalize entry points
      call MAPL_GridCompSetEntryPoint(gc, ESMF_Method_Initialize, Initialize,  _RC)
      call MAPL_GridCompSetEntryPoint(gc, ESMF_Method_Run, Run1, phase_name="Run1", _RC)
      call MAPL_GridCompSetEntryPoint(gc, ESMF_Method_Run, Run2, phase_name="Run2", _RC)

      ! Retrieve the private state
      _GET_NAMED_PRIVATE_STATE(gc, GOCART_State, PRIVATE_STATE, self)

      call MAPL_GridCompGetResource(gc, "wavelengths_for_profile_aop_in_nm", self%wavelengths_profile, _RC)
      call MAPL_GridCompGetResource(gc, "wavelengths_for_vertically_integrated_aop_in_nm", self%wavelengths_vertint, _RC)
      call MAPL_GridCompGetResource(gc, "aerosol_monochromatic_optics_wavelength_in_nm_from_LUT", wavelengths_diagmie, _RC)
      ! pchakrab: TODO - Do we re-implement threading?
      ! call MAPL_GridCompGetResource(gc, "use_threads", use_threads, default=.false., _RC)

      ! Defined UngriddedDim items
      ungrd_wavelengths_profile = UngriddedDim( &
           size(self%wavelengths_profile), &
           name="wavelengths_profile", &
           units="nm")
      ungrd_wavelengths_vertint = UngriddedDim( &
           size(self%wavelengths_vertint), &
           name="wavelengths_vertint", &
           units="nm")

      ! ! Get my internal MAPL_Generic state
      ! call MAPL_GetObjectFromGC (GC, MAPL, _RC)
      ! ! set use_threads
      ! call MAPL%set_use_threads(use_threads)

      ! Get instances to determine what children will be born
      ! IMPORTANT: Active instances are created first
      call MAPL_GridCompGet(gc, hconfig=hconfig, _RC)
      call setup_constituents_(self, hconfig, _RC)
      _ASSERT(.not. (self%NI%n_active > 1), "GOCART supports only one active nitrate instance")

      ! Create children's gridded components and invoke their SetServices
      call create_instances_(self, gc, _RC)

      ! Define EXPORT states

      ! State needed by radiation and moist. Contains aerosols and callbacks
      ! CONNECT child's STATE (e.g. SS_AERO) to import state (SS_AERO)
      ! and then copy import state (SS_AERO) to export state AERO
      do iter = 1, size(self%SS%instances)
         child = self%SS%instances(iter)
         call MAPL_GridCompAddSpec( &
              gc, &
              state_intent=ESMF_STATEINTENT_IMPORT, &
              short_name=child%name//"_AERO", &
              standard_name="aerosol_mass_mixing_ratios_ng",  &
              dims="xyz", &
              vstagger=VERTICAL_STAGGER_CENTER, &
              units="kg kg-1", &
              itemtype=MAPL_STATEITEM_STATE, &
              _RC)
      end do
      do iter = 1, size(self%DU%instances)
         child = self%DU%instances(iter)
         call MAPL_GridCompAddSpec( &
              gc, &
              state_intent=ESMF_STATEINTENT_IMPORT, &
              short_name=child%name//"_AERO", &
              standard_name="aerosol_mass_mixing_ratios_ng",  &
              dims="xyz", &
              vstagger=VERTICAL_STAGGER_CENTER, &
              units="kg kg-1", &
              itemtype=MAPL_STATEITEM_STATE, &
              _RC)
      end do
      call MAPL_GridCompAddSpec( &
           gc, &
           state_intent=ESMF_STATEINTENT_EXPORT, &
           short_name="AERO", &
           standard_name="aerosol_mass_mixing_ratios_ng",  &
           dims="xyz", &
           vstagger=VERTICAL_STAGGER_CENTER, &
           units="kg kg-1", &
           itemtype=MAPL_STATEITEM_STATE, &
           _RC)

      ! Bundle needed by surface for snow albedo modification by aerosol settling and deposition
      ! CONNECT child's BUNDLE (e.g. SS_AERO_DP) to import bundle (SS_AERO_DP)
      ! and then copy import bundle (SS_AERO_DP) to export bundle AERO_DP
      do iter = 1, size(self%SS%instances)
         child = self%SS%instances(iter)
         call MAPL_GridCompAddSpec( &
              gc, &
              state_intent=ESMF_STATEINTENT_IMPORT, &
              short_name=child%name//"_AERO_DP", &
              standard_name="aerosol_deposition_ng",  &
              dims="xy", &
              vstagger=VERTICAL_STAGGER_NONE, &
              units="kg m-2 s-1", &
              itemtype=MAPL_STATEITEM_FIELDBUNDLE, &
              _RC)
      end do
      do iter = 1, size(self%SS%instances)
         child = self%DU%instances(iter)
         call MAPL_GridCompAddSpec( &
              gc, &
              state_intent=ESMF_STATEINTENT_IMPORT, &
              short_name=child%name//"_AERO_DP", &
              standard_name="aerosol_deposition_ng",  &
              dims="xy", &
              vstagger=VERTICAL_STAGGER_NONE, &
              units="kg m-2 s-1", &
              itemtype=MAPL_STATEITEM_FIELDBUNDLE, &
              _RC)
      end do
      call MAPL_GridCompAddSpec( &
           gc, &
           state_intent=ESMF_STATEINTENT_EXPORT, &
           short_name="AERO_DP", &
           standard_name="aerosol_deposition_ng",  &
           dims="xy", &
           vstagger=VERTICAL_STAGGER_NONE, &
           units="kg m-2 s-1", &
           itemtype=MAPL_STATEITEM_FIELDBUNDLE, &
           _RC)


#include "GOCART2G_Export___.h"
#include "GOCART2G_Import___.h"

      ! pchakrab: TODO - NEEDS PORTING - ACTIVATE ONCE SU HAS BEEN PORTED
      ! ! Allow children of Chemistry to connect to these fields
      ! if ((self%SU%instances(1)%is_active)) then
      !    call MAPL_AddExportSpec (GC, SHORT_NAME='PSO4', CHILD_ID=self%SU%instances(1)%id, _RC)
      ! end if

      ! pchakrab: TODO - ACTIVATE ONCE NI HAS BEEN PORTED
      ! ! Add connectivities for Nitrate component
      ! ! Nitrate currently only supports one Nitrate component. Nitrate only
      ! ! uses the first active dust and sea salt instance.
      ! if (size(self%NI%instances) > 0) then
      !    if ((self%DU%instances(1)%is_active)) then
      !       call MAPL_GridCompAddConnectivity( &
      !            gc, &
      !            src_comp=self%DU%instances(1)%name, &
      !            src_names="DU", &
      !            dst_comp=self%NI%instances(1)%name, &
      !            _RC)
      !    end if

      !    if ((self%SS%instances(1)%is_active)) then
      !       call MAPL_GridCompAddConnectivity( &
      !            gc, &
      !            src_comp=self%SS%instances(1)%name, &
      !            src_names="SS", &
      !            dst_comp=self%NI%instances(1)%name, &
      !            _RC)
      !    end if

      !    if ((self%SU%instances(1)%is_active)) then
      !       call MAPL_GridCompAddConnectivity( &
      !            gc, &
      !            src_comp=self%SU%instances(1)%name, &
      !            src_names="SO4", &
      !            dst_comp=self%NI%instances(1)%name, &
      !            _RC)
      !    end if
      ! end if

      ! Connections to Sea Salt's export items
      child_items = &
           "SSEXTTAU, SSSTEXTTAU, SSSCATAU, SSSTSCATAU, " // &
           "SSEXTCOEF, SSEXTCOEFRH20, SSEXTCOEFRH80, " // &
           "SSSCACOEF, SSSCACOEFRH20, SSSCACOEFRH80, " // &
           "SSBCKCOEF, SSEXTT25, SSSCAT25, SSEXTTFM, SSSCATFM, " // &
           "SSANGSTR, SSSMASS, SSSMASS25"
      do iter = 1, size(self%SS%instances)
         child = self%SS%instances(iter)
         if ((child%is_active) .and. (index(child%name, "data") == 0 )) then
            call MAPL_GridCompAddConnectivity( &
                 gc, &
                 src_comp=child%name, &
                 src_names=child_items, &
                 dst_comp="<self>", _RC)
            ! AERO
            call MAPL_GridCompAddConnectivity( &
                 gc, &
                 src_comp=child%name, &
                 src_names=child%name//"_AERO", &
                 dst_comp="<self>", &
                 dst_names=child%name//"_AERO", _RC)
            ! AERO_DP
            call MAPL_GridCompAddConnectivity( &
                 gc, &
                 src_comp=child%name, &
                 src_names=child%name//"_AERO_DP", &
                 dst_comp="<self>", &
                 dst_names=child%name//"_AERO_DP", _RC)
         end if
      end do

      ! Connections to Dust's export items
      child_items= &
           "DUEXTTAU, DUSTEXTTAU, DUSCATAU, DUSTSCATAU, " // &
           "DUEXTCOEF, DUEXTCOEFRH20, DUEXTCOEFRH80, " // &
           "DUSCACOEF, DUSCACOEFRH20, DUSCACOEFRH80, " // &
           "DUBCKCOEF, DUEXTT25, DUSCAT25, DUEXTTFM, DUSCATFM, " // &
           "DUANGSTR, DUSMASS, DUSMASS25"
      do iter = 1, size(self%DU%instances)
         child = self%DU%instances(iter)
         if ((child%is_active) .and. (index(child%name, "data") == 0 )) then
            call MAPL_GridCompAddConnectivity( &
                 gc, &
                 src_comp=child%name, &
                 src_names=child_items, &
                 dst_comp="<self>", _RC)
            ! AERO
            call MAPL_GridCompAddConnectivity( &
                 gc, &
                 src_comp=child%name, &
                 src_names=child%name//"_AERO", &
                 dst_comp="<self>", &
                 dst_names=child%name//"_AERO", _RC)
            ! AERO_DP
            call MAPL_GridCompAddConnectivity( &
                 gc, &
                 src_comp=child%name, &
                 src_names=child%name//"_AERO_DP", &
                 dst_comp="<self>", &
                 dst_names=child%name//"_AERO_DP", _RC)
         end if
      end do

      _RETURN(_SUCCESS)

   end subroutine SetServices

   !BOP
   !IROUTINE: Initialize -- Initialize method for the composite Gridded Component
   !INTERFACE:
   subroutine Initialize(gc, import, export, clock, rc)

      !ARGUMENTS:
      type(ESMF_GridComp) :: gc  ! Gridded component
      type(ESMF_State) :: import ! Import state
      type(ESMF_State) :: export ! Export state
      type(ESMF_Clock) :: clock  ! The clock
      integer, intent(out) :: rc ! Error code

      !DESCRIPTION:  This initializes the GOCART Grid Component. It primarily creates
      !                its exports and births its children.
      !REVISION HISTORY:
      ! 14oct2019   E.Sherman  First attempt at refactoring
      !EOP

      type(ESMF_Geom) :: geom
      type(ESMF_Grid) :: grid
      type(ESMF_State) :: aero
      type(ESMF_FieldBundle) :: aero_dp
      type(ESMF_Info) :: info
      type(GOCART_State), pointer :: self
      character(len=ESMF_MAXSTR), allocatable :: aero_aci_modes(:)
      real :: maxclean, ccntuning
      integer :: im, jm, km, status

      call MAPL_GridCompGet(gc, geom=geom, grid=grid, num_levels=km, _RC)
      call MAPL_GridGet(grid, im=im, jm=jm, _RC)

      ! Get my private state
      _GET_NAMED_PRIVATE_STATE(gc, GOCART_State, PRIVATE_STATE, self)

      ! Fill AERO_RAD, AERO_ACI, and AERO_DP with the children's states
      call ESMF_StateGet(export, "AERO", aero, _RC)
      call ESMF_StateGet(export, "AERO_DP", aero_dp, _RC)

      ! ! Add children's AERO states to GOCART2G's AERO states
      ! ! Only active instances are passed to radiation
      call add_aero_states_(self%DU%instances(:))
      call add_aero_states_(self%SS%instances(:))
      ! call add_aero_states_(self%SU%instances(:))
      ! call add_aero_states_(self%CA%instances(:))
      ! call add_aero_states_(self%NI%instances(:))

      ! Begin AERO_RAD
      ! Add variables to AERO_RAD state. Used in aerosol optics calculations
      call add_aero(aero, label="air_pressure_for_aerosol_optics", label2="PLE", geom=geom, km=km, _RC)
      call add_aero(aero, label="relative_humidity_for_aerosol_optics", label2="RH", geom=geom, km=km, _RC)
      call add_aero(aero, label="extinction_in_air_due_to_ambient_aerosol", label2="EXT", geom=geom, km=km, _RC)
      call add_aero(aero, label="single_scattering_albedo_of_ambient_aerosol", label2="SSA", geom=geom, km=km, _RC)
      call add_aero(aero, label="asymmetry_parameter_of_ambient_aerosol", label2="ASY", geom=geom, km=km, _RC)
      call add_aero( &
           aero, &
           label="monochromatic_extinction_in_air_due_to_ambient_aerosol", label2="monochromatic_EXT", &
           geom=geom, _RC)

      ! Used in get_mixRatioSum
      call add_aero(aero, label="sum_of_internalState_aerosol_DU", label2="aerosolSumDU", geom=geom, km=km, _RC)
      call add_aero(aero, label="sum_of_internalState_aerosol_SS", label2="aerosolSumSS", geom=geom, km=km, _RC)
      call add_aero(aero, label="sum_of_internalState_aerosol_NI", label2="aerosolSumNI", geom=geom, km=km, _RC)
      call add_aero(aero, label="sum_of_internalState_aerosol_CA.oc", label2="aerosolSumCA.oc", geom=geom, km=km, _RC)
      call add_aero(aero, label="sum_of_internalState_aerosol_CA.bc", label2="aerosolSumCA.bc", geom=geom, km=km, _RC)
      call add_aero(aero, label="sum_of_internalState_aerosol_CA.br", label2="aerosolSumCA.br", geom=geom, km=km, _RC)
      call add_aero(aero, label="sum_of_internalState_aerosol_SU", label2="aerosolSumSU", geom=geom, km=km, _RC)

      call ESMF_InfoGetFromHost(aero, info, _RC)
      call ESMF_InfoSet(info, key="band_for_aerosol_optics", value=0, _RC)
      call ESMF_InfoSet(info, key="wavelength_for_aerosol_optics", value=0., _RC)
      call ESMF_InfoSet(info, key="aerosolName", value="", _RC)
      call ESMF_InfoSet(info, key="im", value=im, _RC)
      call ESMF_InfoSet(info, key="jm", value=jm, _RC)
      call ESMF_InfoSet(info, key="km", value=km, _RC)

      ! Attach method to return sum of aerosols. Used in GAAS.
      call ESMF_MethodAdd(aero, label="get_mixRatioSum", userRoutine=get_mixRatioSum, _RC)

      ! Attach method to create a Bundle of aerosol fields. Used in GAAS.
      call ESMF_MethodAdd(aero, label="serialize_bundle", userRoutine=serialize_bundle, _RC)

      ! Attach the monochromatic aerosol optics method. Used in GAAS.
      call ESMF_MethodAdd(aero, label="get_monochromatic_aop", userRoutine=get_monochromatic_aop, _RC)

      ! Attach the aerosol optics method. Used in Radiation.
      call ESMF_MethodAdd(aero, label="run_aerosol_optics", userRoutine=run_aerosol_optics, _RC)

      ! This attribute indicates if the aerosol optics method is implemented or not.
      ! Radiation will not call the aerosol optics method unless this attribute is
      ! explicitly set to true.
      call ESMF_InfoSet(info, key="implements_aerosol_optics_method", value=.true., _RC)

      ! Begin adding necessary aerosol cloud interaction information
      aero_aci_modes =  [ &
           "du001    ", "du002    ", "du003    ", &
           "du004    ", "du005    ",              &
           "ss001    ", "ss002    ", "ss003    ", &
           "sulforg01", "sulforg02", "sulforg03", &
           "bcphilic ", "ocphilic ", "brcphilic"]

      call ESMF_InfoSet(info, key="aerosol_modes", values=aero_aci_modes, _RC)

      ! max mixing ratio before switching to "polluted" size distributions
      call MAPL_GridCompGetResource(gc, "MAXCLEAN", maxclean, default=1.0e-9, _RC)
      call ESMF_InfoSet(info, key="max_q_clean", value=maxclean, _RC)

      ! call ESMF_ConfigGetAttribute(CF, CCNtuning, default=1.8, label="CCNTUNING:", _RC)
      call MAPL_GridCompGetResource(gc, "CCNTUNING", ccntuning, default=1.8, _RC)
      call ESMF_InfoSet(info, key="ccn_tuning", value=ccntuning, _RC)

      ! Add variables to AERO state
      call add_aero(aero, label="air_temperature", label2="T", geom=geom, km=km, _RC)
      call add_aero(aero, label="fraction_of_land_type", label2="FRLAND", geom=geom, _RC)
      call add_aero(aero, label="width_of_aerosol_mode", label2="SIGMA", geom=geom, km=km, _RC)
      call add_aero(aero, label="aerosol_number_concentration", label2="NUM", geom=geom, km=km, _RC)
      call add_aero(aero, label="aerosol_dry_size", label2="DGN", geom=geom, km=km, _RC)
      call add_aero(aero, label="aerosol_density", label2="density", geom=geom, km=km, _RC)
      call add_aero(aero, label="aerosol_hygroscopicity", label2="KAPPA", geom=geom, km=km, _RC)
      call add_aero(aero, label="fraction_of_dust_aerosol", label2="FDUST", geom=geom, km=km, _RC)
      call add_aero(aero, label="fraction_of_soot_aerosol", label2="FSOOT", geom=geom, km=km, _RC)
      call add_aero(aero, label="fraction_of_organic_aerosol", label2="FORGANIC", geom=geom, km=km, _RC)

      ! Attach the aerosol optics method
      call ESMF_MethodAdd(aero, label="aerosol_activation_properties", userRoutine=aerosol_activation_properties, _RC)

      ! call ESMF_StatePrint(aero, _RC)
      _RETURN(_SUCCESS)
      _UNUSED_DUMMY(clock)

   contains

      subroutine add_aero_states_(instances)
         type(Instance), intent(in) :: instances(:)

         type(Instance) :: instance_
         type (ESMF_Field), allocatable :: field_list_(:)
         type(ESMF_FieldBundle) :: bundle_
         type(ESMF_State) :: state_
         integer :: iter, status

         do iter = 1, size(instances)
            instance_ = instances(iter)
            if (.not. instance_%is_active) cycle

            ! AERO state
            call ESMF_StateGet(import, instance_%name//"_AERO", state_, _RC)
            call ESMF_StateAdd(aero, [state_], _RC)

            ! AERO_DP bundle
            if (instance_%name(1:2) /= "NI") then
               call ESMF_StateGet(import, instance_%name//"_AERO_DP", bundle_, _RC)
               call MAPL_FieldBundleGet(bundle_, fieldList=field_list_, _RC)
               call ESMF_FieldBundleAdd(aero_dp, field_list_, _RC)
               deallocate(field_list_)
            end if
         end do

         _RETURN(_SUCCESS)
      end subroutine add_aero_states_

   end subroutine Initialize

   !BOP
   !IROUTINE: RUN -- Run method for GOCART2G
   !INTERFACE:
   subroutine Run1 (gc, import, export, clock, RC)
      !ARGUMENTS:
      type (ESMF_GridComp) :: gc  ! Gridded component
      type (ESMF_State) :: import ! Import state
      type (ESMF_State) :: export ! Export state
      type (ESMF_Clock) :: clock  ! The clock
      integer, intent(out) :: RC  ! Error code:

      !DESCRIPTION: Run method
      !EOP
      character(len=:), allocatable :: child_name
      integer :: num_children, iter, status

      call MAPL_GridCompGet(gc, num_children=num_children, _RC)
      do iter = 1, num_children
         child_name = MAPL_GridCompGetChildName(gc, iter, _RC)
         if ((index(child_name, "data")) /= 0) cycle
         call MAPL_GridCompRunChild(gc, child_name, phase_name="Run1", _RC)
      end do

      _RETURN(_SUCCESS)
      _UNUSED_DUMMY(import)
      _UNUSED_DUMMY(export)
      _UNUSED_DUMMY(clock)
   end subroutine Run1

   !BOP
   !IROUTINE: RUN2 -- Run2 method for GOCART2G component
   !INTERFACE:
   subroutine Run2(gc, import, export, clock, rc)

      !ARGUMENTS:
      type(ESMF_GridComp) :: gc  ! Gridded component
      type(ESMF_State) :: import ! Import state
      type(ESMF_State) :: export ! Export state
      type(ESMF_Clock) :: clock  ! The clock
      integer, intent(out) :: rc ! Error code:

      !DESCRIPTION: This version uses the MAPL\_GenericSetServices. This function sets
      !                the Initialize and Finalize services, as well as allocating
      !EOP

      type(ESMF_Grid) :: grid
      type (GOCART_State), pointer :: self

      ! ! Nitrates - ACG will generate this once we add NI's export states as GOCART's import
      ! real, pointer, dimension(:,:,:) :: niexttau, nistexttau
      ! real, pointer, dimension(:,:,:) :: niscatau, nistscatau
      ! real, pointer, dimension(:,:,:) :: niextt25, niscat25
      ! real, pointer, dimension(:,:,:) :: niexttfm, niscatfm
      ! real, pointer, dimension(:,:,:,:) :: niextcoef, niscacoef
      ! real, pointer, dimension(:,:,:,:) :: niextcoefrh20, niextcoefrh80
      ! real, pointer, dimension(:,:,:,:) :: niscacoefrh20, niscacoefrh80
      ! real, pointer, dimension(:,:,:,:) :: nibckcoef
      ! real, pointer, dimension(:,:)   :: niangstr, nismass, nismass25
      ! real, pointer, dimension(:,:)   :: nh4smass

      ! ! Sulfates - ACG will generate this once we add NI's export states as GOCART's import
      ! real, pointer, dimension(:,:,:) :: suexttau, sustexttau
      ! real, pointer, dimension(:,:,:) :: suscatau, sustscatau
      ! real, pointer, dimension(:,:,:,:) :: suextcoef, suscacoef
      ! real, pointer, dimension(:,:,:,:) :: suextcoefrh20, suextcoefrh80
      ! real, pointer, dimension(:,:,:,:) :: suscacoefrh20, suscacoefrh80
      ! real, pointer, dimension(:,:,:,:) :: subckcoef
      ! real, pointer, dimension(:,:)   :: suangstr, so4smass

      ! ! Carbonaceous - ACG will generate this once we add CA's export states as GOCART's import
      ! real, pointer, dimension(:,:,:) :: bcexttau, bcstexttau, bcscatau, bcstscatau
      ! real, pointer, dimension(:,:,:,:) :: bcextcoef, bcscacoef
      ! real, pointer, dimension(:,:,:,:) :: bcextcoefrh20, bcextcoefrh80
      ! real, pointer, dimension(:,:,:,:) :: bcscacoefrh20, bcscacoefrh80
      ! real, pointer, dimension(:,:,:,:) :: bcbckcoef
      ! real, pointer, dimension(:,:)   :: bcangstr, bcsmass
      ! real, pointer, dimension(:,:,:) :: ocexttau, ocstexttau, ocscatau, ocstscatau
      ! real, pointer, dimension(:,:,:,:) :: ocextcoef, ocscacoef
      ! real, pointer, dimension(:,:,:,:) :: ocextcoefrh20, ocextcoefrh80
      ! real, pointer, dimension(:,:,:,:) :: ocscacoefrh20, ocscacoefrh80
      ! real, pointer, dimension(:,:,:,:) :: ocbckcoef
      ! real, pointer, dimension(:,:)   :: ocangstr, ocsmass
      ! real, pointer, dimension(:,:,:) :: brexttau, brstexttau, brscatau, brstscatau
      ! real, pointer, dimension(:,:,:,:) :: brextcoef, brscacoef
      ! real, pointer, dimension(:,:,:,:) :: brextcoefrh20, brextcoefrh80
      ! real, pointer, dimension(:,:,:,:) :: brscacoefrh20, brscacoefrh80
      ! real, pointer, dimension(:,:,:,:) :: brbckcoef
      ! real, pointer, dimension(:,:)   :: brangstr, brsmass

      ! real, pointer, dimension(:,:,:) :: pso4
      real, allocatable :: tau1(:,:), tau2(:,:)
      real :: c1, c2, c3
      real, parameter :: pi = 3.141529265 ! pchakrab: TODO - use MAPL_PI instead??
      integer :: ind550
      integer :: im, jm
      character(len=:), allocatable :: child_name
      integer :: n, w, num_children, iter, status

#include "GOCART2G_DeclarePointer___.h"

      call MAPL_GridCompGet(gc, num_children=num_children, _RC)

      ! Run zero Klid for children
      do iter = 1, num_children
         child_name = MAPL_GridCompGetChildName(gc, iter, _RC)
         if ((index(child_name, "data")) /= 0) cycle
         call MAPL_GridCompRunChild(gc, child_name, phase_name="Run0", _RC)
      end do

      ! Get private state
      _GET_NAMED_PRIVATE_STATE(gc, GOCART_State, PRIVATE_STATE, self)

#include "GOCART2G_GetPointer___.h"
      if(associated(totexttau)) totexttau = 0.
      if(associated(totstexttau)) totstexttau = 0.
      if(associated(totscatau)) totscatau = 0.
      if(associated(totstscatau)) totstscatau = 0.
      if(associated(totextt25)) totextt25 = 0.
      if(associated(totscat25)) totscat25 = 0.
      if(associated(totexttfm)) totexttfm = 0.
      if(associated(totscatfm)) totscatfm = 0.
      if(associated(totextcoef)) totextcoef = 0.
      if(associated(totextcoefrh20)) totextcoefrh20 = 0.
      if(associated(totextcoefrh80)) totextcoefrh80 = 0.
      if(associated(totscacoef)) totscacoef = 0.
      if(associated(totscacoefrh20)) totscacoefrh20 = 0.
      if(associated(totscacoefrh80)) totscacoefrh80 = 0.
      if(associated(totbckcoef)) totbckcoef = 0.
      if(associated(totabcktoa)) totabcktoa = 0.
      if(associated(totabcksfc)) totabcksfc = 0.
      if(associated(pm)) pm(:,:) = 0.
      if(associated(pm25)) pm25(:,:) = 0.
      if(associated(pm_rh35)) pm_rh35(:,:) = 0.
      if(associated(pm25_rh35)) pm25_rh35(:,:) = 0.
      if(associated(pm_rh50)) pm_rh50(:,:) = 0.
      if(associated(pm25_rh50)) pm25_rh50(:,:) = 0.
      if(associated(pso4tot)) pso4tot(:,:,:) = 0.

      ! Run the children
      do iter = 1, num_children
         child_name = MAPL_GridCompGetChildName(gc, iter, _RC)
         if ((index(child_name, "data")) /= 0) cycle
         call MAPL_GridCompRunChild(gc, child_name, phase_name="Run2", _RC)
      end do

      ! Compute total aerosol diagnostic values for export
      call MAPL_GridCompGet(gc, grid=grid, _RC)
      call MAPL_GridGet(grid, im=im, jm=jm, _RC)
      if(associated(totangstr)) then
         ind550 = 0
         do w = 1, size(self%wavelengths_vertint) ! find index for 550nm to compute total angstrom
            if ((self%wavelengths_vertint(w)*1.e-9 .ge. 5.49e-7) .and. &
                 (self%wavelengths_vertint(w)*1.e-9 .le. 5.51e-7)) then
               ind550 = w
               exit
            end if
         end do
         _ASSERT(ind550 /= 0, "Cannot produce TOTANGSTR variable without 550nm wavelength")
         totangstr = 0.0
         allocate(tau1(im, jm), tau2(im, jm), _STAT)
         tau1(:,:) = tiny(1.0)
         tau2(:,:) = tiny(1.0)
         c1 = -log(470./550.)
         c2 = -log(870./550.)
         c3 = -log(470./870.)
      end if

      ! Dust
      do n = 1, size(self%DU%instances)
         if ((self%DU%instances(n)%is_active) .and. (index(self%DU%instances(n)%name, 'data') == 0 )) then
            if(associated(totexttau)) totexttau = totexttau + duexttau
            if(associated(totstexttau)) totstexttau = totstexttau+dustexttau
            if(associated(totscatau)) totscatau = totscatau+duscatau
            if(associated(totstscatau)) totstscatau = totstscatau+dustscatau
            if(associated(totextt25)) totextt25 = totextt25+duextt25
            if(associated(totscat25)) totscat25 = totscat25+duscat25
            if(associated(totexttfm)) totexttfm = totexttfm+duexttfm
            if(associated(totscatfm)) totscatfm = totscatfm+duscatfm
            if(associated(totextcoef)) totextcoef = totextcoef+duextcoef
            if(associated(totextcoefrh20)) totextcoefrh20 = totextcoefrh20+duextcoefrh20
            if(associated(totextcoefrh80)) totextcoefrh80 = totextcoefrh80+duextcoefrh80
            if(associated(totscacoef)) totscacoef = totscacoef+duscacoef
            if(associated(totscacoefrh20)) totscacoefrh20 = totscacoefrh20+duscacoefrh20
            if(associated(totscacoefrh80)) totscacoefrh80 = totscacoefrh80+duscacoefrh80
            if(associated(totbckcoef)) totbckcoef = totbckcoef+dubckcoef

            if(associated(pm)        .and. associated(dusmass))   pm        = pm        + dusmass
            if(associated(pm25)      .and. associated(dusmass25)) pm25      = pm25      + dusmass25
            if(associated(pm_rh35)   .and. associated(dusmass))   pm_rh35   = pm_rh35   + dusmass
            if(associated(pm25_rh35) .and. associated(dusmass25)) pm25_rh35 = pm25_rh35 + dusmass25
            if(associated(pm_rh50)   .and. associated(dusmass))   pm_rh50   = pm_rh50   + dusmass
            if(associated(pm25_rh50) .and. associated(dusmass25)) pm25_rh50 = pm25_rh50 + dusmass25

            if(associated(totangstr) .and. associated(duexttau) .and. associated(duangstr)) then
               tau1 = tau1 + duexttau(:,:,ind550)*exp(c1*duangstr)
               tau2 = tau2 + duexttau(:,:,ind550)*exp(c2*duangstr)
            end if
         end if
      end do

      ! Sea Salt
      do n = 1, size(self%SS%instances)
         if ((self%SS%instances(n)%is_active) .and. (index(self%SS%instances(n)%name, "data") == 0 )) then
            if(associated(totexttau)) totexttau = totexttau + ssexttau
            if(associated(totstexttau)) totstexttau = totstexttau + ssstexttau
            if(associated(totscatau)) totscatau = totscatau + ssscatau
            if(associated(totstscatau)) totstscatau = totstscatau + ssstscatau
            if(associated(totextt25)) totextt25 = totextt25 + ssextt25
            if(associated(totscat25)) totscat25 = totscat25 + ssscat25
            if(associated(totexttfm)) totexttfm = totexttfm + ssexttfm
            if(associated(totscatfm)) totscatfm = totscatfm + ssscatfm
            if(associated(totextcoef)) totextcoef = totextcoef + ssextcoef
            if(associated(totextcoefrh20)) totextcoefrh20 = totextcoefrh20 + ssextcoefrh20
            if(associated(totextcoefrh80)) totextcoefrh80 = totextcoefrh80 + ssextcoefrh80
            if(associated(totscacoef)) totscacoef = totscacoef + ssscacoef
            if(associated(totscacoefrh20)) totscacoefrh20 = totscacoefrh20 + ssscacoefrh20
            if(associated(totscacoefrh80)) totscacoefrh80 = totscacoefrh80 + ssscacoefrh80
            if(associated(totbckcoef)) totbckcoef = totbckcoef + ssbckcoef

            if(associated(pm)        .and. associated(sssmass))   pm        = pm        + sssmass
            if(associated(pm25)      .and. associated(sssmass25)) pm25      = pm25      + sssmass25
            if(associated(pm_rh35)   .and. associated(sssmass))   pm_rh35   = pm_rh35   + 1.86*sssmass
            if(associated(pm25_rh35) .and. associated(sssmass25)) pm25_rh35 = pm25_rh35 + 1.86*sssmass25
            if(associated(pm_rh50)   .and. associated(sssmass))   pm_rh50   = pm_rh50   + 2.42*sssmass
            if(associated(pm25_rh50) .and. associated(sssmass25)) pm25_rh50 = pm25_rh50 + 2.42*sssmass25

            if(associated(totangstr) .and. associated(ssexttau) .and. associated(ssangstr)) then
               tau1 = tau1 + ssexttau(:,:,ind550) * exp(c1*ssangstr)
               tau2 = tau2 + ssexttau(:,:,ind550) * exp(c2*ssangstr)
            end if
         end if
      end do

      ! ! Nitrates - NOTE! Nitrates currently only support one active instance
      ! do n = 1, size(self%NI%instances)
      !    if ((self%NI%instances(n)%is_active) .and. (index(self%NI%instances(n)%name, 'data') == 0 )) then
      !       call MAPL_GetPointer (gex(self%NI%instances(n)%id), niexttau, 'NIEXTTAU', _RC)
      !       call MAPL_GetPointer (gex(self%NI%instances(n)%id), nistexttau, 'NISTEXTTAU', _RC)
      !       call MAPL_GetPointer (gex(self%NI%instances(n)%id), niscatau, 'NISCATAU', _RC)
      !       call MAPL_GetPointer (gex(self%NI%instances(n)%id), nistscatau, 'NISTSCATAU', _RC)
      !       call MAPL_GetPointer (gex(self%NI%instances(n)%id), niextcoef, 'NIEXTCOEF', _RC)
      !       call MAPL_GetPointer (gex(self%NI%instances(n)%id), niextcoefrh20, 'NIEXTCOEFRH20', _RC)
      !       call MAPL_GetPointer (gex(self%NI%instances(n)%id), niextcoefrh80, 'NIEXTCOEFRH80', _RC)
      !       call MAPL_GetPointer (gex(self%NI%instances(n)%id), niscacoef, 'NISCACOEF', _RC)
      !       call MAPL_GetPointer (gex(self%NI%instances(n)%id), niscacoefrh20, 'NISCACOEFRH20', _RC)
      !       call MAPL_GetPointer (gex(self%NI%instances(n)%id), niscacoefrh80, 'NISCACOEFRH80', _RC)
      !       call MAPL_GetPointer (gex(self%NI%instances(n)%id), nibckcoef, 'NIBCKCOEF', _RC)
      !       call MAPL_GetPointer (gex(self%NI%instances(n)%id), niextt25, 'NIEXTT25', _RC)
      !       call MAPL_GetPointer (gex(self%NI%instances(n)%id), niscat25, 'NISCAT25', _RC)
      !       call MAPL_GetPointer (gex(self%NI%instances(n)%id), niexttfm, 'NIEXTTFM', _RC)
      !       call MAPL_GetPointer (gex(self%NI%instances(n)%id), niscatfm, 'NISCATFM', _RC)
      !       call MAPL_GetPointer (gex(self%NI%instances(n)%id), niangstr, 'NIANGSTR', _RC)

      !       ! Iterate over the wavelengths
      !       do w = 1, size(self%wavelengths_vertint)
      !          if(associated(totexttau) .and. associated(niexttau)) totexttau(:,:,w) = totexttau(:,:,w)+niexttau(:,:,w)
      !          if(associated(totstexttau) .and. associated(nistexttau)) totstexttau(:,:,w) = totstexttau(:,:,w)+nistexttau(:,:,w)
      !          if(associated(totscatau) .and. associated(niscatau)) totscatau(:,:,w) = totscatau(:,:,w)+niscatau(:,:,w)
      !          if(associated(totstscatau) .and. associated(nistscatau)) totstscatau(:,:,w) = totstscatau(:,:,w)+nistscatau(:,:,w)
      !          if(associated(totextt25) .and. associated(niextt25)) totextt25(:,:,w) = totextt25(:,:,w)+niextt25(:,:,w)
      !          if(associated(totscat25) .and. associated(niscat25)) totscat25(:,:,w) = totscat25(:,:,w)+niscat25(:,:,w)
      !          if(associated(totexttfm) .and. associated(niexttfm)) totexttfm(:,:,w) = totexttfm(:,:,w)+niexttfm(:,:,w)
      !          if(associated(totscatfm) .and. associated(niscatfm)) totscatfm(:,:,w) = totscatfm(:,:,w)+niscatfm(:,:,w)
      !       end do

      !       do w = 1, size(self%wavelengths_profile)
      !          if(associated(totextcoef) .and. associated(niextcoef)) totextcoef(:,:,:,w) = totextcoef(:,:,:,w)+niextcoef(:,:,:,w)
      !          if(associated(totextcoefrh20) .and. associated(niextcoefrh20)) totextcoefrh20(:,:,:,w) = totextcoefrh20(:,:,:,w)+niextcoefrh20(:,:,:,w)
      !          if(associated(totextcoefrh80) .and. associated(niextcoefrh80)) totextcoefrh80(:,:,:,w) = totextcoefrh80(:,:,:,w)+niextcoefrh80(:,:,:,w)
      !          if(associated(totscacoef) .and. associated(niscacoef)) totscacoef(:,:,:,w) = totscacoef(:,:,:,w)+niscacoef(:,:,:,w)
      !          if(associated(totscacoefrh20) .and. associated(niscacoefrh20)) totscacoefrh20(:,:,:,w) = totscacoefrh20(:,:,:,w)+niscacoefrh20(:,:,:,w)
      !          if(associated(totscacoefrh80) .and. associated(niscacoefrh80)) totscacoefrh80(:,:,:,w) = totscacoefrh80(:,:,:,w)+niscacoefrh80(:,:,:,w)
      !          if(associated(totbckcoef) .and. associated(nibckcoef)) totbckcoef(:,:,:,w) = totbckcoef(:,:,:,w)+nibckcoef(:,:,:,w)
      !       end do

      !       call MAPL_GetPointer (gex(self%NI%instances(n)%id), nismass,   'NISMASS',   _RC)
      !       call MAPL_GetPointer (gex(self%NI%instances(n)%id), nismass25, 'NISMASS25', _RC)
      !       call MAPL_GetPointer (gex(self%NI%instances(n)%id), nh4smass,  'NH4SMASS',   _RC)
      !       if(associated(pm)        .and. associated(nismass)   .and. associated(nh4smass)) pm        = pm   + nismass   + nh4smass
      !       if(associated(pm25)      .and. associated(nismass25) .and. associated(nh4smass)) pm25      = pm25 + nismass25 + nh4smass
      !       if(associated(pm_rh35)   .and. associated(nismass)   .and. associated(nh4smass)) pm_rh35   = pm_rh35   + 1.33*(nismass   + nh4smass)
      !       if(associated(pm25_rh35) .and. associated(nismass25) .and. associated(nh4smass)) pm25_rh35 = pm25_rh35 + 1.33*(nismass25 + nh4smass)
      !       if(associated(pm_rh50)   .and. associated(nismass)   .and. associated(nh4smass)) pm_rh50   = pm_rh50   + 1.51*(nismass   + nh4smass)
      !       if(associated(pm25_rh50) .and. associated(nismass25) .and. associated(nh4smass)) pm25_rh50 = pm25_rh50 + 1.51*(nismass25 + nh4smass)

      !       if(associated(totangstr) .and. associated(niexttau) .and. associated(niangstr)) then
      !          tau1 = tau1 + niexttau(:,:,ind550)*exp(c1*niangstr)
      !          tau2 = tau2 + niexttau(:,:,ind550)*exp(c2*niangstr)
      !       end if
      !    end if
      ! end do

      ! ! Sulfates
      ! nifactor = 132.14/96.06
      ! if (size(self%NI%instances) > 0) then
      !    if ((self%NI%instances(1)%is_active) .and. (index(self%NI%instances(1)%name, 'data') == 0 )) nifactor = 1.0
      ! end if

      ! do n = 1, size(self%SU%instances)
      !    if ((self%SU%instances(n)%is_active) .and. (index(self%SU%instances(n)%name, 'data') == 0 )) then
      !       call MAPL_GetPointer (gex(self%SU%instances(n)%id), suexttau, 'SUEXTTAU', _RC)
      !       call MAPL_GetPointer (gex(self%SU%instances(n)%id), suextcoef, 'SUEXTCOEF', _RC)
      !       call MAPL_GetPointer (gex(self%SU%instances(n)%id), suextcoefrh20, 'SUEXTCOEFRH20', _RC)
      !       call MAPL_GetPointer (gex(self%SU%instances(n)%id), suextcoefrh80, 'SUEXTCOEFRH80', _RC)
      !       call MAPL_GetPointer (gex(self%SU%instances(n)%id), suscacoef, 'SUSCACOEF', _RC)
      !       call MAPL_GetPointer (gex(self%SU%instances(n)%id), suscacoefrh20, 'SUSCACOEFRH20', _RC)
      !       call MAPL_GetPointer (gex(self%SU%instances(n)%id), suscacoefrh80, 'SUSCACOEFRH80', _RC)
      !       call MAPL_GetPointer (gex(self%SU%instances(n)%id), subckcoef, 'SUBCKCOEF', _RC)
      !       call MAPL_GetPointer (gex(self%SU%instances(n)%id), sustexttau, 'SUSTEXTTAU', _RC)
      !       call MAPL_GetPointer (gex(self%SU%instances(n)%id), suscatau, 'SUSCATAU', _RC)
      !       call MAPL_GetPointer (gex(self%SU%instances(n)%id), sustscatau, 'SUSTSCATAU', _RC)
      !       call MAPL_GetPointer (gex(self%SU%instances(n)%id), suangstr, 'SUANGSTR', _RC)

      !     !   Iterate over the wavelengths
      !       do w = 1, size(self%wavelengths_vertint)
      !          if(associated(totexttau) .and. associated(suexttau)) totexttau(:,:,w) = totexttau(:,:,w)+suexttau(:,:,w)
      !          if(associated(totstexttau) .and. associated(sustexttau)) totstexttau(:,:,w) = totstexttau(:,:,w)+sustexttau(:,:,w)
      !          if(associated(totscatau) .and. associated(suscatau)) totscatau(:,:,w) = totscatau(:,:,w)+suscatau(:,:,w)
      !          if(associated(totstscatau) .and. associated(sustscatau)) totstscatau(:,:,w) = totstscatau(:,:,w)+sustscatau(:,:,w)
      !          if(associated(totextt25) .and. associated(suexttau)) totextt25(:,:,w) = totextt25(:,:,w)+suexttau(:,:,w)
      !          if(associated(totscat25) .and. associated(suscatau)) totscat25(:,:,w) = totscat25(:,:,w)+suscatau(:,:,w)
      !          if(associated(totexttfm) .and. associated(suexttau)) totexttfm(:,:,w) = totexttfm(:,:,w)+suexttau(:,:,w)
      !          if(associated(totscatfm) .and. associated(suscatau)) totscatfm(:,:,w) = totscatfm(:,:,w)+suscatau(:,:,w)
      !       end do

      !       do w = 1, size(self%wavelengths_profile)
      !          if(associated(totextcoef) .and. associated(suextcoef)) totextcoef(:,:,:,w) = totextcoef(:,:,:,w)+suextcoef(:,:,:,w)
      !          if(associated(totextcoefrh20) .and. associated(suextcoefrh20)) totextcoefrh20(:,:,:,w) = totextcoefrh20(:,:,:,w)+suextcoefrh20(:,:,:,w)
      !          if(associated(totextcoefrh80) .and. associated(suextcoefrh80)) totextcoefrh80(:,:,:,w) = totextcoefrh80(:,:,:,w)+suextcoefrh80(:,:,:,w)
      !          if(associated(totscacoef) .and. associated(suscacoef)) totscacoef(:,:,:,w) = totscacoef(:,:,:,w)+suscacoef(:,:,:,w)
      !          if(associated(totscacoefrh20) .and. associated(suscacoefrh20)) totscacoefrh20(:,:,:,w) = totscacoefrh20(:,:,:,w)+suscacoefrh20(:,:,:,w)
      !          if(associated(totscacoefrh80) .and. associated(suscacoefrh80)) totscacoefrh80(:,:,:,w) = totscacoefrh80(:,:,:,w)+suscacoefrh80(:,:,:,w)
      !          if(associated(totbckcoef) .and. associated(subckcoef)) totbckcoef(:,:,:,w) = totbckcoef(:,:,:,w)+subckcoef(:,:,:,w)
      !       end do

      !       call MAPL_GetPointer (gex(self%SU%instances(n)%id), pso4, 'PSO4', _RC)
      !       if(associated(pso4tot) .and. associated(pso4)) pso4tot = pso4tot + pso4

      !       call MAPL_GetPointer (gex(self%SU%instances(n)%id), so4smass, 'SO4SMASS', _RC)
      !       if(associated(so4smass)) then
      !          if(associated(pm)       ) pm        = pm        + nifactor*so4smass
      !          if(associated(pm25)     ) pm25      = pm25      + nifactor*so4smass
      !          if(associated(pm_rh35)  ) pm_rh35   = pm_rh35   + 1.33*nifactor*so4smass
      !          if(associated(pm25_rh35)) pm25_rh35 = pm25_rh35 + 1.33*nifactor*so4smass
      !          if(associated(pm_rh50)  ) pm_rh50   = pm_rh50   + 1.51*nifactor*so4smass
      !          if(associated(pm25_rh50)) pm25_rh50 = pm25_rh50 + 1.51*nifactor*so4smass
      !       end if

      !       if(associated(totangstr) .and. associated(suexttau) .and. associated(suangstr)) then
      !          tau1 = tau1 + suexttau(:,:,ind550)*exp(c1*suangstr)
      !          tau2 = tau2 + suexttau(:,:,ind550)*exp(c2*suangstr)
      !       end if
      !    end if
      ! end do


      ! ! Carbonaceous aerosols
      ! do n = 1, size(self%CA%instances)
      !    if ((self%CA%instances(n)%is_active) .and. (index(self%CA%instances(n)%name, 'data') == 0 ) &
      !         .and. (index(self%CA%instances(n)%name, 'CA.bc') > 0)) then

      !       call MAPL_GetPointer (gex(self%CA%instances(n)%id), bcexttau, 'CA.bcEXTTAU', _RC)
      !       call MAPL_GetPointer (gex(self%CA%instances(n)%id), bcstexttau, 'CA.bcSTEXTTAU', _RC)
      !       call MAPL_GetPointer (gex(self%CA%instances(n)%id), bcscatau, 'CA.bcSCATAU', _RC)
      !       call MAPL_GetPointer (gex(self%CA%instances(n)%id), bcstscatau, 'CA.bcSTSCATAU', _RC)
      !       call MAPL_GetPointer (gex(self%CA%instances(n)%id), bcangstr, 'CA.bcANGSTR', _RC)
      !       call MAPL_GetPointer (gex(self%CA%instances(n)%id), bcextcoef, 'CA.bcEXTCOEF', _RC)
      !       call MAPL_GetPointer (gex(self%CA%instances(n)%id), bcextcoefrh20, 'CA.bcEXTCOEFRH20', _RC)
      !       call MAPL_GetPointer (gex(self%CA%instances(n)%id), bcextcoefrh80, 'CA.bcEXTCOEFRH80', _RC)
      !       call MAPL_GetPointer (gex(self%CA%instances(n)%id), bcscacoef, 'CA.bcSCACOEF', _RC)
      !       call MAPL_GetPointer (gex(self%CA%instances(n)%id), bcscacoefrh20, 'CA.bcSCACOEFRH20', _RC)
      !       call MAPL_GetPointer (gex(self%CA%instances(n)%id), bcscacoefrh80, 'CA.bcSCACOEFRH80', _RC)
      !       call MAPL_GetPointer (gex(self%CA%instances(n)%id), bcbckcoef, 'CA.bcBCKCOEF', _RC)

      !       ! Iterate over the wavelengths
      !       do w = 1, size(self%wavelengths_vertint)
      !          if(associated(totexttau) .and. associated(bcexttau)) totexttau(:,:,w) = totexttau(:,:,w)+bcexttau(:,:,w)
      !          if(associated(totstexttau) .and. associated(bcstexttau)) totstexttau(:,:,w) = totstexttau(:,:,w)+bcstexttau(:,:,w)
      !          if(associated(totscatau) .and. associated(bcscatau)) totscatau(:,:,w) = totscatau(:,:,w)+bcscatau(:,:,w)
      !          if(associated(totstscatau) .and. associated(bcstscatau)) totstscatau(:,:,w) = totstscatau(:,:,w)+bcstscatau(:,:,w)
      !          if(associated(totextt25) .and. associated(bcexttau)) totextt25(:,:,w) = totextt25(:,:,w)+bcexttau(:,:,w)
      !          if(associated(totscat25) .and. associated(bcscatau)) totscat25(:,:,w) = totscat25(:,:,w)+bcscatau(:,:,w)
      !          if(associated(totexttfm) .and. associated(bcexttau)) totexttfm(:,:,w) = totexttfm(:,:,w)+bcexttau(:,:,w)
      !          if(associated(totscatfm) .and. associated(bcscatau)) totscatfm(:,:,w) = totscatfm(:,:,w)+bcscatau(:,:,w)
      !       end do

      !       do w = 1, size(self%wavelengths_profile)
      !          if(associated(totextcoef) .and. associated(bcextcoef)) totextcoef(:,:,:,w) = totextcoef(:,:,:,w)+bcextcoef(:,:,:,w)
      !          if(associated(totextcoefrh20) .and. associated(bcextcoefrh20)) totextcoefrh20(:,:,:,w) = totextcoefrh20(:,:,:,w)+bcextcoefrh20(:,:,:,w)
      !          if(associated(totextcoefrh80) .and. associated(bcextcoefrh80)) totextcoefrh80(:,:,:,w) = totextcoefrh80(:,:,:,w)+bcextcoefrh80(:,:,:,w)
      !          if(associated(totscacoef) .and. associated(bcscacoef)) totscacoef(:,:,:,w) = totscacoef(:,:,:,w)+bcscacoef(:,:,:,w)
      !          if(associated(totscacoefrh20) .and. associated(bcscacoefrh20)) totscacoefrh20(:,:,:,w) = totscacoefrh20(:,:,:,w)+bcscacoefrh20(:,:,:,w)
      !          if(associated(totscacoefrh80) .and. associated(bcscacoefrh80)) totscacoefrh80(:,:,:,w) = totscacoefrh80(:,:,:,w)+bcscacoefrh80(:,:,:,w)
      !          if(associated(totbckcoef) .and. associated(bcbckcoef)) totbckcoef(:,:,:,w) = totbckcoef(:,:,:,w)+bcbckcoef(:,:,:,w)
      !       end do

      !       call MAPL_GetPointer (gex(self%CA%instances(n)%id), bcsmass, 'CA.bcSMASS', _RC)
      !       if(associated(pm)        .and. associated(bcsmass)) pm        = pm        + bcsmass
      !       if(associated(pm25)      .and. associated(bcsmass)) pm25      = pm25      + bcsmass
      !       if(associated(pm_rh35)   .and. associated(bcsmass)) pm_rh35   = pm_rh35   + bcsmass
      !       if(associated(pm25_rh35) .and. associated(bcsmass)) pm25_rh35 = pm25_rh35 + bcsmass
      !       if(associated(pm_rh50)   .and. associated(bcsmass)) pm_rh50   = pm_rh50   + bcsmass
      !       if(associated(pm25_rh50) .and. associated(bcsmass)) pm25_rh50 = pm25_rh50 + bcsmass

      !       if(associated(totangstr) .and. associated(bcexttau) .and. associated(bcangstr)) then
      !          tau1 = tau1 + bcexttau(:,:,ind550)*exp(c1*bcangstr)
      !          tau2 = tau2 + bcexttau(:,:,ind550)*exp(c2*bcangstr)
      !       end if

      !    else if ((self%CA%instances(n)%is_active) .and. (index(self%CA%instances(n)%name, 'data') == 0 ) &
      !         .and. (index(self%CA%instances(n)%name, 'CA.oc') > 0)) then
      !       call MAPL_GetPointer (gex(self%CA%instances(n)%id), ocexttau, 'CA.ocEXTTAU', _RC)
      !       call MAPL_GetPointer (gex(self%CA%instances(n)%id), ocstexttau, 'CA.ocSTEXTTAU', _RC)
      !       call MAPL_GetPointer (gex(self%CA%instances(n)%id), ocscatau, 'CA.ocSCATAU', _RC)
      !       call MAPL_GetPointer (gex(self%CA%instances(n)%id), ocstscatau, 'CA.ocSTSCATAU', _RC)
      !       call MAPL_GetPointer (gex(self%CA%instances(n)%id), ocangstr, 'CA.ocANGSTR', _RC)
      !       call MAPL_GetPointer (gex(self%CA%instances(n)%id), ocextcoef, 'CA.ocEXTCOEF', _RC)
      !       call MAPL_GetPointer (gex(self%CA%instances(n)%id), ocextcoefrh20, 'CA.ocEXTCOEFRH20', _RC)
      !       call MAPL_GetPointer (gex(self%CA%instances(n)%id), ocextcoefrh80, 'CA.ocEXTCOEFRH80', _RC)
      !       call MAPL_GetPointer (gex(self%CA%instances(n)%id), ocscacoef, 'CA.ocSCACOEF', _RC)
      !       call MAPL_GetPointer (gex(self%CA%instances(n)%id), ocscacoefrh20, 'CA.ocSCACOEFRH20', _RC)
      !       call MAPL_GetPointer (gex(self%CA%instances(n)%id), ocscacoefrh80, 'CA.ocSCACOEFRH80', _RC)
      !       call MAPL_GetPointer (gex(self%CA%instances(n)%id), ocbckcoef, 'CA.ocBCKCOEF', _RC)

      !       ! Iterate over the wavelengths
      !       do w = 1, size(self%wavelengths_vertint)
      !          if(associated(totexttau) .and. associated(ocexttau)) totexttau(:,:,w) = totexttau(:,:,w)+ocexttau(:,:,w)
      !          if(associated(totstexttau) .and. associated(ocstexttau)) totstexttau(:,:,w) = totstexttau(:,:,w)+ocstexttau(:,:,w)
      !          if(associated(totscatau) .and. associated(ocscatau)) totscatau(:,:,w) = totscatau(:,:,w)+ocscatau(:,:,w)
      !          if(associated(totstscatau) .and. associated(ocstscatau)) totstscatau(:,:,w) = totstscatau(:,:,w)+ocstscatau(:,:,w)
      !          if(associated(totextt25) .and. associated(ocexttau)) totextt25(:,:,w) = totextt25(:,:,w)+ocexttau(:,:,w)
      !          if(associated(totscat25) .and. associated(ocscatau)) totscat25(:,:,w) = totscat25(:,:,w)+ocscatau(:,:,w)
      !          if(associated(totexttfm) .and. associated(ocexttau)) totexttfm(:,:,w) = totexttfm(:,:,w)+ocexttau(:,:,w)
      !          if(associated(totscatfm) .and. associated(ocscatau)) totscatfm(:,:,w) = totscatfm(:,:,w)+ocscatau(:,:,w)
      !       end do

      !       do w = 1, size(self%wavelengths_profile)
      !          if(associated(totextcoef) .and. associated(ocextcoef)) totextcoef(:,:,:,w) = totextcoef(:,:,:,w)+ocextcoef(:,:,:,w)
      !          if(associated(totextcoefrh20) .and. associated(ocextcoefrh20)) totextcoefrh20(:,:,:,w) = totextcoefrh20(:,:,:,w)+ocextcoefrh20(:,:,:,w)
      !          if(associated(totextcoefrh80) .and. associated(ocextcoefrh80)) totextcoefrh80(:,:,:,w) = totextcoefrh80(:,:,:,w)+ocextcoefrh80(:,:,:,w)
      !          if(associated(totscacoef) .and. associated(ocscacoef)) totscacoef(:,:,:,w) = totscacoef(:,:,:,w)+ocscacoef(:,:,:,w)
      !          if(associated(totscacoefrh20) .and. associated(ocscacoefrh20)) totscacoefrh20(:,:,:,w) = totscacoefrh20(:,:,:,w)+ocscacoefrh20(:,:,:,w)
      !          if(associated(totscacoefrh80) .and. associated(ocscacoefrh80)) totscacoefrh80(:,:,:,w) = totscacoefrh80(:,:,:,w)+ocscacoefrh80(:,:,:,w)
      !          if(associated(totbckcoef) .and. associated(ocbckcoef)) totbckcoef(:,:,:,w) = totbckcoef(:,:,:,w)+ocbckcoef(:,:,:,w)
      !       end do

      !       call MAPL_GetPointer (gex(self%CA%instances(n)%id), ocsmass, 'CA.ocSMASS', _RC)
      !       if(associated(pm)        .and. associated(ocsmass)) pm        = pm        + ocsmass
      !       if(associated(pm25)      .and. associated(ocsmass)) pm25      = pm25      + ocsmass
      !       if(associated(pm_rh35)   .and. associated(ocsmass)) pm_rh35   = pm_rh35   + 1.16*ocsmass  ! needs to be revisited: OCpho + 1.16 OCphi
      !       if(associated(pm25_rh35) .and. associated(ocsmass)) pm25_rh35 = pm25_rh35 + 1.16*ocsmass  !
      !       if(associated(pm_rh50)   .and. associated(ocsmass)) pm_rh50   = pm_rh50   + 1.24*ocsmass  ! needs to be revisited: OCpho + 1.24 OCphi
      !       if(associated(pm25_rh50) .and. associated(ocsmass)) pm25_rh50 = pm25_rh50 + 1.24*ocsmass  !

      !       if(associated(totangstr) .and. associated(ocexttau) .and. associated(ocangstr)) then
      !          tau1 = tau1 + ocexttau(:,:,ind550)*exp(c1*ocangstr)
      !          tau2 = tau2 + ocexttau(:,:,ind550)*exp(c2*ocangstr)
      !       end if

      !    else if ((self%CA%instances(n)%is_active) .and. (index(self%CA%instances(n)%name, 'data') == 0 ) &
      !         .and. (index(self%CA%instances(n)%name, 'CA.br') > 0)) then
      !       call MAPL_GetPointer (gex(self%CA%instances(n)%id), brexttau, 'CA.brEXTTAU', _RC)
      !       call MAPL_GetPointer (gex(self%CA%instances(n)%id), brstexttau, 'CA.brSTEXTTAU', _RC)
      !       call MAPL_GetPointer (gex(self%CA%instances(n)%id), brscatau, 'CA.brSCATAU', _RC)
      !       call MAPL_GetPointer (gex(self%CA%instances(n)%id), brstscatau, 'CA.brSTSCATAU', _RC)
      !       call MAPL_GetPointer (gex(self%CA%instances(n)%id), brangstr, 'CA.brANGSTR', _RC)
      !       call MAPL_GetPointer (gex(self%CA%instances(n)%id), brextcoef, 'CA.brEXTCOEF', _RC)
      !       call MAPL_GetPointer (gex(self%CA%instances(n)%id), brextcoefrh20, 'CA.brEXTCOEFRH20', _RC)
      !       call MAPL_GetPointer (gex(self%CA%instances(n)%id), brextcoefrh80, 'CA.brEXTCOEFRH80', _RC)
      !       call MAPL_GetPointer (gex(self%CA%instances(n)%id), brscacoef, 'CA.brSCACOEF', _RC)
      !       call MAPL_GetPointer (gex(self%CA%instances(n)%id), brscacoefrh20, 'CA.brSCACOEFRH20', _RC)
      !       call MAPL_GetPointer (gex(self%CA%instances(n)%id), brscacoefrh80, 'CA.brSCACOEFRH80', _RC)
      !       call MAPL_GetPointer (gex(self%CA%instances(n)%id), brbckcoef, 'CA.brBCKCOEF', _RC)

      !       ! Iterate over the wavelengths
      !       do w = 1, size(self%wavelengths_vertint)
      !          if(associated(totexttau) .and. associated(brexttau)) totexttau(:,:,w) = totexttau(:,:,w)+brexttau(:,:,w)
      !          if(associated(totstexttau) .and. associated(brstexttau)) totstexttau(:,:,w) = totstexttau(:,:,w)+brstexttau(:,:,w)
      !          if(associated(totscatau) .and. associated(brscatau)) totscatau(:,:,w) = totscatau(:,:,w)+brscatau(:,:,w)
      !          if(associated(totstscatau) .and. associated(brstscatau)) totstscatau(:,:,w) = totstscatau(:,:,w)+brstscatau(:,:,w)
      !          if(associated(totextt25) .and. associated(brexttau)) totextt25(:,:,w) = totextt25(:,:,w)+brexttau(:,:,w)
      !          if(associated(totscat25) .and. associated(brscatau)) totscat25(:,:,w) = totscat25(:,:,w)+brscatau(:,:,w)
      !          if(associated(totexttfm) .and. associated(brexttau)) totexttfm(:,:,w) = totexttfm(:,:,w)+brexttau(:,:,w)
      !          if(associated(totscatfm) .and. associated(brscatau)) totscatfm(:,:,w) = totscatfm(:,:,w)+brscatau(:,:,w)
      !       end do

      !       do w = 1, size(self%wavelengths_profile)
      !          if(associated(totextcoef) .and. associated(brextcoef)) totextcoef(:,:,:,w) = totextcoef(:,:,:,w)+brextcoef(:,:,:,w)
      !          if(associated(totextcoefrh20) .and. associated(brextcoefrh20)) totextcoefrh20(:,:,:,w) = totextcoefrh20(:,:,:,w)+brextcoefrh20(:,:,:,w)
      !          if(associated(totextcoefrh80) .and. associated(brextcoefrh80)) totextcoefrh80(:,:,:,w) = totextcoefrh80(:,:,:,w)+brextcoefrh80(:,:,:,w)
      !          if(associated(totscacoef) .and. associated(brscacoef)) totscacoef(:,:,:,w) = totscacoef(:,:,:,w)+brscacoef(:,:,:,w)
      !          if(associated(totscacoefrh20) .and. associated(brscacoefrh20)) totscacoefrh20(:,:,:,w) = totscacoefrh20(:,:,:,w)+brscacoefrh20(:,:,:,w)
      !          if(associated(totscacoefrh80) .and. associated(brscacoefrh80)) totscacoefrh80(:,:,:,w) = totscacoefrh80(:,:,:,w)+brscacoefrh80(:,:,:,w)
      !          if(associated(totbckcoef) .and. associated(brbckcoef)) totbckcoef(:,:,:,w) = totbckcoef(:,:,:,w)+brbckcoef(:,:,:,w)
      !       end do

      !       call MAPL_GetPointer (gex(self%CA%instances(n)%id), brsmass, 'CA.brSMASS', _RC)
      !       if(associated(pm)        .and. associated(brsmass)) pm        = pm        + brsmass
      !       if(associated(pm25)      .and. associated(brsmass)) pm25      = pm25      + brsmass
      !       if(associated(pm_rh35)   .and. associated(brsmass)) pm_rh35   = pm_rh35   + 1.16*brsmass  ! needs to be revisited: OCpho + 1.16 OCphi
      !       if(associated(pm25_rh35) .and. associated(brsmass)) pm25_rh35 = pm25_rh35 + 1.16*brsmass  !
      !       if(associated(pm_rh50)   .and. associated(brsmass)) pm_rh50   = pm_rh50   + 1.24*brsmass  ! needs to be revisited: OCpho + 1.24 OCphi
      !       if(associated(pm25_rh50) .and. associated(brsmass)) pm25_rh50 = pm25_rh50 + 1.24*brsmass  !

      !       if(associated(totangstr) .and. associated(brexttau) .and. associated(brangstr)) then
      !          tau1 = tau1 + brexttau(:,:,ind550)*exp(c1*brangstr)
      !          tau2 = tau2 + brexttau(:,:,ind550)*exp(c2*brangstr)
      !       end if
      !    end if
      ! end do

      ! Finish calculating totangstr
      if(associated(totangstr)) then
         totangstr = log(tau1/tau2)/c3
      end if

      ! pchakrab - THE FOLLOWING CODE IS DIFFICULT TO TEST
      ! For testing, we activate all exports, but 532nm wavelengths are not provided
      ! Commenting out, for now
      ! ! Calculate the total (molecular + aer) single scattering attenuated backscater coef from the TOA
      ! if(associated(totabcktoa).or.associated(totabcksfc)) then
      !    _ASSERT(associated(totextcoef), "Cannot produce TOTABCKTOA/TOTABCKSFC without totextcoef")
      !    _ASSERT(associated(totbckcoef), "Cannot produce TOTABCKTOA/TOTABCKSFC without totbckcoef")
      !    ind532 = 0
      !    do w = 1, size(self%wavelengths_profile) ! find index for 532nm to compute TBA
      !       if ((self%wavelengths_profile(w)*1.e-9 .ge. 5.31e-7) .and. &
      !            (self%wavelengths_profile(w)*1.e-9 .le. 5.33e-7)) then
      !          ind532 = w
      !          exit
      !       end if
      !    end do
      !    _ASSERT(ind532 /= 0, "Cannot produce TOTBCKCOEF variables without 532nm wavelength")

      !    ! Pressure at layer edges (ple shape (im,jm, km+1) on the edge
      !    i1 = lbound(ple, 1); i2 = ubound(ple, 1)
      !    j1 = lbound(ple, 2); j2 = ubound(ple, 2)
      !    km = ubound(ple, 3) ! km =72 index starts at 0
      !    ! Pressure for each layer
      !    allocate(P(i1:i2,j1:j2,km), __STAT__)
      !    do k = 1, km
      !       P(:,:,k) = 0.5 * (ple(:,:,k-1) + ple(:,:,k))   ! in Pa
      !    enddo

      !    !molecular backscattering cross section for each layer at 532nm: Cair  * P(Pa) / T(K)
      !    !Cair = 4.51944e-9 at 532nm # unit K Pa-1 m-1 sr-1 http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19960051003.pdf
      !    allocate(backscat_mol(i1:i2,j1:j2,km), __STAT__)
      !    backscat_mol = (5.45e-32/1.380648e-23) * (532./550.)**(-4.0)  * P / T
      !    ! tau mol for each layer
      !    allocate(tau_mol_layer(i1:i2,j1:j2,km), delz(i1:i2,j1:j2,km),__STAT__)
      !    delz  = delp / (MAPL_GRAV * airdens)
      !    tau_mol_layer = backscat_mol * 8.* pi /3. * delz

      !    ! tau aer for each layer
      !    allocate(tau_aer_layer(i1:i2,j1:j2,km), __STAT__)
      !    tau_aer_layer = totextcoef(:,:,:,ind532) * delz

      !    allocate(tau_aer(i1:i2,j1:j2), __STAT__)
      !    allocate(tau_mol(i1:i2,j1:j2), __STAT__)

      !    ! TOTAL ABCK TOA
      !    ! top layer
      !    totabcktoa(:,:,1) = (totbckcoef(:,:,1,ind532) + backscat_mol(:,:,1)) * exp(-tau_aer_layer(:,:,1)) * exp(-tau_mol_layer(:,:,1))
      !    ! layer 2 to the layer at the surface(km)
      !    do k = 2, km
      !       tau_aer = 0.
      !       tau_mol = 0. ! for each layer
      !       do kk = 1, k
      !          tau_aer = tau_aer + tau_aer_layer(:,:,kk)
      !          tau_mol = tau_mol + tau_mol_layer(:,:,kk)
      !       enddo
      !       tau_aer = tau_aer + 0.5 *  tau_aer_layer(:,:,k)
      !       tau_mol = tau_mol + 0.5 *  tau_mol_layer(:,:,k)
      !       totabcktoa(:,:,k) = (totbckcoef(:,:,k,ind532) + backscat_mol(:,:,k)) * exp(-tau_aer) * exp(-tau_mol)
      !    enddo

      !    ! TOTAL ABCK SFC
      !    ! bottom layer
      !    totabcksfc(:,:,km) = (totbckcoef(:,:,km,ind532) + backscat_mol(:,:,km)) * exp(-tau_aer_layer(:,:,km)) * exp(-tau_mol_layer(:,:,km))
      !    ! layer 2nd from the surface to the top of the atmoshere (km)
      !    do k = km-1, 1, -1
      !       tau_aer = 0.
      !       tau_mol = 0. ! for each layer
      !       do kk = km, k+1, -1
      !          tau_aer = tau_aer + tau_aer_layer(:,:,kk)
      !          tau_mol = tau_mol + tau_mol_layer(:,:,kk)
      !       enddo
      !       tau_aer = tau_aer + 0.5 *  tau_aer_layer(:,:,k)
      !       tau_mol = tau_mol + 0.5 *  tau_mol_layer(:,:,k)
      !       totabcksfc(:,:,k) = (totbckcoef(:,:,k,ind532) + backscat_mol(:,:,k)) * exp(-tau_aer) * exp(-tau_mol)
      !    enddo

      ! endif ! end of total attenuated backscatter coef calculation

      _RETURN(_SUCCESS)
      _UNUSED_DUMMY(clock)

   end subroutine Run2

   subroutine setup_constituents_(self, hconfig, rc)
      type (GOCART_State), pointer, intent(in) :: self
      type(ESMF_HConfig), intent(in) :: hconfig
      integer, intent(out) :: rc

      type(ESMF_HConfig) :: active_cfg, passive_cfg
      logical :: has_section
      character(len=*), parameter :: ACTIVE_INSTANCES_SECTION = "ACTIVE_INSTANCES"
      character(len=*), parameter :: PASSIVE_INSTANCES_SECTION = "PASSIVE_INSTANCES"
      integer :: status

      self%DU%name = "DU"
      self%SS%name = "SS"
      self%SU%name = "SU"
      self%CA%name = "CA"
      self%NI%name = "NI"

      ! Active instances
      has_section = ESMF_HConfigIsDefined(hconfig, keyString=ACTIVE_INSTANCES_SECTION, _RC)
      _ASSERT(has_section, ACTIVE_INSTANCES_SECTION // " not found")
      active_cfg = ESMF_HConfigCreateAt(hconfig, keyString=ACTIVE_INSTANCES_SECTION, _RC)
      ! Passive instances - checked for each aerosol
      has_section = ESMF_HConfigIsDefined(hconfig, keyString=PASSIVE_INSTANCES_SECTION, _RC)
      _ASSERT(has_section, PASSIVE_INSTANCES_SECTION // " not found")
      passive_cfg = ESMF_HConfigCreateAt(hconfig, keyString=PASSIVE_INSTANCES_SECTION, _RC)

      call setup_instances__(self%DU, active_cfg, passive_cfg, _RC)
      call setup_instances__(self%SS, active_cfg, passive_cfg, _RC)
      call setup_instances__(self%SU, active_cfg, passive_cfg, _RC)
      call setup_instances__(self%CA, active_cfg, passive_cfg, _RC)
      call setup_instances__(self%NI, active_cfg, passive_cfg, _RC)

      _RETURN(_SUCCESS)
   end subroutine setup_constituents_

   subroutine setup_instances__(species, active_cfg, passive_cfg, rc)
      type(Constituent), intent(inout)  :: species
      type(ESMF_HConfig), intent(in) :: active_cfg, passive_cfg
      integer, intent(out) :: rc

      character(len=:), allocatable :: active_instances(:), passive_instances(:)
      integer :: iter, n_active, n_passive, status

      active_instances = ESMF_HConfigAsStringSeq(active_cfg, keyString=species%name, stringLen=ESMF_MAXSTR, _RC)
      n_active = size(active_instances)
      passive_instances = ESMF_HConfigAsStringSeq(passive_cfg, keyString=species%name, stringLen=ESMF_MAXSTR, _RC)
      n_passive = size(passive_instances)
      allocate(species%instances(n_active+n_passive), _STAT)

      ! IMPORTANT: Active instances must be created first! This ordering is necessary for
      ! filing the AERO states that are passed to radiation.
      ! This is achieved by arranging the names of the active instances first.

      ! Fill the instances list with active instances first
      do iter = 1, n_active
         species%instances(iter)%name = trim(active_instances(iter))
         species%instances(iter)%is_active = .true.
      end do
      species%n_active = n_active

      ! Now fill instances list with passive instances
      do iter = 1, n_passive
         species%instances(iter+n_active)%name = trim(passive_instances(iter))
         species%instances(iter+n_active)%is_active = .false.
      end do

      _RETURN(_SUCCESS)
   end subroutine setup_instances__

   ! Creates GOCART2G children. Active instances must be created first. If
   ! additional GOCART2G children are added, this subroutine will need to be updated
   subroutine create_instances_(self, gc, rc)
      type (GOCART_State), pointer, intent(in) :: self
      type (ESMF_GridComp), intent(inout) :: gc
      integer, intent(out) :: rc

      integer :: status

      call add_children__(gc, self%DU, DU2G_SetServices, _RC)
      call add_children__(gc, self%SS, SS2G_SetServices, _RC)
      ! call add_children__(gc, self%SU, SU2G_SetServices, _RC)
      ! call add_children__(gc, self%CA, CA2G_SetServices, _RC)
      ! call add_children__(gc, self%NI, NI2G_SetServices, _RC)

      _RETURN(_SUCCESS)
   end subroutine create_instances_

   subroutine add_children__(gc, species, setservices, rc)
      use mapl3g_UserSetServices, only: user_setservices
      type (ESMF_GridComp), intent(inout) :: gc
      type(Constituent), intent(inout) :: species
      external :: setservices
      integer, intent(out) :: rc

      integer :: iter, status
      character(len=:), allocatable :: child_name, config_file
      type(ESMF_HConfig) :: hconfig

      do iter = 1, size(species%instances)
         child_name = species%instances(iter)%name
         config_file = species%name // "_instance_" // child_name // ".yaml"
         hconfig = ESMF_HConfigCreate(filename=config_file, _RC)
         call MAPL_GridCompAddChild(gc, child_name, user_setservices(setservices), hconfig, _RC)
         call ESMF_HConfigDestroy(hconfig, _RC)
      end do

      _RETURN(_SUCCESS)
   end subroutine add_children__

   subroutine serialize_bundle(state, rc)
      !ARGUMENTS:
      type(ESMF_State) :: state
      integer, intent(out) :: rc

      !Local
      character(len=ESMF_MAXSTR), allocatable :: itemList(:)
      type(ESMF_State) :: child_state
      type(ESMF_StateItem_Flag), allocatable :: itemTypes(:)
      type(ESMF_FieldBundle) :: bundle
      type(ESMF_Grid) :: grid
      type(ESMF_Field) :: field, serializedField
      type(ESMF_Info) :: info

      character(len=ESMF_MAXSTR) :: binIndexstr
      character(len=ESMF_MAXSTR), allocatable :: aeroName(:)

      real, pointer, dimension(:,:,:,:) :: orig_ptr
      real, pointer, dimension(:,:,:) :: ptr3d

      integer :: b, i, j, n, rank, status

      !Description: Callback for AERO_RAD state used in GAAS module to provide a
      !                 serialized ESMF_Bundle of aerosol fields.

      ! Get list of child states within state and add to aeroList
      ! Remember, AERO_RAD contains its children's AERO_RAD states
      call ESMF_StateGet(state, itemCount=n, _RC)
      allocate(itemList(n), _STAT)
      allocate(itemTypes(n), _STAT)
      call ESMF_StateGet(state, itemNameList=itemList, itemTypeList=itemTypes, _RC)

      ! Create empty ESMF_FieldBundle to add Children's aerosol fields to
      bundle = ESMF_FieldBundleCreate(name="serialized_aerosolBundle", _RC)
      call ESMF_StateAdd(state, [bundle], _RC)

      do i = 1, n
         if (itemTypes(i) /= ESMF_STATEITEM_STATE) cycle ! exclude non-states
         call ESMF_StateGet(state, trim(itemList(i)), child_state, _RC)
         call ESMF_InfoGetFromHost(state, info, _RC)
         call ESMF_InfoGet(info, key="internal_variable_name", values=aeroName, _RC)
         do b = 1, size(aeroName)
            call ESMF_StateGet(child_state, trim(aeroName(b)), field, _RC)
            call ESMF_FieldGet(field, rank=rank, _RC)
            select case(rank)
            case(3)
               call MAPL_FieldBundleAdd(bundle, [field], _RC)
            case(4) ! serialize 4d variables to mulitple 3d variables
               call ESMF_FieldGet(field, grid=grid, _RC)
               call MAPL_StateGetPointer(child_state, itemName=trim(aeroName(b)), farrayPtr=orig_ptr, _RC)
               do j = 1, size(orig_ptr, 4)
                  write (binIndexstr, "(I0.3)") j
                  ptr3d => orig_ptr(:,:,:,j)
                  ! pchakrab: TODO, we are sharing data here
                  serializedField = ESMF_FieldCreate( &
                       grid=grid, &
                       datacopyFlag=ESMF_DATACOPY_REFERENCE, &
                       farrayPtr=ptr3d, &
                       name=trim(aeroName(b))//trim(binIndexstr), _RC)
                  ! probably need to add a flag to allow for adding multilple fields of the same name
                  call MAPL_FieldBundleAdd(bundle, [serializedField], _RC)
               end do ! do j
            case default
               _FAIL("rank not supported")
            end select ! select case
         end do ! do b
         deallocate (aeroName, _STAT)
      end do ! do i

      _RETURN(_SUCCESS)
   end subroutine serialize_bundle

   subroutine run_aerosol_optics (state, rc)

      !ARGUMENTS:
      type(ESMF_State) :: state
      integer, intent(out) :: rc

      !Local
      real, dimension(:,:,:), pointer :: ple
      real, dimension(:,:,:), pointer :: rh
      real, dimension(:,:,:), pointer  :: var

      character(len=:), allocatable :: fld_name

      real(kind=8), dimension(:,:,:), pointer :: ext_, ssa_, asy_    ! (lon:,lat:,lev:)
      real(kind=8), dimension(:,:,:), allocatable :: ext,  ssa,  asy ! (lon:,lat:,lev:)

      integer :: i, n, b, j
      integer :: i1, j1, i2, j2, km
      integer :: band
      integer, parameter :: n_bands = 1

      character(len=ESMF_MAXSTR), allocatable :: itemList(:), aeroList(:)
      type(ESMF_State) :: child_state
      real, pointer, dimension(:,:,:) :: as_ptr_3d

      type(ESMF_StateItem_Flag), allocatable :: itemTypes(:)
      type(ESMF_Info) :: info, child_info
      integer :: status

      ! Description: Used in Radiation gridded components to provide aerosol properties

      call ESMF_InfoGetFromHost(state, info, _RC)

      ! Radiation band
      call ESMF_InfoGet(info, key="band_for_aerosol_optics", value=band, _RC)

      ! Relative humidity
      call ESMF_InfoGet(info, key="relative_humidity_for_aerosol_optics", value=fld_name, _RC)
      call MAPL_StateGetPointer(state, itemName=fld_name, farrayPtr=rh, _RC)

      ! Pressure at layer edges
      call ESMF_InfoGet(info, key="air_pressure_for_aerosol_optics", value=fld_name, _RC)
      call MAPL_StateGetPointer(state, itemName=fld_name, farrayPtr=ple, _RC)

      ! TODO: pchakrab - CAREFUL! In MAPL3 land, PLE is (:, :, 1:km+1), instead of (:, :, 0:km)
      i1 = lbound(ple, 1); i2 = ubound(ple, 1)
      j1 = lbound(ple, 2); j2 = ubound(ple, 2)
      km = ubound(ple, 3)

      allocate( &
           ext(i1:i2,j1:j2,km),  &
           ssa(i1:i2,j1:j2,km),  &
           asy(i1:i2,j1:j2,km), _STAT)

      ! Get list of child states within state and add to aeroList
      call ESMF_StateGet(state, itemCount=n, _RC)
      allocate(itemList(n), _STAT)
      allocate(itemTypes(n), _STAT)
      call ESMF_StateGet(state, itemNameList=itemList, itemTypeList=itemTypes, _RC)

      b=0
      do i = 1, n
         if (itemTypes(i) == ESMF_STATEITEM_STATE) then
            b = b + 1
         end if
      end do

      allocate (aeroList(b), _STAT)

      j = 1
      do i = 1, n
         if (itemTypes(i) == ESMF_STATEITEM_STATE) then
            aeroList(j) = trim(itemList(i))
            j = j + 1
         end if
      end do

      ext = 0.0d0
      ssa = 0.0d0
      asy = 0.0d0

      ! Get aerosol optic properties from children
      do i = 1, size(aeroList)
         call ESMF_StateGet(state, trim(aeroList(i)), child_state, _RC)

         call ESMF_InfoGetFromHost(child_state, child_info, _RC)

         ! set RH in child's aero state
         call ESMF_InfoGet(child_info, key="relative_humidity_for_aerosol_optics", value=fld_name, _RC)

         if (fld_name /= "") then
            call MAPL_StateGetPointer(child_state, itemName=fld_name, farrayPtr=as_ptr_3d, _RC)
            as_ptr_3d = rh
         end if

         ! set PLE in child's aero state
         call ESMF_InfoGet(child_info, key="air_pressure_for_aerosol_optics", value=fld_name, _RC)

         if (fld_name /= "") then
            call MAPL_StateGetPointer(child_state, itemName=fld_name, farrayPtr=as_ptr_3d, _RC)
            as_ptr_3d = ple
         end if

         ! set band in child's aero state
         call ESMF_InfoSet(child_info, key="band_for_aerosol_optics", value=band, _RC)

         ! execute the aerosol optics method
         call ESMF_MethodExecute(child_state, label="aerosol_optics", _RC)

         ! Retrieve extinction from each child
         call ESMF_InfoGet(child_info, key="extinction_in_air_due_to_ambient_aerosol", value=fld_name, _RC)
         if (fld_name /= "") then
            call MAPL_StateGetPointer(child_state, itemName=fld_name, farrayPtr=ext_, _RC)
         end if

         ! Retrieve scattering extinction from each child
         call ESMF_InfoGet(child_info, key="single_scattering_albedo_of_ambient_aerosol", value=fld_name, _RC)
         if (fld_name /= "") then
            call MAPL_StateGetPointer(child_state, itemName=fld_name, farrayPtr=ssa_, _RC)
         end if

         ! Retrieve asymetry parameter multiplied by scatering extiction from each child
         call ESMF_InfoGet(child_info, key="asymmetry_parameter_of_ambient_aerosol", value=fld_name, _RC)
         if (fld_name /= "") then
            call MAPL_StateGetPointer(child_state, itemName=fld_name, farrayPtr=asy_, _RC)
         end if

         ! Sum aerosol optic properties from each child
         ext = ext + ext_
         ssa = ssa + ssa_
         asy = asy + asy_
      end do

      call ESMF_InfoGetFromHost(state, info, _RC)

      ! Set ext, ssa, asy to equal the sum of ext, ssa, asy from the children. This is what is passed to radiation.
      call ESMF_InfoGet(info, key="extinction_in_air_due_to_ambient_aerosol", value=fld_name, _RC)
      if (fld_name /= "") then
         call MAPL_StateGetPointer(state, itemName=fld_name, farrayPtr=var, _RC)
         var = ext(:,:,:)
      end if

      call ESMF_InfoGet(info, key="single_scattering_albedo_of_ambient_aerosol", value=fld_name, _RC)
      if (fld_name /= "") then
         call MAPL_StateGetPointer(state, itemName=fld_name, farrayPtr=var, _RC)
         var = ssa(:,:,:)
      end if

      call ESMF_InfoGet(info, key="asymmetry_parameter_of_ambient_aerosol", value=fld_name, _RC)
      if (fld_name /= "") then
         call MAPL_StateGetPointer(state, itemName=fld_name, farrayPtr=var, _RC)
         var = asy(:,:,:)
      end if

      deallocate(ext, ssa, asy, __STAT__)

      _RETURN(_SUCCESS)

   end subroutine run_aerosol_optics

   subroutine aerosol_activation_properties(state, rc)

      ! Arguments
      type(ESMF_State) :: state
      integer, intent(out) :: rc

      ! Local
      character(len=ESMF_MAXSTR) :: mode                ! mode name
      character(len=ESMF_MAXSTR) :: mode_               ! lowercase mode name

      type(ESMF_State) :: child_state
      type(ESMF_Info) :: info

      real, dimension(:,:,:), pointer :: ple            ! pressure at the edges of model layers
      real, dimension(:,:,:), pointer :: temperature    ! air temperature
      real, dimension(:,:), pointer :: f_land           ! fraction of land type in a grid cell

      real, dimension(:,:,:), allocatable :: q          ! aerosol mass mixing ratio
      real, dimension(:,:,:,:), pointer :: ptr_4d       ! aerosol mass mixing ratio (temporary)
      real, dimension(:,:,:), pointer :: ptr_3d         ! aerosol mass mixing ratio (temporary)

      real, dimension(:,:,:), pointer :: num            ! number concentration of aerosol particles
      real, dimension(:,:,:), pointer :: diameter       ! dry size of aerosol
      real, dimension(:,:,:), pointer :: sigma          ! width of aerosol mode
      real, dimension(:,:,:), pointer :: density        ! density of aerosol
      real, dimension(:,:,:), pointer :: hygroscopicity ! hygroscopicity of aerosol
      real, dimension(:,:,:), pointer :: f_dust         ! fraction of dust aerosol
      real, dimension(:,:,:), pointer :: f_soot         ! fraction of soot aerosol
      real, dimension(:,:,:), pointer :: f_organic      ! fraction of organic aerosol

      real :: max_clean          ! max mixing ratio before considered polluted
      real :: ccn_tuning         ! tunes conversion factors for sulfate

      character(len=:), allocatable :: fld_name

      integer :: i2, j2, km
      integer :: b, i, j, n, aerosol_bin
      integer :: varNameLen

      character(len=ESMF_MAXSTR), allocatable :: itemList(:), aeroList(:)
      type(ESMF_StateItem_Flag), allocatable :: itemTypes(:)

      ! auxilliary parameters
      real, parameter :: densSO4 = 1700.0
      real, parameter :: densORG = 1600.0
      real, parameter :: densSS = 2200.0
      real, parameter :: densDU = 1700.0
      real, parameter :: densBC = 1600.0
      real, parameter :: densOC =  900.0
      real, parameter :: densBR =  900.0

      real, parameter :: k_SO4 = 0.65
      real, parameter :: k_ORG = 0.20
      real, parameter :: k_SS = 1.28
      real, parameter :: k_DU = 0.0001
      real, parameter :: k_BC = 0.0001
      real, parameter :: k_OC = 0.0001
      real, parameter :: k_BR = 0.0001

      integer, parameter :: UNKNOWN_AEROSOL_MODE = 2015
      integer :: status

      ! Get list of child states within state and add to aeroList
      call ESMF_StateGet(state, itemCount=n, _RC)
      allocate(itemList(n), _STAT)
      allocate(itemTypes(n), _STAT)
      call ESMF_StateGet(state, itemNameList=itemList, itemTypeList=itemTypes, _RC)

      b=0
      do i = 1, n
         if ((itemTypes(i) == ESMF_StateItem_State) .and. (trim(itemList(i)(1:2)) /= "NI")) then
            b = b + 1
         end if
      end do

      allocate(aeroList(b), _STAT)

      j = 1
      do i = 1, n
         if ((itemTypes(i) == ESMF_StateItem_State) .and. (trim(itemList(i)(1:2)) /= "NI")) then
            aeroList(j) = trim(itemList(i))
            j = j + 1
         end if
      end do

      call ESMF_InfoGetFromHost(state, info, _RC)

      ! Aerosol mode
      call ESMF_InfoGet(info, key="aerosol_mode", value=mode, _RC)

      ! Land fraction
      call ESMF_InfoGet(info, key="fraction_of_land_type", value=fld_name, _RC)
      call MAPL_StateGetPointer(state, itemName=fld_name, farrayPtr=f_land, _RC)

      ! Pressure at layer edges
      call ESMF_InfoGet(info, key="air_pressure_for_aerosol_optics", value=fld_name, _RC)
      call MAPL_StateGetPointer(state, itemName=fld_name, farrayPtr=ple, _RC)

      ! Temperature
      call ESMF_InfoGet(info, key="air_temperature", value=fld_name, _RC)
      call MAPL_StateGetPointer(state, itemName=fld_name, farrayPtr=temperature, _RC)

      i2 = ubound(temperature, 1)
      j2 = ubound(temperature, 2)
      km = ubound(temperature, 3)

      ! Activation activation properties
      call ESMF_InfoGet(info, key="aerosol_number_concentration", value=fld_name, _RC)
      call MAPL_StateGetPointer(state, itemName=fld_name, farrayPtr=num, _RC)

      call ESMF_InfoGet(info, key="aerosol_dry_size", value=fld_name, _RC)
      call MAPL_StateGetPointer(state, itemName=fld_name, farrayPtr=diameter, _RC)

      call ESMF_InfoGet(info, key="width_of_aerosol_mode", value=fld_name, _RC)
      call MAPL_StateGetPointer(state, itemName=fld_name, farrayPtr=sigma, _RC)

      call ESMF_InfoGet(info, key="aerosol_density", value=fld_name, _RC)
      call MAPL_StateGetPointer(state, itemName=fld_name, farrayPtr=density, _RC)

      call ESMF_InfoGet(info, key="aerosol_hygroscopicity", value=fld_name, _RC)
      call MAPL_StateGetPointer(state, itemName=fld_name, farrayPtr=hygroscopicity, _RC)

      call ESMF_InfoGet(info, key="fraction_of_dust_aerosol", value=fld_name, _RC)
      call MAPL_StateGetPointer(state, itemName=fld_name, farrayPtr=f_dust, _RC)

      call ESMF_InfoGet(info, key="fraction_of_soot_aerosol", value=fld_name, _RC)
      call MAPL_StateGetPointer(state, itemName=fld_name, farrayPtr=f_soot, _RC)

      call ESMF_InfoGet(info, key="fraction_of_organic_aerosol", value=fld_name, _RC)
      call MAPL_StateGetPointer(state, itemName=fld_name, farrayPtr=f_organic, _RC)

      ! Sea salt scaling fctor
      call ESMF_InfoGet(info, key="max_q_clean", value=max_clean, _RC)
      call ESMF_InfoGet(info, key="ccn_tuning", value=ccn_tuning, _RC)

      ! Aerosol mass mixing ratios
      mode_ = trim(mode)
      mode_ = ESMF_UtilStringLowerCase(mode_, _RC)

      allocate(q(i2,j2,km), _STAT)
      q = 0.0

      if (index(mode_, "du00") > 0) then ! Dust
         ! dust is mapped one-to-one
         do i = 1, size(aeroList)
            if (index(aeroList(i), "DU") > 0) then
               read (mode_(3:len(mode_)),*) aerosol_bin
               call ESMF_StateGet(state, trim(aeroList(i)), child_state, _RC)
               call MAPL_StateGetPointer(child_state, itemName="DU", farrayPtr=ptr_4d, _RC)
               q = q + ptr_4d(:,:,:,aerosol_bin)
               ptr_3d => ptr_4d(:,:,:,aerosol_bin)

               hygroscopicity = k_DU
               density = densDU
            end if
         end do

      else if (index(mode_, "ss00") > 0) then ! Sea Salt
         ! compute the total mass mixing ratio and impose a tri-modal size distribution
         do i = 1, size(aeroList)
            if (index(aeroList(i), "SS") > 0) then
               call ESMF_StateGet(state, trim(aeroList(i)), child_state, _RC)
               call MAPL_StateGetPointer(child_state, itemName="SS", farrayPtr=ptr_4d, _RC)
               do j = 1, ubound(ptr_4d, 4)
                  q = q + ptr_4d(:,:,:,j)
                  ptr_3d => ptr_4d(:,:,:,j)
               end do

               hygroscopicity = k_SS
               density = densSS
            end if
         end do

      else if (index(mode_, "sulforg") > 0) then ! Sulfate
         hygroscopicity = 0.0
         density = 0.0

         do i = 1, size(aeroList)
            if (index(aeroList(i), "SU") > 0) then
               call ESMF_StateGet(state, trim(aeroList(i)), child_state, _RC)
               call MAPL_StateGetPointer(child_state, itemName="SO4", farrayPtr=ptr_3d, _RC)
               q = q + ptr_3d
               hygroscopicity = k_SO4 * ptr_3d + hygroscopicity
               density = densSO4 * ptr_3d + density
            end if

            if (index(aeroList(i), "CA.oc") > 0) then
               call ESMF_StateGet(state, trim(aeroList(i)), child_state, _RC)
               varNameLen = len_trim(aeroList(i))
               ! the "5" refers to "_AERO", which we want to remove
               ! to get the CA component name (e.g. CA.oc, or CA.oc.data)
               varNameLen = varNameLen - 5
               call MAPL_StateGetPointer(child_state, itemName=aeroList(i)(1:varNameLen)//"philic", farrayPtr=ptr_3d, _RC)
               q = q + ptr_3d
               hygroscopicity = k_ORG * ptr_3d + hygroscopicity
               density = densORG * ptr_3d + density
            end if
         end do

         where (q > 2.0e-12 .and. hygroscopicity > tiny(0.0))
            hygroscopicity = hygroscopicity / q
            hygroscopicity = max(0.001, hygroscopicity)

            density = density / q
            density = min(max(density, densORG), densSO4)
         elsewhere
            hygroscopicity = k_SO4
            density = densSO4
         end where

      else if (index(mode_, "bcphilic") > 0) then ! Black Carbon
         do i = 1, size(aeroList)
            if (index(aeroList(i), "CA.bc") > 0) then
               call ESMF_StateGet(state, trim(aeroList(i)), child_state, _RC)
               varNameLen = len_trim(aeroList(i))
               ! the "5" refers to "_AERO", which we want to remove
               ! to get the CA component name (e.g. CA.bc, or CA.bc.data)
               varNameLen = varNameLen - 5
               call MAPL_StateGetPointer(child_state, itemName=aeroList(i)(1:varNameLen)//"philic", farrayPtr=ptr_3d, _RC)
               q = q + ptr_3d
               hygroscopicity = k_BC
               density = densBC
            end if
         end do

      else if (index(mode_, "ocphilic") > 0) then ! Organic Carbon
         do i = 1, size(aeroList)
            if (index(aeroList(i), "CA.oc") > 0) then
               call ESMF_StateGet(state, trim(aeroList(i)), child_state, _RC)
               varNameLen = len_trim(aeroList(i))
               ! the "5" refers to "_AERO", which we want to remove
               ! to get the CA component name (e.g. CA.oc, or CA.oc.data)
               varNameLen = varNameLen - 5
               call MAPL_StateGetPointer(child_state, itemName=aeroList(i)(1:varNameLen)//"philic", farrayPtr=ptr_3d, _RC)
               q = q + ptr_3d
               hygroscopicity = k_OC
               density = densOC
            end if
         end do

      else if (index(mode_, "brcphilic") > 0) then ! Organic Carbon
         do i = 1, size(aeroList)
            if (index(aeroList(i), "CA.br") > 0) then
               call ESMF_StateGet(state, trim(aeroList(i)), child_state, _RC)
               varNameLen = len_trim(aeroList(i))
               ! the "5" refers to "_AERO", which we want to remove
               ! to get the CA component name (e.g. CA.bc, or CA.bc.data)
               varNameLen = varNameLen - 5
               call MAPL_StateGetPointer(child_state, itemName=aeroList(i)(1:varNameLen)//"philic", farrayPtr=ptr_3d, _RC)
               q = q + ptr_3d
               hygroscopicity = k_BR
               density = densBR
            end if
         end do

      end if !(index(mode_, "du00") > 0) then

      ! Obtain aerosol activation properties of this aerosol mode
      call aap_( &
           mode,               &
           q,                  &
           num,                &
           diameter,           &
           sigma,              &
           f_dust,             &
           f_soot,             &
           f_organic,          &
           density,            &
           ptr_3d,             &
           1, i2, 1, j2, km,   &
           _RC)

      deallocate(q, _STAT)

      _RETURN(_SUCCESS)

   contains

      subroutine aap_(mode, q, num, diameter, sigma, f_dust, f_soot, f_organic, dens_, q_, &
           i1, i2, j1, j2, km, rc)

         integer, intent(in) :: i1, i2                             ! dimension bounds
         integer, intent(in) :: j1, j2                             ! ... // ..
         integer, intent(in) :: km                                 ! ... // ..

         character(len=*), intent(in ) :: mode                     ! name of aerosol mode
         real, intent(in), dimension(i1:i2,j1:j2,km) :: q          ! aerosol mass mixing ratio, kg kg-1
         real, intent(in), dimension(i1:i2,j1:j2,km) :: q_         ! auxiliary mass
         real, intent(in), dimension(i1:i2,j1:j2,km) :: dens_      ! density

         real, intent(out), dimension(i1:i2,j1:j2,km) :: num       ! number concentration of aerosol particles
         real, intent(out), dimension(i1:i2,j1:j2,km) :: diameter  ! dry size of aerosol
         real, intent(out), dimension(i1:i2,j1:j2,km) :: sigma     ! width of aerosol mode
         real, intent(out), dimension(i1:i2,j1:j2,km) :: f_dust    ! fraction of dust aerosol
         real, intent(out), dimension(i1:i2,j1:j2,km) :: f_soot    ! fraction of soot aerosol
         real, intent(out), dimension(i1:i2,j1:j2,km) :: f_organic ! fraction of organic aerosol

         integer, intent(out) :: rc

         ! local
         integer :: status
         character(len=:), allocatable :: mode_

         integer, parameter :: UNKNOWN_AEROSOL_MODE = 2015

         integer :: kinx
         real :: fmassaux, fmassclean
         real, dimension(3) :: TPI, DPGI, SIGI
         real, dimension(3) :: TPIclean, DPGIclean, SIGIclean
         real, dimension(i1:i2,j1:j2,km) :: qaux
         !real, parameter    :: max_clean = 5.0e-7  !max mixing ratio before considered polluted

         mode_ = trim(mode)
         mode_ = ESMF_UtilStringLowerCase(mode_, _RC)

         num       = 0.0
         diameter  = 1.0e-9
         sigma     = log(2.0)
         f_dust    = 0.0
         f_soot    = 0.0
         f_organic = 0.0

         qaux=q !this corrects a bug

         if (index(mode_, "ss00") > 0) then
            TPI  (1) = 230e6          ! num fraction (reduced 091015)
            DPGI (1) = 0.02e-6        ! modal diameter (m)
            SIGI (1) = log(1.6)       ! geometric dispersion (sigma_g)
            ! accumulation
            TPI  (2) = 60.0e6         ! total concentration (# m-3)
            DPGI (2) = 0.071e-6       ! modal diameter (m)
            SIGI (2) = log(2.0)       ! geometric dispersion (sigma_g)
            ! coarse
            TPI  (3) = 3.1e6          ! total concentration (# m-3)
            DPGI (3) = 0.62e-6        ! modal diameter (m)
            SIGI (3) = log(2.7)       ! geometric dispersion (sigma_g)

            fmassaux = 0.0
            do kinx = 1, 3
               fmassaux = (TPI(kinx)*densSS*MAPL_PI*exp(4.5*SIGI(kinx)*SIGI(kinx))*DPGI(kinx)*DPGI(kinx)*DPGI(kinx))/6.0 + fmassaux
            end do
         end if

         if (index(mode_, "sulforg0") > 0) then
            TPI  (1) = 1.06e11        ! num fraction
            DPGI (1) = .014e-6        ! modal diameter (m)
            SIGI (1) = log(1.8)       ! geometric dispersion (sigma_g)
            ! accumulation
            TPI  (2) = 3.2e10         ! total concentration (# m-3)
            DPGI (2) = 0.054e-6       ! modal diameter (m)
            SIGI (2) = log(2.16)      ! geometric dispersion (sigma_g)
            !coarse
            TPI  (3) = 5.4e6          ! total concentration (# m-3)
            DPGI (3) = 0.86e-6        ! modal diameter (m)
            SIGI (3) = log(2.21)      ! geometric dispersion (sigma_g)

            fmassaux = 0.0
            do kinx = 1, 3
               ! density is multiplied below since this is a case of a 3-d field
               fmassaux = (TPI(kinx)*MAPL_PI*exp(4.5*SIGI(kinx)*SIGI(kinx))*DPGI(kinx)*DPGI(kinx)*DPGI(kinx))/6.0 + fmassaux
            end do

            ! clean continental polluted plus org
            ! fine
            TPIclean  (1) = 1.0e9      ! total concentration (# m-3)
            DPGIclean (1) = 0.016e-6   ! modal diameter (m)
            SIGIclean (1) = log(1.6)   ! geometric dispersion (sigma_g)
            ! accumulation
            TPIclean  (2) = 8.0e8      ! total concentration (# m-3)
            DPGIclean (2) = 0.067e-6   ! modal diameter (m)
            SIGIclean (2) = log(2.1)   ! geometric dispersion (sigma_g)
            !Coarse
            TPIclean  (3) = 2.0e6      ! total concentration (# m-3)
            DPGIclean (3) = 0.93e-6    ! modal diameter (m)
            SIGIclean (3) = log(2.2)   ! geometric dispersion (sigma_g)

            fmassclean= 0.0
            do kinx = 1, 3
               fmassclean = (TPIclean(kinx)*MAPL_PI*exp(4.5*SIGIclean(kinx)*SIGIclean(kinx))*DPGIclean(kinx)*DPGIclean(kinx)*DPGIclean(kinx))/6.0 + fmassclean  !
            end do
         end if

         select case(mode_)

         case ("du001")
            sigma    = log(1.8)
            f_dust   = 1.0
            diameter = 1.46e-6
            num      = q / ((MAPL_PI/6.0) * densDU * diameter*diameter*diameter * exp(4.5*sigma*sigma))

         case ("du002")
            sigma    = log(1.8)
            f_dust   = 1.0
            diameter = 2.80e-6
            num      = q / ((MAPL_PI/6.0) * densDU * diameter*diameter*diameter * exp(4.5*sigma*sigma))

         case ("du003")
            sigma    = log(1.8)
            f_dust   = 1.0
            diameter = 4.80e-6
            num      = q / ((MAPL_PI/6.0) * densDU * diameter*diameter*diameter * exp(4.5*sigma*sigma))

         case ("du004")
            sigma    = log(1.8)
            f_dust   = 1.0
            diameter = 9.0e-6
            num      = q / ((MAPL_PI/6.0) * densDU * diameter*diameter*diameter * exp(4.5*sigma*sigma))

         case ("du005")
            sigma    = log(1.8)
            f_dust   = 1.0
            diameter = 16.0e-6
            num      = q / ((MAPL_PI/6.0) * densDU * diameter*diameter*diameter * exp(4.5*sigma*sigma))

         case ("ss001")
            sigma    = SIGI(1)
            diameter = DPGI(1)
            num      = TPI(1) * q / fmassaux

         case ("ss002")
            sigma    = SIGI(2)
            diameter = DPGI(2)
            num      = TPI(2) * q / fmassaux

         case ("ss003")
            sigma    = SIGI(3)
            diameter = DPGI(3)
            num      = TPI(3) * q / fmassaux

         case ("sulforg01")  !different distributions for clean and polluted environments
            where (q > max_clean)
               sigma    = SIGI(1)
               diameter = DPGI(1)
               num      = TPI(1) * qaux*ccn_tuning / (dens_*fmassaux)             ! only sulfate  mass
            elsewhere
               sigma    = SIGIclean(1)
               diameter = DPGIclean(1)
               num      = TPIclean(1) * qaux*ccn_tuning / (dens_*fmassclean)      ! only sulfate
            end where

         case ("sulforg02")
            where (q > max_clean)
               sigma    = SIGI(2)
               diameter = DPGI(2)
               num      = TPI(2) * qaux*ccn_tuning / (dens_*fmassaux)            ! only sulfate mass
            elsewhere
               sigma    = SIGIclean(2)
               diameter = DPGIclean(2)
               num      = TPIclean(2) * qaux*ccn_tuning / (dens_*fmassclean)     ! only sulfate
            end where

         case ("sulforg03")
            where (q > max_clean)
               sigma    = SIGI(3)
               diameter = DPGI(3)
               num      = TPI(3) * qaux*ccn_tuning / (dens_*fmassaux)           ! only sulfate mass
            elsewhere
               sigma    = SIGIclean(3)
               diameter = DPGIclean(3)
               num      = TPIclean(3) * qaux*ccn_tuning / (dens_*fmassclean)    ! only sulfate
            end where

         case ("bcphilic")
            sigma    = log(2.0)
            f_soot   = 1.0
            diameter = 0.0118*2e-6
            num = q / ((MAPL_PI/6.0) * densBC * diameter*diameter*diameter * exp(4.5*sigma*sigma))

         case ("ocphilic")
            sigma     = log(2.2)
            f_organic = 1.0
            diameter  = 0.0212*2.0e-6
            num = q / ((MAPL_PI/6.0) * densOrg * diameter*diameter*diameter * exp(4.5*sigma*sigma))

         case ("brcphilic")
            sigma     = log(2.2)
            f_organic = 1.0
            diameter  = 0.0212*2.0e-6
            num = q / ((MAPL_PI/6.0) * densOrg * diameter*diameter*diameter * exp(4.5*sigma*sigma))

         case default
            _FAIL("Unknown aerosol mode used in the GOCART aerosol activation properties method: "//trim(mode))

         end select


         _RETURN(_SUCCESS)

      end subroutine aap_

   end subroutine aerosol_activation_properties

   subroutine get_monochromatic_aop (state, rc)

      implicit none

      !ARGUMENTS:
      type(ESMF_State) :: state
      integer, intent(out) :: rc

      !Local
      real, dimension(:,:,:), pointer :: ple
      real, dimension(:,:,:), pointer :: rh
      real, dimension(:,:), pointer :: var
      character(len=:), allocatable :: fld_name
      real, dimension(:,:), pointer :: tau_      ! (lon:,lat:,lev:)
      real, dimension(:,:), allocatable :: tau   ! (lon:,lat:,lev:)
      integer :: i, n, b, j
      integer :: i1, j1, i2, j2, km
      real :: wavelength
      character(len=ESMF_MAXSTR), allocatable :: itemList(:), aeroList(:)
      type(ESMF_State) :: child_state
      real, pointer, dimension(:,:,:) :: as_ptr_3d
      type(ESMF_StateItem_Flag), allocatable :: itemTypes(:)
      type(ESMF_Info) :: info, child_info
      integer :: status

      ! Description: Used in GAAS gridded component to provide aerosol properties

      call ESMF_InfoGetFromHost(state, info, _RC)

      ! Radiation band
      call ESMF_InfoGet(info, key="wavelength_for_aerosol_optics", value=wavelength, _RC)

      ! Relative humidity
      call ESMF_InfoGet(info, key="relative_humidity_for_aerosol_optics", value=fld_name, _RC)
      call MAPL_StateGetPointer(state, itemName=fld_name, farrayPtr=rh, _RC)

      ! Pressure at layer edges
      call ESMF_InfoGet(info, key="air_pressure_for_aerosol_optics", value=fld_name, _RC)
      call MAPL_StateGetPointer(state, itemName=fld_name, farrayPtr=ple, _RC)

      ! TODO: pchakrab - CAREFUL! ple in MAPL3 is (:, :, 1:km+1), instead of (:, :, km)
      i1 = lbound(ple, 1); i2 = ubound(ple, 1)
      j1 = lbound(ple, 2); j2 = ubound(ple, 2)
      km = ubound(ple, 3)

      allocate(tau(i1:i2,j1:j2), _STAT)
      tau = 0.0

      ! Get list of child states within state and add to aeroList
      call ESMF_StateGet (state, itemCount=n, _RC)
      allocate(itemList(n), _STAT)
      allocate(itemTypes(n), _STAT)
      call ESMF_StateGet(state, itemNameList=itemList, itemTypeList=itemTypes, _RC)

      b=0
      do i = 1, n
         if (itemTypes(i) == ESMF_StateItem_State) then
            b = b + 1
         end if
      end do

      allocate(aeroList(b), _STAT)

      j = 1
      do i = 1, n
         if (itemTypes(i) == ESMF_StateItem_State) then
            aeroList(j) = trim(itemList(i))
            j = j + 1
         end if
      end do

      ! Get aerosol optic properties from children
      do i = 1, size(aeroList)
         call ESMF_StateGet(state, trim(aeroList(i)), child_state, _RC)
         call ESMF_InfoGetFromHost(child_state, child_info, _RC)

         ! set RH in child's aero state
         call ESMF_InfoGet(child_info, key="relative_humidity_for_aerosol_optics", value=fld_name, _RC)

         if (fld_name /= "") then
            call MAPL_StateGetPointer(child_state, itemName=fld_name, farrayPtr=as_ptr_3d, _RC)
            as_ptr_3d = rh
         end if

         ! set PLE in child's aero state
         call ESMF_InfoGet(child_info, key="air_pressure_for_aerosol_optics", value=fld_name, _RC)

         if (fld_name /= "") then
            call MAPL_StateGetPointer(child_state, itemName=fld_name, farrayPtr=as_ptr_3d, _RC)
            as_ptr_3d = ple
         end if

         ! set wavelength in child's aero state
         call ESMF_InfoSet(child_info, key="wavelength_for_aerosol_optics", value=wavelength, _RC)

         ! execute the aerosol optics method
         call ESMF_MethodExecute(child_state, label="monochromatic_aerosol_optics", _RC)

         ! Retrieve extinction from each child
         call ESMF_InfoGet(child_info, key="monochromatic_extinction_in_air_due_to_ambient_aerosol", value=fld_name, _RC)
         if (fld_name /= "") then
            call MAPL_StateGetPointer(child_state, itemName=fld_name, farrayPtr=tau_, _RC)
         end if

         ! Sum aerosol optic properties from each child
         tau = tau + tau_
      end do

      ! Set ext, ssa, asy to equal the sum of ext, ssa, asy from the children. This is what is passed to radiation.
      call ESMF_InfoGet(info, key="monochromatic_extinction_in_air_due_to_ambient_aerosol", value=fld_name, _RC)
      if (fld_name /= "") then
         call MAPL_StateGetPointer(state, itemName=fld_name, farrayPtr=var, _RC)
         var = tau
      end if

      deallocate(tau, _STAT)

      _RETURN(_SUCCESS)

   end subroutine get_monochromatic_aop

   subroutine get_mixRatioSum (state, rc)

      implicit none

      !ARGUMENTS:
      type(ESMF_State) :: state
      integer, intent(out) :: rc

      !Local
      character(len=ESMF_MAXSTR), allocatable :: itemList(:), aeroList(:)
      character(len=:), allocatable :: aeroName, fld_name

      real, pointer, dimension(:,:,:) :: var
      real, dimension(:,:,:), allocatable :: aeroOut
      type(ESMF_StateItem_Flag), allocatable :: itemTypes(:)
      type(ESMF_Info) :: info

      integer :: b, i, n, j, im, jm, km, status

      ! Description: Used in GAAS gridded component to provide sum of aerosol mixing ratio

      call ESMF_InfoGetFromHost(state, info, _RC)

      call ESMF_InfoGet(info, key="aerosolName", value=aeroName, _RC)
      call ESMF_InfoGet(info, key="im", value=im, _RC)
      call ESMF_InfoGet(info, key="jm", value=jm, _RC)
      call ESMF_InfoGet(info, key="km", value=km, _RC)

      allocate(aeroOut(im,jm,km), _STAT)
      aeroOut = 0.0

      ! Get list of child states within state and add to aeroList
      call ESMF_StateGet(state, itemCount=n, _RC)
      allocate(itemList(n), _STAT)
      allocate(itemTypes(n), _STAT)
      call ESMF_StateGet(state, itemNameList=itemList, itemTypeList=itemTypes, _RC)

      b=0
      do i = 1, n
         if (itemTypes(i) == ESMF_StateItem_State) then
            b = b + 1
         end if
      end do

      allocate(aeroList(b), _STAT)

      j = 1
      do i = 1, n
         if (itemTypes(i) == ESMF_StateItem_State) then
            aeroList(j) = trim(itemList(i))
            j = j + 1
         end if
      end do

      ! Retrieve summed aerosol mixing ratios from active instances
      select case (trim(aeroName))
      case ("dust")
         call getAerosolSum("DU", state, aeroList, aeroOut, _RC)

         call ESMF_InfoGet(info, key="sum_of_internalState_aerosol_DU", value=fld_name, _RC)
         if (fld_name /= "") then
            call MAPL_StateGetPointer(state, itemName=fld_name, farrayPtr=var, _RC)
            var = aeroOut
         end if

      case ("seasalt")
         call getAerosolSum("SS", state, aeroList, aeroOut, _RC)

         call ESMF_InfoGet(info, key="sum_of_internalState_aerosol_SS", value=fld_name, _RC)
         if (fld_name /= "") then
            call MAPL_StateGetPointer(state, itemName=fld_name, farrayPtr=var, _RC)
            var = aeroOut
         end if

      case ("organicCarbon")
         call getAerosolSum("CA.oc", state, aeroList, aeroOut, _RC)

         call ESMF_InfoGet(info, key="sum_of_internalState_aerosol_CA.oc", value=fld_name, _RC)
         if (fld_name /= "") then
            call MAPL_StateGetPointer(state, itemName=fld_name, farrayPtr=var, _RC)
            var = aeroOut
         end if

      case ("blackCarbon")
         call getAerosolSum("CA.bc", state, aeroList, aeroOut, _RC)

         call ESMF_InfoGet(info, key="sum_of_internalState_aerosol_CA.bc", value=fld_name, _RC)
         if (fld_name /= "") then
            call MAPL_StateGetPointer(state, itemName=fld_name, farrayPtr=var, _RC)
            var = aeroOut
         end if

      case ("brownCarbon")
         call getAerosolSum ("CA.br", state, aeroList, aeroOut, _RC)

         call ESMF_InfoGet(info, key="sum_of_internalState_aerosol_CA.br", value=fld_name, _RC)
         if (fld_name /= "") then
            call MAPL_StateGetPointer(state, itemName=fld_name, farrayPtr=var, _RC)
            var = aeroOut
         end if

      case ("sulfate")
         call getAerosolSum("SU", state, aeroList, aeroOut, _RC)

         call ESMF_InfoGet(info, key="sum_of_internalState_aerosol_SU", value=fld_name, _RC)
         if (fld_name /= "") then
            call MAPL_StateGetPointer(state, itemName=fld_name, farrayPtr=var, _RC)
            var = aeroOut
         end if

      case ("nitrate")
         call getAerosolSum("NI", state, aeroList, aeroOut, _RC)

         call ESMF_InfoGet(info, key="sum_of_internalState_aerosol_NI", value=fld_name, _RC)
         if (fld_name /= "") then
            call MAPL_StateGetPointer(state, itemName=fld_name, farrayPtr=var, _RC)
            var = aeroOut
         end if

      case default
         !$omp critical (G2G_2)
         print *,"Invalid aerosolName of "",trim(aeroName), "" in GOCART2G::get_mixRatioSum"
         !$omp end critical (G2G_2)
      end select

      _RETURN(_SUCCESS)

   contains

      subroutine getAerosolSum(aeroToken, state, aeroList, aeroOut, rc)
         !ARGUMENTS:
         character(len=*), intent(in) :: aeroToken
         type(ESMF_State), intent(in) :: state
         character(len=ESMF_MAXSTR), intent(in) :: aeroList(:)
         real, dimension(:,:,:), intent(out) :: aeroOut
         integer, optional, intent(out) :: rc

         !LOCALS:
         integer :: i, endInd
         character(len=ESMF_MAXSTR) :: fld_name
         type(ESMF_State) :: child_state
         type(ESMF_Info) :: child_info
         real, pointer, dimension(:,:,:) :: ptr3d

         endInd = len_trim(aeroToken)

         aeroOut = 0.0
         do i = 1, size(aeroList)
            if (trim(aeroList(i)(1:endInd)) == trim(aeroToken)) then
               call ESMF_StateGet(state, trim(aeroList(i)), child_state, _RC)
               call ESMF_MethodExecute(child_state, label="get_mixR", _RC)
               call ESMF_InfoGetFromHost(child_state, child_info, _RC)
               call ESMF_InfoGet(child_info, key="sum_of_internalState_aerosol", value=fld_name, _RC)
               if (fld_name /= "") then
                  call MAPL_StateGetPointer(child_state, itemName=fld_name, farrayPtr=ptr3d, _RC)
                  aeroOut = aeroOut + ptr3d
               end if
            end if
         end do

         _RETURN(_SUCCESS)
      end subroutine getAerosolSum

   end subroutine get_mixRatioSum

end module GOCART2G_GridCompMod

subroutine SetServices(gc, rc)
   use ESMF
   use GOCART2G_GridCompMod, only : mySetservices=>SetServices
   type(ESMF_GridComp) :: gc
   integer, intent(out) :: rc
   call mySetServices(gc, rc=rc)
end subroutine SetServices
