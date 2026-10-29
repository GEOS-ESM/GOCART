#include "MAPL.h"

!BOP
! !MODULE: DuSsStub_GridCompMod
!
! !DESCRIPTION: A minimal, dedicated stand-in gridded component used only by
!   the GOCART2G regression/nitrate-sa standalone harness. NI2G_GridCompMod's
!   Initialize unconditionally reads custom ESMF Info attributes ('radius',
!   'fnum') off its 'DU'/'SS' import fields. Those attributes are only ever
!   attached by the real DU2G/SS2G components' own Initialize routines
!   (ESMF_InfoSet onto their exported DU/SS fields) -- the generic
!   configurable_gridcomp stand-in used elsewhere in this harness has no
!   mechanism to do that. This tiny component imports 'DU'/'SS' (connected to
!   the DataForNI internal state, i.e. the same import netcdf4 file the
!   harness's other import fields are read from) and re-exports them,
!   attaching 'radius'/'fnum' Info exactly like DU2G/SS2G do, so that NI2G's
!   Initialize can complete when run outside the full GOCART2G parent.
!
!EOP

module mapl_DuSsStubGridComp_mod

   use MAPL
   use ESMF

   implicit none
   private

   public :: SetServices

   integer, parameter :: NBINS = 5

   ! Taken from DU2G_instance_DU.yaml / SS2G_instance_SS.yaml
   ! ('particle_radius_microns' / 'fnum'), which is exactly what the real
   ! DU2G/SS2G components attach to their exported DU/SS fields.
   real, parameter :: DU_RADIUS(NBINS) = [0.73, 1.4, 2.4, 4.5, 8.0]
   real, parameter :: DU_FNUM(NBINS) = [2.45e14, 3.28e13, 6.52e12, 9.89e11, 1.76e11]
   real, parameter :: SS_RADIUS(NBINS) = [0.079, 0.316, 1.119, 2.818, 7.772]
   real, parameter :: SS_FNUM(NBINS) = [3.017e17, 1.085e16, 1.207e14, 9.391e12, 2.922e11]

contains

   subroutine SetServices(gc, rc)
      type(ESMF_GridComp) :: gc
      integer, intent(out) :: rc

      integer :: status
      type(MAPL_UngriddedDim) :: ungrd_bins

      ungrd_bins = MAPL_UngriddedDim(NBINS, name="nbins", units="1")

      ! Imports: connected (via root.yaml) to DataForNI's internal DU/SS,
      ! i.e. the same import netcdf4 file that DataForNI's other fields
      ! (DELP, EMI_NH3_*, ...) are read from.
      call MAPL_GridCompAddSpec(gc, &
           state_intent=ESMF_STATEINTENT_IMPORT, &
           short_name="DU", &
           standard_name="Dust Mixing Ratio all bins", &
           units="kg kg-1", &
           dims="xyz", &
           vertical_stagger=MAPL_VERTICAL_STAGGER_CENTER, &
           ungridded_dim_array=[ungrd_bins], &
           _RC)

      call MAPL_GridCompAddSpec(gc, &
           state_intent=ESMF_STATEINTENT_IMPORT, &
           short_name="SS", &
           standard_name="Sea Salt Mixing Ratio all bins", &
           units="kg kg-1", &
           dims="xyz", &
           vertical_stagger=MAPL_VERTICAL_STAGGER_CENTER, &
           ungridded_dim_array=[ungrd_bins], &
           _RC)

      ! Exports: values copied from the imports above (see Initialize),
      ! with 'radius'/'fnum' Info attached, standing in for DU2G/SS2G.
      call MAPL_GridCompAddSpec(gc, &
           state_intent=ESMF_STATEINTENT_EXPORT, &
           short_name="DU", &
           standard_name="Dust Mixing Ratio all bins", &
           units="kg kg-1", &
           dims="xyz", &
           vertical_stagger=MAPL_VERTICAL_STAGGER_CENTER, &
           ungridded_dim_array=[ungrd_bins], &
           fill_value=0.0, _RC)

      call MAPL_GridCompAddSpec(gc, &
           state_intent=ESMF_STATEINTENT_EXPORT, &
           short_name="SS", &
           standard_name="Sea Salt Mixing Ratio all bins", &
           units="kg kg-1", &
           dims="xyz", &
           vertical_stagger=MAPL_VERTICAL_STAGGER_CENTER, &
           ungridded_dim_array=[ungrd_bins], &
           fill_value=0.0, _RC)

      call MAPL_GridCompSetEntryPoint(gc, ESMF_METHOD_INITIALIZE, Initialize, _RC)
      call MAPL_GridCompSetEntryPoint(gc, ESMF_METHOD_RUN, Run, phase_name="run", _RC)

      _RETURN(_SUCCESS)
   end subroutine SetServices

   subroutine Initialize(gc, import, export, clock, rc)
      type(ESMF_GridComp) :: gc
      type(ESMF_State) :: import
      type(ESMF_State) :: export
      type(ESMF_Clock) :: clock
      integer, intent(out) :: rc

      integer :: status
      type(ESMF_Field) :: field
      type(ESMF_Info) :: field_info
      real, pointer :: du_import(:, :, :, :), du_export(:, :, :, :)
      real, pointer :: ss_import(:, :, :, :), ss_export(:, :, :, :)

      ! Copy the real DU/SS values in from the connected imports (DataForNI's
      ! internal state, read from the same import netcdf4 file as the
      ! harness's other input fields) onto the exports.
      call MAPL_StateGetPointer(import, du_import, "DU", _RC)
      call MAPL_StateGetPointer(export, du_export, "DU", _RC)
      du_export = du_import

      call MAPL_StateGetPointer(import, ss_import, "SS", _RC)
      call MAPL_StateGetPointer(export, ss_export, "SS", _RC)
      ss_export = ss_import

      call ESMF_StateGet(export, itemName="DU", field=field, _RC)
      call ESMF_InfoGetFromHost(field, field_info, _RC)
      call ESMF_InfoSet(field_info, key="radius", values=DU_RADIUS, _RC)
      call ESMF_InfoSet(field_info, key="fnum", values=DU_FNUM, _RC)

      call ESMF_StateGet(export, itemName="SS", field=field, _RC)
      call ESMF_InfoGetFromHost(field, field_info, _RC)
      call ESMF_InfoSet(field_info, key="radius", values=SS_RADIUS, _RC)
      call ESMF_InfoSet(field_info, key="fnum", values=SS_FNUM, _RC)

      _RETURN(_SUCCESS)
      _UNUSED_DUMMY(gc)
      _UNUSED_DUMMY(clock)
   end subroutine Initialize

   recursive subroutine Run(gc, import, export, clock, rc)
      type(ESMF_GridComp) :: gc
      type(ESMF_State) :: import
      type(ESMF_State) :: export
      type(ESMF_Clock) :: clock
      integer, intent(out) :: rc

      ! Nothing to do: DU/SS are copied from the imports once, at Initialize.
      _RETURN(_SUCCESS)
      _UNUSED_DUMMY(gc)
      _UNUSED_DUMMY(import)
      _UNUSED_DUMMY(export)
      _UNUSED_DUMMY(clock)
   end subroutine Run

end module mapl_DuSsStubGridComp_mod

subroutine setServices(gc, rc)
   use ESMF
   use MAPL
   use mapl_DuSsStubGridComp_mod, only: DuSsStub_SetServices => SetServices
   type(ESMF_GridComp) :: gc
   integer, intent(out) :: rc

   integer :: status

   call DuSsStub_SetServices(gc, _RC)

   _RETURN(_SUCCESS)
end subroutine setServices
