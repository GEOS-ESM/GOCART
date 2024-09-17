#include "MAPL_Generic.h"

!=============================================================================
!BOP

! !MODULE: GOCART2G_GridCompMod - The GOCART 2nd Generation Aerosol Grid Component

! !INTERFACE:

module GOCART2G_GridCompMod

! !USES:

   use ESMF
   use MAPL
   use Chem_AeroGeneric

! !Establish the Childen's SetServices
 !-----------------------------------
   use DU2G_GridCompMod,    only   : DU2G_setServices  => SetServices
   use SS2G_GridCompMod,    only   : SS2G_setServices  => SetServices
   use SU2G_GridCompMod,    only   : SU2G_setServices  => SetServices
   use CA2G_GridCompMod,    only   : CA2G_setServices  => SetServices
   use NI2G_GridCompMod,    only   : NI2G_setServices  => SetServices

   implicit none
   private

! !PUBLIC MEMBER FUNCTIONS:
   public  SetServices

  ! Private State
  type :: Instance
     integer :: id = -1
     logical :: is_active
     character(:), allocatable :: name
  end type Instance

  type Constituent
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

  type wrap_
     type (GOCART_State), pointer     :: PTR => null()
  end type wrap_

! !DESCRIPTION:
!
!   {\tt GOCART} is a gridded component from the GOCART model and includes
!  dust, sea salt, sulfates, nitrate, organic and black carbon.

!
!
! !REVISION HISTORY:
!
!  25feb2005  da Silva   First crack.
!  19jul2006  da Silva   First separate GOCART component.
!  14Oct2019  E.Sherman, A.Darmenov, A. da Silva, T. Clune  First attempt at refactoring.
!
!EOP
!============================================================================

contains

!BOP

! !IROUTINE: SetServices -- Sets ESMF services for this component

! !INTERFACE:

  subroutine SetServices (GC, RC)

! !ARGUMENTS:

    type (ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, optional                   :: RC  ! return code

! !DESCRIPTION: This version uses MAPL_GenericSetServices, which sets
!   the Initialize and Finalize services to generic versions. It also
!   allocates our instance of a generic state and puts it in the
!   gridded component (GC). Here we only set the two-stage run method and
!   declare the data services.

! !REVISION HISTORY:

!  14oct2019  Sherman, da Silva, Darmenov, Clune - First attempt at refactoring for ESMF compatibility


!EOP
!============================================================================
!
!   Locals
    character (len=ESMF_MAXSTR)                   :: COMP_NAME
    type (ESMF_Config)                            :: myCF      ! GOCART2G_GridComp.rc
    type (ESMF_Config)                            :: cf        ! universal config
    type (GOCART_State), pointer                  :: self
    type (wrap_)                                  :: wrap

    integer :: n_wavelengths_profile, n_wavelengths_vertint, n_wavelengths_diagmie
    integer, allocatable, dimension(:) :: wavelengths_diagmie
    type (MAPL_MetaComp),       pointer    :: MAPL
    logical :: use_threads

    __Iam__('SetServices')

!****************************************************************************
! Begin...

!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet (GC, name=comp_name, config=cf, __RC__)
    Iam = trim(comp_name)//'::'//'SetServices'

!   Wrap internal state for storing in GC
!   -------------------------------------
    allocate (self, __STAT__)
    wrap%ptr => self

!   Set the Initialize, Run, Finalize entry points
!   ------------------------------------------------
    call MAPL_GridCompSetEntryPoint (GC, ESMF_Method_Initialize,  Initialize,  __RC__)
    call MAPL_GridCompSetEntryPoint (GC, ESMF_Method_Run,  Run1, __RC__)
    call MAPL_GridCompSetEntryPoint (GC, ESMF_Method_Run,  Run2, __RC__)

!   Store internal state in GC
!   --------------------------
    call ESMF_UserCompSetInternalState (GC, 'GOCART_State', wrap, STATUS)
    VERIFY_(STATUS)

    myCF = ESMF_ConfigCreate (__RC__)
    call ESMF_ConfigLoadFile (myCF, 'GOCART2G_GridComp.rc', __RC__)

!   Retrieve wavelengths from GOCART2G_GridComp.rc
    n_wavelengths_profile = ESMF_ConfigGetLen (myCF, label='wavelengths_for_profile_aop_in_nm:', __RC__)
    n_wavelengths_vertint = ESMF_ConfigGetLen (myCF, label='wavelengths_for_vertically_integrated_aop_in_nm:', __RC__)
    n_wavelengths_diagmie = ESMF_ConfigGetLen (myCF, label='aerosol_monochromatic_optics_wavelength_in_nm_from_LUT:', __RC__)

    allocate(self%wavelengths_profile(n_wavelengths_profile), self%wavelengths_vertint(n_wavelengths_vertint), &
             wavelengths_diagmie(n_wavelengths_diagmie), __STAT__)

    call ESMF_ConfigGetAttribute (myCF, self%wavelengths_profile, label='wavelengths_for_profile_aop_in_nm:', __RC__)
    call ESMF_ConfigGetAttribute (myCF, self%wavelengths_vertint, label='wavelengths_for_vertically_integrated_aop_in_nm:', __RC__)
    call ESMF_ConfigGetAttribute (myCF, wavelengths_diagmie, label='aerosol_monochromatic_optics_wavelength_in_nm_from_LUT:', __RC__)

!   Set wavelengths in universal config

    call MAPL_ConfigSetAttribute (cf, self%wavelengths_profile, label='wavelengths_for_profile_aop_in_nm:', __RC__)
    call MAPL_ConfigSetAttribute (cf, self%wavelengths_vertint, label='wavelengths_for_vertically_integrated_aop_in_nm:', __RC__)
    call MAPL_ConfigSetAttribute (cf, wavelengths_diagmie, label='aerosol_monochromatic_optics_wavelength_in_nm_from_LUT:', __RC__)
    call ESMF_ConfigGetAttribute (myCF, use_threads, label='use_threads:', default=.FALSE., __RC__)

!   Get my internal MAPL_Generic state
!   -----------------------------------
    call MAPL_GetObjectFromGC (GC, MAPL, __RC__)
!   set use_threads
    call MAPL%set_use_threads(use_threads)

!   Get instances to determine what children will be born
!   -----------------------------------------------------
    call getInstances_('DU', myCF, species=self%DU, __RC__)
    call getInstances_('SS', myCF, species=self%SS, __RC__)
    call getInstances_('SU', myCF, species=self%SU, __RC__)
    call getInstances_('CA', myCF, species=self%CA, __RC__)
    call getInstances_('NI', myCF, species=self%NI, __RC__)

!   Nitrate currently only supports one active instance
    if (self%NI%n_active > 1) then
       if(mapl_am_i_root()) print*,'WARNING: GOCART can only support one active nitrate instance. Check the RC/GOCART2G_GridComp.rc'
    end if

    call ESMF_ConfigDestroy(myCF, __RC__)

!   Create children's gridded components and invoke their SetServices
!   Active instances are created first
!   -----------------------------------------------------------------
    call createInstances_(self, GC, __RC__)

!   Define EXPORT states

!   This state is needed by radiation and moist. It contains
!   aerosols and callback methods
!   --------------------------------------------------------
    call MAPL_AddExportSpec(GC,                       &
       short_name = 'AERO',                           &
       long_name  = 'aerosol_mass_mixing_ratios_ng',  &
       units      = 'kg kg-1',                        &
       dims       = MAPL_DimsHorzVert,                &
       vlocation  = MAPL_VLocationCenter,             &
       datatype   = MAPL_StateItem, __RC__)

!   This bundle is needed by surface for snow albedo modification
!   by aerosol settling and deposition
!   --------------------------------------------------------
    call MAPL_AddExportSpec(GC,                       &
       short_name = 'AERO_DP',                      &
       long_name  = 'aerosol_deposition_ng',          &
       units      = 'kg m-2 s-1',                     &
       dims       = MAPL_DimsHorzOnly,                &
       datatype   = MAPL_BundleItem, __RC__)


#include "GOCART2G_Export___.h"
#include "GOCART2G_Import___.h"

!   Add connectivities for Nitrate component
!   Nitrate currently only supports one Nitrate component. Nitrate only
!   uses the first active dust and sea salt instance.
    if (size(self%NI%instances) > 0) then
       if ((self%DU%instances(1)%is_active)) then
          call MAPL_AddConnectivity (GC, SHORT_NAME = ["DU"], &
                                     DST_ID=self%NI%instances(1)%id, &
                                     SRC_ID=self%DU%instances(1)%id, __RC__)
       end if

       if ((self%SS%instances(1)%is_active)) then
          call MAPL_AddConnectivity (GC, SHORT_NAME = ["SS"] , &
                                     DST_ID=self%NI%instances(1)%id, &
                                     SRC_ID=self%SS%instances(1)%id, __RC__)
       end if

       if ((self%SU%instances(1)%is_active)) then
          call MAPL_AddConnectivity (GC, SHORT_NAME = ["SO4"] , &
                                     DST_ID=self%NI%instances(1)%id, &
                                     SRC_ID=self%SU%instances(1)%id, __RC__)
       end if
    end if

!   Set generic services
!   ----------------------------------
    call MAPL_GenericSetServices (GC, __RC__)

    RETURN_(ESMF_SUCCESS)

  end subroutine SetServices


!============================================================================
!BOP

! !IROUTINE: Initialize -- Initialize method for the composite Gridded Component

! !INTERFACE:
  subroutine Initialize (GC, import, export, clock, RC)

! !ARGUMENTS:
    type (ESMF_GridComp), intent(inout) :: GC     ! Gridded component
    type (ESMF_State),    intent(inout) :: import ! Import state
    type (ESMF_State),    intent(inout) :: export ! Export state
    type (ESMF_Clock),    intent(inout) :: clock  ! The clock
    integer, optional,    intent(  out) :: RC     ! Error code

! !DESCRIPTION:  This initializes the GOCART Grid Component. It primarily creates
!                its exports and births its children.

! !REVISION HISTORY:
! 14oct2019   E.Sherman  First attempt at refactoring

!EOP
!============================================================================

!   Locals
    character (len=ESMF_MAXSTR)            :: COMP_NAME

    type (MAPL_MetaComp),       pointer    :: MAPL
    type (ESMF_GridComp),       pointer    :: gcs(:)
    type (ESMF_State),          pointer    :: gex(:)
    type (ESMF_Grid)                       :: grid
    type (ESMF_Config)                     :: CF

    type (ESMF_State)                      :: aero
    type (ESMF_FieldBundle)                :: aero_dp

    type (GOCART_State),      pointer      :: self
    type (wrap_)                           :: wrap

    integer                                :: n_modes
    integer, parameter                     :: n_gocart_modes = 14
    integer                                :: dims(3)

    character(len=ESMF_MAXSTR)             :: aero_aci_modes(n_gocart_modes)
    real                                   :: f_aci_seasalt, maxclean, ccntuning
    character(LEN=ESMF_MAXSTR)             :: CLDMICRO

    __Iam__('Initialize')

!****************************************************************************
! Begin...

!   Get the target components name and set-up traceback handle.
!   -----------------------------------------------------------
    call ESMF_GridCompGet (GC, grid=grid, name=COMP_NAME, __RC__)
    Iam = trim(COMP_NAME)//'::'//'Initialize'

    if (mapl_am_i_root()) then
       print *, TRIM(Iam)//': Starting...'
       print *,' '
    end if

!   Get my internal MAPL_Generic state
!   -----------------------------------
    call MAPL_GetObjectFromGC (GC, MAPL, __RC__)

    call MAPL_GridGet ( grid, localCellCountPerDim=dims, __RC__ )

!   Call Generic Initialize
!   ----------------------------------------
    call MAPL_GenericInitialize (GC, import, export, clock, __RC__)

!   Get my internal state
!   ---------------------
    call ESMF_UserCompGetInternalState (GC, 'GOCART_State', wrap, STATUS)
    VERIFY_(STATUS)
    self => wrap%ptr

    CF = ESMF_ConfigCreate (__RC__)
    call ESMF_ConfigLoadFile (CF, 'AGCM.rc', __RC__) ! should the rc file be changed?

!   Get children and their export states from my generic state
!   -----------------------------------------------------------
    call MAPL_Get (MAPL, gcs=gcs, gex=gex, __RC__ )


!   Fill AERO_RAD, AERO_ACI, and AERO_DP with the children's states
!   ---------------------------------------------------------------
    call ESMF_StateGet (export, 'AERO', aero, __RC__)
    call ESMF_StateGet (export, 'AERO_DP', aero_dp, __RC__)


!   Add children's AERO states to GOCART2G's AERO states
!   Only active instances are passed to radiation
!   ------------------------------------------------------
    call add_aero_states_(self%DU%instances(:))
    call add_aero_states_(self%SS%instances(:))
    call add_aero_states_(self%SU%instances(:))
    call add_aero_states_(self%CA%instances(:))
    call add_aero_states_(self%NI%instances(:))

!   Begin AERO_RAD
!   --------------
!   Add variables to AERO_RAD state. Used in aerosol optics calculations
    call add_aero (aero, label='air_pressure_for_aerosol_optics', label2='PLE', &
                   grid=grid, typekind=MAPL_R4, __RC__)
    call add_aero (aero, label='relative_humidity_for_aerosol_optics', label2='RH', &
                   grid=grid, typekind=MAPL_R4, __RC__)
    call add_aero (aero, label='extinction_in_air_due_to_ambient_aerosol', label2='EXT', &
                   grid=grid, typekind=MAPL_R4, __RC__)
    call add_aero (aero, label='single_scattering_albedo_of_ambient_aerosol', label2='SSA', &
                   grid=grid, typekind=MAPL_R4, __RC__)
    call add_aero (aero, label='asymmetry_parameter_of_ambient_aerosol', label2='ASY', &
                   grid=grid, typekind=MAPL_R4, __RC__)
    call add_aero (aero, label='monochromatic_extinction_in_air_due_to_ambient_aerosol', &
                   label2='monochromatic_EXT', grid=grid, typekind=MAPL_R4, __RC__)

!   Used in get_mixRatioSum
    call add_aero (aero, label='sum_of_internalState_aerosol_DU', label2='aerosolSumDU', &
                   grid=grid, typekind=MAPL_R4, __RC__)
    call add_aero (aero, label='sum_of_internalState_aerosol_SS', label2='aerosolSumSS', &
                   grid=grid, typekind=MAPL_R4, __RC__)
    call add_aero (aero, label='sum_of_internalState_aerosol_NI', label2='aerosolSumNI', &
                   grid=grid, typekind=MAPL_R4, __RC__)
    call add_aero (aero, label='sum_of_internalState_aerosol_CA.oc', label2='aerosolSumCA.oc', &
                   grid=grid, typekind=MAPL_R4, __RC__)
    call add_aero (aero, label='sum_of_internalState_aerosol_CA.bc', label2='aerosolSumCA.bc', &
                   grid=grid, typekind=MAPL_R4, __RC__)
    call add_aero (aero, label='sum_of_internalState_aerosol_CA.br', label2='aerosolSumCA.br', &
                   grid=grid, typekind=MAPL_R4, __RC__)
    call add_aero (aero, label='sum_of_internalState_aerosol_SU', label2='aerosolSumSU', &
                   grid=grid, typekind=MAPL_R4, __RC__)

    call ESMF_AttributeSet(aero, name='band_for_aerosol_optics', value=0, __RC__)
    call ESMF_AttributeSet(aero, name='wavelength_for_aerosol_optics', value=0., __RC__)
    call ESMF_AttributeSet(aero, name='aerosolName', value='', __RC__)
    call ESMF_AttributeSet(aero, name='im', value=dims(1), __RC__)
    call ESMF_AttributeSet(aero, name='jm', value=dims(2), __RC__)
    call ESMF_AttributeSet(aero, name='km', value=dims(3), __RC__)

!   Attach method to return sum of aerosols. Used in GAAS.
    call ESMF_MethodAdd (aero, label='get_mixRatioSum', userRoutine=get_mixRatioSum, __RC__)

!   Attach method to create a Bundle of aerosol fields. Used in GAAS.
    call ESMF_MethodAdd (aero, label='serialize_bundle', userRoutine=serialize_bundle, __RC__)

!   Attach the monochromatic aerosol optics method. Used in GAAS.
    call ESMF_MethodAdd (aero, label='get_monochromatic_aop', &
                         userRoutine=get_monochromatic_aop, __RC__)

!   Attach the aerosol optics method. Used in Radiation.
    call ESMF_MethodAdd (aero, label='run_aerosol_optics', userRoutine=run_aerosol_optics, __RC__)

    ! This attribute indicates if the aerosol optics method is implemented or not.
    ! Radiation will not call the aerosol optics method unless this attribute is
    ! explicitly set to true.
    call ESMF_AttributeSet(aero, name='implements_aerosol_optics_method', value=.true., __RC__)

!   Begin adding necessary aerosol cloud interaction information
!   ------------------------------------------------------------
    aero_aci_modes =  (/'du001    ', 'du002    ', 'du003    ', &
                        'du004    ', 'du005    ',              &
                        'ss001    ', 'ss002    ', 'ss003    ', &
                        'sulforg01', 'sulforg02', 'sulforg03', &
                        'bcphilic ', 'ocphilic ', 'brcphilic'/)

    n_modes = size(aero_aci_modes)

    call ESMF_AttributeSet(aero, name='number_of_aerosol_modes', value=n_modes, __RC__)
    call ESMF_AttributeSet(aero, name='aerosol_modes', itemcount=n_modes, valuelist=aero_aci_modes, __RC__)

    ! max mixing ratio before switching to "polluted" size distributions
    call ESMF_ConfigGetAttribute(CF, maxclean, default=1.0e-9, label='MAXCLEAN:', __RC__)
    call ESMF_AttributeSet(aero, name='max_q_clean', value=maxclean, __RC__)

    call ESMF_ConfigGetAttribute(CF, CCNtuning, default=1.8, label='CCNTUNING:', __RC__)
    call ESMF_AttributeSet(aero, name='ccn_tuning', value=CCNtuning, __RC__)

    call ESMF_ConfigGetAttribute( CF, CLDMICRO, Label='CLDMICR_OPTION:',  default="BACM_1M", RC=STATUS)
    call ESMF_AttributeSet(aero, name='cldmicro', value=CLDMICRO, __RC__)

!   Add variables to AERO state
    call add_aero (aero, label='air_temperature', label2='T', grid=grid, typekind=MAPL_R4, __RC__)
    call add_aero (aero, label='fraction_of_land_type', label2='FRLAND', grid=grid, typekind=MAPL_R4, __RC__)
    call add_aero (aero, label='width_of_aerosol_mode', label2='SIGMA', grid=grid, typekind=MAPL_R4, __RC__)
    call add_aero (aero, label='aerosol_number_concentration', label2='NUM', grid=grid, typekind=MAPL_R4, __RC__)
    call add_aero (aero, label='aerosol_dry_size', label2='DGN', grid=grid, typekind=MAPL_R4, __RC__)
    call add_aero (aero, label='aerosol_density', label2='density', grid=grid, typekind=MAPL_R4, __RC__)
    call add_aero (aero, label='aerosol_hygroscopicity', label2='KAPPA', grid=grid, typekind=MAPL_R4, __RC__)
    call add_aero (aero, label='fraction_of_dust_aerosol', label2='FDUST', grid=grid, typekind=MAPL_R4, __RC__)
    call add_aero (aero, label='fraction_of_soot_aerosol', label2='FSOOT', grid=grid, typekind=MAPL_R4, __RC__)
    call add_aero (aero, label='fraction_of_organic_aerosol', label2='FORGANIC', grid=grid, typekind=MAPL_R4, __RC__)

!   Attach the aerosol optics method
    call ESMF_MethodAdd(aero, label='aerosol_activation_properties', userRoutine=aerosol_activation_properties, __RC__)

    RETURN_(ESMF_SUCCESS)

  contains

     subroutine add_aero_states_(instances)
        type(Instance), intent(in) :: instances(:)
        type (ESMF_State)       :: child_state
        type (ESMF_FieldBundle) :: child_bundle
        type (ESMF_Field), allocatable :: fieldList(:)

        integer :: i
        integer :: id
        integer :: fieldCount
        __Iam__('Initialize::add_aero_states_')

        do i = 1, size(instances)
           if (.not. instances(i)%is_active) cycle
           id = instances(i)%id

           call ESMF_GridCompGet (gcs(id), __RC__ )
           call ESMF_StateGet (gex(id), trim(instances(i)%name)//'_AERO', child_state, __RC__)
           call ESMF_StateAdd (aero, [child_state], __RC__)

           if (instances(i)%name(1:2) /= 'NI') then
              call ESMF_StateGet (gex(id), trim(instances(i)%name)//'_AERO_DP', child_bundle, __RC__)
              call ESMF_FieldBundleGet (child_bundle, fieldCount=fieldCount, __RC__)
              allocate (fieldList(fieldCount), __STAT__)
              call ESMF_FieldBundleGet (child_bundle, fieldList=fieldList, __RC__)
              call ESMF_FieldBundleAdd (aero_dp, fieldList, multiflag=.true., __RC__)
              deallocate(fieldList, __STAT__)
           end if
        end do
        RETURN_(ESMF_SUCCESS)
     end subroutine add_aero_states_

 end subroutine Initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!BOP
! !IROUTINE: RUN -- Run method for GOCART2G


! !INTERFACE:

  subroutine Run1 (GC, import, export, clock, RC)

! !ARGUMENTS:
    type (ESMF_GridComp), intent(inout) :: GC     ! Gridded component
    type (ESMF_State),    intent(inout) :: import ! Import state
    type (ESMF_State),    intent(inout) :: export ! Export state
    type (ESMF_Clock),    intent(inout) :: clock  ! The clock
    integer, optional,    intent(  out) :: RC     ! Error code:

! !DESCRIPTION: Run method

!EOP
!============================================================================

!   Locals
    character(len=ESMF_MAXSTR)          :: COMP_NAME
    type (MAPL_MetaComp),      pointer  :: MAPL
    type (ESMF_GridComp),      pointer  :: gcs(:)
    type (ESMF_State),         pointer  :: gim(:)
    type (ESMF_State),         pointer  :: gex(:)
    type (ESMF_State)                   :: internal

    integer                             :: i

    __Iam__('Run1')

!****************************************************************************
! Begin...


!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, __RC__ )
    Iam = trim(COMP_NAME)//'::'//Iam

!   Get my internal MAPL_Generic state
!   -----------------------------------
    call MAPL_GetObjectFromGC ( GC, MAPL, __RC__ )

!   Get parameters from generic state.
!   -----------------------------------
    call MAPL_Get ( MAPL, gcs=gcs, gim=gim, gex=gex, INTERNAL_ESMF_STATE=internal, __RC__ )

!   Run the children
!   -----------------
    do i = 1, size(gcs)
      call ESMF_GridCompRun (gcs(i), importState=gim(i), exportState=gex(i), phase=1, clock=clock, __RC__)
    end do


    RETURN_(ESMF_SUCCESS)

  end subroutine Run1

!============================================================================
!BOP
! !IROUTINE: RUN2 -- Run2 method for GOCART2G component

! !INTERFACE:

 subroutine Run2 (GC, import, export, clock, RC)

! !ARGUMENTS:
    type (ESMF_GridComp), intent(inout) :: GC     ! Gridded component
    type (ESMF_State),    intent(inout) :: import ! Import state
    type (ESMF_State),    intent(inout) :: export ! Export state
    type (ESMF_Clock),    intent(inout) :: clock  ! The clock
    integer, optional,    intent(  out) :: RC     ! Error code:

! !DESCRIPTION: This version uses the MAPL\_GenericSetServices. This function sets
!                the Initialize and Finalize services, as well as allocating

!EOP
!============================================================================

!   Locals
    character(len=ESMF_MAXSTR)          :: COMP_NAME
    type (MAPL_MetaComp),      pointer  :: MAPL
    type (ESMF_GridComp),      pointer  :: gcs(:)
    type (ESMF_State),         pointer  :: gim(:)
    type (ESMF_State),         pointer  :: gex(:)
    type (ESMF_State)                   :: internal
    type (GOCART_State),       pointer  :: self

    type (wrap_)                        :: wrap
    character(len=ESMF_MAXSTR)          :: child_name
    integer                             :: i, n, w
    real, pointer, dimension(:,:)       :: LATS
    real, pointer, dimension(:,:)       :: LONS

    real, pointer, dimension(:,:,:) :: duexttau, dustexttau, &
                                       duscatau, dustscatau, &
                                       duextt25, duscat25, &
                                       duexttfm, duscatfm
    real, pointer, dimension(:,:,:,:) :: duextcoef, duscacoef
    real, pointer, dimension(:,:,:,:) :: duextcoefrh20, duextcoefrh80
    real, pointer, dimension(:,:,:,:) :: duscacoefrh20, duscacoefrh80
    real, pointer, dimension(:,:,:,:) :: dubckcoef
    real, pointer, dimension(:,:)   :: duangstr, dusmass,  &
                                       dusmass25
    real, pointer, dimension(:,:,:) :: ssexttau, ssstexttau, &
                                       ssscatau, ssstscatau, &
                                       ssextt25, ssscat25, &
                                       ssexttfm, ssscatfm
    real, pointer, dimension(:,:,:,:) :: ssextcoef, ssscacoef
    real, pointer, dimension(:,:,:,:) :: ssextcoefrh20, ssextcoefrh80
    real, pointer, dimension(:,:,:,:) :: ssscacoefrh20, ssscacoefrh80
    real, pointer, dimension(:,:,:,:) :: ssbckcoef
    real, pointer, dimension(:,:)   :: ssangstr, sssmass,  &
                                       sssmass25
    real, pointer, dimension(:,:,:) :: niexttau, nistexttau, &
                                       niscatau, nistscatau, &
                                       niextt25, niscat25, &
                                       niexttfm, niscatfm
    real, pointer, dimension(:,:,:,:) :: niextcoef, niscacoef
    real, pointer, dimension(:,:,:,:) :: niextcoefrh20, niextcoefrh80
    real, pointer, dimension(:,:,:,:) :: niscacoefrh20, niscacoefrh80
    real, pointer, dimension(:,:,:,:) :: nibckcoef
    real, pointer, dimension(:,:)   :: niangstr, nismass,  &
                                       nismass25
    real, pointer, dimension(:,:)   :: nh4smass
    real, pointer, dimension(:,:,:) :: suexttau, sustexttau, &
                                       suscatau, sustscatau
    real, pointer, dimension(:,:,:,:) :: suextcoef, suscacoef
    real, pointer, dimension(:,:,:,:) :: suextcoefrh20, suextcoefrh80
    real, pointer, dimension(:,:,:,:) :: suscacoefrh20, suscacoefrh80
    real, pointer, dimension(:,:,:,:) :: subckcoef
    real, pointer, dimension(:,:)   :: suangstr, so4smass
    real, pointer, dimension(:,:,:) :: bcexttau, bcstexttau, bcscatau, bcstscatau
    real, pointer, dimension(:,:,:,:) :: bcextcoef, bcscacoef
    real, pointer, dimension(:,:,:,:) :: bcextcoefrh20, bcextcoefrh80
    real, pointer, dimension(:,:,:,:) :: bcscacoefrh20, bcscacoefrh80
    real, pointer, dimension(:,:,:,:) :: bcbckcoef
    real, pointer, dimension(:,:)   :: bcangstr, bcsmass
    real, pointer, dimension(:,:,:) :: ocexttau, ocstexttau, ocscatau, ocstscatau
    real, pointer, dimension(:,:,:,:) :: ocextcoef, ocscacoef
    real, pointer, dimension(:,:,:,:) :: ocextcoefrh20, ocextcoefrh80
    real, pointer, dimension(:,:,:,:) :: ocscacoefrh20, ocscacoefrh80
    real, pointer, dimension(:,:,:,:) :: ocbckcoef
    real, pointer, dimension(:,:)   :: ocangstr, ocsmass
    real, pointer, dimension(:,:,:) :: brexttau, brstexttau, brscatau, brstscatau
    real, pointer, dimension(:,:,:,:) :: brextcoef, brscacoef
    real, pointer, dimension(:,:,:,:) :: brextcoefrh20, brextcoefrh80
    real, pointer, dimension(:,:,:,:) :: brscacoefrh20, brscacoefrh80
    real, pointer, dimension(:,:,:,:) :: brbckcoef
    real, pointer, dimension(:,:)   :: brangstr, brsmass
    real, pointer, dimension(:,:,:) :: pso4
    real, allocatable               :: tau1(:,:), tau2(:,:)
    real, allocatable               :: backscat_mol(:,:,:)
    real, allocatable               :: P(:,:,:), delz(:,:,:)
    real, allocatable               :: tau_mol_layer(:,:,:), tau_aer_layer(:,:,:)
    real, allocatable               :: tau_mol(:,:), tau_aer(:,:)
    real                            :: c1, c2, c3
    real                            :: nifactor
    real, parameter                 :: pi = 3.141529265
    integer                         :: ind550, ind532
    integer                         :: i1, i2, j1, j2, km, k,kk

#include "GOCART2G_DeclarePointer___.h"

    __Iam__('Run2')

!****************************************************************************
! Begin...

!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, __RC__ )
    Iam = trim(COMP_NAME)//'::'//Iam

!   Get my internal MAPL_Generic state
!   -----------------------------------
    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS )
    VERIFY_(STATUS)

!   Get parameters from generic state.
!   -----------------------------------
    call MAPL_Get ( MAPL, gcs=gcs, gim=gim, gex=gex, INTERNAL_ESMF_STATE=internal, &
                    LONS=LONS, LATS=LATS, __RC__ )

!   Get my internal state
!   ---------------------
    call ESMF_UserCompGetInternalState (GC, 'GOCART_State', wrap, STATUS)
    VERIFY_(STATUS)
    self => wrap%ptr

#include "GOCART2G_GetPointer___.h"

    if(associated(totexttau)) totexttau = 0.
    if(associated(totstexttau)) totstexttau = 0.
    if(associated(totscatau)) totscatau = 0.
    if(associated(totstscatau)) totstscatau = 0.
    if(associated(totextt25)) totextt25 = 0.
    if(associated(totscat25)) totscat25 = 0.
    if(associated(totexttfm)) totexttfm = 0.
    if(associated(totscatfm)) totscatfm = 0.
    if(associated(totextcoef))     totextcoef = 0.
    if(associated(totextcoefrh20)) totextcoefrh20 = 0.
    if(associated(totextcoefrh80)) totextcoefrh80 = 0.
    if(associated(totscacoef))     totscacoef = 0.
    if(associated(totscacoefrh20)) totscacoefrh20 = 0.
    if(associated(totscacoefrh80)) totscacoefrh80 = 0.
    if(associated(totbckcoef))     totbckcoef = 0.
    if(associated(totabcktoa))     totabcktoa = 0.
    if(associated(totabcksfc))     totabcksfc = 0.
    if(associated(pm))        pm(:,:)        = 0.
    if(associated(pm25))      pm25(:,:)      = 0.
    if(associated(pm_rh35))   pm_rh35(:,:)   = 0.
    if(associated(pm25_rh35)) pm25_rh35(:,:) = 0.
    if(associated(pm_rh50))   pm_rh50(:,:)   = 0.
    if(associated(pm25_rh50)) pm25_rh50(:,:) = 0.
    if(associated(pso4tot))   pso4tot(:,:,:) = 0.

!   Run the children
!   -----------------
    do i = 1, size(gcs)
      call ESMF_GridCompGet (gcs(i), NAME=child_name, __RC__ )
      if ((index(child_name, 'data')) == 0) then ! only execute Run2 method if a computational instance
         call ESMF_GridCompRun (gcs(i), importState=gim(i), exportState=gex(i), phase=2, clock=clock, __RC__)
      end if
    end do

!   Compute total aerosol diagnostic values for export
!   --------------------------------------------------
    if(associated(totangstr)) then
    ind550 = 0
       do w = 1, size(self%wavelengths_vertint) ! find index for 550nm to compute total angstrom
          if ((self%wavelengths_vertint(w)*1.e-9 .ge. 5.49e-7) .and. &
              (self%wavelengths_vertint(w)*1.e-9 .le. 5.51e-7)) then
             ind550 = w
             exit
          end if
       end do

       if (ind550 == 0) then
          !$omp critical (G2G_1)
          print*,trim(Iam),' : 550nm wavelengths is not present in GOCART2G_GridComp.rc.',&
                           ' Cannot produce TOTANGSTR variable without 550nm wavelength.'
          !$omp end critical (G2G_1)
          VERIFY_(100)
       end if

       totangstr = 0.0
       allocate(tau1(SIZE(LATS,1), SIZE(LATS,2)), &
                tau2(SIZE(LATS,1), SIZE(LATS,2)), __STAT__)

       tau1(:,:) = tiny(1.0)
       tau2(:,:) = tiny(1.0)
       c1 = -log(470./550.)
       c2 = -log(870./550.)
       c3 = -log(470./870.)
    end if


!   Dust
    do n = 1, size(self%DU%instances)
       if ((self%DU%instances(n)%is_active) .and. (index(self%DU%instances(n)%name, 'data') == 0 )) then
          call MAPL_GetPointer (gex(self%DU%instances(n)%id), duexttau, 'DUEXTTAU', __RC__)
          call MAPL_GetPointer (gex(self%DU%instances(n)%id), dustexttau, 'DUSTEXTTAU', __RC__)
          call MAPL_GetPointer (gex(self%DU%instances(n)%id), duscatau, 'DUSCATAU', __RC__)
          call MAPL_GetPointer (gex(self%DU%instances(n)%id), dustscatau, 'DUSTSCATAU', __RC__)
          call MAPL_GetPointer (gex(self%DU%instances(n)%id), duextcoef, 'DUEXTCOEF', __RC__)
          call MAPL_GetPointer (gex(self%DU%instances(n)%id), duextcoefrh20, 'DUEXTCOEFRH20', __RC__)
          call MAPL_GetPointer (gex(self%DU%instances(n)%id), duextcoefrh80, 'DUEXTCOEFRH80', __RC__)
          call MAPL_GetPointer (gex(self%DU%instances(n)%id), duscacoef, 'DUSCACOEF', __RC__)
          call MAPL_GetPointer (gex(self%DU%instances(n)%id), duscacoefrh20, 'DUSCACOEFRH20', __RC__)
          call MAPL_GetPointer (gex(self%DU%instances(n)%id), duscacoefrh80, 'DUSCACOEFRH80', __RC__)
          call MAPL_GetPointer (gex(self%DU%instances(n)%id), dubckcoef, 'DUBCKCOEF', __RC__)
          call MAPL_GetPointer (gex(self%DU%instances(n)%id), duextt25, 'DUEXTT25', __RC__)
          call MAPL_GetPointer (gex(self%DU%instances(n)%id), duscat25, 'DUSCAT25', __RC__)
          call MAPL_GetPointer (gex(self%DU%instances(n)%id), duexttfm, 'DUEXTTFM', __RC__)
          call MAPL_GetPointer (gex(self%DU%instances(n)%id), duscatfm, 'DUSCATFM', __RC__)
          call MAPL_GetPointer (gex(self%DU%instances(n)%id), duangstr, 'DUANGSTR', __RC__)

      !   Iterate over the wavelengths
          do w = 1, size(self%wavelengths_vertint)
             if(associated(totexttau) .and. associated(duexttau)) totexttau(:,:,w) = totexttau(:,:,w)+duexttau(:,:,w)
             if(associated(totstexttau) .and. associated(dustexttau)) totstexttau(:,:,w) = totstexttau(:,:,w)+dustexttau(:,:,w)
             if(associated(totscatau) .and. associated(duscatau)) totscatau(:,:,w) = totscatau(:,:,w)+duscatau(:,:,w)
             if(associated(totstscatau) .and. associated(dustscatau)) totstscatau(:,:,w) = totstscatau(:,:,w)+dustscatau(:,:,w)
             if(associated(totextt25) .and. associated(duextt25)) totextt25(:,:,w) = totextt25(:,:,w)+duextt25(:,:,w)
             if(associated(totscat25) .and. associated(duscat25)) totscat25(:,:,w) = totscat25(:,:,w)+duscat25(:,:,w)
             if(associated(totexttfm) .and. associated(duexttfm)) totexttfm(:,:,w) = totexttfm(:,:,w)+duexttfm(:,:,w)
             if(associated(totscatfm) .and. associated(duscatfm)) totscatfm(:,:,w) = totscatfm(:,:,w)+duscatfm(:,:,w)
          end do

          do w = 1, size(self%wavelengths_profile)
             if(associated(totextcoef) .and. associated(duextcoef)) totextcoef(:,:,:,w) = totextcoef(:,:,:,w)+duextcoef(:,:,:,w)
             if(associated(totextcoefrh20) .and. associated(duextcoefrh20)) totextcoefrh20(:,:,:,w) = totextcoefrh20(:,:,:,w)+duextcoefrh20(:,:,:,w)
             if(associated(totextcoefrh80) .and. associated(duextcoefrh80)) totextcoefrh80(:,:,:,w) = totextcoefrh80(:,:,:,w)+duextcoefrh80(:,:,:,w)
             if(associated(totscacoef) .and. associated(duscacoef)) totscacoef(:,:,:,w) = totscacoef(:,:,:,w)+duscacoef(:,:,:,w)
             if(associated(totscacoefrh20) .and. associated(duscacoefrh20)) totscacoefrh20(:,:,:,w) = totscacoefrh20(:,:,:,w)+duscacoefrh20(:,:,:,w)
             if(associated(totscacoefrh80) .and. associated(duscacoefrh80)) totscacoefrh80(:,:,:,w) = totscacoefrh80(:,:,:,w)+duscacoefrh80(:,:,:,w)
             if(associated(totbckcoef) .and. associated(dubckcoef)) totbckcoef(:,:,:,w) = totbckcoef(:,:,:,w)+dubckcoef(:,:,:,w)
          end do

          call MAPL_GetPointer (gex(self%DU%instances(n)%id), dusmass,   'DUSMASS',   __RC__)
          call MAPL_GetPointer (gex(self%DU%instances(n)%id), dusmass25, 'DUSMASS25', __RC__)
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

!   Sea Salt
    do n = 1, size(self%SS%instances)
       if ((self%SS%instances(n)%is_active) .and. (index(self%SS%instances(n)%name, 'data') == 0 )) then
          call MAPL_GetPointer (gex(self%SS%instances(n)%id), ssexttau, 'SSEXTTAU', __RC__)
          call MAPL_GetPointer (gex(self%SS%instances(n)%id), ssstexttau, 'SSSTEXTTAU', __RC__)
          call MAPL_GetPointer (gex(self%SS%instances(n)%id), ssscatau, 'SSSCATAU', __RC__)
          call MAPL_GetPointer (gex(self%SS%instances(n)%id), ssstscatau, 'SSSTSCATAU', __RC__)
          call MAPL_GetPointer (gex(self%SS%instances(n)%id), ssextcoef, 'SSEXTCOEF', __RC__)
          call MAPL_GetPointer (gex(self%SS%instances(n)%id), ssextcoefrh20, 'SSEXTCOEFRH20', __RC__)
          call MAPL_GetPointer (gex(self%SS%instances(n)%id), ssextcoefrh80, 'SSEXTCOEFRH80', __RC__)
          call MAPL_GetPointer (gex(self%SS%instances(n)%id), ssscacoef, 'SSSCACOEF', __RC__)
          call MAPL_GetPointer (gex(self%SS%instances(n)%id), ssscacoefrh20, 'SSSCACOEFRH20', __RC__)
          call MAPL_GetPointer (gex(self%SS%instances(n)%id), ssscacoefrh80, 'SSSCACOEFRH80', __RC__)
          call MAPL_GetPointer (gex(self%SS%instances(n)%id), ssbckcoef, 'SSBCKCOEF', __RC__)
          call MAPL_GetPointer (gex(self%SS%instances(n)%id), ssextt25, 'SSEXTT25', __RC__)
          call MAPL_GetPointer (gex(self%SS%instances(n)%id), ssscat25, 'SSSCAT25', __RC__)
          call MAPL_GetPointer (gex(self%SS%instances(n)%id), ssexttfm, 'SSEXTTFM', __RC__)
          call MAPL_GetPointer (gex(self%SS%instances(n)%id), ssscatfm, 'SSSCATFM', __RC__)
          call MAPL_GetPointer (gex(self%SS%instances(n)%id), ssangstr, 'SSANGSTR', __RC__)

      !   Iterate over the wavelengths
          do w = 1, size(self%wavelengths_vertint)
             if(associated(totexttau) .and. associated(ssexttau)) totexttau(:,:,w) = totexttau(:,:,w)+ssexttau(:,:,w)
             if(associated(totstexttau) .and. associated(ssstexttau)) totstexttau(:,:,w) = totstexttau(:,:,w)+ssstexttau(:,:,w)
             if(associated(totscatau) .and. associated(ssscatau)) totscatau(:,:,w) = totscatau(:,:,w)+ssscatau(:,:,w)
             if(associated(totstscatau) .and. associated(ssstscatau)) totstscatau(:,:,w) = totstscatau(:,:,w)+ssstscatau(:,:,w)
             if(associated(totextt25) .and. associated(ssextt25)) totextt25(:,:,w) = totextt25(:,:,w)+ssextt25(:,:,w)
             if(associated(totscat25) .and. associated(ssscat25)) totscat25(:,:,w) = totscat25(:,:,w)+ssscat25(:,:,w)
             if(associated(totexttfm) .and. associated(ssexttfm)) totexttfm(:,:,w) = totexttfm(:,:,w)+ssexttfm(:,:,w)
             if(associated(totscatfm) .and. associated(ssscatfm)) totscatfm(:,:,w) = totscatfm(:,:,w)+ssscatfm(:,:,w)
          end do

          do w = 1, size(self%wavelengths_profile)
             if(associated(totextcoef) .and. associated(ssextcoef)) totextcoef(:,:,:,w) = totextcoef(:,:,:,w)+ssextcoef(:,:,:,w)
             if(associated(totextcoefrh20) .and. associated(ssextcoefrh20)) totextcoefrh20(:,:,:,w) = totextcoefrh20(:,:,:,w)+ssextcoefrh20(:,:,:,w)
             if(associated(totextcoefrh80) .and. associated(ssextcoefrh80)) totextcoefrh80(:,:,:,w) = totextcoefrh80(:,:,:,w)+ssextcoefrh80(:,:,:,w)
             if(associated(totscacoef) .and. associated(ssscacoef)) totscacoef(:,:,:,w) = totscacoef(:,:,:,w)+ssscacoef(:,:,:,w)
             if(associated(totscacoefrh20) .and. associated(ssscacoefrh20)) totscacoefrh20(:,:,:,w) = totscacoefrh20(:,:,:,w)+ssscacoefrh20(:,:,:,w)
             if(associated(totscacoefrh80) .and. associated(ssscacoefrh80)) totscacoefrh80(:,:,:,w) = totscacoefrh80(:,:,:,w)+ssscacoefrh80(:,:,:,w)
             if(associated(totbckcoef) .and. associated(ssbckcoef)) totbckcoef(:,:,:,w) = totbckcoef(:,:,:,w)+ssbckcoef(:,:,:,w)
          enddo

          call MAPL_GetPointer (gex(self%SS%instances(n)%id), sssmass,   'SSSMASS',   __RC__)
          call MAPL_GetPointer (gex(self%SS%instances(n)%id), sssmass25, 'SSSMASS25', __RC__)
          if(associated(pm)        .and. associated(sssmass))   pm        = pm        + sssmass
          if(associated(pm25)      .and. associated(sssmass25)) pm25      = pm25      + sssmass25
          if(associated(pm_rh35)   .and. associated(sssmass))   pm_rh35   = pm_rh35   + 1.86*sssmass
          if(associated(pm25_rh35) .and. associated(sssmass25)) pm25_rh35 = pm25_rh35 + 1.86*sssmass25
          if(associated(pm_rh50)   .and. associated(sssmass))   pm_rh50   = pm_rh50   + 2.42*sssmass
          if(associated(pm25_rh50) .and. associated(sssmass25)) pm25_rh50 = pm25_rh50 + 2.42*sssmass25

          if(associated(totangstr) .and. associated(ssexttau) .and. associated(ssangstr)) then
             tau1 = tau1 + ssexttau(:,:,ind550)*exp(c1*ssangstr)
             tau2 = tau2 + ssexttau(:,:,ind550)*exp(c2*ssangstr)
          end if
       end if
    end do

!   Nitrates - NOTE! Nitrates currently only support one active instance
    do n = 1, size(self%NI%instances)
       if ((self%NI%instances(n)%is_active) .and. (index(self%NI%instances(n)%name, 'data') == 0 )) then
          call MAPL_GetPointer (gex(self%NI%instances(n)%id), niexttau, 'NIEXTTAU', __RC__)
          call MAPL_GetPointer (gex(self%NI%instances(n)%id), nistexttau, 'NISTEXTTAU', __RC__)
          call MAPL_GetPointer (gex(self%NI%instances(n)%id), niscatau, 'NISCATAU', __RC__)
          call MAPL_GetPointer (gex(self%NI%instances(n)%id), nistscatau, 'NISTSCATAU', __RC__)
          call MAPL_GetPointer (gex(self%NI%instances(n)%id), niextcoef, 'NIEXTCOEF', __RC__)
          call MAPL_GetPointer (gex(self%NI%instances(n)%id), niextcoefrh20, 'NIEXTCOEFRH20', __RC__)
          call MAPL_GetPointer (gex(self%NI%instances(n)%id), niextcoefrh80, 'NIEXTCOEFRH80', __RC__)
          call MAPL_GetPointer (gex(self%NI%instances(n)%id), niscacoef, 'NISCACOEF', __RC__)
          call MAPL_GetPointer (gex(self%NI%instances(n)%id), niscacoefrh20, 'NISCACOEFRH20', __RC__)
          call MAPL_GetPointer (gex(self%NI%instances(n)%id), niscacoefrh80, 'NISCACOEFRH80', __RC__)
          call MAPL_GetPointer (gex(self%NI%instances(n)%id), nibckcoef, 'NIBCKCOEF', __RC__)
          call MAPL_GetPointer (gex(self%NI%instances(n)%id), niextt25, 'NIEXTT25', __RC__)
          call MAPL_GetPointer (gex(self%NI%instances(n)%id), niscat25, 'NISCAT25', __RC__)
          call MAPL_GetPointer (gex(self%NI%instances(n)%id), niexttfm, 'NIEXTTFM', __RC__)
          call MAPL_GetPointer (gex(self%NI%instances(n)%id), niscatfm, 'NISCATFM', __RC__)
          call MAPL_GetPointer (gex(self%NI%instances(n)%id), niangstr, 'NIANGSTR', __RC__)

      !   Iterate over the wavelengths
          do w = 1, size(self%wavelengths_vertint)
             if(associated(totexttau) .and. associated(niexttau)) totexttau(:,:,w) = totexttau(:,:,w)+niexttau(:,:,w)
             if(associated(totstexttau) .and. associated(nistexttau)) totstexttau(:,:,w) = totstexttau(:,:,w)+nistexttau(:,:,w)
             if(associated(totscatau) .and. associated(niscatau)) totscatau(:,:,w) = totscatau(:,:,w)+niscatau(:,:,w)
             if(associated(totstscatau) .and. associated(nistscatau)) totstscatau(:,:,w) = totstscatau(:,:,w)+nistscatau(:,:,w)
             if(associated(totextt25) .and. associated(niextt25)) totextt25(:,:,w) = totextt25(:,:,w)+niextt25(:,:,w)
             if(associated(totscat25) .and. associated(niscat25)) totscat25(:,:,w) = totscat25(:,:,w)+niscat25(:,:,w)
             if(associated(totexttfm) .and. associated(niexttfm)) totexttfm(:,:,w) = totexttfm(:,:,w)+niexttfm(:,:,w)
             if(associated(totscatfm) .and. associated(niscatfm)) totscatfm(:,:,w) = totscatfm(:,:,w)+niscatfm(:,:,w)
          end do

          do w = 1, size(self%wavelengths_profile)
             if(associated(totextcoef) .and. associated(niextcoef)) totextcoef(:,:,:,w) = totextcoef(:,:,:,w)+niextcoef(:,:,:,w)
             if(associated(totextcoefrh20) .and. associated(niextcoefrh20)) totextcoefrh20(:,:,:,w) = totextcoefrh20(:,:,:,w)+niextcoefrh20(:,:,:,w)
             if(associated(totextcoefrh80) .and. associated(niextcoefrh80)) totextcoefrh80(:,:,:,w) = totextcoefrh80(:,:,:,w)+niextcoefrh80(:,:,:,w)
             if(associated(totscacoef) .and. associated(niscacoef)) totscacoef(:,:,:,w) = totscacoef(:,:,:,w)+niscacoef(:,:,:,w)
             if(associated(totscacoefrh20) .and. associated(niscacoefrh20)) totscacoefrh20(:,:,:,w) = totscacoefrh20(:,:,:,w)+niscacoefrh20(:,:,:,w)
             if(associated(totscacoefrh80) .and. associated(niscacoefrh80)) totscacoefrh80(:,:,:,w) = totscacoefrh80(:,:,:,w)+niscacoefrh80(:,:,:,w)
             if(associated(totbckcoef) .and. associated(nibckcoef)) totbckcoef(:,:,:,w) = totbckcoef(:,:,:,w)+nibckcoef(:,:,:,w)
          end do

          call MAPL_GetPointer (gex(self%NI%instances(n)%id), nismass,   'NISMASS',   __RC__)
          call MAPL_GetPointer (gex(self%NI%instances(n)%id), nismass25, 'NISMASS25', __RC__)
          call MAPL_GetPointer (gex(self%NI%instances(n)%id), nh4smass,  'NH4SMASS',   __RC__)
          if(associated(pm)        .and. associated(nismass)   .and. associated(nh4smass)) pm        = pm   + nismass   + nh4smass
          if(associated(pm25)      .and. associated(nismass25) .and. associated(nh4smass)) pm25      = pm25 + nismass25 + nh4smass
          if(associated(pm_rh35)   .and. associated(nismass)   .and. associated(nh4smass)) pm_rh35   = pm_rh35   + 1.33*(nismass   + nh4smass)
          if(associated(pm25_rh35) .and. associated(nismass25) .and. associated(nh4smass)) pm25_rh35 = pm25_rh35 + 1.33*(nismass25 + nh4smass)
          if(associated(pm_rh50)   .and. associated(nismass)   .and. associated(nh4smass)) pm_rh50   = pm_rh50   + 1.51*(nismass   + nh4smass)
          if(associated(pm25_rh50) .and. associated(nismass25) .and. associated(nh4smass)) pm25_rh50 = pm25_rh50 + 1.51*(nismass25 + nh4smass)

          if(associated(totangstr) .and. associated(niexttau) .and. associated(niangstr)) then
             tau1 = tau1 + niexttau(:,:,ind550)*exp(c1*niangstr)
             tau2 = tau2 + niexttau(:,:,ind550)*exp(c2*niangstr)
          end if
       end if
    end do

!   Sulfates
    nifactor = 132.14/96.06
    if (size(self%NI%instances) > 0) then
      if ((self%NI%instances(1)%is_active) .and. (index(self%NI%instances(1)%name, 'data') == 0 )) nifactor = 1.0
    end if

    do n = 1, size(self%SU%instances)
       if ((self%SU%instances(n)%is_active) .and. (index(self%SU%instances(n)%name, 'data') == 0 )) then
          call MAPL_GetPointer (gex(self%SU%instances(n)%id), suexttau, 'SUEXTTAU', __RC__)
          call MAPL_GetPointer (gex(self%SU%instances(n)%id), suextcoef, 'SUEXTCOEF', __RC__)
          call MAPL_GetPointer (gex(self%SU%instances(n)%id), suextcoefrh20, 'SUEXTCOEFRH20', __RC__)
          call MAPL_GetPointer (gex(self%SU%instances(n)%id), suextcoefrh80, 'SUEXTCOEFRH80', __RC__)
          call MAPL_GetPointer (gex(self%SU%instances(n)%id), suscacoef, 'SUSCACOEF', __RC__)
          call MAPL_GetPointer (gex(self%SU%instances(n)%id), suscacoefrh20, 'SUSCACOEFRH20', __RC__)
          call MAPL_GetPointer (gex(self%SU%instances(n)%id), suscacoefrh80, 'SUSCACOEFRH80', __RC__)
          call MAPL_GetPointer (gex(self%SU%instances(n)%id), subckcoef, 'SUBCKCOEF', __RC__)
          call MAPL_GetPointer (gex(self%SU%instances(n)%id), sustexttau, 'SUSTEXTTAU', __RC__)
          call MAPL_GetPointer (gex(self%SU%instances(n)%id), suscatau, 'SUSCATAU', __RC__)
          call MAPL_GetPointer (gex(self%SU%instances(n)%id), sustscatau, 'SUSTSCATAU', __RC__)
          call MAPL_GetPointer (gex(self%SU%instances(n)%id), suangstr, 'SUANGSTR', __RC__)

          !   Iterate over the wavelengths
          do w = 1, size(self%wavelengths_vertint)
             if(associated(totexttau) .and. associated(suexttau)) totexttau(:,:,w) = totexttau(:,:,w)+suexttau(:,:,w)
             if(associated(totstexttau) .and. associated(sustexttau)) totstexttau(:,:,w) = totstexttau(:,:,w)+sustexttau(:,:,w)
             if(associated(totscatau) .and. associated(suscatau)) totscatau(:,:,w) = totscatau(:,:,w)+suscatau(:,:,w)
             if(associated(totstscatau) .and. associated(sustscatau)) totstscatau(:,:,w) = totstscatau(:,:,w)+sustscatau(:,:,w)
             if(associated(totextt25) .and. associated(suexttau)) totextt25(:,:,w) = totextt25(:,:,w)+suexttau(:,:,w)
             if(associated(totscat25) .and. associated(suscatau)) totscat25(:,:,w) = totscat25(:,:,w)+suscatau(:,:,w)
             if(associated(totexttfm) .and. associated(suexttau)) totexttfm(:,:,w) = totexttfm(:,:,w)+suexttau(:,:,w)
             if(associated(totscatfm) .and. associated(suscatau)) totscatfm(:,:,w) = totscatfm(:,:,w)+suscatau(:,:,w)
          end do

          do w = 1, size(self%wavelengths_profile)
             if(associated(totextcoef) .and. associated(suextcoef)) totextcoef(:,:,:,w) = totextcoef(:,:,:,w)+suextcoef(:,:,:,w)
             if(associated(totextcoefrh20) .and. associated(suextcoefrh20)) totextcoefrh20(:,:,:,w) = totextcoefrh20(:,:,:,w)+suextcoefrh20(:,:,:,w)
             if(associated(totextcoefrh80) .and. associated(suextcoefrh80)) totextcoefrh80(:,:,:,w) = totextcoefrh80(:,:,:,w)+suextcoefrh80(:,:,:,w)
             if(associated(totscacoef) .and. associated(suscacoef)) totscacoef(:,:,:,w) = totscacoef(:,:,:,w)+suscacoef(:,:,:,w)
             if(associated(totscacoefrh20) .and. associated(suscacoefrh20)) totscacoefrh20(:,:,:,w) = totscacoefrh20(:,:,:,w)+suscacoefrh20(:,:,:,w)
             if(associated(totscacoefrh80) .and. associated(suscacoefrh80)) totscacoefrh80(:,:,:,w) = totscacoefrh80(:,:,:,w)+suscacoefrh80(:,:,:,w)
             if(associated(totbckcoef) .and. associated(subckcoef)) totbckcoef(:,:,:,w) = totbckcoef(:,:,:,w)+subckcoef(:,:,:,w)
          end do

          call MAPL_GetPointer (gex(self%SU%instances(n)%id), pso4, 'PSO4', __RC__)
          if(associated(pso4tot) .and. associated(pso4)) pso4tot = pso4tot + pso4

          call MAPL_GetPointer (gex(self%SU%instances(n)%id), so4smass, 'SO4SMASS', __RC__)
          if(associated(so4smass)) then
             if(associated(pm)       ) pm        = pm        + nifactor*so4smass
             if(associated(pm25)     ) pm25      = pm25      + nifactor*so4smass
             if(associated(pm_rh35)  ) pm_rh35   = pm_rh35   + 1.33*nifactor*so4smass
             if(associated(pm25_rh35)) pm25_rh35 = pm25_rh35 + 1.33*nifactor*so4smass
             if(associated(pm_rh50)  ) pm_rh50   = pm_rh50   + 1.51*nifactor*so4smass
             if(associated(pm25_rh50)) pm25_rh50 = pm25_rh50 + 1.51*nifactor*so4smass
          end if

          if(associated(totangstr) .and. associated(suexttau) .and. associated(suangstr)) then
             tau1 = tau1 + suexttau(:,:,ind550)*exp(c1*suangstr)
             tau2 = tau2 + suexttau(:,:,ind550)*exp(c2*suangstr)
          end if
       end if
    end do


!   Carbonaceous aerosols
    do n = 1, size(self%CA%instances)
       if ((self%CA%instances(n)%is_active) .and. (index(self%CA%instances(n)%name, 'data') == 0 ) &
           .and. (index(self%CA%instances(n)%name, 'CA.bc') > 0)) then

          call MAPL_GetPointer (gex(self%CA%instances(n)%id), bcexttau, 'CA.bcEXTTAU', __RC__)
          call MAPL_GetPointer (gex(self%CA%instances(n)%id), bcstexttau, 'CA.bcSTEXTTAU', __RC__)
          call MAPL_GetPointer (gex(self%CA%instances(n)%id), bcscatau, 'CA.bcSCATAU', __RC__)
          call MAPL_GetPointer (gex(self%CA%instances(n)%id), bcstscatau, 'CA.bcSTSCATAU', __RC__)
          call MAPL_GetPointer (gex(self%CA%instances(n)%id), bcangstr, 'CA.bcANGSTR', __RC__)
          call MAPL_GetPointer (gex(self%CA%instances(n)%id), bcextcoef, 'CA.bcEXTCOEF', __RC__)
          call MAPL_GetPointer (gex(self%CA%instances(n)%id), bcextcoefrh20, 'CA.bcEXTCOEFRH20', __RC__)
          call MAPL_GetPointer (gex(self%CA%instances(n)%id), bcextcoefrh80, 'CA.bcEXTCOEFRH80', __RC__)
          call MAPL_GetPointer (gex(self%CA%instances(n)%id), bcscacoef, 'CA.bcSCACOEF', __RC__)
          call MAPL_GetPointer (gex(self%CA%instances(n)%id), bcscacoefrh20, 'CA.bcSCACOEFRH20', __RC__)
          call MAPL_GetPointer (gex(self%CA%instances(n)%id), bcscacoefrh80, 'CA.bcSCACOEFRH80', __RC__)
          call MAPL_GetPointer (gex(self%CA%instances(n)%id), bcbckcoef, 'CA.bcBCKCOEF', __RC__)

          !   Iterate over the wavelengths
          do w = 1, size(self%wavelengths_vertint)
             if(associated(totexttau) .and. associated(bcexttau)) totexttau(:,:,w) = totexttau(:,:,w)+bcexttau(:,:,w)
             if(associated(totstexttau) .and. associated(bcstexttau)) totstexttau(:,:,w) = totstexttau(:,:,w)+bcstexttau(:,:,w)
             if(associated(totscatau) .and. associated(bcscatau)) totscatau(:,:,w) = totscatau(:,:,w)+bcscatau(:,:,w)
             if(associated(totstscatau) .and. associated(bcstscatau)) totstscatau(:,:,w) = totstscatau(:,:,w)+bcstscatau(:,:,w)
             if(associated(totextt25) .and. associated(bcexttau)) totextt25(:,:,w) = totextt25(:,:,w)+bcexttau(:,:,w)
             if(associated(totscat25) .and. associated(bcscatau)) totscat25(:,:,w) = totscat25(:,:,w)+bcscatau(:,:,w)
             if(associated(totexttfm) .and. associated(bcexttau)) totexttfm(:,:,w) = totexttfm(:,:,w)+bcexttau(:,:,w)
             if(associated(totscatfm) .and. associated(bcscatau)) totscatfm(:,:,w) = totscatfm(:,:,w)+bcscatau(:,:,w)
          end do

          do w = 1, size(self%wavelengths_profile)
             if(associated(totextcoef) .and. associated(bcextcoef)) totextcoef(:,:,:,w) = totextcoef(:,:,:,w)+bcextcoef(:,:,:,w)
             if(associated(totextcoefrh20) .and. associated(bcextcoefrh20)) totextcoefrh20(:,:,:,w) = totextcoefrh20(:,:,:,w)+bcextcoefrh20(:,:,:,w)
             if(associated(totextcoefrh80) .and. associated(bcextcoefrh80)) totextcoefrh80(:,:,:,w) = totextcoefrh80(:,:,:,w)+bcextcoefrh80(:,:,:,w)
             if(associated(totscacoef) .and. associated(bcscacoef)) totscacoef(:,:,:,w) = totscacoef(:,:,:,w)+bcscacoef(:,:,:,w)
             if(associated(totscacoefrh20) .and. associated(bcscacoefrh20)) totscacoefrh20(:,:,:,w) = totscacoefrh20(:,:,:,w)+bcscacoefrh20(:,:,:,w)
             if(associated(totscacoefrh80) .and. associated(bcscacoefrh80)) totscacoefrh80(:,:,:,w) = totscacoefrh80(:,:,:,w)+bcscacoefrh80(:,:,:,w)
             if(associated(totbckcoef) .and. associated(bcbckcoef)) totbckcoef(:,:,:,w) = totbckcoef(:,:,:,w)+bcbckcoef(:,:,:,w)
          end do

          call MAPL_GetPointer (gex(self%CA%instances(n)%id), bcsmass, 'CA.bcSMASS', __RC__)
          if(associated(pm)        .and. associated(bcsmass)) pm        = pm        + bcsmass
          if(associated(pm25)      .and. associated(bcsmass)) pm25      = pm25      + bcsmass
          if(associated(pm_rh35)   .and. associated(bcsmass)) pm_rh35   = pm_rh35   + bcsmass
          if(associated(pm25_rh35) .and. associated(bcsmass)) pm25_rh35 = pm25_rh35 + bcsmass
          if(associated(pm_rh50)   .and. associated(bcsmass)) pm_rh50   = pm_rh50   + bcsmass
          if(associated(pm25_rh50) .and. associated(bcsmass)) pm25_rh50 = pm25_rh50 + bcsmass

          if(associated(totangstr) .and. associated(bcexttau) .and. associated(bcangstr)) then
             tau1 = tau1 + bcexttau(:,:,ind550)*exp(c1*bcangstr)
             tau2 = tau2 + bcexttau(:,:,ind550)*exp(c2*bcangstr)
          end if

       else if ((self%CA%instances(n)%is_active) .and. (index(self%CA%instances(n)%name, 'data') == 0 ) &
                .and. (index(self%CA%instances(n)%name, 'CA.oc') > 0)) then
          call MAPL_GetPointer (gex(self%CA%instances(n)%id), ocexttau, 'CA.ocEXTTAU', __RC__)
          call MAPL_GetPointer (gex(self%CA%instances(n)%id), ocstexttau, 'CA.ocSTEXTTAU', __RC__)
          call MAPL_GetPointer (gex(self%CA%instances(n)%id), ocscatau, 'CA.ocSCATAU', __RC__)
          call MAPL_GetPointer (gex(self%CA%instances(n)%id), ocstscatau, 'CA.ocSTSCATAU', __RC__)
          call MAPL_GetPointer (gex(self%CA%instances(n)%id), ocangstr, 'CA.ocANGSTR', __RC__)
          call MAPL_GetPointer (gex(self%CA%instances(n)%id), ocextcoef, 'CA.ocEXTCOEF', __RC__)
          call MAPL_GetPointer (gex(self%CA%instances(n)%id), ocextcoefrh20, 'CA.ocEXTCOEFRH20', __RC__)
          call MAPL_GetPointer (gex(self%CA%instances(n)%id), ocextcoefrh80, 'CA.ocEXTCOEFRH80', __RC__)
          call MAPL_GetPointer (gex(self%CA%instances(n)%id), ocscacoef, 'CA.ocSCACOEF', __RC__)
          call MAPL_GetPointer (gex(self%CA%instances(n)%id), ocscacoefrh20, 'CA.ocSCACOEFRH20', __RC__)
          call MAPL_GetPointer (gex(self%CA%instances(n)%id), ocscacoefrh80, 'CA.ocSCACOEFRH80', __RC__)
          call MAPL_GetPointer (gex(self%CA%instances(n)%id), ocbckcoef, 'CA.ocBCKCOEF', __RC__)

          !   Iterate over the wavelengths
          do w = 1, size(self%wavelengths_vertint)
             if(associated(totexttau) .and. associated(ocexttau)) totexttau(:,:,w) = totexttau(:,:,w)+ocexttau(:,:,w)
             if(associated(totstexttau) .and. associated(ocstexttau)) totstexttau(:,:,w) = totstexttau(:,:,w)+ocstexttau(:,:,w)
             if(associated(totscatau) .and. associated(ocscatau)) totscatau(:,:,w) = totscatau(:,:,w)+ocscatau(:,:,w)
             if(associated(totstscatau) .and. associated(ocstscatau)) totstscatau(:,:,w) = totstscatau(:,:,w)+ocstscatau(:,:,w)
             if(associated(totextt25) .and. associated(ocexttau)) totextt25(:,:,w) = totextt25(:,:,w)+ocexttau(:,:,w)
             if(associated(totscat25) .and. associated(ocscatau)) totscat25(:,:,w) = totscat25(:,:,w)+ocscatau(:,:,w)
             if(associated(totexttfm) .and. associated(ocexttau)) totexttfm(:,:,w) = totexttfm(:,:,w)+ocexttau(:,:,w)
             if(associated(totscatfm) .and. associated(ocscatau)) totscatfm(:,:,w) = totscatfm(:,:,w)+ocscatau(:,:,w)
          end do

          do w = 1, size(self%wavelengths_profile)
             if(associated(totextcoef) .and. associated(ocextcoef)) totextcoef(:,:,:,w) = totextcoef(:,:,:,w)+ocextcoef(:,:,:,w)
             if(associated(totextcoefrh20) .and. associated(ocextcoefrh20)) totextcoefrh20(:,:,:,w) = totextcoefrh20(:,:,:,w)+ocextcoefrh20(:,:,:,w)
             if(associated(totextcoefrh80) .and. associated(ocextcoefrh80)) totextcoefrh80(:,:,:,w) = totextcoefrh80(:,:,:,w)+ocextcoefrh80(:,:,:,w)
             if(associated(totscacoef) .and. associated(ocscacoef)) totscacoef(:,:,:,w) = totscacoef(:,:,:,w)+ocscacoef(:,:,:,w)
             if(associated(totscacoefrh20) .and. associated(ocscacoefrh20)) totscacoefrh20(:,:,:,w) = totscacoefrh20(:,:,:,w)+ocscacoefrh20(:,:,:,w)
             if(associated(totscacoefrh80) .and. associated(ocscacoefrh80)) totscacoefrh80(:,:,:,w) = totscacoefrh80(:,:,:,w)+ocscacoefrh80(:,:,:,w)
             if(associated(totbckcoef) .and. associated(ocbckcoef)) totbckcoef(:,:,:,w) = totbckcoef(:,:,:,w)+ocbckcoef(:,:,:,w)
          end do

          call MAPL_GetPointer (gex(self%CA%instances(n)%id), ocsmass, 'CA.ocSMASS', __RC__)
          if(associated(pm)        .and. associated(ocsmass)) pm        = pm        + ocsmass
          if(associated(pm25)      .and. associated(ocsmass)) pm25      = pm25      + ocsmass
          if(associated(pm_rh35)   .and. associated(ocsmass)) pm_rh35   = pm_rh35   + 1.16*ocsmass  ! needs to be revisited: OCpho + 1.16 OCphi
          if(associated(pm25_rh35) .and. associated(ocsmass)) pm25_rh35 = pm25_rh35 + 1.16*ocsmass  !
          if(associated(pm_rh50)   .and. associated(ocsmass)) pm_rh50   = pm_rh50   + 1.24*ocsmass  ! needs to be revisited: OCpho + 1.24 OCphi
          if(associated(pm25_rh50) .and. associated(ocsmass)) pm25_rh50 = pm25_rh50 + 1.24*ocsmass  !

          if(associated(totangstr) .and. associated(ocexttau) .and. associated(ocangstr)) then
             tau1 = tau1 + ocexttau(:,:,ind550)*exp(c1*ocangstr)
             tau2 = tau2 + ocexttau(:,:,ind550)*exp(c2*ocangstr)
          end if

       else if ((self%CA%instances(n)%is_active) .and. (index(self%CA%instances(n)%name, 'data') == 0 ) &
                .and. (index(self%CA%instances(n)%name, 'CA.br') > 0)) then
          call MAPL_GetPointer (gex(self%CA%instances(n)%id), brexttau, 'CA.brEXTTAU', __RC__)
          call MAPL_GetPointer (gex(self%CA%instances(n)%id), brstexttau, 'CA.brSTEXTTAU', __RC__)
          call MAPL_GetPointer (gex(self%CA%instances(n)%id), brscatau, 'CA.brSCATAU', __RC__)
          call MAPL_GetPointer (gex(self%CA%instances(n)%id), brstscatau, 'CA.brSTSCATAU', __RC__)
          call MAPL_GetPointer (gex(self%CA%instances(n)%id), brangstr, 'CA.brANGSTR', __RC__)
          call MAPL_GetPointer (gex(self%CA%instances(n)%id), brextcoef, 'CA.brEXTCOEF', __RC__)
          call MAPL_GetPointer (gex(self%CA%instances(n)%id), brextcoefrh20, 'CA.brEXTCOEFRH20', __RC__)
          call MAPL_GetPointer (gex(self%CA%instances(n)%id), brextcoefrh80, 'CA.brEXTCOEFRH80', __RC__)
          call MAPL_GetPointer (gex(self%CA%instances(n)%id), brscacoef, 'CA.brSCACOEF', __RC__)
          call MAPL_GetPointer (gex(self%CA%instances(n)%id), brscacoefrh20, 'CA.brSCACOEFRH20', __RC__)
          call MAPL_GetPointer (gex(self%CA%instances(n)%id), brscacoefrh80, 'CA.brSCACOEFRH80', __RC__)
          call MAPL_GetPointer (gex(self%CA%instances(n)%id), brbckcoef, 'CA.brBCKCOEF', __RC__)

          !   Iterate over the wavelengths
          do w = 1, size(self%wavelengths_vertint)
             if(associated(totexttau) .and. associated(brexttau)) totexttau(:,:,w) = totexttau(:,:,w)+brexttau(:,:,w)
             if(associated(totstexttau) .and. associated(brstexttau)) totstexttau(:,:,w) = totstexttau(:,:,w)+brstexttau(:,:,w)
             if(associated(totscatau) .and. associated(brscatau)) totscatau(:,:,w) = totscatau(:,:,w)+brscatau(:,:,w)
             if(associated(totstscatau) .and. associated(brstscatau)) totstscatau(:,:,w) = totstscatau(:,:,w)+brstscatau(:,:,w)
             if(associated(totextt25) .and. associated(brexttau)) totextt25(:,:,w) = totextt25(:,:,w)+brexttau(:,:,w)
             if(associated(totscat25) .and. associated(brscatau)) totscat25(:,:,w) = totscat25(:,:,w)+brscatau(:,:,w)
             if(associated(totexttfm) .and. associated(brexttau)) totexttfm(:,:,w) = totexttfm(:,:,w)+brexttau(:,:,w)
             if(associated(totscatfm) .and. associated(brscatau)) totscatfm(:,:,w) = totscatfm(:,:,w)+brscatau(:,:,w)
          end do

          do w = 1, size(self%wavelengths_profile)
             if(associated(totextcoef) .and. associated(brextcoef)) totextcoef(:,:,:,w) = totextcoef(:,:,:,w)+brextcoef(:,:,:,w)
             if(associated(totextcoefrh20) .and. associated(brextcoefrh20)) totextcoefrh20(:,:,:,w) = totextcoefrh20(:,:,:,w)+brextcoefrh20(:,:,:,w)
             if(associated(totextcoefrh80) .and. associated(brextcoefrh80)) totextcoefrh80(:,:,:,w) = totextcoefrh80(:,:,:,w)+brextcoefrh80(:,:,:,w)
             if(associated(totscacoef) .and. associated(brscacoef)) totscacoef(:,:,:,w) = totscacoef(:,:,:,w)+brscacoef(:,:,:,w)
             if(associated(totscacoefrh20) .and. associated(brscacoefrh20)) totscacoefrh20(:,:,:,w) = totscacoefrh20(:,:,:,w)+brscacoefrh20(:,:,:,w)
             if(associated(totscacoefrh80) .and. associated(brscacoefrh80)) totscacoefrh80(:,:,:,w) = totscacoefrh80(:,:,:,w)+brscacoefrh80(:,:,:,w)
             if(associated(totbckcoef) .and. associated(brbckcoef)) totbckcoef(:,:,:,w) = totbckcoef(:,:,:,w)+brbckcoef(:,:,:,w)
          end do

          call MAPL_GetPointer (gex(self%CA%instances(n)%id), brsmass, 'CA.brSMASS', __RC__)
          if(associated(pm)        .and. associated(brsmass)) pm        = pm        + brsmass
          if(associated(pm25)      .and. associated(brsmass)) pm25      = pm25      + brsmass
          if(associated(pm_rh35)   .and. associated(brsmass)) pm_rh35   = pm_rh35   + 1.16*brsmass  ! needs to be revisited: OCpho + 1.16 OCphi
          if(associated(pm25_rh35) .and. associated(brsmass)) pm25_rh35 = pm25_rh35 + 1.16*brsmass  !
          if(associated(pm_rh50)   .and. associated(brsmass)) pm_rh50   = pm_rh50   + 1.24*brsmass  ! needs to be revisited: OCpho + 1.24 OCphi
          if(associated(pm25_rh50) .and. associated(brsmass)) pm25_rh50 = pm25_rh50 + 1.24*brsmass  !

          if(associated(totangstr) .and. associated(brexttau) .and. associated(brangstr)) then
             tau1 = tau1 + brexttau(:,:,ind550)*exp(c1*brangstr)
             tau2 = tau2 + brexttau(:,:,ind550)*exp(c2*brangstr)
          end if
       end if
    end do

!   Finish calculating totangstr
    if(associated(totangstr)) then
       totangstr = log(tau1/tau2)/c3
    end if

!  Calculate the total (molecular + aer) single scattering attenuated backscater coef from the TOA
    if(associated(totabcktoa).or.associated(totabcksfc)) then
        if (.not.associated(totextcoef) .or. .not.associated(totbckcoef)) then
             print*,trim(Iam),' : TOTEXTCOEF and TOTBCKCOEF and their children needs to be requested in HISTORY.rc.',&
                           ' Cannot produce TOTABCKTOA or TOTABCKSFC variables without these exports.'
             VERIFY_(100)
        endif

       ind532 = 0
       do w = 1, size(self%wavelengths_profile) ! find index for 532nm to compute TBA
          if ((self%wavelengths_profile(w)*1.e-9 .ge. 5.31e-7) .and. &
              (self%wavelengths_profile(w)*1.e-9 .le. 5.33e-7)) then
             ind532 = w
             exit
          end if
       end do

       if (ind532 == 0) then
          print*,trim(Iam),' : 532nm wavelengths is not present in GOCART2G_GridComp.rc.',&
                           ' Cannot produce TOTBCKCOEF variable without 532nm wavelength.'
          VERIFY_(100)
       end if

        ! Pressure at layer edges (ple shape (im,jm, km+1) on the edge

       i1 = lbound(ple, 1); i2 = ubound(ple, 1)
       j1 = lbound(ple, 2); j2 = ubound(ple, 2)
                            km = ubound(ple, 3) ! km =72 index starts at 0
       ! Pressure for each layer
       allocate(P(i1:i2,j1:j2,km), __STAT__)
       do k = 1, km
           P(:,:,k) = 0.5 * (ple(:,:,k-1) + ple(:,:,k))   ! in Pa
       enddo

      !molecular backscattering cross section for each layer at 532nm: Cair  * P(Pa) / T(K)
      !Cair = 4.51944e-9 at 532nm # unit K Pa-1 m-1 sr-1 http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19960051003.pdf
       allocate(backscat_mol(i1:i2,j1:j2,km), __STAT__)
       backscat_mol = (5.45e-32/1.380648e-23) * (532./550.)**(-4.0)  * P / T
       ! tau mol for each layer
       allocate(tau_mol_layer(i1:i2,j1:j2,km), delz(i1:i2,j1:j2,km),__STAT__)
       delz  = delp / (MAPL_GRAV * airdens)
       tau_mol_layer = backscat_mol * 8.* pi /3. * delz

       ! tau aer for each layer
       allocate(tau_aer_layer(i1:i2,j1:j2,km), __STAT__)
       tau_aer_layer = totextcoef(:,:,:,ind532) * delz

       allocate(tau_aer(i1:i2,j1:j2), __STAT__)
       allocate(tau_mol(i1:i2,j1:j2), __STAT__)

       ! TOTAL ABCK TOA
       ! top layer
       totabcktoa(:,:,1) = (totbckcoef(:,:,1,ind532) + backscat_mol(:,:,1)) * exp(-tau_aer_layer(:,:,1)) * exp(-tau_mol_layer(:,:,1))
       ! layer 2 to the layer at the surface(km)
       do k = 2, km
           tau_aer = 0.
           tau_mol = 0. ! for each layer
           do kk = 1, k
             tau_aer = tau_aer + tau_aer_layer(:,:,kk)
             tau_mol = tau_mol + tau_mol_layer(:,:,kk)
           enddo
           tau_aer = tau_aer + 0.5 *  tau_aer_layer(:,:,k)
           tau_mol = tau_mol + 0.5 *  tau_mol_layer(:,:,k)
           totabcktoa(:,:,k) = (totbckcoef(:,:,k,ind532) + backscat_mol(:,:,k)) * exp(-tau_aer) * exp(-tau_mol)
       enddo

       ! TOTAL ABCK SFC
       ! bottom layer
       totabcksfc(:,:,km) = (totbckcoef(:,:,km,ind532) + backscat_mol(:,:,km)) * exp(-tau_aer_layer(:,:,km)) * exp(-tau_mol_layer(:,:,km))
       ! layer 2nd from the surface to the top of the atmoshere (km)
       do k = km-1, 1, -1
           tau_aer = 0.
           tau_mol = 0. ! for each layer
           do kk = km, k+1, -1
             tau_aer = tau_aer + tau_aer_layer(:,:,kk)
             tau_mol = tau_mol + tau_mol_layer(:,:,kk)
           enddo
           tau_aer = tau_aer + 0.5 *  tau_aer_layer(:,:,k)
           tau_mol = tau_mol + 0.5 *  tau_mol_layer(:,:,k)
           totabcksfc(:,:,k) = (totbckcoef(:,:,k,ind532) + backscat_mol(:,:,k)) * exp(-tau_aer) * exp(-tau_mol)
       enddo

   endif ! end of total attenuated backscatter coef calculation

    RETURN_(ESMF_SUCCESS)

  end subroutine Run2


!===============================================================================

  subroutine getInstances_ (aerosol, myCF, species, rc)

!   Description: Fills the GOCART_State (aka, self%instance_XX) with user
!                defined instances from the GOCART2G_GridComp.rc.

    implicit none

    character (len=*),                intent(in   )  :: aerosol
    type (ESMF_Config),               intent(inout)  :: myCF
    type(Constituent),                intent(inout)  :: species
    integer,                          intent(  out)  :: rc


!   locals
    integer                                          :: i
    integer                                          :: n_active
    integer                                          :: n_passive
    integer                                          :: n_instances
    character (len=ESMF_MAXSTR)                      :: inst_name

    __Iam__('GOCART2G::getInstances_')

!--------------------------------------------------------------------------------------

!   Begin...
    n_active  = ESMF_ConfigGetLen (myCF, label='ACTIVE_INSTANCES_'//trim(aerosol)//':', __RC__)
    n_passive = ESMF_ConfigGetLen (myCF, label='PASSIVE_INSTANCES_'//trim(aerosol)//':', __RC__)
    n_instances = n_active + n_passive
    allocate (species%instances(n_instances), __STAT__)

!   !Fill the instances list with active instances first
    call ESMF_ConfigFindLabel (myCF, 'ACTIVE_INSTANCES_'//trim(aerosol)//':', __RC__)
    do i = 1, n_active
       call ESMF_ConfigGetAttribute (myCF, inst_name, __RC__)
       species%instances(i)%name = inst_name
       species%instances(i)%is_active = .true.
    end do
    species%n_active = n_active

!   !Now fill instances list with passive instances
    call ESMF_ConfigFindLabel (myCF, 'PASSIVE_INSTANCES_'//trim(aerosol)//':', __RC__)
    do i = n_active+1, n_active+n_passive
       call ESMF_ConfigGetAttribute (myCF, inst_name, __RC__)
       species%instances(i)%name = inst_name
       species%instances(i)%is_active = .false.
    end do


    RETURN_(ESMF_SUCCESS)

  end subroutine getInstances_


!====================================================================================
  subroutine createInstances_ (self, GC, rc)

!   Description: Creates GOCART2G children. Active instances must be created first. If
!     additional GOCART2G children are added, this subroutine will need to be updated.

    implicit none

    type (GOCART_State), pointer,            intent(in   )     :: self
    type (ESMF_GridComp),                    intent(inout)     :: GC
    integer,                                 intent(  out)     :: rc

    ! locals
    integer                                                    :: i

    __Iam__('GOCART2G::createInstances_')

!-----------------------------------------------------------------------------------
!   Begin...

!   Active instances must be created first! This ordering is necessary for
!   filing the AERO states that are passed to radiation.
!   This is achieved by arranging the names of the active instances first.

    call addChildren__ (gc, self%DU, setServices=DU2G_setServices, __RC__)
    call addChildren__ (gc, self%SS, setServices=SS2G_setServices, __RC__)
    call addChildren__ (gc, self%CA, setServices=CA2G_setServices, __RC__)
    call addChildren__ (gc, self%SU, setServices=SU2G_setServices, __RC__)
    call addChildren__ (gc, self%NI, setServices=NI2G_setServices, __RC__)

    RETURN_(ESMF_SUCCESS)

    contains

        subroutine addChildren__ (gc, species, setServices, rc)

          type (ESMF_GridComp),            intent(inout)     :: gc
          type(Constituent),               intent(inout)     :: species
          external                                           :: setServices
          integer,                         intent(  out)     :: rc

          ! local
          integer  :: n

          __Iam__('GOCART2G::createInstances_::addChildren__')

          n=size(species%instances)

          do i = 1, n
             species%instances(i)%id = MAPL_AddChild(gc, name=species%instances(i)%name, SS=SetServices, __RC__)
          end do

        RETURN_(ESMF_SUCCESS)

     end subroutine addChildren__

  end subroutine createInstances_

!===================================================================================
  subroutine serialize_bundle (state, rc)

    implicit none

!   !ARGUMENTS:
    type (ESMF_State)                             :: state
    integer,            intent(out)               :: rc

!   !Local
    character (len=ESMF_MAXSTR), allocatable      :: itemList(:)
    type (ESMF_State)                             :: child_state
    type (ESMF_StateItem_Flag), allocatable       :: itemTypes(:)
    type (ESMF_FieldBundle)                       :: bundle
    type (ESMF_Grid)                              :: grid
    type (ESMF_Field)                             :: field, serializedField

    character (len=ESMF_MAXSTR)                   :: binIndexstr
    character (len=ESMF_MAXSTR), allocatable      :: aeroName(:)

    real, pointer, dimension(:,:,:,:)             :: orig_ptr
    real, pointer, dimension(:,:,:)               :: ptr3d

    integer     :: b, i, j, n, rank, nbins

    __Iam__('GOCART2G::serialize_bundle')

!   !Description: Callback for AERO_RAD state used in GAAS module to provide a
!                 serialized ESMF_Bundle of aerosol fields.
!-----------------------------------------------------------------------------------
!   Begin...

!   Get list of child states within state and add to aeroList
!   Remember, AERO_RAD contains its children's AERO_RAD states
!   ----------------------------------------------------------
    call ESMF_StateGet (state, itemCount=n, __RC__)
    allocate (itemList(n), __STAT__)
    allocate (itemTypes(n), __STAT__)
    call ESMF_StateGet (state, itemNameList=itemList, itemTypeList=itemTypes, __RC__)

!  Create empty ESMF_FieldBundle to add Children's aerosol fields to
   bundle = ESMF_FieldBundleCreate(name="serialized_aerosolBundle", __RC__)
   call MAPL_StateAdd(state, bundle, __RC__)

   do i = 1, n
      if (itemTypes(i) /= ESMF_StateItem_State) cycle ! exclude non-states
      call ESMF_StateGet (state, trim(itemList(i)), child_state, __RC__)
      call ESMF_AttributeGet (child_state, name='internal_variable_name', itemCount=nbins, __RC__)
      allocate (aeroName(nbins), __STAT__)
      call ESMF_AttributeGet (child_state, name='internal_variable_name', valueList=aeroName, __RC__)


      do b = 1, size(aeroName)
         call ESMF_StateGet (child_state, trim(aeroName(b)), field, __RC__)
         call ESMF_FieldGet (field, rank=rank, __RC__)

         if (rank == 3) then
            call MAPL_FieldBundleAdd (bundle, field, __RC__)

         else if (rank == 4) then ! serialize 4d variables to mulitple 3d variables
            call ESMF_FieldGet (field, grid=grid, __RC__)
            call MAPL_GetPointer (child_state, orig_ptr, trim(aeroName(b)), __RC__)
            do j = 1, size(orig_ptr, 4)
               write (binIndexstr, '(I0.3)') j
               ptr3d => orig_ptr(:,:,:,j)
               serializedField = ESMF_FieldCreate (grid=grid, datacopyFlag=ESMF_DATACOPY_REFERENCE, &
                                              farrayPtr=ptr3d, name=trim(aeroName(b))//trim(binIndexstr), __RC__)
               call MAPL_FieldBundleAdd (bundle, serializedField, __RC__) ! probably need to add a flag to allow for adding multilple fields of the same name.
            end do ! do j
         end if ! if (rank
      end do ! do b
      deallocate (aeroName, __STAT__)
   end do ! do i

  end subroutine serialize_bundle

!===================================================================================
  subroutine run_aerosol_optics (state, rc)

    implicit none

!   !ARGUMENTS:
    type (ESMF_State)                                :: state
    integer,            intent(out)                  :: rc

!   !Local
    real, dimension(:,:,:), pointer                  :: ple
    real, dimension(:,:,:), pointer                  :: rh
    real, dimension(:,:,:), pointer                  :: var

    character (len=ESMF_MAXSTR)                      :: fld_name

    real(kind=8), dimension(:,:,:),pointer           :: ext_, ssa_, asy_      ! (lon:,lat:,lev:)
    real(kind=8), dimension(:,:,:), allocatable      :: ext,  ssa,  asy       ! (lon:,lat:,lev:)

    integer                                          :: i, n, b, j
    integer                                          :: i1, j1, i2, j2, km
    integer                                          :: band
    integer, parameter                               :: n_bands = 1

    character (len=ESMF_MAXSTR), allocatable         :: itemList(:), aeroList(:)
    type (ESMF_State)                                :: child_state
    real, pointer,     dimension(:,:,:)              :: as_ptr_3d

    type (ESMF_StateItem_Flag), allocatable          :: itemTypes(:)

    __Iam__('GOCART2G::run_aerosol_optics')

!   Description: Used in Radiation gridded components to provide aerosol properties
!-----------------------------------------------------------------------------------
!   Begin...

!   Radiation band
!   --------------
    call ESMF_AttributeGet(state, name='band_for_aerosol_optics', value=band, __RC__)

!   Relative humidity
!   -----------------
    call ESMF_AttributeGet(state, name='relative_humidity_for_aerosol_optics', value=fld_name, __RC__)
    call MAPL_GetPointer(state, RH, trim(fld_name), __RC__)

!   Pressure at layer edges
!   ------------------------
    call ESMF_AttributeGet(state, name='air_pressure_for_aerosol_optics', value=fld_name, __RC__)
    call MAPL_GetPointer(state, PLE, trim(fld_name), __RC__)

    i1 = lbound(ple, 1); i2 = ubound(ple, 1)
    j1 = lbound(ple, 2); j2 = ubound(ple, 2)
                         km = ubound(ple, 3)

    allocate(ext(i1:i2,j1:j2,km),  &
             ssa(i1:i2,j1:j2,km),  &
             asy(i1:i2,j1:j2,km), __STAT__)

!   Get list of child states within state and add to aeroList
!   ---------------------------------------------------------
    call ESMF_StateGet (state, itemCount=n, __RC__)
    allocate (itemList(n), __STAT__)
    allocate (itemTypes(n), __STAT__)
    call ESMF_StateGet (state, itemNameList=itemList, itemTypeList=itemTypes, __RC__)

    b=0
    do i = 1, n
        if (itemTypes(i) == ESMF_StateItem_State) then
            b = b + 1
        end if
    end do

    allocate (aeroList(b), __STAT__)

    j = 1
    do i = 1, n
        if (itemTypes(i) == ESMF_StateItem_State) then
            aeroList(j) = trim(itemList(i))
            j = j + 1
        end if
    end do

    ext = 0.0d0
    ssa = 0.0d0
    asy = 0.0d0

!  ! Get aerosol optic properties from children
   do i = 1, size(aeroList)
        call ESMF_StateGet(state, trim(aeroList(i)), child_state, __RC__)

!       ! set RH in child's aero state
        call ESMF_AttributeGet(child_state, name='relative_humidity_for_aerosol_optics', value=fld_name, __RC__)

        if (fld_name /= '') then
            call MAPL_GetPointer(child_state, as_ptr_3d, trim(fld_name), __RC__)
            as_ptr_3d = rh
        end if

!       ! set PLE in child's aero state
        call ESMF_AttributeGet(child_state, name='air_pressure_for_aerosol_optics', value=fld_name, __RC__)

        if (fld_name /= '') then
            call MAPL_GetPointer(child_state, as_ptr_3d, trim(fld_name), __RC__)
            as_ptr_3d = ple
        end if

!       ! set band in child's aero state
        call ESMF_AttributeSet(child_state, name='band_for_aerosol_optics', value=band, __RC__)

!       ! execute the aerosol optics method
        call ESMF_MethodExecute(child_state, label="aerosol_optics", __RC__)

!       ! Retrieve extinction from each child
        call ESMF_AttributeGet(child_state, name='extinction_in_air_due_to_ambient_aerosol', value=fld_name, __RC__)
        if (fld_name /= '') then
            call MAPL_GetPointer(child_state, ext_, trim(fld_name), __RC__)
        end if

!       ! Retrieve scattering extinction from each child
        call ESMF_AttributeGet(child_state, name='single_scattering_albedo_of_ambient_aerosol', value=fld_name, __RC__)
        if (fld_name /= '') then
            call MAPL_GetPointer(child_state, ssa_, trim(fld_name), __RC__)
        end if

!       ! Retrieve asymetry parameter multiplied by scatering extiction from each child
        call ESMF_AttributeGet(child_state, name='asymmetry_parameter_of_ambient_aerosol', value=fld_name, __RC__)
        if (fld_name /= '') then
            call MAPL_GetPointer(child_state, asy_, trim(fld_name), __RC__)
        end if

!       ! Sum aerosol optic properties from each child
        ext = ext + ext_
        ssa = ssa + ssa_
        asy = asy + asy_

    end do


!   ! Set ext, ssa, asy to equal the sum of ext, ssa, asy from the children. This is what is passed to radiation.
    call ESMF_AttributeGet(state, name='extinction_in_air_due_to_ambient_aerosol', value=fld_name, __RC__)
    if (fld_name /= '') then
        call MAPL_GetPointer(state, var, trim(fld_name), __RC__)
        var = ext(:,:,:)
    end if

    call ESMF_AttributeGet(state, name='single_scattering_albedo_of_ambient_aerosol', value=fld_name, __RC__)
    if (fld_name /= '') then
        call MAPL_GetPointer(state, var, trim(fld_name), __RC__)
        var = ssa(:,:,:)
    end if

    call ESMF_AttributeGet(state, name='asymmetry_parameter_of_ambient_aerosol', value=fld_name, __RC__)
    if (fld_name /= '') then
        call MAPL_GetPointer(state, var, trim(fld_name), __RC__)
        var = asy(:,:,:)
    end if

    deallocate(ext, ssa, asy, __STAT__)



   RETURN_(ESMF_SUCCESS)

  end subroutine run_aerosol_optics

!=====================================================================================================

  subroutine aerosol_activation_properties(state, rc)

    implicit none

!   Arguments
!   ---------
    type(ESMF_State)     :: state
    integer, intent(out) :: rc

!   Local
!   ---------
    character(len=ESMF_MAXSTR)      :: mode              ! mode name
    character(len=ESMF_MAXSTR)      :: mode_             ! lowercase mode name

    type(ESMF_State)                :: child_state

    real, dimension(:,:,:), pointer :: ple               ! pressure at the edges of model layers
    real, dimension(:,:,:), pointer :: temperature       ! air temperature
    real, dimension(:,:),   pointer :: f_land            ! fraction of land type in a grid cell

    real, dimension(:,:,:), pointer :: f                 ! correction factor for sea salt

    real, dimension(:,:,:), allocatable :: q             ! aerosol mass mixing ratio
    real, dimension(:,:,:,:), pointer   :: ptr_4d        ! aerosol mass mixing ratio (temporary)
    real, dimension(:,:,:), pointer     :: ptr_3d        ! aerosol mass mixing ratio (temporary)

    real, dimension(:,:,:), pointer :: num               ! number concentration of aerosol particles
    real, dimension(:,:,:), pointer :: diameter          ! dry size of aerosol
    real, dimension(:,:,:), pointer :: sigma             ! width of aerosol mode
    real, dimension(:,:,:), pointer :: density           ! density of aerosol
    real, dimension(:,:,:), pointer :: hygroscopicity    ! hygroscopicity of aerosol
    real, dimension(:,:,:), pointer :: f_dust            ! fraction of dust aerosol
    real, dimension(:,:,:), pointer :: f_soot            ! fraction of soot aerosol
    real, dimension(:,:,:), pointer :: f_organic         ! fraction of organic aerosol

    real                            :: max_clean          ! max mixing ratio before considered polluted
    real                            :: ccn_tuning         ! tunes conversion factors for sulfate
    character(LEN=ESMF_MAXSTR)      :: cld_micro

    character(len=ESMF_MAXSTR)      :: fld_name

    integer                         :: i2, j2, km
    integer                         :: b, i, j, n, aerosol_bin
    integer                         :: varNameLen

    character (len=ESMF_MAXSTR), allocatable  :: itemList(:), aeroList(:)
    type (ESMF_StateItem_Flag), allocatable   :: itemTypes(:)

!   auxilliary parameters
!   ---------------------
    real, parameter :: densSO4 = 1700.0
    real, parameter :: densORG = 1600.0
    real, parameter :: densSS  = 2200.0
    real, parameter :: densDU  = 1700.0
    real, parameter :: densBC  = 1600.0
    real, parameter :: densOC  =  900.0
    real, parameter :: densBR  =  900.0

    real, parameter :: k_SO4   = 0.65
    real, parameter :: k_ORG   = 0.20
    real, parameter :: k_SS    = 1.28
    real, parameter :: k_DU    = 0.0001
    real, parameter :: k_BC    = 0.0001
    real, parameter :: k_OC    = 0.0001
    real, parameter :: k_BR    = 0.0001

    integer, parameter :: UNKNOWN_AEROSOL_MODE = 2015

    __Iam__('GOCART2G::aerosol_activation_properties')

!   Begin...

!   Get list of child states within state and add to aeroList
!   ---------------------------------------------------------
    call ESMF_StateGet (state, itemCount=n, __RC__)
    allocate (itemList(n), __STAT__)
    allocate (itemTypes(n), __STAT__)
    call ESMF_StateGet (state, itemNameList=itemList, itemTypeList=itemTypes, __RC__)

    b=0
    do i = 1, n
       if ((itemTypes(i) == ESMF_StateItem_State) .and. (trim(itemList(i)(1:2)) /= 'NI')) then
          b = b + 1
       end if
    end do

    allocate (aeroList(b), __STAT__)

    j = 1
    do i = 1, n
       if ((itemTypes(i) == ESMF_StateItem_State) .and. (trim(itemList(i)(1:2)) /= 'NI')) then
          aeroList(j) = trim(itemList(i))
          j = j + 1
       end if
    end do

!   Aerosol mode
!   ------------
    call ESMF_AttributeGet(state, name='aerosol_mode', value=mode, __RC__)

!   Land fraction
!   -------------
    call ESMF_AttributeGet(state, name='fraction_of_land_type', value=fld_name, __RC__)
    call MAPL_GetPointer(state, f_land, trim(fld_name), __RC__)

!   Pressure at layer edges
!   ------------------------
    call ESMF_AttributeGet(state, name='air_pressure_for_aerosol_optics', value=fld_name, __RC__)
    call MAPL_GetPointer(state, ple, trim(fld_name), __RC__)

!   Temperature
!   -----------
    call ESMF_AttributeGet(state, name='air_temperature', value=fld_name, __RC__)
    call MAPL_GetPointer(state, temperature, trim(fld_name), __RC__)

    i2 = ubound(temperature, 1)
    j2 = ubound(temperature, 2)
    km = ubound(temperature, 3)

!   Activation activation properties
!   --------------------------------
    call ESMF_AttributeGet(state, name='aerosol_number_concentration', value=fld_name, __RC__)
    call MAPL_GetPointer(state, num, trim(fld_name), __RC__)

    call ESMF_AttributeGet(state, name='aerosol_dry_size', value=fld_name, __RC__)
    call MAPL_GetPointer(state, diameter, trim(fld_name), __RC__)

    call ESMF_AttributeGet(state, name='width_of_aerosol_mode', value=fld_name, __RC__)
    call MAPL_GetPointer(state, sigma, trim(fld_name), __RC__)

    call ESMF_AttributeGet(state, name='aerosol_density', value=fld_name, __RC__)
    call MAPL_GetPointer(state, density, trim(fld_name), __RC__)

    call ESMF_AttributeGet(state, name='aerosol_hygroscopicity', value=fld_name, __RC__)
    call MAPL_GetPointer(state, hygroscopicity, trim(fld_name), __RC__)

    call ESMF_AttributeGet(state, name='fraction_of_dust_aerosol', value=fld_name, __RC__)
    call MAPL_GetPointer(state, f_dust, trim(fld_name), __RC__)

    call ESMF_AttributeGet(state, name='fraction_of_soot_aerosol', value=fld_name, __RC__)
    call MAPL_GetPointer(state, f_soot, trim(fld_name), __RC__)

    call ESMF_AttributeGet(state, name='fraction_of_organic_aerosol', value=fld_name, __RC__)
    call MAPL_GetPointer(state, f_organic, trim(fld_name), __RC__)

!   Sea salt scaling fctor
!   ----------------------
    call ESMF_AttributeGet(state, name='max_q_clean', value=max_clean, __RC__)
    call ESMF_AttributeGet(state, name='cldmicro', value=cld_micro, __RC__)
    call ESMF_AttributeGet(state, name='ccn_tuning', value=ccn_tuning, __RC__)

!   Aerosol mass mixing ratios
!   --------------------------
    mode_ = trim(mode)
    mode_ = ESMF_UtilStringLowerCase(mode_, __RC__)

    allocate(q(i2,j2,km),  __STAT__)
    q = 0.0

    if (index(mode_, 'du00') > 0) then ! Dust
       ! dust is mapped one-to-one
       do i = 1, size(aeroList)
          if (index(aeroList(i), 'DU') > 0) then
             read (mode_(3:len(mode_)),*) aerosol_bin
             call ESMF_StateGet(state, trim(aeroList(i)), child_state, __RC__)
             call MAPL_GetPointer(child_state, ptr_4d, 'DU', __RC__)
             q = q + ptr_4d(:,:,:,aerosol_bin)
             ptr_3d => ptr_4d(:,:,:,aerosol_bin)

             hygroscopicity = k_DU
             density = densDU
          end if
       end do

    else if (index(mode_, 'ss00') > 0) then ! Sea Salt
       ! compute the total mass mixing ratio and impose a tri-modal size distribution
       do i = 1, size(aeroList)
          if (index(aeroList(i), 'SS') > 0) then
             call ESMF_StateGet(state, trim(aeroList(i)), child_state, __RC__)
             call MAPL_GetPointer(child_state, ptr_4d, 'SS', __RC__)
             do j = 1, ubound(ptr_4d, 4)
               q = q + ptr_4d(:,:,:,j)
               ptr_3d => ptr_4d(:,:,:,j)
             end do

             hygroscopicity = k_SS
             density = densSS
          end if
       end do

    else if (index(mode_, 'sulforg') > 0) then ! Sulfate
       hygroscopicity = 0.0
       density = 0.0

       do i = 1, size(aeroList)
          if (index(aeroList(i), 'SU') > 0) then
             call ESMF_StateGet(state, trim(aeroList(i)), child_state, __RC__)
             call MAPL_GetPointer(child_state, ptr_3d, 'SO4', __RC__)
             q = q + ptr_3d
             hygroscopicity = k_SO4 * ptr_3d + hygroscopicity
             density = densSO4 * ptr_3d + density
          end if

          if (index(aeroList(i), 'CA.oc') > 0) then
             call ESMF_StateGet(state, trim(aeroList(i)), child_state, __RC__)
             varNameLen = len_trim(aeroList(i))
!            the '5' refers to '_AERO', which we want to remove to get the CA component name (e.g. CA.oc, or CA.oc.data)
             varNameLen = varNameLen - 5
             call MAPL_GetPointer(child_state, ptr_3d, aeroList(i)(1:varNameLen)//'philic', __RC__)
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

    else if (index(mode_, 'bcphilic') > 0) then ! Black Carbon
       do i = 1, size(aeroList)
          if (index(aeroList(i), 'CA.bc') > 0) then
             call ESMF_StateGet(state, trim(aeroList(i)), child_state, __RC__)
             varNameLen = len_trim(aeroList(i))
!            the '5' refers to '_AERO', which we want to remove to get the CA component name (e.g. CA.bc, or CA.bc.data)
             varNameLen = varNameLen - 5
             call MAPL_GetPointer(child_state, ptr_3d, aeroList(i)(1:varNameLen)//'philic', __RC__)
             q = q + ptr_3d
             hygroscopicity = k_BC
             density = densBC
          end if
       end do

    else if (index(mode_, 'ocphilic') > 0) then ! Organic Carbon
       do i = 1, size(aeroList)
          if (index(aeroList(i), 'CA.oc') > 0) then
             call ESMF_StateGet(state, trim(aeroList(i)), child_state, __RC__)
             varNameLen = len_trim(aeroList(i))
!            the '5' refers to '_AERO', which we want to remove to get the CA component name (e.g. CA.oc, or CA.oc.data)
             varNameLen = varNameLen - 5
             call MAPL_GetPointer(child_state, ptr_3d, aeroList(i)(1:varNameLen)//'philic', __RC__)
             q = q + ptr_3d
             hygroscopicity = k_OC
             density = densOC
          end if
       end do

    else if (index(mode_, 'brcphilic') > 0) then ! Organic Carbon
       do i = 1, size(aeroList)
          if (index(aeroList(i), 'CA.br') > 0) then
             call ESMF_StateGet(state, trim(aeroList(i)), child_state, __RC__)
             varNameLen = len_trim(aeroList(i))
!            the '5' refers to '_AERO', which we want to remove to get the CA component name (e.g. CA.bc, or CA.bc.data)
             varNameLen = varNameLen - 5
             call MAPL_GetPointer(child_state, ptr_3d, aeroList(i)(1:varNameLen)//'philic', __RC__)
             q = q + ptr_3d
             hygroscopicity = k_BR
             density = densBR
          end if
       end do

    end if !(index(mode_, 'du00') > 0) then

!   Obtain aerosol activation properties of this aerosol mode
!   ---------------------------------------------------------
    call aap_(mode,               &
              q,                  &
              num,                &
              diameter,           &
              sigma,              &
              f_dust,             &
              f_soot,             &
              f_organic,          &
              density,            &
              ptr_3d,             &
              1, i2, 1, j2, km, &
              __RC__)

    deallocate(q, __STAT__)

    RETURN_(ESMF_SUCCESS)

   contains

    subroutine aap_(mode, q, num, diameter, sigma, f_dust, f_soot, f_organic, dens_, q_, &
                    i1, i2, j1, j2, km, rc)

     implicit none

     integer, intent(in) :: i1, i2                                  ! dimension bounds
     integer, intent(in) :: j1, j2                                  ! ... // ..
     integer, intent(in) :: km                                      ! ... // ..

     character(len=*),  intent(in )               :: mode           ! name of aerosol mode
     real, intent(in),  dimension(i1:i2,j1:j2,km) :: q              ! aerosol mass mixing ratio, kg kg-1
     real, intent(in),  dimension(i1:i2,j1:j2,km) :: q_             ! auxiliary mass
     real, intent(in),  dimension(i1:i2,j1:j2,km) :: dens_          ! density

     real, intent(out), dimension(i1:i2,j1:j2,km) :: num            ! number concentration of aerosol particles
     real, intent(out), dimension(i1:i2,j1:j2,km) :: diameter       ! dry size of aerosol
     real, intent(out), dimension(i1:i2,j1:j2,km) :: sigma          ! width of aerosol mode
     real, intent(out), dimension(i1:i2,j1:j2,km) :: f_dust         ! fraction of dust aerosol
     real, intent(out), dimension(i1:i2,j1:j2,km) :: f_soot         ! fraction of soot aerosol
     real, intent(out), dimension(i1:i2,j1:j2,km) :: f_organic      ! fraction of organic aerosol

     integer, intent(out) :: rc                                     ! return code

     ! local
     integer :: STATUS
     character(len=ESMF_MAXSTR) :: mode_
     character(len=ESMF_MAXSTR) :: Iam = 'GOCART::aerosol_activation_properties::aap_()'

     integer, parameter :: UNKNOWN_AEROSOL_MODE = 2015

     integer            :: kinx
     real               :: fmassaux, fmassclean
     real, dimension(3) :: TPI, DPGI, SIGI
     real, dimension(3) :: TPIclean, DPGIclean, SIGIclean
     real, dimension(i1:i2,j1:j2,km) :: qaux
      !real, parameter    :: max_clean = 5.0e-7  !max mixing ratio before considered polluted

     mode_ = trim(mode)
     mode_ = ESMF_UtilStringLowerCase(mode_, __RC__)

     num       = 0.0
     diameter  = 1.0e-9
     sigma     = log(2.0)
     f_dust    = 0.0
     f_soot    = 0.0
     f_organic = 0.0

     qaux=q !this corrects a bug

     if (index(mode_, 'ss00') > 0) then
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

     if (index(mode_, 'sulforg0') > 0) then
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

     case ('du001')
         sigma    = log(1.8)
         f_dust   = 1.0
         diameter = 1.46e-6
         num      = q / ((MAPL_PI/6.0) * densDU * diameter*diameter*diameter * exp(4.5*sigma*sigma))

     case ('du002')
         sigma    = log(1.8)
         f_dust   = 1.0
         diameter = 2.80e-6
         num      = q / ((MAPL_PI/6.0) * densDU * diameter*diameter*diameter * exp(4.5*sigma*sigma))

     case ('du003')
         sigma    = log(1.8)
         f_dust   = 1.0
         diameter = 4.80e-6
         num      = q / ((MAPL_PI/6.0) * densDU * diameter*diameter*diameter * exp(4.5*sigma*sigma))

     case ('du004')
         sigma    = log(1.8)
         f_dust   = 1.0
         diameter = 9.0e-6
         num      = q / ((MAPL_PI/6.0) * densDU * diameter*diameter*diameter * exp(4.5*sigma*sigma))

     case ('du005')
         sigma    = log(1.8)
         f_dust   = 1.0
         diameter = 16.0e-6
         num      = q / ((MAPL_PI/6.0) * densDU * diameter*diameter*diameter * exp(4.5*sigma*sigma))

     case ('ss001')
         sigma    = SIGI(1)
         diameter = DPGI(1)
         num      = TPI(1) * q / fmassaux

     case ('ss002')
         sigma    = SIGI(2)
         diameter = DPGI(2)
         num      = TPI(2) * q / fmassaux

     case ('ss003')
         sigma    = SIGI(3)
         diameter = DPGI(3)
         num      = TPI(3) * q / fmassaux

     case ('sulforg01')  !different distributions for clean and polluted environments
         where (q > max_clean)
             sigma    = SIGI(1)
             diameter = DPGI(1)
             num      = TPI(1) * qaux*ccn_tuning / (dens_*fmassaux)             ! only sulfate  mass
         elsewhere
             sigma    = SIGIclean(1)
             diameter = DPGIclean(1)
             num      = TPIclean(1) * qaux*ccn_tuning / (dens_*fmassclean)      ! only sulfate
         end where

     case ('sulforg02')
         where (q > max_clean)
             sigma    = SIGI(2)
             diameter = DPGI(2)
             num      = TPI(2) * qaux*ccn_tuning / (dens_*fmassaux)            ! only sulfate mass
         elsewhere
             sigma    = SIGIclean(2)
             diameter = DPGIclean(2)
             num      = TPIclean(2) * qaux*ccn_tuning / (dens_*fmassclean)     ! only sulfate
         end where

     case ('sulforg03')
         where (q > max_clean)
             sigma    = SIGI(3)
             diameter = DPGI(3)
             num      = TPI(3) * qaux*ccn_tuning / (dens_*fmassaux)           ! only sulfate mass
         elsewhere
             sigma    = SIGIclean(3)
             diameter = DPGIclean(3)
             num      = TPIclean(3) * qaux*ccn_tuning / (dens_*fmassclean)    ! only sulfate
         end where

     case ('bcphilic')
         sigma    = log(2.0)
         f_soot   = 1.0
         diameter = 0.0118*2e-6
         num = q / ((MAPL_PI/6.0) * densBC * diameter*diameter*diameter * exp(4.5*sigma*sigma))

     case ('ocphilic')
         sigma     = log(2.2)
         f_organic = 1.0
         diameter  = 0.0212*2.0e-6
         num = q / ((MAPL_PI/6.0) * densOrg * diameter*diameter*diameter * exp(4.5*sigma*sigma))

     case ('brcphilic')
         sigma     = log(2.2)
         f_organic = 1.0
         diameter  = 0.0212*2.0e-6
         num = q / ((MAPL_PI/6.0) * densOrg * diameter*diameter*diameter * exp(4.5*sigma*sigma))

     case default
         __raise__(UNKNOWN_AEROSOL_MODE,"Unknown aerosol mode used in the GOCART aerosol activation properties method: "//trim(mode))

     end select


     RETURN_(ESMF_SUCCESS)

    end subroutine aap_

  end subroutine aerosol_activation_properties


!===================================================================================
  subroutine get_monochromatic_aop (state, rc)

    implicit none

!   !ARGUMENTS:
    type (ESMF_State)                                :: state
    integer,            intent(out)                  :: rc

!   !Local
    real, dimension(:,:,:), pointer                  :: ple
    real, dimension(:,:,:), pointer                  :: rh
    real, dimension(:,:), pointer                    :: var

    character (len=ESMF_MAXSTR)                      :: fld_name

    real, dimension(:,:),pointer                     :: tau_      ! (lon:,lat:,lev:)
    real, dimension(:,:), allocatable                :: tau       ! (lon:,lat:,lev:)

    integer                                          :: i, n, b, j
    integer                                          :: i1, j1, i2, j2, km
    real                                             :: wavelength

    character (len=ESMF_MAXSTR), allocatable         :: itemList(:), aeroList(:)
    type (ESMF_State)                                :: child_state
    real, pointer,     dimension(:,:,:)              :: as_ptr_3d

    type (ESMF_StateItem_Flag), allocatable          :: itemTypes(:)

    __Iam__('GOCART2G::get_monochromatic_aop')

!   Description: Used in GAAS gridded component to provide aerosol properties
!-----------------------------------------------------------------------------------
!   Begin...

!   Radiation band
!   --------------
    call ESMF_AttributeGet(state, name='wavelength_for_aerosol_optics', value=wavelength, __RC__)

!   Relative humidity
!   -----------------
    call ESMF_AttributeGet(state, name='relative_humidity_for_aerosol_optics', value=fld_name, __RC__)
    call MAPL_GetPointer(state, RH, trim(fld_name), __RC__)

!   Pressure at layer edges
!   ------------------------
    call ESMF_AttributeGet(state, name='air_pressure_for_aerosol_optics', value=fld_name, __RC__)
    call MAPL_GetPointer(state, PLE, trim(fld_name), __RC__)

    i1 = lbound(ple, 1); i2 = ubound(ple, 1)
    j1 = lbound(ple, 2); j2 = ubound(ple, 2)
                         km = ubound(ple, 3)

    allocate(tau(i1:i2,j1:j2), __STAT__)
    tau = 0.0

!   Get list of child states within state and add to aeroList
!   ---------------------------------------------------------
    call ESMF_StateGet (state, itemCount=n, __RC__)
    allocate (itemList(n), __STAT__)
    allocate (itemTypes(n), __STAT__)
    call ESMF_StateGet (state, itemNameList=itemList, itemTypeList=itemTypes, __RC__)

    b=0
    do i = 1, n
        if (itemTypes(i) == ESMF_StateItem_State) then
            b = b + 1
        end if
    end do

    allocate (aeroList(b), __STAT__)

    j = 1
    do i = 1, n
        if (itemTypes(i) == ESMF_StateItem_State) then
            aeroList(j) = trim(itemList(i))
            j = j + 1
        end if
    end do

!   ! Get aerosol optic properties from children
    do i = 1, size(aeroList)
       call ESMF_StateGet(state, trim(aeroList(i)), child_state, __RC__)

!      ! set RH in child's aero state
       call ESMF_AttributeGet(child_state, name='relative_humidity_for_aerosol_optics', value=fld_name, __RC__)

       if (fld_name /= '') then
          call MAPL_GetPointer(child_state, as_ptr_3d, trim(fld_name), __RC__)
          as_ptr_3d = rh
       end if

!      ! set PLE in child's aero state
       call ESMF_AttributeGet(child_state, name='air_pressure_for_aerosol_optics', value=fld_name, __RC__)

       if (fld_name /= '') then
          call MAPL_GetPointer(child_state, as_ptr_3d, trim(fld_name), __RC__)
          as_ptr_3d = ple
       end if

!      ! set wavelength in child's aero state
       call ESMF_AttributeSet(child_state, name='wavelength_for_aerosol_optics', value=wavelength, __RC__)

!      ! execute the aerosol optics method
       call ESMF_MethodExecute(child_state, label="monochromatic_aerosol_optics", __RC__)

!      ! Retrieve extinction from each child
       call ESMF_AttributeGet(child_state, name='monochromatic_extinction_in_air_due_to_ambient_aerosol', value=fld_name, __RC__)
       if (fld_name /= '') then
          call MAPL_GetPointer(child_state, tau_, trim(fld_name), __RC__)
       end if

!      ! Sum aerosol optic properties from each child
       tau = tau + tau_
    end do

!   ! Set ext, ssa, asy to equal the sum of ext, ssa, asy from the children. This is what is passed to radiation.
    call ESMF_AttributeGet(state, name='monochromatic_extinction_in_air_due_to_ambient_aerosol', value=fld_name, __RC__)
    if (fld_name /= '') then
       call MAPL_GetPointer(state, var, trim(fld_name), __RC__)
       var = tau
    end if

    deallocate(tau, __STAT__)

   RETURN_(ESMF_SUCCESS)

  end subroutine get_monochromatic_aop


!===================================================================================
  subroutine get_mixRatioSum (state, rc)

    implicit none

!   !ARGUMENTS:
    type (ESMF_State)                                :: state
    integer,            intent(out)                  :: rc

!   !Local
    character (len=ESMF_MAXSTR), allocatable         :: itemList(:), aeroList(:)
    character (len=ESMF_MAXSTR)                      :: aeroName
    character (len=ESMF_MAXSTR)                      :: fld_name

    real, pointer, dimension(:,:,:)                  :: var
    real, dimension(:,:,:), allocatable              :: aeroOut
    type (ESMF_StateItem_Flag), allocatable          :: itemTypes(:)

    integer  :: b, i, n, j, im, jm, km

    __Iam__('GOCART2G::get_mixRatioSum')

!   Description: Used in GAAS gridded component to provide sum of aerosol mixing ratio
!--------------------------------------------------------------------------------------
!   Begin...

    call ESMF_AttributeGet(state, name='aerosolName', value=aeroName, __RC__)
    call ESMF_AttributeGet(state, name='im', value=im, __RC__)
    call ESMF_AttributeGet(state, name='jm', value=jm, __RC__)
    call ESMF_AttributeGet(state, name='km', value=km, __RC__)

    allocate(aeroOut(im,jm,km), __STAT__)
    aeroOut = 0.0

!   Get list of child states within state and add to aeroList
!   ---------------------------------------------------------
    call ESMF_StateGet (state, itemCount=n, __RC__)
    allocate (itemList(n), __STAT__)
    allocate (itemTypes(n), __STAT__)
    call ESMF_StateGet (state, itemNameList=itemList, itemTypeList=itemTypes, __RC__)

    b=0
    do i = 1, n
        if (itemTypes(i) == ESMF_StateItem_State) then
            b = b + 1
        end if
    end do

    allocate (aeroList(b), __STAT__)

    j = 1
    do i = 1, n
        if (itemTypes(i) == ESMF_StateItem_State) then
            aeroList(j) = trim(itemList(i))
            j = j + 1
        end if
    end do


!   Retrieve summed aerosol mixing ratios from active instances
    select case (trim(aeroName))
       case ('dust')
          call getAerosolSum ('DU', state, aeroList, aeroOut, __RC__)

          call ESMF_AttributeGet (state, name='sum_of_internalState_aerosol_DU', value=fld_name, __RC__)
          if (fld_name /= '') then
             call MAPL_GetPointer (state, var, trim(fld_name), __RC__)
             var = aeroOut
          end if

       case ('seasalt')
          call getAerosolSum ('SS', state, aeroList, aeroOut, __RC__)

          call ESMF_AttributeGet (state, name='sum_of_internalState_aerosol_SS', value=fld_name, __RC__)
          if (fld_name /= '') then
             call MAPL_GetPointer (state, var, trim(fld_name), __RC__)
             var = aeroOut
          end if

       case ('organicCarbon')
          call getAerosolSum ('CA.oc', state, aeroList, aeroOut, __RC__)

          call ESMF_AttributeGet (state, name='sum_of_internalState_aerosol_CA.oc', value=fld_name, __RC__)
          if (fld_name /= '') then
             call MAPL_GetPointer (state, var, trim(fld_name), __RC__)
             var = aeroOut
          end if

       case ('blackCarbon')
          call getAerosolSum ('CA.bc', state, aeroList, aeroOut, __RC__)

          call ESMF_AttributeGet (state, name='sum_of_internalState_aerosol_CA.bc', value=fld_name, __RC__)
          if (fld_name /= '') then
             call MAPL_GetPointer (state, var, trim(fld_name), __RC__)
             var = aeroOut
          end if

       case ('brownCarbon')
          call getAerosolSum ('CA.br', state, aeroList, aeroOut, __RC__)

          call ESMF_AttributeGet (state, name='sum_of_internalState_aerosol_CA.br', value=fld_name, __RC__)
          if (fld_name /= '') then
             call MAPL_GetPointer (state, var, trim(fld_name), __RC__)
             var = aeroOut
          end if

       case ('sulfate')
          call getAerosolSum ('SU', state, aeroList, aeroOut, __RC__)

          call ESMF_AttributeGet (state, name='sum_of_internalState_aerosol_SU', value=fld_name, __RC__)
          if (fld_name /= '') then
             call MAPL_GetPointer (state, var, trim(fld_name), __RC__)
             var = aeroOut
          end if

       case ('nitrate')
          call getAerosolSum ('NI', state, aeroList, aeroOut, __RC__)

          call ESMF_AttributeGet (state, name='sum_of_internalState_aerosol_NI', value=fld_name, __RC__)
          if (fld_name /= '') then
             call MAPL_GetPointer (state, var, trim(fld_name), __RC__)
             var = aeroOut
          end if

       case default
          !$omp critical (G2G_2)
          print *,"Invalid aerosolName of '",trim(aeroName), "' in GOCART2G::get_mixRatioSum"
          !$omp end critical (G2G_2)
    end select

contains
    subroutine getAerosolSum (aeroToken, state, aeroList, aeroOut, rc)

    implicit none

!   !ARGUMENTS:
    character (len=*), intent(in)                  :: aeroToken
    type (ESMF_State),           intent(in)        :: state
    character (len=ESMF_MAXSTR), intent(in)        :: aeroList(:)
    real, dimension(:,:,:),      intent(out)       :: aeroOut
    integer, optional,           intent(out)       :: rc

!   !LOCALS:
    integer                               :: i, endInd
    character (len=ESMF_MAXSTR)           :: fld_name
    type (ESMF_State)                     :: child_state
    real, pointer, dimension(:,:,:)       :: ptr3d


!   Begin...

    endInd = len_trim(aeroToken)

    aeroOut = 0.0
    do i = 1, size(aeroList)
       if (trim(aeroList(i)(1:endInd)) == trim(aeroToken)) then
          call ESMF_StateGet(state, trim(aeroList(i)), child_state, __RC__)
          call ESMF_MethodExecute(child_state, label="get_mixR", __RC__)
          call ESMF_AttributeGet(child_state, name='sum_of_internalState_aerosol', &
                                 value=fld_name, __RC__)
          if (fld_name /= '') then
             call MAPL_GetPointer(child_state, ptr3d, trim(fld_name), __RC__)
             aeroOut = aeroOut + ptr3d
          end if
       end if
    end do

    end subroutine getAerosolSum

  end subroutine get_mixRatioSum


end module GOCART2G_GridCompMod
