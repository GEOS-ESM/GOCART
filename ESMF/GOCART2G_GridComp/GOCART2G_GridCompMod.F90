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
    type (ESMF_Config)                            :: myCF
    type (GOCART_State), pointer                  :: self
    type (wrap_)                                  :: wrap

    __Iam__('SetServices')

!****************************************************************************
! Begin...

!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet (GC, NAME=COMP_NAME, __RC__)
    Iam = trim(COMP_NAME)//'::'//'SetServices'

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

!   Get instances to determine what children will be born
!   -----------------------------------------------------
    myCF = ESMF_ConfigCreate (__RC__)
    call ESMF_ConfigLoadFile (myCF, 'GOCART2G_GridComp.rc', __RC__)

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

!   Create children`s gridded components and invoke their SetServices
!   Active instances are created first
!   -----------------------------------------------------------------
    call createInstances_(self, GC, __RC__)

!   Define EXPORT states
!   This state is needed by radiation - It will contain 
!   aerosols and aerosol optics
!   --------------------------------------------------------
    call MAPL_AddExportSpec(GC,                       &
       short_name = 'AERO2G_RAD',                         &
       long_name  = 'aerosol_mass_mixing_ratios_ng',  &
       units      = 'kg kg-1',                        &
       dims       = MAPL_DimsHorzVert,                &
       vlocation  = MAPL_VLocationCenter,             &
       datatype   = MAPL_StateItem, __RC__)

!   This state is needed by MOIST - It will contain
!   aerosols
!   --------------------------------------------------------
    call MAPL_AddExportSpec(GC,                       &
       short_name = 'AERO2G_ACI',                     &
       long_name  = 'aerosol_cloud_interaction_ng',   &
       units      = 'kg kg-1',                        &
       dims       = MAPL_DimsHorzVert,                &
       vlocation  = MAPL_VLocationCenter,             &
       datatype   = MAPL_StateItem, __RC__)

!   This bundle is needed by surface for snow albedo modification
!   by aerosol settling and deposition
!   ~~~DEVELOPERS NOTE~~~ Change to StateItem when possible
!                         This will require refactoring Radiation
!   --------------------------------------------------------
    call MAPL_AddExportSpec(GC,                       &
       short_name = 'AERO2G_DP',                      &
       long_name  = 'aerosol_deposition_ng',          &
       units      = 'kg m-2 s-1',                     &
       dims       = MAPL_DimsHorzOnly,                &
       datatype   = MAPL_BundleItem, __RC__)

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
    character (len=ESMF_MAXSTR)              :: COMP_NAME

    type (MAPL_MetaComp),       pointer      :: MAPL
    type (ESMF_GridComp),       pointer      :: gcs(:)
    type (ESMF_State),          pointer      :: gex(:)
    type (ESMF_Grid)                         :: grid

    type (ESMF_State)                        :: aero, aero_aci
    type (ESMF_FieldBundle)                  :: aero_dp

    type (GOCART_State),      pointer        :: self
    type (wrap_)                             :: wrap

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

!   Call Generic Initialize
!   ----------------------------------------
    call MAPL_GenericInitialize (GC, import, export, clock, __RC__)

!   Get my internal state
!   ---------------------
    call ESMF_UserCompGetInternalState (GC, 'GOCART_State', wrap, STATUS)
    VERIFY_(STATUS)
    self => wrap%ptr

!   Get children and their export states from my generic state
!   -----------------------------------------------------------
    call MAPL_Get (MAPL, gcs=gcs, gex=gex, __RC__ )


!   Fill AERO2G_RAD, AERO2G_ACI, and AERO2G_DP with the children's states
!   ------------------------------------------------------------
    call ESMF_StateGet (export, 'AERO2G_RAD' , aero     , __RC__)
    call ESMF_StateGet (export, 'AERO2G_ACI' , aero_aci , __RC__)
    call ESMF_StateGet (export, 'AERO2G_DP'  , aero_dp  , __RC__)


!   Add children's AERO states to GOCART2G's AERO states
!   Only active instances are passed to radiation
!   ------------------------------------------------------
!   Active instances were created before passive instance. Summing the number of 
!   active instances will provide the last index for active instances within gcs
!   -----------------------------------------------------------------------   

!    tot = self%n_DU + self%n_SS + ...

    call add_aero_states_(self%DU%instances(:))
    call add_aero_states_(self%SS%instances(:))
    call add_aero_states_(self%SU%instances(:))
    call add_aero_states_(self%CA%instances(:))
    call add_aero_states_(self%NI%instances(:))


    ! Add variables to AERO_RAD state. Used in aerosol optics calculations
    call add_aero (aero, label='air_pressure_for_aerosol_optics',             label2='PLE', grid=grid, typekind=MAPL_R4, __RC__)
    call add_aero (aero, label='relative_humidity_for_aerosol_optics',        label2='RH',  grid=grid, typekind=MAPL_R4, __RC__)
    call add_aero (aero, label='extinction_in_air_due_to_ambient_aerosol',    label2='EXT', grid=grid, typekind=MAPL_R4, __RC__)
    call add_aero (aero, label='single_scattering_albedo_of_ambient_aerosol', label2='SSA', grid=grid, typekind=MAPL_R4, __RC__)
    call add_aero (aero, label='asymmetry_parameter_of_ambient_aerosol',      label2='ASY', grid=grid, typekind=MAPL_R4, __RC__)

    call ESMF_AttributeSet(aero, name='band_for_aerosol_optics',             value=0,     __RC__)

    ! Attach the aerosol optics method
    call ESMF_MethodAdd (aero, label='run_aerosol_optics', userRoutine=run_aerosol_optics, __RC__)

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
        __Iam__('Initialize::ad_aero_states_')
        
        do i = 1, size(instances)
           if (.not. instances(i)%is_active) cycle
           id = instances(i)%id
           
           call ESMF_GridCompGet (gcs(id), __RC__ )
           
           call ESMF_StateGet (gex(id), trim(instances(i)%name)//'_AERO', child_state, __RC__)
           call ESMF_StateAdd (aero, [child_state], __RC__)
           
           if (instances(i)%name(1:2) /= 'NI') then
              call ESMF_StateGet (gex(id), trim(instances(i)%name)//'_AERO_ACI', child_state, __RC__)
              call ESMF_StateAdd (aero_ACI, [child_state], __RC__)
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
    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS )
    VERIFY_(STATUS)

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

    character(len=ESMF_MAXSTR)          :: child_name
    integer                             :: i

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
    call MAPL_Get ( MAPL, gcs=gcs, gim=gim, gex=gex, INTERNAL_ESMF_STATE=internal, __RC__ )

!   Run the children
!   -----------------
    do i = 1, size(gcs)
      call ESMF_GridCompGet (gcs(i), NAME=child_name, __RC__ )
      if ((index(child_name, 'data')) == 0) then ! only execute Run2 method if a computational instance
         call ESMF_GridCompRun (gcs(i), importState=gim(i), exportState=gex(i), phase=2, clock=clock, __RC__)
      end if
    end do

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


end module GOCART2G_GridCompMod






