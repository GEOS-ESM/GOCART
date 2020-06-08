
#include "MAPL_Generic.h"

!-------------------------------------------------------------------------
!      NASA/GSFC, Global Modeling & Assimilation Office, Code 610.1      !
!-------------------------------------------------------------------------
!BOP
!

! !MODULE:  Aerosol_Callbacks --- Call back methods for aerosol optics.
!                             
!
! !INTERFACE:
!
module  Chem_AeroGeneric

! !USES:
   use ESMF
   use MAPL
!   USE Chem_MieMod2G

   implicit none
   private

!
! !PUBLIC MEMBER FUNCTIONS:
   public add_aero
   public append_to_bundle
   public determine_data_driven
!   public aerosol_optics1
!
! !DESCRIPTION:
!
!  These modules compute aerosol optical properties for GOCART2G.
!
! !REVISION HISTORY:
!
!  23Mar2005 Colarco - Initial code.
!  11Jul2005 da Silva   Standardization.
!  30Dec2019 Sherman, da Silva, Darmenov, Clune - 2nd Gen. Made ESMF compliant.
!                                                 No longer relies on Chem_Reg
!
!EOP
!-------------------------------------------------------------------------
contains


!====================================================================================
  subroutine add_aero (state, label, label2, grid, typekind, ptr, rc)

!   Description: Adds fields to aero state for aerosol optics calcualtions. 

    implicit none

    type (ESMF_State),                          intent(inout)     :: state
    character (len=*),                          intent(in   )     :: label
    character (len=*),                          intent(in   )     :: label2
    type (ESMF_Grid),                           intent(inout)     :: grid
    integer,                                    intent(in   )     :: typekind
    real, pointer, dimension(:,:,:), optional,  intent(in   )     :: ptr
    integer,                                    intent(  out)     :: rc

    ! locals
    type (ESMF_Field)                                             :: field
    character (len=ESMF_MAXSTR)                                   :: field_name

    __Iam__('add_aero')

!----------------------------------------------------------------------------------
!   Begin...

    call ESMF_AttributeSet (state, name=trim(label), value=trim(label2),  __RC__)

    call ESMF_AttributeGet (state, name=trim(label), value=field_name, __RC__)
    if (field_name /= '') then
        field = MAPL_FieldCreateEmpty(trim(field_name), grid, __RC__)

        call MAPL_FieldAllocCommit (field, dims=MAPL_DimsHorzVert, location=MAPL_VLocationEdge, typekind=typekind, hw=0, __RC__)
        call MAPL_StateAdd (state, field, __RC__)
    end if

!   if (field_name /= '') then
!       field = ptr
!       call MAPL_StateAdd (state, field, __RC__)
!   end if

    RETURN_(ESMF_SUCCESS)

  end subroutine add_aero

!=====================================================================================

  subroutine determine_data_driven(COMP_NAME, data_driven, RC)

    !ARGUMENTS:
    integer, optional,               intent(  out)   :: RC          ! Error code:
    character (len=ESMF_MAXSTR),     intent(in   )   :: COMP_NAME
    logical,                         intent(  out)   :: data_driven

    !Local
    integer                                          :: i

!   Description: Determines whether gridded component is data driven or not.

     __Iam__('determine_data_driven')

!   Begin... 

!   Is DU data driven?
!   ------------------
    data_driven = .false.

    i = index(COMP_NAME, 'data')
    if (i > 0) then
      data_driven = .true.
    end if

    RETURN_(ESMF_SUCCESS)

  end subroutine determine_data_driven

!=====================================================================================

  subroutine append_to_bundle(varName, providerState, prefix, bundle, rc)

    implicit none

!   !ARGUMENTS:
    character (len=*),           intent(in   )   :: varName, prefix
    type (ESMF_State),           intent(in   )   :: providerState
    type (ESMF_FieldBundle),     intent(inout)   :: bundle
    integer,                     intent(  out)   :: rc  ! return code

!   !Local
    type (ESMF_Field)                              :: field

!   Description: Adds deposition variables to deposition bundle

     __Iam__('append_to_bundle')

!   Dry deposition
!   ---------------
    call ESMF_StateGet (providerState, trim(prefix)//trim(varName), field, __RC__)
    call MAPL_AllocateCoupling (field, __RC__)
    call MAPL_FieldBundleAdd (bundle, field, __RC__)

    RETURN_(ESMF_SUCCESS)

  end subroutine append_to_bundle

!===================================================================================
#if 0
  subroutine aerosol_optics1(state, rc)

    implicit none

!   !ARGUMENTS:
    type (ESMF_State)                                :: state
    integer,            intent(out)                  :: rc

!   !Local
    integer, parameter                               :: DP=kind(1.0d0)
    real, dimension(:,:,:), pointer                  :: ple, rh
    real(kind=DP), dimension(:,:,:), pointer         :: var
    real, dimension(:,:,:,:), pointer                :: q, q_4d
    integer, allocatable                             :: opaque_self(:)
    type(C_PTR)                                      :: address
    type(SS2G_GridComp), pointer                     :: self

    character (len=ESMF_MAXSTR)                      :: fld_name, int_fld_name
    type(ESMF_Field)                                 :: fld
    character (len=ESMF_MAXSTR)                      :: COMP_NAME

    real(kind=DP), dimension(:,:,:), allocatable     :: ext_s, ssa_s, asy_s  ! (lon:,lat:,lev:)
    real, dimension(:,:,:), allocatable              :: x
    integer                                          :: instance
    integer                                          :: n, nbins, dims(4)
    integer                                          :: i1, j1, i2, j2, km
    integer                                          :: band, offset
    integer, parameter                               :: n_bands = 1

    integer :: i, j, k

    __Iam__('aerosol_optics1')

!   Begin... 
if(mapl_am_i_root()) print*,'aerosol_optics1 is working!'
!   Mie Table instance/index
!   ------------------------
    call ESMF_AttributeGet(state, name='mie_table_instance', value=instance, __RC__)

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
             asy_s(i1:i2, j1:j2, km), &
                 x(i1:i2, j1:j2, km), __STAT__)

    call ESMF_AttributeGet(state, name='internal_varaible_name', value=int_fld_name, __RC__)
    call ESMF_StateGet (state, trim(int_fld_name), field=fld, __RC__) !add as attribute - dont hard code?
    call ESMF_FieldGet (fld, farrayPtr=q, __RC__)

    nbins = size(q,4)

    allocate(q_4d(i1:i2, j1:j2, km, nbins), __STAT__)

    do n = 1, nbins
       do k = 1, km
          x(:,:,k) = ((PLE(:,:,k) - PLE(:,:,k-1))*0.01)*(100./MAPL_GRAV)
          q_4d(:,:,k,n) = x(:,:,k) * q(:,:,k,n)
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

    RETURN_(ESMF_SUCCESS)

  contains

!    subroutine mie_(mie_table, aerosol_names, nb, offset, q, rh, bext_s, bssa_s, basym_s, rc)
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

    __Iam__('SS2G::aerosol_optics::mie_')

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

  end subroutine aerosol_optics1

#endif

!===================================================================================


end module  Chem_AeroGeneric


