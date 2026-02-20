#include "MAPL.h"

!-------------------------------------------------------------------------
!      NASA/GSFC, Global Modeling & Assimilation Office, Code 610.1      !
!-------------------------------------------------------------------------
!BOP
! !MODULE: Chem_AeroGeneric - Utilitarian subroutines used by GOCART2G children.

! !INTERFACE:
module Chem_AeroGeneric

   !USES:
   use ESMF
   use mapl_ErrorHandling, only: MAPL_Verify, MAPL_Assert, MAPL_Return
   use mapl3g_State_API, only: MAPL_StateGetPointer
   use mapl3g_Field_API, only: MAPL_FieldGet, MAPL_FieldCreate
   use mapl3g_FieldBundle_API, only: MAPL_FieldBundleAdd
   use mapl3g_VerticalStaggerLoc, only: VerticalStaggerLoc, VERTICAL_STAGGER_EDGE, VERTICAL_STAGGER_CENTER
   use mapl3g_UngriddedDims, only: UngriddedDims
   ! USE Chem_MieMod2G

   implicit none
   private

   !PUBLIC MEMBER FUNCTIONS:
   public add_aero
   public append_to_bundle
   public determine_data_driven
   public setZeroKlid
   public setZeroKlid4d
   public findKlid
   public get_mixR

   !DESCRIPTION:
   ! These subroutines perform repetitive tasks needed by GOCART2G children.

   !REVISION HISTORY:
   !  March2020 Sherman, da Silva, Darmenov, Clune - created
   !EOP

contains

   subroutine add_aero(state, label, label2, geom, km, typekind, ptr, rc)

      ! Description: Adds fields to aero state for aerosol optics calcualtions.

      type(ESMF_State), intent(inout) :: state
      character(len=*), intent(in) :: label
      character(len=*), intent(in) :: label2
      type(ESMF_Geom), intent(in) :: geom
      integer, optional, intent(in) :: km
      type(ESMF_TypeKind_Flag), optional, intent(in) :: typekind
      real, pointer, dimension(:,:,:), optional, intent(in) :: ptr
      integer, intent(out) :: rc

      ! locals
      type(ESMF_Field) :: field
      type(ESMF_Info) :: info
      character(len=:), allocatable :: field_name
      type(ESMF_TypeKind_Flag) :: typekind_
      integer :: status

      typekind_ = ESMF_TYPEKIND_R4
      if (present(typekind)) typekind_ = typekind

      call ESMF_InfoGetFromHost(state, info, _RC)
      call ESMF_InfoSet(info, key=trim(label), value=trim(label2), _RC)

      field_name = trim(label2)
      if (field_name /= "") then
         if (trim(field_name) == "PLE") then
            _ASSERT(present(km), "missing km for a 3D field")
            field = MAPL_FieldCreate( &
                 geom, typekind_, &
                 name=field_name, &
                 num_levels=km+1, &
                 vert_staggerloc=VERTICAL_STAGGER_EDGE, _RC)
         else if ((trim(field_name) == "FRLAND") .or. (trim(field_name) == "monochromatic_EXT")) then
            field = MAPL_FieldCreate(geom, typekind_, name=field_name, _RC)
         else
            _ASSERT(present(km), "missing km for a 3D field")
            field = MAPL_FieldCreate( &
                 geom, typekind_, &
                 name=field_name, &
                 num_levels=km, &
                 vert_staggerloc=VERTICAL_STAGGER_CENTER, _RC)
         end if
         call ESMF_StateAdd(state, [field], _RC)
      end if

      ! if (field_name /= "") then
      !     field = ptr
      !     call ESMF_StateAdd(state, [field], _RC)
      ! end if
      _UNUSED_DUMMY(ptr)

      _RETURN(_SUCCESS)

   end subroutine add_aero

   recursive subroutine determine_data_driven(comp_name, data_driven, rc)
      !ARGUMENTS:
      character(len=*), intent(in) :: comp_name
      logical, intent(out) :: data_driven
      integer, optional, intent(out) :: rc

      !Local
      integer :: i

      ! Description: Determines whether gridded component is data driven or not.

      data_driven = .false.
      i = index(comp_name, 'data')
      if (i > 0) then
         data_driven = .true.
      end if

      _RETURN(_SUCCESS)
   end subroutine determine_data_driven

   subroutine append_to_bundle(varname, provider_state, prefix, bundle, rc)
      !ARGUMENTS:
      character(len=*), intent(in) :: varname, prefix
      type(ESMF_State), intent(inout) :: provider_state
      type(ESMF_FieldBundle), intent(inout) :: bundle
      integer, intent(out) :: rc

      !Local
      type(ESMF_Field) :: field, field2d
      type(ESMF_Info) :: info
      type(ESMF_Geom), allocatable :: geom
      type(ESMF_TypeKind_Flag) :: typekind
      real, pointer :: orig_ptr(:,:,:)
      real, pointer :: ptr2d(:,:)
      character(len=ESMF_MAXSTR) :: bin_index
      character(:), allocatable :: units, stdname
      type(VerticalStaggerLoc) :: vert_stagger
      integer :: dim_count, iter, status

      ! Description: Adds deposition variables to deposition bundle

      ! Dry deposition
      call ESMF_StateGet(provider_state, trim(prefix)//trim(varname), field, _RC)
      call ESMF_FieldGet(field, dimCount=dim_count, _RC)

      _ASSERT(dim_count==2 .or. dim_count==3, "only 2d and 3d fields are supported")

      select case(dim_count)
      case(2) ! this handles data instances
         call MAPL_FieldBundleAdd(bundle, [field], _RC)

      case(3) ! this handles computational instances
         call MAPL_FieldGet(field, &
              geom=geom, &
              typekind=typekind, &
              units=units, &
              standard_name=stdname, &
              vert_staggerloc=vert_stagger, _RC)
         stdname = stdname(1:index(stdname, "(Bin")-1)
         call MAPL_StateGetPointer(provider_state, itemName=trim(prefix)//trim(varname), farrayPtr=orig_ptr, _RC)

         if ((index(trim(varname), "DU") > 0) .or. (index(trim(varname), "SS") > 0)) then
            do iter = 1, size(orig_ptr, 3)
               write (bin_index, "(A, I0.3)") "", iter
               ptr2d => orig_ptr(:,:,iter)
               field2d = MAPL_FieldCreate( &
                    geom, typekind, &
                    name=trim(varname)//trim(bin_index), &
                    ungridded_dims=UngriddedDims(), &
                    vert_staggerloc=vert_stagger, &
                    units=units, &
                    standard_name=stdname//' Bin '//trim(bin_index), &
                    long_name="unknown", _RC)
               call ESMF_FieldEmptyReset(field2d, status=ESMF_FIELDSTATUS_GEOMSET, _RC)
               call ESMF_FieldEmptyComplete( &
                    field2d, &
                    farray=ptr2d, &
                    indexflag=ESMF_INDEX_DELOCAL, &
                    datacopyflag=ESMF_DATACOPY_REFERENCE, _RC)
               call MAPL_FieldBundleAdd(bundle, [field2d], _RC)
            end do
         end if

         ! if (index(trim(varname), 'SU') > 0) then ! only use SO4, which is the 3rd index
         !    ptr2d => orig_ptr(:,:,3)
         !    field2d = ESMF_FieldCreate(grid=grid, datacopyflag=ESMF_DATACOPY_REFERENCE, farray=ptr2d,&
         !         name=trim(varname)//'003' , indexflag=ESMF_INDEX_DELOCAL, _RC)
         !    call ESMF_AttributeSet(field2d, name='DIMS', value=MAPL_DimsHorzOnly, _RC)
         !    call ESMF_AttributeSet(field2d, name='VLOCATION', value=MAPL_VLocationNone, _RC)
         !    call ESMF_AttributeSet(field2d, name='UNITS', value=units, _RC)
         !    call ESMF_AttributeSet(field2d, name='STANDARD_NAME', value=stdname//' Bin 003', _RC)
         !    call MAPL2_AllocateCoupling(field2d, _RC)
         !    call MAPL2_FieldBundleAdd(bundle, field2d, _RC)
         ! end if

         ! if (index(trim(varname), 'CA.oc') > 0) then
         !    do iter = 1, size(orig_ptr, 3)
         !       write (bin_index,'(A, I0.3)') '', iter
         !       ptr2d => orig_ptr(:,:,iter)
         !       varname_new = 'OC'//varname(6:7)
         !       field2d = ESMF_FieldCreate(grid=grid, datacopyflag=ESMF_DATACOPY_REFERENCE, farray=ptr2d,&
         !            name=trim(varname_new)//trim(bin_index) , indexflag=ESMF_INDEX_DELOCAL, _RC)
         !       call ESMF_AttributeSet(field2d, name='DIMS', value=MAPL_DimsHorzOnly, _RC)
         !       call ESMF_AttributeSet(field2d, name='VLOCATION', value=MAPL_VLocationNone, _RC)
         !       call ESMF_AttributeSet(field2d, name='UNITS', value=units, _RC)
         !       call ESMF_AttributeSet(field2d, name='STANDARD_NAME', value=stdname//' Bin '//trim(bin_index), _RC)
         !       call MAPL2_AllocateCoupling(field2d, _RC)
         !       call MAPL2_FieldBundleAdd(bundle, field2d, _RC)
         !    end do
         ! end if

         ! if (index(trim(varname), 'CA.bc') > 0) then
         !    do i = 1, size(orig_ptr, 3)
         !       write (bin_index,'(A, I0.3)') '', iter
         !       ptr2d => orig_ptr(:,:,iter)
         !       varname_new = 'BC'//varname(6:7)
         !       field2d = ESMF_FieldCreate(grid=grid, datacopyflag=ESMF_DATACOPY_REFERENCE, farray=ptr2d,&
         !            name=trim(varname_new)//trim(bin_index) , indexflag=ESMF_INDEX_DELOCAL, _RC)
         !       call ESMF_AttributeSet(field2d, name='DIMS', value=MAPL_DimsHorzOnly, _RC)
         !       call ESMF_AttributeSet(field2d, name='VLOCATION', value=MAPL_VLocationNone, _RC)
         !       call ESMF_AttributeSet(field2d, name='UNITS', value=units, _RC)
         !       call ESMF_AttributeSet(field2d, name='STANDARD_NAME', value=stdname//' Bin '//trim(bin_index), _RC)
         !       call MAPL2_AllocateCoupling(field2d, _RC)
         !       call MAPL2_FieldBundleAdd(bundle, field2d, _RC)
         !    end do
         ! end if

      case default
         _FAIL("dim_count is other than 2 and 3")
      end select

      _RETURN(_SUCCESS)
   end subroutine append_to_bundle

   !BOP
   !IROUTINE: setZeroKlid
   subroutine setZeroKlid(km, klid, int_ptr)

      !INPUT PARAMETERS:
      integer, intent(in) :: km   ! total model levels
      integer, intent(in) :: klid ! index for pressure level

      !INOUTPUT PARAMETERS:
      real, dimension(:,:,:), intent(inout) :: int_ptr ! aerosol pointer

      !DESCRIPTION: Set values to 0 where above klid
      !REVISION HISTORY:
      ! 25Aug2020 E.Sherman - Written

      !Local Variables
      integer :: k
      !EOP

      do k = 1, km
         if (k < klid) then
            int_ptr(:,:,k) = 0.0
         else if (k >= klid) then
            exit
         end if
      end do

   end subroutine setZeroKlid

   !BOP
   !IROUTINE: setZeroKlid
   subroutine setZeroKlid4d (km, klid, int_ptr)

      !INPUT PARAMETERS:
      integer, intent(in) :: km   ! total model levels
      integer, intent(in) :: klid ! index for pressure level

      !INOUTPUT PARAMETERS:
      real, dimension(:,:,:,:), intent(inout) :: int_ptr ! aerosol pointer

      !DESCRIPTION: Set values to 0 where above klid
      !REVISION HISTORY:
      ! 25Aug2020 E.Sherman - Written

      !Local Variables
      integer :: k, n
      !EOP

      do n = 1, ubound(int_ptr, 4)
         do k = 1, km
            if (k < klid) then
               int_ptr(:,:,k,n) = 0.0
            else if (k >= klid) then
               exit
            end if
         end do
      end do

   end subroutine setZeroKlid4d

   !BOP
   !IROUTINE: findKlid
   subroutine findKlid (klid, plid, ple, rc)

      !INPUT PARAMETERS:
      integer, intent(inout) :: klid ! index for pressure lid
      real, intent(in) :: plid ! pressure lid [hPa]
      real, dimension(:,:,:), intent(in) :: ple  ! air pressure [Pa]

      !OUTPUT PARAMETERS:
      integer, intent(out) :: rc ! return code; 0 - all is good, 1 - bad

      !DESCRIPTION: Finds corresponding vertical index for defined pressure lid
      !REVISION HISTORY:
      ! 25Aug2020 E.Sherman - Written

      !Local Variables
      integer :: k, j, i
      real :: plid_, diff, refDiff
      real, allocatable, dimension(:) :: pres  ! pressure at each model level [Pa]
      !EOP

      klid = 1
      rc = 0

      ! convert from hPa to Pa
      plid_ = plid*100.0

      allocate(pres(ubound(ple,3)))

      ! find pressure at each model level
      do k = 1, ubound(ple,3)
         pres(k) = ple(1,1,k)
      end do

      ! find smallest absolute difference between plid and average pressure at each model level
      refDiff = 150000.0
      do k = 1, ubound(ple,3)
         diff = abs(pres(k) - plid_)
         if (diff < refDiff) then
            klid = k
            refDiff = diff
         end if
      end do

      ! Check to make sure that all pressures at (i,j) were the same
      do j = 1, ubound(ple,2)
         do i = 1, ubound(ple,1)
            if (pres(klid) /= ple(i,j,klid)) then
               rc = 1
               return
            end if
         end do
      end do

   end subroutine findKlid

   !BOP
   !IROUTINE: get_mixR
   subroutine get_mixR (state, rc)

      !ARGUMENTS:
      type (ESMF_State) :: state
      integer, intent(out) :: rc

      !LOCALS:
      type(ESMF_Info) :: info
      real, dimension(:,:,:), pointer :: ptr3d
      real, dimension(:,:,:,:), pointer :: ptr4d
      real, dimension(:,:,:), pointer :: aeroSum
      character(len=:), allocatable :: fld_name
      integer :: i
      character(len=ESMF_MAXSTR), allocatable :: aerosolNames(:)
      integer :: status
      !EOP

      call ESMF_InfoGetFromHost(state, info, _RC)
      call ESMF_InfoGet(info, key="internal_variable_name", values=aerosolNames, _RC)

      ! Zero out previous aerosol sum value so it doesn't keep growing.
      call ESMF_AttributeGet (state, name="sum_of_internalState_aerosol", value=fld_name, _RC)
      if (fld_name /= '') then
         call MAPL_StateGetPointer(state, itemName=fld_name, farrayPtr=aeroSum, _RC)
         aeroSum = 0.0
      end if

      do i = 1, size(aerosolNames)
         if ((aerosolNames(i) == 'DU') .or. (aerosolNames(i) == 'SS')) then
            call MAPL_StateGetPointer(state, itemName=aerosolNames(i), farrayPtr=ptr4d, _RC)
            aeroSum = sum(ptr4d, dim=4) !DU and SS only have 1 internal state variable so no need to =+
         else
            call MAPL_StateGetPointer(state, itemName=aerosolNames(i), farrayPtr=ptr3d, _RC)
            aeroSum = aeroSum + ptr3d
         end if
      end do

   end subroutine get_mixR

end module Chem_AeroGeneric


