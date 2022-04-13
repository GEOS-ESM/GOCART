
#include "MAPL_Generic.h"

!-------------------------------------------------------------------------
!      NASA/GSFC, Global Modeling & Assimilation Office, Code 610.1      !
!-------------------------------------------------------------------------
!BOP
!

! !MODULE: Chem_AeroGeneric - Utilitarian subroutines used by GOCART2G children. 
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
   public setZeroKlid
   public setZeroKlid4d
   public findKlid
   public get_mixR
!
! !DESCRIPTION:
!
!  These subroutines perform repetitive tasks needed by GOCART2G children.
!
! !REVISION HISTORY:
!
!  March2020 Sherman, da Silva, Darmenov, Clune - created
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
       if (trim(field_name) == 'PLE') then
          call MAPL_FieldAllocCommit (field, dims=MAPL_DimsHorzVert, location=MAPL_VLocationEdge, typekind=typekind, hw=0, __RC__)
       else if ((trim(field_name) == 'FRLAND') .or. (trim(field_name) == 'monochromatic_EXT')) then
          call MAPL_FieldAllocCommit(field, dims=MAPL_DimsHorzOnly, location=MAPL_VLocationCenter, typekind=MAPL_R4, hw=0, __RC__)
       else
          call MAPL_FieldAllocCommit (field, dims=MAPL_DimsHorzVert, location=MAPL_VLocationCenter, typekind=typekind, hw=0, __RC__)
       end if
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
    type (ESMF_State),           intent(inout)   :: providerState
    type (ESMF_FieldBundle),     intent(inout)   :: bundle
    integer,                     intent(  out)   :: rc  ! return code

!   !Local
    type (ESMF_Field)     :: field, field2D
    type (ESMF_Grid)      :: grid
    integer               :: dimCount, i
    real, pointer         :: orig_ptr(:,:,:)
    real, pointer         :: ptr2d(:,:)
    character(len=ESMF_MAXSTR)  :: bin_index, varNameNew

!   Description: Adds deposition variables to deposition bundle

     __Iam__('append_to_bundle')

!   Dry deposition
!   ---------------
    call ESMF_StateGet (providerState, trim(prefix)//trim(varName), field, __RC__)
    call MAPL_AllocateCoupling (field, __RC__)
    call ESMF_FieldGet (field, dimCount=dimCount, __RC__)

    if (dimCount == 2) then ! this handles data instances
       call MAPL_FieldBundleAdd (bundle, field, __RC__) 

    else if (dimCount == 3) then ! this handles computational instances
       call ESMF_FieldGet (field, grid=grid, __RC__)
       call MAPL_GetPointer (providerState, orig_ptr, trim(prefix)//trim(varName), __RC__)

       if ((index(trim(varname), 'DU') > 0) .or. (index(trim(varname), 'SS') > 0)) then
          do i = 1, size(orig_ptr, 3)
             write (bin_index,'(A, I0.3)') '', i
             ptr2d => orig_ptr(:,:,i)
             field2D = ESMF_FieldCreate(grid=grid, datacopyflag=ESMF_DATACOPY_REFERENCE, farrayPtr=ptr2d,&
                                        name=trim(varName)//trim(bin_index) , __RC__)
             call MAPL_AllocateCoupling (field2D, __RC__)
             call MAPL_FieldBundleAdd (bundle, field2D, __RC__)
          end do
       end if

       if (index(trim(varname), 'SU') > 0) then ! only use SO4, which is the 3rd index
          ptr2d => orig_ptr(:,:,3)
          field2D = ESMF_FieldCreate(grid=grid, datacopyflag=ESMF_DATACOPY_REFERENCE, farrayPtr=ptr2d,&
                                     name=trim(varName)//'003' , __RC__)
          call MAPL_AllocateCoupling (field2D, __RC__)
          call MAPL_FieldBundleAdd (bundle, field2D, __RC__)
       end if

       if (index(trim(varname), 'CA.oc') > 0) then
          do i = 1, size(orig_ptr, 3)
             write (bin_index,'(A, I0.3)') '', i
             ptr2d => orig_ptr(:,:,i)
             varNameNew = 'OC'//varName(6:7)
             field2D = ESMF_FieldCreate(grid=grid, datacopyflag=ESMF_DATACOPY_REFERENCE, farrayPtr=ptr2d,&
                                        name=trim(varNameNew)//trim(bin_index) , __RC__)
             call MAPL_AllocateCoupling (field2D, __RC__)
             call MAPL_FieldBundleAdd (bundle, field2D, __RC__)
          end do
       end if

       if (index(trim(varname), 'CA.bc') > 0) then
          do i = 1, size(orig_ptr, 3)
             write (bin_index,'(A, I0.3)') '', i
             ptr2d => orig_ptr(:,:,i)
             varNameNew = 'BC'//varName(6:7)
             field2D = ESMF_FieldCreate(grid=grid, datacopyflag=ESMF_DATACOPY_REFERENCE, farrayPtr=ptr2d,&
                                        name=trim(varNameNew)//trim(bin_index) , __RC__)
             call MAPL_AllocateCoupling (field2D, __RC__)
             call MAPL_FieldBundleAdd (bundle, field2D, __RC__)
          end do
       end if

    else if (dimCount > 3) then
       if(mapl_am_i_root()) print*,'Chem_AeroGenric::append_to_bundle does not currently support fields greater than 3 dimensions'
       VERIFY_(824)       
    end if

    RETURN_(ESMF_SUCCESS)

  end subroutine append_to_bundle

!===================================================================================
!BOP
! !IROUTINE: setZeroKlid
   subroutine setZeroKlid(km, klid, int_ptr)

! !USES:
   implicit NONE

! !INPUT PARAMETERS:
   integer, intent(in) :: km   ! total model levels
   integer, intent(in) :: klid ! index for pressure level

! !INOUTPUT PARAMETERS:
   real, dimension(:,:,:), intent(inout) :: int_ptr ! aerosol pointer

! !DESCRIPTION: Set values to 0 where above klid
!
! !REVISION HISTORY:
!
! 25Aug2020 E.Sherman - Written 
!
! !Local Variables
   integer :: k

!EOP
!----------------------------------------------------------------------------------
!  Begin...

    do k = 1, km
       if (k < klid) then
          int_ptr(:,:,k) = 0.0
       else if (k >= klid) then
          exit
       end if
    end do

   end subroutine setZeroKlid
!===================================================================================
!BOP
! !IROUTINE: setZeroKlid
   subroutine setZeroKlid4d (km, klid, int_ptr)

! !USES:
   implicit NONE

! !INPUT PARAMETERS:
   integer, intent(in) :: km   ! total model levels
   integer, intent(in) :: klid ! index for pressure level

! !INOUTPUT PARAMETERS:
   real, dimension(:,:,:,:), intent(inout) :: int_ptr ! aerosol pointer

! !DESCRIPTION: Set values to 0 where above klid
!
! !REVISION HISTORY:
!
! 25Aug2020 E.Sherman - Written 
!
! !Local Variables
   integer :: k, n

!EOP
!----------------------------------------------------------------------------------
!  Begin...

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


!===================================================================================
!BOP
! !IROUTINE: findKlid
   subroutine findKlid (klid, plid, ple, rc)

! !USES:
   implicit NONE

! !INPUT PARAMETERS:
   integer, intent(inout) :: klid ! index for pressure lid
   real, intent(in)       :: plid ! pressure lid [hPa]
   real, dimension(:,:,:), intent(in) :: ple  ! air pressure [Pa]

! !OUTPUT PARAMETERS:
   integer, intent(out) :: rc ! return code; 0 - all is good
!                                            1 - bad

! !DESCRIPTION: Finds corresponding vertical index for defined pressure lid
!
! !REVISION HISTORY:
!
! 25Aug2020 E.Sherman - Written 
!
! !Local Variables
   integer :: k, j, i
   real :: plid_, diff, refDiff
   real, allocatable, dimension(:) :: pres  ! pressure at each model level [Pa]

!EOP
!----------------------------------------------------------------------------------
!  Begin...
   klid = 1
   rc = 0

!  convert from hPa to Pa
   plid_ = plid*100.0

   allocate(pres(ubound(ple,3)))

!  find pressure at each model level
   do k = 1, ubound(ple,3)
      pres(k) = ple(1,1,k)
   end do

!  find smallest absolute difference between plid and average pressure at each model level
   refDiff = 150000.0
   do k = 1, ubound(ple,3)
      diff = abs(pres(k) - plid_)
      if (diff < refDiff) then
         klid = k
         refDiff = diff
      end if
   end do

!  Check to make sure that all pressures at (i,j) were the same
   do j = 1, ubound(ple,2)
      do i = 1, ubound(ple,1)
         if (pres(klid) /= ple(i,j,klid)) then
            rc = 1
            return
         end if
      end do
   end do

   end subroutine findKlid

!===================================================================================
!BOP
! !IROUTINE: get_mixR
  subroutine get_mixR (state, rc)

! !USES:
    implicit none

!   !ARGUMENTS:
    type (ESMF_State)                                :: state
    integer,            intent(out)                  :: rc

!   !LOCALS:
    real, dimension(:,:,:), pointer                  :: ptr3d
    real, dimension(:,:,:,:), pointer                :: ptr4d
    real, dimension(:,:,:), pointer                  :: aeroSum
    character (len=ESMF_MAXSTR)                      :: fld_name
    integer                                          :: aeroN, i
    character (len=ESMF_MAXSTR), allocatable         :: aerosolNames(:)
    integer                                          :: status

!   Begin...

    call ESMF_AttributeGet(state, name='internal_variable_name', itemCount=aeroN, __RC__)
    allocate (aerosolNames(aeroN), __STAT__)
    call ESMF_AttributeGet(state, name='internal_variable_name', valueList=aerosolNames, __RC__)

!   Zero out previous aerosol sum value so it doesn't keep growing.
    call ESMF_AttributeGet (state, name='sum_of_internalState_aerosol', value=fld_name, __RC__)
    if (fld_name /= '') then
       call MAPL_GetPointer (state, aeroSum, trim(fld_name), __RC__)
       aeroSum = 0.0
    end if

    do i = 1, size(aerosolNames)
       if ((aerosolNames(i) == 'DU') .or. (aerosolNames(i) == 'SS')) then
          call MAPL_GetPointer (state, ptr4d, trim(aerosolNames(i)), __RC__)
          aeroSum = sum(ptr4d, dim=4) !DU and SS only have 1 internal state variable so no need to =+
       else
          call MAPL_GetPointer (state, ptr3d, trim(aerosolNames(i)), __RC__)
          aeroSum = aeroSum + ptr3d
       end if
    end do

 end subroutine get_mixR


end module  Chem_AeroGeneric


