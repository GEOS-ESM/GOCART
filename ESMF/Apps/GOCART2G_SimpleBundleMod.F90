!
! This class extends MAPL_SimpleBundle with some convenience methods.
! This would better implemented as an extension of MPL_SimpleBundle.
!
! Arlindo da Silva <arlindo.dasilva@nasa.gov>, March 2022
! V. Buchard July 2022- added new subroutines
!----------------------------------------------------------------------------

#  include "MAPL_Generic.h"

module GOCART2G_SimpleBundleMod

   use ESMF
   use MAPL
!   use RegistryMod
!   use m_StrTemplate

   implicit NONE

   private

!  Inheritted from MAPL_SImpleBundle
!  ---------------------------------
   public MAPL_SimpleBundle
   public MAPL_SimpleBundleCreate
   public MAPL_SimpleBundlePrint
   public MAPL_SimpleBundleRead
   public MAPL_SimpleBundleWrite

!  Defined Here
!  ------------
   public GOCART2G_ESMFBundleCreate
   public GOCART2G_SimpleBundleCreate
   public GOCART2G_SimpleBundleRead
   public GOCART2G_SimpleBundleWrite
   public GOCART2G_ESMFBundleWrite

   integer, parameter :: READ_ONLY=1

CONTAINS

!..........................................................................................

  Function GOCART2G_ESMFBundleCreate (name, CF, rcname, Grid, rc)  result (Bundle)

    type(ESMF_FieldBundle)                    :: Bundle
    character(len=*),           intent(in)     :: name
    type(ESMF_Config),          intent(inout)  :: CF     ! resources
    character(len=*),           intent(in)     :: rcname ! name of resource with variable table
    type(ESMF_Grid),            intent(inout)  :: Grid
    integer, OPTIONAL,          intent(out)    :: rc


!                           ---

    character(len=255),allocatable    :: short_name(:)
    character(len=255),allocatable    :: long_name(:)
    character(len=255),allocatable    :: units(:)
    character(len=255), allocatable   :: table_info(:,:)
    type(ESMF_Field)                  :: Field
    integer :: i, j,im, jm, km, dims(7), nch
    integer :: linecount, columncount
                          __Iam__ ('SimpleBundleCreate')


!   Grid sizes
!   ----------
    call MAPL_GridGet(Grid, localCellCountPerDim = dims, __RC__)
    im = dims(1);  jm = dims(2);  km = dims(3)

!   Create an ESMF Bundle for holding variables
!   -------------------------------------------
    Bundle = ESMF_FieldBundleCreate ( name=name, __RC__ )
    call ESMF_FieldBundleSet ( bundle, grid=Grid, __RC__ )

!   Parse rc and fill in short, long and units
!   ------------------------------------------
    call ESMF_ConfigGetDim(CF, linecount, columncount, label='aop_variables::', __RC__)
    if (columncount /= 3) then
     __raise__(MAPL_RC_ERROR,"columncount should be 3 (short_name, units and long_name)")
    end if
    allocate(table_info(linecount, columncount), __STAT__)
    allocate(short_name(linecount), units(linecount), long_name(linecount), __STAT__)
    call ESMF_ConfigFindLabel(CF, 'aop_variables::', __RC__)
    do i = 1, linecount
       call ESMF_ConfigNextLine(CF, __RC__) 
       do j = 1, columncount
           call ESMF_ConfigGetAttribute(CF, table_info(i,j),__RC__)
       enddo
       short_name(i) = table_info(i,1)
       units(i)      = table_info(i,2)
       long_name(i)  = table_info(i,3)       
    enddo

!   add number of wavelengths nch
!------------------------
    nch  = ESMF_ConfigGetLen(CF,Label='wavelengths_in_nm:',__RC__)
!   Add fields to Bundle
!   --------------------
    do i = 1, size(short_name)
       Field=ESMF_FieldCreate(grid,name=short_name(i),typekind=ESMF_typekind_r4,ungriddedlbound=[1,1], ungriddedubound=[km,nch], __RC__)
       call ESMF_AttributeSet(field, name='ShortName', value=short_name(i),  __RC__)
       call ESMF_AttributeSet(field, name='LongName', value=long_name(i),  __RC__)
       call ESMF_AttributeSet(field, name='Units', value=units(i),  __RC__)

       call MAPL_FieldBundleAdd(Bundle, Field=Field, __RC__)
    end do
   
   
   end Function GOCART2G_ESMFBundleCreate

!..........................................................................................

  Function GOCART2G_SimpleBundleCreate (name, CF, rcname, Grid, rc, &
                               Levs, LevUnits, ptop, delp  )  result (self)

    type(MAPL_SimpleBundle)                    :: self
    character(len=*),           intent(in)     :: name
    type(ESMF_Config),          intent(inout)  :: CF     ! resources
    character(len=*),           intent(in)     :: rcname ! name of resource with variable table
    type(ESMF_Grid),            intent(inout)  :: Grid
    integer, OPTIONAL,          intent(out)    :: rc

                                                ! Vertical coordinates
    real(ESMF_KIND_R4), OPTIONAL,   intent(in) :: Levs(:)       ! Constant levels
    character(len=*),   OPTIONAL,   intent(in) :: LevUnits      ! Level units
                                                ! Lagrangian Control Volume Info
    real(ESMF_KIND_R4), OPTIONAL,   intent(in) ::   ptop        ! top pressure (Pa)
    real(ESMF_KIND_R4), OPTIONAL, pointer, &
                                    intent(in) ::   delp(:,:,:) ! layer thickness (Pa)

!                           ---

    character(len=255),allocatable    :: short_name(:)
    character(len=255),allocatable    :: long_name(:)
    character(len=255),allocatable    :: units(:)
    character(len=255), allocatable   :: table_info(:,:)
    type(ESMF_FieldBundle), pointer :: Bundle
    type(ESMF_Field), pointer       :: Field
    real(ESMF_KIND_R4), pointer :: ptr(:,:,:)
    integer :: i, j,im, jm, km, dims(7)
    integer :: linecount, columncount
                          __Iam__ ('SimpleBundleCreate')


!   Grid sizes
!   ----------
    call MAPL_GridGet(Grid, localCellCountPerDim = dims, __RC__)
    im = dims(1);  jm = dims(2);  km = dims(3)

!   Create an ESMF Bundle for holding variables
!   -------------------------------------------
    allocate(Bundle,__STAT__)
    Bundle = ESMF_FieldBundleCreate ( name=name, __RC__ )
    call ESMF_FieldBundleSet ( bundle, grid=Grid, __RC__ )

!   Parse rc and fill in short, long and units
!   ------------------------------------------
    call ESMF_ConfigGetDim(CF, linecount, columncount, label='aop_variables::', __RC__)
    if (columncount /= 3) then
     __raise__(MAPL_RC_ERROR,"columncount should be 3 (short_name, units and long_name)")
    end if
    allocate(table_info(linecount, columncount), __STAT__)
    allocate(short_name(linecount), units(linecount), long_name(linecount), __STAT__)
    call ESMF_ConfigFindLabel(CF, 'aop_variables::', __RC__)
    do i = 1, linecount
       call ESMF_ConfigNextLine(CF, __RC__) 
       do j = 1, columncount
           call ESMF_ConfigGetAttribute(CF, table_info(i,j),__RC__)
       enddo
       short_name(i) = table_info(i,1)
       units(i)      = table_info(i,2)
       long_name(i)  = table_info(i,3)       
    enddo

!   Add fields to Bundle
!   --------------------
    do i = 1, size(short_name)
       allocate(Field,ptr(im,jm,km), __STAT__)
       ptr = MAPL_UNDEF
       field = ESMF_FieldCreate(grid=Grid, &
               fArrayPtr = ptr, dataCopyFlag=ESMF_DATACOPY_REFERENCE, __RC__)
       call ESMF_AttributeSet(field, name='ShortName', value=short_name(i),  __RC__)
       call ESMF_AttributeSet(field, name='LongName', value=long_name(i),  __RC__)
       call ESMF_AttributeSet(field, name='Units', value=units(i),  __RC__)

       call MAPL_FieldBundleAdd(Bundle, Field=Field, __RC__)
    end do

!   Create the simple bundle
!   ------------------------
    self = MAPL_SimpleBundleCreate(Bundle, Levs=Levs, LevUnits=LevUnits, &
                                   ptop=ptop, delp=delp, __RC__)


   end Function GOCART2G_SimpleBundleCreate


!..........................................................................................

  Function GOCART2G_SimpleBundleRead (CF, rc_name, grid, time, verbose, only_vars, rc ) result(self)
!
!    Variant interface for MAPL_SimpleBundleRead.
!
    type(MAPL_SimpleBundle)                    :: self
    type(ESMF_Config),          intent(inout)  :: CF
    character(len=*),           intent(in)     :: rc_name
    type(ESMF_Grid),            intent(in)     :: Grid
    type(ESMF_Time),            intent(inout)  :: Time
    logical, OPTIONAL,          intent(in)     :: verbose
    character(len=*), optional,  intent(IN)    :: only_vars 
    integer, OPTIONAL,          intent(out)    :: rc
!                                ---
    character(len=256)         :: filename, template, expid, fname
    type(ESMF_Time)            :: Time_
    integer                    :: nymd, nhms, yy, mm, dd, h, m, s, fid, incSecs

    __Iam__ ('SimpleBundleRead')

    call ESMF_ConfigGetAttribute(CF, expid, Label='EXPID:', Default='unknown',__RC__ )
    call ESMF_ConfigGetAttribute(CF, filename, Label=trim(rc_name)//':',  __RC__ )
    fname = trim(rc_name)

       self = MAPL_SimpleBundleRead (filename, fname, Grid, Time, verbose, &
                                     ONLY_VARS=only_vars, expid=expid, __RC__ )

  end Function GOCART2G_SimpleBundleRead

!................................................

  subroutine GOCART2G_ESMFBundleWrite (bundle, CF, rc_name, wavelengths, Time, verbose, rc )
!
!   split the 4D fields (4th dim = wavelength) into nwavelengths bundle to write them 
!   in separate files

    type(ESMF_FieldBundle)                     :: bundle
    type(ESMF_Config),          intent(inout)  :: CF
    character(len=*),           intent(in)     :: rc_name
    real,                       intent(in)     :: wavelengths(:)
    type(ESMF_Clock),           intent(inout)  :: Time
    logical, OPTIONAL,          intent(in)     :: verbose
    integer, OPTIONAL,          intent(out)    :: rc
!                                ---
    character(len=256) :: filename
    integer :: j,i,itemCount,lm, n_wavelengths
    character(len=255), allocatable :: itemNames(:)
    type(ESMF_Field) :: field,new_field
    real, pointer :: ptr3d(:,:,:), ptr4d(:,:,:,:)
    type(ESMF_FieldBundle) :: new_bundle
    character(len=:), allocatable :: output_file
    character(len=4) :: istr
    type(ESMF_Grid) :: grid
    integer :: x1
    
    __Iam__ ('SimpleBundleWrite')

    call ESMF_ConfigGetAttribute(CF, filename, Label=trim(rc_name)//':',  __RC__ )

    call ESMF_FieldBundleGet(bundle,fieldCount=itemCount, __RC__)
    allocate(itemNames(itemCount))
    call ESMF_FieldBundleGet(bundle,fieldNameList=itemNames,grid=grid, __RC__)

    call ESMF_FieldBundleGet(bundle,trim(itemNames(1)),field=field, __RC__)
    call ESMF_FieldGet(field,0,farrayptr=ptr4d, __RC__)
    n_wavelengths = size(ptr4d,4)
    ! sanity check
    if (size(wavelengths) /= n_wavelengths) then
     __raise__(MAPL_RC_ERROR,"mismatch in the wavelengths count")
    end if        
    lm = size(ptr4d,3)

    new_bundle = ESMF_FieldBundleCreate ( __RC__ )
    call ESMF_FieldBundleSet ( new_bundle, grid=Grid, __RC__ )
    
    do i=1,itemCount
       call ESMF_FieldBundleGet(bundle,trim(itemNames(i)),field=field, __RC__)
       new_field = ESMF_FieldCreate(grid,typekind=ESMF_TYPEKIND_R4,ungriddedLBound=[1],ungriddedUBound=[lm],name=trim(itemNames(i)), __RC__)
       call ESMF_AttributeCopy(field,new_field, __RC__)
       call MAPL_FieldBundleAdd(new_bundle,new_field, __RC__)
    enddo

    do j=1,n_wavelengths
       do i=1,itemCount
          call ESMF_FieldBundleGet(bundle,trim(itemNames(i)),field=field, __RC__)
          call ESMF_FieldBundleGet(new_bundle,trim(itemNames(i)),field=new_field, __RC__)
          call ESMF_FieldGet(field,0,farrayptr=ptr4d, __RC__)
          call ESMF_FieldGet(new_field,0,farrayptr=ptr3d, __RC__)
          ptr3d=ptr4d(:,:,:,j)
       enddo
       x1 =  aint(wavelengths(j))
       write(istr,'(I4.4)') x1

       output_file='my_ouput_file_wavelength_'//istr//'.nc4'
       call MAPL_Write_Bundle(new_bundle,Time,output_file, __RC__)
    enddo
  end subroutine GOCART2G_ESMFBundleWrite

!-----------------
  subroutine GOCART2G_SimpleBundleWrite (self, CF, rc_name, Time, verbose, rc )
!
!    Variant interface for MAPL_SimpleBundleWrite.
!
    type(MAPL_SimpleBundle)                    :: self
    type(ESMF_Config),          intent(inout)  :: CF
    character(len=*),           intent(in)     :: rc_name
    type(ESMF_Time),            intent(inout)  :: Time
    logical, OPTIONAL,          intent(in)     :: verbose
    integer, OPTIONAL,          intent(out)    :: rc
!                                ---
    character(len=256) :: filename

    __Iam__ ('SimpleBundleWrite')

    call ESMF_ConfigGetAttribute(CF, filename, Label=trim(rc_name)//':',  __RC__ )

    call MAPL_SimpleBundleWrite ( self, filename, time, verbose, rc )

  end subroutine GOCART2G_SimpleBundleWrite




end module GOCART2G_SimpleBundleMod
