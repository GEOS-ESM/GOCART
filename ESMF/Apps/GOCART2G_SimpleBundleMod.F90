!
! This class extends MAPL_SimpleBundle with some convenience methods.
! This would better implemented as an extension of MPL_SimpleBundle.
!
! Arlindo da Silva <arlindo.dasilva@nasa.gov>, March 2022
!----------------------------------------------------------------------------

#  include "MAPL_Generic.h"

module GOCART2G_SimpleBundleMod

   use ESMF
   use MAPL
   use RegistryMod
   use m_StrTemplate

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
   public GOCART2G_SimpleBundleCreate
   public GOCART2G_SimpleBundleRead
   public GOCART2G_SimpleBundleWrite

   integer, parameter :: READ_ONLY=1

CONTAINS

!..........................................................................................

  Function GOCART2G_SimpleBundleCreate (name, CF, rcname, Grid, rc, &
                               Levs, LevUnits, ptop, delp  )  result (self)

    type(MAPL_SimpleBundle)                    :: self
    character(len=*),           intent(in)     :: name
    type(ESMF_Config), intent(in)              :: CF     ! resources
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

    character(len=*)  :: short_name(:), long_name(:), units(:)
    type(ESMF_FieldBundle), pointer :: Bundle
    type(ESMF_Field), pointer       :: Field
    real(ESMF_KIND_R4), pointer :: ptr(:,:,:)
    integer :: i, im, jm, km, dims(7)

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
    ...
    
    
!   Add fields to Bundle
!   --------------------
    do i = 1, size(field_names)
       allocate(Field,ptr(im,jm,km), __STAT__)
       ptr = MAPL_UNDEF
       field = ESMF_FieldCreate(short(i), grid=Grid, &
            fArrayPtr = ptr, dataCopyFlag=ESMF_DATACOPY_REFERENCE, __RC__)
       ... add units and long attribute ...
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
