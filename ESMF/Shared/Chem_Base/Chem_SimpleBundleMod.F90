!
! This class extends MAPL_SimpleBundle with more Chemistry specific methods.
!
! Arlindo da Silva <arlindo.dasilva@nasa.gov>, November 2010
!----------------------------------------------------------------------------

#  include "MAPL_Generic.h"

module Chem_SimpleBundleMod

   use ESMF
   use MAPL
   use Chem_RegistryMod
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
   public Chem_SimpleBundleCreate
   public Chem_SimpleBundleRead
   public Chem_SimpleBundleWrite

   integer, parameter :: READ_ONLY=1

CONTAINS

!..........................................................................................

  Function Chem_SimpleBundleCreate (name, Reg, Grid, rc, &
                                    Levs, LevUnits,      &
                                    ptop, delp  )           result (self)

    type(MAPL_SimpleBundle)                    :: self
    character(len=*),           intent(in)     :: name
    type(Chem_Registry),        intent(in)     :: Reg
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

    type(ESMF_FieldBundle), pointer :: Bundle
    type(ESMF_Field), pointer       :: Field
    real(ESMF_KIND_R4), pointer :: ptr(:,:,:)
    integer :: i, im, jm, km, dims(7)

                          __Iam__ ('Chem_SimpleBundleCreate')


!   Grid sizes
!   ----------
    call MAPL_GridGet(Grid, localCellCountPerDim = dims, __RC__)
    im = dims(1);  jm = dims(2);  km = dims(3)

!   Create an ESMF Bundle for holding variables in Chem Registry
!   ------------------------------------------------------------
    allocate(Bundle,__STAT__)
    Bundle = ESMF_FieldBundleCreate ( name=name, __RC__ )
    call ESMF_FieldBundleSet ( bundle, grid=Grid, __RC__ )
   
!   Add fields
!   ----------
    do i = 1, Reg%nq
       allocate(Field,ptr(im,jm,km), __STAT__)
       ptr = MAPL_UNDEF
       field = ESMF_FieldCreate(name=Reg%vname(i), grid=Grid, &
               fArrayPtr = ptr, dataCopyFlag=ESMF_DATACOPY_REFERENCE, __RC__)
       call MAPL_FieldBundleAdd(Bundle, Field=Field, __RC__)
    end do

!   Create the simple bundle
!   ------------------------
    self = MAPL_SimpleBundleCreate(Bundle, Levs=Levs, LevUnits=LevUnits, &
                                   ptop=ptop, delp=delp, __RC__)

   end Function Chem_SimpleBundleCreate

!..........................................................................................

  Function Chem_SimpleBundleRead (CF, rc_name, grid, time, verbose, only_vars, rc ) result(self)
!
!    Variant interface for MAPL_SimpleBundleRead.
!
    type(MAPL_SimpleBundle)                    :: self
    type(ESMF_Config),          intent(inout)  :: CF
    character(len=*),           intent(in)     :: rc_name
    type(ESMF_Grid),            intent(in)     :: Grid
    type(ESMF_Time), OPTIONAL,  intent(inout)  :: Time
    logical, OPTIONAL,          intent(in)     :: verbose
    character(len=*), optional,  intent(IN)    :: only_vars 
    integer, OPTIONAL,          intent(out)    :: rc
!                                ---
    character(len=256)         :: filename, template, expid, fname
    type(ESMF_Time)            :: Time_
    integer                    :: nymd, nhms, yy, mm, dd, h, m, s, fid, incSecs

    __Iam__ ('Chem_SimpleBundleRead')

    call ESMF_ConfigGetAttribute(CF, expid, Label='EXPID:', Default='unknown',__RC__ )
    call ESMF_ConfigGetAttribute(CF, filename, Label=trim(rc_name)//':',  __RC__ )
    fname = trim(rc_name)

    if ( present(Time) ) then
       self = MAPL_SimpleBundleRead (filename, fname, Grid, Time, verbose, &
                                     ONLY_VARS=only_vars, expid=expid, __RC__ )
    else
       call GFIO_Open ( filename, READ_ONLY, fid, rc )
       if ( rc /= 0 ) return
       call GetBegDateTime ( fid, nymd, nhms, incSecs, rc )
       if ( rc /= 0 ) return
       call GFIO_Close ( fid, rc )
       if ( rc /= 0 ) return
       yy = nymd/10000; mm = (nymd-yy*10000) / 100; dd = nymd - (10000*yy + mm*100)
       h  = nhms/10000;  m = (nhms- h*10000) / 100;  s = nhms - (10000*h  +  m*100)
       call ESMF_TimeSet(Time_, yy=yy, mm=mm, dd=dd,  h=h,  m=m, s=s)
       self = MAPL_SimpleBundleRead (filename, fname, Grid, Time_, verbose, &
              ONLY_VARS=only_vars, expid=expid, __RC__ )
    end if

  end Function Chem_SimpleBundleRead

!................................................

  subroutine Chem_SimpleBundleWrite (self, CF, rc_name, Time, verbose, rc )
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

    __Iam__ ('Chem_SimpleBundleWrite')

    call ESMF_ConfigGetAttribute(CF, filename, Label=trim(rc_name)//':',  __RC__ )

    call MAPL_SimpleBundleWrite ( self, filename, time, verbose, rc )

  end subroutine Chem_SimpleBundleWrite

end module Chem_SimpleBundleMod
