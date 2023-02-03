#define __SUCCESS__ 0
#define __VERIFY__(x) if(x/=0) then; if(present(rc)) rc=x; return; endif
#define __RC__ rc=status); __VERIFY__(status
#define __STAT__ stat=status); __VERIFY__(status
#define __IOSTAT__ iostat=status); __VERIFY__(status
#define __RETURN__(x) if (present(rc)) rc=x; return
#define __ASSERT__(expr) if(.not. (expr)) then; if (present(rc)) rc=-1; return; endif

#include "MAPL_Generic.h"
!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  CH4_GridCompMod --- CH4 Grid Component Class
!
! !INTERFACE:
!

module  CH4_GridCompMod

! USES:

use ESMF
use MAPL
use Chem_Mod         ! Chemistry Base Class
use Chem_StateMod        ! Chemistry State
use Chem_UtilMod         ! I/O
use m_inpak90        ! Resource file management
!   USE Henrys_law_ConstantsMod, ONLY: get_HenrysLawCts

implicit none

! !public types:
!
private
public  CH4_GridComp       ! Multiple instance CH4 object
public  CH4_GridComp1      ! Single instance CH4 object

!
! !PUBLIC MEMBER FUNCTIONS:
!

public  CH4_GridCompSetServices
public  CH4_GridCompInitialize
public  CH4_GridCompRun
public  CH4_GridCompFinalize

!
! !DESCRIPTION:
!
!  This module implements the CH4 Grid Component.
!
! !REVISION HISTORY:
!
!  24 Jun 2010 Nielsen:  First crack.
!  25 Oct 2012 Nielsen:  Added photolysis.
!  06 Oct 2017 Weir:     Fixing whitespace (STOP USING TABS IN FORTRAN!).
!
!EOP
!-------------------------------------------------------------------------

type CH4_GridComp1

   character(len=ESMF_MAXSTR) :: name        ! generic name of the package
   character(len=ESMF_MAXSTR) :: iname       ! instance name
   character(len=ESMF_MAXSTR) :: rcfilen     ! resource file name
   character(len=ESMF_MAXSTR) :: regionsString   ! Comma-delimited string of regions

   integer :: instance                 ! Instantiation number
   integer :: nymd_oh
   integer :: nymd_ch4
   integer :: BCnymd                   ! Date of last emissions/prodction read
   integer :: n_categ                  ! number of CH4 categories or sources that contri

   real, pointer :: regionMask(:,:) ! regional mask
   real, pointer :: CH4(:,:,:)      ! CH4 mixing ratio mol/mol
   real, pointer :: OHnd(:,:,:)     ! OH number density (cm^{-3})
   real, pointer :: Clnd(:,:,:)     ! Cl number density (cm^{-3})
   real, pointer :: O1Dnd(:,:,:)    ! O(1D) number density (cm^{-3})

   integer                                   :: n_sources
   character(len=ESMF_MAXSTR), allocatable   :: source_categs(:) ! which sources
   real, allocatable                         :: emis(:,:,:)      ! n_sources x i x j
   logical, allocatable                      :: diurnal_fire(:)  ! each source type can decide whether diurnal cycle is to be applied or not
   logical, allocatable                      :: pbl_inject(:)    ! ditto for PBL injection

   logical :: DebugIsOn     ! Run-time debug switch
   logical :: CH4FeedBack   ! Permit increments to CH4 from CH4 + hv => 2H2O + CO
   logical :: H2OFeedBack   ! Permit increments to   Q from CH4 + hv => 2H2O + CO
   logical :: chemistry     ! Should we do chemistry, such as oxidation and photolysis? Default true
   logical :: photolysis    ! Separate flag for photolysis
   logical :: in_total      ! Should this instance be included in a sum to get total CH4?
   real    :: szaCutoff     ! Largest solar zenith angle (degrees) allowed as daytime

end type CH4_GridComp1

type CH4_GridComp
   integer                            ::  n_inst   ! number of instances
   type(CH4_GridComp1), allocatable   ::  gcs(:)   ! instances
end type CH4_GridComp

real, parameter :: radToDeg = 180./MAPL_PI
real, parameter :: mwtCH4   = 16.0422

character(len=ESMF_MAXSTR), allocatable :: instances_in_total(:) ! list of instance names that make up total methane
character(len=ESMF_MAXSTR), allocatable :: all_instance_names(:) ! list of all instances from variable_table_CH4
character(len=*), parameter             :: rcbasen = 'CH4_GridComp'
character(len=1), parameter             :: func_sep = '_' ! The separator between the instance name (CH4wetland) and the function (EM) in exports

contains

subroutine CH4_GridCompSetServices(gc, chemReg, rc)

   type(ESMF_GridComp), intent(inout)  :: gc
   type(Chem_Registry), intent(inout)  :: chemReg
   integer,             intent(out)    :: rc

   character(len=255) :: name
   character(len=ESMF_MAXSTR) :: dummy_str

   integer            :: status
   character(len=*), parameter :: Iam = "CH4::GridCompSetServices"

   integer :: ier,n,i

   type(ESMF_Config) :: cfg

   !  Load resource file
   !  ------------------
   cfg = ESMF_ConfigCreate(__RC__)

   call ESMF_ConfigLoadFile(cfg, trim(rcbasen)//'.rc', __RC__) ! reading CH4_GridComp.rc

   ! Debug :: print what is already in Chem_Registry
   if (MAPL_AM_I_ROOT()) then
      do i = chemReg%i_CH4, chemReg%j_CH4
         write(*,'("    ", a, ": from Chem_Registry, tracer ", i0, " is methane instance ", i0, " with name ", a)') &
            trim(Iam), i, (i-chemReg%i_CH4+1), trim(chemReg%vname(i))
      end do
   end if
   ! End debug

   allocate(all_instance_names(chemReg%i_CH4:chemReg%j_CH4))
   do i = chemReg%i_CH4, chemReg%j_CH4
      call CH4_GridCompSetServices1_(gc, chemReg%vname(i), __RC__)
      all_instance_names(i) = trim(chemReg%vname(i))
   end do

   ! Now read which instances need to be summed to make the total methane
   n = ESMF_ConfigGetLen(cfg, label='CH4_total_components:', __RC__)
   allocate(instances_in_total(n))
   call ESMF_ConfigFindLabel(cfg, 'CH4_total_components:', __RC__)
   do i = 1, n
      call ESMF_ConfigGetAttribute(cfg, dummy_str, __RC__)
      instances_in_total(i) = trim(dummy_str)
   end do

   call MAPL_AddImportSpec(GC,           &
        SHORT_NAME = 'CH4_regionMask',   &
        LONG_NAME  = 'source species'  , &
        UNITS      = '1',                &
        DIMS       = MAPL_DimsHorzOnly,  &
        VLOCATION  = MAPL_VLocationNone, &
        RESTART    = MAPL_RestartSkip,   &
        __RC__)

   call MAPL_AddExportSpec(GC,  &
      SHORT_NAME         = 'CH4',  &
      LONG_NAME          = 'CH4 total mole fraction',  &
      UNITS              = 'mol mol-1', &
      DIMS               = MAPL_DimsHorzVert,    &
      VLOCATION          = MAPL_VLocationCenter,    &
      __RC__)

   call MAPL_AddExportSpec(GC,  &
      SHORT_NAME         = 'CH4DRY',  &
      LONG_NAME          = 'CH4 total dry-air mole fraction',  &
      UNITS              = 'mol mol-1', &
      DIMS               = MAPL_DimsHorzVert,    &
      VLOCATION          = MAPL_VLocationCenter,    &
      __RC__)

   RETURN_(ESMF_SUCCESS)

end subroutine CH4_GridCompSetServices

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  CH4_GridCompInitialize --- Initialize CH4_GridComp
!
! !INTERFACE:
!

subroutine CH4_GridCompInitialize ( gcCH4, w_c, impChem, expChem, nymd, nhms, cdt, rc )

   !USES:

   implicit none

   !INPUT PARAMETERS:

   TYPE(Chem_Bundle), intent(in) :: w_c        ! Chemical tracer fields
   integer, intent(in) :: nymd, nhms           ! time
   real,    intent(in) :: cdt              ! chemical timestep (secs)


   ! OUTPUT PARAMETERS:

   type(CH4_GridComp), intent(inout) :: gcCH4   ! Grid Component
   type(ESMF_State), intent(inout)  :: impChem  ! Import State
   type(ESMF_State), intent(inout)  :: expChem  ! Export State
   integer, intent(out) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 -

   ! !DESCRIPTION: Initializes the CH4 Grid Component. Multiple instance
   !               version.
   !
   ! !REVISION HISTORY:
   !
   !  24 Jun 2010 Nielsen:  First crack.
   !  25 Oct 2012 Nielsen:  Added photolysis.
   !
   !EOP
   !-------------------------------------------------------------------------

   character(len=*), parameter :: Iam = 'CH4::GridCompInitialize'
   character(len=esmf_maxstr) :: name

   integer :: i, n, status, min_i

   gcCH4%n_inst = size(all_instance_names)
   min_i = lbound(all_instance_names, dim=1) ! i_CH4

   !  Next allocate necessary memory
   !  ------------------------------
   allocate ( gcCH4%gcs(gcCH4%n_inst), STAT=status ) ! each instance gets a gcCH4%gcs
   VERIFY_(status)

   !  Record name of each instance
   !  ----------------------------
   do i = 1, gcCH4%n_inst
      gcCH4%gcs(i)%rcfilen = trim(rcbasen)//'.rc' ! Experiment to read all rc keys from one file
      gcCH4%gcs(i)%instance = i              ! instance number
      gcCH4%gcs(i)%iname = trim(all_instance_names(min_i+i-1))
   end do

   !  Next initialize each instance
   !  -----------------------------
   do i = 1, gcCH4%n_inst
      if (MAPL_AM_I_ROOT()) then
         write(*,'("    ", a, ": initializing instance ", i0, ", ", a)') trim(Iam), gcCH4%gcs(i)%instance, TRIM(gcCH4%gcs(i)%iname)
      end if
      call CH4_SingleInstance_ ( CH4_GridCompInitialize1_, i, &
                                gcCH4%gcs(i), w_c, impChem, expChem,  &
                                nymd, nhms, cdt, status )
      VERIFY_(status)
   end do

   !!  Get Henrys Law cts for the parameterized convective wet removal
   !!  -----------------------------------------------------------
   !   CALL get_HenrysLawCts('CH4',c1,c2,c3,c4)
   !   w_c%reg%Hcts(1,w_c%reg%i_CH4 : w_c%reg%j_CH4)=c1
   !   w_c%reg%Hcts(2,w_c%reg%i_CH4 : w_c%reg%j_CH4)=c2
   !   w_c%reg%Hcts(3,w_c%reg%i_CH4 : w_c%reg%j_CH4)=c3
   !   w_c%reg%Hcts(4,w_c%reg%i_CH4 : w_c%reg%j_CH4)=c4

end subroutine CH4_GridCompInitialize

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  CH4_GridCompRun --- Run CH4_GridComp
!
! !INTERFACE:
!

subroutine CH4_GridCompRun ( gcCH4, w_c, impChem, expChem, nymd, nhms, cdt, rc )

   !USES:

   implicit none

   !INPUT PARAMETERS:

   type(Chem_Bundle), intent(in) :: w_c        ! Chemical tracer fields
   integer, intent(in) :: nymd, nhms           ! time
   real,    intent(in) :: cdt              ! chemical timestep (secs)


   !OUTPUT PARAMETERS:

   type(CH4_GridComp), intent(inout) :: gcCH4   ! Grid Component
   type(ESMF_State),   intent(inout) :: impChem ! Import State
   type(ESMF_State),   intent(inout) :: expChem ! Export State
   integer, intent(out)              :: rc      ! Error return code:
                                                !  0 - all is well
                                                !  1 -

   ! !DESCRIPTION: Runs the CH4 Grid Component. Multiple instance
   !               version.
   !
   ! !REVISION HISTORY:
   !
   !  24 Jun 2010 Nielsen:  First crack.
   !  25 Oct 2012 Nielsen:  Added photolysis.
   !
   !EOP
   !-------------------------------------------------------------------------

   character(len=*), parameter :: Iam = 'CH4::GridCompRun'
   INTEGER :: i, status

   real, pointer, dimension(:,:,:) :: CH4_total => null() ! total wet-air methane mole fraction ! Sourish
   real, pointer, dimension(:,:,:) :: CH4_dry => null() ! total dry-air methane mole fraction ! Sourish
   real, pointer, dimension(:,:,:) :: qtot => null()
   integer :: i1, i2, j1, j2, km

   do i = 1, gcCH4%n_inst
      call CH4_SingleInstance_ ( CH4_GridCompRun1_, i, &
                                gcCH4%gcs(i), w_c, impChem, expChem, &
                                nymd, nhms, cdt, status )
      VERIFY_(status)
   end do

   call MAPL_GetPointer ( expChem, CH4_total, 'CH4', __RC__ )
   if (associated(CH4_total)) then
      i1 = w_c%grid%i1
      i2 = w_c%grid%i2
      j1 = w_c%grid%j1
      j2 = w_c%grid%j2
      km = w_c%grid%km
      CH4_total = 0.0
      do i = 1, gcCH4%n_inst
         if (gcCH4%gcs(i)%in_total) then
            CH4_total = CH4_total + w_c%qa(i)%data3d(i1:i2,j1:j2,1:km)
         end if
      end do

      ! We need CH4_total for CH4_dry, so keep it within this if block
      call MAPL_GetPointer(impChem, qtot, 'QTOT', __RC__) ! Sourish
      call MAPL_GetPointer(expChem, CH4_dry, 'CH4DRY', __RC__)
      if (associated(CH4_dry)) &
         CH4_dry(i1:i2,j1:j2,1:km) = CH4_total / (1. - qtot(i1:i2,j1:j2,1:km))
   end if

end subroutine CH4_GridCompRun


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  CH4_GridCompFinalize --- Initialize CH4_GridComp
!
! !INTERFACE:
!

subroutine CH4_GridCompFinalize ( gcCH4, w_c, impChem, expChem, nymd, nhms, cdt, rc )

   implicit none

   ! INPUT PARAMETERS:

   type(Chem_Bundle), intent(in) :: w_c        ! Chemical tracer fields
   integer, intent(in) :: nymd, nhms           ! time
   real,    intent(in) :: cdt              ! chemical timestep (secs)


   ! OUTPUT PARAMETERS:

   type(CH4_GridComp), intent(inout)   :: gcCH4   ! Grid Component
   type(ESMF_State), intent(inout)     :: impChem ! Import State
   type(ESMF_State), intent(inout)     :: expChem ! Export State
   integer, intent(out)                ::  rc     ! Error return code:
                                                  !  0 - all is well
                                                  !  1 -

   ! !DESCRIPTION: Finalizes the CH4 Grid Component. Multiple instance
   !               version.
   !
   ! !REVISION HISTORY:
   !
   !  27Feb2008  da Silva  Introduced multiple instances
   !
   !EOP
   !-------------------------------------------------------------------------

   character(len=*), parameter :: Iam = 'CH4::GridCompFinalize'
   integer :: i, status

   do i = 1, gcCH4%n_inst
      call CH4_SingleInstance_ ( CH4_GridCompFinalize1_, i, &
                                gcCH4%gcs(i), w_c, impChem, expChem, &
                                nymd, nhms, cdt, status )
      VERIFY_(status)
   end do

   deallocate ( gcCH4%gcs, stat=status )
   gcCH4%n_inst = -1

   deallocate(instances_in_total, all_instance_names)

end subroutine CH4_GridCompFinalize

subroutine CH4_GridCompSetServices1_(  gc, iname, rc)

   type(ESMF_GridComp), intent(inout) :: GC
   character(len=*),    intent(in)    :: iname
   integer,             intent(out)   :: rc

   integer :: Status, n_sources, i
   character(len=*), parameter :: Iam = "CH4::GridCompSetServices1_"
   type(ESMF_Config) :: cfg
   character(len=ESMF_MAXSTR) :: emis_varname

   ! Figure out where the emissions are coming from and add imports
   cfg = ESMF_ConfigCreate(__RC__)
   call ESMF_ConfigLoadFile(cfg, trim(rcbasen)//'.rc', __RC__) ! reading CH4_GridComp.rc
   call ESMF_ConfigFindLabel(cfg, 'emissions.'//trim(iname)//':', rc=status)
   if (status .ne. ESMF_SUCCESS) then
      emis_varname = 'emis_'//trim(iname)
      call MAPL_AddImportSpec(GC, &
         SHORT_NAME = trim(emis_varname), &
         LONG_NAME  = 'source species'  , &
         UNITS      = 'kg CH4 m-2 s-1',     &
         DIMS       = MAPL_DimsHorzOnly, &
         VLOCATION  = MAPL_VLocationNone, &
         RESTART    = MAPL_RestartSkip,     &
         __RC__)
   else
      n_sources = ESMF_ConfigGetLen(cfg, label='emissions.'//trim(iname)//':', __RC__)
      call ESMF_ConfigFindLabel(cfg, 'emissions.'//trim(iname)//':', __RC__)
      do i = 1, n_sources
         call ESMF_ConfigGetAttribute(cfg, emis_varname, __RC__)
         call MAPL_AddImportSpec(GC, &
            SHORT_NAME = trim(emis_varname), &
            LONG_NAME  = 'source species'  , &
            UNITS      = 'kg CH4 m-2 s-1',     &
            DIMS       = MAPL_DimsHorzOnly, &
            VLOCATION  = MAPL_VLocationNone, &
            RESTART    = MAPL_RestartSkip,     &
            __RC__)
      end do
   end if

   call MAPL_AddImportSpec(GC, &
      SHORT_NAME = 'CH4_oh', &
      LONG_NAME  = 'source species'  , &
      UNITS      = 'molecules cm-3',     &
      DIMS       = MAPL_DimsHorzVert,    &
      VLOCATION  = MAPL_VLocationCenter, &
      __RC__)
   call MAPL_AddImportSpec(GC,             &
      SHORT_NAME = 'CH4_cl',      &
      LONG_NAME  = 'source species',     &
      UNITS      = 'molecules cm-3',     &
      DIMS       = MAPL_DimsHorzVert,    &
      VLOCATION  = MAPL_VLocationCenter, &
      __RC__)
   call MAPL_AddImportSpec(GC,             &
      SHORT_NAME = 'CH4_o1d',     &
      LONG_NAME  = 'source species',     &
      UNITS      = 'molecules cm-3',     &
      DIMS       = MAPL_DimsHorzVert, &
      VLOCATION  = MAPL_VLocationCenter, &
      __RC__)

   call MAPL_AddExportSpec(GC,  &
      SHORT_NAME         = trim(iname)//func_sep//'EM',  &
      LONG_NAME          = 'Emission of instance '//trim(iname),  &
      UNITS              = 'kg m-2 s-1', &
      DIMS               = MAPL_DimsHorzOnly,    &
      VLOCATION          = MAPL_VLocationNone,    &
      __RC__)
   call MAPL_AddExportSpec(GC,  &
      SHORT_NAME         = trim(iname)//func_sep//'LS',  &
      LONG_NAME          = 'Loss due to oxidation of instance '//trim(iname),  &
      UNITS              = 'kg m-2 s-1', &
      DIMS               = MAPL_DimsHorzOnly,    &
      VLOCATION          = MAPL_VLocationNone,    &
      __RC__)
   call MAPL_AddExportSpec(GC,  &
      SHORT_NAME         = trim(iname)//func_sep//'JL',  &
      LONG_NAME          = 'Loss due to photolysis of instance '//trim(iname),  &
      UNITS              = 'kg m-2 s-1', &
      DIMS               = MAPL_DimsHorzOnly,    &
      VLOCATION          = MAPL_VLocationNone,    &
      __RC__)
   call MAPL_AddExportSpec(GC,  &
      SHORT_NAME         = trim(iname)//func_sep//'SC',  &
      LONG_NAME          = 'Surface layer mole fraction of instance '//trim(iname),  &
      UNITS              = 'mol mol-1 (dry air)', &
      DIMS               = MAPL_DimsHorzOnly,    &
      VLOCATION          = MAPL_VLocationNone,    &
      __RC__)
   call MAPL_AddExportSpec(GC,  &
      SHORT_NAME         = trim(iname)//func_sep//'CL',  &
      LONG_NAME          = 'Column average mole fraction of instance '//trim(iname),  &
      UNITS              = 'mol mol-1 (dry air)', &
      DIMS               = MAPL_DimsHorzOnly,    &
      VLOCATION          = MAPL_VLocationNone,    &
      __RC__)

   RETURN_(ESMF_SUCCESS)

end subroutine CH4_GridCompSetServices1_

!--------------------------------------------------------------------------

!                      Single Instance Methods

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  CH4_GridCompInitialize --- Initialize CH4_GridComp
!
! !INTERFACE:
!

subroutine CH4_GridCompInitialize1_ ( gcCH4, w_c, impChem, expChem, nymd, nhms, cdt, rc )

   implicit none

   ! INPUT PARAMETERS:

   type(Chem_Bundle), intent(in) :: w_c        ! Chemical tracer fields
   integer, intent(in) :: nymd, nhms           ! time
   real,    intent(in) :: cdt              ! chemical timestep (secs)


   ! OUTPUT PARAMETERS:

   type(CH4_GridComp1), intent(inout)  :: gcCH4    ! Grid Component
   type(ESMF_State), intent(inout)     :: impChem  ! Import State
   type(ESMF_State), intent(inout)     :: expChem  ! Export State
   integer, intent(out)                ::  rc      ! Error return code:
                                                   !  0 - all is well
                                                   !  1 -

   ! !DESCRIPTION: Initializes the CH4 Grid Component. It primarily sets
   !               the import state for each active constituent package.
   !
   ! !REVISION HISTORY:
   !
   !  18Sep2003 da Silva  First crack.
   !  31May2005  Nielsen  Mods for 7 CH4 bins, 5 region masks
   !  04Nov2005     Bian  CO tagged to 4 regions
   !                      (global, North America, South America, and Africa)
   !                      for CR-AVE
   !  25Oct2012  Nielsen  Added photolysis.
   !
   !EOP
   !-------------------------------------------------------------------------

   character(len=*), parameter  :: Iam = 'CH4::GridCompInitialize1_'

   character(len=ESMF_MAXSTR)   :: rcfilen, dummy_str
   type(ESMF_Config) :: cfg

   integer :: j, n, comps_total, i, status
   integer :: i1, i2, im, j1, j2, jm, km
   integer :: nTimes, begTime, incSecs
   integer :: nbeg, nend, nymd1, nhms1
   logical :: NoRegionalConstraint

   real :: limitN, limitS
   real, allocatable :: var2D(:,:)
   real, pointer     :: ptr2d(:,:) => null()
   real :: dummy_float
   logical :: dummy_bool
   integer :: dummy_int
   character(len=ESMF_MAXSTR) :: dummy_string

   rcfilen = gcCH4%rcfilen
   gcCH4%name = 'GEOS-5/GOCART Parameterized CH4 Package'
   gcCH4%BCnymd = -1

   !  Initialize local variables
   !  --------------------------
   rc = 0
   i1 = w_c%grid%i1
   i2 = w_c%grid%i2
   im = w_c%grid%im

   j1 = w_c%grid%j1
   j2 = w_c%grid%j2
   jm = w_c%grid%jm

   km = w_c%grid%km

   nbeg  = w_c%reg%i_CH4
   nend  = w_c%reg%j_CH4

   !  It requires 1 bin
   !  -----------------
   if ( nbeg /= nend ) then
      if (MAPL_AM_I_ROOT()) write(*, '(a, ": Must have only 1 bin at the single instance level")') trim(Iam)
      status = 1
      VERIFY_(status)
   end if

   !  Allocate memory, etc
   !  --------------------
   allocate( &
            !gcCH4%eCH4(i1:i2,j1:j2),       &
            gcCH4%regionMask(i1:i2,j1:j2), &
            gcCH4%OHnd(i1:i2,j1:j2,km),    &
            gcCH4%Clnd(i1:i2,j1:j2,km),    &
            gcCH4%O1Dnd(i1:i2,j1:j2,km), STAT=status )
   VERIFY_(status)

   !  Load resource file
   !  ------------------
   cfg = ESMF_ConfigCreate()
   call ESMF_ConfigLoadFile(cfg, trim(rcfilen), __RC__)

   ! Is this instance part of 'total CH4' diagnostic?
   !  -------------------------------------------------
   gcCH4%in_total = .false.
   comps_total = size(instances_in_total)
   do i = 1, comps_total
      dummy_int = index(instances_in_total(i), trim(gcCH4%iname))
      if (dummy_int > 0) gcCH4%in_total = .true.
   end do
   if (MAPL_AM_I_ROOT()) then
      if (gcCH4%in_total) write(*,'("    ", a, ": total methane contains component ", a)') trim(Iam), trim(gcCH4%iname)
      if (.not. gcCH4%in_total) write(*,'("    ", a, ": total methane does not contain component ", a)') trim(Iam), trim(gcCH4%iname)
   end if

   ! Which emission categories will contribute to this instance? By default assign 'emis_CH4instance' if not specified
   !  -------------------------------------------------
   call ESMF_ConfigFindLabel(cfg, 'emissions.'//trim(gcCH4%iname)//':', rc=status)
   if (status .ne. ESMF_SUCCESS) then
      gcCH4%n_sources = 1
      allocate(gcCH4%source_categs(1))
      gcCH4%source_categs(1) = 'emis_'//trim(gcCH4%iname)
      if (MAPL_AM_I_ROOT()) write(*,'("    ", a, ": emission from ", a, " to be added to ", a)') trim(Iam), trim(gcCH4%source_categs(1)), trim(gcCH4%iname)
   else
      gcCH4%n_sources = ESMF_ConfigGetLen(cfg, label='emissions.'//trim(gcCH4%iname)//':', __RC__)
      allocate(gcCH4%source_categs(gcCH4%n_sources))
      call ESMF_ConfigFindLabel(cfg, 'emissions.'//trim(gcCH4%iname)//':', __RC__)
      do i = 1, gcCH4%n_sources
         call ESMF_ConfigGetAttribute(cfg, dummy_str, __RC__)
         gcCH4%source_categs(i) = trim(dummy_str)
         if (MAPL_AM_I_ROOT()) write(*,'("    ", a, ": emission from ", a, " to be added to ", a)') trim(Iam), trim(gcCH4%source_categs(i)), trim(gcCH4%iname)
      end do
   end if
   allocate(gcCH4%emis(gcCH4%n_sources, i1:i2, j1:j2))
   allocate(gcCH4%diurnal_fire(gcCH4%n_sources))
   allocate(gcCH4%pbl_inject(gcCH4%n_sources))
   do i = 1, gcCH4%n_sources
      ! Should we apply a fire-like diurnal cycle for this source?
      call ESMF_ConfigGetAttribute(cfg, value=gcCH4%diurnal_fire(i), label=trim(gcCH4%source_categs(i))//'.diurnal_fire:', &
         default=.false., __RC__)
      if (MAPL_AM_I_ROOT() .and. gcCH4%diurnal_fire(i)) &
         write(*,'("    ", a, ": fire-like diurnal cycle applied for source ", a, " in instance ", a)') trim(Iam), trim(gcCH4%source_categs(i)), trim(gcCH4%iname)
      ! Should we distribute this source throughout the PBL instead of at the surface?
      call ESMF_ConfigGetAttribute(cfg, value=gcCH4%pbl_inject(i), label=trim(gcCH4%source_categs(i))//'.pbl_injection:', &
         default=.false., __RC__)
      if (MAPL_AM_I_ROOT() .and. gcCH4%pbl_inject(i)) &
         write(*,'("    ", a, ": emission of ", a, " in category ", a, " to be spread through the PBL")') trim(Iam), trim(gcCH4%source_categs(i)), trim(gcCH4%iname)
   end do

   !  Maximum allowed solar zenith angle for "daylight"
   !  -------------------------------------------------
   call ESMF_ConfigGetAttribute(cfg, value=dummy_float, label='solar_ZA_cutoff:', __RC__)
   call ESMF_ConfigGetAttribute(cfg, value=gcCH4%szaCutoff, label='solar_ZA_cutoff.'//trim(gcCH4%iname)//':', default=dummy_float, __RC__)

   !  Run-time debug switch
   !  ---------------------
   call ESMF_ConfigGetAttribute(cfg, value=dummy_bool, label='DEBUG:', default=.false., __RC__)
   call ESMF_ConfigGetAttribute(cfg, value=gcCH4%DebugIsOn, label='DEBUG.'//trim(gcCH4%iname)//':', default=dummy_bool, __RC__)

   !  Methane photolysis feedback switch
   !  ----------------------------------
   call ESMF_ConfigGetAttribute(cfg, value=dummy_bool, label='CH4_Feedback:', default=.false., __RC__)
   call ESMF_ConfigGetAttribute(cfg, value=gcCH4%CH4FeedBack, label='CH4_Feedback.'//trim(gcCH4%iname)//':', default=dummy_bool, __RC__)

   !  Water vapor feedback switch
   !  ---------------------------
   call ESMF_ConfigGetAttribute(cfg, value=dummy_bool, label='H2O_Feedback:', default=.false., __RC__)
   call ESMF_ConfigGetAttribute(cfg, value=gcCH4%H2OFeedBack, label='H2O_Feedback.'//trim(gcCH4%iname)//':', default=dummy_bool, __RC__)

   ! Should we do chemistry?
   !  ----------------------
   call ESMF_ConfigGetAttribute(cfg, value=dummy_bool, label='do_chemistry:', default=.true., __RC__)
   call ESMF_ConfigGetAttribute(cfg, value=gcCH4%chemistry, label='do_chemistry.'//trim(gcCH4%iname)//':', default=dummy_bool, __RC__)
   if (MAPL_AM_I_ROOT() .and. (.not. gcCH4%chemistry)) write(*,'("    Chemistry turned off for ", a)') trim(gcCH4%iname)

   ! Should we do photolysis?
   !  ----------------------
   call ESMF_ConfigGetAttribute(cfg, value=dummy_bool, label='do_photolysis:', default=.false., __RC__)
   call ESMF_ConfigGetAttribute(cfg, value=gcCH4%photolysis, label='do_photolysis.'//trim(gcCH4%iname)//':', default=dummy_bool, __RC__)
   if (MAPL_AM_I_ROOT() .and. gcCH4%photolysis) write(*,'("    Photolytic loss turned on for ", a)') trim(gcCH4%iname)

   call MAPL_GetPointer(impChem,ptr2D,'CH4_regionMask', __RC__)

   !  Grab the region string.
   !  -----------------------
   call ESMF_ConfigGetAttribute(cfg, value=dummy_string, label='CH4_regions_indices:', __RC__)
   call ESMF_ConfigGetAttribute(cfg, value=gcCH4%regionsString, label='CH4_regions_indices.'//trim(gcCH4%iname)//':', default=dummy_string, __RC__)

   !  Is this instantiation a global case?
   !  -----------------------------------
   if (gcCH4%regionsString(1:2) == "-1") then
      NoRegionalConstraint = .true.
   else
      select case (ESMF_UtilStringLowerCase(gcCH4%regionsString(1:2)))
      case ("gl")
         NoRegionalConstraint = .true.
      case ("al")
         NoRegionalConstraint = .true.
      case default
         NoRegionalConstraint = .false.
      end select
   end if

   !  Set regionsString to "-1" for the global case
   !  ---------------------------------------------
   if (NoRegionalConstraint) gcCH4%regionsString = "-1"

   if (MAPL_AM_I_ROOT()) then
      if (NoRegionalConstraint) then
         write(*,'(a, ": This instantiation has no regional constraints.")') trim(Iam)
      else
         write(*,'(a, ": This instantiation is regionally constrained.")') trim(Iam)
         write(*,'(a, ": List of region numbers included: ", a)') trim(Iam), trim(gcCH4%regionsString)
      end if
   end if

end subroutine CH4_GridCompInitialize1_

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  CH4_GridCompRun
!
! !INTERFACE:
!

subroutine Calc_emissions ( gcCH4, w_c, impChem, nhms, cdt, flux_3d, rc )

   implicit none

   ! Input variables
   type(CH4_GridComp1), intent(in)  :: gcCH4   ! Grid Component
   type(Chem_Bundle), intent(in)    :: w_c     ! Chemical tracer fields
   type(ESMF_State), intent(inout)  :: impChem ! Import State
   integer, intent(in)              :: nhms    ! time
   real, intent(in)                 :: cdt     ! chemical timestep (secs)

   ! output variables
   integer, intent(out)    :: rc ! Error return code
   real, dimension(:,:,:)  :: flux_3d

   ! local variables
   character(len=*), parameter :: Iam = 'CH4::Calc_emissions'

   integer              ::  i, j, k, kt, minkPBL, i1, i2, j1, j2, km, status
   integer, allocatable :: index(:)

   real, allocatable, dimension(:,:)   :: pblLayer, sfcFlux
   real, allocatable, dimension(:,:,:) :: fPBL, mf_to_be_added, mf_current
   real, pointer, dimension(:,:,:)     :: zle => null()
   real, pointer, dimension(:,:)       :: pblh => null()

   ! local variables and shortcuts
   rc = 0
   i1 = w_c%grid%i1
   i2 = w_c%grid%i2
   j1 = w_c%grid%j1
   j2 = w_c%grid%j2
   km = w_c%grid%km

   allocate(sfcFlux(i1:i2,j1:j2), stat=rc)  ! emissions
   allocate(fPBL(i1:i2,j1:j2,1:km), stat=rc)  ! partitioning of BB

   ! Import PBL height and layer thickness
   call MAPL_GetPointer(impChem, pblh, 'ZPBL', __RC__)
   CALL MAPL_GetPointer(impChem, zle, 'ZLE', __RC__)

   ! Find the layer that contains the PBL
   ! Layer thicknesses are ZLE(:,:,0:km)
   ! ------------------------------------
   allocate(index(0:km), STAT=rc)
   allocate(pblLayer(i1:i2,j1:j2), STAT=rc)
   do j = j1,j2
      do i = i1,i2
         index(0:km) = 0
         where(zle(i,j,0:km) - zle(i,j,km) > pblh(i,j)) index(0:km) = 1
         pblLayer(i,j) = sum(index)
      end do
   end do
   minkPBL = minval(pblLayer)

   ! Determine partitioning fraction based on layer thicknesses
   ! ----------------------------------------------------------
   fPBL(i1:i2,j1:j2,1:km) = 0.
   do j = j1,j2
      do i = i1,i2
         kt = pblLayer(i,j)
         do k = kt,km
            fPBL(i,j,k) = (zle(i,j, k-1) - zle(i,j,k)) / (zle(i,j,kt-1) - zle(i,j,km))
         end do
      end do
   end do

   deallocate(index, pblLayer, STAT=rc)

   ! Establish range of layers on which to work
   ! ------------------------------------------
   kt = minkPBL ! kt means 'k top'

   ! Accumulate total 3D flux to be added in flux_3d
   flux_3d = 0.0
   do i = 1, gcCH4%n_sources
      ! Diurnal cycle needed for this source?
      if (gcCH4%diurnal_fire(i)) then
         call Chem_BiomassDiurnal(sfcFlux(:,:), gcCH4%emis(i,:,:), w_c%grid%lon(:,:)*radToDeg, w_c%grid%lat(:,:)*radToDeg, nhms, cdt)
      else
         sfcFlux(:,:) = gcCH4%emis(i,:,:)
      end if
      ! Do we need to spread this throughout the PBL?
      if (gcCH4%pbl_inject(i)) then
         do k = kt, km
            flux_3d(:,:,k) = flux_3d(:,:,k) + sfcFlux(:,:) * fPBL(:,:,k)
         end do
      else
         flux_3d(:,:,km) = flux_3d(:,:,km) + sfcFlux(:,:)
      end if
   end do

   deallocate(fPBL, sfcFlux, stat=rc)

end subroutine Calc_emissions

subroutine Calc_photolysis_rates(gcCH4, w_c, impChem, photJ, rc)
   !-------------------------------------------------------------------------
   ! Borrowed from meso_phot.F of StratChem, where number densities are cgs [cm^{-3}]

   implicit none

   ! Input variables
   type(CH4_GridComp1), intent(in)     :: gcCH4   ! Grid Component
   type(Chem_Bundle), intent(in)       :: w_c     ! Chemical tracer fields
   type(ESMF_State), intent(inout)     :: impChem ! Import State

   ! Output variables
   integer, intent(out)                :: rc
   real, dimension(:,:,:), intent(out) :: photJ

   real, allocatable :: o2Column(:,:,:), ndwet(:,:,:), cellDepth(:,:,:)
   real, allocatable :: SZARad(:,:), SZADeg(:,:), sinSZA(:,:)
   real, allocatable :: zgrz(:,:), sfaca(:,:), arg(:,:)
   real, pointer, dimension(:,:,:) ::  rhowet => null()
   real, pointer, dimension(:,:,:) ::  ple => null()

   real, parameter :: wavel = 1215.7
   real, parameter :: O2xs  = 1.000E-20
   real, parameter :: CH4xs = 2.000E-17
   real, parameter :: sflux = 4.006E+11

   ! Constants for Chapman function at high solar zenith angle
   ! ---------------------------------------------------------
   real, parameter :: hbar = 6.79
   real, parameter :: zbar = 30.0
   real, parameter :: r0   = 6.371E+03
   real, parameter :: zp   = 60.0

   real, parameter :: d1 = 1.060693
   real, parameter :: d2 = 0.55643831
   real, parameter :: d3 = 1.0619896
   real, parameter :: d4 = 1.7245609
   real, parameter :: d5 = 0.56498823
   real, parameter :: d6 = 0.06651874

   real, parameter :: O2Abv80km = 7.072926E+19 ![cm^{-2}]
   real, parameter :: O2VMR = 0.20946

   real :: b, f, r, s
   integer :: status, i1, i2, j1, j2, km, i, j, k

   character(len=*), parameter :: Iam = "CH4::Calc_photolysis_rates"

   rc = 0
   photJ(:,:,:) = 0

   i1 = w_c%grid%i1
   i2 = w_c%grid%i2
   j1 = w_c%grid%j1
   j2 = w_c%grid%j2
   km = w_c%grid%km

   ! Get pointers
   call MAPL_GetPointer(impChem, rhowet, 'AIRDENS', __RC__)
   call MAPL_GetPointer(impChem, ple, 'PLE', __RC__)

   b = sqrt(0.50*r0/hbar)

   ! O2 overhead number density profile [cm^{-2}]
   ! --------------------------------------------
   allocate(O2Column(i1:i2,j1:j2,1:km), ndwet(i1:i2,j1:j2,km), cellDepth(i1:i2,j1:j2,km), STAT=status)
   VERIFY_(status)

   !  Wet-air number density
   ndwet(i1:i2,j1:j2,1:km) = rhowet(i1:i2,j1:j2,1:km)*MAPL_AVOGAD/MAPL_AIRMW
   !  Cell depth
   do k=1,km
      cellDepth(i1:i2,j1:j2,k) = (ple(i1:i2,j1:j2,k)-ple(i1:i2,j1:j2,k-1)) / (rhowet(i1:i2,j1:j2,k)*MAPL_GRAV)
   end do

   f = O2VMR*5.00E-05
   O2Column(:,:,1) = O2Abv80km+cellDepth(:,:,1)*ndwet(:,:,1)*f

   do k = 2,km
      O2Column(:,:,k) = O2Column(:,:,k-1)+(cellDepth(:,:,k-1)*ndwet(:,:,k-1)+ &
                                     cellDepth(:,:,  k)*ndwet(:,:,  k))*f
   end do

   !IF(gcCH4%DebugIsOn) THEN
      !CALL pmaxmin('CH4: O2Column', O2Column, qmin, qmax, iXj, km,  1. )
   !END IF

   ! Grab some memory
   ! ----------------
   allocate(SZARad(i1:i2,j1:j2), STAT=status)
   VERIFY_(status)
   allocate(SZADeg(i1:i2,j1:j2), STAT=status)
   VERIFY_(status)
   allocate(sinSZA(i1:i2,j1:j2), STAT=status)
   VERIFY_(status)

   where(w_c%cosz(i1:i2,j1:j2) > 1.00)
      SZARad(i1:i2,j1:j2) = 0.00
   elsewhere
      SZARad(i1:i2,j1:j2) = ACOS(w_c%cosz(i1:i2,j1:j2))
   endwhere
   SZADeg(i1:i2,j1:j2) = SZARad(i1:i2,j1:j2)*radToDeg
   sinSZA(i1:i2,j1:j2) = sin(SZARad(i1:i2,j1:j2))

   allocate(zgrz(i1:i2,j1:j2), STAT=status)
   VERIFY_(status)

   where(SZADeg(i1:i2,j1:j2) <= 90.00)
      zgrz(i1:i2,j1:j2) = 1000.00
   elsewhere
      zgrz(i1:i2,j1:j2) = sinSZA(i1:i2,j1:j2)*(zp+r0)-r0
   endwhere

   !IF(gcCH4%DebugIsOn) THEN
      !CALL pmaxmin('CH4: zgrz', zgrz, qmin, qmax, iXj, 1,  1. )
      !CALL pmaxmin('CH4: cosz', w_c%cosz, qmin, qmax, iXj, 1,  1. )
   !END IF

   allocate(sfaca(i1:i2,j1:j2), STAT=status)
   VERIFY_(status)
   sfaca(i1:i2,j1:j2) = 0.00

   ! Chapman function calculation from ACDB 2-D model
   ! ------------------------------------------------
   do j = j1,j2
      do i = i1,i2

         if (SZADeg(i,j) < gcCH4%szaCutoff) then ! Daytime

            if (SZADeg(i,j) < 70.00) then

               sfaca(i,j) = 1.00/w_c%cosz(i,j)

            else if (zgrz(i,j) > 0.00) then

               s = b*abs(w_c%cosz(i,j))

               if (s <= 8.00) then
                  s = (d1+d2*s)/(d3+d4*s+s**2)
               else
                  s = d5/(d6+s)
               end if

               r = b*sqrt(MAPL_PI)
               sfaca(i,j) = r*s

               if (SZADeg(i,j) > 90.00) then
                  sfaca(i,j) = 2.00*r*exp((r0+zbar)*(1.00-sinSZA(i,j))/hbar)-sfaca(i,j)
               end if

            end if

         end if ! Daytime

      end do
   end do

   !IF(gcCH4%DebugIsOn) THEN
      !CALL pmaxmin('CH4: sfaca', sfaca, qmin, qmax, iXj, 1,  1. )
   !END IF

   allocate(arg(i1:i2,j1:j2), STAT=status)
   VERIFY_(status)

   ! At each layer, compute the rate constant, J [s^{-1}], if the sun is up
   ! ----------------------------------------------------------------------
   do k = 1,km

      where(SZADeg(i1:i2,j1:j2) < gcCH4%szaCutoff)
         arg(i1:i2,j1:j2) = O2Column(i1:i2,j1:j2,k)*O2xs*sfaca(i1:i2,j1:j2)
         photJ(i1:i2,j1:j2,k) = sflux*exp(-arg(i1:i2,j1:j2))*CH4xs
      end where

   end do

   !IF(gcCH4%DebugIsOn) THEN
      !CALL pmaxmin('CH4: photJ', photJ, qmin, qmax, iXj, km,  1. )
   !END IF

   deallocate(SZARad, STAT=status)
   VERIFY_(status)
   deallocate(SZADeg, STAT=status)
   VERIFY_(status)
   deallocate(sinSZA, STAT=status)
   VERIFY_(status)
   deallocate(zgrz, STAT=status)
   VERIFY_(status)
   deallocate(sfaca, STAT=status)
   VERIFY_(status)
   deallocate(arg, STAT=status)
   VERIFY_(status)
   deallocate(O2Column, ndwet, cellDepth, STAT=status)
   VERIFY_(status)

end subroutine Calc_photolysis_rates

subroutine CH4_GridCompRun1_ ( gcCH4, w_c, impChem, expChem, nymd, nhms, cdt, rc )

   implicit none

   ! INPUT/OUTPUT PARAMETERS:

   type(CH4_GridComp1), intent(inout)  :: gcCH4      ! Grid Component
   type(Chem_Bundle), intent(inout)    :: w_c        ! Chemical tracer fields

   ! INPUT PARAMETERS:

   type(ESMF_State), intent(inout)     :: impChem    ! Import State
   integer, intent(in)                 :: nymd, nhms ! time
   real,    intent(in)                 :: cdt        ! chemical timestep (secs)

   ! OUTPUT PARAMETERS:

   type(ESMF_State), intent(inout)     :: expChem    ! Export State
   integer, intent(out)                ::  rc        ! Error return code:
                                                     !  0 - all is well
                                                     !  1 -

   ! DESCRIPTION: This routine implements the CH4 Driver for GOCART.
   !
   ! !REVISION HISTORY:
   !
   !  24 Jun 2010 Nielsen:  First crack.
   !  25 Oct 2012 Nielsen:  Added photolysis.
   !
   !EOP
   !-------------------------------------------------------------------------

   character(len=*), parameter :: Iam = 'CH4::GridCompRun1_'

   !  Input fields from fvGCM
   !  -----------------------
   real, pointer, dimension(:,:)   ::  pblh     => null()
   real, pointer, dimension(:,:,:) ::  zle      => null()
   real, pointer, dimension(:,:,:) ::  T        => null()
   real, pointer, dimension(:,:,:) ::  Q        => null()
   real, pointer, dimension(:,:,:) ::  qtot     => null()

   integer :: i1, i2, im, j1, j2, jm, km, idiag, iXj
   integer :: i, j, k, kReverse, n, nbeg, nend
   integer :: nymd1, nhms1, status

   real    :: qmin, qmax

   real, allocatable, dimension(:,:,:) :: rkoh, rkcl, rko1d, rktot, dCH4ox, photJ, dCH4Phot
   real, allocatable, dimension(:,:,:) :: flux_3d, mf_current, mf_to_be_added, delp_dry

   real, pointer, dimension(:,:)   :: ptr2d => null()
   real, pointer, dimension(:,:,:) :: ptr3d => null()
   real, pointer, dimension(:,:)   :: CH4_emis => null() ! EXPORT: CH4 Emission
   !real, pointer, dimension(:,:)   :: CH4PD => null() ! EXPORT: CH4 Chemical Production
   real, pointer, dimension(:,:)   :: CH4_loss => null() ! EXPORT: CH4 Chemical Loss
   real, pointer, dimension(:,:)   :: CH4_surf => null() ! EXPORT: CH4 Surface Concentration
   real, pointer, dimension(:,:)   :: CH4_column => null() ! EXPORT: CH4 Column Burden
   real, pointer, dimension(:,:)   :: CH4_phot => null() ! EXPORT: CH4 Photolytic Loss
   !real, pointer, dimension(:,:,:) :: CH4QP => null() ! EXPORT: H2O tendency from CH4 photolysis

   call MAPL_GetPointer(expChem, CH4_emis,   trim(gcCH4%iname)//func_sep//'EM', __RC__)
   call MAPL_GetPointer(expChem, CH4_loss,   trim(gcCH4%iname)//func_sep//'LS', __RC__)
   call MAPL_GetPointer(expChem, CH4_surf,   trim(gcCH4%iname)//func_sep//'SC', __RC__)
   call MAPL_GetPointer(expChem, CH4_column, trim(gcCH4%iname)//func_sep//'CL', __RC__)
   call MAPL_GetPointer(expChem, CH4_phot,   trim(gcCH4%iname)//func_sep//'JL', __RC__)
   !call MAPL_GetPointer( expChem, CH4PD,  'CH4PD'//trim(gcCH4%iname), __RC__)
   !call MAPL_GetPointer ( expChem, CH4QP,  'CH4QP'//trim(gcCH4%iname), __RC__)

   !  Initialize local variables
   !  --------------------------
   rc = 0
   i1 = w_c%grid%i1
   i2 = w_c%grid%i2
   im = w_c%grid%im

   j1 = w_c%grid%j1
   j2 = w_c%grid%j2
   jm = w_c%grid%jm
   km = w_c%grid%km

   iXj = ( i2 - i1 + 1 ) * ( j2 - j1 + 1 )

   nbeg  = w_c%reg%i_CH4
   nend  = w_c%reg%j_CH4

   !  It requires 1 bin
   !  -----------------
   if (nbeg /= nend) then
      if (MAPL_AM_I_ROOT()) then
         print *, trim(gcCH4%iname)//": Must have only 1 bin at the single instance level"
      end if
      status = 1
      VERIFY_(status)
   end if

   !  Get imports
   !  -----------
   call MAPL_GetPointer(impChem, PBLH,     'ZPBL',    __RC__)
   call MAPL_GetPointer(impChem, ZLE,      'ZLE',     __RC__)
   call MAPL_GetPointer(impChem, T,        'T',       __RC__)
   call MAPL_GetPointer(impChem, Q,        'Q',       __RC__)
   call MAPL_GetPointer(impChem, qtot,     'QTOT',    __RC__)

   if (gcCH4%DebugIsOn) then
      !call pmaxmin('CH4:AREA', cellArea, qmin, qmax, iXj,  1,   1. )
      call pmaxmin('CH4:ZPBL',     pblh, qmin, qmax, iXj, km+1, 1. )
      call pmaxmin('CH4:ZLE',       zle, qmin, qmax, iXj, km+1, 1. )
      call pmaxmin('CH4:T',           T, qmin, qmax, iXj, km,   1. )
      call pmaxmin('CH4:Q',           q, qmin, qmax, iXj, km,   1. )
      call pmaxmin('CH4:QTOT',     qtot, qmin, qmax, iXj, km,   1. ) ! Sourish
      !call pmaxmin('CH4:RHOWET', rhowet, qmin, qmax, iXj, km,   1. )
      !call pmaxmin('CH4:PLE',       ple, qmin, qmax, iXj, km+1, 1. )
   end if

   !  Update CH4 emissions and OH, Cl, and O(1D) number densities
   !  (in molec cm^-3) once each day
   !  ------------------------------------------------------------
   ! Read sources for this instance
   do i = 1, gcCH4%n_sources
      call MAPL_GetPointer(impChem, ptr2d, trim(gcCH4%source_categs(i)), __RC__)
      gcCH4%emis(i,:,:) = ptr2d(:,:)
   end do

   call MAPL_GetPointer(impChem,ptr3d,'CH4_oh', __RC__)
   gcCH4%OHnd=ptr3d

   call MAPL_GetPointer(impChem,ptr3d,'CH4_cl', __RC__)
   gcCH4%Clnd=ptr3d

   call MAPL_GetPointer(impChem,ptr3d,'CH4_o1d', __RC__)
   gcCH4%O1Dnd=ptr3d

   !  Convert number densities from molec cm^-3 to molec m^-3
   !  -------------------------------------------------------
   gcCH4%OHnd(i1:i2,j1:j2,1:km) = gcCH4%OHnd(i1:i2,j1:j2,1:km)*1.00E+06
   gcCH4%Clnd(i1:i2,j1:j2,1:km)  =  gcCH4%Clnd(i1:i2,j1:j2,1:km)*1.00E+06
   gcCH4%O1Dnd(i1:i2,j1:j2,1:km) = gcCH4%O1Dnd(i1:i2,j1:j2,1:km)*1.00E+06

   !  Allocate temporary workspace
   !  ----------------------------
   allocate( &
      rkoh(i1:i2,j1:j2,km), rkcl(i1:i2,j1:j2,km), rko1d(i1:i2,j1:j2,km), rktot(i1:i2,j1:j2,km), dCH4ox(i1:i2,j1:j2,km), &
      STAT=status)
   VERIFY_(status)

!   if (gcCH4%DebugIsOn) then
!      call pmaxmin('CH4: OH Conc',        gcCH4%OHnd, qmin, qmax, iXj, 1,  1. )
!      call pmaxmin('CH4: Cl Conc',        gcCH4%Clnd, qmin, qmax, iXj, 1,  1. )
!      call pmaxmin('CH4: O(1D) Conc',    gcCH4%O1Dnd, qmin, qmax, iXj, 1,  1. )
!      call pmaxmin('CH4: Cell Vol',       cellVolume, qmin, qmax, iXj, km, 1. )
!      call pmaxmin('CH4: Cell Depth',      cellDepth, qmin, qmax, iXj, km, 1. )
!      call pmaxmin('CH4: Wet-air ND',          ndwet, qmin, qmax, iXj, km, 1. )
!   end if

   !  Compute and add surface emissions
   !  ---------------------------------
   allocate(mf_to_be_added(i1:i2,j1:j2,1:km), mf_current(i1:i2,j1:j2,1:km), STAT=rc)
   allocate(flux_3d(i1:i2,j1:j2,1:km), stat=rc)  ! 3D emissions
   call Calc_emissions(gcCH4, w_c, impChem, nhms, cdt, flux_3d, __RC__)

   mf_to_be_added = cdt * flux_3d*MAPL_AIRMW / mwtCH4 / (w_c%delp(i1:i2,j1:j2,1:km)/MAPL_GRAV)
   ! Debug to check if any of the mole fractions will turn negative after adding fluxes
   !mf_current = w_c%qa(nbeg)%data3d(i1:i2,j1:j2,1:km)
   !if (any(mf_current+mf_to_be_added .le. 0.0)) then
      !write(*,'("    For tag ", a, ", cells where adding emissions makes MF negative = ", i0)') &
         !trim(gcCH4%iname), count(mf_current+mf_to_be_added .le. 0.0)
   !end if
   ! End debug

   ! Update mole fraction
   w_c%qa(nbeg)%data3d(i1:i2,j1:j2,1:km) = w_c%qa(nbeg)%data3d(i1:i2,j1:j2,1:km) + mf_to_be_added

   ! Update Surface flux diagnostic for this bin
   ! -------------------------------------------
   if (associated(CH4_emis)) then
      CH4_emis(i1:i2,j1:j2) = sum(flux_3d, dim=3)
   end if

   deallocate(mf_to_be_added, mf_current, flux_3d, STAT=rc)

   if (gcCH4%chemistry) then
      !  Loss rate [m^3 s^-1] for OH + CH4 => CO + products
      rkoh(i1:i2,j1:j2,1:km)  = 2.45E-12*1.00E-06*exp(-1775./T(i1:i2,j1:j2,1:km))
      !  Loss rate [m^3 s^-1] for Cl + CH4 => CO + products
      rkcl(i1:i2,j1:j2,1:km)  = 7.10E-12*1.00E-06*exp(-1270./T(i1:i2,j1:j2,1:km))
      !  Loss rate [m^3 s^-1] for O(1D) + CH4 => CO + products
      rko1d(i1:i2,j1:j2,1:km) = 1.75E-10*1.00E-06
      !  Combine the three loss mechanisms
      rktot = rkoh*gcCH4%OHnd(i1:i2,j1:j2,1:km) + rkcl*gcCH4%Clnd(i1:i2,j1:j2,1:km) + rko1d*gcCH4%O1Dnd(i1:i2,j1:j2,1:km)
      !  Change in CH4 mole fraction due to oxidation
      dCH4ox = -cdt * rktot * w_c%qa(nbeg)%data3d(i1:i2,j1:j2,1:km) ! fwd Euler (1st order)
      w_c%qa(nbeg)%data3d(i1:i2,j1:j2,1:km) = w_c%qa(nbeg)%data3d(i1:i2,j1:j2,1:km) + dCH4ox
   else
      dCH4ox = 0.0
   end if

   if (gcCH4%photolysis) then
      !  Calculate photolytic loss rates, J [s^-1], for CH4 + hv => 2H2O + CO
      allocate(photJ(i1:i2,j1:j2,1:km), dCH4Phot(i1:i2,j1:j2,1:km), STAT=status)
      VERIFY_(STATUS)
      call Calc_photolysis_rates(gcCH4, w_c, impChem, photJ, __RC__)
      !  Michael Long says (email Dec 5, 2022) that photJ is 'per second', i.e., it's a decay rate
      dCH4Phot = -cdt * photJ(i1:i2,j1:j2,1:km) * w_c%qa(nbeg)%data3d(i1:i2,j1:j2,1:km) ! fwd Euler (1st order)
      !  Change the mole fraction accordingly
      w_c%qa(nbeg)%data3d(i1:i2,j1:j2,1:km) = w_c%qa(nbeg)%data3d(i1:i2,j1:j2,1:km) + dCH4Phot
      deallocate(photJ, dCH4Phot, STAT=status)
      VERIFY_(status)
   else
      dCH4Phot = 0.0
   end if

   !  Summarize the loss as an export, vertically integrated
   if (associated(CH4_loss)) then
      CH4_loss(i1:i2,j1:j2) = 0.
      do k = 1, km
         CH4_loss(i1:i2,j1:j2) = CH4_loss(i1:i2,j1:j2) + &
                                 dCH4ox(i1:i2,j1:j2,k) * (mwtCH4/MAPL_AIRMW) * w_c%delp(i1:i2,j1:j2,k)/MAPL_GRAV ! kg/m^2 over this time step
      end do
      CH4_loss(i1:i2,j1:j2) = CH4_loss(i1:i2,j1:j2)/cdt ! convert to kg/m^2/s, then history can export a time-averaged value
   end if
   ! Also summarize photolytic loss
   if (associated(CH4_phot)) THEN
      CH4_phot(i1:i2,j1:j2) = 0.
      do k = 1, km
         CH4_phot(i1:i2,j1:j2) = CH4_phot(i1:i2,j1:j2) + &
                                 dCH4Phot(i1:i2,j1:j2,k) * (mwtCH4/MAPL_AIRMW) * w_c%delp(i1:i2,j1:j2,k)/MAPL_GRAV ! kg/m^2 over this time step
      end do
      CH4_phot(i1:i2,j1:j2) = CH4_phot(i1:i2,j1:j2)/cdt ! convert to kg/m^2/s, then history can export a time-averaged value
   end if

   !  If both feedback switches are on, increment water vapor by
   !  adding two molecules of H2O for each CH4 molecule lost
   !  ----------------------------------------------------------
!   if (associated(CH4_qprod)) CH4_qprod(i1:i2,j1:j2,1:km) = 0.

   if (gcCH4%CH4FeedBack .and. gcCH4%H2OFeedBack) then
      Q(i1:i2,j1:j2,1:km) = Q(i1:i2,j1:j2,1:km) + 2.00*cdt*dCH4Phot(i1:i2,j1:j2,1:km)*MAPL_H2OMW/MAPL_AIRMW

   !     Water vapor tendency [kg kg^-1 s^-1]
   !     ------------------------------------
!      if (associated(CH4_qprod)) then
!         CH4_qprod(i1:i2,j1:j2,1:km) = 2.00*dCH4Phot(i1:i2,j1:j2,1:km)*MAPL_H2OMW/MAPL_AIRMW
!      end if
   end if

   !  Surface concentration (dry air mole fraction)
   !  ----------------------------
   if (associated(CH4_surf)) then
      CH4_surf(i1:i2,j1:j2) = w_c%qa(nbeg)%data3d(i1:i2,j1:j2,km) / (1. - qtot(i1:i2,j1:j2,km))
   end if

   !  Column burden (dry air mole fraction)
   !  ----------------------
   if (associated(CH4_column)) then
      allocate(delp_dry(i1:i2,j1:j2,1:km))
      do k = 1, km
         delp_dry(:,:,k) = w_c%delp(i1:i2,j1:j2,k) * (1. - qtot(i1:i2,j1:j2,k))
      end do
      CH4_column(i1:i2,j1:j2) = sum(w_c%qa(nbeg)%data3d(i1:i2,j1:j2,1:km)*w_c%delp(i1:i2,j1:j2,1:km), dim=3)/sum(delp_dry, dim=3)
      deallocate(delp_dry)
   end if

!   IF(gcCH4%DebugIsOn) THEN
!      IF(ASSOCIATED(CH4_emis))    CALL pmaxmin(    'CH4: emis', CH4_emis,    qmin, qmax, iXj,  1, 1. )
!      IF(ASSOCIATED(CH4_loss))    CALL pmaxmin(    'CH4: loss', CH4_loss,    qmin, qmax, iXj,  1, 1. )
!      IF(ASSOCIATED(CH4_prod))    CALL pmaxmin(    'CH4: prod', CH4_prod,    qmin, qmax, iXj,  1, 1. )
!      IF(ASSOCIATED(CH4_column))  CALL pmaxmin(  'CH4: column', CH4_column,  qmin, qmax, iXj,  1, 1. )
!      IF(ASSOCIATED(CH4_surface)) CALL pmaxmin( 'CH4: surface', CH4_surface, qmin, qmax, iXj,  1, 1. )
!      IF(ASSOCIATED(CH4_qprod))   CALL pmaxmin('CH4: qprod',    CH4_qprod,   qmin, qmax, iXj, km, 1. )
!      IF(ASSOCIATED(CH4_phot))    CALL pmaxmin('CH4: dch4phot', dCH4phot,    qmin, qmax, iXj, km, 1. )
!      IF(ASSOCIATED(CH4_dry))     CALL pmaxmin(     'CH4: dry', CH4_dry,     qmin, qmax, iXj, km, 1. )
!   END IF


   !  Housekeeping
   !  ------------
   deallocate(rkoh, rkcl, rko1d, rktot, dCH4ox, STAT=status)
   VERIFY_(status)

end subroutine CH4_GridCompRun1_

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  CH4_GridCompFinalize --- The Chem Driver
!
! !INTERFACE:
!

subroutine CH4_GridCompFinalize1_ ( gcCH4, w_c, impChem, expChem, nymd, nhms, cdt, rc )

   implicit none

   ! INPUT/OUTPUT PARAMETERS:

   type(CH4_GridComp1), intent(inout) :: gcCH4   ! Grid Component

   ! INPUT PARAMETERS:

   type(Chem_Bundle), intent(in)  :: w_c      ! Chemical tracer fields
   integer, intent(in) :: nymd, nhms          ! time
   real,    intent(in) :: cdt             ! chemical timestep (secs)


   ! OUTPUT PARAMETERS:

   type(ESMF_State), intent(inout) :: impChem   ! Import State
   type(ESMF_State), intent(inout) :: expChem   ! Import State
   integer, intent(out) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 -

   ! DESCRIPTION: This routine finalizes this Grid Component.
   !
   ! REVISION HISTORY:
   !
   !  18Sep2003 da Silva  First crack.
   !
   !EOP
   !-------------------------------------------------------------------------

   character(len=*), parameter :: Iam = 'CH4::GridCompFinalize1_'
   rc = 0

   deallocate(gcCH4%regionMask, gcCH4%OHnd, gcCH4%Clnd, gcCH4%O1Dnd, STAT=rc)
   VERIFY_(rc)

   deallocate(gcCH4%emis, gcCH4%source_categs, gcCH4%diurnal_fire, gcCH4%pbl_inject, stat=rc)
   VERIFY_(rc)

end subroutine CH4_GridCompFinalize1_

end module CH4_GridCompMod

!-----------------------------------------------------------------------

!                     Single Instance Wrapper

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  CH4_SingleInstance_ --- Runs single instance of method
!
! !INTERFACE:
!
subroutine CH4_SingleInstance_( Method_, instance, gcCH4, w_c, impChem, expChem, nymd, nhms, cdt, rc )

   ! USES:

   Use CH4_GridCompMod
   Use ESMF
   Use MAPL
   Use Chem_Mod

   implicit none

   ! INPUT PARAMETERS:

   !  Input "function pointer"
   !  -----------------------
   interface
      subroutine Method_ (gc, w, imp, exp, ymd, hms, dt, rcode )
         Use CH4_GridCompMod
         Use ESMF
         Use MAPL
         Use Chem_Mod
         type(CH4_GridComp1), intent(inout)  :: gc
         type(Chem_Bundle),   intent(in)     :: w
         type(ESMF_State),    intent(inout)  :: imp
         type(ESMF_State),    intent(inout)  :: exp
         integer,             intent(in)     :: ymd, hms
         real,                intent(in)     :: dt
         integer,             intent(out)    :: rcode
      end subroutine Method_
   end interface

   integer, intent(in)           :: instance   ! instance number

   type(Chem_Bundle), intent(inout) :: w_c     ! Chemical tracer fields
   integer, intent(in) :: nymd, nhms           ! time
   real,    intent(in) :: cdt              ! chemical timestep (secs)


   ! OUTPUT PARAMETERS:

   type(CH4_GridComp1), intent(inout) :: gcCH4    ! Grid Component
   type(ESMF_State), intent(inout)  :: impChem  ! Import State
   type(ESMF_State), intent(inout)  :: expChem  ! Export State
   integer, intent(out) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 -

   ! !DESCRIPTION: Finalizes the CH4 Grid Component. Multiple instance
   !               version.
   !
   ! !REVISION HISTORY:
   !
   !  27Feb2008  da Silva  Introduced multiple instances
   !
   !EOP
   !-------------------------------------------------------------------------

   integer :: n_CH4, i_CH4, j_CH4

   ! Save overall CH4 indices
   ! -----------------------
   n_CH4 = w_c%reg%n_CH4
   i_CH4 = w_c%reg%i_CH4
   j_CH4 = w_c%reg%j_CH4

   ! Customize indices for this particular instance
   ! ----------------------------------------------
   w_c%reg%n_CH4 = 1
   w_c%reg%i_CH4 = i_CH4 + instance - 1
   w_c%reg%j_CH4 = i_CH4 + instance - 1

   ! Execute the instance method
   ! ---------------------------
   call Method_ ( gcCH4, w_c, impChem, expChem, &
                 nymd, nhms, cdt, rc )

   ! Restore the overall CH4 indices
   ! ------------------------------
   w_c%reg%n_CH4 = n_CH4
   w_c%reg%i_CH4 = i_CH4
   w_c%reg%j_CH4 = j_CH4

end subroutine CH4_SingleInstance_

!-----------------------------------------------------------------------
