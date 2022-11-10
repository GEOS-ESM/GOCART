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

TYPE CH4_GridComp1

   character(LEN=ESMF_MAXSTR) :: name        ! generic name of the package
   character(LEN=ESMF_MAXSTR) :: iname       ! instance name
   character(LEN=ESMF_MAXSTR) :: rcfilen     ! resource file name
   character(LEN=ESMF_MAXSTR) :: regionsString   ! Comma-delimited string of regions
   character(LEN=ESMF_MAXSTR) :: CH4Source   ! Source name on emission file (CH4_ANIMLS, for example)
   character(len=256), allocatable :: categ_names(:) ! each instance can be a sum of multiple categories

   integer :: instance                 ! Instantiation number
   integer :: nymd_oh
   integer :: nymd_ch4
   integer :: BCnymd                   ! Date of last emissions/prodction read
   integer :: n_categ                  ! number of CH4 categories or sources that contri

   REAL, POINTER :: regionMask(:,:)    ! regional mask
   REAL, POINTER ::  CH4(:,:,:)    ! CH4 mixing ratio mol/mol
   REAL, POINTER :: OHnd(:,:,:)    ! OH number density (cm^{-3})
   REAL, POINTER ::     Clnd(:,:,:)    ! Cl number density (cm^{-3})
   REAL, POINTER ::    O1Dnd(:,:,:)    ! O(1D) number density (cm^{-3})

!   REAL, POINTER :: eCH4_wetland(:,:)    ! kgCH4/m2/s, Earth surface
!   REAL, POINTER :: eCH4_industrial(:,:) ! kgCH4/m2/s, Earth surface
!   REAL, POINTER :: eCH4_extract(:,:)    ! kgCH4/m2/s, Earth surface
!   REAL, POINTER :: eCH4_transport(:,:)  ! kgCH4/m2/s, Earth surface
!   REAL, POINTER :: eCH4_agwaste(:,:)    ! kgCH4/m2/s, PBL (before diurnal)
!   REAL, POINTER :: eCH4_onat(:,:)       ! kgCH4/m2/s, PBL
!   REAL, POINTER :: eCH4_fire(:,:)       ! kgCH4/m2/s, Earth surface
   real, allocatable :: eCH4(:,:)         ! trying out one category per instance (kg CH4/m2/s)
   real, allocatable :: eCH4_d(:,:)       ! needed if we want diurnal cycle

   LOGICAL :: DebugIsOn     ! Run-time debug switch
   LOGICAL :: CH4FeedBack   ! Permit increments to CH4 from CH4 + hv => 2H2O + CO
   LOGICAL :: H2OFeedBack   ! Permit increments to   Q from CH4 + hv => 2H2O + CO
   logical :: chemistry     ! Should we do chemistry, such as oxidation and photolysis? Default true

   REAL :: szaCutoff        ! Largest solar zenith angle (degrees) allowed as daytime

END TYPE CH4_GridComp1

TYPE CH4_GridComp
   INTEGER                            ::  n_inst   ! number of instances
   TYPE(CH4_GridComp1), allocatable   ::  gcs(:)   ! instances
END TYPE CH4_GridComp

real, parameter :: radToDeg = 180./MAPL_PI
real, parameter :: mwtCH4   = 16.0422

CONTAINS

subroutine CH4_GridCompSetServices(  gc, chemReg, rc)
   type(ESMF_GridComp), intent(INOUT) :: GC
   type(Chem_Registry), intent(INOUT) :: chemReg
   integer,             intent(OUT  ) :: rc

   CHARACTER(LEN=255) :: rcbasen = 'CH4_GridComp'
   CHARACTER(LEN=255) :: name

   integer            :: status
   character(len=ESMF_MAXSTR) :: Iam

   integer :: ier,n,i

   type(ESMF_Config) :: cfg

   Iam = "CH4_GridCompSetServices"

   !  Load resource file
   !  ------------------
   cfg = ESMF_ConfigCreate(rc=status)
   VERIFY_(STATUS)
   call ESMF_ConfigLoadFile(cfg,TRIM(rcbasen)//'.rc',rc=status) ! reading CH4_GridComp.rc
   VERIFY_(STATUS)

   !  Parse resource file
   !  -------------------
   n = ESMF_ConfigGetLen(cfg,label='CH4_instances:',rc=status)
   VERIFY_(STATUS)

   !  We cannot have fewer instances than the number of
   !   CH4 bins in the registry (it is OK to have less, though)
   !  --------------------------------------------------------
   if ( n .LT. chemReg%n_CH4 ) then ! nbins_CH4 > number of instances in CH4_GridComp.rc
        rc = 35
        return
   else if ( n .GT. chemReg%n_CH4 ) then
        if (MAPL_AM_I_ROOT()) &
        PRINT *, TRIM(Iam)//": Bins = ",chemReg%n_CH4," of ",n," expected."
   end if
   n = min(n,chemReg%n_CH4 )

   !write(*,'(a, " :: on line ", i0, ", rc = ", i0)') trim(Iam), __LINE__, rc

   !  Record name of each instance
   !  ----------------------------
   call ESMF_ConfigFindLabel(cfg,'CH4_instances:',rc=status) ! these could be called 'Basu', 'Weir', etc.
   VERIFY_(STATUS)
   do i = 1, n
      call ESMF_ConfigGetAttribute(cfg,name,rc=status)
      VERIFY_(STATUS)
                                            ! resource file name
      !IF(TRIM(name) == "full" .or. trim(name) == "total") THEN
       !name = " "              ! blank instance name for full (1)
      !ELSE
       !name = TRIM(name)       ! instance name for others
      !END IF
      call CH4_GridCompSetServices1_(gc,chemReg,name,rc=status)
      VERIFY_(STATUS)
   end do

   call MAPL_AddImportSpec(GC,           &
        SHORT_NAME = 'CH4_regionMask',   &
        LONG_NAME  = 'source species'  , &
        UNITS      = '1',                &
        DIMS       = MAPL_DimsHorzOnly,  &
        VLOCATION  = MAPL_VLocationNone, &
        RESTART    = MAPL_RestartSkip,   &
        RC         = STATUS)
   VERIFY_(STATUS)
! Sourish Basu
   call MAPL_AddExportSpec(GC,  &
      SHORT_NAME         = 'CH4',  &
      LONG_NAME          = 'CH4 total mole fraction',  &
      UNITS              = 'mol mol-1', &
      DIMS               = MAPL_DimsHorzVert,    &
      VLOCATION          = MAPL_VLocationCenter,    &
      RC=STATUS  )
   _VERIFY(STATUS)

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
   INTEGER, INTENT(IN) :: nymd, nhms           ! time
   REAL,    INTENT(IN) :: cdt              ! chemical timestep (secs)


   ! OUTPUT PARAMETERS:

   TYPE(CH4_GridComp), intent(inout) :: gcCH4   ! Grid Component
   TYPE(ESMF_State), intent(inout)  :: impChem  ! Import State
   TYPE(ESMF_State), intent(inout)  :: expChem  ! Export State
   INTEGER, intent(out) ::  rc                  ! Error return code:
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

   CHARACTER(LEN=*), PARAMETER :: Iam = 'CH4_GridCompInitialize'
   CHARACTER(LEN=ESMF_MAXSTR) :: rcbasen = 'CH4_GridComp'
   CHARACTER(LEN=ESMF_MAXSTR) :: name

   INTEGER :: i, n, status
   !   REAL :: c1,c2,c3,c4

   !  Load resource file
   !  ------------------
   CALL I90_loadf ( TRIM(rcbasen)//'.rc', status ) ! read from CH4_GridComp.rc (must be in current folder?)
   VERIFY_(status)

   !  Parse resource file
   !  -------------------
   CALL I90_label ( 'CH4_instances:', status ) ! e.g., 'wetland industrial extract transport agwaste onat fire total'
   VERIFY_(status)

   !  First determine how many instances we have
   !  ------------------------------------------
   n = 0
   status = 0

   DO WHILE ( status == 0 )
      CALL I90_gtoken( name, status )
      IF(status == 0) n = n + 1
   END DO
   IF ( n == 0 ) THEN
      status = 1
      VERIFY_(status)
   END IF

   !  We cannot have fewer instances than the number of
   !   CH4 bins in the registry (it is OK to have less, though)
   !  --------------------------------------------------------
   ! SB :: Don't understand this yet, will get to it later
   IF ( n < w_c%reg%n_CH4 ) THEN
      status = 1
      VERIFY_(status)
   ELSE IF ( n >= w_c%reg%n_CH4 ) THEN
      IF(MAPL_AM_I_ROOT()) PRINT *, TRIM(Iam)//": Bins = ",w_c%reg%n_CH4," of ",n," expected."
   END IF
   n = min(n,w_c%reg%n_CH4 )
   gcCH4%n_inst = n

   !  Next allocate necessary memory
   !  ------------------------------
   ALLOCATE ( gcCH4%gcs(n), STAT=status ) ! each instance gets a gcCH4%gcs
   VERIFY_(status)

   !  Record name of each instance
   !  ----------------------------
   CALL I90_label ( 'CH4_instances:', status )
   VERIFY_(status)
   DO i = 1, n
      CALL I90_gtoken( name, status )
      VERIFY_(status)
      gcCH4%gcs(i)%rcfilen = trim(rcbasen)//'.rc' ! Experiment to read all rc keys from one file
      gcCH4%gcs(i)%instance = i              ! instance number
      !IF(TRIM(name) == "full" .or. trim(name) == "total") THEN
         !gcCH4%gcs(i)%iname = " "              ! blank instance name for full (1) ! why?
      !ELSE
         !gcCH4%gcs(i)%iname = TRIM(name)       ! instance name for others
      !END IF
      gcCH4%gcs(i)%iname = TRIM(name)
   END DO

   !  Next initialize each instance
   !  -----------------------------
   DO i = 1, gcCH4%n_inst
      IF(MAPL_AM_I_ROOT()) THEN
         PRINT *," "
         PRINT *,TRIM(Iam)//": Initializing instance ",TRIM(gcCH4%gcs(i)%iname)," [",gcCH4%gcs(i)%instance,"]"
      END IF
      CALL CH4_SingleInstance_ ( CH4_GridCompInitialize1_, i, &
                                gcCH4%gcs(i), w_c, impChem, expChem,  &
                                nymd, nhms, cdt, status )
      VERIFY_(status)
   END DO

   !!  Get Henrys Law cts for the parameterized convective wet removal
   !!  -----------------------------------------------------------
   !   CALL get_HenrysLawCts('CH4',c1,c2,c3,c4)
   !   w_c%reg%Hcts(1,w_c%reg%i_CH4 : w_c%reg%j_CH4)=c1
   !   w_c%reg%Hcts(2,w_c%reg%i_CH4 : w_c%reg%j_CH4)=c2
   !   w_c%reg%Hcts(3,w_c%reg%i_CH4 : w_c%reg%j_CH4)=c3
   !   w_c%reg%Hcts(4,w_c%reg%i_CH4 : w_c%reg%j_CH4)=c4

   !  All done
   !  --------
   CALL I90_FullRelease( status )
   VERIFY_(status)

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

   TYPE(Chem_Bundle), intent(in) :: w_c        ! Chemical tracer fields
   INTEGER, INTENT(IN) :: nymd, nhms           ! time
   REAL,    INTENT(IN) :: cdt              ! chemical timestep (secs)


   !OUTPUT PARAMETERS:

   TYPE(CH4_GridComp), INTENT(INOUT) :: gcCH4   ! Grid Component
   TYPE(ESMF_State), INTENT(INOUT)  :: impChem  ! Import State
   TYPE(ESMF_State), INTENT(INOUT)  :: expChem  ! Export State
   INTEGER, INTENT(OUT) ::  rc                  ! Error return code:
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

   CHARACTER(LEN=ESMF_MAXSTR), PARAMETER :: Iam = 'CH4_GridCompRun'
   INTEGER :: i, status

   DO i = 1, gcCH4%n_inst
      CALL CH4_SingleInstance_ ( CH4_GridCompRun1_, i, &
                                gcCH4%gcs(i), w_c, impChem, expChem, &
                                nymd, nhms, cdt, status )
      VERIFY_(status)
   END DO

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

   !USES:

   IMPLICIT NONE

   ! INPUT PARAMETERS:

   TYPE(Chem_Bundle), intent(in) :: w_c        ! Chemical tracer fields
   INTEGER, INTENT(IN) :: nymd, nhms           ! time
   REAL,    INTENT(IN) :: cdt              ! chemical timestep (secs)


   ! OUTPUT PARAMETERS:

   TYPE(CH4_GridComp), INTENT(INOUT) :: gcCH4   ! Grid Component
   TYPE(ESMF_State), INTENT(INOUT)  :: impChem  ! Import State
   TYPE(ESMF_State), INTENT(INOUT)  :: expChem  ! Export State
   INTEGER, INTENT(OUT) ::  rc                  ! Error return code:
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

   CHARACTER(LEN=ESMF_MAXSTR), PARAMETER :: Iam = 'CH4_GridCompFinalize'
   INTEGER :: i, status

   DO i = 1, gcCH4%n_inst
      CALL CH4_SingleInstance_ ( CH4_GridCompFinalize1_, i, &
                                gcCH4%gcs(i), w_c, impChem, expChem, &
                                nymd, nhms, cdt, status )
      VERIFY_(status)
   END DO

   DEALLOCATE ( gcCH4%gcs, stat=status )
   gcCH4%n_inst = -1

end subroutine CH4_GridCompFinalize

subroutine CH4_GridCompSetServices1_(  gc, chemReg, iname, rc)
   type(ESMF_GridComp), intent(INOUT) :: GC
   type(Chem_Registry), intent(INOUT) :: chemReg
   character(len=*),    intent(IN   ) :: iname
   integer,             intent(OUT  ) :: rc

   integer :: Status
   character(len=ESMF_MAXSTR) :: Iam

   Iam ="CH4_GridCOmpSetServices1_"

   call MAPL_AddImportSpec(GC, &
      SHORT_NAME = 'CH4_'//trim(iname), &
      LONG_NAME  = 'source species'  , &
      UNITS      = 'kg CH4 m-2 s-1',     &
      DIMS       = MAPL_DimsHorzOnly, &
      VLOCATION  = MAPL_VLocationNone, &
      RESTART    = MAPL_RestartSkip,     &
      RC         = STATUS)
!   call MAPL_AddImportSpec(GC, &
!      SHORT_NAME = 'CH4_industr'//iname, &
!      LONG_NAME  = 'source species'  , &
!      UNITS      = 'kg CH4 m-2 s-1',     &
!      DIMS       = MAPL_DimsHorzOnly, &
!      VLOCATION  = MAPL_VLocationNone, &
!      RESTART    = MAPL_RestartSkip,     &
!      RC         = STATUS)
!   call MAPL_AddImportSpec(GC,             &
!      SHORT_NAME = 'CH4_agwaste'//iname, &
!      LONG_NAME  = 'source species'  ,   &
!      UNITS      = 'kg CH4 m-2 s-1',     &
!      DIMS       = MAPL_DimsHorzOnly,    &
!      VLOCATION  = MAPL_VLocationNone,   &
!      RESTART    = MAPL_RestartSkip,     &
!      RC         = STATUS)
!   call MAPL_AddImportSpec(GC,             &
!      SHORT_NAME = 'CH4_bioburn'//iname, &
!      LONG_NAME  = 'source species'  ,   &
!      UNITS      = 'kg CH4 m-2 s-1',     &
!      DIMS       = MAPL_DimsHorzOnly,    &
!      VLOCATION  = MAPL_VLocationNone,   &
!      RESTART    = MAPL_RestartSkip,     &
!      RC         = STATUS)
!   call MAPL_AddImportSpec(GC,             &
!      SHORT_NAME = 'CH4_biofuel'//iname, &
!      LONG_NAME  = 'source species'  ,   &
!      UNITS      = 'kg CH4 m-2 s-1',     &
!      DIMS       = MAPL_DimsHorzOnly,    &
!      VLOCATION  = MAPL_VLocationNone,   &
!      RESTART    = MAPL_RestartSkip,     &
!      RC         = STATUS)
!   call MAPL_AddImportSpec(GC,             &
!      SHORT_NAME = 'CH4_minnatl'//iname, &
!      LONG_NAME  = 'source species'  ,   &
!      UNITS      = 'kg CH4 m-2 s-1',     &
!      DIMS       = MAPL_DimsHorzOnly,    &
!      VLOCATION  = MAPL_VLocationNone,   &
!      RESTART    = MAPL_RestartSkip,     &
!      RC         = STATUS)
!   call MAPL_AddImportSpec(GC,             &
!      SHORT_NAME = 'CH4_wetland'//iname, &
!      LONG_NAME  = 'source species',     &
!      UNITS      = 'kg CH4 m-2 s-1',     &
!      DIMS       = MAPL_DimsHorzOnly,    &
!      VLOCATION  = MAPL_VLocationNone,   &
!      RESTART    = MAPL_RestartSkip,     &
!      RC         = STATUS)
   call MAPL_AddImportSpec(GC, &
      SHORT_NAME = 'CH4_oh', &
      LONG_NAME  = 'source species'  , &
      UNITS      = 'molecules cm-3',     &
      DIMS       = MAPL_DimsHorzVert,    &
      VLOCATION  = MAPL_VLocationCenter, &
      RC         = STATUS)
   call MAPL_AddImportSpec(GC,             &
      SHORT_NAME = 'CH4_cl',      &
      LONG_NAME  = 'source species',     &
      UNITS      = 'molecules cm-3',     &
      DIMS       = MAPL_DimsHorzVert,    &
      VLOCATION  = MAPL_VLocationCenter, &
      RC         = STATUS)
   call MAPL_AddImportSpec(GC,             &
      SHORT_NAME = 'CH4_o1d',     &
      LONG_NAME  = 'source species',     &
      UNITS      = 'molecules cm-3',     &
      DIMS       = MAPL_DimsHorzVert, &
      VLOCATION  = MAPL_VLocationCenter, &
      RC         = STATUS)
!   call MAPL_AddImportSpec(GC, &
!      SHORT_NAME = 'CH4_oh'//iname, &
!      LONG_NAME  = 'source species'  , &
!      UNITS      = 'molecules cm-3',     &
!      DIMS       = MAPL_DimsHorzVert,    &
!      VLOCATION  = MAPL_VLocationCenter, &
!      RC         = STATUS)
!   call MAPL_AddImportSpec(GC,             &
!      SHORT_NAME = 'CH4_cl'//iname,      &
!      LONG_NAME  = 'source species',     &
!      UNITS      = 'molecules cm-3',     &
!      DIMS       = MAPL_DimsHorzVert,    &
!      VLOCATION  = MAPL_VLocationCenter, &
!      RC         = STATUS)
!   call MAPL_AddImportSpec(GC,             &
!      SHORT_NAME = 'CH4_o1d'//iname,     &
!      LONG_NAME  = 'source species',     &
!      UNITS      = 'molecules cm-3',     &
!      DIMS       = MAPL_DimsHorzVert, &
!      VLOCATION  = MAPL_VLocationCenter, &
!      RC         = STATUS)
   VERIFY_(STATUS)

   call MAPL_AddExportSpec(GC,  &
      SHORT_NAME         = 'CH4EM'//trim(iname),  &
      LONG_NAME          = 'CH4 Emission '//trim(iname),  &
      UNITS              = 'kg m-2 s-1', &
      DIMS               = MAPL_DimsHorzOnly,    &
      VLOCATION          = MAPL_VLocationNone,    &
      RC=STATUS  )
   _VERIFY(STATUS)

   call MAPL_AddExportSpec(GC,  &
      SHORT_NAME         = 'CH4LS'//trim(iname),  &
      LONG_NAME          = 'CH4 Loss '//trim(iname),  &
      UNITS              = 'kg m-2 s-1', &
      DIMS               = MAPL_DimsHorzOnly,    &
      VLOCATION          = MAPL_VLocationNone,    &
      RC=STATUS  )
   _VERIFY(STATUS)

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

   ! USES:

   IMPLICIT NONE

   ! INPUT PARAMETERS:

   TYPE(Chem_Bundle), intent(in) :: w_c        ! Chemical tracer fields
   INTEGER, INTENT(IN) :: nymd, nhms           ! time
   REAL,    INTENT(IN) :: cdt              ! chemical timestep (secs)


   ! OUTPUT PARAMETERS:

   TYPE(CH4_GridComp1), INTENT(INOUT) :: gcCH4  ! Grid Component
   TYPE(ESMF_State), INTENT(INOUT)  :: impChem  ! Import State
   TYPE(ESMF_State), INTENT(INOUT)  :: expChem  ! Export State
   INTEGER, INTENT(OUT) ::  rc                  ! Error return code:
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

   CHARACTER(LEN=*), PARAMETER  :: Iam = 'CH4_GridCompInitialize1_'

   CHARACTER(LEN=ESMF_MAXSTR)   :: rcfilen
   type(ESMF_Config) :: cfg

   INTEGER :: j, n, status
   INTEGER :: i1, i2, im, j1, j2, jm, km
   INTEGER :: nTimes, begTime, incSecs
   INTEGER :: nbeg, nend, nymd1, nhms1
   LOGICAL :: NoRegionalConstraint

   REAL :: limitN, limitS
   REAL, ALLOCATABLE :: var2D(:,:)
   real, pointer     :: ptr2d(:,:) => null()
   real :: dummy_float
   logical :: dummy_bool
   character(LEN=ESMF_MAXSTR) :: dummy_string

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
      IF(MAPL_AM_I_ROOT()) PRINT *,TRIM(Iam)//": Must have only 1 bin at the single instance level"
      status = 1
      VERIFY_(status)
   end if

   !  Allocate memory, etc
   !  --------------------
   ALLOCATE ( gcCH4%eCH4(i1:i2,j1:j2),     &
!            gcCH4%eCH4_ind(i1:i2,j1:j2),   &
!            gcCH4%eCH4_wetl(i1:i2,j1:j2),  &
!            gcCH4%eCH4_agw(i1:i2,j1:j2),   &
!            gcCH4%eCH4_bb(i1:i2,j1:j2),    &
!            gcCH4%eCH4_bb_(i1:i2,j1:j2),   &
!            gcCH4%eCH4_bf(i1:i2,j1:j2),    &
!            gcCH4%eCH4_mnat(i1:i2,j1:j2),  &
            gcCH4%regionMask(i1:i2,j1:j2), &
            gcCH4%OHnd(i1:i2,j1:j2,km),    &
            gcCH4%Clnd(i1:i2,j1:j2,km),    &
            gcCH4%O1Dnd(i1:i2,j1:j2,km), STAT=status )
   VERIFY_(status)

   !  Load resource file
   !  ------------------
   cfg = ESMF_ConfigCreate()
   call ESMF_ConfigLoadFile(cfg, trim(rcfilen), rc=status)
   VERIFY_(status)

   !  Maximum allowed solar zenith angle for "daylight"
   !  -------------------------------------------------
   call ESMF_ConfigGetAttribute(cfg, value=dummy_float, label='solar_ZA_cutoff:', rc=status)
   VERIFY_(status)
   call ESMF_ConfigGetAttribute(cfg, value=gcCH4%szaCutoff, label='solar_ZA_cutoff.'//trim(gcCH4%iname)//':', default=dummy_float, rc=status)
   VERIFY_(status)

   !  Run-time debug switch
   !  ---------------------
   call ESMF_ConfigGetAttribute(cfg, value=dummy_bool, label='DEBUG:', default=.false., rc=status)
   VERIFY_(status)
   call ESMF_ConfigGetAttribute(cfg, value=gcCH4%DebugIsOn, label='DEBUG.'//trim(gcCH4%iname)//':', default=dummy_bool, rc=status)
   VERIFY_(status)

   !  Methane photolysis feedback switch
   !  ----------------------------------
   call ESMF_ConfigGetAttribute(cfg, value=dummy_bool, label='CH4_Feedback:', default=.false., rc=status)
   VERIFY_(status)
   call ESMF_ConfigGetAttribute(cfg, value=gcCH4%CH4FeedBack, label='CH4_Feedback.'//trim(gcCH4%iname)//':', default=dummy_bool, rc=status)
   VERIFY_(status)

   !  Water vapor feedback switch
   !  ---------------------------
   call ESMF_ConfigGetAttribute(cfg, value=dummy_bool, label='H2O_Feedback:', default=.false., rc=status)
   VERIFY_(status)
   call ESMF_ConfigGetAttribute(cfg, value=gcCH4%H2OFeedBack, label='H2O_Feedback.'//trim(gcCH4%iname)//':', default=dummy_bool, rc=status)
   VERIFY_(status)

   ! Should we do chemistry?
   !  ----------------------
   call ESMF_ConfigGetAttribute(cfg, value=dummy_bool, label='do_chemistry:', default=.true., rc=status)
   VERIFY_(STATUS)
   call ESMF_ConfigGetAttribute(cfg, value=gcCH4%chemistry, label='do_chemistry.'//trim(gcCH4%iname)//':', default=dummy_bool, rc=status)
   VERIFY_(STATUS)
   if (MAPL_AM_I_ROOT() .and. (.not. gcCH4%chemistry)) write(*,'("Chemistry turned off for ", a)') trim(gcCH4%iname)

   call MAPL_GetPointer(impChem,ptr2D,'CH4_regionMask',rc=status)
   VERIFY_(STATUS)

   !  Grab the region string.
   !  -----------------------
   call ESMF_ConfigGetAttribute(cfg, value=dummy_string, label='CH4_regions_indices:', rc=status)
   VERIFY_(status)
   call ESMF_ConfigGetAttribute(cfg, value=gcCH4%regionsString, label='CH4_regions_indices.'//trim(gcCH4%iname)//':', default=dummy_string, rc=status)
   VERIFY_(status)

   !  Is this instantiation a global case?
   !  -----------------------------------
   IF(gcCH4%regionsString(1:2) == "-1") THEN
      NoRegionalConstraint = .TRUE.
   ELSE
      SELECT CASE (ESMF_UtilStringLowerCase(gcCH4%regionsString(1:2)))
      CASE ("gl")
         NoRegionalConstraint = .TRUE.
      CASE ("al")
         NoRegionalConstraint = .TRUE.
      CASE DEFAULT
         NoRegionalConstraint = .FALSE.
      END SELECT
   END IF

   !  Set regionsString to "-1" for the global case
   !  ---------------------------------------------
   IF(NoRegionalConstraint) gcCH4%regionsString = "-1"

   IF(MAPL_AM_I_ROOT()) THEN
      IF(NoRegionalConstraint) THEN
         PRINT *,TRIM(Iam)//": This instantiation has no regional constraints."
      ELSE
         PRINT *,TRIM(Iam)//": This instantiation is regionally constrained."
         PRINT *,TRIM(Iam)//": List of region numbers included: ",trim(gcCH4%regionsString)
      END IF
   END IF

   !  Use instance name as key to CH4 emission source
   !  -----------------------------------------------
   gcCH4%CH4Source = "CH4_"//trim(gcCH4%iname)

   return

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

subroutine CH4_GridCompRun1_ ( gcCH4, w_c, impChem, expChem, nymd, nhms, cdt, rc )

   !USES:

   IMPLICIT NONE

   ! INPUT/OUTPUT PARAMETERS:

   TYPE(CH4_GridComp1), INTENT(INOUT) :: gcCH4   ! Grid Component
   TYPE(Chem_Bundle), INTENT(INOUT) :: w_c  ! Chemical tracer fields

   ! INPUT PARAMETERS:

   TYPE(ESMF_State), INTENT(inout) :: impChem    ! Import State
   INTEGER, INTENT(IN) :: nymd, nhms          ! time
   REAL,    INTENT(IN) :: cdt             ! chemical timestep (secs)


   ! OUTPUT PARAMETERS:

   TYPE(ESMF_State), intent(inout) :: expChem     ! Export State
   INTEGER, INTENT(OUT) ::  rc                  ! Error return code:
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

   CHARACTER(LEN=*), PARAMETER :: Iam = 'CH4_GridCompRun1_'

   !  Input fields from fvGCM
   !  -----------------------
   REAL, POINTER, DIMENSION(:,:)   ::  cellArea => null()
   REAL, POINTER, DIMENSION(:,:)   ::  pblh     => null()
   REAL, POINTER, DIMENSION(:,:,:) ::  zle      => null()
   REAL, POINTER, DIMENSION(:,:,:) ::  T        => null()
   REAL, POINTER, DIMENSION(:,:,:) ::  Q        => null()
   REAL, POINTER, DIMENSION(:,:,:) ::  qtot     => null()
   REAL, POINTER, DIMENSION(:,:,:) ::  rhowet   => null()
   REAL, POINTER, DIMENSION(:,:,:) ::  ple      => null()

   INTEGER :: i1, i2, im, j1, j2, jm, km, idiag, iXj
   INTEGER :: i, j, k, kReverse, n, nbeg, nend
   INTEGER :: nymd1, nhms1, status

   REAL    :: qmin, qmax

   REAL, ALLOCATABLE :: cellDepth(:,:,:), cellVolume(:,:,:)
   REAL, ALLOCATABLE ::      rkoh(:,:,:),       rkcl(:,:,:),    rko1d(:,:,:),   rktot(:,:,:)
   REAL, ALLOCATABLE ::         p(:,:,:),      ndwet(:,:,:)
   REAL, ALLOCATABLE ::    dCH4ox(:,:,:),      photJ(:,:,:), dCH4Phot(:,:,:)

   CHARACTER(LEN=256) :: CH4Source

   REAL, POINTER :: ptr2d(:,:)   => null()
   REAL, POINTER :: ptr3d(:,:,:) => null()

#define EXPORT   expChem
#define iNAME    TRIM(gcCH4%iname)

#define CH4EM    CH4_emis
!#define CH4CL    CH4_column
!#define CH4SC    CH4_surface
!#define CH4PD    CH4_prod
#define CH4LS    CH4_loss
!#define CH4JL    CH4_phot
!#define CH4QP    CH4_qprod
!#define CH4DRY   CH4_dry ! Sourish

!#include "CH4_GetPointer___.h"
   real, pointer, dimension(:,:)   :: CH4EM ! EXPORT: CH4 Emission
   real, pointer, dimension(:,:)   :: CH4PD ! EXPORT: CH4 Chemical Production
   real, pointer, dimension(:,:)   :: CH4LS ! EXPORT: CH4 Chemical Loss
   real, pointer, dimension(:,:)   :: CH4SC ! EXPORT: CH4 Surface Concentration
   real, pointer, dimension(:,:)   :: CH4CL ! EXPORT: CH4 Column Burden
   real, pointer, dimension(:,:,:) :: CH4JL ! EXPORT: CH4 Photolytic Loss
   real, pointer, dimension(:,:,:) :: CH4QP ! EXPORT: H2O tendency from CH4 photolysis
   real, pointer, dimension(:,:,:) :: CH4DRY ! EXPORT: CH4_dry_air_mole_fraction
   real, pointer, dimension(:,:,:) :: CH4_for_rad ! come up with a better name later ! Sourish

   call MAPL_GetPointer ( EXPORT, CH4EM,  'CH4EM'//iNAME, RC=STATUS )
   _VERIFY(STATUS)
   !call MAPL_GetPointer ( EXPORT, CH4PD,  'CH4PD'//iNAME, RC=STATUS )
   !_VERIFY(STATUS)
   call MAPL_GetPointer ( EXPORT, CH4LS,  'CH4LS'//iNAME, RC=STATUS )
   _VERIFY(STATUS)
   !call MAPL_GetPointer ( EXPORT, CH4SC,  'CH4SC'//iNAME, RC=STATUS )
   !_VERIFY(STATUS)
   !call MAPL_GetPointer ( EXPORT, CH4CL,  'CH4CL'//iNAME, RC=STATUS )
   !_VERIFY(STATUS)
   !call MAPL_GetPointer ( EXPORT, CH4JL,  'CH4JL'//iNAME, RC=STATUS )
   !_VERIFY(STATUS)
   !call MAPL_GetPointer ( EXPORT, CH4QP,  'CH4QP'//iNAME, RC=STATUS )
   !_VERIFY(STATUS)
   !call MAPL_GetPointer ( EXPORT, CH4DRY,  'CH4DRY'//iNAME, RC=STATUS ) ! Sourish
   !_VERIFY(STATUS) ! Sourish


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
         print *, trim(iNAME)//": Must have only 1 bin at the single instance level"
      end if
      status = 1
      VERIFY_(status)
   end if

   !  Get imports
   !  -----------
   CALL MAPL_GetPointer(impChem, PBLH,     'ZPBL',        RC=status)
   VERIFY_(status)
   CALL MAPL_GetPointer(impChem, ZLE,      'ZLE',         RC=status)
   VERIFY_(status)
   CALL MAPL_GetPointer(impChem, T,    'T',           RC=status)
   VERIFY_(status)
   CALL MAPL_GetPointer(impChem, Q,    'Q',           RC=status)
   VERIFY_(status)
   !CALL MAPL_GetPointer(impChem, qtot,     'QTOT',        RC=status) ! Sourish
   !VERIFY_(status) ! Sourish
   CALL MAPL_GetPointer(impChem, rhowet,   'AIRDENS',     RC=status)
   VERIFY_(status)
   CALL MAPL_GetPointer(impChem, cellArea, 'AREA',        RC=status)
   VERIFY_(status)
   CALL MAPL_GetPointer(impChem, ple,      'PLE',         RC=status)
   VERIFY_(status)

   if (gcCH4%DebugIsOn) then
      call pmaxmin('CH4:AREA', cellArea, qmin, qmax, iXj,  1,   1. )
      call pmaxmin('CH4:ZPBL',     pblh, qmin, qmax, iXj, km+1, 1. )
      call pmaxmin('CH4:ZLE',       zle, qmin, qmax, iXj, km+1, 1. )
      call pmaxmin('CH4:T',           T, qmin, qmax, iXj, km,   1. )
      call pmaxmin('CH4:Q',           q, qmin, qmax, iXj, km,   1. )
      !call pmaxmin('CH4:QTOT',     qtot, qmin, qmax, iXj, km,   1. ) ! Sourish
      call pmaxmin('CH4:RHOWET', rhowet, qmin, qmax, iXj, km,   1. )
      call pmaxmin('CH4:PLE',       ple, qmin, qmax, iXj, km+1, 1. )
   endif

   !  Update CH4 emissions and OH, Cl, and O(1D) number densities
   !  (in molec cm^-3) once each day
   !  ------------------------------------------------------------

   ! All sources
   call MAPL_GetPointer(impChem, ptr2d, 'CH4_'//trim(iNAME), rc=status)
   VERIFY_(STATUS)
   gcCH4%eCH4 = ptr2d


   !  Background OH, Cl, and O(1D) for loss term
   !  ------------------------------------------
!   call MAPL_GetPointer(impChem,ptr3d,'CH4_oh'//trim(iNAME),rc=status)
!   VERIFY_(STATUS)
!   gcCH4%OHnd=ptr3d

!   call MAPL_GetPointer(impChem,ptr3d,'CH4_cl'//trim(iNAME),rc=status)
!   VERIFY_(STATUS)
!   gcCH4%Clnd=ptr3d

!   call MAPL_GetPointer(impChem,ptr3d,'CH4_o1d'//trim(iNAME),rc=status)
!   VERIFY_(STATUS)
!   gcCH4%O1Dnd=ptr3d

   call MAPL_GetPointer(impChem,ptr3d,'CH4_oh', rc=status)
   VERIFY_(STATUS)
   gcCH4%OHnd=ptr3d

   call MAPL_GetPointer(impChem,ptr3d,'CH4_cl', rc=status)
   VERIFY_(STATUS)
   gcCH4%Clnd=ptr3d

   call MAPL_GetPointer(impChem,ptr3d,'CH4_o1d', rc=status)
   VERIFY_(STATUS)
   gcCH4%O1Dnd=ptr3d

   !  Convert number densities from molec cm^-3 to molec m^-3
   !  -------------------------------------------------------
   gcCH4%OHnd(i1:i2,j1:j2,1:km) = gcCH4%OHnd(i1:i2,j1:j2,1:km)*1.00E+06
   gcCH4%Clnd(i1:i2,j1:j2,1:km)  =  gcCH4%Clnd(i1:i2,j1:j2,1:km)*1.00E+06
   gcCH4%O1Dnd(i1:i2,j1:j2,1:km) = gcCH4%O1Dnd(i1:i2,j1:j2,1:km)*1.00E+06

   !  Allocate temporary workspace
   !  ----------------------------
   allocate( p(i1:i2,j1:j2,km), &
         ndwet(i1:i2,j1:j2,km), &
          rkoh(i1:i2,j1:j2,km), &
          rkcl(i1:i2,j1:j2,km), &
         rko1d(i1:i2,j1:j2,km), &
         rktot(i1:i2,j1:j2,km), &
        dCH4ox(i1:i2,j1:j2,km), &
    cellVolume(i1:i2,j1:j2,km), &
     cellDepth(i1:i2,j1:j2,km), STAT=status)
   VERIFY_(status)

   !  Layer mean pressures
   !  --------------------
   do k=1,km
      p(i1:i2,j1:j2,k) = (ple(i1:i2,j1:j2,k-1)+ple(i1:i2,j1:j2,k))*0.50
   end do

   !  Wet-air number density
   !  ----------------------
   ndwet(i1:i2,j1:j2,1:km) = rhowet(i1:i2,j1:j2,1:km)*MAPL_AVOGAD/MAPL_AIRMW

   !  Cell depth and volume
   !  ---------------------
   do k=1,km
      cellDepth(i1:i2,j1:j2,k)  = (ple(i1:i2,j1:j2,k)-ple(i1:i2,j1:j2,k-1)) &
                                / (rhowet(i1:i2,j1:j2,k)*MAPL_GRAV)
      cellVolume(i1:i2,j1:j2,k) = cellArea(i1:i2,j1:j2)*cellDepth(i1:i2,j1:j2,k)
   end do

!   if (gcCH4%DebugIsOn) then
!      call pmaxmin('CH4: eCH4_ind',   gcCH4%eCH4_ind, qmin, qmax, iXj, 1,  1. )
!      call pmaxmin('CH4: eCH4_agw',   gcCH4%eCH4_agw, qmin, qmax, iXj, 1,  1. )
!      call pmaxmin('CH4: eCH4_bb',     gcCH4%eCH4_bb, qmin, qmax, iXj, 1,  1. )
!      call pmaxmin('CH4: eCH4_bf',     gcCH4%eCH4_bf, qmin, qmax, iXj, 1,  1. )
!      call pmaxmin('CH4: eCH4_wetl', gcCH4%eCH4_wetl, qmin, qmax, iXj, 1,  1. )
!      call pmaxmin('CH4: eCH4_mnat', gcCH4%eCH4_mnat, qmin, qmax, iXj, 1,  1. )
!      call pmaxmin('CH4: OH Conc',        gcCH4%OHnd, qmin, qmax, iXj, 1,  1. )
!      call pmaxmin('CH4: Cl Conc',        gcCH4%Clnd, qmin, qmax, iXj, 1,  1. )
!      call pmaxmin('CH4: O(1D) Conc',    gcCH4%O1Dnd, qmin, qmax, iXj, 1,  1. )
!      call pmaxmin('CH4: Cell Vol',       cellVolume, qmin, qmax, iXj, km, 1. )
!      call pmaxmin('CH4: Cell Depth',      cellDepth, qmin, qmax, iXj, km, 1. )
!      call pmaxmin('CH4: Wet-air ND',          ndwet, qmin, qmax, iXj, km, 1. )
!   end if

   !  Compute and add surface emissions
   !  ---------------------------------
   call CH4_Emission(rc)

   !  Loss rate [m^3 s^-1] for OH + CH4 => CO + products
   !  --------------------------------------------------
   rkoh(i1:i2,j1:j2,1:km)  = 2.45E-12*1.00E-06*exp(-1775./T(i1:i2,j1:j2,1:km))

   !  Loss rate [m^3 s^-1] for Cl + CH4 => CO + products
   !  --------------------------------------------------
   rkcl(i1:i2,j1:j2,1:km)  = 7.10E-12*1.00E-06*exp(-1270./T(i1:i2,j1:j2,1:km))

   !  Loss rate [m^3 s^-1] for O(1D) + CH4 => CO + products
   !  -----------------------------------------------------
   rko1d(i1:i2,j1:j2,1:km) = 1.75E-10*1.00E-06

   !  Change in CH4 mole fraction due to oxidation
   !  --------------------------------------------
!   dCH4ox(i1:i2,j1:j2,1:km) = w_c%qa(nbeg)%data3d(i1:i2,j1:j2,1:km)                         &
!                              * (    rkoh(i1:i2,j1:j2,1:km)*gcCH4%OHnd(i1:i2,j1:j2,1:km)    &
!                                  +  rkcl(i1:i2,j1:j2,1:km)*gcCH4%Clnd(i1:i2,j1:j2,1:km)    &
!                                  + rko1d(i1:i2,j1:j2,1:km)*gcCH4%O1Dnd(i1:i2,j1:j2,1:km) )

    if (gcCH4%chemistry) then
        rktot = rkoh*gcCH4%OHnd(i1:i2,j1:j2,1:km) + rkcl*gcCH4%Clnd(i1:i2,j1:j2,1:km) + rko1d*gcCH4%O1Dnd(i1:i2,j1:j2,1:km)
        dCH4ox = w_c%qa(nbeg)%data3d(i1:i2,j1:j2,1:km) * (exp(-cdt*rktot)-1.0)
        !w_c%qa(nbeg)%data3d(i1:i2,j1:j2,1:km) = w_c%qa(nbeg)%data3d(i1:i2,j1:j2,1:km) + dCH4ox
        !w_c%qa(nbeg)%data3d(i1:i2,j1:j2,1:km) = merge(w_c%qa(nbeg)%data3d(i1:i2,j1:j2,1:km), 0.0, w_c%qa(nbeg)%data3d(i1:i2,j1:j2,1:km) > 0.0)
        w_c%qa(nbeg)%data3d(i1:i2,j1:j2,1:km) = w_c%qa(nbeg)%data3d(i1:i2,j1:j2,1:km) * exp(-cdt*rktot) ! avoid negatives
    else
        dCH4ox = 0.0
    end if

   !  Vertically integrated CH4 loss due to oxidation (only)
   !  ------------------------------------------------------
   if (associated(CH4_loss)) then
      CH4_loss(i1:i2,j1:j2) = 0.
      do k = 1,km
         CH4_loss(i1:i2,j1:j2) = CH4_loss(i1:i2,j1:j2) &
                               + dCH4ox(i1:i2,j1:j2,k) * (mwtCH4/MAPL_AIRMW) * w_c%delp(i1:i2,j1:j2,k)/MAPL_GRAV ! kg/m^2 over this time step
      end do
      CH4_loss(i1:i2,j1:j2) = CH4_loss(i1:i2,j1:j2)/cdt ! convert to kg/m^2/s, then history can export a time-averaged value
   end if

   !  CH4 production (none)
   !  ---------------------
!   if (associated(CH4_prod)) CH4_prod(i1:i2,j1:j2) = 0.

   !  Calculate photolytic loss rates, J [s^-1], for CH4 + hv => 2H2O + CO
   !  Notice that J and the losses are always computed. However, the setting
   !  of the feedback switch(es) determines if the increments are actually applied
   !  ----------------------------------------------------------------------------
   ALLOCATE(photJ(i1:i2,j1:j2,1:km), dCH4Phot(i1:i2,j1:j2,1:km), STAT=status)
   VERIFY_(STATUS)

   photJ(:,:,:) = 0.
   CALL getJRates(status)
   VERIFY_(status)

   !  Change in CH4 number density [m^-3 s^-1] due to photolysis
   !  ----------------------------------------------------------
   dCH4Phot(i1:i2,j1:j2,1:km) = photJ(i1:i2,j1:j2,1:km)*w_c%qa(nbeg)%data3d(i1:i2,j1:j2,1:km)

!   if (associated(CH4_phot)) THEN
!      CH4_phot(i1:i2,j1:j2,1:km) = dCH4Phot(i1:i2,j1:j2,1:km)*ndwet(i1:i2,j1:j2,1:km)
!   endif

!   !  Increment the CH4 mole fraction due to oxidation
!   !  ------------------------------------------------
!   w_c%qa(nbeg)%data3d(i1:i2,j1:j2,1:km) = w_c%qa(nbeg)%data3d(i1:i2,j1:j2,1:km) - cdt*dCH4ox(i1:i2,j1:j2,1:km)

   !  Increment the CH4 mole fraction due to photolysis when the switch is on
   !  -----------------------------------------------------------------------
   if (gcCH4%CH4FeedBack .and. gcCH4%chemistry) then
      w_c%qa(nbeg)%data3d(i1:i2,j1:j2,1:km) = w_c%qa(nbeg)%data3d(i1:i2,j1:j2,1:km) - cdt*dCH4Phot(i1:i2,j1:j2,1:km)
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

   !  Surface concentration [ppbv]
   !  ----------------------------
!   if (associated(CH4_surface)) then
!      CH4_surface(i1:i2,j1:j2) = w_c%qa(nbeg)%data3d(i1:i2,j1:j2,km)*1.00E+09
!   end if

   !  Column burden [kg m-2]
   !  ----------------------
!   if (associated(CH4_column)) then
!      CH4_column(i1:i2,j1:j2) = 0.
!      do k = 1,km
!         CH4_column(i1:i2,j1:j2) = CH4_column(i1:i2,j1:j2)                              &
!                                + w_c%qa(nbeg)%data3d(i1:i2,j1:j2,k)*mwtCH4/MAPL_AIRMW &
!                                           * w_c%delp(i1:i2,j1:j2,k)/MAPL_GRAV
!      end do
!   end if

!! Sourish Basu
   !!  Dry-air mole fraction [mol mol-1]
   !!  ---------------------------------
   !if (associated(CH4_dry)) then
      !CH4_dry(i1:i2,j1:j2,1:km) = w_c%qa(nbeg)%data3d(i1:i2,j1:j2,1:km) / (1. - qtot(i1:i2,j1:j2,1:km))
   !end if

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

! Sourish Basu
   if (trim(iNAME) == "total") then
      call MAPL_GetPointer ( EXPORT, CH4_for_rad, 'CH4', RC=STATUS )
      _VERIFY(STATUS)
      if (associated(CH4_for_rad)) then
         CH4_for_rad(i1:i2,j1:j2,1:km) = w_c%qa(nbeg)%data3d(i1:i2,j1:j2,1:km)
      end if
   end if

   !  Housekeeping
   !  ------------
   DEALLOCATE(p, rkoh, rkcl, rko1d, rktot, dCH4ox, cellDepth, cellVolume, STAT=status)
   VERIFY_(status)
   DEALLOCATE(photJ, dCH4Phot, STAT=status)
   VERIFY_(status)

   RETURN

   CONTAINS

   !-------------------------------------------------------------------------
   subroutine CH4_Emission ( rc )

      IMPLICIT NONE

      ! OUTPUT VARIABLES
      INTEGER, INTENT(OUT) :: rc ! Error return code

      ! LOCAL VARIABLES
      CHARACTER(LEN=*), PARAMETER :: myname = 'CH4_Emission'

      INTEGER ::  i, j, k, kt, minkPBL
      INTEGER, ALLOCATABLE :: index(:)

      REAL, ALLOCATABLE, dimension(:,:) :: pblLayer, sfcFlux, myMask, mf_to_be_added, mf_current
      REAL, ALLOCATABLE :: fPBL(:,:,:)

      rc = 0

      allocate(sfcFlux(i1:i2,j1:j2),      STAT=rc)      ! emissions
      allocate( myMask(i1:i2,j1:j2),      STAT=rc)      ! region mask
      allocate(   fPBL(i1:i2,j1:j2,1:km), STAT=rc)      ! partitioning of BB
      allocate(mf_to_be_added(i1:i2,j1:j2), mf_current(i1:i2,j1:j2), STAT=rc)

!      ! Apply biomass burning diurnal cycle if desired
!      ! ----------------------------------------------
!   SB: Commented out for now because I'm not sure how the daily average is preserved in successive calls
!      if (w_c%diurnal_bb) then
!         allocate(gcCH4%eCH4_d(i1:i2,j1:j2))
!         gcCH4%eCH4_d = gcCH4%eCH4 ! temp storage, put the daily average in eCH4_d
!!         gcCH4%eCH4_bb_(:,:) = gcCH4%eCH4_bb(:,:)

!!         call Chem_BiomassDiurnal ( gcCH4%eCH4_bb, gcCH4%eCH4_bb_,   &
!!                                w_c%grid%lon(:,:)*radToDeg,      &
!!                                w_c%grid%lat(:,:)*radToDeg, nhms, cdt )
!         call Chem_BiomassDiurnal ( gcCH4%eCH4, gcCH4%eCH4_d,   &
!                                w_c%grid%lon(:,:)*radToDeg,      &
!                                w_c%grid%lat(:,:)*radToDeg, nhms, cdt )
!         deallocate(gcCH4%eCH4_d)
!      end if

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
      kt = minkPBL

!      Layer: do k = kt,km

!         !    Emissions: Weighted biomass burning
!         !    -----------------------------------
!         sfcFlux(i1:i2,j1:j2) = gcCH4%eCH4_BB(i1:i2,j1:j2)*fPBL(i1:i2,j1:j2,k)

!         ! Add other emission components when in surface layer
!         ! ---------------------------------------------------
!         if (k == km) then
!            sfcFlux(i1:i2,j1:j2) = sfcFlux(i1:i2,j1:j2) &
!                                    + gcCH4%eCH4_bf(i1:i2,j1:j2)   &
!                                    + gcCH4%eCH4_ind(i1:i2,j1:j2)  &
!                                    + gcCH4%eCH4_agw(i1:i2,j1:j2)  &
!                                    + gcCH4%eCH4_mnat(i1:i2,j1:j2) &
!                                    + gcCH4%eCH4_wetl(i1:i2,j1:j2)
!         end if

!      !    Update CH4 at this level
!      !    ------------------------
!         w_c%qa(nbeg)%data3d(i1:i2,j1:j2,k) = w_c%qa(nbeg)%data3d(i1:i2,j1:j2,k)           &
!                                             + cdt * sfcFlux(i1:i2,j1:j2)*MAPL_AIRMW/mwtCH4 &
!                                                   / (w_c%delp(i1:i2,j1:j2,k)/MAPL_GRAV)
!      end do Layer

      ! debug
      mf_current = w_c%qa(nbeg)%data3d(i1:i2,j1:j2,km)
      mf_to_be_added = cdt * gcCH4%eCH4(i1:i2,j1:j2)*MAPL_AIRMW / mwtCH4 / (w_c%delp(i1:i2,j1:j2,km)/MAPL_GRAV)
      if (any(mf_current+mf_to_be_added .le. 0.0)) then
         !do j = j1,j2
            !do i = i1,i2
               !if (mf_current(i,j)+mf_to_be_added(i,j) .le. 0.0) then
                  !write(*,'("NEG MF :: Tag ", a, " at cell i=", i4, ", j=", i4, " has mf = ", es14.7, ", after adding ", es14.7, " will end up with negative MF")') &
                     !trim(gcCH4%iname), i, j, mf_current(i,j), mf_to_be_added(i,j)
               !end if
            !end do
         !end do
         write(*,'("For tag ", a, ", pixels where adding emissions makes MF negative = ", i7)') trim(gcCH4%iname), count(mf_current+mf_to_be_added .le. 0.0)
      end if
      ! end debug

      w_c%qa(nbeg)%data3d(i1:i2,j1:j2,km) = w_c%qa(nbeg)%data3d(i1:i2,j1:j2,km) + &
            cdt * gcCH4%eCH4(i1:i2,j1:j2)*MAPL_AIRMW / mwtCH4 / (w_c%delp(i1:i2,j1:j2,km)/MAPL_GRAV)


      ! Update Surface flux diagnostic for this bin
      ! -------------------------------------------
      if (associated(CH4_emis)) then
         CH4_emis(i1:i2,j1:j2) = gcCH4%eCH4(i1:i2,j1:j2)
!!         CH4_emis(i1:i2,j1:j2) = gcCH4%eCH4_bb(i1:i2,j1:j2)   + gcCH4%eCH4_bf(i1:i2,j1:j2)   &
!!                           + gcCH4%eCH4_ind(i1:i2,j1:j2)  + gcCH4%eCH4_agw(i1:i2,j1:j2)  &
!!                           + gcCH4%eCH4_mnat(i1:i2,j1:j2) + gcCH4%eCH4_wetl(i1:i2,j1:j2)
      end if

      deallocate(fPBL, myMask, sfcFlux, mf_to_be_added, mf_current, STAT=rc)

      return

   end subroutine CH4_Emission

   !-------------------------------------------------------------------------
   ! Borrowed from meso_phot.F of StratChem, where number densities are cgs [cm^{-3}]
   SUBROUTINE getJRates(rc)

      IMPLICIT NONE

      INTEGER, INTENT(out) :: rc

      REAL, ALLOCATABLE :: o2Column(:,:,:)
      REAL, ALLOCATABLE :: SZARad(:,:)
      REAL, ALLOCATABLE :: SZADeg(:,:)
      REAL, ALLOCATABLE :: sinSZA(:,:)
      REAL, ALLOCATABLE :: zgrz(:,:)
      REAL, ALLOCATABLE :: sfaca(:,:)
      REAL, ALLOCATABLE :: arg(:,:)

      REAL, PARAMETER :: wavel = 1215.7
      REAL, PARAMETER :: O2xs  = 1.000E-20
      REAL, PARAMETER :: CH4xs = 2.000E-17
      REAL, PARAMETER :: sflux = 4.006E+11

      ! Constants for Chapman function at high solar zenith angle
      ! ---------------------------------------------------------
      REAL, PARAMETER :: hbar = 6.79
      REAL, PARAMETER :: zbar = 30.0
      REAL, PARAMETER :: r0   = 6.371E+03
      REAL, PARAMETER :: zp   = 60.0

      REAL, PARAMETER :: d1 = 1.060693
      REAL, PARAMETER :: d2 = 0.55643831
      REAL, PARAMETER :: d3 = 1.0619896
      REAL, PARAMETER :: d4 = 1.7245609
      REAL, PARAMETER :: d5 = 0.56498823
      REAL, PARAMETER :: d6 = 0.06651874

      REAL, PARAMETER :: O2Abv80km = 7.072926E+19 ![cm^{-2}]
      REAL, PARAMETER :: O2VMR = 0.20946

      REAL :: b
      REAL :: f
      REAL :: r
      REAL :: s

      INTEGER :: status

      CHARACTER(LEN=ESMF_MAXSTR) :: Iam = "CH4::getJRates"

      rc = 0
      b = SQRT(0.50*r0/hbar)

      ! O2 overhead number density profile [cm^{-2}]
      ! --------------------------------------------
      ALLOCATE(O2Column(i1:i2,j1:j2,1:km), STAT=status)
      VERIFY_(status)

      f = O2VMR*5.00E-05
      O2Column(:,:,1) = O2Abv80km+cellDepth(:,:,1)*ndwet(:,:,1)*f

      DO k = 2,km
         O2Column(:,:,k) = O2Column(:,:,k-1)+(cellDepth(:,:,k-1)*ndwet(:,:,k-1)+ &
                                        cellDepth(:,:,  k)*ndwet(:,:,  k))*f
      END DO

      IF(gcCH4%DebugIsOn) THEN
         CALL pmaxmin('CH4: O2Column', O2Column, qmin, qmax, iXj, km,  1. )
      END IF

      ! Grab some memory
      ! ----------------
      ALLOCATE(SZARad(i1:i2,j1:j2), STAT=status)
      VERIFY_(status)
      ALLOCATE(SZADeg(i1:i2,j1:j2), STAT=status)
      VERIFY_(status)
      ALLOCATE(sinSZA(i1:i2,j1:j2), STAT=status)
      VERIFY_(status)

      WHERE(w_c%cosz(i1:i2,j1:j2) > 1.00)
         SZARad(i1:i2,j1:j2) = 0.00
      ELSEWHERE
         SZARad(i1:i2,j1:j2) = ACOS(w_c%cosz(i1:i2,j1:j2))
      ENDWHERE
      SZADeg(i1:i2,j1:j2) = SZARad(i1:i2,j1:j2)*radToDeg
      sinSZA(i1:i2,j1:j2) = SIN(SZARad(i1:i2,j1:j2))

      ALLOCATE(zgrz(i1:i2,j1:j2), STAT=status)
      VERIFY_(status)

      WHERE(SZADeg(i1:i2,j1:j2) <= 90.00)
         zgrz(i1:i2,j1:j2) = 1000.00
      ELSEWHERE
         zgrz(i1:i2,j1:j2) = sinSZA(i1:i2,j1:j2)*(zp+r0)-r0
      ENDWHERE

      IF(gcCH4%DebugIsOn) THEN
         CALL pmaxmin('CH4: zgrz', zgrz, qmin, qmax, iXj, 1,  1. )
         CALL pmaxmin('CH4: cosz', w_c%cosz, qmin, qmax, iXj, 1,  1. )
      END IF

      ALLOCATE(sfaca(i1:i2,j1:j2), STAT=status)
      VERIFY_(status)
      sfaca(i1:i2,j1:j2) = 0.00

      ! Chapman function calculation from ACDB 2-D model
      ! ------------------------------------------------
      DO j = j1,j2
         DO i = i1,i2

            Daytime: IF(SZADeg(i,j) < gcCH4%szaCutoff) THEN

               IF(SZADeg(i,j) < 70.00) THEN

                  sfaca(i,j) = 1.00/w_c%cosz(i,j)

               ELSE IF(zgrz(i,j) > 0.00) THEN

                  s = b*ABS(w_c%cosz(i,j))

                  IF(s <= 8.00) THEN
                     s = (d1+d2*s)/(d3+d4*s+s**2)
                  ELSE
                     s = d5/(d6+s)
                  END IF

                  r = b*SQRT(MAPL_PI)
                  sfaca(i,j) = r*s

                  IF(SZADeg(i,j) > 90.00) THEN
                     sfaca(i,j) = 2.00*r*EXP((r0+zbar)*(1.00-sinSZA(i,j))/hbar)-sfaca(i,j)
                  END IF

               END IF

            END IF Daytime

         END DO
      END DO

      IF(gcCH4%DebugIsOn) THEN
         CALL pmaxmin('CH4: sfaca', sfaca, qmin, qmax, iXj, 1,  1. )
      END IF

      ALLOCATE(arg(i1:i2,j1:j2), STAT=status)
      VERIFY_(status)

      ! At each layer, compute the rate constant, J [s^{-1}], if the sun is up
      ! ----------------------------------------------------------------------
      DO k = 1,km

         WHERE(SZADeg(i1:i2,j1:j2) < gcCH4%szaCutoff)
            arg(i1:i2,j1:j2) = O2Column(i1:i2,j1:j2,k)*O2xs*sfaca(i1:i2,j1:j2)
            photJ(i1:i2,j1:j2,k) = sflux*EXP(-arg(i1:i2,j1:j2))*CH4xs
         END WHERE

      END DO

      IF(gcCH4%DebugIsOn) THEN
         CALL pmaxmin('CH4: photJ', photJ, qmin, qmax, iXj, km,  1. )
      END IF

      DEALLOCATE(SZARad, STAT=status)
      VERIFY_(status)
      DEALLOCATE(SZADeg, STAT=status)
      VERIFY_(status)
      DEALLOCATE(sinSZA, STAT=status)
      VERIFY_(status)
      DEALLOCATE(zgrz, STAT=status)
      VERIFY_(status)
      DEALLOCATE(sfaca, STAT=status)
      VERIFY_(status)
      DEALLOCATE(arg, STAT=status)
      VERIFY_(status)
      DEALLOCATE(O2Column, STAT=status)
      VERIFY_(status)

      return

   end subroutine getJRates

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

   ! USES:

   IMPLICIT NONE

   ! INPUT/OUTPUT PARAMETERS:

   TYPE(CH4_GridComp1), INTENT(INOUT) :: gcCH4   ! Grid Component

   ! INPUT PARAMETERS:

   TYPE(Chem_Bundle), INTENT(IN)  :: w_c      ! Chemical tracer fields
   INTEGER, INTENT(IN) :: nymd, nhms          ! time
   REAL,    INTENT(IN) :: cdt             ! chemical timestep (secs)


   ! OUTPUT PARAMETERS:

   TYPE(ESMF_State), INTENT(INOUT) :: impChem   ! Import State
   TYPE(ESMF_State), INTENT(INOUT) :: expChem   ! Import State
   INTEGER, INTENT(OUT) ::  rc                  ! Error return code:
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

   CHARACTER(LEN=*), PARAMETER :: Iam = 'CH4_GridCompFinalize1_'
   rc = 0

   DEALLOCATE(gcCH4%eCH4,                                         &
!              gcCH4%eCH4_ind, gcCH4%eCH4_wetl,  gcCH4%eCH4_agw,  &
!              gcCH4%eCH4_bb,  gcCH4%eCH4_bb_,   gcCH4%eCH4_mnat, &
!              gcCH4%eCH4_bf,  gcCH4%regionMask, gcCH4%OHnd,      &
!              gcCH4%Clnd,     gcCH4%O1Dnd,                       &
               gcCH4%regionMask, gcCH4%OHnd, gcCH4%Clnd, gcCH4%O1Dnd, &
               STAT=rc)
   VERIFY_(rc)

   return

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

   IMPLICIT NONE

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

   TYPE(Chem_Bundle), intent(inout) :: w_c     ! Chemical tracer fields
   INTEGER, INTENT(IN) :: nymd, nhms           ! time
   REAL,    INTENT(IN) :: cdt              ! chemical timestep (secs)


   ! OUTPUT PARAMETERS:

   TYPE(CH4_GridComp1), INTENT(INOUT) :: gcCH4    ! Grid Component
   TYPE(ESMF_State), INTENT(INOUT)  :: impChem  ! Import State
   TYPE(ESMF_State), INTENT(INOUT)  :: expChem  ! Export State
   INTEGER, INTENT(OUT) ::  rc                  ! Error return code:
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

   integer n_CH4, i_CH4, j_CH4

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
