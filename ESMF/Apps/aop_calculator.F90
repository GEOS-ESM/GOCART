#  include "MAPL_Generic.h"

      program ext_calculator

!-----------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1   !
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: ext_calculator --- Extinction calculator
!
! !INTERFACE:
!
!      Usage:  ext_calculator.x 
!
! !USES:
!
      use  ESMF
      use  MAPL_Mod
      use  GOCART2G_SimpleBundleMod
      use  GOCART2G_MieMod
      use  GOCART2G_AodMod

      use  m_die
 
      implicit NONE

! !DESCRIPTION: This is a parallel version of the 3D AOD Calculator.
!
! !REVISION HISTORY:
!
!  18Jun2011  da Silva  Derived from mpana_aod.x
!
!EOP
!-----------------------------------------------------------------------

      character(len=*), parameter :: myname = 'ext_calculator'

!     Local variables
!     ---------------
      integer :: rc
      integer :: nymd=0, nhms=0
      integer :: yy, mm, dd, h, m, s

      logical :: verbose = .TRUE.

!     Control variables and obervations
!     ---------------------------------
      type (MAPL_SimpleBundle)   :: q_f              ! aerosol mixing ratio 
      type (MAPL_SimpleBundle)   :: y_f              ! extinction parameters

      integer                          :: n_species    ! number of species
      type (GOCART2G_Mie), pointer     :: MieTables(:) ! (n_species) Mie Tables, etc


      integer                          :: n_tracers  ! number of tracers
      type (GOCART2G_Mie), pointer     :: Mie(:)     ! (n_tracers) Mie Tables, etc

!     Basic ESMF objects
!     ------------------
      type(ESMF_Config)    :: CF       ! Resource file
      type(ESMF_Grid)      :: etaGrid  ! Eta Grid (lon, lat, eta)
      type(ESMF_Time)      :: Time
      type(ESMF_VM)        :: VM

      integer :: Nx, Ny                        ! Layout
      integer :: IM_World, JM_World, LM_WORLD  ! Global Grid dimensions

      call Main()

CONTAINS

!...............................................................................................

     Subroutine Main()

                                   __Iam__('ext_calculator')

!   Initialize the ESMF. For performance reasons, it is important
!    to turn OFF ESMF's automatic logging feature
!   -------------------------------------------------------------
    call ESMF_Initialize (LogKindFlag=ESMF_LOGKIND_NONE, VM=VM, __RC__)
    call ESMF_CalendarSetDefault ( ESMF_CALKIND_GREGORIAN, __RC__ )

    if ( MAPL_am_I_root() ) then
       print *
       print *, '     --------------------------------------'
       print *, '             3D Extinction Calculator'
       print *, '     --------------------------------------'
       print *
    end if

!                                     -------------------
!                                       ESMF Grid, Etc
!                                     -------------------

!   Load resources
!   --------------
    CF = ESMF_ConfigCreate(__RC__)
    call ESMF_ConfigLoadFile(CF, fileName='ext.rc', __RC__)

!   World grid dimensions and layout
!   --------------------------------
    call ESMF_ConfigGetAttribute(CF, IM_World, Label='IM_World:',  __RC__ )
    call ESMF_ConfigGetAttribute(CF, JM_World, Label='JM_World:',  __RC__ )
    call ESMF_ConfigGetAttribute(CF, LM_World, Label='LM_World:',  __RC__ )
    call ESMF_ConfigGetAttribute(CF, Nx,       Label='Layout_Nx:', __RC__ )
    call ESMF_ConfigGetAttribute(CF, Ny,       Label='Layout_Ny:', __RC__ )
    call ESMF_ConfigGetAttribute(CF, verbose,  Label='verbose:',   __RC__ )

!   Create global lat/lon grid
!   --------------------------
    etaGrid = MAPL_LatLonGridCreate (Name='etaGrid',        &
                                     Nx = Nx, Ny = Ny,      &
                                     IM_World = IM_World,   &
                                     JM_World = JM_World,   &
                                     LM_World = LM_World,   &
                                    __RC__ )

!   Validate grid
!   -------------
    call ESMF_GridValidate(etaGrid,__RC__)

!   Get date/time from CF
!   ---------------------
    call ESMF_ConfigGetAttribute(CF, nymd, Label='nymd:', __RC__ )
    call ESMF_ConfigGetAttribute(CF, nhms, Label='nhms:', __RC__ )

!   Create ESMF Time
!   ----------------
    yy = nymd/10000; mm = (nymd-yy*10000) / 100; dd = nymd - (10000*yy + mm*100)
    h  = nhms/10000;  m = (nhms - h*10000) / 100;  s = nhms - (10000*h  +  m*100)
    call ESMF_TimeSet(Time, yy=yy, mm=mm, dd=dd,  h=h,  m=m, s=s)

!                                     -------------------
!                                     Gridded Background
!                                     -------------------


!     Read aerosol mixing ratio
!     -------------------------
      q_f = GOCART2G_SimpleBundleRead (CF, 'aer_filename', etaGrid, &
                                       time=Time, verbose=verbose, __RC__ )

!     Associate mixing ratio tracers with corresponding MieTable
!     ----------------------------------------------------------
      n_tracers = q_f%n3d
      allocate(Mie(n_tracers),__STAT__)
      do i = 1, n_tracers
         j = getTable__(CF, species, q_f%r3d%name)
         Mie(i) => MieTables(j)
      end do
      
!     Load Mie tables
!     ---------------
      n_species = ESMF_ConfigGetLen(CF,'Species:', __RC__)
      - alocate and read mie tables

!     Create SimpleBundle for output fields
!     -------------------------------------      
      y_f = GOCART2G_SimpleBundleCreate ('ext', CF, 'aop_variables', etaGrid, __RC__ )
      
!     Perform Mie calculation
!     -----------------------
      call GOCART2G_AopCalculator3D (y_f, q_f, Mie, verbose, __RC__)

      call MAPL_SimpleBundlePrint(q_f)
      call MAPL_SimpleBundlePrint(y_f)

!     Write file with AOD/Extinction output
!     ------------------------------------
      call GOCART2G_SimpleBundleWrite (y_f, CF, 'ext_filename', Time, __RC__ )

!     All done
!     --------
      call ESMF_Finalize ( __RC__ )

     end subroutine Main


     integer function getTable__(CF, species, name)
     .... write code ....
     end function getTable__
       
    end program ext_calculator



