#  include "MAPL_Generic.h"

      program reff_calculator

!-----------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1   !
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: ana_aod:  AOD Analysis Application
!
! !INTERFACE:
!
!      Usage:  reff_calculator.x 
!
! !USES:
!
      use  ESMF
      use  MAPL_Mod
      use  Chem_SimpleBundleMod
      use  Chem_RegistryMod
      use  Chem_MieMod
      use  Chem_AodMod

      use  m_die
      use  m_FileResolv
      use  m_StrTemplate
      use  m_inpak90
 
      implicit NONE

! !DESCRIPTION: Reads a chem bundle file and writes a similar with the effective radius.
!               This is a MPI application.
!
! !REVISION HISTORY:
!
!  13jan2011  da Silva  Derived from mpana_aod.
!
!EOP
!-----------------------------------------------------------------------

      character(len=*), parameter :: myname = 'reff_calculator'

!     Local variables
!     ---------------
      integer :: rc
      integer :: nymd=0, nhms=0
      integer :: yy, mm, dd, h, m, s

      logical :: verbose = .TRUE.

!     Control variables and obervations
!     ---------------------------------
      type (MAPL_SimpleBundle)   :: q_f          ! Gridded concentration background 

      type (Chem_Mie)      :: Mie                ! Mie Tables, etc
      type (Chem_Registry) :: aerReg             ! Registry with many species

!     Basic ESMF objects
!     ------------------
      type(ESMF_Config)    :: CF       ! Resource file
      type(ESMF_Grid)      :: etaGrid  ! Eta Grid (lon, lat, eta)
      type(ESMF_Time)      :: Time
      type(ESMF_VM)        :: VM

      integer :: Nx, Ny         ! Layout
      integer :: IM_World, JM_World, LM_WORLD ! Global Grid dimensions

      call Main()

CONTAINS

!...............................................................................................
      Subroutine Main()

                                   __Iam__('reff_calculator')

!   Initialize the ESMF. For performance reasons, it is important
!    to turn OFF ESMF's automatic logging feature
!   -------------------------------------------------------------
    call ESMF_Initialize (logKindFlag=ESMF_LOGKIND_NONE, VM=VM, __RC__)
    call ESMF_CalendarSetDefault ( ESMF_CALKIND_GREGORIAN, __RC__ )

    if ( MAPL_am_I_root() ) then
       print *
       print *, '     --------------------------------------'
       print *, '     ana_aod - the aod analysis application'
       print *, '     --------------------------------------'
       print *
    end if

!                                     -------------------
!                                       ESMF Grid, Etc
!                                     -------------------

!   Load resources
!   --------------
    CF = ESMF_ConfigCreate(__RC__)
    call ESMF_ConfigLoadFile(CF, fileName='reff.rc', __RC__)

!   World grid dimensions and layout
!   --------------------------------
    call ESMF_ConfigGetAttribute(CF, IM_World, Label='IM_World:',  __RC__ )
    call ESMF_ConfigGetAttribute(CF, JM_World, Label='JM_World:',  __RC__ )
    call ESMF_ConfigGetAttribute(CF, LM_World, Label='LM_World:',  __RC__ )
    call ESMF_ConfigGetAttribute(CF, Nx,       Label='Layout_Nx:', __RC__ )
    call ESMF_ConfigGetAttribute(CF, Ny,       Label='Layout_Ny:', __RC__ )
    call ESMF_ConfigGetAttribute(CF, verbose,  Label='verbose:',   __RC__ )

!   Create global grids:
!   -------------------
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

!     Registry
!     --------
      aerReg = Chem_RegistryCreate ( rc, 'Chem_AerRegistry.rc' )
      if ( rc == 0 ) then
         if ( MAPL_AM_I_ROOT() ) then
            call Chem_RegistryPrint(aerReg)
         end if
      else
         call die(myname,'could not read Chem Registry for specifies')
      end if

!     Read aerosol mixing ratio
!     -------------------------
      q_f = Chem_SimpleBundleRead (CF, 'aer_bkg_filename', etaGrid, &
                                   time=Time, verbose=verbose, __RC__ )
      Mie = Chem_MieCreate(CF, chemReg=aerReg, __RC__)

!     Calcular rEff
!     -------------
      call MAPL_SimpleBundlePrint(q_f)
      call Chem_rEffCalculator (q_f, q_f, Mie, verbose, __RC__)
      call MAPL_SimpleBundlePrint(q_f)

!     Write rEff file
!     ---------------
      call Chem_SimpleBundleWrite (q_f, CF, 'aer_reff_filename', Time, __RC__ )

!     All done
!     --------
      call ESMF_Finalize ( __RC__ )

     end subroutine Main

    end program reff_calculator



