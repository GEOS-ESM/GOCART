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
      use  MAPL
      use  MAPL_LatLonGridFactoryMod 
      use  GOCART2G_SimpleBundleMod
      use  GOCART2G_MieMod
      use  GOCART2G_AopMod
      use  mpi
      implicit NONE

! !DESCRIPTION: This is a parallel version of the 3D AOD Calculator.
!
! !REVISION HISTORY:
!
!  18Jun2011  da Silva  Derived from mpana_aod.x
!  Jun2022    da Silva/Buchard  Updated for GOCART2G
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
      type (MAPL_SimpleBundle)   :: q_f              ! input: aerosol mixing ratio 
      type (ESMF_FieldBundle)   :: y_f              ! output: extinction parameters

      integer                           :: n_species    ! number of species
      integer                           :: n_tracers  ! number of tracers
      
      type (MieTable), allocatable       :: Mie(:)       ! (n_tracers) Mie Tables, etc

      type (GOCART2G_Mie), allocatable, target    :: MieTables(:)     ! (n_species) Mie Tables, etc

!     Basic ESMF objects
!     ------------------
      type(ESMF_Config)       :: CF       ! Resource file
      type(ESMF_Grid)         :: etaGrid  ! Eta Grid (lon, lat, eta)
      type(ESMF_Time)         :: Time
      type(ESMF_Clock)        :: Clock
      type(ESMF_timeinterval) :: Tinterval
      type(ESMF_VM)           :: VM
      type(ServerManager) ::io_server
      integer :: Nx, Ny                        ! Layout
      integer :: IM_World, JM_World, LM_WORLD  ! Global Grid dimensions
      integer :: i, j, n_wav, iTable, nbins
      integer ::  iBin 
      real, allocatable :: wavelengths(:)     
      integer, pointer :: bin_number(:)
      character(len=255), allocatable :: species_name(:)
      character(len=255), allocatable :: tmp(:)
      character(len=255) :: MieFile
 
      call Main()

CONTAINS

!...............................................................................................

     Subroutine Main()

                                   __Iam__('ext_calculator')

!   Initialize the ESMF. For performance reasons, it is important
!    to turn OFF ESMF's automatic logging feature
!   -------------------------------------------------------------
    call ESMF_Initialize (LogKindFlag=ESMF_LOGKIND_MULTI, VM=VM, __RC__)
    call ESMF_CalendarSetDefault ( ESMF_CALKIND_GREGORIAN, __RC__ )
    call MAPL_Initialize(__RC__) 
    call io_server%initialize(mpi_comm_world)
       
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
    call ESMF_ConfigLoadFile(CF, fileName='aop_calculator.rc', __RC__)

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

    etaGrid = grid_manager%make_grid(                                                 &
                   LatLonGridFactory(im_world=IM_World, jm_world=JM_World, lm=LM_World, &
                   nx=NX, ny=NY, pole='PC', dateline= 'DC', rc = status))
    VERIFY_(status)


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
    call ESMF_TimeIntervalSet(Tinterval, h=1)
    clock = ESMF_ClockCreate (timeStep=Tinterval, &
                                startTime=Time, __RC__ )

!                                     -------------------
!                                     Gridded Background
!                                     -------------------


!     Read aerosol mixing ratio
!     -------------------------
      q_f = GOCART2G_SimpleBundleRead (CF, 'aer_filename', etaGrid, &
                                       time=Time, verbose=verbose, __RC__ )
      
!     Load Mie tables
!     ---------------
      n_wav = ESMF_ConfigGetLen(CF,Label='wavelengths_in_nm:',__RC__)
      allocate(wavelengths(n_wav),__STAT__)
      call ESMF_ConfigGetAttribute(CF, wavelengths, Label='wavelengths_in_nm:', __RC__)  


      n_species = ESMF_ConfigGetLen(CF, Label='Species:', __RC__)
      allocate(species_name(n_species), __STAT__)
      call ESMF_ConfigGetAttribute(CF, species_name, Label='Species:', __RC__)

      allocate(MieTables(n_species), __STAT__)
      do i = 1, n_species
          nbins = ESMF_ConfigGetLen(CF, Label=trim(species_name(i))//':',__RC__)
          allocate(tmp(nbins), __STAT__)  ! note: arg(1) is miefile then bin name
          call ESMF_ConfigGetAttribute(CF, tmp, Label=trim(species_name(i))//':', __RC__)
          MieFile = tmp(1)
          MieTables(i) = GOCART2G_Mie(MieFile, wavelengths*1e-9, __RC__)
          deallocate(tmp)
      end do
      
!     Associate mixing ratio tracers with corresponding MieTable
!
!     ----------------------------------------------------------
      n_tracers = q_f%n3d
      allocate(Mie(n_tracers),__STAT__)
      allocate(bin_number(n_tracers), __STAT__)
      do i = 1, n_tracers

         call getTable__(CF,  q_f%r3(i)%name, iTable, iBin)

         if ( iTable>0 ) then
            Mie(i)%Table => MieTables(iTable)
            Mie(i)%bin_number = iBin
         else
            Mie(i)%Table => null()
         end if
      end do
      
!     Create SimpleBundle for output fields!     -------------------------------------      
!     -------------------------------------      
!      y_f = GOCART2G_SimpleBundleCreate ('ext', CF, 'aop_variables', etaGrid, __RC__ )
      y_f = GOCART2G_ESMFBundleCreate ('ext', CF, 'aop_variables', etaGrid, __RC__ )
      
!     Perform Mie calculation
!     -----------------------
      call GOCART2G_AopCalculator3D (y_f, q_f, Mie, verbose, __RC__)

      call MAPL_SimpleBundlePrint(q_f)
!      call MAPL_SimpleBundlePrint(y_f)

!     Write file with AOP output
!     ---------------------------
!      call GOCART2G_SimpleBundleWrite (y_f, CF, 'ext_filename', Time, __RC__ )
      call GOCART2G_ESMFBundleWrite (y_f, CF, 'ext_filename', wavelengths, Clock, __RC__ )

!     All done
!     --------
      call io_server%finalize()
      call MAPL_Finalize( __RC__ )
      call ESMF_Finalize( __RC__ )
      

     end subroutine Main

     subroutine getTable__(CF,  qname,  ind_table, bin_number)

        implicit NONE

        type (ESMF_config), intent(inout) :: CF
        character(len=*), intent(in)      :: qname               
        integer,            intent(out)   :: ind_table
        integer,            intent(out)   :: bin_number

        integer :: i, j, n_species, nbins, STATUS
        character(len=255), allocatable :: species_name(:)
        character(len=255), allocatable :: bins_name(:)

        n_species = ESMF_ConfigGetLen(CF, Label='Species:', __RC__)
        allocate(species_name(n_species), __STAT__)
        call ESMF_ConfigGetAttribute(CF, species_name, Label='Species:', __RC__)


          ind_table = 0
          bin_number = 0
          do i = 1, n_species
          nbins = ESMF_ConfigGetLen(CF, Label=trim(species_name(i))//':',__RC__)
          allocate(bins_name(nbins), __STAT__)  ! note: arg(1) is miefile then bin name
          call ESMF_ConfigGetAttribute(CF, bins_name, Label=trim(species_name(i))//':', __RC__)
          do j = 2, nbins   ! first index is Miefile
              if(ESMF_UtilStringUpperCase(trim(qname)) == ESMF_UtilStringUpperCase(trim(bins_name(j)))) then
                 ind_table = i
                 bin_number = j - 1
              endif
                            
           enddo
           deallocate(bins_name)
        enddo
     
     end subroutine  getTable__
end program ext_calculator 



