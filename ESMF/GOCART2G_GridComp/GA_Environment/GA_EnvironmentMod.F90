#include "MAPL_Generic.h"

module GA_EnvironmentMod

   use ESMF
   use MAPL
   use GOCART2G_MieMod
   use mapl3g_generic, only: MAPL_GridCompGetResource

   implicit none
   private

   public :: GA_Environment

   type :: GA_Environment
       type(GOCART2G_Mie)     :: rad_Mie, diag_Mie
       real, allocatable      :: radius(:)      ! particle effective radius [um]
       real, allocatable      :: rhop(:)        ! soil class density [kg m-3]
       real, allocatable      :: fscav(:)       ! scavenging efficiency
!      logical                :: scav_byColdCloud ! new flag example
       real, allocatable      :: molwght(:)     ! molecular weight            !NOT UNIVERSAL ONLY FOR GASES,
       real, allocatable      :: fnum(:)        ! number of particles per kg mass
       real, allocatable      :: fwet_ice(:)    ! large scale wet removal scaling factor for ice
       real, allocatable      :: fwet_snow(:)   ! large scale wet removal scaling factor for snow
       real, allocatable      :: fwet_rain(:)   ! large scale wet removal scaling factor for rain
       real                   :: washout_tuning ! tuning factor for washout process (1 by default)
       real                   :: wet_radius_thr ! wet radius threshold [um]
       integer                :: rhFlag
       integer                :: nbins
       integer                :: km             ! vertical grid dimension
       real                   :: CDT            ! chemistry timestep (secs)
       integer                :: instance       ! data or computational instance
       real                   :: plid           ! pressure lid [hPa]
       integer                :: klid           ! vertical index of pressure lid
       character(:), allocatable :: wet_removal_scheme     ! name of wet removal scheme
       character(:), allocatable :: settling_scheme        ! settling option (1 - use SettlingSolver, 2 - use SettlingSolverUFS)
       real, allocatable      :: wavelengths_profile(:) ! wavelengths for profile aop [nm]
       real, allocatable      :: wavelengths_vertint(:) ! wavelengths for vertically integrated aop [nm]
    contains
       procedure :: load_from_config
    end type GA_Environment

 contains

    subroutine load_from_config(self, gc, rc)
       class(GA_Environment), intent(inout) :: self
       type(ESMF_GridComp), intent(inout) :: gc
       integer, optional, intent(out) :: rc

       ! Local variables
       character(:), allocatable :: wet_removal_scheme
       character(:), allocatable :: settling_scheme
       real, allocatable :: ones(:)
       integer :: nbins, status

       call MAPL_GridCompGetResource(gc, "nbins", self%nbins, _RC)
       nbins = self%nbins

       call MAPL_GridCompGetResource(gc, "particle_radius_microns", self%radius, _RC)
       call MAPL_GridCompGetResource(gc, "particle_density", self%rhop, _RC)
       call MAPL_GridCompGetResource(gc, "fscav", self%fscav, _RC)
       call MAPL_GridCompGetResource(gc, "molecular_weight", self%molwght, _RC)
       call MAPL_GridCompGetResource(gc, "fnum", self%fnum, _RC)
       call MAPL_GridCompGetResource(gc, "pressure_lid_in_hPa", self%plid, _RC)

       call MAPL_GridCompGetResource(gc, "wet_radius_thr", self%wet_radius_thr, default=0.05, _RC)
       call MAPL_GridCompGetResource(gc, "washout_tuning", self%washout_tuning, default=1.0, _RC)
       call MAPL_GridCompGetResource(gc, "wet_removal_scheme", wet_removal_scheme, default="gocart", _RC)
       call MAPL_GridCompGetResource(gc, "settling_scheme", settling_scheme, default='gocart', _RC)
       self%wet_removal_scheme = ESMF_UtilStringLowerCase(wet_removal_scheme, _RC)
       self%settling_scheme = ESMF_UtilStringLowerCase(settling_scheme, _RC)

       allocate(ones(nbins), source=1.0)
       call MAPL_GridCompGetResource(gc, "fwet_ice", self%fwet_ice, default=ones, _RC)
       call MAPL_GridCompGetResource(gc, "fwet_snow", self%fwet_snow, default=ones, _RC)
       call MAPL_GridCompGetResource(gc, "fwet_rain", self%fwet_rain, default=ones, _RC)

       call MAPL_GridCompGetResource(gc, "wavelengths_for_profile_aop_in_nm", self%wavelengths_profile, _RC)
       call MAPL_GridCompGetResource(gc, "wavelengths_for_vertically_integrated_aop_in_nm", self%wavelengths_vertint, _RC)

       call validate_config(self, _RC)

    end subroutine load_from_config

    subroutine validate_config(self, rc)
       class(GA_Environment), intent(in) :: self
       integer, optional, intent(out) :: rc

       !   Local variables
       character(len=255) :: msg

       !   Validate config
       !   ---------------
       !
       !   * Rainout efficiency
       write(msg,'(5(2x,g20.8))') self%fwet_ice
       call ESMF_LogWrite("GA: config: fwet_ice: "//msg)
       _ASSERT_RC(any(abs(self%fwet_ice - 1.) <= 1.), "Error. Rainout efficiency (fwet) must be between 0. and 1.", ESMF_RC_VAL_OUTOFRANGE)
       write(msg,'(5(2x,g20.8))') self%fwet_snow
       call ESMF_LogWrite("GA: config: fwet_snow: "//msg)
       _ASSERT_RC(any(abs(self%fwet_snow - 1.) <= 1.), "Error. Rainout efficiency (fwet) must be between 0. and 1.", ESMF_RC_VAL_OUTOFRANGE)
       write(msg,'(5(2x,g20.8))') self%fwet_rain
       call ESMF_LogWrite("GA: config: fwet_rain: "//msg)
       _ASSERT_RC(any(abs(self%fwet_rain - 1.) <= 1.), "Error. Rainout efficiency (fwet) must be between 0. and 1.", ESMF_RC_VAL_OUTOFRANGE)

       !   * Wet removal scheme
       _ASSERT_RC(any(self%wet_removal_scheme == [character(len=7) :: 'gocart','ufs']), "Error. Unallowed wet removal scheme: "//trim(self%wet_removal_scheme)//". Allowed: gocart, ufs", ESMF_RC_NOT_IMPL)
       _ASSERT_RC(any(self%settling_scheme == [character(len=7) :: 'gocart','ufs']), "Error. Unallowed settling option: "//trim(self%settling_scheme)//". Allowed: gocart, ufs", ESMF_RC_NOT_IMPL)

    end subroutine validate_config

end module GA_EnvironmentMod
