#include "MAPL_Generic.h"

module GA_EnvironmentMod

   use ESMF
   use MAPL
   use GOCART2G_MieMod
   use mapl3g_HConfig_get, only: HConfigParams, MAPL_HConfigGet
   use pflogger, only: logger_t => logger

   implicit none
   private

   public :: GA_Environment

   type :: GA_Environment
       type(GOCART2G_Mie)     :: rad_Mie, diag_Mie
       real, allocatable      :: radius(:)      ! particle effective radius [um]
       real, allocatable      :: rhop(:)        ! soil class density [kg m-3]
       real, allocatable      :: fscav(:)       ! scavenging efficiency
       ! logical              :: scav_byColdCloud ! new flag example
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
       real, allocatable      :: wavelengths_profile(:) ! wavelengths for profile aop [nm]
       real, allocatable      :: wavelengths_vertint(:) ! wavelengths for vertically integrated aop [nm]
    contains
       procedure :: load_from_config
    end type GA_Environment

 contains

    subroutine load_from_config(self, hconfig, logger, rc)
       class(GA_Environment), intent(inout) :: self
       type(ESMF_HConfig), intent(in) :: hconfig
       class(logger_t), pointer, intent(in) :: logger
       integer, optional, intent(out) :: rc

       !   Local variables
       character(:), allocatable :: wet_removal_scheme
       type(HConfigParams) :: params
       real, allocatable :: ones(:)
       integer :: nbins, status

       !   Get nbins from cfg
       self%nbins = ESMF_HConfigAsI4(hconfig, keyString="nbins", _RC)
       ! call ESMF_ConfigGetAttribute (cfg, self%nbins, label='nbins:', __RC__)
       nbins = self%nbins

       ! n_wavelengths_profile = ESMF_ConfigGetLen (universal_cfg, label='wavelengths_for_profile_aop_in_nm:', __RC__)
       ! n_wavelengths_vertint = ESMF_ConfigGetLen (universal_cfg, label='wavelengths_for_vertically_integrated_aop_in_nm:', __RC__)

       !   Parse config file into private internal state
       !   ----------------------------------------------
       ! allocate(self%radius(nbins), self%rhop(nbins), self%fscav(nbins), self%molwght(nbins), &
       !          self%fnum(nbins), self%fwet_ice(nbins), self%fwet_snow(nbins), self%fwet_rain(nbins), &
       !          self%wavelengths_profile(n_wavelengths_profile), &
       !          self%wavelengths_vertint(n_wavelengths_vertint), &
       !          __STAT__)

       self%radius = ESMF_HConfigAsR4Seq(hconfig, keyString="particle_radius_microns", _RC)
       self%rhop = ESMF_HConfigAsR4Seq(hconfig, keyString="particle_density", _RC)
       self%fscav = ESMF_HConfigAsR4Seq(hconfig, keyString="fscav", _RC)
       self%molwght = ESMF_HConfigAsR4Seq(hconfig, keyString="molecular_weight", _RC)
       self%fnum = ESMF_HConfigAsR4Seq(hconfig, keyString="fnum", _RC)
       self%plid = ESMF_HConfigAsR4(hconfig, keyString="pressure_lid_in_hPa", _RC)
       ! call ESMF_ConfigGetAttribute (cfg, self%radius,     label='particle_radius_microns:', __RC__)
       ! call ESMF_ConfigGetAttribute (cfg, self%rhop,       label='particle_density:', __RC__)
       ! call ESMF_ConfigGetAttribute (cfg, self%fscav,      label='fscav:', __RC__)
       ! call ESMF_ConfigGetAttribute (cfg, self%molwght,    label='molecular_weight:', __RC__)
       ! call ESMF_ConfigGetAttribute (cfg, self%fnum,       label='fnum:', __RC__)
       ! call ESMF_ConfigGetAttribute (cfg, self%plid,       label='pressure_lid_in_hPa:', __RC__)
       
       params = HConfigParams(hconfig, "wet_radius_thr", logger=logger)
       call MAPL_HConfigGet(params, self%wet_radius_thr, 0.05, _RC)
       params = HConfigParams(hconfig, "washout_tuning", logger=logger)
       call MAPL_HConfigGet(params, self%washout_tuning, 1.0, _RC)
       params = HConfigParams(hconfig, "wet_removal_scheme", logger=logger)
       call MAPL_HConfigGet(params, wet_removal_scheme, "gocart", _RC)
       self%wet_removal_scheme = ESMF_UtilStringLowerCase(wet_removal_scheme, _RC)
       ! call ESMF_ConfigGetAttribute (cfg, self%wet_radius_thr,  label='wet_radius_thr:', default=0.05, __RC__)
       ! call ESMF_ConfigGetAttribute (cfg, self%washout_tuning,  label='washout_tuning:', default=1.0, __RC__)
       ! call ESMF_ConfigGetAttribute (cfg, wet_removal_scheme, label='wet_removal_scheme:', default='gocart', __RC__)
       ! self%wet_removal_scheme = ESMF_UtilStringLowerCase(trim(wet_removal_scheme), __RC__)

       allocate(ones(nbins), source=1.0)
       params = HConfigParams(hconfig, "fwet_ice", logger=logger)
       call MAPL_HConfigGet(params, self%fwet_ice, ones, _RC)
       params = HConfigParams(hconfig, "fwet_snow", logger=logger)
       call MAPL_HConfigGet(params, self%fwet_snow, ones, _RC)
       params = HConfigParams(hconfig, "fwet_rain", logger=logger)
       call MAPL_HConfigGet(params, self%fwet_rain, ones, _RC)
       ! call ESMF_ConfigGetAttribute (cfg, self%fwet_ice,   label='fwet_ice:', default=1.0, __RC__)
       ! call ESMF_ConfigGetAttribute (cfg, self%fwet_snow,  label='fwet_snow:', default=1.0, __RC__)
       ! call ESMF_ConfigGetAttribute (cfg, self%fwet_rain,  label='fwet_rain:', default=1.0, __RC__)

       self%wavelengths_profile = ESMF_HConfigAsR4Seq(hconfig, keyString="wavelengths_for_profile_aop_in_nm", _RC)
       self%wavelengths_vertint = ESMF_HConfigAsR4Seq(hconfig, keyString="wavelengths_for_vertically_integrated_aop_in_nm", _RC)

       call validate_config(self, __RC__)

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

    end subroutine validate_config

end module GA_EnvironmentMod
