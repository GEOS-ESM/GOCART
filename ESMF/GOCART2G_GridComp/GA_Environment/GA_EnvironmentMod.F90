#include "MAPL_Generic.h"

module GA_EnvironmentMod

   use ESMF
   use MAPL
   use GOCART2G_MieMod

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
       integer                :: rhFlag
       integer                :: nbins
       integer                :: km             ! vertical grid dimension
       real                   :: CDT            ! chemistry timestep (secs)
       integer                :: instance       ! data or computational instance
       real                   :: plid           ! pressure lid [hPa]
       integer                :: klid           ! vertical index of pressure lid
       real, allocatable      :: wavelengths_profile(:) ! wavelengths for profile aop [nm]
       real, allocatable      :: wavelengths_vertint(:) ! wavelengths for vertically integrated aop [nm]
    contains
       procedure :: load_from_config
    end type GA_Environment


    !LOCALS
     integer :: status
     integer :: nbins
     integer :: n_wavelengths_profile, n_wavelengths_vertint, n_channels

 contains


    subroutine load_from_config(self, cfg, universal_cfg, rc)
       class(GA_Environment), intent(inout) :: self
       type(ESMF_Config), intent(inout) :: cfg
       type(ESMF_Config), intent(inout) :: universal_cfg
       integer, optional, intent(out) :: rc

       !   Get nbins from cfg
       call ESMF_ConfigGetAttribute (cfg, self%nbins, label='nbins:', __RC__)
       nbins = self%nbins

       n_wavelengths_profile = ESMF_ConfigGetLen (universal_cfg, label='wavelengths_for_profile_aop_in_nm:', __RC__)
       n_wavelengths_vertint = ESMF_ConfigGetLen (universal_cfg, label='wavelengths_for_vertically_integrated_aop_in_nm:', __RC__)

       !   Parse config file into private internal state
       !   ----------------------------------------------
       allocate(self%radius(nbins), self%rhop(nbins), self%fscav(nbins), self%molwght(nbins), &
                self%fnum(nbins), self%wavelengths_profile(n_wavelengths_profile), &
                self%wavelengths_vertint(n_wavelengths_vertint), &
                __STAT__)
       
       call ESMF_ConfigGetAttribute (cfg, self%radius,     label='particle_radius_microns:', __RC__)
       call ESMF_ConfigGetAttribute (cfg, self%rhop,       label='particle_density:', __RC__)
       call ESMF_ConfigGetAttribute (cfg, self%fscav,      label='fscav:', __RC__)
       call ESMF_ConfigGetAttribute (cfg, self%molwght,    label='molecular_weight:', __RC__)
       call ESMF_ConfigGetAttribute (cfg, self%fnum,       label='fnum:', __RC__)
       call ESMF_ConfigGetAttribute (cfg, self%rhFlag,     label='rhFlag:', __RC__)
       call ESMF_ConfigGetAttribute (cfg, self%plid,       label='pressure_lid_in_hPa:', __RC__)
       call ESMF_ConfigGetAttribute (universal_cfg, self%wavelengths_profile, label='wavelengths_for_profile_aop_in_nm:', __RC__)
       call ESMF_ConfigGetAttribute (universal_cfg, self%wavelengths_vertint, &
                                     label='wavelengths_for_vertically_integrated_aop_in_nm:', __RC__)

    end subroutine load_from_config

end module GA_EnvironmentMod
