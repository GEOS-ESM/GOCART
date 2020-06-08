
!=============================================================================
!
! $Id$
!
! CODE DEVELOPER
!   John Tannahill, LLNL
!   jrt@llnl.gov
!
! FILE
!   gmi_emiss_constants.h
!
! Last modified: 03/29/04
!     -  incorporated Bryan Duncan's suggested values.
! DESCRIPTION
!   This include file contains constants for the Harvard biogenic and soil
!   emissions.
! HISTORY
!   - August 12, 2005 * Jules Kouatchou
!     Remove "NLANDHAR" as a parameter in this file and add it as a namelist
!     variable. This is done because we want the user to determine its value
!     at run time. If not, the user would have to modify its value at
!     compilation whenever there is a change of model resolution.
!
!   Mar 30, 2017: Moved this file from GmiInclude/ to Chem_Shared/ for TR
!
!=============================================================================


!c!?  Should be fixed sometime.
!     -------------------------------------------------------------
!     NOTE THAT NEED TO MATCH THE NLANDHAR PARAMETER BELOW WITH THE
!     RESOLUTION YOU ARE RUNNING AT!
!     -------------------------------------------------------------

!!      integer, parameter ::
!!     &  NLANDHAR = 3920      ! number of grids on land for 2x2.5 resolution
!!     &  NLANDHAR = 1118      ! number of grids on land for 4x5   resolution
!        integer :: NLANDHAR

!        common /emiss_NLANDHAR/ NLANDHAR


      integer, parameter ::  &
     &  NPOLY    =   20,       & ! number of coefficients for polynomial fits
     &  NPULSE   =    3,       & ! number of types of pulsing
     &  NSOIL    =   11,       & ! defined soil types
     &  NTYPE    =   15,       & ! maximum number of vegetation types in any
                             ! grid box
     &  NVEGTYPE =   74      ! maximum number of surface types (Olson)


!     ----------------------------------------------------------
!     MWTCARBON : molecular weight of Carbon (g/mol)
!     MWTCO     : molecular weight of CO     (g/mol)
!     MWTMONOT  : molecular weight of monoterpene (C10H16) (g/mol)
!     ----------------------------------------------------------

      real*8, parameter ::  &
     &  MWTCARBON = 12.01d0,  &
     &  MWTCO     = 28.01d0,  &
     &  MWTMONOT  = 136.00d0


!     ---------------------------------------------
!     ATOMSC_PER_MOLECCO    : atomsC/molec CO
!     ATOMSC_PER_MOLECISOP  : atomsC/molec isoprene
!     ATOMSC_PER_MOLECMONOT : atomsC/molec monoterpene
!     ATOMSC_PER_MOLECPRPE  : atomsC/molec proprene
!     ---------------------------------------------

      real*8, parameter ::  &
     &  ATOMSC_PER_MOLECCO    = 1.0d0,  &
     &  ATOMSC_PER_MOLECISOP  = 5.0d0,  &
     &  ATOMSC_PER_MOLECMONOT = 10.0d0,  &
     &  ATOMSC_PER_MOLECPRPE  = 3.0d0


!     ----------------------------------------------------------------------
!     BIOSCAL : factor for biogenic source of propene, scaled from isoprene
!               (atomsCalkenes / atomsCisoprene)
!
!     We need to scale the isoprene flux to get the biogenic emissions of
!     alkenes (probably OK for summertime conditions).
!
!     The scaling factor is based on work by Allen Goldstein.  His values
!     indicate emission ratios of ethene:propene:butene=4:2:1 (on a per
!     molecule basis), with total emissions approx. equal to 10% of isoprene
!     emissions (again, on a per molecule basis):
!
!       (10 molec alkenes  / 100 molec  isoprene) *
!       ( 1 molec isoprene /   5 atomsC isoprene) *
!       ( 3      molec butene + propene   / 7 molec total alkenes) *
!       ( 3.3333 atomsCbutene+propene mix / 1 molec butene+propene mix) =
!
!         0.0286 atomsCbutene+propene / atomsCisoprene
!
!     Note that 3.3333 atomsC/molec is the weighted average for this mix.
!     Note also that this factor now excludes ethene.
!     ----------------------------------------------------------------------

      real*8, parameter ::  &
     &  BIOSCAL = 0.0286d0


!     ------------------------------------------------------------------
!     ICO_FAC_ISOP : factor for biogenic source of CO from methanol
!                    oxidation, scaled from isoprene
!
!     We need to scale the isoprene flux to get the CH3OH (methanol)
!     flux.  Currently, the annual isoprene flux in GEOS-CHEM is
!     ~397 tg C.
!
!     Daniel Jacob recommends a flux of 100 tg/yr CO from CH3OH
!     oxidation based on Singh et al., 2000 (JGR 105, 3795-3805), who
!     estimate a global methanol source of 122 tg/yr, of which most
!     (75 tg/yr) is "primary biogenic".  He also recommends for now,
!     that the CO flux from CH3OH oxidation be scaled from monthly mean
!     isoprene flux.
!
!     To get CO from methanol oxidation, we must therefore multiply
!     the isoprene flux by the following scale factor:
!
!       (100 tg CO from CH3OH oxidation / 380 tg C from isoprene flux) *
!       ( 60 g C/mole isoprene          /  68 g isoprene/mole)
!     ------------------------------------------------------------------

      real*8, parameter ::  &
     &  ICO_FAC_ISOP = (100.0d0 / 380.0d0) * (60.0d0 / 68.0d0)


!     ------------------------------------------------------------------
!     ICO_FAC_MONOT : factor for biogenic source of CO from monoterpene
!                     oxidation.
!
!     Assume the production of CO from monoterpenes is instantaneous
!     even though the lifetime of intermediate species may be on the
!     order of hours or days.  This assumption will likely cause CO
!     from monoterpene oxidation to be too high in the box in which
!     the monoterpene is emitted.
!
!     The CO yield here is taken from:
!
!     Hatakeyama et al., JGR, Vol. 96, p. 947-958 (1991) =>
!       "The ultimate yield of CO from the tropospheric oxidation of
!       terpenes (including both O3 and OH reactions) was estimated to
!       be 20% on the carbon number basis."  (They studied alpha- &
!       beta-pinene.)
!
!     Vinckier et al., Fresenius Env. Bull., Vol. 7, p.361-368 (1998) =>
!       "R(CO) = 1.8+/-0.3"  ((1.8/10) is about 20%.)
!     ------------------------------------------------------------------

      real*8, parameter ::  &
     &  ICO_FAC_MONOT =  &
     &    0.2d0 *  &
     &    (ATOMSC_PER_MOLECMONOT / ATOMSC_PER_MOLECCO) *  &
     &    (MWTCO / MWTMONOT)


      real*8, parameter  ::  &
     &  PULSE_DECAY(NPULSE) =    & ! values from Yienger & Levy
     &    (/ 0.805d0, 0.384d0, 0.208d0 /)

      real*8, parameter  ::  &
     &  PULSE_FAC(NPULSE)   =    & ! initial pulsing factor following a precip.
                               ! event
     &    (/ 5.0d0, 10.0d0, 15.0d0 /)


      real*8, parameter  ::  &
     &  SOIL_AD(NSOIL)      =    & ! dry biome coefficient
     &    (/ 0.0d0,  8.6d0,  0.22d0, 0.40d0,  0.22d0, 1.44d0,  &
     &       2.65d0, 2.65d0, 2.65d0, 0.003d0, 0.37d0 /)

      real*8, parameter  ::  &
     &  SOIL_AW(NSOIL)      =    & ! wet biome coefficient
     &    (/ 0.0d0,  2.6d0,  0.03d0, 0.06d0,  0.03d0, 0.17d0,  &
     &       0.36d0, 0.36d0, 0.36d0, 0.003d0, 0.05d0 /)


      real*8, parameter  ::  &
     &  SOIL_EXT(NSOIL)     =    & ! canopy wind extinction coefficient
     &    (/ 0.1d0, 4.0d0, 4.0d0, 4.0d0, 4.0d0, 2.0d0,  &
     &       1.0d0, 2.0d0, 2.0d0, 0.5d0, 0.1d0 /)


      real*8, parameter  ::  &
     &  SOIL_T1(NSOIL)      =    & ! first  coefficient used to convert from
                               ! surface temperture to soil temperature
     &    (/ 0.0d0,  0.84d0, 0.84d0, 0.84d0, 0.84d0, 0.66d0, 0.66d0,  &
     &       1.03d0, 1.03d0, 0.92d0, 0.66d0 /)

      real*8, parameter  ::  &
     &  SOIL_T2(NSOIL)      =    & ! second coefficient used to convert from
                               ! surface temperture to soil temperature
     &    (/ 0.0d0, 3.6d0, 3.6d0, 3.6d0, 3.6d0, 8.8d0, 8.8d0,  &
     &       2.9d0, 2.9d0, 4.4d0, 8.8d0 /)

