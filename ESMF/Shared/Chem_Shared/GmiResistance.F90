
!=============================================================================
!
! $Id$
!
! CODE DEVELOPER
!   Dan Bergmann, LLNL
!   dbergmann@llnl.gov
!
! FILE
!   GmiResistance.F90
!
! ROUTINES
!   Canopy_Resistance
!   Surface_Resistance
!
! NOTE
!   Mar 30, 2017: Moved and renamed GMI_GridComp/GmiDeposition/resistance.F90
!                 to Chem_Shared/GmiResistance.F90, for use by TR
!
!=============================================================================


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Canopy_Resistance
!
! DESCRIPTION
!   This routine calculates bulk surface resistance of the canopy from the
!   network of resistances in parallel and in series.
!
! ARGUMENTS
!   rdc     : tbd
!   rix     :
!   press   :
!   tempk1  :
!   f01     :
!   hstar1  :
!   xmw1    :
!   rac1    :
!   rclo1   :
!   rcls1   :
!   rgso1   :
!   rgss1   :
!   rlu1    :
!   rsurfc1 :
!
!-----------------------------------------------------------------------------

      subroutine Canopy_Resistance  &
     &  (rdc, rix, press, tempk1, f01, hstar1, xmw1,  &
     &   rac1, rclo1, rcls1, rgso1, rgss1, rlu1, rsurfc1)

      implicit none

#     include "gmi_phys_constants.h"


!     ----------------------
!     Argument declarations.
!     ----------------------

      real*8  :: rdc
      real*8  :: rix
      real*8  :: press
      real*8  :: tempk1
      real*8  :: f01
      real*8  :: hstar1
      real*8  :: xmw1
      real*8  :: rac1
      real*8  :: rclo1
      real*8  :: rcls1
      real*8  :: rgso1
      real*8  :: rgss1
      real*8  :: rlu1
      real*8  :: rsurfc1


!     -----------------------
!     Parameter declarations.
!     -----------------------

      real*8, parameter ::  &
     &  BIG1 = 1.0d25


!     ----------------------
!     Function declarations.
!     ----------------------

      real*8, external  :: Diffg


!     ----------------------
!     Variable declarations.
!     ----------------------

      real*8  :: dtmp1, dtmp2, dtmp3, dtmp4

      real*8  :: rclx
      real*8  :: rgsx
      real*8  :: rixx
      real*8  :: rluxx


!     ----------------
!     Begin execution.
!     ----------------

      rixx =  rix * ( Diffg (tempk1, press, MWTH2O*KGPG) /  &
     &                Diffg (tempk1, press, xmw1  *KGPG) ) +  &
     &        1.0d0 / (hstar1 / 3000.0d0 + 100.0d0 * f01)

      if (rlu1 < 9999.0d0) then
        rluxx = rlu1 / (hstar1 / 1.0d+05 + f01)
      else
        rluxx = BIG1
      end if

!     ---------------------------------------------------------------------
!     To prevent virtually zero resistance to species with huge HSTAR, such
!     as HNO3, a minimum value of rluxx needs to be set.  The rationality
!     of the existence of such a minimum is demonstrated by the observed
!     relationship between Vd(NOy-NOx) and ustar in Munger et. al. [1996];
!     Vd(HNO3) never exceeds 2 cm*s^-1 in observations.  The corresponding
!     minimum resistance is 50 s*m^-1.  This correction was introduced by
!     J.Y. Liang on 7/9/95.
!     ---------------------------------------------------------------------

      rluxx = Max (rluxx, 50.0d0)


      rgsx =  &
     &  1.0d0 /  &
     &  (hstar1 / 1.0d+05 / rgss1 + f01 / rgso1)

      rclx =  &
     &  1.0d0 /  &
     &  (hstar1 / 1.0d+05 / rcls1 + f01 / rclo1)


!     -----------------------------------------------------------------------
!     Get the bulk surface resistance of the canopy, rsurfc, from the network
!     of resistances in parallel and in series (Fig. 1 of Wesely [1989]).
!     -----------------------------------------------------------------------

      dtmp1 = 1.0d0 / rixx
      dtmp2 = 1.0d0 / rluxx
      dtmp3 = 1.0d0 / (rac1 + rgsx)
      dtmp4 = 1.0d0 / (rdc + rclx)

      rsurfc1 = 1.0d0 / (dtmp1 + dtmp2 + dtmp3 + dtmp4)


      return

      end subroutine Canopy_Resistance


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Surface_Resistace
!
! DESCRIPTION
!   This routine tbd
!
! ARGUMENTS
!   idep1   : tbd
!   rac1    :
!   rclo1   :
!   rcls1   :
!   rgso1   :
!   rgss1   :
!   ri1     :
!   rlu1    :
!   rt      :
!   tempc1  :
!   cfrac1  :
!   radiat1 :
!   suncos1 :
!   xlai1   :
!   rix     :
!
!-----------------------------------------------------------------------------

      subroutine Surface_Resistance  &
     &  (idep1, rac1, rclo1, rcls1, rgso1, rgss1, ri1, rlu1,  &
     &   rt, tempc1, cfrac1, radiat1, suncos1, xlai1, rix)

      implicit none

#     include "gmi_emiss_constants.h"
#     include "gmi_drydep_data.h"


!     ----------------------
!     Argument declarations.
!     ----------------------

      integer :: idep1
      real*8  :: rac1
      real*8  :: rclo1
      real*8  :: rcls1
      real*8  :: rgso1
      real*8  :: rgss1
      real*8  :: ri1
      real*8  :: rlu1
      real*8  :: rt
      real*8  :: tempc1
      real*8  :: cfrac1
      real*8  :: radiat1
      real*8  :: suncos1
      real*8  :: xlai1
      real*8  :: rix


!     -----------------------
!     Parameter declarations.
!     -----------------------

      real*8, parameter ::  &
     &  BIG1 = 1.0d6

      real*8, parameter ::  &
     &  BIG2 = 1.0d25


!     ----------------------
!     Function declarations.
!     ----------------------

      real*8, external  :: Biofit


!     ----------------------
!     Variable declarations.
!     ----------------------

      real*8  :: gfaci, gfact

      real*8  :: xdrycoeff(NPOLY)


!     ----------------
!     Begin execution.
!     ----------------

!     ----------------------------------------------------------------------
!     Read the internal resistance ri (minimum stomatal resistance for water
!     vapor, per unit area of leaf) from the IRI array; a "9999" value means
!     no deposition to stomata, so we impose a very large value for ri.
!     ----------------------------------------------------------------------

      ri1 = iri(idep1)

      if (ri1 >= 9999.0d0) ri1 = 1.0d25


!     ----------------------------------------------------------------------
!     Cuticular resistances IRLU read in from drydep table are per unit area
!     of leaf; divide them by the leaf area index to get a cuticular
!     resistance for the bulk canopy.  If IRLU is "9999", it means there are
!     no cuticular surfaces on which to deposit, so we impose a very large
!     value for rlu.
!     ----------------------------------------------------------------------

      if ((irlu(idep1) >= 9999) .or. (xlai1 <= 0.0d0)) then

        rlu1 = 1.0d25

      else

        rlu1 = irlu(idep1)
        rlu1 = (rlu1 / xlai1) + rt

      end if


!     ----------------------------------------------------------
!     The following are the remaining resistances for the Wesely
!     resistance-in-series model for a surface canopy
!     (see Atmos. Environ. paper, Fig.1).
!     ----------------------------------------------------------

      rac1  = Max (irac(idep1), 1)
      if (rac1  >= 9999.0d0) rac1  = BIG2

      rgss1 = irgss(idep1)
      rgss1 = Max (rgss1+rt, 1.0d0)

      rgso1 = irgso(idep1)
      rgso1 = Max (rgso1+rt, 1.0d0)
      if (rgso1 >= 9999.0d0) rgso1 = BIG2

      rcls1 = ircls(idep1)
      rcls1 = rcls1 + rt
      if (rcls1 >= 9999.0d0) rcls1 = BIG2

      rclo1 = irclo(idep1)
      rclo1 = rclo1 + rt
      if (rclo1 >= 9999.0d0) rclo1 = BIG2


!     ----------------------------------------------------------------------
!     Adjust stomatal resistances for insolation and temperature =>
!
!     Temperature adjustment is from Wesely [1989], equation (3).
!
!     Light adjustment by the function Biofit is described by Wang [1996].
!     It combines:
!       - local dependence of stomal resistance on the intensity I of light
!         impinging the leaf; this is expressed as a mutliplicative
!         factor I/(I+b) to the stomatal resistance where b = 50 W*m^-2
!         (equation (7) of Baldocchi et. al. [1987]);
!       - radiative transfer of direct and diffuse radiation in the
!         canopy using equations (12)-(16) from Guenther et. al. [1995];
!       - separate accounting of sunlit and shaded leaves using
!         equation (12) of Guenther et. al. [1995];
!       - partitioning of the radiation at the top of the canopy into direct
!         and diffuse components using a parameterization to results from
!         an atmospheric radiative transfer model [Wang, 1996].
!     The dependent variables of the function Biofit are the leaf area
!     index (xylai), the cosine of zenith angle (suncos) and the fractional
!     cloud cover (cfrac).  The factor gfaci integrates the light
!     dependence over the canopy depth; sp even though ri is input per
!     unit area of leaf, it need not be scaled by lai to yield a bulk
!     canopy value because that is already done in the gfaci formulation.
!     ----------------------------------------------------------------------

      if (ri1 >= 9999.0d0) then

        rix = ri1

      else

        if ((tempc1 > 0.0d0) .and. (tempc1 < 40.0d0)) then
          gfact = 400.0d0 / tempc1 / (40.0d0 - tempc1)
        else
          gfact = 100.0d0
        end if

        if ((radiat1 > 1.0d-5) .and. (xlai1 > 0.0d0)) then
          xdrycoeff(:) = DRYCOEFF(:)
          gfaci = 1.0d0 /  Biofit (cfrac1, suncos1, xlai1, xdrycoeff)
        else
          gfaci = 100.0d0
        end if

        rix = ri1 * gfact * gfaci

      end if

      return

      end subroutine Surface_Resistance

