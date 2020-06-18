
!=============================================================================
!
! $Id$
!
! CODE DEVELOPER
!   John Tannahill, LLNL
!   jrt@llnl.gov
!
!   Original code from:
!     Harvard tropospheric emissions module for 3D applications;
!       by Yuhang Wang, Gerry Gardner, and Prof. Daniel Jacob
!       of Harvard University (Release V1.0)
!
! FILE
!   emiss_utils.F
!
! ROUTINES
!   Biofit
!   Sun_Param
!   Diffg
!
!=============================================================================


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Biofit
!
! DESCRIPTION
!   This routine calculates the light correction.
!
!   Light adjustment by the function Biofit is described by Wang [1996].
!   It combines:
!     * local dependence of stomal resistance on the intensity I of light
!       impinging the leaf; this is expressed as a mutliplicative factor
!       I/(I+b) to the stomatal resistance where b = 50 (W*m^-2)
!       (equation (7) of Baldocchi et al. [1987])
!     * radiative transfer of direct and diffuse radiation in the canopy
!       using equations (12)-(16) from Guenther et al. [1995]
!     * separate accounting of sunlit and shaded leaves using equation (12)
!       of Guenther et al. [1995]
!     * partitioning of the radiation at the top of the canopy into direct
!       and diffuse components using a parameterization of results from an
!       atmospheric radiative transfer model [Wang, 1996].
!
! ARGUMENTS
!   cloud_frac1 : fractional cloud cover
!   suncos1     : cosine of the solar zenith angle
!   xlai1       : leaf area index of land type for current month
!   coeff       : factor that integrates the light dependence over the canopy
!                 depth; "sp" even though "ri" is input per unit area of leaf
!
!-----------------------------------------------------------------------------

      function Biofit  &
     &  (cloud_frac1, suncos1, xlai1, coeff)

      implicit none

#     include "gmi_emiss_constants.h"

!     ----------------------
!     Argument declarations.
!     ----------------------

      real*8  :: cloud_frac1
      real*8  :: suncos1
      real*8  :: xlai1
      real*8  :: coeff(NPOLY)


!     ----------------------
!     Function declarations.
!     ----------------------

      real*8  :: Biofit


!     -----------------------
!     Parameter declarations.
!     -----------------------

      integer, parameter ::  &
     &  BIOFIT_DIM = 4   ! Biofit dimension


!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: ii, nn
      integer :: k0, k1, k2, k3
      integer :: spdim

      real*8  :: bigterm(NPOLY)

      real*8  :: spterm (BIOFIT_DIM-1)

      real*8  :: term   (BIOFIT_DIM)


!     ----------------
!     Begin execution.
!     ----------------

      term(1) = 1.0d0
      term(2) = xlai1
      term(3) = suncos1
      term(4) = cloud_frac1


      spdim = BIOFIT_DIM - 1


      do ii = 1, spdim
        spterm(ii) = term(ii+1)
      end do


!     ==============
      call Sun_Param  &
!     ==============
     &  (spdim, spterm)


      do ii = 1, spdim
        term(ii+1) = spterm(ii)
      end do


      k0 = 0

      do k3 = 1, BIOFIT_DIM
        do k2 = k3, BIOFIT_DIM
          do k1 = k2, BIOFIT_DIM

            k0 = k0 + 1

            bigterm(k0) = term(k1) * term(k2) * term(k3)

          end do
        end do
      end do


      Biofit = 0.0d0

      do nn = 1, NPOLY

        Biofit = Biofit + (coeff(nn) * bigterm(nn))

      end do

      Biofit = Max (Biofit, 0.1d0)


      return

      end


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Sun_Param
!
! DESCRIPTION
!   This routine scales and constrains xlai, suncos, and cloud_frac; called
!   by Biofit.
!
! ARGUMENTS
!   spdim  : spterm dimension
!   spterm : array of terms containing xlai, suncos, cloud_frac
!
!-----------------------------------------------------------------------------

      subroutine Sun_Param  &
     &  (spdim, spterm)

      use GmiPrintError_mod, only : GmiPrintError

      implicit none

!     ----------------------
!     Argument declarations.
!     ----------------------

      integer :: spdim
      real*8  :: spterm(spdim)


!     -----------------------
!     Parameter declarations.
!     -----------------------

      integer, parameter ::  &
     &  SPDIM_PARAM = 3

      integer, parameter ::  &
     &  ND (SPDIM_PARAM) =       & ! scaling factor for each variable
     &    (/ 55, 20, 11 /)

      real*8,  parameter ::  &
     &  XHI(SPDIM_PARAM) =       & ! maximum for each variable
     &    (/ 11.0d0, 1.0d0, 1.0d0 /)


!     ----------------------
!     Variable declarations.
!     ----------------------

      character (len=75) :: err_msg

      integer :: ii

      real*8  :: rnd
      real*8  :: xlow          ! minimum for each variable


!     ----------------
!     Begin execution.
!     ----------------

      if (spdim /= SPDIM_PARAM) then
        err_msg = 'spdim/SPDIM_PARAM problem in Sun_Param.'
        call GmiPrintError  &
     &    (err_msg, .true., 2, spdim, SPDIM_PARAM, 0, 0.0d0, 0.0d0)
      end if


      do ii = 1, spdim

        spterm(ii) = Min (spterm(ii), XHI(ii))

        if (ii /= spdim) then

          rnd  = ND (ii)

          xlow = XHI(ii) / rnd

        else

          xlow = 0.0d0

        end if

        spterm(ii) = Max (spterm(ii), xlow)

        spterm(ii) = spterm(ii) / XHI(ii)

      end do


      return

      end


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Diffg
!
! DESCRIPTION
!   This routine calculates the molecular diffusivity of a gas in air
!   (m^2/s).
!
!   The molecular radius of air is given in a table on p. 479 of Levine
!   [1988]; the table also gives radii for some other molecules.  Rather
!   than using a specific molecular radius a generic value is used for all
!   molecules, which is good enough in terms of calculating the diffusivity
!   as long as the molecule is not too big.
!
! ARGUMENTS
!   tk    : temperature (degK)
!   press : pressure    (Pa)
!   xm    : molecular weight of gas (kg)
!
!-----------------------------------------------------------------------------

      function Diffg  &
     &  (tk, press, xm)

      implicit none

#     include "gmi_phys_constants.h"

!     ----------------------
!     Argument declarations.
!     ----------------------

      real*8  :: tk
      real*8  :: press
      real*8  :: xm


!     ----------------------
!     Function declarations.
!     ----------------------

      real*8  :: Diffg


!     -----------------------
!     Parameter declarations.
!     -----------------------

      real*8, parameter ::  &
     &  RADAIR = 1.2d-10,    & ! hard-sphere molecular radii of air (m)
     &  RADX   = 1.5d-10   ! hard-sphere molecular radii of the diffusing
                           ! gas (m)


!     ----------------------
!     Variable declarations.
!     ----------------------

      real*8  :: airden
      real*8  :: diam
      real*8  :: frpath
      real*8  :: speed
      real*8  :: zz


!     ----------------
!     Begin execution.
!     ----------------

      airden = press * AVOGAD / (GAS_CONST_J * tk)


!     -------------------------------------------------------------------
!     Calculate the mean free path for gas X in air:  eq. 8.5 of Seinfeld
!     [1986]; diam is the collision diameter for gas X with air.
!     -------------------------------------------------------------------

      zz = xm / (MWTAIR * KGPG)

      diam = RADX + RADAIR

      frpath = 1.0d0 / (GMI_PI * Sqrt (1.0d0 + zz ) *  &
     &         airden * (diam * diam))


!     -------------------------------------------------------------
!     Calculate average speed of gas X; eq. 15.47 of Levine [1988].
!     -------------------------------------------------------------

      speed = Sqrt (8.0d0 * GAS_CONST_J * tk / (GMI_PI * xm))


!     --------------------------------------------------------------------
!     Calculate diffusion coefficient of gas X in air; eq. 8.9 of Seinfeld
!     [1986].
!     --------------------------------------------------------------------

      Diffg = (3.0d0 * GMI_PI / 32.0d0) *  &
     &        (1.0d0 + zz) * frpath * speed


      return

      end

