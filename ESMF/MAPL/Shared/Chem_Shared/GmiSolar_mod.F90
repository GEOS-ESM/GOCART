module GmiSolar_mod

      implicit none

      private
      public  :: CalcSolarDeclination
      public  :: CalcSolarZenithAngle
      public  :: CalcCosSolarZenithAngle
      public  :: computeSolarZenithAngle_Photolysis

      contains

!=============================================================================
!
! CODE DEVELOPER
!   John Tannahill, LLNL (Original code from Keith Grant, LLNL)
!   jrt@llnl.gov
!
!=============================================================================
!
! Mar 30, 2017: Moved this file from GmiSupportingModules/ to Chem_Shared/ for TR
!
!-----------------------------------------------------------------------------
!
! ROUTINE
!   CalcSolarDeclination
!
! DESCRIPTION
!   Given the Julian day, this routine calculates the solar declination and
!   the square of the ratio of the mean Earth-Sun distance to the current
!   Earth-Sun distance.   Refer to:  Paltridge, G.W., and C.M.R. Platt, 1976:
!   "Radiative Processes in Meteorology and Climatology", Elsevier, pp. 57-63.
!
! ARGUMENTS
!   julday  : Julian day counting ranging from 0 on Jan. 1st to 364 on
!             Dec. 31st
!   decl    : solar declination (deg)
!   rdistsq : the square of the ratio of the mean Earth-Sun distance
!             to the current Earth-Sun distance
!
!-----------------------------------------------------------------------------

      subroutine CalcSolarDeclination  &
     &  (julday, decl, rdistsq)

      implicit none

#     include "gmi_phys_constants.h"


!     ----------------------
!     Argument declarations.
!     ----------------------

      real*8, intent(in )  :: julday
      real*8, intent(out)  :: decl
      real*8, intent(out)  :: rdistsq

!     -----------------------
!     Parameter declarations.
!     -----------------------

      real*8, parameter :: TWO_PI = 2.0d0 * GMI_PI

!     ----------------------
!     Variable declarations.
!     ----------------------

      real*8  :: theta           ! year angle on specified day
      real*8  :: theta2, theta3

!     ----------------
!     Begin execution.
!     ----------------

!c?   Leap year?
      theta  = (TWO_PI * julday) / 365.0d0

      theta2 = 2.0d0 * theta
      theta3 = 3.0d0 * theta

      decl =  &
     &  (360.0d0 / TWO_PI) *  &
     &  (0.006918d0 -  &
     &   0.399912d0 * Cos (theta)  + 0.070257d0 * Sin (theta)  -  &
     &   0.006758d0 * Cos (theta2) + 0.000907d0 * Sin (theta2) -  &
     &   0.002697d0 * Cos (theta3) + 0.001480d0 * Sin (theta3))

      rdistsq =  &
     &  1.000110d0  +  &
     &  0.034221d0  * Cos (theta)  + 0.001280d0 * Sin (theta)  +  &
     &  0.000719d0  * Cos (theta2) + 0.000077d0 * Sin (theta2)

      return

      end subroutine CalcSolarDeclination

!-----------------------------------------------------------------------------
!
! ROUTINE
!   Solrza
!
! DESCRIPTION
!   Given a Greenwich time, a solar declination and reference latitude and
!   longitudes,  this routine returns a corresponding list of cosines of solar
!   zenith angles.  Refer to:  Paltridge, G.W., and C.M.R. Platt, 1976:
!   "Radiative Processes in Meteorology and Climatology", Elsevier, p. 62.
!
! ARGUMENTS
!   time   : Greenwich time since Jan 1, counting zero from midnight (days)
!   decl   : solar declination (deg)
!   lat    : latitude   (deg)
!   lon    : longitudes (deg)
!   nn     : number of longitudes for which to calculate cosines of the
!            solar zenith angle
!   cossza : cosines of the solar zenith angle (output)
!
!-----------------------------------------------------------------------------

      subroutine CalcSolarZenithAngle  &
     &  (time, decl, lat, lon, nn, cossza)

      implicit none

#     include "gmi_phys_constants.h"

!     ----------------------
!     Argument declarations.
!     ----------------------

      integer, intent(in ) :: nn
      real*8 , intent(in ) :: time
      real*8 , intent(in ) :: decl
      real*8 , intent(in ) :: lat
      real*8 , intent(in ) :: lon   (nn)
      real*8 , intent(out) :: cossza(nn)

!     -----------------------
!     Parameter declarations.
!     -----------------------

      real*8,  parameter :: TWO_PI = 2.0d0 * GMI_PI
      real*8,  parameter :: PI_180 = TWO_PI / 360.0d0

!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: ii

      real*8  :: cosha

!     ----------------
!     Begin execution.
!     ----------------

      do ii = 1, nn

!       -------------------------------------------------------------
!       Calculate the cosine of the hour angle which is referenced to
!       noon.
!       -------------------------------------------------------------

        cosha =  &
     &    Cos (TWO_PI *  &
     &         (Mod (time + (lon(ii) / 360.d0), 1.0d0) - 0.5d0))


!       -----------------------------------------------
!       Calculate the cosine of the solar zenith angle.
!       -----------------------------------------------

        cossza(ii) =  &
     &    Sin (PI_180 * lat) * Sin (PI_180 * decl) +  &
     &    Cos (PI_180 * lat) * Cos (PI_180 * decl) * cosha

      end do


      return

      end subroutine CalcSolarZenithAngle

!-----------------------------------------------------------------------------
!
! ROUTINE
!  CalcCosSolarZenithAngle 
!
! DESCRIPTION
!   Given a Greenwich time, a solar declination and reference latitude and
!   longitudes,  this routine returns a corresponding list of cosines of solar
!   zenith angles.  Refer to:  Paltridge, G.W., and C.M.R. Platt, 1976:
!   "Radiative Processes in Meteorology and Climatology", Elsevier, p. 62.
!
! ARGUMENTS
!   time   : Greenwich time since Jan 1, counting zero from midnight (days)
!   latDeg : latitude   (deg)
!   lonDeg : longitudes (deg)
!   i1,i2  : Dimensions of latDeg, lonDeg
!    j1,j2
!   cossza : cosines of the solar zenith angle (output)
!
!-----------------------------------------------------------------------------

      subroutine CalcCosSolarZenithAngle  &
     &  (time, latDeg, lonDeg, cosSolarZenithAngle, i1, i2, j1, j2)

      implicit none

#     include "gmi_phys_constants.h"

!     ----------------------
!     Argument declarations.
!     ----------------------

      integer, intent(in ) :: i1, i2, j1, j2
      real*8 , intent(in ) :: time
      real*8 , intent(in ) :: latDeg(i1:i2,j1:j2)
      real*8 , intent(in ) :: lonDeg(i1:i2,j1:j2) 
      real*8 , intent(out) :: cosSolarZenithAngle(i1:i2,j1:j2)

!     -----------------------
!     Parameter declarations.
!     -----------------------

      real*8,  parameter :: TWO_PI = 2.0d0 * GMI_PI
      real*8,  parameter :: PI_180 = TWO_PI / 360.0d0

!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: ii, ij
      real*8  :: decl, theta, theta2, theta3
      real*8  :: cosha(i1:i2,j1:j2), coslat, sinlat

!     ----------------
!     Begin execution.
!     ----------------

!c?   Leap year?
      theta  = (TWO_PI * time) / 365.0d0

      theta2 = 2.0d0 * theta
      theta3 = 3.0d0 * theta

!     Calculate the solar declination (deg)

      decl =  &
     &  (360.0d0 / TWO_PI) *  &
     &  (0.006918d0 -  &
     &   0.399912d0 * Cos (theta)  + 0.070257d0 * Sin (theta)  -  &
     &   0.006758d0 * Cos (theta2) + 0.000907d0 * Sin (theta2) -  &
     &   0.002697d0 * Cos (theta3) + 0.001480d0 * Sin (theta3))

!     -------------------------------------------------------------------
!     Calculate the cosine of the hour angle which is referenced to noon.
!     -------------------------------------------------------------------

      cosha (i1:i2,j1:j2) =  Cos(TWO_PI * (Mod (time + (lonDeg(i1:i2,j1:j2) / 360.d0), 1.0d0) - 0.5d0))

!     -----------------------------------------------
!     Calculate the cosine of the solar zenith angle.
!     -----------------------------------------------

      do ij = j1, j2
       do ii = i1, i2
        sinlat =  Sin (PI_180 * latDeg(ii,ij)) * Sin (PI_180 * decl)
        coslat =  Cos (PI_180 * latDeg(ii,ij)) * Cos (PI_180 * decl)
        cosSolarZenithAngle(ii,ij) =  sinlat + coslat*cosha(ii,ij)
       end do
      end do

      return

      end subroutine CalcCosSolarZenithAngle
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: computeSolarZenithAngle_Photolysis
!
! !INTERFACE:
!
      function computeSolarZenithAngle_Photolysis (jday, time_sec, &
                   fastj_offset_sec, latDeg, lonDeg, i1, i2, j1, j2) &
                   result(this_)
!
! !INPUT PARAMETERS:
      integer :: i1, i2, j1, j2
      integer :: jday             ! day of year (1-365)
      real*8  :: time_sec         ! current time in seconds
      real*8  :: fastj_offset_sec ! offset from tau at which to do photolysis (s)
      real*8  :: latDeg(i1:i2,j1:j2)    ! latitude (deg)
      real*8  :: lonDeg(i1:i2,j1:j2)    ! longitude (deg)
!
! !RETURNED VALUE
      real*8  :: this_(i1:i2,j1:j2)
!
! !DESCRIPTION:
!  Computes the solar zenith angle used in the photolysis package.
!
! !DEFINED PARAMETERS:
      real*8, parameter :: locPI = 3.141592653589793D0
      real*8, parameter :: PI180 = locPI / 180.0d0
!
! !LOCAL VARIABLES:
    REAL*8 :: sindec, soldek, cosdec
    REAL*8 :: sinlat(i1:i2,j1:j2), sollat(i1:i2,j1:j2), coslat(i1:i2,j1:j2)
    REAL*8 :: cosz(i1:i2,j1:j2)
    real*8 :: tau, timej, loct
    integer        :: ii, ij
!EOP
!------------------------------------------------------------------------------
!BOC
      tau   = time_sec         / 3600.0d0
      timej = fastj_offset_sec / 3600.0d0
  
      sindec=0.3978d0*sin(0.9863d0*(dble(jday)-80.d0)*PI180)
      soldek=asin(sindec)
      cosdec=cos(soldek)
      sinlat(i1:i2,j1:j2)=sin(latDeg(i1:i2,j1:j2)*PI180)
      sollat(i1:i2,j1:j2)=asin(sinlat(i1:i2,j1:j2))
      coslat(i1:i2,j1:j2)=cos(sollat(i1:i2,j1:j2))

      do ij = j1, j2
       do ii = i1, i2
         loct            = (((tau+timej)*15.d0)-180.d0)*PI180 + (lonDeg(ii,ij)*PI180)
         cosz(ii,ij)  = cosdec*coslat(ii,ij)*cos(loct) + sindec*sinlat(ii,ij)
         this_(ii,ij) = acos(cosz(ii,ij))/PI180
       enddo
      enddo

      end function computeSolarZenithAngle_Photolysis
!EOC
!------------------------------------------------------------------------------
end module GmiSolar_mod
