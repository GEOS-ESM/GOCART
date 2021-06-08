module solarZenithAngle_mod

      IMPLICIT NONE
      
      PRIVATE
      PUBLIC  :: computeSolarZenithAngle

CONTAINS

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: computeSolarZenithAngle
!
! !INTERFACE:
!
      function computeSolarZenithAngle (jday, time_sec, &
                   offset_sec, latDeg, lonDeg, i1, i2, j1, j2) &
                   result(this_)
!
! !INPUT PARAMETERS:
      integer :: i1, i2, j1, j2
      integer :: jday             ! day of year (1-365)
      real*8  :: time_sec         ! current time in seconds
      real*8  :: offset_sec       ! offset from tau at which to do photolysis (s)
      real*8  :: latDeg(j1:j2)    ! latitude (deg)
      real*8  :: lonDeg(i1:i2)    ! longitude (deg)
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
    REAL*8 :: sinlat(j1:j2), sollat(j1:j2), coslat(j1:j2)
    REAL*8 :: cosz(i1:i2, j1:j2)
    real*8 :: tau, timej, loct
    integer        :: ii
!EOP
!------------------------------------------------------------------------------
!BOC
      tau   = time_sec         / 3600.0d0
      timej = offset_sec / 3600.0d0
  
      sindec=0.3978d0*sin(0.9863d0*(dble(jday)-80.d0)*PI180)
      soldek=asin(sindec)
      cosdec=cos(soldek)
      sinlat(j1:j2)=sin(latDeg(j1:j2)*PI180)
      sollat(j1:j2)=asin(sinlat(j1:j2))
      coslat(j1:j2)=cos(sollat(j1:j2))
!
      do ii = i1, i2
         loct            = (((tau+timej)*15.d0)-180.d0)*PI180 + (lonDeg(ii)*PI180)
         cosz(ii,j1:j2)  = cosdec*coslat(j1:j2)*cos(loct) + sindec*sinlat(j1:j2)
         this_(ii,j1:j2) = acos(cosz(ii,j1:j2))/PI180
      enddo

      end function computeSolarZenithAngle
!EOC
!------------------------------------------------------------------------------
end module solarZenithAngle_mod
