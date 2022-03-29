   subroutine SulfateUpdateOxidants (nymd_current, nhms_current, lonRad, latRad, &
                                     rhoa, km, cdt, nymd_last, &
                                     undefval, radToDeg, nAvogadro, pi, airMolWght, &
                                     oh_clim, no3_clim, h2o2_clim, &
                                     xoh, xno3, xh2o2, recycle_h2o2, rc)
! !USES:
   implicit NONE

! !INPUT PARAMETERS:
   integer, intent(in)    :: nymd_current, &   ! current model NYMD
                             nhms_current      ! current model NHMS
   real, dimension(:,:), intent(in)   :: lonRad, latRad ! model grid lon and lat
   real, dimension(:,:,:), intent(in) :: rhoa           ! layer air density [kg/m^3]
   integer, intent(in)    :: km         ! number of model levels
   real, intent(in)       :: cdt        ! chemistry model time-step
   integer, intent(inout) :: nymd_last  ! NYMD of last emission update
   real, intent(in)       :: undefval   ! value for undefined values
   real, intent(in)       :: radToDeg   ! radian to degrees conversion
   real, intent(in)       :: nAvogadro  ! Avogadro's number [molecules per mole of air]
   real, intent(in)       :: pi         ! pi constant
   real, intent(in)       :: airMolWght ! air molecular weight [kg/Kmole]
   real, pointer, dimension(:,:,:) :: oh_clim, &   ! climatological OH
                                      no3_clim, &  ! climatological NO3
                                      h2o2_clim  ! climatological H2O2
   real, dimension(:,:,:), intent(inout) :: xoh, xno3, xh2o2 ! returned oxidant values
   logical, intent(inout) :: recycle_h2o2

! !OUTPUT PARAMETERS:
  integer, optional, intent(out)   :: rc    ! Error return code:
                                            !  0 - all is well

! !DESCRIPTION: Update Oxidant Fields for Sulfate
!               We have 3 oxidant fields (OH, NO3, H2O2) which may come
!               from either a climatological file or from interactive GMI.
!               IF from climatology, update (reset) values from climatology
!               if necessary (e.g., for a new day) and set to current values
!               needed by chemistry.
!               IF from GMI read as is

!
! !REVISION HISTORY:
! ???        ???       - Legacy code
! 23July2020 E.Sherman - ported/refactored for use in process library.
!

! !Local Variables
   integer :: i, j, k, jday
   real    :: qmax, xhour, xhouruse
   real, dimension(:,:), allocatable  :: cossza, sza
   real, dimension(:,:), allocatable  :: tcosz, tday, tnight
   integer :: n, ndystep
   integer :: i1=1, j1=1, i2, j2

!EOP
!-------------------------------------------------------------------------
!  Begin...

    i2 = size(rhoa,1)
    j2 = size(rhoa,2)

    allocate(cossza(i1:i2,j1:j2), sza(i1:i2,j1:j2), tcosz(i1:i2,j1:j2), &
             tday(i1:i2,j1:j2), tnight(i1:i2,j1:j2))

! Update emissions/production if necessary (daily)
!  -----------------------------------------------
!   Oxidant fields
!   The expectation here is that OH is being read in the form
!   volume mixing ratio from a file (so, like GMI would provide).
!   Below, in the scaling by solar zenith angle, we convert from
!   VMR to # cm-3 expected by the chemistry.
    where(1.01*oh_clim(i1:i2,j1:j2,1:km) > undefval) oh_clim(i1:i2,j1:j2,1:km) = 0.
    where(     oh_clim(i1:i2,j1:j2,1:km) < 0       ) oh_clim(i1:i2,j1:j2,1:km) = 0.

    where(1.01*no3_clim(i1:i2,j1:j2,1:km) > undefval) no3_clim(i1:i2,j1:j2,1:km) = 0.
    where(     no3_clim(i1:i2,j1:j2,1:km) < 0       ) no3_clim(i1:i2,j1:j2,1:km) = 0.

    where(1.01*h2o2_clim(i1:i2,j1:j2,1:km) > undefval) h2o2_clim(i1:i2,j1:j2,1:km) = 0.
    where(     h2o2_clim(i1:i2,j1:j2,1:km) < 0       ) h2o2_clim(i1:i2,j1:j2,1:km) = 0.

!   The first time through the reads we will save the h2o2 monthly
!   average in the instantaneous field
!   ---------------------------------
    if (nymd_last == nymd_current) then
       xh2o2 = h2o2_clim
       nymd_last = nymd_current
    end if

!   Find the day number of the year and hour (needed for later doing sza)
!   ----------------------------------
    jday = idaynum(nymd_current)
    xhour = (  real(nhms_current/10000)*3600. &
             + real(mod(nhms_current,10000)/100)*60. &
             + real(mod(nhms_current,100)) &
             ) / 3600.

!   Recycle H2O2 to input on 3 hour boundaries if not coupled to GMI
!   ----------------------------------
    if (recycle_h2o2) then
       xh2o2 = h2o2_clim
       recycle_h2o2 = .false.
    end if

!   If not getting instantaneous values from GMI, update for time of day.
!   ---------------------------------------------------------------------
!   OH
    xoh = oh_clim
    cossza(:,:) = 0.

!   Want to find the sum of the cos(sza) for use in scaling OH diurnal variation
!   tcosz is the sum of cossza over the whole day
!   tday is the time of day spent in light
!   Requires integrating over future times, so cannot use w_c%cosz
    xHourUse = xHour
    ndystep = 86400. / cdt
    tcosz(:,:) = 0.
    tday(:,:) = 0.
    do n = 1, ndystep
       call szangle(jday,xHourUse,lonRad,latRad,PI,radToDeg,sza,cossza, i2, j2)
       tcosz = tcosz + cossza
       xHourUse = xHourUse + cdt/3600.
       if(xHourUse .gt. 24.) xHourUse = xHourUse - 24.
!      Find the daylight portion of the day
       do j = j1, j2
          do i = i1, i2
             if(cossza(i,j) .gt. 0.) tday(i,j) = tday(i,j) + cdt
          end do
       end do
    end do

!   Find the cos(sza) now for use in scaling OH and NO3
    call szangle(jday,xHour,lonRad,latRad,PI,radToDeg,sza,cossza, i2, j2)

    tnight(i1:i2,j1:j2) = (86400.-tday(i1:i2,j1:j2))

    do k = 1, km
       where (tcosz(i1:i2,j1:j2) > 0)
          xoh(i1:i2,j1:j2,k) = oh_clim(i1:i2,j1:j2,k)*(86400./cdt)*cossza(i1:i2,j1:j2) / tcosz(i1:i2,j1:j2)
       elsewhere
          xoh(i1:i2,j1:j2,k) = 0.00
       end where
    end do
    where(xoh(i1:i2,j1:j2,1:km) < 0.00) xoh(i1:i2,j1:j2,1:km) = 0.00

!   To go from volume mixing ratio to # cm-3 (expected in chemistry)
!   include the following line
    xoh = xoh * 1000.*rhoa / airMolWght * nAvogadro * 1.e-6

!   NO3
    xno3 = no3_clim
    cossza(:,:) = 0.
    call szangle(jday,xHour,lonRad,latRad,PI,radToDeg,sza,cossza, i2, j2)

!   If there is daylight then no3 is small (assume zero) and the
!   average is distributed only over the night time portion

    do k=1,km
       where(cossza(i1:i2,j1:j2) > 0 .OR. tnight(i1:i2,j1:j2) < tiny(1.0))
          xno3(i1:i2,j1:j2,k) = 0.00
       elsewhere
          xno3(i1:i2,j1:j2,k) = no3_clim(i1:i2,j1:j2,k) * 86400./ tnight(i1:i2,j1:j2)
       end where
    end do

    __RETURN__(__SUCCESS__)
   end subroutine SulfateUpdateOxidants
