   subroutine SeasaltEmission ( rLow, rUp, method, u10m, v10m, ustar, pi, &
                                memissions, nemissions, rc )

! !DESCRIPTION: Calculates the seasalt mass emission flux every timestep.
!  The particular method (algorithm) used for the calculation is based
!  on the value of "method" passed on input.  Mostly these algorithms are
!  a function of wind speed and particle size (nominally at 80% RH).
!  Routine is called once for each size bin, passing in the edge radii
!  "rLow" and "rUp" (in dry radius, units of um).  Returned in the emission
!  mass flux [kg m-2 s-1].  A sub-bin assumption is made to break (possibly)
!  large size bins into a smaller space.
!
! !USES:

  implicit NONE

! !INPUT PARAMETERS:

   real, intent(in)             :: rLow, rUp   ! Dry particle bin edge radii [um]
   real, intent(in)             :: u10m(:,:)   ! 10-meter eastward wind [m s-1]
   real, intent(in)             :: v10m(:,:)   ! 10-m northward wind [m s-1]
   real, target, intent(in)     :: ustar(:,:)  ! friction velocity [m s-1]
   integer, intent(in)          :: method      ! Algorithm to use
   real, intent(in)             :: pi          ! pi constant

! !INOUTPUT PARAMETERS:
   real, dimension(:,:), intent(inout) :: memissions      ! Mass Emissions Flux [kg m-2 s-1]
   real, dimension(:,:), intent(inout) :: nemissions      ! Number Emissions Flux [# m-2 s-1]

! !OUTPUT PARAMETERS:
   integer, intent(out)          :: rc              ! Error return code:
                                                    !  0 - all is well
                                                    !  1 -
! !Local Variables
   integer       :: ir
   real, pointer :: w(:,:)                          ! Intermediary wind speed [m s-1]
   real          :: r, dr                           ! sub-bin radius spacing (dry, um)
   real          :: rwet, drwet                     ! sub-bin radius spacing (rh=80%, um)
   real          :: aFac, bFac, scalefac, rpow, exppow, wpow
   real, allocatable, dimension(:,:), target  :: w10m  ! 10-m wind speed [m s-1]

! !CONSTANTS
   real, parameter    :: r80fac = 1.65     ! ratio of radius(RH=0.8)/radius(RH=0.) [Gerber]
   real, parameter    :: rhop = 2200.      ! dry seasalt density [kg m-3]
!   real, parameter    :: pi = 3.1415       ! ratio of circumference to diameter of circle
   integer, parameter :: nr = 10                    ! Number of (linear) sub-size bins

   character(len=*), parameter :: myname = 'SeasaltEmission'

!EOP
!-------------------------------------------------------------------------
!  Begin...

   rc = __SUCCESS__

!  Define 10-m wind speed
   allocate(w10m, mold=u10m)
   w10m = sqrt(u10m*u10m + v10m*v10m)
!  Define the sub-bins (still in dry radius)
   dr = (rUp - rLow)/nr
   r  = rLow + 0.5*dr

!  Loop over size bins
   nemissions = 0.
   memissions = 0.

   do ir = 1, nr

    rwet  = r80fac * r
    drwet = r80fac * dr

    select case(method)

     case(1)  ! Gong 2003
      aFac     = 4.7*(1.+30.*rwet)**(-0.017*rwet**(-1.44))
      bFac     = (0.433-log10(rwet))/0.433
      scalefac = 1.
      rpow     = 3.45
      exppow   = 1.607
      wpow     = 3.41
      w        => w10m

     case(2)  ! Gong 1997
      aFac     = 3.
      bFac     = (0.380-log10(rwet))/0.650
      scalefac = 1.
      rpow     = 1.05
      exppow   = 1.19
      wpow     = 3.41
      w        => w10m

     case(3)  ! GEOS5 2012
      aFac     = 4.7*(1.+30.*rwet)**(-0.017*rwet**(-1.44))
      bFac     = (0.433-log10(rwet))/0.433
      scalefac = 33.0e3
      rpow     = 3.45
      exppow   = 1.607
      wpow     = 3.41 - 1.
      w        => ustar

     case default
      print *, 'GOCART2G_Process.F90 - SeasaltEmission - missing algorithm method'
      rc = __FAIL__
      return

    end select

!   Number emissions flux (# m-2 s-1)
    nemissions = nemissions + SeasaltEmissionGong( rwet, drwet, w, scalefac, aFac, bFac, rpow, exppow, wpow )

!   Mass emissions flux (kg m-2 s-1)
    scalefac = scalefac * 4./3.*pi*rhop*r**3.*1.e-18
    memissions = memissions + SeasaltEmissionGong( rwet, drwet, w, scalefac, aFac, bFac, rpow, exppow, wpow )

    r = r + dr

   end do

   deallocate(w10m)
  end subroutine SeasaltEmission 
