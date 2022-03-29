   subroutine wetRadius (radius, rhop, rh, flag, radius_wet, rhop_wet, rc)

! !USES:
   implicit NONE

! !INPUT PARAMETERS:
   real, intent(in)  :: radius    ! dry radius [m]
   real, intent(in)  :: rhop      ! dry density [kg m-3]
   real, intent(in)  :: rh        ! relative humidity [0-1]
   integer           :: flag      ! 1 (Fitzgerald, 1975)
                                  ! 2 (Gerber, 1985)

! !OUTPUT PARAMETERS:
   real, intent(out) :: radius_wet ! humidified radius [m]
   real, intent(out) :: rhop_wet   ! wet density [kg m-3]
   integer, intent(out) :: rc

! !Local Variables
   real :: sat, rcm, rrat
   real, parameter ::  rhow = 1000.  ! Density of water [kg m-3]

!  The following parameters relate to the swelling of seasalt like particles
!  following Fitzgerald, Journal of Applied Meteorology, 1975.
   real, parameter :: epsilon = 1.   ! soluble fraction of deliqeuscing particle
   real, parameter :: alphaNaCl = 1.35
   real :: alpha, alpha1, alpharat, beta, theta, f1, f2

!  parameter from Gerber 1985 (units require radius in cm, see rcm)
   real, parameter :: c1=0.7674, c2=3.079, c3=2.573e-11, c4=-1.424

!EOP
!------------------------------------------------------------------------------------
!  Begin...

   rc = __SUCCESS__

!  Default is to return radius as radius_wet, rhop as rhop_wet
   radius_wet = radius
   rhop_wet   = rhop

!  Make sure saturation ratio (RH) is sensible
   sat = max(rh,tiny(1.0)) ! to avoid zero FPE

!  Fitzgerald Scheme
   if(flag .eq. 1 .and. sat .ge. 0.80) then
!     parameterization blows up for RH > 0.995, so set that as max
!     rh needs to be scaled 0 - 1
      sat = min(0.995,sat)
!     Calculate the alpha and beta parameters for the wet particle
!     relative to amonium sulfate
      beta = exp( (0.00077*sat) / (1.009-sat) )
      if(sat .le. 0.97) then
         theta = 1.058
      else
         theta = 1.058 - (0.0155*(sat-0.97)) /(1.02-sat**1.4)
      endif
      alpha1 = 1.2*exp( (0.066*sat) / (theta-sat) )
      f1 = 10.2 - 23.7*sat + 14.5*sat**2.
      f2 = -6.7 + 15.5*sat - 9.2*sat**2.
      alpharat = 1. - f1*(1.-epsilon) - f2*(1.-epsilon**2.)
      alpha = alphaNaCl * (alpha1*alpharat)
!     radius_wet is the radius of the wet particle
      radius_wet = alpha * radius**beta
      rrat       = (radius/radius_wet)**3.
      rhop_wet   = rrat*rhop + (1.-rrat)*rhow
   elseif(flag .eq. 2) then   ! Gerber
      sat = min(0.995,sat)
      rcm = radius*100.
      radius_wet = 0.01 * (c1*rcm**c2 / (c3*rcm**c4-alog10(sat)) &
                           + rcm**3.)**(1./3.)
      rrat       = (radius/radius_wet)**3.
      rhop_wet   = rrat*rhop + (1.-rrat)*rhow
   endif

 end subroutine wetRadius
