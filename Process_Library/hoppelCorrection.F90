   subroutine hoppelCorrection (radius, rhop, rh, dz, ustar, rhFlag, &
                                airdens, t, grav, karman, fhoppel, rc)

! !USES:
   implicit NONE

! !INPUT PARAMETERS:
   real, intent(in)     :: radius    ! dry radius [m]
   real, intent(in)     :: rhop      ! dry density [kg m-3]
   integer, intent(in)  :: rhFlag    ! 1 (Fitzgerald, 1975)
                                     ! 2 (Gerber, 1985)
   real, dimension(:,:), intent(in)  :: rh    ! relative humidity [0-1]
   real, dimension(:,:), intent(in)  :: dz    ! surface layer height [m]
   real, dimension(:,:), intent(in)  :: ustar ! surface velocity scale [m s-1]
   real, dimension(:,:), intent(in)  :: airdens ! air density [kg/m^3]s
   real, dimension(:,:), intent(in)  :: t  ! temperature [k]
   real, intent(in)  :: grav    ! gravity [m/sec^2]
   real, intent(in)  :: karman  ! Von Karman constant [unitless]


! !INOUTPUT PARAMETERS:
   real, dimension(:,:), intent(inout) :: fhoppel

! !OUTPUT PARAMETERS:
   integer, intent(out) :: rc

! !Local Variables
   real    :: radius_wet ! humidified radius [m]
   real    :: rhop_wet   ! wet density [kg m-3]
   real    :: diff_coef
   real, allocatable, dimension(:,:) ::  vsettle
   integer :: i, j

   integer :: status


!EOP
!------------------------------------------------------------------------------------
!  Begin..

   rc = __SUCCESS__
   fhoppel = 1.0
   allocate(vsettle, mold=rh)

   do j = 1, ubound(rh,2)
      do i = 1, ubound(rh,1)
         call wetRadius (radius, rhop, rh(i,j), rhFlag, &
                         radius_wet, rhop_wet, __RC_NO_OPT__)
         call Chem_CalcVsettle2Gorig (radius_wet, rhop_wet, airdens(i,j), t(i,j), &
                                      GRAV, diff_coef, vsettle(i,j))
         fhoppel(i,j) = (10./dz(i,j)) ** (vsettle(i,j)/KARMAN/ustar(i,j))
      end do
   end do


   deallocate(vsettle)

   end subroutine hoppelCorrection
