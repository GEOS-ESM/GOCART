   subroutine weibullDistribution (gweibull, weibullFlag, u10m, v10m, rc)

! !USES:
   implicit NONE

! !INPUT/OUTPUT PARAMETERS:
   real(kind=DP), dimension(:,:), intent(inout)    :: gweibull

! !INPUT PARAMETERS:
   logical, intent(in)                    :: weibullFlag
   real, dimension(:,:), intent(in)       :: u10m
   real, dimension(:,:), intent(in)       :: v10m


! !OUTPUT PARAMETERS:
   integer, optional, intent(out) :: rc

! !Descrption: The Weibull distribution correction ends up being a multiplicative constant
!  (g) times our present source function (see Eq. 12 in Fan & Toon, 2011 and notes for
!  (9/22/11). This constant is derived from the incomplete and complete forms of the gamma
!  function, hence the utilities pasted below.  The Weibull function and shape
!  parameters (k, c) assumed are from Justus 1978.

!EOP

! !Local Variables
   real(kind=DP)                 :: a, c, k, wt, x
   real(kind=DP), dimension(:,:), allocatable :: wm
   integer     :: i, j

   integer :: status

!EOP
!-------------------------------------------------------------------------
!  Begin...

   gweibull = 1.0

   allocate(wm(ubound(u10m, 1),ubound(u10m, 2)))
   wm = sqrt(u10m**2 + v10m**2)   ! mean wind speed
   wt = 4.d0                      ! a threshold (Fan & Toon, 2011)

   if (weibullFlag) then
       gweibull = 0.0

   do j = 1, ubound(u10m, 2)
      do i = 1, ubound(u10m, 1)
         if (wm(i,j) > 0.01) then
            k = 0.94d0 * sqrt(wm(i,j))         ! Weibull shape parameter
            c = wm(i,j) / gamma(1.d0 + 1.d0/k) ! Weibull shape parameter
            x = (wt / c) ** k
            a = 3.41d0 / k + 1.d0
            gweibull(i,j)  = (c / wm(i,j))**3.41d0 * igamma(a,x,__RC__)
         end if
      end do ! i
   end do ! j
   endif

   deallocate(wm)

   __RETURN__(__SUCCESS__)
   end subroutine weibullDistribution
