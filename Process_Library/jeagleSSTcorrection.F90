   subroutine jeagleSSTcorrection(sstEmisFlag, fsstemis, ts, rc)

! !USES:
  implicit NONE

! !INPUT/OUTPUT PARAMETERS:
  real, dimension(:,:), intent(inout) :: fsstemis     !

! !INPUT PARAMETERS:
   integer, intent(in)                       :: sstEmisFlag  ! 1 or 2
   real, dimension(:,:), intent(in)          :: ts  ! surface temperature (K)

! !OUTPUT PARAMETERS:
   integer, optional, intent(out) :: rc
!EOP

! !Local Variables
   real, allocatable, dimension(:,:) :: tskin_c
!EOP
!-------------------------------------------------------------------------
!  Begin...

   fsstemis = 1.0

   if (sstemisFlag == 1) then          ! SST correction folowing Jaegle et al. 2011
      fsstemis = 0.0

      allocate( tskin_c, mold=fsstemis )
      tskin_c  = ts - 273.15
      fsstemis = (0.3 + 0.1*tskin_c - 0.0076*tskin_c**2 + 0.00021*tskin_c**3)

      where(fsstemis < 0.0) fsstemis = 0.0

      deallocate( tskin_c )
   else if (sstemisFlag == 2) then     ! GEOS5 SST correction
      fsstemis = 0.0

      allocate( tskin_c, mold=fsstemis )
      tskin_c  = ts - 273.15

      where(tskin_c < -0.1) tskin_c = -0.1    ! temperature range (0, 36) C
      where(tskin_c > 36.0) tskin_c = 36.0    !

      fsstemis = (-1.107211 -0.010681*tskin_c -0.002276*tskin_c**2 + 60.288927*1.0/(40.0 - tskin_c))
      where(fsstemis < 0.0) fsstemis = 0.0
      where(fsstemis > 7.0) fsstemis = 7.0

      deallocate( tskin_c )
   end if

   __RETURN__(__SUCCESS__)
   end subroutine jeagleSSTcorrection
