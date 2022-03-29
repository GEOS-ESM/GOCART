   subroutine phobicTophilic (aerosol_phobic, aerosol_philic, aerosol_toHydrophilic, &
                              km, cdt, grav, delp, rc)

! !USES:
   implicit NONE

! !INPUT PARAMETERS:
   integer, intent(in)   :: km   ! total model level
   real, intent(in)      :: cdt  ! chemistry model time-step [sec]
   real, intent(in)      :: grav ! [m/sec^2]
   real, dimension(:,:,:), intent(in)  :: delp  ! pressure thickness [Pa]

! !INOUTPUT PARAMETERS:
   real, dimension(:,:,:), intent(inout)  :: aerosol_phobic   ! OCphobic [kg kg-1]
   real, dimension(:,:,:), intent(inout)  :: aerosol_philic   ! OCphilic [kg kg-1]
   real, dimension(:,:), pointer   :: aerosol_toHydrophilic ! OCHYPHIL [kg m-2 s-1]
! !OUTPUT PARAMETERS:
   integer, optional, intent(out) :: rc

! !Local Variables
   integer :: i, j, k
   real :: qUpdate, delq

!EOP
!------------------------------------------------------------------------------------
!  Begin...

   if(associated(aerosol_toHydrophilic)) aerosol_toHydrophilic = 0.0

   do k = 1, km
    do j = 1, ubound(delp, 2)
     do i = 1, ubound(delp, 1)
      qUpdate = aerosol_phobic(i,j,k)*exp(-4.63e-6*cdt)
      qUpdate = max(qUpdate,1.e-32)
      delq = max(0.,aerosol_phobic(i,j,k)-qUpdate)
      aerosol_phobic(i,j,k) = qUpdate
      aerosol_philic(i,j,k) = aerosol_philic(i,j,k)+delq
      if(associated(aerosol_toHydrophilic)) &
       aerosol_toHydrophilic(i,j) = aerosol_toHydrophilic(i,j) &
        + delq*delp(i,j,k)/grav/cdt
     end do
    end do
   end do

   __RETURN__(__SUCCESS__)
  end subroutine phobicTophilic
