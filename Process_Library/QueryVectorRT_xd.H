      integer                              :: ch
      integer, allocatable, dimension(:)   :: irh, irhp1, isnap
      real, allocatable, dimension(:)      :: rh_, arh
      real, allocatable, dimension(:,:,:)  :: pmomIn
      integer :: i, j
      real, allocatable :: momIn1(:,:,:), momIn2(:,:,:)

      if ( present(rc) ) rc = 0

      ch = this%getChannel(wavelength)

!     Now linearly interpolate the input table for the requested aerosol and
!     channel; rh is the relative humidity.

      call this%Query( wavelength, bin, q_mass, rh, &
                 tau,  ssa, gasym, bext,bsca, bbck, &
                 reff, p11, p22,   gf,  rhop, rhod, &
                 vol,  area, refr, refi, rc )

!     Now map the input RH to the high resolution hash table for RH

      rh_ = reshape(max(min(rh,0.99),0.), [product(shape(rh))])
      isnap = int((rh_+0.001)*1000.)
      where(isnap .lt. 1) isnap = 1
      arh   = this%rha( isnap )
      irh   = this%rhi( isnap )
      irhp1 = irh+1
      where(irhp1 .gt. this%nrh) irhp1 = this%nrh

!     Now linearly interpolate the input table for the requested aerosol and
!     channel; rh is the relative humidity.

      allocate(momIn1(size(irh),size(this%pmom,4),size(this%pmom,5)))
      allocate(momIn2(size(irh),size(this%pmom,4),size(this%pmom,5)))
      allocate(pmomIn(size(irh),size(this%pmom,4),size(this%pmom,5)))
      momIn1 = this%pmom(irh   ,ch, bin,:,:)
      momIn2 = this%pmom(irhp1 ,ch, bin,:,:)
      do j = 1, size(momIn1,3)
         do i = 1, size(momIn1,2)
            pmomIn(:,i,j)  =  momIn1(:,i,j) * (1.-arh) &
                            + momIn2(:,i,j) * arh
         enddo
      enddo

      pmom  = reshape(pmomIn, [shape(q_mass),size(pmomIn,2),size(pmomIn,3)])
