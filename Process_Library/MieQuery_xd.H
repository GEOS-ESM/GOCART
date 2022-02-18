      integer                      :: ch
      integer, allocatable, dimension(:)   :: irh, irhp1, isnap
      real, allocatable, dimension(:)      :: rh_, arh
      real, allocatable, dimension(:)      :: bextIn, bscaIn, bbckIn, gasymIn, p11In, p22In, &
                                              gfIn, rhopIn, rhodIn, volIn, areaIn, &
                                              refrIn, refiIn, rEffIn
      real, allocatable, dimension(:,:,:) :: pmomIn
      integer :: i, j

      if ( present(rc) ) rc = 0

      ch = this%getChannel(wavelength, rc=rc)

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

include "MieQuery.H"

      if (present(pmom)) then
         allocate(pmomIn(size(irh),size(this%pmom,4),size(this%pmom,5)))
         do j = 1, size(this%pmom,5)
            do i = 1, size(this%pmom,4)
               pmomIn(:,i,j)  =  this%pmom(irh   ,ch,bin,i,j) * (1.-arh) &
                               + this%pmom(irhp1 ,ch,bin,i,j) * arh
            enddo
         enddo
      endif

!     Fill the requested outputs
      if(present(tau  )) tau   = reshape(bextIn,shape(q_mass)) * q_mass
      if(present(ssa  )) ssa   = reshape(bscaIn/bextIn,shape(q_mass))
      if(present(bext )) bext  = reshape(bextIn,shape(q_mass))
      if(present(bsca )) bsca  = reshape(bscaIn,shape(q_mass))
      if(present(bbck )) bbck  = reshape(bbckIn,shape(q_mass))
      if(present(gasym)) gasym = reshape(gasymIn,shape(q_mass))
      if(present(rEff )) rEff  = reshape(1.E6*rEffIn,shape(q_mass)) ! convert to microns
      if(present(pmom )) pmom  = reshape(pmomIn, [shape(q_mass),size(pmomIn,2),size(pmomIn,3)])
      if(present(p11  )) p11   = reshape(p11In,shape(q_mass))
      if(present(p22  )) p22   = reshape(p22In, shape(q_mass))
      if(present(gf   )) gf    = reshape(gfIn, shape(q_mass))
      if(present(rhop )) rhop  = reshape(rhopIn, shape(q_mass))
      if(present(rhod )) rhod  = reshape(rhodIn, shape(q_mass))
      if(present(vol ))  vol   = reshape(volIn, shape(q_mass))
      if(present(area )) area  = reshape(areaIn, shape(q_mass))
      if(present(refr )) refr  = reshape(refrIn, shape(q_mass))
      if(present(refi )) refi  = reshape(refiIn,shape(q_mass))
