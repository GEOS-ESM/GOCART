      integer                      :: ICHANNEL, TYPE
      integer, allocatable, dimension(:)   :: irh, irhp1, isnap
      real, allocatable, dimension(:)      :: rhUse, arh
      real, allocatable, dimension(:)      :: bextIn, bscaIn, bbckIn, gasymIn, p11In, p22In, &
                                              gfIn, rhopIn, rhodIn, volIn, areaIn, &
                                              refrIn, refiIn, rEffIn

      if ( present(rc) ) rc = 0

      ICHANNEL = nint(CHANNEL)
      TYPE = idx

!     Now map the input RH to the high resolution hash table for RH

      rhUse = reshape(max(rh,0.), [product(shape(rh))])
      rhUse = reshape(min(rh,0.99), [product(shape(rh))])
      isnap = int((rhUse+0.001)*1000.)
      where(isnap .lt. 1) isnap = 1
      arh   = this%rha( isnap )
      irh   = this%rhi( isnap )
      irhp1 = irh+1
      where(irhp1 .gt. this%nrh) irhp1 = this%nrh

!     Now linearly interpolate the input table for the requested aerosol and
!     channel; rh is the relative humidity.

include "MieQuery.H"

!     Fill the requested outputs
      if(present(tau  )) tau   = reshape(bextIn,shape(q_mass)) * q_mass
      if(present(ssa  )) ssa   = reshape(bscaIn/bextIn,shape(q_mass))
      if(present(bext )) bext  = reshape(bextIn,shape(q_mass))
      if(present(bsca )) bsca  = reshape(bscaIn,shape(q_mass))
      if(present(bbck )) bbck  = reshape(bbckIn,shape(q_mass))
      if(present(gasym)) gasym = reshape(gasymIn,shape(q_mass))
      if(present(rEff )) rEff  = reshape(1.E6*rEffIn,shape(q_mass)) ! convert to microns
      if(present(p11  )) p11   = reshape(p11In,shape(q_mass))
      if(present(p22  )) p22   = reshape(p22In, shape(q_mass))
      if(present(gf   )) gf    = reshape(gfIn, shape(q_mass))
      if(present(rhop )) rhop  = reshape(rhopIn, shape(q_mass))
      if(present(rhod )) rhod  = reshape(rhodIn, shape(q_mass))
      if(present(vol ))  vol   = reshape(volIn, shape(q_mass))
      if(present(area )) area  = reshape(areaIn, shape(q_mass))
      if(present(refr )) refr  = reshape(refrIn, shape(q_mass))
      if(present(refi )) refi  = reshape(refiIn,shape(q_mass))
