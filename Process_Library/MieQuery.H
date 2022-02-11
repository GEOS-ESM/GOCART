      if(present(bext) .or. present(tau) .or. present(ssa) ) then
         bextIn =   this%bext(irh  ,ichannel,TYPE) * (1.-arh) &
                  + this%bext(irhp1,ichannel,TYPE) * arh
      endif

      if(present(bsca) .or. present(ssa) ) then
         bscaIn =   this%bsca(irh  ,ichannel,TYPE) * (1.-arh) &
                  + this%bsca(irhp1,ichannel,TYPE) * arh
      endif

      if(present(bbck)) then
         bbckIn =   this%bbck(irh  ,ichannel,TYPE) * (1.-arh) &
                  + this%bbck(irhp1,ichannel,TYPE) * arh
      endif

      if(present(gasym)) then
         gasymIn =  this%g(irh  ,ichannel,TYPE) * (1.-arh) &
                  + this%g(irhp1,ichannel,TYPE) * arh
      endif

      if(present(rEff) ) then
         rEffin =     this%rEff(irh  ,TYPE) * (1.-arh) &
                  + this%rEff(irhp1,TYPE) * arh
      endif

!      if(present(pmom)) then
!         pmom(:,:) = TABLE%pmom(irh  ,ichannel,TYPE,:,:) * (1.-arh) &
!                   + TABLE%pmom(irhp1,ichannel,TYPE,:,:) * arh
!      endif

!      if(present(pmom)) then
!         call Chem_MieQueryByIntWithpmom(this, idx, channel, q_mass, rh, pmom)
!      endif


      if(present(p11) ) then
         p11In =   this%pback(irh  ,ichannel,TYPE,1) * (1.-arh) &
                 + this%pback(irhp1,ichannel,TYPE,1) * arh
      endif

      if(present(p22) ) then
         p22In =   this%pback(irh  ,ichannel,TYPE,5) * (1.-arh) &
                 + this%pback(irhp1,ichannel,TYPE,5) * arh
      endif

      if(present(gf) ) then
         gfIn =     this%gf(irh  ,TYPE) * (1.-arh) &
                  + this%gf(irhp1,TYPE) * arh
      endif

      if(present(rhod) ) then
         rhodIn =   this%rhod(1  ,TYPE)
      endif

      if(present(vol) ) then
         volIn  =   this%vol(irh  ,TYPE) * (1.-arh) &
                  + this%vol(irhp1,TYPE) * arh
      endif

      if(present(area) ) then
         areaIn  =   this%area(irh  ,TYPE) * (1.-arh) &
                  + this%area(irhp1,TYPE) * arh
      endif

      if(present(refr) .or. present(tau) .or. present(ssa) ) then
         refrIn =   this%refr(irh  ,ichannel,TYPE) * (1.-arh) &
                  + this%refr(irhp1,ichannel,TYPE) * arh
      endif

      if(present(refi) .or. present(tau) .or. present(ssa) ) then
         refiIn =   this%refi(irh  ,ichannel,TYPE) * (1.-arh) &
                  + this%refi(irhp1,ichannel,TYPE) * arh
      endif
