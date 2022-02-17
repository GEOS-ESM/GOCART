      if(present(bext) .or. present(tau) .or. present(ssa) ) then
         bextIn =   this%bext(irh  ,ch,bin) * (1.-arh) &
                  + this%bext(irhp1,ch,bin) * arh
      endif

      if(present(bsca) .or. present(ssa) ) then
         bscaIn =   this%bsca(irh  ,ch,bin) * (1.-arh) &
                  + this%bsca(irhp1,ch,bin) * arh
      endif

      if(present(bbck)) then
         bbckIn =   this%bbck(irh  ,ch,bin) * (1.-arh) &
                  + this%bbck(irhp1,ch,bin) * arh
      endif

      if(present(gasym)) then
         gasymIn =  this%g(irh  ,ch,bin) * (1.-arh) &
                  + this%g(irhp1,ch,bin) * arh
      endif

      if(present(rEff) ) then
         rEffin =     this%rEff(irh  ,bin) * (1.-arh) &
                  + this%rEff(irhp1,bin) * arh
      endif

      if(present(p11) ) then
         p11In =   this%pback(irh  ,ch,bin,1) * (1.-arh) &
                 + this%pback(irhp1,ch,bin,1) * arh
      endif

      if(present(p22) ) then
         p22In =   this%pback(irh  ,ch,bin,5) * (1.-arh) &
                 + this%pback(irhp1,ch,bin,5) * arh
      endif

      if(present(gf) ) then
         gfIn =     this%gf(irh  ,bin) * (1.-arh) &
                  + this%gf(irhp1,bin) * arh
      endif

      if(present(rhod) ) then
         rhodIn =   this%rhod(1  ,bin)
      endif

      if(present(vol) ) then
         volIn  =   this%vol(irh  ,bin) * (1.-arh) &
                  + this%vol(irhp1,bin) * arh
      endif

      if(present(area) ) then
         areaIn  =   this%area(irh  ,bin) * (1.-arh) &
                  + this%area(irhp1,bin) * arh
      endif

      if(present(refr) .or. present(tau) .or. present(ssa) ) then
         refrIn =   this%refr(irh  ,ch,bin) * (1.-arh) &
                  + this%refr(irhp1,ch,bin) * arh
      endif

      if(present(refi) .or. present(tau) .or. present(ssa) ) then
         refiIn =   this%refi(irh  ,ch,bin) * (1.-arh) &
                  + this%refi(irhp1,ch,bin) * arh
      endif
