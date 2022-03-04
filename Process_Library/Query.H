#if (BYWAVE_ == 1)
#define ARG_ wavelength
#endif
#if (BYCHANNEL_ == 1)
#define ARG_ ch
#endif

#if (RANK_ == 0)

#if (BYWAVE_ == 1)
#define SUB_ QueryByWavelength_0d
#endif

#if (BYCHANNEL_ == 1)
#define SUB_ QueryByChannel_0d
#endif

#define DIMENSIONS_ 
#define PMOMDIMENSIONS_ (:,:)
#define DECLAREDIM_ 
#define DECLAREPMOMDIM_ ,allocatable,dimention(:,:)

#endif

#if (RANK_ == 1)

#if (BYWAVE_ == 1)
#define SUB_ QueryByWavelength_1d
#endif

#if (BYCHANNEL_ == 1)
#define SUB_ QueryByChannel_1d
#endif

#define DIMENSIONS_ (:)
#define PMOMDIMENSIONS_ (:,:,:)
#define DECLAREDIM_ ,allocatable,dimension(:) 
#define DECLAREPMOMDIM_ ,allocatable,dimention(:,:,:)

#endif

#if (RANK_ == 2)

#if (BYWAVE_ == 1)
#define SUB_ QueryByWavelength_2d
#endif

#if (BYCHANNEL_ == 1)
#define SUB_ QueryByChannel_2d
#endif

#define DIMENSIONS_ (:,:)
#define PMOMDIMENSIONS_ (:,:,:,:)
#define DECLAREDIM_ ,allocatable,dimension(:) 
#define DECLAREPMOMDIM_ ,allocatable,dimention(:,:,:)
#endif

#if (RANK_ == 3)

#if (BYWAVE_ == 1)
#define SUB_ QueryByWavelength_3d
#endif

#if (BYCHANNEL_ == 1)
#define SUB_ QueryByChannel_3d
#endif

#define DIMENSIONS_ (:,:)
#define PMOMDIMENSIONS_ (:,:,:,:)
#define DECLAREDIM_ ,allocatable,dimension(:) 
#define DECLAREPMOMDIM_ ,allocatable,dimention(:,:,:)
#endif

 subroutine SUB_ ( this, ARG_ , bin, q_mass, rh,   &
                   tau, ssa, gasym, bext, bsca, bbck,   &
                   reff,pmom, p11, p22, gf, rhop, rhod, &
                   vol, area, refr, refi, rc )

! !INPUT PARAMETERS:

     class (GOCART2G_Mie),   intent(in ) :: this
#if (BYWAVE_ ==1)
     real,                   intent(in ) :: wavelength ! wave length [m] 
#endif
#if (BYCHANNEL_ ==1)
     integer,                intent(in ) :: ch         ! wave length [m] 
#endif

     integer,                intent(in ) :: bin        ! bin number
     real,                   intent(in ) :: q_mass DIMENSIONS_     ! aerosol mass [kg/m2],
     real,                   intent(in ) :: rh DIMENSIONS_         ! relative himidity

! !OUTPUT PARAMETERS:

     real,    optional,      intent(out) :: tau DIMENSIONS_        ! aerol extinction optical depth
     real,    optional,      intent(out) :: ssa DIMENSIONS_        ! single scattering albedo
     real,    optional,      intent(out) :: gasym DIMENSIONS_      ! asymmetry parameter
     real,    optional,      intent(out) :: bext DIMENSIONS_       ! mass extinction efficiency [m2 (kg dry mass)-1]
     real,    optional,      intent(out) :: bsca DIMENSIONS_       ! mass scattering efficiency [m2 (kg dry mass)-1]
     real,    optional,      intent(out) :: bbck DIMENSIONS_       ! mass backscatter efficiency [m2 (kg dry mass)-1]
     real,    optional,      intent(out) :: reff DIMENSIONS_       ! effective radius (micron)
     real,    optional,      intent(out) :: pmom PMOMDIMENSIONS_   ! moments of phase function 
     real,    optional,      intent(out) :: p11 DIMENSIONS_        ! P11 phase function at backscatter
     real,    optional,      intent(out) :: p22 DIMENSIONS_        ! P22 phase function at backscatter
     real,    optional,      intent(out) :: gf DIMENSIONS_         ! Growth factor (ratio of wet to dry radius)
     real,    optional,      intent(out) :: rhop DIMENSIONS_       ! Wet particle density [kg m-3]
     real,    optional,      intent(out) :: rhod DIMENSIONS_       ! Dry particle density [kg m-3]
     real,    optional,      intent(out) :: vol DIMENSIONS_        ! Wet particle volume [m3 kg-1]
     real,    optional,      intent(out) :: area DIMENSIONS_       ! Wet particle cross section [m2 kg-1]
     real,    optional,      intent(out) :: refr DIMENSIONS_       ! Wet particle real part of ref. index
     real,    optional,      intent(out) :: refi DIMENSIONS_       ! Wet particle imag. part of ref. index
     integer, optional,      intent(out) :: rc         ! error code


     integer                        :: status
#if (BYWAVE_ == 1)
     integer                        :: ch
#endif

     integer  DECLAREDIM_           :: irh, irhp1, isnap
     real     DECLAREDIM_           :: rh_, arh
     real     DECLAREDIM_           :: bextIn, bscaIn, bbckIn, gasymIn, p11In, p22In, &
                                       gfIn, rhopIn, rhodIn, volIn, areaIn, &
                                       refrIn, refiIn, rEffIn
     real  DECLAREPMOMDIM_          :: pmomIn
     integer :: i, j

     character(len=*), parameter  :: Iam = 'SUB_'

     if ( present(rc) ) rc = 0

#if (BYWAVE_ == 1)
    ch = this%getChannel(wavelength, rc=rc)
    if (present(rc))
      if (rc /=0) return 
    endif
#endif

#if (RANK_ == 0)
      rh_ = max(min(rh,0.99),0.)
      isnap = int((rh_+0.001)*1000.)
      where(isnap .lt. 1) isnap = 1
      arh   = this%rha( isnap )
      irh   = this%rhi( isnap )
      irhp1 = irh+1
      if(irhp1 .gt. this%nrh) irhp1 = this%nrh
#else
      rh_ = reshape(max(min(rh,0.99),0.), [product(shape(rh))])
      isnap = int((rh_+0.001)*1000.)
      where(isnap .lt. 1) isnap = 1
      arh   = this%rha( isnap )
      irh   = this%rhi( isnap )
      irhp1 = irh+1
      where(irhp1 .gt. this%nrh) irhp1 = this%nrh
#endif 

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
#if (RANK_ == 0)
      if (present(pmom)) then
        pmomIn  = this%pmom(irh  ,ch, bin,:,:) * (1.-arh) &
                + this%pmom(irhp1,ch, bin,:,:) * arh
      endif
#else
      if (present(pmom)) then
         allocate(pmomIn(size(irh),size(this%pmom,4),size(this%pmom,5)))
         do j = 1, size(this%pmom,5)
            do i = 1, size(this%pmom,4)
               pmomIn(:,i,j)  =  this%pmom(irh   ,ch,bin,i,j) * (1.-arh) &
                               + this%pmom(irhp1 ,ch,bin,i,j) * arh
            enddo
         enddo
      endif
#endif

#if (RANK_ == 0)

     if(present(tau  )) tau   = bextIn * q_mass
     if(present(ssa  )) ssa   = bscaIn/bextIn
     if(present(bext )) bext  = bextIn
     if(present(bsca )) bsca  = bscaIn
     if(present(bbck )) bbck  = bbckIn
     if(present(gasym)) gasym = gasymIn
     if(present(rEff )) rEff  = 1.E6*rEffIn
     if(present(pmom )) pmom  = pmomIn
     if(present(p11  )) p11   = p11In
     if(present(p22  )) p22   = p22In
     if(present(gf   )) gf    = gfIn
     if(present(rhop )) rhop  = rhopIn
     if(present(rhod )) rhod  = rhodIn
     if(present(vol ))  vol   = volIn
     if(present(area )) area  = areaIn
     if(present(refr )) refr  = refrIn
     if(present(refi )) refi  = refiIn
#else
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
#endif

  end subroutine SUB_

#undef RANK_
#undef BYWAVE_
#undef BYCHANNEL_
