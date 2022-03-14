#define IDENTITY_(a) a

#if (RANK_ == 1)
#  define SUFIX_ _1d
#  define __DIMS__ : 
#elif (RANK_ == 2)
#  define SUFIX_ _2d
#  define __DIMS__ :,: 
#elif (RANK_ == 3)
#  define SUFIX_ _3d
#  define __DIMS__ :,:,: 
#endif

#define __INTERP2__(table,rh,ah,bin,sh) reshape((table(rh,bin) * (1-ah) + table(rh+1,bin)*ah), shape(sh))
#define __INTERP3__(table,rh,ah,ch,bin,sh) reshape((table(rh,ch,bin) * (1-ah) + table(rh+1,ch,bin)*ah), shape(sh))
#define __INTERP5__(table,rh,ah,ch,bin,k,l,sh) reshape((table(rh,ch,bin,k,l) * (1-ah) + table(rh+1,ch,bin,k,l)*ah), shape(sh))

#define SUBCHANNEL_ IDENTITY_(QueryByChannel)IDENTITY_(SUFIX_)
#define SUBWAVELENGTH_ IDENTITY_(QueryByWavelength)IDENTITY_(SUFIX_)

!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  GOCART2G_Query --- Return Tau, SSA, etc 
!
!
! !INTERFACE:
!

 subroutine SUBCHANNEL_ ( this, channel , bin, q_mass, rh,   &
                   tau, ssa, gasym, bext, bsca, bbck,   &
                   reff, pmom, p11, p22, gf, rhop, rhod, &
                   vol, area, refr, refi, rc )

! !INPUT PARAMETERS:

     class (GOCART2G_Mie),   intent(in ) :: this
     integer,                intent(in ) :: channel    ! wave length [m] 
     integer,                intent(in ) :: bin        ! bin number
     real,                   intent(in ) :: q_mass(__DIMS__)      ! aerosol mass [kg/m2],
     real,                   intent(in ) :: rh    (__DIMS__)      ! relative himidity

! !OUTPUT PARAMETERS:

     real,    optional,      intent(out) :: tau    (__DIMS__)      ! aerol extinction optical depth
     real,    optional,      intent(out) :: ssa    (__DIMS__)      ! single scattering albedo
     real,    optional,      intent(out) :: gasym  (__DIMS__)      ! asymmetry parameter
     real,    optional,      intent(out) :: bext   (__DIMS__)      ! mass extinction efficiency [m2 (kg dry mass)-1]
     real,    optional,      intent(out) :: bsca   (__DIMS__)      ! mass scattering efficiency [m2 (kg dry mass)-1]
     real,    optional,      intent(out) :: bbck   (__DIMS__)      ! mass backscatter efficiency [m2 (kg dry mass)-1]
     real,    optional,      intent(out) :: reff   (__DIMS__)      ! effective radius (micron)
     real,    optional,      intent(out) :: pmom   (__DIMS__,:,:)  ! moments of phase function 
     real,    optional,      intent(out) :: p11    (__DIMS__)      ! P11 phase function at backscatter
     real,    optional,      intent(out) :: p22    (__DIMS__)      ! P22 phase function at backscatter
     real,    optional,      intent(out) :: gf     (__DIMS__)      ! Growth factor (ratio of wet to dry radius)
     real,    optional,      intent(out) :: rhop   (__DIMS__)      ! Wet particle density [kg m-3]
     real,    optional,      intent(out) :: rhod   (__DIMS__)      ! Dry particle density [kg m-3]
     real,    optional,      intent(out) :: vol    (__DIMS__)      ! Wet particle volume [m3 kg-1]
     real,    optional,      intent(out) :: area   (__DIMS__)      ! Wet particle cross section [m2 kg-1]
     real,    optional,      intent(out) :: refr   (__DIMS__)      ! Wet particle real part of ref. index
     real,    optional,      intent(out) :: refi   (__DIMS__)      ! Wet particle imag. part of ref. index
     integer, optional,      intent(out) :: rc                      ! error code

! !DESCRIPTION:
!
!   Returns requested parameters from the Mie tables, as a function 
!   of species, relative humidity, and channel
!
!  Notes: Needs some checking, and I still force an interpolation step

!
! !REVISION HISTORY:
!
!  03Mar2022 Tom Clune refactoring
!  23Mar2005 Colarco
!  11Jul2005 da Silva   Standardization.
!
!EOP
!-------------------------------------------------------------------------

     integer                        :: status

     integer, allocatable, dimension(:)       :: irh, qrh
     real, allocatable, dimension(:)          :: rh_, arh
     real, allocatable, dimension(__DIMS__)   :: bext_, bsca_
     integer :: i, j

     character(len=*), parameter  :: Iam = 'Query'

     if ( present(rc) ) rc = 0

     rh_ = reshape(min(max(rh,0.),0.99),[size(rh)])    ! no (super) saturation 
     qrh = max(1,nint((rh_+0.001)*1000.))              ! quantized RH
     arh = this%rha(qrh)                               ! interpolation weight in table
     irh = this%rhi(qrh)                               ! index in table for each RH vcalue
     where (irh==this%nrh )
       irh = this%nrh-1
       arh = 1.0
     end where

     if(present(bext) .or. present(tau) .or. present(ssa) ) then
       bext_ =  __INTERP3__(this%bext, irh, arh, channel, bin, bext)
     endif

     if(present(bsca) .or. present(ssa) ) then
       bsca_ =  __INTERP3__(this%bsca, irh, arh, channel, bin, bsca)
     endif

     if(present(bbck)) then
       bbck =  __INTERP3__(this%bbck, irh, arh, channel, bin, bbck)
     endif

     if(present(gasym)) then
       gasym =  __INTERP3__(this%g, irh, arh, channel, bin, gasym)
     endif

     if(present(rEff)) then
       rEff = 1.E6* __INTERP2__(this%rEff, irh, arh, bin, rEff)
     endif

     if(present(p11)) then
       p11 =  __INTERP3__(this%p11, irh, arh, channel, bin, p11)
     endif

     if(present(p22)) then
       p22 =  __INTERP3__(this%p22, irh, arh, channel, bin, p22)
     endif

     if(present(gf)) then
       gf =  __INTERP2__(this%gf, irh, arh, bin, gf)
     endif

     if(present(rhod)) then
       rhod =  __INTERP2__(this%rhod, irh, arh, bin, rhod)
     endif

     if(present(vol)) then
       vol =  __INTERP2__(this%vol, irh, arh, bin, vol)
     endif

     if(present(area)) then
       area =  __INTERP2__(this%area, irh, arh, bin, area)
     endif

     if(present(refr)) then
       refr =  __INTERP3__(this%refr, irh, arh, channel, bin, refr)
     endif

     if(present(refi)) then
       refi =  __INTERP3__(this%refi, irh, arh, channel, bin, refi)
     endif

     if (present(pmom)) then
       do j = 1, size(this%pmom,5)
          do i = 1, size(this%pmom,4)
             pmom(__DIMS__,i,j)  =  __INTERP5__(this%pmom, irh, arh, channel, bin, i, j, rh)
          enddo
       enddo
     endif

     if(present(tau  )) tau   = bext_ * q_mass
     if(present(ssa  )) ssa   = bsca_/bext_
     if(present(bext )) bext  = bext_
     if(present(bsca )) bsca  = bsca_

  end subroutine SUBCHANNEL_


!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  GOCART2G_Query --- Return Tau, SSA, etc
!
!
! !INTERFACE:
!

 subroutine SUBWAVELENGTH_ ( this, wavelength, bin, q_mass, rh,   &
                   tau, ssa, gasym, bext, bsca, bbck,   &
                   reff, pmom, p11, p22, gf, rhop, rhod, &
                   vol, area, refr, refi, rc )

! !INPUT PARAMETERS:

     class (GOCART2G_Mie),   intent(in ) :: this
     real,                   intent(in ) :: wavelength ! wave length [m]
     integer,                intent(in ) :: bin        ! bin number
     real,                   intent(in ) :: q_mass(__DIMS__)      ! aerosol mass [kg/m2],
     real,                   intent(in ) :: rh    (__DIMS__)      ! relative himidity

! !OUTPUT PARAMETERS:

     real,    optional,      intent(out) :: tau    (__DIMS__)      ! aerol extinction optical depth
     real,    optional,      intent(out) :: ssa    (__DIMS__)      ! single scattering albedo
     real,    optional,      intent(out) :: gasym  (__DIMS__)      ! asymmetry parameter
     real,    optional,      intent(out) :: bext   (__DIMS__)      ! mass extinction efficiency [m2 (kg dry mass)-1]
     real,    optional,      intent(out) :: bsca   (__DIMS__)      ! mass scattering efficiency [m2 (kg dry mass)-1]
     real,    optional,      intent(out) :: bbck   (__DIMS__)      ! mass backscatter efficiency [m2 (kg dry mass)-1]
     real,    optional,      intent(out) :: reff   (__DIMS__)      ! effective radius (micron)
     real,    optional,      intent(out) :: pmom   (__DIMS__,:,:)  ! moments of phase function
     real,    optional,      intent(out) :: p11    (__DIMS__)      ! P11 phase function at backscatter
     real,    optional,      intent(out) :: p22    (__DIMS__)      ! P22 phase function at backscatter
     real,    optional,      intent(out) :: gf     (__DIMS__)      ! Growth factor (ratio of wet to dry radius)
     real,    optional,      intent(out) :: rhop   (__DIMS__)      ! Wet particle density [kg m-3]
     real,    optional,      intent(out) :: rhod   (__DIMS__)      ! Dry particle density [kg m-3]
     real,    optional,      intent(out) :: vol    (__DIMS__)      ! Wet particle volume [m3 kg-1]
     real,    optional,      intent(out) :: area   (__DIMS__)      ! Wet particle cross section [m2 kg-1]
     real,    optional,      intent(out) :: refr   (__DIMS__)      ! Wet particle real part of ref. index
     real,    optional,      intent(out) :: refi   (__DIMS__)      ! Wet particle imag. part of ref. index
     integer, optional,      intent(out) :: rc                      ! error code

! !DESCRIPTION:
!
!   Returns requested parameters from the Mie tables, as a function
!   of species, relative humidity, and channel
!
!  Notes: Needs some checking, and I still force an interpolation step

!
! !REVISION HISTORY:
!
!  03Mar2022 Tom Clune refactoring
!  23Mar2005 Colarco
!  11Jul2005 da Silva   Standardization.
!
!EOP
!-------------------------------------------------------------------------

     integer                        :: status
     integer                        :: channel

     character(len=*), parameter  :: Iam = 'QueryByWavelength'

     if ( present(rc) ) rc = 0

    channel = this%getChannel(wavelength, rc=rc)
    if (present(rc)) then
      if (rc /=0) return
    endif

    call this%Query(channel, bin, q_mass, rh,                   &
                   tau, ssa, gasym, bext, bsca, bbck,           &
                   reff, pmom, p11, p22, gf, rhop, rhod, &
                   vol, area, refr, refi, rc )

 end subroutine SUBWAVELENGTH_

#undef SUBWAVELENGTH_
#undef SUBCHANNEL_
#undef __INTERP2__
#undef __INTERP3__
#undef __INTERP5__
#undef __DIMS__
#undef SUFIX_
#undef IDENTITY_

