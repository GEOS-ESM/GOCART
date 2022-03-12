#include "overload.macro"

!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  GOCART2G_Query --- Return Tau, SSA, etc 
!
!
! !INTERFACE:
!

 subroutine SUB_ ( this, ARG_ , bin, q_mass, rh,   &
                   tau, ssa, gasym, bext, bsca, bbck,   &
                   reff,pmom, p11, p22, gf, rhop, rhod, &
                   vol, area, refr, refi, rc )

! !INPUT PARAMETERS:

     class (GOCART2G_Mie),   intent(in ) :: this
#if defined(BYWAVE_)
     real,                   intent(in ) :: wavelength ! wave length [m] 
#endif
#if defined(BYCHANNEL_)
     integer,                intent(in ) :: channel    ! wave length [m] 
#endif
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
#if defined(BYWAVE_)
     integer                        :: channel
#endif

     integer, allocatable, dimension(__DIMS__)  :: irh
     integer, allocatable, dimension(:)         :: qrh
     real, allocatable, dimension(__DIMS__)     :: rh_, arh
     real, allocatable, dimension(__DIMS__)     :: bext_, bsca_
     integer :: i, j

     character(len=*), parameter  :: Iam = 'Query'

     if ( present(rc) ) rc = 0

#if defined(BYWAVE_)
    channel = this%getChannel(wavelength, rc=rc)
    if (present(rc)) then
      if (rc /=0) return 
    endif
#endif

    rh_ = min(max(rh,0.),0.99)                                  ! no (super) saturation 
    qrh = reshape(max(1,nint((rh_+0.001)*1000.)), [size(rh_)])  ! quantized RH
    arh = reshape(this%rha(qrh), shape(rh_))                    ! interpolation weight in table
    irh = reshape(this%rhi(qrh), shape(rh_))                    ! index in table for each RH vcalue
    where (irh==this%nrh )
       irh = this%nrh-1
       arh = 1.0
    end where

    if(present(bext) .or. present(tau) .or. present(ssa) ) then
       bext_ =  this%interp(this%bext(channel,bin), irh, arh)
    endif

    if(present(bsca) .or. present(ssa) ) then
       bsca_ =  this%interp(this%bsca(channel,bin), irh, arh)
    endif

    if(present(bbck)) then
       bbck =  this%interp(this%bbck(channel,bin), irh, arh)
    endif

    if(present(gasym)) then
       gasym =  this%interp(this%g(channel,bin), irh, arh)
    endif

    if(present(rEff)) then
       rEff = 1.E6 * this%interp(this%rEff(bin), irh, arh)
    endif

    if(present(p11)) then
       p11 =  this%interp(this%pback(channel,bin,1), irh, arh)
    endif

    if(present(p22)) then
       p22 =  this%interp(this%pback(channel,bin,5), irh, arh)
    endif

    if(present(gf)) then
       gf =  this%interp(this%gf(bin), irh, arh)
    endif

    if(present(rhod)) then
       rhod = this%interp(this%rhod(bin), irh, arh)
    endif

    if(present(vol)) then
       vol = this%interp(this%vol(bin), irh, arh)
    endif

    if(present(area)) then
       area = this%interp(this%area(bin), irh, arh)
    endif

    if(present(refr)) then
       refr = this%interp(this%refr(channel,bin), irh, arh)
    endif

    if(present(refi)) then
       refi = this%interp(this%refi(channel,bin), irh, arh)
    endif

    if (present(pmom)) then
       do j = 1, size(this%pmom,4)
          do i = 1, size(this%pmom,3)
             pmom(__DIMS__,i,j) = this%interp(this%pmom(channel,bin,i,j),irh,arh)
          enddo
       enddo
    endif

    if(present(tau  )) tau   = bext_ * q_mass
    if(present(ssa  )) ssa   = bsca_/bext_
    if(present(bext )) bext  = bext_
    if(present(bsca )) bsca  = bsca_

  end subroutine SUB_

#include "undef_overload.macro"
