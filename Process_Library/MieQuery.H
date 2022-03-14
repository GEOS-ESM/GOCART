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

     integer, allocatable, dimension(:)       :: irh, qrh
     real, allocatable, dimension(:)          :: rh_, arh
     real, allocatable, dimension(__DIMS__)   :: bext_, bsca_
     integer :: i, j

     character(len=*), parameter  :: Iam = 'Query'

     if ( present(rc) ) rc = 0

#if defined(BYWAVE_)
    channel = this%getChannel(wavelength, rc=rc)
    if (present(rc)) then
      if (rc /=0) return 
    endif
#endif

    rh_ = reshape(min(max(rh,0.),0.99),[size(rh)])    ! no (super) saturation 
    qrh = max(1,nint((rh_+0.001)*1000.))              ! quantized RH
    arh = this%rha(qrh)                               ! interpolation weight in table
    irh = this%rhi(qrh)                               ! index in table for each RH vcalue
    where (irh==this%nrh )
       irh = this%nrh-1
       arh = 1.0
    end where

    if(present(bext) .or. present(tau) .or. present(ssa) ) then
       bext_ =  __INTERP__(this%bext, irh, arh, channel, bin, bext)
    endif

    if(present(bsca) .or. present(ssa) ) then
       bsca_ =  __INTERP__(this%bsca, irh, arh, channel, bin, bsca)
    endif

    if(present(bbck)) then
       bbck =  __INTERP__(this%bbck, irh, arh, channel, bin, bbck)
    endif

    if(present(gasym)) then
       gasym =  __INTERP__(this%g, irh, arh, channel, bin, gasym)
    endif

    if(present(rEff)) then
       rEff = 1.E6* __INTERP1__(this%rEff, irh, arh, bin, rEff)
    endif

    if(present(p11)) then
       p11 =  __INTERPBACK__(this%pback, irh, arh, channel, bin, 1, p11)
    endif

    if(present(p22)) then
       p22 =  __INTERPBACK__(this%pback, irh, arh, channel, bin, 5, p22)
    endif

    if(present(gf)) then
       gf =  __INTERP1__(this%gf, irh, arh, bin, gf)
    endif

    if(present(rhod)) then
       rhod =  __INTERP1__(this%rhod, irh, arh, bin, rhod)
    endif

    if(present(vol)) then
       vol =  __INTERP1__(this%vol, irh, arh, bin, vol)
    endif

    if(present(area)) then
       area =  __INTERP1__(this%area, irh, arh, bin, area)
    endif

    if(present(refr)) then
       refr =  __INTERP__(this%refr, irh, arh, channel, bin, refr)
    endif

    if(present(refi)) then
       refi =  __INTERP__(this%refi, irh, arh, channel, bin, refi)
    endif

    if (present(pmom)) then
       do j = 1, size(this%pmom,5)
          do i = 1, size(this%pmom,4)
             pmom(__DIMS__,i,j)  =  __INTERPMOM__(this%pmom, irh, arh, channel, bin, i, j, rh)
          enddo
       enddo
    endif

    if(present(tau  )) tau   = bext_ * q_mass
    if(present(ssa  )) ssa   = bsca_/bext_
    if(present(bext )) bext  = bext_
    if(present(bsca )) bsca  = bsca_

  end subroutine SUB_

#include "undef_overload.macro"
