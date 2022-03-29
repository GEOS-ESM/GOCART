   subroutine DMSemission (km, cdt, grav, tmpu, u10m, v10m, oro, delp, &
                           fMassDMS, dmso_conc, dms, SU_emis, ndms, rc)

! !USES:
   implicit NONE

! !INPUT PARAMETERS:
   integer, intent(in) :: km  ! number model layers, and number of species respectively
   real, intent(in)    :: cdt ! model time step [seconds]
   real, intent(in)    :: grav ! gravity [m sec-1]
   real, pointer, dimension(:,:,:), intent(in)  :: tmpu  ! temperature [K]
   real, pointer, dimension(:,:), intent(in)    :: u10m  ! 10-m u-wind component [m s-1]
   real, pointer, dimension(:,:), intent(in)    :: v10m  ! 10-m v-wind component [m s-1]
   real, pointer, dimension(:,:), intent(in)    :: oro   ! orography flag
   real, pointer, dimension(:,:,:), intent(in)  :: delp  ! pressure thickness [Pa]
   real, dimension(:,:), intent(in) :: dmso_conc ! DMS source [1]
   integer, intent(in) :: ndms      ! index of DMS relative to other sulfate tracers
   real, intent(in)    :: fMassDMS  ! gram molecular weight of DMS


! !INOUT PARAMETERS:
   real, dimension(:,:,:), intent(inout)  :: dms ! dms [kg kg-1]
   real, pointer, dimension(:,:,:), intent(inout)  :: SU_emis   ! SU emissions, kg/m2/s

! !OUTPUT PARAMETERS:
   integer, optional, intent(out)   :: rc    ! Error return code:
                                             !  0 - all is well

! !DESCRIPTION: Computes the dust emissions for one time step
!
! !REVISION HISTORY:
!
! 11Feb2020 E.Sherman - First attempt at refactor
!

! !Local Variables
   integer         :: i, j, k
   integer         :: i1=1, j1=1, i2, j2

   real, dimension(:,:), allocatable :: srcDMS
   real :: sCO2, schmidt, w10m, akw, sst

!EOP
!-------------------------------------------------------------------------
!  Begin

!  Add in the DMS source
!  ---------------------
!  DMS emissions go into the lowest model layer only
!  The transfer of DMS from the ocean surface to the atmosphere is
!  a function of surface temperature and wind speed.
!  For now we use the lowest atmospheric temperature (really want SST)
!  and the 10-m wind speed.
!  This code follows from GOCART with the following notes:
!  :the Schmidt number for CO2 is assumed to be 600
!  :the Schmidt number of DMSo follows Saltzman et al., 1993
!  :the Schmidt number dependence breaks for high SST
!  :following www.knmi.nl/~velthove/TM/input we introduce a maximum
!   temperature of 28 C for the calculation
!  :the w10m dependence is from Liss and Merlivat (1986)
!  All this needs some thorough checking!

    i2 = size(tmpu,1)
    j2 = size(tmpu,2)

    allocate(srcDMS(i2,j2))
    srcDMS = 0.

    k = km
    sCO2 = 600.
    do j = j1, j2
       do i = i1, i2
          sst = tmpu(i,j,k)-273.15
          if(sst .gt. 28.) sst = 28.
!         only valid for ocean and warm enough temperatures
          if( (oro(i,j) /= OCEAN) .or. (sst .lt. -20.)) cycle
            schmidt = 2764.0 - 147.12*sst + 3.726*(sst**2.) - 0.038*(sst**3.)
!           w10m is the 10-m wind speed in m s-1
            w10m = sqrt(u10m(i,j)**2. + v10m(i,j)**2.)
          if(w10m .le. 3.6) then
             akw = 0.17*w10m*((sCO2/schmidt)**0.667)
          else if (w10m .le. 13.) then
             akw = (2.85*w10m - 9.65)*sqrt(sCO2/schmidt)
          else
             akw = (5.90*w10m - 49.3)*sqrt(sCO2/schmidt)
          endif
!         This parameterization has put akw in units cm hr-1 -> goto m s-1
          akw = akw/100./3600.
!         DMSo concentration is nMol/L
!         Want to put the source into units of kg m-2 s-1
          srcDMS(i,j) = akw * (fmassDMS/1000.)*(dmso_conc(i,j)*1.e-9/1.e-3)
          dms(i,j,k) =  dms(i,j,k) + srcDMS(i,j)*cdt*grav/delp(i,j,k)
       end do
    end do

    if( associated(SU_emis )) SU_emis(:,:,ndms) = srcDMS


      __RETURN__(__SUCCESS__)
   end subroutine DMSemission
