   subroutine NIthermo (km, klid, cdt, grav, delp, rhoa, tmpu, rh, fMassHNO3, fMassAir, &
                        SO4, NH3, NO3an1, NH4a, xhno3, &
                        NI_pno3aq, NI_pnh4aq, NI_pnh3aq, rc)


! !USES:
   implicit NONE

! !INPUT PARAMETERS:
   integer, intent(in) :: km    ! total model levels
   integer, intent(in) :: klid   ! index for pressure lid
   real, intent(in)    :: cdt   ! model time step [s]
   real, intent(in)    :: grav  ! gravity [m s-2]
   real, dimension(:,:,:), intent(in)  :: delp  ! pressure thickness [Pa]
   real, dimension(:,:,:), intent(in)  :: rhoa   ! Layer air density [kg m-3]
   real, dimension(:,:,:), intent(in)  :: tmpu   ! Layer temperature [K]
   real, dimension(:,:,:), intent(in)  :: rh     ! relative humidity [0-1]
   real, intent(in)  :: fMassHNO3   ! gram molecular weight of HNO3
   real, intent(in)  :: fMassAir    ! gram molecular weight of air

! !INOUTPUT PARAMETERS:
   real, dimension(:,:,:), intent(inout)  :: SO4    ! Sulphate aerosol [kg kg-1]
   real, dimension(:,:,:), intent(inout)  :: NH3    ! Ammonia (NH3, gas phase) [kg kg-1]
   real, dimension(:,:,:), intent(inout)  :: NO3an1 ! Nitrate size bin 001 [kg kg-1]
   real, dimension(:,:,:), intent(inout)  :: NH4a   ! Ammonium ion (NH4+, aerosol phase) [kg kg-1]
   real, dimension(:,:,:), intent(inout)  :: xhno3  ! buffer for NITRATE_HNO3 [kg m-2 sec-1]
   real, pointer, dimension(:,:), intent(inout) :: NI_pno3aq ! Nitrate Production from Aqueous Chemistry [kg m-2 s-1]
   real, pointer, dimension(:,:), intent(inout) :: NI_pnh4aq ! Ammonium Production from Aqueous Chemistry [kg m-2 s-1]
   real, pointer, dimension(:,:), intent(inout) :: NI_pnh3aq ! Ammonia Change from Aqueous Chemistry [kg m-2 s-1]

! !OUTPUT PARAMETERS:
   integer, optional, intent(out) :: rc                   ! Error return code:


! !DESCRIPTION: Prepares variables and calls the RPMARES (thermodynamics module)
!
! !REVISION HISTORY:
!
! Aug2020 E.Sherman - Refactored for process library
!

! !Local Variables
   real   :: fmmr_to_conc
   real(kind=DP) :: SO4_, GNO3, GNH3, RH_, TEMP, ASO4, AHSO4, AH2O, ANO3, ANH4
   integer :: k, j, i

   integer :: status

!EOP
!-------------------------------------------------------------------------
!  Begin...

   do k = klid, km
    do j = 1, ubound(tmpu,2)
     do i = 1, ubound(tmpu,1)

!     Conversion of mass mixing ratio to concentration (ug m-3)
      fmmr_to_conc = 1.e9 * rhoa(i,j,k)

!     Unit conversion for input to thermodynamic module
!     Per grid box call to RPMARES thermodynamic module
!     We do not presently treat chemistry of sulfate completely,
!     hence we ignore terms for ASO4, AHSO4, AH2O, and we do
!     not update SO4 on output from RPMARES.
!     At present we are importing HNO3 from offline file, so we
!     do not update on return.
      SO4_  = 1.d-32
      SO4_  = max(1.d-32,SO4(i,j,k) * fmmr_to_conc)
      GNO3  = max(1.d-32,xhno3(i,j,k) * fMassHNO3 / fMassAir * fmmr_to_conc)
      GNH3  = max(1.d-32,NH3(i,j,k)  * fmmr_to_conc)
      RH_   = rh(i,j,k)
      TEMP  = tmpu(i,j,k)
      ASO4  = 1.d-32
      AHSO4 = 1.d-32
      ANO3  = max(1.d-32,NO3an1(i,j,k) * fmmr_to_conc)
      AH2O  = 1.d-32
      ANH4  = max(1.d-32,NH4a(i,j,k) * fmmr_to_conc)

!print*,'GOCART2G NIthermo TEST 2'

      call RPMARES (  SO4_, GNO3,  GNH3, RH_,  TEMP, &
                      ASO4, AHSO4, ANO3, AH2O, ANH4, __RC__ )

!     Diagnostic terms
      if(associated(NI_pno3aq)) &
       NI_pno3aq(i,j) = NI_pno3aq(i,j) &
        + (ANO3 / fmmr_to_conc - NO3an1(i,j,k)) &
          * delp(i,j,k)/grav/cdt
      if(associated(NI_pnh4aq)) &
       NI_pnh4aq(i,j) = NI_pnh4aq(i,j) &
        + (ANH4 / fmmr_to_conc - NH4a(i,j,k)) &
          * delp(i,j,k)/grav/cdt
      if(associated(NI_pnh3aq)) &
       NI_pnh3aq(i,j) = NI_pnh3aq(i,j) &
        + (GNH3 / fmmr_to_conc - NH3(i,j,k)) &
          * delp(i,j,k)/grav/cdt

!     Unit conversion back on return from thermodynamic module
      NH3(i,j,k)    = GNH3 / fmmr_to_conc
      NO3an1(i,j,k) = ANO3 / fmmr_to_conc
      NH4a(i,j,k)   = ANH4 / fmmr_to_conc
      xhno3(i,j,k) = max(1.d-32, GNO3 * fMassAir / fMassHNO3 / fmmr_to_conc)

     enddo
    enddo
   enddo


   __RETURN__(__SUCCESS__)
   end subroutine NIthermo
