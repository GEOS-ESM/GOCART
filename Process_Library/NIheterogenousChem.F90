   subroutine NIheterogenousChem (NI_phet, xhno3, UNDEF, AVOGAD, AIRMW, PI, RUNIV, rhoa, tmpu, relhum, delp, &
                                  DU, SS, rmedDU, rmedSS, fnumDU, fnumSS, nbinsDU, nbinsSS, &
                                  km, klid, cdt, grav, fMassHNO3, fMassNO3, nNO3an1, nNO3an2, &
                                  nNO3an3, HNO3_conc, HNO3_sfcmass, HNO3_colmass, rc)


! !DESCRIPTION: Nitrogen heterogeneous chemistry - Optimized by A. Darmenov

! !USES:
   implicit NONE

! !INPUT PARAMETERS:
   real, intent(in)                    :: UNDEF          ! set an undefined value (MAPL_UNDEF)
   real, intent(in)                    :: AVOGAD         ! Avogadro's number [1/kmol]
   real, intent(in)                    :: AIRMW          ! molecular weight of air [kg/kmol]
   real, intent(in)                    :: PI             ! pi constant
   real, intent(in)                    :: RUNIV          ! ideal gas constant [J/(Kmole*K)]
   real, dimension(:,:,:), intent(in)  :: rhoa           ! Layer air density [kg/m^3]
   real, dimension(:,:,:), intent(in)  :: tmpu           ! Layer temperature [K]
   real, dimension(:,:,:), intent(in)  :: relhum         ! relative humidity [1]
   real, dimension(:,:,:), intent(in)  :: delp           ! pressure thickness [Pa]
   real, pointer, dimension(:,:,:,:), intent(in) :: DU   ! dust aerosol [kg/kg]
   real, pointer, dimension(:,:,:,:), intent(in) :: SS   ! sea salt aerosol [kg/kg]
   real, dimension(:) ,intent(in)      :: rmedDU         ! dust aerosol radius [um]
   real, dimension(:) ,intent(in)      :: rmedSS         ! sea salt aerosol radius [um]
   real, dimension(:) ,intent(in)      :: fnumDU         ! number of dust particles per kg mass
   real, dimension(:) ,intent(in)      :: fnumSS         ! number of sea salt particles per kg mass
   integer, intent(in)                 :: nbinsDU        ! number of dust bins
   integer, intent(in)                 :: nbinsSS        ! number of sea salt bins
   integer, intent(in)                 :: km             ! number of model levels
   integer, intent(in)                 :: klid           ! index for pressure lid
   real, intent(in)                    :: cdt            ! chemistry model timestep (sec)
   real, intent(in)                    :: grav           ! gravity (m/sec)
   real, intent(in)                    :: fMassHNO3      ! gram molecular weight
   real, intent(in)                    :: fMassNO3       ! gram molecular weight

! !INOUTPUT PARAMETERS:
   real, pointer, dimension(:,:,:), intent(inout)  :: NI_phet   ! Nitrate Production from Het Chem [kg/(m^2 sec)]
   real, dimension(:,:,:), intent(inout)  :: xhno3     ! buffer for NITRATE_HNO3 [kg/(m^2 sec)]
   real, pointer, dimension(:,:,:), intent(inout)  :: HNO3_conc ! Nitric Acid Mass Concentration [kg/m^3]
   real, pointer, dimension(:,:), intent(inout)    :: HNO3_sfcmass ! Nitric Acid Surface Mass Concentration [kg/m^3]
   real, pointer, dimension(:,:), intent(inout)    :: HNO3_colmass ! Nitric Acid Column Mass Density [kg/m^3]
   real, pointer, dimension(:,:,:), intent(inout)  :: nNO3an1 ! Nitrate bin 1 [kg/kg]
   real, pointer, dimension(:,:,:), intent(inout)  :: nNO3an2 ! Nitrate bin 2 [kg/kg]
   real, pointer, dimension(:,:,:), intent(inout)  :: nNO3an3 ! Nitrate bin 3 [kg/kg]

! !OUTPUT PARAMETERS:
   integer, optional, intent(out) :: rc

! !Local Variables
   real, dimension(:,:,:), allocatable :: kan1, kan2, kan3, kan
   real, dimension(:,:,:), allocatable :: deltahno3

   integer :: i1, j1, i2, j2, n, i, j, k


!
! !REVISION HISTORY:
!
! ???? Optimized NI Het Chem - A. Darmenov
! 15Dec2020 - Ported to process library - E. Sherman

!EOP
!------------------------------------------------------------------------------------
!  Begin..

!  Heterogeneous chemistry
!  -----------------------
!  Heterogeneous chemistry wants to know about GOCART dust and sea
!  salt tracers.  This code is not at the moment generalized as it
!  seems very wedded to the traditional GOCART arrangement (5 dust,
!  5 sea salt) and the particulars of the nitrate aerosol arrangement.

   j1 = lbound(tmpu, 2)
   j2 = ubound(tmpu, 2)
   i1 = lbound(tmpu, 1)
   i2 = ubound(tmpu, 1)

   allocate(kan1, mold=tmpu)
   allocate(kan2, mold=tmpu)
   allocate(kan3, mold=tmpu)
   allocate(kan, mold=tmpu)

   kan1 = 0.0
   kan2 = 0.0
   kan3 = 0.0
   kan  = UNDEF

   DUST_HETEROGENOUS_CHEM: if (associated(DU)) then
      DUST_REACTION_RATES: do n = 1, nbinsDU
         kan = 0.0
         call HNO3_reaction_rate(i1, i2, j1, j2, km, klid, &
                                 rmedDU(n), fnumDU(n), &
                                 rhoa, tmpu, relhum, DU(:,:,:,n), kan, &
                                 AVOGAD, AIRMW, PI, RUNIV, fMassHNO3)

         select case(n)
            case (1)
               kan1 = kan1 + kan
            case (2)
               kan2 = kan2 + kan
            case (3)
               kan2 = kan2 + kan
            case (4)
               kan3 = kan3 + kan
            case (5)
               kan3 = kan3 + kan
         end select

      end do DUST_REACTION_RATES
   end if DUST_HETEROGENOUS_CHEM


   SALT_HETEROGENOUS_CHEM: if (associated(SS)) then
      SALT_REACTION_RATES: do n = 1, nbinsSS
         kan = 0.0
         call SSLT_reaction_rate(i1, i2, j1, j2, km, klid, &
                                 rmedSS(n), fnumSS(n), &
                                 rhoa, tmpu, relhum, SS(:,:,:,n), kan, &
                                 AVOGAD, AIRMW, PI, RUNIV, fMassHNO3)

         select case(n)
            case (1)
               kan1 = kan1 + kan
            case (2)
               kan1 = kan1 + kan
            case (3)
               kan2 = kan2 + kan
            case (4)
               kan2 = kan2+ kan
            case (5)
               kan3 = kan3 + kan
         end select

      end do SALT_REACTION_RATES
   end if SALT_HETEROGENOUS_CHEM

!  Compute the nitric acid loss (but don't actually update)
   kan = max(0.0, (kan1 + kan2 + kan3))

   call apportion_reaction_rate(i1, i2, j1, j2, km, kan1, kan)
   call apportion_reaction_rate(i1, i2, j1, j2, km, kan2, kan)
   call apportion_reaction_rate(i1, i2, j1, j2, km, kan3, kan)

   allocate(deltahno3, mold=kan)
   deltahno3 = xhno3 * fMassHNO3 / AIRMW * (1.0 - exp(-kan*cdt))
   deltahno3 = max(0.0, deltahno3)

   xhno3 = xhno3 - deltahno3 * AIRMW / fMassHNO3

   nNO3an1 = nNO3an1 + kan1 * deltahno3 * fMassNO3 / fMassHNO3
   nNO3an2 = nNO3an2 + kan2 * deltahno3 * fMassNO3 / fMassHNO3
   nNO3an3 = nNO3an3 + kan3 * deltahno3 * fMassNO3 / fMassHNO3

   if(associated(NI_phet)) then
      NI_phet(:,:,1) = (1.0 / (grav*cdt)) * sum(kan1*deltahno3*delp, dim=3)
      NI_phet(:,:,2) = (1.0 / (grav*cdt)) * sum(kan2*deltahno3*delp, dim=3)
      NI_phet(:,:,3) = (1.0 / (grav*cdt)) * sum(kan3*deltahno3*delp, dim=3)
   end if

!  Output diagnostic HNO3
!  ----------------------
!  Calculate the HNO3 mass concentration
   if( associated(HNO3_conc) ) then
      HNO3_conc = xhno3 * fMassHNO3 / AIRMW * rhoa
   endif
!  Calculate the HNO3 surface mass concentration
   if( associated(HNO3_sfcmass) ) then
      HNO3_sfcmass(i1:i2,j1:j2) = xhno3(i1:i2,j1:j2,km) * fMassHNO3 / AIRMW * rhoa(i1:i2,j1:j2,km)
   endif
!  Calculate the HNO3 column loading
   if( associated(HNO3_colmass) ) then
      HNO3_colmass(i1:i2,j1:j2) = 0.
      do k = klid, km
        HNO3_colmass(i1:i2,j1:j2) &
         =   HNO3_colmass(i1:i2,j1:j2) + xhno3(i1:i2,j1:j2,k)*delp(i1:i2,j1:j2,k)/grav
      end do
   endif

   __RETURN__(__SUCCESS__)
   end subroutine NIheterogenousChem
