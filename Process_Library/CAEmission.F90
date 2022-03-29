   subroutine CAEmission (mie, km, nbins, cdt, grav, prefix, ratPOM, eAircraftfuel, aircraft_fuel_src, &
                           aviation_lto_src, aviation_cds_src, aviation_crs_src, &
                           fHydrophobic, pblh, tmpu, rhoa, rh, aerosolPhilic, aerosolPhobic, &
                           delp, aviation_layers, &
                            biomass_src, terpene_src, eocant1_src, eocant2_src, oc_ship_src, biofuel_src, &
                           OC_emis, OC_emisAN, OC_emisBB, OC_emisBF, OC_emisBG, rc )

! !USES:

  implicit NONE

! !INPUT PARAMETERS:
   type(GOCART2G_Mie),  intent(in) :: mie        ! mie table
   integer, intent(in) :: km     ! total model levels
   integer, intent(in) :: nbins  ! number of aerosol size bins
   real, intent(in)    :: cdt    ! chemistry model time-step [sec]
   real, intent(in)    :: grav   ! gravity [m/sec^2]
   character(len=2), intent(in)  :: prefix ! varaible name prefix
   real, intent(in)    :: ratPOM
   real, intent(in)    :: eAircraftFuel ! Aircraft emission factor: go from kg fuel to kg C
   real, dimension(:), intent(in)  :: aviation_layers ! Heights [m] of LTO, CDS and CRS aviation emissions layers
   real, pointer, dimension(:,:), intent(in)    :: pblh  ! PBL height [m]
   real, pointer, dimension(:,:,:), intent(in)  :: tmpu  ! temperature [K]
   real, pointer, dimension(:,:,:), intent(in)  :: rhoa  ! air density [kg m-3]
   real, pointer, dimension(:,:,:), intent(in)  :: rh    ! relative humidity [1]
   real, pointer, dimension(:,:,:), intent(in)  :: delp  ! pressure level thickness [Pa]
   real, dimension(:,:,:), intent(in) :: aircraft_fuel_src ! aircraft fuel source [1]
   real, dimension(:,:), intent(in) :: aviation_cds_src ! Climb/Descent aircraft fuel emission [1]
   real, dimension(:,:), intent(in) :: aviation_crs_src ! Cruise aircraft fuel emission [1]
   real, dimension(:,:), intent(in) :: aviation_lto_src ! Landing/Take-off aircraft fuel emission [1]
   real, dimension(:,:), intent(in) :: biomass_src
   real, dimension(:,:), intent(in) :: terpene_src
   real, dimension(:,:), intent(in) :: eocant1_src  ! anthropogenic emissions
   real, dimension(:,:), intent(in) :: eocant2_src  ! anthropogenic emissions
   real, dimension(:,:), intent(in) :: oc_ship_src  ! ship emissions
   real, dimension(:,:), intent(in) :: biofuel_src  ! biofuel emissions
   real, intent(in) :: fHydrophobic

! !OUTPUT PARAMETERS:
   real, dimension(:,:,:), intent(inout) :: aerosolPhobic
   real, dimension(:,:,:), intent(inout) :: aerosolPhilic
   real, pointer, dimension(:,:,:)  :: OC_emis  ! OC emissions, kg/m2/s
   real, pointer, dimension(:,:)  :: OC_emisAN  ! OC emissions, kg/m2/s
   real, pointer, dimension(:,:)  :: OC_emisBB  ! OC emissions, kg/m2/s
   real, pointer, dimension(:,:)  :: OC_emisBF  ! OC emissions, kg/m2/s
   real, pointer, dimension(:,:)  :: OC_emisBG  ! OC emissions, kg/m2/s
   integer, optional, intent(out) :: rc         ! Error return code:
                                                !  0 - all is well
                                                !  1 -
   character(len=*), parameter :: myname = 'CAEmission'

! !DESCRIPTION: Updates the OC concentration with emissions every timestep
!
! !REVISION HISTORY:
!
!  06Nov2003, Colarco
!  Based on Ginoux
!    June2020 E.Sherman - moved to process library
!EOP
!-------------------------------------------------------------------------

! !Local Variables
   integer  ::  i, j, k, n, ios, ijl
   integer  ::  n1, n2
   integer  :: i1=1, i2, j1=1, j2
!  pressure at 100m, 500m, & PBLH
   real, dimension(:,:), allocatable :: p100, p500, pPBL
   real, dimension(:,:), allocatable :: p0, z0, ps
   real :: p1, z1, dz, delz, delp_, f100, f500, fPBL, fBot
   real :: qmax, qmin, eBiofuel, eBiomass, eTerpene, eAnthro

   real, dimension(:,:), allocatable :: factor, srcHydrophobic, srcHydrophilic
   real, dimension(:,:), allocatable :: srcBiofuel, srcBiomass, srcAnthro, srcBiogenic
   real                         :: srcTmp, zpbl, maxAll

   real, dimension(:,:,:), allocatable :: emis_aviation
   real, dimension(:,:,:), allocatable :: srcAviation
   real                            :: z_lto_bot, z_lto_top
   real                            :: z_cds_bot, z_cds_top
   real                            :: z_crs_bot, z_crs_top

   real, dimension(:,:), allocatable          :: f_bb_        ! scaling factor for BB emissions based on maximum allowed exttau
   real, dimension(:,:), allocatable          :: exttau_bb_   ! increment of exttau due to BB during the current time step
   real, allocatable, dimension(:,:,:,:) :: qa_bb_       ! increment of qa due to BB during the current time step (nbins,i1:i2,j1:j2:km) 
                                                         ! W.Jiang note, changed to (i1:i2,j1:j2,km,nbins) for efficiency
   real                                  :: cutoff_bb_exttau
   integer                               :: idx
   integer                               :: ilam550, status
   real                                  :: wavelength550
   real, dimension(:,:,:), allocatable   :: tau
   character(len=255)                    :: qname
   real, parameter                       :: max_bb_exttau = 30.0

!  Indices for point emissions
   real, dimension(km)          :: point_column_emissions

!  Source function terms for SOA from Anthropogenic VOCs
   real :: srcSOAanthro = 0.0
!  Initialize local variables
!  --------------------------
   i2 = size(rhoa,1)
   j2 = size(rhoa,2)
   ijl = ( i2 - i1 + 1 ) * ( j2 - j1 + 1 )

   allocate(factor(i2,j2), srcHydrophobic(i2,j2), srcHydrophilic(i2,j2), srcBiofuel(i2,j2), &
            srcBiomass(i2,j2), srcAnthro(i2,j2), srcBiogenic(i2,j2), f_bb_(i2,j2), exttau_bb_(i2,j2))

!  Emission factors scaling from source files to desired mass quantity
   eBiomass = ratPOM
   eBiofuel = ratPOM
   eTerpene = ratPOM
   eAnthro  = ratPOM

!  Zero diagnostic accumulators
     if(associated(OC_emis)) OC_emis = 0.0
     if(associated(OC_emisAN)) OC_emisAN = 0.0
     if(associated(OC_emisBF)) OC_emisBF = 0.0
     if(associated(OC_emisBB)) OC_emisBB = 0.0
     if(associated(OC_emisBG)) OC_emisBG = 0.0

!  Distribute aircraft emissions from LTO, CDS and CRS layers
!  ----------------------------------------------------------
   z_lto_bot = max(1e-3, aviation_layers(1))
   z_lto_top = max(2e-3, aviation_layers(2))

   z_cds_bot = max(2e-3, aviation_layers(2))
   z_cds_top = max(3e-3, aviation_layers(3))

   z_crs_bot = max(3e-3, aviation_layers(3))
   z_crs_top = max(4e-3, aviation_layers(4))

   allocate(emis_aviation, mold=tmpu)
   allocate(srcAviation, mold=tmpu)
   emis_aviation = 0.0
   srcAviation   = 0.0

   call distribute_aviation_emissions(delp, rhoa, z_lto_bot, z_lto_top, aviation_lto_src, emis_aviation, i1, i2, j1, j2, km, grav)
   srcAviation = srcAviation + emis_aviation

   call distribute_aviation_emissions(delp, rhoa, z_cds_bot, z_cds_top, aviation_cds_src, emis_aviation, i1, i2, j1, j2, km, grav)
   srcAviation = srcAviation + emis_aviation

   call distribute_aviation_emissions(delp, rhoa, z_crs_bot, z_crs_top, aviation_crs_src, emis_aviation, i1, i2, j1, j2, km, grav)
   srcAviation = srcAviation + emis_aviation


!  Determine surface pressure
!  AMS Note: pass this in
!  --------------------------
   allocate(ps, mold=pblh)
   allocate(p0, mold=pblh)
   allocate(z0, mold=pblh)
   allocate(p100, mold=pblh)
   allocate(p500, mold=pblh)
   allocate(pPBL, mold=pblh)
   ps = 0.0
   do k = 1, km
    ps(i1:i2,j1:j2) = ps(i1:i2,j1:j2) + delp(i1:i2,j1:j2,k)
   end do

!  Find the pressure of the 100m, 500m, and PBLH altitudes
!  AMS Note: this could be greatly simplified by using ze/zm and having a
!      generic routine from the bottom up with an early exit condition
!  -----------------------------------------------------------------------
   p0 = ps
   z0(i1:i2,j1:j2) = 0.
   do k = km, 1, -1
    do j = j1, j2
     do i = i1, i2
      p1 = p0(i,j) - delp(i,j,k)
      dz = delp(i,j,k)/rhoa(i,j,k)/grav
      z1 = z0(i,j)+dz
      if(z0(i,j) .lt. 100 .and. z1 .ge. 100.) then
       delz = z1-100.
       delp_ = delz*rhoa(i,j,k)*grav
       p100(i,j) = p1+delp_
      endif
      if(z0(i,j) .lt. 500 .and. z1 .ge. 500.) then
       delz = z1-500.
       delp_ = delz*rhoa(i,j,k)*grav
       p500(i,j) = p1+delp_
      endif
      zpbl = max ( pblh(i,j), 100. )
      if(z0(i,j) .lt. zpbl .and. z1 .ge. zpbl ) then
       delz = z1-zpbl
       delp_ = delz*rhoa(i,j,k)*grav
       pPBL(i,j) = p1+delp_
      endif
      p0(i,j) = p1
      z0(i,j) = z1
     end do
    end do
   end do

!   Limit biomass burning emissions
!   -------------------------------
    allocate(qa_bb_(i1:i2,j1:j2,km,nbins))
    qa_bb_ = 0.0

    p0 = ps
K_LOOP_BB: do k = km, 1, -1

!   First determine emissions for this layer
!   ----------------------------------------
    maxAll = 0.0
    do j = j1, j2
     do i = i1, i2

      p1 = p0(i,j) - delp(i,j,k)

!     Pressure @ PBL height
!     ---------------------
      fPBL = 0.
      if(p1 .ge. pPBL(i,j)) fPBL = delp(i,j,k)/(ps(i,j)-pPBL(i,j))
      if(p1 .lt. pPBL(i,j) .and. p0(i,j) .ge. pPBL(i,j)) fPBL = (p0(i,j)-pPBL(i,j))/(ps(i,j)-pPBL(i,j))

!     Sources by class in kg m-2 s-1
!     ------------------------------
      srcBiomass(i,j)  = fPBL * eBiomass * biomass_src(i,j)

      srcHydrophobic(i,j) =     fHydrophobic  * srcBiomass(i,j)
      srcHydrophilic(i,j) = (1.-fHydrophobic) * srcBiomass(i,j)

!     Update pressure of lower level
!     ------------------------------
      p0(i,j) = p1

     end do ! i
    end do  ! j


!   Update concentrations at this layer
!   The "1" element is hydrophobic
!   The "2" element is hydrophilic
!   -----------------------------------
    factor = cdt * grav / delp(:,:,k)

    qa_bb_(:,:,k,1) = factor * srcHydrophobic
    qa_bb_(:,:,k,2) = factor * srcHydrophilic

   end do K_LOOP_BB

!   Get the wavelength indices
!   --------------------------
!   Must provide ilam550 for AOT calculation
    ilam550 = mie%getChannel(5.50e-7)
    if (ilam550 <=0) ilam550 = 1
    wavelength550 = mie%getWavelength(ilam550, __RC__)

!  Calculate the extinction and/or scattering AOD

   exttau_bb_(i1:i2,j1:j2) = 0.0
   allocate(tau(i1:i2,j1:j2,km), source = 0.)
   do n = 1, nbins
!     Select the name for species and the index
     call mie%Query(wavelength550, n,           &
              qa_bb_(:,:,:,n)*delp(:,:,:)/grav, &
              rh, tau=tau, __RC__)
     do k = 1, km
!        Integrate in the vertical
        exttau_bb_(:,:) = exttau_bb_(:,:) + tau(:,:,k)
     enddo
   enddo  ! nbins

   f_bb_ = 1.0
   cutoff_bb_exttau = (cdt / (24 * 3600.0)) * max_bb_exttau

   do j = j1, j2
    do i = i1, i2
     if (exttau_bb_(i,j) > cutoff_bb_exttau) then
      f_bb_(i,j) = cutoff_bb_exttau / exttau_bb_(i,j)
     end if
    enddo
   enddo

   deallocate(qa_bb_, tau)

!  Now update the tracer mixing ratios with the aerosol sources
!  ------------------------------------------------------------
   p0 = ps
K_LOOP: do k = km, 1, -1

!   First determine emissions for this layer
!   ----------------------------------------
    maxAll = 0.0
    do j = j1, j2
     do i = i1, i2

      p1 = p0(i,j) - delp(i,j,k)

!     Pressure @ 100m
!     ---------------
      f100 = 0.
      if(p1 .ge. p100(i,j)) f100 = delp(i,j,k)/(ps(i,j)-p100(i,j))
      if(p1 .lt. p100(i,j) .and. p0(i,j) .ge. p100(i,j)) &
       f100 = (p0(i,j)-p100(i,j))/(ps(i,j)-p100(i,j))

!     Pressure @ 500m
!     ---------------
      f500 = 0.
      if ( p0(i,j) .ge. p100(i,j) .and. p1 .lt. p100(i,j) .and. p1 .ge. p500(i,j)) &
       f500 = (p100(i,j)-p1)/(p100(i,j)-p500(i,j))
      if(p0(i,j) .lt. p100(i,j) .and. p1 .ge. p500(i,j)) &
       f500 = delp(i,j,k)/(p100(i,j)-p500(i,j))
      if(p0(i,j) .ge. p500(i,j) .and. p1 .lt. p500(i,j)) &
       f500 = (p0(i,j)-p500(i,j))/(p100(i,j)-p500(i,j))

!     Pressure @ PBL height
!     ---------------------
      fPBL = 0.
      if(p1 .ge. pPBL(i,j)) fPBL = delp(i,j,k)/(ps(i,j)-pPBL(i,j))
      if(p1 .lt. pPBL(i,j) .and. p0(i,j) .ge. pPBL(i,j)) &
       fPBL = (p0(i,j)-pPBL(i,j))/(ps(i,j)-pPBL(i,j))

!     Terpene is tree-top emission; only add in bottom layer
!     ------------------------------------------------------
      if ( k .eq. km ) then
         fBot = 1.0
      else
         fBot = 0.0
      end if

!     Sources by class in kg m-2 s-1
!     ------------------------------
      srcBiofuel(i,j)  = f100 * eBiofuel * biofuel_src(i,j)
      srcAnthro(i,j)   = f100 * eAnthro  * eocant1_src(i,j) &
                       + f500 * eAnthro  * eocant2_src(i,j) &
                       + f100 * eAnthro  * oc_ship_src(i,j) &
                       +        eAnthro  * srcAviation(i,j,k) &
                       +        eAnthro  * eAircraftFuel * aircraft_fuel_src(i,j,k)
      if ((prefix == 'OC') .or. (prefix == 'BR')) then
         srcBiomass(i,j)  = fPBL * eBiomass * biomass_src(i,j) * f_bb_(i,j)
      else
         srcBiomass(i,j)  = fPBL * eBiomass * biomass_src(i,j)
      end if

      srcBiogenic(i,j) = fBot * eTerpene * terpene_src(i,j) !Black carbon has no biogenic source. Should be zeros.

      srcTmp = srcBiofuel(i,j) + srcAnthro(i,j) + srcBiomass(i,j)

      srcHydrophobic(i,j) =     fHydrophobic  * srcTmp
      srcHydrophilic(i,j) = (1.-fHydrophobic) * srcTmp + srcBiogenic(i,j)

!     Update pressure of lower level
!     ------------------------------
      p0(i,j) = p1

     end do ! i
    end do  ! j

!   Update concentrations at this layer
!   The "1" element is hydrophobic
!   The "2" element is hydrophilic
!   -----------------------------------
    factor = cdt * grav / delp(:,:,k)

    aerosolPhobic(:,:,k) = aerosolPhobic(:,:,k) &
                             + factor * srcHydrophobic

    aerosolPhilic(:,:,k) = aerosolPhilic(:,:,k) &
                             + factor * srcHydrophilic

!   Fill in diagnostics if requested
!   --------------------------------
    if ( associated(OC_emis)) &
                    OC_emis(:,:,1) = OC_emis(:,:,1) + srcHydrophobic

    if ( associated(OC_emis)) &
                    OC_emis(:,:,2) = OC_emis(:,:,2) + srcHydrophilic

    if ( associated(OC_emisBF)) &
                    OC_emisBF  = OC_emisBF  + srcBiofuel

    if ( associated(OC_emisBB)) &
                    OC_emisBB  = OC_emisBB  + srcBiomass

    if ( associated(OC_emisAN)) &
                    OC_emisAN  = OC_emisAN  + srcAnthro

    if ( associated(OC_emisBG)) &
                    OC_emisBG  = OC_emisBG + srcBiogenic
   end do K_LOOP

   __RETURN__(__SUCCESS__)
   end subroutine CAEmission
