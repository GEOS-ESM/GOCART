  subroutine SulfateDistributeEmissions ( km, nbins, cdt, grav, nymd, nhms, &
                                          fMassSO4, fMassSO2, fSO4ant, eAircraftFuel, &
                                          nSO2, nSO4, &
                                          so2anthro_l1_src, so2anthro_l2_src, &
                                          so2biomass_src, dmso_conc, &
                                          so2ship_src, so4ship_src, &
                                          aircraft_fuel_src, &
                                          so2, so4, &
                                          oro, u10m, v10m, hghte, pblh, &
                                          tmpu, rhoa, delp, nVolc, &
                                          SU_emis, SU_SO4eman, SU_SO2eman, SU_SO2embb, &
!                                          maskString, gridMask, &
                                          aviation_layers,   &
                                          aviation_lto_src, &
                                          aviation_cds_src, &
                                          aviation_crs_src, rc)

! !USES:
   implicit NONE

! !INPUT PARAMETERS:
   integer, intent(in) :: km, nbins ! number model layers, and number of species respectively
   real, intent(in)    :: cdt, grav ! model time, and gravity respectively
   integer, intent(in) :: nymd, nhms
   real, intent(in)    :: fMassSO4  ! gram molecular weight of SO4
   real, intent(in)    :: fMassSO2  ! gram molecular weight of SO2
   real, intent(in)    :: fSO4ant ! Fraction of anthropogenic emissions that are SO4
   integer, intent(in) :: nSO2  ! index of SO2 relative to other sulfate tracers
   integer, intent(in) :: nSO4  ! index of SO2 relative to other sulfate tracers
   real, intent(in)    :: eAircraftFuel ! Aircraft emission factor: go from kg fuel to kg SO2
   real, dimension(:,:), intent(in) :: so2anthro_l1_src ! anthropogenic source surface[1]
   real, dimension(:,:), intent(in) :: so2anthro_l2_src ! anthropogenic source [1]
   real, dimension(:,:), intent(in) :: so2biomass_src ! biomass burning source [1]
   real, dimension(:,:), intent(in) :: dmso_conc ! DMS source [1]
   real, dimension(:,:), intent(in) :: so2ship_src ! SO2 ship emissions [1]
   real, dimension(:,:), intent(in) :: so4ship_src ! SO4 ship emissions [1]
   real, dimension(:,:,:), intent(in) :: aircraft_fuel_src ! aircraft fuel source [1]

   real, pointer, dimension(:,:), intent(in)    :: oro   ! orography flag
   real, pointer, dimension(:,:), intent(in)    :: u10m  ! 10-m u-wind component [m s-1]
   real, pointer, dimension(:,:), intent(in)    :: v10m  ! 10-m v-wind component [m s-1]
   real, pointer, dimension(:,:,:), intent(in)  :: hghte ! top of layer geopotential height [m]
   real, pointer, dimension(:,:), intent(in)    :: pblh
   real, pointer, dimension(:,:,:), intent(in)  :: tmpu  ! temperature [K]
   real, pointer, dimension(:,:,:), intent(in)  :: rhoa  ! Layer air density [kg/m^3]
   real, pointer, dimension(:,:,:), intent(in)  :: delp  ! pressure thickness [Pa]
   integer, intent(in) :: nVolc     ! number of volcanic emissions
   real, dimension(:), intent(in)  :: aviation_layers ! Heights [m] of LTO, CDS and CRS aviation emissions layers
   real, dimension(:,:), intent(in) :: aviation_cds_src ! Climb/Descent aircraft fuel emission [1]
   real, dimension(:,:), intent(in) :: aviation_crs_src ! Cruise aircraft fuel emission [1]
   real, dimension(:,:), intent(in) :: aviation_lto_src ! Landing/Take-off aircraft fuel emission [1]

! !INOUTPUT PARAMETERS:
   real, dimension(:,:,:), intent(inout)  :: so2, so4 ! Sulfate  internal state varaibles [kg/kg]
   real, pointer, dimension(:,:,:)  :: SU_emis      ! SU emissions, kg/m2/s
   real, pointer, dimension(:,:)  :: SU_SO4eman  ! SO4 anthro emissions, kg/m2/s
   real, pointer, dimension(:,:)  :: SU_SO2eman  ! SO2 anthro emissions, kg/m2/s
   real, pointer, dimension(:,:)  :: SU_SO2embb  ! SO2 bioburn emissions, kg/m2/s

!  OUTPUT PARAMETERS:
   integer, optional, intent(out)   :: rc    ! Error return code:
                                             !  0 - all is well

! !DESCRIPTION: SulfateDistributeEmissions - Adds sulfate source emission for one timestep
!               We have emissions from 4 sources, which are distributed
!               differently in the vertical
!               1) biomass burning - uniformly mixed in PBL (SO2)
!               2) anthropogenic l1 - emitted into lowest 100 m (SO2,SO4)
!               3) anthropogenic l2 - emitted into 100 - 500 m levels (SO2,SO4)
!               4) volcanic emissions
!               Additionally have a source of DMS from transfer from seawater
!               into lowest model layer
!               Consider factors in conversion: we estimate that 5% of sulfur
!               from anthropogenic sources (by mass) goes directly to SO4.
!
! !REVISION HISTORY:
!
!  06Nov2003, Colarco -  Based on Ginoux
!  17July2020, Sherman - Refactored for GOCART2G. Only uses intrinsic Fortran

! !Local Variables
   integer  ::  i, j, k
   integer  :: i1=1, j1=1, i2, j2

   real, dimension(:,:), allocatable :: srcSO2, srcSO4, srcDMS, srcSO4anthro, &
                                        srcSO2anthro, srcSO2bioburn
   real, allocatable, dimension(:,:)    :: hsurf

   real :: p1, z1, dz, deltaz, deltap, f100, f500, fPblh
   real :: zpbl
                          ! pressure at 100m, 500m, & PBLH
   real, dimension(:,:), allocatable :: p100, p500, pPblh, p0, z0, ps

   real, dimension(:,:,:), allocatable :: emis_aviation
   real, dimension(:,:,:), allocatable :: srcAviation
   real  :: z_lto_bot, z_lto_top
   real  :: z_cds_bot, z_cds_top
   real  :: z_crs_bot, z_crs_top

!EOP
!-------------------------------------------------------------------------
!  Begin

   i2 = size(rhoa,1)
   j2 = size(rhoa,2)
   allocate(hsurf(i1:i2,j1:j2))
   hsurf = hghte(i1:i2,j1:j2,km)

   allocate(srcSO2(i2,j2), srcSO4(i2,j2), srcDMS(i2,j2), srcSO4anthro(i2,j2), &
            srcSO2anthro(i2,j2), srcSO2bioburn(i2,j2))

!  Initialize local variables
!  --------------------------
   srcSO2 = 0.0
   srcSO4 = 0.0
   srcDMS = 0.0

   if ((nVolc <= 0) .and. associated(SU_emis)) SU_emis = 0.0 !SU_emis is usually set to zero in SUvolcanicEmissions.
!                                               !If there are no volcanic emissions, we need to set it to zero here.
   if (associated(SU_SO4eman)) SU_SO4eman = 0.0
   if (associated(SU_SO2eman)) SU_SO2eman = 0.0
   if (associated(SU_SO2embb)) SU_SO2embb = 0.0

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

   call distribute_aviation_emissions(delp, rhoa, z_lto_bot, z_lto_top, aviation_lto_src, &
                                      emis_aviation, i1, i2, j1, j2, km, grav)
   srcAviation = srcAviation + emis_aviation

   call distribute_aviation_emissions(delp, rhoa, z_cds_bot, z_cds_top, aviation_cds_src, &
                                      emis_aviation, i1, i2, j1, j2, km, grav)
   srcAviation = srcAviation + emis_aviation

   call distribute_aviation_emissions(delp, rhoa, z_crs_bot, z_crs_top, aviation_crs_src, &
                                      emis_aviation, i1, i2, j1, j2, km, grav)
   srcAviation = srcAviation + emis_aviation

!  Find the pressure of the 100m, 500m, and PBLH altitudes
   allocate(ps, mold=pblh)
   allocate(p0, mold=pblh)
   allocate(z0, mold=pblh)
   allocate(p100, mold=pblh)
   allocate(p500, mold=pblh)
   allocate(pPblh, mold=pblh)

   ps = 0.0
   do k = 1, km
    ps(i1:i2,j1:j2) = ps(i1:i2,j1:j2) + delp(i1:i2,j1:j2,k)
   end do
   p0 = ps
   z0(i1:i2,j1:j2) = 0.
   do k = km, 1, -1
    do j = j1, j2
     do i = i1, i2
      p1 = p0(i,j) - delp(i,j,k)
      dz = delp(i,j,k)/rhoa(i,j,k)/grav
      z1 = z0(i,j)+dz
      if(z0(i,j) .lt. 100 .and. z1 .ge. 100.) then
       deltaz = z1-100.
       deltap = deltaz*rhoa(i,j,k)*grav
       p100(i,j) = p1+deltap
      endif
      if(z0(i,j) .lt. 500 .and. z1 .ge. 500.) then
       deltaz = z1-500.
       deltap = deltaz*rhoa(i,j,k)*grav
       p500(i,j) = p1+deltap
      endif
      zpbl = max ( pblh(i,j), 100. )
      if(z0(i,j) .lt. zpbl .and. z1 .ge. zpbl) then
       deltaz = z1-zpbl
       deltap = deltaz*rhoa(i,j,k)*grav
       pPblh(i,j) = p1+deltap
      endif
      p0(i,j) = p1
      z0(i,j) = z1
     end do
    end do
   end do

!  Now update the tracer mixing ratios with the aerosol sources
   p0 = ps
   z0 = hsurf
   do k = km, 1, -1

    do j = j1, j2
     do i = i1, i2

      p1 = p0(i,j) - delp(i,j,k)
      f100 = 0.
      if(p1 .ge. p100(i,j)) f100 = delp(i,j,k)/(ps(i,j)-p100(i,j))
      if(p1 .lt. p100(i,j) .and. p0(i,j) .ge. p100(i,j)) &
       f100 = (p0(i,j)-p100(i,j))/(ps(i,j)-p100(i,j))

      f500 = 0.
      if(p0(i,j) .ge. p100(i,j) .and. p1 .lt. p100(i,j) .and. p1 .ge. p500(i,j)) &
       f500 = (p100(i,j)-p1)/(p100(i,j)-p500(i,j))
      if(p0(i,j) .lt. p100(i,j) .and. p1 .ge. p500(i,j)) &
       f500 = delp(i,j,k)/(p100(i,j)-p500(i,j))
      if(p0(i,j) .ge. p500(i,j) .and. p1 .lt. p500(i,j)) &
       f500 = (p0(i,j)-p500(i,j))/(p100(i,j)-p500(i,j))

      fPblh = 0.
      if(p1 .ge. pPblh(i,j)) fPblh = delp(i,j,k)/(ps(i,j)-pPblh(i,j))
      if(p1 .lt. pPblh(i,j) .and. p0(i,j) .ge. pPblh(i,j)) &
       fPblh = (p0(i,j)-pPblh(i,j))/(ps(i,j)-pPblh(i,j))

!     All source from files specified in kg SO2 m-2 s-1 (unless filename
!     indicates otherwise!).
      srcSO4anthro(i,j) = fSO4ant * fMassSO4/fMassSO2 * &
                (   f100 * so2anthro_l1_src(i,j) &
                  + f500 * so2anthro_l2_src(i,j)  )
      srcSO2anthro(i,j) = (1.-fSO4ant) * &
                (   f100 * so2anthro_l1_src(i,j) &
                  + f500 * so2anthro_l2_src(i,j)  )

      srcSO2bioburn(i,j) = fPblh*so2biomass_src(i,j)

!     Add the ship emissions to anthro
      srcSO2anthro(i,j) = srcSO2anthro(i,j) + f100*so2ship_src(i,j)
      srcSO4anthro(i,j) = srcSO4anthro(i,j) + f100*so4ship_src(i,j)

!     Add the aircraft fuel emissions to anthro SO2
      srcSO2anthro(i,j) = srcSO2anthro(i,j) + &
       eAircraftFuel * aircraft_fuel_src(i,j,k)

      srcSO2anthro(i,j) = srcSO2anthro(i,j) + srcAviation(i,j,k)

      srcSO4(i,j) = srcSO4anthro(i,j)
      srcSO2(i,j) = srcSO2anthro(i,j)+srcSO2bioburn(i,j)

      so2(i,j,k)  =   so2(i,j,k) + srcSO2(i,j)*cdt*grav/delp(i,j,k)
      so4(i,j,k)  =   so4(i,j,k) + srcSO4(i,j)*cdt*grav/delp(i,j,k)

      p0(i,j) = p1

     end do ! i
    end do  ! j

    if (associated(SU_emis)) SU_emis(:,:,nSO2) = SU_emis(:,:,nSO2) + srcSO2
    if (associated(SU_emis)) SU_emis(:,:,nSO4) = SU_emis(:,:,nSO4) + srcSO4
    if (associated(SU_SO4eman)) SU_SO4eman = SU_SO4eman + srcSO4anthro
    if (associated(SU_SO2eman)) SU_SO2eman = SU_SO2eman + srcSO2anthro
    if (associated(SU_SO2embb)) SU_SO2embb = SU_SO2embb + srcSO2bioburn

   end do ! k

   __RETURN__(__SUCCESS__)
  end subroutine SulfateDistributeEmissions
