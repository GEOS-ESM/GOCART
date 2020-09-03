#include "MAPL_Generic.h"
#include "unused_dummy.H"
!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  SulfateChemDriverMod.F90 --- Calculate the sulfate aerosol
!                                        chemistry
!
! !INTERFACE:
!

   module  SulfateChemDriverMod

! !USES:

   USE ESMF
   USE MAPL
   use m_StrTemplate
   use m_die, only: die

   use Chem_Mod
   use Chem_ConstMod, only: grav, undefval => undef, &
                            airMolWght => airmw        ! Constants !
   use Chem_UtilMod
   use DryDepositionMod

   use m_mpout

!   use pflogger

   implicit none

! !PUBLIC TYPES:
!
   PRIVATE

!
! !PUBLIC MEMBER FUNCTIONS:
!

   PUBLIC  SulfateChemDriverGOCART
   PUBLIC  SulfateUpdateEmissions
   PUBLIC  SulfateUpdateOxidants
   PUBLIC  SulfateDistributeEmissions

  real, parameter :: OCEAN=0.0, LAND = 1.0, SEA_ICE = 2.0

  real, parameter :: pi = 3.1415, rearth = 6.37e6
  real, parameter :: radToDeg = 57.2957795
  real, parameter :: nAvogadro  = 6.022e23 ! molecules per mole of air

! gram molecular weights of species
  real, parameter :: fMassSulfur = 32., fMassSO2 = 64., fMassSO4 = 96., &
                     fMassDMS = 62., fMassMSA = 96.

! relative position of sulfate tracers
  integer, parameter :: nDMS = 1, &
                        nSO2 = 2, &
                        nSO4 = 3, &
                        nMSA = 4


!
! !DESCRIPTION:
!
!  This module implements the Dust Emission calculations
!
! !REVISION HISTORY:
!
!  29Dec2009 Colarco    First crack!
!
!EOP
!-------------------------------------------------------------------------
CONTAINS
!
! GOCART-based sulfate chemistry driver
!
! !ROUTINE:  SulfateChemDriverGOCART
!
! !INTERFACE:
!
   subroutine SulfateChemDriverGOCART ( i1, i2, j1, j2, km, nbeg, &
                                        nbins, cdt, nymd, nhms, &
                                        lonRad, latRad, &
                                        dms, so2, so4, msa, &
                                        xoh, xno3, xh2o2, &
                                        u, v, delp, tmpu, cloud, rhoa, hghte, &
                                        ustar, shflux, oro, pblh, z0h, &
                                        SU_dep, SU_PSO2, SU_PMSA, &
                                        SU_PSO4, SU_PSO4g, SU_PSO4aq, &     ! 2d diagnostics
                                        pso2, pmsa, pso4, pso4g, pso4aq,  & ! 3d diagnostics
                                        rc)


! !USES:

  implicit NONE

! !INPUT PARAMETERS:

   integer, intent(in)                 :: i1, i2, j1, j2, km, nbeg, nbins, &
                                          nymd, nhms
   real, intent(in)                    :: cdt
   real, pointer, dimension(:,:,:)     :: dms, so2, so4, msa
   real, pointer, dimension(:,:,:)     :: xoh, xno3, xh2o2
   real, pointer, dimension(:,:,:)     :: u, v, delp, tmpu, cloud, rhoa, hghte
   real, pointer, dimension(:,:)       :: ustar, shflux, oro, pblh, z0h
   real, pointer, dimension(:,:)       :: lonRad, latRad

! !OUTPUT PARAMETERS:
   type(Chem_Array), intent(inout)     :: SU_dep(nbins)  ! Mass lost by deposition
                                                         ! to surface, kg/m2/s
!  chemical production terms d(mixing ratio) /s
   type(Chem_Array), intent(inout)     :: su_pSO2, su_pMSA, su_PSO4, su_pSO4g, su_pSO4aq
   type(Chem_Array), intent(inout)     :: pSO2, pMSA, pSO4, pSO4g, pSO4aq 

   integer, intent(out)                :: rc          ! Error return code:
                                                      !  0 - all is well
                                                      !  1 - 
   character(len=*), parameter         :: myname = 'SU_ChemDrv'

! !DESCRIPTION: Updates the SU concentration due to chemistry
!  The SU grid component is currently established with 4 different
!  species (bins) following this convection:
!   1) DMS
!   2) SO2
!   3) SO4
!   4) MSA
!  Accordingly we have 4 chemical cycles to follow through, which are
!  sub-subroutines under this one.
!  The chemistry is a function of OH, NO3, and H2O2 concentrations
!  as well as DMS, SO2, SO4, MSA concentrations.  It is also a function
!  of solar zenith angle and temperature.  We pass in temperature.  SZA
!  will be a function of time of day and lat/lon.  For now we simply add
!  this to the grid component before calculating it.  I bet this is
!  somewhere else in the model.
!
! !REVISION HISTORY:
!
!  06Nov2003, Colarco
!  Based on Ginoux
!
!EOP
!-------------------------------------------------------------------------

! !Local Variables
   real    :: cossza(i1:i2,j1:j2), sza(i1:i2,j1:j2)
   integer :: k, n, jday
   real    :: pSO2_DMS(i1:i2,j1:j2,1:km), pMSA_DMS(i1:i2,j1:j2,1:km), &
              pSO4g_SO2(i1:i2,j1:j2,1:km), pSO4aq_SO2(i1:i2,j1:j2,1:km)
   real    :: drydepositionfrequency(i1:i2,j1:j2)
   real    :: xhour
   integer :: ijl, ijkl
#ifdef DEBUG
   real :: qmin, qmax
#endif
   _UNUSED_DUMMY(nbeg)
   _UNUSED_DUMMY(u)
   _UNUSED_DUMMY(v)

   ijl  = ( i2 - i1 + 1 ) * ( j2 - j1 + 1 )
   ijkl = ijl * km

!  Reset the production terms
   pSO2_DMS(i1:i2,j1:j2,1:km) = 0.
   pMSA_DMS(i1:i2,j1:j2,1:km) = 0.
   pSO4g_SO2(i1:i2,j1:j2,1:km) = 0.
   pSO4aq_SO2(i1:i2,j1:j2,1:km) = 0.
   if( associated(su_pSO2%data2d) )   su_pSO2%data2d(i1:i2,j1:j2) = 0.
   if( associated(su_pMSA%data2d) )   su_pMSA%data2d(i1:i2,j1:j2) = 0.
   if( associated(su_pSO4%data2d) )   su_pSO4%data2d(i1:i2,j1:j2) = 0.
   if( associated(su_pSO4g%data2d) )  su_pSO4g%data2d(i1:i2,j1:j2) = 0.
   if( associated(su_pSO4aq%data2d) ) su_pSO4aq%data2d(i1:i2,j1:j2) = 0.
   if( associated(pSO2%data3d) )      pSO2%data3d(i1:i2,j1:j2,1:km) = 0.
   if( associated(pMSA%data3d) )      pMSA%data3d(i1:i2,j1:j2,1:km) = 0.
   if( associated(pSO4%data3d) )      pSO4%data3d(i1:i2,j1:j2,1:km) = 0.
   if( associated(pSO4g%data3d) )     pSO4g%data3d(i1:i2,j1:j2,1:km) = 0.
   if( associated(pSO4aq%data3d) )    pSO4aq%data3d(i1:i2,j1:j2,1:km) = 0.


!  Find the cossza
!  ----------------------------------
   jday = idaynum(nymd)
   xhour = (  real(nhms/10000)*3600. &
            + real(mod(nhms,10000)/100)*60. &
            + real(mod(nhms,100)) &
           ) / 3600.
   call szangle(jday, xHour, lonRad, latRad, sza, cossza, i1, i2, j1, j2)

!  Reset the dry deposition fluxes & frequencies
   do n = 1, nbins
    if( associated(su_dep(n)%data2d) ) su_dep(n)%data2d(i1:i2,j1:j2) = 0.0
   end do
   call DryDepositionGOCART( i1, i2, j1, j2, km, &
                             tmpu, rhoa, hghte, oro, ustar, &
                             pblh, shflux, z0h, drydepositionfrequency, rc )

if(mapl_am_i_root()) print*,'SU sum(drydepositionfrequency) = ',sum(drydepositionfrequency)

!  Now call the chemistry packages...
!  ----------------------------------

!  DMS source and oxidation to SO2 and MSA
   call SU_ChemDrv_DMS( i1, i2, j1, j2, km, cdt, xoh, xno3, rhoa, &
                        dms, delp, tmpu, cossza, pSO2_DMS, pMSA_DMS, &
                        drydepositionfrequency, su_dep(nDMS), rc)

   if( associated(pSO2%data3d) ) &
     pSO2%data3d(i1:i2,j1:j2,1:km) = pSO2_DMS(i1:i2,j1:j2,1:km)
   if( associated(su_pSO2%data2d)) then
     do k = 1, km
      su_pSO2%data2d(i1:i2,j1:j2) &
        =   su_pSO2%data2d(i1:i2,j1:j2) &
          + pSO2_DMS(i1:i2,j1:j2,k)*delp(i1:i2,j1:j2,k)/grav
     enddo
   endif

   if( associated(pMSA%data3d) ) &
     pMSA%data3d(i1:i2,j1:j2,1:km) = pMSA_DMS(i1:i2,j1:j2,1:km)
   if( associated(su_pMSA%data2d)) then
     do k = 1, km
      su_pMSA%data2d(i1:i2,j1:j2) &
        =   su_pMSA%data2d(i1:i2,j1:j2) &
          + pMSA_DMS(i1:i2,j1:j2,k)*delp(i1:i2,j1:j2,k)/grav
     enddo
   endif

!  SO2 source and oxidation to SO4
   call SU_ChemDrv_SO2( i1, i2, j1, j2, km, cdt, xoh, xh2o2, rhoa, &
                        so2, delp, tmpu, cloud, pSO2_DMS, pSO4g_SO2, &
                        pSO4aq_SO2, drydepositionfrequency, oro, su_dep(nSO2), rc)
   if( associated(pSO4g%data3d) ) &
     pSO4g%data3d(i1:i2,j1:j2,1:km) = pSO4g_SO2(i1:i2,j1:j2,1:km)
   if( associated(su_pSO4g%data2d)) then
     do k = 1, km
      su_pSO4g%data2d(i1:i2,j1:j2) &
        =   su_pSO4g%data2d(i1:i2,j1:j2) &
          + pSO4g_SO2(i1:i2,j1:j2,k)*delp(i1:i2,j1:j2,k)/grav
     enddo
   endif

   if( associated(pSO4aq%data3d) ) &
     pSO4aq%data3d(i1:i2,j1:j2,1:km) = pSO4aq_SO2(i1:i2,j1:j2,1:km)
   if( associated(su_pSO4aq%data2d)) then
     do k = 1, km
      su_pSO4aq%data2d(i1:i2,j1:j2) &
        =   su_pSO4aq%data2d(i1:i2,j1:j2) &
          + pSO4aq_SO2(i1:i2,j1:j2,k)*delp(i1:i2,j1:j2,k)/grav
     enddo
   endif

   if( associated(pSO4%data3d) ) &
     pSO4%data3d(i1:i2,j1:j2,1:km) = pSO4g_SO2(i1:i2,j1:j2,1:km) + pSO4aq_SO2(i1:i2,j1:j2,1:km)
   if( associated(su_pSO4%data2d)) then
     do k = 1, km
      su_pSO4%data2d(i1:i2,j1:j2) &
        =   su_pSO4%data2d(i1:i2,j1:j2) &
          + pSO4g_SO2(i1:i2,j1:j2,k)*delp(i1:i2,j1:j2,k)/grav   &
          + pSO4aq_SO2(i1:i2,j1:j2,k)*delp(i1:i2,j1:j2,k)/grav
     enddo
   endif



!  SO4 source and loss
   call SU_ChemDrv_SO4( i1, i2, j1, j2, km, cdt, so4, delp, &
                        pSO4g_SO2, pSO4aq_SO2, drydepositionfrequency, su_dep(nSO4), rc)


!  MSA source and loss
   call SU_ChemDrv_MSA( i1, i2, j1, j2, km, cdt, msa, delp, &
                        pMSA_DMS, drydepositionfrequency, su_dep(nMSA), rc)

#ifdef DEBUG
   if(associated(su_pso2%data2d)) call pmaxmin('SU: su_pso2',su_pso2%data2d,qmin,qmax,ijl,1,1.)
   if(associated(su_pmsa%data2d)) call pmaxmin('SU: su_pmsa',su_pmsa%data2d,qmin,qmax,ijl,1,1.)
   if(associated(su_pso4g%data2d)) call pmaxmin('SU: su_pso4g',su_pso4g%data2d,qmin,qmax,ijl,1,1.)
   if(associated(su_pso4aq%data2d)) call pmaxmin('SU: su_pso4aq',su_pso4aq%data2d,qmin,qmax,ijl,1,1.)
   call pmaxmin('SU:  pSO4g_SO2',  pSO4g_SO2, qmin, qmax, ijl, km, 1. )
   call pmaxmin('SU: pSO4aq_SO2', pSO4aq_SO2, qmin, qmax, ijl, km, 1. )
   call pmaxmin('SU:   pSO2_DMS',   pSO2_DMS, qmin, qmax, ijl, km, 1. )
   call pmaxmin('SU:   pMSA_DMS',   pMSA_DMS, qmin, qmax, ijl, km, 1. )
#endif

   rc = 0

   end subroutine SulfateChemDriverGOCART


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  SU_ChemDrv_DMS
!  NOTE: This is the DMS oxidation subroutine.

!
! !INTERFACE:
!

   subroutine SU_ChemDrv_DMS( i1, i2, j1, j2, km, cdt, xoh, xno3, rhoa, &
                              qa, delp, tmpu, cossza, pSO2_DMS, pMSA_DMS, &
                              drydepf, fluxOut, rc) 

! !USES:

  implicit NONE

! !INPUT PARAMETERS:

   integer, intent(in)              :: i1, i2, j1, j2, km
   real, intent(in)                 :: cdt
   real, pointer                    :: qa(:,:,:)
   real, pointer                    :: delp(:,:,:)
   real, pointer, dimension(:,:,:)  :: tmpu, rhoa
   real, intent(in)                 :: drydepf(i1:i2,j1:j2)
   real, intent(in)                 :: xoh(i1:i2,j1:j2,km), xno3(i1:i2,j1:j2,km)
   real, intent(in)                 :: cossza(i1:i2,j1:j2)


! !OUTPUT PARAMETERS:

   type(Chem_Array), intent(inout)  :: fluxout     ! Mass lost by deposition
                                                   ! to surface, kg/m2/s
   real, intent(out)                :: pSO2_DMS(i1:i2,j1:j2,km), pMSA_DMS(i1:i2,j1:j2,km)
   integer, intent(out)             :: rc          ! Error return code:
                                                   !  0 - all is well
                                                   !  1 - 
   character(len=*), parameter :: myname = 'SU_ChemDrv_DMS'

! !DESCRIPTION: Computes the production of SO2 and MSA due to DMS oxidation
!
!   R1:    DMS + OH  -> a*SO2 + b*MSA                OH addition channel
!          k1 = { 1.7d-42*exp(7810/T)*[O2] / (1+5.5e-31*exp(7460/T)*[O2] }
!          a = 0.75, b = 0.25
!
!   R2:    DMS + OH  ->   SO2 + ...                  OH abstraction channel
!          k2 = 1.2e-11*exp(-260/T)
!
!      DMS_OH = DMS0 * exp(-(r1+r2)*NDT1)
!          where DMS0 is the DMS concentration at the beginning,
!          r1 = k1*[OH], r2 = k2*[OH]
!
!   R3:    DMS + NO3 ->   SO2 + ...
!          k3 = 1.9e-13*exp(500/T)
!
!      DMS = DMS_OH * exp(-r3*NDT1)
!          where r3 = k3*[NO3]
!
!   R4:    DMS + X   ->   SO2 + ...
!          assume to be at the rate of DMS+OH and DMS+NO3 combined.
!
!   The production of SO2 and MSA here, PSO2_DMS and PMSA_DMS, are saved
!   for use in CHEM_SO2 and CHEM_MSA subroutines as a source term.  They
!   are in unit of MixingRatio/second.
!
! !REVISION HISTORY:
!
!  06Nov2003, Colarco
!  Based on Ginoux
!
!EOP
!-------------------------------------------------------------------------

! !Local Variables
   integer :: i, j, k
   real*8  :: Fx, b, eff
   real*8  :: rk1, rk2, rk3, rk4
   real*8  :: tk, o2, oh, no3, air
   real*8  :: dms, dms0, dms_oh

   data Fx  / 1.0 /
   data b   / 0.25 /
   data eff / 1. /

   _UNUSED_DUMMY(delp)
   _UNUSED_DUMMY(drydepf)

!  spatial loop 
   do k = 1, km
    do j = j1, j2
     do i = i1, i2

      rk1 = 0.
      rk2 = 0.
      rk3 = 0.
      rk4 = 0.

      tk  = tmpu(i,j,k)
      oh  = xoh(i,j,k)
!     air molecules in # cm-3
      air = 1000.*rhoa(i,j,k) / airMolWght * nAvogadro * 1.e-6
!     oxygen molecules in # cm-3
      o2 = 0.21 * air
!     no3 -> go from volume mixing ratio to # cm-3
      no3 = xno3(i,j,k) * air

!     initial DMS concentration (kg kg-1)
      dms0 = qa(i,j,k)
      dms0 = max(dms0,tiny(dms0))

!     1 & 2) DMS + OH: RK1 = addition, RK2 = abstraction
      if(oh .gt. 0.) then
       rk1 = (1.7d-42 * exp(7810./tk) * o2) / &
             (1. + 5.5e-31 * exp(7460./tk) * o2) * oh
       rk2 = (1.2e-11 * exp(-260./tk)) * oh
      endif

!     3) DMS + NO3: only happens at night
      if(cossza(i,j) .le. 0.) then
       rk3 = (1.9e-13 * exp(500./tk)) * no3
      endif

!     Now do the DMS loss
      dms_oh = dms0   * exp( -(rk1+rk2)* Fx * cdt)
      dms    = dms_oh * exp( -(rk3)    * Fx * cdt)

!     SO2 and MSA production terms
!     MSA is formed from the DMS+OH addition step
!     Production should go as mass mixing ratio change in MSA
      if( (rk1+rk2) .eq. 0.) then
       pMSA_DMS(i,j,k) = 0.
      else
       pMSA_DMS(i,j,k) =  (dms0 - dms_oh) * b*rk1/((rk1+rk2)*Fx) * eff &
                         * (fMassMSA/fMassDMS) / cdt
      endif

!     Everything else goes into SO2 formation step
      pSO2_DMS(i,j,k) = ( dms0 - dms - &
                          pMSA_DMS(i,j,k)*cdt*(fMassDMS/fMassMSA) &
                        ) * (fMassSO2/fMassDMS) / cdt


!     4) Dry deposition of DMS (not in GOCART?)
!      if(k .eq. km) rk4 = drydepf(i,j)
!      dms0 = dms
!      dms  = dms0 * exp(-rk4*cdt)
!      dms    = max(dms,1.e-32)

!     Update the mass mixing ratio and the dry deposition flux out of DMS
      dms    = max(dms,tiny(dms))
      qa(i,j,k) = dms

     end do ! i
    end do  ! j
    if(k .eq. km .and. associated(fluxout%data2d) ) fluxout%data2d = 0.
   end do   ! k


   rc = 0

   end subroutine SU_ChemDrv_DMS


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  SU_ChemDrv_SO2
!  NOTE: This is the SO2 oxidation subroutine.

!
! !INTERFACE:
!

   subroutine SU_ChemDrv_SO2( i1, i2, j1, j2, km, cdt, xoh, xh2o2, rhoa,&
                              qa, delp, tmpu, cloud, pSO2_DMS, pSO4g_SO2, &
                              pSO4aq_SO2, drydepf, oro, fluxOut, rc) 

! !USES:

  implicit NONE

! !INPUT PARAMETERS:

   integer, intent(in)       :: i1, i2, j1, j2, km
   real, intent(in)          :: cdt
   real, pointer                    :: qa(:,:,:)
   real, pointer                    :: delp(:,:,:)
   real, pointer, dimension(:,:,:)  :: tmpu, cloud, rhoa
   real, pointer, dimension(:,:)    :: oro
   real, intent(in)                 :: drydepf(i1:i2,j1:j2)
   real, intent(in)                 :: pSO2_DMS(i1:i2,j1:j2,km)
   real, intent(inout)              :: xoh(i1:i2,j1:j2,km), xh2o2(i1:i2,j1:j2,km)

! !OUTPUT PARAMETERS:

   type(Chem_Array), intent(inout)  :: fluxout     ! Mass lost by deposition
                                                   ! to surface, kg/m2/s
   real, intent(out)                :: pSO4g_SO2(i1:i2,j1:j2,km)
   real, intent(out)                :: pSO4aq_SO2(i1:i2,j1:j2,km)
   integer, intent(out)             :: rc          ! Error return code:
                                                   !  0 - all is well
                                                   !  1 - 
   character(len=*), parameter :: myname = 'SU_ChemDrv_SO2'

! !DESCRIPTION: Computes the concentration of SO2 and production of SO4
!
!  SO2 production:
!    DMS + OH, DMS + NO3 (saved in SU_ChemDrv_DMS)
!
!  SO2 loss:
!    SO2 + OH  -> SO4
!    SO2       -> drydep
!    SO2 + H2O2 or O3 (aq) -> SO4
!
!  SO2 = SO2_0 * exp(-bt)
!      + PSO2_DMS*dt/bt * [1-exp(-bt)]
!    where b is the sum of the reaction rate of SO2 + OH and the dry
!    deposition rate of SO2, PSO2_DMS is SO2 production from DMS in
!    MixingRatio/timestep.
!
!  If there is cloud in the gridbox (fraction = fc), then the aqueous
!  phase chemistry also takes place in cloud. The amount of SO2 oxidized
!  by H2O2 in cloud is limited by the available H2O2; the rest may be
!  oxidized due to additional chemistry, e.g, reaction with O3 or O2
!  (catalyzed by trace metal).
!
! !REVISION HISTORY:
!
!  15Jul2010, Colarco - modularized
!  06Nov2003, Colarco
!  Based on Ginoux
!
!EOP
!-------------------------------------------------------------------------

! !Local Variables
   integer :: i, j, k
   real*8  :: rk1, rk2, rk, rkt, f1
   real*8  :: L1, L2, Ld, SO2, SO2_cd, fc, fMR
   real*8  :: oh, h2o2, SO20, tk, air, k0, ki, kk
   real, dimension(i1:i2,j1:j2) :: fout

   data ki / 1.5e-12 /

!  Conversion of SO2 mmr to SO2 vmr
   fMR = airMolWght / fMassSO2

!  Initialize flux variable   
   fout = 0.

!  spatial loop 
   do k = 1, km
    do j = j1, j2
     do i = i1, i2

      rk1 = 0.
      rk2 = 0.
      L1  = 0.
      L2  = 0.
      Ld  = 0.

      tk   = tmpu(i,j,k)
      oh   = xoh(i,j,k)
      h2o2 = max(xh2o2(i,j,k),tiny(xh2o2(i,j,k)))

!     air molecules in # cm-3
      air  = 1000.*rhoa(i,j,k) / airMolWght * nAvogadro * 1.e-6
!     1) SO2 + OH(g) in s-1
      k0 = 3.0e-31 * (300./tk)**4.3
      kk = k0 * air / ki
      f1 = (1. + (log10(kk))**2.)**(-1.)
      rk1 = ( (k0*air/(1.+kk)) * 0.6**f1) * oh

!     2) rk2 is the loss of SO2 due to dry deposition.
      if(k .eq. km) then
!      drydepf calculated for aerosol
!      follow Walcek: ocean drydepf_so2 = 10*drydepf_aer
!      or if land drydepf_so2 = 3*drydepf_aer
       if(oro(i,j) .eq. OCEAN) then
        rk2 = 10.*drydepf(i,j)
       else
        rk2 = 3.*drydepf(i,j)
       endif
!       rk2 = drydepf(i,j)
      else
       rk2 = 0.
      endif

      rk = (rk1 + rk2)
      rkt = rk*cdt

!     Update the SO2 concentration
!     Originally this was solved like a simple exponential solution
!     after Jacobson eq. 13.38, which is more accurate but not mass
!     conserving.  We've already timesplit everything, so accuracy is
!     out to lunch, and I'd prefer to conserve mass.

!     initial SO2 concentration (kg kg-1) after adding source
      SO20 = qa(i,j,k) + pSO2_DMS(i,j,k)*cdt
      SO20 = max(SO20,tiny(SO20))

      if(rk .gt. 0.) then
       SO2_cd =  SO20 * exp(-rkt)
       L1     = (SO20 - SO2_cd) * rk1/rk
       if(k .eq. km) then
        Ld    = (SO20 - SO2_cd) * rk2/rk
        fout(i,j) = Ld * delp(i,j,km)/grav/cdt
       else
        Ld    = 0.
       endif
      else
       SO2_cd = SO20
       L1     = 0.
      endif


!     Update SO2 concentration after cloud chemistry, if it occurs
      fc = cloud(i,j,k)
      if(fc .gt. 0. .and. SO2_cd .gt. 0. .and. tk .gt. 258.) then
!      Check on H2O2 vmr -> is SO2 vmr greater?
       if(fMr * SO2_cd .gt. h2o2) then
        fc = fc*(h2o2/(fMR*SO2_cd))
        h2o2 = h2o2*(1.-cloud(i,j,k))
       else
        h2o2 = h2o2*(1. - cloud(i,j,k)*(fMR*SO2_cd)/h2o2)
       endif
       SO2 = SO2_cd*(1.-fc)
!      aqueous loss rate (mixing ratio/timestep)
       L2 = SO2_cd * fc
      else
       SO2 = SO2_cd
       L2 = 0.
      endif

!     Ideally you would update the H2O2 mixing ratio at this point
!     and then reset it periodically
      xh2o2(i,j,k) = max(h2o2,tiny(h2o2))

      SO2 = max(SO2,tiny(SO2))
      qa(i,j,k) = SO2
      pSO4g_SO2(i,j,k) = L1 * (fMassSO4/fMassSO2) / cdt
      pSO4aq_SO2(i,j,k) = L2 * (fMassSO4/fMassSO2) / cdt

     end do
    end do
   end do

   if( associated(fluxout%data2d) ) fluxout%data2d = fout

   rc = 0

   end subroutine SU_ChemDrv_SO2

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  SU_ChemDrv_SO4 - Do SU cycle chemistry following GOCART
!  NOTE: This is the SO4 production by SO2 oxidation subroutine.

!
! !INTERFACE:
!

   subroutine SU_ChemDrv_SO4( i1, i2, j1, j2, km, cdt, qa, delp, &
                              pSO4g_SO2, pSO4aq_SO2, drydepf, fluxOut, rc) 

! !USES:

  implicit NONE

! !INPUT PARAMETERS:

   integer, intent(in)       :: i1, i2, j1, j2, km
   real, intent(in)          :: cdt
   real, pointer             :: qa(:,:,:)
   real, pointer             :: delp(:,:,:)
   real, intent(in)          :: drydepf(i1:i2,j1:j2)
   real, intent(in)          :: pSO4g_SO2(i1:i2,j1:j2,km)
   real, intent(in)          :: pSO4aq_SO2(i1:i2,j1:j2,km)

! !OUTPUT PARAMETERS:

   type(Chem_Array), intent(inout)  :: fluxout     ! Mass lost by deposition
                                                   ! to surface, kg/m2/s
   integer, intent(out)             :: rc          ! Error return code:
                                                   !  0 - all is well
                                                   !  1 - 
   character(len=*), parameter :: myname = 'SU_ChemDrv_SO4'

! !DESCRIPTION:
!
!  SO4 production:
!    The only production term is due to SO2 oxidation.
!    SO4 = SO4_0 * exp(-kt) + pSO4_SO2/kt * (1.-exp(-kt))
!     where k is the dry deposition
!
! !REVISION HISTORY:
!
!  15Jul2010, Colarco - Modularized
!  06Nov2003, Colarco
!  Based on Ginoux
!
!EOP
!-------------------------------------------------------------------------

! !Local Variables
   integer :: i, j, k
   real*8  :: rk, rkt, Ld
   real*8  :: SO4, SO40, pSO4
   real, dimension(i1:i2,j1:j2) :: fout


!  Initialize flux variable
   fout = 0.

!  spatial loop 
   do k = 1, km
    do j = j1, j2
     do i = i1, i2

      pSO4 = pSO4g_SO2(i,j,k)+pSO4aq_SO2(i,j,k)

!     initial SO4 concentration (kg kg-1)
      SO40 = qa(i,j,k)
      SO40 = max(SO40,tiny(SO40))

!     Update the SO4 concentration
!     Originally this was solved like a simple exponential solution
!     after Jacobson eq. 13.38, which is more accurate but not mass
!     conserving.  We've already timesplit everything, so accuracy is
!     out to lunch, and I'd prefer to conserve mass.
!     RK is the dry deposition frequency
      if(k .eq. km) then
       RK = drydepf(i,j)
       RKT = RK*cdt
       SO4 = (SO40 + pSO4*cdt) * exp(-rkt)
       Ld  = (SO40 - SO4 + pSO4*cdt)
       fout(i,j) = Ld * delp(i,j,km)/grav/cdt
      else
       SO4 = SO40 + pSO4*cdt
       Ld = 0.
      endif

      SO4 = max(SO4,tiny(SO4))
      qa(i,j,k) = SO4

     end do
    end do
   end do

   if( associated(fluxout%data2d) ) fluxout%data2d = fout

   rc = 0

   end subroutine SU_ChemDrv_SO4


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  SU_ChemDrv_MSA - MSA production by DMS oxidation

!
! !INTERFACE:
!

   subroutine SU_ChemDrv_MSA( i1, i2, j1, j2, km, cdt, qa, delp, &
                              pMSA_DMS, drydepf, fluxOut, rc) 

! !USES:

  implicit NONE

! !INPUT PARAMETERS:

   integer, intent(in)       :: i1, i2, j1, j2, km
   real, intent(in)          :: cdt
   real, pointer             :: qa(:,:,:)
   real, pointer             :: delp(:,:,:)
   real, intent(in)          :: drydepf(i1:i2,j1:j2)
   real, intent(in)          :: pMSA_DMS(i1:i2,j1:j2,km)

! !OUTPUT PARAMETERS:

   type(Chem_Array), intent(inout)  :: fluxout  ! Mass lost by deposition
                                                ! to surface, kg/m2/s
   integer, intent(out)             :: rc       ! Error return code:
                                                !  0 - all is well
                                                !  1 - 
   character(len=*), parameter :: myname = 'SU_ChemDrv_MSA'

! !DESCRIPTION: 
!
!  MSA production:
!    The only production term is due to DMS oxidation.
!    MSA = MSA_0 * exp(-kt) + pSO4_SO2/kt * (1.-exp(-kt))
!     where k is the dry deposition
!
! !REVISION HISTORY:
!
!  15Jul2010, Colarco -- modularized
!  06Nov2003, Colarco
!  Based on Ginoux
!
!EOP
!-------------------------------------------------------------------------

! !Local Variables
   integer :: i, j, k
   real*8  :: rk, rkt, Ld
   real*8  :: MSA, MSA0
   real, dimension(i1:i2,j1:j2) :: fout

!  Initialize flux variable
   fout = 0.

!  spatial loop 
   do k = 1, km
    do j = j1, j2
     do i = i1, i2

!     initial MSA concentration (kg kg-1)
      MSA0 = qa(i,j,k)
      MSA0 = max(MSA0,tiny(MSA0))

!     Update the MSA concentration
!     Originally this was solved like a simple exponential solution
!     after Jacobson eq. 13.38, which is more accurate but not mass
!     conserving.  We've already timesplit everything, so accuracy is
!     out to lunch, and I'd prefer to conserve mass.
!     RK is the dry deposition frequency
      if(k .eq. km) then
       RK = drydepf(i,j)
       RKT = RK*cdt
       MSA = (MSA0 + pMSA_DMS(i,j,k)*cdt) * exp(-rkt)
       Ld  = (MSA0 + pMSA_DMS(i,j,k)*cdt - MSA)
       fout(i,j) = Ld * delp(i,j,km)/grav/cdt
      else
       MSA = MSA0 + pMSA_DMS(i,j,k)*cdt
       Ld = 0.
      endif

      MSA = max(MSA,tiny(MSA))
      qa(i,j,k) = MSA

     end do
    end do
   end do

   if( associated(fluxout%data2d) ) fluxout%data2d = fout

   rc = 0

   end subroutine SU_ChemDrv_MSA

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  SulfateUpdateOxidants - Update Oxidant Fields for Sulfate
!             We have 3 oxidant fields (OH, NO3, H2O2) which may come
!             from either a climatological file or from interactive GMI.
!             IF from climatology, update (reset) values from climatology
!             if necessary (e.g., for a new day) and set to current values
!             needed by chemistry.
!             IF from GMI read as is
!
! !INTERFACE:
   subroutine SulfateUpdateOxidants ( impChem,iName, i1, i2, im, j1, j2, jm, km, cdt, &
                                      using_GMI_OH, using_GMI_NO3, &
                                      using_GMI_H2O2, &
                                      GMI_OHmr, GMI_NO3mr, GMI_H2O2mr, &
                                      nymd_current, nhms_current, &
                                      grid, lonRad, latRad, &
                                      rhoa, &
                                      nymd_last, &
                                      oh_clim, no3_clim, h2o2_clim, &
                                      xoh, xno3, xh2o2, recycle_h2o2 )

!
!  Input
   type(ESMF_State), intent(inout) :: impChem
   character(len=*), intent(in   ) :: iName
   integer, intent(in)    :: i1, i2, im, j1, j2, jm, km
   real, intent(in)       :: cdt
   integer, intent(in)    :: nymd_current, &   ! current model NYMD
                             nhms_current      ! current model NHMS
   integer, intent(inout) :: nymd_last         ! NYMD of last emission update
   logical, intent(inout) :: recycle_h2o2      ! triggers recycling of H2O2
   logical, intent(in)    :: using_GMI_OH, &
                             using_GMI_NO3, &
                             using_GMI_H2O2
   type(ESMF_Grid)  :: grid
   real, pointer, dimension(:,:) :: lonRad, latRad

   real, pointer, dimension(:,:,:) :: oh_clim, &
                                      no3_clim, &
                                      h2o2_clim, &
                                      xoh, xno3, xh2o2, &
                                      GMI_OHmr, &
                                      GMI_H2O2mr, &
                                      GMI_NO3mr, &
                                      rhoa

!  Local
   integer :: i, j, k, jday
   integer :: STATUS
   real    :: qmax, xhour, xhouruse
   real    :: cossza(i1:i2,j1:j2), sza(i1:i2,j1:j2)
   real    :: tcosz(i1:i2,j1:j2), tday(i1:i2,j1:j2), tnight(i1:i2,j1:j2)
   integer :: n, ndystep
   real, pointer :: ptr3d(:,:,:) => null()

   _UNUSED_DUMMY(im)
   _UNUSED_DUMMY(jm)
   _UNUSED_DUMMY(grid)


! Update emissions/production if necessary (daily)
!  -----------------------------------------------
!   Oxidant fields
!   The expectation here is that OH is being read in the form
!   volume mixing ratio from a file (so, like GMI would provide).
!   Below, in the scaling by solar zenith angle, we convert from
!   VMR to # cm-3 expected by the chemistry.
    IF(.NOT. using_GMI_OH) THEN
     call MAPL_GetPointer(impChem,ptr3d,"SU_OH"//trim(iName),rc=status)
     oh_clim = ptr3d

     where(1.01*oh_clim(i1:i2,j1:j2,1:km) > undefval) oh_clim(i1:i2,j1:j2,1:km) = 0.
     where(     oh_clim(i1:i2,j1:j2,1:km) < 0       ) oh_clim(i1:i2,j1:j2,1:km) = 0.
    END IF


    IF(.NOT. using_GMI_NO3) THEN
     call MAPL_GetPointer(impChem,ptr3d,"SU_NO3"//trim(iName),rc=status)
     no3_clim = ptr3d

     where(1.01*no3_clim(i1:i2,j1:j2,1:km) > undefval) no3_clim(i1:i2,j1:j2,1:km) = 0.
     where(     no3_clim(i1:i2,j1:j2,1:km) < 0       ) no3_clim(i1:i2,j1:j2,1:km) = 0.
    END IF

    IF(.NOT. using_GMI_H2O2) THEN
     call MAPL_GetPointer(impChem,ptr3d,"SU_H2O2"//trim(iName),rc=status)
     h2o2_clim = ptr3d

     where(1.01*h2o2_clim(i1:i2,j1:j2,1:km) > undefval) h2o2_clim(i1:i2,j1:j2,1:km) = 0.
     where(     h2o2_clim(i1:i2,j1:j2,1:km) < 0       ) h2o2_clim(i1:i2,j1:j2,1:km) = 0.
    END IF

!   The first time through the reads we will save the h2o2 monthly
!   average in the instantaneous field
!   ---------------------------------
    if (nymd_last == nymd_current .and. (.not. using_GMI_H2O2)) then
     xh2o2 = h2o2_clim
     nymd_last = nymd_current
if(mapl_am_i_root()) print*,'SU nymd_oxidants == nymd_current'
    end if 


!  Find the day number of the year and hour (needed for later doing sza)
!  ----------------------------------
   jday = idaynum(nymd_current)
   xhour = (  real(nhms_current/10000)*3600. &
            + real(mod(nhms_current,10000)/100)*60. &
            + real(mod(nhms_current,100)) &
           ) / 3600.

!  Recycle H2O2 to input on 3 hour boundaries if not coupled to GMI
!  ----------------------------------
   if (.NOT. using_GMI_H2O2 .and. recycle_h2o2) then
        xh2o2 = h2o2_clim
        recycle_h2o2 = .false.
   end if

!  If not getting instantaneous values from GMI, update for time of day.
!  ---------------------------------------------------------------------
!  OH
   if( .not. using_GMI_OH) then
    xoh = oh_clim
    cossza(:,:) = 0.

!   Want to find the sum of the cos(sza) for use in scaling OH diurnal variation
!   tcosz is the sum of cossza over the whole day
!   tday is the time of day spent in light
!   Requires integrating over future times, so cannot use w_c%cosz
    xHourUse = xHour
    ndystep = 86400. / cdt
    tcosz(:,:) = 0.
    tday(:,:) = 0.
    do n = 1, ndystep
     call szangle(jday, xHourUse, lonRad, latRad, sza, cossza, i1, i2, j1, j2)
     tcosz = tcosz + cossza
     xHourUse = xHourUse + cdt/3600.
     if(xHourUse .gt. 24.) xHourUse = xHourUse - 24.
!    Find the daylight portion of the day
     do j = j1, j2
      do i = i1, i2
       if(cossza(i,j) .gt. 0.) tday(i,j) = tday(i,j) + cdt
      end do
     end do
    end do

!   Find the cos(sza) now for use in scaling OH and NO3
    call szangle(jday,xHour,lonRad,latRad,sza,cossza,i1,i2,j1,j2)

    tnight(i1:i2,j1:j2) = (86400.-tday(i1:i2,j1:j2))

    DO k=1,km
     WHERE(tcosz(i1:i2,j1:j2) > 0)
      xoh(i1:i2,j1:j2,k) = oh_clim(i1:i2,j1:j2,k)*(86400./cdt)*cossza(i1:i2,j1:j2) / tcosz(i1:i2,j1:j2)
     ELSEWHERE
      xoh(i1:i2,j1:j2,k) = 0.00
     END WHERE
    END DO
    WHERE(xoh(i1:i2,j1:j2,1:km) < 0.00) xoh(i1:i2,j1:j2,1:km) = 0.00
   endif
!  To go from volume mixing ratio to # cm-3 (expected in chemistry)
!  include the following line
   xoh = xoh * 1000.*rhoa / airMolWght * nAvogadro * 1.e-6

!  NO3
   IF(.NOT. using_GMI_NO3) THEN
    xno3 = no3_clim
    cossza(:,:) = 0.
    call szangle(jday,xHour,lonRad,latRad,sza,cossza,i1,i2,j1,j2)

!   If there is daylight then no3 is small (assume zero) and the
!   average is distributed only over the night time portion

    DO k=1,km
     WHERE(cossza(i1:i2,j1:j2) > 0 .OR. tnight(i1:i2,j1:j2) < tiny(1.0))
      xno3(i1:i2,j1:j2,k) = 0.00
     ELSEWHERE
      xno3(i1:i2,j1:j2,k) = no3_clim(i1:i2,j1:j2,k) * 86400./ tnight(i1:i2,j1:j2)
     END WHERE
    END DO
   END IF

!  If doing GMI, grab oxidants from GMICHEM if the pointers were found.
!  Note: OH units must be cm^{-3}.
!  --------------------------------------------------------------------
   IF( using_GMI_NO3) THEN
    xno3(i1:i2,j1:j2,1:km)  = GMI_NO3mr(i1:i2,j1:j2,1:km)
    WHERE(xno3(i1:i2,j1:j2,1:km) < 0.00) xno3(i1:i2,j1:j2,1:km) = 0.00
   END IF

   IF(using_GMI_H2O2) THEN
    xh2o2(i1:i2,j1:j2,1:km) = GMI_H2O2mr(i1:i2,j1:j2,1:km)
    WHERE(xh2o2(i1:i2,j1:j2,1:km) < 0.00) xh2o2(i1:i2,j1:j2,1:km) = 0.00
   END IF

   IF(  using_GMI_OH) THEN
    qmax = 17.01/airMolWght
    xoh(i1:i2,j1:j2,1:km) =  GMI_OHmr(i1:i2,j1:j2,1:km)* &
                             nAvogadro / airMolWght * 1000.* &
                             rhoa(i1:i2,j1:j2,1:km)*1.00E-06
    WHERE(xoh(i1:i2,j1:j2,1:km) < 0.00) xoh(i1:i2,j1:j2,1:km) = 0.00
   END IF




   end subroutine SulfateUpdateOxidants




!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  SulfateUpdateEmissions - Update Sulfate Emissions
!             We have emissions from 4 sources, which are distributed
!             differently in the vertical
!             1) biomass burning - uniformly mixed in PBL (SO2)
!             2) anthropogenic l1 - emitted into lowest 100 m (SO2,SO4)
!             3) anthropogenic l2 - emitted into 100 - 500 m levels (SO2,SO4)
!             4) volcanic emissions
!             Additionally have a source of DMS from transfer from seawater
!             into lowest model layer
!             Consider factors in conversion: we estimate that 5% of sulfur
!             from anthropogenic sources (by mass) goes directly to SO4.
!
! !INTERFACE:
   subroutine SulfateUpdateEmissions ( impChem, iName, i1, i2, im, j1, j2, jm, km, cdt, &
                                       nymd_current, nhms_current, &
                                       grid, lonRad, latRad, &
                                       nymd_last, &
                                       diurnal_bb, &
                                       so2biomass_src, so2biomass_src_, &
                                       so2anthro_l1_src, &
                                       so2anthro_l2_src, &
                                       so2ship_src, &
                                       so4ship_src, &
                                       dmso_conc, &
                                       aircraft_fuel_src, &
                                       aviation_lto_src, &
                                       aviation_cds_src, &
                                       aviation_crs_src, &
                                       volcano_srcfilen, &
                                       nvolc, vLat, vLon, vElev, vCloud, vSO2, vStart, vEnd, &
                                       doing_NEI, nei_hour, nei_year, nei_srcfilen, &
                                       nei_lon, nei_lat, lons, lats, &
                                       maskString, gridMask, &
                                       rc)

!
!  Input
   type(ESMF_State), intent(inout) :: impChem
   character(len=*), intent(in   ) :: iName
   integer :: i1, i2, im, j1, j2, jm, km
   real    :: cdt
   integer :: nymd_current, &   ! current model NYMD
              nhms_current, &   ! current model NHMS
              nymd_last         ! NYMD of last emission update
   logical :: diurnal_bb
   type(ESMF_Grid)  :: grid
   real, pointer, dimension(:,:) :: lonRad, latRad

   character(len=*) :: volcano_srcfilen

   real, pointer, dimension(:,:) :: so2biomass_src, &
                                    so2biomass_src_, &
                                    so2anthro_l1_src, &
                                    so2anthro_l2_src, &
                                    so2ship_src, &
                                    so4ship_src, &
                                    dmso_conc

   real, pointer, dimension(:,:,:) :: aircraft_fuel_src
   real, pointer, dimension(:)     :: vLat, vLon, vElev, vCloud, vSO2
   integer, pointer, dimension(:)  :: vStart, vEnd
   real, pointer, dimension(:)     :: vLatE   => null(), &
                                      vLonE   => null(), &
                                      vElevE  => null(), &
                                      vCloudE => null(), &
                                      vSO2E   => null()

   real, pointer, dimension(:)     :: vLatC   => null(), &
                                      vLonC   => null(), &
                                      vElevC  => null(), &
                                      vCloudC => null(), &
                                      vSO2C   => null()

   integer                         :: nVolc, nVolcE, nVolcC
   logical                         :: useVolcanicDailyTables = .false.

!  Input parameters for NEI08 handling
!  Note: these should are being defined as optional until the CARMA API is updated.
   logical, intent(in), optional :: doing_NEI
   integer, intent(inout), optional :: nei_hour, nei_year
   character(len=*), intent(in), optional :: nei_srcfilen(2)
   real, intent(in), optional :: nei_lon(2), nei_lat(2), lons(:,:), lats(:,:)

!  Emissions from aviation sector
!  ------------------------------
   real, dimension(:,:), intent(inout) :: aviation_lto_src, &
                                          aviation_cds_src, &
                                          aviation_crs_src

   integer, intent(out) :: rc

!  Optional parameters to be passed for masking separate instances
   character(len=*), OPTIONAL, intent(in) :: maskString             !Delimited string of integers
   real, OPTIONAL, intent(in)             :: gridMask(i1:i2,j1:j2)  !Grid mask (NOTE: No ghosting) 

! Local
  integer :: ijl, ijkl
  integer :: i, ios
  real, pointer :: ptr2d(:,:)   => null()
  real, pointer :: ptr3d(:,:,:) => null()
  integer       :: status

! Workspace for NEI emissions
! ---------------------------
  real, pointer, dimension(:,:)         ::  nei_src1, nei_src2
  integer                               ::  hour, year, nei_nymd
#ifdef DEBUG
   real :: qmin, qmax
#endif


  character(len=32) :: Iam

  _UNUSED_DUMMY(maskString)
  _UNUSED_DUMMY(gridMask)

  Iam = 'SulfateUpdateEmissions'

  ijl  = ( i2 - i1 + 1 ) * ( j2 - j1 + 1 )
  ijkl = ijl * km

! Special handling of the volcanic sources, which are provided only in a table
! based format at the present time.  The key string searched on here is
! "volcanic_".  If present this refers to the daily tables.  Anything else
! points to the data statement tables in this piece of Fortran code.  The 
! special case is "/dev/null" which is to run without volcanoes.
  if(index(volcano_srcfilen,'volcanic_') .ne. 0) useVolcanicDailyTables = .true.
  

!   Anthropogenic and ship emissions are now outside of UpdateEmiss loop to 
!   account for diurnal variability in HEMCO emissions (ckeller, 2/27/16)

!   Anthropogenic emissions
!   -----------------------
    call MAPL_GetPointer(impChem,ptr2d,"SU_ANTHROL1"//trim(iName),rc=status)
    VERIFY_(STATUS)
    so2anthro_l1_src = ptr2d

    call MAPL_GetPointer(impChem,ptr2d,"SU_ANTHROL2"//trim(iName),rc=status)
    VERIFY_(STATUS)
    so2anthro_l2_src = ptr2d

!   Ship based emissions of SO2 and SO4
!   -----------------------------------
    call MAPL_GetPointer(impChem,ptr2d,"SU_SHIPSO2"//trim(iName),rc=status)
    VERIFY_(STATUS)
    so2ship_src = ptr2d

    call MAPL_GetPointer(impChem,ptr2d,"SU_SHIPSO4"//trim(iName),rc=status)
    VERIFY_(STATUS)
    so4ship_src = ptr2d

! Update emissions/production if necessary (daily)
!  -----------------------------------------------
if(mapl_am_i_root()) print*,'SU nymd_last = ', nymd_last
if(mapl_am_i_root()) print*,'SU nymd_current = ', nymd_current

   UpdateEmiss: if(nymd_last .ne. nymd_current) then
    nymd_last = nymd_current

if(mapl_am_i_root()) print*,'SU inside UpdateEmiss'

!   Biomass Burning -- select on known inventories
!   ----------------------------------------------
    call MAPL_GetPointer(impChem,ptr2d,"SU_BIOMASS"//trim(iName),rc=status)
    so2biomass_src = ptr2d

!   Save read in emissions if doing diurnal BB
!   ------------------------------------------
    if ( diurnal_bb ) then
         so2biomass_src_(:,:) = so2biomass_src(:,:)
    end if

!   DMS concentrations (from climatology)
!   -------------------------------------
    call MAPL_GetPointer(impChem,ptr2d,"SU_DMSO"//trim(iName),rc=status)
    VERIFY_(STATUS)
    dmso_conc = ptr2d

!   Aircraft fuel source
!   --------------------
    call MAPL_GetPointer(impChem,ptr3d,"SU_AIRCRAFT"//trim(iName),rc=status)
    VERIFY_(STATUS)
    aircraft_fuel_src = ptr3d

!   Aviation-LTO emissions
!   ----------------------
    call MAPL_GetPointer(impChem,ptr2d,'SU_AVIATION_LTO'//iNAME,rc=status)
    VERIFY_(STATUS)
    aviation_lto_src = ptr2d

    call MAPL_GetPointer(impChem,ptr2d,'SU_AVIATION_CDS'//iNAME,rc=status)
    VERIFY_(STATUS)
    aviation_cds_src = ptr2d

    call MAPL_GetPointer(impChem,ptr2d,'SU_AVIATION_CRS'//iNAME,rc=status)
    VERIFY_(STATUS)
    aviation_crs_src = ptr2d

!   As a safety check, where values are undefined set to 0
    where(1.01*so2biomass_src(i1:i2,j1:j2)  > undefval)     so2biomass_src(i1:i2,j1:j2) = 0.
    where(1.01*dmso_conc(i1:i2,j1:j2)    > undefval)        dmso_conc(i1:i2,j1:j2) = 0.
    where(1.01*so2anthro_l1_src(i1:i2,j1:j2)    > undefval) so2anthro_l1_src(i1:i2,j1:j2) = 0.
    where(1.01*so2anthro_l2_src(i1:i2,j1:j2)    > undefval) so2anthro_l2_src(i1:i2,j1:j2) = 0.
    where(1.01*so2ship_src(i1:i2,j1:j2) > undefval)         so2ship_src(i1:i2,j1:j2) = 0.
    where(1.01*so4ship_src(i1:i2,j1:j2) > undefval)         so4ship_src(i1:i2,j1:j2) = 0.
    where(1.01*aircraft_fuel_src(i1:i2,j1:j2,1:km) > undefval ) &
               aircraft_fuel_src(i1:i2,j1:j2,1:km) = 0.
    where(1.01*aviation_lto_src(i1:i2,j1:j2) > undefval )   aviation_lto_src(i1:i2,j1:j2) = 0.
    where(1.01*aviation_cds_src(i1:i2,j1:j2) > undefval )   aviation_cds_src(i1:i2,j1:j2) = 0.    
    where(1.01*aviation_crs_src(i1:i2,j1:j2) > undefval )   aviation_crs_src(i1:i2,j1:j2) = 0.

#ifdef DEBUG
    call pmaxmin('SU: so2biomass_src  ', so2biomass_src, qmin, qmax, ijl,1, 1.  )
    call pmaxmin('SU: so2anthro_l1_src', so2anthro_l1_src,   qmin, qmax, ijl,1, 1.  )
    call pmaxmin('SU: so2anthro_l2_src', so2anthro_l2_src,   qmin, qmax, ijl,1, 1.  )
    call pmaxmin('SU: so2ship_src     ', so2ship_src, qmin, qmax, ijl,1, 1.  )
    call pmaxmin('SU: so4ship_src     ', so4ship_src, qmin, qmax, ijl,1, 1.  )
    call pmaxmin('SU: DMSO_conc       ', dmso_conc,   qmin, qmax, ijl,1, 1.  )
    call pmaxmin('SU: fuel            ', aircraft_fuel_src, qmin, qmax, ijl,km, 1. )
    call pmaxmin('SU: so2_aviation_lto', aviation_lto_src,  qmin, qmax, ijl,1, 1.  )
    call pmaxmin('SU: so2_aviation_cds', aviation_cds_src,  qmin, qmax, ijl,1, 1.  )
    call pmaxmin('SU: so2_aviation_crs', aviation_crs_src,  qmin, qmax, ijl,1, 1.  )
#endif

!   Volcanic emissions
!   ------------------
!   Two possible sources of volcanic emissions
!    1) daily tables referenced as external text files
!    2) contained inventory of explosive and non-explosive
!       volcanoes
!   What is provided at this point is the number of
!   volcanoes, their locations, elevations, and SO2 emissions.
!   Added some additional variables "vStart" and "vEnd" which
!   for now are the start and end time of the eruption.  For hard
!   wired data tables the emissions are assumed to hold over entire
!   day.  With "useVolcanicDailyTables" there are optionally two
!   columns in each line of table that provide start and end hour
!   of event.  In all cases we are assuming the value of vSO2 is
!   kg SO2 s-1 over the *duration of the event*.  The duration is
!   24 hours unless otherwise specified.  So when updating the daily
!   tables put in numbers accordingly.

    if( useVolcanicDailyTables ) then
     call GetVolcDailyTables( nymd_current, volcano_srcfilen, &
                              nVolc, vLat, vLon, vElev, vCloud, vSO2, vStart, vEnd )
!   Read from the previous inventory of non-explosive volcanoes
!   Special handling to partition is possible (and /dev/null handling)
    else
     if(volcano_srcfilen(1:4) == 'cont') then
        call GetVolcContinuous( nVolcC, vLatC, vLonC, vElevC, vCloudC, vSO2C )
        nVolcE=0
     else if(volcano_srcfilen(1:4) == 'expl') then
        call GetVolcExplosive( nymd_current, nVolcE, vLatE, vLonE, vElevE, vCloudE, vSO2E )
        nVolcC=0
     else
        call GetVolcExplosive( nymd_current, nVolcE, vLatE, vLonE, vElevE, vCloudE, vSO2E )
        call GetVolcContinuous( nVolcC, vLatC, vLonC, vElevC, vCloudC, vSO2C )
     endif

!    combine the two data sets
!    If previous instance of volcano point data tables exist, deallocate it
!    to get the correct number of elements
     if(associated(vLat))    deallocate(vLat, stat=ios)
     if(associated(vLon))    deallocate(vLon, stat=ios)
     if(associated(vSO2))    deallocate(vSO2, stat=ios)
     if(associated(vElev))   deallocate(vElev, stat=ios)
     if(associated(vCloud))  deallocate(vCloud, stat=ios)
     if(associated(vStart))  deallocate(vStart, stat=ios)
     if(associated(vEnd))    deallocate(vEnd, stat=ios)

!    Allocate space for the volcanoes
     nVolc = nVolcE + nVolcC
     allocate(vLat(nvolc), vLon(nvolc), &
              vSO2(nvolc), vElev(nvolc), &
              vCloud(nvolc), vStart(nvolc), vEnd(nvolc), &
              stat=ios)
     if(nVolc > 0) then
      if(nVolcE > 0) then
       do i = 1, nVolcE
        vLat(i) = vLatE(i)
        vLon(i) = vLonE(i)
        vElev(i) = vElevE(i)
        vCloud(i) = vCloudE(i)
        vSO2(i) = vSO2E(i)
       end do
      end if
      if(nVolcC > 0) then
       do i = 1, nVolcC
        vLat(i+nVolcE) = vLatC(i)
        vLon(i+nVolcE) = vLonC(i)
        vElev(i+nVolcE) = vElevC(i)
        vCloud(i+nVolcE) = vCloudC(i)
        vSO2(i+nVolcE) = vSO2C(i)
       end do
      end if
     endif

!    For these tables vStart and vEnd are not provided, so we assume
!    eruption is throughout day and set to default values
     vStart = -1
     vEnd   = -1

!    Clean Up
     if(associated(vLatC))    deallocate(vLatC, stat=ios)
     if(associated(vLonC))    deallocate(vLonC, stat=ios)
     if(associated(vSO2C))    deallocate(vSO2C, stat=ios)
     if(associated(vElevC))   deallocate(vElevC, stat=ios)
     if(associated(vCloudC))  deallocate(vCloudC, stat=ios)
     if(associated(vLatE))    deallocate(vLatE, stat=ios)
     if(associated(vLonE))    deallocate(vLonE, stat=ios)
     if(associated(vSO2E))    deallocate(vSO2E, stat=ios)
     if(associated(vElevE))   deallocate(vElevE, stat=ios)
     if(associated(vCloudE))  deallocate(vCloudE, stat=ios)
!    Special possible case
     if(volcano_srcfilen(1:9) == '/dev/null') nvolc = 0
    endif

!   For volcanos, check value of vStart and vEnd.  Set to be
!   vStart = 000000 if default (=-1) is provided
!   vEnd   = 240000 if default (=-1) is provided
    where(vStart < 0) vStart = 000000
    where(vEnd < 0)   vEnd   = 240000


  endif UpdateEmiss

! Apply dirunal emissions to BB
! ------------------------------------------
  if ( diurnal_bb ) then
       call Chem_BiomassDiurnal ( so2biomass_src, so2biomass_src_,   &
                                  lonRad*radToDeg, latRad*radToDeg, nhms_current, cdt )      
if(mapl_am_i_root()) print*,'SU inside sum(so2biomass_src) = ',sum(so2biomass_src)
if(mapl_am_i_root()) print*,'SU inside sum(so2biomass_src_) = ',sum(so2biomass_src_)
  end if

!  Apply NEI emissions over North America if so desired
!  ----------------------------------------------------
   if ( present(doing_nei) ) then
    if (doing_NEI) then

       hour = nhms_current/10000
       year = nymd_current/10000

       if ( hour /= nei_hour ) then

            allocate(nei_src1(i1:i2,j1:j2),nei_src2(i1:i2,j1:j2),stat=ios)

!           Handle SO2
!           ----------
            call MAPL_GetPointer(impChem,ptr2d,'SU_NEI_SRC1',rc=status)
            VERIFY_(STATUS)
            call MAPL_GetPointer(impChem,ptr2d,'SU_NEI_SRC2',rc=status)
            VERIFY_(STATUS)

            WHERE ( (lons.ge.nei_lon(1)) .AND. &
                    (lons.le.nei_lon(2)) .AND. &
                    (lats.ge.nei_lat(1)) .AND. &
                    (lats.le.nei_lat(2))   )

                    so2anthro_l1_src = nei_src1
                    so2anthro_l2_src = nei_src2

            end where

#ifdef DEBUG
            call pmaxmin('SO2: nei_bot', nei_src1, qmin, qmax, ijl,1, 1. )
            call pmaxmin('SO2: nei_src', nei_src2, qmin, qmax, ijl,1, 1. )
#endif

            nei_hour = hour ! only update NEI once hourly
            deallocate(nei_src1,nei_src2)

       end if ! time to update NEI

    end if ! doing NEI
   end if ! present(doing_nei)

  rc = 0

  end subroutine SulfateUpdateEmissions


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  SulfateDistributeEmissions - Adds sulfate source emission for one timestep
!             We have emissions from 4 sources, which are distributed
!             differently in the vertical
!             1) biomass burning - uniformly mixed in PBL (SO2)
!             2) anthropogenic l1 - emitted into lowest 100 m (SO2,SO4)
!             3) anthropogenic l2 - emitted into 100 - 500 m levels (SO2,SO4)
!             4) volcanic emissions
!             Additionally have a source of DMS from transfer from seawater
!             into lowest model layer
!             Consider factors in conversion: we estimate that 5% of sulfur
!             from anthropogenic sources (by mass) goes directly to SO4.
!
! !INTERFACE:
!

   subroutine SulfateDistributeEmissions ( i1, i2, j1, j2, km, nbins, cdt, nymd, nhms, &
                                           fSO4ant, eAircraftFuel, &
                                           so2anthro_l1_src, so2anthro_l2_src, &
                                           so2biomass_src, dmso_conc, &
                                           so2ship_src, so4ship_src, &
                                           aircraft_fuel_src, &
                                           nvolc, vlat, vlon, velev, vcloud, &
                                           vso2, vstart, vend, &
                                           dms, so2, so4, &
                                           oro, u10m, v10m, hsurf, hghte, pblh, &
                                           tmpu, rhoa, delp, &
                                           cell_area, grid, &
                                           SU_emis, &
                                           SU_SO4eman, SU_SO2eman, SU_SO2embb, &
                                           SU_SO2emvn, SU_SO2emve, &
                                           rc, maskString, gridMask, &
                                           aviation_layers,   &
                                           aviation_lto_src, &
                                           aviation_cds_src, &
                                           aviation_crs_src)

! !USES:

  implicit NONE

! !INPUT PARAMETERS:

   integer, intent(in)              :: i1, i2, j1, j2, km, nbins, nvolc
   integer, intent(in)              :: nymd, nhms
   real, intent(in)                 :: cdt, fSO4ant, eAircraftFuel
   real, pointer, dimension(:,:)    :: so2anthro_l1_src, so2anthro_l2_src, &
                                       so2biomass_src, dmso_conc, &
                                       so2ship_src, so4ship_src
   real, pointer, dimension(:,:,:)  :: aircraft_fuel_src
   real, pointer, dimension(:,:)    :: cell_area
   type(ESMF_Grid), intent(inout)   :: Grid  ! ESMF Grid
   real, pointer, dimension(:)      :: vlat, vlon, velev, vcloud, vso2
   integer, pointer, dimension(:)   :: vstart, vend
   real, pointer, dimension(:,:)    :: oro, u10m, v10m, pblh, hsurf
   real, pointer, dimension(:,:,:)  :: tmpu, rhoa, hghte, delp

! !OUTPUT PARAMETERS:

   real, pointer, dimension(:,:,:)  :: dms, so2, so4
   type(Chem_Array), intent(inout)  :: SU_emis(nbins)  ! SU emissions, kg/m2/s
   type(Chem_Array), intent(inout)  :: SU_SO4eman  ! SO4 anthro emissions, kg/m2/s
   type(Chem_Array), intent(inout)  :: SU_SO2eman  ! SO2 anthro emissions, kg/m2/s
   type(Chem_Array), intent(inout)  :: SU_SO2embb  ! SO2 bioburn emissions, kg/m2/s
   type(Chem_Array), intent(inout)  :: SU_SO2emvn  ! SO2 volcanic (non-explosive) emissions, kg/m2/s
   type(Chem_Array), intent(inout)  :: SU_SO2emve  ! SO2 volcanic (explosive) emissions, kg/m2/s
   integer, intent(out)             :: rc    ! Error return code:
                                             !  0 - all is well
                                             !  1 - 

!  Optional parameters to be passed for masking separate instances
   character(len=*), OPTIONAL, intent(in) :: maskString             !Delimited string of integers
   real, OPTIONAL, intent(in)             :: gridMask(i1:i2,j1:j2)  !Grid mask (NOTE: No ghosting) 

   real, optional, intent(in)                         :: aviation_layers(4)   ! heights of LTO, CDS and CRS layers
   real, optional, dimension(i1:i2,j1:j2), intent(in) :: aviation_lto_src     ! SO2 Aviation-LTO
   real, optional, dimension(i1:i2,j1:j2), intent(in) :: aviation_cds_src     ! SO2 Aviation-CDS
   real, optional, dimension(i1:i2,j1:j2), intent(in) :: aviation_crs_src     ! SO2 Aviation-CRS

   character(len=*), parameter :: myname = 'SU_Emission'

! !DESCRIPTION: Updates the SU concentration with emissions every timestep
!
! !REVISION HISTORY:
!
!  06Nov2003, Colarco
!  Based on Ginoux
!
!EOP
!-------------------------------------------------------------------------

! !Local Variables
   integer  ::  i, j, k, n, ios
   real :: p1, z1, dz, deltaz, deltap, f100, f500, fPblh
   real :: sCO2, schmidt, w10m, akw, sst
   real :: zpbl

                                 ! pressure at 100m, 500m, & PBLH
   real, dimension(i1:i2,j1:j2) :: p100, p500, pPblh  
   real, dimension(i1:i2,j1:j2) :: p0, z0, ps


   real, dimension(i1:i2,j1:j2) :: srcSO2
   real, dimension(i1:i2,j1:j2) :: srcSO4
   real, dimension(i1:i2,j1:j2) :: srcDMS  
   real, dimension(i1:i2,j1:j2) :: srcSO4anthro
   real, dimension(i1:i2,j1:j2) :: srcSO2anthro
   real, dimension(i1:i2,j1:j2) :: srcSO2bioburn
   real, dimension(i1:i2,j1:j2) :: srcSO2volc
   real, dimension(i1:i2,j1:j2) :: srcSO2volce
   real, dimension(i1:i2,j1:j2) :: so2srcvolc 

   integer :: it
   real :: hup, hlow, dzvolc
   real :: deltaSO2v, so2volcano
   integer :: ijl, ijkl

real :: deltaSO2v_sum

!  Handle masking of volcanic sources
    logical :: doingMasking
    INTEGER, ALLOCATABLE :: regionNumbers(:),flag(:)
    INTEGER, ALLOCATABLE :: mask(:,:)

! Indices for volcanic sources
    integer :: iVolc(nvolc), jVolc(nvolc)

!  Aviation
   real, dimension(i1:i2,j1:j2,km) :: emis_aviation
   real, dimension(i1:i2,j1:j2,km) :: srcAviation
   real                            :: z_lto_bot, z_lto_top
   real                            :: z_cds_bot, z_cds_top
   real                            :: z_crs_bot, z_crs_top
#ifdef DEBUG
   real :: qmin, qmax
#endif

!class (logger), pointer :: lgr
!lgr => logging%get_logger('volcanic_emissions')

   _UNUSED_DUMMY(nymd)

   ijl  = ( i2 - i1 + 1 ) * ( j2 - j1 + 1 )
   ijkl = ijl * km

!  Initialize local variables
!  --------------------------
   srcSO2 = 0.0
   srcSO4 = 0.0
   srcDMS = 0.0
   srcSO2volc = 0.0
   srcSO2volce = 0.0
   so2srcvolc  = 0.0

   do n = 1, nbins
    if( associated(SU_emis(n)%data2d) ) SU_emis(n)%data2d(i1:i2,j1:j2) = 0.0
   end do
   if( associated(SU_SO4eman%data2d)) SU_SO4eman%data2d(i1:i2,j1:j2) = 0.0
   if( associated(SU_SO2eman%data2d)) SU_SO2eman%data2d(i1:i2,j1:j2) = 0.0
   if( associated(SU_SO2embb%data2d)) SU_SO2embb%data2d(i1:i2,j1:j2) = 0.0
   if( associated(SU_SO2emvn%data2d)) SU_SO2emvn%data2d(i1:i2,j1:j2) = 0.0
   if( associated(SU_SO2emve%data2d)) SU_SO2emve%data2d(i1:i2,j1:j2) = 0.0

!  Distribute aircraft emissions from LTO, CDS and CRS layers
!  ----------------------------------------------------------
   z_lto_bot = max(1e-3, aviation_layers(1))
   z_lto_top = max(2e-3, aviation_layers(2))

   z_cds_bot = max(2e-3, aviation_layers(2))
   z_cds_top = max(3e-3, aviation_layers(3))

   z_crs_bot = max(3e-3, aviation_layers(3))
   z_crs_top = max(4e-3, aviation_layers(4))

   emis_aviation = 0.0
   srcAviation   = 0.0

   call distribute_aviation_emissions(delp, rhoa, z_lto_bot, z_lto_top, aviation_lto_src, emis_aviation, i1, i2, j1, j2, km)
   srcAviation = srcAviation + emis_aviation

   call distribute_aviation_emissions(delp, rhoa, z_cds_bot, z_cds_top, aviation_cds_src, emis_aviation, i1, i2, j1, j2, km)
   srcAviation = srcAviation + emis_aviation

   call distribute_aviation_emissions(delp, rhoa, z_crs_bot, z_crs_top, aviation_crs_src, emis_aviation, i1, i2, j1, j2, km)
   srcAviation = srcAviation + emis_aviation

!  Find the pressure of the 100m, 500m, and PBLH altitudes
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

    if( associated(SU_emis(nSO2)%data2d) ) &
                   SU_emis(nSO2)%data2d =  SU_emis(nSO2)%data2d + srcSO2
    if( associated(SU_emis(nSO4)%data2d) ) &
                   SU_emis(nSO4)%data2d =  SU_emis(nSO4)%data2d + srcSO4
    if( associated(SU_SO4eman%data2d) ) &
                   SU_SO4eman%data2d    =  SU_SO4eman%data2d    + srcSO4anthro
    if( associated(SU_SO2eman%data2d) ) &
                   SU_SO2eman%data2d    =  SU_SO2eman%data2d    + srcSO2anthro
    if( associated(SU_SO2embb%data2d) ) &
                   SU_SO2embb%data2d    =  SU_SO2embb%data2d    + srcSO2bioburn

#ifdef DEBUG
   if ( k >= km-1 ) then
      call pmaxmin('SU: srcSO2        ', srcSO2 , qmin, qmax, ijl, 1, 1. )
      call pmaxmin('SU: srcSO4        ', srcSO4 , qmin, qmax, ijl, 1, 1. )
      call pmaxmin('SU: srcSO4anthro  ', srcSO4anthro , qmin, qmax, ijl, 1, 1. )
      call pmaxmin('SU: srcSO2anthro  ', srcSO2anthro , qmin, qmax, ijl, 1, 1. )
      call pmaxmin('SU: srcSO2bioburn ', srcSO2bioburn , qmin, qmax, ijl, 1, 1. )
   end if
#endif

   end do ! k

!  Create mask for volcanic emissions
!  When masking, both the mask and the string
!  of integers (region numbers) must be present
!  --------------------------------------------
   IF ( (PRESENT(maskString) .AND. .NOT. PRESENT(gridMask)) .OR.   &
        (PRESENT(gridMask) .AND. .NOT. PRESENT(maskString))      ) &
      CALL die ( myname, ": Both gridMask and maskString must be specified." )

   IF(PRESENT(gridMask)) THEN
      IF(TRIM(maskString) == "-1") THEN
       doingMasking = .FALSE.
      ELSE
       doingMasking = .TRUE.
      END IF
   ELSE
    doingMasking = .FALSE.
   END IF

!  Masking initialization
!  ----------------------
   SetMask: IF(doingMasking) THEN

    k = 32
    ALLOCATE(regionNumbers(k),flag(k),mask(i1:i2,j1:j2),STAT=ios)
    IF ( ios /= 0 ) CALL die ( myname, ": Cannot allocate for masking.")

!  Obtain region numbers from delimited list of integers
!  -----------------------------------------------------
    regionNumbers(:) = 0
    CALL Chem_UtilExtractIntegers(maskString,k,regionNumbers,RC=ios)
    IF ( ios /= 0 ) CALL die ( myname, ": Unable to extract integers for regionNumbers.")

!  How many integers were found?
!  -----------------------------
    flag(:) = 1
    WHERE(regionNumbers(:) == 0) flag(:) = 0
    k = SUM(flag)
    DEALLOCATE(flag,STAT=ios)
    IF ( ios /= 0 ) CALL die ( myname, ": Cannot dallocate flag.")

!  Set local mask to 1 where gridMask matches each integer (within precision!) 
!  ---------------------------------------------------------------------------
    mask(i1:i2,j1:j2) = 0
    DO ios=1,k
     WHERE(regionNumbers(ios)-0.01 <= gridMask .AND. &
           gridMask <= regionNumbers(ios)+0.01) mask = 1
    END DO

   END IF SetMask

!  Add the volcanic source
!  -----------------------
!  Note: the model lat and lon are wired in radians
!  but the dx and dy are in degrees
!  We add here the regional mask checking to exclude
!  any points that are not wanted for a particular
!  ensemble member calling this routine.

!  Point source volcanos (loop over each volcano)
   srcSO2volc(:,:) = 0.
   srcSO2volce(:,:) = 0.
   z0 = hghte(:,:,km)

   if(nvolc > 0) then

 !    Get indices for volcanic emissions
 !    ----------------------------------
      call MAPL_GetHorzIJIndex(nvolc,iVolc,jVolc,Grid=Grid,lon=vLon/radToDeg,lat=vLat/radToDeg,rc=rc)

      if ( rc /= 0 ) call die(myname,'cannot get indices for volcanic emissions')

!     Loop over all volcanoes in the database
      do it = 1, nvolc

deltaSO2v_sum = 0.

         i = iVolc(it)
         j = jVolc(it)
         
!        Skip this volcano?
!        ------------------
         if ( i<1 .OR. j<1 ) cycle ! volcano not in sub-domain
!         if(doingMasking) then
!            if( mask(i,j) == 0 ) cycle
!         end if

!        Check time against time range of eruption
!        -----------------------------------------
         if(nhms < vStart(it) .or. nhms >= vEnd(it)) cycle

         so2volcano = 0.

!        Emissions per volcano
!        -------------------------------------------------------------------------------
         if(cell_area(i,j) .gt. 1.) then
            so2volcano = vSO2(it) /cell_area(i,j)     ! to kg SO2/sec/m2
            so2volcano = max(so2volcano,tiny(so2volcano))
         endif

!        Distribute in the vertical
!        Database provides altitude of top of volcano cone (vElev) and altitude
!        of plume top (vCloud).  If vCloud != vElev then distribute emissions
!        in top 1/3 of column extending from vElev to vCloud (case of explosive
!        eruption), else put emissions in grid cell containing vElev (degassing)
!        --------------------------
         hup  = vCloud(it)
         hlow = vElev(it)
         if (hup .ne. hlow) then
            hlow = hup - (hup-hlow)/3.
         endif

!        Diagnostic - sum of volcanos
!        ----------------------------
         if (hup .eq. hlow) then
            srcSO2volc(i,j) = srcSO2volc(i,j) + so2volcano
         else
            srcSO2volce(i,j) = srcSO2volce(i,j) + so2volcano
         endif

         dzvolc = hup-hlow
         do k = km, 1, -1
            z1 = hghte(i,j,k-1) ! geopotential altitude at gridbox top
            dz = z1-z0(i,j)     ! thickness of gridbox
            deltaSO2v = 0.

!           Volcano is above this level
!           ---------------------------
            if(z1 .lt. hlow) then
               z0(i,j) = z1
               cycle
            end if

!           Volcano is below this level (except at surface)
!           -----------------------------------------------
            if(z0(i,j) .gt. hup .and. k .ne. km) then
               z0(i,j) = z1
               cycle
            end if

!           Volcano is in this level
!           ------------------------
            if( (k .eq. km .and. z0(i,j) .gt. hup) .or. &     ! below surface
                 (z0(i,j) .le. hlow .and. z1 .ge. hup) ) then ! in level
               deltaSO2v = so2volcano

!           Volcano only partly in level                       ! Cell:
!           ----------------------------
            else if (z0(i,j) .lt. hlow .and. z1 .lt. hup) then ! has bottom of cloud
               deltaSO2v = (z1-hlow)/dzvolc*so2volcano
 
            else if (z0(i,j) .gt. hlow .and. z1 .gt. hup) then ! has top of cloud
               deltaSO2v = (hup-z0(i,j))/dzvolc*so2volcano
 
            else                                               ! is filled with cloud
               deltaSO2v = dz/dzvolc*so2volcano
            end if

            z0(i,j) = z1
            so2(i,j,k) = so2(i,j,k) + deltaSO2v*cdt*grav/delp(i,j,k)

deltaSO2v_sum = deltaSO2v_sum + deltaSO2v

      end do ! k
!call lgr%debug('legacy emissions at %g0 %g0 : %g25.17', vLat(it), vLon(it), deltaSO2v_sum)

   enddo     ! it

   endif


!if(mapl_am_i_root()) print*,'SU deltaSO2v_sum = ',deltaSO2v_sum

!  Diagnostics -- this is really the point defined volcanos
   if(associated(SU_SO2emve%data2d)) then
      SU_SO2emve%data2d = srcSO2volce
   endif
   if(associated(SU_SO2emvn%data2d)) then
      SU_SO2emvn%data2d = srcSO2volc
if(mapl_am_i_root()) print*,'SU inside sum(SO2EMVN) = ',sum(SU_SO2emvn%data2d)
   endif
   if( associated(SU_emis(nSO2)%data2d) ) &
                  SU_emis(nSO2)%data2d =  SU_emis(nSO2)%data2d + srcSO2volc + srcSO2volce

!  Clean up volcano masking function
   IF(doingMasking) THEN
    DEALLOCATE(regionNumbers, mask, STAT=ios)
    IF ( ios /= 0 ) CALL die ( myname, ': Cannot deallocate masking tape.')
   END IF

#ifdef DEBUG
      call pmaxmin('SU: srcSO2volcExp ', srcSO2volc , qmin, qmax, ijl, 1, 1. )
#endif

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
   k = km
   sCO2 = 600.
   do j = j1, j2
    do i = i1, i2 
     sst = tmpu(i,j,k)-273.15
     if(sst .gt. 28.) sst = 28.
!    only valid for ocean and warm enough temperatures
     if( (oro(i,j) /= OCEAN) .or. (sst .lt. -20.)) cycle
     schmidt = 2764.0 - 147.12*sst + 3.726*(sst**2.) - 0.038*(sst**3.)
!    w10m is the 10-m wind speed in m s-1
     w10m = sqrt(u10m(i,j)**2. + v10m(i,j)**2.)
     if(w10m .le. 3.6) then
      akw = 0.17*w10m*((sCO2/schmidt)**0.667)
     else if (w10m .le. 13.) then
      akw = (2.85*w10m - 9.65)*sqrt(sCO2/schmidt)
     else
      akw = (5.90*w10m - 49.3)*sqrt(sCO2/schmidt)
     endif
!    This parameterization has put akw in units cm hr-1 -> goto m s-1
     akw = akw/100./3600.
!    DMSo concentration is nMol/L
!    Want to put the source into units of kg m-2 s-1
     srcDMS(i,j) = akw * (fmassDMS/1000.)*(dmso_conc(i,j)*1.e-9/1.e-3)
     dms(i,j,k) =  dms(i,j,k) + srcDMS(i,j)*cdt*grav/delp(i,j,k)
     end do
   end do

   if( associated(SU_emis(nDMS)%data2d) ) SU_emis(nDMS)%data2d =  srcDMS

#ifdef DEBUG
   call pmaxmin('SU: srcDMS        ', srcDMS , qmin, qmax, ijl, 1, 1. )
#endif

   rc = 0

contains
   subroutine distribute_aviation_emissions(delp, rhoa, z_bot, z_top, emissions_layer, emissions, i1, i2, j1, j2, km)

    implicit none

    integer, intent(in) :: i1, i2, j1, j2, km

    real, dimension(:,:,:), intent(in) :: delp
    real, dimension(:,:,:), intent(in) :: rhoa
    real, dimension(:,:),   intent(in) :: emissions_layer
    real, intent(in)                   :: z_bot
    real, intent(in)                   :: z_top
    real, dimension(:,:,:), intent(out):: emissions
    
!   local
    integer :: i, j, k
    integer :: k_bot, k_top
    real    :: z_
    real, dimension(km) :: z, dz, w_
    
    do j = j1, j2
        do i = i1, i2
            ! find level height
            z = 0.0
            z_= 0.0 

            do k = km, 1, -1
                dz(k) = delp(i,j,k)/rhoa(i,j,k)/grav
                z_    = z_ + dz(k)
                z(k)  = z_
            end do

            ! find the bottom level
            do k = km, 1, -1
                if (z(k) >= z_bot) then
                    k_bot = k
                    exit
                end if
            end do
            
            ! find the top level
            do k = k_bot, 1, -1
                if (z(k) >= z_top) then
                    k_top = k
                    exit
                end if
            end do

            ! find the weights
            w_ = 0

!           if (k_top > k_bot) then
!               need to bail - something went wrong here
!           end if

            if (k_bot .eq. k_top) then
                w_(k_bot) = z_top - z_bot
            else
                do k = k_bot, k_top, -1
                    if ((k < k_bot) .and. (k > k_top)) then
                        w_(k) = dz(k)
                    else
                        if (k == k_bot) then
                            w_(k) = (z(k) - z_bot)
                        end if

                        if (k == k_top) then
                            w_(k) = z_top - (z(k)-dz(k))
                        end if
                    end if
                end do
            end if
           
            ! distribute emissions in the vertical 
            emissions(i,j,:) = (w_ / sum(w_)) * emissions_layer(i,j)
        end do 
    end do

    end subroutine distribute_aviation_emissions

   end subroutine SulfateDistributeEmissions


!-------------------------------------------------------------------------
!     NASA/GSFC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  GetVolcDailyTables - Get pointwise SO2 and altitude of volcanoes
!                                  from a daily file data base

!
! !INTERFACE:
!

   subroutine GetVolcDailyTables( nymd, volcnon_srcfilen, &
                                  nVolcPts, vLat, vLon, vElev, vCloud, vSO2, vStart, vEnd )
			  
! !USES:

  implicit NONE

! Data for volcanic emissions comes from the daily inventory of all
! volcanos (as represented by the text tables).  We return all the
! volcanic emissions (as points, per volcano).

  integer, intent(in)            :: nymd
  character(len=255)             :: volcnon_srcfilen
  integer                        :: nVolcPts
  real, pointer, dimension(:)    :: vLat, vLon, vElev, vCloud, vSO2
  integer, pointer, dimension(:) :: vStart, vEnd
  integer :: i, j, nLines, nCols, rc, STATUS, nymd1, nhms1, ios
  character(len=255) :: fname
  type(ESMF_Config)  :: cf
  real, pointer, dimension(:) :: vData

! If previous instance of volcano point data tables exist, deallocate it
! to get the correct number of elements
  if(associated(vLat))    deallocate(vLat, stat=ios)
  if(associated(vLon))    deallocate(vLon, stat=ios)
  if(associated(vSO2))    deallocate(vSO2, stat=ios)
  if(associated(vElev))   deallocate(vElev, stat=ios)
  if(associated(vCloud))  deallocate(vCloud, stat=ios)
  if(associated(vStart))  deallocate(vStart, stat=ios)
  if(associated(vEnd))    deallocate(vEnd, stat=ios)

! Daily files (e.g., from AEROCOM)
! --------------------------------
! Note: Volcanic emissions in these files are in mass of Sulfur
  nymd1 = nymd
  nhms1 = 120000
  call StrTemplate ( fname, volcnon_srcfilen, xid='unknown', &
                     nymd=nymd1, nhms=nhms1 )
  cf = ESMF_ConfigCreate()
  call ESMF_ConfigLoadFile(cf, fileName=trim(fname), rc=STATUS )
  call ESMF_ConfigGetDim(cf, nLines, nCols, LABEL='volcano::', rc=STATUS )
  nVolcPts = nLines
  allocate(vData(nCols), vLat(nLines), vLon(nLines), &
           vSO2(nLines), vElev(nLines), vStart(nLines), &
           vEnd(nLines), vCloud(nLines), stat=ios)
  vStart = -1
  vEnd   = -1
  call ESMF_ConfigFindLabel(cf, 'volcano::',rc=STATUS)
     do i = 1, nLines
      call ESMF_ConfigNextLine(cf, rc=rc)
      do j = 1, nCols
       call ESMF_ConfigGetAttribute(cf, vData(j), default=-1.)
      end do
      vLat(i)    = vData(1)
      vLon(i)    = vData(2)
      vSO2(i)    = vData(3) * fMassSO2 / fMassSulfur
      vElev(i)   = vData(4)
      vCloud(i)  = vData(5)
      if(nCols >= 6) vStart(i)  = vData(6)
      if(nCols >= 7) vEnd(i)    = vData(7)
  end do

  call ESMF_ConfigDestroy(cf)
  deallocate(vData, stat=ios)

  end subroutine GetVolcDailyTables

!-------------------------------------------------------------------------
!     NASA/GSFC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  GetVolcExplosive - Get pointwise SO2 and altitude of explosive
!                                volcanos from a data base

!
! !INTERFACE:
!

   subroutine GetVolcExplosive( nymd, nVolcExp, &
                                vLatP, vLonP, vElevP, vCloudP, vSO2P )
			  
! !USES:

  implicit NONE

! Description
! Data for volcanic explosions provided by Thomas Diehl.  I have converted
! from kt SO2 event-1 to kt SO2 day-1 in the data table.  Following that I
! convert to kg SO2 s-1 needed in emissions.  Assumption is uniform 
! injection over 24 hour period.
! What is returned is the number of volcanoes and an array of locations,
! elevations, and SO2 amounts.

  integer                       :: nymd
  real, pointer, dimension(:)   :: vLatP, vLonP, vElevP, vCloudP, vSO2P
  integer         :: it, nv, ios, nVolcExp

! database parameters and data
  integer, parameter :: nvolc = 349
  integer :: startday(nvolc), endday(nvolc)
  real    :: so2exp(nvolc), vlon(nvolc), vlat(nvolc), & 
             velev(nvolc), celev(nvolc)

! If previous instance of volcano point data tables exist, deallocate it
! to get the correct number of elements
  if(associated(vLatP))    deallocate(vLatP, stat=ios)
  if(associated(vLonP))    deallocate(vLonP, stat=ios)
  if(associated(vSO2P))    deallocate(vSO2P, stat=ios)
  if(associated(vElevP))   deallocate(vElevP, stat=ios)
  if(associated(vCloudP))  deallocate(vCloudP, stat=ios)

  data startday / &
    20051231, 20051231, 20051231, 20061231, 20001015, 20061231,   &
    20000301, 20021019, 20010717, 20061231, 20010905, 20001220,   &
    20050619, 20061231, 20020615, 19990815, 20000425, 20010525,   &
    19990815, 19991109, 19990719, 20000315, 20001214, 19990915,   &
    19990224, 19990815, 19990712, 19990225, 19990418, 19990915,   &
    20031028, 20000615, 19990329, 19990402, 19990701, 19990417,   &
    19990412, 19990915, 19990421, 19990527, 19990915, 20000615,   &
    19990514, 20001209, 20000415, 19990915, 19990630, 20000319,   &
    19990628, 19991023, 19990720, 20010729, 19990807, 20061231,   &
    19991115, 19991020, 19991116, 20010805, 20000302, 19991229,   &
    20000315, 20000127, 20000210, 20000208, 20000214, 20000304,   &
    20000224, 20000226, 20000229, 20000306, 20000905, 20000403,   &
    20000326, 20000518, 20010915, 20001029, 20010818, 20030831,   &
    20010215, 20000610, 20001030, 20000604, 20001113, 20000818,   &
    20001018, 20001015, 20000831, 20001104, 20000721, 20000922,   &
    20010705, 20010115, 20000820, 20000823, 20000827, 20000904,   &
    20000910, 20001108, 20000926, 20061231, 20000930, 20001101,   &
    20010721, 20010416, 20001129, 20001130, 20010115, 20001215,   &
    20001218, 20001222, 20040705, 20010428, 20010808, 20021124,   &
    20010429, 20010218, 20010302, 20010219, 20010219, 20010415,   &
    20010405, 20010404, 20010605, 20010425, 20010429, 20031122,   &
    20010501, 20010503, 20011209, 20020827, 20011003, 20010619,   &
    20010707, 20010809, 20010915, 20010809, 20010730, 20030712,   &
    20010806, 20061231, 20010828, 20011115, 20011019, 20011115,   &
    20011005, 20061231, 20011031, 20061231, 20011126, 20020106,   &
    20021026, 20061231, 20020116, 20020521, 20020116, 20020117,   &
    20020203, 20020422, 20020515, 20020315, 20020609, 20040415,   &
    20020715, 20021006, 20030409, 20021216, 20020617, 20020607,   &
    20020825, 20021115, 20020725, 20020802, 20020815, 20020802,   &
    20020915, 20020820, 20021103, 20021015, 20020925, 20020926,   &
    20020929, 20040217, 20021106, 20021011, 20021207, 20021012,   &
    20021026, 20021027, 20021030, 20021103, 20021103, 20021120,   &
    20030128, 20021202, 20021112, 20021114, 20021219, 20021116,   &
    20021118, 20021203, 20030110, 20021128, 20021228, 20030401,   &
    20030101, 20040408, 20030228, 20030418, 20031015, 20030528,   &
    20030723, 20031109, 20030514, 20030416, 20031010, 20030417,   &
    20030703, 20030506, 20030510, 20030513, 20030523, 20030523,   &
    20040325, 20030602, 20061231, 20040110, 20030901, 20030712,   &
    20030608, 20030609, 20040614, 20030616, 20030723, 20030714,   &
    20030715, 20031008, 20030801, 20031002, 20030915, 20030901,   &
    20030904, 20030912, 20031115, 20031011, 20040328, 20031209,   &
    20040114, 20050215, 20040127, 20040205, 20040214, 20040224,   &
    20040915, 20040502, 20040925, 20050805, 20040414, 20050405,   &
    20041008, 20041003, 20040528, 20050222, 20040517, 20040526,   &
    20040607, 20040815, 20040912, 20040608, 20040609, 20040624,   &
    20041024, 20040624, 20040916, 20040704, 20050207, 20050911,   &
    20040730, 20040805, 20061231, 20040914, 20050315, 20040915,   &
    20040915, 20041209, 20041005, 20061231, 20041212, 20061231,   &
    20041110, 20041104, 20061231, 20041111, 20041115, 20041123,   &
    20041125, 20041225, 20041126, 20041127, 20041219, 20041209,   &
    20041213, 20041227, 20041220, 20050127, 20050214, 20061231,   &
    20050407, 20061231, 20050525, 20050128, 20061231, 20050216,   &
    20050331, 20050227, 20050225, 20050223, 20050407, 20050701,   &
    20050406, 20050903, 20050718, 20050518, 20050414, 20061231,   &
    20050418, 20061231, 20050718, 20050504, 20050529, 20061231,   &
    20050616, 20051007, 20050815, 20051112, 20051104, 20050929,   &
    20051001, 20051005, 20061231, 20051117, 20051030, 20061231,   &
    20051208, 20061231, 20061231, 20051201, 20061231, 20051222,   &
    20061231  &
    /
  data endday / &
    20051231, 20051231, 20051231, 20061231, 20001015, 20061231,   &
    20000301, 20021019, 20010717, 20061231, 20010905, 20001220,   &
    20050619, 20061231, 20020615, 19990815, 20000425, 20010525,   &
    19990815, 19991109, 19990719, 20000315, 20001214, 19990915,   &
    19990224, 19990815, 19990712, 19990225, 19990418, 19990915,   &
    20031028, 20000615, 19990329, 19990402, 19990701, 19990417,   &
    19990412, 19990915, 19990421, 19990527, 19990915, 20000615,   &
    19990514, 20001209, 20000415, 19990915, 19990630, 20000319,   &
    19990628, 19991023, 19990720, 20010729, 19990807, 20061231,   &
    19991115, 19991020, 19991116, 20010805, 20000302, 19991229,   &
    20000315, 20000127, 20000210, 20000208, 20000214, 20000304,   &
    20000224, 20000226, 20000229, 20000306, 20000905, 20000403,   &
    20000326, 20000518, 20010915, 20001029, 20010818, 20030831,   &
    20010215, 20000610, 20001030, 20000604, 20001113, 20000818,   &
    20001018, 20001015, 20000831, 20001104, 20000721, 20000922,   &
    20010705, 20010115, 20000820, 20000823, 20000827, 20000904,   &
    20000910, 20001108, 20000926, 20061231, 20000930, 20001101,   &
    20010721, 20010416, 20001129, 20001130, 20010115, 20001215,   &
    20001218, 20001222, 20040705, 20010428, 20010808, 20021124,   &
    20010429, 20010218, 20010302, 20010219, 20010219, 20010415,   &
    20010405, 20010404, 20010605, 20010425, 20010429, 20031122,   &
    20010501, 20010503, 20011209, 20020827, 20011003, 20010619,   &
    20010707, 20010809, 20010915, 20010809, 20010730, 20030712,   &
    20010806, 20061231, 20010828, 20011115, 20011019, 20011115,   &
    20011005, 20061231, 20011031, 20061231, 20011126, 20020106,   &
    20021026, 20061231, 20020116, 20020521, 20020116, 20020117,   &
    20020203, 20020422, 20020515, 20020315, 20020609, 20040415,   &
    20020715, 20021006, 20030409, 20021216, 20020617, 20020607,   &
    20020825, 20021115, 20020725, 20020802, 20020815, 20020802,   &
    20020915, 20020820, 20021103, 20021015, 20020925, 20020926,   &
    20020929, 20040217, 20021106, 20021011, 20021207, 20021012,   &
    20021026, 20021027, 20021030, 20021103, 20021103, 20021120,   &
    20030128, 20021202, 20021112, 20021114, 20021219, 20021116,   &
    20021118, 20021203, 20030110, 20021128, 20021228, 20030401,   &
    20030101, 20040408, 20030228, 20030418, 20031015, 20030528,   &
    20030723, 20031109, 20030514, 20030416, 20031010, 20030417,   &
    20030703, 20030506, 20030510, 20030513, 20030523, 20030523,   &
    20040325, 20030602, 20061231, 20040110, 20030901, 20030712,   &
    20030608, 20030609, 20040614, 20030616, 20030723, 20030714,   &
    20030715, 20031008, 20030801, 20031002, 20030915, 20030901,   &
    20030904, 20030912, 20031115, 20031011, 20040328, 20031209,   &
    20040114, 20050215, 20040127, 20040205, 20040214, 20040224,   &
    20040915, 20040502, 20040925, 20050805, 20040414, 20050405,   &
    20041008, 20041003, 20040528, 20050222, 20040517, 20040526,   &
    20040607, 20040815, 20040912, 20040608, 20040609, 20040624,   &
    20041024, 20040624, 20040916, 20040704, 20050207, 20050911,   &
    20040730, 20040805, 20061231, 20040914, 20050315, 20040915,   &
    20040915, 20041209, 20041005, 20061231, 20041212, 20061231,   &
    20041110, 20041104, 20061231, 20041111, 20041115, 20041123,   &
    20041125, 20041225, 20041126, 20041127, 20041219, 20041209,   &
    20041213, 20041227, 20041220, 20050127, 20050214, 20061231,   &
    20050407, 20061231, 20050525, 20050128, 20061231, 20050216,   &
    20050331, 20050227, 20050225, 20050223, 20050407, 20050701,   &
    20050406, 20050903, 20050718, 20050518, 20050414, 20061231,   &
    20050418, 20061231, 20050718, 20050504, 20050529, 20061231,   &
    20050616, 20051007, 20050815, 20051112, 20051104, 20050929,   &
    20051001, 20051005, 20061231, 20051117, 20051030, 20061231,   &
    20051208, 20061231, 20061231, 20051201, 20061231, 20051222,   &
    20061231  &
    /
  data so2exp / &
       0.004,    0.004,    0.008,    0.001,    0.011,    0.022,   &
       0.031,    0.004,    0.044,    0.042,    0.008,    0.063,   &
       0.001,    0.035,    0.001,    0.036,    0.177,    0.112,   &
       0.047,    0.012,    0.038,    0.001,    0.041,    0.062,   &
       0.044,    0.089,    0.108,   17.000,    1.513,    0.092,   &
       0.068,    0.005,   30.952,   30.952,    0.183,   30.952,   &
       1.700,    0.110,   21.000,    1.513,    0.016,    0.041,   &
     190.000,    0.030,    0.051,    0.002,    0.141,    0.423,   &
       2.250,    1.959,   16.000,    0.038,    5.667,    0.223,   &
       0.006,    2.250,    3.000,    0.006,    0.022,    0.750,   &
       0.279,   43.333,   43.333,    2.833,    9.167,    9.500,   &
      17.000,  250.000,  250.000,  250.000,    0.628,    0.708,   &
       1.308,    0.038,    0.032,    0.086,    0.036,    0.095,   &
       2.434,   46.429,    0.015,    8.500,    1.319,    1.250,   &
       0.155,    0.024,    0.362,    0.155,    8.500,    0.298,   &
       0.007,    0.002,   11.500,    1.250,    0.375,   20.583,   &
       0.225,    0.034,    1.250,    0.001,   17.143,   17.143,   &
       0.008,    0.014,    0.354,    9.000,    0.354,   23.000,   &
       0.041,   10.000,    0.089,    0.041,    0.540,    0.025,   &
       1.065,    8.219,   28.000,    8.219,   17.000,    8.219,   &
      11.017,   21.111,    0.315,    0.750,    4.000,    0.041,   &
      15.000,    1.065,    0.011,    0.036,    0.002,    9.583,   &
       7.037,    0.708,    0.039,    6.389,   33.000,    0.038,   &
       3.000,    0.006,   17.000,    0.004,    0.078,    0.043,   &
       2.250,    0.001,    2.250,    0.061,    2.250,    0.607,   &
       0.007,    0.009,   15.833,    0.891,    2.250,   30.000,   &
      10.556,    0.193,    0.177,    2.250,    0.274,    0.003,   &
       0.258,    0.385,    0.053,    0.011,    0.118,    2.250,   &
       0.236,    0.002,   12.264,   90.000,    0.125,   17.000,   &
      12.264,    1.889,    0.230,    0.102,  120.000,  120.000,   &
     120.000,    0.034,    8.387,    2.250,    0.039,    2.250,   &
       1.211,    8.500,   10.000,   10.000,   10.000,   15.278,   &
       1.211,    0.385,    0.436,    5.000,    0.436,   36.111,   &
      36.111,   36.111,    8.696,    2.250,    5.000,    0.170,   &
       1.700,    0.036,    0.385,    0.031,    0.070,   11.236,   &
       0.016,    0.009,    0.288,    2.125,    0.094,    2.250,   &
       0.221,    0.250,    1.797,   38.000,    1.797,    0.321,   &
       0.007,   36.000,    0.385,    0.841,    0.179,    1.797,   &
      12.778,   12.778,    0.108,   12.778,    0.061,    3.400,   &
      33.333,    0.038,   16.429,    0.266,    0.125,   16.000,   &
       0.095,  115.000,    0.005,    0.062,    0.015,    2.250,   &
       2.250,    0.288,    2.125,    0.281,    0.750,    1.889,   &
       0.080,    1.885,    0.083,    0.035,    5.667,    0.225,   &
       0.096,    0.258,   52.581,    0.001,    1.125,   17.000,   &
      52.581,    0.227,    0.022,    1.000,   30.000,    1.000,   &
       0.136,  190.000,    0.224,    2.250,    0.556,    0.005,   &
      20.000,   17.000,    0.003,    0.170,    0.012,   17.000,   &
     100.000,    0.170,    3.400,    0.021,    1.620,    0.021,   &
       0.751,  625.000,    0.022,    8.000,    0.450,    0.751,   &
      55.000,    0.531,    0.751,    7.000,    0.751,    0.225,   &
      15.000,    1.620,   40.000,    0.751,    0.405,    0.024,   &
       0.218,    0.024,    0.140,  140.000,    0.751,    0.895,   &
       0.279,    0.102,    4.444,    2.250,    0.083,    0.175,   &
      70.000,    0.225,    0.173,    0.061,    2.250,    0.027,   &
       5.667,    0.027,    0.187,  115.000,   38.235,    0.029,   &
       2.250,    0.177,    0.708,    0.157,    0.038,   28.750,   &
     115.000,    0.062,    0.088,    0.515,  277.778,    0.039,   &
       7.667,    0.042,    1.625,    5.667,    0.006,    0.321,   &
       0.006  &
    /
  data velev / &
        1185,     5230,     3676,     3794,     1330,     1222,   &
        2552,     2968,     3350,     2960,      688,     1536,   &
        1334,     3850,     2847,      704,     1413,     4784,   &
         321,     1807,      915,     1023,     5426,     1325,   &
         799,      813,     4835,     2882,     2857,     3800,   &
        1784,     1717,     4095,     4095,     1703,     4095,   &
        3283,     2891,     2857,     2857,     3428,     1745,   &
        4317,     3763,     1061,     1018,      799,     2462,   &
        2799,     2631,      915,      915,      728,     3283,   &
        5023,     2334,     5023,     5023,      635,     1700,   &
         704,     3058,     3058,     4835,     3058,     2631,   &
         799,     1491,     1491,     1491,      321,     2891,   &
        2882,     4276,      737,     5967,     1580,     1784,   &
        2745,     4095,      813,     1807,     2631,      815,   &
        2997,     3428,     2462,     2882,     5592,     4835,   &
        2552,      990,      815,      815,     1952,      815,   &
        2799,     1131,      815,     1750,     2334,     2334,   &
         704,      851,     2329,     2329,     2329,     5426,   &
        5426,     5426,      799,     5426,     2462,      815,   &
        2334,     1730,     3058,     1730,      321,     1730,   &
        3058,     2631,     2891,      635,     5426,     5426,   &
        2334,     2334,     1745,     3800,     1325,     1413,   &
        2631,     3350,      813,     2882,      915,      915,   &
        5023,     5023,     2334,      990,      161,     2597,   &
        2741,     1370,     2552,     1536,     4784,     2882,   &
        3350,     3763,     2631,     1807,     2130,     3470,   &
        3470,     1816,     1580,      833,     4835,     4784,   &
         704,     3470,     1330,     1745,     2552,     4276,   &
        3332,      990,     3058,     3058,     2799,      140,   &
        3058,      394,     2334,     2507,      725,      725,   &
         725,      688,     3470,     2462,     4784,     1703,   &
        3350,     5592,     3350,     3350,     3350,     3562,   &
        3350,     3470,     2665,     2665,     2665,     2631,   &
        2631,     2631,     3562,     2435,     3470,     1580,   &
        2882,     4835,     3470,     2568,      704,     3470,   &
        2435,     3350,     2462,     3125,     2334,     4784,   &
        1816,      990,      790,      790,      790,     1807,   &
        2847,      790,     3470,     2631,     1703,      790,   &
        1413,     1413,     2745,     1413,     1745,     1592,   &
         915,      915,     2882,     1715,     3212,     1784,   &
        1784,     1580,      635,     2462,     1807,     5592,   &
        1592,     2882,     1330,     1703,     3350,      833,   &
        2507,      915,      704,     1784,     2334,      790,   &
        3332,     2631,     3058,     1325,     2857,     5426,   &
        3058,     1320,     2462,     2329,     2329,     2329,   &
        3800,     1230,     1703,      635,     4276,     2552,   &
        2060,     2891,     2847,     2568,     3350,     1413,   &
        2568,     2568,     3726,     2549,     1784,      799,   &
        1807,     1725,     3562,     1807,     1061,     1807,   &
        1807,     1330,     1807,     1807,     1807,      815,   &
        1784,     1784,     1807,     1807,     2507,     5426,   &
        4835,      688,     2435,     1807,     1807,     1156,   &
        1413,     1703,     2631,     1533,     1816,     2334,   &
         790,      790,     2597,      815,     1592,      915,   &
        2361,     1330,     1784,     5592,     1476,      354,   &
        2381,     1730,     3332,     1700,     2507,     1442,   &
        2381,      990,     2631,      564,     1490,     1413,   &
        2361,     4276,     1496,     2882,     1252,     3350,   &
        1784  &
    /
  data celev / &
        9000,     9000,     9000,     6794,     9000,     1772,   &
        9000,     5968,     9000,     3510,     3688,     9000,   &
        1884,     9000,     3397,     3704,     9000,     9000,   &
        3321,     9000,     9000,     1073,     9000,     4325,   &
        1349,     3813,     7835,     5882,     9000,     6800,   &
        9000,     2267,     7095,     7095,     4703,     7095,   &
        6283,     5891,     9000,     9000,     3978,     4745,   &
        4867,     6763,     4061,     1068,     1349,     9000,   &
        3349,     3181,     9000,     9000,     3728,    18000,   &
        8023,     2884,     8023,     8023,     1185,     2250,   &
        3704,     6058,     6058,     7835,     6058,     3181,   &
        3799,     9000,     9000,     9000,     9000,     5891,   &
        5882,     4826,     3737,     8967,     4580,     9000,   &
        5745,     7095,     1363,     4807,     3181,     9000,   &
        5997,     3978,     5462,     5882,     8592,     7835,   &
        3102,     1040,     9000,     9000,     2502,     9000,   &
        3349,     1681,     9000,     2300,    18000,    18000,   &
        1254,     1401,     5329,     5329,     5329,     9000,   &
        9000,     9000,     9000,     9000,     9000,     3815,   &
        9000,    18000,     6058,    18000,     3321,    18000,   &
        6058,     3181,     5891,     1185,     9000,     9000,   &
        9000,     9000,     2295,     6800,     1375,     9000,   &
        3181,     6350,     1363,     9000,     9000,     9000,   &
        8023,     8023,     5334,     1040,      711,     3147,   &
        3291,     1920,     3102,     9000,     5334,     5882,   &
        3900,     6763,     3181,     9000,     2680,     4020,   &
        4020,     4816,     4580,     1383,     7835,     5334,   &
        3704,     6470,     4330,     2295,     3102,     4826,   &
        6332,     1040,     6058,     6058,     3349,     3140,   &
        6058,     3394,     5334,     3057,    18000,    18000,   &
       18000,     3688,     6470,     3012,     5334,     2253,   &
        9000,     8592,     3900,     9000,     3900,    18000,   &
        9000,     6470,     5665,     5665,     5665,     5631,   &
        5631,     5631,    18000,     2985,     6470,     4580,   &
        5882,     7835,     6470,     3118,     3704,     6470,   &
        2985,     3900,     5462,     6125,     5334,     5334,   &
        4816,     1040,     9000,     9000,     9000,     2357,   &
        3397,     9000,     6470,     3181,     4703,     9000,   &
        9000,     9000,     2795,     9000,     2295,     4592,   &
        9000,     9000,     9000,     4715,     3762,     9000,   &
        9000,     9000,      685,     2512,     2357,     6142,   &
        2142,     9000,     4330,     2253,     3900,     3833,   &
        5507,     9000,     3704,     4784,     5334,     9000,   &
        6332,     2681,     6058,     1375,     3407,     8426,   &
        6058,     4320,     3012,     5329,     5329,     5329,   &
        6800,     1780,     4703,     1185,     9000,     3102,   &
        2110,     5891,     3397,     5568,     3900,     4413,   &
        5568,     5568,     6726,     5549,     9000,     3799,   &
       18000,     9000,     6562,    18000,     1611,    18000,   &
       18000,     4330,    18000,    18000,    18000,     1365,   &
        9000,     9000,    18000,    18000,     5507,     8426,   &
        7835,     3688,     5435,    18000,    18000,     4156,   &
        4413,     2253,     2681,     2083,     2366,     5334,   &
        9000,     9000,     5597,     1365,     2142,     3915,   &
        5361,     4330,     4784,     9000,     4476,     3354,   &
        2931,     4730,     6332,     4700,     3057,     9000,   &
        9000,     1040,     2681,     3564,     9000,     4413,   &
        9000,     7276,     4496,     5882,     1802,     3900,   &
        2334  &
    /
  data vlon / &
     127.880,  281.659,  112.920,  167.170,  148.420,  204.708,   &
     269.399,  110.442,   15.004,   35.902,  152.203,  159.430,   &
     168.120,  256.380,  288.070,  130.308,  168.346,  281.402,   &
     177.180,  145.061,  297.820,  332.680,  261.378,  127.642,   &
     129.716,  105.423,  160.638,  160.587,  196.030,  101.264,   &
     125.400,  115.375,    9.170,    9.170,  122.775,    9.170,   &
     161.360,  100.473,  196.030,  196.030,  109.208,  272.996,   &
     215.980,  269.120,  273.155,  123.590,  129.716,  123.685,   &
     114.242,   55.713,  297.820,  297.820,  273.298,  161.360,   &
     281.558,  151.330,  281.558,  281.558,  273.839,  274.378,   &
     130.308,   29.200,   29.200,  160.638,   29.200,   55.713,   &
     129.716,  340.300,  340.300,  340.300,  177.180,  100.473,   &
     160.587,  282.630,  140.843,  288.150,  124.792,  124.725,   &
      73.513,    9.170,  105.423,  145.061,   55.713,  139.529,   &
     288.830,  109.208,  123.685,  160.587,  292.270,  160.638,   &
     269.399,  333.550,  139.529,  139.529,  102.620,  139.529,   &
     114.242,  140.681,  139.529,  155.195,  151.330,  151.330,   &
     130.308,  165.800,  112.950,  112.950,  112.950,  261.378,   &
     261.378,  261.378,  129.716,  261.378,  123.685,  139.529,   &
     151.330,  190.056,   29.200,  190.056,  177.180,  190.056,   &
      29.200,   55.713,  100.473,  273.839,  261.378,  261.378,   &
     151.330,  151.330,  272.996,  101.264,  127.642,  168.346,   &
      55.713,   15.004,  105.423,  160.587,  297.820,  297.820,   &
     281.558,  281.558,  151.330,  333.550,  141.290,  100.679,   &
     158.830,  333.670,  269.399,  159.430,  281.402,  160.587,   &
      15.004,  269.120,   55.713,  145.061,  271.731,   29.250,   &
      29.250,  155.458,  124.792,  168.370,  160.638,  281.402,   &
     130.308,   29.250,  148.420,  272.996,  269.399,  282.630,   &
     114.042,  333.550,   29.200,   29.200,  114.242,  148.121,   &
      29.200,  140.306,  151.330,  200.620,  125.425,  125.425,   &
     125.425,  152.203,   29.250,  123.685,  281.402,  122.775,   &
      15.004,  292.270,   15.004,   15.004,   15.004,  282.344,   &
      15.004,   29.250,  107.730,  107.730,  107.730,   55.713,   &
      55.713,   55.713,  282.344,  123.132,   29.250,  124.792,   &
     160.587,  160.638,   29.250,  138.526,  130.308,   29.250,   &
     123.132,   15.004,  123.685,  288.271,  151.330,  281.402,   &
     155.458,  333.550,  145.670,  145.670,  145.670,  145.061,   &
     288.070,  145.670,   29.250,   55.713,  122.775,  145.670,   &
     168.346,  168.346,   73.513,  168.346,  272.996,  131.106,   &
     297.820,  297.820,  160.587,  127.325,  288.623,  124.725,   &
     124.725,  124.792,  273.839,  123.685,  145.061,  292.270,   &
     131.106,  160.587,  148.420,  122.450,   15.004,  168.370,   &
     200.620,  297.820,  130.308,  125.400,  151.330,  145.670,   &
     114.042,   55.713,   29.200,  127.642,  196.030,  261.378,   &
      29.200,  125.500,  123.685,  112.950,  112.950,  112.950,   &
     101.264,   37.750,  122.450,  273.839,  282.630,  269.399,   &
     347.720,  100.473,  288.070,  138.526,   15.004,  168.346,   &
     138.526,  138.526,  116.470,  237.820,  124.725,  129.716,   &
     145.061,  342.670,  282.344,  145.061,  273.155,  145.061,   &
     145.061,  148.420,  145.061,  145.061,  145.061,  139.529,   &
     124.725,  124.725,  145.061,  145.061,  200.620,  261.378,   &
     160.638,  152.203,  123.132,  145.061,  145.061,  156.020,   &
     168.346,  122.450,   55.713,  185.846,  155.458,  151.330,   &
     145.670,  145.670,  100.679,  139.529,  131.106,  297.820,   &
      43.380,  148.420,  124.725,  292.270,  268.450,   93.858,   &
     270.370,  190.056,  114.042,  274.378,  200.620,   40.480,   &
     270.370,  333.550,   55.713,  150.030,  268.830,  168.346,   &
      43.380,  282.630,  167.830,  160.587,  206.570,   15.004,   &
     124.725  &
    /
  data vlat / &
       1.680,   -2.002,   -8.108,  -77.530,   -5.525,   19.425,   &
      14.381,   -7.542,   37.734,   -2.751,   -4.271,   54.050,   &
     -16.250,   19.514,  -39.420,   30.789,  -16.507,   -0.171,   &
     -37.520,   -4.100,   16.720,   38.730,   19.023,    1.475,   &
      29.635,   -6.102,   56.057,   55.978,   54.756,   -1.814,   &
       2.780,   -8.242,    4.203,    4.203,   -8.530,    4.203,   &
      56.653,   -0.381,   54.756,   54.756,   -7.242,   12.702,   &
      62.000,   14.473,   12.602,   -8.540,   29.635,   13.257,   &
      -8.058,  -21.229,   16.720,   16.720,   12.506,   56.653,   &
      -1.467,   -5.050,   -1.467,   -1.467,   11.984,   11.538,   &
      30.789,   -1.408,   -1.408,   56.057,   -1.408,  -21.229,   &
      29.635,   63.980,   63.980,   63.980,  -37.520,   -0.381,   &
      55.978,    1.220,   42.541,  -15.780,    1.358,    1.108,   &
     -53.106,    4.203,   -6.102,   -4.100,  -21.229,   34.079,   &
     -37.850,   -7.242,   13.257,   55.978,  -23.370,   56.057,   &
      14.381,  -57.780,   34.079,   34.079,   -3.520,   34.079,   &
      -8.058,   42.061,   34.079,   -6.140,   -5.050,   -5.050,   &
      30.789,  -10.380,   -7.942,   -7.942,   -7.942,   19.023,   &
      19.023,   19.023,   29.635,   19.023,   13.257,   34.079,   &
      -5.050,   52.825,   -1.408,   52.825,  -37.520,   52.825,   &
      -1.408,  -21.229,   -0.381,   11.984,   19.023,   19.023,   &
      -5.050,   -5.050,   12.702,   -1.814,    1.475,  -16.507,   &
     -21.229,   37.734,   -6.102,   55.978,   16.720,   16.720,   &
      -1.467,   -1.467,   -5.050,  -57.780,   24.754,   -0.978,   &
      53.255,  -58.420,   14.381,   54.050,   -0.171,   55.978,   &
      37.734,   14.473,  -21.229,   -4.100,   13.434,   -1.520,   &
      -1.520,   50.325,    1.358,  -16.680,   56.057,   -0.171,   &
      30.789,   -1.520,   -5.525,   12.702,   14.381,    1.220,   &
      -8.125,  -57.780,   -1.408,   -1.408,   -8.058,   -5.520,   &
      -1.408,   30.480,   -5.050,   56.170,    2.280,    2.280,   &
       2.280,   -4.271,   -1.520,   13.257,   -0.171,   -8.530,   &
      37.734,  -23.370,   37.734,   37.734,   37.734,   -0.077,   &
      37.734,   -1.520,   -7.320,   -7.320,   -7.320,  -21.229,   &
     -21.229,  -21.229,   -0.077,   10.412,   -1.520,    1.358,   &
      55.978,   56.057,   -1.520,   36.403,   30.789,   -1.520,   &
      10.412,   37.734,   13.257,  -38.692,   -5.050,   -0.171,   &
      50.325,  -57.780,   16.350,   16.350,   16.350,   -4.100,   &
     -39.420,   16.350,   -1.520,  -21.229,   -8.530,   16.350,   &
     -16.507,  -16.507,  -53.106,  -16.507,   12.702,   32.881,   &
      16.720,   16.720,   55.978,    0.800,  -36.863,    1.108,   &
       1.108,    1.358,   11.984,   13.257,   -4.100,  -23.370,   &
      32.881,   55.978,   -5.525,   -8.670,   37.734,  -16.680,   &
      56.170,   16.720,   30.789,    2.780,   -5.050,   16.350,   &
      -8.125,  -21.229,   -1.408,    1.475,   54.756,   19.023,   &
      -1.408,    3.670,   13.257,   -7.942,   -7.942,   -7.942,   &
      -1.814,  -46.900,   -8.670,   11.984,    1.220,   14.381,   &
     -37.092,   -0.381,  -39.420,   36.403,   37.734,  -16.507,   &
      36.403,   36.403,   -8.420,   46.200,    1.108,   29.635,   &
      -4.100,   64.420,   -0.077,   -4.100,   12.602,   -4.100,   &
      -4.100,   -5.525,   -4.100,   -4.100,   -4.100,   34.079,   &
       1.108,    1.108,   -4.100,   -4.100,   56.170,   19.023,   &
      56.057,   -4.271,   10.412,   -4.100,   -4.100,   50.680,   &
     -16.507,   -8.670,  -21.229,   52.381,   50.325,   -5.050,   &
      16.350,   16.350,   -0.978,   34.079,   32.881,   16.720,   &
     -11.750,   -5.525,    1.108,  -23.370,   -0.370,   12.278,   &
      13.853,   52.825,   -8.125,   11.538,   56.170,   12.600,   &
      13.853,  -57.780,  -21.229,   -5.450,   -0.830,  -16.507,   &
     -11.750,    1.220,  -15.400,   55.978,   59.363,   37.734,   &
       1.108  &
    /

!  Reorient the longitudes for GEOS-5
   where(vlon > 180.) vlon = vlon-360.

!  Count the number of volcanoes on your given day
   nv = 0
   do it = 1, nvolc
    if(nymd .lt. startday(it) .or. nymd .gt. endday(it)) cycle
    nv = nv + 1
   end do

!  Allocate space for the volcanoes
   allocate(vLatP(nv), vLonP(nv), &
            vSO2P(nv), vElevP(nv), &
            vCloudP(nv), stat=ios)

!  Accumulate the volcanoes
   nv = 0
   do it = 1, nvolc
    if(nymd .lt. startday(it) .or. nymd .gt. endday(it)) cycle
    nv = nv + 1
    vLatP(nv) = vlat(it)
    vLonP(nv) = vlon(it)
    vSO2P(nv) = so2exp(it) * 1.e6 / 86400.   ! to kg SO2/sec
    vElevP(nv) = velev(it)
    vCloudP(nv) = celev(it)
   enddo

   nVolcExp = nv

  end subroutine GetVolcExplosive

!-------------------------------------------------------------------------
!     NASA/GSFC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  GetVolcContinuous - Get pointwise SO2 and altitude of 
!                                 continuous outgassing volcanos from 
!                                 a data base

!
! !INTERFACE:
!

   subroutine GetVolcContinuous( nVolcC, &
                                 vLatP, vLonP, vElevP, vCloudP, vSO2P )
			  
! !USES:

  implicit NONE

! Description
! Data for outgassing volcanos provided by Thomas Diehl.  Data table is
! Mg SO2 day-1 and I convert to kg SO2 s-1 needed in emissions.  Assumption
! is continuous emissions throughout day.
! What is returned is the number of volcanoes and an array of locations,
! elevations, and SO2 amounts.

  real, pointer, dimension(:)   :: vLatP, vLonP, vElevP, vCloudP, vSO2P
  integer                       :: it, ios, nVolcC

! database parameters and data
  integer, parameter :: nvolc = 47
  real    :: vso2(nvolc), vlon(nvolc), vlat(nvolc), & 
             velev(nvolc), celev(nvolc)

! If previous instance of volcano point data tables exist, deallocate it
! to get the correct number of elements
  if(associated(vLatP))    deallocate(vLatP, stat=ios)
  if(associated(vLonP))    deallocate(vLonP, stat=ios)
  if(associated(vSO2P))    deallocate(vSO2P, stat=ios)
  if(associated(vElevP))   deallocate(vElevP, stat=ios)
  if(associated(vCloudP))  deallocate(vCloudP, stat=ios)

  data vso2 / &
        730,           44,         4000,           21,           16, &
        520,          920,          690,          480,         3300, &
        900,           75,           58,          140,           14, &
        370,          530,          570,         1900,          130, &
         76,          140,          370,          270,           56, &
         68,            3,           48,           22,         1027, &
        140,          230,          640,          510,           20, &
         20,          590,           84,           73,          790, &
        110,          500,         1900,          650,         2400, &
          3,           79 &
    /

  data velev / &
        926,          500,         3350,          613,         2890, &
        321,         1807,         1330,         2334,         1750, &
        361,         2084,         3432,         2911,         2329, &
       1565,         2462,          717,         1117,         1359, &
       1592,         1788,         2560,          758,          731, &
       1124,         1860,         1252,         3053,         1222, &
       4100,         3772,         3763,         2552,         2365, &
       1950,         1745,         1010,         1258,          635, &
       1657,         2708,         5321,         4276,         5592, &
       1920,         3794 &
   /
  data celev / &
        926,          500,         3350,          613,         2890, &
        321,         1807,         1330,         2334,         1750, &
        361,         2084,         3432,         2911,         2329, &
       1565,         2462,          717,         1117,         1359, &
       1592,         1788,         2560,          758,          731, &
       1124,         1860,         1252,         3053,         1222, &
       4100,         3772,         3763,         2552,         2365, &
       1950,         1745,         1010,         1258,          635, &
       1657,         2708,         5321,         4276,         5592, &
       1920,         3794 &
    /
  data vlon / &
    15.2130,      14.9620,      15.0040,      40.6700,      35.9020, &
   177.1800,     145.0610,     148.4200,     151.3300,     155.1950, &
   169.4250,     107.6000,     109.2080,     110.4420,     112.9500, &
   124.0500,     123.6850,     130.3080,     130.6570,     130.2940, &
   131.1060,     131.2510,     138.5260,     139.3980,     140.8430, &
   148.8430,    -155.3610,    -153.4300,    -153.0900,    -155.2920, &
  -103.6200,     -91.5520,     -90.8800,     -90.6010,     -89.6300, &
   -89.6330,     -87.0040,     -86.8450,     -86.5400,     -86.1610, &
   -84.7030,     -84.2330,     -75.3220,     -77.3700,     -67.7300, &
   -16.7200,     167.1700 &
    /
  data vlat / &
    38.7890,      38.4040,      37.7340,      13.6000,      -2.7510, &
   -37.5200,      -4.1000,      -5.5250,      -5.0500,      -6.1400, &
   -19.5200,      -6.7700,      -7.2420,      -7.5420,      -7.9420, &
    12.7700,      13.2570,      30.7890,      31.5850,      32.7570, &
    32.8810,      33.0830,      36.4030,      34.7210,      42.5410, &
    45.3870,      58.1720,      59.3630,      60.0320,      19.4250, &
    19.5140,      14.7560,      14.4730,      14.3810,      13.8530, &
    13.8130,      12.7020,      12.6020,      12.4220,      11.9840, &
    10.4630,      10.2000,       4.8950,       1.2200,     -23.3700, &
    64.6500,     -77.5300 &
    /

!  Allocate space for the volcanoes
   allocate(vLatP(nvolc), vLonP(nvolc), &
            vSO2P(nvolc), vElevP(nvolc), &
            vCloudP(nvolc), stat=ios)

!  Accumulate the volcanoes
   do it = 1, nvolc
    vLatP(it) = vlat(it)
    vLonP(it) = vlon(it)
    vSO2P(it) = vso2(it) * 1000. / 86400.  ! to kg SO2/sec
    vElevP(it) = velev(it)
    vCloudP(it) = celev(it)
   enddo
   nVolcC = nvolc

  end subroutine GetVolcContinuous

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  idaynum -- given nymd compute the day number of the year
!
! Colarco, July 29, 2004

   integer function idaynum (nymd)
   integer :: nymd, yyyy, mm, dd, imon, isleapyr
   integer :: ndays(12)

   data ndays /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/

   yyyy = nymd / 10000
   mm = mod(nymd,10000) / 100
   dd = mod(nymd,100)

!  Is it a leap year?
   isleapyr = 0
   if(mod(yyyy,4) .eq. 0) then
    isleapyr = 1
    if(mod(yyyy,100) .eq. 0) then
     isleapyr = 0
     if(mod(yyyy,400) .eq. 0) then
      isleapyr = 1
     endif
    endif
   endif

!  What day number
   idaynum = 0
   if(mm .eq. 1) then
    idaynum = dd
   else
    do imon = 1, mm-1
     if(imon .eq. 2 .and. isleapyr .eq. 1) then
      idaynum = idaynum+29
     else
      idaynum = idaynum + ndays(imon)
     endif
    enddo
    idaynum = idaynum + dd
   endif

   return
   end function idaynum



!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  szangle -- given locations and hour find the sza
!                        from GOCART (source?)
!
! Colarco, July 29, 2004

   subroutine szangle(jday,xhour,lonRad,latRad,sza,cossza,i1,i2,j1,j2)

   real, pointer, dimension(:,:) :: lonRad, latRad
   integer :: jday, i1, i2, j1, j2, i, j
   real :: a0, a1, a2, a3, b1, b2, b3, r, dec
   real :: pi, timloc, ahr, xlon, rlat, xHour
   real :: cossza(i1:i2,j1:j2), sza(i1:i2,j1:j2)
   data pi / 3.1415926 /

   a0 = 0.006918
   a1 = 0.399912
   a2 = 0.006758
   a3 = 0.002697
   b1 = 0.070257
   b2 = 0.000907
   b3 = 0.000148
   r  = 2.*pi*float(jday-1)/365. ! where jday is day # of the year

!  dec is the solar declination in radians
   dec = a0 - a1*cos(   r) + b1*sin(   r) &
            - a2*cos(2.*r) + b2*sin(2.*r) &
            - a3*cos(3.*r) + b3*sin(3.*r)

   do j = j1, j2
     do i = i1, i2
!    timloc is the local time in hours
     xlon = lonRad(i,j)*radToDeg
     timloc = xhour + xlon/15.
     if(timloc .lt. 0.)  timloc = timloc+24.
     if(timloc .gt. 24.) timloc = timloc-24.
!    ahr is the hour angle in radians
     ahr = abs(timloc - 12.)*15.*pi/180.

      rlat = latRad(i,j)
      cossza(i,j) =   sin(rlat)*sin(dec) &
                    + cos(rlat)*cos(dec)*cos(ahr)
      cossza(i,j)    = min(max(cossza(i,j),-1.0),1.0) !ALT make sure cos stays between -1.0 and 1.0
      sza(i,j)    = acos(cossza(i,j)) * radToDeg
      if(cossza(i,j) .lt. 0.) cossza(i,j) = 0.
     end do
   end do

   end subroutine szangle



end module 
