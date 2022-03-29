   subroutine SU_Wet_Removal ( km, nbins, klid, cdt, kin, grav, airMolWght, delp, fMassSO4, fMassSO2, &
                               h2o2_int, ple, rhoa, precc, precl, pfllsan, pfilsan, tmpu, &
                               nDMS, nSO2, nSO4, nMSA, DMS, SO2, SO4, MSA, &
                               fluxout, pSO4_colflux, pSO4wet_colflux, &
                               pso4, pso4wet, rc )


! !USES:
   implicit NONE

! !INPUT PARAMETERS:

   integer, intent(in) :: km, nbins  ! number of model levels and number of species respectively
   integer, intent(in) :: klid  ! index for pressure lid
   real, intent(in)    :: cdt   ! chemisty model timestep
   logical, intent(in) :: KIN   ! true for aerosol
   real, intent(in)    :: grav  ! gravity [m/sec]
   real, intent(in)    :: airMolWght ! air molecular weight [kg]
   real, dimension(:,:,:), intent(in) :: delp   ! pressure thickness [Pa]
   real, intent(in) :: fMassSO4, fMassSO2
   real, dimension(:,:,:) :: h2o2_int
   real, pointer, dimension(:,:,:), intent(in) :: ple     ! level edge air pressure
   real, pointer, dimension(:,:,:), intent(in) :: rhoa    ! air density, [kg m-3]
   real, pointer, dimension(:,:), intent(in)   :: precc   ! total convective precip, [mm day-1]
   real, pointer, dimension(:,:), intent(in)   :: precl   ! total large-scale prec,  [mm day-1]
   real, pointer, dimension(:,:,:), intent(in) :: pfllsan ! 3D flux of liquid nonconvective precipitation [kg/(m^2 sec)]
   real, pointer, dimension(:,:,:), intent(in) :: pfilsan ! 3D flux of ice nonconvective precipitation [kg/(m^2 sec)]
   real, pointer, dimension(:,:,:), intent(in) :: tmpu    ! temperature, [K]
   integer, intent(in) :: nDMS, nSO2, nSO4, nMSA !index position of sulfates
   real, dimension(:,:,:), intent(inout) :: DMS ! [kg/kg]
   real, dimension(:,:,:), intent(inout) :: SO2 ! [kg/kg]
   real, dimension(:,:,:), intent(inout) :: SO4 ! [kg/kg]
   real, pointer, dimension(:,:,:), intent(inout) :: MSA ! [kg/kg]

! !OUTPUT PARAMETERS:
   real, pointer, dimension(:,:,:),intent(inout) :: fluxout
   real, pointer, dimension(:,:),intent(inout)  :: pSO4_colflux
   real, pointer, dimension(:,:),intent(inout)  :: pSO4wet_colflux
   real, pointer, dimension(:,:,:),intent(inout) :: pso4
   real, pointer, dimension(:,:,:),intent(inout) :: pso4wet
   integer, optional, intent(out)   :: rc    ! Error return code:
                                            !  0 - all is well


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
!
! !Local Variables
   integer :: status

   integer :: i1=1, j1=1, i2, j2
   integer :: dims(3)

   integer  ::  i, j, k, iit, n, LH, kk, ios
!   real :: pdog(i1:i2,j1:j2,km)      ! air mass factor dp/g [kg m-2]
   real, dimension(:,:,:), allocatable :: pdog           ! air mass factor dp/g [kg m-2]
   real*8 :: Td_ls, Td_cv              ! ls and cv timescales [s]
   real*8 :: pls, pcv, pac             ! ls, cv, tot precip [mm day-1]
   real*8 :: qls(km), qcv(km)          ! ls, cv portion of moisture tendency [kg m-3 s-1]
   real*8 :: qmx, qd, A                ! temporary variables on moisture
   real*8 :: F, B, BT                  ! temporary variables on cloud, freq.
   real*8, allocatable :: fd(:,:)      ! flux across layers [kg m-2]
   real*8, allocatable :: dpfli(:,:,:) !
   real*8, allocatable :: DC(:)        ! scavenge change in mass mixing ratio
!   real :: c_h2o(i1:i2,j1:j2,km), cldliq(i1:i2,j1:j2,km), cldice(i1:i2,j1:j2,km)
   real, dimension(:,:,:), allocatable :: c_h2o, cldliq, cldice
   real, parameter :: kb = 1.3807e-23 ! Boltzmann constant [kg m2 s-1 K-1 mol-1]

!  Rain parameters (from where?)
   real, parameter :: B0_ls = 1.0e-4
   real, parameter :: F0_ls = 1.0
   real, parameter :: XL_ls = 5.0e-4
   real, parameter :: B0_cv = 1.5e-3
   real, parameter :: F0_cv = 0.3
   real, parameter :: XL_cv = 2.0e-3
   real, parameter :: one = 1.0, zero = 0.0

   integer :: nbins_ = 0 ! nbins needs to be redefined in case MSA is not being computed

!  Conversion of SO2 mmr to SO2 vmr (since H2O2 is carried around like
!  a volume mixing ratio)
   real*8 :: fmr, SO2Soluble
   fMR = airMolWght / fMassSO2

!EOP
!-------------------------------------------------------------------------
!  Begin

   rc = __SUCCESS__

   allocate(c_h2o, mold=rhoa)
   allocate(cldliq, mold=rhoa)
   allocate(cldice, mold=rhoa)
!  Initialize local variables
!  --------------------------
!  c_h2o, cldliq, and cldice are respectively intended to be the
!  water mixing ratio (liquid or vapor?, in or out of cloud?)
!  cloud liquid water mixing ratio
!  cloud ice water mixing ratio
   c_h2o  = (10d0**(-2663.5d0/tmpu(:,:,:) + 12.537d0 ) ) /  &
                   (ple(:,:,0:km-1)+ple(:,:,1:km)) /2d0
   cldliq = 0.d0
   where(tmpu >  248.) cldliq = 1.d-6 * ( ( tmpu - 248.d0) / 20.d0 )
   where(tmpu >= 268.) cldliq = 1.d-6
   cldice = 1.d-6 - cldliq

   dims = shape(rhoa)
   i2 = dims(1); j2 = dims(2)

!  check if doing MSA and define nbins_ accordingly
!   if (associated(MSA)) then
!      nbins_ = nbins
!   else
!      nbins_ = nbins - 1
!   end if

   do n = 1, nbins
    if( associated(fluxout)) fluxout(:,:,n) = 0.0
   end do
   if( associated(pso4wet_colflux)) pso4wet_colflux(i1:i2,j1:j2) = 0.
   if( associated(pso4wet)) pso4wet(i1:i2,j1:j2,1:km) = 0.

!  Allocate the dynamic arrays
   allocate(fd(km,nbins),__STAT__)
   allocate(dc(nbins),__STAT__)
   allocate(dpfli(i1:i2, j1:j2, km),__STAT__)

!  Duration of rain: ls = model timestep, cv = 1800 s (<= cdt)
   Td_ls = cdt
   Td_cv = 1800.

!  Accumulate the 3-dimensional arrays of rhoa and pdog
   pdog = delp/grav

   dpfli = pfllsan(:,:,1:km)-pfllsan(:,:,0:km-1)+pfilsan(:,:,1:km)-pfilsan(:,:,0:km-1)

!  Loop over spatial indices
   do j = j1, j2
    do i = i1, i2

!    Check for total precipitation amount
!    Assume no precip in column if precl+precc = 0
     pac = precl(i,j) + precc(i,j)
     if(pac .le. 0.) goto 100
     pls = precl(i,j)
     pcv = precc(i,j)

!    Initialize the precipitation fields
     qls(:)  = 0.
     qcv(:)  = 0.
     fd(:,:) = 0.
     Dc(:)   = 0.

!    Find the highest model layer experiencing rainout.  Assumes no
!    scavenging if T < 258 K
     LH = 0
     do k = klid, km
      if(dpfli(i,j,k) .gt. 0. .and. tmpu(i,j,k) .gt. 258.) then
       LH = k
       goto 15
      endif
     end do
 15  continue
     if(LH .lt. 1) goto 100

     do k = LH, km
      qls(k) = dpfli(i,j,k)/pdog(i,j,k)*rhoa(i,j,k)
     end do

!    Loop over vertical to do the scavenging!
     do k = LH, km

!-----------------------------------------------------------------------------
!   (1) LARGE-SCALE RAINOUT:
!       Tracer loss by rainout = TC0 * F * exp(-B*dt)
!         where B = precipitation frequency,
!               F = fraction of grid box covered by precipitating clouds.
!       We assume that tracer scavenged by rain is falling down to the
!       next level, where a fraction could be re-evaporated to gas phase
!       if Qls is less then 0 in that level.
!-----------------------------------------------------------------------------
      if (qls(k) .gt. tiny(0.)) then
       F  = F0_ls / (1. + F0_ls*B0_ls*XL_ls/(qls(k)*cdt/Td_ls))
       B  = B0_ls/F0_ls +1./(F0_ls*XL_ls/qls(k))
       BT = B * Td_ls
       if (BT.gt.10.) BT = 10.               !< Avoid overflow >
!      What is the soluble amount of SO2?
       SO2Soluble = min(fmr*SO2(i,j,k),h2o2_int(i,j,k)*one)/fmr
       if(SO2Soluble .lt. 0.) SO2Soluble = 0.

!      Adjust SU amounts
       DC(nDMS) = 0.
       DC(nSO2) = SO2Soluble * F * (1.-exp(-BT))
       DC(nSO4) = SO4(i,j,k) * F * (1.-exp(-BT))
       if (associated(MSA)) then
          DC(nMSA) = MSA(i,j,k) * F * (1.-exp(-BT))
       endif

!      Adjust H2O2 concentration in cloudy portion of cell
       if(fmr*SO2(i,j,k) .gt. h2o2_int(i,j,k)) then
!        gcSU%h2o2_int(i,j,k) = max(zero,(1.-F)*gcSU%h2o2_int(i,j,k))
! GOCART removes all
        h2o2_int(i,j,k) = 0.
       else
        h2o2_int(i,j,k) &
          = h2o2_int(i,j,k) - F*fmr*SO2(i,j,k)
       endif

       do n = 1, nbins
        if (DC(n).lt.0.) DC(n) = 0.
       end do

       call updateAerosol(DMS(i,j,k), DC(nDMS))
       call updateAerosol(SO2(i,j,k), DC(nSO2))
       call updateAerosol(SO4(i,j,k), DC(nSO4))
       if (associated(MSA)) then
          call updateAerosol(MSA(i,j,k), DC(nMSA))
       end if

!      Flux down:  unit is kg m-2
!      Formulated in terms of production in the layer.  In the revaporation step
!      we consider possibly adding flux from above...
       do n = 1, nbins
        Fd(k,n) = DC(n) * pdog(i,j,k)
       end do

      end if                                    ! if Qls > 0  >>>

!-----------------------------------------------------------------------------
! * (2) LARGE-SCALE WASHOUT:
! *     Occurs when rain at this level is less than above.
!-----------------------------------------------------------------------------
      if(k .gt. LH .and. qls(k) .ge. 0.) then
       if(qls(k) .lt. qls(k-1)) then
!       Find a maximum F overhead until the level where Qls<0.
        Qmx   = 0.
        do kk = k-1,LH,-1
         if (Qls(kk).gt.0.) then
          Qmx = max(Qmx,Qls(kk))
         else
          goto 333
         end if
        end do

 333    continue
        F = F0_ls / (1. + F0_ls*B0_ls*XL_ls/(Qmx*cdt/Td_ls))
        if (F.lt.0.01) F = 0.01

!-----------------------------------------------------------------------------
!  The following is to convert Q(k) from kgH2O/m3/sec to mm/sec in order
!  to use the Harvard formula.  Convert back to mixing ratio by multiplying
!  by rhoa.  Multiply by pdog gives kg/m2/s of precip.  Divide by density
!  of water (=1000 kg/m3) gives m/s of precip and multiply by 1000 gives
!  units of mm/s (omit the multiply and divide by 1000).
!-----------------------------------------------------------------------------

        Qd = Qmx /rhoa(i,j,k)*pdog(i,j,k)
        if (Qd.ge.50.) then
         B = 0.
        else
         B = Qd * 0.1
        end if
        BT = B * cdt
        if (BT.gt.10.) BT = 10.

!      What is the soluble amount of SO2?
       SO2Soluble = min(fmr*SO2(i,j,k),h2o2_int(i,j,k)*one)/fmr
       if(SO2Soluble .lt. 0.) SO2Soluble = 0.

!      Adjust SU amounts
       DC(nDMS) = 0.
       DC(nSO2) = SO2Soluble * F * (1.-exp(-BT))
       DC(nSO4) = SO4(i,j,k) * F * (1.-exp(-BT))
       if (associated(MSA)) then
          DC(nMSA) = MSA(i,j,k) * F * (1.-exp(-BT))
       end if

!       Adjust H2O2 concentration in cloudy portion of cell
        if(fmr*SO2(i,j,k) .gt. h2o2_int(i,j,k)) then
         h2o2_int(i,j,k) = max(zero,(one-F)*h2o2_int(i,j,k))
!  GOCART removes all
!         gcSU%h2o2_int(i,j,k) = 0.
        else
         h2o2_int(i,j,k) &
           = h2o2_int(i,j,k) - F*fmr*SO2(i,j,k)
        endif

        do n = 1, nbins
         if (DC(n).lt.0.) DC(n) = 0.
        end do

       call updateAerosol(DMS(i,j,k), DC(nDMS))
       call updateAerosol(SO2(i,j,k), DC(nSO2))
       call updateAerosol(SO4(i,j,k), DC(nSO4))
       if (associated(MSA)) then
          call updateAerosol(MSA(i,j,k), DC(nMSA))
       end if

        do n = 1, nbins
         if( associated(fluxout) ) then
          fluxout(i,j,n) = fluxout(i,j,n)+DC(n)*pdog(i,j,k)/cdt
         endif
        end do

       end if
      end if                                    ! if ls washout  >>>

!-----------------------------------------------------------------------------
!  (3) CONVECTIVE RAINOUT:
!      Tracer loss by rainout = dd0 * F * exp(-B*dt)
!        where B = precipitation frequency,
!              F = fraction of grid box covered by precipitating clouds.
!-----------------------------------------------------------------------------
      if (qcv(k) .gt. 0.) then
       F  = F0_cv / (1. + F0_cv*B0_cv*XL_cv/(Qcv(k)*cdt/Td_cv))
       B  = B0_cv
       BT = B * Td_cv
       if (BT.gt.10.) BT = 10.               !< Avoid overflow >

!      Adjust SO2 for H2O2 oxidation
       SO2Soluble = min(fmr*SO2(i,j,k),h2o2_int(i,j,k)*one)/fmr
       if(SO2Soluble .lt. 0.) SO2Soluble = 0.

!      Adjust SU amounts
       DC(nDMS) = 0.
       DC(nSO2) = SO2Soluble * F * (1.-exp(-BT))
       DC(nSO4) = SO4(i,j,k) * F * (1.-exp(-BT))
       if (associated(MSA)) then
          DC(nMSA) = MSA(i,j,k) * F * (1.-exp(-BT))
       end if
       DC(nSO4) = 0.
       if (associated(MSA)) then
          DC(nMSA) = 0.
       end if

!      Adjust H2O2 concentration in cloudy portion of cell
       if(fmr*SO2(i,j,k) .gt. h2o2_int(i,j,k)) then
        h2o2_int(i,j,k) = max(zero,(one-F)*h2o2_int(i,j,k))
       else
        h2o2_int(i,j,k) &
          = h2o2_int(i,j,k) - F*fmr*SO2(i,j,k)
       endif

       do n = 1, nbins
        if (DC(n).lt.0.) DC(n) = 0.
       end do

       call updateAerosol(DMS(i,j,k), DC(nDMS))
       call updateAerosol(SO2(i,j,k), DC(nSO2))
       call updateAerosol(SO4(i,j,k), DC(nSO4))
       if (associated(MSA)) then
          call updateAerosol(MSA(i,j,k), DC(nMSA))
       end if

!      Flux down:  unit is kg m-2
!      Formulated in terms of production in the layer.  In the revaporation step
!      we consider possibly adding flux from above...
       do n = 1, nbins
        Fd(k,n) = Fd(k,n) + DC(n)*pdog(i,j,k)
       end do

      end if                                  ! if Qcv > 0   >>>

!-----------------------------------------------------------------------------
!  (4) CONVECTIVE WASHOUT:
!      Occurs when rain at this level is less than above.
!-----------------------------------------------------------------------------

      if (k.gt.LH .and. Qcv(k).ge.0.) then
       if (Qcv(k).lt.Qcv(k-1)) then
!-----  Find a maximum F overhead until the level where Qls<0.
        Qmx   = 0.
        do kk = k-1, LH, -1
         if (Qcv(kk).gt.0.) then
          Qmx = max(Qmx,Qcv(kk))
         else
          goto 444
         end if
        end do

 444    continue
        F = F0_cv / (1. + F0_cv*B0_cv*XL_cv/(Qmx*cdt/Td_cv))
        if (F.lt.0.01) F = 0.01
!-----------------------------------------------------------------------------
!  The following is to convert Q(k) from kgH2O/m3/sec to mm/sec in order
!  to use the Harvard formula.  Convert back to mixing ratio by multiplying
!  by rhoa.  Multiply by pdog gives kg/m2/s of precip.  Divide by density
!  of water (=1000 kg/m3) gives m/s of precip and multiply by 1000 gives
!  units of mm/s (omit the multiply and divide by 1000).
!-----------------------------------------------------------------------------

        Qd = Qmx / rhoa(i,j,k)*pdog(i,j,k)
        if (Qd.ge.50.) then
         B = 0.
        else
         B = Qd * 0.1
        end if
        BT = B * cdt
        if (BT.gt.10.) BT = 10.

!       Adjust SO2 for H2O2 oxidation
        SO2Soluble = min(fmr*SO2(i,j,k),h2o2_int(i,j,k)*one)/fmr
        if(SO2Soluble .lt. 0.) SO2Soluble = 0.

!       Adjust SU amounts
        DC(nDMS) = 0.
        DC(nSO2) = SO2Soluble * F * (1.-exp(-BT))
! Sulfate scavenged in moist
!        DC(nSO4) = w_c%qa(n1+nSO4-1)%data3d(i,j,k) * F * (1.-exp(-BT))
!        DC(nMSA) = w_c%qa(n1+nMSA-1)%data3d(i,j,k) * F * (1.-exp(-BT))
        DC(nSO4) = 0.
        if (associated(MSA)) then
           DC(nMSA) = 0.
        end if

!       Adjust H2O2 concentration in cloudy portion of cell
        if(fmr*SO2(i,j,k) .gt. h2o2_int(i,j,k)) then
         h2o2_int(i,j,k) = max(zero,(one-F)*h2o2_int(i,j,k))
        else
         h2o2_int(i,j,k) &
           = h2o2_int(i,j,k) - F*fmr*SO2(i,j,k)
        endif

        do n = 1, nbins
         if (DC(n).lt.0.) DC(n) = 0.
        end do

       call updateAerosol(DMS(i,j,k), DC(nDMS))
       call updateAerosol(SO2(i,j,k), DC(nSO2))
       call updateAerosol(SO4(i,j,k), DC(nSO4))
       if (associated(MSA)) then
          call updateAerosol(MSA(i,j,k), DC(nMSA))
       end if

        do n = 1, nbins
         if( associated(fluxout) ) then
          fluxout(i,j,n) = fluxout(i,j,n)+DC(n)*pdog(i,j,k)/cdt
         endif
        end do

       end if
      end if                                    ! if cv washout  >>>

!-----------------------------------------------------------------------------
!  (5) RE-EVAPORATION.  Assume that SO2 is re-evaporated as SO4 since it
!      has been oxidized by H2O2 at the level above.
!-----------------------------------------------------------------------------
!     Add in the flux from above, which will be subtracted if reevaporation occurs
      if(k .gt. LH) then
       do n = 1, nbins
        Fd(k,n) = Fd(k,n) + Fd(k-1,n)
       end do

!      Is there evaporation in the currect layer?
       if (dpfli(i,j,k) .lt. 0.) then
!       Fraction evaporated = H2O(k)evap / H2O(next condensation level).
        if (dpfli(i,j,k-1) .gt. 0.) then

          A =  abs(  dpfli(i,j,k) /  dpfli(i,j,k-1)  )
          if (A .gt. 1.) A = 1.

!         Adjust tracer in the level
!         For the SO2 tracer we do not allow re-evaporation.
!         We compute DC(nSO2) solely to add this to DC(nSO4) and to remove
!         from Fd(k,nSO2)
!         Instead, the SO2 gets re-evaporated to the SO4 bin because of
!         previous H2O2 oxidation

          DC(nDMS) = 0.
          DC(nSO2) = Fd(k-1,nSO2) / pdog(i,j,k) * A
          DC(nSO4) = Fd(k-1,nSO4) / pdog(i,j,k) * A
          if (associated(MSA)) then
             DC(nMSA) = Fd(k-1,nMSA) / pdog(i,j,k) * A
          end if

          do n = 1, nbins
           if (DC(n).lt.0.) DC(n) = 0.
          end do

          if (associated(MSA)) then
             MSA(i,j,k) = MSA(i,j,k) + DC(nMSA)
          end if
!         SO2 gets added to SO4, but remember to remove the SO2 from FD!
          SO4(i,j,k) = SO4(i,j,k) + DC(nSO4) + DC(nSO2)*fMassSO4/fMassSO2
          if( associated(pso4wet_colflux)) &
             pso4wet_colflux(i,j) = pso4wet_colflux(i,j) &
              + DC(nSO2)*fMassSO4/fMassSO2 / cdt * delp(i,j,k)/grav
          if( associated(pso4wet) ) &
             pso4wet(i,j,k) = DC(nSO2)*fMassSO4/fMassSO2 / cdt

          if( associated(pso4_colflux)) &
             pso4_colflux(i,j) = pso4_colflux(i,j) &
              + DC(nSO2)*fMassSO4/fMassSO2 / cdt * delp(i,j,k)/grav
          if( associated(pso4) ) &
             pso4(i,j,k) = pso4(i,j,k) + DC(nSO2)*fMassSO4/fMassSO2 / cdt


!         Adjust the flux out of the bottom of the layer--remove SO2 here!
          DMS(i,j,k) = max(DMS(i,j,k),tiny(1.0))
          Fd(k,nDMS) = Fd(k,nDMS) - DC(nDMS)*pdog(i,j,k)
          SO2(i,j,k) = max(SO2(i,j,k),tiny(1.0))
          Fd(k,nSO2) = Fd(k,nSO2) - DC(nSO2)*pdog(i,j,k)
          SO4(i,j,k) = max(SO4(i,j,k),tiny(1.0))
          Fd(k,nSO4) = Fd(k,nSO4) - DC(nSO4)*pdog(i,j,k)
          if (associated(MSA)) then
             MSA(i,j,k) = max(MSA(i,j,k),tiny(1.0))
             Fd(k,nMSA) = Fd(k,nMSA) - DC(nMSA)*pdog(i,j,k)
          end if
        endif
       endif                                   ! if -moistq < 0
      endif
     end do  ! k

     do n = 1, nbins
      if( associated(fluxout) ) then
       fluxout(i,j,n) = fluxout(i,j,n)+Fd(km,n)/cdt
      endif
     end do

 100 continue
    end do   ! i
   end do    ! j

   deallocate(fd,DC,dpfli,stat=ios)


   __RETURN__(__SUCCESS__)
   contains
     subroutine updateAerosol (aerosol, DC)

     ! !USES:
     implicit NONE
     ! !INPUT PARAMETERS:
      real, intent(inout) :: aerosol
      real*8, intent(in)  :: DC

        aerosol = aerosol - DC
        if (aerosol .lt. 1.0E-32) aerosol = 1.0E-32

      end subroutine updateAerosol

   end subroutine SU_Wet_Removal
