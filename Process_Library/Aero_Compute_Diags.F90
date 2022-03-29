   subroutine Aero_Compute_Diags (mie, km, klid, nbegin, nbins, rlow, rup, &
                                  wavelengths_profile, wavelengths_vertint, aerosol, &
                                  grav, tmpu, rhoa, rh, u, v, delp, ple,tropp, &
                                  sfcmass, colmass, mass, exttau, stexttau, scatau, stscatau,&
                                  sfcmass25, colmass25, mass25, exttau25, scatau25, &
                                  fluxu, fluxv, conc, extcoef, scacoef, &
                                  exttaufm, scataufm, angstrom, aerindx, NO3nFlag, rc )

! !USES:

  implicit NONE

! !INPUT PARAMETERS:
   type(GOCART2G_Mie),  intent(in) :: mie        ! mie table
   integer, intent(in) :: km, nbegin, nbins
   integer,    intent(in)    :: klid   ! index for pressure lid
   real, optional, dimension(:), intent(in)    :: rlow   ! bin radii - low bounds
   real, optional, dimension(:), intent(in)    :: rup    ! bin radii - upper bounds
   real, dimension(:), intent(in)    :: wavelengths_profile
   real, dimension(:), intent(in)    :: wavelengths_vertint
   real, dimension(:,:,:,:), intent(in) :: aerosol     !
   real, intent(in) :: grav
   real, pointer, dimension(:,:,:), intent(in) :: tmpu  ! temperature [K]
   real, pointer, dimension(:,:,:), intent(in) :: rhoa  ! air density [kg/m^3]
   real, pointer, dimension(:,:,:), intent(in) :: delp  ! pressure thickness [Pa]
   real, pointer, dimension(:,:,:), intent(in) :: rh    ! relative humidity [1]
   real, pointer, dimension(:,:,:), intent(in) :: u     ! east-west wind [m/s]
   real, pointer, dimension(:,:,:), intent(in) :: v     ! north-south wind [m/s]
   real, pointer, dimension(:,:,:), intent(in) :: ple   ! level edge air pressure [Pa]
   real, pointer, dimension(:,:), intent(in)   :: tropp ! tropopause pressure [Pa]
   logical, optional, intent(in)               :: NO3nFlag

! !OUTPUT PARAMETERS:
!  Total mass
   real, optional, dimension(:,:), intent(inout)   :: sfcmass   ! sfc mass concentration kg/m3
   real, optional, dimension(:,:), intent(inout)   :: colmass   ! col mass density kg/m2
   real, pointer, dimension(:,:,:), intent(inout) :: mass      ! 3d mass mixing ratio kg/kg
   real, pointer, dimension(:,:,:), intent(inout) :: conc      ! 3d mass concentration, kg/m3
!  Total optical properties
   real, optional, dimension(:,:,:), intent(inout)   :: exttau    ! ext. AOT at 550 nm
   real, optional, dimension(:,:,:), intent(inout)   :: stexttau  ! stratospheric ext. AOT at 550 nm
   real, optional, dimension(:,:,:), intent(inout)   :: scatau    ! sct. AOT at 550 nm
   real, optional, dimension(:,:,:), intent(inout)   :: stscatau  ! stratospheric sct. AOT at 550 nm
   real, optional, dimension(:,:), intent(inout)   :: sfcmass25 ! sfc mass concentration kg/m3 (pm2.5)
   real, optional, dimension(:,:), intent(inout)   :: colmass25 ! col mass density kg/m2 (pm2.5)
   real, optional, dimension(:,:,:), intent(inout) :: mass25    ! 3d mass mixing ratio kg/kg (pm2.5)
   real, optional, dimension(:,:,:), intent(inout)   :: exttau25  ! ext. AOT at 550 nm (pm2.5)
   real, optional, dimension(:,:,:), intent(inout)   :: scatau25  ! sct. AOT at 550 nm (pm2.5)
   real, optional, dimension(:,:),  intent(inout)  :: aerindx   ! TOMS UV AI
   real, optional, dimension(:,:), intent(inout)   :: fluxu     ! Column mass flux in x direction
   real, optional, dimension(:,:), intent(inout)   :: fluxv     ! Column mass flux in y direction
   real, optional, dimension(:,:,:,:), intent(inout) :: extcoef   ! 3d ext. coefficient, 1/m
   real, optional, dimension(:,:,:,:), intent(inout) :: scacoef   ! 3d scat.coefficient, 1/m
   real, optional, dimension(:,:,:), intent(inout)   :: exttaufm  ! fine mode (sub-micron) ext. AOT at 550 nm
   real, optional, dimension(:,:,:), intent(inout)   :: scataufm  ! fine mode (sub-micron) sct. AOT at 550 nm
   real, optional, dimension(:,:), intent(inout)   :: angstrom  ! 470-870 nm Angstrom parameter
   integer, optional, intent(out)   :: rc        ! Error return code:
                                                 !  0 - all is well
                                                 !  1 -

! !DESCRIPTION: Calculates some simple 2d diagnostics from the dust fields
!
! !REVISION HISTORY:
!
!  16APR2004, Colarco
!  11MAR2010, Nowottnick
!  11AUG2020, E.Sherman - refactored to work for multiple aerosols

! !Local Variables
   character(len=*), parameter :: myname = 'Aero_Compute_Diags'
   integer :: i, j, k, n, w, ios, status
   integer :: i1 =1, i2, j1=1, j2
   integer :: ilam470, ilam870
   real, allocatable, dimension(:,:,:) :: tau, ssa
!   real :: fPMfm(nbins)  ! fraction of bin with particles diameter < 1.0 um
!   real :: fPM25(nbins)  ! fraction of bin with particles diameter < 2.5 um
   real, dimension(:), allocatable :: fPMfm  ! fraction of bin with particles diameter < 1.0 um
   real, dimension(:), allocatable :: fPM25  ! fraction of bin with particles diameter < 2.5 um
   logical :: do_angstrom
   real, dimension(:,:), allocatable :: tau470, tau870
   logical   :: NO3nFlag_ !local version of the input

!EOP
!-------------------------------------------------------------------------
!  Begin...

   if( present(NO3nFlag) ) then
      NO3nFlag_ = NO3nFlag
   else
      NO3nFlag_ = .false.
   end if

!  Initialize local variables
!  --------------------------
   i2 = size(rhoa,1)
   j2 = size(rhoa,2)
   allocate(fPMfm(nbins))
   allocate(fPM25(nbins))

!  Get the wavelength indices
!  --------------------------

   ilam470 = mie%getChannel(4.70e-7)
   if(ilam470 <= 0) ilam470 = 0

   ilam870 = mie%getChannel(8.70e-7)
   if(ilam870 <= 0) ilam870 = 0

!  Determine if going to do Angstrom parameter calculation
!  -------------------------------------------------------
   do_angstrom = .false.
!  If both 470 and 870 channels provided (and not the same) then
!  possibly will do Angstrom parameter calculation
   if(ilam470 .ne. 0 .and. &
      ilam870 .ne. 0 .and. &
      ilam470 .ne. ilam870) do_angstrom = .true.

   if( present(angstrom) .and. do_angstrom ) then
      allocate(tau470(i1:i2,j1:j2), tau870(i1:i2,j1:j2))
   end if

!  Compute the fine mode (sub-micron) and PM2.5 bin-wise fractions
!  ------------------------------------
   if (present(rlow) .and. present(rup)) then
      call Aero_Binwise_PM_Fractions(fPMfm, 0.50, rlow, rup, nbins)   ! 2*r < 1.0 um
      call Aero_Binwise_PM_Fractions(fPM25, 1.25, rlow, rup, nbins)   ! 2*r < 2.5 um
   end if

   if (present(aerindx))  aerindx = 0.0  ! for now

!  Calculate the diagnostic variables if requested
!  -----------------------------------------------
!  Calculate the surface mass concentration
   if( present(sfcmass) ) then
      sfcmass(i1:i2,j1:j2) = 0.
      do n = nbegin, nbins
         sfcmass(i1:i2,j1:j2) &
              =   sfcmass(i1:i2,j1:j2) &
              + aerosol(i1:i2,j1:j2,km,n)*rhoa(i1:i2,j1:j2,km)
      end do
   endif
   if( present(sfcmass25) ) then
      sfcmass25(i1:i2,j1:j2) = 0.
      do n = nbegin, nbins
         sfcmass25(i1:i2,j1:j2) &
              =   sfcmass25(i1:i2,j1:j2) &
              + aerosol(i1:i2,j1:j2,km,n)*rhoa(i1:i2,j1:j2,km)*fPM25(n)
      end do
   endif

!  Calculate the aerosol column loading
   if( present(colmass) ) then
      colmass(i1:i2,j1:j2) = 0.
      do n = nbegin, nbins
       do k = klid, km
        colmass(i1:i2,j1:j2) &
         =   colmass(i1:i2,j1:j2) &
           + aerosol(i1:i2,j1:j2,k,n)*delp(i1:i2,j1:j2,k)/grav
       end do
      end do
   endif
   if( present(colmass25) ) then
      colmass25(i1:i2,j1:j2) = 0.
      do n = nbegin, nbins
         do k = klid, km
            colmass25(i1:i2,j1:j2) &
             = colmass25(i1:i2,j1:j2) &
             + aerosol(i1:i2,j1:j2,k,n)*delp(i1:i2,j1:j2,k)/grav*fPM25(n)
       end do
      end do
   endif

!  Calculate the total mass concentration
   if( associated(conc) ) then
      conc(i1:i2,j1:j2,1:km) = 0.
      do n = nbegin, nbins
         conc(i1:i2,j1:j2,1:km) &
             = conc(i1:i2,j1:j2,1:km) &
             + aerosol(i1:i2,j1:j2,1:km,n)*rhoa(i1:i2,j1:j2,1:km)
      end do
   endif

!  Calculate the total mass mixing ratio
   if( associated(mass) ) then
      mass(i1:i2,j1:j2,1:km) = 0.
      do n = nbegin, nbins
       mass(i1:i2,j1:j2,1:km) &
         =   mass(i1:i2,j1:j2,1:km) &
           + aerosol(i1:i2,j1:j2,1:km,n)
      end do
   endif
   if( present(mass25) ) then
      mass25(i1:i2,j1:j2,1:km) = 0.
      do n = nbegin, nbins
       mass25(i1:i2,j1:j2,1:km) &
         =   mass25(i1:i2,j1:j2,1:km) &
           + aerosol(i1:i2,j1:j2,1:km,n)*fPM25(n)
      end do
   endif

!  Calculate the column mass flux in x direction
   if( present(fluxu) ) then
      fluxu(i1:i2,j1:j2) = 0.
      do n = nbegin, nbins
         do k = klid, km
           fluxu(i1:i2,j1:j2) &
            = fluxu(i1:i2,j1:j2) &
            + aerosol(i1:i2,j1:j2,k,n)*delp(i1:i2,j1:j2,k)/grav*u(i1:i2,j1:j2,k)
         end do
      end do
   endif

!  Calculate the column mass flux in y direction
   if( present(fluxv) ) then
      fluxv(i1:i2,j1:j2) = 0.
      do n = nbegin, nbins
         do k = klid, km
           fluxv(i1:i2,j1:j2) &
           = fluxv(i1:i2,j1:j2) &
           + aerosol(i1:i2,j1:j2,k,n)*delp(i1:i2,j1:j2,k)/grav*v(i1:i2,j1:j2,k)
         end do
      end do
   endif

   allocate(tau(i1:i2,j1:j2,km),source = 0.)
   allocate(ssa(i1:i2,j1:j2,km),source = 0.)
!  Calculate the extinction and/or scattering AOD
   if( present(extcoef)  .or. &
       present(scacoef) ) then

      if( present(extcoef) ) extcoef = 0.
      if( present(scacoef) ) scacoef = 0.

      do n = nbegin, nbins
        do w = 1, size(wavelengths_profile)
          call mie%Query(wavelengths_profile(w),n,   &
                         aerosol(:,:,:,n)*delp/grav, &
                         rh, tau=tau, ssa=ssa, __RC__)
!         Calculate the total ext. and scat. coefficients
          if ( present(extcoef) ) then
             extcoef(:,:,:,w) = extcoef(:,:,:,w) + &
                               tau * (grav * rhoa / delp)
          endif
          if ( present(scacoef) ) then
             scacoef(:,:,:,w) = scacoef(:,:,:,w) + &
                               ssa * tau * (grav * rhoa / delp)
          endif
       enddo !wavelengths_profile
      enddo !nbins
    end if !present(extcoef)...

   if( present(exttau) .or. &
       present(stexttau) .or. &
       present(scatau) .or. &
       present(stscatau) .or. &
       present(exttau25) .or. &
       present(exttaufm) .or. &
       present(scatau25) .or. &
       present(scataufm) ) then

      if( present(exttau) ) exttau = 0.
      if( present(stexttau) ) stexttau = 0.
      if( present(scatau) ) scatau = 0.
      if( present(stscatau) ) stscatau = 0.

      if( present(exttau25) ) exttau25 = 0.
      if( present(scatau25) ) scatau25 = 0.

      if( present(exttaufm) ) exttaufm = 0.
      if( present(scataufm) ) scataufm = 0.

      do w = 1, size(wavelengths_vertint)
        do n = nbegin, nbins
           call mie%Query(wavelengths_vertint(w), n,  &
                          aerosol(:,:,:,n)*delp/grav, &
                          rh, tau=tau, ssa=ssa, __RC__)
           do k = klid, km
!             Integrate in the vertical
              if( present(exttau) ) exttau(:,:,w) = exttau(:,:,w) + tau(:,:,k)
              if( present(stexttau) ) then
                 where (ple(:,:,k) .le. tropp) 
                    stexttau(:,:,w) = stexttau(:,:,w) + tau(:,:,k)
                 elsewhere(ple(:,:,k) .gt. tropp .and. ple(:,:,k-1) .lt. tropp) 
                    stexttau(:,:,w) = stexttau(:,:,w) + log(tropp/ple(:,:,k-1))/log(ple(:,:,k)/ple(:,:,k-1))*tau(:,:,k)
                 endwhere
              endif

              if( present(exttaufm) ) then
                 if( NO3nFlag_ ) then
                    exttaufm(:,:,w) = exttaufm(:,:,w) + tau(:,:,k)
                 else
                    exttaufm(:,:,w) = exttaufm(:,:,w) + tau(:,:,k) * fPMfm(n)
                 end if
              end if

              if( present(exttau25) ) then
                 if( NO3nFlag_ ) then
                    exttau25(:,:,w) = exttau25(:,:,w) + tau(:,:,k)
                 else
                    exttau25(:,:,w) = exttau25(:,:,w) + tau(:,:,k) * fPM25(n)
                 end if
              end if

              if( present(scatau) ) scatau(:,:,w) = scatau(:,:,w) + tau(:,:,k)*ssa(:,:,k)
              if( present(stscatau) ) then
                 where (ple(:,:,k) .le. tropp) 
                    stscatau(:,:,w) = stscatau(:,:,w) + tau(:,:,k)*ssa(:,:,k)
                 elsewhere(ple(:,:,k) .gt. tropp .and. ple(:,:,k-1) .lt. tropp) 
                    stscatau(:,:,w) = stscatau(:,:,w) + log(tropp/ple(:,:,k-1))/log(ple(:,:,k)/ple(:,:,k-1))*tau(:,:,k)*ssa(:,:,k)
                 endwhere
              endif
              if( present(scataufm) ) then
                 if( NO3nFlag_ ) then
                    scataufm(:,:,w) = scataufm(:,:,w) + tau(:,:,k) * ssa(:,:,k)
                 else
                    scataufm(:,:,w) = scataufm(:,:,w) + tau(:,:,k) * ssa(:,:,k) * fPMfm(n)
                 end if
              end if

              if( present(scatau25) ) then
                 if( NO3nFlag_ ) then
                    scatau25(:,:,w) = scatau25(:,:,w) + tau(:,:,k) * ssa(:,:,k)
                 else
                    scatau25(:,:,w) = scatau25(:,:,w) + tau(:,:,k) * ssa(:,:,k) * fPM25(n)
                 end if
              end if
           enddo !k
        enddo !wavelengths_vertint
      enddo !nbins
   endif !present(exttau)...

!  Calculate the 470-870 Angstrom parameter
   if( present(angstrom) .and. do_angstrom ) then

      angstrom(i1:i2,j1:j2) = 0.
!     Set tau to small number by default
      tau470(i1:i2,j1:j2) = tiny(1.0)
      tau870(i1:i2,j1:j2) = tiny(1.0)

      do n = nbegin, nbins

!       Select the name for species
        call mie%Query(4.70E-7, n,                 &
                       aerosol(:,:,:,n)*delp/grav, &
                       rh, tau=tau, __RC__)
        do k = klid, km
          tau470 = tau470 + tau(:,:,k)
        enddo

        call mie%Query(8.70E-7, n,                 &
                       aerosol(:,:,:,n)*delp/grav, &
                       rh, tau=tau, __RC__)
        do k = klid, km
          tau870 = tau870 + tau(:,:,k)
        enddo
      enddo  ! nbins

      angstrom(i1:i2,j1:j2) = &
        -log(tau470(i1:i2,j1:j2)/tau870(i1:i2,j1:j2)) / &
         log(470./870.)
   endif
   deallocate(tau,ssa)
   __RETURN__(__SUCCESS__)
   end subroutine Aero_Compute_Diags
