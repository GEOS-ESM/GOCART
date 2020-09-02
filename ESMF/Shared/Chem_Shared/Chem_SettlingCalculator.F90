! Colarco, January 2009
! Code to run offline the Chem_SettlingMod Module from GEOS-5
! to test how it works (if at all!).

! NOTE: There are filenames coded below that need to be
! interfaced to properly: requires an input file
! of the chemistry (this is specified on the input line)
! but also needs some ancillary meteorology which is
! currently hard-coded below.  Should be fixed.

  program Chem_SettleCalc

  use ESMF                  ! Needed to initialize MPI
  use Chem_Mod              ! Chemistry Base Class
  use Chem_StateMod         ! Chemistry State
  use Chem_SettlingMod      ! Settling
  use Chem_ConstMod, only: grav, von_karman, cpd     ! Constants !
  use Chem_UtilMod          ! I/O
  use m_inpak90             ! Resource file management
  use m_die, only: die

  implicit none

! Variables
  integer :: i, j, k, n, iarg, iargc, argc, fid
  integer, parameter :: READ_ONLY=1
  integer :: nymd, nhms, flag
  integer :: nbins, i1, i2, j1, j2, ig, jg, im, jm, km, ijl, ijkl
  integer, dimension(32) :: ier
  type(Chem_Bundle)   :: w_c      ! chemistry bundle used for RH and pressure
  type(Chem_Bundle)   :: w_vsettle ! bundle to hold output
  type(Chem_Registry) :: regInp
  character(len=255) :: argv, infile, rcfile, myname
  character(len=2)   :: XX
  real :: qmin, qmax, cdt
  real, pointer, dimension(:) :: radius, & ! particle effective radius [um]
                                 rhop      ! particle density [kg m-3]
  real, pointer, dimension(:,:,:) :: pm, &   ! layer midpoint air pressure [Pa]
                                     rhoa, & ! layer midpoint air density [kg m-3]
                                     t, &    ! layer midpoint temperature [K]
                                     rh, &   ! layer relative humidity
                                     delp, & ! layer pressure thickness [Pa]
                                     hght    ! layer mid-point geopotential height [m]
  real, pointer, dimension(:,:)   :: hsurf
  real, pointer, dimension(:,:,:) :: hghte   ! layer edge geopotential height [m] 
  type(Chem_Array), pointer, dimension(:)  :: fluxout, vsettle

!
  ier(:) = 0
  myname = 'Chem_SettlingCalc'

! This is needed to set up the MPI for the pmaxmin calls below and in
! Chem_SettlingMod
  call ESMF_Initialize( rc = ier(0))

! Parse the command line (see usage() below)
  argc = iargc()
  iarg = 0
  rcfile = ''
  do i = 0, 32767
   iarg = iarg+1
   if(iarg .gt. argc) exit
   call GetArg(iarg, argv)
   select case(argv)
    case ("-t")
     if(iarg+1 .gt. argc) call usage()
     iarg = iarg+1
     call GetArg(iarg, rcfile)
    case default
     infile = argv
   end select
  end do


! Provide a resource file (e.g., SS_GridComp.rc) which has the
! particle bins and sizes.  If no resource file is provided, 
! then default conditions are for 1 bin of 1 um effective radius
! and 1000 kg m-3 density.
! The value of "flag" is used to define what relative humidity
! correction to apply to the particles.  0: none, 1: Fitzgerald SS, 
! 2: Gerber SS
  flag = 0
  if(rcfile .eq. '') then
   nbins = 1
   allocate(radius(nbins), rhop(nbins), stat = ier(1))
   if(ier(1) /= 0) call die(myname,'cannot allocate radius, rhop')

   radius(1) = 1.e-6
   rhop(1)   = 1000.

  else
!  Read the resource file provided
!  Works only for SS or DU
   if(rcfile(1:2) .eq. 'DU' .or. &
      rcfile(1:2) .eq. 'SS') then
    XX = trim(rcfile(1:2))
   else
    call die(myname,'do not know what to do with rcfile '//rcfile)
   endif

   call i90_loadf ( rcfile, ier(1) )
   if(ier(1) /= 0) call die(myname,'cannot open rcfile '//rcfile)
   call i90_label ( 'number_'//XX//'_bins:', ier(1) )
   nbins = i90_gint ( ier(2) )
   if(any(ier(:) /= 0)) call die(myname,'cannot find number of bins in '//rcfile)

   allocate(radius(nbins), rhop(nbins), stat = ier(1))
   if(ier(1) /= 0) call die(myname,'cannot allocate radius, rhop')

!  Particle radius
!  ---------------
   call i90_label ( 'particle_radius:', ier(1) )
   do n = 1, nbins
      radius(n)    = i90_gfloat ( ier(n+1) )*1.e-6
   end do
   if ( any(ier(1:nbins+1) /= 0) ) call die(myname,'cannot get radii')

!  Dry Particle Density
!  ---------------
   call i90_label ( XX//'_density:', ier(1) )
   do n = 1, nbins
      rhop(n)      = i90_gfloat ( ier(n+1) )
   end do
   if ( any(ier(1:nbins+1) /= 0) ) call die(myname,'cannot get rhop')
!                          -------

!  If Seasalt, check for RH flag
!  -----------------------------
   if(XX .eq. 'SS') then
    call i90_label ( 'rhFlag:', ier(1) )
    flag = i90_gint(ier(2))
    if ( any(ier(1:nbins+1) /= 0) ) call die(myname,'cannot get rhflag')
   endif
  endif

  regInp = Chem_RegistryCreate(ier(1),'Chem_Registry.rc')
  call Chem_BundleRead(infile, nymd, nhms, w_c, ier(2), ChemReg=regInp)
  i1 = w_c%grid%i1
  i2 = w_c%grid%i2
  ig = w_c%grid%ig
  j1 = w_c%grid%j1
  j2 = w_c%grid%j2
  jg = w_c%grid%jg
  im = w_c%grid%im
  jm = w_c%grid%jm
  km = w_c%grid%km
  ijl  = ( i2 - i1 + 1 ) * ( j2 - j1 + 1 )
  ijkl = ijl * km
  allocate(pm(i1:i2,j1:j2,1:km), t(i1:i2,j1:j2,1:km), &
           rhoa(i1:i2,j1:j2,1:km), rh(i1:i2,j1:j2,1:km), &
           delp(i1:i2,j1:j2,1:km), hght(i1:i2,j1:j2,1:km), stat = ier(1))
  allocate(hsurf(i1:i2,j1:j2), stat=ier(2))
  allocate(hghte(i1:i2,j1:j2,1:km+1), stat = ier(3))
  if(any(ier(:) /= 0)) call die(myname,'cannot allocate atmosphere')

! Allocate space for output arrays
  allocate(fluxout(nbins), vsettle(nbins), stat=ier(1))
  do n = 1, nbins
   allocate(fluxout(n)%data2d(i1:i2,j1:j2), &
            vsettle(n)%data3d(i1:i2,j1:j2,1:km), stat=ier(n+1))
  enddo
  if( any(ier(1:nbins+1) /= 0)) call die(myname,'cannot allocate output arrays')

! Temperature, RH, DELP, HGHT
  call GFIO_Open('d520_fp.tavg3d_dyn_v.20090101_0000z.hdf', &
                 READ_ONLY, fid, ier(1))
  call GFIO_GetVar(fid,'T', nymd, nhms, im, jm, 1, km, t, ier(2))
  call GFIO_GetVar(fid,'RH', nymd, nhms, im, jm, 1, km, rh, ier(3))
  call GFIO_GetVar(fid,'DELP', nymd, nhms, im, jm, 1, km, delp, ier(4))
  call GFIO_GetVar(fid,'HGHT', nymd, nhms, im, jm, 1, km, hght, ier(5))
  call GFIO_Close(fid,ier(6))

! HGHTE
  call GFIO_Open('d520_fp.tavg3d_met_e.20090101_0000z.hdf', &
                 READ_ONLY, fid, ier(1))
  call GFIO_GetVar(fid,'HGHTE', nymd, nhms, im, jm, 1, km+1, hghte, ier(2))
  call GFIO_Close(fid,ier(3))

! Loop over vertical coordinate to get pressures
  pm(i1:i2,j1:j2,1) = w_c%grid%ptop + 0.5*delp(i1:i2,j1:j2,1)
  do k = 2, km
   pm(i1:i2,j1:j2,k) = pm(i1:i2,j1:j2,k-1) &
        +0.5*(delp(i1:i2,j1:j2,k)+delp(i1:i2,j1:j2,k-1))
  end do

! Air Density
  rhoa = pm/287./t

! Surface Height
  hsurf = hghte(:,:,km+1)

!  call pmaxmin('delp ', delp, qmin, qmax, ijkl, 1, 1.)
!  call pmaxmin('hght ', hght, qmin, qmax, ijkl, 1, 1.)
!  call pmaxmin('rh   ', rh, qmin, qmax, ijkl, 1, 1.)
!  call pmaxmin('t    ', t, qmin, qmax, ijkl, 1, 1.)
!  call pmaxmin('rhoa ', rhoa, qmin, qmax, ijkl, 1, 1.)

  cdt = 1800.

!  Create and initialize output chemistry bundle
   call Chem_BundleCreate(regInp, &
                          i1, i2, ig, im, &
                          j1, j2, jg, jm, km, &
                          w_vsettle, ier(1))

   if(ier(1) /= 0) call die(myname, 'cannot create bundle')
   w_vsettle%delp = w_c%delp
   w_vsettle%rh = w_c%rh
   do n = 1, regInp%nq
    w_vsettle%qa(n)%data3d = 0.0
   end do

  call Chem_Settling ( i1, i2, j1, j2, km, &
                       w_c%reg%i_SS, w_c%reg%j_SS, nbins, flag, &
                       radius, rhop, cdt, w_c, &
                       t, rhoa, hsurf, hghte, fluxout, ier(1), &
                       vsettle)

  do n = 1, nbins
   w_vsettle%qa(n)%data3d = vsettle(n)%data3d
  end do

  call Chem_BundleWrite('d520_fp.tavg3d_vsettle_v.20090101.hdf', &
                         nymd, nhms, 0, w_vsettle, ier(1), verbose=.true., new=.true.)

  do n = 1, nbins
   deallocate(fluxout(n)%data2d, vsettle(n)%data3d, stat=ier(n))
  enddo
  deallocate(fluxout, vsettle, stat=ier(nbins+1))
  if( any(ier(1:nbins+1) /= 0)) call die(myname,'cannot deallocate output arrays')

  deallocate(radius, rhop, stat=ier(1))
  if(ier(1) /= 0) call die(myname,'cannot deallocate radius, rhop')
  deallocate(pm, rhoa, t, rh, delp, hght, hghte, hsurf, stat=ier(1))
  if(ier(1) /= 0) call die(myname,'cannot deallocate atmosphere')



! ----------------------------------------------------------------------------
  contains

  subroutine usage()
  print *
  print *,'Usage: '
  print *,'  Chem_SettleCalc.x [-t rcfile ] infile'
  print *
  print *, 'where'
  print *
  print *, '-t rcfile    resource file specifying particle radii, density'
  print *, '-v           request verbose output'
  print *, 'infile       mandatory input aer_v file'
  print *
  call exit(1)
  end subroutine usage

end program Chem_SettleCalc
