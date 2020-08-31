! This little program will demonstrate the skeleton of what needs to be
! done to compute the parameters needed for the SW radiative transfer 
! calculations.
! Requires a chem.eta file

  program chem_miecalc

  use m_die, only: die
  use Chem_RegistryMod
  use Chem_BundleMod

! +++++++++++++++++++++++++++++++++++++++++++++++++++
! This will be mandatory for the SW band calculations
  use Chem_MieMod
! ---------------------------------------------------


  implicit none

  character(len=*), parameter :: myname = 'chem_mie'
  type(Chem_Registry) :: reg      ! chemistry registry
  type(Chem_Bundle)   :: w_c      ! chemistry bundle
  type(Chem_Bundle)   :: w_tau    ! tau chemistry bundle
  type(Chem_Bundle)   :: w_ssa    ! ssa chemistry bundle
  type(Chem_Bundle)   :: w_gasym  ! asymmetry parameter bundle
  type(Chem_Mie)      :: mie_tables
  real :: channel                 ! requested band number to get properties of
  real :: channels(8)             ! bands wanted from table
  integer :: i, j, k, im, jm, km, iq
  integer :: nymd, nhms, timidx, freq, rc, ier
  integer iarg, iargc, argc, lenfile
  real :: tau, ssa, gasym
  real, pointer :: rh(:,:,:)
  character(len=255) :: infile, outfile, filename, rcfile, argv
  character(len=255) :: which(5)

  data channels /1., 2., 3., 4., 5., 6., 7., 8./


! Parse the command line (see usage() below)
  argc = iargc()
  if(argc .lt. 1) call usage()
  iarg = 0
  outfile = 'test'
  rcfile  = 'AodBands_Registry.rc'
  do i = 0, 32767
   iarg = iarg+1
   if(iarg .gt. argc) exit
   call GetArg(iarg, argv)
   select case(argv)
    case ("-o")
     if(iarg+1 .gt. argc) call usage()
     iarg = iarg+1
     call GetArg(iarg, outfile)
    case ("-t")
     if(iarg+1 .gt. argc) call usage()
     iarg = iarg+1
     call GetArg(iarg, rcfile)
    case default
     infile = argv
   end select
  end do
  rcfile = trim(rcfile)
  infile = trim(infile)
  outfile = trim(outfile)
  lenfile = len(trim(outfile))

! Read the chemistry registry
! ---------------------------
  reg = Chem_RegistryCreate(ier)
  if(ier /= 0) call die(myname, 'cannot create registry')

! Read the chemistry bundle from the infile
! -------------------------------------------------
  call Chem_BundleRead(infile, nymd, nhms, w_c, rc, freq=freq, ChemReg=reg)

! Create the tau bundle
! -------------------------------------------------
  im = w_c%grid%im
  jm = w_c%grid%jm
  km = w_c%grid%km
  reg = Chem_RegistryCreate(rc,rcfile=rcfile)
  call Chem_BundleCreate(reg, im, jm, km, w_tau, ier)
  if(ier /= 0) call die(myname, 'cannot create tau bundle')
  call Chem_BundleCreate(reg, im, jm, km, w_ssa, ier)
  if(ier /= 0) call die(myname, 'cannot create ssa bundle')
  call Chem_BundleCreate(reg, im, jm, km, w_gasym, ier)
  if(ier /= 0) call die(myname, 'cannot create gasym bundle')

! Fill the bundles
! -------------------------------------------------
  w_tau%delp = w_c%delp
  w_tau%rh = w_c%rh
  w_ssa%delp = w_c%delp
  w_ssa%rh = w_c%rh
  w_gasym%delp = w_c%delp
  w_gasym%rh = w_c%rh

! +++++++++++++++++++++++++++++++++++++++++++++++++++
! Read the Mie Tables
! -------------------
  mie_tables = Chem_MieCreate(rcfile, rc)  

! Compute what you want
! ---------------------
! Need to clean up the functionality.
! You must pass the following things:
!  1) the mie_tables
!  2) tracer name
!  3) requested channel (band) by number (e.g., 1, 2, 3, ...)
!  4) tracer mass in gridbox, 
!     = mixing  ratio * delp / grav
!  5) relative humidity (scaled 0 - 1)
! You get out ier (error code) and optionally
! tau, ssa, and gasym of band
! Loop should be over space, constituent, and band
! Need to add functionality, but skeletally this is correct
! Only works now dust bin 1 and channel 1

  channel = 1.
  iq = w_c%reg%i_DU
  do k = 1, km
   do j = 1, jm
    do i = 1, im

!    Get the parameters from the mie tables 
!    Note that the RH from the Chem Bundle is 0 - 100; want 0 - 1 (fraction)
!    so divide by 100
     call Chem_MieQuery(mie_tables, 'du001', channel, &
                         w_c%q(i,j,k,iq)*w_c%delp(i,j,k)/9.81, &
                         w_c%rh(i,j,k), &
                         tau=tau, ssa=ssa, gasym=gasym,rc=ier)

!    Fill in the values
     w_tau%q(i,j,k,1)   = tau
     w_ssa%q(i,j,k,1)   = ssa
     w_gasym%q(i,j,k,1) = gasym

    enddo
   enddo
  enddo

! Destroy the Mie Tables
! ----------------------
  call Chem_MieDestroy(mie_tables, rc)  


! ---------------------------------------------------


! Write the AOD
! -------------------------------------------------
  filename = trim(outfile(1:lenfile)//'*.nc4')

  filename = trim(outfile(1:lenfile)//'.tau.nc4')
  call Chem_BundleWrite( filename, nymd, nhms, 0, w_tau, rc, &
                         verbose=.true.)
  filename = trim(outfile(1:lenfile)//'.ssa.nc4')
  call Chem_BundleWrite( filename, nymd, nhms, 0, w_ssa, rc, &
                         verbose=.true.)

  filename = trim(outfile(1:lenfile)//'.gasym.nc4')
  call Chem_BundleWrite( filename, nymd, nhms, 0, w_gasym, rc, &
                         verbose=.true.)

! ----------------------------------------------------------------------------
  contains

  subroutine usage()
  print *
  print *,'Usage: '
  print *,'  chem_mie.x [-o outfile -t rcfile ] infile'
  print *
  print *, 'where'
  print *
  print *, '-o outfile   output file containing AOD'
  print *, '-t rcfile    resource file specifying channels for AOD calc'
  print *, 'infile       mandatory c_rst file'
  print *
  call exit(1)
  end subroutine usage

end
