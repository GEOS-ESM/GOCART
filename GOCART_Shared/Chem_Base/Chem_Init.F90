! Create an initialization for the chemistry bundle
! - get the chemistry registry
! - create the chemistry bundle
! !REVISION HISTORY:
!   31Mar2005 - Todling - iuic was not used before defined

  program chem_initialize

  use m_ioutil, only: luavail
  use m_die, only: die
  use Chem_InitMod
  use Chem_RegistryMod
  use Chem_BundleMod

  implicit none

! Hard-wired c55 grid parameters
  integer, parameter :: im = 288, jm = 181, km = 55
  real :: ak(km+1), bk(km+1)

  real, parameter :: grav = 9.81, RAIR = 287., CP = 1004., PREF = 1.e5
  character(len=*), parameter :: myname = 'chem_ini'
  type(Chem_Registry) :: reg      ! chemistry registry
  type(Chem_Bundle)   :: w_c      ! chemistry bundle
  integer :: i,j,k,n,nq,prec      ! working variables
  integer :: iuic, nstep, nymd, nhms
  integer :: ier                  ! error checking
  real :: ps(im,jm), pc(im,jm,km)                   ! Pressures
  real :: zc(im,jm,km), dz(im,jm,km)                ! vertical coordinates
  real :: u(im,jm,km), v(im,jm,km)                  ! velocity
  real :: pt(im,jm,km)                              ! potential temperature  
  real, parameter :: ptop = 1.                      ! pressure at top of atm.
  real, pointer :: delp(:,:,:), q(:,:,:,:)
  integer iarg, iargc, argc
  character(len=255) :: chemfile, dynfile, argv, dummystr

  data ak / 1.00000,       2.00000,       3.27000,                     &
            4.75850,       6.60000,       8.93450,                     &
           11.97030,      15.94950,      21.13490,                     &
           27.85260,      36.50410,      47.58060,                     &
           61.67790,      79.51340,     101.94420,                     &
          130.05080,     165.07920,     208.49720,                     &
          262.02120,     327.64330,     407.65670,                     &
          504.68050,     621.68000,     761.98390,                     &
          929.29430,    1127.68880,    1364.33920,                     &
         1645.70720,    1979.15540,    2373.03610,                     &
         2836.78160,    3380.99550,    4017.54170,                     &
         4764.39320,    5638.79380,    6660.33770,                     &
         7851.22980,    9236.56610,   10866.34270,                     &
        12783.70000,   15039.30000,   17693.00000,                     &
        20119.20876,   21686.49129,   22436.28749,                     &
        22388.46844,   21541.75227,   19873.78342,                     &
        17340.31831,   13874.44006,   10167.16551,                     &
         6609.84274,    3546.59643,    1270.49390,                     &
            0.00000,       0.00000   /

  data bk  /0.00000,       0.00000,       0.00000,           &
            0.00000,       0.00000,       0.00000,           & 
            0.00000,       0.00000,       0.00000,           &
            0.00000,       0.00000,       0.00000,           &
            0.00000,       0.00000,       0.00000,           &
            0.00000,       0.00000,       0.00000,           &
            0.00000,       0.00000,       0.00000,           &
            0.00000,       0.00000,       0.00000,           &
            0.00000,       0.00000,       0.00000,           &
            0.00000,       0.00000,       0.00000,           &
            0.00000,       0.00000,       0.00000,           &
            0.00000,       0.00000,       0.00000,           &
            0.00000,       0.00000,       0.00000,           &
            0.00000,       0.00000,       0.00000,           &
            0.00696,       0.02801,       0.06372,           &
            0.11503,       0.18330,       0.27033,           &
            0.37844,       0.51046,       0.64271,           &
            0.76492,       0.86783,       0.94329,           &
            0.98511,       1.00000  /


! Parse the command line (see usage() below)
  argc = iargc()
  if(argc .lt. 1) call usage()
  iarg = 0
  chemfile = 'chem_ini.c_rst'
  do i = 0, 32767
   iarg = iarg+1
   if(iarg .gt. argc) exit
   call GetArg(iarg, argv)
   select case(argv)
    case ("-o")
     if(iarg+1 .gt. argc) call usage()
     iarg = iarg+1
     call GetArg(iarg, chemfile)
    case ("-grid")
     if(iarg+1 .gt. argc) call usage()
     iarg = iarg+1
     call GetArg(iarg, dummystr)
    case default
     dynfile = argv
   end select
  end do
  if ( trim(chemfile) .eq. 'chem_ini.c_rst' ) then 
   i = index ( dynfile, 'd_rst' )
   if ( i .gt. 1 ) then
    chemfile = dynfile(1:i-1) // 'c_rst' // dynfile(i+5:)
   end if
  end if


! Read the chemistry registry
! ---------------------------
  reg = Chem_RegistryCreate(ier)
  if(ier /= 0) call die(myname, 'cannot create registry')
  nq = reg%nq


! Create the chemistry bundle to be filled in later
! -------------------------------------------------
  call Chem_BundleCreate1PE_(reg, im, jm, km, w_c, ier)
  if(ier /= 0) call die(myname, 'cannot create bundle')

  delp => w_c%delp
  q => w_c%q


! Read in some dynamics fields from a restart file in order to
! fill in delp and to calculate vertical coordinates in a useful
! fashion.
! --------------------------------------------------------------
  iuic = luavail()
  open(iuic,file=dynfile,  &
            form='unformatted',access='sequential')
  read(iuic) nstep,nymd,nhms
  read(iuic) ps, delp, u, v, pt
  do n = 1, 2
   read(iuic,END=9999) (((q(i,j,k,n),i=1,im),j=1,jm),k=1,km)
  end do
9999 continue
  close(iuic)


! Calculate the pressures and vertical grid coordinates from the
! hydrostatic equation.  Use ak and bk to calculate grid mid-point
! pressures pc [Pa].  Use the potential temperature read in along with the
! pressures to calculate level thickness dz [m] with hydrostatic 
! equation.  Integrate level thickness to get mid-point altitudes zc [m].
  k = km
  do j = 1, jm
   do i = 1, im
    pc(i,j,k) = 0.5*(ak(k)+ak(k+1)+(bk(k)+bk(k+1))*ps(i,j))
    dz(i,j,k) = RAIR*(pt(i,j,k)+273.)*delp(i,j,k)/GRAV/pc(i,j,k) &
                    *(PREF/pc(i,j,k))**(-RAIR/CP)
    zc(i,j,k) = 0.5*dz(i,j,k)
   end do
  end do

  do k = km-1, 1, -1
   do j = 1, jm
    do i = 1, im
     pc(i,j,k) = 0.5*(ak(k)+ak(k+1)+(bk(k)+bk(k+1))*ps(i,j))
     dz(i,j,k) = RAIR*(pt(i,j,k)+273.)*delp(i,j,k)/GRAV/pc(i,j,k) &
                     *(PREF/pc(i,j,k))**(-RAIR/CP)
     zc(i,j,k) = zc(i,j,k+1) + 0.5*(dz(i,j,k+1)+dz(i,j,k))
    end do
   end do
  end do


! Fill the bundle with initial conditions
  if(reg%doing_CO) call initial_q_('CO',reg%i_CO, reg%j_CO)
  if(reg%doing_DU) call initial_q_('DU',reg%i_DU, reg%j_DU)
  if(reg%doing_SS) call initial_q_('SS',reg%i_SS, reg%j_SS)
  if(reg%doing_SU) call initial_q_('SU',reg%i_SU, reg%j_SU)
  if(reg%doing_BC) call initial_q_('BC',reg%i_BC, reg%j_BC)
  if(reg%doing_OC) call initial_q_('OC',reg%i_OC, reg%j_OC)
  if(reg%doing_XX) call initial_q_('XX',reg%i_XX, reg%j_XX)

! Stat the bundle
  call Chem_BundleStat(6, w_c, ier)
  if(ier /= 0) call die(myname, 'cannot stat bundle')


! Write the bundle
  prec = 0
  call Chem_BundleWrite(chemfile, nymd, nhms, prec, w_c, ier, verbose=.true.)
  if(ier /=0 ) call die(myname,'cannot write bundle')
 

! ----------------------------------------------------------------------------
  contains

  subroutine initial_q_(name,n1,n2)

  character(len=*) :: name
  type(Chem_Init)     :: init     ! initialization registry
  integer :: n1, n2, ier
  real :: lat, lat_min, lat_del, lon, lon_min, lon_del
  real, parameter :: REARTH=6370000., GRAV = 9.81
  real :: pi, sqrttwopi, mass, delm, area
  real :: x, y, ang_x, ang_y, angarg
  real :: amp, lat0, lon0, z0, rx, ry, rz

  pi = 4.*atan(1.)
  mass = 0.

  init = Chem_InitResource(name, ier)
  if(ier /= 0) call die(name, 'cannot read init resource')

  amp = init%amp
  lon0 = 2.*pi*init%lon0/360.
  lat0 = 2.*pi*init%lat0/360.
  z0 = init%z0 * 1000.
  rx = init%rx * 1000.
  ry = init%ry * 1000.
  rz = init%rz * 1000.

! for a positive amplitude I assume that that is a total amount of mass
! in kg I want distributed in an ellipse (assumes consituent is in mass
! mixing ratio units).  Integrate ellipse for total mass of air.  Divide
! constituent mass by air mass to arrive at uniform mass mixing ratio.
  if(amp .gt. 0.) then
!  first find the grid cells in the ellipse and integrate to get mass
!  of air inside
   do k = 1, km
    do j = 1, jm
     lat = 2.*pi*(w_c%grid%lat_min + (j-1)*w_c%grid%lat_del)/360.
     area =  2.*pi*REARTH*cos(lat)/360.*w_c%grid%lon_del &
           * 2.*pi*REARTH/360.*w_c%grid%lat_del
     do i = 1, im
      lon = 2.*pi*(w_c%grid%lon_min+(i-1)*w_c%grid%lon_del)/360.
!     Use angarg to check for round-off error
      angarg = sin(lat0)**2.+cos(lat0)**2.*cos(lon0-lon)
      if(angarg .ge. 1.) angarg = 1.d0
      ang_x = acos(angarg)
      angarg = sin(lat)*sin(lat0)+cos(lat)*cos(lat0)
      if(angarg .ge. 1.) angarg = 1.d0
      ang_y = acos(angarg)
      x = REARTH*ang_x
      y = REARTH*ang_y
      if( (x/rx)**2.+(y/ry)**2.+((zc(i,j,k)-z0)/rz)**2. .lt. 1) then  
       mass = mass+delp(i,j,k)*area/GRAV
      endif
     end do
    end do
   end do

   delm = amp/mass

   do k = 1, km
    do j = 1, jm
     lat = 2.*pi*(w_c%grid%lat_min + (j-1)*w_c%grid%lat_del)/360.
     do i = 1, im
      lon = 2.*pi*(w_c%grid%lon_min+(i-1)*w_c%grid%lon_del)/360.
!     Use angarg to check for round-off error
      angarg = sin(lat0)**2.+cos(lat0)**2.*cos(lon0-lon)
      if(angarg .ge. 1.) angarg = 1.d0
      ang_x = acos(angarg)
      angarg = sin(lat)*sin(lat0)+cos(lat)*cos(lat0)
      if(angarg .ge. 1.) angarg = 1.d0
      ang_y = acos(angarg)
      x = REARTH*ang_x
      y = REARTH*ang_y
      if( (x/rx)**2.+(y/ry)**2.+((zc(i,j,k)-z0)/rz)**2. .lt. 1) then  
       do n = n1, n2
        q(i,j,k,n) = delm
       end do
      endif
     end do
    end do
   end do

  endif   


!!! THIS PART IS NOT QUITE CORRECT

! for a negative amplitude I assume that that the absolute value of the 
! amplitude is a desired maximum mixing ratio for the species and that
! it will be distributed in space in some sort of gaussian distribution
! of lat, lon, and z.
  if(amp .le. 0.) then 
   sqrttwopi = sqrt(2.*pi)
   do k = 1, km
    do j = 1, jm
     lat = w_c%grid%lat_min + (j-1)*w_c%grid%lat_del
     do i = 1, im
      lon = w_c%grid%lon_min+(i-1)*w_c%grid%lon_del
      do n = n1, n2
       q(i,j,k,n) = -amp*exp(-(1./rx**2.)*(lon-lon0)**2.) &
                        *exp(-(1./ry**2.)*(lat-lat0)**2.) &
                        *exp(-(1./rz**2.)*(zc(i,j,k)-z0)**2.)
      end do
     end do
    end do
   end do
  endif

  end subroutine initial_q_


! ---------------------------------------------------------------------------
  subroutine usage()
  print *
  print *,'Usage: '
  print *,'  chem_init.x [-o chemfile -grid ] dynfile'
  print *
  print *, 'where'
  print *
  print *, '-o chemfile   output chemistry initialization'
  print *, '              file. Default: same as dynfile'
  print *, '              with substring "d_rst" replaced'
  print *, '              with "c_rst"'
  print *, '-grid grtype  a55, b55, c55, etc. (inactive now!)'
  print *, 'dynfile       mandatory d_rst file'
  print *
  call exit(1)
  end subroutine usage

end
