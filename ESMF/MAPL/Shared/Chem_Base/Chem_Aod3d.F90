! Colarco, June 2006
! Offline calculator of 3D aerosol optical properties
! 1) Read a registry determining aerosol optical properties to compute
! 2) Read a chem bundle of aerosol distribution from a file
! 3) Read an output registry naming the output chem bundle
! 4) Compute

  program chem_aodcalc

  use m_die, only: die
  use Chem_MieMod
  use Chem_RegistryMod
  use Chem_BundleMod
  use m_fpe, only: isnan

  implicit none

  real, parameter :: grav = 9.80616
  real, parameter :: rgas = 287.04

  character(len=*), parameter :: myname = 'chem_aod'
  type(Chem_Mie)      :: mie_tables
  type(Chem_Registry) :: regInp      ! chemistry registry
  type(Chem_Registry) :: regOut      ! chemistry registry
  type(Chem_Bundle)   :: w_c         ! input aerosol chemistry bundle 
  type(Chem_Bundle)   :: w_tau       ! output total tau chemistry bundle
  type(Chem_Bundle)   :: w_taudu     ! output du tau chemistry bundle
  type(Chem_Bundle)   :: w_tauant    ! output anthro tau chemistry bundle
  type(Chem_Bundle)   :: w_tauss     ! output seasalt tau chemistry bundle
  type(Chem_Bundle)   :: w_tauoc     ! output organic carbon tau chemistry bundle
  type(Chem_Bundle)   :: w_taubc     ! output black carbon tau chemistry bundle
  type(Chem_Bundle)   :: w_taucc     ! output total carbon tau chemistry bundle
  type(Chem_Bundle)   :: w_tausu     ! output sulfate tau chemistry bundle
  integer :: i, j, k, im, jm, km, idx, n, ndx
  integer :: i1, i2, ig, j1, j2, jg, ik, iq, kk
  integer :: nymd, nhms, freq, rc, ier, fid
  integer :: idxTable
  integer :: iarg, iargc, argc, lenfile
  logical :: doing_dust, doing_anthro, doing_ss, doing_oc, doing_bc, doing_su, doing_cc
  logical :: doing_geos4
  logical :: doing_dry
  logical :: new
  logical :: found_airdensfile, found_hghtefile
  real    :: channel, tau_, ssa_, bbck_, bext_, taulev, gasym_, scaleRH, maxRH, &
             qMass, p11_, p22_, vol_, area_, refr_, refi_
  integer :: itau, iext, issa, immr, ibbck, &
             iabck0, iabck1, ietob, igasym, idepol, &
             ivol, iarea, ireff, irefr, irefi
  integer, parameter :: READ_ONLY = 1
  real, pointer :: delz(:,:,:), t(:,:,:), q(:,:,:), hghte(:,:,:)
  real, pointer :: airdens(:,:,:) => null()
  character(len=255) :: infile, outfile, filename, airdensfile, hghtefile, rcfile, argv
  character(len=8)   :: datestr
  character(len=11)  :: chstr
  character(len=3)   :: chnstr

! Parse the command line (see usage() below)
  argc = iargc()
  if(argc .lt. 1) call usage()
  iarg = 0
  outfile = 'tau3d.nc4'
  rcfile  = 'Aod_CALIPSO.rc'
  doing_dust   = .false.
  doing_ss     = .false.
  doing_bc     = .false.
  doing_oc     = .false.
  doing_cc     = .false.
  doing_su     = .false.
  doing_anthro = .false.
  doing_geos4  = .false.
  doing_dry    = .false.
  found_airdensfile = .false.
  found_hghtefile = .false.
  do i = 0, 32767
   iarg = iarg+1
   if(iarg .gt. argc) exit
   call GetArg(iarg, argv)
   select case(argv)
    case ("-geos4")
     doing_geos4 = .true.
    case ("-dryaer")
     doing_dry = .true.
    case ("-dust")
     doing_dust = .true.
    case ("-ss")
     doing_ss = .true.
    case ("-su")
     doing_su = .true.
    case ("-bc")
     doing_bc = .true.
    case ("-oc")
     doing_oc = .true.
    case ("-cc")
     doing_cc = .true.
    case ("-anthro")
     doing_anthro = .true.
    case ("-airdensfile")
     if(iarg+1 .gt. argc) call usage()
     iarg = iarg+1
     call GetArg(iarg, airdensfile)
     found_airdensfile = .true.
    case ("-hghtefile")
     if(iarg+1 .gt. argc) call usage()
     iarg = iarg+1
     call GetArg(iarg, hghtefile)
     found_hghtefile = .true.
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
  if(found_hghtefile) hghtefile = trim(hghtefile)
  if(found_airdensfile) then
   airdensfile = trim(airdensfile)
  else
   airdensfile = infile
  endif

! Scaling of Relative Humidity
! Input optics files (e.g., optics_XX.nc4) have a fractional RH
! coordinate (that is, RH varies 0 - 1 in the file, as in GEOS-5).
! In GEOS-4, however, RH in the chem.eta file is represented as
! a percentage, so it varies 0 - 100%.  For compatibility we
! introduce a flag "-geos4" on execution.  If present, the RH
! from the chem.eta file is divided by 100 on input to the Mie
! calculator.
! If you requested to do the calculation like the aerosols were
! dry, then we set the RH like it is 0%
  scaleRH = 1.
  if(doing_geos4) scaleRH = 1. / 100.
  if(doing_dry)   scaleRH = 0.

! Hardwired: Create the output bundle
! -----------------------------------
  regOut = Chem_RegistryCreate(ier, rcfile)
  if(ier /= 0) call die(myname, 'cannot create output registry')

! Check the output registry for desired output fields
! ---------------------------------------------------
  itau   = -1
  iext   = -1
  issa   = -1
  immr   = -1
  ietob  = -1
  ibbck  = -1
  iabck0 = -1
  iabck1 = -1
  igasym = -1
  idepol = -1
  ivol   = -1
  iarea  = -1
  ireff  = -1
  irefr  = -1
  irefi  = -1
  do iq = 1, regOut%nq
   if(trim(regOut%vname(iq)) .eq. 'tau')        itau   = iq
   if(trim(regOut%vname(iq)) .eq. 'extinction') iext   = iq
   if(trim(regOut%vname(iq)) .eq. 'ssa')        issa   = iq
   if(trim(regOut%vname(iq)) .eq. 'mmr')        immr   = iq
   if(trim(regOut%vname(iq)) .eq. 'ext2back')   ietob  = iq
   if(trim(regOut%vname(iq)) .eq. 'attback0')   iabck0 = iq
   if(trim(regOut%vname(iq)) .eq. 'aback_sfc')  iabck0 = iq ! alias
   if(trim(regOut%vname(iq)) .eq. 'attback1')   iabck1 = iq
   if(trim(regOut%vname(iq)) .eq. 'aback_toa')  iabck1 = iq ! alias
   if(trim(regOut%vname(iq)) .eq. 'backscat')   ibbck  = iq
   if(trim(regOut%vname(iq)) .eq. 'gasym')      igasym = iq
   if(trim(regOut%vname(iq)) .eq. 'depol')      idepol = iq
   if(trim(regOut%vname(iq)) .eq. 'depol')      idepol = iq
   if(trim(regOut%vname(iq)) .eq. 'volume')     ivol   = iq
   if(trim(regOut%vname(iq)) .eq. 'area')       iarea  = iq
   if(trim(regOut%vname(iq)) .eq. 'reff')       ireff  = iq
   if(trim(regOut%vname(iq)) .eq. 'refreal')    irefr  = iq
   if(trim(regOut%vname(iq)) .eq. 'refimag')    irefi  = iq
  enddo

! At this point should do some checking that certain fields are
! provided.
  if(issa .gt. 0   .and. itau .lt. 0)                      ier = 1
  if(igasym .gt. 0 .and. (itau .lt. 0 .or. issa .lt. 0))   ier = 2
  if(ibbck .gt. 0  .and. itau .lt. 0)                      ier = 3
  if(ietob .gt. 0  .and. (itau .lt. 0 .or. ibbck .lt. 0))  ier = 4
  if((iabck0 .gt. 0 .or. iabck1 .gt. 0) .and. &
     (itau .lt. 0 .or.ibbck .lt. 0) )                      ier = 5
  if(idepol .gt. 0 .and. (itau .lt. 0 .or. issa .lt. 0))   ier = 6
  if(ireff .gt. 0 .and. (ivol .lt. 0 .or. iarea .lt. 0))   ier = 7
  if(irefr .gt. 0 .and. ivol .lt. 0)                       ier = 8
  if(irefi .gt. 0 .and. ivol .lt. 0)                       ier = 8

  if(ier /= 0) then
     print *,'----------------------------------------------------'
     print *,'ier = ', ier
     print *,'issa, itau, igasym, ibbck, ietob, iabck0, iabck1, idepol, ireff, irefr, irefi = ', &
            issa, itau, igasym, ibbck, ietob, iabck0, iabck1, idepol, ireff, irefr, irefi
     call die(myname,'inconsistency in output registry')
  end if

! Hardwired: Read the input chemistry registry from Chem_Registry.rc
! This registry file describes the input chemistry bundle being 
! operated on.
! ------------------------------------------------------------------
  regInp = Chem_RegistryCreate(ier,'Chem_MieRegistry.rc')
  if(ier /= 0) call die(myname, 'cannot create registry')

! Hardwired: use the Chem_MieMod function to create the Mie tables
! ------------------------------------------------------------------
  mie_tables = Chem_MieCreate(rcfile,ier)

! Call subroutine to check the number of times and frequency on
! input file; output will have same information
! ------------------------------------------------------------------
  call check_infile  ! returns ndx and freq

  new = .true.
  do idx = 1, ndx

!  Read the input chemistry bundle from the infile
!  -------------------------------------------------
   call Chem_BundleRead(infile, nymd, nhms, w_c, rc, freq=freq, &
                        ChemReg=regInp, timidx=idx)

   i1 = w_c%grid%i1
   i2 = w_c%grid%i2
   ig = w_c%grid%ig
   j1 = w_c%grid%j1
   j2 = w_c%grid%j2
   jg = w_c%grid%jg
   im = w_c%grid%im
   jm = w_c%grid%jm
   km = w_c%grid%km

!  Layer thickness in cartesian sense is needed for certain extrinsic
!  properties, specifically any calculations requiring backscatter
!  or extinction.  If these fields are requested then we need to allocate
!  space for the layer thickness ("delz").  The value of delz is
!  filled in either by (1) reading the layer edge heights from a file (the
!  "hghtefile") or (2) by reading the air density from a file (the
!  "airdensfile") and computing via hydrostatic relation 
!  (i.e., dz = dp/g/rhoa) or (3) by computing directly from provided 
!  thermodynamic variables (e.g., t, qv).
   if(iext .gt. 0 .or. ibbck .gt. 0) then 

!   Allocate space to contain the thickness of the layer
    allocate(delz(i1:i2,j1:j2,1:km), stat=ier)
    if(ier /= 0) call die(myname,'could not allocate space for delz')

!   If "hghtefile" is provided let's get delz that way
    if(found_hghtefile) then
     allocate(hghte(i1:i2,j1:j2,1:km+1), stat=ier)
     if(ier /= 0) call die(myname,'could not allocate space for HGHTE')
     call GFIO_Open ( hghtefile, READ_ONLY, fid, ier )
     if(ier /= 0) call die(myname,'could not open '//hghtefile//' to read HGHTE')
     call GFIO_GetVar ( fid, 'HGHTE', nymd, nhms,     &
                        im, jm, 1, km+1, hghte, ier )
     if(ier /= 0) call die(myname,'could not read HGHTE from '//hghtefile)
     call GFIO_Close ( fid, ier)
     if(ier /= 0) call die(myname,'could not close infile for HGHTE read')
     do ik = 1, km
      delz(:,:,ik) = hghte(:,:,ik)-hghte(:,:,ik+1)
     end do
     deallocate(hghte,stat=ier)
     if(ier /= 0) call die(myname,'could not deallocate space for HGHTE')

    else
!    Get the air density, either from aerosol file or externally
     if(trim(airdensfile) .eq. trim(infile)) then
      airdens => delz
!     Get the air density from the aerosol file
      call GFIO_Open ( airdensfile, READ_ONLY, fid, ier )
      if(ier /= 0) call die(myname,'could not open '//airdensfile//' to read airdens')
      call GFIO_GetVar ( fid, 'AIRDENS', nymd, nhms,     &
                         im, jm, 1, km, airdens, ier )
      if(ier /= 0) then
         call GFIO_GetVar ( fid, 'airdens', nymd, nhms,     &
                            im, jm, 1, km, airdens, ier )
         if(ier /= 0) call die(myname,'could not read airdens from '//airdensfile)
      end if
      call GFIO_Close ( fid, ier)
      if(ier /= 0) call die(myname,'could not close infile for airdens read')
     else
!     Get air density from an external file
      call GFIO_Open ( airdensfile, READ_ONLY, fid, ier )
      if(ier /= 0) call die(myname,'could not open '//airdensfile//' to read airdens')
      call GFIO_GetVar ( fid, 'AIRDENS', nymd, nhms,     &
                         im, jm, 1, km, delz, ier )
      if ( ier /= 0 ) then
            call GFIO_GetVar ( fid, 'airdens', nymd, nhms,     &
                               im, jm, 1, km, delz, ier )
      end if
      if(ier .eq. 0) then
       airdens => delz
       call GFIO_Close ( fid, ier)
       if(ier /= 0) call die(myname,'could not close infile for airdens read')
      else
!      If the air density is not on the file try to compute delz directly
!      from the thermodynamic variables (need T, QV)
       allocate(t(i1:i2,j1:j2,1:km), q(i1:i2,j1:j2,1:km), stat=ier)
       if(ier /= 0) call die(myname,'could not allocate space to compute airdens')
       call GFIO_GetVar ( fid, 'T', nymd, nhms,     &
                          im, jm, 1, km, t, ier )
       if(ier /= 0) call die(myname,'could not read temperature '//airdensfile)
       call GFIO_GetVar ( fid, 'QV', nymd, nhms,     &
                          im, jm, 1, km, q, ier )
       if(ier /= 0) call die(myname,'could not read specific humidity from '//airdensfile)
       call GFIO_Close ( fid, ier)
       if(ier /= 0) call die(myname,'could not close infile for airdens read')
       do j = j1, j2
        do i = i1, i2
         call delz_
        enddo
       enddo
       deallocate(t,q,stat=ier)
       if(ier /= 0) call die(myname,'cloud not deallocate t,q')
      endif
     endif

!    If you actually read air density, compute dz
     if(associated(airdens)) then
      delz = w_c%delp/grav/airdens
     endif
    endif

!   Turn delz from m to km
    delz = delz/1000.

   endif

   print *, 'Computing AOD for ', nymd, nhms

!  Check that the RH seems sane
   if(.not. doing_geos4 .and. .not. doing_dry) then
    maxrh = maxval(w_c%rh)
    if(maxrh .gt. 2.) then
     print *, 'Maximum RH value = ', maxRH
     print *, 'Should you have chosen "-geos4" as a command line option?'
    endif
   endif

!  ==================================================================================
!  The enclosed bundle of code selects on what calculation we run and what is written
   do ik = 1, mie_tables%nch

    channel = mie_tables%channels(ik)

!   Create and initialize output chemistry bundle
    call create_species_bundle(w_tau)

!   If doing species
    if(doing_dust) call create_species_bundle(w_taudu)
    if(doing_anthro) call create_species_bundle(w_tauant)
    if(doing_ss) call create_species_bundle(w_tauss)
    if(doing_su) call create_species_bundle(w_tausu)
    if(doing_oc) call create_species_bundle(w_tauoc)
    if(doing_bc) call create_species_bundle(w_taubc)
    if(doing_cc) call create_species_bundle(w_taucc)


    do iq = 1, mie_tables%nq

!    Sanity check: Does the mie_table name = chem_registry name?
     if(trim(mie_tables%vname(iq)) .ne. trim(w_c%reg%vname(iq)) ) &
       call die(myname, 'mie_tables and chem_registry vname mismatch')
!    Cases
     idxTable = Chem_MieQueryIdx(mie_tables,mie_tables%vname(iq),rc)

     if(idxTable .ne. -1) then
      do k = 1, km
      do j = 1, jm
      do i = 1, im

!     fix in case NaN on input of qa
#ifndef sysAIX
      if(isnan(w_c%qa(iq)%data3d(i,j,k))) &
        w_c%qa(iq)%data3d(i,j,k) = tiny(w_c%qa(iq)%data3d(i,j,k))
#endif
      qMass = w_c%qa(iq)%data3d(i,j,k)*w_c%delp(i,j,k)/grav
      call Chem_MieQuery(mie_tables, idxTable, 1.*ik, &
                         qMass, &
                         w_c%rh(i,j,k) * scaleRH, tau=tau_, ssa=ssa_, &
                         bbck=bbck_, bext=bext_, gasym=gasym_, &
                         p11=p11_, p22=p22_, vol=vol_, area=area_, &
                         refr=refr_, refi=refi_)

!     Fill in the total values
!     Note the weighting of the ssa, backscatter, and e_to_b ratio
      call fill(w_tau)

      if(doing_dust .and. (iq .ge. w_c%reg%i_du .and. iq .le. w_c%reg%j_du)) then
       call fill(w_taudu)
      endif

      if(doing_su .and. (iq .ge. w_c%reg%i_su .and. iq .le. w_c%reg%j_su)) then
       call fill(w_tausu)
      endif

      if(doing_bc .and. (iq .ge. w_c%reg%i_bc .and. iq .le. w_c%reg%j_bc)) then
       call fill(w_taubc)
      endif

      if(doing_oc .and. (iq .ge. w_c%reg%i_oc .and. iq .le. w_c%reg%j_oc)) then
       call fill(w_tauoc)
      endif

      if(doing_ss .and. (iq .ge. w_c%reg%i_ss .and. iq .le. w_c%reg%j_ss)) then
       call fill(w_tauss)
      endif

      if(doing_anthro .and. &
         ( (iq .ge. w_c%reg%i_oc .and. iq .le. w_c%reg%j_oc) .or. &
           (iq .ge. w_c%reg%i_bc .and. iq .le. w_c%reg%j_bc) .or. &
           (iq .ge. w_c%reg%i_su .and. iq .le. w_c%reg%j_su) )      ) then
       call fill(w_tauant)
      endif

      if(doing_cc .and. &
         ( (iq .ge. w_c%reg%i_oc .and. iq .le. w_c%reg%j_oc) .or. &
           (iq .ge. w_c%reg%i_bc .and. iq .le. w_c%reg%j_bc) )      ) then
       call fill(w_taucc)
      endif

      enddo  ! i
      enddo  ! j
      enddo  ! k
     endif   
    enddo    ! iq

!    call Chem_RegistryPrint(regout)

!   Normalize the ssa calculation
    call normal(w_tau)
    if(doing_dust)   call normal(w_taudu)
    if(doing_ss)     call normal(w_tauss)
    if(doing_su)     call normal(w_tausu)
    if(doing_bc)     call normal(w_taubc)
    if(doing_oc)     call normal(w_tauoc)
    if(doing_cc)     call normal(w_taucc)
    if(doing_anthro) call normal(w_tauant)

!   Compute the attentuated backscatter terms
    call backscatter(w_tau)
    if(doing_dust) call backscatter(w_taudu)
    if(doing_ss)   call backscatter(w_tauss)
    if(doing_su)   call backscatter(w_tausu)
    if(doing_oc)   call backscatter(w_tauoc)
    if(doing_bc)   call backscatter(w_taubc)
    if(doing_cc)   call backscatter(w_taucc)
    if(doing_anthro) call backscatter(w_tauant)

!   Write the Chem_Bundle out
    write(datestr,'(i8.8)') nymd
    chstr = ''
    write(chnstr,'(i3.3)') ik
    if(mie_tables%nch > 1) chstr = '.channel'//chnstr
    filename = trim(outfile(1:lenfile)//trim(chstr))
    w_tau%rh = w_c%rh
    call Chem_BundleWrite( filename, nymd, nhms, 0, w_tau, rc, &
                           verbose=.true., new=new, freq=freq)
    if(doing_dust) then 
     filename = trim(outfile(1:lenfile)//trim(chstr)//'.dust')
     w_taudu%rh = w_c%rh
     call Chem_BundleWrite( filename, nymd, nhms, 0, w_taudu, rc, &
                            verbose=.true., new=new, freq=freq)
    endif
    if(doing_ss) then 
     filename = trim(outfile(1:lenfile)//trim(chstr)//'.ss')
     w_tauss%rh = w_c%rh
     call Chem_BundleWrite( filename, nymd, nhms, 0, w_tauss, rc, &
                            verbose=.true., new=new, freq=freq)
    endif
    if(doing_su) then 
     filename = trim(outfile(1:lenfile)//trim(chstr)//'.su')
     w_tausu%rh = w_c%rh
     call Chem_BundleWrite( filename, nymd, nhms, 0, w_tausu, rc, &
                            verbose=.true., new=new, freq=freq)
    endif
    if(doing_oc) then 
     filename = trim(outfile(1:lenfile)//trim(chstr)//'.oc')
     w_tauoc%rh = w_c%rh
     call Chem_BundleWrite( filename, nymd, nhms, 0, w_tauoc, rc, &
                            verbose=.true., new=new, freq=freq)
    endif
    if(doing_bc) then 
     filename = trim(outfile(1:lenfile)//trim(chstr)//'.bc')
     w_taubc%rh = w_c%rh
     call Chem_BundleWrite( filename, nymd, nhms, 0, w_taubc, rc, &
                            verbose=.true., new=new, freq=freq)
    endif
    if(doing_cc) then 
     filename = trim(outfile(1:lenfile)//trim(chstr)//'.cc')
     w_taucc%rh = w_c%rh
     call Chem_BundleWrite( filename, nymd, nhms, 0, w_taucc, rc, &
                            verbose=.true., new=new, freq=freq)
    endif
    if(doing_anthro) then 
     filename = trim(outfile(1:lenfile)//trim(chstr)//'.anthro')
     w_tauant%rh = w_c%rh
     call Chem_BundleWrite( filename, nymd, nhms, 0, w_tauant, rc, &
                            verbose=.true., new=new, freq=freq)
    endif


!   Clean up pointers
!   -----------------
    call Chem_BundleDestroy(w_tau,rc)
    if(doing_anthro) call Chem_BundleDestroy(w_tauant,rc)
    if(doing_dust)   call Chem_BundleDestroy(w_taudu,rc)
    if(doing_ss)     call Chem_BundleDestroy(w_tauss,rc)
    if(doing_su)     call Chem_BundleDestroy(w_tausu,rc)
    if(doing_oc)     call Chem_BundleDestroy(w_tauoc,rc)
    if(doing_bc)     call Chem_BundleDestroy(w_taubc,rc)
    if(doing_cc)     call Chem_BundleDestroy(w_taucc,rc)


   enddo  ! channels

!  Clean up allocation of space for delz if doing extinction
   if(iext .gt. 0 .or. ibbck .gt. 0) deallocate(delz)


!  ==================================================================================

   call Chem_BundleDestroy(w_c,rc)
   new = .false.
  enddo   ! idx (time increment in input file)

! Destroy Mie tables
  call Chem_MieDestroy(mie_tables,ier)
  call Chem_RegistryDestroy(regInp, rc)
  call Chem_RegistryDestroy(regOut, rc)

! ----------------------------------------------------------------------------
  contains

  subroutine usage()
  print *
  print *,'Usage: '
  print *,'  Chem_Aod3d.x [-dust -anthro '
  print *,'                -o outfile -t rcfile ] infile'
  print *
  print *, 'where'
  print *
  print *, '-geos4       to specify that the relative humidity of input file'
  print *, '             varies 0 - 100 instead of 0 - 1 as in GEOS-5'
  print *, '-dryaer      to specify to ignore the relative humidity in the'
  print *, '             input file; compute all properties like RH = 0%'
  print *, '-dust        compute for dust separately'
  print *, '-ss          compute for seasalt separately'
  print *, '-su          compute for sulfate separately'
  print *, '-bc          compute for black carbon separately'
  print *, '-oc          compute for organic carbon separately'
  print *, '-cc          compute for total carbon separately'
  print *, '-anthro      compute for SU, BC, and OC separately'
  print *, '-o outfile   output file header'
  print *, '-t rcfile    resource file specifying channels for AOD calc'
  print *, '-airdensfile airdensfilename   filename providing variable AIRDENS'
  print *, '                               if omitted, use infile'
  print *, '-hghtefile   hghtefilename     filename providing variable HGHTE'
  print *, '                               supercedes airdensfile if both provided'
  print *, 'infile       mandatory c_rst file'
  print *
  call exit(1)
  end subroutine usage


  subroutine create_species_bundle(this)
  type(Chem_Bundle) :: this

  call Chem_BundleCreate(regOut, &
                         i1, i2, ig, im, &
                         j1, j2, jg, jm, km, &
                         this, ier)
  if(ier /= 0) call die(myname, 'cannot create bundle')
  this%delp = w_c%delp
  this%rh   = 0.
! NB: as a hack I may use this%rh as a storage variable;
!     I'll reset at the end
  do n = 1, regOut%nq
   this%qa(n)%data3d = 0.0
  end do
  end subroutine create_species_bundle


  subroutine fill(this)
  type(Chem_Bundle) :: this
      if(itau .gt. 0)   this%qa(itau)%data3d(i,j,k) &
                         = this%qa(itau)%data3d(i,j,k) + tau_
      if(iext .gt. 0)   this%qa(iext)%data3d(i,j,k) &
                         = this%qa(iext)%data3d(i,j,k) + tau_/delz(i,j,k)
      if(issa .gt. 0)   this%qa(issa)%data3d(i,j,k) &
                         = this%qa(issa)%data3d(i,j,k) + ssa_*tau_
      if(igasym .gt. 0) this%qa(igasym)%data3d(i,j,k) &
                         = this%qa(igasym)%data3d(i,j,k) + gasym_*ssa_*tau_
      if(immr .gt. 0)   this%qa(immr)%data3d(i,j,k) &
                         = this%qa(immr)%data3d(i,j,k) + w_c%qa(iq)%data3d(i,j,k)
      if(ibbck .gt. 0)  this%qa(ibbck)%data3d(i,j,k) &
                         = this%qa(ibbck)%data3d(i,j,k) + bbck_*qMass/delz(i,j,k)
      if(idepol .gt. 0) then
        this%qa(idepol)%data3d(i,j,k) &
         = this%qa(idepol)%data3d(i,j,k) + (p11_-p22_)*ssa_*tau_
!       See how we are using this%rh; don't forget!
        this%rh(i,j,k) = this%rh(i,j,k) + (p11_+p22_)*ssa_*tau_
      endif

      if(ivol .gt. 0) then
        this%qa(ivol)%data3d(i,j,k) = this%qa(ivol)%data3d(i,j,k) + vol_*qMass/delz(i,j,k)
      endif

      if(iarea .gt. 0) then
        this%qa(iarea)%data3d(i,j,k) = this%qa(iarea)%data3d(i,j,k) + area_*qMass/delz(i,j,k)
      endif

      if(irefr .gt. 0) then
        this%qa(irefr)%data3d(i,j,k) = this%qa(irefr)%data3d(i,j,k) + refr_*vol_*qMass/delz(i,j,k)
      endif

      if(irefi .gt. 0) then
        this%qa(irefi)%data3d(i,j,k) = this%qa(irefi)%data3d(i,j,k) + refi_*vol_*qMass/delz(i,j,k)
      endif

  end subroutine fill


  subroutine normal(this)
  type(Chem_Bundle) :: this
!     Check floating underflow in tau
      do k = 1, km
       do j = 1, jm
        do i = 1, im
         if(itau .gt. 0) this%qa(itau)%data3d(i,j,k) = max(     this%qa(itau)%data3d(i,j,k), &
                                                           tiny(this%qa(itau)%data3d(i,j,k)))
         if(iext .gt. 0) this%qa(iext)%data3d(i,j,k) = max(     this%qa(iext)%data3d(i,j,k), &
                                                           tiny(this%qa(iext)%data3d(i,j,k)))
         if(ibbck .gt. 0) this%qa(ibbck)%data3d(i,j,k) = max(     this%qa(ibbck)%data3d(i,j,k), &
                                                           tiny(this%qa(ibbck)%data3d(i,j,k)))
         if(issa .gt. 0) this%qa(issa)%data3d(i,j,k) = max(     this%qa(issa)%data3d(i,j,k), &
                                                           tiny(this%qa(issa)%data3d(i,j,k)))
         if(issa .gt. 0)   this%qa(issa)%data3d(i,j,k) = this%qa(issa)%data3d(i,j,k)/this%qa(itau)%data3d(i,j,k)
         if(igasym .gt. 0) this%qa(igasym)%data3d(i,j,k) = max(     this%qa(igasym)%data3d(i,j,k), &
                                                           tiny(this%qa(igasym)%data3d(i,j,k)))
         if(igasym .gt. 0) this%qa(igasym)%data3d(i,j,k) = this%qa(igasym)%data3d(i,j,k) &
                                                    /(this%qa(issa)%data3d(i,j,k)*this%qa(itau)%data3d(i,j,k))
         if(ietob .gt. 0) this%qa(ietob)%data3d(i,j,k) = this%qa(itau)%data3d(i,j,k)/delz(i,j,k)/this%qa(ibbck)%data3d(i,j,k)
         if(ietob .gt. 0) this%qa(ietob)%data3d(i,j,k) = max(     this%qa(ietob)%data3d(i,j,k), &
                                                           tiny(this%qa(ietob)%data3d(i,j,k)))
         if(idepol .gt. 0) this%qa(idepol)%data3d(i,j,k) = this%qa(idepol)%data3d(i,j,k)/max(this%rh(i,j,k),tiny(this%rh))
         if(idepol .gt. 0) this%qa(idepol)%data3d(i,j,k) = max(     this%qa(idepol)%data3d(i,j,k), &
                                                           tiny(this%qa(idepol)%data3d(i,j,k)))
         if(ivol .gt. 0) this%qa(ivol)%data3d(i,j,k) = max(     this%qa(ivol)%data3d(i,j,k), &
                                                           tiny(this%qa(ivol)%data3d(i,j,k)))
         if(iarea .gt. 0) this%qa(iarea)%data3d(i,j,k) = max(     this%qa(iarea)%data3d(i,j,k), &
                                                           tiny(this%qa(iarea)%data3d(i,j,k)))
         if(ireff .gt. 0) this%qa(ireff)%data3d(i,j,k) = 3./4.*this%qa(ivol)%data3d(i,j,k)/this%qa(iarea)%data3d(i,j,k)
         if(ireff .gt. 0) this%qa(ireff)%data3d(i,j,k) = max(     this%qa(ireff)%data3d(i,j,k), &
                                                           tiny(this%qa(ireff)%data3d(i,j,k)))
         if(irefr .gt. 0) this%qa(irefr)%data3d(i,j,k) = this%qa(irefr)%data3d(i,j,k)/this%qa(ivol)%data3d(i,j,k)
         if(irefr .gt. 0) this%qa(irefr)%data3d(i,j,k) = max(     this%qa(irefr)%data3d(i,j,k), &
                                                           tiny(this%qa(irefr)%data3d(i,j,k)))
         if(irefi .gt. 0) this%qa(irefi)%data3d(i,j,k) = this%qa(irefi)%data3d(i,j,k)/this%qa(ivol)%data3d(i,j,k)
         if(irefi .gt. 0) this%qa(irefi)%data3d(i,j,k) = max(     this%qa(irefi)%data3d(i,j,k), &
                                                           tiny(this%qa(irefi)%data3d(i,j,k)))
        enddo
       enddo
      enddo
  end subroutine normal


  subroutine backscatter(this)
  type(Chem_Bundle) :: this

!  Attenuated backscatter from space (I think, right, k = 1 is TOA)
   if(iabck1 .gt. 0) then 
    this%qa(iabck1)%data3d(:,:,1) = this%qa(ibbck)%data3d(:,:,1)*exp(-this%qa(itau)%data3d(:,:,1))
    do k = 2, km
     do j = 1, jm
      do i = 1, im
       taulev = 0.
       do kk = 1, k-1
        taulev = taulev + this%qa(itau)%data3d(i,j,kk)
       enddo
       taulev = taulev + 0.5 * this%qa(itau)%data3d(i,j,k)
       this%qa(iabck1)%data3d(i,j,k) = this%qa(ibbck)%data3d(i,j,k)*exp(-2.*taulev)
      enddo
     enddo
    enddo
   endif

!  Attenuated backscatter from surface (I think, right, k = 1 is TOA)
   if(iabck0 .gt. 0) then 
    this%qa(iabck0)%data3d(:,:,km) = this%qa(ibbck)%data3d(:,:,km)*exp(-this%qa(itau)%data3d(:,:,km))
    do k = km-1, 1, -1
     do j = 1, jm
      do i = 1, im
       taulev = 0.
       do kk = km, k+1, -1
        taulev = taulev + this%qa(itau)%data3d(i,j,kk)
       enddo
       taulev = taulev + 0.5 * this%qa(itau)%data3d(i,j,k)
       this%qa(iabck0)%data3d(i,j,k) = this%qa(ibbck)%data3d(i,j,k)*exp(-2.*taulev)
      enddo
     enddo
    enddo
   endif
  end subroutine backscatter


  subroutine delz_
   implicit none

   real :: Tv, mixr, pe(km+1)
   real, parameter :: ptop = 1.
!  Construct edge pressures
!  ------------------------
   pe(1) = ptop
   do k = 2, km + 1
      pe(k) = pe(k-1) + w_c%delp(i,j,k-1)
   end do

!  Construct mid-layer pressures and layer thickness
!  -------------------------------------------------
   do k = 1, km
      mixr = q(i,j,k) / ( 1.0 - q(i,j,k) ) ! Mixing ratio from specific humidity
      Tv = T(i,j,k) * ( 1 + 0.61 * mixr )
      delz(i,j,k)  = Rgas * Tv * log(pe(k+1)/pe(k)) / grav
   end do

  end subroutine delz_


  subroutine check_infile
   ! This will die non-gracefully if any problems reading file
   integer fid, err, hour, min, timInc
   integer l_im, l_jm, l_km, l_lm, l_nvars, l_ngatts, incsecs, yyyymmdd_beg, hhmmss_beg
   integer, parameter :: READ_ONLY = 1

   !  Open the file
   !  -------------
   call GFIO_Open ( infile, READ_ONLY, fid, err )
   call GFIO_DimInquire ( fid, l_im, l_jm, l_km, l_lm, l_nvars, l_ngatts, err)
   call GetBegDateTime ( fid, yyyymmdd_beg, hhmmss_beg, incSecs, err )
   timInc = 000100   ! default: 1 minute
   if ( l_lm .ge. 1 ) then   !ams: changed lm.gt.1 to lm.ge.1
           hour = incSecs/3600
           if (hour == 0) hour=1
           min = mod(incSecs,3600*hour)/60
           timInc = incSecs/3600*10000 + mod(incSecs,3600)/60*100 + mod(incSecs,60)
   end if
   call GFIO_close ( fid, err )
   ndx = l_lm
   freq = timInc
  end subroutine check_infile

end
