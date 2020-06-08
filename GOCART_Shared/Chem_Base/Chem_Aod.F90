! Compute the Aerosol Optical Thickness for an input Chem Bundle
! - get the chemistry registry
! - pass the chem bundle to the AOD routine

  program Chem_Aod

  use m_die, only: die
  use Chem_MieMod
  use Chem_RegistryMod
  use Chem_BundleMod

  implicit none

  character(len=*), parameter :: myname = 'chem_aod'
  type(Chem_Mie)      :: mie_tables
  type(Chem_Registry) :: regInp   ! chemistry registry
  type(Chem_Bundle)   :: w_c      ! chemistry bundle
  type(Chem_Bundle)   :: w_tau    ! tau chemistry bundle
  type(Chem_Bundle)   :: w_tauabs ! tau absorption chemistry bundle
  type(Chem_Bundle)   :: w_ssa    ! ssa chemistry bundle
  integer :: i, j, k, im, jm, km, idx
  integer :: i1, i2, ig, j1, j2, jg, ik, iq
  integer :: nymd, nhms, freq, rc, ier
  integer :: idxTable
  integer :: out_fid
  integer iarg, iargc, argc, lenfile
  logical :: only_taod
  logical :: doing_tauabs2d
  logical :: doing_totext2d
  logical :: doing_ssa2d
  logical :: doing_geos4
  logical :: doing_dry   ! if true, calculate like rh = 0%
  logical :: new, verbose
  real :: channel, tau_, ssa_, scalerh, maxRH
  character(len=255) :: infile, outfile, filename, rcfile, argv
  character(len=255) :: filename0
  character(len=14)  :: datestr
  character(len=8)   :: yyyymmddstr
  character(len=4)   :: hhnnstr

! Parse the command line (see usage() below)
  argc = iargc()
  if(argc .lt. 1) call usage()
  iarg = 0
  outfile = 'chem_aod'
  rcfile  = 'Aod_Registry.rc'
  only_taod=.false.
  doing_tauabs2d = .false.
  doing_totext2d = .false.
  doing_ssa2d = .false.
  doing_geos4 = .false.
  verbose = .false.
  doing_dry = .false.
  do i = 0, 32767
   iarg = iarg+1
   if(iarg .gt. argc) exit
   call GetArg(iarg, argv)
   select case(argv)
    case ("-geos4")
     doing_geos4 = .true.
    case ("-dryaer")
     doing_dry = .true.
    case ("-v")
     verbose = .true.
    case ("-tauabs2d")
     doing_tauabs2d = .true.
    case ("-totext2d")
     doing_totext2d = .true.
!_RT doing_tauabs2d = .true.
    case ("-ssa2d")
     doing_ssa2d = .true.
    case ("-only_taod")
     only_taod = .true.
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

! Hardwired: Read the input chemistry registry from Chem_Registry.rc
! This registry file describes the input chemistry bundle being 
! operated on.
! ------------------------------------------------------------------
  regInp = Chem_RegistryCreate(ier,'Chem_MieRegistry.rc')
  if(ier /= 0) call die(myname, 'cannot create registry')
  if(verbose) call Chem_RegistryPrint(regInp)

! Hardwired: use the Chem_MieMod function to create the Mie tables
! -------------------------------------------------------------------------
  mie_tables = Chem_MieCreate(rcfile,ier)

! Hardwired: we know the chem bundle files contain four time steps per file
! We will loop over all the times in the file and write them out
! -------------------------------------------------------------------------
  new = .true.
  do idx = 1, 1

!  Read the chemistry bundle from the infile
!  -------------------------------------------------
   call Chem_BundleRead(infile, nymd, nhms, w_c, rc, freq=freq, &
                        ChemReg=regInp, timidx=idx)

   print *, 'Computing AOD for ', nymd, nhms

!  Check the RH seems sane
   if(.not. doing_geos4 .and. .not. doing_dry) then
    maxrh = maxval(w_c%rh)
    if(maxrh .gt. 2.) then
     print *, 'Maximum RH value = ', maxRH
     print *, 'Should you have chosen "-geos4" as a command line option?'
    endif
   endif

!  ==================================================================================
!  The enclosed bundle of code selects on what calculation we run and what is written

!  Simple for now: only do tau2d
!  Create the output Chem_Bundle
   i1 = w_c%grid%i1
   i2 = w_c%grid%i2
   ig = w_c%grid%ig
   j1 = w_c%grid%j1
   j2 = w_c%grid%j2
   jg = w_c%grid%jg
   im = w_c%grid%im
   jm = w_c%grid%jm
   km = mie_tables%nch

   call Chem_BundleCreate(regInp, &
                          i1, i2, ig, im, &
                          j1, j2, jg, jm, km, &
                          w_tau, ier, &
                          lev=mie_tables%channels, levUnits="m")
   if(doing_ssa2d .or. doing_tauabs2d) then
    call Chem_BundleCreate(regInp, &
                           i1, i2, ig, im, &
                           j1, j2, jg, jm, km, &
                           w_ssa, ier, &
                           lev=mie_tables%channels, levUnits="m")
   endif

   if(doing_tauabs2d) then
    call Chem_BundleCreate(regInp, &
                           i1, i2, ig, im, &
                           j1, j2, jg, jm, km, &
                           w_tauabs, ier, &
                           lev=mie_tables%channels, levUnits="m")
   endif

   if(ier /= 0) call die(myname, 'cannot create tau2d bundle')

   do ik = 1, km
    channel = mie_tables%channels(ik)

    do iq = 1, mie_tables%nq
     idxTable = Chem_MieQueryIdx(mie_tables,mie_tables%vname(iq),rc)

     if(idxTable .ne. -1) then
      do k = 1, w_c%grid%km
      do j = 1, jm
      do i = 1, im
         if ( doing_ssa2d .or. doing_tauabs2d ) then
            call Chem_MieQuery(mie_tables, idxTable, float(ik), &
                            w_c%qa(iq)%data3d(i,j,k)*w_c%delp(i,j,k)/9.80616, &
                            w_c%rh(i,j,k) * scaleRH, tau=tau_, ssa=ssa_)
         else
            call Chem_MieQuery(mie_tables, idxTable, float(ik), &
                            w_c%qa(iq)%data3d(i,j,k)*w_c%delp(i,j,k)/9.80616, &
                            w_c%rh(i,j,k) * scaleRH, tau=tau_)
         endif
         w_tau%qa(iq)%data3d(i,j,ik) = w_tau%qa(iq)%data3d(i,j,ik) + tau_
         if(doing_ssa2d) w_ssa%qa(iq)%data3d(i,j,ik) = w_ssa%qa(iq)%data3d(i,j,ik) + tau_*ssa_
         if(doing_tauabs2d) w_tauabs%qa(iq)%data3d(i,j,ik) = &
                            w_tauabs%qa(iq)%data3d(i,j,ik) + (1.-ssa_)*tau_
      enddo
      enddo
      enddo
     endif
    enddo

   enddo

 if(doing_ssa2d) then
      do iq = 1, mie_tables%nq
       idxTable = Chem_MieQueryIdx(mie_tables,mie_tables%vname(iq),rc)
       if(idxTable .ne. -1) then
         w_ssa%qa(iq)%data3d = w_ssa%qa(iq)%data3d / w_tau%qa(iq)%data3d
       endif
      end do
 end if

!  Write the Chem_Bundle out
   write(yyyymmddstr,'(i8.8)') nymd
   write(hhnnstr,'(i4.4)') nhms/100
   datestr = yyyymmddstr//'_'//hhnnstr//'z'


  if(doing_totext2d) then 
    if (only_taod) then
       filename0 = trim(outfile)
    else
       filename0 = trim(outfile(1:lenfile)//'.ext_Nc.'//datestr//'.nc4')
    endif
    call cmp_totext(filename0,w_tau,nymd,nhms,0,1,out_fid,freq=freq)
  endif

  if(.not.only_taod) then 

    filename = trim(outfile(1:lenfile)//'.taod_Nc.'//datestr//'.nc4')
    call Chem_BundleWrite( filename, nymd, nhms, 0, w_tau, rc, &
                          verbose=verbose, new=new, freq=freq)

   if(doing_ssa2d) then 
      filename = trim(outfile(1:lenfile)//'.tssa_Nc.'//datestr//'.nc4')
      call Chem_BundleWrite( filename, nymd, nhms, 0, w_ssa, rc, &
                             verbose=verbose, new=new, freq=freq)
   endif

   if(doing_tauabs2d) then 

    if(doing_totext2d) then 
     filename0 = trim(outfile(1:lenfile)//'.ext_Nc.'//datestr//'.nc4')
     call cmp_totext(filename0,w_tauabs,nymd,nhms,0,2,out_fid,freq=freq)
    endif

    filename = trim(outfile(1:lenfile)//'.aaod_Nc.'//datestr//'.nc4')
    call Chem_BundleWrite( filename, nymd, nhms, 0, w_tauabs, rc, &
                           verbose=verbose, new=new, freq=freq)
   endif
  endif ! .not.only_taod

!  ==================================================================================

!  Don't overwrite the file
!  ------------------------
   new = .false.

  enddo   ! idx (time increment in input file)

! Destroy Mie tables
  call Chem_BundleDestroy(w_tau, rc)
  if(doing_ssa2d) call Chem_BundleDestroy(w_ssa, rc)
  if(doing_tauabs2d) call Chem_BundleDestroy(w_tauabs, rc)
  call Chem_MieDestroy(mie_tables,ier)
  call Chem_RegistryDestroy(regInp, rc)

! ----------------------------------------------------------------------------
  contains

!-------------------------------------------------------------------------
!      NASA/GSFC Global Modeling & Assimilation Office, Code 900.3       !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  cmp_totext --- Computes Total AOD and ABS_AOD.
! 
! !INTERFACE:
!


  subroutine cmp_totext( outFile,w_c,nymd,nhms,prec,iflag,out_fid,freq)
!
! !USES:
!
  use m_chars, only: uppercase
  implicit NONE

!
! !INPUT PARAMETERS: 
!
   integer, intent(in), OPTIONAL  :: freq    ! time frequency (HHMMSS) for
   integer,             intent(in)   :: prec ! precision
                                             ! multiple instance files
                                             ! (default: 060000)

  type(Chem_Bundle), intent(in)     :: w_c    ! chemical bundle
  integer      :: im_e,jm_e,km_e,i,j,k,iq,iflag

! !OUTPUT PARAMETER:

     real, allocatable   :: total_ext(:,:,:)
     real, allocatable   :: valid_range(:,:), packing_range(:,:)
     real, allocatable   :: lat_e(:), lon_e(:), lev_e(:),kmVar_e(:)
     character(len=*)    :: outFile            ! Output file name
     integer             :: out_fid,err4
     integer             :: nymd,nhms,timeinc
     integer             :: nVars_new,rc
     character(len=64)   :: outVars(2)  ! output variable names (nVars
     character(len=64)   :: outUnits(2)  ! Units of output variables (nVars)
     character(len=256)  :: out_title(2) ! output title
     character(len=256)  :: source, contact, levunits,title
!
! !DESCRIPTION: Computes Total AOD & Total Absorption AOD.
!               
! !REVISION HISTORY: 
!
!  16Nov2011 Ravi      Initial code 
!
!EOP
!------------------------------------------------------------------
  
      im_e = w_c%grid%im
      jm_e = w_c%grid%jm
      km_e = w_c%grid%km

     allocate(total_ext(im_e,jm_e,km_e),stat=err4)

     total_ext = 0.0
     do iq = 1, w_c%reg%nq
       if(trim(uppercase(w_c%reg%vname(iq)(1:2))) == 'DU'       .or. &
          trim(uppercase(w_c%reg%vname(iq)(1:2))) == 'SS'       .or. &
          trim(uppercase(w_c%reg%vname(iq)))      == 'SO4'      .or. &
          trim(uppercase(w_c%reg%vname(iq)))      == 'NO3AN1'   .or. &
          trim(uppercase(w_c%reg%vname(iq)))      == 'NO3AN2'   .or. &
          trim(uppercase(w_c%reg%vname(iq)))      == 'NO3AN3'   .or. &
          trim(uppercase(w_c%reg%vname(iq)))      == 'BCPHOBIC' .or. &
          trim(uppercase(w_c%reg%vname(iq)))      == 'BCPHILIC' .or. &
          trim(uppercase(w_c%reg%vname(iq)))      == 'OCPHOBIC' .or. &
          trim(uppercase(w_c%reg%vname(iq)))      == 'OCPHILIC') then

          print *, 'w_c%reg%vname(iq) = ', trim(w_c%reg%vname(iq))

          do k = 1,km_e
           do j = 1,jm_e
            do i = 1,im_e
             if(w_c%qa(iq)%data3d(i,j,k) < 1.e+10) then
               total_ext(i,j,k) = total_ext(i,j,k) + w_c%qa(iq)%data3d(i,j,k)
             endif
            end do
           end do
          end do
       endif
     end do

        nvars_new = 2
        if(iflag == 1) then
          allocate ( packing_range(2,nVars_new), valid_range(2,nVars_new), stat=err4 )
          allocate ( kmVar_e(nvars_new), stat=err4 )
          outVars(1)   = 'taod'
          out_title(1) = 'Total Aerosol Optical Depth'
          outVars(2)   = 'aaod'
          out_title(2) = 'Total Absorption Aerosol Optical Depth'
          outUnits(1)  = 'kg/kg'
          outUnits(2)  = 'kg/kg'
          source  = 'Data Assimilation Office, NASA/GSFC'
          contact = 'data@gmao.gsfc.nasa.gov'
          title   = "unknown"
          levunits = 'm'
          kmVar_e = km_e

!         Cannot handle cubed sphere for now
!         ----------------------------------
          if ( w_c%grid%cubed_sphere ) then
             call die('chem_aod','cannot yet handle cubed sphere')
          end if
          allocate ( lat_e(jm_e), lon_e(im_e), lev_e(km_e),stat=err4 )
          lat_e = w_c%grid%lat(1,:)
          lon_e = w_c%grid%lon(:,1)
          lev_e = w_c%grid%lev

          do j = 1, nVars_new
           do i = 1, 2
            valid_range(i,j) = w_c%missing_value
            packing_range(i,j) = w_c%missing_value
           end do
          end do

          if ( present(freq) ) then
             timeinc = freq
          else
             timeinc = 060000
          end if



          call GFIO_Create ( outFile, title, source, contact, w_c%missing_value, &
                             im_e, jm_e, km_e, lon_e, lat_e, Lev_e, levunits,    &
                             nymd,nhms,timeinc,                                  &
                             nVars_new, outVars, out_title, outUnits,            &
                             kmVar_e,valid_range,packing_range,prec,             &
                             out_fid, rc )

          if ( rc /= 0 )  call die (myname, 'wrong in GFIO_Create')
          deallocate ( packing_range, valid_range,lat_e,lon_e,lev_e,kmVar_e)
        endif

        call GFIO_PutVar (out_fid,outVars(iflag),nymd,nhms,  &
                          im_e, jm_e, 1, km_e, total_ext,rc )
        if ( rc /= 0 ) call die (myname, 'something wrong in GFIO_PutVarT for 3D file')

        deallocate(total_ext)

  end subroutine cmp_totext
! -----------------------------------------------------

  subroutine usage()
  print *
  print *,'Usage: '
  print *,'  Chem_Aod.x [-tauabs2d -ssa2d '
  print *,'              -o outfile -t rcfile ] infile'
  print *
  print *, 'where'
  print *
  print *, '-geos4       to specify that the relative humidity of input file'
  print *, '             varies 0 - 100 instead of 0 - 1 as in GEOS-5'
  print *, '-dryaer      to specify to ignore the relative humidity in the'
  print *, '             input file; compute all properties like RH = 0%'
  print *, '-tauabs2d    request column integrated absorption aerosol optical thickness'
  print *, '-totext2d    request total aerosol optical thickness and '
  print *, '             column integrated absorption aerosol optical thickness'
  print *, '             (-tauabs2d default.)'
  print *, '-ssa2d       request column integrated single scattering albedo'
  print *, '-o expid     filename will look like expid.ext_nx.YYYYMMDD_HHNNz.nc4'
  print *, '             filename for totext2d will look like expid.ext_nc.YYYYMMDD_HHNNz.nc4'
  print *, '-t rcfile    resource file specifying channels for AOD calc'
  print *, '-v           request verbose output'
  print *, 'infile       mandatory input aer_v file'
  print *
  call exit(1)
  end subroutine usage

end
