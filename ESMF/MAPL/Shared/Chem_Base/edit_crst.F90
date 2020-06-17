!

! Simple program to read a chem bundle, modify its values and write it back.
!
 
     program edit_crst

     use Chem_RegistryMod
     use Chem_BundleMod

     character(len=255) :: in_file, out_file 
     character(len=255) :: in_add_file
     character(len=255) :: dir
     real, pointer :: lon(:), lat(:), pm(:,:,:)
     real, allocatable :: co2_tmp(:,:,:) 

     type(Chem_Registry) :: reg_in, reg_out
     type(Chem_Registry) :: reg
     type(Chem_Bundle)   :: w_in, w_out
     type(Chem_Bundle)   :: w

     integer :: k,ier,nhms,nymd,im,jm,km
     integer :: ii, jj, i, j

!     real, pointer :: q(:,:,:,:)


!    Look in 
!    /output0/dao_ops/GEOS-4_AVE_Houston/fvcm_ave_01/rs/Y2005/M06
     in_file  = '/nobackup2/colarco/aqua_gfedv2.c_rst'

!    Input file from CTM run
     in_add_file  = '/nobackup2/colarco/init_co2_20000630.c32.dat.biged'

!    Output file to create
     out_file = '/nobackup2/colarco/aqua_gfedv2.c_rst.co2'

     reg_in  = Chem_RegistryCreate ( ier, rcfile = 'Chem_Registry_old.rc' )
     reg_out = Chem_RegistryCreate ( ier, rcfile = 'Chem_Registry_new.rc' )
     call Chem_registryprint(reg_out)
     if ( ier /= 0 ) then
        print *,'oops, error'
        call exit(1)
     end if


!    Read initial chem bundle
!    ------------------------
     call Chem_BundleRead ( trim(in_file), nymd, nhms, w_in, ier, &
                            chemReg = reg_in )
     print *, 'nymd, nhms = ', nymd, nhms

!    Fill in the static portion of the outgoing chem bundle
!    ------------------------------------------------------
     im = w_in%grid%im
     jm = w_in%grid%jm
     km = w_in%grid%km
     call Chem_BundleCreate ( reg_out, im, jm, km, &
                              w_out, ier )
     call Chem_RegistryPrint(w_out%reg)
     w_out%delp = w_in%delp
     w_out%rh = w_in%rh
     w_out%qa = w_in%qa
     w_out%q(:,:,:,1) = w_In%q(:,:,:,1)


!    Now map the constituent values
!    Hard-wired for all constituents to assume all tracers have same
!    dimensions, except CO and CO2 treated separately below
     if(reg_out%doing_DU) then
      ii = reg_out%i_DU
      jj = reg_out%j_DU
      if(reg_in%doing_DU) then
       i = reg_in%i_DU
       j = reg_in%j_DU
       w_out%q(:,:,:,ii:jj) = w_in%q(:,:,:,i:j)
      else
       w_out%q(:,:,:,ii:jj) = 0.
      endif
     endif
     if(reg_out%doing_OC) then
      ii = reg_out%i_OC
      jj = reg_out%j_OC
      if(reg_in%doing_OC) then
       i = reg_in%i_OC
       j = reg_in%j_OC
       w_out%q(:,:,:,ii:jj) = w_in%q(:,:,:,i:j)
      else
       w_out%q(:,:,:,ii:jj) = 0.
      endif
     endif
     if(reg_out%doing_BC) then
      ii = reg_out%i_BC
      jj = reg_out%j_BC
      if(reg_in%doing_BC) then
       i = reg_in%i_BC
       j = reg_in%j_BC
       w_out%q(:,:,:,ii:jj) = w_in%q(:,:,:,i:j)
      else
       w_out%q(:,:,:,ii:jj) = 0.
      endif
     endif
     if(reg_out%doing_SU) then
      ii = reg_out%i_SU
      jj = reg_out%j_SU
      if(reg_in%doing_SU) then
       i = reg_in%i_SU
       j = reg_in%j_SU
       w_out%q(:,:,:,ii:jj) = w_in%q(:,:,:,i:j)
      else
       w_out%q(:,:,:,ii:jj) = 0.
      endif
     endif
     if(reg_out%doing_SS) then
      ii = reg_out%i_SS
      jj = reg_out%j_SS
      if(reg_in%doing_SS) then
       i = reg_in%i_SS
       j = reg_in%j_SS
       w_out%q(:,:,:,ii:jj) = w_in%q(:,:,:,i:j)
      else
       w_out%q(:,:,:,ii:jj) = 0.
      endif
     endif
     if(reg_out%doing_O3) then
      ii = reg_out%i_O3
      jj = reg_out%j_O3
      if(reg_in%doing_O3) then
       i = reg_in%i_O3
       j = reg_in%j_O3
       w_out%q(:,:,:,ii:jj) = w_in%q(:,:,:,i:j)
      else
       w_out%q(:,:,:,ii:jj) = 0.
      endif
     endif
     if(reg_out%doing_SC) then
      ii = reg_out%i_SC
      jj = reg_out%j_SC
      if(reg_in%doing_SC) then
       i = reg_in%i_SC
       j = reg_in%j_SC
       w_out%q(:,:,:,ii:jj) = w_in%q(:,:,:,i:j)
      else
       w_out%q(:,:,:,ii:jj) = 0.
      endif
     endif
     if(reg_out%doing_AC) then
      ii = reg_out%i_AC
      jj = reg_out%j_AC
      if(reg_in%doing_AC) then
       i = reg_in%i_AC
       j = reg_in%j_AC
       w_out%q(:,:,:,ii:jj) = w_in%q(:,:,:,i:j)
      else
       w_out%q(:,:,:,ii:jj) = 0.
      endif
     endif
     if(reg_out%doing_XX) then
      ii = reg_out%i_XX
      jj = reg_out%j_XX
      if(reg_in%doing_XX) then
       i = reg_in%i_XX
       j = reg_in%j_XX
       w_out%q(:,:,:,ii:jj) = w_in%q(:,:,:,i:j)
      else
       w_out%q(:,:,:,ii:jj) = 0.
      endif
     endif


!    Set the CO tracers to 50 ppbv
     if(reg_out%doing_CO) then
      ii = reg_out%i_CO
      jj = reg_out%j_CO
      w_out%q(:,:,:,ii:jj) = 50.e-9
     endif

!    Get and fill in the CO2 tracers
!    For the first tracer use the UCTM output
!    For the others set to 350 ppmv
     if(reg_out%doing_CO2) then
      ii = reg_out%i_CO2
      jj = reg_out%j_CO2
      allocate(co2_tmp(im,jm,km))
      open(10,file=TRIM(in_add_file),form='unformatted', &
              status='old',action='read')
      read(10) co2_tmp
      close(10)
      k = ii
      w_out%q(:,:,:,k) = co2_tmp(:,:,:)
      if(jj .gt. ii) then
       do k=ii+1, jj
        w_out%q(:,:,:,k) = 350.e-6
       enddo
      endif
      print *, 'q=',w_out%q(1,1,1,ii:jj)
      deallocate(co2_tmp)
     endif

!    Write restart
!    ------------
     call Chem_BundleWrite ( trim(out_file), nymd, nhms, 1, w_out, ier )

!    Read it back in
     reg = Chem_RegistryCreate ( ier, rcfile = 'Chem_Registry.rc' )
     call Chem_registryprint(reg)
     call Chem_BundleCreate ( reg, im, jm, km, w, ier )
     call Chem_BundleRead ( trim(out_file), nymd, nhms, w, ier, &
                            chemReg = reg )


     end program edit_crst


