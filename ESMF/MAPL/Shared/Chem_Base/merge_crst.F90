!
!  Merge 2 Chem_Bundles. Quick and dirty.
!
 
     program merge_crst

     use Chem_RegistryMod
     use Chem_BundleMod

     character(len=255) :: in_file(2), in_reg(2), out_file 
     integer :: ier

     type(Chem_Registry) :: reg_in(2), reg_out
     type(Chem_Registry) :: reg
     type(Chem_Bundle)   :: w_in(2), w_out

     integer :: k,nhms,nymd,nymd2,nhms2,im,jm,km
     integer :: i, j, i_in, i_out, j_in, j_out

!    Take Aerosols and global CO from CRAVE run
!    ------------------------------------------
     in_file(1) = '/output/dao_ops/GEOS-4_CRAVE/a_flk_04C/rs/Y2006/M02/a_flk_04C.rst.chem.20060221_15z.bin'

     in_reg(1) = '/output/dao_ops/GEOS-4_CRAVE/a_flk_04C/run/Chem_Registry.rc'

!    Take O3 from Eric's file
!    ------------------------
     in_file(2) = '/nobackup1/enielsen/fvchem/INTEX2006/c55/recycle/c55.c_rst.20060101'
     in_reg(2) = '/nobackup1/enielsen/fvchem/INTEX2006/c55/run/Chem_Registry.rc'

!    Output file to create
     out_file = '/nobackup1/dasilva/rs4intex/20060221_15/a_flk_04C.c_rst'

!    Create registries
!    -----------------
     reg_in(1)  = Chem_RegistryCreate ( ier, rcfile = in_reg(1) )
     call Chem_RegistryPrint(reg_in(1))
     reg_in(2)  = Chem_RegistryCreate ( ier, rcfile = in_reg(2) )
     call Chem_RegistryPrint(reg_in(2))
     reg_out = Chem_RegistryCreate ( ier, rcfile = 'Chem_Registry.rc' )
     call Chem_RegistryPrint(reg_out)
     if ( ier /= 0 ) then
        print *,'oops, error'
        call exit(1)
     end if

!    Read in source bundles
!    ------------------------
     call Chem_BundleRead ( trim(in_file(1)), nymd, nhms, w_in(1), ier, &
                            chemReg = reg_in(1) )
     print *, 'Read first bundle at on ', nymd, nhms
     call Chem_BundleRead ( trim(in_file(2)), nymd2, nhms2, w_in(2), ier, &
                            chemReg = reg_in(2) )
     print *, 'nymd, nhms = ', nymd2, nhms2

!    Fill in the static portion of the outgoing chem bundle
!    ------------------------------------------------------
     im = w_in(1)%grid%im
     jm = w_in(1)%grid%jm
     km = w_in(1)%grid%km
     call Chem_BundleCreate ( reg_out, im, jm, km, &
                              w_out, ier )
     call Chem_RegistryPrint(w_out%reg)

!    Most things come from the CRAVE run
!    -----------------------------------
     w_out%delp = w_in(1)%delp
     w_out%rh = w_in(1)%rh
     w_out%qa = w_in(1)%qa
     w_out%grid = w_in(1)%grid

     w_out% q = 0.0  ! start with a clean slate

!    H2O
!    ---
     i_in  = reg_in(1)%i_H2O; j_in  = reg_in(1)%j_H2O;
     i_out = reg_out%i_H2O;   j_out =   reg_out%j_H2O;
     w_out%q(:,:,:,i_out:j_out) = w_in(1)%q(:,:,:,i_in:j_in)

!    DU
!    ---
     i_in  = reg_in(1)%i_DU; j_in  = reg_in(1)%j_DU;
     i_out = reg_out%i_DU;   j_out =   reg_out%j_DU;
     w_out%q(:,:,:,i_out:j_out) = w_in(1)%q(:,:,:,i_in:j_in)

!    SS
!    ---
     i_in  = reg_in(1)%i_SS; j_in  = reg_in(1)%j_SS;
     i_out = reg_out%i_SS;   j_out =   reg_out%j_SS;
     w_out%q(:,:,:,i_out:j_out) = w_in(1)%q(:,:,:,i_in:j_in)

!    BC
!    ---
     i_in  = reg_in(1)%i_BC; j_in  = reg_in(1)%j_BC;
     i_out = reg_out%i_BC;   j_out =   reg_out%j_BC;
     w_out%q(:,:,:,i_out:j_out) = w_in(1)%q(:,:,:,i_in:j_in)

!    OC
!    ---
     i_in  = reg_in(1)%i_OC; j_in  = reg_in(1)%j_OC;
     i_out = reg_out%i_OC;   j_out =   reg_out%j_OC;
     w_out%q(:,:,:,i_out:j_out) = w_in(1)%q(:,:,:,i_in:j_in)

!    SU
!    ---
     i_in  = reg_in(1)%i_SU; j_in  = reg_in(1)%j_SU;
     i_out = reg_out%i_SU;   j_out =   reg_out%j_SU;
     w_out%q(:,:,:,i_out:j_out) = w_in(1)%q(:,:,:,i_in:j_in)

!    Get O3 from Eric
!    ----------------
     i_in  = reg_in(2)%i_O3; j_in  = reg_in(2)%j_O3;
     i_out = reg_out%i_O3;   j_out =   reg_out%j_O3;
     w_out%q(:,:,:,i_out:j_out) = w_in(2)%q(:,:,:,i_in:j_in)

!    CO from Eric
!    ------------
     i_in  = reg_in(2)%i_CO; j_in  = reg_in(2)%j_CO;
     i_out = reg_out%i_CO;   j_out =   reg_out%j_CO;
     w_out%q(:,:,:,i_out:j_out) = w_in(2)%q(:,:,:,i_in:j_in)


!    Write output file with date from Ops
!    ------------------------------------
     call Chem_BundleWrite ( trim(out_file), nymd, nhms, 0, w_out, ier )

     print *, 'all done'

     end program merge_crst


