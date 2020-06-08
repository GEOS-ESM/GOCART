!
! 1) First version added rough estimates SO2 for the Montserrat volcano on 
!    May 20th, 2006. 
! 2) On May 15, SO2 pattern had not moved west as much as seen by OMI. Here
!    we apply a simple relocation of that feature to where OMI reported it.
!
!
 
     program edit_crst

     use Chem_RegistryMod
     use Chem_BundleMod

     character(len=255) :: in_file, out_file 
     character(len=255) :: in_add_file
     character(len=255) :: dir, time
     real, pointer :: lon(:), lat(:), pm(:,:,:)
     real, allocatable :: co2_tmp(:,:,:), tmp(:,:,:)

     type(Chem_Registry) :: reg_in, reg_out
     type(Chem_Registry) :: reg
     type(Chem_Bundle)   :: w_in, w_out
     type(Chem_Bundle)   :: w

     integer :: k,ier,nhms,nymd,im,jm,km
     integer :: ii, jj, i, j, is, js, i1, i2, i1r, i2r, j1, j2, k1, k2, i_so2
     real :: q_so2

!     real, pointer :: q(:,:,:,:)

     time = '20060520_15z' ! eruption restart
     time = '20060525_21z' ! first relocation

!    Look in 
!    /output0/dao_ops/GEOS-4_AVE_Houston/fvcm_ave_01/rs/Y2005/M06
     in_file  = 'a_flk_04C.rst.chem.'//trim(time)//'.bin'

!    Output file to create
     out_file = 'a_flk_04C.c_rst.'//time

     reg_in  = Chem_RegistryCreate ( ier, rcfile = 'Chem_Registry.rc' )
     reg_out = Chem_RegistryCreate ( ier, rcfile = 'Chem_Registry.rc' )
!     call Chem_registryprint(reg_out)
     if ( ier /= 0 ) then
        print *,'oops, error'
        call exit(1)
     end if


!    Read initial chem bundle
!    ------------------------
     call Chem_BundleRead ( trim(in_file), nymd, nhms, w_in, ier, &
                            chemReg = reg_in )
     print *, 'nymd, nhms = ', nymd, nhms

!    Alter the SO2 amounts (see /home/dasilva/out/montserrat.m on calculon)
!    ----------------------------------------------------------------------
#ifdef MONTSERRAT_ERUPTION

     k1 = 38 
     k2 = 40
     is = -49
     if ( is <= 0 ) is = w_in%grid%im + is 
     js = 108
     print *, 'Nearest gridpoints to Montserrat:'
     print *, '   Levels: ', w_in%grid%lev(k1:k2)
     print *, ' Lon, lat: ', w_in%grid%lon(is)-360., w_in%grid%lat(js)

     q_so2 = 1.4262e-06 ! Based on Nick estimated 55K tons of SO2.
     i_so2 = reg_in%i_su + 1
     print *, 'Variable name: ', reg_in%vname(i_so2)

     w_in%q(is,js,k1:k2,i_so2) = q_so2
     
#else

     k1 = 55 - 21 + 1
     k2 = 55 - 16 + 1

     i1 = w_in%grid%im - 82  ! source
     i2 = w_in%grid%im - 67  ! source
     i1r = i1 - 14 ! destination
     i2r = i2 - 14 ! destination

     j1 = 98
     j2 = 107

     allocate ( tmp(i1:i2,j1:j2,k1:k2) )

     print *, 'Nearest gridpoints to feature:'
     print *, '   Levels: ', k1, k2, w_in%grid%lev(k1:k2)
     print *, '     Lats: ', w_in%grid%lat(j1:j2)
     print *, ' Src Lons: ', w_in%grid%lon(i1)-360., w_in%grid%lon(i2)-360.
     print *, ' Dst Lons: ', w_in%grid%lon(i1r)-360., w_in%grid%lon(i2r)-360.

     i_so2 = reg_in%i_su + 1
     print *, 'Variable name: ', reg_in%vname(i_so2)

     tmp = w_in%q(i1r:i2r,j1:j2,k1:k2,i_so2) 
     w_in%q(i1r:i2r,j1:j2,k1:k2,i_so2) = w_in%q(i1:i2,j1:j2,k1:k2,i_so2) 
     w_in%q(i1:i2,j1:j2,k1:k2,i_so2) = tmp
     
#endif


!    Write restart
!    ------------
     call Chem_BundleWrite ( trim(out_file), nymd, nhms, 1, w_in, ier )

     end program edit_crst


