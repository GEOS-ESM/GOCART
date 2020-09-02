!
! Unit tester for Chem_Registry.F90
!

   Program ut_Bundle

   use m_die, only: die
   use Chem_RegistryMod 
   use Chem_BundleMod

   implicit NONE


   character(len=*), parameter ::  myname = 'ut_Bundle'
   type(Chem_Registry) :: reg
   type(Chem_Bundle) :: w_c
   integer ier, i, j, k, l, prec, nymd, nhms
   integer :: im = 72, jm = 46, km = 18
   real d2r, factor, lat, lon, dp
   integer, external :: system

!  No rc file name provided
!  ------------------------
   reg = Chem_RegistryCreate ( ier )
   if ( ier /= 0 ) call die ( myname, 'cannot create registry' )
   call reg_print_ ( reg )

!  Create Bundle - No memory alloc
!  -------------------------------
   call Chem_BundleCreate1PE_ ( reg, im, jm, km, w_c, ier )
   if ( ier /= 0 ) call die ( myname, 'cannot create bundle' )

!  Uninitialized bundle
!  --------------------
   print *, 'Uninitialized bundle'
   call Chem_BundleStat ( 6, w_c, ier )
   if ( ier /= 0 ) call die ( myname, 'cannot stat bundle' )

!  Put data in it
!  --------------
   call fill_bundle_()

!  Initialized bundle
!  ------------------
   print *, 'Initialized bundle'
   call Chem_BundleStat ( 6, w_c, ier )
   if ( ier /= 0 ) call die ( myname, 'cannot stat bundle' )

! Write out bundle
! ----------------
  prec = 0
  nymd = 19600205
  nhms = 120000
  ier = system ( "/bin/rm -rf bundle.hdf" )
  call  Chem_BundleWrite ( 'bundle.hdf', nymd, nhms, prec, w_c, ier, &
                           verbose = .true. )
   if ( ier /= 0 ) then
      print *, 'ier = ', ier
      call die ( myname, 'cannot write bundle' )
   end if
   print *, 'ncdumping bundle1...'
   ier = system ( "hdfdump bundle.hdf > bundle1.ncl" )
   ier = system ( "/bin/mv bundle.hdf bundle1.hdf" )

!  Destroy bundle
!  --------------
   call Chem_BundleDestroy ( w_c, ier )
   if ( ier /= 0 ) call die ( myname, 'cannot destroy bundle' )

!  Next read bundle from file
!  --------------------------
   call Chem_BundleRead ( 'bundle1.hdf', nymd, nhms, w_c, ier, &
                          timidx=0, chemReg=reg )   
   if ( ier /= 0 ) call die ( myname, 'cannot read bundle' )

   print *, 'Read bundle'
   call Chem_BundleStat ( 6, w_c, ier )
   if ( ier /= 0 ) call die ( myname, 'cannot stat bundle' )

!  Write out file and diff them
!  ----------------------------
   ier = system ( "/bin/rm -rf bundle.hdf" )
    call  Chem_BundleWrite ( 'bundle.hdf', nymd, nhms, prec, w_c, ier, &
                           verbose = .true. )
   if ( ier /= 0 ) then
      print *, 'ier = ', ier
      call die ( myname, 'cannot write bundle' )
   end if
   print *, 'ncdumping bundle2...'
   ier = system ( "hdfdump bundle.hdf > bundle2.ncl" )
   ier = system ( "/bin/mv bundle.hdf bundle2.hdf" )

   print *, 'differences are...'
    ier = system ( "diff bundle1.ncl bundle2.ncl" )

!  Clean up mess
!  -------------
   call Chem_RegistryDestroy ( reg, ier )
   if ( ier /= 0 ) call die ( myname, 'cannot destroy registry' )
   call Chem_BundleDestroy ( w_c, ier )
   if ( ier /= 0 ) call die ( myname, 'cannot destroy bundle' )


! ..................................................................

   contains

   subroutine reg_print_ ( reg )
   type(Chem_Registry) :: reg

   integer i

   print *
   print *, '-----------------------------------------------------------'
   print *
   print *, 'Total number of tracers: ', reg%nq
!   print *, 'Number of fixed tracers: ', reg%nf
   do i = 1, reg%nq
      print *
      print *, '    Tracer: ', i
      print *, 'Short Name: ', trim(reg%vname(i))
      print *, 'Long  Name: ', trim(reg%vtitle(i))
      print *, '     Units: ', trim(reg%vunits(i))
   end do

   print *

  if ( reg%doing_H2O ) then
   print *, 'Tracer H2O: ', reg%doing_H2O, reg%n_H2O, reg%i_H2O, reg%j_H2O
  end if
  if ( reg%doing_O3 ) then
   print *, 'Tracer  O3: ', reg%doing_O3, reg%n_O3, reg%i_O3, reg%j_O3
  end if
  if ( reg%doing_CO ) then
   print *, 'Tracer  CO: ', reg%doing_CO, reg%n_CO, reg%i_CO, reg%j_CO
  end if
  if ( reg%doing_DU ) then
   print *, 'Tracer  DU: ', reg%doing_DU, reg%n_DU, reg%i_DU, reg%j_DU
  end if
  if ( reg%doing_SS ) then
   print *, 'Tracer  SS: ', reg%doing_SS, reg%n_SS, reg%i_SS, reg%j_SS
  end if
  if ( reg%doing_SU ) then
   print *, 'Tracer  SU: ', reg%doing_SU, reg%n_SU, reg%i_SU, reg%j_SU
  end if
  if ( reg%doing_BC ) then
   print *, 'Tracer  BC: ', reg%doing_BC, reg%n_BC, reg%i_BC, reg%j_BC
  end if
  if ( reg%doing_OC ) then
   print *, 'Tracer  OC: ', reg%doing_OC, reg%n_OC, reg%i_OC, reg%j_OC
  end if
  if ( reg%doing_XX ) then
   print *, 'Tracer  XX: ', reg%doing_XX, reg%n_XX, reg%i_XX, reg%j_XX
  end if

  print *

  end subroutine reg_print_

  subroutine fill_bundle_()

  integer kk

   d2r = 3.1415 / 180.
   w_c%grid%ptop = 1.0
   dp = 100000 / km
   do k = 1, km
      kk = 1 + (5  * k ) / km
      do j = 1, jm
         lat = w_c%grid%lat_min + (j-1) * w_c%grid%lat_del 
         do i = 1, im
            lon = w_c%grid%lon_min + (i-1) * w_c%grid%lon_del 
            factor = (1 - sin(kk*d2r*lat) * cos(d2r*kk*lon) )
            w_c%delp(i,j,k) = dp + factor
            w_c%rh(i,j,k) = 10. + (k-1.) * 80. / ( km - 1. ) + factor
            do l = 1, reg%nq
!               w_c%q(i,j,k,l) = l * k * factor
               w_c%qa(l)%data3d(i,j,k) = l * k * factor
            end do
         end do
      end do
   end do

   print *
   do k = 1, km
      print *, 'k, delp = ', k, minval(w_c%delp(1:im,1:jm,k)), &
                             maxval(w_c%delp(1:im,1:jm,k))
   end do

   print *
   do k = 1, km
      print *, 'k, rh  = ', k, minval(w_c%rh(1:im,1:jm,k)), &
                             maxval(w_c%delp(1:im,1:jm,k))
   end do

 end subroutine fill_bundle_


end Program ut_Bundle
