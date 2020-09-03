!
! Unit test for Chem_UtilMPread()
!

PROGRAM ut_mpread

#if defined (SPMD)
   use mod_comm, only : gid, mp_init, mp_exit, y_decomp
#endif
   use Chem_UtilMod, only: CHem_UtilMPread

   implicit NONE

   integer, parameter :: im = 288, jm = 181, km = 55, nq = 1
   integer :: i1, i2, ig=0, j1, j2, jg=0
   integer :: jnp=jm, nl=km, jfirst=1, jlast=jm, kfirst=1, klast=km
   integer nymd, nhms
   character(len=255) :: filen

#if !defined (SPMD)
    integer :: gid = 0
#endif

   real, allocatable :: src(:,:)


#if defined(SPMD)
   call mp_init()
   call y_decomp( jnp, nl, jfirst, jlast, kfirst, klast, gid)
   if ( gid .eq. 0 ) then
        print *, 'mod_comm initialized'
        print *, 'gid, jnp, nl, jfirst, jlast, kfirst, klast = '
   end if
#endif

   print *, 'gid: ', gid, jnp, nl, jfirst, jlast, kfirst, klast

   allocate ( src(im,jfirst:jlast) )

!!!   filen = '/share/fvchem/ginoux.DU_src.sfc_288x181.clm.hdf'
   filen = 'ginoux.DU_src.sfc_288x181.clm.hdf'
   nymd = 20021219
   nhms = 0
   i1 = 1
   i2 = im
   ig = 0 
   j1 = jfirst
   j2 = jlast
   jg = 0 
   
   call Chem_UtilMPread ( filen, 'du_src', nymd, nhms, &
                          i1, i2, ig, im, j1, j2, jg, jm, km, &
                          var2d = src )

   print *, 'gid ', gid, ':', minval(src), maxval(src)

   deallocate ( src )

#if defined(SPMD)
   call mp_exit( )
#endif

end Program ut_mpread
