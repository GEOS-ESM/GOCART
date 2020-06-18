  program ut_diurnal_bb

     use Chem_UtilMod
     implicit NONE

!     integer, parameter :: im = 720, jm=361
!     integer, parameter :: im = 288, jm=181
     integer, parameter :: im = 144, jm=91
     integer :: nhms, i, j, hh, mm, n
     real :: dlon, dlat, cdt
     real :: bb(im,jm), bb_(im,jm), lons(im), lats(jm)

     bb_ = 1.
     dlon = 360. / im
     dlat = 180. / (jm-1)

     print *, 'dlon = ', dlon
     print *, 'dlat = ', dlat

     do i = 1, im
        lons(i) = -180. + (i-1) * dlon
     end do
     do j = 1, jm
        lats(j) = -90. + (j-1) * dlat
     end do

     open(20,file='gfed2.bin',form='unformatted',status='old')
     read(20) bb_
     close(20)

     n = 0
     cdt = 30 * 60. ! 15 min in secs

     open(10,file='diurnal2.bin',form='unformatted')
     do hh = 0, 23
        do mm = 0, 30, 30
           n = n + 1
           nhms = hh * 10000 + mm*100
           print *, 'nhms = ', nhms, ' --- bb = ', minval(bb), maxval(bb), n
           call Chem_BiomassDiurnal ( bb, bb_, lons, lats, nhms, cdt)
           write(10) bb
        end do
     end do
     close(10)

   end program ut_diurnal_bb
