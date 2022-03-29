   function Chem_UtilIdow(nymd) result (idow)
     implicit NONE
     integer, intent(in) :: nymd
     integer :: idow ! day of the week: Sun=1, Mon=2, etc.
     integer :: y, m, d
     integer, parameter :: t(0:11) = (/ 0, 3, 2, 5, 0, 3, 5, 1, 4, 6, 2, 4 /)
     y = nymd / 10000
     m = (nymd - y*10000)/100
     d = nymd - (y*10000 + m*100)
     if ( m<3 ) then
        y = y - 1
     end if
     idow = 1+mod(y + y/4 - y/100 + y/400 + t(m-1) + d,7)
     return
   end function Chem_UtilIdow
