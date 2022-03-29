   function Chem_UtilCdow(nymd) result (cdow)
     implicit NONE
     integer, intent(in) :: nymd
     character(len=3) :: cdow ! day of the week: Sun, Mon, etc.
     character(len=3) :: cday(7) = (/ 'Sun','Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat' /)
     cdow = cday(Chem_UtilIdow(nymd))
     return
   end function Chem_UtilCdow
