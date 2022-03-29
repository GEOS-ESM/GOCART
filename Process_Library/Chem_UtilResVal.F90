   function Chem_UtilResVal( im_World, jm_World, res_value, rc ) result (val)

! !USES:

   implicit NONE

   real :: val                                ! resolution dependent value

! !INPUT/OUTPUT PARAMETERS:
   integer, intent(in) :: im_World, jm_World  ! number of global grid cells
   real,    intent(in) :: res_value(:)        ! array with the resolution dependent values:
                                              ! the 'a', 'b', ..., 'e' resolution values have
                                              ! indexes 1, 2, ..., 5.

! !OUTPUT PARAMETERS:
   integer, intent(inout) :: rc               ! return code


! !DESCRIPTION:
!
! !REVISION HISTORY:
!
! 13 Feb2012   Anton Darmenov  First crack.
! 25 Oct2012   Anton Darmenov  Added support for FV3 resolutions.
! 19 Aug2020   E. Sherman - moved from Chem_UtilMod.F90 to process library
!
!EOP
!-------------------------------------------------------------------------
       character(len=*), parameter :: Iam = 'Chem_UtilResVal'

       integer            :: i_res
       integer, parameter :: res_a = 1  ! 'a' to 'e' resolution indexes
       integer, parameter :: res_b = 2  !
       integer, parameter :: res_c = 3  !
       integer, parameter :: res_d = 4  !
       integer, parameter :: res_e = 5  !
       integer, parameter :: res_f = 6  !

       i_res = 0

       if ((im_World < 1) .or. (jm_World < 1)) then
!           call die(Iam, 'incorrect model resolution')
           print*,'GOCART2G_Process::Chem_UtilResVal - incorrect model resolution'
           return
       end if

       if (jm_World == 6*im_World) then
           if (im_World <= 24) then
               i_res = res_a
           else if (im_World <=  48) then
               i_res = res_b
           else if (im_World <=  90) then
               i_res = res_c
           else if (im_World <= 180) then
               i_res = res_d
           else if (im_World <= 360) then
               i_res = res_e
           else if (im_World <= 720) then
               i_res = res_f
           else
               i_res = res_f
           end if
       else
           if ((im_World <= 72) .and. (jm_World <= 46)) then
               i_res = res_a
           else if ((im_World <=  144) .and. (jm_World <=  91)) then
               i_res = res_b
           else if ((im_World <=  288) .and. (jm_World <= 181)) then
               i_res = res_c
           else if ((im_World <=  576) .and. (jm_World <= 361)) then
               i_res = res_d
           else if ((im_World <= 1152) .and. (jm_World <= 721)) then
               i_res = res_e
           else if ((im_World <= 2304) .and. (jm_World <=1441)) then
               i_res = res_f
           else
               i_res = res_f
           end if


       end if

       if ((i_res < 1) .or. (i_res > size(res_value))) then
           val = 0.0
           rc  = __FAIL__
       else
           val = res_value(i_res)
           rc  = __SUCCESS__
       end if

   end function Chem_UtilResVal
