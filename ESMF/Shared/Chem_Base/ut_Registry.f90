!
! Unit tester for Chem_Registry.F90
!

   Program ut_Registry

!!!   use m_die, only: die
   use Chem_RegistryMod 

   implicit NONE


   character(len=*), parameter ::  myname = 'ut_Registry'
   type(Chem_Registry) :: reg
   integer ier

!  No rc file name provided
!  ------------------------
   reg = Chem_RegistryCreate ( ier )
   if ( ier /= 0 ) call die ( myname, 'cannot create registry' )
   CALL Chem_RegistryPrint ( reg )
   call Chem_RegistryDestroy ( reg, ier )
   if ( ier /= 0 ) call die ( myname, 'cannot destroy registry' )

!  No rc file name provided
!  ------------------------
   reg = Chem_RegistryCreate ( ier, 'Chem_Registry.rc' )
   if ( ier /= 0 ) call die ( myname, 'cannot create registry' )
   CALL Chem_RegistryPrint ( reg )
   call Chem_RegistryDestroy ( reg, ier )
   if ( ier /= 0 ) call die ( myname, 'cannot destroy registry' )

! ..................................................................

   end Program ut_Registry

   subroutine die(name,msg)
     character(len=*) name, msg
     print *, trim(name)//': '//trim(msg)
     call exit(7)
   end subroutine die
