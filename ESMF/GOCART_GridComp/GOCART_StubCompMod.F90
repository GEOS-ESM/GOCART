!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1     !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: GOCART_GridCompMod - The GOCART Aerosol Grid Component
!
! !INTERFACE:
!
   Module GOCART_GridCompMod
!
! !USES:
!

   use ESMF

   implicit none
   private
!
! !PUBLIC MEMBER FUNCTIONS:

   public SetServices
!
! !DESCRIPTION: 
!
!   This is a stub for the {\tt GOCART} gridded component.
!
!EOP
!-------------------------------------------------------------------------

CONTAINS


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: SetServices --- Sets IRF services for GOCART Grid Component
!
! !INTERFACE:

   subroutine SetServices ( GC, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, optional                  :: RC  ! return code

! !DESCRIPTION: Sets Initialize, Run and Finalize services. 
!
! !REVISION HISTORY:
!
!  25feb2005  da Silva  First crack.
!
!EOP
!-------------------------------------------------------------------------

  if ( present(rc) ) rc = 0
  
  end subroutine SetServices

end Module GOCART_GridCompMod
