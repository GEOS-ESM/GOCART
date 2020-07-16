   
!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!
#include "MAPL_Generic.h"
   
MODULE Reciever_GridCompMod
!
! !USES:
!
      use ESMF
      use MAPL

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:

      PUBLIC SetServices

      include "mpif.h"

   contains

!BOP

! !IROUTINE: SetServices -- Sets ESMF services for this component

! !INTERFACE:

      subroutine SetServices ( GC, RC )

! !ARGUMENTS:

         type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
         integer,             intent(  OUT) :: RC  ! return code

! !DESCRIPTION:  The SetServices registers the Radiation component

!EOP

!=============================================================================
!
! ErrLog Variables

         character(len=ESMF_MAXSTR)              :: IAm
         integer                                 :: STATUS
         character(len=ESMF_MAXSTR)              :: COMP_NAME

         type(ESMF_Config)          :: cf

!=============================================================================

! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------

         Iam = 'SetServices'
         call ESMF_GridCompGet( GC, NAME=COMP_NAME, CONFIG=CF, RC=STATUS )
         _VERIFY(STATUS)
         Iam = trim(COMP_NAME) // "::" // Iam

!   Set the Initialize, Run, Finalize entry points
!   ----------------------------------------------
         call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_INITIALIZE,  Initialize_, RC=status)
         _VERIFY(STATUS)
         call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,   Run_, RC=status)
         _VERIFY(STATUS)


         call MAPL_AddImportSpec(GC,&
               short_name='var1', &
               long_name='na' , &
               units = 'na', &
               dims = MAPL_DimsHorzOnly, &
               vlocation = MAPL_VLocationNone, rc=status)

!   Generic Set Services
!   --------------------
         call MAPL_GenericSetServices ( GC, rc=status)
         _VERIFY(STATUS)

         _RETURN(ESMF_SUCCESS)

      end subroutine SetServices

!BOP
!
! !IROUTINE:  Initialize_ --- Initialize RUT
!
! !INTERFACE:
!

      SUBROUTINE Initialize_ ( GC, IMPORT, EXPORT, CLOCK, rc )

! !USES:

         implicit NONE

! !INPUT PARAMETERS:

         type(ESMF_Clock),  intent(inout) :: CLOCK     ! The clock

! !OUTPUT PARAMETERS:

         type(ESMF_GridComp), intent(inout) :: GC      ! Grid Component
         type(ESMF_State), intent(inout) :: IMPORT     ! Import State
         type(ESMF_State), intent(inout) :: EXPORT     ! Export State
         integer, intent(out)            :: rc         ! Error return code:
!  0 - all is well
!  1 - 

! !DESCRIPTION: This is a simple ESMF wrapper.
!
! !REVISION HISTORY:
!
!  29 Sep 2011  Anton Darmenov   Cloned from the ExtData GC code
!
!EOP
!-------------------------------------------------------------------------


         type(ESMF_Config)           :: CF          ! Universal Config 
         character(len=ESMF_MAXSTR)  :: Iam
         integer                     :: status
         character(len=ESMF_MAXSTR)  :: comp_name

!  Get my name and set-up traceback handle
!  ---------------------------------------
         Iam = "Initialize_"
         call ESMF_GridCompGet( GC, name=comp_name,rc=status)
         _VERIFY(status)
         Iam = trim(comp_name) // '::' // trim(Iam)

         call MAPL_GridCreate(GC, rc=status)
         _VERIFY(STATUS)

!  Initialize MAPL Generic
!  -----------------------
         call MAPL_GenericInitialize ( GC, IMPORT, EXPORT, clock, RC=status)
         _VERIFY(STATUS)

         _RETURN(ESMF_SUCCESS)

      END SUBROUTINE Initialize_

      SUBROUTINE Run_ ( GC, IMPORT, EXPORT, CLOCK, rc )

! !USES:

         implicit NONE

! !INPUT PARAMETERS:

         type(ESMF_Clock),  intent(inout) :: CLOCK     ! The clock

! !OUTPUT PARAMETERS:

         type(ESMF_GridComp), intent(inout)  :: GC     ! Grid Component
         type(ESMF_State), intent(inout) :: IMPORT     ! Import State
         type(ESMF_State), intent(inout) :: EXPORT     ! Export State
         integer, intent(out) ::  rc                   ! Error return code:

         character(len=ESMF_MAXSTR)    :: Iam
         integer                       :: STATUS
         character(len=ESMF_MAXSTR)    :: comp_name
         real, pointer :: ptr2d(:,:)
         integer :: localPet
         type(ESMF_VM) :: vm

!  Get my name and set-up traceback handle
!  ---------------------------------------
         Iam = "Run_"
         call ESMF_GridCompGet( GC, name=comp_name, rc=status )
         _VERIFY(STATUS)
         Iam = trim(comp_name) // '::' // trim(Iam)

         call ESMF_VMGetCurrent(vm,rc=status)
         _VERIFY(status)
         call ESMF_VMGet(vm,localPet=localPet,rc=status)
         _VERIFY(status)
         call MAPL_GetPointer(import,ptr2d,'var1',rc=status)
         _VERIFY(status)
         print*,'The value of the pointer is ',minval(ptr2d),maxval(ptr2d)

!  --------
         _RETURN(ESMF_SUCCESS)

      END SUBROUTINE Run_

end module Reciever_GridCompMod

