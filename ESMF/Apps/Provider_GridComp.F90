   
!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!
#include "MAPL_Generic.h"
   
MODULE Provider_GridCompMod
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

         print*, "Provider start Set Services"

         Iam = 'SetServices'
         call ESMF_GridCompGet( GC, NAME=COMP_NAME, CONFIG=CF, RC=STATUS )
         _VERIFY(STATUS)
         Iam = trim(COMP_NAME) // "::" // Iam

!   Set the Initialize, Run, Finalize entry points
!   ----------------------------------------------
         print*, "Provider set initialize"
         call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_INITIALIZE,  Initialize_, RC=status)
         _VERIFY(STATUS)
         print*, "Provider set run"
         call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,   Run_, RC=status)
         _VERIFY(STATUS)


         print*, "Provider set export"
         call MAPL_AddExportSpec(GC,&
               short_name='var1', &
               long_name='na' , &
               units = 'na', &
               dims = MAPL_DimsHorzOnly, &
               vlocation = MAPL_VLocationNone, rc=status)

!   Generic Set Services
!   --------------------
         print*, "Provider set Generic Set Services"
         call MAPL_GenericSetServices ( GC, rc=status)
         _VERIFY(STATUS)

         print*, "Provider finish Set Services"

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

         type(ESMF_Grid) :: grid

!  Get my name and set-up traceback handle
!  ---------------------------------------

         print*, "Provider start Initialize"

         call ESMF_GridCompGet(GC,name=comp_name,rc=status)
         _VERIFY(status)
         Iam = "Initialize_"
         Iam = trim(comp_name) // '::' // trim(Iam)

         print*, "Provider MAPL_GridCreate"
         call MAPL_GridCreate(GC, rc=status)
         _VERIFY(STATUS)

!  Initialize MAPL Generic
!  -----------------------
         print*, "Provider Generic initialize"
         call MAPL_GenericInitialize ( GC, IMPORT, EXPORT, clock, RC=status)
         _VERIFY(STATUS)
         print*, "Provider Force allocate"
         call ForceAllocation(Export,rc=status)
         _VERIFY(STATUS)

         print*, "Provider finish Initialize"

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

!  Get my name and set-up traceback handle
!  ---------------------------------------

         print*, "Provider start Run"

         Iam = "Run_"
         call ESMF_GridCompGet( GC, name=comp_name, rc=status )
         _VERIFY(STATUS)
         Iam = trim(comp_name) // '::' // trim(Iam)

         print*,"Provider set export value"
         call MAPL_GetPointer(export,ptr2d,'var1',rc=status)
         _VERIFY(status)
         ptr2d = ptr2d + 1.0

         print*, "Provider finsh Run"

         _RETURN(ESMF_SUCCESS)

      END SUBROUTINE Run_

      subroutine ForceAllocation(state,rc)
         type(ESMF_State), intent(inout) :: state
         integer, optional, intent(out) :: rc
       
         integer :: status
         character(len=*), parameter :: Iam=__FILE__//"::ForceAllocation"
  
         real, pointer :: ptr3d(:,:,:)
         real, pointer :: ptr2d(:,:)
         integer       :: ii
         integer :: itemcount,dims
         character(len=ESMF_MAXSTR), allocatable :: NameList(:)
         type (ESMF_StateItem_Flag), allocatable :: itemTypeList(:)
         type(ESMF_Field) :: Field

         call ESMF_StateGet(State,itemcount=itemCount,__RC__)
         allocate(NameList(itemCount),stat=status)
         _VERIFY(status)
         allocate(itemTypeList(itemCount),stat=status)
         _VERIFY(status)
         call ESMF_StateGet(State,itemNameList=NameList,itemTypeList=itemTypeList,__RC__)
         if (itemCount == 0) then
            _RETURN(ESMF_SUCCESS)
         end if
         do ii=1,itemCount
            if (itemTypeList(ii)==ESMF_STATEITEM_FIELD) then
               call ESMF_StateGet(State,trim(nameList(ii)),field,__RC__)
               call ESMF_AttributeGet(field,name='DIMS',value=dims,__RC__)
               if (dims==MAPL_DimsHorzOnly) then
                  call MAPL_GetPointer(state,ptr2d,trim(nameList(ii)),alloc=.true.,__RC__)
                  ptr2d = 17.0
               else if (dims==MAPL_DimsHorzVert) then
                  call MAPL_GetPointer(state,ptr3d,trim(nameList(ii)),alloc=.true.,__RC__)
               end if
            end if
         enddo
         _RETURN(ESMF_SUCCESS)

      end subroutine ForceAllocation

end module Provider_GridCompMod

