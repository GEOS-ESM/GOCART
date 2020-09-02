#include "unused_dummy.H"
!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  Chem_StateMod --- Replacement for Chem_State under ESMF
!
! !INTERFACE:
!

   module  Chem_StateMod

! !USES:

   use ESMF
   use ESMFL_Mod, only: ESMFL_StateGetPointerToData

   use Chem_Mod           ! Chemistry Base Class
   use mod_diag           ! fvGCM diagnostics
   use m_die
   implicit NONE

!
! !PUBLIIC MEMBER FUNCTIONS:
!

   PUBLIC  Chem_StateSetNeeded
   PUBLIC  Chem_StateGetArray2D
   PUBLIC  Chem_StateGetArray3D

!
! !DESCRIPTION:
!
!  This module provides a few methods for managing data arrays from an
!  ESMF state. This is a replacing for the orifinal, pre-ESMF, Chem_State
!  class.
!
! !REVISION HISTORY:
!
!  27Feb2005 da Silva  First crack.
!
!EOP
!-------------------------------------------------------------------------

   integer, parameter :: XY = 2, XYZ = 3 ! all that is needed for now

!  Holds GEOS-4 variable catalogue (not instantiable by definition)
!  ---------------------------------------------------------------
   logical :: diaglist_init = .false.
   type(diag_type), save :: diag(pdiag)

CONTAINS

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Chem_StateSetNeeded - Activate/deactivate fields in State
!
! !INTERFACE:
!

   subroutine Chem_StateSetNeeded ( This, which, needed, rc )
  
! !USES:

  implicit NONE

! !INPUT/OUTPUT PARAMETERS:

   type(ESMF_State), intent(inout) :: This  ! Import State for Chem_GridComp

! !INPUT PARAMETERS:

   integer, intent(in) :: which  ! which field to activate/deactivate
   logical, intent(in) :: needed ! whether it is needed or not

! !OUTPUT PARAMETERS:

   integer, intent(out) ::  rc                ! Error return code:
                                              !  0 - all is well
                                              !  1 -
 
! !DESCRIPTION: This routine set a field as needed or not.
!
! !REVISION HISTORY:
!
!  18Sep2003 da Silva  First crack.
!
!EOP
!-------------------------------------------------------------------------

!
!  Nothing to do in GEOS-5
!
   _UNUSED_DUMMY(this)
   _UNUSED_DUMMY(which)
   _UNUSED_DUMMY(needed)

   rc = 0

 end subroutine Chem_StateSetNeeded

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Chem_StateGetArray2D
!
! !INTERFACE:
!

   subroutine Chem_StateGetArray2D ( This, which, array, rc )

! !USES:

  implicit NONE

! !INPUT/OUTPUT PARAMETERS:

   type(ESMF_State), intent(inout) :: This  ! Import State for Chem_GridComp

! !INPUT PARAMETERS:

   integer, intent(in) :: which  ! which field to activate/deactivate

! !OUTPUT PARAMETERS:

   real, pointer :: array(:,:)
   integer, intent(out) ::  rc                ! Error return code:
                                              !  0 - all is well
                                              !  1 -
 
! !DESCRIPTION: This routine assign a pointer to a 2D array in the state.
!
! !REVISION HISTORY:
!
!  18Sep2003 da Silva  First crack.
!
!EOP
!-------------------------------------------------------------------------

   character(len=80) :: vname

   call GetVarName_ ( which, vname )

   call ESMFL_StateGetPointerToData ( This, array, vname)
   rc =0


 end subroutine Chem_StateGetArray2D


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Chem_StateGetArray3D
!
! !INTERFACE:
!

   subroutine Chem_StateGetArray3D ( This, which, array, rc )

! !USES:

  implicit NONE

! !INPUT/OUTPUT PARAMETERS:

   type(ESMF_State), intent(inout) :: This  ! Import State for Chem_GridComp

! !INPUT PARAMETERS:

   integer, intent(in) :: which  ! which field to activate/deactivate

! !OUTPUT PARAMETERS:

   real, pointer :: array(:,:,:)
   integer, intent(out) ::  rc                ! Error return code:
                                              !  0 - all is well
                                              !  1 -
 
! !DESCRIPTION: This routine assign a pointer to a 3D array in the state.
!
! !REVISION HISTORY:
!
!  18Sep2003 da Silva  First crack.
!
!EOP
!-------------------------------------------------------------------------

   character(len=80) :: vname

   call GetVarName_ ( which, vname )

   call ESMFL_StateGetPointerToData ( This, array, vname)
   rc =0
   if (which == iAIRDENS) array = 1.2 ! THIS is ugly and VERY WRONG!!!!

 end subroutine Chem_StateGetArray3D

!.................. private auxiliary routines .................

 subroutine GetVarName_ ( which, vname )

   implicit NONE
   integer, intent(in) :: which
   character(len=80), intent(out) :: vname 

!  Make sure GEOS-4 catalogue is initialized
!  -----------------------------------------
   if ( .not. diaglist_init ) then
      call diaglist ( diag )
   end if

   if ( which < 1 .OR. which > pdiag ) &
     call die ( 'Chem_stateMod::GetVarName_', 'invalid variable index' )

   vname = trim ( diag(which)%name )
   
 end subroutine GetVarName_


 end module Chem_StateMod

