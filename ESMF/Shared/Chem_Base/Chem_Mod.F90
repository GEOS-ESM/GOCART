!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  Chem_Mod --- Chemistry Base Class
!
! !INTERFACE:
!

   module  Chem_Mod

! !USES:

   use Chem_RegistryMod
   use Chem_ArrayMod
   use Chem_BundleMod
   use Chem_aodMod
   implicit NONE

   PUBLIC  ! All of the above

!
! !DESCRIPTION:
!
!  This module implements a base class for the GMAO Chemistry class. This
!  initial class is intended to serve as a stop gap before an ESMF 
!  implementation is adopted. 
!
!
! !REVISION HISTORY:
!
!  04May2003 da Silva  First crack.
!
!EOP
!-------------------------------------------------------------------------

 end module Chem_Mod

