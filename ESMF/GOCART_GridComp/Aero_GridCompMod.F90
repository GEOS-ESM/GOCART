#include "MAPL_Generic.h"
!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  Aero_GridCompMod --- Legacy GOCART GridComponent
!
! !INTERFACE:
!

   module  Aero_GridCompMod

! !USES:

   use ESMF
   use MAPL

   use Chem_Mod              ! Chemistry Base Class
   use Chem_StateMod         ! Chemistry State
   use Chem_MieMod           ! Aerosol LU Tables

   Use Chem_UtilMod, only: pmaxmin

   use O3_GridCompMod        ! Ozone
   use CO_GridCompMod        ! Carbon monoxide
   use CO2_GridCompMod       ! Carbon dioxide
   use CFC_GridCompMod       ! CFCs
   use Rn_GridCompMod        ! Radon
   use CH4_GridCompMod       ! Methane

   implicit none

! !PUBLIC TYPES:
!
   PRIVATE
   PUBLIC  Aero_GridComp       ! The Legacy GOCART Object

!
! !PUBLIC MEMBER FUNCTIONS:
!

   PUBLIC  Aero_GridCompSetServices
   PUBLIC  Aero_GridCompInitialize
   PUBLIC  Aero_GridCompRun1
   PUBLIC  Aero_GridCompRun2
   PUBLIC  Aero_GridCompFinalize

!
! !DESCRIPTION:
!
!  This module implements the (pre-ESMF) GOCART Grid Component. This is
!  a composite component which delegates the real work to its 
!  sub-components.
!
! !REVISION HISTORY:
!
!  16Sep2003 da Silva  First crack.
!  24Jan2004 da Silva  Added expChem/cdt to interfaces.
!  24Mar2005 da Silva  Requires RH and saves it under w_c%rh
!  29Mar2005 da Silva  Initializes AOD LUTs.
!  18Oct2005 da Silva  Added CO2.
!  24Jul2006 da Silva  Adapted from Chem_GridComp.
!
!EOP
!-------------------------------------------------------------------------

  type Aero_GridComp
        character(len=255) :: name
        type(Chem_Mie), pointer :: mie_tables
        type(O3_GridComp)  :: gcO3
        type(CO_GridComp)  :: gcCO
        type(CO2_GridComp) :: gcCO2
        type(CFC_GridComp) :: gcCFC
        type(Rn_GridComp)  :: gcRn
        type(CH4_GridComp) :: gcCH4
  end type Aero_GridComp

CONTAINS

   subroutine Aero_GridCompSetServices(GC,chemReg,rc)

   type(ESMF_GridComp), intent(INOUT) :: GC
   type(Chem_Registry), intent(INOUT) :: chemReg
   integer,             intent(OUT  ) :: rc

   integer             :: status
   character(len=ESMF_MAXSTR) :: Iam

   Iam = "Aero_GridCompSetServices"
 
!  Carbon Monoxide
!  ---------------
   if ( chemReg%doing_CO ) then
#     include "CO_ExportSpec___.h"
      call CO_GridCompSetServices(GC,chemReg, __RC__)
   end if

   if ( chemReg%doing_CO2 ) then
#     include "CO2_ExportSpec___.h"
      call CO2_GridCompSetServices(GC, chemReg, __RC__)
   end if

   if ( chemReg%doing_CFC ) then
#     include "CFC_ExportSpec___.h"
      call CFC_GridCompSetServices(GC, chemReg, __RC__)
   end if

   if ( chemReg%doing_O3 ) then
#       include "O3_ExportSpec___.h"
        call O3_GridCompSetServices(GC, chemReg, __RC__)
   endif

   if ( chemReg%doing_Rn ) then
#     include "Rn_ExportSpec___.h"
      call Rn_GridCompSetServices(GC, chemReg, __RC__)
   endif

   if ( chemReg%doing_CH4 ) then
#     include "CH4_ExportSpec___.h"
      call CH4_GridCompSetServices(GC, chemReg, __RC__)
   endif

   RETURN_(ESMF_SUCCESS)

   end subroutine Aero_GridCompSetServices

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Aero_GridCompInitialize --- Initialize Aero_GridComp
!
! !INTERFACE:
!

   subroutine Aero_GridCompInitialize ( gcThis, w_c, gc, impChem, expChem, &
                                        nymd, nhms, cdt, data_driven, rc )

! !USES:

  implicit NONE

! !INPUT PARAMETERS:

   type(Chem_Bundle), intent(inout) :: w_c        ! Chemical tracer fields      
   type(ESMF_GridComp), intent(inout) :: gc
   integer, intent(in) :: nymd, nhms              ! time
   real, intent(in)    :: cdt                     ! chemistry timestep (secs)
   logical, intent(in) :: data_driven             ! GOCART data instance flag


! !OUTPUT PARAMETERS:

   type(Aero_GridComp), intent(out) :: gcThis     ! Grid Component
   type(ESMF_State), intent(inout)  :: impChem    ! Import State
   type(ESMF_State), intent(inout)  :: expChem    ! Export State
   integer, intent(out) ::  rc                    ! Error return code:
                                                  !  0 - all is well
                                                  !  1 -

! !DESCRIPTION: Initializes the GOCART Grid Component. It primarily sets
!               the import state for each active constituent package.
!
! !REVISION HISTORY:
!
!  18Sep2003 da Silva  First crack.
!
!EOP
!-------------------------------------------------------------------------

   character(len=*), parameter :: Iam = 'Aero_GridCompInit'

   gcThis%name = 'Composite Constituent Package'


!  Require Relative Humidity
!  -------------------------
   call Chem_StateSetNeeded ( impChem, iRELHUM, .true., rc )
   if ( rc /= 0 ) then
      if (MAPL_AM_I_ROOT()) print *, Iam//': failed StateSetNeeded'
      return
   end if

!  Ozone & friends
!  ---------------
   if ( w_c%reg%doing_O3 ) then
      call O3_GridCompInitialize ( gcThis%gcO3, w_c, impChem, expChem, &
                                   nymd, nhms, cdt, rc )
      if ( rc /= 0 ) then  
           if (MAPL_AM_I_ROOT()) print *, Iam//': O3 failed to initialize ', rc
           rc = 1000 + rc
           return
      end if
   end if

!  Carbon Monoxide
!  ---------------
   if ( w_c%reg%doing_CO ) then
      if (.not. data_driven) then
         call CO_GridCompInitialize ( gcThis%gcCO, w_c, impChem, expChem, &
                                      nymd, nhms, cdt, rc )
         if ( rc /= 0 ) then  
            if (MAPL_AM_I_ROOT()) print *, Iam//': CO failed to initialize ', rc
            rc = 2000 + rc
            return
         end if
      end if
   end if !

!  Carbon Dioxide
!  ---------------
   if ( w_c%reg%doing_CO2 ) then
      if (.not. data_driven) then
         call CO2_GridCompInitialize ( gcThis%gcCO2, w_c, impChem, expChem, &
                                       nymd, nhms, cdt, rc )
         if ( rc /= 0 ) then  
            if (MAPL_AM_I_ROOT()) print *, Iam//': CO2 failed to initialize ', rc
            rc = 2500 + rc
            return
         end if
      end if   
   end if

!  CFCs
!  ----
   if ( w_c%reg%doing_CFC ) then
      call CFC_GridCompInitialize ( gcThis%gcCFC, w_c, impChem, expChem, &
                                    nymd, nhms, cdt, rc )
      if ( rc /= 0 ) then  
           if (MAPL_AM_I_ROOT()) print *, Iam//': CFC failed to initialize ', rc
           rc = 8000 + rc
           return
      end if
   end if

!  Radon
!  -----
   if ( w_c%reg%doing_Rn ) then
      call Rn_GridCompInitialize ( gcThis%gcRn, w_c, impChem, expChem, &
                                    nymd, nhms, cdt, rc )
      if ( rc /= 0 ) then  
           if (MAPL_AM_I_ROOT()) print *, Iam//': Rn failed to initialize ', rc
           rc = 8500 + rc
           return
      end if
   end if

!  Methane
!  -------
   if ( w_c%reg%doing_CH4 ) then
      call CH4_GridCompInitialize ( gcThis%gcCH4, w_c, impChem, expChem, &
                                    nymd, nhms, cdt, rc )
      if ( rc /= 0 ) then  
           if (MAPL_AM_I_ROOT()) print *, Iam//': CH4 failed to initialize ', rc
           rc = 8800 + rc
           return
      end if
   end if


   call print_init_()

   return

CONTAINS

   subroutine print_init_()

   integer :: i1, i2, j1, j2, ijl, km, n
   real :: qmin, qmax

   i1 = w_c%grid%i1; i2 = w_c%grid%i2
   j1 = w_c%grid%j1; j2 = w_c%grid%j2
   km = w_c%grid%km
   ijl  = ( i2 - i1 + 1 ) * ( j2 - j1 + 1 )

#ifdef DEBUG
   do n = w_c%reg%i_GOCART, w_c%reg%j_GOCART
      call pmaxmin('Init::'//trim(w_c%reg%vname(n)), &
                   w_c%qa(n)%data3d(i1:i2,j1:j2,1:km), qmin, qmax, &
                   ijl, km, 1. )
   end do
#endif

   end subroutine print_init_

   end subroutine Aero_GridCompInitialize




!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Aero_GridCompRun1 --- The GOCART Driver 
!
! !INTERFACE:
!

   subroutine Aero_GridCompRun1 ( gcThis, w_c, gc, impChem, expChem, &
                                  nymd, nhms, cdt, rc )

! !USES:

  implicit NONE

! !INPUT/OUTPUT PARAMETERS:

   type(Aero_GridComp), intent(inout) :: gcThis   ! Grid Component
   type(Chem_Bundle), intent(inout)   :: w_c      ! Chemical tracer fields   
   type(ESMF_GridComp), intent(inout) :: gc

! !INPUT PARAMETERS:

   type(ESMF_State), intent(inout) :: impChem  ! Import State
   integer, intent(in) :: nymd, nhms           ! time
   real, intent(in)    :: cdt                  ! chemistry timestep (secs)


! !OUTPUT PARAMETERS:

   type(ESMF_State), intent(inout) :: expChem    ! Export State
   integer, intent(out) :: rc                    ! Error return code:
                                                 !  0 - all is well
                                                 !  1 -
   type(MAPL_MetaComp), pointer :: state
 
! !DESCRIPTION: This routine implements the so-called GOCART Driver. That 
!               is, adds chemical tendencies to each of the constituents,
!  Note: water wapor, the first constituent is not considered a chemical
!  constituents. 
!
! !REVISION HISTORY:
!
!  18Sep2003 da Silva  First crack.
!  19May2005 da Silva  Phased execution option.
!
!EOP
!-------------------------------------------------------------------------

   call MAPL_GetObjectFromGC(gc,state,rc)


   return

 end subroutine Aero_GridCompRun1



!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Aero_GridCompRun --- The GOCART Driver 
!
! !INTERFACE:
!

   subroutine Aero_GridCompRun2 ( gcThis, w_c, gc, impChem, expChem, &
                                  run_alarm, nymd, nhms, cdt, rc )

! !USES:

  implicit NONE

! !INPUT/OUTPUT PARAMETERS:

   type(Aero_GridComp), intent(inout) :: gcThis   ! Grid Component
   type(Chem_Bundle), intent(inout)   :: w_c      ! Chemical tracer fields   
   type(ESMF_GridComp), intent(inout) :: gc

! !INPUT PARAMETERS:

   type(ESMF_State), intent(inout) :: impChem  ! Import State
   logical             :: run_alarm            ! run alarm 
   integer, intent(in) :: nymd, nhms           ! time
   real, intent(in)    :: cdt                  ! chemistry timestep (secs)


! !OUTPUT PARAMETERS:

   type(ESMF_State), intent(inout) :: expChem    ! Export State
   integer, intent(out) :: rc                    ! Error return code:
                                                 !  0 - all is well
                                                 !  1 -
   type(MAPL_MetaComp), pointer :: state
 
! !DESCRIPTION: This routine implements the so-called GOCART Driver. That 
!               is, adds chemical tendencies to each of the constituents,
!  Note: water wapor, the first constituent is not considered a chemical
!  constituents. 
!
! !REVISION HISTORY:
!
!  18Sep2003 da Silva  First crack.
!  19May2005 da Silva  Phased execution option.
!
!EOP
!-------------------------------------------------------------------------

   call MAPL_GetObjectFromGC(gc,state)
!  Ozone & friends
!  ---------------
   if ( w_c%reg%doing_O3 ) then
      call MAPL_TimerOn(state,"O3")
      call O3_GridCompRun( gcThis%gcO3, w_c, impChem, expChem, &
                           nymd, nhms, cdt, rc )
      call MAPL_TimerOff(state,"O3")
      if ( rc /= 0 ) then  
           rc = 1000 + rc
           return
      end if
   end if

!  Carbon Monoxide
!  ---------------
   if ( w_c%reg%doing_CO ) then
      call MAPL_TimerOn(state,"CO")
      call CO_GridCompRun ( gcThis%gcCO, w_c, impChem, expChem, &
                            nymd, nhms, cdt, rc )
      call MAPL_TimerOff(state,"CO")
       if ( rc /= 0 ) then  
           rc = 2000 + rc
           return
      end if
   end if

!  Carbon Dioxide
!  --------------
   if ( w_c%reg%doing_CO2 ) then
      call MAPL_TimerOn(state,"CO2")
      call CO2_GridCompRun ( gcThis%gcCO2, w_c, impChem, expChem, &
                             nymd, nhms, cdt, rc )
      call MAPL_TimerOff(state,"CO2")
      if ( rc /= 0 ) then  
           rc = 2500 + rc
           return
      end if
   end if

!  CFCs
!  ----
   if ( w_c%reg%doing_CFC ) then
      call MAPL_TimerOn(state,"CFC")
      call CFC_GridCompRun ( gcThis%gcCFC, w_c, impChem, expChem, &
                             nymd, nhms, cdt, rc )
      call MAPL_TimerOff(state,"CFC")
      if ( rc /= 0 ) then  
           rc = 8000 + rc
           return
      end if
   end if

!  Radon
!  -----
   if ( w_c%reg%doing_Rn ) then
      call Rn_GridCompRun ( gcThis%gcRn, w_c, impChem, expChem, &
                            nymd, nhms, cdt, rc )
      if ( rc /= 0 ) then  
           rc = 8500 + rc
           return
      end if
   end if

!  Methane
!  -------
   if ( w_c%reg%doing_CH4 ) then
      call MAPL_TimerOn(state,"CH4")
      call CH4_GridCompRun ( gcThis%gcCH4, w_c, impChem, expChem, &
                             nymd, nhms, cdt, rc )
      call MAPL_TimerOff(state,"CH4")
      if ( rc /= 0 ) then  
           rc = 8800 + rc
           return
      end if
   end if

   return

 end subroutine Aero_GridCompRun2

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Aero_GridCompFinalize --- Finalizes the GOCART Component
!
! !INTERFACE:
!

   subroutine Aero_GridCompFinalize ( gcThis, w_c, impChem, expChem, &
                                      nymd, nhms, cdt, rc )

! !USES:

  implicit NONE

! !INPUT/OUTPUT PARAMETERS:

   type(Aero_GridComp), intent(inout) :: gcThis   ! Grid Component

! !INPUT PARAMETERS:

   type(Chem_Bundle), intent(in)  :: w_c      ! Chemical tracer fields   
   integer, intent(in) :: nymd, nhms          ! time
   real, intent(in) :: cdt                    ! chemistry timestep (secs)

! !OUTPUT PARAMETERS:

   type(ESMF_State), intent(inout) :: impChem ! Import State
   type(ESMF_State), intent(inout) :: expChem   ! Export State
   integer, intent(out) ::  rc                  ! Error return code:
                                                !  0 - all is well
                                                !  1 -
 
! !DESCRIPTION: This routine finalizes this Grid Component.
!
! !REVISION HISTORY:
!
!  18Sep2003 da Silva  First crack.
!
!EOP
!-------------------------------------------------------------------------


!  Ozone & friends
!  ---------------
   if ( w_c%reg%doing_O3 ) then
      call O3_GridCompFinalize ( gcThis%gcO3, w_c, impChem, expChem, &
                                 nymd, nhms, cdt, rc )
      if ( rc /= 0 ) then  
           rc = 1000 + rc
           return
      end if
   end if

!  Carbon Monoxide
!  ---------------
   if ( w_c%reg%doing_CO ) then
      call CO_GridCompFinalize ( gcThis%gcCO, w_c, impChem, expChem, &
                                 nymd, nhms, cdt, rc )
      if ( rc /= 0 ) then  
           rc = 2000 + rc
           return
      end if
   end if

!  Carbon Dioxide
!  --------------
   if ( w_c%reg%doing_CO2 ) then
      call CO2_GridCompFinalize ( gcThis%gcCO2, w_c, impChem, expChem, &
                                  nymd, nhms, cdt, rc )
      if ( rc /= 0 ) then  
           rc = 2500 + rc
           return
      end if
   end if

!  CFCs
!  ----
   if ( w_c%reg%doing_CFC ) then
      call CFC_GridCompFinalize ( gcThis%gcCFC, w_c, impChem, expChem, &
                                  nymd, nhms, cdt, rc )
      if ( rc /= 0 ) then  
           rc = 8000 + rc
           return
      end if
   end if

!  Radon
!  -----
   if ( w_c%reg%doing_Rn ) then
      call Rn_GridCompFinalize ( gcThis%gcRn, w_c, impChem, expChem, &
                                  nymd, nhms, cdt, rc )
      if ( rc /= 0 ) then  
           rc = 8500 + rc
           return
      end if
   end if

!  Methane
!  -------
   if ( w_c%reg%doing_CH4 ) then
      call CH4_GridCompFinalize ( gcThis%gcCH4, w_c, impChem, expChem, &
                                  nymd, nhms, cdt, rc )
      if ( rc /= 0 ) then  
           rc = 8800 + rc
           return
      end if
   end if


   return

 end subroutine Aero_GridCompFinalize

 end module Aero_GridCompMod

