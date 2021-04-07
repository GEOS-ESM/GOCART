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
   use BC_GridCompMod        ! Black Carbon
   use DU_GridCompMod        ! Dust
   use OC_GridCompMod        ! Organic Carbon
   use BRC_GridCompMod       ! Brown Carbon
   use SS_GridCompMod        ! Sea Salt
   use SU_GridCompMod        ! Sulfates
   use NI_GridCompMod        ! Nitrate
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
        type(DU_GridComp)  :: gcDU
        type(SS_GridComp)  :: gcSS
        type(BC_GridComp)  :: gcBC
        type(OC_GridComp)  :: gcOC
        type(BRC_GridComp) :: gcBRC
        type(SU_GridComp)  :: gcSU
        type(CFC_GridComp) :: gcCFC
        type(Rn_GridComp)  :: gcRn
        type(CH4_GridComp) :: gcCH4
        type(NI_GridComp)  :: gcNI
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
   if ( chemReg%doing_DU ) then
#     include "DU_ExportSpec___.h"
      call DU_GridCompSetServices(GC, chemReg, __RC__)
   end if

   if ( chemReg%doing_CO ) then
#     include "CO_ExportSpec___.h"
      call CO_GridCompSetServices(GC, chemReg, __RC__)
   end if

   if ( chemReg%doing_BC ) then
#     include "BC_ExportSpec___.h"
      call BC_GridCompSetServices(GC, chemReg, __RC__)
   end if

   if ( chemReg%doing_OC ) then
#     include "OC_ExportSpec___.h"
      call OC_GridCompSetServices(GC, chemReg, __RC__)
   end if

   if ( chemReg%doing_BRC ) then
#     include "BRC_ExportSpec___.h"
      call BRC_GridCompSetServices(GC, chemReg, __RC__)
   end if

   if ( chemReg%doing_SS ) then
#     include "SS_ExportSpec___.h"
      call SS_GridCompSetServices(GC, chemReg, __RC__)
   end if

   if ( chemReg%doing_SU ) then
#     include "SU_ExportSpec___.h"
      call SU_GridCompSetServices(GC, chemReg, __RC__)
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

   if ( chemReg%doing_NI ) then
#     include "NI_ExportSpec___.h"
      call NI_GridCompSetServices(GC, chemReg, __RC__)
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

   type(MAPL_MetaComp), pointer :: state

   gcThis%name = 'Composite Constituent Package'


   call MAPL_GetObjectFromGC(gc, state, RC)
   VERIFY_(RC)


!  Require Relative Humidity
!  -------------------------
   call Chem_StateSetNeeded ( impChem, iRELHUM, .true., rc )
   if ( rc /= 0 ) then
      if (MAPL_AM_I_ROOT()) print *, Iam//': failed StateSetNeeded'
      return
   end if

!  Initialize AOD tables
!  ---------------------
   if ( w_c%reg%doing_DU .or. w_c%reg%doing_SS .or. w_c%reg%doing_SU .or. &
        w_c%reg%doing_BC .or. w_c%reg%doing_OC .or. w_c%reg%doing_NI .or. &
        w_C%reg%doing_BRC ) then
        allocate ( gcThis%mie_tables, stat = rc )
        if ( rc /= 0 ) then
           if (MAPL_AM_I_ROOT()) print *, Iam//': cannot allocate Mie tables'
           return
        end if
!       Here we are assuming the mie tables to use
!       In future, pass a registry file
        gcThis%mie_tables = Chem_MieCreate ( 'Aod_Registry.rc', rc )
        if ( rc /= 0 ) then
           if (MAPL_AM_I_ROOT()) print *, Iam//': MieCreate failed ', rc
           return
        end if
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

!  Dust
!  ----
   if ( w_c%reg%doing_DU ) then
      gcThis%gcDU%mie_tables => gcThis%mie_tables

      if (.not. data_driven) then
         call DU_GridCompInitialize ( gcThis%gcDU, w_c, impChem, expChem, state, &
                                      nymd, nhms, cdt, rc )
         if ( rc /= 0 ) then  
            if (MAPL_AM_I_ROOT()) print *, Iam//': DU failed to initialize ', rc
            rc = 3000 + rc
            return
         end if
      end if
   end if

!  Sea Salt
!  --------
   if ( w_c%reg%doing_SS ) then
      gcThis%gcSS%mie_tables => gcThis%mie_tables

      if (.not. data_driven) then
         call SS_GridCompInitialize ( gcThis%gcSS, w_c, impChem, expChem, state, &
                                      nymd, nhms, cdt, rc )
         if ( rc /= 0 ) then  
            if (MAPL_AM_I_ROOT()) print *, Iam//': SS failed to initialize ', rc
            rc = 4000 + rc
            return
         end if
      end if
   end if

!  Black carbon
!  ------------
   if ( w_c%reg%doing_BC ) then
      gcThis%gcBC%mie_tables => gcThis%mie_tables

      if (.not. data_driven) then
         call BC_GridCompInitialize ( gcThis%gcBC, w_c, impChem, expChem, state, &
                                      nymd, nhms, cdt, rc )
         if ( rc /= 0 ) then  
            if (MAPL_AM_I_ROOT()) print *, Iam//': BC failed to initialize ', rc
            rc = 5000 + rc
            return
         end if
      end if
   end if

!  Organic Carbon
!  --------------
   if ( w_c%reg%doing_OC ) then
      gcThis%gcOC%mie_tables => gcThis%mie_tables

      if (.not. data_driven) then
         call OC_GridCompInitialize ( gcThis%gcOC, w_c, impChem, expChem, state, &
                                      nymd, nhms, cdt, rc )
         if ( rc /= 0 ) then  
            if (MAPL_AM_I_ROOT()) print *, Iam//': OC failed to initialize ', rc
            rc = 6000 + rc
            return
         end if
      end if   
   end if

!  Brown Carbon
!  ------------
   if ( w_c%reg%doing_BRC ) then
      gcThis%gcBRC%mie_tables => gcThis%mie_tables

      if (.not. data_driven) then
         call BRC_GridCompInitialize ( gcThis%gcBRC, w_c, impChem, expChem, state, &
                                      nymd, nhms, cdt, rc )
         if ( rc /= 0 ) then  
            if (MAPL_AM_I_ROOT()) print *, Iam//': BRC failed to initialize ', rc
            rc = 6000 + rc
            return
         end if
      end if   
   end if

!  Sulfates
!  --------
   if ( w_c%reg%doing_SU ) then
      gcThis%gcSU%mie_tables => gcThis%mie_tables

      if (.not. data_driven) then
         call SU_GridCompInitialize ( gcThis%gcSU, w_c, impChem, expChem, state, &
                                      nymd, nhms, cdt, rc )
         if ( rc /= 0 ) then  
            if (MAPL_AM_I_ROOT()) print *, Iam//': SU failed to initialize ', rc
            rc = 7000 + rc
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

!  Nitrate
!  -------
   if ( w_c%reg%doing_NI ) then
      gcThis%gcNI%mie_tables => gcThis%mie_tables

      if (.not. data_driven) then
         call NI_GridCompInitialize ( gcThis%gcNI, w_c, impChem, expChem, state, &
                                      nymd, nhms, cdt, rc )
         if ( rc /= 0 ) then  
            if (MAPL_AM_I_ROOT()) print *, Iam//': NI failed to initialize ', rc
            rc = 8900 + rc
            return
         end if
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
   character(len=ESMF_MAXSTR)   :: Iam
 
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

   Iam = 'Aero_GridCompRun1'

   call MAPL_GetObjectFromGC(gc, state, RC)
   VERIFY_(RC)

!  Dust
!  ----
   if ( w_c%reg%doing_DU ) then
      call MAPL_TimerOn(state,"DU")
      call DU_GridCompRun1 ( gcThis%gcDU, w_c, impChem, expChem, state, &
                             nymd, nhms, cdt, rc )
      call MAPL_TimerOff(state,"DU")
      if ( rc /= 0 ) then  
           rc = 3000 + rc
           return
      end if
   end if

!  Sea Salt
!  --------
   if ( w_c%reg%doing_SS ) then
      call MAPL_TimerOn(state,"SS")
      call SS_GridCompRun1 ( gcThis%gcSS, w_c, impChem, expChem, state, &
                             nymd, nhms, cdt, rc )
      call MAPL_TimerOff(state,"SS")
      if ( rc /= 0 ) then  
           rc = 4000 + rc
           return
      end if
   end if

!  Black carbon
!  ------------
   if ( w_c%reg%doing_BC ) then
      call MAPL_TimerOn(state,"BC")
      call BC_GridCompRun1 ( gcThis%gcBC, w_c, impChem, expChem, state, &
                            nymd, nhms, cdt, rc )
      call MAPL_TimerOff(state,"BC")
      if ( rc /= 0 ) then  
           rc = 5000 + rc
           return
      end if
   end if

!  Organic Carbon
!  --------------
   if ( w_c%reg%doing_OC ) then
      call MAPL_TimerOn(state,"OC")
      call OC_GridCompRun1 ( gcThis%gcOC, w_c, impChem, expChem, state, &
                            nymd, nhms, cdt, rc )
      call MAPL_TimerOff(state,"OC")
      if ( rc /= 0 ) then  
           rc = 6000 + rc
           return
      end if
   end if

!  Brown Carbon
!  ------------
   if ( w_c%reg%doing_BRC ) then
      call MAPL_TimerOn(state,"BRC")
      call BRC_GridCompRun1 ( gcThis%gcBRC, w_c, impChem, expChem, state, &
                            nymd, nhms, cdt, rc )
      call MAPL_TimerOff(state,"BRC")
      if ( rc /= 0 ) then  
           rc = 6100 + rc
           return
      end if
   end if

!  Sulfates
!  --------
   if ( w_c%reg%doing_SU ) then
      call MAPL_TimerOn(state,"SU")
      call SU_GridCompRun1 ( gcThis%gcSU, w_c, impChem, expChem, state, &
                            nymd, nhms, cdt, rc )
      call MAPL_TimerOff(state,"SU")
      if ( rc /= 0 ) then  
           rc = 7000 + rc
           return
      end if
   end if

!  Nitrate
!  -------
   if ( w_c%reg%doing_NI ) then
      call MAPL_TimerOn(state,"NI")
      call NI_GridCompRun1 ( gcThis%gcNI, w_c, impChem, expChem, state, &
                             nymd, nhms, cdt, rc )
      call MAPL_TimerOff(state,"NI")
      if ( rc /= 0 ) then  
           rc = 8900 + rc
           return
      end if
   end if


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
   character(len=ESMF_MAXSTR)   :: Iam
 
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
   
   Iam = 'Aero_GridCompRun2'

   call MAPL_GetObjectFromGC(gc, state, RC)
   VERIFY_(RC)

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

!  Dust
!  ----
   if ( w_c%reg%doing_DU ) then
      call MAPL_TimerOn(state,"DU")
      call DU_GridCompRun2 ( gcThis%gcDU, w_c, impChem, expChem, state, &
                             run_alarm, nymd, nhms, cdt, rc )
      call MAPL_TimerOff(state,"DU")
      if ( rc /= 0 ) then  
           rc = 3000 + rc
           return
      end if
   end if

!  Sea Salt
!  --------
   call MAPL_TimerOn(state,"SS")
   if ( w_c%reg%doing_SS ) then
      call SS_GridCompRun2 ( gcThis%gcSS, w_c, impChem, expChem, state, &
                             run_alarm, nymd, nhms, cdt, rc )
      if ( rc /= 0 ) then  
           rc = 4000 + rc
           return
      end if
   end if
   call MAPL_TimerOff(state,"SS")

!  Black carbon
!  ------------
   if ( w_c%reg%doing_BC ) then
      call MAPL_TimerOn(state,"BC")
      call BC_GridCompRun2 ( gcThis%gcBC, w_c, impChem, expChem, state, &
                             run_alarm, nymd, nhms, cdt, rc )
      call MAPL_TimerOff(state,"BC")
      if ( rc /= 0 ) then  
           rc = 5000 + rc
           return
      end if
   end if

!  Organic Carbon
!  --------------
   if ( w_c%reg%doing_OC ) then
      call MAPL_TimerOn(state,"OC")
      call OC_GridCompRun2 ( gcThis%gcOC, w_c, impChem, expChem, state, &
                             run_alarm, nymd, nhms, cdt, rc )
      call MAPL_TimerOff(state,"OC")
      if ( rc /= 0 ) then  
           rc = 6000 + rc
           return
      end if
   end if

!  Brown Carbon
!  ------------
   if ( w_c%reg%doing_BRC ) then
      call MAPL_TimerOn(state,"BRC")
      call BRC_GridCompRun2 ( gcThis%gcBRC, w_c, impChem, expChem, state, &
                             run_alarm, nymd, nhms, cdt, rc )
      call MAPL_TimerOff(state,"BRC")
      if ( rc /= 0 ) then  
           rc = 6100 + rc
           return
      end if
   end if

!  Sulfates
!  --------
   if ( w_c%reg%doing_SU ) then
      call MAPL_TimerOn(state,"SU")
      call SU_GridCompRun2 ( gcThis%gcSU, w_c, impChem, expChem, state, &
                             run_alarm, nymd, nhms, cdt, rc )
      call MAPL_TimerOff(state,"SU")
      if ( rc /= 0 ) then  
           rc = 7000 + rc
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

!  Nitrate
!  -------
   if ( w_c%reg%doing_NI ) then
      call MAPL_TimerOn(state,"NI")
      call NI_GridCompRun2 ( gcThis%gcNI, w_c, impChem, expChem, state, &
                             run_alarm, nymd, nhms, cdt, rc )
      call MAPL_TimerOff(state,"NI")
      if ( rc /= 0 ) then  
           rc = 8900 + rc
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

   subroutine Aero_GridCompFinalize ( gcThis, w_c, gc, impChem, expChem, &
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
   type(ESMF_GridComp), intent(inout) :: gc
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

   TYPE(MAPL_MetaComp), pointer :: state
   character(len=*), parameter  :: Iam = 'Aero_GridCompFinalize'

   call MAPL_GetObjectFromGC(gc, state, RC)
   VERIFY_(RC)
   
 

!  Finalize AOD tables
!  -------------------
   if ( w_c%reg%doing_DU .or. w_c%reg%doing_SS .or. w_c%reg%doing_SU .or. &
        w_c%reg%doing_BC .or. w_c%reg%doing_OC .or. w_c%reg%doing_NI .or. &
        w_c%reg%doing_BRC ) then
        call Chem_MieDestroy ( gcThis%mie_tables, rc )
        if ( rc /= 0 ) return
        deallocate ( gcThis%mie_tables, stat = rc )
        if ( rc /= 0 ) return
   end if


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

!  Dust
!  ----
   if ( w_c%reg%doing_DU ) then
      call DU_GridCompFinalize ( gcThis%gcDU, w_c, impChem, expChem, state, &
                                 nymd, nhms, cdt, rc )
      if ( rc /= 0 ) then  
           rc = 3000 + rc
           return
      end if
   end if

!  Sea Salt
!  --------
   if ( w_c%reg%doing_SS ) then
      call SS_GridCompFinalize ( gcThis%gcSS, w_c, impChem, expChem, state, &
                                 nymd, nhms, cdt, rc )
      if ( rc /= 0 ) then  
           rc = 4000 + rc
           return
      end if
   end if

!  Black carbon
!  ------------
   if ( w_c%reg%doing_BC ) then
      call BC_GridCompFinalize ( gcThis%gcBC, w_c, impChem, expChem, state, &
                                 nymd, nhms, cdt, rc )
      if ( rc /= 0 ) then  
           rc = 5000 + rc
           return
      end if
   end if

!  Organic Carbon
!  --------------
   if ( w_c%reg%doing_OC ) then
      call OC_GridCompFinalize ( gcThis%gcOC, w_c, impChem, expChem, state, &
                                 nymd, nhms, cdt, rc )
      if ( rc /= 0 ) then  
           rc = 6000 + rc
           return
      end if
   end if

!  Brown Carbon
!  ------------
   if ( w_c%reg%doing_BRC ) then
      call BRC_GridCompFinalize ( gcThis%gcBRC, w_c, impChem, expChem, state, &
                                 nymd, nhms, cdt, rc )
      if ( rc /= 0 ) then  
           rc = 6000 + rc
           return
      end if
   end if

!  Sulfates
!  --------
   if ( w_c%reg%doing_SU ) then
      call SU_GridCompFinalize ( gcThis%gcSU, w_c, impChem, expChem, state, &
                                 nymd, nhms, cdt, rc )
      if ( rc /= 0 ) then  
           rc = 7000 + rc
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

!  Nitrate
!  -------
   if ( w_c%reg%doing_NI ) then
      call NI_GridCompFinalize ( gcThis%gcNI, w_c, impChem, expChem, state, &
                                  nymd, nhms, cdt, rc )
      if ( rc /= 0 ) then  
           rc = 8900 + rc
           return
      end if
   end if

   return

 end subroutine Aero_GridCompFinalize

 end module Aero_GridCompMod

