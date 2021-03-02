#include "MAPL_Generic.h"

!!! TO DO: Please revise Prologues!!!!

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  CO_GridCompMod --- CO Grid Component Class
!
! !INTERFACE:
!

   MODULE  CO_GridCompMod

! !USES:

   USE ESMF
   USE MAPL

   USE Chem_Mod                        ! Chemistry Base Class
   USE Chem_StateMod                   ! Chemistry State
   USE Chem_ConstMod, only: grav
   USE Chem_UtilMod                    ! I/O

   USE m_inpak90                       ! Resource file management
   USE m_die, ONLY: die
   USE Henrys_law_ConstantsMod, ONLY: get_HenrysLawCts

   IMPLICIT NONE

! !PUBLIC TYPES:
!
   PRIVATE
   PUBLIC  CO_GridComp                 ! Multiple instance CO object 
   PUBLIC  CO_GridComp1                ! Single instance CO object

!
! !PUBLIC MEMBER FUNCTIONS:
!

   PUBLIC  CO_GridCompSetServices
   PUBLIC  CO_GridCompInitialize
   PUBLIC  CO_GridCompRun
   PUBLIC  CO_GridCompFinalize

!
! !DESCRIPTION:
!
!  This module implements the (pre-ESMF) CO Grid Component. 
!
! !REVISION HISTORY:
!
!  16Sep2003 da Silva  First crack.
!  31May2005  Nielsen  Mods for 7 CO bins, 5 region masks
!  31May2005 da Silva  Seperate file for biomass emissions; option for
!                       daily templatable files
!  31May2005 da Silva  Moved reading of region mask to init, specified 
!                      fixed time.
!  17Oct2005     Bian  add biogenic emission and CH4 oxidation, two options 
!                      for updating emissions 
!  19dec2005 da Silva  Activated 3D diags for output
!  14Apr2006     Bian  Add CO tagged to fossil fuel, biofuel, biomass burning
!                      and biogenic
!    Oct2006     Bian  Evaluate total and tagged CO performace in GEOS4 system 
!                      with emissions and oxident fields described in 
!                      Bian et al., [2007]. The observations included GMD ground
!                      surface and aircraft measurements, TRACE-P aircraft 
!                      measurements, and satellite MOPITT and AIRS retrieves. 
!  01Aug2006 da Silva  Extensions for GEOS-5.
!  10Mar2008 da Silva  Multiple instances for ARCTAS.
!  18Mar2011  Nielsen  Simplified PBL partitioning for biomass burning emissions  
!  12May2015 Thompson  Bring units into state (Fix for gfortran)
!
!EOP
!-------------------------------------------------------------------------

  TYPE CO_GridComp1

        CHARACTER(LEN=255) :: name            ! generic name of the package
        CHARACTER(LEN=255) :: iname           ! instance name
        CHARACTER(LEN=255) :: rcfilen         ! resource file name

        INTEGER :: instance                   ! instance number

        REAL, POINTER :: eCO_bioburn_(:,:)   ! molec/cm2/s  (before diurnal)
        REAL, POINTER :: eCO_bioburn(:,:)    ! molec/cm2/s
        REAL, POINTER :: eCO_biofuel(:,:)    ! molec/cm2/s
        REAL, POINTER :: eCO_fosfuel(:,:)    ! molec/cm2/s
        REAL, POINTER ::     eCO_iso(:,:)    ! mgC/m2/s, Earth surface
        REAL, POINTER ::     eCO_mon(:,:)    ! mgC/m2/s, Earth surface
        REAL, POINTER ::     eCO_mtn(:,:)    ! mgC/m2/s, Earth surface

        REAL, POINTER ::       CH4(:,:,:)    ! CH4 mixing ratio (ppb)
        REAL, POINTER ::      OHnd(:,:,:)    ! OH number density (#/cm3 or mol mol-1)

        REAL, POINTER ::   COsfcFlux(:,:)    ! CO surface flux kg m^-2 s^-1

        LOGICAL :: DBG                       ! Run-time debug switch
        LOGICAL :: doingBB = .true.          ! Switch to consider biomass burning
        CHARACTER(LEN=ESMF_MAXSTR) :: units_oh ! Units for OH 
        CHARACTER(LEN=ESMF_MAXSTR) :: units_ff ! Units for fossil fuel
        CHARACTER(LEN=ESMF_MAXSTR) :: units_bf ! Units for biofuel

  END TYPE CO_GridComp1

  TYPE CO_GridComp
     INTEGER                     ::  n        ! number of instances 
     TYPE(CO_GridComp1), POINTER ::  gcs(:)   ! instances
  END TYPE CO_GridComp

  REAL, PARAMETER :: radToDeg = 57.2957795

CONTAINS

   subroutine CO_GridCompSetServices(  gc, chemReg, rc)
   type(ESMF_GridComp), intent(INOUT) :: GC
   type(Chem_Registry), intent(INOUT) :: chemReg
   integer,             intent(OUT  ) :: rc

   CHARACTER(LEN=255) :: rcbasen = 'CO_GridComp'
   CHARACTER(LEN=255) :: name

   integer            :: status
   character(len=ESMF_MAXSTR) :: Iam

   integer :: ier,n,i

   type(ESMF_Config) :: cfg

   Iam = "CO_GridCompSetServices"

!  Load resource file
!  ------------------
   cfg = ESMF_ConfigCreate(rc=status)
   VERIFY_(STATUS)
   call ESMF_ConfigLoadFile(cfg,TRIM(rcbasen)//'.rc',rc=status)
   VERIFY_(STATUS)

!  Parse resource file
!  -------------------
   n = ESMF_ConfigGetLen(cfg,label='CO_instances:',rc=status)
   VERIFY_(STATUS)

!  We cannot have fewer instances than the number of
!   CO bins in the registry (it is OK to have less, though)
!  --------------------------------------------------------
   if ( n .LT. chemReg%n_CO ) then
        rc = 35
        return
   else if ( n .GT. chemReg%n_CO ) then
        if (MAPL_AM_I_ROOT()) &
        print *, trim(Iam)// &
                 ': fewer CO bins than possible CO instances: ',&
                 n, chemReg%n_CO
   end if
   n = min(n,chemReg%n_CO )

!  Record name of each instance
!  ----------------------------
   call ESMF_ConfigFindLabel(cfg,'CO_instances:',rc=status)
   VERIFY_(STATUS)
   do i = 1, n
      call ESMF_ConfigGetAttribute(cfg,name,rc=status)
      VERIFY_(STATUS)
                                            ! resource file name
      IF(TRIM(name) == "full" ) THEN
       name = " "              ! blank instance name for full (1)
      ELSE
       name = TRIM(name)       ! instance name for others
      END IF
      call CO_GridCompSetServices1_(gc,chemReg,name,rc=status)
      VERIFY_(STATUS)
   end do

   RETURN_(ESMF_SUCCESS)

   end subroutine CO_GridCompSetServices

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  CO_GridCompInitialize --- Initialize CO_GridComp
!
! !INTERFACE:
!

   subroutine CO_GridCompInitialize ( gcCO, w_c, impChem, expChem, &
                                      nymd, nhms, cdt, rc )

! !USES:

  IMPLICIT NONE

! !INPUT PARAMETERS:

   TYPE(Chem_Bundle), intent(in) :: w_c          ! Chemical tracer fields      
   INTEGER, INTENT(IN) :: nymd, nhms             ! time
   REAL,    INTENT(IN) :: cdt                    ! chemical timestep (secs)


! !OUTPUT PARAMETERS:

   TYPE(CO_GridComp), INTENT(INOUT) :: gcCO      ! Grid Component
   TYPE(ESMF_State), INTENT(INOUT)  :: impChem   ! Import State
   TYPE(ESMF_State), INTENT(INOUT)  :: expChem   ! Export State
   INTEGER, INTENT(OUT) ::  rc                   ! Error return code:
                                                 !  0 - all is well
                                                 !  1 - 

! !DESCRIPTION: Initializes the CO Grid Component. Multiple instance
!               version.
!
! !REVISION HISTORY:
!
!  27Feb2008  da Silva  Introduced multiple instances
!
!EOP
!-------------------------------------------------------------------------

   CHARACTER(LEN=*), PARAMETER :: myname = 'CO_GridCompInitialize'
   CHARACTER(LEN=255) :: rcbasen = 'CO_GridComp'
   CHARACTER(LEN=255) :: name
   
   integer i, ier, n
   REAL :: c1,c2,c3,c4

!  Load resource file
!  ------------------
   CALL I90_loadf ( TRIM(rcbasen)//'.rc', ier )
   if ( ier .NE. 0 ) then
      rc = 10
      return
   end if

!  Parse resource file
!  -------------------
   CALL I90_label ( 'CO_instances:', ier )
   if ( ier .NE. 0 ) then
      rc = 20
      return
   end if

!  First determine how many instances we have
!  ------------------------------------------   
   n = 0
   do while ( ier .EQ. 0 )
      CALL I90_gtoken( name, ier )
      n = n + 1
   end do
   if ( n .EQ. 0 ) then
      rc = 30
      return
   end if
   
!  We cannot have fewer instances than the number of
!   CO bins in the registry (it is OK to have less, though)
!  --------------------------------------------------------
   if ( n .LT. w_c%reg%n_CO ) then
        rc = 35
        return
   else if ( n .GT. w_c%reg%n_CO ) then
        if (MAPL_AM_I_ROOT()) &
        print *, trim(myname)// &
                 ': fewer CO bins than possible CO instances: ',&
                 n, w_c%reg%n_CO
   end if
   n = min(n,w_c%reg%n_CO )
   gcCO%n = n

!  Next allocate necessary memory
!  ------------------------------
   allocate ( gcCO%gcs(n), stat=ier )    
   if ( ier .NE. 0 ) then
      rc = 40
      return
   end if

!  Record name of each instance
!  ----------------------------
   CALL I90_label ( 'CO_instances:', ier )
   do i = 1, n
      CALL I90_gtoken( name, ier )
      if ( ier .NE. 0 ) then
         rc = 40
         return
      end if
                                            ! resource file name
      gcCO%gcs(i)%rcfilen = trim(rcbasen)//'---'//trim(name)//'.rc'
      gcCO%gcs(i)%instance = i              ! instance number 
      IF(TRIM(name) == "full" ) THEN
       gcCO%gcs(i)%iname = " "              ! blank instance name for full (1)
      ELSE
       gcCO%gcs(i)%iname = TRIM(name)       ! instance name for others
      END IF
   end do    

!  Next initialize each instance
!  -----------------------------
   do i = 1, gcCO%n
      IF(MAPL_AM_I_ROOT()) THEN
       PRINT *," "
       PRINT *,myname,": Initializing instance ",TRIM(gcCO%gcs(i)%iname)," [",gcCO%gcs(i)%instance,"]"
      END IF
      call CO_SingleInstance_ ( CO_GridCompInitialize1_, i, &
                                gcCO%gcs(i), w_c, impChem, expChem,  &
                                nymd, nhms, cdt, ier )
      if ( ier .NE. 0 ) then
         rc = 1000+ier
         return
      end if
   end do

!  Get Henrys Law cts for the parameterized convective wet removal
!  -----------------------------------------------------------
   CALL get_HenrysLawCts('CO',c1,c2,c3,c4)  
   w_c%reg%Hcts(1,w_c%reg%i_CO : w_c%reg%j_CO)=c1
   w_c%reg%Hcts(2,w_c%reg%i_CO : w_c%reg%j_CO)=c2
   w_c%reg%Hcts(3,w_c%reg%i_CO : w_c%reg%j_CO)=c3
   w_c%reg%Hcts(4,w_c%reg%i_CO : w_c%reg%j_CO)=c4


!  All done
!  --------
   CALL I90_FullRelease( ier )
   IF( ier /= 0 ) THEN
    PRINT *,myname,": I90_FullRelease not successful."
    rc = 40
   END IF


 end subroutine CO_GridCompInitialize

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  CO_GridCompRun --- Run CO_GridComp
!
! !INTERFACE:
!

   subroutine CO_GridCompRun ( gcCO, w_c, impChem, expChem, &
                                      nymd, nhms, cdt, rc )

! !USES:

  IMPLICIT NONE

! !INPUT PARAMETERS:

   TYPE(Chem_Bundle), intent(in) :: w_c          ! Chemical tracer fields      
   INTEGER, INTENT(IN) :: nymd, nhms             ! time
   REAL,    INTENT(IN) :: cdt                    ! chemical timestep (secs)


! !OUTPUT PARAMETERS:

   TYPE(CO_GridComp), INTENT(INOUT) :: gcCO      ! Grid Component
   TYPE(ESMF_State), INTENT(INOUT)  :: impChem   ! Import State
   TYPE(ESMF_State), INTENT(INOUT)  :: expChem   ! Export State
   INTEGER, INTENT(OUT) ::  rc                   ! Error return code:
                                                 !  0 - all is well
                                                 !  1 - 

! !DESCRIPTION: Runs the CO Grid Component. Multiple instance
!               version.
!
! !REVISION HISTORY:
!
!  27Feb2008  da Silva  Introduced multiple instances
!
!EOP
!-------------------------------------------------------------------------

   integer i, ier

   do i = 1, gcCO%n
      call CO_SingleInstance_ ( CO_GridCompRun1_, i, &
                                gcCO%gcs(i), w_c, impChem, expChem, &
                                nymd, nhms, cdt, ier )
      if ( ier .NE. 0 ) then
         rc = i * 1000+ier
         return
      end if
   end do

 end subroutine CO_GridCompRun


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  CO_GridCompFinalize --- Initialize CO_GridComp
!
! !INTERFACE:
!

   subroutine CO_GridCompFinalize ( gcCO, w_c, impChem, expChem, &
                                      nymd, nhms, cdt, rc )

! !USES:

  IMPLICIT NONE

! !INPUT PARAMETERS:

   TYPE(Chem_Bundle), intent(in) :: w_c          ! Chemical tracer fields      
   INTEGER, INTENT(IN) :: nymd, nhms             ! time
   REAL,    INTENT(IN) :: cdt                    ! chemical timestep (secs)


! !OUTPUT PARAMETERS:

   TYPE(CO_GridComp), INTENT(INOUT) :: gcCO      ! Grid Component
   TYPE(ESMF_State), INTENT(INOUT)  :: impChem   ! Import State
   TYPE(ESMF_State), INTENT(INOUT)  :: expChem   ! Export State
   INTEGER, INTENT(OUT) ::  rc                   ! Error return code:
                                                 !  0 - all is well
                                                 !  1 - 

! !DESCRIPTION: Finalizes the CO Grid Component. Multiple instance
!               version.
!
! !REVISION HISTORY:
!
!  27Feb2008  da Silva  Introduced multiple instances
!
!EOP
!-------------------------------------------------------------------------

   integer i, ier

   do i = 1, gcCO%n
      call CO_SingleInstance_ ( CO_GridCompFinalize1_, i, &
                                gcCO%gcs(i), w_c, impChem, expChem, &
                                nymd, nhms, cdt, ier )
      if ( ier .NE. 0 ) then
         rc = i * 1000+ier
         return
      end if
   end do

   deallocate ( gcCO%gcs, stat=ier )    
   gcCO%n = -1

 end subroutine CO_GridCompFinalize

 subroutine CO_GridCompSetServices1_(  gc, chemReg, iname, rc)
 type(ESMF_GridComp), intent(INOUT) :: GC
 type(Chem_Registry), intent(INOUT) :: chemReg
 character(len=*),    intent(IN   ) :: iname
 integer,             intent(OUT  ) :: rc

 integer :: Status
 character(len=ESMF_MAXSTR) :: Iam

 Iam ="CO_GridCOmpSetServices1_"

  call MAPL_AddImportSpec(GC,             &
       SHORT_NAME = 'CO_OH'//iname,       &
       LONG_NAME  = 'source species'  ,   &
       UNITS      = '1',                  &
       DIMS       = MAPL_DimsHorzVert,    &
       VLOCATION  = MAPL_VLocationCenter, &
       RESTART    = MAPL_RestartSkip,     &
       RC         = STATUS)
  VERIFY_(STATUS)
  call MAPL_AddImportSpec(GC,             &
       SHORT_NAME = 'CO_CH4'//iname,      &
       LONG_NAME  = 'source species'  ,   &
       UNITS      = '1',                  &
       DIMS       = MAPL_DimsHorzVert,    &
       VLOCATION  = MAPL_VLocationCenter, &
       RESTART    = MAPL_RestartSkip,     &
       RC         = STATUS)
  VERIFY_(STATUS)
  call MAPL_AddImportSpec(GC,           &
       SHORT_NAME = 'CO_BF'//iname,     &
       LONG_NAME  = 'source species'  , &
       UNITS      = '1',                &
       DIMS       = MAPL_DimsHorzOnly,  &
       VLOCATION  = MAPL_VLocationNone, &
       RESTART    = MAPL_RestartSkip,     &
       RC         = STATUS)
  VERIFY_(STATUS)
  call MAPL_AddImportSpec(GC,           &
       SHORT_NAME = 'CO_FS'//iname,     &
       LONG_NAME  = 'source species'  , &
       UNITS      = '1',                &
       DIMS       = MAPL_DimsHorzOnly,  &
       VLOCATION  = MAPL_VLocationNone, &
       RESTART    = MAPL_RestartSkip,     &
       RC         = STATUS)
  VERIFY_(STATUS)
  call MAPL_AddImportSpec(GC,           &
       SHORT_NAME = 'CO_ISOP'//iname,   &
       LONG_NAME  = 'source species'  , &
       UNITS      = '1',                &
       DIMS       = MAPL_DimsHorzOnly,  &
       VLOCATION  = MAPL_VLocationNone, &
       RESTART    = MAPL_RestartSkip,     &
       RC         = STATUS)
  VERIFY_(STATUS)
  call MAPL_AddImportSpec(GC,           &
       SHORT_NAME = 'CO_NVOC'//iname,   &
       LONG_NAME  = 'source species'  , &
       UNITS      = '1',                &
       DIMS       = MAPL_DimsHorzOnly,  &
       VLOCATION  = MAPL_VLocationNone, &
       RESTART    = MAPL_RestartSkip,     &
       RC         = STATUS)
  VERIFY_(STATUS)
  call MAPL_AddImportSpec(GC, &
       SHORT_NAME = 'CO_TERP'//iname,   &
       LONG_NAME  = 'source species'  , &
       UNITS      = '1',                &
       DIMS       = MAPL_DimsHorzOnly,  &
       VLOCATION  = MAPL_VLocationNone, &
       RESTART    = MAPL_RestartSkip,     &
       RC         = STATUS)
  VERIFY_(STATUS)
  call MAPL_AddImportSpec(GC,           &
       SHORT_NAME = 'CO_BIOMASS'//iname,&
       LONG_NAME  = 'source species'  , &
       UNITS      = '1',                &
       DIMS       = MAPL_DimsHorzOnly,  &
       VLOCATION  = MAPL_VLocationNone, &
       RESTART    = MAPL_RestartSkip,     &
       RC         = STATUS)
  VERIFY_(STATUS)

  RETURN_(ESMF_SUCCESS) 

 end subroutine CO_GridCompSetServices1_

!--------------------------------------------------------------------------

!                      Single Instance Methods

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  CO_GridCompInitialize --- Initialize CO_GridComp
!
! !INTERFACE:
!

   subroutine CO_GridCompInitialize1_ ( gcCO, w_c, impChem, expChem, &
                                        nymd, nhms, cdt, rc )

! !USES:

  IMPLICIT NONE

! !INPUT PARAMETERS:

   TYPE(Chem_Bundle), intent(in) :: w_c          ! Chemical tracer fields      
   INTEGER, INTENT(IN) :: nymd, nhms             ! time
   REAL,    INTENT(IN) :: cdt                    ! chemical timestep (secs)


! !OUTPUT PARAMETERS:

   TYPE(CO_GridComp1), INTENT(INOUT) :: gcCO     ! Grid Component
   TYPE(ESMF_State), INTENT(INOUT)  :: impChem   ! Import State
   TYPE(ESMF_State), INTENT(INOUT)  :: expChem   ! Export State
   INTEGER, INTENT(OUT) ::  rc                   ! Error return code:
                                                 !  0 - all is well
                                                 !  1 - 

! !DESCRIPTION: Initializes the CO Grid Component. It primarily sets
!               the import state for each active constituent package.
!
! !REVISION HISTORY:
!
!  18Sep2003 da Silva  First crack.
!  31May2005  Nielsen  Mods for 7 CO bins, 5 region masks
!  04Nov2005     Bian  CO tagged to 4 regions 
!                      (global, North America, South America, and Africa)
!                      for CR-AVE
!
!EOP
!-------------------------------------------------------------------------

   CHARACTER(LEN=*), PARAMETER :: myname = 'CO_GridCompInitialize1'

   CHARACTER(LEN=255) :: rcfilen 

   INTEGER :: ios, j, n
   INTEGER, ALLOCATABLE :: ier(:)
   INTEGER :: i1, i2, im, j1, j2, jm, km
   INTEGER :: nTimes, begTime, incSecs
   INTEGER :: nbeg, nend, nymd1, nhms1

   REAL :: limitN, limitS
   REAL, ALLOCATABLE :: var2d(:,:)

   rcfilen = gcCO%rcfilen
   gcCO%name = 'GEOS-5/GOCART Parameterized CO Package'

!  Initialize local variables
!  --------------------------
   rc = 0
   i1 = w_c%grid%i1
   i2 = w_c%grid%i2
   im = w_c%grid%im

   j1 = w_c%grid%j1
   j2 = w_c%grid%j2
   jm = w_c%grid%jm

   km = w_c%grid%km

   nbeg  = w_c%reg%i_CO
   nend  = w_c%reg%j_CO

!  It requires 1 bin
!  -----------------
   if ( nbeg /= nend ) then
      IF(MAPL_AM_I_ROOT()) PRINT *,myname,": Must have only 1 bin at the single instance level"
      rc = 1
      return 
   end if

!  Allocate memory, etc
!  --------------------
   CALL INIT_()
   IF ( rc /= 0 ) RETURN

!  Load resource file
!  ------------------
   CALL I90_loadf ( TRIM(rcfilen), ier(1) )
   IF ( ier(1) .NE. 0 ) THEN
      CALL final_(10)
      RETURN
   END IF
   ier(:)=0


!  Run-time debug switch
!  ---------------------
   call i90_label ( 'DEBUG:', ier(1) )
   n = i90_gint ( ier(2) )
   if(n /= 0) then
    gcco%dbg = .true.
   else
    gcco%dbg = .false.
   end if

!  Possibly OH is in a different unit
!  Allowable choices are: "mol/mol" or "mol mol-1"
!  or else behavior is as though input in 
!  "molecules cm-3"
!  Should this checking be done now in ExtData?
!  --------------------------------------------
   CALL i90_label ('units_oh:',ier(3) )
   if (ier(3) /= 0) then
      gcCO%units_oh = " "
      ier(3) = 0
   else
      CALL I90_gtoken(gcCO%units_oh,ier(4) )
   end if

   IF( ANY( ier(:) /= 0 ) ) THEN
    CALL final_(21)
    RETURN
   END IF

!  Set the initial CO surface fluxes to zero
!  -----------------------------------------
   gcCO%COsfcFlux(i1:i2,j1:j2) = 0.00

   DEALLOCATE(ier)

   RETURN

CONTAINS

   SUBROUTINE init_()

   INTEGER ios, nerr
   nerr = 128
   ALLOCATE ( gcCO%eCO_bioburn(i1:i2,j1:j2), & 
              gcCO%eCO_bioburn_(i1:i2,j1:j2), & 
              gcCO%eCO_biofuel(i1:i2,j1:j2), &
              gcCO%eCO_fosfuel(i1:i2,j1:j2), &
              gcCO%COsfcFlux(i1:i2,j1:j2), &
              gcCO%eCO_iso(i1:i2,j1:j2), &
              gcCO%eCO_mon(i1:i2,j1:j2), &
              gcCO%eCO_mtn(i1:i2,j1:j2), &
              gcCO%CH4(i1:i2,j1:j2,km), &
              gcCO%OHnd(i1:i2,j1:j2,km), &
              ier(nerr),STAT=ios )
   ier = 0
   IF ( ios /= 0 ) rc = 100
   END SUBROUTINE init_

   SUBROUTINE final_(ierr)
   INTEGER :: ierr
   INTEGER ios
   DEALLOCATE ( gcCO%eCO_bioburn, gcCO%eCO_biofuel, gcCO%eCO_fosfuel, & 
                gcCO%eCO_bioburn_, &
                gcCO%COsfcFlux, gcCO%eCO_iso, gcCO%eCO_mon, &
                gcCO%eCO_mtn, gcCO%CH4, gcCO%OHnd, &
                ier, STAT=ios )
   CALL I90_release()
   rc = ierr
   END SUBROUTINE final_

 END SUBROUTINE CO_GridCompInitialize1_

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  CO_GridCompRun --- The Chem Driver 
!
! !INTERFACE:
!

   SUBROUTINE CO_GridCompRun1_ ( gcCO, w_c, impChem, expChem, &
                                 nymd, nhms, cdt, rc )

#define MR_PBL
!#define SFLUX_PBL

! !USES:

  IMPLICIT NONE

! !INPUT/OUTPUT PARAMETERS:

   TYPE(CO_GridComp1), INTENT(INOUT) :: gcCO     ! Grid Component
   TYPE(Chem_Bundle), INTENT(INOUT) :: w_c       ! Chemical tracer fields   

! !INPUT PARAMETERS:

   TYPE(ESMF_State), INTENT(inout) :: impChem    ! Import State
   INTEGER, INTENT(IN) :: nymd, nhms             ! time
   REAL,    INTENT(IN) :: cdt                    ! chemical timestep (secs)


! !OUTPUT PARAMETERS:

   TYPE(ESMF_State), intent(inout) :: expChem    ! Export State
   INTEGER, INTENT(OUT) ::  rc                   ! Error return code:
                                                 !  0 - all is well
                                                 !  1 -
 
! !DESCRIPTION: This routine implements the CO Driver for INTEX. That 
!               is, adds chemical tendencies to each of the constituents,
!  Note: water wapor, the first constituent is not considered a chemical
!  constituents.
!
! !REVISION HISTORY:
!
!  18Sep2003 da Silva  First crack.
!  31May2005  Nielsen  Mods for 7 tags, 5 regions
!  04Nov2005     Bian  CO tagged to 4 regions    
!  13Apr2005     Bian  CO tagged to emissions    
!
!EOP
!-------------------------------------------------------------------------

   CHARACTER(LEN=*), PARAMETER :: myname = 'CO_GridCompRun'
   CHARACTER(LEN=*), PARAMETER :: Iam = myname

!  Input fields from fvGCM
!  -----------------------
   REAL, POINTER, DIMENSION(:,:)   ::  pblh  => null()
   REAL, POINTER, DIMENSION(:,:,:) ::  T     => null()
   REAL, POINTER, DIMENSION(:,:,:) ::  rhoa  => null()
   REAL, POINTER, DIMENSION(:,:,:) ::  zle   => null()

   INTEGER :: i1, i2, im, j1, j2, jm, km, ios, idiag, iXj
   INTEGER :: i, j, k, kReverse, n, nbeg, nend
   INTEGER :: nymd1, nhms1, ier(8)
   integer :: iregWant

   REAL, PARAMETER :: nsuba=6.022E+26
   REAL, PARAMETER :: mwtAir=28.97
   REAL, PARAMETER :: mwtCO=28.01
   REAL, PARAMETER :: rstar=8.3143E+03
   REAL, PARAMETER :: rpstd=1.00E-05

   REAL    :: qmin, qmax, toMass, c2co
   REAL    :: fiso, fmtn, fmon

   REAL, ALLOCATABLE :: CH4nd(:,:,:)
   REAL, ALLOCATABLE :: OHnd(:,:,:)
   REAL, ALLOCATABLE :: pe(:,:,:),p(:,:,:),nd(:,:,:)
   REAL, ALLOCATABLE :: rkoh(:,:,:),rkch4(:,:,:)

   real, pointer, dimension(:,:,:) :: ptr3d => null()
   real, pointer, dimension(:,:)   :: ptr2d => null()

#define EXPORT   expChem
#define iNAME    TRIM(gcCO%iname)

#define COEM     CO_emis
#define COCL     CO_column
#define COSC     CO_surface
#define COPD     CO_prod
#define COLS     CO_loss

   integer :: STATUS

#include "CO_GetPointer___.h"

!  Initialize local variables
!  --------------------------
   rc = 0
   i1 = w_c%grid%i1
   i2 = w_c%grid%i2
   im = w_c%grid%im

   j1 = w_c%grid%j1
   j2 = w_c%grid%j2
   jm = w_c%grid%jm
   km = w_c%grid%km

   iXj = ( i2 - i1 + 1 ) * ( j2 - j1 + 1 )

   nbeg  = w_c%reg%i_CO
   nend  = w_c%reg%j_CO

!  It requires 1 bin
!  -----------------
   if ( nbeg /= nend ) then
      IF(MAPL_AM_I_ROOT()) PRINT *,myname,": Must have only 1 bin at the single instance level"
      rc = 1
      return 
   end if

!  Conversion factor, molecules CO cm^-2 s^-1 to kg CO m^-2 s^-1
!  -------------------------------------------------------------
   toMass = 1.00E+04*mwtCO/nsuba
   c2co   = mwtCO/12.00

!  Update emissions and OH number density once each day.
!  The latter appears to be in molecules cm^-3.
!  -----------------------------------------------------
   
!   Selections based on biomass burning emission set chosen
!   Currently, parse on:
!    harvard   -> molecules cm-2 s-1, need to convert to mass
!    modisfire -> kg CO m-2 s-1
!    else      -> based on dry matter consumed

!   Biomass Burning -- select on known inventories
!   ----------------------------------------------

!   Harvard biomass burning climatology, is in molecules cm^-2 s^-1
!   ---------------------------------------------------------------
    call MAPL_GetPointer(impChem, ptr2d, 'CO_BIOMASS'//iNAME,rc=status)
    VERIFY_(STATUS)
    gcCO%eCO_bioburn = ptr2d

! Background OH, for loss term
! ----------------------------
    call MAPL_GetPointer(impChem, ptr3d, 'CO_OH'//iNAME,rc=status)
    VERIFY_(STATUS)
    gcCO%OHnd = ptr3d

! Background CH4, for source term.
! NOTE: Return zeroes in all but the global instantiation.
! --------------------------------------------------------
    call MAPL_GetPointer(impChem, ptr3d, 'CO_CH4'//iNAME,rc=status)
    VERIFY_(STATUS)
    gcCO%CH4 = ptr3d

! Biofuel source
! --------------
    call MAPL_GetPointer(impChem, ptr2d, 'CO_BF'//iNAME,rc=status)
    VERIFY_(STATUS)
    gcCO%eCO_biofuel = ptr2d

! Fossil fuel source
! ------------------
    call MAPL_GetPointer(impChem, ptr2d, 'CO_FS'//iNAME,rc=status)
    VERIFY_(STATUS)
    gcCO%eCO_fosfuel = ptr2d

! Isoprene source
! ---------------
    call MAPL_GetPointer(impChem, ptr2d, 'CO_ISOP'//iNAME,rc=status)
    VERIFY_(STATUS)
    gcCO%eCO_iso = ptr2d

! VOC source
! ----------
    call MAPL_GetPointer(impChem, ptr2d, 'CO_NVOC'//iNAME,rc=status)
    VERIFY_(STATUS)
    gcCO%eCO_mon = ptr2d

! Monoterpene source
! ------------------
    call MAPL_GetPointer(impChem, ptr2d, 'CO_TERP'//iNAME,rc=status)
    VERIFY_(STATUS)
    gcCO%eCO_mtn = ptr2d

   IF(gcCO%DBG) THEN
    CALL pmaxmin('CO: eCO_bioburn', gcCO%eCO_bioburn, qmin, qmax, iXj,1, 1. )
    CALL pmaxmin('CO: eCO_biofuel', gcCO%eCO_biofuel, qmin, qmax, iXj,1, 1. )
    CALL pmaxmin('CO: eCO_fosfuel', gcCO%eCO_fosfuel, qmin, qmax, iXj,1, 1. )
    CALL pmaxmin('CO: eCO_iso',     gcCO%eCO_iso,     qmin, qmax, iXj,1, 1. )
    CALL pmaxmin('CO: eCO_mon',     gcCO%eCO_mon,     qmin, qmax, iXj,1, 1. )
    CALL pmaxmin('CO: eCO_mtn',     gcCO%eCO_mtn,     qmin, qmax, iXj,1, 1. )
   END IF

!   Save this in case we need to apply diurnal cycle
!   ------------------------------------------------
   if ( w_c%diurnal_bb ) then
        gcCO%eCO_bioburn_(:,:) = gcCO%eCO_bioburn(:,:)
   end if

!  Apply diurnal cycle if so desired
!  ---------------------------------
   if ( w_c%diurnal_bb ) then
      call Chem_BiomassDiurnal ( gcCO%eCO_bioburn, gcCO%eCO_bioburn_,   &
                                 w_c%grid%lon(:,:)*radToDeg, &
                                 w_c%grid%lat(:,:)*radToDeg, nhms, cdt )      
   end if


!  Allocate temporary workspace
!  ----------------------------
   allocate ( pe(i1:i2,j1:j2,km+1), p(i1:i2,j1:j2,km), nd(i1:i2,j1:j2,km), &
              rkoh(i1:i2,j1:j2,km), rkch4(i1:i2,j1:j2,km), &
              CH4nd(i1:i2,j1:j2,km), OHnd(i1:i2,j1:j2,km), stat = ios )

   if ( ios /= 0 ) then
      rc = 3
      return
   end if

!  Layer interface pressures
!  -------------------------
   pe(i1:i2,j1:j2,1)=w_c%grid%ptop
   DO k=2,km+1
    pe(i1:i2,j1:j2,k)=pe(i1:i2,j1:j2,k-1)+w_c%delp(i1:i2,j1:j2,k-1)
   END DO

!  Layer mean pressures
!  --------------------
   DO k=1,km
    p(i1:i2,j1:j2,k)=(pe(i1:i2,j1:j2,k)+pe(i1:i2,j1:j2,k+1))*0.50
   END DO
 
!  Get imports
!  -----------
   call MAPL_GetPointer( impChem, pblh,  'ZPBL',    rc=ier(1) ) 
   call MAPL_GetPointer( impChem, T,     'T',       rc=ier(2) ) 
   call MAPL_GetPointer( impChem, rhoa,  'AIRDENS', rc=ier(3) ) 
   call MAPL_GetPointer( impChem, zle,   'ZLE',     rc=ier(4) ) 

   if ( any(ier(1:4) /= 0) ) then
        rc = 10
        return
   end if

   IF(gcCO%DBG) THEN
    CALL pmaxmin('CO: pblh', pblh, qmin, qmax, iXj,    1, 1. )
    CALL pmaxmin('CO:    T',    T, qmin, qmax, iXj,   km, 1. )
    CALL pmaxmin('CO: rhoa', rhoa, qmin, qmax, iXj,   km, 1. )
    CALL pmaxmin('CO:  zle',  zle, qmin, qmax, iXj, km+1, 1. )
   END IF

!  Number density
!  --------------
   nd(i1:i2,j1:j2,1:km)= nsuba*p(i1:i2,j1:j2,1:km)/ &
                        (rstar*t(i1:i2,j1:j2,1:km))

!  CH4 number density.  CH4 on file is in mole fraction.
!  -----------------------------------------------------
   CH4nd(i1:i2,j1:j2,1:km)=gcCO%CH4(i1:i2,j1:j2,1:km)* &
                                 nd(i1:i2,j1:j2,1:km)

!  OH number density. Handle mole fraction or number density.
!  ----------------------------------------------------------
   if ( (trim(gcCO%units_oh) .eq. 'mol/mol') .or. (trim(gcCO%units_oh) .eq. 'mol mol-1') ) then
       OHnd(i1:i2,j1:j2,1:km) = gcCO%OHnd(i1:i2,j1:j2,1:km)* &
                                       nd(i1:i2,j1:j2,1:km)
   else
       ! assume that units are 'molecules cm-3' and convert to 'molecules m^-3'
       OHnd(i1:i2,j1:j2,1:km) = gcCO%OHnd(i1:i2,j1:j2,1:km)*1.00E+06
   end if    

!  Clear surface flux array
!  ------------------------
   gcCO%COsfcFlux(i1:i2,j1:j2) = 0.0

#if defined( MR_PBL)
!  Emissions, direct update of mixing ratio 
!  ----------------------------------------
   call CO_Emission(rc)
#endif

!  Convert carbon monoxide from mole fraction to number density
!  ------------------------------------------------------------
   w_c%qa(nbeg)%data3d(i1:i2,j1:j2,1:km) = &
                  w_c%qa(nbeg)%data3d(i1:i2,j1:j2,1:km)*nd(i1:i2,j1:j2,1:km)

!  Loss due to OH.
!  ---------------
   rkoh(i1:i2,j1:j2,1:km) = &
                  1.5E-13*1.00E-06*(1.00+0.6*p(i1:i2,j1:j2,1:km)*rpstd)
   n = gcCO%instance 
   IF(ASSOCIATED(CO_loss)) CO_loss(i1:i2,j1:j2) = 0.
   DO k = 1, km

    IF(ASSOCIATED(CO_loss)) CO_loss(i1:i2,j1:j2) = CO_loss(i1:i2,j1:j2) &
       + w_c%qa(nbeg)%data3d(i1:i2,j1:j2,k)*rkoh(i1:i2,j1:j2,k) &
       * OHnd(i1:i2,j1:j2,k)/nd(i1:i2,j1:j2,k) &
       * mwtCO/mwtAir*w_c%delp(i1:i2,j1:j2,k)/grav

    w_c%qa(nbeg)%data3d(i1:i2,j1:j2,k) = &
         w_c%qa(nbeg)%data3d(i1:i2,j1:j2,k)*(1.00-cdt* &
         rkoh(i1:i2,j1:j2,k)*OHnd(i1:i2,j1:j2,k))

   END DO ! Next layer, k

!  CH4 production
!  --------------
   rkch4(i1:i2,j1:j2,1:km)= &
                 2.45e-12*1.00E-06*exp(-1775./t(i1:i2,j1:j2,1:km))
   n = gcCO%instance 
   IF(ASSOCIATED(CO_prod)) CO_prod(i1:i2,j1:j2) = 0.
   DO k = 1, km

    IF(ASSOCIATED(CO_prod)) CO_prod(i1:i2,j1:j2) = CO_prod(i1:i2,j1:j2) &
       + rkch4(i1:i2,j1:j2,k)*OHnd(i1:i2,j1:j2,k)*CH4nd(i1:i2,j1:j2,k) &
       / nd(i1:i2,j1:j2,k) * mwtCO/mwtAir*w_c%delp(i1:i2,j1:j2,k)/grav

    w_c%qa(nbeg)%data3d(i1:i2,j1:j2,k)=w_c%qa(nbeg)%data3d(i1:i2,j1:j2,k)+cdt* &
         rkch4(i1:i2,j1:j2,k)*OHnd(i1:i2,j1:j2,k)* &
         CH4nd(i1:i2,j1:j2,k)

   END DO ! Next layer, k

!  Return to mole fraction
!  -----------------------
   w_c%qa(nbeg)%data3d(i1:i2,j1:j2,1:km)=w_c%qa(nbeg)%data3d(i1:i2,j1:j2,1:km)/nd(i1:i2,j1:j2,1:km)

!  Surface concentration in ppbv
!  -----------------------------
   n = gcCO%instance 
    if(associated(CO_surface)) &
      CO_surface(i1:i2,j1:j2) = w_c%qa(nbeg)%data3d(i1:i2,j1:j2,km)*1.e9

!  Column burden in kg m-2
!  -----------------------
   n = gcCO%instance 
    if(associated(CO_column)) then
     CO_column(i1:i2,j1:j2) = 0.
     do k = 1, km
      CO_column(i1:i2,j1:j2) &
       =   CO_column(i1:i2,j1:j2) &
         +   w_c%qa(nbeg)%data3d(i1:i2,j1:j2,k)*mwtCO/mwtAir &
           * w_c%delp(i1:i2,j1:j2,k)/grav
     enddo
    endif

!  CO Surface Emission Flux in kg m-2 s-1
!  --------------------------------------
    n = gcCO%instance 
    if(associated(CO_emis)) &
         CO_emis(i1:i2,j1:j2) = gcCO%COsfcFlux(i1:i2,j1:j2)

   IF(gcCO%DBG) THEN
     n = gcCO%instance 
     if(associated(CO_emis)) &
     CALL pmaxmin('CO: emis', CO_emis(i1:i2,j1:j2), qmin, qmax, &
                   iXj,1, 1. )
     if(associated(CO_loss)) &
     CALL pmaxmin('CO: loss', CO_loss(i1:i2,j1:j2), qmin, qmax, &
                   iXj,1, 1. )
     if(associated(CO_prod)) &
     CALL pmaxmin('CO: prod', CO_prod(i1:i2,j1:j2), qmin, qmax, &
                   iXj,1, 1. )
     if(associated(CO_column)) &
     CALL pmaxmin('CO: column', CO_column(i1:i2,j1:j2), qmin, qmax, &
                   iXj,1, 1. )
     if(associated(CO_surface)) &
     CALL pmaxmin('CO: surface', CO_surface(i1:i2,j1:j2), qmin, qmax,&
                   iXj,1, 1. )
   END IF

!  Housekeeping
!  ------------
   DEALLOCATE(nd,p,pe,rkoh,rkch4,CH4nd,OHnd,STAT=ier(1))

   RETURN

CONTAINS
!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
! !DESCRIPTION: Updates the CO concentration with emissions every timestep
!
! !IROUTINE:  CO_Emission - Adds emissions for CO for one timestep
!             We have emissions from 4 sources, which are distributed
!             differently in the vertical
!             1) fossil fuel - emitted at surface
!             2) biofuel sources - emitted at surface 
!             3) biomass burning - uniformly mixed in PBL
!             4) biogenic - emitted at surface
!                           include: isoprene, converting factor 0.15
!                                    terpene,  converting factor 0.2
!                                    nvoc,     converting factor 0.2
! !REVISION HISTORY:
!
!  17Oct2005, Bian!
!  14Apr2006, Bian: Add indirect NMHC from FF (0.20), BF (0.19), BB (0.11)
!                   Add seasonality for FF
!                   Modify FF & BF over Asia region (1.39) for Streets' data
!  18Mar2011, Nielsen: Simplified PBL partitioning for biomass burning emissions   
!
! !INTERFACE:
!
!EOP
!-------------------------------------------------------------------------
   SUBROUTINE CO_Emission ( rc )
!-------------------------------------------------------------------------

! !USES:

  IMPLICIT NONE

! !INPUT PARAMETERS:

! !OUTPUT PARAMETERS:

   INTEGER, INTENT(OUT) :: rc  ! Error return code

! !LOCAL VARIABLES

   CHARACTER(LEN=*), PARAMETER :: myname = 'CO_Emission'

   INTEGER :: i, j, k, kt, minkPBL
   INTEGER, ALLOCATABLE :: index(:)

   REAL, PARAMETER :: mwtAir=28.97
   REAL, PARAMETER :: mwtCO=28.01
   REAL, ALLOCATABLE :: pblLayer(:,:),sfcFlux(:,:),fPBL(:,:,:)

   rc    = 0

! Grab some memory for manipulating surface fluxes
! ------------------------------------------------
   ALLOCATE(sfcFlux(i1:i2,j1:j2),STAT=ios)

! Biomass burning
! ---------------
   BioBurn: IF(gcCO%doingBB) THEN

    sfcFlux(i1:i2,j1:j2)=gcCO%eCO_bioburn(i1:i2,j1:j2)
    gcCO%COsfcFlux(i1:i2,j1:j2)=gcCO%COsfcFlux(i1:i2,j1:j2)+sfcFlux(i1:i2,j1:j2)

! Find the layer that contains the PBL.
! Layer thicknesses are ZLE(:,:,0:km).
! -------------------------------------
    ALLOCATE(index(0:km),STAT=ios)
    ALLOCATE(pblLayer(i1:i2,j1:j2),STAT=ios)
    DO j=j1,j2
     DO i=i1,i2
      index(0:km)=0
      WHERE(zle(i,j,0:km)-zle(i,j,km) > pblh(i,j)) index(0:km)=1
      pblLayer(i,j)=SUM(index)
     END DO
    END DO
    DEALLOCATE(index,STAT=ios)
    minkPBL=MINVAL(pblLayer)

! Determine partitioning fraction based on layer thicknesses
! ----------------------------------------------------------
    ALLOCATE(fPBL(i1:i2,j1:j2,1:km),STAT=ios)
    fPBL(i1:i2,j1:j2,1:km)=0.00
    DO j=j1,j2
     DO i=i1,i2
      kt=pblLayer(i,j)
      DO k=kt,km
       fPBL(i,j,k)=(zle(i,j,k-1)-zle(i,j,k))/(zle(i,j,kt-1)-zle(i,j,km))
      END DO
     END DO
    END DO

! Partition surface flux into layers within the PBL
! -------------------------------------------------
    DO k=minkPBL,km
     w_c%qa(nbeg)%data3d(i1:i2,j1:j2,k) = w_c%qa(nbeg)%data3d(i1:i2,j1:j2,k)+ &
                                          sfcFlux(i1:i2,j1:j2)*fPBL(i1:i2,j1:j2,k)*cdt* &
                                          (mwtAir/mwtCO)/(w_c%delp(i1:i2,j1:j2,k)/grav)
    END DO

! Release memory
! --------------
    DEALLOCATE(fPBL,STAT=ios)
    DEALLOCATE(pblLayer,STAT=ios)

   END IF BioBurn

! Biogenic
! --------
   sfcFlux(i1:i2,j1:j2) = gcCO%eCO_iso(i1:i2,j1:j2)+gcCO%eCO_mon(i1:i2,j1:j2)+gcCO%eCO_mtn(i1:i2,j1:j2)
   gcCO%COsfcFlux(i1:i2,j1:j2) = gcCO%COsfcFlux(i1:i2,j1:j2)+sfcFlux(i1:i2,j1:j2)
   w_c%qa(nbeg)%data3d(i1:i2,j1:j2,km) = w_c%qa(nbeg)%data3d(i1:i2,j1:j2,km)+sfcFlux(i1:i2,j1:j2)*cdt* &
                                         (mwtAir/mwtCO)/(w_c%delp(i1:i2,j1:j2,km)/grav)
! Fossil fuel and biofuel
! -----------------------
   sfcFlux(i1:i2,j1:j2) = gcCO%eCO_fosfuel(i1:i2,j1:j2)+gcCO%eCO_biofuel(i1:i2,j1:j2)
   gcCO%COsfcFlux(i1:i2,j1:j2) = gcCO%COsfcFlux(i1:i2,j1:j2)+sfcFlux(i1:i2,j1:j2)
   w_c%qa(nbeg)%data3d(i1:i2,j1:j2,km) = w_c%qa(nbeg)%data3d(i1:i2,j1:j2,km)+sfcFlux(i1:i2,j1:j2)*cdt* &
                                        (mwtAir/mwtCO)/(w_c%delp(i1:i2,j1:j2,km)/grav)
! Release memory
! --------------
   DEALLOCATE(sfcFlux,STAT=ios)

   RETURN
   END SUBROUTINE CO_Emission

 END SUBROUTINE CO_GridCompRun1_

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  CO_GridCompFinalize --- The Chem Driver 
!
! !INTERFACE:
!

   SUBROUTINE CO_GridCompFinalize1_ ( gcCO, w_c, impChem, expChem, &
                                      nymd, nhms, cdt, rc )

! !USES:

  IMPLICIT NONE

! !INPUT/OUTPUT PARAMETERS:

   TYPE(CO_GridComp1), INTENT(INOUT) :: gcCO     ! Grid Component

! !INPUT PARAMETERS:

   TYPE(Chem_Bundle), INTENT(IN)  :: w_c         ! Chemical tracer fields   
   INTEGER, INTENT(IN) :: nymd, nhms             ! time
   REAL,    INTENT(IN) :: cdt                    ! chemical timestep (secs)


! !OUTPUT PARAMETERS:

   TYPE(ESMF_State), INTENT(INOUT) :: impChem    ! Import State
   TYPE(ESMF_State), INTENT(INOUT) :: expChem    ! Import State
   INTEGER, INTENT(OUT) ::  rc                   ! Error return code:
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

   CHARACTER(LEN=*), PARAMETER :: myname = 'CO_GridCompFinalize'
   INTEGER :: ios

   DEALLOCATE ( gcCO%eCO_bioburn, gcCO%eCO_biofuel, gcCO%eCO_fosfuel, & 
                gcCO%COsfcFlux, gcCO%eCO_iso, gcCO%eCO_mon, &
                gcCO%eCO_mtn, gcCO%CH4, gcCO%OHnd, &
                STAT=ios )
   rc = 0
   IF ( ios /= 0 ) rc = 1

   RETURN

 END SUBROUTINE CO_GridCompFinalize1_

 END MODULE CO_GridCompMod

!-----------------------------------------------------------------------

!                     Single Instance Wrapper

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  CO_SingleInstance_ --- Runs single instance of method
!
! !INTERFACE:
!
  subroutine CO_SingleInstance_ ( Method_, instance, &
                                  gcCO, w_c, impChem, expChem, &
                                  nymd, nhms, cdt, rc )

! !USES:

  Use CO_GridCompMod
  Use ESMF
  Use MAPL
  Use Chem_Mod 

  IMPLICIT NONE

! !INPUT PARAMETERS:

!  Input "function pointer"
!  -----------------------
   interface 
     subroutine Method_ (gc, w, imp, exp, ymd, hms, dt, rcode )
       Use CO_GridCompMod
       Use ESMF
       Use MAPL
       Use Chem_Mod 
       type(CO_GridComp1),  intent(inout)  :: gc
       type(Chem_Bundle),   intent(in)     :: w
       type(ESMF_State),    intent(inout)  :: imp
       type(ESMF_State),    intent(inout)  :: exp
       integer,             intent(in)     :: ymd, hms
       real,                intent(in)     :: dt
       integer,             intent(out)    :: rcode
     end subroutine Method_
   end interface

   integer, intent(in)           :: instance     ! instance number

   TYPE(Chem_Bundle), intent(inout) :: w_c       ! Chemical tracer fields      
   INTEGER, INTENT(IN) :: nymd, nhms             ! time
   REAL,    INTENT(IN) :: cdt                    ! chemical timestep (secs)


! !OUTPUT PARAMETERS:

   TYPE(CO_GridComp1), INTENT(INOUT) :: gcCO     ! Grid Component
   TYPE(ESMF_State), INTENT(INOUT)  :: impChem   ! Import State
   TYPE(ESMF_State), INTENT(INOUT)  :: expChem   ! Export State
   INTEGER, INTENT(OUT) ::  rc                   ! Error return code:
                                                 !  0 - all is well
                                                 !  1 - 

! !DESCRIPTION: Finalizes the CO Grid Component. Multiple instance
!               version.
!
! !REVISION HISTORY:
!
!  27Feb2008  da Silva  Introduced multiple instances
!
!EOP
!-------------------------------------------------------------------------

  integer n_CO, i_CO, j_CO

! Save overall CO indices
! -----------------------
  n_CO = w_c%reg%n_CO
  i_CO = w_c%reg%i_CO
  j_CO = w_c%reg%j_CO
  
! Customize indices for this particular instance
! ----------------------------------------------
  w_c%reg%n_CO = 1
  w_c%reg%i_CO = i_CO + instance - 1
  w_c%reg%j_CO = i_CO + instance - 1
  
! Execute the instance method
! ---------------------------
  call Method_ ( gcCO, w_c, impChem, expChem, &
                 nymd, nhms, cdt, rc )

! Restore the overall CO indices
! ------------------------------
  w_c%reg%n_CO = n_CO
  w_c%reg%i_CO = i_CO
  w_c%reg%j_CO = j_CO

  end subroutine CO_SingleInstance_

!-----------------------------------------------------------------------
