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
   USE Chem_UtilMod                    ! I/O

!  bweir: for photolysis
   USE ESMF_CFIOFileMOD
   USE MAPL_CFIOMOD

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

        REAL, POINTER ::       CH4(:,:,:)    ! CH4 mixing ratio (mol/mol)
        REAL, POINTER ::      OHnd(:,:,:)    ! OH number density (#/cm3)
        REAL, POINTER ::      Clnd(:,:,:)    ! Cl number density (#/cm3)
        REAL, POINTER ::     O1Dnd(:,:,:)    ! O(1D) number density (#/cm3)

        REAL, POINTER :: COsfcFlux(:,:)      ! CO surface flux kg m^-2 s^-1

        LOGICAL :: DBG                       ! Run-time debug switch
        LOGICAL :: doingBB = .true.          ! Switch to consider biomass burning
        CHARACTER(LEN=ESMF_MAXSTR) :: units_oh ! Units for OH 

!       Photolysis tables (bweir: from StratChem)
!       -----------------
        INTEGER :: numphoto
        INTEGER :: nxdo
        INTEGER :: nlam
        INTEGER :: nsza
        INTEGER :: numo3
        INTEGER :: nts
        INTEGER :: aqsize

        REAL, POINTER :: sdat(:,:,:,:)
        REAL, POINTER :: o2jdat(:,:,:)
        REAL, POINTER :: sza_tab(:)
        REAL, POINTER :: o3_tab(:,:)
        REAL, POINTER :: xtab(:,:,:)
        REAL, POINTER :: CH2O_aq(:)
        REAL, POINTER :: rlam(:)
  END TYPE CO_GridComp1

  TYPE CO_GridComp
     INTEGER                     ::  n        ! number of instances 
     TYPE(CO_GridComp1), POINTER ::  gcs(:)   ! instances
  END TYPE CO_GridComp

! REAL, PARAMETER :: radToDeg = 57.2957795
  REAL, PARAMETER :: radToDeg = 180./MAPL_PI
  REAL, PARAMETER :: mwtCO    = 28.0104

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
!  CO bins in the registry (it is OK to have less, though)
!  -------------------------------------------------------
   if ( n .LT. chemReg%n_CO ) then
        rc = 35
        return
   else if ( n .GT. chemReg%n_CO ) then
        if (MAPL_AM_I_ROOT()) &
        print *, trim(Iam)// &
                 ': fewer CO bins than possible CO instances: ',&
                 n, chemReg%n_CO
   end if
   n = min(n, chemReg%n_CO)

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
   
   integer :: i, ier, n
   real :: c1, c2, c3, c4

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
!  CO bins in the registry (it is OK to have less, though)
!  -------------------------------------------------------
   if ( n .LT. w_c%reg%n_CO ) then
        rc = 35
        return
   else if ( n .GT. w_c%reg%n_CO ) then
        if (MAPL_AM_I_ROOT()) &
        print *, trim(myname)// &
                 ': fewer CO bins than possible CO instances: ',&
                 n, w_c%reg%n_CO
   end if
   n = min(n, w_c%reg%n_CO)
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
   call get_HenrysLawCts('CO', c1, c2, c3, c4)
   w_c%reg%Hcts(1,w_c%reg%i_CO:w_c%reg%j_CO) = c1
   w_c%reg%Hcts(2,w_c%reg%i_CO:w_c%reg%j_CO) = c2
   w_c%reg%Hcts(3,w_c%reg%i_CO:w_c%reg%j_CO) = c3
   w_c%reg%Hcts(4,w_c%reg%i_CO:w_c%reg%j_CO) = c4

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
       LONG_NAME  = 'source species',     &
       UNITS      = '1',                  &
       DIMS       = MAPL_DimsHorzVert,    &
       VLOCATION  = MAPL_VLocationCenter, &
       RESTART    = MAPL_RestartSkip,     &
       RC         = STATUS)
  VERIFY_(STATUS)
  call MAPL_AddImportSpec(GC,             &
       SHORT_NAME = 'CO_Cl'//iname,       &
       LONG_NAME  = 'source species',     &
       UNITS      = '1',                  &
       DIMS       = MAPL_DimsHorzVert,    &
       VLOCATION  = MAPL_VLocationCenter, &
       RESTART    = MAPL_RestartSkip,     &
       RC         = STATUS)
  VERIFY_(STATUS)
  call MAPL_AddImportSpec(GC,             &
       SHORT_NAME = 'CO_O1D'//iname,      &
       LONG_NAME  = 'source species',     &
       UNITS      = '1',                  &
       DIMS       = MAPL_DimsHorzVert,    &
       VLOCATION  = MAPL_VLocationCenter, &
       RESTART    = MAPL_RestartSkip,     &
       RC         = STATUS)
  VERIFY_(STATUS)
  call MAPL_AddImportSpec(GC,             &
       SHORT_NAME = 'CO_CH4'//iname,      &
       LONG_NAME  = 'source species',     &
       UNITS      = '1',                  &
       DIMS       = MAPL_DimsHorzVert,    &
       VLOCATION  = MAPL_VLocationCenter, &
       RESTART    = MAPL_RestartSkip,     &
       RC         = STATUS)
  VERIFY_(STATUS)
  call MAPL_AddImportSpec(GC,           &
       SHORT_NAME = 'CO_BF'//iname,     &
       LONG_NAME  = 'source species',   &
       UNITS      = '1',                &
       DIMS       = MAPL_DimsHorzOnly,  &
       VLOCATION  = MAPL_VLocationNone, &
       RESTART    = MAPL_RestartSkip,   &
       RC         = STATUS)
  VERIFY_(STATUS)
  call MAPL_AddImportSpec(GC,           &
       SHORT_NAME = 'CO_FS'//iname,     &
       LONG_NAME  = 'source species',   &
       UNITS      = '1',                &
       DIMS       = MAPL_DimsHorzOnly,  &
       VLOCATION  = MAPL_VLocationNone, &
       RESTART    = MAPL_RestartSkip,   &
       RC         = STATUS)
  VERIFY_(STATUS)
  call MAPL_AddImportSpec(GC,           &
       SHORT_NAME = 'CO_ISOP'//iname,   &
       LONG_NAME  = 'source species',   &
       UNITS      = '1',                &
       DIMS       = MAPL_DimsHorzOnly,  &
       VLOCATION  = MAPL_VLocationNone, &
       RESTART    = MAPL_RestartSkip,   &
       RC         = STATUS)
  VERIFY_(STATUS)
  call MAPL_AddImportSpec(GC,           &
       SHORT_NAME = 'CO_NVOC'//iname,   &
       LONG_NAME  = 'source species',   &
       UNITS      = '1',                &
       DIMS       = MAPL_DimsHorzOnly,  &
       VLOCATION  = MAPL_VLocationNone, &
       RESTART    = MAPL_RestartSkip,   &
       RC         = STATUS)
  VERIFY_(STATUS)
  call MAPL_AddImportSpec(GC,           &
       SHORT_NAME = 'CO_TERP'//iname,   &
       LONG_NAME  = 'source species',   &
       UNITS      = '1',                &
       DIMS       = MAPL_DimsHorzOnly,  &
       VLOCATION  = MAPL_VLocationNone, &
       RESTART    = MAPL_RestartSkip,   &
       RC         = STATUS)
  VERIFY_(STATUS)
  call MAPL_AddImportSpec(GC,           &
       SHORT_NAME = 'CO_BIOMASS'//iname,&
       LONG_NAME  = 'source species',   &
       UNITS      = '1',                &
       DIMS       = MAPL_DimsHorzOnly,  &
       VLOCATION  = MAPL_VLocationNone, &
       RESTART    = MAPL_RestartSkip,   &
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
   CHARACTER(LEN=*), PARAMETER :: Iam = myname

   CHARACTER(LEN=255) :: rcfilen 

   INTEGER :: ios, j, n
   INTEGER, ALLOCATABLE :: ier(:)
   INTEGER :: i1, i2, im, j1, j2, jm, km
   INTEGER :: nTimes, begTime, incSecs
   INTEGER :: nbeg, nend, nymd1, nhms1

   REAL :: limitN, limitS
   REAL, ALLOCATABLE :: var2d(:,:)

!  Photolysis (bweir: from StratChem)
   CHARACTER(LEN=ESMF_MAXSTR) :: fnphoto

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

!  Possibly oxidants are in a different unit
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

!  Read photolysis tables (bweir: from StratChem)
!  ----------------------
   CALL I90_label('photolysisFile:', ios)
   IF( ios == 0 ) THEN
      gcCO%numphoto = 55
      CALL I90_Gtoken(fnphoto, ios)
      VERIFY_(ios)
      CALL readPhotTables(trim(fnphoto), ios)
      VERIFY_(ios)
   ELSE
      gcCO%numphoto = 0
   END IF

   RETURN

CONTAINS

   SUBROUTINE init_()

   INTEGER ios, nerr
   nerr = 128
   ALLOCATE ( gcCO%eCO_bioburn(i1:i2,j1:j2),  & 
              gcCO%eCO_bioburn_(i1:i2,j1:j2), & 
              gcCO%eCO_biofuel(i1:i2,j1:j2),  &
              gcCO%eCO_fosfuel(i1:i2,j1:j2),  &
              gcCO%COsfcFlux(i1:i2,j1:j2),    &
              gcCO%eCO_iso(i1:i2,j1:j2),      &
              gcCO%eCO_mon(i1:i2,j1:j2),      &
              gcCO%eCO_mtn(i1:i2,j1:j2),      &
              gcCO%CH4(i1:i2,j1:j2,km),       &
              gcCO%OHnd(i1:i2,j1:j2,km),      &
              gcCO%Clnd(i1:i2,j1:j2,km),      &
              gcCO%O1Dnd(i1:i2,j1:j2,km),     &
              ier(nerr), STAT=ios )
   IF ( ios /= 0 ) rc = 100
   END SUBROUTINE init_

   SUBROUTINE final_(ierr)
   INTEGER :: ierr
   INTEGER ios
   DEALLOCATE ( gcCO%eCO_bioburn, gcCO%eCO_biofuel, gcCO%eCO_fosfuel,  & 
                gcCO%eCO_bioburn_, gcCO%COsfcFlux, gcCO%eCO_iso,       &
                gcCO%eCO_mon, gcCO%eCO_mtn, gcCO%CH4, gcCO%OHnd,       &
                gcCO%Clnd, gcCO%O1Dnd, ier, STAT=ios )
!  Photolysis (bweir: from StratChem)
   DEALLOCATE ( gcCO%sdat, gcCO%o2jdat, gcCO%o3_tab, gcCO%xtab,        & 
                gcCO%sza_tab, gcCO%CH2O_aq, gcCO%rlam, ier, STAT=ios )
   CALL I90_release()
   rc = ierr
   END SUBROUTINE final_

   SUBROUTINE readPhotTables(fileName, rc)

   IMPLICIT NONE

!  Read tables for photolysis in GOCART ... from a NetCDF file
!
!  Input parameters:
!
   CHARACTER(LEN=*), INTENT(IN) :: fileName
!
!  Output parameters:
!
   INTEGER, INTENT(OUT) :: rc
!
!  Restrictions:
!  ASSERT that the number of pressure layers in the dataset equals km.
!
!  REVISION HISTORY:
!  Nielsen     11 May 2012: First crack.
!  Weir        29 Jan 2021: Pilferd from StratChem
!-----------------------------------------------------------------------

  CHARACTER(LEN=ESMF_MAXSTR) :: Iam = "CO::readPhotTables"

  TYPE(ESMF_VM) :: vm

  INTEGER :: comm, info, unit, status
  INTEGER :: dimid, i, n

  INTEGER :: length

  INTEGER, PARAMETER :: nD = 7
  CHARACTER(LEN=ESMF_MAXSTR) :: dimName(nD)= (/"nsza  ", "numO3 ", "layers", &
                                               "nlam  ", "nts   ", "nxdo  ", "aqsize" /)

  INTEGER, PARAMETER :: nV = 7
  CHARACTER(LEN=ESMF_MAXSTR) :: varName(nV)= (/"sza    ", &
                        "lambda ", "O3TAB  ",  "SDAT   ", &
                        "O2JDAT ", "XTAB   ",  "CH2O_AQ" /)
  rc = 0

! Grab the virtual machine
! ------------------------
  CALL ESMF_VMGetCurrent(vm, RC=status)
  VERIFY_(status)

  CALL ESMF_VMGet(vm, MPICOMMUNICATOR=comm, rc=status)
  VERIFY_(status)

#ifdef H5_HAVE_PARALLEL

  CALL MPI_Info_create(info, status)
  VERIFY_(status)
  CALL MPI_Info_set(info, "romio_cb_read", "automatic", status)
  VERIFY_(status)

#ifdef NETCDF_NEED_NF_MPIIO
  status = NF_OPEN_PAR(TRIM(fileName), IOR(NF_NOWRITE,NF_MPIIO), comm, info, unit)
#else
  status = NF_OPEN_PAR(TRIM(fileName), NF_NOWRITE, comm, info, unit)
#endif

#else

  IF(MAPL_AM_I_ROOT(vm)) THEN 
   status = NF_OPEN(TRIM(fileName), NF_NOWRITE, unit)

#endif

   IF(status /= NF_NOERR) THEN
    PRINT *,'Error opening file ',TRIM(fileName), status
    PRINT *, NF_STRERROR(status)
    VERIFY_(status)
   END IF

   DO i = 1,nD

    status = NF_INQ_DIMID(unit, TRIM(dimName(i)), dimid)
    IF(status /= NF_NOERR) THEN
     PRINT *,"Error inquiring dimension ID for ", TRIM(dimName(i)), status
     PRINT *, NF_STRERROR(status)
     VERIFY_(status)
    END IF

    status = NF_INQ_DIMLEN(unit, dimid, n)
    IF(status /= NF_NOERR) THEN
     PRINT *,"Error inquiring  dimension length for ", TRIM(dimName(i)), status
     PRINT *, NF_STRERROR(status)
    END IF

    SELECT CASE (i)
     CASE (1)
      gcCO%nsza = n
     CASE (2)
      gcCO%numO3 = n
     CASE (3)
      ASSERT_(n == km)
     CASE (4)
      gcCO%nlam = n
     CASE (5)
      gcCO%nts = n
     CASE (6)
      gcCO%nxdo = n
     CASE (7)
      gcCO%aqsize = n
     CASE DEFAULT
    END SELECT

   END DO

#ifndef H5_HAVE_PARALLEL

  END IF ! MAPL_AM_I_ROOT

  CALL MAPL_CommsBcast(vm, gcCO%nsza, 1, 0, RC=status)
  VERIFY_(status)
  CALL MAPL_CommsBcast(vm, gcCO%numO3, 1, 0, RC=status)
  VERIFY_(status)
  CALL MAPL_CommsBcast(vm, gcCO%nlam, 1, 0, RC=status)
  VERIFY_(status)
  CALL MAPL_CommsBcast(vm, gcCO%nts, 1, 0, RC=status)
  VERIFY_(status)
  CALL MAPL_CommsBcast(vm, gcCO%nxdo, 1, 0, RC=status)
  VERIFY_(status)
  CALL MAPL_CommsBcast(vm, gcCO%aqSize, 1, 0, RC=status)
  VERIFY_(status)

#endif

  ALLOCATE(gcCO%sdat(gcCO%nsza,gcCO%numo3,km,gcCO%nlam), STAT=status)
  VERIFY_(status)
  ALLOCATE(gcCO%o2jdat(gcCO%nsza,gcCO%numo3,km), STAT=status)
  VERIFY_(status)
  ALLOCATE(gcCO%o3_tab(gcCO%numo3,km), STAT=status)
  VERIFY_(status)
  ALLOCATE(gcCO%xtab(gcCO%nlam,gcCO%nxdo,gcCO%nts), STAT=status)
  VERIFY_(status)
  ALLOCATE(gcCO%sza_tab(gcCO%nsza), STAT=status)
  VERIFY_(status)
  ALLOCATE(gcCO%CH2O_aq(gcCO%aqSize), STAT=status)
  VERIFY_(status)
  ALLOCATE(gcCO%rlam(gcCO%nlam), STAT=status)
  VERIFY_(status)

#ifndef H5_HAVE_PARALLEL

  IF(MAPL_AM_I_ROOT()) THEN

#endif

   DO i = 1,nV

    status = NF_INQ_VARID(unit, TRIM(varName(i)), n)
    IF(status /= NF_NOERR) THEN
     PRINT *,"Error getting varid for ", TRIM(varName(i)), status
     PRINT *, NF_STRERROR(status)
     VERIFY_(status)
    END IF

    SELECT CASE (i)
     CASE (1)
      status = NF_GET_VAR_REAL(unit, n, gcCO%sza_tab)
     CASE (2)
      status = NF_GET_VAR_REAL(unit, n, gcCO%rlam)
     CASE (3)
      status = NF_GET_VAR_REAL(unit, n, gcCO%o3_tab)
     CASE (4)
      status = NF_GET_VAR_REAL(unit, n, gcCO%sdat)
     CASE (5)
      status = NF_GET_VAR_REAL(unit, n, gcCO%o2jdat)
     CASE (6)
      status = NF_GET_VAR_REAL(unit, n, gcCO%xtab)
     CASE (7)
      status = NF_GET_VAR_REAL(unit, n, gcCO%CH2O_aq)
     CASE DEFAULT
    END SELECT

    IF(status /= NF_NOERR) THEN
     PRINT *,"Error getting values for ", TRIM(varName(i)), status
     PRINT *, NF_STRERROR(status)
     VERIFY_(status)
    END IF

   END DO

#ifdef H5_HAVE_PARALLEL

   CALL MPI_Info_free(info, status)
   VERIFY_(status)

#else

  END IF ! MAPL_AM_I_ROOT

  length = SIZE(gcCO%sza_tab)
  CALL MPI_Bcast(gcCO%sza_tab, length, MPI_REAL, 0, comm, status)
  VERIFY_(status)

  length = SIZE(gcCO%rlam)
  CALL MPI_Bcast(gcCO%rlam, length, MPI_REAL, 0, comm, status)
  VERIFY_(status)

  length = SIZE(gcCO%o3_tab)
  CALL MPI_Bcast(gcCO%o3_tab, length, MPI_REAL, 0, comm, status)
  VERIFY_(status)

  length = SIZE(gcCO%sdat)
  CALL MPI_Bcast(gcCO%sdat, length, MPI_REAL, 0, comm, status)
  VERIFY_(status)

  length = SIZE(gcCO%o2jdat)
  CALL MPI_Bcast(gcCO%o2jdat, length, MPI_REAL, 0, comm, status)
  VERIFY_(status)

  length = SIZE(gcCO%xtab)
  CALL MPI_Bcast(gcCO%xtab, length, MPI_REAL, 0, comm, status)
  VERIFY_(status)

  CALL MAPL_CommsBcast(vm, gcCO%CH2O_aq, gcCO%aqsize, 0, RC=status)
  VERIFY_(status)

#endif

  status = NF_CLOSE(unit)
  VERIFY_(status)

  RETURN
 END SUBROUTINE readPhotTables

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
   REAL, POINTER, DIMENSION(:,:)   ::  pblh   => null()
   REAL, POINTER, DIMENSION(:,:,:) ::  T      => null()
   REAL, POINTER, DIMENSION(:,:,:) ::  rhowet => null()
   REAL, POINTER, DIMENSION(:,:,:) ::  zle    => null()
   REAL, POINTER, DIMENSION(:,:,:) ::  qtot   => null()

   INTEGER :: i1, i2, im, j1, j2, jm, km, ios, idiag, iXj
   INTEGER :: i, j, k, kReverse, n, nbeg, nend
   INTEGER :: nymd1, nhms1, ier(8)
   integer :: iregWant

   REAL    :: qmin, qmax
   REAL    :: fiso, fmtn, fmon

   REAL, ALLOCATABLE :: pe(:,:,:), p(:,:,:), ndwet(:,:,:)
   REAL, ALLOCATABLE :: rkoh(:,:,:), rkch4_oh(:,:,:)
   REAL, ALLOCATABLE :: rkch4_cl(:,:,:), rkch4_o1d(:,:,:)

!  Photolysis (bweir: from StratChem, but aj is SINGLE)
!  ----------
   REAL, ALLOCATABLE :: photJ(:,:,:), dCOPhot(:,:,:)
   REAL, ALLOCATABLE :: aj(:)
   REAL    :: szan

   real, pointer, dimension(:,:,:) :: ptr3d => null()
   real, pointer, dimension(:,:)   :: ptr2d => null()

   REAL, ALLOCATABLE, DIMENSION(:,:) :: psdry, psco

#define EXPORT   expChem
#define iNAME    TRIM(gcCO%iname)

#define COEM     CO_emis
#define COSC     CO_surface
#define COCL     CO_column
#define CODRY	 CO_dry
#define COSD     CO_surfdry
#define COCD     CO_coldry
#define COPD     CO_prod
#define COLS     CO_loss
#define COJP	 CO_phot

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
   if (nbeg /= nend) then
      if (MAPL_AM_I_ROOT()) print *, myname, ": Must have only 1 bin at the single instance level"
      rc = 1
      return 
   endif

!  Biomass Burning
!  ---------------
   call MAPL_GetPointer(impChem, ptr2d, 'CO_BIOMASS'//iNAME, rc=status)
   VERIFY_(STATUS)
   gcCO%eCO_bioburn = ptr2d

!  Biofuel source
!  --------------
   call MAPL_GetPointer(impChem, ptr2d, 'CO_BF'//iNAME, rc=status)
   VERIFY_(STATUS)
   gcCO%eCO_biofuel = ptr2d

!  Fossil fuel source
!  ------------------
   call MAPL_GetPointer(impChem, ptr2d, 'CO_FS'//iNAME, rc=status)
   VERIFY_(STATUS)
   gcCO%eCO_fosfuel = ptr2d

!  Background OH, for loss term
!  ----------------------------
   call MAPL_GetPointer(impChem, ptr3d, 'CO_OH'//iNAME, rc=status)
   VERIFY_(STATUS)
   gcCO%OHnd = ptr3d

!  Background Cl, for loss term
!  ----------------------------
   call MAPL_GetPointer(impChem, ptr3d, 'CO_Cl'//iNAME, rc=status)
   VERIFY_(STATUS)
   gcCO%Clnd = ptr3d

!  Background O1D, for loss term
!  ----------------------------
   call MAPL_GetPointer(impChem, ptr3d, 'CO_O1D'//iNAME, rc=status)
   VERIFY_(STATUS)
   gcCO%O1Dnd = ptr3d

!  Background CH4, for source term
!  NOTE: Return zeroes in all but the global instantiation
!  NOTE: Not sure this NOTE is true anymore
!  -------------------------------------------------------
   call MAPL_GetPointer(impChem, ptr3d, 'CO_CH4'//iNAME, rc=status)
   VERIFY_(STATUS)
   gcCO%CH4 = ptr3d

!  Isoprene source
!  ---------------
   call MAPL_GetPointer(impChem, ptr2d, 'CO_ISOP'//iNAME, rc=status)
   VERIFY_(STATUS)
   gcCO%eCO_iso = ptr2d

!  VOC source
!  ----------
   call MAPL_GetPointer(impChem, ptr2d, 'CO_NVOC'//iNAME, rc=status)
   VERIFY_(STATUS)
   gcCO%eCO_mon = ptr2d

!  Monoterpene source
!  ------------------
   call MAPL_GetPointer(impChem, ptr2d, 'CO_TERP'//iNAME, rc=status)
   VERIFY_(STATUS)
   gcCO%eCO_mtn = ptr2d

   if (gcCO%DBG) then
      call pmaxmin('CO: eCO_bioburn', gcCO%eCO_bioburn, qmin, qmax, iXj,1, 1. )
      call pmaxmin('CO: eCO_biofuel', gcCO%eCO_biofuel, qmin, qmax, iXj,1, 1. )
      call pmaxmin('CO: eCO_fosfuel', gcCO%eCO_fosfuel, qmin, qmax, iXj,1, 1. )
      call pmaxmin('CO: eCO_iso',     gcCO%eCO_iso,     qmin, qmax, iXj,1, 1. )
      call pmaxmin('CO: eCO_mon',     gcCO%eCO_mon,     qmin, qmax, iXj,1, 1. )
      call pmaxmin('CO: eCO_mtn',     gcCO%eCO_mtn,     qmin, qmax, iXj,1, 1. )
   endif

!  Apply diurnal cycle if so desired
!  ---------------------------------
   if (w_c%diurnal_bb) then
      gcCO%eCO_bioburn_(:,:) = gcCO%eCO_bioburn(:,:)

      call Chem_BiomassDiurnal(gcCO%eCO_bioburn, gcCO%eCO_bioburn_,   &
                               w_c%grid%lon(:,:)*radToDeg,            &
                               w_c%grid%lat(:,:)*radToDeg, nhms, cdt)      
   endif

!  Allocate temporary workspace
!  ----------------------------
   allocate(pe(i1:i2,j1:j2,km+1), p(i1:i2,j1:j2,km), ndwet(i1:i2,j1:j2,km), &
            rkoh(i1:i2,j1:j2,km), rkch4_oh(i1:i2,j1:j2,km),                 &
            rkch4_cl(i1:i2,j1:j2,km), rkch4_o1d(i1:i2,j1:j2,km), stat = ios )

   if (ios /= 0) then
      rc = 3
      return
   endif

!  Layer interface pressures
!  -------------------------
   pe(i1:i2,j1:j2,1) = w_c%grid%ptop
   do k=2,km+1
      pe(i1:i2,j1:j2,k) = pe(i1:i2,j1:j2,k-1) + w_c%delp(i1:i2,j1:j2,k-1)
   enddo

!  Layer mean pressures
!  --------------------
   do k=1,km
      p(i1:i2,j1:j2,k) = (pe(i1:i2,j1:j2,k)+pe(i1:i2,j1:j2,k+1))*0.50
   enddo
 
!  Get imports
!  -----------
   call MAPL_GetPointer( impChem, pblh,   'ZPBL',    rc=ier(1) ) 
   call MAPL_GetPointer( impChem, T,      'T',       rc=ier(2) ) 
   call MAPL_GetPointer( impChem, rhowet, 'AIRDENS', rc=ier(3) ) 
   call MAPL_GetPointer( impChem, zle,    'ZLE',     rc=ier(4) ) 
   call MAPL_GetPointer( impChem, qtot,   'QTOT',    rc=ier(5) ) 

   if (any(ier(1:5) /= 0)) then
      rc = 10
      return
   endif

   if (gcCO%DBG) then
      call pmaxmin('CO:PBLH',     pblh, qmin, qmax, iXj,    1, 1. )
      call pmaxmin('CO:T',           T, qmin, qmax, iXj,   km, 1. )
      call pmaxmin('CO:RHOWET', rhowet, qmin, qmax, iXj,   km, 1. )
      call pmaxmin('CO:ZLE',       zle, qmin, qmax, iXj, km+1, 1. )
      call pmaxmin('CO:QTOT',     qtot, qmin, qmax, iXj,   km, 1. )
   endif

!  Wet-air number density
!  ----------------------
   ndwet(i1:i2,j1:j2,1:km) = rhowet(i1:i2,j1:j2,1:km)*MAPL_AVOGAD/MAPL_AIRMW

!  Handle mole fraction or number density units of oxidants
!  --------------------------------------------------------
   if (trim(gcCO%units_oh) == 'mol/mol' .or. trim(gcCO%units_oh) == 'mol mol-1') then
       gcCO%OHnd(i1:i2,j1:j2,1:km)  = gcCO%OHnd(i1:i2,j1:j2,1:km) &
                                    * ndwet(i1:i2,j1:j2,1:km)
       gcCO%Clnd(i1:i2,j1:j2,1:km)  = gcCO%Clnd(i1:i2,j1:j2,1:km) &
                                    * ndwet(i1:i2,j1:j2,1:km)
       gcCO%O1Dnd(i1:i2,j1:j2,1:km) = gcCO%O1Dnd(i1:i2,j1:j2,1:km) &
                                    * ndwet(i1:i2,j1:j2,1:km)
   else
!      Otherwise, assume units are molec cm^-3 and convert to molec m^-3
       gcCO%OHnd(i1:i2,j1:j2,1:km)  =  gcCO%OHnd(i1:i2,j1:j2,1:km)*1.00E+06
       gcCO%Clnd(i1:i2,j1:j2,1:km)  =  gcCO%Clnd(i1:i2,j1:j2,1:km)*1.00E+06
       gcCO%O1Dnd(i1:i2,j1:j2,1:km) = gcCO%O1Dnd(i1:i2,j1:j2,1:km)*1.00E+06
   end if

!  Loss due to OH
!  --------------
   rkoh(i1:i2,j1:j2,1:km) = 1.50E-13*1.00E-06*(1.00+0.60E-05*p(i1:i2,j1:j2,1:km))

   if (associated(CO_loss)) then
      CO_loss(i1:i2,j1:j2) = 0.
      do k = 1,km
         CO_loss(i1:i2,j1:j2) = CO_loss(i1:i2,j1:j2)                                        &
                              + w_c%qa(nbeg)%data3d(i1:i2,j1:j2,k)*mwtCO/MAPL_AIRMW         &
                                         *     rkoh(i1:i2,j1:j2,k)*gcCO%OHnd(i1:i2,j1:j2,k) &
                                         * w_c%delp(i1:i2,j1:j2,k)/MAPL_GRAV
      enddo
   endif

!  Production due to CH4
!  ---------------------
   rkch4_oh(i1:i2,j1:j2,1:km)  = 2.45E-12*1.00E-06*exp(-1775./t(i1:i2,j1:j2,1:km))
   rkch4_cl(i1:i2,j1:j2,1:km)  = 7.10E-12*1.00E-06*exp(-1270./t(i1:i2,j1:j2,1:km))
   rkch4_o1d(i1:i2,j1:j2,1:km) = 1.75E-10*1.00E-06

   if (associated(CO_prod)) then
      CO_prod(i1:i2,j1:j2) = 0.
      do k = 1,km
         CO_prod(i1:i2,j1:j2) = CO_prod(i1:i2,j1:j2)                                    &
                              + (    rkch4_oh(i1:i2,j1:j2,k)* gcCO%OHnd(i1:i2,j1:j2,k)  &
                                  +  rkch4_cl(i1:i2,j1:j2,k)* gcCO%Clnd(i1:i2,j1:j2,k)  &
                                  + rkch4_o1d(i1:i2,j1:j2,k)*gcCO%O1Dnd(i1:i2,j1:j2,k)) &
                                * gcCO%CH4(i1:i2,j1:j2,k)*mwtCO/MAPL_AIRMW              &
                                * w_c%delp(i1:i2,j1:j2,k)/MAPL_GRAV
      enddo
   endif

!  Calculate photolytic loss rates, J [s^-1] for
!     CH4 + hv => 2H2O + CO
!     CO2 + hv => CO + ???
!  Notice that J and the losses are always computed. However, the setting 
!  of the feedback switch(es) determines if the increments are actually applied
!  ----------------------------------------------------------------------------
   allocate(photJ(i1:i2,j1:j2,1:km), dCOPhot(i1:i2,j1:j2,1:km), STAT=status)
   VERIFY_(STATUS)

   photJ(i1:i2,j1:j2,1:km) = 0.
   dCOPhot(i1:i2,j1:j2,1:km) = 0.

!  Change in CO due to CH4 photolysis
!  ----------------------------------
!  call getJRates(status)
!  VERIFY_(status)
!
!  dCOPhot(i1:i2,j1:j2,1:km) = photJ(i1:i2,j1:j2,1:km) * CH4(i1:i2,j1:j2,1:km)

!  Change in CO due to CO2 photolysis
!  ----------------------------------
   if (gcCO%numphoto > 0) then
      photJ(i1:i2,j1:j2,1:km) = 0.
      allocate(aj(gcCO%numphoto), STAT=status)
      VERIFY_(STATUS)

      do k = 1,km
         do j = j1,j2
            do i = i1,i2
               szan = 0.
               if (w_c%cosz(i,j) <= 1.) szan = acos(w_c%cosz(i,j))
!              bweir: Using 0 for O3 (FIXME)
               call jcalc4(km-k+1, szan, 0., p(i,j,k), t(i,j,k), aj, gcCO)
               photJ(i,j,k) = aj(12)
            enddo
         enddo
      enddo
   endif

!  bweir: Using 400 ppm for CO2 (FIXME)
!  dCOPhot(i1:i2,j1:j2,1:km) = dCOPhot(i1:i2,j1:j2,1:km) + photJ(i1:i2,j1:j2,1:km) * CO2(i1:i2,j1:j2,1:km)
   dCOPhot(i1:i2,j1:j2,1:km) = dCOPhot(i1:i2,j1:j2,1:km) + photJ(i1:i2,j1:j2,1:km) * 400.e-6

!  Photolysis
!  ----------
   if (associated(CO_phot)) then
      CO_phot(i1:i2,j1:j2,1:km) = dCOPhot(i1:i2,j1:j2,1:km)
   endif

!  Decrement the CO mole fraction due to oxidation 
!  -----------------------------------------------
   w_c%qa(nbeg)%data3d(i1:i2,j1:j2,1:km) = w_c%qa(nbeg)%data3d(i1:i2,j1:j2,1:km) &
         - cdt * rkoh(i1:i2,j1:j2,1:km)*gcCO%OHnd(i1:i2,j1:j2,1:km)              &
               * w_c%qa(nbeg)%data3d(i1:i2,j1:j2,1:km)

!  bweir: just checking
   w_c%qa(nbeg)%data3d(i1:i2,j1:j2,1:km) = max(w_c%qa(nbeg)%data3d(i1:i2,j1:j2,1:km), 0.)

!  Compute and add surface emissions
!  ---------------------------------
   gcCO%COsfcFlux(i1:i2,j1:j2) = 0.
   call CO_Emission(rc)

!  Increment the CO mole fraction due to photolysis
!  ------------------------------------------------
   w_c%qa(nbeg)%data3d(i1:i2,j1:j2,1:km) = w_c%qa(nbeg)%data3d(i1:i2,j1:j2,1:km) + cdt*dCOPhot(i1:i2,j1:j2,1:km)

!  Increment the CO mole fraction due to production
!  ------------------------------------------------
   w_c%qa(nbeg)%data3d(i1:i2,j1:j2,1:km) = w_c%qa(nbeg)%data3d(i1:i2,j1:j2,1:km) &
         + cdt * (    rkch4_oh(i1:i2,j1:j2,1:km)* gcCO%OHnd(i1:i2,j1:j2,1:km)    &
                   +  rkch4_cl(i1:i2,j1:j2,1:km)* gcCO%Clnd(i1:i2,j1:j2,1:km)    &
                   + rkch4_o1d(i1:i2,j1:j2,1:km)*gcCO%O1Dnd(i1:i2,j1:j2,1:km))   &
               * gcCO%CH4(i1:i2,j1:j2,1:km)

!  Surface concentration [ppbv]
!  ----------------------------
   if (associated(CO_surface)) then
      CO_surface(i1:i2,j1:j2) = w_c%qa(nbeg)%data3d(i1:i2,j1:j2,km)*1.e9
   endif

!  Column burden [kg m-2]
!  ----------------------
   if (associated(CO_column)) then
      CO_column(i1:i2,j1:j2) = 0.
      do k = 1, km
         CO_column(i1:i2,j1:j2) = CO_column(i1:i2,j1:j2)                                &
                                + w_c%qa(nbeg)%data3d(i1:i2,j1:j2,k)*mwtCO/MAPL_AIRMW * &
                                             w_c%delp(i1:i2,j1:j2,k)/MAPL_GRAV
     enddo
   endif

!  Dry-air mole fraction
!  ---------------------
   if (associated(CO_dry)) then
      CO_dry(i1:i2,j1:j2,1:km) = w_c%qa(nbeg)%data3d(i1:i2,j1:j2,1:km) &
                                        / (1. - qtot(i1:i2,j1:j2,1:km))
   endif

!  Dry-air surface concentration [mol mol-1]
!  -----------------------------------------
   if (associated(CO_surfdry)) then
      CO_surfdry(i1:i2,j1:j2) = w_c%qa(nbeg+n-1)%data3d(i1:i2,j1:j2,km)*(1. - qtot(i1:i2,j1:j2,km))
   endif

!  Dry-air column average [mol mol-1]
!  ----------------------------------
   allocate(psdry(i1:i2,j1:j2), stat=status)    ! dry-air surface pressure
   allocate(psco( i1:i2,j1:j2), stat=status)    ! co      surface pressure

   if (associated(CO_coldry)) then
      psdry(i1:i2,j1:j2) = 0.
      psco( i1:i2,j1:j2) = 0.
      do k = 1,km
         psdry(i1:i2,j1:j2) = psdry(i1:i2,j1:j2) + w_c%delp(i1:i2,j1:j2,k)*(1. - qtot(i1:i2,j1:j2,k))
         psco( i1:i2,j1:j2) = psco( i1:i2,j1:j2) + w_c%delp(i1:i2,j1:j2,k)*w_c%qa(nbeg+n-1)%data3d(i1:i2,j1:j2,k)
      enddo
      CO_coldry(i1:i2,j1:j2) = psco(i1:i2,j1:j2)/psdry(i1:i2,j1:j2)
   endif

   deallocate(psdry, psco)

!  CO Surface Emission Flux in kg m-2 s-1
!  --------------------------------------
   if (associated(CO_emis)) CO_emis(i1:i2,j1:j2) = gcCO%COsfcFlux(i1:i2,j1:j2)

   if (gcCO%DBG) then
      if (associated(CO_emis))    call pmaxmin('CO: emis',       CO_emis, qmin, qmax, iXj,  1, 1. )
      if (associated(CO_loss))    call pmaxmin('CO: loss',       CO_loss, qmin, qmax, iXj,  1, 1. )
      if (associated(CO_prod))    call pmaxmin('CO: prod',       CO_prod, qmin, qmax, iXj,  1, 1. )
      if (associated(CO_phot))    call pmaxmin('CO: phot',       CO_phot, qmin, qmax, iXj,  1, 1. )
      if (associated(CO_column))  call pmaxmin('CO: column',   CO_column, qmin, qmax, iXj,  1, 1. )
      if (associated(CO_surface)) call pmaxmin('CO: surface', CO_surface, qmin, qmax, iXj,  1, 1. )
      if (associated(CO_dry))     call pmaxmin('CO: dry',         CO_dry, qmin, qmax, iXj, km, 1. )
   endif

!  Housekeeping
!  ------------
   deallocate(ndwet, p, pe, rkoh, rkch4_oh, rkch4_cl, rkch4_o1d, STAT=ier(1))
   deallocate(photJ, dCOPhot, aj, STAT=ier(1))

   return

contains
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

   REAL, ALLOCATABLE :: pblLayer(:,:), sfcFlux(:,:), fPBL(:,:,:)

   rc = 0

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
                                          (MAPL_AIRMW/mwtCO)/(w_c%delp(i1:i2,j1:j2,k)/MAPL_GRAV)
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
                                         (MAPL_AIRMW/mwtCO)/(w_c%delp(i1:i2,j1:j2,km)/MAPL_GRAV)
! Fossil fuel and biofuel
! -----------------------
   sfcFlux(i1:i2,j1:j2) = gcCO%eCO_fosfuel(i1:i2,j1:j2)+gcCO%eCO_biofuel(i1:i2,j1:j2)
   gcCO%COsfcFlux(i1:i2,j1:j2) = gcCO%COsfcFlux(i1:i2,j1:j2)+sfcFlux(i1:i2,j1:j2)
   w_c%qa(nbeg)%data3d(i1:i2,j1:j2,km) = w_c%qa(nbeg)%data3d(i1:i2,j1:j2,km)+sfcFlux(i1:i2,j1:j2)*cdt* &
                                        (MAPL_AIRMW/mwtCO)/(w_c%delp(i1:i2,j1:j2,km)/MAPL_GRAV)
! Release memory
! --------------
   DEALLOCATE(sfcFlux,STAT=ios)

   RETURN
   END SUBROUTINE CO_Emission

!-------------------------------------------------------------------------
!  Borrowed from meso_phot.F of StratChem, where number densities are cgs [cm^{-3}]
   SUBROUTINE getJRates(rc)

   IMPLICIT NONE

   INTEGER, INTENT(out) :: rc

   REAL, ALLOCATABLE :: o2Column(:,:,:)
   REAL, ALLOCATABLE :: SZARad(:,:)
   REAL, ALLOCATABLE :: SZADeg(:,:)
   REAL, ALLOCATABLE :: sinSZA(:,:)
   REAL, ALLOCATABLE :: zgrz(:,:)
   REAL, ALLOCATABLE :: sfaca(:,:)
   REAL, ALLOCATABLE :: arg(:,:)

   REAL, PARAMETER :: wavel = 1215.7
   REAL, PARAMETER :: O2xs  = 1.000E-20
   REAL, PARAMETER :: CH4xs = 2.000E-17
   REAL, PARAMETER :: sflux = 4.006E+11

!  Constants for Chapman function at high solar zenith angle
!  ---------------------------------------------------------
   REAL, PARAMETER :: hbar = 6.79
   REAL, PARAMETER :: zbar = 30.0
   REAL, PARAMETER :: r0   = 6.371E+03
   REAL, PARAMETER :: zp   = 60.0

   REAL, PARAMETER :: d1 = 1.060693
   REAL, PARAMETER :: d2 = 0.55643831
   REAL, PARAMETER :: d3 = 1.0619896
   REAL, PARAMETER :: d4 = 1.7245609
   REAL, PARAMETER :: d5 = 0.56498823
   REAL, PARAMETER :: d6 = 0.06651874

   REAL, PARAMETER :: O2ABV80KM = 7.072926E+19 ![cm^{-2}]
   REAL, PARAMETER :: O2VMR = 0.20946

!  bweir: Hardcoding because I don't want to brick old RC files (FIXME)
   REAL, PARAMETER :: SZACUTOFF = 70.0

   REAL :: b, r, s
   REAL :: fO2VMR

   INTEGER :: status

   CHARACTER(LEN=ESMF_MAXSTR) :: Iam = "CO::getJRates"

   rc = 0
   b = SQRT(0.50*r0/hbar)

!  O2 overhead number density profile [cm^{-2}]
!  --------------------------------------------
   ALLOCATE(O2Column(i1:i2,j1:j2,1:km), STAT=status)
   VERIFY_(status)

   fO2VMR = O2VMR*5.00E-05
!  O2Column(:,:,1) = O2ABV80KM+cellDepth(:,:,1)*ndwet(:,:,1)*fO2VMR
   O2Column(i1:i2,j1:j2,1) = O2ABV80KM + w_c%delp(i1:i2,j1:j2,1)*MAPL_AVOGAD/(MAPL_AIRMW*MAPL_GRAV)*fO2VMR

   DO k = 2,km
!     O2Column(:,:,k) = O2Column(:,:,k-1)+(cellDepth(:,:,k-1)*ndwet(:,:,k-1)+ &
!                                          cellDepth(:,:,  k)*ndwet(:,:,  k))*fO2VMR
      O2Column(i1:i2,j1:j2,k) = O2Column(i1:i2,j1:j2,k-1) &
                              + (w_c%delp(i1:i2,j1:j2,k-1) + w_c%delp(i1:i2,j1:j2,k)) &
                                *MAPL_AVOGAD/(MAPL_AIRMW*MAPL_GRAV)*fO2VMR
   END DO

   IF(gcCO%dbg) THEN
      CALL pmaxmin('CO: O2Column', O2Column, qmin, qmax, iXj, km,  1. )
   END IF

!  Grab some memory
!  ----------------
   ALLOCATE(SZARad(i1:i2,j1:j2), STAT=status)
   VERIFY_(status)
   ALLOCATE(SZADeg(i1:i2,j1:j2), STAT=status)
   VERIFY_(status)
   ALLOCATE(sinSZA(i1:i2,j1:j2), STAT=status)
   VERIFY_(status)

   WHERE(w_c%cosz(i1:i2,j1:j2) > 1.00)
      SZARad(i1:i2,j1:j2) = 0.00
   ELSEWHERE
      SZARad(i1:i2,j1:j2) = ACOS(w_c%cosz(i1:i2,j1:j2))
   ENDWHERE
   SZADeg(i1:i2,j1:j2) = SZARad(i1:i2,j1:j2)*radToDeg
   sinSZA(i1:i2,j1:j2) = SIN(SZARad(i1:i2,j1:j2))

   ALLOCATE(zgrz(i1:i2,j1:j2), STAT=status)
   VERIFY_(status)

   WHERE(SZADeg(i1:i2,j1:j2) <= 90.00)
      zgrz(i1:i2,j1:j2) = 1000.00
   ELSEWHERE
      zgrz(i1:i2,j1:j2) = sinSZA(i1:i2,j1:j2)*(zp+r0)-r0
   ENDWHERE

   IF(gcCO%dbg) THEN
      CALL pmaxmin('CO: zgrz',     zgrz, qmin, qmax, iXj, 1,  1. )
      CALL pmaxmin('CO: cosz', w_c%cosz, qmin, qmax, iXj, 1,  1. )
   END IF

   ALLOCATE(sfaca(i1:i2,j1:j2), STAT=status)
   VERIFY_(status)
   sfaca(i1:i2,j1:j2) = 0.00

!  Chapman function calculation from ACDB 2-D model
!  ------------------------------------------------
   DO j = j1,j2
      DO i = i1,i2
         Daytime: IF(SZADeg(i,j) < SZACUTOFF) THEN
            IF(SZADeg(i,j) < 70.00) THEN
               sfaca(i,j) = 1.00/w_c%cosz(i,j)
            ELSE IF(zgrz(i,j) > 0.00) THEN
               s = b*ABS(w_c%cosz(i,j))
               IF(s <= 8.00) THEN
                  s = (d1+d2*s)/(d3+d4*s+s**2)
               ELSE
                  s = d5/(d6+s)
               ENDIF

               r = b*SQRT(MAPL_PI)
               sfaca(i,j) = r*s

               IF(SZADeg(i,j) > 90.00) THEN
                  sfaca(i,j) = 2.00*r*EXP((r0+zbar)*(1.00-sinSZA(i,j))/hbar)-sfaca(i,j)
               ENDIF
            ENDIF
         ENDIF Daytime
      ENDDO
   ENDDO

   IF(gcCO%dbg) THEN
      CALL pmaxmin('CO: sfaca', sfaca, qmin, qmax, iXj, 1,  1. )
   END IF

   ALLOCATE(arg(i1:i2,j1:j2), STAT=status)
   VERIFY_(status)

!  At each layer, compute the rate constant, J [s^{-1}], if the sun is up
!  ----------------------------------------------------------------------
   DO k = 1,km
      WHERE(SZADeg(i1:i2,j1:j2) < SZACUTOFF)
         arg(i1:i2,j1:j2) = O2Column(i1:i2,j1:j2,k)*O2xs*sfaca(i1:i2,j1:j2)
         photJ(i1:i2,j1:j2,k) = sflux*EXP(-arg(i1:i2,j1:j2))*CH4xs
      ENDWHERE
   ENDDO

   IF(gcCO%dbg) THEN
      CALL pmaxmin('CO: photJ', photJ, qmin, qmax, iXj, km,  1. )
   ENDIF

   DEALLOCATE(SZARad, STAT=status)
   VERIFY_(status)
   DEALLOCATE(SZADeg, STAT=status)
   VERIFY_(status)
   DEALLOCATE(sinSZA, STAT=status)
   VERIFY_(status)
   DEALLOCATE(zgrz, STAT=status)
   VERIFY_(status)
   DEALLOCATE(sfaca, STAT=status)
   VERIFY_(status)
   DEALLOCATE(arg, STAT=status)
   VERIFY_(status)
   DEALLOCATE(O2Column, STAT=status)
   VERIFY_(status)

   RETURN
   END SUBROUTINE getJRates

   SUBROUTINE interp_s(k,sza,o3column,s,jo2,gcCO)
! ----------------------------------------------------------------------------
! NAME:
!   interp_s
!
! PURPOSE:
!   Interpolate S values for each wavelength in table to specified O3
!   column and zenith angle
!
! INPUTS:
!   k         Current layer number
!   szaRad    Solar zenith angle [radians]
!   o3column  Overhead o3 column value [cm^{-2}]
!   gcCO      The GOCART::CO grid component, which contains
!     sza_tab Solar zenith angle table
!     o3_tab  Overhead O3 values table
!     sdat    Radiative source function 
!     o2jdat  Table of J(O2) values
!
! OUTPUTS:
!   s         S value for each wavelength at current k, interpolated to
!               the given o3column and sza
!   jo2       J(O2) values interpolated as above
!
! 
! PROCEDURE:
!   Bi-linear interpolation, for sza > 94 s=0, for O3 out of range use min/max
!
! MODIFICATION HISTORY: 
!   25 Aug 1993  Kawa
!   10 Jul 1996  Kawa    For 28 levels and to handle J(O2) separately
!   11 May 2012  Nielsen Accomodation for GEOS-5 FV cubed release
!   30 Jan 2021  Weir    Copied from StratChem
! ----------------------------------------------------------------------------
       
   IMPLICIT NONE

   TYPE(CO_GridComp1), INTENT(IN) :: gcCO   ! Grid Component

   INTEGER, INTENT(IN) :: k
   REAL, INTENT(IN) :: sza, o3column 
   REAL, INTENT(OUT) :: s(gcCO%nlam), jo2

   INTEGER :: ijj, ik, ikk, ikkm, il, is
   REAL :: omt, omu, t, u
   REAL, PARAMETER :: PI = 3.14159265

! For each input solar zenith angle, find the first element of gcCO%sza_tab that 
! is greater.  Use this element and previous one to determine the interpolated value.
! -----------------------------------------------------------------------------------
   DO is = 1,gcCO%nsza
      ijj = is 
      IF(gcCO%sza_tab(is) > sza) EXIT 
   ENDDO
      
! Zenith angle test       
! -----------------
   IF(sza > gcCO%sza_tab(gcCO%nsza)) THEN
!     Cell is dark, set s and jo2=0        
!     -----------------------------
      s(1:gcCO%nlam) = 0.
      jo2 = 0.
   ELSE  
!     Cell is illuminated     
!     -------------------
      t = (sza-gcCO%sza_tab(ijj-1))/(gcCO%sza_tab(ijj)-gcCO%sza_tab(ijj-1))
      omt = 1.-t
         
! For each overhead O3 column, find the first element in gcCO%o3_tab that is
! greater. Use this element and previous one to determine the interpolated value.
! -------------------------------------------------------------------------------
      DO is = 1,gcCO%numo3
         ikk = is 
         IF(gcCO%o3_tab(is,k) > o3column) EXIT
      ENDDO

      ikkm = ikk-1 
      IF(ikk > 1 .AND. o3column <= gcCO%o3_tab(gcCO%numo3,k)) THEN
         u = (o3column-gcCO%o3_tab(ikkm,k))/(gcCO%o3_tab(ikk,k)-gcCO%o3_tab(ikkm,k))
         omu = 1.-u

! Do bilinear interpolation for each wavelength.
! ----------------------------------------------
         DO il = 1,gcCO%nlam       
            s(il) = omt*omu*gcCO%sdat(ijj-1,ikkm,k,il)+t*omu*gcCO%sdat(ijj,ikkm,k,il)+ &
                    t*u*gcCO%sdat(ijj,ikk,k,il)+omt*u*gcCO%sdat(ijj-1,ikk,k,il)
         ENDDO
         jo2 = omt*omu*gcCO%o2jdat(ijj-1,ikkm,k)+t*omu*gcCO%o2jdat(ijj,ikkm,k)+ &
               t*u*gcCO%o2jdat(ijj,ikk,k)+omt*u*gcCO%o2jdat(ijj-1,ikk,k)
    
! Extrapolate ahead of table
! --------------------------
      ELSE IF (ikk == 1) THEN
         DO il = 1,gcCO%nlam
            s(il) = omt*gcCO%sdat(ijj-1,1,k,il)+t*gcCO%sdat(ijj,1,k,il)
         ENDDO
         jo2 = omt*gcCO%o2jdat(ijj-1,1,k)+t*gcCO%o2jdat(ijj,1,k)

! Extrapolate beyond table
! ------------------------
      ELSE
         DO il = 1,gcCO%nlam
            s(il) = omt*gcCO%sdat(ijj-1,gcCO%numo3,k,il)+t*gcCO%sdat(ijj,gcCO%numo3,k,il)
         END DO 
         jo2 = omt*gcCO%o2jdat(ijj-1,gcCO%numo3,k)+t*gcCO%o2jdat(ijj,gcCO%numo3,k)
      ENDIF  
   ENDIF
      
   RETURN
   END SUBROUTINE interp_s

   SUBROUTINE jcalc4(k,szan,o3column,press,kel,aj,gcCO)
! ---------------------------------------------------------------------------------
! NAME: jcalc4
! PURPOSE:
!   Calculate photolysis rates
! INPUT:
!   k         Current layer number
!   levels    Number of layers
!   szan      Solar zenith angle (radians)
!   o3column  Overhead O3 values
!   press     Mid-layer pressure (hPa)
!   kel       Mid-layer temperature (K)
! OUTPUT:
!   aj        Array of photolysis rates
! RESTRICTIONS:
!   Currently set up for 23-J set (see var gcCO%nxdo)
! REQUIRED ROUTINES:
!   interp_s
! MODIFICATION HISTORY: 
!   26 Aug 1993 Kawa    Created
!   23 Nov 1993 Kawa    Remade xtab to do multiplication by solar flux beforehand 
!                        and removed inputs.
!   25 Feb 1994         Add 3 additional Js, incl N2O
!   18 Sep 1995         Add 2 additional Js, up to 22, and do CH2O special
!   13 May 1996 Crum    Removed fossils, move toward Fortran 90
!   10 Jul 1996         Modified to handle J(O2) separately and use 28 levels
!    1 Apr 2009 Nielsen GEOS-5 form with standardized SC_GridComp interface.
!    1 Jun 2009 Nielsen Updated to JPL 2006
!   12 Dec 2010 Nielsen Updated to JPL 2010 following Luke Oman's testing.
!   11 May 2012 Nielsen Accomodation for GEOS-5 FV cubed release
!    3 Jun 2015 Liang   Updated to the new 50-slot table with addition of halons,
!                       HCFCs, and 5 VSLSs
!                       numphoto is now updated to 52
!   30 Jan 2021 Weir    Copied from StratChem
!
! WARNING: Photolysis reaction rate numbers 38-42 are calculated in MESO_PHOT.
! ---------------------------------------------------------------------------------
   IMPLICIT NONE
   INTEGER, PARAMETER :: DBL = KIND(0.00D+00)

   TYPE(CO_GridComp1), INTENT(INOUT) :: gcCO   ! Grid Component

   INTEGER, INTENT(IN) :: k
   REAL, INTENT(IN) :: szan, o3column, press, kel
!  REAL(KIND=DBL), INTENT(OUT) :: aj(gcCO%numphoto)
!  bweir: demoted to single
   REAL, INTENT(OUT) :: aj(gcCO%numphoto)

   INTEGER :: ilam,indt,ix

   REAL :: alpha300, alphat, jo2, rjs(gcCO%nxdo), q1, q2, r1mq1
   REAL :: s(gcCO%nlam), sx(2,gcCO%nlam), tfac, wvl

! Start with a clean slate
! ------------------------
   aj(1:gcCO%numphoto) = 0.

! Interpolate radiative flux function values to model conditions
! --------------------------------------------------------------
   CALL interp_s(k,szan,o3column,s,jo2,gcCO)
   indt = kel-148.5
   indt = MAX(1,indt)
   indt = MIN(indt,200)

! Preliminaries for CH2O quantum yield dependence on m, T, wavelength
! -------------------------------------------------------------------
   tfac = (kel-80.0)/80.0

   DO ilam=1,gcCO%nlam
      ZeroS: IF(s(ilam) == 0.) THEN
         sx(1,ilam) = 0.00
         sx(2,ilam) = 0.00
      ELSE 

         wvl = gcCO%rlam(ilam)*0.10

         IF(wvl < 250.00) THEN
            q1 = 0.24
         ELSE IF(wvl >= 339.00) THEN
            q1 = 0.00
         ELSE
            q1 = gcCO%CH2O_aq(1) + gcCO%CH2O_aq(2)*wvl         + &
                                   gcCO%CH2O_aq(3)*wvl*wvl     + &
                                   gcCO%CH2O_aq(4)*wvl*wvl*wvl + &
                                   gcCO%CH2O_aq(5)*wvl*wvl*wvl*wvl
         ENDIF

         r1mq1 = 1./(1.-q1)

         IF(wvl < 330.00) THEN
            q2 = gcCO%xtab(ilam,22,indt)
         ELSE IF(wvl > 360.00) THEN
            q2 = 0.00
         ELSE
            alpha300 = 1.00E-03*(1./gcCO%xtab(ilam,22,1)-r1mq1)
            alphat = alpha300*(1.+0.05*(wvl-329.)*((300.-kel)/80.))
            q2 = 1.00/(r1mq1+alphat*press)
         ENDIF

         IF(wvl .LT. 250.00) q2=0.5

         sx(2,ilam) = s(ilam)*gcCO%xtab(ilam,21,indt)*q2
         sx(1,ilam) = s(ilam)*gcCO%xtab(ilam,21,indt)*q1
      ENDIF ZeroS
   ENDDO

! J(BrONO2) through J(OCLO)
! -------------------------
   DO ix=1,14
      rjs(ix) = 0.

      DO ilam=1,gcCO%nlam
         rjs(ix) = rjs(ix)+s(ilam)*gcCO%xtab(ilam,ix,indt)
      ENDDO
   ENDDO

! J(O2)
! -----
   rjs(15) = jo2

! J(O3_O1D) through J(N2O)
! ------------------------
   DO ix=16,20
      rjs(ix) = 0.

      DO ilam=1,gcCO%nlam
         rjs(ix) = rjs(ix)+s(ilam)*gcCO%xtab(ilam,ix,indt)
      ENDDO
   ENDDO

! J(CH2O)
! -------
   rjs(21) = 0.
   rjs(22) = 0.
   DO ilam=1,gcCO%nlam
      rjs(21) = rjs(21)+sx(1,ilam)
      rjs(22) = rjs(22)+sx(2,ilam)
   ENDDO

! J(CO2 -> CO + O) through xH1211
! -------------------------------
   DO ix=23,gcCO%nxdo
      rjs(ix) = 0.

      DO ilam=1,gcCO%nlam
         rjs(ix) = rjs(ix)+s(ilam)*gcCO%xtab(ilam,ix,indt)
      ENDDO
   ENDDO
               
! ---------------------------------------------------------------
! Order photolysis rates to match order in full chemistry model.  
! Sort rjs into CTM photolysis rate array, aj.  Order of rjs:
!
!  1-J(BrONO2)
!  2-J(BrO)
!  3-J(Cl2O2)
!  4-J(ClONO2)
!  5-J(H2O2)
!  6-J(HCl)
!  7-J(HNO3)
!  8-J(HO2NO2)
!  9-J(HOCl)
! 10-J(N2O5)
! 11-J(NO2)
! 12-J(NO3_NO)
! 13-J(NO3_NO2)
! 14-J(OClO)
! 15-J(O2)
! 16-J(O3_O1D)
! 17-J(O3_3P)
! 18-J(HOBr)
! 19-J(CH3OOH)
! 20-J(N2O)
! 21-J(CH2O_HCO)
! 22-J(CH2O_CO)
! 23-J(CO2 -> CO + O)
! 24-xCFC-11
! 25-xCFC-12
! 26-xCCl4
! 27-xCH3CCl3
! 28-xHCFC-22
! 29-xCFC-113
! 30-xCH3Cl
! 31-xCH3Br
! 32-xH1301
! 33-xH1211 
! 34-xH1202
! 35-xH2402
! 36-xCHBr3
! 37-xCH2Br2
! 38-xCH2ClBr
! 39-xCHClBr2
! 40-xCHCl2Br
! 41-xHCFC-141b
! 42-xHCFC-142b
! 43-xCFC-114 
! 44-xCFC-115
! 45-xOCS
! 46-
! 47-
! 48-
! 49-
! 50-
! ---------------------------------------------------------------
! Solar cycle goes here when ready  
!     aj( 1) = rjs(15)*gcCO%s_cycle(3,gcCO%iscyr)
! ----------------------------------------------------------------
   aj( 1) = rjs(15)
   aj( 2) = rjs(16)
   aj( 3) = rjs(17)
! H2O
! ---
   aj( 4) = 0.
   aj( 5) = rjs(13)
   aj( 6) = rjs(7)
   aj( 7) = rjs(11)
   aj( 8) = rjs(5)
   aj( 9) = rjs(10)
   aj(10) = rjs(21)
   aj(11) = rjs(22)
   aj(12) = rjs(23)
   aj(13) = rjs(19)
   aj(14) = rjs(20)
   aj(15) = rjs(4)
   aj(16) = 0.
   aj(17) = rjs(12)
   aj(18) = rjs(6)
   aj(19) = 0.

! CH3Br(20) H1301(21) H12_24(22)
! ------------------------------
   aj(20) = rjs(31)
   aj(21) = rjs(32)
   aj(22) = rjs(33)
   aj(23) = rjs(9)
   aj(24) = rjs(8)
   aj(25) = rjs(18)
   aj(26) = 0.
   aj(27) = rjs(2)
   aj(28) = rjs(1)

! F11(29) F12(30) CCl4(31) CHCCl3(32) HCFC(33) F113(34) CH3Cl(35)
! ---------------------------------------------------------------
   aj(29) = rjs(24)
   aj(30) = rjs(25)
   aj(31) = rjs(26)
   aj(32) = rjs(27)
   aj(33) = rjs(28)
   aj(34) = rjs(29)
   aj(35) = rjs(30)
   aj(36) = rjs(3)
   aj(37) = rjs(14)

! ------------------------------------------
! WARNING: Photolysis reaction rate
! numbers 38-42 are calculated in MESO_PHOT.
! ------------------------------------------
! Add aj(43) which is J(Cl2O2) for partitioning but not Ox loss 
! which is aj(36). In lookup table J(Cl2O2) is J*qy where qy is 0.8 
! so multiply by 1.25 to equal J and used in part.F and partest.F

   aj(43) = rjs(3)*1.25

! QingLiang -- 06/03/2015
! CHBr3(44) CH2Br2(45) CH2BrCl(46) CHBrCl2(47) CHBr2Cl(48)
   aj(44) = rjs(36)
   aj(45) = rjs(37)
   aj(46) = rjs(38)
   aj(47) = rjs(39)
   aj(48) = rjs(40)

! QingLiang -- 06/03/2015
! Add two new halons: H-1202 (49) H2402 (50) 
! and two new HCFCs: HCFC-141b (51) HCFC-142b (52) 
   aj(49) = rjs(34)
   aj(50) = rjs(35)
   aj(51) = rjs(41)
   aj(52) = rjs(42)

! QingLiang -- 02/05/2016
! Add CFC-114 and CFC-115
! Add OCS for GOCART module
   aj(53) = rjs(43)
   aj(54) = rjs(44)
   aj(55) = rjs(45)
!  aj(53) = rjs(34)
!  aj(54) = rjs(34)
!  aj(55) = rjs(34)

   RETURN
   END SUBROUTINE jcalc4

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
                gcCO%eCO_mtn, gcCO%CH4, gcCO%OHnd, STAT=ios )
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
