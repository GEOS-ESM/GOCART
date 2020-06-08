#include "MAPL_Generic.h"

!-------------------------------------------------------------------------
!                       NASA/GSFC, GMAO, Code 610.1                      !
!-------------------------------------------------------------------------
!BOP
!
!
! !MODULE:  Bundle_IncrementMod --- Computes tracer increments and puts them into a bundle
! 
! !INTERFACE:
!

   MODULE  Bundle_IncrementMod

! !USES:



!
! !DESCRIPTION: This module can be used to compute model variable increments if the
!               following conditions are met:
!                  1) A non-increment, or "parent", bundle exists (e.g. TR, TRADV, MTR)
!                  2) The increment is computed as: (X_t2 - X_t1)/(t2-t1)
!                For new increment bundles the user must do at least the following:
!                  1) Define an increment bundle in Set Services of a gridded component.
!                  2) add an "IncType" and make appropriate changes to the code for the new "IncType".

! NOTE!!! - The MTRI (moist increments) are NOT mass weighted. If DYCORE uses mass weighting in 
!           its solution, the moist increment will be incorrect.


!EOP
!-------------------------------------------------------------------------

!BOC

   USE ESMF
   USE MAPL

   IMPLICIT NONE
   PRIVATE
 
! !PUBLIC MEMBER FUNCTIONS: 
   PUBLIC Initialize_IncBundle_init
   PUBLIC Initialize_IncBundle_run
   PUBLIC Compute_IncBundle

   type(ESMF_Field)                          :: field, TempField
   character(len=ESMF_MAXSTR)                :: fieldname
   character(len=ESMF_MAXSTR)                :: org_bundle       ! Original bundle name
   character(len=ESMF_MAXSTR)                :: inc_bundle       ! Increment bundle name
   integer                                   :: i, ppos, NQ

   integer, parameter, public                :: DYNinc      = 1
   integer, parameter, public                :: H2Oinc      = 2
   integer, parameter, public                :: CHMinc      = 3
   integer, parameter, public                :: MTRIinc     = 4
   integer, parameter, public                :: CHMincR2    = 5
   integer, parameter, public                :: MTRIincCTM  = 6
   integer, parameter, public                :: TRIincCTM   = 7


CONTAINS

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!BOP
!
!  !ROUTINE:  Initialize_IncBundle_init - Initialize increment bundle with fields from a 
!             "parent" bundle. This must be called in the Initialize method of a gridded component.

 SUBROUTINE Initialize_IncBundle_init(GC, state1, state2, IncType, RC)

  IMPLICIT NONE
  
  ! ARGUMENTS
  type(ESMF_GridComp),        intent(in)        :: GC          ! Gridded component 
  type(ESMF_State),           intent(in)        :: state1      ! Original bundle state
  type(ESMF_State),           intent(in)        :: state2      ! Increment bundle state
  integer,                    intent(in)        :: IncType     ! Increment bundle type
  integer, optional,          intent(  out)     :: RC          ! Error code


  ! TYPES and VARIABLES
  integer                                       :: STATUS
  type(ESMF_FieldBundle)                        :: BUNDLE, BUNDLEi, BUNDLEemiss
  type(ESMF_Config)                             :: cf               ! AGCM.rc
  character(len=ESMF_MAXSTR)                    :: IAm, valueOld
  character(len=ESMF_MAXSTR)                    :: incEmiss_bundle  ! Increment chemistry emissions bundle name
  character(len=ESMF_MAXSTR)                    :: longname         ! longname metadata description
  character(len=2)                              :: suffix           ! suffix appended to variable name
  integer                                       :: nCols
  character(len=ESMF_MAXSTR), allocatable       :: NAMES(:)

  Iam = "Initialize_IncBundle_init"

! ============================================================================

! Begin...
    if (IncType == DYNinc) then
      org_bundle = 'TRADV'
      inc_bundle = 'TRADVI'
      suffix = 'ID'
      longname = '_due_to_dynamics'
    else if (IncType == H2Oinc) then
      org_bundle = 'TRADV'
      inc_bundle = 'H2ORTRI'
      suffix = 'IW'
      longname = '_due_to_water_rescaling'
    else if (IncType == CHMinc) then
      org_bundle = 'CHEM_TRACERS'
      inc_bundle = 'CHEMTRI'
      incEmiss_bundle = 'CHEMTRIr1'
      suffix = 'IC'
      longname = '_due_to_chemistry'
    else if (IncType == MTRIinc) then
      org_bundle = 'MTR'
      inc_bundle = 'MTRI'
      suffix = 'IM'
      longname = '_due_to_moist_processes'
    else if (IncType == MTRIincCTM) then
      org_bundle = 'ConvTR'
      inc_bundle = 'MTRI'
      suffix = 'IM'
      longname = '_due_to_moist_processes'
    else if (IncType == TRIincCTM) then
      org_bundle = 'DiffTR'
      inc_bundle = 'TRI'
      suffix = 'IT'
      longname = '_due_to_turbulence_processes'
    end if

    call ESMF_GridCompGet ( GC, config=cf, RC=STATUS )
    VERIFY_(STATUS)

    call ESMF_ConfigGetDim (cf, NQ, nCols, label=(trim(inc_bundle)//'_increments::'), rc=STATUS)

    if (NQ > 0) then
      call ESMF_ConfigFindLabel (cf, (trim(inc_bundle)//'_increments::'), rc=STATUS)
      VERIFY_(STATUS)

      allocate (NAMES(NQ), stat=STATUS)
      VERIFY_(STATUS)

      do i = 1, NQ
        call ESMF_ConfigNextLine(cf, rc=STATUS)
        VERIFY_(STATUS)
        call ESMF_ConfigGetAttribute(cf, NAMES(i), rc=STATUS)
        VERIFY_(STATUS)
      enddo

! Fill the increments bundle with fields from "parent" bundle 
!------------------------------------------------------------------
      call ESMF_StateGet(state1, org_bundle, BUNDLE,   rc=STATUS)
      VERIFY_(STATUS)

      call ESMF_StateGet(state2, inc_bundle, BUNDLEi, rc=STATUS)
      VERIFY_(STATUS)
  
      if (IncType == CHMinc) then
        call ESMF_StateGet(state2, incEmiss_bundle, BUNDLEemiss, rc=STATUS)
        VERIFY_(STATUS)
      end if

      do i = 1, NQ
        call ESMF_FieldBundleGet (BUNDLE, NAMES(i), field=field, rc=status)
        if (status/=0) then
          if (mapl_am_i_root()) print*, trim(NAMES(i)),' is not valid. It likely does not exist in ',trim(org_bundle)
          VERIFY_(23) 
        end if

        call ESMF_FieldGet (field, name=fieldname, RC=STATUS)
        VERIFY_(STATUS)

        TempField = MAPL_FieldCreate (field, name=(trim(fieldname)//suffix) ,DoCopy=.true., rc=status)
        VERIFY_(STATUS)
        call MAPL_FieldBundleAdd (BUNDLEi, TempField, rc=status)
        VERIFY_(STATUS)
     
        if (IncType == CHMinc) then
          TempField = MAPL_FieldCreate (field, name=(trim(fieldname)//suffix//'emiss') ,DoCopy=.true., rc=status)
          VERIFY_(STATUS)
          call MAPL_FieldBundleAdd (BUNDLEemiss, TempField, rc=status)
          VERIFY_(STATUS)
        end if
      end do

      deallocate (NAMES, stat=STATUS)

! Set Species Attributes
!----------------------------------
      if ((IncType==DYNinc) .OR. (IncType==CHMinc) .OR. (IncType==H2Oinc).OR. (IncType == MTRIinc) &
        .OR. (IncType == MTRIincCTM) .OR. (IncType == TRIincCTM)) then
        do i = 1, NQ
          call ESMF_FieldBundleGet (BUNDLEi, fieldIndex=i, field=field, rc=STATUS )
          VERIFY_(STATUS)
          call ESMF_FieldGet (field, name=fieldname, RC=STATUS)
          VERIFY_(STATUS)
 
          if (fieldname==('AOADAYS'//suffix)) then
            call ESMF_AttributeSET (field, name='UNITS', value='days s-1', rc=status)
            VERIFY_(STATUS)
          else
            call ESMF_AttributeGET (field, name='UNITS', value=valueOld, rc=status)
            VERIFY_(STATUS)
            call ESMF_AttributeSET (field, name='UNITS', value=trim(valueOld)//' s-1', rc=status)
            VERIFY_(STATUS)
          end if

          ppos = len(trim(fieldname))
          call ESMF_AttributeSET (field, name='LONG_NAME', value=('tendency_of_'//fieldname(1:ppos-2)//trim(longname)), rc=status)
          VERIFY_(STATUS)
        end do
      end if

      if (IncType == CHMinc) then
        do i = 1, NQ
          call ESMF_FieldBundleGet (BUNDLEemiss, fieldIndex=i, field=field, rc=STATUS )
          VERIFY_(STATUS)
          call ESMF_FieldGet (field, name=fieldname, RC=STATUS)
          VERIFY_(STATUS)

          if (fieldname==('AOADAYS'//suffix//'emiss')) then
            call ESMF_AttributeSET (field, name='UNITS', value='days s-1', rc=status)
            VERIFY_(STATUS)
          else
            call ESMF_AttributeGET (field, name='UNITS', value=valueOld, rc=status)
            VERIFY_(STATUS)
            call ESMF_AttributeSET (field, name='UNITS', value=trim(valueOld)//' s-1', rc=status)
            VERIFY_(STATUS)
          end if

          ppos = len(trim(fieldname))
          call ESMF_AttributeSET (field, name='LONG_NAME', value=('chemistry_tendency_of_'//fieldname(1:ppos-7)//'_from_emissions'), rc=status)
          VERIFY_(STATUS)
        end do
      end if ! (IncType == 'CHMinc')
    end if ! NQ > 0

    RETURN_(ESMF_SUCCESS)
 END SUBROUTINE Initialize_IncBundle_init
!==============================================================================================================

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 610.1 GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
!  !ROUTINE:  Initialize_IncBundle_run - Re-initialize increment bundle with data from the 
!  !          "parent" bundle within the Run method. The "parent" bundle is the non-increment 
!              bundle (e.g. TR, MTR, TRADV)

 SUBROUTINE Initialize_IncBundle_run(state1, state2, IncType, RC)

  IMPLICIT NONE

  ! ARGUMENTS
  type(ESMF_State),            intent(in)            :: state1           ! Original bundle state
  type(ESMF_State),            intent(in)            :: state2           ! Increment bundle state
  integer,                     intent(in)            :: IncType          ! Increment bundle type
  integer, optional,           intent(  out)         :: RC               ! Error code

  ! TYPES and VARIABLES
  integer                                            :: STATUS
  type(ESMF_FieldBundle)                             :: BUNDLE, BUNDLEi
  character(len=ESMF_MAXSTR)                         :: IAm
  character(len=ESMF_MAXSTR), allocatable            :: NAMES(:)
  real                                               :: DT
  real, dimension(:,:,:), pointer                    :: org_ptr, inc_ptr

  Iam = "Initialize_IncBundle_run"

! ============================================================================

! Begin...
    if (IncType == DYNinc) then
      org_bundle = 'TRADV'
      inc_bundle = 'TRADVI'
    else if (IncType == H2Oinc) then
      org_bundle = 'TRADV'
      inc_bundle = 'H2ORTRI'
    else if (IncType == CHMinc) then
      org_bundle = 'CHEM_TRACERS'
      inc_bundle = 'CHEMTRIr1'
    else if (IncType == CHMincR2) then
      org_bundle = 'CHEM_TRACERS'
      inc_bundle = 'CHEMTRI'
    else if (IncType == MTRIinc) then
      org_bundle = 'MTR'
      inc_bundle = 'MTRI'
    else if (IncType == MTRIincCTM) then
      org_bundle = 'ConvTR'
      inc_bundle = 'MTRI'
    else if (IncType == TRIincCTM) then
      org_bundle = 'DiffTR'
      inc_bundle = 'TRI'
    end if


!  !Initialize increment bundle in Run method before the child is called
!  !--------------------------------------------------------------------
    call ESMF_StateGet (state2, inc_bundle, BUNDLEi, rc=STATUS)
    VERIFY_(STATUS)
    call ESMF_StateGet (state1, org_bundle, BUNDLE, rc=STATUS)
    VERIFY_(STATUS)
    call ESMF_FieldBundleGet (BUNDLEi, fieldCount=NQ, rc=STATUS )
    VERIFY_(STATUS)

!  !Check if there is anything in the bundle.
    if (NQ > 0) then
      allocate (NAMES(NQ), stat=STATUS)
      VERIFY_(STATUS)

      call ESMF_FieldBundleGet(BUNDLEi, fieldNameList=NAMES, rc=STATUS)
      VERIFY_(STATUS)

!    !Get increment data pointer and initialize value
      do i = 1, NQ
        ppos = len(trim(NAMES(i)))
        if (IncType == CHMinc) then
          call ESMFL_BundleGetPointerToData (BUNDLE, trim(NAMES(i)(1:ppos-7)), org_ptr, rc=status)
          VERIFY_(STATUS)
        else
          call ESMFL_BundleGetPointerToData (BUNDLE, trim(NAMES(i)(1:ppos-2)), org_ptr, rc=status)
          VERIFY_(STATUS)
        end if

        call ESMFL_BundleGetPointerToData (BUNDLEi, trim(NAMES(i)), inc_ptr, rc=status)
        VERIFY_(STATUS)

        inc_ptr = org_ptr
      end do
      deallocate(NAMES)
    end if ! NQ > 0

    RETURN_(ESMF_SUCCESS)
 END SUBROUTINE Initialize_IncBundle_run
!=======================================================================================
!
!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 610.1 GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
!  !ROUTINE:  Compute_IncBundle - Compute the increment after the child has run.

 SUBROUTINE Compute_IncBundle(state1, state2, IncType, META, RC)

  IMPLICIT NONE

  ! ARGUMENTS
  type(ESMF_State),               intent(in)         :: state1           ! Original bundle state
  type(ESMF_State),               intent(in)         :: state2           ! Increment bundle state
  integer,                        intent(in)         :: IncType          ! Increment bundle type
  type(MAPL_MetaComp), pointer,   intent(in)         :: META
  integer, optional,              intent(  out)      :: RC               ! Error code


  ! TYPES and VARIABLES
  integer                                            :: STATUS
  type(ESMF_FieldBundle)                             :: BUNDLE, BUNDLEi, BUNDLEemiss
  character(len=ESMF_MAXSTR)                         :: IAm
  character(len=ESMF_MAXSTR)                         :: inc_emiss_bundle  ! Increment chemistry emissions bundle name
  character(len=ESMF_MAXSTR), allocatable            :: NAMES(:)
  real                                               :: DT
  real, dimension(:,:,:), pointer                    :: org_ptr, inc_ptr, inc_emiss_ptr

  Iam = "Compute_IncBundle"
! ============================================================================
! Begin...
    if (IncType == DYNinc) then
      org_bundle = 'TRADV'
      inc_bundle = 'TRADVI'
    else if (IncType == H2Oinc) then
      org_bundle = 'TRADV'
      inc_bundle = 'H2ORTRI'
    else if (IncType == CHMinc) then
      org_bundle = 'CHEM_TRACERS'
      inc_bundle = 'CHEMTRIr1'
    else if (IncType == CHMincR2) then
      org_bundle = 'CHEM_TRACERS'
      inc_bundle = 'CHEMTRI'
      inc_emiss_bundle = 'CHEMTRIr1'
    else if (IncType == MTRIinc) then
      org_bundle = 'MTR'
      inc_bundle = 'MTRI'
    else if (IncType == MTRIincCTM) then
      org_bundle = 'ConvTR'
      inc_bundle = 'MTRI'
    else if (IncType == TRIincCTM) then
      org_bundle = 'DiffTR'
      inc_bundle = 'TRI'
    end if

    call ESMF_StateGet (state2, inc_bundle, BUNDLEi, rc=STATUS)
    VERIFY_(STATUS)
    call ESMF_FieldBundleGet (BUNDLEi, fieldCount=NQ, rc=STATUS )
    VERIFY_(STATUS)

!  !Check if there is anything in the bundle.
    if (NQ > 0) then
      call ESMF_StateGet (state1, org_bundle, BUNDLE, rc=STATUS)
      VERIFY_(STATUS)

      if (IncType == CHMincR2) then
        call ESMF_StateGet (state2, inc_emiss_bundle, BUNDLEemiss, rc=STATUS)
        VERIFY_(STATUS)
      end if

      call MAPL_GetResource(META, DT, label="RUN_DT:", RC=STATUS)
      VERIFY_(STATUS)
      allocate (NAMES(NQ), stat=STATUS)
      VERIFY_(STATUS)

      call ESMF_FieldBundleGet(BUNDLEi, fieldNameList=NAMES, rc=STATUS)
      VERIFY_(STATUS)

!    !Get pointers to data
      do i = 1, NQ
        ppos = len(trim(NAMES(i)))
        if (IncType == CHMinc) then
          call ESMFL_BundleGetPointerToData (BUNDLE, trim(NAMES(i)(1:ppos-7)), org_ptr, rc=status)
          VERIFY_(STATUS)
          call ESMFL_BundleGetPointerToData (BUNDLEi, trim(NAMES(i)), inc_ptr, rc=status)
          VERIFY_(STATUS)
!        end if
        else if (IncType == CHMincR2) then
          call ESMFL_BundleGetPointerToData (BUNDLEemiss, trim(NAMES(i))//'emiss', inc_emiss_ptr, rc=status)
          VERIFY_(STATUS)
          call ESMFL_BundleGetPointerToData (BUNDLE, trim(NAMES(i)(1:ppos-2)), org_ptr, rc=status)
          VERIFY_(STATUS)
          call ESMFL_BundleGetPointerToData (BUNDLEi, trim(NAMES(i)), inc_ptr, rc=status)
          VERIFY_(STATUS)
        else
          call ESMFL_BundleGetPointerToData (BUNDLE, trim(NAMES(i)(1:ppos-2)), org_ptr, rc=status)
          VERIFY_(STATUS)
          call ESMFL_BundleGetPointerToData (BUNDLEi, trim(NAMES(i)), inc_ptr, rc=status)
          VERIFY_(STATUS)
        end if

!      !Compute increment and update pointer
        if (IncType == CHMincR2) then
          inc_ptr = ((org_ptr-inc_ptr)/DT) + inc_emiss_ptr
        else
          inc_ptr = (org_ptr-inc_ptr)/DT
        end if
      end do
      deallocate(NAMES)
    end if ! NQ > 0
 
    RETURN_(ESMF_SUCCESS)
 END SUBROUTINE Compute_IncBundle


END MODULE Bundle_IncrementMod
