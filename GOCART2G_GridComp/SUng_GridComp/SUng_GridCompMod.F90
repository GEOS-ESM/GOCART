#include "MAPL_Generic.h"

!=============================================================================
!BOP

! !MODULE: SUng_GridCompMod - GOCART refactoring of the SU gridded component 

! !INTERFACE:

module SUng_GridCompMod

!  !USES:
   USE ESMF
   USE MAPL_Mod

   implicit none
   private

! !PUBLIC MEMBER FUNCTIONS:
   PUBLIC  SetServices
!   PUBLIC  Initialize


! !DESCRIPTION: This module implements GOCARTS' Sulfate (SU) Gridded Component.

! !REVISION HISTORY:
! 29Oct2019  E.Sherman  First attempt at refactoring.

!EOP
!===========================================================================

contains


!BOP

! !IROUTINE: SetServices 

! !INTERFACE:

  subroutine SetServices (GC, RC)

! !ARGUMENTS:

    type (ESMF_GridComp), intent(INOUT)   :: GC  ! gridded component
    integer,              intent(  OUT)   :: RC  ! return code


! !DESCRIPTION: This version uses the {\tt GEOS\_GenericSetServices}, which sets
!   the Initialize and Finalize services to generic versions. It also
!   allocates our instance of a generic state and puts it in the 
!   gridded component (GC). Here we only set the two-stage run method and
!   declare the data services.
! \newline

! !REVISION HISTORY: 
! 29oct2019   E.Sherman  First attempt at refactoring

!EOP

!****************************************************************************
!
! ErrLog Variables

    character (len=ESMF_MAXSTR)                 :: IAm
    integer                                     :: STATUS
    character (len=ESMF_MAXSTR)                 :: COMP_NAME


    type (ESMF_Config)                          :: cfg

    character (len=ESMF_MAXSTR)                 :: field_name

    integer                                     :: n, i, n_bins
    real                                        :: DEFVAL
    logical                                     :: data_driven = .true.

    !development testing variables - to be deleted
    real, dimension(:,:), pointer       :: ptr_test

!****************************************************************************
!   Begin...

!   Get my name and set-up traceback handle
!   ---------------------------------------
    Iam = 'SetServices'
    call ESMF_GridCompGet (GC, NAME=COMP_NAME, __RC__)
    Iam = trim(COMP_NAME) // Iam

if (mapl_am_i_root()) print*,'GOCARTng SUng SetServices BEGIN'           ! for testing - to be deleted
if (mapl_am_i_root()) print*,'GOCARTng SU COMP_NAME = ', trim(COMP_NAME) ! for testing - to be deleted


!   Load resource file and get number of bins 
!   -------------------------------------------
    cfg = ESMF_ConfigCreate (__RC__)
    call ESMF_ConfigLoadFile (cfg, 'SUng_GridComp_'//trim(COMP_NAME)//'.rc', RC=STATUS)
    if (STATUS /= 0) then
        if (mapl_am_i_root()) print*,'SUngGridComp_'//trim(COMP_NAME)//'.rc does not exist! loading SUng_GridComp_SU.rc instead'
        call ESMF_ConfigLoadFile (cfg, 'SUng_GridComp_SU.rc', __RC__)
    end if

    call ESMF_ConfigGetAttribute (cfg, n_bins, label='bins:', __RC__)

!   Set entry points
!   ------------------------
!    call MAPL_GridCompSetEntryPoint (GC, ESMF_METHOD_INITIALIZE,  Initialize, __RC__)
!    call MAPL_GridCompSetEntryPoint (GC, ESMF_METHOD_RUN, Run1, __RC__)

!   Is SU data driven?
!   ------------------
!    call data_driven_ (COMP_NAME, data_driven, __RC__)


!   INTERNAL STATE
!   ---------------
!   Default INTERNAL state values
!   -----------------------------
    DEFVAL = 0.0

!   Aerosol Tracers to be transported
!   ---------------------------------
    call MAPL_AddInternalSpec(GC,                         &
      SHORT_NAME = trim(COMP_NAME)//'DMS',                &
      LONG_NAME  = 'Dimethylsulphide',                    &
      UNITS      = 'kg kg-1',                             &
      RESTART    = MAPL_RestartOptional,                  &
      DEFAULT    = DEFVAL,                                &
      DIMS       = MAPL_DimsHorzVert,                     &
      VLOCATION  = MAPL_VLocationCenter, __RC__)

    call MAPL_AddInternalSpec(GC,                         &
      SHORT_NAME = trim(COMP_NAME)//'SO2',                &
      LONG_NAME  = 'Sulphur dioxide',                     &
      UNITS      = 'kg kg-1',                             &
      RESTART    = MAPL_RestartOptional,                  &
      DEFAULT    = DEFVAL,                                &
      DIMS       = MAPL_DimsHorzVert,                     &
      VLOCATION  = MAPL_VLocationCenter, __RC__)

    call MAPL_AddInternalSpec(GC,                         &
      SHORT_NAME = trim(COMP_NAME)//'SO4',                &
      LONG_NAME  = 'Sulphate aerosol',                    &
      UNITS      = 'kg kg-1',                             &
      RESTART    = MAPL_RestartOptional,                  &
      DEFAULT    = DEFVAL,                                &
      DIMS       = MAPL_DimsHorzVert,                     &
      VLOCATION  = MAPL_VLocationCenter, __RC__)

    call MAPL_AddInternalSpec(GC,                         &
      SHORT_NAME = trim(COMP_NAME)//'MSA',                &
      LONG_NAME  = 'Methanesulphonic acid',               &
      UNITS      = 'kg kg-1',                             &
      RESTART    = MAPL_RestartOptional,                  &
      DEFAULT    = DEFVAL,                                &
      DIMS       = MAPL_DimsHorzVert,                     &
      VLOCATION  = MAPL_VLocationCenter, __RC__)







!   EXPORT STATE
!   -------------
    if (.not. data_driven) then
#       include "SUng_ExportSpec___.h"
    end if

!   This state holds fields needed by radiation
!   ---------------------------------------------
    call MAPL_AddExportSpec (GC,                             &
       SHORT_NAME = trim(COMP_NAME)//'_AERO',                &
       LONG_NAME  = 'dust_aerosols_from_'//trim(COMP_NAME),  &
       UNITS      = 'kg kg-1',                               &
       DIMS       = MAPL_DimsHorzVert,                       &
       VLOCATION  = MAPL_VLocationCenter,                    &
       DATATYPE   = MAPL_StateItem, __RC__)

!   This state is needed by MOIST - It will contain aerosols
!   ----------------------------------------------------------
    call MAPL_AddExportSpec (GC,                                              &
       SHORT_NAME = trim(COMP_NAME)//'_AERO_ACI',                             &
       LONG_NAME  = 'aerosol_cloud_interaction_dust_from_'//trim(COMP_NAME),  &
       UNITS      = 'kg kg-1',                                                &
       DIMS       = MAPL_DimsHorzVert,                                        &
       VLOCATION  = MAPL_VLocationCenter,                                     &
       DATATYPE   = MAPL_StateItem, __RC__)

!   This bundle is needed by surface for snow albedo modification
!   by aerosol settling and deposition
!   DEVELOPMENT NOTE - Change to StateItem in future
!   ---------------------------------------------------------------
    call MAPL_AddExportSpec (GC,                                  &
       SHORT_NAME = trim(COMP_NAME)//'_AERO_DP',                  &
       LONG_NAME  = 'aerosol_deposition_from_'//trim(COMP_NAME),  &
       UNITS      = 'kg m-2 s-1',                                 &
       DIMS       = MAPL_DimsHorzOnly,                            &
       DATATYPE   = MAPL_BundleItem, __RC__)


if (mapl_am_i_root()) print*,'GOCARTng SUng SetServices END' ! for testing - to be deleted

!   Set generic services
!   ----------------------------------
    call MAPL_GenericSetServices (GC, __RC__)

    RETURN_(ESMF_SUCCESS)

  end subroutine SetServices



















end module SUng_GridCompMod

