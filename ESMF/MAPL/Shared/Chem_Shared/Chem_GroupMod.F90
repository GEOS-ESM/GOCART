#include "MAPL_Generic.h"

!! 
!! #define CHEM_INFO
!! #define DEBUG_STUFF
!! 


!---------------------------------------------------------------------------
!      NASA/GSFC, Atmospheric Chemistry and Dynamics Lab,  Code 614        !
!---------------------------------------------------------------------------
!BOP
!
! !MODULE:  Chem_GroupMod --- Support for Bry and Cly chemical families
!
! !INTERFACE:
!

   MODULE  Chem_GroupMod

! !USES:

   USE ESMF
   USE MAPL

   IMPLICIT NONE

! !PUBLIC TYPES:
!
   PRIVATE

!
! !PUBLIC MEMBER FUNCTIONS:
!

   PUBLIC  Init_GMI_Chem_Groups
   PUBLIC  Init_GCC_Chem_Groups

   PUBLIC  Pack_Chem_Groups
   PUBLIC  Unpack_Chem_Groups

   PRIVATE Set_GMI_Globals
   PRIVATE Set_GCC_Globals

   LOGICAL   :: GMI_groups_active = .FALSE.
   LOGICAL   :: GCC_groups_active = .FALSE.

   REAL, PARAMETER  :: mw_BR = 79.904
   REAL, PARAMETER  :: mw_CL = 35.45
   REAL, PARAMETER  :: mw_N  = 14.007
   REAL, PARAMETER  :: mw_O  = 15.999
   REAL, PARAMETER  :: mw_H  =  1.008

   REAL, PARAMETER  :: mw_BR2     = (2*mw_BR                 )
   REAL, PARAMETER  :: mw_BRCL    = (  mw_BR + mw_CL         )
   REAL, PARAMETER  :: mw_BRNO2   = (  mw_BR + mw_N  + 2*mw_O)
   REAL, PARAMETER  :: mw_BRNO3   = (  mw_BR + mw_N  + 3*mw_O)
   REAL, PARAMETER  :: mw_BRO     = (  mw_BR + mw_O          )
   REAL, PARAMETER  :: mw_HBR     = (  mw_BR + mw_H          )
   REAL, PARAMETER  :: mw_HOBR    = (  mw_BR + mw_H  +   mw_O)

   REAL, PARAMETER  :: mw_CL2     = (2*mw_CL                  )
   REAL, PARAMETER  :: mw_CL2O2   = (2*mw_CL + 2*mw_O         )
   REAL, PARAMETER  :: mw_CLNO2   = (  mw_CL +   mw_N + 2*mw_O)
   REAL, PARAMETER  :: mw_CLNO3   = (  mw_CL +   mw_N + 3*mw_O)
   REAL, PARAMETER  :: mw_CLO     = (  mw_CL +   mw_O         )
   REAL, PARAMETER  :: mw_CLOO    = (  mw_CL + 2*mw_O         )
   REAL, PARAMETER  :: mw_HCL     = (  mw_CL +   mw_H         )
   REAL, PARAMETER  :: mw_HOCL    = (  mw_CL +   mw_H +   mw_O)
   REAL, PARAMETER  :: mw_OCLO    = (  mw_CL + 2*mw_O         )

   REAL, PARAMETER  :: mw_NO      = (  mw_N  +   mw_O         )
   REAL, PARAMETER  :: mw_NO2     = (  mw_N  + 2*mw_O         )
   REAL, PARAMETER  :: mw_N2O5    = (2*mw_N  + 5*mw_O         )
   REAL, PARAMETER  :: mw_HNO3    = (  mw_H  +   mw_N + 3*mw_O)

! For Chem Groups
! ---------------
   integer, parameter ::                   NUMGRP  =  2
   character(len=10)  ::  chem_group_names(NUMGRP)  ! See Set_GMI_Globals and Set_GCC_Globals

   integer, parameter ::      MAXGRP_ELEM = 10      ! Max number of species per family

   INTEGER   :: sgrp_elem_map(MAXGRP_ELEM, NUMGRP)  ! See Set_GMI_Globals and Set_GCC_Globals
   REAL      ::      sgrp_fac(MAXGRP_ELEM, NUMGRP)  ! See Set_GMI_Globals and Set_GCC_Globals



   TYPE ::  ChemGroupElemPtr
     REAL, DIMENSION(:,:,:), POINTER :: p
   END TYPE ChemGroupElemPtr

   INTEGER, PARAMETER      ::                   GMI_GROUP_SPECIES_COUNT = 18
   TYPE (ChemGroupElemPtr) ::                 x(GMI_GROUP_SPECIES_COUNT)  ! vector of pointers
   character(len=10)       :: gmi_group_species(GMI_GROUP_SPECIES_COUNT)  ! See Set_GMI_Globals

   INTEGER, PARAMETER      ::                   GCC_GROUP_SPECIES_COUNT = 22
   TYPE (ChemGroupElemPtr) ::                 y(GCC_GROUP_SPECIES_COUNT)  ! vector of pointers
   character(len=10)       :: gcc_group_species(GCC_GROUP_SPECIES_COUNT)  ! See Set_GCC_Globals

  ! See Set_GMI_Globals for a subset of these:
  ! See Set_GCC_Globals for a subset of these:
  ! These indices are relative to vector x (for GMI) or vector y (for GCC)
  integer :: iBR
  integer :: iBR2
  integer :: iBRCL
  integer :: iBRO
  integer :: iBRONO2
  integer :: iBRNO2
  integer :: iBRNO3
  integer :: iHBR
  integer :: iHOBR

  integer :: iCL
  integer :: iCL2
  integer :: iCLO
  integer :: iCLOO
  integer :: iCL2O2
  integer :: iCLONO2  ! GMI version
  integer :: iCLNO2
  integer :: iCLNO3   ! GCC version
  integer :: iHCL
  integer :: iHOCL
  integer :: iOCLO

  integer :: iN2O5
  integer :: iNO
  integer :: iNO2
  integer :: iHNO3


!
! !DESCRIPTION:
!
!  This module provides types and subroutines pertaining to chemical families
!
! !REVISION HISTORY:
!
!  31Jul2018 - Manyin, first crack
!  21Nov2018 - Manyin: Resolve Cl from BrCl (and N from ClONO2) within cell,
!                      and then distribute remainder within column
!
!EOP

!-------------------------------------------------------------------------
CONTAINS

!---------------------------------------------------------------------------
! NASA/GSFC, Atmospheric Chemistry and Dynamics Lab,  Code 614             !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  Init_GMI_Chem_Groups - One time initialization
!
! !INTERFACE:

 SUBROUTINE Init_GMI_Chem_Groups()

  IMPLICIT NONE

   GMI_groups_active = .TRUE.

  RETURN
 END SUBROUTINE Init_GMI_Chem_Groups

!---------------------------------------------------------------------------
! NASA/GSFC, Atmospheric Chemistry and Dynamics Lab,  Code 614             !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  Init_GCC_Chem_Groups - One time initialization
!
! !INTERFACE:

 SUBROUTINE Init_GCC_Chem_Groups()

  IMPLICIT NONE

   GCC_groups_active = .TRUE.

  RETURN
 END SUBROUTINE Init_GCC_Chem_Groups

!---------------------------------------------------------------------------
! NASA/GSFC, Atmospheric Chemistry and Dynamics Lab,  Code 614             !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  Set_GMI_Globals
!
! !INTERFACE:

 SUBROUTINE Set_GMI_Globals()

  IMPLICIT NONE

  chem_group_names(1) = 'Bry'
  chem_group_names(2) = 'Cly'

  iBR         = 1   ; gmi_group_species(iBR)     = 'Br    '
  iBRCL       = 2   ; gmi_group_species(iBRCL)   = 'BrCl  '
  iBRO        = 3   ; gmi_group_species(iBRO)    = 'BrO   '
  iBRONO2     = 4   ; gmi_group_species(iBRONO2) = 'BrONO2'
  iHBR        = 5   ; gmi_group_species(iHBR)    = 'HBr   '
  iHOBR       = 6   ; gmi_group_species(iHOBR)   = 'HOBr  '

  iCL         = 7   ; gmi_group_species(iCL)     = 'Cl    '
  iCL2        = 8   ; gmi_group_species(iCL2)    = 'Cl2   '
  iCLO        = 9   ; gmi_group_species(iCLO)    = 'ClO   '
  iCL2O2      = 10  ; gmi_group_species(iCL2O2)  = 'Cl2O2 '
  iCLONO2     = 11  ; gmi_group_species(iCLONO2) = 'ClONO2'
  iHCL        = 12  ; gmi_group_species(iHCL)    = 'HCl   '
  iHOCL       = 13  ; gmi_group_species(iHOCL)   = 'HOCl  '
  iOCLO       = 14  ; gmi_group_species(iOCLO)   = 'OClO  '

  iN2O5       = 15  ; gmi_group_species(iN2O5)   = 'N2O5  '
  iNO         = 16  ; gmi_group_species(iNO)     = 'NO    '
  iNO2        = 17  ; gmi_group_species(iNO2)    = 'NO2   '
  iHNO3       = 18  ; gmi_group_species(iHNO3)   = 'HNO3  ' 

  ! Species indices and factors for Bry
  ! Factors account for the units being mol/mol
  ! (indices are relative to the vector x)
  sgrp_elem_map( 1,1) = iBR      ; sgrp_fac( 1,1) = 1.0
  sgrp_elem_map( 2,1) = iBRCL    ; sgrp_fac( 2,1) = 1.0
  sgrp_elem_map( 3,1) = iBRO     ; sgrp_fac( 3,1) = 1.0
  sgrp_elem_map( 4,1) = iBRONO2  ; sgrp_fac( 4,1) = 1.0
  sgrp_elem_map( 5,1) = iHBR     ; sgrp_fac( 5,1) = 1.0
  sgrp_elem_map( 6,1) = iHOBR    ; sgrp_fac( 6,1) = 1.0
  sgrp_elem_map( 7,1) = 0        ; sgrp_fac( 7,1) = 0.0
  sgrp_elem_map( 8,1) = 0        ; sgrp_fac( 8,1) = 0.0
  sgrp_elem_map( 9,1) = 0        ; sgrp_fac( 9,1) = 0.0
  sgrp_elem_map(10,1) = 0        ; sgrp_fac(10,1) = 0.0

  ! Species indices and factors for Cly
  ! Factors account for the units being mol/mol
  ! (indices are relative to the vector x)
  sgrp_elem_map( 1,2) = iCL      ; sgrp_fac( 1,2) = 1.0
  sgrp_elem_map( 2,2) = iCL2     ; sgrp_fac( 2,2) = 2.0
  sgrp_elem_map( 3,2) = iCLO     ; sgrp_fac( 3,2) = 1.0
  sgrp_elem_map( 4,2) = iCL2O2   ; sgrp_fac( 4,2) = 2.0
  sgrp_elem_map( 5,2) = iCLONO2  ; sgrp_fac( 5,2) = 1.0
  sgrp_elem_map( 6,2) = iHCL     ; sgrp_fac( 6,2) = 1.0
  sgrp_elem_map( 7,2) = iHOCL    ; sgrp_fac( 7,2) = 1.0
  sgrp_elem_map( 8,2) = iOCLO    ; sgrp_fac( 8,2) = 1.0
  sgrp_elem_map( 9,2) = 0        ; sgrp_fac( 9,2) = 0.0
  sgrp_elem_map(10,2) = 0        ; sgrp_fac(10,2) = 0.0

  RETURN
 END SUBROUTINE Set_GMI_Globals


!---------------------------------------------------------------------------
! NASA/GSFC, Atmospheric Chemistry and Dynamics Lab,  Code 614             !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  Set_GCC_Globals
!
! !INTERFACE:

 SUBROUTINE Set_GCC_Globals()

  IMPLICIT NONE

  chem_group_names(1) = 'TRC_Bry'
  chem_group_names(2) = 'TRC_Cly'

  iBR         = 1   ; gcc_group_species(iBR)     = 'TRC_Br    '
  iBR2        = 2   ; gcc_group_species(iBR2)    = 'TRC_Br2   '
  iBRCL       = 3   ; gcc_group_species(iBRCL)   = 'TRC_BrCl  '
  iBRNO2      = 4   ; gcc_group_species(iBRNO2)  = 'TRC_BrNO2 '
  iBRNO3      = 5   ; gcc_group_species(iBRNO3)  = 'TRC_BrNO3 '
  iBRO        = 6   ; gcc_group_species(iBRO)    = 'TRC_BrO   '
  iHBR        = 7   ; gcc_group_species(iHBR)    = 'TRC_HBr   '
  iHOBR       = 8   ; gcc_group_species(iHOBR)   = 'TRC_HOBr  '

  iCL         = 9   ; gcc_group_species(iCL)     = 'TRC_Cl    '
  iCL2        = 10  ; gcc_group_species(iCL2)    = 'TRC_Cl2   '
  iCL2O2      = 11  ; gcc_group_species(iCL2O2)  = 'TRC_Cl2O2 '
  iCLNO2      = 12  ; gcc_group_species(iCLNO2)  = 'TRC_ClNO2 '
  iCLNO3      = 13  ; gcc_group_species(iCLNO3)  = 'TRC_ClNO3 '
  iCLO        = 14  ; gcc_group_species(iCLO)    = 'TRC_ClO   '
  iCLOO       = 15  ; gcc_group_species(iCLOO)   = 'TRC_ClOO  '
  iHCL        = 16  ; gcc_group_species(iHCL)    = 'TRC_HCl   '
  iHOCL       = 17  ; gcc_group_species(iHOCL)   = 'TRC_HOCl  '
  iOCLO       = 18  ; gcc_group_species(iOCLO)   = 'TRC_OClO  '

  iN2O5       = 19  ; gcc_group_species(iN2O5)   = 'TRC_N2O5  '
  iNO         = 20  ; gcc_group_species(iNO)     = 'TRC_NO    '
  iNO2        = 21  ; gcc_group_species(iNO2)    = 'TRC_NO2   '
  iHNO3       = 22  ; gcc_group_species(iHNO3)   = 'TRC_HNO3  ' 

  ! Species indices and factors for Bry
  ! Factors account for the units being kg/kg
  ! (indices are relative to the vector y)
  sgrp_elem_map( 1,1) = iBR      ; sgrp_fac( 1,1) =   mw_BR / mw_BR
  sgrp_elem_map( 2,1) = iBR2     ; sgrp_fac( 2,1) = 2*mw_BR / mw_BR2
  sgrp_elem_map( 3,1) = iBRCL    ; sgrp_fac( 3,1) =   mw_BR / mw_BRCL
  sgrp_elem_map( 4,1) = iBRNO2   ; sgrp_fac( 4,1) =   mw_BR / mw_BRNO2
  sgrp_elem_map( 5,1) = iBRNO3   ; sgrp_fac( 5,1) =   mw_BR / mw_BRNO3
  sgrp_elem_map( 6,1) = iBRO     ; sgrp_fac( 6,1) =   mw_BR / mw_BRO
  sgrp_elem_map( 7,1) = iHBR     ; sgrp_fac( 7,1) =   mw_BR / mw_HBR
  sgrp_elem_map( 8,1) = iHOBR    ; sgrp_fac( 8,1) =   mw_BR / mw_HOBR
  sgrp_elem_map( 9,1) = 0        ; sgrp_fac( 9,1) =   0.0
  sgrp_elem_map(10,1) = 0        ; sgrp_fac(10,1) =   0.0

  ! Species indices and factors for Cly
  ! Factors account for the units being kg/kg
  ! (indices are relative to the vector y)
  sgrp_elem_map( 1,2) = iCL      ; sgrp_fac( 1,2) =   mw_CL / mw_CL
  sgrp_elem_map( 2,2) = iCL2     ; sgrp_fac( 2,2) = 2*mw_CL / mw_CL2
  sgrp_elem_map( 3,2) = iCL2O2   ; sgrp_fac( 3,2) = 2*mw_CL / mw_CL2O2
  sgrp_elem_map( 4,2) = iCLNO2   ; sgrp_fac( 4,2) =   mw_CL / mw_CLNO2
  sgrp_elem_map( 5,2) = iCLNO3   ; sgrp_fac( 5,2) =   mw_CL / mw_CLNO3
  sgrp_elem_map( 6,2) = iCLO     ; sgrp_fac( 6,2) =   mw_CL / mw_CLO
  sgrp_elem_map( 7,2) = iCLOO    ; sgrp_fac( 7,2) =   mw_CL / mw_CLOO
  sgrp_elem_map( 8,2) = iHCL     ; sgrp_fac( 8,2) =   mw_CL / mw_HCL
  sgrp_elem_map( 9,2) = iHOCL    ; sgrp_fac( 9,2) =   mw_CL / mw_HOCL
  sgrp_elem_map(10,2) = iOCLO    ; sgrp_fac(10,2) =   mw_CL / mw_OCLO

  RETURN
 END SUBROUTINE Set_GCC_Globals



!---------------------------------------------------------------------------
! NASA/GSFC, Atmospheric Chemistry and Dynamics Lab,  Code 614             !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  Pack_Chem_Groups - Initialize the family entries in the Advection Bundle
!
! !INTERFACE:

 SUBROUTINE Pack_Chem_Groups(state)

  IMPLICIT NONE

! !ARGUMENTS
  TYPE(ESMF_State), INTENT(in) :: state

   character(len=ESMF_MAXSTR)            :: IAm

   REAL, POINTER, DIMENSION(:,:,:)       :: qq1  ! for chem groups
   INTEGER                               :: ig,im,imsgrp,i
   REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: qq2
   TYPE(ESMF_FieldBundle)                :: TRADV_BUNDLE

   INTEGER                               :: STATUS, RC
   INTEGER                               :: i1,i2,j1,j2,k1,k2

   Iam = "Pack_Chem_Groups"


    IF ( GMI_groups_active ) THEN

      ! Store in Bry:  atoms of Br / molecules of air
      ! Store in Cly:  atoms of Cl / molecules of air

      CALL Set_GMI_Globals()

      CALL ESMF_StateGet(state, 'TRADV', TRADV_BUNDLE, __RC__ )

      DO i=1,GMI_GROUP_SPECIES_COUNT
        CALL ESMFL_BundleGetPointerToData(TRADV_BUNDLE, 'GMICHEM::'//TRIM(gmi_group_species(i)), x(i)%p, __RC__ )
      ENDDO

      i1 = LBOUND(x(1)%p,1); i2 = UBOUND(x(1)%p,1)
      j1 = LBOUND(x(1)%p,2); j2 = UBOUND(x(1)%p,2)
      k1 = LBOUND(x(1)%p,3); k2 = UBOUND(x(1)%p,3)

      allocate( qq2(i1:i2,j1:j2,k1:k2), __STAT__ )

      DO ig=1,NUMGRP
        CALL ESMFL_BundleGetPointerToData(TRADV_BUNDLE, 'GMICHEM::'//TRIM(chem_group_names(ig)), qq1, __RC__ )

!       PACK GROUP  using REAL*8 if possible

        ! Accumulate in qq2 which is double precision
        qq2(:,:,:) = 0.0

        do im = 1, MAXGRP_ELEM

          imsgrp = sgrp_elem_map(im,ig)

          if (imsgrp > 0) then
            qq2(:,:,:) = qq2(:,:,:) + x(imsgrp)%p(:,:,:) * sgrp_fac(im,ig)
          end if

        end do

        qq1 = qq2

      ENDDO

      deallocate(qq2, __STAT__ )

    ENDIF

    IF ( GCC_groups_active ) THEN

      ! Store in TRC_Bry:  kg of Br / kg of air
      ! Store in TRC_Cly:  kg of Cl / kg of air

      CALL Set_GCC_Globals()

      CALL ESMF_StateGet(state, 'TRADV', TRADV_BUNDLE, __RC__ )

      DO i=1,GCC_GROUP_SPECIES_COUNT
        CALL ESMFL_BundleGetPointerToData(TRADV_BUNDLE, 'GEOSCHEMCHEM::'//TRIM(gcc_group_species(i)), y(i)%p, __RC__ )
      ENDDO

      i1 = LBOUND(y(1)%p,1); i2 = UBOUND(y(1)%p,1)
      j1 = LBOUND(y(1)%p,2); j2 = UBOUND(y(1)%p,2)
      k1 = LBOUND(y(1)%p,3); k2 = UBOUND(y(1)%p,3)

      allocate( qq2(i1:i2,j1:j2,k1:k2), __STAT__ )

      DO ig=1,NUMGRP
        CALL ESMFL_BundleGetPointerToData(TRADV_BUNDLE, 'GEOSCHEMCHEM::'//TRIM(chem_group_names(ig)), qq1, __RC__ )

!       PACK GROUP  using REAL*8 if possible

        ! Accumulate in qq2 which is double precision
        qq2(:,:,:) = 0.0

        do im = 1, MAXGRP_ELEM

          imsgrp = sgrp_elem_map(im,ig)

          if (imsgrp > 0) then
            qq2(:,:,:) = qq2(:,:,:) + y(imsgrp)%p(:,:,:) * sgrp_fac(im,ig)
          end if

        end do

        qq1 = qq2

      ENDDO

      deallocate(qq2, __STAT__ )

    ENDIF

  RETURN
 END SUBROUTINE Pack_Chem_Groups

!---------------------------------------------------------------------------
! NASA/GSFC, Atmospheric Chemistry and Dynamics Lab,  Code 614             !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  Unpack_Chem_Groups - Modify species based on the family tracers
!
! !INTERFACE:

#ifdef CHEM_INFO
 SUBROUTINE Unpack_Chem_Groups(state, PLE, AREA, Q_separate, bry_ratio, aBRCL, aCL2, aOCLO, aCL2O2, aCLO, aHCL, aHOCL, zBRCL, zCL2, zOCLO, zCL2O2, zCLO, zHCL, zHOCL, zBRY, aCLY, zCLY )
#else
 SUBROUTINE Unpack_Chem_Groups(state, PLE, AREA, Q_separate)
#endif

  IMPLICIT NONE

! !ARGUMENTS
   type(ESMF_State), intent(in) :: state
   REAL*4, POINTER, DIMENSION(:,:,:), INTENT(IN)           :: PLE
   REAL*4, POINTER, DIMENSION(:,:),   INTENT(IN)           :: AREA
   REAL*4, POINTER, DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: Q_separate  ! water vapor  [kg vapor / kg moist air]
                                                                          ! option to provide Q separate from the advection bundle,
                                                                          ! otherwise we get Q from that bundle (TRADV)
#ifdef CHEM_INFO
   REAL*4, POINTER, DIMENSION(:,:,:)  :: bry_ratio, aBRCL, aCL2, aOCLO, aCL2O2, aCLO, aHCL, aHOCL
   REAL*4, POINTER, DIMENSION(:,:,:)  ::            zBRCL, zCL2, zOCLO, zCL2O2, zCLO, zHCL, zHOCL
   REAL*4, POINTER, DIMENSION(:,:,:)  ::            zBRY,  aCLY, zCLY
#endif


   character(len=ESMF_MAXSTR)            :: IAm

   REAL, POINTER, DIMENSION(:,:,:)       :: qq1  ! for chem groups
   INTEGER                               :: ig,im,imsgrp,i
   REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: group_factor
   REAL*8, ALLOCATABLE, DIMENSION(:,:,:) ::   BRCL_pre_adjust
   REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: CLONO2_pre_adjust
   REAL*8, ALLOCATABLE, DIMENSION(:,:,:) ::  CLNO2_pre_adjust
   REAL*8, ALLOCATABLE, DIMENSION(:,:,:) ::  CLNO3_pre_adjust

   REAL*8, POINTER,     DIMENSION(:,:,:) :: XXX_diff
   REAL*8, POINTER,     DIMENSION(:,:,:) :: BRCL_diff   ! alias of XXX_diff
   REAL*8, POINTER,     DIMENSION(:,:,:) :: CLNOX_diff  ! alias of XXX_diff
   REAL*8, POINTER,     DIMENSION(:,:,:) :: CLONO2_diff ! alias of XXX_diff

   REAL*8, POINTER,     DIMENSION(:,:,:) :: XXX_sum
   REAL*8, POINTER,     DIMENSION(:,:,:) :: CL_sum      ! alias of XXX_sum
   REAL*8, POINTER,     DIMENSION(:,:,:) :: N_sum       ! alias of XXX_sum

   REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: qq2
   REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: mass_AIR ! kg moist air
   REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: kmol_per_kg_AIR ! conversion for moist air
   REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: frac     ! scaling fraction
   REAL*8                                :: minloc_out(3)
   TYPE(ESMF_FieldBundle)                :: TRADV_BUNDLE
   REAL*4,     POINTER, DIMENSION(:,:,:) :: Q    ! water vapor  [kg vapor / kg moist air]

   LOGICAL, ALLOCATABLE, DIMENSION(:)    :: in_play  ! whether gridbox can accept more adjustment

   INTEGER                               :: STATUS, RC
   INTEGER                               :: i1,i2,j1,j2,k1,k2
   INTEGER                               :: k

   REAL, PARAMETER                       :: MAXFRAC = 0.99   ! Do not reduce more than this
   REAL, PARAMETER                       :: VERY_SMALL = 1.0e-25

   Iam = "Unpack_Chem_Groups"

   IF ( GMI_groups_active ) THEN

!     IF (MAPL_AM_I_ROOT()) print*,'Unpacking GMI groups Bry and Cly'

      CALL Set_GMI_Globals()

      i1 = LBOUND(x(1)%p,1); i2 = UBOUND(x(1)%p,1)
      j1 = LBOUND(x(1)%p,2); j2 = UBOUND(x(1)%p,2)
      k1 = LBOUND(x(1)%p,3); k2 = UBOUND(x(1)%p,3)

      allocate(     group_factor(i1:i2,j1:j2,k1:k2), &
                 BRCL_pre_adjust(i1:i2,j1:j2,k1:k2), &
               CLONO2_pre_adjust(i1:i2,j1:j2,k1:k2), &
                        mass_AIR(i1:i2,j1:j2,k1:k2), &
                 kmol_per_kg_AIR(i1:i2,j1:j2,k1:k2), &
                            frac(i1:i2,j1:j2,k1:k2), &
                        XXX_diff(i1:i2,j1:j2,k1:k2), &
                         XXX_sum(i1:i2,j1:j2,k1:k2), &
                                     in_play(k1:k2), &
                             qq2(i1:i2,j1:j2,k1:k2),   __STAT__ )

      do k=k1,k2
        mass_AIR(:,:,k) = (PLE(:,:,k)-PLE(:,:,k-1)) * AREA(:,:) / MAPL_GRAV
      end do

      CALL ESMF_StateGet(state, 'TRADV', TRADV_BUNDLE, __RC__ )

      IF ( PRESENT(Q_separate) ) THEN
        Q => Q_separate
      ELSE
        CALL ESMFL_BundleGetPointerToData(TRADV_BUNDLE, 'Q', Q, __RC__ )
      ENDIF

      kmol_per_kg_AIR(:,:,:) = (1.0 - Q(:,:,:))/MAPL_AIRMW + Q(:,:,:)/MAPL_H2OMW

      DO ig=1,NUMGRP
        CALL ESMFL_BundleGetPointerToData(TRADV_BUNDLE, 'GMICHEM::'//TRIM(chem_group_names(ig)), qq1, __RC__ )
!
!  BASED ON INCLUDE FILE  vvv
!

      qq2(:,:,:)                   = 0.0d0

      select case (ig)
        case (1)

!.... Bry

#ifdef CHEM_INFO
  IF ( ASSOCIATED(zBRY)  )   zBRY  = qq1
  IF ( ASSOCIATED(aBRCL) )   aBRCL = x(iBRCL)%p(:,:,:)
#endif

        BRCL_pre_adjust(:,:,:)     = x(iBRCL)%p(:,:,:)

          do im                    = 1, MAXGRP_ELEM
            imsgrp                 = sgrp_elem_map(im,ig)
            if (imsgrp > 0)  &
     &        qq2(:,:,:)           = qq2(:,:,:) +  &
     &                               x(imsgrp)%p(:,:,:) *  &
     &                               sgrp_fac(im,ig)
          end do

          where ( qq2(:,:,:) > VERY_SMALL )
            group_factor(:,:,:)      = qq1(:,:,:) / qq2(:,:,:)
          elsewhere
            group_factor(:,:,:)      = -999.0
          end where

          do im                    = 1, MAXGRP_ELEM
            imsgrp                 = sgrp_elem_map(im,ig)
            if (imsgrp > 0)  &
     &        x(imsgrp)%p(:,:,:)   = x(imsgrp)%p(:,:,:) *  &
     &                               group_factor(:,:,:)
          end do

!.... In the very rare cases where the species sum is tiny, divide up Bry
          where( group_factor(:,:,:) == -999.0 )
            x(iBR    )%p(:,:,:) = qq1(:,:,:) / 6.0
            x(iBRCL  )%p(:,:,:) = qq1(:,:,:) / 6.0
            x(iBRO   )%p(:,:,:) = qq1(:,:,:) / 6.0
            x(iBRONO2)%p(:,:,:) = qq1(:,:,:) / 6.0
            x(iHBR   )%p(:,:,:) = qq1(:,:,:) / 6.0
            x(iHOBR  )%p(:,:,:) = qq1(:,:,:) / 6.0
          end where

#ifdef CHEM_INFO
  IF ( ASSOCIATED(zBRCL) )   zBRCL   = x(iBRCL)%p(:,:,:)
  IF ( ASSOCIATED(bry_ratio) )   bry_ratio = group_factor
#endif

        case (2)

!.... Cly  -- NOTE:  family does NOT include BrCl

#ifdef CHEM_INFO
        IF ( ASSOCIATED(aCLY) )   aCLY = qq1
#endif

        CLONO2_pre_adjust(:,:,:)   = x(iCLONO2)%p(:,:,:)

        do im                      = 1, MAXGRP_ELEM
          imsgrp                   = sgrp_elem_map(im,ig)
          if (imsgrp > 0)  &
     &      qq2(:,:,:)             = qq2(:,:,:) +  &
     &                               x(imsgrp)%p(:,:,:) *  &
     &                               sgrp_fac(im,ig)
        end do

        where ( qq2(:,:,:) > VERY_SMALL )
          group_factor(:,:,:)      = qq1(:,:,:) / qq2(:,:,:)
        elsewhere
          group_factor(:,:,:)      = -999.0
        end where

        do im                      = 1, MAXGRP_ELEM
          imsgrp                   = sgrp_elem_map(im,ig)
          if (imsgrp > 0)  &
     &      x(imsgrp)%p(:,:,:)     = x(imsgrp)%p(:,:,:) *  &
     &                               group_factor(:,:,:)
        end do

!.... In the very rare cases where the species sum is tiny, divide up Cly
          where( group_factor(:,:,:) == -999.0 )
            x(iCL    )%p(:,:,:) = (qq1(:,:,:) / 10.0)
            x(iCL2   )%p(:,:,:) = (qq1(:,:,:) / 10.0) * 2.0
            x(iCLO   )%p(:,:,:) = (qq1(:,:,:) / 10.0)
            x(iCL2O2 )%p(:,:,:) = (qq1(:,:,:) / 10.0) * 2.0
            x(iCLONO2)%p(:,:,:) = (qq1(:,:,:) / 10.0)
            x(iHCL   )%p(:,:,:) = (qq1(:,:,:) / 10.0)
            x(iHOCL  )%p(:,:,:) = (qq1(:,:,:) / 10.0)
            x(iOCLO  )%p(:,:,:) = (qq1(:,:,:) / 10.0)
          end where

      end select

      ENDDO

!!!
!!!  Compensate for the Cl changed in BrCl
!!!

#ifdef CHEM_INFO
  IF ( ASSOCIATED(aCL2) )   aCL2   = x(iCL2)%p(:,:,:)
  IF ( ASSOCIATED(aOCLO) )  aOCLO  = x(iOCLO)%p(:,:,:)
  IF ( ASSOCIATED(aCL2O2) ) aCL2O2 = x(iCL2O2)%p(:,:,:)
  IF ( ASSOCIATED(aCLO) )   aCLO   = x(iCLO)%p(:,:,:)
  IF ( ASSOCIATED(aHCL) )   aHCL   = x(iHCL)%p(:,:,:)
  IF ( ASSOCIATED(aHOCL) )  aHOCL  = x(iHOCL)%p(:,:,:)
#endif

        BRCL_diff => XXX_diff
          CL_sum  => XXX_sum

!.... Compute delta BrCl [mol_BRCL / mol_AIR] = [mol_CL / mol_AIR]
        BRCL_diff(:,:,:) = x(iBRCL)%p(:,:,:) - BRCL_pre_adjust(:,:,:)

!.... Compute   CL_sum   [mol_CL   / mol_AIR]
        CL_sum(:,:,:) = 2 * x(iCL2  )%p(:,:,:) +  &
     &                      x(iOCLO )%p(:,:,:) +  &
     &                  2 * x(iCL2O2)%p(:,:,:) +  &
     &                      x(iCLO  )%p(:,:,:) +  &
     &                      x(iHCL  )%p(:,:,:) +  &
     &                      x(iHOCL )%p(:,:,:)

!print*,'CL_sum min = ', MINVAL(CL_sum)

!.... In gridboxes where BrCl has been reduced, increase Cly species
        where ( BRCL_diff(:,:,:) < 0.0 )
          where ( CL_sum(:,:,:) > VERY_SMALL )

            frac(:,:,:) = -1.0 * BRCL_diff(:,:,:) / CL_sum(:,:,:)  ! Note frac is positive
                                                                   ! Multiply by  (1 + frac)

            x(iCL2  )%p(:,:,:) = x(iCL2  )%p(:,:,:) * ( 1.0 + frac(:,:,:) )
            x(iOCLO )%p(:,:,:) = x(iOCLO )%p(:,:,:) * ( 1.0 + frac(:,:,:) )
            x(iCL2O2)%p(:,:,:) = x(iCL2O2)%p(:,:,:) * ( 1.0 + frac(:,:,:) )
            x(iCLO  )%p(:,:,:) = x(iCLO  )%p(:,:,:) * ( 1.0 + frac(:,:,:) )
            x(iHCL  )%p(:,:,:) = x(iHCL  )%p(:,:,:) * ( 1.0 + frac(:,:,:) )
            x(iHOCL )%p(:,:,:) = x(iHOCL )%p(:,:,:) * ( 1.0 + frac(:,:,:) )

          elsewhere

            ! Add to species (subtract a negative):
            x(iCL2  )%p(:,:,:) = x(iCL2  )%p(:,:,:) - BRCL_diff(:,:,:) / 12.0
            x(iOCLO )%p(:,:,:) = x(iOCLO )%p(:,:,:) - BRCL_diff(:,:,:) /  6.0
            x(iCL2O2)%p(:,:,:) = x(iCL2O2)%p(:,:,:) - BRCL_diff(:,:,:) / 12.0
            x(iCLO  )%p(:,:,:) = x(iCLO  )%p(:,:,:) - BRCL_diff(:,:,:) /  6.0
            x(iHCL  )%p(:,:,:) = x(iHCL  )%p(:,:,:) - BRCL_diff(:,:,:) /  6.0
            x(iHOCL )%p(:,:,:) = x(iHOCL )%p(:,:,:) - BRCL_diff(:,:,:) /  6.0

          end where

        end where

!.... In gridboxes where BrCl has increased, reduce Cly species

        call GMI_reduce_Cly_vmr ( )

!.... Clean up
        nullify(BRCL_diff, CL_sum)

#ifdef CHEM_INFO
  IF ( ASSOCIATED(zCL2) )   zCL2   = x(iCL2)%p(:,:,:)
  IF ( ASSOCIATED(zOCLO) )  zOCLO  = x(iOCLO)%p(:,:,:)
  IF ( ASSOCIATED(zCL2O2) ) zCL2O2 = x(iCL2O2)%p(:,:,:)
  IF ( ASSOCIATED(zCLO) )   zCLO   = x(iCLO)%p(:,:,:)
  IF ( ASSOCIATED(zHCL) )   zHCL   = x(iHCL)%p(:,:,:)
  IF ( ASSOCIATED(zHOCL) )  zHOCL  = x(iHOCL)%p(:,:,:)
#endif

!!!
!!!  Done handling Cl change from BrCl
!!!

!.... Account for changes in NOy reservoir ClONO2

        CLONO2_diff => XXX_diff
        N_sum       => XXX_sum

!.... Compute delta ClONO2 [mol_CLONO2 / mol_AIR] = [mol_N / mol_AIR]
        CLONO2_diff(:,:,:) = x(iCLONO2)%p(:,:,:) - CLONO2_pre_adjust(:,:,:)

!.... Compute N_sum     [mol_N / mol_AIR]
        N_sum(:,:,:) = 2 * x(iN2O5)%p(:,:,:) +  &
                           x(iNO  )%p(:,:,:) +  &
                           x(iNO2 )%p(:,:,:) +  &
                           x(iHNO3)%p(:,:,:)

!print*,'N_sum min = ', MINVAL(N_sum)

!.... In gridboxes where ClONO2 has been reduced, increase N species
        where ( CLONO2_diff(:,:,:) < 0.0 )
          where ( N_sum(:,:,:) > VERY_SMALL )

            frac(:,:,:) = -1.0 * CLONO2_diff(:,:,:) / N_sum(:,:,:)  ! Note frac is positive
                                                                    ! Multiply by  (1 + frac)

            x(iN2O5)%p(:,:,:) = x(iN2O5)%p(:,:,:) * ( 1.0 + frac(:,:,:) )
            x(iNO  )%p(:,:,:) = x(iNO  )%p(:,:,:) * ( 1.0 + frac(:,:,:) )
            x(iNO2 )%p(:,:,:) = x(iNO2 )%p(:,:,:) * ( 1.0 + frac(:,:,:) )
            x(iHNO3)%p(:,:,:) = x(iHNO3)%p(:,:,:) * ( 1.0 + frac(:,:,:) )

          elsewhere

            ! Add to species (subtract a negative):
            x(iN2O5)%p(:,:,:) = x(iN2O5)%p(:,:,:) - CLONO2_diff(:,:,:) / 8.0
            x(iNO  )%p(:,:,:) = x(iNO  )%p(:,:,:) - CLONO2_diff(:,:,:) / 4.0
            x(iNO2 )%p(:,:,:) = x(iNO2 )%p(:,:,:) - CLONO2_diff(:,:,:) / 4.0
            x(iHNO3)%p(:,:,:) = x(iHNO3)%p(:,:,:) - CLONO2_diff(:,:,:) / 4.0

          end where
        end where


!.... In gridboxes where ClONO2 has increased, reduce N species

        call GMI_reduce_N_vmr ( )

!.... Clean up
        nullify(CLONO2_diff, N_sum)

!
!  ^^^^^^^^^^^^^^^^^^^^^^^^^^
!

#ifdef CHEM_INFO
      IF ( ASSOCIATED(zCLY) ) THEN

        qq2(:,:,:) = 0.0
        ig = 2

        do im                      = 1, MAXGRP_ELEM
          imsgrp                   = sgrp_elem_map(im,ig)
          if (imsgrp > 0)  &
     &      qq2(:,:,:)             = qq2(:,:,:) +  &
     &                               x(imsgrp)%p(:,:,:) *  &
     &                               sgrp_fac(im,ig)
        end do

        zCLY = qq2

      ENDIF
#endif

      deallocate(group_factor, BRCL_pre_adjust, CLONO2_pre_adjust, mass_AIR,      &
                 kmol_per_kg_AIR, in_play, frac, XXX_diff, XXX_sum, qq2, __STAT__ )

   ENDIF


   IF ( GCC_groups_active ) THEN

!     IF (MAPL_AM_I_ROOT()) print*,'Unpacking GCC groups TRC_Bry and TRC_Cly'

      CALL Set_GCC_Globals()

      i1 = LBOUND(y(1)%p,1); i2 = UBOUND(y(1)%p,1)
      j1 = LBOUND(y(1)%p,2); j2 = UBOUND(y(1)%p,2)
      k1 = LBOUND(y(1)%p,3); k2 = UBOUND(y(1)%p,3)

      allocate(     group_factor(i1:i2,j1:j2,k1:k2), &
                 BRCL_pre_adjust(i1:i2,j1:j2,k1:k2), &
                CLNO2_pre_adjust(i1:i2,j1:j2,k1:k2), &
                CLNO3_pre_adjust(i1:i2,j1:j2,k1:k2), &
                        mass_AIR(i1:i2,j1:j2,k1:k2), &
                            frac(i1:i2,j1:j2,k1:k2), &
                        XXX_diff(i1:i2,j1:j2,k1:k2), &
                         XXX_sum(i1:i2,j1:j2,k1:k2), &
                                     in_play(k1:k2), &
                             qq2(i1:i2,j1:j2,k1:k2),   __STAT__ )

      do k=k1,k2
        mass_AIR(:,:,k) = (PLE(:,:,k)-PLE(:,:,k-1)) * AREA(:,:) / MAPL_GRAV
      end do

      CALL ESMF_StateGet(state, 'TRADV', TRADV_BUNDLE, __RC__ )

      DO ig=1,NUMGRP
        CALL ESMFL_BundleGetPointerToData(TRADV_BUNDLE, 'GEOSCHEMCHEM::'//TRIM(chem_group_names(ig)), qq1, __RC__ )
!
!  BASED ON INCLUDE FILE  vvv
!

      qq2(:,:,:)                   = 0.0d0

      select case (ig)
        case (1)

!.... Bry

#ifdef CHEM_INFO
  IF ( ASSOCIATED(zBRY)  )   zBRY  = qq1
  IF ( ASSOCIATED(aBRCL) )   aBRCL = y(iBRCL)%p(:,:,:)
#endif

        BRCL_pre_adjust(:,:,:)     = y(iBRCL)%p(:,:,:)

          do im                    = 1, MAXGRP_ELEM
            imsgrp                 = sgrp_elem_map(im,ig)
            if (imsgrp > 0)  &
     &        qq2(:,:,:)           = qq2(:,:,:) +  &
     &                               y(imsgrp)%p(:,:,:) *  &
     &                               sgrp_fac(im,ig)
          end do

          group_factor(:,:,:)      = qq1(:,:,:) / qq2(:,:,:)

          do im                    = 1, MAXGRP_ELEM
            imsgrp                 = sgrp_elem_map(im,ig)
            if (imsgrp > 0)  &
     &        y(imsgrp)%p(:,:,:)   = y(imsgrp)%p(:,:,:) *  &
     &                               group_factor(:,:,:)
          end do

#ifdef CHEM_INFO
  IF ( ASSOCIATED(zBRCL) )   zBRCL   = y(iBRCL)%p(:,:,:)
  IF ( ASSOCIATED(bry_ratio) )   bry_ratio = group_factor
#endif

        case (2)

!  aCLY is initilized below, so z-a shows affect of BrCl adjustment
!#ifdef CHEM_INFO
!        IF ( ASSOCIATED(aCLY) )   aCLY = qq1
!#endif

!.... Cly  -- NOTE:  family does NOT include BrCl

        CLNO2_pre_adjust(:,:,:)    = y(iCLNO2)%p(:,:,:)
        CLNO3_pre_adjust(:,:,:)    = y(iCLNO3)%p(:,:,:)

        do im                      = 1, MAXGRP_ELEM
          imsgrp                   = sgrp_elem_map(im,ig)
          if (imsgrp > 0)  &
     &      qq2(:,:,:)             = qq2(:,:,:) +  &
     &                               y(imsgrp)%p(:,:,:) *  &
     &                               sgrp_fac(im,ig)
        end do

        group_factor(:,:,:)        = qq1(:,:,:) / qq2(:,:,:)

        do im                      = 1, MAXGRP_ELEM
          imsgrp                   = sgrp_elem_map(im,ig)
          if (imsgrp > 0)  &
     &      y(imsgrp)%p(:,:,:)     = y(imsgrp)%p(:,:,:) *  &
     &                               group_factor(:,:,:)
        end do

      end select

      ENDDO

!!!!!!!!!!!!!!!!!!!!!!!
!!!  Compensate for the Cl changed in BrCl

#ifdef CHEM_INFO
  IF ( ASSOCIATED(aCL2) )   aCL2   = y(iCL2)%p(:,:,:)
  IF ( ASSOCIATED(aOCLO) )  aOCLO  = y(iOCLO)%p(:,:,:)
  IF ( ASSOCIATED(aCL2O2) ) aCL2O2 = y(iCL2O2)%p(:,:,:)
  IF ( ASSOCIATED(aCLO) )   aCLO   = y(iCLO)%p(:,:,:)
  IF ( ASSOCIATED(aHCL) )   aHCL   = y(iHCL)%p(:,:,:)
  IF ( ASSOCIATED(aHOCL) )  aHOCL  = y(iHOCL)%p(:,:,:)
#endif


        BRCL_diff => XXX_diff
          CL_sum  => XXX_sum

!.... Compute delta BrCl [kmol_BrCl / kg_AIR] = [kmol_Cl / kg_AIR]
        BRCL_diff(:,:,:) = (y(iBRCL)%p(:,:,:) - BRCL_pre_adjust(:,:,:)) / mw_BRCL

!.... Compute CL_sum     [kmol_Cl   / kg_AIR]
        CL_sum(:,:,:) = 2 * y(iCL2  )%p(:,:,:) / mw_CL2    +  &
                            y(iOCLO )%p(:,:,:) / mw_OCLO   +  &
                        2 * y(iCL2O2)%p(:,:,:) / mw_CL2O2  +  &
                            y(iCLO  )%p(:,:,:) / mw_CLO    +  &
                            y(iHCL  )%p(:,:,:) / mw_HCL    +  &
                            y(iHOCL )%p(:,:,:) / mw_HOCL

!print*,'CL_sum min = ', MINVAL(CL_sum)

!.... In gridboxes where BrCl has been reduced, increase Cly species
        where ( BRCL_diff(:,:,:) < 0.0 )
          where ( CL_sum(:,:,:) > VERY_SMALL )

            frac(:,:,:) = -1.0 * BRCL_diff(:,:,:) / CL_sum(:,:,:)  ! Note frac is positive
                                                                   ! Multiply by  (1 + frac)

            y(iCL2  )%p(:,:,:) = y(iCL2  )%p(:,:,:) * ( 1.0 + frac(:,:,:) )
            y(iOCLO )%p(:,:,:) = y(iOCLO )%p(:,:,:) * ( 1.0 + frac(:,:,:) )
            y(iCL2O2)%p(:,:,:) = y(iCL2O2)%p(:,:,:) * ( 1.0 + frac(:,:,:) )
            y(iCLO  )%p(:,:,:) = y(iCLO  )%p(:,:,:) * ( 1.0 + frac(:,:,:) )
            y(iHCL  )%p(:,:,:) = y(iHCL  )%p(:,:,:) * ( 1.0 + frac(:,:,:) )
            y(iHOCL )%p(:,:,:) = y(iHOCL )%p(:,:,:) * ( 1.0 + frac(:,:,:) )

          elsewhere

            ! Add to species (subtract a negative):
            y(iCL2  )%p(:,:,:) = y(iCL2  )%p(:,:,:) - BRCL_diff(:,:,:) / 12.0  * mw_CL2
            y(iOCLO )%p(:,:,:) = y(iOCLO )%p(:,:,:) - BRCL_diff(:,:,:) /  6.0  * mw_OCLO
            y(iCL2O2)%p(:,:,:) = y(iCL2O2)%p(:,:,:) - BRCL_diff(:,:,:) / 12.0  * mw_CL2O2
            y(iCLO  )%p(:,:,:) = y(iCLO  )%p(:,:,:) - BRCL_diff(:,:,:) /  6.0  * mw_CLO
            y(iHCL  )%p(:,:,:) = y(iHCL  )%p(:,:,:) - BRCL_diff(:,:,:) /  6.0  * mw_HCL
            y(iHOCL )%p(:,:,:) = y(iHOCL )%p(:,:,:) - BRCL_diff(:,:,:) /  6.0  * mw_HOCL

          end where

        end where


!.... In gridboxes where BrCl has increased, reduce Cly species

        call GCC_reduce_Cly_mmr ( )

!.... Clean up
        nullify(BRCL_diff, CL_sum)

#ifdef CHEM_INFO
  IF ( ASSOCIATED(zCL2) )   zCL2   = y(iCL2)%p(:,:,:)
  IF ( ASSOCIATED(zOCLO) )  zOCLO  = y(iOCLO)%p(:,:,:)
  IF ( ASSOCIATED(zCL2O2) ) zCL2O2 = y(iCL2O2)%p(:,:,:)
  IF ( ASSOCIATED(zCLO) )   zCLO   = y(iCLO)%p(:,:,:)
  IF ( ASSOCIATED(zHCL) )   zHCL   = y(iHCL)%p(:,:,:)
  IF ( ASSOCIATED(zHOCL) )  zHOCL  = y(iHOCL)%p(:,:,:)
#endif

!!!  Done handling Cl change from BrCl
!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!
!!!  Account for changes in NOy reservoirs ClNO2 and ClNO3

        CLNOX_diff => XXX_diff
        N_sum      => XXX_sum

!.... Compute delta ClNOX [(kmol_CLNO2 + kmol_CLNO3) / kg_AIR] = [kmol_N / kg_AIR]
        CLNOX_diff(:,:,:) = (y(iCLNO2)%p(:,:,:) - CLNO2_pre_adjust(:,:,:)) / mw_CLNO2 + &
                            (y(iCLNO3)%p(:,:,:) - CLNO3_pre_adjust(:,:,:)) / mw_CLNO3

!.... Compute N_sum     [kmol_N / kg_AIR]
        N_sum(:,:,:) = 2 * y(iN2O5)%p(:,:,:) / mw_N2O5  +  &
                           y(iNO  )%p(:,:,:) / mw_NO    +  &
                           y(iNO2 )%p(:,:,:) / mw_NO2   +  &
                           y(iHNO3)%p(:,:,:) / mw_HNO3

!print*,'N_sum min = ', MINVAL(N_sum)

!.... In gridboxes where ClNOX has been reduced, increase N species
        where ( CLNOX_diff(:,:,:) < 0.0 )
          where ( N_sum(:,:,:) > VERY_SMALL )

            frac(:,:,:) = -1.0 * CLNOX_diff(:,:,:) / N_sum(:,:,:)  ! Note frac is positive
                                                                   ! Multiply by  (1 + frac)

            y(iN2O5)%p(:,:,:) = y(iN2O5)%p(:,:,:) * ( 1.0 + frac(:,:,:) )
            y(iNO  )%p(:,:,:) = y(iNO  )%p(:,:,:) * ( 1.0 + frac(:,:,:) )
            y(iNO2 )%p(:,:,:) = y(iNO2 )%p(:,:,:) * ( 1.0 + frac(:,:,:) )
            y(iHNO3)%p(:,:,:) = y(iHNO3)%p(:,:,:) * ( 1.0 + frac(:,:,:) )

          elsewhere

            ! Add to species (subtract a negative):
            y(iN2O5)%p(:,:,:) = y(iN2O5)%p(:,:,:) - CLNOX_diff(:,:,:) / 8.0 * mw_N2O5
            y(iNO  )%p(:,:,:) = y(iNO  )%p(:,:,:) - CLNOX_diff(:,:,:) / 4.0 * mw_NO
            y(iNO2 )%p(:,:,:) = y(iNO2 )%p(:,:,:) - CLNOX_diff(:,:,:) / 4.0 * mw_NO2
            y(iHNO3)%p(:,:,:) = y(iHNO3)%p(:,:,:) - CLNOX_diff(:,:,:) / 4.0 * mw_HNO3

          end where

        end where


!.... In gridboxes where ClNOX has increased, reduce N species

        call GCC_reduce_N_mmr ( )

!.... Clean up
        nullify(CLNOX_diff, N_sum)

!!!  Done handling N change from ClNO2 and ClNO3
!!!!!!!!!!!!!!!!!!!!!!!

!
!  ^^^^^^^^^^^^^^^^^^^^^^^^^^
!

#ifdef CHEM_INFO
      IF ( ASSOCIATED(zCLY) ) THEN

        qq2(:,:,:) = 0.0
        ig = 2

        do im                      = 1, MAXGRP_ELEM
          imsgrp                   = sgrp_elem_map(im,ig)
          if (imsgrp > 0)  &
     &      qq2(:,:,:)             = qq2(:,:,:) +  &
     &                               y(imsgrp)%p(:,:,:) *  &
     &                               sgrp_fac(im,ig)
        end do

        zCLY = qq2

      ENDIF
#endif

      deallocate(group_factor, BRCL_pre_adjust, CLNO2_pre_adjust, CLNO3_pre_adjust,  &
                 mass_AIR, in_play, frac, XXX_diff, XXX_sum, qq2, __STAT__ )

   ENDIF

  RETURN


  contains

!
!   Where BrCl has increased, we want to reduce Cl to maintain conservation.
!   To distribute the Cl loss proportionally across 6 species, we compute
!   the fraction that is to be lost (guard against very low values
!   in all 6 species which would make the fraction too high or undefined).
!   Once we have the fraction by which to reduce the field, do not let
!   the field go negative; if mass loss is un-accounted for in a gridbox,
!   distribute it within the column.
!

!.... In gridboxes where BrCl has increased, reduce Cly species
    SUBROUTINE GMI_reduce_Cly_vmr ( )

       integer :: loop_count, i, j, k
       real*8  :: f, extra, vsum_EXTRA, col_mixrat, box_mixrat

        do i = i1,i2
          do j = j1,j2

            vsum_EXTRA = 0.0d0  !  extra within column  [kmol_CL]

            !! Guard against (1-frac) getting too small, or going negative

            !!
            !! First pass through the column; make sure that in_play(k) is initialized for all k
            !!
            do k = k1,k2

              if ( BRCL_diff(i,j,k) <= 0.0 ) then

                !! There was addition (not reduction) in this gridbox;
                !! permit reduction in the next phase (column work)
                !! if there is enough Cly
                if ( CL_sum(i,j,k) < VERY_SMALL ) then
                  in_play(k) = .FALSE.
                else
                  in_play(k) = .TRUE.
                end if

                CYCLE

              end if

              if ( CL_sum(i,j,k) < VERY_SMALL ) then

                !! We should not divide  BRCL_diff  by  CL_sum
                !! Consider all of the amount to be extra   [kmol_CL / kmol_AIR]
                extra = BRCL_diff(i,j,k)

                !! Convert to   [kmol_CL]
                extra = extra * kmol_per_kg_AIR(i,j,k) * mass_AIR(i,j,k)

                vsum_EXTRA = vsum_EXTRA + extra

                in_play(k) = .FALSE.

                CYCLE

              end if

              !! At this point, BRCL_diff > 0, and CL_sum is not too small

              f =   BRCL_diff(i,j,k) / CL_sum(i,j,k)  ! Note frac is positive
                                                      ! Multiply by  (1 - frac)

              !! In gridboxes where (1-frac) gets too small or negative,
              !! accumulate the extra mass of the species
              if ( f >  MAXFRAC ) then

                !! First compute total in gridbox    [kmol_CL / kmol_AIR]
                extra = CL_sum(i,j,k)

                !! Convert to remaining fraction  [kmol_CL]
                extra = extra * ( f - MAXFRAC ) * kmol_per_kg_AIR(i,j,k) * mass_AIR(i,j,k)

                vsum_EXTRA = vsum_EXTRA + extra

                !! Apply f=MAXFRAC
                x(iCL2  )%p(i,j,k) = x(iCL2  )%p(i,j,k) * ( 1.0 - MAXFRAC )
                x(iOCLO )%p(i,j,k) = x(iOCLO )%p(i,j,k) * ( 1.0 - MAXFRAC )
                x(iCL2O2)%p(i,j,k) = x(iCL2O2)%p(i,j,k) * ( 1.0 - MAXFRAC )
                x(iCLO  )%p(i,j,k) = x(iCLO  )%p(i,j,k) * ( 1.0 - MAXFRAC )
                x(iHCL  )%p(i,j,k) = x(iHCL  )%p(i,j,k) * ( 1.0 - MAXFRAC )
                x(iHOCL )%p(i,j,k) = x(iHOCL )%p(i,j,k) * ( 1.0 - MAXFRAC )

                in_play(k) = .FALSE.

              else   ! ( f > 0.0  .AND.  f <= MAXFRAC )

                x(iCL2  )%p(i,j,k) = x(iCL2  )%p(i,j,k) * ( 1.0 - f )
                x(iOCLO )%p(i,j,k) = x(iOCLO )%p(i,j,k) * ( 1.0 - f )
                x(iCL2O2)%p(i,j,k) = x(iCL2O2)%p(i,j,k) * ( 1.0 - f )
                x(iCLO  )%p(i,j,k) = x(iCLO  )%p(i,j,k) * ( 1.0 - f )
                x(iHCL  )%p(i,j,k) = x(iHCL  )%p(i,j,k) * ( 1.0 - f )
                x(iHOCL )%p(i,j,k) = x(iHOCL )%p(i,j,k) * ( 1.0 - f )

                !! If fraction is small we will allow further reduction
                !! in the gridbox; we expect that to happen at most 2 or 3 times
                if (f < 0.70) then
                  in_play(k) = .TRUE.
                else
                  in_play(k) = .FALSE.
                end if

              end if

            end do


            !! Distribute any extra mass of Cly species throughout the entire column

            loop_count = 0

            !!
            !! Additional passes through the column, in which
            !! the remainder amount is divided among
            !! all column gridboxes that are still 'in play'
            !!

            DO WHILE ( vsum_EXTRA > 0.0d0 .AND. ANY(in_play(:)) )

              !! Mixing ratio to be subtracted from every box in
              !! the column that is still in play    [kmol_CL / kg_AIR]
              col_mixrat = vsum_EXTRA / SUM(mass_AIR(i,j,:), MASK=in_play(:) )

              vsum_EXTRA = 0.0d0

              do k = k1,k2

                if ( .NOT. in_play(k) ) CYCLE

                !!
                !! Compute the fraction
                !!

                !! First compute total in gridbox   [kmol_CL / kmol_AIR]
                box_mixrat = 2 * x(iCL2  )%p(i,j,k) +  &
                                 x(iOCLO )%p(i,j,k) +  &
                             2 * x(iCL2O2)%p(i,j,k) +  &
                                 x(iCLO  )%p(i,j,k) +  &
                                 x(iHCL  )%p(i,j,k) +  &
                                 x(iHOCL )%p(i,j,k)

                !! Convert to [kmol_CL / kg_AIR]
                box_mixrat = box_mixrat * kmol_per_kg_AIR(i,j,k)

                f = col_mixrat/box_mixrat

                !!
                !! Apply the fraction
                !!

                if ( f > MAXFRAC ) then

                  !! First compute total in gridbox   [kmol_CL / kmol_AIR]
                  extra = 2 * x(iCL2  )%p(i,j,k) +  &
                              x(iOCLO )%p(i,j,k) +  &
                          2 * x(iCL2O2)%p(i,j,k) +  &
                              x(iCLO  )%p(i,j,k) +  &
                              x(iHCL  )%p(i,j,k) +  &
                              x(iHOCL )%p(i,j,k)

                  !! Convert to remaining fraction  [kmol_CL]
                  extra = extra * ( f - MAXFRAC ) * kmol_per_kg_AIR(i,j,k) * mass_AIR(i,j,k)

                  vsum_EXTRA = vsum_EXTRA + extra

                  !! Apply f=MAXFRAC
                  x(iCL2  )%p(i,j,k) = x(iCL2  )%p(i,j,k) * ( 1.0 - MAXFRAC )
                  x(iOCLO )%p(i,j,k) = x(iOCLO )%p(i,j,k) * ( 1.0 - MAXFRAC )
                  x(iCL2O2)%p(i,j,k) = x(iCL2O2)%p(i,j,k) * ( 1.0 - MAXFRAC )
                  x(iCLO  )%p(i,j,k) = x(iCLO  )%p(i,j,k) * ( 1.0 - MAXFRAC )
                  x(iHCL  )%p(i,j,k) = x(iHCL  )%p(i,j,k) * ( 1.0 - MAXFRAC )
                  x(iHOCL )%p(i,j,k) = x(iHOCL )%p(i,j,k) * ( 1.0 - MAXFRAC )

                  in_play(k) = .FALSE.

                else

                  x(iCL2  )%p(i,j,k) = x(iCL2  )%p(i,j,k) * ( 1.0 - f )
                  x(iOCLO )%p(i,j,k) = x(iOCLO )%p(i,j,k) * ( 1.0 - f )
                  x(iCL2O2)%p(i,j,k) = x(iCL2O2)%p(i,j,k) * ( 1.0 - f )
                  x(iCLO  )%p(i,j,k) = x(iCLO  )%p(i,j,k) * ( 1.0 - f )
                  x(iHCL  )%p(i,j,k) = x(iHCL  )%p(i,j,k) * ( 1.0 - f )
                  x(iHOCL )%p(i,j,k) = x(iHOCL )%p(i,j,k) * ( 1.0 - f )

                  !! If fraction is small we will allow further reduction
                  !! in the gridbox; we expect that to happen at most 2 or 3 times
                  if (f < 0.70) then
                    in_play(k) = .TRUE.
                  else
                    in_play(k) = .FALSE.
                  end if

                end if

              end do

              loop_count = loop_count+1

            END DO

            IF ( loop_count > 0 ) print*,'CL LOOP COUNT = ', loop_count, COUNT(in_play(:))

            if ( vsum_EXTRA > 0.0d0 ) then
              print*,'Cannot conserve Cl !!'
            end if

          end do
        end do

    END SUBROUTINE GMI_reduce_Cly_vmr

!.... In gridboxes where ClONO2 has increased, reduce N species
    SUBROUTINE GMI_reduce_N_vmr ( )

       integer :: loop_count, i, j, k
       real*8  :: f, extra, vsum_EXTRA, col_mixrat, box_mixrat

#ifdef DEBUG_STUFF
       real :: the_min
       integer:: ii,jj,kk

       the_min = 1000.0
#endif

        do i = i1,i2
          do j = j1,j2

            vsum_EXTRA = 0.0d0  !  extra within column  [kmol_N]

            !! Guard against (1-frac) getting too small, or going negative

            !!
            !! First pass through the column; make sure that in_play(k) is initialized for all k
            !!
            do k = k1,k2

#ifdef DEBUG_STUFF
       IF ( N_sum(i,j,k) < the_min ) THEN
         the_min = N_sum(i,j,k)
         ii = i
         jj = j
         kk = k
       ENDIF
#endif

              if ( CLONO2_diff(i,j,k) <= 0.0 ) then

                !! There was addition (not reduction) in this gridbox;
                !! permit reduction in the next phase (column work)
                in_play(k) = .TRUE.

                CYCLE

              end if

              if ( N_sum(i,j,k) < VERY_SMALL ) then

                !! We should not divide  CLONO2_diff  by  N_sum
                !! Consider all of the amount to be extra   [kmol_N / kmol_AIR]
                extra = CLONO2_diff(i,j,k)

                !! Convert to   [kmol_N]
                extra = extra * kmol_per_kg_AIR(i,j,k) * mass_AIR(i,j,k)

                vsum_EXTRA = vsum_EXTRA + extra

                in_play(k) = .FALSE.

                CYCLE

              end if

              f =   CLONO2_diff(i,j,k) / N_sum(i,j,k)  ! Note frac is positive
                                                       ! Multiply by  (1 - frac)

              !! In gridboxes where (1-frac) gets too small or negative,
              !! accumulate the extra mass of the species
              if ( f >  MAXFRAC ) then

                !! First compute total in gridbox    [kmol_N / kmol_AIR]
                extra = N_sum(i,j,k)

                !! Convert to remaining fraction  [kmol_N]
                extra = extra * ( f - MAXFRAC ) * kmol_per_kg_AIR(i,j,k) * mass_AIR(i,j,k)

                vsum_EXTRA = vsum_EXTRA + extra

                !! Apply f=MAXFRAC
                x(iN2O5)%p(i,j,k) = x(iN2O5)%p(i,j,k) * ( 1.0 - MAXFRAC )
                x(iNO  )%p(i,j,k) = x(iNO  )%p(i,j,k) * ( 1.0 - MAXFRAC )
                x(iNO2 )%p(i,j,k) = x(iNO2 )%p(i,j,k) * ( 1.0 - MAXFRAC )
                x(iHNO3)%p(i,j,k) = x(iHNO3)%p(i,j,k) * ( 1.0 - MAXFRAC )

                in_play(k) = .FALSE.

              else   ! ( f > 0.0  .AND.  f <= MAXFRAC )

                x(iN2O5)%p(i,j,k) = x(iN2O5)%p(i,j,k) * ( 1.0 - f )
                x(iNO  )%p(i,j,k) = x(iNO  )%p(i,j,k) * ( 1.0 - f )
                x(iNO2 )%p(i,j,k) = x(iNO2 )%p(i,j,k) * ( 1.0 - f )
                x(iHNO3)%p(i,j,k) = x(iHNO3)%p(i,j,k) * ( 1.0 - f )

                !! If fraction is small we will allow further reduction
                !! in the gridbox; we expect that to happen at most 2 or 3 times
                if (f < 0.70) then
                  in_play(k) = .TRUE.
                else
                  in_play(k) = .FALSE.
                end if

              end if

            end do


            !! Distribute any extra mass of N species throughout the entire column

            loop_count = 0

            !!
            !! Additional passes through the column, in which
            !! the remainder amount is divided among
            !! all column gridboxes that are still 'in play'
            !!

            DO WHILE ( vsum_EXTRA > 0.0d0 .AND. ANY(in_play(:)) )

              !! Mixing ratio to be subtracted from every box in
              !! the column that is still in play    [kmol_N / kg_AIR]
              col_mixrat = vsum_EXTRA / SUM(mass_AIR(i,j,:), MASK=in_play(:) )

              vsum_EXTRA = 0.0d0

              do k = k1,k2

                if ( .NOT. in_play(k) ) CYCLE

                !!
                !! Compute the fraction
                !!

                !! First compute total in gridbox    [kmol_N / kmol_AIR]
                box_mixrat = 2 * x(iN2O5)%p(i,j,k) +  &
                                 x(iNO  )%p(i,j,k) +  &
                                 x(iNO2 )%p(i,j,k) +  &
                                 x(iHNO3)%p(i,j,k)

                !! Convert to [kmol_N / kg_AIR]
                box_mixrat = box_mixrat * kmol_per_kg_AIR(i,j,k)

                f = col_mixrat/box_mixrat

                !!
                !! Apply the fraction
                !!

                if ( f > MAXFRAC ) then

                  !! First compute total in gridbox    [kmol_N / kmol_AIR]
                  extra = 2 * x(iN2O5)%p(i,j,k) / mw_N2O5 +  &
                              x(iNO  )%p(i,j,k) / mw_NO   +  &
                              x(iNO2 )%p(i,j,k) / mw_NO2  +  &
                              x(iHNO3)%p(i,j,k) / mw_HNO3

                  !! Convert to remaining fraction  [kmol_N]
                  extra = extra * ( f - MAXFRAC ) * kmol_per_kg_AIR(i,j,k) * mass_AIR(i,j,k)

                  vsum_EXTRA = vsum_EXTRA + extra

                  !! Apply f=MAXFRAC
                  x(iN2O5)%p(i,j,k) = x(iN2O5)%p(i,j,k) * ( 1.0 - MAXFRAC )
                  x(iNO  )%p(i,j,k) = x(iNO  )%p(i,j,k) * ( 1.0 - MAXFRAC )
                  x(iNO2 )%p(i,j,k) = x(iNO2 )%p(i,j,k) * ( 1.0 - MAXFRAC )
                  x(iHNO3)%p(i,j,k) = x(iHNO3)%p(i,j,k) * ( 1.0 - MAXFRAC )

                  in_play(k) = .FALSE.

                else

                  x(iN2O5)%p(i,j,k) = x(iN2O5)%p(i,j,k) * ( 1.0 - f )
                  x(iNO  )%p(i,j,k) = x(iNO  )%p(i,j,k) * ( 1.0 - f )
                  x(iNO2 )%p(i,j,k) = x(iNO2 )%p(i,j,k) * ( 1.0 - f )
                  x(iHNO3)%p(i,j,k) = x(iHNO3)%p(i,j,k) * ( 1.0 - f )

                  !! If fraction is small we will allow further reduction
                  !! in the gridbox; we expect that to happen at most 2 or 3 times
                  if (f < 0.70) then
                    in_play(k) = .TRUE.
                  else
                    in_play(k) = .FALSE.
                  end if

                end if

              end do

              loop_count = loop_count+1

            END DO

            IF ( loop_count > 0 ) print*,'N LOOP COUNT = ', loop_count, COUNT(in_play(:))

            if ( vsum_EXTRA > 0.0d0 ) then
              print*,'Cannot conserve N!!'
            end if

          end do
        end do

#ifdef DEBUG_STUFF
        IF ( the_min < 1.0e-15 ) print*,'N_sum min = ',the_min, ' ijk ', ii, jj, kk
#endif

    END SUBROUTINE GMI_reduce_N_vmr

!.... In gridboxes where BrCl has increased, reduce Cly species
    SUBROUTINE GCC_reduce_Cly_mmr ( )

       integer :: loop_count, i, j, k
       real*8  :: f, extra, vsum_EXTRA, col_mixrat, box_mixrat

        do i = i1,i2
          do j = j1,j2

            vsum_EXTRA = 0.0d0  !  extra within column  [kmol_CL]

            !! Guard against (1-frac) getting too small, or going negative

            !!
            !! First pass through the column; make sure that in_play(k) is initialized for all k
            !!
            do k = k1,k2

              if ( BRCL_diff(i,j,k) <= 0.0 ) then

                !! There was addition (not reduction) in this gridbox;
                !! permit reduction in the next phase (column work)
                !! if there is enough Cly
                if ( CL_sum(i,j,k) < VERY_SMALL ) then
                  in_play(k) = .FALSE.
                else
                  in_play(k) = .TRUE.
                end if

                CYCLE

              end if

              if ( CL_sum(i,j,k) < VERY_SMALL ) then

                !! We should not divide  BRCL_diff  by  CL_sum
                !! Consider all of the amount to be extra   [kmol_CL / kg_AIR]
                extra = BRCL_diff(i,j,k)

                !! Convert to   [kmol_CL]
                extra = extra * mass_AIR(i,j,k)

                vsum_EXTRA = vsum_EXTRA + extra

                in_play(k) = .FALSE.

                CYCLE

              end if

              !! At this point, BRCL_diff > 0, and CL_sum is not too small

              f =   BRCL_diff(i,j,k) / CL_sum(i,j,k)  ! Note frac is positive
                                                      ! Multiply by  (1 - frac)

              !! In gridboxes where (1-frac) gets too small or negative,
              !! accumulate the extra mass of the species
              if ( f >  MAXFRAC ) then

                !! First compute total in gridbox    [kmol_CL / kg_AIR]
                extra = CL_sum(i,j,k)

                !! Convert to remaining fraction  [kmol_CL]
                extra = extra * ( f - MAXFRAC ) * mass_AIR(i,j,k)

                vsum_EXTRA = vsum_EXTRA + extra

                !! Apply f=MAXFRAC
                y(iCL2  )%p(i,j,k) = y(iCL2  )%p(i,j,k) * ( 1.0 - MAXFRAC )
                y(iOCLO )%p(i,j,k) = y(iOCLO )%p(i,j,k) * ( 1.0 - MAXFRAC )
                y(iCL2O2)%p(i,j,k) = y(iCL2O2)%p(i,j,k) * ( 1.0 - MAXFRAC )
                y(iCLO  )%p(i,j,k) = y(iCLO  )%p(i,j,k) * ( 1.0 - MAXFRAC )
                y(iHCL  )%p(i,j,k) = y(iHCL  )%p(i,j,k) * ( 1.0 - MAXFRAC )
                y(iHOCL )%p(i,j,k) = y(iHOCL )%p(i,j,k) * ( 1.0 - MAXFRAC )

                in_play(k) = .FALSE.

              else   ! ( f > 0.0  .AND.  f <= MAXFRAC )

                y(iCL2  )%p(i,j,k) = y(iCL2  )%p(i,j,k) * ( 1.0 - f )
                y(iOCLO )%p(i,j,k) = y(iOCLO )%p(i,j,k) * ( 1.0 - f )
                y(iCL2O2)%p(i,j,k) = y(iCL2O2)%p(i,j,k) * ( 1.0 - f )
                y(iCLO  )%p(i,j,k) = y(iCLO  )%p(i,j,k) * ( 1.0 - f )
                y(iHCL  )%p(i,j,k) = y(iHCL  )%p(i,j,k) * ( 1.0 - f )
                y(iHOCL )%p(i,j,k) = y(iHOCL )%p(i,j,k) * ( 1.0 - f )

                !! If fraction is small we will allow further reduction
                !! in the gridbox; we expect that to happen at most 2 or 3 times
                if (f < 0.70) then
                  in_play(k) = .TRUE.
                else
                  in_play(k) = .FALSE.
                end if

              end if

            end do


            !! Distribute any extra mass of Cly species throughout the entire column

#ifdef CHEM_INFO
      ! aCLY and zCLY will show the change IN THE COLUMN, due to BrCl adjustment
      IF ( ASSOCIATED(aCLY) ) THEN

        qq2(i,j,:) = 0.0
        ig = 2

        do im                      = 1, MAXGRP_ELEM
          imsgrp                   = sgrp_elem_map(im,ig)
          if (imsgrp > 0)  &
     &      qq2(i,j,:)             = qq2(i,j,:) +  &
     &                               y(imsgrp)%p(i,j,:) *  &
     &                               sgrp_fac(im,ig)
        end do

        aCLY(i,j,:) = qq2(i,j,:)

      ENDIF
#endif

            loop_count = 0

            !!
            !! Additional passes through the column, in which
            !! the remainder amount is divided among
            !! all column gridboxes that are still 'in play'
            !!

            DO WHILE ( vsum_EXTRA > 0.0d0 .AND. ANY(in_play(:)) )

              !! Mixing ratio to be subtracted from every box in
              !! the column that is still in play    [kmol_CL / kg_AIR]
              col_mixrat = vsum_EXTRA / SUM(mass_AIR(i,j,:), MASK=in_play(:) )

              vsum_EXTRA = 0.0d0

              do k = k1,k2

                if ( .NOT. in_play(k) ) CYCLE

                !!
                !! Compute the fraction
                !!

                !! First compute total in gridbox   [kmol_CL / kg_AIR]
                box_mixrat = 2 * y(iCL2  )%p(i,j,k) / mw_CL2    +  &
                                 y(iOCLO )%p(i,j,k) / mw_OCLO   +  &
                             2 * y(iCL2O2)%p(i,j,k) / mw_CL2O2  +  &
                                 y(iCLO  )%p(i,j,k) / mw_CLO    +  &
                                 y(iHCL  )%p(i,j,k) / mw_HCL    +  &
                                 y(iHOCL )%p(i,j,k) / mw_HOCL

                f = col_mixrat/box_mixrat

                !!
                !! Apply the fraction
                !!

                if ( f > MAXFRAC ) then

                  !! First compute total in gridbox   [kmol_CL / kg_AIR]
                  extra = 2 * y(iCL2  )%p(i,j,k) / mw_CL2    +  &
                              y(iOCLO )%p(i,j,k) / mw_OCLO   +  &
                          2 * y(iCL2O2)%p(i,j,k) / mw_CL2O2  +  &
                              y(iCLO  )%p(i,j,k) / mw_CLO    +  &
                              y(iHCL  )%p(i,j,k) / mw_HCL    +  &
                              y(iHOCL )%p(i,j,k) / mw_HOCL

                  !! Convert to remaining fraction  [kmol_CL]
                  extra = extra * ( f - MAXFRAC ) * mass_AIR(i,j,k)

                  vsum_EXTRA = vsum_EXTRA + extra

                  !! Apply f=MAXFRAC
                  y(iCL2  )%p(i,j,k) = y(iCL2  )%p(i,j,k) * ( 1.0 - MAXFRAC )
                  y(iOCLO )%p(i,j,k) = y(iOCLO )%p(i,j,k) * ( 1.0 - MAXFRAC )
                  y(iCL2O2)%p(i,j,k) = y(iCL2O2)%p(i,j,k) * ( 1.0 - MAXFRAC )
                  y(iCLO  )%p(i,j,k) = y(iCLO  )%p(i,j,k) * ( 1.0 - MAXFRAC )
                  y(iHCL  )%p(i,j,k) = y(iHCL  )%p(i,j,k) * ( 1.0 - MAXFRAC )
                  y(iHOCL )%p(i,j,k) = y(iHOCL )%p(i,j,k) * ( 1.0 - MAXFRAC )

                  in_play(k) = .FALSE.

                else

                  y(iCL2  )%p(i,j,k) = y(iCL2  )%p(i,j,k) * ( 1.0 - f )
                  y(iOCLO )%p(i,j,k) = y(iOCLO )%p(i,j,k) * ( 1.0 - f )
                  y(iCL2O2)%p(i,j,k) = y(iCL2O2)%p(i,j,k) * ( 1.0 - f )
                  y(iCLO  )%p(i,j,k) = y(iCLO  )%p(i,j,k) * ( 1.0 - f )
                  y(iHCL  )%p(i,j,k) = y(iHCL  )%p(i,j,k) * ( 1.0 - f )
                  y(iHOCL )%p(i,j,k) = y(iHOCL )%p(i,j,k) * ( 1.0 - f )

                  !! If fraction is small we will allow further reduction
                  !! in the gridbox; we expect that to happen at most 2 or 3 times
                  if (f < 0.70) then
                    in_play(k) = .TRUE.
                  else
                    in_play(k) = .FALSE.
                  end if

                end if

              end do

              loop_count = loop_count+1

            END DO

            IF ( loop_count > 0 ) print*,'CL LOOP COUNT = ', loop_count, COUNT(in_play(:))

            if ( vsum_EXTRA > 0.0d0 ) then
              print*,'Cannot conserve Cl !!'
            end if

          end do
        end do

    END SUBROUTINE GCC_reduce_Cly_mmr

!.... In gridboxes where ClNOX has increased, reduce N species
    SUBROUTINE GCC_reduce_N_mmr ( )

       integer :: loop_count, i, j, k
       real*8  :: f, extra, vsum_EXTRA, col_mixrat, box_mixrat

        do i = i1,i2
          do j = j1,j2

            vsum_EXTRA = 0.0d0  !  extra within column  [kmol_N]

            !! Guard against (1-frac) getting too small, or going negative

            !!
            !! First pass through the column; make sure that in_play(k) is initialized for all k
            !!
            do k = k1,k2

              if ( CLNOX_diff(i,j,k) <= 0.0 ) then

                !! There was addition (not reduction) in this gridbox;
                !! permit reduction in the next phase (column work)
                in_play(k) = .TRUE.

                CYCLE

              end if

              if ( N_sum(i,j,k) < VERY_SMALL ) then

                !! We should not divide  CLNOX_diff  by  N_sum
                !! Consider all of the amount to be extra   [kmol_N / kg_AIR]
                extra = CLNOX_diff(i,j,k)

                !! Convert to   [kmol_N]
                extra = extra * mass_AIR(i,j,k)

                vsum_EXTRA = vsum_EXTRA + extra

                in_play(k) = .FALSE.

                CYCLE

              end if

              f =   CLNOX_diff(i,j,k) / N_sum(i,j,k)  ! Note frac is positive
                                                      ! Multiply by  (1 - frac)

              !! In gridboxes where (1-frac) gets too small or negative,
              !! accumulate the extra mass of the species
              if ( f >  MAXFRAC ) then

                !! First compute total in gridbox    [kmol_N / kg_AIR]
                extra = N_sum(i,j,k)

                !! Convert to remaining fraction  [kmol_N]
                extra = extra * ( f - MAXFRAC ) * mass_AIR(i,j,k)

                vsum_EXTRA = vsum_EXTRA + extra

                !! Apply f=MAXFRAC
                y(iN2O5)%p(i,j,k) = y(iN2O5)%p(i,j,k) * ( 1.0 - MAXFRAC )
                y(iNO  )%p(i,j,k) = y(iNO  )%p(i,j,k) * ( 1.0 - MAXFRAC )
                y(iNO2 )%p(i,j,k) = y(iNO2 )%p(i,j,k) * ( 1.0 - MAXFRAC )
                y(iHNO3)%p(i,j,k) = y(iHNO3)%p(i,j,k) * ( 1.0 - MAXFRAC )

                in_play(k) = .FALSE.

              else   ! ( f > 0.0  .AND.  f <= MAXFRAC )

                y(iN2O5)%p(i,j,k) = y(iN2O5)%p(i,j,k) * ( 1.0 - f )
                y(iNO  )%p(i,j,k) = y(iNO  )%p(i,j,k) * ( 1.0 - f )
                y(iNO2 )%p(i,j,k) = y(iNO2 )%p(i,j,k) * ( 1.0 - f )
                y(iHNO3)%p(i,j,k) = y(iHNO3)%p(i,j,k) * ( 1.0 - f )

                !! If fraction is small we will allow further reduction
                !! in the gridbox; we expect that to happen at most 2 or 3 times
                if (f < 0.70) then
                  in_play(k) = .TRUE.
                else
                  in_play(k) = .FALSE.
                end if

              end if

            end do


            !! Distribute any extra mass of N species throughout the entire column

            loop_count = 0

            !!
            !! Additional passes through the column, in which
            !! the remainder amount is divided among
            !! all column gridboxes that are still 'in play'
            !!

            DO WHILE ( vsum_EXTRA > 0.0d0 .AND. ANY(in_play(:)) )

              !! Mixing ratio to be subtracted from every box in
              !! the column that is still in play    [kmol_N / kg_AIR]
              col_mixrat = vsum_EXTRA / SUM(mass_AIR(i,j,:), MASK=in_play(:) )

              vsum_EXTRA = 0.0d0

              do k = k1,k2

                if ( .NOT. in_play(k) ) CYCLE

                !!
                !! Compute the fraction
                !!

                !! First compute total in gridbox    [kmol_N / kg_AIR]
                box_mixrat = 2 * y(iN2O5)%p(i,j,k) / mw_N2O5 +  &
                                 y(iNO  )%p(i,j,k) / mw_NO   +  &
                                 y(iNO2 )%p(i,j,k) / mw_NO2  +  &
                                 y(iHNO3)%p(i,j,k) / mw_HNO3

                f = col_mixrat/box_mixrat

                !!
                !! Apply the fraction
                !!

                if ( f > MAXFRAC ) then

                  !! First compute total in gridbox    [kmol_N / kg_AIR]
                  extra = 2 * y(iN2O5)%p(i,j,k) / mw_N2O5 +  &
                              y(iNO  )%p(i,j,k) / mw_NO   +  &
                              y(iNO2 )%p(i,j,k) / mw_NO2  +  &
                              y(iHNO3)%p(i,j,k) / mw_HNO3

                  !! Convert to remaining fraction  [kmol_N]
                  extra = extra * ( f - MAXFRAC ) * mass_AIR(i,j,k)

                  vsum_EXTRA = vsum_EXTRA + extra

                  !! Apply f=MAXFRAC
                  y(iN2O5)%p(i,j,k) = y(iN2O5)%p(i,j,k) * ( 1.0 - MAXFRAC )
                  y(iNO  )%p(i,j,k) = y(iNO  )%p(i,j,k) * ( 1.0 - MAXFRAC )
                  y(iNO2 )%p(i,j,k) = y(iNO2 )%p(i,j,k) * ( 1.0 - MAXFRAC )
                  y(iHNO3)%p(i,j,k) = y(iHNO3)%p(i,j,k) * ( 1.0 - MAXFRAC )

                  in_play(k) = .FALSE.

                else

                  y(iN2O5)%p(i,j,k) = y(iN2O5)%p(i,j,k) * ( 1.0 - f )
                  y(iNO  )%p(i,j,k) = y(iNO  )%p(i,j,k) * ( 1.0 - f )
                  y(iNO2 )%p(i,j,k) = y(iNO2 )%p(i,j,k) * ( 1.0 - f )
                  y(iHNO3)%p(i,j,k) = y(iHNO3)%p(i,j,k) * ( 1.0 - f )

                  !! If fraction is small we will allow further reduction
                  !! in the gridbox; we expect that to happen at most 2 or 3 times
                  if (f < 0.70) then
                    in_play(k) = .TRUE.
                  else
                    in_play(k) = .FALSE.
                  end if

                end if

              end do

              loop_count = loop_count+1

            END DO

            IF ( loop_count > 0 ) print*,'N LOOP COUNT = ', loop_count, COUNT(in_play(:))

            if ( vsum_EXTRA > 0.0d0 ) then
              print*,'Cannot conserve N!!'
            end if

          end do
        end do

    END SUBROUTINE GCC_reduce_N_mmr

 END SUBROUTINE Unpack_Chem_Groups

 END MODULE Chem_GroupMod
