   MODULE  cblock_CO_budget_mod

!
! !DESCRIPTION:
!
!  This module replaces CMN_CO_BUDGET, storing common blocks and parameters
!
! !REVISION HISTORY:
!
!  13October2015 Manyin  First crack.
!

   IMPLICIT NONE

!
!*****************************************************************************
!  CMN_CO_BUDGET contains variables for computing CO budgets in the
!  CO simulation with parameterized OH. (bnd, bmy, 4/19/00)
!*****************************************************************************
!
      ! NTALLY is the number of budget types for the TCO array
      INTEGER, PARAMETER :: NTALLY=12

      ! TCO is the diagnostic array for cross tropopause fluxes
!KNJR      REAL*8             :: TCO(IIPAR,JJPAR,LLPAR,NTALLY)
!KNJR      COMMON /COb55/        TCO

      ! FMOL_CO = the molecular weight of CO in kg
      REAL*8, PARAMETER  :: FMOL_CO = 28d-3

      ! XNUMOL_CO = Molecules CO / kg CO
      REAL*8, PARAMETER  :: XNUMOL_CO = 6.022D+23 / FMOL_CO

      ! TCVV_CO = 28.97 [g/mole air] / 28.0 [g/mole CO]
      REAL*8, PARAMETER  :: TCVV_CO = 28.97d0 / 28.0d0


!-------------------------------------------------------------------------
CONTAINS
!------------------------------------------------------------------------------
END MODULE  cblock_CO_budget_mod
