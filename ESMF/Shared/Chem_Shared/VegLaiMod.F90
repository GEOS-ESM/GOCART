#include "MAPL_Generic.h"

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  VegLaiMod --- Vegetation Index and Leaf Area Index
!
! !INTERFACE:
!

   MODULE  VegLaiMod

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

   PUBLIC  Decode_Land_Types
   PUBLIC  Decode_XLAI

!
! !DESCRIPTION:
!
!  This module decodes land-use and LAI info from 3D arrays.
!  The Olson land types are used.
!
! !REVISION HISTORY:
!
!  16May2016 - Manyin, first crack
!
!EOP

!-------------------------------------------------------------------------
CONTAINS

!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS !
!            Atmospheric Chemistry and Dynamics Lab,  Code 614             !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  Decode_Land_Types - Decode Olson Land Types from a 3D array
!
! !INTERFACE:

 SUBROUTINE Decode_Land_Types(PTR3D, NTYPE, i1, i2, j1, j2, km, ireg, iuse, iland, rc)

  IMPLICIT NONE

! !INPUT PARAMETERS:
  REAL,    POINTER, INTENT(IN),    DIMENSION(:,:,:) :: PTR3D
  INTEGER,          INTENT(IN)                      :: NTYPE      ! max number of land types in a gridbox
  INTEGER,          INTENT(IN)                      :: i1, i2, j1, j2, km

! !OUTPUT PARAMETERS
  INTEGER, POINTER, INTENT(INOUT), DIMENSION(:,:)   :: ireg       ! number of land types in each grid square
  INTEGER, POINTER, INTENT(INOUT), DIMENSION(:,:,:) :: iuse       ! fraction of grid box area occupied by land type
  INTEGER, POINTER, INTENT(INOUT), DIMENSION(:,:,:) :: iland      ! land type id in grid square for ireg land types

  INTEGER,          INTENT(OUT) :: rc
!
! !DESCRIPTION:
!
!  This routine ingests land use mil fractions from Michael Manyin's NetCDF files, which are 
!  set up to facilitate mapping to the cubed sphere. The incoming 3-D field actually consists of
!  72 2-D fractional coverage fields, one for each Olson land type. (There are 74 land types, but
!  a few such as 13 and 14 are unused, so in the file we store the coverage for type 73 in 13, and
!  the coverage for type 74 in 14.)  Fractional coverage is a value between 0 and 1000.
!  We impose the constraint that, for any gridbox, there can be at most NTYPE land types stored.
!  So we test each coverage fraction to be sure it satisfies the threshold VEG_MIN_VALID, and
!  we only track info for those entries. If the number of qualifying entries in a gridbox
!  exceeds NTPYE, the mil fractions are sorted in ascending order by MAPL_Sort.
!  The "upper" NTYPE values are then copied into the output arrays. If the number of
!  qualifying land types is less than or equal to NTYPE, there is no need to sort, and all the
!  values are copied. 
!
!  NOTE: The fraction coverage fields are constant in the current Olson scheme.
!
!EOP
!---------------------------------------------------------------------------
  REAL, PARAMETER :: VEG_MIN_VALID = 1.0

  CHARACTER(LEN=ESMF_MAXSTR) :: IAm
  INTEGER :: STATUS
  INTEGER :: i,j,k,m,ic
  INTEGER, ALLOCATABLE :: landNum(:)
  INTEGER, ALLOCATABLE :: milFrac(:)

  rc = 0
  IAm = "Decode_Land_Types"

    ireg(:,:)   = 0
    iuse(:,:,:) = 0
   iland(:,:,:) = 0

   ! Allocate the max number of possible entries
   ALLOCATE( landNum(km), __STAT__ )
   ALLOCATE( milFrac(km), __STAT__ )

   DO j = j1, j2
    DO i = i1, i2

     ic = 0
     landNum(:) = -1  ! in case we need to sort (ic > NTYPE)
     milFrac(:) = -1  ! in case we need to sort (ic > NTYPE)
     DO k = 1, km

      IF(PTR3D(i,j,k) >= VEG_MIN_VALID) THEN
       ic = ic+1
       m = k
       IF(k == 13) m = 73
       IF(k == 14) m = 74
! Expecting land type number 0-73 but read in as 1-74 
!      landNum(ic) = m
       landNum(ic) = m - 1
       milFrac(ic) = INT(PTR3D(i,j,k)+0.0001)
      END IF

     END DO

     IF(ic > NTYPE) THEN

      CALL MAPL_Sort(milFrac(1:ic),landNum(1:ic))

      iland(i,j,1:NTYPE) = landNum(ic:ic-(NTYPE-1):-1)
       iuse(i,j,1:NTYPE) = milFrac(ic:ic-(NTYPE-1):-1)
       ireg(i,j)         = NTYPE

     ELSE

      iland(i,j,1:ic) = landNum(1:ic)
       iuse(i,j,1:ic) = milFrac(1:ic)
       ireg(i,j)      = ic

     END IF

    END DO
   END DO

   DEALLOCATE(landNum, __STAT__ )
   DEALLOCATE(milFrac, __STAT__ )

  RETURN
 END SUBROUTINE Decode_Land_Types


!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS !
!            Atmospheric Chemistry and Dynamics Lab,  Code 614             !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  Decode_XLAI - Decode Leaf Area Indices from a 3D array
!
! !INTERFACE:

 SUBROUTINE Decode_XLAI(PTR3D, NTYPE, i1, i2, j1, j2, km, ireg, iuse, iland, xlai, rc)

  IMPLICIT NONE

! !INPUT PARAMETERS:
  REAL,    POINTER, INTENT(IN),    DIMENSION(:,:,:) :: PTR3D
  INTEGER,          INTENT(IN)                      :: NTYPE
  INTEGER,          INTENT(IN)                      :: i1, i2, j1, j2, km
  INTEGER, POINTER, INTENT(IN),    DIMENSION(:,:)   :: ireg       ! number of land types in a grid square
  INTEGER, POINTER, INTENT(IN),    DIMENSION(:,:,:) :: iuse       ! fraction of grid box area occupied by land type
  INTEGER, POINTER, INTENT(IN),    DIMENSION(:,:,:) :: iland      ! land type id in grid square for ireg land types

! !OUTPUT PARAMETERS
  REAL*8 , POINTER, INTENT(INOUT), DIMENSION(:,:,:) :: xlai       ! leaf area index of land type

  INTEGER,          INTENT(OUT) :: rc
!
! !DESCRIPTION:
!
!  This routine ingests Leaf Area Index fractions from Michael Manyin's NetCDF files, which are 
!  set up to facilitate mapping to the cubed sphere. The incoming 3-D field actually consists of
!  72  2-D  LAI fields, one for each Olson land type. (There are 74 land types, but
!  a few such as 13 and 14 are unused, so in the file we store the LAI for type 73 in 13, and
!  the coverage for type 74 in 14.)  LAI is a floating point fraction.  It varies monthly.
!
!EOP
!---------------------------------------------------------------------------
  CHARACTER(LEN=ESMF_MAXSTR) :: IAm
  INTEGER :: STATUS
  INTEGER :: i,j,k,ic

  rc = 0
  IAm = "Decode_XLAI"

   DO j = j1,j2
    DO i = i1,i2

     DO ic = 1,ireg(i,j)
! Add back 1 for correct k indices 
!     k = iland(i,j,ic)
      k = iland(i,j,ic) + 1
      IF(k > 72) k = k-60  ! There are 74 land types, but some are unused so 73->13, 74->14
      IF(k < 1 .OR. k > 72) THEN
        print*,'DECODE_XLAI bad value of k=',k
      ENDIF
      xlai(i,j,ic) = PTR3D(i,j,k)
     END DO

     DO ic = ireg(i,j)+1,NTYPE
      xlai(i,j,ic) = 0.0d0
     END DO

    END DO
   END DO

  RETURN
 END SUBROUTINE Decode_XLAI


 END MODULE VegLaiMod
