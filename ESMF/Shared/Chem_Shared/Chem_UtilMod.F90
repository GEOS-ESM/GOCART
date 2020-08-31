#include "MAPL_Generic.h"
#include "unused_dummy.H"
#if 0
#define min(x,y) amin(real(x),real(y))
#define MIN(x,y) AMIN(real(x),real(y))
#define max(x,y) amax(real(x),real(y))
#define MAX(x,y) AMAX(real(x),real(y))
#endif

#define DEALOC_(A) if(associated(A)) then; A=0; call MAPL_DeAllocNodeArray(A,rc=STATUS); if(STATUS==MAPL_NoShm) deallocate(A, stat=STATUS); VERIFY_(STATUS); NULLIFY(A); endif

#define DEALOC2_(A) if(associated(A)) then; deallocate(A, stat=STATUS); VERIFY_(STATUS); NULLIFY(A); endif

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  Chem_UtilMod --- Assorted Utilities for fvChem
!
! !INTERFACE:
!

   module  Chem_UtilMod

! !USES:

   use ESMF
   use MAPL

   use Chem_Mod                  ! Chemistry Base Class
   use mod_diag                  ! fvGCM diagnostics
   use m_die
   use m_StrTemplate             ! string templates
   use m_chars, only: uppercase

   implicit NONE

!
! !PUBLIIC MEMBER FUNCTIONS:
!

   PRIVATE
   PUBLIC  Chem_UtilNegFiller       ! Fills negative values in a column
   PUBLIC  Chem_UtilTroppFixer      ! Repairs tropopause pressure bad values
   PUBLIC  Chem_UtilGetTimeInfo     ! Time info on file
   PUBLIC  Chem_UtilExtractIntegers ! Extract integers from a character-delimited string
   PUBLIC  Chem_BiomassDiurnal      ! Biomass burning diurnal cycle
   PUBLIC  Chem_UtilResVal          ! Finds resolution dependent value that corresponds 
                                    ! to the model resolution
   PUBLIC  Chem_UtilIdow            ! Integer day of the week: Sun=1, Mon=2, etc.
   PUBLIC  Chem_UtilCdow            ! String day of the week: 'Sun', 'Mon', etc.

   PUBLIC  Chem_UtilPointEmissions  ! From a provided list returns a table of
                                    ! pointwise emissions (e.g., for volcanoes
                                    ! or wildfires)

   PUBLIC tick      ! GEOS-4 stub
   PUBLIC mcalday   ! GEOS-4 stub
   PUBLIC pmaxmin   ! functional
   PUBLIC zenith    ! GEOS-4 stub

!
! !DESCRIPTION:
!
!  This module implements assorted odds & ends for fvChem.
!
! !REVISION HISTORY:
!
!  29oct2003  da Silva  First crack.
!  16aug2005  da Silva  Introduced scatter from MAPL_CommsMod.
!
!EOP
!-------------------------------------------------------------------------

   interface pmaxmin
     module procedure pmaxmin2d
     module procedure pmaxmin3d
   end interface

CONTAINS

#ifdef USE_MAPL_MPREAD

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Chem_UtilMPread --- Reads fields from file and distribute
!
! !INTERFACE:
!
   subroutine Chem_UtilMPread_g5 ( filen, varn, nymd, nhms, &
                                i1, i2, ig, im, j1, j2, jg, jm, km, &
                                var2d, var3d, cyclic, grid )

! !USES:

  implicit NONE

! !INPUT PARAMETERS:

  character(len=*), intent(in) :: filen   ! GFIO compatible file name
  character(len=*), intent(in) :: varn    ! variable name
  integer, intent(in)          :: nymd, nhms ! date/time

                                          ! Distributed grid info:
  integer,      intent(in)     :: i1, i2  !   local  zonal indices
  integer,      intent(in)     :: ig      !   zonal ghosting
  integer,      intent(in)     :: im      !   global zonal dimension
  integer,      intent(in)     :: j1, j2  !   local meridional indices
  integer,      intent(in)     :: jg      !   meridional ghosting
  integer,      intent(in)     :: jm      !   global zonal dimension
  integer,      intent(in)     :: km      !   vertical dimension

  logical, OPTIONAL, intent(in) :: cyclic ! whether time dimension is periodic
 
                                          ! ESMF Grid; this is required
                                          !  in GEOS-5 under ESMF
  type(ESMF_Grid), OPTIONAL, intent(in) :: grid 

! !OUTPUT PARAMETERS:

  real, OPTIONAL, intent(out), target :: var2d(i1-ig:i2+ig,j1-jg:j2+jg)
  real, OPTIONAL, intent(out), target :: var3d(i1-ig:i2+ig,j1-jg:j2+jg,km)

! !DESCRIPTION: 
!
! !REVISION HISTORY:
!
!  15Aug2006  da Silva  It is now a simple wrap around CFIOReadArray().
!
!EOP
!-------------------------------------------------------------------------

    character(len=*), parameter ::  Iam = 'Chem_UtilMPread_g5'
    type(ESMF_TIME)             :: time
    integer                     :: yy, mm, dd, h, m, s, rc, STATUS

    real, pointer               :: ptr2(:,:), ptr3(:,:,:)

    
!   Convert time to ESMF format
!   ---------------------------
    call parseIntTime_ ( nymd, yy, mm, dd )
    call parseIntTime_ ( nhms,  h,  m,  s )
    call ESMF_TimeSet(time, yy=yy, mm=mm, dd=dd,  h=h,  m=m, s=s, rc=status )
    if ( status /= 0 ) call die(Iam,'failed to convert time')

!   Read either 2D or 3D array
!   --------------------------
    if ( .not. present(grid) ) &
       call die ( Iam,'when running under the ESMF "grid" must be specified' )

    if ( present(var2d) ) then

         ptr2 => var2d
         call MAPL_CFIORead ( varn, filen, time, grid, ptr2, rc=STATUS,  &
                              verbose = .true., time_is_cyclic=cyclic,   &
                              time_interp = .true. )

    else if ( present(var3d) ) then

         ptr3 => var3d
         call MAPL_CFIORead ( varn, filen, time, grid, ptr3, rc=STATUS,  &
                              verbose = .true., time_is_cyclic=cyclic,   &
                              time_interp = .true. )

    else

       call die ( Iam,'either "var2d" or "var3d" must be specified' )

    end if

    if ( status /= 0 ) call die(Iam,'cannot read '//trim(varn))


CONTAINS
    subroutine parseIntTime_ ( hhmmss, hour, min, sec )      
      integer, intent(in)  :: hhmmss
      integer, intent(out) :: hour, min, sec 
      hour = hhmmss / 10000
      min  = mod(hhmmss,10000)/100
      sec  = mod(hhmmss,100)
    end subroutine parseIntTime_

  end subroutine Chem_UtilMPread_g5

#endif


subroutine Chem_UtilGetTimeInfo ( fname, begDate, begTime, nTimes, incSecs )

  implicit NONE
  character(len=*), intent(in) :: fname    ! GFIO/CFIO filename
  integer, intent(out) :: begDate, begTime ! initial time/date on file
                                           ! given as YYYYMMDD and HHMMSS
  integer, intent(out) :: nTimes           ! number of time steps on file
  integer, intent(out) :: incSecs          ! time steps in seconds
 

  integer :: READ_ONLY=1
  integer :: fid, rc, im,jm,km,nvars,ngatts
  character(len=*), parameter ::  myname = 'Chem_UtilGetTimeInfo'
  
! Special case
! -----------
  if ( fname(1:9) == '/dev/null' ) then    
     begDate=0; begTime=0;  nTimes=0; incSecs=0
     return
  end if

! Open file
! ---------
  call GFIO_Open ( fname, READ_ONLY, fid, rc )
  if ( rc /= 0 ) call die(myname, 'Unable to open '// trim(fname))
  if ( rc .ne. 0 ) then
     begDate=-1; begTime=-1;  nTimes=-1; incSecs=-1
     return
  endif

! Get dimension sizes
! -------------------
  call GFIO_DimInquire (fid,im,jm,km,nTimes,nvars,ngatts,rc)
  if ( rc /= 0 ) call die(myname, 'Unable to inquire about dimension sizes for file '// trim(fname))
  if ( rc .ne. 0 ) then
     begDate=-1; begTime=-1;  nTimes=-1; incSecs=-1
     return
  endif

! Get initial time/timestep
! -------------------------
  call GetBegDateTime ( fid, begDate, begTime, incSecs, rc )
  if ( rc .ne. 0 ) then
     begDate=-1; begTime=-1;  nTimes=-1; incSecs=-1
     return
  endif

  call GFIO_close(fid,rc)

end subroutine Chem_UtilGetTimeInfo

!............................... geos4 stubs  ..........................

! Parallelized utility routine for computing/printing
! max/min of an input array
!
      subroutine pmaxmin3d ( qname, a, pmin, pmax, im, jt, fac )
      implicit none
      character*(*)  qname
      integer im, jt
      real :: a(:,:,:)
      real pmax, pmin
      real fac                     ! multiplication factor
      call pmaxmin2d ( qname, reshape(a,(/ im, jt /)), &
                             pmin, pmax, im, jt, fac )
      
      end subroutine pmaxmin3d

      subroutine pmaxmin2d ( qname, a, pmin, pmax, im, jt, fac )

      implicit none

      character*(*)  qname
      integer im, jt
      real a(im,jt)
      real pmax, pmin
      real fac                     ! multiplication factor

      integer :: i, j, two=2

      real qmin(jt), qmax(jt)
      real pm1(2)
      real pm_res(2)
      type(ESMF_VM) :: vm

      character(len=32) :: name
      integer :: status


      call ESMF_VmGetCurrent(vm=vm, rc=status)
!$omp parallel do private(i, j, pmax, pmin)

      do j=1,jt
         pmax = a(1,j)
         pmin = a(1,j)
         do i=2,im
            pmax = max(pmax, a(i,j))
            pmin = min(pmin, a(i,j))
         enddo
         qmax(j) = pmax
         qmin(j) = pmin
      enddo
!
! Now find max/min of amax/amin
!
      pmax = qmax(1)
      pmin = qmin(1)
      do j=2,jt
         pmax = max(pmax, qmax(j))
         pmin = min(pmin, qmin(j))
      enddo

      pm1(1) = pmax
      pm1(2) = -pmin
      call MAPL_CommsAllReduceMax(vm, sendbuf=pm1, recvbuf=pm_res, cnt=two, RC=status)
      pmax=pm_res(1)
      pmin=-pm_res(2)
     
      if ( fac /= 0.0 ) then  ! trick to prevent printing
         if ( MAPL_am_I_root() ) then
            name = '            '
            name(1:len(qname)) = qname
            write(*,*) name, ' max = ', pmax*fac, ' min = ', pmin*fac
            return
         end if
      end if

      return

    end subroutine pmaxmin2d


      function leap_year(ny)
!
! Determine if year ny is a leap year
!
! Author: S.-J. Lin
      implicit none
      logical leap_year
      integer ny
      integer ny00

!
! No leap years prior to 0000
!
      parameter ( ny00 = 0000 )   ! The threshold for starting leap-year 

      if( ny >= ny00 ) then
         if( mod(ny,100) == 0. .and. mod(ny,400) == 0. ) then
             leap_year = .true.
         elseif( mod(ny,4) == 0. .and. mod(ny,100) /= 0.  ) then
             leap_year = .true.
         else
             leap_year = .false.
         endif
      else
          leap_year = .false.
      endif

      return 
    end function leap_year


      integer FUNCTION INCYMD (NYMD,M)

!  PURPOSE
!     INCYMD:  NYMD CHANGED BY ONE DAY
!     MODYMD:  NYMD CONVERTED TO JULIAN DATE
!  DESCRIPTION OF PARAMETERS
!     NYMD     CURRENT DATE IN YYMMDD FORMAT
!     M        +/- 1 (DAY ADJUSTMENT)

      integer nymd, m, ny, nm, nd

      INTEGER NDPM(12)
      DATA    NDPM /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
!!!      logical leap_year

      NY = NYMD / 10000
      NM = MOD(NYMD,10000) / 100
      ND = MOD(NYMD,100) + M

      IF (ND.EQ.0) THEN
      NM = NM - 1
      IF (NM.EQ.0) THEN
          NM = 12
          NY = NY - 1
      ENDIF
      ND = NDPM(NM)
      IF (NM.EQ.2 .AND. leap_year(NY))  ND = 29
      ENDIF

      IF (ND.EQ.29 .AND. NM.EQ.2 .AND. leap_year(ny))  GO TO 20

      IF (ND.GT.NDPM(NM)) THEN
      ND = 1
      NM = NM + 1
      IF (NM.GT.12) THEN
          NM = 1
          NY = NY + 1
      ENDIF
      ENDIF

   20 CONTINUE
      INCYMD = NY*10000 + NM*100 + ND
      RETURN
    END FUNCTION INCYMD

      subroutine tick (nymd, nhms, ndt)

! Input:
      integer ndt                     ! TIME-STEP
! Inpuit/Output:
      integer nymd                    ! CURRENT YYYYMMDD
      integer nhms                    ! CURRENT HHMMSS
!!!      integer incymd

! Revision:   S.-J. Lin Mar 2000
       integer nsecf, nhmsf, n, nsec

       NSECF(N)   = N/10000*3600 + MOD(N,10000)/100* 60 + MOD(N,100)
       NHMSF(N)   = N/3600*10000 + MOD(N,3600 )/ 60*100 + MOD(N, 60)

       NSEC = NSECF(NHMS) + ndt

       IF (NSEC.GT.86400)  THEN
           DO WHILE (NSEC.GT.86400)
              NSEC = NSEC - 86400
              NYMD = INCYMD (NYMD,1)
           ENDDO
       ENDIF

       IF (NSEC.EQ.86400)  THEN
           NSEC = 0
           NYMD = INCYMD (NYMD,1)
       ENDIF

       IF (NSEC .LT. 0)  THEN
           DO WHILE (NSEC .LT. 0)
               NSEC = 86400 + NSEC
               NYMD = INCYMD (NYMD,-1)
           ENDDO
        ENDIF

          NHMS = NHMSF (NSEC)
      return
    end subroutine tick

      subroutine mcalday(nymd, nhms, calday)
      implicit none

! input:
      integer nymd
      integer nhms
! Output:
      real calday                    ! Julian day (1 to 366 for non-leap year)
                                     ! Julian day (-1 to -367 for   leap year)
! Local:
      real tsec
      integer n, nsecf, m, mm
      integer dd, ds
      integer days(12)
      integer ny

      data days /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
      nsecf(n)  = n/10000*3600 + mod(n,10000)/100* 60 + mod(n,100)

      ny = nymd / 10000
      mm = mod(nymd, 10000) / 100 
      dd = mod(nymd,   100)

      ds = dd -1

      if( mm .ne. 1) then
      do m=1, mm-1
         if( m.eq.2  .and. leap_year(ny) ) then 
             ds = ds + 29
         else
             ds = ds + days(m)
         endif
      enddo
      endif

      tsec = ds * 86400 + nsecf(nhms)

      calday = tsec / 86400.  + 1.
      if( leap_year(ny) ) calday = -calday

      return
    end subroutine mcalday

      subroutine zenith(calday  ,dodiavg ,clat    ,coszrs  )

!
! Input arguments
!
      real calday              ! Calendar day, including fraction
      logical dodiavg          ! true => do diurnal averaging
      real clat                ! Current latitude (radians)
!
! Output arguments
!
      real coszrs(*)       ! Cosine solar zenith angle
!
!---------------------------Local variables-----------------------------
!
      _UNUSED_DUMMY(calday)
      _UNUSED_DUMMY(dodiavg)
      _UNUSED_DUMMY(clat)
      _UNUSED_DUMMY(coszrs(1)) ! assumed-size

      call die ('zenith','stub only, please do not call zenith_()' )

    end subroutine zenith

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Chem_UtilNegFiller --- Negative Filler
!
! !INTERFACE:
!
   subroutine Chem_UtilNegFiller ( q, delp, in, jn, qmin )

! !USES:

  implicit NONE

! !INPUT/OUTPUT PARAMETERS:

  integer                    :: in, jn             ! number of local lon/lat 
  real, pointer              :: delp(:,:,:)
  real, OPTIONAL, intent(in) :: qmin

! !OUTPUT PARAMETERS:

  real, pointer :: q(:,:,:)     ! 3D tracer

! !DESCRIPTION: 
!
! !REVISION HISTORY: Makes sure tracer has no negative values. This is
!                    a "flat tax" algorithm: first negative values are
!  replaced with tiny() or user specified value. Then profiles are rescaled 
!  to preserve column mass, whenever possible. No mass conservation is
!  imposed when the initial column mass is negative or zero.
!
!  18May2007  da Silva  First crack.
!
!EOP
!-------------------------------------------------------------------------

  real :: mass_1(in,jn), mass_2(in,jn), ratio(in,jn)
  real :: qmin_
#ifdef DEBUG
  real :: vmax, vmin
#endif
  integer :: k, k1, k2, km

  k1 = lbound(q,3)
  k2 = ubound(q,3)
  km = k2 - k1 + 1

! Unless specified, minimum is the smallest positive float, not zero
! ------------------------------------------------------------------
  if ( present(qmin) ) then
     qmin_ = qmin
  else  
     qmin_ = tiny(1.0)
  end if

#ifdef DEBUG
  call pmaxmin ( 'NegFill:  q_beg', q, vmax, vmin, in*jn, km, 1. )
#endif

! Column mass before fixer
! ------------------------
  mass_1 = sum ( delp * q, 3 )

! Cap q
! -----
  where ( q < qmin_ ) q = qmin_

! Enforce conservation of column mass
! -----------------------------------
  mass_2 = sum ( delp * q, 3 )
  where ( (mass_2 /= mass_1) .AND. (mass_1 > 0.0) )
          ratio = mass_1 / mass_2
  elsewhere
          ratio = 1.0
  end where

! Next correct q in each layer
! ----------------------------
  do k = k1, k2
     where ( ratio /= 1.0 )
             q(:,:,k) = ratio * q(:,:,k)
     end where
  end do

#ifdef DEBUG
  call pmaxmin ( 'NegFill: mass_1', mass_1, vmax, vmin, in*jn,1, 1. )
  call pmaxmin ( 'NegFill: mass_2', mass_2, vmax, vmin, in*jn,1, 1. )
  call pmaxmin ( 'NegFill:  ratio',  ratio, vmax, vmin, in*jn,1, 1. )
  call pmaxmin ( 'NegFill:  q_end', q, vmax, vmin, in*jn, km, 1. )
#endif

end subroutine Chem_UtilNegFiller

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Chem_UtilTroppFixer - Repair tropopause pressure bad values
!
! !INTERFACE:
!
  SUBROUTINE Chem_UtilTroppFixer(im, jm, tropp, threshold, verbose, newTropp, rc)

! !USES:

! !USES:

  IMPLICIT NONE

! !INPUT/OUTPUT PARAMETERS:

  INTEGER, INTENT(IN)		:: im          !Index range, longitude
  INTEGER, INTENT(IN)		:: jm          !Index range, latitude
  REAL, POINTER , INTENT(INOUT)	:: tropp(:,:)  !Tropopause pressure (Pa)
  REAL, OPTIONAL, INTENT(IN)    :: threshold   !User-supplied threshold (Pa).  If
                                               ! not present, default is 110000 Pa.
  LOGICAL, OPTIONAL, INTENT(IN) :: verbose     !Write message when bad value is 
                                               ! encountered. DEBUG directive turns
                                               ! on the message even if verbose is not
                                               ! present or if verbose = .FALSE.
  INTEGER, INTENT(OUT)		:: rc          !Return code: 0=OK, 1=Unable to repair


! !OUTPUT PARAMETERS:

  REAL, OPTIONAL, INTENT(OUT)   :: newTropp(1:im,1:jm)
                                               !If present, fill with repaired
                                               ! tropopause pressures (Pa).  Otherwise
                                               ! tropp is overwritten.

! !DESCRIPTION: 
!
! Replace bad values in tropopause pressure array with an average of nearby good
! values.  Working outward from the afflicted cell, nearby cells are scanned for 
! the presence of good values. In the first scan, all adjacent cells, which can
! number up to eight, are considered.  If at least one of the cells has a valid
! pressure, the scanning is terminated.  Otherwise, the "radius" of the scan is
! increased by one, and up to 24 cells are considered, and so on, until, at the
! extreme, all cells on the current processor fall under consideration.
!
! After the scanning is done, the bad value is replaced with the average the 
! valid pressures found by the scan.  Thus, the accuracy of the replaced value
! drops rapidly as the scan expands outward.  At the same time, the only case
! in which a valid pressure will not be found to replace a bad one is if all 
! pressures on the current processor are invalid.  In this case the return 
! code is set to 1.
!
! If newTropp is not present, then the input array, tropp, is overwritten.  If
! it is present, it is filled, even if there are no invalid pressures in tropp.
!
! No scanning is done if all tropopause pressures are valid.
!
! Assumptions/bugs:
!
! Bad values are high, but not Indef, and the primary purpose of this routine
! is to repair the case where GEOS-5 fails to find the tropopause and assigns
! MAPL_UNDEF as the tropopause pressure.  We recommend using the blended 
! tropopause values for tropp, because the frequency of bad values is quite
! rare compared to the unblended case.
!
! !REVISION HISTORY: 
!
!  25Jan2008  Nielsen  Initial code and testing.
!
!EOP
!-------------------------------------------------------------------------
  CHARACTER(LEN=*), PARAMETER :: myName = 'Chem_UtilTroppFixer'

  LOGICAL :: tellMe

  INTEGER :: i,ier,j,m
  INTEGER :: ie,iw,jn,js

  INTEGER, ALLOCATABLE :: mask(:,:)
  INTEGER, ALLOCATABLE :: mx(:)
  REAL, ALLOCATABLE :: p(:,:)
  
  REAL :: badValue, r

  rc = 0

! Determine verbosity, letting the DEBUG 
! directive override local specification
! --------------------------------------
  tellMe = .FALSE.
  IF(PRESENT(verbose)) THEN
   IF(verbose) tellMe = .TRUE.
  END IF
#ifdef DEBUG
  tellMe = .TRUE.
#endif

! Set the bad value to 110000 Pa (1100 hPa)
! -----------------------------------------
  IF(PRESENT(threshold)) THEN
   badValue = threshold
  ELSE
   badValue = 1.10E+05 !Pa, 1100 hPa
  END IF

! There may be no bad values ...
! ------------------------------
  IF( ALL( tropp(1:im,1:jm) < badValue ) ) THEN
   IF(PRESENT(newTropp)) newTropp(1:im,1:jm) = tropp(1:im,1:jm)
   RETURN
  END IF

! ... or there is at least one bad value
! --------------------------------------
  ALLOCATE(mask(1:im,1:jm),STAT=ier)
  ALLOCATE(p(1:im,1:jm),STAT=ier)
  ALLOCATE(mx(4),STAT=ier)

! Loop over each cell
! -------------------
  DO j=1,jm
   DO i=1,im

! Invalid pressure found at cell(i,j)
! -----------------------------------
    IF(tropp(i,j) >= badValue) THEN

! Determine maximum "radius" of search
! ------------------------------------
     mx(1) = im-i
     mx(2) = i-1
     mx(3) = jm-j
     mx(4) = j-1

! Start search
! ------------
     DO m=1,MAXVAL(mx)

! Clear the mask
! --------------
      mask(1:im,1:jm) = 0

! Range of search
! ---------------
      iw = MAX( 1,i-m)
      ie = MIN(im,i+m)
      js = MAX( 1,j-m)
      jn = MIN(jm,j+m)

! Set mask to one for cells in range of search
! --------------------------------------------
       mask(iw:ie,js:jn) = 1

! Set mask back to zero for cells in range
! of search that have invalid pressures.
! ----------------------------------------
       WHERE(tropp(iw:ie,js:jn) >= badValue) mask(iw:ie,js:jn) = 0

! One valid pressure is enough ...
! --------------------------------
       IF(SUM(MASK) >= 1) EXIT

! ... or "radius" of search needs to be extended
! ----------------------------------------------
     END DO

! Repair bad value at cell(i,j) with average 
! of valid pressures found in range of search
! -------------------------------------------
     r = SUM(tropp,mask == 1)
     p(i,j) = r/(1.00*SUM(mask))

! For debugging
! -------------
!    IF(tellMe) THEN
!     WRITE(*,FMT="(A,': ',ES12.5,' becomes ',ES12.5,' Pa [',I4,2X,I4,']')") &
!           TRIM(myName),tropp(i,j),p(i,j),m,SUM(mask)
!    END IF

    ELSE

! Input pressure at cell(i,j) was valid
! -------------------------------------
     p(i,j) = tropp(i,j)

    END IF

! Next cell
! ---------
   END DO
  END DO

! Clean up
! --------
  DEALLOCATE(mask,STAT=ier)
  DEALLOCATE(mx,STAT=ier)

! If all cells have bad values, then 
! register a failure, but continue.
! ----------------------------------
  IF( ANY( p(1:im,1:jm) >= badValue ) ) THEN
   PRINT *, myName,": WARNING Unable to fix bad tropopause pressure(s)"
   rc = 1
  END IF
  
! Overwrite input or fill output array
! ------------------------------------
  IF(PRESENT(newTropp)) THEN
   newTropp(1:im,1:jm) = p(1:im,1:jm)
  ELSE
   tropp(1:im,1:jm) = p(1:im,1:jm)
  END IF

! Clean up some more
! ------------------
  DEALLOCATE(p,STAT=ier)

  RETURN
  END SUBROUTINE Chem_UtilTroppFixer

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  Chem_UtilExtractIntegers - Extract integers from a delimited string
!
! !INTERFACE:
!
  SUBROUTINE Chem_UtilExtractIntegers(string,iSize,iValues,delimiter,verbose,fillValue,rc)

! !USES:

  IMPLICIT NONE

! !INPUT/OUTPUT PARAMETERS:

  CHARACTER(LEN=*), INTENT(IN)   :: string	  ! Character-delimited string of integers
  INTEGER, INTENT(IN)            :: iSize
  INTEGER, INTENT(INOUT)         :: iValues(iSize)! Space allocated for extracted integers
  CHARACTER(LEN=*), OPTIONAL     :: delimiter     ! 1-character delimiter
  LOGICAL, OPTIONAL, INTENT(IN)  :: verbose	  ! Let me know iValues as they are found. 
                                		  ! DEBUG directive turns on the message even 
                                		  ! if verbose is not present or if 
                                		  ! verbose = .FALSE.
  INTEGER, OPTIONAL, INTENT(IN)  :: fillValue     ! unfilled iValue entries get this value (default=0)
  INTEGER, OPTIONAL, INTENT(OUT) :: rc            ! Return code

! !DESCRIPTION: 
!
!  Extract integers from a character-delimited string, for example, "-1,45,256,7,10".  In the context
!  of Chem_Util, this is provided for determining the numerically indexed regions over which an 
!  emission might be applied.
!
!  In multiple passes, the string is parsed for the delimiter, and the characters up to, but not
!  including the delimiter are taken as consecutive digits of an integer.  A negative sign ("-") is
!  allowed.  After the first pass, each integer and its trailing delimiter are lopped of the head of
!  the (local copy of the) string, and the process is started over.
!
!  The default delimiter is a comma (",").
!
!  "Unfilled" iValues get set to fillValue.
!  
!  Return codes:
!  1 Zero-length string.
!  2 iSize needs to be increased.
!
!  Assumptions/bugs:
!
!  A non-zero return code does not stop execution.
!  Allowed numerals are: 0,1,2,3,4,5,6,7,8,9.
!  A delimiter must be separated from another delimiter by at least one numeral.
!  The delimiter cannot be a numeral or a negative sign.
!  The character following a negative sign must be an allowed numeral.
!  The first character must be an allowed numeral or a negative sign.
!  The last character must be an allowed numeral.
!  The blank character (" ") cannot serve as a delimiter.
!
!  Examples of strings that will work:
!  "1"
!  "-1"
!  "-1,2004,-3"
!  "1+-2+3"
!  "-1A100A5"
!
!  Examples of strings that will not work:
!  "1,--2,3"
!  "1,,2,3"
!  "1,A,3"
!  "1,-,2"
!  "1,2,3,4,"
!  "+1"
!  "1 3 6"
!
! !REVISION HISTORY: 
!
!  29Feb2008  Nielsen  Initial code and testing.
!  25Jul2014  Manyin   Added fillValue arg
!
!EOP
!-------------------------------------------------------------------------
 CHARACTER(LEN=*), PARAMETER :: myName = 'Chem_UtilExtractIntegers'

 INTEGER :: base,count,i,iDash,last,lenStr
 INTEGER :: multiplier,pos,posDelim,sign
 CHARACTER(LEN=255) :: str
 CHARACTER(LEN=1) :: char,delimChar
 LOGICAL :: Done
 LOGICAL :: tellMe

! Initializations
! ---------------
 rc = 0
 count = 1
 Done = .FALSE.
 iValues(:) = 0
 base = ICHAR("0")
 iDash = ICHAR("-")

! Determine verbosity, letting the DEBUG 
! directive override local specification
! --------------------------------------
  tellMe = .FALSE.
  IF(PRESENT(verbose)) THEN
   IF(verbose) tellMe = .TRUE.
 END IF
#ifdef DEBUG
  tellMe = .TRUE.
#endif

! Check for zero-length string
! ----------------------------
 lenStr = LEN_TRIM(string)
 IF(lenStr == 0) THEN
  rc = 1
  PRINT *,myname,": ERROR - Found zero-length string."
  RETURN
 END IF

! Default delimiter is a comma
! ----------------------------
 delimChar = ","
 IF(PRESENT(delimiter)) delimChar(1:1) = delimiter(1:1)

! Work on a local copy
! --------------------
 str = TRIM(string)

! One pass for each delimited integer
! -----------------------------------
 Parse: DO

  lenStr = LEN_TRIM(str)

! Parse the string for the delimiter
! ----------------------------------
  posDelim = INDEX(TRIM(str),TRIM(delimChar))
  IF(tellMe) PRINT *,myname,": Input string is >",TRIM(string),"<"

! If the delimiter does not exist,
! one integer remains to be extracted.
! ------------------------------------
  IF(posDelim == 0) THEN
   Done = .TRUE.
   last = lenStr
  ELSE
   last = posDelim-1
  END IF
  multiplier = 10**last

! Examine the characters of this integer
! --------------------------------------
  Extract: DO pos=1,last

   char = str(pos:pos)
   i = ICHAR(char)

! Account for a leading "-"
! -------------------------
   IF(pos == 1) THEN 
    IF(i == iDash) THEN
     sign = -1
    ELSE
     sign = 1
    END IF
   END IF

! "Power" of 10 for this character
! --------------------------------
   multiplier = multiplier/10

   IF(pos == 1 .AND. sign == -1) CYCLE Extract

! Integer comes from remaining characters
! ---------------------------------------
   i = (i-base)*multiplier
   iValues(count) = iValues(count)+i
   IF(pos == last) THEN
    iValues(count) = iValues(count)*sign
    IF(tellMe) PRINT *,myname,":Integer number ",count," is ",iValues(count)
   END IF
  
  END DO Extract

  IF(Done) THEN
    IF(PRESENT(fillValue) .AND. (count < iSize)) iValues((count+1):iSize) = fillValue
    EXIT
  END IF

! Lop off the leading integer and try again
! -----------------------------------------
  str(1:lenStr-posDelim) = str(posDelim+1:lenStr)
  str(lenStr-posDelim+1:255) = " "
  count = count+1

! Check size
! ----------
  IF(count > iSize) THEN
   rc = 2
   PRINT *,myname,": ERROR - iValues does not have enough elements."
  END IF

 END DO Parse

 RETURN
END SUBROUTINE Chem_UtilExtractIntegers

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This is a candidate for ESMFL, here for dependency reasons
!  

  subroutine GridGetLatLons_ ( grid, lons, lats )

    implicit NONE
    type(ESMF_Grid) :: grid

    real, pointer   :: lons(:), lats(:)

!                     ---

    character(len=*), parameter :: Iam = 'GridGetLatLons'

    real(KIND=8), pointer  :: R8D2(:,:)
    real, pointer          :: lons2d(:,:), lats2d(:,:)
    real, pointer          :: LONSLocal(:,:), LATSlocal(:,:)
    integer                :: IM_WORLD, JM_WORLD, dims(3), STATUS, RC

!                          ----

!      Get world dimensions
!      --------------------
       call MAPL_GridGet(GRID, globalCellCountPerDim=DIMS, RC=STATUS)

       IM_WORLD = dims(1)
       JM_WORLD = dims(2)

!      Allocate memory for output if necessary
!      ---------------------------------------
       if ( .not. associated(lons) ) then
            allocate(lons(IM_WORLD), stat=STATUS)
            VERIFY_(status)
       else
            if(size(LONS,1) /= IM_WORLD) STATUS = 1
            VERIFY_(status)
       end if
       if ( .not. associated(lats) ) then
            allocate(lats(JM_WORLD), stat=STATUS)
            VERIFY_(status)
       else
            if(size(LATS,1) /= JM_WORLD) STATUS = 1
            VERIFY_(status)
       end if

!      Local work space
!      ----------------
       allocate(LONS2d(IM_WORLD,JM_WORLD), LATS2d(IM_WORLD,JM_WORLD), &
                STAT=status)             
       VERIFY_(status)
       LONS2d=0
       LATS2d=0

!      Get the local longitudes and gather them into a global array
!      ------------------------------------------------------------
       call ESMF_GridGetCoord(grid, localDE=0, coordDim=1, &
             staggerloc=ESMF_STAGGERLOC_CENTER, &
             datacopyFlag = ESMF_DATACOPY_REFERENCE,       &
             farrayPtr=R8D2, rc=status)

       allocate(LONSLOCAL(size(R8D2,1),size(R8D2,2)), STAT=status)             
       VERIFY_(status)

       LONSLOCAL = R8D2*(180/MAPL_PI)

       call ArrayGather(LONSLOCAL, LONS2D, GRID, RC=STATUS)

!      Get the local longitudes and gather them into a global array
!      ------------------------------------------------------------
       call ESMF_GridGetCoord(grid, localDE=0, coordDim=2, &
             staggerloc=ESMF_STAGGERLOC_CENTER, &
             datacopyFlag = ESMF_DATACOPY_REFERENCE,       &
             farrayPtr=R8D2, rc=status)

       allocate(LATSLOCAL(size(R8D2,1),size(R8D2,2)), STAT=status)             
       VERIFY_(status)

       LATSlocal = R8D2*(180/MAPL_PI)

       call ArrayGather(LATSLOCAL, LATS2D, GRID, RC=STATUS)
       VERIFY_(STATUS)

!      Return 1D arrays
!      ----------------
       LONS = LONS2D(:,1)
       LATS = LATS2D(1,:)

       DEALLOCATE(LONSLOCAL, LATSLOCAL, LONS2d, LATS2d )
        
     end subroutine GridGetLatLons_

!..........................................................................................


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  Chem_BiomassDiurnal - Applies diurnal cycle to biomass emissions.
!
! !INTERFACE:
!

     subroutine Chem_BiomassDiurnal ( Eout, Ein, lons, lats, nhms, cdt)

! !USES:

  IMPLICIT NONE

! !ARGUMENTS:

       real, intent(out)   :: Eout(:,:) ! Emissions valid at NHMS
       real, intent(in)    :: Ein(:,:)  ! Daily-mean emissions
       real, intent(in)    :: lons(:,:) ! Latitudes in degrees
       real, intent(in)    :: lats(:,:) ! Latitudes in degrees
       integer, intent(in) :: nhms
       real, intent(in)    :: cdt       ! time step in seconds

! !DESCRIPTION: 
!
!      Applies diurnal cycle to biomass emissions.       
!
! !DESCRIPTION:
!
!  This module implements assorted odds & ends for fvChem.
!
! !REVISION HISTORY:
!
!  13nov2009  da Silva  First crack.
!
!EOP
!-------------------------------------------------------------------------

!      Hardwired diurnal cycle (multiplied by 100)
!      These numbers were derived from GOES-12
!      fire counts for 2003-2007.
!      -------------------------------------------
       integer, parameter :: N = 240
       real,    parameter :: DT = 86400. / N

!      Apply flat diurnal cycle for boreal forests as a 
!      temporary solution to prevent very high aerosol
!      optical depth during the day
       real,    parameter :: Boreal(N) = 1.0
!      real,    parameter :: Boreal(N) = &
!      (/ 0.0277, 0.0292, 0.0306, 0.0318, 0.0327, 0.0335, &
!         0.0340, 0.0342, 0.0341, 0.0338, 0.0333, 0.0326, &
!         0.0316, 0.0305, 0.0292, 0.0278, 0.0263, 0.0248, &
!         0.0233, 0.0217, 0.0202, 0.0187, 0.0172, 0.0158, &
!         0.0145, 0.0133, 0.0121, 0.0110, 0.0100, 0.0091, &
!         0.0083, 0.0075, 0.0068, 0.0062, 0.0056, 0.0051, &
!         0.0046, 0.0042, 0.0038, 0.0035, 0.0032, 0.0030, &
!         0.0028, 0.0026, 0.0025, 0.0024, 0.0024, 0.0024, &
!         0.0024, 0.0026, 0.0027, 0.0030, 0.0033, 0.0036, &
!         0.0041, 0.0046, 0.0052, 0.0060, 0.0069, 0.0079, &
!         0.0090, 0.0104, 0.0119, 0.0137, 0.0157, 0.0180, &
!         0.0205, 0.0235, 0.0268, 0.0305, 0.0346, 0.0393, &
!         0.0444, 0.0502, 0.0565, 0.0634, 0.0711, 0.0794, &
!         0.0884, 0.0982, 0.1087, 0.1201, 0.1323, 0.1453, &
!         0.1593, 0.1742, 0.1900, 0.2069, 0.2249, 0.2439, &
!         0.2642, 0.2858, 0.3086, 0.3329, 0.3587, 0.3860, &
!         0.4149, 0.4455, 0.4776, 0.5115, 0.5470, 0.5840, &
!         0.6227, 0.6628, 0.7043, 0.7470, 0.7908, 0.8355, &
!         0.8810, 0.9271, 0.9735, 1.0200, 1.0665, 1.1126, &
!         1.1580, 1.2026, 1.2460, 1.2880, 1.3282, 1.3664, &
!         1.4023, 1.4356, 1.4660, 1.4933, 1.5174, 1.5379, &
!         1.5548, 1.5679, 1.5772, 1.5826, 1.5841, 1.5818, &
!         1.5758, 1.5661, 1.5529, 1.5365, 1.5169, 1.4944, &
!         1.4693, 1.4417, 1.4119, 1.3801, 1.3467, 1.3117, &
!         1.2755, 1.2383, 1.2003, 1.1616, 1.1225, 1.0832, &
!         1.0437, 1.0044, 0.9653, 0.9265, 0.8882, 0.8504, &
!         0.8134, 0.7771, 0.7416, 0.7070, 0.6734, 0.6407, &
!         0.6092, 0.5787, 0.5493, 0.5210, 0.4939, 0.4680, &
!         0.4433, 0.4197, 0.3974, 0.3763, 0.3565, 0.3380, &
!         0.3209, 0.3051, 0.2907, 0.2777, 0.2662, 0.2561, &
!         0.2476, 0.2407, 0.2352, 0.2313, 0.2289, 0.2279, &
!         0.2283, 0.2300, 0.2329, 0.2369, 0.2417, 0.2474, &
!         0.2536, 0.2602, 0.2670, 0.2738, 0.2805, 0.2869, &
!         0.2927, 0.2979, 0.3024, 0.3059, 0.3085, 0.3101, &
!         0.3107, 0.3102, 0.3087, 0.3061, 0.3026, 0.2983, &
!         0.2931, 0.2871, 0.2806, 0.2735, 0.2659, 0.2579, &
!         0.2497, 0.2412, 0.2326, 0.2240, 0.2153, 0.2066, &
!         0.1979, 0.1894, 0.1809, 0.1726, 0.1643, 0.1562, &
!         0.1482, 0.1404, 0.1326, 0.1250, 0.1175, 0.1101, &
!         0.1028, 0.0956, 0.0886, 0.0818, 0.0751, 0.0687 /)       
       real,    parameter :: NonBoreal(N) = &
       (/ 0.0121, 0.0150, 0.0172, 0.0185, 0.0189, 0.0184, &
          0.0174, 0.0162, 0.0151, 0.0141, 0.0133, 0.0126, &
          0.0121, 0.0117, 0.0115, 0.0114, 0.0114, 0.0116, &
          0.0120, 0.0126, 0.0133, 0.0142, 0.0151, 0.0159, &
          0.0167, 0.0174, 0.0180, 0.0184, 0.0187, 0.0189, &
          0.0190, 0.0190, 0.0191, 0.0192, 0.0192, 0.0193, &
          0.0194, 0.0194, 0.0193, 0.0192, 0.0190, 0.0187, &
          0.0185, 0.0182, 0.0180, 0.0178, 0.0177, 0.0176, &
          0.0174, 0.0172, 0.0169, 0.0166, 0.0162, 0.0158, &
          0.0153, 0.0149, 0.0144, 0.0138, 0.0132, 0.0126, &
          0.0118, 0.0109, 0.0101, 0.0092, 0.0085, 0.0081, &
          0.0080, 0.0083, 0.0091, 0.0102, 0.0117, 0.0135, &
          0.0157, 0.0182, 0.0210, 0.0240, 0.0273, 0.0308, &
          0.0345, 0.0387, 0.0432, 0.0483, 0.0540, 0.0606, &
          0.0683, 0.0775, 0.0886, 0.1022, 0.1188, 0.1388, &
          0.1625, 0.1905, 0.2229, 0.2602, 0.3025, 0.3500, &
          0.4031, 0.4623, 0.5283, 0.6016, 0.6824, 0.7705, &
          0.8650, 0.9646, 1.0676, 1.1713, 1.2722, 1.3662, &
          1.4491, 1.5174, 1.5685, 1.6014, 1.6173, 1.6200, &
          1.6150, 1.6082, 1.6040, 1.6058, 1.6157, 1.6353, &
          1.6651, 1.7045, 1.7513, 1.8024, 1.8541, 1.9022, &
          1.9429, 1.9738, 1.9947, 2.0072, 2.0132, 2.0141, &
          2.0096, 1.9994, 1.9829, 1.9604, 1.9321, 1.8977, &
          1.8562, 1.8052, 1.7419, 1.6646, 1.5738, 1.4734, &
          1.3693, 1.2676, 1.1724, 1.0851, 1.0052, 0.9317, &
          0.8637, 0.8004, 0.7414, 0.6862, 0.6348, 0.5871, &
          0.5434, 0.5037, 0.4682, 0.4368, 0.4097, 0.3864, &
          0.3667, 0.3499, 0.3355, 0.3231, 0.3123, 0.3029, &
          0.2944, 0.2862, 0.2773, 0.2670, 0.2547, 0.2402, &
          0.2238, 0.2061, 0.1882, 0.1712, 0.1562, 0.1434, &
          0.1332, 0.1251, 0.1189, 0.1141, 0.1103, 0.1071, &
          0.1043, 0.1018, 0.0996, 0.0979, 0.0968, 0.0964, &
          0.0966, 0.0970, 0.0973, 0.0970, 0.0959, 0.0938, &
          0.0909, 0.0873, 0.0831, 0.0784, 0.0732, 0.0676, &
          0.0618, 0.0565, 0.0521, 0.0491, 0.0475, 0.0473, &
          0.0480, 0.0492, 0.0504, 0.0514, 0.0519, 0.0521, &
          0.0520, 0.0517, 0.0513, 0.0510, 0.0507, 0.0507, &
          0.0508, 0.0512, 0.0515, 0.0518, 0.0519, 0.0518, &
          0.0513, 0.0506, 0.0496, 0.0482, 0.0465, 0.0443, &
          0.0418, 0.0387, 0.0351, 0.0310, 0.0263, 0.0214 /)

!      Fixed normalization factors; a more accurate normalization would take
!      in consideration longitude and time step
!      ---------------------------------------------------------------------
       real*8, save :: fBoreal = -1., fNonBoreal = -1
       real,   save :: fDT=-1

       integer :: hh, mm, ss, ndt, i, j, k
       integer :: NN
       real :: secs, secs_local, aBoreal, aNonBoreal, alpha

!                              -----

!      Normalization factor depends on timestep
!      ----------------------------------------
       if ( fDT /= cdt ) then
            fBoreal = 0.0
            fNonBoreal = 0.0
            NN = 0
            ndt = max(1,nint(cdt/DT))

            do k = 1, N, ndt
               NN = NN + 1
               fBoreal    = fBoreal    + Boreal(k)
               fNonBoreal = fNonBoreal + NonBoreal(k)
            end do

            fBoreal    = fBoreal / NN
            fnonBoreal = fnonBoreal / NN
            fDT = cdt ! so it recalculates only if necessary
       end if


!      Find number of secs since begining of the day (GMT)
!      ---------------------------------------------------
       hh = nhms/10000
       mm = (nhms - 10000*hh) / 100
       ss = nhms - 10000*hh - 100*mm
       secs = 3600.*hh + 60.*mm + ss

!      Apply factors depending on latitude
!      -----------------------------------
       do j = lbound(Ein,2), ubound(Ein,2)
         do i = lbound(Ein,1), ubound(Ein,1)

!            Find corresponding index in hardwired diurnal cycle
!            240 = 24 * 60 * 60 secs / 360 deg
!            ---------------------------------------------------
             secs_local = secs + 240. * lons(i,j)
             k = 1 + mod(nint(secs_local/DT),N)
             if ( k < 1 ) k = N + k

!            Apply diurnal cycle
!            -------------------
             aBoreal = Boreal(k) / fBoreal 
             aNonBoreal = NonBoreal(k) / fNonBoreal

                if ( lats(i,j) >= 50. ) then
                   Eout(i,j) = aBoreal    * Ein(i,j)
                else if ( lats(i,j) >= 30. ) then
                   alpha = (lats(i,j) - 30. ) / 20.
                   Eout(i,j) = (1-alpha) * aNonBoreal * Ein(i,j) + &
                                  alpha  * aBoreal    * Ein(i,j)
                else                  
                   Eout(i,j) = aNonBoreal * Ein(i,j)
                end if
          end do
       end do

     end subroutine Chem_BiomassDiurnal


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Chem_UtilResVal --- returns resolution dependent value
!
! !INTERFACE:
!
   function Chem_UtilResVal( im_World, jm_World, res_value, rc ) result (val)

! !USES:

   implicit NONE

   real :: val                                ! resolution dependent value

! !INPUT/OUTPUT PARAMETERS:
   integer, intent(in) :: im_World, jm_World  ! number of global grid cells
   real,    intent(in) :: res_value(:)        ! array with the resolution dependent values:
                                              ! the 'a', 'b', ..., 'e' resolution values have 
                                              ! indexes 1, 2, ..., 5.

! !OUTPUT PARAMETERS:
   integer, intent(inout) :: rc               ! return code


! !DESCRIPTION: 
!
! !REVISION HISTORY:
!
! 13 Feb2012   Anton Darmenov  First crack.
! 25 Oct2012   Anton Darmenov  Added support for FV3 resolutions.
!
!EOP
!-------------------------------------------------------------------------
       character(len=*), parameter :: Iam = 'Chem_UtilResVal'

       integer            :: i_res       

       integer, parameter :: res_a = 1  ! 'a' to 'e' resolution indexes
       integer, parameter :: res_b = 2  !
       integer, parameter :: res_c = 3  !
       integer, parameter :: res_d = 4  !
       integer, parameter :: res_e = 5  !
       integer, parameter :: res_f = 6  !

       i_res = 0

       if ((im_World < 1) .or. (jm_World < 1)) then
           call die(Iam, 'incorrect model resolution')
       end if

       if (jm_World == 6*im_World) then
           if (im_World <= 24) then
               i_res = res_a
           else if (im_World <=  48) then
               i_res = res_b
           else if (im_World <=  90) then
               i_res = res_c
           else if (im_World <= 180) then
               i_res = res_d
           else if (im_World <= 360) then
               i_res = res_e
           else if (im_World <= 720) then
               i_res = res_f
           else
               i_res = res_f
           end if
       else
           if ((im_World <= 72) .and. (jm_World <= 46)) then
               i_res = res_a
           else if ((im_World <=  144) .and. (jm_World <=  91)) then
               i_res = res_b
           else if ((im_World <=  288) .and. (jm_World <= 181)) then
               i_res = res_c
           else if ((im_World <=  576) .and. (jm_World <= 361)) then
               i_res = res_d
           else if ((im_World <= 1152) .and. (jm_World <= 721)) then
               i_res = res_e
           else if ((im_World <= 2304) .and. (jm_World <=1441)) then
               i_res = res_f
           else
               i_res = res_f
           end if
       end if 

       if ((i_res < 1) .or. (i_res > size(res_value))) then
           val = 0.0
           rc  = 42
       else
           val = res_value(i_res)
           rc  = 0
       end if

   end function Chem_UtilResVal

   function Chem_UtilIdow(nymd) result (idow)
     implicit NONE
     integer, intent(in) :: nymd
     integer :: idow ! day of the week: Sun=1, Mon=2, etc.
     integer :: y, m, d
     integer, parameter :: t(0:11) = (/ 0, 3, 2, 5, 0, 3, 5, 1, 4, 6, 2, 4 /)
     y = nymd / 10000
     m = (nymd - y*10000)/100
     d = nymd - (y*10000 + m*100)
     if ( m<3 ) then
        y = y - 1
     end if
     idow = 1+mod(y + y/4 - y/100 + y/400 + t(m-1) + d,7)
     return 
   end function Chem_UtilIdow

   function Chem_UtilCdow(nymd) result (cdow)
     implicit NONE
     integer, intent(in) :: nymd
     character(len=3) :: cdow ! day of the week: Sun, Mon, etc.
     character(len=3) :: cday(7) = (/ 'Sun','Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat' /)
     cdow = cday(Chem_UtilIdow(nymd))
     return 
   end function Chem_UtilCdow


!  Chem_UtilPointEmissions
!  Colarco, February 9, 2015
!  Given a text file of point wise emissions (see example) return a table with
!  the emissions location (lat,lon), altitude (bottom, top), amount (kg/s over
!  duration of event), and (optionally) the NHMS start and end times of
!  emissions.  Useful for events, like volcanic eruptions or individual fires.
!  This is inspired by GetVolcDailyTables in SulfateChemDriverMod.F90
!
!  Table format
!  Here's an example from volcanic tables (remove ! from lines for functionality
!###  latitude (-90,90), longitude (-180,180), amount [kg/s], 
!###  base elevation [m], top altitude [m], 
!###  (optional) begin time [HHMMSS], (optional) end time [HHMMSS]
!source::
!50.170 6.850 3.587963e-03 600. 600.
!::
!
!  Arguments
!  Input
!   nymd         -- integer YYYYMMDD for file
!   filetemplate -- grads-like filename template filled in with nymd
!  Output
!   nPts         -- number of events in file
!   vLat         -- latitude (one per event...)
!   vLon         -- longitude
!   vBase        -- base altitude (e.g., bottom of plume)
!   vTop         -- top altitude (e.g., top of plume)
!   vEmis        -- emission flux (e.g., kg s-1 of species)
!   vStart       -- HHMMSS to start emissions (optional)
!   vEnd         -- HHMMSS to end emissions (optional)

   subroutine Chem_UtilPointEmissions( nymd, filetemplate, &
                                       nPts, vLat, vLon, vBase, vTop, vEmis, vStart, vEnd )
			  
  implicit NONE

  integer, intent(in)            :: nymd
  character(len=255)             :: filetemplate
  integer                        :: nPts
  real, pointer, dimension(:)    :: vLat, vLon, vTop, vBase, vEmis
  integer, pointer, dimension(:) :: vStart, vEnd
  integer :: i, j, nLines, nCols, rc, STATUS, nymd1, nhms1, ios
  character(len=255) :: fname
  type(ESMF_Config)  :: cf
  real, pointer, dimension(:) :: vData

! If previous instance of volcano point data tables exist, deallocate it
! to get the correct number of elements
  if(associated(vLat))    deallocate(vLat, stat=ios)
  if(associated(vLon))    deallocate(vLon, stat=ios)
  if(associated(vEmis))   deallocate(vEmis, stat=ios)
  if(associated(vBase))   deallocate(vBase, stat=ios)
  if(associated(vTop))    deallocate(vTop, stat=ios)
  if(associated(vStart))  deallocate(vStart, stat=ios)
  if(associated(vEnd))    deallocate(vEnd, stat=ios)


! Assumes files provided daily (or less frequently)
! -------------------------------------------------
  nymd1 = nymd
  nhms1 = 120000
  call StrTemplate ( fname, filetemplate, xid='unknown', &
                     nymd=nymd1, nhms=nhms1 )
  cf = ESMF_ConfigCreate()
  call ESMF_ConfigLoadFile(cf, fileName=trim(fname), rc=STATUS )
  call ESMF_ConfigGetDim(cf, nLines, nCols, LABEL='source::', rc=STATUS )
  nPts = nLines
  allocate(vData(nCols), vLat(nLines), vLon(nLines), &
           vEmis(nLines), vBase(nLines), vStart(nLines), &
           vEnd(nLines), vTop(nLines), stat=ios)
  vStart = -1
  vEnd   = -1
  call ESMF_ConfigFindLabel(cf, 'source::',rc=STATUS)
     do i = 1, nLines
      call ESMF_ConfigNextLine(cf, rc=rc)
      do j = 1, nCols
       call ESMF_ConfigGetAttribute(cf, vData(j), default=-1.)
      end do
      vLat(i)    = vData(1)
      vLon(i)    = vData(2)
      vEmis(i)   = vData(3)
      vBase(i)   = vData(4)
      vTop(i)    = vData(5)
      if(nCols >= 6) vStart(i)  = vData(6)
      if(nCols >= 7) vEnd(i)    = vData(7)
  end do
! Check value of vStart and vEnd.  Set to be
! vStart = 000000 if default (=-1) is provided
! vEnd   = 240000 if default (=-1) is provided
  where(vStart < 0) vStart = 000000
  where(vEnd < 0)   vEnd   = 240000

  call ESMF_ConfigDestroy(cf)
  deallocate(vData, stat=ios)

  end subroutine Chem_UtilPointEmissions

 end module Chem_UtilMod

