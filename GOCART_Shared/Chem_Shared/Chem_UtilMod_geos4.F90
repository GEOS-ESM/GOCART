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

#if defined(SPMD)
   use mod_comm, only: gid           ! fvGCM communication library
#ifdef FORCE_R4
   use mod_comm, only: mp_scatter2d => mp_scatter2d_r4 ! to cope lf95 bug 
   use mod_comm, only: mp_scatter4d => mp_scatter4d_r4
#else
   use mod_comm, only: mp_scatter2d  ! fvGCM communication library
   use mod_comm, only: mp_scatter4d  ! fvGCM communication library
#endif
#endif
   use Chem_Mod           ! Chemistry Base Class
   use mod_diag           ! fvGCM diagnostics
   use m_die
   use m_StrTemplate

   implicit NONE

!
! !PUBLIIC MEMBER FUNCTIONS:
!

   PRIVATE
   PUBLIC  Chem_UtilMPread

!
! !DESCRIPTION:
!
!  This module implements assorted odds & ends for fvChem.
!
! !REVISION HISTORY:
!
!  29oct2003 da Silva  First crack.
!
!EOP
!-------------------------------------------------------------------------

CONTAINS

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Chem_UtilMPread --- Reads fields from file and distribute
!
! !INTERFACE:
!
   subroutine Chem_UtilMPread ( filen, varn, nymd, nhms, &
                                i1, i2, ig, im, j1, j2, jg, jm, km, &
! ++PRC
                                grid, &
! --PRC
                                var2d, var3d, cyclic )

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
  integer, OPTIONAL, intent(in) :: grid   ! not need in GEOS-4


! !OUTPUT PARAMETERS:

  real, OPTIONAL, intent(out)   :: var2d(i1-ig:i2+ig,j1-jg:j2+jg)
  real, OPTIONAL, intent(out)   :: var3d(i1-ig:i2+ig,j1-jg:j2+jg,km)
  logical, OPTIONAL, intent(in) :: cyclic ! whether time dimension is periodic

! !DESCRIPTION: 
!
! !REVISION HISTORY:
!
!  28oct2003 da Silva  First crack.
!  03ayg2004 da Silva  Uses GetVarT for time interpolation
!  18nov2004 da Silva  Added cyclic option for climatological files.
!  31may2005 da Silva  Template expansion.
!
!EOP
!-------------------------------------------------------------------------

    character(len=*), parameter ::  myname = 'Chem_UtilMPread'
    logical :: tcyclic

    integer  :: READ_ONLY=1, nokm=0
    real, allocatable :: v2d(:,:), v3d(:,:,:)
    integer :: fid, rc, ios
    character(len=255) :: fname

#if !defined (SPMD)
    integer :: gid = 0
#endif

!   Consistency check
!   -----------------
    if ( .not. ( present(var2d) .or. present(var3d) ) ) then
       call die ( myname, 'missing var2d or var3d' )
    else if ( present(var2d) .and. present(var3d) ) then
       call die ( myname, 'either var2d or var3d, but not both' )
    end if
    if ( i1 /=1 .or. i2 /=im .or. ig /= 0 ) &
       call die ( myname, 'fvgcm only allows distributed latitudes' )

    if ( present(cyclic) ) then
        tcyclic = cyclic
    else
        tcyclic = .false. ! by default time dimension is not periodic
    end if

!   Expand templates
!   ----------------
    if ( index(filen,'%') .gt. 0 ) then
         call StrTemplate ( fname, filen, xid='unknown', &
                            nymd=nymd, nhms=nhms )
    else
         fname = filen
    end if


!   Read file
!   ---------
    if ( gid .eq. 0 ) then

#ifdef DEBUG
       print *, myname // ': Opening GFIO file ' // trim(fname)
#endif

!      Open file, get first time
!      -------------------------
       call GFIO_Open ( fname, READ_ONLY, fid, rc )
       if ( rc .ne. 0 ) then
          call die(myname,'cannot open GFIO file '//trim(fname))
       end if

!      Read global array
!      -----------------
#if defined(SPMD)
       if ( present(var2d) ) then
#ifdef DEBUG
            print *, myname // ':    reading variable ' // trim(varn)
#endif
            allocate(v2d(im,jm),stat=ios)
            if ( ios /= 0 ) call die ( myname, 'cannot allocate v2d' )
            call GFIO_GetVarT1 ( fid, trim(varn), nymd, nhms, im, jm, nokm, 1, &
                               v2d, rc, tcyclic, fid )
            if ( rc .ne. 0 ) call die(myname,'cannot read '//trim(varn) )
       else
            allocate(v3d(im,jm,km),stat=ios)
            if ( ios /= 0 ) call die ( myname, 'cannot allocate v3d' )
            call GFIO_GetVarT1 ( fid, trim(varn), nymd, nhms, im, jm, 1, km, &
                               v3d, rc, tcyclic, fid )
            if ( rc .ne. 0 ) call die(myname,'cannot read '//trim(varn) )
       end if
#else
       if ( present(var2d) ) then
#ifdef DEBUG
            print *, myname // ':    reading variable ' // trim(varn)
#endif
            call GFIO_GetVarT1 ( fid, trim(varn), nymd, nhms, im, jm, nokm, 1, &
                               var2d, rc, tcyclic, fid )
            if ( rc .ne. 0 ) call die(myname,'cannot read '//trim(varn) )
       else
            if ( ios /= 0 ) call die ( myname, 'cannot allocate v3d' )
            call GFIO_GetVarT1 ( fid, trim(varn), nymd, nhms, im, jm, 1, km, &
                               var3d, rc, tcyclic, fid )
            if ( rc .ne. 0 ) call die(myname,'cannot read '//trim(varn) )
       end if
#endif

!      Close file
!      ----------
       call GFIO_Close ( fid, rc )
#ifdef DEBUG
       print *, myname // ': Closing GFIO file ' // trim(fname)
#endif

    end if ! masterproc

!   Distribute data
!   ---------------
#if defined(SPMD)
    if ( present(var2d) ) then
         call mp_scatter2d(v2d, var2d, im,  jm, j1, j2, 0) 
         if ( gid == 0 ) then 
              deallocate(v2d,stat=ios)
              if ( ios /= 0 ) call die ( myname, 'cannot deallocate v2d' )
         end if
    else
         call mp_scatter4d ( v3d, var3d, im,  jm, km, 1, j1, j2,   &
                             1, km, jg, jg, 0)
         if ( gid == 0 ) then
              deallocate(v3d,stat=ios)
              if ( ios /= 0 ) call die ( myname, 'cannot deallocate v3d' )
         end if
    end if
#endif

!   All done
!   --------
    return

end subroutine Chem_UtilMPread

 end module Chem_UtilMod

