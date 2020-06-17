#include "unused_dummy.H"
!-------------------------------------------------------------------------
!      NASA/GSFC, Global Modeling & Assimilation Office, Code 900.3      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  Chem_BundleMod.F90 --- Implements Chemical Bundle Class
! 
! !INTERFACE:
!
      MODULE  Chem_BundleMod
            
! !USES:

      Use ESMF
      Use Chem_RegistryMod
      Use Chem_ArrayMod
      Use m_chars, only: uppercase
      Use MAPL, only: MAPL_UNDEF

      Implicit NONE

! !PUBLIC TYPES:
!
      PRIVATE 
      PUBLIC  Chem_Bundle           ! chemical bundle type
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  Chem_Grid             ! grid definition
      PUBLIC  Chem_BundleCreate     ! initializes chemical bundle
      PUBLIC  Chem_BundleDestroy    ! "cleans" chemical bundle
      PUBLIC  Chem_BundleWrite      ! writes Chem_Bundle to  file
      PUBLIC  Chem_BundleRead       ! reads Chem_Bundle from file
      PUBLIC  Chem_BundleStat       ! Prints vital stats by level
      PUBLIC  Chem_BundleSetPtr     ! Fix internal pointers

!
! !DESCRIPTION: This module implements the chemical bundle class.
!               It is loosely based on the fvDAS dynamics class m\_dyn.
!
! !REVISION HISTORY: 
!
! 18sep2001  da Silva  Initial code.
! 29jul2003  da Silva  Revised for ESMF consistency
! 16oct2003  da Silva  Introduced ghosting
! 15jan2003  da Silva  Removed ptop,delp,q as arguments; added skipAlloc
! 29Mar2005  da Silva  Revised RH implementation: before it was always
!                      allocated, not it is treated as the other variables.
! 11Aug2005  da Silva  Added optional levels during creation
! 12Aug2005  da Silva  Fixed bug in SetPtr routne (now it checks whether
!                      array is associated before pointing to it).
! 15Jan2008  Nielsen   Added lon_min=-180.0 for GEOS5.
! 16Jun2016  Buchard   Added option to do concentration instead of MR if needed (read airdens) 
!EOP
!-------------------------------------------------------------------------

   real, parameter ::  missing_val = MAPL_UNDEF ! hardwire this for now
   integer, parameter :: nch = 256

!   Grid
!   ----
    type Chem_Grid

!     Zonal grid
!     ----------
      integer       :: i1, i2, iml               ! local indices
      integer       :: ig                        ! ghosting
      integer       :: im                        ! global dimension
      integer       :: iLeft                     ! i1's index on global grid
      real, pointer :: lon(:,:) => null()        ! longitudes (deg)

!     Meridional grid
!     ---------------
      integer       :: j1, j2, jml               ! local indices
      integer       :: jg                        ! ghosting
      integer       :: jm                        ! global dimension
      real, pointer :: lat(:,:) => null()        ! latitudes (deg)

!     Vertical grid
!     -------------
      integer       :: km
      real, pointer :: lev(:) => null()
      character(len=nch) :: levUnits 
      real          :: ptop          ! Top pressure [Pa]

!     Horizontal gridbox area
!     -----------------------
      real, pointer :: cell_area(:,:) => null()

!     Cubed sphere or not
!     -------------------
      logical :: Cubed_Sphere = .FALSE.

    end type Chem_Grid

!   Chemical vector
!   ---------------
    type Chem_Bundle

!     Chem Registry
!     -------------
      type(Chem_Registry) :: reg

!     Grid
!     ----
      type(Chem_Grid) :: grid

      real, pointer   :: cosz(:,:) => null()  ! cosine solar zenith angle
      real, pointer   :: sinz(:,:) => null()  !   sine solar zenith angle

      type(ESMF_Grid) :: grid_esmf

!     Whether this class allocated the memory for q, delp
!     ---------------------------------------------------
      logical :: did_allocation = .false.
      logical :: has_rh = .false.            ! for backward compatibility
      logical :: diurnal_bb = .false.        ! whether using diurnal biomass burning
      logical :: do_concentration = .false.  ! using concentration: MR * airdens instead of MR alone
      real    :: missing_value = MAPL_UNDEF

!     Tracer array
!     ------------
      real, pointer :: delp(:,:,:) => null()! Layer thickness [Pa] (not ghosted)
      real, pointer :: rh(:,:,:) => null()  ! Layer thickness [Pa] (not ghosted)
      real, pointer :: airdens(:,:,:) => null() ! Air density

      type(chem_array), pointer :: qa(:) => null()
                                            ! access 4D array in q as a 
                                            ! collection of 3D arrays; used
                                            ! for gradually removing the 4D
                                            ! arrays     

!     Two calendar elements (from ESMF)
!     ---------------------------------
      LOGICAL :: isLeapYear
      REAL :: dayOfYear

    end type Chem_Bundle

!  Interfaces
!  ----------
   PUBLIC  Chem_BundleCreate_    ! initializes chemical bundle
   PUBLIC  Chem_BundleCreate1PE_ ! initializes chemical bundle
   interface Chem_BundleCreate
      module procedure Chem_BundleCreate_    ! Distributed data
      module procedure Chem_BundleCreate1PE_ ! Undistributed data (1PE)
   end interface
 
   CONTAINS

!-------------------------------------------------------------------------
!      NASA/GSFC Global Modeling & Assimilation Office, Code 900.3       !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Chem_BundleCreate_ --- Creates Chemical Bundle
! 
! !INTERFACE:
!
  subroutine  Chem_BundleCreate_ ( reg,                 &
                                   i1, i2, ig, im,      &
                                   j1, j2, jg, jm, km,  &
                                   w_c, rc,             &
                                   skipAlloc, lat, lon, &
                                   cell_area,           &
                                   lev, levUnits, ptop, &
                                   do_Conc )  ! Optional
!
! !USES:
!
  implicit NONE
!
! !INPUT PARAMETERS: 
!

  type(Chem_registry)        :: reg             ! Chemical Registry

                                                ! Distributed grid info:
  integer,      intent(in)   :: i1, i2          !   local  zonal indices
  integer,      intent(in)   :: ig              !   zonal ghosting
  integer,      intent(in)   :: im              !   global zonal dimension
  integer,      intent(in)   :: j1, j2          !   local meridional indices
  integer,      intent(in)   :: jg              !   meridional ghosting
  integer,      intent(in)   :: jm              !   global zonal dimension
  integer,      intent(in)   :: km              !   vertical dimension

  logical, OPTIONAL, intent(in) :: skipAlloc    ! Do not allocate arrays
  real,    OPTIONAL, intent(in) :: lon(i1:i2,j1:j2) ! longitude in degrees
  real,    OPTIONAL, intent(in) :: lat(i1:i2,j1:j2) ! latitude in degrees
  real,    OPTIONAL, pointer    :: cell_area(:,:) ! grid box area
  real,    OPTIONAL, intent(in) :: lev(1:km)    ! levels
  character(len=*), OPTIONAL, intent(in) :: levUnits ! level units
  real,    OPTIONAL, intent(in) :: ptop         ! top pressure in Pa
  logical, OPTIONAL, intent(in) :: do_Conc ! Do_concentration: In case you need MR*airdens instead of MR alone
!
! !OUTPUT PARAMETERS:
!
  type(Chem_Bundle), intent (out) :: w_c   ! Chemical bundle

  integer, intent(out)            :: rc    ! error return code:
                                           !  0 - all is well
                                           !  1 - already allocated
                                           !  2 - could not allocate memory
                                           !  3 - invalid dimensions

!
! !DESCRIPTION: Creates a Chemical Bundle, allocating the necessary memory or
!               optionally associating internal pointers with model declared 
!  flat arrays.
!
! !REVISION HISTORY: 
!
!  20Sep2001 da Silva  Initial code based on m_dyn.
!  05Sep2003 da Silva  Revised ESMF-like design.
!  14mar2003 da Silva  Added rh, lat, lon 
!  16Jun2016 Buchard   Added do_concentration option
!EOP
!-------------------------------------------------------------------------

     character(len=*), parameter ::  myname = 'Chem_BundleCreate_'

     integer err, i, j, n, nq, ios, ios1, ios2, ios3, ios4
     logical :: do_allocation
     logical :: do_concentration
     real*8 :: delta

!    Sanity check
!    ------------
     rc = 0
     nq = reg%nq
     if ( im<1 .or. jm<1 .or. km<1 .or. nq<1) then
          rc = -3
          return
     endif
              
     w_c%reg = reg
     w_c%missing_value = MAPL_UNDEF

!    Whether or not we allocate memory for arrays
!    --------------------------------------------
     if ( present(skipAlloc) ) then
          do_allocation = .not. skipAlloc
     else
          do_allocation = .true.
     end if

!   Whether or not read airdens to be able to work with concentration instead  of MR
!   -----------------------------------------------
     if ( present(do_Conc) ) then
          do_concentration = do_Conc
     else
          do_concentration = .false.
     end if



!    Initialize dimensional attributes
!    ---------------------------------
     w_c%grid%i1 = i1; w_c%grid%i2 = i2; w_c%grid%ig = ig; w_c%grid%im = im
     w_c%grid%iml = i2 - i1 + 1 
     w_c%grid%j1 = j1; w_c%grid%j2 = j2; w_c%grid%jg = jg; w_c%grid%jm = jm
     w_c%grid%jml = j2 - j1 + 1 
     w_c%grid%km = km

!    Detect cubed sphere for sanity checks latter
!    --------------------------------------------
     if ( jm == im * 6 ) then
          w_c%grid%Cubed_Sphere = .TRUE.
     else
          w_c%grid%Cubed_Sphere = .FALSE.
     end if

!    Horizontal grid (hardwire A-grid for now)
!    -----------------------------------------
     if ( present(ptop) ) then
          w_c%grid%ptop =  ptop
     else
          w_c%grid%ptop =  1.0 ! Pa: reasonable default 
     endif

!    Save lat/lons
!    -------------
     allocate ( w_c%grid%lon(i1:i2,j1:j2), w_c%grid%lat(i1:i2,j1:j2), &
                stat = ios ) ! 
     if ( ios /= 0 ) then
        rc = 2
        return
     end if
     if ( present(lon) ) then
          w_c%grid%lon = lon
     else
          !ALT w_c%grid%lon = MAPL_UNDEF
          delta = 360.0d0/im
          do i = 1, im
            w_c%grid%lon(i,:) = -180.0d0 + (i-1)*delta
          end do
     end if

     if ( present(lat) ) then
          w_c%grid%lat = lat
     else
          !ALT w_c%grid%lat = MAPL_UNDEF
          if(jm==1) then
            delta = 0.0d0
          else
            delta = 180.0d0/(jm-1)
          endif
          do j = 1, jm
            w_c%grid%lat(:,j) = -90.0d0 + (j-1)*delta
          end do
     end if

     if ( present(cell_area) ) then ! will be left unallocated otherwise
          w_c%grid%cell_area => cell_area
     end if

     if ( present(lev) ) then
        allocate ( w_c%grid%lev(1:km), stat = ios )
        if ( ios /= 0 ) then
           rc = 3
           return
        end if
        w_c%grid%lev = lev 
     end if
     if ( present(levUnits) ) then
          w_c%grid%levUnits = levUnits
     else
          w_c%grid%levUnits = 'none'
     end if

     w_c%did_allocation = .false.
     w_c%has_rh = .false.       ! will be set to TRUE when set
     w_c%do_concentration = .false.

     allocate(w_c%qa(nq),stat=ios2)

     if ( do_allocation ) then

        allocate(w_c%delp(i1:i2,j1:j2,km),stat=ios )
        do n = 1, nq
           allocate(w_c%qa(n)%data3d(i1-ig:i2+ig,j1-jg:j2+jg,km),stat=ios1 )
           ios2 = ios2 + ios1
        end do
        allocate(w_c%rh(i1:i2,j1:j2,km),stat=ios3 )
            
        if ( ios + ios2 + ios3 == 0 ) then 
             w_c%did_allocation = .true.
             w_c%delp = 0.0
             do n = 1, nq
                w_c%qa(n)%data3d = 0.0
             end do
             w_c%rh = 0.0             
        else
             w_c%did_allocation = .false.
             rc = 4
        end if

        if (do_concentration ) then
           w_c%do_concentration = .true.
           allocate(w_c%airdens(i1:i2,j1:j2,km),stat=ios4 )
           if ( ios4 == 0 ) w_c%airdens = 0.0
        endif

         
     end if

!    Set array of pointers: may be null() if no allocation took place
!    ----------------------------------------------------------------
     call Chem_BundleSetPtr ( w_c, rc ) 

  
   end subroutine Chem_BundleCreate_

!-------------------------------------------------------------------------
!      NASA/GSFC Global Modeling & Assimilation Office, Code 900.3       !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Chem_BundleCreate1PE_ --- Creates Chemical Bundle (1 PE)
! 
! !INTERFACE:
!
  subroutine  Chem_BundleCreate1PE_ ( reg, im, jm, km, &
                                      w_c, rc,             &
                                      skipAlloc, ptop  )  ! Optional
!
! !USES:
!
  implicit NONE
!
! !INPUT PARAMETERS: 
!

  type(Chem_registry)        :: reg             ! Chemical Registry

                                                ! Distributed grid info:
  integer,      intent(in)   :: im              !   global zonal dimension
  integer,      intent(in)   :: jm              !   global zonal dimension
  integer,      intent(in)   :: km              !   vertical dimension



  logical, OPTIONAL, intent(in) :: skipAlloc    ! Do not allocate arrays
  real,    OPTIONAL, intent(in) :: ptop         ! top pressure in Pa

!
! !OUTPUT PARAMETERS:
!
  type(Chem_Bundle), intent (out) :: w_c   ! Chemical bundle
  integer, intent(out)            :: rc    ! error return code:
                                           !  0 - all is well
                                           !  1 - already allocated
                                           !  2 - could not allocate memory
                                           !  3 - invalid dimensions
                                           !  
!
! !DESCRIPTION: Creates a Chemical Bundle, allocating the necessary memory or
!               optionally associating internal pointers with model declared 
!  flat array. All optional arguments must be specified, or NOT at all.
!  This interface is for global, undistributed Bundles.
!
!
! !REVISION HISTORY: 
!
!  20Sep2001 da Silva  Initial code based on m_dyn.
!  05Sep2003 da Silva  Revised ESMF-like design.
!
!EOP
!-------------------------------------------------------------------------

    character(len=*), parameter ::  myname = 'Chem_BundleCreate1PE_'

    integer :: i1, i2, ig, j1, j2, jg
    logical :: skipAllocation = .false.
    real :: p_top = 1.0


    if ( present(skipAlloc) ) then
       skipAllocation = skipAlloc
    end if
    if ( present(ptop) ) then
         p_top = ptop
    end if

    i1 = 1;  i2 = im; ig = 0
    j1 = 1;  j2 = jm; jg = 0 

    call Chem_BundleCreate_ ( reg, i1, i2, ig, im, j1, j2, jg, jm, km, &
                                    w_c, rc, skipAlloc=skipAllocation, &
                                    ptop = p_top )

  end subroutine Chem_BundleCreate1PE_

!-------------------------------------------------------------------------
!      NASA/GSFC Global Modeling & Assimilation Office, Code 900.3       !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Chem_BundleDestroy --- Deallocates memory used by chemical state
! 
! !INTERFACE:
!
  subroutine  Chem_BundleDestroy ( w_c, rc )
!
! !USES:
!
  implicit NONE
!
! !INPUT/OUTPUT PARAMETERS: 
!
  type(Chem_Bundle), intent (inout) :: w_c   ! chemical bundle

! !OUTPUT PARAMETERS:

  integer, intent(out)          ::  rc     ! Error return code:
                                           !  0 - all is well
                                           !  1 - 
! !DESCRIPTION: 
!
!  Deallocates memory used by chemical bundle.
!
! !REVISION HISTORY: 
!
!  20Jul1999 da Silva  Initial code.
!  26oct1999 da Silva  Added hs_stdv, ts, lwi, a, b
!  16Jun2016 Buchard   Added do_concentration option
!EOP
!-------------------------------------------------------------------------

   integer n, ier

   rc = 0

   if ( w_c%did_allocation ) then

      if ( associated(w_c%delp) ) deallocate(w_c%delp, stat=ier)
      if ( associated(w_c%rh)  )  deallocate(w_c%rh, stat=ier)

       
     
    
      do n = 1, w_c%reg%nq 
         if ( associated(w_c%qa(n)%data3d)  ) &
              deallocate(w_c%qa(n)%data3d, stat=ier)
      end do
      if(associated(w_c%grid%lon)) deallocate( w_c%grid%lon )
      if(associated(w_c%grid%lat)) deallocate( w_c%grid%lat )
      if(associated(w_c%grid%lev)) deallocate( w_c%grid%lev )
      if(associated(w_c%qa))       deallocate( w_c%qa )


      if ( w_c%do_concentration ) then
          if ( associated(w_c%airdens)  )  deallocate(w_c%airdens, stat=ier)
      endif

      
   else

      if ( associated(w_c%delp) ) nullify(w_c%delp)
      if ( associated(w_c%rh) )   nullify(w_c%rh)
            
      do n = 1, w_c%reg%nq 
         if ( associated(w_c%qa(n)%data3d)  ) nullify(w_c%qa(n)%data3d)
      end do
      deallocate( w_c%grid%lon, w_c%grid%lat, w_c%grid%lev,w_c%qa, stat=ier)  

      if ( w_c%do_concentration ) then
          if ( associated(w_c%airdens)  )  nullify(w_c%airdens)
      endif
 
   end if


   w_c%reg%nq = -1

  end subroutine Chem_BundleDestroy

!-------------------------------------------------------------------------
!      NASA/GSFC Global Modeling & Assimilation Office, Code 900.3       !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Chem_BundleWrite --- writes out single instance of chemical state
! 
! !INTERFACE:
!
  subroutine  Chem_BundleWrite ( fname, nymd, nhms, prec, w_c, rc, &
                                 verbose, new, freq, gfio_prec )   ! optional
!
! !USES:
!
  implicit NONE

!
! !INPUT PARAMETERS: 
!
  character(len=*),    intent(in)   :: fname  ! output file name
  integer,             intent(in)   :: nymd   ! Date: year-month-day
  integer,             intent(in)   :: nhms   ! Time: hour-min-sec
  integer,             intent(in)   :: prec   ! precision:
                                              ! 0 = 32 bits
                                              ! 1 = 64 bits
  type(Chem_Bundle), intent(in)     :: w_c    ! chemical bundle

  logical, intent(in), OPTIONAL     :: verbose ! if true, send log to stdout
  logical, intent(in), OPTIONAL     :: new     ! create new file even if it
                                               ! already exists.
  integer, intent(in), OPTIONAL     :: freq    ! time frequency (HHMMSS) for
                                               ! multiple instance files
                                               ! (default: 060000)
  integer, OPTIONAL,   intent(in)   :: gfio_prec    ! specify user precision

!
! !OUTPUT PARAMETERS:
!

  integer, intent(out)              :: rc    ! error return code:
                                             !  0 - all is well
                                             !  >0 - errors
!
! !DESCRIPTION: Writes a GFIO file with one or more instances of the
!               chemical tracer bundle {\tt w\_f} valid at a given time
!  to a file named {\tt fname}. The file is created or opened, written to and
!  closed upon completion. 
!               
! !REVISION HISTORY: 
!
!  20Sep2002 da Silva  Initial code based dyn_put()
!  04Nov2015 Todling/Buchard add underlying gfio precision as option
!  16Jun2016 Buchard Even if do_concentration option is True, the file contains only MR
!
!EOP
!-------------------------------------------------------------------------

   character(len=nch)              :: title, source, contact, levunits

   real,    allocatable :: lat(:), lon(:), lev(:)
   character(len=nch), allocatable :: vname(:), vtitle(:), vunits(:)
   real,    allocatable :: valid_range(:,:), packing_range(:,:)
   integer, allocatable :: kmvar(:)

   real(8),allocatable :: lon8(:),lat8(:),lev8(:)
   real(8),allocatable :: valid_range8(:,:), packing_range8(:,:)
   real(8),allocatable,dimension(:,:,:)::rank3
   integer :: i1, i2, j1, j2, im, jm, km, nq, nvars
   real    :: p, ptop
   integer :: i, j, k, timeinc
   integer :: fid, err
   integer :: gfio_prec_


   logical verb, creating, fexists

   integer, parameter :: READ_WRITE = 0

   rc = 0

!  Short hand for dimensions
!  -------------------------
   i1 = w_c%grid%i1; j1 = w_c%grid%j1
   i2 = w_c%grid%i2; j2 = w_c%grid%j2
   im = w_c%grid%im; jm = w_c%grid%jm
   km = w_c%grid%km; nq = w_c%reg%nq
   nvars = nq + 2            ! delp, RH and tracers 

!  Cannot handle cubed sphere
!  --------------------------
   if ( w_c%grid%cubed_sphere ) then
      rc = 1
      return
   end if

!  No chemical tracers, nothing to do
!  ----------------------------------
   if ( nvars .lt. 2 ) return  

!  Requires global arrays for now
!  ------------------------------
   if ( i1 /= 1 .or. i2 /= im .or. j1 /= 1 .or. j2 /= jm ) then
      rc = 1
      return
   end if

!  Handle optional parameters
!  --------------------------
   if ( present(verbose) ) then
        verb = verbose
   else                    
        verb = .false.
   end if

   if ( present(new) ) then
      creating = new
   else
      creating = .false.
   end if

   gfio_prec_ = 32 ! default precision
   if ( present(gfio_prec) ) then
      gfio_prec_= gfio_prec
   end if


! Check whether file exists
! -------------------------
  inquire ( file=trim(fname), exist=fexists )
  if ( .not. fexists ) creating = .true.       ! must create then

!  Allocate local work space
!  -------------------------
   rc = 0
   call init_ ( err )
   if ( err .ne. 0 ) then
        call clean_()
        rc = 2
        return
   end if


!  Create coordinate variables
!  ---------------------------

   if ( gfio_prec_==64 ) then
      lat8 = w_c%grid%lat(1,:)
      lon8 = w_c%grid%lon(:,1)
   else
      lat = w_c%grid%lat(1,:)
      lon = w_c%grid%lon(:,1)
   endif

!  Vertical coordinates: fake something for GrADS sake
!  ---------------------------------------------------
  

  if ( associated(w_c%grid%lev) ) then
      lev =  w_c%grid%lev
      levUnits = w_c%grid%levUnits
  else
      ptop = w_c%grid%ptop
      i = im / 2
      j = jm / 2
      p = ptop + 0.5 * w_c%delp(i,j,1)
      lev(1) = p
      do k = 2, km
         p = p + 0.5 * ( w_c%delp(i,j,k-1) + w_c%delp(i,j,k) )
         lev(k) = p               
      end do
      lev(1:km) = lev(1:km) / 100.
      levunits = 'hPa'
   end if
   if ( gfio_prec_==64 ) then
      lev8 = lev
    endif
   
   kmvar = km

!  Global metadata
!  ---------------
   title   = 'fvCCM Chemical Bundle (Lagrangian Control Coordinates) v2.0'
   source  = 'Data Assimilation Office, NASA/GSFC'
   contact = 'data@gmao.gsfc.nasa.gov'

!  For now, do not exercise packing/valid range feature
!  -----------------------------------------------------
   if ( gfio_prec_== 64 ) then
      do j = 1, nvars
         do i = 1, 2
            valid_range8(i,j) = w_c%missing_value
            packing_range8(i,j) = w_c%missing_value
         end do
      end do
   else
      do j = 1, nvars
         do i = 1, 2
            valid_range(i,j) = w_c%missing_value
            packing_range(i,j) = w_c%missing_value
         end do
      end do
   endif
!  Time attributes
!  ---------------
   if ( present(freq) ) then
      timeinc = freq
   else
      timeinc = 060000      
   end if

!  Create new GFIO file ...
!  ------------------------
   if ( creating ) then
      if (verb) print *, '	[] creating GFIO file ', trim(fname)
        if ( gfio_prec_ == 64 ) then
         call GFIO_Create ( fname, title, source, contact, &
                            real(w_c%missing_value,8),  &
                            im, jm, km, lon8, lat8, lev8, levunits,      & 
                            nymd, nhms, timeinc,                         &
                            nvars, vname, vtitle, vunits,                &
                            kmvar, valid_range8, packing_range8, prec,   &
                            fid, err )
       else
        call GFIO_Create ( fname, title, source, contact, &
                       w_c%missing_value,  &
                       im, jm, km, lon, lat, lev, levunits,         & 
                       nymd, nhms, timeinc,                         &
                       nvars, vname, vtitle, vunits,                &
                       kmvar, valid_range, packing_range, prec,     &
                       fid, err )
       endif
      
!  ... or open existing GFIO file ...
!  ----------------------------------
   else

    if (verb) print *, '	[] opening GFIO file ', trim(fname)

    call GFIO_Open ( fname, READ_WRITE, fid, err )

   end if

    if ( err .ne. 0 ) then
         rc = 3
         call clean_()
         return
    end if

!   Write the data to GFIO file
!   ---------------------------
    if  ( gfio_prec_ == 64 ) then
        if (verb) print *, '	[] writing ', trim(vname(1))
           rank3 = w_c%delp
           call GFIO_PutVar ( fid, vname(1), nymd, nhms, im, jm, 1, km, &
                       rank3, err )
           if ( err .ne. 0 ) rc = 101

        if (verb) print *, '	[] writing ', trim(vname(2))
           rank3 = w_c%rh
           call GFIO_PutVar ( fid, vname(2), nymd, nhms, im, jm, 1, km, &
                       rank3, err )
           if ( err .ne. 0 ) rc = 102

        do i = 3, nvars
        if (verb) print *, '	[] writing ', trim(vname(i))
           rank3 = w_c%qa(i-2)%data3d(:,:,:)
           call GFIO_PutVar ( fid, vname(i), nymd, nhms, im, jm, 1, km, &
                          rank3, err )
           if ( err .ne. 0 ) rc = 100 + i-1

           if (w_c%do_concentration) then
           rank3 = w_c%qa(i-2)%data3d(:,:,:)/w_c%airdens
           call GFIO_PutVar ( fid, vname(i), nymd, nhms, im, jm, 1, km, &
                          rank3, err )
           if ( err .ne. 0 ) rc = 100 + i-1
           endif
        end do
        
    else

        if (verb) print *, '	[] writing ', trim(vname(1))
           call GFIO_PutVar ( fid, vname(1), nymd, nhms, im, jm, 1, km, &
                       w_c%delp, err )
           if ( err .ne. 0 ) rc = 101
   
        if (verb) print *, '	[] writing ', trim(vname(2))
           call GFIO_PutVar ( fid, vname(2), nymd, nhms, im, jm, 1, km, &
                       w_c%rh, err )
           if ( err .ne. 0 ) rc = 102

        do i = 3, nvars
        if (verb) print *, '	[] writing ', trim(vname(i))
           call GFIO_PutVar ( fid, vname(i), nymd, nhms, im, jm, 1, km, &
                          w_c%qa(i-2)%data3d(:,:,:), err )
           if ( err .ne. 0 ) rc = 100 + i-1

           if (w_c%do_concentration) then
           call GFIO_PutVar ( fid, vname(i), nymd, nhms, im, jm, 1, km, &
                          w_c%qa(i-2)%data3d(:,:,:)/w_c%airdens, err )
           if ( err .ne. 0 ) rc = 100 + i-1
           endif 

        end do

    endif

!   Now save vertical grid info as attributes
!   -----------------------------------------
    if ( creating ) then
       call GFIO_PutRealAtt ( fid, 'ptop',   1, w_c%grid%ptop, prec, err )
       if ( err .ne. 0 ) rc = 201
    end if

    if (verb) print *, '	[] closing GFIO file ', trim(fname)
    call GFIO_close ( fid, err )

!  Clean up
!  --------
   call clean_()

!  All done
!  --------
   return

  CONTAINS

     subroutine init_ ( err )       ! allocates local memory
     integer err, err1, err2, err3, err4
     allocate ( lat(jm), lon(im), lev(km),stat=err1 )
     allocate ( kmvar(nvars), vname(nvars), vtitle(nvars), stat=err2 )
     allocate ( vunits(nvars), stat=err3 )
     allocate ( packing_range(2,nvars), valid_range(2,nvars), stat=err4 )
     if (gfio_prec_==64) then
        allocate(lon8(im))
        allocate(lat8(jm))
        allocate(lev8(km))
        allocate(packing_range8(2,nvars))
        allocate(valid_range8(2,nvars))
        allocate(rank3(im,jm,km))
     end if

     err = err1 + err2 + err3 + err4 
     if ( err /= 0 ) return
     vname(1) = 'delp'
     vtitle(1) = 'Pressure Thickness'
     vunits(1) = 'hPa'
     vname(2) = 'rh'
     vtitle(2) = 'Relative Humidity'
     vunits(2) = 'percent'
     do i = 1, nq
        vname(i+2)  = w_c%reg%vname(i)
        vtitle(i+2) = w_c%reg%vtitle(i)
        vunits(i+2) = w_c%reg%vunits(i)
     end do
     end subroutine init_

     subroutine clean_()             ! de-allocates local memory
     deallocate ( lat, lon, lev, kmvar, vname, vtitle, vunits, stat=err )
     deallocate ( valid_range, packing_range, stat=err )
     if (gfio_prec_==64) then
        deallocate(rank3)
        deallocate(lon8)
        deallocate(lat8)
        deallocate(lev8)
        deallocate(packing_range8)
        deallocate(valid_range8)
     endif

     end subroutine clean_

  end Subroutine Chem_BundleWrite

!-------------------------------------------------------------------------
!      NASA/GSFC Global Modeling & Assimilation Office, Code 900.3       !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Chem_BundleRead --- reads single instance of chemical state
! 
! !INTERFACE:
!
  subroutine  Chem_BundleRead ( fname, nymd, nhms, w_c, rc, &
                                timidx, freq, chemReg, gfio_prec, do_Conc )         ! optional
!
! !USES:
!
  implicit NONE
!
! !INPUT PARAMETERS: 
!
  character(len=*),    intent(in)   :: fname  ! output file name
  integer, OPTIONAL,   intent(in)   :: timidx ! time index; by default
                                              ! last time index is returned
  integer, OPTIONAL,   intent(in)   :: gfio_prec   ! specify user precision
                                                   ! 32 or 64 bits
  logical, OPTIONAL,   intent(in)   :: do_Conc     ! True if you want MR*airdens
                                                   
                                              
  type(Chem_Registry), OPTIONAL, intent(in) :: chemReg ! Chemistry Registry
!
! !INPUT/OUTPUT PARAMETERS:
!
  integer,             intent(inout)   :: nymd  ! Date: year-month-day
  integer,             intent(inout)   :: nhms  ! Time: hour-min-sec
                                                ! Note: unless timidx=0,
                                                !       nymd/nhms will be
                                                !       an output parameter.
   

!
! !OUTPUT PARAMETERS:
!
  type(Chem_Bundle),   intent(inout)   :: w_c   ! chemical bundle
  integer,             intent(out)   :: rc    ! error return code:
                                              !   0 - all is well
                                              !  >0 - errors
  integer, OPTIONAL,   intent(out)   :: freq  ! time frequency on file

!
! !DESCRIPTION: This routine reads GFIO files with one or more instances
!               of the chemical bundle {\tt w\_f} from a file named 
!  {\tt fname}. The file is opened, read from and closed upon completion.
!  By default, this routine returns the last time written to the file.
!  The optional parameter {\tt timidx} allows the retrieval of times
!  other than the last; use routine {\tt GFIO\_DimInquire()} to determine
!  how many times have been written to the file. When {\tt timidx=0}
!  is specified, (nymd,nhms) is used to specify the actual date required
!  (this is the only instance when nymd/nhms are input parameters).
!  When the OPTIONAL parameter chemReg only tracers in specified in it
!  are read. When chemReg it is not specified, the chemical registry is
!  read from file.
!
! !REVISION HISTORY: 
!
!  20Sep2001 da Silva  Initial code.
!  04Nov2015 Todling/Buchard - add underlying gfio precision as option
!                              all bug fix in lat/lon handle to BundleCreate
!
!EOP
!-------------------------------------------------------------------------

   character(len=nch)              :: title, source, contact, levunits
   character(len=nch), allocatable :: vname(:), vtitle(:), vunits(:)
   character(len=nch) :: vname_

   real,    allocatable :: lat(:), lon(:), lev(:)
   real,    allocatable :: lat2d(:,:), lon2d(:,:)
   real,    allocatable :: valid_range(:,:), packing_range(:,:)
   integer, allocatable :: kmvar(:), yyyymmdd(:), hhmmss(:), ivar(:)

   integer :: im, jm, km, lm, nvars
   integer :: i1, i2, ig, j1, j2, jg
   integer :: l, timinc, i, j, n
   real    :: amiss, buf(1), buf8(1)

   integer, parameter :: READ_ONLY = 1
   integer :: gfio_prec_
   integer :: fid, err, ngatts

   real(8) amiss8
   real(8),allocatable :: lon8(:),lat8(:),lev8(:)
   real(8),allocatable :: valid_range8(:,:), packing_range8(:,:)
   real(8),allocatable,dimension(:,:,:)::rank3
   type(Chem_registry) :: Reg
   logical :: all_upper ! whether all variables are upper case

   rc = 0
   gfio_prec_ = 32 ! default precision
   if ( present(gfio_prec) ) then
      gfio_prec_= gfio_prec
   end if

!  Open the file
!  -------------
   call GFIO_Open ( fname, READ_ONLY, fid, err )
   if ( err .ne. 0 ) then
      rc = 1
      call clean_()
      return
   end if

!  Get dimensions
!  --------------
   call GFIO_DimInquire ( fid, im, jm, km, lm, nvars, ngatts, err)
   if ( err .ne. 0 ) then
      call clean_()
      rc = 2
   end if   
   call init_ ( err )
   if ( err .ne. 0 ) then
      call clean_()
      rc = 3
   end if

!  For now must specify chemReg
!  ----------------------------
   if ( present(chemReg) ) then
      reg = chemReg
   else
      rc = 4 ! for now, later Reg will come from file
      return
   end if

!  Get file attributes
!  -------------------
   if ( gfio_prec_== 64 ) then
      call GFIO_Inquire ( fid, im, jm, km, lm, nvars,     &
                          title, source, contact, amiss8, &
                          lon8, lat8, lev8, levunits,     &
                          yyyymmdd, hhmmss, timinc,       &
                          vname, vtitle, vunits,          &
                          kmvar, valid_range8, packing_range8, err )
      amiss=amiss8;lon=lon8;lat=lat8;lev=lev8
   else
      call GFIO_Inquire ( fid, im, jm, km, lm, nvars,     &
                          title, source, contact, amiss,  &
                          lon, lat, lev, levunits,        &
                          yyyymmdd, hhmmss, timinc,       &
                          vname, vtitle, vunits,          &
                          kmvar, valid_range , packing_range, err )
   end if
   if ( err .ne. 0 ) then
      call clean_()
      rc = 5
   end if

   if ( present(freq) ) then
        freq = timinc
   end if

!  Loop over variables and detect whether we have the new
!  GEOS-5 convention where all variables are upper case.
!  ------------------------------------------------------
   all_upper = .TRUE.
   do n = 1, nvars
      if ( trim(vname(n)) /= uppercase(trim(vname(n))) ) &
           all_upper = .FALSE.
   end do

   if ( all_upper ) &
        print *, "Chem_BundleRead: Using GEOS-5 all upercase mode"

!  Pick time to return
!  -------------------
   if ( present(timidx) ) then
        if ( timidx .eq. 0 ) then
             continue  ! nothing to do, nymd/nhms set on input
        else if ( timidx .lt. 0 .or. timidx .gt. lm ) then
           call clean_()
           rc = 6
           return
        else
           nymd = yyyymmdd(timidx)
           nhms = hhmmss(timidx)
        end if
   else 
      nymd = yyyymmdd(lm)   
      nhms = hhmmss(lm)
   end if

   if ( nvars < 2 ) then
      rc = 7
      call clean_()
      return
   end if


!  Allocate memory if necessary 
!  ----------------------------
   i1 = 1;  i2 = im; ig = 0
   j1 = 1;  j2 = jm; jg = 0 
   allocate( lat2d(i1:i2,j1:j2), lon2d(i1:i2,j1:j2), stat=err)
   if(err .ne. 0) then
    call clean_()
    rc = 11
    return
   endif
   do i = i1, i2
    lat2d(i,:) = lat
   enddo
   do j = j1, j2
    lon2d(:,j) = lon
   enddo
   call Chem_BundleCreate_ ( reg, i1, i2, ig, im, j1, j2, jg, jm, km, &
                             w_c, err,    &
                             lat=lat2d, lon=lon2d, lev=lev,  &
                             levUnits = levUnits, do_Conc = do_Conc)
   deallocate(lat2d, lon2d, stat=err)
   if(err .ne. 0) then
    call clean_()
    rc = 12
    return
   endif
   if ( err .eq. 1 ) then                     ! already allocated
        if ( w_c%grid%im .ne. im  .OR.  &
             w_c%grid%jm .ne. jm  .OR.  &
             w_c%grid%km .ne. km  .OR.  &
             w_c%reg%nq .ne. reg%nq  ) then
             rc = 8                           ! current size not compatible
             call clean_()
             return
        end if
   else if ( err .ne. 0 ) then
        rc = 9
        call clean_()
        return
   end if

!  Verify that variables in register are on file
!  ---------------------------------------------
   if ( reg%nq > nvars-1 ) then
        rc = 200
        return
   end if
   do j = 1, reg%nq
      ivar(j) = -1
      if ( all_upper ) then
           vname_ = uppercase(trim(reg%vname(j)))
      else
           vname_ = trim(reg%vname(j))
      end if
      do i = 1, nvars
         if ( trim(vname_) .eq. trim(vname(i)) ) then
              ivar(j) = i
         end if
      end do
      if ( ivar(j) < 1 ) then
         print *, 'Missing variable: ', trim(reg%vname(j))
         rc = 10
         call clean_()
         return
      end if
   end do
      

!  retrieve the variables
!  ----------------------
   if ( gfio_prec_== 64 ) then
      if ( all_upper ) then 
         call GFIO_GetVar ( fid, 'DELP', nymd, nhms,     &
                            im, jm, 1, km, rank3, err )
      else
         call GFIO_GetVar ( fid, 'delp', nymd, nhms,     &
                            im, jm, 1, km, rank3, err )
      endif
      if ( err .ne. 0 ) then
         rc = 101
      else
         w_c%delp=rank3
      end if

      call GFIO_GetVar ( fid, 'RH', nymd, nhms,     &
                         im, jm, 1, km, rank3, err )
      if ( err .ne. 0 ) then
         call GFIO_GetVar ( fid, 'rh', nymd, nhms,     &
                            im, jm, 1, km, rank3, err )
      end if
      if ( err .ne. 0 ) then
           w_c%rh = w_c%missing_value
        w_c%has_rh = .false.
      else
           w_c%has_rh = .true.  ! for backward compatibility
           w_c%rh = rank3
      end if

      if (w_c%do_concentration) then 
         if (all_upper) then 
            call GFIO_GetVar ( fid, 'AIRDENS', nymd, nhms,     &
                            im, jm, 1, km, rank3, err )
         else
            call GFIO_GetVar ( fid, 'airdens', nymd, nhms,     &
                            im, jm, 1, km, rank3, err )
         endif
      if ( err .ne. 0 ) then
         rc = 101
      else
         w_c%airdens=rank3
      end if
      endif         
   else
      if ( all_upper ) then 
         call GFIO_GetVar ( fid, 'DELP', nymd, nhms,     &
                            im, jm, 1, km, w_c%delp, err )
      else
         call GFIO_GetVar ( fid, 'delp', nymd, nhms,     &
                            im, jm, 1, km, w_c%delp, err )
      endif

      if ( err .ne. 0 )                                        rc = 101
      call GFIO_GetVar ( fid, 'RH', nymd, nhms,     &
                         im, jm, 1, km, w_c%rh, err )
      if ( err .ne. 0 ) then
         call GFIO_GetVar ( fid, 'rh', nymd, nhms,     &
                            im, jm, 1, km, w_c%rh, err )
      end if
      if ( err .ne. 0 ) then
           w_c%rh = w_c%missing_value
        w_c%has_rh = .false.
      else
           w_c%has_rh = .true.  ! for backward compatibility
      end if
      
      if (w_c%do_concentration) then
         if ( all_upper ) then 
            call GFIO_GetVar ( fid, 'AIRDENS', nymd, nhms,     &
                            im, jm, 1, km, w_c%airdens, err )
         else
            call GFIO_GetVar ( fid, 'airdens', nymd, nhms,     &
                            im, jm, 1, km, w_c%airdens, err )
         endif
      endif
   end if ! <prec_>
   do n = 1, reg%nq
        l = ivar(n)
        if ( all_upper ) then
             vname_ = uppercase(trim(vname(l)))
        else
             vname_ = trim(vname(l))
        end if
        if ( gfio_prec_== 64 ) then

           if ( w_c%do_concentration ) then

              call GFIO_GetVar ( fid, vname_, nymd, nhms,         & 
                              im, jm, 1, km, rank3, err )
               w_c%qa(n)%data3d(:,:,:)=rank3*w_c%airdens

           else
              call GFIO_GetVar ( fid, vname_, nymd, nhms,         & 
                              im, jm, 1, km, rank3, err )
              w_c%qa(n)%data3d(:,:,:)=rank3
           endif

        else

           if ( w_c%do_concentration ) then

              call GFIO_GetVar ( fid, vname_, nymd, nhms,         & 
                              im, jm, 1, km, w_c%qa(n)%data3d(:,:,:)*w_c%airdens, err )
           else
              call GFIO_GetVar ( fid, vname_, nymd, nhms,         & 
                              im, jm, 1, km, w_c%qa(n)%data3d(:,:,:), err )
           endif                   
        endif
        if ( err .ne. 0 )                                  rc = 100 + l
   end do


!  Retrieve vertical grid attributes
!  ---------------------------------
    if ( gfio_prec_==64 ) then
       call GFIO_GetRealAtt ( fid, 'ptop',   1, buf8, err )
       w_c%grid%ptop = buf8(1)
    else
       call GFIO_GetRealAtt ( fid, 'ptop',   1, buf , err )
       w_c%grid%ptop = buf(1)
    endif
    if ( err .ne. 0 ) then
       w_c%grid%ptop = 1. ! do not fuss about this
       rc = 0
    end if

!   Close GFIO file
!   ---------------
    call GFIO_close ( fid, err )

!  All done
!  --------
   call clean_()
   return

  CONTAINS

     subroutine init_ ( err )       ! allocates local memory
     integer err, err1, err2, err3, err4
     allocate( lat(jm), lon(im), lev(km), yyyymmdd(lm), stat=err1 )
     allocate( hhmmss(lm),vname(nvars), vunits(nvars), stat=err2 )
     allocate( vtitle(nvars),kmvar(nvars),valid_range(2,nvars),stat=err3 )
     allocate( packing_range(2,nvars),ivar(nvars),stat=err4) 
     if (gfio_prec_==64) then
        allocate(lon8(im))
        allocate(lat8(jm))
        allocate(lev8(km))
        allocate(packing_range8(2,nvars))
        allocate(valid_range8(2,nvars))
        allocate(rank3(im,jm,km))
     end if
     err = err1 + err2 + err3 + err4 
     end subroutine init_

     subroutine clean_()             ! de-allocates local memory
     deallocate (lat,lon,lev, yyyymmdd, hhmmss, vname, vunits, stat=err)
     deallocate (vtitle,kmvar,valid_range,packing_range, ivar, stat=err)
     if (gfio_prec_==64) then
        deallocate(rank3)
        deallocate(lon8)
        deallocate(lat8)
        deallocate(lev8)
        deallocate(packing_range8)
        deallocate(valid_range8)
     endif
     end subroutine clean_

   end subroutine Chem_BundleRead

!-------------------------------------------------------------------------
!      NASA/GSFC Global Modeling & Assimilation Office, Code 900.3       !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Chem_BundleSetPtr - Make sure internal pointers are set
! 
! !INTERFACE:
!
  subroutine  Chem_BundleSetPtr ( w_c, rc )
!
! !USES:
!
  implicit NONE
!
! !INPUT/OUTPUT PARAMETERS: 
!
  type(Chem_Bundle), intent (inout) :: w_c   ! chemical bundle

! !OUTPUT PARAMETERS:

  integer, intent(out)          ::  rc     ! Error return code:
                                           !  0 - all is well
                                           !  1 - 
! !DESCRIPTION: 
!
!  Make sure the internal array of pointers points to the current
!  memory location of the 4D q array
!
! !REVISION HISTORY: 
!
!  22Jul2005 da Silva  Initial code.
!
!EOP
!-------------------------------------------------------------------------

   _UNUSED_DUMMY(w_c)
   rc = 0

   return

   end subroutine Chem_BundleSetPtr

!-------------------------------------------------------------------------
!      NASA/GSFC Global Modeling & Assimilation Office, Code 900.3       !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Chem_BundleStat --- Prints vital stats of Chem_Bundle
!
! !INTERFACE:
!
  subroutine Chem_BundleStat ( lu, w_c, rc )
!
! !USES:
!
  implicit NONE
!
! !INPUT PARAMETERS:
!
  integer, intent(in)               :: lu     ! FORTRAN unit number for ASCII output
  type(Chem_Bundle), intent(in)     :: w_c    ! dynamics state vector

!
! !OUTPUT PARAMETERS:
!

  integer, intent(out), optional    :: rc    ! error return code:
                                             !  0 - all is well
                                             !  >0 - errors
!
! !DESCRIPTION: This routine prints basic stats about a dynamics state vector.
!
! !REVISION HISTORY:
!
!  15sep2003 da Silva  Adapted from dyn_stat
!
!EOP
!-------------------------------------------------------------------------

       integer :: i1, i2, im, j1, j2, jm, km, nq, l,  ios, k
       real, allocatable :: levs(:)
       real :: amiss

       if ( present(rc) ) rc = 0
       i1 = w_c%grid%i1; j1 = w_c%grid%j1
       i2 = w_c%grid%i2; j2 = w_c%grid%j2
       im = w_c%grid%im; jm = w_c%grid%jm
       km = w_c%grid%km; nq = w_c%reg%nq
       amiss = w_c%missing_value
       allocate ( levs(km), stat = ios )
       if ( ios .ne. 0 ) then
          if ( present(rc) ) rc = 1
          return
       end if

       do k = 1, km
          levs(k) = k
       end do

!!!   do k = 1, km
!!!      print *, 'k, delp = ', k, minval(w_c%delp(i1:i2,j1:j2,k)), &
!!!                             maxval(w_c%delp(i1:i2,j1:j2,k))
!!!   end do
!!!   do k = 1, km
!!!      print *, 'k,   rh = ', k, minval(w_c%rh(i1:i2,j1:j2,k)), &
!!!                                maxval(w_c%rh(i1:i2,j1:j2,k))
!!!   end do

       write(lu,*) 'Chem_BundleStat: i1, i2, im ', i1, i2, im
       write(lu,*) 'Chem_BundleStat: j1, j2, jm ', j1, j2, jm
       write(lu,*) 'Chem_BundleStat: km, ptop: ', km, w_c%grid%ptop

       call GDSTAT_ (lu,im,jm,km,w_c%delp, &
                     levs,'PRES','lev ',amiss, &
                     'delp', 1 )

       call GDSTAT_ (lu,im,jm,km,w_c%rh, &
                     levs,'PRES','lev ',amiss, &
                     'rh', 1 )

       do l = 1, nq
         call GDSTAT_ (lu,im,jm,km,w_c%qa(l)%data3d(i1:i2,j1:j2,1:km), &
                       levs,'PRES','lev',amiss, &
                       w_c%reg%vtitle(l), 1 )
       end do

       deallocate ( levs )

CONTAINS

        subroutine GDSTAT_ (lu,mx,my,mz,a,h,atype,htype,amiss,header,inc)

!       Print statistics of one 3-d variable. This is from the PSAS library.
!       It is reproduced here to avoid unnecessary dependencies.

        implicit none

        integer lu              ! Output unit
        integer mx,my,mz        ! Array sizes
        real a(mx,my,mz)        ! The array
        real h(mz)              ! The argument(levels)
        character(*), intent(in) :: atype       ! Type of the variable
        character(*), intent(in) :: htype       ! Typf of the levels
        real amiss              ! missing value flag of a
        character*(*) header    ! A header message
        integer inc             ! order of the listing

        integer i,j,k
        integer kfr,kto,kinc
        integer imx,imn,jmx,jmn
        integer knt
        real amx,amn
        real avg,dev,d
        logical first

!       ..A practical value for the magnitude of the fraction of a real
!       number.

        real rfrcval
        parameter(rfrcval=1.e-5)

        character(len=nch) dash

!       ..function


        logical spv
        real aspv
        spv(aspv)=abs((aspv-amiss)/amiss).le.rfrcval

        do i = 1, nch
           dash(i:i) = '-'
        end do
        i = len(trim(header))
        write(lu,'(//a)') trim(header)
        write(lu,'(a/)')  dash(1:i)
        if(htype.eq.'PRES') then
          write(lu,'(a,3x,a,2x,a,5x,a,6x,a,9x,a,15x,a)') 'lvl','mbar', &
            'count','mean','stdv','maxi','mini'
        elseif(htype.eq.'HGHT') then
          write(lu,'(a,2x,a,2x,a,5x,a,6x,a,9x,a,15x,a)') 'lvl','meter', &
            'count','mean','stdv','maxi','mini'
        elseif(htype.eq.'TEMP') then
          write(lu,'(a,4x,a,4x,a,5x,a,6x,a,9x,a,15x,a)') 'lvl','K', &
            'count','mean','stdv','maxi','mini'
        else
          write(lu,'(a,4x,a,4x,a,5x,a,6x,a,9x,a,15x,a)') 'lvl',htype, &
            'count','mean','stdv','maxi','mini'
        endif

!       ..Check the order of the listing, increase or decrease
        if(inc.ge.0) then
          kfr=1
          kto=mz
          kinc=1
        else
          kfr=mz
          kto=1
          kinc=-1
        endif

        do k=kfr,kto,kinc
          knt=0
          avg=0.
          do j=1,my
            do i=1,mx
              if(.not.spv(a(i,j,k))) then
                knt=knt+1
                avg=avg+a(i,j,k)
              endif
            end do
          end do
          avg=avg/max(1,knt)

          dev=0.
          do j=1,my
            do i=1,mx
              if(.not.spv(a(i,j,k))) then
                d=a(i,j,k)-avg
                dev=dev+d*d
              endif
            end do
          end do
          dev=sqrt(dev/max(1,knt-1))

          amx=a(1,1,k)
          amn=a(1,1,k)
          first=.true.
          do j=1,my
            do i=1,mx
              if(.not.spv(a(i,j,k))) then
                if(first) then
                  imx=i
                  imn=i
                  jmx=j
                  jmn=j
                  amx=a(imx,jmx,k)
                  amn=a(imn,jmn,k)
                  first=.false.
                else
                  if(a(i,j,k).gt.amx) then
                    amx=a(i,j,k)
                    imx=i
                    jmx=j
                  endif
                  if(a(i,j,k).lt.amn) then
                    amn=a(i,j,k)
                    imn=i
                    jmn=j
                  endif
                endif
              endif
            end do
          end do

          if(atype.eq.'RELH') then
            avg=avg*100.
            dev=dev*100.
            amx=amx*100.
            amn=amn*100.
          endif

          if(htype.eq.'PRES'.or.htype.eq.'HGHT') then
            if(atype.eq.'HGHT'.or.atype.eq.'STRM') then
              write(lu,'(i3,2i7,2i10,'// &
                '2(i10,a,i3,a,i3,a))') &
                k,nint(h(k)),knt,nint(avg),nint(dev), &
                nint(amx),'(',imx,',',jmx,')', &
                nint(amn),'(',imn,',',jmn,')'
            elseif(atype.eq.'TEMP'.or.atype.eq.'PRES'.or. &
              atype.eq.'WIND'.or.atype.eq.'%REH'.or. &
              atype.eq.'RELH'.or.atype.eq.'MIXR') then
              write(lu,'(i3,2i7,2f10.2,'// &
                '2(f10.2,a,i3,a,i3,a))') &
                k,nint(h(k)),knt,avg,dev, &
                amx,'(',imx,',',jmx,')', &
                amn,'(',imn,',',jmn,')'
            elseif(atype.eq.'NORM') then
              write(lu,'(i3,2i7,2f10.4,'// &
                '2(f10.4,a,i3,a,i3,a))') &
                k,nint(h(k)),knt,avg,dev, &
                amx,'(',imx,',',jmx,')', &
                amn,'(',imn,',',jmn,')'
            else
              write(lu,'(i3,2i7,1p,2e10.3e1,0p,'// &
                '2(1p,e10.3e1,0p,a,i3,a,i3,a))') &
                k,nint(h(k)),knt,avg,dev, &
                amx,'(',imx,',',jmx,')', &
                amn,'(',imn,',',jmn,')'
            endif

          elseif(htype.eq.'TEMP') then
            if(atype.eq.'HGHT'.or.atype.eq.'STRM') then
              write(lu,'(i3,f7.2,i7,2i10,'// &
                '2(i10,a,i3,a,i3,a))') &
                k,h(k),knt,nint(avg),nint(dev), &
                nint(amx),'(',imx,',',jmx,')', &
                nint(amn),'(',imn,',',jmn,')'
            elseif(atype.eq.'TEMP'.or.atype.eq.'PRES'.or. &
              atype.eq.'WIND'.or.atype.eq.'%REH'.or. &
              atype.eq.'RELH'.or.atype.eq.'MIXR') then
              write(lu,'(i3,f7.2,i7,2f10.2,'// &
                '2(f10.2,a,i3,a,i3,a))') &
                k,h(k),knt,avg,dev, &
                amx,'(',imx,',',jmx,')', &
                amn,'(',imn,',',jmn,')'
            elseif(atype.eq.'NORM') then
              write(lu,'(i3,f7.2,i7,2f10.4,'// &
                '2(f10.4,a,i3,a,i3,a))') &
                k,h(k),knt,avg,dev,&
                amx,'(',imx,',',jmx,')', &
                amn,'(',imn,',',jmn,')'
            else
              write(lu,'(i3,f7.2,i7,1p,2e10.3e1,0p,'// &
                '2(1p,e10.3e1,0p,a,i3,a,i3,a))') &
                k,h(k),knt,avg,dev, &
                amx,'(',imx,',',jmx,')', &
                amn,'(',imn,',',jmn,')'
            endif

          else
            if(atype.eq.'HGHT'.or.atype.eq.'STRM') then
              write(lu,'(i3,1p,e10.3e1,0p,i7,2i10,'// &
                '2(i10,a,i3,a,i3,a))') &
                k,h(k),knt,nint(avg),nint(dev), &
                nint(amx),'(',imx,',',jmx,')', &
                nint(amn),'(',imn,',',jmn,')'
            elseif(atype.eq.'TEMP'.or.atype.eq.'PRES'.or. &
              atype.eq.'WIND'.or.atype.eq.'%REH'.or. &
              atype.eq.'RELH'.or.atype.eq.'MIXR') then
              write(lu,'(i3,1p,e10.3e1,0p,i7,2f10.2,'// &
                '2(f10.2,a,i3,a,i3,a))') &
                k,h(k),knt,avg,dev, &
                amx,'(',imx,',',jmx,')', &
                amn,'(',imn,',',jmn,')'
            elseif(atype.eq.'NORM') then
              write(lu,'(i3,1p,e10.3e1,0p,i7,2f10.4,'// &
                '2(f10.4,a,i3,a,i3,a))') &
                k,h(k),knt,avg,dev, &
                amx,'(',imx,',',jmx,')', &
                amn,'(',imn,',',jmn,')'
            else
              write(lu,'(i3,1p,e10.3e1,0p,i7,1p,2e10.3e1,0p,'// &
                '2(1p,e10.3e1,0p,a,i3,a,i3,a))') &
                k,h(k),knt,avg,dev, &
                amx,'(',imx,',',jmx,')', &
                amn,'(',imn,',',jmn,')'
            endif

          endif
        end do          ! k=kfr,kto,kinc

      end  subroutine GDSTAT_

    end subroutine Chem_BundleStat


#ifdef FUTURE_SPMD

!
! The routines below are incomplete.
!

!-------------------------------------------------------------------------
!      NASA/GSFC Global Modeling & Assimilation Office, Code 900.3       !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Chem_BundleGather --- Gather bundle on master PE
!
! !INTERFACE:
!
  subroutine Chem_BundleGather ( w_c, w_g, rc )
!
! !USES:
!
!  use mod_comm
  implicit NONE
!
! !INPUT PARAMETERS:
!

  type(Chem_Bundle), intent(inout)     :: w_c    ! distributed bundle

!
! !OUTPUT PARAMETERS:
!

  type(Chem_Bundle), intent(inout)     :: w_g    ! global bundle

  integer, intent(out), optional    :: rc    ! error return code:
                                             !  0 - all is well
                                             !  >0 - errors
!
! !DESCRIPTION: This routine gather a chem bundle on single PE.
!
! !REVISION HISTORY:
!
!  01mar2005  da Silva  Initial version.
!
!EOP
!-------------------------------------------------------------------------

   type(Chem_Registry) :: reg
   integer :: i1, i2, ig
   integer :: j1, j2, jg
   integer :: km, nq

!  Short hand for dimensions
!  -------------------------
   i1 = w_c%grid%i1; j1 = w_c%grid%j1
   i2 = w_c%grid%i2; j2 = w_c%grid%j2
   ig = w_c%grid%i2; jg = w_c%grid%j2
   im = w_c%grid%im; jm = w_c%grid%jm
   km = w_c%grid%km; nq = w_c%reg%nq

   reg = w_c%reg

   rc = 0

!  Create global bundle on master
!  -----------------------------
   if ( gid .eq. 0 ) then
        call Chem_BundleCreate1PE_ ( reg, im, jm, km, w_c, rc )
        if ( rc /= 0 ) return
   end if

!   Gather from all subdomains to GID=0
!  ------------------------------------
   call mp_gather4d(delp, delptmp, im, jm, km, 1, jfirst, jlast, 1, km, 0, 0, 0)
end subroutine Chem_BundleGather

!-------------------------------------------------------------------------
!      NASA/GSFC Global Modeling & Assimilation Office, Code 900.3       !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Chem_BundleScatter --- Distributes a chem bundle
!
! !INTERFACE:
!
  subroutine Chem_BundleScatter ( w_g, w_c, rc )
!
! !USES:
!
  use mod_comm
  implicit NONE
!
! !INPUT PARAMETERS:
!

  type(Chem_Bundle), intent(inout)     :: w_c    ! distributed bundle

!
! !OUTPUT PARAMETERS:
!

  type(Chem_Bundle), intent(inout)     :: w_g    ! global bundle

  integer, intent(out), optional    :: rc    ! error return code:
                                             !  0 - all is well
                                             !  >0 - errors
!
! !DESCRIPTION: This routine distributes a chem bundle.
!
! !REVISION HISTORY:
!
!  01mar2005  da Silva  Initial version.
!
!EOP
!-------------------------------------------------------------------------

end subroutine Chem_BundleScatter

#endif

 end MODULE Chem_BundleMod
