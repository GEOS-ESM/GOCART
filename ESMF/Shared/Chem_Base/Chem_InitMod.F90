!-------------------------------------------------------------------------
!      NASA/GSFC, Global Modeling & Assimilation Office, Code 900.3      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  Chem_InitMod --- Chemistry Initialization Class
!
! !INTERFACE:
!

   module  chem_initmod

! !USES:

   use m_inpak90            ! resource file management (MPEU)
   use m_die, only: die
   implicit none

! !PUBLIC TYPES:
!
   PRIVATE
   PUBLIC  Chem_Init     ! Keeps track of desired initial conditions for
                         ! tracers to be specified in the Chem_Registry
                           
!
! !PUBLIIC MEMBER FUNCTIONS:
!
   PUBLIC  Chem_InitResource  ! Query the resource file

!
! !DESCRIPTION:
!
!  This module implements an initialization registry for (chemical)
!  constituents.
!  This initial class is intended to serve as a stop gap before an ESMF 
!  implementation is adopted. 
!
!
! !REVISION HISTORY:
!
!  21Oct2003 Colarco
!
!EOP
!-------------------------------------------------------------------------

  integer, parameter :: nch = 255

! Registry
! --------
  type Chem_Init

!    Amplitude
!    ---------
     real :: amp, lat0, lon0, z0, rx, ry, rz

  end type Chem_Init

CONTAINS

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Chem_InitCreate --- Cronstruct Chemisty Initialization
!
! !INTERFACE:
!

  Function Chem_InitResource ( name, rc, rcfile )

  implicit none
  type(Chem_Init) Chem_InitResource

! !USES:


! !INPUT PARAMETERS:

   character(len=*) :: name
   character(len=*), OPTIONAL :: rcfile  ! Resource file name; default is
                                         ! 'Chem_IC.rc'

! !OUTPUT PARAMETERS:

   integer, intent(out) ::  rc            ! Error return code:
                                          !  0 - all is well
                                          !  1 - 

! !DESCRIPTION:
!
!
! !REVISION HISTORY:
!
!  21Oct2003 Colarco
!
!EOP
!-------------------------------------------------------------------------

   character(len=*), parameter ::  myname = 'Chem_InitResource'

   type(Chem_Init) :: this
   character(len=255) :: rcfilen = 'Chem_Init.rc'
   integer :: ier

   rc = 0
                
!  ------------------------------------------------------
!  Parse resource file to see desired tracer parameters
!  ------------------------------------------------------
   if ( present(rcfile) ) rcfilen = trim(rcfile)
   call i90_loadf ( rcfilen, ier )
   if ( ier .ne. 0 ) call die(myname, 'could not read rc file '// &
        trim(rcfile) )
   call parserc_ ( name, this%amp, this%lat0, this%lon0, this%z0, &
                         this%ry, this%rx, this%rz)
   call I90_Release()

!  All done
!  --------
   Chem_InitResource = this
   
   return 

!                 -----------------------------
!                 Internal Constructor Routines
!                 -----------------------------

   CONTAINS

!     parses rc file
      subroutine parserc_ ( name, amp, lat0, lon0, z0, ry, rx, rz ) 
!     -------------------
      character(len=*), intent(in) :: name
      real, intent(out) :: amp, lat0, lon0, z0, ry, rx, rz

      character(len=255) :: answer
      integer ier

!     Defaults
!     --------
      lat0 = 0.
      lon0 = 90.
      z0 = 20.
      ry = 500.
      rx = 5000.
      rz = 5.

!     What amplitude?
!     -------------------------
      call i90_label ( 'amp_'//trim(name)//':', ier )
      if ( ier .eq. 0 ) then
         amp = i90_gint ( ier )
      end if
      if ( ier .ne. 0 ) then
         call die ( myname, 'cannot determine AMP for '//trim(name) )
      end if

!     What latitude?
!     -------------------------
      call i90_label ( 'lat0_'//trim(name)//':', ier )
      if ( ier .eq. 0 ) then
         lat0 = i90_gint ( ier )
      end if
      if ( ier .ne. 0 ) then
         call die ( myname, 'cannot determine LAT for '//trim(name) )
      end if

!     What longitude?
!     -------------------------
      call i90_label ( 'lon0_'//trim(name)//':', ier )
      if ( ier .eq. 0 ) then
         lon0 = i90_gint ( ier )
      end if
      if ( ier .ne. 0 ) then
         call die ( myname, 'cannot determine LON for '//trim(name) )
      end if

!     What z?
!     -------------------------
      call i90_label ( 'z0_'//trim(name)//':', ier )
      if ( ier .eq. 0 ) then
         z0 = i90_gint ( ier )
      end if
      if ( ier .ne. 0 ) then
         call die ( myname, 'cannot determine Z for '//trim(name) )
      end if

!     What rlat?
!     -------------------------
      call i90_label ( 'ry_'//trim(name)//':', ier )
      if ( ier .eq. 0 ) then
         ry = i90_gint ( ier )
      end if
      if ( ier .ne. 0 ) then
         call die ( myname, 'cannot determine RLAT for '//trim(name) )
      end if

!     What rlon?
!     -------------------------
      call i90_label ( 'rx_'//trim(name)//':', ier )
      if ( ier .eq. 0 ) then
         rx = i90_gint ( ier )
      end if
      if ( ier .ne. 0 ) then
         call die ( myname, 'cannot determine RLON for '//trim(name) )
      end if

!     What rz?
!     -------------------------
      call i90_label ( 'rz_'//trim(name)//':', ier )
      if ( ier .eq. 0 ) then
         rz = i90_gint ( ier )
      end if
      if ( ier .ne. 0 ) then
         call die ( myname, 'cannot determine RZ for '//trim(name) )
      end if

      end subroutine parserc_ 
      

 end Function Chem_InitResource


 end module Chem_InitMod

