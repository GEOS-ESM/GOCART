!BOP
!
! !MODULE:  GOCART2G_MieMod --- Reader for aerosol mie tables
!
! !INTERFACE:
!

module GOCART2G_MieMod
! !USES:
   use netcdf
   implicit none

! !PUBLIC TYPES:
!
   private
   public  GOCART2G_Mie             ! Holds Chem_MieTable and other information

!
! !PUBLIC MEMBER FUNCTIONS:
!
!
! !DESCRIPTION:
!
!  This module read the mie aerosol tables.
!
! !REVISION HISTORY:
!
!  23Mar2005 Colarco - Initial code.
!  31Mar2005 Todling - Declared netcdf nf_ routines as external (OSF1) 
!                      Removed # from include netcdf.inc
!  16Feb2022 da Silva  Refactoring from old Chem_MieTable and Chem_MieMod;
!                      these functionalities have been merged.   
!   
!EOP
!-------------------------------------------------------------------------
!
! Mie LUT table
! Will be reduced from input files to the desired channels
! --------
   integer, parameter :: NRH_BINS = 991

   type :: RH_Mie
      real, allocatable :: rh(:)
   end type RH_Mie

   type GOCART2G_Mie
      
      private
      character(len=:), allocatable :: table_name
      integer :: nch             ! number of channels in table (replacement of nlamfda)
      integer :: nrh             ! number of RH values in table
      integer :: nbin            ! number of size bins in table
      integer :: nMom            ! number of moments of phase function
      integer :: nPol            ! number of elements of scattering phase matrix

                                            ! c=channel, r=rh, b=bin, m=moments, p=nPol
      real, allocatable  :: wavelengths(:)  ! (c) wavelengths [m]
      real, allocatable  :: rh(:)           ! (r) RH values   [fraction]
      type(RH_Mie), allocatable  :: reff(:)       ! (b) effective radius [m]
      type(RH_Mie), allocatable  :: bext(:,:)     ! (c,b) bext values [m2 kg-1]
      type(RH_Mie), allocatable  :: bsca(:,:)     ! (c,b) bsca values [m2 kg-1]
      type(RH_Mie), allocatable  :: bbck(:,:)     ! (c,b) bbck values [m2 kg-1]
      type(RH_Mie), allocatable  :: g(:,:)        ! (c,b) asymmetry parameter
      type(RH_Mie), allocatable  :: pback(:,:,:)  ! (c,b,p) Backscatter phase function
      type(RH_Mie), allocatable  :: pmom(:,:,:,:) ! (c,b,m,p) moments of phase function
      type(RH_Mie), allocatable  :: gf(:)         ! (b) hygroscopic growth factor
      type(RH_Mie), allocatable  :: rhop(:)       ! (b) wet particle density [kg m-3]
      type(RH_Mie), allocatable  :: rhod(:)       ! (b) wet particle density [kg m-3]
      type(RH_Mie), allocatable  :: vol(:)        ! (b) wet particle volume [m3 kg-1]
      type(RH_Mie), allocatable  :: area(:)       ! (b) wet particle cross section [m2 kg-1]
      type(RH_Mie), allocatable  :: refr(:,:)     ! (c,b) real part of refractive index
      type(RH_Mie), allocatable  :: refi(:,:)     ! (c,b) imaginary part of refractive index

      integer            :: rhi(NRH_BINS)   ! pointer to rh LUT
      real               :: rha(NRH_BINS)   ! slope on rh LUT
      
   CONTAINS
      
      procedure :: QueryByWavelength_1d
      procedure :: QueryByWavelength_2d
      procedure :: QueryByWavelength_3d
      procedure :: QueryByChannel_1d
      procedure :: QueryByChannel_2d
      procedure :: QueryByChannel_3d
      generic   :: Query => QueryByWavelength_1d, &
                            QueryByWavelength_2d, &
                            QueryByWavelength_3d, &
                            QueryByChannel_1d,    &
                            QueryByChannel_2d,    &
                            QueryByChannel_3d
      procedure :: getChannel
      procedure :: getWavelength
   end type GOCART2G_Mie

   interface GOCART2G_Mie
      module procedure GOCART2G_MieCreate
   end interface GOCART2G_Mie

CONTAINS

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  GOCART2G_MieCreate  --- Construct Chemistry Registry
!
!
! !DESCRIPTION:
!
!
! !REVISION HISTORY:
!
!  09Mar2005 da Silva  API, prologues.
!

  type(GOCART2G_Mie) function GOCART2G_MieCreate ( rcfile, wavelengths, nmom, rc ) result (this)

! !INPUT PARAMETERS:

     character(len=*), intent(in) :: rcfile  ! Mie table file name
     real, intent(in) :: wavelengths(:)
     integer, optional, intent(in) :: nmom

! !OUTPUT PARAMETERS:

     integer, intent(out) ::  rc          ! Error return code:
                                          !  0 - all is well
! !DESCRIPTION:
!
!   Fills in the Mie table
!
! !REVISION HISTORY:
!
!  23Mar2005 Colarco
!
!EOP
!-------------------------------------------------------------------------
     character(len=*), parameter ::  myname = 'GOCART2G_MieCreate'
     integer :: nch
     integer :: ncid, idimid, ivarid, n, i, j, ip1
     integer :: nch_table, nrh_table, nbin_table, nmom_table, nPol_table
!    Tables are hard-wired as single precision
     real, allocatable ::   channels_table(:),    rh_table(:), reff_table(:,:),    &
                            bext_table(:,:,:),    bsca_table(:,:,:),               &
                            bbck_table(:,:,:),    g_table(:,:,:),                  &
                            pmom_table(:,:,:,:,:),pback_table(:,:,:,:),            &
                            gf_table(:,:),        rhop_table(:,:), rhod_table(:,:),&
                            vol_table(:,:),       area_table(:,:),                 &
                            refr_table(:,:,:),    refi_table(:,:,:)

     real :: yerr
     integer :: nmom_, imom, ipol
     integer :: status

#define NF_VERIFY_(expr) rc = expr; if (rc /= 0) return
#define __NF_STAT__ stat=status); NF_VERIFY_(status

     rc = 0
     this%table_name = rcfile
     nch = size(wavelengths)

!    Whether or not doing phase function
!    -----------------------------------
     if ( present(nmom) ) then
        nmom_ = nmom
     else
        nmom_ = 0
     end if

!    Set up nPol_table for reading backscatter phase function
!    This will get overwritten if pmoments is requested
!     --------------------------------------------------------
     nPol_table = 6

!     Open the table and get the dimensions
!     -------------------------------------
      rc = nf90_open(this%table_name, NF90_NOWRITE, ncid)
      IF ( rc /= 0 ) THEN
        print *, 'nf90_open '//this%table_name//'  RETURN CODE=', rc
        NF_VERIFY_(rc)
      END IF

!     RH
!     --
      NF_VERIFY_(nf90_inq_dimid(ncid,'rh',idimid))
      NF_VERIFY_(nf90_inquire_dimension(ncid,idimid,len=nrh_table))

!     Channels
!     --------
      NF_VERIFY_(nf90_inq_dimid(ncid,'lambda',idimid))
      NF_VERIFY_(nf90_inquire_dimension(ncid,idimid,len=nch_table))

!     Dry Effective radius
!     --------------------
      NF_VERIFY_(nf90_inq_dimid(ncid,'radius',idimid))
      NF_VERIFY_(nf90_inquire_dimension(ncid,idimid,len=nbin_table))

!     Moments of phase function
!     -------------------------
      if ( nmom_ > 0 ) then
         NF_VERIFY_(nf90_inq_dimid(ncid,'nMom',idimid))
         NF_VERIFY_(nf90_inquire_dimension(ncid,idimid,len=nmom_table))
         if ( nmom_ > nmom_table ) then
!            rc = 99
            print*,'Error: nmom_ > nmom_table, see:'//myname
            NF_VERIFY_(1)
         end if
         NF_VERIFY_(nf90_inq_dimid(ncid,'nPol',idimid))
         NF_VERIFY_(nf90_inquire_dimension(ncid,idimid,len=nPol_table))
      endif

!     Get the table contents
!     -------------------------------------

      allocate(channels_table(nch_table), __NF_STAT__)
      allocate(rh_table(nrh_table), __NF_STAT__)
      allocate(reff_table(nrh_table,nbin_table), __NF_STAT__)
      allocate(bext_table(nch_table,nrh_table,nbin_table), __NF_STAT__)
      allocate(bsca_table(nch_table,nrh_table,nbin_table), __NF_STAT__)
      allocate(bbck_table(nch_table,nrh_table,nbin_table),  __NF_STAT__)
      allocate(g_table(nch_table,nrh_table,nbin_table), stat = rc )
      allocate(pback_table(nch_table,nrh_table,nbin_table,nPol_table),  __NF_STAT__)
      allocate(gf_table(nrh_table,nbin_table), __NF_STAT__)
      allocate(rhop_table(nrh_table,nbin_table), __NF_STAT__)
      allocate(rhod_table(nrh_table,nbin_table), __NF_STAT__)
      allocate(vol_table(nrh_table,nbin_table), __NF_STAT__)
      allocate(area_table(nrh_table,nbin_table), __NF_STAT__)
      allocate(refr_table(nch_table,nrh_table,nbin_table), __NF_STAT__)
      allocate(refi_table(nch_table,nrh_table,nbin_table), __NF_STAT__)

      if ( nmom_ > 0 ) then
         allocate(pmom_table(nch_table,nrh_table,nbin_table,nmom_table,nPol_table), __NF_STAT__)
      end if
      NF_VERIFY_(nf90_inq_varid(ncid,'lambda',ivarid))
      NF_VERIFY_(nf90_get_var(ncid,ivarid,channels_table))
      NF_VERIFY_(nf90_inq_varid(ncid,'rEff',ivarid))
      NF_VERIFY_(nf90_get_var(ncid,ivarid,reff_table))
      NF_VERIFY_(nf90_inq_varid(ncid,'bext',ivarid))
      NF_VERIFY_(nf90_get_var(ncid,ivarid,bext_table))
      NF_VERIFY_(nf90_inq_varid(ncid,'bsca',ivarid))
      NF_VERIFY_(nf90_get_var(ncid,ivarid,bsca_table))
      NF_VERIFY_(nf90_inq_varid(ncid,'bbck',ivarid))
      NF_VERIFY_(nf90_get_var(ncid,ivarid,bbck_table))
      NF_VERIFY_(nf90_inq_varid(ncid,'g',ivarid))
      NF_VERIFY_(nf90_get_var(ncid,ivarid,g_table))
      NF_VERIFY_(nf90_inq_varid(ncid,'rh',ivarid))
      NF_VERIFY_(nf90_get_var(ncid,ivarid,rh_table))

      ! TODO: we need to look at these NF90_NOERR checks
!     Get the backscatter phase function values
      rc = nf90_inq_varid(ncid,'pback',ivarid)
      if(rc .ne. NF90_NOERR) then   ! pback not in table, fill in dummy variable
        pback_table = 1.
      else
        NF_VERIFY_(nf90_get_var(ncid,ivarid,pback_table))
      endif

      if ( nmom_ > 0 ) then
         NF_VERIFY_(nf90_inq_varid(ncid,'pmom',ivarid))
         NF_VERIFY_(nf90_get_var(ncid,ivarid,pmom_table))
      end if

!     Aerosol optical properties not necessarily stored in all versions of the tables
!     ----------------------
!     Particle growth factor
      rc = nf90_inq_varid(ncid,'growth_factor',ivarid)
      if(rc .ne. NF90_NOERR) then   ! not in table, fill in dummy variable
        gf_table = -999.
      else
        NF_VERIFY_(nf90_get_var(ncid,ivarid,gf_table))
      endif

!     Wet particle density
      rc = nf90_inq_varid(ncid,'rhop',ivarid)
      if(rc .ne. NF90_NOERR) then   ! not in table, fill in dummy variable
        rhop_table = -999.
      else
        NF_VERIFY_(nf90_get_var(ncid,ivarid,rhop_table))
      endif

!     Dry particle density (will be pulled from wet particle radius)
      rc = nf90_inq_varid(ncid,'rhop',ivarid)
      if(rc .ne. NF90_NOERR) then   ! not in table, fill in dummy variable
        rhod_table = -999.
      else
        rc = nf90_get_var(ncid,ivarid,rhod_table)
        do i = 1, nrh_table
          rhod_table(i,:) = rhod_table(1,:)
        enddo
        if (rc /=0) return
      endif

!     Wet particle real part of refractive index
      rc = nf90_inq_varid(ncid,'refreal',ivarid)
      if(rc .ne. NF90_NOERR) then   ! not in table, fill in dummy variable
        refr_table = -999.
      else
        NF_VERIFY_(nf90_get_var(ncid,ivarid,refr_table))
      endif

!     Wet particle imaginary part of refractive index (ensure positive)
      rc = nf90_inq_varid(ncid,'refimag',ivarid)
      if(rc .ne. NF90_NOERR) then   ! not in table, fill in dummy variable
        refi_table = -999.
      else
        NF_VERIFY_(nf90_get_var(ncid,ivarid,refi_table))
        refi_table = abs(refi_table)
      endif

!     Wet particle volume [m3 kg-1]
!     Ratio of wet to dry volume is gf^3, hence the following
      vol_table = gf_table**3 / rhod_table

!     Wet particle cross sectional area [m2 kg-1]
!     Assume area is volume divided by (4./3.*reff)
      area_table = vol_table / (4./3.*reff_table)

!     Close the table file
!     -------------------------------------
      NF_VERIFY_(nf90_close(ncid))


!     Setup the table to be returned
!     -------------------------------------
      this%nch  = nch
      this%nrh  = nrh_table
      this%nbin = nbin_table
      this%nMom = nmom_
      this%nPol = nPol_table

!     Preserve the full RH structure of the input table
      this%rh = rh_table ! assignment does allocation

!     Insert the requested channels in the output table
      this%wavelengths = wavelengths

      call fill_table2d(this%reff, reff_table)
      !Insert growth factor
      call fill_table2d(this%gf,   gf_table)
      !Wet particle density [kg m-3]
      call fill_table2d(this%rhop, rhop_table)
      !Dry particle density [kg m-3]
      call fill_table2d(this%rhod, rhod_table)
      !Volume [m3 kg-1]
      call fill_table2d(this%vol,  vol_table)
      !Area [m2 kg-1]
      call fill_table2d(this%area, area_table)


      call fill_table3d(this%bext, bext_table, channels_table, wavelengths, yerr)
      call fill_table3d(this%bsca, bsca_table, channels_table, wavelengths, yerr)
      call fill_table3d(this%bbck, bbck_table, channels_table, wavelengths, yerr)
      call fill_table3d(this%g,       g_table, channels_table, wavelengths, yerr)
      call fill_table3d(this%refi, refi_table, channels_table, wavelengths, yerr)
      call fill_table3d(this%refr, refr_table, channels_table, wavelengths, yerr)

      call fill_table4d(this%pback, pback_table, channels_table, wavelengths, yerr)

      if ( nmom_ > 0 ) then
         call fill_table5d(this%pmom, pmom_table, channels_table, wavelengths, yerr)
      endif

!     Now we do a mapping of the RH from the input table to some high
!     resolution representation.  This is to spare us the need to
!     do a full-up interpolation later on.
!     RH input from the table is scaled 0 - 0.99
!     We resolve the map to 0 - 0.990 in steps of 0.001 (991 total steps)
      do j = 1, NRH_BINS
         do i = this%nrh, 1, -1
            if ((j-1) .ge. int(this%rh(i)*1000)) then
               ip1 = i + 1
               this%rhi(j) = i
               if (ip1 .gt. this%nrh) then
                  this%rha(j) = 0.
               else
                  this%rha(j) =  ( (j-1)/1000. - this%rh(i)) &
                                /  ( this%rh(ip1)- this%rh(i))
               endif
               exit
            endif
         enddo
!       print *, j, this%rhi(j), this%rha(j), this%rh(this%rhi(j))
      enddo

      return

  contains

     subroutine fill_table2d(table, array)
       type(RH_Mie), allocatable, intent(out) :: table(:)
       real, intent(in) :: array(:,:)
       integer :: j  
       associate (nbin => this%nbin, nrh => this%nrh)
         allocate(table(nbin))
         do j = 1, nbin
            table(j)%rh = array(:,j)
         enddo
       end associate
     end subroutine

     subroutine fill_table3d(table, array, channels_table, wavelengths, yerr)
       type(RH_Mie), allocatable, intent(out) :: table(:,:)
       real, intent(in) :: array(:,:,:)
       real, intent(in) :: channels_table(:), wavelengths(:)
       real, intent(out) :: yerr
       integer :: j,n

       associate (nbin => this%nbin, nch => this%nch, nrh => this%nrh)
         allocate(table(nch, nbin), __NF_STAT__)
         do j = 1, nbin
            do n = 1, nch
               allocate (table(n,j)%rh(nrh), __NF_STAT__)
               call polint(channels_table, array(:,:,j),nch, &
                     this%wavelengths(n),table(n,j),yerr)
            enddo
         enddo
       end associate
     end subroutine fill_table3d

     subroutine fill_table4d(table, array, channels_table, wavelengths, yerr)
       type(RH_Mie), allocatable, intent(out) :: table(:,:,:)
       real, intent(in) :: array(:,:,:,:)
       real, intent(in) :: channels_table(:), wavelengths(:)
       real, intent(out) :: yerr
       integer :: j,n, ipol

       associate (nbin => this%nbin, nch => this%nch, nrh => this%nrh, npol=> this%npol)
         allocate(table(nch, nbin, npol), __NF_STAT__)
         do j = 1, nbin
            do n = 1, nch
               do ipol = 1, npol
                  allocate (table(n,j,ipol)%rh(nrh), __NF_STAT__)
                  call polint(channels_table, array(:,:,j,ipol),nch, &
                              wavelengths(n),table(n,j,ipol),yerr)
               enddo
            enddo
         enddo
       end associate
     end subroutine fill_table4d

     subroutine fill_table5d(table, array, channels_table, wavelengths, yerr)
       type(RH_Mie), allocatable, intent(out) :: table(:,:,:,:)
       real, intent(in) :: array(:,:,:,:,:)
       real, intent(in) :: channels_table(:), wavelengths(:)
       real, intent(out) :: yerr
       integer :: j, n, ipol, imom

       associate (nbin => this%nbin, nch => this%nch, nrh => this%nrh, npol=> this%npol, nmom => this%nmom)
         allocate(table(nch, nbin, nmom, npol), __NF_STAT__)
         do j = 1, nbin
            do n = 1, nch
               do imom = 1, nmom
                  do ipol = 1, npol
                     allocate(table(n,j,imom,ipol)%rh(nrh), __NF_STAT__)
                     call polint(channels_table, array(:,:,j,imom, ipol),nch, &
                                 wavelengths(n),table(n,j,imom,ipol),yerr)
                  enddo
               enddo
            enddo
         enddo
       end associate
     end subroutine fill_table5d

     subroutine polint(x,y,n,xWant,yWant,yErr)
       integer :: n
       !  recall, table hard-wired single precision
       real   :: x(:),y(:,:)
       real   :: xWant, yErr
       type(RH_Mie) :: yWant

       !  given array x(n) of independent variables and array y(n) of dependent
       !  variables, compute the linear interpolated result yWant at xWant and return
       !  with a dummy error estimate yErr.  Hacked up from Numerical Recipes Chapter 3

       integer :: i, j, k
       real    :: dx, slope
       character(len=255) :: msg

       do k = 1, size(y,2)

          !  on out of bounds, set i to lower or upper limit
          i = 0
          if(xWant .lt. x(1)) then
          !   write(msg,*) "in polint, wanted: ", xWant, ", got lower bound: ", x(1)
          !   call warn(myname,msg)
          !    if (mapl_am_i_root()) print *,'in polint(), wnted: ', xWant, ', got lower bound: ', x(1)
            i = 1
          endif
          if(xWant .gt. x(n)) then
          !    write(msg,*) "in polint, wanted: ", xWant, ", got upper bound: ", x(n)
          !    call warn(myname,msg)
          !    if (mapl_am_i_root()) print *,'in polint(), wnted: ', xWant, ', got upper bound: ', x(n)

            i = n
          endif

          !  if i is still zero find i less than xWant
          if (i .eq. 0) then
             do j = 1, n
                if(xWant .ge. x(j)) i = j
             enddo
          endif

          !  slope
          if (i .eq. n) then
            slope = 0.
          else
            slope = (y(i+1,k)-y(i,k)) / (x(i+1)-x(i))
          endif
          dx = xWant - x(i)
          yWant%rh(k) = y(i,k) + slope*dx

          yErr = 0.
        enddo
     end subroutine polint

   end function GOCART2G_MieCreate
   
!
! QueryByWave subroutines
!
#define BYWAVE_ 

#define RANK_ 1
#include "MieQuery.H"
#undef RANK_

#define RANK_ 2
#include "MieQuery.H"
#undef RANK_

#define RANK_ 3
#include "MieQuery.H"
#undef RANK_

#undef BYWAVE_

!
! QueryByChannel subroutines
!

#define BYCHANNEL_ 

#define RANK_ 1
#include "MieQuery.H"
#undef RANK_

#define RANK_ 2
#include "MieQuery.H"
#undef RANK_

#define RANK_ 3
#include "MieQuery.H"
#undef RANK_

#undef BYCHANNEL_

  integer function getChannel(this, wavelength, rc) result (channel)
     class (GOCART2G_Mie), intent(in) :: this
     real, intent(in) :: wavelength
     integer, optional, intent(out) :: rc
     real, parameter :: w_tol = 1.e-9
     integer :: i

     channel = -1
     do i = 1, this%nch
       if (abs(this%wavelengths(i)-wavelength) <= w_tol) then
          channel = i
          exit
       endif
    enddo

    if (present(rc)) rc = 0

    if (channel < 0) then
       !$omp critical
       print*, "wavelength of ",wavelength, " is an invalid value."
       !$omp end critical
       if (present(rc)) rc = -1
    endif

  end function getChannel

  real function getWavelength(this, channel, rc) result (wavelength)
     class (GOCART2G_Mie), intent(in) :: this
     integer, intent(in) :: channel
     integer, optional, intent(out) :: rc
     real, parameter :: w_tol = 1.e-9
     integer :: i

     if (present(rc)) rc = 0

     if (channel <=0 .or. channel > this%nch ) then
       !$omp critical
       print*, "The channel of ",channel, " is an invalid channel number."
       !$omp end critical
       if (present(rc)) rc = -1
       wavelength = -1. ! meanlingless nagative
       return
     endif

     wavelength = this%wavelengths(channel)

  end function getWavelength

  elemental real function interp(table, irh, arh)
    type(RH_Mie), intent(in) :: table
    integer, intent(in) :: irh
    real, intent(in) :: arh
    interp = sum(table%rh(irh:irh+1) * [(1-arh),arh])
  end function interp

end module GOCART2G_MieMod
