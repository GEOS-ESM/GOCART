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

   type GOCART2G_Mie
      
      character(len=:), allocatable :: table_name
      integer :: nch             ! number of channels in table (replacement of nlamfda)
      integer :: nrh             ! number of RH values in table
      integer :: nbin            ! number of size bins in table
      integer :: nMom            ! number of moments of phase function
      integer :: nPol            ! number of elements of scattering phase matrix

                                            ! c=channel, r=rh, b=bin, m=moments, p=nPol
      real, pointer  :: wavelengths(:) => Null()  ! (c) wavelengths [m]
      real, pointer  :: rh(:) => Null()           ! (r) RH values   [fraction]
      real, pointer  :: reff(:,:) => Null()       ! (r,b) effective radius [m]
      real, pointer  :: bext(:,:,:) => Null()     ! (r,c,b) bext values [m2 kg-1]
      real, pointer  :: bsca(:,:,:) => Null()     ! (r,c,b) bsca values [m2 kg-1]
      real, pointer  :: bbck(:,:,:) => Null()     ! (r,c,b) bbck values [m2 kg-1]
      real, pointer  :: g(:,:,:) => Null()        ! (r,c,b) asymmetry parameter
!ams  real, pointer  :: pback(:,:,:,:) => Null()  ! (r,c,b,p) Backscatter phase function
      real, pointer  :: p11(:,:,:) => Null()      ! (r,c,b) Backscatter phase function, index 1 
      real, pointer  :: p22(:,:,:) => Null()      ! (r,c,b) Backscatter phase function, index 5
      real, pointer  :: pmom(:,:,:,:,:) => Null() ! (r,c,b,m,p) moments of phase function
      real, pointer  :: gf(:,:) => Null()         ! (r,b) hygroscopic growth factor
      real, pointer  :: rhop(:,:) => Null()       ! (r,b) wet particle density [kg m-3]
      real, pointer  :: rhod(:,:) => Null()       ! (r,b) wet particle density [kg m-3]
      real, pointer  :: vol(:,:) => Null()        ! (r,b) wet particle volume [m3 kg-1]
      real, pointer  :: area(:,:) => Null()       ! (r,b) wet particle cross section [m2 kg-1]
      real, pointer  :: refr(:,:,:) => Null()     ! (r,c,b) real part of refractive index
      real, pointer  :: refi(:,:,:) => Null()     ! (r,c,b) imaginary part of refractive index

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

  type(GOCART2G_Mie) function GOCART2G_MieCreate ( MieFile, wavelengths, nmom, rc ) result (this)

! !INPUT PARAMETERS:

     character(len=*), intent(in) :: MieFile  ! Mie table file name
     real, optional,   intent(in) :: wavelengths(:)
     integer, optional,intent(in) :: nmom

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

     real, allocatable ::   channels_table(:),    rh_table(:), reff_table(:,:),    &
                            bext_table(:,:,:),    bsca_table(:,:,:),               &
                            bbck_table(:,:,:),    g_table(:,:,:),                  &
                            pmom_table(:,:,:,:,:),pback_table(:,:,:,:),            &
                            gf_table(:,:),        rhop_table(:,:), rhod_table(:,:),&
                            vol_table(:,:),       area_table(:,:),                 &
                            refr_table(:,:,:),    refi_table(:,:,:)

     real, pointer  :: pback(:,:,:,:)  ! (r,c,b,p) Backscatter phase function
     
     real :: yerr
     integer :: nmom_, imom, ipol
     integer :: status

#define NF_VERIFY_(expr) rc = expr; if (rc /= 0) return
#define __NF_STAT__ stat=status); NF_VERIFY_(status

     rc = 0
     this%table_name = MieFile

      
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

      if (present(wavelengths) ) then
        nch = size(wavelengths)
      else
        nch = nch_table
     end if
      
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

      allocate (this%rh(this%nrh), __NF_STAT__)
      allocate (this%reff(this%nrh,this%nbin), __NF_STAT__)
      allocate (this%bext(this%nrh,this%nch,this%nbin), __NF_STAT__)
      allocate (this%bsca(this%nrh,this%nch,this%nbin), __NF_STAT__)
      allocate (this%bbck(this%nrh,this%nch,this%nbin), __NF_STAT__)
      allocate (this%g(this%nrh,this%nch,this%nbin),    __NF_STAT__)
      allocate (pback(this%nrh,this%nch,this%nbin,this%nPol),    __NF_STAT__)
      if ( nmom_ > 0 ) then
         allocate (this%pmom(this%nrh,this%nch,this%nbin,this%nMom,this%nPol),    __NF_STAT__)
      end if
      allocate (this%gf(this%nrh,this%nbin),    __NF_STAT__)
      allocate (this%rhop(this%nrh,this%nbin),    __NF_STAT__)
      allocate (this%rhod(this%nrh,this%nbin),    __NF_STAT__)
      allocate (this%vol(this%nrh,this%nbin),    __NF_STAT__)
      allocate (this%area(this%nrh,this%nbin),    __NF_STAT__)
      allocate (this%refr(this%nrh,this%nch,this%nbin), __NF_STAT__)
      allocate (this%refi(this%nrh,this%nch,this%nbin), __NF_STAT__)
      allocate (this%wavelengths(this%nch), __NF_STAT__)
      allocate (this%p11(this%nrh,this%nch,this%nbin),    __NF_STAT__)
      allocate (this%p22(this%nrh,this%nch,this%nbin),    __NF_STAT__)

!     Preserve the full RH structure of the input table
      this%rh = rh_table 

!     Insert rEff (moist effective radius)
      this%reff = reff_table

!     Insert growth factor
      this%gf = gf_table

!     Wet particle density [kg m-3]
      this%rhop = rhop_table

!     Dry particle density [kg m-3]
      this%rhod = rhod_table

!     Volume [m3 kg-1]
      this%vol  = vol_table

!     Area [m2 kg-1]
      this%area = area_table

!     Insert the requested channels in the output table
      if ( present(wavelengths) ) then
         this%wavelengths = wavelengths
      else
         this%wavelengths = channels_table 
      endif

!     Now we linearly interpolate the input table to the output table grid
!     of requested channels
      if ( present(wavelengths) ) then
         do j = 1, this%nbin
            do i = 1, this%nrh
               do n = 1, this%nch
                  call polint(channels_table,bext_table(:,i,j),nch_table, &
                       this%wavelengths(n),this%bext(i,n,j),yerr)
                  call polint(channels_table,bsca_table(:,i,j),nch_table, &
                       this%wavelengths(n),this%bsca(i,n,j),yerr)
                  call polint(channels_table,bbck_table(:,i,j),nch_table, &
                       this%wavelengths(n),this%bbck(i,n,j),yerr)
                  call polint(channels_table,g_table(:,i,j),nch_table,    &
                       this%wavelengths(n),this%g(i,n,j),yerr)
                  call polint(channels_table,refr_table(:,i,j),nch_table, &
                       this%wavelengths(n),this%refr(i,n,j),yerr)
                  call polint(channels_table,refi_table(:,i,j),nch_table, &
                       this%wavelengths(n),this%refi(i,n,j),yerr)
                  do ipol = 1, this%nPol
                      call polint(channels_table,pback_table(:,i,j,ipol),nch_table,    &
                             this%wavelengths(n),pback(i,n,j,ipol),yerr)
                  end do
                  if ( nmom_ > 0 ) then
                     do imom = 1, this%nMom
                        do ipol = 1, this%nPol
                           call polint(channels_table,pmom_table(:,i,j,imom,ipol),nch_table, &
                                this%wavelengths(n),this%pmom(i,n,j,imom,ipol),yerr)
                        enddo
                     enddo
                  endif
               enddo
            enddo
         enddo
      else !(no wavelength)
         !swap the order
         this%bext = reshape(bext_table, [nrh_table, nch, nbin_table],order =[2,1,3])
         this%bsca = reshape(bsca_table, [nrh_table, nch, nbin_table],order =[2,1,3])
         this%bbck = reshape(bbck_table, [nrh_table, nch, nbin_table],order =[2,1,3])
         this%g    = reshape(   g_table, [nrh_table, nch, nbin_table],order =[2,1,3])
         this%refr = reshape(refr_table, [nrh_table, nch, nbin_table],order =[2,1,3])
         this%refi = reshape(refi_table, [nrh_table, nch, nbin_table],order =[2,1,3])
         pback     = reshape(pback_table,[nrh_table, nch, nbin_table, npol_table],order =[2,1,3,4])
         if ( nmom_ > 0 ) then
           this%pmom = reshape(pmom_table,[nrh_table,nch, nbin_table, nmom_, npol_table], order = [2,1,3,4,5])
         endif
      endif

!     Pick p11, p12
      this%p11 = pback(:,:,:,1)
      this%p22 = pback(:,:,:,5)

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
      enddo

      return

  contains

     subroutine polint(x,y,n,xWant,yWant,yErr)
       integer :: n
       !  recall, table hard-wired single precision
       real   :: x(n),y(n)
       real   :: xWant, yWant, yErr

       !  given array x(n) of independent variables and array y(n) of dependent
       !  variables, compute the linear interpolated result yWant at xWant and return
       !  with a dummy error estimate yErr.  Hacked up from Numerical Recipes Chapter 3

       integer :: i, j
       real    :: dx, slope
       character(len=255) :: msg

       !  on out of bounds, set i to lower or upper limit
       i = 0
       if(xWant .lt. x(1)) then
         i = 1
       endif
       if(xWant .gt. x(n)) then
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
         slope = (y(i+1)-y(i)) / (x(i+1)-x(i))
       endif
       dx = xWant - x(i)
       yWant = y(i) + slope*dx

       yErr = 0.

     end subroutine polint

   end function GOCART2G_MieCreate
   
!
! Query subroutines
!

#define RANK_ 1
#include "MieQuery.H"
#undef RANK_

#define RANK_ 2
#include "MieQuery.H"
#undef RANK_

#define RANK_ 3
#include "MieQuery.H"
#undef RANK_

  integer function getChannel(this, wavelength, rc) result (ch)
     class (GOCART2G_Mie), intent(in) :: this
     real, intent(in) :: wavelength
     integer, optional, intent(out) :: rc
     real, parameter :: w_tol = 1.e-9
     integer :: i

     ch = -1
     do i = 1, this%nch
       if (abs(this%wavelengths(i)-wavelength) <= w_tol) then
          ch = i
          exit
       endif
    enddo

    if (present(rc)) rc = 0

    if (ch < 0) then
       !$omp critical (GetCha)
       print*, "wavelength of ",wavelength, " is an invalid value."
       !$omp end critical (GetCha)
       if (present(rc)) rc = -1
    endif

  end function getChannel

  real function getWavelength(this, ith_channel, rc) result (wavelength)
     class (GOCART2G_Mie), intent(in) :: this
     integer, intent(in) :: ith_channel
     integer, optional, intent(out) :: rc
     real, parameter :: w_tol = 1.e-9
     integer :: i

     if (present(rc)) rc = 0

     if (ith_channel <=0 .or. ith_channel > this%nch ) then
       !$omp critical (GetWav)
       print*, "The channel of ",ith_channel, " is an invalid channel number."
       !$omp end critical (GetWav)
       if (present(rc)) rc = -1
       wavelength = -1. ! meanlingless nagative
       return
     endif
     
     wavelength = this%wavelengths(ith_channel)

  end function getWavelength

end module GOCART2G_MieMod
