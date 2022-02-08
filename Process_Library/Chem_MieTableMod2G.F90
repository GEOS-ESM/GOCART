!BOP
!
! !MODULE:  Chem_MieTableMod --- Reader for aerosol mie tables
!
! !INTERFACE:
!

   module Chem_MieTableMod2G

! !USES:
   implicit none
   include "netcdf.inc"    ! Required for Mie tables stored as NCDF files

! !PUBLIC TYPES:
!
   private
   public  Chem_MieTable        ! Holds Mie Lookup Tables
   public  Chem_Mie             ! Holds Chem_MieTable and other information

!
! !PUBLIC MEMBER FUNCTIONS:
!
   public  Chem_MieTableCreate  ! Constructor 
   public  Chem_MieTableRead    ! Read the mie table from the file
   public  Chem_MieQuery
   public  Chem_MieTableIndex
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
!EOP
!-------------------------------------------------------------------------
!
! Mie LUT table
! Will be reduced from input files to the desired channels
! --------

  type Chem_MieTable
     character(len=255) :: mietablename
     integer :: nlambda         ! number of wavelengths in table
     integer :: nrh             ! number of RH values in table
     integer :: nbin            ! number of size bins in table
     integer :: nMom            ! number of moments of phase function
     integer :: nPol            ! number of elements of scattering phase matrix
     real, pointer    :: lambda(:) => null()       ! wavelengths [m]
     real, pointer    :: rh(:) => null()           ! RH values   [fraction]
     real, pointer    :: reff(:,:) => null()       ! effective radius [m]
     real, pointer    :: bext(:,:,:) => null()     ! bext values [m2 kg-1]
     real, pointer    :: bsca(:,:,:) => null()     ! bsca values [m2 kg-1]
     real, pointer    :: bbck(:,:,:) => null()     ! bbck values [m2 kg-1]
     real, pointer    :: g(:,:,:) => null()        ! asymmetry parameter
     real, pointer    :: pback(:,:,:,:) => null()  ! Backscatter phase function
     real, pointer    :: pmom(:,:,:,:,:) => null( )  ! moments of phase function
     real, pointer    :: gf(:,:) => null()         ! hygroscopic growth factor
     real, pointer    :: rhop(:,:) => null()       ! wet particle density [kg m-3]
     real, pointer    :: rhod(:,:) => null()       ! wet particle density [kg m-3]
     real, pointer    :: vol(:,:) => null()        ! wet particle volume [m3 kg-1]
     real, pointer    :: area(:,:) => null()       ! wet particle cross section [m2 kg-1]
     real, pointer    :: refr(:,:,:) => null()     ! real part of refractive index
     real, pointer    :: refi(:,:,:) => null()     ! imaginary part of refractive index

     integer          :: rhi(991)        ! pointer to rh map
     real             :: rha(991)        ! slope on rh map
  end type Chem_MieTable

  type Chem_Mie
     integer :: nch                               ! number of channels
     integer :: nMom=0                            ! number of moments (phase function)
     integer :: nPol=0                            ! number of moments (phase function)
     real, pointer    :: channels(:)              ! wavelengths

     character(len=255) :: rcfile
     character(len=255) :: optics_file
                                                  ! mie tables -- dim(nch,nrh,nbin)
     type(Chem_MieTable), pointer :: mie_aerosol  => null()

     integer :: nq                                ! number of tracers
     character(len=255), pointer  :: vname(:)  => null()      ! possibly remove lines 67-71 vname,vindex,vtable
     integer, pointer             :: vindex(:) => null()
     type(Chem_MieTable), pointer :: vtable(:) => null()
                                                  ! mapping of vtable for given idx
     type(Chem_MieTable), pointer :: vtableUse => null()
  end type Chem_Mie

# ifndef HAS_NETCDF3
  external nf_open, nf_inq_dimid, nf_inq_dimlen, nf_inq_varid, &
           nf_get_var_double, nf_close
#endif


  interface Chem_MieQuery
     module procedure Chem_MieQueryByInt
     module procedure Chem_MieQueryByIntWithpmom  !possibly remove
  end interface


CONTAINS

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Chem_MieTableCreate --- Construct Chemistry Registry
!
! !INTERFACE:
!

  Function Chem_MieTableCreate ( rcfile, rc )

  implicit none
  type(Chem_MieTable) Chem_MieTableCreate

! !INPUT PARAMETERS:

   character(len=*) :: rcfile  ! Mie table file name

! !OUTPUT PARAMETERS:

   integer, intent(out) ::  rc            ! Error return code:
                                          !  0 - all is well
                                          !  1 - 

! !DESCRIPTION:
!
!
! !REVISION HISTORY:
!
!  09Mar2005 da Silva  API, prologues.
!
!EOP
!-------------------------------------------------------------------------
   character(len=*), parameter ::  myname = 'Chem_MieTableCreate'

   type(Chem_MieTable) :: this

!                           _Iam_("Chem_MieTableCreate")

   rc = 0

   this%mietablename = rcfile

!  Note: The actual allocation is done when reading because dimensions are
!        read from file

!  All done
!  --------
   Chem_MieTableCreate = this

   return

 end Function Chem_MieTableCreate


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Chem_MieTableRead --- Read and fill in the Mie table, interpolated
!                              to the requested channels
!
! !INTERFACE:
!
   SUBROUTINE Chem_MieTableRead ( this, nch, channels, rc, nmom )

! !INPUT PARAMETERS:

   IMPLICIT none
   TYPE(Chem_MieTable), intent(inout)       :: this
   integer, intent(in) :: nch               ! number of channels to interpolate table to
   real, intent(in)    :: channels(:)       ! channels to interpolate table to
   integer, OPTIONAL, intent(in) :: nmom    ! number of moments to keep (default=0)
   integer, intent(out) :: rc   ! return code


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

      character(len=*), parameter :: myname = 'Chem_MieTableRead'

      integer :: ncid, idimid, ivarid, n, i, j, ip1
      integer :: nch_table, nrh_table, nbin_table, nmom_table, nPol_table
!     Tables are hard-wired as single precision
      real*8, pointer :: channels_table(:), rh_table(:), reff_table(:,:), &
                         bext_table(:,:,:), bsca_table(:,:,:), &
                         bbck_table(:,:,:), g_table(:,:,:), &
                         pmom_table(:,:,:,:,:), pback_table(:,:,:,:), &
                         gf_table(:,:), rhop_table(:,:), rhod_table(:,:), &
                         vol_table(:,:), area_table(:,:), &
                         refr_table(:,:,:), refi_table(:,:,:)

      real :: yerr
      integer :: nmom_, imom, ipol
      integer :: status

#define NF_VERIFY_(expr) rc = expr; if (rc /= 0) return
#define __NF_STAT__ stat=status); NF_VERIFY_(status

      rc = 0

!     Whether or not doing phase function
!     -----------------------------------
      if ( present(nmom) ) then
         nmom_ = nmom
      else
         nmom_ = 0
      end if

!     Set up nPol_table for reading backscatter phase function
!     This will get overwritten if pmoments is requested
!     --------------------------------------------------------
      nPol_table = 6

!     Open the table and get the dimensions
!     -------------------------------------
      rc = nf_open(this%mietablename, NF_NOWRITE, ncid)
      IF ( rc /= 0 ) THEN
        print *, 'nf_open '//this%mietablename//'  RETURN CODE=', rc
        NF_VERIFY_(rc)
      END IF

!     RH
!     --
      NF_VERIFY_(nf_inq_dimid(ncid,'rh',idimid))
      NF_VERIFY_(nf_inq_dimlen(ncid,idimid,nrh_table))

!     Channels
!     --------
      NF_VERIFY_(nf_inq_dimid(ncid,'lambda',idimid))
      NF_VERIFY_(nf_inq_dimlen(ncid,idimid,nch_table))

!     Dry Effective radius
!     --------------------
      NF_VERIFY_(nf_inq_dimid(ncid,'radius',idimid))
      NF_VERIFY_(nf_inq_dimlen(ncid,idimid,nbin_table))

!     Moments of phase function
!     -------------------------
      if ( nmom_ > 0 ) then
         NF_VERIFY_(nf_inq_dimid(ncid,'nMom',idimid))
         NF_VERIFY_(nf_inq_dimlen(ncid,idimid,nmom_table))
         if ( nmom_ > nmom_table ) then
!            rc = 99
            print*,'Error: nmom_ > nmom_table, see:'//myname
            NF_VERIFY_(1)
         end if
         NF_VERIFY_(nf_inq_dimid(ncid,'nPol',idimid))
         NF_VERIFY_(nf_inq_dimlen(ncid,idimid,nPol_table))
      endif

!     Get the table contents
!     -------------------------------------
!      allocate ( channels_table(nch_table), rh_table(nrh_table), &
!                bext_table(nch_table,nrh_table,nbin_table), &
!                bsca_table(nch_table,nrh_table,nbin_table), &
!                bbck_table(nch_table,nrh_table,nbin_table), &
!                g_table(nch_table,nrh_table,nbin_table), stat = rc )

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
      NF_VERIFY_(nf_inq_varid(ncid,'lambda',ivarid))
      NF_VERIFY_(nf_get_var_double(ncid,ivarid,channels_table))
      NF_VERIFY_(nf_inq_varid(ncid,'rEff',ivarid))
      NF_VERIFY_(nf_get_var_double(ncid,ivarid,reff_table))
      NF_VERIFY_(nf_inq_varid(ncid,'bext',ivarid))
      NF_VERIFY_(nf_get_var_double(ncid,ivarid,bext_table))
      NF_VERIFY_(nf_inq_varid(ncid,'bsca',ivarid))
      NF_VERIFY_(nf_get_var_double(ncid,ivarid,bsca_table))
      NF_VERIFY_(nf_inq_varid(ncid,'bbck',ivarid))
      NF_VERIFY_(nf_get_var_double(ncid,ivarid,bbck_table))
      NF_VERIFY_(nf_inq_varid(ncid,'g',ivarid))
      NF_VERIFY_(nf_get_var_double(ncid,ivarid,g_table))
      NF_VERIFY_(nf_inq_varid(ncid,'rh',ivarid))
      NF_VERIFY_(nf_get_var_double(ncid,ivarid,rh_table))

      ! TODO: we need to look at these NF_NOERR checks
!     Get the backscatter phase function values
      rc = nf_inq_varid(ncid,'pback',ivarid)
      if(rc .ne. NF_NOERR) then   ! pback not in table, fill in dummy variable
        pback_table = 1.
      else
        NF_VERIFY_(nf_get_var_double(ncid,ivarid,pback_table))
      endif

      if ( nmom_ > 0 ) then
         NF_VERIFY_(nf_inq_varid(ncid,'pmom',ivarid))
         NF_VERIFY_(nf_get_var_double(ncid,ivarid,pmom_table))
      end if

!     Aerosol optical properties not necessarily stored in all versions of the tables
!     ----------------------
!     Particle growth factor
      rc = nf_inq_varid(ncid,'growth_factor',ivarid)
      if(rc .ne. NF_NOERR) then   ! not in table, fill in dummy variable
        gf_table = -999.
      else
        NF_VERIFY_(nf_get_var_double(ncid,ivarid,gf_table))
      endif

!     Wet particle density
      rc = nf_inq_varid(ncid,'rhop',ivarid)
      if(rc .ne. NF_NOERR) then   ! not in table, fill in dummy variable
        rhop_table = -999.
      else
        NF_VERIFY_(nf_get_var_double(ncid,ivarid,rhop_table))
      endif

!     Dry particle density (will be pulled from wet particle radius)
      rc = nf_inq_varid(ncid,'rhop',ivarid)
      if(rc .ne. NF_NOERR) then   ! not in table, fill in dummy variable
        rhod_table = -999.
      else
        rc = nf_get_var_double(ncid,ivarid,rhod_table)
        do i = 1, nrh_table
          rhod_table(i,:) = rhod_table(1,:)
        enddo
        if (rc /=0) return
      endif

!     Wet particle real part of refractive index
      rc = nf_inq_varid(ncid,'refreal',ivarid)
      if(rc .ne. NF_NOERR) then   ! not in table, fill in dummy variable
        refr_table = -999.
      else
        NF_VERIFY_(nf_get_var_double(ncid,ivarid,refr_table))
      endif

!     Wet particle imaginary part of refractive index (ensure positive)
      rc = nf_inq_varid(ncid,'refimag',ivarid)
      if(rc .ne. NF_NOERR) then   ! not in table, fill in dummy variable
        refi_table = -999.
      else
        NF_VERIFY_(nf_get_var_double(ncid,ivarid,refi_table))
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
      NF_VERIFY_(nf_close(ncid))


!     Setup the table to be returned
!     -------------------------------------
      this%nlambda = nch
      this%nrh = nrh_table
      this%nbin = nbin_table
      this%nMom = nmom_
!      if ( nmom_ > 0 ) this%nPol = nPol_table
      this%nPol = nPol_table

!      allocate ( this%lambda(this%nLambda), this%rh(this%nrh), &
!                 this%bext(this%nLambda,this%nrh,this%nbin),   &
!                 this%bsca(this%nLambda,this%nrh,this%nbin),   &
!                 this%bbck(this%nLambda,this%nrh,this%nbin),   &
!                 this%g(this%nLambda,this%nrh,this%nbin),      &
!                 stat = rc )

      allocate (this%lambda(this%nLambda), __NF_STAT__)
      allocate (this%rh(this%nrh), __NF_STAT__)
      allocate (this%reff(this%nrh,this%nbin), __NF_STAT__)
      allocate (this%bext(this%nrh,this%nLambda,this%nbin), __NF_STAT__)
      allocate (this%bsca(this%nrh,this%nLambda,this%nbin), __NF_STAT__)
      allocate (this%bbck(this%nrh,this%nLambda,this%nbin), __NF_STAT__)
      allocate (this%g(this%nrh,this%nLambda,this%nbin),    __NF_STAT__)
      allocate (this%pback(this%nrh,this%nLambda,this%nbin,this%nPol),    __NF_STAT__)
      if ( nmom_ > 0 ) then
         allocate (this%pmom(this%nrh,this%nLambda,this%nbin,this%nMom,this%nPol),    __NF_STAT__)
      end if
      allocate (this%gf(this%nrh,this%nbin),    __NF_STAT__)
      allocate (this%rhop(this%nrh,this%nbin),    __NF_STAT__)
      allocate (this%rhod(this%nrh,this%nbin),    __NF_STAT__)
      allocate (this%vol(this%nrh,this%nbin),    __NF_STAT__)
      allocate (this%area(this%nrh,this%nbin),    __NF_STAT__)
      allocate (this%refr(this%nrh,this%nLambda,this%nbin), __NF_STAT__)
      allocate (this%refi(this%nrh,this%nLambda,this%nbin), __NF_STAT__)

!     Preserve the full RH structure of the input table
      this%rh(:) = rh_table(:)

!     Insert the requested channels in the output table
      this%lambda(:) = channels(:)

!     Insert rEff (moist effective radius)
      this%reff(:,:) = reff_table(:,:)

!     Now we linearly interpolate the input table to the output table grid
!     of requested channels
      do j = 1, this%nbin
       do i = 1, this%nrh
        do n = 1, this%nlambda
         call polint(channels_table,bext_table(:,i,j),nch_table, &
                     this%lambda(n),this%bext(i,n,j),yerr)
         call polint(channels_table,bsca_table(:,i,j),nch_table, &
                     this%lambda(n),this%bsca(i,n,j),yerr)
         call polint(channels_table,bbck_table(:,i,j),nch_table, &
                     this%lambda(n),this%bbck(i,n,j),yerr)
         call polint(channels_table,g_table(:,i,j),nch_table, &
                     this%lambda(n),this%g(i,n,j),yerr)
         call polint(channels_table,refr_table(:,i,j),nch_table, &
                     this%lambda(n),this%refr(i,n,j),yerr)
         call polint(channels_table,refi_table(:,i,j),nch_table, &
                     this%lambda(n),this%refi(i,n,j),yerr)
         do ipol = 1, this%nPol
                  call polint(channels_table,pback_table(:,i,j,ipol),nch_table, &
                       this%lambda(n),this%pback(i,n,j,ipol),yerr)
         end do
         if ( nmom_ > 0 ) then
            do imom = 1, this%nMom
               do ipol = 1, this%nPol
                  call polint(channels_table,pmom_table(:,i,j,imom,ipol),nch_table, &
                       this%lambda(n),this%pmom(i,n,j,imom,ipol),yerr)
               end do
            end do
         end if
        enddo
       enddo
      enddo

!     Insert growth factor
      this%gf(:,:) = gf_table(:,:)

!     Wet particle density [kg m-3]
      this%rhop(:,:) = rhop_table(:,:)

!     Dry particle density [kg m-3]
      this%rhod(:,:) = rhod_table(:,:)

!     Volume [m3 kg-1]
      this%vol(:,:) = vol_table(:,:)

!     Area [m2 kg-1]
      this%area(:,:) = area_table(:,:)

!     Now we do a mapping of the RH from the input table to some high
!     resolution representation.  This is to spare us the need to
!     do a full-up interpolation later on.
!     RH input from the table is scaled 0 - 0.99
!     We resolve the map to 0 - 0.990 in steps of 0.001 (991 total steps)
      do j = 1, 991
       do i = this%nrh, 1, -1
        if( (j-1) .ge. int(this%rh(i)*1000)) then
         ip1 = i + 1
         this%rhi(j) = i
         if(ip1 .gt. this%nrh) then
          this%rha(j) = 0.
         else
          this%rha(j) =   ( (j-1)/1000. - this%rh(i)) &
                       /  ( this%rh(ip1)- this%rh(i))
         endif
         exit
        endif
       enddo
!       print *, j, this%rhi(j), this%rha(j), this%rh(this%rhi(j))
      enddo

!      deallocate (channels_table, rh_table, bext_table, bsca_table, &
!                  bbck_table, g_table, stat = rc )

      deallocate (channels_table, __NF_STAT__)
      deallocate (rh_table, __NF_STAT__)
      deallocate (reff_table, __NF_STAT__)
      deallocate (bext_table, __NF_STAT__)
      deallocate (bsca_table, __NF_STAT__)
      deallocate (bbck_table, __NF_STAT__)
      deallocate (g_table, __NF_STAT__)
      deallocate (pback_table, __NF_STAT__)
      if ( nmom_ > 0 ) then
         deallocate (pmom_table, __NF_STAT__)
      endif
      deallocate (gf_table, __NF_STAT__)
      deallocate (rhop_table, __NF_STAT__)
      deallocate (rhod_table, __NF_STAT__)
      deallocate (vol_table, __NF_STAT__)
      deallocate (area_table, __NF_STAT__)
      deallocate (refr_table, __NF_STAT__)
      deallocate (refi_table, __NF_STAT__)

return

contains

   subroutine polint(x,y,n,xWant,yWant,yErr)
   integer :: n
!  recall, table hard-wired single precision
   real*8 :: x(n),y(n)
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
   if(i .eq. 0) then
    do j = 1, n
     if(xWant .ge. x(j)) i = j
    enddo
   endif

!  slope
   if(i .eq. n) then
    slope = 0.
   else
    slope = (y(i+1)-y(i)) / (x(i+1)-x(i))
   endif
   dx = xWant - x(i)
   yWant = y(i) + slope*dx

   yErr = 0.

   return
   end subroutine polint

END SUBROUTINE Chem_MieTableRead
!----------------------------------------------------------------------------

!BOP
!
! !IROUTINE:  Chem_MieQueryByInt --- Return Tau, SSA, etc 
!
!
! !INTERFACE:
!
   impure elemental subroutine Chem_MieQueryByInt ( this, idx, channel, q_mass, rh,     &
                                   tau, ssa, gasym, bext, bsca, bbck,  &
                                   reff, p11, p22, gf, rhop, rhod, &
                                   vol, area, refr, refi, rc )

! !INPUT PARAMETERS:

   type(Chem_Mie), target, intent(in ) :: this
   integer,                intent(in ) :: idx     ! variable index on Chem_Mie
   real,                   intent(in ) :: channel ! channel number 
   real,                   intent(in ) :: q_mass  ! aerosol mass [kg/m2],
   real,                   intent(in ) :: rh      ! relative himidity

! !OUTPUT PARAMETERS:

   real,    optional,      intent(out) :: tau   ! aerol extinction optical depth
   real,    optional,      intent(out) :: ssa   ! single scattering albedo
   real,    optional,      intent(out) :: gasym ! asymmetry parameter
   real,    optional,      intent(out) :: bext  ! mass extinction efficiency [m2 (kg dry mass)-1]
   real,    optional,      intent(out) :: bsca  ! mass scattering efficiency [m2 (kg dry mass)-1]
   real,    optional,      intent(out) :: bbck  ! mass backscatter efficiency [m2 (kg dry mass)-1]
   real,    optional,      intent(out) :: reff  ! effective radius (micron)
   real,    optional,      intent(out) :: p11   ! P11 phase function at backscatter
   real,    optional,      intent(out) :: p22   ! P22 phase function at backscatter
   real,    optional,      intent(out) :: gf    ! Growth factor (ratio of wet to dry radius)
   real,    optional,      intent(out) :: rhop  ! Wet particle density [kg m-3]
   real,    optional,      intent(out) :: rhod  ! Dry particle density [kg m-3]
   real,    optional,      intent(out) :: vol   ! Wet particle volume [m3 kg-1]
   real,    optional,      intent(out) :: area  ! Wet particle cross section [m2 kg-1]
   real,    optional,      intent(out) :: refr  ! Wet particle real part of ref. index
   real,    optional,      intent(out) :: refi  ! Wet particle imag. part of ref. index
   integer, optional,      intent(out) :: rc    ! error code

! !DESCRIPTION:
!
!   Returns requested parameters from the Mie tables, as a function 
!   of species, relative humidity, and channel
!
!  Notes: Needs some checking, and I still force an interpolation step

!
! !REVISION HISTORY:
!
!  23Mar2005 Colarco
!  11Jul2005 da Silva   Standardization.
!
!EOP
!-------------------------------------------------------------------------


      integer                      :: ICHANNEL, TYPE
      integer                      :: irh, irhp1, isnap
      real                         :: rhUse, arh
      real                         :: bextIn, bscaIn, bbckIn, gasymIn, p11In, p22In, &
                                      gfIn, rhopIn, rhodIn, volIn, areaIn, &
                                      refrIn, refiIn
      type(Chem_MieTable), pointer :: TABLE

      character(len=*), parameter  :: Iam = 'Chem_MieQuery'

      if ( present(rc) ) rc = 0

      ICHANNEL = nint(CHANNEL)
!      TABLE => this%vtableUse
      Table => this%mie_aerosol
      TYPE = idx

!      ASSERT_(TYPE>0)
!      ASSERT_(ICHANNEL>=LBOUND(TABLE%bext,1))
!      ASSERT_(ICHANNEL<=UBOUND(TABLE%bext,1))

!     Now map the input RH to the high resolution hash table for RH
      rhUse = max(rh,0.)
      rhUse = min(rh,0.99)
      isnap = int((rhUse+0.001)*1000.)
      if(isnap .lt. 1) isnap = 1
      arh   = TABLE%rha( isnap )
      irh   = TABLE%rhi( isnap )
      irhp1 = irh+1
      if(irhp1 .gt. TABLE%nrh) irhp1 = TABLE%nrh

!     Now linearly interpolate the input table for the requested aerosol and
!     channel; rh is the relative humidity.

      if(present(bext) .or. present(tau) .or. present(ssa) ) then
         bextIn =   TABLE%bext(irh  ,ichannel,TYPE) * (1.-arh) &
                  + TABLE%bext(irhp1,ichannel,TYPE) * arh
      endif

      if(present(bsca) .or. present(ssa) ) then
         bscaIn =   TABLE%bsca(irh  ,ichannel,TYPE) * (1.-arh) &
                  + TABLE%bsca(irhp1,ichannel,TYPE) * arh
      endif

      if(present(bbck)) then
         bbckIn =   TABLE%bbck(irh  ,ichannel,TYPE) * (1.-arh) &
                  + TABLE%bbck(irhp1,ichannel,TYPE) * arh
      endif

      if(present(gasym)) then
         gasymIn =  TABLE%g(irh  ,ichannel,TYPE) * (1.-arh) &
                  + TABLE%g(irhp1,ichannel,TYPE) * arh
      endif

      if(present(rEff) ) then
         rEff =     TABLE%rEff(irh  ,TYPE) * (1.-arh) &
                  + TABLE%rEff(irhp1,TYPE) * arh
         rEff = 1.E6 * rEff ! convert to microns
      endif

!      if(present(pmom)) then
!         pmom(:,:) = TABLE%pmom(irh  ,ichannel,TYPE,:,:) * (1.-arh) &
!                   + TABLE%pmom(irhp1,ichannel,TYPE,:,:) * arh
!      endif

!      if(present(pmom)) then
!         call Chem_MieQueryByIntWithpmom(this, idx, channel, q_mass, rh, pmom)
!      endif


      if(present(p11) ) then
         p11In =   TABLE%pback(irh  ,ichannel,TYPE,1) * (1.-arh) &
                 + TABLE%pback(irhp1,ichannel,TYPE,1) * arh
      endif

      if(present(p22) ) then
         p22In =   TABLE%pback(irh  ,ichannel,TYPE,5) * (1.-arh) &
                 + TABLE%pback(irhp1,ichannel,TYPE,5) * arh
      endif

      if(present(gf) ) then
         gfIn =     TABLE%gf(irh  ,TYPE) * (1.-arh) &
                  + TABLE%gf(irhp1,TYPE) * arh
      endif

      if(present(rhod) ) then
         rhodIn =   TABLE%rhod(1  ,TYPE)
      endif

      if(present(vol) ) then
         volIn  =   TABLE%vol(irh  ,TYPE) * (1.-arh) &
                  + TABLE%vol(irhp1,TYPE) * arh
      endif

      if(present(area) ) then
         areaIn  =   TABLE%area(irh  ,TYPE) * (1.-arh) &
                  + TABLE%area(irhp1,TYPE) * arh
      endif

      if(present(refr) .or. present(tau) .or. present(ssa) ) then
         refrIn =   TABLE%refr(irh  ,ichannel,TYPE) * (1.-arh) &
                  + TABLE%refr(irhp1,ichannel,TYPE) * arh
      endif

      if(present(refi) .or. present(tau) .or. present(ssa) ) then
         refiIn =   TABLE%refi(irh  ,ichannel,TYPE) * (1.-arh) &
                  + TABLE%refi(irhp1,ichannel,TYPE) * arh
      endif

!     Fill the requested outputs
      if(present(tau  )) tau   = bextIn * q_mass
      if(present(ssa  )) ssa   = bscaIn/bextIn
      if(present(bext )) bext  = bextIn
      if(present(bsca )) bsca  = bscaIn
      if(present(bbck )) bbck  = bbckIn
      if(present(gasym)) gasym = gasymIn
      if(present(p11  )) p11   = p11In
      if(present(p22  )) p22   = p22In
      if(present(gf   )) gf    = gfIn
      if(present(rhop )) rhop  = rhopIn
      if(present(rhod )) rhod  = rhodIn
      if(present(vol ))  vol   = volIn
      if(present(area )) area  = areaIn
      if(present(refr )) refr  = refrIn
      if(present(refi )) refi  = refiIn

!  All Done
!----------
  end subroutine Chem_MieQueryByInt


!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Chem_MieQueryByIntWithpmom --- Return Tau, SSA, etc 
!
!
! !INTERFACE:
!
   subroutine Chem_MieQueryByIntWithpmom ( this, idx, channel, q_mass, rh,     &
                                   tau, ssa, gasym, bext, bsca, bbck,  &
                                   reff, pmom, p11, p22, gf, rhop, rhod, &
                                   vol, area, refr, refi, rc )

! !INPUT PARAMETERS:

   type(Chem_Mie), target, intent(in ) :: this
   integer,                intent(in ) :: idx     ! variable index on Chem_Mie
   real,                   intent(in ) :: channel ! channel number
   real,                   intent(in ) :: q_mass  ! aerosol mass [kg/m2],
   real,                   intent(in ) :: rh      ! relative himidity

! !OUTPUT PARAMETERS:

   real,    optional,      intent(out) :: tau   ! aerol extinction optical depth
   real,    optional,      intent(out) :: ssa   ! single scattering albedo
   real,    optional,      intent(out) :: gasym ! asymmetry parameter
   real,    optional,      intent(out) :: bext  ! mass extinction efficiency [m2 (kg dry mass)-1]
   real,    optional,      intent(out) :: bsca  ! mass scattering efficiency [m2 (kg dry mass)-1]
   real,    optional,      intent(out) :: bbck  ! mass backscatter efficiency [m2 (kg dry mass)-1]
   real,    optional,      intent(out) :: reff  ! effective radius (micron)
   real,                   intent(out) :: pmom(:,:)
   real,    optional,      intent(out) :: p11   ! P11 phase function at backscatter
   real,    optional,      intent(out) :: p22   ! P22 phase function at backscatter
   real,    optional,      intent(out) :: gf    ! Growth factor (ratio of wet to dry radius)
   real,    optional,      intent(out) :: rhop  ! Wet particle density [kg m-3]
   real,    optional,      intent(out) :: rhod  ! Dry particle density [kg m-3]
   real,    optional,      intent(out) :: vol   ! Wet particle volume [m3 kg-1]
   real,    optional,      intent(out) :: area  ! Wet particle cross section [m2 kg-1]
   real,    optional,      intent(out) :: refr  ! Wet particle real part of ref. index
   real,    optional,      intent(out) :: refi  ! Wet particle imag. part of ref. index
   integer, optional,      intent(out) :: rc    ! error code

! !DESCRIPTION:
!
!   Returns requested parameters from the Mie tables, as a function 
!   of species, relative humidity, and channel
!
!  Notes: Needs some checking, and I still force an interpolation step

!
! !REVISION HISTORY:
!
!  23Mar2005 Colarco
!  11Jul2005 da Silva   Standardization.
!
!EOP
!-------------------------------------------------------------------------


      integer                      :: ICHANNEL, TYPE
      integer                      :: irh, irhp1, isnap
      real                         :: rhUse, arh
      type(Chem_MieTable), pointer :: TABLE

      character(len=*), parameter  :: Iam = 'Chem_MieQueryByIntWithpmom'

      if ( present(rc) ) rc = 0

      ICHANNEL = nint(CHANNEL)
      TABLE => this%vtableUse
      TYPE = idx

!     Now map the input RH to the high resolution hash table for RH
      rhUse = max(rh,0.)
      rhUse = min(rh,0.99)
      isnap = int((rhUse+0.001)*1000.)
      if(isnap .lt. 1) isnap = 1
      arh   = TABLE%rha( isnap )
      irh   = TABLE%rhi( isnap )
      irhp1 = irh+1
      if(irhp1 .gt. TABLE%nrh) irhp1 = TABLE%nrh

!     Now linearly interpolate the input table for the requested aerosol and
!     channel; rh is the relative humidity.
    call Chem_MieQuery ( this, idx, channel, q_mass, rh, &
                         tau, ssa, gasym, bext, bsca, bbck, &
                         reff, p11, p22, gf, rhop, rhod, &
                         vol, area, refr, refi, rc )
    NF_VERIFY_(rc)

         pmom(:,:) = TABLE%pmom(irh  ,ichannel,TYPE,:,:) * (1.-arh) &
                   + TABLE%pmom(irhp1,ichannel,TYPE,:,:) * arh



!  All Done
!----------
  end subroutine Chem_MieQueryByIntWithpmom

  function Chem_MieTableIndex(wavelength, channels, rc) result (mieTable_index)
    real :: mieTable_index
    real, intent(in) :: wavelength
    real, intent(in) :: channels(:)
    integer, optional, intent(out) :: rc
    real, parameter :: width = 1.e-9
    integer :: i
    
    mieTable_index = -1.
    do i = 1, size(channels)
       if (channels[i]-width <= wavelength .and. wavelength <= channels[i]+width) then
          mieTable_index = i
          exit
       endif
    enddo
    if (present(rc)) rc = 0
    if (mieTable_index < 0)
       print*, "wavelength of ",wavelength, " is an invalid value."
       if (present(rc)) rc = -1
    endif

  end function Chem_MieTableIndex

 end module Chem_MieTableMod2G

