
#include "MAPL_Generic.h"

!-------------------------------------------------------------------------
!      NASA/GSFC, Global Modeling & Assimilation Office, Code 900.3      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  GOCART2G_AopMod --- Aerosol Optical Depth Calculator
!
! !INTERFACE:
!

   module  GOCART2G_AopMod

! !USES:


   use ESMF
   use MAPL

   use GOCART2G_MieMod

   implicit none

!
! !PUBLIC MEMBER FUNCTIONS:
!
   PUBLIC  GOCART2G_AopCalculator3d  ! Calculates 3D AOP profiles

   type MieTable
        type (GOCART2G_Mie), pointer :: Table => null()     !  Mie Table structure, etc
        integer                      :: bin_number          !  bin number for each specie
   end type MieTable
!
! !DESCRIPTION:
!
!  This module implements an Aerosol Optical Property calculator.
!
! !REVISION HISTORY:
!
!  29Nov2010 da Silva  Revamped to use Simple Bundle construct.
!  Jun2022   Buchard/da Silva  Use with GOCART2G.
!
!EOP
!-------------------------------------------------------------------------


CONTAINS


!-----------------------------------------------------------------------------
!BOP
 
! !IROUTINE: GOCART2G_AopCalculator3D ---- Compute 3D Aerosol Extinction
!
! !INTERFACE:
!
  subroutine GOCART2G_AopCalculator3d ( y, x, Mie, verbose, rc ) 


! !ARGUMENTS:

  type(ESMF_FieldBundle), intent(inout) :: y             ! Aerosol extinction
  type(MAPL_SimpleBundle), intent(in)    :: x             ! Aerosol mixing ratio
  type(MieTable),          intent(inout) :: Mie(:)        ! Mie tables for each tracer

  logical, OPTIONAL,       intent(in)    :: verbose
  integer, OPTIONAL,       intent(out)   :: rc

! !DESCRIPTION: Return Aerosol Extinction bundle.
!
!EOP
!-----------------------------------------------------------------------------
  real(ESMF_KIND_R4), allocatable :: delc(:,:,:), delm(:,:,:), delz(:,:,:), rh(:,:,:), tau(:,:,:), bck(:,:,:), aod_(:,:,:,:)
  real(ESMF_KIND_R4), pointer :: ptrtau(:,:,:,:), ptrext(:,:,:,:), ptrbck(:,:,:,:)
  type(ESMF_Field) :: Field
  integer :: iq, im, jm, km, i, j, k, n, idxTable, nch
  integer :: iTau, iExt, iBck, iRH, iRHO
  logical :: verbose_
  real, allocatable :: aod(:,:)
  logical :: has_tau, has_ext, has_bck
              __Iam__("GOCART2G_ExtCalculator")

  if ( present(verbose) ) then
     verbose_ = verbose .and. MAPL_AM_I_ROOT()
  else
     verbose_ = .FALSE.
  end if

  im = size(x%coords%lons,1)
  jm = size(x%coords%lons,2)
  km = size(x%coords%levs)

  iRH = MAPL_SimpleBundleGetIndex(x,'rh',rank=3,quiet=.true.,rc=STATUS)
  if ( STATUS /= 0 ) then
     iRH = MAPL_SimpleBundleGetIndex(x,'RH',rank=3,__RC__)
  end if

  iRHO = MAPL_SimpleBundleGetIndex(x,'airdens',rank=3,quiet=.true.,rc=STATUS)
  if ( STATUS /= 0 ) then
     iRHO = MAPL_SimpleBundleGetIndex(x,'AIRDENS',rank=3,__RC__)
  end if
 
  call ESMF_FieldBundleGet(y, 'tau', ispresent=has_tau, __RC__)
  call ESMF_FieldBundleGet(y, 'extinction', ispresent=has_ext, __RC__)
  call ESMF_FieldBundleGet(y, 'backscat', ispresent=has_bck, __RC__)

  if (has_tau) then
      call ESMF_FieldBundleGet(y, 'tau', Field=Field, __RC__)
      call ESMF_FieldGet(Field,0,Farrayptr=ptrtau,__RC__)
  endif
  if (has_ext) then
      call ESMF_FieldBundleGet(y, 'extinction', Field=Field, __RC__)
      call ESMF_FieldGet(Field,0,Farrayptr=ptrext,__RC__)
  endif
  if (has_bck) then
      call ESMF_FieldBundleGet(y, 'backscat', Field=Field, __RC__)
      call ESMF_FieldGet(Field,0,Farrayptr=ptrbck,__RC__)
  endif
!  iTau = MAPL_SimpleBundleGetIndex(y,'tau',       rank=3,quiet=.true.,rc=STATUS)
!  iExt = MAPL_SimpleBundleGetIndex(y,'extinction',rank=3,quiet=.true.,rc=STATUS)
!  iBck = MAPL_SimpleBundleGetIndex(y,'backscat',  rank=3,quiet=.true.,rc=STATUS)

! Consistency check
! -----------------
  if ( .not. associated(x%coords%lcv%delp) ) then
     __raise__(MAPL_RC_ERROR,"x%coords%lcv must have delp set")
  end if
!  if ( any (Mie(1)%Table%nch /= 1) ) then
!     __raise__(MAPL_RC_ERROR,"for now, Mie must have a single channel")
!  end if
  n = 1

! Loop over aerosol species
! -------------------------
!  if ( iTau>0 ) y%r3(iTau)%q = 0.0
!  if ( iExt>0 ) y%r3(iExt)%q = 0.0
!  if ( iBck>0 ) y%r3(iBck)%q = 0.0
  nch = size(Mie(2)%Table%wavelengths)
  allocate(delm(im,jm,km), delc(im,jm,km), delz(im,jm,km),rh(im,jm,km))
  allocate(tau(im,jm,km),bck(im,jm,km), aod_(im,jm,km,nch))
  aod_ = 0
!  ptrtau = 0
  
     do iq =1, x%n3d
     print*, 'iq species', iq, x%r3(iq)%name
        if (associated(Mie(iq)%Table)) then
                 do n = 1 , size(Mie(iq)%Table%wavelengths)
                    delm = x%r3(iq)%q * x%coords%lcv%delp / MAPL_GRAV
                    delc = x%r3(iq)%q * x%r3(iRHO)%q ! concentration
                    delz = x%coords%lcv%delp / (MAPL_GRAV*x%r3(iRHO)%q)
                     rh = x%r3(iRH)%q
                 
                    ! print*, 'Table name', Mie(iq)%Table%table_name, Mie(iq)%bin_number, Mie(iq)%Table%wavelengths(n)
                     call Mie(iq)%Table%Query(Mie(iq)%Table%wavelengths(n), Mie(iq)%bin_number, delm, rh, tau=tau, bbck=bck) 
                     print*, 'aod lev', iq, x%r3(iq)%name, tau(100,100,72), maxval(tau(:,:,72))
                     aod_(:,:,:,n) = aod_(:,:,:,n) + tau
                     print*, 'aod_', aod_(100,100,72,1), aod_(100,100,72,2)
                 if ( has_tau ) ptrtau(:,:,:,n) = ptrtau(:,:,:,n) + tau
                 if ( has_ext ) ptrext(:,:,:,n) = ptrext(:,:,:,n) + tau / delz
                 if ( has_bck ) ptrbck(:,:,:,n) = ptrbck(:,:,:,n) + bck * delc
                 enddo
         endif
     enddo
       print*, 'max val aod lev 72' maxval(ptrtau(:,:,72,2)
!     print*, 'tau total over species', y%r3(iTau)%q(100,100,:), y%r3(iTau)%q(130,200,:)
!     print*, 'tau total over species', shape(y%r3(iTau)%q)
     allocate(aod(im,jm))
     aod = 0
  
     do j = 1, 72
        aod = aod + aod_(:,:,j,2)
     enddo
     print*, 'aod', aod(100,100), aod(130,200), maxval(aod)
!  enddo



!  do iq = 1, x%n3d

!     idxTable = GOCART2G_MieQueryIdx(Mie,x%r3(iq)%name,rc)
!     if(idxTable == -1) cycle
!     if ( rc/=0 ) then
!        __raise__(MAPL_RC_ERROR,"cannot get Mie index for "//trim(x%r3(iq)%name))
!     end if

!     if (verbose_) &
!          print *, '[+] Adding '//trim(x%r3(iq)%name)//' contribution'

!    Loop over x, y, z
!    -----------------
!     do k = 1, km
!    delm = x%r3(iq)%q(i,j,k) * x%coords%lcv%delp(i,j,k) / MAPL_GRAV
!                 delc = x%r3(iq)%q(i,j,k) * x%r3(iRHO)%q(i,j,k) ! concentration
!                 delz = x%coords%lcv%delp(i,j,k) / (MAPL_GRAV*x%r3(iRHO)%q(i,j,k))
!                   rh = x%r3(iRH)%q(i,j,k)     
!           do j = 1, jm
!           do i = 1, im

!                 delm = x%r3(iq)%q(i,j,k) * x%coords%lcv%delp(i,j,k) / MAPL_GRAV
!                 delc = x%r3(iq)%q(i,j,k) * x%r3(iRHO)%q(i,j,k) ! concentration
!                 delz = x%coords%lcv%delp(i,j,k) / (MAPL_GRAV*x%r3(iRHO)%q(i,j,k))
!                   rh = x%r3(iRH)%q(i,j,k)

!                 call GOCART2G_MieQuery(Mie, idxTable, float(n), delm, rh, tau=tau, bbck=bck )

!                 if ( iTau>0 ) y%r3(iTau)%q(i,j,k) = y%r3(iTau)%q(i,j,k) + tau
!
!                 if ( iExt>0 ) y%r3(iExt)%q(i,j,k) = y%r3(iExt)%q(i,j,k) + tau / delz

!                 if ( iBck>0 ) y%r3(iBck)%q(i,j,k) = y%r3(iBck)%q(i,j,k) + bck * delc

!           end do ! longitudes
!        end do ! latitudes
!     end do ! levels

!  end do ! aerosol tracers


  if (verbose_) &
       print *, '[x] All done!'
        
  rc = 0

  end subroutine GOCART2G_AopCalculator3D

 end module GOCART2G_AopMod

