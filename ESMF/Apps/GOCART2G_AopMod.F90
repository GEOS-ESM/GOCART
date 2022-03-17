
#include "MAPL_Generic.h"

!-------------------------------------------------------------------------
!      NASA/GSFC, Global Modeling & Assimilation Office, Code 900.3      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  Chem_AodMod --- Aerosol Optical Depth Calculator
!
! !INTERFACE:
!

   module  Chem_AodMod

! !USES:


   use ESMF
   use MAPL

   use Chem_MieTableMod
   use Chem_MieMod
   use Chem_RegistryMod

   implicit none

!
! !PUBLIC MEMBER FUNCTIONS:
!
   PUBLIC  Chem_AodCalculator  ! Calculates AOD given mixing ratios
   PUBLIC  Chem_rEffCalculator ! Calculates Aerosol Effective Radius

!
! !DESCRIPTION:
!
!  This module implements an Aerosol Optical Depth calculator.
!
! !REVISION HISTORY:
!
!  29Nov2010 da Silva  Revamped to use Simple Bundle construct.
!  12Jan2011 da Silva  Reff capability.
!
!EOP
!-------------------------------------------------------------------------

  interface Chem_AodCalculator
     Module Procedure Chem_AodCalculator2d
  end interface Chem_AodCalculator

CONTAINS

!-----------------------------------------------------------------------------
!BOP
 
! !IROUTINE: Chem_AodCalculator ---- Mie Calculator given Aerosol Mixing Ratio
!
! !INTERFACE:
!
  subroutine Chem_AodCalculator2d ( y, x, Mie, verbose, rc ) 


! !ARGUMENTS:

  type(MAPL_SimpleBundle), intent(inout) :: y   ! Already created Bundle
  type(MAPL_SimpleBundle), intent(in)    :: x   ! Aerosol concentrations
  type(Chem_Mie),          intent(inout) :: Mie ! Mie tables, etc

  logical, OPTIONAL,       intent(in)    :: verbose
  integer, OPTIONAL,       intent(out)   :: rc

! !DESCRIPTION: Compute AOD given aerosol mixing ratio.
!
!EOP
!-----------------------------------------------------------------------------

  integer iq, idxTable, im, jm, km, iRH, i, j, k, n
  real(ESMF_KIND_R4) :: tau, delm, rh
  logical :: verbose_

              __Iam__("Chem_AodCalculator")

  if ( present(verbose) ) then
     verbose_ = verbose .and. MAPL_AM_I_ROOT()
  else
     verbose_ = .FALSE.
  end if

  im = size(x%coords%lons,1)
  jm = size(x%coords%lons,2)
  km = size(x%coords%levs)

  iRH = MAPL_SimpleBundleGetIndex(x,'RH',rank=3,rc=STATUS,quiet=.TRUE.)
  if ( STATUS /= 0 ) then
     iRH = MAPL_SimpleBundleGetIndex(x,'RH2',rank=3,__RC__)
  end if

! Consistency check
! -----------------
  if ( Mie%nch /= size(y%coords%levs) ) then
     __raise__(MAPL_RC_ERROR,"mieTable/y has inconsistent number of channels")
  end if
  if ( .not. associated(x%coords%lcv%delp) ) then
     __raise__(MAPL_RC_ERROR,"x%coords%lcv must have delp set")
  end if

! Initialize output arrays to zero
! --------------------------------
  y%r3(1)%q = 0.0

! Loop over aerosol species
! -------------------------
  do iq = 1, x%n3d

     idxTable = Chem_MieQueryIdx(Mie,x%r3(iq)%name,rc)

     if (idxTable == -1) then
         if (verbose_) &
          print *, '[-] Skipping '//trim(x%r3(iq)%name)//' contribution'

         cycle
     end if

     if ( rc/=0 ) then
        __raise__(MAPL_RC_ERROR,"cannot get Mie index for "//trim(x%r3(iq)%name))
     end if

     if (verbose_) &
          print *, '[+] Adding '//trim(x%r3(iq)%name)//' contribution'

!    Loop over channel, x, y, z
!    --------------------------
     do n = 1, Mie%nch
        do k = 1, km
           do j = 1, jm
              do i = 1, im

                 delm = x%r3(iq)%q(i,j,k) * x%coords%lcv%delp(i,j,k) / MAPL_GRAV
                   rh = x%r3(iRH)%q(i,j,k)

                 call Chem_MieQuery(Mie, idxTable, float(n), delm, rh, tau=tau)

                 y%r3(1)%q(i,j,n) = y%r3(1)%q(i,j,n) + tau

             end do ! longitudes
          end do ! latitudes
       end do ! levels
    end do ! channels

  end do ! aerosol tracers


  if (verbose_) &
       print *, '[x] All done!'
        
  rc = 0

  end subroutine Chem_AodCalculator2d

!-----------------------------------------------------------------------------
!BOP
 
! !IROUTINE: Chem_rEffCalculator ---- Compute Aerosol Effective Radius
!
! !INTERFACE:
!
  subroutine Chem_rEffCalculator ( y, x, Mie, verbose, rc ) 


! !ARGUMENTS:

  type(MAPL_SimpleBundle), intent(inout) :: y   ! Already created Bundle
  type(MAPL_SimpleBundle), intent(in)    :: x   ! Aerosol concentrations
  type(Chem_Mie),          intent(inout) :: Mie ! Mie tables, etc

  logical, OPTIONAL,       intent(in)    :: verbose
  integer, OPTIONAL,       intent(out)   :: rc

! !DESCRIPTION: Return Aerosol Effective Radius.
!
!EOP
!-----------------------------------------------------------------------------

  integer iq, idxTable, im, jm, km, iRH, i, j, k
  real(ESMF_KIND_R4) :: rEff, delm, rh
  logical :: verbose_

              __Iam__("Chem_rEffCalculator")

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

! Consistency check
! -----------------
  if ( .not. associated(x%coords%lcv%delp) ) then
     __raise__(MAPL_RC_ERROR,"x%coords%lcv must have delp set")
  end if

  delm = 0.0

! Loop over aerosol species
! -------------------------
  do iq = 1, x%n3d

     idxTable = Chem_MieQueryIdx(Mie,x%r3(iq)%name,rc)
     if(idxTable == -1) cycle
     if ( rc/=0 ) then
        __raise__(MAPL_RC_ERROR,"cannot get Mie index for "//trim(x%r3(iq)%name))
     end if

     if (verbose_) &
          print *, '[+] Adding '//trim(x%r3(iq)%name)//' contribution'

!    Loop over x, y, z
!    -----------------
     do k = 1, km
        do j = 1, jm
           do i = 1, im
              
              rh = x%r3(iRH)%q(i,j,k)

              call Chem_MieQuery(Mie, idxTable, 1.0, delm, rh, rEff=rEff)

              y%r3(iq)%q(i,j,k) = rEff


           end do ! longitudes
        end do ! latitudes
     end do ! levels

  end do ! aerosol tracers


  if (verbose_) &
       print *, '[x] All done!'
        
  rc = 0

  end subroutine Chem_rEffCalculator

!-----------------------------------------------------------------------------
!BOP
 
! !IROUTINE: Chem_ConcCalculator ---- Compute Aerosol Effective Radius
!
! !INTERFACE:
!
  subroutine Chem_ConcCalculator ( y, x, Mie, verbose, rc ) 


! !ARGUMENTS:

  type(MAPL_SimpleBundle), intent(inout) :: y   ! Already created Bundle
  type(MAPL_SimpleBundle), intent(in)    :: x   ! Aerosol concentrations
  type(Chem_Mie),          intent(inout) :: Mie ! Mie tables, etc

  logical, OPTIONAL,       intent(in)    :: verbose
  integer, OPTIONAL,       intent(out)   :: rc

! !DESCRIPTION: Return Aerosol Concentration bundle.
!
!EOP
!-----------------------------------------------------------------------------

  integer iq, im, jm, km, i, j, k, idxTable
  real(ESMF_KIND_R4) :: delm
  logical :: verbose_

              _Iam_("Chem_ConcCalculator")

  if ( present(verbose) ) then
     verbose_ = verbose .and. MAPL_AM_I_ROOT()
  else
     verbose_ = .FALSE.
  end if

  im = size(x%coords%lons,1)
  jm = size(x%coords%lons,2)
  km = size(x%coords%levs)

! Consistency check
! -----------------
  if ( .not. associated(x%coords%lcv%delp) ) then
     __raise__(MAPL_RC_ERROR,"x%coords%lcv must have delp set")
  end if

! Loop over aerosol species
! -------------------------
  do iq = 1, x%n3d

     idxTable = Chem_MieQueryIdx(Mie,x%r3(iq)%name,rc)
     if(idxTable == -1) cycle
     if ( rc/=0 ) then
        __raise__(MAPL_RC_ERROR,"cannot get Mie index for "//trim(x%r3(iq)%name))
     end if

     if (verbose_) &
          print *, '[+] Adding '//trim(x%r3(iq)%name)//' contribution'

!    Loop over x, y, z
!    -----------------
     do k = 1, km
        do j = 1, jm
           do i = 1, im
              
              delm = x%r3(iq)%q(i,j,k) * x%coords%lcv%delp(i,j,k) / MAPL_GRAV

              y%r3(iq)%q(i,j,k) = delm

           end do ! longitudes
        end do ! latitudes
     end do ! levels

  end do ! aerosol tracers


  if (verbose_) &
       print *, '[x] All done!'
        
  rc = 0

  end subroutine Chem_ConcCalculator

!-----------------------------------------------------------------------------
!BOP
 
! !IROUTINE: Chem_ExtCalculator ---- Compute 3D Aerosol Extinction
!
! !INTERFACE:
!
  subroutine Chem_ExtCalculator ( y, x, Mie, verbose, rc ) 


! !ARGUMENTS:

  type(MAPL_SimpleBundle), intent(inout) :: y   ! Aerosol extinction
  type(MAPL_SimpleBundle), intent(in)    :: x   ! Aerosol mixing ratio
  type(Chem_Mie),          intent(inout) :: Mie ! Mie tables, etc

  logical, OPTIONAL,       intent(in)    :: verbose
  integer, OPTIONAL,       intent(out)   :: rc

! !DESCRIPTION: Return Aerosol Extinction bundle.
!
!EOP
!-----------------------------------------------------------------------------

  integer :: iq, im, jm, km, i, j, k, n, idxTable
  integer :: iTau, iExt, iBck, iRH, iRHO
  real(ESMF_KIND_R4) :: delc, delm, delz, rh, tau, bck
  logical :: verbose_

              __Iam__("Chem_ExtCalculator")

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

  iTau = MAPL_SimpleBundleGetIndex(y,'tau',       rank=3,quiet=.true.,rc=STATUS)
  iExt = MAPL_SimpleBundleGetIndex(y,'extinction',rank=3,quiet=.true.,rc=STATUS)
  iBck = MAPL_SimpleBundleGetIndex(y,'backscat',  rank=3,quiet=.true.,rc=STATUS)

! Consistency check
! -----------------
  if ( .not. associated(x%coords%lcv%delp) ) then
     __raise__(MAPL_RC_ERROR,"x%coords%lcv must have delp set")
  end if
  if ( Mie%nch /= 1 ) then
     __raise__(MAPL_RC_ERROR,"for now, Mie must have a single channel")
  end if
  n = 1

! Loop over aerosol species
! -------------------------
  if ( iTau>0 ) y%r3(iTau)%q = 0.0
  if ( iExt>0 ) y%r3(iExt)%q = 0.0
  if ( iBck>0 ) y%r3(iBck)%q = 0.0

  do iq = 1, x%n3d

     idxTable = Chem_MieQueryIdx(Mie,x%r3(iq)%name,rc)
     if(idxTable == -1) cycle
     if ( rc/=0 ) then
        __raise__(MAPL_RC_ERROR,"cannot get Mie index for "//trim(x%r3(iq)%name))
     end if

     if (verbose_) &
          print *, '[+] Adding '//trim(x%r3(iq)%name)//' contribution'

!    Loop over x, y, z
!    -----------------
     do k = 1, km
        do j = 1, jm
           do i = 1, im

                 delm = x%r3(iq)%q(i,j,k) * x%coords%lcv%delp(i,j,k) / MAPL_GRAV
                 delc = x%r3(iq)%q(i,j,k) * x%r3(iRHO)%q(i,j,k) ! concentration
                 delz = x%coords%lcv%delp(i,j,k) / (MAPL_GRAV*x%r3(iRHO)%q(i,j,k))
                   rh = x%r3(iRH)%q(i,j,k)

                 call Chem_MieQuery(Mie, idxTable, float(n), delm, rh, tau=tau, bbck=bck )

                 if ( iTau>0 ) y%r3(iTau)%q(i,j,k) = y%r3(iTau)%q(i,j,k) + tau

                 if ( iExt>0 ) y%r3(iExt)%q(i,j,k) = y%r3(iExt)%q(i,j,k) + tau / delz

                 if ( iBck>0 ) y%r3(iBck)%q(i,j,k) = y%r3(iBck)%q(i,j,k) + bck * delc

           end do ! longitudes
        end do ! latitudes
     end do ! levels

  end do ! aerosol tracers


  if (verbose_) &
       print *, '[x] All done!'
        
  rc = 0

  end subroutine Chem_ExtCalculator

 end module Chem_AodMod

