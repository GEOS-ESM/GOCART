
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

!
! !DESCRIPTION:
!
!  This module implements an Aerosol Optical Property calculator.
!
! !REVISION HISTORY:
!
!  29Nov2010 da Silva  Revamped to use Simple Bundle construct.
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

  type(MAPL_SimpleBundle), intent(inout) :: y          ! Aerosol extinction
  type(MAPL_SimpleBundle), intent(in)    :: x          ! Aerosol mixing ratio
  type(GOCART2G_Mie),          intent(inout) :: Mie(:) ! Mie tables for each tracer

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

     idxTable = GOCART2G_MieQueryIdx(Mie,x%r3(iq)%name,rc)
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

                 call GOCART2G_MieQuery(Mie, idxTable, float(n), delm, rh, tau=tau, bbck=bck )

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

  end subroutine GOCART2G_ExtCalculator

 end module GOCART2G_AopMod

