#include "unused_dummy.H"
! Module: ConvectionMod --- mimic the convective scavenging algorithm from the 
!                                 offline GOCART CTM

  module ConvectionMod

! USES
  use Chem_Mod
  use m_die

! PUBLIC

! Description
! TC is tracer mixing ratio (Mass Mixing Ratio)
! Note that I change the sense of this so that rather than calling all aerosols
! at once (as in offline CTM) I call separately from each component (DU, SS, etc.)
! Assumption is that H2O2 passed in is in units of mass mixing ratio

! Parameters
 real*8, parameter :: kc = 5.0e-3      ! conversion rate of cloud condensate to precipation [s-1]
 logical, parameter :: lsadirect = .false.

! ALT: for the time being we need a simple mechanism to disable
!      convection if GF was chosen in Moist 
!      (otherwise we would be doing it twice) 
 logical, private :: doing_convection = .true.
 public enable_convection
 public disable_convection
 
CONTAINS
  subroutine enable_convection
    doing_convection = .true.
  end subroutine enable_convection

  subroutine disable_convection
    doing_convection = .false.
  end subroutine disable_convection

  subroutine zflip (varin, varout, km)
! reorder a variable in the vertical
  real*8, dimension(:,:,:) :: varin, varout
  integer*4                :: km, k

  do k = 1, km
   varout(:,:,k) = varin(:,:,km-k+1)
  enddo

end subroutine zflip


! ----------------------------------------------------------------------------------
! OMIT sulfate stuff for now
  SUBROUTINE convection(i1, i2, j1, j2, km, n1, n2, dt_30m, aero_type, kin, &
                        tc, cldmas, dtrain, area, delz, delp, vud, &
                        airmass, airmol, tmpu, ple, bcnv, h2o2 )

  IMPLICIT NONE

  INTEGER, INTENT(IN)    :: i1, i2, j1, j2, km, n1, n2, dt_30m
  character(len=*)     :: aero_type
  REAL*8,    INTENT(INOUT) :: tc(i1:i2,j1:j2,km,n1:n2)
!  REAL*8,    INTENT(INOUT) :: cldscv(i1:i2,j1:j2,km,n1:n2), cldso2(i1:i2,j1:j2,km)
!  REAL*8,    INTENT(INOUT) :: cldso4(i1:i2,j1:j2,km), cldmsa(i1:i2,j1:j2,km)
  REAL*8                    :: cldscv(i1:i2,j1:j2,km,n1:n2), cldso2(i1:i2,j1:j2,km)
  REAL*8                    :: cldso4(i1:i2,j1:j2,km), cldmsa(i1:i2,j1:j2,km)
!  REAL*8,    INTENT(INOUT) :: tcnv(i1:i2,j1:j2,n1:n2)
!  REAL*8,    INTENT(INOUT) :: wet_conv_in(i1:i2,j1:j2,km,n1:n2)
  REAL*8,    INTENT(INOUT) :: airmass(i1:i2,j1:j2,km)
  REAL*8,    DIMENSION(i1:i2,j1:j2,km),   INTENT(IN) :: vud
  REAL*8,    DIMENSION(i1:i2,j1:j2,km+1), INTENT(IN) :: cldmas
  REAL*8,    DIMENSION(i1:i2,j1:j2,km),   INTENT(IN) :: dtrain
  REAL*8,    DIMENSION(i1:i2,j1:j2),      INTENT(IN) :: area
  REAL*8,    DIMENSION(i1:i2,j1:j2,km),   INTENT(IN) :: delz, delp
  REAL*8,    DIMENSION(i1:i2,j1:j2,km),   INTENT(IN) :: airmol
  REAL*8,    DIMENSION(i1:i2,j1:j2,km),   INTENT(IN) :: tmpu
  REAL*8,    DIMENSION(i1:i2,j1:j2,km+1), INTENT(IN) :: ple
  REAL*8,    INTENT(OUT)   :: bcnv(i1:i2,j1:j2,n1:n2)
  REAL*8,    INTENT(INOUT), optional :: h2o2(i1:i2,j1:j2,km)
  LOGICAL, INTENT(INOUT)   :: KIN  ! true for aerosol

  REAL*8    :: tc1(i1:i2,j1:j2,km,n1:n2), f(i1:i2,j1:j2,km,n1:n2)
  REAL*8    :: cldmas_tmp(i1:i2,j1:j2,km), so2loss
! epsilon: A very small positive number   [unitless]
  REAL*8,  PARAMETER   :: EPSILON = 1.0E-32
  REAL*8,  PARAMETER   :: R = 8.2057d-2  ! universal gas constant [L*atm/moles/K]
  REAL*8,  PARAMETER   :: INV_T0 = 1d0 / 298d0
  REAL*8,  PARAMETER   :: conv_NH3 = 5.69209978831d-1 ! 0.6*SQRT(0.9) for ice to gas ratio
  REAL*8  :: kg, Kstar298, H298_R, I2G, L2G, C_TOT, F_L, F_I
  INTEGER :: n, i, j, l, NSO2, NSO4, NMSA
  REAL*8,    DIMENSION(i1:i2,j1:j2,km) :: c_h2o
  REAL*8,    DIMENSION(i1:i2,j1:j2,km) :: cldliq
  REAL*8,    DIMENSION(i1:i2,j1:j2,km) :: cldice


!srf----------------------
!ALT: a better protection to not do convection is Moist chooses GF
  if (.not. doing_convection) then
     bcnv(:,:,:) = 0.0
     RETURN
  end if
!srf----------------------


!  Initialize local variables
!  --------------------------
!  c_h2o, cldliq, and cldice are respectively intended to be the 
!  water mixing ratio (liquid or vapor?, in or out of cloud?)
!  cloud liquid water mixing ratio
!  cloud ice water mixing ratio
   c_h2o  = (10d0**(-2663.5d0/tmpu(:,:,:) + 12.537d0 ) ) /  &
                   (ple(:,:,1:km)+ple(:,:,2:km+1)) /2d0   
   cldliq = 0.d0
   where(tmpu >  248.) cldliq = 1.d-6 * ( ( tmpu - 248.d0) / 20.d0 )
   where(tmpu >= 268.) cldliq = 1.d-6
   cldice = 1.d-6 - cldliq

  ! executable statements

  tc1(:,:,:,:) = tc(:,:,:,:)
  !if (MAPL_AM_I_ROOT()) print *, 'hbian convection tmpu =', tmpu(i1,j1,1), tmpu(i1,j1,km)
  
! compute the fraction of tracer scavenged in convective cloud updrafts
  f = 0.0
  kg = 0d0

  DO n = n1, n2
     if (TRIM(aero_type) .eq. 'nitrate' .and. n .eq. n1 )  kin = .false.  ! treat NH3 as a gas tracer
     if (TRIM(aero_type) .eq. 'nitrate' .and. n .gt. n1 )  kin = .true.   ! treat others as aerosol

     if (kin) then
        CALL f_aerosol(i1, i2, j1, j2, km, kc, f(:,:,:,n), delz, vud )
     else
        ! gas tracer NH3
        if (TRIM(aero_type) .eq. 'nitrate' .and. n .eq. n1 )  then
        ! values adopted in Umich/IMPACT and GMI, effective Herry's law coefficient at pH = 5
           Kstar298 = 1.05d6
           H298_R = -4.2d3
        endif
           DO L = 2, KM
           DO J = j1, j2
           DO I = i1, i2
              ! ice to gas ratio
              if ( c_h2o(i,j,l) > 0.d0) then
                 I2G = (cldice(i,j,l) / c_h2o(i,j,l)) * conv_NH3
              else
                 I2G = 0.d0
              endif
              L2G = cldliq(i,j,l) * R * tmpu(i,j,l) * &
                      Kstar298 * EXP( -H298_R * ( ( 1d0 / tmpu(i,j,l) ) - INV_T0 ) )
              ! fraction of NH3 in liquid & ice phases
              C_TOT = 1d0 + L2G + I2G
              F_L = L2G / C_TOT
              F_I = I2G / C_TOT
              ! compute kg, the retention factor for liquid NH3 is 0 at T < 248K and
              ! 0.05 at 248K < T < 268K
              if (tmpu(i,j,l) >=268d0) then
                 kg = kc * ( F_L+F_I )
              elseif ( (248d0 < tmpu(i,j,l)) .and. (tmpu(i,j,l) < 268d0) ) then
                 kg = kc * ( (0.05*F_L)+F_I )
              else
                 kg = kc * F_I
              endif
              if(kg > 0.d0 .and. vud(i,j,l) > 1.e-14) &
               f(i,j,l,n) = 1.0 - EXP( -kg * delz(i,j,l) / vud(i,j,l) )
           ENDDO
           ENDDO
           ENDDO
     endif

!    Special treatment for DMS and SO2 if aero_type is "sulfur"
     if(trim(aero_type) .eq. 'sulfur') then

        if(.not.present(h2o2)) call die ('GOCARTConvectionMod.F90', &
                                         'missing required H2O2 for sulfur') 

        if(n .eq. n1)    f(:,:,:,n1) = 0.0    ! DMS
        !if(n .eq. n1+1)  f(:,:,:,n1+1) = 0.0  ! SO2 for now is not scavenged   
#undef PRC
!#ifdef PRC
        if(n .eq. n1+1) then                  ! SO2 requires special handling

        !==============================================================
        ! Coupled full chemistry/aerosol simulation:
        ! Use the wet scavenging formula of Chin et al [1996], 
        ! such that a soluble fraction of SO2 is limited by the
        ! availability of H2O2 in the precipitating grid box. 
        ! Scavenge the soluble SO2 at the same rate as the sulfate.
        ! Update H2O2_sav and SO2_sav for use in RAINOUT, WASHOUT
        !==============================================================
        DO L = 2, KM
           DO J = j1, j2
              DO I = i1, i2

                 ! Make sure to deplete H2O2s the same as SO2s.
                 ! (dkh, rjp, bmy, 11/17/05)
                 ! based on GEOS-Chem. tq, 01/09
                 
                 IF ( tc1(i,j,l,n) > epsilon ) THEN
                    
                    ! limit f
                    so2loss  = MIN( h2o2(i,j,l), tc1(i,j,l,n) )
                    f(i,j,l,n) = f(i,j,l,n) * so2loss / tc1(i,j,l,n)
                    f(i,j,l,n) = MAX(f(i,j,l,n), 0.0)
                    
                    ! update saved h2o2 concentration
                    h2o2(i,j,l) = h2o2(i,j,l) - ( tc1(i,j,l,n) * f(i,j,l,n) )
                    h2o2(i,j,l) = MAX( h2o2(i,j,l), epsilon )
                    
                 ELSE
                    
                    ! set f = 0 if so2 < epsilon (dkh, rjp, bmy, 11/17/05)
                    f(i,j,l,n) = 0.d0
                    
                 END IF
                 
              ENDDO
           ENDDO
        ENDDO
        endif                            ! SO2
!#endif
     endif                               ! sulfur

  ENDDO  ! n

! if tracer is type "carbon" then set coefficient to 0 for hydrophobic
! implementing QQ Wang's change by Huisheng Bian (4/24/2015)
! not scavenging  BCn1 (hydrophobic) when T > 258 K
  if(trim(aero_type) .eq. 'OC') f(:,:,:,n1) = 0d0

! suppress scavenging most aerosols at cold T except BCn1 (hydrophobic), dust, and HNO3
  if (trim(aero_type) .eq. 'BC') then
     where (tmpu >= 258.d0)
        f(:,:,:,n1) = 0.d0
     end where
  end if

  if (trim(aero_type) .eq. 'BC') then
     where (tmpu < 258.d0)
        f(:,:,:,n2) = 0.d0
     end where
  end if

  if (trim(aero_type) .eq. 'OC'       .or. &
      trim(aero_type) .eq. 'sea_salt' .or. &
      trim(aero_type) .eq. 'sulfur'   .or. &
      trim(aero_type) .eq. 'seasalt'  .or. &
      trim(aero_type) .eq. 'sulfate'  .or. &
      trim(aero_type) .eq. 'nitrate'  .or. &
      trim(aero_type) .eq. 'bromine'  .or. &
      trim(aero_type) .eq. 'NH3'      .or. &
      trim(aero_type) .eq. 'NH4a') then

      do n = n1, n2
         where (tmpu < 258.d0 )
            f(:,:,:,n) = 0.d0
         endwhere
      end do

  end if


  ! re-index for routine cldcnv
  cldmas_tmp(:,:,1:km) = cldmas(:,:,2:km+1)

  ! internal time step for the convection routine is 300s
  CALL cldcnv(i1, i2, j1, j2, km, n1, n2, dt_30m, aero_type, &
              tc, f, airmass, area, cldmas_tmp, dtrain, delz, delp)

  ! -- Mass balance
  SELECT CASE(TRIM(aero_type))

  CASE('sulfur')
     NSO2 = n1+1
     NSO4 = n1+2
     NMSA = n1+3

     cldso2 = 0.0d0
     cldso4 = 0.0d0
     cldmsa = 0.0d0

     DO l = 1,km
        DO j = j1, j2
           DO i = i1, i2
              cldso2(i,j,l) = cldso2(i,j,l) + &
                  (tc(i,j,l,NSO2) - tc1(i,j,l,NSO2)) * airmass(i,j,l)
              cldso4(i,j,l) = cldso4(i,j,l) + &
                  (tc(i,j,l,NSO4) - tc1(i,j,l,NSO4)) * airmass(i,j,l)
              cldmsa(i,j,l) = cldmsa(i,j,l) + &
                  (tc(i,j,l,NMSA) - tc1(i,j,l,NMSA)) * airmass(i,j,l)
           END DO
        END DO
     END DO
     
     DO n = n1, n2
        DO i = i1, i2
           DO j = j1, j2
              bcnv(i,j,n) = 0.0
              DO l = 1,km
                 IF (tc(i,j,l,n) < 0.0) tc(i,j,l,n) = 1.0E-32
                 ! kg tracer
                 bcnv(i,j,n) = bcnv(i,j,n) + &
                     (tc(i,j,l,n) - tc1(i,j,l,n)) *airmass(i,j,l)
              END DO
           END DO
        END DO
     END DO

  CASE('co')
     cldscv = 0.0d0
     DO n = n1, n2
        DO l = 1,km
           DO j = j1, j2
              DO i = i1, i2
                 cldscv(i,j,l,n) = cldscv(i,j,l,n) + &
 !                wet_conv_in(i,j,l,n) = wet_conv_in(i,j,l,n) + &
                      (tc(i,j,l,n) - tc1(i,j,l,n)) * airmol(i,j,l)
              END DO
           END DO
        END DO
     END DO

     bcnv(:,:,:) = 0.0
     DO n = n1, n2
        DO i = i1, i2
           DO j = j1, j2
              DO l = 1,km
                 IF (tc(i,j,l,n) < 0.0) tc(i,j,l,n) = 1.0E-32
                 bcnv(i,j,n) = bcnv(i,j,n) + &
                      (tc(i,j,l,n) - tc1(i,j,l,n)) * airmol(i,j,l)
              END DO
           END DO
        END DO
     END DO

  CASE DEFAULT

!!$     DO n = n1, n2
!!$        DO l = 1,km
!!$           DO j = j1, j2
!!$              DO i = i1, i2
!!$                 wet_conv_in(i,j,l,n) = wet_conv_in(i,j,l,n) + &
!!$                      (tc(i,j,l,n) - tc1(i,j,l,n)) * airmass(i,j,l)
!!$              END DO
!!$           END DO
!!$        END DO
!!$     END DO

     DO n = n1, n2
        DO i = i1, i2
           DO j = j1, j2
              bcnv(i,j,n) = 0.0
              DO l = 1,km
                 IF (tc(i,j,l,n) < 0.0) tc(i,j,l,n) = 1.0E-32
                 bcnv(i,j,n) = bcnv(i,j,n) + &
                      (tc(i,j,l,n) - tc1(i,j,l,n)) * airmass(i,j,l)
              END DO
           END DO
        END DO
     END DO

  END SELECT

!  tcnv(:,:,:) = tcnv(:,:,:) + bcnv(:,:,:)

END SUBROUTINE convection



! ----------------------------------------------------------------------------------
! set_vud
  SUBROUTINE set_vud(i1, i2, j1, j2, km, frlake, frocean, frseaice, cldmas, qccu, &
                     airmass, delz, area, vud)

    INTEGER, INTENT(IN)  :: i1, i2, j1, j2, km
    REAL*8,    INTENT(IN), DIMENSION(i1:i2,j1:j2)       :: frlake, frocean, frseaice
    REAL*8,    INTENT(IN), DIMENSION(i1:i2,j1:j2)       :: area
    REAL*8,    INTENT(IN), DIMENSION(i1:i2,j1:j2,km)    :: qccu, airmass, delz
    REAL*8,    INTENT(IN), DIMENSION(i1:i2,j1:j2,km+1)  :: cldmas
    REAL*8,    INTENT(OUT) :: vud(i1:i2,j1:j2,km)

    REAL*8, PARAMETER :: max_vud=100 ! maximum updraft velocity [m/s]
    REAL*8 :: water(i1:i2,j1:j2), dvud(i1:i2,j1:j2,km)
    INTEGER :: i, j, k

    ! executable statements

    !==============================================================
    ! Compute vud -- 5 m/s over oceans, 10 m/s over land
    ! Assume vud is the same at all altitudes; the array can be 2-D
    !==============================================================
!    WHERE ((frlake + frocean - frseaice) >= 0.5) ! water 
!       dvud = 5.0
!    ELSEWHERE ! land (including permanent ice) and sea ice
!       dvud = 10.0
!    END WHERE

    water = frlake + frocean - frseaice
    WHERE (water < 0.0) water = 0.0
    WHERE (water > 1.0) water = 1.0
    ! Compute updraft velocity as a weighted average over water and land+ice. 
    DO k = 1,km
!    dvud(:,:,k) = water(:,:)*5.0 + MAX((1.0-water(:,:)),0.0)*10.0
       dvud(:,:,k) = water(:,:)*5.0 + (1.0-water(:,:))*10.0
    END DO

    ! compute updraft velocity from cldmas=rho*vud

    DO k = 1,km-1
       DO j = j1,j2
          DO i = i1,i2
             IF (qccu(i,j,k) >= TINY(0.0)) THEN
                vud(i,j,k) = cldmas(i,j,k+1)/(qccu(i,j,k)*airmass(i,j,k)) * &
                             area(i,j) * delz(i,j,k)
             ELSE
               vud(i,j,k) = dvud(i,j,k)
            END IF
         END DO
      END DO
   END DO

   vud(:,:,km) = 0.0

    ! What should be used as threshold value here? 100 m/s?
!   WHERE (vud > max_vud) vud = max_vud
   WHERE (vud > max_vud) vud = dvud 

  END SUBROUTINE set_vud


! ----------------------------------------------------------------------------------
!  SUBROUTINE COMPUTE_F( i1, i2 ,j1 ,j2, km, n, aero_type, F, bxheight, vud, tc1, h2o2, kc)
  SUBROUTINE COMPUTE_F( i1, i2 ,j1 ,j2, km, n, aero_type, F, bxheight, vud, tc1, kc)
!
!******************************************************************************
!  Subroutine COMPUTE_F computes F, the fraction of soluble tracer lost by 
!  scavenging in convective cloud updrafts. (hyl, bmy, djj, 2/23/00, 7/20/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) N    (INTEGER) : Tracer number
!
!  Arguments as Output:
!  ============================================================================
!  (2 ) F    (REAL*8)  : Fraction of tracer scavenged in cloud updraft [0-1]
!
!  References (see above for full citations):
!  ===========================================================================
!  (1 ) Jacob et al, 2000
!  (2 ) Chin et al, 1996
!
!  NOTES:
!  (1 ) Currently works computes scavenging fractions for either full
!        chemistry simulation (NSRCX == 3) or Rn-Pb-Be chemistry simulation
!        (NSRCX == 1).  Set the scavenging fraction to zero for other
!        simulations which do not carry soluble tracers. (bmy, 3/2/00)
!  (2 ) Need to call INIT_SCAV to initialize the Vud, C_H2O, CLDLIQ, 
!        and CLDICE fields once per timestep. (bmy, 2/23/00)
!  (3 ) For aerosols only: now apply Eq. 2 for all temperatures.  Also
!        use the distance between the grid box centers in Eq. 2.  Updated
!        comments and made some cosmetic changes (hyl, bmy, 6/18/01)
!  (4 ) Remove IREF, JREF -- these are obsolete.  T is now dimensioned
!        (I1:I2,J1:J2,KM).  T(IREF,JREF,L) is now T(I,J,L). (bmy, 9/27/01)
!  (5 ) Removed obsolete code from 9/01 (bmy, 10/23/01)
!  (6 ) Fix 2 bugs for aerosol scavenging in Rn-Pb-Be simulation: 
!        (a) set F(:,:,1) = 0 since we don't do any scavenging there.  
!        (b) DO L = 2, KM to avoid any subscript range out of bounds 
!        errors (rjp, hyl, bmy, 1/10/02)
!  (7 ) Now set F=0 in the first level for all tracers.  Also now
!        compute the distance between grid box centers and use that in 
!        in Eq. 10 from Jacob et al, 2000 to compute F. (hyl, bmy, 1/24/02)
!  (8 ) Eliminated obsolete code from 1/02 (bmy, 2/27/02)
!  (9 ) Now reference T from "dao_mod.f" instead of from "CMN".  Also reference
!        BXHEIGHT from "dao_mod.f" instead of from "CMN_NOX".  Now bundled
!        into "wetscav_mod.f".  Now references IDTHNO3, IDTH2O2, etc, from
!        F90 module "tracerid_mod.f".  Added internal routines F_AEROSOL
!        and GET_ISOL.  Rewritten so that we don't duplicate code for 
!        different chemistry simulations. (bmy, 1/17/03)
!  (10) Now compute F for SO2 in the same way for both fullchem and offline 
!        simulations (rjp, bmy, 3/23/03)
!  (11) Added slots for carbon aerosol & dust tracers.  Now modified internal
!        routine GET_ISOL so it's not hardwired anymore. (rjp, bmy, 4/5/04)
!  (12) Added slots for sea salt aerosol tracers (rjp, bec, bmy, 4/20/04)
!  (13) Added slots for secondary organic aerosol tracers (rjp, bmy, 7/13/04)
!  (14) Remove reference to CMN, it's not needed.  Made internal routine
!        F_AEROSOL a module procedure rather than an internal routine to
!        COMPUTE_F in order to facilitate parallelization on the Altix.  Also
!        now pass all arguments explicitly to F_AEROSOL. (bmy, 7/20/04)
!******************************************************************************
!
      ! References to F90 modules

!  USE mo_control, ONLY: lsadirect
!  USE mo_tracer,  ONLY: NBC1, NOC1, NBC2, NOC2, NDMS, NSO2, NSO4, NMSA, &
!                        NCO, NCOSA, NCOEA, NCOEU, NCONA, NCOOR, &
!                        NAVOC, NBVOC, NCOAVOC, NCOBVOC, NCOCH4, &
!                        NCOSM, NCOAF, NCOAU, NCOFF, NCOBF, NCOBB, NCOBI, &
!                        aero_type

  IMPLICIT NONE

  ! Arguments
  INTEGER, INTENT(IN)    :: i1, i2, j1, j2, km
  character(len=*)     :: aero_type
  INTEGER, INTENT(IN)    :: n
  REAL*8,    INTENT(IN)    :: bxheight(I1:I2,J1:J2,KM), vud(i1:i2,j1:j2,km)
  REAL*8,    INTENT(IN)    :: kc
!  REAL*8,    INTENT(INOUT) :: tc1(i1:i2,j1:j2,km), h2o2(i1:i2,j1:j2,km)
  REAL*8,    INTENT(INOUT) :: tc1(i1:i2,j1:j2,km)
  REAL*8,    INTENT(OUT)   :: f(I1:I2,J1:J2,KM)

  ! Local variables 
  
  ! Kc is the conversion rate from cloud condensate to precip [s^-1]
!  REAL*8, PARAMETER    :: KC   = 5.0E-3
  
  ! CONV = 0.6 * SQRT( 1.9 ), used for the ice to gas ratio for H2O2
  REAL*8, PARAMETER    :: CONV = 8.27042925126E-1
  !  epsilon: A very small positive number   [unitless]
  REAL*8,  PARAMETER   :: EPSILON = 1.0E-32

  !=================================================================
  ! COMPUTE_F begins here!
  !
  ! For aerosol tracers, compute F with internal routine F_AEROSOL.
  !
  ! ISOL = tracer index for the ND38 diagnostic.  Values are:
  !
  ! Tracer   Rn-Pb-Be run   Fullchem run   Offline sulfate run
  ! ------   ------------   ------------   -------------------
  !  210Pb         1              -                -
  !  7Be           2              -                -
  !  HNO3          -              1                -
  !  H2O2          -              2                7 
  !  CH2O          -              3                -
  !  MP            -              4                -
  !  SO2           -              5                1
  !  SO4           -              6                2
  !  MSA           -              7                3
  !  NH3           -              8                4
  !  NH4           -              9                5
  !  NIT           -             10                6
  !=================================================================
  
  _UNUSED_DUMMY(tc1)
  _UNUSED_DUMMY(n)

  SELECT CASE (TRIM(aero_type))

  CASE ('sulfur')

#undef PRC
#ifdef PRC

!    SO2
     IF ( n == NSO2 ) THEN

        !---------------------------
        ! SO2 (aerosol)
        !---------------------------
        
        ! Compute fraction of SO2 scavenged
        
        CALL f_aerosol(i1, i2, j1, j2, km, kc, f, bxheight, vud )
        
        !==============================================================
        ! Coupled full chemistry/aerosol simulation:
        ! Use the wet scavenging formula of Chin et al [1996], 
        ! such that a soluble fraction of SO2 is limited by the
        ! availability of H2O2 in the precipitating grid box. 
        ! Scavenge the soluble SO2 at the same rate as the sulfate.
        ! Update H2O2_sav and SO2_sav for use in RAINOUT, WASHOUT
        !==============================================================
        DO L = 2, KM
           DO J = j1, j2
              DO I = i1, i2

                 ! Make sure to deplete H2O2s the same as SO2s.
                 ! (dkh, rjp, bmy, 11/17/05)
                 ! based on GEOS-Chem. tq, 01/09
                 
                 IF ( tc1(i,j,l) > epsilon ) THEN
                    
                    ! limit f
                    so2loss  = MIN( h2o2(i,j,l), tc1(i,j,l) )
                    f(i,j,l) = f(i,j,l) * so2loss / tc1(i,j,l)
                    f(i,j,l) = MAX(f(i,j,l), 0.0)
                    
                    ! update saved h2o2 concentration
                    h2o2(i,j,l) = h2o2(i,j,l) - ( tc1(i,j,l) * f(i,j,l) )
                    h2o2(i,j,l) = MAX( h2o2(i,j,l), epsilon )
                    
                 ELSE
                    
                    ! set f = 0 if so2 < epsilon (dkh, rjp, bmy, 11/17/05)
                    f(i,j,l) = 0.0
                    
                 END IF
                 
              ENDDO
           ENDDO
        ENDDO

     ELSE IF ( n == NSO4 ) THEN
     
        !----------------------------
        ! SO4 (aerosol)
        !----------------------------

        CALL f_aerosol(i1, i2, j1, j2, km, kc, f, bxheight, vud ) 

     ELSE IF ( n == NMSA ) THEN
        
        !---------------------------
        ! MSA (aerosol)
        !---------------------------

        CALL f_aerosol(i1, i2, j1, j2, km, kc, f, bxheight, vud )

     ELSE  IF ( n == NDMS) THEN 
        
        !---------------------------
        ! DMS (aerosol)
        !---------------------------

        !----------------------------
        ! Insoluble tracer, set F=0
        !----------------------------

        F(:,:,:) = 0.0
        
     ENDIF
#endif
F = 0.0
     
  CASE ('carbon')

! PRC rewrite this a bit
!     IF ( n == n2 ) THEN

        !----------------------------
        ! HYDROPHILIC (aerosol)
        !----------------------------
     
        CALL f_aerosol(i1, i2, j1, j2, km, kc, f, bxheight, vud )
        
        
!     ELSE IF ( n == n1 ) THEN
        
!        !----------------------------
!        ! HYDROPHOBIC (aerosol)
!        !----------------------------

!        ! Force not to be lost in convective updraft for now
!        F    = 0.0

!     END IF

  CASE ('dust')

     !----------------------------
     ! DUST (aerosol) (all dust bins)
     !----------------------------
     
     CALL f_aerosol(i1, i2, j1, j2, km, kc, f, bxheight, vud )
     
  CASE ('sea_salt')
     
     !----------------------------
     ! seasalt aerosol (accum mode and coarse mode)
     !----------------------------
     
     CALL f_aerosol(i1, i2, j1, j2, km, kc, f, bxheight, vud )

  CASE ('co')

     IF ( lsadirect ) THEN
     
        CALL f_aerosol(i1, i2, j1, j2, km, kc, f, bxheight, vud ) 

     ELSE

        !----------------------------
        ! ALL CO ARE HYDROPHOBIC 
        !----------------------------
        
        ! Force not to be lost in convective updraft for now
        F    = 0.0
        
     END IF

     
  END SELECT
  
  ! Return to calling program

END SUBROUTINE COMPUTE_F


! ----------------------------------------------------------------------------------
  SUBROUTINE f_aerosol( i1, i2, j1, j2, km, kc, f, bxheight, vud) 
!
!******************************************************************************
!  Subroutine F_AEROSOL returns the fraction of aerosol scavenged in updrafts
!  (bmy, 11/7/02, 7/20/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) KC (REAL*8) : Conversion rate from cloud condensate to precip [s^-1]
!
!  Arguments as Output:
!  ============================================================================
!  (2 ) F  (REAL*8) : Fraction of aerosol scavenged in updrafts [unitless]
!
!  NOTES:
!  (1 ) Split off 
!******************************************************************************
!
  ! References to F90 modules
      
  IMPLICIT NONE

  ! Arguments
  INTEGER, INTENT(IN)  :: i1, i2, j1, j2, km
  REAL*8,    INTENT(IN)  :: kc
  REAL*8,    INTENT(IN)  :: bxheight(i1:i2,j1:j2,km), vud(i1:i2,j1:j2,km)
  REAL*8,    INTENT(OUT) :: f(i1:i2,j1:j2,km)
  
  ! Local variables
  INTEGER             :: i, j, l
  
  !=================================================================
  ! F_AEROSOL begins here!
  !
  ! Aerosol tracers are 100% in the cloud condensate phase, so 
  ! we set K = Kc, and compute F accordingly (cf Jacob et al 2000 )    
  !=================================================================
  
  ! Turn off scavenging in the first level by setting F = 0
!!!  f(:,:,1) = 0.0
  
  ! Apply scavenging in levels 2 and higher
!!  DO l = 2, km

  DO l = 1, km-1
     DO j = j1, j2
        DO i = i1, i2
           ! Distance between grid box centers [m]
!           tmp = 0.5 * ( bxheight(i,j,l-1) + bxheight(i,j,l) ) 
       
           ! (Eq. 2, Jacob et al, 2000, with K = Kc)
!           f(i,j,l) = 1.0 - EXP( -kc * tmp / vud(i,j,l) )
           if(vud(i,j,l) > 1.e-14) &
            f(i,j,l) = 1.0 - EXP( -kc * bxheight(i,j,l) / vud(i,j,l) )
           
        END DO
     END DO
  END DO
  
  ! Return to calling program
END SUBROUTINE f_aerosol

! ----------------------------------------------------------------------------------
  SUBROUTINE cldcnv(i1, i2, j1, j2, km, n1, n2, dt_conv, aero_type, &
                    q, f, airmass, area, cldmas, dtrn, delz, delp)

! ============================================================================
!                                                                           
!   This is the cumulus transport module for 3D GEOS-CTM.                   
!   Author: Shian-Jiann Lin, Code 910.3, NASA/GSFC,                         
!   Feb 12, 1997.                                                           
!   Version 3, Detrainment and Entrainment are considered.                  
!   The algorithm reduces to that of version 2 if Dtrn = 0.                 
!                                                                           
!   Q:      tracer mixing ratio.                                            
!   CLDMAS: cloud mass flux in kg/(s m**2); must be positive definite.      
!   DTRN:   Detrainment rate in kg/(s m**2); must be positive definite.     
!   NDT:    Large scale advection time step (SEC.).                         
!                                                                           
!                                                                           
!                                                                           
!                ^                                                          
!                |                                                          
!                |                                                          
!      ------ cldmas(k) -----------                                         
!                                                                           
!          q(k), dtrn(k)   (layer k)                                        
!                                                                           
!      ----- cldmas(k-1) -----------                                        
!                ^                                                          
!                |                                                          
!                |                                                          
!                                                                           
! ============================================================================
!                                                                           
!   cldcnv modified by Bob Yantosca and Mian Chin for wet scavenging of     
!   soluble species in cloud updraft.                                       
!   Modified by Thomas Diehl for usage with GEOS-5 data.                    
!                                                                           
! ============================================================================

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: i1, i2, j1, j2, km, n1, n2, dt_conv
  character(len=*)      :: aero_type
  REAL*8, INTENT(IN), DIMENSION(i1:i2,j1:j2,km) :: airmass, cldmas, dtrn
  REAL*8, INTENT(IN)    :: delz(i1:i2,j1:j2,km), delp(i1:i2,j1:j2,km), area(i1:i2,j1:j2) 
  REAL*8, INTENT(IN)    :: f(i1:i2,j1:j2,km,n1:n2)
!  REAL*8, INTENT(INOUT) :: q(i1:i2,j1:j2,km,n1:n2), h2o2(i1:i2,j1:j2,km)
  REAL*8, INTENT(INOUT) :: q(i1:i2,j1:j2,km,n1:n2)

  REAL*8    :: bmass(i1:i2,j1:j2,km), qb(i1:i2,j1:j2), mb(i1:i2,j1:j2), qc(i1:i2,j1:j2)
  REAL*8    :: tdt, qc_pres, cmout, entrn, delq
  REAL*8    :: term_1, term_2, term_3, term_4, tsum
  INTEGER :: nsteps, ktop, ic, istep, i, j, k
  REAL*8, PARAMETER :: tiny = 1.0E-14 

  ! executable statements

! ============================================================================
!   Define active convective region,from J = JS(outh) to J = JN(orth),      
!   and to level K = KTOP.                                                  
!                                                                           
!   Polar regions are too cold to have moist convection.                    
!   (Dry convection should be done elsewhere.)                              
! ============================================================================

  ktop = km - 1
!  jump = (jmx - 1)/20
!  js = 1 + jump
!  jn = jmx - js + 1
  
!write(*,*) js, jn
!write(*,*) minval(q), maxval(q)
!write(*,*) minval(delp), maxval(delp)

! ============================================================================
!   Internal time step for convective mixing is 300 sec.                    
! ============================================================================

! Use fixed internal convection time step of 300s
! nstep: number of internal time steps to reach dt_conv, which is the external
  ! time interval for the convection process

  _UNUSED_DUMMY(aero_type)
  _UNUSED_DUMMY(delz)

  nsteps = dt_conv / 300
  nsteps = MAX(nsteps,1)
  tdt    = REAL(dt_conv) / REAL(nsteps)

! ============================================================================
!   Compute air mass of each grid box (I,J,K).                              
! ============================================================================
  
  DO k = 1,km
     DO j = j1, j2
        DO i = i1, i2
           bmass(i,j,k) = airmass(i,j,k)/area(i,j)
        END DO
     END DO
  END DO

  tracer_loop: DO ic = n1, n2
     time_loop:  DO istep = 1,nsteps             

! ============================================================================
!   (1) Below cloud base                                                    
!                                                                           
!   If Cloud Mass Flux exists at (I,J,2), then compute QB.                  
!   QB is "weighted average" mixing ratio below the cloud base.             
!   QB is used to compute QC, which is defined as:                          
!                                                                           
!   QC =  ( Total mass of tracer below cloud base +                         
!           Subsidence into cloud base from above )                         
!        -------------------------------------------------------            
!                 Total air mass below cloud base                           
!                                                                           
!   MB is the total mass of air below the cloud base.                       
! ============================================================================

        j_loop_1: DO j = j1, j2
           i_loop_1: DO i = i1, i2

              IF (cldmas(i,j,2) > tiny) THEN
                 
                 qb(i,j) = &
                      (q(i,j,1,ic)*delp(i,j,1) + q(i,j,2,ic)*delp(i,j,2)) / &
                      (delp(i,j,1) + delp(i,j,2)) 

                 ! alternative:
                 ! use delz as weight 
!                 qb(i,j) = (q(i,j,1,ic)*delz(i,j,1) + &
!                            q(i,j,2,ic)*delz(i,j,2)) / &
!                            (delz(i,j,1) + delz(i,j,2))
 
                 mb(i,j) =  bmass(i,j,1) + bmass(i,j,2)

                 qc(i,j) = &
                      ( mb(i,j)       * qb(i,j) + &
                        cldmas(i,j,2) * q(i,j,3,ic)   * tdt) / &
                      ( mb(i,j)       + cldmas(i,j,2) * tdt)
 
! ============================================================================
!   Compute net change in mixing ratio.                                     
!                                                                           
!   DQ = QB - QC is the total mass to be transported out of the cloud base. 
!   Changes below cloud base are proportional to the background mass.       
!                                                                           
!   Subtract DQ from Q(*,*,K=1,*) and from Q(*,*,K=2,*), but do not make    
!   Q(*,*,K=1,*) or Q(*,*,K=2,*) negative.                                  
! ============================================================================

! modification for now based on GEOS-CHEM
! should be revisited ...

!                 dq = qb(i,j) - qc(i,j)
!                 IF (dq > q(i,j,1,ic) .OR. dq > q(i,j,2,ic)) THEN
!                    q(i,j,2,ic) = qc(i,j)
!                    q(i,j,1,ic) = qc(i,j)
!                 ELSE
!                    q(i,j,2,ic) = q(i,j,2,ic) - dq
!                    q(i,j,1,ic) = q(i,j,1,ic) - dq
!                 END IF

                 q(i,j,2,ic) = qc(i,j)
                 q(i,j,1,ic) = qc(i,j)

! ============================================================================
!   If there is no Cloud mass flux, set QC = Q(K=3) at this I,J location    
! ============================================================================
              ELSE
                 qc(i,j) = q(i,j,3,ic)                 
              END IF
           END DO i_loop_1
        END DO j_loop_1

! ============================================================================
!   (2) Cloud interior mixing                                               
! ============================================================================

        k_loop: DO k = 3,ktop

           j_loop_2: DO j = j1, j2
              i_loop_2: DO i = i1, i2

! ============================================================================
!   If there is cloud mass flux at this location, do the convective         
!   transport.                                                              
!                                                                            
!   QC_PRES = amount of QC preserved against wet scavenging                 
!   CMOUT = air mass flowing out of cloud at level K                        
!   ENTRN = air mass flowing into cloud at level K                          
!                                                                            
!   If Entrainment >= 0 then compute the new value of QC(I,J):              
!                                                                            
!                  CLDMAS(K-1)*QC_PRES + ENTRN(K)*Q(K)                        
!      QC(I,J) =  ---------------------------------------                   
!                     CLDMAS(I,J,K) + DTRN(I,J,K)                           
!
!              =   tracer mass coming in from below      (i.e. level K-1) + 
!                  tracer mass coming in from this level (i.e. level K)
!                 -----------------------------------------------------------
!                             total mass coming into cloud

!                                                                            
!   Entrainment must be >= 0 (since we cannot have a negative flux of air   
!   into the cloud).  This condition is strong enough to ensure that        
!   CMOUT > 0 and will prevent floating-point exception.                    
! ============================================================================

                 IF (cldmas(i,j,k-1) > tiny) THEN

                    ! Soluble species are scavenged during entrainment.
                    qc_pres = qc(i,j) * (1.0 - f(i,j,k,ic))
                    cmout = cldmas(i,j,k) + dtrn(i,j,k)
                    entrn = cmout - cldmas(i,j,k-1)
 
                    IF (entrn >= 0.0) THEN
                       qc(i,j) = ( cldmas(i,j,k-1) * qc_pres + &
                                   entrn           * q(i,j,k,ic) ) / &
                                   cmout
                    endif
 
! ============================================================================
!   The cumulus transport above the cloud base is done as follows:          
!      C_k-1  = cloud air mass flux from level k-1 to level k                   
!      C_k    = cloud air mass flux from level k to level k+1
!      QC_k-1 = mixing ratio of tracer INSIDE CLOUD in level k-1               
!      QC_k   = mixing ratio of tracer INSIDE CLOUD in level k                  
!      Q_k    = mixing ratio of tracer in level k                              
!      Q_k+1  = mixing ratio of tracer in level k+1                            
!                                                                           
!                        (       Cloud        )                             
!                        (                    )                             
!                        (                    )3)      C_k * Q_k+1          
!             k+1        (         ^          )            |                
!            ------------(---------|----------)------------|--------        
!                        (         |          )            V                
!                        (2)   C_k * QC_k     )                             
!                        (                    )                             
!             k          (                    )                             
!                        (         ^          )4)    C_k-1 * Q_k            
!                        (         |          )            |                
!            ------------(---------|----------)------------|----------      
!             k-1        (         |          )            |                
!                        (1) C_k-1 * QC_k-1   )            V                
!                        (         * AP       )                             
!                        (                    )                             
!                        (                    )                             
!                                                                           
!   There are 4 terms that contribute to mass flow in and out of level k:   
!                                                                           
!   1) C_k-1 * QC_PRES = tracer convected from level k-1 to level k       
!   2) C_k   * QC_k    = tracer convected from level k   to level k+1     
!   3) C_k   * Q_k+1   = tracer subsiding from level k+1 to level k       
!   4) C_k-1 * Q_k     = tracer subsiding from level k   to level k-1     
!                                                                           
!   Therefore the change in tracer concentration is given by                
!      DELQ = (Term 1) - (Term 2) + (Term 3) - (Term 4)             
!
!   and Q(I,J,K,IC) = Q(I,J,K,IC) + DELQ        
!                                                                           
! ============================================================================
                    
                    term_1   =  cldmas(i,j,k-1) * qc_pres
                    term_2   = -cldmas(i,j,k  ) * qc(i,j       )
                    term_3   =  cldmas(i,j,k  ) * q (i,j,k+1,ic)
                    term_4   = -cldmas(i,j,k-1) * q (i,j,k,  ic)
                    
                    tsum = term_1 + term_2 + term_3 + term_4 
                    
                    delq = ( tdt / bmass(i,j,k) ) * tsum
                    
                    q(i,j,k,ic) = q(i,j,k,ic) + delq

                    ! prevent concentrations from being negative
                    IF (q(i,j,k,ic) < 1.0E-32) q(i,j,k,ic) = 1.0E-32

                 ELSE
! ============================================================================
!   No cloud transport if cloud mass flux < TINY; Change Qc to q            
! ============================================================================

                    qc(i,j) = q(i,j,k,ic)
                    
                    !--------------------------------------------------
                    ! FIX FOR GEOS-5 MET FIELDS!
                    !
                    ! Bug fix for the cloud base layer, which is not 
                    ! necessarily in the boundary layer, and for 
                    ! GEOS-5, there could be "secondary convection 
                    ! plumes - one in the PBL and another one not.
                    !
                    ! NOTE: TERM_2 and TERM_3 are the same terms as described
                    ! in the above section.
                    !
                    ! (swu, 08/13/2007)
                    !--------------------------------------------------
                    IF ( cldmas(i,j,k) > tiny ) THEN 
                       ! Tracer convected from K -> K+1 
                       term_2   = -cldmas(i,j,k  ) * qc(i,j)
                       
                       ! Tracer subsiding from K+1 -> K 
                       term_3   =  cldmas(i,j,k  ) * q (i,j,k+1,ic)
                       
                       ! Change in tracer concentration
                       delq = ( tdt / bmass(i,j,k) ) * (term_2 + term_3)
                                              
                       ! Add change in tracer to Q array
                       q(i,j,k,ic) = q(i,j,k,ic) + delq
                       
                       ! prevent concentrations from being negative
                       IF (q(i,j,k,ic) < 1.0E-32) q(i,j,k,ic) = 1.0E-32   
                    END IF

                 END IF

              END DO i_loop_2
           END DO j_loop_2
        END DO k_loop

     END DO time_loop
  END DO tracer_loop

END SUBROUTINE cldcnv

end module
