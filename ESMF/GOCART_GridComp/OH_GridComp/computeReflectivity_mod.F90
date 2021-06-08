!-------------------------------------------------------------------------
!BOP
!
! !MODULE: computeReflectivity_mod
!
! !INTERFACE:
!
module computeReflectivity_mod

     USE cblock_size_mod         ! "CMN_SIZE"  ! Size parameters
     USE cblock_CO_mod           ! "CMN_CO"    ! CO arrays

     IMPLICIT NONE
!
! !PUBLIC MEMBER FUNCTIONS:
!
     PRIVATE
     PUBLIC  :: avgrefl

     REAL*8,  SAVE, ALLOCATABLE :: Ravg  (:,:,:)
     INTEGER, SAVE, ALLOCATABLE :: Rcount(:,:,:)
!
! !DESCRIPTION:
!
! !AUTHOR:
! Bryan Duncan
! Jules Kouatchou, Jules.Kouatchou-1@nasa.gov
!
!EOP
!------------------------------------------------------------------------------
CONTAINS
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: avgrefl
!
! !INTERFACE:
!
      SUBROUTINE avgrefl(OPTD, SZA, LPAUSE, Ravga, Ravgb, &
                         nhms, IIPAR, JJPAR, LLPAR)
!*****************************************************************************

      INTEGER, intent(in) :: IIPAR, JJPAR, LLPAR
      INTEGER, intent(in) :: nhms
      REAL*8 , intent(in) :: SZA(IIPAR,JJPAR)
      INTEGER, intent(in) :: LPAUSE(IIPAR,JJPAR)
      REAL*8 , intent(in) :: OPTD(IIPAR, JJPAR, LLPAR)        ! bottom-up

      REAL*8 , intent(inOut) :: Ravga(IIPAR, JJPAR, LLPAR)    ! bottom-up
      REAL*8 , intent(inOut) :: Ravgb(IIPAR, JJPAR, LLPAR)    ! bottom-up

      INTEGER I,J,K,L,M, MOUTPUTON
      INTEGER IJLOOP,IREF,JREF,LHR0,NFNUM
      REAL*8 :: Rout(IIPAR,JJPAR,LLPAR), tempout

      LOGICAL, SAVE :: FIRST = .true.
     
      NFNUM = 65

!*****************************************************************************

!  NEED TO SEE IF DAYLIGHT AVERAGE VS 6-HR AVG AROUND
!  NOON PRODUCES DIFFERENT RESULTS!!!!!!! 

!*****************************************************************************
!
! Designed to convert optical depths to reflectivities
!   and calculate a daylight average which will be used
!   for the CO/OH parameterization option (see SR chemco).
!   Cannot average optical depths, but can average reflectivities.
! Created by Bryan Duncan (1/99).
!
! The "daylight average" can be misleading since, for OH,
!   the optical depths around the solar apex (i.e., noon)
!   are more important in determining its concentration than
!   around dawn or sunset.  Therefore, areas with morning
!   fog (such as marine environments) or areas affected by
!   tropical convection may have a "skewed" daylight average
!   reflectivities. (see Rossow and Schiffer, ISCCP Cloud Data 
!   Products, p.15)
!
! NCHEM is chemical time step (i.e., 24 hours).
! NDT   is meteorological data time step (i.e., 6 hours).
! NTDT  is dynamical time step (i.e., 1/2 hours).
!
!*****************************************************************************
! The annual mean tropopause is stored in the LPAUSE array 
! (from header file "CMN").  LPAUSE is defined such that: 
!
! Levels            1 <= L <= LPAUSE(I,J) - 1 are tropospheric
!         LPAUSE(I,J) <= L <= LLPAR           are stratospheric
!
! (bmy, 4/18/00)  
!*****************************************************************************
!
      IF (FIRST) THEN
         ALLOCATE(Ravg(IIPAR, JJPAR, LLPAR))
         Ravg(:,:,:) = 0.0d0
         ALLOCATE(Rcount(IIPAR, JJPAR, LLPAR))
         Rcount(:,:,:) = 0
         FIRST = .FALSE.
      END IF
!=================================================================
!       
! SR get_refl converts optical depths to reflectivities
!   and then to flux fractions for a given solar zenith angle.
!
      Rout(:,:,:)=0.
      CALL get_refl(OPTD, Rout, SZA, LPAUSE, IIPAR, JJPAR, LLPAR)
     
!
! Sum reflectivities in Ravga() over daylight hours.
! If Rout = -999 then dark, so don't include in sum.
!
      WHERE (Rout(:,:,:).GE.0.) Ravg(:,:,:) = Ravg(:,:,:)+ Rout(:,:,:)

      IF (nhms == 0) THEN 
!
! Average Ravg over NCHEM time step.
!
         WHERE(Rcount(:,:,:).GT.0) &
     &        Ravg(:,:,:)=Ravg(:,:,:)/REAL(Rcount(:,:,:))
!
! If Rcount = 0 and Ravg < 0 then polar night!
!
         WHERE(Rcount(:,:,:).LE.0) Ravg(:,:,:) = -999.
! 
! Convert average reflectivities to flux fractions above and below.
!     
         Ravga(:,:,:)=0.
         Ravgb(:,:,:)=0.
!
! above:
!
         WHERE(Ravg(:,:,:).GT.0.) Ravga(:,:,:)=Ravg(:,:,:)
         WHERE(Ravg(:,:,:).GT.-1..AND.Ravg(:,:,:).LT.0.) Ravga(:,:,:)=0.
         WHERE(Ravg(:,:,:).LT.-1.) Ravga(:,:,:)=-999.
!
! below:
!
         WHERE(Ravga(:,:,:).LT.0.) Ravgb(:,:,:)=-999.

         ! Now L goes up to the annual mean tropopause (bmy, 4/17/00)       
         DO J = 1, JJPAR
            DO I = 1, IIPAR
               DO L = 1, LPAUSE(I,J) - 1
                  IF (Ravga(I,J,L).GT.0.) Ravgb(I,J,L)=Ravg(I,J,1)- &
     &                 Ravga(I,J,L)
               ENDDO
            ENDDO
         ENDDO

         WHERE(Ravgb(:,:,:).GT.-1.AND.Ravgb(:,:,:).LT.0.) &
     &        Ravgb(:,:,:)=0.
!
! EY DEBUG BEGIN
     	 MOUTPUTON=0
         IF(MOUTPUTON.EQ.1) THEN
            open(NFNUM,file='bryanoutput',status='unknown')
            rewind(NFNUM)
            
            DO I=1,IIPAR
            DO J=1,JJPAR
            ! Now L goes up to the annual mean tropopause (bmy, 4/17/00)
            DO L = 1, LPAUSE(I,J) - 1
               write(NFNUM,*)Ravg(I,J,L)
            ENDDO
            ENDDO
            ENDDO

            DO I=1,IIPAR
            DO J=1,JJPAR
            ! Now L goes up to the annual mean tropopause (bmy, 4/17/00)
            DO L = 1, LPAUSE(I,J) - 1
               write(NFNUM,*)Ravga(I,J,L)
            ENDDO
            ENDDO
            ENDDO

            DO I=1,IIPAR
            DO J=1,JJPAR
            ! Now L goes up to the annual mean tropopause (bmy, 4/17/00)
            DO L = 1, LPAUSE(I,J) - 1
               write(NFNUM,*)Ravgb(I,J,L)
            ENDDO
            ENDDO
            ENDDO

            DO I=1,IIPAR
            DO J=1,JJPAR
            ! Now L goes up to the annual mean tropopause (bmy, 4/17/00)
            DO L = 1, LPAUSE(I,J) - 1 
               if(Ravg(I,J,L).gt.0.) THen
                  tempout=(Ravgb(I,J,L)+Ravga(I,J,L))/Ravg(I,J,1)
               else
                  tempout=-999
               endif
               write(NFNUM,*)tempout
            ENDDO
            ENDDO
            ENDDO
  	    write(NFNUM,*), 'OPTD'
            ! Now L goes up to the annual mean tropopause (bmy, 4/17/00)
            DO I=1,IIPAR
            DO J=1,JJPAR
            DO L = 1, LPAUSE(I,J) - 1
               write(NFNUM,*)OPTD(L,I,J)
            ENDDO
            ENDDO
            ENDDO

            DO I=1,IIPAR
            DO J=1,JJPAR
            DO L = 1, LPAUSE(I,J)
               write(NFNUM,*)Rcount(I,J,L)
            ENDDO
            ENDDO
            ENDDO
               
            CLOSE(NFNUM)
         ENDIF ! END DEBUG PRINTOUTS
! EY DEBUG END
         Rcount(:,:,:)=0.
         Ravg(:,:,:)=0.

      ENDIF ! END 24 hour average loop
!
      RETURN

      END SUBROUTINE avgrefl
!EOC
!------------------------------------------------------------------------------
!BOP

      SUBROUTINE get_refl(OPTD, Rout, SZA, LPAUSE, IIPAR, JJPAR, LLPAR)
!*****************************************************************************
! Using optical depth and the cos(SZA), use reflectivity
!   table to find the reflectivity above each point.  
!   To calculate the flux fraction at a point:
!    flux fraction=(flux above)/(flux total from TOA to surface)
! Adapated & modified from chem1d by Bryan Duncan (1/99).
!*****************************************************************************
! The annual mean tropopause is stored in the LPAUSE array 
! (from header file "CMN").  LPAUSE is defined such that: 
!
! Levels            1 <= L <= LPAUSE(I,J) - 1 are tropospheric
!         LPAUSE(I,J) <= L <= LLPAR           are stratospheric
!
! We now use LPAUSE instead of NSKIPL to denote the strat/trop boundary. 
! (bmy, 4/18/00)  
!*****************************************************************************
      IMPLICIT NONE

      INTEGER, intent(in) :: IIPAR, JJPAR, LLPAR
      INTEGER, intent(in) :: LPAUSE(IIPAR,JJPAR)
      REAL*8,  intent(in) :: SZA(IIPAR,JJPAR)
      REAL*8,  intent(in)  :: OPTD(IIPAR,JJPAR,LLPAR)
      REAL*8,  intent(inOut) :: Rout(IIPAR,JJPAR,LLPAR)
!     
      INTEGER NMAX,MMAX,JA,KA
      PARAMETER (NMAX=5,MMAX=5)
      REAL*8 YMTMP(MMAX),YNTMP(NMAX)
      INTEGER ICOL,IINT,JINT,IROW,II,JJ,Q
      INTEGER I,J,K,L
      REAL*8 Y,DY,RFLCVAL,PI,TORAD,BCOSSZA(IIPAR,JJPAR), &
     &     REFINT(3,3),OPTDINT(3),COSTHETAINT(3), &
     &     PTRFLCABOVE,PTRFLCBELOW,PTRFLC, &
     &     PTOPTDEPTH(IIPAR,JJPAR,LLPAR), &
     &     COSTHETA(NREFLCOL),WGHTS(NREFLCOL)
      DATA COSTHETA /0.0469,0.2308,0.5,0.7692,0.9531/
      DATA WGHTS/0.1185,0.2393,0.2844,0.2393,0.1185/
!     
!  To calculate the flux fraction above (Ravga) each point,
!    first sum up the optical depths above. (Optical depths
!    are additive.)  Make the assumption that the point never
!    lies inside a cloud.
!  If the optical depth is greater than 90, then set it =90,
!  which is the maximum optical depth on the chart.
!  Compute cos(SZA) for the point.
!  If the cos(SZA) is beyond the limits of the reflectivity
!  table, set it equal to the limit.
!
!  
      PI=3.14
      TORAD=PI/180.
      BCOSSZA(:,:)=0.
      PTOPTDEPTH(:,:,:)=0.
!
      DO I=1,IIPAR
         DO J=1,JJPAR
            DO L = 1, LPAUSE(I,J) - 1
!  EY SEE HERE!!!!!
               PTOPTDEPTH(I,J,L)=SUM(OPTD(I,J,L+1:LLPAR))
               IF(PTOPTDEPTH(I,J,L).GT.90.) PTOPTDEPTH(I,J,L)=90.
!     
            ENDDO
         ENDDO
      ENDDO
!  
      BCOSSZA(:,:)= COS(SZA(:,:)*TORAD)
      WHERE(BCOSSZA(:,:).LT.0.0469)BCOSSZA(:,:)=0.0469
      WHERE(BCOSSZA(:,:).GT.0.9531)BCOSSZA(:,:)=0.9531
!
!*****************************************************************************
! Begin looping over IIPAR and JJPAR
!*****************************************************************************
!
!  Now we have both optical depth and cos(SZA), and can use cms'
!  reflectivity table to find the reflectivity of point.  note:
!  need to interpolate in two directions in the reflectivity table.
!
      DO I=1,IIPAR
         DO J=1,JJPAR
            DO L = 1, LPAUSE(I,J) - 1
!
!  If the solar zenith angle is > 85, then set PTRFLC=0.      
!
               IF (SZA(I,J).GT.85.) THEN
                  PTRFLC=-999.
                  GOTO 200 
               ENDIF
!
!  Calculated reflectivity is stored in PTRFLC.
!
               PTRFLC=0.
!
!  If the solar zenith angle is < 85, then update Rcount.
!
               Rcount(I,J,L)=Rcount(I,J,L)+1
!
!  Loop across columns, looking for the correct cos(SZA) values
!  to interpolate between.  note:  the columns correspond to the
!  quadrature points.
!     
               DO ICOL=1,NREFLCOL-1
                  IF ((BCOSSZA(I,J).GE.COSTHETA(ICOL)).AND. &
     &                (BCOSSZA(I,J).LE.COSTHETA(ICOL+1))) THEN 
                      JINT=ICOL
                      GOTO 230
                  ENDIF
               ENDDO
!     
 230           IF ((BCOSSZA(I,J)-COSTHETA(NREFLCOL-1).LE.&
                  (COSTHETA(NREFLCOL)-BCOSSZA(I,J))).AND.(JINT.GT.1)) THEN
                  JINT=JINT-1
               ENDIF
               IF(JINT.GT.3) JINT=3
!     
! Look down columns for correct optical depth to interpolate between.
!
               DO IROW=1,NREFL-1
                  IF ((PTOPTDEPTH(I,J,L).GE.OPTDEPTH(IROW)).AND. &
     &                (PTOPTDEPTH(I,J,L).LE.OPTDEPTH(IROW+1))) THEN
                     IINT=IROW
                     GOTO 330
                  ENDIF 
               ENDDO
!
 330           IF (((PTOPTDEPTH(I,J,L)-OPTDEPTH(NREFL-1)).LE. &
     &              (OPTDEPTH(NREFL)-PTOPTDEPTH(I,J,L))).AND. &
     &              (IINT.GT.1)) THEN 
                  IINT=IINT-1
               ENDIF
               IF (IINT.GT.(NREFL-2)) THEN
                  IINT=NREFL-2
               ENDIF 
!     
               DO II=1,3
                  OPTDINT(II)=OPTDEPTH(IINT+II-1)
                  COSTHETAINT(II)=COSTHETA(JINT+II-1)
                  DO JJ=1,3
                     REFINT(II,JJ)=RFLC(IINT+II-1,JINT+JJ-1)
                  ENDDO
               ENDDO
!     
! ****************************************************************************
! Interpolation
!  Calculates an interpolated function value y
!  and an accuracy indication dy.
!
               DO JA=1,3
                  DO KA=1,3
                     YNTMP(KA)=REFINT(JA,KA)
                  END DO
            
                  CALL POLINT(COSTHETAINT,YNTMP,3,BCOSSZA(I,J),YMTMP(JA),DY)
               END DO
            
               CALL POLINT(OPTDINT,YMTMP,3,PTOPTDEPTH(I,J,L),Y,DY)
!
! ****************************************************************************
!
               RFLCVAL=Y
!
! Now use the quadrature weights to get the final reflectivity.
!
               DO 120 Q=1,5
 120              PTRFLC=PTRFLC+WGHTS(Q)*RFLCVAL
!
 200           CONTINUE
!*****************************************************************************
! Reflectivity above is stored in RAVGA.
!  The total reflectivity from the surface to TOA is stored in
!  ROUT(I,J,1).
               Rout(I,J,L)=PTRFLC

!*****************************************************************************
            ENDDO
         ENDDO
      ENDDO
!*****************************************************************************
!
!      print*, ' TOTAL reflectivity to TOA ', Rout(:,:,1) 
      RETURN
      END SUBROUTINE get_refl

end module computeReflectivity_mod
