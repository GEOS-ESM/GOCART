!--------------------------------------------------exi-----------------------
!      NASA GSFC, Global Modeling & Assimilation Office, Code 900.3      !
!-------------------------------------------------------------------------
!BOP
!
! !PROGRAM:  SCChem_Init.F90 --- Establishes the inital Chem_Bundle

PROGRAM Chem_initialize

! !USES:

  USE m_die, only: die
  USE Chem_InitMod
  USE Chem_RegistryMod
  USE Chem_BundleMod

  IMPLICIT none

! !DESCRIPTION:  Create a c55 initialization for the SC Code 916 chemistry 
!                bundle under fvGCM version 1_4r2 using
!                  - c55 d_rst file for dynamics from dao_ops
!                  - b55 d_rst file for chemistry from "old" (no SC_GridComp) 
!                    runs.
!
!                The new c55 d_rst file contains the dynamics and specific 
!                humidity from the dao_ops c55 d_rst file, plus the ozone from
!                the ozone from the chemistry.
!                
!                The c55 Chem_Bundle contains 52 species that have been
!                very roughly interpolated from 2x2.5 to 1x1.25, plus RO3OX
!                set to -1.
!
! !USAGE:  Typical command line lists output chem bundle, input dao_ops dynamics,
!          and input chem916 d_rst in the following manner:
!
!    SCChem_Init.x -o /gmao/model/enielsen/c55.rst.chem.19650101_00z.hdf \
!     -d /output/dao_ops/GEOS-4.0.3/a_llk_04/rs/Y2004/M11/a_llk_04.rst.lcv.20041101_00z.bin
!     -c /output/enielsen/chem916/1_4r2/b55/rs/Y1964/M11/b55.rst.lcv.19641101_00z.bin
!
!          For now the output d_rst file is hard-coded.  Sorry.
!
! !REVISION HISTORY: 
!
! 18sep2001  da Silva  Initial code. (?)
! 10Oct2004  Nielsen   Adapted for SCGridComp
! 28Nov2004  Nielsen   Rewritten for getting IC for Polar-AVE mission.
!
!EOP
!-------------------------------------------------------------------------
! Hard-wired c55 grid parameters

  INTEGER, PARAMETER :: im = 288, jm = 181, km = 55

! Hard-wired b55 grid parameters

  INTEGER, PARAMETER :: imb = 144, jmb = 91

! This sting identifies this procedure

  CHARACTER(LEN=*), PARAMETER :: myname = 'SCChem_Init'

! Registry and bundle

  TYPE(Chem_Registry) :: reg	  ! chemistry registry
  TYPE(Chem_Bundle)   :: w_c	  ! chemistry bundle

! Local variables

  CHARACTER(LEN=255) :: argv, chemfile, dummystr, dynfile, b55chemfile

! 1 x 1.25
  REAL :: ps(im,jm)
  REAL :: specie(im,jm,km),u(im,jm,km)
  REAL :: v(im,jm,km),pt(im,jm,km),ro3ox(im,jm,km)

! 2 x 2.5
  REAL :: psb(imb,jmb),delpb(imb,jmb,km)
  REAL :: specieb(imb,jmb,km),ub(imb,jmb,km),qb(imb,jmb,km)
  REAL :: vb(imb,jmb,km),ptb(imb,jmb,km),ro3oxb(imb,jmb,km)

  INTEGER :: argc,i,iarg,iargc,ic,ier
  INTEGER :: j,k,n,nhms,nq,nstep,nymd,prec
  INTEGER, EXTERNAL :: system
  INTEGER, PARAMETER :: iuic = 12

  REAL, POINTER :: delp(:,:,:), q(:,:,:,:)

! Parse the command line [see usage() below] for file names
! of input dynamics restart and output chemistry bundle.
! ---------------------------------------------------------
  argc = iargc()
  if(argc .lt. 1) call usage()
  iarg = 0
  chemfile = 'chem_ini.c_rst'
  do i = 0, 32767
   iarg = iarg+1
   if(iarg .gt. argc) exit
   call GetArg(iarg, argv)
   select case(argv)
    case ("-o")
     if(iarg+1 .gt. argc) call usage()
     iarg = iarg+1
     call GetArg(iarg, chemfile)
    case ("-d")
     if(iarg+1 .gt. argc) call usage()
     iarg = iarg+1
     call GetArg(iarg, dynfile)
    case ("-c")
     if(iarg+1 .gt. argc) call usage()
     iarg = iarg+1
     call GetArg(iarg, b55chemfile)
    case default
     call usage()
   end select
  end do
  if ( trim(chemfile) .eq. 'chem_ini.c_rst' ) then 
   i = index ( dynfile, 'd_rst' )
   if ( i .gt. 1 ) then
    chemfile = dynfile(1:i-1) // 'c_rst' // dynfile(i+5:)
   end if
  end if

! Read the chemistry registry
! ---------------------------
  reg = Chem_RegistryCreate( ier )
  IF(ier /= 0) CALL die(myname, 'Cannot create registry')
  CALL Chem_RegistryPrint( reg )
  nq = reg%nq

! Create the chemistry bundle to be filled in later
! -------------------------------------------------
  CALL Chem_BundleCreate1PE_(reg, im, jm, km, w_c, ier)
  IF(ier /= 0) CALL die(myname, 'Cannot create bundle')

  delp => w_c%delp
  q => w_c%q

! Read the chem916 b55 d_rst file, which contains both dynamics
! and chemistry
! -------------------------------------------------------------
  OPEN(UNIT=iuic,FILE=b55chemfile,FORM='unformatted', &
       ACCESS='sequential',STATUS='old',ACTION='read')

  WRITE(*,FMT="(' CHEM_INIT: Reading b55 ',A,' ...')") TRIM(b55chemfile)
  READ(iuic) nstep,nymd,nhms
  READ(iuic) psb, delpb, ub, vb, ptb

  WRITE(*,FMT="(2X,'PS       b55 MAX=',1PE12.4,'   MIN=',1PE12.4)") &
        MAXVAL(psb),MINVAL(psb)
  WRITE(*,FMT="(2X,'U        b55 MAX=',1PE12.4,'   MIN=',1PE12.4)") &
        MAXVAL(ub),MINVAL(ub)
  WRITE(*,FMT="(2X,'V        b55 MAX=',1PE12.4,'   MIN=',1PE12.4)") &
        MAXVAL(vb),MINVAL(vb)
  WRITE(*,FMT="(2X,'PT       b55 MAX=',1PE12.4,'   MIN=',1PE12.4)") &
        MAXVAL(ptb),MINVAL(ptb)

! Read ozone-to-odd oxygen ratio
! ------------------------------
  READ(iuic) ro3oxb
  WRITE(*,FMT="(2X,'RO3OX    b55 MAX=',1PE12.4,'   MIN=',1PE12.4)") &
        MAXVAL(ro3oxb),MINVAL(ro3oxb)

! Acquire specific humidity
! -------------------------
  READ(iuic) specieb
  qb(:,:,:) = specieb(:,:,:)
  WRITE(*,FMT="(2X,'SPEC HUM b55 MAX=',1PE12.4,'   MIN=',1PE12.4)") &
        MAXVAL(specieb),MINVAL(specieb)

  CALL make1x125(imb,jmb,im,jm,km,qb,specie)
  q(:,:,:,1)=specie(:,:,:)

! Spin through the rest of the species
! ------------------------------------
  DO ic=1,52

   READ(iuic) specieb
   qb(:,:,:) = specieb(:,:,:)
   WRITE(*,FMT="(2X,A8,' b55 MAX=',1PE12.4,'   MIN=',1PE12.4)") &
         reg%vname(ic+2),MAXVAL(specieb),MINVAL(specieb)
   CALL make1x125(imb,jmb,im,jm,km,qb,specie)
   q(:,:,:,ic+2)=specie(:,:,:)
  
  END DO
  
  CLOSE(iuic)

! Just for fun, set the "operational ozone," which we do not use, to
! the sum of the stratospheric and tropospheric parts.  NOTE: vmr!
! ------------------------------------------------------------------
   q(:,:,:,2)=q(:,:,:,37)+q(:,:,:,38)
  
! Add ozone-to-odd oxygen ration at the end
! -----------------------------------------
   q(:,:,:,nq) = -1

! Read in the dynamics fields from dao_ops c55 d_rst file
! --------------------------------------------------------
  OPEN(UNIT=iuic,FILE=dynfile,FORM='unformatted', &
       ACCESS='sequential',STATUS='old',ACTION='read')

  WRITE(*,FMT="(' CHEM_INIT: Reading dao_ops ',A,' ...')") TRIM(dynfile)
  READ(iuic) nstep,nymd,nhms
  READ(iuic) ps, delp, u, v, pt

  WRITE(*,FMT="(2X,'PS        MAX=',1PE12.4,'   MIN=',1PE12.4)") &
        MAXVAL(ps),MINVAL(ps)
  WRITE(*,FMT="(2X,'U         MAX=',1PE12.4,'   MIN=',1PE12.4)") &
        MAXVAL(u),MINVAL(u)
  WRITE(*,FMT="(2X,'V         MAX=',1PE12.4,'   MIN=',1PE12.4)") &
        MAXVAL(v),MINVAL(v)
  WRITE(*,FMT="(2X,'PT        MAX=',1PE12.4,'   MIN=',1PE12.4)") &
        MAXVAL(pt),MINVAL(pt)

! Acquire specific humidity
! -------------------------
  READ(iuic) specie
  q(:,:,:,1) = specie(:,:,:)
  WRITE(*,FMT="(2X,'SPEC HUM  MAX=',1PE12.4,'   MIN=',1PE12.4)") &
        MAXVAL(specie),MINVAL(specie)

! Stat the bundle
! ---------------
!  call Chem_BundleStat(6, w_c, ier)
  if(ier /= 0) call die(myname, 'Cannot stat bundle')

! Write the bundle
! ----------------
  prec = 0
  ier = system("/bin/rm -rf "//chemfile)
  call Chem_BundleWrite(chemfile, nymd, nhms, prec, w_c, ier, verbose=.true.)
  if(ier /=0 ) call die(myname,'Cannot write bundle')

! Write a "new" dynamics restart file.  It must have SPHU and O3
! --------------------------------------------------------------
  dynfile='/gmao/model/enielsen/c55.rst.dyn.20041101_00z.bin'
  OPEN(UNIT=iuic,FILE=TRIM(dynfile),FORM='unformatted', &
       ACCESS='sequential',STATUS='new',ACTION='write')
  WRITE(*,FMT="(' CHEM_INIT: Writing ',A,' ...')") TRIM(dynfile)
  WRITE(iuic) nstep,nymd,nhms
  WRITE(iuic) ps, delp, u, v, pt
  specie(:,:,:)=q(:,:,:,1)
  WRITE(iuic) specie
  specie(:,:,:)=q(:,:,:,2)
  WRITE(iuic) specie
  CLOSE(iuic)
  
! Try to read the newly created bundle
! ------------------------------------
  WRITE(*,FMT="(' CHEM_INIT: Test reading ',A,' ...')") TRIM(chemfile)
  CALL Chem_BundleRead ( TRIM(chemfile), nymd, nhms, w_c, ier, &
                         timidx=0, chemReg=reg )
  WRITE(*,FMT="(' CHEM_INIT: Return code is ',I3)") ier
  
  DO ic=1,55
   WRITE(*,FMT="(2X,A8,' MAX=',1PE12.4,'   MIN=',1PE12.4)") &
         reg%vname(ic),MAXVAL(q(:,:,:,ic)),MINVAL(q(:,:,:,ic))
  END DO
  
CONTAINS

  SUBROUTINE usage()
  PRINT *
  PRINT *,'Usage: '
  PRINT *,'  chem_init.x -o chemfile -d dynfile -c b55chemfile'
  PRINT *
  PRINT *, 'where'
  PRINT *
  PRINT *, '-o chemfile   Output chemistry initialization'
  PRINT *, '		  file name. Default: same as dynfile'
  PRINT *, '		  with substring "d_rst" replaced'
  PRINT *, '		  with "c_rst"'
  PRINT *, '-d dynfile	  Mandatory c55 dao_ops d_rst file name'
  PRINT *, 'b55chemfile	  Mandatory b55 chemistry input file name'
  PRINT *
  CALL exit(1)
  END SUBROUTINE usage
  
  SUBROUTINE make1x125(imb,jmb,im,jm,km,qb,q)
  IMPLICIT none
  INTEGER, INTENT(IN) :: imb,jmb,im,jm,km
  REAL, INTENT(IN) :: qb(imb,jmb,km)
  REAL, INTENT(OUT) :: q(im,jm,km)
  INTEGER :: i,j,k
  REAL :: r(imb,jm,km)
  
  DO k=1,km
   DO i=1,imb
   
    DO j=1,jmb-1
    r(i,2*j-1,k) = qb(i,j,k)
    r(i,2*j  ,k) = (qb(i,j,k)+qb(i,j+1,k))*0.50
    END DO
    r(i,jm,k) = qb(i,jmb,k)

   END DO
  END DO
  
  DO k=1,km
   DO j=1,jm

    DO i=1,imb
    q(2*i-1,j,k) = r(i,j,k)
    q(2*i  ,j,k) = (r(i,j,k)+r(i,j,k))*0.50
    END DO
    q(im   ,j,k) = (r(imb,j,k)+r(1,j,k))*0.50

   END DO
  END DO

  END SUBROUTINE make1x125

END PROGRAM Chem_initialize
