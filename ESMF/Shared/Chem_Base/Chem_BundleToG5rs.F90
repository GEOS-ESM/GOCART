!-------------------------------------------------------------------------
!           NASA GSFC, Global Modeling & Assimilation Office             !
!-------------------------------------------------------------------------
!BOP
!
! !PROGRAM:  rst_hdf_bin.F90 --- Convert an hdf chem restart file to binary,
!                                primarily for obtaining GEOS5/Eros restarts
!                                from GEOS4 fvchem chemistry restart bundles.

PROGRAM rst_hdf_bin

! !USES:

  USE m_die, only: die
  USE Chem_InitMod
  USE Chem_RegistryMod

  IMPLICIT none

! !DESCRIPTION:  Create a binary aerochem restart from an existing HDF
!                restart.
!
! !USAGE:  Typical command line lists the existing HDF file and the 
!          output binary file:
!
!          rst_hdf_bin.x -i b72.rst.chem.19650101_00z.hdf \
!                        -o aerochem_internal_restart
!
!          The Chem_Registry.rc must be applicable to the existing HDF 
!          dataset.  However, the user may supply a list of species that
!          are not to be included in the output file.  See below.
!
!          The binary file is a simple concatenation of 3-D arrays for
!          each specie.
!
!          The program contains no provision for changing horizontal or
!          vertical discretization.
!
! !REVISION HISTORY: 
!
!  15 Jun 2006 Nielsen    First version.
!  28 Sep 2006 da Silva   Major overhaul
!
!EOP
!-------------------------------------------------------------------------
! Grid parameters for existing chem file


   integer, parameter :: READ_ONLY = 1

! This string identifies this procedure

  CHARACTER(LEN=*), PARAMETER :: myname = 'Chem_BundleToG5rs'

! Registry and bundle

  TYPE(Chem_Registry) :: reg	  ! chemistry registry

! Local variables

  CHARACTER(LEN=255) :: argv, inFile, outFile

  REAL, allocatable :: g4(:,:,:)
  real*4, allocatable :: g5(:,:)

  INTEGER :: argc,iargc,ier,k
  INTEGER :: n,nhms,nq,nymd
  INTEGER, PARAMETER :: lu=12
  REAL :: gmin

  INTEGER :: fid, im, jm, km, lm, nvars, ngatts
  INTEGER :: imh, n1, n2, n3, n4, m1, m2, m3, m4


! Parse the command line [see usage() below] for file names
! ---------------------------------------------------------
  argc = iargc()
  if(argc .lt. 4) call usage()
  CALL GetArg(1, inFile)
  CALL GetArg(2, outFile)
  CALL GetArg(3, argv); read(argv,*) nymd
  CALL GetArg(4, argv); read(argv,*) nhms

! Read the chemistry registry
! ---------------------------
  reg = Chem_RegistryCreate( ier )
  IF(ier /= 0) CALL die(myname, 'Cannot create registry')
  CALL Chem_RegistryPrint( reg )
  nq = reg%nq


! Open the GFIO file
! ------------------
  print *, myname //': reading ' // trim(inFile)
  call GFIO_Open ( trim(inFile), READ_ONLY, fid, ier )
  if ( ier /= 0 ) call die(myname,'cannot read '//trim(inFile), ier )

  call GFIO_DimInquire ( fid, im, jm, km, lm, nvars, ngatts, ier)
  if ( ier /= 0 ) call die(myname,'cannot query '//trim(inFile), ier )

  allocate(g4(im,jm,km),g5(im,jm),stat=ier)
  if ( ier /= 0 ) call die(myname,'cannot allocate buffer', ier )

! Read the bundle
! ---------------
  imh = im / 2
  n1 = 1;  n2 = imh;      n3 = n2+1;  n4 = im
  m1 = 1;  m2 = n4-n3+1;  m3 = m2+1;  m4 = im
  open(lu,file=trim(outFile), form='unformatted')
  gmin = tiny(1.0)  ! a very small number
  do n = 1, nq

    call gfio_GetVar ( fid, trim(reg%vname(n)), nymd, nhms, &
                       im, jm, 1, km, g4, ier )
    if ( ier /= 0 ) call die(myname,'cannot read variable', ier )

!   Lower cap for tracers: the smallest float point number
!   ------------------------------------------------------
    where ( g4 < gmin )
            g4 = gmin
    end where

    print 10, trim(Reg%vname(n))//':', minval(g4), maxval(g4), &
              count (g4==gmin), ' points capped'
10 format(1x,a15,1p2e12.3,i8,a)

    do k = 1, km
       g5(m1:m2,:) = g4(n3:n4,:,k)
       g5(m3:m4,:) = g4(n1:n2,:,k)
       write(lu) g5
    end do
  
 end do

 call gfio_close (fid, ier )
 close(lu)
 deallocate(g4,g5)

CONTAINS

 SUBROUTINE usage()
  PRINT *
  PRINT *,'Usage: '
  PRINT *,'  Chem_BundleToG5rs.x inFile outFile nymd nhms'
  PRINT *
  PRINT *,'where both arguments are required and'
  PRINT *
  PRINT *,'  inFile  Input chemistry restart file, HDF format'
  PRINT *,' outFile  Output restart file, binary'
  PRINT *
  CALL exit(1)
 END SUBROUTINE usage

END PROGRAM rst_hdf_bin
