!-------------------------------------------------------------------------
!      NASA/GSFC, Global Modeling & Assimilation Office, Code 900.3      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  Chem_RegistryMod --- Chemistry Registry Class (Stub Version)
!
! !INTERFACE:
!

   module  Chem_RegistryMod

! !USES:

   implicit none

! !PUBLIC TYPES:
!
   PRIVATE
   PUBLIC  Chem_Registry     ! Keeps track of which constituents are active,
                             !  to become internal state of Grid Component
                           
!
! !PUBLIIC MEMBER FUNCTIONS:
!
   PUBLIC  Chem_RegistryCreate   ! Constructor from RC file
   PUBLIC  Chem_RegistryDestroy  ! Destructor
   PUBLIC  Chem_RegistryPrint    ! Prints a summary of the Chemistry registry

   PUBLIC  Chem_RegistrySetIsGOCART ! Whether variable belongs to GOCART 

!
! !DESCRIPTION:
!
!  This module implements a registry for (chemical) constituents.
!  This initial class is intended to serve as a stop gap before an ESMF 
!  implementation is adopted. 
!
!
! !REVISION HISTORY:
!
!  21Oct2009 da Silva  Derived from real Chem_Registry for use at NCEP.
!
!EOP
!-------------------------------------------------------------------------

  integer, parameter :: nch = 255

! Registry
! --------
  type Chem_Registry

     integer :: nq    ! Total number of tracers 

!    Fixed Tracers
!    -------------
     logical :: doing_H2O = .true.  ! water vapor
     logical :: doing_O3  = .true.  ! ozone

     character(len=nch) :: units_H2O = 'kg/kg'
     character(len=nch) :: units_O3  = 'ppmv'

     integer :: n_H2O = 1, i_H2O = 1, j_H2O = 1
     integer :: n_O3  = 1, i_O3  = 2, j_O3  = 2

!    Floating Tracers
!    ----------------
     logical :: doing_CO    ! carbon monoxide
     logical :: doing_CO2   ! carbon dioxide
     logical :: doing_DU    ! mineral dust
     logical :: doing_SS    ! sea salt
     logical :: doing_SU    ! sulfates
     logical :: doing_CFC   ! CFCs
     logical :: doing_BC    ! black carbon
     logical :: doing_OC    ! organic carbon
     logical :: doing_BRC   ! brown carbon
     logical :: doing_Rn    ! radon
     logical :: doing_CH4   ! Methane
     logical :: doing_SC    ! stratospheric chemistry
     logical :: doing_XX    ! ancillary data
     logical :: doing_AC    ! auto chem 
     logical :: doing_PC    ! Parameterized Chemistry (GEOS-5)
     logical :: doing_GMI   ! GMI Chemistry (GEOS-5)
     logical :: doing_OCS   ! ACHEM chemistry (OCS)
     logical :: doing_NI    ! Nitrate
     logical :: doing_TR    ! passive tracers

!    Number of bins and tracer index ranges for each constituent:
!        n_TT - number of bins for tracer TT (n_TT = j_TT - i_TT + 1)
!        i_TT - first index for tracer TT
!        j_TT - last  index for tracer TT
!    -----------------------------------------------------------
     integer :: n_CO, i_CO, j_CO     ! carbon monoxide
     integer :: n_CO2,i_CO2,j_CO2    ! carbon dioxide
     integer :: n_DU, i_DU, j_DU     ! mineral dust
     integer :: n_SS, i_SS, j_SS     ! sea salt
     integer :: n_SU, i_SU, j_SU     ! sulfates
     integer :: n_CFC,i_CFC,j_CFC    ! CFCs
     integer :: n_BC, i_BC, j_BC     ! black carbon
     integer :: n_OC, i_OC, j_OC     ! organic carbon
     integer :: n_BRC, i_BRC, j_BRC  ! brown carbon
     integer :: n_Rn, i_Rn, j_Rn     ! radon
     integer :: n_CH4,i_CH4,j_CH4    ! Methane
     integer :: n_SC, i_SC, j_SC     ! stratospheric chemistry
     integer :: n_XX, i_XX, j_XX     ! ancillary data
     integer :: n_AC, i_AC, j_AC     ! auto chem
     integer :: n_PC, i_PC, j_PC     ! parameterized chemistry (GEOS-5)
     integer :: n_GMI, i_GMI, j_GMI  ! GMI chemistry (GEOS-5)
     integer :: n_OCS, i_OCS, j_OCS  ! OCS chemistry (ACHEM)
     integer :: n_NI, i_NI, j_NI     ! Nitrate
     integer :: n_TR, i_TR, j_TR     ! passive tracers

!    GEOS-5 Short-hands: all combined tracers from CO to OC
!    ------------------------------------------------------
     logical :: doing_GOCART  
     integer :: n_GOCART, i_GOCART, j_GOCART

!    Tracer units
!    ------------
     character(len=nch) :: units_CO    ! carbon monoxide
     character(len=nch) :: units_CO2   ! carbon dioxide
     character(len=nch) :: units_DU    ! mineral dust
     character(len=nch) :: units_SS    ! sea salt
     character(len=nch) :: units_SU    ! sulfates
     character(len=nch) :: units_CFC   ! CFCs
     character(len=nch) :: units_BC    ! black carbon
     character(len=nch) :: units_OC    ! organic carbon
     character(len=nch) :: units_BRC   ! brown carbon
     character(len=nch) :: units_Rn    ! radon
     character(len=nch) :: units_CH4   ! Methane
     character(len=nch) :: units_SC    ! stratospheric chemistry
     character(len=nch) :: units_XX    ! ancillary data
     character(len=nch) :: units_AC    ! auto chem
     character(len=nch) :: units_PC    ! parameterized chemistry (GEOS-5)
     character(len=nch) :: units_GMI   ! GMI chemistry (GEOS-5)
     character(len=nch) :: units_OCS   ! OCS chemistry (ACHEM)
     character(len=nch) :: units_NI    ! Nitrate
     character(len=nch) :: units_TR    ! passive tracers

!    CF Style metadata
!    -----------------
     character(len=nch), pointer :: vname(:)   ! (nq), variable short name
     character(len=nch), pointer :: vtitle(:)  ! (nq), variable long  name
     character(len=nch), pointer :: vunits(:)  ! (nq), variable units

!    Tracer transport properties
!    ---------------------------
!!!  logical, pointer :: advect(:)  ! (nq), whether to advect it
!!!  logical, pointer :: diffuse(:) ! (nq), whether to diffuse it
!    Set (or not) from component resource files
     real, pointer    :: fscav(:)   ! (nq), scavenging coefficient
     real, pointer    :: rhop(:)    ! (nq), dry particle mass density [kg m-3]
     real, pointer    :: molwght(:) ! (nq), molecular weight [kg mole-1]
     real, pointer    :: rlow(:)    ! (nq), lower edge of particle size bin [m]
     real, pointer    :: rup(:)     ! (nq), upper edge of particle size bin [m]
     real, pointer    :: rmed(:)    ! (nq), particle bin number median radius [m]
     real, pointer    :: sigma(:)   ! (nq), particle lognormal width
     real, pointer    :: fNum(:)    ! (nq), ratio of particle number to mass

  end type Chem_Registry

CONTAINS

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Chem_RegistryCreate --- Construct Chemistry Registry
!
! !INTERFACE:
!

  Function Chem_RegistryCreate ( rc, rcfile )

  implicit none
  type(Chem_Registry) Chem_RegistryCreate 

! !USES:

! !INPUT PARAMETERS:

   character(len=*), OPTIONAL :: rcfile  ! Resource file name; default is
                                         ! 'Chem_Registry.rc'

! !OUTPUT PARAMETERS:

   integer, intent(out) ::  rc            ! Error return code:
                                          !  0 - all is well
                                          !  1 - 

! !DESCRIPTION:  Disable all tracers.
!
!
! !REVISION HISTORY:
!
!  21Oct2009 da Silva  Derived from real Chem_Registry for use at NCEP.
!
!EOP
!-------------------------------------------------------------------------

   character(len=*), parameter ::  myname = 'Chem_RegistryCreate'

   type(Chem_Registry) :: this
   character(len=255) :: rcfilen
   integer :: nq, ios, ier, n
   logical, allocatable :: isGOCART(:)

   rc = 0
                
!  Defaults
!  --------
   nq = 0

!  ------------------------------------------------------
!  Parse resource file to see which tracers are on
!     defines doing_xx, and n_xx for each tracer
!  ------------------------------------------------------
   call parserc_ ( 'H2O', this%doing_H2O, this%n_H2O, this%units_H2O )
   call parserc_ ( 'O3', this%doing_O3, this%n_O3, this%units_O3 )
   call parserc_ ( 'CO', this%doing_CO, this%n_CO, this%units_CO )
   call parserc_ ( 'CO2', this%doing_CO2, this%n_CO2, this%units_CO2 )
   call parserc_ ( 'DU', this%doing_DU, this%n_DU, this%units_DU )
   call parserc_ ( 'SS', this%doing_SS, this%n_SS, this%units_SS )
   call parserc_ ( 'SU', this%doing_SU, this%n_SU, this%units_SU )
   call parserc_ ( 'CFC', this%doing_CFC, this%n_CFC, this%units_CFC )
   call parserc_ ( 'BC', this%doing_BC, this%n_BC, this%units_BC )
   call parserc_ ( 'OC', this%doing_OC, this%n_OC, this%units_OC )
   call parserc_ ( 'BRC', this%doing_BRC, this%n_BRC, this%units_BRC )
   call parserc_ ( 'Rn', this%doing_Rn, this%n_Rn, this%units_Rn )
   call parserc_ ( 'CH4', this%doing_CH4, this%n_CH4, this%units_CH4 )
   call parserc_ ( 'SC', this%doing_SC, this%n_SC, this%units_SC )
   call parserc_ ( 'GMI', this%doing_GMI, this%n_GMI, this%units_GMI )
   call parserc_ ( 'XX', this%doing_XX, this%n_XX, this%units_XX )
   call parserc_ ( 'AC', this%doing_AC, this%n_AC, this%units_AC )
   call parserc_ ( 'PC', this%doing_PC, this%n_PC, this%units_PC )
   call parserc_ ( 'OCS', this%doing_OCS, this%n_OCS, this%units_OCS )
   call parserc_ ( 'NI', this%doing_NI, this%n_NI, this%units_NI )
   call parserc_ ( 'TR', this%doing_TR, this%n_TR, this%units_TR )

!  Set internal indices
!  --------------------
   call setidx_ ( this%doing_H2O, this%n_H2O, this%i_H2O, this%j_H2O )
   call setidx_ ( this%doing_O3, this%n_O3, this%i_O3, this%j_O3 )
   call setidx_ ( this%doing_CO, this%n_CO, this%i_CO, this%j_CO )
   call setidx_ ( this%doing_CO2, this%n_CO2, this%i_CO2, this%j_CO2 )
   call setidx_ ( this%doing_DU, this%n_DU, this%i_DU, this%j_DU )
   call setidx_ ( this%doing_SS, this%n_SS, this%i_SS, this%j_SS )
   call setidx_ ( this%doing_SU, this%n_SU, this%i_SU, this%j_SU )
   call setidx_ ( this%doing_CFC, this%n_CFC, this%i_CFC, this%j_CFC )
   call setidx_ ( this%doing_BC, this%n_BC, this%i_BC, this%j_BC )
   call setidx_ ( this%doing_OC, this%n_OC, this%i_OC, this%j_OC )
   call setidx_ ( this%doing_BRC, this%n_BRC, this%i_BRC, this%j_BRC )
   call setidx_ ( this%doing_Rn, this%n_Rn, this%i_Rn, this%j_Rn )
   call setidx_ ( this%doing_CH4, this%n_CH4, this%i_CH4, this%j_CH4 )
   call setidx_ ( this%doing_SC, this%n_SC, this%i_SC, this%j_SC )
   call setidx_ ( this%doing_GMI, this%n_GMI, this%i_GMI, this%j_GMI )
   call setidx_ ( this%doing_XX, this%n_XX, this%i_XX, this%j_XX )
   call setidx_ ( this%doing_AC, this%n_AC, this%i_AC, this%j_AC )
   call setidx_ ( this%doing_PC, this%n_PC, this%i_PC, this%j_PC )
   call setidx_ ( this%doing_OCS, this%n_OCS, this%i_OCS, this%j_OCS )
   call setidx_ ( this%doing_NI, this%n_NI, this%i_NI, this%j_NI )
   call setidx_ ( this%doing_TR, this%n_TR, this%i_TR, this%j_TR )

!  Allocate memory in registry
!  ---------------------------
   this%nq = nq
   allocate ( this%vname(nq), this%vtitle(nq), this%vunits(nq), &
              this%fscav(nq), this%rhop(nq), this%molwght(nq), &
              this%rlow(nq), this%rup(nq), this%rmed(nq), &
              this%sigma(nq), this%fNum(nq), stat=ios )
   if ( ios /= 0 ) then
        rc = 2
        return 
   end if

   this%fscav    = 0.0 ! no scavanging by default
   this%rhop     = -1. ! default
   this%molwght  = -1. ! default
   this%rlow     = -1. ! default
   this%rup      = -1. ! default
   this%rmed     = -1. ! default
   this%sigma    = -1. ! default
   this%fNum     = -1. ! default

!  Fill in CF metadata
!  -------------------
   call setmeta_ ( this%doing_H2O, 'q ', 'Specific Humidity', & 
                   this%units_H2O, this%i_H2O, this%j_H2O )
   call setmeta_ ( this%doing_O3,  'o3', 'Ozone Mixing Ratio', &
                   this%units_O3, this%i_O3, this%j_O3 )
   call setmeta_ ( this%doing_CO,  'CO', 'Carbon Monoxide Mixing Ratio', &
                   this%units_CO,  this%i_CO, this%j_CO )
   call setmeta_ ( this%doing_CO2,  'CO2', 'Carbon Dioxide Mixing Ratio', &
                   this%units_CO2,  this%i_CO2, this%j_CO2 )
   call setmeta_ ( this%doing_DU,  'du', 'Dust Mixing Ratio', &
                   this%units_DU,  this%i_DU, this%j_DU )
   call setmeta_ ( this%doing_SS,  'ss', 'Sea Salt Mixing Ratio', & 
                   this%units_SS,  this%i_SS, this%j_SS )
   call setmeta_ ( this%doing_SU,  'su', 'Surfates Mixing Ratio', &
                   this%units_SU,  this%i_SU, this%j_SU )
   call setmeta_ ( this%doing_CFC,  'CFC', 'CFC-12 (CCl2F2) Mixing Ratio', &
                   this%units_CFC,  this%i_CFC, this%j_CFC )
   call setmeta_ ( this%doing_BC,  'bc', 'Black Carbon Mixing Ratio', &
                   this%units_BC,  this%i_BC, this%j_BC )
   call setmeta_ ( this%doing_OC,  'oc', 'Organic Carbon Mixing Ratio', &
                   this%units_OC,  this%i_OC, this%j_OC )
   call setmeta_ ( this%doing_BRC,  'brc', 'Brown Carbon Mixing Ratio', &
                   this%units_BRC,  this%i_BRC, this%j_BRC )
   call setmeta_ ( this%doing_Rn,  'Rn', 'Radon Mixing Ratio', &
                   this%units_Rn,  this%i_Rn, this%j_Rn )
   call setmeta_ ( this%doing_CH4,  'CH4', 'Methane Mixing Ratio', &
                   this%units_CH4,  this%i_CH4, this%j_CH4 )
   call setmeta_ ( this%doing_SC,  'sc', 'Stratosperic Chemistry Species', &
                   this%units_SC,  this%i_SC, this%j_SC )
   call setmeta_ ( this%doing_GMI,  'GMI', 'GMI Chemistry', &
                   this%units_GMI,  this%i_GMI, this%j_GMI )
   call setmeta_ ( this%doing_XX,  'xx', 'Ancillary Data', &
                   this%units_XX,  this%i_XX, this%j_XX )
   call setmeta_ ( this%doing_AC,  'ac', 'Auto Chemistry Species', &
                   this%units_AC,  this%i_AC, this%j_AC )
   call setmeta_ ( this%doing_PC,  'pc', 'Parameterized Chemistry', &
                   this%units_PC,  this%i_PC, this%j_PC )
   call setmeta_ ( this%doing_OCS,  'ocs', 'Carbonyl Sulfide', &
                   this%units_OCS,  this%i_OCS, this%j_OCS )
   call setmeta_ ( this%doing_NI,  'ni', 'Nitrate', &
                   this%units_NI,  this%i_NI, this%j_NI )
   call setmeta_ ( this%doing_TR,  'TR', 'Passive Tracers', &
                   this%units_TR,  this%i_TR, this%j_TR )
		   
!  Set indices for the GOCART family: from CO to OC
!  ------------------------------------------------
   allocate ( isGOCART(nq), stat=ios )
   if ( ios /= 0 ) then
        rc = 3
        return
   end if
   call Chem_RegistrySetIsGOCART ( this, isGOCART, nq )   
   if ( any(isGOCART) ) then
      this%doing_GOCART = .true.
      do n = 1, nq
         if ( isGOCART(n) ) then
            this%i_GOCART = n
            exit
         end if
      end do
      do n = nq, 1, -1
         if ( isGOCART(n) ) then
            this%j_GOCART = n
            exit
         end if
      end do
      this%n_GOCART = this%j_GOCART - this%i_GOCART + 1 
   else
      this%doing_GOCART = .false.
      this%n_GOCART =  0
      this%i_GOCART = -1
      this%j_GOCART = -2
   end if
   deallocate ( isGOCART )

!  All done
!  --------
   Chem_RegistryCreate = this
   
   return 

!                 -----------------------------
!                 Internal Constructor Routines
!                 -----------------------------

   CONTAINS

      subroutine parserc_ ( name, doing_it, n_tt, units ) ! parses rc file
!     -------------------
      character(len=*), intent(in) :: name
      logical, intent(out)  :: doing_it
      integer, intent(out)  :: n_tt  ! number of bins for tracer 
      character(len=*), intent(out) :: units

      character(len=255) :: answer
      integer ier

!     Defaults
!     --------
      doing_it = .false.
      n_tt = -1 
      units = 'unknown'

      end subroutine parserc_ 
      

      subroutine setidx_ ( doing_it, n_tt, i_tt, j_tt ) ! set tracer indices
!     ------------------
      logical, intent(in)  :: doing_it
      integer, intent(in)  :: n_tt  ! number of bins for tracer 
      integer, intent(out) :: i_tt  ! first tracer index
      integer, intent(out) :: j_tt  ! last tracer index
         i_tt = -1
         j_tt = -2
      end subroutine setidx_


      subroutine setmeta_ ( doing_it, vname, vtitle, vunits, i_tt, j_tt )
!     -------------------
      logical, intent(in) :: doing_it
      character(len=*), intent(in) :: vname, vtitle, vunits
      integer, intent(in) :: i_tt, j_tt
      integer i, nbins, ibin, ier, n
      character(len=3) :: cbin
      character(len=255) :: token
      character(len=255) :: uvname

      return

     end subroutine setmeta_

 end Function Chem_RegistryCreate

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Chem_RegistryDestroy --- Destruct Chemisty Registry
!
! !INTERFACE:
!
  subroutine Chem_RegistryDestroy ( this, rc )

! !USES:

  implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(Chem_Registry), intent(inout) :: this

! !OUTPUT PARAMETERS:

  integer, intent(out)          ::  rc     ! Error return code:
                                           !  0 - all is well
                                           !  1 - 

! !DESCRIPTION: Destructor for Chemistry Registry object.
!
! !REVISION HISTORY:
!
!  22Jul2003 da Silva  First crack.
!  03Sep2004 da Silva  Added stratospheric chemistry hooks.
!
!EOP
!-------------------------------------------------------------------------
   integer ios

   rc = 0
   this%nq = -1 
   this%doing_H2O = .false.
   this%doing_O3  = .false.
   this%doing_CO = .false.    ! carbon monoxide
   this%doing_CO2 = .false.   ! carbon dioxide
   this%doing_DU = .false.    ! mineral dust
   this%doing_SS = .false.    ! sea salt
   this%doing_SU = .false.    ! sulfates
   this%doing_CFC = .false.   ! CFCs
   this%doing_BC = .false.    ! black carbon
   this%doing_OC = .false.    ! organic carbon
   this%doing_BRC = .false.   ! brown carbon
   this%doing_Rn = .false.    ! radon
   this%doing_CH4 = .false.   ! Methane
   this%doing_SC = .false.    ! stratospheric chemistry
   this%doing_AC = .false.    ! stratospheric chemistry
   this%doing_XX = .false.    ! ancillary data
   this%doing_PC = .false.    ! parameterized chemistry (GEOS-5)
   this%doing_OCS = .false.   ! ACHEM chemistry (OCS)
   this%doing_NI = .false.    ! Nitrate
   this%doing_GMI = .false.   ! GMI chemistry (GEOS-5)
   this%doing_TR = .false.    ! passive tracers
   deallocate ( this%vname, this%vtitle, this%vunits, this%fscav, &
                this%rhop, this%molwght, this%rlow, this%rup, this%rmed, &
                this%sigma, this%fNum, stat=ios )
   if ( ios /= 0 ) then
        rc = 1
        return 
   end if

end subroutine Chem_RegistryDestroy 

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Chem_RegistryPrint --- Print summary of Chemistry Registry
!
! !INTERFACE:
!
   SUBROUTINE Chem_RegistryPrint( reg )

! !USES:

! !INPUT PARAMETERS:

   IMPLICIT none
   TYPE(Chem_Registry) :: reg

! !OUTPUT PARAMETERS:


! !DESCRIPTION:
!
!   Prints summary of Chemistry Registry
!
! !REVISION HISTORY:
!
!  22Jul2003 da Silva  First crack.
!   9Dec2004 Nielsen   Enhancements.
!
!EOP
!-------------------------------------------------------------------------

   CHARACTER (LEN=255) :: ActiveList
   
   PRINT *
   PRINT *,'--------------------------------------------------------'
   PRINT *,'|          Summary of the Chemistry Registry           |'
   PRINT *,'|               from Chem_RegistryPrint                |'
   PRINT *,'--------------------------------------------------------'
   WRITE(*,FMT="(' ','       Number of species: ',I3)") reg%nq

   ActiveList = '  '
   IF ( reg%doing_H2O ) ActiveList = TRIM(ActiveList)//'  H2O'
   IF ( reg%doing_O3 ) ActiveList = TRIM(ActiveList)//'  O3'
   IF ( reg%doing_CO ) ActiveList = TRIM(ActiveList)//'  CO'
   IF ( reg%doing_CO2 ) ActiveList = TRIM(ActiveList)//'  CO2'
   IF ( reg%doing_DU ) ActiveList = TRIM(ActiveList)//'  DU'
   IF ( reg%doing_SS ) ActiveList = TRIM(ActiveList)//'  SS'
   IF ( reg%doing_SU ) ActiveList = TRIM(ActiveList)//'  SU'
   IF ( reg%doing_CFC ) ActiveList = TRIM(ActiveList)//'  CFC'
   IF ( reg%doing_BC ) ActiveList = TRIM(ActiveList)//'  BC'
   IF ( reg%doing_BRC ) ActiveList = TRIM(ActiveList)//'  BRC'
   IF ( reg%doing_OC ) ActiveList = TRIM(ActiveList)//'  OC'
   IF ( reg%doing_Rn ) ActiveList = TRIM(ActiveList)//'  Rn'
   IF ( reg%doing_CH4 ) ActiveList = TRIM(ActiveList)//'  CH4'
   IF ( reg%doing_SC ) ActiveList = TRIM(ActiveList)//'  SC'
   IF ( reg%doing_GMI ) ActiveList = TRIM(ActiveList)//'  GMI'
   IF ( reg%doing_AC ) ActiveList = TRIM(ActiveList)//'  AC'
   IF ( reg%doing_XX ) ActiveList = TRIM(ActiveList)//'  XX'
   IF ( reg%doing_PC ) ActiveList = TRIM(ActiveList)//'  PC'
   IF ( reg%doing_OCS ) ActiveList = TRIM(ActiveList)//'  OCS'
   IF ( reg%doing_NI ) ActiveList = TRIM(ActiveList)//'  NI'
   IF ( reg%doing_TR ) ActiveList = TRIM(ActiveList)//'  TR'
   
   PRINT *
   PRINT *, 'Active chemistry components:',TRIM(ActiveList)
   
   IF ( reg%doing_H2O ) CALL reg_prt_( 'H2O', reg%n_H2O, reg%i_H2O, reg%j_H2O )
   IF ( reg%doing_O3 ) CALL reg_prt_( 'O3', reg%n_O3, reg%i_O3, reg%j_O3 )
   IF ( reg%doing_CO ) CALL reg_prt_( 'CO', reg%n_CO, reg%i_CO, reg%j_CO ) 
   IF ( reg%doing_CO2 ) CALL reg_prt_( 'CO2', reg%n_CO2, reg%i_CO2, reg%j_CO2 )
   IF ( reg%doing_DU ) CALL reg_prt_( 'DU', reg%n_DU, reg%i_DU, reg%j_DU )
   IF ( reg%doing_SS ) CALL reg_prt_( 'SS', reg%n_SS, reg%i_SS, reg%j_SS )
   IF ( reg%doing_SU ) CALL reg_prt_( 'SU', reg%n_SU, reg%i_SU, reg%j_SU )
   IF ( reg%doing_CFC ) CALL reg_prt_( 'CFC', reg%n_CFC, reg%i_CFC, reg%j_CFC )
   IF ( reg%doing_BC ) CALL reg_prt_( 'BC', reg%n_BC, reg%i_BC, reg%j_BC )
   IF ( reg%doing_OC ) CALL reg_prt_( 'OC', reg%n_OC, reg%i_OC, reg%j_OC )
   IF ( reg%doing_BRC ) CALL reg_prt_( 'BRC', reg%n_BRC, reg%i_BRC, reg%j_BRC )
   IF ( reg%doing_Rn ) CALL reg_prt_( 'Rn', reg%n_Rn, reg%i_Rn, reg%j_Rn )
   IF ( reg%doing_CH4 ) CALL reg_prt_( 'CH4', reg%n_CH4, reg%i_CH4, reg%j_CH4 )
   IF ( reg%doing_SC ) CALL reg_prt_( 'SC', reg%n_SC, reg%i_SC, reg%j_SC )
   IF ( reg%doing_GMI ) CALL reg_prt_( 'GMI', reg%n_GMI, reg%i_GMI, reg%j_GMI )
   IF ( reg%doing_AC ) CALL reg_prt_( 'AC', reg%n_AC, reg%i_AC, reg%j_AC )
   IF ( reg%doing_XX ) CALL reg_prt_( 'XX', reg%n_XX, reg%i_XX, reg%j_XX )
   IF ( reg%doing_PC ) CALL reg_prt_( 'PC', reg%n_PC, reg%i_PC, reg%j_PC )
   IF ( reg%doing_OCS ) CALL reg_prt_( 'OCS', reg%n_OCS, reg%i_OCS, reg%j_OCS )
   IF ( reg%doing_NI ) CALL reg_prt_( 'NI', reg%n_NI, reg%i_NI, reg%j_NI )
   IF ( reg%doing_TR ) CALL reg_prt_( 'TR', reg%n_TR, reg%i_TR, reg%j_TR )

   IF ( reg%doing_GOCART ) & 
        CALL reg_prt_( 'GOCART is a COMPOSITE and', &
             reg%n_GOCART, reg%i_GOCART, reg%j_GOCART )

   PRINT *

!  PRINT *,'----- End of the summary of the Chemistry Registry -----'
!  PRINT *
   
  101 FORMAT(/,'       Number of species: ',I3)

   RETURN
   
   CONTAINS
   
      SUBROUTINE reg_prt_ ( compName, n, i1, i2 )
      
      IMPLICIT none
      CHARACTER(LEN=*), INTENT(IN) :: compName
      INTEGER, INTENT(IN) :: n, i1, i2
      INTEGER :: i
      CHARACTER(LEN=7) :: string
      
      string = 'species'
      IF( n == 1 ) string='specie '
      
      WRITE(*,101) TRIM(compName),n,string
      DO i = i1, i2
       WRITE(*,201) i,TRIM(reg%vname(i)),TRIM(reg%vunits(i)),TRIM(reg%vtitle(i))
      END DO

  101 FORMAT(/,' Component ',A,' has ',I3,' ',A7,/, &
	     ' No ',2X,'  Name  ',2X,'   Units  ',2X,'Description',/, &
	     ' ---',2X,'--------',2X,'----------',2X,'-----------')
  201 FORMAT(' ',I3,2X,A8,2X,A10,2X,A)
          
  END SUBROUTINE reg_prt_
  
END SUBROUTINE Chem_RegistryPrint

   subroutine Chem_RegistrySetIsGOCART ( chemReg, isGOCART, nq )
     type(Chem_registry), intent(in) :: chemReg
     integer, intent(in)  :: nq ! total number of tracers in registry 
     logical, intent(out) :: isGOCART(nq)
     !                     -----
     isGOCART = .false.
     if ( chemReg%doing_O3 )  isGOCART(chemReg%i_O3 :chemReg%j_O3)  = .true.
     if ( chemReg%doing_CO )  isGOCART(chemReg%i_CO :chemReg%j_CO)  = .true.
     if ( chemReg%doing_CO2 ) isGOCART(chemReg%i_CO2:chemReg%j_CO2) = .true.
     if ( chemReg%doing_DU )  isGOCART(chemReg%i_DU :chemReg%j_DU)  = .true.
     if ( chemReg%doing_SS )  isGOCART(chemReg%i_SS :chemReg%j_SS)  = .true.
     if ( chemReg%doing_BC )  isGOCART(chemReg%i_BC :chemReg%j_BC)  = .true.
     if ( chemReg%doing_OC )  isGOCART(chemReg%i_OC :chemReg%j_OC)  = .true.
     if ( chemReg%doing_BRC ) isGOCART(chemReg%i_BRC :chemReg%j_BRC)  = .true.
     if ( chemReg%doing_SU )  isGOCART(chemReg%i_SU :chemReg%j_SU)  = .true.
     if ( chemReg%doing_CFC ) isGOCART(chemReg%i_CFC:chemReg%j_CFC) = .true.
     if ( chemReg%doing_Rn )  isGOCART(chemReg%i_Rn :chemReg%j_Rn)  = .true.
     if ( chemReg%doing_CH4 ) isGOCART(chemReg%i_CH4:chemReg%j_CH4) = .true.
     if ( chemReg%doing_NI )  isGOCART(chemReg%i_NI :chemReg%j_NI)  = .true.
   end subroutine Chem_RegistrySetIsGOCART

 end module Chem_RegistryMod

