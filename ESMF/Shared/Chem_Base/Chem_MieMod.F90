! $Id$

#include "MAPL_Exceptions.h"

!-------------------------------------------------------------------------
!      NASA/GSFC, Global Modeling & Assimilation Office, Code 900.3      !
!-------------------------------------------------------------------------
!BOP
!

! !MODULE:  Chem_MieMod --- Load and manipulate Mie tables
!
! !INTERFACE:
!

   module  Chem_MieMod

! !USES:

   use Chem_MieTableMod
   use Chem_RegistryMod
   use m_chars, only : uppercase
   use m_die, only: die
   use m_inpak90

#if defined(GEOS5)
   use ESMF
   use MAPL
#endif

   implicit none

! !PUBLIC TYPES:
!
   private
   public  Chem_Mie                ! Holds Mie Lookup Tables
                           
!
! !PUBLIC MEMBER FUNCTIONS:
!
   public  Chem_MieCreate          ! Constructor 
   public  Chem_MieDestroy         ! Destructor
   public  Chem_MieQuery           ! Query the Mie table to return parameters (qname interface)
   public  Chem_MieQueryTauList    
   public  Chem_MieQueryAllBand3D
   public  Chem_MieQueryAllBand4D
   public  Chem_MieQueryIdx        ! Query the index of the mie table given the qname

!
! !DESCRIPTION:
!
!  This module read the mie aerosol tables.
!
! !REVISION HISTORY:
!
!  23Mar2005 Colarco - Initial code.
!  11Jul2005 da Silva   Standardization.
!
!EOP
!-------------------------------------------------------------------------

! Mie LUT table
! Will be reduced from input files to the desired channels
! --------
  type Chem_Mie
!     private
     integer :: nch                               ! number of channels
     integer :: nMom=0                            ! number of moments (phase function)
     integer :: nPol=0                            ! number of moments (phase function)
     real, pointer    :: channels(:)              ! wavelengths

     character(len=255) :: rcfile
     character(len=255) :: du_optics_file
     character(len=255) :: ss_optics_file
     character(len=255) :: bc_optics_file
     character(len=255) :: oc_optics_file
     character(len=255) :: brc_optics_file
     character(len=255) :: su_optics_file
     character(len=255) :: ni_optics_file
     character(len=255) :: sm_optics_file

                                           ! mie tables -- dim(nch,nrh,nbin)
     type(Chem_MieTable), pointer :: mie_DU  => null()
     type(Chem_MieTable), pointer :: mie_SS  => null()
     type(Chem_MieTable), pointer :: mie_BC  => null()
     type(Chem_MieTable), pointer :: mie_OC  => null()
     type(Chem_MieTable), pointer :: mie_BRC => null()
     type(Chem_MieTable), pointer :: mie_SU  => null()
     type(Chem_MieTable), pointer :: mie_NI  => null()
     type(Chem_MieTable), pointer :: mie_SM  => null()

     integer :: nq                                ! number of tracers
     character(len=255), pointer  :: vname(:)  => null()
     integer, pointer             :: vindex(:) => null()
     type(Chem_MieTable), pointer :: vtable(:) => null()
                                          ! mapping of vtable for given idx
     type(Chem_MieTable), pointer :: vtableUse => null()   

  end type Chem_Mie

  interface Chem_MieCreate
     module procedure Chem_MieCreateFromCF
     module procedure Chem_MieCreateFromRC
  end interface
  
  interface Chem_MieQuery
     module procedure Chem_MieQueryByInt
     module procedure Chem_MieQueryByChar
  end interface


contains

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Chem_MieCreate --- Construct Mie LUTs from RC File
!
! !INTERFACE:
!

  Function Chem_MieCreateFromRC ( rcfile, rc, chemReg ) result(this)

  implicit none

! !INPUT PARAMETERS:

   character(len=*), intent(in)            :: rcfile  ! Mie table file name
   type(Chem_Registry), target, optional, intent(in) :: chemReg ! Optional chemReg

! !OUTPUT PARAMETERS:

   integer, intent(out) ::  rc            ! Error return code:
                                          !  0 - all is well
                                          !  1 - 

! !DESCRIPTION: Creates a Mie object. It reads Mie Tables and associates 
!               with a regulr Chem Registry describing aerosol mixing
!  ratios. Then chemReg is not specified, "Chem_MieRegistry.rc" or
!  "Chem_Registry.rc" are in, if available.
!
! !REVISION HISTORY:
!
!  09Mar2005 da Silva  API, prologues.
!  21Mar2010 da Silva  Added optional argument chemReg.
!
!EOP
!-------------------------------------------------------------------------

   character(len=*), parameter ::  myname = 'Chem_MieCreate'

   type(Chem_Mie) :: this
   type(Chem_Registry), pointer :: reg
   integer        :: ios, n, iq, iiq
   real, pointer  :: rh_table(:), lambda_table(:), &
                     bext(:,:,:), bsca(:,:,:), reff(:,:), gf(:,:), &
                     rhop(:,:), rhod(:)
   logical :: fexists
   character(len=255) :: name

   rc = 0

!  NOTE: when rc is mandatory it is not cool do call die(); in this case
!        the user should do the error trapping

!  Get the Chem Registry: optionally, uses a private name: Chem_MieRegistry

   if ( present(chemReg) ) then
      reg => ChemReg
   else
      allocate(reg)
      inquire ( file='Chem_MieRegistry.rc', exist=fexists )
      if ( fexists ) then
         reg = Chem_RegistryCreate(rc,rcfile='Chem_MieRegistry.rc')
         if ( rc /= 0 ) call die(myname, 'Cannot read Chem_MieRegistry.rc' )
      else
         reg = Chem_RegistryCreate(rc,rcfile='Chem_Registry.rc')
         if ( rc /= 0 ) call die(myname, 'Cannot read Chem_Registry.rc' )
      end if
   end if

!  Set up the hash table to map the Chem Registry to the
!  Mie tables
!  -----------------------------------------------------
   this%nq = reg%nq
   allocate(this%vname( this%nq) )
   allocate(this%vindex(this%nq) )
   allocate(this%vtable(this%nq) )
   do iq = 1, this%nq
    this%vindex(iq) = -1
    this%vname(iq)  = reg%vname(iq)
   enddo

   this%rcfile = rcfile

!   Load the resource file
!   ----------------------
    call i90_loadf ( rcfile, ios )
    if ( ios /= 0 ) call die(myname, 'could not read rc file '// &
         trim(rcfile) )


!   Set the number of channels to calculate over
!   --------------------------------------------
    call i90_label ( 'n_channels:', ios )
    if ( ios /= 0 ) then
     call die(myname, 'could not find channel number request')
    else
     this%nch = i90_gint ( ios )
     if ( ios /= 0 ) call die(myname,'could not parse number of channels')
    end if

!   Set the number of moments
!   -------------------------
    call i90_label ( 'n_moments:', ios )
    if ( ios /= 0 ) then
     this%nmom = 0
    else
     this%nmom = i90_gint ( ios )
     if ( ios /= 0 ) call die(myname,'could not parse number of moments')
    end if

!   Set the channels to calculate over
!   ----------------------------------
    allocate( this%channels(this%nch), stat = ios )
    call i90_label ( 'r_channels:', ios )
    if ( ios /= 0 ) then
     call die(myname, 'could not find channel number request')
    else
     do n = 1, this%nch
      this%channels(n) = i90_gfloat ( ios )
      if ( ios /= 0 ) call die(myname,'could not parse channels')
     enddo
    end if

!   Logic needs to be placed so that you check the mie tables against
!   the chem registry (bin size, species, etc.)  For now assume
!   they are all right.
!   -----------------------------------------------------------------
    call i90_label ( 'filename_optical_properties_DU:', ios )
     if ( ios /= 0 ) then
      call die(myname, 'could not parse DU filename label')
     else
      call i90_gtoken ( this%du_optics_file, ios )
      if ( ios /= 0 ) call die(myname,'could not parse DU filename')
     end if

    call i90_label ( 'filename_optical_properties_SS:', ios )
     if ( ios /= 0 ) then
      call die(myname, 'could not parse SS filename label')
     else
      call i90_gtoken ( this%ss_optics_file, ios )
      if ( ios /= 0 ) call die(myname,'could not parse SS filename')
     end if

    call i90_label ( 'filename_optical_properties_BC:', ios )
     if ( ios /= 0 ) then
      call die(myname, 'could not parse BC filename label')
     else
      call i90_gtoken ( this%bc_optics_file, ios )
      if ( ios /= 0 ) call die(myname,'could not parse BC filename')
     end if

    call i90_label ( 'filename_optical_properties_OC:', ios )
     if ( ios /= 0 ) then
      call die(myname, 'could not parse OC filename label')
     else
      call i90_gtoken ( this%oc_optics_file, ios )
      if ( ios /= 0 ) call die(myname,'could not parse OC filename')
     end if

    call i90_label ( 'filename_optical_properties_SU:', ios )
     if ( ios /= 0 ) then
      call die(myname, 'could not parse SU filename label')
     else
      call i90_gtoken ( this%su_optics_file, ios )
      if ( ios /= 0 ) call die(myname,'could not parse SU filename')
     end if

    call i90_label ( 'filename_optical_properties_NI:', ios )
     if ( ios /= 0 ) then
      call die(myname, 'could not parse NI filename label')
     else
      call i90_gtoken ( this%ni_optics_file, ios )
      if ( ios /= 0 ) call die(myname,'could not parse NI filename')
     end if

    call i90_label ( 'filename_optical_properties_BRC:', ios )
     if ( ios /= 0 ) then
      call die(myname, 'could not parse BRC filename label')
     else
      call i90_gtoken ( this%brc_optics_file, ios )
      if ( ios /= 0 ) call die(myname,'could not parse BRC filename')
     end if


!   Close resource file
    call I90_Release()


!  Allocate and fill Mie Table
!  ---------------------------
   allocate(this%mie_DU, this%mie_SS, this%mie_SU, stat = rc )
   if ( rc /= 0 ) return
   allocate(this%mie_BC, this%mie_OC, stat = rc )
   if ( rc /= 0 ) return
   allocate(this%mie_NI, this%mie_BRC, stat = rc )
   if ( rc /= 0 ) return
   this%mie_DU = Chem_MieTableCreate(this%du_optics_file, rc)
   if ( rc /= 0 ) call die(myname, 'could not create table for dust')
   this%mie_SS = Chem_MieTableCreate(this%ss_optics_file, rc)
   if ( rc /= 0 ) call die(myname, 'could not create table for sea salt')
   this%mie_SU = Chem_MieTableCreate(this%su_optics_file, rc)
   if ( rc /= 0 ) call die(myname, 'could not create table for sulfates')
   this%mie_OC = Chem_MieTableCreate(this%oc_optics_file, rc)
   if ( rc /= 0 ) call die(myname, 'could not create table for organic carbon')
   this%mie_BC = Chem_MieTableCreate(this%bc_optics_file, rc)
   if ( rc /= 0 ) call die(myname, 'could not create table for black carbon')
   this%mie_NI = Chem_MieTableCreate(this%ni_optics_file, rc)
   if ( rc /= 0 ) call die(myname, 'could not create table for nitrates')
   this%mie_BRC = Chem_MieTableCreate(this%brc_optics_file, rc)
   if ( rc /= 0 ) call die(myname, 'could not create table for brown carbon')

   call Chem_MieTableRead(this%mie_DU,this%nch,this%channels,rc,nmom=this%nmom)
   if ( rc /= 0 ) call die(myname, 'could not read table for dust')
   call Chem_MieTableRead(this%mie_SS,this%nch,this%channels,rc,nmom=this%nmom)
   if ( rc /= 0 ) call die(myname, 'could not read table for sea salt')
   call Chem_MieTableRead(this%mie_SU,this%nch,this%channels,rc,nmom=this%nmom)
   if ( rc /= 0 ) call die(myname, 'could not read table for sulfates')
   call Chem_MieTableRead(this%mie_OC,this%nch,this%channels,rc,nmom=this%nmom)
   if ( rc /= 0 ) call die(myname, 'could not read table for organic carbon')
   call Chem_MieTableRead(this%mie_BC,this%nch,this%channels,rc,nmom=this%nmom)
   if ( rc /= 0 ) call die(myname, 'could not read table for black carbon')
   call Chem_MieTableRead(this%mie_NI,this%nch,this%channels,rc,nmom=this%nmom)
   if ( rc /= 0 ) call die(myname, 'could not read table for nitrates')
   call Chem_MieTableRead(this%mie_BRC,this%nch,this%channels,rc,nmom=this%nmom)
   if ( rc /= 0 ) call die(myname, 'could not read table for brown carbon')

   this%nPol = this%mie_DU%nPol

!  Now map the mie tables to the hash table for the registry
!  This part is hard-coded for now!
!  ---------------------------------------------------------
   if(reg%doing_DU) then
    do iq = reg%i_DU, reg%j_DU
     this%vindex(iq) = iq-reg%i_DU + 1
     this%vtable(iq)  = this%mie_DU
    enddo
   endif
   if(reg%doing_SS) then
    do iq = reg%i_SS, reg%j_SS
     this%vindex(iq) = iq-reg%i_SS + 1
     this%vtable(iq)  = this%mie_SS
    enddo
   endif
   if(reg%doing_OC) then
    do iq = reg%i_OC, reg%j_OC
     this%vindex(iq) = iq-reg%i_OC + 1
     this%vtable(iq)  = this%mie_OC
    enddo
   endif
   if(reg%doing_BRC) then
    do iq = reg%i_BRC, reg%j_BRC
     this%vindex(iq) = iq-reg%i_BRC + 1
     this%vtable(iq)  = this%mie_BRC
    enddo
   endif
   if(reg%doing_BC) then
    do iq = reg%i_BC, reg%j_BC
     this%vindex(iq) = iq-reg%i_BC + 1
     this%vtable(iq)  = this%mie_BC
    enddo
   endif
   if(reg%doing_SU) then
    iiq = 0
    do iq = reg%i_SU, reg%j_SU
     name = trim(this%vname(iq))
!    Only sulfate aerosol species have entries in the Mie table
     if(name(1:3) == 'so4' .or. name(1:3) == 'SO4' .or. &
        name(1:3) == 'sul' .or. name(1:3) == 'SUL') then
        iiq = iiq + 1
        this%vindex(iq) = iiq
        this%vtable(iq) = this%mie_SU
     endif
    enddo
   endif
   if(reg%doing_NI) then
    iiq = 0
    do iq = reg%i_NI, reg%j_NI
     name = trim(this%vname(iq))
!    Only nitrate aerosol species have entries in the Mie table
     if(name(1:3) == 'no3' .or. name(1:3) == 'NO3') then
        iiq = iiq + 1
        this%vindex(iq) = iiq
        this%vtable(iq) = this%mie_NI
     endif
    enddo
   endif

!  All done
!  --------
   if ( .not. present(chemReg) ) then
      call Chem_RegistryDestroy(reg,rc)
      if ( rc /= 0 ) return
   end if

   return 

 end Function Chem_MieCreateFromRC


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Chem_MieCreate --- Construct Mie LUTs from CF object
!
! !INTERFACE:
!

  function Chem_MieCreateFromCF ( cf, rc, chemReg ) result(this)

#if !defined(GEOS5)

  integer, intent(in)  :: cf
  integer, intent(out) :: rc
  type(Chem_Mie) this

#else

! !INPUT PARAMETERS:

   type(ESMF_Config) :: cf  ! Mie table file name
   type(Chem_Registry), target, optional, intent(in) :: chemReg ! Optional chemReg

! !OUTPUT PARAMETERS:

   type(Chem_Mie) this
   integer, intent(out) ::  rc

! !DESCRIPTION:
!
!     This routine creates a LUT object from an ESMF configuration
!  attribute CF. This routine is usually called from GEOS-5.
!
!  IMPORTANT: Does not yet handle the phase function!!!!
!
! !REVISION HISTORY:
!
!  09Mar2005 da Silva  API, prologues.
!
!EOP
!-------------------------------------------------------------------------

   type(Chem_Registry), pointer :: reg
   integer        :: iq, rcs(32)
   integer        :: i, itick, iiq
   real, pointer  :: rh_table(:), lambda_table(:), &
                     bext(:,:,:), bsca(:,:,:), reff(:,:), gf(:,:), &
                     rhop(:,:), rhod(:)
   character(len=255) :: reg_filename, name

                      __Iam__('Chem_MieCreateFromCF')


!  We need a Chem Registry to map a variable name into 
!  the relevant Mie Table, mostly for efficiency reasons.
!  ------------------------------------------------------
   if ( present(chemReg) ) then

      reg => chemReg ! User supplied registry

!  Load the Chem Registry from Chem_MieRegistry.rc
!  -----------------------------------------------
   else

      allocate(reg,__STAT__)
      call ESMF_ConfigGetAttribute( CF, reg_filename, Label="CHEM_REGISTRY_FILENAME:" , &
                                    default='Chem_MieRegistry.rc', __RC__ )
      reg = Chem_RegistryCreate(rc,rcfile=reg_filename)
      if ( rc /= 0 ) return 

   end if

!  Set up the hash table to map the variable names to the
!  corresponding Mie Table
!  -----------------------------------------------------
   this%nq = reg%nq
   allocate(this%vname(this%nq), this%vindex(this%nq), stat=rc  )
   if ( rc /= 0 ) return 
   allocate(this%vtable(this%nq), stat=rc )
   if ( rc /= 0 ) return 
   do iq = 1, this%nq
    this%vindex(iq) = -1
    this%vname(iq)  = reg%vname(iq)
   enddo

!  Get file names for the optical tables
!  -------------------------------------
   call ESMF_ConfigGetAttribute( CF, this%du_optics_file, Label="DU_OPTICS:" , &
                                 default='ExtData/g5chem/x/opticsBands_DU.nc4', &
                                 __RC__ )
   call ESMF_ConfigGetAttribute( CF, this%ss_optics_file, Label="SS_OPTICS:" , &
                                 default='ExtData/g5chem/x/opticsBands_SS.nc4', &
                                 __RC__ )
   call ESMF_ConfigGetAttribute( CF, this%su_optics_file, Label="SU_OPTICS:" , &
                                 default='ExtData/g5chem/x/opticsBands_SU.nc4', &
                                 __RC__ )
   call ESMF_ConfigGetAttribute( CF, this%oc_optics_file, Label="OC_OPTICS:" , &
                                 default='ExtData/g5chem/x/opticsBands_OC.nc4', &
                                 __RC__ )
   call ESMF_ConfigGetAttribute( CF, this%bc_optics_file, Label="BC_OPTICS:" , &
                                 default='ExtData/g5chem/x/opticsBands_BC.nc4', &
                                 __RC__ )
   call ESMF_ConfigGetAttribute( CF, this%ni_optics_file, Label="NI_OPTICS:" , &
                                 default='ExtData/g5chem/x/opticsBands_NI.nc4', &
                                 __RC__ )
   call ESMF_ConfigGetAttribute( CF, this%brc_optics_file, Label="BRC_OPTICS:" , &
                                 default='ExtData/g5chem/x/opticsBands_BRC.nc4', &
                                 __RC__ )
   call ESMF_ConfigGetAttribute( CF, this%nch           , Label= "NUM_BANDS:" , &
                                 default=18, __RC__)

   allocate ( this%channels(this%nch), stat=rc )
   if ( rc /= 0 ) return 

   call ESMF_ConfigGetAttribute( CF, this%channels       , Label= "BANDS:" , &
                                 count=this%nch, rc=rc )

!  If there is no BAND definition on CF, make something up
!  -------------------------------------------------------
   if(rc /= ESMF_SUCCESS) then
      do i=1,this%nch
         this%channels(i) = i
      end do
   end if

   allocate(this%mie_DU, this%mie_SS, this%mie_SU, this%mie_BRC, &
            this%mie_BC, this%mie_OC, this%mie_NI, stat=rc)
   if ( rc /= 0 ) return 

   this%mie_DU  = Chem_MieTableCreate(this%du_optics_file, __RC__ )
   this%mie_SS  = Chem_MieTableCreate(this%ss_optics_file, __RC__ )
   this%mie_SU  = Chem_MieTableCreate(this%su_optics_file, __RC__ )
   this%mie_OC  = Chem_MieTableCreate(this%oc_optics_file, __RC__ )
   this%mie_BC  = Chem_MieTableCreate(this%bc_optics_file, __RC__ )
   this%mie_NI  = Chem_MieTableCreate(this%ni_optics_file, __RC__ )
   this%mie_BRC = Chem_MieTableCreate(this%brc_optics_file, __RC__ )

   call Chem_MieTableRead(this%mie_DU,this%nch,this%channels, __RC__)
   call Chem_MieTableRead(this%mie_SS,this%nch,this%channels, __RC__)
   call Chem_MieTableRead(this%mie_SU,this%nch,this%channels, __RC__)
   call Chem_MieTableRead(this%mie_OC,this%nch,this%channels, __RC__)
   call Chem_MieTableRead(this%mie_BC,this%nch,this%channels, __RC__)
   call Chem_MieTableRead(this%mie_NI,this%nch,this%channels, __RC__)
   call Chem_MieTableRead(this%mie_BRC,this%nch,this%channels, __RC__)

!  Now map the mie tables to the hash table for the registry
!  This part is hard-coded for now!
!  ---------------------------------------------------------
   if(reg%doing_DU) then
    do iq = reg%i_DU, reg%j_DU
     this%vname(iq)  = reg%vname(iq)
     this%vindex(iq) = iq-reg%i_DU + 1
     this%vtable(iq)  = this%mie_DU
    enddo
   endif
   if(reg%doing_SS) then
    do iq = reg%i_SS, reg%j_SS
     this%vindex(iq) = iq-reg%i_SS + 1
     this%vtable(iq)  = this%mie_SS
    enddo
   endif
   if(reg%doing_OC) then
    do iq = reg%i_OC, reg%j_OC
     this%vindex(iq) = iq-reg%i_OC + 1
     this%vtable(iq)  = this%mie_OC
    enddo
   endif
   if(reg%doing_BRC) then
    do iq = reg%i_BRC, reg%j_BRC
     this%vindex(iq) = iq-reg%i_BRC + 1
     this%vtable(iq)  = this%mie_BRC
    enddo
   endif
   if(reg%doing_BC) then
    do iq = reg%i_BC, reg%j_BC
     this%vindex(iq) = iq-reg%i_BC + 1
     this%vtable(iq)  = this%mie_BC
    enddo
   endif
   if(reg%doing_SU) then
    iiq = 0
    do iq = reg%i_SU, reg%j_SU
     name = trim(this%vname(iq))
!    Only sulfate aerosol species have entries in the Mie table
     if(name(1:3) == 'so4' .or. name(1:3) == 'SO4' .or. &
        name(1:3) == 'sul' .or. name(1:3) == 'SUL') then
        iiq = iiq + 1
        this%vindex(iq) = iiq
        this%vtable(iq) = this%mie_SU
     endif
    enddo
   endif
   if(reg%doing_NI) then
    iiq = 0
    do iq = reg%i_NI, reg%j_NI
     name = trim(this%vname(iq))
!    Only nitrate aerosol species have entries in the Mie table
     if(name(1:3) == 'no3' .or. name(1:3) == 'NO3') then
        iiq = iiq + 1
        this%vindex(iq) = iiq
        this%vtable(iq) = this%mie_NI
     endif
    enddo
   endif

!  All done
!  --------
   if ( .not. present(chemReg) ) then
      call Chem_RegistryDestroy(reg,__RC__)
      deallocate(reg)
   end if

#endif
   
   return

 end function Chem_MieCreateFromCF

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Chem_MieDestroy --- Destruct Mie Table
!
! !INTERFACE:
!
  subroutine Chem_MieDestroy ( this, rc )

! !USES:

! !INPUT/OUTPUT PARAMETERS:

   type(Chem_Mie), intent(inout) :: this

! !OUTPUT PARAMETERS:

  integer, intent(out) ::  rc              ! Error return code:
                                           !  0 - all is well
                                           !  1 - 


! !DESCRIPTION: Destructor for AOD object.
!
! !REVISION HISTORY:
!
!  23Mar2005 Colarco
!
!EOP
!-------------------------------------------------------------------------

   character(len=*), parameter :: Iam = 'Chem_MieDestroy'

   call Chem_MieTableDestroy(this%mie_DU, rc=rc)
   if ( rc /= 0 ) return
   call Chem_MieTableDestroy(this%mie_SS, rc=rc)
   if ( rc /= 0 ) return
   call Chem_MieTableDestroy(this%mie_SU, rc=rc)
   if ( rc /= 0 ) return
   call Chem_MieTableDestroy(this%mie_OC, rc=rc)
   if ( rc /= 0 ) return
   call Chem_MieTableDestroy(this%mie_BC, rc=rc)
   if ( rc /= 0 ) return
   call Chem_MieTableDestroy(this%mie_NI, rc=rc)
   if ( rc /= 0 ) return
   call Chem_MieTableDestroy(this%mie_BRC, rc=rc)
   if ( rc /= 0 ) return

   if ( associated(this%channels) )  deallocate(this%channels, stat=rc)
   if ( rc /= 0 ) return
   if ( associated(this%vname) )     deallocate(this%vname, stat=rc)
   if ( rc /= 0 ) return
   if ( associated(this%vindex) )    deallocate(this%vindex, stat=rc)
   if ( rc /= 0 ) return
   if ( associated(this%vtable) )    deallocate(this%vtable, stat=rc)
   if ( rc /= 0 ) return
   if ( associated(this%mie_DU) )    deallocate(this%mie_DU, stat=rc)
   if ( rc /= 0 ) return
   if ( associated(this%mie_SS) )    deallocate(this%mie_SS, stat=rc)
   if ( rc /= 0 ) return
   if ( associated(this%mie_SS) )    deallocate(this%mie_SS, stat=rc)
   if ( rc /= 0 ) return
   if ( associated(this%mie_OC) )    deallocate(this%mie_OC, stat=rc)
   if ( rc /= 0 ) return
   if ( associated(this%mie_BC) )    deallocate(this%mie_BC, stat=rc)
   if ( rc /= 0 ) return
   if ( associated(this%mie_NI) )    deallocate(this%mie_NI, stat=rc)
   if ( rc /= 0 ) return
   if ( associated(this%mie_BRC) )   deallocate(this%mie_BRC, stat=rc)
   if ( rc /= 0 ) return

end subroutine Chem_MieDestroy 

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Chem_MieQueryIdx --- Return the index of the mie table given
!                                  a qname requested
!
!
! !INTERFACE:
!
   Function Chem_MieQueryIdx ( this, qname, rc ) result(idx)

   implicit none

! !INPUT PARAMETERS:

   type(Chem_Mie), intent(inout) :: this   ! Input mie table structure
   character(len=*), intent(in)  :: qname  ! Variable name to find in table, e.g., du001

! !OUTPUT PARAMETERS:

   integer, optional, intent(out) ::  rc ! Error return code:
                                         !  0 - all is well
                                         !  1 - 
! !DESCRIPTION:
!
!
! !REVISION HISTORY:
!
!   24Apr2006, PRC
!
!EOP
!-------------------------------------------------------------------------
      character(len=255) :: NAME
      integer            :: idx         ! Index number in Mie table of qname
      integer            :: iq, i

!     Find the right table for this aerosol from its name

      NAME = trim(qname)

!     Remove qualifier from variable name: GOCART::du001 --> du001
!     ------------------------------------------------------------
      i = index(NAME,'::')
      if ( i > 0 ) then
         NAME = NAME(i+2:)
      end if

      idx = -1
      do iq = 1, this%nq
       if(uppercase(trim((NAME))) .eq. uppercase(trim(this%vname(iq)))) then
        idx = this%vindex(iq)
        this%vtableUse => this%vtable(iq)
        exit
       endif
      enddo

      if(present(rc)) then
         if(idx .eq. -1) then
            rc = 1
         else
            rc = 0
         end if
      end if

      return

  end Function Chem_MieQueryIdx

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Chem_MieQuery --- Return Tau, SSA, etc (scalar version)
!
!
! !INTERFACE:
!
   subroutine Chem_MieQueryByInt ( this, idx, channel, q_mass, rh,     &
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
   real,    optional,      intent(out) :: pmom(:,:)
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

      character(len=*), parameter  :: Iam = 'Chem_MieQueryByInt'

      if ( present(rc) ) rc = 0

      ICHANNEL = nint(CHANNEL)
      TABLE => this%vtableUse
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

      if(present(pmom)) then
         pmom(:,:) = TABLE%pmom(irh  ,ichannel,TYPE,:,:) * (1.-arh) &
                   + TABLE%pmom(irhp1,ichannel,TYPE,:,:) * arh
      endif

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

      return

 end subroutine Chem_MieQueryByInt


   subroutine Chem_MieQueryTauList ( this, idx, channel, q_mass, rh, tau, rc )


   type(Chem_Mie), target, intent(in ) :: this     
   integer,                intent(in ) :: idx     ! variable index on Chem_Mie
   real,                   intent(in ) :: channel ! channel number
   real,                   intent(in ) :: q_mass(:)  ! aerosol mass [kg/m2],
   real,                   intent(in ) :: rh(:)      ! relative himidity
   real,                   intent(out) :: tau(:)   ! aerol optical depth
   integer, optional,      intent(out) :: rc    ! error code

!-------------------------------------------------------------------------

      integer                      :: ICHANNEL, i
      integer                      :: irh, irhp1, isnap
      real                         :: arh
      type(Chem_MieTable), pointer :: TABLE

      character(len=*), parameter  :: Iam = 'Chem_MieQueryList'

      if ( present(rc) ) rc = 0

      ICHANNEL = nint(CHANNEL)
      TABLE => this%vtableUse

!     Now map the input RH to the high resolution hash table for RH

      do i=1,size(tau)
         arh = rh(i)
         if(arh > .99) arh = .99
         if(arh < 0.0) arh = 0.0

         isnap = int((arh+0.001)*1000.)
         if(isnap .lt. 1) isnap = 1

         arh   = TABLE%rha( isnap )
         irh   = TABLE%rhi( isnap )
         irhp1 = irh+1

         if(irhp1 .gt. TABLE%nrh) irhp1 = TABLE%nrh

         tau(i) = (  TABLE%bext(irh  ,ichannel,idx ) * (1.-arh) &
                  +  TABLE%bext(irhp1,ichannel,idx ) *     arh  )*q_mass(i)
      enddo

!  All Done
!----------

      return

    end subroutine Chem_MieQueryTauList


   subroutine Chem_MieQueryByChar( this, idx, channel, q_mass, rh,     &
                                   tau, ssa, gasym, bext, bsca, bbck,  &
                                   rEff, pmom, p11, p22, rc )

!  ! INPUT parameters
   type(Chem_Mie), target, intent(in ) :: this     
   character(*),           intent(in ) :: idx     ! variable index on Chem_Mie
   real,                   intent(in ) :: channel ! channel number
   real,                   intent(in ) :: q_mass  ! aerosol mass [kg/m2],
   real,                   intent(in ) :: rh      ! relative himidity

!  ! OUTPUT Parameters
   real,    optional,      intent(out) :: tau   ! aerol extinction optical depth
   real,    optional,      intent(out) :: ssa   ! single scattering albedo
   real,    optional,      intent(out) :: gasym ! asymmetry parameter
   real,    optional,      intent(out) :: bext  ! mass extinction efficiency [m2 (kg dry mass)-1]
   real,    optional,      intent(out) :: bsca  ! mass scattering efficiency [m2 (kg dry mass)-1]
   real,    optional,      intent(out) :: bbck  ! mass backscatter efficiency [m2 (kg dry mass)-1]
   real,    optional,      intent(out) :: reff  ! effective radius (micron)
   real,    optional,      intent(out) :: pmom(:,:)
   real,    optional,      intent(out) :: p11   ! P11 phase function at backscatter
   real,    optional,      intent(out) :: p22   ! P22 phase function at backscatter
   integer, optional,      intent(out) :: rc    ! error code

   integer :: iq, i

   character(len=*), parameter  :: Iam = 'Chem_MieQueryByChar'
   character(len=255) :: NAME

   if ( present(rc) ) rc = 0

!  Remove qualifier from variable name: GOCART::du001 --> du001
!   ------------------------------------------------------------
   NAME = trim(idx)
   i = index(NAME,'::')
   if ( i > 0 ) then
      NAME = NAME(i+2:)
   end if

   do iq = 1, this%nq
      if( uppercase(trim(NAME)) == uppercase(trim(this%vname(iq)))) then
         call  Chem_MieQueryByInt( this, iq, channel, q_mass, rh,     &
                             tau, ssa, gasym, bext, bsca, bbck, &
                             rEff, pmom, p11, p22, rc=rc )
         if ( rc /= 0 ) return
      endif
   enddo

 end subroutine Chem_MieQueryByChar

   subroutine Chem_MieQueryAllBand3D  ( this, idx, nbands, offset, q_mass, rh, &
                                        tau, ssa, asy, rc )


   type(Chem_Mie), target, intent(in ) :: this     
   integer,                intent(in ) :: idx     ! variable index on Chem_Mie
   integer,                intent(in ) :: nbands  ! number of bands
   integer,                intent(in ) :: offset  ! offset of bands
   real,                   intent(in ) :: q_mass(:,:)  ! aerosol mass [kg/m2],
   real,                   intent(in ) :: rh(:,:)      ! relative humidity
   real,                   intent(out) :: tau(:,:,:)   ! aerosol optical depth
   real,                   intent(out) :: ssa(:,:,:)   ! single scattering albedo
   real,                   intent(out) :: asy(:,:,:)   ! asymmetry parameter
   integer, optional,      intent(out) :: rc    ! error code

!-------------------------------------------------------------------------

      integer                      :: ICHANNEL, i, k, STATUS
      integer                      :: irh, irhp1, isnap
      real                         :: arh
      type(Chem_MieTable), pointer :: TABLE
      integer                      :: II, KK
      real                         :: bextIn, bscaIn, gasymIn

      real, allocatable            :: bext_band(:)
      real, allocatable            :: bsca_band(:)
      real, allocatable            :: gasy_band(:)

      character(len=*), parameter  :: Iam = 'Chem_MieQueryAllBand3D'

      if ( present(rc) ) rc = 0

      TABLE => this%vtableUse

      II = size(rh,1)
      KK = size(rh,2)

      allocate(bext_band(TABLE%nrh), source=0.0, __STAT__)
      allocate(bsca_band(TABLE%nrh), source=0.0, __STAT__)
      allocate(gasy_band(TABLE%nrh), source=0.0, __STAT__)

!     Now map the input RH to the high resolution hash table for RH

      do ichannel=1,nbands

         bext_band = TABLE%bext(:,ichannel+offset,idx)
         bsca_band = TABLE%bsca(:,ichannel+offset,idx)
         gasy_band =    TABLE%g(:,ichannel+offset,idx)

         do k=1,KK
            do i=1,II

               arh = rh(i,k)
               arh = max(arh,0.0)
               arh = min(arh,0.99)
      
               isnap = int((arh+0.001)*1000.)
               if(isnap .lt. 1) isnap = 1
      
               arh   = TABLE%rha( isnap )
               irh   = TABLE%rhi( isnap )
               irhp1 = irh+1
      
               if(irhp1 .gt. TABLE%nrh) irhp1 = TABLE%nrh

               bextIn  = bext_band(irh  ) * (1.-arh) &
                       + bext_band(irhp1) *     arh
      
               bscaIn  = bsca_band(irh  ) * (1.-arh) &
                       + bsca_band(irhp1) *     arh
         
               gasymIn = gasy_band(irh  ) * (1.-arh) &
                       + gasy_band(irhp1) *     arh
      
!              Fill the requested outputs

               tau(i,k,ichannel) = bextIn * q_mass(i,k)
               ssa(i,k,ichannel) = bscaIn/bextIn
               asy(i,k,ichannel) = gasymIn

            enddo
         enddo
      enddo

      deallocate(bext_band)
      deallocate(bsca_band)
      deallocate(gasy_band)

!  All Done
!----------

      return

    end subroutine Chem_MieQueryAllBand3D

   subroutine Chem_MieQueryAllBand4D  ( this, idx, nbands, offset, q_mass, rh, &
                                        tau, ssa, asy, rc )


   type(Chem_Mie), target, intent(in ) :: this     
   integer,                intent(in ) :: idx     ! variable index on Chem_Mie
   integer,                intent(in ) :: nbands  ! number of bands
   integer,                intent(in ) :: offset  ! offset of bands
   real,                   intent(in ) :: q_mass(:,:,:)  ! aerosol mass [kg/m2],
   real,                   intent(in ) :: rh(:,:,:)      ! relative humidity
   real,                   intent(out) :: tau(:,:,:,:)   ! aerosol optical depth
   real,                   intent(out) :: ssa(:,:,:,:)   ! single scattering albedo
   real,                   intent(out) :: asy(:,:,:,:)   ! asymmetry parameter
   integer, optional,      intent(out) :: rc    ! error code

!-------------------------------------------------------------------------

      integer                      :: ICHANNEL, i, j, k, STATUS
      integer                      :: irh, irhp1, isnap
      real                         :: arh
      type(Chem_MieTable), pointer :: TABLE
      integer                      :: II, JJ, KK
      real                         :: bextIn, bscaIn, gasymIn

      real, allocatable            :: bext_band(:)
      real, allocatable            :: bsca_band(:)
      real, allocatable            :: gasy_band(:)

      character(len=*), parameter  :: Iam = 'Chem_MieQueryAllBand4D'

      if ( present(rc) ) rc = 0

      TABLE => this%vtableUse

      II = size(rh,1)
      JJ = size(rh,2)
      KK = size(rh,3)

      allocate(bext_band(TABLE%nrh), source=0.0, __STAT__)
      allocate(bsca_band(TABLE%nrh), source=0.0, __STAT__)
      allocate(gasy_band(TABLE%nrh), source=0.0, __STAT__)

!     Now map the input RH to the high resolution hash table for RH

      do ichannel=1,nbands

         bext_band = TABLE%bext(:,ichannel+offset,idx)
         bsca_band = TABLE%bsca(:,ichannel+offset,idx)
         gasy_band =    TABLE%g(:,ichannel+offset,idx)

         do k=1,KK
            do j=1,JJ
               do i=1,II

                  arh = rh(i,j,k)
                  arh = max(arh,0.0)
                  arh = min(arh,0.99)
         
                  isnap = int((arh+0.001)*1000.)
                  if(isnap .lt. 1) isnap = 1
         
                  arh   = TABLE%rha( isnap )
                  irh   = TABLE%rhi( isnap )
                  irhp1 = irh+1
         
                  if(irhp1 .gt. TABLE%nrh) irhp1 = TABLE%nrh

                  bextIn  = bext_band(irh  ) * (1.-arh) &
                          + bext_band(irhp1) *     arh
         
                  bscaIn  = bsca_band(irh  ) * (1.-arh) &
                          + bsca_band(irhp1) *     arh
            
                  gasymIn = gasy_band(irh  ) * (1.-arh) &
                          + gasy_band(irhp1) *     arh
         
!                 Fill the requested outputs

                  tau(i,j,k,ichannel) = bextIn * q_mass(i,j,k)
                  ssa(i,j,k,ichannel) = bscaIn/bextIn
                  asy(i,j,k,ichannel) = gasymIn

               enddo
            enddo
         enddo
      enddo

      deallocate(bext_band)
      deallocate(bsca_band)
      deallocate(gasy_band)

!  All Done
!----------

      return

    end subroutine Chem_MieQueryAllBand4D

 end module Chem_MieMod

