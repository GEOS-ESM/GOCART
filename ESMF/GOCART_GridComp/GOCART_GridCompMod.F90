#include "MAPL_Generic.h"

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1     !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: GOCART_GridCompMod - The GOCART Aerosol Grid Component
!
! !INTERFACE:
!
   Module GOCART_GridCompMod
!
! !USES:
!
   use ESMF
   use MAPL

   use Chem_Mod              ! Chemistry Base Class
   use Chem_UtilMod, only: Chem_UtilNegFiller
   use Aero_GridCompMod      ! Parent Aerosol component with IRF methods but no SetServices()

   use ConvectionMod, only: Disable_Convection
   
   use mod_network , only: network_type
   use mod_kinds, only: ik, rk


   implicit none
   private

   type(Chem_Mie), dimension(2), save :: gocartMieTable
   integer, parameter :: instanceComputational = 1
   integer, parameter :: instanceData          = 2

   character(len=*), parameter :: H2O2_RECYCLE_ALARM = 'GOCART::RECYCLE_H2O2'
   character(len=*), parameter :: HNO3_RECYCLE_ALARM = 'GOCART::RECYCLE_HNO3'

!
! !PUBLIC MEMBER FUNCTIONS:

   public SetServices
!
! !DESCRIPTION: 
!
!   {\tt GOCART} is a gridded component from the GOCART model and includes 
!  dust, sea salt, sulfates, organic and black carbon. In addition, we
!  also include closely related components for CO and CO2 with relatively
!  simple parameterization of the chemical processes, but sharing
!  consistent emissions with the aerosols.
!
!  This code derives from the pre-ESMF Chem component from GEOS-4. This
!  GEOS-4 Chem "component" used ESMF like constructs (Chem component class, 
!  import/export states, etc) but no ESMF specific data types because of 
!  an odd incompatibility with the fvGCM code (the so-called 
!  {\tt oldworld} library. Unlike GEOS-4, the Stratospheric Chemistry
!  component is treated separately here.
!
! !REVISION HISTORY:
!
!  25feb2005  da Silva  First crack.
!  19jul2006  da Silva  First separate GOCART component.
!
!EOP
!-------------------------------------------------------------------------

  type GOCART_State
     private
     type(Chem_Registry), pointer :: chemReg => null()
     type(Aero_GridComp), pointer :: gcChem  => null()
     type(Chem_Bundle), pointer   :: w_c     => null()
     logical                      :: data_driven = .false.
   end type GOCART_State

  type GOCART_WRAP
     type (GOCART_State), pointer :: PTR => null()
  end type GOCART_WRAP

CONTAINS


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: SetServices --- Sets IRF services for GOCART Grid Component
!
! !INTERFACE:

   subroutine SetServices ( GC, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, optional                  :: RC  ! return code

! !DESCRIPTION: Sets Initialize, Run and Finalize services. 
!
! !REVISION HISTORY:
!
!  25feb2005  da Silva  First crack.
!
!EOP
!-------------------------------------------------------------------------


!   ErrLog Variables
!   ----------------
    character(len=ESMF_MAXSTR)      :: IAm = 'SetServices'
    integer                         :: STATUS
    character(len=ESMF_MAXSTR)      :: COMP_NAME

!   Local derived type aliases
!   --------------------------
    type (ESMF_Config)              :: CF
    type (GOCART_State), pointer    :: state   ! internal, that is
    type (GOCART_wrap)              :: wrap
    type(Chem_Registry), pointer    :: r

    integer                         :: n, nq
    integer                         :: DO_CO2CNNEE
    character(len=ESMF_MAXSTR)      :: FRIENDLIES
    character(len=ESMF_MAXSTR)      :: AEROFRIENDLY
    character(len=ESMF_MAXSTR)      :: providerName
    character(len=ESMF_MAXSTR)      :: short_name
    real                            :: DEFVAL
    real                            :: DEFVAL_CO2

    character(len=ESMF_MAXSTR)      :: field_name
    character(len=ESMF_MAXSTR)      :: chem_registry_file

!                              ------------

    
!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, CONFIG=CF, __RC__ )

    Iam = TRIM(COMP_NAME) // '::' // 'SetServices'

!   Wrap internal state for storing in GC; rename legacyState
!   -------------------------------------
    allocate ( state, __STAT__ )
    wrap%ptr => state

!   Is the component data driven
!   ----------------------------
    state%data_driven = IsDataDrivenGC_(GC, __RC__)

!   Start by loading the Chem Registry
!   ----------------------------------
    allocate ( state%chemReg )

    if (state%data_driven) then
        state%chemReg = Chem_RegistryCreate(STATUS, rcfile='GOCARTdata_AerRegistry.rc')
        VERIFY_(STATUS)
    else
       call ESMF_ConfigGetAttribute(cf, chem_registry_file, label = "Chem_Registry_File:", &
            default = "Chem_Registry.rc", rc = status)
       VERIFY_(status)
       state%chemReg = Chem_RegistryCreate(STATUS, rcfile=chem_registry_file)
        VERIFY_(STATUS)
    end if    

    r => state%chemReg   ! short hand
    

!                       ------------------------
!                       ESMF Functional Services
!                       ------------------------

!   Set the Initialize, Run, Finalize entry points
!   ----------------------------------------------
     if ( r%doing_GOCART ) then

        if(MAPL_AM_I_ROOT()) then
         print *, trim(Iam)//': ACTIVE'
         print *,' '
        end if

        call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_INITIALIZE,  Initialize_, __RC__ )
        call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,         Run1_,       __RC__ )
        call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,         Run2_,       __RC__ )

        call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_FINALIZE,    Finalize_,   __RC__ )
        
!       Store internal state in GC
!       --------------------------
        call ESMF_UserCompSetInternalState ( GC, 'GOCART_state', wrap, STATUS )
        VERIFY_(STATUS)

     else

        if (MAPL_AM_I_ROOT()) then
         print *, trim(Iam)//': NOT ACTIVE, defaulting to Generic No-op stubs'
        end if 

        call MAPL_GenericSetServices ( GC, __RC__ )
        RETURN_(ESMF_SUCCESS)

     endif


!                         ------------------
!                         GEOS Data Services
!                         ------------------

!  NOTE: For now, always define import state to avoid breaking connectivities.

!!BOS
!
! !IMPORT STATE:

GOCARTdata_IMPORTS: if (state%data_driven) then

!   Pressure at layer edges
!   -----------------------
    call MAPL_AddImportSpec(GC,                            &
       SHORT_NAME = 'PLE',                                 &
       LONG_NAME  = 'air_pressure',                        &
       UNITS      = 'Pa',                                  &
       DIMS       = MAPL_DimsHorzVert,                     &
       VLOCATION  = MAPL_VLocationEdge,                    &
       RESTART    = MAPL_RestartSkip,     __RC__)

!   Pressure thickness
!   ------------------
    call MAPL_AddImportSpec(GC,                            &
       SHORT_NAME = 'DELP',                                &
       LONG_NAME  = 'pressure_thickness',                  &
       UNITS      = 'Pa',                                  &
       DIMS       = MAPL_DimsHorzVert,                     &
       VLOCATION  = MAPL_VLocationCenter,                  &
       RESTART    = MAPL_RestartSkip,     __RC__)

!   RH: is between 0 and 1
!   ----------------------
    call MAPL_AddImportSpec(GC,                            &
       SHORT_NAME = 'RH2',                                 &
       LONG_NAME  = 'Rel_Hum_after_moist',                 &
       UNITS      = '1',                                   &
       DIMS       = MAPL_DimsHorzVert,                     &
       VLOCATION  = MAPL_VLocationCenter,                  &
       RESTART    = MAPL_RestartSkip,     __RC__)

!   T
!   -
    call MAPL_AddImportSpec(GC,                            &
       SHORT_NAME = 'T',                                   &
       LONG_NAME  = 'air_temperature',                     &
       UNITS      = 'K',                                   &
       DIMS       = MAPL_DimsHorzVert,                     &
       VLOCATION  = MAPL_VLocationCenter,                  &
       RESTART    = MAPL_RestartSkip,     __RC__)

    nq = r%nq     ! total number of chemical tracers
       
!   Loop over all constituents on registry
!   --------------------------------------
    do n = r%i_GOCART, r%j_GOCART
             
!       3D mass mixing ratios
!       ---------------------
        call MAPL_AddImportSpec(GC,                        &
           SHORT_NAME = 'clim'//trim(r%vname(n)),          &
           LONG_NAME  = r%vtitle(n),                       &
           UNITS      = r%vunits(n),                       &
           DIMS       = MAPL_DimsHorzVert,                 &
           VLOCATION  = MAPL_VLocationCenter,              &
           RESTART    = MAPL_RestartSkip, __RC__)
    end do


!   2D deposition fluxes
!   --------------------
    IMPORT_DUST_DEP_FLUXES: if (r%doing_DU) then
        do n = 1, r%n_DU
            ! dry deposition            
            write (field_name, '(A, I0.3)') 'DUDP', n

            call MAPL_AddImportSpec(GC,                    &
               SHORT_NAME = 'clim'//field_name,            &
               LONG_NAME  = r%vtitle(n),                   &
               UNITS      = r%vunits(n),                   &
               DIMS       = MAPL_DimsHorzOnly,             &
               VLOCATION  = MAPL_VLocationCenter,          &
               RESTART    = MAPL_RestartSkip, __RC__)
    
            ! wet deposition    
            write (field_name, '(A, I0.3)') 'DUWT', n

            call MAPL_AddImportSpec(GC,                    &
               SHORT_NAME = 'clim'//field_name,            &
               LONG_NAME  = r%vtitle(n),                   &
               UNITS      = r%vunits(n),                   &
               DIMS       = MAPL_DimsHorzOnly,             &
               VLOCATION  = MAPL_VLocationCenter,          &
               RESTART    = MAPL_RestartSkip, __RC__)

            ! gravitational settling
            write (field_name, '(A, I0.3)') 'DUSD', n

            call MAPL_AddImportSpec(GC,                    &
               SHORT_NAME = 'clim'//field_name,            &
               LONG_NAME  = r%vtitle(n),                   &
               UNITS      = r%vunits(n),                   &
               DIMS       = MAPL_DimsHorzOnly,             &
               VLOCATION  = MAPL_VLocationCenter,          &  
               RESTART    = MAPL_RestartSkip, __RC__)

            ! convective scavenging
            write (field_name, '(A, I0.3)') 'DUSV', n

            call MAPL_AddImportSpec(GC,                    &
               SHORT_NAME = 'clim'//field_name,            &
               LONG_NAME  = r%vtitle(n),                   &
               UNITS      = r%vunits(n),                   &
               DIMS       = MAPL_DimsHorzOnly,             &
               VLOCATION  = MAPL_VLocationCenter,          &  
               RESTART    = MAPL_RestartSkip, __RC__)
        end do
    end if IMPORT_DUST_DEP_FLUXES


!   Black Carbon
!   ------------
    IMPORT_BC_DEP_FLUXES: if (r%doing_BC) then
        do n = 1, r%n_BC
            ! dry deposition
            write (field_name, '(A, I0.3)') 'BCDP', n

            call MAPL_AddImportSpec(GC,                    &
               SHORT_NAME = 'clim'//field_name,            &
               LONG_NAME  = r%vtitle(n),                   &
               UNITS      = r%vunits(n),                   &
               DIMS       = MAPL_DimsHorzOnly,             &
               VLOCATION  = MAPL_VLocationCenter,          &
               RESTART    = MAPL_RestartSkip, __RC__)

            ! wet deposition    
            write (field_name, '(A, I0.3)') 'BCWT', n

            call MAPL_AddImportSpec(GC,                    &
               SHORT_NAME  = 'clim'//field_name,           &
               LONG_NAME   = r%vtitle(n),                  &
               UNITS       = r%vunits(n),                  &
               DIMS        = MAPL_DimsHorzOnly,            &
               VLOCATION   = MAPL_VLocationCenter,         &
               RESTART     = MAPL_RestartSkip,   __RC__)

            ! gravitational settling
            write (field_name, '(A, I0.3)') 'BCSD', n

            call MAPL_AddImportSpec(GC,                    &
               SHORT_NAME = 'clim'//field_name,            &
               LONG_NAME  = r%vtitle(n),                   &
               UNITS      = r%vunits(n),                   &
               DIMS       = MAPL_DimsHorzOnly,             &
               VLOCATION  = MAPL_VLocationCenter,          &  
               RESTART    = MAPL_RestartSkip, __RC__)

            ! convective scavenging
            write (field_name, '(A, I0.3)') 'BCSV', n

            call MAPL_AddImportSpec(GC,                    &
               SHORT_NAME = 'clim'//field_name,            &
               LONG_NAME  = r%vtitle(n),                   &
               UNITS      = r%vunits(n),                   &
               DIMS       = MAPL_DimsHorzOnly,             &
               VLOCATION  = MAPL_VLocationCenter,          &  
               RESTART    = MAPL_RestartSkip, __RC__)
        end do
    end if IMPORT_BC_DEP_FLUXES

!   Organic Carbon
!   --------------
    IMPORT_OC_DEP_FLUXES: if (r%doing_OC) then
        do n = 1, r%n_OC
            ! dry deposition            
            write (field_name, '(A, I0.3)') 'OCDP', n

            call MAPL_AddImportSpec(GC,                    &
               SHORT_NAME = 'clim'//field_name,            &
               LONG_NAME  = r%vtitle(n),                   &
               UNITS      = r%vunits(n),                   &
               DIMS       = MAPL_DimsHorzOnly,             &
               VLOCATION  = MAPL_VLocationCenter,          &
               RESTART    = MAPL_RestartSkip,   __RC__)

            ! wet deposition    
            write (field_name, '(A, I0.3)') 'OCWT', n

            call MAPL_AddImportSpec(GC,                    &
               SHORT_NAME = 'clim'//field_name,            &
               LONG_NAME  = r%vtitle(n),                   &
               UNITS      = r%vunits(n),                   &
               DIMS       = MAPL_DimsHorzOnly,             &
               VLOCATION  = MAPL_VLocationCenter,          &
               RESTART    = MAPL_RestartSkip,   __RC__)

            ! gravitational settling
            write (field_name, '(A, I0.3)') 'OCSD', n

            call MAPL_AddImportSpec(GC,                    &
               SHORT_NAME = 'clim'//field_name,            &
               LONG_NAME  = r%vtitle(n),                   &
               UNITS      = r%vunits(n),                   &
               DIMS       = MAPL_DimsHorzOnly,             &
               VLOCATION  = MAPL_VLocationCenter,          &  
               RESTART    = MAPL_RestartSkip, __RC__)

            ! convective scavenging
            write (field_name, '(A, I0.3)') 'OCSV', n

            call MAPL_AddImportSpec(GC,                    &
               SHORT_NAME = 'clim'//field_name,            &
               LONG_NAME  = r%vtitle(n),                   &
               UNITS      = r%vunits(n),                   &
               DIMS       = MAPL_DimsHorzOnly,             &
               VLOCATION  = MAPL_VLocationCenter,          &  
               RESTART    = MAPL_RestartSkip, __RC__)
        end do
    end if IMPORT_OC_DEP_FLUXES

!   Sulfate
!   --------
    IMPORT_SU_DEP_FLUXES: if (r%doing_SU) then
        do n = 1, r%n_SU
            ! dry deposition
            write (field_name, '(A, I0.3)') 'SUDP', n

            call MAPL_AddImportSpec(GC,                    &
               SHORT_NAME = 'clim'//field_name,            &
               LONG_NAME  = r%vtitle(n),                   &
               UNITS      = r%vunits(n),                   &
               DIMS       = MAPL_DimsHorzOnly,             &
               VLOCATION  = MAPL_VLocationCenter,          &
               RESTART    = MAPL_RestartSkip, __RC__)

            ! wet deposition    
            write (field_name, '(A, I0.3)') 'SUWT', n

            call MAPL_AddImportSpec(GC,                    &
               SHORT_NAME  = 'clim'//field_name,           &
               LONG_NAME   = r%vtitle(n),                  &
               UNITS       = r%vunits(n),                  &
               DIMS        = MAPL_DimsHorzOnly,            &
               VLOCATION   = MAPL_VLocationCenter,         &
               RESTART     = MAPL_RestartSkip,   __RC__)

            ! gravitational settling
            write (field_name, '(A, I0.3)') 'SUSD', n

            call MAPL_AddImportSpec(GC,                    &
               SHORT_NAME = 'clim'//field_name,            &
               LONG_NAME  = r%vtitle(n),                   &
               UNITS      = r%vunits(n),                   &
               DIMS       = MAPL_DimsHorzOnly,             &
               VLOCATION  = MAPL_VLocationCenter,          &  
               RESTART    = MAPL_RestartSkip, __RC__)

            ! convective scavenging
            write (field_name, '(A, I0.3)') 'SUSV', n

            call MAPL_AddImportSpec(GC,                    &
               SHORT_NAME = 'clim'//field_name,            &
               LONG_NAME  = r%vtitle(n),                   &
               UNITS      = r%vunits(n),                   &
               DIMS       = MAPL_DimsHorzOnly,             &
               VLOCATION  = MAPL_VLocationCenter,          &  
               RESTART    = MAPL_RestartSkip, __RC__)
        end do
    end if IMPORT_SU_DEP_FLUXES

!   Sea Salt
!   --------
    IMPORT_SS_DEP_FLUXES: if (r%doing_SS) then
        do n = 1, r%n_SS
            ! dry deposition
            write (field_name, '(A, I0.3)') 'SSDP', n

            call MAPL_AddImportSpec(GC,                    &
               SHORT_NAME = 'clim'//field_name,            &
               LONG_NAME  = r%vtitle(n),                   &
               UNITS      = r%vunits(n),                   &
               DIMS       = MAPL_DimsHorzOnly,             &
               VLOCATION  = MAPL_VLocationCenter,          &
               RESTART    = MAPL_RestartSkip, __RC__)

            ! wet deposition    
            write (field_name, '(A, I0.3)') 'SSWT', n

            call MAPL_AddImportSpec(GC,                    &
               SHORT_NAME  = 'clim'//field_name,           &
               LONG_NAME   = r%vtitle(n),                  &
               UNITS       = r%vunits(n),                  &
               DIMS        = MAPL_DimsHorzOnly,            &
               VLOCATION   = MAPL_VLocationCenter,         &
               RESTART     = MAPL_RestartSkip,   __RC__)

            ! gravitational settling
            write (field_name, '(A, I0.3)') 'SSSD', n

            call MAPL_AddImportSpec(GC,                    &
               SHORT_NAME = 'clim'//field_name,            &
               LONG_NAME  = r%vtitle(n),                   &
               UNITS      = r%vunits(n),                   &
               DIMS       = MAPL_DimsHorzOnly,             &
               VLOCATION  = MAPL_VLocationCenter,          &  
               RESTART    = MAPL_RestartSkip, __RC__)

            ! convective scavenging
            write (field_name, '(A, I0.3)') 'SSSV', n

            call MAPL_AddImportSpec(GC,                    &
               SHORT_NAME = 'clim'//field_name,            &
               LONG_NAME  = r%vtitle(n),                   &
               UNITS      = r%vunits(n),                   &
               DIMS       = MAPL_DimsHorzOnly,             &
               VLOCATION  = MAPL_VLocationCenter,          &  
               RESTART    = MAPL_RestartSkip, __RC__)
        end do
    end if IMPORT_SS_DEP_FLUXES

else

!   3-D Quantities
!   --------------    

!   GMICHEM species
!   ---------------
    GMI_on: if (state%chemReg%doing_GMI) then

     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME = 'OH',                                 &
        LONG_NAME  = 'Hydroxyl_radical',                   &
        UNITS      = 'mol/mol',                            &
        DIMS       = MAPL_DimsHorzVert,                    &
        VLOCATION  = MAPL_VLocationCenter,                 &
        RESTART    = MAPL_RestartSkip,     __RC__)

     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME = 'CH4',                                &
        LONG_NAME  = 'Methane',                            &
        UNITS      = 'mol/mol',                            &
        DIMS       = MAPL_DimsHorzVert,                    &
        VLOCATION  = MAPL_VLocationCenter,                 &
        RESTART    = MAPL_RestartSkip,     __RC__)

     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME = 'H2O2',                               &
        LONG_NAME  = 'Hydrogen_peroxide',                  &
        UNITS      = 'mol/mol',                            &
        DIMS       = MAPL_DimsHorzVert,                    &
        VLOCATION  = MAPL_VLocationCenter,                 &
        RESTART    = MAPL_RestartSkip,     __RC__)

     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME = 'NO3',                                &
        LONG_NAME  = 'Nitrogen_trioxide',                  &
        UNITS      = 'mol/mol',                            &
        DIMS       = MAPL_DimsHorzVert,                    &
        VLOCATION  = MAPL_VLocationCenter,                 &
        RESTART    = MAPL_RestartSkip,     __RC__)

    end if GMI_on

!    Pressure at layer edges
!    -----------------------
    call MAPL_AddImportSpec(GC,                            &
       SHORT_NAME = 'PLE',                                 &
       LONG_NAME  = 'air_pressure',                        &
       UNITS      = 'Pa',                                  &
       DIMS       = MAPL_DimsHorzVert,                     &
       VLOCATION  = MAPL_VLocationEdge,                    &
       RESTART    = MAPL_RestartSkip,     __RC__) 

!   Pressure thickness
!   ------------------
    call MAPL_AddImportSpec(GC,                            &
       SHORT_NAME = 'DELP',                                &
       LONG_NAME  = 'pressure_thickness',                  &
       UNITS      = 'Pa',                                  &
       DIMS       = MAPL_DimsHorzVert,                     &
       VLOCATION  = MAPL_VLocationCenter,                  &
       RESTART    = MAPL_RestartSkip,     __RC__)

!   Height at the edges
!   -------------------
    call MAPL_AddImportSpec(GC,                            &
       SHORT_NAME = 'ZLE',                                 &
       LONG_NAME  = 'geopotential_height',                 &
       UNITS      = 'm',                                   &
       DIMS       = MAPL_DimsHorzVert,                     &
       VLOCATION  = MAPL_VLocationEdge,                    &
       RESTART    = MAPL_RestartSkip,     __RC__)

!    AIRDENS: moist air density
!    --------------------------
     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME = 'AIRDENS',                            &
        LONG_NAME  = 'moist_air_density',                  &
        UNITS      = 'kg/m^3',                             &
        DIMS       = MAPL_DimsHorzVert,                    &
        VLOCATION  = MAPL_VLocationCenter,                 &
        RESTART    = MAPL_RestartSkip,    __RC__)

!    AIRDENS_DRY: dry air density
!    ----------------------------
     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME = 'AIRDENS_DRYP',                       &
        LONG_NAME  = 'partial_dry_air_density',            &
        UNITS      = 'kg dry m-3 tot',                     &
        DIMS       = MAPL_DimsHorzVert,                    &
        VLOCATION  = MAPL_VLocationCenter,                 &
        RESTART    = MAPL_RestartSkip,     __RC__)

!    CLOUD
!    -----
     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME = 'FCLD'  ,                             &
        LONG_NAME  = 'Cloud fraction for radiation',       &
        UNITS      = '1',                                  &
        DIMS       = MAPL_DimsHorzVert,                    &
        VLOCATION  = MAPL_VLocationCenter,                 &
        RESTART    = MAPL_RestartSkip,     __RC__)


!   T
!   -
    call MAPL_AddImportSpec(GC,                            &
       SHORT_NAME = 'T',                                   &
       LONG_NAME  = 'air_temperature',                     &
       UNITS      = 'K',                                   &
       DIMS       = MAPL_DimsHorzVert,                     &
       VLOCATION  = MAPL_VLocationCenter,                  &
       RESTART    = MAPL_RestartSkip,     __RC__)

!   U
!   -
    call MAPL_AddImportSpec(GC,                            &
       SHORT_NAME = 'U',                                   &
       LONG_NAME  = 'eastward_wind',                       &
       UNITS      = 'm s-1',                               &
       DIMS       = MAPL_DimsHorzVert,                     &
       VLOCATION  = MAPL_VLocationCenter,                  &
       RESTART    = MAPL_RestartSkip,     __RC__)

!   V
!   -
    call MAPL_AddImportSpec(GC,                            &
       SHORT_NAME = 'V',                                   &
       LONG_NAME  = 'northward_wind',                      &
       UNITS      = 'm s-1',                               &
       DIMS       = MAPL_DimsHorzVert,                     &
       VLOCATION  = MAPL_VLocationCenter,                  &
       RESTART    = MAPL_RestartSkip,     __RC__)

!   Q
!   -
    call MAPL_AddImportSpec(GC,                            &
       SHORT_NAME = 'Q',                                   &
       LONG_NAME  = 'specific_humidity',                   &
       UNITS      = 'kg kg-1',                             &
       DIMS       = MAPL_DimsHorzVert,                     &
       VLOCATION  = MAPL_VLocationCenter,                  &
       RESTART    = MAPL_RestartSkip,     __RC__)  

!   QITOT + QITOT
!   -------------
    call MAPL_AddImportSpec(GC,                            &
       SHORT_NAME = 'QCTOT',                               &
       LONG_NAME  = 'mass_fraction_of_total_cloud_water',  &
       UNITS      = 'kg kg-1',                             &
       DIMS       = MAPL_DimsHorzVert,                     &
       VLOCATION  = MAPL_VLocationCenter, __RC__)

!   RH: is between 0 and 1
!   ----------------------
    call MAPL_AddImportSpec(GC,                            &
       SHORT_NAME = 'RH2',                                 &
       LONG_NAME  = 'Rel_Hum_after_moist',                 &
       UNITS      = '1',                                   &
       DIMS       = MAPL_DimsHorzVert,                     &
       VLOCATION  = MAPL_VLocationCenter,                  &
       RESTART    = MAPL_RestartSkip,     __RC__)

!   PFI_LAAN: this is nonconvective precipition
!   ----------------------------------------
     call MAPL_AddImportSpec(GC,                           &
         SHORT_NAME='PFI_LSAN',                            &
         LONG_NAME ='3D_flux_of_ice_nonconvective_precipitation',  &
         UNITS     ='kg/m2/s',                             &
         DIMS      = MAPL_DimsHorzVert,                    &
         VLOCATION = MAPL_VLocationEdge,                   &
         RESTART   = MAPL_RestartSkip,    __RC__)

!    PFL_LAAN: this is nonconvective precipition
!    ----------------------------------------
     call MAPL_AddImportSpec(GC,                           &
         SHORT_NAME='PFL_LSAN',                            &
         LONG_NAME ='3D_flux_of_liquid_nonconvective_precipitation',  &
         UNITS     ='kg/m2/s',                             &
         DIMS      = MAPL_DimsHorzVert,                    &
         VLOCATION = MAPL_VLocationEdge,                   &
         RESTART   = MAPL_RestartSkip,    __RC__)

!   Ozone from PCHEM for CFC-12 photolysis
!   --------------------------------------
  CFC_on: if (state%chemReg%doing_CFC) then
    call MAPL_AddImportSpec(GC,                            &
       SHORT_NAME = 'O3',                                  &
       LONG_NAME  = 'ozone_mass_mixing_ratio',             &
       UNITS      = 'kg/kg',                               &
       DIMS       = MAPL_DimsHorzVert,                     &
       VLOCATION  = MAPL_VLocationCenter,                  &
       RESTART    = MAPL_RestartSkip,     __RC__)
  end if CFC_on

!    2-D Quantities
!    --------------    

!    TROPP - Connectivity from SDYN to PHYS is TROPP_BLENDED to TROPP
!    ----------------------------------------------------------------
     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME = 'TROPP',                              &
        LONG_NAME  = 'tropopause_pressure_based_on_blended_estimate', &
        UNITS      = 'Pa',                                 &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationNone,                   &
        RESTART    = MAPL_RestartSkip,    __RC__)

!    LWI
!    ---
     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME = 'LWI',                                &
        LONG_NAME  = 'land-ocean-ice_mask',                &
        UNITS      = '1',                                  &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationNone,                   &
      ! RESTART    = MAPL_RestartSkip,                     &
                                          __RC__)

!    PBL 
!    ---
     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME = 'ZPBL',                               &
        LONG_NAME  = 'Planetary boundary layer height',    &
        UNITS      = 'm',                                  &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationNone,                   &
        RESTART    = MAPL_RestartSkip,    __RC__)

!    FRACLAKE
!    --------
     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME = 'FRLAKE',                             &
        LONG_NAME  = 'fraction_of_lake',                   &
        UNITS      = '1',                                  &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationNone,                   &
        RESTART    = MAPL_RestartSkip,    __RC__)

!    FROCEAN
!    -------
     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME = 'FROCEAN',                            &
        LONG_NAME  = 'fraction_of_ocean',                  &
        UNITS      = '1',                                  &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationNone,                   &
        RESTART    = MAPL_RestartSkip,    __RC__)

!    FRACI
!    -----
     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME = 'FRACI',                              &
        LONG_NAME  = 'ice_covered_fraction_of_tile',       &
        UNITS      = '1',                                  &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationNone,                   &
      ! RESTART    = MAPL_RestartSkip,                     &
                                          __RC__)

!    GWETTOP
!    -------
     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME = 'WET1',                               &
        LONG_NAME  = 'surface_soil_wetness',               &
        UNITS      = '1',                                  &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationNone, __RC__)

!    LAI
!    ---
     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME = 'LAI',                                &
        LONG_NAME  = 'leaf_area_index',                    &
        UNITS      = '1',                                  &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationNone,                   &
        RESTART    = MAPL_RestartSkip,   __RC__)

!    This could be useful,  but it is not needed now
!    -----------------------------------------------
     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME = 'GRN',                                &
        LONG_NAME  = 'greeness_fraction',                  &
        UNITS      = '1',                                  &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationNone,                   &
        RESTART    = MAPL_RestartSkip,   __RC__)

!    PRECC: I hope this is defined over oceans
!    -----------------------------------------
     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME = 'CN_PRCP',                            &
        LONG_NAME  = 'Surface Conv. rain flux needed by land', &
        UNITS      = 'kg/m^2/s',                           &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationNone,                   &
        RESTART    = MAPL_RestartSkip,   __RC__)

!    PRECL: Non-convective precip, provided by Cinderella
!    ----------------------------------------------------
     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME = 'NCN_PRCP',                           &
        LONG_NAME  = 'Non-convective precipitation',       &
        UNITS      = 'kg/m^2/s',                           &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationNone,                   &
        RESTART    = MAPL_RestartSkip,   __RC__)

!    PS: from where???
!    -----------------
     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME = 'PS',                                 &
        LONG_NAME  = 'surface_pressure',                   &
        UNITS      = 'Pa',                                 &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationNone,                   &
        RESTART    = MAPL_RestartSkip,   __RC__)

!    SHFX (pos is up) - why not evap, Ri, ???
!    ----------------------------------------
     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME = 'SH',                                 &
        LONG_NAME  = 'sensible_heat_flux_from_turbulence', &
        UNITS      = 'W m-2',                              &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationNone,                   &
        RESTART    = MAPL_RestartSkip,   __RC__)

!    TA -- Surface Air Temperature
!    ----
     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME = 'TA',                                 &
        LONG_NAME  = 'surface_temperature_from_surface',   &
        UNITS      = 'K',                                  &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationNone,                   &
        RESTART    = MAPL_RestartSkip,   __RC__)

!    DZ -- Surface Mid-layer Height
!    ----
     call MAPL_AddImportSpec(GC,                           &
        LONG_NAME  = 'surface_layer_height',               &
        UNITS      = 'm',                                  &
        SHORT_NAME = 'DZ',                                 &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationNone,                   &
        RESTART    = MAPL_RestartSkip,   __RC__)

!    TSOIL1, from SURFACE
!    --------------------
     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME = 'TSOIL1',                             &
        LONG_NAME  = 'soil_temperatures_layer_1',          &
        UNITS      = 'K',                                  &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationNone,                   &
        RESTART    = MAPL_RestartSkip,   __RC__)

!    U10M
!    ----
     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME = 'U10M',                               &
        LONG_NAME  = '10-meter_eastward_wind',             &
        UNITS      = 'm s-1',                              &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationNone,                   &
        RESTART    = MAPL_RestartSkip,   __RC__)

!    V10M
!    ----
     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME = 'V10M',                               &
        LONG_NAME  = '10-meter_northward_wind',            &
        UNITS      = 'm s-1',                              &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationNone,                   &
        RESTART    = MAPL_RestartSkip,   __RC__)

!    U10N
!    ----
     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME = 'U10N',                               &
        LONG_NAME  = 'equivalent_neutral_10-meter_eastward_wind', &
        UNITS      = 'm s-1',                              &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationNone,                   &
        RESTART    = MAPL_RestartSkip,   __RC__)

!    V10N
!    ----
     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME = 'V10N',                               &
        LONG_NAME  = 'equivalent_neutral_10-meter_northward_wind', &
        UNITS      = 'm s-1',                              &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationNone,                   &
        RESTART    = MAPL_RestartSkip,   __RC__)

!    USTAR
!    -----
     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME = 'USTAR',                              &
        LONG_NAME  = 'surface_velocity_scale',             &
        UNITS      = 'm s-1',                              &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationNone,                   &
        RESTART    = MAPL_RestartSkip,   __RC__)

!    Z0H
!    ---
     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME = 'Z0H',                                &
        LONG_NAME  = 'surface_roughness_for_heat',         &
        UNITS      = 'm',                                  &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationNone,                   &
        RESTART    = MAPL_RestartSkip,   __RC__)

!    Cell area
!    ---------
     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME = 'AREA',                               &
        LONG_NAME  = 'agrid_cell_area',                    &
        UNITS      = 'm^2',                                &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationNone,                   &
        RESTART    = MAPL_RestartSkip,   __RC__)

     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME = 'TS',                                 &
        LONG_NAME  = 'surface skin temperature',           &
        UNITS      = 'K',                                  &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationNone,                   &
        RESTART    = MAPL_RestartSkip,   __RC__)

     call MAPL_AddImportSpec(GC,                           & 
         SHORT_NAME = 'CNV_MFD',                           & 
         LONG_NAME  = 'detraining_mass_flux',              &
         UNITS      = 'kg m-2 s-1',                        &    
         DIMS       = MAPL_DimsHorzVert,                   &  
         VLOCATION  = MAPL_VLocationCenter,                &
         RESTART    = MAPL_RestartSkip,  __RC__)

     call MAPL_AddImportSpec(GC,                           &                  
         SHORT_NAME = 'CNV_MFC',                           & 
         LONG_NAME  = 'cumulative_mass_flux',              &
         UNITS      = 'kg m-2 s-1',                        &
         DIMS       = MAPL_DimsHorzVert,                   &
         VLOCATION  = MAPL_VLocationEdge,                  &
         RESTART    = MAPL_RestartSkip,  __RC__)

     call MAPL_AddImportSpec(GC,                           &
         SHORT_NAME = 'QLCN',                              &
         LONG_NAME  = 'mass_fraction_of_convective_cloud_liquid_water', &
         UNITS      = 'kg kg-1',                           &
         DIMS       = MAPL_DimsHorzVert,                   &
         VLOCATION  = MAPL_VLocationCenter,                &
         RESTART    = MAPL_RestartSkip,  __RC__)

     call MAPL_AddImportSpec(GC,                           &
         SHORT_NAME = 'QICN',                              &
         LONG_NAME  = 'mass_fraction_of_convective_cloud_ice_water', &
         UNITS      = 'kg kg-1',                           &
         DIMS       = MAPL_DimsHorzVert,                   &
         VLOCATION  = MAPL_VLocationCenter,                &
         RESTART    = MAPL_RestartSkip,  __RC__)

     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME  = 'SWNDSRF',                           &
        LONG_NAME   = 'surface_net_downward_shortwave_flux', &
        UNITS       = 'W m^-2',                            &
        DIMS        = MAPL_DimsHorzOnly,                   &
        VLOCATION   = MAPL_VLocationNone,                  &
        RESTART     = MAPL_RestartSkip,  __RC__)
     
     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME = 'CLDTT',                              &
        LONG_NAME  = 'total_cloud_area_fraction',          &
        UNITS      = '1',                                  &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationNone,                   &
        RESTART    = MAPL_RestartSkip,   __RC__)


     call ESMF_ConfigGetAttribute(CF, DO_CO2CNNEE, label='USE_CNNEE:', default=0, __RC__)

     if (DO_CO2CNNEE == 1) then
        call MAPL_AddImportSpec(GC,                           &
             SHORT_NAME = 'CNNEE',                              &
             LONG_NAME  = 'CN_net_ecosystem_exchange',          &
             UNITS      = 'kg m-2 s-1',                         &
             DIMS       = MAPL_DimsHorzOnly,                    &
             VLOCATION  = MAPL_VLocationNone, __RC__)
     endif
     
end if GOCARTdata_IMPORTS


if ( r%doing_GOCART ) then


! !INTERNAL STATE:

!
!  NOTES: 
!    1)  vtitle as it stands is as the CF definition of long name.
!        I may need to create a "standard name" in chemReg and pass
!        this to GEOS Generic
!    2)  Host model MUST provide convective transport as well
!

    nq = r%nq     ! total number of chemical tracers
    
! Get BOOTSTRAP Default Values for GOCART INTERNAL
! ------------------------------------------------
    CALL ESMF_ConfigGetAttribute(CF, DEFVAL_CO2, Default=380.0e-6, Label='DEFVAL_CO2:', __RC__)

! Is GOCART providing O3 to the ANALYSIS bundle?
! ----------------------------------------------
    CALL ESMF_ConfigGetAttribute(CF, providerName, Default="PCHEM", Label='ANALYSIS_OX_PROVIDER:', __RC__)

! r%doing_O3 must be TRUE if the ANALYSIS_OX_PROVIDER is GOCART.
! --------------------------------------------------------------
    IF (providerName == 'GOCART' .AND. .NOT. r%doing_O3) THEN
     IF(MAPL_AM_I_ROOT()) THEN
      PRINT *," "
      PRINT *,TRIM(Iam)//": Set doing_O3 to yes in Chem_Registry.rc if GOCART "
      PRINT *,"             is the ANALYSIS_OX_PROVIDER."
      PRINT *," "
     END IF
     STATUS = 1
     VERIFY_(STATUS)
    END IF
     
!   Loop over all constituents on registry
!   --------------------------------------
    do n = r%i_GOCART, r%j_GOCART 

       if (state%data_driven) then 
          FRIENDLIES = trim(COMP_NAME)

          call ESMF_ConfigGetAttribute(CF, AEROFRIENDLY, Label='AERO_FRIENDLIES:', default=FRIENDLIES, __RC__)
       else
          if (trim(r%vname(n)) == 'OX' .and. trim(providerName) == 'GOCART') then
             FRIENDLIES = 'ANALYSIS:DYNAMICS:TURBULENCE:MOIST'
          else
             FRIENDLIES = 'DYNAMICS:TURBULENCE:MOIST'
          end if

!         Set aerosol friendly attribute to MOIST as function of Convective Parameterization
!         ----------------------------------------------------------------------------------

          short_name = ESMF_UtilStringUpperCase(trim(r%vname(n)))
          if ( short_name(1:2) .eq. 'DU'    .or. &
               short_name(1:2) .eq. 'SS'    .or. &
               short_name(1:2) .eq. 'OC'    .or. &
               short_name(1:3) .eq. 'BRC'   .or. &
               short_name(1:2) .eq. 'BC'    .or. &
               short_name(1:3) .eq. 'DMS'   .or. &
               short_name(1:3) .eq. 'SO2'   .or. &
               short_name(1:3) .eq. 'SO4'   .or. &
               short_name(1:3) .eq. 'MSA'   .or. &
               short_name(1:3) .eq. 'NH3'   .or. &
               short_name(1:4) .eq. 'NH4A'  .or. &
               short_name(1:5) .eq. 'NO3AN' ) then

             FRIENDLIES = 'DYNAMICS:TURBULENCE:MOIST'
             call ESMF_ConfigGetAttribute(CF, AEROFRIENDLY, Label='AERO_FRIENDLIES:', default=trim(FRIENDLIES), __RC__)

             if (index(trim(FRIENDLIES), 'MOIST') > 0)  call Disable_Convection
          endif

       end if ! data or computational GC

                                         DEFVAL = 0.0
       if ( short_name(1:3) .eq. 'CO2' ) DEFVAL = DEFVAL_CO2

!      Aerosol Tracers to be transported
!      ---------------------------------

       call MAPL_AddInternalSpec(GC,               &
          SHORT_NAME = trim(COMP_NAME)//'::'//trim(r%vname(n)), &
          LONG_NAME  = r%vtitle(n),                &
          UNITS      = r%vunits(n),                &     
          FRIENDLYTO = FRIENDLIES,                 &
          RESTART    = MAPL_RestartOptional,       &
          DEFAULT    = DEFVAL,                     &
          DIMS       = MAPL_DimsHorzVert,          &
          VLOCATION  = MAPL_VLocationCenter, __RC__)

    end do

!   This state is needed by radiation - It will contain 
!   aerosols and aerosol optics
!   --------------------------------------------------------
    call MAPL_AddExportSpec(GC,                    &
       SHORT_NAME = 'AERO',                        &
       LONG_NAME  = 'aerosol_mass_mixing_ratios',  &
       UNITS      = 'kg kg-1',                     &
       DIMS       = MAPL_DimsHorzVert,             &
       VLOCATION  = MAPL_VLocationCenter,          &
       DATATYPE   = MAPL_StateItem, __RC__)

!   This state is needed by MOIST - It will contain
!   aerosols
!   --------------------------------------------------------
    call MAPL_AddExportSpec(GC,                    &
       SHORT_NAME = 'AERO_ACI',                    &
       LONG_NAME  = 'aerosol_cloud_interaction',   &
       UNITS      = 'kg kg-1',                     &
       DIMS       = MAPL_DimsHorzVert,             &
       VLOCATION  = MAPL_VLocationCenter,          &
       DATATYPE   = MAPL_StateItem, __RC__)

!   This bundle is needed by surface for snow albedo modification
!   by aerosol settling and deposition
!   --------------------------------------------------------
    call MAPL_AddExportSpec(GC,                    &
       SHORT_NAME = 'AERO_DP',                     &
       LONG_NAME  = 'aerosol_deposition',          &
       UNITS      = 'kg m-2 s-1',                  &
       DIMS       = MAPL_DimsHorzOnly,             &
       DATATYPE   = MAPL_BundleItem, __RC__)

!   Export RH and DELP used in GOCART calculations
!   ----------------------------------------------
    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='RH2',                                         &
         LONG_NAME ='relative_humidity_after_moist',               &
         UNITS     ='1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'DELP',                                      &
         LONG_NAME  = 'pressure_thickness',                        &
         UNITS      = 'Pa',                                        &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

!   Diagnostic Exports over all aerosol tracers
!   -------------------------------------------

GOCART_COMPUTATIONAL_EXPORTS: if (.not. state%data_driven) then

    call MAPL_AddExportSpec(GC,               &
       SHORT_NAME = 'PSO4TOT',             &
       LONG_NAME  = 'Total Sulfate Produced in GOCART', &
       UNITS      = 'kg m-2 s-1',             &
       DIMS       = MAPL_DimsHorzVert,        &
       VLOCATION  = MAPL_VLocationCenter, __RC__)

    call MAPL_AddExportSpec(GC,               &
       SHORT_NAME = 'TOTEXTTAU',              &
       LONG_NAME  = 'Total Aerosol Extinction AOT [550 nm]', &
       UNITS      = '1',                      &
       DIMS       = MAPL_DimsHorzOnly,        &
       VLOCATION  = MAPL_VLocationNone, __RC__)

    call MAPL_AddExportSpec(GC,               &
       SHORT_NAME = 'TOTSCATAU',              &
       LONG_NAME  = 'Total Aerosol Scattering AOT [550 nm]', &
       UNITS      = '1',                      &
       DIMS       = MAPL_DimsHorzOnly,        &
       VLOCATION  = MAPL_VLocationNone, __RC__)

    call MAPL_AddExportSpec(GC,               &
       SHORT_NAME = 'TOTEXTT25',              &
       LONG_NAME  = 'Total Aerosol Extinction AOT [550 nm] - PM2.5', &
       UNITS      = '1',                      &
       DIMS       = MAPL_DimsHorzOnly,        &
       VLOCATION  = MAPL_VLocationNone, __RC__)

    call MAPL_AddExportSpec(GC,               &
       SHORT_NAME = 'TOTSCAT25',              &
       LONG_NAME  = 'Total Aerosol Scattering AOT [550 nm] - PM2.5', &
       UNITS      = '1',                      &
       DIMS       = MAPL_DimsHorzOnly,        &
       VLOCATION  = MAPL_VLocationNone, __RC__)

    call MAPL_AddExportSpec(GC,               &
       SHORT_NAME = 'TOTEXTTFM',              &
       LONG_NAME  = 'Total Aerosol Extinction AOT [550 nm] - PM1.0', &
       UNITS      = '1',                      &
       DIMS       = MAPL_DimsHorzOnly,        &
       VLOCATION  = MAPL_VLocationNone, __RC__)

    call MAPL_AddExportSpec(GC,               &
       SHORT_NAME = 'TOTSCATFM',              &
       LONG_NAME  = 'Total Aerosol Scattering AOT [550 nm] - PM1.0', &
       UNITS      = '1',                      &
       DIMS       = MAPL_DimsHorzOnly,        &
       VLOCATION  = MAPL_VLocationNone, __RC__)

    call MAPL_AddExportSpec(GC,               &
       SHORT_NAME = 'TOTANGSTR',              &
       LONG_NAME  = 'Total Aerosol Angstrom parameter [470-870 nm]', &
       UNITS      = '1',                      &
       DIMS       = MAPL_DimsHorzOnly,        &
       VLOCATION  = MAPL_VLocationNone, __RC__)

    call MAPL_AddExportSpec(GC,               &
       SHORT_NAME = 'PM',                     &
       LONG_NAME  = 'Total reconstructed PM', &
       UNITS      = 'kg m-3',                 &
       DIMS       = MAPL_DimsHorzOnly,        &
       VLOCATION  = MAPL_VLocationNone, __RC__)

    call MAPL_AddExportSpec(GC,               &
       SHORT_NAME = 'PM_RH35',                &
       LONG_NAME  = 'Total reconstructed PM(RH=35%)', &
       UNITS      = 'kg m-3',                 &
       DIMS       = MAPL_DimsHorzOnly,        &
       VLOCATION  = MAPL_VLocationNone, __RC__)

    call MAPL_AddExportSpec(GC,               &
       SHORT_NAME = 'PM_RH50',                &
       LONG_NAME  = 'Total reconstructed PM(RH=50%)', &
       UNITS      = 'kg m-3',                 &
       DIMS       = MAPL_DimsHorzOnly,        &
       VLOCATION  = MAPL_VLocationNone, __RC__)

    call MAPL_AddExportSpec(GC,               &
       SHORT_NAME = 'PM25',                   &
       LONG_NAME  = 'Total reconstructed PM2.5', &
       UNITS      = 'kg m-3',                 &
       DIMS       = MAPL_DimsHorzOnly,        &
       VLOCATION  = MAPL_VLocationNone, __RC__)

    call MAPL_AddExportSpec(GC,               &
       SHORT_NAME = 'PM25_RH35',              &
       LONG_NAME  = 'Total reconstructed PM2.5(RH=35%)', &
       UNITS      = 'kg m-3',                 &
       DIMS       = MAPL_DimsHorzOnly,        &
       VLOCATION  = MAPL_VLocationNone, __RC__)

    call MAPL_AddExportSpec(GC,               &
       SHORT_NAME = 'PM25_RH50',              &
       LONG_NAME  = 'Total reconstructed PM2.5(RH=50%)', &
       UNITS      = 'kg m-3',                 &
       DIMS       = MAPL_DimsHorzOnly,        &
       VLOCATION  = MAPL_VLocationNone, __RC__)

end if GOCART_COMPUTATIONAL_EXPORTS


!   Call Legacy Set Services
!   ----------------------
    if (.not.state%data_driven) then
       call Aero_GridCompSetServices ( gc, r, __RC__)
    end if


!!EOS

!   Set the Profiling timers
!   ------------------------
    call MAPL_TimerAdd ( GC, name = 'RUN',        __RC__ )
    call MAPL_TimerAdd ( GC, name = 'INITIALIZE', __RC__ )
    call MAPL_TimerAdd ( GC, name = 'FINALIZE',   __RC__ )
    call MAPL_TimerAdd ( GC, name = 'AERO1',   __RC__ )
    call MAPL_TimerAdd ( GC, name = 'AERO2',   __RC__ )

    call MAPL_TimerAdd (GC, name = 'SS', __RC__)
    call MAPL_TimerAdd (GC, name = 'O3', __RC__)
    call MAPL_TimerAdd (GC, name = 'DU', __RC__)
    call MAPL_TimerAdd (GC, name = 'BC', __RC__)
    call MAPL_TimerAdd (GC, name = 'OC', __RC__)
    call MAPL_TimerAdd (GC, name = 'SU', __RC__)
    call MAPL_TimerAdd (GC, name = 'CO', __RC__)
    call MAPL_TimerAdd (GC, name = 'CO2', __RC__)
    call MAPL_TimerAdd (GC, name = 'NI', __RC__)
    call MAPL_TimerAdd (GC, name = 'BRC', __RC__)
    call MAPL_TimerAdd (GC, name = 'CH4', __RC__)
    call MAPL_TimerAdd (GC, name = 'CFC', __RC__)

end if ! doing GOCART


!   Generic Set Services
!   --------------------
    call MAPL_GenericSetServices ( GC, __RC__ )

!   All done
!   --------

    RETURN_(ESMF_SUCCESS)
  
  end subroutine SetServices


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Initialize_ --- Initialize Aero_GridComp (ESMF)
!
! !INTERFACE:
!

   subroutine Initialize_ ( gc, impChem, expChem, clock, rc )

! !USES:

   implicit NONE

! !INPUT PARAMETERS:

   type(ESMF_Clock),  intent(inout) :: clock      ! The clock

! !OUTPUT PARAMETERS:

   type(ESMF_GridComp), intent(inout)  :: gc      ! Grid Component
   type(ESMF_State), intent(inout) :: impChem     ! Import State
   type(ESMF_State), intent(inout) :: expChem     ! Export State
   integer, intent(out) ::  rc                    ! Error return code:
                                                  !  0 - all is well
                                                  !  1 - 

! !DESCRIPTION: This is a simple ESMF wrapper.
!
! !REVISION HISTORY:
!
!  27Feb2005 da Silva  First crack.
!
!EOP
!-------------------------------------------------------------------------

!  ErrLog Variables
!  ----------------
   character(len=ESMF_MAXSTR)      :: IAm = 'Initialize_'
   integer                         :: STATUS
   character(len=ESMF_MAXSTR)      :: COMP_NAME

   type(Chem_Registry), pointer    :: chemReg
   type(Aero_GridComp), pointer    :: gcChem      ! Grid Component
   type(Chem_Bundle), pointer      :: w_c         ! Chemical tracer fields     
   integer                         :: nymd, nhms  ! time
   real                            :: cdt         ! chemistry timestep (secs)
   integer                         :: hdt         ! model     timestep (secs)
   integer                         :: rft

   type(ESMF_Grid)                 :: grid       
 
   integer                         :: i1=1, i2, ig=0, im  ! dist grid indices
   integer                         :: j1=1, j2, jg=0, jm  ! dist grid indices
   integer                         :: km, nq              ! dist grid indices
   integer                         :: n, dims(3), l

   type(ESMF_Config)               :: CF
   character(len=ESMF_MAXSTR)      :: diurnal_bb

   type(MAPL_MetaComp), pointer    :: ggState      ! GEOS Generic State
   type(GOCART_state),  pointer    :: myState       ! GOCART state
   type(ESMF_State)                :: internal
   type(ESMF_Field)                :: field
   type(ESMF_Field)                :: fld
   type(ESMF_FieldBundle)          :: bundle
   type(ESMF_State)                :: aero
   type(ESMF_FieldBundle)          :: aero_state_aerosols
   type(ESMF_State)                :: aero_aci
   type(ESMF_FieldBundle)          :: aero_aci_aerosols
   character(len=ESMF_MAXSTR)      :: fld_name
   integer                         :: n_aerosols
   integer                         :: n_modes
   integer, parameter              :: n_gocart_modes = 13
   !character(len=ESMF_MAXSTR) :: aero_aci_modes(n_gocart_modes)
   character(len=ESMF_MAXSTR), dimension(13),  parameter :: aero_aci_modes = (/'du001    ', 'du002    ', 'du003    ', &
                        'du004    ', 'du005    ',              &
                        'ss001    ', 'ss002    ', 'ss003    ', &
                        'sulforg01', 'sulforg02', 'sulforg03', &
                        'bcphilic ', 'ocphilic '/)
                        
   character(len=ESMF_MAXSTR), dimension(7), parameter :: aero_aci_modes_mam =  (/'NUM_A_ACC', 'NUM_A_AIT', 'NUM_A_CDU', 'NUM_A_CSS', &
                             'NUM_A_FDU', 'NUM_A_FSS', 'NUM_A_PCM'/)
                                             
   character(len=ESMF_MAXSTR)      :: short_name
   real                            :: f_aci_seasalt, maxclean, ccntuning
   character(LEN=ESMF_MAXSTR)      :: CLDMICRO
   
   logical                         :: USE_MAMNET
     
   type(MAPL_VarSpec), pointer     :: InternalSpec(:)
   integer                         :: instance

   real(ESMF_KIND_R4), pointer, dimension(:,:) :: LATS
   real(ESMF_KIND_R4), pointer, dimension(:,:) :: LONS
   real(ESMF_KIND_R4), pointer, dimension(:,:) :: CELL_AREA

   type(ESMF_Calendar)     :: calendar
   type(ESMF_Time)         :: currentTime
   type(ESMF_Alarm)        :: alarm_H2O2
   type(ESMF_Alarm)        :: alarm_HNO3
   type(ESMF_Time)         :: ringTime
   type(ESMF_TimeInterval) :: ringInterval
   integer                 :: year, month, day, hh, mm, ss

   type(ESMF_State)           :: providerState
   character(len=ESMF_MAXSTR) :: prefix
  
  real(ESMF_KIND_R4),  dimension(4) :: Vect_Hcts

!  Get my name and set-up traceback handle
!  ---------------------------------------
   call ESMF_GridCompGet( GC, NAME=COMP_NAME, CONFIG=CF, __RC__ )

   Iam = trim(COMP_NAME) // '::' // 'Initialize_'

   if (MAPL_AM_I_ROOT()) then
      print *, TRIM(Iam)//': Starting...'
      print *,' '
   end if

!  Get my internal MAPL_Generic state
!  -----------------------------------
   call MAPL_GetObjectFromGC ( GC, ggState, __RC__ )

   call MAPL_TimerOn(ggState, 'INITIALIZE')

!  Compute proper REFERENCE_TIME if not explicitly specified
!  ---------------------------------------------------------
   call MAPL_GetResource( ggState, hdt, Label='RUN_DT:',                                   __RC__ )
   call MAPL_GetResource( ggState, cdt, Label='GOCART_DT:',             default=real(hdt), __RC__ )
   call MAPL_GetResource( ggState, rft, Label='GOCART_REFERENCE_TIME:', default=-999,      __RC__ )

   if (rft == -999 .and. int(cdt) /= hdt) then
       hh  = int(  hdt/3600          )
       mm  = int( (hdt-hh*3600)/60   )
       ss  = int(  hdt-hh*3600-mm*60 )
       rft = hh*10000 + mm*100 + ss

       call MAPL_ConfigSetAttribute(cf, value=rft, Label="GOCART_REFERENCE_TIME:", __RC__ )
       if (MAPL_AM_I_ROOT()) write(*,"(21x,'Re-Setting GOCART_REFERENCE_TIME: ',i6.6)") rft
   endif


!  Initialize GEOS Generic
!  ------------------------
   call MAPL_GenericInitialize ( gc, impChem, expChem, clock, __RC__ )

   call MAPL_TimerOn(ggState, 'TOTAL')

!  Get pre-ESMF parameters from gc and clock
!  -----------------------------------------
   call extract_ ( gc, clock, chemReg, gcChem, w_c, nymd, nhms, cdt, STATUS, state=myState )
   VERIFY_(STATUS)


!  Create Chem Bundle
!  ------------------
   call ESMF_GridCompGet ( GC, grid=grid, __RC__ )

   call MAPL_GridGet ( grid, globalCellCountPerDim=DIMS, __RC__ )

   im = dims(1)
   jm = dims(2)
   nq = chemReg%nq

   call ESMF_GridGet(GRID, localDE = 0, &
                           staggerloc = ESMF_STAGGERLOC_CENTER, &
                           computationalCount = DIMS, __RC__)

!  Associate the Internal State fields with our legacy state 
!  ---------------------------------------------------------
   call MAPL_Get ( ggSTATE, INTERNALSPEC = InternalSpec, &
                            INTERNAL_ESMF_STATE = internal, &
                            LONS = LONS, &
                            LATS = LATS, __RC__ )

! A-Grid cell area
! ----------------
  if (myState%data_driven) then
      CELL_AREA => null()
  else
      call MAPL_GetPointer(impChem, NAME='AREA', ptr=CELL_AREA, __RC__)
  end if 

! Local sizes of three dimensions
!--------------------------------
   i2 = dims(1)
   j2 = dims(2)
   km = dims(3)

!  Initalize the legacy state but do not allocate memory for arrays
!  ----------------------------------------------------------------
   call Chem_BundleCreate_ ( chemReg, i1, i2, ig, im, j1, j2, jg, jm, km,  &
                             w_c, lon=LONS(:,:), lat=LATS(:,:), cell_area=CELL_AREA, &
                             skipAlloc=.true., __RC__ )

   w_c%grid_esmf = grid  ! Will need this for I/O later

!   Check whether to de-activate diurnal biomass burning (default is *on*)
!   ----------------------------------------------------------------------
    call ESMF_ConfigGetAttribute(CF, diurnal_bb, label='DIURNAL_BIOMASS_BURNING:', default='yes', __RC__)

    if ( diurnal_bb(1:3) .eq. 'yes' .or. &
         diurnal_bb(1:3) .eq. 'YES' .or. &
         diurnal_bb(1:3) .eq. 'Yes' ) then
         if (MAPL_AM_I_ROOT()) print *, trim(Iam)//': Diurnal Biomass Burning is ON'
         w_c%diurnal_bb = .true.
    else
         if (MAPL_AM_I_ROOT()) print *, trim(Iam)//': Diurnal Biomass Burning is OFF'
         w_c%diurnal_bb = .false.
    endif

!  Allocate these because they are not friendly tracers
!  ----------------------------------------------------
   allocate(w_c%delp(i1:i2,j1:j2,km), w_c%rh(i1:i2,j1:j2,km), __STAT__)

   _ASSERT( size(InternalSpec) == chemReg%n_GOCART, 'needs informative message' )

   do L = 1, size(InternalSpec)

      call MAPL_VarSpecGet(InternalSpec(L), SHORT_NAME=short_name, __RC__)

      N = chemReg%i_GOCART + L - 1
      call MAPL_GetPointer(internal, NAME=short_name, ptr=w_c%qa(N)%data3d, __RC__)

   end do

#ifdef PRINT_STATES

   if (MAPL_AM_I_ROOT()) then
       print *, trim(Iam)//': INTERNAL State during Initialize():' 
       call ESMF_StatePrint ( internal )
       print *, trim(Iam)//': IMPORT   State during Initialize():'
       call ESMF_StatePrint ( impChem  )
       print *, trim(Iam)//': EXPORT   State during Initialize():' 
       call ESMF_StatePrint ( expChem  )
    end if

#endif

!   Call Legacy Initialize
!   ----------------------
    call Aero_GridCompInitialize ( gcChem, w_c, gc, impChem, expChem, &
                                   nymd, nhms, cdt, myState%data_driven, STATUS )
    VERIFY_(STATUS)


!   Create H2O2 and HNO3 alarms
!   ---------------------------
    if (.not. myState%data_driven .and. w_c%reg%doing_SU) then

        call ESMF_ClockGet(clock, calendar=calendar, currTime=currentTime, RC=STATUS)
        VERIFY_(STATUS)

        call ESMF_TimeGet(currentTime, YY=year, MM=month, DD=day, H=hh, M=mm, S=ss, RC=STATUS)
        VERIFY_(STATUS)
        
        call ESMF_TimeSet(ringTime, YY=year, MM=month, DD=day, H=0, M=0, S=0, RC=STATUS)
        VERIFY_(STATUS)

        call ESMF_TimeIntervalSet(ringInterval, H=3, calendar=calendar, RC=STATUS)
        VERIFY_(STATUS)

        do while (ringTime < currentTime) 
            ringTime = currentTime + ringInterval
        end do

        alarm_H2O2 = ESMF_AlarmCreate(Clock        = clock,        &
                                      Name         = trim(H2O2_RECYCLE_ALARM), &
                                      RingInterval = ringInterval, &
                                      RingTime     = currentTime,  &
                                      Enabled      = .true.   ,    &
                                      Sticky       = .false.  ,    &
                                      RC           = STATUS)
        VERIFY_(STATUS)
    end if


    if (.not. myState%data_driven .and. w_c%reg%doing_NI) then
        call ESMF_ClockGet(clock, calendar=calendar, currTime=currentTime, RC=STATUS)
        VERIFY_(STATUS)

        call ESMF_TimeGet(currentTime, YY=year, MM=month, DD=day, H=hh, M=mm, S=ss, RC=STATUS)
        VERIFY_(STATUS)
        
        call ESMF_TimeSet(ringTime, YY=year, MM=month, DD=day, H=0, M=0, S=0, RC=STATUS)
        VERIFY_(STATUS)

        call ESMF_TimeIntervalSet(ringInterval, H=3, calendar=calendar, RC=STATUS)
        VERIFY_(STATUS)

        do while (ringTime < currentTime) 
            ringTime = currentTime + ringInterval
        end do

        alarm_HNO3 = ESMF_AlarmCreate(Clock        = clock,        &
                                      Name         = trim(HNO3_RECYCLE_ALARM), &
                                      RingInterval = ringInterval, &
                                      RingTime     = currentTime,  &
                                      Enabled      = .true.   ,    &
                                      Sticky       = .false.  ,    &
                                      RC           = STATUS)
        VERIFY_(STATUS)
    end if


!   Only at this point we have the scavenging coefficients filled,
!    so annotate the convection friendly internal state
!   Note: Move this to AddInternalSpec but first we need to have
!         the subcomponents as bonafide ESMF components
!-srf added Henrys law constants
!   --------------------------------------------------------------
    do n = ChemReg%i_GOCART, ChemReg%j_GOCART 
       call ESMF_StateGet(internal, trim(COMP_NAME)//'::'//trim(ChemReg%vname(n)), field, __RC__)
       call ESMF_AttributeSet(field, NAME='ScavengingFractionPerKm', VALUE=ChemReg%fscav(n), __RC__)
       Vect_Hcts(1:4)= ChemReg%hcts(1:4,n) 
  
       call ESMF_AttributeSet(field, 'SetofHenryLawCts', Vect_Hcts, __RC__)
    end do

!   Now that the internal state is nice and ready, add its contents and
!   attach aerosol optics method to the AERO state needed by radiation
!   ---------------------------------------------------------------------
    call ESMF_StateGet(expChem, 'AERO', aero, __RC__)

    ! This attribute indicates if the aerosol optics method is implemented or not. 
    ! Radiation will not call the aerosol optics method unless this attribute is 
    ! explicitly set to true.
    call ESMF_AttributeSet(aero, name='implements_aerosol_optics_method', value=.true., __RC__)

    aero_state_aerosols = ESMF_FieldBundleCreate(name='AEROSOLS', __RC__)
    call MAPL_StateAdd(aero, aero_state_aerosols, __RC__)

    do n = ChemReg%i_GOCART, ChemReg%j_GOCART

        short_name = ESMF_UtilStringUpperCase(trim(ChemReg%vname(n)))

        if ( short_name .eq. 'DU001'    .or. &
             short_name .eq. 'DU002'    .or. &
             short_name .eq. 'DU003'    .or. &
             short_name .eq. 'DU004'    .or. &
             short_name .eq. 'DU005'    .or. &
             short_name .eq. 'SS001'    .or. &
             short_name .eq. 'SS002'    .or. &
             short_name .eq. 'SS003'    .or. &
             short_name .eq. 'SS004'    .or. &
             short_name .eq. 'SS005'    .or. &
             short_name .eq. 'NO3AN1'   .or. &
             short_name .eq. 'NO3AN2'   .or. &
             short_name .eq. 'NO3AN3'   .or. &
             short_name .eq. 'OCPHOBIC' .or. &
             short_name .eq. 'OCPHILIC' .or. &
             short_name .eq. 'BRCPHOBIC' .or. &
             short_name .eq. 'BRCPHILIC' .or. &
             short_name .eq. 'BCPHOBIC' .or. &
             short_name .eq. 'BCPHILIC' .or. &
             short_name .eq. 'SO4'      .or. &
             short_name .eq. 'SO4V'     )    &
        then
           call ESMF_StateGet(INTERNAL,                     &
                              trim(COMP_NAME) // '::'//     &
                              trim(ChemReg%vname(n)),       &
                              FIELD, __RC__ )

           fld = MAPL_FieldCreate(FIELD, name=ChemReg%vname(n), __RC__)
           call MAPL_FieldBundleAdd(aero_state_aerosols, fld, __RC__)
        end if
    end do
    
    call ESMF_FieldBundleGet(aero_state_aerosols, fieldCount=n_aerosols, __RC__)

    if (n_aerosols > 0) then
        
        if (myState%data_driven) then
            instance = instanceData
        else
            instance = instanceComputational
        end if

        gocartMieTable(instance) = Chem_MieCreate(CF, __RC__)

        ! Mie Table instance/index
        call ESMF_AttributeSet(aero, name='mie_table_instance', value=instance, __RC__)

        ! state of the atmosphere
        call ESMF_AttributeSet(aero, name='air_pressure_for_aerosol_optics',             value='PLE', __RC__)
        call ESMF_AttributeSet(aero, name='relative_humidity_for_aerosol_optics',        value='RH',  __RC__)
        call ESMF_AttributeSet(aero, name='cloud_area_fraction_for_aerosol_optics',      value='',    __RC__) ! 'cloud_area_fraction_in_atmosphere_layer_for_aerosol_optics'

        ! aerosol optics
        call ESMF_AttributeSet(aero, name='band_for_aerosol_optics',                     value=0,     __RC__)
        call ESMF_AttributeSet(aero, name='extinction_in_air_due_to_ambient_aerosol',    value='EXT', __RC__)
        call ESMF_AttributeSet(aero, name='single_scattering_albedo_of_ambient_aerosol', value='SSA', __RC__)
        call ESMF_AttributeSet(aero, name='asymmetry_parameter_of_ambient_aerosol',      value='ASY', __RC__)

        ! add PLE to aero state
        call ESMF_AttributeGet(aero, name='air_pressure_for_aerosol_optics', value=fld_name, __RC__)
        if (fld_name /= '') then
            fld = MAPL_FieldCreateEmpty(trim(fld_name), w_c%grid_esmf, __RC__)

            call MAPL_FieldAllocCommit(fld, dims=MAPL_DimsHorzVert, location=MAPL_VLocationEdge, typekind=MAPL_R4, hw=0, __RC__)
            call MAPL_StateAdd(aero, fld, __RC__)
        end if

        ! add RH to Aero state
        call ESMF_AttributeGet(aero, name='relative_humidity_for_aerosol_optics', value=fld_name, __RC__)
        if (fld_name /= '') then
            fld = MAPL_FieldCreateEmpty(trim(fld_name), w_c%grid_esmf, __RC__)

            call MAPL_FieldAllocCommit(fld, dims=MAPL_DimsHorzVert, location=MAPL_VLocationCenter, typekind=MAPL_R4, hw=0, __RC__)
            call MAPL_StateAdd(aero, fld, __RC__)
        end if

        ! add EXT to aero state
        call ESMF_AttributeGet(aero, name='extinction_in_air_due_to_ambient_aerosol', value=fld_name, __RC__)
        if (fld_name /= '') then 
            fld = MAPL_FieldCreateEmpty(trim(fld_name), w_c%grid_esmf, __RC__)

            call MAPL_FieldAllocCommit(fld, dims=MAPL_DimsHorzVert, location=MAPL_VLocationCenter, typekind=MAPL_R4, hw=0, __RC__)            
            call MAPL_StateAdd(aero, fld, __RC__)
        end if

        ! add SSA to aero state
        call ESMF_AttributeGet(aero, name='single_scattering_albedo_of_ambient_aerosol', value=fld_name, __RC__)
        if (fld_name /= '') then
            fld = MAPL_FieldCreateEmpty(trim(fld_name), w_c%grid_esmf, __RC__)

            call MAPL_FieldAllocCommit(fld, dims=MAPL_DimsHorzVert, location=MAPL_VLocationCenter, typekind=MAPL_R4, hw=0, __RC__)
            call MAPL_StateAdd(aero, fld, __RC__)
        end if

        ! add ASY to aero state
        call ESMF_AttributeGet(aero, name='asymmetry_parameter_of_ambient_aerosol', value=fld_name, RC=STATUS)
        if (fld_name /= '') then 
            fld = MAPL_FieldCreateEmpty(trim(fld_name), w_c%grid_esmf, __RC__)

            call MAPL_FieldAllocCommit(fld, dims=MAPL_DimsHorzVert, location=MAPL_VLocationCenter, typekind=MAPL_R4, hw=0, __RC__)
            call MAPL_StateAdd(aero, fld, __RC__)
        end if
       
        ! attach the aerosol optics method
        call ESMF_MethodAdd(aero, label='aerosol_optics', userRoutine=aerosol_optics, __RC__)

    end if

#ifdef PRINT_STATES
    if (MAPL_AM_I_ROOT()) then
        print *, trim(Iam)//': AERO State during Initialize():'
        call ESMF_StatePrint(aero, nestedFlag=.true., __RC__)
    end if
#endif


!   Now that the internal state is nice and ready, add its contents and
!   attach aerosol-cloud interaction method to the AERO_ACI state needed by moist
!   ---------------------------------------------------------------------
    call ESMF_StateGet(expChem, 'AERO_ACI', aero_aci, __RC__)

    ! This attribute indicates if the aerosol optics method is implemented or not. 
    ! Radiation will not call the aerosol optics method unless this attribute is 
    ! explicitly set to true.
    call ESMF_AttributeSet(aero_aci, name='implements_aerosol_activation_properties_method', value=.true., __RC__)

    aero_aci_aerosols = ESMF_FieldBundleCreate(name='AEROSOLS', __RC__)
    call MAPL_StateAdd(aero_aci, aero_aci_aerosols, __RC__)

    do n = ChemReg%i_GOCART, ChemReg%j_GOCART 
        short_name = ESMF_UtilStringUpperCase(trim(ChemReg%vname(n)))

        if ( short_name .eq. 'DU001'    .or. &
             short_name .eq. 'DU002'    .or. &
             short_name .eq. 'DU003'    .or. &
             short_name .eq. 'DU004'    .or. &
             short_name .eq. 'DU005'    .or. &
             short_name .eq. 'SS001'    .or. &
             short_name .eq. 'SS002'    .or. &
             short_name .eq. 'SS003'    .or. &
             short_name .eq. 'SS004'    .or. &
             short_name .eq. 'SS005'    .or. &
!!           short_name .eq. 'NO3AN1'   .or. &
!!           short_name .eq. 'NO3AN2'   .or. &
!!           short_name .eq. 'NO3AN3'   .or. &
             short_name .eq. 'OCPHOBIC' .or. &
             short_name .eq. 'OCPHILIC' .or. &
             short_name .eq. 'BCPHOBIC' .or. &
             short_name .eq. 'BCPHILIC' .or. &
             short_name .eq. 'SO4'      .or. &
             short_name .eq. 'SO4V'     )    &
        then
            call ESMF_StateGet(INTERNAL,                     &
                               trim(COMP_NAME) // '::'//     &
                               trim(ChemReg%vname(n)),       &
                               FIELD, __RC__ )

            fld = MAPL_FieldCreate(FIELD, name=ChemReg%vname(n), __RC__)
            call MAPL_FieldBundleAdd(aero_aci_aerosols, fld, __RC__)
        end if
    end do

    !
    ! NOTE: The implementation of aerosol-cloud interaction in case of GOCART aerosols
    ! treats every aerosol tracer as a distinctive aerosol mode.
    !

    ! Following the aerosol-cloud-interaction state protocol, next steps are:  
    !  - attach a list with the aerosol modes
    !  - attach required met fields
    !  - attach method that computes the aerosol activation properties
    
    call ESMF_ConfigGetAttribute(CF, USE_MAMNET, default=.FALSE., label='USE_MAMNET:', __RC__)
    
    call ESMF_FieldBundleGet(aero_aci_aerosols, fieldCount=n_aerosols, __RC__)


    if (.not. USE_MAMNET) then 
        n_modes =  13    	
    else  	
        n_modes =  7
    end if                     

   !n_modes = size(aero_aci_modes)


    if (n_modes > 0 .and. n_aerosols > 0) then
      
        call ESMF_AttributeSet(aero_aci, name='number_of_aerosol_modes', value=n_modes, __RC__)
       
       
       if (.not. USE_MAMNET) then 
        call ESMF_AttributeSet(aero_aci, name='aerosol_modes', itemcount=n_modes, valuelist=aero_aci_modes, __RC__)
       else
         call ESMF_AttributeSet(aero_aci, name='aerosol_modes', itemcount=n_modes, valuelist=aero_aci_modes_mam, __RC__)
       end if  
   
        ! met fields and land fraction
        call ESMF_AttributeSet(aero_aci, name='air_pressure',                 value='PLE',      __RC__)
        call ESMF_AttributeSet(aero_aci, name='air_density',              value='AIRDEN',        __RC__)
        call ESMF_AttributeSet(aero_aci, name='air_temperature',              value='T',        __RC__)
        call ESMF_AttributeSet(aero_aci, name='fraction_of_land_type',        value='FRLAND',   __RC__)

        ! max mixing ratio before switching to "polluted" size distributions
        call ESMF_ConfigGetAttribute(CF, maxclean, default=1.0e-9, label='MAXCLEAN:', __RC__)
        call ESMF_AttributeSet(aero_aci, name='max_q_clean', value=maxclean, __RC__)
        
        call ESMF_ConfigGetAttribute(CF, CCNtuning, default=1.8, label='CCNTUNING:', __RC__)        
        call ESMF_AttributeSet(aero_aci, name='ccn_tuning', value=CCNtuning, __RC__)
       
        call ESMF_ConfigGetAttribute( CF, CLDMICRO, Label='CLDMICRO:',  default="1MOMENT", RC=STATUS)
        call ESMF_AttributeSet(aero_aci, name='cldmicro', value=CLDMICRO, __RC__)

        ! scaling factor for sea salt
        if(adjustl(CLDMICRO)=="2MOMENT") then
          call ESMF_ConfigGetAttribute(CF, f_aci_seasalt, default=4.0, label='SS_SCALE:', __RC__)
          call ESMF_AttributeSet(aero_aci, name='seasalt_scaling_factor', value=f_aci_seasalt, __RC__)
        else
          ! scaling factor for sea salt
          call ESMF_ConfigGetAttribute(CF, f_aci_seasalt, default=14.0, label='SS_SCALE:', __RC__)
          call ESMF_AttributeSet(aero_aci, name='seasalt_scaling_factor', value=f_aci_seasalt, __RC__)
        endif

        ! aerosol activation properties
        call ESMF_AttributeSet(aero_aci, name='width_of_aerosol_mode',        value='SIGMA',    __RC__)
        call ESMF_AttributeSet(aero_aci, name='aerosol_number_concentration', value='NUM',      __RC__)
        call ESMF_AttributeSet(aero_aci, name='aerosol_dry_size',             value='DGN',      __RC__)
        call ESMF_AttributeSet(aero_aci, name='aerosol_density',              value='density',  __RC__)
        call ESMF_AttributeSet(aero_aci, name='aerosol_hygroscopicity',       value='KAPPA',    __RC__)
        call ESMF_AttributeSet(aero_aci, name='fraction_of_dust_aerosol',     value='FDUST',    __RC__)
        call ESMF_AttributeSet(aero_aci, name='fraction_of_soot_aerosol',     value='FSOOT',    __RC__)
        call ESMF_AttributeSet(aero_aci, name='fraction_of_organic_aerosol',  value='FORGANIC', __RC__)


        ! add PLE to ACI state
        call ESMF_AttributeGet(aero_aci, name='air_pressure', value=fld_name, __RC__)
        if (fld_name /= '') then
            fld = MAPL_FieldCreateEmpty(trim(fld_name), w_c%grid_esmf, __RC__)

            call MAPL_FieldAllocCommit(fld, dims=MAPL_DimsHorzVert, location=MAPL_VLocationEdge, typekind=MAPL_R4, hw=0, __RC__)
            call MAPL_StateAdd(aero_aci, fld, __RC__)
        end if
        
        
        ! add AIRDENSITY to ACI state
        call ESMF_AttributeGet(aero_aci, name='air_density', value=fld_name, __RC__)
        if (fld_name /= '') then
            fld = MAPL_FieldCreateEmpty(trim(fld_name), w_c%grid_esmf, __RC__)

            call MAPL_FieldAllocCommit(fld, dims=MAPL_DimsHorzVert, location=MAPL_VLocationCenter, typekind=MAPL_R4, hw=0, __RC__)
            call MAPL_StateAdd(aero_aci, fld, __RC__)
        end if

        ! add T to ACI state
        call ESMF_AttributeGet(aero_aci, name='air_temperature', value=fld_name, __RC__)
        if (fld_name /= '') then
            fld = MAPL_FieldCreateEmpty(trim(fld_name), w_c%grid_esmf, __RC__)

            call MAPL_FieldAllocCommit(fld, dims=MAPL_DimsHorzVert, location=MAPL_VLocationCenter, typekind=MAPL_R4, hw=0, __RC__)
            call MAPL_StateAdd(aero_aci, fld, __RC__)
        end if

        ! add FRLAND to ACI state
        call ESMF_AttributeGet(aero_aci, name='fraction_of_land_type', value=fld_name, __RC__)
        if (fld_name /= '') then
            fld = MAPL_FieldCreateEmpty(trim(fld_name), w_c%grid_esmf, __RC__)

            call MAPL_FieldAllocCommit(fld, dims=MAPL_DimsHorzOnly, location=MAPL_VLocationCenter, typekind=MAPL_R4, hw=0, __RC__)
            call MAPL_StateAdd(aero_aci, fld, __RC__)
        end if


        ! add aerosol activation properties to ACI state
        call ESMF_AttributeGet(aero_aci, name='width_of_aerosol_mode', value=fld_name, __RC__)
        if (fld_name /= '') then
            fld = MAPL_FieldCreateEmpty(trim(fld_name), w_c%grid_esmf, __RC__)

            call MAPL_FieldAllocCommit(fld, dims=MAPL_DimsHorzVert, location=MAPL_VLocationCenter, typekind=MAPL_R4, hw=0, __RC__)
            call MAPL_StateAdd(aero_aci, fld, __RC__)
        end if

        call ESMF_AttributeGet(aero_aci, name='aerosol_number_concentration', value=fld_name, __RC__)
        if (fld_name /= '') then
            fld = MAPL_FieldCreateEmpty(trim(fld_name), w_c%grid_esmf, __RC__)

            call MAPL_FieldAllocCommit(fld, dims=MAPL_DimsHorzVert, location=MAPL_VLocationCenter, typekind=MAPL_R4, hw=0, __RC__)
            call MAPL_StateAdd(aero_aci, fld, __RC__)
        end if

        call ESMF_AttributeGet(aero_aci, name='aerosol_dry_size', value=fld_name, __RC__)
        if (fld_name /= '') then
            fld = MAPL_FieldCreateEmpty(trim(fld_name), w_c%grid_esmf, __RC__)

            call MAPL_FieldAllocCommit(fld, dims=MAPL_DimsHorzVert, location=MAPL_VLocationCenter, typekind=MAPL_R4, hw=0, __RC__)
            call MAPL_StateAdd(aero_aci, fld, __RC__)
        end if

        call ESMF_AttributeGet(aero_aci, name='aerosol_density', value=fld_name, __RC__)
        if (fld_name /= '') then
            fld = MAPL_FieldCreateEmpty(trim(fld_name), w_c%grid_esmf, __RC__)

            call MAPL_FieldAllocCommit(fld, dims=MAPL_DimsHorzVert, location=MAPL_VLocationCenter, typekind=MAPL_R4, hw=0, __RC__)
            call MAPL_StateAdd(aero_aci, fld, __RC__)
        end if

        call ESMF_AttributeGet(aero_aci, name='aerosol_hygroscopicity', value=fld_name, __RC__)
        if (fld_name /= '') then
            fld = MAPL_FieldCreateEmpty(trim(fld_name), w_c%grid_esmf, __RC__)

            call MAPL_FieldAllocCommit(fld, dims=MAPL_DimsHorzVert, location=MAPL_VLocationCenter, typekind=MAPL_R4, hw=0, __RC__)
            call MAPL_StateAdd(aero_aci, fld, __RC__)
        end if

        call ESMF_AttributeGet(aero_aci, name='fraction_of_dust_aerosol', value=fld_name, __RC__)
        if (fld_name /= '') then
            fld = MAPL_FieldCreateEmpty(trim(fld_name), w_c%grid_esmf, __RC__)

            call MAPL_FieldAllocCommit(fld, dims=MAPL_DimsHorzVert, location=MAPL_VLocationCenter, typekind=MAPL_R4, hw=0, __RC__)
            call MAPL_StateAdd(aero_aci, fld, __RC__)
        end if

        call ESMF_AttributeGet(aero_aci, name='fraction_of_soot_aerosol', value=fld_name, __RC__)
        if (fld_name /= '') then
            fld = MAPL_FieldCreateEmpty(trim(fld_name), w_c%grid_esmf, __RC__)

            call MAPL_FieldAllocCommit(fld, dims=MAPL_DimsHorzVert, location=MAPL_VLocationCenter, typekind=MAPL_R4, hw=0, __RC__)
            call MAPL_StateAdd(aero_aci, fld, __RC__)
        end if

        call ESMF_AttributeGet(aero_aci, name='fraction_of_organic_aerosol', value=fld_name, __RC__)
        if (fld_name /= '') then
            fld = MAPL_FieldCreateEmpty(trim(fld_name), w_c%grid_esmf, __RC__)

            call MAPL_FieldAllocCommit(fld, dims=MAPL_DimsHorzVert, location=MAPL_VLocationCenter, typekind=MAPL_R4, hw=0, __RC__)
            call MAPL_StateAdd(aero_aci, fld, __RC__)
        end if

        ! attach the aerosol optics method
        call ESMF_MethodAdd(aero_aci, label='aerosol_activation_properties', userRoutine=aerosol_activation_properties, __RC__)
        
        call ESMF_MethodAdd(aero_aci, label='MAMnet', userRoutine=MAMnet, __RC__)
    end if


!   Add settling and deposition to the AERO_DP bundle
!   -------------------------------------------------
    call ESMF_StateGet(expChem, 'AERO_DP', bundle, __RC__ )

!   If using GOCART.data, the data is provided in the import
!   state via ExtData versus the actual GOCART children
!   -------------------------------------------------------- 
    if ( myState%data_driven ) then
       providerState = impChem
       prefix = 'clim'
    else
       providerState = expChem
       prefix = ''
    end if

!   Dust
!   ----
    if ( ChemReg%doing_DU ) then

!      Dry deposition
!      --------------
       call AddFromExportToBundle_(providerState, bundle, trim(prefix)//'DUDP001', __RC__)
       call AddFromExportToBundle_(providerState, bundle, trim(prefix)//'DUDP002', __RC__)
       call AddFromExportToBundle_(providerState, bundle, trim(prefix)//'DUDP003', __RC__)
       call AddFromExportToBundle_(providerState, bundle, trim(prefix)//'DUDP004', __RC__)
       call AddFromExportToBundle_(providerState, bundle, trim(prefix)//'DUDP005', __RC__)

!      Wet deposition (Convective scavenging)
!      --------------------------------------
       call AddFromExportToBundle_(providerState, bundle, trim(prefix)//'DUSV001', __RC__)
       call AddFromExportToBundle_(providerState, bundle, trim(prefix)//'DUSV002', __RC__)
       call AddFromExportToBundle_(providerState, bundle, trim(prefix)//'DUSV003', __RC__)
       call AddFromExportToBundle_(providerState, bundle, trim(prefix)//'DUSV004', __RC__)
       call AddFromExportToBundle_(providerState, bundle, trim(prefix)//'DUSV005', __RC__)

!      Wet deposition
!      --------------
       call AddFromExportToBundle_(providerState, bundle, trim(prefix)//'DUWT001', __RC__)
       call AddFromExportToBundle_(providerState, bundle, trim(prefix)//'DUWT002', __RC__)
       call AddFromExportToBundle_(providerState, bundle, trim(prefix)//'DUWT003', __RC__)
       call AddFromExportToBundle_(providerState, bundle, trim(prefix)//'DUWT004', __RC__)
       call AddFromExportToBundle_(providerState, bundle, trim(prefix)//'DUWT005', __RC__)

!      Gravitational Settling
!      ----------------------
       call AddFromExportToBundle_(providerState, bundle, trim(prefix)//'DUSD001', __RC__)
       call AddFromExportToBundle_(providerState, bundle, trim(prefix)//'DUSD002', __RC__)
       call AddFromExportToBundle_(providerState, bundle, trim(prefix)//'DUSD003', __RC__)
       call AddFromExportToBundle_(providerState, bundle, trim(prefix)//'DUSD004', __RC__)
       call AddFromExportToBundle_(providerState, bundle, trim(prefix)//'DUSD005', __RC__)

    end if

!   Black Carbon
!   ------------
    if ( ChemReg%doing_BC ) then

!      Dry deposition
!      --------------
       call AddFromExportToBundle_(providerState, bundle, trim(prefix)//'BCDP001', __RC__)
       call AddFromExportToBundle_(providerState, bundle, trim(prefix)//'BCDP002', __RC__)

!      Wet deposition (Convective scavenging)
!      --------------------------------------
       call AddFromExportToBundle_(providerState, bundle, trim(prefix)//'BCSV001', __RC__)
       call AddFromExportToBundle_(providerState, bundle, trim(prefix)//'BCSV002', __RC__)

!      Wet deposition
!      --------------
       call AddFromExportToBundle_(providerState, bundle, trim(prefix)//'BCWT001', __RC__)
       call AddFromExportToBundle_(providerState, bundle, trim(prefix)//'BCWT002', __RC__)

!      Gravitational Settling
!      ----------------------
       call AddFromExportToBundle_(providerState, bundle, trim(prefix)//'BCSD001', __RC__)
       call AddFromExportToBundle_(providerState, bundle, trim(prefix)//'BCSD002', __RC__)

    end if

!   Organic Carbon
!   --------------
    if ( ChemReg%doing_OC ) then

!      Dry deposition
!      --------------
       call AddFromExportToBundle_(providerState, bundle, trim(prefix)//'OCDP001', __RC__)
       call AddFromExportToBundle_(providerState, bundle, trim(prefix)//'OCDP002', __RC__)

!      Wet deposition (Convective scavenging)
!      --------------------------------------
       call AddFromExportToBundle_(providerState, bundle, trim(prefix)//'OCSV001', __RC__)
       call AddFromExportToBundle_(providerState, bundle, trim(prefix)//'OCSV002', __RC__)

!      Wet deposition
!      --------------
       call AddFromExportToBundle_(providerState, bundle, trim(prefix)//'OCWT001', __RC__)
       call AddFromExportToBundle_(providerState, bundle, trim(prefix)//'OCWT002', __RC__)

!      Gravitational Settling
!      ----------------------
       call AddFromExportToBundle_(providerState, bundle, trim(prefix)//'OCSD001', __RC__)
       call AddFromExportToBundle_(providerState, bundle, trim(prefix)//'OCSD002', __RC__)

    end if

!   Sulfate (SO4; only aerosol component; bin 003)
!   ----------------------------------------------
    if ( ChemReg%doing_SU ) then

!      Dry deposition
!      --------------
       call AddFromExportToBundle_(providerState, bundle, trim(prefix)//'SUDP003', __RC__)

!      Wet deposition (Convective scavenging)
!      --------------
       call AddFromExportToBundle_(providerState, bundle, trim(prefix)//'SUSV003', __RC__)

!      Wet deposition
!      --------------
       call AddFromExportToBundle_(providerState, bundle, trim(prefix)//'SUWT003', __RC__)

!      Gravitational Settling
!      ----------------------
       call AddFromExportToBundle_(providerState, bundle, trim(prefix)//'SUSD003', __RC__)

    end if

!   Sea Salt
!   --------
    if ( ChemReg%doing_SS ) then

!      Dry deposition
!      --------------
       call AddFromExportToBundle_(providerState, bundle, trim(prefix)//'SSDP001', __RC__)
       call AddFromExportToBundle_(providerState, bundle, trim(prefix)//'SSDP002', __RC__)
       call AddFromExportToBundle_(providerState, bundle, trim(prefix)//'SSDP003', __RC__)
       call AddFromExportToBundle_(providerState, bundle, trim(prefix)//'SSDP004', __RC__)
       call AddFromExportToBundle_(providerState, bundle, trim(prefix)//'SSDP005', __RC__)

!      Wet deposition (Convective scavenging)
!      --------------------------------------
       call AddFromExportToBundle_(providerState, bundle, trim(prefix)//'SSSV001', __RC__)
       call AddFromExportToBundle_(providerState, bundle, trim(prefix)//'SSSV002', __RC__)
       call AddFromExportToBundle_(providerState, bundle, trim(prefix)//'SSSV003', __RC__)
       call AddFromExportToBundle_(providerState, bundle, trim(prefix)//'SSSV004', __RC__)
       call AddFromExportToBundle_(providerState, bundle, trim(prefix)//'SSSV005', __RC__)

!      Wet deposition
!      --------------
       call AddFromExportToBundle_(providerState, bundle, trim(prefix)//'SSWT001', __RC__)
       call AddFromExportToBundle_(providerState, bundle, trim(prefix)//'SSWT002', __RC__)
       call AddFromExportToBundle_(providerState, bundle, trim(prefix)//'SSWT003', __RC__)
       call AddFromExportToBundle_(providerState, bundle, trim(prefix)//'SSWT004', __RC__)
       call AddFromExportToBundle_(providerState, bundle, trim(prefix)//'SSWT005', __RC__)

!      Gravitational Settling
!      ----------------------
       call AddFromExportToBundle_(providerState, bundle, trim(prefix)//'SSSD001', __RC__)
       call AddFromExportToBundle_(providerState, bundle, trim(prefix)//'SSSD002', __RC__)
       call AddFromExportToBundle_(providerState, bundle, trim(prefix)//'SSSD003', __RC__)
       call AddFromExportToBundle_(providerState, bundle, trim(prefix)//'SSSD004', __RC__)
       call AddFromExportToBundle_(providerState, bundle, trim(prefix)//'SSSD005', __RC__)

    end if

#ifdef PRINT_STATES

   if (MAPL_AM_I_ROOT()) then
       print *, trim(Iam)//': AERO_DP Bundle during Initialize():' 
       call ESMF_FieldBundlePrint ( bundle )
   end if

#endif

    call MAPL_TimerOff(ggState, 'TOTAL')
    call MAPL_TimerOff(ggState, 'INITIALIZE')

    RETURN_(ESMF_SUCCESS)

CONTAINS

       subroutine AddFromExportToBundle_(STATE, BUNDLE, NAME, RC)
         type(ESMF_State)       :: STATE
         type(ESMF_FieldBundle) :: BUNDLE
         CHARACTER(LEN=*)       :: NAME
         integer, optional      :: RC
         type(ESMF_Field) :: FIELD
                   __Iam__('AddFromExportToBundle_')
         call ESMF_StateGet( STATE, NAME, FIELD, __RC__ )
         call MAPL_AllocateCoupling( FIELD, __RC__ )
         call MAPL_FieldBundleAdd ( BUNDLE, FIELD, __RC__ )
         RETURN_(ESMF_SUCCESS)
       end subroutine AddFromExportToBundle_

   end subroutine Initialize_


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Run1_ --- Runs Aero_GridComp (ESMF)
!
! !INTERFACE:
!

   subroutine Run1_ ( gc, impChem, expChem, clock, rc )

! !USES:

  implicit NONE

! !INPUT PARAMETERS:

   type(ESMF_Clock),  intent(inout) :: clock      ! The clock

! !OUTPUT PARAMETERS:

   type(ESMF_GridComp), intent(inout)  :: gc      ! Grid Component
   type(ESMF_State), intent(inout) :: impChem     ! Import State
   type(ESMF_State), intent(inout) :: expChem     ! Export State
   integer, intent(out) ::  rc                    ! Error return code:
                                                  !  0 - all is well
                                                  !  1 - 

! !DESCRIPTION: This is a simple ESMF wrapper.
!
! !REVISION HISTORY:
!
!  27Feb2005 da Silva  First crack.
!
!EOP
!-------------------------------------------------------------------------

!  ErrLog Variables
!  ----------------
   character(len=ESMF_MAXSTR)      :: IAm
   integer                         :: STATUS
   character(len=ESMF_MAXSTR)      :: COMP_NAME

   type(Chem_Registry), pointer    :: chemReg
   type(Aero_GridComp), pointer    :: gcChem      ! Grid Component
   type(Chem_Bundle), pointer      :: w_c         ! Chemical tracer fields     
   integer                         :: nymd, nhms  ! time
   real                            :: cdt         ! chemistry timestep (secs)
   real                            :: hdt         ! heartbeat time step (secs)
   real, pointer                   :: var(:,:,:)
   integer                         :: n

   type(ESMF_Config)               :: CF

   type(MAPL_MetaComp), pointer    :: ggState      ! GEOS Generic State
   type(ESMF_Alarm)                :: ALARM

   real(ESMF_KIND_R4), pointer, dimension(:,:) :: LATS
   real(ESMF_KIND_R4), pointer, dimension(:,:) :: LONS

   type (MAPL_SunOrbit)            :: ORBIT
   real, allocatable, target       :: ZTH(:,:)    ! can be R8
   real(ESMF_KIND_R4), allocatable :: r4ZTH(:,:)
   real(ESMF_KIND_R4), allocatable :: SLR(:,:)

   real, pointer                   :: rh2(:,:,:)
   integer                         :: in, jn

   type(GOCART_state), pointer     :: myState


!                               ---

!  Get my name and set-up traceback handle
!  ---------------------------------------
   call ESMF_GridCompGet( GC, NAME=COMP_NAME, CONFIG=CF, __RC__ )
   Iam = trim(COMP_NAME) // '::' // 'Run1_'

!  Get my internal MAPL_Generic state
!  -----------------------------------
   call MAPL_GetObjectFromGC ( GC, ggState, __RC__)

   call MAPL_TimerOn(ggState, 'TOTAL')
   call MAPL_TimerOn(ggState, 'RUN')

!  Get parameters from generic state.
!  ----------------------------------
   call MAPL_Get(ggState, LONS=LONS, LATS=LATS, ORBIT=ORBIT, RUNALARM=ALARM, __RC__)

!  Get heartbeat time step 
!  -----------------------
   call MAPL_GetResource(ggState, hdt, label='RUN_DT:', __RC__)
 
!  Get pre-ESMF parameters from gc and clock
!  -----------------------------------------
   call extract_ ( gc, clock, chemReg, gcChem, w_c, nymd, nhms, cdt, STATUS, state=myState )
   VERIFY_(STATUS)

!  Until all the gas phase species handle GOCART_DT correctly, we must run at the heartbeat
!  ----------------------------------------------------------------------------------------
!  Assume that DT is always an integral number of seconds
!  Add a fraction to both (and then truncate to int), to avoid cases like 900 /= 899.999999
   _ASSERT(abs(cdt-hdt) < 0.1, 'Implementation of GOCART_DT is problematic; set GOCART_DT = HEARTBEAT_DT')

   allocate(r4ZTH(SIZE(LATS,1), SIZE(LATS,2)), __STAT__)
   allocate(  ZTH(SIZE(LATS,1), SIZE(LATS,2)), __STAT__)
   allocate(  SLR(SIZE(LATS,1), SIZE(LATS,2)), __STAT__)

!  Update solar zenith angle
!  --------------------------
   call MAPL_SunGetInsolation(LONS, LATS, ORBIT, r4ZTH, SLR, CLOCK=CLOCK, __RC__)

!  Set pointers for sine/cosine zenith angle
!  -----------------------------------------

!  w_c%sinz => ...
   ZTH = r4ZTH
   w_c%cosz => zth

!  Fill in delp
!  ------------
   call MAPL_GetPointer ( impChem, var, 'DELP', __RC__ )
   w_c%delp = var(:,:,:)

!  Fill in RH
!  ----------
   call MAPL_GetPointer ( impChem, rh2, 'RH2', __RC__ )
   w_c%rh = rh2

!  Make sure tracers remain positive
!  ---------------------------------
   in = size(w_c%delp,1);   jn = size(w_c%delp,2)
   do n = ChemReg%i_GOCART, ChemReg%j_GOCART 
      call Chem_UtilNegFiller ( w_c%qa(n)%data3d, w_c%delp, in, jn, &
                                qmin=tiny(1.0) )
   end do

!  Call pre-ESMF version: runs at the heartbeat
!  --------------------------------------------
   call MAPL_TimerOn(ggState,'AERO1')
   call Aero_GridCompRun1 ( gcChem, w_c, gc, impChem, expChem, &
                            nymd, nhms, hdt, STATUS )
   VERIFY_(STATUS)
   call MAPL_TimerOff(ggState,'AERO1')

   deallocate(SLR,   __STAT__)
   deallocate(ZTH,   __STAT__)
   deallocate(r4ZTH, __STAT__)

   call MAPL_TimerOff(ggState, 'RUN')
   call MAPL_TimerOff(ggState, 'TOTAL')

   RETURN_(ESMF_SUCCESS)

   end subroutine Run1_


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Run_ --- Runs Aero_GridComp (ESMF)
!
! !INTERFACE:
!

   subroutine Run2_ ( gc, impChem, expChem, clock, rc )

! !USES:

  implicit NONE

! !INPUT PARAMETERS:

   type(ESMF_Clock),  intent(inout) :: clock      ! The clock

! !OUTPUT PARAMETERS:

   type(ESMF_GridComp), intent(inout)  :: gc      ! Grid Component
   type(ESMF_State), intent(inout) :: impChem     ! Import State
   type(ESMF_State), intent(inout) :: expChem     ! Export State
   integer, intent(out) ::  rc                    ! Error return code:
                                                  !  0 - all is well
                                                  !  1 - 

! !DESCRIPTION: This is a simple ESMF wrapper.
!
! !REVISION HISTORY:
!
!  27Feb2005 da Silva  First crack.
!
!EOP
!-------------------------------------------------------------------------

!  ErrLog Variables
!  ----------------
   character(len=ESMF_MAXSTR)      :: IAm
   integer                         :: STATUS
   character(len=ESMF_MAXSTR)      :: COMP_NAME

   type(Chem_Registry), pointer    :: chemReg
   type(Aero_GridComp), pointer    :: gcChem      ! Grid Component
   type(Chem_Bundle), pointer      :: w_c         ! Chemical tracer fields     
   integer                         :: nymd, nhms  ! time
   real                            :: hdt         ! heartbeat timestep (secs)
   real                            :: cdt         ! chemistry timestep (secs)
   real, pointer                   :: var(:,:,:)
   integer                         :: n

   type(ESMF_Config)               :: CF

   type(MAPL_MetaComp), pointer    :: ggState      ! GEOS Generic State
   type(ESMF_Alarm)                :: ALARM

   real(ESMF_KIND_R4), pointer, dimension(:,:) :: LATS
   real(ESMF_KIND_R4), pointer, dimension(:,:) :: LONS

   type (MAPL_SunOrbit)            :: ORBIT
   real, allocatable, target       :: ZTH(:,:)    ! can be R8
   real(ESMF_KIND_R4), allocatable :: r4ZTH(:,:)
   real(ESMF_KIND_R4), allocatable :: SLR(:,:)

   real, pointer                   :: rh2(:,:,:)
   integer                         :: in, jn

   type(ESMF_State)                :: internal
   type(GOCART_state), pointer     :: myState
   real, pointer, dimension(:,:,:) :: ptr3d_int
   real, pointer, dimension(:,:,:) :: ptr3d_imp
   
   logical                         :: run_alarm
   logical                         :: alarm_is_ringing


!  Diagnostics
   real, pointer, dimension(:,:)   :: totexttau, totscatau, &
                                      totextt25, totscat25, &
                                      totexttfm, totscatfm, &
                                      totangstr 
   real, pointer, dimension(:,:)   :: pm,        pm25,      &
                                      pm_rh35,   pm25_rh35, &
                                      pm_rh50,   pm25_rh50
   real, pointer, dimension(:,:)   :: duexttau, duscatau, &
                                      duextt25, duscat25, &
                                      duexttfm, duscatfm, &
                                      duangstr, dusmass,  &
                                      dusmass25
   real, pointer, dimension(:,:)   :: ssexttau, ssscatau, &
                                      ssextt25, ssscat25, &
                                      ssexttfm, ssscatfm, &
                                      ssangstr, sssmass,  &
                                      sssmass25
   real, pointer, dimension(:,:)   :: niexttau, niscatau, &
                                      niextt25, niscat25, &
                                      niexttfm, niscatfm, &
                                      niangstr, nismass,  &
                                      nismass25
   real, pointer, dimension(:,:)   :: nh4smass
   real, pointer, dimension(:,:)   :: suexttau, suscatau, &
                                      suangstr, susmass
   real, pointer, dimension(:,:)   :: suexttauvolc, suscatauvolc, &
                                      suangstrvolc, susmassvolc
   real, pointer, dimension(:,:)   :: bcexttau, bcscatau, &
                                      bcangstr, bcsmass
   real, pointer, dimension(:,:)   :: ocexttau, ocscatau, &
                                      ocangstr, ocsmass
   real, pointer, dimension(:,:)   :: brcexttau, brcscatau, &
                                      brcangstr, brcsmass
   real, pointer, dimension(:,:,:) :: rh2x, delpx
   real, pointer, dimension(:,:,:) :: pso4, pso4v, pso4t
   real, allocatable               :: tau1(:,:), tau2(:,:)
   real                            :: c1, c2, c3

!                               ---

!  Get my name and set-up traceback handle
!  ---------------------------------------
   call ESMF_GridCompGet( GC, NAME=COMP_NAME, CONFIG=CF, __RC__ )
   Iam = trim(COMP_NAME) // '::' // 'Run2_'

!  Get my internal MAPL_Generic state
!  -----------------------------------
   call MAPL_GetObjectFromGC ( GC, ggState, __RC__)

   call MAPL_TimerOn(ggState, 'TOTAL')
   call MAPL_TimerOn(ggState, 'RUN')

!  Get heartbeat time step 
!  -----------------------
   call MAPL_GetResource(ggState, hdt, label='RUN_DT:', __RC__)

!  Is time to recycle H2O2 and HNO3?
!  ---------------------------------
   call extract_ ( gc, clock, chemReg, gcChem, w_c, nymd, nhms, cdt, STATUS, state=myState )
   VERIFY_(STATUS)

   if (.not. myState%data_driven) then

       if (w_c%reg%doing_SU) then
           call ESMF_ClockGetAlarm(clock, trim(H2O2_RECYCLE_ALARM), alarm, __RC__)

           alarm_is_ringing = ESMF_AlarmIsRinging(alarm, __RC__)

           if (alarm_is_ringing) then
               do n = 1, gcChem%gcSU%n
                   if (.not. gcChem%gcSU%gcs(n)%using_GMI_H2O2) then
                       gcChem%gcSU%gcs(n)%recycle_H2O2 = .true.
                   end if
               end do

               call ESMF_AlarmRingerOff(alarm, __RC__)
           end if
       end if

       if (w_c%reg%doing_NI) then
           call ESMF_ClockGetAlarm(clock, trim(HNO3_RECYCLE_ALARM), alarm, __RC__)

           alarm_is_ringing = ESMF_AlarmIsRinging(alarm, __RC__)

           if (alarm_is_ringing) then
               do n = 1, gcChem%gcNI%n
                   gcChem%gcNI%gcs(n)%recycle_HNO3 = .true.
               end do

               call ESMF_AlarmRingerOff(alarm, __RC__)
           end if
       end if

   end if

!  Get parameters from generic state.
!  ----------------------------------
   call MAPL_Get(ggState, LONS=LONS, LATS=LATS, ORBIT=ORBIT, RUNALARM=ALARM, __RC__)

!  Get pre-ESMF parameters from gc and clock
!  -----------------------------------------
   call extract_ ( gc, clock, chemReg, gcChem, w_c, nymd, nhms, cdt, STATUS, state=myState )
   VERIFY_(STATUS)

   if (myState%data_driven) then
       
       call MAPL_Get ( ggState, INTERNAL_ESMF_STATE=internal, __RC__ )

       do n = chemReg%i_GOCART, chemReg%j_GOCART
           call MAPL_GetPointer ( internal, NAME=trim(COMP_NAME)//'::'//trim(chemReg%vname(n)), ptr=ptr3d_int, __RC__ )
           call MAPL_GetPointer ( impChem,  NAME='clim'//trim(chemReg%vname(n)), ptr=ptr3d_imp, __RC__ )
             
           ptr3d_int = ptr3d_imp
       end do

       call MAPL_TimerOff(ggState, 'RUN')
       call MAPL_TimerOff(ggState, 'TOTAL')

       RETURN_(ESMF_SUCCESS)
   end if 


   allocate(r4ZTH(SIZE(LATS,1), SIZE(LATS,2)), __STAT__)
   allocate(  ZTH(SIZE(LATS,1), SIZE(LATS,2)), __STAT__)
   allocate(  SLR(SIZE(LATS,1), SIZE(LATS,2)), __STAT__)

!  Update solar zenith angle
!  --------------------------
   call MAPL_SunGetInsolation(LONS, LATS, ORBIT, r4ZTH, SLR, CLOCK=CLOCK, __RC__)

!  Set pointers for sine/cosine zenith angle
!  -----------------------------------------

!  w_c%sinz => ...
   ZTH = r4ZTH
   w_c%cosz => zth

!  Fill in delp
!  ------------
   call MAPL_GetPointer ( impChem, var, 'DELP', __RC__ )
   w_c%delp = var(:,:,:)

!  Fill in RH
!  ----------
   call MAPL_GetPointer ( impChem, rh2, 'RH2', __RC__ )
   w_c%rh = rh2

!  Make sure tracers remain positive
!  ---------------------------------
   in = size(w_c%delp,1);   jn = size(w_c%delp,2)
   do n = ChemReg%i_GOCART, ChemReg%j_GOCART 
      call Chem_UtilNegFiller ( w_c%qa(n)%data3d, w_c%delp, in, jn, &
                                qmin=tiny(1.0) )
   end do

!  Call pre-ESMF version
!  ---------------------
   run_alarm = ESMF_AlarmIsRinging(ALARM, RC=STATUS)

   call MAPL_TimerOn(ggState,'AERO2')
   call Aero_GridCompRun2 ( gcChem, w_c, gc, impChem, expChem, &
                            run_alarm, nymd, nhms, cdt, STATUS )
   VERIFY_(STATUS)
   call MAPL_TimerOff(ggState,'AERO2')

   if (run_alarm) then
      call ESMF_AlarmRingerOff(ALARM, __RC__)
   end if

!  Get the diagnostics
   call MAPL_GetPointer (expChem, totexttau, 'TOTEXTTAU', __RC__)
   call MAPL_GetPointer (expChem, totscatau, 'TOTSCATAU', __RC__)
   call MAPL_GetPointer (expChem, totextt25, 'TOTEXTT25', __RC__)
   call MAPL_GetPointer (expChem, totscat25, 'TOTSCAT25', __RC__)
   call MAPL_GetPointer (expChem, totexttfm, 'TOTEXTTFM', __RC__)
   call MAPL_GetPointer (expChem, totscatfm, 'TOTSCATFM', __RC__)
   call MAPL_GetPointer (expChem, totangstr, 'TOTANGSTR', __RC__)

   ! dry PM
   call MAPL_GetPointer (expChem, pm25,      'PM25',      __RC__)
   call MAPL_GetPointer (expChem, pm,        'PM',        __RC__)
   call MAPL_GetPointer (expChem, rh2x,      'RH2',       __RC__)
   call MAPL_GetPointer (expChem, delpx,     'DELP',      __RC__)

   ! PM at RH=35%
   call MAPL_GetPointer (expChem, pm25_rh35, 'PM25_RH35', __RC__)
   call MAPL_GetPointer (expChem, pm_rh35,   'PM_RH35',   __RC__)

   ! PM at RH=50%
   call MAPL_GetPointer (expChem, pm25_rh50, 'PM25_RH50', __RC__)
   call MAPL_GetPointer (expChem, pm_rh50,   'PM_RH50',   __RC__)

   ! Sulfate produced (SO4-) in GOCART sulfur chemistry [kg m-2 s-1]
   call MAPL_GetPointer (expChem, pso4t,     'PSO4TOT',   __RC__)

   if(associated(totexttau)) totexttau(:,:) = 0.
   if(associated(totscatau)) totscatau(:,:) = 0.
   if(associated(totextt25)) totextt25(:,:) = 0.
   if(associated(totscat25)) totscat25(:,:) = 0.
   if(associated(totexttfm)) totexttfm(:,:) = 0.
   if(associated(totscatfm)) totscatfm(:,:) = 0.

   if(associated(pm))        pm(:,:)        = 0.
   if(associated(pm25))      pm25(:,:)      = 0.
   if(associated(pm_rh35))   pm_rh35(:,:)   = 0.
   if(associated(pm25_rh35)) pm25_rh35(:,:) = 0.
   if(associated(pm_rh50))   pm_rh50(:,:)   = 0.
   if(associated(pm25_rh50)) pm25_rh50(:,:) = 0.

   if(associated(rh2x))      rh2x           = w_c%rh
   if(associated(delpx))     delpx          = w_c%delp

   if(associated(pso4t))     pso4t(:,:,:)   = 0.

   if(w_c%reg%doing_du) then
     call MAPL_GetPointer (expChem, duexttau, 'DUEXTTAU', __RC__)
     call MAPL_GetPointer (expChem, duscatau, 'DUSCATAU', __RC__)
     call MAPL_GetPointer (expChem, duextt25, 'DUEXTT25', __RC__)
     call MAPL_GetPointer (expChem, duscat25, 'DUSCAT25', __RC__)
     call MAPL_GetPointer (expChem, duexttfm, 'DUEXTTFM', __RC__)
     call MAPL_GetPointer (expChem, duscatfm, 'DUSCATFM', __RC__)
     call MAPL_GetPointer (expChem, duangstr, 'DUANGSTR', __RC__)
     if(associated(totexttau) .and. associated(duexttau)) totexttau = totexttau+duexttau
     if(associated(totscatau) .and. associated(duscatau)) totscatau = totscatau+duscatau
     if(associated(totextt25) .and. associated(duextt25)) totextt25 = totextt25+duextt25
     if(associated(totscat25) .and. associated(duscat25)) totscat25 = totscat25+duscat25
     if(associated(totexttfm) .and. associated(duexttfm)) totexttfm = totexttfm+duexttfm
     if(associated(totscatfm) .and. associated(duscatfm)) totscatfm = totscatfm+duscatfm

     call MAPL_GetPointer (expChem, dusmass,   'DUSMASS',   __RC__)
     call MAPL_GetPointer (expChem, dusmass25, 'DUSMASS25', __RC__)
     if(associated(pm)        .and. associated(dusmass))   pm        = pm        + dusmass
     if(associated(pm25)      .and. associated(dusmass25)) pm25      = pm25      + dusmass25
     if(associated(pm_rh35)   .and. associated(dusmass))   pm_rh35   = pm_rh35   + dusmass
     if(associated(pm25_rh35) .and. associated(dusmass25)) pm25_rh35 = pm25_rh35 + dusmass25
     if(associated(pm_rh50)   .and. associated(dusmass))   pm_rh50   = pm_rh50   + dusmass
     if(associated(pm25_rh50) .and. associated(dusmass25)) pm25_rh50 = pm25_rh50 + dusmass25
   endif

   if(w_c%reg%doing_ss) then
     call MAPL_GetPointer (expChem, ssexttau, 'SSEXTTAU', __RC__)
     call MAPL_GetPointer (expChem, ssscatau, 'SSSCATAU', __RC__)
     call MAPL_GetPointer (expChem, ssextt25, 'SSEXTT25', __RC__)
     call MAPL_GetPointer (expChem, ssscat25, 'SSSCAT25', __RC__)
     call MAPL_GetPointer (expChem, ssexttfm, 'SSEXTTFM', __RC__)
     call MAPL_GetPointer (expChem, ssscatfm, 'SSSCATFM', __RC__)
     call MAPL_GetPointer (expChem, ssangstr, 'SSANGSTR', __RC__)
     if(associated(totexttau) .and. associated(ssexttau)) totexttau = totexttau+ssexttau
     if(associated(totscatau) .and. associated(ssscatau)) totscatau = totscatau+ssscatau
     if(associated(totextt25) .and. associated(ssextt25)) totextt25 = totextt25+ssextt25
     if(associated(totscat25) .and. associated(ssscat25)) totscat25 = totscat25+ssscat25
     if(associated(totexttfm) .and. associated(ssexttfm)) totexttfm = totexttfm+ssexttfm
     if(associated(totscatfm) .and. associated(ssscatfm)) totscatfm = totscatfm+ssscatfm

     call MAPL_GetPointer (expChem, sssmass,   'SSSMASS',   __RC__)
     call MAPL_GetPointer (expChem, sssmass25, 'SSSMASS25', __RC__)
     if(associated(pm)        .and. associated(sssmass))   pm        = pm        + sssmass
     if(associated(pm25)      .and. associated(sssmass25)) pm25      = pm25      + sssmass25
     if(associated(pm_rh35)   .and. associated(sssmass))   pm_rh35   = pm_rh35   + 1.86*sssmass
     if(associated(pm25_rh35) .and. associated(sssmass25)) pm25_rh35 = pm25_rh35 + 1.86*sssmass25
     if(associated(pm_rh50)   .and. associated(sssmass))   pm_rh50   = pm_rh50   + 2.42*sssmass
     if(associated(pm25_rh50) .and. associated(sssmass25)) pm25_rh50 = pm25_rh50 + 2.42*sssmass25
   endif

   if(w_c%reg%doing_ni) then
     call MAPL_GetPointer (expChem, niexttau, 'NIEXTTAU', __RC__)
     call MAPL_GetPointer (expChem, niscatau, 'NISCATAU', __RC__)
     call MAPL_GetPointer (expChem, niextt25, 'NIEXTT25', __RC__)
     call MAPL_GetPointer (expChem, niscat25, 'NISCAT25', __RC__)
     call MAPL_GetPointer (expChem, niexttfm, 'NIEXTTFM', __RC__)
     call MAPL_GetPointer (expChem, niscatfm, 'NISCATFM', __RC__)
     call MAPL_GetPointer (expChem, niangstr, 'NIANGSTR', __RC__)
     if(associated(totexttau) .and. associated(niexttau)) totexttau = totexttau+niexttau
     if(associated(totscatau) .and. associated(niscatau)) totscatau = totscatau+niscatau
     if(associated(totextt25) .and. associated(niextt25)) totextt25 = totextt25+niextt25
     if(associated(totscat25) .and. associated(niscat25)) totscat25 = totscat25+niscat25
     if(associated(totexttfm) .and. associated(niexttfm)) totexttfm = totexttfm+niexttfm
     if(associated(totscatfm) .and. associated(niscatfm)) totscatfm = totscatfm+niscatfm

     call MAPL_GetPointer (expChem, nismass,   'NISMASS',   __RC__)
     call MAPL_GetPointer (expChem, nismass25, 'NISMASS25', __RC__)
     call MAPL_GetPointer (expChem, nh4smass,  'NH4SMASS',  __RC__)
     if(associated(pm)        .and. associated(nismass)   .and. associated(nh4smass)) pm        = pm   + nismass   + nh4smass
     if(associated(pm25)      .and. associated(nismass25) .and. associated(nh4smass)) pm25      = pm25 + nismass25 + nh4smass
     if(associated(pm_rh35)   .and. associated(nismass)   .and. associated(nh4smass)) pm_rh35   = pm_rh35   + 1.33*(nismass   + nh4smass)
     if(associated(pm25_rh35) .and. associated(nismass25) .and. associated(nh4smass)) pm25_rh35 = pm25_rh35 + 1.33*(nismass25 + nh4smass)
     if(associated(pm_rh50)   .and. associated(nismass)   .and. associated(nh4smass)) pm_rh50   = pm_rh50   + 1.51*(nismass   + nh4smass)
     if(associated(pm25_rh50) .and. associated(nismass25) .and. associated(nh4smass)) pm25_rh50 = pm25_rh50 + 1.51*(nismass25 + nh4smass)
   endif

   if(w_c%reg%doing_su) then
     call MAPL_GetPointer (expChem, suexttau, 'SUEXTTAU', __RC__)
     call MAPL_GetPointer (expChem, suscatau, 'SUSCATAU', __RC__)
     call MAPL_GetPointer (expChem, suangstr, 'SUANGSTR', __RC__)
     call MAPL_GetPointer (expChem, suexttauvolc, 'SUEXTTAUvolc', __RC__)
     call MAPL_GetPointer (expChem, suscatauvolc, 'SUSCATAUvolc', __RC__)
     call MAPL_GetPointer (expChem, suangstrvolc, 'SUANGSTRvolc', __RC__)
     if(associated(totexttau) .and. associated(suexttau)) totexttau = totexttau+suexttau
     if(associated(totscatau) .and. associated(suscatau)) totscatau = totscatau+suscatau
     if(associated(totextt25) .and. associated(suexttau)) totextt25 = totextt25+suexttau
     if(associated(totscat25) .and. associated(suscatau)) totscat25 = totscat25+suscatau
     if(associated(totexttfm) .and. associated(suexttau)) totexttfm = totexttfm+suexttau
     if(associated(totscatfm) .and. associated(suscatau)) totscatfm = totscatfm+suscatau
!    Volcanic tracer is additive if present (should it go into PM2.5 as below?)
     if(associated(totexttau) .and. associated(suexttauvolc)) totexttau = totexttau+suexttauvolc
     if(associated(totscatau) .and. associated(suscatauvolc)) totscatau = totscatau+suscatauvolc
     if(associated(totextt25) .and. associated(suexttauvolc)) totextt25 = totextt25+suexttauvolc
     if(associated(totscat25) .and. associated(suscatauvolc)) totscat25 = totscat25+suscatauvolc
     if(associated(totexttfm) .and. associated(suexttauvolc)) totexttfm = totexttfm+suexttauvolc
     if(associated(totscatfm) .and. associated(suscatauvolc)) totscatfm = totscatfm+suscatauvolc

!    Sulfate production by GOCART chemistry (SO4-, kg m-2 s-1)
!    Philosophy here is there is always "full" sulfate instance
!    This "full" instance may be accompanied by "volc" instance
!    in which volcanoes are broken out from rest of "full" and 
!    the sulfur is additive.  Other tagged instances do not get
!    added.
     call MAPL_GetPointer (expChem, pso4,   'PSO4',       __RC__)
     call MAPL_GetPointer (expChem, pso4v,  'PSO4volc',   __RC__)
     if(associated(pso4t)) then
      if(associated(pso4))   pso4t = pso4t   + pso4
      if(associated(pso4v))  pso4t = pso4t   + pso4v
     endif

     call MAPL_GetPointer (expChem, susmass, 'SO4SMASS', __RC__)
     if (w_c%reg%doing_ni) then
         if(associated(pm)        .and. associated(susmass)) pm        = pm        + susmass
         if(associated(pm25)      .and. associated(susmass)) pm25      = pm25      + susmass
         if(associated(pm_rh35)   .and. associated(susmass)) pm_rh35   = pm_rh35   + 1.33*susmass
         if(associated(pm25_rh35) .and. associated(susmass)) pm25_rh35 = pm25_rh35 + 1.33*susmass
         if(associated(pm_rh50)   .and. associated(susmass)) pm_rh50   = pm_rh50   + 1.51*susmass
         if(associated(pm25_rh50) .and. associated(susmass)) pm25_rh50 = pm25_rh50 + 1.51*susmass
     else
         if(associated(pm)        .and. associated(susmass)) pm        = pm        + (132.14/96.06)*susmass
         if(associated(pm25)      .and. associated(susmass)) pm25      = pm25      + (132.14/96.06)*susmass
         if(associated(pm_rh35)   .and. associated(susmass)) pm_rh35   = pm_rh35   + 1.33*(132.14/96.06)*susmass
         if(associated(pm25_rh35) .and. associated(susmass)) pm25_rh35 = pm25_rh35 + 1.33*(132.14/96.06)*susmass
         if(associated(pm_rh50)   .and. associated(susmass)) pm_rh50   = pm_rh50   + 1.51*(132.14/96.06)*susmass
         if(associated(pm25_rh50) .and. associated(susmass)) pm25_rh50 = pm25_rh50 + 1.51*(132.14/96.06)*susmass
     endif
   endif

   if(w_c%reg%doing_bc) then
     call MAPL_GetPointer (expChem, bcexttau, 'BCEXTTAU', __RC__)
     call MAPL_GetPointer (expChem, bcscatau, 'BCSCATAU', __RC__)
     call MAPL_GetPointer (expChem, bcangstr, 'BCANGSTR', __RC__)
     if(associated(totexttau) .and. associated(bcexttau)) totexttau = totexttau+bcexttau
     if(associated(totscatau) .and. associated(bcscatau)) totscatau = totscatau+bcscatau
     if(associated(totextt25) .and. associated(bcexttau)) totextt25 = totextt25+bcexttau
     if(associated(totscat25) .and. associated(bcscatau)) totscat25 = totscat25+bcscatau
     if(associated(totexttfm) .and. associated(bcexttau)) totexttfm = totexttfm+bcexttau
     if(associated(totscatfm) .and. associated(bcscatau)) totscatfm = totscatfm+bcscatau

     call MAPL_GetPointer (expChem, bcsmass, 'BCSMASS', __RC__)
     if(associated(pm)        .and. associated(bcsmass)) pm        = pm        + bcsmass
     if(associated(pm25)      .and. associated(bcsmass)) pm25      = pm25      + bcsmass
     if(associated(pm_rh35)   .and. associated(bcsmass)) pm_rh35   = pm_rh35   + bcsmass
     if(associated(pm25_rh35) .and. associated(bcsmass)) pm25_rh35 = pm25_rh35 + bcsmass
     if(associated(pm_rh50)   .and. associated(bcsmass)) pm_rh50   = pm_rh50   + bcsmass
     if(associated(pm25_rh50) .and. associated(bcsmass)) pm25_rh50 = pm25_rh50 + bcsmass
   endif

   if(w_c%reg%doing_oc) then
     call MAPL_GetPointer (expChem, ocexttau, 'OCEXTTAU', __RC__)
     call MAPL_GetPointer (expChem, ocscatau, 'OCSCATAU', __RC__)
     call MAPL_GetPointer (expChem, ocangstr, 'OCANGSTR', __RC__)
     if(associated(totexttau) .and. associated(ocexttau)) totexttau = totexttau+ocexttau
     if(associated(totscatau) .and. associated(ocscatau)) totscatau = totscatau+ocscatau
     if(associated(totextt25) .and. associated(ocexttau)) totextt25 = totextt25+ocexttau
     if(associated(totscat25) .and. associated(ocscatau)) totscat25 = totscat25+ocscatau
     if(associated(totexttfm) .and. associated(ocexttau)) totexttfm = totexttfm+ocexttau
     if(associated(totscatfm) .and. associated(ocscatau)) totscatfm = totscatfm+ocscatau

     call MAPL_GetPointer (expChem, ocsmass, 'OCSMASS', __RC__)
     if(associated(pm)        .and. associated(ocsmass)) pm        = pm        + ocsmass
     if(associated(pm25)      .and. associated(ocsmass)) pm25      = pm25      + ocsmass
     if(associated(pm_rh35)   .and. associated(ocsmass)) pm_rh35   = pm_rh35   + 1.16*ocsmass  ! needs to be revisited: OCpho + 1.16 OCphi
     if(associated(pm25_rh35) .and. associated(ocsmass)) pm25_rh35 = pm25_rh35 + 1.16*ocsmass  ! 
     if(associated(pm_rh50)   .and. associated(ocsmass)) pm_rh50   = pm_rh50   + 1.24*ocsmass  ! needs to be revisited: OCpho + 1.24 OCphi
     if(associated(pm25_rh50) .and. associated(ocsmass)) pm25_rh50 = pm25_rh50 + 1.24*ocsmass  !
   endif

   if(w_c%reg%doing_brc) then
     call MAPL_GetPointer (expChem, brcexttau, 'BRCEXTTAU', __RC__)
     call MAPL_GetPointer (expChem, brcscatau, 'BRCSCATAU', __RC__)
     call MAPL_GetPointer (expChem, brcangstr, 'BRCANGSTR', __RC__)
     if(associated(totexttau) .and. associated(brcexttau)) totexttau = totexttau+brcexttau
     if(associated(totscatau) .and. associated(brcscatau)) totscatau = totscatau+brcscatau
     if(associated(totextt25) .and. associated(brcexttau)) totextt25 = totextt25+brcexttau
     if(associated(totscat25) .and. associated(brcscatau)) totscat25 = totscat25+brcscatau
     if(associated(totexttfm) .and. associated(brcexttau)) totexttfm = totexttfm+brcexttau
     if(associated(totscatfm) .and. associated(brcscatau)) totscatfm = totscatfm+brcscatau

     call MAPL_GetPointer (expChem, brcsmass, 'BRCSMASS', __RC__)
     if(associated(pm)        .and. associated(brcsmass)) pm        = pm        + brcsmass
     if(associated(pm25)      .and. associated(brcsmass)) pm25      = pm25      + brcsmass
     if(associated(pm_rh35)   .and. associated(brcsmass)) pm_rh35   = pm_rh35   + 1.16*brcsmass  ! needs to be revisited: BRCpho + 1.16 BRCphi
     if(associated(pm25_rh35) .and. associated(brcsmass)) pm25_rh35 = pm25_rh35 + 1.16*brcsmass  ! 
     if(associated(pm_rh50)   .and. associated(brcsmass)) pm_rh50   = pm_rh50   + 1.24*brcsmass  ! needs to be revisited: BRCpho + 1.24 BRCphi
     if(associated(pm25_rh50) .and. associated(brcsmass)) pm25_rh50 = pm25_rh50 + 1.24*brcsmass  !
   endif

   if(associated(totangstr)) then
    totangstr(:,:) = 0.0

    allocate(tau1(SIZE(LATS,1), SIZE(LATS,2)), &
             tau2(SIZE(LATS,1), SIZE(LATS,2)), __STAT__)

    tau1(:,:) = tiny(1.0)
    tau2(:,:) = tiny(1.0)
    c1 = -log(470./550.)
    c2 = -log(870./550.)
    c3 = -log(470./870.)
    if(w_c%reg%doing_du .and. associated(duexttau) .and. associated(duangstr)) then
     tau1 = tau1 + duexttau*exp(c1*duangstr)
     tau2 = tau2 + duexttau*exp(c2*duangstr)
    endif
    if(w_c%reg%doing_ss .and. associated(ssexttau) .and. associated(ssangstr)) then
     tau1 = tau1 + ssexttau*exp(c1*ssangstr)
     tau2 = tau2 + ssexttau*exp(c2*ssangstr)
    endif
    if(w_c%reg%doing_ni .and. associated(niexttau) .and. associated(niangstr)) then
     tau1 = tau1 + niexttau*exp(c1*niangstr)
     tau2 = tau2 + niexttau*exp(c2*niangstr)
    endif
    if(w_c%reg%doing_su .and. associated(suexttau) .and. associated(suangstr)) then
     tau1 = tau1 + suexttau*exp(c1*suangstr)
     tau2 = tau2 + suexttau*exp(c2*suangstr)
    endif
    if(w_c%reg%doing_su .and. associated(suexttauvolc) .and. associated(suangstrvolc)) then
     tau1 = tau1 + suexttauvolc*exp(c1*suangstrvolc)
     tau2 = tau2 + suexttauvolc*exp(c2*suangstrvolc)
    endif
    if(w_c%reg%doing_bc .and. associated(bcexttau) .and. associated(bcangstr)) then
     tau1 = tau1 + bcexttau*exp(c1*bcangstr)
     tau2 = tau2 + bcexttau*exp(c2*bcangstr)
    endif
    if(w_c%reg%doing_oc .and. associated(ocexttau) .and. associated(ocangstr)) then
     tau1 = tau1 + ocexttau*exp(c1*ocangstr)
     tau2 = tau2 + ocexttau*exp(c2*ocangstr)
    endif
    if(w_c%reg%doing_brc .and. associated(brcexttau) .and. associated(brcangstr)) then
     tau1 = tau1 + brcexttau*exp(c1*brcangstr)
     tau2 = tau2 + brcexttau*exp(c2*brcangstr)
    endif
    totangstr = log(tau1/tau2)/c3

    deallocate(tau1, tau2, __STAT__)
   endif

   deallocate(SLR,   __STAT__)
   deallocate(ZTH,   __STAT__)
   deallocate(r4ZTH, __STAT__)

   call MAPL_TimerOff(ggState, 'RUN')
   call MAPL_TimerOff(ggState, 'TOTAL')

   RETURN_(ESMF_SUCCESS)

   end subroutine Run2_


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Finalize_ --- Finalize Aero_GridComp (ESMF)
!
! !INTERFACE:
!

   subroutine Finalize_ ( gc, impChem, expChem, clock, rc )

! !USES:

  implicit NONE

! !INPUT PARAMETERS:

   type(ESMF_Clock),  intent(inout) :: clock      ! The clock

! !OUTPUT PARAMETERS:

   type(ESMF_GridComp), intent(inout)  :: gc      ! Grid Component
   type(ESMF_State), intent(inout) :: impChem     ! Import State
   type(ESMF_State), intent(inout) :: expChem     ! Export State
   integer, intent(out) ::  rc                    ! Error return code:
                                                  !  0 - all is well
                                                  !  1 - 

! !DESCRIPTION: This is a simple ESMF wrapper.
!
! !REVISION HISTORY:
!
!  27Feb2005 da Silva  First crack.
!
!EOP
!-------------------------------------------------------------------------


!  ErrLog Variables
!  ----------------
   character(len=ESMF_MAXSTR)      :: IAm = 'Finalize_'
   integer                         :: STATUS
   character(len=ESMF_MAXSTR)      :: COMP_NAME

   type(Chem_Registry), pointer    :: chemReg
   type(Aero_GridComp), pointer    :: gcChem      ! Grid Component
   type(Chem_Bundle), pointer      :: w_c         ! Chemical tracer fields     
   integer                         :: nymd, nhms  ! time
   real                            :: cdt         ! chemistry timestep (secs)
   type(MAPL_MetaComp), pointer    :: ggState     ! GEOS Generic State
   type(GOCART_state), pointer     :: state

!  Get my name and set-up traceback handle
!  ---------------------------------------
   call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
   VERIFY_(STATUS)
   Iam = trim(COMP_NAME) // '::' // 'Finalize_'

!  Get my internal MAPL_Generic state
!  -----------------------------------
   call MAPL_GetObjectFromGC ( GC, ggState, __RC__)

   call MAPL_TimerON(ggState, 'TOTAL')
   call MAPL_TimerON(ggState, 'FINALIZE')

!  Get pre-ESMF parameters from gc and clock
!  -----------------------------------------
   call extract_ ( gc, clock, chemReg, gcChem, w_c, nymd, nhms, cdt, STATUS, &
                   state = state )
   VERIFY_(STATUS)

!  Call pre-ESMF version
!  ---------------------
   call Aero_GridCompFinalize ( gcChem, w_c, impChem, expChem, &
                                nymd, nhms, cdt, STATUS )
   VERIFY_(STATUS)

!  Destroy Chem_Bundle
!  -------------------
   call Chem_BundleDestroy ( w_c, STATUS )
   VERIFY_(STATUS)

!  Destroy Chem_Registry
!  ---------------------
   call Chem_RegistryDestroy ( chemReg, STATUS ) 
   VERIFY_(STATUS)

!  Destroy Legacy state
!  --------------------
   deallocate ( state%chemReg, state%gcChem, state%w_c, stat = STATUS )
   VERIFY_(STATUS)

   call MAPL_TimerOff(ggState, 'FINALIZE')
   call MAPL_TimerOff(ggState, 'TOTAL')

!  Finalize GEOS Generic
!  ---------------------
!ALT: do not deallocate "foreign objects"
   call MAPL_GenericFinalize ( gc, impChem, expChem, clock, __RC__ )

   RETURN_(ESMF_SUCCESS)

   end subroutine Finalize_


!.......................................................................

    subroutine extract_ ( gc, clock, chemReg, gcChem, w_c, nymd, nhms, cdt, &
                          rc, state )

    type(ESMF_GridComp), intent(INout)  :: gc
    type(ESMF_Clock), intent(in)     :: clock
    type(Chem_Registry), pointer     :: chemReg
    type(Aero_GridComp), pointer     :: gcChem
    type(Chem_Bundle), pointer       :: w_c
    integer, intent(out)             :: nymd, nhms
    real, intent(out)                :: cdt
    integer, intent(out)             :: rc
    type(MAPL_MetaComp), pointer            :: ggState
    type(GOCART_state), pointer, optional   :: state


    type(GOCART_state), pointer    :: myState

!   ErrLog Variables
!   ----------------
    character(len=ESMF_MAXSTR)      :: IAm
    integer                         :: STATUS
    character(len=ESMF_MAXSTR)      :: COMP_NAME

    type(ESMF_Alarm)                :: ALARM
    type(ESMF_TimeInterval)         :: RingInterval

    type(ESMF_Time)      :: TIME
    type(ESMF_Config)    :: CF
    type(GOCART_Wrap)    :: wrap
    integer              :: IYR, IMM, IDD, IHR, IMN, ISC
    real(ESMF_KIND_R8)   :: dt_r8

!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // '::' // 'extract_'

    rc = 0

!  Get my internal MAPL_Generic state
!  -----------------------------------
   call MAPL_GetObjectFromGC ( GC, ggState, __RC__ )


!   Get my internal state
!   ---------------------
    call ESMF_UserCompGetInternalState(gc, 'GOCART_state', WRAP, STATUS)
    VERIFY_(STATUS)
    myState => wrap%ptr
    if ( present(state) ) then
         state => wrap%ptr
    end if

!   This is likely to be allocated during initialize only
!   -----------------------------------------------------
    if ( .not. associated(myState%chemReg) ) then
         allocate ( myState%chemReg, stat=STATUS )
         VERIFY_(STATUS)
    end if
    if ( .not. associated(myState%gcChem) ) then
         allocate ( myState%gcChem, stat=STATUS )
         VERIFY_(STATUS)
    end if
    if ( .not. associated(myState%w_c) ) then
         allocate ( myState%w_c, stat=STATUS )
         VERIFY_(STATUS)
    end if

    chemReg => myState%chemReg
    gcChem  => myState%gcChem
    w_c     => myState%w_c

!   Get the configuration
!   ---------------------
    call ESMF_GridCompGet ( GC, CONFIG = CF, RC=STATUS )
    VERIFY_(STATUS)

!   Get time step
!   -------------
    call MAPL_Get(ggState, RUNALARM=ALARM, __RC__ )
    call ESMF_AlarmGet(ALARM, ringInterval=RingInterval, __RC__)

    call ESMF_TimeIntervalGet(RingInterval, s_r8=dt_r8, __RC__)
    cdt = real(dt_r8)


!   Need code to extract nymd(20050205), nhms(120000) from clock
!   ------------------------------------------

    call ESMF_ClockGet(CLOCK,currTIME=TIME,rc=STATUS)
    VERIFY_(STATUS)

    call ESMF_TimeGet(TIME ,YY=IYR, MM=IMM, DD=IDD, H=IHR, M=IMN, S=ISC, rc=STATUS)
    VERIFY_(STATUS)

    call MAPL_PackTime(NYMD,IYR,IMM,IDD)
    call MAPL_PackTime(NHMS,IHR,IMN,ISC)

    RETURN_(ESMF_SUCCESS)

   end subroutine extract_


logical function isDataDrivenGC_(gc, rc)

   implicit none

   type(ESMF_GridComp), intent(INout) :: gc
   integer, intent(out)               :: rc
 
!  local 
   character(len=ESMF_MAXSTR)         :: IAm
   integer                            :: STATUS

   integer                            :: i
   character(len=ESMF_MAXSTR)         :: comp_name
   character(len=*), parameter        :: modifier = '.data'

   call ESMF_GridCompGet(gc, name=comp_name, __RC__)   
   i = index(trim(comp_name), trim(modifier), back=.true.)
 
   if (i > 0) then 
       ! lets be strict
       if (comp_name(i:) == modifier) then
           isDataDrivenGC_ = .true.
       else
           isDataDrivenGC_ = .false.
       end if
   else
       isDataDrivenGC_ = .false.
   end if

   RETURN_(ESMF_SUCCESS)

end function isDataDrivenGC_


subroutine aerosol_optics(state, rc)

  implicit none

! Arguments
! ---------
  type(ESMF_State)     :: state
  integer, intent(out) :: rc


! Local
! ---------
  integer                                 :: n_aerosols
  character(len=ESMF_MAXSTR), allocatable :: aerosol_names(:)
  type(ESMF_FieldBundle)                  :: aerosols

  real, dimension(:,:,:), pointer         :: ple
  real, dimension(:,:,:), pointer         :: rh
  real, dimension(:,:,:), pointer         :: var
  real, dimension(:,:,:), pointer         :: q
  real, dimension(:,:,:,:), pointer       :: q_4d

  real, dimension(:,:,:), allocatable     :: dp, f_p

  character(len=ESMF_MAXSTR)              :: fld_name
  type(ESMF_Field)                        :: fld

  real, dimension(:,:,:,:), allocatable   :: ext, ssa, asy  ! (lon:,lat:,lev:,band:)

  integer                                 :: n
  integer                                 :: i1, j1, i2, j2, km

  integer                                 :: band, offset

  integer                                 :: instance

  integer                                 :: STATUS
  character(len=ESMF_MAXSTR)              :: Iam

  integer, parameter                      :: n_bands = 1

  real    :: x
  integer :: i, j, k

  Iam = 'GOCART::aerosol_optics()'


! Mie Table instance/index
! ------------------------
  call ESMF_AttributeGet(state, name='mie_table_instance', value=instance, __RC__)

! Radiation band
! --------------
  band = 0
  call ESMF_AttributeGet(state, name='band_for_aerosol_optics', value=band, __RC__)
  offset = band - n_bands

! Pressure at layer edges 
! ------------------------
  call ESMF_AttributeGet(state, name='air_pressure_for_aerosol_optics', value=fld_name, __RC__)
  call MAPL_GetPointer(state, ple, trim(fld_name), __RC__)

  i1 = lbound(ple, 1); i2 = ubound(ple, 1)
  j1 = lbound(ple, 2); j2 = ubound(ple, 2)
                       km = ubound(ple, 3)

! Relative humidity
! -----------------
  call ESMF_AttributeGet(state, name='relative_humidity_for_aerosol_optics', value=fld_name, __RC__)
  call MAPL_GetPointer(state, rh, trim(fld_name), __RC__)

  i1 = lbound(rh, 1); i2 = ubound(rh, 1)
  j1 = lbound(rh, 2); j2 = ubound(rh, 2)
                      km = ubound(rh, 3)
  
  call ESMF_StateGet(state, 'AEROSOLS', aerosols, __RC__)
  call ESMF_FieldBundleGet(aerosols, fieldCount=n_aerosols, __RC__)

  allocate(aerosol_names(n_aerosols), __STAT__)
 
  call ESMF_FieldBundleGet(aerosols, itemorderflag=ESMF_ITEMORDER_ADDORDER, &
                                     FieldNameList=aerosol_names, __RC__)
 
  allocate(ext(i1:i2,j1:j2,km,n_bands), &
           ssa(i1:i2,j1:j2,km,n_bands), &
           asy(i1:i2,j1:j2,km,n_bands), __STAT__)

  allocate(q_4d(i1:i2,j1:j2,km,n_aerosols), __STAT__)

#if (0)
  allocate(dp(i1:i2,j1:j2,km), f_p(i1:i2,j1:j2,km), __STAT__)

  dp  = ple(:,:,1:km) - ple(:,:,0:km-1)
  f_p = dp / MAPL_GRAV

  do n = 1, n_aerosols
      call ESMF_FieldBundleGet(aerosols, trim(aerosol_names(n)), field=fld, __RC__)
      call ESMF_FieldGet(fld, farrayPtr=q, __RC__)

      q_4d(:,:,:,n) = f_p * q
  end do

  call ESMF_AttributeGet(state, name='mie_table_instance', value=instance, __RC__)
  call mie_(gocartMieTable(instance, aerosol_names, n_bands, offset, q_4d, rh, ext, ssa, asy, __RC__)

  deallocate(dp, f_p, __STAT__)
#else
  do n = 1, n_aerosols
      call ESMF_FieldBundleGet(aerosols, trim(aerosol_names(n)), field=fld, __RC__)
      call ESMF_FieldGet(fld, farrayPtr=q, __RC__)

      do k = 1, km
          do j = j1, j2
              do i = i1, i2
                  x = ((PLE(i,j,k) - PLE(i,j,k-1))*0.01)*(100./MAPL_GRAV)
                  q_4d(i,j,k,n) = x * q(i,j,k)
              end do
          end do
      end do
  end do

  call mie_(gocartMieTable(instance), aerosol_names, n_bands, offset, q_4d, rh, ext, ssa, asy, __RC__)
#endif
  
  call ESMF_AttributeGet(state, name='extinction_in_air_due_to_ambient_aerosol', value=fld_name, __RC__)
  if (fld_name /= '') then 
      call MAPL_GetPointer(state, var, trim(fld_name), __RC__)
      var = ext(:,:,:,1)
  end if

  call ESMF_AttributeGet(state, name='single_scattering_albedo_of_ambient_aerosol', value=fld_name, __RC__)
  if (fld_name /= '') then 
      call MAPL_GetPointer(state, var, trim(fld_name), __RC__)
      var = ssa(:,:,:,1)
  end if

  call ESMF_AttributeGet(state, name='asymmetry_parameter_of_ambient_aerosol', value=fld_name, __RC__)
  if (fld_name /= '') then 
      call MAPL_GetPointer(state, var, trim(fld_name), __RC__)
      var = asy(:,:,:,1)
  end if

  deallocate(aerosol_names, ext, ssa, asy, q_4d, __STAT__)

  RETURN_(ESMF_SUCCESS)

contains 

    subroutine mie_(mie_table, aerosol, nb, offset, q, rh, ext, ssa, asy, rc)
     
     implicit none

     type(Chem_Mie),    intent(inout):: mie_table    ! mie table
     character(len=*),  intent(in )  :: aerosol(:)   ! list of aerosols
     integer,           intent(in )  :: nb           ! number of bands
     integer,           intent(in )  :: offset       ! bands offset 
     real,              intent(in )  :: q(:,:,:,:)   ! aerosol mass mixing ratio, kg kg-1
     real,              intent(in )  :: rh(:,:,:)    ! relative humidity

     real,              intent(out)  :: ext(:,:,:,:) ! extinction
     real,              intent(out)  :: ssa(:,:,:,:) ! SSA
     real,              intent(out)  :: asy(:,:,:,:) ! asymmetry parameter

     integer,           intent(out)  :: rc

     ! local
     integer :: STATUS
     character(len=ESMF_MAXSTR) :: Iam='aerosol_optics::mie_' 

     integer :: l, idx, na

     real(kind=8) :: ext_(size(ext,1),size(ext,2),size(ext,3),size(ext,4))
     real(kind=8) :: ssa_(size(ext,1),size(ext,2),size(ext,3),size(ext,4))
     real(kind=8) :: asy_(size(ext,1),size(ext,2),size(ext,3),size(ext,4))

     na = size(aerosol)

     _ASSERT(na == size(q,4), 'needs informative message')

     ext_ = 0.0d0
     ssa_ = 0.0d0
     asy_ = 0.0d0

     do l = 1, na
        idx = Chem_MieQueryIdx(mie_table, trim(aerosol(l)), __RC__)

        call Chem_MieQueryAllBand4D(mie_table, idx, nb, offset, q(:,:,:,l), rh, ext, ssa, asy, __RC__)

        ext_ = ext_ +          ext     ! total extinction
        ssa_ = ssa_ +     (ssa*ext)    ! total scattering
        asy_ = asy_ + asy*(ssa*ext)    ! sum of (asy * sca)
     end do

     ext = ext_
     ssa = ssa_
     asy = asy_

     RETURN_(ESMF_SUCCESS)

    end subroutine mie_

end subroutine aerosol_optics


subroutine aerosol_activation_properties(state, rc)

  implicit none

! Arguments
! ---------
  type(ESMF_State)     :: state
  integer, intent(out) :: rc


! Local
! ---------
  character(len=ESMF_MAXSTR)      :: mode              ! mode name
  character(len=ESMF_MAXSTR)      :: mode_             ! lowercase mode name 
  type(ESMF_FieldBundle)          :: aerosols          ! field bundle containing the aerosol mass mixing ratios

  real, dimension(:,:,:), pointer :: ple               ! pressure at the edges of model layers
  real, dimension(:,:,:), pointer :: temperature       ! air temperature
  real, dimension(:,:),   pointer :: f_land            ! fraction of land type in a grid cell

  real, dimension(:,:,:), pointer :: f                 ! correction factor for sea salt

  real, dimension(:,:,:), pointer :: q                 ! aerosol mass mixing ratio
  real, dimension(:,:,:), pointer :: q_                ! aerosol mass mixing ratio (temporary)

  real, dimension(:,:,:), pointer :: num               ! number concentration of aerosol particles 
  real, dimension(:,:,:), pointer :: diameter          ! dry size of aerosol
  real, dimension(:,:,:), pointer :: sigma             ! width of aerosol mode
  real, dimension(:,:,:), pointer :: density           ! density of aerosol
  real, dimension(:,:,:), pointer :: hygroscopicity    ! hygroscopicity of aerosol 
  real, dimension(:,:,:), pointer :: f_dust            ! fraction of dust aerosol
  real, dimension(:,:,:), pointer :: f_soot            ! fraction of soot aerosol 
  real, dimension(:,:,:), pointer :: f_organic         ! fraction of organic aerosol

  real                            :: ss_scale          ! sea salt scaling factor
  real                            :: max_clean          ! max mixing ratio before considered polluted
  real                            :: ccn_tuning         ! tunes conversion factors for sulfate
  character(LEN=ESMF_MAXSTR)      :: cld_micro
  
  character(len=ESMF_MAXSTR)      :: fld_name
  type(ESMF_Field)                :: fld

  integer                         :: i1, j1, i2, j2, km

  integer                         :: STATUS
  character(len=ESMF_MAXSTR)      :: Iam

! auxilliary parameters
! ---------------------
  real, parameter :: densSO4 = 1700.0
  real, parameter :: densORG = 1600.0
  real, parameter :: densSS  = 2200.0
  real, parameter :: densDU  = 1700.0
  real, parameter :: densBC  = 1600.0
  real, parameter :: densOC  =  900.0
  real, parameter :: densBRC =  900.0

  real, parameter :: k_SO4   = 0.65
  real, parameter :: k_ORG   = 0.20
  real, parameter :: k_SS    = 1.28
  real, parameter :: k_DU    = 0.0001
  real, parameter :: k_BC    = 0.0001
  real, parameter :: k_OC    = 0.0001
  real, parameter :: k_BRC   = 0.0001

  integer, parameter :: UNKNOWN_AEROSOL_MODE = 2015 


  Iam = 'GOCART::aerosol_activation_properties()'


! Aerosol mode
! ------------
  call ESMF_AttributeGet(state, name='aerosol_mode', value=mode, __RC__)

! Land fraction 
! -------------
  call ESMF_AttributeGet(state, name='fraction_of_land_type', value=fld_name, __RC__)
  call MAPL_GetPointer(state, f_land, trim(fld_name), __RC__)

! Pressure at layer edges 
! ------------------------
  call ESMF_AttributeGet(state, name='air_pressure', value=fld_name, __RC__)
  call MAPL_GetPointer(state, ple, trim(fld_name), __RC__)

! Temperature
! -----------
  call ESMF_AttributeGet(state, name='air_temperature', value=fld_name, __RC__)
  call MAPL_GetPointer(state, temperature, trim(fld_name), __RC__)

  i1 = lbound(temperature, 1); i2 = ubound(temperature, 1)
  j1 = lbound(temperature, 2); j2 = ubound(temperature, 2)
                               km = ubound(temperature, 3)

! Activation activation properties
! --------------------------------
  call ESMF_AttributeGet(state, name='aerosol_number_concentration', value=fld_name, __RC__)
  call MAPL_GetPointer(state, num, trim(fld_name), __RC__)

  call ESMF_AttributeGet(state, name='aerosol_dry_size', value=fld_name, __RC__)
  call MAPL_GetPointer(state, diameter, trim(fld_name), __RC__)

  call ESMF_AttributeGet(state, name='width_of_aerosol_mode', value=fld_name, __RC__)
  call MAPL_GetPointer(state, sigma, trim(fld_name), __RC__)
  
  call ESMF_AttributeGet(state, name='aerosol_density', value=fld_name, __RC__)
  call MAPL_GetPointer(state, density, trim(fld_name), __RC__)

  call ESMF_AttributeGet(state, name='aerosol_hygroscopicity', value=fld_name, __RC__)
  call MAPL_GetPointer(state, hygroscopicity, trim(fld_name), __RC__)

  call ESMF_AttributeGet(state, name='fraction_of_dust_aerosol', value=fld_name, __RC__)
  call MAPL_GetPointer(state, f_dust, trim(fld_name), __RC__)

  call ESMF_AttributeGet(state, name='fraction_of_soot_aerosol', value=fld_name, __RC__)
  call MAPL_GetPointer(state, f_soot, trim(fld_name), __RC__)

  call ESMF_AttributeGet(state, name='fraction_of_organic_aerosol', value=fld_name, __RC__)
  call MAPL_GetPointer(state, f_organic, trim(fld_name), __RC__)

! Sea salt scaling fctor
! ----------------------
  call ESMF_AttributeGet(state, name='seasalt_scaling_factor', value=ss_scale, __RC__)
  call ESMF_AttributeGet(state, name='max_q_clean', value=max_clean, __RC__)
  call ESMF_AttributeGet(state, name='cldmicro', value=cld_micro, __RC__)
  call ESMF_AttributeGet(state, name='ccn_tuning', value=ccn_tuning, __RC__)

! Aerosol mass mixing ratios
! --------------------------
  mode_ = trim(mode)
  mode_ = ESMF_UtilStringLowerCase(mode_, __RC__)

  call ESMF_StateGet(state, 'AEROSOLS', aerosols, __RC__) !GOCART state  

  allocate(q(i1:i2,j1:j2,km), __STAT__)
  q = 0.0

  hygroscopicity = 0.01
  density = 2200.0

  if (index(mode_, 'du00') > 0) then  
      ! dust is mapped one-to-one
      call ESMF_FieldBundleGet(aerosols, trim(mode), field=fld, __RC__)
      call ESMF_FieldGet(fld, farrayPtr=q_, __RC__)         
      q = q_
      hygroscopicity = k_DU
      density = densDU
      
  else if (index(mode_, 'ss00') > 0) then
      ! compute the total mass mixing ratio and impose a tri-modal size distribution
      call ESMF_FieldBundleGet(aerosols, 'ss001', field=fld, __RC__)
      call ESMF_FieldGet(fld, farrayPtr=q_, __RC__)
      q = q + q_

      call ESMF_FieldBundleGet(aerosols, 'ss002', field=fld, __RC__)
      call ESMF_FieldGet(fld, farrayPtr=q_, __RC__)
      q = q + q_

      call ESMF_FieldBundleGet(aerosols, 'ss003', field=fld, __RC__)
      call ESMF_FieldGet(fld, farrayPtr=q_, __RC__)
      q = q + q_

      call ESMF_FieldBundleGet(aerosols, 'ss004', field=fld, __RC__)
      call ESMF_FieldGet(fld, farrayPtr=q_, __RC__)
      q = q + q_

      call ESMF_FieldBundleGet(aerosols, 'ss005', field=fld, __RC__)
      call ESMF_FieldGet(fld, farrayPtr=q_, __RC__)
      q = q + q_

      ! temperature correction over the ocean
      allocate(f(i1:i2,j1:j2, km), __STAT__)
      call ocean_correction_(f, f_land, temperature(i1:i2,j1:j2,km), ss_scale, i1, i2, j1, j2, km)

      ! apply the correction factor
      q = f * q
      deallocate(f, __STAT__)

      hygroscopicity = k_SS
      density = densSS
       
  else if (index(mode_, 'sulforg') > 0) then
      hygroscopicity = 0.0
      density = 0.0

      !internally mixed organics and sulfate
      call ESMF_FieldBundleGet(aerosols, 'SO4', field=fld, __RC__)
      call ESMF_FieldGet(fld, farrayPtr=q_, __RC__)
      q = q + q_
      hygroscopicity = k_SO4*q_ + hygroscopicity
      density = densSO4*q_ + density
    
      call ESMF_FieldBundleGet(aerosols, 'OCphilic', field=fld, __RC__)
      call ESMF_FieldGet(fld, farrayPtr=q_, __RC__)
      q = q + q_
      hygroscopicity = k_ORG*q_ + hygroscopicity
      density = densORG*q_ + density

      where (q > 2.0e-12 .and. hygroscopicity > tiny(0.0))
          hygroscopicity = hygroscopicity / q
          hygroscopicity = max(0.001, hygroscopicity)

          density = density / q
          density = min(max(density, densORG), densSO4)
      elsewhere
          hygroscopicity = k_SO4
          density = densSO4
      end where

      ! required by the aap_(...)
      if(adjustl(cld_micro)/="2MOMENT") then ! maintained for compatibility with the single moment
      
       call ESMF_FieldBundleGet(aerosols, 'SO4', field=fld, __RC__)
       call ESMF_FieldGet(fld, farrayPtr=q_, __RC__)  ! only use the mass of sulfate to make the conversion
     
      end if 

  else if (index(mode_, 'bcphilic') > 0) then
      call ESMF_FieldBundleGet(aerosols, 'BCphilic', field=fld, __RC__)
      call ESMF_FieldGet(fld, farrayPtr=q_, __RC__)
      q = q_
      hygroscopicity = k_BC
      density = densBC
 
  else if (index(mode_, 'ocphilic') > 0) then !this does not activate into droplets, only relevant for ice nuc.
      call ESMF_FieldBundleGet(aerosols, 'OCphilic', field=fld, __RC__)
      call ESMF_FieldGet(fld, farrayPtr=q_, __RC__)         
      q = q_
      hygroscopicity = k_OC
      density = densOC

  else
      __raise__(UNKNOWN_AEROSOL_MODE, "Unknown aerosol mode used in the GOCART aerosol activation properties method: "//trim(mode))

  end if

! Obtain aerosol activation properties of this aerosol mode
! ---------------------------------------------------------
  call aap_(mode,               &
            q,                  &
            num,                &
            diameter,           &
            sigma,              &        
            f_dust,             &
            f_soot,             &
            f_organic,          &
            density,            &
            q_,                 &
            i1, i2, j1, j2, km, &
            __RC__)

  deallocate(q, __STAT__)

  RETURN_(ESMF_SUCCESS)

contains

    subroutine aap_(mode, q, num, diameter, sigma, f_dust, f_soot, f_organic, dens_, q_, &
                    i1, i2, j1, j2, km, rc)
     
     implicit none

     integer, intent(in) :: i1, i2                                  ! dimension bounds
     integer, intent(in) :: j1, j2                                  ! ... // ..
     integer, intent(in) :: km                                      ! ... // ..

     character(len=*),  intent(in )               :: mode           ! name of aerosol mode
     real, intent(in),  dimension(i1:i2,j1:j2,km) :: q              ! aerosol mass mixing ratio, kg kg-1
     real, intent(in),  dimension(i1:i2,j1:j2,km) :: q_             ! auxiliary mass
     real, intent(in),  dimension(i1:i2,j1:j2,km) :: dens_          ! density


     real, intent(out), dimension(i1:i2,j1:j2,km) :: num            ! number concentration of aerosol particles 
     real, intent(out), dimension(i1:i2,j1:j2,km) :: diameter       ! dry size of aerosol
     real, intent(out), dimension(i1:i2,j1:j2,km) :: sigma          ! width of aerosol mode  
     real, intent(out), dimension(i1:i2,j1:j2,km) :: f_dust         ! fraction of dust aerosol
     real, intent(out), dimension(i1:i2,j1:j2,km) :: f_soot         ! fraction of soot aerosol 
     real, intent(out), dimension(i1:i2,j1:j2,km) :: f_organic      ! fraction of organic aerosol

     integer, intent(out) :: rc                                     ! return code

     ! local
     integer :: STATUS
     character(len=ESMF_MAXSTR) :: mode_
     character(len=ESMF_MAXSTR) :: Iam = 'GOCART::aerosol_activation_properties::aap_()'

     integer, parameter :: UNKNOWN_AEROSOL_MODE = 2015

     integer            :: kinx
     real               :: fmassaux, fmassclean
     real, dimension(3) :: TPI, DPGI, SIGI
     real, dimension(3) :: TPIclean, DPGIclean, SIGIclean
     real, dimension(i1:i2,j1:j2,km) :: qaux
      !real, parameter    :: max_clean = 5.0e-7  !max mixing ratio before considered polluted
 

 

     mode_ = trim(mode)
     mode_ = ESMF_UtilStringLowerCase(mode_, __RC__)

     num       = 0.0
     diameter  = 1.0e-9
     sigma     = log(2.0)
     f_dust    = 0.0
     f_soot    = 0.0
     f_organic = 0.0

      if(adjustl(cld_micro)=="2MOMENT") then
        qaux=q !this corrects a bug
      else
        qaux  =  q_ !keep it to get zero diff with the single moment
        max_clean = 5.0e-7
        ccn_tuning = 1.0
      end if


     if (index(mode_, 'ss00') > 0) then
       if(adjustl(cld_micro)=="2MOMENT") then
         TPI  (1) = 230e6          ! num fraction (reduced 091015)        
       else       
         TPI  (1) = 100e6          ! num fraction (reduced 091015)                   
       end if
     
         DPGI (1) = 0.02e-6        ! modal diameter (m)
         SIGI (1) = log(1.6)       ! geometric dispersion (sigma_g)
         ! accumulation
         TPI  (2) = 60.0e6         ! total concentration (# m-3)
         DPGI (2) = 0.071e-6       ! modal diameter (m)
         SIGI (2) = log(2.0)       ! geometric dispersion (sigma_g)
         ! coarse
         TPI  (3) = 3.1e6          ! total concentration (# m-3)
         DPGI (3) = 0.62e-6        ! modal diameter (m)
         SIGI (3) = log(2.7)       ! geometric dispersion (sigma_g)
          
         fmassaux = 0.0
         do kinx = 1, 3
             fmassaux = (TPI(kinx)*densSS*MAPL_PI*exp(4.5*SIGI(kinx)*SIGI(kinx))*DPGI(kinx)*DPGI(kinx)*DPGI(kinx))/6.0 + fmassaux
         end do
     end if

     if (index(mode_, 'sulforg0') > 0) then      
         TPI  (1) = 1.06e11        ! num fraction
         DPGI (1) = .014e-6        ! modal diameter (m)
         SIGI (1) = log(1.8)       ! geometric dispersion (sigma_g)
         ! accumulation
         TPI  (2) = 3.2e10         ! total concentration (# m-3)
         DPGI (2) = 0.054e-6       ! modal diameter (m)
         SIGI (2) = log(2.16)      ! geometric dispersion (sigma_g)
         !coarse
         TPI  (3) = 5.4e6          ! total concentration (# m-3)
         DPGI (3) = 0.86e-6        ! modal diameter (m)
         SIGI (3) = log(2.21)      ! geometric dispersion (sigma_g)

         fmassaux = 0.0
         do kinx = 1, 3
             ! density is multiplied below since this is a case of a 3-d field
             fmassaux = (TPI(kinx)*MAPL_PI*exp(4.5*SIGI(kinx)*SIGI(kinx))*DPGI(kinx)*DPGI(kinx)*DPGI(kinx))/6.0 + fmassaux
         end do

         ! clean continental polluted plus org
         ! fine
         TPIclean  (1) = 1.0e9      ! total concentration (# m-3)
         DPGIclean (1) = 0.016e-6   ! modal diameter (m)
         SIGIclean (1) = log(1.6)   ! geometric dispersion (sigma_g)      
         ! accumulation
         TPIclean  (2) = 8.0e8      ! total concentration (# m-3)
         DPGIclean (2) = 0.067e-6   ! modal diameter (m)
         SIGIclean (2) = log(2.1)   ! geometric dispersion (sigma_g) 
         !Coarse
         TPIclean  (3) = 2.0e6      ! total concentration (# m-3)
         DPGIclean (3) = 0.93e-6    ! modal diameter (m)
         SIGIclean (3) = log(2.2)   ! geometric dispersion (sigma_g)
              
         fmassclean= 0.0
         do kinx = 1, 3
             fmassclean = (TPIclean(kinx)*MAPL_PI*exp(4.5*SIGIclean(kinx)*SIGIclean(kinx))*DPGIclean(kinx)*DPGIclean(kinx)*DPGIclean(kinx))/6.0 + fmassclean  !
         end do 
     end if 
           


     select case(mode_)

     case ('du001')
         sigma    = log(1.8)
         f_dust   = 1.0
         diameter = 1.46e-6
         num      = q / ((MAPL_PI/6.0) * densDU * diameter*diameter*diameter * exp(4.5*sigma*sigma))

     case ('du002')
         sigma    = log(1.8)
         f_dust   = 1.0
         diameter = 2.80e-6
         num      = q / ((MAPL_PI/6.0) * densDU * diameter*diameter*diameter * exp(4.5*sigma*sigma))

     case ('du003')
         sigma    = log(1.8)
         f_dust   = 1.0
         diameter = 4.80e-6
         num      = q / ((MAPL_PI/6.0) * densDU * diameter*diameter*diameter * exp(4.5*sigma*sigma))

     case ('du004')
         sigma    = log(1.8)
         f_dust   = 1.0
         diameter = 9.0e-6
         num      = q / ((MAPL_PI/6.0) * densDU * diameter*diameter*diameter * exp(4.5*sigma*sigma))

     case ('du005')
         sigma    = log(1.8)
         f_dust   = 1.0
         diameter = 16.0e-6
         num      = q / ((MAPL_PI/6.0) * densDU * diameter*diameter*diameter * exp(4.5*sigma*sigma))

     case ('ss001')
         sigma    = SIGI(1)
         diameter = DPGI(1)
         num      = TPI(1) * q / fmassaux

     case ('ss002')
         sigma    = SIGI(2)
         diameter = DPGI(2)
         num      = TPI(2) * q / fmassaux

     case ('ss003')
         sigma    = SIGI(3)
         diameter = DPGI(3)
         num      = TPI(3) * q / fmassaux
    
     case ('sulforg01')  !different distributions for clean and polluted environments
         where (q > max_clean)          
             sigma    = SIGI(1)
             diameter = DPGI(1)
             num      = TPI(1) * qaux*ccn_tuning / (dens_*fmassaux)             ! only sulfate  mass
         elsewhere 
             sigma    = SIGIclean(1)
             diameter = DPGIclean(1)
             num      = TPIclean(1) * qaux*ccn_tuning / (dens_*fmassclean)      ! only sulfate 
         end where 

     case ('sulforg02')
         where (q > max_clean)
             sigma    = SIGI(2)
             diameter = DPGI(2)
             num      = TPI(2) * qaux*ccn_tuning / (dens_*fmassaux)            ! only sulfate mass
         elsewhere
             sigma    = SIGIclean(2)
             diameter = DPGIclean(2)        
             num      = TPIclean(2) * qaux*ccn_tuning / (dens_*fmassclean)     ! only sulfate
         end where 
   
     case ('sulforg03')
         where (q > max_clean)
             sigma    = SIGI(3)
             diameter = DPGI(3)        
             num      = TPI(3) * qaux*ccn_tuning / (dens_*fmassaux)           ! only sulfate mass
         elsewhere 
             sigma    = SIGIclean(3)
             diameter = DPGIclean(3)
             num      = TPIclean(3) * qaux*ccn_tuning / (dens_*fmassclean)    ! only sulfate
         end where 

     case ('bcphilic')
         sigma    = log(2.0)
         f_soot   = 1.0
         diameter = 0.0118*2e-6
         num = q / ((MAPL_PI/6.0) * densBC * diameter*diameter*diameter * exp(4.5*sigma*sigma))

     case ('ocphilic')
         sigma     = log(2.2)
         f_organic = 1.0
         diameter  = 0.0212*2.0e-6
         num = q / ((MAPL_PI/6.0) * densOrg * diameter*diameter*diameter * exp(4.5*sigma*sigma))

     case default
         __raise__(UNKNOWN_AEROSOL_MODE,"Unknown aerosol mode used in the GOCART aerosol activation properties method: "//trim(mode))

     end select


     RETURN_(ESMF_SUCCESS)

    end subroutine aap_


    subroutine ocean_correction_(f, f_land, t_air_sfc, ss_scale, i1, i2, j1, j2, km)
     
     implicit none

     integer, intent(in) :: i1, i2                               ! dimension bounds
     integer, intent(in) :: j1, j2                               ! ... // ..
     integer, intent(in) :: km                                   ! ... // ..

     real, intent(in ), dimension(i1:i2,j1:j2) :: f_land         ! fraction of land
     real, intent(in ), dimension(i1:i2,j1:j2) :: t_air_sfc      ! air temperature in the surface model layer
     real, intent(in )                         :: ss_scale       ! scaling factor for sea salt at low T

     real, intent(out), dimension(i1:i2,j1:j2, km) :: f          ! correction factor

     ! local
     integer :: i, j
     real    :: usurf

     f = 1.0

     do j = j1, j2
         do i = i1, i2
             if (f_land(i,j) < 0.1) then  !ocean

                 if(adjustl(cld_micro) .ne."2MOMENT") then
                    usurf = max(min((t_air_sfc(i,j) - 285.0) / 2.0, 10.0), -10.0) !smooth transition around some T value		   		   	      	      
                 else
                    usurf = max(min((t_air_sfc(i,j) - 285.0) / 2.0, 30.0), -30.0) !smooth transition around some T value
                 end if 
                 usurf = min(ss_scale / (1.0 + exp(usurf)), 20.0)

                 f(i,j,:) = (1.0 + usurf)
             end if
         end do
     end do
          
    end subroutine ocean_correction_

end subroutine aerosol_activation_properties



subroutine MAMnet(state, rc)

    implicit none

    ! Arguments
    ! ---------
    type(ESMF_State)     :: state
    integer, intent(out) :: rc

    ! Local
    ! ---------
    type(ESMF_FieldBundle)          :: aerosols          ! field bundle containing the aerosol mass mixing ratios

    real, dimension(:,:,:), pointer :: ple               ! pressure at the edges of model layers
    real, dimension(:,:,:), pointer :: temperature       ! air temperature
    real, dimension(:,:,:), pointer :: airden            ! air density
    real, dimension(:,:,:), pointer :: q, q_             ! aerosol mass mixing ratio

    real, dimension(:,:,:), pointer :: num               ! number concentration of aerosol particles 
    real, dimension(:,:,:), pointer :: diameter          ! dry size of aerosol
    real, dimension(:,:,:), pointer :: sigma             ! width of aerosol mode
    real, dimension(:,:,:), pointer :: density           ! density of aerosol
    real, dimension(:,:,:), pointer :: hygroscopicity    ! hygroscopicity of aerosol 
    real, dimension(:,:,:), pointer :: f_dust            ! fraction of dust aerosol
    real, dimension(:,:,:), pointer :: f_soot            ! fraction of soot aerosol 
    real, dimension(:,:,:), pointer :: f_organic         ! fraction of organic aerosol
    
    real, dimension(:,:,:, :),  pointer :: num_4d               ! number concentration of aerosol particles 
    real, dimension(:,:,:, :), pointer :: diameter_4d          ! dry size of aerosol
    real, dimension(:,:,:, :), pointer :: sigma_4d             ! width of aerosol mode
    real, dimension(:,:,:, :), pointer :: density_4d           ! density of aerosol
    real, dimension(:,:,:, :), pointer :: hygroscopicity_4d    ! hygroscopicity of aerosol 
    real, dimension(:,:,:, :), pointer :: f_dust_4d            ! fraction of dust aerosol
    real, dimension(:,:,:, :), pointer :: f_soot_4d            ! fraction of soot aerosol 
    real, dimension(:,:,:, :), pointer :: f_organic_4d         ! fraction of organic aerosol
    
    real, dimension(:,:,:,:), pointer :: mass         ! MAM species 
    real, dimension(:,:,:), pointer :: modal_mass         ! MAM species
    real, dimension(:,:,:), pointer :: modal_number        ! MAM modal numbers

    character(len=ESMF_MAXSTR)      :: fld_name
    type(ESMF_Field)                :: fld

    integer                         :: i1, j1, i2, j2, km

    integer                         :: STATUS
    character(len=ESMF_MAXSTR)      :: Iam

    ! auxilliary parameters
    ! ---------------------
    real, parameter :: dens_SO4 = 1700.0
    real, parameter :: dens_ORG = 1600.0
    real, parameter :: dens_SS  = 2200.0
    real, parameter :: dens_DU  = 1700.0
    real, parameter :: dens_BC  = 1600.0

    real, parameter :: k_SO4   = 0.65
    real, parameter :: k_ORG   = 0.20
    real, parameter :: k_SS    = 1.28
    real, parameter :: k_DU    = 0.01
    real, parameter :: k_BC    = 0.0001
    
    real, parameter :: sig_acc   = 1.8
    real, parameter :: sig_ait   = 1.6
    real, parameter :: sig_cdu   = 1.8
    real, parameter :: sig_css   = 2.0
    real, parameter :: sig_fdu   = 1.8
    real, parameter :: sig_fss   = 2.0
    real, parameter :: sig_pcm   = 1.6

    integer, parameter :: UNKNOWN_AEROSOL_MODE = 2015 
    
    character(len=ESMF_MAXSTR)      :: mam_species(24)


    real, parameter, dimension(7) :: means =  (/-14.35511, -14.725009, -14.91018, -15.612967, -15.061058, 2.3847961, -1.353214/)
    integer :: samples, n_species, n_modes, im, jm, sx 
    
    ! Neural Network parameterization of aerosol number concentration  (Breen et al. 2022)
    ! Two neural nets work in series:
    ! mass_net (Gocart total mass, T, P [samples, 7*72] to MAM mass species [samples, 24*72])    
    ! number_net (MAM species [samples, 24*72] to MAM modal number conc. [samples, 7*72])   
    !mass_net input is log10 transformed, then removing the precalculated mean 
    !mass net output is number-like using precalculated modal diameters and feeds to number net
    ! number_net output is log10-transformed
     
    type (network_type) :: mass_net
    type (network_type) :: number_net
    
    real, allocatable, dimension(:,:) :: massnet_in
    real, allocatable, dimension(:,:) :: massnet_out
    real, allocatable, dimension(:,:) :: numnet_out
    real, allocatable, dimension(:,:) :: q_reshaped
    
    real :: scale_factor
    

    Iam = 'MAMnet'
    
    mam_species = (/'SU_A_ACC', 'SOA_A_ACC', 'SS_A_ACC', 'POM_A_ACC', 'BC_A_ACC', 'AMM_A_ACC', &
                   'SU_A_AIT', 'SOA_A_AIT', 'SS_A_AIT', 'AMM_A_AIT', &
                   'SU_A_CDU', 'DU_A_CDU', 'AMM_A_CDU', &
                   'SU_A_CSS', 'SS_A_CSS', 'AMM_A_CSS', &
                   'SU_A_FDU', 'DU_A_FDU', 'AMM_A_FDU', &
                   'SU_A_FSS', 'SS_A_FSS', 'AMM_A_FSS', &
                   'POM_A_PCM', 'BC_A_PCM'/)
                 
    n_species  = size(mam_species)             
     
    !initialize the neural networks
    
    call mass_net%load('./massnet_weights.rc')
    call number_net%load('./numbernet_weights.rc') 


    ! Temperature
    ! -----------
    call ESMF_AttributeGet(state, name='air_temperature', value=fld_name, __RC__)
    call MAPL_GetPointer(state, temperature, trim(fld_name), __RC__)

    ! dry air density 
    ! ------------------------
    call ESMF_AttributeGet(state, name='air_density', value=fld_name, __RC__)
    call MAPL_GetPointer(state, airden, trim(fld_name), __RC__)
    
    ! Aerosol properties
    ! --------------------------------    

    call ESMF_AttributeGet(state, name='aerosol_dry_size', value=fld_name, __RC__)
    call MAPL_GetPointer(state, diameter, trim(fld_name), __RC__)
    
    call ESMF_AttributeGet(state, name='aerosol_number_concentration', value=fld_name, __RC__)
    call MAPL_GetPointer(state, num, trim(fld_name), __RC__)

    call ESMF_AttributeGet(state, name='width_of_aerosol_mode', value=fld_name, __RC__)
    call MAPL_GetPointer(state, sigma, trim(fld_name), __RC__)

    call ESMF_AttributeGet(state, name='aerosol_density', value=fld_name, __RC__)
    call MAPL_GetPointer(state, density, trim(fld_name), __RC__)

    call ESMF_AttributeGet(state, name='aerosol_hygroscopicity', value=fld_name, __RC__)
    call MAPL_GetPointer(state, hygroscopicity, trim(fld_name), __RC__)

    call ESMF_AttributeGet(state, name='fraction_of_dust_aerosol', value=fld_name, __RC__)
    call MAPL_GetPointer(state, f_dust, trim(fld_name), __RC__)

    call ESMF_AttributeGet(state, name='fraction_of_soot_aerosol', value=fld_name, __RC__)
    call MAPL_GetPointer(state, f_soot, trim(fld_name), __RC__)

    call ESMF_AttributeGet(state, name='fraction_of_organic_aerosol', value=fld_name, __RC__)
    call MAPL_GetPointer(state, f_organic, trim(fld_name), __RC__)


    i1 = lbound(temperature, 1); i2 = ubound(temperature, 1)
    j1 = lbound(temperature, 2); j2 = ubound(temperature, 2)
    km = ubound(temperature, 3)
    
    !if (km .ne. 72) __raise__(MAMnet only works for 72 levels) !for now the NNs always expect 72 levels
    
    call ESMF_StateGet(state, 'AEROSOLS', aerosols, __RC__) !GOCART state  

    allocate(q(i1:i2,j1:j2,km), __STAT__)
    allocate(q_(i1:i2,j1:j2,km), __STAT__)
    q = 0.0 
    q_ =0.0        
    
    
    !!!!!! make the feature set for mass_net
    !!!!!! sulfate

    call ESMF_FieldBundleGet(aerosols, 'SO4', field=fld, __RC__)
    call ESMF_FieldGet(fld, farrayPtr=q, __RC__)
    
    
    im =  size(q, 1)
    jm =  size(q, 2)
    km =  72 !size(q, 3)
    
    samples = im*jm! each column is an example
           
    allocate(q_reshaped(samples, km))
    km =  72*7 
    allocate(massnet_in(samples, km)) 
    allocate(numnet_out(samples, km))
    km =  72*24
    allocate(massnet_out(samples, km))   

    km  = 72
    
    !allocate (mass(im, jm, km, n_species))

    

    i1 =  1
    i2 = 72 
    where(q> 0.0)
        q =  log10(q) - means(1)       
    end where 
        q_reshaped  = reshape(q, shape(q_reshaped))    
    
    massnet_in(:, i1:i2) =  q_reshaped 
    
    
   ! print *, 'massnet_in'
   ! print *, massnet_in
    !!!!!!! sea salt
    
    q =  0.0

    call ESMF_FieldBundleGet(aerosols, 'ss001', field=fld, __RC__)
    call ESMF_FieldGet(fld, farrayPtr=q_, __RC__)
    q = q + q_

    call ESMF_FieldBundleGet(aerosols, 'ss002', field=fld, __RC__)
    call ESMF_FieldGet(fld, farrayPtr=q_, __RC__)
    q = q + q_

    call ESMF_FieldBundleGet(aerosols, 'ss003', field=fld, __RC__)
    call ESMF_FieldGet(fld, farrayPtr=q_, __RC__)
    q = q + q_

    call ESMF_FieldBundleGet(aerosols, 'ss004', field=fld, __RC__)
    call ESMF_FieldGet(fld, farrayPtr=q_, __RC__)
    q = q + q_

    call ESMF_FieldBundleGet(aerosols, 'ss005', field=fld, __RC__)
    call ESMF_FieldGet(fld, farrayPtr=q_, __RC__)
    q = q + q_

    i1 = 73  
    i2 = 144 
    where(q> 0.0)
        q =  log10(q) - means(2)        
    end where 
    q_reshaped  = reshape(q, shape(q_reshaped))    
    massnet_in(:, i1:i2) =  q_reshaped 
     
    !!organics
    q = 0.0  
    call ESMF_FieldBundleGet(aerosols, 'OCphilic', field=fld, __RC__)
    call ESMF_FieldGet(fld, farrayPtr=q_, __RC__)
    q = q + q_
    
    call ESMF_FieldBundleGet(aerosols, 'OCphobic', field=fld, __RC__)
    call ESMF_FieldGet(fld, farrayPtr=q_, __RC__)
    q = q + q_
  
  
    i1 = 145  
    i2 = 216 
    where(q > 0.0)
        q =  log10(q) - means(3)        
    end where 
    q_reshaped  = reshape(q, shape(q_reshaped))    
    massnet_in(:, i1:i2) =  q_reshaped 
    
    
    !!black carbon
    q = 0.0  
    call ESMF_FieldBundleGet(aerosols, 'BCphilic', field=fld, __RC__)
    call ESMF_FieldGet(fld, farrayPtr=q_, __RC__)
    q = q + q_
    
    call ESMF_FieldBundleGet(aerosols, 'BCphobic', field=fld, __RC__)
    call ESMF_FieldGet(fld, farrayPtr=q_, __RC__)
    q = q + q_
  
  
    i1 = 217  
    i2 = 288 
    where(q > 0.0)
        q =  log10(q) - means(4)        
    end where
    q_reshaped  = reshape(q, shape(q_reshaped))     
    massnet_in(:, i1:i2) =  q_reshaped 
    
    !!!!!!! dust
    
    q =  0.0

    call ESMF_FieldBundleGet(aerosols, 'du001', field=fld, __RC__)
    call ESMF_FieldGet(fld, farrayPtr=q_, __RC__)
    q = q + q_

    call ESMF_FieldBundleGet(aerosols, 'du002', field=fld, __RC__)
    call ESMF_FieldGet(fld, farrayPtr=q_, __RC__)
    q = q + q_

    call ESMF_FieldBundleGet(aerosols, 'du003', field=fld, __RC__)
    call ESMF_FieldGet(fld, farrayPtr=q_, __RC__)
    q = q + q_

    call ESMF_FieldBundleGet(aerosols, 'du004', field=fld, __RC__)
    call ESMF_FieldGet(fld, farrayPtr=q_, __RC__)
    q = q + q_

    call ESMF_FieldBundleGet(aerosols, 'du005', field=fld, __RC__)
    call ESMF_FieldGet(fld, farrayPtr=q_, __RC__)
    q = q + q_

    i1 = 289  
    i2 = 360 
    where(q> 0.0)
        q =  log10(q) - means(5)        
    end where 
    q_reshaped  = reshape(q, shape(q_reshaped))    
    massnet_in(:, i1:i2) =  q_reshaped 
    
    !! temperature 
    
    i1 = 361  
    i2 = 432
     
    q =  0.0 
    where(temperature > 0.0)
        q =  log10(temperature) - means(6)        
    end where  
    q_reshaped  = reshape(q, shape(q_reshaped))   
    massnet_in(:, i1:i2) =  q_reshaped 
    
    !! air density 
    
    i1 = 433  
    i2 = 504
     
    q =  0.0 
    where(airden > 0.0)
        q =  log10(airden) - means(7)        
    end where
    q_reshaped  = reshape(q, shape(q_reshaped))     
    massnet_in(:, i1:i2) =  q_reshaped
    
    !!!!!! run the Neural Networks
    do sx  = 1, samples
    	!******call mass_net
    	massnet_out(sx, :) = mass_net%output(massnet_in(sx, :)) ![samples, 7*24]
    	!****** call number_net    
    	numnet_out(sx, :) = number_net%output(massnet_out(sx, :)) ![samples, 7*72]
    end do    
    
    
    
    n_modes =  7 
    km = 72
    allocate(num_4d(im, jm, km, n_modes))
    allocate(sigma_4d(im, jm, km, n_modes)) 
    allocate(density_4d(im, jm, km, n_modes)) 
    allocate(hygroscopicity_4d(im, jm, km, n_modes))
    allocate( diameter_4d(im, jm, km, n_modes))
    allocate(f_dust_4d(im, jm, km, n_modes))
    allocate(f_soot_4d(im, jm, km, n_modes))
    allocate(f_organic_4d(im, jm, km, n_modes))
    allocate(modal_mass(im, jm, km))
    allocate(modal_number(im, jm, km))
     
    allocate(mass(im, jm, km, n_species))        
    !!!map numbernet_out to each mode
    mass = reshape(massnet_out,  (/im, jm, km, n_species/)) ![im, jm, 72, 24]
    num_4d =  reshape(numnet_out, (/im, jm, km, n_modes/)) ![im, jm, 72, 7]    
    num_4d =  0.0
    !return number concentration to #/Kg
    where ((num_4d .le. 18.0) .and. (num_4d .gt. -3.0))
     num_4d =  10.0**num_4d 
    end where 
    
    !!!!!!return mixing ratios to Kg/Kg and fill up the other properties
   
    mass  = 10.0**mass !remove log scaling
    
    !set default values
    sigma_4d =  num_4d*0.0 + sig_acc
    density_4d =  num_4d*0.0 + dens_SO4
    hygroscopicity_4d =  num_4d*0.0  + k_SO4
    diameter_4d =  num_4d*0.0 + 1.0e-9
    f_dust_4d  =  num_4d*0.0 
    f_soot_4d =  num_4d*0.0
    f_organic_4d =  num_4d*0.0
 
    !!!!!!!!!!!! mode 1: accumulation, species 1 to 6
    !'SU_A_ACC', 'SOA_A_ACC', 'SS_A_ACC', 'POM_A_ACC', 'BC_A_ACC', 'AMM_A_ACC'
    scale_factor  = (((0.26+0.056)/2.0)**3.0)* 1600.0 * 8.0 * 1.0e-18 ! scaling using during training 
    
    i1 = 1
    i2 = 6
    mass(:, :, :, i1:i2)  =  mass(:, :, :, i1:i2)*scale_factor
    modal_mass =  sum(mass(:, :, :, i1:i2), 4)
    modal_number = num_4d(:, :, :,1)
    
    sigma_4d(:, :, :, 1) = sig_acc
        
    where ((modal_mass .gt. 0.0 ) .and.  (modal_number .gt. 0.0 )) 
       hygroscopicity_4d(:, :, :, 1) =  (mass(:, :, :, 1)*k_SO4 + mass(:, :, :, 2)*k_ORG + mass(:, :, :, 3)*k_SS  &
      					 + mass(:, :, :, 4)*k_ORG + mass(:, :, :, 5)*k_BC +  mass(:, :, :, 6)*k_SO4)/modal_mass !approximated by mass fraction
       diameter_4d(:, :, :, 1) =   ((modal_mass/modal_number)**(1./3.))*exp(-3.0*log(sig_acc)*log(sig_acc))          
       density_4d(:, :, :, 1) =  (mass(:, :, :, 1)*dens_SO4 + mass(:, :, :, 2)*dens_ORG + mass(:, :, :, 3)*dens_SS  &
       							+ mass(:, :, :, 4)*dens_ORG + mass(:, :, :, 5)*dens_BC +  mass(:, :, :, 6)*dens_SO4) /modal_mass 
       f_organic_4d(:, :, :, 1) =  (mass(:, :, :, 2) + mass(:, :, :, 4)) /modal_mass
       f_soot_4d(:, :, :, 1) =  mass(:, :, :, 5) /modal_mass
    end where
    
    diameter_4d(:, :, :, 1)  = min(max(diameter_4d(:, :, :, 1), 0.056e-6), 0.26e-6) 
    
    
    !!!!!!!!!!!!! mode 2: aitken, species 7 to 10
    ! 'SU_A_AIT', 'SOA_A_AIT', 'SS_A_AIT', 'AMM_A_AIT',
    scale_factor  = (((0.052+0.015)/2.0)**3.0)* 1600.0 * 8.0 * 1.0e-18 ! scaling using durinhg training    
    i1 = 7
    i2 = 10
    mass(:, :, :, i1:i2)  =  mass(:, :, :, i1:i2)*scale_factor
    modal_mass =  sum(mass(:, :, :, i1:i2), 4)
    modal_number = num_4d(:, :, :, 2)
    
    sigma_4d(:, :, :, 2) = sig_ait
    
    where ((modal_mass .gt. 0.0 ) .and.  (modal_number .gt. 0.0 )) 
       hygroscopicity_4d(:, :, :, 2) =  (mass(:, :, :, 7)*k_SO4 + mass(:, :, :, 8)*k_ORG + mass(:, :, :, 9)*k_SS  &
      								 +   mass(:, :, :, 10*k_SO4))/modal_mass    
       diameter_4d(:, :, :, 2) =   ((modal_mass/modal_number)**(1./3.))*exp(-3.0*log(sig_ait)*log(sig_ait)) 
       density_4d(:, :, :, 2) =  (mass(:, :, :, 7)*dens_SO4 + mass(:, :, :, 8)*dens_ORG + mass(:, :, :, 9)*dens_SS  &
      							 +  mass(:, :, :, 10)*dens_SO4) /modal_mass 
       f_organic_4d(:, :, :, 2) =  mass(:, :, :, 8) /modal_mass
    end where
    diameter_4d(:, :, :, 2)  = min(max(diameter_4d(:, :, :, 2), 0.015e-6), 0.052e-6) 
 
   !!!!!!!!!!!!!!!!! mode 3: coarse dust, species 11 to 13
    !'SU_A_CDU', 'DU_A_CDU', 'AMM_A_CDU',
    scale_factor  = (((2.75+0.59)/2.0)**3.0)* 1600.0 * 8.0 * 1.0e-18 ! scaling using durinhg training =    
    i1 = 11
    i2 = 13
    mass(:, :, :, i1:i2)  =  mass(:, :, :, i1:i2)*scale_factor
    modal_mass =  sum(mass(:, :, :, i1:i2), 4)
    modal_number = num_4d(:, :, :, 3)
    
    sigma_4d(:, :, :, 3) = sig_cdu
    
    where ((modal_mass .gt. 0.0 ) .and.  (modal_number .gt. 0.0 )) 
       hygroscopicity_4d(:, :, :, 3) =  (mass(:, :, :, 11)*k_SO4 + mass(:, :, :, 12)*k_DU + mass(:, :, :, 13)*k_SO4)/modal_mass  
       diameter_4d(:, :, :, 3) =   ((modal_mass/modal_number)**(1./3.))*exp(-3.0*log(sig_cdu)*log(sig_cdu))    
       density_4d(:, :, :, 3) =  (mass(:, :, :, 11)*dens_SO4 + mass(:, :, :, 12)*dens_DU +  mass(:, :, :, 13)*dens_SO4) /modal_mass    
       f_dust_4d(:, :, :, 3) =  mass(:, :, :, 12) /modal_mass
    end where
    
    diameter_4d(:, :, :, 3)  = min(max(diameter_4d(:, :, :, 3), 0.59e-6), 2.75e-6) 
     
    !!!!!!! mode 4: coarse sea salt, species 14 to 16
    !'SU_A_CSS', 'SS_A_CSS', 'AMM_A_CSS',,
    scale_factor  = (((3.7+0.63)/2.0)**3.0)* 1600.0 * 8.0 * 1.0e-18 ! scaling using durinhg training    
    i1 = 14
    i2 = 16
    mass(:, :, :, i1:i2)  =  mass(:, :, :, i1:i2)*scale_factor
    modal_mass =  sum(mass(:, :, :, i1:i2), 4)
    modal_number = num_4d(:, :, :, 4)
    
    sigma_4d(:, :, :, 4) = sig_css
    
    where ((modal_mass .gt. 0.0 ) .and.  (modal_number .gt. 0.0 )) 
       hygroscopicity_4d(:, :, :, 4) =  (mass(:, :, :, 14)*k_SO4 + mass(:, :, :, 15)*k_SS + mass(:, :, :, 16)*k_SO4)/modal_mass 
       diameter_4d(:, :, :, 4) =   ((modal_mass/modal_number)**(1./3.))*exp(-3.0*log(sig_css)*log(sig_css))  
       density_4d(:, :, :, 4) =  (mass(:, :, :, 14)*dens_SO4 + mass(:, :, :, 15)*dens_SS +  mass(:, :, :, 16)*dens_SO4) /modal_mass 
    end where     
    diameter_4d(:, :, :, 4) =   min(max(diameter_4d(:, :, :, 3), 0.63e-6), 3.7e-6)  
    
    
    !!!!!!!!!!!! mode 5: fine dust, species 17 to 19
    !'SU_A_FDU', 'DU_A_FDU', 'AMM_A_FDU',
    scale_factor  = (((0.62+0.14)/2.0)**3.0)* 1600.0 * 8.0 * 1.0e-18 ! scaling using durinhg training     
    i1 = 17
    i2 = 19
    mass(:, :, :, i1:i2)  =  mass(:, :, :, i1:i2)*scale_factor
    modal_mass =  sum(mass(:, :, :, i1:i2), 4)
    modal_number = num_4d(:, :, :, 5)
    
    sigma_4d(:, :, :, 5) = sig_fdu
    
    where ((modal_mass .gt. 0.0 ) .and.  (modal_number .gt. 0.0 )) 
       hygroscopicity_4d(:, :, :, 5) =  (mass(:, :, :, 17)*k_SO4 + mass(:, :, :, 18)*k_DU + mass(:, :, :, 19)*k_SO4)/modal_mass 
       diameter_4d(:, :, :, 5) =   ((modal_mass/modal_number)**(1./3.))*exp(-3.0*log(sig_fdu)*log(sig_fdu))  
       density_4d(:, :, :, 5) =  (mass(:, :, :, 17)*dens_SO4 + mass(:, :, :, 18)*dens_DU +  mass(:, :, :, 19)*dens_SO4) /modal_mass 
       f_dust_4d(:, :, :, 5) =  mass(:, :, :, 18) /modal_mass
    end where
     diameter_4d(:, :, :, 5)  = min(max(diameter_4d(:, :, :, 5), 0.14e-6), 0.62e-6)    
   
    !!!!!!!! mode 6: fine sea salt, species 20 to 22
    !'SU_A_FSS', 'SS_A_FSS', 'AMM_A_FSS',
    scale_factor  = (((0.56+0.095)/2.0)**3.0)* 1600.0 * 8.0 * 1.0e-18 ! scaling using durinhg training   
    i1 = 20
    i2 = 22
    mass(:, :, :, i1:i2)  =  mass(:, :, :, i1:i2)*scale_factor
    modal_mass =  sum(mass(:, :, :, i1:i2), 4)
    modal_number = num_4d(:, :, :, 6)
    
    sigma_4d(:, :, :, 6) = sig_fss
    
    where ((modal_mass .gt. 0.0 ) .and.  (modal_number .gt. 0.0 )) 
       hygroscopicity_4d(:, :, :, 6) =  (mass(:, :, :, 20)*k_SO4 + mass(:, :, :, 21)*k_SS + mass(:, :, :, 22)*k_SO4)/modal_mass 
       diameter_4d(:, :, :, 6) =   ((modal_mass/modal_number)**(1./3.))*exp(-3.0*log(sig_fss)*log(sig_fss))  
       density_4d(:, :, :, 6) =  (mass(:, :, :, 20)*dens_SO4 + mass(:, :, :, 21)*dens_SS +  mass(:, :, :, 22)*dens_SO4) /modal_mass 
    end where    
        diameter_4d(:, :, :, 6)  = min(max(diameter_4d(:, :, :, 6), 0.095e-6), 0.56e-6) 

    !!!!!!! mode 7: primary organic matter, species 23 to 24
    !'POM_A_PCM', 'BC_A_PCM'
    scale_factor  = (((0.13+0.039)/2.0)**3.0)* 1600.0 * 8.0 * 1.0e-18 ! scaling using durinhg training 
    
    i1 = 23
    i2 = 24
    mass(:, :, :, i1:i2)  =  mass(:, :, :, i1:i2)*scale_factor
    modal_mass =  sum(mass(:, :, :, i1:i2), 4)
    modal_number = num_4d(:, :, :, 7)
    
    sigma_4d(:, :, :, 7) = sig_pcm
    
    where ((modal_mass .gt. 0.0 ) .and.  (modal_number .gt. 0.0 )) 
       hygroscopicity_4d(:, :, :, 7) =  (mass(:, :, :, 23)*k_ORG + mass(:, :, :, 24)*k_BC)/modal_mass 
       diameter_4d(:, :, :, 7) =   ((modal_mass/modal_number)**(1./3.))*exp(-3.0*log(sig_pcm)*log(sig_pcm)) 
       density_4d(:, :, :, 7) =  (mass(:, :, :, 23)*dens_ORG + mass(:, :, :, 24)*dens_BC) /modal_mass
       f_organic_4d(:, :, :, 7) =  mass(:, :, :, 23) /modal_mass
       f_soot_4d(:, :, :, 7) =  mass(:, :, :, 24) /modal_mass
    end where
    diameter_4d(:, :, :, 7)  = min(max(diameter_4d(:, :, :, 7), 0.039e-6), 0.13e-6) 
    
    
    !Trick to pass it as a 3d array
    
    km =  72*n_modes
    num =  reshape(num_4d, (/im, jm, km/))
    hygroscopicity =  reshape(hygroscopicity_4d, (/im, jm, km/))
    diameter =  reshape(diameter_4d, (/im, jm, km/))
    density =  reshape(density_4d, (/im, jm, km/))
    f_dust =  reshape(f_dust_4d, (/im, jm, km/))
    f_organic =  reshape(f_organic_4d, (/im, jm, km/))
    f_soot =  reshape(f_soot_4d, (/im, jm, km/))

   deallocate(q_reshaped)
   deallocate(massnet_in) 
   deallocate(numnet_out)
   deallocate(massnet_out)   
   deallocate (mass)
   deallocate (num_4d)
   deallocate(q)
   deallocate(q_)
   
   deallocate(sigma_4d) 
   deallocate(density_4d) 
   deallocate(hygroscopicity_4d)
   deallocate( diameter_4d)
   deallocate(f_dust_4d)
   deallocate(f_soot_4d)
   deallocate(f_organic_4d)
   deallocate(modal_mass)
   deallocate(modal_number)
      
    
   
   RETURN_(ESMF_SUCCESS)
         
end subroutine MAMnet 
     
end module GOCART_GridCompMod

