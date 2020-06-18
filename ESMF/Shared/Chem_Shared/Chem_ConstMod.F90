!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  Chem_ConstMod --- Defines physical constants
! 
! !INTERFACE:
!

   MODULE  Chem_ConstMod

! !USES:

   Implicit NONE

!
! !DESCRIPTION: This module defines physical constants used throughout the
!               system. Since it is highly desirable that model and analysis
!  share the same constants, it is advised that this module be customized
!  for each model.
!
! !REVISION HISTORY: 
!
!  10oct1999  da Silva  Based on function getcon() from the GEOS GCM. 
!
!EOP
!-------------------------------------------------------------------------

!BOC
                                                                                
!     Computational constants                                                       
!     -----------------------                                                       
      real, parameter ::  vecmax = 65535.5 
      real, parameter ::  caltoj = 4184.   
      real, parameter ::  undef  = 1.e15   ! missing value
                                                                                
!     Astronomical constants                                                        
!     ----------------------                                                        
      real, parameter ::  obliquity       = 23.45
      real, parameter ::  perihelion      = 102.    
      real, parameter ::  eccentricity    = 0.0167   
      real, parameter ::  radius_earth    = 6371e3
      real, parameter ::  vernal_equinox  = 80.5    
      real, parameter ::  summer_solstice = 176.5   
      real, parameter ::  s0              = 1365.0  
                                                                                
!     Terrestrial constants                                                         
!     ---------------------                                                         
!!!   real, parameter ::  grav   = 9.81    
      real, parameter ::  grav   = 9.80616
      real, parameter ::  srfprs = 984.7   
      real, parameter ::  pimean = 984.7   
      real, parameter ::  pstd   = 1000.0  
      real, parameter ::  tstd   = 280.0   
      real, parameter ::  sday   = 86400.0 
      real, parameter ::  ssalb  = 0.99    
      real, parameter ::  co2    = 330.0   
                                                                                
!     Thermodynamic constants                                                       
!     -----------------------                                                       
      real, parameter ::  cpd    = 1004.16 
      real, parameter ::  cpm    = 1004.64
      real, parameter ::  cpv    = 1869.46 
      real, parameter ::  alhl   = 2.499e6 ! Latent heat consensation
      real, parameter ::  alhs   = 2.845e6 ! Latent heat sublimation
      real, parameter ::  stfbol = 5.67e-8 
      real, parameter ::  airmw  = 28.97   ! molecular weight of air
      real, parameter ::  h2omw  = 18.01   ! molecular weight of H2O
      real, parameter ::  runiv  = 8314.3  
!!!   real, parameter ::  rgas   = runiv/airmw
      real, parameter ::  rgas   = 287.04
      real, parameter ::  rvap   = runiv/h2omw
!!!   real, parameter ::  kappa  = rgas/cpd   
      real, parameter ::  kappa  = rgas/cpm   
      real, parameter ::  heatw  = 597.2      
      real, parameter ::  heati  = 680.0      
      real, parameter ::  tice   = 273.16  ! Freezing point    
      real, parameter ::  zvir   = 4.61e2/rgas - 1.
                                      


!     Lapse rate
!     ----------
      real, parameter ::  gamma = 6.5E-3  ! 6.5 Kelvin / Km
                                          
!     turbulence constants                                                          
!     --------------------                                                          
      real, parameter ::  von_karman = 0.4        
                                                                                
!     moisture constants                                                            
!     ------------------                                                            
      real, parameter ::  eps     = 0.622      
      real, parameter ::  virtcon = 0.609      
      real, parameter ::  epsfac  = eps*heatw/rgas*caltoj 
                                                                                
!EOC

     end MODULE Chem_ConstMod
