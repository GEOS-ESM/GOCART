!!!! #include <params.h>

      SUBROUTINE diaglist(diag)

      USE mod_diag
 
! ----------------------------------------------------------------------
!
!  This routine contains the initialization of the diagnostic name,
!  description, and unit.
!
!     type diag_type
!          character*8          :: name         ! name for diagnostic
!          character*16         :: unit         ! unit for diagnostic fields
!          character*80         :: desc         ! description for diagnostic fields
!          integer              :: pick         ! is the diag needed for output 1:yes 0:no
!          logical              :: counted      ! if true, diagnostic is counted
!          integer              :: vdim         ! number of levels 1 or km
!          integer              :: nlist        ! number of output streams
!          integer, pointer     :: fldloc(:,:)  ! location in diagnostic buffer
!          integer, pointer     :: count(:)     ! counter
!          logical, pointer     :: alt(:)       ! if true use alternate name/unit for output
!          character*8          :: aname        ! Alternate name for diagnostic
!          character*16         :: aunit        ! Alternate unit for diagnostic fields
!          character*80         :: adesc        ! description for diagnostic fields
!          real*4               :: convfac      ! conversion factor if alternate unit different than primary
!     endtype
!
! ----------------------------------------------------------------------
!  Variable Declaration 

      type (diag_type)    	:: diag(pdiag)	! diagnostic attributes  (see above)
      integer              	:: n            ! do loop counter

      do n = 1, pdiag
        diag(n)%convfac = 1.
        diag(n)%aname   = '        '
        diag(n)%aunit   = '                '
        diag(n)%adesc   = '                                        '// &
                          '                                        '
        diag(n)%counted = .false.
      enddo


      diag(iALBEDO  )%name = 'ALBEDO  '
      diag(iALBEDO  )%desc = 'Surface Albedo'
      diag(iALBEDO  )%unit = 'fraction'
      diag(iALBEDO  )%aname= 'ALBEDO  '
      diag(iALBEDO  )%adesc= 'Surface Albedo'
      diag(iALBEDO  )%counted = .true.

      diag(iALDIF   )%name = 'ALDIF   '
      diag(iALDIF   )%desc = 'Albedo: longwave, diffuse '
      diag(iALDIF   )%unit = 'fraction'
      diag(iALDIF   )%aname= 'ALBNIRDF'
      diag(iALDIF   )%adesc= 'Diffuse Beam NIR Surface Albedo'

      diag(iALDIR   )%name = 'ALDIR   '
      diag(iALDIR   )%desc = 'Albedo: longwave, direct '
      diag(iALDIR   )%unit = 'fraction'
      diag(iALDIR   )%aname= 'ALBNIRDR'
      diag(iALDIR   )%adesc= 'Direct Beam NIR Surface Albedo'

      diag(iASDIF   )%name = 'ASDIF   '
      diag(iASDIF   )%desc = 'Albedo: shortwave, diffuse '
      diag(iASDIF   )%unit = 'fraction'
      diag(iASDIF   )%aname= 'ALBVISDF'
      diag(iASDIF   )%adesc= 'Diffuse Beam VIS Surface Albedo'

      diag(iASDIR   )%name = 'ASDIR   '
      diag(iASDIR   )%desc = 'Albedo: shortwave, direct '
      diag(iASDIR   )%unit = 'fraction'
      diag(iASDIR   )%aname= 'ALBVISDR'
      diag(iASDIR   )%adesc= 'Direct Beam VIS Surface Albedo'

      diag(iBMA     )%name = 'BMA     '
      diag(iBMA     )%desc = 'Bulk moisture avaliability '
      diag(iBMA     )%unit = 'fraction'

      diag(iBULKTS  )%name = 'BULKTS  '
      diag(iBULKTS  )%desc = 'Bulk surface temperatur (average of tile temperature)'
      diag(iBULKTS  )%unit = 'K'

      diag(iCAPEMX  )%name = 'CAPEMX  '
      diag(iCAPEMX  )%desc = 'Maximum CAPE'
      diag(iCAPEMX  )%unit = 'J/kg'

      diag(iCLDHGH  )%name = 'CLDHGH  '
      diag(iCLDHGH  )%desc = 'Vertically-integrated, random overlap, ' &
                         // 'high cloud amount'
      diag(iCLDHGH  )%unit = 'fraction'
      diag(iCLDHGH  )%aname= 'CLDHI   '
      diag(iCLDHGH  )%adesc= 'High-Level (above 400 hPa) Cloud Fraction' 

      diag(iCLDLOW  )%name = 'CLDLOW  '
      diag(iCLDLOW  )%desc = 'Vertically-integrated, random overlap, ' &
                         // 'low cloud amount'
      diag(iCLDLOW  )%unit = 'fraction'
      diag(iCLDLOW  )%adesc= 'Low-Level (1000-700 hPa) Cloud Fraction' 

      diag(iCLDMED  )%name = 'CLDMED  '
      diag(iCLDMED  )%desc = 'Vertically-integrated, random overlap, ' &
                         // 'mid cloud amount'
      diag(iCLDMED  )%unit = 'fraction'
      diag(iCLDMED  )%aname= 'CLDMID  '
      diag(iCLDMED  )%adesc= 'Mid-Level (700-400 hPa) Cloud Fraction'

      diag(iCLDPRS  )%name = 'CLDPRS  '
      diag(iCLDPRS  )%desc = 'Cloud Top Pressure (when cloudy)'
      diag(iCLDPRS  )%unit = 'Pa'
      diag(iCLDPRS  )%aunit= 'hPa'
      diag(iCLDPRS  )%convfac = 0.01
      diag(iCLDPRS  )%counted = .true.

      diag(iCLDTMP  )%name = 'CLDTMP  '
      diag(iCLDTMP  )%desc = 'Cloud Top Temperature (when cloudy)'
      diag(iCLDTMP  )%unit = 'K'
      diag(iCLDTMP  )%counted = .true.

      diag(iCLDTOT  )%name = 'CLDTOT  '
      diag(iCLDTOT  )%desc = 'Vertically-integrated, random overlap, ' &
                         // 'total cloud cover'
      diag(iCLDTOT  )%unit = 'fraction'
      diag(iCLDTOT  )%aname= 'CLDFRC  '
      diag(iCLDTOT  )%adesc= '2-D Total Cloud Fraction'

      diag(iCNVCLD  )%name = 'CNVCLD  '
      diag(iCNVCLD  )%desc = 'Random overlap total convective cloud amount'
      diag(iCNVCLD  )%unit = 'fraction'

      diag(iEMSFC   )%name = 'EMSFC   '
      diag(iEMSFC   )%desc = 'Bulk surface emissivity'
      diag(iEMSFC   )%unit = 'fraction'

      diag(iFLNS    )%name = 'FLNS    '
      diag(iFLNS    )%desc = 'Net longwave flux at surface'
      diag(iFLNS    )%unit = 'W/m2'
      diag(iFLNS    )%aname= 'RADLWG  '
      diag(iFLNS    )%adesc= 'Net Upward Longwave Flux at the Ground'

      diag(iFLNSC   )%name = 'FLNSC   '
      diag(iFLNSC   )%desc = 'Clear sky net longwave flux at surface'
      diag(iFLNSC   )%unit = 'W/m2'
      diag(iFLNSC   )%aname= 'LWGCLR  '
      diag(iFLNSC   )%adesc= 'Clear Sky Net Longwave Flux at the Ground'

      diag(iFLNT    )%name = 'FLNT    '
      diag(iFLNT    )%desc = 'Net longwave flux at top'
      diag(iFLNT    )%unit = 'W/m2'
      diag(iFLNT    )%aname= 'OLR     '
      diag(iFLNT    )%adesc= 'Outgoing longwave radiation'

      diag(iFLNTC   )%name = 'FLNTC   '
      diag(iFLNTC   )%desc = 'Clear sky net longwave flux at top'
      diag(iFLNTC   )%unit = 'W/m2'
      diag(iFLNTC   )%aname= 'OLRCLR  '
      diag(iFLNTC   )%adesc= 'Clear sky outgoing longwave radiation'

      diag(iFRACLAKE)%name = 'FRACLAKE'
      diag(iFRACLAKE)%desc = 'Lake fraction'
      diag(iFRACLAKE)%unit = 'fraction'

      diag(iFRACVEG )%name = 'FRACVEG '
      diag(iFRACVEG )%desc = 'Vegetation fraction'
      diag(iFRACVEG )%unit = 'fraction'

      diag(iFSDS    )%name = 'FSDS    '
      diag(iFSDS    )%desc = 'Flux shortwave downwelling surface'
      diag(iFSDS    )%unit = 'W/m2'

      diag(iFSNS    )%name = 'FSNS    '
      diag(iFSNS    )%desc = 'Net solar flux at surface'
      diag(iFSNS    )%unit = 'W/m2'
      diag(iFSNS    )%aname= 'RADSWG  '
      diag(iFSNS    )%adesc= 'Net Downward Shortwave Flux at the Ground'

      diag(iFSNSC   )%name = 'FSNSC   '
      diag(iFSNSC   )%desc = 'Clear sky net solar flux at surface '
      diag(iFSNSC   )%unit = 'W/m2'
      diag(iFSNSC   )%aname= 'SWGCLR  '
      diag(iFSNSC   )%adesc= 'Clear Sky Net Downward SW Radiation at the Ground'

      diag(iFSNT    )%name = 'FSNT    '
      diag(iFSNT    )%desc = 'Net solar flux at top'
      diag(iFSNT    )%unit = 'W/m2'

      diag(iFSNTC   )%name = 'FSNTC   '
      diag(iFSNTC   )%desc = 'Clear sky net solar flux at top'
      diag(iFSNTC   )%unit = 'W/m2'

      diag(iGWETROOT)%name = 'GWETROOT'
      diag(iGWETROOT)%desc = 'Root zone soil wetness'
      diag(iGWETROOT)%unit = 'fraction'
      diag(iGWETROOT)%aname= 'GWETROOT'

      diag(iGWETTC )%name = 'GWETTC '
      diag(iGWETTC )%desc = 'Total column soil layer wetness'
      diag(iGWETTC )%unit = 'fraction'
      diag(iGWETTC )%aname= 'GWETTC '

      diag(iGWETTOP )%name = 'GWETTOP '
      diag(iGWETTOP )%desc = 'Top soil layer wetness'
      diag(iGWETTOP )%unit = 'fraction'
      diag(iGWETTOP )%aname= 'GWETTOP '

      diag(iH300    )%name = 'H300    '
      diag(iH300    )%desc = '300 hPa Geopotential height'
      diag(iH300    )%unit = 'm'

      diag(iH500    )%name = 'H500    '
      diag(iH500    )%desc = '500 hPa Geopotential height'
      diag(iH500    )%unit = 'm'

      diag(iHKBETA  )%name = 'HKBETA  '
      diag(iHKBETA  )%desc = 'Overshoot parameter in Hack scheme'
      diag(iHKBETA  )%unit = 'fraction'

      diag(iHKETA   )%name = 'HKETA   '
      diag(iHKETA   )%desc = 'Mass flux without overshoot in Hack scheme'
      diag(iHKETA   )%unit = 'kg/m2 s'

      diag(iHSURF   )%name = 'HSURF   '
      diag(iHSURF   )%desc = 'Surface height'
      diag(iHSURF   )%unit = 'm'
      diag(iHSURF   )%aname= 'PHIS    '
      diag(iHSURF   )%aunit= 'm'
      diag(iHSURF   )%adesc= 'Surface geopotential'

      diag(iHTLCL   )%name = 'HTLCL   '
      diag(iHTLCL   )%desc = ' Height above surface at LCL level'
      diag(iHTLCL   )%unit = 'm'

      diag(iHTMMSE  )%name = 'HTMMSE  '
      diag(iHTMMSE  )%desc = ' Height above surface at maximum moist' &
                                //' static energy level'
      diag(iHTMMSE  )%unit = 'm'

      diag(iLAI     )%name = 'LAI     '
      diag(iLAI     )%desc = 'Leaf area index'
      diag(iLAI     )%unit = 'm2/m2'

      diag(iDTG     )%name = 'DTG     '
      diag(iDTG     )%desc = 'Change in ground temperature'
      diag(iDTG     )%unit = 'K/s'

      diag(iLHFX    )%name = 'LHFX    '
      diag(iLHFX    )%desc = 'Surface latent heat flux'
      diag(iLHFX    )%unit = 'W/m2'
      diag(iLHFX    )%aname= 'EFLUX   '
      diag(iLHFX    )%adesc= 'Latent Heat Flux (pos.upwrd)'

      diag(iLWSH    )%name = 'LWSH    '
      diag(iLWSH    )%desc = 'Liquid water scale height'
      diag(iLWSH    )%unit = 'm'

      diag(iO3DU    )%name = 'O3DU    '
      diag(iO3DU    )%desc = 'Total Column Ozone'
      diag(iO3DU    )%unit = 'Dobson Unit'

      diag(iORO     )%name = 'ORO     '
      diag(iORO     )%desc = 'Surface type flag'
      diag(iORO     )%unit = 'flag'
      diag(iORO     )%aname= 'SURFTYPE'

      diag(iOSR     )%name = 'OSR     '
      diag(iOSR     )%desc = 'Outgoing Shortwave Radiation'
      diag(iOSR     )%unit = 'W/m2'

      diag(iOSRCLR  )%name = 'OSRCLR  '
      diag(iOSRCLR  )%desc = 'Clear Sky Outgoing Shortwave Radiation'
      diag(iOSRCLR  )%unit = 'W/m2'

      diag(iPARDF   )%name = 'PARDF   '
      diag(iPARDF   )%desc = 'Diffuse photosynthetically active radiation' &
                       // ' (0.35-0.70 um)'
      diag(iPARDF   )%unit = 'W/m2'
      diag(iPARDF   )%aname= 'PARDF   '

      diag(iPARDR   )%name = 'PARDR   '
      diag(iPARDR   )%desc = 'Direct photosynthetically active radiation'  &
                       //' (0.35-0.70 um)'
      diag(iPARDR   )%unit = 'W/m2'
      diag(iPARDR   )%aname= 'PARDR   '

      diag(iPBLH    )%name = 'PBLH    '
      diag(iPBLH    )%desc = 'Planetary boundary layer height'
      diag(iPBLH    )%unit = 'm'
      diag(iPBLH    )%aname= 'PBL     '

      diag(iPREACC  )%name = 'PREACC  '
      diag(iPREACC  )%desc = 'Total precipitation rate'
      diag(iPREACC  )%unit = 'mm/day'

      diag(iPRECC   )%name = 'PRECC   '
      diag(iPRECC   )%desc = 'Convective precipitation rate'
      diag(iPRECC   )%unit = 'mm/day'
      diag(iPRECC   )%aname= 'PRECON  '
      diag(iPRECC   )%adesc= 'Convective Precipitation'

      diag(iPRECL   )%name = 'PRECL   '
      diag(iPRECL   )%desc = 'Large-scale precipitation rate'
      diag(iPRECL   )%unit = 'mm/day'

      diag(iPRECL_RH)%name = 'PRECL_RH'
      diag(iPRECL_RH)%desc = 'Precipitation rate due to limit RH'
      diag(iPRECL_RH)%unit = 'mm/day'

      diag(iQ10M    )%name = 'Q10M    '
      diag(iQ10M    )%desc = '10 meter Specific humidity'
      diag(iQ10M    )%unit = 'kg/kg'
      diag(iQ10M    )%adesc= 'Specific Humidity Interpolated to 10 Meters'
      diag(iQ10M    )%aunit= 'g/kg'
      diag(iQ10M    )%convfac = 1000.0

      diag(iQ2M     )%name = 'Q2M     '
      diag(iQ2M     )%desc = '2 meter Specific humidity'
      diag(iQ2M     )%unit = 'kg/kg'
      diag(iQ2M     )%adesc= 'Specific Humidity Interpolated to 2 Meters'
      diag(iQ2M     )%aunit= 'g/kg'
      diag(iQ2M     )%convfac = 1000.0

      diag(iQFLX    )%name = 'QFLX    '
      diag(iQFLX    )%desc = 'Surface water flux'
      diag(iQFLX    )%unit = 'kg/m2/s'
      diag(iQFLX    )%aname= 'EVAP    '
      diag(iQFLX    )%aunit= 'mm/day'
      diag(iQFLX    )%adesc= 'Surface Evaporation'
      diag(iQFLX    )%convfac = 86400.0

      diag(iQPERT   )%name = 'QPERT   '
      diag(iQPERT   )%desc = 'Perturbation specific humidity '  &
                         // '(eddies in PBL)'
      diag(iQPERT   )%unit = 'kg/kg'

      diag(iSHFX    )%name = 'SHFX    '
      diag(iSHFX    )%desc = 'Surface sensible heat flux'
      diag(iSHFX    )%unit = 'W/m2'
      diag(iSHFX    )%aname= 'HFLUX   '
      diag(iSHFX    )%adesc= 'Sensible Heat Flux (pos.upwrd)'

      diag(iSLP     )%name = 'SLP     '
      diag(iSLP     )%desc = 'Sea level pressure'
      diag(iSLP     )%unit = 'Pa'
      diag(iSLP     )%aunit= 'hPa'
      diag(iSLP     )%convfac = 0.01

      diag(iSNOWDP  )%name = 'SNOWDP  '
      diag(iSNOWDP  )%desc = 'Snow depth'
      diag(iSNOWDP  )%unit = 'm'
      diag(iSNOWDP  )%aname= 'SNOWDPTH'
      diag(iSNOWDP  )%adesc= 'Snow Depth (mm)'
      diag(iSNOWDP  )%aunit= 'mm'
      diag(iSNOWDP  )%convfac = 1000.0

      diag(iSNOWH   )%name = 'SNOWH   '
      diag(iSNOWH   )%desc = 'Water equivalent snow depth'
      diag(iSNOWH   )%unit = 'm'
      diag(iSNOWH   )%aname= 'SNOW    '
      diag(iSNOWH   )%adesc= 'Snow Depth (mm water equivalent)'
      diag(iSNOWH   )%aunit= 'mm'
      diag(iSNOWH   )%convfac = 1000.0

      diag(iSOILWC1 )%name = 'SOILWC1 '
      diag(iSOILWC1 )%desc = 'Total column soil water content for soil tiles'
      diag(iSOILWC1 )%unit = 'mm'

      diag(iSOILWC2 )%name = 'SOILWC2 '
      diag(iSOILWC2 )%desc = 'Total column soil water content for water budget'
      diag(iSOILWC2 )%unit = 'mm'

      diag(iSOLIN   )%name = 'SOLIN    '
      diag(iSOLIN   )%desc = 'Solar insolation'
      diag(iSOLIN   )%unit = 'W/m2'
      diag(iSOLIN   )%aname= 'RADSWT  '
      diag(iSOLIN   )%adesc= 'Incident Shortwave Radiation at TOA'

      diag(iSRFRAD  )%name = 'SRFRAD  '
      diag(iSRFRAD  )%desc = 'Net radiative forcing at surface (net SW + downward LW)'
      diag(iSRFRAD  )%unit = 'W/m2'

      diag(iSRUNOFF )%name = 'SRUNOFF '
      diag(iSRUNOFF )%desc = 'Surface runoff'
      diag(iSRUNOFF )%unit = 'mm/s'

      diag(iSURFP   )%name = 'SURFP   '
      diag(iSURFP   )%desc = 'Surface pressure'
      diag(iSURFP   )%unit = 'Pa'
      diag(iSURFP   )%aname= 'PS      '
      diag(iSURFP   )%aunit= 'hPa'
      diag(iSURFP   )%convfac = 0.01  

      diag(iT10M    )%name = 'T10M    '
      diag(iT10M    )%desc = '10 meter temperature'
      diag(iT10M    )%unit = 'K'
      diag(iT10M    )%adesc= 'Temperature Interpolated to 10 Meters'

      diag(iT200    )%name = 'T200    '
      diag(iT200    )%desc = '200 hPa temperature'
      diag(iT200    )%unit = 'K'

      diag(iT2M     )%name = 'T2M     '
      diag(iT2M     )%desc = '2 meter temperature'
      diag(iT2M     )%unit = 'K'
      diag(iT2M     )%adesc= 'Temperature Interpolated to 2 Meters'

      diag(iT850    )%name = 'T850    '
      diag(iT850    )%desc = '850 hPa temperature'
      diag(iT850    )%unit = 'K'
      diag(iT850    )%counted = .true.

      diag(iTAUCLI  )%name = 'TAUCLI  '
      diag(iTAUCLI  )%desc = 'Cloud Optical Depth Ice'
      diag(iTAUCLI  )%unit = 'unitless'
      diag(iTAUCLI  )%counted = .true.

      diag(iTAUCLW  )%name = 'TAUCLW  '
      diag(iTAUCLW  )%desc = 'Cloud Optical Depth Water'
      diag(iTAUCLW  )%unit = 'unitless'
      diag(iTAUCLW  )%counted = .true.

      diag(iTAUGWX  )%name = 'TAUGWX  '
      diag(iTAUGWX  )%desc = 'East-west gravity wave drag surface stress'
      diag(iTAUGWX  )%unit = 'N/m2'
      diag(iTAUGWX  )%aname= 'GWDUS   '
      diag(iTAUGWX  )%adesc= 'Zonal Wind Gravity Wave Surface Stress'

      diag(iTAUGWY  )%name = 'TAUGWY  '
      diag(iTAUGWY  )%desc = 'North-south gravity wave drag surface stress'
      diag(iTAUGWY  )%unit = 'N/m2'
      diag(iTAUGWY  )%aname= 'GWDVS   '
      diag(iTAUGWY  )%adesc= 'Meridional Wind Gravity Wave Surface Stress'

      diag(iTAUX    )%name = 'TAUX    '
      diag(iTAUX    )%desc = 'X-component (east-west) of surface stress'
      diag(iTAUX    )%unit = 'N/m2'
      diag(iTAUX    )%aname= 'UFLUX   '
      diag(iTAUX    )%adesc= 'Zonal Wind Surface Stress'

      diag(iTAUY    )%name = 'TAUY    '
      diag(iTAUY    )%desc = 'Y-component (north-south) of surface stress'
      diag(iTAUY    )%unit = 'N/m2'
      diag(iTAUY    )%aname= 'VFLUX   '
      diag(iTAUY    )%adesc= 'Meridional Wind Surface Stress'

      diag(iTHICK   )%name = 'THICK    '
      diag(iTHICK   )%desc = 'Thickness of the 500 hPa to 1000 hPa Layer (5400m = freezing line)'
      diag(iTHICK   )%unit = 'm'

      diag(iTLAKE1  )%name = 'TLAKE1  '
      diag(iTLAKE1  )%desc = 'Top layer lake temperature'
      diag(iTLAKE1  )%unit = 'K'

      diag(iTPERT   )%name = 'TPERT   '
      diag(iTPERT   )%desc = 'Perturbation temperature (eddies in PBL)'
      diag(iTPERT   )%unit = 'K'

      diag(iTQ      )%name = 'TQ      '
      diag(iTQ      )%desc = 'Total precipitable water'
      diag(iTQ      )%unit = 'kg/m2'
      diag(iTQ      )%aname= 'TPW     '
      diag(iTQ      )%aunit= 'g/cm2'
      diag(iTQ      )%convfac = 0.1

      diag(iTRAD    )%name = 'TRAD    '
      diag(iTRAD    )%desc = 'Surface brightness temperature (average of tile brightness temperature'
      diag(iTRAD    )%unit = 'K'

      diag(iTROPP   )%name = 'TROPP   '
      diag(iTROPP   )%desc = 'Tropopause pressure'
      diag(iTROPP   )%unit = 'Pa'
      diag(iTROPP   )%aunit= 'hPa'
      diag(iTROPP   )%convfac = 0.01
      diag(iTROPP   )%counted = .true.

      diag(iTROPQ   )%name = 'TROPQ   '
      diag(iTROPQ   )%desc = 'Tropopause specific humidity'
      diag(iTROPQ   )%unit = 'kg/kg'
      diag(iTROPQ   )%aunit= 'g/kg'
      diag(iTROPQ   )%convfac = 1000.
      diag(iTROPQ   )%counted = .true.

      diag(iTROPT   )%name = 'TROPT   '
      diag(iTROPT   )%desc = 'Tropopause temperature'
      diag(iTROPT   )%unit = 'K'
      diag(iTROPT   )%counted = .true.

      diag(iTSKIN   )%name = 'TSKIN   '
      diag(iTSKIN   )%desc = 'Surface skin temperature for CERES (derived from tile IR flux'
      diag(iTSKIN   )%unit = 'K'

      diag(iTSLAKE  )%name = 'TSLAKE  '
      diag(iTSLAKE  )%desc = 'Lake skin temperature'
      diag(iTSLAKE  )%unit = 'K'

      diag(iTSOIL1  )%name = 'TSOIL1  '
      diag(iTSOIL1  )%desc = 'Top layer soil temperature'
      diag(iTSOIL1  )%unit = 'K'

      diag(iTVEG    )%name = 'TVEG    '
      diag(iTVEG    )%desc = 'Bulk vegetation temperature'
      diag(iTVEG    )%unit = 'K'

      diag(iU10M    )%name = 'U10M    '
      diag(iU10M    )%desc = '10 meter U wind'
      diag(iU10M    )%unit = 'm/s'
      diag(iU10M    )%adesc= 'Zonal Wind Interpolated to 10 Meters'

      diag(iU200    )%name = 'U200    '
      diag(iU200    )%desc = '200 hPa U wind'
      diag(iU200    )%unit = 'm/s'

      diag(iU2M     )%name = 'U2M     '
      diag(iU2M     )%desc = '2 meter U wind'
      diag(iU2M     )%unit = 'm/s'
      diag(iU2M     )%adesc= 'Zonal Wind Interpolated to 2 Meters'

      diag(iU500    )%name = 'U500    '
      diag(iU500    )%desc = '500 hPa U wind'
      diag(iU500    )%unit = 'm/s'

      diag(iU850    )%name = 'U850    '
      diag(iU850    )%desc = '850 hPa U wind'
      diag(iU850    )%unit = 'm/s'
      diag(iU850    )%counted = .true.

      diag(iUSTAR   )%name = 'USTAR   '
      diag(iUSTAR   )%desc = 'Surface friction velocity'
      diag(iUSTAR   )%unit = 'm/s'
      diag(iUSTAR   )%adesc= 'Friction velocity'

      diag(iV10M    )%name = 'V10M    '
      diag(iV10M    )%desc = '10 meter V wind'
      diag(iV10M    )%unit = 'm/s'
      diag(iV10M    )%adesc= 'Meridional Wind Interpolated to 10 Meters'

      diag(iV200    )%name = 'V200    '
      diag(iV200    )%desc = '200 hPa V wind'
      diag(iV200    )%unit = 'm/s'

      diag(iV2M     )%name = 'V2M     '
      diag(iV2M     )%desc = '2 meter V wind'
      diag(iV2M     )%unit = 'm/s'
      diag(iV2M     )%adesc= 'Meridional Wind Interpolated to 2 Meters'

      diag(iV500    )%name = 'V500    '
      diag(iV500    )%desc = '500 hPa V wind'
      diag(iV500    )%unit = 'm/s'
 
      diag(iV850    )%name = 'V850    '
      diag(iV850    )%desc = '850 hPa V wind'
      diag(iV850    )%unit = 'm/s'
      diag(iV850    )%counted = .true.

      diag(iVAVET   )%name = 'VAVET   '
      diag(iVAVET   )%desc = 'Vertically averaged temperature'
      diag(iVAVET   )%unit = 'K'

      diag(iVAVEU   )%name = 'VAVEU   '
      diag(iVAVEU   )%desc = 'Vertically averaged U Wind'
      diag(iVAVEU   )%unit = 'm/s'

      diag(iVAVEUQ  )%name = 'VAVEUQ  '
      diag(iVAVEUQ  )%desc = 'Vertically averaged U Wind * specific humidity'
      diag(iVAVEUQ  )%unit = 'm/s*kg/kg'
      diag(iVAVEUQ  )%aunit= 'm/s*g/kg'
      diag(iVAVEUQ  )%convfac = 1000.0 

      diag(iVAVEUT  )%name = 'VAVEUT  '
      diag(iVAVEUT  )%desc = 'Vertically averaged U Wind * temperature'
      diag(iVAVEUT  )%unit = 'Km/s'

      diag(iVAVEV   )%name = 'VAVEV   '
      diag(iVAVEV   )%desc = 'Vertically averaged V Wind'
      diag(iVAVEV   )%unit = 'm/s'

      diag(iVAVEVQ  )%name = 'VAVEVQ  '
      diag(iVAVEVQ  )%desc = 'Vertically averaged V Wind * specific humidity'
      diag(iVAVEVQ  )%unit = 'm/s*kg/kg'
      diag(iVAVEVQ  )%aunit= 'm/s*g/kg'
      diag(iVAVEVQ  )%convfac = 1000.0 

      diag(iVAVEVT  )%name = 'VAVEVT  '
      diag(iVAVEVT  )%desc = 'Vertically averaged V Wind * temperature'
      diag(iVAVEVT  )%unit = 'Km/s'

      diag(iWATCAN  )%name = 'WATCAN  '
      diag(iWATCAN  )%desc = 'Water on canopy'
      diag(iWATCAN  )%unit = 'mm'

      diag(iWBALLAKE)%name = 'WBALLAKE'
      diag(iWBALLAKE)%desc = 'Water balance term for lake, wet land, and glacier'
      diag(iWBALLAKE)%unit = 'mm/s'

      diag(iZ0H     )%name = 'Z0H     '
      diag(iZ0H     )%desc = 'Roughness length, sensible heat'
      diag(iZ0H     )%unit = 'm'

      diag(iZ0M     )%name = 'Z0M     '
      diag(iZ0M     )%desc = 'Roughness length, momentum'
      diag(iZ0M     )%unit = 'm'

      diag(iZMMB    )%name = 'ZMMB    '
      diag(iZMMB    )%desc = 'Cloud base mass flux from Z&M scheme'
      diag(iZMMB    )%unit = 'kg/m2/s'

      diag(iZMPR    )%name = 'ZMPR    '
      diag(iZMPR    )%desc = 'Precipitation from Z&M scheme'
      diag(iZMPR    )%unit = 'mm/day'

      diag(iZPD     )%name = 'ZPD     '
      diag(iZPD     )%desc = 'Displacement height'
      diag(iZPD     )%unit = 'm'

#ifdef FVCHEM

      diag(iDUEM001 )%name = 'DUEM001 '
      diag(iDUEM001 )%desc = 'Dust Emission Bin 1'
      diag(iDUEM001 )%unit = 'kg/m2/s'

      diag(iDUEM002 )%name = 'DUEM002 '
      diag(iDUEM002 )%desc = 'Dust Emission Bin 2'
      diag(iDUEM002 )%unit = 'kg/m2/s'

      diag(iDUEM003 )%name = 'DUEM003 '
      diag(iDUEM003 )%desc = 'Dust Emission Bin 3'
      diag(iDUEM003 )%unit = 'kg/m2/s'

      diag(iDUEM004 )%name = 'DUEM004 '
      diag(iDUEM004 )%desc = 'Dust Emission Bin 4'
      diag(iDUEM004 )%unit = 'kg/m2/s'

      diag(iDUEM005 )%name = 'DUEM005 '
      diag(iDUEM005 )%desc = 'Dust Emission Bin 5'
      diag(iDUEM005 )%unit = 'kg/m2/s'

      diag(iDUEM006 )%name = 'DUEM006 '
      diag(iDUEM006 )%desc = 'Dust Emission Bin 6'
      diag(iDUEM006 )%unit = 'kg/m2/s'

      diag(iDUEM007 )%name = 'DUEM007 '
      diag(iDUEM007 )%desc = 'Dust Emission Bin 7'
      diag(iDUEM007 )%unit = 'kg/m2/s'

      diag(iDUEM008 )%name = 'DUEM008 '
      diag(iDUEM008 )%desc = 'Dust Emission Bin 8'
      diag(iDUEM008 )%unit = 'kg/m2/s'

      diag(iDUSD001 )%name = 'DUSD001 '
      diag(iDUSD001 )%desc = 'Dust Sedimentation Bin 1'
      diag(iDUSD001 )%unit = 'kg/m2/s'

      diag(iDUSD002 )%name = 'DUSD002 '
      diag(iDUSD002 )%desc = 'Dust Sedimentation Bin 2'
      diag(iDUSD002 )%unit = 'kg/m2/s'

      diag(iDUSD003 )%name = 'DUSD003 '
      diag(iDUSD003 )%desc = 'Dust Sedimentation Bin 3'
      diag(iDUSD003 )%unit = 'kg/m2/s'

      diag(iDUSD004 )%name = 'DUSD004 '
      diag(iDUSD004 )%desc = 'Dust Sedimentation Bin 4'
      diag(iDUSD004 )%unit = 'kg/m2/s'

      diag(iDUSD005 )%name = 'DUSD005 '
      diag(iDUSD005 )%desc = 'Dust Sedimentation Bin 5'
      diag(iDUSD005 )%unit = 'kg/m2/s'

      diag(iDUSD006 )%name = 'DUSD006 '
      diag(iDUSD006 )%desc = 'Dust Sedimentation Bin 6'
      diag(iDUSD006 )%unit = 'kg/m2/s'

      diag(iDUSD007 )%name = 'DUSD007 '
      diag(iDUSD007 )%desc = 'Dust Sedimentation Bin 7'
      diag(iDUSD007 )%unit = 'kg/m2/s'

      diag(iDUSD008 )%name = 'DUSD008 '
      diag(iDUSD008 )%desc = 'Dust Sedimentation Bin 8'
      diag(iDUSD008 )%unit = 'kg/m2/s'

      diag(iDUDP001 )%name = 'DUDP001 '
      diag(iDUDP001 )%desc = 'Dust Dry Deposition Bin 1'
      diag(iDUDP001 )%unit = 'kg/m2/s'

      diag(iDUDP002 )%name = 'DUDP002 '
      diag(iDUDP002 )%desc = 'Dust Dry Deposition Bin 2'
      diag(iDUDP002 )%unit = 'kg/m2/s'

      diag(iDUDP003 )%name = 'DUDP003 '
      diag(iDUDP003 )%desc = 'Dust Dry Deposition Bin 3'
      diag(iDUDP003 )%unit = 'kg/m2/s'

      diag(iDUDP004 )%name = 'DUDP004 '
      diag(iDUDP004 )%desc = 'Dust Dry Deposition Bin 4'
      diag(iDUDP004 )%unit = 'kg/m2/s'

      diag(iDUDP005 )%name = 'DUDP005 '
      diag(iDUDP005 )%desc = 'Dust Dry Deposition Bin 5'
      diag(iDUDP005 )%unit = 'kg/m2/s'

      diag(iDUDP006 )%name = 'DUDP006 '
      diag(iDUDP006 )%desc = 'Dust Dry Deposition Bin 6'
      diag(iDUDP006 )%unit = 'kg/m2/s'

      diag(iDUDP007 )%name = 'DUDP007 '
      diag(iDUDP007 )%desc = 'Dust Dry Deposition Bin 7'
      diag(iDUDP007 )%unit = 'kg/m2/s'

      diag(iDUDP008 )%name = 'DUDP008 '
      diag(iDUDP008 )%desc = 'Dust Dry Deposition Bin 8'
      diag(iDUDP008 )%unit = 'kg/m2/s'

      diag(iDUWT001 )%name = 'DUWT001 '
      diag(iDUWT001 )%desc = 'Dust Wet Deposition Bin 1'
      diag(iDUWT001 )%unit = 'kg/m2/s'

      diag(iDUWT002 )%name = 'DUWT002 '
      diag(iDUWT002 )%desc = 'Dust Wet Deposition Bin 2'
      diag(iDUWT002 )%unit = 'kg/m2/s'

      diag(iDUWT003 )%name = 'DUWT003 '
      diag(iDUWT003 )%desc = 'Dust Wet Deposition Bin 3'
      diag(iDUWT003 )%unit = 'kg/m2/s'

      diag(iDUWT004 )%name = 'DUWT004 '
      diag(iDUWT004 )%desc = 'Dust Wet Deposition Bin 4'
      diag(iDUWT004 )%unit = 'kg/m2/s'

      diag(iDUWT005 )%name = 'DUWT005 '
      diag(iDUWT005 )%desc = 'Dust Wet Deposition Bin 5'
      diag(iDUWT005 )%unit = 'kg/m2/s'

      diag(iDUWT006 )%name = 'DUWT006 '
      diag(iDUWT006 )%desc = 'Dust Wet Deposition Bin 6'
      diag(iDUWT006 )%unit = 'kg/m2/s'

      diag(iDUWT007 )%name = 'DUWT007 '
      diag(iDUWT007 )%desc = 'Dust Wet Deposition Bin 7'
      diag(iDUWT007 )%unit = 'kg/m2/s'

      diag(iDUWT008 )%name = 'DUWT008 '
      diag(iDUWT008 )%desc = 'Dust Wet Deposition Bin 8'
      diag(iDUWT008 )%unit = 'kg/m2/s'

      diag(iDUSV001 )%name = 'DUSV001 '
      diag(iDUSV001 )%desc = 'Dust Convective Scavenging Bin 1'
      diag(iDUSV001 )%unit = 'kg/m2/s'

      diag(iDUSV002 )%name = 'DUSV002 '
      diag(iDUSV002 )%desc = 'Dust Convective Scavenging Bin 2'
      diag(iDUSV002 )%unit = 'kg/m2/s'

      diag(iDUSV003 )%name = 'DUSV003 '
      diag(iDUSV003 )%desc = 'Dust Convective Scavenging Bin 3'
      diag(iDUSV003 )%unit = 'kg/m2/s'

      diag(iDUSV004 )%name = 'DUSV004 '
      diag(iDUSV004 )%desc = 'Dust Convective Scavenging Bin 4'
      diag(iDUSV004 )%unit = 'kg/m2/s'

      diag(iDUSV005 )%name = 'DUSV005 '
      diag(iDUSV005 )%desc = 'Dust Convective Scavenging Bin 5'
      diag(iDUSV005 )%unit = 'kg/m2/s'

      diag(iDUSV006 )%name = 'DUSV006 '
      diag(iDUSV006 )%desc = 'Dust Convective Scavenging Bin 6'
      diag(iDUSV006 )%unit = 'kg/m2/s'

      diag(iDUSV007 )%name = 'DUSV007 '
      diag(iDUSV007 )%desc = 'Dust Convective Scavenging Bin 7'
      diag(iDUSV007 )%unit = 'kg/m2/s'

      diag(iDUSV008 )%name = 'DUSV008 '
      diag(iDUSV008 )%desc = 'Dust Convective Scavenging Bin 8'
      diag(iDUSV008 )%unit = 'kg/m2/s'

      diag(iDUSMASS )%name = 'DUSMASS '
      diag(iDUSMASS )%desc = 'Dust Surface Mass Concentration'
      diag(iDUSMASS )%unit = 'kg/m3'

      diag(iDUCMASS )%name = 'DUCMASS '
      diag(iDUCMASS )%desc = 'Dust Column Mass Density'
      diag(iDUCMASS )%unit = 'kg/m2'

      diag(iDUSMASS1 )%name = 'DUSMASS1 '
      diag(iDUSMASS1 )%desc = 'Dust sub-micron Surface Mass Concentration'
      diag(iDUSMASS1 )%unit = 'kg/m3'

      diag(iDUCMASS1 )%name = 'DUCMASS1 '
      diag(iDUCMASS1 )%desc = 'Dust sub-micron Column Mass Density'
      diag(iDUCMASS1 )%unit = 'kg/m2'

      diag(iDUEXTTAU )%name = 'DUEXTTAU '
      diag(iDUEXTTAU )%desc = 'Dust Extinction AOT [550 nm]'
      diag(iDUEXTTAU )%unit = 'unitless'

      diag(iDUSCATAU )%name = 'DUSCATAU '
      diag(iDUSCATAU )%desc = 'Dust Scattering AOT [550 nm]'
      diag(iDUSCATAU )%unit = 'unitless'

      diag(iDUAERIDX )%name = 'DUAERIDX '
      diag(iDUAERIDX )%desc = 'Dust TOMS UV Aerosol Index'
      diag(iDUAERIDX )%unit = 'unitless'

      diag(iDUSM25 )%name = 'DUSM25 '
      diag(iDUSM25 )%desc = 'Dust Surface Mass Concentration (PM2.5)'
      diag(iDUSM25 )%unit = 'kg/m3'

      diag(iDUCM25 )%name = 'DUCM25 '
      diag(iDUCM25 )%desc = 'Dust Column Mass Density (PM2.5)'
      diag(iDUCM25 )%unit = 'kg/m2'

      diag(iDUEXTT25 )%name = 'DUEXTT25 '
      diag(iDUEXTT25 )%desc = 'Dust Extinction AOT [550 nm] (PM2.5)'
      diag(iDUEXTT25 )%unit = 'unitless'

      diag(iDUSCAT25 )%name = 'DUSCAT25 '
      diag(iDUSCAT25 )%desc = 'Dust Scattering AOT [550 nm] (PM2.5)'
      diag(iDUSCAT25 )%unit = 'unitless'

      diag(iSSEM001 )%name = 'SSEM001 '
      diag(iSSEM001 )%desc = 'Seasalt Emission Bin 1'
      diag(iSSEM001 )%unit = 'kg/m2/s'

      diag(iSSEM002 )%name = 'SSEM002 '
      diag(iSSEM002 )%desc = 'Seasalt Emission Bin 2'
      diag(iSSEM002 )%unit = 'kg/m2/s'

      diag(iSSEM003 )%name = 'SSEM003 '
      diag(iSSEM003 )%desc = 'Seasalt Emission Bin 3'
      diag(iSSEM003 )%unit = 'kg/m2/s'

      diag(iSSEM004 )%name = 'SSEM004 '
      diag(iSSEM004 )%desc = 'Seasalt Emission Bin 4'
      diag(iSSEM004 )%unit = 'kg/m2/s'

      diag(iSSEM005 )%name = 'SSEM005 '
      diag(iSSEM005 )%desc = 'Seasalt Emission Bin 5'
      diag(iSSEM005 )%unit = 'kg/m2/s'

      diag(iSSEM006 )%name = 'SSEM006 '
      diag(iSSEM006 )%desc = 'Seasalt Emission Bin 6'
      diag(iSSEM006 )%unit = 'kg/m2/s'

      diag(iSSEM007 )%name = 'SSEM007 '
      diag(iSSEM007 )%desc = 'Seasalt Emission Bin 7'
      diag(iSSEM007 )%unit = 'kg/m2/s'

      diag(iSSEM008 )%name = 'SSEM008 '
      diag(iSSEM008 )%desc = 'Seasalt Emission Bin 8'
      diag(iSSEM008 )%unit = 'kg/m2/s'

      diag(iSSSD001 )%name = 'SSSD001 '
      diag(iSSSD001 )%desc = 'Seasalt Sedimentation Bin 1'
      diag(iSSSD001 )%unit = 'kg/m2/s'

      diag(iSSSD002 )%name = 'SSSD002 '
      diag(iSSSD002 )%desc = 'Seasalt Sedimentation Bin 2'
      diag(iSSSD002 )%unit = 'kg/m2/s'

      diag(iSSSD003 )%name = 'SSSD003 '
      diag(iSSSD003 )%desc = 'Seasalt Sedimentation Bin 3'
      diag(iSSSD003 )%unit = 'kg/m2/s'

      diag(iSSSD004 )%name = 'SSSD004 '
      diag(iSSSD004 )%desc = 'Seasalt Sedimentation Bin 4'
      diag(iSSSD004 )%unit = 'kg/m2/s'

      diag(iSSSD005 )%name = 'SSSD005 '
      diag(iSSSD005 )%desc = 'Seasalt Sedimentation Bin 5'
      diag(iSSSD005 )%unit = 'kg/m2/s'

      diag(iSSSD006 )%name = 'SSSD006 '
      diag(iSSSD006 )%desc = 'Seasalt Sedimentation Bin 6'
      diag(iSSSD006 )%unit = 'kg/m2/s'

      diag(iSSSD007 )%name = 'SSSD007 '
      diag(iSSSD007 )%desc = 'Seasalt Sedimentation Bin 7'
      diag(iSSSD007 )%unit = 'kg/m2/s'

      diag(iSSSD008 )%name = 'SSSD008 '
      diag(iSSSD008 )%desc = 'Seasalt Sedimentation Bin 8'
      diag(iSSSD008 )%unit = 'kg/m2/s'

      diag(iSSDP001 )%name = 'SSDP001 '
      diag(iSSDP001 )%desc = 'Seasalt Dry Deposition Bin 1'
      diag(iSSDP001 )%unit = 'kg/m2/s'

      diag(iSSDP002 )%name = 'SSDP002 '
      diag(iSSDP002 )%desc = 'Seasalt Dry Deposition Bin 2'
      diag(iSSDP002 )%unit = 'kg/m2/s'

      diag(iSSDP003 )%name = 'SSDP003 '
      diag(iSSDP003 )%desc = 'Seasalt Dry Deposition Bin 3'
      diag(iSSDP003 )%unit = 'kg/m2/s'

      diag(iSSDP004 )%name = 'SSDP004 '
      diag(iSSDP004 )%desc = 'Seasalt Dry Deposition Bin 4'
      diag(iSSDP004 )%unit = 'kg/m2/s'

      diag(iSSDP005 )%name = 'SSDP005 '
      diag(iSSDP005 )%desc = 'Seasalt Dry Deposition Bin 5'
      diag(iSSDP005 )%unit = 'kg/m2/s'

      diag(iSSDP006 )%name = 'SSDP006 '
      diag(iSSDP006 )%desc = 'Seasalt Dry Deposition Bin 6'
      diag(iSSDP006 )%unit = 'kg/m2/s'

      diag(iSSDP007 )%name = 'SSDP007 '
      diag(iSSDP007 )%desc = 'Seasalt Dry Deposition Bin 7'
      diag(iSSDP007 )%unit = 'kg/m2/s'

      diag(iSSDP008 )%name = 'SSDP008 '
      diag(iSSDP008 )%desc = 'Seasalt Dry Deposition Bin 8'
      diag(iSSDP008 )%unit = 'kg/m2/s'

      diag(iSSWT001 )%name = 'SSWT001 '
      diag(iSSWT001 )%desc = 'Seasalt Wet Deposition Bin 1'
      diag(iSSWT001 )%unit = 'kg/m2/s'

      diag(iSSWT002 )%name = 'SSWT002 '
      diag(iSSWT002 )%desc = 'Seasalt Wet Deposition Bin 2'
      diag(iSSWT002 )%unit = 'kg/m2/s'

      diag(iSSWT003 )%name = 'SSWT003 '
      diag(iSSWT003 )%desc = 'Seasalt Wet Deposition Bin 3'
      diag(iSSWT003 )%unit = 'kg/m2/s'

      diag(iSSWT004 )%name = 'SSWT004 '
      diag(iSSWT004 )%desc = 'Seasalt Wet Deposition Bin 4'
      diag(iSSWT004 )%unit = 'kg/m2/s'

      diag(iSSWT005 )%name = 'SSWT005 '
      diag(iSSWT005 )%desc = 'Seasalt Wet Deposition Bin 5'
      diag(iSSWT005 )%unit = 'kg/m2/s'

      diag(iSSWT006 )%name = 'SSWT006 '
      diag(iSSWT006 )%desc = 'Seasalt Wet Deposition Bin 6'
      diag(iSSWT006 )%unit = 'kg/m2/s'

      diag(iSSWT007 )%name = 'SSWT007 '
      diag(iSSWT007 )%desc = 'Seasalt Wet Deposition Bin 7'
      diag(iSSWT007 )%unit = 'kg/m2/s'

      diag(iSSWT008 )%name = 'SSWT008 '
      diag(iSSWT008 )%desc = 'Seasalt Wet Deposition Bin 8'
      diag(iSSWT008 )%unit = 'kg/m2/s'

      diag(iSSSV001 )%name = 'SSSV001 '
      diag(iSSSV001 )%desc = 'Seasalt Convective Scavenging Bin 1'
      diag(iSSSV001 )%unit = 'kg/m2/s'

      diag(iSSSV002 )%name = 'SSSV002 '
      diag(iSSSV002 )%desc = 'Seasalt Convective Scavenging Bin 2'
      diag(iSSSV002 )%unit = 'kg/m2/s'

      diag(iSSSV003 )%name = 'SSSV003 '
      diag(iSSSV003 )%desc = 'Seasalt Convective Scavenging Bin 3'
      diag(iSSSV003 )%unit = 'kg/m2/s'

      diag(iSSSV004 )%name = 'SSSV004 '
      diag(iSSSV004 )%desc = 'Seasalt Convective Scavenging Bin 4'
      diag(iSSSV004 )%unit = 'kg/m2/s'

      diag(iSSSV005 )%name = 'SSSV005 '
      diag(iSSSV005 )%desc = 'Seasalt Convective Scavenging Bin 5'
      diag(iSSSV005 )%unit = 'kg/m2/s'

      diag(iSSSV006 )%name = 'SSSV006 '
      diag(iSSSV006 )%desc = 'Seasalt Convective Scavenging Bin 6'
      diag(iSSSV006 )%unit = 'kg/m2/s'

      diag(iSSSV007 )%name = 'SSSV007 '
      diag(iSSSV007 )%desc = 'Seasalt Convective Scavenging Bin 7'
      diag(iSSSV007 )%unit = 'kg/m2/s'

      diag(iSSSV008 )%name = 'SSSV008 '
      diag(iSSSV008 )%desc = 'Seasalt Convective Scavenging Bin 8'
      diag(iSSSV008 )%unit = 'kg/m2/s'

      diag(iSSSMASS )%name = 'SSSMASS '
      diag(iSSSMASS )%desc = 'Seasalt Surface Mass Concentration'
      diag(iSSSMASS )%unit = 'kg/m3'

      diag(iSSCMASS )%name = 'SSCMASS '
      diag(iSSCMASS )%desc = 'Seasalt Column Mass Density'
      diag(iSSCMASS )%unit = 'kg/m2'

      diag(iSSEXTTAU )%name = 'SSEXTTAU '
      diag(iSSEXTTAU )%desc = 'Seasalt Extinction AOT [550 nm]'
      diag(iSSEXTTAU )%unit = 'unitless'

      diag(iSSSCATAU )%name = 'SSSCATAU '
      diag(iSSSCATAU )%desc = 'Seasalt Scattering AOT [550 nm]'
      diag(iSSSCATAU )%unit = 'unitless'

      diag(iSSSM25 )%name = 'SSSM25 '
      diag(iSSSM25 )%desc = 'Seasalt Surface Mass Concentration (PM2.5)'
      diag(iSSSM25 )%unit = 'kg/m3'

      diag(iSSCM25 )%name = 'SSCM25 '
      diag(iSSCM25 )%desc = 'Seasalt Column Mass Density (PM2.5)'
      diag(iSSCM25 )%unit = 'kg/m2'

      diag(iSSEXTT25 )%name = 'SSEXTT25 '
      diag(iSSEXTT25 )%desc = 'Seasalt Extinction AOT [550 nm] (PM2.5)'
      diag(iSSEXTT25 )%unit = 'unitless'

      diag(iSSSCAT25 )%name = 'SSSCAT25 '
      diag(iSSSCAT25 )%desc = 'Seasalt Scattering AOT [550 nm] (PM2.5)'
      diag(iSSSCAT25 )%unit = 'unitless'

      diag(iBCEM001 )%name = 'BCEM001 '
      diag(iBCEM001 )%desc = 'Black Carbon Emission Bin 1'
      diag(iBCEM001 )%unit = 'kg/m2/s'

      diag(iBCEM002 )%name = 'BCEM002 '
      diag(iBCEM002 )%desc = 'Black Carbon Emission Bin 2'
      diag(iBCEM002 )%unit = 'kg/m2/s'

      diag(iBCEM003 )%name = 'BCEM003 '
      diag(iBCEM003 )%desc = 'Black Carbon Emission Bin 3'
      diag(iBCEM003 )%unit = 'kg/m2/s'

      diag(iBCEM004 )%name = 'BCEM004 '
      diag(iBCEM004 )%desc = 'Black Carbon Emission Bin 4'
      diag(iBCEM004 )%unit = 'kg/m2/s'

      diag(iBCEM005 )%name = 'BCEM005 '
      diag(iBCEM005 )%desc = 'Black Carbon Emission Bin 5'
      diag(iBCEM005 )%unit = 'kg/m2/s'

      diag(iBCEM006 )%name = 'BCEM006 '
      diag(iBCEM006 )%desc = 'Black Carbon Emission Bin 6'
      diag(iBCEM006 )%unit = 'kg/m2/s'

      diag(iBCEM007 )%name = 'BCEM007 '
      diag(iBCEM007 )%desc = 'Black Carbon Emission Bin 7'
      diag(iBCEM007 )%unit = 'kg/m2/s'

      diag(iBCEM008 )%name = 'BCEM008 '
      diag(iBCEM008 )%desc = 'Black Carbon Emission Bin 8'
      diag(iBCEM008 )%unit = 'kg/m2/s'

      diag(iBCDP001 )%name = 'BCDP001 '
      diag(iBCDP001 )%desc = 'Black Carbon Deposition Bin 1'
      diag(iBCDP001 )%unit = 'kg/m2/s'

      diag(iBCDP002 )%name = 'BCDP002 '
      diag(iBCDP002 )%desc = 'Black Carbon Deposition Bin 2'
      diag(iBCDP002 )%unit = 'kg/m2/s'

      diag(iBCDP003 )%name = 'BCDP003 '
      diag(iBCDP003 )%desc = 'Black Carbon Deposition Bin 3'
      diag(iBCDP003 )%unit = 'kg/m2/s'

      diag(iBCDP004 )%name = 'BCDP004 '
      diag(iBCDP004 )%desc = 'Black Carbon Deposition Bin 4'
      diag(iBCDP004 )%unit = 'kg/m2/s'

      diag(iBCDP005 )%name = 'BCDP005 '
      diag(iBCDP005 )%desc = 'Black Carbon Deposition Bin 5'
      diag(iBCDP005 )%unit = 'kg/m2/s'

      diag(iBCDP006 )%name = 'BCDP006 '
      diag(iBCDP006 )%desc = 'Black Carbon Deposition Bin 6'
      diag(iBCDP006 )%unit = 'kg/m2/s'

      diag(iBCDP007 )%name = 'BCDP007 '
      diag(iBCDP007 )%desc = 'Black Carbon Deposition Bin 7'
      diag(iBCDP007 )%unit = 'kg/m2/s'

      diag(iBCDP008 )%name = 'BCDP008 '
      diag(iBCDP008 )%desc = 'Black Carbon Deposition Bin 8'
      diag(iBCDP008 )%unit = 'kg/m2/s'

      diag(iBCWT001 )%name = 'BCWT001 '
      diag(iBCWT001 )%desc = 'Black Carbon Wet Deposition Bin 1'
      diag(iBCWT001 )%unit = 'kg/m2/s'

      diag(iBCWT002 )%name = 'BCWT002 '
      diag(iBCWT002 )%desc = 'Black Carbon Wet Deposition Bin 2'
      diag(iBCWT002 )%unit = 'kg/m2/s'

      diag(iBCWT003 )%name = 'BCWT003 '
      diag(iBCWT003 )%desc = 'Black Carbon Wet Deposition Bin 3'
      diag(iBCWT003 )%unit = 'kg/m2/s'

      diag(iBCWT004 )%name = 'BCWT004 '
      diag(iBCWT004 )%desc = 'Black Carbon Wet Deposition Bin 4'
      diag(iBCWT004 )%unit = 'kg/m2/s'

      diag(iBCWT005 )%name = 'BCWT005 '
      diag(iBCWT005 )%desc = 'Black Carbon Wet Deposition Bin 5'
      diag(iBCWT005 )%unit = 'kg/m2/s'

      diag(iBCWT006 )%name = 'BCWT006 '
      diag(iBCWT006 )%desc = 'Black Carbon Wet Deposition Bin 6'
      diag(iBCWT006 )%unit = 'kg/m2/s'

      diag(iBCWT007 )%name = 'BCWT007 '
      diag(iBCWT007 )%desc = 'Black Carbon Wet Deposition Bin 7'
      diag(iBCWT007 )%unit = 'kg/m2/s'

      diag(iBCWT008 )%name = 'BCWT008 '
      diag(iBCWT008 )%desc = 'Black Carbon Wet Deposition Bin 8'
      diag(iBCWT008 )%unit = 'kg/m2/s'

      diag(iBCSV001 )%name = 'BCSV001 '
      diag(iBCSV001 )%desc = 'Black Carbon Convective Scavenging Bin 1'
      diag(iBCSV001 )%unit = 'kg/m2/s'

      diag(iBCSV002 )%name = 'BCSV002 '
      diag(iBCSV002 )%desc = 'Black Carbon Convective Scavenging Bin 2'
      diag(iBCSV002 )%unit = 'kg/m2/s'

      diag(iBCSV003 )%name = 'BCSV003 '
      diag(iBCSV003 )%desc = 'Black Carbon Convective Scavenging Bin 3'
      diag(iBCSV003 )%unit = 'kg/m2/s'

      diag(iBCSV004 )%name = 'BCSV004 '
      diag(iBCSV004 )%desc = 'Black Carbon Convective Scavenging Bin 4'
      diag(iBCSV004 )%unit = 'kg/m2/s'

      diag(iBCSV005 )%name = 'BCSV005 '
      diag(iBCSV005 )%desc = 'Black Carbon Convective Scavenging Bin 5'
      diag(iBCSV005 )%unit = 'kg/m2/s'

      diag(iBCSV006 )%name = 'BCSV006 '
      diag(iBCSV006 )%desc = 'Black Carbon Convective Scavenging Bin 6'
      diag(iBCSV006 )%unit = 'kg/m2/s'

      diag(iBCSV007 )%name = 'BCSV007 '
      diag(iBCSV007 )%desc = 'Black Carbon Convective Scavenging Bin 7'
      diag(iBCSV007 )%unit = 'kg/m2/s'

      diag(iBCSV008 )%name = 'BCSV008 '
      diag(iBCSV008 )%desc = 'Black Carbon Convective Scavenging Bin 8'
      diag(iBCSV008 )%unit = 'kg/m2/s'

      diag(iBCSMASS )%name = 'BCSMASS '
      diag(iBCSMASS )%desc = 'Black Carbon Surface Mass Concentration'
      diag(iBCSMASS )%unit = 'kg/m3'

      diag(iBCCMASS )%name = 'BCCMASS '
      diag(iBCCMASS )%desc = 'Black Carbon Column Mass Density'
      diag(iBCCMASS )%unit = 'kg/m2'

      diag(iBCEXTTAU )%name = 'BCEXTTAU '
      diag(iBCEXTTAU )%desc = 'Black Carbon Extinction AOT [550 nm]'
      diag(iBCEXTTAU )%unit = 'unitless'

      diag(iBCSCATAU )%name = 'BCSCATAU '
      diag(iBCSCATAU )%desc = 'Black Carbon Scattering AOT [550 nm]'
      diag(iBCSCATAU )%unit = 'unitless'

      diag(iBCEMAN )%name = 'BCEMAN '
      diag(iBCEMAN )%desc = 'Black Carbon Anthropogenic Emissions'
      diag(iBCEMAN )%unit = 'kg/m2/s'

      diag(iBCEMBB )%name = 'BCEMBB '
      diag(iBCEMBB )%desc = 'Black Carbon Biomass Burning Emissions'
      diag(iBCEMBB )%unit = 'kg/m2/s'

      diag(iBCEMBF )%name = 'BCEMBF '
      diag(iBCEMBF )%desc = 'Black Carbon Biofuel Emissions'
      diag(iBCEMBF )%unit = 'kg/m2/s'

      diag(iBCHYPHIL )%name = 'BCHYPHIL '
      diag(iBCHYPHIL )%desc = 'Black Carbon Conversion of Hydrophobic to Hydrophilic'
      diag(iBCHYPHIL )%unit = 'kg/m2/s'

      diag(iOCEM001 )%name = 'OCEM001 '
      diag(iOCEM001 )%desc = 'Organic Carbon Emission Bin 1'
      diag(iOCEM001 )%unit = 'kg/m2/s'

      diag(iOCEM002 )%name = 'OCEM002 '
      diag(iOCEM002 )%desc = 'Organic Carbon Emission Bin 2'
      diag(iOCEM002 )%unit = 'kg/m2/s'

      diag(iOCEM003 )%name = 'OCEM003 '
      diag(iOCEM003 )%desc = 'Organic Carbon Emission Bin 3'
      diag(iOCEM003 )%unit = 'kg/m2/s'

      diag(iOCEM004 )%name = 'OCEM004 '
      diag(iOCEM004 )%desc = 'Organic Carbon Emission Bin 4'
      diag(iOCEM004 )%unit = 'kg/m2/s'

      diag(iOCEM005 )%name = 'OCEM005 '
      diag(iOCEM005 )%desc = 'Organic Carbon Emission Bin 5'
      diag(iOCEM005 )%unit = 'kg/m2/s'

      diag(iOCEM006 )%name = 'OCEM006 '
      diag(iOCEM006 )%desc = 'Organic Carbon Emission Bin 6'
      diag(iOCEM006 )%unit = 'kg/m2/s'

      diag(iOCEM007 )%name = 'OCEM007 '
      diag(iOCEM007 )%desc = 'Organic Carbon Emission Bin 7'
      diag(iOCEM007 )%unit = 'kg/m2/s'

      diag(iOCEM008 )%name = 'OCEM008 '
      diag(iOCEM008 )%desc = 'Organic Carbon Emission Bin 8'
      diag(iOCEM008 )%unit = 'kg/m2/s'

      diag(iOCDP001 )%name = 'OCDP001 '
      diag(iOCDP001 )%desc = 'Organic Carbon Deposition Bin 1'
      diag(iOCDP001 )%unit = 'kg/m2/s'

      diag(iOCDP002 )%name = 'OCDP002 '
      diag(iOCDP002 )%desc = 'Organic Carbon Deposition Bin 2'
      diag(iOCDP002 )%unit = 'kg/m2/s'

      diag(iOCDP003 )%name = 'OCDP003 '
      diag(iOCDP003 )%desc = 'Organic Carbon Deposition Bin 3'
      diag(iOCDP003 )%unit = 'kg/m2/s'

      diag(iOCDP004 )%name = 'OCDP004 '
      diag(iOCDP004 )%desc = 'Organic Carbon Deposition Bin 4'
      diag(iOCDP004 )%unit = 'kg/m2/s'

      diag(iOCDP005 )%name = 'OCDP005 '
      diag(iOCDP005 )%desc = 'Organic Carbon Deposition Bin 5'
      diag(iOCDP005 )%unit = 'kg/m2/s'

      diag(iOCDP006 )%name = 'OCDP006 '
      diag(iOCDP006 )%desc = 'Organic Carbon Deposition Bin 6'
      diag(iOCDP006 )%unit = 'kg/m2/s'

      diag(iOCDP007 )%name = 'OCDP007 '
      diag(iOCDP007 )%desc = 'Organic Carbon Deposition Bin 7'
      diag(iOCDP007 )%unit = 'kg/m2/s'

      diag(iOCDP008 )%name = 'OCDP008 '
      diag(iOCDP008 )%desc = 'Organic Carbon Deposition Bin 8'
      diag(iOCDP008 )%unit = 'kg/m2/s'

      diag(iOCWT001 )%name = 'OCWT001 '
      diag(iOCWT001 )%desc = 'Organic Carbon Wet Deposition Bin 1'
      diag(iOCWT001 )%unit = 'kg/m2/s'

      diag(iOCWT002 )%name = 'OCWT002 '
      diag(iOCWT002 )%desc = 'Organic Carbon Wet Deposition Bin 2'
      diag(iOCWT002 )%unit = 'kg/m2/s'

      diag(iOCWT003 )%name = 'OCWT003 '
      diag(iOCWT003 )%desc = 'Organic Carbon Wet Deposition Bin 3'
      diag(iOCWT003 )%unit = 'kg/m2/s'

      diag(iOCWT004 )%name = 'OCWT004 '
      diag(iOCWT004 )%desc = 'Organic Carbon Wet Deposition Bin 4'
      diag(iOCWT004 )%unit = 'kg/m2/s'

      diag(iOCWT005 )%name = 'OCWT005 '
      diag(iOCWT005 )%desc = 'Organic Carbon Wet Deposition Bin 5'
      diag(iOCWT005 )%unit = 'kg/m2/s'

      diag(iOCWT006 )%name = 'OCWT006 '
      diag(iOCWT006 )%desc = 'Organic Carbon Wet Deposition Bin 6'
      diag(iOCWT006 )%unit = 'kg/m2/s'

      diag(iOCWT007 )%name = 'OCWT007 '
      diag(iOCWT007 )%desc = 'Organic Carbon Wet Deposition Bin 7'
      diag(iOCWT007 )%unit = 'kg/m2/s'

      diag(iOCWT008 )%name = 'OCWT008 '
      diag(iOCWT008 )%desc = 'Organic Carbon Wet Deposition Bin 8'
      diag(iOCWT008 )%unit = 'kg/m2/s'

      diag(iOCSV001 )%name = 'OCSV001 '
      diag(iOCSV001 )%desc = 'Organic Carbon Convective Scavenging Bin 1'
      diag(iOCSV001 )%unit = 'kg/m2/s'

      diag(iOCSV002 )%name = 'OCSV002 '
      diag(iOCSV002 )%desc = 'Organic Carbon Convective Scavenging Bin 2'
      diag(iOCSV002 )%unit = 'kg/m2/s'

      diag(iOCSV003 )%name = 'OCSV003 '
      diag(iOCSV003 )%desc = 'Organic Carbon Convective Scavenging Bin 3'
      diag(iOCSV003 )%unit = 'kg/m2/s'

      diag(iOCSV004 )%name = 'OCSV004 '
      diag(iOCSV004 )%desc = 'Organic Carbon Convective Scavenging Bin 4'
      diag(iOCSV004 )%unit = 'kg/m2/s'

      diag(iOCSV005 )%name = 'OCSV005 '
      diag(iOCSV005 )%desc = 'Organic Carbon Convective Scavenging Bin 5'
      diag(iOCSV005 )%unit = 'kg/m2/s'

      diag(iOCSV006 )%name = 'OCSV006 '
      diag(iOCSV006 )%desc = 'Organic Carbon Convective Scavenging Bin 6'
      diag(iOCSV006 )%unit = 'kg/m2/s'

      diag(iOCSV007 )%name = 'OCSV007 '
      diag(iOCSV007 )%desc = 'Organic Carbon Convective Scavenging Bin 7'
      diag(iOCSV007 )%unit = 'kg/m2/s'

      diag(iOCSV008 )%name = 'OCSV008 '
      diag(iOCSV008 )%desc = 'Organic Carbon Convective Scavenging Bin 8'
      diag(iOCSV008 )%unit = 'kg/m2/s'

      diag(iOCSMASS )%name = 'OCSMASS '
      diag(iOCSMASS )%desc = 'Organic Carbon Surface Mass Concentration'
      diag(iOCSMASS )%unit = 'kg/m3'

      diag(iOCCMASS )%name = 'OCCMASS '
      diag(iOCCMASS )%desc = 'Organic Carbon Column Mass Density'
      diag(iOCCMASS )%unit = 'kg/m2'

      diag(iOCEXTTAU )%name = 'OCEXTTAU '
      diag(iOCEXTTAU )%desc = 'Organic Carbon Extinction AOT [550 nm]'
      diag(iOCEXTTAU )%unit = 'unitless'

      diag(iOCSCATAU )%name = 'OCSCATAU '
      diag(iOCSCATAU )%desc = 'Organic Carbon Scattering AOT [550 nm]'
      diag(iOCSCATAU )%unit = 'unitless'

      diag(iOCEMAN )%name = 'OCEMAN '
      diag(iOCEMAN )%desc = 'Organic Carbon Anthropogenic Emissions'
      diag(iOCEMAN )%unit = 'kg/m2/s'

      diag(iOCEMBB )%name = 'OCEMBB '
      diag(iOCEMBB )%desc = 'Organic Carbon Biomass Burning Emissions'
      diag(iOCEMBB )%unit = 'kg/m2/s'

      diag(iOCEMBF )%name = 'OCEMBF '
      diag(iOCEMBF )%desc = 'Organic Carbon Biofuel Emissions'
      diag(iOCEMBF )%unit = 'kg/m2/s'

      diag(iOCEMBG )%name = 'OCEMBG '
      diag(iOCEMBG )%desc = 'Organic Carbon Biogenic Emissions'
      diag(iOCEMBG )%unit = 'kg/m2/s'

      diag(iOCHYPHIL )%name = 'OCHYPHIL '
      diag(iOCHYPHIL )%desc = 'Organic Carbon Conversion of Hydrophobic to Hydrophilic'
      diag(iOCHYPHIL )%unit = 'kg/m2/s'

      diag(iSUEM001 )%name = 'SUEM001 '
      diag(iSUEM001 )%desc = 'Sulfate Emission Bin 1'
      diag(iSUEM001 )%unit = 'kg/m2/s'

      diag(iSUEM002 )%name = 'SUEM002 '
      diag(iSUEM002 )%desc = 'Sulfate Emission Bin 2'
      diag(iSUEM002 )%unit = 'kg/m2/s'

      diag(iSUEM003 )%name = 'SUEM003 '
      diag(iSUEM003 )%desc = 'Sulfate Emission Bin 3'
      diag(iSUEM003 )%unit = 'kg/m2/s'

      diag(iSUEM004 )%name = 'SUEM004 '
      diag(iSUEM004 )%desc = 'Sulfate Emission Bin 4'
      diag(iSUEM004 )%unit = 'kg/m2/s'

      diag(iSUEM005 )%name = 'SUEM005 '
      diag(iSUEM005 )%desc = 'Sulfate Emission Bin 5'
      diag(iSUEM005 )%unit = 'kg/m2/s'

      diag(iSUEM006 )%name = 'SUEM006 '
      diag(iSUEM006 )%desc = 'Sulfate Emission Bin 6'
      diag(iSUEM006 )%unit = 'kg/m2/s'

      diag(iSUEM007 )%name = 'SUEM007 '
      diag(iSUEM007 )%desc = 'Sulfate Emission Bin 7'
      diag(iSUEM007 )%unit = 'kg/m2/s'

      diag(iSUEM008 )%name = 'SUEM008 '
      diag(iSUEM008 )%desc = 'Sulfate Emission Bin 8'
      diag(iSUEM008 )%unit = 'kg/m2/s'

      diag(iSUDP001 )%name = 'SUDP001 '
      diag(iSUDP001 )%desc = 'Sulfate Deposition Bin 1'
      diag(iSUDP001 )%unit = 'kg/m2/s'

      diag(iSUDP002 )%name = 'SUDP002 '
      diag(iSUDP002 )%desc = 'Sulfate Deposition Bin 2'
      diag(iSUDP002 )%unit = 'kg/m2/s'

      diag(iSUDP003 )%name = 'SUDP003 '
      diag(iSUDP003 )%desc = 'Sulfate Deposition Bin 3'
      diag(iSUDP003 )%unit = 'kg/m2/s'

      diag(iSUDP004 )%name = 'SUDP004 '
      diag(iSUDP004 )%desc = 'Sulfate Deposition Bin 4'
      diag(iSUDP004 )%unit = 'kg/m2/s'

      diag(iSUDP005 )%name = 'SUDP005 '
      diag(iSUDP005 )%desc = 'Sulfate Deposition Bin 5'
      diag(iSUDP005 )%unit = 'kg/m2/s'

      diag(iSUDP006 )%name = 'SUDP006 '
      diag(iSUDP006 )%desc = 'Sulfate Deposition Bin 6'
      diag(iSUDP006 )%unit = 'kg/m2/s'

      diag(iSUDP007 )%name = 'SUDP007 '
      diag(iSUDP007 )%desc = 'Sulfate Deposition Bin 7'
      diag(iSUDP007 )%unit = 'kg/m2/s'

      diag(iSUDP008 )%name = 'SUDP008 '
      diag(iSUDP008 )%desc = 'Sulfate Deposition Bin 8'
      diag(iSUDP008 )%unit = 'kg/m2/s'

      diag(iSUWT001 )%name = 'SUWT001 '
      diag(iSUWT001 )%desc = 'Sulfate Wet Deposition Bin 1'
      diag(iSUWT001 )%unit = 'kg/m2/s'

      diag(iSUWT002 )%name = 'SUWT002 '
      diag(iSUWT002 )%desc = 'Sulfate Wet Deposition Bin 2'
      diag(iSUWT002 )%unit = 'kg/m2/s'

      diag(iSUWT003 )%name = 'SUWT003 '
      diag(iSUWT003 )%desc = 'Sulfate Wet Deposition Bin 3'
      diag(iSUWT003 )%unit = 'kg/m2/s'

      diag(iSUWT004 )%name = 'SUWT004 '
      diag(iSUWT004 )%desc = 'Sulfate Wet Deposition Bin 4'
      diag(iSUWT004 )%unit = 'kg/m2/s'

      diag(iSUWT005 )%name = 'SUWT005 '
      diag(iSUWT005 )%desc = 'Sulfate Wet Deposition Bin 5'
      diag(iSUWT005 )%unit = 'kg/m2/s'

      diag(iSUWT006 )%name = 'SUWT006 '
      diag(iSUWT006 )%desc = 'Sulfate Wet Deposition Bin 6'
      diag(iSUWT006 )%unit = 'kg/m2/s'

      diag(iSUWT007 )%name = 'SUWT007 '
      diag(iSUWT007 )%desc = 'Sulfate Wet Deposition Bin 7'
      diag(iSUWT007 )%unit = 'kg/m2/s'

      diag(iSUWT008 )%name = 'SUWT008 '
      diag(iSUWT008 )%desc = 'Sulfate Wet Deposition Bin 8'
      diag(iSUWT008 )%unit = 'kg/m2/s'

      diag(iSUSV001 )%name = 'SUSV001 '
      diag(iSUSV001 )%desc = 'Sulfate Convective Scavenging Bin 1'
      diag(iSUSV001 )%unit = 'kg/m2/s'

      diag(iSUSV002 )%name = 'SUSV002 '
      diag(iSUSV002 )%desc = 'Sulfate Convective Scavenging Bin 2'
      diag(iSUSV002 )%unit = 'kg/m2/s'

      diag(iSUSV003 )%name = 'SUSV003 '
      diag(iSUSV003 )%desc = 'Sulfate Convective Scavenging Bin 3'
      diag(iSUSV003 )%unit = 'kg/m2/s'

      diag(iSUSV004 )%name = 'SUSV004 '
      diag(iSUSV004 )%desc = 'Sulfate Convective Scavenging Bin 4'
      diag(iSUSV004 )%unit = 'kg/m2/s'

      diag(iSUSV005 )%name = 'SUSV005 '
      diag(iSUSV005 )%desc = 'Sulfate Convective Scavenging Bin 5'
      diag(iSUSV005 )%unit = 'kg/m2/s'

      diag(iSUSV006 )%name = 'SUSV006 '
      diag(iSUSV006 )%desc = 'Sulfate Convective Scavenging Bin 6'
      diag(iSUSV006 )%unit = 'kg/m2/s'

      diag(iSUSV007 )%name = 'SUSV007 '
      diag(iSUSV007 )%desc = 'Sulfate Convective Scavenging Bin 7'
      diag(iSUSV007 )%unit = 'kg/m2/s'

      diag(iSUSV008 )%name = 'SUSV008 '
      diag(iSUSV008 )%desc = 'Sulfate Convective Scavenging Bin 8'
      diag(iSUSV008 )%unit = 'kg/m2/s'

      diag(iSUSO2SMASS )%name = 'SO2SMASS '
      diag(iSUSO2SMASS )%desc = 'SO2 Surface Mass Concentration'
      diag(iSUSO2SMASS )%unit = 'kg/m3'

      diag(iSUSO2CMASS )%name = 'SO2CMASS '
      diag(iSUSO2CMASS )%desc = 'SO2 Column Mass Density'
      diag(iSUSO2CMASS )%unit = 'kg/m2'

      diag(iSUSO4SMASS )%name = 'SO4SMASS '
      diag(iSUSO4SMASS )%desc = 'SO4 Surface Mass Concentration'
      diag(iSUSO4SMASS )%unit = 'kg/m3'

      diag(iSUSO4CMASS )%name = 'SO4CMASS '
      diag(iSUSO4CMASS )%desc = 'SO4 Column Mass Density'
      diag(iSUSO4CMASS )%unit = 'kg/m2'

      diag(iSUDMSSMASS )%name = 'DMSSMASS '
      diag(iSUDMSSMASS )%desc = 'SO2 Surface Mass Concentration'
      diag(iSUDMSSMASS )%unit = 'kg/m3'

      diag(iSUDMSCMASS )%name = 'DMSCMASS '
      diag(iSUDMSCMASS )%desc = 'SO2 Column Mass Density'
      diag(iSUDMSCMASS )%unit = 'kg/m2'

      diag(iSUPSO2 )%name = 'SUPSO2 '
      diag(iSUPSO2 )%desc = 'SO2 Production from DMS Oxidation column integrated'
      diag(iSUPSO2 )%unit = 'kg/m2/s'

      diag(iSUPSO4g )%name = 'SUPSO4g '
      diag(iSUPSO4g )%desc = 'SO4 Production from gas-phase SO2 Oxidation column integrated'
      diag(iSUPSO4g )%unit = 'kg/m2/s'

      diag(iSUPSO4aq )%name = 'SUPSO4aq '
      diag(iSUPSO4aq )%desc = 'SO4 Production from aqueous SO2 Oxidation column integrated'
      diag(iSUPSO4aq )%unit = 'kg/m2/s'

      diag(iSUPSO4wet )%name = 'SUPSO4wt '
      diag(iSUPSO4wet )%desc = 'SO4 Production from aqueous SO2 Oxidation (wet dep) column integrated'
      diag(iSUPSO4wet )%unit = 'kg/m2/s'

      diag(iSUPMSA )%name = 'SUPMSA '
      diag(iSUPMSA )%desc = 'MSA Production from DMS Oxidation column integrated'
      diag(iSUPMSA )%unit = 'kg/m2/s'

      diag(iSUEXTTAU )%name = 'SUEXTTAU '
      diag(iSUEXTTAU )%desc = 'SO4 Extinction AOT [550 nm]'
      diag(iSUEXTTAU )%unit = 'unitless'

      diag(iSUSCATAU )%name = 'SUSCATAU '
      diag(iSUSCATAU )%desc = 'SO4 Scattering AOT [550 nm]'
      diag(iSUSCATAU )%unit = 'unitless'

      diag(iSUEMSO4AN )%name = 'SO4EMAN '
      diag(iSUEMSO4AN )%desc = 'SO4 Anthropogenic Emissions'
      diag(iSUEMSO4AN )%unit = 'kg/m2/s'

      diag(iSUEMSO2AN )%name = 'SO2EMAN '
      diag(iSUEMSO2AN )%desc = 'SO2 Anthropogenic Emissions'
      diag(iSUEMSO2AN )%unit = 'kg/m2/s'

      diag(iSUEMSO2BB )%name = 'SO2EMBB '
      diag(iSUEMSO2BB )%desc = 'SO2 Biomass Burning Emissions'
      diag(iSUEMSO2BB )%unit = 'kg/m2/s'

      diag(iSUEMSO2VN )%name = 'SO2EMVN '
      diag(iSUEMSO2VN )%desc = 'SO2 Non-explosive Volcanic Emissions'
      diag(iSUEMSO2VN )%unit = 'kg/m2/s'

      diag(iSUEMSO2VE )%name = 'SO2EMVE '
      diag(iSUEMSO2VE )%desc = 'SO2 Explosive Volanic Emissions'
      diag(iSUEMSO2VE )%unit = 'kg/m2/s'

      diag(iN2OFLX )%name = 'N2OFLX'
      diag(iN2OFLX )%desc = 'Derived N2O surface flux'
      diag(iN2OFLX )%unit = 'kg/s'

      diag(iCH4FLX )%name = 'CH4FLX'
      diag(iCH4FLX )%desc = 'Derived CH4 surface flux'
      diag(iCH4FLX )%unit = 'kg/s'

      diag(iF11FLX )%name = 'F11FLX'
      diag(iF11FLX )%desc = 'Derived F11 surface flux'
      diag(iF11FLX )%unit = 'kg/s'

      diag(iF12FLX )%name = 'F12FLX'
      diag(iF12FLX )%desc = 'Derived F12 surface flux'
      diag(iF12FLX )%unit = 'kg/s'

      diag(iF113FLX )%name = 'F113FLX'
      diag(iF113FLX )%desc = 'Derived F113 surface flux'
      diag(iF113FLX )%unit = 'kg/s'

      diag(iHCFCFLX )%name = 'HCFCFLX'
      diag(iHCFCFLX )%desc = 'Derived HCFC surface flux'
      diag(iHCFCFLX )%unit = 'kg/s'

      diag(iCCL4FLX )%name = 'CCL4FLX'
      diag(iCCL4FLX )%desc = 'Derived CCL4 surface flux'
      diag(iCCL4FLX )%unit = 'kg/s'

      diag(iCH3CCL3FLX )%name = 'MCFFLX'
      diag(iCH3CCL3FLX )%desc = 'Derived CH3CCL3 surface flux'
      diag(iCH3CCL3FLX )%unit = 'kg/s'

      diag(iCH3CLFLX )%name = 'CH3CLFLX'
      diag(iCH3CLFLX )%desc = 'Derived CH3CL surface flux'
      diag(iCH3CLFLX )%unit = 'kg/s'

      diag(iCH3BRFLX )%name = 'CH3BRFLX'
      diag(iCH3BRFLX )%desc = 'Derived CH3BR surface flux'
      diag(iCH3BRFLX )%unit = 'kg/s'

      diag(iH1301FLX )%name = 'H1301FLX'
      diag(iH1301FLX )%desc = 'Derived H1301 surface flux'
      diag(iH1301FLX )%unit = 'kg/s'

      diag(iH12_24FLX )%name = 'H1224FLX'
      diag(iH12_24FLX )%desc = 'Derived H12_24 surface flux'
      diag(iH12_24FLX )%unit = 'kg/s'

      diag(iCOSSZA )%name = 'COSSZA '
      diag(iCOSSZA )%desc = 'Cosine of SZA'
      diag(iCOSSZA )%unit = 'unitless'

      diag(iTCOSZ )%name = 'TCOSZ '
      diag(iTCOSZ )%desc = 'Integrated cosine of SZA'
      diag(iTCOSZ )%unit = 'unitless'

      diag(iXOH )%name = 'XOH '
      diag(iXOH )%desc = 'OH surface concentration'
      diag(iXOH )%unit = '# cm-3'

      diag(iXNO3 )%name = 'XNO3 '
      diag(iXNO3 )%desc = 'NO3 surface mixing ratio'
      diag(iXNO3 )%unit = 'unitless'

      diag(iXH2O2 )%name = 'XH2O2 '
      diag(iXH2O2 )%desc = 'H2O2 surface mixing ratio'
      diag(iXH2O2 )%unit = 'unitless'

      diag(iCOEM001 )%name = 'COEM001 '
      diag(iCOEM001 )%desc = 'Carbon Monoxide Emission Bin 1'
      diag(iCOEM001 )%unit = 'kg/m2/s'

      diag(iCOEM002 )%name = 'COEM002 '
      diag(iCOEM002 )%desc = 'Carbon Monoxide Emission Bin 2'
      diag(iCOEM002 )%unit = 'kg/m2/s'

      diag(iCOEM003 )%name = 'COEM003 '
      diag(iCOEM003 )%desc = 'Carbon Monoxide Emission Bin 3'
      diag(iCOEM003 )%unit = 'kg/m2/s'

      diag(iCOEM004 )%name = 'COEM004 '
      diag(iCOEM004 )%desc = 'Carbon Monoxide Emission Bin 4'
      diag(iCOEM004 )%unit = 'kg/m2/s'

      diag(iCOEM005 )%name = 'COEM005 '
      diag(iCOEM005 )%desc = 'Carbon Monoxide Emission Bin 5'
      diag(iCOEM005 )%unit = 'kg/m2/s'

      diag(iCOEM006 )%name = 'COEM006 '
      diag(iCOEM006 )%desc = 'Carbon Monoxide Emission Bin 6'
      diag(iCOEM006 )%unit = 'kg/m2/s'

      diag(iCOEM007 )%name = 'COEM007 '
      diag(iCOEM007 )%desc = 'Carbon Monoxide Emission Bin 7'
      diag(iCOEM007 )%unit = 'kg/m2/s'

      diag(iCOEM008 )%name = 'COEM008 '
      diag(iCOEM008 )%desc = 'Carbon Monoxide Emission Bin 8'
      diag(iCOEM008 )%unit = 'kg/m2/s'

      diag(iCOLS001 )%name = 'COLS001 '
      diag(iCOLS001 )%desc = 'Carbon Monoxide Chemical Loss Bin 1'
      diag(iCOLS001 )%unit = 'kg/m2/s'

      diag(iCOLS002 )%name = 'COLS002 '
      diag(iCOLS002 )%desc = 'Carbon Monoxide Chemical Loss Bin 2'
      diag(iCOLS002 )%unit = 'kg/m2/s'

      diag(iCOLS003 )%name = 'COLS003 '
      diag(iCOLS003 )%desc = 'Carbon Monoxide Chemical Loss Bin 3'
      diag(iCOLS003 )%unit = 'kg/m2/s'

      diag(iCOLS004 )%name = 'COLS004 '
      diag(iCOLS004 )%desc = 'Carbon Monoxide Chemical Loss Bin 4'
      diag(iCOLS004 )%unit = 'kg/m2/s'

      diag(iCOLS005 )%name = 'COLS005 '
      diag(iCOLS005 )%desc = 'Carbon Monoxide Chemical Loss Bin 5'
      diag(iCOLS005 )%unit = 'kg/m2/s'

      diag(iCOLS006 )%name = 'COLS006 '
      diag(iCOLS006 )%desc = 'Carbon Monoxide Chemical Loss Bin 6'
      diag(iCOLS006 )%unit = 'kg/m2/s'

      diag(iCOLS007 )%name = 'COLS007 '
      diag(iCOLS007 )%desc = 'Carbon Monoxide Chemical Loss Bin 7'
      diag(iCOLS007 )%unit = 'kg/m2/s'

      diag(iCOLS008 )%name = 'COLS008 '
      diag(iCOLS008 )%desc = 'Carbon Monoxide Chemical Loss Bin 8'
      diag(iCOLS008 )%unit = 'kg/m2/s'

      diag(iCOPD001 )%name = 'COPD001 '
      diag(iCOPD001 )%desc = 'Carbon Monoxide Chemical Production Bin 1'
      diag(iCOPD001 )%unit = 'kg/m2/s'

      diag(iCOPD002 )%name = 'COPD002 '
      diag(iCOPD002 )%desc = 'Carbon Monoxide Chemical Production Bin 2'
      diag(iCOPD002 )%unit = 'kg/m2/s'

      diag(iCOPD003 )%name = 'COPD003 '
      diag(iCOPD003 )%desc = 'Carbon Monoxide Chemical Production Bin 3'
      diag(iCOPD003 )%unit = 'kg/m2/s'

      diag(iCOPD004 )%name = 'COPD004 '
      diag(iCOPD004 )%desc = 'Carbon Monoxide Chemical Production Bin 4'
      diag(iCOPD004 )%unit = 'kg/m2/s'

      diag(iCOPD005 )%name = 'COPD005 '
      diag(iCOPD005 )%desc = 'Carbon Monoxide Chemical Production Bin 5'
      diag(iCOPD005 )%unit = 'kg/m2/s'

      diag(iCOPD006 )%name = 'COPD006 '
      diag(iCOPD006 )%desc = 'Carbon Monoxide Chemical Production Bin 6'
      diag(iCOPD006 )%unit = 'kg/m2/s'

      diag(iCOPD007 )%name = 'COPD007 '
      diag(iCOPD007 )%desc = 'Carbon Monoxide Chemical Production Bin 7'
      diag(iCOPD007 )%unit = 'kg/m2/s'

      diag(iCOPD008 )%name = 'COPD008 '
      diag(iCOPD008 )%desc = 'Carbon Monoxide Chemical Production Bin 8'
      diag(iCOPD008 )%unit = 'kg/m2/s'

      diag(iCOCL001 )%name = 'COCL001 '
      diag(iCOCL001 )%desc = 'Carbon Monoxide Column Burden Bin 1'
      diag(iCOCL001 )%unit = 'kg/m2'

      diag(iCOCL002 )%name = 'COCL002 '
      diag(iCOCL002 )%desc = 'Carbon Monoxide Column Burden Bin 2'
      diag(iCOCL002 )%unit = 'kg/m2'

      diag(iCOCL003 )%name = 'COCL003 '
      diag(iCOCL003 )%desc = 'Carbon Monoxide Column Burden Bin 3'
      diag(iCOCL003 )%unit = 'kg/m2'

      diag(iCOCL004 )%name = 'COCL004 '
      diag(iCOCL004 )%desc = 'Carbon Monoxide Column Burden Bin 4'
      diag(iCOCL004 )%unit = 'kg/m2'

      diag(iCOCL005 )%name = 'COCL005 '
      diag(iCOCL005 )%desc = 'Carbon Monoxide Column Burden Bin 5'
      diag(iCOCL005 )%unit = 'kg/m2'

      diag(iCOCL006 )%name = 'COCL006 '
      diag(iCOCL006 )%desc = 'Carbon Monoxide Column Burden Bin 6'
      diag(iCOCL006 )%unit = 'kg/m2'

      diag(iCOCL007 )%name = 'COCL007 '
      diag(iCOCL007 )%desc = 'Carbon Monoxide Column Burden Bin 7'
      diag(iCOCL007 )%unit = 'kg/m2'

      diag(iCOCL008 )%name = 'COCL008 '
      diag(iCOCL008 )%desc = 'Carbon Monoxide Column Burden Bin 8'
      diag(iCOCL008 )%unit = 'kg/m2'

      diag(iCOSC001 )%name = 'COSC001 '
      diag(iCOSC001 )%desc = 'Carbon Monoxide Surface Concentration Bin 1'
      diag(iCOSC001 )%unit = 'ppbv'

      diag(iCOSC002 )%name = 'COSC002 '
      diag(iCOSC002 )%desc = 'Carbon Monoxide Surface Concentration Bin 2'
      diag(iCOSC002 )%unit = 'ppbv'

      diag(iCOSC003 )%name = 'COSC003 '
      diag(iCOSC003 )%desc = 'Carbon Monoxide Surface Concentration Bin 3'
      diag(iCOSC003 )%unit = 'ppbv'

      diag(iCOSC004 )%name = 'COSC004 '
      diag(iCOSC004 )%desc = 'Carbon Monoxide Surface Concentration Bin 4'
      diag(iCOSC004 )%unit = 'ppbv'

      diag(iCOSC005 )%name = 'COSC005 '
      diag(iCOSC005 )%desc = 'Carbon Monoxide Surface Concentration Bin 5'
      diag(iCOSC005 )%unit = 'ppbv'

      diag(iCOSC006 )%name = 'COSC006 '
      diag(iCOSC006 )%desc = 'Carbon Monoxide Surface Concentration Bin 6'
      diag(iCOSC006 )%unit = 'ppbv'

      diag(iCOSC007 )%name = 'COSC007 '
      diag(iCOSC007 )%desc = 'Carbon Monoxide Surface Concentration Bin 7'
      diag(iCOSC007 )%unit = 'ppbv'

      diag(iCOSC008 )%name = 'COSC008 '
      diag(iCOSC008 )%desc = 'Carbon Monoxide Surface Concentration Bin 8'
      diag(iCOSC008 )%unit = 'ppbv'

      diag(iCO2EM001 )%name = 'CO2EM001 '
      diag(iCO2EM001 )%desc = 'Carbon Dioxide Emission Bin 1'
      diag(iCO2EM001 )%unit = 'kg/m2/s'

      diag(iCO2EM002 )%name = 'CO2EM002 '
      diag(iCO2EM002 )%desc = 'Carbon Dioxide Emission Bin 2'
      diag(iCO2EM002 )%unit = 'kg/m2/s'

      diag(iCO2EM003 )%name = 'CO2EM003 '
      diag(iCO2EM003 )%desc = 'Carbon Dioxide Emission Bin 3'
      diag(iCO2EM003 )%unit = 'kg/m2/s'

      diag(iCO2EM004 )%name = 'CO2EM004 '
      diag(iCO2EM004 )%desc = 'Carbon Dioxide Emission Bin 4'
      diag(iCO2EM004 )%unit = 'kg/m2/s'

      diag(iCO2EM005 )%name = 'CO2EM005 '
      diag(iCO2EM005 )%desc = 'Carbon Dioxide Emission Bin 5'
      diag(iCO2EM005 )%unit = 'kg/m2/s'

      diag(iCO2EM006 )%name = 'CO2EM006 '
      diag(iCO2EM006 )%desc = 'Carbon Dioxide Emission Bin 6'
      diag(iCO2EM006 )%unit = 'kg/m2/s'

      diag(iCO2EM007 )%name = 'CO2EM007 '
      diag(iCO2EM007 )%desc = 'Carbon Dioxide Emission Bin 7'
      diag(iCO2EM007 )%unit = 'kg/m2/s'

      diag(iCO2EM008 )%name = 'CO2EM008 '
      diag(iCO2EM008 )%desc = 'Carbon Dioxide Emission Bin 8'
      diag(iCO2EM008 )%unit = 'kg/m2/s'

      diag(iCO2CL001 )%name = 'CO2CL001 '
      diag(iCO2CL001 )%desc = 'Carbon Dioxide Column Burden Bin 1'
      diag(iCO2CL001 )%unit = 'kg/m2'

      diag(iCO2CL002 )%name = 'CO2CL002 '
      diag(iCO2CL002 )%desc = 'Carbon Dioxide Column Burden Bin 2'
      diag(iCO2CL002 )%unit = 'kg/m2'

      diag(iCO2CL003 )%name = 'CO2CL003 '
      diag(iCO2CL003 )%desc = 'Carbon Dioxide Column Burden Bin 3'
      diag(iCO2CL003 )%unit = 'kg/m2'

      diag(iCO2CL004 )%name = 'CO2CL004 '
      diag(iCO2CL004 )%desc = 'Carbon Dioxide Column Burden Bin 4'
      diag(iCO2CL004 )%unit = 'kg/m2'

      diag(iCO2CL005 )%name = 'CO2CL005 '
      diag(iCO2CL005 )%desc = 'Carbon Dioxide Column Burden Bin 5'
      diag(iCO2CL005 )%unit = 'kg/m2'

      diag(iCO2CL006 )%name = 'CO2CL006 '
      diag(iCO2CL006 )%desc = 'Carbon Dioxide Column Burden Bin 6'
      diag(iCO2CL006 )%unit = 'kg/m2'

      diag(iCO2CL007 )%name = 'CO2CL007 '
      diag(iCO2CL007 )%desc = 'Carbon Dioxide Column Burden Bin 7'
      diag(iCO2CL007 )%unit = 'kg/m2'

      diag(iCO2CL008 )%name = 'CO2CL008 '
      diag(iCO2CL008 )%desc = 'Carbon Dioxide Column Burden Bin 8'
      diag(iCO2CL008 )%unit = 'kg/m2'

      diag(iCO2SC001 )%name = 'CO2SC001 '
      diag(iCO2SC001 )%desc = 'Carbon Dioxide Surface Concentration Bin 1'
      diag(iCO2SC001 )%unit = 'ppmv'

      diag(iCO2SC002 )%name = 'CO2SC002 '
      diag(iCO2SC002 )%desc = 'Carbon Dioxide Surface Concentration Bin 2'
      diag(iCO2SC002 )%unit = 'ppmv'

      diag(iCO2SC003 )%name = 'CO2SC003 '
      diag(iCO2SC003 )%desc = 'Carbon Dioxide Surface Concentration Bin 3'
      diag(iCO2SC003 )%unit = 'ppmv'

      diag(iCO2SC004 )%name = 'CO2SC004 '
      diag(iCO2SC004 )%desc = 'Carbon Dioxide Surface Concentration Bin 4'
      diag(iCO2SC004 )%unit = 'ppmv'

      diag(iCO2SC005 )%name = 'CO2SC005 '
      diag(iCO2SC005 )%desc = 'Carbon Dioxide Surface Concentration Bin 5'
      diag(iCO2SC005 )%unit = 'ppmv'

      diag(iCO2SC006 )%name = 'CO2SC006 '
      diag(iCO2SC006 )%desc = 'Carbon Dioxide Surface Concentration Bin 6'
      diag(iCO2SC006 )%unit = 'ppmv'

      diag(iCO2SC007 )%name = 'CO2SC007 '
      diag(iCO2SC007 )%desc = 'Carbon Dioxide Surface Concentration Bin 7'
      diag(iCO2SC007 )%unit = 'ppmv'

      diag(iCO2SC008 )%name = 'CO2SC008 '
      diag(iCO2SC008 )%desc = 'Carbon Dioxide Surface Concentration Bin 8'
      diag(iCO2SC008 )%unit = 'ppmv'


#endif

! ... 3-D diagnostic fields (in alphabetic order)
!

      diag(iAIRDENS )%name = 'AIRDENS '
      diag(iAIRDENS )%desc = 'Air Density'
      diag(iAIRDENS )%unit = 'Kg m-3'

      diag(iCAPE    )%name = 'CAPE    '
      diag(iCAPE    )%desc = 'CAPE'
      diag(iCAPE    )%unit = 'J/kg'

      diag(iCGS     )%name = 'CGS     '
      diag(iCGS     )%desc = 'Counter-gradient coefficient on surface '  &
                         // 'kinematic fluxes'
      diag(iCGS     )%unit = 's/m2'
 
      diag(iCLDLWP  )%name = 'CLDLWP  '
      diag(iCLDLWP  )%desc = 'Actual cloud liquid water path length'
      diag(iCLDLWP  )%unit = 'gram/m2'

      diag(iCLOUD   )%name = 'CLOUD   '
      diag(iCLOUD   )%desc = 'Cloud fraction'
      diag(iCLOUD   )%unit = 'fraction'
      diag(iCLOUD   )%aname= 'CLDTOT  '
      diag(iCLOUD   )%adesc= '3-D Total Cloud Fraction'

      diag(iCLOUDUP )%name = 'CLOUDUP '
      diag(iCLOUDUP )%desc = 'Cloud fraction during omega < 0.'
      diag(iCLOUDUP )%unit = 'fraction'
 
      diag(iCMFDQ   )%name = 'CMFDQ   '
      diag(iCMFDQ   )%desc = 'q tendency - Hack convection'
      diag(iCMFDQ   )%unit = 'kg/kg/s'

      diag(iCMFDQR2 )%name = 'CMFDQR2 '
      diag(iCMFDQR2 )%desc = 'Rain production rate - Hack convection'
      diag(iCMFDQR2 )%unit = 'kg/kg/s'

      diag(iCMFDT   )%name = 'CMFDT   '
      diag(iCMFDT   )%desc = 'T tendency - Hack convection'
      diag(iCMFDT   )%unit = 'K/s'

      diag(iCMFDTR  )%name = 'CMFDTR  '
      diag(iCMFDTR  )%desc = 'Detrainment mass flux - Hack convection'
      diag(iCMFDTR  )%unit = 'Pa/s'

      diag(iCMFETR  )%name = 'CMFETR  '
      diag(iCMFETR  )%desc = 'Entrainment mass flux - Hack convection'
      diag(iCMFETR  )%unit = 'Pa/s'

      diag(iCMFMC   )%name = 'CMFMC   '
      diag(iCMFMC   )%desc = 'Total Moist convection mass flux'
      diag(iCMFMC   )%unit = 'Pa/s'
      diag(iCMFMC   )%aname= 'CLDMAS  '
      diag(iCMFMC   )%adesc= 'Cloud mass flux'
      diag(iCMFMC   )%aunit= 'kg/m2/s'
      diag(iCMFMC   )%convfac = 1./9.80616	! should use constant rga
!     diag(iCMFMC   )%convfac = rga       	! include comcon when f90

      diag(iCMFMC2  )%name = 'CMFMC2  '
      diag(iCMFMC2  )%desc = 'Hack Moist convection mass flux'
      diag(iCMFMC2  )%unit = 'Pa/s'

      diag(iCONVCLD  )%name = 'CONVCLD  '
      diag(iCONVCLD  )%desc = 'Convective cloud amount'
      diag(iCONVCLD  )%unit = 'fraction'

      diag(iDCAFDT  )%name = 'DCAFDT  '
      diag(iDCAFDT  )%desc = 'T Tendency - Dry convective adjustment'
      diag(iDCAFDT  )%unit = 'K/s'

      diag(iDIABDT  )%name = 'DIABDT  '
      diag(iDIABDT  )%desc = 'T Tendency - Total adiabatic (physics)'
      diag(iDIABDT  )%unit = 'K/s'                   
      diag(iDIABDT  )%aname= 'DIABT   '
      diag(iDIABDT  )%adesc= 'Temperature Tendency due to Diabatic Forcing'
      diag(iDIABDT  )%aunit= 'K/day'                   
      diag(iDIABDT  )%convfac = 86400.0

      diag(iDQCOND  )%name = 'DQCOND  '
      diag(iDQCOND  )%desc = 'Q tendency - moist physics'
      diag(iDQCOND  )%unit = 'kg/kg/s'
      diag(iDQCOND  )%aname= 'MOISTQ  '
      diag(iDQCOND  )%adesc= 'Specific Humidity Tendency due to Moist Processes'
      diag(iDQCOND  )%aunit= 'g/kg/day'
      diag(iDQCOND  )%convfac = 86400000.0

      diag(iDQPBLCG )%name = 'DQPBLCG '
      diag(iDQPBLCG )%desc = 'Q tendency - PBL counter gradient'
      diag(iDQPBLCG )%unit = 'kg/kg/s'

      diag(iDQRL    )%name = 'DQRL    '
      diag(iDQRL    )%desc = 'Rain production rate - large-scale'
      diag(iDQRL    )%unit = 'kg/kg/s'

      diag(iDRUNOFF )%name = 'DRUNOFF '
      diag(iDRUNOFF )%desc = 'Sub-surface runoff'
      diag(iDRUNOFF )%unit = 'mm/s'

      diag(iDTCOND  )%name = 'DTCOND  '
      diag(iDTCOND  )%desc = 'T tendency - moist physics'
      diag(iDTCOND  )%unit = 'K/s'
      diag(iDTCOND  )%aname= 'MOISTT  '
      diag(iDTCOND  )%adesc= 'Temperature Tendency due to Moist Processes'
      diag(iDTCOND  )%aunit= 'K/day'
      diag(iDTCOND  )%convfac = 86400.0

      diag(iDTPBLCG )%name = 'DTPBLCG '
      diag(iDTPBLCG )%desc = 'T tendency - PBL counter gradient'
      diag(iDTPBLCG )%unit = 'K/s'
 
      diag(iDTRAIN  )%name = 'DTRAIN  '
      diag(iDTRAIN  )%desc = 'Detrainment Cloud Mass Flux'
      diag(iDTRAIN  )%unit = 'Pa/s'
      diag(iDTRAIN  )%aunit= 'kg/m2/s'
      diag(iDTRAIN  )%convfac = 1./9.80616	! should use constant rga
!     diag(iDTRAIN  )%convfac = rga       	! include comcon when f90

      diag(iDTV     )%name = 'DTV     '
      diag(iDTV     )%desc = 'T vertical diffusion'
      diag(iDTV     )%unit = 'K/s'
      diag(iDTV     )%aname= 'TURBT   '
      diag(iDTV     )%adesc= 'Temperature Tendency due to Turbulence'
      diag(iDTV     )%aunit= 'K/day'
      diag(iDTV     )%convfac = 86400.0

      diag(iDUV     )%name = 'DUV     '
      diag(iDUV     )%desc = 'U tendency from vertical diffusion'
      diag(iDUV     )%unit = 'm/s2'
      diag(iDUV     )%aname= 'TURBU   '
      diag(iDUV     )%adesc= 'Zonal Wind Tendency due to Turbulence'
      diag(iDUV     )%aunit= 'm/s/day'
      diag(iDUV     )%convfac = 86400.

      diag(iDVV     )%name = 'DVV     '
      diag(iDVV     )%desc = 'V tendency from vertical diffusion'
      diag(iDVV     )%unit = 'm/s2'
      diag(iDVV     )%aname= 'TURBV   '
      diag(iDVV     )%adesc= 'Meridional Wind Tendency due to Turbulence'
      diag(iDVV     )%aunit= 'm/s/day'
      diag(iDVV     )%convfac = 86400.

      diag(iEFFCLD  )%name = 'EFFCLD  '
      diag(iEFFCLD  )%desc = 'Effective cloud fraction'
      diag(iEFFCLD  )%unit = 'fraction'

      diag(iEVAPL   )%name = 'EVAPL   '
      diag(iEVAPL   )%desc = 'Large-scale evaporation'
      diag(iEVAPL   )%unit = 'kg/kg/s'

      diag(iH       )%name = 'H       '
      diag(iH       )%desc = 'Geopotential height at mid-layer'
      diag(iH       )%unit = 'm'
      diag(iH       )%aname= 'HGHT    '

      diag(iHGHTE   )%name = 'HGHTE   '
      diag(iHGHTE   )%desc = 'Geopotential height at layer top'
      diag(iHGHTE   )%unit = 'm'

      diag(iKVH     )%name = 'KVH     '
      diag(iKVH     )%desc = 'Vertical diffusion diffusivities '  &
                           // '(heat/moisture)'
      diag(iKVH     )%unit = 'm2/s'
      diag(iKVH     )%aname= 'KH      '

      diag(iKVM     )%name = 'KVM     '
      diag(iKVM     )%desc = 'Eddy diffusivity for momentum'
      diag(iKVM     )%unit = 'm2/s'

      diag(iO3VMR   )%name = 'O3VMR   '
      diag(iO3VMR   )%desc = 'Ozone used in radiative transfer'
      diag(iO3VMR   )%unit = 'mol/mol'

      diag(iOMEGA   )%name = 'OMEGA   '
      diag(iOMEGA   )%desc = 'Vertical pressure velocity'
      diag(iOMEGA   )%unit = 'Pa/s'
      diag(iOMEGA   )%adesc= 'Vertical velocity'
      diag(iOMEGA   )%aunit= 'hPa/day'
      diag(iOMEGA   )%convfac = 864.0 

      diag(iPV      )%name = 'PV      '
      diag(iPV      )%desc = 'Ertels potential vorticity'
      diag(iPV      )%unit = 'm2/(kg*sec)'

      diag(iQ       )%name = 'Q       '
      diag(iQ       )%desc = 'Specific humidity'
      diag(iQ       )%unit = 'kg/kg'
      diag(iQ       )%aname= 'SPHU    '
      diag(iQ       )%aunit= 'g/kg'
      diag(iQ       )%convfac = 1000.0

      diag(iQC      )%name = 'QC      '
      diag(iQC      )%desc = 'Specific humidity of cloud condensate'
      diag(iQC      )%unit = 'kg/kg'

      diag(iQRL     )%name = 'QRL     '
      diag(iQRL     )%desc = 'Longwave heating rate'
      diag(iQRL     )%unit = 'K/s'
      diag(iQRL     )%aname= 'RADLW   '
      diag(iQRL     )%adesc= 'Temperature Tendency due to Longwave Radiation'
      diag(iQRL     )%aunit= 'K/day'
      diag(iQRL     )%convfac = 86400.0

      diag(iQRS     )%name = 'QRS     '
      diag(iQRS     )%desc = 'Shortwave heating rate'
      diag(iQRS     )%unit = 'K/s'
      diag(iQRS     )%aname= 'RADSW   '
      diag(iQRS     )%adesc= 'Temperature Tendency due to Shortwave Radiation'
      diag(iQRS     )%aunit= 'K/day'
      diag(iQRS     )%convfac = 86400.0

      diag(iRAYFDT  )%name = 'RAYFDT  '
      diag(iRAYFDT  )%desc = 'T Tendency - Rayleigh friction'
      diag(iRAYFDT  )%unit = 'K/s'
      diag(iRAYFDT  )%aname= 'RFT     '
      diag(iRAYFDT  )%adesc= 'Temperature Tendency due to Rayleigh Friction'
      diag(iRAYFDT  )%aunit= 'K/day'
      diag(iRAYFDT  )%convfac = 86400.0

      diag(iRELHUM  )%name = 'RELHUM  '
      diag(iRELHUM  )%desc = 'Relative Humidity after cloud physics'
      diag(iRELHUM  )%unit = '%'

      diag(iRHCLR   )%name = 'RHCLR   '
      diag(iRHCLR   )%desc = 'Relative Humidity in clear region'
      diag(iRHCLR   )%unit = '%'

      diag(iRNEVPDQ )%name = 'RNEVPDQ '
      diag(iRNEVPDQ )%desc = 'Q Tendency - Rain evaporation'
      diag(iRNEVPDQ )%unit = 'kg/kg/s'
      diag(iRNEVPDQ )%aname= 'DQLS    '
      diag(iRNEVPDQ )%adesc= 'Specific Humidity Tendency due to Stratiform Processes'
      diag(iRNEVPDQ )%aunit= 'g/kg/day'
      diag(iRNEVPDQ )%convfac = 86400000.0

      diag(iRNEVPDT )%name = 'RNEVPDT '
      diag(iRNEVPDT )%desc = 'T Tendency - Rain evaporation'
      diag(iRNEVPDT )%unit = 'K/s'         
      diag(iRNEVPDT )%aname= 'DTLS    '
      diag(iRNEVPDT )%adesc= 'Temperature Tendency due to Stratiform Processes'
      diag(iRNEVPDT )%aunit= 'K/day'         
      diag(iRNEVPDT )%convfac = 86400.       

      diag(iSETLWP  )%name = 'SETLWP  '
      diag(iSETLWP  )%desc = 'Specified liquid water path lengths'
      diag(iSETLWP  )%unit = 'gram/m2'

      diag(iSTRATCLD)%name = 'STRATCLD'
      diag(iSTRATCLD)%desc = 'Stratiform cloud amount'
      diag(iSTRATCLD)%unit = 'fraction'

      diag(iT       )%name = 'T       '
      diag(iT       )%desc = 'Temperature'
      diag(iT       )%unit = 'K'
      diag(iT       )%aname= 'TMPU    '

      diag(iTKE     )%name = 'TKE     '
      diag(iTKE     )%desc = 'Turbulent kinetic energy'
      diag(iTKE     )%unit = '(m/s)2'

      diag(iTTMGW   )%name = 'TTMGW   '
      diag(iTTMGW   )%desc = 'T tendency - gravity wave drag'
      diag(iTTMGW   )%unit = 'K/s'
      diag(iTTMGW   )%aname= 'GWDT    '
      diag(iTTMGW   )%adesc= 'Temperature Tendency due to Gravity Wave Drag'
      diag(iTTMGW   )%aunit= 'K/day'
      diag(iTTMGW   )%convfac = 86400.0

      diag(iU       )%name = 'U       '
      diag(iU       )%desc = 'U wind'
      diag(iU       )%unit = 'm/s'
      diag(iU       )%aname= 'UWND    '
      diag(iU       )%adesc= 'Zonal Wind'

      diag(iUQ      )%name = 'UQ      '
      diag(iUQ      )%desc = 'U wind * specific humidity'
      diag(iUQ      )%unit = 'm/s*kg/kg'

      diag(iUT      )%name = 'UT      '
      diag(iUT      )%desc = 'U wind * temperature'
      diag(iUT      )%unit = 'm/s*K'

      diag(iUTGW    )%name = 'UTGW    '
      diag(iUTGW    )%desc = 'U tendency - gravity wave drag'
      diag(iUTGW    )%unit = 'm/s2'
      diag(iUTGW    )%aname= 'GWDU    '
      diag(iUTGW    )%aunit= 'm/s/day'
      diag(iUTGW    )%adesc= 'Zonal Wind Tendency due to Gravity Wave Drag'
      diag(iUTGW    )%convfac = 86400. 

      diag(iUU      )%name = 'UU      '
      diag(iUU      )%desc = 'U wind * U wind'
      diag(iUU      )%unit = 'm/s*m/s'

      diag(iUV      )%name = 'UV      '
      diag(iUV      )%desc = 'U wind * V wind'
      diag(iUV      )%unit = 'm/s*m/s'

      diag(iV       )%name = 'V       '
      diag(iV       )%desc = 'V wind'
      diag(iV       )%unit = 'm/s'
      diag(iV       )%aname= 'VWND    '
      diag(iV       )%adesc= 'Meridional Wind'

      diag(iVD01    )%name = 'VD01    '
      diag(iVD01    )%desc = 'Vertical diffusion tendency of water vapor'
      diag(iVD01    )%unit = 'kg/kg/s'
      diag(iVD01    )%aname= 'TURBQ   '
      diag(iVD01    )%adesc= 'Specific Humidity Tendency due to Turbulence'
      diag(iVD01    )%aunit = 'g/kg/day'
      diag(iVD01    )%convfac = 86400000.0

      diag(iVQ      )%name = 'VQ      '
      diag(iVQ      )%desc = 'V wind * specific humidity'
      diag(iVQ      )%unit = 'm/s*kg/kg'

      diag(iVT      )%name = 'VT      '
      diag(iVT      )%desc = 'V wind * temperature'
      diag(iVT      )%unit = 'm/s*K'

      diag(iVTGW    )%name = 'VTGW    '
      diag(iVTGW    )%desc = 'V tendency - gravity wave drag'
      diag(iVTGW    )%unit = 'm/s2'
      diag(iVTGW    )%aname= 'GWDV    '
      diag(iVTGW    )%adesc= 'Meridional Wind Tendency due to Gravity Wave Drag'
      diag(iVTGW    )%aunit= 'm/s/day'
      diag(iVTGW    )%convfac = 86400. 

      diag(iVV      )%name = 'VV      '
      diag(iVV      )%desc = 'V wind * V wind'
      diag(iVV      )%unit = 'm/s*m/s'

      diag(iZMCME   )%name = 'ZMCME   '
      diag(iZMCME   )%desc = 'Condensation - evaporation from Z&M scheme'
      diag(iZMCME   )%unit = 'kg/kg/s'

      diag(iZMDLF   )%name = 'ZMDLF   '
      diag(iZMDLF   )%desc = 'Detrainment of cloud water from Z&M scheme'
      diag(iZMDLF   )%unit = 'kg/kg/s'

      diag(iZMDQ    )%name = 'ZMDQ    '
      diag(iZMDQ    )%desc = 'q tendency - Zhang-McFarlane convection'
      diag(iZMDQ    )%unit = 'kg/kg/s'

      diag(iZMDQR   )%name = 'ZMDQR   '
      diag(iZMDQR   )%desc = 'Rain production rate - Z&M convection'
      diag(iZMDQR   )%unit = 'kg/kg/s'

      diag(iZMDT    )%name = 'ZMDT    '
      diag(iZMDT    )%desc = 'T tendency - Zhang-McFarlane convection'
      diag(iZMDT    )%unit = 'K/s'

      diag(iZMDU    )%name = 'ZMDU    '
      diag(iZMDU    )%desc = 'Updraft detrainment mass flux - Z&M'  &
                                  //' convection'
      diag(iZMDU    )%unit = 'Pa/s'

      diag(iZMED    )%name = 'ZMED    '
      diag(iZMED    )%desc = 'Downdraft entrainment mass flux - Z&M'  &
                                  //' convection'
      diag(iZMED    )%unit = 'Pa/s'

      diag(iZMEPS   )%name = 'ZMEPS   '
      diag(iZMEPS   )%desc = 'Fractional entrainment - Z&M convection'
      diag(iZMEPS   )%unit = '1/s'

      diag(iZMEU    )%name = 'ZMEU    '
      diag(iZMEU    )%desc = 'Updraft entrainment mass flux - Z&M'  &
                                  //' convection'
      diag(iZMEU    )%unit = 'Pa/s'

      diag(iZMEVP   )%name = 'ZMEVP   '
      diag(iZMEVP   )%desc = 'downdraft evaporation from Z&M convection'
      diag(iZMEVP   )%unit = 'kg/kg/s'

      diag(iZMMD    )%name = 'ZMMD    '
      diag(iZMMD    )%desc = 'Downdraft mass flux - Z&M convection'
      diag(iZMMD    )%unit = 'Pa/s'

      diag(iZMMU    )%name = 'ZMMU    '
      diag(iZMMU    )%desc = 'Updraft mass flux - Z&M convection'
      diag(iZMMU    )%unit = 'pa/s'

      diag(iZMPFLX  )%name = 'ZMPFLX  '
      diag(iZMPFLX  )%desc = 'Precipitation flux - Z&M convection'
      diag(iZMPFLX  )%unit = 'kg/m2/s'

      diag(iZMQL    )%name = 'ZMQL    '
      diag(iZMQL    )%desc = 'Cloud water in updraft - Z&M convection'
      diag(iZMQL    )%unit = 'kg/kg'

#ifdef FVCHEM
      diag(iBR      )%name = 'BR      '
      diag(iBR      )%desc = 'Atomic bromine'
      diag(iBR      )%unit = 'mol/mol'

      diag(iBRCL    )%name = 'BRCL    '
      diag(iBRCL    )%desc = 'Bromine chloride'
      diag(iBRCL    )%unit = 'mol/mol'

      diag(iBRO     )%name = 'BRO     '
      diag(iBRO     )%desc = 'Bromine monoxide'
      diag(iBRO     )%unit = 'mol/mol'

      diag(iBRONO2  )%name = 'BRONO2  '
      diag(iBRONO2  )%desc = 'Bromine nitrate'
      diag(iBRONO2  )%unit = 'mol/mol'

      diag(iBRX     )%name = 'BRX     '
      diag(iBRX     )%desc = 'Odd bromine'
      diag(iBRX     )%unit = 'mol/mol'

      diag(iCCL4    )%name = 'CCL4    '
      diag(iCCL4    )%desc = 'Carbon tetrachloride'
      diag(iCCL4    )%unit = 'mol/mol'

      diag(iCH2O    )%name = 'CH2O    '
      diag(iCH2O    )%desc = 'Formaldehyde'
      diag(iCH2O    )%unit = 'mol/mol'

      diag(iCH3BR   )%name = 'CH3BR   '
      diag(iCH3BR   )%desc = 'Methyl bromide'
      diag(iCH3BR   )%unit = 'mol/mol'

      diag(iCH3CCL3 )%name = 'CH3CCL3 '
      diag(iCH3CCL3 )%desc = 'Methyl chloroform'
      diag(iCH3CCL3 )%unit = 'mol/mol'

      diag(iCH3CL   )%name = 'CH3CL   '
      diag(iCH3CL   )%desc = 'Methyl chloride'
      diag(iCH3CL   )%unit = 'mol/mol'
      
      diag(iCH3O2   )%name = 'CH3O2   '
      diag(iCH3O2   )%desc = 'Methyl peroxide'
      diag(iCH3O2   )%unit = 'mol/mol'
      
      diag(iCH3OOH  )%name = 'CH3OOH  '
      diag(iCH3OOH  )%desc = 'Methyl hydroperoxide'
      diag(iCH3OOH  )%unit = 'mol/mol'

      diag(iCH4     )%name = 'CH4     '
      diag(iCH4     )%desc = 'Methane'
      diag(iCH4     )%unit = 'mol/mol'

      diag(iCL      )%name = 'CL      '
      diag(iCL      )%desc = 'Atomic chlorine'
      diag(iCL      )%unit = 'mol/mol'

      diag(iCL2     )%name = 'CL2     '
      diag(iCL2     )%desc = 'Molecular chlorine'
      diag(iCL2     )%unit = 'mol/mol'

      diag(iCL2O2   )%name = 'CL2O2   '
      diag(iCL2O2   )%desc = 'Dichlorine peroxide'
      diag(iCL2O2   )%unit = 'mol/mol'

      diag(iCLO     )%name = 'CLO     '
      diag(iCLO     )%desc = 'Chlorine monoxide'
      diag(iCLO     )%unit = 'mol/mol'

      diag(iCLONO2  )%name = 'CLONO2  '
      diag(iCLONO2  )%desc = 'Chlorine nitrate'
      diag(iCLONO2  )%unit = 'mol/mol'

      diag(iCLX     )%name = 'CLX     '
      diag(iCLX     )%desc = 'Odd chlorine'
      diag(iCLX     )%unit = 'mol/mol'

! --------------------------------------------------------------------
!               8 CO regions and types for INTEX-B 2006
! --------------------------------------------------------------------

      diag(iCO      )%name = 'CO      '
      diag(iCO      )%desc = 'Global Carbon Monoxide'
      diag(iCO      )%unit = 'mol/mol'

      diag(iCONOAMAN)%name = 'CONOAMAN'
      diag(iCONOAMAN)%desc = 'North American anthropogenic CO'
      diag(iCONOAMAN)%unit = 'mol/mol'

      diag(iCOCEAMAN)%name = 'COCEAMAN'
      diag(iCOCEAMAN)%desc = 'Central American anthropogenic CO'
      diag(iCOCEAMAN)%unit = 'mol/mol'

      diag(iCOWHBB  )%name = 'COWHBB  '
      diag(iCOWHBB  )%desc = 'Western Hemisphere biomass burning CO'
      diag(iCOWHBB  )%unit = 'mol/mol'

      diag(iCOASIAAN)%name = 'COASIAAN'
      diag(iCOASIAAN)%desc = 'Asian anthropogenic CO'
      diag(iCOASIAAN)%unit = 'mol/mol'

      diag(iCOASNBB )%name = 'COASNBB '
      diag(iCOASNBB )%desc = 'Northern Asia biomass burning CO'
      diag(iCOASNBB )%unit = 'mol/mol'

      diag(iCOASSBB )%name = 'COASSBB '
      diag(iCOASSBB )%desc = 'Southern Asia biomass burning CO'
      diag(iCOASSBB )%unit = 'mol/mol'

      diag(iCOFDAN  )%name = 'COFDAN  '
      diag(iCOFDAN  )%desc = 'Mexico City anthropogenic CO'
      diag(iCOFDAN  )%unit = 'mol/mol'

! --------------------------------------------------------------------

      diag(iCOFF    )%name = 'COFF    '
      diag(iCOFF    )%desc = 'CO tagged to Fossil Fuel'
      diag(iCOFF    )%unit = 'mol/mol'

      diag(iCOBF    )%name = 'COBF    '
      diag(iCOBF    )%desc = 'CO tagged to Biofuel'
      diag(iCOBF    )%unit = 'mol/mol'

      diag(iCOBB    )%name = 'COBB    '
      diag(iCOBB    )%desc = 'CO tagged to Biomass Burning'
      diag(iCOBB    )%unit = 'mol/mol'

      diag(iCOBI    )%name = 'COBI    '
      diag(iCOBI    )%desc = 'CO tagged to Biogenic'
      diag(iCOBI    )%unit = 'mol/mol'

      diag(iCONAMERI)%name = 'COnam   '
      diag(iCONAMERI)%desc = 'North American CO'
      diag(iCONAMERI)%unit = 'mol/mol'

      diag(iCOSAMERI)%name = 'COsam   '
      diag(iCOSAMERI)%desc = 'South American CO'
      diag(iCOSAMERI)%unit = 'mol/mol'

      diag(iCOAFRICA)%name = 'COafr   '
      diag(iCOAFRICA)%desc = 'African CO'
      diag(iCOAFRICA)%unit = 'mol/mol'

      diag(iCO2     )%name = 'CO2     '
      diag(iCO2     )%desc = 'Carbon Dioxide'
      diag(iCO2     )%unit = 'mol/mol'

      diag(iCO2NAMER)%name = 'CO2nam  '
      diag(iCO2NAMER)%desc = 'North American CO2'
      diag(iCO2NAMER)%unit = 'mol/mol'

      diag(iCO2SAMER)%name = 'CO2sam  '
      diag(iCO2SAMER)%desc = 'South American CO2'
      diag(iCO2SAMER)%unit = 'mol/mol'

      diag(iCO2AFRIC)%name = 'CO2afr  '
      diag(iCO2AFRIC)%desc = 'African CO2'
      diag(iCO2AFRIC)%unit = 'mol/mol'

      diag(iF11     )%name = 'F11     '
      diag(iF11     )%desc = 'CFC-11 (CCl3F)'
      diag(iF11     )%unit = 'mol/mol'

      diag(iF113    )%name = 'F113    '
      diag(iF113    )%desc = 'CFC-113 (CCl2FCClF2)'
      diag(iF113    )%unit = 'mol/mol'

      diag(iF12     )%name = 'F12     '
      diag(iF12     )%desc = 'CFC-12 (CCl2F2)'
      diag(iF12     )%unit = 'mol/mol'

      diag(iH12_24  )%name = 'H12_24  '
      diag(iH12_24  )%desc = 'Halon 12_24'
      diag(iH12_24  )%unit = 'mol/mol'

      diag(iH1301   )%name = 'H1301   '
      diag(iH1301   )%desc = 'Halon 1301 (CBrF3)'
      diag(iH1301   )%unit = 'mol/mol'

      diag(iH2O2    )%name = 'H2O2    '
      diag(iH2O2    )%desc = 'Hydrogen peroxide'
      diag(iH2O2    )%unit = 'mol/mol'

      diag(iH2OCOND )%name = 'H2OCOND '
      diag(iH2OCOND )%desc = 'Condensed water vapor in chemistry'
      diag(iH2OCOND )%unit = 'mol/mol'

      diag(iHATOMIC )%name = 'HATOMIC '
      diag(iHATOMIC )%desc = 'Atomic hydrogen'
      diag(iHATOMIC )%unit = 'mol/mol'

      diag(iHBR     )%name = 'HBR     '
      diag(iHBR     )%desc = 'Hydrogen bromide'
      diag(iHBR     )%unit = 'mol/mol'

      diag(iHCFC    )%name = 'HCFC    '
      diag(iHCFC    )%desc = 'HCFC'
      diag(iHCFC    )%unit = 'mol/mol'

      diag(iHCL     )%name = 'HCL     '
      diag(iHCL     )%desc = 'Hydrochloric acid'
      diag(iHCL     )%unit = 'mol/mol'

      diag(iHNO3    )%name = 'HNO3    '
      diag(iHNO3    )%desc = 'Nitric acid'
      diag(iHNO3    )%unit = 'mol/mol'

      diag(iHNO3COND)%name = 'HNO3COND'
      diag(iHNO3COND)%desc = 'Condensed nitric acid'
      diag(iHNO3COND)%unit = 'mol/mol'

      diag(iHO2     )%name = 'HO2     '
      diag(iHO2     )%desc = 'Hydroperoxyl radical'
      diag(iHO2     )%unit = 'mol/mol'

      diag(iHO2NO2  )%name = 'HO2NO2  '
      diag(iHO2NO2  )%desc = 'Peroxynitric acid'
      diag(iHO2NO2  )%unit = 'mol/mol'

      diag(iHOBR    )%name = 'HOBR    '
      diag(iHOBR    )%desc = 'Hypobromous acid'
      diag(iHOBR    )%unit = 'mol/mol'

      diag(iHOCL    )%name = 'HOCL    '
      diag(iHOCL    )%desc = 'Hypochlorous acid'
      diag(iHOCL    )%unit = 'mol/mol'

      diag(iN       )%name = 'N       '
      diag(iN       )%desc = 'Atomic nitrogen'
      diag(iN       )%unit = 'mol/mol'

      diag(iN2O     )%name = 'N2O     '
      diag(iN2O     )%desc = 'Nitrous oxide'
      diag(iN2O     )%unit = 'mol/mol'

      diag(iN2O5    )%name = 'N2O5    '
      diag(iN2O5    )%desc = 'Dinitrogen pentoxide'
      diag(iN2O5    )%unit = 'mol/mol'

      diag(iNO      )%name = 'NO      '
      diag(iNO      )%desc = 'Nitric oxide'
      diag(iNO      )%unit = 'mol/mol'

      diag(iNO2     )%name = 'NO2     '
      diag(iNO2     )%desc = 'Nitrogen dioxide'
      diag(iNO2     )%unit = 'mol/mol'

      diag(iNO3     )%name = 'NO3     '
      diag(iNO3     )%desc = 'Nitrogen trioxide'
      diag(iNO3     )%unit = 'mol/mol'

      diag(iNOX     )%name = 'NOX     '
      diag(iNOX     )%desc = 'Odd nitrogen'
      diag(iNOX     )%unit = 'mol/mol'

      diag(iO1D     )%name = 'O1D     '
      diag(iO1D     )%desc = 'Atomic oxygen in the first excited state'
      diag(iO1D     )%unit = 'mol/mol'

      diag(iO3CHEM  )%name = 'O3CHEM  '
      diag(iO3CHEM  )%desc = 'Ozone volume mixing ratio from chemistry'
      diag(iO3CHEM  )%unit = 'mol/mol'

      diag(iO3P     )%name = 'O3P     '
      diag(iO3P     )%desc = 'Atomic oxygen in the ground state'
      diag(iO3P     )%unit = 'mol/mol'

      diag(iO3PARAM )%name = 'O3PARAM '
      diag(iO3PARAM )%desc = 'Parameterized ozone'
      diag(iO3PARAM )%unit = 'mol/mol'

      diag(iOCLO    )%name = 'OCLO    '
      diag(iOCLO    )%desc = 'Chlorine dioxide'
      diag(iOCLO    )%unit = 'mol/mol'

      diag(iOH      )%name = 'OH      '
      diag(iOH      )%desc = 'Hydroxyl radical'
      diag(iOH      )%unit = 'mol/mol'

      diag(iOX      )%name = 'OX      '
      diag(iOX      )%desc = 'Odd oxygen from Parameterized Chemistry'
      diag(iOX      )%unit = 'mol/mol'

      diag(iOXSTRAT )%name = 'OXSTRAT '
      diag(iOXSTRAT )%desc = 'Odd oxygen from Stratospheric Chemistry'
      diag(iOXSTRAT )%unit = 'mol/mol'

      diag(iOXTROP  )%name = 'OXTROP  '
      diag(iOXTROP  )%desc = 'Tropospheric odd oxygen'
      diag(iOXTROP  )%unit = 'mol/mol'

      diag(iDUMASS )%name = 'DUMASS '
      diag(iDUMASS )%desc = 'Dust 3D Mass Mixing Ratio'
      diag(iDUMASS )%unit = 'kg/kg'

      diag(iDUMASS25 )%name = 'DUMASS25 '
      diag(iDUMASS25 )%desc = 'Dust 3D Mass Mixing Ratio (PM2.5)'
      diag(iDUMASS25 )%unit = 'kg/kg'

      diag(iDUMASS1 )%name = 'DUMASS1 '
      diag(iDUMASS1 )%desc = 'Dust sub-micron 3D Mass Mixing Ratio'
      diag(iDUMASS1 )%unit = 'kg/kg'

      diag(iSSMASS )%name = 'SSMASS '
      diag(iSSMASS )%desc = 'Seasalt 3D Mass Mixing Ratio'
      diag(iSSMASS )%unit = 'kg/kg'

      diag(iSSMASS25 )%name = 'SSMASS25 '
      diag(iSSMASS25 )%desc = 'Seasalt 3D Mass Mixing Ratio (PM2.5)'
      diag(iSSMASS25 )%unit = 'kg/kg'

      diag(iBCMASS )%name = 'BCMASS '
      diag(iBCMASS )%desc = 'Black Carbon 3D Mass Mixing Ratio'
      diag(iBCMASS )%unit = 'kg/kg'

      diag(iOCMASS )%name = 'OCMASS '
      diag(iOCMASS )%desc = 'Organic Carbon 3D Mass Mixing Ratio'
      diag(iOCMASS )%unit = 'kg/kg'

      diag(iSO4MASS )%name = 'SO4MASS '
      diag(iSO4MASS )%desc = 'SO4 Aerosol 3D Mass Mixing Ratio'
      diag(iSO4MASS )%unit = 'kg/kg'

      diag(iPSO2    )%name = 'PSO2    '
      diag(iPSO2    )%desc = 'Chemical production of SO2 from DMS ox'
      diag(iPSO2    )%unit = 'kg/m2/s'

      diag(iPMSA    )%name = 'PMSA    '
      diag(iPMSA    )%desc = 'Chemical production of MSA from DMS ox'
      diag(iPMSA    )%unit = 'kg/m2/s'

      diag(iPSO4g    )%name = 'PSO4g    '
      diag(iPSO4g    )%desc = 'Chemical production of SO4 from SO2 ox (gas)'
      diag(iPSO4g    )%unit = 'kg/m2/s'

      diag(iPSO4aq  )%name = 'PSO4aq    '
      diag(iPSO4aq  )%desc = 'Chemical production of SO4 from SO2 ox (aqueous)'
      diag(iPSO4aq  )%unit = 'kg/m2/s'

      diag(iPSO4wet  )%name = 'PSO4wet    '
      diag(iPSO4wet  )%desc = 'Chemical production of SO4 from SO2 ox (wet dep)'
      diag(iPSO4wet  )%unit = 'kg/m2/s'

      diag(iQ4AGE   )%name = 'Q4AGE   '
      diag(iQ4AGE   )%desc = 'Surface souce gas for computation of age-of-air'
      diag(iQ4AGE   )%unit = 'kg/kg'

#endif

      RETURN
      END SUBROUTINE diaglist
