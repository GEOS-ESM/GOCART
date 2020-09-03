      MODULE mod_diag

! !REVISION HISTORY:
!
!  16Jul2003  Norris  Added a few comments. Will now allow a
!                     -ve history freq (-1) for non-output diags.
!  18Apr2005 Nielsen  Removed iTRACER1-iTRACER5 and iCTRCR1-iCTRCR5, and
!                     added i-values for SC_GridComp species.

      implicit none
!
!    -----------------------------------------------------------
!     Refer to ../README file for how to add or delete a field.
!    -----------------------------------------------------------
!
! ... 2-D diagnostic fields
!
      integer  iALBEDO
      integer  iALDIF
      integer  iALDIR
      integer  iASDIF
      integer  iASDIR           
      integer  iBMA           
      integer  iBULKTS        
      integer  iCAPEMX
      integer  iCLDHGH
      integer  iCLDLOW
      integer  iCLDMED
      integer  iCLDPRS
      integer  iCLDTMP
      integer  iCLDTOT
      integer  iCNVCLD
      integer  iDRUNOFF
      integer  iDTG
      integer  iEMSFC
      integer  iFLNS
      integer  iFLNSC
      integer  iFLNT
      integer  iFLNTC
      integer  iFRACLAKE
      integer  iFRACVEG
      integer  iFSDS
      integer  iFSNS
      integer  iFSNSC
      integer  iFSNT
      integer  iFSNTC
      integer  iGWETROOT
      integer  iGWETTC
      integer  iGWETTOP
      integer  iH300 
      integer  iH500 
      integer  iHSURF
      integer  iHTLCL
      integer  iHTMMSE
      integer  iLAI
      integer  iLHFX
      integer  iLWSH
      integer  iO3DU
      integer  iORO
      integer  iOSR
      integer  iOSRCLR
      integer  iPARDF
      integer  iPARDR
      integer  iPBLH
      integer  iPREACC
      integer  iPRECC
      integer  iPRECL
      integer  iPRECL_RH
      integer  iQ10M
      integer  iQ2M
      integer  iQFLX
      integer  iQPERT
      integer  iSHFX
      integer  iSLP
      integer  iSNOWDP
      integer  iSNOWH
      integer  iSOILWC1
      integer  iSOILWC2
      integer  iSOLIN
      integer  iSRFRAD
      integer  iSRUNOFF
      integer  iSURFP
      integer  iT10M
      integer  iT200 
      integer  iT2M
      integer  iT850 
      integer  iTAUGWX
      integer  iTAUGWY
      integer  iTAUX
      integer  iTAUY
      integer  iTHICK
      integer  iTLAKE1
      integer  iTPERT
      integer  iTQ
      integer  iTRAD
      integer  iTROPP
      integer  iTROPQ
      integer  iTROPT
      integer  iTSKIN
      integer  iTSLAKE
      integer  iTSOIL1
      integer  iTVEG
      integer  iU10M
      integer  iU200 
      integer  iU2M
      integer  iU500
      integer  iU850 
      integer  iUSTAR
      integer  iV10M
      integer  iV200 
      integer  iV2M
      integer  iV500
      integer  iV850 
      integer  iVAVET 
      integer  iVAVEU 
      integer  iVAVEUQ
      integer  iVAVEUT
      integer  iVAVEV 
      integer  iVAVEVQ
      integer  iVAVEVT
      integer  iWATCAN
      integer  iWBALLAKE
      integer  iZ0H 
      integer  iZ0M 
      integer  iZMMB 
      integer  iZMPR 
      integer  iZPD 
!
! ... 3-D diagnostic fields
!

      integer  iAIRDENS
      integer  iCAPE
      integer  iCGS
      integer  iCLDLWP
      integer  iCLOUD
      integer  iCLOUDUP
      integer  iCMFDQ
      integer  iCMFDQR2
      integer  iCMFDT
      integer  iCMFDTR
      integer  iCMFETR
      integer  iCMFMC
      integer  iCMFMC2
      integer  iCONVCLD
      integer  iDCAFDT
      integer  iDIABDT
      integer  iDQCOND
      integer  iDQPBLCG
      integer  iDQRL
      integer  iDTCOND
      integer  iDTPBLCG
      integer  iDTRAIN 
      integer  iDTV
      integer  iDUV
      integer  iDVV
      integer  iEFFCLD
      integer  iEVAPL
      integer  iH
      integer  iHGHTE
      integer  iHKBETA
      integer  iHKETA
      integer  iKVH
      integer  iKVM
      integer  iO3VMR
      integer  iOMEGA
      integer  iPV     
      integer  iQ
      integer  iQC
      integer  iQRL
      integer  iQRS
      integer  iRAYFDT
      integer  iRELHUM
      integer  iRHCLR
      integer  iRNEVPDQ
      integer  iRNEVPDT
      integer  iSETLWP
      integer  iSTRATCLD
      integer  iT
      integer  iTAUCLI
      integer  iTAUCLW
      integer  iTKE
      integer  iTTMGW
      integer  iU
      integer  iUQ
      integer  iUT
      integer  iUTGW
      integer  iUU
      integer  iUV
      integer  iV
      integer  iVD01
      integer  iVQ
      integer  iVT
      integer  iVTGW
      integer  iVV
      integer  iZMCME
      integer  iZMDLF
      integer  iZMDQ
      integer  iZMDQR
      integer  iZMDT
      integer  iZMDU
      integer  iZMED
      integer  iZMEPS
      integer  iZMEU
      integer  iZMEVP
      integer  iZMMD
      integer  iZMMU
      integer  iZMPFLX
      integer  iZMQL

      integer  pd2d     ! Total number of 2-D diagnostic fields
      integer  pd3d     ! Total number of 3-D diagnostic fields
      integer  pdiag    ! Total number of diagnostic fields
!
! ... 2-D diagnosis fields
!
      parameter (iALBEDO   =            1)
      parameter (iALDIF    = iALBEDO  + 1)
      parameter (iALDIR    = iALDIF   + 1)
      parameter (iASDIF    = iALDIR   + 1)
      parameter (iASDIR    = iASDIF   + 1)
      parameter (iBMA      = iASDIR   + 1)
      parameter (iBULKTS   = iBMA     + 1)
      parameter (iCAPEMX   = iBULKTS  + 1)
      parameter (iCLDHGH   = iCAPEMX  + 1)
      parameter (iCLDLOW   = iCLDHGH  + 1)
      parameter (iCLDMED   = iCLDLOW  + 1)
      parameter (iCLDPRS   = iCLDMED  + 1)
      parameter (iCLDTMP   = iCLDPRS  + 1)
      parameter (iCLDTOT   = iCLDTMP  + 1)
      parameter (iCNVCLD   = iCLDTOT  + 1)
      parameter (iDRUNOFF  = iCNVCLD  + 1)
      parameter (iDTG      = iDRUNOFF + 1)
      parameter (iEMSFC    = iDTG     + 1)
      parameter (iFLNS     = iEMSFC   + 1)
      parameter (iFRACLAKE = iFLNS    + 1)
      parameter (iFRACVEG  = iFRACLAKE+ 1)
      parameter (iFLNSC    = iFRACVEG + 1)
      parameter (iFLNT     = iFLNSC   + 1)
      parameter (iFLNTC    = iFLNT    + 1)
      parameter (iFSDS     = iFLNTC   + 1)
      parameter (iFSNS     = iFSDS    + 1)
      parameter (iFSNSC    = iFSNS    + 1)
      parameter (iFSNT     = iFSNSC   + 1)
      parameter (iFSNTC    = iFSNT    + 1)
      parameter (iGWETROOT = iFSNTC   + 1)
      parameter (iGWETTC   = iGWETROOT+ 1)
      parameter (iGWETTOP  = iGWETTC  + 1)
      parameter (iH300     = iGWETTOP + 1)
      parameter (iH500     = iH300    + 1)
      parameter (iHSURF    = iH500    + 1)
      parameter (iHTLCL    = iHSURF   + 1)
      parameter (iHTMMSE   = iHTLCL   + 1)
      parameter (iLAI      = iHTMMSE  + 1)
      parameter (iLHFX     = iLAI     + 1)
      parameter (iLWSH     = iLHFX    + 1)
      parameter (iO3DU     = iLWSH    + 1)
      parameter (iORO      = iO3DU    + 1)
      parameter (iOSR      = iORO     + 1)
      parameter (iOSRCLR   = iOSR     + 1)
      parameter (iPARDF    = iOSRCLR  + 1)
      parameter (iPARDR    = iPARDF   + 1)
      parameter (iPBLH     = iPARDR   + 1)
      parameter (iPREACC   = iPBLH    + 1)
      parameter (iPRECC    = iPREACC  + 1)
      parameter (iPRECL    = iPRECC   + 1)
      parameter (iPRECL_RH = iPRECL   + 1)
      parameter (iQ10M     = iPRECL_RH+ 1)
      parameter (iQ2M      = iQ10M    + 1)
      parameter (iQFLX     = iQ2M     + 1)
      parameter (iQPERT    = iQFLX    + 1)
      parameter (iSHFX     = iQPERT   + 1)
      parameter (iSLP      = iSHFX    + 1)
      parameter (iSNOWDP   = iSLP     + 1)
      parameter (iSNOWH    = iSNOWDP  + 1)
      parameter (iSOILWC1  = iSNOWH   + 1)
      parameter (iSOILWC2  = iSOILWC1 + 1)
      parameter (iSOLIN    = iSOILWC2 + 1)
      parameter (iSRFRAD   = iSOLIN   + 1)
      parameter (iSRUNOFF  = iSRFRAD  + 1)
      parameter (iSURFP    = iSRUNOFF + 1)
      parameter (iT10M     = iSURFP   + 1)
      parameter (iT200     = iT10M    + 1)
      parameter (iT2M      = iT200    + 1)
      parameter (iT850     = iT2M     + 1)
      parameter (iTAUGWX   = iT850    + 1)
      parameter (iTAUGWY   = iTAUGWX  + 1)
      parameter (iTAUX     = iTAUGWY  + 1)
      parameter (iTAUY     = iTAUX    + 1)
      parameter (iTHICK    = iTAUY    + 1)
      parameter (iTLAKE1   = iTHICK   + 1)
      parameter (iTPERT    = iTLAKE1  + 1)
      parameter (iTQ       = iTPERT   + 1)
      parameter (iTRAD     = iTQ      + 1)
      parameter (iTROPP    = iTRAD    + 1)
      parameter (iTROPQ    = iTROPP   + 1)
      parameter (iTROPT    = iTROPQ   + 1)
      parameter (iTSKIN    = iTROPT   + 1)
      parameter (iTSLAKE   = iTSKIN   + 1)
      parameter (iTSOIL1   = iTSLAKE  + 1)
      parameter (iTVEG     = iTSOIL1  + 1)
      parameter (iU10M     = iTVEG    + 1)
      parameter (iU200     = iU10M    + 1)
      parameter (iU2M      = iU200    + 1)
      parameter (iU500     = iU2M     + 1)
      parameter (iU850     = iU500    + 1)
      parameter (iUSTAR    = iU850    + 1)
      parameter (iV10M     = iUSTAR   + 1)
      parameter (iV200     = iV10M    + 1)
      parameter (iV2M      = iV200    + 1)
      parameter (iV500     = iV2M     + 1)
      parameter (iV850     = iV500    + 1)
      parameter (iVAVET    = iV850    + 1)
      parameter (iVAVEU    = iVAVET   + 1)
      parameter (iVAVEUQ   = iVAVEU   + 1)
      parameter (iVAVEUT   = iVAVEUQ  + 1)
      parameter (iVAVEV    = iVAVEUT  + 1)
      parameter (iVAVEVQ   = iVAVEV   + 1)
      parameter (iVAVEVT   = iVAVEVQ  + 1)
      parameter (iWATCAN   = iVAVEVT  + 1)
      parameter (iWBALLAKE = iWATCAN  + 1)
      parameter (iZ0H      = iWBALLAKE+ 1)
      parameter (iZ0M      = iZ0H     + 1)
      parameter (iZMMB     = iZ0M     + 1)
      parameter (iZMPR     = iZMMB    + 1)
      parameter (iZPD      = iZMPR    + 1)

#ifdef FVCHEM
      integer, parameter :: iDUEM001 = iZPD + 1
      integer, parameter :: iDUEM002 = iDUEM001 + 1
      integer, parameter :: iDUEM003 = iDUEM002 + 1
      integer, parameter :: iDUEM004 = iDUEM003 + 1
      integer, parameter :: iDUEM005 = iDUEM004 + 1
      integer, parameter :: iDUEM006 = iDUEM005 + 1
      integer, parameter :: iDUEM007 = iDUEM006 + 1
      integer, parameter :: iDUEM008 = iDUEM007 + 1

      integer, parameter :: iDUSD001 = iDUEM008 + 1
      integer, parameter :: iDUSD002 = iDUSD001 + 1
      integer, parameter :: iDUSD003 = iDUSD002 + 1
      integer, parameter :: iDUSD004 = iDUSD003 + 1
      integer, parameter :: iDUSD005 = iDUSD004 + 1
      integer, parameter :: iDUSD006 = iDUSD005 + 1
      integer, parameter :: iDUSD007 = iDUSD006 + 1
      integer, parameter :: iDUSD008 = iDUSD007 + 1

      integer, parameter :: iDUDP001 = iDUSD008 + 1
      integer, parameter :: iDUDP002 = iDUDP001 + 1
      integer, parameter :: iDUDP003 = iDUDP002 + 1
      integer, parameter :: iDUDP004 = iDUDP003 + 1
      integer, parameter :: iDUDP005 = iDUDP004 + 1
      integer, parameter :: iDUDP006 = iDUDP005 + 1
      integer, parameter :: iDUDP007 = iDUDP006 + 1
      integer, parameter :: iDUDP008 = iDUDP007 + 1

      integer, parameter :: iDUWT001 = iDUDP008 + 1
      integer, parameter :: iDUWT002 = iDUWT001 + 1
      integer, parameter :: iDUWT003 = iDUWT002 + 1
      integer, parameter :: iDUWT004 = iDUWT003 + 1
      integer, parameter :: iDUWT005 = iDUWT004 + 1
      integer, parameter :: iDUWT006 = iDUWT005 + 1
      integer, parameter :: iDUWT007 = iDUWT006 + 1
      integer, parameter :: iDUWT008 = iDUWT007 + 1

      integer, parameter :: iDUSV001 = iDUWT008 + 1
      integer, parameter :: iDUSV002 = iDUSV001 + 1
      integer, parameter :: iDUSV003 = iDUSV002 + 1
      integer, parameter :: iDUSV004 = iDUSV003 + 1
      integer, parameter :: iDUSV005 = iDUSV004 + 1
      integer, parameter :: iDUSV006 = iDUSV005 + 1
      integer, parameter :: iDUSV007 = iDUSV006 + 1
      integer, parameter :: iDUSV008 = iDUSV007 + 1

      integer, parameter :: iDUSMASS = iDUSV008 + 1
      integer, parameter :: iDUCMASS = iDUSMASS + 1
      integer, parameter :: iDUSMASS1 = iDUCMASS + 1
      integer, parameter :: iDUCMASS1 = iDUSMASS1 + 1
      integer, parameter :: iDUEXTTAU = iDUCMASS1 + 1
      integer, parameter :: iDUSCATAU = iDUEXTTAU + 1
      integer, parameter :: iDUAERIDX = iDUSCATAU + 1

!     Dust submicron diagnostics (column integral)
      integer, parameter :: iDUSM25 = iDUAERIDX + 1
      integer, parameter :: iDUCM25 = iDUSM25 + 1
      integer, parameter :: iDUEXTT25 = iDUCM25 + 1
      integer, parameter :: iDUSCAT25 = iDUEXTT25 + 1

      integer, parameter :: iSSEM001 = iDUSCAT25 + 1
      integer, parameter :: iSSEM002 = iSSEM001 + 1
      integer, parameter :: iSSEM003 = iSSEM002 + 1
      integer, parameter :: iSSEM004 = iSSEM003 + 1
      integer, parameter :: iSSEM005 = iSSEM004 + 1
      integer, parameter :: iSSEM006 = iSSEM005 + 1
      integer, parameter :: iSSEM007 = iSSEM006 + 1
      integer, parameter :: iSSEM008 = iSSEM007 + 1

      integer, parameter :: iSSSD001 = iSSEM008 + 1
      integer, parameter :: iSSSD002 = iSSSD001 + 1
      integer, parameter :: iSSSD003 = iSSSD002 + 1
      integer, parameter :: iSSSD004 = iSSSD003 + 1
      integer, parameter :: iSSSD005 = iSSSD004 + 1
      integer, parameter :: iSSSD006 = iSSSD005 + 1
      integer, parameter :: iSSSD007 = iSSSD006 + 1
      integer, parameter :: iSSSD008 = iSSSD007 + 1

      integer, parameter :: iSSDP001 = iSSSD008 + 1
      integer, parameter :: iSSDP002 = iSSDP001 + 1
      integer, parameter :: iSSDP003 = iSSDP002 + 1
      integer, parameter :: iSSDP004 = iSSDP003 + 1
      integer, parameter :: iSSDP005 = iSSDP004 + 1
      integer, parameter :: iSSDP006 = iSSDP005 + 1
      integer, parameter :: iSSDP007 = iSSDP006 + 1
      integer, parameter :: iSSDP008 = iSSDP007 + 1

      integer, parameter :: iSSWT001 = iSSDP008 + 1
      integer, parameter :: iSSWT002 = iSSWT001 + 1
      integer, parameter :: iSSWT003 = iSSWT002 + 1
      integer, parameter :: iSSWT004 = iSSWT003 + 1
      integer, parameter :: iSSWT005 = iSSWT004 + 1
      integer, parameter :: iSSWT006 = iSSWT005 + 1
      integer, parameter :: iSSWT007 = iSSWT006 + 1
      integer, parameter :: iSSWT008 = iSSWT007 + 1

      integer, parameter :: iSSSV001 = iSSWT008 + 1
      integer, parameter :: iSSSV002 = iSSSV001 + 1
      integer, parameter :: iSSSV003 = iSSSV002 + 1
      integer, parameter :: iSSSV004 = iSSSV003 + 1
      integer, parameter :: iSSSV005 = iSSSV004 + 1
      integer, parameter :: iSSSV006 = iSSSV005 + 1
      integer, parameter :: iSSSV007 = iSSSV006 + 1
      integer, parameter :: iSSSV008 = iSSSV007 + 1

      integer, parameter :: iSSSMASS = iSSSV008 + 1
      integer, parameter :: iSSCMASS = iSSSMASS + 1
      integer, parameter :: iSSEXTTAU = iSSCMASS + 1
      integer, parameter :: iSSSCATAU = iSSEXTTAU + 1

!     Seasalt submicron diagnostics (column integral)
      integer, parameter :: iSSSM25 = iSSSCATAU + 1
      integer, parameter :: iSSCM25 = iSSSM25 + 1
      integer, parameter :: iSSEXTT25 = iSSCM25 + 1
      integer, parameter :: iSSSCAT25 = iSSEXTT25 + 1

      integer, parameter :: iBCEM001 = iSSSCAT25 + 1
      integer, parameter :: iBCEM002 = iBCEM001 + 1
      integer, parameter :: iBCEM003 = iBCEM002 + 1
      integer, parameter :: iBCEM004 = iBCEM003 + 1
      integer, parameter :: iBCEM005 = iBCEM004 + 1
      integer, parameter :: iBCEM006 = iBCEM005 + 1
      integer, parameter :: iBCEM007 = iBCEM006 + 1
      integer, parameter :: iBCEM008 = iBCEM007 + 1

      integer, parameter :: iBCDP001 = iBCEM008 + 1
      integer, parameter :: iBCDP002 = iBCDP001 + 1
      integer, parameter :: iBCDP003 = iBCDP002 + 1
      integer, parameter :: iBCDP004 = iBCDP003 + 1
      integer, parameter :: iBCDP005 = iBCDP004 + 1
      integer, parameter :: iBCDP006 = iBCDP005 + 1
      integer, parameter :: iBCDP007 = iBCDP006 + 1
      integer, parameter :: iBCDP008 = iBCDP007 + 1

      integer, parameter :: iBCWT001 = iBCDP008 + 1
      integer, parameter :: iBCWT002 = iBCWT001 + 1
      integer, parameter :: iBCWT003 = iBCWT002 + 1
      integer, parameter :: iBCWT004 = iBCWT003 + 1
      integer, parameter :: iBCWT005 = iBCWT004 + 1
      integer, parameter :: iBCWT006 = iBCWT005 + 1
      integer, parameter :: iBCWT007 = iBCWT006 + 1
      integer, parameter :: iBCWT008 = iBCWT007 + 1

      integer, parameter :: iBCSV001 = iBCWT008 + 1
      integer, parameter :: iBCSV002 = iBCSV001 + 1
      integer, parameter :: iBCSV003 = iBCSV002 + 1
      integer, parameter :: iBCSV004 = iBCSV003 + 1
      integer, parameter :: iBCSV005 = iBCSV004 + 1
      integer, parameter :: iBCSV006 = iBCSV005 + 1
      integer, parameter :: iBCSV007 = iBCSV006 + 1
      integer, parameter :: iBCSV008 = iBCSV007 + 1

      integer, parameter :: iBCSMASS = iBCSV008 + 1
      integer, parameter :: iBCCMASS = iBCSMASS + 1
      integer, parameter :: iBCEXTTAU = iBCCMASS + 1
      integer, parameter :: iBCSCATAU = iBCEXTTAU + 1
      integer, parameter :: iBCEMBF = iBCSCATAU+1
      integer, parameter :: iBCEMBB = iBCEMBF+1
      integer, parameter :: iBCEMAN = iBCEMBB+1
      integer, parameter :: iBCHYPHIL = iBCEMAN+1

      integer, parameter :: iOCEM001 = iBCHYPHIL + 1
      integer, parameter :: iOCEM002 = iOCEM001 + 1
      integer, parameter :: iOCEM003 = iOCEM002 + 1
      integer, parameter :: iOCEM004 = iOCEM003 + 1
      integer, parameter :: iOCEM005 = iOCEM004 + 1
      integer, parameter :: iOCEM006 = iOCEM005 + 1
      integer, parameter :: iOCEM007 = iOCEM006 + 1
      integer, parameter :: iOCEM008 = iOCEM007 + 1

      integer, parameter :: iOCDP001 = iOCEM008 + 1
      integer, parameter :: iOCDP002 = iOCDP001 + 1
      integer, parameter :: iOCDP003 = iOCDP002 + 1
      integer, parameter :: iOCDP004 = iOCDP003 + 1
      integer, parameter :: iOCDP005 = iOCDP004 + 1
      integer, parameter :: iOCDP006 = iOCDP005 + 1
      integer, parameter :: iOCDP007 = iOCDP006 + 1
      integer, parameter :: iOCDP008 = iOCDP007 + 1

      integer, parameter :: iOCWT001 = iOCDP008 + 1
      integer, parameter :: iOCWT002 = iOCWT001 + 1
      integer, parameter :: iOCWT003 = iOCWT002 + 1
      integer, parameter :: iOCWT004 = iOCWT003 + 1
      integer, parameter :: iOCWT005 = iOCWT004 + 1
      integer, parameter :: iOCWT006 = iOCWT005 + 1
      integer, parameter :: iOCWT007 = iOCWT006 + 1
      integer, parameter :: iOCWT008 = iOCWT007 + 1

      integer, parameter :: iOCSV001 = iOCWT008 + 1
      integer, parameter :: iOCSV002 = iOCSV001 + 1
      integer, parameter :: iOCSV003 = iOCSV002 + 1
      integer, parameter :: iOCSV004 = iOCSV003 + 1
      integer, parameter :: iOCSV005 = iOCSV004 + 1
      integer, parameter :: iOCSV006 = iOCSV005 + 1
      integer, parameter :: iOCSV007 = iOCSV006 + 1
      integer, parameter :: iOCSV008 = iOCSV007 + 1

      integer, parameter :: iOCSMASS = iOCSV008 + 1
      integer, parameter :: iOCCMASS = iOCSMASS + 1
      integer, parameter :: iOCEXTTAU = iOCCMASS + 1
      integer, parameter :: iOCSCATAU = iOCEXTTAU + 1
      integer, parameter :: iOCEMBF = iOCSCATAU + 1
      integer, parameter :: iOCEMBB = iOCEMBF + 1
      integer, parameter :: iOCEMAN = iOCEMBB + 1
      integer, parameter :: iOCEMBG = iOCEMAN + 1
      integer, parameter :: iOCHYPHIL = iOCEMBG + 1

      integer, parameter :: iSUEM001 = iOCHYPHIL + 1
      integer, parameter :: iSUEM002 = iSUEM001 + 1
      integer, parameter :: iSUEM003 = iSUEM002 + 1
      integer, parameter :: iSUEM004 = iSUEM003 + 1
      integer, parameter :: iSUEM005 = iSUEM004 + 1
      integer, parameter :: iSUEM006 = iSUEM005 + 1
      integer, parameter :: iSUEM007 = iSUEM006 + 1
      integer, parameter :: iSUEM008 = iSUEM007 + 1

      integer, parameter :: iSUDP001 = iSUEM008 + 1
      integer, parameter :: iSUDP002 = iSUDP001 + 1
      integer, parameter :: iSUDP003 = iSUDP002 + 1
      integer, parameter :: iSUDP004 = iSUDP003 + 1
      integer, parameter :: iSUDP005 = iSUDP004 + 1
      integer, parameter :: iSUDP006 = iSUDP005 + 1
      integer, parameter :: iSUDP007 = iSUDP006 + 1
      integer, parameter :: iSUDP008 = iSUDP007 + 1

      integer, parameter :: iSUWT001 = iSUDP008 + 1
      integer, parameter :: iSUWT002 = iSUWT001 + 1
      integer, parameter :: iSUWT003 = iSUWT002 + 1
      integer, parameter :: iSUWT004 = iSUWT003 + 1
      integer, parameter :: iSUWT005 = iSUWT004 + 1
      integer, parameter :: iSUWT006 = iSUWT005 + 1
      integer, parameter :: iSUWT007 = iSUWT006 + 1
      integer, parameter :: iSUWT008 = iSUWT007 + 1

      integer, parameter :: iSUSV001 = iSUWT008 + 1
      integer, parameter :: iSUSV002 = iSUSV001 + 1
      integer, parameter :: iSUSV003 = iSUSV002 + 1
      integer, parameter :: iSUSV004 = iSUSV003 + 1
      integer, parameter :: iSUSV005 = iSUSV004 + 1
      integer, parameter :: iSUSV006 = iSUSV005 + 1
      integer, parameter :: iSUSV007 = iSUSV006 + 1
      integer, parameter :: iSUSV008 = iSUSV007 + 1

      integer, parameter :: iSUSO2SMASS = iSUSV008 + 1
      integer, parameter :: iSUSO2CMASS = iSUSO2SMASS + 1
      integer, parameter :: iSUSO4SMASS = iSUSO2CMASS + 1
      integer, parameter :: iSUSO4CMASS = iSUSO4SMASS + 1
      integer, parameter :: iSUDMSSMASS = iSUSO4CMASS + 1
      integer, parameter :: iSUDMSCMASS = iSUDMSSMASS + 1
      integer, parameter :: iSUPSO2     = iSUDMSCMASS + 1
      integer, parameter :: iSUPSO4g    = iSUPSO2 + 1
      integer, parameter :: iSUPSO4aq   = iSUPSO4g + 1
      integer, parameter :: iSUPMSA     = iSUPSO4aq + 1
      integer, parameter :: iSUPSO4wet  = iSUPMSA + 1
      integer, parameter :: iSUEXTTAU   = iSUPSO4wet + 1
      integer, parameter :: iSUSCATAU = iSUEXTTAU + 1
      integer, parameter :: iSUEMSO4AN = iSUSCATAU+1
      integer, parameter :: iSUEMSO2AN = iSUEMSO4AN+1
      integer, parameter :: iSUEMSO2BB = iSUEMSO2AN+1
      integer, parameter :: iSUEMSO2VN = iSUEMSO2BB+1
      integer, parameter :: iSUEMSO2VE = iSUEMSO2VN+1

      INTEGER, PARAMETER :: iN2OFLX = iSUEMSO2VE + 1
      INTEGER, PARAMETER :: iCH4FLX = iN2OFLX + 1
      INTEGER, PARAMETER :: iF11FLX = iCH4FLX + 1
      INTEGER, PARAMETER :: iF12FLX = iF11FLX + 1
      INTEGER, PARAMETER :: iF113FLX = iF12FLX + 1
      INTEGER, PARAMETER :: iHCFCFLX = iF113FLX + 1
      INTEGER, PARAMETER :: iCCL4FLX = iHCFCFLX + 1
      INTEGER, PARAMETER :: iCH3CCL3FLX = iCCL4FLX + 1
      INTEGER, PARAMETER :: iCH3CLFLX = iCH3CCL3FLX + 1
      INTEGER, PARAMETER :: iCH3BRFLX = iCH3CLFLX + 1
      INTEGER, PARAMETER :: iH1301FLX = iCH3BRFLX + 1
      INTEGER, PARAMETER :: iH12_24FLX = iH1301FLX + 1

      integer, parameter :: iCOSSZA = iH12_24FLX + 1
      integer, parameter :: iTCOSZ = iCOSSZA + 1
      integer, parameter :: iXOH = iTCOSZ + 1
      integer, parameter :: iXNO3 = iXOH + 1
      integer, parameter :: iXH2O2 = iXNO3 + 1

      integer, parameter :: iCOEM001 = iXH2O2 + 1
      integer, parameter :: iCOEM002 = iCOEM001 + 1
      integer, parameter :: iCOEM003 = iCOEM002 + 1
      integer, parameter :: iCOEM004 = iCOEM003 + 1
      integer, parameter :: iCOEM005 = iCOEM004 + 1
      integer, parameter :: iCOEM006 = iCOEM005 + 1
      integer, parameter :: iCOEM007 = iCOEM006 + 1
      integer, parameter :: iCOEM008 = iCOEM007 + 1

      integer, parameter :: iCOLS001 = iCOEM008 + 1
      integer, parameter :: iCOLS002 = iCOLS001 + 1
      integer, parameter :: iCOLS003 = iCOLS002 + 1
      integer, parameter :: iCOLS004 = iCOLS003 + 1
      integer, parameter :: iCOLS005 = iCOLS004 + 1
      integer, parameter :: iCOLS006 = iCOLS005 + 1
      integer, parameter :: iCOLS007 = iCOLS006 + 1
      integer, parameter :: iCOLS008 = iCOLS007 + 1

      integer, parameter :: iCOPD001 = iCOLS008 + 1
      integer, parameter :: iCOPD002 = iCOPD001 + 1
      integer, parameter :: iCOPD003 = iCOPD002 + 1
      integer, parameter :: iCOPD004 = iCOPD003 + 1
      integer, parameter :: iCOPD005 = iCOPD004 + 1
      integer, parameter :: iCOPD006 = iCOPD005 + 1
      integer, parameter :: iCOPD007 = iCOPD006 + 1
      integer, parameter :: iCOPD008 = iCOPD007 + 1

      integer, parameter :: iCOCL001 = iCOPD008 + 1
      integer, parameter :: iCOCL002 = iCOCL001 + 1
      integer, parameter :: iCOCL003 = iCOCL002 + 1
      integer, parameter :: iCOCL004 = iCOCL003 + 1
      integer, parameter :: iCOCL005 = iCOCL004 + 1
      integer, parameter :: iCOCL006 = iCOCL005 + 1
      integer, parameter :: iCOCL007 = iCOCL006 + 1
      integer, parameter :: iCOCL008 = iCOCL007 + 1

      integer, parameter :: iCOSC001 = iCOCL008 + 1
      integer, parameter :: iCOSC002 = iCOSC001 + 1
      integer, parameter :: iCOSC003 = iCOSC002 + 1
      integer, parameter :: iCOSC004 = iCOSC003 + 1
      integer, parameter :: iCOSC005 = iCOSC004 + 1
      integer, parameter :: iCOSC006 = iCOSC005 + 1
      integer, parameter :: iCOSC007 = iCOSC006 + 1
      integer, parameter :: iCOSC008 = iCOSC007 + 1

      integer, parameter :: iCO2EM001 = iCOSC008 + 1
      integer, parameter :: iCO2EM002 = iCO2EM001 + 1
      integer, parameter :: iCO2EM003 = iCO2EM002 + 1
      integer, parameter :: iCO2EM004 = iCO2EM003 + 1
      integer, parameter :: iCO2EM005 = iCO2EM004 + 1
      integer, parameter :: iCO2EM006 = iCO2EM005 + 1
      integer, parameter :: iCO2EM007 = iCO2EM006 + 1
      integer, parameter :: iCO2EM008 = iCO2EM007 + 1

      integer, parameter :: iCO2CL001 = iCO2EM008 + 1
      integer, parameter :: iCO2CL002 = iCO2CL001 + 1
      integer, parameter :: iCO2CL003 = iCO2CL002 + 1
      integer, parameter :: iCO2CL004 = iCO2CL003 + 1
      integer, parameter :: iCO2CL005 = iCO2CL004 + 1
      integer, parameter :: iCO2CL006 = iCO2CL005 + 1
      integer, parameter :: iCO2CL007 = iCO2CL006 + 1
      integer, parameter :: iCO2CL008 = iCO2CL007 + 1

      integer, parameter :: iCO2SC001 = iCO2CL008 + 1
      integer, parameter :: iCO2SC002 = iCO2SC001 + 1
      integer, parameter :: iCO2SC003 = iCO2SC002 + 1
      integer, parameter :: iCO2SC004 = iCO2SC003 + 1
      integer, parameter :: iCO2SC005 = iCO2SC004 + 1
      integer, parameter :: iCO2SC006 = iCO2SC005 + 1
      integer, parameter :: iCO2SC007 = iCO2SC006 + 1
      integer, parameter :: iCO2SC008 = iCO2SC007 + 1

      parameter (pd2d      = iCO2SC008 )

#else
      parameter (pd2d      = iZPD )
#endif

!
! ... 3-D diagnosis fields
!
      parameter (iAIRDENS  = pd2d     + 1)
      parameter (iCAPE     = iAIRDENS + 1)
      parameter (iCGS      = iCAPE    + 1)
      parameter (iCLDLWP   = iCGS     + 1)
      parameter (iCLOUD    = iCLDLWP  + 1)
      parameter (iCLOUDUP  = iCLOUD   + 1)
      parameter (iCMFDQ    = iCLOUDUP + 1)
      parameter (iCMFDQR2  = iCMFDQ   + 1)
      parameter (iCMFDT    = iCMFDQR2 + 1)
      parameter (iCMFDTR   = iCMFDT   + 1)
      parameter (iCMFETR   = iCMFDTR  + 1)
      parameter (iCMFMC    = iCMFETR  + 1)
      parameter (iCMFMC2   = iCMFMC   + 1)
      parameter (iCONVCLD  = iCMFMC2  + 1)
      parameter (iDCAFDT   = iCONVCLD + 1)
      parameter (iDIABDT   = iDCAFDT  + 1)
      parameter (iDQCOND   = iDIABDT  + 1)
      parameter (iDQPBLCG  = iDQCOND  + 1)
      parameter (iDQRL     = iDQPBLCG + 1)
      parameter (iDTCOND   = iDQRL    + 1)
      parameter (iDTPBLCG  = iDTCOND  + 1)
      parameter (iDTRAIN   = iDTPBLCG + 1)
      parameter (iDTV      = iDTRAIN  + 1)
      parameter (iDUV      = iDTV     + 1)
      parameter (iDVV      = iDUV     + 1)
      parameter (iEFFCLD   = iDVV     + 1)
      parameter (iEVAPL    = iEFFCLD  + 1)
      parameter (iH        = iEVAPL   + 1)
      parameter (iHGHTE    = iH       + 1)
      parameter (iHKBETA   = iHGHTE   + 1)
      parameter (iHKETA    = iHKBETA  + 1)
      parameter (iKVH      = iHKETA   + 1)
      parameter (iKVM      = iKVH     + 1)
      parameter (iO3VMR    = iKVM     + 1)
      parameter (iOMEGA    = iO3VMR   + 1)
      parameter (iPV       = iOMEGA   + 1)
      parameter (iQ        = iPV      + 1)
      parameter (iQC       = iQ       + 1)
      parameter (iQRL      = iQC      + 1)
      parameter (iQRS      = iQRL     + 1)
      parameter (iRAYFDT   = iQRS     + 1)
      parameter (iRELHUM   = iRAYFDT  + 1)
      parameter (iRHCLR    = iRELHUM  + 1)
      parameter (iRNEVPDQ  = iRHCLR   + 1)
      parameter (iRNEVPDT  = iRNEVPDQ + 1)
      parameter (iSETLWP   = iRNEVPDT + 1)
      parameter (iSTRATCLD = iSETLWP  + 1)
      parameter (iT        = iSTRATCLD+ 1)
      parameter (iTAUCLI   = iT       + 1)
      parameter (iTAUCLW   = iTAUCLI  + 1)
      parameter (iTKE      = iTAUCLW  + 1)
      parameter (iTTMGW    = iTKE     + 1)
      parameter (iU        = iTTMGW   + 1)
      parameter (iUQ       = iU       + 1)
      parameter (iUT       = iUQ      + 1)
      parameter (iUTGW     = iUT      + 1)
      parameter (iUU       = iUTGW    + 1)
      parameter (iUV       = iUU      + 1)
      parameter (iV        = iUV      + 1)
      parameter (iVD01     = iV       + 1)
      parameter (iVQ       = iVD01    + 1)
      parameter (iVT       = iVQ      + 1)
      parameter (iVTGW     = iVT      + 1)
      parameter (iVV       = iVTGW    + 1)
      parameter (iZMCME    = iVV      + 1)
      parameter (iZMDLF    = iZMCME   + 1)
      parameter (iZMDQ     = iZMDLF   + 1)
      parameter (iZMDQR    = iZMDQ    + 1)
      parameter (iZMDT     = iZMDQR   + 1)
      parameter (iZMDU     = iZMDT    + 1)
      parameter (iZMED     = iZMDU    + 1)
      parameter (iZMEPS    = iZMED    + 1)
      parameter (iZMEU     = iZMEPS   + 1)
      parameter (iZMEVP    = iZMEU    + 1)
      parameter (iZMMD     = iZMEVP   + 1)
      parameter (iZMMU     = iZMMD    + 1)
      parameter (iZMPFLX   = iZMMU    + 1)
      parameter (iZMQL     = iZMPFLX  + 1)

#ifdef FVCHEM
      INTEGER, PARAMETER :: iBR      = iZMQL     + 1
      INTEGER, PARAMETER :: iBRCL    = iBR + 1
      INTEGER, PARAMETER :: iBRO     = iBRCL + 1
      INTEGER, PARAMETER :: iBRONO2  = iBRO + 1
      INTEGER, PARAMETER :: iBRX     = iBRONO2 + 1
      INTEGER, PARAMETER :: iCCL4    = iBRX + 1
      INTEGER, PARAMETER :: iCH2O    = iCCL4 + 1
      INTEGER, PARAMETER :: iCH3BR   = iCH2O + 1
      INTEGER, PARAMETER :: iCH3CCL3 = iCH3BR + 1
      INTEGER, PARAMETER :: iCH3CL   = iCH3CCL3 + 1
      INTEGER, PARAMETER :: iCH3O2   = iCH3CL + 1
      INTEGER, PARAMETER :: iCH3OOH  = iCH3O2 + 1
      INTEGER, PARAMETER :: iCH4     = iCH3OOH + 1
      INTEGER, PARAMETER :: iCL      = iCH4 + 1
      INTEGER, PARAMETER :: iCL2     = iCL + 1
      INTEGER, PARAMETER :: iCL2O2   = iCL2 + 1
      INTEGER, PARAMETER :: iCLO     = iCL2O2 + 1
      INTEGER, PARAMETER :: iCLONO2  = iCLO + 1
      INTEGER, PARAMETER :: iCLX     = iCLONO2 + 1

! --------------------------------------------------------------------
!          8 CO regions and types for INTEX-B 2006
! --------------------------------------------------------------------

      INTEGER, PARAMETER :: iCO      = iCLX + 1
      INTEGER, PARAMETER :: iCONOAMAN= iCO + 1
      INTEGER, PARAMETER :: iCOCEAMAN= iCONOAMAN + 1
      INTEGER, PARAMETER :: iCOWHBB  = iCOCEAMAN + 1
      INTEGER, PARAMETER :: iCOASIAAN= iCOWHBB + 1
      INTEGER, PARAMETER :: iCOASNBB = iCOASIAAN + 1
      INTEGER, PARAMETER :: iCOASSBB = iCOASNBB + 1
      INTEGER, PARAMETER :: iCOFDAN  = iCOASSBB + 1

! --------------------------------------------------------------------

      INTEGER, PARAMETER :: iCOFF    = iCOFDAN  + 1
      INTEGER, PARAMETER :: iCOBF    = iCOFF + 1
      INTEGER, PARAMETER :: iCOBB    = iCOBF + 1
      INTEGER, PARAMETER :: iCOBI    = iCOBB + 1
      INTEGER, PARAMETER :: iCONAMERI= iCOBI  + 1
      INTEGER, PARAMETER :: iCOSAMERI= iCONAMERI + 1
      INTEGER, PARAMETER :: iCOAFRICA= iCOSAMERI + 1
      INTEGER, PARAMETER :: iCO2     = iCOAFRICA + 1
      INTEGER, PARAMETER :: iCO2NAMER= iCO2 + 1
      INTEGER, PARAMETER :: iCO2SAMER= iCO2NAMER + 1
      INTEGER, PARAMETER :: iCO2AFRIC= iCO2SAMER + 1
      INTEGER, PARAMETER :: iF11     = iCO2AFRIC + 1
      INTEGER, PARAMETER :: iF113    = iF11 + 1
      INTEGER, PARAMETER :: iF12     = iF113 + 1
      INTEGER, PARAMETER :: iH12_24  = iF12 + 1
      INTEGER, PARAMETER :: iH1301   = iH12_24 + 1
      INTEGER, PARAMETER :: iH2O2    = iH1301 + 1
      INTEGER, PARAMETER :: iH2OCOND = iH2O2 + 1
      INTEGER, PARAMETER :: iHATOMIC = iH2OCOND + 1
      INTEGER, PARAMETER :: iHBR     = iHATOMIC + 1
      INTEGER, PARAMETER :: iHCFC    = iHBR + 1
      INTEGER, PARAMETER :: iHCL     = iHCFC + 1
      INTEGER, PARAMETER :: iHNO3    = iHCL + 1
      INTEGER, PARAMETER :: iHNO3COND= iHNO3 + 1
      INTEGER, PARAMETER :: iHO2     = iHNO3COND + 1
      INTEGER, PARAMETER :: iHO2NO2  = iHO2 + 1
      INTEGER, PARAMETER :: iHOBR    = iHO2NO2 + 1
      INTEGER, PARAMETER :: iHOCL    = iHOBR + 1
      INTEGER, PARAMETER :: iN       = iHOCL + 1
      INTEGER, PARAMETER :: iN2O     = iN + 1
      INTEGER, PARAMETER :: iN2O5    = iN2O + 1
      INTEGER, PARAMETER :: iNO      = iN2O5 + 1
      INTEGER, PARAMETER :: iNO2     = iNO + 1
      INTEGER, PARAMETER :: iNO3     = iNO2 + 1
      INTEGER, PARAMETER :: iNOX     = iNO3 + 1
      INTEGER, PARAMETER :: iO1D     = iNOX + 1
      INTEGER, PARAMETER :: iO3CHEM  = iO1D + 1
      INTEGER, PARAMETER :: iO3P     = iO3CHEM + 1
      INTEGER, PARAMETER :: iO3PARAM = iO3P + 1
      INTEGER, PARAMETER :: iOCLO    = iO3PARAM + 1
      INTEGER, PARAMETER :: iOH      = iOCLO + 1
      INTEGER, PARAMETER :: iOX      = iOH + 1
      INTEGER, PARAMETER :: iOXSTRAT = iOX + 1
      INTEGER, PARAMETER :: iOXTROP  = iOXSTRAT + 1

      integer, parameter :: iDUMASS  = iOXTROP + 1
      integer, parameter :: iDUMASS1  = iDUMASS + 1
      integer, parameter :: iDUMASS25  = iDUMASS1 + 1
      integer, parameter :: iSSMASS  = iDUMASS25 + 1
      integer, parameter :: iSSMASS25  = iSSMASS + 1
      integer, parameter :: iBCMASS  = iSSMASS25 + 1

      integer, parameter :: iOCMASS  = iBCMASS + 1
      integer, parameter :: iSO4MASS = iOCMASS + 1
      integer, parameter :: iPSO2    = iSO4MASS + 1
      integer, parameter :: iPSO4g   = iPSO2 + 1
      integer, parameter :: iPSO4aq  = iPSO4g + 1
      integer, parameter :: iPMSA    = iPSO4aq + 1
      integer, parameter :: iPSO4wet = iPMSA + 1
      integer, parameter :: iQ4AGE   = iPSO4wet + 1
      parameter (pdiag     = iQ4AGE)
#else
      parameter (pdiag     = iZMQL)
#endif

      parameter (pd3d      = pdiag  - pd2d)

      integer maxlist
      parameter (maxlist = 30)	! maximum number of output files

!      real*4 			:: undef = 1.e25
      real 			:: undef = 1.e25

! There is one such diagnostic descriptor (diag_type) for each diagnostic
! which can be output using outfld. Each diagnostic can be output to multiple
! output streams (up to nlist of them, the index of % count and % alt and the 
! second index of % fldloc) and for each output stream different types of
! the diagnostic may be output, currently averaged, instantaneous, maximum 
! and minimum (the first index of % fldloc).

      type diag_type
           character*8 		:: name        	! name for diagnostic
           character*16 	:: unit        	! unit for diagnostic fields
           character*80 	:: desc        	! description for diagnostic fields
           integer              :: pick         ! is the diag needed for output 1:yes 0:no
           logical      	:: counted      ! if true, diagnostic is counted 
           integer		:: vdim		! number of levels 1 or km
           integer              :: nlist	! number of output streams
           integer, pointer 	:: fldloc(:,:)	! location in diagnostic buffer (negative if not on)
           integer, pointer 	:: count(:)	! counter 
           logical, pointer 	:: alt(:)	! if true use alternate name/unit for output
           character*8 		:: aname       	! Alternate name for diagnostic
           character*16 	:: aunit       	! Alternate unit for diagnostic fields
           character*80         :: adesc        ! description for diagnostic fields
           real*4       	:: convfac     	! conversion factor if alternate unit different than primary
      endtype 

! Each output stream (history file) has a hist_type descriptor. Each stream has one
! frequency and any number of subset regions of the model grid (nreg of them, the array
! index in most of the fields). A positive frequency gives the period of periodic output.
! A zero frequency is used for segment end output. As of 7-16-03 we allow a negative
! frequency (-1) to indicate a set of non-output diagnostics. Such diagnostics (usually
! instantaneous) are maintained in the same way as a normal diagnostics but they are
! never output to file. They can be accessed through getdiag (). 

      type hist_type
           character*80, pointer :: name(:)	! Output Filenames for each region
           integer             	:: frequency    ! frequency for output stream HHMMSS
           character*1          :: timestamp    ! Time stamp for output stream Centered, Forward, Backward
           character*3,pointer  :: type(:)  	! file type for output hdf or bin, hdf is default
           integer, pointer	:: imin(:)	! longitude minimum bounds for each region
           integer, pointer	:: imax(:)	! longitude maximum bounds for each region
           integer, pointer	:: jmin(:)	! latitude minimum bounds for each region
           integer, pointer	:: jmax(:)	! latitude maximum bounds for each region
           integer, pointer	:: ktop(:)	! vertical top bounds for each region
           integer, pointer	:: kbot(:)	! vertical bottom bounds for each region
           integer, pointer 	:: nrec(:)      ! diag record number (for binary not hdf output)
           integer, pointer 	:: lun(:)       ! logical unit numbers (binary, grads ctl reuses lun)
           integer          	:: nreg    	! Number of regions for the output stream
      endtype 

      logical doingDiag

      END MODULE mod_diag

