MODULE Henrys_law_ConstantsMod
  IMPLICIT NONE
  !--- : ak0(ispc), dak(ispc),  hstar(ispc), dhr(ispc)
  !--- corrh=1.+ak0(ispc)*exp(dak(ispc)*tcorr)/hplus  
  !--- hplus = 1.175E-4  - for cloud water. pH is assumed to be 3.93: pH=3.93 =>hplus=10**(-pH)  
  !--- tcorr = 1./temp - 1./298.15
  !--- fct   = 1.e-3 * rgas * temp  
  !--- henry_coef =  hstar(ispc)* exp(dhr(ispc)*tcorr) * fct * corrh
  
  INTEGER,PARAMETER :: nspecies_HL=051
  REAL   ,PARAMETER :: notfound = -1.
  
  type Hcts_vars
   real :: hstar,dhr,ak0,dak
  END TYPE Hcts_vars
  type (Hcts_vars) :: Hcts(nspecies_HL)

  !Name of species 
  CHARACTER(LEN=8),PARAMETER,DIMENSION(nspecies_HL) :: spc_name=(/ &
      'O3  ' & !001
     ,'H2O2' & !002
     ,'NO  ' & !003
     ,'NO2 ' & !004
     ,'NO3 ' & !005
     ,'N2O5' & !006
     ,'HONO' & !007
     ,'HNO3' & !008
     ,'HNO4' & !009
     ,'SO2 ' & !010
     ,'SULF' & !011
     ,'CO  ' & !012
     ,'CO2 ' & !013
     ,'N2  ' & !014
     ,'O2  ' & !015
     ,'H2O ' & !016
     ,'H2  ' & !017
     ,'O3P ' & !018
     ,'O1D ' & !019
     ,'HO  ' & !020
     ,'HO2 ' & !021
     ,'CH4 ' & !022
     ,'ETH ' & !023
     ,'ALKA' & !024
     ,'ALKE' & !025
     ,'BIO ' & !026
     ,'ARO ' & !027
     ,'HCHO' & !028
     ,'ALD ' & !029
     ,'KET ' & !030
     ,'CRBO' & !031
     ,'ONIT' & !032
     ,'PAN ' & !033
     ,'OP1 ' & !034
     ,'OP2 ' & !035
     ,'ORA1' & !036
     ,'ORA2' & !037
     ,'MO2 ' & !038
     ,'AKAP' & !039
     ,'AKEP' & !040
     ,'BIOP' & !041
     ,'PHO ' & !042
     ,'ADD ' & !043
     ,'AROP' & !044
     ,'CBOP' & !045
     ,'OLN ' & !046
     ,'XO2 ' & !047
     ,'DMS ' & !048
     ,'NH3 ' & !049
     ,'CFC ' & !050
     ,'N2O ' & !050
   /)
  
  
  !Number of each specie   
  INTEGER,PARAMETER :: O3  =001
  INTEGER,PARAMETER :: H2O2=002
  INTEGER,PARAMETER :: NO  =003
  INTEGER,PARAMETER :: NO2 =004
  INTEGER,PARAMETER :: NO3 =005
  INTEGER,PARAMETER :: N2O5=006
  INTEGER,PARAMETER :: HONO=007
  INTEGER,PARAMETER :: HNO3=008
  INTEGER,PARAMETER :: HNO4=009
  INTEGER,PARAMETER :: SO2 =010
  INTEGER,PARAMETER :: SULF=011
  INTEGER,PARAMETER :: CO  =012
  INTEGER,PARAMETER :: CO2 =013
  INTEGER,PARAMETER :: N2  =014
  INTEGER,PARAMETER :: O2  =015
  INTEGER,PARAMETER :: H2O =016
  INTEGER,PARAMETER :: H2  =017
  INTEGER,PARAMETER :: O3P =018
  INTEGER,PARAMETER :: O1D =019
  INTEGER,PARAMETER :: HO  =020
  INTEGER,PARAMETER :: HO2 =021
  INTEGER,PARAMETER :: CH4 =022
  INTEGER,PARAMETER :: ETH =023
  INTEGER,PARAMETER :: ALKA=024
  INTEGER,PARAMETER :: ALKE=025
  INTEGER,PARAMETER :: BIO =026
  INTEGER,PARAMETER :: ARO =027
  INTEGER,PARAMETER :: HCHO=028
  INTEGER,PARAMETER :: ALD =029
  INTEGER,PARAMETER :: KET =030
  INTEGER,PARAMETER :: CRBO=031
  INTEGER,PARAMETER :: ONIT=032
  INTEGER,PARAMETER :: PAN =033
  INTEGER,PARAMETER :: OP1 =034
  INTEGER,PARAMETER :: OP2 =035
  INTEGER,PARAMETER :: ORA1=036
  INTEGER,PARAMETER :: ORA2=037
  INTEGER,PARAMETER :: MO2 =038
  INTEGER,PARAMETER :: AKAP=039
  INTEGER,PARAMETER :: AKEP=040
  INTEGER,PARAMETER :: BIOP=041
  INTEGER,PARAMETER :: PHO =042
  INTEGER,PARAMETER :: ADD =043
  INTEGER,PARAMETER :: AROP=044
  INTEGER,PARAMETER :: CBOP=045
  INTEGER,PARAMETER :: OLN =046
  INTEGER,PARAMETER :: XO2 =047
  INTEGER,PARAMETER :: DMS =048
  INTEGER,PARAMETER :: NH3 =049
  INTEGER,PARAMETER :: CFC =050
  INTEGER,PARAMETER :: N2O =051
 
  
  
  
!     HENRYS LAW COEFFICIENTS
!     Henrys law coefficient
!     [KH298]=mole/(l atm)
!     Referencias em R. Sander (1999)
!     Compilation of Henry Law Constants 
!     for Inorganic and Organic Species 
!     of Potential Importance in 
!     Environmental Chemistry (Version 3) 
!     http://www.henrys-law.org 
!     * indica artigos nao encontrados nesse endereço eletronico
  REAL,PARAMETER,DIMENSION(nspecies_HL) :: hstar=(/&
    1.10E-2              ,   & ! O3 - 001
    8.30E+4              ,   & ! H2O2 - 002
    1.90E-3              ,   & ! NO - 003
    1.20E-2              ,   & ! NO2 - 004
    6.1E-01              ,   & ! NO3 - 005
    2.1E+00              ,   & ! N2O5 - 006
    5.00E+1              ,   & ! HONO - 007
    2.10E+5              ,   & ! HNO3 - 008
    1.20E+4              ,   & ! HNO4 - 009
    1.40E+0              ,   & ! SO2 - 010
    2.10E+5              ,   & ! SULF - 011
    9.90E-4              ,   & ! CO - 012
    3.6E-02              ,   & ! CO2 - 013
    6.1E-04              ,   & ! N2 - 014
    1.3E-03              ,   & ! O2 - 015
    0.0E+00              ,   & ! H2O - 016
    7.8E-04              ,   & ! H2 - 017
    0.00E+0              ,   & ! O3P - 018
    0.00E+0              ,   & ! O1D - 019
    3.00E+1              ,   & ! HO - 020
    5.70E+3              ,   & ! HO2 - 021
    1.40E-3              ,   & ! CH4 - 022
    1.90E-3              ,   & ! ETH - 023
    1.00E-3              ,   & ! ALKA - 024
    5.00E-3              ,   & ! ALKE - 025
    2.80E-2              ,   & ! BIO - 026
    1.73E-1              ,   & ! ARO - 027
    3.20E+3              ,   & ! HCHO - 028
    1.40E+1              ,   & ! ALD - 029
    3.00E+1              ,   & ! KET - 030
    2.1E+05              ,   & ! CRBO - 031
    1.00E+0              ,   & ! ONIT - 032
    3.60E+0              ,   & ! PAN - 033
    3.10E+2              ,   & ! OP1 - 034
    3.40E+2              ,   & ! OP2 - 035
    8.90E+3              ,   & ! ORA1 - 036
    4.10E+3              ,   & ! ORA2 - 037
    2.00E+3              ,   & ! MO2 - 038
    0.0E+00              ,   & ! AKAP - 039
    0.0E+00              ,   & ! AKEP - 040
    0.0E+00              ,   & ! BIOP - 041
    0.0E+00              ,   & ! PHO - 042
    0.0E+00              ,   & ! ADD - 043
    0.0E+00              ,   & ! AROP - 044
    1.14E+1              ,   & ! CBOP - 045
    0.0E+00              ,   & ! OLN - 046
    0.0E+00              ,   & ! XO2 - 047
    5.6E-01              ,   & ! DMS - 048
    5.9E+01              ,   & ! NH3 - 048
    -1.                  ,   & ! CFC - 048
    2.4E-02                  & ! N2O - 051
    /)
  
  
  
!     -DH/R (for temperature correction)
!     [-DH/R]=K
!     Referencias em R. Sander (1999)
!     Compilation of Henry Law Constants
!     for Inorganic and Organic Species 
!     of Potential Importance in 
!     Environmental Chemistry (Version 3)
!     http://www.henrys-law.org 
  REAL,PARAMETER,DIMENSION(nspecies_HL) :: dhr=(/&
    2400.         ,   & ! O3 - 001
    7400.         ,   & ! H2O2 - 002
    1400.         ,   & ! NO - 003
    2500.         ,   & ! NO2 - 004
    2000.         ,   & ! NO3 - 005
    3400.         ,   & ! N2O5 - 006
    4900.         ,   & ! HONO - 007
    8700.         ,   & ! HNO3 - 008
    6900.         ,   & ! HNO4 - 009
    2900.         ,   & ! SO2 - 010
    0.            ,   & ! SULF - 011
    1300.         ,   & ! CO - 012
    2200.         ,   & ! CO2 - 013
    1300.         ,   & ! N2 - 014
    1500.         ,   & ! O2 - 015
    0.            ,   & ! H2O - 016
    500.          ,   & ! H2 - 017
    0.            ,   & ! O3P - 018
    0.            ,   & ! O1D - 019
    4500.         ,   & ! HO - 020
    5900.         ,   & ! HO2 - 021
    1600.         ,   & ! CH4 - 022
    2300.         ,   & ! ETH - 023
    2700.         ,   & ! ALKA - 024
    3000.         ,   & ! ALKE - 025
    0.            ,   & ! BIO - 026
    4045.         ,   & ! ARO - 027
    6800.         ,   & ! HCHO - 028
    5600.         ,   & ! ALD - 029
    4600.         ,   & ! KET - 030
    5300.         ,   & ! CRBO - 031
    5800.         ,   & ! ONIT - 032
    6500.         ,   & ! PAN - 033
    5200.         ,   & ! OP1 - 034
    6000.         ,   & ! OP2 - 035
    5700.         ,   & ! ORA1 - 036
    6300.         ,   & ! ORA2 - 037
    6600.         ,   & ! MO2 - 038
    0.            ,   & ! AKAP - 039
    0.            ,   & ! AKEP - 040
    0.            ,   & ! BIOP - 041
    0.            ,   & ! PHO - 042
    0.            ,   & ! ADD - 043
    0.            ,   & ! AROP - 044
    0.            ,   & ! CBOP - 045
    0.            ,   & ! OLN - 046
    0.            ,   & ! XO2 - 047
    3500.         ,   & ! DMS - 048
    4200.         ,   & ! NH3 - 048
    -1.           ,   & ! CFC - 048
    2700.             & ! N2O - 048
    /)
  
  
  REAL,PARAMETER,DIMENSION(nspecies_HL) :: weight=(/&
    48.  ,   & ! O3 - 001
    34.  ,   & ! H2O2 - 002
    30.  ,   & ! NO - 003
    46.  ,   & ! NO2 - 004
    62.  ,   & ! NO3 - 005
    108. ,   & ! N2O5 - 006
    47.  ,   & ! HONO - 007
    63.  ,   & ! HNO3 - 008
    79.  ,   & ! HNO4 - 009
    64.  ,   & ! SO2 - 010
    98.  ,   & ! SULF - 011
    28.  ,   & ! CO - 012
    44.  ,   & ! CO2 - 013
    28.  ,   & ! N2 - 014
    32.  ,   & ! O2 - 015
    18.  ,   & ! H2O - 016
    2.   ,   & ! H2 - 017
    16.  ,   & ! O3P - 018
    16.  ,   & ! O1D - 019
    17.  ,   & ! HO - 020
    33.  ,   & ! HO2 - 021
    16.  ,   & ! CH4 - 022
    30.  ,   & ! ETH - 023
    61.6 ,   & ! ALKA - 024
    33.0 ,   & ! ALKE - 025
    68.  ,   & ! BIO - 026
    97.9 ,   & ! ARO - 027
    30.  ,   & ! HCHO - 028
    44.  ,   & ! ALD - 029
    72.  ,   & ! KET - 030
    68.6 ,   & ! CRBO - 031
    119. ,   & ! ONIT - 032
    122. ,   & ! PAN - 033
    48.  ,   & ! OP1 - 034
    62.  ,   & ! OP2 - 035
    46.  ,   & ! ORA1 - 036
    60.  ,   & ! ORA2 - 037
    47.  ,   & ! MO2 - 038
    102. ,   & ! AKAP - 039
    88.4 ,   & ! AKEP - 040
    117. ,   & ! BIOP - 041
    107. ,   & ! PHO - 042
    107. ,   & ! ADD - 043
    151. ,   & ! AROP - 044
    85.4 ,   & ! CBOP - 045
    136. ,   & ! OLN - 046
    44.  ,   & ! XO2 - 047
    62.13,   & ! DMS - 048
    17.03,   & ! NH3 - 048
    -1.  ,   & ! CFC - 048
    44.      & ! CFC - 048
   /)
  
  
!    ACID DISSOCIATION CONSTANT AT 298K 
!     [mole/liter of liquid water]
!     Referencias: Barth et al. JGR 112, D13310 2007
!     Martell and Smith, 1976, Critical stability
!     vol1-4 Plenum Press New York
  REAL,PARAMETER,DIMENSION(nspecies_HL) :: ak0=(/&
    0.00E+00     ,   & ! O3 - 001
    2.20E-12     ,   & ! H2O2 - 002
    0.00E+00     ,   & ! NO - 003
    0.00E+00     ,   & ! NO2 - 004
    0.00E+00     ,   & ! NO3 - 005
    0.00E+00     ,   & ! N2O5 - 006
    7.10E-04     ,   & ! HONO - 007
    1.54E+01     ,   & ! HNO3 - 008
    0.00E+00     ,   & ! HNO4 - 009
    1.30E-02     ,   & ! SO2 - 010
    1.00E-02     ,   & ! SULF - 011
    0.00E+00     ,   & ! CO - 012
    4.50E-07     ,   & ! CO2 - 013
    0.00E+00     ,   & ! N2 - 014
    0.00E+00     ,   & ! O2 - 015
    0.00E+00     ,   & ! H2O - 016
    0.00E+00     ,   & ! H2 - 017
    0.00E+00     ,   & ! O3P - 018
    0.00E+00     ,   & ! O1D - 019
    0.00E+00     ,   & ! HO - 020
    3.50E-05     ,   & ! HO2 - 021
    0.00E+00     ,   & ! CH4 - 022
    0.00E+00     ,   & ! ETH - 023
    0.00E+00     ,   & ! ALKA - 024
    0.00E+00     ,   & ! ALKE - 025
    0.00E+00     ,   & ! BIO - 026
    0.00E+00     ,   & ! ARO - 027
    0.00E+00     ,   & ! HCHO - 028
    0.00E+00     ,   & ! ALD - 029
    0.00E+00     ,   & ! KET - 030
    0.00E+00     ,   & ! CRBO - 031
    0.00E+00     ,   & ! ONIT - 032
    0.00E+00     ,   & ! PAN - 033
    0.00E+00     ,   & ! OP1 - 034
    0.00E+00     ,   & ! OP2 - 035
    1.80E-04     ,   & ! ORA1 - 036
    1.75E-05     ,   & ! ORA2 - 037
    0.00E+00     ,   & ! MO2 - 038
    0.00E+00     ,   & ! AKAP - 039
    0.00E+00     ,   & ! AKEP - 040
    0.00E+00     ,   & ! BIOP - 041
    0.00E+00     ,   & ! PHO - 042
    0.00E+00     ,   & ! ADD - 043
    0.00E+00     ,   & ! AROP - 044
    0.00E+00     ,   & ! CBOP - 045
    0.00E+00     ,   & ! OLN - 046
    0.00E+00     ,   & ! XO2 - 047
    0.00E+00     ,   & ! DMS - 048
    0.00E+00     ,   & ! NH3 - 049
    0.00E+00     ,   & ! NH3 - 049
    0.00E+00         & ! CFC - 050
   /)
  
!     Temperature correction factor for
!     acid dissociation constants
!     [K]
!     Referencias: Barth et al. JGR 112, D13310 2007
  REAL,PARAMETER,DIMENSION(nspecies_HL) :: dak=(/&
    0.         ,   & ! O3 - 001
    -3700.     ,   & ! H2O2 - 002
    0.         ,   & ! NO - 003
    0.         ,   & ! NO2 - 004
    0.         ,   & ! NO3 - 005
    0.         ,   & ! N2O5 - 006
    0.         ,   & ! HONO - 007
    0.         ,   & ! HNO3 - 008
    0.         ,   & ! HNO4 - 009
    2000.      ,   & ! SO2 - 010
    0.         ,   & ! SULF - 011
    0.         ,   & ! CO - 012
    -1000.     ,   & ! CO2 - 013
    0.         ,   & ! N2 - 014
    0.         ,   & ! O2 - 015
    0.         ,   & ! H2O - 016
    0.         ,   & ! H2 - 017
    0.         ,   & ! O3P - 018
    0.         ,   & ! O1D - 019
    0.         ,   & ! HO - 020
    0.         ,   & ! HO2 - 021
    0.         ,   & ! CH4 - 022
    0.         ,   & ! ETH - 023
    0.         ,   & ! ALKA - 024
    0.         ,   & ! ALKE - 025
    0.         ,   & ! BIO - 026
    0.         ,   & ! ARO - 027
    0.         ,   & ! HCHO - 028
    0.         ,   & ! ALD - 029
    0.         ,   & ! KET - 030
    0.         ,   & ! CRBO - 031
    0.         ,   & ! ONIT - 032
    0.         ,   & ! PAN - 033
    0.         ,   & ! OP1 - 034
    0.         ,   & ! OP2 - 035
    -1500.     ,   & ! ORA1 - 036
    0.         ,   & ! ORA2 - 037
    0.         ,   & ! MO2 - 038
    0.         ,   & ! AKAP - 039
    0.         ,   & ! AKEP - 040
    0.         ,   & ! BIOP - 041
    0.         ,   & ! PHO - 042
    0.         ,   & ! ADD - 043
    0.         ,   & ! AROP - 044
    0.         ,   & ! CBOP - 045
    0.         ,   & ! OLN - 046
    0.         ,   & ! XO2 - 047
    0.         ,   & ! DMS - 048
    0.         ,   & ! NH3 - 049
    0.         ,   & ! NH3 - 049
    0.             & ! CFC - 050
    /)
CONTAINS
!---------------------------------------------------------------------------------------------------
   SUBROUTINE get_HenrysLawCts(name,c1,c2,c3,c4)  
     IMPLICIT NONE
     character(len=*), intent(in) :: name
     real, intent(out):: c1,c2,c3,c4
     integer :: l,found
     found = 0 
loop2: DO l = 1,nspecies_HL
        IF(TRIM(spc_name(l)) == TRIM(name)) then 
          c1  = hstar(l)
          c2  =   dhr(l)
          c3  =   ak0(l)
          c4  =   dak(l)   
          found = 1	  
          EXIT loop2
	ENDIF
       enddo loop2
       IF(found == 0) then 
          c1  = notfound
          c2  = notfound
          c3  = notfound
          c4  = notfound
       ENDIF
   END SUBROUTINE get_HenrysLawCts  
  
END MODULE Henrys_law_ConstantsMod
