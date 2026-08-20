SUBROUTINE AER_EQSAM4CLIM(xPP,xTT,xRH,xAD,PLP,PRP,xsPM,xaPM,xPMs,xPMT,xRHO,xVOL,xPH,xHp,xGF,&
                     & xWH2O,xYPa,xYPs,xYMa,xYMs,xYG,Dd,lEQskip,xNO3_2NA,xNO3_2CA,xNO3_2OTH,& 
                     & xAERPH, xCLOUDPH, xRAINPH, xAERWAT)
! 
 USE PARKIND1  ,ONLY : JPIM     ,JPRB, JPRD
 USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
!
! For the purposes of the CAMS2_35 Framework Agreement, the EQSAM4clim module, which is
! owned by ResearchConcept io GmbH, will not be treated as a Deliverable but will be
! treated as background  intellectual property and any improvements made to EQSAM4clim during the  
! Term of the Framework Agreement will be treated as improvement intellectual property. Each EQSAM4clim
! sub-routine in the IFS shall include a clear header identifying the ownership. All changes, made by
! the Contractor or any Sub-contractor to IFS code, including the interfaces between IFS routines and
! EQSAM4clim routines, shall be treated as Deliverables.
!____________________________________________________________________________________________
 IMPLICIT NONE
!____________________________________________________________________________________________
! WRITTEN BY SWEN METZGER, 2012-2022
!
! Copyright 2012-2021+. All rights reserved.
!
! Contact, w.r.t. EQSAM4clim modeling:
! Swen Metzger <sm@researchconcepts.io>
! ResearchConcepts io GmbH, HRB 717519
! https://www.researchconcepts.io
!
! Metzger, S., B. Steil, L. Xu, J. E. Penner, and J. Lelieveld, 
! New representation of water activity based on a single solute specific constant to 
! parameterize the hygroscopic growth of aerosols in atmospheric models, 
! Atmos. Chem. Phys., 12, 5429-5446, https://doi.org/10.5194/acp-12-5429-2012, 2012; 
!
! Metzger, S., B. Steil, M. Abdelkader, K. Klingmueller, L. Xu, J.E. Penner, C. Fountoukis, 
! A. Nenes, and J. Lelieveld; Aerosol Water Parameterization: A single parameter framework; 
! Atmos. Chem. Phys., 16, 7213-7237, https://doi.org/10.5194/acp-16-7213-2016, 2016; 
!
! S Metzger et al., EUMETSAT ITT 15/210839 Validation Report, 12/2016:
! "Comparison of Metop - PMAp Version 2 AOD Products using Model Data"
! https://www.eumetsat.int/PMAP
!
! Metzger, S., Abdelkader, M., Steil, B., and Klingmueller, K.: 
! Aerosol water parameterization: long-term evaluation and importance for climate studies, 
! Atmos. Chem. Phys., 18, 16747-16774, https://doi.org/10.5194/acp-18-16747-2018, 2018.
!
! Svetlana Tsyro and Swen Metzger, contribution to the EMEP Status Report 1/2019
! "Transboundary particulate matter, photo-oxidants, acidifying and eutrophying components"
! Joint MSC-W & CCC & CEIP Report, Chapter 9 (p133): EQSAM4clim 
! Report available from https://emep.int/publ/emep2019_publications.html (pdf 31 MB).
!
! B Koo, S Metzger, P Vennam, C Emery, G Wilson, G Yarwood, 
! "Comparing the ISORROPIA and EQSAM Aerosol Thermodynamic Options in CAMx", 
! International Technical Meeting on Air Pollution Modelling and its Application, 
! ITM 2018: Air Pollution Modeling and its Application XXVI pp 93-98, 2020.
! https://link.springer.com/chapter/10.1007/978-3-030-22055-6_16.
!
! - diagnostic output:
! xAD   = aerosol domain         [0-5]
! xsPM  = aerosol dry     mass   [umol/m3(air)]
! xaPM  = aerosol aqueous mass   [umol/m3(air)]
! xPMs  = aerosol dry     mass   [ ug/m3(air)]
! xPMt  = aerosol total   mass   [ ug/m3(air)]
! xRHO  = aerosol density        [  g/cm3    ]
! xVOL  = aerosol volume         [cm3/m3(air)]
! xHp   = aerosol H+             [mol/m3(air)]
! xPH   = aerosol pH             [-]
! xGF   = growth  factor         [-]
! - input/output:
! xWH2O = aerosol water   mass   [kg/m3(air)]
! xYPa  = cations (aqueous phase)[mol/m3(air)]
! xYPs  = cations (solid phase)  [mol/m3(air)]
!    (1)=NH4+, (2)=Na+,  (3)=K+,(4)=Ca++,(5)=Mg++ 
! xYMa  = anions  (aqueous phase)[mol/m3(air)]
! xYMs  = anions  (solid phase)  [mol/m3(air)]
!    (1)=SO4--,(2)=HSO4-,(3)=NO3-,(4)=Cl-
! xYG   = gases   (gaseous phase)[mol/m3(air)]
!    (1)=HNO3, (2)=HCl,  (3)=NH3, (4)=H2SO4 
! - input:
! xPP    = pressure              [Pa] - optional
! xTT    = temperature           [K]
! xRH    = relative humidity     [0-1]
! PLP    = cloud water content   [kg/m3(air)]
! PRP    = rain water content    [kg/m3(air)]
! Dd     = aerosol dry diameter  [m]
! lo     = array switch dimension [-] - optional
! imask  = array to skip values   [-] - optional (only for 3D)
!
! Last update:
! - September 2021 
!                  -- single / double precison mode
!                  -- H+ / pH parameterization revised
! - October 2021
!                  -- added I/O documentation
!                  -- some changes to comments
!                  -- removed imask array from input
!                  -- tested in single / double precison mode
!                  -- reduced dimensions to simplify code structure
!                  -> H+ / pH parameterization adition based on EQSAM2:
!                  -- aqueous uptake of NH3 and dissociation into NH4+
!                  -- added pressure to interface for 3-D applications (test only)
!                  -- added pre-existing water to the pH eqn for 3-D applications
! - December 2021
!                  -- added rain and cloud water content for pH calculation
!                  -- Coarse mode nitrate (SR)
! - February 2022
!                  -- added optional organic acid neutralization from EQSAM2/EQSAM3
!                     new switch: lorganics lumped organics (considering CH2O2)
!                     Formic acid, i.e.,    CH2O2 = H(COOH), and their salts:
!                     Ammonium  formate,  NH4CHO2 = H(COONH4)
!                     Sodium    formate,   NaCHO2 = H(COONa)
!                     Potassium formate,    KCHO2 = H(COOK)
!                     Calcium   formate,  Ca(CHO2)2
!                     Magnesium formate,  Mg(CHO2)2
! - March 2022
!                  -- added I/O and merged aqueous phase chemistry from TM5_WETCHEM
!                  -- testing size distribution dependency on mineral cations/speciation
! Notes:
!                  -> missing processes relevant for H+ (so far only in EQSAM2/EQSAM3):
!                  -> aqueous uptake of aerosol precursor gases (H2SO4/NH3/HNO3/HCl)
!                  -> dissociation of H2SO4/HSO4-/SO4--, NH3/NH4+, HNO3/NO3-, HCl/Cl-
!                  => to be merged with aqueous phase chemistry from TM5_WETCHEM
!____________________________________________________________________________________________
 CHARACTER(LEN=*),PARAMETER :: modstr  = 'EQSAM4clim'    ! name of module
 CHARACTER(LEN=*),PARAMETER :: modver  = 'v12'           ! module version
 CHARACTER(LEN=*),PARAMETER :: moddat  = '15Mar2022'     ! last modification date
!____________________________________________________________________________________________
! ' Hydrogen |   H2O           H2SO4          HNO3            HCl            NH3      HCOOH '
! '    H+    |    1              2             3               4              5         31  '
! '  Index   |  eqh2o          eqhsa          eqhna          eqhca          eqxam     eqhof '
! '----------|------------------------------------------------------------------------------'
! ' Ammonium |(NH4)3H(SO4)2  (NH4)2SO4       NH4HSO4        NH4NO3          NH4Cl    NH4CHO2'
! '   NH4+   |    6              7             8               9             10         32  '
! '  Index   |  eqalc          eqasu          eqahs          eqano          eqacl     eqaof '
! '----------|------------------------------------------------------------------------------'
! ' Sodium   | Na3H(SO4)2     Na2SO4          NaHSO4         NaNO3          NaCl      NaCHO2'
! '   Na+    |   11             12            13              14             15         33  '
! '  Index   |  eqslc          eqssu          eqshs          eqsno          eqscl     eqsof '
! '----------|------------------------------------------------------------------------------'
! ' Potassium|  K3H(SO4)2      K2SO4           KHSO4          KNO3           KCl      KCHO2 '
! '    K+    |   16             17            18              19             20         34  '
! '  Index   |  eqplc          eqpsu          eqphs          eqpno          eqpcl     eqpof '
! '----------|------------------------------------------------------------------------------'
! ' Calcium  |   ---           CaSO4          ---           Ca(NO3)2        CaCl2  Ca(CHO2)2'
! '   Ca++   |   21             22            23              24             25         35  '
! '  Index   |  eqc01          eqcsu          eqc02          eqcno          eqccl     eqcof '
! '----------|------------------------------------------------------------------------------'
! ' Magnesium|   ---           MgSO4          ---           Mg(NO3)2        MgCl2  Mg(CHO2)2'
! '   Mg++   |   26             27            28              29             30         36  '
! '  Index   |  eqm01          eqmsu          eqm02          eqmno          eqmcl     empof '
!____________________________________________________________________________________________
! INTEGER,PARAMETER                            :: JPIM = SELECTED_INT_KIND(9)       ! SP
! INTEGER,PARAMETER                            :: JPRB = SELECTED_REAL_KIND(6,37)   ! SP
!INTEGER,PARAMETER                            :: JPRB = SELECTED_REAL_KIND(12,307) ! DP
!______________________________________________
!- CAMS2_35 ->
!REAL(KIND=JPRB)                               :: ZHOOK_HANDLE
 REAL(KIND=JPRB),INTENT(IN)                    :: Dd,PLP,PRP
 REAL(KIND=JPRB),DIMENSION(4),INTENT(INOUT)    :: xYG
 REAL(KIND=JPRB),DIMENSION(4),INTENT(INOUT)    :: xYMa,xYMs
 REAL(KIND=JPRB),DIMENSION(5),INTENT(INOUT)    :: xYPa,xYPs
 REAL(KIND=JPRB),INTENT(IN)                    :: xPP,xTT, xRH
 REAL(KIND=JPRB),INTENT(INOUT)                 :: xWH2O
 REAL(KIND=JPRB),INTENT(OUT)                   :: xsPM,xPMs,xaPM,xPMt,xPH,xVOL,xRHO,xGF,xHp,xAD
 REAL(KIND=JPRB),INTENT(OUT)                   :: xNO3_2NA,xNO3_2CA,xNO3_2OTH,xAERPH
 REAL(KIND=JPRB),INTENT(OUT)                   :: xCLOUDPH,xRAINPH, xAERWAT
 REAL(KIND=JPRB)                               :: eqPHaer, eqPHcloud, eqPHrain
 LOGICAL,INTENT(INOUT)                         :: leqskip
!- CAMS2_35 <-
!______________________________________________
REAL(KIND=JPRB),PARAMETER                            :: REALZERO=tiny(0._JPRB),ZERO=0._JPRB,ONE=1._JPRB
REAL(KIND=JPRB),PARAMETER                            :: TINYX=1.e-15_JPRB,eqT0=298.15_JPRB ![K]
REAL(KIND=JPRB),PARAMETER                            :: R=8.314409_JPRB            ![J/mol/K]
REAL(KIND=JPRB),PARAMETER                            :: eqR=R/101325._JPRB   ![atm*m^3/mol/K]
REAL(KIND=JPRB),PARAMETER                            :: sigma=0.0761_JPRB    ![J/m^2]
REAL(KIND=JPRB),PARAMETER                            :: RHMAX=0.98_JPRB      ![-]
REAL(KIND=JPRB),PARAMETER                            :: RHMIN=0.01_JPRB      ![-]
!______________________________________________
 LOGICAL,PARAMETER                            :: lke         =.FALSE. ! TRUE = Kelvin effect, else Ke=1
 LOGICAL,PARAMETER                            :: lvola       =.TRUE.  ! TRUE = semi-volatiles, else none
 LOGICAL,PARAMETER                            :: lmixs       =.FALSE. ! TRUE = mixed solution Keq, else pure
 LOGICAL,PARAMETER                            :: lrhdm       =.TRUE.  ! TRUE = mixed solution RHD, else RHD only
 LOGICAL,PARAMETER                            :: lHSO4       =.TRUE.  ! TRUE = SO4/HSO4 partioning, else SO4 only
 LOGICAL,PARAMETER                            :: lH2SO4gas   =.FALSE. ! TRUE = residual gaseous H2SO4, else aerosol
 LOGICAL,PARAMETER                            :: lmetastable =.TRUE. ! TRUE = aqueous phase only, else + solids 
 LOGICAL,PARAMETER                            :: lorganics   =.FALSE. ! TRUE= lumped organics (considering CH2O2)
 LOGICAL,PARAMETER                            :: laqchem     =.FALSE. ! TRUE= merged aqueous phase chem of TM5_WETCHEM
!______________________________________________
INTEGER(KIND=JPIM),PARAMETER :: jMs=1
INTEGER(KIND=JPIM),PARAMETER :: jDs=2
INTEGER(KIND=JPIM),PARAMETER :: jZa=3
INTEGER(KIND=JPIM),PARAMETER :: jWs=4
INTEGER(KIND=JPIM),PARAMETER :: jNs=5
INTEGER(KIND=JPIM),PARAMETER :: jNi=6
INTEGER(KIND=JPIM),PARAMETER :: jRHD=7
INTEGER(KIND=JPIM),PARAMETER :: jRHDc=8
INTEGER(KIND=JPIM),PARAMETER :: jSOL=30
INTEGER(KIND=JPIM),PARAMETER :: jVOL=2
INTEGER(KIND=JPIM),PARAMETER :: jGAS=4
INTEGER(KIND=JPIM),PARAMETER :: jKEQ=2
INTEGER(KIND=JPIM),PARAMETER :: jAP=1
INTEGER(KIND=JPIM),PARAMETER :: jDP=2
INTEGER(KIND=JPIM),PARAMETER :: jGP=3
INTEGER(KIND=JPIM),PARAMETER :: jSP=1
INTEGER(KIND=JPIM),PARAMETER :: jSM=2
INTEGER(KIND=JPIM),PARAMETER :: jSS=3
!______________________________________________
INTEGER(KIND=JPIM),PARAMETER :: mciam =  1  ! ammonium        -  NH4+
INTEGER(KIND=JPIM),PARAMETER :: mciso =  2  ! sodium          -  Na+
INTEGER(KIND=JPIM),PARAMETER :: mcipo =  3  ! potassium       -  K+
INTEGER(KIND=JPIM),PARAMETER :: mcica =  4  ! calcium         -  Ca++
INTEGER(KIND=JPIM),PARAMETER :: mcimg =  5  ! magnesium       -  Mg++
INTEGER(KIND=JPIM),PARAMETER :: mcati =  5
!______________________________________________
INTEGER(KIND=JPIM),PARAMETER :: mdumm =  0  ! dummy
INTEGER(KIND=JPIM),PARAMETER :: maisu =  1  ! sulfate         -  SO4--
INTEGER(KIND=JPIM),PARAMETER :: maihs =  2  ! bisulfate       -  HSO4-
INTEGER(KIND=JPIM),PARAMETER :: maino =  3  ! nitrate         -  NO3-
INTEGER(KIND=JPIM),PARAMETER :: maicl =  4  ! chloride        -  Cl-
INTEGER(KIND=JPIM),PARAMETER :: manio =  4
!_______________________________________________
INTEGER(KIND=JPIM)                      :: II,IJ,Ip,Im,Is
!______________________________________________
CHARACTER(len= 5),  DIMENSION(0:jSOL)   :: ceqsolute
INTEGER(KIND=JPIM), DIMENSION(0:jSOL,5) :: ieqsolute
INTEGER(KIND=JPIM), DIMENSION(1:jVOL)   :: ieqvola
INTEGER(KIND=JPIM), DIMENSION(1:jGAS)   :: ieqgases
!_______________________________________________
LOGICAL,            DIMENSION(0:jSOL)   :: leqsolute
!______________________________________________
REAL(KIND=JPRB)                         :: T0T,AKW,X1,X2,X3,XZ,KEQ,Mw,Dw,COEF
!______________________________________________
REAL(KIND=JPRB),DIMENSION(0:mcati)  :: eqxpz,eqxpm
REAL(KIND=JPRB),DIMENSION(0:manio)  :: eqxmz,eqxmm
REAL(KIND=JPRB),DIMENSION(0:jSOL)   :: eqxyZ,eqxyA,eqxyB,eqxyC
REAL(KIND=JPRB),DIMENSION(8,0:jSOL) :: eqxy
!______________________________________________
INTEGER(KIND=JPIM) :: NZeq,IWATeq,IDeq
!______________________________________________
INTEGER(KIND=JPIM),DIMENSION(0:jSOL) :: ieqsalt
!______________________________________________
REAL(KIND=JPRB) :: eqP,eqTT,eqRH,eqsPM,eqPMs,eqaPM,eqPMt,eqWH2O,eqPH,eqGF,YY,YZ
REAL(KIND=JPRB) :: eqVOL,eqRHO,eqXPi,eqXMi,eqHPLUS,XRHDMIN,XRHDMAX,XRHDIFF,TSO4
REAL(KIND=JPRB) :: eqMolm,eqMs,eqNs,eqWs,eqNZm,eqRHDM,aqWH2O,tWH2O,lWH2O!,eqR,Dd
!_______________________________________________
REAL(KIND=JPRB),DIMENSION(mcati,2)   :: eqyp
REAL(KIND=JPRB),DIMENSION(manio,2)   :: eqym
REAL(KIND=JPRB),DIMENSION(jGAS)      :: eqyg
REAL(KIND=JPRB),DIMENSION(0:jSOL,3)  :: eqys
REAL(KIND=JPRB),DIMENSION(0:jSOL)    :: eqxyRHD
REAL(KIND=JPRB),DIMENSION(0:jSOL)    :: eqxyXS,eqxyMol,eqxyKe,eqxyGF
REAL(KIND=JPRB),DIMENSION(0:mcati)   :: eqxp
REAL(KIND=JPRB),DIMENSION(0:manio)   :: eqxm
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!______________________________________________
INTEGER(KIND=JPIM),PARAMETER ::     eqdum =  0, &
        &  eqh2o =  1,    eqhsa =  2,    eqhna =  3,    eqhca =  4,    eqxam = 5,&
        &  eqalc =  6,    eqasu =  7,    eqahs =  8,    eqano =  9,    eqacl = 10,&
        &  eqslc = 11,    eqssu = 12,    eqshs = 13,    eqsno = 14,    eqscl = 15,&
        &  eqplc = 16,    eqpsu = 17,    eqphs = 18,    eqpno = 19,    eqpcl = 20,&
        &  eqc01 = 21,    eqcsu = 22,    eqc02 = 23,    eqcno = 24,    eqccl = 25,&
        &  eqm01 = 26,    eqmsu = 27,    eqm02 = 28,    eqmno = 29,    eqmcl =30
!____________________________________________________________________________________________
INTEGER(KIND=JPIM)                :: ieqsu(mcati)=(/eqasu,eqssu,eqpsu,eqcsu,eqmsu/)
INTEGER(KIND=JPIM)                :: ieqhs(mcati)=(/eqahs,eqshs,eqphs,mdumm,mdumm/)
INTEGER(KIND=JPIM)                :: ieqno(mcati)=(/eqano,eqsno,eqpno,eqcno,eqmno/)
INTEGER(KIND=JPIM)                :: ieqcl(mcati)=(/eqacl,eqscl,eqpcl,eqccl,eqmcl/)
!____________________________________________________________________________________________
INTEGER(KIND=JPIM)                :: ieqams(manio)=(/eqasu,eqahs,eqano,eqacl/)
INTEGER(KIND=JPIM)                :: ieqsos(manio)=(/eqssu,eqshs,eqsno,eqscl/)
INTEGER(KIND=JPIM)                :: ieqpos(manio)=(/eqpsu,eqphs,eqpno,eqpcl/)
INTEGER(KIND=JPIM)                :: ieqcas(manio)=(/eqcsu,mdumm,eqcno,eqccl/)
INTEGER(KIND=JPIM)                :: ieqmgs(manio)=(/eqmsu,mdumm,eqmno,eqmcl/)
!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('AER_EQSAM4CLIM',0,ZHOOK_HANDLE) ! CAMS2_35
!______________________________________________
!____________________________________________________________________________________________
!REAL(KIND=JPRB),DIMENSION(4)    :: xYG
!REAL(KIND=JPRB),DIMENSION(4)    :: xYMa,xYMs
!REAL(KIND=JPRB),DIMENSION(5)    :: xYPa,xYPs
!REAL(KIND=JPRB)                 :: xPP,xTT,xRH,xAD,xsPM,xPMs,xaPM,xPMt,xWH2O,xPH,xVOL,xRHO,xGF,xHp
!____________________________________________________________________________________________
! cation name     |            mciam      mciso      mcipo      mcica       mcimg
eqxpz(0:mcati)=(/ 1._JPRB,  1.0_JPRB,    1.0_JPRB,    1.0_JPRB,     2.0_JPRB,     2.0_JPRB     /)
eqxpm(0:mcati)=(/ 1._JPRB, 0.01805_JPRB, 0.02299_JPRB, 0.039098_JPRB, 0.040080_JPRB, 0.024305_JPRB /)
!____________________________________________________________________________________________
! anion  name     |            maisu      maihs      maino      maicl
eqxmz(0:manio)=(/ 1._JPRB,  2.0_JPRB,    1.0_JPRB,    1.0_JPRB,    1.0_JPRB /)
eqxmm(0:manio)=(/ 1._JPRB,0.09607_JPRB,0.09708_JPRB,0.062010_JPRB,0.03545_JPRB/)
!____________________________________________________________________________________________
!____________________________________________________________________________________________
ieqvola (1:jVOL)=(/eqacl,eqano/)
ieqgases(1:jGAS)=(/eqhna,eqhca,eqxam,eqhsa/)
!____________________________________________________________________________________________
leqsolute(0:jSOL)=(/.FALSE.,&
         .FALSE.,       .FALSE.,       .FALSE.,       .FALSE.,       .FALSE.,&
         .TRUE.,        .TRUE.,        .TRUE.,        .TRUE.,        .TRUE.,&
         .FALSE.,       .TRUE.,        .TRUE.,        .TRUE.,        .TRUE.,&
         .FALSE.,       .TRUE.,        .TRUE.,        .TRUE.,        .TRUE.,&
         .FALSE.,       .TRUE.,        .FALSE.,       .TRUE.,        .TRUE.,&
         .FALSE.,       .TRUE.,        .FALSE.,       .TRUE.,        .TRUE./)
!____________________________________________________________________________________________
ieqsolute(0:jSOL,jSS)=(/mdumm,&
          eqh2o,         eqhsa,         eqhna,         eqhca,         eqxam,&
          eqalc,         eqasu,         eqahs,         eqano,         eqacl,&
          eqslc,         eqssu,         eqshs,         eqsno,         eqscl,&
          eqplc,         eqpsu,         eqphs,         eqpno,         eqpcl,&
          eqc01,         eqcsu,         eqc02,         eqcno,         eqccl,&
          eqm01,         eqmsu,         eqm02,         eqmno,         eqmcl/)
!____________________________________________________________________________________________
ieqsolute(0:jSOL,jSP)=(/mdumm,&
          mdumm,         mdumm,         mdumm,         mdumm,         mciam,&
          mciam,         mciam,         mciam,         mciam,         mciam,&
          mciso,         mciso,         mciso,         mciso,         mciso,&
          mcipo,         mcipo,         mcipo,         mcipo,         mcipo,&
          mdumm,         mcica,         mcica,         mcica,         mcica,&
          mdumm,         mcimg,         mdumm,         mcimg,         mcimg/)
!____________________________________________________________________________________________
ieqsolute(0:jSOL,jSM)=(/mdumm,&
          mdumm,         maisu,         maino,         maicl,         mdumm,&
          maisu,         maisu,         maihs,         maino,         maicl,&
          maisu,         maisu,         maihs,         maino,         maicl,&
          maisu,         maisu,         maihs,         maino,         maicl,&
          mdumm,         maisu,         mdumm,         maino,         maicl,&
          mdumm,         maisu,         mdumm,         maino,         maicl/)
!____________________________________________________________________________________________
ceqsolute(0:jSOL)=(/ 'eqdum',                                             &
         'eqh2o',       'eqhsa',       'eqhna',       'eqhca',       'eqxam',&
         'eqalc',       'eqasu',       'eqahs',       'eqano',       'eqacl',&
         'eqslc',       'eqssu',       'eqshs',       'eqsno',       'eqscl',&
         'eqplc',       'eqpsu',       'eqphs',       'eqpno',       'eqpcl',&
         'eqc01',       'eqcsu',       'eqc02',       'eqcno',       'eqccl',&
         'eqm01',       'eqmsu',       'eqm02',       'eqmno',       'eqmcl'/)
!____________________________________________________________________________________________
!        Ms   [kg/mol],&
eqxy(jMs,0:jSOL)=(/  ZERO, &
       0.018020_JPRB,   0.098090_JPRB,   0.063020_JPRB,   0.036460_JPRB,   0.017040_JPRB,&
       0.247300_JPRB,   0.132170_JPRB,   0.115130_JPRB,   0.080060_JPRB,   0.053500_JPRB,&
       0.262120_JPRB,   0.142050_JPRB,   0.120070_JPRB,   0.085000_JPRB,   0.058440_JPRB,&
       0.310444_JPRB,   0.174266_JPRB,   0.136178_JPRB,   0.101108_JPRB,   0.074548_JPRB,&
       0.000000_JPRB,   0.136150_JPRB,   0.000000_JPRB,   0.164100_JPRB,   0.110980_JPRB,&
       0.000000_JPRB,   0.120375_JPRB,   0.000000_JPRB,   0.148325_JPRB,   0.095205_JPRB/)
!____________________________________________________________________________________________
!        Ds   [kg/m^3],&
eqxy(jDs,0:jSOL)=(/  ZERO, &
     997.000000_JPRB,1830.000000_JPRB,1513.000000_JPRB,1490.000000_JPRB, 696.000000_JPRB,&
    1775.000000_JPRB,1770.000000_JPRB,1780.000000_JPRB,1720.000000_JPRB,1519.000000_JPRB,&
    2565.000000_JPRB,2700.000000_JPRB,2430.000000_JPRB,2260.000000_JPRB,2170.000000_JPRB,&
    2490.000000_JPRB,2660.000000_JPRB,2320.000000_JPRB,2110.000000_JPRB,1988.000000_JPRB,&
    1700.000000_JPRB,2960.000000_JPRB,1700.000000_JPRB,2500.000000_JPRB,2150.000000_JPRB,&
    1000.000000_JPRB,2660.000000_JPRB,1000.000000_JPRB,2300.000000_JPRB,2325.000000_JPRB/)
!____________________________________________________________________________________________
!        Za        [-],&
eqxy(jZa,0:jSOL)=(/  ZERO, &
       1.000000_JPRB,   2.000000_JPRB,   1.000000_JPRB,   1.000000_JPRB,   1.000000_JPRB,&
       3.000000_JPRB,   2.000000_JPRB,   1.000000_JPRB,   1.000000_JPRB,   1.000000_JPRB,&
       3.000000_JPRB,   2.000000_JPRB,   1.000000_JPRB,   1.000000_JPRB,   1.000000_JPRB,&
       3.000000_JPRB,   2.000000_JPRB,   1.000000_JPRB,   1.000000_JPRB,   1.000000_JPRB,&
       3.000000_JPRB,   2.000000_JPRB,   3.000000_JPRB,   2.000000_JPRB,   2.000000_JPRB,&
       1.000000_JPRB,   2.000000_JPRB,   1.000000_JPRB,   2.000000_JPRB,   2.000000_JPRB/)
!____________________________________________________________________________________________
!        Ws (T_o)  [%],&
eqxy(jWs,0:jSOL)=(/  ZERO, &
       0.000000_JPRB,  70.000000_JPRB,  25.000000_JPRB,  15.000000_JPRB,  30.000000_JPRB,&
      53.300000_JPRB,  43.310000_JPRB,  76.000000_JPRB,  68.050000_JPRB,  28.340000_JPRB,&
      44.060000_JPRB,  21.940000_JPRB,  66.180000_JPRB,  47.700000_JPRB,  26.470000_JPRB,&
       0.000000_JPRB,  10.710000_JPRB,  33.600000_JPRB,  27.690000_JPRB,  26.230000_JPRB,&
       0.000000_JPRB,   0.210000_JPRB,   0.000000_JPRB,  59.020000_JPRB,  44.840000_JPRB,&
       0.000000_JPRB,  26.310000_JPRB,   0.000000_JPRB,  41.590000_JPRB,  35.900000_JPRB/)
!____________________________________________________________________________________________
!        nu_s     [-],&
eqxy(jNs,0:jSOL)=(/  ZERO, &
       2.000000_JPRB,   3.000000_JPRB,   2.000000_JPRB,   2.000000_JPRB,   1.000000_JPRB,&
       5.000000_JPRB,   3.000000_JPRB,   2.000000_JPRB,   2.000000_JPRB,   2.000000_JPRB,&
       5.000000_JPRB,   3.000000_JPRB,   2.000000_JPRB,   2.000000_JPRB,   2.000000_JPRB,&
       5.000000_JPRB,   3.000000_JPRB,   2.000000_JPRB,   2.000000_JPRB,   2.000000_JPRB,&
       5.000000_JPRB,   2.000000_JPRB,   5.000000_JPRB,   3.000000_JPRB,   3.000000_JPRB,&
       1.000000_JPRB,   2.000000_JPRB,   1.000000_JPRB,   3.000000_JPRB,   3.000000_JPRB/)
!____________________________________________________________________________________________
!        nu_i     [-],&
eqxy(jNi,0:jSOL)=(/  ZERO, &
       1.000000_JPRB,   1.761309_JPRB,   1.793689_JPRB,   2.681333_JPRB,   0.527416_JPRB,&
       1.616356_JPRB,   1.274822_JPRB,   1.253573_JPRB,   1.051480_JPRB,   1.243054_JPRB,&
       1.000000_JPRB,   1.278762_JPRB,   1.293906_JPRB,   1.160345_JPRB,   1.358377_JPRB,&
       1.000000_JPRB,   1.286445_JPRB,   1.308499_JPRB,   1.014102_JPRB,   1.256989_JPRB,&
       1.000000_JPRB,   1.271828_JPRB,   1.000000_JPRB,   1.586562_JPRB,   2.024869_JPRB,&
       1.000000_JPRB,   1.435281_JPRB,   1.000000_JPRB,   1.878693_JPRB,   2.107772_JPRB/)
!____________________________________________________________________________________________
!        RHD (T_o) [-],&
eqxy(jRHD,0:jSOL)=(/  ZERO, &
       1.000000_JPRB,   1.000000_JPRB,   1.000000_JPRB,   1.000000_JPRB,   1.000000_JPRB,&
       0.690000_JPRB,   0.799700_JPRB,   0.400000_JPRB,   0.618300_JPRB,   0.771000_JPRB,&
       0.390000_JPRB,   0.930000_JPRB,   0.520000_JPRB,   0.737900_JPRB,   0.752800_JPRB,&
       1.000000_JPRB,   0.975000_JPRB,   0.860000_JPRB,   0.924800_JPRB,   0.842600_JPRB,&
       1.000000_JPRB,   0.990000_JPRB,   1.000000_JPRB,   0.490600_JPRB,   0.283000_JPRB,&
       1.000000_JPRB,   0.861300_JPRB,   1.000000_JPRB,   0.540000_JPRB,   0.328400_JPRB/)
!____________________________________________________________________________________________
!        T-coef.   [-],&
eqxy(jRHDc,0:jSOL)=(/  ZERO, &
       0.000000_JPRB,   0.000000_JPRB,   0.000000_JPRB,   0.000000_JPRB,   0.000000_JPRB,&
     186.000000_JPRB,  80.000000_JPRB, 384.000000_JPRB, 852.000000_JPRB, 239.000000_JPRB,&
       0.000000_JPRB,  80.000000_JPRB, -45.000000_JPRB, 304.000000_JPRB,  25.000000_JPRB,&
       0.000000_JPRB,  35.600000_JPRB,   0.000000_JPRB,   0.000000_JPRB, 159.000000_JPRB,&
       0.000000_JPRB,   0.000000_JPRB,   0.000000_JPRB, 509.400000_JPRB, 551.100000_JPRB,&
       0.000000_JPRB,-714.450000_JPRB,   0.000000_JPRB, 230.200000_JPRB,  42.230000_JPRB/)
!____________________________________________________________________________________________
! GLOBAL INITIALIZATION
!______________________________________
!    lke         = leqm(1)
!    lvola       = leqm(2)
!    lmixs       = leqm(3)
!    lrhdm       = leqm(4)
!    lHSO4       = leqm(5)
!    lH2SO4gas   = leqm(6)
!    lmetastable = leqm(7)
!    laqchem     = leqm(8)
!    lorganics   = leqm(9)
!   WRITE(*,*) ' lke=',lke,' lvola=',lvola,' lmixs=',lmixs,' lrhdm=',lrhdm,&
!              ' lHSO4=',lHSO4,' lH2SO4gas=',lH2SO4gas,' lmetastable=',lmetastable,&
!              ' laqchem=',laqchem,' lorganics=',lorganics
!______________________________________
      xHp=ZERO
      xPH=ZERO
      eqXPi=ZERO
      eqXMi=ZERO
      eqxp (0)=ZERO
      eqxp (1)=ZERO
      eqxp (2)=ZERO
      eqxp (3)=ZERO
      eqxp (4)=ZERO
      eqxp (5)=ZERO
      eqxm (0)=ZERO
      eqxm (1)=ZERO
      eqxm (2)=ZERO
      eqxm (3)=ZERO
      eqxm (4)=ZERO
      eqyg (1)=ZERO
      eqyg (2)=ZERO
      eqyg (3)=ZERO
      eqyg (4)=ZERO
      lEQskip=.FALSE.
!______________________________________
! INPUT
   ! Pressure [Pa]
      eqP=xPP
   !  eqR=R/eqP               !   [m^3/mol/K] - test sensitivity
   ! Aerosol water   mass     [ ug/m3(air)]
!     aqWH2O=xWH2O*1.E-9_JPRB ! [ ug/m3(air)]
      aqWH2O=xWH2O            ! [ kg/m3(air)]
   ! Temperature            [K]
      eqTT  =xTT  
   ! Relative humidity      [0-1]
      eqRH  =MAX(RHMIN,MIN(xRH,RHMAX))  
   ! CATIONS [MOL/M^3(air)]
   DO IJ=1,mcati
      eqxp (IJ)    =xYPa(IJ)
      eqXPi=eqXPi+eqxp(IJ)*eqxpz(IJ)
   END DO ! mcati
   ! ANIONS [MOL/M^3(air)]
   DO IJ=1,manio
      eqxm (IJ)    =xYMa(IJ)
      eqXMi=eqXMi+eqxm(IJ)*eqxmz(IJ)
   END DO ! manio
   IF(eqXPi+eqXMi < TINYX) THEN
      lEQskip=.TRUE.
      IF (LHOOK) CALL DR_HOOK('AER_EQSAM4CLIM',1,ZHOOK_HANDLE) ! CAMS2_35
      RETURN
   END IF
!______________________________________
! DOMAINS
      IWATeq=2 ! SOLID+LIQUID AEROSOL
      IF(lmetastable.OR.aqWH2O > REALZERO) &
      IWATeq=1 ! METASTABLE AEROSOL
      TSO4 = eqxm(maihs)+eqxm(maisu)
      X1 = eqxmz (maihs)*eqxm(maihs)
      X2 = eqxmz (maisu)
      IF(lHSO4.AND.eqXPi>ZERO.AND.eqXPi<X1+X2) THEN
         Is=eqasu
         IF(eqxp(mciam) >= eqxm(maisu)*0.5_JPRB) Is=eqahs
         IF(eqxp(mcipo) >= eqxm(maisu)*0.5_JPRB) Is=eqphs
         IF(eqxp(mciso) >= eqxm(maisu)*0.5_JPRB) Is=eqshs
         IF(eqRH >= eqxy(jRHD,Is)) X2=1.95_JPRB
      END IF
      X2 =            X2*eqxm(maisu)
      X3 = eqxpz (mciam)*eqxp(mciam)
      IDeq = 0 ! SULFATE POOR
      IF(eqXPi < TINYX .AND. TSO4 > TINYX) THEN
         IDeq = 5 ! ONLY SULFURIC ACID
         eqxm(maisu)=eqxm(maihs)+eqxm(maisu)
         eqxm(maihs)=ZERO
      ELSE IF(eqXPi > TINYX .AND. eqXPi < TSO4) THEN
         IDeq = 4 ! SULFATE VERY RICH
         eqxm(maihs)=eqxm(maihs)+eqxm(maisu)
         eqxm(maisu)=ZERO
      ELSE IF(eqXPi >= TSO4 .AND. eqXPi < X1+X2) THEN
         IDeq = 3 ! SULFATE RICH
         eqxm(maisu)=eqxm(maihs)+eqxm(maisu)
         eqxm(maihs)=ZERO
         XZ = MAX(ZERO,X1+X2-eqXPi)
         eqxm(maihs)=eqxm(maihs)+XZ
         eqxm(maisu)=eqxm(maisu)-XZ
      ELSE IF(eqXPi >= X1+X2 .AND. eqXPi-X3 < TSO4) THEN
         IDeq = 2 ! SULFATE NEUTRAL
      ELSE IF(eqXPi-X3 >= TSO4) THEN
         IDeq = 1 ! SULFATE POOR / MINERAL CATION RICH
      END IF
      IF(IDeq <= 2) THEN
         eqxm(maisu)=eqxm(maisu)+eqxm(maihs)
         eqxm(maihs)=ZERO
      END IF
!______________________________________
! LOCAL INITIALIZATION
   Mw=eqxy(jMs,eqh2o)
   Dw=eqxy(jDs,eqh2o)
   YY     =ZERO
   eqWH2O =ZERO
   eqMolm =ZERO
   eqMs   =ZERO
   eqNs   =ZERO
   eqWs   =ZERO
   eqNZm  =ZERO
   NZeq   =0
   eqHPLUS=ZERO
   eqsPM  =ZERO
   eqaPM  =ZERO
   eqPMs  =ZERO
   eqPMt  =ZERO
   eqVOL  =ZERO
   eqGF   =ONE
   eqRHO  =ONE
   eqRHDM =ONE
   eqXPi  =ZERO
   eqXMi  =ZERO
   XRHDMAX=ZERO
   XRHDMIN=ONE
   XRHDIFF=ONE
   DO IJ=1,mcati
      eqyp(IJ,1) = ZERO
      eqyp(IJ,2) = ZERO
   END DO
   DO IJ=1,manio
      eqym(IJ,1) = ZERO
      eqym(IJ,2) = ZERO
   END DO
   DO IS=0,jSOL
      eqxyA(Is)=ONE
      eqxyB(Is)=ZERO
      eqxyC(Is)=ZERO
      eqxyZ(Is)=ZERO
      ieqsalt(Is)=0
      eqxyGF (Is)=ONE
      eqxyRHD(Is)=ONE
      eqxyKe  (Is)=ONE
      eqxyXS  (Is)=ONE
      eqxyMol (Is)=ZERO
      eqys(Is,jAP)=ZERO
      eqys(Is,jDP)=ZERO
      eqys(Is,jGP)=ZERO
   END DO
   eqPH     = -999._JPRB
   eqPHaer  = -999._JPRB
   eqPHcloud= -999._JPRB
   eqPHrain = -999._JPRB
!______________________________________
! NEUTRALIZATION/REACTION ORDER
      IF(IDeq < 3) THEN
          ieqsalt( 1)=eqcsu ; ieqsalt( 2)=eqmsu ; ieqsalt( 3)=eqpsu
          ieqsalt( 4)=eqssu ; ieqsalt( 5)=eqasu ; ieqsalt( 6)=eqcno
          ieqsalt( 7)=eqmno ; ieqsalt( 8)=eqpno ; ieqsalt( 9)=eqsno
          ieqsalt(10)=eqano ; ieqsalt(11)=eqccl ; ieqsalt(12)=eqmcl
          ieqsalt(13)=eqpcl ; ieqsalt(14)=eqscl ; ieqsalt(15)=eqacl
      END IF
      IF(IDeq == 3) THEN
          ieqsalt( 1)=eqcsu ; ieqsalt( 2)=eqmsu
          ieqsalt( 3)=eqpsu ; ieqsalt( 4)=eqphs
          ieqsalt( 5)=eqssu ; ieqsalt( 6)=eqshs
          ieqsalt( 7)=eqasu ; ieqsalt( 8)=eqahs
      END IF
      IF(IDeq == 4) THEN
          ieqsalt( 1)=eqcsu ; ieqsalt( 2)=eqmsu
          ieqsalt( 3)=eqphs ; ieqsalt( 4)=eqshs
          ieqsalt( 5)=eqahs ; ieqsalt( 6)=eqhsa
          ieqsalt( 7)=eqalc
      END IF
      IF(IDeq == 5) THEN
          ieqsalt( 1)=eqhsa
          ieqsalt( 2)=eqalc
      END IF
!______________________________________
! SOLUTE MOLALITY
   DO IJ=1,jSOL
         Is=ieqsalt(IJ)
         IF(.NOT.leqsolute(Is)) CYCLE
         eqxyB(Is)=ZERO
         eqxyZ(Is)=eqxy(jNi,Is)
         eqxyMol(Is)=((eqxyKe(Is)/eqRH-eqxyA(Is))/MW/eqxyZ(Is))**(ONE/eqxyZ(Is))
         !1st
         eqxyGF(Is) = (eqxy(jDs,Is)/(eqxy(jMs,Is)*Dw*eqxyMol(Is))+1._JPRB)**(1._JPRB/3._JPRB)
         IF(lke)    eqxyKe(Is) = exp(4._JPRB*Mw*sigma/(R*eqTT*Dw*eqxyGF(Is)*Dd))
         eqxyXS(Is) = ONE/(ONE/(eqxy(jMs,Is)*eqxyMol(Is))+ONE)
         eqxyB   (Is) = eqxyXS(Is)**(ONE/(ONE+eqxyZ(Is)+eqxyXS(Is)))
         eqxyC   (Is) = ((eqxyKe(Is)/eqRH-eqxyA(Is))/MW/eqxyZ(Is))**(ONE/eqxyZ(Is))
         IF(eqxyB(Is) < eqxyC(Is)) eqxyMol(Is) = eqxyC(Is)-eqxyB(Is)
         !2nd
         eqxyGF(Is) = (eqxy(jDs,Is)/(eqxy(jMs,Is)*Dw*eqxyMol(Is))+1._JPRB)**(1._JPRB/3._JPRB)
         IF(lke)    eqxyKe(Is) = exp(4._JPRB*Mw*sigma/(R*eqTT*Dw*eqxyGF(Is)*Dd))
         eqxyXS(Is) = ONE/(ONE/(eqxy(jMs,Is)*eqxyMol(Is))+ONE)
         eqxyB   (Is) = eqxyXS(Is)**(ONE/(ONE+eqxyZ(Is)+eqxyXS(Is)))
         eqxyC   (Is) = ((eqxyKe(Is)/eqRH-eqxyA(Is))/MW/eqxyZ(Is))**(ONE/eqxyZ(Is))
         IF(eqxyB(Is) < eqxyC(Is)) eqxyMol(Is) = eqxyC(Is)-eqxyB(Is)
         !3rd
         eqxyGF(Is) = (eqxy(jDs,Is)/(eqxy(jMs,Is)*Dw*eqxyMol(Is))+1._JPRB)**(1._JPRB/3._JPRB)
         IF(lke)    eqxyKe(Is) = exp(4._JPRB*Mw*sigma/(R*eqTT*Dw*eqxyGF(Is)*Dd))
         eqxyXS(Is) = ONE/(ONE/(eqxy(jMs,Is)*eqxyMol(Is))+ONE)
         eqxyB   (Is) = eqxyXS(Is)**(ONE/(ONE+eqxyZ(Is)+eqxyXS(Is)))
         eqxyC   (Is) = ((eqxyKe(Is)/eqRH-eqxyA(Is))/MW/eqxyZ(Is))**(ONE/eqxyZ(Is))
         IF(eqxyB(Is) < eqxyC(Is)) eqxyMol(Is) = eqxyC(Is)-eqxyB(Is)
         eqxyXS(Is) = ONE/(ONE/(eqxy(jMs,Is)*eqxyMol(Is))+ONE)
      END DO
!______________________________________
! RELATIVE HUMIDITY OF DELIQUESCENCE
   DO IJ=1,jSOL
         Is=ieqsalt(IJ)
         IF(IWATeq == 2) THEN
            IF(.NOT.leqsolute(Is)) CYCLE
            IF(eqxy(jRHD,Is) >= ONE) CYCLE
            eqxyRHD (Is)=eqxy(jRHD,Is)*EXP(eqxy(jRHDc,Is)*(ONE/eqTT-ONE/eqT0))
            eqxyRHD (Is)=eqxyRHD(Is)*eqxyKe(Is)
            eqxyRHD (Is)=MAX(ZERO,MIN(eqxyRHD(Is),0.999_JPRB))
         ELSE
            eqxyRHD (Is)=eqRH
         END IF
   END DO
!______________________________________
! EQUILIBRIUM
   DO IJ=1,jSOL
         Is=ieqsalt(IJ)
         IF(.NOT.leqsolute(Is)) CYCLE
         Ip=ieqsolute(Is,jSP)
         Im=ieqsolute(Is,jSM)
         IF(Ip == 0 .OR. Im == 0 .OR. IDeq==5) CYCLE
         IF(eqxp(Ip)*eqxm(Im) > REALZERO) THEN
            NZeq = NZeq + 1
            XZ = MAX(ZERO, MIN(eqxpz(Ip)*eqxp(Ip),eqxmz(Im )*eqxm(Im)))
            eqxp(Ip)   = MAX(ZERO,eqxp(Ip)-XZ/eqxpz(Ip))
            eqxm(Im)   = MAX(ZERO,eqxm(Im)-XZ/eqxmz(Im))
            eqys(Is,jAP) = eqys(Is,jAP)+XZ/eqxy(jZa,Is)
            eqNs=eqNs+eqys(Is,jAP)
         ELSE
            eqxyRHD(Is) = ONE
         END IF
   END DO
!______________________________________
   IF(lrhdm) THEN
      DO IJ=1,jSOL
         Is=ieqsalt(IJ)
         IF(.NOT.leqsolute(Is)) CYCLE
         IF(eqys(Is,jAP) > TINYX .AND. IWATeq == 2) THEN
            eqNZm =eqNZm+ONE
            eqMs  =eqMs+eqxy(jMs,Is)
            XZ       =(ONE/(100._JPRB/eqxy(jWs,Is)-ONE))/eqxy(jMs,Is)
            eqMolm=eqMolm+XZ
            eqWs  =MAX(0.1_JPRB,MIN(ONE,ONE/(ONE/(eqMolm*eqMs)+ONE)))
            YY    =ONE/MAX(ZERO,(0.25_JPRB*log(eqWs)+ONE))
            eqRHDM=ONE/(ONE+Mw*YY*(eqMolm)**YY)
         END IF
      END DO
   END IF
!______________________________________
! GAS/AEROSOL PARTITIONING OF SEMI-VOLATILES
   IF(lvola) THEN
   ! semi-volatiles
   DO IJ=1,jVOL
      Is=ieqvola(IJ)
      IF(.NOT.leqsolute(Is)) CYCLE
      Ip=ieqsolute(Is,jSP)
      Im=ieqsolute(Is,jSM)
         IF(IDeq >2) CYCLE
         XZ = eqys(Is,jAP)
         IF(XZ > TINYX) THEN
            YY=ONE
            ! TEMPERATURE DEPENDENT EQUILIBRIUM CONSTANTS
            T0T=eqT0/eqTT
            COEF=ONE+LOG(T0T)-T0T
            X1=ONE/(ONE/(eqxy(jMs,Is)*eqxyMol(Is))+ONE)
            X2=MAX(ZERO,MIN(2._JPRB*X1**2._JPRB,2._JPRB))
            XRHDMIN=eqxyRHD(Is)
            IF(Is == eqacl) THEN
               ! NH4CL(S) <==> NH3(G) + HCL(G)   [ppb^2]
               X3   = 1.086E-16_JPRB ! [ppb^2]
               X3   = X3*EXP(-71.00_JPRB*(T0T-ONE)+2.400_JPRB*COEF)
               KEQ  = X3/(eqR*eqTT)/(eqR*eqTT) ! [(mol^2/m^3(air))^2]
            END IF
            IF(Is == eqano) THEN
               ! NH4NO3(S) <==> NH3(G) + HNO3(G)
               ! ISORROPIA2
               X3   = 5.746E-17_JPRB ! [ppb^2]
               X3   = X3*EXP(-74.38_JPRB*(T0T-ONE)+6.120_JPRB*COEF)
               ! Mozurkewich (1993)
               !X3   = 4.199E-17_JPRB ! [ppb^2]
               !X3   = X3*EXP(-74.7351_JPRB*(T0T-ONE)+6.025_JPRB*COEF)
               ! SEQUILIB
               !X3   = 2.985e-17_JPRB ! [ppb^2]
               !X3   = X3*EXP(-75.11_JPRB*(T0T-ONE)+13.460_JPRB*COEF)
               KEQ   = X3/(eqR*eqTT)/(eqR*eqTT) ! [(mol^2/m^3(air))^2]
            END IF
            COEF=X2
            IF(lmixs.AND.TSO4>eqXPi) THEN
!           IF(lmixs.AND.TSO4>TINYX) THEN
               YY=(XZ/(XZ+3._JPRB*TSO4))**0.8_JPRB
               IF(lrhdm.AND.IWATeq==2) &
               XRHDMIN=eqRHDM*YY**0.25_JPRB+eqxyRHD(Is)   *(ONE-YY**0.25_JPRB)
!              XRHDMIN=eqRHDM*YY**0.25_JPRB+eqxyRHD(eqasu)*(ONE-YY**0.25_JPRB)
               IF(eqRH>=XRHDMIN) THEN
                  COEF=X2*YY
                  XZ=XZ*(ONE-COEF)
               END IF
            END IF
            IF(eqRH<XRHDMIN) THEN
               COEF=ONE
            ELSE
               IF(Is==eqacl) KEQ=KEQ*6.0_JPRB
            END IF
            KEQ=KEQ*COEF
            X1=eqxp(Ip)+eqxm(Im) ! [mol/m3(air)]
            X2=SQRT(X1*X1+4._JPRB*KEQ)
            X3=0.5_JPRB*(-X1+X2)
            X3=MIN(XZ,X3)
            eqxp(Ip)=eqxp(Ip)+X3
            eqxm(Im)=eqxm(Im)+X3
            eqys(Is,jAP)=MAX(0._JPRB,eqys(Is,jAP)-X3)
         END IF
   END DO
   END IF
!______________________________________
! LIQUID/SOLID PARTITIONING
   DO IJ=1,jSOL
         Is=ieqsalt(IJ)
         IF(.NOT.leqsolute(Is)) CYCLE
         IF(eqys(Is,jAP) < REALZERO) CYCLE
         IF(IWATeq == 2) THEN
            XRHDIFF=ONE
            XRHDMAX=eqxyRHD(Is)
            IF(lrhdm.AND.eqNZm>ONE.AND.Is/=eqcsu.AND.Is/=eqpsu) THEN
               XRHDMIN=eqRHDM
               YZ=eqys(Is,jAP)/eqNs
               XRHDMAX=XRHDMIN*YZ**0.25_JPRB+XRHDMAX*(ONE-YZ**0.25_JPRB)
            ELSE
               XRHDMIN=XRHDMAX
            END IF
            IF(eqRH < XRHDMAX) THEN
               IF(eqRH> XRHDMIN   .AND.   XRHDMIN<XRHDMAX) &
               XRHDIFF=(XRHDMAX-eqRH)/(XRHDMAX-XRHDMIN)
               eqys(Is,jDP) = MAX(0._JPRB,eqys(Is,jDP) + eqys(Is,jAP)*XRHDIFF)
               eqys(Is,jAP) = MAX(0._JPRB,eqys(Is,jAP) * (ONE-XRHDIFF))
            END IF
         END IF
   END DO
!______________________________________
! AEROSOL WATER (AW) [KG/M^3(AIR)]
   DO IJ=1,jSOL
         Is=ieqsalt(IJ)
         IF(.NOT.leqsolute(Is)) CYCLE
         IF(eqxyMol(Is) > REALZERO) &
         eqWH2O = eqWH2O + eqys(Is,jAP)/eqxyMol(Is)
   END DO
!______________________________________
! Remaining H+/OH- [MOL]
   DO IJ=1,mcati
         eqXPi=eqXPi+eqxp(IJ)*eqxpz(IJ)
!        WRITE(*,*) 'eqXPi= ',IJ,eqXPi,eqxp(IJ)
   END DO
   DO IJ=1,manio
         eqXMi=eqXMi+eqxm(IJ)*eqxmz(IJ)
!        WRITE(*,*) 'eqXMi= ',IJ,eqXMi,eqxm(IJ)
   END DO
!______________________________________
! RESIDUAL GASES
   DO IJ=1,jGAS
      Is=ieqgases(IJ)
      Ip=ieqsolute(Is,jSP)
      Im=ieqsolute(Is,jSM)
         XZ=eqxm(Im)
         IF(Is == eqxam) XZ=eqxp(Ip)
         IF(Is == eqhsa) XZ=eqxm(maisu)+eqxm(maihs)
         IF(XZ < REALZERO) CYCLE
         IF(Is == eqhna) THEN
            eqyg(1)  = XZ
            eqxm(Im) = ZERO
         ELSE IF(Is == eqhca) THEN
            eqyg(2)  = XZ
            eqxm(Im) = ZERO
         ELSE IF(Is == eqxam) THEN
            eqyg(3)  = XZ
            eqxp(Ip) = ZERO
         ELSE IF(Is == eqhsa) THEN
            II=eqalc
            IF(eqxyMol(II) > REALZERO) &
            eqWH2O = eqWH2O + XZ/eqxyMol(II)
            IF(lH2SO4gas) THEN
               eqyg(4)  = XZ
               eqxm(maisu) = ZERO
               eqxm(maihs) = ZERO
            END IF
         END IF
   END DO
!______________________________________
!  add previous AW to current time step
   tWH2O=aqWH2O & ! pre-existing AW
      & +eqWH2O & ! aerosol water
      & +PLP    & ! cloud   water
      & +PRP      ! rain    water
!tWH2O=1._JPRB ! test / remove
!tWH2O=1.E-6_JPRB! test / remove
!______________________________________
      ! AQUEOUS PHASE PROPERTIES
      eqGF = ONE
       IF(tWH2O > 0.1e-12_JPRB) THEN
         T0T=eqT0/eqTT
         COEF=ONE+LOG(T0T)-T0T
         KEQ  = 1.010E-14_JPRB*EXP(-22.52_JPRB*(T0T-ONE) + 26.920_JPRB*COEF)
         AKW  = KEQ**0.5_JPRB
         YY   = AKW/Mw/(1._JPRB-eqRH**2._JPRB)
         XZ   = YY-eqXPi+eqXMi
         IF(XZ < ZERO) THEN
            eqHPLUS = tWH2O/10._JPRB**(7._JPRB+log(-XZ/1.E-4_JPRB))
         ELSE IF(XZ == ZERO) THEN
            eqHPLUS = YY*1.E-6_JPRB
         ELSE
            eqHPLUS = XZ*1.E-6_JPRB
         END IF
!______________________________________
         eqPH  = -LOG10(eqHPLUS/tWH2O)    ! total aerosol pH
         lWH2O = (aqWH2O+eqWH2O)
         IF(lWH2O > 0.1e-12_JPRB) &
         eqPHaer = -LOG10(eqHPLUS/lWH2O)  ! aerosol water pH
          IF(PLP > 0.1e-12_JPRB) &
         eqPHcloud = -LOG10(eqHPLUS/PLP)  ! cloud water pH 
         IF(PRP > 0.1e-12_JPRB) &
         eqPHrain  = -LOG10(eqHPLUS/PRP)  ! rain water pH 
!______________________________________
!         IF(laqchem) THEN
         ! Test: needs to be merged with TM5_WETCHEM 
!           call TM5_WETCHEM2!()
!          ! EQUILIBRIUM CONSTANT NH3(g) <==> NH3(aq) [mol/kg atm]
!          X1 = 57.639_JPRB*EXP(13.79_JPRB*(T0T-ONE)-5.393_JPRB*COEF)
!          ! DISS. EQUILIBRIUM NH3(aq) <==> NH4(aq) + OH(aq) [mol/kg]
!          X2  = 1.805e-5_JPRB*EXP(-1.50_JPRB*(T0T-ONE)+26.920_JPRB*COEF)
!          XZ = X1*X2/AKW/tWH2O
!          X1 = XZ + eqHPLUS
!          X2 = XZ + eqyg(3)
!          X3 = ROOTF(X1,-X2)
!          XZ = MAX(ZERO,MIN(X3,eqyg(3)))
!          eqHPLUS = eqHPLUS - XZ
!          IF (eqHPLUS < ZERO) eqHPLUS = 10._JPRB**(7._JPRB+log(-eqHPLUS))
!          eqyg(3) = eqyg(3) - XZ
!          eqxp(1) = eqxp(1) + XZ
! !        WRITE(*,*) 'NH3= ',KEQ,XZ,eqyg(3),eqXPi,X3
!         END IF
!______________________________________
    ENDIF
!______________________________________
! OUTPUT
   DO IJ=1,mcati
        eqyp (IJ,1)=eqyp(IJ,1)+eqxp(IJ)
        eqaPM=eqaPM+eqxp(IJ)
        eqPMt=eqPMt+eqxp(IJ)*eqxpm(IJ)
   END DO
   DO IJ=1,manio
        eqym (IJ,1)=eqym(IJ,1)+eqxm(IJ)
        eqaPM=eqaPM+eqxm(IJ)
        eqPMt=eqPMt+eqxm(IJ)*eqxmm(IJ)
   END DO
   DO IJ=1,mcati
         eqym(maisu,1) = eqym(maisu,1) + eqys(ieqsu(IJ),jAP)*eqxy(jZa,ieqsu(IJ))/eqxmz(ieqsolute(ieqsu(IJ),jSM))
         eqym(maisu,2) = eqym(maisu,2) + eqys(ieqsu(IJ),jDP)*eqxy(jZa,ieqsu(IJ))/eqxmz(ieqsolute(ieqsu(IJ),jSM))
         eqym(maihs,1) = eqym(maihs,1) + eqys(ieqhs(IJ),jAP)*eqxy(jZa,ieqhs(IJ))/eqxmz(ieqsolute(ieqhs(IJ),jSM))
         eqym(maihs,2) = eqym(maihs,2) + eqys(ieqhs(IJ),jDP)*eqxy(jZa,ieqhs(IJ))/eqxmz(ieqsolute(ieqhs(IJ),jSM))
         eqym(maino,1) = eqym(maino,1) + eqys(ieqno(IJ),jAP)*eqxy(jZa,ieqno(IJ))/eqxmz(ieqsolute(ieqno(IJ),jSM))
         eqym(maino,2) = eqym(maino,2) + eqys(ieqno(IJ),jDP)*eqxy(jZa,ieqno(IJ))/eqxmz(ieqsolute(ieqno(IJ),jSM))
         eqym(maicl,1) = eqym(maicl,1) + eqys(ieqcl(IJ),jAP)*eqxy(jZa,ieqcl(IJ))/eqxmz(ieqsolute(ieqcl(IJ),jSM))
         eqym(maicl,2) = eqym(maicl,2) + eqys(ieqcl(IJ),jDP)*eqxy(jZa,ieqcl(IJ))/eqxmz(ieqsolute(ieqcl(IJ),jSM))
   END DO
   DO IJ=1,manio
         eqyp(mciam,1) = eqyp(mciam,1) + eqys(ieqams(IJ),jAP)*eqxy(jZa,ieqams(IJ))/eqxpz(ieqsolute(ieqams(IJ),jSP))
         eqyp(mciam,2) = eqyp(mciam,2) + eqys(ieqams(IJ),jDP)*eqxy(jZa,ieqams(IJ))/eqxpz(ieqsolute(ieqams(IJ),jSP))
         eqyp(mciso,1) = eqyp(mciso,1) + eqys(ieqsos(IJ),jAP)*eqxy(jZa,ieqsos(IJ))/eqxpz(ieqsolute(ieqsos(IJ),jSP))
         eqyp(mciso,2) = eqyp(mciso,2) + eqys(ieqsos(IJ),jDP)*eqxy(jZa,ieqsos(IJ))/eqxpz(ieqsolute(ieqsos(IJ),jSP))
         eqyp(mcipo,1) = eqyp(mcipo,1) + eqys(ieqpos(IJ),jAP)*eqxy(jZa,ieqpos(IJ))/eqxpz(ieqsolute(ieqpos(IJ),jSP))
         eqyp(mcipo,2) = eqyp(mcipo,2) + eqys(ieqpos(IJ),jDP)*eqxy(jZa,ieqpos(IJ))/eqxpz(ieqsolute(ieqpos(IJ),jSP))
         eqyp(mcica,1) = eqyp(mcica,1) + eqys(ieqcas(IJ),jAP)*eqxy(jZa,ieqcas(IJ))/eqxpz(ieqsolute(ieqcas(IJ),jSP))
         eqyp(mcica,2) = eqyp(mcica,2) + eqys(ieqcas(IJ),jDP)*eqxy(jZa,ieqcas(IJ))/eqxpz(ieqsolute(ieqcas(IJ),jSP))
         eqyp(mcimg,1) = eqyp(mcimg,1) + eqys(ieqmgs(IJ),jAP)*eqxy(jZa,ieqmgs(IJ))/eqxpz(ieqsolute(ieqmgs(IJ),jSP))
         eqyp(mcimg,2) = eqyp(mcimg,2) + eqys(ieqmgs(IJ),jDP)*eqxy(jZa,ieqmgs(IJ))/eqxpz(ieqsolute(ieqmgs(IJ),jSP))
   END DO
   ! CAM2_35: Coarse mode nitrate (SR, 12/2021)
   xNO3_2NA=eqys(eqsno,jDP)*eqxy(jZa,eqsno)/eqxpz(ieqsolute(eqsno,jSP))+ &
          & eqys(eqsno,jAP)*eqxy(jZa,eqsno)/eqxpz(ieqsolute(eqsno,jSP))
   xNO3_2CA=eqys(eqcno,jDP)*eqxy(jZa,eqcno)/eqxpz(ieqsolute(eqcno,jSP))+ &
          & eqys(eqcno,jAP)*eqxy(jZa,eqcno)/eqxpz(ieqsolute(eqcno,jSP))
   xNO3_2OTH=eqys(eqpno,jDP)*eqxy(jZa,eqpno)/eqxpz(ieqsolute(eqpno,jSP))+ &
          &  eqys(eqpno,jAP)*eqxy(jZa,eqpno)/eqxpz(ieqsolute(eqpno,jSP))+ &
          &  eqys(eqmno,jDP)*eqxy(jZa,eqmno)/eqxpz(ieqsolute(eqmno,jSP))+ &
          &  eqys(eqmno,jAP)*eqxy(jZa,eqmno)/eqxpz(ieqsolute(eqmno,jSP))
   ! PARTICULATE MATTER
   DO Is=1,jSOL
      IF(.NOT.leqsolute(Is)) CYCLE
         eqsPM=eqsPM+eqys(Is,jDP)*eqxy(jNs,Is)
         eqaPM=eqaPM+eqys(Is,jAP)*eqxy(jNs,Is)
         eqPMt=eqPMt+eqys(Is,jAP)*eqxy(jMs,Is)
         eqPMs=eqPMs+eqys(Is,jDP)*eqxy(jMs,Is)
         eqVOL=eqVOL+eqys(Is,jAP)*eqxy(jMs,Is)/eqxy(jDs,Is)
         eqVOL=eqVOL+eqys(Is,jDP)*eqxy(jMs,Is)/eqxy(jDs,Is)
   END DO
      ! TOTAL PM      [KG/M^3(AIR)]
      eqPMt=eqPMt+eqPMs
      ! TOTAL VOLUME  [M^3/M^3(AIR)]
      ! TOTAL DENSITY [KG/M^3]
      IF(eqVOL > REALZERO) &
      eqRHO=eqPMt/eqVOL
      ! Aerosol Growth Factor [-]
      IF(eqPMT > 0.1e-12_JPRB) THEN
        eqGF =(eqRHO/Dw*eqWH2O/eqPMt+ONE)**(1._JPRB/3._JPRB)
      ENDIF
      ! AEROSOL WATER  [UG/M^3(AIR)]
!     xWH2O  = aqWH2O*1.E9_JPRB !  input total liquid water content (previous time-step)
!     xWH2O  = eqWH2O*1.E9_JPRB !  equilbrium aerosol water content (current  time-step)
!     xWH2O  =  tWH2O*1.E9_JPRB ! output total liquid water content (current  time-step)
      ! AEROSOL WATER  [KG/M^3(AIR)]
!     xWH2O  = aqWH2O           !  input total liquid water content (previous time-step)
!     xWH2O  = eqWH2O           !  equilbrium aerosol water content (current  time-step)
      xWH2O  =  tWH2O           ! output total liquid water content (current  time-step)
      ! TOTAL PM      [UG/M^3(AIR)]
      eqPMt=eqPMt*1.E9_JPRB
      ! DRY PM        [UG/M^3(AIR)]
      eqPMs=eqPMs*1.E9_JPRB
      ! DRY PM        [UMOL/M^3(AIR)]
      eqsPM=eqsPM*1.E6_JPRB
      ! AQUEOUS PM    [UMOL/M^3(AIR)]
      eqaPM=eqaPM*1.E6_JPRB
      ! total AQUEOUS H+ [MOL/Kg(H2O)]
      IF(eqPH > -999._JPRB) THEN
        eqHPLUS = 10._JPRB**(-eqPH)
      ELSE
        eqHPLUS=ZERO
      END IF
      ! DENSITY [G/CM^3]
      eqRHO=eqRHO*1.e-3_JPRB
      eqVOL=eqPMt/eqRHO
      xAD  =float(IDeq)
      xPMt =eqPMt 
      xPMs =eqPMs 
      xsPM =eqsPM 
      xaPM =eqaPM 
      xRHO =eqRHO 
      xVOL =eqVOL 
      xPH  =eqPH  
      xGF  =eqGF  
      xHp  =eqHPLUS
   DO IJ=1,mcati
       xYPa(IJ)=eqyp(IJ,1)
       xYPs(IJ)=eqyp(IJ,2)
   END DO
   DO IJ=1,manio
       xYMa(IJ)=eqym(IJ,1)
       xYMs(IJ)=eqym(IJ,2)
   END DO
   DO IJ=1,jGAS
      xYG(IJ)=eqyg(IJ)
   END DO
!- CAMS2_35 ->
   xAERPH    = MIN(MAX(eqPHaer,0._JPRB),8._JPRB)
   xCLOUDPH  = MIN(MAX(eqPHcloud,0._JPRB),8._JPRB)
   xRAINPH   = MIN(MAX(eqPHrain,0._JPRB),8._JPRB)
   xAERWAT   = aqWH2O+eqWH2O
   IF (LHOOK) CALL DR_HOOK('AER_EQSAM4CLIM',1,ZHOOK_HANDLE)
!- CAMS2_35 <-
!____________________________________________________________________________________________
CONTAINS

FUNCTION ROOTF (A,B)
USE PARKIND1  ,ONLY : JPRB, JPRD
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
IMPLICIT NONE

REAL(KIND=JPRD),INTENT (IN) :: A,B
REAL(KIND=JPRD)             :: ROOTF

IF (LHOOK) CALL DR_HOOK('AER_EQSAM4CLIM:ROOTF',0,ZHOOK_HANDLE)

ROOTF = 0.50_JPRB * (-A + SQRT(A*A - 4.0_JPRB*B))

IF (LHOOK) CALL DR_HOOK('AER_EQSAM4CLIM:ROOTF',1,ZHOOK_HANDLE)
END FUNCTION ROOTF
!____________________________________________________________________________________________
END SUBROUTINE AER_EQSAM4CLIM

