!****************** SUBROUTINE RRTM_ECRT_140GP **************************

SUBROUTINE RRTM_ECRT_140GP&
 & (YDDIMV,YDERAD,KLON,KIDIA,KFDIA, KLEV, KCLD,&
 & PRAER, PAER , PAPH , PAP,&
 & PTS  , PTH  , PT,&
 & P_ZEMIS, P_ZEMIW,&
 & PQ   , PCO2, PCH4, PN2O, PNO2, PC11, PC12, PC22, PCL4, POZN, PCLDF, PTAUCLD, PTCLEAR,&
 & P_CLDFRAC,P_TAUCLD,P_COLDRY, P_WBRODL,P_WKL,P_WX,&
 & P_TAUAERL,PAVEL,P_TAVEL,PZ,P_TZ,P_TBOUND,K_NLAYERS,P_SEMISS)  

!     Reformatted for F90 by JJMorcrette, ECMWF, 980714

!     Read in atmospheric profile from ECMWF radiation code, and prepare it
!     for use in RRTM.  Set other RRTM input parameters.  Values are passed
!     back through existing RRTM arrays and commons.

!- Modifications

!     2000-05-15 Deborah Salmond  Speed-up
!     JJMorcrette 20071015  3D fields of CO2, CH4, N2O and NO2
!        NEC           25-Oct-2007 Optimisations
!     2008-08-25 Yann Seity : add  missing cases NOVLP > 4
!     PBechtold+NSemane        09-Jul-2012 Gravity
!     201305 ABozzo PWBRODL,O2
!      R. El Khatib 12-Aug-2016 optimization


USE YOMDIMV  , ONLY : TDIMV
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE PARRRTM  , ONLY : JPBAND, JPXSEC, JPINPX  
USE YOERAD   , ONLY : TERAD
USE YOMCST   , ONLY : RG
USE YOMDYNCORE,ONLY : RPLRG

IMPLICIT NONE

!------------------------------Arguments--------------------------------

TYPE(TDIMV)       ,INTENT(IN)    :: YDDIMV
TYPE(TERAD)       ,INTENT(INOUT) :: YDERAD
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON! Number of atmospheres (longitudes) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV! Number of atmospheric layers 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KCLD(KIDIA:KFDIA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRAER(6,6) ! LW absorption coefficients
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAER(KLON,6,KLEV) ! Aerosol optical thickness
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPH(KLON,KLEV+1) ! Interface pressures (Pa)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAP(KLON,KLEV) ! Layer pressures (Pa)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTS(KLON) ! Surface temperature (K)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTH(KLON,KLEV+1) ! Interface temperatures (K)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KLON,KLEV) ! Layer temperature (K)
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_ZEMIS(KLON) ! Non-window surface emissivity
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_ZEMIW(KLON) ! Window surface emissivity
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQ(KLON,KLEV) ! H2O specific humidity (mmr)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCO2(KLON,KLEV) ! CO2 mass mixing ratio
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCH4(KLON,KLEV) ! CH4 mass mixing ratio
REAL(KIND=JPRB)   ,INTENT(IN)    :: PN2O(KLON,KLEV) ! N2O mass mixing ratio
REAL(KIND=JPRB)   ,INTENT(IN)    :: PNO2(KLON,KLEV) ! NO2 mass mixing ratio
REAL(KIND=JPRB)   ,INTENT(IN)    :: PC11(KLON,KLEV) ! CFC11 mass mixing ratio
REAL(KIND=JPRB)   ,INTENT(IN)    :: PC12(KLON,KLEV) ! CFC12 mass mixing ratio
REAL(KIND=JPRB)   ,INTENT(IN)    :: PC22(KLON,KLEV) ! CFC22 mass mixing ratio
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCL4(KLON,KLEV) ! CCL4  mass mixing ratio
REAL(KIND=JPRB)   ,INTENT(IN)    :: POZN(KLON,KLEV) ! O3 mass mixing ratio
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCLDF(KLON,KLEV) ! Cloud fraction
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTAUCLD(KLON,KLEV,JPBAND) ! Cloud optical depth
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTCLEAR (KIDIA:KFDIA)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_CLDFRAC(KIDIA:KFDIA,YDDIMV%NFLEVG) ! Cloud fraction
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_TAUCLD(KIDIA:KFDIA,YDDIMV%NFLEVG,JPBAND) ! Spectral optical thickness
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_COLDRY(KIDIA:KFDIA,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_WBRODL(KIDIA:KFDIA,YDDIMV%NFLEVG) ! broadening gas column density (mol/cm2)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_WKL(KIDIA:KFDIA,JPINPX,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_WX(KIDIA:KFDIA,JPXSEC,YDDIMV%NFLEVG) ! Amount of trace gases
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_TAUAERL(KIDIA:KFDIA,YDDIMV%NFLEVG,JPBAND) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PAVEL(KIDIA:KFDIA,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_TAVEL(KIDIA:KFDIA,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PZ(KIDIA:KFDIA,0:YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_TZ(KIDIA:KFDIA,0:YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_TBOUND(KIDIA:KFDIA) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: K_NLAYERS(KIDIA:KFDIA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_SEMISS(KIDIA:KFDIA,JPBAND) 
!      real rch4                       ! CH4 mass mixing ratio
!      real rn2o                       ! N2O mass mixing ratio
!      real rcfc11                     ! CFC11 mass mixing ratio
!      real rcfc12                     ! CFC12 mass mixing ratio
!      real rcfc22                     ! CFC22 mass mixing ratio
!      real rccl4                      ! CCL4  mass mixing ratio
!- from AER
!- from PROFILE             
!- from SURFACE             
REAL(KIND=JPRB) :: ZAMD                  ! Effective molecular weight of dry air (g/mol)
REAL(KIND=JPRB) :: ZAMW                  ! Molecular weight of water vapor (g/mol)
REAL(KIND=JPRB) :: ZAMCO2                ! Molecular weight of carbon dioxide (g/mol)
REAL(KIND=JPRB) :: ZAMO                  ! Molecular weight of ozone (g/mol)
REAL(KIND=JPRB) :: ZAMCH4                ! Molecular weight of methane (g/mol)
REAL(KIND=JPRB) :: ZAMN2O                ! Molecular weight of nitrous oxide (g/mol)
REAL(KIND=JPRB) :: ZAMC11                ! Molecular weight of CFC11 (g/mol) - CFCL3
REAL(KIND=JPRB) :: ZAMC12                ! Molecular weight of CFC12 (g/mol) - CF2CL2
REAL(KIND=JPRB) :: ZAMC22                ! Molecular weight of CFC22 (g/mol) - CF2CL2
REAL(KIND=JPRB) :: ZAMCL4                ! Molecular weight of CCL4  (g/mol) - CCl4
REAL(KIND=JPRB) :: ZAVGDRO               ! Avogadro's number (molecules/mole)
REAL(KIND=JPRB) :: ZGRAVIT               ! Gravitational acceleration (cm/sec2)

REAL(KIND=JPRB) :: ZSUMMOL

! Atomic weights for conversion from mass to volume mixing ratios; these
!  are the same values used in ECRT to assure accurate conversion to vmr
DATA ZAMD   /  28.970_JPRB    /
DATA ZAMW   /  18.0154_JPRB   /
DATA ZAMCO2 /  44.011_JPRB    /
DATA ZAMO   /  47.9982_JPRB   /
DATA ZAMCH4 /  16.043_JPRB    /
DATA ZAMN2O /  44.013_JPRB    /
DATA ZAMC11 / 137.3686_JPRB   /
DATA ZAMC12 / 120.9140_JPRB   /
DATA ZAMC22 /  86.4690_JPRB   /
DATA ZAMCL4 / 153.8230_JPRB   /
DATA ZAVGDRO/ 6.02214E23_JPRB /

INTEGER(KIND=JPIM) :: IATM, IMOL, IXMAX, J1, J2, IAE, JB, JK, JLEV
INTEGER(KIND=JPIM), PARAMETER :: I_NMOL =7
INTEGER(KIND=JPIM) :: JLON

REAL(KIND=JPRB) :: ZAMM, ZCLDLY(KIDIA:KFDIA,KLEV), ZCLOUD(KIDIA:KFDIA), ZEPSEC
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ***

! *** mji
! Initialize all molecular amounts and aerosol optical depths to zero here, 
! then pass ECRT amounts into RRTM arrays below.

!      DATA ZWKL /MAXPRDW*0.0/
!      DATA ZWX  /MAXPROD*0.0/

! Activate cross section molecules:
!     NXMOL     - number of cross-sections input by user
!     IXINDX(I) - index of cross-section molecule corresponding to Ith
!                 cross-section specified by user
!                 = 0 -- not allowed in RRTM
!                 = 1 -- CCL4
!                 = 2 -- CFC11
!                 = 3 -- CFC12
!                 = 4 -- CFC22
!      DATA KXMOL  /2/
!      DATA KXINDX /0,2,3,0,31*0/

!      NXMOL=KXMOL

IF (LHOOK) CALL DR_HOOK('RRTM_ECRT_140GP',0,ZHOOK_HANDLE)
ASSOCIATE(NFLEVG=>YDDIMV%NFLEVG, &
 & NOVLP=>YDERAD%NOVLP, RAER=>PRAER)
ZGRAVIT=(RG/RPLRG)*1.E2_JPRB

! IXINDX(J1)=0
!IXINDX(2)=2
!IXINDX(3)=3

!     Set parameters needed for RRTM execution:
IATM    = 0
!      IXSECT  = 1
!      NUMANGS = 0
!      IOUT    = -1
IXMAX   = 4

!     Bands 6,7,8 are considered the 'window' and allowed to have a
!     different surface emissivity (as in ECMWF).  Eli wrote this part....

DO J2=1,5
  P_SEMISS(KIDIA:KFDIA,J2)  = P_ZEMIS(KIDIA:KFDIA)
ENDDO
DO J2=6,8
  P_SEMISS(KIDIA:KFDIA,J2)  = P_ZEMIW(KIDIA:KFDIA)
ENDDO
DO J2=9,16
  P_SEMISS(KIDIA:KFDIA,J2)  = P_ZEMIS(KIDIA:KFDIA)
ENDDO

!     Set surface temperature.  

P_TBOUND(KIDIA:KFDIA) = PTS(KIDIA:KFDIA)
K_NLAYERS(KIDIA:KFDIA) = KLEV

!     Install ECRT arrays into RRTM arrays for pressure, temperature,
!     and molecular amounts.  Pressures are converted from Pascals
!     (ECRT) to mb (RRTM).  H2O, CO2, O3 and trace gas amounts are 
!     converted from mass mixing ratio to volume mixing ratio.  CO2
!     converted with same dry air and CO2 molecular weights used in 
!     ECRT to assure correct conversion back to the proper CO2 vmr.
!     The dry air column COLDRY (in molec/cm2) is calculated from 
!     the level pressures PZ (in mb) based on the hydrostatic equation
!     and includes a correction to account for H2O in the layer.  The
!     molecular weight of moist air (amm) is calculated for each layer.
!     Note: RRTM levels count from bottom to top, while the ECRT input
!     variables count from the top down and must be reversed here.

DO JLEV=0,KLEV
  P_TZ(KIDIA:KFDIA,JLEV) = PTH(KIDIA:KFDIA,KLEV-JLEV+1)
ENDDO
PZ(KIDIA:KFDIA,0) = PAPH(KIDIA:KFDIA,KLEV+1)/100._JPRB
DO JLEV=1,KLEV
  JK=KLEV-JLEV+1
  DO JLON=KIDIA,KFDIA
    PZ(JLON,JLEV) = PAPH(JLON,JK)/100._JPRB
    PAVEL(JLON,JLEV) = PAP(JLON,JK)/100._JPRB
    P_WKL(JLON,1,JLEV) = PQ(JLON,JK)  *ZAMD/ZAMW
    ZAMM = (1-P_WKL(JLON,1,JLEV))*ZAMD + P_WKL(JLON,1,JLEV)*ZAMW
    P_COLDRY(JLON,JLEV) = (PZ(JLON,JLEV-1)-PZ(JLON,JLEV))*1.E3_JPRB*ZAVGDRO/(ZGRAVIT*ZAMM*(1+P_WKL(JLON,1,JLEV)))
!- Set cross section molecule amounts from ECRT; convert to vmr
    P_WX(JLON,1,JLEV) = PCL4(JLON,JK) * ZAMD/ZAMCL4
    P_WX(JLON,1,JLEV) = P_COLDRY(JLON,JLEV) * P_WX(JLON,1,JLEV) * 1.E-20_JPRB
    P_WX(JLON,2,JLEV) = PC11(JLON,JK) * ZAMD/ZAMC11
    P_WX(JLON,2,JLEV) = P_COLDRY(JLON,JLEV) * P_WX(JLON,2,JLEV) * 1.E-20_JPRB
    P_WX(JLON,3,JLEV) = PC12(JLON,JK) * ZAMD/ZAMC12
    P_WX(JLON,3,JLEV) = P_COLDRY(JLON,JLEV) * P_WX(JLON,3,JLEV) * 1.E-20_JPRB
    P_WX(JLON,4,JLEV) = PC22(JLON,JK) * ZAMD/ZAMC22
    P_WX(JLON,4,JLEV) = P_COLDRY(JLON,JLEV) * P_WX(JLON,4,JLEV) * 1.E-20_JPRB
    P_WKL(JLON,2,JLEV) = PCO2(JLON,JK)*ZAMD/ZAMCO2
    P_WKL(JLON,3,JLEV) = POZN(JLON,JK)*ZAMD/ZAMO
    P_WKL(JLON,4,JLEV) = PN2O(JLON,JK)*ZAMD/ZAMN2O
    P_WKL(JLON,5,JLEV) = 0.0_JPRB
    P_WKL(JLON,6,JLEV) = PCH4(JLON,JK)*ZAMD/ZAMCH4

!- Fill RRTM aerosol arrays with operational ECMWF aerosols,
!  do the mixing and distribute over the 16 spectral intervals

    IAE=1
    P_TAUAERL(JLON,JLEV, 1)=&
   & RAER(IAE,1)*PAER(JLON,1,JK)+RAER(IAE,2)*PAER(JLON,2,JK)&
   & +RAER(IAE,3)*PAER(JLON,3,JK)+RAER(IAE,4)*PAER(JLON,4,JK)&
   & +RAER(IAE,5)*PAER(JLON,5,JK)+RAER(IAE,6)*PAER(JLON,6,JK)  
    IAE=2
    P_TAUAERL(JLON,JLEV, 3)=&
   & RAER(IAE,1)*PAER(JLON,1,JK)+RAER(IAE,2)*PAER(JLON,2,JK)&
   & +RAER(IAE,3)*PAER(JLON,3,JK)+RAER(IAE,4)*PAER(JLON,4,JK)&
   & +RAER(IAE,5)*PAER(JLON,5,JK)+RAER(IAE,6)*PAER(JLON,6,JK)  
    IAE=3
    P_TAUAERL(JLON,JLEV, 6)=&
   & RAER(IAE,1)*PAER(JLON,1,JK)+RAER(IAE,2)*PAER(JLON,2,JK)&
   & +RAER(IAE,3)*PAER(JLON,3,JK)+RAER(IAE,4)*PAER(JLON,4,JK)&
   & +RAER(IAE,5)*PAER(JLON,5,JK)+RAER(IAE,6)*PAER(JLON,6,JK)  
    IAE=4
    P_TAUAERL(JLON,JLEV, 7)=&
   & RAER(IAE,1)*PAER(JLON,1,JK)+RAER(IAE,2)*PAER(JLON,2,JK)&
   & +RAER(IAE,3)*PAER(JLON,3,JK)+RAER(IAE,4)*PAER(JLON,4,JK)&
   & +RAER(IAE,5)*PAER(JLON,5,JK)+RAER(IAE,6)*PAER(JLON,6,JK)  
    IAE=5
    P_TAUAERL(JLON,JLEV,10)=&
   & RAER(IAE,1)*PAER(JLON,1,JK)+RAER(IAE,2)*PAER(JLON,2,JK)&
   & +RAER(IAE,3)*PAER(JLON,3,JK)+RAER(IAE,4)*PAER(JLON,4,JK)&
   & +RAER(IAE,5)*PAER(JLON,5,JK)+RAER(IAE,6)*PAER(JLON,6,JK)  
  ENDDO
  P_TAVEL(KIDIA:KFDIA,JLEV) = PT(KIDIA:KFDIA,JK)
  P_TAUAERL(KIDIA:KFDIA,JLEV, 2)=P_TAUAERL(KIDIA:KFDIA,JLEV, 1)
  P_TAUAERL(KIDIA:KFDIA,JLEV, 4)=P_TAUAERL(KIDIA:KFDIA,JLEV, 3)
  P_TAUAERL(KIDIA:KFDIA,JLEV, 5)=P_TAUAERL(KIDIA:KFDIA,JLEV, 3)
  P_TAUAERL(KIDIA:KFDIA,JLEV, 8)=P_TAUAERL(KIDIA:KFDIA,JLEV, 6)
  P_TAUAERL(KIDIA:KFDIA,JLEV, 9)=P_TAUAERL(KIDIA:KFDIA,JLEV, 6)
  P_TAUAERL(KIDIA:KFDIA,JLEV,11)=P_TAUAERL(KIDIA:KFDIA,JLEV,10)
  P_TAUAERL(KIDIA:KFDIA,JLEV,12)=P_TAUAERL(KIDIA:KFDIA,JLEV,10)
  P_TAUAERL(KIDIA:KFDIA,JLEV,13)=P_TAUAERL(KIDIA:KFDIA,JLEV,10)
  P_TAUAERL(KIDIA:KFDIA,JLEV,14)=P_TAUAERL(KIDIA:KFDIA,JLEV,10)
  P_TAUAERL(KIDIA:KFDIA,JLEV,15)=P_TAUAERL(KIDIA:KFDIA,JLEV,10)
  P_TAUAERL(KIDIA:KFDIA,JLEV,16)=P_TAUAERL(KIDIA:KFDIA,JLEV,10)
  P_WX(KIDIA:KFDIA,5:,JLEV)=0.0_JPRB
  P_WKL(KIDIA:KFDIA,7,JLEV) = 0.209488_JPRB
  DO J1=8,35
    P_WKL(KIDIA:KFDIA,J1,JLEV)=0.0_JPRB 
  ENDDO

!- Here, all molecules in WKL and WX are in volume mixing ratio; convert to
!  molec/cm2 based on COLDRY for use in RRTM

  DO JLON=KIDIA,KFDIA
!CDIR UNROLL=I_NMOL
    ZSUMMOL = 0.0_JPRB
!AB broadening gases
    DO IMOL = 2, I_NMOL
      ZSUMMOL = ZSUMMOL + P_WKL(JLON,IMOL,JLEV)
    ENDDO
    P_WBRODL(JLON,JLEV) = P_COLDRY(JLON,JLEV) * (1._JPRB - ZSUMMOL)
    DO IMOL = 1, I_NMOL
      P_WKL(JLON,IMOL,JLEV) = P_COLDRY(JLON,JLEV) * P_WKL(JLON,IMOL,JLEV)
    ENDDO  
  
! DO IX = 1,JPXSEC
! IF (IXINDX(IX)  /=  0) THEN
!     WX(IXINDX(IX),L) = COLDRY(L) * WX(IX,L) * 1.E-20_JPRB
! ENDIF
! ENDDO  

  ENDDO
ENDDO

!- Approximate treatment for various cloud overlaps
ZEPSEC=1.E-03_JPRB

DO JK=1,KLEV
  DO JLON=KIDIA,KFDIA
    IF (PCLDF(JLON,JK) > ZEPSEC) THEN
      ZCLDLY(JLON,JK)=PCLDF(JLON,JK)
    ELSE
      ZCLDLY(JLON,JK)=0.0_JPRB
    ENDIF
  ENDDO
ENDDO

IF ((NOVLP == 1).OR.(NOVLP ==6).OR.(NOVLP ==8)) THEN

  PTCLEAR(KIDIA:KFDIA)=1.0_JPRB
  ZCLOUD(KIDIA:KFDIA)=0.0_JPRB
  DO JK=1,KLEV
    DO JLON=KIDIA,KFDIA
      PTCLEAR(JLON)=PTCLEAR(JLON)&
       & *(1.0_JPRB-MAX( ZCLDLY(JLON,JK) , ZCLOUD(JLON) ))&
       & /(1.0_JPRB-MIN( ZCLOUD(JLON) , 1.0_JPRB-ZEPSEC ))  
      PTCLEAR(JLON)=1.0_JPRB-(1.0_JPRB - PTCLEAR(JLON))
      ZCLOUD(JLON) = ZCLDLY(JLON,JK)
    ENDDO
  ENDDO

ELSEIF ((NOVLP == 2).OR.(NOVLP ==7)) THEN

  ZCLOUD(KIDIA:KFDIA)=0.0_JPRB
  DO JK=1,KLEV
    DO JLON=KIDIA,KFDIA
      ZCLOUD(JLON) = MAX( ZCLDLY(JLON,JK) , ZCLOUD(JLON) )
    ENDDO
  ENDDO
  DO JLON=KIDIA,KFDIA
    PTCLEAR(JLON)=1.0_JPRB-ZCLOUD(JLON)
  ENDDO

ELSEIF ((NOVLP == 3).OR.(NOVLP ==5)) THEN

  PTCLEAR(KIDIA:KFDIA)=1.0_JPRB
  DO JK=1,KLEV
    DO JLON=KIDIA,KFDIA
      PTCLEAR(JLON)=1.0_JPRB-(1.0_JPRB - (PTCLEAR(JLON) * (1.0_JPRB-ZCLDLY(JLON,JK))))
    ENDDO
  ENDDO

ELSEIF (NOVLP == 4) THEN

  PTCLEAR(KIDIA:KFDIA)=1.0_JPRB

ENDIF

! Transfer cloud fraction and cloud optical depth to RRTM arrays; 
! invert array index for PCLDF to go from bottom to top for RRTM
!- clear-sky column
!- cloudy column
!   The diffusivity factor (Savijarvi, 1997) on the cloud optical 
!   thickness TAUCLD has already been applied in RADLSW


DO JLON=KIDIA,KFDIA
  IF (PTCLEAR(JLON)  >  1.0_JPRB-ZEPSEC) THEN
    KCLD(JLON)=0
  ELSE
    KCLD(JLON)=1
  ENDIF
ENDDO

DO JLEV = 1, KLEV
  DO JLON=KIDIA,KFDIA
    P_CLDFRAC(JLON,JLEV) = PCLDF(JLON,JLEV)*KCLD(JLON)
  ENDDO
ENDDO

DO JB=1,JPBAND
  DO JLEV = 1, KLEV
    DO JLON=KIDIA,KFDIA
      P_TAUCLD(JLON,JLEV,JB) = PTAUCLD(JLON,JLEV,JB)*KCLD(JLON)
    ENDDO
  ENDDO
ENDDO

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('RRTM_ECRT_140GP',1,ZHOOK_HANDLE)
END SUBROUTINE RRTM_ECRT_140GP
