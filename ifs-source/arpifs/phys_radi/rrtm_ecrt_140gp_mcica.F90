SUBROUTINE RRTM_ECRT_140GP_MCICA&
 &(YDDIMV, YDEAERATM,YDERAD,YGFL,KIDIA, KFDIA, KLON, KLEV, KCOLS ,&
 &PRAER, PAER , PAPH , PAP , PAERTAUL, PAERASYL, PAEROMGL,&
 &PTS  , PTH  , PT  ,&
 &PEMIS, PEMIW,&
 &PQ   , PCO2 , PCH4, PN2O  , PNO2, PC11, PC12, PC22, PCL4, POZN, PCLDF  , PTAUCLDI,&
 &PCLDFRAC, PTAUCLD, PCOLDRY, PWBRODL, PWKL, PWX ,&
 &  PTAUAERL, PAVEL  , PTAVEL , PZ  , PTZ , PTBOUND, PSEMISS , KREFLECT)  

!----compiled for Cray with -h nopattern----

!     Reformatted for F90 by JJMorcrette, ECMWF, 980714

!     Read in atmospheric profile from ECMWF radiation code, and prepare it
!     for use in RRTM.  Set other RRTM input parameters.  Values are passed
!     back through existing RRTM arrays and commons.

!- Modifications

!     2000-05-15 Deborah Salmond  Speed-up
!     JJMorcrette 20050110  McICA version
!        NEC           25-Oct-2007 Optimisations
!     PBechtold+NSemane        09-Jul-2012 Gravity
!     201305 ABozzo added PWBRODL and O2 mix. ratio as in the latest version of RRTMG-LW
!     JJMorcrette 20130805 climatol.stratos.aero. with MACC-derived aero.climatol.

USE YOMDIMV  , ONLY : TDIMV
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST   , ONLY : RG
USE PARRRTM  , ONLY : JPBAND, JPXSEC, JPINPX  
USE YOERAD   , ONLY : TERAD
USE YOEAERATM, ONLY : TEAERATM
USE YOM_YGFL , ONLY : TYPE_GFLD
USE YOMDYNCORE,ONLY : RPLRG

!------------------------------Arguments--------------------------------

IMPLICIT NONE

TYPE(TDIMV)       ,INTENT(IN)    :: YDDIMV
TYPE(TEAERATM)    ,INTENT(INOUT) :: YDEAERATM
TYPE(TERAD)       ,INTENT(INOUT) :: YDERAD
TYPE(TYPE_GFLD)   ,INTENT(INOUT) :: YGFL
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON! Number of atmospheres (longitudes) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV! Number of atmospheric layers 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA, KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KCOLS

REAL(KIND=JPRB)   ,INTENT(IN)    :: PRAER(6,6) ! LW absorption coefficients
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAER(KLON,6,KLEV) ! Aerosol optical thickness
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAERTAUL(KLON,KLEV,16), PAERASYL(KLON,KLEV,16), PAEROMGL(KLON,KLEV,16)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPH(KLON,KLEV+1) ! Interface pressures (Pa)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAP(KLON,KLEV) ! Layer pressures (Pa)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTS(KLON) ! Surface temperature (K)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTH(KLON,KLEV+1) ! Interface temperatures (K)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KLON,KLEV) ! Layer temperature (K)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEMIS(KLON) ! Non-window surface emissivity
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEMIW(KLON) ! Window surface emissivity
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

REAL(KIND=JPRB)   ,INTENT(IN)    :: PCLDF(KLON,KCOLS,KLEV)    ! Cloud fraction
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTAUCLDI(KLON,KLEV,KCOLS) ! Cloud optical depth

REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCLDFRAC(KIDIA:KFDIA,KCOLS,YDDIMV%NFLEVG)   ! Cloud fraction
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTAUCLD(KIDIA:KFDIA,YDDIMV%NFLEVG,KCOLS)    ! Spectral optical thickness

REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCOLDRY(KIDIA:KFDIA,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PWBRODL(KIDIA:KFDIA,YDDIMV%NFLEVG) ! broadening gas column density (mol/cm2)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PWKL(KIDIA:KFDIA,JPINPX,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PWX(KIDIA:KFDIA,JPXSEC,YDDIMV%NFLEVG) ! Amount of trace gases
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTAUAERL(KIDIA:KFDIA,YDDIMV%NFLEVG,JPBAND) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PAVEL(KIDIA:KFDIA,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTAVEL(KIDIA:KFDIA,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PZ(KIDIA:KFDIA,0:YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTZ(KIDIA:KFDIA,0:YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTBOUND(KIDIA:KFDIA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSEMISS(KIDIA:KFDIA,JPBAND) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KREFLECT(KIDIA:KFDIA) 

!      real rch4                       ! CH4 mass mixing ratio
!      real rn2o                       ! N2O mass mixing ratio
!      real rcfc11                     ! CFC11 mass mixing ratio
!      real rcfc12                     ! CFC12 mass mixing ratio
!      real rcfc22                     ! CFC22 mass mixing ratio
!      real rccl4                      ! CCl4  mass mixing ratio
!- from AER
!- from PROFILE             
!- from SURFACE             
REAL(KIND=JPRB) :: ZTAUAER(5)
REAL(KIND=JPRB) :: ZAMD                  ! Effective molecular weight of dry air (g/mol)
REAL(KIND=JPRB) :: ZAMW                  ! Molecular weight of water vapor (g/mol)
REAL(KIND=JPRB) :: ZAMCO2                ! Molecular weight of carbon dioxide (g/mol)
REAL(KIND=JPRB) :: ZAMO                  ! Molecular weight of ozone (g/mol)
REAL(KIND=JPRB) :: ZAMCH4                ! Molecular weight of methane (g/mol)
REAL(KIND=JPRB) :: ZAMN2O                ! Molecular weight of nitrous oxide (g/mol)
REAL(KIND=JPRB) :: ZAMC11                ! Molecular weight of CFC11 (g/mol) - CFCL3
REAL(KIND=JPRB) :: ZAMC12                ! Molecular weight of CFC12 (g/mol) - CF2CL2
REAL(KIND=JPRB) :: ZAMC22                ! Molecular weight of CFC22 (g/mol) - CHF2CL
REAL(KIND=JPRB) :: ZAMCL4                ! Molecular weight of CCl4  (g/mol) - CCL4
REAL(KIND=JPRB) :: ZAVGDRO               ! Avogadro's number (molecules/mole)
REAL(KIND=JPRB) :: ZGRAVIT               ! Gravitational acceleration (cm/s**2)

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

INTEGER(KIND=JPIM) :: IATM, JMOL, IXMAX, J1, J2, IAE, IKL, JK, JCOLS, JL, JLW
INTEGER(KIND=JPIM) :: ITMOL, INXMOL

REAL(KIND=JPRB) :: ZAMM

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ***

! *** mji
! Initialize all molecular amounts and aerosol optical depths to zero here, 
! then pass ECRT amounts into RRTM arrays below.

!      DATA ZWKL /MAXPRDW*0.0/
!      DATA ZWX  /MAXPROD*0.0/
!      DATA KREFLECT /0/

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

!      IREFLECT=KREFLECT
!      NXMOL=KXMOL

IF (LHOOK) CALL DR_HOOK('RRTM_ECRT_140GP_MCICA',0,ZHOOK_HANDLE)
ASSOCIATE(NACTAERO=>YGFL%NACTAERO, &
 & NFLEVG=>YDDIMV%NFLEVG, &
 & LAERCSTR=>YDEAERATM%LAERCSTR, LAERRRTM=>YDEAERATM%LAERRRTM, &
 & LAERADCLI=>YDERAD%LAERADCLI, NAERMACC=>YDERAD%NAERMACC, &
 & NSPMAPL=>YDERAD%NSPMAPL, RAER=>PRAER)
ZGRAVIT=(RG/RPLRG)*1.E2_JPRB

DO JL = KIDIA, KFDIA
  KREFLECT(JL)=0
  INXMOL=2
ENDDO

!DO J1=1,35
! IXINDX(J1)=0
DO J2=1,KLEV
  DO J1=1,35
    DO JL = KIDIA, KFDIA
      PWKL(JL,J1,J2)=0.0_JPRB 
    ENDDO
  ENDDO
ENDDO
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
DO JL = KIDIA, KFDIA
  PSEMISS(JL,1)  = PEMIS(JL)
  PSEMISS(JL,2)  = PEMIS(JL)
  PSEMISS(JL,3)  = PEMIS(JL)
  PSEMISS(JL,4)  = PEMIS(JL)
  PSEMISS(JL,5)  = PEMIS(JL)
  PSEMISS(JL,6)  = PEMIW(JL)
  PSEMISS(JL,7)  = PEMIW(JL)
  PSEMISS(JL,8)  = PEMIW(JL)
  PSEMISS(JL,9)  = PEMIS(JL)
  PSEMISS(JL,10) = PEMIS(JL)
  PSEMISS(JL,11) = PEMIS(JL)
  PSEMISS(JL,12) = PEMIS(JL)
  PSEMISS(JL,13) = PEMIS(JL)
  PSEMISS(JL,14) = PEMIS(JL)
  PSEMISS(JL,15) = PEMIS(JL)
  PSEMISS(JL,16) = PEMIS(JL)

!     Set surface temperature.  

  PTBOUND(JL) = PTS(JL)

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
!     variables count from the top down and must be reversed 
  ITMOL = 7
  PZ(JL,0) = PAPH(JL,KLEV+1)/100._JPRB
  PTZ(JL,0) = PTH(JL,KLEV+1)
ENDDO

  DO JK = 1, KLEV
DO JL = KIDIA, KFDIA
    PAVEL(JL,JK) = PAP(JL,KLEV-JK+1)/100._JPRB
    PTAVEL(JL,JK) = PT(JL,KLEV-JK+1)
    PZ(JL,JK) = PAPH(JL,KLEV-JK+1)/100._JPRB
    PTZ(JL,JK) = PTH(JL,KLEV-JK+1)
    PWKL(JL,1,JK) = PQ(JL,KLEV-JK+1)*ZAMD/ZAMW
    PWKL(JL,2,JK) = PCO2(JL,KLEV-JK+1)*ZAMD/ZAMCO2
    PWKL(JL,3,JK) = POZN(JL,KLEV-JK+1)*ZAMD/ZAMO
    PWKL(JL,4,JK) = PN2O(JL,KLEV-JK+1)*ZAMD/ZAMN2O
    PWKL(JL,6,JK) = PCH4(JL,KLEV-JK+1)*ZAMD/ZAMCH4
    PWKL(JL,7,JK) = 0.209488_JPRB
    ZAMM = (1.0_JPRB-PWKL(JL,1,JK))*ZAMD + PWKL(JL,1,JK)*ZAMW
    PCOLDRY(JL,JK) = (PZ(JL,JK-1)-PZ(JL,JK))*1.E3_JPRB*ZAVGDRO/(ZGRAVIT*ZAMM*(1.0_JPRB+PWKL(JL,1,JK)))
ENDDO
  ENDDO

!- If prognostic aerosols with proper RRTM optical properties, 
!  or alternatively if MACC-aerosol climatology is to be used, 
!  fill the RRTM aerosol arrays without. Stratospheric/tropospheric
!  background types are added in aer_rrtm when LAERCSRT=F. If 
!  LAERCSRT=T, only the climatological stratospheric background
!  is consedered and this is added to the sullphates type in
!  radlswr.
!  AB added multiplication by (1-omega) to PAERTAUL 
!  to compute the absorption optical thickness
!  (otherwise the extinction is used)

IF (LAERRRTM .OR. LAERADCLI) THEN
    DO JLW=1,16
      DO JK=1,KLEV
        IKL=KLEV-JK+1
          DO JL=KIDIA,KFDIA
            PTAUAERL(JL,JK,JLW)=PAERTAUL(JL,IKL,JLW)*(1.0_JPRB-PAEROMGL(JL,IKL,JLW))
          ENDDO
      ENDDO
    ENDDO

ELSEIF (.NOT.LAERRRTM .AND. .NOT.LAERADCLI) THEN

!- Fill RRTM aerosol arrays with operational ECMWF (Tegen) aerosols,
!  do the mixing and distribute over the 16 spectral intervals

  DO JK=1,KLEV
    IKL=KLEV-JK+1
    DO JL = KIDIA, KFDIA
    IAE=1
    ZTAUAER(IAE) =&
   & RAER(IAE,1)*PAER(JL,1,IKL)+RAER(IAE,2)*PAER(JL,2,IKL)&
   & +RAER(IAE,3)*PAER(JL,3,IKL)+RAER(IAE,4)*PAER(JL,4,IKL)&
   & +RAER(IAE,5)*PAER(JL,5,IKL)+RAER(IAE,6)*PAER(JL,6,IKL)  
    PTAUAERL(JL,JK, 1)=ZTAUAER(1)
    PTAUAERL(JL,JK, 2)=ZTAUAER(1)
    IAE=2
    ZTAUAER(IAE) =&
   & RAER(IAE,1)*PAER(JL,1,IKL)+RAER(IAE,2)*PAER(JL,2,IKL)&
   & +RAER(IAE,3)*PAER(JL,3,IKL)+RAER(IAE,4)*PAER(JL,4,IKL)&
   & +RAER(IAE,5)*PAER(JL,5,IKL)+RAER(IAE,6)*PAER(JL,6,IKL)  
    PTAUAERL(JL,JK, 3)=ZTAUAER(2)
    PTAUAERL(JL,JK, 4)=ZTAUAER(2)
    PTAUAERL(JL,JK, 5)=ZTAUAER(2)
    IAE=3
    ZTAUAER(IAE) =&
   & RAER(IAE,1)*PAER(JL,1,IKL)+RAER(IAE,2)*PAER(JL,2,IKL)&
   & +RAER(IAE,3)*PAER(JL,3,IKL)+RAER(IAE,4)*PAER(JL,4,IKL)&
   & +RAER(IAE,5)*PAER(JL,5,IKL)+RAER(IAE,6)*PAER(JL,6,IKL)  
    PTAUAERL(JL,JK, 6)=ZTAUAER(3)
    PTAUAERL(JL,JK, 8)=ZTAUAER(3)
    PTAUAERL(JL,JK, 9)=ZTAUAER(3)
    IAE=4
    ZTAUAER(IAE) =&
   & RAER(IAE,1)*PAER(JL,1,IKL)+RAER(IAE,2)*PAER(JL,2,IKL)&
   & +RAER(IAE,3)*PAER(JL,3,IKL)+RAER(IAE,4)*PAER(JL,4,IKL)&
   & +RAER(IAE,5)*PAER(JL,5,IKL)+RAER(IAE,6)*PAER(JL,6,IKL)  
    PTAUAERL(JL,JK, 7)=ZTAUAER(4)
    IAE=5
    ZTAUAER(IAE) =&
   & RAER(IAE,1)*PAER(JL,1,IKL)+RAER(IAE,2)*PAER(JL,2,IKL)&
   & +RAER(IAE,3)*PAER(JL,3,IKL)+RAER(IAE,4)*PAER(JL,4,IKL)&
   & +RAER(IAE,5)*PAER(JL,5,IKL)+RAER(IAE,6)*PAER(JL,6,IKL)  
    PTAUAERL(JL,JK,10)=ZTAUAER(5)
    PTAUAERL(JL,JK,11)=ZTAUAER(5)
    PTAUAERL(JL,JK,12)=ZTAUAER(5)
    PTAUAERL(JL,JK,13)=ZTAUAER(5)
    PTAUAERL(JL,JK,14)=ZTAUAER(5)
    PTAUAERL(JL,JK,15)=ZTAUAER(5)
    PTAUAERL(JL,JK,16)=ZTAUAER(5)
    ENDDO
  ENDDO
ENDIF

  DO J2=1,KLEV
    DO J1=1,JPXSEC
DO JL = KIDIA, KFDIA
      PWX(JL,J1,J2)=0.0_JPRB
ENDDO
    ENDDO
  ENDDO

  DO JK = 1, KLEV
DO JL = KIDIA, KFDIA
!- Set cross section molecule amounts from ECRT; convert to vmr
    PWX(JL,1,JK) = PCL4(JL,KLEV-JK+1) * ZAMD/ZAMCL4
    PWX(JL,2,JK) = PC11(JL,KLEV-JK+1) * ZAMD/ZAMC11
    PWX(JL,3,JK) = PC12(JL,KLEV-JK+1) * ZAMD/ZAMC12
    PWX(JL,4,JK) = PC22(JL,KLEV-JK+1) * ZAMD/ZAMC22
    PWX(JL,1,JK) = PCOLDRY(JL,JK) * PWX(JL,1,JK) * 1.E-20_JPRB
    PWX(JL,2,JK) = PCOLDRY(JL,JK) * PWX(JL,2,JK) * 1.E-20_JPRB
    PWX(JL,3,JK) = PCOLDRY(JL,JK) * PWX(JL,3,JK) * 1.E-20_JPRB
    PWX(JL,4,JK) = PCOLDRY(JL,JK) * PWX(JL,4,JK) * 1.E-20_JPRB

!- Here, all molecules in WKL and WX are in volume mixing ratio; convert to
!  molec/cm2 based on COLDRY for use in RRTM

!CDIR UNROLL=6
ZSUMMOL = 0.0_JPRB
!AB broadening gases
    DO JMOL = 2, ITMOL
      ZSUMMOL = ZSUMMOL + PWKL(JL,JMOL,JK)
    ENDDO
    PWBRODL(JL,JK) = PCOLDRY(JL,JK) * (1._JPRB - ZSUMMOL)
    DO JMOL = 1, ITMOL
      PWKL(JL,JMOL,JK) = PCOLDRY(JL,JK) * PWKL(JL,JMOL,JK)
    ENDDO    
ENDDO
  ENDDO

!- McICA: No overlap; simple copy of optical thickness; layer cloud cover is 0. or 1.

  DO JK=1,KLEV
    DO JCOLS=1,KCOLS
DO JL = KIDIA, KFDIA
      PCLDFRAC(JL,JCOLS,JK)=PCLDF(JL,JCOLS,JK)
      PTAUCLD(JL,JK,JCOLS) =PTAUCLDI(JL,JK,JCOLS)
ENDDO
    ENDDO
  ENDDO

!     ------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('RRTM_ECRT_140GP_MCICA',1,ZHOOK_HANDLE)
END SUBROUTINE RRTM_ECRT_140GP_MCICA
