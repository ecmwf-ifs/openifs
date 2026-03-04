SUBROUTINE RRTM_PREPARE_GASES &
 &( KIDIA, KFDIA, KLON, KLEV, &
 &  PAPH , PAP , &
 &  PTH  , PT  , &
 &  PQ   , PCO2 , PCH4, PN2O  , PNO2, PC11, PC12, PC22, PCL4, POZN, &
 &  PCOLDRY, PWBRODL, PWKL, PWX , &
 &  PAVEL  , PTAVEL , PZ  , PTZ , KREFLECT)  

!----compiled for Cray with -h nopattern----

!     Prepare the units of the gas concentrations for the longwave
!     RRTM gas absorption model.  This file is adapted from
!     rrtm_ecrt_140gp_mcica.F90, written mainly by Jean-Jacques
!     Morcrette.

!- Original
!     2015-07-15  Robin Hogan

!- Modifications

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST   , ONLY : RG
USE PARRRTM  , ONLY : JPXSEC, JPINPX  
USE YOMDYNCORE,ONLY : RPLRG

!------------------------------Arguments--------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KLON! Number of atmospheres (longitudes) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV! Number of atmospheric layers 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA, KFDIA 

REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPH(KLON,KLEV+1) ! Interface pressures (Pa)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAP(KLON,KLEV) ! Layer pressures (Pa)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTH(KLON,KLEV+1) ! Interface temperatures (K)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KLON,KLEV) ! Layer temperature (K)
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

REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCOLDRY(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PWBRODL(KIDIA:KFDIA,KLEV) ! broadening gas column density (mol/cm2)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PWKL(KIDIA:KFDIA,JPINPX,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PWX(KIDIA:KFDIA,JPXSEC,KLEV) ! Amount of trace gases
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PAVEL(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTAVEL(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PZ(KIDIA:KFDIA,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTZ(KIDIA:KFDIA,0:KLEV) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KREFLECT(KIDIA:KFDIA) 

!      real rch4                       ! CH4 mass mixing ratio
!      real rn2o                       ! N2O mass mixing ratio
!      real rcfc11                     ! CFC11 mass mixing ratio
!      real rcfc12                     ! CFC12 mass mixing ratio
!      real rcfc22                     ! CFC22 mass mixing ratio
!      real rccl4                      ! CCl4  mass mixing ratio
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
REAL(KIND=JPRB) :: ZAMC22                ! Molecular weight of CFC22 (g/mol) - CHF2CL
REAL(KIND=JPRB) :: ZAMCL4                ! Molecular weight of CCl4  (g/mol) - CCL4
REAL(KIND=JPRB) :: ZAVGDRO               ! Avogadro's number (molecules/mole)
REAL(KIND=JPRB) :: ZGRAVIT               ! Gravitational acceleration (cm/s**2)

REAL(KIND=JPRB) :: ZSUMMOL

! Atomic weights for conversion from mass to volume mixing ratios; these
!  are the same values used in ECRT to assure accurate conversion to vmr
data ZAMD   /  28.970_JPRB    /
data ZAMW   /  18.0154_JPRB   /
data ZAMCO2 /  44.011_JPRB    /
data ZAMO   /  47.9982_JPRB   /
data ZAMCH4 /  16.043_JPRB    /
data ZAMN2O /  44.013_JPRB    /
data ZAMC11 / 137.3686_JPRB   /
data ZAMC12 / 120.9140_JPRB   /
data ZAMC22 /  86.4690_JPRB   /
data ZAMCL4 / 153.8230_JPRB   /
data ZAVGDRO/ 6.02214E23_JPRB /

INTEGER(KIND=JPIM) :: IATM, JMOL, IXMAX, J1, J2, JK, JL
INTEGER(KIND=JPIM) :: ITMOL, INXMOL

REAL(KIND=JPRB) :: ZAMM

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ***

! *** mji
! Initialize all molecular amounts to zero here, 
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

IF (LHOOK) CALL DR_HOOK('RRTM_PREPARE_GASES',0,ZHOOK_HANDLE)

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

DO JL = KIDIA, KFDIA
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

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('RRTM_PREPARE_GASES',1,ZHOOK_HANDLE)

END SUBROUTINE RRTM_PREPARE_GASES
