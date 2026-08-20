! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE AER_CGROWTH_H13 &
  &( KIDIA, KFDIA, KLON, KLEV, &
  &  PAERPHO, PITAERPHO, PTSPHY, &
  &  PQP, PRHO, & 
  &  POHOK, PO3OK, & 
  &  PTAERPHI, PTAERPHO,PRATEOXD, PRATECC )

!*** * AER_CGROWTH_H13* - GROWTH OF O.M. AND B.C. AEROSOLS following Huang13
!parameterization

!**   INTERFACE.
!     ----------
!          *AER_CGROWTH_H13* IS CALLED FROM *AER_PHY3*.
! INPUTS:
! -------
! PAERPHO (KLON,KLEV)  : HYDROPHOBIC SPECIES (OM/BC) MMR        (kg/kg)
! PITAERPHO (KLON,KLEV)  : HYDROPHOBIC SPECIES (OM/BC) TENDENCY   (kg/kg.s-1)
! PQP(KLON,KLEV)       : FULL-LEVEL HUMIDITY (W. DYN.TEND.) (kg kg-1)
! PRHO(KLON,KLEV)      : AIR DENSITY (kg m-3)
! POHOK(KLON,KLEV)     : OH MMR (kg/kg)
! PO3OK(KLON,KLEV)     : O3 MMR (kg/kg)
!
! OUTUTS:
! -------
! PTAERPHI (KLON,KLEV)  : HYDROPHILIC SPECIES (OM/BC) TENDENCY   (kg/kg.s-1)
! PTAERPHO (KLON,KLEV)  : HYDROPHOBIC SPECIES (OM/BC) TENDENCY   (kg/kg.s-1)
! PRATEOXD (KLON,KLEV)  : OXIDATION RATE    (s-1)
! PRATECC (KLON,KLEV)   : COAGULATION/CONDENSATION RATE    (s-1)


!     AUTHOR.
!     -------
!        THIERRY ELIAS AND SAMUEL REMY

!     SOURCE.
!     -------

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 2022-03-22
!-----------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMCST    , ONLY : RG,  RTT,  RPI,  RNAVO,  RMV

IMPLICIT NONE

!-----------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV

REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSPHY
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAERPHO(KLON,KLEV) ,PITAERPHO(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQP(KLON,KLEV), PRHO(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: POHOK(KLON,KLEV), PO3OK(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTAERPHI(KLON,KLEV), PTAERPHO(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRATEOXD(KLON,KLEV),PRATECC(KLON,KLEV)


!*       0.5   LOCAL VARIABLES
!              ---------------

INTEGER(KIND=JPIM) :: JL, JK
REAL(KIND=JPRB) :: ZCOEF, ZAERCONV
REAL(KIND=JPRB) :: ZFCT_SH
REAL(KIND=JPRB) :: ZO3_MOLEC(KLON,KLEV), ZOH_MOLEC(KLON,KLEV),ZH2O_MOLEC(KLON,KLEV)
REAL(KIND=JPRB) :: ZK_CC, ZALPHA, ZBETA
REAL(KIND=JPRB) :: ZK_OXD, ZSHIELD, ZNUM, ZDEN, ZK_INF, ZK_O3, ZK_H2O
REAL(KIND=JPRB), PARAMETER :: ZMWO3  = 48.E-3_JPRB
REAL(KIND=JPRB), PARAMETER :: ZMWOH  = 17.E-3_JPRB


REAL(KIND=JPHOOK) :: ZHOOK_HANDLE


!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('AER_CGROWTH_H13',0,ZHOOK_HANDLE)

! param from Huang et al. 2013:
ZALPHA = 5.8E-7_JPRB
ZBETA = 4.6E-12_JPRB

ZK_INF = 0.015_JPRB     ! s-1
ZK_O3 = 2.8e-13_JPRB    !   cm3
ZK_H2O = 2.1e-17_JPRB    !   cm3
ZSHIELD = 0.01_JPRB   !   = lambda He et al 2016 Eq. 15
ZK_OXD=0._JPRB
ZK_CC=0._JPRB

DO JK=1,KLEV
  DO JL=KIDIA,KFDIA

    ZOH_MOLEC(JL,JK) = POHOK(JL,JK) * PRHO(JL,JK) / ZMWOH * RNAVO * 1E-6_JPRB    !  molecules cm-3
    ZO3_MOLEC(JL,JK) = PO3OK(JL,JK) * PRHO(JL,JK) / ZMWO3 * RNAVO * 1E-6_JPRB    !  molecules cm-3
    IF (PQP(JL,JK) > 1.E-9_JPRB) THEN
      ZFCT_SH = 1._JPRB/((1._JPRB/PQP(JL,JK)) - 1._JPRB)
    ELSE
      ZFCT_SH=0._JPRB
    ENDIF
    ZH2O_MOLEC(JL,JK) = ZFCT_SH * PRHO(JL,JK) * RNAVO / RMV * 1E-6_JPRB !molecules / cm3

    
    
    ZK_CC = ZBETA * ZOH_MOLEC(JL,JK) + ZALPHA   !  Huang et al. 2013 Eq. 6

    ZNUM = ZK_INF * ZK_O3  * ZO3_MOLEC(JL,JK)    !  Huang et al. 2013 Eq. 1
    ZDEN = 1._JPRB + ZK_O3 * ZO3_MOLEC(JL,JK) + ZK_H2O * ZH2O_MOLEC(JL,JK)
    IF (ZDEN > 0._JPRB) THEN
       ZK_OXD = ZSHIELD * ZNUM / ZDEN   ! s-1   !  He et al 2016 Eq. 15
    ENDIF

    ZCOEF = ZK_OXD + ZK_CC      !  Huang et al. Eq. 8, He et al 2016 Eq. 17
    PRATEOXD(JL,JK)=ZK_OXD
    PRATECC(JL,JK)=ZK_CC

    !   not changed from AER_CGROWTH
    ZAERCONV = (PAERPHO(JL,JK) + PTSPHY * PITAERPHO(JL,JK)) * ZCOEF
    PTAERPHI(JL,JK) =  ZAERCONV     !   outputs
    PTAERPHO(JL,JK) = -ZAERCONV    
  ENDDO
ENDDO
!-----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('AER_CGROWTH_H13',1,ZHOOK_HANDLE)

END SUBROUTINE AER_CGROWTH_H13
