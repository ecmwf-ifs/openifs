! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE CUENTR &
 & ( YDECUMF, YDSPP_CONFIG,&
 & KIDIA,     KFDIA,    KLON,     KLEV,&
 & KK,        KCBOT,    KTYPE,&
 & LDCUM,     LDWORK,&
 & PQSEN,     PAPH,     PGEOH,&
 & PMFU,      PGP2DSPP, PDMFEN,   PDMFDE )

!          M.TIEDTKE         E.C.M.W.F.     12/89
!          P.BECHTOLD        E.C.M.W.F.     06/07
!          P.BECHTOLD        E.C.M.W.F.     03/12

!          PURPOSE.
!          --------
!          THIS ROUTINE CALCULATES ENTRAINMENT/DETRAINMENT RATES
!          FOR UPDRAFTS IN CUMULUS PARAMETERIZATION

!          INTERFACE
!          ---------

!          THIS ROUTINE IS CALLED FROM *CUASC*.
!          INPUT ARE ENVIRONMENTAL VALUES T,Q ETC
!          AND UPDRAFT VALUES T,Q ETC
!          IT RETURNS ENTRAINMENT/DETRAINMENT RATES

!          METHOD.
!          --------
!          TURBULENT ENTRAINMENT IS SIMULATED BY A CONSTANT
!          MULTIPLIED BY A VERTICAL SCALING FUNCTION

!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):

!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KLEV*         NUMBER OF LEVELS
!    *KK*           CURRENT LEVEL
!    *KCBOT*        CLOUD BASE LEVEL

!    INPUT PARAMETERS (LOGICAL):

!    *LDCUM*        FLAG: .TRUE. FOR CONVECTIVE POINTS

!    INPUT PARAMETERS (REAL):

!    *PQSEN*        SATURATION SPEC. HUMIDITY                   KG/KG
!    *PAPH*         PROVISIONAL PRESSURE ON HALF LEVELS          PA
!    *PGEOH*        PROVISIONAL GEOPOTENTIAL ON HALF LEVELS      PA
!    *PMFU*         MASSFLUX IN UPDRAFTS                        KG/(M2*S)
!    *PGP2DSPP*     Standard stochastic variable (mean=0, SD=1)

!    OUTPUT PARAMETERS (REAL):

!    *PDMFEN*       ENTRAINMENT RATE                            KG/(M2*S)
!    *PDMFDE*       DETRAINMENT RATE                            KG/(M2*S)

!          EXTERNALS
!          ---------
!          NONE
!
!
!    MODIFICATIONS.
!    --------------
!
!    M. Leutbecher & S.-J. Lock (Jan 2016) Introduced SPP scheme (LSPP)
!    M. Leutbecher (Oct 2020) SPP abstraction
!----------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMCST   , ONLY : RG
USE YOECUMF  , ONLY : TECUMF
USE SPP_MOD     , ONLY : TSPP_CONFIG
USE SPP_GEN_MOD , ONLY : SPP_PERT

IMPLICIT NONE

TYPE(TECUMF)      ,INTENT(INOUT) :: YDECUMF
TYPE(TSPP_CONFIG) ,INTENT(IN)    :: YDSPP_CONFIG
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KK 
INTEGER(KIND=JPIM),INTENT(IN)    :: KCBOT(KLON) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTYPE(KLON) 
LOGICAL           ,INTENT(IN)    :: LDCUM(KLON) 
LOGICAL           ,INTENT(IN)    :: LDWORK 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQSEN(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPH(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOH(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMFU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGP2DSPP(KLON,YDSPP_CONFIG%SM%NRFTOTAL)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDMFEN(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDMFDE(KLON) 

LOGICAL ::  LLO1

INTEGER(KIND=JPIM) :: JL

! A bunch of SPP variables that are saved between calls for efficiency reasons
LOGICAL            :: LLPERT_DETRPEN     ! SPP perturbation on?
INTEGER(KIND=JPIM) :: IPDETRPEN ! SPP random field pointer
INTEGER(KIND=JPIM) :: IPN  ! SPP perturbation pointer
TYPE(SPP_PERT)     :: PN1  ! SPP pertn. config. for RTAU

REAL(KIND=JPRB) :: ZDZ, ZENTR(KLON), ZMF, ZRG, ZXDETRPEN
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!----------------------------------------------------------------------

!*    1.           CALCULATE ENTRAINMENT AND DETRAINMENT RATES
!                  -------------------------------------------

IF (LHOOK) CALL DR_HOOK('CUENTR',0,ZHOOK_HANDLE)
ASSOCIATE(DETRPEN=>YDECUMF%DETRPEN)
IF(LDWORK) THEN

  ZRG=1.0_JPRB/RG

  ! prepare SPP perturbations (just once)
  IF (YDSPP_CONFIG%LSPP) THEN
    IPN = YDSPP_CONFIG%PPTR%DETRPEN
    LLPERT_DETRPEN = IPN > 0
    IF (LLPERT_DETRPEN) THEN
      PN1=YDSPP_CONFIG%SM%PN(IPN)
      IPDETRPEN = PN1%MP
    ENDIF
  ELSE
    LLPERT_DETRPEN=.FALSE.
  ENDIF

  DO JL=KIDIA,KFDIA
    PDMFEN(JL)=0.0_JPRB
    PDMFDE(JL)=0.0_JPRB
    ZENTR(JL)=0.0
  ENDDO

!*    1.1          SPECIFY ENTRAINMENT RATES
!                  -------------------------

  DO JL=KIDIA,KFDIA
    IF(LDCUM(JL)) THEN
      ZDZ=(PGEOH(JL,KK)-PGEOH(JL,KK+1))*ZRG
      ZMF=PMFU(JL,KK+1)*ZDZ
      LLO1=KK < KCBOT(JL)
      IF(LLO1) THEN
        IF (LLPERT_DETRPEN) THEN
          ZXDETRPEN=DETRPEN*EXP(PN1%MU(1)+PN1%XMAG(1)*PGP2DSPP(JL, IPDETRPEN))
        ELSE
          ZXDETRPEN=DETRPEN
        ENDIF
        PDMFEN(JL)=ZENTR(JL)*ZMF
        PDMFDE(JL)=ZXDETRPEN*ZMF
      ENDIF
    ENDIF
  ENDDO

ENDIF

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CUENTR',1,ZHOOK_HANDLE)
END SUBROUTINE CUENTR
