SUBROUTINE SRTM_SETCOEF&
 & (YDDIMV, KIDIA   , KFDIA    , KLEV    ,&
 & PAVEL   , PTAVEL   ,&
 & PCOLDRY , PWKL     ,&
 & KLAYTROP,&
 & PCOLCH4  , PCOLCO2 , PCOLH2O , PCOLMOL  , PCOLO2 , PCOLO3 ,&
 & PFORFAC , PFORFRAC , KINDFOR , PSELFFAC, PSELFFRAC, KINDSELF ,&
 & PFAC00  , PFAC01   , PFAC10  , PFAC11  ,&
 & KJP     , KJT      , KJT1    , PRMU0&
 & )  

!     J. Delamere, AER, Inc. (version 2.5, 02/04/01)

!     Modifications:
!     JJMorcrette 030224   rewritten / adapted to ECMWF F90 system
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        D.Salmond  31-Oct-2007 Vector version in the style of RRTM from Meteo France & NEC

!     Purpose:  For a given atmosphere, calculate the indices and
!     fractions related to the pressure and temperature interpolations.

USE YOMDIMV  , ONLY : TDIMV
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOESRTWN , ONLY : PREFLOG, TREF
!!  USE YOESWN  , ONLY : NDBUG

IMPLICIT NONE

!-- Input arguments

TYPE(TDIMV)       ,INTENT(IN)    :: YDDIMV
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA, KFDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAVEL(KIDIA:KFDIA,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTAVEL(KIDIA:KFDIA,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCOLDRY(KIDIA:KFDIA,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWKL(KIDIA:KFDIA,35,YDDIMV%NFLEVG) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KLAYTROP(KIDIA:KFDIA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCOLCH4(KIDIA:KFDIA,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCOLCO2(KIDIA:KFDIA,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCOLH2O(KIDIA:KFDIA,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCOLMOL(KIDIA:KFDIA,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCOLO2(KIDIA:KFDIA,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCOLO3(KIDIA:KFDIA,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFORFAC(KIDIA:KFDIA,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFORFRAC(KIDIA:KFDIA,YDDIMV%NFLEVG) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KINDFOR(KIDIA:KFDIA,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSELFFAC(KIDIA:KFDIA,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSELFFRAC(KIDIA:KFDIA,YDDIMV%NFLEVG) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KINDSELF(KIDIA:KFDIA,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFAC00(KIDIA:KFDIA,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFAC01(KIDIA:KFDIA,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFAC10(KIDIA:KFDIA,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFAC11(KIDIA:KFDIA,YDDIMV%NFLEVG) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KJP(KIDIA:KFDIA,YDDIMV%NFLEVG) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KJT(KIDIA:KFDIA,YDDIMV%NFLEVG) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KJT1(KIDIA:KFDIA,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRMU0(KIDIA:KFDIA) 
!-- Output arguments

!-- local integers

INTEGER(KIND=JPIM) :: I_NLAYERS, JK, JL, JP1

!-- local reals

REAL(KIND=JPRB) :: Z_STPFAC, Z_PLOG
REAL(KIND=JPRB) :: Z_FP, Z_FT, Z_FT1, Z_WATER, Z_SCALEFAC
REAL(KIND=JPRB) :: Z_FACTOR, Z_CO2REG, Z_COMPFP
!REAL(KIND=JPRB) :: Z_TBNDFRAC, Z_T0FRAC
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE




IF (LHOOK) CALL DR_HOOK('SRTM_SETCOEF',0,ZHOOK_HANDLE)
ASSOCIATE(NFLEVG=>YDDIMV%NFLEVG)

Z_STPFAC = 296._JPRB/1013._JPRB
I_NLAYERS = KLEV

DO JL = KIDIA, KFDIA
  IF (PRMU0(JL) > 0.0_JPRB) THEN
    KLAYTROP(JL)  = 0
  ENDIF
ENDDO

DO JK = 1, I_NLAYERS
  DO JL = KIDIA, KFDIA
    IF (PRMU0(JL) > 0.0_JPRB) THEN
      !        Find the two reference pressures on either side of the
      !        layer pressure.  Store them in JP and JP1.  Store in FP the
      !        fraction of the difference (in ln(pressure)) between these
      !        two values that the layer pressure lies.

      Z_PLOG = LOG(PAVEL(JL,JK))
      KJP(JL,JK) = INT(36. - 5*(Z_PLOG+0.04))
      IF (KJP(JL,JK) < 1) THEN
        KJP(JL,JK) = 1
      ELSEIF (KJP(JL,JK) > 58) THEN
        KJP(JL,JK) = 58
      ENDIF
      JP1 = KJP(JL,JK) + 1
      Z_FP = 5. * (PREFLOG(KJP(JL,JK)) - Z_PLOG)

      !        Determine, for each reference pressure (JP and JP1), which
      !        reference temperature (these are different for each  
      !        reference pressure) is nearest the layer temperature but does
      !        not exceed it.  Store these indices in JT and JT1, resp.
      !        Store in FT (resp. FT1) the fraction of the way between JT
      !        (JT1) and the next highest reference temperature that the 
      !        layer temperature falls.

      KJT(JL,JK) = INT(3. + (PTAVEL(JL,JK)-TREF(KJP(JL,JK)))/15.)
      IF (KJT(JL,JK) < 1) THEN
        KJT(JL,JK) = 1
      ELSEIF (KJT(JL,JK) > 4) THEN
        KJT(JL,JK) = 4
      ENDIF
      Z_FT = ((PTAVEL(JL,JK)-TREF(KJP(JL,JK)))/15.) - REAL(KJT(JL,JK)-3)
      KJT1(JL,JK) = INT(3. + (PTAVEL(JL,JK)-TREF(JP1))/15.)
      IF (KJT1(JL,JK) < 1) THEN
        KJT1(JL,JK) = 1
      ELSEIF (KJT1(JL,JK) > 4) THEN
        KJT1(JL,JK) = 4
      ENDIF
      Z_FT1 = ((PTAVEL(JL,JK)-TREF(JP1))/15.) - REAL(KJT1(JL,JK)-3)

      Z_WATER = PWKL(JL,1,JK)/PCOLDRY(JL,JK)
      Z_SCALEFAC = PAVEL(JL,JK) * Z_STPFAC / PTAVEL(JL,JK)

      !        If the pressure is less than ~100mb, perform a different
      !        set of species interpolations.

!      IF (Z_PLOG <= 4.56) GO TO 5300
      IF (KJP(JL,JK) >= 13) GO TO 5300
      KLAYTROP(JL) =  KLAYTROP(JL) + 1

      !        Set up factors needed to separately include the water vapor
      !        foreign-continuum in the calculation of absorption coefficient.

      PFORFAC(JL,JK) = Z_SCALEFAC / (1.+Z_WATER)
      Z_FACTOR = (332.0-PTAVEL(JL,JK))/36.0
      KINDFOR(JL,JK) = MIN(2, MAX(1, INT(Z_FACTOR)))
      PFORFRAC(JL,JK) = Z_FACTOR - REAL(KINDFOR(JL,JK))

      !        Set up factors needed to separately include the water vapor
      !        self-continuum in the calculation of absorption coefficient.

      PSELFFAC(JL,JK) = Z_WATER * PFORFAC(JL,JK)
      Z_FACTOR = (PTAVEL(JL,JK)-188.0)/7.2
      KINDSELF(JL,JK) = MIN(9, MAX(1, INT(Z_FACTOR)-7))
      PSELFFRAC(JL,JK) = Z_FACTOR - REAL(KINDSELF(JL,JK) + 7)

      !        Calculate needed column amounts.

      PCOLH2O(JL,JK) = 1.E-20 * PWKL(JL,1,JK)
      PCOLCO2(JL,JK) = 1.E-20 * PWKL(JL,2,JK)
      PCOLO3(JL,JK) = 1.E-20 * PWKL(JL,3,JK)
      !         COLO3(LAY) = 0.
      !         COLO3(LAY) = colo3(lay)/1.16
      PCOLCH4(JL,JK) = 1.E-20 * PWKL(JL,6,JK)
      PCOLO2(JL,JK) = 1.E-20 * PWKL(JL,7,JK)
      PCOLMOL(JL,JK) = 1.E-20 * PCOLDRY(JL,JK) + PCOLH2O(JL,JK)
      !         colco2(lay) = 0.
      !         colo3(lay) = 0.
      !         colch4(lay) = 0.
      !         colo2(lay) = 0.
      !         colmol(lay) = 0.
      IF (PCOLCO2(JL,JK) == 0.) PCOLCO2(JL,JK) = 1.E-32 * PCOLDRY(JL,JK)
      IF (PCOLCH4(JL,JK) == 0.) PCOLCH4(JL,JK) = 1.E-32 * PCOLDRY(JL,JK)
      IF (PCOLO2(JL,JK) == 0.) PCOLO2(JL,JK) = 1.E-32 * PCOLDRY(JL,JK)
      !        Using E = 1334.2 cm-1.
      Z_CO2REG = 3.55E-24 * PCOLDRY(JL,JK)
      GO TO 5400

      !        Above LAYTROP.
5300  CONTINUE

      !        Set up factors needed to separately include the water vapor
      !        foreign-continuum in the calculation of absorption coefficient.

      PFORFAC(JL,JK) = Z_SCALEFAC / (1.+Z_WATER)
      Z_FACTOR = (PTAVEL(JL,JK)-188.0)/36.0
      KINDFOR(JL,JK) = 3
      PFORFRAC(JL,JK) = Z_FACTOR - 1.0

      !        Calculate needed column amounts.

      PCOLH2O(JL,JK) = 1.E-20 * PWKL(JL,1,JK)
      PCOLCO2(JL,JK) = 1.E-20 * PWKL(JL,2,JK)
      PCOLO3(JL,JK)  = 1.E-20 * PWKL(JL,3,JK)
      PCOLCH4(JL,JK) = 1.E-20 * PWKL(JL,6,JK)
      PCOLO2(JL,JK)  = 1.E-20 * PWKL(JL,7,JK)
      PCOLMOL(JL,JK) = 1.E-20 * PCOLDRY(JL,JK) + PCOLH2O(JL,JK)
      IF (PCOLCO2(JL,JK) == 0.) PCOLCO2(JL,JK) = 1.E-32 * PCOLDRY(JL,JK)
      IF (PCOLCH4(JL,JK) == 0.) PCOLCH4(JL,JK) = 1.E-32 * PCOLDRY(JL,JK)
      IF (PCOLO2(JL,JK) == 0.) PCOLO2(JL,JK)  = 1.E-32 * PCOLDRY(JL,JK)
      Z_CO2REG = 3.55E-24 * PCOLDRY(JL,JK)

      PSELFFAC(JL,JK) =0.0_JPRB
      PSELFFRAC(JL,JK)=0.0_JPRB
      KINDSELF(JL,JK) = 0

5400  CONTINUE

      !        We have now isolated the layer ln pressure and temperature,
      !        between two reference pressures and two reference temperatures 
      !        (for each reference pressure).  We multiply the pressure 
      !        fraction FP with the appropriate temperature fractions to get 
      !        the factors that will be needed for the interpolation that yields
      !        the optical depths (performed in routines TAUGBn for band n).

      Z_COMPFP = 1. - Z_FP
      PFAC10(JL,JK) = Z_COMPFP * Z_FT
      PFAC00(JL,JK) = Z_COMPFP * (1. - Z_FT)
      PFAC11(JL,JK) = Z_FP * Z_FT1
      PFAC01(JL,JK) = Z_FP * (1. - Z_FT1)

      !  IF (NDBUG.LE.3) THEN
      !    print 9000,LAY,LAYTROP,JP(LAY),JT(LAY),JT1(LAY),TAVEL(LAY) &
      !      &,FAC00(LAY),FAC01(LAY),FAC10(LAY),FAC11(LAY) &
      !      &,COLMOL(LAY),COLCH4(LAY),COLCO2(LAY),COLH2O(LAY) &
      !      &,COLO2(LAY),COLO3(LAY),SELFFAC(LAY),SELFFRAC(LAY) &
      !      &,FORFAC(LAY),FORFRAC(LAY),INDSELF(LAY),INDFOR(LAY)
9000  FORMAT(1X,2I3,3I4,F6.1,4F7.2,12E9.2,2I5)
      !  ENDIF

    ENDIF
  ENDDO
ENDDO

!----------------------------------------------------------------------- 
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SRTM_SETCOEF',1,ZHOOK_HANDLE)
END SUBROUTINE SRTM_SETCOEF
