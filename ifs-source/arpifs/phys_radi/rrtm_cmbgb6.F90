!***************************************************************************
SUBROUTINE RRTM_CMBGB6
!***************************************************************************

!     BAND 6:  820-980 cm-1 (low - H2O; high - nothing)
!     ABozzo 201306 updated to rrtmg v4.85
!     band 6:  820-980 cm-1 (low key - h2o; low minor - co2)
!                           (high key - nothing; high minor - cfc11, cfc12)
!***************************************************************************

! Parameters
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOERRTO6 , ONLY : KAO     ,SELFREFO   , FORREFO, FRACREFAO  ,&
 & KAO_MCO2 ,CFC11ADJO,CFC12O  
USE YOERRTA6 , ONLY : KA      ,SELFREF    , FORREF, FRACREFA   ,&
 & KA_MCO2  ,CFC11ADJ ,CFC12  
USE YOERRTRWT, ONLY : RWGT
USE YOERRTFTR, ONLY : NGC      ,NGS      ,NGN      

IMPLICIT NONE

INTEGER(KIND=JPIM) :: IGC, IPR, IPRSM, JP, JT

REAL(KIND=JPRB) :: Z_SUMF, Z_SUMK, Z_SUMK2, Z_SUMK3
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('RRTM_CMBGB6',0,ZHOOK_HANDLE)
DO JT = 1,5
  DO JP = 1,13
    IPRSM = 0
    DO IGC = 1,NGC(6)
      Z_SUMK = 0.0_JPRB
      DO IPR = 1, NGN(NGS(5)+IGC)
        IPRSM = IPRSM + 1

        Z_SUMK = Z_SUMK + KAO(JT,JP,IPRSM)*RWGT(IPRSM+80)
      ENDDO

      KA(JT,JP,IGC) = Z_SUMK
    ENDDO
  ENDDO
ENDDO

DO JT = 1,19
    IPRSM = 0
    DO IGC = 1,NGC(6)
        Z_SUMK = 0.0_JPRB
        DO IPR = 1, NGN(NGS(5)+IGC)
            IPRSM = IPRSM + 1
            Z_SUMK = Z_SUMK + KAO_MCO2(JT,IPRSM)*RWGT(IPRSM+80)
        ENDDO
        KA_MCO2(JT,IGC) = Z_SUMK
    ENDDO
ENDDO

DO JT = 1,10
  IPRSM = 0
  DO IGC = 1,NGC(6)
    Z_SUMK = 0.0_JPRB
    DO IPR = 1, NGN(NGS(5)+IGC)
      IPRSM = IPRSM + 1

      Z_SUMK = Z_SUMK + SELFREFO(JT,IPRSM)*RWGT(IPRSM+80)
    ENDDO

    SELFREF(JT,IGC) = Z_SUMK
  ENDDO
ENDDO

DO JT = 1,4
   IPRSM = 0
   DO IGC = 1,NGC(6)
      Z_SUMK = 0.0_JPRB
      DO IPR = 1, NGN(NGS(5)+IGC)
         IPRSM = IPRSM + 1
         Z_SUMK = Z_SUMK + FORREFO(JT,IPRSM)*RWGT(IPRSM+80)
      ENDDO
      FORREF(JT,IGC) = Z_SUMK
   ENDDO
ENDDO

IPRSM = 0
DO IGC = 1,NGC(6)
  Z_SUMF = 0.0_JPRB
  Z_SUMK2= 0.0_JPRB
  Z_SUMK3= 0.0_JPRB
  DO IPR = 1, NGN(NGS(5)+IGC)
    IPRSM = IPRSM + 1

    Z_SUMF = Z_SUMF + FRACREFAO(IPRSM)
    Z_SUMK2= Z_SUMK2+ CFC11ADJO(IPRSM)*RWGT(IPRSM+80)
    Z_SUMK3= Z_SUMK3+ CFC12O(IPRSM)*RWGT(IPRSM+80)
  ENDDO

  FRACREFA(IGC) = Z_SUMF
  CFC11ADJ(IGC) = Z_SUMK2
  CFC12(IGC) = Z_SUMK3
ENDDO


IF (LHOOK) CALL DR_HOOK('RRTM_CMBGB6',1,ZHOOK_HANDLE)
END SUBROUTINE RRTM_CMBGB6
