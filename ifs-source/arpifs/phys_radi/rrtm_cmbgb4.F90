!***************************************************************************
SUBROUTINE RRTM_CMBGB4
!***************************************************************************

!     BAND 4:  630-700 cm-1 (low - H2O,CO2; high - O3,CO2)
!     ABozzo 201306 updated to rrtmg v4.85
!***************************************************************************

! Parameters
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOERRTO4 , ONLY : KAO     ,KBO     ,SELFREFO   , FORREFO, FRACREFAO  ,FRACREFBO
USE YOERRTA4 , ONLY : KA      ,KB      ,SELFREF    , FORREF, FRACREFA   ,FRACREFB
USE YOERRTRWT, ONLY : RWGT
USE YOERRTFTR, ONLY : NGC      ,NGS      ,NGN      

IMPLICIT NONE

INTEGER(KIND=JPIM) :: IGC, IPR, IPRSM, JN, JP, JT

REAL(KIND=JPRB) :: Z_SUMF, Z_SUMK
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('RRTM_CMBGB4',0,ZHOOK_HANDLE)
DO JN = 1,9
  DO JT = 1,5
    DO JP = 1,13
      IPRSM = 0
      DO IGC = 1,NGC(4)
        Z_SUMK = 0.0_JPRB
        DO IPR = 1, NGN(NGS(3)+IGC)
          IPRSM = IPRSM + 1

          Z_SUMK = Z_SUMK + KAO(JN,JT,JP,IPRSM)*RWGT(IPRSM+48)
        ENDDO

        KA(JN,JT,JP,IGC) = Z_SUMK
      ENDDO
    ENDDO
  ENDDO
ENDDO
DO JN = 1,5
  DO JT = 1,5
    DO JP = 13,59
      IPRSM = 0
      DO IGC = 1,NGC(4)
        Z_SUMK = 0.0_JPRB
        DO IPR = 1, NGN(NGS(3)+IGC)
          IPRSM = IPRSM + 1

          Z_SUMK = Z_SUMK + KBO(JN,JT,JP,IPRSM)*RWGT(IPRSM+48)
        ENDDO

        KB(JN,JT,JP,IGC) = Z_SUMK
      ENDDO
    ENDDO
  ENDDO
ENDDO

DO JT = 1,10
  IPRSM = 0
  DO IGC = 1,NGC(4)
    Z_SUMK = 0.0_JPRB
    DO IPR = 1, NGN(NGS(3)+IGC)
      IPRSM = IPRSM + 1

      Z_SUMK = Z_SUMK + SELFREFO(JT,IPRSM)*RWGT(IPRSM+48)
    ENDDO

    SELFREF(JT,IGC) = Z_SUMK
  ENDDO
ENDDO

DO JT = 1,4
   IPRSM = 0
   DO IGC = 1,NGC(4)
     Z_SUMK = 0.0_JPRB
     DO IPR = 1, NGN(NGS(3)+IGC)
       IPRSM = IPRSM + 1
       Z_SUMK = Z_SUMK + FORREFO(JT,IPRSM)*RWGT(IPRSM+48)
     ENDDO
     FORREF(JT,IGC) = Z_SUMK
   ENDDO
ENDDO

DO JP = 1,9
  IPRSM = 0
  DO IGC = 1,NGC(4)
    Z_SUMF = 0.0_JPRB
    DO IPR = 1, NGN(NGS(3)+IGC)
      IPRSM = IPRSM + 1

      Z_SUMF = Z_SUMF + FRACREFAO(IPRSM,JP)
    ENDDO

    FRACREFA(IGC,JP) = Z_SUMF
  ENDDO
ENDDO

DO JP = 1,5
  IPRSM = 0
  DO IGC = 1,NGC(4)
    Z_SUMF = 0.0_JPRB
    DO IPR = 1, NGN(NGS(3)+IGC)
      IPRSM = IPRSM + 1

      Z_SUMF = Z_SUMF + FRACREFBO(IPRSM,JP)
    ENDDO

    FRACREFB(IGC,JP) = Z_SUMF
  ENDDO
ENDDO


IF (LHOOK) CALL DR_HOOK('RRTM_CMBGB4',1,ZHOOK_HANDLE)
END SUBROUTINE RRTM_CMBGB4
