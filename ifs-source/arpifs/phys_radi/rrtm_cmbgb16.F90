!***************************************************************************
SUBROUTINE RRTM_CMBGB16
!***************************************************************************

!     BAND 16:  2600-3000 cm-1 (low - H2O,CH4; high - nothing)
!     ABozzo 201306 updated to rrtmg v4.85
!     band 16:  2600-3250 cm-1 (low key- h2o,ch4; high key - ch4)
!***************************************************************************

! Parameters
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOERRTO16, ONLY : KAO,KBO   ,SELFREFO,FORREFO  ,FRACREFAO,FRACREFBO
USE YOERRTA16, ONLY : KA,KB     ,SELFREF,FORREF    ,FRACREFA,FRACREFB
USE YOERRTRWT, ONLY : RWGT
USE YOERRTFTR, ONLY : NGC      ,NGS      ,NGN      

IMPLICIT NONE

INTEGER(KIND=JPIM) :: IGC, IPR, IPRSM, JN, JP, JT

REAL(KIND=JPRB) :: Z_SUMF, Z_SUMK
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('RRTM_CMBGB16',0,ZHOOK_HANDLE)
DO JN = 1,9
  DO JT = 1,5
    DO JP = 1,13
      IPRSM = 0
      DO IGC = 1,NGC(16)
        Z_SUMK = 0.0_JPRB
        DO IPR = 1, NGN(NGS(15)+IGC)
          IPRSM = IPRSM + 1

          Z_SUMK = Z_SUMK + KAO(JN,JT,JP,IPRSM)*RWGT(IPRSM+240)
        ENDDO

        KA(JN,JT,JP,IGC) = Z_SUMK
      ENDDO
    ENDDO
  ENDDO
ENDDO

DO JT = 1,5
   DO JP = 13,59
      IPRSM = 0
      DO IGC = 1,NGC(16)
         Z_SUMK = 0.0_JPRB
         DO IPR = 1, NGN(NGS(15)+IGC)
            IPRSM = IPRSM + 1
            Z_SUMK = Z_SUMK + KBO(JT,JP,IPRSM)*RWGT(IPRSM+240)
         ENDDO
         KB(JT,JP,IGC) = Z_SUMK
      ENDDO
   ENDDO
ENDDO


DO JT = 1,10
  IPRSM = 0
  DO IGC = 1,NGC(16)
    Z_SUMK = 0.0_JPRB
    DO IPR = 1, NGN(NGS(15)+IGC)
      IPRSM = IPRSM + 1

      Z_SUMK = Z_SUMK + SELFREFO(JT,IPRSM)*RWGT(IPRSM+240)
    ENDDO

    SELFREF(JT,IGC) = Z_SUMK
  ENDDO
ENDDO

DO JT = 1,4
   IPRSM = 0
   DO IGC = 1,NGC(16)
      Z_SUMK = 0.0_JPRB
      DO IPR = 1, NGN(NGS(15)+IGC)
         IPRSM = IPRSM + 1
         Z_SUMK = Z_SUMK + FORREFO(JT,IPRSM)*RWGT(IPRSM+240)
      ENDDO
      FORREF(JT,IGC) = Z_SUMK
   ENDDO
ENDDO

IPRSM = 0
DO IGC = 1,NGC(16)
   Z_SUMF = 0.0_JPRB
   DO IPR = 1, NGN(NGS(15)+IGC)
      IPRSM = IPRSM + 1
      Z_SUMF = Z_SUMF + FRACREFBO(IPRSM)
   ENDDO
   FRACREFB(IGC) = Z_SUMF
ENDDO

DO JP = 1,9
  IPRSM = 0
  DO IGC = 1,NGC(16)
    Z_SUMF = 0.0_JPRB
    DO IPR = 1, NGN(NGS(15)+IGC)
      IPRSM = IPRSM + 1

      Z_SUMF = Z_SUMF + FRACREFAO(IPRSM,JP)
    ENDDO

    FRACREFA(IGC,JP) = Z_SUMF
  ENDDO
ENDDO

 
IF (LHOOK) CALL DR_HOOK('RRTM_CMBGB16',1,ZHOOK_HANDLE)
END SUBROUTINE RRTM_CMBGB16
