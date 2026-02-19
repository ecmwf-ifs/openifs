!***************************************************************************
SUBROUTINE RRTM_CMBGB13
!***************************************************************************

!     BAND 13:  2080-2250 cm-1 (low - H2O,N2O; high - nothing)
!     ABozzo 201306 updated to rrtmg v4.85
!     band 13:  2080-2250 cm-1 (low key - h2o,n2o; high minor - o3 minor)
!***************************************************************************

! Parameters
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOERRTO13, ONLY : KAO     ,SELFREFO, FORREFO   ,FRACREFAO, FRACREFBO, &
                     & KAO_MCO2, KAO_MCO, KBO_MO3
USE YOERRTA13, ONLY : KA      ,SELFREF, FORREF    ,FRACREFA, FRACREFB, &
                     & KA_MCO2, KA_MCO, KB_MO3
USE YOERRTRWT, ONLY : RWGT
USE YOERRTFTR, ONLY : NGC      ,NGS      ,NGN      

IMPLICIT NONE

INTEGER(KIND=JPIM) :: IGC, IPR, IPRSM, JN, JP, JT

REAL(KIND=JPRB) :: Z_SUMF, Z_SUMK, Z_SUMK1, Z_SUMK2
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE


IF (LHOOK) CALL DR_HOOK('RRTM_CMBGB13',0,ZHOOK_HANDLE)
DO JN = 1,9
  DO JT = 1,5
    DO JP = 1,13
      IPRSM = 0
      DO IGC = 1,NGC(13)
        Z_SUMK = 0.0_JPRB
        DO IPR = 1, NGN(NGS(12)+IGC)
          IPRSM = IPRSM + 1

          Z_SUMK = Z_SUMK + KAO(JN,JT,JP,IPRSM)*RWGT(IPRSM+192)
        ENDDO

        KA(JN,JT,JP,IGC) = Z_SUMK
      ENDDO
    ENDDO
  ENDDO
ENDDO

DO JN = 1,9
   DO JT = 1,19
      IPRSM = 0
      DO IGC = 1,NGC(13)
        Z_SUMK1 = 0.0_JPRB
        Z_SUMK2 = 0.0_JPRB
         DO IPR = 1, NGN(NGS(12)+IGC)
            IPRSM = IPRSM + 1
            Z_SUMK1 = Z_SUMK1 + KAO_MCO2(JN,JT,IPRSM)*RWGT(IPRSM+192)
            Z_SUMK2 = Z_SUMK2 + KAO_MCO(JN,JT,IPRSM)*RWGT(IPRSM+192)
         ENDDO
         KA_MCO2(JN,JT,IGC) = Z_SUMK1
         KA_MCO(JN,JT,IGC) = Z_SUMK2
      ENDDO
   ENDDO
ENDDO

DO JT = 1,19
   IPRSM = 0
   DO IGC = 1,NGC(13)
      Z_SUMK = 0.0_JPRB
      DO IPR = 1, NGN(NGS(12)+IGC)
         IPRSM = IPRSM + 1
         Z_SUMK = Z_SUMK + KBO_MO3(JT,IPRSM)*RWGT(IPRSM+192)
      ENDDO
      KB_MO3(JT,IGC) = Z_SUMK
   ENDDO
ENDDO


DO JT = 1,10
  IPRSM = 0
  DO IGC = 1,NGC(13)
    Z_SUMK = 0.0_JPRB
    DO IPR = 1, NGN(NGS(12)+IGC)
      IPRSM = IPRSM + 1

      Z_SUMK = Z_SUMK + SELFREFO(JT,IPRSM)*RWGT(IPRSM+192)
    ENDDO

    SELFREF(JT,IGC) = Z_SUMK
  ENDDO
ENDDO

DO JT = 1,4
   IPRSM = 0
   DO IGC = 1,NGC(13)
      Z_SUMK = 0.0_JPRB
      DO IPR = 1, NGN(NGS(12)+IGC)
         IPRSM = IPRSM + 1
         Z_SUMK = Z_SUMK + FORREFO(JT,IPRSM)*RWGT(IPRSM+192)
      ENDDO
      FORREF(JT,IGC) = Z_SUMK
   ENDDO
ENDDO

IPRSM = 0
DO IGC = 1,NGC(13)
   Z_SUMF = 0.0_JPRB
   DO IPR = 1, NGN(NGS(12)+IGC)
      IPRSM = IPRSM + 1
      Z_SUMF = Z_SUMF + FRACREFBO(IPRSM)
   ENDDO
   FRACREFB(IGC) = Z_SUMF
ENDDO

DO JP = 1,9
  IPRSM = 0
  DO IGC = 1,NGC(13)
    Z_SUMF = 0.0_JPRB
    DO IPR = 1, NGN(NGS(12)+IGC)
      IPRSM = IPRSM + 1

      Z_SUMF = Z_SUMF + FRACREFAO(IPRSM,JP)
    ENDDO

    FRACREFA(IGC,JP) = Z_SUMF
  ENDDO
ENDDO

IF (LHOOK) CALL DR_HOOK('RRTM_CMBGB13',1,ZHOOK_HANDLE)
END SUBROUTINE RRTM_CMBGB13
