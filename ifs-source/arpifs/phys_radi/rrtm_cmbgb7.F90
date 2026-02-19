!***************************************************************************
SUBROUTINE RRTM_CMBGB7
!***************************************************************************

!     BAND 7:  980-1080 cm-1 (low - H2O,O3; high - O3)
!     ABozzo 201306 updated to rrtmg v4.85
!     band 7:  980-1080 cm-1 (low key - h2o,o3; low minor - co2)
!                            (high key - o3; high minor - co2)
!***************************************************************************

! Parameters
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOERRTO7 , ONLY : KAO     ,KBO     ,SELFREFO   ,FORREFO, FRACREFAO  ,&
 & FRACREFBO,  KAO_MCO2     ,KBO_MCO2  
USE YOERRTA7 , ONLY : KA      ,KB      ,SELFREF    ,FORREF, FRACREFA   ,&
 & FRACREFB,  KA_MCO2     ,KB_MCO2     
USE YOERRTRWT, ONLY : RWGT
USE YOERRTFTR, ONLY : NGC      ,NGS      ,NGN      

IMPLICIT NONE

INTEGER(KIND=JPIM) :: IGC, IPR, IPRSM, JN, JP, JT

REAL(KIND=JPRB) :: Z_SUMF, Z_SUMK
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('RRTM_CMBGB7',0,ZHOOK_HANDLE)
DO JN = 1,9
  DO JT = 1,5
    DO JP = 1,13
      IPRSM = 0
      DO IGC = 1,NGC(7)
        Z_SUMK = 0.0_JPRB
        DO IPR = 1, NGN(NGS(6)+IGC)
          IPRSM = IPRSM + 1

          Z_SUMK = Z_SUMK + KAO(JN,JT,JP,IPRSM)*RWGT(IPRSM+96)
        ENDDO

        KA(JN,JT,JP,IGC) = Z_SUMK
      ENDDO
    ENDDO
  ENDDO
ENDDO
DO JT = 1,5
  DO JP = 13,59
    IPRSM = 0
    DO IGC = 1,NGC(7)
      Z_SUMK = 0.0_JPRB
      DO IPR = 1, NGN(NGS(6)+IGC)
        IPRSM = IPRSM + 1

        Z_SUMK = Z_SUMK + KBO(JT,JP,IPRSM)*RWGT(IPRSM+96)
      ENDDO

      KB(JT,JP,IGC) = Z_SUMK
    ENDDO
  ENDDO
ENDDO

DO JN = 1,9
      DO JT = 1,19
         IPRSM = 0
         DO IGC = 1,NGC(7)
            Z_SUMK = 0.0_JPRB
            DO IPR = 1, NGN(NGS(6)+IGC)
               IPRSM = IPRSM + 1
               Z_SUMK = Z_SUMK + KAO_MCO2(JN,JT,IPRSM)*RWGT(IPRSM+96)
            ENDDO
            KA_MCO2(JN,JT,IGC) = Z_SUMK
         ENDDO
       ENDDO
ENDDO

DO JT = 1,19
      IPRSM = 0
      DO IGC = 1,NGC(7)
         Z_SUMK = 0.0_JPRB
         DO IPR = 1, NGN(NGS(6)+IGC)
            IPRSM = IPRSM + 1
            Z_SUMK = Z_SUMK + KBO_MCO2(JT,IPRSM)*RWGT(IPRSM+96)
         ENDDO
         KB_MCO2(JT,IGC) = Z_SUMK
      ENDDO
ENDDO


DO JT = 1,10
  IPRSM = 0
  DO IGC = 1,NGC(7)
    Z_SUMK = 0.0_JPRB
    DO IPR = 1, NGN(NGS(6)+IGC)
      IPRSM = IPRSM + 1

      Z_SUMK = Z_SUMK + SELFREFO(JT,IPRSM)*RWGT(IPRSM+96)
    ENDDO

    SELFREF(JT,IGC) = Z_SUMK
  ENDDO
ENDDO

DO JT = 1,4
      IPRSM = 0
      DO IGC = 1,NGC(7)
         Z_SUMK = 0.0_JPRB
         DO IPR = 1, NGN(NGS(6)+IGC)
            IPRSM = IPRSM + 1
            Z_SUMK = Z_SUMK + FORREFO(JT,IPRSM)*RWGT(IPRSM+96)
         ENDDO
         FORREF(JT,IGC) = Z_SUMK
      ENDDO
ENDDO

DO JP = 1,9
  IPRSM = 0
  DO IGC = 1,NGC(7)
    Z_SUMF = 0.0_JPRB
    DO IPR = 1, NGN(NGS(6)+IGC)
      IPRSM = IPRSM + 1

      Z_SUMF = Z_SUMF + FRACREFAO(IPRSM,JP)
    ENDDO

    FRACREFA(IGC,JP) = Z_SUMF
  ENDDO
ENDDO

IPRSM = 0
DO IGC = 1,NGC(7)
  Z_SUMF = 0.0_JPRB
  DO IPR = 1, NGN(NGS(6)+IGC)
    IPRSM = IPRSM + 1

    Z_SUMF = Z_SUMF + FRACREFBO(IPRSM)
  ENDDO

  FRACREFB(IGC) = Z_SUMF
ENDDO


IF (LHOOK) CALL DR_HOOK('RRTM_CMBGB7',1,ZHOOK_HANDLE)
END SUBROUTINE RRTM_CMBGB7
