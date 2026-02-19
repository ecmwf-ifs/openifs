!***************************************************************************
SUBROUTINE RRTM_CMBGB3
!***************************************************************************

!     BAND 3:  500-630 cm-1 (low - H2O,CO2; high - H2O,CO2)
!      ABozzo 200130517 updated to rrtmg_lw_v4.85:
!     band 3:  500-630 cm-1 (low key - h2o,co2; low minor - n2o)
!                           (high key - h2o,co2; high minor - n2o)
!***************************************************************************

! Parameters
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOERRTO3 , ONLY : KAO     ,KBO     ,SELFREFO   ,FRACREFAO  ,&
 & FRACREFBO  ,FORREFO    ,KAO_MN2O   ,KBO_MN2O  
USE YOERRTA3 , ONLY : KA      ,KB      ,SELFREF    ,FRACREFA   ,&
 & FRACREFB   ,FORREF    ,KA_MN2O   ,KB_MN2O  
USE YOERRTRWT, ONLY : RWGT
USE YOERRTFTR, ONLY : NGC      ,NGS      ,NGN      

IMPLICIT NONE

INTEGER(KIND=JPIM) :: IGC, IPR, IPRSM, JN, JP, JT

REAL(KIND=JPRB) :: Z_SUMF, Z_SUMK
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('RRTM_CMBGB3',0,ZHOOK_HANDLE)
DO JN = 1,9
  DO JT = 1,5
    DO JP = 1,13
      IPRSM = 0
      DO IGC = 1,NGC(3)
        Z_SUMK = 0.0_JPRB
        DO IPR = 1, NGN(NGS(2)+IGC)
          IPRSM = IPRSM + 1

          Z_SUMK = Z_SUMK + KAO(JN,JT,JP,IPRSM)*RWGT(IPRSM+32)
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
      DO IGC = 1,NGC(3)
        Z_SUMK = 0.0_JPRB
        DO IPR = 1, NGN(NGS(2)+IGC)
          IPRSM = IPRSM + 1

          Z_SUMK = Z_SUMK + KBO(JN,JT,JP,IPRSM)*RWGT(IPRSM+32)
        ENDDO

        KB(JN,JT,JP,IGC) = Z_SUMK
      ENDDO
    ENDDO
  ENDDO
ENDDO

    DO JN = 1,9
         DO JT = 1,19
            IPRSM = 0
            DO IGC = 1,NGC(3)
              Z_SUMK = 0.
               DO IPR = 1, NGN(NGS(2)+IGC)
                  IPRSM = IPRSM + 1
                  Z_SUMK = Z_SUMK + KAO_MN2O(JN,JT,IPRSM)*RWGT(IPRSM+32)
               ENDDO
               KA_MN2O(JN,JT,IGC) = Z_SUMK
            ENDDO
         ENDDO
      ENDDO

      DO JN = 1,5
         DO JT = 1,19
            IPRSM = 0
            DO IGC = 1,NGC(3)
              Z_SUMK = 0.
               DO IPR = 1, NGN(NGS(2)+IGC)
                  IPRSM = IPRSM + 1
                  Z_SUMK = Z_SUMK + KBO_MN2O(JN,JT,IPRSM)*RWGT(IPRSM+32)
               ENDDO
               KB_MN2O(JN,JT,IGC) = Z_SUMK
            ENDDO
         ENDDO
      ENDDO



DO JT = 1,10
  IPRSM = 0
  DO IGC = 1,NGC(3)
    Z_SUMK = 0.0_JPRB
    DO IPR = 1, NGN(NGS(2)+IGC)
      IPRSM = IPRSM + 1
      Z_SUMK = Z_SUMK + SELFREFO(JT,IPRSM)*RWGT(IPRSM+32)
    ENDDO
    SELFREF(JT,IGC) = Z_SUMK
  ENDDO
ENDDO

      DO JT = 1,4
         IPRSM = 0
         DO IGC = 1,NGC(3)
            Z_SUMK = 0.
            DO IPR = 1, NGN(NGS(2)+IGC)
               IPRSM = IPRSM + 1
               Z_SUMK = Z_SUMK + FORREFO(JT,IPRSM)*RWGT(IPRSM+32)
            ENDDO
            FORREF(JT,IGC) = Z_SUMK
         ENDDO
      ENDDO

      DO JP = 1,9
         IPRSM = 0
         DO IGC = 1,NGC(3)
            Z_SUMF = 0.
            DO IPR = 1, NGN(NGS(2)+IGC)
               IPRSM = IPRSM + 1
               Z_SUMF = Z_SUMF + FRACREFAO(IPRSM,JP)
            ENDDO
            FRACREFA(IGC,JP) = Z_SUMF
         ENDDO
      ENDDO



DO JP = 1,5
  IPRSM = 0
  DO IGC = 1,NGC(3)
    Z_SUMF = 0.0_JPRB
    DO IPR = 1, NGN(NGS(2)+IGC)
      IPRSM = IPRSM + 1

      Z_SUMF = Z_SUMF + FRACREFBO(IPRSM,JP)
    ENDDO

    FRACREFB(IGC,JP) = Z_SUMF
  ENDDO
ENDDO


IF (LHOOK) CALL DR_HOOK('RRTM_CMBGB3',1,ZHOOK_HANDLE)
END SUBROUTINE RRTM_CMBGB3
