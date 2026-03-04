!***************************************************************************
SUBROUTINE RRTM_CMBGB5
!***************************************************************************

!     BAND 5:  700-820 cm-1 (low - H2O,CO2; high - O3,CO2)
!     ABozzo 201306 updated to rrtmg v4.85
!     band 5:  700-820 cm-1 (low key - h2o,co2; low minor - o3, ccl4)
!                           (high key - o3,co2)
!***************************************************************************

! Parameters
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOERRTO5 , ONLY : KAO     ,KBO     ,SELFREFO   ,FORREFO, FRACREFAO  ,&
 & FRACREFBO, CCL4O, KAO_MO3
USE YOERRTA5 , ONLY : KA      ,KB      ,SELFREF    ,FORREF, FRACREFA   ,&
 & FRACREFB , CCL4, KA_MO3  
USE YOERRTRWT, ONLY : RWGT
USE YOERRTFTR, ONLY : NGC      ,NGS      ,NGN      

IMPLICIT NONE

INTEGER(KIND=JPIM) :: IGC, IPR, IPRSM, JN, JP, JT

REAL(KIND=JPRB) :: Z_SUMF, Z_SUMK
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('RRTM_CMBGB5',0,ZHOOK_HANDLE)
DO JN = 1,9
  DO JT = 1,5
    DO JP = 1,13
      IPRSM = 0
      DO IGC = 1,NGC(5)
        Z_SUMK = 0.0_JPRB
        DO IPR = 1, NGN(NGS(4)+IGC)
          IPRSM = IPRSM + 1

          Z_SUMK = Z_SUMK + KAO(JN,JT,JP,IPRSM)*RWGT(IPRSM+64)
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
      DO IGC = 1,NGC(5)
        Z_SUMK = 0.0_JPRB
        DO IPR = 1, NGN(NGS(4)+IGC)
          IPRSM = IPRSM + 1

          Z_SUMK = Z_SUMK + KBO(JN,JT,JP,IPRSM)*RWGT(IPRSM+64)
        ENDDO

        KB(JN,JT,JP,IGC) = Z_SUMK
      ENDDO
    ENDDO
  ENDDO
ENDDO

 DO JN = 1,9
    DO JT = 1,19
       IPRSM = 0
       DO IGC = 1,NGC(5)
          Z_SUMK = 0.0_JPRB
          DO IPR = 1, NGN(NGS(4)+IGC)
               IPRSM = IPRSM + 1
               Z_SUMK = Z_SUMK + KAO_MO3(JN,JT,IPRSM)*RWGT(IPRSM+64)
          ENDDO
          KA_MO3(JN,JT,IGC) = Z_SUMK
       ENDDO
    ENDDO
ENDDO



DO JT = 1,10
  IPRSM = 0
  DO IGC = 1,NGC(5)
    Z_SUMK = 0.0_JPRB
    DO IPR = 1, NGN(NGS(4)+IGC)
      IPRSM = IPRSM + 1

      Z_SUMK = Z_SUMK + SELFREFO(JT,IPRSM)*RWGT(IPRSM+64)
    ENDDO

    SELFREF(JT,IGC) = Z_SUMK
  ENDDO
ENDDO

DO JT = 1,4
   IPRSM = 0
   DO IGC = 1,NGC(5)
      Z_SUMK = 0.0_JPRB
      DO IPR = 1, NGN(NGS(4)+IGC)
         IPRSM = IPRSM + 1
         Z_SUMK = Z_SUMK + FORREFO(JT,IPRSM)*RWGT(IPRSM+64)
      ENDDO
      FORREF(JT,IGC) = Z_SUMK
   ENDDO
ENDDO



DO JP = 1,9
  IPRSM = 0
  DO IGC = 1,NGC(5)
    Z_SUMF = 0.0_JPRB
    DO IPR = 1, NGN(NGS(4)+IGC)
      IPRSM = IPRSM + 1

      Z_SUMF = Z_SUMF + FRACREFAO(IPRSM,JP)
    ENDDO

    FRACREFA(IGC,JP) = Z_SUMF
  ENDDO
ENDDO

DO JP = 1,5
  IPRSM = 0
  DO IGC = 1,NGC(5)
    Z_SUMF = 0.0_JPRB
    DO IPR = 1, NGN(NGS(4)+IGC)
      IPRSM = IPRSM + 1

      Z_SUMF = Z_SUMF + FRACREFBO(IPRSM,JP)
    ENDDO

    FRACREFB(IGC,JP) = Z_SUMF
  ENDDO
ENDDO

IPRSM = 0
DO IGC = 1,NGC(5)
  Z_SUMK = 0.0_JPRB
  DO IPR = 1, NGN(NGS(4)+IGC)
    IPRSM = IPRSM + 1

    Z_SUMK = Z_SUMK + CCL4O(IPRSM)*RWGT(IPRSM+64)
  ENDDO

  CCL4(IGC) = Z_SUMK
ENDDO


IF (LHOOK) CALL DR_HOOK('RRTM_CMBGB5',1,ZHOOK_HANDLE)
END SUBROUTINE RRTM_CMBGB5
