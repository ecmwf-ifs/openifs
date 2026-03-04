!***************************************************************************
SUBROUTINE RRTM_CMBGB11
!***************************************************************************

!     BAND 11:  1480-1800 cm-1 (low - H2O; high - H2O)
!     ABozzo updated to rrtmg v4.85
!     band 11:  1480-1800 cm-1 (low - h2o; low minor - o2)
!                              (high key - h2o; high minor - o2)
!***************************************************************************

! Parameters
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOERRTO11, ONLY : KAO     ,KBO     ,SELFREFO,FORREFO    ,FRACREFAO ,FRACREFBO, &
                    & KAO_MO2, KBO_MO2
USE YOERRTA11, ONLY : KA      ,KB      ,SELFREF,FORREF     ,FRACREFA  ,FRACREFB, &
                    & KA_MO2, KB_MO2
USE YOERRTRWT, ONLY : RWGT
USE YOERRTFTR, ONLY : NGC      ,NGS      ,NGN      

IMPLICIT NONE

INTEGER(KIND=JPIM) :: IGC, IPR, IPRSM, JP, JT

REAL(KIND=JPRB) :: Z_SUMF1, Z_SUMF2, Z_SUMK,Z_SUMK1,Z_SUMK2
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('RRTM_CMBGB11',0,ZHOOK_HANDLE)
DO JT = 1,5
  DO JP = 1,13
    IPRSM = 0
    DO IGC = 1,NGC(11)
      Z_SUMK = 0.0_JPRB
      DO IPR = 1, NGN(NGS(10)+IGC)
        IPRSM = IPRSM + 1

        Z_SUMK = Z_SUMK + KAO(JT,JP,IPRSM)*RWGT(IPRSM+160)
      ENDDO

      KA(JT,JP,IGC) = Z_SUMK
    ENDDO
  ENDDO
ENDDO
DO JT = 1,5
  DO JP = 13,59
    IPRSM = 0
    DO IGC = 1,NGC(11)
      Z_SUMK = 0.0_JPRB
      DO IPR = 1, NGN(NGS(10)+IGC)
        IPRSM = IPRSM + 1

        Z_SUMK = Z_SUMK + KBO(JT,JP,IPRSM)*RWGT(IPRSM+160)
      ENDDO

      KB(JT,JP,IGC) = Z_SUMK
    ENDDO
  ENDDO
ENDDO

      DO JT = 1,19
         IPRSM = 0
         DO IGC = 1,NGC(11)
            Z_SUMK1 = 0.0_JPRB
            Z_SUMK2 = 0.0_JPRB
            DO IPR = 1, NGN(NGS(10)+IGC)
               IPRSM = IPRSM + 1
               Z_SUMK1 = Z_SUMK1 + KAO_MO2(JT,IPRSM)*RWGT(IPRSM+160)
               Z_SUMK2 = Z_SUMK2 + KBO_MO2(JT,IPRSM)*RWGT(IPRSM+160)
            ENDDO
            KA_MO2(JT,IGC) = Z_SUMK1
            KB_MO2(JT,IGC) = Z_SUMK2
         ENDDO
      ENDDO


DO JT = 1,10
  IPRSM = 0
  DO IGC = 1,NGC(11)
    Z_SUMK = 0.0_JPRB
    DO IPR = 1, NGN(NGS(10)+IGC)
      IPRSM = IPRSM + 1

      Z_SUMK = Z_SUMK + SELFREFO(JT,IPRSM)*RWGT(IPRSM+160)
    ENDDO

    SELFREF(JT,IGC) = Z_SUMK
  ENDDO
ENDDO

DO JT = 1,4
         IPRSM = 0
         DO IGC = 1,NGC(11)
            Z_SUMK = 0.0_JPRB
            DO IPR = 1, NGN(NGS(10)+IGC)
               IPRSM = IPRSM + 1
               Z_SUMK = Z_SUMK + FORREFO(JT,IPRSM)*RWGT(IPRSM+160)
            ENDDO
            FORREF(JT,IGC) = Z_SUMK
         ENDDO
ENDDO

IPRSM = 0
DO IGC = 1,NGC(11)
         Z_SUMF1= 0.0_JPRB
         Z_SUMF2= 0.0_JPRB
   DO IPR = 1, NGN(NGS(10)+IGC)
            IPRSM = IPRSM + 1
            Z_SUMF1= Z_SUMF1+ FRACREFAO(IPRSM)
            Z_SUMF2= Z_SUMF2+ FRACREFBO(IPRSM)
   ENDDO
         FRACREFA(IGC) = Z_SUMF1
         FRACREFB(IGC) = Z_SUMF2
 ENDDO

IF (LHOOK) CALL DR_HOOK('RRTM_CMBGB11',1,ZHOOK_HANDLE)
END SUBROUTINE RRTM_CMBGB11
