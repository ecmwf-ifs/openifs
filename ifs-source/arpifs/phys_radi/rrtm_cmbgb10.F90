!***************************************************************************
SUBROUTINE RRTM_CMBGB10
!***************************************************************************

!     BAND 10:  1390-1480 cm-1 (low - H2O; high - H2O)
!     ABozzo updated to rrtmg v4.85
!***************************************************************************

! Parameters
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOERRTO10, ONLY : KAO     ,KBO      ,FRACREFAO   ,FRACREFBO, SELFREFO,FORREFO
USE YOERRTA10, ONLY : KA      ,KB       ,FRACREFA    ,FRACREFB, SELFREF,FORREF
USE YOERRTRWT, ONLY : RWGT
USE YOERRTFTR, ONLY : NGC      ,NGS      ,NGN      

IMPLICIT NONE

INTEGER(KIND=JPIM) :: IGC, IPR, IPRSM, JP, JT

REAL(KIND=JPRB) :: Z_SUMF1, Z_SUMF2, Z_SUMK
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('RRTM_CMBGB10',0,ZHOOK_HANDLE)
DO JT = 1,5
  DO JP = 1,13
    IPRSM = 0
    DO IGC = 1,NGC(10)
      Z_SUMK = 0.0_JPRB
      DO IPR = 1, NGN(NGS(9)+IGC)
        IPRSM = IPRSM + 1

        Z_SUMK = Z_SUMK + KAO(JT,JP,IPRSM)*RWGT(IPRSM+144)
      ENDDO

      KA(JT,JP,IGC) = Z_SUMK
    ENDDO
  ENDDO
ENDDO
DO JT = 1,5
  DO JP = 13,59
    IPRSM = 0
    DO IGC = 1,NGC(10)
      Z_SUMK = 0.0_JPRB
      DO IPR = 1, NGN(NGS(9)+IGC)
        IPRSM = IPRSM + 1

        Z_SUMK = Z_SUMK + KBO(JT,JP,IPRSM)*RWGT(IPRSM+144)
      ENDDO

      KB(JT,JP,IGC) = Z_SUMK
    ENDDO
  ENDDO
ENDDO

  DO JT = 1,10
         IPRSM = 0
         DO IGC = 1,NGC(10)
            Z_SUMK = 0.0_JPRB
            DO IPR = 1, NGN(NGS(9)+IGC)
               IPRSM = IPRSM + 1
               Z_SUMK = Z_SUMK + SELFREFO(JT,IPRSM)*RWGT(IPRSM+144)
            ENDDO
            SELFREF(JT,IGC) = Z_SUMK
         ENDDO
      ENDDO

      DO JT = 1,4
         IPRSM = 0
         DO IGC = 1,NGC(10)
            Z_SUMK = 0.0_JPRB
            DO IPR = 1, NGN(NGS(9)+IGC)
               IPRSM = IPRSM + 1
               Z_SUMK = Z_SUMK + FORREFO(JT,IPRSM)*RWGT(IPRSM+144)
            ENDDO
            FORREF(JT,IGC) = Z_SUMK
         ENDDO
      ENDDO

      IPRSM = 0
      DO IGC = 1,NGC(10)
         Z_SUMF1= 0.0_JPRB
         Z_SUMF2= 0.0_JPRB
         DO IPR = 1, NGN(NGS(9)+IGC)
            IPRSM = IPRSM + 1
            Z_SUMF1= Z_SUMF1+ FRACREFAO(IPRSM)
            Z_SUMF2= Z_SUMF2+ FRACREFBO(IPRSM)
         ENDDO
         FRACREFA(IGC) = Z_SUMF1
         FRACREFB(IGC) = Z_SUMF2
      ENDDO


IF (LHOOK) CALL DR_HOOK('RRTM_CMBGB10',1,ZHOOK_HANDLE)
END SUBROUTINE RRTM_CMBGB10
