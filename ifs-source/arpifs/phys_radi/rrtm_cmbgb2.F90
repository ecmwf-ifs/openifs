!***************************************************************************
SUBROUTINE RRTM_CMBGB2
!***************************************************************************

!     BAND 2:  250-500 cm-1 (low - H2O; high - H2O)
! ABozzo May 2013 updated to last version of rrtmg
!     band 2:  350-500 cm-1 (low key - h2o; high key - h2o)
!***************************************************************************

! Parameters
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOERRTO2 , ONLY : KAO     ,KBO     ,SELFREFO   ,FRACREFAO  ,&
 & FRACREFBO  ,FORREFO  
USE YOERRTA2 , ONLY : KA      ,KB      ,SELFREF    ,FRACREFA   ,&
 & FRACREFB   ,FORREF       
USE YOERRTRWT, ONLY : RWGT
USE YOERRTFTR, ONLY : NGC      ,NGS      ,NGN      

IMPLICIT NONE

INTEGER(KIND=JPIM) :: IGC, IPR, IPRSM, JP, JT

REAL(KIND=JPRB) :: Z_SUMF, Z_SUMK
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('RRTM_CMBGB2',0,ZHOOK_HANDLE)
DO JT = 1,5
  DO JP = 1,13
    IPRSM = 0
    DO IGC = 1,NGC(2)
      Z_SUMK = 0.0_JPRB
      DO IPR = 1, NGN(NGS(1)+IGC)
        IPRSM = IPRSM + 1

        Z_SUMK = Z_SUMK + KAO(JT,JP,IPRSM)*RWGT(IPRSM+16)
      ENDDO

      KA(JT,JP,IGC) = Z_SUMK
    ENDDO
  ENDDO
  DO JP = 13,59
    IPRSM = 0
    DO IGC = 1,NGC(2)
      Z_SUMK = 0.0_JPRB
      DO IPR = 1, NGN(NGS(1)+IGC)
        IPRSM = IPRSM + 1

        Z_SUMK = Z_SUMK + KBO(JT,JP,IPRSM)*RWGT(IPRSM+16)
      ENDDO
!               KBC(JT,JP,IGC) = SUMK
      KB(JT,JP,IGC) = Z_SUMK
    ENDDO
  ENDDO
ENDDO

DO JT = 1,10
  IPRSM = 0
  DO IGC = 1,NGC(2)
    Z_SUMK = 0.0_JPRB
    DO IPR = 1, NGN(NGS(1)+IGC)
      IPRSM = IPRSM + 1

      Z_SUMK = Z_SUMK + SELFREFO(JT,IPRSM)*RWGT(IPRSM+16)
    ENDDO

    SELFREF(JT,IGC) = Z_SUMK
  ENDDO
ENDDO

DO JT = 1,4
   IPRSM = 0
   DO IGC = 1,NGC(2)
      Z_SUMK = 0.0_JPRB
        DO IPR = 1, NGN(NGS(1)+IGC)
           IPRSM = IPRSM + 1
           Z_SUMK = Z_SUMK + FORREFO(JT,IPRSM)*RWGT(IPRSM+16)
        ENDDO
      FORREF(JT,IGC) = Z_SUMK
   ENDDO
ENDDO


IPRSM = 0
DO IGC = 1,NGC(2)
  Z_SUMK = 0.0_JPRB
  Z_SUMF = 0.0_JPRB
  DO IPR = 1, NGN(NGS(1)+IGC)
    IPRSM = IPRSM + 1

    Z_SUMK = Z_SUMK + FRACREFAO(IPRSM)
    Z_SUMF = Z_SUMF + FRACREFBO(IPRSM)
  ENDDO

  FRACREFA(IGC) = Z_SUMK
  FRACREFB(IGC) = Z_SUMF
ENDDO

!DO JP = 1,13
!  DO IGC = 1,NGC(2)

!    FREFA(NGS(1)+IGC,JP) = FRACREFA(IGC,JP)
!  ENDDO
!ENDDO
!DO JP = 2,13
!  DO IGC = 1,NGC(2)

!    FREFADF(NGS(1)+IGC,JP) = FRACREFA(IGC,JP-1) -FRACREFA(IGC,JP)
!  ENDDO
!ENDDO
!DO IGC = 1,NGC(2)

!  FREFB(NGS(1)+IGC,1) = FRACREFB(IGC)
!ENDDO

IF (LHOOK) CALL DR_HOOK('RRTM_CMBGB2',1,ZHOOK_HANDLE)
END SUBROUTINE RRTM_CMBGB2
