!***************************************************************************
SUBROUTINE RRTM_CMBGB1
!***************************************************************************

!  The subroutines CMBGB1->CMBGB16 input the absorption coefficient
!  data for each band, which are defined for 16 g-points and 16 spectral
!  bands. The data are combined with appropriate weighting following the
!  g-point mapping arrays specified in RRTMINIT.  Plank fraction data
!  in arrays FRACREFA and FRACREFB are combined without weighting.  All
!  g-point reduced data are put into new arrays for use in RRTM.

!  BAND 1:  10-250 cm-1 (low - H2O; high - H2O)

! ABozzo may 2013 update to the last version of rrtmg

!band 1:  10-350 cm-1 (low key - h2o; low minor - n2)
!                       (high key - h2o; high minor - n2)
!***************************************************************************

! Parameters
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOERRTO1 , ONLY : KAO, KBO, SELFREFO, FORREFO, FRACREFAO,FRACREFBO,KAO_MN2, KBO_MN2
USE YOERRTA1 , ONLY : KA , KB , SELFREF , FORREF , FRACREFA ,FRACREFB,  KA_MN2, KB_MN2
USE YOERRTRWT, ONLY : RWGT
USE YOERRTFTR, ONLY : NGC      ,NGN      

IMPLICIT NONE

INTEGER(KIND=JPIM) :: IGC, IPR, IPRSM, JP, JT

REAL(KIND=JPRB) :: Z_SUMF1, Z_SUMF2, Z_SUMK, Z_SUMK1, Z_SUMK2
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('RRTM_CMBGB1',0,ZHOOK_HANDLE)
DO JT = 1,5
  DO JP = 1,13
    IPRSM = 0
    DO IGC = 1,NGC(1)
      Z_SUMK = 0.0_JPRB
      DO IPR = 1, NGN(IGC)
        IPRSM = IPRSM + 1

        Z_SUMK = Z_SUMK + KAO(JT,JP,IPRSM)*RWGT(IPRSM)
      ENDDO

      KA(JT,JP,IGC) = Z_SUMK
    ENDDO
  ENDDO
  DO JP = 13,59
    IPRSM = 0
    DO IGC = 1,NGC(1)
      Z_SUMK = 0.0_JPRB
      DO IPR = 1, NGN(IGC)
        IPRSM = IPRSM + 1

        Z_SUMK = Z_SUMK + KBO(JT,JP,IPRSM)*RWGT(IPRSM)
      ENDDO

      KB(JT,JP,IGC) = Z_SUMK
    ENDDO
  ENDDO
ENDDO

DO JT = 1,10
  IPRSM = 0
  DO IGC = 1,NGC(1)
    Z_SUMK = 0.0_JPRB
    DO IPR = 1, NGN(IGC)
      IPRSM = IPRSM + 1

      Z_SUMK = Z_SUMK + SELFREFO(JT,IPRSM)*RWGT(IPRSM)
    ENDDO

    SELFREF(JT,IGC) = Z_SUMK
  ENDDO
ENDDO

DO JT = 1,4
         IPRSM = 0
         DO IGC = 1,NGC(1)
            Z_SUMK = 0.0_JPRB
            DO IPR = 1, NGN(IGC)
               IPRSM = IPRSM + 1
               Z_SUMK = Z_SUMK + FORREFO(JT,IPRSM)*RWGT(IPRSM)
            ENDDO
            FORREF(JT,IGC) = Z_SUMK
         ENDDO
      ENDDO

      DO JT = 1,19
         IPRSM = 0
         DO IGC = 1,NGC(1)
            Z_SUMK1 = 0.0_JPRB
            Z_SUMK2 = 0.0_JPRB
            DO IPR = 1, NGN(IGC)
               IPRSM = IPRSM + 1
               Z_SUMK1 = Z_SUMK1 + KAO_MN2(JT,IPRSM)*RWGT(IPRSM)
               Z_SUMK2 = Z_SUMK2 + KBO_MN2(JT,IPRSM)*RWGT(IPRSM)
            ENDDO
            KA_MN2(JT,IGC) = Z_SUMK1
            KB_MN2(JT,IGC) = Z_SUMK2
         ENDDO
      ENDDO

IPRSM = 0
DO IGC = 1,NGC(1)
  Z_SUMF1 = 0.0_JPRB
  Z_SUMF2 = 0.0_JPRB
  DO IPR = 1, NGN(IGC)
    IPRSM = IPRSM + 1

    Z_SUMF1= Z_SUMF1+ FRACREFAO(IPRSM)
    Z_SUMF2= Z_SUMF2+ FRACREFBO(IPRSM)
  ENDDO

  FRACREFA(IGC) = Z_SUMF1
  FRACREFB(IGC) = Z_SUMF2
ENDDO


IF (LHOOK) CALL DR_HOOK('RRTM_CMBGB1',1,ZHOOK_HANDLE)
END SUBROUTINE RRTM_CMBGB1
