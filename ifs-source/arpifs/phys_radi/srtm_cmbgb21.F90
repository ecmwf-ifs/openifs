SUBROUTINE SRTM_CMBGB21

!     BAND 21:  6150-7700 cm-1 (low - H2O,CO2; high - H2O,CO2)
!-----------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM , JPRB
USE YOMHOOK   ,ONLY : LHOOK, DR_HOOK, JPHOOK

USE YOESRTM  , ONLY : NGN
USE YOESRTWN , ONLY : NGC, NGS, RWGT
!USE YOESRTWN , ONLY : NGC, NGS, NGN, RWGT
USE YOESRTA21, ONLY : KA, KB, SELFREF, FORREF, SFLUXREF, &
                    & KAC, KBC, SELFREFC, FORREFC, SFLUXREFC

IMPLICIT NONE

! Local variables
INTEGER(KIND=JPIM) :: JN, JT, JP, IGC, IPR, IPRSM
REAL(KIND=JPRB)    :: ZSUMK, ZSUMF

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SRTM_CMBGB21',0,ZHOOK_HANDLE)

DO JN = 1,9
  DO JT = 1,5
    DO JP = 1,13
      IPRSM = 0
      DO IGC = 1,NGC(6)
        ZSUMK = 0.
        DO IPR = 1, NGN(NGS(5)+IGC)
          IPRSM = IPRSM + 1
          ZSUMK = ZSUMK + KA(JN,JT,JP,IPRSM)*RWGT(IPRSM+80)
        ENDDO
        KAC(JN,JT,JP,IGC) = ZSUMK
      ENDDO
    ENDDO
  ENDDO
ENDDO

DO JN = 1,5
  DO JT = 1,5
    DO JP = 13,59
      IPRSM = 0
      DO IGC = 1,NGC(6)
        ZSUMK = 0.
        DO IPR = 1, NGN(NGS(5)+IGC)
          IPRSM = IPRSM + 1
          ZSUMK = ZSUMK + KB(JN,JT,JP,IPRSM)*RWGT(IPRSM+80)
        ENDDO
        KBC(JN,JT,JP,IGC) = ZSUMK
      ENDDO
    ENDDO
  ENDDO
ENDDO

DO JT = 1,10
  IPRSM = 0
  DO IGC = 1,NGC(6)
    ZSUMK = 0.
    DO IPR = 1, NGN(NGS(5)+IGC)
      IPRSM = IPRSM + 1
      ZSUMK = ZSUMK + SELFREF(JT,IPRSM)*RWGT(IPRSM+80)
    ENDDO
    SELFREFC(JT,IGC) = ZSUMK
  ENDDO
ENDDO

DO JT = 1,4
  IPRSM = 0
  DO IGC = 1,NGC(6)
    ZSUMK = 0.
    DO IPR = 1, NGN(NGS(5)+IGC)
      IPRSM = IPRSM + 1
      ZSUMK = ZSUMK + FORREF(JT,IPRSM)*RWGT(IPRSM+80)
    ENDDO
    FORREFC(JT,IGC) = ZSUMK
  ENDDO
ENDDO

DO JP = 1,9
  IPRSM = 0
  DO IGC = 1,NGC(6)
    ZSUMF = 0.
    DO IPR = 1, NGN(NGS(5)+IGC)
      IPRSM = IPRSM + 1
      ZSUMF = ZSUMF + SFLUXREF(IPRSM,JP)
    ENDDO
    SFLUXREFC(IGC,JP) = ZSUMF
  ENDDO
ENDDO

!     -----------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SRTM_CMBGB21',1,ZHOOK_HANDLE)
END SUBROUTINE SRTM_CMBGB21

