SUBROUTINE SRTM_CMBGB16


!  Original version:       Michael J. Iacono; July, 1998
!  Revision for RRTM_SW:   Michael J. Iacono; November, 2002
!  Revision for RRTMG_SW:  Michael J. Iacono; December, 2003

!  The subroutines CMBGB16->CMBGB29 input the absorption coefficient
!  data for each band, which are defined for 16 g-points and 14 spectral
!  bands. The data are combined with appropriate weighting following the
!  g-point mapping arrays specified in RRTMG_SW_INIT.  Solar source 
!  function data in array SFLUXREF are combined without weighting.  All
!  g-point reduced data are put into new arrays for use in RRTMG_SW.

!  BAND 16:  2600-3250 cm-1 (low key- H2O,CH4; high key - CH4)

!-----------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM , JPRB
USE YOMHOOK   ,ONLY : LHOOK, DR_HOOK, JPHOOK

USE YOESRTM  , ONLY : NGN
USE YOESRTWN , ONLY : NGC, RWGT
!USE YOESRTWN , ONLY : NGC, NGN, RWGT
USE YOESRTA16, ONLY : KA, KB, SELFREF, FORREF, SFLUXREF, &
                    & KAC, KBC, SELFREFC, FORREFC, SFLUXREFC

IMPLICIT NONE

! Local variables
INTEGER(KIND=JPIM) :: JN, JT, JP, IGC, IPR, IPRSM
REAL(KIND=JPRB)    :: ZSUMK, ZSUMF

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SRTM_CMBGB16',0,ZHOOK_HANDLE)

DO JN = 1,9
  DO JT = 1,5
    DO JP = 1,13
      IPRSM = 0
      DO IGC = 1,NGC(1)
        ZSUMK = 0.
        DO IPR = 1, NGN(IGC)
          IPRSM = IPRSM + 1
          ZSUMK = ZSUMK + KA(JN,JT,JP,IPRSM)*RWGT(IPRSM)
        ENDDO
        KAC(JN,JT,JP,IGC) = ZSUMK
      ENDDO
    ENDDO
  ENDDO
ENDDO

DO JT = 1,5
  DO JP = 13,59
    IPRSM = 0
    DO IGC = 1,NGC(1)
      ZSUMK = 0.
      DO IPR = 1, NGN(IGC)
        IPRSM = IPRSM + 1
        ZSUMK = ZSUMK + KB(JT,JP,IPRSM)*RWGT(IPRSM)
      ENDDO
      KBC(JT,JP,IGC) = ZSUMK
    ENDDO
  ENDDO
ENDDO

DO JT = 1,10
  IPRSM = 0
  DO IGC = 1,NGC(1)
    ZSUMK = 0.
    DO IPR = 1, NGN(IGC)
      IPRSM = IPRSM + 1
      ZSUMK = ZSUMK + SELFREF(JT,IPRSM)*RWGT(IPRSM)
    ENDDO
    SELFREFC(JT,IGC) = ZSUMK
  ENDDO
ENDDO

DO JT = 1,3
  IPRSM = 0
  DO IGC = 1,NGC(1)
    ZSUMK = 0.
    DO IPR = 1, NGN(IGC)
      IPRSM = IPRSM + 1
      ZSUMK = ZSUMK + FORREF(JT,IPRSM)*RWGT(IPRSM)
    ENDDO
    FORREFC(JT,IGC) = ZSUMK
  ENDDO
ENDDO

IPRSM = 0
DO IGC = 1,NGC(1)
  ZSUMF = 0.
  DO IPR = 1, NGN(IGC)
    IPRSM = IPRSM + 1
    ZSUMF = ZSUMF + SFLUXREF(IPRSM)
  ENDDO
  SFLUXREFC(IGC) = ZSUMF
ENDDO

!     -----------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SRTM_CMBGB16',1,ZHOOK_HANDLE)
END SUBROUTINE SRTM_CMBGB16

