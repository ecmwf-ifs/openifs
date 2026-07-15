      SUBROUTINE CTCOR (KIJS, KIJL, F, CTR)

! ----------------------------------------------------------------------

!**** *CTCOR* - COMPUTATION OF CREST-TROUGH CORRELATION

!*    PURPOSE.
!     --------

!       COMPUTE THE CREST-TROUGH CORRELATION (only over the discretised part of the spectrum)

!**   INTERFACE.
!     ----------

!       *CALL* *CTCOR (KIJS, KIJL, F, CTR)*
!              *KIJS*    - INDEX OF FIRST GRIDPOINT
!              *KIJL*    - INDEX OF LAST GRIDPOINT
!              *F*       - SPECTRUM.
!              *CTR*     - CREST-TROUGH CORRELATION.

!     METHOD.
!     -------

!       NONE.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       HAEFNER D, GEMMRICH J, JOCHUM J, 2021: FOWD: a free ocean wave dataset for data mining and machine learning

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : FR       ,DFIM ,  DFIMFR
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : ZMISS    ,PI

      USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      REAL(KIND=JWRB), INTENT(IN) :: F(KIJS:KIJL,NANG,NFRE)
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(OUT) :: CTR 


      INTEGER(KIND=JWIM) :: IJ, K, M
      REAL(KIND=JWRB) :: FR1M1, ZARG, ZAMP
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: EM, ZT1 
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: ZRHO, ZLAM
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL, NFRE) :: TEMP

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('CTCOR',0,ZHOOK_HANDLE)


      DO IJ=KIJS,KIJL
        EM(IJ) = 0.0_JWRB
        ZT1(IJ) = 0.0_JWRB
        ZRHO(IJ) = 0.0_JWRB
        ZLAM(IJ) = 0.0_JWRB
      ENDDO

      DO M=1,NFRE
        DO IJ=KIJS,KIJL
          TEMP(IJ,M) = 0.0_JWRB
        ENDDO
        DO K=1,NANG
          DO IJ=KIJS,KIJL
            TEMP(IJ,M) = TEMP(IJ,M)+F(IJ,K,M)
          ENDDO
        ENDDO
        DO IJ=KIJS,KIJL
          EM(IJ) = EM(IJ)+DFIM(M)*TEMP(IJ,M)
          ZT1(IJ) = ZT1(IJ)+DFIMFR(M)*TEMP(IJ,M)
        ENDDO
      ENDDO


      FR1M1 = 1.0_JWRB/FR(1)

      DO IJ=KIJS,KIJL
        IF(ZT1(IJ) > 0.0_JWRB ) THEN
          ZT1(IJ) = EM(IJ)/ZT1(IJ)
          ZT1(IJ) = MIN(ZT1(IJ),FR1M1)
        ELSE
          ZT1(IJ) = 0.0_JWRB
        ENDIF
      ENDDO


!    COMPUTE THE CREST-TROUGH CORRELATION (only over the discretised part of the spectrum)

      DO M=1,NFRE
        DO IJ=KIJS,KIJL
          ZARG = PI*FR(M)*ZT1(IJ)
          ZAMP = DFIM(M)*TEMP(IJ,M)
          ZRHO(IJ) = ZRHO(IJ) + ZAMP*COS(ZARG)
          ZLAM(IJ) = ZLAM(IJ) + ZAMP*SIN(ZARG)
        ENDDO
      ENDDO

      DO IJ=KIJS,KIJL
        IF(EM(IJ) > 0.0_JWRB ) THEN
          CTR(IJ) = SQRT( ZRHO(IJ)**2 + ZLAM(IJ)**2) / EM(IJ)
        ELSE
          CTR(IJ) = ZMISS
        ENDIF
      ENDDO

      IF (LHOOK) CALL DR_HOOK('CTCOR',1,ZHOOK_HANDLE)

      END SUBROUTINE CTCOR
