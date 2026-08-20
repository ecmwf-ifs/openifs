! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE SEBTMEAN (KIJS, KIJL, FL1, TB, TT, EBT)

! ----------------------------------------------------------------------

!**** *SEBTMEAN* - COMPUTATION OF SPECTRAL VARIANCE FOR ALL WAVES
!                  WITH  TB <= PERIODS <= TT
!                  IF NEEDED A HIGH FREQUENCY f**-5 TAIL IS ASSUMED !!!!

!*    PURPOSE.
!     --------

!       TO COMPUTE TOTAL ENERGY AT EACH GRID POINT.

!**   INTERFACE.
!     ----------

!       *CALL* *SEBTMEAN(KIJS, KIJL, FL1, TB, TT, EBT)*
!          *KIJS*    - INDEX OF FIRST GRIDPOINT
!          *KIJL*    - INDEX OF LAST GRIDPOINT
!          *FL1*     - SPECTRUM.
!          *TB*      - BOTTOM PERIOD
!          *TT*      - TOP PERIOD
!          *EBT*     - MEAN VARIANCE

!     METHOD.
!     -------

!       NONE.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : FR, FR5, DFIM, DELTH, FLOGSPRDM1
      USE YOWPARAM , ONLY : NANG, NFRE
      USE YOWPCONS , ONLY : EPSMIN

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL, NANG, NFRE), INTENT(IN) :: FL1
      REAL(KIND=JWRB), INTENT(IN) :: TB, TT
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(OUT) :: EBT


      INTEGER(KIND=JWIM) :: IJ, M, K, MCUTB, MCUTT

      REAL(KIND=JWRB) :: FCUTB, FCUTT, DFCUT, FBOT, FTOP, ZW
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(NFRE) :: DFIMLOC
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: TEMP

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('SEBTMEAN',0,ZHOOK_HANDLE)

      DFIMLOC(:)=DFIM(:)

      FBOT=1.0_JWRB/MAX(TT,EPSMIN)
      FCUTB=MAX(FR(1),MIN(FBOT,FR(NFRE)))
      FBOT=MAX(FBOT,FR(NFRE))

      FTOP=1.0_JWRB/MAX(TB,EPSMIN)
      FCUTT=MAX(FR(1),MIN(FTOP,FR(NFRE)))
      FTOP=MAX(FTOP,FR(NFRE))

      MCUTB=NINT(LOG10(FCUTB/FR(1))*FLOGSPRDM1)+1
      MCUTB=MIN(MAX(1,MCUTB),NFRE)
      DFCUT=0.5_JWRB*MAX(0.0_JWRB,FCUTB-0.5_JWRB*(FR(MAX(1,MCUTB-1))+FR(MCUTB)))
      DFIMLOC(MCUTB)=MIN(DFCUT*DELTH,DFIM(MCUTB))
      DFIMLOC(MCUTB)=DFIM(MCUTB)-DFIMLOC(MCUTB)

      MCUTT=NINT(LOG10(FCUTT/FR(1))*FLOGSPRDM1)+1
      MCUTT=MIN(MAX(1,MCUTT),NFRE)
      DFCUT=0.5_JWRB*MAX(0.0_JWRB,FCUTT-0.5_JWRB*(FR(MAX(1,MCUTT-1))+FR(MCUTT)))
      DFIMLOC(MCUTT)=MIN(DFCUT*DELTH,DFIM(MCUTT))

      IF(FCUTB == FCUTT) MCUTT=MCUTB-1

!*    1. INITIALISE ENERGY ARRAY.
!        ------------------------

      DO IJ=KIJS,KIJL
        EBT(IJ) = EPSMIN 
      ENDDO

! ----------------------------------------------------------------------

!*    2. INTEGRATE OVER FREQUENCIES AND DIRECTION.
!        -----------------------------------------

      DO M=MCUTB,MCUTT
        K=1
        DO IJ=KIJS,KIJL
          TEMP(IJ) = FL1(IJ,K,M)
        ENDDO
        DO K=2,NANG
          DO IJ=KIJS,KIJL
            TEMP(IJ) = TEMP(IJ)+FL1(IJ,K,M)
          ENDDO
        ENDDO
        DO IJ=KIJS,KIJL
          EBT(IJ) = EBT(IJ)+DFIMLOC(M)*TEMP(IJ)
        ENDDO
      ENDDO

!     CHECK IF A THE REQUESTED FREQUENCIES ARE ABOVE FR(NFRE)
      IF ( FBOT < FTOP ) THEN

        ZW = 0.25_JWRB * DELTH * FR5(NFRE) * ( 1.0_JWRB/FBOT**4 - 1.0_JWRB/FTOP**4 )
        K=1
        DO IJ=KIJS,KIJL
          TEMP(IJ) = FL1(IJ,K,NFRE)
        ENDDO
        DO K=2,NANG
          DO IJ=KIJS,KIJL
            TEMP(IJ) = TEMP(IJ)+FL1(IJ,K,NFRE)
          ENDDO
        ENDDO
        DO IJ=KIJS,KIJL
          EBT(IJ) = EBT(IJ) + ZW*TEMP(IJ)
        ENDDO

      ENDIF 

      IF (LHOOK) CALL DR_HOOK('SEBTMEAN',1,ZHOOK_HANDLE)

      END SUBROUTINE SEBTMEAN
