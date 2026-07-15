! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE JONSWAP (ALPHAJ, ZGAMMA, SA, SB, FP, KIJS, KIJL, ET)

! ---------------------------------------------------------------------

!**** *JONSWAP* - ROUTINE TO COMPUTE THE 1-D JONSWAP SPECTRUM.

!     

!*    PURPOSE.
!     --------

!     TO COMPUTE THE 1-D JONSWAP SPECTRUM. 

!**   INTERFACE.
!     ----------

!       *CALL* *JONSWAP (ALPHAJ, ZGAMMA, SA, SB, FP, KIJS, KIJL, ET)*
!        *ALPHAJ*  -  OVERALL ENERGY LEVEL OF JONSWAP SPECTRUM
!        *ZGAMMA*  -  OVERSHOOT FACTOR.
!        *SA*      -  LEFT PEAK WIDTH.
!        *SB*      -  RIGHT PEAK WIDTH.
!        *FP*      -  PEAK FREQUENCY.
!        *KIJS*    -  FIRST POINT IN BLOCK.
!        *KIJL*    -  LAST  POINT IN BLOCK.
!        *ET*      -  JOSWAP SPECTRUM AT DISCRETISED FRQUENCY BINS.


!     METHOD.
!     -------

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : FR
      USE YOWPARAM , ONLY : NFRE
      USE YOWPCONS , ONLY : G        ,ZPI      ,ZPI4GM2

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      REAL(KIND=JWRB), INTENT(IN)  :: ZGAMMA, SA, SB
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: ALPHAJ, FP
      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE), INTENT(OUT) :: ET


      INTEGER(KIND=JWIM) :: M, IJ

      REAL(KIND=JWRB) :: SIGMA, G2ZPI4FRH5M
      REAL(KIND=JWRB) :: FRH, EARG, FJON, FMPF, FJONH
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('JONSWAP',0,ZHOOK_HANDLE)

      DO M=1,NFRE
        FRH = FR(M)
        G2ZPI4FRH5M=1.0_JWRB/(FRH**5*ZPI4GM2)
        DO IJ=KIJS,KIJL
          IF (ALPHAJ(IJ) /= 0.0_JWRB .AND. FP(IJ) /= 0.0_JWRB) THEN
            IF (FRH > FP(IJ)) THEN
              SIGMA = SB
            ELSE
              SIGMA = SA
            ENDIF
            EARG = 0.5_JWRB*((FRH-FP(IJ)) / (SIGMA*FP(IJ)))**2
            EARG = MIN(EARG,50.0_JWRB)
            FJON = ZGAMMA**EXP(-EARG)
            FMPF = 1.25_JWRB*(FP(IJ)/FRH)**4
            FMPF = MIN(FMPF,50.0_JWRB)
            FJONH = EXP(-FMPF)
            ET(IJ,M) = ALPHAJ(IJ)*G2ZPI4FRH5M*FJONH*FJON
          ELSE
            ET(IJ,M) = 0.0_JWRB
          ENDIF
        ENDDO
      ENDDO

      IF (LHOOK) CALL DR_HOOK('JONSWAP',1,ZHOOK_HANDLE)

      END SUBROUTINE JONSWAP
