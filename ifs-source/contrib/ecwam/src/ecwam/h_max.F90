! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE H_MAX(C3,C4,XNSLC,KIJS,KIJL,AA,BB,HMAXN,SIG_HM)
 
!***  DETERMINE EXPECTED MAXIMUM WAVE HEIGHT, NORMALISED WITH
!     SIGNIFICANT WAVE HEIGHT.
 
!     PURPOSE:
!     -------
 
!     DETERMINE EXPECTED MAXIMUM ENVELOPE WAVE HEIGHT WHICH IS OBTAINED
!     FROM EXPECTED MAXIMUM WAVE ENERGY
!                 /
!         <E_MAX>=| DE E p_MAX(E)            
!                 /
!     WHERE P_MAX(E) IS DERIVED IN JANSSEN (2015) AND BASED ON WORK DONE BY
!     MONTINA ET AL (2009). MAXIMUM WAVE HEIGHT THEN FOLLOWS FROM
!     <H_MAX>=SQRT(<E_MAX>/2. 

!     VARIABLE       TYPE         PURPOSE
!     --------       ----         -------
 
!     C_3            REAL         SKEWNESS
!     C_4            REAL         KURTOSIS
!     XNSLC          REAL         NUMBER OF SIGNIFICANT LEVEL CROSSINGS
!     KIJS            INTEGER      FIRST INDEX
!     KIJL            INTEGER      LAST INDEX
!     AA             REAL         FIRST PARAMETER OF RESIDORI PDF
!     BB             REAL         =2(AA+1)
 
!     OUTPUT:
!     ------

!     HMAXN          REAL         MAXIMUM WAVE HEIGHT NORMALIZED WITH
!                                 SIGNIFICANT WAVE HEIGHT H_S. 
!     SIG_HM         REAL         WIDTH OF HMAX DISTRIBUTION 

!     AUTHOR:
!     ------
!     P.A.E.M. JANSSEN, NOVEMBER 2017
 
!----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWPCONS , ONLY : PI

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

!----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: C3, C4, XNSLC
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(OUT) :: AA, BB, HMAXN, SIG_HM


      INTEGER(KIND=JWIM), PARAMETER :: NITER = 5
      INTEGER(KIND=JWIM) :: IJ, I

!     H_C   FIRST GUESS OF EXPECTED NORMALISED MAX WAVE HEIGHT
      REAL(KIND=JWRB), PARAMETER :: H_C= 2._JWRB 
      REAL(KIND=JWRB), PARAMETER :: H_C_MIN= 1._JWRB 
      REAL(KIND=JWRB), PARAMETER :: H_C_MAX= 4._JWRB 

      REAL(KIND=JWRB), PARAMETER :: TWOSQRT6=2._JWRB*SQRT(6._JWRB)
      REAL(KIND=JWRB), PARAMETER :: GAM = 0.5772_JWRB
      REAL(KIND=JWRB), PARAMETER :: EB = 10._JWRB
      REAL(KIND=JWRB), PARAMETER :: FLOGMIN = 0.1_JWRB
      REAL(KIND=JWRB), PARAMETER :: AA_MAX = 1000._JWRB

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB) :: ZEPSILON
      REAL(KIND=JWRB) :: TWOG1, G2, AE, BE, F, Z0, EMIN, EMAX, EVAL
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL):: E, BBM1, DFNORMA

!----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('H_MAX',0,ZHOOK_HANDLE)

!*    1. DETERMINE EXPECTED MAXIMUM WAVE HEIGHT.
!     -----------------------------------------

      ZEPSILON=10._JWRB*EPSILON(ZEPSILON)
      TWOG1 = -2._JWRB*GAM
      G2 = GAM**2+PI**2/6._JWRB
      AE = 0.5_JWRB*EB*(EB-2._JWRB)
      BE = 0.5_JWRB*EB*(EB**2-6._JWRB*EB+6._JWRB)

      EMIN = 2._JWRB*H_C_MIN**2
      EMAX = 2._JWRB*H_C_MAX**2
      EVAL = 2._JWRB*H_C**2
      DO IJ=KIJS,KIJL
        E(IJ) = EVAL
      ENDDO

      DO IJ=KIJS,KIJL
        DFNORMA(IJ)=C4(IJ)*AE+C3(IJ)**2*BE
        IF (XNSLC(IJ) > 0._JWRB .AND. ABS(DFNORMA(IJ)).GT.ZEPSILON) THEN
          F  = LOG(MAX(1._JWRB+DFNORMA(IJ),FLOGMIN))

          AA(IJ) = MIN(((EB-F)**2-2._JWRB*EB)/(2._JWRB*F),AA_MAX)
          BB(IJ) = 2._JWRB*(1._JWRB+AA(IJ))
          BBM1(IJ) = 1._JWRB/(BB(IJ)+ZEPSILON*SIGN(1.0_JWRB,BB(IJ)))

          DO I=1,NITER
            Z0 = LOG(XNSLC(IJ)*SQRT(0.5_JWRB*E(IJ)))
            E(IJ) = (G2-TWOG1*(AA(IJ)+Z0)+(2._JWRB*AA(IJ)+Z0)*Z0)*BBM1(IJ)
            E(IJ)=MIN(MAX(E(IJ),EMIN),EMAX)
          ENDDO
          HMAXN(IJ) = SQRT(0.5_JWRB*E(IJ))
          SIG_HM(IJ) = PI*HMAXN(IJ)/(TWOSQRT6*(Z0+0.5_JWRB*GAM))
        ELSE
          AA(IJ) = 0._JWRB
          BB(IJ) = 2._JWRB*(1._JWRB+AA(IJ))
          HMAXN(IJ) = H_C_MIN
          SIG_HM(IJ)= 0._JWRB
        ENDIF
      ENDDO

      IF (LHOOK) CALL DR_HOOK('H_MAX',1,ZHOOK_HANDLE)

      END SUBROUTINE H_MAX
