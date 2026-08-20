! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE TAUT_Z0(KIJS, KIJL, IUSFG,          &
&                  HALP, UTOP, UDIR, TAUW, TAUWDIR, RNFAC, &
&                  USTAR, Z0, Z0B, CHRNCK)

! ----------------------------------------------------------------------

!**** *TAUT_Z0* - COMPUTATION OF TOTAL STRESS AND ROUGHNESS LENGTH SCALE.


!**   INTERFACE.
!     ----------

!       *CALL* *TAUT_Z0(KIJS, KIJL, IUSFG, FL1, WAVNUM,
!                       UTOP, UDIR, TAUW, TAUWDIR, RNFAC,
!                       USTAR, Z0, Z0B, CHRNCK)
!          *KIJS*    - INDEX OF FIRST GRIDPOINT
!          *KIJL*    - INDEX OF LAST GRIDPOINT
!          *IUSFG*   - IF = 1 THEN USE THE FRICTION VELOCITY (US) AS FIRST GUESS in TAUT_Z0
!                           0 DO NOT USE THE FIELD USTAR
!          *FL1*     - 2D-SPECTRA
!          *WAVNUM*  - WAVE NUMBER
!          *HALP*    - 1/2 PHILLIPS PARAMETER
!          *UTOP*    - WIND SPEED AT REFERENCE LEVEL XNLEV
!          *UDIR*    - WIND SPEED DIRECTION AT REFERENCE LEVEL XNLEV
!          *TAUW*    - WAVE STRESS.
!          *TAUWDIR* - WAVE STRESS DIRECTION.
!          *RNFAC*   - WIND DEPENDENT FACTOR USED IN THE GROWTH RENORMALISATION.
!          *USTAR*   - FRICTION VELOCITY
!          *Z0*      - ROUGHNESS LENGTH
!          *Z0B*     - BACKGROUND ROUGHNESS LENGTH
!          *CHRNCK*  - CHARNOCK COEFFICIENT

!     METHOD.
!     -------

!       A STEADY STATE WIND PROFILE IS ASSUMED.
!       THE WIND STRESS IS COMPUTED USING THE ROUGHNESS LENGTH

!                  Z1=Z0/SQRT(1-TAUW/TAU)

!       WHERE Z0 IS THE CHARNOCK RELATION , TAUW IS THE WAVE-
!       INDUCED STRESS AND TAU IS THE TOTAL STRESS.
!       WE SEARCH FOR STEADY-STATE SOLUTIONS FOR WHICH TAUW/TAU < 1.

!       IT WAS EXTENDED TO INCLUDE THE GRAVITY-CAPILLARY MODEL FOR THE CALCULATION
!       OF THE BACKGROUND ROUGHNESS.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       FOR QUASILINEAR EFFECT SEE PETER A.E.M. JANSSEN,1990.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUP  , ONLY : LLCAPCHNK, LLGCBZ0
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : G, GM1, EPSUS, EPSMIN, ACD, BCD, CDMAX
      USE YOWPHYS  , ONLY : XKAPPA, XNLEV, RNU, RNUM, ALPHA, ALPHAMIN, ALPHAMAX, &
     &                      ANG_GC_A, ANG_GC_B, ANG_GC_C
      USE YOWTABL  , ONLY : EPS1 

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "chnkmin.intfb.h"
#include "stress_gc.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL, IUSFG
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: HALP, UTOP, UDIR, TAUW, TAUWDIR, RNFAC
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(INOUT) :: USTAR
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(OUT) :: Z0, Z0B, CHRNCK


      INTEGER(KIND=JWIM), PARAMETER :: NITER=17

      REAL(KIND=JWRB), PARAMETER :: TWOXMP1=3.0_JWRB

      INTEGER(KIND=JWIM) :: IJ, ITER
      INTEGER(KIND=JWIM) :: IFRPH

      ! Cd and Z0 from Hersbach 2010, ECMWF Tech Memo (without the viscous part)
!     CD = ACDLIN + BCDLIN*SQRT(PCHAR) * U10
      REAL(KIND=JWRB), PARAMETER :: ACDLIN=0.0008_JWRB
      REAL(KIND=JWRB), PARAMETER :: BCDLIN=0.00047_JWRB
      REAL(KIND=JWRB) :: ALPHAGM1

      REAL(KIND=JWRB), PARAMETER :: Z0MIN = 0.000001_JWRB
      REAL(KIND=JWRB) :: PCE_GC
      REAL(KIND=JWRB) :: Z0MINRST
      REAL(KIND=JWRB) :: CHARNOCK_MIN
      REAL(KIND=JWRB) :: COSDIFF 
      REAL(KIND=JWRB) :: ZCHAR
      REAL(KIND=JWRB) :: US2TOTAUW, USMAX
      REAL(KIND=JWRB) :: XLOGXL, XKUTOP, XOLOGZ0
      REAL(KIND=JWRB) :: USTOLD, USTNEW, TAUOLD, TAUNEW, X, F, DELF, CDFG
      REAL(KIND=JWRB) :: USTM1, Z0TOT, Z0CH, Z0VIS, HZ0VISO1MX, ZZ
      REAL(KIND=JWRB) :: CONST, TAUV, DEL
      REAL(KIND=JWRB) :: RNUEFF, RNUKAPPAM1
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: ALPHAOG, XMIN
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: W1
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: TAUWACT, TAUWEFF 
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: ANG_GC, ALPHAP, TAUUNR

      LOGICAL,  DIMENSION(KIJS:KIJL) :: LLCOSDIFF

! ----------------------------------------------------------------------

!     INLINE FUNCTION.
!     ----------------

!     Simple empirical fit to model drag coefficient
      REAL(KIND=JWRB) :: CDM, U10

      CDM(U10) = MAX(MIN(0.0006_JWRB+0.00008_JWRB*U10, 0.001_JWRB+0.0018_JWRB*EXP(-0.05_JWRB*(U10-33._JWRB))),0.001_JWRB)

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('TAUT_Z0',0,ZHOOK_HANDLE)

      XLOGXL=LOG(XNLEV)
      US2TOTAUW=1.0_JWRB+EPS1

!     ONLY take the contribution of TAUW that is in the wind direction
      DO IJ = KIJS, KIJL
        COSDIFF = COS(UDIR(IJ)-TAUWDIR(IJ))
        TAUWACT(IJ) = MAX(TAUW(IJ)*COSDIFF, EPSMIN )
        LLCOSDIFF(IJ) = (COSDIFF > 0.9_JWRB )
      ENDDO

!  USING THE CG MODEL:
IF (LLGCBZ0) THEN

      IF (LLCAPCHNK) THEN
        DO IJ=KIJS,KIJL
          CHARNOCK_MIN = CHNKMIN(UTOP(IJ))
          ALPHAOG(IJ) = CHARNOCK_MIN*GM1
        ENDDO
      ELSE
        ALPHAOG(KIJS:KIJL)= 0.0_JWRB
      ENDIF

      DO IJ = KIJS, KIJL
        USMAX = MAX(-0.21339_JWRB + 0.093698_JWRB*UTOP(IJ) -0.0020944_JWRB*UTOP(IJ)**2 + 5.5091E-5_JWRB*UTOP(IJ)**3, 0.03_JWRB)
        TAUWEFF(IJ) = MIN(TAUWACT(IJ)*US2TOTAUW, USMAX**2 )
      ENDDO

      RNUEFF = 0.04_JWRB*RNU

      RNUKAPPAM1 = RNUEFF/XKAPPA

      PCE_GC = 0.001_JWRB * IUSFG + (1-IUSFG) * 0.005_JWRB

      IF (IUSFG == 0 ) THEN
        ALPHAGM1 = ALPHA*GM1
        DO IJ = KIJS, KIJL
          IF ( UTOP(IJ) < 1.0_JWRB ) THEN
            CDFG = 0.002_JWRB
          ELSEIF ( LLCOSDIFF(IJ) ) THEN
            X = MIN(TAUWACT(IJ)/MAX(USTAR(IJ),EPSUS)**2,0.99_JWRB)
            ZCHAR = MIN( ALPHAGM1 * USTAR(IJ)**2 / SQRT(1.0_JWRB - X), 0.05_JWRB*EXP(-0.05_JWRB*(UTOP(IJ)-35._JWRB)) )
            ZCHAR = MIN(ZCHAR,ALPHAMAX)
            CDFG = ACDLIN + BCDLIN*SQRT(ZCHAR) * UTOP(IJ)
          ELSE
            CDFG = CDM(UTOP(IJ))
          ENDIF
          USTAR(IJ) = UTOP(IJ)*SQRT(CDFG)
        ENDDO
      ENDIF

      DO IJ = KIJS, KIJL
        W1(IJ) = 0.85_JWRB - 0.05_JWRB*( TANH(10.0_JWRB*(UTOP(IJ)-5.0_JWRB)) + 1.0_JWRB )
      ENDDO

      DO IJ = KIJS, KIJL
        XKUTOP = XKAPPA * UTOP(IJ)

        USTOLD = USTAR(IJ)
        TAUOLD = USTOLD**2

        DO ITER=1,NITER
!         Z0 IS DERIVED FROM THE NEUTRAL LOG PROFILE: UTOP = (USTAR/XKAPPA)*LOG((XNLEV+Z0)/Z0)
          Z0(IJ) = MAX(XNLEV/(EXP(MIN(XKUTOP/USTOLD, 50.0_JWRB))-1.0_JWRB), Z0MIN)
          ! Viscous kinematic stress nu_air * dU/dz at z=0 of the neutral log profile reduced by factor 25 (0.04)
          TAUV = RNUKAPPAM1*USTOLD/Z0(IJ)

          ANG_GC(IJ) = ANG_GC_A + ANG_GC_B * TANH(ANG_GC_C * TAUOLD)

          TAUUNR(IJ) = STRESS_GC(ANG_GC(IJ), USTAR(IJ), Z0(IJ), Z0MIN, HALP(IJ), RNFAC(IJ))

!         TOTAL kinematic STRESS:
          TAUNEW = TAUWEFF(IJ) + TAUV + TAUUNR(IJ)
          USTNEW = SQRT(TAUNEW)
          USTAR(IJ) = W1(IJ)*USTOLD+(1.0_JWRB-W1(IJ))*USTNEW

!         CONVERGENCE ?
          DEL = USTAR(IJ)-USTOLD
          IF (ABS(DEL) < PCE_GC*USTAR(IJ)) EXIT 
          TAUOLD = USTAR(IJ)**2
          USTOLD = USTAR(IJ)
        ENDDO
        ! protection just in case there is no convergence
        IF (ITER > NITER ) THEN
          CDFG = CDM(UTOP(IJ))
          USTAR(IJ) = UTOP(IJ)*SQRT(CDFG)
          Z0MINRST = USTAR(IJ)**2 * ALPHA*GM1
          Z0(IJ) = MAX(XNLEV/(EXP(XKUTOP/USTAR(IJ))-1.0_JWRB), Z0MINRST)
          Z0B(IJ) = Z0MINRST
        ELSE
          Z0(IJ) = MAX(XNLEV/(EXP(XKUTOP/USTAR(IJ))-1.0_JWRB), Z0MIN)
          Z0B(IJ) = Z0(IJ)*SQRT(TAUUNR(IJ)/TAUOLD)
        ENDIF

!       Refine solution
        X = TAUWEFF(IJ)/TAUOLD

        IF (X < 0.99_JWRB) THEN
          USTOLD = USTAR(IJ)
          TAUOLD = MAX(USTOLD**2,TAUWEFF(IJ))

          DO ITER=1,NITER
            X = MIN(TAUWEFF(IJ)/TAUOLD, 0.99_JWRB)
            USTM1 = 1.0_JWRB/MAX(USTOLD,EPSUS)
            !!!! Limit how small z0 could become
            !!!! This is a bit of a compromise to limit very low Charnock for intermediate high winds (15 -25 m/s)
            !!!! It is not ideal !!!
            Z0(IJ) = MAX(XNLEV/(EXP(MIN(XKUTOP/USTOLD, 50.0_JWRB))-1.0_JWRB), Z0MIN)

            TAUUNR(IJ) = STRESS_GC(ANG_GC(IJ), USTOLD, Z0(IJ), Z0MIN, HALP(IJ), RNFAC(IJ))

            Z0B(IJ) = MAX( Z0(IJ)*SQRT(TAUUNR(IJ)/TAUOLD), ALPHAOG(IJ)*TAUOLD)
            Z0VIS = RNUM*USTM1
            HZ0VISO1MX = 0.5_JWRB*Z0VIS/(1.0_JWRB-X)
            Z0(IJ) = HZ0VISO1MX+SQRT(HZ0VISO1MX**2+Z0B(IJ)**2/(1.0_JWRB-X))

            XOLOGZ0= 1.0_JWRB/(XLOGXL-LOG(Z0(IJ)))
            F = USTOLD-XKUTOP*XOLOGZ0
            ZZ = 2.0_JWRB*USTM1*(3.0_JWRB*Z0B(IJ)**2+0.5_JWRB*Z0VIS*Z0(IJ)-Z0(IJ)**2) &
&                / (2.0_JWRB*Z0(IJ)**2*(1.0_JWRB-X)-Z0VIS*Z0(IJ))

            DELF= 1.0_JWRB-XKUTOP*XOLOGZ0**2*ZZ
            IF (DELF /= 0.0_JWRB) USTAR(IJ) = USTOLD-F/DELF

!           CONVERGENCE ?
            DEL = USTAR(IJ)-USTOLD

            IF (ABS(DEL) < PCE_GC*USTAR(IJ)) EXIT 
            USTOLD = USTAR(IJ)
            TAUOLD = MAX(USTOLD**2,TAUWEFF(IJ))
          ENDDO
          ! protection just in case there is no convergence
          IF (ITER > NITER ) THEN
            CDFG = CDM(UTOP(IJ))
            USTAR(IJ) = UTOP(IJ)*SQRT(CDFG)
            Z0MINRST = USTAR(IJ)**2 * ALPHA*GM1
            Z0(IJ) = MAX(XNLEV/(EXP(XKUTOP/USTAR(IJ))-1.0_JWRB), Z0MINRST)
            Z0B(IJ) = Z0MINRST
            CHRNCK(IJ) = MAX(G*Z0(IJ)/USTAR(IJ)**2, ALPHAMIN)
          ELSE
            CHRNCK(IJ) = MAX( G*(Z0B(IJ)/SQRT(1.0_JWRB-X))/MAX(USTAR(IJ),EPSUS)**2, ALPHAMIN)
          ENDIF

        ELSE
          USTM1 = 1.0_JWRB/MAX(USTAR(IJ), EPSUS)
          Z0VIS = RNUM*USTM1
          CHRNCK(IJ) = MAX(G*(Z0(IJ)-Z0VIS) * USTM1**2, ALPHAMIN)
        ENDIF

      ENDDO


ELSE

      DO IJ = KIJS, KIJL
        TAUWEFF(IJ) = TAUWACT(IJ)*US2TOTAUW
      ENDDO

      IF (LLCAPCHNK) THEN
        DO IJ=KIJS,KIJL
          CHARNOCK_MIN = CHNKMIN(UTOP(IJ))
          XMIN(IJ) = 0.15_JWRB*(ALPHA-CHARNOCK_MIN)
          ALPHAOG(IJ) = CHARNOCK_MIN*GM1
        ENDDO
      ELSE
        DO IJ=KIJS,KIJL
          XMIN(IJ)= 0.0_JWRB
          ALPHAOG(IJ)= ALPHA*GM1
        ENDDO
      ENDIF

      DO IJ=KIJS,KIJL
        XKUTOP = XKAPPA * UTOP(IJ)

        USTOLD = (1-IUSFG)*UTOP(IJ)*SQRT(MIN(ACD+BCD*UTOP(IJ),CDMAX)) + IUSFG*USTAR(IJ)
        TAUOLD = MAX(USTOLD**2,TAUWEFF(IJ))
        USTAR(IJ) = SQRT(TAUOLD)
        USTM1 = 1.0_JWRB/MAX(USTAR(IJ),EPSUS) 

        DO ITER=1,NITER
          X = MAX(TAUWACT(IJ)/TAUOLD,XMIN(IJ))
          Z0CH = ALPHAOG(IJ)*TAUOLD/SQRT(1.0_JWRB-X)
          Z0VIS = RNUM*USTM1
          Z0TOT = Z0CH+Z0VIS

          XOLOGZ0= 1.0_JWRB/(XLOGXL-LOG(Z0TOT))
          F = USTAR(IJ)-XKUTOP*XOLOGZ0
          ZZ = USTM1*(Z0CH*(2.0_JWRB-TWOXMP1*X)/(1.0_JWRB-X)-Z0VIS)/Z0TOT
          DELF= 1.0_JWRB-XKUTOP*XOLOGZ0**2*ZZ

          IF (DELF /= 0.0_JWRB) USTAR(IJ) = USTAR(IJ)-F/DELF
          TAUNEW = MAX(USTAR(IJ)**2,TAUWEFF(IJ))
          USTAR(IJ) = SQRT(TAUNEW)
          IF (TAUNEW == TAUOLD) EXIT
          USTM1 = 1.0_JWRB/MAX(USTAR(IJ),EPSUS)
          TAUOLD = TAUNEW
        ENDDO

        Z0(IJ) = Z0CH
        Z0B(IJ) = ALPHAOG(IJ)*TAUOLD
        CHRNCK(IJ) = MAX(G*Z0(IJ)*USTM1**2, ALPHAMIN)

      ENDDO

ENDIF

IF (LHOOK) CALL DR_HOOK('TAUT_Z0',1,ZHOOK_HANDLE)

END SUBROUTINE TAUT_Z0
