! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE WSIGSTAR (KIJS, KIJL, WSWAVE, UFRIC, Z0M, WSTAR, SIG_N)
! ----------------------------------------------------------------------

!**** *WSIGSTAR* - COMPUTATION OF THE RELATIVE STANDARD DEVIATION OF USTAR.

!*    PURPOSE.
!     ---------

!     COMPUTES THE STANDARD DEVIATION OF USTAR DUE TO SMALL SCALE GUSTINESS
!     RELATIVE TO USTAR

!**   INTERFACE.
!     ----------

!     *CALL* *WSIGSTAR (KIJS, KIJL, WSWAVE, UFRIC, Z0M, WSTAR, SIG_N)
!             *KIJS*   - INDEX OF FIRST GRIDPOINT.
!             *KIJL*   - INDEX OF LAST GRIDPOINT.
!             *WSWAVE* - 10M WIND SPEED (m/s).
!             *UFRIC*  - NEW FRICTION VELOCITY IN M/S.
!             *Z0M*    - ROUGHNESS LENGTH IN M.
!             *WSTAR*  - FREE CONVECTION VELOCITY SCALE (M/S).
!             *SIG_N*  - ESTINATED RELATIVE STANDARD DEVIATION OF USTAR.

!     METHOD.
!     -------

!     USE PANOFSKY (1991) TO EXPRESS THE STANDARD DEVIATION OF U10 IN TERMS
!     USTAR AND  w* THE CONVECTIVE VELOCITY SCALE.
!     (but with the background gustiness set to 0.)
!     and USTAR=SQRT(Cd)*U10 to DERIVE THE STANDARD DEVIATION OF USTAR.
!     WITH CD=A+B*U10 (see below).

!     REFERENCE.
!     ----------

!     SEE SECTION 3.2.1 OF THE WAM DOCUMENTATION.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUP  , ONLY : LLGCBZ0
      USE YOWPCONS , ONLY : G, EPSUS, ACDLIN, BCDLIN 
      USE YOWPHYS  , ONLY : XKAPPA, RNUM, ALPHAMIN, ALPHAMAX
      USE YOWWIND  , ONLY : WSPMIN

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: WSWAVE, UFRIC, Z0M, WSTAR
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(OUT) :: SIG_N

      INTEGER(KIND=JWIM) :: IJ

      REAL(KIND=JWRB), PARAMETER :: BG_GUST = 0.0_JWRB ! NO BACKGROUND GUSTINESS (S0 12. IS NOT USED)
      REAL(KIND=JWRB), PARAMETER :: ONETHIRD = 1.0_JWRB/3.0_JWRB
      REAL(KIND=JWRB), PARAMETER :: SIG_NMAX = 0.9_JWRB ! MAX OF RELATIVE STANDARD DEVIATION OF USTAR 

      REAL(KIND=JWRB), PARAMETER :: LOG10 = LOG(10.0_JWRB)
      REAL(KIND=JWRB), PARAMETER :: C1 = 1.03E-3_JWRB
      REAL(KIND=JWRB), PARAMETER :: C2 = 0.04E-3_JWRB
      REAL(KIND=JWRB), PARAMETER :: P1 = 1.48_JWRB
      REAL(KIND=JWRB), PARAMETER :: P2 = -0.21_JWRB

      REAL(KIND=JWRB) :: ZCHAR, C_D, DC_DDU, SIG_CONV
      REAL(KIND=JWRB) :: XKAPPAD, U10, C2U10P1, U10P2
      REAL(KIND=JWRB) :: BCD, U10M1, ZN, Z0VIS
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE


! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('WSIGSTAR',0,ZHOOK_HANDLE)

      IF (LLGCBZ0) THEN
        ZN = RNUM

        DO IJ=KIJS,KIJL
          U10M1=1.0_JWRB/MAX(WSWAVE(IJ),WSPMIN)
          ! CHARNOCK:
          Z0VIS = ZN/MAX(UFRIC(IJ),EPSUS)
          ZCHAR=G*(Z0M(IJ)-Z0VIS)/MAX(UFRIC(IJ)**2,EPSUS)
          ZCHAR=MAX(MIN(ZCHAR,ALPHAMAX),ALPHAMIN)

          BCD = BCDLIN*SQRT(ZCHAR)
          C_D = ACDLIN + BCD * WSWAVE(IJ)
          DC_DDU = BCD
          SIG_CONV = 1.0_JWRB + 0.5_JWRB*WSWAVE(IJ)/C_D * DC_DDU
          SIG_N(IJ) = MIN(SIG_NMAX, SIG_CONV * U10M1*(BG_GUST*UFRIC(IJ)**3 + &
     &                    0.5_JWRB*XKAPPA*WSTAR(IJ)**3)**ONETHIRD )
        ENDDO

       ELSE
        ZN = 0.0_JWRB


!!! for consistency I have kept the old method, even though the new method above could be used,
!!! but until LLGCBZ0 is the default, keep the old scheme whe it is not...
!
!       IN THE FOLLOWING U10 IS ESTIMATED ASSUMING EVERYTHING IS
!       BASED ON U*
!
        XKAPPAD=1.0_JWRB/XKAPPA
        DO IJ=KIJS,KIJL
          U10 = UFRIC(IJ)*XKAPPAD*(LOG10-LOG(Z0M(IJ)))
          U10 = MAX(U10,WSPMIN)
          U10M1=1.0_JWRB/U10
          C2U10P1=C2*U10**P1
          U10P2=U10**P2
          C_D = (C1 + C2U10P1)*U10P2
          DC_DDU = (P2*C1+(P1+P2)*C2U10P1)*U10P2*U10M1
          SIG_CONV = 1.0_JWRB + 0.5_JWRB*U10/C_D*DC_DDU
          SIG_N(IJ) = MIN(SIG_NMAX, SIG_CONV * U10M1*(BG_GUST*UFRIC(IJ)**3 + &
     &                    0.5_JWRB*XKAPPA*WSTAR(IJ)**3)**ONETHIRD )
        ENDDO

      ENDIF

      IF (LHOOK) CALL DR_HOOK('WSIGSTAR',1,ZHOOK_HANDLE)

      END SUBROUTINE WSIGSTAR
