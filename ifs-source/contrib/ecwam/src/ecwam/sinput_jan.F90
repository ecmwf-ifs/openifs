! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE SINPUT_JAN (NGST, LLSNEG, KIJS, KIJL, FL1 , &
     &                       WAVNUM, CINV, XK2CG,            &
     &                       WSWAVE, UFRIC, Z0M,     &
     &                       COSWDIF, SINWDIF2,              &
     &                       RAORW, WSTAR, RNFAC,            &
                             FLD, SL, SPOS, XLLWS)
! ----------------------------------------------------------------------

!**** *SINPUT_JAN* - COMPUTATION OF INPUT SOURCE FUNCTION.

!     P.A.E.M. JANSSEN    KNMI      AUGUST    1990

!     OPTIMIZED BY : H. GUENTHER

!     MODIFIED BY : 
!       J-R BIDLOT NOVEMBER 1995
!       J-R BIDLOT FEBRUARY 1996-97
!       J-R BIDLOT FEBRUARY 1999 : INTRODUCE ICALL AND NCALL
!       P.A.E.M. JANSSEN MAY 2000 : INTRODUCE GUSTINESS
!       J-R BIDLOT FEBRUARY 2001 : MAKE IT FULLY IMPLICIT BY ONLY
!                                  USING NEW STRESS AND ROUGHNESS. 
!       S. ABDALLA OCTOBER 2001:  INTRODUCTION OF VARIABLE AIR
!                                 DENSITY AND STABILITY-DEPENDENT 
!                                 WIND GUSTINESS
!       P.A.E.M. JANSSEN OCTOBER 2008: INTRODUCE DAMPING WHEN WAVES ARE 
!                                      RUNNING FASTER THAN THE WIND.
!       J-R BIDLOT JANUARY 2013: SHALLOW WATER FORMULATION.

!*    PURPOSE.
!     ---------

!       COMPUTE INPUT SOURCE FUNCTION AND STORE ADDITIVELY INTO NET
!       SOURCE FUNCTION ARRAY, ALSO COMPUTE FUNCTIONAL DERIVATIVE OF
!       INPUT SOURCE FUNCTION.
!
!       GUSTINESS IS INTRODUCED FOLL0WING THE APPROACH OF JANSSEN(1986),
!       USING A GAUSS-HERMITE APPROXIMATION SUGGESTED BY MILES(1997).
!       IN THE PRESENT VERSION ONLY TWO HERMITE POLYNOMIALS ARE UTILISED
!       IN THE EVALUATION OF THE PROBABILITY INTEGRAL. EXPLICITELY ONE THEN
!       FINDS:
!
!             <GAMMA(X)> = 0.5*( GAMMA(X(1+SIG)) + GAMMA(X(1-SIG)) )
!
!       WHERE X IS THE FRICTION VELOCITY AND SIG IS THE RELATIVE GUSTINESS
!       LEVEL.

!**   INTERFACE.
!     ----------

!     *CALL* *SINPUT_JAN (NGST, LLSNEG, KIJS, KIJL, FL1,
!    &                    WAVNUM, CINV, XK2CG,
!    &                    WDWAVE, WSWAVE, UFRIC, Z0M,
!    &                    COSWDIF, SINWDIF2,
!    &                    RAORW, WSTAR, RNFAC,
!    &                    FLD, SL, SPOS, XLLWS)
!         *NGST* - IF = 1 THEN NO GUSTINESS PARAMETERISATION
!                - IF = 2 THEN GUSTINESS PARAMETERISATION
!         *LLSNEG- IF TRUE THEN THE NEGATIVE SINPUT (SWELL DAMPING) WILL BE COMPUTED
!         *KIJS* - INDEX OF FIRST GRIDPOINT.
!         *KIJL* - INDEX OF LAST GRIDPOINT.
!          *FL1* - SPECTRUM.
!       *WAVNUM* - WAVE NUMBER.
!         *CINV* - INVERSE PHASE VELOCITY.
!       *XK2CG*  - (WAVNUM)**2 * GROUP SPPED.
!       *WDWAVE* - WIND DIRECTION IN RADIANS IN OCEANOGRAPHIC
!                  NOTATION (POINTING ANGLE OF WIND VECTOR,
!                  CLOCKWISE FROM NORTH).
!        *UFRIC* - FRICTION VELOCITY IN M/S.
!        *Z0M*   - ROUGHNESS LENGTH IN M.
!      *COSWDIF* - COS(TH(K)-WDWAVE(IJ))
!     *SINWDIF2* - SIN(TH(K)-WDWAVE(IJ))**2
!        *RAORW* - RATIO AIR DENSITY TO WATER DENSITY
!        *RNFAC* - WIND DEPENDENT FACTOR USED IN THE GROWTH RENORMALISATION.
!        *WSTAR* - FREE CONVECTION VELOCITY SCALE (M/S).
!          *FLD* - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE.
!           *SL* - TOTAL SOURCE FUNCTION ARRAY.
!         *SPOS* - ONLY POSITIVE PART OF INPUT SOURCE FUNCTION ARRAY.
!        *XLLWS* - 1 WHERE SINPUT IS POSITIVE.


!     METHOD.
!     -------

!       SEE REFERENCE.

!     EXTERNALS.
!     ----------

!       WSIGSTAR. 

!     MODIFICATIONS
!     -------------

!     - REMOVAL OF CALL TO CRAY SPECIFIC FUNCTIONS EXPHF AND ALOGHF
!       BY THEIR STANDARD FORTRAN EQUIVALENT EXP and ALOGHF
!     - MODIFIED TO MAKE INTEGRATION SCHEME FULLY IMPLICIT
!     - INTRODUCTION OF VARIABLE AIR DENSITY
!     - INTRODUCTION OF WIND GUSTINESS

!     REFERENCE.
!     ----------

!       P. JANSSEN, J.P.O., 1989.
!       P. JANSSEN, J.P.O., 1991

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUP  , ONLY : LLNORMAGAM
      USE YOWFRED  , ONLY : ZPIFR    ,DELTH    ,TH
      USE YOWFRED  , ONLY : FR       ,TH       , ZPIFR
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : G        ,GM1      ,ZPI  , EPSUS
      USE YOWPHYS  , ONLY : ZALP     ,XKAPPA, BETAMAXOXKAPPA2
      USE YOWSTAT  , ONLY : IDAMPING
      USE YOWTEST  , ONLY : IU06

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "abort1.intfb.h"
#include "wsigstar.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: NGST
      LOGICAL, INTENT(IN) :: LLSNEG
      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(IN) :: FL1
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE), INTENT(IN) :: WAVNUM, CINV, XK2CG
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: WSWAVE, UFRIC, Z0M
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG), INTENT(IN) :: COSWDIF, SINWDIF2
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: RAORW, WSTAR, RNFAC
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(OUT) :: FLD, SL, SPOS
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(OUT) :: XLLWS


      INTEGER(KIND=JWIM) :: IJ, IG, K, M
      INTEGER(KIND=JWIM) :: IGST

      REAL(KIND=JWRB) :: CONST1, CONST3, XKAPPAD
      REAL(KIND=JWRB) :: CONSTN
      REAL(KIND=JWRB) :: ZNZ
      REAL(KIND=JWRB) :: X, ZLOG, ZLOG2X, ZBETA
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(NGST) :: WSIN
      REAL(KIND=JWRB), DIMENSION(NFRE) :: CONST
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: ZTANHKD 
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: SIG_N
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: CNSN
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: SUMF, SUMFSIN2 
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: CSTRNFAC
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE) :: XNGAMCONST
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NGST) :: GAMNORMA ! ! RENORMALISATION FACTOR OF THE GROWTH RATE
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NGST) :: SIGDEV ,US, Z0, UCN, ZCN
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NGST) :: USTPM1
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NGST) :: XVD, UCND, CONST3_UCN2
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG) :: UFAC1, UFAC2
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG) :: TEMPD
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NGST) :: GAM0

      LOGICAL, DIMENSION(KIJS:KIJL,NANG) :: LZ

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('SINPUT_JAN',0,ZHOOK_HANDLE)

      CONST1   = BETAMAXOXKAPPA2 
      CONST3   = 2.0_JWRB*XKAPPA/CONST1  ! SEE IDAMPING
      XKAPPAD  = 1.E0_JWRB/XKAPPA

      CONST3 = IDAMPING*CONST3

      CONSTN = DELTH/(XKAPPA*ZPI)

!*    1. PRECALCULATED ANGULAR DEPENDENCE.
!        ---------------------------------

      DO K=1,NANG
        DO IJ=KIJS,KIJL
          IF (COSWDIF(IJ,K) > 0.01_JWRB) THEN
            LZ(IJ,K) = .TRUE.
            TEMPD(IJ,K) = XKAPPA/COSWDIF(IJ,K)
          ELSE
            LZ(IJ,K) = .FALSE.
            TEMPD(IJ,K) = XKAPPA 
          ENDIF
        ENDDO
      ENDDO

      IF (LLNORMAGAM) THEN
        DO IJ=KIJS,KIJL
          CSTRNFAC(IJ) = CONSTN * RNFAC(IJ) / RAORW(IJ)
        ENDDO
        DO M=1,NFRE
          DO IJ=KIJS,KIJL
            XNGAMCONST(IJ,M) = CSTRNFAC(IJ)*XK2CG(IJ,M)
          ENDDO
        ENDDO

      ELSE
        XNGAMCONST(KIJS:KIJL,:) = 0.0_JWRB
      ENDIF


!     ESTIMATE THE STANDARD DEVIATION OF GUSTINESS.
      IF (NGST > 1) CALL WSIGSTAR (KIJS, KIJL, WSWAVE, UFRIC, Z0M, WSTAR, SIG_N)


!     DEFINE WHERE SINPUT WILL BE EVALUATED IN RELATIVE TERM WRT USTAR
!     DEFINE ALSO THE RELATIVE WEIGHT OF EACH.

      IF (NGST == 1) THEN
        WSIN(1) = 1.0_JWRB 
        DO IJ=KIJS,KIJL
          SIGDEV(IJ,1) = 1.0_JWRB
        ENDDO
      ELSE IF (NGST == 2) THEN
        WSIN(1) = 0.5_JWRB 
        WSIN(2) = 0.5_JWRB 
        DO IJ=KIJS,KIJL
          SIGDEV(IJ,1) = 1.0_JWRB-SIG_N(IJ)
          SIGDEV(IJ,2) = 1.0_JWRB+SIG_N(IJ)
        ENDDO
      ELSE
         WRITE (IU06,*) '**************************************'
         WRITE (IU06,*) '*    FATAL ERROR                     *'
         WRITE (IU06,*) '*    ===========                     *'
         WRITE (IU06,*) '* IN SINPUT_JAN: NGST > 2            *'
         WRITE (IU06,*) '* NGST = ', NGST
         WRITE (IU06,*) '* PROGRAM ABORTS.   PROGRAM ABORTS.  *'
         WRITE (IU06,*) '*                                    *'
         WRITE (IU06,*) '**************************************'
         CALL ABORT1
      ENDIF


      IF (NGST == 1) THEN
        DO IJ=KIJS,KIJL
          US(IJ,1) = UFRIC(IJ)
          Z0(IJ,1) = Z0M(IJ)
        ENDDO
      ELSE
        DO IGST=1,NGST
          DO IJ=KIJS,KIJL
            US(IJ,IGST) = UFRIC(IJ)*SIGDEV(IJ,IGST)
            Z0(IJ,IGST) = Z0M(IJ)
          ENDDO
        ENDDO
      ENDIF

      DO IGST=1,NGST
        DO IJ=KIJS,KIJL
          USTPM1(IJ,IGST) = 1.0_JWRB/MAX(US(IJ,IGST),EPSUS)
        ENDDO
      ENDDO

! ----------------------------------------------------------------------

!*    2. LOOP OVER FREQUENCIES.
!        ----------------------

      IF ( .NOT. LLNORMAGAM) THEN
        GAMNORMA(KIJS:KIJL,:) = 1.0_JWRB
      ENDIF

      IF ( .NOT. LLSNEG) THEN
        UFAC2(KIJS:KIJL,:) = 0.0_JWRB
      ENDIF

      DO M=1,NFRE

        CONST(M)=ZPIFR(M)*CONST1

!*      PRECALCULATE FREQUENCY DEPENDENCE.
!       ----------------------------------

        DO IJ=KIJS,KIJL
          ZTANHKD(IJ) = ZPIFR(M)**2/(G*WAVNUM(IJ,M)) 
        ENDDO

        DO IJ=KIJS,KIJL
          CNSN(IJ) = CONST(M)*ZTANHKD(IJ)*RAORW(IJ)
        ENDDO

        DO IGST=1,NGST
          DO IJ=KIJS,KIJL
            UCN(IJ,IGST) = US(IJ,IGST)*CINV(IJ,M) + ZALP
            CONST3_UCN2(IJ,IGST) = CONST3*UCN(IJ,IGST)**2
            UCND(IJ,IGST) = 1.0_JWRB/ UCN(IJ,IGST)
            ZCN(IJ,IGST)  = LOG(WAVNUM(IJ,M)*Z0(IJ,IGST))
            XVD(IJ,IGST) = 1.0_JWRB/(-US(IJ,IGST)*XKAPPAD*ZCN(IJ,IGST)*CINV(IJ,M))
          ENDDO
        ENDDO

!*    2.1 LOOP OVER DIRECTIONS.
!         ---------------------


!       WIND INPUT:
        DO K=1,NANG

          DO IJ=KIJS,KIJL
            XLLWS(IJ,K,M)= 0.0_JWRB
          ENDDO

          DO IGST=1,NGST
            DO IJ=KIJS,KIJL
              IF (LZ(IJ,K)) THEN
                ZLOG = ZCN(IJ,IGST) + TEMPD(IJ,K)*UCND(IJ,IGST)
                IF (ZLOG < 0.0_JWRB) THEN
                  X=COSWDIF(IJ,K)*UCN(IJ,IGST)
                  ZLOG2X=ZLOG*ZLOG*X
                  GAM0(IJ,K,IGST) = ZLOG2X*ZLOG2X*EXP(ZLOG) * CNSN(IJ)
                  XLLWS(IJ,K,M)= 1.0_JWRB
                ELSE
                  GAM0(IJ,K,IGST) = 0.0_JWRB
                ENDIF
              ELSE
                GAM0(IJ,K,IGST) = 0.0_JWRB
              ENDIF
            ENDDO
          ENDDO

        ENDDO


        IF (LLNORMAGAM) THEN

          DO IGST=1,NGST

            SUMF(KIJS:KIJL) = 0.0_JWRB
            SUMFSIN2(KIJS:KIJL) = 0.0_JWRB
            DO K=1,NANG
              DO IJ=KIJS,KIJL
                SUMF(IJ) = SUMF(IJ) + GAM0(IJ,K,IGST)*FL1(IJ,K,M)
                SUMFSIN2(IJ) = SUMFSIN2(IJ) + GAM0(IJ,K,IGST)*FL1(IJ,K,M)*SINWDIF2(IJ,K)
              ENDDO
            ENDDO

            DO IJ=KIJS,KIJL
              ZNZ = XNGAMCONST(IJ,M)*USTPM1(IJ,IGST)
              GAMNORMA(IJ,IGST) = (1.0_JWRB + ZNZ*SUMFSIN2(IJ)) / (1.0_JWRB + ZNZ*SUMF(IJ))
            ENDDO

          ENDDO

        ENDIF


        DO K=1,NANG
          DO IJ=KIJS,KIJL
            UFAC1(IJ,K) = WSIN(1)*GAM0(IJ,K,1)*GAMNORMA(IJ,1)
          ENDDO
          DO IGST=2,NGST
            DO IJ=KIJS,KIJL
              UFAC1(IJ,K) = UFAC1(IJ,K) + WSIN(IGST)*GAM0(IJ,K,IGST)*GAMNORMA(IJ,IGST)
            ENDDO
          ENDDO
        ENDDO

        IF (LLSNEG) THEN
!         SWELL DAMPING:
          DO K=1,NANG
            DO IGST=1,1
              DO IJ=KIJS,KIJL
                ZBETA = CONST3_UCN2(IJ,IGST)*(COSWDIF(IJ,K)-XVD(IJ,IGST))
                UFAC2(IJ,K) = WSIN(IGST)*ZBETA
              ENDDO
            ENDDO
            DO IGST=2,NGST
              DO IJ=KIJS,KIJL
                ZBETA = CONST3_UCN2(IJ,IGST)*(COSWDIF(IJ,K)-XVD(IJ,IGST))
                UFAC2(IJ,K) = UFAC2(IJ,K)+WSIN(IGST)*ZBETA
              ENDDO
            ENDDO
          ENDDO
        ENDIF

!*    2.2 ADDING INPUT SOURCE TERM TO NET SOURCE FUNCTION.
!         ------------------------------------------------

        DO K=1,NANG
          DO IJ=KIJS,KIJL
            FLD(IJ,K,M) = UFAC1(IJ,K) + UFAC2(IJ,K)*CNSN(IJ)
            SPOS(IJ,K,M) = UFAC1(IJ,K)*FL1(IJ,K,M)
            SL(IJ,K,M) = FLD(IJ,K,M)*FL1(IJ,K,M)
          ENDDO
        ENDDO

      ENDDO

      IF (LHOOK) CALL DR_HOOK('SINPUT_JAN',1,ZHOOK_HANDLE)

      END SUBROUTINE SINPUT_JAN
