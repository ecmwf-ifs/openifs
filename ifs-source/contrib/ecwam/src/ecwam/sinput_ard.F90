! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE SINPUT_ARD (NGST, LLSNEG, KIJS, KIJL, FL1, &
 &                     WAVNUM, CINV, XK2CG,           &
 &                     WDWAVE, WSWAVE, UFRIC, Z0M,    &
 &                     COSWDIF, SINWDIF2,             &
 &                     RAORW, WSTAR, RNFAC,           &
 &                     FLD, SL, SPOS, XLLWS)
! ----------------------------------------------------------------------

!**** *SINPUT_ARD* - COMPUTATION OF INPUT SOURCE FUNCTION.


!*    PURPOSE.
!     ---------

!       COMPUTE THE WIND INPUT SOURCE TRERM BASED ON ARDHUIN ET AL. 2010.

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

!     *CALL* *SINPUT_ARD (NGST, LLSNEG, KIJS, KIJL, FL1,
!    &                    WAVNUM, CINV, XK2CG,
!    &                    WSWAVE, WDWAVE, UFRIC, Z0M,
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
!        *UFRIC* - NEW FRICTION VELOCITY IN M/S.
!        *Z0M* - ROUGHNESS LENGTH IN M.
!      *COSWDIF* - COS(TH(K)-WDWAVE(IJ))
!     *SINWDIF2* - SIN(TH(K)-WDWAVE(IJ))**2
!        *RAORW* - RATIO AIR DENSITY TO WATER DENSITY.
!        *WSTAR* - FREE CONVECTION VELOCITY SCALE (M/S).
!        *RNFAC* - WIND DEPENDENT FACTOR USED IN THE GROWTH RENORMALISATION.
!          *FLD* - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE.
!           *SL* - TOTAL SOURCE FUNCTION ARRAY.
!         *SPOS* - POSITIVE SOURCE FUNCTION ARRAY.
!       *XLLWS*  - = 1 WHERE SINPUT IS POSITIVE

!     METHOD.
!     -------

!       SEE REFERENCE.


! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUP  , ONLY : LLCAPCHNK,LLNORMAGAM
      USE YOWFRED  , ONLY : FR       ,TH       ,DFIM     ,COSTH  ,SINTH, ZPIFR, DELTH
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : G        ,GM1      ,EPSMIN, EPSUS, ZPI
      USE YOWPHYS  , ONLY : ZALP     ,TAUWSHELTER, XKAPPA, BETAMAXOXKAPPA2,    &
     &                      RNU      ,RNUM, &
     &                      SWELLF   ,SWELLF2  ,SWELLF3  ,SWELLF4  , SWELLF5, &
     &                      SWELLF6  ,SWELLF7  ,SWELLF7M1, Z0RAT   ,Z0TUBMAX , &
     &                      ABMIN  ,ABMAX
      USE YOWTEST  , ONLY : IU06
      USE YOWTABL  , ONLY : IAB      ,SWELLFT

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "wsigstar.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: NGST
      LOGICAL, INTENT(IN) :: LLSNEG
      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(IN) :: FL1
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE), INTENT(IN) :: WAVNUM, CINV, XK2CG
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: WDWAVE, WSWAVE, UFRIC, Z0M
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: RAORW, WSTAR, RNFAC
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG), INTENT(IN) :: COSWDIF, SINWDIF2
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(OUT) :: FLD, SL, SPOS
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(OUT) :: XLLWS


      INTEGER(KIND=JWIM) :: IJ, K, M, IND, IGST

      REAL(KIND=JWRB) :: CONSTN 
      REAL(KIND=JWRB) :: AVG_GST, ABS_TAUWSHELTER 
      REAL(KIND=JWRB) :: CONST1
      REAL(KIND=JWRB) :: ZNZ
      REAL(KIND=JWRB) :: X1, X2, ZLOG, ZLOG1, ZLOG2, ZLOG2X, XV1, XV2, ZBETA1, ZBETA2
      REAL(KIND=JWRB) :: XI, X, DELI1, DELI2
      REAL(KIND=JWRB) :: FU, FUD, NU_AIR,SMOOTH, HFTSWELLF6, Z0TUB
      REAL(KIND=JWRB) :: FAC_NU_AIR, FACM1_NU_AIR
      REAL(KIND=JWRB) :: ARG, DELABM1
      REAL(KIND=JWRB) :: TAUPX, TAUPY
      REAL(KIND=JWRB) :: DSTAB2
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(NFRE) :: CONST, SIG, SIGM1, SIG2, COEF, COEF5, DFIM_SIG2
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: CONSTF, CONST11, CONST22
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: Z0VIS, Z0NOZ, FWW
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: PVISC, PTURB
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: ZCN
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: SIG_N, UORBT, AORB, TEMP, RE, RE_C, ZORB
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: CNSN, SUMF, SUMFSIN2
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: CSTRNFAC
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: FLP_AVG, SLP_AVG
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: ROGOROAIR, AIRD_PVISC
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NGST) :: XSTRESS, YSTRESS, FLP, SLP
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NGST) :: USG2, TAUX, TAUY, USTP, USTPM1, USDIRP, UCN
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NGST) :: UCNZALPD
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE) :: XNGAMCONST
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NGST) :: GAMNORMA ! ! RENORMALISATION FACTOR OF THE GROWTH RATE
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: DSTAB1, TEMP1, TEMP2
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NGST) :: COSLP, GAM0, DSTAB

      LOGICAL :: LTAUWSHELTER
! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SINPUT_ARD',0,ZHOOK_HANDLE)

      AVG_GST = 1.0_JWRB/NGST
      CONST1  = BETAMAXOXKAPPA2
      CONSTN = DELTH/(XKAPPA*ZPI) 

      ABS_TAUWSHELTER=ABS(TAUWSHELTER)
      IF (ABS_TAUWSHELTER == 0.0_JWRB ) THEN
        LTAUWSHELTER = .FALSE.
      ELSE
        LTAUWSHELTER = .TRUE.
      ENDIF

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

! ----------------------------------------------------------------------

      DO M=1,NFRE
        SIG(M) = ZPIFR(M)
        SIGM1(M) = 1.0_JWRB/SIG(M)
        SIG2(M) = SIG(M)**2
        DFIM_SIG2(M)=DFIM(M)*SIG2(M)
        CONST(M)=SIG(M)*CONST1
      ENDDO


      IF (LLSNEG) THEN
!!!!  only for the negative sinput
         NU_AIR = RNU
         FACM1_NU_AIR = 4.0_JWRB/NU_AIR

         FAC_NU_AIR = RNUM

         FU=ABS(SWELLF3)
         FUD=SWELLF2
         DELABM1= REAL(IAB)/(ABMAX-ABMIN)


!       computation of Uorb and Aorb
        DO IJ=KIJS,KIJL
          UORBT(IJ) = EPSMIN 
          AORB(IJ) = EPSMIN
        ENDDO

        DO M=1,NFRE
          COEF(M) =-SWELLF*16._JWRB*SIG2(M)/G
          COEF5(M)=-SWELLF5*2._JWRB*SQRT(2._JWRB*NU_AIR*SIG(M))

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
            UORBT(IJ) = UORBT(IJ)+DFIM_SIG2(M)*TEMP(IJ)
            AORB(IJ) = AORB(IJ)+DFIM(M)*TEMP(IJ)
          ENDDO
        ENDDO

        DO IJ=KIJS,KIJL
          UORBT(IJ) = 2.0_JWRB*SQRT(UORBT(IJ))  ! this is the significant orbital amplitude
          AORB(IJ)  = 2.0_JWRB*SQRT(AORB(IJ))   ! this 1/2 Hs
          RE(IJ)    = FACM1_NU_AIR*UORBT(IJ)*AORB(IJ) ! this is the Reynolds number 
          Z0VIS(IJ) = FAC_NU_AIR/MAX(UFRIC(IJ),0.0001_JWRB)
          Z0TUB = Z0RAT*MIN(Z0TUBMAX,Z0M(IJ))
          Z0NOZ(IJ) = MAX(Z0VIS(IJ),Z0TUB)
          ZORB(IJ)  = AORB(IJ)/Z0NOZ(IJ)

!         compute fww
          XI=(LOG10(MAX(ZORB(IJ),3.0_JWRB))-ABMIN)*DELABM1
          IND  = MIN (IAB-1, INT(XI))
          DELI1= MIN (1.0_JWRB ,XI-REAL(IND,JWRB))
          DELI2= 1.0_JWRB - DELI1
          FWW(IJ)= SWELLFT(IND)*DELI2+SWELLFT(IND+1)*DELI1
          TEMP2(IJ) = FWW(IJ)*UORBT(IJ)
        ENDDO

!       Define the critical Reynolds number
        IF ( SWELLF6 == 1.0_JWRB) THEN
          DO IJ=KIJS,KIJL
            RE_C(IJ) = SWELLF4
          ENDDO
        ELSE
          HFTSWELLF6=1.0_JWRB-SWELLF6
          DO IJ=KIJS,KIJL
            RE_C(IJ) = SWELLF4*(2.0_JWRB/AORB(IJ))**HFTSWELLF6
          ENDDO
        ENDIF

!       Swell damping weight between viscous and turbulent boundary layer
        IF (SWELLF7 > 0.0_JWRB) THEN
          DO IJ=KIJS,KIJL
            SMOOTH=0.5_JWRB*TANH((RE(IJ)-RE_C(IJ))*SWELLF7M1)
            PTURB(IJ)=0.5_JWRB+SMOOTH
            PVISC(IJ)=0.5_JWRB-SMOOTH
          ENDDO
        ELSE
          DO IJ=KIJS,KIJL
            IF (RE(IJ) <= RE_C(IJ)) THEN
              PTURB(IJ)=0.0_JWRB
              PVISC(IJ)=0.5_JWRB
            ELSE
              PTURB(IJ)=0.5_JWRB
              PVISC(IJ)=0.0_JWRB
            ENDIF
          ENDDO
        ENDIF

        DO IJ=KIJS,KIJL
          AIRD_PVISC(IJ) = PVISC(IJ)*RAORW(IJ)
        ENDDO

      ENDIF



! Initialisation

      IF (NGST == 1) THEN
        DO IJ=KIJS,KIJL
          USTP(IJ,1) = UFRIC(IJ)
        ENDDO
      ELSE IF (NGST == 2) THEN
        DO IJ=KIJS,KIJL
          USTP(IJ,1) = UFRIC(IJ)*(1.0_JWRB+SIG_N(IJ))
          USTP(IJ,2) = UFRIC(IJ)*(1.0_JWRB-SIG_N(IJ))
        ENDDO
      ELSE
         WRITE (IU06,*) '**************************************'
         WRITE (IU06,*) '*    FATAL ERROR                     *'
         WRITE (IU06,*) '*    ===========                     *'
         WRITE (IU06,*) '* IN SINPUT_ARD: NGST > 2            *'
         WRITE (IU06,*) '* NGST = ', NGST
         WRITE (IU06,*) '* PROGRAM ABORTS.   PROGRAM ABORTS.  *'
         WRITE (IU06,*) '*                                    *'
         WRITE (IU06,*) '**************************************'
         CALL ABORT1
      ENDIF

      DO IGST=1,NGST
        DO IJ=KIJS,KIJL
          USTPM1(IJ,IGST) = 1.0_JWRB/MAX(USTP(IJ,IGST),EPSUS)
        ENDDO
      ENDDO

      IF (LTAUWSHELTER) THEN
        DO IGST=1,NGST
          DO IJ=KIJS,KIJL
            XSTRESS(IJ,IGST)=0.0_JWRB
            YSTRESS(IJ,IGST)=0.0_JWRB
            USG2(IJ,IGST)=USTP(IJ,IGST)**2
            TAUX(IJ,IGST)=USG2(IJ,IGST)*SIN(WDWAVE(IJ))
            TAUY(IJ,IGST)=USG2(IJ,IGST)*COS(WDWAVE(IJ))
          ENDDO
        ENDDO

        DO IJ=KIJS,KIJL
          ROGOROAIR(IJ) = G/RAORW(IJ)
        ENDDO

      ELSE
        DO IGST=1,NGST
          DO K=1,NANG
            DO IJ=KIJS,KIJL
              COSLP(IJ,K,IGST) = COSWDIF(IJ,K)
            ENDDO
          ENDDO
        ENDDO
      ENDIF


!*    2. MAIN LOOP OVER FREQUENCIES.
!        ---------------------------

      IF ( .NOT. LLNORMAGAM) THEN
        GAMNORMA(KIJS:KIJL,:) = 1.0_JWRB
      ENDIF

      IF ( .NOT. LLSNEG) THEN
        DSTAB(KIJS:KIJL,:,:) = 0.0_JWRB
      ENDIF

      DO M=1,NFRE

        IF (LTAUWSHELTER) THEN
          DO IGST=1,NGST
            DO IJ=KIJS,KIJL
              TAUPX=TAUX(IJ,IGST)-ABS_TAUWSHELTER*XSTRESS(IJ,IGST)
              TAUPY=TAUY(IJ,IGST)-ABS_TAUWSHELTER*YSTRESS(IJ,IGST)
              USDIRP(IJ,IGST)=ATAN2(TAUPX,TAUPY)
              USTP(IJ,IGST)=(TAUPX**2+TAUPY**2)**0.25_JWRB
              USTPM1(IJ,IGST)=1.0_JWRB/MAX(USTP(IJ,IGST),EPSUS)
            ENDDO
          ENDDO

          DO IGST=1,NGST
            DO K=1,NANG
              DO IJ=KIJS,KIJL
                COSLP(IJ,K,IGST) = COS(TH(K)-USDIRP(IJ,IGST))
              ENDDO
            ENDDO
          ENDDO

          DO IJ=KIJS,KIJL
            CONSTF(IJ) = ROGOROAIR(IJ)*CINV(IJ,M)*DFIM(M)
          ENDDO
        ENDIF


!*      PRECALCULATE FREQUENCY DEPENDENCE.
!       ----------------------------------

        DO IGST=1,NGST
          DO IJ=KIJS,KIJL
            UCN(IJ,IGST) = USTP(IJ,IGST)*CINV(IJ,M)
            UCNZALPD(IJ,IGST) = XKAPPA/(UCN(IJ,IGST) + ZALP)
          ENDDO
        ENDDO
        DO IJ=KIJS,KIJL
          ZCN(IJ) = LOG(WAVNUM(IJ,M)*Z0M(IJ))
          CNSN(IJ) = CONST(M)*RAORW(IJ)
        ENDDO

!*    2.1 LOOP OVER DIRECTIONS.
!         ---------------------

        DO K=1,NANG
          DO IJ=KIJS,KIJL
            XLLWS(IJ,K,M)=0.0_JWRB
          ENDDO
        ENDDO

        DO IGST=1,NGST
          DO K=1,NANG
            DO IJ=KIJS,KIJL
              IF (COSLP(IJ,K,IGST) > 0.01_JWRB) THEN
                X    = COSLP(IJ,K,IGST)*UCN(IJ,IGST)
                ZLOG = ZCN(IJ) + UCNZALPD(IJ,IGST)/COSLP(IJ,K,IGST)
                IF (ZLOG < 0.0_JWRB) THEN
                  ZLOG2X=ZLOG*ZLOG*X
                  GAM0(IJ,K,IGST) = EXP(ZLOG)*ZLOG2X*ZLOG2X * CNSN(IJ)
                  XLLWS(IJ,K,M) = 1.0_JWRB
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


        IF (LLSNEG) THEN
!       SWELL DAMPING:
          DO IJ=KIJS,KIJL
            DSTAB1(IJ) = COEF5(M)*AIRD_PVISC(IJ)*WAVNUM(IJ,M)
            TEMP1(IJ) = COEF(M)*RAORW(IJ)
          ENDDO

          DO IGST=1,NGST
            DO K=1,NANG
              DO IJ=KIJS,KIJL
                DSTAB2 = TEMP1(IJ)*(TEMP2(IJ)+(FU+FUD*COSLP(IJ,K,IGST))*USTP(IJ,IGST))
                DSTAB(IJ,K,IGST) = DSTAB1(IJ)+PTURB(IJ)*DSTAB2
              ENDDO
            ENDDO
          ENDDO
        ENDIF


!*    2.2 UPDATE THE SHELTERING STRESS (in any),
!         AND THEN ADDING INPUT SOURCE TERM TO NET SOURCE FUNCTION.
!         ---------------------------------------------------------

        DO K=1,NANG

          DO IGST=1,NGST
            DO IJ=KIJS,KIJL
              ! SLP: only the positive contributions
              SLP(IJ,IGST) =  GAM0(IJ,K,IGST) * GAMNORMA(IJ,IGST)
              FLP(IJ,IGST) = SLP(IJ,IGST)+DSTAB(IJ,K,IGST)
            ENDDO
          ENDDO

          DO IGST=1,NGST
            DO IJ=KIJS,KIJL
              SLP(IJ,IGST) = SLP(IJ,IGST)*FL1(IJ,K,M)
            ENDDO
          ENDDO

          IF (LTAUWSHELTER) THEN
            DO IJ=KIJS,KIJL
              CONST11(IJ)=CONSTF(IJ)*SINTH(K)
              CONST22(IJ)=CONSTF(IJ)*COSTH(K)
            ENDDO
            DO IGST=1,NGST
              DO IJ=KIJS,KIJL
                XSTRESS(IJ,IGST)=XSTRESS(IJ,IGST)+SLP(IJ,IGST)*CONST11(IJ)
                YSTRESS(IJ,IGST)=YSTRESS(IJ,IGST)+SLP(IJ,IGST)*CONST22(IJ)
              ENDDO
            ENDDO
          ENDIF

          IGST=1
            DO IJ=KIJS,KIJL
              SLP_AVG(IJ) = SLP(IJ,IGST)
              FLP_AVG(IJ) = FLP(IJ,IGST)
            ENDDO
          DO IGST=2,NGST
            DO IJ=KIJS,KIJL
              SLP_AVG(IJ) = SLP_AVG(IJ)+SLP(IJ,IGST)
              FLP_AVG(IJ) = FLP_AVG(IJ)+FLP(IJ,IGST)
            ENDDO
          ENDDO

          DO IJ=KIJS,KIJL
            SPOS(IJ,K,M) = AVG_GST*SLP_AVG(IJ)
          ENDDO

          DO IJ=KIJS,KIJL
            FLD(IJ,K,M) = AVG_GST*FLP_AVG(IJ)
            SL(IJ,K,M) = FLD(IJ,K,M)*FL1(IJ,K,M)
          ENDDO

        ENDDO

      ENDDO ! END LOOP OVER FREQUENCIES

IF (LHOOK) CALL DR_HOOK('SINPUT_ARD',1,ZHOOK_HANDLE)

END SUBROUTINE SINPUT_ARD
