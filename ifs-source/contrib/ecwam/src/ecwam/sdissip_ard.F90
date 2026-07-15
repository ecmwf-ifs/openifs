! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE SDISSIP_ARD (KIJS, KIJL, FL1, FLD, SL,          &
     &                        WAVNUM, CGROUP, XK2CG,             &
     &                        UFRIC, COSWDIF, RAORW)
! ----------------------------------------------------------------------

!**** *SDISSIP_ARD* - COMPUTATION OF DISSIPATION SOURCE FUNCTION.

!     LOTFI AOUF       METEO FRANCE 2013
!     FABRICE ARDHUIN  IFREMER  2013


!*    PURPOSE.
!     --------
!       COMPUTE DISSIPATION SOURCE FUNCTION AND STORE ADDITIVELY INTO
!       NET SOURCE FUNCTION ARRAY. ALSO COMPUTE FUNCTIONAL DERIVATIVE
!       OF DISSIPATION SOURCE FUNCTION.

!**   INTERFACE.
!     ----------

!       *CALL* *SDISSIP_ARD (KIJS, KIJL, FL1, FLD,SL,*
!                            WAVNUM, CGROUP, XK2CG,
!                            UFRIC, COSWDIF, RAORW)*
!          *KIJS*   - INDEX OF FIRST GRIDPOINT
!          *KIJL*   - INDEX OF LAST GRIDPOINT
!          *FL1*    - SPECTRUM.
!          *FLD*    - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
!          *SL*     - TOTAL SOURCE FUNCTION ARRAY
!          *WAVNUM* - WAVE NUMBER
!          *CGROUP* - GROUP SPEED
!          *XK2CG*  - (WAVE NUMBER)**2 * GROUP SPEED
!          *UFRIC*  - FRICTION VELOCITY IN M/S.
!          *RAORW*  - RATIO AIR DENSITY TO WATER DENSITY
!          *COSWDIF*-  COS(TH(K)-WDWAVE(IJ))


!     METHOD.
!     -------

!       SEE REFERENCES.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       ARDHUIN et AL. JPO DOI:10.1175/20110JPO4324.1


! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : FR      , TH     ,ZPIFR   ,FRATIO    ,DELTH
      USE YOWPCONS , ONLY : G        ,ZPI
      USE YOWPARAM , ONLY : NANG    ,NFRE
      USE YOWPHYS  , ONLY : SDSBR   ,ISDSDTH ,ISB     ,IPSAT    ,      &
&                  SSDSC2  , SSDSC4, SSDSC6,  MICHE, SSDSC3, SSDSBRF1, &
&                  BRKPBCOEF ,SSDSC5, NSDSNTH, NDIKCUMUL,              &
&                  INDICESSAT, SATWEIGHTS

      USE YOMHOOK  , ONLY : LHOOK   ,DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL

      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(IN) :: FL1
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(INOUT) :: FLD, SL
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE), INTENT(IN) :: WAVNUM, CGROUP, XK2CG 
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: UFRIC, RAORW 
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL, NANG), INTENT(IN) :: COSWDIF 


      INTEGER(KIND=JWIM) :: IJ, K, M, I, J, M2, K2, KK, NANGD
      INTEGER(KIND=JWIM), DIMENSION(NANG) :: KKD

      REAL(KIND=JWRB) :: TPIINV, TPIINVH, TMP01, TMP02, TMP03
      REAL(KIND=JWRB) :: EPSR, SSDSC6M1, ZCOEF, ZCOEFM1
      REAL(KIND=JWRB) :: XLOGDFRTH, BRLAMBDA

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE


      REAL(KIND=JWRB), DIMENSION(NFRE) :: SIG, SSDSC2_SIG 
      REAL(KIND=JWRB), DIMENSION(0:NANG/2) :: COSDTH
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: FACTURB
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE) :: FACSAT, FACWTRB, TEMP1 
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE) :: BTH0 !saturation spectrum
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE) :: C_, C_C, DSIP, TRPZ_DSIP
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE) :: BTH !saturation spectrum 
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE) :: TEMP2
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE) :: D
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE) :: SCUMUL 
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE) :: RENEWALFREQ
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,0:NANG/2,NFRE) :: WCUMUL

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('SDISSIP_ARD',0,ZHOOK_HANDLE)

      ! INITIALISATION

      EPSR=SQRT(SDSBR)

      TPIINV = 1.0_JWRB/ZPI
      TPIINVH= 0.5_JWRB*TPIINV


      TMP03 = 1.0_JWRB/(SDSBR*MICHE)

      SSDSC6M1=1._JWRB-SSDSC6

      DO M=1,NFRE
        SIG(M) = ZPIFR(M)
        SSDSC2_SIG(M)=SSDSC2*SIG(M)
      END DO

      DO M=1, NFRE
        DO IJ=KIJS,KIJL
          FACSAT(IJ,M) = WAVNUM(IJ,M)*TPIINV*XK2CG(IJ,M)
        ENDDO
      ENDDO

      ! COMPUTE SATURATION SPECTRUM
      DO M=1,NFRE
        DO IJ=KIJS,KIJL
          BTH0(IJ,M) = 0.0_JWRB
        ENDDO
        DO K=1,NANG
          DO IJ=KIJS,KIJL
            BTH(IJ,K,M) = 0.0_JWRB
          ENDDO
        ENDDO
      ENDDO

      DO M=1,NFRE
        DO K=1,NANG
          ! integrates in directional sector
          DO K2=1,NSDSNTH*2+1
            KK=INDICESSAT(K,K2)
            DO IJ=KIJS,KIJL
              BTH(IJ,K,M) = BTH(IJ,K,M) + SATWEIGHTS(K,K2)*FL1(IJ,KK,M)
            ENDDO
          ENDDO
          DO IJ=KIJS,KIJL
            BTH(IJ,K,M)=BTH(IJ,K,M)*FACSAT(IJ,M)
            BTH0(IJ,M)=MAX(BTH0(IJ,M), BTH(IJ,K,M))
          ENDDO
        ENDDO
      ENDDO


      ! SATURATION TERM

      DO  M=1,NFRE
        ZCOEF = SSDSC2_SIG(M)*SSDSC6
        DO IJ=KIJS,KIJL
          TEMP1(IJ,M) = ZCOEF * (MAX(0._JWRB, BTH0(IJ,M)*TMP03-SSDSC4))**IPSAT
        ENDDO
      ENDDO
      DO  M=1,NFRE
        ZCOEFM1 = SSDSC2_SIG(M)*SSDSC6M1
        DO K=1,NANG
          DO IJ=KIJS,KIJL
            D(IJ,K,M) = TEMP1(IJ,M) + ZCOEFM1 * (MAX(0._JWRB, BTH(IJ,K,M)*TMP03-SSDSC4))**IPSAT
          ENDDO
        ENDDO
      ENDDO


      ! CUMULATIVE TERM
      IF (SSDSC3 /= 0.0_JWRB) THEN

        NANGD=NANG/2
        XLOGDFRTH=LOG(FRATIO)*DELTH
!       l(k,th)=1/(2*pi²)= the breaking crest density
        BRLAMBDA=BRKPBCOEF/(2.0_JWRB*ZPI**2)
        TMP02 = SSDSC3*BRLAMBDA

!       NDIKCUMUL is the  integer difference in frequency bands
!       between the "large breakers" and short "wiped-out waves"
!!! wrong !!???        NDIKCUMUL = NINT(SSDSBRF1/(FRATIO-1.))
        NDIKCUMUL = NINT(-LOG(SSDSBRF1)/LOG(FRATIO))

        DO KK=0,NANGD
          COSDTH(KK)=COS(KK*DELTH)
        ENDDO

        DO M=1,NFRE
          DO IJ=KIJS,KIJL
            C_(IJ,M)=SIG(M)/WAVNUM(IJ,M)
            C_C(IJ,M)=C_(IJ,M)**2
            DSIP(IJ,M)=TMP02*SIG(M)*XLOGDFRTH/CGROUP(IJ,M) !  coef*dtheta*dk = coef*dtheta*dsigma/cg
          ENDDO
        ENDDO

        DO M2=1,NFRE-NDIKCUMUL
          DO IJ=KIJS,KIJL
            IF (BTH0(IJ,M2) > SDSBR) THEN
              TEMP1(IJ,M2)=1.0_JWRB
            ELSE
              TEMP1(IJ,M2)=0.0_JWRB
            ENDIF
          ENDDO
        ENDDO
        DO M2=1,NFRE-NDIKCUMUL
          DO K2=1,NANG
            DO IJ=KIJS,KIJL
              SCUMUL(IJ,K2,M2)=TEMP1(IJ,M2)*(MAX(SQRT(BTH(IJ,K2,M2))-EPSR,0.0_JWRB))**2
            ENDDO
          ENDDO
        ENDDO

        DO M=NDIKCUMUL+1,NFRE
          DO K=1,NANG
            DO IJ=KIJS,KIJL
              RENEWALFREQ(IJ,K,M)=0.0_JWRB
            ENDDO
          ENDDO
        ENDDO


        DO M=NDIKCUMUL+1,NFRE

          IF (M-NDIKCUMUL >= 3) THEN
            DO IJ=KIJS,KIJL
              TRPZ_DSIP(IJ,1)=0.5_JWRB*DSIP(IJ,1)
            ENDDO
            DO M2=2,M-NDIKCUMUL-1
              DO IJ=KIJS,KIJL
                TRPZ_DSIP(IJ,M2)=DSIP(IJ,M2)
              ENDDO
            ENDDO
            DO IJ=KIJS,KIJL
              TRPZ_DSIP(IJ,M-NDIKCUMUL)=0.5_JWRB*DSIP(IJ,M-NDIKCUMUL)
            ENDDO
          ELSE
            DO M2=1,M-NDIKCUMUL
              DO IJ=KIJS,KIJL
                TRPZ_DSIP(IJ,M2)=DSIP(IJ,M2)
              ENDDO
            ENDDO
          ENDIF

          DO M2=1,M-NDIKCUMUL
            DO KK=0,NANGD
              DO IJ=KIJS,KIJL
                WCUMUL(IJ,KK,M2)=SQRT(ABS(C_C(IJ,M)+C_C(IJ,M2)-2.0_JWRB*C_(IJ,M)*C_(IJ,M2)*COSDTH(KK)))*TRPZ_DSIP(IJ,M2)
              ENDDO
            ENDDO
          ENDDO

          DO K=1,NANG
          ! Correction of saturation level for shallow-water kinematics
          ! Cumulative effect based on lambda   (breaking probability is
          ! the expected rate of sweeping by larger breaking waves)

            DO K2=1,NANG
              KKD(K2)=ABS(K2-K)
              IF ( KKD(K2) > NANGD) KKD(K2)=KKD(K2)-NANGD
            ENDDO

            DO M2=1,M-NDIKCUMUL
              DO K2=1,NANG
                KK=KKD(K2)
                DO IJ=KIJS,KIJL
                  ! Integrates over frequencies M2 and directions K2 to 
                  ! Integration is performed from M2=1 to a frequency lower than M: IK-NDIKCUMUL
                  RENEWALFREQ(IJ,K,M)=RENEWALFREQ(IJ,K,M)+ WCUMUL(IJ,KK,M2)*SCUMUL(IJ,K2,M2)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO

        DO M=NDIKCUMUL+1,NFRE
          DO K=1,NANG
            DO IJ=KIJS,KIJL
              D(IJ,K,M)= D(IJ,K,M) + RENEWALFREQ(IJ,K,M)
            ENDDO
          ENDDO
        ENDDO
      ENDIF


!     WAVE-TURBULENCE INTERACTION TERM
      IF (SSDSC5 /= 0.0_JWRB) THEN
        TMP01 = 2._JWRB*SSDSC5/G
        DO IJ=KIJS,KIJL
          FACTURB(IJ) = TMP01*RAORW(IJ)*UFRIC(IJ)*UFRIC(IJ)
        ENDDO
        DO M=1, NFRE
          DO IJ=KIJS,KIJL
            FACWTRB(IJ,M) = SIG(M)*WAVNUM(IJ,M)*FACTURB(IJ)
          ENDDO
          DO K=1,NANG
            DO IJ=KIJS,KIJL
              D(IJ,K,M)= D(IJ,K,M)- FACWTRB(IJ,M)*COSWDIF(IJ,K)
            ENDDO
          ENDDO
        ENDDO
      ENDIF


      ! ADD ALL CONTRIBUTIONS TO SOURCE TERM
      DO  M=1, NFRE
        DO K=1, NANG
          DO IJ=KIJS,KIJL
            SL(IJ,K,M) = SL(IJ,K,M)+D(IJ,K,M)*FL1(IJ,K,M)
            FLD(IJ,K,M) = FLD(IJ,K,M)+D(IJ,K,M)
          ENDDO
        ENDDO
      ENDDO

      IF (LHOOK) CALL DR_HOOK('SDISSIP_ARD',1,ZHOOK_HANDLE)

      END SUBROUTINE SDISSIP_ARD
