! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

      SUBROUTINE INIT_SDISS_ARDH
! ----------------------------------------------------------------------

!**** *INIT_SDISS_ARDH* - INITIALISATION FOR SDISS_ARD

!     LOTFI AOUF       METEO FRANCE 2013
!     FABRICE ARDHUIN  IFREMER  2013


!*    PURPOSE.
!     --------

!**   INTERFACE.
!     ----------

!       *CALL* *INIT_SDISS_ARDH

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

      USE YOWFRED  , ONLY : FR, TH, FRATIO, DELTH, COSTH, SINTH
      USE YOWPCONS , ONLY : RAD     ,G        ,ZPI
      USE YOWPARAM , ONLY : NANG    ,NFRE
      USE YOWPHYS  , ONLY : SDSBR   ,ISDSDTH ,ISB     ,IPSAT    ,      &
&                  SSDSC2  , SSDSC4,  SSDSC6,  MICHE, SSDSC3, SSDSBRF1,&
&                  BRKPBCOEF  ,SSDSC5, NSDSNTH, NDIKCUMUL,             &
&                  INDICESSAT, SATWEIGHTS, CUMULW

      USE YOWSHAL  , ONLY : NDEPTH  ,TFAK , TCGOND

      USE YOMHOOK  , ONLY : LHOOK   ,DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM) :: JD, K, M, I_INT, J_INT, M2, KK, NANGD

      REAL(KIND=JWRB) :: TPIINV, TMP01, TMP02
      REAL(KIND=JWRB) :: DELTH_TRUNC, DELTH_LOC
      REAL(KIND=JWRB) :: DTURB
      REAL(KIND=JWRB) :: XLOGDFRTH
      REAL(KIND=JWRB) :: BRLAMBDA
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(0:NANG/2) :: COSDTH
      REAL(KIND=JWRB), DIMENSION(NFRE) :: SIG, C_, C_C, CGM1, DSIP, TRPZ_DSIP 

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('INIT_SDISS_ARDH',0,ZHOOK_HANDLE)

      TPIINV = 1.0_JWRB/ZPI

      NANGD=NANG/2

      XLOGDFRTH=LOG(FRATIO)*DELTH

!     COMPUTE SATWEIGHTS

!     l(k,th)=1/(2*pi²)= the breaking crest density
      BRLAMBDA=BRKPBCOEF/(2.0_JWRB*ZPI**2)

      TMP02 = SSDSC3*BRLAMBDA

      NSDSNTH  = MIN(NINT(ISDSDTH*RAD/(DELTH)),NANGD-1)
      DELTH_TRUNC=(TH(1)+ISDSDTH*RAD)-(TH(1+NSDSNTH)-0.5_JWRB*DELTH)
      DELTH_TRUNC=MAX(0.0_JWRB, MIN(DELTH_TRUNC,DELTH))

      IF (ALLOCATED(INDICESSAT)) DEALLOCATE(INDICESSAT)
      ALLOCATE(INDICESSAT(NANG,NSDSNTH*2+1))
      IF (ALLOCATED(SATWEIGHTS)) DEALLOCATE(SATWEIGHTS)
      ALLOCATE(SATWEIGHTS(NANG,NSDSNTH*2+1))

      DO K=1,NANG
        DO I_INT=K-NSDSNTH, K+NSDSNTH
          J_INT=I_INT
          IF (I_INT < 1)  J_INT=I_INT+NANG
          IF (I_INT > NANG) J_INT=I_INT-NANG
          INDICESSAT(K,I_INT-(K-NSDSNTH)+1)=J_INT

          IF (I_INT == K-NSDSNTH .OR. I_INT == K+NSDSNTH) THEN
            DELTH_LOC=DELTH_TRUNC
          ELSE
            DELTH_LOC=DELTH
          ENDIF
          SATWEIGHTS(K,I_INT-(K-NSDSNTH)+1)=DELTH_LOC*COS(TH(K)-TH(J_INT))**ISB
        END DO
      END DO

!     COMPUTE CUMULW
      IF (ALLOCATED(CUMULW)) DEALLOCATE(CUMULW)
      ALLOCATE(CUMULW(NDEPTH,0:NANGD,NFRE,NFRE))

      IF (SSDSC3 /= 0.0_JWRB) THEN
!       NDIKCUMUL is the  integer difference in frequency bands
!       between the "large breakers" and short "wiped-out waves"
!!! wrong !!???        NDIKCUMUL = NINT(SSDSBRF1/(FRATIO-1.))
        NDIKCUMUL = NINT(-LOG(SSDSBRF1)/LOG(FRATIO))

        DO KK=0,NANGD
          COSDTH(KK)=COS(KK*DELTH)
        ENDDO

        DO M=1,NFRE
          SIG(M) = ZPI*FR(M)
        ENDDO

        DO JD=1,NDEPTH

          DO M=1,NFRE
            C_(M)=SIG(M)/TFAK(JD,M)
            CGM1(M)=1.0_JWRB/TCGOND(JD,M)
          ENDDO

          DO M=1,NFRE
            C_C(M)=C_(M)*C_(M)
            DSIP(M)=TMP02*SIG(M)*XLOGDFRTH*CGM1(M) !  coef*dtheta*dk = coef*dtheta*dsigma/cg
          ENDDO

          DO M=NDIKCUMUL+1,NFRE

            IF (M-NDIKCUMUL >= 3) THEN
              TRPZ_DSIP(1)=0.5_JWRB*DSIP(1)
              DO M2=2,M-NDIKCUMUL-1
                TRPZ_DSIP(M2)=DSIP(M2)
              ENDDO
              TRPZ_DSIP(M-NDIKCUMUL)=0.5_JWRB*DSIP(M-NDIKCUMUL)
            ELSE
              DO M2=1,M-NDIKCUMUL
                TRPZ_DSIP(M2)=DSIP(M2)
              ENDDO
            ENDIF

            DO M2=1,M-NDIKCUMUL
              DO KK=0,NANGD
                CUMULW(JD,KK,M2,M)=SQRT(ABS(C_C(M)+C_C(M2)-2.0_JWRB*C_(M)*C_(M2)*COSDTH(KK)))*TRPZ_DSIP(M2)
              ENDDO 
            ENDDO
          ENDDO

        ENDDO ! JD

      ELSE
        CUMULW(:,:,:,:) = 0.0_JWRB
      ENDIF

      IF (LHOOK) CALL DR_HOOK('INIT_SDISS_ARDH',1,ZHOOK_HANDLE)

      END SUBROUTINE INIT_SDISS_ARDH
