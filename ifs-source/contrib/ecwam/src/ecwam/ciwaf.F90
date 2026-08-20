! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE CIWAF (KIJS, KIJL, CGROUP, CICOVER, CITHICK, CIWA)

! ----------------------------------------------------------------------

!**** *CIWAF* - COMPUTE SEA ICE WAVE ATTENUATION FACTORS DUE TO ICE FLOES
!               SCATTERING.

!*    PURPOSE.
!     --------

!       CIWAF COMPUTES SEA ICE WAVE ATTENUATION FACTORS

!**   INTERFACE.
!     ----------

!       *CALL* *CIWAF (KIJS, KIJL, CGROUP, CICOVER, CITHICK, CIWA)

!          *KIJS*     - INDEX OF FIRST POINT.
!          *KIJL*     - INDEX OF LAST POINT.
!          *CGROUP*   - GROUP SPEED.
!          *CICOVER*  - SEA ICE COVER.
!          *CITHICK*  - SEA ICE THICKNESS. 
!          *CIWA*     - SEA ICE WAVE ATTENUATION FACTOR. 

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!     REFERENCES.                                                       
!     -----------  

!     KOHOUT AND MEYLAN, 2008: JGR, 113, doi:10.1029/2007JC004424


! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : FR
      USE YOWICE   , ONLY : NICT   ,NICH     ,TICMIN   ,HICMIN   ,      &
     &              DTIC   ,DHIC   ,CIDEAC   ,CIBLOCK
      USE YOWPARAM , ONLY : NFRE
      USE YOWPCONS , ONLY : EPSMIN
      USE YOWSTAT  , ONLY : IDELT

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL 
      REAL(KIND=JWRB),DIMENSION(KIJS:KIJL,NFRE), INTENT(IN) :: CGROUP
      REAL(KIND=JWRB),DIMENSION(KIJS:KIJL), INTENT(IN) :: CICOVER
      REAL(KIND=JWRB),DIMENSION(KIJS:KIJL), INTENT(IN) :: CITHICK
      REAL(KIND=JWRB),DIMENSION(KIJS:KIJL,NFRE), INTENT(OUT) :: CIWA


      INTEGER(KIND=JWIM) :: ICM, I, MAXICM
      INTEGER(KIND=JWIM) :: M, IJ
      INTEGER(KIND=JWIM) :: IT, IT1, IH, IH1 

      REAL(KIND=JWRB) :: CIDMIN, CIDMAX, CIDMEAN, CIFRGL, CIFRGMT
      REAL(KIND=JWRB) :: SD, SN 
      REAL(KIND=JWRB) :: TW, X
      REAL(KIND=JWRB) :: A, B, C
      REAL(KIND=JWRB) :: CIDEAC_INT, WT, WT1, WH, WH1 
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB),DIMENSION(KIJS:KIJL) :: DINV
      REAL(KIND=JWRB),DIMENSION(KIJS:KIJL,NFRE) :: ALP 

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('CIWAF',0,ZHOOK_HANDLE)

!     following  Dumont et al. (2011), eqn (13):
      ! sea ice fragility
      CIFRGL=0.955_JWRB
      ! sea ice foes minimum size (m)
      CIDMIN=20.0_JWRB
      ! sea ice fragmentation number (>2)
      CIFRGMT=2.0_JWRB

      A=200.0_JWRB
      C=300.0_JWRB
!     the maximum number of fragmentation step is set
!     for the limit value for sea ice = 0 (CIDMAX=A)
!     to insure that CIDMEAN will be increasing with increasing sea ice cover
      MAXICM=INT(LOG(A/CIDMIN)/LOG(CIFRGMT))
      DO IJ=KIJS,KIJL
        IF(CITHICK(IJ) > 0.0_JWRB) THEN
          ! sea ice foes maxmimum size (m)
!!!       testing making it function of sea ice cover
          CIDMAX=A+C*CICOVER(IJ)
          ICM=MIN(INT(LOG(CIDMAX/CIDMIN)/LOG(CIFRGMT)),MAXICM)
          SN=0.0_JWRB
          SD=0.0_JWRB
           DO I=0,ICM
            X=(CIFRGMT**2*CIFRGL)**I
            SN=SN+X*CIDMAX/CIFRGMT**I
            SD=SD+X
          ENDDO
        CIDMEAN=SN/SD
        DINV(IJ)=1.0_JWRB/CIDMEAN
        ELSE
          DINV(IJ)=CIDMIN
        ENDIF
      ENDDO


      DO M=1,NFRE
        TW=1.0_JWRB/FR(M)
        IT=FLOOR((TW-TICMIN)/DTIC+1)
        IT=MAX(1,MIN(IT,NICT))
        IT1=IT+1
        IT1=MAX(1,MIN(IT1,NICT))
        WT1=MAX(MIN(1.0_JWRB,(TW-(TICMIN+(IT-1)*DTIC))/DTIC),0.0_JWRB)
        WT=1.0_JWRB-WT1
        DO IJ=KIJS,KIJL
          IF(CITHICK(IJ) > 0.0_JWRB) THEN
            IH=FLOOR((CITHICK(IJ)-HICMIN)/DHIC+1)
            IH=MAX(1,MIN(IH,NICH))
            IH1=IH+1
            IH1=MAX(1,MIN(IH1,NICH))
            WH1=MAX(MIN(1._JWRB,(CITHICK(IJ)-(HICMIN+(IH-1)*DHIC))/DHIC),0.0_JWRB)
            WH=1.0_JWRB-WH1
            CIDEAC_INT=WT*(WH*CIDEAC(IT,IH)+ WH1*CIDEAC(IT,IH1)) +      &
     &                WT1*(WH*CIDEAC(IT1,IH)+WH1*CIDEAC(IT1,IH1))
!!!            ALP(IJ,M)=CICOVER(IJ)*CIDEAC_INT*DINV(IJ)
            ALP(IJ,M)=CICOVER(IJ)*EXP(CIDEAC_INT)*DINV(IJ)
          ELSE
            ALP(IJ,M)=0.0_JWRB
          ENDIF
        ENDDO
      ENDDO

      DO M=1,NFRE
        DO IJ=KIJS,KIJL
          X=ALP(IJ,M)*CGROUP(IJ,M)*IDELT
          IF(X.LT.EPSMIN) THEN
            CIWA(IJ,M)=1.0_JWRB
          ELSE IF(CICOVER(IJ) > CIBLOCK) THEN
            CIWA(IJ,M)=0.0_JWRB
          ELSE
            CIWA(IJ,M)=EXP(-MIN(X,50.0_JWRB))
          ENDIF
        ENDDO
      ENDDO

      IF (LHOOK) CALL DR_HOOK('CIWAF',1,ZHOOK_HANDLE)

      END SUBROUTINE CIWAF
