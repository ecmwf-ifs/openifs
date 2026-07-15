! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE MSWELL (KIJS, KIJL, IFROMIJ, JFROMIJ, NXS, NXE, NYS, NYE, FIELDG, FL1)
! ----------------------------------------------------------------------

!**** *MSWELL* - MAKES START SWELL FIELDS FOR WAMODEL.

!     J. BIDLOT     ECMWF    

!*    PURPOSE.
!     --------

!       TO GENERATE WAMODEL START FIELDS.

!**   INTERFACE.
!     ----------

!   *CALL* *MSWELL (KIJS, KIJL, BLK2LOC, NXS, NXE, NYS, NYE, FIELDG, FL1)
!      *BLK2LOC*             POINTERS FROM LOCAL GRID POINTS TO 2-D MAP
!      *NXS:NXE*  FIRST DIMENSION OF FIELDG
!      *NYS:NYE*  SECOND DIMENSION OF FIELDG
!      *FIELDG*   INPUT FORCING FIELDS ON THE WAVE MODEL GRID
!      *FL1*      REAL      2-D SPECTRUM FOR EACH GRID POINT 

!     METHOD.
!     -------

!       NLOC SWELL SYSTEMS ARE SPECIFIED (SEE BELOW), CENTERED AT
!       DIFFERENT (XLON0, YLAT0) COORDINATES.
!
!       EACH SWELL SPECTRUM AT (XLON0, YLAT0) IS SPECIFIED AS
!       F(OMEGA,THETA)=F1(OMEGA)*Q(THETA)
!       WHERE
!
!       F1(OMEGA)=(N+1)*m0*(OMEGA_p**N/OMEGA**(N+1))*EXP[-((N+1)/N)*(OMEGA_p/OMEGA)**N]
!       WITH
!          OMEGA=2PI*f ; f: wave frequency
!          m0 the zero spectral moment = Hs**2/(2PI) ; Hs the wave height (to be specified) 
!          OMEGA_p=2PI*f_p  ; f_p the peak frequency 
!          N=5
!
!       Q(THETA)=(8/3PI)*COS(THETA-THETA0)**4 for ABS(THETA-THATA0) >  PI/2 
!       Q(THETA)=0                            for ABS(THETA-THATA0) <= PI/2 
!       WITH
!           THETA the wave direction
!           THETA0 the mean wave direction of the swell system (to be specified)
!
!       THE SPECTRUM AT EACH (XLON0,YLAT0) is SPREAD ON THE SPHERE BY MULTIPLYING
!       F(OMEGA,THETA) BY SPEAD(XLON,YLAT)
!       WHERE
!       SPREAD(XLON,YLAT)=EXP[-2*DIST/L}]
!       WITH
!          DIST THE DISTANCE BETWEEN (XLON0,YLAT0)  and (XLON,YLAT)
!          L the correlation distance of the initial disturbance to be specified (in m)
!
!       NOTE THAT XLON,YLAT,XLON0,YLAT0 WILL SPECIFIED IN DEGREE BUT ARE CONVERTED 
!       IN RADIAN FOR CALCULATIONS.
!

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------

!       I. V. LAVRENOV, 2003: WIND-WAVES IN OCEANS, DYNAMICS AND NUMERICAl SIMILATIONS
!       SPRINGER-VERLAG, 376pp.
!       see pages 52-53

!       NONE.

! ----------------------------------------------------------------------

!
!       F1(OMEGA)=(N+1)*m0*(OMEGA_p**N/OMEGA**(N+1))*EXP[-((N+1)/N)*(OMEGA_p/OMEGA)**N]
!       WITH
!          OMEGA=2PI*f ; f: wave frequency
!          m0 the zero spectral moment = Hs**2/(2PI) ; Hs the wave height (to be specified) 
!          OMEGA_p=2PI*f_p  ; f_p the peak frequency 
!          N=5
!
!       Q(THETA)=(8/3PI)*COS(THETA-THETA0)**4 for ABS(THETA-THATA0) >  PI/2 
!       Q(THETA)=0                            for ABS(THETA-THATA0) <= PI/2 
!       WITH
!           THETA the wave direction
!           THETA0 the mean wave direction of the swell system (to be specified)
!
!       THE SPECTRUM AT EACH (XLON0,YLAT0) is SPREAD ON THE SPHERE BY MULTIPLYING
!       F(OMEGA,THETA) BY SPEAD(XLON,YLAT)
!       WHERE
!       SPREAD(XLON,YLAT)=EXP[-2*DIST/L}]
!       WITH
!          DIST THE DISTANCE BETWEEN (XLON0,YLAT0)  and (XLON,YLAT)
!          L the correlation distance of the initial disturbance to be specified (in m)
!
!       NOTE THAT XLON,YLAT,XLON0,YLAT0 WILL SPECIFIED IN DEGREE BUT ARE CONVERTED 
!       IN RADIAN FOR CALCULATIONS.
!

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------

!       I. V. LAVRENOV, 2003: WIND-WAVES IN OCEANS, DYNAMICS AND NUMERICAl SIMILATIONS
!       SPRINGER-VERLAG, 376pp.
!       see pages 52-53

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWDRVTYPE  , ONLY : WVGRIDLOC, FORCING_FIELDS

      USE YOWFRED  , ONLY : FR       ,TH
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : ZPI      ,RAD      ,R       ,ZMISS
      USE YOWSPHERE, ONLY : SPHERICAL_COORDINATE_DISTANCE

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      INTEGER(KIND=JWIM), DIMENSION(KIJS:KIJL), INTENT(IN) :: IFROMIJ, JFROMIJ
      INTEGER(KIND=JWIM), INTENT(IN) :: NXS, NXE, NYS, NYE
      TYPE(FORCING_FIELDS), INTENT(IN) :: FIELDG
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL, NANG, NFRE), INTENT(OUT) :: FL1


      INTEGER(KIND=JWIM), PARAMETER :: NLOC=4 ! TOTAL NUMBER OF SWELL SYSTEMS
      INTEGER(KIND=JWIM) :: IJ, K, M, ILOC, IX, JY
      INTEGER(KIND=JWIM) :: NSP
      INTEGER(KIND=JWIM), DIMENSION(NLOC) :: KLOC, MLOC

      REAL(KIND=JWRB) :: CQ0, COSDIRMAX, COSDIR
      REAL(KIND=JWRB) :: OMEGA, OMEGADIFFMIN, OMEGADIFF
      REAL(KIND=JWRB) :: E0, CEX, CS0, S0, SPRD
      REAL(KIND=JWRB), DIMENSION(NLOC) :: H0, OMEGAP, XL
      REAL(KIND=JWRB), DIMENSION(NLOC) :: THETA0, COSLAT2
      REAL(KIND=JWRB), DIMENSION(NANG) :: Q0
      REAL(KIND=JWRB), DIMENSION(NLOC,NANG,NFRE):: FL0

      REAL(KIND=JWRU) :: XLO, YLA, DIST
      REAL(KIND=JWRU), DIMENSION(NLOC) :: YLAT0, XLON0

!----------------------------------------------------------------------

        CQ0=16.0_JWRB/(3*ZPI)

        DO M=1,NFRE
          DO K=1,NANG
            DO IJ = KIJS, KIJL
              FL1(IJ,K,M)=0.0_JWRB
            ENDDO
          ENDDO
        ENDDO

!       DEFINE THE SWELL SYSTEMS
        H0(1)=2.0_JWRB
        THETA0(1)=135.0_JWRB
        OMEGAP(1)=0.3117_JWRB
        XL(1)=250000.0_JWRB
        YLAT0(1)=47.0_JWRU
        XLON0(1)=165.0_JWRU

        H0(2)=2.0_JWRB
        THETA0(2)=90.0_JWRB
        OMEGAP(2)=0.3117_JWRB
        XL(2)=200000.0_JWRB
        YLAT0(2)=-50.0_JWRU
        XLON0(2)=20.0_JWRU

        H0(3)=2.0_JWRB
        THETA0(3)=180.0_JWRB
        OMEGAP(3)=0.3117_JWRB
        XL(3)=200000.0_JWRB
        YLAT0(3)=35.0_JWRU
        XLON0(3)=331.0_JWRU

        H0(4)=2.0_JWRB
        THETA0(4)=45.0_JWRB
        OMEGAP(4)=0.3117_JWRB
        XL(4)=150000.0_JWRB
        YLAT0(4)=52.0_JWRU
        XLON0(4)=329.0_JWRU

        NSP=5


        DO ILOC=1,NLOC
          THETA0(ILOC)=RAD*THETA0(ILOC)
          COSLAT2(ILOC)=COS(RAD*YLAT0(ILOC))**2
        ENDDO

        DO ILOC=1,NLOC
          KLOC(ILOC)=1
          COSDIRMAX=COS(TH(K)-THETA0(1))
          DO K=2,NANG
            COSDIR=COS(TH(K)-THETA0(ILOC))
            IF (COSDIRMAX < COSDIR) THEN
              KLOC(ILOC)=K
              COSDIRMAX=COSDIR
            ENDIF
          ENDDO
          MLOC(ILOC)=1
          OMEGA=ZPI*FR(1)
          OMEGADIFFMIN=ABS(OMEGA-OMEGAP(ILOC))
          DO M=2,NFRE
            OMEGA=ZPI*FR(M)
            OMEGADIFF=ABS(OMEGA-OMEGAP(ILOC))
            IF (OMEGADIFF < OMEGADIFFMIN) THEN
              MLOC(ILOC)=M
              OMEGADIFFMIN=OMEGADIFF
            ENDIF
          ENDDO
        ENDDO

        DO ILOC=1,NLOC
          DO K=1,NANG
            COSDIR=COS(TH(K)-THETA0(ILOC))
            IF (COSDIR > 0.0_JWRB) THEN
              Q0(K) = CQ0*COSDIR**4
            ELSE
              Q0(K) = 0. 
            ENDIF
          ENDDO
          E0=H0(ILOC)**2/16.0_JWRB
          CEX=REAL(NSP+1,JWRB)/REAL(NSP,JWRB)
          CS0=(NSP+1)*E0*OMEGAP(ILOC)**NSP
          DO M=1,NFRE
            OMEGA=ZPI*FR(M)
            S0=(CS0/OMEGA**(NSP+1))*EXP(-CEX*(OMEGAP(ILOC)/OMEGA)**NSP)
            IF (S0 < 0.001_JWRB) S0=0.0_JWRB
            DO K=1,NANG
              FL0(ILOC,K,M)= Q0(K)*S0
            ENDDO
          ENDDO
        ENDDO

        DO IJ = KIJS, KIJL
          IX = IFROMIJ(IJ)
          JY = JFROMIJ(IJ)
          XLO=FIELDG%XLON(IX,JY)
          YLA=FIELDG%YLAT(IX,JY)
          IF (YLA == ZMISS .OR. XLO == ZMISS) CYCLE 

          DO ILOC=1,NLOC
            CALL SPHERICAL_COORDINATE_DISTANCE(XLON0(ILOC),XLO,YLAT0(ILOC),YLA,DIST)
            DIST=2*R*DIST/XL(ILOC)
            IF (DIST < 10.0_JWRU) THEN
              SPRD=EXP(-DIST)
              DO M=1,NFRE
                DO K=1,NANG
                  FL1(IJ,K,M)=FL1(IJ,K,M)+FL0(ILOC,K,M)*SPRD 
                ENDDO
              ENDDO
            ENDIF
          ENDDO

        ENDDO

      END SUBROUTINE MSWELL
