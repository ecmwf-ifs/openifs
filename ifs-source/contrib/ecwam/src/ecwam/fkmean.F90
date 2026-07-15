! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE FKMEAN (KIJS, KIJL, FL1, WAVNUM,   &
     &                   EM, FM1, F1, AK, XK)

! ----------------------------------------------------------------------

!**** *FKMEAN* - COMPUTATION OF MEAN FREQUENCIES AT EACH GRID POINT
!                AND MEAN WAVE NUMBER (based in  sqrt(k)*F moment) .
!                COMPUTATION OF THE MEAN WAVE ENERGY WAS ALSO
!                ADDED SUCH THAT A CALL TO FKMEAN DOES NOT NEED


!*    PURPOSE.
!     --------

!       COMPUTE MEAN FREQUENCIES AND WAVE NUMBER AT EACH GRID POINT.

!**   INTERFACE.
!     ----------

!       *CALL* *FKMEAN (KIJS, KIJL, FL1, WAVNUM, EM, FM1, F1, AK, XK)*
!              *KIJS*    - LOCAL INDEX OF FIRST GRIDPOINT
!              *KIJL*    - LOCAL INDEX OF LAST GRIDPOINT
!              *FL1*     - SPECTRUM.
!              *WAVNUM*  - WAVE NUMBER.
!              *EM*      - MEAN WAVE ENERGY
!              *FM1*     - MEAN WAVE FREQUENCY BASED ON (1/f)*FL1 INTEGRATION
!              *F1*      - MEAN WAVE FREQUENCY BASED ON f*FL1 INTEGRATION
!              *AK*      - MEAN WAVE NUMBER  BASED ON sqrt(1/k)*FL1 INTGRATION
!                          ONLY FOR SHALLOW WATER RUNS.
!!!                        AK IS STILL NEEDED IN SNONLIN !!!!
!!!                        IF THE OLD FORMULATION IS USED.
!              *XK*      - MEAN WAVE NUMBER  BASED ON sqrt(k)*FL1 INTEGRATION
!                          ONLY FOR SHALLOW WATER RUNS.

!     METHOD.
!     -------

!       NONE.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : FR       ,DFIM     ,DFIMOFR  ,DFIMFR   ,    &
     &            DELTH    ,WETAIL   ,FRTAIL   ,WP1TAIL
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : G        ,ZPI      ,EPSMIN

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL

      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(IN) :: FL1
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE), INTENT(IN) :: WAVNUM 

      REAL(KIND=JWRB),DIMENSION(KIJS:KIJL), INTENT(OUT) :: EM, FM1, F1, AK, XK


      INTEGER(KIND=JWIM) :: IJ, M, K
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB) :: DELT25, COEFM1, COEF1, COEFA, COEFX, SQRTK
      REAL(KIND=JWRB),DIMENSION(KIJS:KIJL) :: TEMPA, TEMPX,  TEMP2

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('FKMEAN',0,ZHOOK_HANDLE)


!*    1. INITIALISE MEAN FREQUENCY ARRAY AND TAIL FACTOR.
!        ------------------------------------------------

      DO IJ=KIJS,KIJL
        EM(IJ) = EPSMIN
        FM1(IJ)= EPSMIN
        F1(IJ) = EPSMIN
        AK(IJ) = EPSMIN
        XK(IJ) = EPSMIN
      ENDDO

      DELT25 = WETAIL*FR(NFRE)*DELTH
      COEFM1 = FRTAIL*DELTH
      COEF1 = WP1TAIL*DELTH*FR(NFRE)**2
      COEFA = COEFM1*SQRT(G)/ZPI
      COEFX = COEF1*(ZPI/SQRT(G))

!*    2. INTEGRATE OVER FREQUENCIES AND DIRECTIONS.
!        ------------------------------------------

!*    2.2 SHALLOW WATER INTEGRATION.
!         --------------------------

        DO M=1,NFRE
          DO IJ=KIJS,KIJL
            SQRTK=SQRT(WAVNUM(IJ,M))
            TEMPA(IJ) = DFIM(M)/SQRTK
            TEMPX(IJ) = SQRTK*DFIM(M)
          ENDDO
          K=1
          DO IJ=KIJS,KIJL
            TEMP2(IJ) = FL1(IJ,K,M) 
          ENDDO
          DO K=2,NANG
            DO IJ=KIJS,KIJL
              TEMP2(IJ) = TEMP2(IJ)+FL1(IJ,K,M)
            ENDDO
          ENDDO
          DO IJ=KIJS,KIJL
            EM(IJ) = EM(IJ)+DFIM(M)*TEMP2(IJ)
            FM1(IJ)= FM1(IJ)+DFIMOFR(M)*TEMP2(IJ)
            F1(IJ) = F1(IJ)+DFIMFR(M)*TEMP2(IJ)
            AK(IJ) = AK(IJ)+TEMPA(IJ)*TEMP2(IJ)
            XK(IJ) = XK(IJ)+TEMPX(IJ)*TEMP2(IJ)
          ENDDO
        ENDDO

!*      ADD TAIL CORRECTION TO MEAN FREQUENCY AND
!*      NORMALIZE WITH TOTAL ENERGY.
        DO IJ=KIJS,KIJL
          EM(IJ) = EM(IJ)+DELT25*TEMP2(IJ)
          FM1(IJ) = FM1(IJ)+COEFM1*TEMP2(IJ)
          FM1(IJ) = EM(IJ)/FM1(IJ)
          F1(IJ) = F1(IJ)+COEF1*TEMP2(IJ)
          F1(IJ) = F1(IJ)/EM(IJ)
          AK(IJ) = AK(IJ)+COEFA*TEMP2(IJ)
          AK(IJ) = (EM(IJ)/AK(IJ))**2
          XK(IJ) = XK(IJ)+COEFX*TEMP2(IJ)
          XK(IJ) = (XK(IJ)/EM(IJ))**2
        ENDDO

      IF (LHOOK) CALL DR_HOOK('FKMEAN',1,ZHOOK_HANDLE)

      END SUBROUTINE FKMEAN
