! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE WEFLUX (KIJS, KIJL, FL1, CGROUP,         &
     &                   NFRE, NANG, DFIM, DELTH,         &
     &                   COSTH, SINTH,                    &
     &                   WEFMAG, WEFDIR)
!
! ----------------------------------------------------------------------
!
!**** *WEFLUX* - COMPUTATION OF WAVE ENERGY FLUX AT EACH GRID POINT


!*    PURPOSE.
!     --------

!       COMPUTE WAVE ENERGY FLUX MAGNITUDE AND MEAN DIRECTION AT EACH GRID POINT.

!**   INTERFACE.
!     ----------

!       *CALL* *WEFLUX* (KIJS, KIJL, FL1, CGROUP,
!                        NFRE, NANG, DFIM, DELTH,
!                        COSTH, SINTH,
!                        WEFMAG, WEFDIR)
!              *KIJS*    - INDEX OF FIRST GRIDPOINT
!              *KIJL*    - INDEX OF LAST GRIDPOINT
!              *FL1*     - SPECTRUM.
!              *CGROUP*  - GROUP VELOCITIES
!              *NFRE*    - NUMBER OF FREQUENCIES
!              *NANG*    - NUMBER OF DIRECTIONS
!              *DFIM*    - FREQUENCY-DIRECTION INCREMENT
!              *DELTH*   - DIRECTIONAL INCREMENT
!              *COSTH*   - COS(THETA)
!              *SINTH*   - SIN(THETA)
!              *WEFMAG*  - WAVE ENERGY FLUX MAGNITUDE (W/m) 
!              *WEFDIR*  - WAVE ENERGY FLUX MEAN DIRECTION (radian oceanographic convention) 
!

!     METHOD.
!     -------

!       NONE.

!     EXTERNALS.
!     ----------

!       NONE.
!
!     REFERENCE.
!     ----------
!
!       NONE.
!
! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
 
      USE YOWPCONS , ONLY : G        ,ZPI      ,ROWATER   ,EPSMIN
      USE YOWFRED  , ONLY : FRTAIL

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
 
! ----------------------------------------------------------------------
 
      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(IN) :: FL1
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE), INTENT(IN) :: CGROUP
      INTEGER(KIND=JWIM), INTENT(IN) :: NFRE, NANG
      REAL(KIND=JWRB), INTENT(IN) :: DELTH
      REAL(KIND=JWRB), DIMENSION(NFRE), INTENT(IN) :: DFIM
      REAL(KIND=JWRB), DIMENSION(NANG), INTENT(IN) :: COSTH, SINTH
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(OUT) :: WEFMAG, WEFDIR


      INTEGER(KIND=JWIM) :: IJ, M, K

      REAL(KIND=JWRB) :: ROG, FCG, DELT
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: TEMP, TEMPX, TEMPY
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: WEFX, WEFY 
!
! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('WEFLUX',0,ZHOOK_HANDLE)

!*    1. INITIALISE MEAN FREQUENCY ARRAY AND TAIL FACTOR.
!        ------------------------------------------------

      DO IJ=KIJS,KIJL
        WEFMAG(IJ)= 0.0_JWRB
        WEFX(IJ)  = 0.0_JWRB
        WEFY(IJ)  = 0.0_JWRB
      ENDDO

      ROG = ROWATER*G  
      DELT = FRTAIL*DELTH*G/(2.0_JWRB*ZPI)
!
!*    2. INTEGRATE OVER FREQUENCIES AND DIRECTIONS.
!        ------------------------------------------

      DO M=1,NFRE
         K=1
         DO IJ=KIJS,KIJL
           FCG = FL1(IJ,K,M)*CGROUP(IJ,M)
           TEMP(IJ) = FCG
           TEMPX(IJ) = FCG*SINTH(K) 
           TEMPY(IJ) = FCG*COSTH(K) 
         ENDDO
         DO K=2,NANG
            DO IJ=KIJS,KIJL
              FCG = FL1(IJ,K,M)*CGROUP(IJ,M)
              TEMP(IJ) = TEMP(IJ)+FCG
              TEMPX(IJ) = TEMPX(IJ)+FCG*SINTH(K)
              TEMPY(IJ) = TEMPY(IJ)+FCG*COSTH(K)
            ENDDO
         ENDDO
         DO IJ=KIJS,KIJL
           WEFMAG(IJ)  = WEFMAG(IJ)+DFIM(M)*TEMP(IJ)
           WEFX(IJ)  = WEFX(IJ)+DFIM(M)*TEMPX(IJ)
           WEFY(IJ)  = WEFY(IJ)+DFIM(M)*TEMPY(IJ)
         ENDDO
      ENDDO


!*    3. ADD TAIL CORRECTION
!        -------------------

       K=1
       DO IJ=KIJS,KIJL
         FCG = FL1(IJ,K,NFRE)
         TEMP(IJ) = FCG
         TEMPX(IJ) = FCG*SINTH(K) 
         TEMPY(IJ) = FCG*COSTH(K) 
       ENDDO
       DO K=2,NANG
          DO IJ=KIJS,KIJL
            FCG = FL1(IJ,K,NFRE)
            TEMP(IJ) = TEMP(IJ)+FCG
            TEMPX(IJ) = TEMPX(IJ)+FCG*SINTH(K)
            TEMPY(IJ) = TEMPY(IJ)+FCG*COSTH(K)
          ENDDO
       ENDDO
       DO IJ=KIJS,KIJL
         WEFMAG(IJ) = WEFMAG(IJ)+DELT*TEMP(IJ)
         WEFX(IJ)  = WEFX(IJ)+DELT*TEMPX(IJ)
         WEFY(IJ)  = WEFY(IJ)+DELT*TEMPY(IJ)
       ENDDO

!*    4. MULTIPLY MAGNITUDE BY ROG:
!        -------------------------
      DO IJ=KIJS,KIJL
        WEFMAG(IJ)  = ROG*WEFMAG(IJ)
      ENDDO

!*    5. FLUX MEAN DIRECTION
!        -------------------

      DO IJ=KIJS,KIJL
        IF (WEFY(IJ) == 0.0_JWRB) WEFY(IJ) = EPSMIN
      ENDDO

      DO IJ=KIJS,KIJL
        WEFDIR(IJ)  = ATAN2(WEFX(IJ),WEFY(IJ)) 
      ENDDO

      DO IJ=KIJS,KIJL
        IF (WEFDIR(IJ) < 0.0_JWRB) WEFDIR(IJ) = WEFDIR(IJ) + ZPI
      ENDDO

      IF (LHOOK) CALL DR_HOOK('WEFLUX',1,ZHOOK_HANDLE)

      END SUBROUTINE WEFLUX
