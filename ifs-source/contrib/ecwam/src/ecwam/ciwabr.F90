! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE CIWABR (KIJS, KIJL, CICOVER, FL1, WAVNUM, CGROUP, CIWAB)

! ----------------------------------------------------------------------

!**** *CIWABR* - COMPUTE SEA ICE WAVE ATTENUATION FACTORS DUE TO ICE FLOES
!                BOTTOM FRICTION.

!*    PURPOSE.
!     --------

!       CIWABR COMPUTES SEA ICE WAVE ATTENUATION FACTORS DUE TO ICE FLOES
!              BOTTOM FRICTION.

!**   INTERFACE.
!     ----------

!       *CALL* *CIWABR (KIJS,KIJL,CICOVER,FL1,CIWAB)

!          *KIJS*     - INDEX OF FIRST POINT.
!          *KIJL*     - INDEX OF LAST POINT.
!          *CICOVER*  -SEA ICE COVER.
!          *FL1*      -ENERGY SPECTRUM. 
!          *CIWAB*    -SEA ICE WAVE ATTENUATION FACTOR DUE TO ICE FLOE BOTTOM FRICTION 

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!     REFERENCES.                                                       
!     -----------  

!     KOHOUT A., M. MEYLAN, D PLEW, 2011: ANNALS OF GLACIOLOGY, 2011. 


! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : FR       ,DFIM     , DELTH
      USE YOWICE   , ONLY : LICERUN  ,LMASKICE , CDICWA
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : G        ,ZPI      ,ZPI4GM2    ,EPSMIN
      USE YOWSTAT  , ONLY : IDELT

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL 
      REAL(KIND=JWRB),DIMENSION(KIJS:KIJL), INTENT(IN) :: CICOVER
      REAL(KIND=JWRB),DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(IN) :: FL1
      REAL(KIND=JWRB),DIMENSION(KIJS:KIJL,NFRE), INTENT(IN) :: WAVNUM, CGROUP 
      REAL(KIND=JWRB),DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(OUT) :: CIWAB


      INTEGER(KIND=JWIM) :: K, M, IJ
      REAL(KIND=JWRB) :: EWH 
      REAL(KIND=JWRB) :: X, ALP
      REAL(KIND=JWRB),DIMENSION(NFRE) :: XK2 
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('CIWABR',0,ZHOOK_HANDLE)

      IF ( .NOT. LICERUN .OR. LMASKICE ) THEN

        DO M=1,NFRE
          DO K=1,NANG
            DO IJ=KIJS,KIJL
              CIWAB(IJ,K,M)=1.0_JWRB
            ENDDO
          ENDDO
        ENDDO

      ELSE

          DO M=1,NFRE
            DO K=1,NANG
              DO IJ=KIJS,KIJL
                EWH=4.0_JWRB*SQRT(MAX(EPSMIN,FL1(IJ,K,M)*DFIM(M)))
                XK2(M)=WAVNUM(IJ,M)**2
                ALP=CDICWA*XK2(M)*EWH
                X=ALP*CGROUP(IJ,M)*IDELT
                CIWAB(IJ,K,M)=1.0_JWRB-CICOVER(IJ)*(1.0_JWRB-EXP(-MIN(X,50.0_JWRB)))
              ENDDO
            ENDDO
          ENDDO

      ENDIF

      IF (LHOOK) CALL DR_HOOK('CIWABR',1,ZHOOK_HANDLE)

      END SUBROUTINE CIWABR
