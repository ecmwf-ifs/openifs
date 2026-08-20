! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE CIMSSTRN(KIJS, KIJL, FL1, WAVNUM, DEPTH, CITHICK, STRN)

! ----------------------------------------------------------------------

!**** *CIMSSTRN* - COMPUTATION OF THE MEAN SQUARE WAVE STRAIN IN SEA ICE.

!     J. BIDLOT  ECMWF  JANUARY 2013.

!*    PURPOSE.
!     --------

!       COMPUTES MEAN SQUARE WAVE STRAIN AT EACH GRID POINT.

!**   INTERFACE.
!     ----------

!       *CALL* *CIMSSTRN (KIJS, KIJL, FL1, WAVNUM, DEPTH, CITHICK, STRN)*
!              *KIJS*    - INDEX OF FIRST GRIDPOINT
!              *KIJL*    - INDEX OF LAST GRIDPOINT
!              *FL1*     - SPECTRUM.
!              *WAVNUM*  - OPEN WATER WAVE NUMBER
!              *DEPTH*   - WATER DEPTH
!              *CITHICK* - SEA ICE THICKNESS
!              *STRN*    - MEAN SQUARE WAVE STRAIN IN ICE (OUTPUT).

!     METHOD.
!     -------

!      !!! IT ASSUMES SO DEFAULT SETTING FOR THE MECHANICAL PROPERTIES OF
!          THE SEA ICE (SEE AKI_ICE) !!!!!!!

!       NONE.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWICE   , ONLY : FLMIN
      USE YOWFRED  , ONLY : FR       ,DFIM     ,DELTH
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : G        ,ZPI      ,ROWATER

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "aki_ice.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(IN) :: FL1
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE), INTENT(IN) :: WAVNUM
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: DEPTH 
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: CITHICK 
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(OUT) :: STRN


      INTEGER(KIND=JWIM) :: IJ, M, K
      REAL(KIND=JWRB) :: F1LIM 
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: XKI, E, SUME 
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('CIMSSTRN',0,ZHOOK_HANDLE)

!*    1. INITIALISE
!        ----------

      F1LIM=FLMIN/DELTH

      DO IJ=KIJS,KIJL
        STRN(IJ) = 0.0_JWRB
      ENDDO

! ----------------------------------------------------------------------

!*    2. INTEGRATE OVER FREQUENCIES AND DIRECTIONS.
!        ------------------------------------------

      DO M=1,NFRE
        DO IJ=KIJS,KIJL
          XKI(IJ)=AKI_ICE(G,WAVNUM(IJ,M),DEPTH(IJ),ROWATER,CITHICK(IJ))
          E(IJ)=0.5_JWRB*CITHICK(IJ)*XKI(IJ)**3/WAVNUM(IJ,M)
        ENDDO

        DO IJ=KIJS,KIJL
          SUME(IJ) = 0.0_JWRB
        ENDDO
        DO K=1,NANG
          DO IJ=KIJS,KIJL
            SUME(IJ) = SUME(IJ)+FL1(IJ,K,M)
          ENDDO
        ENDDO

        DO IJ=KIJS,KIJL
          IF (SUME(IJ) > F1LIM) THEN
            STRN(IJ) = STRN(IJ)+E(IJ)**2*SUME(IJ)*DFIM(M)
          ENDIF
        ENDDO

      ENDDO

      IF (LHOOK) CALL DR_HOOK('CIMSSTRN',1,ZHOOK_HANDLE)

      END SUBROUTINE CIMSSTRN
