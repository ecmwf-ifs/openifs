! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE ALPHAP_TAIL(KIJS, KIJL, FL1, ALPHAP)

! ----------------------------------------------------------------------

!**** *ALPHAP_TAIL* - COMPUTATION OF PHILLIPS PARAMETER using the last frequency bins


!**   INTERFACE.
!     ----------

!       *CALL* *ALPHAP_TAIL(KIJS, KIJL, FL1, ALPHAP)
!          *KIJS*    - INDEX OF FIRST GRIDPOINT
!          *KIJL*    - INDEX OF LAST GRIDPOINT
!          *FL1*     - SPECTRA
!          *ALPHAP*  - PHILLIPS PARAMETER 

!     METHOD.
!     -------

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : TH       ,FR5      ,DELTH
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : G        ,ZPI      ,ZPI4GM2

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(IN) :: FL1
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(OUT) :: ALPHAP

      INTEGER(KIND=JWIM) :: IJ, K

      REAL(KIND=JWRB) :: CONST
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('ALPHAP_TAIL',0,ZHOOK_HANDLE)

!     COMPUTE THE PHILLIPS PARAMETER
      CONST = DELTH*ZPI4GM2*FR5(NFRE)
      ALPHAP(:) = 0.0_JWRB
      DO K = 1, NANG
        DO IJ = KIJS, KIJL
          ALPHAP(IJ) = ALPHAP(IJ) + CONST*FL1(IJ,K,NFRE)
        ENDDO
      ENDDO

IF (LHOOK) CALL DR_HOOK('ALPHAP_TAIL',1,ZHOOK_HANDLE)

END SUBROUTINE ALPHAP_TAIL
