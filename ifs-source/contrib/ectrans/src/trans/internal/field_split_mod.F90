! (C) Copyright 2001- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE FIELD_SPLIT_MOD
CONTAINS
SUBROUTINE FIELD_SPLIT(KBLK,KF_GP,KKF_UV_G,KVSETUV,KVSETSC,&
 & KSTUV_G,KENUV_G,KF_UV_G,KSTSC_G,KENSC_G,KF_SCALARS_G,&
 & KSTUV,KENUV,KF_UV,KSTSC,KENSC,KF_SCALARS)

!**** *FIELD_SPLIT* - Split fields

!     Purpose.
!     --------
!        Split fields

!**   Interface.
!     ----------
!     CALL FIELD_SPLIT(...)

!     Explicit arguments :
!     --------------------
!     KBLK          - block number
!     KF_GP         - total number of output gridpoint fields
!     KKF_UV_G      - global number of spectral u-v fields
!     KVSETUV       - IVSETUV from SHUFFLE
!     KVSETSC       - IVSETUV from SHUFFLE

!     All the following output arguments are quantities for THIS packet.
!     KSTUV_G       -
!     KENUV_G       -
!     KF_UV_G       -
!     KSTSC_G       -
!     KENSC_G       -
!     KF_SCALARS_G  -
!     KSTUV         -
!     KENUV         -
!     KF_UV         -
!     KSTSC         -
!     KENSC         -
!     KF_SCALARS    -

!     Externals.  NONE
!     ----------

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 01-01-03

!     ------------------------------------------------------------------
USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE TPM_GEN         ,ONLY : NPROMATR
!USE TPM_TRANS
USE TPM_DISTR       ,ONLY : MYSETV, NPRTRV
!

IMPLICIT NONE

! Declaration of arguments

INTEGER(KIND=JPIM),INTENT(IN)  :: KBLK,KF_GP,KKF_UV_G
INTEGER(KIND=JPIM), INTENT(IN)  :: KVSETUV(:),KVSETSC(:)
INTEGER(KIND=JPIM), INTENT(OUT) :: KSTUV_G,KENUV_G,KF_UV_G,KSTSC_G,KENSC_G,KF_SCALARS_G
INTEGER(KIND=JPIM), INTENT(OUT) :: KSTUV,KENUV,KF_UV,KSTSC,KENSC,KF_SCALARS

! Local variables

INTEGER(KIND=JPIM) :: ISTF,IENF,J

!     ------------------------------------------------------------------

ISTF  = (KBLK-1)*NPROMATR+1
IENF  = MIN(KBLK*NPROMATR,KF_GP)

KSTUV_G = (KBLK-1)*NPROMATR/2+1
KENUV_G = MIN(KBLK*NPROMATR/2,KKF_UV_G)
IF(ISTF > 2*KKF_UV_G) KSTUV_G = KENUV_G+1
KF_UV_G = KENUV_G-KSTUV_G+1
KSTSC_G = MAX(ISTF-2*KKF_UV_G,1)
KENSC_G = MAX(IENF-2*KKF_UV_G,0)
KF_SCALARS_G = KENSC_G-KSTSC_G+1

! Spectral fields distributed over fields

IF(NPRTRV > 1) THEN
  KF_UV = 0
  KSTUV = 1
  DO J=1,KSTUV_G-1
    IF(KVSETUV(J) == MYSETV) THEN
      KSTUV = KSTUV+1
    ENDIF
  ENDDO
  KENUV = KSTUV-1
  DO J=KSTUV_G,KENUV_G
    IF(KVSETUV(J) == MYSETV) THEN
      KF_UV = KF_UV+1
      KENUV = KENUV+1
    ENDIF
  ENDDO
  KF_SCALARS = 0
  KSTSC = 1
  DO J=1,KSTSC_G-1
    IF(KVSETSC(J) == MYSETV) THEN
      KSTSC =KSTSC+1
    ENDIF
  ENDDO
  KENSC = KSTSC-1
  DO J=KSTSC_G,KENSC_G
    IF(KVSETSC(J) == MYSETV) THEN
      KF_SCALARS = KF_SCALARS+1
      KENSC = KENSC+1
    ENDIF
  ENDDO
ELSE

  ! Spectral fields not distributed over fields

  KF_UV = KF_UV_G
  KSTUV = KSTUV_G
  KENUV = KENUV_G
  KF_SCALARS = KF_SCALARS_G
  KSTSC = KSTSC_G
  KENSC = KENSC_G
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE FIELD_SPLIT
END MODULE FIELD_SPLIT_MOD
