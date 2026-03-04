! (C) Copyright 2001- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE SHUFFLE_MOD
CONTAINS
SUBROUTINE SHUFFLE(KF_UV_G,KF_SCALARS_G,KSHFUV_G,KIVSETUV,KSHFSC_G,KIVSETSC,&
 & KVSETUV,KVSETSC)

!**** *SHUFFLE* - Re-shuffle fields for load balancing

!     Purpose.
!     --------
!        Re-shuffle fields for load balancing if NPRTRV>1. Note that the
!      relative order of the local spectral fields has to maintained.

!**   Interface.
!     ----------
!     CALL SHUFFLE(...)

!     Explicit arguments :
!     --------------------
!     KF_UV_G      - global number of spectral u-v fields
!     KF_SCALARS_G - global number of scalar spectral fields
!     KSHFUV_G     - reshuffling index for uv fields
!     KIVSETUV     - reshuffled KVSETUV
!     KSHFSC_G     - reshuffling index for scalar fields
!     KIVSETSC     - reshuffled KVSETSC
!     KVSETUV(:)  - indicating which 'b-set' in spectral space owns a
!                   vor/div field. Equivalant to NBSETLEV in the IFS.
!                   The length of KVSETUV should be the GLOBAL number
!                   of u/v fields which is the dimension of u and v releated
!                   fields in grid-point space.
!     KVESETSC(:) - indicating which 'b-set' in spectral space owns a
!                   scalar field. As for KVSETUV this argument is required
!                   if the total number of processors is greater than
!                   the number of processors used for distribution in
!                   spectral wave space.

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

!USE TPM_GEN
!USE TPM_TRANS
USE TPM_DISTR       ,ONLY : NPRTRV
!

IMPLICIT NONE

! Declaration of arguments

INTEGER(KIND=JPIM), INTENT(IN)  :: KF_UV_G,KF_SCALARS_G
INTEGER(KIND=JPIM), INTENT(OUT) :: KSHFUV_G(:),KSHFSC_G(:)
INTEGER(KIND=JPIM), INTENT(OUT) :: KIVSETUV(:),KIVSETSC(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETUV(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETSC(:)

INTEGER(KIND=JPIM) :: IHELP(MAX(KF_UV_G,KF_SCALARS_G),NPRTRV),IHELPC(NPRTRV)
INTEGER(KIND=JPIM) :: IDW,J

!     ------------------------------------------------------------------

IF(NPRTRV > 1) THEN
  IHELP(:,:) = 0
  IHELPC(:)  = 0
  DO J=1,KF_UV_G
    IHELPC(KVSETUV(J)) = IHELPC(KVSETUV(J))+1
    IHELP(IHELPC(KVSETUV(J)),KVSETUV(J)) = J
  ENDDO
  IDW = KF_UV_G+1
  DO
    DO J=NPRTRV,1,-1
      IF(IHELPC(J) > 0) THEN
        IDW = IDW-1
        KSHFUV_G(IDW) = IHELP(IHELPC(J),J)
        IHELPC(J) =IHELPC(J)-1
      ENDIF
    ENDDO
    IF(IDW == 1) EXIT
  ENDDO

  IHELP(:,:) = 0
  IHELPC(:)  = 0
  DO J=1,KF_SCALARS_G
    IHELPC(KVSETSC(J)) = IHELPC(KVSETSC(J))+1
    IHELP(IHELPC(KVSETSC(J)),KVSETSC(J)) = J
  ENDDO
  IDW = KF_SCALARS_G+1
  DO
    DO J=NPRTRV,1,-1
      IF(IHELPC(J) > 0) THEN
        IDW = IDW-1
        KSHFSC_G(IDW) = IHELP(IHELPC(J),J)
        IHELPC(J) =IHELPC(J)-1
      ENDIF
    ENDDO
    IF(IDW == 1) EXIT
  ENDDO

  DO J=1,KF_UV_G
    KIVSETUV(J) = KVSETUV(KSHFUV_G(J))
  ENDDO
  DO J=1,KF_SCALARS_G
    KIVSETSC(J) = KVSETSC(KSHFSC_G(J))
  ENDDO
ELSE
  DO J=1,KF_UV_G
    KSHFUV_G(J) = J
    KIVSETUV(J) = 1
  ENDDO
  DO J=1,KF_SCALARS_G
    KSHFSC_G(J) = J
    KIVSETSC(J) = 1
  ENDDO
ENDIF

!     ------------------------------------------------------------------
  
END SUBROUTINE SHUFFLE
END MODULE SHUFFLE_MOD
