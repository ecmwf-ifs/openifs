! (C) Copyright 2015- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE VORDIV_TO_UV(PSPVOR,PSPDIV,PSPU,PSPV,KSMAX,KVSETUV)

!**** *VORDIV_TO_UV* - Convert spectral vorticity and divergence to spectral U (u*cos(theta)) and V (v*cos(theta).

!     Purpose.
!     --------
!        Interface routine for  Convert spectral vorticity and divergence to spectral U  and V 

!**   Interface.
!     ----------
!     CALL VORDIV_TO_UV(...)

!     Explicit arguments :
!     --------------------
!     PSPVOR(:,:) - spectral vorticity (input)
!     PSPDIV(:,:) - spectral divergence (input)
!     PSPU(:,:)   - spectral U (u*cos(theta) (output)
!     PSPV(:,:)   - spectral V (v*cos(theta) (output)
!     KSMAX       - spectral resolution (input)
!     KVSETUV(:)  - Optionally indicating which 'b-set' in spectral space owns a
!                   vor/div field. Equivalant to NBSETLEV in the IFS.
!                   The length of KVSETUV should be the GLOBAL number
!                   of u/v fields which is the dimension of u and v releated
!                   fields in grid-point space.

!     Method.
!     -------

!     Externals.  SET_RESOL   - set resolution
!     ----------  VD2UV_CTL   - control vordiv to uv

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 15-06-15


!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

!ifndef INTERFACE

USE TPM_GEN         ,ONLY : NERR, NOUT,MSETUP0
USE TPM_DISTR       ,ONLY : D, NPRTRV, MYSETV

USE SET_RESOL_MOD   ,ONLY : SET_RESOL
USE VD2UV_CTL_MOD ,ONLY   : VD2UV_CTL
USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

!endif INTERFACE

IMPLICIT NONE

! Declaration of arguments

REAL(KIND=JPRB), INTENT(IN) :: PSPVOR(:,:)
REAL(KIND=JPRB), INTENT(IN) :: PSPDIV(:,:)
REAL(KIND=JPRB), INTENT(OUT) :: PSPU(:,:)
REAL(KIND=JPRB), INTENT(OUT) :: PSPV(:,:)
INTEGER(KIND=JPIM) , INTENT(IN) :: KSMAX
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETUV(:)

!ifndef INTERFACE

! Local varaibles
INTEGER(KIND=JPIM) :: IUBOUND(4),J
INTEGER(KIND=JPIM) :: IF_UV,IF_UV_G,IRESOL,IDGL
LOGICAL :: LTMP_SETUP0
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "setup_trans0.h"
#include "setup_trans.h"
#include "trans_release.h"
#include "trans_end.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('VORDIV_TO_UV',0,ZHOOK_HANDLE)

!CALL GSTATS(XXXX,0)

IF(MSETUP0 == 0) THEN
  CALL SETUP_TRANS0()
  LTMP_SETUP0 = .TRUE.
ELSE
  LTMP_SETUP0 = .FALSE.
ENDIF
IDGL = 2 ! It doesn't matter as long as it's a positive even number
CALL SETUP_TRANS(KSMAX,IDGL,LDSPSETUPONLY=.TRUE.,KRESOL=IRESOL)
CALL SET_RESOL(IRESOL)


! Set defaults

IF_UV = 0
IF_UV_G = 0
! Decide requirements

IF(PRESENT(KVSETUV)) THEN
  IF_UV_G = UBOUND(KVSETUV,1)
  DO J=1,IF_UV_G
    IF(KVSETUV(J) > NPRTRV .OR. KVSETUV(J) < 1) THEN
      WRITE(NERR,*) 'VORDIV_TO_UV:KVSETUV(J) > NPRTRV ',J,KVSETUV(J),NPRTRV
      CALL ABORT_TRANS('VORDIV_TO_UV:KVSETUV TOO LONG OR CONTAINS VALUES OUTSIDE RANGE')
    ENDIF
    IF(KVSETUV(J) == MYSETV) THEN
      IF_UV = IF_UV+1
    ENDIF
  ENDDO
ELSE
  IF_UV = UBOUND(PSPVOR,1)
  IF_UV_G = IF_UV
ENDIF

! Consistency checks

IF (IF_UV > 0) THEN
  IF(UBOUND(PSPVOR,1) < IF_UV) THEN
    WRITE(NERR,*)'VORDIV_TO_UV : UBOUND(PSPVOR,1) < IF_UV ',UBOUND(PSPVOR,1),IF_UV
    CALL ABORT_TRANS('VORDIV_TO_UV : PSPVOR TOO SHORT')
  ENDIF
  IF(UBOUND(PSPDIV,1) < IF_UV) THEN
    WRITE(NERR,*)'VORDIV_TO_UV : UBOUND(PSPDIV,1) < IF_UV ',UBOUND(PSPDIV,1),IF_UV
    CALL ABORT_TRANS('VORDIV_TO_UV : PSPDIV TOO SHORT')
  ENDIF
  IF(UBOUND(PSPU,1) < IF_UV) THEN
    WRITE(NERR,*)'VORDIV_TO_UV : UBOUND(PSPU,1) < IF_UV ',UBOUND(PSPU,1),IF_UV
    CALL ABORT_TRANS('VORDIV_TO_UV : PSPU TOO SHORT')
  ENDIF
  IF(UBOUND(PSPV,1) < IF_UV) THEN
    WRITE(NERR,*)'VORDIV_TO_UV : UBOUND(PSPV,1) < IF_UV ',UBOUND(PSPV,1),IF_UV
    CALL ABORT_TRANS('VORDIV_TO_UV : PSPV TOO SHORT')
  ENDIF
ENDIF


IF(NPRTRV >1) THEN
  IF(IF_UV > 0 .AND. .NOT. PRESENT(KVSETUV)) THEN
    WRITE(NERR,*)'NPRTRV >1 AND IF_UV > 0 AND NOT PRESENT(KVSETUV)',&
                 &NPRTRV,IF_UV
    CALL ABORT_TRANS('VORDIV_TO_UV: SPECIFY VERTICAL SPECTRAL DISTRIBUTION!')
  ENDIF
ENDIF

!CALL GSTATS(XXXX,1)

!     ------------------------------------------------------------------

! Perform transform

CALL VD2UV_CTL(IF_UV,PSPVOR,PSPDIV,PSPU,PSPV)

CALL TRANS_RELEASE(IRESOL)
IF (LTMP_SETUP0) THEN
  CALL TRANS_END()
ENDIF

IF (LHOOK) CALL DR_HOOK('VORDIV_TO_UV',1,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

!endif INTERFACE

END SUBROUTINE VORDIV_TO_UV

