! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE WAM_MULTIO_MOD
!
!== wam_multio_mod.F90 ==
!
! A module to encapsulate MULTIO-calls in order to call Dr.Hook
!
! Author: Willem Deconinck, 9-Aug-2021
!

USE EC_PARKIND,  ONLY : JPIM
USE YOMHOOK,     ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOWABORT,    ONLY : WAM_ABORT

IMPLICIT NONE
PRIVATE

PUBLIC :: WAM_MULTIO_FLUSH
PUBLIC :: WAM_MULTIO_WRITE

CONTAINS

! ------------------------------------------------------------------

SUBROUTINE WAM_MULTIO_FLUSH(KRET)
INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT)   :: KRET

#ifdef WAM_HAVE_MULTIO
INTEGER(KIND=JPIM), EXTERNAL :: IMULTIO_FLUSH
#endif

INTEGER(KIND=JPIM) :: IRET
CHARACTER(LEN=256) :: CDMSG
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('WAM_MULTIO_FLUSH',0,ZHOOK_HANDLE)

#ifdef WAM_HAVE_MULTIO
IRET = IMULTIO_FLUSH()
#else
IRET = -1
#endif
IF(PRESENT(KRET)) THEN
  KRET = IRET
ELSEIF(IRET /= 0) THEN
#ifdef WAM_HAVE_MULTIO
    WRITE(CDMSG,'(A,I0)') 'IMULTIO_FLUSH FAILED. Return code: ',IRET
    CALL WAM_ABORT('WAM_MULTIO_FLUSH: '//TRIM(CDMSG))
#else
    CALL WAM_ABORT('WAM_MULTIO_FLUSH: wam is not compiled with multio support')
#endif
ENDIF

IF (LHOOK) CALL DR_HOOK('WAM_MULTIO_FLUSH',1,ZHOOK_HANDLE)
END SUBROUTINE WAM_MULTIO_FLUSH

! ------------------------------------------------------------------

SUBROUTINE WAM_MULTIO_WRITE(KGRIB,KLEN,KRET)
INTEGER(KIND=JPIM),           INTENT(IN)    :: KLEN, KGRIB(:)
INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT)   :: KRET

#ifdef WAM_HAVE_MULTIO
INTEGER(KIND=JPIM), EXTERNAL :: IMULTIO_WRITE
#endif

INTEGER(KIND=JPIM) :: IRET
CHARACTER(LEN=256) :: CDMSG
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('WAM_MULTIO_WRITE',0,ZHOOK_HANDLE)

#ifdef WAM_HAVE_MULTIO
IRET = IMULTIO_WRITE(KGRIB(:), KLEN)
#else
IRET = -1
#endif
IF(PRESENT(KRET)) THEN
  KRET = IRET
ELSEIF(IRET /= 0) THEN
#ifdef WAM_HAVE_MULTIO
  WRITE(CDMSG,'(A,I0,A,I0)') 'IMULTIO_WRITE(KLEN=', KLEN,') FAILED. Return code: ',IRET
  CALL WAM_ABORT('WAM_MULTIO_WRITE: '//TRIM(CDMSG))
#else
  CALL WAM_ABORT('WAM_MULTIO_WRITE: wam is not compiled with multio support')
#endif
ENDIF

IF (LHOOK) CALL DR_HOOK('WAM_MULTIO_WRITE',1,ZHOOK_HANDLE)
END SUBROUTINE WAM_MULTIO_WRITE

! ------------------------------------------------------------------

END MODULE WAM_MULTIO_MOD
