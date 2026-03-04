! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction
! 
! (C) Copyright 1989- Meteo-France.
! 

SUBROUTINE SUIOSTREAM
USE PARKIND1,            ONLY : JPIM, JPRB
USE ALGORITHM_STATE_MOD, ONLY : SET_IOSTREAM_FDB
USE YOMHOOK,             ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCT0,              ONLY : LFDBOP
USE YOMMP0,              ONLY : NPROC, NWRTOUT
USE YOMVAR,              ONLY : LFDBERR
USE IOSTREAM_MIX,        ONLY : SETUP_IOSTREAM, Y_IOSTREAM_FDB
IMPLICIT NONE

INTEGER(KIND=JPIM) :: IOP,JROC
INTEGER(KIND=JPIM) :: IOPROCS(NPROC)
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUIOSTREAM',0,ZHOOK_HANDLE)

IF (LFDBOP .OR. LFDBERR) THEN
  IOP=0
  DO JROC=1,NPROC,NWRTOUT
    IOP = IOP+1
    IOPROCS(IOP) = JROC
  ENDDO
  CALL SETUP_IOSTREAM(Y_IOSTREAM_FDB,CDTYPE='FDB',CDMODE='w',&
   & KIOPROCS=IOPROCS(1:IOP))
  CALL SET_IOSTREAM_FDB(.TRUE.)
ENDIF

IF (LHOOK) CALL DR_HOOK('SUIOSTREAM',1,ZHOOK_HANDLE)
END SUBROUTINE SUIOSTREAM
