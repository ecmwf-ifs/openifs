! (C) Copyright 2005- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE DR_HOOK_UTIL_MULTI(LDHOOK,CDNAME,KCASE,PKEY,KPKEY,CDFILENAME,KSIZEINFO)
USE EC_PARKIND  ,ONLY : JPIM     ,JPRD
USE OML_MOD,ONLY : OML_GET_MAX_THREADS,OML_MY_THREAD
IMPLICIT NONE
LOGICAL,INTENT(INOUT)       :: LDHOOK
CHARACTER(LEN=*),INTENT(IN) :: CDNAME,CDFILENAME
INTEGER(KIND=JPIM),INTENT(IN) :: KPKEY, KCASE,KSIZEINFO
REAL(KIND=JPRD),INTENT(INOUT) :: PKEY(KPKEY)

LOGICAL,SAVE :: LL_FIRST_TIME = .TRUE.
REAL(KIND=JPRD) :: ZDUMMY
INTEGER(KIND=JPIM) :: IMYTID, ISILENT, IMAXTH

#include "dr_hook_util.h"

! -----------------------------------------------------------------

IF (.NOT.LDHOOK) RETURN
IF (LL_FIRST_TIME) THEN
  LL_FIRST_TIME = .FALSE.
  CALL DR_HOOK_UTIL(LDHOOK,'',-1,ZDUMMY,'',-1_JPIM)

  ! Approximately the very first OpenMP-loop
  IMAXTH = OML_GET_MAX_THREADS()
  ! trapfpe setting also for slave threads -- was missing
  !$OMP PARALLEL DO SCHEDULE(STATIC,1) PRIVATE(IMYTID,ISILENT) IF (IMAXTH > 1)
  DO IMYTID=1,IMAXTH
     ISILENT = 1 ! no verbosity
     IF (IMYTID == IMAXTH) ISILENT = 0 ! be verbose with the last thread
     CALL TRAPFPE_SLAVE_THREADS(ISILENT) ! see drhook.c; does not do anything on master thread
  ENDDO ! IMYTID=1,IMAXTH
  !$OMP END PARALLEL DO
ENDIF

!$OMP PARALLEL DO SCHEDULE(STATIC,1) PRIVATE(IMYTID)
DO IMYTID=1,KPKEY
  IF (KCASE == 0) THEN
    CALL C_DRHOOK_START(CDNAME, IMYTID, PKEY(IMYTID), CDFILENAME, KSIZEINFO)
  ELSE IF (KCASE == 1) THEN
    CALL C_DRHOOK_END  (CDNAME, IMYTID, PKEY(IMYTID), CDFILENAME, KSIZEINFO)
  ENDIF
ENDDO ! IMYTID=1,KPKEY
!$OMP END PARALLEL DO

END SUBROUTINE DR_HOOK_UTIL_MULTI
