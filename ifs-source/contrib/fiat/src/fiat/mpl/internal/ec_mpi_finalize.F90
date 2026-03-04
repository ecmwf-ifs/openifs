! (C) Copyright 2005- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE EC_MPI_FINALIZE(KERROR,LDCALLFINITO,LDMEMINFO,CALLER)
USE EC_PARKIND, ONLY : JPIM
USE MPL_MPIF
IMPLICIT NONE
INTEGER(KIND=JPIM), INTENT(OUT) :: KERROR
LOGICAL, INTENT(IN) :: LDCALLFINITO
LOGICAL, INTENT(IN) :: LDMEMINFO
CHARACTER(LEN=*), INTENT(IN) :: CALLER
LOGICAL :: LLINIT, LLFIN, LLNOTMPIWORLD
INTEGER(KIND=JPIM) :: IERR, ICOMM
INTEGER(KIND=JPIM) :: NCOMM_MEMINFO
COMMON /cmn_meminfo/ NCOMM_MEMINFO
#include "ec_meminfo.intfb.h"
#include "dr_hook_end.intfb.h"
KERROR = 0
IF (LDCALLFINITO) THEN !*** common MPI_Finalize()
  CALL MPI_INITIALIZED(LLINIT,IERR)
  IF (LLINIT .AND. IERR == 0) THEN
    CALL MPI_FINALIZED(LLFIN,IERR)
    IF (.NOT.LLFIN .AND. IERR == 0) THEN
      LLNOTMPIWORLD = (NCOMM_MEMINFO /= 0 .and. NCOMM_MEMINFO /= MPI_COMM_WORLD)
      IF (LLNOTMPIWORLD) THEN
        ICOMM = NCOMM_MEMINFO
      ELSE
        ICOMM = MPI_COMM_WORLD
      ENDIF

      IF( LDMEMINFO ) CALL EC_MEMINFO(-1,"ec_mpi_finalize:"//caller,ICOMM,KBARR=1,KIOTASK=-1,KCALL=1)

      CALL DR_HOOK_END() ! Make sure DrHook output is produced before MPI_Finalize (in case it fails)
      CALL MPI_BARRIER(ICOMM,IERR)
      IF (LLNOTMPIWORLD) THEN
        ! CALL MPI_COMM_FREE(NCOMM_MEMINFO,IERR)
        NCOMM_MEMINFO = 0
      ENDIF
      CALL MPI_FINALIZE(KERROR)
    ENDIF
  ENDIF
ENDIF
END SUBROUTINE EC_MPI_FINALIZE
