! (C) Copyright 2005- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE DR_HOOK_PROCINFO(KMYPROC, KNPROC)
USE EC_PARKIND  ,ONLY : JPIM
USE MPL_MPIF
IMPLICIT NONE
INTEGER(KIND=JPIM),INTENT(OUT) :: KMYPROC, KNPROC
LOGICAL :: LMPI_INITIALIZED
INTEGER(KIND=JPIM) :: IERROR

KMYPROC=-1
KNPROC=0
CALL MPI_INITIALIZED(LMPI_INITIALIZED,IERROR)
IF( LMPI_INITIALIZED ) THEN
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,KNPROC,IERROR)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,KMYPROC,IERROR)
  KMYPROC = KMYPROC+1 ! 1-based in IFS context
ENDIF
END SUBROUTINE DR_HOOK_PROCINFO
