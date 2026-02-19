! (C) Copyright 2005- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

#define DR_HOOK_ASSERT_MPI_INITITALIZED 1
  !! TEMPORARY !!
  !! DR_HOOK will abort when MPI is not initialized.
  !! DR_HOOK used to initialize MPI via MPL_INIT, but no longer.

SUBROUTINE DR_HOOK_INIT()
  !! Initialises DR_HOOK
  !! Also calls
  !! - OML_INIT to e.g. save the OMP_NUM_THREADS environment variable
  !! - EC_ARGS to save command-line arguments in case the main program is Fortran
  !! Environment variable "DR_HOOK=1" will enable DR_HOOK
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT
  USE OML_MOD       ,ONLY : OML_MAX_THREADS, OML_INIT
  USE EC_ARGS_MOD   ,ONLY : EC_ARGS
  USE YOMHOOK       ,ONLY : LHOOK ! True by default
  IMPLICIT NONE
  LOGICAL,SAVE :: LL_FIRST_TIME = .TRUE.
  CHARACTER(LEN=512) :: CLENV
  INTEGER(KIND=C_INT) :: IMAX_THREADS

  IF (LL_FIRST_TIME) THEN
    LL_FIRST_TIME = .FALSE.

    CALL GET_ENVIRONMENT_VARIABLE('DR_HOOK',CLENV)
    IF ( CLENV == ' ' .OR. CLENV == '0'  .OR. &
       & CLENV == 'false' .OR. CLENV == 'FALSE' .OR. &
       & CLENV == 'off'   .OR. CLENV == 'OFF' .OR. &
       & CLENV == 'no'    .OR. CLENV == 'NO' ) THEN
       LHOOK = .FALSE.
    ENDIF
    CALL OML_INIT()

#if DR_HOOK_ASSERT_MPI_INITITALIZED
    CALL DR_HOOK_ASSERT_MPI_INITIALIZED_()
#endif

    CALL EC_ARGS()

    IF (.NOT. LHOOK) RETURN
    
    IMAX_THREADS = OML_MAX_THREADS()
    CALL C_DRHOOK_INIT('',IMAX_THREADS)
      !! First argument (progname) is empty ==> c_drhook_init will retrieve progname via ec_args itself
  ENDIF

CONTAINS

#if DR_HOOK_ASSERT_MPI_INITITALIZED
SUBROUTINE DR_HOOK_ASSERT_MPI_INITIALIZED_()
  LOGICAL :: LMPI_REQUIRED
  INTEGER :: ILEN
  INTEGER(KIND=C_INT) :: IERR
  LOGICAL :: LMPI_INITIALIZED
  INTEGER, PARAMETER :: NVARS = 5
  CHARACTER(LEN=32), DIMENSION(NVARS) :: CMPIRUN_DETECT
  CHARACTER(LEN=4) :: CLENV_DR_HOOK_ASSERT_MPI_INITIALIZED
  INTEGER :: IVAR

#include "mpif.h"
#include "abor1.intfb.h"

  ! Environment variables that are set when mpirun, srun, aprun, ... are used (see eckit/mpi/Comm.cc)
  CMPIRUN_DETECT(1) = 'OMPI_COMM_WORLD_SIZE'  ! OpenMPI
  CMPIRUN_DETECT(2) = 'ALPS_APP_PE'           ! Cray PE
  CMPIRUN_DETECT(3) = 'PMI_SIZE'              ! Intel
  CMPIRUN_DETECT(4) = 'SLURM_NTASKS'          ! Slurm

  LMPI_REQUIRED = .FALSE.
  DO IVAR=1,NVARS
    CALL GET_ENVIRONMENT_VARIABLE(NAME=TRIM(CMPIRUN_DETECT(IVAR)),LENGTH=ILEN)
    IF (ILEN > 0) THEN
      LMPI_REQUIRED = .TRUE.
      EXIT ! break
    ENDIF
  ENDDO

  IF (LMPI_REQUIRED) THEN
    CALL GET_ENVIRONMENT_VARIABLE(NAME="DR_HOOK_ASSERT_MPI_INITIALIZED", VALUE=CLENV)
    IF ( CLENV == '0'  .OR. &
       & CLENV == 'false' .OR. CLENV == 'FALSE' .OR. &
       & CLENV == 'off'   .OR. CLENV == 'OFF' .OR. &
       & CLENV == 'no'    .OR. CLENV == 'NO' ) THEN
      LMPI_REQUIRED = .FALSE.
    ENDIF
  ENDIF
  IF (LMPI_REQUIRED) THEN
    CALL MPI_INITIALIZED(LMPI_INITIALIZED,IERR)
    IF( IERR /= 0 ) THEN
      CALL ABOR1FL( "dr_hook_init.F90", __LINE__, &
        & "DR_HOOK: MPI_INITIALIZED failed" )
    ENDIF
    IF( .NOT.LMPI_INITIALIZED ) THEN
      CALL ABOR1FL( "dr_hook_init.F90", __LINE__, &
        & "DR_HOOK no longer calls MPL_INIT. Please initialize MPI (or MPL) before first DR_HOOK call."//NEW_LINE('A')//&
        & "This assertion can be disabled with environment: DR_HOOK_ASSERT_MPI_INITIALIZED=0" )
    ENDIF
  ENDIF
END SUBROUTINE
#endif
END SUBROUTINE
