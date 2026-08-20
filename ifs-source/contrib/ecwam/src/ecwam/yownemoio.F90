! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE YOWNEMOIO
IMPLICIT NONE
PUBLIC :: LNEMOIOSERVER
PUBLIC :: WAMININEMOIO
PUBLIC :: WAMENDNEMOIO

LOGICAL :: LNEMOIOSERVER = .FALSE.

CONTAINS

SUBROUTINE WAMININEMOIO(LNEMOIO)
!
!**** *WAMININEMOIO*  - Initialize NEMO IO server (if present)
!
!     Purpose.
!     --------
!     Get a IFS MPI communicator to use when NEMO IO server is active.
!
!**   Interface.
!     ----------
!       *CALL*  *WAMININEMOIO*
!
!     Input:
!     -----
!       LNEMOIO (optional, default=FALSE)     to enable NEMO IO SERVER
!           NOTE that NEMOIOSERVER environment variable can override LNEMOIO
!
!     Output:
!     ------
!       MPLUSERCOMM in MPL_MODULE is updated.
!       LNEMOIOSERVER module variable is updated
!
!     Method:
!     ------
!       NEMOIO usage is controlled by environment variables
!
!     Externals:
!     ---------
!       GETENV - Get enviroment variables
!
!     Reference:
!     ---------
!
!     Author:
!     -------
!       K. Mogensen, ECMWF
!
!     Modifications.
!     --------------
!
! -----------------------------------------------------------

USE PARKIND_WAVE, ONLY : JWIM
USE MPL_MODULE,   ONLY : LMPLUSERCOMM, MPLUSERCOMM, MPI_COMM_WORLD

! -----------------------------------------------------------

IMPLICIT NONE

LOGICAL, INTENT(IN), OPTIONAL :: LNEMOIO
LOGICAL :: LLNEMOIO
INTEGER(KIND=JWIM) :: ILOCAL_COMM, ICOMM    ! Local communicators

INTEGER, EXTERNAL :: EC_MPIRANK

LNEMOIOSERVER=.FALSE.

LLNEMOIO = .FALSE.
IF (PRESENT(LNEMOIO)) LLNEMOIO = LNEMOIO

CALL UPDATE_FROM_ENVIRONMENT_VARIABLE('NEMOIOSERVER',LLNEMOIO)

#ifdef WITH_NEMO

IF (LLNEMOIO) THEN
  CALL NEMOGCMCOUP_INIT_IOSERVER( ILOCAL_COMM, LNEMOIOSERVER )
  ICOMM=MPI_COMM_WORLD
  CALL NEMOGCMCOUP_INIT_IOSERVER_2( ICOMM )
  LMPLUSERCOMM=.TRUE.
  MPLUSERCOMM=ILOCAL_COMM
ENDIF

#endif

IF (EC_MPIRANK() == 0) THEN ! Access MPI rank before MPI_INIT is called
  WRITE(6,'(A,L)') "*** WAMININEMOIO: LNEMOIOSERVER = ",LNEMOIOSERVER
ENDIF

CONTAINS

      SUBROUTINE UPDATE_FROM_ENVIRONMENT_VARIABLE(CNAME,LLNEMOIO)
            USE YOWABORT, ONLY : WAM_ABORT
            IMPLICIT NONE

            CHARACTER(LEN=*), INTENT(IN)    :: CNAME
            LOGICAL,          INTENT(INOUT) :: LLNEMOIO

            CHARACTER(LEN=128) :: CLNEMOIO
            INTEGER(KIND=JWIM) :: ISTATUS

            INTEGER, EXTERNAL :: EC_MPIRANK

            CALL GET_ENVIRONMENT_VARIABLE(NAME=TRIM(CNAME),VALUE=CLNEMOIO,STATUS=ISTATUS)

            IF (ISTATUS == 0) THEN ! Variable exists
              SELECT CASE (TRIM(CLNEMOIO))
                CASE ('no')
                CASE ('OFF')
                CASE ('0')
                CASE ('F')
                  LLNEMOIO    =.FALSE.
                CASE ('yes')
                CASE ('ON')
                CASE ('1')
                CASE ('T')
                  LLNEMOIO    =.TRUE.
                CASE DEFAULT
                  CALL WAM_ABORT('Invalid value of NEMOIOSERVER environment variable.', &
           &              __FILENAME__, __LINE__)
              END SELECT
              IF (EC_MPIRANK() == 0) THEN ! Access MPI rank before MPI_INIT is called
                  WRITE(6,'(5A,L)') "*** WAMININEMOIO: ${",TRIM(CNAME),"}=",TRIM(CLNEMOIO)
              ENDIF
            ENDIF

      END SUBROUTINE UPDATE_FROM_ENVIRONMENT_VARIABLE

END SUBROUTINE WAMININEMOIO

! ----------------------------------------------------------------------------------------------------------------------

SUBROUTINE WAMENDNEMOIO
!
!**** *WAMENDNEMOIO*  - End
!
!     Purpose.
!     --------
!     Call MPI finalize in case of IO servers.
!
!**   Interface.
!     ----------
!       *CALL*  *WAMENDNEMOIO*
!
!     Input:
!     -----
!
!     Output:
!     ------
!
!     Method:
!     ------
!       NEMOIO usage is controlled by environment variables
!
!     Externals:
!     ---------
!
!     Reference:
!     ---------
!
!     Author:
!     -------
!       K. Mogensen, ECMWF
!
!     Modifications.
!     --------------
!
!     -----------------------------------------------------------

IMPLICIT NONE

#ifdef WITH_NEMO

IF (LNEMOIOSERVER) THEN
    CALL NEMOGCMCOUP_END_IOSERVER
ENDIF

#endif

END SUBROUTINE WAMENDNEMOIO

END MODULE YOWNEMOIO
