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

SUBROUTINE ININEMOIO(LDNEMOIO,LDNEMOIOSERVER,KREQUIRED,KPROVIDED,LDMPI1)
!
!**** *ININEMOIO*  - Initialize NEMO IO server (if present)
!
!     Purpose.
!     --------
!     Get a IFS MPI communicator to use when NEMO IO server is active.
!
!**   Interface.
!     ----------
!       *CALL*  *ININEMOIO*
!
!     Input:
!     -----
!
!     Output:
!     ------
!       MPLUSERCOMM in MPL_MODULE is updated.
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
!     -----------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM 
USE MPL_MODULE
IMPLICIT NONE

LOGICAL, INTENT(OUT) :: LDNEMOIO, LDNEMOIOSERVER
INTEGER, INTENT(INOUT) :: KREQUIRED,KPROVIDED
LOGICAL, INTENT(IN) :: LDMPI1
CHARACTER(LEN=128) :: CLNEMOIO              ! String to store NEMOIO env.
INTEGER(KIND=JPIM) :: ILOCAL_COMM           ! Local communicator 
#ifdef WITH_NEMO

! Get the enviroment variables

!IF (LHOOK) CALL DR_HOOK('ININEMOIO',0,ZHOOK_HANDLE)

CALL GET_ENVIRONMENT_VARIABLE('NEMOIOSERVER',CLNEMOIO)

! Set NEMOIO settings depending of the NEMOIO environment variables
SELECT CASE (TRIM(CLNEMOIO))
! If NEMOIO is not set don't use NEMOIO.
CASE ('')
   LDNEMOIO    =.FALSE.
! If NEMOIO is "no" don't use NEMOIO.
CASE ('no')
   LDNEMOIO    =.FALSE.
! If NEMOIO is "yes" use NEMOIO for real.
CASE ('yes') 
   LDNEMOIO    =.TRUE.
CASE DEFAULT
   CALL ABOR1('Invalid value of NEMOIO environment variable in ininemoio.F90')
END SELECT

IF (LDNEMOIO) THEN
  CALL NEMOGCMCOUP_INIT_IOSERVER( ILOCAL_COMM, LDNEMOIOSERVER, &
    & KREQUIRED, KPROVIDED, LDMPI1 )
  LMPLUSERCOMM=.TRUE.
  MPLUSERCOMM=ILOCAL_COMM
ELSE
  LDNEMOIOSERVER=.FALSE.
ENDIF

#else

LDNEMOIO=.FALSE.
LDNEMOIOSERVER=.FALSE.
!IF (LHOOK) CALL DR_HOOK('ININEMOIO',1,ZHOOK_HANDLE)

#endif

END SUBROUTINE ININEMOIO
