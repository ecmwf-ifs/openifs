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

SUBROUTINE ENDNEMOIO
!
!**** *ENDNEMOIO*  - End
!
!     Purpose.
!     --------
!     Call MPI finalize in case of IO servers.
!
!**   Interface.
!     ----------
!       *CALL*  *ENDNEMOIO*
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

CHARACTER(LEN=128) :: CLNEMOIO              ! String to store NEMOIO env.
LOGICAL :: LLNEMOIO 

#ifdef WITH_NEMO

!IF (LHOOK) CALL DR_HOOK('ENDNEMOIO',0,ZHOOK_HANDLE)

CALL GET_ENVIRONMENT_VARIABLE('NEMOIOSERVER',CLNEMOIO)

! Set NEMOIO settings depending of the NEMOIO environment variables
SELECT CASE (TRIM(CLNEMOIO))
! If NEMOIO is not set don't use NEMOIO.
CASE ('')
   LLNEMOIO    =.FALSE.
! If NEMOIO is "no" don't use NEMOIO.
CASE ('no')
   LLNEMOIO    =.FALSE.
! If NEMOIO is "yes" use NEMOIO for real.
CASE ('yes') 
   LLNEMOIO    =.TRUE.
CASE DEFAULT
   CALL ABOR1('Invalid value of NEMOIO environment variable in ininemoio.F90')
END SELECT

IF(LLNEMOIO) CALL NEMOGCMCOUP_END_IOSERVER

!IF (LHOOK) CALL DR_HOOK('ENDNEMOIO',1,ZHOOK_HANDLE)

#endif

END SUBROUTINE ENDNEMOIO
