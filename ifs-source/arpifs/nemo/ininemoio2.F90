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

SUBROUTINE ININEMOIO2(KNEMOCOMM)
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

INTEGER(KIND=JPIM), INTENT(IN) :: KNEMOCOMM

#ifdef WITH_NEMO

CALL NEMOGCMCOUP_INIT_IOSERVER_2( KNEMOCOMM )

#endif

END SUBROUTINE ININEMOIO2
