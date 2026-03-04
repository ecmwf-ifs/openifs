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

SUBROUTINE GETNEMO1WAY()
!
!**** *GETNEMO1WAY*  - .
!
!     Purpose.
!     --------
!       UPDATE FIELDS IN NEMO IF 1WAY COUPLING
!
!**   Interface.
!     ----------
!       *CALL*  *GETNEMO1WAY*
!
!     Input:
!     -----
!
!     Output:
!     ------
!
!     Method:
!     ------
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
!     -----------------------------------------------------------
   
USE PARKIND1 , ONLY : JPRD, JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN   , ONLY : NULOUT
USE MPL_MODULE, ONLY : MPL_COMM
USE YOMMP0   , ONLY : MYPROC, NPROC

IMPLICIT NONE

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     -----------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GETNEMO1WAY',0,ZHOOK_HANDLE)

#ifdef WITH_NEMO

CALL NEMOGCMCOUP_GET_1WAY( MYPROC-1, NPROC, MPL_COMM )

#else

CALL ABOR1('ININEMO: COMPILED WITHOUT WITH_NEMO')

#endif

IF (LHOOK) CALL DR_HOOK('GETNEMO1WAY',1,ZHOOK_HANDLE)

!     -----------------------------------------------------------

END SUBROUTINE GETNEMO1WAY

