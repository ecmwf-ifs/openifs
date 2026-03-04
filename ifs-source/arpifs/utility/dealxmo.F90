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

SUBROUTINE DEALXMO(YDGFL,YDGMV)

!**** *DEALXMO* - Deallocate GMV,GMVS,GFL,GFLSLP,GFLPT,GFLPC

!     Purpose.
!     -------
!           Deallocate GMV,GMVS,GFL,GFLSLP,GFLPT,GFLPC

!**   Interface.
!**   ---------
!**         DEALXMO is called from routine CNT1

!     Author.
!     -------
!        J. HAGUE *ECMWF*

!     Modifications.
!     --------------
!        Original :  07-Dec-2004
!     T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!------------------------------------------------------


USE YOMGFL   , ONLY : TGFL
USE YOMGMV   , ONLY : TGMV
USE PARKIND1 , ONLY : JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
IMPLICIT NONE

TYPE(TGFL) , INTENT(INOUT) :: YDGFL
TYPE(TGMV) , INTENT(INOUT) :: YDGMV
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('DEALXMO',0,ZHOOK_HANDLE)

IF (ALLOCATED(YDGMV%GMV))    DEALLOCATE(YDGMV%GMV)
IF (ALLOCATED(YDGMV%GMVS))   DEALLOCATE(YDGMV%GMVS)
IF (ALLOCATED(YDGFL%GFL))    DEALLOCATE(YDGFL%GFL)
IF (ALLOCATED(YDGFL%GFLSLP)) DEALLOCATE(YDGFL%GFLSLP)
IF (ALLOCATED(YDGFL%GFLPT))  DEALLOCATE(YDGFL%GFLPT)
IF (ALLOCATED(YDGFL%GFLPC))  DEALLOCATE(YDGFL%GFLPC)

IF (LHOOK) CALL DR_HOOK('DEALXMO',1,ZHOOK_HANDLE)
END SUBROUTINE DEALXMO

