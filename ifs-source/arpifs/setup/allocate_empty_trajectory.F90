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

SUBROUTINE ALLOCATE_EMPTY_TRAJECTORY(YDDIM,YDRIP)

!    Purpose:
!    --------
!     ??????

!    Author:
!    -------
!     ??????

!    Modifications:
!    --------------
!      K. Yessad (July 2014): Move some variables.
!      R. El Khatib 25-Oct-2018 workaround for DFI
! ----------------------------------------------------------------------------

USE YOMDIM  , ONLY : TDIM
USE PARKIND1, ONLY : JPIM, JPRB
USE YOMHOOK,  ONLY : LHOOK, DR_HOOK, JPHOOK

USE YOMRIP,   ONLY : TRIP
USE YOMTRAJ,  ONLY : TRAJEC, MSTEPTRAJW, MSTEPTRAJR

USE YOMINI, ONLY : LDFI
USE YOMVAR, ONLY : LJCDFI
USE YOMDFI, ONLY : NSTDFI, NSTDFIA
! ----------------------------------------------------------------------------

IMPLICIT NONE

TYPE(TDIM) , INTENT(IN) :: YDDIM
TYPE(TRIP)  ,INTENT(INOUT):: YDRIP

INTEGER(KIND=JPIM) :: JT, ISTOP
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ALLOCATE_EMPTY_TRAJECTORY',0,ZHOOK_HANDLE)
ASSOCIATE(NGPBLKS=>YDDIM%NGPBLKS, &
 & NSTART=>YDRIP%NSTART, NSTOP=>YDRIP%NSTOP)
! ----------------------------------------------------------------------------

IF (LDFI .OR. LJCDFI) THEN
  ISTOP=MAX(NSTOP,2*MAX(NSTDFI,NSTDFIA))
ELSE
  ISTOP=NSTOP
ENDIF

! Allocate trajectory timers and initialize them
IF (.NOT. ALLOCATED(MSTEPTRAJW)) THEN
  ALLOCATE(MSTEPTRAJW(NSTART:ISTOP))
  MSTEPTRAJW(:)=-99
ENDIF
IF (.NOT. ALLOCATED(MSTEPTRAJR)) THEN
  ALLOCATE(MSTEPTRAJR(NSTART:ISTOP))
  MSTEPTRAJR(:)=-99
ENDIF

! Allocate trajectory
! -------------------
  ALLOCATE (TRAJEC(0:ISTOP))

  ! Allocation only done for zero time-step
  ALLOCATE(TRAJEC(0)%CST(NGPBLKS))
  ALLOCATE(TRAJEC(0)%SLAG(NGPBLKS))
  ALLOCATE(TRAJEC(0)%PHYS(NGPBLKS))

  ! Now make the other time levels pointing to the 0 one.
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JT)
  DO JT=1,ISTOP
    TRAJEC(JT)%CST  => TRAJEC(0)%CST(:)
    TRAJEC(JT)%SLAG => TRAJEC(0)%SLAG(:)
    TRAJEC(JT)%PHYS => TRAJEC(0)%PHYS(:)
  ENDDO
!$OMP END PARALLEL DO

! ----------------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('ALLOCATE_EMPTY_TRAJECTORY',1,ZHOOK_HANDLE)
END SUBROUTINE ALLOCATE_EMPTY_TRAJECTORY
