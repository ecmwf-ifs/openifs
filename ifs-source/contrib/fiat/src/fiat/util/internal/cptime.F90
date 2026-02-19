! (C) Copyright 2005- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

!RJ: make interfaceable; generalization
SUBROUTINE CPTIME(PVCP,PTCP)
USE EC_PARKIND, ONLY : JPRD, JPIM
IMPLICIT NONE
REAL(KIND=JPRD) :: PVCP
REAL(KIND=JPRD) :: PTCP
!
#if defined (NEWTIMER)
! this routine should work better with OpenMP
! But doesn't work on Cray - and in any case does not return 
! CPU time for all threads combined 
INTEGER(KIND=JPIM),SAVE :: IFIRST=1
INTEGER(KIND=JPIM),SAVE :: KFIRST,KTPS
INTEGER(KIND=JPIM) :: KTICK
!     Usage of Fortran95 intrinsic function system_clock for ELAPSED time,
!     thus taking into account the parallelism if inside an open-mp region. REK
IF(IFIRST.EQ.1) THEN
  IFIRST=0
  CALL SYSTEM_CLOCK(KFIRST,KTPS)
  PVCP=0.0_JPRD
  PTCP=PVCP
ELSE
  CALL SYSTEM_CLOCK(KTICK)
  PVCP=REAL(KTICK-KFIRST,KIND=JPRD)/REAL(KTPS,KIND=JPRD)
  PTCP=PVCP
ENDIF
#else
INTEGER(KIND=JPIM),SAVE :: IFIRST=0
REAL(KIND=JPRD),SAVE :: ZFIRST
REAL(KIND=JPRD) :: ZSEC
!     Usage of Fortran95 intrinsic function for CPU timing.
IF(IFIRST.EQ.0) THEN
  IFIRST=1
  CALL CPU_TIME(ZFIRST)
  PVCP=0.0_JPRD
  PTCP=PVCP
ELSE
  CALL CPU_TIME(ZSEC)
  PVCP=ZSEC-ZFIRST
  PTCP=PVCP
ENDIF
#endif
!
RETURN
END SUBROUTINE CPTIME
