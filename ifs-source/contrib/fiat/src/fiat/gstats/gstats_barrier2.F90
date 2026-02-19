! (C) Copyright 2005- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE GSTATS_BARRIER2(KNUM)

USE EC_PARKIND  ,ONLY : JPIM

USE YOMGSTATS, ONLY : LBARRIER_STATS2,NBAR_PTR,NBAR2
USE MPL_MODULE  , ONLY : MPL_BARRIER

IMPLICIT NONE

INTEGER(KIND=JPIM) :: KNUM
INTEGER(KIND=JPIM) :: INUM

IF(LBARRIER_STATS2)THEN
  IF(NBAR_PTR(KNUM) == 0) THEN
    INUM=NBAR2
    NBAR2=NBAR2+1
    NBAR_PTR(KNUM)=INUM
  ENDIF
  INUM=NBAR_PTR(KNUM)
  CALL GSTATS(INUM,0)
  CALL MPL_BARRIER()
  CALL GSTATS(INUM,1)
ENDIF

END SUBROUTINE GSTATS_BARRIER2

