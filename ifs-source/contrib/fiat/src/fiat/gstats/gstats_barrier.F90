! (C) Copyright 2005- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE GSTATS_BARRIER(KNUM)

USE EC_PARKIND  ,ONLY : JPIM

USE YOMGSTATS, ONLY : LBARRIER_STATS
USE MPL_MODULE  , ONLY : MPL_BARRIER

IMPLICIT NONE

INTEGER(KIND=JPIM) :: KNUM

IF(LBARRIER_STATS)THEN
  CALL GSTATS(KNUM,0)
  CALL MPL_BARRIER()
  CALL GSTATS(KNUM,1)
ENDIF

END SUBROUTINE GSTATS_BARRIER

