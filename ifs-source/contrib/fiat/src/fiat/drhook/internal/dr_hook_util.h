! (C) Copyright 2005- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

INTERFACE
SUBROUTINE DR_HOOK_UTIL(LDHOOK,CDNAME,KCASE,PKEY,CDFILENAME,KSIZEINFO)
USE EC_PARKIND  ,ONLY : JPIM     ,JPRD
IMPLICIT NONE
LOGICAL, INTENT(INOUT)      :: LDHOOK
CHARACTER(LEN=*),INTENT(IN) :: CDNAME,CDFILENAME
INTEGER(KIND=JPIM),INTENT(IN) :: KCASE,KSIZEINFO
REAL(KIND=JPRD),INTENT(INOUT) :: PKEY
END SUBROUTINE DR_HOOK_UTIL
END INTERFACE
