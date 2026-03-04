! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE ETRANSDIRH(YDGEOMETRY,YDGFL,YDGMV,CDCONF,KNFTHER,YDSP)
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMGFL , ONLY : TGFL
USE YOMGMV , ONLY : TGMV
use parkind1 , only:&
 & jpim
USE SPECTRAL_FIELDS_MOD
TYPE(GEOMETRY), INTENT(IN) :: YDGEOMETRY
TYPE(TGFL) , INTENT(INOUT) :: YDGFL
TYPE(TGMV) , INTENT(INOUT) :: YDGMV
CHARACTER(LEN=1),INTENT(IN) :: CDCONF
INTEGER(KIND=JPIM) ,INTENT(IN) :: KNFTHER
TYPE(SPECTRAL_FIELD), INTENT(INOUT) :: YDSP
call abor1("etransdirh.F90 should never be called with OpenIFS - EXITING")
END SUBROUTINE ETRANSDIRH
