! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE ESPCM(YDGEOMETRY,YDMODEL,YDFIELDS,CDCONF,LDESPCL,LDIAU,KCPLSPEC2,KCPLSPEC2V)
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE TYPE_MODEL, ONLY : MODEL
use parkind1 , only:&
 & jpim
USE FIELDS_MOD, ONLY: FIELDS
TYPE(GEOMETRY), INTENT(IN) :: YDGEOMETRY
TYPE(MODEL),INTENT(INOUT):: YDMODEL
TYPE(FIELDS),INTENT(INOUT):: YDFIELDS
CHARACTER(LEN=1), INTENT(IN) :: CDCONF
LOGICAL, INTENT(IN) :: LDESPCL
LOGICAL, INTENT(IN) :: LDIAU
INTEGER(KIND=JPIM), INTENT(IN) :: KCPLSPEC2
INTEGER(KIND=JPIM), INTENT(IN) :: KCPLSPEC2V
call abor1("espcm.F90 should never be called with OpenIFS - EXITING")
END SUBROUTINE ESPCM
