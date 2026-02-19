! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE ECOMMSPNORM(YDGEOMETRY,KIOMASTER,PMEAN,YDSP)
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 ,ONLY : JPIM ,JPRB
USE SPECTRAL_FIELDS_DATA, ONLY: SPECTRAL_FIELD
USE MPL_MODULE
TYPE(GEOMETRY), INTENT(IN) :: YDGEOMETRY
INTEGER(KIND=JPIM),INTENT(IN) :: KIOMASTER
REAL(KIND=JPRB) ,INTENT(OUT) :: PMEAN(YDGEOMETRY%YRDIMV%NFLEVG)
TYPE(SPECTRAL_FIELD), INTENT(IN) :: YDSP
call abor1("ecommspnorm.F90 should never be called with OpenIFS - EXITING")
END SUBROUTINE ECOMMSPNORM
