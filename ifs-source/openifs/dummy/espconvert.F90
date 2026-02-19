! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE ESPCONVERT(YDGEOMETRY,LDMODEL_TO_FILE,YDSP)
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE SPECTRAL_FIELDS_DATA, ONLY: SPECTRAL_FIELD
TYPE(GEOMETRY), INTENT(IN) :: YDGEOMETRY
LOGICAL ,INTENT(IN) :: LDMODEL_TO_FILE
TYPE(SPECTRAL_FIELD), INTENT(INOUT) :: YDSP
call abor1("espconvert.F90 should never be called with OpenIFS - EXITING")
END SUBROUTINE ESPCONVERT
