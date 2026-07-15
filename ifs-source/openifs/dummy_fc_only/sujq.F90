! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SUJQ(YDERRMOD,YDGEOMETRY,YDSPCTLMODERR)
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMJQ
USE SPECTRAL_FIELDS_MOD, ONLY : SPECTRAL_FIELD
TYPE(ERRMOD_STRUCT), INTENT(INOUT) :: YDERRMOD
TYPE(GEOMETRY), INTENT(IN) :: YDGEOMETRY
TYPE(SPECTRAL_FIELD), INTENT(IN) :: YDSPCTLMODERR
call abor1("SUJQ should never be called with this build configuration - EXITING")
END SUBROUTINE SUJQ
