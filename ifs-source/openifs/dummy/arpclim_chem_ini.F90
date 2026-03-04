! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE ARPCLIM_CHEM_INI(YDGEOMETRY,YGFL,YDCHEM)
USE YOM_YGFL, ONLY : TYPE_GFLD
USE GEOMETRY_MOD, ONLY : GEOMETRY
USE YOMCHEM , ONLY : TCHEM
TYPE(TYPE_GFLD),INTENT(INOUT) :: YGFL
TYPE(GEOMETRY) ,INTENT(IN) :: YDGEOMETRY
TYPE(TCHEM), INTENT(INOUT) :: YDCHEM
call abor1("arpclim_chem_ini.F90 should never be called with OpenIFS - EXITING")
END SUBROUTINE ARPCLIM_CHEM_INI
