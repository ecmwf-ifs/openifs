! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SUEDYN(YDEDYN,YDEGEO,YDGEM)
USE YOMGEM , ONLY : TGEM
USE YEMDYN , ONLY : TEDYN
USE YEMGEO , ONLY : TEGEO
TYPE(TEDYN), TARGET,INTENT(INOUT) :: YDEDYN
TYPE(TEGEO), INTENT(INOUT) :: YDEGEO
TYPE(TGEM) , INTENT(INOUT) :: YDGEM
call abor1("suedyn.F90 should never be called with OpenIFS - EXITING")
END SUBROUTINE SUEDYN
