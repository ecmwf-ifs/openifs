! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE DEALMOD(YDRADF,YDSLPHY,YDSPPT,YDSPPT_CONFIG)
USE YOMSLPHY , ONLY : TSLPHY
USE YOMRADF , ONLY : TRADF
USE YOMSPSDT , ONLY : TSPPT_CONFIG, TSPPT_DATA
TYPE(TRADF), INTENT(INOUT) :: YDRADF
TYPE(TSLPHY), INTENT(INOUT) :: YDSLPHY
TYPE(TSPPT_DATA), INTENT(INOUT) :: YDSPPT
TYPE(TSPPT_CONFIG), INTENT(INOUT) :: YDSPPT_CONFIG
call abor1("dealmod.F90 should never be called with OpenIFS - EXITING")
END SUBROUTINE DEALMOD
