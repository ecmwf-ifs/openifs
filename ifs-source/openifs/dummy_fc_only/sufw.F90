! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SUFW(YDLAP,YDDIM,YDEGEO,YDELAP)
USE YOMLAP , ONLY : TLAP
USE YOMDIM , ONLY : TDIM
USE YEMGEO , ONLY : TEGEO
USE YEMLAP , ONLY : TLEP
TYPE(TLAP) , INTENT(IN) :: YDLAP
TYPE(TDIM) , INTENT(IN) :: YDDIM
TYPE(TEGEO) , INTENT(IN) :: YDEGEO
TYPE(TLEP) , INTENT(IN) :: YDELAP
call abor1("SUFW should never be called with this build configuration - EXITING")
END SUBROUTINE SUFW
