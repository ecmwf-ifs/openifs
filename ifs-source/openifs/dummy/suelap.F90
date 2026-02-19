! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SUELAP(YDDIM,YDEDIM,YDLAP,YDLEP,YDEGEO)
USE YOMDIM , ONLY : TDIM
USE YEMGEO , ONLY : TEGEO
USE YOMLAP , ONLY : TLAP
USE YEMLAP , ONLY : TLEP
USE YEMDIM , ONLY : TEDIM
TYPE(TDIM) , INTENT(INOUT) :: YDDIM
TYPE(TEDIM), INTENT(INOUT) :: YDEDIM
TYPE(TLAP), INTENT(INOUT) :: YDLAP
TYPE(TLEP), INTENT(INOUT) :: YDLEP
TYPE(TEGEO), INTENT(INOUT) :: YDEGEO
call abor1("suelap.F90 should never be called with OpenIFS - EXITING")
END SUBROUTINE SUELAP
