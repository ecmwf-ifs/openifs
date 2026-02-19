! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SUEMODJK(YDSTA,YDDIMV,YDDIMF)
USE YOMDIMV , ONLY : TDIMV
USE YOMDIMF , ONLY : TDIMF
USE YOMSTA , ONLY : TSTA
TYPE(TSTA), INTENT(INOUT) :: YDSTA
TYPE(TDIMV), INTENT(INOUT) :: YDDIMV
TYPE(TDIMF) ,INTENT(INOUT) :: YDDIMF
call abor1("suemodjk.F90 should never be called with OpenIFS - EXITING")
END SUBROUTINE SUEMODJK
