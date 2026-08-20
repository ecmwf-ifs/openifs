! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SUDFI(YDEPHY,YDRIP,YDPHY,LDTENC)
USE YOMRIP , ONLY : TRIP
USE YOMPHY , ONLY : TPHY
USE YOEPHY , ONLY : TEPHY
TYPE(TEPHY),INTENT(INOUT) :: YDEPHY
TYPE(TPHY),INTENT(INOUT) :: YDPHY
TYPE(TRIP),INTENT(INOUT) :: YDRIP
LOGICAL,INTENT(IN) :: LDTENC
call abor1("SUDFI should never be called with this build configuration - EXITING")
END SUBROUTINE SUDFI
