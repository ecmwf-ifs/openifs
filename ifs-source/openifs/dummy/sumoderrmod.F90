! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SUMODERRMOD(YDMODERR,YDRIP,CDNAMELIST)
USE YOMRIP, ONLY : TRIP
USE YOMMODERRMOD, ONLY : TMODERR
USE YOMMODERRMOD, ONLY : TMODERR
TYPE(TMODERR), INTENT(INOUT) :: YDMODERR
TYPE(TRIP), INTENT(IN) :: YDRIP
CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: CDNAMELIST
call abor1("sumoderrmod.F90 should never be called with OpenIFS - EXITING")
END SUBROUTINE SUMODERRMOD
