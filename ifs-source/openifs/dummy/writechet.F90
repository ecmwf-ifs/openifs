! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE WRITECHET(YDRIP,KTSNO)
use parkind1 , only:&
 & jpim
USE YOMRIP , ONLY : TRIP
TYPE(TRIP) ,INTENT(INOUT) :: YDRIP
INTEGER(KIND=JPIM),INTENT(IN) :: KTSNO
call abor1("writechet.F90 should never be called with OpenIFS - EXITING")
END SUBROUTINE WRITECHET
