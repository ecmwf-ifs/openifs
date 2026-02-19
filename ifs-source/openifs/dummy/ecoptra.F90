! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE ECOPTRA(YDDIMV,YDRIP,CDCONF,YDSP)
USE YOMDIMV , ONLY : TDIMV
USE YOMRIP , ONLY : TRIP
USE SPECTRAL_FIELDS_DATA, ONLY: SPECTRAL_FIELD
TYPE(TDIMV), INTENT(INOUT) :: YDDIMV
TYPE(TRIP) ,INTENT(INOUT) :: YDRIP
CHARACTER(LEN=1) ,INTENT(IN) :: CDCONF
TYPE(SPECTRAL_FIELD), INTENT(IN) :: YDSP
call abor1("ecoptra.F90 should never be called with OpenIFS - EXITING")
END SUBROUTINE ECOPTRA
