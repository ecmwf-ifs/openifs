! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE OPEN_OUTPUT_LFA(YDPHY2)
USE YOMPHY2 , ONLY : TPHY2
TYPE(TPHY2),INTENT(INOUT):: YDPHY2
call abor1("open_output_lfa.F90 should never be called with OpenIFS - EXITING")
END SUBROUTINE OPEN_OUTPUT_LFA
