! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SUJBWAVELET0(YDDIM,YD_JB_STRUCT,CDFILE)
USE YOMDIM , ONLY : TDIM
USE YOMJG, ONLY : TYPE_JB_STRUCT
TYPE(TDIM) , INTENT(IN) :: YDDIM
TYPE(TYPE_JB_STRUCT), INTENT(INOUT) :: YD_JB_STRUCT
CHARACTER(LEN=*) , INTENT(IN) :: CDFILE
call abor1("sujbwavelet0.F90 should never be called with OpenIFS - EXITING")
END SUBROUTINE SUJBWAVELET0
