! (C) Copyright 2005- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

INTERFACE

SUBROUTINE ABOR1(CDTEXT)
CHARACTER(LEN=*), INTENT(IN) :: CDTEXT
END SUBROUTINE ABOR1

SUBROUTINE ABOR1FL(CDFILE, KLINENUM, CDTEXT)
USE EC_PARKIND  ,ONLY : JPIM
CHARACTER(LEN=*), INTENT(IN) :: CDFILE,CDTEXT
INTEGER(KIND=JPIM), INTENT(IN) :: KLINENUM
END SUBROUTINE ABOR1FL

SUBROUTINE ABOR1_EXCEPTION_HANDLER()
END SUBROUTINE ABOR1_EXCEPTION_HANDLER

SUBROUTINE SET_ABOR1_EXCEPTION_HANDLER() BIND(C,name="set_abor1_exception_handler")
END SUBROUTINE SET_ABOR1_EXCEPTION_HANDLER

END INTERFACE
