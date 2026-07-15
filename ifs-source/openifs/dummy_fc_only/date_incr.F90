! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE DATE_INCR (KDATE_IN, KSEC_INCR, KDATE_OUT, KSEC_OUT)
use parkind1, only:&
 & jpim
INTEGER (KIND=JPIM), INTENT(IN) :: KDATE_IN
INTEGER (KIND=JPIM), INTENT(IN) :: KSEC_INCR
INTEGER (KIND=JPIM), INTENT(OUT) :: KDATE_OUT
INTEGER (KIND=JPIM), INTENT(OUT) :: KSEC_OUT
call abor1("DATE_INCR should never be called with this build configuration - EXITING")
END SUBROUTINE DATE_INCR
