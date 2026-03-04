! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

MODULE ORDER_INDEPENDENT_SUMMATION_MOD

! Dummy version for single column model

USE PARKIND1  ,ONLY : JPIM     ,JPRB

PUBLIC

!  single column model only needs ORDER_INDEP_GLOBAL_SUM

 CONTAINS

FUNCTION ORDER_INDEP_GLOBAL_SUM(PIN,KNG)

IMPLICIT NONE

REAL(KIND=JPRB) :: ORDER_INDEP_GLOBAL_SUM
REAL(KIND=JPRB), INTENT(IN) :: PIN(:)
INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KNG

PRINT *,'DUMMY ORDER_INDEP_GLOBAL_SUM CALLED'
ORDER_INDEP_GLOBAL_SUM = SUM(PIN)

END FUNCTION

END MODULE