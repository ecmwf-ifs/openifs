! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SUEDIM(YDDIM,YDEDIM,KSEFRE,KSPECG,KSUPERSEDE)
USE YOMDIM , ONLY : TDIM
use parkind1 , only:&
 & jpim
USE YEMDIM , ONLY : TEDIM
TYPE(TDIM) , INTENT(INOUT) :: YDDIM
TYPE(TEDIM), INTENT(INOUT), TARGET :: YDEDIM
INTEGER(KIND=JPIM), INTENT(OUT) :: KSEFRE
INTEGER(KIND=JPIM), INTENT(OUT) :: KSPECG
INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: KSUPERSEDE
call abor1("suedim.F90 should never be called with OpenIFS - EXITING")
END SUBROUTINE SUEDIM
