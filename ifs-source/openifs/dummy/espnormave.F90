! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE ESPNORMAVE(KLEVXG,PNORM,PMEAN,PAVE)
USE PARKIND1 ,ONLY : JPIM ,JPRB
INTEGER(KIND=JPIM),INTENT(IN) :: KLEVXG
REAL(KIND=JPRB) ,INTENT(INOUT) :: PNORM(KLEVXG)
REAL(KIND=JPRB) ,INTENT(IN) :: PMEAN(KLEVXG)
REAL(KIND=JPRB) ,INTENT(OUT) :: PAVE
call abor1("espnormave.F90 should never be called with OpenIFS - EXITING")
END SUBROUTINE ESPNORMAVE
