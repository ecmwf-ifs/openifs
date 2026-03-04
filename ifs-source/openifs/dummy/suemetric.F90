! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SUEMETRIC(YDDIM,YDEDIM,YDLEP,KSMAX,PMET,PMETDER,PMETKE)
USE YOMDIM , ONLY : TDIM
USE PARKIND1 ,ONLY : JPIM ,JPRB
USE YEMLAP , ONLY : TLEP
USE YEMDIM , ONLY : TEDIM
TYPE(TDIM) , INTENT(IN) :: YDDIM
TYPE(TEDIM), INTENT(IN) :: YDEDIM
TYPE(TLEP), INTENT(IN) :: YDLEP
INTEGER(KIND=JPIM),INTENT(OUT) :: KSMAX
REAL(KIND=JPRB) ,INTENT(OUT) :: PMET(YDDIM%NSPECG)
REAL(KIND=JPRB) ,INTENT(OUT) :: PMETDER(YDDIM%NSPECG)
REAL(KIND=JPRB) ,INTENT(OUT) :: PMETKE(YDDIM%NSPECG)
call abor1("suemetric.F90 should never be called with OpenIFS - EXITING")
END SUBROUTINE SUEMETRIC
