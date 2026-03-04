! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE MODGRIN(YDGEM,YDMCC,YDEPHY,KFIELDS,KPAR,PFPDGG)
USE YOMGEM , ONLY : TGEM
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOEPHY , ONLY : TEPHY
USE YOMMCC , ONLY : TMCC
TYPE(TGEM) ,INTENT(IN) :: YDGEM
TYPE(TEPHY) ,INTENT(INOUT) :: YDEPHY
TYPE(TMCC) ,INTENT(INOUT) :: YDMCC
INTEGER(KIND=JPIM),INTENT(IN) :: KFIELDS
INTEGER(KIND=JPIM),INTENT(OUT) :: KPAR(:)
REAL(KIND=JPRB) ,INTENT(OUT) :: PFPDGG(:,:)
call abor1("modgrin.F90 should never be called with OpenIFS - EXITING")
END SUBROUTINE MODGRIN
