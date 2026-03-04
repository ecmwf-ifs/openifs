! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction
! 
! (C) Copyright 1989- Meteo-France.
! 

INTERFACE
SUBROUTINE SUAFN1SFX(CLFILEPGD2,KULOUT,KMASK_CNT,KMASK_BAS,CDSFXMASK,KPRE_CNT,CDNOMA,KMASK)

USE PARKIND1, ONLY : JPRB, JPIM
USE YOMHOOK, ONLY : LHOOK, DR_HOOK
IMPLICIT NONE
CHARACTER(LEN=*),    INTENT(IN)    :: CLFILEPGD2(:)
INTEGER (KIND=JPIM), INTENT(IN)    :: KULOUT
INTEGER (KIND=JPIM), INTENT(OUT)   :: KMASK_CNT
INTEGER (KIND=JPIM), INTENT(OUT)   :: KMASK_BAS
CHARACTER(LEN=*), ALLOCATABLE, INTENT (OUT) :: CDSFXMASK(:)
INTEGER (KIND=JPIM), INTENT(OUT)   :: KPRE_CNT
CHARACTER(LEN=*),    INTENT(OUT)   :: CDNOMA(:)
INTEGER (KIND=JPIM), INTENT(OUT)   :: KMASK(:)
END SUBROUTINE SUAFN1SFX
END INTERFACE
