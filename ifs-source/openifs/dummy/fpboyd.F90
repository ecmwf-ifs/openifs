! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE FPBOYD(PFPBSCAL,KFPXFLD,YDFPGEO_DEP,YDFPGIND,YDFPGEO,YDFPUSERGEO,PBUF_DEP,PBUF)
USE PARKIND1 ,ONLY : JPIM ,JPRB
USE TYPE_FPUSERGEO, ONLY : TFPUSERGEO
USE YOMFPGEO, ONLY : TFPGEO
USE YOMFPGIND, ONLY : TFPGIND
REAL(KIND=JPRB) ,INTENT(IN) :: PFPBSCAL
INTEGER(KIND=JPIM), INTENT(IN) :: KFPXFLD
TYPE (TFPGEO), INTENT(IN) :: YDFPGEO_DEP
TYPE (TFPGIND), INTENT(IN) :: YDFPGIND
TYPE (TFPGEO), INTENT(IN) :: YDFPGEO
TYPE (TFPUSERGEO) ,INTENT(IN) :: YDFPUSERGEO(:)
REAL(KIND=JPRB) ,INTENT(IN) :: PBUF_DEP(:,:,:)
REAL(KIND=JPRB) ,INTENT(OUT) :: PBUF(:,:,:)
call abor1("fpboyd.F90 should never be called with OpenIFS - EXITING")
END SUBROUTINE FPBOYD
