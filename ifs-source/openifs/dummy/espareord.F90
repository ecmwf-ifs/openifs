! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE ESPAREORD(YDDIM,YDEDIM,YDLEP,KFLEV,PSPFILE,PSPBUF,LD_FILE_TO_MODEL)
USE YOMDIM , ONLY : TDIM
USE PARKIND1 ,ONLY : JPIM ,JPRB
USE YEMDIM , ONLY : TEDIM
USE YEMLAP , ONLY : TLEP
TYPE(TDIM) , INTENT(IN) :: YDDIM
TYPE(TEDIM), INTENT(IN) :: YDEDIM
TYPE(TLEP), INTENT(IN) :: YDLEP
INTEGER(KIND=JPIM),INTENT(IN) :: KFLEV
REAL(KIND=JPRB) ,INTENT(INOUT) :: PSPFILE(YDDIM%NSEFRE,KFLEV)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PSPBUF(YDDIM%NSPEC2G,KFLEV)
LOGICAL ,INTENT(IN) :: LD_FILE_TO_MODEL
call abor1("espareord.F90 should never be called with OpenIFS - EXITING")
END SUBROUTINE ESPAREORD
