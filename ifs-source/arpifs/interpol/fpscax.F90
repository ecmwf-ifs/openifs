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

SUBROUTINE FPSCAX(YDFPSTRUCT,KPROMA,KSTART,KEND,KLATA,KLONA,KS0)

!**** *FPSCAX*  - Compute the index of points in the SL halo
!             

!     Author. 
!     ------- 
!      Philippe Marguinaud *METEO FRANCE* 
!      Original : 01-10-2014

!     Description.
!     ------------
!     This routine is very similar to FPSCAW, expect that we
!     compute indexes for all points of the halo in a column 
!     near the target point


USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE EINT_MOD, ONLY : SL_STRUCT

IMPLICIT NONE

TYPE(SL_STRUCT),   INTENT(IN)  :: YDFPSTRUCT
INTEGER(KIND=JPIM),INTENT(IN)  :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)  :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)  :: KEND 
INTEGER(KIND=JPIM),INTENT(IN)  :: KLATA(KPROMA)
INTEGER(KIND=JPIM),INTENT(IN)  :: KLONA(KPROMA,2*YDFPSTRUCT%NSLWIDE)
INTEGER(KIND=JPIM),INTENT(OUT) :: KS0(KPROMA,2*YDFPSTRUCT%NSLWIDE)

INTEGER(KIND=JPIM) :: ILA, ILAG, ILO, JLO
INTEGER(KIND=JPIM) :: JLEN

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

IF (LHOOK) CALL DR_HOOK('FPSCAX',0,ZHOOK_HANDLE)

IF (YDFPSTRUCT%NSLWIDE  == 0) THEN
  CALL ABOR1 ('FPSCAX: INVALID YDFPSTRUCT%NSLWIDE')
ELSE
  DO JLEN=KSTART,KEND
    DO JLO = 1, 2*YDFPSTRUCT%NSLWIDE
      ILAG=KLATA(JLEN)+YDFPSTRUCT%NFRSTLOFF+(JLO-YDFPSTRUCT%NSLWIDE)
      ILA=MIN(MAX(ILAG,YDFPSTRUCT%NDGUNG),YDFPSTRUCT%NDGUXG)-YDFPSTRUCT%NFRSTLOFF
      ILO=MIN(MAX(KLONA(JLEN,JLO),YDFPSTRUCT%NDLUNG-1),YDFPSTRUCT%NDLUXG-1)
      KS0(JLEN,JLO)=YDFPSTRUCT%NSLOFF(ILA)+YDFPSTRUCT%NSLEXT(ILO,ILA)
    ENDDO
  ENDDO
ENDIF

IF (LHOOK) CALL DR_HOOK('FPSCAX',1,ZHOOK_HANDLE)

END SUBROUTINE FPSCAX

