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

SUBROUTINE LEGTRI (PSIN,KCP,KDIM,PALP)

!**** *LEGTRI* - *LEGENDRE FUNCTIONS FOR A TRIANGULAR TRUNCATION.

!     J.F.GELEYN     E.C.M.W.F.     03/06/82.

!     PURPOSE.
!     --------

!          THIS ROUTINE COMPUTES THE VALUES *PALP* FOR THE ARGUMENT
!     *PSIN* OF THE NORMALISED *LEGENDRE ASSOCIATED FUNCTIONS IN THE
!     ORDER ((JN1=JM1,KCP),JM1=1,KCP) FOR JN=JN1-1 AND JM=JM1-1 .

!**   INTERFACE.
!     ----------

!          *LEGTRI* IS CALLED FROM *RADMOD*.
!          THERE ARE THREE DUMMY ARGUMENTS: *PSIN* IS THE SINE OF
!     LATITUDE.
!                                           *KCP* IS ONE PLUS THE LIMIT
!     WAVE NUMBER.
!                                           *PALP* IS THE ARRAY OF THE
!     RESULTS.

!     METHOD.
!     -------

!          SIMPLE RECURENCE FORMULA.

!     EXTERNALS.
!     ----------

!          NONE.

!     REFERENCE.
!     ----------

!          NONE.

USE PARKIND1  ,ONLY : JPIM     ,JPRB, JPRD
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KDIM 
REAL(KIND=JPRD)   ,INTENT(IN)    :: PSIN 
INTEGER(KIND=JPIM),INTENT(IN)    :: KCP 
REAL(KIND=JPRD)   ,INTENT(OUT)   :: PALP(KDIM) 
INTEGER(KIND=JPIM) :: IC, ICP, II, IM, IM2, JM1, JN

REAL(KIND=JPRD) :: Z2M, ZCOS, ZE1, ZE2, ZF1M, ZF2M, ZM, ZN, ZN2, ZRE1, ZSIN
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

!*         1.     PRELIMINARY SETTING.
!                 ----------- --------

IF (LHOOK) CALL DR_HOOK('LEGTRI',0,ZHOOK_HANDLE)
ZSIN=PSIN
ICP=KCP

!     ------------------------------------------------------------------

!*         2.     COMPUTATIONS.
!                 -------------

IC=ICP-1
ZCOS=SQRT(1.0_JPRD-ZSIN**2)
II=2
PALP(1)=1.0_JPRD
ZF1M=SQRT(3._JPRD)
PALP(2)=ZF1M*ZSIN
DO JM1=1,ICP
  IM=JM1-1
  ZM=IM
  Z2M=ZM+ZM
  ZRE1=SQRT(Z2M+3._JPRD)
  ZE1=1.0_JPRD/ZRE1
  IF(IM /= 0) THEN
    ZF2M=ZF1M*ZCOS/SQRT(Z2M)
    ZF1M=ZF2M*ZRE1
    II=II+1
    PALP(II)=ZF2M
    IF(IM == IC) CYCLE
    II=II+1
    PALP(II)=ZF1M*ZSIN
    IF(JM1 == IC) CYCLE
  ENDIF
  IM2=IM+2
  DO JN=IM2,IC
    ZN=JN
    ZN2=ZN**2
    ZE2=SQRT((4._JPRD*ZN2-1.0_JPRD)/(ZN2-ZM**2))
    II=II+1
    PALP(II)=ZE2*(ZSIN*PALP(II-1)-ZE1*PALP(II-2))
    ZE1=1.0_JPRD/ZE2
  ENDDO
ENDDO

!     ------------------------------------------------------------------

!*         3.     RETURN.
!                 -------

IF (LHOOK) CALL DR_HOOK('LEGTRI',1,ZHOOK_HANDLE)
END SUBROUTINE LEGTRI

