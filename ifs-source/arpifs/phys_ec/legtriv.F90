! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE LEGTRIV (PSIN,KCP,KDIM,PALP,KLON)

!**** *LEGTRIV* - *LEGENDRE FUNCTIONS FOR A TRIANGULAR TRUNCATION.

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

!     MODIFICATIONS.
!     --------------
!     J.HAGUE      01/07/03   Vector version of LEGTRI created for bit
!                             reproducibility
!        M.Hamrud      01-Oct-2003 CY28 Cleaning

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KDIM 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSIN(KLON) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KCP 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PALP(KDIM,KLON) 
INTEGER(KIND=JPIM) :: IC, ICP, II, IM, IM2, JM1, JN, JL

!REAL_B :: ZE1
REAL(KIND=JPRB) :: ZE1(KCP,KCP)
!REAL_B :: ZE2
REAL(KIND=JPRB) :: ZE2(KCP,KCP)

REAL(KIND=JPRB) :: Z2M, ZCOS, ZF1M, ZF2M, ZM, ZN, ZN2, ZRE1, ZSIN
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

!*         1.     PRELIMINARY SETTING.
!                 ----------- --------

IF (LHOOK) CALL DR_HOOK('LEGTRIV',0,ZHOOK_HANDLE)
ICP=KCP
IC=ICP-1
DO JM1=1,ICP
  IM=JM1-1
  ZM=IM
  Z2M=ZM+ZM
  ZRE1=SQRT(Z2M+3._JPRB)
  IM2=IM+2
  ZE1(IM2-1,JM1)=1.0_JPRB/ZRE1
  DO JN=IM2,IC
    ZN=JN
    ZN2=ZN**2
    ZE2(JN,JM1)=SQRT((4._JPRB*ZN2-1.0_JPRB)/(ZN2-ZM**2))
    ZE1(JN,JM1)=1.0_JPRB/ZE2(JN,JM1)
  ENDDO
ENDDO

DO JL=1,KLON

  ZSIN=PSIN(JL)
  ICP=KCP

!     ------------------------------------------------------------------

!*         2.     COMPUTATIONS.
!                 -------------

  IC=ICP-1
  ZCOS=SQRT(1.0_JPRB-ZSIN**2)
  II=2
  PALP(1,JL)=1.0_JPRB
  ZF1M=SQRT(3._JPRB)
  PALP(2,JL)=ZF1M*ZSIN
  DO JM1=1,ICP
    IM=JM1-1
    ZM=IM
    Z2M=ZM+ZM
    ZRE1=SQRT(Z2M+3._JPRB)
!   ZE1=_ONE_/ZRE1
    IF(IM == 0) GO TO 201
    ZF2M=ZF1M*ZCOS/SQRT(Z2M)
    ZF1M=ZF2M*ZRE1
    II=II+1
    PALP(II,JL)=ZF2M
    IF(IM == IC) GO TO 203
    II=II+1
    PALP(II,JL)=ZF1M*ZSIN
    IF(JM1 == IC) GO TO 203
    201 CONTINUE
    IM2=IM+2
    DO JN=IM2,IC
!     ZN=JN
!     ZN2=ZN**2
!     ZE2=SQRT((4._JPRB*ZN2-_ONE_)/(ZN2-ZM**2))
      II=II+1
!     PALP(II,JL)=ZE2*(ZSIN*PALP(II-1,JL)-ZE1*PALP(II-2,JL))
      PALP(II,JL)=ZE2(JN,JM1)*(ZSIN*PALP(II-1,JL)-ZE1(JN-1,JM1)*PALP(II-2,JL)) 
!     ZE1=_ONE_/ZE2(JN,JM1)
    ENDDO
    203 CONTINUE
  ENDDO

ENDDO

!     ------------------------------------------------------------------

!*         3.     RETURN.
!                 -------

IF (LHOOK) CALL DR_HOOK('LEGTRIV',1,ZHOOK_HANDLE)
END SUBROUTINE LEGTRIV

