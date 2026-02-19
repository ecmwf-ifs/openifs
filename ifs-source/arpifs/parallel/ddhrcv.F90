! (C) Copyright 1989- Meteo-France.

SUBROUTINE DDHRCV(PRCV,KLEN,YD_IOSTREAM,YD_IOREQUEST)

!**** *ddhrcv * - Routine to input ddh restart data and distribute

!     Purpose.
!     --------
!     To input ddh restart data and distribute

!**   Interface.
!     ----------
!        *CALL* *ddhrcv *

!        Explicit arguments :
!        --------------------
!        PRCV   : data to send
!        KLEN   : length of data to receive
!        KUNIT  : unit to read restart data

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.   NONE.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      MPP Group *ECMWF*
!      Original : 96-08-27

!     Modifications.
!     --------------
!      Modified : 02-06-10 by D.Dent - Additional file for GP fields
!                          LFASTRES option and improved DDH handling
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      M.Hamrud  13-Oct-2008 Use IOSTREAM, one file per task, only FASTRES supported
!     ------------------------------------------------------------------

USE PARKIND1     , ONLY : JPIM     ,JPRB
USE YOMHOOK      , ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOMLUN       , ONLY : NULERR
USE IOSTREAM_MIX , ONLY : IO_GET,TYPE_IOSTREAM , TYPE_IOREQUEST

IMPLICIT NONE

INTEGER(KIND=JPIM)  , INTENT(IN)    :: KLEN
REAL(KIND=JPRB)     , INTENT(OUT)   :: PRCV(KLEN)
TYPE(TYPE_IOSTREAM) , INTENT(INOUT) :: YD_IOSTREAM
TYPE(TYPE_IOREQUEST), INTENT(INOUT) :: YD_IOREQUEST


INTEGER(KIND=JPIM) :: ILENRR, IBUF(1)  
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

IF (LHOOK) CALL DR_HOOK('DDHRCV',0,ZHOOK_HANDLE)

IF( KLEN > 0 )THEN
  CALL IO_GET(YD_IOSTREAM,YD_IOREQUEST,KR1=IBUF)
  ILENRR=IBUF(1)
  IF( KLEN /= ILENRR )THEN
    WRITE(NULERR,*) 'ILENRR =',ILENRR,' KLEN= ',KLEN
    CALL ABOR1('DDHRCV: WRONG LENGTH, INTERNAL ERROR.1')
  ENDIF
  CALL IO_GET(YD_IOSTREAM,YD_IOREQUEST,PR1=PRCV)
ENDIF


IF (LHOOK) CALL DR_HOOK('DDHRCV',1,ZHOOK_HANDLE)
END SUBROUTINE DDHRCV
