! (C) Copyright 1989- Meteo-France.

SUBROUTINE DDHSND(PSND,KLEN,YD_IOSTREAM,YD_IOREQUEST)

!**** *ddhsnd * - Routine to gather DDH restart data and output

!     Purpose.
!     --------
!     To gather ddh restart data and write out

!**   Interface.
!     ----------
!        *CALL* *ddhsnd *

!        Explicit arguments :
!        --------------------
!        PSND   : data to send
!        KLEN   : length of data to send
!        KUNIT  : unit to write restart data

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
!      M.Hamrud   13-Oct-2008 Use IOSTREAM, one file per task, only FASTRES supported
!     ------------------------------------------------------------------

USE PARKIND1     , ONLY : JPIM     ,JPRB
USE YOMHOOK      , ONLY : LHOOK,   DR_HOOK, JPHOOK
USE IOSTREAM_MIX , ONLY : IO_PUT,TYPE_IOSTREAM , TYPE_IOREQUEST

IMPLICIT NONE

INTEGER(KIND=JPIM)  ,INTENT(IN)    :: KLEN
REAL(KIND=JPRB)     ,INTENT(IN)    :: PSND(KLEN)
TYPE(TYPE_IOSTREAM) ,INTENT(INOUT) :: YD_IOSTREAM
TYPE(TYPE_IOREQUEST),INTENT(INOUT) :: YD_IOREQUEST

INTEGER(KIND=JPIM) :: IBUF(1)  
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('DDHSND',0,ZHOOK_HANDLE)

IF( KLEN > 0 )THEN
  IBUF(1)=KLEN
  CALL IO_PUT(YD_IOSTREAM,YD_IOREQUEST,KR1=IBUF)
  CALL IO_PUT(YD_IOSTREAM,YD_IOREQUEST,PR1=PSND)
ENDIF

IF (LHOOK) CALL DR_HOOK('DDHSND',1,ZHOOK_HANDLE)
END SUBROUTINE DDHSND
