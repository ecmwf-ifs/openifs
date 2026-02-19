! (C) Copyright 2005- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE MPL_WAIT_MOD

!**** MPL_WAIT Waits for completion

!     Purpose.
!     --------
!     Returns control when the operation(s) identified by the request
!     is completed.
!     Normally used in conjunction with non-blocking buffering type

!**   Interface.
!     ----------
!        CALL MPL_WAIT

!        Input required arguments :
!        -------------------------
!           KREQUEST -  array or scalar containing
!                       Communication request(s)
!                       as provided by MPL_RECV or MPL_SEND

!        Input optional arguments :
!        -------------------------
!           CDSTRING -  Character string for ABORT messages
!                       used when KERROR is not provided

!        Output required arguments :
!        -------------------------
!           none

!        Output optional arguments :
!        -------------------------
!           KOUNT    -  must be the same size and shape as KREQUEST
!                       contains number of items sent/received
!           KBYTES   -  number of bytes in a single element in all KOUNTs
!                       *must* be supplied with KOUNT
!                       KBYTES normally determited robustly by MPL_BYTES
!           KERROR   -  return error code.     If not supplied, 
!                       MPL_WAIT aborts when an error is detected.
!     Author.
!     -------
!        D.Dent, M.Hamrud     ECMWF

!     Modifications.
!     --------------
!        Original: 2000-09-01
!        J. Hague: 2005-04-25  WAITALL replaced by WAIT loop
!      F. Vana  05-Mar-2015  Support for single precision
!  S. Saarinen  17-Feb-2017  Removed PBUF argument (not realy needed)
!                            KREQUEST must be INOUT (as per MPI_Wait)
!                            MPL_WAITS calls MPI_Waitall unless MPI1
!  S. Saarinen  01-Mar-2017  Added KBYTES

!     ------------------------------------------------------------------

USE EC_PARKIND  ,ONLY : JPRD, JPIM, JPRM, JPIB

USE MPL_MPIF
USE MPL_DATA_MODULE
USE MPL_MESSAGE_MOD

IMPLICIT NONE

PRIVATE

INTERFACE MPL_WAIT
MODULE PROCEDURE MPL_WAITS, MPL_WAIT1
END INTERFACE

PUBLIC MPL_WAIT

CONTAINS 


SUBROUTINE MPL_WAITS(KREQUEST,KOUNT,KBYTES,KERROR,CDSTRING)

#ifdef USE_8_BYTE_WORDS
  USE MPI4TO8, ONLY : &
    MPI_WAITALL => MPI_WAITALL8, MPI_GET_COUNT => MPI_GET_COUNT8, &
    MPI_WAIT => MPI_WAIT8
#endif

INTEGER(KIND=JPIM),INTENT(INOUT)  :: KREQUEST(:)
CHARACTER*(*),INTENT(IN),OPTIONAL :: CDSTRING
INTEGER(KIND=JPIM),INTENT(OUT),OPTIONAL :: KOUNT(SIZE(KREQUEST))
INTEGER(KIND=JPIM),INTENT(IN),OPTIONAL :: KBYTES
INTEGER(KIND=JPIM),INTENT(OUT),OPTIONAL :: KERROR
INTEGER(KIND=JPIM) :: IWAITERR,ICOUNTERR,JL,IREQLEN,ICOUNT,IW
INTEGER(KIND=JPIM) :: IWAIT_STATUS(MPI_STATUS_SIZE,SIZE(KREQUEST))
LOGICAL :: LLABORT
LLABORT=.TRUE.
IWAITERR=0
ICOUNTERR=0

IF(MPL_NUMPROC < 1) CALL MPL_MESSAGE( &
  & CDMESSAGE='MPL_WAITS: MPL NOT INITIALISED ',LDABORT=LLABORT) 

IREQLEN=SIZE(KREQUEST)
#ifndef MPI1
CALL MPI_WAITALL(IREQLEN,KREQUEST,IWAIT_STATUS,IWAITERR)
#else
DO JL=1,IREQLEN
  CALL MPI_WAIT(KREQUEST(JL),IWAIT_STATUS(1,JL),IW)
  IWAITERR=MAX(IWAITERR,IW)
ENDDO
#endif

IF(PRESENT(KOUNT))THEN
  IF (.not.PRESENT(KBYTES)) THEN
    CALL MPL_MESSAGE( &
    & CDMESSAGE='MPL_WAIT: KBYTES MUST BE PRESENT WITH KOUNT ', &
    & CDSTRING=CDSTRING,LDABORT=LLABORT)
  ENDIF
  IF(SIZE(KOUNT) /= IREQLEN) THEN
    CALL MPL_MESSAGE( &
    & CDMESSAGE='MPL_WAIT: KOUNT AND KREQUEST INCONSISTENT ', &
    & CDSTRING=CDSTRING,LDABORT=LLABORT)
  ENDIF
  DO JL=1,IREQLEN
    CALL MPI_GET_COUNT(IWAIT_STATUS(1,JL),INT(MPI_BYTE),KOUNT(JL),ICOUNTERR)
    KOUNT(JL) = KOUNT(JL) / KBYTES
  ENDDO
ENDIF

IF(PRESENT(KERROR))THEN
  KERROR=IWAITERR+ICOUNTERR
ELSE IF(IWAITERR /= 0) THEN
  CALL MPL_MESSAGE(IWAITERR,'MPL_WAITS_WAITING',CDSTRING,LDABORT=LLABORT)
ELSE IF(ICOUNTERR /= 0) THEN
  CALL MPL_MESSAGE(ICOUNTERR,'MPL_WAITS_COUNT',CDSTRING,LDABORT=LLABORT)
ENDIF

RETURN
END SUBROUTINE MPL_WAITS


SUBROUTINE MPL_WAIT1(KREQUEST,KOUNT,KBYTES,KERROR,CDSTRING)

#ifdef USE_8_BYTE_WORDS
  USE MPI4TO8, ONLY : &
    MPI_WAITALL => MPI_WAITALL8, MPI_GET_COUNT => MPI_GET_COUNT8
#endif

INTEGER(KIND=JPIM),INTENT(INOUT)  :: KREQUEST
CHARACTER*(*),INTENT(IN),OPTIONAL :: CDSTRING
INTEGER(KIND=JPIM),INTENT(OUT),OPTIONAL :: KOUNT
INTEGER(KIND=JPIM),INTENT(IN),OPTIONAL :: KBYTES
INTEGER(KIND=JPIM),INTENT(OUT),OPTIONAL :: KERROR
INTEGER(KIND=JPIM) :: IWAITERR,ICOUNTERR,JL,IREQLEN,ICOUNT
INTEGER(KIND=JPIM) :: IWAIT_STATUS(MPI_STATUS_SIZE)
LOGICAL :: LLABORT
LLABORT=.TRUE.
IWAITERR=0
ICOUNTERR=0

IF(MPL_NUMPROC < 1) CALL MPL_MESSAGE( &
  & CDMESSAGE='MPL_WAIT: MPL NOT INITIALISED ',LDABORT=LLABORT) 

CALL MPI_WAIT(KREQUEST,IWAIT_STATUS,IWAITERR)

IF(PRESENT(KOUNT))THEN
  IF (.not.PRESENT(KBYTES)) THEN
    CALL MPL_MESSAGE( &
    & CDMESSAGE='MPL_WAIT: KBYTES MUST BE PRESENT WITH KOUNT ', &
    & CDSTRING=CDSTRING,LDABORT=LLABORT)
  ENDIF
  CALL MPI_GET_COUNT(IWAIT_STATUS(1),INT(MPI_BYTE),KOUNT,ICOUNTERR)
  KOUNT = KOUNT / KBYTES
ENDIF

IF(PRESENT(KERROR))THEN
  KERROR=IWAITERR+ICOUNTERR
ELSE IF(IWAITERR /= 0) THEN
  CALL MPL_MESSAGE(IWAITERR,'MPL_WAIT_WAITING',CDSTRING,LDABORT=LLABORT)
ELSE IF(ICOUNTERR /= 0) THEN
  CALL MPL_MESSAGE(ICOUNTERR,'MPL_WAIT_COUNT',CDSTRING,LDABORT=LLABORT)
ENDIF

RETURN
END SUBROUTINE MPL_WAIT1

END MODULE MPL_WAIT_MOD
