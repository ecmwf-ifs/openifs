! (C) Copyright 2005- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE MPL_PROBE_MOD

!**** MPL_PROBE - Check for incoming message

!     Purpose.
!     --------
!     Look for existence of an incoming message.

!**   Interface.
!     ----------
!        CALL MPL_PROBE

!        Input required arguments :
!        -------------------------
!           none

!        Input optional arguments :
!        -------------------------
!           KSOURCE  -  rank of process sending the message
!                       default is MPI_ANY_SOURCE
!           KTAG     -  tag of incoming message
!                       default is MPI_ANY_TAG
!           KCOMM    -  Communicator number if different from MPI_COMM_WORLD 
!           LDWAIT   -  = TRUE : waits for a message to be available
!                       = FALSE: return immediately and set
!                                LDFLAG to indicate if a message exists
!           CDSTRING -  Character string for ABORT messages
!                       used when KERROR is not provided

!        Output required arguments :
!        -------------------------
!           none

!        Output optional arguments :
!        -------------------------
!           KERROR   -  return error code.     If not supplied, 
!                       MPL_PROBE aborts when an error is detected.
!           LDFLAG   -  must be supplied if LDWAIT=false
!                       = TRUE if a message exists
!     Author.
!     -------
!        D.Dent, M.Hamrud     ECMWF

!     Modifications.
!     --------------
!        Original: 2000-09-01
!        P. Marguinaud : 01-Jan-2011 : Extends original interface with 
!                                      KCOUNT,KRECVTAG,KFROM (same meaning as
!                                      in all MPL_* routines)

!     ------------------------------------------------------------------

USE EC_PARKIND  ,ONLY : JPIM
USE OML_MOD   ,ONLY : OML_MY_THREAD

USE MPL_MPIF
USE MPL_DATA_MODULE
USE MPL_MESSAGE_MOD

IMPLICIT NONE

PUBLIC MPL_PROBE

PRIVATE

!--- Moved into subroutine to make thrreadsafe----
! INTEGER(KIND=JPIM) :: IRECV_STATUS(MPI_STATUS_SIZE)
! INTEGER(KIND=JPIM) :: ICOMM,ITAG,ISOURCE,IERROR
! LOGICAL :: LLWAIT,LLABORT=.TRUE.

CONTAINS

SUBROUTINE MPL_PROBE(KSOURCE,KTAG,KCOMM,LDWAIT,LDFLAG,CDSTRING,KERROR,KCOUNT,KRECVTAG,KFROM)


#ifdef USE_8_BYTE_WORDS
  USE MPI4TO8, ONLY : &
    MPI_PROBE => MPI_PROBE8, MPI_IPROBE => MPI_IPROBE8
#endif


INTEGER(KIND=JPIM),INTENT(IN), OPTIONAL  :: KSOURCE,KTAG,KCOMM
INTEGER(KIND=JPIM),INTENT(OUT),OPTIONAL  :: KERROR
LOGICAL,INTENT(IN), OPTIONAL  :: LDWAIT
LOGICAL,INTENT(OUT),OPTIONAL  :: LDFLAG
CHARACTER*(*),INTENT(IN),OPTIONAL :: CDSTRING
INTEGER(KIND=JPIM),INTENT(OUT),OPTIONAL  :: KCOUNT, KRECVTAG, KFROM

INTEGER(KIND=JPIM) :: IRECV_STATUS(MPI_STATUS_SIZE)
INTEGER(KIND=JPIM) :: ICOMM,ITAG,ISOURCE,IERROR
LOGICAL :: LLWAIT,LLABORT=.TRUE.
INTEGER(KIND=JPIM) :: ITID
ITID = OML_MY_THREAD()
IF(MPL_NUMPROC < 1) CALL MPL_MESSAGE( &
  & CDMESSAGE='MPL_PROBE: MPL NOT INITIALISED ',LDABORT=LLABORT) 

IF(PRESENT(KCOMM)) THEN
  ICOMM=KCOMM
ELSE
  ICOMM=MPL_COMM_OML(ITID)
ENDIF
IF(PRESENT(KSOURCE)) THEN
  ISOURCE=KSOURCE-1
ELSE
  ISOURCE=MPI_ANY_SOURCE
ENDIF
IF(PRESENT(KTAG)) THEN
  ITAG=KTAG
ELSE
  ITAG=MPI_ANY_TAG
ENDIF

IF(PRESENT(LDWAIT)) THEN
  LLWAIT=LDWAIT
ELSE
  LLWAIT=.TRUE.
ENDIF

IF(LLWAIT) THEN
  CALL MPI_PROBE(ISOURCE,ITAG,ICOMM,IRECV_STATUS,IERROR)
  IF (IERROR == 0) THEN
    IF (PRESENT (KCOUNT))   CALL MPI_GET_COUNT (IRECV_STATUS, MPI_CHARACTER, KCOUNT, IERROR) 
    IF (PRESENT (KRECVTAG)) KRECVTAG = IRECV_STATUS(MPI_TAG)
    IF (PRESENT (KFROM))    KFROM    = IRECV_STATUS(MPI_SOURCE)+1
  ENDIF
ELSE
  IF(PRESENT(LDFLAG)) THEN
    CALL MPI_IPROBE(ISOURCE,ITAG,ICOMM,LDFLAG,IRECV_STATUS,IERROR)
    IF (IERROR == 0 .AND. LDFLAG) THEN
      IF (PRESENT (KCOUNT))   CALL MPI_GET_COUNT (IRECV_STATUS, MPI_CHARACTER, KCOUNT, IERROR)
      IF (PRESENT (KRECVTAG)) KRECVTAG = IRECV_STATUS(MPI_TAG)
      IF (PRESENT (KFROM))    KFROM    = IRECV_STATUS(MPI_SOURCE)+1
    ENDIF
  ELSE
    CALL MPL_MESSAGE(IERROR,'MPL_PROBE: MUST PROVIDE LDFLAG ',CDSTRING, &
                    & LDABORT=LLABORT)
  ENDIF
ENDIF
IF(PRESENT(KERROR)) THEN
  KERROR=IERROR
ELSE
  IF(IERROR /= 0 ) CALL MPL_MESSAGE(IERROR,'MPL_PROBE',CDSTRING,LDABORT=LLABORT)
ENDIF
  
RETURN
END SUBROUTINE MPL_PROBE

END MODULE MPL_PROBE_MOD
