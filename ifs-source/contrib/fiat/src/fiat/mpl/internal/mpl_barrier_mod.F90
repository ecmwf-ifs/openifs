! (C) Copyright 2005- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE MPL_BARRIER_MOD

!**** MPL_BARRIER - Barrier synchronisation

!     Purpose.
!     --------
!     Blocks the caller until all group members have called it.

!**   Interface.
!     ----------
!        CALL MPL_BARRIER

!        Input required arguments :
!        -------------------------
!           none

!        Input optional arguments :
!        -------------------------
!           KCOMM    -  Communicator number if different from MPI_COMM_WORLD 
!                       or from that established as the default 
!                       by an MPL communicator routine
!           CDSTRING -  Character string for ABORT messages
!                       used when KERROR is not provided

!        Output required arguments :
!        -------------------------
!           none

!        Output optional arguments :
!        -------------------------
!           KERROR   -  return error code.     If not supplied, 
!                       MPL_BARRIER aborts when an error is detected.
!     Author.
!     -------
!        D.Dent, M.Hamrud     ECMWF

!     Modifications.
!     --------------
!        Original: 2000-09-01
!        Threadsafe: 2004-12-15, J.Hague

!     ------------------------------------------------------------------

USE EC_PARKIND  ,ONLY : JPIM
USE OML_MOD   ,ONLY : OML_MY_THREAD
USE MPL_DATA_MODULE
USE MPL_MESSAGE_MOD

IMPLICIT NONE


PRIVATE

LOGICAL :: LLABORT=.TRUE.

PUBLIC MPL_BARRIER

CONTAINS

SUBROUTINE MPL_BARRIER(KCOMM,CDSTRING,KERROR)


#ifdef USE_8_BYTE_WORDS
  USE MPI4TO8, ONLY : &
    MPI_BARRIER => MPI_BARRIER8
#endif


INTEGER(KIND=JPIM),INTENT(IN),OPTIONAL   :: KCOMM
INTEGER(KIND=JPIM),INTENT(OUT),OPTIONAL  :: KERROR
CHARACTER*(*),INTENT(IN),OPTIONAL :: CDSTRING
INTEGER :: ICOMM,IERROR,ITID
IERROR = 0
ITID = OML_MY_THREAD()
IF(MPL_NUMPROC < 1) CALL MPL_MESSAGE(CDSTRING=CDSTRING,&
  & CDMESSAGE='MPL_BARRIER: MPL NOT INITIALISED ',LDABORT=LLABORT)
 
IF(PRESENT(KCOMM)) THEN
  ICOMM=KCOMM
ELSE
  ICOMM=MPL_COMM_OML(ITID)
ENDIF

CALL MPI_BARRIER(ICOMM,IERROR)

IF(PRESENT(KERROR)) THEN
  KERROR=IERROR
ELSE
  IF(IERROR /= 0 ) CALL MPL_MESSAGE(IERROR,'MPL_BARRIER',CDSTRING,LDABORT=LLABORT)
ENDIF
  
RETURN
END SUBROUTINE MPL_BARRIER

END MODULE MPL_BARRIER_MOD
