! (C) Copyright 2005- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE MPL_MESSAGE_MOD

!**** MPL_MESSAGE - Prints message

!     Purpose.
!     --------
!     Creates an ASCII message for printing and optionally aborts

!**   Interface.
!     ----------
!        CALL MPL_MESSAGE

!        Input required arguments :
!        -------------------------
!           CDMESSAGE-  character string for message

!        Input optional arguments :
!        -------------------------
!           KERROR   -  Error number
!           CDSTRING -  Optional additional message 
!                       prepended to CDMESSAGE
!           LDABORT  -  forces ABORT if true

!        Output required arguments :
!        -------------------------
!           none

!        Output optional arguments :
!        -------------------------
!           none
!     Author.
!     -------
!        D.Dent, M.Hamrud     ECMWF

!     Modifications.
!     --------------
!        Original: 2000-09-01

!     ------------------------------------------------------------------

USE MPL_MPIF
USE MPL_DATA_MODULE
USE MPL_ABORT_MOD

USE EC_PARKIND  ,ONLY : JPIM

PRIVATE

PUBLIC MPL_MESSAGE

CONTAINS 

SUBROUTINE MPL_MESSAGE(KERROR,CDMESSAGE,CDSTRING,LDABORT)


#ifdef USE_8_BYTE_WORDS
  USE MPI4TO8, ONLY : &
    MPI_ERROR_STRING => MPI_ERROR_STRING8
#endif


IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN),OPTIONAL     :: KERROR
CHARACTER*(*),INTENT(IN)          :: CDMESSAGE
CHARACTER*(*),INTENT(IN),OPTIONAL :: CDSTRING
LOGICAL      ,INTENT(IN),OPTIONAL :: LDABORT

CHARACTER*(MPI_MAX_ERROR_STRING)  :: CLMPI_ERROR
CHARACTER*12                      :: CLERROR
INTEGER(KIND=JPIM)                :: IRESULTLEN,IERROR

IF(PRESENT(KERROR)) THEN
  WRITE(CLERROR,'(I0)') KERROR
ELSE
  CLERROR=' '
ENDIF
IF(PRESENT(CDSTRING)) THEN
  WRITE(MPL_UNIT,'(4(1X,A),I0)') CDSTRING,CDMESSAGE,CLERROR, &
                               & ' FROM PROCESSOR ',MPL_RANK
ELSE
  WRITE(MPL_UNIT,'(3(1X,A),I0)') CDMESSAGE,CLERROR, &
                               & ' FROM PROCESSOR ',MPL_RANK
ENDIF

IF(PRESENT(KERROR)) THEN
  CALL MPI_ERROR_STRING(KERROR,CLMPI_ERROR,IRESULTLEN,IERROR)
  WRITE(MPL_UNIT,'(1X,2A,I0)') CLMPI_ERROR(1:IRESULTLEN),' in processor ',MPL_RANK
ELSE
  CLMPI_ERROR=' '
  IRESULTLEN=1
ENDIF

IF(PRESENT(LDABORT)) THEN
  IF(LDABORT) THEN
    WRITE(0,'(1X,2A,I0)') CLMPI_ERROR(1:IRESULTLEN),' in processor ',MPL_RANK
    CALL MPL_ABORT('ABORT')
  ENDIF
ENDIF
RETURN
END SUBROUTINE MPL_MESSAGE

END MODULE MPL_MESSAGE_MOD
