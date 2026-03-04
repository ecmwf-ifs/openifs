! (C) Copyright 2005- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE MPL_BUFFER_METHOD_MOD

!**** MPL_BUFFER_METHOD Establish message passing default method

!     Purpose.
!     --------
!     Setup the message passing buffering 
!     by allocating an attached buffer if required.

!**   Interface.
!     ----------
!        CALL MPL_BUFFER_METHOD

!        Input required arguments :
!        -------------------------
!           KMP_TYPE -  buffering type
!                       possible values are :
!                       JP_BLOCKING_STANDARD, JP_BLOCKING_BUFFERED
!                       defined as parameters in MPL_DATA_MODULE

!        Input optional arguments :
!        -------------------------
!           KMBX_SIZE - Size (in bytes) of attached buffer 
!                       if KMP_TYPE=JP_BLOCKING_BUFFERED
!           KPROCIDS  - array of processor ids
!           LDINFO   -  = .TRUE.  : Print informative msgs from MPL_INIT (default) 
!                       = .FALSE. : Do not print

!        Output required arguments :
!        -------------------------
!           none

!        Output optional arguments :
!        -------------------------
!           KERROR   -  return error code.     If not supplied, 
!                       MPL_BUFFER_METHOD aborts when an error is detected.
!     Author.
!     -------
!        D.Dent, M.Hamrud     ECMWF

!     Modifications.
!     --------------
!        Original: 2000-09-01

!     ------------------------------------------------------------------

USE MPL_DATA_MODULE
USE MPL_MESSAGE_MOD

IMPLICIT NONE

PRIVATE
PUBLIC MPL_BUFFER_METHOD

CONTAINS 

SUBROUTINE MPL_BUFFER_METHOD(KMP_TYPE,KMBX_SIZE,KERROR,KPROCIDS,LDINFO)


#ifdef USE_8_BYTE_WORDS
  USE MPI4TO8, ONLY : &
    MPI_BUFFER_DETACH => MPI_BUFFER_DETACH8, MPI_BUFFER_ATTACH => MPI_BUFFER_ATTACH8
#endif


USE EC_PARKIND  ,ONLY : JPIM

INTEGER(KIND=JPIM),INTENT(IN) ::  KMP_TYPE
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN)  :: KMBX_SIZE
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN)  :: KPROCIDS(:)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(OUT) :: KERROR
LOGICAL,INTENT(IN),OPTIONAL :: LDINFO
INTEGER(KIND=JPIM) :: IMBX_DEFAULT_SIZE = 1000000
INTEGER(KIND=JPIM) :: IBUFFMPI,IERROR,ILEN
LOGICAL :: LLABORT=.TRUE., LLINFO

IF(MPL_NUMPROC < 1) CALL MPL_MESSAGE( &
  & CDMESSAGE='MPL_BUFFER_METHOD: MPL NOT INITIALISED ',LDABORT=LLABORT) 

IF (ALLOCATED(MPL_ATTACHED_BUFFER)) THEN
  CALL MPI_BUFFER_DETACH(MPL_ATTACHED_BUFFER,MPL_MBX_SIZE,IERROR)
  DEALLOCATE(MPL_ATTACHED_BUFFER)
ENDIF

IF(PRESENT(LDINFO)) THEN
  LLINFO = LDINFO
ELSE
  LLINFO = .TRUE.
ENDIF

IF(KMP_TYPE == JP_BLOCKING_STANDARD) THEN
  IBUFFMPI=MPL_MBX_SIZE
ELSE IF(KMP_TYPE == JP_BLOCKING_BUFFERED) THEN
  IBUFFMPI=KMBX_SIZE
  IF(IBUFFMPI == 0) IBUFFMPI=IMBX_DEFAULT_SIZE
!    convert to bytes
  ILEN = (IBUFFMPI-1)/JP_ATTACHED_BUFFER_BYTES+1
  ALLOCATE(MPL_ATTACHED_BUFFER(ILEN))
#ifdef OPS_COMPILE
  IERROR = 0
#else
  CALL MPI_BUFFER_ATTACH(MPL_ATTACHED_BUFFER,IBUFFMPI,IERROR)
#endif
  IF(PRESENT(KERROR)) THEN
    KERROR=IERROR
  ELSE
    IF( IERROR /= 0 )THEN
      CALL MPL_MESSAGE(IERROR,'MPL_BUFFER_METHOD ','MPI_BUFFER_ATTACH ERROR',LDABORT=LLABORT)
    ENDIF
  ENDIF
ELSE
!    invalid type
  IF(PRESENT(KERROR)) THEN
    KERROR=1
  ELSE
    CALL MPL_MESSAGE(KMP_TYPE,'MPL_BUFFER_METHOD','INVALID KMP_TYPE=',LDABORT=LLABORT)
  ENDIF
ENDIF

MPL_MBX_SIZE=IBUFFMPI
MPL_METHOD=KMP_TYPE

IF (MPL_RANK == 1) THEN
  IF (LLINFO) WRITE(MPL_UNIT,'(A,I2,I12)') 'MPL_BUFFER_METHOD: ',MPL_METHOD,MPL_MBX_SIZE
ENDIF

IF(PRESENT(KPROCIDS)) THEN
  IF(SIZE(KPROCIDS) < MPL_NUMPROC) THEN
    CALL MPL_MESSAGE(CDMESSAGE='MPL_BUFFER_METHOD: KPROCIDS NOT CORRECT',LDABORT=LLABORT)
  ELSE
    MPL_IDS=KPROCIDS
  ENDIF
ENDIF

RETURN
END SUBROUTINE MPL_BUFFER_METHOD

END MODULE MPL_BUFFER_METHOD_MOD
