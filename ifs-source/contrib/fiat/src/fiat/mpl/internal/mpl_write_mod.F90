! (C) Copyright 2005- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE MPL_WRITE_MOD

!
!     Purpose.  write to an MPIIO file 
!     --------
!
!
!     Interface.
!     ----------
!        call mpl_write(...)
!
!        Explicit arguments :
!        --------------------
!
!        input arguments:
!        kfptr   - handle for file
!        kop     - requested operation 
!        kbuf    - buffer containing data to be written
!        klen    - length of buffer in words
!        input/output arguements:
!        kreq    - request handle for non-blocking operations
!        output arguments:
!        kerror  - error code
!
!        Implicit arguments :
!        --------------------
!
!     Method.
!     -------
!     MPL supports 4 styles of MPIIO
!
!     kop = 1    -  Blocking, non collective, shared file pointer
!                   using MPI_FILE_WRITE_SHARED,
!                         MPI_FILE_READ_SHARED
!     kop = 2    -  Blocking, collective, ordered, shared file pointer
!                   using MPI_FILE_WRITE_ORDERED,
!                         MPI_FILE_READ_ORDERED
!     kop = 3    -  Non Blocking, non collective, shared file pointer
!                   using MPI_FILE_IWRITE_SHARED,
!                         MPI_FILE_IREAD_SHARED
!                   and MPI_WAIT
!     kop = 4    -  Non Blocking, collective, ordered, shared file pointer
!                   using MPI_FILE_WRITE_ORDERED_BEGIN/END,
!                         MPI_FILE_READ_ORDERED_BEGIN/END
!
!     Externals.
!     ----------
!
!     Reference.
!     ----------
!        none yet
!
!     Author.
!     -------
!        G.Mozdzynski
!
!     Modifications.
!     --------------
!        Original : 2000-12-08 (Based on MPE_WRITE)
!        R. EL Khatib 24-May-2011 Change ifdef MPI2 into ifndef MPI1
!     -----------------------------------------------------------------
!
USE EC_PARKIND  ,ONLY : JPIM ,JPRM

USE MPL_MPIF
USE MPL_DATA_MODULE
USE MPL_MESSAGE_MOD
USE MPL_IOINIT_MOD

IMPLICIT NONE

INTERFACE MPL_WRITE
MODULE PROCEDURE MPL_WRITE_INT,MPL_WRITE_REAL8
END INTERFACE

PRIVATE
PUBLIC MPL_WRITE

CONTAINS

SUBROUTINE MPL_WRITE_INT(KFPTR,KOP,KBUF,KLEN,KREQ,KERROR)


#ifdef USE_8_BYTE_WORDS
  USE MPI4TO8, ONLY : &
    MPI_FILE_WRITE_SHARED => MPI_FILE_WRITE_SHARED8, &
    MPI_FILE_WRITE_ORDERED => MPI_FILE_WRITE_ORDERED8, &
    MPI_FILE_IWRITE_SHARED => MPI_FILE_IWRITE_SHARED8, &
    MPI_FILE_WRITE_ORDERED_BEGIN => MPI_FILE_WRITE_ORDERED_BEGIN8, &
    MPI_WAIT => MPI_WAIT8, &
    MPI_FILE_WRITE_ORDERED_END => MPI_FILE_WRITE_ORDERED_END8
#endif

INTEGER(KIND=JPIM),INTENT(IN) :: KFPTR,KOP,KLEN
INTEGER(KIND=JPIM),INTENT(OUT) :: KERROR
INTEGER(KIND=JPIM) KBUF(:)
INTEGER(KIND=JPIM) KREQ


INTEGER(KIND=JPIM) STATUS(MPI_STATUS_SIZE)
!
#ifndef MPI1

!     -----------------------------------------------------------------
!
!     1.    Preamble 
!           --------

IF( MPL_RANK > MPL_NUMIO ) THEN
  KERROR = -1
  RETURN
ENDIF

IF( KOP >= 1.AND.KOP <= 4 ) THEN

  IF( KOP /= MPL_IOP ) THEN
    KERROR = -1
    RETURN
  ENDIF

ENDIF
!     -----------------------------------------------------------------
!
!     2.    Check style and take appropriate action 
!           ---------------------------------------


IF( KOP == 1 ) THEN

! blocking write, non collective, shared file pointer
        
  CALL MPI_FILE_WRITE_SHARED(KFPTR,&
                           & KBUF,&
                           & KLEN,&
                           & MPI_INTEGER,&
                           & STATUS,&
                           & KERROR)

ELSEIF( KOP == 2 ) THEN

! blocking write, collective, ordered with shared file pointer

  CALL MPI_FILE_WRITE_ORDERED(KFPTR,&
                            & KBUF,&
                            & KLEN,&
                            & MPI_INTEGER,&
                            & STATUS,&
                            & KERROR)

ELSEIF( KOP == 3 ) THEN

! non blocking write, non collective, shared file pointer

  CALL MPI_FILE_IWRITE_SHARED(KFPTR,&
                            & KBUF,&
                            & KLEN,&
                            & MPI_INTEGER,&
                            & KREQ,&
                            & KERROR)

ELSEIF( KOP == 4 ) THEN

! non blocking write, collective, ordered with shared file pointer

  CALL MPI_FILE_WRITE_ORDERED_BEGIN(KFPTR,&
                                  & KBUF,&
                                  & KLEN,&
                                  & MPI_INTEGER,&
                                  & KERROR)

ELSEIF( KOP == 5 ) THEN

  CALL MPI_WAIT(KREQ,&
              & STATUS,&
              & KERROR )

ELSEIF( KOP == 6 ) THEN

  CALL MPI_FILE_WRITE_ORDERED_END(KFPTR,&
                                & KBUF,&
                                & STATUS,&
                                & KERROR)
 
ELSE

  KERROR =-1
  RETURN

ENDIF

#else

CALL ABOR1('MPL_WRITE_INT not built with MPI2')

#endif
!
!     -----------------------------------------------------------------
RETURN
END SUBROUTINE MPL_WRITE_INT

SUBROUTINE MPL_WRITE_REAL8(KFPTR,KOP,PBUF,KLEN,KREQ,KERROR)


#ifdef USE_8_BYTE_WORDS
  USE MPI4TO8, ONLY : &
    MPI_FILE_WRITE_SHARED => MPI_FILE_WRITE_SHARED8, &
    MPI_FILE_WRITE_ORDERED => MPI_FILE_WRITE_ORDERED8, &
    MPI_FILE_IWRITE_SHARED => MPI_FILE_IWRITE_SHARED8, &
    MPI_FILE_WRITE_ORDERED_BEGIN => MPI_FILE_WRITE_ORDERED_BEGIN8, &
    MPI_WAIT => MPI_WAIT8, &
    MPI_FILE_WRITE_ORDERED_END => MPI_FILE_WRITE_ORDERED_END8
#endif

INTEGER(KIND=JPIM),INTENT(IN) :: KFPTR,KOP,KLEN
INTEGER(KIND=JPIM),INTENT(OUT) :: KERROR
REAL(KIND=JPRM) PBUF(:)
INTEGER(KIND=JPIM) KREQ


INTEGER(KIND=JPIM) STATUS(MPI_STATUS_SIZE)
!
#ifndef MPI1

!     -----------------------------------------------------------------
!
!     1.    Preamble 
!           --------

IF( MPL_RANK > MPL_NUMIO ) THEN
  KERROR = -1
  RETURN
ENDIF

IF( KOP >= 1.AND.KOP <= 4 ) THEN

  IF( KOP /= MPL_IOP ) THEN
    KERROR = -1
    RETURN
  ENDIF

ENDIF
!     -----------------------------------------------------------------
!
!     2.    Check style and take appropriate action 
!           ---------------------------------------


IF( KOP == 1 ) THEN

! blocking write, non collective, shared file pointer
        
  CALL MPI_FILE_WRITE_SHARED(KFPTR,&
                           & PBUF,&
                           & KLEN,&
                           & MPI_REAL8,&
                           & STATUS,&
                           & KERROR)

ELSEIF( KOP == 2 ) THEN

! blocking write, collective, ordered with shared file pointer

  CALL MPI_FILE_WRITE_ORDERED(KFPTR,&
                            & PBUF,&
                            & KLEN,&
                            & MPI_REAL8,&
                            & STATUS,&
                            & KERROR)

ELSEIF( KOP == 3 ) THEN

! non blocking write, non collective, shared file pointer

  CALL MPI_FILE_IWRITE_SHARED(KFPTR,&
                            & PBUF,&
                            & KLEN,&
                            & MPI_REAL8,&
                            & KREQ,&
                            & KERROR)

ELSEIF( KOP == 4 ) THEN

! non blocking write, collective, ordered with shared file pointer

  CALL MPI_FILE_WRITE_ORDERED_BEGIN(KFPTR,&
                                  & PBUF,&
                                  & KLEN,&
                                  & MPI_REAL8,&
                                  & KERROR)

ELSEIF( KOP == 5 ) THEN

  CALL MPI_WAIT(KREQ,&
              & STATUS,&
              & KERROR )

ELSEIF( KOP == 6 ) THEN

  CALL MPI_FILE_WRITE_ORDERED_END(KFPTR,&
                                & PBUF,&
                                & STATUS,&
                                & KERROR)
 
ELSE

  KERROR =-1
  RETURN

ENDIF

#else

CALL ABOR1('MPL_WRITE_REAL8 not built with MPI2')

#endif

!
!     -----------------------------------------------------------------
RETURN
END SUBROUTINE MPL_WRITE_REAL8

END MODULE MPL_WRITE_MOD
