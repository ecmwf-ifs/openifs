! (C) Copyright 2005- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE MPL_CLOSE_MOD
!
!     Purpose.  close an MPIIO file 
!     --------
!
!
!     Interface.
!     ----------
!        call mpl_close(...)
!
!        Explicit arguments :
!        --------------------
!
!        input arguments:
!        kfptr   - handle for file pointer
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
!        Original : 2000-12-08 (Based on MPE_CLOSE)
!        R. EL Khatib 24-May-2011 Change ifdef MPI2 into ifndef MPI1
!     -----------------------------------------------------------------
!
USE EC_PARKIND  ,ONLY : JPIM

USE MPL_DATA_MODULE
USE MPL_MESSAGE_MOD
USE MPL_IOINIT_MOD

IMPLICIT NONE

PRIVATE
PUBLIC MPL_CLOSE

CONTAINS

SUBROUTINE MPL_CLOSE(KFPTR,KERROR)


#ifdef USE_8_BYTE_WORDS
  USE MPI4TO8, ONLY : &
    MPI_FILE_CLOSE => MPI_FILE_CLOSE8
#endif

INTEGER(KIND=JPIM),INTENT(IN) :: KFPTR
INTEGER(KIND=JPIM),INTENT(OUT) :: KERROR
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
!
!     -----------------------------------------------------------------
!
!     1.    Close the File
!           --------------
CALL MPI_FILE_CLOSE(KFPTR,KERROR)
!
!     -----------------------------------------------------------------
#else

CALL ABOR1('MPI_CLOSE not built with MPI2')

#endif
!
RETURN
END SUBROUTINE MPL_CLOSE

END MODULE MPL_CLOSE_MOD
