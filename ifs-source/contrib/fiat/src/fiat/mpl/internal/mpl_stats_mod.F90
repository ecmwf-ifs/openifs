! (C) Copyright 2005- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE MPL_STATS_MOD

PRIVATE
PUBLIC :: MPL_STATSINIT, MPL_STATSON, MPL_STATSREAD, MPL_SENDSTATS, MPL_RECVSTATS

CONTAINS

SUBROUTINE MPL_STATSINIT

!**** MPL_STATSINIT - Initialise collection of mpl statistics

!     Purpose.
!     --------
!     Initialises the mpl_stats package 

!**   Interface.
!     ----------
!        CALL MPL_STATSINIT

!        Input required arguments :
!        -------------------------
!           none

!        Input optional arguments :
!        -------------------------
!           none

!        Output required arguments :
!        -------------------------
!           none

!        Output optional arguments :
!        -------------------------
!           none 
 
!     Author.
!     -------
!        P.Towers             ECMWF

!     Modifications.
!     --------------
!        Original: 2011-04-06
!      F. Vana  05-Mar-2015  Support for single precision

!     ------------------------------------------------------------------

USE EC_PARKIND  ,ONLY : JPIM

USE YOMMPLSTATS

IMPLICIT NONE

INTEGER(KIND=JPIM) :: ITHR,OMP_GET_MAX_THREADS

LMPLSTATS=.TRUE.

ITHR = 1
!$ ITHR = OMP_GET_MAX_THREADS()
ITHR = ITHR-1

ALLOCATE(MPLSENDBYTES(0:ITHR))
ALLOCATE(MPLRECVBYTES(0:ITHR))
ALLOCATE(MPLSENDNUM(0:ITHR))
ALLOCATE(MPLRECVNUM(0:ITHR))
  
MPLSENDBYTES(:) = 0
MPLRECVBYTES(:) = 0
MPLSENDNUM(:) = 0
MPLRECVNUM(:) = 0

RETURN
END SUBROUTINE MPL_STATSINIT

SUBROUTINE MPL_STATSON(SENDNUM,SENDBYTES,RECVNUM,RECVBYTES)

!**** MPL_STATSON - Reset mpl statistics counters

!     Purpose.
!     --------
!     Returns the mpl statistics counter values
!     and sets them back to zero
!     non zero returned values correspond to messages that have
!     been sent/received outside of a GSTATS MPL region

!**   Interface.
!     ----------
!        CALL MPL_STATSON(SENDNUM,SENDBYTES,RECVNUM,RECVBYTES)

!        Input required arguments :
!        -------------------------
!           none

!        Input optional arguments :
!        -------------------------
!           none

!        Output required arguments :
!        -------------------------
!           SENDNUM     - number of unknown messages sent
!           SENDBYTES   - number of unknown bytes sent
!           RECVNUM     - number of unknown messages received
!           RECVBYTES   - number of unknown bytes received

!        Output optional arguments :
!        -------------------------
!           none 
 
!     Author.
!     -------
!        P.Towers             ECMWF

!     Modifications.
!     --------------
!        Original: 2011-04-06

!     ------------------------------------------------------------------

USE EC_PARKIND  ,ONLY : JPRD, JPIM

USE YOMMPLSTATS

IMPLICIT NONE

REAL(KIND=JPRD), INTENT(OUT)    :: SENDBYTES,RECVBYTES
INTEGER(KIND=JPIM), INTENT(OUT) :: SENDNUM,RECVNUM


 SENDBYTES = SUM(MPLSENDBYTES(:))
 RECVBYTES = SUM(MPLRECVBYTES(:))
 SENDNUM   = SUM(MPLSENDNUM(:))
 RECVNUM   = SUM(MPLRECVNUM(:))
 
 MPLSENDBYTES(:)=0.0_JPRD
 MPLRECVBYTES(:)=0.0_JPRD
 MPLSENDNUM(:)=0
 MPLRECVNUM(:)=0
 
RETURN
END SUBROUTINE MPL_STATSON

SUBROUTINE MPL_STATSREAD(SENDNUM,SENDBYTES,RECVNUM,RECVBYTES)

!**** MPL_STATSREAD - read mpl statistics counters

!     Purpose.
!     --------
!     returns the mpl statistics counter values
!     before setting them back to zero

!**   Interface.
!     ----------
!        CALL MPL_STATSREAD(SENDNUM,SENDBYTES,RECVNUM,RECVBYTES)

!        Input required arguments :
!        -------------------------
!           none

!        Input optional arguments :
!        -------------------------
!           none

!        Output required arguments :
!        -------------------------
!           SENDNUM     - number of messages sent
!           SENDBYTES   - number of bytes sent
!           RECVNUM     - number of messages received
!           RECVBYTES   - number of bytes received

!        Output optional arguments :
!        -------------------------
!           none 
 
!     Author.
!     -------
!        P.Towers             ECMWF

!     Modifications.
!     --------------
!        Original: 2011-04-06

!     ------------------------------------------------------------------

USE EC_PARKIND  ,ONLY : JPRD, JPIM

USE YOMMPLSTATS

IMPLICIT NONE

REAL(KIND=JPRD), INTENT(OUT)    :: SENDBYTES,RECVBYTES
INTEGER(KIND=JPIM), INTENT(OUT) :: SENDNUM,RECVNUM

  SENDNUM=SUM(MPLSENDNUM(:))
  RECVNUM=SUM(MPLRECVNUM(:))
  SENDBYTES=SUM(MPLSENDBYTES(:))
  RECVBYTES=SUM(MPLRECVBYTES (:))

  MPLSENDNUM(:)=0
  MPLRECVNUM(:)=0
  MPLSENDBYTES(:)=0.0_JPRD
  MPLRECVBYTES(:)=0.0_JPRD

RETURN
END SUBROUTINE MPL_STATSREAD

SUBROUTINE MPL_SENDSTATS(ICOUNT,ITYPE)

!**** MPL_SENDSTATS - collect mpl send statistics

!     Purpose.
!     --------
!     counts the number of messages and volume sent

!**   Interface.
!     ----------
!        CALL MPL_SENDSTATS(ICOUNT,ITYPE)

!        Input required arguments :
!        -------------------------
!           ICOUNT   - The number of elements sent
!           ITYPE    - The type of an element 

!        Input optional arguments :
!        -------------------------
!           none

!        Output required arguments :
!        -------------------------
!           none

!        Output optional arguments :
!        -------------------------
!           none 
 
!     Author.
!     -------
!        P.Towers             ECMWF

!     Modifications.
!     --------------
!        Original: 2011-04-06

!     ------------------------------------------------------------------

USE EC_PARKIND  ,ONLY : JPIM

USE MPL_DATA_MODULE
USE MPL_MESSAGE_MOD
USE YOMMPLSTATS

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)            :: ICOUNT
INTEGER(KIND=JPIM),INTENT(IN)            :: ITYPE

INTEGER(KIND=JPIM)  ISIZE,IERR,ITH,OMP_GET_THREAD_NUM

ITH = 0
!$ ITH = OMP_GET_THREAD_NUM()

MPLSENDNUM(ITH) = MPLSENDNUM(ITH) + 1

CALL MPI_TYPE_SIZE(ITYPE,ISIZE,IERR)

MPLSENDBYTES(ITH)=MPLSENDBYTES(ITH) + FLOAT(ISIZE * ICOUNT)

  
RETURN
END SUBROUTINE MPL_SENDSTATS

SUBROUTINE MPL_RECVSTATS(ICOUNT,ITYPE)

!**** MPL_RECVSTATS - collect mpl recv statistics

!     Purpose.
!     --------
!     counts the number of messages and volume received

!**   Interface.
!     ----------
!        CALL MPL_RECVSTATS(ICOUNT,ITYPE)

!        Input required arguments :
!        -------------------------
!           ICOUNT   - The number of elements received
!           ITYPE    - The type of an element 

!        Input optional arguments :
!        -------------------------
!           none

!        Output required arguments :
!        -------------------------
!           none

!        Output optional arguments :
!        -------------------------
!           none 
 
!     Author.
!     -------
!        P.Towers             ECMWF

!     Modifications.
!     --------------
!        Original: 2011-04-06

!     ------------------------------------------------------------------

USE EC_PARKIND  ,ONLY : JPIM

USE MPL_DATA_MODULE
USE MPL_MESSAGE_MOD
USE YOMMPLSTATS

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)            :: ICOUNT
INTEGER(KIND=JPIM),INTENT(IN)            :: ITYPE

INTEGER(KIND=JPIM) ISIZE,IERR,ITH,OMP_GET_THREAD_NUM

ITH = 0
!$ ITH = OMP_GET_THREAD_NUM()

MPLRECVNUM(ITH) = MPLRECVNUM(ITH) + 1

CALL MPI_TYPE_SIZE(ITYPE,ISIZE,IERR)

MPLRECVBYTES(ITH)=MPLRECVBYTES(ITH) + FLOAT(ISIZE * ICOUNT)
  
RETURN
END SUBROUTINE MPL_RECVSTATS

END MODULE MPL_STATS_MOD
