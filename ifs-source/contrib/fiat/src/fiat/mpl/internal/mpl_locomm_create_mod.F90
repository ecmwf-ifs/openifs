! (C) Copyright 2005- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE MPL_LOCOMM_CREATE_MOD

!**** MPL_LOCOMM_CREATE Create a new communicator

!     Purpose.
!     --------
!     Create a new communicator from lowest N tasks in MPI_COMM_WORLD 
!     and set as default
                                                                                
!**   Interface.
!     ----------
!        CALL MPL_LOCOMM_CREATE

!        Input required arguments :
!        -------------------------
!           N     -  Number of tasks in New Communicator

!        Input optional arguments :
!        -------------------------

!        Output required arguments :
!        -------------------------
!           KCOMM -  New Communicator

!        Output optional arguments :
!        -------------------------
!           MPL_LOCOMM_CREATE aborts when an error is detected.
 
!     Author.
!     -------
!        J.Hague        

!     Modifications.
!     --------------
!        Original: 21/07/2003

!     ------------------------------------------------------------------

USE EC_PARKIND  ,ONLY : JPIM

USE MPL_MPIF
USE MPL_DATA_MODULE
USE MPL_MESSAGE_MOD

IMPLICIT NONE

PRIVATE
PUBLIC MPL_LOCOMM_CREATE

CONTAINS 


SUBROUTINE MPL_LOCOMM_CREATE(N,KCOMM)


#ifdef USE_8_BYTE_WORDS
  USE MPI4TO8, ONLY : &
    MPI_COMM_GROUP => MPI_COMM_GROUP8, MPI_GROUP_INCL => MPI_GROUP_INCL8, MPI_COMM_CREATE => MPI_COMM_CREATE8
#endif


INTEGER(KIND=JPIM),INTENT(OUT) :: KCOMM
INTEGER(KIND=JPIM),INTENT(IN)  :: N

INTEGER(KIND=JPIM) :: IRANK(N)
INTEGER(KIND=JPIM) :: J, IER, IGROUP, MPI_GROUP_WORLD
LOGICAL   :: LLABORT=.TRUE.

DO J=1,N
  IRANK(J)=J-1
ENDDO

CALL MPI_COMM_GROUP(MPL_COMM,MPI_GROUP_WORLD,IER)
IF (IER/=0) CALL MPL_MESSAGE(IER,'MPL_LOCOMM_CREATE: MPI_COMM_GROUP',LDABORT=LLABORT)

CALL MPI_GROUP_INCL(MPI_GROUP_WORLD,N,IRANK,IGROUP,IER)
IF (IER/=0) CALL MPL_MESSAGE(IER,'MPL_LOCOMM_CREATE: MPI_GROUP_INCL',LDABORT=LLABORT)

CALL MPI_COMM_CREATE(MPL_COMM,IGROUP,KCOMM,IER)
IF (IER/=0) CALL MPL_MESSAGE(IER,'MPL_LOCOMM_CREATE: MPI_COMM_CREATE',LDABORT=LLABORT)

RETURN

END SUBROUTINE MPL_LOCOMM_CREATE

END MODULE MPL_LOCOMM_CREATE_MOD
