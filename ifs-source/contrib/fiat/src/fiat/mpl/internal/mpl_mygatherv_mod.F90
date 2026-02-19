! (C) Copyright 2005- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE MPL_MYGATHERV_MOD

USE EC_PARKIND  ,ONLY : JPRD, JPIM

USE MPL_MPIF
USE MPL_DATA_MODULE
USE MPL_MESSAGE_MOD

IMPLICIT NONE
PRIVATE
PUBLIC MPL_MYGATHERV

LOGICAL :: LLABORT=.TRUE.

CONTAINS

! ------------------------------------------------------------------
SUBROUTINE MPL_MYGATHERV(PSEND,KSEND,PRECV,KRECV,KDISPL,KROOT,KCOMM)


#ifdef USE_8_BYTE_WORDS
  USE MPI4TO8, ONLY : &
    MPI_GATHERV => MPI_GATHERV8
#endif


REAL(KIND=JPRD), INTENT(IN)  :: PSEND(:)
REAL(KIND=JPRD), INTENT(OUT) :: PRECV(:)
INTEGER(KIND=JPIM), INTENT(IN) :: KSEND, KRECV(:), KDISPL(:)
INTEGER(KIND=JPIM), INTENT(IN) :: KROOT, KCOMM
INTEGER(KIND=JPIM) :: IERR

CALL MPI_GATHERV(PSEND,KSEND,INT(MPI_REAL8), &
               & PRECV,KRECV,KDISPL,INT(MPI_REAL8),KROOT-1,KCOMM,IERR)

IF (IERR/=0) CALL MPL_MESSAGE(IERR,'MPL_MYGATHERV',LDABORT=LLABORT)

END SUBROUTINE MPL_MYGATHERV
! ------------------------------------------------------------------

END MODULE MPL_MYGATHERV_MOD
