! (C) Copyright 2005- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE MPL_GROUPS

!     Purpose.
!     --------
!       Use MPI groups for easier to read code (and more efficient
!       communications, at least on IBM).

!     Author.
!     -------
!       Y. Tremolet

!     Modifications.
!     --------------
!       Original: 02-03-13
! ------------------------------------------------------------------

! --- *NOT* THREAD SAFE YET ---

USE EC_PARKIND  ,ONLY : JPIM

USE MPL_MPIF
USE MPL_DATA_MODULE
USE MPL_MESSAGE_MOD

IMPLICIT NONE
PRIVATE
PUBLIC MPL_COMM_GRID, MPL_ALL_LEVS_COMM, MPL_ALL_MS_COMM, &
     & MPL_GROUPS_CREATE, MPL_CART_RANK, MPL_CART_COORDS

INTEGER(KIND=JPIM) :: MPL_COMM_GRID, MPL_ALL_LEVS_COMM, MPL_ALL_MS_COMM, &
           & MPL_GP_GRID
LOGICAL,SAVE :: LGROUPSETUP=.FALSE.
CONTAINS
! ------------------------------------------------------------------

SUBROUTINE MPL_GROUPS_CREATE(KPROCW, KPROCV)


#ifdef USE_8_BYTE_WORDS
  USE MPI4TO8, ONLY : &
    MPI_CART_CREATE => MPI_CART_CREATE8, MPI_COMM_GROUP => MPI_COMM_GROUP8, &
    MPI_CART_SUB => MPI_CART_SUB8
#endif


IMPLICIT NONE
INTEGER(KIND=JPIM), INTENT(IN) :: KPROCW, KPROCV

INTEGER(KIND=JPIM) :: IDIMS(2), IERR
LOGICAL :: LTORUS(2), LDIMS(2), LREORDER

IF(LGROUPSETUP) RETURN

IDIMS(1)=KPROCW
IDIMS(2)=KPROCV
LTORUS(1)=.FALSE.
LTORUS(2)=.FALSE.
LREORDER=.FALSE.

CALL MPI_CART_CREATE(MPL_COMM_OML(1), 2, IDIMS, LTORUS, LREORDER, &
                   & MPL_COMM_GRID, IERR)
IF (IERR/=0) CALL MPL_MESSAGE(IERR,'MPL_GROUPS_CREATE: MPI_CART_CREATE')

CALL MPI_COMM_GROUP(MPL_COMM_GRID, MPL_GP_GRID, IERR)
IF (IERR/=0) CALL MPL_MESSAGE(IERR,'MPL_GROUPS_CREATE: mpi_comm_group')

! Group all levels for same Ms
! ----------------------------
LDIMS(1)=.FALSE.
LDIMS(2)=.TRUE.
CALL MPI_CART_SUB(MPL_COMM_GRID, LDIMS, MPL_ALL_LEVS_COMM, IERR)
IF (IERR/=0) CALL MPL_MESSAGE(IERR,'MPL_GROUPS_CREATE: mpi_cart_sub 1')

! Group all Ms for same levels
! ----------------------------
LDIMS(1)=.TRUE.
LDIMS(2)=.FALSE.
CALL MPI_CART_SUB(MPL_COMM_GRID, LDIMS, MPL_ALL_MS_COMM, IERR)
IF (IERR/=0) CALL MPL_MESSAGE(IERR,'MPL_GROUPS_CREATE: mpi_cart_sub 2')

LGROUPSETUP=.TRUE.

END SUBROUTINE MPL_GROUPS_CREATE

! ------------------------------------------------------------------

FUNCTION MPL_CART_RANK(KPROCW, KPROCV)

#ifdef USE_8_BYTE_WORDS
  USE MPI4TO8, ONLY : &
    MPI_CART_RANK => MPI_CART_RANK8
#endif

IMPLICIT NONE
INTEGER(KIND=JPIM), INTENT(IN)  :: KPROCW, KPROCV
INTEGER(KIND=JPIM) :: MPL_CART_RANK

INTEGER(KIND=JPIM) :: IDIMS(2), IPROC, IERR

IDIMS(1)=KPROCW-1
IDIMS(2)=KPROCV-1

CALL MPI_CART_RANK(MPL_COMM_GRID, IDIMS, IPROC, IERR)
IF (IERR/=0) CALL MPL_MESSAGE(IERR,'MPL_CART_RANK: mpi_cart_rank')

MPL_CART_RANK=IPROC+1

END FUNCTION MPL_CART_RANK

! ------------------------------------------------------------------

SUBROUTINE MPL_CART_COORDS(KPROC, KPROCW, KPROCV)

#ifdef USE_8_BYTE_WORDS
  USE MPI4TO8, ONLY : &
    MPI_CART_COORDS => MPI_CART_COORDS8
#endif


IMPLICIT NONE
INTEGER(KIND=JPIM), INTENT(IN)   :: KPROC
INTEGER(KIND=JPIM), INTENT(OUT)  :: KPROCW, KPROCV

INTEGER(KIND=JPIM) :: IDIMS(2), IPROC, IERR

IPROC=KPROC-1

CALL MPI_CART_COORDS(MPL_COMM_GRID, IPROC, 2, IDIMS, IERR)
IF (IERR/=0) CALL MPL_MESSAGE(IERR,'MPL_CART_COORDS: mpi_cart_coords')

KPROCW=IDIMS(1)+1
KPROCV=IDIMS(2)+1

END SUBROUTINE MPL_CART_COORDS

! ------------------------------------------------------------------

END MODULE MPL_GROUPS
