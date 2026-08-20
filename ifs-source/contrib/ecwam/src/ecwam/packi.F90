! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE PACKI (NIN, NOUT, NDIM, IFLAG, IPACK)

! ----------------------------------------------------------------------

!**** *PACKI* - PACKS A INTEGER ARRAY.

!     R. PORTZ     MPI         15/01/1991

!*    PURPOSE.
!     --------

!       TO REMOVE FLAGGED POINTS FROM AN INTEGER ARRAY.

!**   INTERFACE.
!     ----------

!       *CALL* *PACKI (NIN, NOUT, NDIM, IFLAG, IPACK)*
!          *NIN*   -  NUMBER OP POINTS IN ARRAYS BEFORE PACKING.
!          *NOUT*  -  NUMBER OF POINTS IN ARRAYES AFTER PACKING.
!          *NDIM*  -  DIMENSION OF ARRAYS.
!          *IFLAG* -  FLAG ARRAY.
!          *IPACK* -  ARRAY TO BE PACKED / PACKED ARRAY AT OUTPUT.

!     METHOD.
!     -------

!       NONE.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: NIN, NDIM
      INTEGER(KIND=JWIM), INTENT(OUT) :: NOUT
      INTEGER(KIND=JWIM), DIMENSION(NDIM), INTENT(INOUT) :: IFLAG, IPACK


      INTEGER(KIND=JWIM) :: K, I

! ----------------------------------------------------------------------


!*    1. LOOP OVER POINTS.
!        -----------------

      K = 0
      DO I = 1, NIN
        IF (IFLAG(I).EQ.0) THEN
          K = K + 1
        ELSE
          IPACK(I-K) = IPACK(I)
        ENDIF
      ENDDO

!*    2. NUMBER OF POINTS LEFT.
!        ----------------------

      NOUT = NIN - K

      END SUBROUTINE PACKI
