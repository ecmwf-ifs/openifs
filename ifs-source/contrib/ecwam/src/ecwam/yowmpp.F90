! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      MODULE YOWMPP
 
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      IMPLICIT NONE

!*    ** *MPP* - MESSAGE PASSING PARAMETERS.

      INTEGER(KIND=JWIM) :: IRANK 
      INTEGER(KIND=JWIM) :: NPROC 
      INTEGER(KIND=JWIM) :: NPREVIOUS 
      INTEGER(KIND=JWIM) :: NNEXT 
      INTEGER(KIND=JWIM) :: NINF 
      INTEGER(KIND=JWIM) :: NSUP 
      INTEGER(KIND=JWIM) :: MPMAXLENGTH 
      INTEGER(KIND=JWIM) :: KTAG 
      INTEGER(KIND=JWIM) :: NPRECR 
      INTEGER(KIND=JWIM) :: NPRECI 

!*     VARIABLE.     TYPE.   PURPOSE.
!      ---------     -----   --------
!      *IRANK*       INTEGER RANK OF THE CURRENT PROCESS (IN MPP MODE) 
!      *NPROC*       INTEGER TOTAL NUMBER OF PROCESSOR USED   
!      *NPREVIOUS*   INTEGER RANK OF PREVIOUS PROCESS
!      *NNEXT*       INTEGER RANK OF NEXT PROCESS
!      *NINF*        INTEGER SMALLEST GRID POINT INDEX USED BY THE PE
!      *NSUP*        INTEGER LARGEST GRID POINT INDEX USED BY THE PE
!      *MPMAXLENGTH* INTEGER MAXIMUM LENGTH OF SCALAR BLOCK MESSAGE
!      *KTAG*        INTEGER GLOBAL TAG FOR MESSAGE PASSING
!      *NPRECR*      INTEGER NUMBER OF BYTES TO REPRESENT REAL NUMBERS
!      *NPRECI*      INTEGER NUMBER OF BYTES TO REPRESENT INTEGERS

! ----------------------------------------------------------------------
      END MODULE YOWMPP
