! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      MODULE YOWINDN

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      IMPLICIT NONE

!*    ** *INDNL* - INDICES AND WEIGHTS USED IN THE COMPUTATION
!                  OF THE NONLINEAR TRANSFER RATE.

      INTEGER(KIND=JWIM), PARAMETER  :: NINL = 5
      INTEGER(KIND=JWIM), PARAMETER  :: NRNL = 25
      INTEGER(KIND=JWIM)             :: KFRH 
      INTEGER(KIND=JWIM)             :: MFRSTLW
      INTEGER(KIND=JWIM)             :: MLSTHG

      INTEGER(KIND=JWIM), ALLOCATABLE :: IKP(:) 
      INTEGER(KIND=JWIM), ALLOCATABLE :: IKP1(:) 
      INTEGER(KIND=JWIM), ALLOCATABLE :: IKM(:) 
      INTEGER(KIND=JWIM), ALLOCATABLE :: IKM1(:) 
      INTEGER(KIND=JWIM), ALLOCATABLE :: K1W(:,:) 
      INTEGER(KIND=JWIM), ALLOCATABLE :: K2W(:,:) 
      INTEGER(KIND=JWIM), ALLOCATABLE :: K11W(:,:) 
      INTEGER(KIND=JWIM), ALLOCATABLE :: K21W(:,:) 
      INTEGER(KIND=JWIM), ALLOCATABLE :: INLCOEF(:,:)

      REAL(KIND=JWRB), ALLOCATABLE :: AF11(:) 
      REAL(KIND=JWRB), ALLOCATABLE :: FKLAP(:) 
      REAL(KIND=JWRB), ALLOCATABLE :: FKLAP1(:) 
      REAL(KIND=JWRB), ALLOCATABLE :: FKLAM(:) 
      REAL(KIND=JWRB), ALLOCATABLE :: FKLAM1(:) 
      REAL(KIND=JWRB)              :: ACL1 
      REAL(KIND=JWRB)              :: ACL2 
      REAL(KIND=JWRB)              :: CL11 
      REAL(KIND=JWRB)              :: CL21 
      REAL(KIND=JWRB)              :: DAL1 
      REAL(KIND=JWRB)              :: DAL2 
      REAL(KIND=JWRB), ALLOCATABLE :: FRH(:)
      REAL(KIND=JWRB), ALLOCATABLE :: FTRF(:)
      REAL(KIND=JWRB), ALLOCATABLE :: RNLCOEF(:,:)

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   -------
!      *NINL*      INTEGER   SIZE OF INLCOEF
!      *NRNL*      INTEGER   SIZE OF RNLCOEF
!      *KFRH*      INTEGER   SIZE OF FRH
!      *MFRSTLW*   INTEGER   INDEX OF FIRST EXTRA LOW FREQUENCY FOR SNL
!      *MLSTHG*    INTEGER   INDEX OF LAST EXTRA HIGH FREQUENCY FOR SNL
!      *IKP*       INTEGER   FREQUENCY INDEX ARRAY FOR STORING ENERGY
!                            TRANSFER INCREMENTS INTO BINS, WAVE NO. 3.
!      *IKP1*      INTEGER   IKP+1.
!      *IKM*       INTEGER   FREQUENCY INDEX ARRAY FOR STORING ENERGY
!                            TRANSFER INCREMENTS INTO BINS, WAVE NO. 4.
!      *IKM1*      INTEGER   IKM+1
!      *K1W*       INTEGER   ANGULAR INDEX ARRAY FOR STORING ENERGY
!                            TRANSFER INCREMENTS INTO BINS, WAVE NO. 3.
!      *K11W*      INTEGER   K1W(.,1)-1, K1W(.,2)+1.
!      *K2W*       INTEGER   ANGULAR INDEX ARRAY FOR STORING ENERGY
!                            TRANSFER INCREMENTS INTO BINS, WAVE NO. 4.
!      *K21W*      INTEGER   K2W(.,1)+1, K2W(.,2)-1.
!      *INLCOEF*   INTEGER   ARRAY USED TO STORE ALL FREQUENCY DEPENDENT
!                            INDICES FOUND IN SNONLIN
!      *AF11*      REAL      WEIGHTS FOR DISCRETE APPROXIMATION OF NONL
!                            TRANSFER (AT PRESENT ONE TERM ONLY SET TO
!                            3000). MULTIPLIED BY FREQUENCIES **11.
!      *FKLAP*     REAL      WEIGHT IN FREQUENCY GRID FOR INTERPOLATION,
!                            WAVE NO. 3 ("1+LAMBDA" TERM).
!      *FKLAP1*    REAL      1-FKLAP.
!      *FKLAM*     REAL      WEIGHT IN FREQUENCY GRID FOR INTERPOLATION,
!                            WAVE NO. 4 ("1-LAMBDA" TERM).
!      *ACL1*      REAL      WEIGHT IN ANGULAR GRID FOR INTERPOLATION,
!                            WAVE NO. 3 ("1+LAMBDA" TERM).
!      *ACL2*      REAL      WEIGHT IN ANGULAR GRID FOR INTERPOLATION,
!                            WAVE NO. 4 ("1-LAMBDA" TERM).
!      *CL11*      REAL      1.-ACL1.
!      *CL21*      REAL      1.-ACL2.
!      *DAL1*      REAL      1./ACL1.
!      *DAL2*      REAL      1./ACL2.
!      *FRH*       REAL      TAIL FREQUENCY RATION **5
!      *FTRF*      REAL      FRONT TAIL REDUCTIOn FACTOR USED TO A SPECTRAL
!                            TAIL IN FRONT OF THE FIRST DISCRETISED FREQUENCY
!      *RNLCOEF*   REAL      ARRAY USED TO STORE ALL FREQUENCY DEPENDENT
!                            COEFFICIENT FOUND IN SNONLIN

! ----------------------------------------------------------------------
      END MODULE YOWINDN
