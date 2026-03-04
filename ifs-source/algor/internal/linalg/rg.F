! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction
! 
! (C) Copyright 1989- Meteo-France.
! 

      SUBROUTINE RG(NM,N,PA,PWR,PWI,MATZ,PZ,KIV1,PFV1,IERR)
!AUTOPROMOTE
      USE PARKIND1, ONLY : JPIM, JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
!
      IMPLICIT NONE
!
      INTEGER(KIND=JPIM) :: N
      INTEGER(KIND=JPIM) :: NM
      INTEGER(KIND=JPIM) :: IS1,IS2
      INTEGER(KIND=JPIM) :: IERR
      INTEGER(KIND=JPIM) :: MATZ
      REAL(KIND=JPRB) :: PA(NM,N)
      REAL(KIND=JPRB) :: PWR(N)
      REAL(KIND=JPRB) :: PWI(N)
      REAL(KIND=JPRB) :: PZ(NM,N)
      REAL(KIND=JPRB) :: PFV1(N)
      INTEGER(KIND=JPIM) :: KIV1(N)
!
!     this subroutine calls the recommended sequence of
!     subroutines from the eigensystem subroutine package (eispack)
!     to find the eigenvalues and eigenvectors (if desired)
!     of a real general matrix.
!
!     on input
!
!        nm  must be set to the row dimension of the two-dimensional
!        array parameters as declared in the calling program
!        dimension statement.
!
!        n  is the order of the matrix pa.
!
!        pa contains the real general matrix.
!
!        matz  is an integer variable set equal to zero if
!        only eigenvalues are desired.  otherwise it is set to
!        any non-zero integer for both eigenvalues and eigenvectors.
!
!     on output
!
!        pwr and  pwi contain the real and imaginary parts,
!        respectively, of the eigenvalues.  complex conjugate
!        pairs of eigenvalues appear consecutively with the
!        eigenvalue having the positive imaginary part first.
!
!        pz contains the real and imaginary parts of the eigenvectors
!        if matz is not zero.  if the j-th eigenvalue is real, the
!        j-th column of pz  contains its eigenvector.  if the j-th
!        eigenvalue is complex with positive imaginary part, the
!        j-th and (j+1)-th columns of pz  contain the real and
!        imaginary parts of its eigenvector.  the conjugate of this
!        vector is the eigenvector for the conjugate eigenvalue.
!
!        ierr  is an integer output variable set equal to an error
!           completion code described in the documentation for hqr
!           and hqr2.  the normal completion code is zero.
!
!        kiv1  and  pfv1  are temporary storage arrays.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('RG',0,ZHOOK_HANDLE)
      IF (N .LE. NM) GO TO 10
      IERR = 10 * N
      GO TO 50
!
   10 CALL  BALANC(NM,N,PA,IS1,IS2,PFV1)
      CALL  ELMHES(NM,N,IS1,IS2,PA,KIV1)
      IF (MATZ .NE. 0) GO TO 20
!     .......... find eigenvalues only ..........
      CALL  HQR(NM,N,IS1,IS2,PA,PWR,PWI,IERR)
      GO TO 50
!     .......... find both eigenvalues and eigenvectors ..........
   20 CALL  ELTRAN(NM,N,IS1,IS2,PA,KIV1,PZ)
      CALL  HQR2(NM,N,IS1,IS2,PA,PWR,PWI,PZ,IERR)
      IF (IERR .NE. 0) GO TO 50
      CALL  BALBAK(NM,N,IS1,IS2,PFV1,N,PZ)
   50 CONTINUE
      IF (LHOOK) CALL DR_HOOK('RG',1,ZHOOK_HANDLE)
      RETURN
      ENDSUBROUTINE RG
