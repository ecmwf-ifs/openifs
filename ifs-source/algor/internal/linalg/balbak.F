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

      SUBROUTINE BALBAK(NM,N,LOW,IGH,SCALE,M,PZ)
!AUTOPROMOTE
      USE PARKIND1, ONLY : JPIM, JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
!
      IMPLICIT NONE
!
      INTEGER(KIND=JPIM) :: NM
      INTEGER(KIND=JPIM) :: N
      INTEGER(KIND=JPIM) :: LOW
      INTEGER(KIND=JPIM) :: IGH
      INTEGER(KIND=JPIM) :: M
      INTEGER(KIND=JPIM) :: I,J,K,II
      REAL(KIND=JPRB) :: SCALE(N)
      REAL(KIND=JPRB) :: PZ(NM,M)
      REAL(KIND=JPRB) :: S
!
!     this subroutine is a translation of the algol procedure balbak,
!     num. math. 13, 293-304(1969) by parlett and reinsch.
!     handbook for auto. comp., vol.ii-linear algebra, 315-326(1971).
!
!     this subroutine forms the eigenvectors of a real general
!     matrix by back transforming those of the corresponding
!     balanced matrix determined by  balanc.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        low and igh are integers determined by  balanc.
!
!        scale contains information determining the permutations
!          and scaling factors used by  balanc.
!
!        m is the number of columns of Pz to be back transformed.
!
!        Pz contains the real and imaginary parts of the eigen-
!          vectors to be back transformed in its first m columns.
!
!     on output
!
!        Pz contains the real and imaginary parts of the
!          transformed eigenvectors in its first m columns.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('BALBAK',0,ZHOOK_HANDLE)
      IF (M .EQ. 0) GO TO 200
      IF (IGH .EQ. LOW) GO TO 120
!
      DO 110 I = LOW, IGH
         S = SCALE(I)
!     .......... left hand eigenvectors are back transformed
!                if the foregoing statement is replaced by
!                s=1.0d0/scale(i). ..........
         DO 100 J = 1, M
  100    PZ(I,J) = PZ(I,J) * S
!
  110 CONTINUE
!     ......... for i=low-1 step -1 until 1,
!               igh+1 step 1 until n do -- ..........
  120 DO 140 II = 1, N
         I = II
         IF (I .GE. LOW .AND. I .LE. IGH) GO TO 140
         IF (I .LT. LOW) I = LOW - II
         K = SCALE(I)
         IF (K .EQ. I) GO TO 140
!
         DO 130 J = 1, M
            S = PZ(I,J)
            PZ(I,J) = PZ(K,J)
            PZ(K,J) = S
  130    CONTINUE
!
  140 CONTINUE
!
  200 CONTINUE
      IF (LHOOK) CALL DR_HOOK('BALBAK',1,ZHOOK_HANDLE)
      RETURN
      ENDSUBROUTINE BALBAK
