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

      SUBROUTINE ELMHES(NM,N,LOW,IGH,A,KINT)
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
      INTEGER(KIND=JPIM) :: I,J,M,LA,KP1,MM1,MP1
      REAL(KIND=JPRB) :: A(NM,N)
      REAL(KIND=JPRB) :: X,Y
      INTEGER(KIND=JPIM) :: KINT(IGH)
!
!     this subroutine is a translation of the algol procedure elmhes,
!     num. math. 12, 349-368(1968) by martin and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).
!
!     given a real general matrix, this subroutine
!     reduces a submatrix situated in rows and columns
!     low through igh to upper hessenberg form by
!     stabilized elementary similarity transformations.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        low and igh are integers determined by the balancing
!          subroutine  balanc.  if  balanc  has not been used,
!          set low=1, igh=n.
!
!        a contains the input matrix.
!
!     on output
!
!        a contains the hessenberg matrix.  the multipliers
!          which were used in the reduction are stored in the
!          remaining triangle under the hessenberg matrix.
!
!        kint contains information on the rows and columns
!          interchanged in the reduction.
!          only elements low through igh are used.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('ELMHES',0,ZHOOK_HANDLE)
      LA = IGH - 1
      KP1 = LOW + 1
      IF (LA .LT. KP1) GO TO 200
!
      DO 180 M = KP1, LA
         MM1 = M - 1
         X = 0.0_JPRB
         I = M
!
         DO 100 J = M, IGH
            IF (ABS(A(J,MM1)) .LE. ABS(X)) GO TO 100
            X = A(J,MM1)
            I = J
  100    CONTINUE
!
         KINT(M) = I
         IF (I .EQ. M) GO TO 130
!     .......... interchange rows and columns of a ..........
         DO 110 J = MM1, N
            Y = A(I,J)
            A(I,J) = A(M,J)
            A(M,J) = Y
  110    CONTINUE
!
         DO 120 J = 1, IGH
            Y = A(J,I)
            A(J,I) = A(J,M)
            A(J,M) = Y
  120    CONTINUE
!     .......... end interchange ..........
  130    IF (X .EQ. 0.0_JPRB) GO TO 180
         MP1 = M + 1
!
         DO 160 I = MP1, IGH
            Y = A(I,MM1)
            IF (Y .EQ. 0.0_JPRB) GO TO 160
            Y = Y / X
            A(I,MM1) = Y
!
            DO 140 J = M, N
  140       A(I,J) = A(I,J) - Y * A(M,J)
!
            DO 150 J = 1, IGH
  150       A(J,M) = A(J,M) + Y * A(J,I)
!
  160    CONTINUE
!
  180 CONTINUE
!
  200 CONTINUE
      IF (LHOOK) CALL DR_HOOK('ELMHES',1,ZHOOK_HANDLE)
      RETURN
      ENDSUBROUTINE ELMHES
