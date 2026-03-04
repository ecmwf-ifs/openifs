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

      SUBROUTINE BALANC(NM,N,A,LOW,IGH,SCALE)
!AUTOPROMOTE
      USE PARKIND1, ONLY : JPIM, JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
!
      IMPLICIT NONE
!
      INTEGER(KIND=JPIM) :: NM
      INTEGER(KIND=JPIM) :: N
      INTEGER(KIND=JPIM) :: IGH
      INTEGER(KIND=JPIM) :: LOW
      INTEGER(KIND=JPIM) :: I,J,K,L,M,JJ,IEXC
!     DOUBLE PRECISION A(NM,N),SCALE(N)
!     DOUBLE PRECISION C,F,G,R,S,B2,RADIX
      REAL(KIND=JPRB) :: A(NM,N)
      REAL(KIND=JPRB) :: SCALE(N)
      REAL(KIND=JPRB) :: C,F,G,R,S,B2,RADIX
      LOGICAL :: NOCONV
!
!     this subroutine is a translation of the algol procedure balance,
!     num. math. 13, 293-304(1969) by parlett and reinsch.
!     handbook for auto. comp., vol.ii-linear algebra, 315-326(1971).
!
!     this subroutine balances a real matrix and isolates
!     eigenvalues whenever possible.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        a contains the input matrix to be balanced.
!
!     on output
!
!        a contains the balanced matrix.
!
!        low and igh are two integers such that a(i,j)
!          is equal to zero if
!           (1) i is greater than j and
!           (2) j=1,...,low-1 or i=igh+1,...,n.
!
!        scale contains information determining the
!           permutations and scaling factors used.
!
!     suppose that the principal submatrix in rows low through igh
!     has been balanced, that p(j) denotes the index interchanged
!     with j during the permutation step, and that the elements
!     of the diagonal matrix used are denoted by d(i,j).  then
!        scale(j) = p(j),    for j = 1,...,low-1
!                 = d(j,j),      j = low,...,igh
!                 = p(j)         j = igh+1,...,n.
!     the order in which the interchanges are made is n to igh+1,
!     then 1 to low-1.
!
!     note that 1 is returned for igh if igh is zero formally.
!
!     the algol procedure exc contained in balance appears in
!     balanc  in line.  (note that the algol roles of identifiers
!     k,l have been reversed.)
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('BALANC',0,ZHOOK_HANDLE)
      RADIX = 16.0_JPRB
!
      B2 = RADIX * RADIX
      K = 1
      L = N
      GO TO 100
!     .......... in-line procedure for row and
!                column exchange ..........
   20 SCALE(M) = J
      IF (J .EQ. M) GO TO 50
!
      DO 30 I = 1, L
         F = A(I,J)
         A(I,J) = A(I,M)
         A(I,M) = F
   30 CONTINUE
!
      DO 40 I = K, N
         F = A(J,I)
         A(J,I) = A(M,I)
         A(M,I) = F
   40 CONTINUE
!
   50 GO TO (80,130), IEXC
!     .......... search for rows isolating an eigenvalue
!                and push them down ..........
   80 IF (L .EQ. 1) GO TO 280
      L = L - 1
!     .......... for j=l step -1 until 1 do -- ..........
  100 DO 120 JJ = 1, L
         J = L + 1 - JJ
!
         DO 110 I = 1, L
            IF (I .EQ. J) GO TO 110
            IF (A(J,I) .NE. 0.0_JPRB) GO TO 120
  110    CONTINUE
!
         M = L
         IEXC = 1
         GO TO 20
  120 CONTINUE
!
      GO TO 140
!     .......... search for columns isolating an eigenvalue
!                and push them left ..........
  130 K = K + 1
!
  140 DO 170 J = K, L
!
         DO 150 I = K, L
            IF (I .EQ. J) GO TO 150
            IF (A(I,J) .NE. 0.0_JPRB) GO TO 170
  150    CONTINUE
!
         M = K
         IEXC = 2
         GO TO 20
  170 CONTINUE
!     .......... now balance the submatrix in rows k to l ..........
      DO 180 I = K, L
  180 SCALE(I) = 1.0_JPRB
!     .......... iterative loop for norm reduction ..........
  190 NOCONV = .FALSE.
!
      DO 270 I = K, L
         C = 0.0_JPRB
         R = 0.0_JPRB
!
         DO 200 J = K, L
            IF (J .EQ. I) GO TO 200
            C = C + ABS(A(J,I))
            R = R + ABS(A(I,J))
  200    CONTINUE
!     .......... guard against zero c or r due to underflow ..........
         IF (C .EQ. 0.0_JPRB .OR. R .EQ. 0.0_JPRB) GO TO 270
         G = R / RADIX
         F = 1.0_JPRB
         S = C + R
  210    IF (C .GE. G) GO TO 220
         F = F * RADIX
         C = C * B2
         GO TO 210
  220    G = R * RADIX
  230    IF (C .LT. G) GO TO 240
         F = F / RADIX
         C = C / B2
         GO TO 230
!     .......... now balance ..........
  240    IF ((C + R) / F .GE. 0.95_JPRB * S) GO TO 270
         G = 1.0_JPRB / F
         SCALE(I) = SCALE(I) * F
         NOCONV = .TRUE.
!
         DO 250 J = K, N
  250    A(I,J) = A(I,J) * G
!
         DO 260 J = 1, L
  260    A(J,I) = A(J,I) * F
!
  270 CONTINUE
!
      IF (NOCONV) GO TO 190
!
  280 LOW = K
      IGH = L
      IF (LHOOK) CALL DR_HOOK('BALANC',1,ZHOOK_HANDLE)
      ENDSUBROUTINE BALANC
