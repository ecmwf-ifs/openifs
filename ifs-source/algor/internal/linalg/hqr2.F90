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

      SUBROUTINE HQR2(NM,N,LOW,IGH,H,WR,WI,PZ,IERR)
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
      INTEGER(KIND=JPIM) :: IERR
      INTEGER(KIND=JPIM) :: I,J,K,L,M,EN,II,JJ,LL,MM,NA,NN
      INTEGER(KIND=JPIM) :: ITN,ITS,MP2,ENM2
      REAL(KIND=JPRB) :: H(NM,N)
      REAL(KIND=JPRB) :: WR(N)
      REAL(KIND=JPRB) :: WI(N)
      REAL(KIND=JPRB) :: PZ(NM,N)
      REAL(KIND=JPRB) :: P,Q,R,S,T,W,X,Y,RA,SA,VI,VR,ZZ,NORM,TST1,TST2
      LOGICAL :: NOTLAS
!
!     this subroutine is a translation of the algol procedure hqr2,
!     num. math. 16, 181-204(1970) by peters and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).
!
!     this subroutine finds the eigenvalues and eigenvectors
!     of a real upper hessenberg matrix by the qr method.  the
!     eigenvectors of a real general matrix can also be found
!     if  elmhes  and  eltran  or  orthes  and  ortran  have
!     been used to reduce this general matrix to hessenberg form
!     and to accumulate the similarity transformations.
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
!        h contains the upper hessenberg matrix.
!
!        Pz contains the transformation matrix produced by  eltran
!          after the reduction by  elmhes, or by  ortran  after the
!          reduction by  orthes, if performed.  if the eigenvectors
!          of the hessenberg matrix are desired, Pz must contain the
!          identity matrix.
!
!     on output
!
!        h has been destroyed.
!
!        wr and wi contain the real and imaginary parts,
!          respectively, of the eigenvalues.  the eigenvalues
!          are unordered except that complex conjugate pairs
!          of values appear consecutively with the eigenvalue
!          having the positive imaginary part first.  if an
!          error exit is made, the eigenvalues should be correct
!          for indices ierr+1,...,n.
!
!        Pz contains the real and imaginary parts of the eigenvectors.
!          if the i-th eigenvalue is real, the i-th column of Pz
!          contains its eigenvector.  if the i-th eigenvalue is complex
!          with positive imaginary part, the i-th and (i+1)-th
!          columns of Pz contain the real and imaginary parts of its
!          eigenvector.  the eigenvectors are unnormalized.  if an
!          error exit is made, none of the eigenvectors has been found.
!
!        ierr is set to
!          zero       for normal return,
!          j          if the limit of 30*n iterations is exhausted
!                     while the j-th eigenvalue is being sought.
!
!     calls cdiv for complex division.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('HQR2',0,ZHOOK_HANDLE)
      IERR = 0
      NORM = 0.0_JPRB
      K = 1
!     .......... store roots isolated by balanc
!                and compute matrix norm ..........
      DO 50 I = 1, N
!
         DO 40 J = K, N
   40    NORM = NORM + ABS(H(I,J))
!
         K = I
         IF (I .GE. LOW .AND. I .LE. IGH) GO TO 50
         WR(I) = H(I,I)
         WI(I) = 0.0_JPRB
   50 CONTINUE
!
      EN = IGH
      T = 0.0_JPRB
      ITN = 30*N
!     .......... search for next eigenvalues ..........
   60 IF (EN .LT. LOW) GO TO 340
      ITS = 0
      NA = EN - 1
      ENM2 = NA - 1
!     .......... look for single small sub-diagonal element
!                for l=en step -1 until low do -- ..........
   70 DO 80 LL = LOW, EN
         L = EN + LOW - LL
         IF (L .EQ. LOW) GO TO 100
         S = ABS(H(L-1,L-1)) + ABS(H(L,L))
         IF (S .EQ. 0.0_JPRB) S = NORM
         TST1 = S
         TST2 = TST1 + ABS(H(L,L-1))
         IF (TST2 .EQ. TST1) GO TO 100
   80 CONTINUE
!     .......... form shift ..........
  100 X = H(EN,EN)
      IF (L .EQ. EN) GO TO 270
      Y = H(NA,NA)
      W = H(EN,NA) * H(NA,EN)
      IF (L .EQ. NA) GO TO 280
      IF (ITN .EQ. 0) GO TO 1000
      IF (ITS .NE. 10 .AND. ITS .NE. 20) GO TO 130
!     .......... form exceptional shift ..........
      T = T + X
!
      DO 120 I = LOW, EN
  120 H(I,I) = H(I,I) - X
!
      S = ABS(H(EN,NA)) + ABS(H(NA,ENM2))
      X = 0.75_JPRB * S
      Y = X
      W = -0.4375_JPRB * S * S
  130 ITS = ITS + 1
      ITN = ITN - 1
!     .......... look for two consecutive small
!                sub-diagonal elements.
!                for m=en-2 step -1 until l do -- ..........
      DO 140 MM = L, ENM2
         M = ENM2 + L - MM
         ZZ = H(M,M)
         R = X - ZZ
         S = Y - ZZ
         P = (R * S - W) / H(M+1,M) + H(M,M+1)
         Q = H(M+1,M+1) - ZZ - R - S
         R = H(M+2,M+1)
         S = ABS(P) + ABS(Q) + ABS(R)
         P = P / S
         Q = Q / S
         R = R / S
         IF (M .EQ. L) GO TO 150
         TST1 = ABS(P)*(ABS(H(M-1,M-1)) + ABS(ZZ) + ABS(H(M+1,M+1)))
         TST2 = TST1 + ABS(H(M,M-1))*(ABS(Q) + ABS(R))
         IF (TST2 .EQ. TST1) GO TO 150
  140 CONTINUE
!
  150 MP2 = M + 2
!
      DO 160 I = MP2, EN
         H(I,I-2) = 0.0_JPRB
         IF (I .EQ. MP2) GO TO 160
         H(I,I-3) = 0.0_JPRB
  160 CONTINUE
!     .......... double qr step involving rows l to en and
!                columns m to en ..........
      DO 260 K = M, NA
         NOTLAS = K .NE. NA
         IF (K .EQ. M) GO TO 170
         P = H(K,K-1)
         Q = H(K+1,K-1)
         R = 0.0_JPRB
         IF (NOTLAS) R = H(K+2,K-1)
         X = ABS(P) + ABS(Q) + ABS(R)
         IF (X .EQ. 0.0E0) GO TO 260
         P = P / X
         Q = Q / X
         R = R / X
  170    S = SIGN(SQRT(P*P+Q*Q+R*R),P)
         IF (K .EQ. M) GO TO 180
         H(K,K-1) = -S * X
         GO TO 190
  180    IF (L .NE. M) H(K,K-1) = -H(K,K-1)
  190    P = P + S
         X = P / S
         Y = Q / S
         ZZ = R / S
         Q = Q / P
         R = R / P
         IF (NOTLAS) GO TO 225
!     .......... row modification ..........
         DO 200 J = K, N
            P = H(K,J) + Q * H(K+1,J)
            H(K,J) = H(K,J) - P * X
            H(K+1,J) = H(K+1,J) - P * Y
  200    CONTINUE
!
         J = MIN(EN,K+3)
!     .......... column modification ..........
         DO 210 I = 1, J
            P = X * H(I,K) + Y * H(I,K+1)
            H(I,K) = H(I,K) - P
            H(I,K+1) = H(I,K+1) - P * Q
  210    CONTINUE
!     .......... accumulate transformations ..........
         DO 220 I = LOW, IGH
            P = X * PZ(I,K) + Y * PZ(I,K+1)
            PZ(I,K) = PZ(I,K) - P
            PZ(I,K+1) = PZ(I,K+1) - P * Q
  220    CONTINUE
         GO TO 255
  225    CONTINUE
!     .......... row modification ..........
         DO 230 J = K, N
            P = H(K,J) + Q * H(K+1,J) + R * H(K+2,J)
            H(K,J) = H(K,J) - P * X
            H(K+1,J) = H(K+1,J) - P * Y
            H(K+2,J) = H(K+2,J) - P * ZZ
  230    CONTINUE
!
         J = MIN(EN,K+3)
!     .......... column modification ..........
         DO 240 I = 1, J
            P = X * H(I,K) + Y * H(I,K+1) + ZZ * H(I,K+2)
            H(I,K) = H(I,K) - P
            H(I,K+1) = H(I,K+1) - P * Q
            H(I,K+2) = H(I,K+2) - P * R
  240    CONTINUE
!     .......... accumulate transformations ..........
         DO 250 I = LOW, IGH
            P = X * PZ(I,K) + Y * PZ(I,K+1) + ZZ * PZ(I,K+2)
            PZ(I,K) = PZ(I,K) - P
            PZ(I,K+1) = PZ(I,K+1) - P * Q
            PZ(I,K+2) = PZ(I,K+2) - P * R
  250    CONTINUE
  255    CONTINUE
!
  260 CONTINUE
!
      GO TO 70
!     .......... one root found ..........
  270 H(EN,EN) = X + T
      WR(EN) = H(EN,EN)
      WI(EN) = 0.0_JPRB
      EN = NA
      GO TO 60
!     .......... two roots found ..........
  280 P = (Y - X) / 2.0_JPRB
      Q = P * P + W
      ZZ = SQRT(ABS(Q))
      H(EN,EN) = X + T
      X = H(EN,EN)
      H(NA,NA) = Y + T
      IF (Q .LT. 0.0_JPRB) GO TO 320
!     .......... real pair ..........
      ZZ = P + SIGN(ZZ,P)
      WR(NA) = X + ZZ
      WR(EN) = WR(NA)
      IF (ZZ .NE. 0.0_JPRB) WR(EN) = X - W / ZZ
      WI(NA) = 0.0_JPRB
      WI(EN) = 0.0_JPRB
      X = H(EN,NA)
      S = ABS(X) + ABS(ZZ)
      P = X / S
      Q = ZZ / S
      R = SQRT(P*P+Q*Q)
      P = P / R
      Q = Q / R
!     .......... row modification ..........
      DO 290 J = NA, N
         ZZ = H(NA,J)
         H(NA,J) = Q * ZZ + P * H(EN,J)
         H(EN,J) = Q * H(EN,J) - P * ZZ
  290 CONTINUE
!     .......... column modification ..........
      DO 300 I = 1, EN
         ZZ = H(I,NA)
         H(I,NA) = Q * ZZ + P * H(I,EN)
         H(I,EN) = Q * H(I,EN) - P * ZZ
  300 CONTINUE
!     .......... accumulate transformations ..........
      DO 310 I = LOW, IGH
         ZZ = PZ(I,NA)
         PZ(I,NA) = Q * ZZ + P * PZ(I,EN)
         PZ(I,EN) = Q * PZ(I,EN) - P * ZZ
  310 CONTINUE
!
      GO TO 330
!     .......... complex pair ..........
  320 WR(NA) = X + P
      WR(EN) = X + P
      WI(NA) = ZZ
      WI(EN) = -ZZ
  330 EN = ENM2
      GO TO 60
!     .......... all roots found.  backsubstitute to find
!                vectors of upper triangular form ..........
  340 IF (NORM .EQ. 0.0_JPRB) GO TO 1001
!     .......... for en=n step -1 until 1 do -- ..........
      DO 800 NN = 1, N
         EN = N + 1 - NN
         P = WR(EN)
         Q = WI(EN)
         NA = EN - 1
         IF (Q) 710, 600, 800
!     .......... real vector ..........
  600    M = EN
         H(EN,EN) = 1.0_JPRB
         IF (NA .EQ. 0) GO TO 800
!     .......... for i=en-1 step -1 until 1 do -- ..........
         DO 700 II = 1, NA
            I = EN - II
            W = H(I,I) - P
            R = 0.0_JPRB
!
            DO 610 J = M, EN
  610       R = R + H(I,J) * H(J,EN)
!
            IF (WI(I) .GE. 0.0_JPRB) GO TO 630
            ZZ = W
            S = R
            GO TO 700
  630       M = I
            IF (WI(I) .NE. 0.0_JPRB) GO TO 640
            T = W
            IF (T .NE. 0.0_JPRB) GO TO 635
               TST1 = NORM
               T = TST1
  632          T = 0.01_JPRB * T
               TST2 = NORM + T
               IF (TST2 .GT. TST1) GO TO 632
  635       H(I,EN) = -R / T
            GO TO 680
!     .......... solve real equations ..........
  640       X = H(I,I+1)
            Y = H(I+1,I)
            Q = (WR(I) - P) * (WR(I) - P) + WI(I) * WI(I)
            T = (X * S - ZZ * R) / Q
            H(I,EN) = T
            IF (ABS(X) .LE. ABS(ZZ)) GO TO 650
            H(I+1,EN) = (-R - W * T) / X
            GO TO 680
  650       H(I+1,EN) = (-S - Y * T) / ZZ
!
!     .......... overflow control ..........
  680       T = ABS(H(I,EN))
            IF (T .EQ. 0.0_JPRB) GO TO 700
            TST1 = T
            TST2 = TST1 + 1.0_JPRB/TST1
            IF (TST2 .GT. TST1) GO TO 700
            DO 690 J = I, EN
               H(J,EN) = H(J,EN)/T
  690       CONTINUE
!
  700    CONTINUE
!     .......... end real vector ..........
         GO TO 800
!     .......... complex vector ..........
  710    M = NA
!     .......... last vector component chosen imaginary so that
!                eigenvector matrix is triangular ..........
         IF (ABS(H(EN,NA)) .LE. ABS(H(NA,EN))) GO TO 720
         H(NA,NA) = Q / H(EN,NA)
         H(NA,EN) = -(H(EN,EN) - P) / H(EN,NA)
         GO TO 730
  720    CALL CDIV(0.0_JPRB,-H(NA,EN),H(NA,NA)-P,Q,H(NA,NA),H(NA,EN))
  730    H(EN,NA) = 0.0_JPRB
         H(EN,EN) = 1.0_JPRB
         ENM2 = NA - 1
         IF (ENM2 .EQ. 0) GO TO 800
!     .......... for i=en-2 step -1 until 1 do -- ..........
         DO 795 II = 1, ENM2
            I = NA - II
            W = H(I,I) - P
            RA = 0.0_JPRB
            SA = 0.0_JPRB
!
            DO 760 J = M, EN
               RA = RA + H(I,J) * H(J,NA)
               SA = SA + H(I,J) * H(J,EN)
  760       CONTINUE
!
            IF (WI(I) .GE. 0.0E0) GO TO 770
            ZZ = W
            R = RA
            S = SA
            GO TO 795
  770       M = I
            IF (WI(I) .NE. 0.0_JPRB) GO TO 780
            CALL CDIV(-RA,-SA,W,Q,H(I,NA),H(I,EN))
            GO TO 790
!     .......... solve complex equations ..........
  780       X = H(I,I+1)
            Y = H(I+1,I)
            VR = (WR(I) - P) * (WR(I) - P) + WI(I) * WI(I) - Q * Q
            VI = (WR(I) - P) * 2.0_JPRB * Q
            IF (VR .NE. 0.0_JPRB .OR. VI .NE. 0.0_JPRB) GO TO 784
               TST1 = NORM * (ABS(W) + ABS(Q) + ABS(X)                  &
     &                      + ABS(Y) + ABS(ZZ))
               VR = TST1
  783          VR = 0.01_JPRB * VR
               TST2 = TST1 + VR
               IF (TST2 .GT. TST1) GO TO 783
  784       CALL CDIV(X*R-ZZ*RA+Q*SA,X*S-ZZ*SA-Q*RA,VR,VI,              &
     &                H(I,NA),H(I,EN))
            IF (ABS(X) .LE. ABS(ZZ) + ABS(Q)) GO TO 785
            H(I+1,NA) = (-RA - W * H(I,NA) + Q * H(I,EN)) / X
            H(I+1,EN) = (-SA - W * H(I,EN) - Q * H(I,NA)) / X
            GO TO 790
  785       CALL CDIV(-R-Y*H(I,NA),-S-Y*H(I,EN),ZZ,Q,                   &
     &                H(I+1,NA),H(I+1,EN))
!
!     .......... overflow control ..........
  790       T = MAX(ABS(H(I,NA)), ABS(H(I,EN)))
            IF (T .EQ. 0.0_JPRB) GO TO 795
            TST1 = T
            TST2 = TST1 + 1.0_JPRB/TST1
            IF (TST2 .GT. TST1) GO TO 795
            DO 792 J = I, EN
               H(J,NA) = H(J,NA)/T
               H(J,EN) = H(J,EN)/T
  792       CONTINUE
!
  795    CONTINUE
!     .......... end complex vector ..........
  800 CONTINUE
!     .......... end back substitution.
!                vectors of isolated roots ..........
      DO 840 I = 1, N
         IF (I .GE. LOW .AND. I .LE. IGH) GO TO 840
!
         DO 820 J = I, N
  820    PZ(I,J) = H(I,J)
!
  840 CONTINUE
!     .......... multiply by transformation matrix to give
!                vectors of original full matrix.
!                for j=n step -1 until low do -- ..........
      DO 880 JJ = LOW, N
         J = N + LOW - JJ
         M = MIN(J,IGH)
!
         DO 880 I = LOW, IGH
            ZZ = 0.0_JPRB
!
            DO 860 K = LOW, M
  860       ZZ = ZZ + PZ(I,K) * H(K,J)
!
            PZ(I,J) = ZZ
  880 CONTINUE
!
      GO TO 1001
!     .......... set error -- all eigenvalues have not
!                converged after 30*n iterations ..........
 1000 IERR = EN
 1001 CONTINUE
      IF (LHOOK) CALL DR_HOOK('HQR2',1,ZHOOK_HANDLE)
      RETURN
      ENDSUBROUTINE HQR2
