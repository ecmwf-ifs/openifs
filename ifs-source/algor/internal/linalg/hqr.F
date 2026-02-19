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

      SUBROUTINE HQR(NM,N,LOW,IGH,H,WR,WI,IERR)
!AUTOPROMOTE
      USE PARKIND1, ONLY : JPIM, JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
!  RESTORED CORRECT INDICES OF LOOPS (200,210,230,240). (9/29/89 BSG)
!
      IMPLICIT NONE
!
      INTEGER(KIND=JPIM) :: NM
      INTEGER(KIND=JPIM) :: N
      INTEGER(KIND=JPIM) :: LOW
      INTEGER(KIND=JPIM) :: IGH
      INTEGER(KIND=JPIM) :: IERR
      INTEGER(KIND=JPIM) :: I,J,K,L,M,EN,LL,MM,NA,ITN,ITS,MP2,ENM2
      REAL(KIND=JPRB) :: H(NM,N)
      REAL(KIND=JPRB) :: WR(N)
      REAL(KIND=JPRB) :: WI(N)
      REAL(KIND=JPRB) :: P,Q,R,S,T,W,X,Y,ZZ,NORM,TST1,TST2
      LOGICAL :: NOTLAS
!
!     this subroutine is a translation of the algol procedure hqr,
!     num. math. 14, 219-231(1970) by martin, peters, and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 359-371(1971).
!
!     this subroutine finds the eigenvalues of a real
!     upper hessenberg matrix by the qr method.
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
!        h contains the upper hessenberg matrix.  information about
!          the transformations used in the reduction to hessenberg
!          form by  elmhes  or  orthes, if performed, is stored
!          in the remaining triangle under the hessenberg matrix.
!
!     on output
!
!        h has been destroyed.  therefore, it must be saved
!          before calling  hqr  if subsequent calculation and
!          back transformation of eigenvectors is to be performed.
!
!        wr and wi contain the real and imaginary parts,
!          respectively, of the eigenvalues.  the eigenvalues
!          are unordered except that complex conjugate pairs
!          of values appear consecutively with the eigenvalue
!          having the positive imaginary part first.  if an
!          error exit is made, the eigenvalues should be correct
!          for indices ierr+1,...,n.
!
!        ierr is set to
!          zero       for normal return,
!          j          if the limit of 30*n iterations is exhausted
!                     while the j-th eigenvalue is being sought.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated september 1989.
!
!     ------------------------------------------------------------------
!
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('HQR',0,ZHOOK_HANDLE)
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
   60 IF (EN .LT. LOW) GO TO 1001
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
         IF (X .EQ. 0.0_JPRB) GO TO 260
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
         DO 200 J = K, EN
            P = H(K,J) + Q * H(K+1,J)
            H(K,J) = H(K,J) - P * X
            H(K+1,J) = H(K+1,J) - P * Y
  200    CONTINUE
!
         J = MIN(EN,K+3)
!     .......... column modification ..........
         DO 210 I = L, J
            P = X * H(I,K) + Y * H(I,K+1)
            H(I,K) = H(I,K) - P
            H(I,K+1) = H(I,K+1) - P * Q
  210    CONTINUE
         GO TO 255
  225    CONTINUE
!     .......... row modification ..........
         DO 230 J = K, EN
            P = H(K,J) + Q * H(K+1,J) + R * H(K+2,J)
            H(K,J) = H(K,J) - P * X
            H(K+1,J) = H(K+1,J) - P * Y
            H(K+2,J) = H(K+2,J) - P * ZZ
  230    CONTINUE
!
         J = MIN(EN,K+3)
!     .......... column modification ..........
         DO 240 I = L, J
            P = X * H(I,K) + Y * H(I,K+1) + ZZ * H(I,K+2)
            H(I,K) = H(I,K) - P
            H(I,K+1) = H(I,K+1) - P * Q
            H(I,K+2) = H(I,K+2) - P * R
  240    CONTINUE
  255    CONTINUE
!
  260 CONTINUE
!
      GO TO 70
!     .......... one root found ..........
  270 WR(EN) = X + T
      WI(EN) = 0.0_JPRB
      EN = NA
      GO TO 60
!     .......... two roots found ..........
  280 P = (Y - X) / 2.0_JPRB
      Q = P * P + W
      ZZ = SQRT(ABS(Q))
      X = X + T
      IF (Q .LT. 0.0_JPRB) GO TO 320
!     .......... real pair ..........
      ZZ = P + SIGN(ZZ,P)
      WR(NA) = X + ZZ
      WR(EN) = WR(NA)
      IF (ZZ .NE. 0.0_JPRB) WR(EN) = X - W / ZZ
      WI(NA) = 0.0_JPRB
      WI(EN) = 0.0_JPRB
      GO TO 330
!     .......... complex pair ..........
  320 WR(NA) = X + P
      WR(EN) = X + P
      WI(NA) = ZZ
      WI(EN) = -ZZ
  330 EN = ENM2
      GO TO 60
!     .......... set error -- all eigenvalues have not
!                converged after 30*n iterations ..........
 1000 IERR = EN
 1001 CONTINUE
      IF (LHOOK) CALL DR_HOOK('HQR',1,ZHOOK_HANDLE)
      RETURN
      ENDSUBROUTINE HQR
