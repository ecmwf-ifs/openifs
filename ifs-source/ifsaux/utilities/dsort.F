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

*DECK DSORT
      SUBROUTINE DSORT (DX, DY, N, KFLAG)
C***BEGIN PROLOGUE  DSORT
C***PURPOSE  SORT AN ARRAY AND OPTIONALLY MAKE THE SAME INTERCHANGES IN
C            AN AUXILIARY ARRAY.  THE ARRAY MAY BE SORTED IN INCREASING
C            OR DECREASING ORDER.  A SLIGHTLY MODIFIED QUICKSORT
C            ALGORITHM IS USED.
C***LIBRARY   SLATEC
C***CATEGORY  N6A2B
C***TYPE      DOUBLE PRECISION (SSORT-S, DSORT-D, ISORT-I)
C***KEYWORDS  SINGLETON QUICKSORT, SORT, SORTING
C***AUTHOR  JONES, R. E., (SNLA)
C           WISNIEWSKI, J. A., (SNLA)
C***DESCRIPTION
C
C   DSORT SORTS ARRAY DX AND OPTIONALLY MAKES THE SAME INTERCHANGES IN
C   ARRAY DY.  THE ARRAY DX MAY BE SORTED IN INCREASING ORDER OR
C   DECREASING ORDER.  A SLIGHTLY MODIFIED QUICKSORT ALGORITHM IS USED.
C
C   DESCRIPTION OF PARAMETERS
C      DX - ARRAY OF VALUES TO BE SORTED   (USUALLY ABSCISSAS)
C      DY - ARRAY TO BE (OPTIONALLY) CARRIED ALONG
C      N  - NUMBER OF VALUES IN ARRAY DX TO BE SORTED
C      KFLAG - CONTROL PARAMETER
C            =  2  MEANS SORT DX IN INCREASING ORDER AND CARRY DY ALONG.
C            =  1  MEANS SORT DX IN INCREASING ORDER (IGNORING DY)
C            = -1  MEANS SORT DX IN DECREASING ORDER (IGNORING DY)
C            = -2  MEANS SORT DX IN DECREASING ORDER AND CARRY DY ALONG.
C
C***REFERENCES  R. C. SINGLETON, ALGORITHM 347, AN EFFICIENT ALGORITHM
C                 FOR SORTING WITH MINIMAL STORAGE, COMMUNICATIONS OF
C                 THE ACM, 12, 3 (1969), PP. 185-187.
C***ROUTINES CALLED  XERMSG
C***REVISION HISTORY  (YYMMDD)
C   761101  DATE WRITTEN
C   761118  MODIFIED TO USE THE SINGLETON QUICKSORT ALGORITHM.  (JAW)
C   890531  CHANGED ALL SPECIFIC INTRINSICS TO GENERIC.  (WRB)
C   890831  MODIFIED ARRAY DECLARATIONS.  (WRB)
C   891009  REMOVED UNREFERENCED STATEMENT LABELS.  (WRB)
C   891024  CHANGED CATEGORY.  (WRB)
C   891024  REVISION DATE FROM VERSION 3.2
C   891214  PROLOGUE CONVERTED TO VERSION 4.0 FORMAT.  (BAB)
C   900315  CALLS TO XERROR CHANGED TO CALLS TO XERMSG.  (THJ)
C   901012  DECLARED ALL VARIABLES; CHANGED X,Y TO DX,DY; CHANGED
C           CODE TO PARALLEL SSORT. (M. MCCLAIN)
C   920501  REFORMATTED THE REFERENCES SECTION.  (DWL, WRB)
C   920519  CLARIFIED ERROR MESSAGES.  (DWL)
C   920801  DECLARATIONS SECTION REBUILT AND CODE RESTRUCTURED TO USE
C           IF-THEN-ELSE-ENDIF.  (RWC, WRB)
C***END PROLOGUE  DSORT
C     .. SCALAR ARGUMENTS ..
      INTEGER KFLAG, N
C     .. ARRAY ARGUMENTS ..
      DOUBLE PRECISION DX(*), DY(*)
C     .. LOCAL SCALARS ..
      DOUBLE PRECISION R, T, TT, TTY, TY
      INTEGER I, IJ, J, K, KK, L, M, NN
C     .. LOCAL ARRAYS ..
      INTEGER IL(21), IU(21)
C     .. EXTERNAL SUBROUTINES ..
C      EXTERNAL XERMSG
C     .. INTRINSIC FUNCTIONS ..
      INTRINSIC ABS, INT
C***FIRST EXECUTABLE STATEMENT  DSORT
      NN = N
      IF (NN .LT. 1) THEN
          WRITE(6,*)  
     +   "The number of values to be sorted is not positive."

C         CALL XERMSG ('SLATEC', 'DSORT',
C     +      'The number of values to be sorted is not positive.', 1, 1)

         RETURN
      ENDIF
C
      KK = ABS(KFLAG)
      IF (KK.NE.1 .AND. KK.NE.2) THEN

          WRITE(6,*) 
     +    "The sort control parameter, K is bad."
C         CALL XERMSG ('SLATEC', 'DSORT',
C     +      'The sort control parameter, K, is not 2, 1, -1, or -2.', 2,
C     +      1)
         RETURN
      ENDIF
C
C     ALTER ARRAY DX TO GET DECREASING ORDER IF NEEDED
C
      IF (KFLAG .LE. -1) THEN
         DO 10 I=1,NN
            DX(I) = -DX(I)
   10    CONTINUE
      ENDIF
C
      IF (KK .EQ. 2) GO TO 100
C
C     SORT DX ONLY
C
      M = 1
      I = 1
      J = NN
      R = 0.375D0
C
   20 IF (I .EQ. J) GO TO 60
      IF (R .LE. 0.5898437D0) THEN
         R = R+3.90625D-2
      ELSE
         R = R-0.21875D0
      ENDIF
C
   30 K = I
C
C     SELECT A CENTRAL ELEMENT OF THE ARRAY AND SAVE IT IN LOCATION T
C
      IJ = I + INT((J-I)*R)
      T = DX(IJ)
C
C     IF FIRST ELEMENT OF ARRAY IS GREATER THAN T, INTERCHANGE WITH T
C
      IF (DX(I) .GT. T) THEN
         DX(IJ) = DX(I)
         DX(I) = T
         T = DX(IJ)
      ENDIF
      L = J
C
C     IF LAST ELEMENT OF ARRAY IS LESS THAN THAN T, INTERCHANGE WITH T
C
      IF (DX(J) .LT. T) THEN
         DX(IJ) = DX(J)
         DX(J) = T
         T = DX(IJ)
C
C        IF FIRST ELEMENT OF ARRAY IS GREATER THAN T, INTERCHANGE WITH T
C
         IF (DX(I) .GT. T) THEN
            DX(IJ) = DX(I)
            DX(I) = T
            T = DX(IJ)
         ENDIF
      ENDIF
C
C     FIND AN ELEMENT IN THE SECOND HALF OF THE ARRAY WHICH IS SMALLER
C     THAN T
C
   40 L = L-1
      IF (DX(L) .GT. T) GO TO 40
C
C     FIND AN ELEMENT IN THE FIRST HALF OF THE ARRAY WHICH IS GREATER
C     THAN T
C
   50 K = K+1
      IF (DX(K) .LT. T) GO TO 50
C
C     INTERCHANGE THESE ELEMENTS
C
      IF (K .LE. L) THEN
         TT = DX(L)
         DX(L) = DX(K)
         DX(K) = TT
         GO TO 40
      ENDIF
C
C     SAVE UPPER AND LOWER SUBSCRIPTS OF THE ARRAY YET TO BE SORTED
C
      IF (L-I .GT. J-K) THEN
         IL(M) = I
         IU(M) = L
         I = K
         M = M+1
      ELSE
         IL(M) = K
         IU(M) = J
         J = L
         M = M+1
      ENDIF
      GO TO 70
C
C     BEGIN AGAIN ON ANOTHER PORTION OF THE UNSORTED ARRAY
C
   60 M = M-1
      IF (M .EQ. 0) GO TO 190
      I = IL(M)
      J = IU(M)
C
   70 IF (J-I .GE. 1) GO TO 30
      IF (I .EQ. 1) GO TO 20
      I = I-1
C
   80 I = I+1
      IF (I .EQ. J) GO TO 60
      T = DX(I+1)
      IF (DX(I) .LE. T) GO TO 80
      K = I
C
   90 DX(K+1) = DX(K)
      K = K-1
      IF (T .LT. DX(K)) GO TO 90
      DX(K+1) = T
      GO TO 80
C
C     SORT DX AND CARRY DY ALONG
C
  100 M = 1
      I = 1
      J = NN
      R = 0.375D0
C
  110 IF (I .EQ. J) GO TO 150
      IF (R .LE. 0.5898437D0) THEN
         R = R+3.90625D-2
      ELSE
         R = R-0.21875D0
      ENDIF
C
  120 K = I
C
C     SELECT A CENTRAL ELEMENT OF THE ARRAY AND SAVE IT IN LOCATION T
C
      IJ = I + INT((J-I)*R)
      T = DX(IJ)
      TY = DY(IJ)
C
C     IF FIRST ELEMENT OF ARRAY IS GREATER THAN T, INTERCHANGE WITH T
C
      IF (DX(I) .GT. T) THEN
         DX(IJ) = DX(I)
         DX(I) = T
         T = DX(IJ)
         DY(IJ) = DY(I)
         DY(I) = TY
         TY = DY(IJ)
      ENDIF
      L = J
C
C     IF LAST ELEMENT OF ARRAY IS LESS THAN T, INTERCHANGE WITH T
C
      IF (DX(J) .LT. T) THEN
         DX(IJ) = DX(J)
         DX(J) = T
         T = DX(IJ)
         DY(IJ) = DY(J)
         DY(J) = TY
         TY = DY(IJ)
C
C        IF FIRST ELEMENT OF ARRAY IS GREATER THAN T, INTERCHANGE WITH T
C
         IF (DX(I) .GT. T) THEN
            DX(IJ) = DX(I)
            DX(I) = T
            T = DX(IJ)
            DY(IJ) = DY(I)
            DY(I) = TY
            TY = DY(IJ)
         ENDIF
      ENDIF
C
C     FIND AN ELEMENT IN THE SECOND HALF OF THE ARRAY WHICH IS SMALLER
C     THAN T
C
  130 L = L-1
      IF (DX(L) .GT. T) GO TO 130
C
C     FIND AN ELEMENT IN THE FIRST HALF OF THE ARRAY WHICH IS GREATER
C     THAN T
C
  140 K = K+1
      IF (DX(K) .LT. T) GO TO 140
C
C     INTERCHANGE THESE ELEMENTS
C
      IF (K .LE. L) THEN
         TT = DX(L)
         DX(L) = DX(K)
         DX(K) = TT
         TTY = DY(L)
         DY(L) = DY(K)
         DY(K) = TTY
         GO TO 130
      ENDIF
C
C     SAVE UPPER AND LOWER SUBSCRIPTS OF THE ARRAY YET TO BE SORTED
C
      IF (L-I .GT. J-K) THEN
         IL(M) = I
         IU(M) = L
         I = K
         M = M+1
      ELSE
         IL(M) = K
         IU(M) = J
         J = L
         M = M+1
      ENDIF
      GO TO 160
C
C     BEGIN AGAIN ON ANOTHER PORTION OF THE UNSORTED ARRAY
C
  150 M = M-1
      IF (M .EQ. 0) GO TO 190
      I = IL(M)
      J = IU(M)
C
  160 IF (J-I .GE. 1) GO TO 120
      IF (I .EQ. 1) GO TO 110
      I = I-1
C
  170 I = I+1
      IF (I .EQ. J) GO TO 150
      T = DX(I+1)
      TY = DY(I+1)
      IF (DX(I) .LE. T) GO TO 170
      K = I
C
  180 DX(K+1) = DX(K)
      DY(K+1) = DY(K)
      K = K-1
      IF (T .LT. DX(K)) GO TO 180
      DX(K+1) = T
      DY(K+1) = TY
      GO TO 170
C
C     CLEAN UP
C
  190 IF (KFLAG .LE. -1) THEN
         DO 200 I=1,NN
            DX(I) = -DX(I)
  200    CONTINUE
      ENDIF
      RETURN
      END SUBROUTINE DSORT
