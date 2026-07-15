! (C) Copyright 1998- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

!
!     SUBROUTINE 'FFT992' - MULTIPLE FAST REAL PERIODIC TRANSFORM
!
!     Author: Clive Temperton, January 1998
!
!     This routine is a modernized and enhanced version of FFT991
!         - Cray directives and ancient Fortran constructs removed
!         - "vector chopping" removed
!         - WORK array is now dynamically allocated
!         - stride in WORK array is now always 1
!
!     REAL TRANSFORM OF LENGTH N PERFORMED BY REMOVING REDUNDANT
!     OPERATIONS FROM COMPLEX TRANSFORM OF LENGTH N
!
!     A IS THE ARRAY CONTAINING INPUT & OUTPUT DATA
!     TRIGS IS A PREVIOUSLY PREPARED LIST OF TRIG FUNCTION VALUES
!     IFAX IS A PREVIOUSLY PREPARED LIST OF FACTORS OF N
!     INC IS THE INCREMENT WITHIN EACH DATA 'VECTOR'
!         (E.G. INC=1 FOR CONSECUTIVELY STORED DATA)
!     JUMP IS THE INCREMENT BETWEEN THE START OF EACH DATA VECTOR
!     N IS THE LENGTH OF THE DATA VECTORS
!     LOT IS THE NUMBER OF DATA VECTORS
!     ISIGN = +1 FOR TRANSFORM FROM SPECTRAL TO GRIDPOINT
!           = -1 FOR TRANSFORM FROM GRIDPOINT TO SPECTRAL
!
!     ORDERING OF COEFFICIENTS:
!         A(0),B(0),A(1),B(1),A(2),B(2),...,A(N/2),B(N/2)
!         WHERE B(0)=B(N/2)=0; (N+2) LOCATIONS REQUIRED
!
!     ORDERING OF DATA:
!         X(0),X(1),X(2),...,X(N-1), 0 , 0 ; (N+2) LOCATIONS REQUIRED
!
!     VECTORIZATION IS ACHIEVED BY DOING THE TRANSFORMS IN PARALLEL
!
!     N MUST BE COMPOSED OF FACTORS 2,3 & 5 BUT DOES NOT HAVE TO BE EVEN
!
!     DEFINITION OF TRANSFORMS:
!     -------------------------
!
!     ISIGN=+1: X(J)=SUM(K=0,...,N-1)(C(K)*EXP(2*I*J*K*PI/N))
!         WHERE C(K)=A(K)+I*B(K) AND C(N-K)=A(K)-I*B(K)
!
!     ISIGN=-1: A(K)=(1/N)*SUM(J=0,...,N-1)(X(J)*COS(2*J*K*PI/N))
!               B(K)=-(1/N)*SUM(J=0,...,N-1)(X(J)*SIN(2*J*K*PI/N))
!
#ifdef MATHKEISAN
! MathKeisan is a scientific library optimized for NEC (www.mathkeisan.com)

      SUBROUTINE FFT992(A,TRIGS_,IFAX_,INC,JUMP,N,LOT,ISIGN)
!AUTOPROMOTE
      USE PARKIND1, ONLY : JPIM, JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
      IMPLICIT NONE
      INTEGER(KIND=JPIM) :: N
      REAL(KIND=JPRB) :: A(*)
      REAL(KIND=JPRB) :: TRIGS_(N)
      INTEGER(KIND=JPIM) :: IFAX_(10)

      INTEGER(KIND=JPIM) :: INC
      INTEGER(KIND=JPIM) :: JUMP
      INTEGER(KIND=JPIM) :: LOT
      INTEGER(KIND=JPIM) :: ISIGN

      REAL(KIND=JPRB),ALLOCATABLE,DIMENSION(:),SAVE :: WORK , TRIGS
      INTEGER(KIND=JPIM),SAVE :: IFAX (32)


      INTEGER(KIND=JPIM), SAVE ::  N_OLD=-1
      INTEGER(KIND=JPIM), SAVE ::  LOT_OLD=-1

!$OMP threadprivate(ifax,n_old,lot_old,trigs,work)


      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('FFT992',0,ZHOOK_HANDLE)
      IF (N .NE. N_OLD) THEN

         IF( ALLOCATED( WORK ) ) DEALLOCATE( WORK )
         IF( ALLOCATED( TRIGS ) ) DEALLOCATE( TRIGS )

         ALLOCATE(WORK(3*N*LOT))
         ALLOCATE(TRIGS(2*N))

         CALL DFTFAX ( N, IFAX, TRIGS )

         N_OLD = N
         LOT_OLD = LOT

      ELSE

       IF (LOT .GT. LOT_OLD) THEN

         IF( ALLOCATED( WORK ) ) DEALLOCATE( WORK )
         ALLOCATE(WORK(3*N*LOT))
         LOT_OLD = LOT

       ENDIF

      ENDIF

      CALL DFFTMLT ( A, WORK, TRIGS, IFAX, INC, JUMP, N, LOT, ISIGN )

      IF (LHOOK) CALL DR_HOOK('FFT992',1,ZHOOK_HANDLE)
      RETURN

      END SUBROUTINE FFT992
#else
!
!     SUBROUTINE 'FFT992' - MULTIPLE FAST REAL PERIODIC TRANSFORM
!
!     Author: Clive Temperton, January 1998
!
!     This routine is a modernized and enhanced version of FFT991
!         - Cray directives and ancient Fortran constructs removed
!         - "vector chopping" removed
!         - WORK array is now dynamically allocated
!         - stride in WORK array is now always 1
!
!     REAL TRANSFORM OF LENGTH N PERFORMED BY REMOVING REDUNDANT
!     OPERATIONS FROM COMPLEX TRANSFORM OF LENGTH N
!
!     A IS THE ARRAY CONTAINING INPUT & OUTPUT DATA
!     TRIGS IS A PREVIOUSLY PREPARED LIST OF TRIG FUNCTION VALUES
!     IFAX IS A PREVIOUSLY PREPARED LIST OF FACTORS OF N
!     INC IS THE INCREMENT WITHIN EACH DATA 'VECTOR'
!         (E.G. INC=1 FOR CONSECUTIVELY STORED DATA)
!     JUMP IS THE INCREMENT BETWEEN THE START OF EACH DATA VECTOR
!     N IS THE LENGTH OF THE DATA VECTORS
!     LOT IS THE NUMBER OF DATA VECTORS
!     ISIGN = +1 FOR TRANSFORM FROM SPECTRAL TO GRIDPOINT
!           = -1 FOR TRANSFORM FROM GRIDPOINT TO SPECTRAL
!
!     ORDERING OF COEFFICIENTS:
!         A(0),B(0),A(1),B(1),A(2),B(2),...,A(N/2),B(N/2)
!         WHERE B(0)=B(N/2)=0; (N+2) LOCATIONS REQUIRED
!
!     ORDERING OF DATA:
!         X(0),X(1),X(2),...,X(N-1), 0 , 0 ; (N+2) LOCATIONS REQUIRED
!
!     VECTORIZATION IS ACHIEVED BY DOING THE TRANSFORMS IN PARALLEL
!
!     N MUST BE COMPOSED OF FACTORS 2,3 & 5 BUT DOES NOT HAVE TO BE EVEN
!
!     DEFINITION OF TRANSFORMS:
!     -------------------------
!
!     ISIGN=+1: X(J)=SUM(K=0,...,N-1)(C(K)*EXP(2*I*J*K*PI/N))
!         WHERE C(K)=A(K)+I*B(K) AND C(N-K)=A(K)-I*B(K)
!
!     ISIGN=-1: A(K)=(1/N)*SUM(J=0,...,N-1)(X(J)*COS(2*J*K*PI/N))
!               B(K)=-(1/N)*SUM(J=0,...,N-1)(X(J)*SIN(2*J*K*PI/N))
!
      SUBROUTINE FFT992(A,TRIGS,IFAX,INC,JUMP,N,LOT,ISIGN)
!disabled for now. REK.!DEC$ OPTIMIZE:3
!AUTOPROMOTE
      USE PARKIND1, ONLY : JPIM, JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
      IMPLICIT NONE
!      
      INTEGER(KIND=JPIM) :: NU
      INTEGER(KIND=JPIM) :: N
      INTEGER(KIND=JPIM) :: INC
      INTEGER(KIND=JPIM) :: JBASE
      INTEGER(KIND=JPIM) :: JUMP
      INTEGER(KIND=JPIM) :: J,JJ,JUMPA
      INTEGER(KIND=JPIM) :: LOT
      INTEGER(KIND=JPIM) :: K,LA,NFAX
      INTEGER(KIND=JPIM) :: ISIGN
      INTEGER(KIND=JPIM) :: I,IA,IBASE,IERR,IFAC,IGO,II,INCA,IX

      REAL(KIND=JPRB) :: A(*)
      REAL(KIND=JPRB) :: TRIGS(N)
      INTEGER(KIND=JPIM) :: IFAX(10)
!     Dynamically allocated work array:
      REAL(KIND=JPRB) :: WORK(N*LOT+1)
      LOGICAL :: LIPL
!
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('FFT992',0,ZHOOK_HANDLE)
      NFAX=IFAX(1)
      IF (ISIGN.EQ.+1) THEN
!
!     ISIGN=+1, SPECTRAL TO GRIDPOINT TRANSFORM
!     -----------------------------------------
!
        I=1
!OCL NOVREC
!DEC$ IVDEP
        DO J=1,LOT
          A(I+INC)=0.5_JPRB*A(I)
          I=I+JUMP
        ENDDO
        IF (MOD(N,2).EQ.0) THEN
          I=N*INC+1
!OCL NOVREC
!DEC$ IVDEP
          DO J=1,LOT
            A(I)=0.5_JPRB*A(I)
            I=I+JUMP
          ENDDO
        ENDIF
!
        IA=INC+1
        LA=1
        IGO=+1
!
        DO K=1,NFAX
          IFAC=IFAX(K+1)
          IERR=-1
          IF (K.EQ.NFAX.AND.NFAX.GT.2.AND.IGO.EQ.+1) THEN
            LIPL=.TRUE.
          ELSE
            LIPL=.FALSE.
          ENDIF
          IF (INC.EQ.1.AND.JUMP.LT.(2*N).AND.                           &
     &        K.GT.1.AND.K.LT.(NFAX-MOD(NFAX,2))) THEN
            INCA=LOT
            JUMPA=1
          ELSE
            INCA=INC
            JUMPA=JUMP
          ENDIF
          IF (IGO.EQ.+1) THEN
!DEC$ FORCEINLINE
           CALL RPASSF(A(IA),A(IA+LA*INCA),WORK(1),WORK(IFAC*LA*LOT+1), &
     &                  TRIGS,INCA,LOT,JUMPA,1,LOT,N,IFAC,LA,IERR,LIPL)
          ELSE
!DEC$ FORCEINLINE
           CALL RPASSF(WORK(1),WORK(LA*LOT+1),A(IA),A(IA+IFAC*LA*INCA), &
     &                  TRIGS,LOT,INCA,1,JUMPA,LOT,N,IFAC,LA,IERR,LIPL)
          ENDIF
          IF (IERR.NE.0) THEN
            IF (IERR.EQ.2) WRITE(6,901) IFAC
            IF (IERR.EQ.3) WRITE(6,902) IFAC
            IF (LHOOK) CALL DR_HOOK('FFT992',1,ZHOOK_HANDLE)
            RETURN
          ENDIF
          LA=IFAC*LA
          IGO=-IGO
          IA=1
        ENDDO
!
!     IF NECESSARY, COPY RESULTS BACK TO A
!     ------------------------------------
        IF (NFAX.EQ.1) THEN
          IBASE=1
          JBASE=1
          DO JJ=1,N
            I=IBASE
            J=JBASE
            DO II=1,LOT
              A(J)=WORK(I)
              I=I+1
              J=J+JUMP
            ENDDO
            IBASE=IBASE+LOT
            JBASE=JBASE+INC
          ENDDO
        ENDIF
!
!     FILL IN ZEROS AT END
!     --------------------
        IX=N*INC+1
!OCL NOVREC
!DEC$ IVDEP
        DO J=1,LOT
          A(IX)=0.0_JPRB
          A(IX+INC)=0.0_JPRB
          IX=IX+JUMP
        ENDDO
!
      ELSEIF (ISIGN.EQ.-1) THEN
!
!     ISIGN=-1, GRIDPOINT TO SPECTRAL TRANSFORM
!     -----------------------------------------
        IA=1
        LA=N
        IGO=+1
!
        DO K=1,NFAX
          IFAC=IFAX(NFAX+2-K)
          LA=LA/IFAC
          IERR=-1
          IF (K.EQ.1.AND.NFAX.GT.2.AND.MOD(NFAX,2).EQ.1) THEN
            LIPL=.TRUE.
          ELSE
            LIPL=.FALSE.
          ENDIF
          IF (INC.EQ.1.AND.JUMP.LT.(2*N).AND.                           &
     &        K.GT.(1+MOD(NFAX,2)).AND.K.LT.NFAX) THEN
            INCA=LOT
            JUMPA=1
          ELSE
            INCA=INC
            JUMPA=JUMP
          ENDIF
          IF (IGO.EQ.+1) THEN
!DEC$ FORCEINLINE
           CALL QPASSF(A(IA),A(IA+IFAC*LA*INCA),WORK(1),WORK(LA*LOT+1), &
     &                  TRIGS,INCA,LOT,JUMPA,1,LOT,N,IFAC,LA,IERR,LIPL)
          ELSE
!DEC$ FORCEINLINE
           CALL QPASSF(WORK(1),WORK(IFAC*LA*LOT+1),A(IA),A(IA+LA*INCA), &
     &                  TRIGS,LOT,INCA,1,JUMPA,LOT,N,IFAC,LA,IERR,LIPL)
          ENDIF
          IF (IERR.NE.0) THEN
            IF (IERR.EQ.2) WRITE(6,901) IFAC
            IF (IERR.EQ.3) WRITE(6,902) IFAC
            IF (LHOOK) CALL DR_HOOK('FFT992',1,ZHOOK_HANDLE)
            RETURN
          ENDIF
          IF (LIPL) THEN
            IA=1
          ELSE
            IGO=-IGO
            IA=INC+1
          ENDIF
        ENDDO
!
!     IF NECESSARY, COPY RESULTS BACK TO A
!     ------------------------------------
        IF (NFAX.EQ.1) THEN
          IBASE=1
          JBASE=INC+1
          DO JJ=1,N
            I=IBASE
            J=JBASE
            DO II=1,LOT
              A(J)=WORK(I)
              I=I+1
              J=J+JUMP
            ENDDO
            IBASE=IBASE+LOT
            JBASE=JBASE+INC
          ENDDO
        ENDIF
!
!     SHIFT A(0) & FILL IN ZERO IMAG PARTS
!     ------------------------------------
        IX=1
!OCL NOVREC
!DEC$ IVDEP
        DO J=1,LOT
          A(IX)=A(IX+INC)
          A(IX+INC)=0.0_JPRB
          IX=IX+JUMP
        ENDDO
        IF (MOD(N,2).EQ.0) THEN
          IX=(N+1)*INC+1
          DO J=1,LOT
            A(IX)=0.0_JPRB
            IX=IX+JUMP
          ENDDO
        ENDIF
!
      ENDIF
!
!     FORMAT STATEMENTS FOR ERROR MESSAGES:
  901 FORMAT(' FACTOR =',I3,' NOT CATERED FOR')
  902 FORMAT(' FACTOR =',I3,' ONLY CATERED FOR IF LA*IFAC=N')
!
      IF (LHOOK) CALL DR_HOOK('FFT992',1,ZHOOK_HANDLE)

      CONTAINS
!     SUBROUTINE 'RPASSF' - PERFORMS ONE PASS THROUGH DATA AS PART
!     OF MULTIPLE REAL FFT (FOURIER SYNTHESIS) ROUTINE
!
!     A IS FIRST REAL INPUT VECTOR
!         EQUIVALENCE B(1) WITH A (LA*INC1+1)
!     C IS FIRST REAL OUTPUT VECTOR
!         EQUIVALENCE D(1) WITH C(IFAC*LA*INC2+1)
!     TRIGS IS A PRECALCULATED LIST OF SINES & COSINES
!     INC1 IS THE ADDRESSING INCREMENT FOR A
!     INC2 IS THE ADDRESSING INCREMENT FOR C
!     INC3 IS THE INCREMENT BETWEEN INPUT VECTORS A
!     INC4 IS THE INCREMENT BETWEEN OUTPUT VECTORS C
!     LOT IS THE NUMBER OF VECTORS
!     N IS THE LENGTH OF THE VECTORS
!     IFAC IS THE CURRENT FACTOR OF N
!     LA IS THE PRODUCT OF PREVIOUS FACTORS
!     IERR IS AN ERROR INDICATOR:
!              0 - PASS COMPLETED WITHOUT ERROR
!              1 - LOT GREATER THAN 64
!              2 - IFAC NOT CATERED FOR
!              3 - IFAC ONLY CATERED FOR IF LA=N/IFAC
!     LIPL=.T. => RESULTS ARE RETURNED TO INPUT ARRAY
!              (ONLY VALID IF LA=N/IFAC, I.E. ON LAST PASS)
!
!-----------------------------------------------------------------------
!
      SUBROUTINE RPASSF(A,B,C,D,TRIGS,INC1,INC2,INC3,INC4,LOT,N,IFAC,   &
     &    LA,IERR,LIPL)
!AUTOPROMOTE
!
      USE PARKIND1, ONLY : JPIM, JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
!
      IMPLICIT NONE
!
      INTEGER(KIND=JPIM) :: N
      REAL(KIND=JPRB) :: A(*)
      REAL(KIND=JPRB) :: B(*)
      REAL(KIND=JPRB) :: C(*)
      REAL(KIND=JPRB) :: D(*)
      REAL(KIND=JPRB) :: TRIGS(N)
      REAL(KIND=JPRB) :: A10,A11,A20,A21
      REAL(KIND=JPRB) :: B10,B11,B20,B21
      REAL(KIND=JPRB) :: C1,C2,C3,C4,C5
      REAL(KIND=JPRB) :: S1,S2,S3,S4,S5
      REAL(KIND=JPRB) :: SIN36,SIN45,SIN60,SIN72
      REAL(KIND=JPRB) :: SSIN36,SSIN45,SSIN60,SSIN72
      REAL(KIND=JPRB) :: QRT5,QQRT5
      REAL(KIND=JPRB) :: T1,T2,T3,T4,T5,T6,T7
      INTEGER(KIND=JPIM) :: IERR
      INTEGER(KIND=JPIM) :: INC1
      INTEGER(KIND=JPIM) :: INC2
      INTEGER(KIND=JPIM) :: INC3
      INTEGER(KIND=JPIM) :: INC4
      INTEGER(KIND=JPIM) :: LOT
      INTEGER(KIND=JPIM) :: IFAC
      INTEGER(KIND=JPIM) :: LA
      INTEGER(KIND=JPIM) :: INC21,IINK,IJK,ILOT,ILA
      INTEGER(KIND=JPIM) :: I,IA,IB,IBAD,IBASE,IC,ID,IE,IF
      INTEGER(KIND=JPIM) :: J,JA,JB,JBASE,JC,JD,JE,JF,JG,JH,JINK,JUMP
      INTEGER(KIND=JPIM) :: K,KB,KC,KD,KE,KF,KSTOP
      INTEGER(KIND=JPIM) :: L,M
      LOGICAL :: LIPL
!
      DATA SIN36/0.587785252292473_JPRB/,SIN72/0.951056516295154_JPRB/, &
     &    QRT5/0.559016994374947_JPRB/,SIN60/0.866025403784437_JPRB/
!
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      M=N/IFAC
      IINK=LA*INC1
      JINK=LA*INC2
      JUMP=(IFAC-1)*JINK
      KSTOP=(N-IFAC)/(2*IFAC)
!
      IBASE=0
      JBASE=0
      IBAD=0
!
!     Increase the vector length by fusing the loops if the
!     data layout is appropriate:
      IF (INC1.EQ.LOT.AND.INC2.EQ.LOT.AND.INC3.EQ.1.AND.INC4.EQ.1) THEN
        ILA=1
        ILOT=LA*LOT
        INC21=LA*LOT
      ELSE
        ILA=LA
        ILOT=LOT
        INC21=INC2
      ENDIF
!
      IF (IFAC.EQ.2) THEN
!
!     CODING FOR FACTOR 2
!     -------------------
  200 CONTINUE
      IA=1
      IB=IA+(2*M-LA)*INC1
      JA=1
      JB=JA+JINK
!
      IF (LA.NE.M) THEN
!
      DO 220 L=1,ILA
      I=IBASE
      J=JBASE
!OCL NOVREC
!DEC$ IVDEP
!NEC$ ivdep
      DO 210 IJK=1,ILOT
      C(JA+J)=A(IA+I)+A(IB+I)
      C(JB+J)=A(IA+I)-A(IB+I)
      I=I+INC3
      J=J+INC4
  210 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC21
  220 CONTINUE
      IA=IA+IINK
      IINK=2*IINK
      IB=IB-IINK
      IBASE=0
      JBASE=JBASE+JUMP
      JUMP=2*JUMP+JINK
!
      IF (IA.LT.IB) THEN
      DO 250 K=LA,KSTOP,LA
      KB=K+K
      C1=TRIGS(KB+1)
      S1=TRIGS(KB+2)
      IBASE=0
      DO 240 L=1,ILA
      I=IBASE
      J=JBASE
!OCL NOVREC
!DEC$ IVDEP
!NEC$ ivdep
      DO 230 IJK=1,ILOT
      C(JA+J)=A(IA+I)+A(IB+I)
      D(JA+J)=B(IA+I)-B(IB+I)
      C(JB+J)=C1*(A(IA+I)-A(IB+I))-S1*(B(IA+I)+B(IB+I))
      D(JB+J)=S1*(A(IA+I)-A(IB+I))+C1*(B(IA+I)+B(IB+I))
      I=I+INC3
      J=J+INC4
  230 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC21
  240 CONTINUE
      IA=IA+IINK
      IB=IB-IINK
      JBASE=JBASE+JUMP
  250 CONTINUE
      ENDIF
!
      IF (IA.EQ.IB) THEN
      IBASE=0
      DO 280 L=1,ILA
      I=IBASE
      J=JBASE
!OCL NOVREC
!DEC$ IVDEP
!NEC$ ivdep
      DO 270 IJK=1,ILOT
      C(JA+J)=A(IA+I)
      C(JB+J)=-B(IA+I)
      I=I+INC3
      J=J+INC4
  270 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC21
  280 CONTINUE
      ENDIF
!
      ELSE                !!! Case LA=M
      IF (LIPL) THEN
        DO 294 L=1,ILA
        I=IBASE
!OCL NOVREC
!NEC$ ivdep
        DO 292 IJK=1,ILOT
        T1=2.0*(A(IA+I)-A(IB+I))
        A(IA+I)=2.0_JPRB*(A(IA+I)+A(IB+I))
        A(IB+I)=T1
        I=I+INC3
  292   CONTINUE
        IBASE=IBASE+INC1
  294   CONTINUE
      ELSE
        DO 298 L=1,ILA
        I=IBASE
        J=JBASE
!OCL NOVREC
!DEC$ IVDEP
!NEC$ ivdep
        DO 296 IJK=1,ILOT
        C(JA+J)=2.0_JPRB*(A(IA+I)+A(IB+I))
        C(JB+J)=2.0_JPRB*(A(IA+I)-A(IB+I))
        I=I+INC3
        J=J+INC4
  296   CONTINUE
        IBASE=IBASE+INC1
        JBASE=JBASE+INC21
  298   CONTINUE
      ENDIF
      ENDIF
!
      ELSEIF (IFAC.EQ.3) THEN
!
!     CODING FOR FACTOR 3
!     -------------------
  300 CONTINUE
      IA=1
      IB=IA+(2*M-LA)*INC1
      IC=IB
      JA=1
      JB=JA+JINK
      JC=JB+JINK
!
      IF (LA.NE.M) THEN
!
      DO 320 L=1,ILA
      I=IBASE
      J=JBASE
!OCL NOVREC
!DEC$ IVDEP
!NEC$ ivdep
      DO 310 IJK=1,ILOT
      C(JA+J)=A(IA+I)+A(IB+I)
      C(JB+J)=(A(IA+I)-0.5_JPRB*A(IB+I))-(SIN60*(B(IB+I)))
      C(JC+J)=(A(IA+I)-0.5_JPRB*A(IB+I))+(SIN60*(B(IB+I)))
      I=I+INC3
      J=J+INC4
  310 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC21
  320 CONTINUE
      IA=IA+IINK
      IINK=2*IINK
      IB=IB+IINK
      IC=IC-IINK
      JBASE=JBASE+JUMP
      JUMP=2*JUMP+JINK
!
      IF (IA.LT.IC) THEN
      DO 350 K=LA,KSTOP,LA
      KB=K+K
      KC=KB+KB
      C1=TRIGS(KB+1)
      S1=TRIGS(KB+2)
      C2=TRIGS(KC+1)
      S2=TRIGS(KC+2)
      IBASE=0
      DO 340 L=1,ILA
      I=IBASE
      J=JBASE
!OCL NOVREC
!DEC$ IVDEP
!NEC$ ivdep
      DO 330 IJK=1,ILOT
      C(JA+J)=A(IA+I)+(A(IB+I)+A(IC+I))
      D(JA+J)=B(IA+I)+(B(IB+I)-B(IC+I))
      C(JB+J)=                                                          &
     &    C1*((A(IA+I)-0.5_JPRB*(A(IB+I)+A(IC+I)))-                     &
     &   (SIN60*(B(IB+I)+B(IC+I))))                                     &
     &   -S1*((B(IA+I)-0.5_JPRB*(B(IB+I)-B(IC+I)))+                     &
     &   (SIN60*(A(IB+I)-A(IC+I))))
      D(JB+J)=                                                          &
     &    S1*((A(IA+I)-0.5_JPRB*(A(IB+I)+A(IC+I)))-                     &
     &   (SIN60*(B(IB+I)+B(IC+I))))                                     &
     &   +C1*((B(IA+I)-0.5_JPRB*(B(IB+I)-B(IC+I)))+                     &
     &   (SIN60*(A(IB+I)-A(IC+I))))
      C(JC+J)=                                                          &
     &    C2*((A(IA+I)-0.5_JPRB*(A(IB+I)+A(IC+I)))+                     &
     &   (SIN60*(B(IB+I)+B(IC+I))))                                     &
     &   -S2*((B(IA+I)-0.5_JPRB*(B(IB+I)-B(IC+I)))-                     &
     &   (SIN60*(A(IB+I)-A(IC+I))))
      D(JC+J)=                                                          &
     &    S2*((A(IA+I)-0.5_JPRB*(A(IB+I)+A(IC+I)))+                     &
     &   (SIN60*(B(IB+I)+B(IC+I))))                                     &
     &   +C2*((B(IA+I)-0.5_JPRB*(B(IB+I)-B(IC+I)))-                     &
     &   (SIN60*(A(IB+I)-A(IC+I))))
      I=I+INC3
      J=J+INC4
  330 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC21
  340 CONTINUE
      IA=IA+IINK
      IB=IB+IINK
      IC=IC-IINK
      JBASE=JBASE+JUMP
  350 CONTINUE
      ENDIF
!
      IF (IA.EQ.IC) THEN
      IBASE=0
      DO 380 L=1,ILA
      I=IBASE
      J=JBASE
!OCL NOVREC
!DEC$ IVDEP
!NEC$ ivdep
      DO 370 IJK=1,ILOT
      C(JA+J)=A(IA+I)+A(IB+I)
      C(JB+J)=(0.5_JPRB*A(IA+I)-A(IB+I))-(SIN60*B(IA+I))
      C(JC+J)=-(0.5_JPRB*A(IA+I)-A(IB+I))-(SIN60*B(IA+I))
      I=I+INC3
      J=J+INC4
  370 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC21
  380 CONTINUE
      ENDIF
!
      ELSE                !!! Case LA=M
      SSIN60=2.0*SIN60
      IF (LIPL) THEN
        DO 394 L=1,ILA
        I=IBASE
!OCL NOVREC
!NEC$ ivdep
        DO 392 IJK=1,ILOT
        T1=(2.0_JPRB*A(IA+I)-A(IB+I))-(SSIN60*B(IB+I))
        T2=(2.0_JPRB*A(IA+I)-A(IB+I))+(SSIN60*B(IB+I))
        A(IA+I)=2.0_JPRB*(A(IA+I)+A(IB+I))
        A(IB+I)=T1
        B(IB+I)=T2
        I=I+INC3
  392   CONTINUE
        IBASE=IBASE+INC1
  394   CONTINUE
      ELSE
        DO 398 L=1,ILA
        I=IBASE
        J=JBASE
!OCL NOVREC
!DEC$ IVDEP
!NEC$ ivdep
        DO 396 IJK=1,ILOT
        C(JA+J)=2.0_JPRB*(A(IA+I)+A(IB+I))
        C(JB+J)=(2.0_JPRB*A(IA+I)-A(IB+I))-(SSIN60*B(IB+I))
        C(JC+J)=(2.0_JPRB*A(IA+I)-A(IB+I))+(SSIN60*B(IB+I))
        I=I+INC3
        J=J+INC4
  396   CONTINUE
        IBASE=IBASE+INC1
        JBASE=JBASE+INC21
  398   CONTINUE
      ENDIF
      ENDIF
!
      ELSEIF (IFAC.EQ.4) THEN
!
!     CODING FOR FACTOR 4
!     -------------------
  400 CONTINUE
      IA=1
      IB=IA+(2*M-LA)*INC1
      IC=IB+2*M*INC1
      ID=IB
      JA=1
      JB=JA+JINK
      JC=JB+JINK
      JD=JC+JINK
!
      IF (LA.NE.M) THEN
!
      DO 420 L=1,ILA
      I=IBASE
      J=JBASE
!OCL NOVREC
!DEC$ IVDEP
!NEC$ ivdep
      DO 410 IJK=1,ILOT
      C(JA+J)=(A(IA+I)+A(IC+I))+A(IB+I)
      C(JB+J)=(A(IA+I)-A(IC+I))-B(IB+I)
      C(JC+J)=(A(IA+I)+A(IC+I))-A(IB+I)
      C(JD+J)=(A(IA+I)-A(IC+I))+B(IB+I)
      I=I+INC3
      J=J+INC4
  410 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC21
  420 CONTINUE
      IA=IA+IINK
      IINK=2*IINK
      IB=IB+IINK
      IC=IC-IINK
      ID=ID-IINK
      JBASE=JBASE+JUMP
      JUMP=2*JUMP+JINK
!
      IF (IB.LT.IC) THEN
      DO 450 K=LA,KSTOP,LA
      KB=K+K
      KC=KB+KB
      KD=KC+KB
      C1=TRIGS(KB+1)
      S1=TRIGS(KB+2)
      C2=TRIGS(KC+1)
      S2=TRIGS(KC+2)
      C3=TRIGS(KD+1)
      S3=TRIGS(KD+2)
      IBASE=0
      DO 440 L=1,ILA
      I=IBASE
      J=JBASE
!OCL NOVREC
!DEC$ IVDEP
!NEC$ ivdep
      DO 430 IJK=1,ILOT
      C(JA+J)=(A(IA+I)+A(IC+I))+(A(IB+I)+A(ID+I))
      D(JA+J)=(B(IA+I)-B(IC+I))+(B(IB+I)-B(ID+I))
      C(JC+J)=                                                          &
     &    C2*((A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I)))                      &
     &   -S2*((B(IA+I)-B(IC+I))-(B(IB+I)-B(ID+I)))
      D(JC+J)=                                                          &
     &    S2*((A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I)))                      &
     &   +C2*((B(IA+I)-B(IC+I))-(B(IB+I)-B(ID+I)))
      C(JB+J)=                                                          &
     &    C1*((A(IA+I)-A(IC+I))-(B(IB+I)+B(ID+I)))                      &
     &   -S1*((B(IA+I)+B(IC+I))+(A(IB+I)-A(ID+I)))
      D(JB+J)=                                                          &
     &    S1*((A(IA+I)-A(IC+I))-(B(IB+I)+B(ID+I)))                      &
     &   +C1*((B(IA+I)+B(IC+I))+(A(IB+I)-A(ID+I)))
      C(JD+J)=                                                          &
     &    C3*((A(IA+I)-A(IC+I))+(B(IB+I)+B(ID+I)))                      &
     &   -S3*((B(IA+I)+B(IC+I))-(A(IB+I)-A(ID+I)))
      D(JD+J)=                                                          &
     &    S3*((A(IA+I)-A(IC+I))+(B(IB+I)+B(ID+I)))                      &
     &   +C3*((B(IA+I)+B(IC+I))-(A(IB+I)-A(ID+I)))
      I=I+INC3
      J=J+INC4
  430 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC21
  440 CONTINUE
      IA=IA+IINK
      IB=IB+IINK
      IC=IC-IINK
      ID=ID-IINK
      JBASE=JBASE+JUMP
  450 CONTINUE
      ENDIF
!
      IF (IB.EQ.IC) THEN
      IBASE=0
      SIN45=SQRT(0.5_JPRB)
      DO 480 L=1,ILA
      I=IBASE
      J=JBASE
!OCL NOVREC
!DEC$ IVDEP
!NEC$ ivdep
      DO 470 IJK=1,ILOT
      C(JA+J)=A(IA+I)+A(IB+I)
      C(JB+J)=SIN45*((A(IA+I)-A(IB+I))-(B(IA+I)+B(IB+I)))
      C(JC+J)=B(IB+I)-B(IA+I)
      C(JD+J)=-SIN45*((A(IA+I)-A(IB+I))+(B(IA+I)+B(IB+I)))
      I=I+INC3
      J=J+INC4
  470 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC21
  480 CONTINUE
      ENDIF
!
      ELSE                !!! Case LA=M
      IF (LIPL) THEN
        DO 494 L=1,ILA
        I=IBASE
!OCL NOVREC
!NEC$ ivdep
        DO 492 IJK=1,ILOT
        T1=2.0_JPRB*((A(IA+I)-A(IC+I))-B(IB+I))
        T2=2.0_JPRB*((A(IA+I)+A(IC+I))-A(IB+I))
        T3=2.0_JPRB*((A(IA+I)-A(IC+I))+B(IB+I))
        A(IA+I)=2.0_JPRB*((A(IA+I)+A(IC+I))+A(IB+I))
        A(IB+I)=T1
        B(IB+I)=T2
        A(IC+I)=T3
        I=I+INC3
  492   CONTINUE
        IBASE=IBASE+INC1
  494   CONTINUE
      ELSE
        DO 498 L=1,ILA
        I=IBASE
        J=JBASE
!OCL NOVREC
!DEC$ IVDEP
!NEC$ ivdep
        DO 496 IJK=1,ILOT
        C(JA+J)=2.0_JPRB*((A(IA+I)+A(IC+I))+A(IB+I))
        C(JB+J)=2.0_JPRB*((A(IA+I)-A(IC+I))-B(IB+I))
        C(JC+J)=2.0_JPRB*((A(IA+I)+A(IC+I))-A(IB+I))
        C(JD+J)=2.0_JPRB*((A(IA+I)-A(IC+I))+B(IB+I))
        I=I+INC3
        J=J+INC4
  496   CONTINUE
        IBASE=IBASE+INC1
        JBASE=JBASE+INC21
  498   CONTINUE
      ENDIF
      ENDIF
!
      ELSEIF (IFAC.EQ.5) THEN
!
!     CODING FOR FACTOR 5
!     -------------------
  500 CONTINUE
      IA=1
      IB=IA+(2*M-LA)*INC1
      IC=IB+2*M*INC1
      ID=IC
      IE=IB
      JA=1
      JB=JA+JINK
      JC=JB+JINK
      JD=JC+JINK
      JE=JD+JINK
!
      IF (LA.NE.M) THEN
!
      DO 520 L=1,ILA
      I=IBASE
      J=JBASE
!OCL NOVREC
!DEC$ IVDEP
!NEC$ ivdep
      DO 510 IJK=1,ILOT
      C(JA+J)=A(IA+I)+(A(IB+I)+A(IC+I))
      C(JB+J)=((A(IA+I)-0.25_JPRB*(A(IB+I)+A(IC+I)))+                   &
     &      QRT5*(A(IB+I)-A(IC+I)))-(SIN72*B(IB+I)+SIN36*B(IC+I))
      C(JC+J)=((A(IA+I)-0.25_JPRB*(A(IB+I)+A(IC+I)))-                   &
     &      QRT5*(A(IB+I)-A(IC+I)))-(SIN36*B(IB+I)-SIN72*B(IC+I))
      C(JD+J)=((A(IA+I)-0.25_JPRB*(A(IB+I)+A(IC+I)))-                   &
     &      QRT5*(A(IB+I)-A(IC+I)))+(SIN36*B(IB+I)-SIN72*B(IC+I))
      C(JE+J)=((A(IA+I)-0.25_JPRB*(A(IB+I)+A(IC+I)))+                   &
     &      QRT5*(A(IB+I)-A(IC+I)))+(SIN72*B(IB+I)+SIN36*B(IC+I))
      I=I+INC3
      J=J+INC4
  510 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC21
  520 CONTINUE
      IA=IA+IINK
      IINK=2*IINK
      IB=IB+IINK
      IC=IC+IINK
      ID=ID-IINK
      IE=IE-IINK
      JBASE=JBASE+JUMP
      JUMP=2*JUMP+JINK
!
      IF (IB.LT.ID) THEN
      DO 550 K=LA,KSTOP,LA
      KB=K+K
      KC=KB+KB
      KD=KC+KB
      KE=KD+KB
      C1=TRIGS(KB+1)
      S1=TRIGS(KB+2)
      C2=TRIGS(KC+1)
      S2=TRIGS(KC+2)
      C3=TRIGS(KD+1)
      S3=TRIGS(KD+2)
      C4=TRIGS(KE+1)
      S4=TRIGS(KE+2)
      IBASE=0
      DO 540 L=1,ILA
      I=IBASE
      J=JBASE
!OCL NOVREC
!DEC$ IVDEP
!NEC$ ivdep
      DO 530 IJK=1,ILOT
!
      A10=(A(IA+I)-0.25_JPRB*((A(IB+I)+A(IE+I))+(A(IC+I)+A(ID+I))))     &
     &    +QRT5*((A(IB+I)+A(IE+I))-(A(IC+I)+A(ID+I)))
      A20=(A(IA+I)-0.25_JPRB*((A(IB+I)+A(IE+I))+(A(IC+I)+A(ID+I))))     &
     &    -QRT5*((A(IB+I)+A(IE+I))-(A(IC+I)+A(ID+I)))
      B10=(B(IA+I)-0.25_JPRB*((B(IB+I)-B(IE+I))+(B(IC+I)-B(ID+I))))     &
     &    +QRT5*((B(IB+I)-B(IE+I))-(B(IC+I)-B(ID+I)))
      B20=(B(IA+I)-0.25_JPRB*((B(IB+I)-B(IE+I))+(B(IC+I)-B(ID+I))))     &
     &    -QRT5*((B(IB+I)-B(IE+I))-(B(IC+I)-B(ID+I)))
      A11=SIN72*(B(IB+I)+B(IE+I))+SIN36*(B(IC+I)+B(ID+I))
      A21=SIN36*(B(IB+I)+B(IE+I))-SIN72*(B(IC+I)+B(ID+I))
      B11=SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I))
      B21=SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I))
!
      C(JA+J)=A(IA+I)+((A(IB+I)+A(IE+I))+(A(IC+I)+A(ID+I)))
      D(JA+J)=B(IA+I)+((B(IB+I)-B(IE+I))+(B(IC+I)-B(ID+I)))
      C(JB+J)=C1*(A10-A11)-S1*(B10+B11)
      D(JB+J)=S1*(A10-A11)+C1*(B10+B11)
      C(JE+J)=C4*(A10+A11)-S4*(B10-B11)
      D(JE+J)=S4*(A10+A11)+C4*(B10-B11)
      C(JC+J)=C2*(A20-A21)-S2*(B20+B21)
      D(JC+J)=S2*(A20-A21)+C2*(B20+B21)
      C(JD+J)=C3*(A20+A21)-S3*(B20-B21)
      D(JD+J)=S3*(A20+A21)+C3*(B20-B21)
!
      I=I+INC3
      J=J+INC4
  530 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC21
  540 CONTINUE
      IA=IA+IINK
      IB=IB+IINK
      IC=IC+IINK
      ID=ID-IINK
      IE=IE-IINK
      JBASE=JBASE+JUMP
  550 CONTINUE
      ENDIF
!
      IF (IB.EQ.ID) THEN
      IBASE=0
      DO 580 L=1,ILA
      I=IBASE
      J=JBASE
!OCL NOVREC
!DEC$ IVDEP
!NEC$ ivdep
      DO 570 IJK=1,ILOT
      C(JA+J)=(A(IA+I)+A(IB+I))+A(IC+I)
      C(JB+J)=(QRT5*(A(IA+I)-A(IB+I))+                                  &
     &     (0.25_JPRB*(A(IA+I)+A(IB+I))-A(IC+I)))                       &
     &    -(SIN36*B(IA+I)+SIN72*B(IB+I))
      C(JE+J)=-(QRT5*(A(IA+I)-A(IB+I))+                                 &
     &     (0.25_JPRB*(A(IA+I)+A(IB+I))-A(IC+I)))                       &
     &    -(SIN36*B(IA+I)+SIN72*B(IB+I))
      C(JC+J)=(QRT5*(A(IA+I)-A(IB+I))-                                  &
     &     (0.25_JPRB*(A(IA+I)+A(IB+I))-A(IC+I)))                       &
     &    -(SIN72*B(IA+I)-SIN36*B(IB+I))
      C(JD+J)=-(QRT5*(A(IA+I)-A(IB+I))-                                 &
     &     (0.25_JPRB*(A(IA+I)+A(IB+I))-A(IC+I)))                       &
     &    -(SIN72*B(IA+I)-SIN36*B(IB+I))
      I=I+INC3
      J=J+INC4
  570 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC21
  580 CONTINUE
      ENDIF
!
      ELSE               !!! Case LA=M
      QQRT5=2.0*QRT5
      SSIN36=2.0*SIN36
      SSIN72=2.0*SIN72
      IF (LIPL) THEN
        DO 594 L=1,ILA
        I=IBASE
!OCL NOVREC
!NEC$ ivdep
        DO 592 IJK=1,ILOT
        T1=(2.0_JPRB*(A(IA+I)-0.25_JPRB*(A(IB+I)+A(IC+I)))              &
     &    +QQRT5*(A(IB+I)-A(IC+I)))-(SSIN72*B(IB+I)+SSIN36*B(IC+I))
        T2=(2.0_JPRB*(A(IA+I)-0.25_JPRB*(A(IB+I)+A(IC+I)))              &
     &    -QQRT5*(A(IB+I)-A(IC+I)))-(SSIN36*B(IB+I)-SSIN72*B(IC+I))
        T3=(2.0_JPRB*(A(IA+I)-0.25_JPRB*(A(IB+I)+A(IC+I)))              &
     &    -QQRT5*(A(IB+I)-A(IC+I)))+(SSIN36*B(IB+I)-SSIN72*B(IC+I))
        T4=(2.0_JPRB*(A(IA+I)-0.25_JPRB*(A(IB+I)+A(IC+I)))              &
     &    +QQRT5*(A(IB+I)-A(IC+I)))+(SSIN72*B(IB+I)+SSIN36*B(IC+I))
        A(IA+I)=2.0_JPRB*(A(IA+I)+(A(IB+I)+A(IC+I)))
        A(IB+I)=T1
        B(IB+I)=T2
        A(IC+I)=T3
        B(IC+I)=T4
        I=I+INC3
  592   CONTINUE
        IBASE=IBASE+INC1
  594   CONTINUE
      ELSE
        DO 598 L=1,ILA
        I=IBASE
        J=JBASE
!OCL NOVREC
!DEC$ IVDEP
!NEC$ ivdep
        DO 596 IJK=1,ILOT
        C(JA+J)=2.0_JPRB*(A(IA+I)+(A(IB+I)+A(IC+I)))
        C(JB+J)=(2.0_JPRB*(A(IA+I)-0.25_JPRB*(A(IB+I)+A(IC+I)))         &
     &    +QQRT5*(A(IB+I)-A(IC+I)))-(SSIN72*B(IB+I)+SSIN36*B(IC+I))
        C(JC+J)=(2.0_JPRB*(A(IA+I)-0.25_JPRB*(A(IB+I)+A(IC+I)))         &
     &    -QQRT5*(A(IB+I)-A(IC+I)))-(SSIN36*B(IB+I)-SSIN72*B(IC+I))
        C(JD+J)=(2.0_JPRB*(A(IA+I)-0.25_JPRB*(A(IB+I)+A(IC+I)))         &
     &    -QQRT5*(A(IB+I)-A(IC+I)))+(SSIN36*B(IB+I)-SSIN72*B(IC+I))
        C(JE+J)=(2.0_JPRB*(A(IA+I)-0.25_JPRB*(A(IB+I)+A(IC+I)))         &
     &    +QQRT5*(A(IB+I)-A(IC+I)))+(SSIN72*B(IB+I)+SSIN36*B(IC+I))
        I=I+INC3
        J=J+INC4
  596   CONTINUE
        IBASE=IBASE+INC1
        JBASE=JBASE+INC21
  598   CONTINUE
      ENDIF
      ENDIF
!
      ELSEIF (IFAC.EQ.6) THEN
!
!     CODING FOR FACTOR 6
!     -------------------
  600 CONTINUE
      IA=1
      IB=IA+(2*M-LA)*INC1
      IC=IB+2*M*INC1
      ID=IC+2*M*INC1
      IE=IC
      IF=IB
      JA=1
      JB=JA+JINK
      JC=JB+JINK
      JD=JC+JINK
      JE=JD+JINK
      JF=JE+JINK
!
      IF (LA.NE.M) THEN
!
      DO 620 L=1,ILA
      I=IBASE
      J=JBASE
!OCL NOVREC
!DEC$ IVDEP
!NEC$ ivdep
      DO 610 IJK=1,ILOT
      C(JA+J)=(A(IA+I)+A(ID+I))+(A(IB+I)+A(IC+I))
      C(JD+J)=(A(IA+I)-A(ID+I))-(A(IB+I)-A(IC+I))
      C(JB+J)=((A(IA+I)-A(ID+I))+0.5_JPRB*(A(IB+I)-A(IC+I)))            &
     &    -(SIN60*(B(IB+I)+B(IC+I)))
      C(JF+J)=((A(IA+I)-A(ID+I))+0.5_JPRB*(A(IB+I)-A(IC+I)))            &
     &    +(SIN60*(B(IB+I)+B(IC+I)))
      C(JC+J)=((A(IA+I)+A(ID+I))-0.5_JPRB*(A(IB+I)+A(IC+I)))            &
     &    -(SIN60*(B(IB+I)-B(IC+I)))
      C(JE+J)=((A(IA+I)+A(ID+I))-0.5_JPRB*(A(IB+I)+A(IC+I)))            &
     &    +(SIN60*(B(IB+I)-B(IC+I)))
      I=I+INC3
      J=J+INC4
  610 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC21
  620 CONTINUE
      IA=IA+IINK
      IINK=2*IINK
      IB=IB+IINK
      IC=IC+IINK
      ID=ID-IINK
      IE=IE-IINK
      IF=IF-IINK
      JBASE=JBASE+JUMP
      JUMP=2*JUMP+JINK
!
      IF (IC.LT.ID) THEN
      DO 650 K=LA,KSTOP,LA
      KB=K+K
      KC=KB+KB
      KD=KC+KB
      KE=KD+KB
      KF=KE+KB
      C1=TRIGS(KB+1)
      S1=TRIGS(KB+2)
      C2=TRIGS(KC+1)
      S2=TRIGS(KC+2)
      C3=TRIGS(KD+1)
      S3=TRIGS(KD+2)
      C4=TRIGS(KE+1)
      S4=TRIGS(KE+2)
      C5=TRIGS(KF+1)
      S5=TRIGS(KF+2)
      IBASE=0
      DO 640 L=1,ILA
      I=IBASE
      J=JBASE
!OCL NOVREC
!DEC$ IVDEP
!NEC$ ivdep
      DO 630 IJK=1,ILOT
!
      A11= (A(IE+I)+A(IB+I))+(A(IC+I)+A(IF+I))
      A20=(A(IA+I)+A(ID+I))-0.5_JPRB*A11
      A21=SIN60*((A(IE+I)+A(IB+I))-(A(IC+I)+A(IF+I)))
      B11= (B(IB+I)-B(IE+I))+(B(IC+I)-B(IF+I))
      B20=(B(IA+I)-B(ID+I))-0.5_JPRB*B11
      B21=SIN60*((B(IB+I)-B(IE+I))-(B(IC+I)-B(IF+I)))
!
      C(JA+J)=(A(IA+I)+A(ID+I))+A11
      D(JA+J)=(B(IA+I)-B(ID+I))+B11
      C(JC+J)=C2*(A20-B21)-S2*(B20+A21)
      D(JC+J)=S2*(A20-B21)+C2*(B20+A21)
      C(JE+J)=C4*(A20+B21)-S4*(B20-A21)
      D(JE+J)=S4*(A20+B21)+C4*(B20-A21)
!
      A11=(A(IE+I)-A(IB+I))+(A(IC+I)-A(IF+I))
      B11=(B(IE+I)+B(IB+I))-(B(IC+I)+B(IF+I))
      A20=(A(IA+I)-A(ID+I))-0.5_JPRB*A11
      A21=SIN60*((A(IE+I)-A(IB+I))-(A(IC+I)-A(IF+I)))
      B20=(B(IA+I)+B(ID+I))+0.5_JPRB*B11
      B21=SIN60*((B(IE+I)+B(IB+I))+(B(IC+I)+B(IF+I)))
!
      C(JD+J)=                                                          &
     &  C3*((A(IA+I)-A(ID+I))+A11)-S3*((B(IA+I)+B(ID+I))-B11)
      D(JD+J)=                                                          &
     &  S3*((A(IA+I)-A(ID+I))+A11)+C3*((B(IA+I)+B(ID+I))-B11)
      C(JB+J)=C1*(A20-B21)-S1*(B20-A21)
      D(JB+J)=S1*(A20-B21)+C1*(B20-A21)
      C(JF+J)=C5*(A20+B21)-S5*(B20+A21)
      D(JF+J)=S5*(A20+B21)+C5*(B20+A21)
!
      I=I+INC3
      J=J+INC4
  630 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC21
  640 CONTINUE
      IA=IA+IINK
      IB=IB+IINK
      IC=IC+IINK
      ID=ID-IINK
      IE=IE-IINK
      IF=IF-IINK
      JBASE=JBASE+JUMP
  650 CONTINUE
      ENDIF
!
      IF (IC.EQ.ID) THEN
      IBASE=0
      DO 680 L=1,ILA
      I=IBASE
      J=JBASE
!OCL NOVREC
!DEC$ IVDEP
!NEC$ ivdep
      DO 670 IJK=1,ILOT
      C(JA+J)=A(IB+I)+(A(IA+I)+A(IC+I))
      C(JD+J)=B(IB+I)-(B(IA+I)+B(IC+I))
      C(JB+J)=(SIN60*(A(IA+I)-A(IC+I)))-                                &
     &        (0.5_JPRB*(B(IA+I)+B(IC+I))+B(IB+I))
      C(JF+J)=-(SIN60*(A(IA+I)-A(IC+I)))-                               &
     &        (0.5_JPRB*(B(IA+I)+B(IC+I))+B(IB+I))
      C(JC+J)=SIN60*(B(IC+I)-B(IA+I))+                                  &
     &        (0.5_JPRB*(A(IA+I)+A(IC+I))-A(IB+I))
      C(JE+J)=SIN60*(B(IC+I)-B(IA+I))-                                  &
     &        (0.5_JPRB*(A(IA+I)+A(IC+I))-A(IB+I))
      I=I+INC3
      J=J+INC4
  670 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC21
  680 CONTINUE
      ENDIF
!
      ELSE                 !!! Case LA=M
      SSIN60=2.0_JPRB*SIN60
      IF (LIPL) THEN
        DO 694 L=1,ILA
        I=IBASE
!OCL NOVREC
!NEC$ ivdep
        DO 692 IJK=1,ILOT
        T1=(2.0_JPRB*(A(IA+I)-A(ID+I))+(A(IB+I)-A(IC+I)))               &
     &    -(SSIN60*(B(IB+I)+B(IC+I)))
        T5=(2.0_JPRB*(A(IA+I)-A(ID+I))+(A(IB+I)-A(IC+I)))               &
     &    +(SSIN60*(B(IB+I)+B(IC+I)))
        T2=(2.0_JPRB*(A(IA+I)+A(ID+I))-(A(IB+I)+A(IC+I)))               &
     &    -(SSIN60*(B(IB+I)-B(IC+I)))
        T4=(2.0_JPRB*(A(IA+I)+A(ID+I))-(A(IB+I)+A(IC+I)))               &
     &    +(SSIN60*(B(IB+I)-B(IC+I)))
        T3=(2.0_JPRB*(A(IA+I)-A(ID+I)))-(2.0_JPRB*(A(IB+I)-A(IC+I)))
        A(IA+I)=(2.0_JPRB*(A(IA+I)+A(ID+I)))+                           &
     &          (2.0_JPRB*(A(IB+I)+A(IC+I)))
        A(IB+I)=T1
        B(IB+I)=T2
        A(IC+I)=T3
        B(IC+I)=T4
        A(ID+I)=T5
        I=I+INC3
  692   CONTINUE
        IBASE=IBASE+INC1
  694   CONTINUE
      ELSE
        DO 698 L=1,ILA
        I=IBASE
        J=JBASE
!OCL NOVREC
!DEC$ IVDEP
!NEC$ ivdep
        DO 696 IJK=1,ILOT
        C(JA+J)=(2.0_JPRB*(A(IA+I)+A(ID+I)))+                           &
     &          (2.0_JPRB*(A(IB+I)+A(IC+I)))
        C(JD+J)=(2.0_JPRB*(A(IA+I)-A(ID+I)))-                           &
     &          (2.0_JPRB*(A(IB+I)-A(IC+I)))
        C(JB+J)=(2.0_JPRB*(A(IA+I)-A(ID+I))+(A(IB+I)-A(IC+I)))          &
     &    -(SSIN60*(B(IB+I)+B(IC+I)))
        C(JF+J)=(2.0_JPRB*(A(IA+I)-A(ID+I))+(A(IB+I)-A(IC+I)))          &
     &    +(SSIN60*(B(IB+I)+B(IC+I)))
        C(JC+J)=(2.0_JPRB*(A(IA+I)+A(ID+I))-(A(IB+I)+A(IC+I)))          &
     &    -(SSIN60*(B(IB+I)-B(IC+I)))
        C(JE+J)=(2.0_JPRB*(A(IA+I)+A(ID+I))-(A(IB+I)+A(IC+I)))          &
     &    +(SSIN60*(B(IB+I)-B(IC+I)))
        I=I+INC3
        J=J+INC4
  696   CONTINUE
        IBASE=IBASE+INC1
        JBASE=JBASE+INC21
  698   CONTINUE
      ENDIF
      ENDIF
!
      ELSEIF (IFAC.EQ.8) THEN
!
!     CODING FOR FACTOR 8
!     -------------------
  800 CONTINUE
      IF (LA.NE.M) THEN
        IBAD=3
      ELSE
      IA=1
      IB=IA+LA*INC1
      IC=IB+2*LA*INC1
      ID=IC+2*LA*INC1
      IE=ID+2*LA*INC1
      JA=1
      JB=JA+JINK
      JC=JB+JINK
      JD=JC+JINK
      JE=JD+JINK
      JF=JE+JINK
      JG=JF+JINK
      JH=JG+JINK
      SSIN45=SQRT(2.0_JPRB)
!
      IF (LIPL) THEN
        DO 820 L=1,ILA
        I=IBASE
!OCL NOVREC
!NEC$ ivdep
        DO 810 IJK=1,ILOT
        T2=2.0_JPRB*(((A(IA+I)+A(IE+I))-A(IC+I))-(B(IB+I)-B(ID+I)))
        T6=2.0_JPRB*(((A(IA+I)+A(IE+I))-A(IC+I))+(B(IB+I)-B(ID+I)))
        T1=2.0_JPRB*((A(IA+I)-A(IE+I))-B(IC+I))                         &
     &    +SSIN45*((A(IB+I)-A(ID+I))-(B(IB+I)+B(ID+I)))
        T5=2.0_JPRB*((A(IA+I)-A(IE+I))-B(IC+I))                         &
     &    -SSIN45*((A(IB+I)-A(ID+I))-(B(IB+I)+B(ID+I)))
        T3=2.0_JPRB*((A(IA+I)-A(IE+I))+B(IC+I))                         &
     &    -SSIN45*((A(IB+I)-A(ID+I))+(B(IB+I)+B(ID+I)))
        T7=2.0_JPRB*((A(IA+I)-A(IE+I))+B(IC+I))                         &
     &    +SSIN45*((A(IB+I)-A(ID+I))+(B(IB+I)+B(ID+I)))
        T4=2.0_JPRB*(((A(IA+I)+A(IE+I))+A(IC+I))-(A(IB+I)+A(ID+I)))
        A(IA+I)=2.0_JPRB*(((A(IA+I)+A(IE+I))+A(IC+I))+(A(IB+I)+A(ID+I)))
        A(IB+I)=T1
        B(IB+I)=T2
        A(IC+I)=T3
        B(IC+I)=T4
        A(ID+I)=T5
        B(ID+I)=T6
        A(IE+I)=T7
        I=I+INC3
  810   CONTINUE
        IBASE=IBASE+INC1
  820   CONTINUE
      ELSE
        DO 840 L=1,ILA
        I=IBASE
        J=JBASE
!OCL NOVREC
!DEC$ IVDEP
!NEC$ ivdep
        DO 830 IJK=1,ILOT
        C(JA+J)=2.0_JPRB*(((A(IA+I)+A(IE+I))+A(IC+I))+(A(IB+I)+A(ID+I)))
        C(JE+J)=2.0_JPRB*(((A(IA+I)+A(IE+I))+A(IC+I))-(A(IB+I)+A(ID+I)))
        C(JC+J)=2.0_JPRB*(((A(IA+I)+A(IE+I))-A(IC+I))-(B(IB+I)-B(ID+I)))
        C(JG+J)=2.0_JPRB*(((A(IA+I)+A(IE+I))-A(IC+I))+(B(IB+I)-B(ID+I)))
        C(JB+J)=2.0_JPRB*((A(IA+I)-A(IE+I))-B(IC+I))                    &
     &    +SSIN45*((A(IB+I)-A(ID+I))-(B(IB+I)+B(ID+I)))
        C(JF+J)=2.0_JPRB*((A(IA+I)-A(IE+I))-B(IC+I))                    &
     &    -SSIN45*((A(IB+I)-A(ID+I))-(B(IB+I)+B(ID+I)))
        C(JD+J)=2.0_JPRB*((A(IA+I)-A(IE+I))+B(IC+I))                    &
     &    -SSIN45*((A(IB+I)-A(ID+I))+(B(IB+I)+B(ID+I)))
        C(JH+J)=2.0_JPRB*((A(IA+I)-A(IE+I))+B(IC+I))                    &
     &    +SSIN45*((A(IB+I)-A(ID+I))+(B(IB+I)+B(ID+I)))
        I=I+INC3
        J=J+INC4
  830   CONTINUE
        IBASE=IBASE+INC1
        JBASE=JBASE+INC21
  840   CONTINUE
      ENDIF
!
      ENDIF
!
      ELSE
!
      IBAD=2       !!! Illegal factor
!
      ENDIF
!
!     RETURN
!     ------
  900 CONTINUE
      IERR=IBAD
      ENDSUBROUTINE RPASSF

!     SUBROUTINE 'QPASSF' - PERFORMS ONE PASS THROUGH DATA AS PART
!     OF MULTIPLE REAL FFT (FOURIER ANALYSIS) ROUTINE
!
!     A IS FIRST REAL INPUT VECTOR
!         EQUIVALENCE B(1) WITH A(IFAC*LA*INC1+1)
!     C IS FIRST REAL OUTPUT VECTOR
!         EQUIVALENCE D(1) WITH C(LA*INC2+1)
!     TRIGS IS A PRECALCULATED LIST OF SINES & COSINES
!     INC1 IS THE ADDRESSING INCREMENT FOR A
!     INC2 IS THE ADDRESSING INCREMENT FOR C
!     INC3 IS THE INCREMENT BETWEEN INPUT VECTORS A
!     INC4 IS THE INCREMENT BETWEEN OUTPUT VECTORS C
!     LOT IS THE NUMBER OF VECTORS
!     N IS THE LENGTH OF THE VECTORS
!     IFAC IS THE CURRENT FACTOR OF N
!     LA = N/(PRODUCT OF FACTORS USED SO FAR)
!     IERR IS AN ERROR INDICATOR:
!              0 - PASS COMPLETED WITHOUT ERROR
!              1 - LOT GREATER THAN 64
!              2 - IFAC NOT CATERED FOR
!              3 - IFAC ONLY CATERED FOR IF LA=N/IFAC
!     LIPL=.T. => RESULTS ARE RETURNED TO INPUT ARRAY
!              (ONLY VALID IF LA=N/IFAC, I.E. ON FIRST PASS)
!
!-----------------------------------------------------------------------
!
      SUBROUTINE QPASSF(A,B,C,D,TRIGS,INC1,INC2,INC3,INC4,LOT,N,IFAC,   &
     &    LA,IERR,LIPL)
!AUTOPROMOTE
      USE PARKIND1, ONLY : JPIM, JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
!
      IMPLICIT NONE
!
      INTEGER(KIND=JPIM) :: N
      REAL(KIND=JPRB) :: A(*)
      REAL(KIND=JPRB) :: B(*)
      REAL(KIND=JPRB) :: C(*)
      REAL(KIND=JPRB) :: D(*)
      REAL(KIND=JPRB) :: TRIGS(N)
      REAL(KIND=JPRB) :: A0,A1,A2,A3,A4,A5,A6,A10,A11,A20,A21
      REAL(KIND=JPRB) :: B0,B1,B2,B3,B4,B5,B6,B10,B11,B20,B21
      REAL(KIND=JPRB) :: C1,C2,C3,C4,C5
      REAL(KIND=JPRB) :: S1,S2,S3,S4,S5
      REAL(KIND=JPRB) :: T1,T2,T3,T4,T5,T6,T7
      REAL(KIND=JPRB) :: Z
      REAL(KIND=JPRB) :: QRT5,SIN36,SIN45,SIN60,SIN72
      REAL(KIND=JPRB) :: ZQRT5,ZSIN36,ZSIN45,ZSIN60,ZSIN72
      INTEGER(KIND=JPIM) :: IERR
      INTEGER(KIND=JPIM) :: INC1
      INTEGER(KIND=JPIM) :: INC2
      INTEGER(KIND=JPIM) :: INC3
      INTEGER(KIND=JPIM) :: INC4
      INTEGER(KIND=JPIM) :: LOT
      INTEGER(KIND=JPIM) :: IFAC
      INTEGER(KIND=JPIM) :: LA
      INTEGER(KIND=JPIM) :: IINK,IJK,ILOT
      INTEGER(KIND=JPIM) :: I,IA,IB,IBAD,IBASE,IC,ID,IE,IF,IG,IH
      INTEGER(KIND=JPIM) :: IJUMP,ILA,INC11
      INTEGER(KIND=JPIM) :: J,JA,JB,JC,JD,JE,JBASE,JF,JINK
      INTEGER(KIND=JPIM) :: K,KB,KC,KD,KE,KF,KSTOP
      INTEGER(KIND=JPIM) :: L,M
      LOGICAL :: LIPL
!
      DATA SIN36/0.587785252292473_JPRB/,SIN72/0.951056516295154_JPRB/, &
     &    QRT5/0.559016994374947_JPRB/,SIN60/0.866025403784437_JPRB/
!
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      M=N/IFAC
      IINK=LA*INC1
      JINK=LA*INC2
      IJUMP=(IFAC-1)*IINK
      KSTOP=(N-IFAC)/(2*IFAC)
!
      IBASE=0
      JBASE=0
      IBAD=0
!
!     Increase the vector length by fusing the loops if the
!     data layout is appropriate:
      IF (INC1.EQ.LOT.AND.INC2.EQ.LOT.AND.INC3.EQ.1.AND.INC4.EQ.1) THEN
        ILA=1
        ILOT=LA*LOT
        INC11=LA*LOT
      ELSE
        ILA=LA
        ILOT=LOT
        INC11=INC1
      ENDIF

!
      IF (IFAC.EQ.2) THEN
!
!     CODING FOR FACTOR 2
!     -------------------
  200 CONTINUE
      IA=1
      IB=IA+IINK
      JA=1
      JB=JA+(2*M-LA)*INC2
!
      IF (LA.NE.M) THEN
!
      DO 220 L=1,ILA
      I=IBASE
      J=JBASE
!OCL NOVREC
!DEC$ IVDEP
!NEC$ ivdep
      DO 210 IJK=1,ILOT
      C(JA+J)=A(IA+I)+A(IB+I)
      C(JB+J)=A(IA+I)-A(IB+I)
      I=I+INC3
      J=J+INC4
  210 CONTINUE
      IBASE=IBASE+INC11
      JBASE=JBASE+INC2
  220 CONTINUE
      JA=JA+JINK
      JINK=2*JINK
      JB=JB-JINK
      IBASE=IBASE+IJUMP
      IJUMP=2*IJUMP+IINK
!
      IF (JA.LT.JB) THEN
      DO 250 K=LA,KSTOP,LA
      KB=K+K
      C1=TRIGS(KB+1)
      S1=TRIGS(KB+2)
      JBASE=0
      DO 240 L=1,ILA
      I=IBASE
      J=JBASE
!OCL NOVREC
!DEC$ IVDEP
!NEC$ ivdep
      DO 230 IJK=1,ILOT
      C(JA+J)=A(IA+I)+(C1*A(IB+I)+S1*B(IB+I))
      C(JB+J)=A(IA+I)-(C1*A(IB+I)+S1*B(IB+I))
      D(JA+J)=(C1*B(IB+I)-S1*A(IB+I))+B(IA+I)
      D(JB+J)=(C1*B(IB+I)-S1*A(IB+I))-B(IA+I)
      I=I+INC3
      J=J+INC4
  230 CONTINUE
      IBASE=IBASE+INC11
      JBASE=JBASE+INC2
  240 CONTINUE
      IBASE=IBASE+IJUMP
      JA=JA+JINK
      JB=JB-JINK
  250 CONTINUE
      ENDIF
!
      IF (JA.EQ.JB) THEN
      JBASE=0
      DO 280 L=1,ILA
      I=IBASE
      J=JBASE
!OCL NOVREC
!DEC$ IVDEP
!NEC$ ivdep
      DO 270 IJK=1,ILOT
      C(JA+J)=A(IA+I)
      D(JA+J)=-A(IB+I)
      I=I+INC3
      J=J+INC4
  270 CONTINUE
      IBASE=IBASE+INC11
      JBASE=JBASE+INC2
  280 CONTINUE
      ENDIF
!
      ELSE                !!! Case LA=M
      Z=1.0_JPRB/REAL(N,KIND=JPRB)
      IF (LIPL) THEN
        DO 294 L=1,ILA
        I=IBASE
!OCL NOVREC
!DEC$ IVDEP
!NEC$ ivdep
        DO 292 IJK=1,ILOT
        T1=Z*(A(IA+I)-A(IB+I))
        A(IA+I)=Z*(A(IA+I)+A(IB+I))
        A(IB+I)=T1
        I=I+INC3
  292   CONTINUE
        IBASE=IBASE+INC11
  294   CONTINUE
      ELSE
        DO 298 L=1,ILA
        I=IBASE
        J=JBASE
!OCL NOVREC
!DEC$ IVDEP
!NEC$ ivdep
        DO 296 IJK=1,ILOT
        C(JA+J)=Z*(A(IA+I)+A(IB+I))
        C(JB+J)=Z*(A(IA+I)-A(IB+I))
        I=I+INC3
        J=J+INC4
  296   CONTINUE
        IBASE=IBASE+INC11
        JBASE=JBASE+INC2
  298   CONTINUE
      ENDIF
      ENDIF
!
      ELSEIF (IFAC.EQ.3) THEN
!
!     CODING FOR FACTOR 3
!     -------------------
  300 CONTINUE
      IA=1
      IB=IA+IINK
      IC=IB+IINK
      JA=1
      JB=JA+(2*M-LA)*INC2
      JC=JB
!
      IF (LA.NE.M) THEN
!
      DO 320 L=1,ILA
      I=IBASE
      J=JBASE
!OCL NOVREC
!DEC$ IVDEP
!NEC$ ivdep
      DO 310 IJK=1,ILOT
      C(JA+J)=A(IA+I)+(A(IB+I)+A(IC+I))
      C(JB+J)=A(IA+I)-0.5_JPRB*(A(IB+I)+A(IC+I))
      D(JB+J)=SIN60*(A(IC+I)-A(IB+I))
      I=I+INC3
      J=J+INC4
  310 CONTINUE
      IBASE=IBASE+INC11
      JBASE=JBASE+INC2
  320 CONTINUE
      JA=JA+JINK
      JINK=2*JINK
      JB=JB+JINK
      JC=JC-JINK
      IBASE=IBASE+IJUMP
      IJUMP=2*IJUMP+IINK
!
      IF (JA.LT.JC) THEN
      DO 350 K=LA,KSTOP,LA
      KB=K+K
      KC=KB+KB
      C1=TRIGS(KB+1)
      S1=TRIGS(KB+2)
      C2=TRIGS(KC+1)
      S2=TRIGS(KC+2)
      JBASE=0
      DO 340 L=1,ILA
      I=IBASE
      J=JBASE
!OCL NOVREC
!DEC$ IVDEP
!NEC$ ivdep
      DO 330 IJK=1,ILOT
      A1=(C1*A(IB+I)+S1*B(IB+I))+(C2*A(IC+I)+S2*B(IC+I))
      B1=(C1*B(IB+I)-S1*A(IB+I))+(C2*B(IC+I)-S2*A(IC+I))
      A2=A(IA+I)-0.5_JPRB*A1
      B2=B(IA+I)-0.5_JPRB*B1
      A3=SIN60*((C1*A(IB+I)+S1*B(IB+I))-(C2*A(IC+I)+S2*B(IC+I)))
      B3=SIN60*((C1*B(IB+I)-S1*A(IB+I))-(C2*B(IC+I)-S2*A(IC+I)))
      C(JA+J)=A(IA+I)+A1
      D(JA+J)=B(IA+I)+B1
      C(JB+J)=A2+B3
      D(JB+J)=B2-A3
      C(JC+J)=A2-B3
      D(JC+J)=-(B2+A3)
      I=I+INC3
      J=J+INC4
  330 CONTINUE
      IBASE=IBASE+INC11
      JBASE=JBASE+INC2
  340 CONTINUE
      IBASE=IBASE+IJUMP
      JA=JA+JINK
      JB=JB+JINK
      JC=JC-JINK
  350 CONTINUE
      ENDIF
!
      IF (JA.EQ.JC) THEN
      JBASE=0
      DO 380 L=1,ILA
      I=IBASE
      J=JBASE
!OCL NOVREC
!DEC$ IVDEP
!NEC$ ivdep
      DO 370 IJK=1,ILOT
      C(JA+J)=A(IA+I)+0.5_JPRB*(A(IB+I)-A(IC+I))
      D(JA+J)=-SIN60*(A(IB+I)+A(IC+I))
      C(JB+J)=A(IA+I)-(A(IB+I)-A(IC+I))
      I=I+INC3
      J=J+INC4
  370 CONTINUE
      IBASE=IBASE+INC11
      JBASE=JBASE+INC2
  380 CONTINUE
      ENDIF
!
      ELSE                 !!! Case LA=M
      Z=1.0_JPRB/REAL(N,KIND=JPRB)
      ZSIN60=Z*SIN60
      IF (LIPL) THEN
        DO 394 L=1,ILA
        I=IBASE
!OCL NOVREC
!DEC$ IVDEP
!NEC$ ivdep
        DO 392 IJK=1,ILOT
        T1=Z*(A(IA+I)-0.5_JPRB*(A(IB+I)+A(IC+I)))
        T2=ZSIN60*(A(IC+I)-A(IB+I))
        A(IA+I)=Z*(A(IA+I)+(A(IB+I)+A(IC+I)))
        A(IB+I)=T1
        A(IC+I)=T2
        I=I+INC3
  392   CONTINUE
        IBASE=IBASE+INC11
  394   CONTINUE
      ELSE
        DO 398 L=1,ILA
        I=IBASE
        J=JBASE
!OCL NOVREC
!DEC$ IVDEP
!NEC$ ivdep
        DO 396 IJK=1,ILOT
        C(JA+J)=Z*(A(IA+I)+(A(IB+I)+A(IC+I)))
        C(JB+J)=Z*(A(IA+I)-0.5_JPRB*(A(IB+I)+A(IC+I)))
        D(JB+J)=ZSIN60*(A(IC+I)-A(IB+I))
        I=I+INC3
        J=J+INC4
  396   CONTINUE
        IBASE=IBASE+INC11
        JBASE=JBASE+INC2
  398   CONTINUE
      ENDIF
      ENDIF
!
      ELSEIF (IFAC.EQ.4) THEN
!
!     CODING FOR FACTOR 4
!     -------------------
  400 CONTINUE
      IA=1
      IB=IA+IINK
      IC=IB+IINK
      ID=IC+IINK
      JA=1
      JB=JA+(2*M-LA)*INC2
      JC=JB+2*M*INC2
      JD=JB
!
      IF (LA.NE.M) THEN
!
      DO 420 L=1,ILA
      I=IBASE
      J=JBASE
!OCL NOVREC
!DEC$ IVDEP
!NEC$ ivdep
      DO 410 IJK=1,ILOT
      C(JA+J)=(A(IA+I)+A(IC+I))+(A(IB+I)+A(ID+I))
      C(JC+J)=(A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I))
      C(JB+J)=A(IA+I)-A(IC+I)
      D(JB+J)=A(ID+I)-A(IB+I)
      I=I+INC3
      J=J+INC4
  410 CONTINUE
      IBASE=IBASE+INC11
      JBASE=JBASE+INC2
  420 CONTINUE
      JA=JA+JINK
      JINK=2*JINK
      JB=JB+JINK
      JC=JC-JINK
      JD=JD-JINK
      IBASE=IBASE+IJUMP
      IJUMP=2*IJUMP+IINK
!
      IF (JB.LT.JC) THEN
      DO 450 K=LA,KSTOP,LA
      KB=K+K
      KC=KB+KB
      KD=KC+KB
      C1=TRIGS(KB+1)
      S1=TRIGS(KB+2)
      C2=TRIGS(KC+1)
      S2=TRIGS(KC+2)
      C3=TRIGS(KD+1)
      S3=TRIGS(KD+2)
      JBASE=0
      DO 440 L=1,ILA
      I=IBASE
      J=JBASE
!OCL NOVREC
!DEC$ IVDEP
!NEC$ ivdep
      DO 430 IJK=1,ILOT
      A0=A(IA+I)+(C2*A(IC+I)+S2*B(IC+I))
      A2=A(IA+I)-(C2*A(IC+I)+S2*B(IC+I))
      A1=(C1*A(IB+I)+S1*B(IB+I))+(C3*A(ID+I)+S3*B(ID+I))
      A3=(C1*A(IB+I)+S1*B(IB+I))-(C3*A(ID+I)+S3*B(ID+I))
      B0=B(IA+I)+(C2*B(IC+I)-S2*A(IC+I))
      B2=B(IA+I)-(C2*B(IC+I)-S2*A(IC+I))
      B1=(C1*B(IB+I)-S1*A(IB+I))+(C3*B(ID+I)-S3*A(ID+I))
      B3=(C1*B(IB+I)-S1*A(IB+I))-(C3*B(ID+I)-S3*A(ID+I))
      C(JA+J)=A0+A1
      C(JC+J)=A0-A1
      D(JA+J)=B0+B1
      D(JC+J)=B1-B0
      C(JB+J)=A2+B3
      C(JD+J)=A2-B3
      D(JB+J)=B2-A3
      D(JD+J)=-(B2+A3)
      I=I+INC3
      J=J+INC4
  430 CONTINUE
      IBASE=IBASE+INC11
      JBASE=JBASE+INC2
  440 CONTINUE
      IBASE=IBASE+IJUMP
      JA=JA+JINK
      JB=JB+JINK
      JC=JC-JINK
      JD=JD-JINK
  450 CONTINUE
      ENDIF
!
      IF (JB.EQ.JC) THEN
      SIN45=SQRT(0.5_JPRB)
      JBASE=0
      DO 480 L=1,ILA
      I=IBASE
      J=JBASE
!OCL NOVREC
!DEC$ IVDEP
!NEC$ ivdep
      DO 470 IJK=1,ILOT
      C(JA+J)=A(IA+I)+SIN45*(A(IB+I)-A(ID+I))
      C(JB+J)=A(IA+I)-SIN45*(A(IB+I)-A(ID+I))
      D(JA+J)=-A(IC+I)-SIN45*(A(IB+I)+A(ID+I))
      D(JB+J)=A(IC+I)-SIN45*(A(IB+I)+A(ID+I))
      I=I+INC3
      J=J+INC4
  470 CONTINUE
      IBASE=IBASE+INC11
      JBASE=JBASE+INC2
  480 CONTINUE
      ENDIF
!
      ELSE              !!! Case LA=M
      Z=1.0_JPRB/REAL(N,KIND=JPRB)
      IF (LIPL) THEN
        DO 494 L=1,ILA
        I=IBASE
!OCL NOVREC
!DEC$ IVDEP
!NEC$ ivdep
        DO 492 IJK=1,ILOT
        T1=Z*(A(IA+I)-A(IC+I))
        T3=Z*(A(ID+I)-A(IB+I))
        T2=Z*((A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I)))
        A(IA+I)=Z*((A(IA+I)+A(IC+I))+(A(IB+I)+A(ID+I)))
        A(IB+I)=T1
        A(IC+I)=T2
        A(ID+I)=T3
        I=I+INC3
  492   CONTINUE
        IBASE=IBASE+INC11
  494   CONTINUE
      ELSE
        DO 498 L=1,ILA
        I=IBASE
        J=JBASE
!OCL NOVREC
!DEC$ IVDEP
!NEC$ ivdep
        DO 496 IJK=1,ILOT
        C(JA+J)=Z*((A(IA+I)+A(IC+I))+(A(IB+I)+A(ID+I)))
        C(JC+J)=Z*((A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I)))
        C(JB+J)=Z*(A(IA+I)-A(IC+I))
        D(JB+J)=Z*(A(ID+I)-A(IB+I))
        I=I+INC3
        J=J+INC4
  496   CONTINUE
        IBASE=IBASE+INC11
        JBASE=JBASE+INC2
  498   CONTINUE
      ENDIF
      ENDIF
!
      ELSEIF (IFAC.EQ.5) THEN
!
!     CODING FOR FACTOR 5
!     -------------------
  500 CONTINUE
      IA=1
      IB=IA+IINK
      IC=IB+IINK
      ID=IC+IINK
      IE=ID+IINK
      JA=1
      JB=JA+(2*M-LA)*INC2
      JC=JB+2*M*INC2
      JD=JC
      JE=JB
!
      IF (LA.NE.M) THEN
!
      DO 520 L=1,ILA
      I=IBASE
      J=JBASE
!OCL NOVREC
!DEC$ IVDEP
!NEC$ ivdep
      DO 510 IJK=1,ILOT
      A1=A(IB+I)+A(IE+I)
      A3=A(IB+I)-A(IE+I)
      A2=A(IC+I)+A(ID+I)
      A4=A(IC+I)-A(ID+I)
      A5=A(IA+I)-0.25_JPRB*(A1+A2)
      A6=QRT5*(A1-A2)
      C(JA+J)=A(IA+I)+(A1+A2)
      C(JB+J)=A5+A6
      C(JC+J)=A5-A6
      D(JB+J)=-SIN72*A3-SIN36*A4
      D(JC+J)=-SIN36*A3+SIN72*A4
      I=I+INC3
      J=J+INC4
  510 CONTINUE
      IBASE=IBASE+INC11
      JBASE=JBASE+INC2
  520 CONTINUE
      JA=JA+JINK
      JINK=2*JINK
      JB=JB+JINK
      JC=JC+JINK
      JD=JD-JINK
      JE=JE-JINK
      IBASE=IBASE+IJUMP
      IJUMP=2*IJUMP+IINK
!
      IF (JB.LT.JD) THEN
      DO 550 K=LA,KSTOP,LA
      KB=K+K
      KC=KB+KB
      KD=KC+KB
      KE=KD+KB
      C1=TRIGS(KB+1)
      S1=TRIGS(KB+2)
      C2=TRIGS(KC+1)
      S2=TRIGS(KC+2)
      C3=TRIGS(KD+1)
      S3=TRIGS(KD+2)
      C4=TRIGS(KE+1)
      S4=TRIGS(KE+2)
      JBASE=0
      DO 540 L=1,ILA
      I=IBASE
      J=JBASE
!OCL NOVREC
!DEC$ IVDEP
!NEC$ ivdep
      DO 530 IJK=1,ILOT
      A1=(C1*A(IB+I)+S1*B(IB+I))+(C4*A(IE+I)+S4*B(IE+I))
      A3=(C1*A(IB+I)+S1*B(IB+I))-(C4*A(IE+I)+S4*B(IE+I))
      A2=(C2*A(IC+I)+S2*B(IC+I))+(C3*A(ID+I)+S3*B(ID+I))
      A4=(C2*A(IC+I)+S2*B(IC+I))-(C3*A(ID+I)+S3*B(ID+I))
      B1=(C1*B(IB+I)-S1*A(IB+I))+(C4*B(IE+I)-S4*A(IE+I))
      B3=(C1*B(IB+I)-S1*A(IB+I))-(C4*B(IE+I)-S4*A(IE+I))
      B2=(C2*B(IC+I)-S2*A(IC+I))+(C3*B(ID+I)-S3*A(ID+I))
      B4=(C2*B(IC+I)-S2*A(IC+I))-(C3*B(ID+I)-S3*A(ID+I))
      A5=A(IA+I)-0.25_JPRB*(A1+A2)
      A6=QRT5*(A1-A2)
      B5=B(IA+I)-0.25_JPRB*(B1+B2)
      B6=QRT5*(B1-B2)
      A10=A5+A6
      A20=A5-A6
      B10=B5+B6
      B20=B5-B6
      A11=SIN72*B3+SIN36*B4
      A21=SIN36*B3-SIN72*B4
      B11=SIN72*A3+SIN36*A4
      B21=SIN36*A3-SIN72*A4
      C(JA+J)=A(IA+I)+(A1+A2)
      C(JB+J)=A10+A11
      C(JE+J)=A10-A11
      C(JC+J)=A20+A21
      C(JD+J)=A20-A21
      D(JA+J)=B(IA+I)+(B1+B2)
      D(JB+J)=B10-B11
      D(JE+J)=-(B10+B11)
      D(JC+J)=B20-B21
      D(JD+J)=-(B20+B21)
      I=I+INC3
      J=J+INC4
  530 CONTINUE
      IBASE=IBASE+INC11
      JBASE=JBASE+INC2
  540 CONTINUE
      IBASE=IBASE+IJUMP
      JA=JA+JINK
      JB=JB+JINK
      JC=JC+JINK
      JD=JD-JINK
      JE=JE-JINK
  550 CONTINUE
      ENDIF
!
      IF (JB.EQ.JD) THEN
      JBASE=0
      DO 580 L=1,ILA
      I=IBASE
      J=JBASE
!OCL NOVREC
!DEC$ IVDEP
!NEC$ ivdep
      DO 570 IJK=1,ILOT
      A1=A(IB+I)+A(IE+I)
      A3=A(IB+I)-A(IE+I)
      A2=A(IC+I)+A(ID+I)
      A4=A(IC+I)-A(ID+I)
      A5=A(IA+I)+0.25_JPRB*(A3-A4)
      A6=QRT5*(A3+A4)
      C(JA+J)=A5+A6
      C(JB+J)=A5-A6
      C(JC+J)=A(IA+I)-(A3-A4)
      D(JA+J)=-SIN36*A1-SIN72*A2
      D(JB+J)=-SIN72*A1+SIN36*A2
      I=I+INC3
      J=J+INC4
  570 CONTINUE
      IBASE=IBASE+INC11
      JBASE=JBASE+INC2
  580 CONTINUE
      ENDIF
!
      ELSE                !!! Case LA=M
      Z=1.0_JPRB/REAL(N,KIND=JPRB)
      ZQRT5=Z*QRT5
      ZSIN36=Z*SIN36
      ZSIN72=Z*SIN72
      IF (LIPL) THEN
        DO 594 L=1,ILA
        I=IBASE
!OCL NOVREC
!DEC$ IVDEP
!NEC$ ivdep
        DO 592 IJK=1,ILOT
        A1=A(IB+I)+A(IE+I)
        A3=A(IB+I)-A(IE+I)
        A2=A(IC+I)+A(ID+I)
        A4=A(IC+I)-A(ID+I)
        A5=Z*(A(IA+I)-0.25_JPRB*(A1+A2))
        A6=ZQRT5*(A1-A2)
        A(IA+I)=Z*(A(IA+I)+(A1+A2))
        A(IB+I)=A5+A6
        A(ID+I)=A5-A6
        A(IC+I)=-ZSIN72*A3-ZSIN36*A4
        A(IE+I)=-ZSIN36*A3+ZSIN72*A4
        I=I+INC3
  592   CONTINUE
        IBASE=IBASE+INC11
  594   CONTINUE
      ELSE
        DO 598 L=1,ILA
        I=IBASE
        J=JBASE
!OCL NOVREC
!DEC$ IVDEP
!NEC$ ivdep
        DO 596 IJK=1,ILOT
        A1=A(IB+I)+A(IE+I)
        A3=A(IB+I)-A(IE+I)
        A2=A(IC+I)+A(ID+I)
        A4=A(IC+I)-A(ID+I)
        A5=Z*(A(IA+I)-0.25_JPRB*(A1+A2))
        A6=ZQRT5*(A1-A2)
        C(JA+J)=Z*(A(IA+I)+(A1+A2))
        C(JB+J)=A5+A6
        C(JC+J)=A5-A6
        D(JB+J)=-ZSIN72*A3-ZSIN36*A4
        D(JC+J)=-ZSIN36*A3+ZSIN72*A4
        I=I+INC3
        J=J+INC4
  596   CONTINUE
        IBASE=IBASE+INC11
        JBASE=JBASE+INC2
  598   CONTINUE
      ENDIF
      ENDIF
!
      ELSEIF (IFAC.EQ.6) THEN
!
!     CODING FOR FACTOR 6
!     -------------------
  600 CONTINUE
      IA=1
      IB=IA+IINK
      IC=IB+IINK
      ID=IC+IINK
      IE=ID+IINK
      IF=IE+IINK
      JA=1
      JB=JA+(2*M-LA)*INC2
      JC=JB+2*M*INC2
      JD=JC+2*M*INC2
      JE=JC
      JF=JB
!
      IF (LA.NE.M) THEN
!
      DO 620 L=1,ILA
      I=IBASE
      J=JBASE
!OCL NOVREC
!DEC$ IVDEP
!NEC$ ivdep
      DO 610 IJK=1,ILOT
      A11=(A(IC+I)+A(IF+I))+(A(IB+I)+A(IE+I))
      C(JA+J)=(A(IA+I)+A(ID+I))+A11
      C(JC+J)=(A(IA+I)+A(ID+I)-0.5_JPRB*A11)
      D(JC+J)=SIN60*((A(IC+I)+A(IF+I))-(A(IB+I)+A(IE+I)))
      A11=(A(IC+I)-A(IF+I))+(A(IE+I)-A(IB+I))
      C(JB+J)=(A(IA+I)-A(ID+I))-0.5_JPRB*A11
      D(JB+J)=SIN60*((A(IE+I)-A(IB+I))-(A(IC+I)-A(IF+I)))
      C(JD+J)=(A(IA+I)-A(ID+I))+A11
      I=I+INC3
      J=J+INC4
  610 CONTINUE
      IBASE=IBASE+INC11
      JBASE=JBASE+INC2
  620 CONTINUE
      JA=JA+JINK
      JINK=2*JINK
      JB=JB+JINK
      JC=JC+JINK
      JD=JD-JINK
      JE=JE-JINK
      JF=JF-JINK
      IBASE=IBASE+IJUMP
      IJUMP=2*IJUMP+IINK
!
      IF (JC.LT.JD) THEN
      DO 650 K=LA,KSTOP,LA
      KB=K+K
      KC=KB+KB
      KD=KC+KB
      KE=KD+KB
      KF=KE+KB
      C1=TRIGS(KB+1)
      S1=TRIGS(KB+2)
      C2=TRIGS(KC+1)
      S2=TRIGS(KC+2)
      C3=TRIGS(KD+1)
      S3=TRIGS(KD+2)
      C4=TRIGS(KE+1)
      S4=TRIGS(KE+2)
      C5=TRIGS(KF+1)
      S5=TRIGS(KF+2)
      JBASE=0
      DO 640 L=1,ILA
      I=IBASE
      J=JBASE
!OCL NOVREC
!DEC$ IVDEP
!NEC$ ivdep
      DO 630 IJK=1,ILOT
      A1=C1*A(IB+I)+S1*B(IB+I)
      B1=C1*B(IB+I)-S1*A(IB+I)
      A2=C2*A(IC+I)+S2*B(IC+I)
      B2=C2*B(IC+I)-S2*A(IC+I)
      A3=C3*A(ID+I)+S3*B(ID+I)
      B3=C3*B(ID+I)-S3*A(ID+I)
      A4=C4*A(IE+I)+S4*B(IE+I)
      B4=C4*B(IE+I)-S4*A(IE+I)
      A5=C5*A(IF+I)+S5*B(IF+I)
      B5=C5*B(IF+I)-S5*A(IF+I)
      A11=(A2+A5)+(A1+A4)
      A20=(A(IA+I)+A3)-0.5_JPRB*A11
      A21=SIN60*((A2+A5)-(A1+A4))
      B11=(B2+B5)+(B1+B4)
      B20=(B(IA+I)+B3)-0.5_JPRB*B11
      B21=SIN60*((B2+B5)-(B1+B4))
      C(JA+J)=(A(IA+I)+A3)+A11
      D(JA+J)=(B(IA+I)+B3)+B11
      C(JC+J)=A20-B21
      D(JC+J)=A21+B20
      C(JE+J)=A20+B21
      D(JE+J)=A21-B20
      A11=(A2-A5)+(A4-A1)
      A20=(A(IA+I)-A3)-0.5_JPRB*A11
      A21=SIN60*((A4-A1)-(A2-A5))
      B11=(B5-B2)-(B4-B1)
      B20=(B3-B(IA+I))-0.5_JPRB*B11
      B21=SIN60*((B5-B2)+(B4-B1))
      C(JB+J)=A20-B21
      D(JB+J)=A21-B20
      C(JD+J)=A11+(A(IA+I)-A3)
      D(JD+J)=B11+(B3-B(IA+I))
      C(JF+J)=A20+B21
      D(JF+J)=A21+B20
      I=I+INC3
      J=J+INC4
  630 CONTINUE
      IBASE=IBASE+INC11
      JBASE=JBASE+INC2
  640 CONTINUE
      IBASE=IBASE+IJUMP
      JA=JA+JINK
      JB=JB+JINK
      JC=JC+JINK
      JD=JD-JINK
      JE=JE-JINK
      JF=JF-JINK
  650 CONTINUE
      ENDIF
!
      IF (JC.EQ.JD) THEN
      JBASE=0
      DO 680 L=1,ILA
      I=IBASE
      J=JBASE
!OCL NOVREC
!DEC$ IVDEP
!NEC$ ivdep
      DO 670 IJK=1,ILOT
      C(JA+J)=(A(IA+I)+0.5_JPRB*(A(IC+I)-A(IE+I)))+                     &
     &         SIN60*(A(IB+I)-A(IF+I))
      D(JA+J)=-(A(ID+I)+0.5_JPRB*(A(IB+I)+A(IF+I)))-                    &
     &         SIN60*(A(IC+I)+A(IE+I))
      C(JB+J)=A(IA+I)-(A(IC+I)-A(IE+I))
      D(JB+J)=A(ID+I)-(A(IB+I)+A(IF+I))
      C(JC+J)=(A(IA+I)+0.5_JPRB*(A(IC+I)-A(IE+I)))-                     &
     &         SIN60*(A(IB+I)-A(IF+I))
      D(JC+J)=-(A(ID+I)+0.5_JPRB*(A(IB+I)+                              &
     &         A(IF+I)))+SIN60*(A(IC+I)+A(IE+I))
      I=I+INC3
      J=J+INC4
  670 CONTINUE
      IBASE=IBASE+INC11
      JBASE=JBASE+INC2
  680 CONTINUE
      ENDIF
!
      ELSE                !!! Case LA=M
      Z=1.0_JPRB/REAL(N,KIND=JPRB)
      ZSIN60=Z*SIN60
      IF (LIPL) THEN
        DO 694 L=1,ILA
        I=IBASE
!OCL NOVREC
!DEC$ IVDEP
!NEC$ ivdep
        DO 692 IJK=1,ILOT
        A11=(A(IC+I)-A(IF+I))+(A(IE+I)-A(IB+I))
        T1=Z*((A(IA+I)-A(ID+I))-0.5_JPRB*A11)
        T5=Z*((A(IA+I)-A(ID+I))+A11)
        T2=ZSIN60*((A(IE+I)-A(IB+I))-(A(IC+I)-A(IF+I)))
        T4=ZSIN60*((A(IC+I)+A(IF+I))-(A(IB+I)+A(IE+I)))
        A11=(A(IC+I)+A(IF+I))+(A(IB+I)+A(IE+I))
        T3=Z*((A(IA+I)+A(ID+I))-0.5_JPRB*A11)
        A(IA+I)=Z*((A(IA+I)+A(ID+I))+A11)
        A(IB+I)=T1
        A(IC+I)=T2
        A(ID+I)=T3
        A(IE+I)=T4
        A(IF+I)=T5
        I=I+INC3
  692   CONTINUE
        IBASE=IBASE+INC11
  694   CONTINUE
      ELSE
        DO 698 L=1,ILA
        I=IBASE
        J=JBASE
!OCL NOVREC
!DEC$ IVDEP
!NEC$ ivdep
        DO 696 IJK=1,ILOT
        A11=(A(IC+I)+A(IF+I))+(A(IB+I)+A(IE+I))
        C(JA+J)=Z*((A(IA+I)+A(ID+I))+A11)
        C(JC+J)=Z*((A(IA+I)+A(ID+I))-0.5_JPRB*A11)
        D(JC+J)=ZSIN60*((A(IC+I)+A(IF+I))-(A(IB+I)+A(IE+I)))
        A11=(A(IC+I)-A(IF+I))+(A(IE+I)-A(IB+I))
        C(JB+J)=Z*((A(IA+I)-A(ID+I))-0.5_JPRB*A11)
        D(JB+J)=ZSIN60*((A(IE+I)-A(IB+I))-(A(IC+I)-A(IF+I)))
        C(JD+J)=Z*((A(IA+I)-A(ID+I))+A11)
        I=I+INC3
        J=J+INC4
  696   CONTINUE
        IBASE=IBASE+INC11
        JBASE=JBASE+INC2
  698   CONTINUE
      ENDIF
      ENDIF
!
      ELSEIF (IFAC.EQ.8) THEN
!
!     CODING FOR FACTOR 8
!     -------------------
  800 CONTINUE
      IF (LA.NE.M) THEN
        IBAD=3
      ELSE
      IA=1
      IB=IA+IINK
      IC=IB+IINK
      ID=IC+IINK
      IE=ID+IINK
      IF=IE+IINK
      IG=IF+IINK
      IH=IG+IINK
      JA=1
      JB=JA+LA*INC2
      JC=JB+2*M*INC2
      JD=JC+2*M*INC2
      JE=JD+2*M*INC2
      Z=1.0_JPRB/REAL(N,KIND=JPRB)
      ZSIN45=Z*SQRT(0.5_JPRB)
!
      IF (LIPL) THEN
        DO 820 L=1,ILA
        I=IBASE
!OCL NOVREC
!DEC$ IVDEP
!NEC$ ivdep
        DO 810 IJK=1,ILOT
        T3=Z*((A(IA+I)+A(IE+I))-(A(IC+I)+A(IG+I)))
        T4=Z*((A(ID+I)+A(IH+I))-(A(IB+I)+A(IF+I)))
        T1=Z*(A(IA+I)-A(IE+I))                                          &
     &    +ZSIN45*((A(IH+I)-A(ID+I))-(A(IF+I)-A(IB+I)))
        T5=Z*(A(IA+I)-A(IE+I))                                          &
     &    -ZSIN45*((A(IH+I)-A(ID+I))-(A(IF+I)-A(IB+I)))
        T2=ZSIN45*((A(IH+I)-A(ID+I))+(A(IF+I)-A(IB+I)))                 &
     &    +Z*(A(IG+I)-A(IC+I))
        T6=ZSIN45*((A(IH+I)-A(ID+I))+(A(IF+I)-A(IB+I)))                 &
     &    -Z*(A(IG+I)-A(IC+I))
        T7=Z*(((A(IA+I)+A(IE+I))+(A(IC+I)+A(IG+I)))-                    &
     &    ((A(ID+I)+A(IH+I))+(A(IB+I)+A(IF+I))))
        A(IA+I)=Z*(((A(IA+I)+A(IE+I))+(A(IC+I)+A(IG+I)))+               &
     &    ((A(ID+I)+A(IH+I))+(A(IB+I)+A(IF+I))))
        A(IB+I)=T1
        A(IC+I)=T2
        A(ID+I)=T3
        A(IE+I)=T4
        A(IF+I)=T5
        A(IG+I)=T6
        A(IH+I)=T7
        I=I+INC3
  810   CONTINUE
        IBASE=IBASE+INC11
  820   CONTINUE
      ELSE
        DO 840 L=1,ILA
        I=IBASE
        J=JBASE
!OCL NOVREC
!DEC$ IVDEP
        DO 830 IJK=1,ILOT
        C(JA+J)=Z*(((A(IA+I)+A(IE+I))+(A(IC+I)+A(IG+I)))+               &
     &    ((A(ID+I)+A(IH+I))+(A(IB+I)+A(IF+I))))
        C(JE+J)=Z*(((A(IA+I)+A(IE+I))+(A(IC+I)+A(IG+I)))-               &
     &    ((A(ID+I)+A(IH+I))+(A(IB+I)+A(IF+I))))
        C(JC+J)=Z*((A(IA+I)+A(IE+I))-(A(IC+I)+A(IG+I)))
        D(JC+J)=Z*((A(ID+I)+A(IH+I))-(A(IB+I)+A(IF+I)))
        C(JB+J)=Z*(A(IA+I)-A(IE+I))                                     &
     &    +ZSIN45*((A(IH+I)-A(ID+I))-(A(IF+I)-A(IB+I)))
        C(JD+J)=Z*(A(IA+I)-A(IE+I))                                     &
     &    -ZSIN45*((A(IH+I)-A(ID+I))-(A(IF+I)-A(IB+I)))
        D(JB+J)=ZSIN45*((A(IH+I)-A(ID+I))+(A(IF+I)-A(IB+I)))            &
     &    +Z*(A(IG+I)-A(IC+I))
        D(JD+J)=ZSIN45*((A(IH+I)-A(ID+I))+(A(IF+I)-A(IB+I)))            &
     &    -Z*(A(IG+I)-A(IC+I))
        I=I+INC3
        J=J+INC4
  830   CONTINUE
        IBASE=IBASE+INC11
        JBASE=JBASE+INC2
  840   CONTINUE
      ENDIF
!
      ENDIF
!
      ELSE
!
        IBAD=2        !!! Illegal factor
!
      ENDIF
!
!     RETURN
!     ------
  900 CONTINUE
      IERR=IBAD
      ENDSUBROUTINE QPASSF

      ENDSUBROUTINE FFT992
#endif
