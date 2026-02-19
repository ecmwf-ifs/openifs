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

#if defined(RS6K)
      SUBROUTINE SGEMMX(NRP,NCP,M,                                      &
     &                  ONE,SA,IAC,IAR,SB,IBC,IBR,                      &
     &                  ZER,SC,ICC,ICR)

!AUTOPROMOTE
      USE PARKIND1, ONLY : JPRD, JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
      USE, INTRINSIC :: ieee_exceptions

      REAL(KIND=JPRB) :: SA(*),SB(*),SC(*)
      REAL(KIND=JPRB) :: A(NRP*M),B(M*NCP),C(NRP*NCP)
      INTEGER(KIND=8):: PA,PB,PC
      CHARACTER*1 CA1,CA2,CB1,CB2,CC1,CC2
      LOGICAL,PARAMETER :: LLDOUBLE = (JPRB == JPRD)
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('SGEMMX',0,ZHOOK_HANDLE)
!
!   Check if ONE and ZERO have valid values, else ABORT
!
      IF (ONE.NE.1.0D0.OR.ZER.NE.0.0D0)                                 &
     &     CALL ABOR1(' INVALID ARGUMENTS IN SGEMMX CALL')

      IF(NCP.EQ.0) THEN
        IF (LHOOK) CALL DR_HOOK('SGEMMX',1,ZHOOK_HANDLE)
        RETURN
      ENDIF

      PA=LOC(A)
      PB=LOC(B)
      PC=LOC(C)
      IF(IAC.EQ.1) THEN
        PA=LOC(SA)
        IAL=IAR
        CA1="N"
        CA2="T"
      ENDIF
      IF(IAR.EQ.1) THEN
        PA=LOC(SA)
        IAL=IAC
        CA1="N"
        CA2="T"
      ENDIF
      IF(IBC.EQ.1) THEN
        PB=LOC(SB)
        IBL=IBR
        CB1="N"
        CB2="T"
      ENDIF
      IF(IBR.EQ.1) THEN
        PB=LOC(SB)
        IBL=IBC
        CB1="N"
        CB2="T"
      ENDIF
      IF(ICC.EQ.1 .AND. ICR.EQ.NRP) PC=LOC(SC)
!     write(0,*) "IAC,IBC,ICC,IAR,IBR,ICR,NRP,NCP,M=",
!    &            IAC,IBC,ICC,IAR,IBR,ICR,NRP,NCP,M
!     write(0,*) "loc(a),pa",loc(a),pa
!     write(0,*) "loc(b),pb",loc(b),pb
!     write(0,*) "loc(c),pc",loc(c),pc

!     write(0,'(a,i12,9i8)')
!    &  "SGEMMX: NRP*NCP*M,NRP,NCP,M,IAC,IAR,IBC,IBR,ICC,ICR=",
!    &           NRP*NCP*M,NRP,NCP,M,IAC,IAR,IBC,IBR,ICC,ICR

      IF(PA.EQ.LOC(A)) THEN
!       write(0,*) "copying a"
        IF(IAR.GT.IAC) THEN
          DO J=1,M
            DO I=1,NRP
              A(I+(J-1)*NRP)=SA(1+(I-1)*IAC+(J-1)*IAR)
            ENDDO
          ENDDO
          IAL=NRP
          CA1="N"
          CA2="T"
        ELSE
          DO I=1,NRP
            DO J=1,M
              A(J+(I-1)*M)=SA(1+(I-1)*IAC+(J-1)*IAR)
            ENDDO
          ENDDO
          IAL=M
          CA1="N"
          CA2="T"
        ENDIF
      ENDIF
!
      IF(PB.EQ.LOC(B)) THEN
!       write(0,*) "copying b"
        IF(IBR.GT.IBC) THEN
          DO K=1,NCP
            DO J=1,M
              B(J+(K-1)*M)=SB(1+(J-1)*IBC+(K-1)*IBR)
            ENDDO
          ENDDO
          IBL=M
          CB1="N"
          CB2="T"
        ELSE
          DO J=1,M
            DO K=1,NCP
              B(K+(J-1)*NCP)=SB(1+(J-1)*IBC+(K-1)*IBR)
            ENDDO
          ENDDO
          IBL=NCP
          CB1="N"
          CB2="T"
        ENDIF
      ENDIF

!     CALL DGEMM('N','N',NRP,NCP,M,ONE,A,NRP,B,M,ZER,C,NRP)

      CALL SGEMMXZ(A,B)

      IF(ICR.GT.ICC) THEN
        IF(IBR.GT.IBC) THEN
          IF(IAR.GT.IAC) THEN
            IF (LLDOUBLE) THEN
              CALL DGEMUL(%VAL(PA),IAL,CA1,                             &
     &                  %VAL(PB),IBL,CB1,%VAL(PC),NRP,NRP,M,NCP)
            ELSE
              CALL SGEMUL(%VAL(PA),IAL,CA1,                             &
     &                  %VAL(PB),IBL,CB1,%VAL(PC),NRP,NRP,M,NCP)
            ENDIF
          ELSE
            IF (LLDOUBLE) THEN
              CALL DGEMUL(%VAL(PA),IAL,CA2,                             &
     &                  %VAL(PB),IBL,CB1,%VAL(PC),NRP,NRP,M,NCP)
            ELSE
              CALL SGEMUL(%VAL(PA),IAL,CA2,                             &
     &                  %VAL(PB),IBL,CB1,%VAL(PC),NRP,NRP,M,NCP)
            ENDIF
          ENDIF
        ELSE
          IF(IAR.GT.IAC) THEN
            IF (LLDOUBLE) THEN
              CALL DGEMUL(%VAL(PA),IAL,CA1,                             &
     &                  %VAL(PB),IBL,CB2,%VAL(PC),NRP,NRP,M,NCP)
            ELSE
              CALL SGEMUL(%VAL(PA),IAL,CA1,                             &
     &                  %VAL(PB),IBL,CB2,%VAL(PC),NRP,NRP,M,NCP)
            ENDIF
          ELSE
            IF (LLDOUBLE) THEN
              CALL DGEMUL(%VAL(PA),IAL,CA2,                             &
     &                  %VAL(PB),IBL,CB2,%VAL(PC),NRP,NRP,M,NCP)
            ELSE
              CALL SGEMUL(%VAL(PA),IAL,CA2,                             &
     &                  %VAL(PB),IBL,CB2,%VAL(PC),NRP,NRP,M,NCP)
            ENDIF
          ENDIF
        ENDIF
      ELSE
        IF(IBR.GT.IBC) THEN
          IF(IAR.GT.IAC) THEN
            IF (LLDOUBLE) THEN
              CALL DGEMUL(%VAL(PB),IBL,CB2,                             &
     &                  %VAL(PA),IAL,CA2,%VAL(PC),NCP,NCP,M,NRP)
            ELSE
              CALL SGEMUL(%VAL(PB),IBL,CB2,                             &
     &                  %VAL(PA),IAL,CA2,%VAL(PC),NCP,NCP,M,NRP)
            ENDIF
          ELSE
            IF (LLDOUBLE) THEN
              CALL DGEMUL(%VAL(PB),IBL,CB2,                             &
     &                  %VAL(PA),IAL,CA1,%VAL(PC),NCP,NCP,M,NRP)
            ELSE
              CALL SGEMUL(%VAL(PB),IBL,CB2,                             &
     &                  %VAL(PA),IAL,CA1,%VAL(PC),NCP,NCP,M,NRP)
            ENDIF
          ENDIF
        ELSE
          IF(IAR.GT.IAC) THEN
            IF (LLDOUBLE) THEN
              CALL DGEMUL(%VAL(PB),IBL,CB1,                             &
     &                  %VAL(PA),IAL,CA2,%VAL(PC),NCP,NCP,M,NRP)
            ELSE
              CALL SGEMUL(%VAL(PB),IBL,CB1,                             &
     &                  %VAL(PA),IAL,CA2,%VAL(PC),NCP,NCP,M,NRP)
            ENDIF
          ELSE
            IF (LLDOUBLE) THEN
              CALL DGEMUL(%VAL(PB),IBL,CB1,                             &
     &                  %VAL(PA),IAL,CA1,%VAL(PC),NCP,NCP,M,NRP)
            ELSE
              CALL SGEMUL(%VAL(PB),IBL,CB1,                             &
     &                  %VAL(PA),IAL,CA1,%VAL(PC),NCP,NCP,M,NRP)
            ENDIF
          ENDIF
        ENDIF
      ENDIF

      IF(PC.EQ.LOC(C)) THEN
!       write(0,*) "copying c"
        IF(ICR.GT.ICC) THEN
          DO K=1,NCP
            DO I=1,NRP
              SC(1+(I-1)*ICC+(K-1)*ICR)=C(I+NRP*(K-1))
            ENDDO
          ENDDO
        ELSE
          DO I=1,NRP
            DO K=1,NCP
              SC(1+(I-1)*ICC+(K-1)*ICR)=C(K+NCP*(I-1))
            ENDDO
          ENDDO
        ENDIF
      ENDIF

!     do i=1,nrp
!       do k=1,ncp
!         write(0,*) SC(1+(I-1)*ICC+(K-1)*ICR)
!       enddo
!     enddo

      IF (LHOOK) CALL DR_HOOK('SGEMMX',1,ZHOOK_HANDLE)

      RETURN
      END SUBROUTINE SGEMMX

      SUBROUTINE SGEMMXZ(A,B)
      END SUBROUTINE SGEMMXZ
#elif defined(NECSX) && defined(BLAS) && !defined(USE_MATMUL)
      SUBROUTINE SGEMMX(NRP,NCP,M,                                      &
     &                  ONE,SA,IAC,IAR,SB,IBC,IBR,                      &
     &                  ZER,SC,ICC,ICR)

!AUTOPROMOTE
      USE PARKIND1, ONLY : JPRD, JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
      use, intrinsic :: ieee_exceptions

!   C=A*B  , with arbitrary increments IAC,IAR, IBC,IBR, ICC,ICR
!
!  It will only work if ONE=1.0 and ZER=0.0
!
      LOGICAL :: LL_HALT_INVALID
#ifdef WITH_IEEE_HALT
      LOGICAL, PARAMETER :: LL_IEEE_HALT = .TRUE.
#else
      LOGICAL, PARAMETER :: LL_IEEE_HALT = .FALSE.
#endif
      LOGICAL,PARAMETER :: LLDOUBLE = (JPRB == JPRD)
      REAL(KIND=JPRB) :: SA(*),SB(*),SC(*)
      REAL(KIND=JPRB) :: A(NRP,M),B(M,NCP),C(NRP,NCP)
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('SGEMMX',0,ZHOOK_HANDLE)
!
!   Check if ONE and ZERO have valid values, else ABORT
!
      IF (ONE.NE.1.0D0.OR.ZER.NE.0.0D0)                                 &
     &     CALL ABOR1(' INVALID ARGUMENTS IN SGEMMX CALL')
!
      IF(NCP.EQ.0) THEN
        IF (LHOOK) CALL DR_HOOK('SGEMMX',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
!
      DO 1000 J=1,M
!cdir nodep
      DO 1000 I=1,NRP
        A(I,J)=SA(1+(I-1)*IAC+(J-1)*IAR)
1000  CONTINUE
!
      DO 1100 K=1,NCP
!cdir nodep
      DO 1100 J=1,M
        B(J,K)=SB(1+(J-1)*IBC+(K-1)*IBR)
1100  CONTINUE
!
      IF (LLDOUBLE) THEN
         CALL DGEMM('N','N',NRP,NCP,M,ONE,A,NRP,B,M,ZER,C,NRP)
      ELSE
         IF (LL_IEEE_HALT) THEN
            call ieee_get_halting_mode(ieee_invalid,LL_HALT_INVALID)
            if (LL_HALT_INVALID)                                           & 
     &         call ieee_set_halting_mode(ieee_invalid,.false.)
         ENDIF
         CALL SGEMM('N','N',NRP,NCP,M,ONE,A,NRP,B,M,ZER,C,NRP)
         if (LL_IEEE_HALT .and. LL_HALT_INVALID)                           &  
     &      call ieee_set_halting_mode(ieee_invalid,.true.)
      ENDIF
!
      DO 1200 K=1,NCP
!cdir nodep
      DO 1200 I=1,NRP
        SC(1+(I-1)*ICC+(K-1)*ICR)=C(I,K)
1200  CONTINUE
      IF (LHOOK) CALL DR_HOOK('SGEMMX',1,ZHOOK_HANDLE)
      RETURN
      END SUBROUTINE SGEMMX
#elif defined(BLAS) && !defined(USE_MATMUL)
!OCL  NOUNROLL,NOPREEX,NOEVAL
      SUBROUTINE SGEMMX(NRP,NCP,M,                                      &
     &                  ONE,SA,IAC,IAR,SB,IBC,IBR,                      &
     &                  ZER,SC,ICC,ICR)

!AUTOPROMOTE
      USE PARKIND1, ONLY : JPRD, JPIM, JPIB, JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
      use, intrinsic :: ieee_exceptions
      IMPLICIT NONE

      INTEGER(KIND=JPIM) :: NRP
      INTEGER(KIND=JPIM) :: NCP
      INTEGER(KIND=JPIM) :: M
      INTEGER(KIND=JPIM) :: IAC
      INTEGER(KIND=JPIM) :: IAL
      INTEGER(KIND=JPIM) :: IAR
      INTEGER(KIND=JPIM) :: IBC
      INTEGER(KIND=JPIM) :: IBR
      INTEGER(KIND=JPIM) :: ICC
      INTEGER(KIND=JPIM) :: ICR
      INTEGER(KIND=JPIM) :: I,J,K,KIAL,IBL
      REAL(KIND=JPRB) :: SA(*)
      REAL(KIND=JPRB) :: SB(*)
      REAL(KIND=JPRB) :: SC(*)
      REAL(KIND=JPRB) :: ONE
      REAL(KIND=JPRB) :: ZER
      REAL(KIND=JPRB) :: A(NRP*M),B(M*NCP),C(NRP*NCP)
      INTEGER(KIND=JPIB) :: PA
      INTEGER(KIND=JPIB) :: PB
      INTEGER(KIND=JPIB) :: PC
      CHARACTER(LEN=1) :: CA1,CA2,CB1,CB2,CC1,CC2
      LOGICAL, PARAMETER :: LLDOUBLE = (JPRB == JPRD)
      LOGICAL :: LL_HALT_INVALID
#ifdef WITH_IEEE_HALT
      LOGICAL, PARAMETER :: LL_IEEE_HALT = .TRUE.
#else
      LOGICAL, PARAMETER :: LL_IEEE_HALT = .FALSE.
#endif
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('SGEMMX',0,ZHOOK_HANDLE)
!
!   Check if ONE and ZERO have valid values, else ABORT
!
      IF (ONE.NE.1.0_JPRB.OR.ZER.NE.0.0_JPRB)                           &
     &     CALL ABOR1(' INVALID ARGUMENTS IN SGEMMX CALL')

      IF(NCP.EQ.0) THEN
        IF (LHOOK) CALL DR_HOOK('SGEMMX',1,ZHOOK_HANDLE)
        RETURN
      ENDIF

      PA=LOC(A)
      PB=LOC(B)
      PC=LOC(C)
      IF(IAC.EQ.1) THEN
        PA=LOC(SA)
        IAL=IAR
        CA1="N"
        CA2="T"
      ENDIF
      IF(IAR.EQ.1) THEN
        PA=LOC(SA)
        IAL=IAC
        CA1="N"
        CA2="T"
      ENDIF
      IF(IBC.EQ.1) THEN
        PB=LOC(SB)
        IBL=IBR
        CB1="N"
        CB2="T"
      ENDIF
      IF(IBR.EQ.1) THEN
        PB=LOC(SB)
        IBL=IBC
        CB1="N"
        CB2="T"
      ENDIF
      IF(ICC.EQ.1 .AND. ICR.EQ.NRP) PC=LOC(SC)

      IF(PA.EQ.LOC(A)) THEN
!       write(0,*) "copying a"
        IF(IAR.GT.IAC) THEN
!$OMP PARALLEL DO PRIVATE(i,j)
          DO J=1,M
!OCL  NOVREC
            DO I=1,NRP
              A(I+(J-1)*NRP)=SA(1+(I-1)*IAC+(J-1)*IAR)
            ENDDO
          ENDDO
!$OMP END PARALLEL DO
          IAL=NRP
          CA1="N"
          CA2="T"
        ELSE
!$OMP PARALLEL DO PRIVATE(i,j)
          DO I=1,NRP
!OCL  NOVREC
            DO J=1,M
              A(J+(I-1)*M)=SA(1+(I-1)*IAC+(J-1)*IAR)
            ENDDO
          ENDDO
!$OMP END PARALLEL DO
          IAL=M
          CA1="N"
          CA2="T"
        ENDIF
      ENDIF
!
      IF(PB.EQ.LOC(B)) THEN
!       write(0,*) "copying b"
        IF(IBR.GT.IBC) THEN
!$OMP PARALLEL DO PRIVATE(k,j)
          DO K=1,NCP
!OCL  NOVREC
            DO J=1,M
              B(J+(K-1)*M)=SB(1+(J-1)*IBC+(K-1)*IBR)
            ENDDO
          ENDDO
!$OMP END PARALLEL DO
          IBL=M
          CB1="N"
          CB2="T"
        ELSE
!$OMP PARALLEL DO PRIVATE(k,j)
          DO J=1,M
!OCL  NOVREC
            DO K=1,NCP
              B(K+(J-1)*NCP)=SB(1+(J-1)*IBC+(K-1)*IBR)
            ENDDO
          ENDDO
!$OMP END PARALLEL DO
          IBL=NCP
          CB1="N"
          CB2="T"
        ENDIF
      ENDIF

!     CALL DGEMM('N','N',NRP,NCP,M,ONE,A,NRP,B,M,ZER,C,NRP)

      CALL SGEMMXZ(A,B)

      IF(ICR.GT.ICC) THEN
        IF(IBR.GT.IBC) THEN
          IF(IAR.GT.IAC) THEN
            IF (LLDOUBLE) THEN
               CALL DGEMM(CA1,CB1,NRP,NCP,M,ONE,                           &
     &              %VAL(PA),IAL,%VAL(PB),IBL,ZER,%VAL(PC),NRP)
            ELSE
               IF (LL_IEEE_HALT) THEN
                  call ieee_get_halting_mode(ieee_invalid,LL_HALT_INVALID)
                  if (LL_HALT_INVALID)                                        & 
     &               call ieee_set_halting_mode(ieee_invalid,.false.)
               ENDIF
               CALL SGEMM(CA1,CB1,NRP,NCP,M,ONE,                           &
     &              %VAL(PA),IAL,%VAL(PB),IBL,ZER,%VAL(PC),NRP)
               if (LL_IEEE_HALT .and. LL_HALT_INVALID)                        &  
     &            call ieee_set_halting_mode(ieee_invalid,.true.)
            ENDIF
          ELSE
            IF (LLDOUBLE) THEN
               CALL DGEMM(CA2,CB1,NRP,NCP,M,ONE,                          &
     &             %VAL(PA),IAL,%VAL(PB),IBL,ZER,%VAL(PC),NRP)
            ELSE
               IF (LL_IEEE_HALT) THEN
                  call ieee_get_halting_mode(ieee_invalid,LL_HALT_INVALID)
                  if (LL_HALT_INVALID)                                        &  
     &               call ieee_set_halting_mode(ieee_invalid,.false.)
               ENDIF
               CALL SGEMM(CA2,CB1,NRP,NCP,M,ONE,                           &
     &              %VAL(PA),IAL,%VAL(PB),IBL,ZER,%VAL(PC),NRP)
               if (LL_IEEE_HALT .and. LL_HALT_INVALID)                        &  
     &            call ieee_set_halting_mode(ieee_invalid,.true.)
            ENDIF
          ENDIF
        ELSE
          IF(IAR.GT.IAC) THEN
            IF (LLDOUBLE) THEN
               CALL DGEMM(CA1,CB2,NRP,NCP,M,ONE,                           &
     &              %VAL(PA),IAL,%VAL(PB),IBL,ZER,%VAL(PC),NRP)
            ELSE
               IF (LL_IEEE_HALT) THEN
                  call ieee_get_halting_mode(ieee_invalid,LL_HALT_INVALID)
                  if (LL_HALT_INVALID)                                        & 
     &               call ieee_set_halting_mode(ieee_invalid,.false.)
               ENDIF
               CALL SGEMM(CA1,CB2,NRP,NCP,M,ONE,                           &
     &              %VAL(PA),IAL,%VAL(PB),IBL,ZER,%VAL(PC),NRP)
               if (LL_IEEE_HALT .and. LL_HALT_INVALID)                        & 
     &            call ieee_set_halting_mode(ieee_invalid,.true.)
            ENDIF
          ELSE
            IF (LLDOUBLE) THEN
               CALL DGEMM(CA2,CB2,NRP,NCP,M,ONE,                           &
     &              %VAL(PA),IAL,%VAL(PB),IBL,ZER,%VAL(PC),NRP)
            ELSE
               IF (LL_IEEE_HALT) THEN
                  call ieee_get_halting_mode(ieee_invalid,LL_HALT_INVALID)
                  if (LL_HALT_INVALID)                                        &  
     &               call ieee_set_halting_mode(ieee_invalid,.false.)
               ENDIF
               CALL SGEMM(CA2,CB2,NRP,NCP,M,ONE,                           &
     &              %VAL(PA),IAL,%VAL(PB),IBL,ZER,%VAL(PC),NRP)
               if (LL_IEEE_HALT .and. LL_HALT_INVALID)                        &  
     &            call ieee_set_halting_mode(ieee_invalid,.true.)
            ENDIF
          ENDIF
        ENDIF
      ELSE
        IF(IBR.GT.IBC) THEN
          IF(IAR.GT.IAC) THEN
            IF (LLDOUBLE) THEN
               CALL DGEMM(CB2,CA2,NCP,NRP,M,ONE,                           &
     &              %VAL(PB),IBL,%VAL(PA),IAL,ZER,%VAL(PC),NCP)
            ELSE
               IF (LL_IEEE_HALT) THEN
                  call ieee_get_halting_mode(ieee_invalid,LL_HALT_INVALID)
                  if (LL_HALT_INVALID)                                        &  
     &               call ieee_set_halting_mode(ieee_invalid,.false.)
               ENDIF
               CALL SGEMM(CB2,CA2,NCP,NRP,M,ONE,                           &
     &              %VAL(PB),IBL,%VAL(PA),IAL,ZER,%VAL(PC),NCP)
               if (LL_IEEE_HALT .and. LL_HALT_INVALID)                        &  
     &            call ieee_set_halting_mode(ieee_invalid,.true.)
            ENDIF
          ELSE
            IF (LLDOUBLE) THEN
               CALL DGEMM(CB2,CA1,NCP,NRP,M,ONE,                           &
     &              %VAL(PB),IBL,%VAL(PA),IAL,ZER,%VAL(PC),NCP)
            ELSE
               IF (LL_IEEE_HALT) THEN
                  call ieee_get_halting_mode(ieee_invalid,LL_HALT_INVALID)
                  if (LL_HALT_INVALID)                                        &  
     &               call ieee_set_halting_mode(ieee_invalid,.false.)
               ENDIF
               CALL SGEMM(CB2,CA1,NCP,NRP,M,ONE,                           &
     &              %VAL(PB),IBL,%VAL(PA),IAL,ZER,%VAL(PC),NCP)
               if (LL_IEEE_HALT .and. LL_HALT_INVALID)                        &  
     &            call ieee_set_halting_mode(ieee_invalid,.true.)
            ENDIF
          ENDIF
        ELSE
          IF(IAR.GT.IAC) THEN
            IF (LLDOUBLE) THEN
               CALL DGEMM(CB1,CA2,NCP,NRP,M,ONE,                           &
     &              %VAL(PB),IBL,%VAL(PA),IAL,ZER,%VAL(PC),NCP)
            ELSE
               IF (LL_IEEE_HALT) THEN
                  call ieee_get_halting_mode(ieee_invalid,LL_HALT_INVALID)
                  if (LL_HALT_INVALID)                                        &  
     &               call ieee_set_halting_mode(ieee_invalid,.false.)
               ENDIF
               CALL SGEMM(CB1,CA2,NCP,NRP,M,ONE,                           &
     &              %VAL(PB),IBL,%VAL(PA),IAL,ZER,%VAL(PC),NCP)
               if (LL_IEEE_HALT .and. LL_HALT_INVALID)                        &  
     &            call ieee_set_halting_mode(ieee_invalid,.true.)
            ENDIF
          ELSE
            IF (LLDOUBLE) THEN
               CALL DGEMM(CB1,CA1,NCP,NRP,M,ONE,                           &
     &              %VAL(PB),IBL,%VAL(PA),IAL,ZER,%VAL(PC),NCP)
            ELSE
               IF (LL_IEEE_HALT) THEN
                  call ieee_get_halting_mode(ieee_invalid,LL_HALT_INVALID)
                  if (LL_HALT_INVALID)                                        &  
     &               call ieee_set_halting_mode(ieee_invalid,.false.)
               ENDIF
               CALL SGEMM(CB1,CA1,NCP,NRP,M,ONE,                           &
     &              %VAL(PB),IBL,%VAL(PA),IAL,ZER,%VAL(PC),NCP)
               if (LL_IEEE_HALT .and. LL_HALT_INVALID)                        &  
     &            call ieee_set_halting_mode(ieee_invalid,.true.)
            ENDIF
          ENDIF
        ENDIF
      ENDIF

      IF(PC.EQ.LOC(C)) THEN
!       write(0,*) "copying c"
        IF(ICR.GT.ICC) THEN
!$OMP PARALLEL DO PRIVATE(k,i)
          DO K=1,NCP
!OCL  NOVREC
            DO I=1,NRP
              SC(1+(I-1)*ICC+(K-1)*ICR)=C(I+NRP*(K-1))
            ENDDO
          ENDDO
!$OMP END PARALLEL DO
        ELSE
!$OMP PARALLEL DO PRIVATE(k,i)
          DO I=1,NRP
!OCL  NOVREC
            DO K=1,NCP
              SC(1+(I-1)*ICC+(K-1)*ICR)=C(K+NCP*(I-1))
            ENDDO
          ENDDO
!$OMP END PARALLEL DO
        ENDIF
      ENDIF

      IF (LHOOK) CALL DR_HOOK('SGEMMX',1,ZHOOK_HANDLE)

      RETURN
      ENDSUBROUTINE SGEMMX

      SUBROUTINE SGEMMXZ(A,B)
      USE PARKIND1, ONLY : JPRD, JPRB
      IMPLICIT NONE
!
      REAL(KIND=JPRB) :: A(*)
      REAL(KIND=JPRB) :: B(*)
!RJ: missing code here?
      END SUBROUTINE SGEMMXZ
#elif defined(USE_MATMUL)
      SUBROUTINE SGEMMX(NRP,NCP,M,                                      &
     &                  ONE,SA,IAC,IAR,SB,IBC,IBR,                      &
     &                  ZER,SC,ICC,ICR)

!AUTOPROMOTE
      USE PARKIND1, ONLY : JPRD, JPIM, JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
!
!  Fortran subroutine equivalent to SGEMMX using MATMUL function
!  Written by F. Vana (august 2004)
!  Corrected on 2006-January-11
!
!
!   C=A*B  , with arbitrary increments IAC,IAR, IBC,IBR, ICC,ICR
!
!  It will only work if ONE=1.0 and ZER=0.0
!
!      IMPLICIT LOGICAL (L)
      IMPLICIT NONE
!
      INTEGER(KIND=JPIM) :: NRP
      INTEGER(KIND=JPIM) :: NCP
      INTEGER(KIND=JPIM) :: M
      INTEGER(KIND=JPIM) :: IAC
      INTEGER(KIND=JPIM) :: IAR
      INTEGER(KIND=JPIM) :: IBC
      INTEGER(KIND=JPIM) :: IBR
      INTEGER(KIND=JPIM) :: ICC
      INTEGER(KIND=JPIM) :: ICR
      REAL(KIND=JPRB) :: SA(*)
      REAL(KIND=JPRB) :: SB(*)
      REAL(KIND=JPRB) :: SC(*)
      REAL(KIND=JPRB) :: A(NRP,M),B(M,NCP),C(NRP,NCP)
      REAL(KIND=JPRB) :: ONE
      REAL(KIND=JPRB) :: ZER
      INTEGER(KIND=JPIM) :: I,J,K

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('SGEMMX',0,ZHOOK_HANDLE)
!
!   Check if ONE and ZERO have valid values, else ABORT
!
      IF(ONE.NE.1.0_JPRB.OR.ZER.NE.0.0_JPRB)                            &
     &     CALL ABOR1(' INVALID ARGUMENTS IN SGEMMX CALL')
!
      IF(NCP.EQ.0) THEN
        IF (LHOOK) CALL DR_HOOK('SGEMMX',1,ZHOOK_HANDLE)
        RETURN
      ENDIF

!
!  Preparation for matrix multiplication
!
!$OMP PARALLEL PRIVATE(i,j,k)
!OMP DO
!cdir OUTERUNROLL=4
      DO 1000 J=1,M
!cdir NODEP
      DO 1000 I=1,NRP
        A(I,J)=SA(1+(I-1)*IAC+(J-1)*IAR)
1000  CONTINUE
!OMP END DO
!
!OMP DO
!cdir OUTERUNROLL=4
      DO 1100 K=1,NCP
!cdir NODEP
      DO 1100 J=1,M
        B(J,K)=SB(1+(J-1)*IBC+(K-1)*IBR)
1100  CONTINUE
!OMP END DO
!$OMP END PARALLEL
!
!  Multiply matrices A and B
!
      C=MATMUL(A,B)
!
!  Final data transfer
!
!$OMP PARALLEL DO PRIVATE(i,k)
!cdir OUTERUNROLL=4
      DO 1200 K=1,NCP
!cdir NODEP
      DO 1200 I=1,NRP
        SC(1+(I-1)*ICC+(K-1)*ICR)=C(I,K)
1200  CONTINUE
!$OMP END PARALLEL DO

      IF (LHOOK) CALL DR_HOOK('SGEMMX',1,ZHOOK_HANDLE)

      RETURN
      END SUBROUTINE SGEMMX
#elif defined(VPP)
      SUBROUTINE SGEMMX(NRP,NCP,M,                                      &
     & ONE,SA,IAC,IAR,SB,IBC,IBR,                                       &
     & ZER,SC,ICC,ICR)

!AUTOPROMOTE
      USE PARKIND1, ONLY : JPRD, JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
!
!  Fortran subroutine equivalent to SGEMMX using SSLII dvmggm routine
!
!   C=A*B  , with arbitrary increments IAC,IAR, IBC,IBR, ICC,ICR
!
!  It will only work if ONE=1.0 and ZER=0.0
!
      REAL(KIND=JPRB) :: SA(*),SB(*),SC(*)
!     REAL(KIND=JPRB) :: A(NRP,m),B(m,NCP),C(NRP,NCP)
!
!   For maximum performance in dvmggm columns of matrices must
!   always start on a 64 byte boundary.
!
      REAL(KIND=JPRB),ALLOCATABLE :: A(:)

      INTEGER APADLEN,BPADLEN,CPADLEN
      INTEGER AOFFSET,BOFFSET,COFFSET
      INTEGER OFFSET

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('SGEMMX',0,ZHOOK_HANDLE)
!
!   Check if ONE and ZERO have valid values, else ABORT
!
      IF (ONE.NE.1.0D0.OR.ZER.NE.0.0D0) THEN
           WRITE(0,*) ' INVALID ARGUMENTS IN SGEMMX CALL'
           STOP
      ENDIF
!

      IF(NCP.EQ.0) THEN
        IF (LHOOK) CALL DR_HOOK('SGEMMX',1,ZHOOK_HANDLE)
        RETURN
      ENDIF

      NSIZE = 8 + (NRP+8)*M + (M+8)*NCP + (NRP+8)*NCP

      ALLOCATE(A(NSIZE))

      AOFFSET = MOD(LOC(A),64)/8
      APADLEN = ((NRP+8)/8)*8


      BOFFSET= AOFFSET + APADLEN*M
      BPADLEN = ((M+8)/8)*8

      COFFSET = BOFFSET + BPADLEN*NCP
      CPADLEN = APADLEN


!
      IF(ICC.EQ.1)THEN

      DO J=1,M
        OFFSET = AOFFSET + (J-1)*APADLEN
!OCL  NOVREC
        DO I=1,NRP
          A(OFFSET + I)=SA(1+(I-1)*IAC+(J-1)*IAR)
        ENDDO
      ENDDO
!
      DO K=1,NCP
        OFFSET = BOFFSET + (K-1)*BPADLEN
!OCL  NOVREC
        DO J=1,M
          A(OFFSET + J)=SB(1+(J-1)*IBC+(K-1)*IBR)
        ENDDO
      ENDDO

      CALL DVMGGM(A(AOFFSET+1),APADLEN,A(BOFFSET+1),BPADLEN,            &
     &            SC,ICR,NRP,M,NCP,ICON)
      IF(ICON.NE.0) STOP 'icon error at dvmggm SGEMMX.f'


      ELSE IF(IAC.EQ.1)THEN

!
      DO K=1,NCP
        OFFSET = BOFFSET + (K-1)*BPADLEN
!OCL  NOVREC
        DO J=1,M
          A(OFFSET + J)=SB(1+(J-1)*IBC+(K-1)*IBR)
        ENDDO
      ENDDO

      CALL DVMGGM(SA,IAR,A(BOFFSET+1),BPADLEN,                          &
     &            A(COFFSET+1),CPADLEN,NRP,M,NCP,ICON)
      IF(ICON.NE.0) STOP 'icon error at dvmggm SGEMMX.f'


      DO K=1,NCP
        OFFSET = COFFSET + (K-1)*CPADLEN
!OCL  NOVREC
        DO I=1,NRP
          SC(1+(I-1)*ICC+(K-1)*ICR)=A(OFFSET +I)
        ENDDO
      ENDDO

      ELSE

      DO J=1,M
        OFFSET = AOFFSET + (J-1)*APADLEN
!OCL  NOVREC
        DO I=1,NRP
          A(OFFSET + I)=SA(1+(I-1)*IAC+(J-1)*IAR)
        ENDDO
      ENDDO
!
      DO K=1,NCP
        OFFSET = BOFFSET + (K-1)*BPADLEN
!OCL  NOVREC
        DO J=1,M
          A(OFFSET + J)=SB(1+(J-1)*IBC+(K-1)*IBR)
        ENDDO
      ENDDO

      CALL DVMGGM(A(AOFFSET+1),APADLEN,A(BOFFSET+1),BPADLEN,            &
     &            A(COFFSET+1),CPADLEN,NRP,M,NCP,ICON)
      IF(ICON.NE.0) STOP 'icon error at dvmggm SGEMMX.f'


      DO K=1,NCP
        OFFSET = COFFSET + (K-1)*CPADLEN
!OCL  NOVREC
        DO I=1,NRP
          SC(1+(I-1)*ICC+(K-1)*ICR)=A(OFFSET +I)
        ENDDO
      ENDDO

      ENDIF

      DEALLOCATE(A)

      IF (LHOOK) CALL DR_HOOK('SGEMMX',1,ZHOOK_HANDLE)

      RETURN
      END SUBROUTINE SGEMMX
#else
      SUBROUTINE SGEMMX(NRP,NCP,M,                                      &
     &                  ONE,SA,IAC,IAR,SB,IBC,IBR,                      &
     &                  ZER,SC,ICC,ICR)

!AUTOPROMOTE
       USE PARKIND1, ONLY : JPRD, JPIM, JPRB
       USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
!
!  Fortran subroutine equivalent to SGEMMX, made from CRAY MXMA(3SCI)
!  See page 170 in CRAY Res. SR-2081 5.0
!
!   C=A*B  , with arbitrary increments IAC,IAR, IBC,IBR, ICC,ICR
!
!  It will only work if ONE=1.0 and ZER=0.0
!
!      IMPLICIT LOGICAL (L)
      IMPLICIT NONE
!
      INTEGER(KIND=JPIM) :: NRP
      INTEGER(KIND=JPIM) :: NCP
      INTEGER(KIND=JPIM) :: M
      INTEGER(KIND=JPIM) :: IAC
      INTEGER(KIND=JPIM) :: IAR
      INTEGER(KIND=JPIM) :: IBC
      INTEGER(KIND=JPIM) :: IBR
      INTEGER(KIND=JPIM) :: ICC
      INTEGER(KIND=JPIM) :: ICR
      REAL(KIND=JPRB) :: SA(*)
      REAL(KIND=JPRB) :: SB(*)
      REAL(KIND=JPRB) :: SC(*)
      REAL(KIND=JPRB) :: ONE
      REAL(KIND=JPRB) :: ZER
      INTEGER(KIND=JPIM) :: I,J,K
!
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('SGEMMX',0,ZHOOK_HANDLE)
!
!   Check if ONE and ZERO have valid values, else ABORT
!
      IF(ONE.NE.1.0_JPRB.OR.ZER.NE.0.0_JPRB)                            &
     &     CALL ABOR1(' INVALID ARGUMENTS IN SGEMMX CALL')
!  Initialize product
      DO 120 K=1,NCP
       DO 110 I=1,NRP
        SC(1+(I-1)*ICC+(K-1)*ICR) = 0.0_JPRB
!       (C(I,K) := 0.)
110    CONTINUE
120   CONTINUE
!  Multiply matrices from SA and SB
      DO 230 K=1,NCP
       DO 220 J=1,M
        DO 210 I=1,NRP
         SC(1+(I-1)*ICC+(K-1)*ICR)                                      &
     &     =SC(1+(I-1)*ICC+(K-1)*ICR)                                   &
     &      +SA(1+(I-1)*IAC+(J-1)*IAR)                                  &
     &       *SB(1+(J-1)*IBC+(K-1)*IBR)
!        (C(I,K):=C(I,K)+A(I,J)*B(J,K))
210     CONTINUE
220    CONTINUE
230   CONTINUE

      IF (LHOOK) CALL DR_HOOK('SGEMMX',1,ZHOOK_HANDLE)

      RETURN
      END SUBROUTINE SGEMMX
#endif
