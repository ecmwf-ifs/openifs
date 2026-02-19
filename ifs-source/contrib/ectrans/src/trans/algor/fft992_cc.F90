! (C) Copyright 1998- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE FFT992_CC (A, KINC, KJUMP, KN, KLOT, KSIGN)
! 
! Perform complex transforms with FFT992 like interface
!
! For KSIGN=-1 (Real -> Complex) call after of FFT992
! For KSIGN=1  (Complex -> Real) call before FFT992
!
USE PARKIND1, ONLY : JPIM, JPRB
IMPLICIT NONE
REAL(KIND=JPRB),INTENT(INOUT):: A(*)
INTEGER(KIND=JPIM),INTENT(IN):: KINC,KJUMP,KN,KLOT,KSIGN
REAL(KIND=JPRB),ALLOCATABLE :: ZWORK(:,:)
REAL(KIND=JPRB),ALLOCATABLE :: ZWORK1(:,:)
REAL(KIND=JPRB),ALLOCATABLE :: ZWORK2(:,:)
REAL(KIND=JPRB) :: ZN
INTEGER(KIND=JPIM) :: NH, NLOTH, I1, INCD, JLOT, I2, I1P, I2P, J, I1N, I2N

NH=KN/2
NLOTH=KLOT/2
I1=1
INCD=KINC*2
ZN=SQRT(REAL(KN,JPRB))


IF( KSIGN==-1)THEN

  IF( KJUMP /= 1 )THEN
    ALLOCATE(ZWORK(1:2*KLOT,0:2*KN))
    DO JLOT=1,NLOTH
      I2=I1+1
      I1P=I1+KINC
      I2P=I2+KINC
      DO J=0,NH
        ZWORK(1,J)=A(I1+INCD*J)-A(I2P+INCD*J)
        ZWORK(2,J)=A(I2+INCD*J)+A(I1P+INCD*J)
        ZWORK(1,KN-J)=A(I1+INCD*J)+A(I2P+INCD*J)
        ZWORK(2,KN-J)=A(I2+INCD*J)-A(I1P+INCD*J)
      ENDDO
      ! M normalization requires sqrt(M) ??
      DO J=0,KN-1
        A(I1+KINC*J)=ZWORK(1,J)*ZN
        A(I2+KINC*J)=ZWORK(2,J)*ZN
     ENDDO
      I1=I1+KJUMP*2
    ENDDO
    DEALLOCATE(ZWORK)
  ELSE
    ALLOCATE(ZWORK1(1:NLOTH,0:KN))
    ALLOCATE(ZWORK2(1:NLOTH,0:KN))
    DO J=0,NH
      DO JLOT=1,NLOTH
        ZWORK1(JLOT,   J)=A((JLOT-1)*2+1+INCD*J)-A((JLOT-1)*2+2+KINC+INCD*J)
        ZWORK2(JLOT,   J)=A((JLOT-1)*2+2+INCD*J)+A((JLOT-1)*2+1+KINC+INCD*J)
        ZWORK1(JLOT,KN-J)=A((JLOT-1)*2+1+INCD*J)+A((JLOT-1)*2+2+KINC+INCD*J)
        ZWORK2(JLOT,KN-J)=A((JLOT-1)*2+2+INCD*J)-A((JLOT-1)*2+1+KINC+INCD*J)
      ENDDO
    ENDDO
    DO J=0,KN-1
      DO JLOT=1,NLOTH
        A((JLOT-1)*2+1+J*KINC)=ZWORK1(JLOT,J)*ZN
      ENDDO
      DO JLOT=1,NLOTH
        A((JLOT-1)*2+2+J*KINC)=ZWORK2(JLOT,J)*ZN
      ENDDO
    ENDDO
    DEALLOCATE(ZWORK1,ZWORK2)
  ENDIF

ELSE

  IF( KJUMP /= 1 )THEN
    ALLOCATE(ZWORK(1:2*KLOT,0:2*KN))
    DO JLOT=1,NLOTH
      I2=I1+1
      I1N=I1+KN*KINC
      I2N=I2+KN*KINC
      DO J=1,NH-1
        ZWORK(1,2*J)=0.5D0*(A(I1+KINC*J)+A(I1N-KINC*J))
        ZWORK(2,2*J)=0.5D0*(A(I2+KINC*J)+A(I2N-KINC*J))
        ZWORK(1,2*J+1)=0.5D0*(A(I2+KINC*J)-A(I2N-KINC*J))
        ZWORK(2,2*J+1)=0.5D0*(A(I1N-KINC*J)-A(I1+KINC*J))
      ENDDO
      ZWORK(1,0)=A(I1)
      ZWORK(2,0)=A(I2)
      ZWORK(1,1)=0.0D0
      ZWORK(2,1)=0.0D0
      ZWORK(1,KN)=A(I1N)
      ZWORK(2,KN)=A(I2N)
      ZWORK(1,KN+1)=0.0D0
      ZWORK(2,KN+1)=0.0D0
      DO J=0,KN+1
        A(I1+KINC*J)=ZWORK(1,J)
        A(I2+KINC*J)=ZWORK(2,J)
      ENDDO
      I1=I1+KJUMP*2
    ENDDO
    DEALLOCATE(ZWORK)
  ELSE
    ALLOCATE(ZWORK1(1:NLOTH,0:KN+1))
    ALLOCATE(ZWORK2(1:NLOTH,0:KN+1))
    DO J=1,NH-1
      DO JLOT=1,NLOTH
        ZWORK1(JLOT,2*J  )=0.5D0*(A((JLOT-1)*2+1+KINC*J)+A((JLOT-1)*2+1+(KN-J)*KINC))
        ZWORK2(JLOT,2*J  )=0.5D0*(A((JLOT-1)*2+2+KINC*J)+A((JLOT-1)*2+2+(KN-J)*KINC))
        ZWORK1(JLOT,2*J+1)=0.5D0*(A((JLOT-1)*2+2+KINC*J)-A((JLOT-1)*2+2+(KN-J)*KINC))
        ZWORK2(JLOT,2*J+1)=0.5D0*(A((JLOT-1)*2+1+(KN-J)*KINC)-A((JLOT-1)*2+1+KINC*J))
        ZWORK1(JLOT,0    )=A((JLOT-1)*2+1)
        ZWORK2(JLOT,0    )=A((JLOT-1)*2+2)
        ZWORK1(JLOT,1    )=0.0D0
        ZWORK2(JLOT,1    )=0.0D0
        ZWORK1(JLOT,KN   )=A((JLOT-1)*2+1+KN*KINC)
        ZWORK2(JLOT,KN   )=A((JLOT-1)*2+2+KN*KINC)
        ZWORK1(JLOT,KN+1 )=0.0D0
        ZWORK2(JLOT,KN+1 )=0.0D0
      ENDDO
    ENDDO
    DO J=0,KN+1
      DO JLOT=1,NLOTH
        A((JLOT-1)*2+1+KINC*J)=ZWORK1(JLOT,J)
        A((JLOT-1)*2+2+KINC*J)=ZWORK2(JLOT,J)
      ENDDO
    ENDDO
    DEALLOCATE(ZWORK1,ZWORK2)
  ENDIF

ENDIF

RETURN
END SUBROUTINE FFT992_CC
