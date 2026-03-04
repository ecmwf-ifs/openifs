! (C) Copyright 2005- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

!OPTIONS NOOPT
MODULE STRHANDLER_MOD

USE EC_PARKIND  ,ONLY : JPIM, JPRM, JPRD

IMPLICIT NONE

PRIVATE

PUBLIC :: TOLOWER, TOUPPER, EXPAND_STRING
PUBLIC :: SADJUSTL, SADJUSTR
PUBLIC :: STRANSFER

INTERFACE STRANSFER
MODULE PROCEDURE &
  STRANSFER_R8_TO_STR, STRANSFER_STR_TO_R8, &
& STRANSFER_R4_TO_STR, STRANSFER_STR_TO_R4
END INTERFACE

CONTAINS

FUNCTION STRANSFER_R8_TO_STR(SOURCE, MOLD) RESULT(C)
REAL(KIND=JPRD) , INTENT(IN) :: SOURCE
CHARACTER(LEN=*), INTENT(IN) :: MOLD
CHARACTER(LEN=8) :: C
CALL ECMWF_TRANSFER(C,MIN(8,LEN(MOLD)),SOURCE,8)
END FUNCTION STRANSFER_R8_TO_STR


FUNCTION STRANSFER_STR_TO_R8(SOURCE, MOLD) RESULT(Z)
CHARACTER(LEN=*), INTENT(IN) :: SOURCE
REAL(KIND=JPRD) , INTENT(IN) :: MOLD
REAL(KIND=JPRD) :: Z
CALL ECMWF_TRANSFER(Z,8,SOURCE,LEN(SOURCE))
END FUNCTION STRANSFER_STR_TO_R8


FUNCTION STRANSFER_R4_TO_STR(SOURCE, MOLD) RESULT(C)
REAL(KIND=JPRM) , INTENT(IN) :: SOURCE
CHARACTER(LEN=*), INTENT(IN) :: MOLD
CHARACTER(LEN=4) :: C
CALL ECMWF_TRANSFER(C,MIN(4,LEN(MOLD)),SOURCE,4)
END FUNCTION STRANSFER_R4_TO_STR


FUNCTION STRANSFER_STR_TO_R4(SOURCE, MOLD) RESULT(Z)
CHARACTER(LEN=*), INTENT(IN) :: SOURCE
REAL(KIND=JPRM) , INTENT(IN) :: MOLD
REAL(KIND=JPRM) :: Z
CALL ECMWF_TRANSFER(Z,4,SOURCE,LEN(SOURCE))
END FUNCTION STRANSFER_STR_TO_R4


FUNCTION SADJUSTL(S) RESULT(C)
CHARACTER(LEN=*), INTENT(IN) :: S
CHARACTER(LEN=MAX(1,LEN(S))) C
C = ' '
IF (LEN(S) > 0) THEN
  IF (S /= ' ') C = ADJUSTL(S)
ENDIF
END FUNCTION SADJUSTL


FUNCTION SADJUSTR(S) RESULT(C)
CHARACTER(LEN=*), INTENT(IN) :: S
CHARACTER(LEN=MAX(1,LEN(S))) C
C = ' '
IF (LEN(S) > 0) THEN
  IF (S /= ' ') C = ADJUSTR(S)
ENDIF
END FUNCTION SADJUSTR


SUBROUTINE TOLOWER(CDS)
CHARACTER(LEN=*), INTENT(INOUT) :: CDS
INTEGER(KIND=JPIM), PARAMETER :: ICH_A = ICHAR('a')
INTEGER(KIND=JPIM), PARAMETER :: ICHA  = ICHAR('A')
INTEGER(KIND=JPIM), PARAMETER :: ICHZ  = ICHAR('Z')
INTEGER(KIND=JPIM) :: I, ICH, NEW_ICH
CHARACTER(LEN=1) CH
DO I=1,LEN(CDS)
  CH = CDS(I:I)
  ICH = ICHAR(CH)
  IF ( ICH >= ICHA .AND. ICH <= ICHZ ) THEN
    NEW_ICH = ICH + (ICH_A - ICHA)
    CH = CHAR(NEW_ICH)
    CDS(I:I) = CH
  ENDIF
ENDDO
END SUBROUTINE TOLOWER


SUBROUTINE TOUPPER(CDS)
CHARACTER(LEN=*), INTENT(INOUT) :: CDS
INTEGER(KIND=JPIM), PARAMETER :: ICH_A = ICHAR('A')
INTEGER(KIND=JPIM), PARAMETER :: ICHA  = ICHAR('a')
INTEGER(KIND=JPIM), PARAMETER :: ICHZ  = ICHAR('z')
INTEGER(KIND=JPIM) :: I, ICH, NEW_ICH
CHARACTER(LEN=1) CH
DO I=1,LEN(CDS)
  CH = CDS(I:I)
  ICH = ICHAR(CH)
  IF ( ICH >= ICHA .AND. ICH <= ICHZ ) THEN
    NEW_ICH = ICH + (ICH_A - ICHA)
    CH = CHAR(NEW_ICH)
    CDS(I:I) = CH
  ENDIF
ENDDO
END SUBROUTINE TOUPPER


SUBROUTINE EXPAND_STRING(&
     &MYPROC,               &! %p
     &nproc,                &! %n
     &timestep,             &! %t
     &max_timestep,&
     &s)                   ! %s

INTEGER(KIND=JPIM), INTENT(IN)          :: MYPROC, NPROC
INTEGER(KIND=JPIM), INTENT(IN)          :: TIMESTEP, MAX_TIMESTEP
CHARACTER(LEN=*), INTENT(INOUT) :: S(:)
CHARACTER(LEN=2*LEN(S))  T
CHARACTER(LEN=2*LEN(S)) TT
INTEGER(KIND=JPIM) :: I, J, JJ, LOC_P, LEN_T, N
INTEGER(KIND=JPIM) :: NDIGS(4), NUM(4)
CHARACTER(LEN=6) FMT(4)

N = SIZE(S)

IF (N < 1) RETURN

!*    Setup output formats
NUM(1) = MYPROC
NUM(2) = MAX(NPROC,MYPROC)
NUM(3) = N
NUM(4) = MAX(MAX_TIMESTEP,TIMESTEP)

!*    Count number of digits in each integer
DO J=1,4
  NDIGS(J) = 1
  IF (NUM(J) /= 0) THEN
    NDIGS(J) = 1 + LOG10(DBLE(ABS(NUM(J))))
    IF (NUM(J) < 0) NDIGS(J) = NDIGS(J) + 1 ! Room for minus sign
  ENDIF
  NDIGS(J) = MIN(NDIGS(J),9)   ! Max 9 digits supported; i.e. '999999999'
  WRITE(FMT(J),'("(i",i1,")")') NDIGS(J)
ENDDO


!*    Expand fields '%s', '%p', '%n' and '%t' with their values


!*    A special treatment with the sequence numbering
IF (N>1) THEN
  LOC_P = INDEX(S(1),'%s')
  IF (LOC_P > 0) THEN
    S(2:) = S(1)
  ENDIF
ENDIF

DO I=1,N
  T = ADJUSTL(S(I))//' '
  LOC_P = INDEX(T,'%')

  IF (LOC_P > 0) THEN
    LEN_T = LEN_TRIM(T)
    J = LOC_P
    TT(:J-1) = T(:J-1)
    TT(J:) = ' '
    JJ = J-1

    DO WHILE (J <= LEN_T)
      IF (T(J:J) == '%') THEN
        J = J + 1
        IF (J <= LEN_T) THEN
          SELECT CASE ( T(J:J) )
          CASE ( 'p' )   ! myproc
          WRITE(TT(JJ+1:JJ+NDIGS(1)),FMT(1)) MYPROC
          JJ = JJ + NDIGS(1)
          CASE ( 'n' )   ! nproc
          WRITE(TT(JJ+1:JJ+NDIGS(2)),FMT(2)) NPROC
          JJ = JJ + NDIGS(2)
          CASE ( 's' )   ! sequence number i=[1..n]
          WRITE(TT(JJ+1:JJ+NDIGS(3)),FMT(3)) I
          JJ = JJ + NDIGS(3)
          CASE ( 't' )   ! timestep
          WRITE(TT(JJ+1:JJ+NDIGS(4)),FMT(4)) TIMESTEP
          JJ = JJ + NDIGS(4)
          CASE DEFAULT
          TT(JJ+1:JJ+2) = '%'//T(J:J)
          JJ = JJ + 2
          END SELECT
        ELSE
          TT(JJ+1:JJ+1) = '%'
          JJ = JJ + 1
        ENDIF
      ELSE
        TT(JJ+1:JJ+1) = T(J:J)
        JJ = JJ + 1
      ENDIF
      J = J + 1
    ENDDO

    T = ADJUSTL(TT)

!*   Get also rid of any blanks in the middle of the string

    LEN_T = LEN_TRIM(T)
    J = 1
    DO WHILE (J < LEN_T)
      IF (T(J:J) == ' ') THEN
        T(J:) = T(J+1:)
        LEN_T = LEN_TRIM(T)
      ELSE
        J = J + 1
      ENDIF
    ENDDO

  ENDIF

  S(I) = T
ENDDO

END SUBROUTINE EXPAND_STRING

END MODULE STRHANDLER_MOD
