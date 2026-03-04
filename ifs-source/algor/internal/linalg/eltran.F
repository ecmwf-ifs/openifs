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

      SUBROUTINE ELTRAN(NM,N,LOW,IGH,PA,KINT,PZ)
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
      INTEGER(KIND=JPIM) :: I,J,KL,MM,MP,MP1
      REAL(KIND=JPRB) :: PA(NM,IGH)
      REAL(KIND=JPRB) :: PZ(NM,N)
      INTEGER(KIND=JPIM) :: KINT(IGH)
!
!     this subroutine is a translation of the algol procedure elmtrans,
!     num. math. 16, 181-204(1970) by peters and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).
!
!     this subroutine accumulates the stabilized elementary
!     similarity transformations used in the reduction of a
!     real general matrix to upper hessenberg form by  elmhes.
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
!        pa contains the multipliers which were used in the
!          reduction by  elmhes  in its lower triangle
!          below the subdiagonal.
!
!        kint contains information on the rows and columns
!          interchanged in the reduction by  elmhes.
!          only elements low through igh are used.
!
!     on output
!
!        pz contains the transformation matrix produced in the
!          reduction by  elmhes.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
!     .......... initialize pz to identity matrix ..........
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('ELTRAN',0,ZHOOK_HANDLE)
      DO 80 J = 1, N
!
         DO 60 I = 1, N
   60    PZ(I,J) = 0.0_JPRB
!
         PZ(J,J) = 1.0_JPRB
   80 CONTINUE
!
      KL = IGH - LOW - 1
      IF (KL .LT. 1) GO TO 200
!     .......... for mp=igh-1 step -1 until low+1 do -- ..........
      DO 140 MM = 1, KL
         MP = IGH - MM
         MP1 = MP + 1
!
         DO 100 I = MP1, IGH
  100    PZ(I,MP) = PA(I,MP-1)
!
         I = KINT(MP)
         IF (I .EQ. MP) GO TO 140
!
         DO 130 J = MP, IGH
            PZ(MP,J) = PZ(I,J)
            PZ(I,J) = 0.0_JPRB
  130    CONTINUE
!
         PZ(I,MP) = 1.0_JPRB
  140 CONTINUE
!
  200 CONTINUE
      IF (LHOOK) CALL DR_HOOK('ELTRAN',1,ZHOOK_HANDLE)
      RETURN
      ENDSUBROUTINE ELTRAN
