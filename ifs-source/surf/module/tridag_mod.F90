MODULE TRIDAG_MOD

CONTAINS
SUBROUTINE TRIDAG(A,B,C,R,U,N,KIDIA,KFDIA,KLON,LLINCR)

USE PARKIND1 ,ONLY : JPIM     ,JPRB

! (C) Copyright 2005- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN) :: N, KIDIA, KFDIA, KLON
REAL(KIND=JPRB), DIMENSION(KLON,N), INTENT(IN) :: A, B, C, R
REAL(KIND=JPRB), DIMENSION(KLON,N), INTENT(INOUT) :: U
LOGICAL, INTENT(IN) :: LLINCR

INTEGER(KIND=JPIM) :: J, JL

REAL(KIND=JPRB), DIMENSION(KLON) :: BET
REAL(KIND=JPRB), DIMENSION(KLON,N) :: GAM, U0

!----------------------------------------------------------------------
 
IF (LLINCR) THEN
   DO J=1,N
     DO JL=KIDIA,KFDIA
       U0(JL,J) = 0.0_JPRB
     ENDDO
   ENDDO
ELSE
   DO J=1,N
     DO JL=KIDIA,KFDIA
       U0(JL,J) = U(JL,J)
     ENDDO
   ENDDO
ENDIF

DO JL=KIDIA,KFDIA
  BET(JL) = B(JL,1)
  U(JL,1) = R(JL,1)/BET(JL)
ENDDO
      
DO J=2,N
  DO JL=KIDIA,KFDIA
    GAM(JL,J) = C(JL,J-1)/BET(JL)
    BET(JL) = B(JL,J)-A(JL,J)*GAM(JL,J)
    U(JL,J) = (R(JL,J)-A(JL,J)*U(JL,J-1))/BET(JL)
  ENDDO
ENDDO

DO J=N-1,1,-1
  DO JL=KIDIA,KFDIA
    U(JL,J) = U(JL,J)-GAM(JL,J+1)*U(JL,J+1)
  ENDDO
ENDDO

DO J=1,N
  DO JL=KIDIA,KFDIA
    U(JL,J) = U0(JL,J) + U(JL,J)
  ENDDO
ENDDO


END SUBROUTINE TRIDAG

END MODULE TRIDAG_MOD
