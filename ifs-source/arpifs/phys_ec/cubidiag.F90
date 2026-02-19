! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE CUBIDIAG &
 & ( KIDIA, KFDIA, KLON, KLEV,&
 & KCTOP, LD_LCUMASK,&
 & PA,    PB,   PR,   PU )  

!          P. Bechtold         E.C.M.W.F.     07/03

!          PURPOSE.
!          --------
!          SOLVES BIDIAGONAL SYSTEM
!          FOR IMPLICIT SOLUTION OF ADVECTION EQUATION

!          INTERFACE
!          ---------

!          THIS ROUTINE IS CALLED FROM *CUDUDV* AND CUDTDQ.
!          IT RETURNS UPDATED VALUE OF QUANTITY

!          METHOD.
!          --------
!          NUMERICAL RECIPES (Cambridge Press)
!          DERIVED FROM TRIDIAGONAL ALGORIHM WITH C=0.
!          (ONLY ONE FORWARD SUBSTITUTION NECESSARY)
!          M  x  U  = R
!          ELEMENTS OF MATRIX M ARE STORED IN VECTORS A, B, C

!          (  B(kctop-1)    C(kctop-1)    0          0        )
!          (  A(kctop)      B(kctop)      C(kctop)   0        )
!          (  0             A(jk)         B(jk)      C(jk)    )
!          (  0             0             A(klev)    B(klev)  )

!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----

!    INPUT PARAMETERS (INTEGER):

!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KLEV*         NUMBER OF LEVELS
!    *KCTOP*        CLOUD TOP LEVELS

!    INPUT PARAMETERS (REAL):

!    *PA, PB*       VECTORS CONTAINING DIAGONAL ELEMENTS 
!    *PR*           RHS VECTOR CONTAINING "CONSTANTS"

!    OUTPUT PARAMETERS (REAL):

!    *PU*            SOLUTION VECTOR = UPDATED VALUE OF QUANTITY

!          EXTERNALS
!          ---------
!          NONE

!----------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

IMPLICIT NONE

!     DUMMY INTEGER 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KCTOP(KLON) 
LOGICAL           ,INTENT(IN)    :: LD_LCUMASK(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PA(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PB(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PR(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PU(KLON,KLEV) 
!     DUMMY REALS
!     DUMMY LOGICALS
!     LOCALS
INTEGER(KIND=JPIM) :: JK, JL
REAL(KIND=JPRB) :: ZBET
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CUBIDIAG',0,ZHOOK_HANDLE)

PU(:,:)=0.0_JPRB

! Forward Substitution

DO JK = 2, KLEV
  DO JL = KIDIA,KFDIA
    IF ( LD_LCUMASK(JL,JK) ) THEN
      IF ( JK==KCTOP(JL)-1 ) THEN
        ZBET      =1.0_JPRB/(PB(JL,JK)+1.E-35_JPRB)
        PU(JL,JK) = PR(JL,JK) * ZBET
      ELSEIF ( JK>KCTOP(JL)-1 ) THEN
        ZBET      = 1.0_JPRB/(PB(JL,JK) + 1.E-35_JPRB)
        PU(JL,JK) =(PR(JL,JK)-PA(JL,JK)*PU(JL,JK-1))*ZBET
      ENDIF
    ENDIF
  ENDDO
ENDDO

IF (LHOOK) CALL DR_HOOK('CUBIDIAG',1,ZHOOK_HANDLE)
END SUBROUTINE CUBIDIAG
