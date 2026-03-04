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

SUBROUTINE TRIDIA(KN,KSYS,KFIRST,KEND,KTYP,PM,PRHS,PSOL)

!**** *TRIDIA*   SOLVES A NUMBER OF TRIDIAGONAL LINEAR SYSTEMS

!     If KTYP=1: solves PM*PSOL = PRHS(js) for js=kfirst,kend
!                (PM does not depend on "js", cf. former TRIDIAVSPL)
!     If KTYP=2: solves PM(js)*PSOL = PRHS(js) for js=kfirst,kend
!                (PM depends on "js", cf. former TRIDIALCZ)

!     Remark: this routine does something similar to routine SGTSL, but:
!      - it has been optimised to invert several tridiagonal linear systems
!        (SGTSL inverts only one tridiagonal linear system).
!      - it does not the additional tests and calculations done in SGTSL
!        (output quantity "kinfo").
!      - it is F90-norm compliant.

!        Explicit arguments :
!        --------------------
!            KN         : Dimension of the systems.                    (in)
!            KSYS       : Number of systems to be solved.              (in)
!            KFIRST     : First system to be solved.                   (in)
!            KEND       : Last system to be solved.                    (in)
!            KTYP       : Type of matrix inversion (see above)         (in)
!            PM         : Tridiagonal matrices of the systems          (inout)
!            PRHS       : RHS of the systems                           (in)
!            PSOL       : Solutions of the systems                     (out)

!     Author.
!     -------
!           OLIVIER TALAGRAND  -  16 SEPTEMBER 1988

!     Modifications.
!     --------------
!  MODIFICATION: ROBERTO BUIZZA - 29 JUL 92
!  K. Yessad (Oct 2008): merge of TRIDIALCZ and TRIDIAVSPL.
! -----------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

! -----------------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KN
INTEGER(KIND=JPIM),INTENT(IN)    :: KSYS
INTEGER(KIND=JPIM),INTENT(IN)    :: KFIRST
INTEGER(KIND=JPIM),INTENT(IN)    :: KEND
INTEGER(KIND=JPIM),INTENT(IN)    :: KTYP
REAL(KIND=JPRB)   ,INTENT(IN)    :: PM(1+(KTYP-1)*(KSYS-1),KN,-1:1)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRHS(KSYS,KN)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSOL(KSYS,KN)

! -----------------------------------------------------------------------------

REAL(KIND=JPRB) :: ZRHS(KSYS,KN)
REAL(KIND=JPRB) :: ZM(1+(KTYP-1)*(KSYS-1),KN,-1:1)
INTEGER(KIND=JPIM) :: J, JS

REAL(KIND=JPRB) :: ZDEN
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! -----------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('TRIDIA',0,ZHOOK_HANDLE)
! -----------------------------------------------------------------------------

!     1.   INITIALIZATIONS.
!          ----------------

!     1.1 FILL IN EXTREMITIES OF SECONDARY DIAGONALS

ZM=PM

IF (KTYP == 1) THEN
  ZM(1,1,-1) = 0.0_JPRB
  ZM(1,KN,1) = 0.0_JPRB
ELSEIF (KTYP == 2) THEN
  ZM(KFIRST:KEND,1,-1) = 0.0_JPRB
  ZM(KFIRST:KEND,KN,1) = 0.0_JPRB
ENDIF

!     1.2 Copy PRHS

ZRHS(KFIRST:KEND,1:KN)=PRHS(KFIRST:KEND,1:KN)

! -----------------------------------------------------------------------------

!     2.   ASCENDING LOOP.
!          ---------------

IF (KTYP == 1) THEN
  DO JS = KFIRST,KEND
    PSOL(JS,1) = -ZM(1,1,1)/ZM(1,1,0)
    ZRHS(JS,1) =  ZRHS(JS,1)/ZM(1,1,0)
  ENDDO
  DO J = 2,KN
    DO JS = KFIRST,KEND
      ZDEN = 1.0_JPRB/(ZM(1,J,-1)*PSOL(JS,J-1) + ZM(1,J,0))
      PSOL(JS,J) = -ZM(1,J,1)*ZDEN
      ZRHS(JS,J) = (ZRHS(JS,J) - ZRHS(JS,J-1)*ZM(1,J,-1))*ZDEN
    ENDDO
  ENDDO
ELSEIF (KTYP == 2) THEN
  DO JS = KFIRST,KEND
    PSOL(JS,1) = -ZM(JS,1,1)/ZM(JS,1,0)
    ZRHS(JS,1) =  ZRHS(JS,1)/ZM(JS,1,0)
  ENDDO
  DO J = 2,KN
    DO JS = KFIRST,KEND
      ZDEN = 1.0_JPRB/(ZM(JS,J,-1)*PSOL(JS,J-1) + ZM(JS,J,0))
      PSOL(JS,J) = -ZM(JS,J,1)*ZDEN
      ZRHS(JS,J) = (ZRHS(JS,J) - ZRHS(JS,J-1)*ZM(JS,J,-1))*ZDEN
    ENDDO
  ENDDO
ENDIF

! -----------------------------------------------------------------------------

!     3.   DESCENDING LOOP.
!          ----------------

PSOL(KFIRST:KEND,KN)=ZRHS(KFIRST:KEND,KN)
DO J = KN-1,1,-1
  DO JS = KFIRST,KEND
    PSOL(JS,J) = ZRHS(JS,J) + PSOL(JS,J)*PSOL(JS,J+1)
  ENDDO
ENDDO

! -----------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('TRIDIA',1,ZHOOK_HANDLE)
END SUBROUTINE TRIDIA
