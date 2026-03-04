! (C) Copyright 2000- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE DIST_GRID(PGPG,KPROMA,KFDISTG,KFROM,KRESOL,PGP,KSORT)

!**** *DIST_GRID* - Distribute global gridpoint array among processors

!     Purpose.
!     --------
!        Interface routine for distributing gridpoint array

!**   Interface.
!     ----------
!     CALL DIST_GRID(...)

!     Explicit arguments :
!     --------------------
!     PGPG(:,:) - Global spectral array
!     KFDISTG     - Global number of fields to be distributed
!     KPROMA      - required blocking factor for gridpoint input
!     KFROM(:)    - Processor resposible for distributing each field
!     KRESOL      - resolution tag  which is required ,default is the
!                   first defined resulution (input)
!     PGP(:,:)  - Local spectral array
!     KSORT (:) - Re-order fields on output
!
!     Method.
!     -------

!     Externals.  SET_RESOL      - set resolution
!     ----------  DIST_GRID_CTL  - control routine

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-03-03
!    P.Marguinaud : 10-10-14 Add KSORT

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

!ifndef INTERFACE

USE TPM_GEN         ,ONLY : NERR, NOUT
!USE TPM_DIM
USE TPM_DISTR       ,ONLY : D, MYPROC, NPROC

USE SET_RESOL_MOD   ,ONLY : SET_RESOL
USE DIST_GRID_CTL_MOD ,ONLY : DIST_GRID_CTL
USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

!endif INTERFACE

IMPLICIT NONE

! Declaration of arguments

REAL(KIND=JPRB)    ,OPTIONAL, INTENT(IN)  :: PGPG(:,:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KPROMA
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KFDISTG
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KFROM(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KRESOL
REAL(KIND=JPRB)             , INTENT(OUT) :: PGP(:,:,:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KSORT (:)

!ifndef INTERFACE

INTEGER(KIND=JPIM) :: IFSEND,J,IUBOUND(3),IPROMA,IGPBLKS
REAL(KIND=JPHOOK)  :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('DIST_GRID',0,ZHOOK_HANDLE)
! Set current resolution
CALL SET_RESOL(KRESOL)

IPROMA = D%NGPTOT
IF(PRESENT(KPROMA)) THEN
  IPROMA = KPROMA
ENDIF
IGPBLKS = (D%NGPTOT-1)/IPROMA+1

IF(UBOUND(KFROM,1) < KFDISTG) THEN
 CALL ABORT_TRANS('DIST_GRID: KFROM TOO SHORT!')
ENDIF
IFSEND = 0
DO J=1,KFDISTG
  IF(KFROM(J) < 1 .OR. KFROM(J) > NPROC) THEN
    WRITE(NERR,*) 'DIST_GRID:ILLEGAL KFROM VALUE',KFROM(J),J
    CALL ABORT_TRANS('DIST_GRID:ILLEGAL KFROM VALUE')
  ENDIF
  IF(KFROM(J) == MYPROC) IFSEND = IFSEND+1
ENDDO

IUBOUND=UBOUND(PGP)
IF(IUBOUND(1) < IPROMA) THEN
  WRITE(NOUT,*)'DIST_GRID:FIRST DIM. OF PGP TOO SMALL ',IUBOUND(1),IPROMA
  CALL ABORT_TRANS('DIST_GRID:FIRST DIMENSION OF PGP TOO SMALL ')
ENDIF
IF(IUBOUND(2) < KFDISTG) THEN
  WRITE(NOUT,*)'DIST_GRID:SEC. DIM. OF PGP TOO SMALL ',IUBOUND(2),KFDISTG
  CALL ABORT_TRANS('DIST_GRID:SECOND DIMENSION OF PGP TOO SMALL ')
ENDIF
IF(IUBOUND(3) < IGPBLKS) THEN
  WRITE(NOUT,*)'DIST_GRID:THIRD DIM. OF PGP TOO SMALL ',IUBOUND(3),IGPBLKS
  CALL ABORT_TRANS('DIST_GRID:THIRD DIMENSION OF PGP TOO SMALL ')
ENDIF

IF(IFSEND > 0) THEN
  IF(.NOT.PRESENT(PGPG)) THEN
    CALL ABORT_TRANS('DIST_GRID:PGPG MISSING')
  ENDIF
  IF(UBOUND(PGPG,1) < D%NGPTOTG) THEN
    CALL ABORT_TRANS('DIST_GRID:FIRST DIMENSION OF PGPG TOO SMALL')
  ENDIF
 IF(UBOUND(PGPG,2) < IFSEND) THEN
    CALL ABORT_TRANS('DIST_GRID:FIRST DIMENSION OF PGPG TOO SMALL')
  ENDIF
ENDIF

IF (PRESENT (KSORT)) THEN
  IF (UBOUND (KSORT, 1) /= UBOUND (PGP, 2)) THEN
    CALL ABORT_TRANS('DIST_GRID: DIMENSION MISMATCH KSORT, PGP')
  ENDIF
ENDIF

CALL DIST_GRID_CTL(PGPG,KFDISTG,IPROMA,KFROM,PGP,KSORT)

IF (LHOOK) CALL DR_HOOK('DIST_GRID',1,ZHOOK_HANDLE)
!endif INTERFACE

!     ------------------------------------------------------------------

END SUBROUTINE DIST_GRID

