! (C) Copyright 2000- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE GATH_GRID_32(PGPG,KPROMA,KFGATHG,KTO,KRESOL,PGP)

!**** *GATH_GRID_32* - Gather global gridpoint array from processors

!     Purpose.
!     --------
!        Interface routine for gathering gripoint array

!**   Interface.
!     ----------
!     CALL GATH_GRID_32(...)

!     Explicit arguments : 
!     -------------------- 
!     PGPG(:,:)   - Global gridpoint array
!     KFGATHG     - Global number of fields to be gathered
!     KPROMA      - blocking factor for gridpoint input
!     KTO(:)      - Processor responsible for gathering each field
!     KRESOL      - resolution tag  which is required ,default is the
!                   first defined resulution (input)
!     PGP(:,:,:)  - Local spectral array
!
!     Method.
!     -------

!     Externals.  SET_RESOL   - set resolution
!     ----------  GATH_GRID_32_CTL -  control routine

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-03-03

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB   ,JPRM

!ifndef INTERFACE

USE TPM_GEN
USE TPM_DIM
USE TPM_DISTR

USE SET_RESOL_MOD
USE GATH_GRID_32_CTL_MOD
USE ABORT_TRANS_MOD
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

!endif INTERFACE

IMPLICIT NONE

! Declaration of arguments

REAL(KIND=JPRM)    ,OPTIONAL, INTENT(OUT) :: PGPG(:,:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KPROMA
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KFGATHG
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KTO(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KRESOL
REAL(KIND=JPRM)             , INTENT(IN)  :: PGP(:,:,:)

!ifndef INTERFACE

INTEGER(KIND=JPIM) :: IFRECV,J,IUBOUND(3),IPROMA,IGPBLKS
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GATH_GRID_32',0,ZHOOK_HANDLE)
! Set current resolution
CALL SET_RESOL(KRESOL)

IPROMA = D%NGPTOT
IF(PRESENT(KPROMA)) THEN
  IPROMA = KPROMA
ENDIF
IGPBLKS = (D%NGPTOT-1)/IPROMA+1


IF(UBOUND(KTO,1) < KFGATHG) THEN
 CALL ABORT_TRANS('GATH_GRID_32: KTO TOO SHORT!')
ENDIF
 
IFRECV = 0
DO J=1,KFGATHG
  IF(KTO(J) < 1 .OR. KTO(J) > NPROC) THEN
    WRITE(NERR,*) 'GATH_GRID_32:ILLEGAL KTO VALUE',KTO(J),J
    CALL ABORT_TRANS('GATH_GRID_32:ILLEGAL KTO VALUE')
  ENDIF
  IF(KTO(J) == MYPROC) IFRECV = IFRECV+1
ENDDO

IUBOUND=UBOUND(PGP)
IF(IUBOUND(1) < IPROMA) THEN
  WRITE(NOUT,*)'GATH_GRID_32:FIRST DIM. OF PGP TOO SMALL ',IUBOUND(1),IPROMA
  CALL ABORT_TRANS('GATH_GRID_32:FIRST DIMENSION OF PGP TOO SMALL ')
ENDIF
IF(IUBOUND(2) < KFGATHG) THEN
  WRITE(NOUT,*)'GATH_GRID_32:SEC. DIM. OF PGP TOO SMALL ',IUBOUND(2),KFGATHG
  CALL ABORT_TRANS('GATH_GRID_32:SECOND DIMENSION OF PGP TOO SMALL ')
ENDIF
IF(IUBOUND(3) < IGPBLKS) THEN
  WRITE(NOUT,*)'GATH_GRID_32:THIRD DIM. OF PGP TOO SMALL ',IUBOUND(3),IGPBLKS
  CALL ABORT_TRANS('GATH_GRID_32:THIRD DIMENSION OF PGP TOO SMALL ')
ENDIF

IF(IFRECV > 0) THEN
  IF(.NOT.PRESENT(PGPG)) THEN
    CALL ABORT_TRANS('GATH_GRID_32:PGPG MISSING')
  ENDIF
  IF(UBOUND(PGPG,1) < D%NGPTOTG) THEN
    CALL ABORT_TRANS('GATH_GRID_32:FIRST DIMENSION OF PGPG TOO SMALL')
  ENDIF 
 IF(UBOUND(PGPG,2) < IFRECV) THEN
    CALL ABORT_TRANS('GATH_GRID_32:SECOND DIMENSION OF PGPG TOO SMALL')
  ENDIF
ENDIF

CALL GATH_GRID_32_CTL(PGPG,KFGATHG,IPROMA,KTO,PGP)

IF (LHOOK) CALL DR_HOOK('GATH_GRID_32',1,ZHOOK_HANDLE)
!endif INTERFACE

!     ------------------------------------------------------------------

END SUBROUTINE GATH_GRID_32

