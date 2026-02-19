! (C) Copyright 2000- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE INIGPTR_MOD
CONTAINS
SUBROUTINE INIGPTR(KGPTRSEND,KGPTRRECV)

!     Compute tables to assist GP to/from Fourier space transpositions

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE TPM_GEN         ,ONLY : NOUT
USE TPM_DISTR       ,ONLY : D, NPRTRNS
USE TPM_TRANS       ,ONLY : NGPBLKS, NPROMA
USE EQ_REGIONS_MOD  ,ONLY : MY_REGION_EW
USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS
!

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(OUT) :: KGPTRSEND(2,NGPBLKS,NPRTRNS)
INTEGER(KIND=JPIM),INTENT(OUT) :: KGPTRRECV(NPRTRNS)

INTEGER(KIND=JPIM) :: IBLOCK,IROF,IBFIRST,IPROCLAST,IPROC,IFIRST,ILAST,IBLAST
INTEGER(KIND=JPIM) :: JGL,JBL,JPRTRNS,JBLKS
!     Compute tables to assist GP to/from Fourier space transpositions


KGPTRSEND(:,:,:)=0
IBLOCK=1
IROF=1
IBFIRST=1
IPROCLAST=D%NPROCL(D%NFRSTLOFF+1)
DO JGL=1,D%NDGL_GP
  !   Find processor which deals with this latitude in Fourier distribution
  IPROC=D%NPROCL(D%NFRSTLOFF+JGL)
  IF(IPROC > NPRTRNS) THEN
    WRITE(NOUT,'(A,I8)')&
     &' INIGPTR ERROR : exceeding processor limit ',NPRTRNS
    CALL ABORT_TRANS(' INIGPTR ERROR : exceeding processor limit ')
  ENDIF

  !           for each latitude on this processor, find first and last points
  !           for each NPROMA chunk, for each destination processor
  IF(IPROC /= IPROCLAST) THEN
    IF(IROF > 1) THEN
      KGPTRSEND(1,IBLOCK,IPROCLAST)=IBFIRST
      KGPTRSEND(2,IBLOCK,IPROCLAST)=IROF-1
    ENDIF
    IF(IROF <= NPROMA) IBFIRST=IROF
    IPROCLAST=IPROC
  ENDIF
  IFIRST=D%NSTA(D%NPTRFLOFF+JGL,MY_REGION_EW)
  ILAST =IFIRST + D%NONL(D%NPTRFLOFF+JGL,MY_REGION_EW) -1
  DO JBL=IFIRST,ILAST
    IF(IROF == NPROMA) THEN
      IBLAST=IROF
      KGPTRSEND(1,IBLOCK,IPROC)=IBFIRST
      KGPTRSEND(2,IBLOCK,IPROC)=IBLAST
      IF(IBLOCK < NGPBLKS) IBLOCK=IBLOCK+1
      IROF=0
      IBFIRST=1
    ENDIF
    IROF=IROF+1
  ENDDO
ENDDO
IF(IROF /= 1.AND.IROF /= IBFIRST) THEN
!           non-empty residual block after last latitude line
  IBLAST=IROF-1
  KGPTRSEND(1,IBLOCK,IPROC)=IBFIRST
  KGPTRSEND(2,IBLOCK,IPROC)=IBLAST
ENDIF
!         sum up over blocks
KGPTRRECV(:)=0
DO JPRTRNS=1,NPRTRNS
  DO JBLKS=1,NGPBLKS
    IF(KGPTRSEND(1,JBLKS,JPRTRNS) > 0) THEN
      KGPTRRECV(JPRTRNS)=KGPTRRECV(JPRTRNS)+&
       &KGPTRSEND(2,JBLKS,JPRTRNS)-KGPTRSEND(1,JBLKS,JPRTRNS)+1
    ENDIF
  ENDDO
ENDDO

END SUBROUTINE INIGPTR
END MODULE INIGPTR_MOD
