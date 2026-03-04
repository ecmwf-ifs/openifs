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

SUBROUTINE EXTFPNORM(YDFPGEO,PGP,KDOM,KDOMPTR,KSIZEG,KNFD,KFLDOFF,KNFG,KBFA,KIOMASTER,PAVE,PMIN,PMAX)

!**** *EXTFPNORM*  - COMPUTE GRIDPOINT NORMS - FULLPOS

!     PURPOSE.
!     --------
!        To compute  a set of gridpoint norms of the post-processed fields 
!        and send the results to the processor in charge of outputs.

!**   INTERFACE.
!     ----------
!       *CALL* *EXTFPNORM*

!        EXPLICIT ARGUMENTS
!        ------------------
!     PGP(:,:,:) - gridpoint fields (input)
!                  PGP is  dimensioned (NPROMA,IFIELDS,NGPBLKS) where
!                  NPROMA is the blocking factor, IFIELDS the total number
!                  of fields and NGPBLKS the number of NPROMA blocks.
!     KNFD       - Number of fields in the current chunk for each task
!     KFLDOFF    - Fields offset after distribution of chunk of fields among tasks.
!     KNFG       - Total number of fields in the current chunk
!     KBFA       - Absolute first field to be extracted
!     KIOMASTER  - MPI task in charge of receiving the norms
!     KDOM       - Number of subdomains for each field
!     KDOMPTR    - subdomains pointers for each field
!     KSIZEG     - size of subdomains
!     PAVE       - average (output) dim 1 = domains ; dim 2 = fields
!     PMIN       - minimum (output) dim 1 = domains ; dim 2 = fields
!     PMAX       - maximum (output) dim 1 = domains ; dim 2 = fields

!        IMPLICIT ARGUMENTS
!        ------------------
!        See modules above.

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     AUTHOR.
!     -------
!      RYAD EL KHATIB *METEO-FRANCE*
!      ORIGINAL : 18-Nov-2015

!     Modifications.
!     --------------
!     R. El Khatib 02-aug-2016 Mask (disabled for now)
!     --------------------------------------------------------------------------

USE PARKIND1    , ONLY : JPIM, JPRB
USE YOMHOOK     , ONLY : LHOOK, DR_HOOK, JPHOOK

USE YOMMP0      , ONLY : NPROC, MYPROC, NPRCIDS
USE YOMFPGEO, ONLY : TFPGEO
USE MPL_MODULE  , ONLY : MPL_SEND, MPL_RECV
USE IOSTREAM_MIX, ONLY : RMDI

!     --------------------------------------------------------------------------

IMPLICIT NONE

TYPE (TFPGEO),  INTENT(IN) :: YDFPGEO
REAL(KIND=JPRB),    INTENT(IN) :: PGP(:,:,:)
INTEGER(KIND=JPIM), INTENT(IN) :: KDOM(:)
INTEGER(KIND=JPIM), INTENT(IN) :: KDOMPTR(:,:)
INTEGER(KIND=JPIM), INTENT(IN) :: KSIZEG(:)
INTEGER(KIND=JPIM), INTENT(IN) :: KNFD(NPROC)
INTEGER(KIND=JPIM), INTENT(IN) :: KFLDOFF(NPROC)
INTEGER(KIND=JPIM), INTENT(IN) :: KNFG
INTEGER(KIND=JPIM), INTENT(IN) :: KBFA
INTEGER(KIND=JPIM), INTENT(IN) :: KIOMASTER
REAL(KIND=JPRB),INTENT(OUT)  :: PAVE(:,:)
REAL(KIND=JPRB),INTENT(OUT)  :: PMIN(:,:)
REAL(KIND=JPRB),INTENT(OUT)  :: PMAX(:,:)

!     --------------------------------------------------------------------------

! IFLDG   : global field address

INTEGER(KIND=JPIM) :: IDOM
INTEGER(KIND=JPIM) :: JFLD, JD, JROC
INTEGER(KIND=JPIM) :: IPTR(SIZE(KDOMPTR,DIM=1))
INTEGER(KIND=JPIM) :: IFLDG
REAL(KIND=JPRB) :: Z1SLEN(SIZE(KDOMPTR,DIM=1),KNFD(MYPROC))
REAL(KIND=JPRB) :: ZNORML(3,SIZE(KDOMPTR,DIM=1),SIZE(PGP,DIM=2))
REAL(KIND=JPRB) :: ZREALG(YDFPGEO%NFPRGPG,KNFD(MYPROC))
LOGICAL         :: LLMASK(YDFPGEO%NFPRGPG,KNFD(MYPROC))

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     --------------------------------------------------------------------------

#include "extfpf.intfb.h"

!     --------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('EXTFPNORM',0,ZHOOK_HANDLE)
ASSOCIATE(NFPRGPG=>YDFPGEO%NFPRGPG, NFPEND=>YDFPGEO%NFPEND, NFPRGPNUM=>YDFPGEO%NFPRGPNUM, &
 & NFPRGPIND=>YDFPGEO%NFPRGPIND)
!     --------------------------------------------------------------------------

IF (SIZE(KDOM) /= SIZE(PGP,DIM=2)) CALL ABOR1('FPNORMS: SIZES OF KDOM AND PGP(DIM=2) MISMATCH')
IF (SIZE(KDOMPTR,DIM=1) /= SIZE(KSIZEG)) CALL ABOR1('EXTFPNORM: SIZES OF KDOMPTR(DIM=1) AND SIZE(KSIZEG) MISMATCH')
IF (SIZE(KDOMPTR,DIM=2) /= SIZE(PGP,DIM=2)) CALL ABOR1('EXTFPNORM: SIZES OF KDOMPTR(DIM=2) AND PGP(DIM=2) MISMATCH')


!*       2.2 TRANSPOSITION FROM GRIDPOINT DISTRIBUTION TO FIELD DISTRIBUTION

CALL EXTFPF(KNFG,KNFD,KFLDOFF,KBFA,SIZE(PGP,DIM=1),SIZE(PGP,DIM=2),SIZE(PGP,DIM=3), &
 & NFPEND,NFPRGPG,NFPRGPNUM,NFPRGPIND,PGP,ZREALG)


!*       2.3 COMPUTE NORMS (DISTRIBUTED)

LLMASK(:,:)=.TRUE.

IPTR(1)=1
DO JD=2,SIZE(KSIZEG)
  IPTR(JD)=IPTR(JD-1)+KSIZEG(JD-1)
ENDDO
DO JFLD=1,KNFD(MYPROC)
  DO JD=1,SIZE(KSIZEG)
    LLMASK(IPTR(JD):IPTR(JD)+KSIZEG(JD)-1,JFLD) = LLMASK(IPTR(JD):IPTR(JD)+KSIZEG(JD)-1,JFLD) .AND. &
      & (ZREALG(IPTR(JD):IPTR(JD)+KSIZEG(JD)-1,JFLD) /= RMDI)
    Z1SLEN(JD,JFLD)=1.0_JPRB/REAL(COUNT(LLMASK(IPTR(JD):IPTR(JD)+KSIZEG(JD)-1,JFLD)),JPRB)
  ENDDO
ENDDO
IF (MYPROC == KIOMASTER) THEN
  DO JFLD=1,KNFD(MYPROC)
    IFLDG=KBFA-1+KFLDOFF(MYPROC)+JFLD
    DO JD=1,KDOM(IFLDG)
      IDOM=KDOMPTR(JD,IFLDG)
      PAVE(JD,KFLDOFF(MYPROC)+JFLD)=SUM   (ZREALG(IPTR(IDOM):IPTR(IDOM)+KSIZEG(IDOM)-1,JFLD), &
       & MASK=LLMASK(IPTR(IDOM):IPTR(IDOM)+KSIZEG(IDOM)-1,JFLD))*Z1SLEN(IDOM,JFLD)
      PMIN(JD,KFLDOFF(MYPROC)+JFLD)=MINVAL(ZREALG(IPTR(IDOM):IPTR(IDOM)+KSIZEG(IDOM)-1,JFLD), &
       & MASK=LLMASK(IPTR(IDOM):IPTR(IDOM)+KSIZEG(IDOM)-1,JFLD))
      PMAX(JD,KFLDOFF(MYPROC)+JFLD)=MAXVAL(ZREALG(IPTR(IDOM):IPTR(IDOM)+KSIZEG(IDOM)-1,JFLD), &
       & MASK=LLMASK(IPTR(IDOM):IPTR(IDOM)+KSIZEG(IDOM)-1,JFLD))
    ENDDO
  ENDDO
  DO JROC=1,NPROC
    IF (JROC /= KIOMASTER) THEN
      CALL MPL_RECV(ZNORML,KTAG=0,KSOURCE=NPRCIDS(JROC),CDSTRING='EXTFPNORM:')  
      DO JFLD=1,KNFD(JROC)
        IFLDG=KBFA-1+KFLDOFF(JROC)+JFLD
        DO JD=1,KDOM(IFLDG)
          PAVE(JD,KFLDOFF(JROC)+JFLD)=ZNORML(1,JD,JFLD)
          PMIN(JD,KFLDOFF(JROC)+JFLD)=ZNORML(2,JD,JFLD)
          PMAX(JD,KFLDOFF(JROC)+JFLD)=ZNORML(3,JD,JFLD)
        ENDDO
      ENDDO
    ENDIF
  ENDDO
ELSE
  DO JFLD=1,KNFD(MYPROC)
    IFLDG=KBFA-1+KFLDOFF(MYPROC)+JFLD
    DO JD=1,KDOM(IFLDG)
      IDOM=KDOMPTR(JD,IFLDG)
      ZNORML(1,JD,JFLD)=SUM   (ZREALG(IPTR(IDOM):IPTR(IDOM)+KSIZEG(IDOM)-1,JFLD), &
       & MASK=LLMASK(IPTR(IDOM):IPTR(IDOM)+KSIZEG(IDOM)-1,JFLD))*Z1SLEN(IDOM,JFLD)
      ZNORML(2,JD,JFLD)=MINVAL(ZREALG(IPTR(IDOM):IPTR(IDOM)+KSIZEG(IDOM)-1,JFLD), &
       & MASK=LLMASK(IPTR(IDOM):IPTR(IDOM)+KSIZEG(IDOM)-1,JFLD))
      ZNORML(3,JD,JFLD)=MAXVAL(ZREALG(IPTR(IDOM):IPTR(IDOM)+KSIZEG(IDOM)-1,JFLD), &
       & MASK=LLMASK(IPTR(IDOM):IPTR(IDOM)+KSIZEG(IDOM)-1,JFLD))
    ENDDO
  ENDDO
  CALL MPL_SEND(ZNORML,KTAG=0,KDEST=NPRCIDS(KIOMASTER),CDSTRING='EXTFPNORM:')
ENDIF

!     --------------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('EXTFPNORM',1,ZHOOK_HANDLE)
END SUBROUTINE EXTFPNORM
