! (C) Copyright 1995- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE TRMTOL_MOD

CONTAINS
SUBROUTINE TRMTOL(PFBUF_IN,PFBUF,KFIELD)

!**** *trmtol * - transposition in Fourier space

!     Purpose.
!     --------
!              Transpose Fourier buffer data from partitioning
!              over wave numbers to partitioning over latitudes.
!              It is called between direct FFT and direct Legendre
!              transform.
!              This routine is the inverse of TRLTOM.


!**   Interface.
!     ----------
!        *call* *trmtol(...)*

!        Explicit arguments : PFBUF  - Fourier coefficient buffer. It is
!        --------------------          used for both input and output.
!                             KFIELD - Number of fields communicated

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        MPP Group *ECMWF*

!     Modifications.
!     --------------
!        Original : 95-10-01
!        Modified : 97-06-17 G. Mozdzynski - control MPI mailbox use
!                                            (NCOMBFLEN) for nphase.eq.1
!        Modified : 99-05-28  D.Salmond - Optimise copies.
!        Modified : 00-02-02  M.Hamrud  - Remove NPHASE
!        D.Salmond : 01-11-23 LIMP_NOOLAP Option for non-overlapping message
!                             passing and buffer packing
!        G.Mozdzynski: 08-01-01 Cleanup
!        Y.Seity   : 07-08-31 add barrien synchronisation under LSYNC_TRANS
!     ------------------------------------------------------------------


USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE MPL_MODULE  ,ONLY : MPL_ALLTOALLV, MPL_BARRIER, MPL_ALL_MS_COMM, MPL_WAIT, JP_NON_BLOCKING_STANDARD

USE TPM_DISTR       ,ONLY : D, MTAGML, MYSETW, NPRTRW, NPROC
USE TPM_GEN         ,ONLY : LSYNC_TRANS

!USE SET2PE_MOD
!USE MYSENDSET_MOD
!USE MYRECVSET_MOD
!USE ABORT_TRANS_MOD
!

IMPLICIT NONE


INTEGER(KIND=JPIM),INTENT(IN)    :: KFIELD
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFBUF(:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFBUF_IN(:)

INTEGER(KIND=JPIM) :: ILENS(NPRTRW),IOFFS(NPRTRW),ILENR(NPRTRW),IOFFR(NPRTRW)

INTEGER(KIND=JPIM) :: ITAG, J, ILEN, ISTA
 
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE_BAR
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE_BAR2

REAL(KIND=JPRB)    :: ZDUM(1)
INTEGER(KIND=JPIM) :: IREQ


!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('TRMTOL',0,ZHOOK_HANDLE)


ITAG = MTAGML

DO J=1,NPRTRW
  ILENS(J) = D%NLTSFTB(J)*KFIELD
  IOFFS(J) = D%NSTAGT0B(J)*KFIELD
  ILENR(J) = D%NLTSGTB(J)*KFIELD
  IOFFR(J) = D%NSTAGT0B(D%MSTABF(J))*KFIELD
ENDDO

IF(NPROC > 1) THEN
  IF (LHOOK) CALL DR_HOOK('TRMTOL_BAR',0,ZHOOK_HANDLE_BAR)
  CALL GSTATS_BARRIER(764)
  IF (LHOOK) CALL DR_HOOK('TRMTOL_BAR',1,ZHOOK_HANDLE_BAR)
! IF (LSYNC_TRANS) THEN
!   CALL MPL_BARRIER(CDSTRING='TRMTOL')
! ENDIF

  CALL GSTATS(807,0)
  CALL MPL_ALLTOALLV(PSENDBUF=PFBUF_IN,KSENDCOUNTS=ILENS,&
   & PRECVBUF=PFBUF,KRECVCOUNTS=ILENR,KSENDDISPL=IOFFS,KRECVDISPL=IOFFR,&
   & KCOMM=MPL_ALL_MS_COMM,CDSTRING='TRMTOL:')
!Faster on Cray - because of peculiarity of their MPICH
! CALL MPL_ALLTOALLV(PSENDBUF=PFBUF_IN,KSENDCOUNTS=ILENS,&
!  & PRECVBUF=PFBUF,KRECVCOUNTS=ILENR,KSENDDISPL=IOFFS,KRECVDISPL=IOFFR,&
!  & KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=IREQ,&
!  & KCOMM=MPL_ALL_MS_COMM,CDSTRING='TRMTOL:')
! CALL MPL_WAIT(KREQUEST=IREQ,CDSTRING='TRMTOL: WAIT')

  CALL GSTATS(807,1)
  IF (LHOOK) CALL DR_HOOK('TRMTOL_BAR2',0,ZHOOK_HANDLE_BAR2)
  CALL GSTATS_BARRIER2(764)
  IF (LHOOK) CALL DR_HOOK('TRMTOL_BAR2',1,ZHOOK_HANDLE_BAR2)
ELSE
  ILEN = D%NLTSGTB(MYSETW)*KFIELD
  ISTA = D%NSTAGT0B(MYSETW)*KFIELD+1
  CALL GSTATS(1608,0)
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(J)
  DO J=ISTA,ISTA+ILEN-1
    PFBUF(J) = PFBUF_IN(J)
  ENDDO
!$OMP END PARALLEL DO
  CALL GSTATS(1608,1)
ENDIF


IF (LHOOK) CALL DR_HOOK('TRMTOL',1,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

END SUBROUTINE TRMTOL
END MODULE TRMTOL_MOD
