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

SUBROUTINE SLCOMM2A(YDDIM,YDSL,KFIXFLD,LDEXCLUDE,KMASK2,KFLDSLB,LDINC,PB1A,LDSLCOMM,LDTLEVOL)
!****-------------------------------------------------------------------
!**** *SLCOMM2A* - Interpolation buffer communication
!****-------------------------------------------------------------------
!     Purpose.
!     --------
!           This routine is called to perform the inter-processor 
!     communication required for schemes needing horizontal interpolations.

!**   Interface.
!     ----------
!        *CALL* *SLCOMM2A

!        Explicit arguments :
!        --------------------
!        INPUT:
!          YDSL           - SL_STRUCT definition
!          KFLDSLB        - number of fields in interpolation buffer.
!          KFIXFLD        - description of 'fixed' fields. 
!          LDEXCLUDE      - .true. - exclude fields in range specified by KFIXFLD
!                           .false.- include fields in range specified by KFIXFLD
!          KMASK2         - Masks for on demand-comms
!          LDINC          - increment flag for interpolation buffer
!        INPUT/OUTPUT:
!          PB1A           - interpolation buffer.
!          LDSLCOMM       - (optional)flag indicating whether this processor needs to
!                           communicate with a specific processor.

!        Implicit arguments :  None.
!        --------------------

!     Method.
!     -------
!       Copy of SLCOMM2 for 'on demand' Semi-Lagrangian Communications
!       Special version which only communicates data necessary for stencil

!     Externals.
!     ----------
!        MPL_BARRIER,  MPL_SEND, MPL_RECV, MPL_WAIT
!                                   - Communication Interface

!        Is called by SCAN2M, SCAN2MTL, CPCLIMO, CPCLIMI, FPMODPREC

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Deborah Salmond and John Hague

!     Modifications
!     -------------
!        01-11-23  Deborah Salmond  LIMP_NOOLAP Option for non-overlapping 
!                                    message passing and buffer packing
!        26-05-02  Deborah Salmond  Memory save for LIMP_NOOLAP
!        02-10-01  G.Mozdzynski     support for radiation on-demand comms
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        G.Mozdzynski  10-Nov-2005 Optimised OpenMP parallelisation
!        M.Jidane 09-04-2006 : Cleaning CY31
!        G.Mozdzynski  01-Jan-2008 Cleanup
!        G.Mozdzynski  02-Jan-2009 use non-blocking recv and send
!        G.Mozdzynski  27-Jul-2009 Optimise SL communications (LDSLCOMM)
!        G.Mozdzynski (Jan 2011): OOPS cleaning, use of derived type SL_STRUCT
!        G. Mozdzynski (Aug 2011): support higher order interpolation
!        R. El Khatib 13-Jun-2016 Optimize (automatic arrays + overlap send/recv with pack/unpack)
!     ------------------------------------------------------------------

USE YOMDIM   , ONLY : TDIM
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN_IFSAUX, ONLY : NULOUT
USE MPL_MODULE
! arp/ifs dependencies to be solved later.
USE YOMMP0   , ONLY : NPRCIDS, LSLDEBUG
USE YOMTAG   , ONLY : MTAGSLAG
USE YOMVAR   , ONLY : LSLADREP
USE ALGORITHM_STATE_MOD  , ONLY : GET_NSIM4D
USE EINT_MOD , ONLY : SL_STRUCT

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TDIM)        ,INTENT(IN)    :: YDDIM
TYPE(SL_STRUCT)   ,INTENT(INOUT) :: YDSL
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDSLB 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFIXFLD(2) 
LOGICAL           ,INTENT(IN)    :: LDEXCLUDE 
INTEGER(KIND=JPIM),INTENT(IN)    :: KMASK2 (YDSL%NASLB1+YDDIM%NSTENCILWIDE*2) 
LOGICAL           ,INTENT(IN)    :: LDINC 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PB1A(YDSL%NASLB1,KFLDSLB) 
LOGICAL,INTENT(INOUT), OPTIONAL  :: LDSLCOMM(YDSL%NSLPROCS)
LOGICAL,INTENT(IN)   , OPTIONAL  :: LDTLEVOL

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: IRECVMAP(0:YDSL%NSLRPT,YDSL%NSLPROCS)
INTEGER(KIND=JPIM) :: ISENDMAP(0:YDSL%NSLSPT,YDSL%NSLPROCS)
INTEGER(KIND=JPIM) :: ISEND_POS(YDSL%NSLPROCS), IRECV_POS(YDSL%NSLPROCS)
INTEGER(KIND=JPIM) :: IREQ_RCV(YDSL%NSLPROCS), IREQ_SND(YDSL%NSLPROCS), IREQINDX (YDSL%NSLPROCS)

INTEGER(KIND=JPIM) :: INUMFIELDS, ISENDCOUNT, IRECVCOUNT
INTEGER(KIND=JPIM) :: ISENDTOT (YDSL%NSLPROCS), IRECVTOT (YDSL%NSLPROCS)
INTEGER(KIND=JPIM) :: I_MAPSENDTOT, I_MAPRECVTOT
INTEGER(KIND=JPIM) :: J, IBEG, IEND, ILEN, IPT, IRECVPROC, ISENDPROC, ITAG, JJ, JNUM, ITOT
INTEGER(KIND=JPIM) :: IRRCV, IRSND

LOGICAL :: LLSLCOMM (YDSL%NSLPROCS)
LOGICAL :: LLCOMMDBGL

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "check_sl_struct.intfb.h"
#include "abor1.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SLCOMM2A',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

CALL GSTATS(1813,0)

ITAG=MTAGSLAG

DO J=1,YDSL%NSLPROCS
  IF(PRESENT(LDSLCOMM))THEN
    IF (.NOT. PRESENT(LDTLEVOL))THEN
      CALL ABOR1('Calling SLCOMM2A with LDSLCOMM but without LDTLEVOL: not permitted') 
    ENDIF
    IF((GET_NSIM4D() > 0).AND.(.NOT.LDTLEVOL)) THEN
      LLSLCOMM(J)=LDSLCOMM(J)
    ELSE
      LLSLCOMM(J)=.TRUE.
    ENDIF
    IF( LSLADREP ) LLSLCOMM(J)=.TRUE.
  ELSE
    LLSLCOMM(J)=.TRUE.
  ENDIF
ENDDO

IF(LDINC)THEN
  CALL ABOR1('SLCOMM2A: NOT DONE FOR LDINC')
ENDIF
CALL GSTATS(1813,1)

IF(LDEXCLUDE)THEN
  INUMFIELDS=KFLDSLB-(KFIXFLD(2)-KFIXFLD(1)+1)
ELSE
  INUMFIELDS=KFIXFLD(2)-KFIXFLD(1)+1
ENDIF

LLCOMMDBGL=.FALSE.

!     Communicate maps

CALL GSTATS_BARRIER(758)
CALL GSTATS(502,0)

IRRCV=0
DO J=1,YDSL%NSLPROCS
  ISENDTOT(J)=YDSL%NSLSENDPTR(J+1)-YDSL%NSLSENDPTR(J)
  ISENDMAP(0,J)=0
  IF( LLSLCOMM(J) )THEN
    IRRCV=IRRCV+1
    IREQINDX(IRRCV)=J
    IRECVPROC=YDSL%NSLCOMM(J)
    ILEN=ISENDTOT(J)
    CALL MPL_RECV(ISENDMAP(0:ILEN,J),KSOURCE=NPRCIDS(IRECVPROC),KTAG=ITAG,&
     & KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=IREQ_RCV(IRRCV),&
     & CDSTRING='SLCOMM2A:SLCOMM2A RECV MAP ' )  
  ENDIF
ENDDO

IRSND=0
DO J=1,YDSL%NSLPROCS
  ITOT=0
  IBEG=YDSL%NSLRECVPTR(J)
  IEND=YDSL%NSLRECVPTR(J+1)-1
  DO JNUM=IBEG,IEND
    IPT=YDSL%NSLRECVPOS(JNUM)
    IRECVMAP(JNUM-IBEG+1,J)= KMASK2(IPT)
    IF(KMASK2(IPT) == 1) ITOT=ITOT+1
  ENDDO
  IRECVTOT(J)=YDSL%NSLRECVPTR(J+1)-YDSL%NSLRECVPTR(J)
  IRECVMAP(0,J)=ITOT
  IF(IRECVMAP(0,J) == 0) IRECVTOT(J)=0 
  IF( LLSLCOMM(J) )THEN
    IRSND=IRSND+1
    ISENDPROC=YDSL%NSLCOMM(J)
    ILEN=IRECVTOT(J)
    CALL MPL_SEND(IRECVMAP(0:ILEN,J),KDEST=NPRCIDS(ISENDPROC),KTAG=ITAG,&
     & KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=IREQ_SND(IRSND),&
     & CDSTRING='SLCOMM2A:SLCOMM2A SEND MAP ' )  
  ENDIF
ENDDO
IRECV_POS(1)=0
DO J=2,YDSL%NSLPROCS
  IRECV_POS(J)=IRECV_POS(J-1)+IRECVMAP(0,J-1)*INUMFIELDS+1
ENDDO
I_MAPRECVTOT=IRECV_POS(YDSL%NSLPROCS)+IRECVMAP(0,YDSL%NSLPROCS)*INUMFIELDS

IF( LLCOMMDBGL )THEN
  DO J=1,YDSL%NSLPROCS
    WRITE(NULOUT,'("SLCOMM2A: J,YDSL%NSLCOMM(J),ISENDTOT(J),",&
     & "IRECVTOT(J)=",4I6)')J,YDSL%NSLCOMM(J),ISENDTOT(J),IRECVTOT(J)
    CALL FLUSH(NULOUT)
  ENDDO
ENDIF

ISEND_POS(1)=0
IF( IRRCV > 0 )THEN
  CALL MPL_WAIT(KREQUEST=IREQ_RCV(1:IRRCV),CDSTRING='SLCOMM2A:SLCOMM2A WAIT FOR MAP RECEIVES')  
ENDIF
DO J=1,YDSL%NSLPROCS
  IF(ISENDMAP(0,J) == 0) ISENDTOT(J)=0 
ENDDO
DO J=2,YDSL%NSLPROCS
  ISEND_POS(J)=ISEND_POS(J-1)+ISENDMAP(0,J-1)*INUMFIELDS+1
ENDDO
I_MAPSENDTOT=ISEND_POS(YDSL%NSLPROCS)+ISENDMAP(0,YDSL%NSLPROCS)*INUMFIELDS

ISENDCOUNT=MAXVAL(ISENDTOT(:))
IRECVCOUNT=MAXVAL(IRECVTOT(:))

IF( PRESENT(LDSLCOMM) .AND. PRESENT(LDTLEVOL)) THEN
  IF ( (GET_NSIM4D()==0 .OR. LDTLEVOL) .AND. .NOT.LSLADREP )THEN
    DO J=1,YDSL%NSLPROCS
      LDSLCOMM(J)=ISENDTOT(J)>0 .OR. IRECVTOT(J)>0
    ENDDO
  ENDIF
ENDIF

IF( IRSND > 0 )THEN
  CALL MPL_WAIT(KREQUEST=IREQ_SND(1:IRSND),CDSTRING='SLCOMM2A:SLCOMM2A WAIT FOR MAP SENDS')  
ENDIF

CALL GSTATS(502,1)
CALL GSTATS_BARRIER2(758)

!    Communicate halos

CALL SLCOMM2A_INT

IF (LSLDEBUG) CALL CHECK_SL_STRUCT(YDSL,'SLCOMM2A')

IF (LHOOK) CALL DR_HOOK('SLCOMM2A',1,ZHOOK_HANDLE)


CONTAINS

SUBROUTINE SLCOMM2A_INT

USE PARKIND1 , ONLY : JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK

INTEGER(KIND=JPIM) :: IPTA_SND(ISENDCOUNT)
INTEGER(KIND=JPIM) :: IPTA_RCV(IRECVCOUNT)
INTEGER(KIND=JPIM) :: I_KPOSJ (KFLDSLB)

INTEGER(KIND=JPIM) :: IBEG, IEND, ILEN, IPT, IRECVPROC, ISENDPROC, ITAG, J, JJ, JFLD, JNUM, ITOT, I_KPOS
INTEGER(KIND=JPIM) :: IRRCV, IRSND , JM, JNR, INR, JPOS

! SS/2018: Reducing stack pressure

REAL(KIND=JPRB), ALLOCATABLE :: ZSENDBUF(:) ! 0:I_MAPSENDTOT
REAL(KIND=JPRB), ALLOCATABLE :: ZRECVBUF(:) ! 0:I_MAPRECVTOT

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE


IF (LHOOK) CALL DR_HOOK('SLCOMM2A:SLCOMM2A_INT',0,ZHOOK_HANDLE)

ALLOCATE(ZSENDBUF(0:I_MAPSENDTOT))
ALLOCATE(ZRECVBUF(0:I_MAPRECVTOT))

ITAG=MTAGSLAG+1

!       Post Receive loop.........................................................

CALL GSTATS_BARRIER(749)
CALL GSTATS(512,0)

IRRCV=0
DO J=1,YDSL%NSLPROCS
  IF( IRECVTOT(J) > 0 .AND. LLSLCOMM(J) )THEN
    IRRCV=IRRCV+1
    IREQINDX(IRRCV)=J
    ISENDPROC=YDSL%NSLCOMM(J)
    IF( LLCOMMDBGL )THEN
      WRITE(NULOUT,'("SLCOMM2A: RECEIVING ",I6," FROM ",I6,&
       & " ITAG=",I9)')ISENDPROC,ITAG  
      CALL FLUSH(NULOUT)
    ENDIF
    ILEN=IRECVMAP(0,J)*INUMFIELDS
    CALL MPL_RECV(ZRECVBUF(IRECV_POS(J):IRECV_POS(J)+ILEN),&
     & KSOURCE=NPRCIDS(ISENDPROC),&
     & KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=IREQ_RCV(IRRCV),&
     & KTAG=ITAG,CDSTRING='SLCOMM2A:SLCOMM2A_INT RECV' )  
  ENDIF
ENDDO

CALL GSTATS(512,1)
CALL GSTATS_BARRIER2(749)

!       Pack & send  loop.............................................................

CALL GSTATS(1121,0)

IRSND=0
DO J=1,YDSL%NSLPROCS
  IF( ISENDTOT(J) > 0 .AND. LLSLCOMM(J) )THEN
    IBEG=YDSL%NSLSENDPTR(J)
    IEND=YDSL%NSLSENDPTR(J+1)-1
    JM=0
    DO JNUM=IBEG,IEND
      IF(ISENDMAP((JNUM-YDSL%NSLSENDPTR(J)+1),J) == 1)THEN
        JM=JM+1
        IPTA_SND(JM)=YDSL%NSLSENDPOS(JNUM)
      ENDIF
    ENDDO
    I_KPOS=0
    IF(LDEXCLUDE)THEN
      DO JFLD=1,KFLDSLB
        IF(JFLD < KFIXFLD(1).OR.JFLD > KFIXFLD(2)) THEN
          I_KPOSJ(JFLD)=I_KPOS
          I_KPOS=I_KPOS+JM
        ENDIF
      ENDDO
    ELSE
      DO JFLD=KFIXFLD(1),KFIXFLD(2)
        I_KPOSJ(JFLD)=I_KPOS
        I_KPOS=I_KPOS+JM
      ENDDO
    ENDIF
    JPOS=I_KPOS
    ZSENDBUF(ISEND_POS(J))=ISENDTOT(J)

!$OMP PARALLEL PRIVATE(JFLD,I_KPOS,JJ)
    IF(LDEXCLUDE)THEN
!$OMP DO SCHEDULE(STATIC)
      DO JFLD=1,KFIXFLD(1)-1
        I_KPOS=ISEND_POS(J)+I_KPOSJ(JFLD)
        DO JJ=1,JM
          ZSENDBUF(I_KPOS+JJ)=PB1A(IPTA_SND(JJ),JFLD)
        ENDDO
      ENDDO
!$OMP END DO NOWAIT
!$OMP DO SCHEDULE(STATIC)
      DO JFLD=KFIXFLD(2)+1,KFLDSLB
        I_KPOS=ISEND_POS(J)+I_KPOSJ(JFLD)
        DO JJ=1,JM
          ZSENDBUF(I_KPOS+JJ)=PB1A(IPTA_SND(JJ),JFLD)
        ENDDO
      ENDDO
!$OMP END DO
    ELSE
!$OMP DO SCHEDULE(STATIC)
      DO JFLD=KFIXFLD(1),KFIXFLD(2)
        I_KPOS=ISEND_POS(J)+I_KPOSJ(JFLD)
        DO JJ=1,JM
          ZSENDBUF(I_KPOS+JJ)=PB1A(IPTA_SND(JJ),JFLD)
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
!$OMP END PARALLEL
    IRECVPROC=YDSL%NSLCOMM(J)
    IF( LLCOMMDBGL )THEN
      WRITE(NULOUT,'("SLCOMM2A: SENDING ",I6," TO ",I6,&
       & " ITAG=",I9)')ISENDTOT(J),IRECVPROC,ITAG
      CALL FLUSH(NULOUT)
    ENDIF
    IRSND=IRSND+1
    CALL MPL_SEND(ZSENDBUF(ISEND_POS(J):ISEND_POS(J)+JPOS),&
     & KDEST=NPRCIDS(IRECVPROC),&
     & KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=IREQ_SND(IRSND),&
     & KTAG=ITAG,CDSTRING='SLCOMM2A:SLCOMM2A_INT SEND ' )  
  ENDIF
ENDDO

CALL GSTATS(1121,1)

!     Receive and Unpack loop.............................................................

CALL GSTATS(1118,0)

DO JNR=1,IRRCV
  CALL MPL_WAITANY(KREQUEST=IREQ_RCV(1:IRRCV),KINDEX=INR,CDSTRING='SLCOMM2A:SLCOMM2A_INT: WAIT FOR ANY RECV')
  J=IREQINDX(INR)
  IBEG=YDSL%NSLRECVPTR(J)
  IEND=YDSL%NSLRECVPTR(J)+ZRECVBUF(IRECV_POS(J))-1
! Extract communication buffer data and place in right position of PB1A.
! Polar mimics are handled by the routine slextpol.F
  JM=0
  DO JNUM=IBEG,IEND
    IF(IRECVMAP((JNUM-YDSL%NSLRECVPTR(J)+1),J) == 1)THEN
      JM=JM+1
      IPTA_RCV(JM)=YDSL%NSLRECVPOS(JNUM)
    ENDIF
  ENDDO
  I_KPOS=0
  IF(LDEXCLUDE)THEN
    DO JFLD=1,KFLDSLB
      IF(JFLD < KFIXFLD(1).OR.JFLD > KFIXFLD(2)) THEN
        I_KPOSJ(JFLD)=I_KPOS
        I_KPOS=I_KPOS+JM
      ENDIF
    ENDDO
  ELSE
    DO JFLD=KFIXFLD(1),KFIXFLD(2)
      I_KPOSJ(JFLD)=I_KPOS
      I_KPOS=I_KPOS+JM
    ENDDO
  ENDIF
  JPOS=I_KPOS
!$OMP PARALLEL PRIVATE(JFLD,I_KPOS,JJ)
  IF(LDEXCLUDE)THEN
!$OMP DO SCHEDULE(STATIC)
    DO JFLD=1,KFIXFLD(1)-1
      I_KPOS=IRECV_POS(J)+I_KPOSJ(JFLD)
      DO JJ=1,JM
        PB1A(IPTA_RCV(JJ),JFLD)=ZRECVBUF(I_KPOS+JJ)
      ENDDO
    ENDDO
!$OMP END DO NOWAIT
!$OMP DO SCHEDULE(STATIC)
    DO JFLD=KFIXFLD(2)+1,KFLDSLB
      I_KPOS=IRECV_POS(J)+I_KPOSJ(JFLD)
      DO JJ=1,JM
        PB1A(IPTA_RCV(JJ),JFLD)=ZRECVBUF(I_KPOS+JJ)
      ENDDO
    ENDDO
!$OMP END DO
  ELSE
!$OMP DO SCHEDULE(STATIC)
    DO JFLD=KFIXFLD(1),KFIXFLD(2)
      I_KPOS=IRECV_POS(J)+I_KPOSJ(JFLD)
      DO JJ=1,JM
        PB1A(IPTA_RCV(JJ),JFLD)=ZRECVBUF(I_KPOS+JJ)
      ENDDO
    ENDDO
!$OMP END DO
  ENDIF
!$OMP END PARALLEL
ENDDO

CALL GSTATS(1118,1)

!       Wait for all outstanding sends

IF(IRSND > 0) THEN
  CALL MPL_WAIT(KREQUEST=IREQ_SND(1:IRSND),CDSTRING='SLCOMM2A:SLCOMM2A_INT WAIT FOR SENDS AND RECEIVES')  
ENDIF

DEALLOCATE(ZRECVBUF)
DEALLOCATE(ZSENDBUF)

IF (LHOOK) CALL DR_HOOK('SLCOMM2A:SLCOMM2A_INT',1,ZHOOK_HANDLE)

END SUBROUTINE SLCOMM2A_INT

END SUBROUTINE SLCOMM2A
