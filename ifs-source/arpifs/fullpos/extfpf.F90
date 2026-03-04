! (C) Copyright 1989- Meteo-France.

SUBROUTINE EXTFPF(KNFG,KNFD,KFLDOFF,KBFA,KFPROMA,KFIELDS,KBLOCKS,KLAST,KGPG,KGPD,KGPIND,PGP,PREALG)

!**** *EXTFPF* - extracts gridpoint fields in the framework of Fullpos

!     Purpose.
!     --------
!       To transpose a gridpoint sliced and distributed set of fields to a 
!       field-distributed set of global gridpoint fields.

!**   Interface.
!     ----------
!        *CALL* *EXTFPF*

!        Explicit arguments :
!        --------------------

!        KNFG   : Total number of fields in the current chunk
!        KNFD   : number of fields in the current chunk for each MPI task (in output arrays)
!        KFLDOFF: offset to the first local field for each MPI task
!        KBFA   : Absolute first field to be extracted
!        KFPROMA: size of cache blocks (in input array)
!        KFIELDS: Total number of fields (in input array)
!        KBLOCKS: Number of cache blocks (in input array)
!        KLAST  : last rank of each block (in input array)
!        KGPG   : global number of gridpoints (in output array)
!        KGPD   : number of gridpoints for each MPI task (in input arrays)
!        KGPIND : KGPIND(i,j)=k means : i-th gridpoint on j-th processor is k-th global gp
!        PGP    : sliced and gridpoint-distributed fields (input)
!        PREALG : field-distributed gridpoints (output)

!        Implicit arguments :
!        --------------------
!           See modules above.

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
!      01-Dec-2015 R. El Khatib *METEO-FRANCE* merge and optimization of the former extfpf + diwrgrfp

!     Modifications.
!     --------------
!      R. El Khatib 12-Aug-2016 waitany
!     -------------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMMP0   , ONLY : MYPROC, NPRCIDS, NPROC
USE MPL_MODULE, ONLY : MPL_SEND, MPL_RECV, MPL_WAIT, JP_NON_BLOCKING_STANDARD, MPL_WAITANY

!     -------------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN)  :: KNFG 
INTEGER(KIND=JPIM), INTENT(IN)  :: KNFD(NPROC) 
INTEGER(KIND=JPIM), INTENT(IN)  :: KFLDOFF(NPROC) 
INTEGER(KIND=JPIM), INTENT(IN)  :: KBFA
INTEGER(KIND=JPIM), INTENT(IN)  :: KFPROMA
INTEGER(KIND=JPIM), INTENT(IN)  :: KFIELDS
INTEGER(KIND=JPIM), INTENT(IN)  :: KBLOCKS
INTEGER(KIND=JPIM), INTENT(IN)  :: KLAST(KBLOCKS)
INTEGER(KIND=JPIM), INTENT(IN)  :: KGPG 
INTEGER(KIND=JPIM), INTENT(IN)  :: KGPD(NPROC) 
INTEGER(KIND=JPIM), INTENT(IN)  :: KGPIND(MAXVAL(KGPD),NPROC) 
REAL(KIND=JPRB),    INTENT(IN)  :: PGP(KFPROMA,KFIELDS,KBLOCKS) 
REAL(KIND=JPRB),    INTENT(OUT) :: PREALG(KGPG,KNFD(MYPROC))

!     -------------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JROC, IRSND, IOFFPROC(NPROC), JI, JFLD, JFP, JBL
INTEGER(KIND=JPIM) :: ISIZERCV(NPROC) ! Size of received messages from each task
INTEGER(KIND=JPIM) :: IREQRCV(NPROC), IREQSND(NPROC) ! maximum possible number of recv/send requests
INTEGER(KIND=JPIM) :: IREQINDX(NPROC), JNR, INR, IRRCV

REAL(KIND=JPRB) :: ZSENDBUF(KGPD(MYPROC),KNFG) ! Send buffer
REAL(KIND=JPRB) :: ZRECVBUF(KGPG*KNFD(MYPROC)) ! Recv buffer

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     -------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('EXTFPF',0,ZHOOK_HANDLE)

!     -------------------------------------------------------------------------

!*     1. NON-BLOCKING PEER-TO-PEER COMMUNICATIONS
!         ----------------------------------------

!      1.1 Post recv first

IOFFPROC(1)=0
ISIZERCV(1)=KGPD(1)*KNFD(MYPROC)
DO JROC = 2, NPROC
  ISIZERCV(JROC)=KGPD(JROC)*KNFD(MYPROC)
  IOFFPROC(JROC)=IOFFPROC(JROC-1)+ISIZERCV(JROC-1)
ENDDO
IRRCV=0
DO JROC = 1, NPROC
  IF (ISIZERCV(JROC) > 0 .AND. JROC /= MYPROC) THEN
    IRRCV=IRRCV+1
    IREQINDX(IRRCV)=JROC
    CALL MPL_RECV(ZRECVBUF(IOFFPROC(JROC)+1:IOFFPROC(JROC)+ISIZERCV(JROC)), &
     & KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=IREQRCV(IRRCV), &
     & KOUNT=ISIZERCV(JROC),KSOURCE=NPRCIDS(JROC),KTAG=JROC,CDSTRING='EXTFPF:')  
  ENDIF
ENDDO

!      1.2 pack and send chunks of fields

IRSND=0
DO JROC = 1, NPROC
  IF (KGPD(MYPROC)*KNFD(JROC) > 0 .AND. JROC /= MYPROC) THEN
    IRSND=IRSND+1
!$OMP PARALLEL DO PRIVATE(JFLD,JBL,JFP) COLLAPSE(2)
    DO JFLD=KFLDOFF(JROC)+1,KFLDOFF(JROC)+KNFD(JROC)
      DO JBL=1,KBLOCKS
        DO JFP=1,KLAST(JBL)
          ZSENDBUF(JFP+KFPROMA*(JBL-1),JFLD)=PGP(JFP,KBFA-1+JFLD,JBL)
        ENDDO
      ENDDO
    ENDDO
!$OMP END PARALLEL DO
    CALL MPL_SEND(ZSENDBUF(:,KFLDOFF(JROC)+1:KFLDOFF(JROC)+KNFD(JROC)), &
     & KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=IREQSND(IRSND), &
     & KTAG=MYPROC,KDEST=NPRCIDS(JROC),CDSTRING='EXTFPF:')  
  ENDIF
ENDDO

!      1.3 Local contribution

!$OMP PARALLEL DO PRIVATE(JFLD,JBL,JFP) COLLAPSE(2)
DO JFLD=1,KNFD(MYPROC)
  DO JBL=1,KBLOCKS
    DO JFP=1,KLAST(JBL)
      PREALG(KGPIND(JFP+KFPROMA*(JBL-1),MYPROC),JFLD)=PGP(JFP,KBFA-1+KFLDOFF(MYPROC)+JFLD,JBL)
    ENDDO
  ENDDO
ENDDO
!$OMP END PARALLEL DO

!      1.4 Wait for receive and unpack

DO JNR=1,IRRCV
  CALL MPL_WAITANY(KREQUEST=IREQRCV(1:IRRCV),KINDEX=INR,CDSTRING='EXTFPF: WAIT FOR ANY RECVS')
  JROC=IREQINDX(INR)
!$OMP PARALLEL DO PRIVATE(JFLD,JI)
  DO JFLD=1,KNFD(MYPROC)
    DO JI=1,KGPD(JROC)
      PREALG(KGPIND(JI,JROC),JFLD)=ZRECVBUF(JI+(JFLD-1)*KGPD(JROC)+IOFFPROC(JROC))
    ENDDO
  ENDDO
!$OMP END PARALLEL DO
ENDDO

!      1.5 Wait completion of sends

IF (IRSND > 0) THEN
  CALL MPL_WAIT(KREQUEST=IREQSND(1:IRSND),CDSTRING='EXTFPF: WAIT FOR SENDS')
ENDIF

!     -------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('EXTFPF',1,ZHOOK_HANDLE)
END SUBROUTINE EXTFPF
