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

SUBROUTINE SLCOMM_OBS(YDSL,KFLDSLB,LDINC,PB1A)

!****-------------------------------------------------------------------
!**** *SLCOMM_OBS* - Interpolation buffer communication
!                 Particular version used in SL scheme when "on demand"
!                 communications are activated.
!****-------------------------------------------------------------------
!     Purpose.
!     --------
!           This routine is called to perform the inter-processor 
!     communication required for schemes needing horizontal interpolations.

!**   Interface.
!     ----------
!        *CALL* *SLCOMM_OBS

!        Explicit arguments :
!        --------------------
!        INPUT:
!          YDSL           - SL_STRUCT definition
!          KFLDSLB        - number of fields in interpolation buffer.
!          LDINC          - increment flag for interpolation buffer
!                           communicate with a specific processor.
!        INPUT/OUTPUT:
!          PB1A           - interpolation buffer.

!        Implicit arguments :  None.
!        --------------------

!     Method.
!     -------
!     The method used is to send and receive data in blocks of a fixed
!     number of grid columns each of size KFLDSLB words. In this method
!     we will only send the next block to a particular processor if 
!     there is no outstanding block to receive from that processor. 
!     We adopt this approach to both minimize the message buffer 
!     requirements for interpolation buffer communication and also to perform
!     dynamic communication. Dynamic communication simply means that we
!     only continue communicating with those processors that are in 
!     SLCOMM_OBS and are sending us data. By probing for messages we avoid 
!     waiting for any processor on a receive.

!     FORMAT of message:
!     (1)   first point number (index into NSLSENDPOS table)
!     (2)   last  point number (index into NSLSENDPOS table)
!     (3)   first fixed field  (index into PB1A field dimension)
!     (4)   last  fixed field  (index into PB1A field dimension)
!     (5)   position of dynamic table (located after the field data)
!           = 0 if no dynamic table

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      MPP Group *ECMWF*
!      Original : 97-06-04

!     Modifications.
!     --------------
!      01-11-23  Deborah Salmond  LIMP_NOOLAP Option for non-overlapping 
!                                 message passing and buffer packing
!      26-05-02  Deborah Salmond  Memory save for LIMP_NOOLAP
!      02-10-01  G.Mozdzynski     support for radiation on-demand comms
!      03-01-07  Deborah Salmond  Tidy up
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      G.Mozdzynski  10-Nov-2005 Optimised OpenMP parallelisation
!      G.Mozdzynski  01-Jan-2008 Cleanup
!      Y. Seity  : add barrier under LSYNC_SLCOM key
!      K. Yessad (Oct 2008) : unified SLCOMM (=merge SLCOMM+SLCOMM1) + clean.
!      G.Mozdzynski  02-Jan-2009 use non-blocking recv and send
!      G. Mozdzynski (Jan 2011): OOPS cleaning, use of derived type SL_STRUCT
!      R. El Khatib 08-Jun-2016 Optimize : use automatic arrays and overlap send/recv with pack/unpack
!     ------------------------------------------------------------------

USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN_IFSAUX, ONLY : NULOUT
USE MPL_MODULE   , ONLY : MPL_SEND, MPL_RECV, MPL_WAIT, MPL_BARRIER, &
 &                        JP_NON_BLOCKING_STANDARD, MPL_WAITANY
! arp/ifs dependencies to be solved later.
USE YOMMP0       , ONLY : NPRCIDS, LSYNC_SLCOM
USE YOMTAG       , ONLY : MTAGSLAG
USE YOMVAR       , ONLY : LSLADREP
USE EINT_MOD     , ONLY : SL_STRUCT

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(SL_STRUCT)   ,INTENT(IN)    :: YDSL
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDSLB 
LOGICAL           ,INTENT(IN)    :: LDINC 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PB1A(KFLDSLB,YDSL%NASLB1) 

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: ISENDTOT (YDSL%NSLPROCS)
INTEGER(KIND=JPIM) :: IRECVTOT (YDSL%NSLPROCS)
INTEGER(KIND=JPIM) :: ISEND_POS(YDSL%NSLPROCS)
INTEGER(KIND=JPIM) :: IRECV_POS(YDSL%NSLPROCS)
INTEGER(KIND=JPIM) :: I_MSENDTOT, I_MRECVTOT
INTEGER(KIND=JPIM) :: INUMFIELDS, IJR, IJS, J

LOGICAL            :: LLSWITCH

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SLCOMM_OBS',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

!     1/ Preliminary calculations.


  
!     2/ Mark all processors as having no outstanding receives,
!        nothing sent or received.

IJR=0
IJS=0
I_MSENDTOT=0
I_MRECVTOT=0
INUMFIELDS=KFLDSLB


DO J=1,YDSL%NSLPROCS
  ISENDTOT(J)=YDSL%NSLSENDPTR(J+1)-YDSL%NSLSENDPTR(J)
  IRECVTOT(J)=YDSL%NSLRECVPTR(J+1)-YDSL%NSLRECVPTR(J)
ENDDO
DO J=1,YDSL%NSLPROCS
  IF(ISENDTOT(J) > 0 ) THEN
    IJS=IJS+1
    ISEND_POS(J)=I_MSENDTOT*INUMFIELDS
    I_MSENDTOT=I_MSENDTOT+ISENDTOT(J)+1
  ENDIF
  IF(IRECVTOT(J) > 0 ) THEN
    IJR=IJR+1
    IRECV_POS(J)=I_MRECVTOT*INUMFIELDS
    I_MRECVTOT=I_MRECVTOT+IRECVTOT(J)+1
  ENDIF
ENDDO

CALL SLCOMM_INT(YDSL,YDSL%NSLRPT,YDSL%NSLSPT,YDSL%NSLRECVPOS,YDSL%NSLSENDPOS,&
 & YDSL%NSLRECVPTR,YDSL%NSLSENDPTR,&
 & KFLDSLB,LDINC,PB1A)


IF (LHOOK) CALL DR_HOOK('SLCOMM_OBS',1,ZHOOK_HANDLE)

CONTAINS

SUBROUTINE SLCOMM_INT(YDSL,KSLRPT,KSLSPT,KSLRECVPOS,KSLSENDPOS,&
 & KSLRECVPTR,KSLSENDPTR,KFLDSLB,LDINC,PB1A)

USE PARKIND1 , ONLY : JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
TYPE(SL_STRUCT),   INTENT(IN)    :: YDSL
INTEGER(KIND=JPIM),INTENT(IN)    :: KSLRPT
INTEGER(KIND=JPIM),INTENT(IN)    :: KSLSPT
INTEGER(KIND=JPIM),INTENT(IN)    :: KSLRECVPOS(KSLRPT)
INTEGER(KIND=JPIM),INTENT(IN)    :: KSLSENDPOS(KSLSPT)
INTEGER(KIND=JPIM),INTENT(IN)    :: KSLRECVPTR(YDSL%NSLPROCS+1)
INTEGER(KIND=JPIM),INTENT(IN)    :: KSLSENDPTR(YDSL%NSLPROCS+1)
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDSLB 
LOGICAL           ,INTENT(IN)    :: LDINC 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PB1A(KFLDSLB,YDSL%NASLB1) 

REAL(KIND=JPRB) :: ZSENDBUF(0:I_MSENDTOT*INUMFIELDS)
REAL(KIND=JPRB) :: ZRECVBUF(0:I_MRECVTOT*INUMFIELDS)

INTEGER(KIND=JPIM) :: IREQ_RCV (YDSL%NSLPROCS)
INTEGER(KIND=JPIM) :: IREQINDX (YDSL%NSLPROCS)
INTEGER(KIND=JPIM) :: IREQ_SND (YDSL%NSLPROCS)
INTEGER(KIND=JPIM) :: IJPOS    (YDSL%NSLPROCS)
INTEGER(KIND=JPIM) :: ISEND    (YDSL%NSLPROCS)

INTEGER(KIND=JPIM) :: IBEG, IEND, ILEN, IRECVPROC, ISENDPROC, ITAG, IRRCV, IRSND
INTEGER(KIND=JPIM) :: J, JJ, JFLD, I_KPOS, INR, JNR

LOGICAL :: LLCOMMDBGL
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('SLCOMM_OBS:SLCOMM_INT',0,ZHOOK_HANDLE)

!     3/ Determine the maximum number of grid columns that can be communicated.

ITAG=MTAGSLAG

LLCOMMDBGL=.FALSE.

IF( LLCOMMDBGL )THEN
  DO J=1,YDSL%NSLPROCS
    WRITE(NULOUT,'("SLCOMM: J,YDSL%NSLCOMM(J),ISENDTOT(J),",&
     & "IRECVTOT(J)=",4I6)')J,YDSL%NSLCOMM(J),ISENDTOT(J),IRECVTOT(J)  
    CALL FLUSH(NULOUT)
  ENDDO
ENDIF

!     4/ post recv loop.

CALL GSTATS_BARRIER(760)

CALL GSTATS(509,0)

IF (LSYNC_SLCOM) THEN
  CALL MPL_BARRIER(CDSTRING='SLCOMM_OBS:SLCOMM_INT')
ENDIF

IRRCV=0
DO J=1,YDSL%NSLPROCS
  IF( IRECVTOT(J) > 0 )THEN
    IRRCV=IRRCV+1
    IREQINDX(IRRCV)=J
    ISENDPROC=YDSL%NSLCOMM(J)
    IF( LLCOMMDBGL )THEN
      WRITE(NULOUT,'("SLCOMM: RECEIVING FROM ",I6," ITAG=",I9)')ISENDPROC,ITAG  
      CALL FLUSH(NULOUT)
    ENDIF
    ILEN=IRECVTOT(J)*INUMFIELDS
    CALL MPL_RECV(ZRECVBUF(IRECV_POS(J):IRECV_POS(J)+ILEN),&
     & KSOURCE=NPRCIDS(ISENDPROC),&
     & KMP_TYPE=JP_NON_BLOCKING_STANDARD,&
     & KREQUEST=IREQ_RCV(IRRCV),&
     & KTAG=ITAG,&
     & CDSTRING='SLCOMM_OBS:SLCOMM_INT RECV')
  ENDIF
ENDDO

CALL GSTATS(509,1)
CALL GSTATS_BARRIER2(760)


!     5/ Pack and send loop.

CALL GSTATS(1115,0)

IRSND=0
DO J=1,YDSL%NSLPROCS
  IF ( ISENDTOT(J) > 0 )THEN
    IRECVPROC=YDSL%NSLCOMM(J)
    IBEG=KSLSENDPTR(J)
    IEND=KSLSENDPTR(J+1)-1
    ILEN=IEND-IBEG+1
    IJPOS(J)=KFLDSLB*ILEN
    ISEND(J)=NPRCIDS(IRECVPROC)
    ZSENDBUF(ISEND_POS(J))=ISENDTOT(J)
    IF( LLCOMMDBGL )THEN
      WRITE(NULOUT,'("SLCOMM: SENDING ",I6," TO ",I6," ITAG=",I9)')ISENDTOT(J),IRECVPROC,ITAG  
      CALL FLUSH(NULOUT)
    ENDIF
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JFLD,I_KPOS,JJ)
    DO JJ=IBEG,IEND
      DO JFLD=1,KFLDSLB
        I_KPOS=ISEND_POS(J)+(JFLD-1)*ILEN
        ZSENDBUF(I_KPOS+JJ-IBEG+1)=PB1A(JFLD,KSLSENDPOS(JJ))
      ENDDO
    ENDDO
!$OMP END PARALLEL DO
    IRSND=IRSND+1
    CALL MPL_SEND(ZSENDBUF(ISEND_POS(J):ISEND_POS(J)+IJPOS(J)),&
     & KDEST=ISEND(J),&
     & KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=IREQ_SND(IRSND),&
     & KTAG=ITAG,CDSTRING='SLCOMM_OBS:SLCOMM_INT SEND' )
  ENDIF
ENDDO

CALL GSTATS(1115,1)



!     6/ Receive and Unpack loop:

!     Extract communication buffer data and place in right position of PB1A.
!     Polar mimics are handled by the routine (E)SLEXTPOL.

CALL GSTATS(1116,0)


IF (LDINC) THEN
  ! Do not use mpi_waitany because we want the order of summation to remain the same for bitwise reproducibility reason
  DO JNR=1,IRRCV
    CALL MPL_WAIT(KREQUEST=IREQ_RCV(JNR),CDSTRING='SLCOMM_OBS:SLCOMM_INT WAIT FOR RECVS')  
    J=IREQINDX(JNR)
    IBEG=KSLRECVPTR(J)
    IEND=KSLRECVPTR(J)+ZRECVBUF(IRECV_POS(J))-1
    ILEN=IEND-IBEG+1
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JFLD,I_KPOS,JJ)
    DO JJ=IBEG,IEND
      DO JFLD=1,KFLDSLB
        I_KPOS=IRECV_POS(J)+(JFLD-1)*ILEN
        PB1A(JFLD,KSLRECVPOS(JJ))=PB1A(JFLD,KSLRECVPOS(JJ))+ZRECVBUF(I_KPOS+JJ-IBEG+1)
      ENDDO
    ENDDO
!$OMP END PARALLEL DO
  ENDDO
ELSE
  DO JNR=1,IRRCV
    CALL MPL_WAITANY(KREQUEST=IREQ_RCV(1:IRRCV),KINDEX=INR,CDSTRING='SLCOMM_OBS:SLCOMM_INT: WAIT FOR ANY RECV')
    J=IREQINDX(INR)
    IBEG=KSLRECVPTR(J)
    IEND=KSLRECVPTR(J)+ZRECVBUF(IRECV_POS(J))-1
    ILEN=IEND-IBEG+1
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JFLD,I_KPOS,JJ)
    DO JJ=IBEG,IEND
      DO JFLD=1,KFLDSLB
        I_KPOS=IRECV_POS(J)+(JFLD-1)*ILEN
        PB1A(JFLD,KSLRECVPOS(JJ))=ZRECVBUF(I_KPOS+JJ-IBEG+1)
      ENDDO
    ENDDO
!$OMP END PARALLEL DO
  ENDDO
ENDIF

CALL GSTATS(1116,1)

!     7/ Wait for all outstanding sends.

IF(IRSND > 0) THEN
  CALL MPL_WAIT(KREQUEST=IREQ_SND(1:IRSND),CDSTRING='SLCOMM_OBS:SLCOMM_INT WAIT FOR SENDS')
ENDIF

IF (LHOOK) CALL DR_HOOK('SLCOMM_OBS:SLCOMM_INT',1,ZHOOK_HANDLE)
END SUBROUTINE SLCOMM_INT

END SUBROUTINE SLCOMM_OBS

