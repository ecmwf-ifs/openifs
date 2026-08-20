! (C) Copyright 2001- Aron Roland (Roland & Partner, Germany).
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE OUTPUT_STRUCT
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWSPHERE, only : SPHERICAL_COORDINATE_DISTANCE
      USE YOWCOUT, ONLY : JPPFLAG, IPFGTBL
      USE YOWPCONS , ONLY : ZMISS
      INTEGER(KIND=JWIM), dimension(:), pointer :: out_recv_rqst
      INTEGER(KIND=JWIM), dimension(:), pointer :: out_send_rqst
      INTEGER(KIND=JWIM), dimension(:,:), pointer :: out_recv_stat
      INTEGER(KIND=JWIM), dimension(:,:), pointer :: out_send_stat
      INTEGER(KIND=JWIM), dimension(:), pointer :: block_type
      INTEGER(KIND=JWIM) :: NbProcOut, NbProcOutRed
      INTEGER(KIND=JWIM) :: nbPointCovered
      TYPE CONTAINER_ARR
      real(KIND=JWRB), allocatable :: ARR(:,:)
      END TYPE CONTAINER_ARR
      TYPE(CONTAINER_ARR), allocatable :: ARR_OUT_SEND_C(:)
      real(KIND=JWRB), allocatable :: ARR_OUT_RECV(:,:)
      real(KIND=JWRB), allocatable :: ARR_OUT_RECV_SEP(:,:)
      INTEGER(KIND=JWIM), allocatable :: ListIEfind(:)
      INTEGER(KIND=JWIM), allocatable :: ListICTpos(:)
      INTEGER(KIND=JWIM), allocatable :: ListICTposRev(:)
      INTEGER(KIND=JWIM), allocatable :: LocalPosICT(:)
      INTEGER(KIND=JWIM), allocatable :: ListICT_to_ProcOut(:)
      INTEGER(KIND=JWIM), allocatable :: ListProcOut(:)
      INTEGER(KIND=JWIM), allocatable :: ListProcOutRev(:)
      INTEGER(KIND=JWIM), allocatable :: ListStartIDXout(:)
      INTEGER(KIND=JWIM), allocatable :: MapFD_locglob(:)
      INTEGER(KIND=JWIM), allocatable :: MapIPexp_IP(:)
      INTEGER(KIND=JWIM), allocatable :: NbEntries(:)
      INTEGER(KIND=JWIM), allocatable :: MAP_CovToExp(:,:)
      INTEGER(KIND=JWIM), allocatable :: IXarr(:), IYarr(:)
      real(KIND=JWRB), allocatable :: WIarr(:,:)
      logical, allocatable :: ApplyZMISS(:)
      INTEGER(KIND=JWIM) :: MPI_EXCH_STR
      INTEGER(KIND=JWIM) :: MaxLen
      INTEGER(KIND=JWIM) :: NIBLO_FD_EXP

PUBLIC :: SET_UP_ARR_OUT_RECV, &
     &    INITIAL_OUTPUT_INITS, &
     &    INITIAL_OUTPUT_INITS_NEXTGEN, &
     &    FIND_ELE, &
     &    INTELEMENT_IPOL, &
     &    IPELEMENT_CLOSEST, &
     &    INTELEMENT

INTERFACE SET_UP_ARR_OUT_RECV 
  MODULE PROCEDURE SET_UP_ARR_OUT_RECV
END INTERFACE

INTERFACE INITIAL_OUTPUT_INITS
  MODULE PROCEDURE INITIAL_OUTPUT_INITS
END INTERFACE

INTERFACE INITIAL_OUTPUT_INITS_NEXTGEN
  MODULE PROCEDURE INITIAL_OUTPUT_INITS_NEXTGEN
END INTERFACE

INTERFACE FIND_ELE 
  MODULE PROCEDURE FIND_ELE
END INTERFACE

INTERFACE INTELEMENT_IPOL
  MODULE PROCEDURE INTELEMENT_IPOL
END INTERFACE

INTERFACE IPELEMENT_CLOSEST  
  MODULE PROCEDURE IPELEMENT_CLOSEST 
END INTERFACE

INTERFACE INTELEMENT  
  MODULE PROCEDURE INTELEMENT
END INTERFACE

      CONTAINS
      SUBROUTINE SET_UP_ARR_OUT_RECV(IJS, IJL, TEMP, NFLDPPE)
      USE yowDatapool, only: comm
      USE YOWMPP   , ONLY : IRANK, NPROC
      USE YOW_RANK_GLOLOC, ONLY : MyRankGlobal
      USE yowpd, only: MNE=>ne, INE, MNP=>npa, NP_RES=>np
      USE YOWUNPOOL, only : NIBLO_FD
      USE yownodepool, only : iplg
      IMPLICIT NONE
      INTEGER(KIND=JWIM), intent(in) :: IJS, IJL
      REAL(KIND=JWRB), intent(IN) :: TEMP(IJS:IJL,JPPFLAG)
      INTEGER(KIND=JWIM), intent(in) :: NFLDPPE(nproc)
      INTEGER(KIND=JWIM) :: IP, eProcOut, fProcOut, iProcOut, jProcOut
      REAL(KIND=JWRB) :: eWI, eVal, eValIns
      INTEGER(KIND=JWIM) :: I, ICTpos, IEfind, idx
      INTEGER(KIND=JWIM) :: iStart, iEnd, pos, ierr, iProc
      INTEGER(KIND=JWIM) :: IPR, IJ
      INTEGER(KIND=JWIM) :: iNode, iNodeGlob, IPexp
      INTEGER(KIND=JWIM) :: eLen, iVar, NbEnt
      INTEGER(KIND=JWIM) :: nbZmiss
      logical, save :: FirstTime = .TRUE.
      INTEGER(KIND=JWIM) :: idxHS, posHS
      INTEGER(KIND=JWIM) :: IPglob
#ifdef DEBUG
      INTEGER(KIND=JWIM) :: iLen
#endif
#ifdef DEBUG
      WRITE(740+MyRankGlobal,*) 'JWRU=', JWRU
      WRITE(740+MyRankGlobal,*) 'JWRB=', JWRB
      WRITE(740+MyRankGlobal,*) 'NbProcOut=', NbProcOut
      WRITE(740+MyRankGlobal,*) 'IJS / IJL=', IJS, IJL
      FLUSH(740+MyRankGlobal)
#endif
      jProcOut=0
      DO iProcOut=1,NbProcOut
#ifdef DEBUG
        WRITE(740+MyRankGlobal,*) 'iProcOut=', iProcOut
        FLUSH(740+MyRankGlobal)
#endif
        eProcOut=ListProcOut(iProcOut)
#ifdef DEBUG
        WRITE(740+MyRankGlobal,*) 'eProcOut=', eProcOut
        FLUSH(740+MyRankGlobal)
#endif
        IF (eProcOut .ne. IRANK) THEN
          jProcOut = jProcOut + 1
          ARR_OUT_SEND_C(jProcOut) % ARR(:,:) = 0
#ifdef DEBUG
          WRITE(740+MyRankGlobal,*) 'Doing the exchange'
          WRITE(740+MyRankGlobal,*) 'IJS=', IJS
          FLUSH(740+MyRankGlobal)
#endif
          iStart=ListStartIDXout(iProcOut)
          iEnd=ListStartIDXout(iProcOut+1)-1
          eLen=iEnd + 1 - iStart
          DO iNode=1,nbPointCovered
            iNodeGlob = MapFD_locglob(iNode)
            IEfind=ListIEfind(iNode)
            ApplyZMISS(:)=.FALSE.
            DO I=1,3
              IP=INE(I,IEfind)
              IF (IP .le. NP_RES) THEN
                IJ= IJS - 1 + IP
                DO idx=iStart,iEnd
                  ICTpos=ListICTpos(idx)
                  pos=idx + 1 - iStart
                  eWI=WIarr(I,iNode)
                  IF (eWI .gt. 0) THEN
                    eValIns=TEMP(IJ,ICTpos)
#ifdef DEBUG
                    fProcOut = ListICT_to_ProcOut(1)
                    idxHS = ListICTposRev(1)
                    posHS = idxHS + 1 - iStart
                    IF ((MapIPexp_IP(iNodeGlob) .eq. 211).and.(fProcOut .eq. eProcOut).and.(posHS.eq.pos)) THEN
                      IPglob=iplg(IP)
                      WRITE(740+MyRankGlobal,*) 'IP=', 211, ' 1: IPexp=', iNodeGlob
                      WRITE(740+MyRankGlobal,*) 'I=', I
                      WRITE(740+MyRankGlobal,*) 'IJ=', IJ
                      WRITE(740+MyRankGlobal,*) 'eWI=', eWI, ' IPglob=', IPglob
                      WRITE(740+MyRankGlobal,*) 'eValIns=', eValIns, ' IPglob=', IPglob
                      FLUSH(740+MyRankGlobal)
                    END IF
#endif
                    IF (eValIns .eq. ZMISS) THEN
                      ApplyZMISS(pos)=.TRUE.
                    ENDIF
                    eVal = ARR_OUT_SEND_C(jProcOut) % ARR(pos,iNode)
                    eVal = eVal + eWI*eValIns
                    ARR_OUT_SEND_C(jPRocOut) % ARR(pos,iNode) = eVal
                  END IF
                END DO
              END IF
            END DO
#ifdef DEBUG
            fProcOut = ListICT_to_ProcOut(10)
            IF ((iNode .eq. 1187).and.(fProcOut.eq.eProcOut)) THEN
              idx = ListICTposRev(10)
              pos = idx + 1 - iStart
              WRITE(740+MyRankGlobal,*) 'U10NEW: spec 1187 = ', ARR_OUT_SEND_C(jPRocOut) % ARR(pos,iNode)
            END IF
#endif
            DO pos=1,eLen
              IF (ApplyZMISS(pos)) THEN
                ARR_OUT_SEND_C(jProcOut) % ARR(pos,iNode) = ZMISS
              END IF
            END DO
          END DO
#ifdef DEBUG
          fProcOut = ListICT_to_ProcOut(10)
          WRITE(740+MyRankGlobal,*) 'fProcOut=', fProcOut
          WRITE(740+MyRankGlobal,*) 'ZMISS=', ZMISS
          WRITE(740+MyRankGlobal,*) 'iStart=', iStart, ' iEnd=', iEnd
          WRITE(740+MyRankGlobal,*) 'iProcOut=', iProcOut
          IF (fProcOut .eq. eProcOut) THEN
            idx = ListICTposRev(10)
            ICTpos=ListICTpos(idx)
            pos = idx + 1 - iStart
            WRITE(740+MyRankGlobal,*) 'ICTpos=', ICTpos
            WRITE(740+MyRankGlobal,*) 'idx=', idx, ' pos=', pos
            WRITE(740+MyRankGlobal,*) 'U10NEW: min(TEMP)        =', minval(TEMP(1:NP_RES,10))
            WRITE(740+MyRankGlobal,*) 'U10NEW: max(TEMP)        =', maxval(TEMP(1:NP_RES,10))
            WRITE(740+MyRankGlobal,*) 'U10NEW: min(ARR_OUT_SEND)=', minval(ARR_OUT_SEND_C(jProcOut) % ARR(pos,:))
            WRITE(740+MyRankGlobal,*) 'U10NEW: max(ARR_OUT_SEND)=', maxval(ARR_OUT_SEND_C(jProcOut) % ARR(pos,:))
          END IF
          WRITE(740+MyRankGlobal,*) 'Before call to MPI_ISEND'
          WRITE(740+MyRankGlobal,*) 'MPI_EXCH_STR=', MPI_EXCH_STR
          WRITE(740+MyRankGlobal,*) 'Maxlen=', Maxlen
          FLUSH(740+MyRankGlobal)
#endif
          CALL MPI_ISEND(ARR_OUT_SEND_C(jProcOut) % ARR,Maxlen*nbPointCovered, &
&                        MPI_EXCH_STR,eProcOut-1,1020,comm,out_send_rqst(jProcOut),ierr)
#ifdef DEBUG         
          WRITE(740+MyRankGlobal,*) 'After call to MPI_ISEND'
          FLUSH(740+MyRankGlobal)
#endif
        END IF
      END DO
#ifdef DEBUG         
      WRITE(740+MyRankGlobal,*) 'Second part of our operations'
      FLUSH(740+MyRankGlobal)
#endif
      IF (NFLDPPE(IRANK).GT.0) THEN
#ifdef DEBUG
        WRITE(740+MyRankGlobal,*) 'Before MPI_IRECV loop'
        FLUSH(740+MyRankGlobal)
#endif
        iProc=0
        ARR_OUT_RECV_SEP=0
        DO IPR=1,NPROC
          IF (IPR .ne. IRANK) THEN
            iProc=iProc + 1
            CALL MPI_IRECV(ARR_OUT_RECV_SEP,1,block_type(IPR), IPR-1,1020,comm,out_recv_rqst(iProc),ierr)
          END IF
        END DO
#ifdef DEBUG
        WRITE(740+MyRankGlobal,*) 'After MPI_IRECV loop'
        FLUSH(740+MyRankGlobal)
#endif
        iProcOut=ListProcOutRev(IRANK)
        iStart=ListStartIDXout(iProcOut)
        iEnd=ListStartIDXout(iProcOut+1)-1
        eLen=iEnd + 1 - iStart
        DO iNode=1,nbPointCovered
          iNodeGlob=MapFD_locglob(iNode)
          IEfind=ListIEfind(iNode)
          ApplyZMISS(:)=.FALSE.
          DO I=1,3
            IP=INE(I,IEfind)
            IF (IP .le. NP_RES) THEN
              DO idx=iStart,iEnd
                ICTpos=ListICTpos(idx)
                pos=idx + 1 - iStart
                eWI=WIarr(I,iNode)
                IF (eWI .gt. 0) THEN
                  eValIns=TEMP(IP,ICTpos)
#ifdef DEBUG
                  fProcOut = ListICT_to_ProcOut(1)
                  idxHS = ListICTposRev(1)
                  posHS = idxHS + 1 - iStart
                  IF ((MapIPexp_IP(iNodeGlob) .eq. 211).and.(fProcOut .eq. IRANK).and.(posHS.eq.pos)) THEN
                    WRITE(740+MyRankGlobal,*) 'IP=', 211, ' 2: IPexp=', iNodeGlob
                    WRITE(740+MyRankGlobal,*) 'I=', I
                    WRITE(740+MyRankGlobal,*) 'IJ=', IJ
                    WRITE(740+MyRankGlobal,*) 'eValIns=', eValIns, ' IPglob=', IPglob
                    FLUSH(740+MyRankGlobal)
                  END IF
#endif
                  IF (eValIns .eq. ZMISS) THEN
                    ApplyZMISS(pos)=.TRUE.
                  ENDIF
                  eVal=ARR_OUT_RECV_SEP(pos,iNodeGlob)
                  eVal = eVal + eWI*eValIns
                  ARR_OUT_RECV_SEP(pos,iNodeGlob) = eVal
                END IF
              END DO
            END IF
          END DO
          DO pos=1,eLen
            IF (ApplyZMISS(pos)) THEN
              ARR_OUT_RECV_SEP(pos,iNodeGlob) = ZMISS
            END IF
          END DO
        END DO
#ifdef DEBUG         
        WRITE(740+MyRankGlobal,*) 'After ARR_OUT_RECV_SEP assignation'
        FLUSH(740+MyRankGlobal)
#endif
        IF (nproc .gt. 1) THEN
          CALL MPI_WAITALL(nproc-1, out_recv_rqst, out_recv_stat,ierr)
        END IF
#ifdef DEBUG         
        WRITE(740+MyRankGlobal,*) 'After WAITALL for recv_rqst'
        FLUSH(740+MyRankGlobal)
        fProcOut = ListICT_to_ProcOut(10)
        WRITE(740+MyRankGlobal,*) 'fProcOut=', fProcOut
        WRITE(740+MyRankGlobal,*) 'size(ARR_OUT_RECV_SEP,1)=', size(ARR_OUT_RECV_SEP,1)
        IF (fProcOut .eq. IRANK) THEN
          idx = ListICTposRev(10)
          pos = idx + 1 - iStart
          WRITE(740+MyRankGlobal,*) 'idx=', idx, ' pos=', pos
          WRITE(740+MyRankGlobal,*) 'U10NEW: min(ARR_OUT_RECV_SEP)=', minval(ARR_OUT_RECV_SEP(pos,:))
          WRITE(740+MyRankGlobal,*) 'U10NEW: max(ARR_OUT_RECV_SEP)=', maxval(ARR_OUT_RECV_SEP(pos,:))
        END IF
        FLUSH(740+MyRankGlobal)
#endif
      END IF
      IF (NbProcOutRed .gt. 0) THEN
        CALL MPI_WAITALL(NbProcOutRed, out_send_rqst, out_send_stat,ierr)
      END IF
#ifdef DEBUG         
      WRITE(740+MyRankGlobal,*) 'After WAITALL for send_rqst'
      FLUSH(740+MyRankGlobal)
#endif
      IF (NFLDPPE(IRANK).GT.0) THEN
        ARR_OUT_RECV=0
        iProcOut=ListProcOutRev(IRANK)
        iStart=ListStartIDXout(iProcOut)
        iEnd=ListStartIDXout(iProcOut+1)-1
        eLen=iEnd + 1 - iStart
        DO IP=1,NIBLO_FD
          NbEnt=NbEntries(IP)
#ifdef DEBUG
          IF (IP .eq. 211) THEN
            WRITE(740+MyRankGlobal,*) 'IP=', IP, ' NbEnt=', NbEnt
            FLUSH(740+MyRankGlobal)
          END IF
#endif
          IF (NbEnt .eq. 1) THEN
            IPexp = MAP_CovToExp(1,IP)
#ifdef DEBUG
            IF (IP .eq. 211) THEN
              WRITE(740+MyRankGlobal,*) '1: IPexp=', IPexp
              FLUSH(740+MyRankGlobal)
            END IF
#endif
            ARR_OUT_RECV(:,IP) = ARR_OUT_RECV_SEP(:,IPexp)
          ELSE
            ApplyZMISS(:) = .FALSE.
            DO I=1,NbEnt
              IPexp = MAP_CovToExp(I,IP)
#ifdef DEBUG
              IF (IP .eq. 211) THEN
                WRITE(740+MyRankGlobal,*) '2: I=', I, ' IPexp=', IPexp
                FLUSH(740+MyRankGlobal)
              END IF
#endif
              DO iVar=1,eLen
                eValIns=ARR_OUT_RECV_SEP(iVar,IPexp)
                IF (eValIns .eq. ZMISS) THEN
                  ApplyZMISS(iVar)=.TRUE.
                END IF
                ARR_OUT_RECV(iVar,IP) = ARR_OUT_RECV(iVar,IP) + eValIns
              END DO
            END DO
            DO iVar=1,eLen
              IF (ApplyZMISS(iVar)) THEN
                ARR_OUT_RECV(iVar,IP) = ZMISS
              END IF
            END DO
          END IF
        END DO
#ifdef DEBUG
        WRITE(740+MyRankGlobal,*) 'size(ARR_OUT_RECV,1)=', size(ARR_OUT_RECV,1)
        WRITE(740+MyRankGlobal,*) 'size(ARR_OUT_RECV,2)=', size(ARR_OUT_RECV,2)
        eLen=iEnd + 1 - iStart
        DO iLen=1,eLen
          idx=iStart - 1 + iLen
          ICTpos=ListICTpos(idx)
          WRITE(740+MyRankGlobal,*) 'iLen=', iLen, 'ICT=', ICTpos, 'sum=', sum(ARR_OUT_RECV_SEP(iLen,:))
        END DO
        WRITE(740+MyRankGlobal,*) 'ZMISS = ', ZMISS
        FLUSH(740+MyRankGlobal)
        eProcOut=ListICT_to_ProcOut(10)
        WRITE(740+MyRankGlobal,*) 'WIND: eProcOut=', eProcOut
        FLUSH(740+MyRankGlobal)
        IF (eProcOut .eq. IRANK) THEN
          WRITE(740+MyRankGlobal,*) 'After test'
          FLUSH(740+MyRankGlobal)
          pos=LocalPosICT(10)
          WRITE(740+MyRankGlobal,*) 'pos=', pos
          FLUSH(740+MyRankGlobal)
          WRITE(740+MyRankGlobal,*) 'U10NEW: sum(ARR_OUT_RECV(10))=', sum(ARR_OUT_RECV(pos,:))
          WRITE(740+MyRankGlobal,*) 'U10NEW: min(ARR_OUT_RECV(10))=', minval(ARR_OUT_RECV(pos,:))
          WRITE(740+MyRankGlobal,*) 'U10NEW: max(ARR_OUT_RECV(10))=', maxval(ARR_OUT_RECV(pos,:))
          DO iNode=1,NIBLO_FD
            WRITE(740+MyRankGlobal,*) 'U10NEW: iNode=', iNode, ' A_O_R=', ARR_OUT_RECV(pos,iNode)
          END DO
          FLUSH(740+MyRankGlobal)
        END IF
        eProcOut=ListICT_to_ProcOut(1)
        WRITE(740+MyRankGlobal,*) 'HS: eProcOut=', eProcOut
        FLUSH(740+MyRankGlobal)
        IF ((eProcOut .eq. IRANK).and. FirstTime) THEN
          WRITE(740+MyRankGlobal,*) 'HS: After test'
          FLUSH(740+MyRankGlobal)
          pos=LocalPosICT(1)
          nbZmiss=0
          WRITE(740+MyRankGlobal,*) 'pos=', pos
          FLUSH(740+MyRankGlobal)
          WRITE(740+MyRankGlobal,*) 'NIBLO_FD=', NIBLO_FD
          FLUSH(740+MyRankGlobal)
          DO iNode=1,NIBLO_FD
!            WRITE(740+MyRankGlobal,*) 'iNode=', iNode, 'NIBLO_FD=', NIBLO_FD
!            FLUSH(740+MyRankGlobal)
!            WRITE(740+MyRankGlobal,*) 'Before'
!            FLUSH(740+MyRankGlobal)
            IF (ARR_OUT_RECV(pos,iNode) .eq. ZMISS) THEN
              nbZmiss = nbZmiss + 1
!              WRITE(740+MyRankGlobal,*) 'print, step 1'
!              FLUSH(740+MyRankGlobal)
!              WRITE(940+MyRankGlobal,*) 'iNode=', iNode, ' nbZmiss=', nbZmiss
!              WRITE(740+MyRankGlobal,*) 'print, step 2'
!              FLUSH(740+MyRankGlobal)
!              WRITE(740+MyRankGlobal,*) 'iNode=', iNode, ' nbZmiss=', nbZmiss
!              WRITE(740+MyRankGlobal,*) 'print, step 3'
!              FLUSH(740+MyRankGlobal)
            END IF
!            WRITE(740+MyRankGlobal,*) 'After'
!            FLUSH(740+MyRankGlobal)
          END DO
!          WRITE(740+MyRankGlobal,*) 'After loop'
!          FLUSH(740+MyRankGlobal)
!          IF (nbZmiss .gt. 0) THEN
!            FLUSH(940+MyRankGlobal)
!          END IF
          WRITE(740+MyRankGlobal,*) 'size(ARR_OUT_RECV,2)=', size(ARR_OUT_RECV,2)
          IF (size(ARR_OUT_RECV,2) .gt. 220) THEN
            WRITE(740+MyRankGlobal,*) 'ARR_OUT_RECV(pos,211)=', ARR_OUT_RECV(pos,211)
            WRITE(740+MyRankGlobal,*) 'HS: nbZmiss=', nbZmiss
            FLUSH(740+MyRankGlobal)
          END IF
!          IF (NIBLO_FD .gt. 0) THEN
!            STOP
!          END IF
        END IF
        FirstTime = .FALSE.
        WRITE(740+MyRankGlobal,*) 'NIBLO_FD=', NIBLO_FD
        FLUSH(740+MyRankGlobal)
#endif
      END IF
      END SUBROUTINE SET_UP_ARR_OUT_RECV
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INITIAL_OUTPUT_INITS
!     PREPARE OUTPUT OF UNSTRUCTURED GRID TO STRUCTURED GRID (as read in)
      USE YOWPD, ONLY: NE_GLOBAL, NP_GLOBAL, INE_GLOBAL, NODES_GLOBAL
      USE YOWUNPOOL, ONLY : IE_OUTPTS
      USE YOWMAP   , ONLY : AMOWEP, AMONOP, XDELLA, ZDELLO, IPER, NGX, NGY
      IMPLICIT NONE
      INTEGER(KIND=JWIM) :: IX, IY
      REAL(KIND=JWRU) :: XLA, XLO
      IF(ALLOCATED(IE_OUTPTS)) DEALLOCATE(IE_OUTPTS)
      ALLOCATE(IE_OUTPTS(NGX,NGY))
      !AR: add openmp here but better it would be if you would have only one loop over ngx*ngy
      DO IY=1,NGY
        XLA = REAL(AMONOP-REAL(IY-1)*XDELLA,JWRU)
        DO IX=1,NGX
          XLO = REAL(AMOWEP+REAL(IX-1)*ZDELLO(IY),JWRU)
          !!debile
!!! need to find out unstructured grid left longitude (I assume -180 for now)
          IF(XLO.GE.180._JWRU) XLO=max(XLO-IPER*360._JWRU,-180._JWRU)
          CALL FIND_ELE(NE_GLOBAL, NP_GLOBAL, NP_GLOBAL, INE_GLOBAL, NODES_GLOBAL, XLO, XLA, IE_OUTPTS(IX,IY))
        ENDDO
      ENDDO
      END SUBROUTINE INITIAL_OUTPUT_INITS
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INITIAL_OUTPUT_INITS_NEXTGEN
!     PREPARE OUTPUT OF UNSTRUCTURED GRID TO STRUCTURED GRID (as read in)
      USE MPL_MPIF
      USE yowpd, only: MNE=>ne, INE, MNP=>npa, NP_RES=>np
      USE yownodepool, ONLY : iplg, ipgl
      USE YOWMPP   , ONLY : IRANK, NPROC
      USE yowDatapool, only: comm
      USE yownodepool, only : nodes, nodes_global, t_Node
      USE YOW_RANK_GLOLOC, ONLY : MyRankGlobal
      USE YOWUNPOOL, only : NIBLO_FD
      USE YOWMAP   , ONLY : AMOWEP, AMONOP, XDELLA, ZDELLO, IPER, NGX, NGY
      IMPLICIT NONE
      INTEGER(KIND=JWIM) :: IX, IY
      REAL(JWRB) :: XLA, XLO
      INTEGER(KIND=JWIM) :: istatus(MPI_STATUS_SIZE)
      real(KIND=JWRB), allocatable :: XLOmat(:,:), XLAmat(:,:)
      INTEGER(KIND=JWIM), allocatable :: eInt(:)
      INTEGER(KIND=JWIM), allocatable :: ListPointCovered(:)
      INTEGER(KIND=JWIM), allocatable :: IE_OUTPTS_LOC(:,:)
      INTEGER(KIND=JWIM), allocatable :: POSblock(:,:)
      INTEGER(KIND=JWIM), allocatable :: MapXY_iNode(:,:)
      INTEGER(KIND=JWIM), allocatable :: dspl_send(:), LandCovered(:,:)
      real(KIND=JWRU) :: X(3), Y(3)
      real(KIND=JWRU) :: eWI(3)
      INTEGER(KIND=JWIM), allocatable :: ListFirst(:)
      INTEGER(KIND=JWIM), allocatable :: ListNbVarOut(:)
      INTEGER(KIND=JWIM), allocatable :: List_dspl_send(:)
      INTEGER(KIND=JWIM), allocatable :: dspl_sendMPI(:)
      INTEGER(KIND=JWIM), allocatable :: CoveringMap(:,:)
      INTEGER(KIND=JWIM), allocatable :: PointStatus(:)
      INTEGER(KIND=JWIM) :: sumSIZ, ePos, eProcOut
      INTEGER(KIND=JWIM) :: ICT, idx1, idx2, idx, IPglob, IPR, len
      INTEGER(KIND=JWIM) :: idxPosBlock, I, iProc
      INTEGER(KIND=JWIM) :: ierr, IP, IEfind
      INTEGER(KIND=JWIM) :: iProcOut, istat
      INTEGER(KIND=JWIM) :: nbPointCoveredLoc
      INTEGER(KIND=JWIM) :: TotalOutVar
      INTEGER(KIND=JWIM) :: iNode, IXY
      INTEGER(KIND=JWIM) :: TheIE, iNodeExp
      INTEGER(KIND=JWIM) :: Ifound, eFirst
      INTEGER(KIND=JWIM) :: I1, I2, iProc1, iProc2
      INTEGER(KIND=JWIM) :: eEnt
      INTEGER(KIND=JWIM) :: iVarOut
      logical :: HasOutputNode
      TYPE(t_Node), allocatable :: nodes_loc(:)
      REAL(KIND=JWRU) :: eX, eY
      IF (JWRB .eq. 4) THEN
        MPI_EXCH_STR = MPI_REAL4
      ELSE
        MPI_EXCH_STR = MPI_REAL8
      END IF
      !
      ! Creating the arrays to determine which nodes belong to the list of points
      !
      ALLOCATE(IE_OUTPTS_LOC(NGX,NGY))
      !AR: add openmp here but better it would be if you would have only one loop over ngx*ngy
      ALLOCATE(XLOmat(NGX,NGY), XLAmat(NGX,NGY), stat=istat)
      nbPointCovered=0
      allocate(nodes_loc(MNP), stat=istat)
      DO IP=1,MNP
        nodes_loc(IP) = nodes_global(iplg(IP))
      END DO
      DO IY=1,NGY
        XLA = AMONOP-REAL(IY-1)*XDELLA
        DO IX=1,NGX
          XLO = AMOWEP+REAL(IX-1)*ZDELLO(IY)
          !!debile
!!! need to find out unstructured grid left longitude (I assume -180 for now)
          IF(XLO.GE.180.d0) XLO=max(XLO-IPER*360.d0,-180.d0)
          XLOmat(IX,IY)=XLO
          XLAmat(IX,IY)=XLA
          CALL FIND_ELE(MNE, MNP, NP_RES, INE, nodes_loc, DBLE(XLO), DBLE(XLA), TheIE)
          IE_OUTPTS_LOC(IX,IY) = TheIE
          IF (TheIE .ne. -1) THEN
            nbPointCovered = nbPointCovered + 1
            HasOutputNode=.FALSE.
            DO I=1,3
              IP=INE(I,TheIE)
              IF (IP .le. NP_RES) THEN
                HasOutputNode=.TRUE.
              END IF
            END DO
#ifdef DEBUG
            IF (HasOutputNode .eqv. .FALSE.) THEN
              Print *, 'One assumption is proven false'
              Print *, 'There are triangles with all nodes being ghost !'
              STOP
            END IF
#endif
          END IF
        ENDDO
      ENDDO
      deallocate(nodes_loc)
#ifdef DEBUG
      WRITE(740+MyRankGlobal,*) 'XLO : min=', minval(XLOmat), ' max=', maxval(XLOmat)
      WRITE(740+MyRankGlobal,*) 'XLA : min=', minval(XLAmat), ' max=', maxval(XLAmat)
      FLUSH(740+MyRankGlobal)
#endif
      !
      ! We need to manage the creation of ZSENDBUF arrays for data sending
      ! Goal is to have
      !
      allocate(WIarr(3,nbPointCovered), ListIEfind(nbPointCovered), POSblock(NGX,NGY))
      idx=0
      idxPosBlock=0
      DO IY=1,NGY
        DO IX=1,NGX
          IEfind=IE_OUTPTS_LOC(IX,IY)
          IF (IEfind .ne. -1) THEN
            idx=idx+1
            ListIEfind(idx)=IEfind
            XLO=XLOmat(IX,IY)
            XLA=XLAmat(IX,IY)
            DO I=1,3
              IP=INE(I,IEfind)
              X(I)=REAL(nodes_global(iplg(IP)) % X, JWRU)
              Y(I)=REAL(nodes_global(iplg(IP)) % Y, JWRU)
            END DO
            eX=REAL(XLO, JWRU)
            eY=REAL(XLA, JWRU)
            CALL INTELEMENT_COEFF(X, Y, eX, eY, eWI)
            WIarr(:,idx) = REAL(eWI, JWRB)
          END IF
          POSblock(IX,IY)=NGX*(IY-1) + IX
        END DO
      END DO
      idx=0
      allocate(dspl_send(nbPointCovered))
      DO IY=1,NGY
        DO IX=1,NGX
          IF (IE_OUTPTS_LOC(IX,IY) .ne. -1) THEN
            idx=idx+1
            idxPosBlock=POSblock(IX,IY)
            dspl_send(idx)=idxPosBlock
          END IF
        END DO
      END DO
      allocate(ListPointCovered(nproc))
      allocate(eInt(1))
      IF (IRANK .eq. 1) THEN
        ListPointCovered(1)=nbPointCovered
        DO iProc=2,nproc
          CALL MPI_RECV(eInt,1,MPI_INTEGER,iProc-1, 53, comm, istatus, ierr)
          ListPointCovered(iProc) = eInt(1)
        END DO
        DO iProc=2,nproc
          CALL MPI_SEND(ListPointCovered,nproc,MPI_INTEGER, iProc-1, 54, comm, ierr)
        END DO
      ELSE
        eInt(1)=nbPointCovered
        CALL MPI_SEND(eInt,1,MPI_INTEGER, 0, 53, comm, ierr)
        CALL MPI_RECV(ListPointCovered,nproc,MPI_INTEGER,0, 54, comm, istatus, ierr)
      END IF
      deallocate(eInt)
#ifdef DEBUG
      IF (minval(ListPointCovered) .eq. 0) THEN
        WRITE(740+MyRankGlobal,*) 'We have ListPointCovered = 0 proven wrong'
        WRITE(740+MyRankGlobal,*) 'Be mindful'
        FLUSH(740+MyRankGlobal)
      END IF
#endif
      !
      ! ListFirst
      !
      allocate(ListFirst(nproc))
      ListFirst=0
      DO iProc=2,nproc
        ListFirst(iProc)=ListFirst(iProc-1) + ListPointCovered(iProc-1)
      END DO
      !
      ! Accumulating all the dspl_send arrays
      !
      sumSIZ = sum(ListPointCovered)
      NIBLO_FD_EXP = sumSIZ
      allocate(List_dspl_send(sumSIZ))
      IF (IRANK .eq. 1) THEN
        idx=0
        DO IP=1,nbPointCovered
          idx=idx+1
          List_dspl_send(idx)=dspl_send(IP)
        END DO
        DO iProc=2,nproc
          nbPointCoveredLoc=ListPointCovered(iProc)
          allocate(eInt(nbPointCoveredLoc))
          CALL MPI_RECV(eInt,nbPointCoveredLoc,MPI_INTEGER,iProc-1, 55, comm, istatus, ierr)
          DO IP=1,nbPointCoveredLoc
            idx=idx+1
            List_dspl_send(idx)=eInt(IP)
          END DO
          deallocate(eInt)
        END DO
        DO iProc=2,nproc
          CALL MPI_SEND(List_dspl_send,sumSIZ,MPI_INTEGER, iProc-1, 56, comm, ierr)
        END DO
      ELSE
        CALL MPI_SEND(dspl_send,nbPointCovered,MPI_INTEGER, 0, 55, comm, ierr)
        CALL MPI_RECV(List_dspl_send,sumSIZ,MPI_INTEGER,0, 56, comm, istatus, ierr)
      END IF
      !
      ! Determining the land covered by the points
      !
      allocate(LandCovered(NGX,NGY), stat=istat)
      LandCovered=0
      DO idx=1,sumSIZ
        IXY=List_dspl_send(idx)
        IX = 1 + MOD(IXY-1,NGX)
        IY = 1 + (IXY - IX)/NGX
#ifdef DEBUG          
        IF (POSblock(IX,IY) .ne. IXY) THEN
          Print *, 'Coherency error'
          Print *, 'idx=', idx
          Print *, 'IX / IY / IXY / POS=', IX, IY, IXY, POSblock(IX,IY)
          STOP
        END IF
#endif
        LandCovered(IX,IY)=1
      END DO
      !
      ! Determination of total number of sea points
      !
      NIBLO_FD=0
      DO IX=1,NGX
        DO IY=1,NGY
          IF (LandCovered(IX,IY) .eq. 1) THEN
            NIBLO_FD = NIBLO_FD + 1
          END IF
        END DO
      END DO
#ifdef DEBUG          
      WRITE(740+MyRankGlobal,*) 'NIBLO_FD=', NIBLO_FD
      FLUSH(740+MyRankGlobal)
#endif
      !
      ! construction of the mapping array global
      !
      allocate(IXarr(NIBLO_FD), IYarr(NIBLO_FD), stat=istat)
      allocate(MapXY_iNode(NGX,NGY), stat=istat)
!      WRITE(740+MyRankGlobal,*) 'allocated(MapXY_iNode)=', allocated(MapXY_iNode)
!      WRITE(740+MyRankGlobal,*) 'NGX=', NGX
!      WRITE(740+MyRankGlobal,*) 'NGY=', NGY
!      WRITE(740+MyRankGlobal,*) 'size(MapXY_iNode)=', size(MapXY_iNode)
      MapXY_iNode=-1
      idx=0
      DO IX=1,NGX
        DO IY=1,NGY
          IF (LandCovered(IX,IY) .eq. 1) THEN
            idx=idx+1
            MapXY_iNode(IX,IY)=idx
            IXarr(idx)=IX
            IYarr(idx)=IY
          END IF
        END DO
      END DO
      !
      ! allocation of array output
      !
      allocate(ListNbVarOut(nproc), stat=istat)
      ListNbVarOut = 0
      DO ICT=1,JPPFLAG
        IPR=IPFGTBL(ICT)
#ifdef DEBUG
        WRITE(740+MyRankGlobal,*) 'ICT=', ICT, ' IPR=', IPR
        FLUSH(740+MyRankGlobal)
#endif
        IF (IPR .ne. 0) ListNbVarOut(IPR) = ListNbVarOut(IPR) + 1
      END DO
      MaxLen=maxval(ListNbVarOut)
      allocate(ApplyZMISS(MaxLen))
      TotalOutVar=sum(ListNbVarOut)
      NbProcOut=0
      DO IPR=1,nproc
#ifdef DEBUG
        WRITE(740+MyRankGlobal,*) 'IPR=', IPR, ' NbVarOut=', ListNbVarOut(IPR)
        FLUSH(740+MyRankGlobal)
#endif
        IF (ListNbVarOut(IPR) .gt. 0) NbProcOut = NbProcOut + 1 
      END DO
#ifdef DEBUG
      WRITE(740+MyRankGlobal,*) 'NbProcOut = ', NbProcOut
      FLUSH(740+MyRankGlobal)
#endif
      allocate(ListProcOut(NbProcOut), ListProcOutRev(nproc), stat=istat)
      ListProcOutRev=-1
      idx=0
      DO IPR=1,nproc
        IF (ListNbVarOut(IPR) .gt. 0) THEN
          idx=idx+1
          ListProcOut(idx) = IPR
          ListProcOutRev(IPR) = idx
        END IF
      END DO
#ifdef DEBUG
      DO iProcOut=1,NbProcOut
        eProcOut=ListProcOut(iProcOut)
        WRITE(740+MyRankGlobal,*) 'iProcOut=', iProcOut, ' eProcOut=', eProcOut
        FLUSH(740+MyRankGlobal)
      END DO
      DO IPR=1,nproc
        WRITE(740+MyRankGlobal,*) 'IPR=', IPR, ' ListProcOutRev=', ListProcOutRev(IPR)
        FLUSH(740+MyRankGlobal)
      END DO
#endif
      allocate(ListStartIDXout(NbProcOut+1), ListICTpos(TotalOutVar), ListICTposRev(JPPFLAG), stat=istat)
      allocate(LocalPosICT(JPPFLAG), ListICT_to_ProcOut(JPPFLAG), stat=istat)
      ListICTposRev = -1
      LocalPosICT=0
      ListStartIDXout(1)=1
      ListICT_to_ProcOut = -1
      idx1=0
      DO iProcOut=1,NbProcOut
        eProcOut=ListProcOut(iProcOut)
        len=ListNbVarOut(eProcOut)
        ListStartIDXout(iProcOut+1) = ListStartIDXout(iProcOut) + len
        idx2=0
        DO ICT=1,JPPFLAG
          IF (IPFGTBL(ICT) .eq. eProcOut) THEN
            idx2=idx2 + 1
            idx1=idx1 + 1
            ListICTpos(idx1) = ICT
            ListICTposRev(ICT) = idx1
            LocalPosICT(ICT) = idx2
            ListICT_to_ProcOut(ICT) = eProcOut
          END IF
        END DO
#ifdef DEBUG
        IF (idx2 .ne. len) THEN
          Print *, 'inconsistency between idx2 and len'
          STOP
        END IF
#endif
      END DO
#ifdef DEBUG
      DO iVarOut=1,TotalOutVar
        WRITE(740+MyRankGlobal,*) 'iVarOut=', iVarOut, ' idx=', ListICTpos(iVarOut)
        FLUSH(740+MyRankGlobal)
      END DO
      DO iProcOut=1,NbProcOut+1
        WRITE(740+MyRankGlobal,*) 'iProcOut=', iProcOut, ' start=', ListStartIDXout(iProcOut)
        FLUSH(740+MyRankGlobal)
      END DO
      DO ICT=1,JPPFLAG
        WRITE(740+MyRankGlobal,*) 'ICT=', ICT, ' eProcOut=', ListICT_to_ProcOut(ICT)
        FLUSH(740+MyRankGlobal)
      END DO
      IF (idx1 .ne. TotalOutVar) THEN
        Print *, 'inconsistency between idx1 and TotalOutVar'
        STOP
      END IF
#endif
      !
      ! determination of land covering
      !
      allocate(CoveringMap(3,NIBLO_FD), stat=istat)
      CoveringMap=0
      DO iProc=1,nproc
        eFirst=ListFirst(iProc)
        len=ListPointCovered(iProc)
        DO idx=1,len
          IXY=List_dspl_send(idx + eFirst)
          IX = 1 + MOD(IXY-1,NGX)
          IY = 1 + (IXY - IX)/NGX
          iNode = MapXY_iNode(IX,IY)
          ePos=1
          DO I=1,3
            IF (CoveringMap(I,iNode) .ne. 0) THEN
              ePos=I+1
            END IF
          END DO
#ifdef DEBUG          
          IF (ePos .eq. 4) THEN
            Print *, 'ePos is larger than we expected. Possible bug'
            STOP
          END IF
#endif       
          CoveringMap(ePos,iNode)=iProc
        END DO
      END DO
      DO iNode=1,NIBLO_FD
        len=0
        DO I=1,3
          IF (CoveringMap(I,iNode) .ne. 0) THEN
            len=I
          END IF
        END DO
#ifdef DEBUG          
        IF (len .eq. 0) THEN
          Print *, 'We should have len=0'
          STOP
        END IF
#endif       
        DO I1=1,len-1
           DO I2=I1+1,len
            iProc1=CoveringMap(I1,iNode)
            iProc2=CoveringMap(I2,iNode)
            IF (iProc2 .lt. iProc1) THEN
              CoveringMap(I1,iNode)=iProc2
              CoveringMap(I2,iNode)=iProc1
            END IF
          END DO
        END DO
      END DO
      !
      ! Computation of reverse maps
      !
      allocate(MAP_CovToExp(3,NIBLO_FD), MapIPexp_IP(NIBLO_FD_EXP), stat=istat)
      MAP_CovToExp=0
      idx=0
      DO IP=1,NIBLO_FD
        DO I=1,3
          IF (CoveringMap(I,IP) .ne. 0) THEN
            idx=idx+1
            MAP_CovToExp(I,IP) = idx
            MapIPexp_IP(idx) = IP
          END IF
        END DO
      END DO
      !
      ! Computation of number of entries per output point
      !
      allocate(NbEntries(NIBLO_FD), stat=istat)
      DO IP=1,NIBLO_FD
        eEnt=0
        DO I=1,3
          IF (CoveringMap(I,IP) .ne. 0) THEN
            eEnt=I
          END IF
        END DO
        NbEntries(IP) = eEnt
      END DO
      !
      ! allocation of the MPI stuff
      !
      allocate(out_recv_rqst(nproc-1), out_recv_stat(MPI_STATUS_SIZE,nproc-1), stat=istat)
      NbProcOutRed=0
      DO iProcOut=1,NbProcOut
        eProcOut = ListProcOut(iProcOut)
        IF (eProcOut .ne. IRANK) THEN
          NbProcOutRed = NbProcOutRed + 1
        END IF
      END DO
      allocate(out_send_rqst(NbProcOutRed), out_send_stat(MPI_STATUS_SIZE,NbProcOutRed), stat=istat)
      !
      ! Construction of single mapping local global for FD
      !
      allocate(MapFD_locglob(nbPointCovered))
      DO IP=1,nbPointCovered
        IXY=dspl_send(IP)
        IX = 1 + MOD(IXY-1,NGX)
        IY = 1 + (IXY - IX)/NGX
        iNode = MapXY_iNode(IX,IY)
        Ifound=-1
        DO I=1,3
          IF (CoveringMap(I,iNode) .eq. IRANK) THEN
            Ifound=I
          END IF
        END DO
#ifdef DEBUG          
        IF (Ifound .eq. -1) THEN
          Print *, 'Error, we should have iFound non zero'
          STOP
        END IF
#endif
        iNodeExp=MAP_CovToExp(Ifound, iNode)
        MapFD_locglob(IP) = iNodeExp
      END DO
      !
      ! Construction of the MPI arrays
      !
      allocate(block_type(nproc), stat=istat)
#ifdef DEBUG
      allocate(PointStatus(sumSIZ), stat=istat)
      PointStatus=0
#endif      
      DO iProc=1,nproc
        nbPointCoveredLoc=ListPointCovered(iProc)
        allocate(dspl_sendMPI(nbPointCoveredLoc))
        DO IP=1,nbPointCoveredLoc
          IXY=List_dspl_send(IP + ListFirst(iProc))
          IX = 1 + MOD(IXY-1,NGX)
          IY = 1 + (IXY - IX)/NGX
          iNode = MapXY_iNode(IX,IY)
          Ifound=-1
          DO I=1,3
            IF (CoveringMap(I,iNode) .eq. iProc) THEN
              Ifound=I
            END IF
          END DO
#ifdef DEBUG          
          IF (Ifound .eq. -1) THEN
            Print *, 'Error, we should have iFound non zero'
            STOP
          END IF
#endif
          iNodeExp=MAP_CovToExp(Ifound, iNode)
#ifdef DEBUG
          IF (PointStatus(iNodeExp) .ne. 0) THEN
            Print *, 'Node has already been assigned'
            Print *, ' val=', PointStatus(iNodeExp)
            Print *, 'iNodeExp=', iNodeExp
            Print *, 'sumSIZ=', sumSIZ
            STOP
          END IF
          PointStatus(iNodeExp) = iProc
          IF ((iNodeExp .gt. sumSIZ).or.(iNodeExp .lt. 1)) THEN
            Print *, 'iNodeExp is too large or too small'
            Print *, 'iNodeExp=', iNodeExp
            Print *, 'sumSIZ=', sumSIZ
            STOP
          END IF
          IF (iNodeExp .eq. 1) THEN
            WRITE(740+MyRankGlobal,*) 'iProc=', iProc, ' dspIP=', IP, ' iNodeExp=', iNodeExp
            FLUSH(740+MyRankGlobal)
          END IF
#endif
          dspl_sendMPI(IP)=Maxlen * (iNodeExp - 1)
        END DO
#ifdef DEBUG
        WRITE(740+MyRankGlobal,*) 'MPI_EXCH_STR=', MPI_EXCH_STR
        FLUSH(740+MyRankGlobal)
#endif
        call mpi_type_create_indexed_block(nbPointCoveredLoc,Maxlen,dspl_sendMPI,MPI_EXCH_STR,block_type(iProc), ierr)
        call mpi_type_commit(block_type(iProc), ierr)
        deallocate(dspl_sendMPI)
      END DO
#ifdef DEBUG
      DO IP=1,sumSIZ
        IF (PointStatus(IP) .eq. 0) THEN
          Print *, 'point IP has not been reached'
          STOP
        END IF
      END DO
      deallocate(PointStatus)
#endif
      !
      ! deallocations
      !
      deallocate(List_dspl_send)
      deallocate(dspl_send)
      deallocate(IE_OUTPTS_LOC, XLOmat, XLAmat)
      deallocate(POSblock)
      deallocate(LandCovered)
      deallocate(ListPointCovered)
      deallocate(MapXY_iNode)
      deallocate(ListFirst)
      !
      ! allocation
      !
#ifdef DEBUG
      WRITE(740+MyRankGlobal,*) 'Before alloc ARR_OUT_SEND'
      FLUSH(740+MyRankGlobal)
#endif
      allocate(ARR_OUT_SEND_C(NbProcOutRed))
      DO iProc=1,NbProcOutRed
        allocate(ARR_OUT_SEND_C(iProc) % ARR(Maxlen, nbPointCovered))
      END DO
      IF (ListNbVarOut(IRANK) .gt. 0) THEN
        allocate(ARR_OUT_RECV_SEP(MaxLen, sumSIZ))
        allocate(ARR_OUT_RECV(MaxLen, NIBLO_FD))
      END IF
#ifdef DEBUG
      WRITE(740+MyRankGlobal,*) 'After alloc ARR_OUT_RECV'
      FLUSH(740+MyRankGlobal)
#endif
      END SUBROUTINE INITIAL_OUTPUT_INITS_NEXTGEN
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE FIND_ELE(MNEloc, MNPloc, NP_RESloc, INEloc, NODESloc, Xo,Yo,IE )

!     ! Xo, Yo is the point of interest
!     ! IE is the global element id, which contains the point Xo, Yo
!
!      use yowpd, only: MNE=>ne_global, INE=>INE_GLOBAL, NODES=>nodes_global
!
      USE yowNodepool, only : t_Node
      IMPLICIT NONE

! Dummy arguments
!
      INTEGER(KIND=JWIM), intent(in) :: MNEloc, MNPloc, NP_RESloc
      INTEGER(KIND=JWIM), intent(in) :: INEloc(3,MNEloc)
      TYPE(t_Node) :: NODESloc(MNPloc)
      REAL(KIND=JWRU), INTENT(IN) :: Xo , Yo
      INTEGER(KIND=JWIM), INTENT(INOUT) :: IE 
!
! Local variables
!
      REAL(KIND=JWRU), SAVE   :: xi, xj, xk, yi, yj, yk, dx, dy, f, Xo8, Yo8
      REAL(KIND=JWRU), SAVE   :: xmax , xmin , ymax , ymin
      INTEGER(KIND=JWIM), SAVE  :: if0 , if1 , ijk , k , ki , kj , kk , l
      INTEGER(KIND=JWIM), SAVE  :: i , i0 , idx , ielem
      

      REAL(KIND=JWRU):: THR = TINY(1.0_JWRU)

      DATA idx/0/ , i0/0/ , ielem/1/ , if0/0/ , if1/0/
!
!     Laengste Kannte (DX,DY) bestimmen
!     Dieser Programmabschnitt wird nur beim ersten Aufruf
!     durchlaufen!
!
      Xo8 = REAL(Xo,JWRU)
      Yo8 = REAL(Yo,JWRU)

      IF ( idx/=1 ) THEN
        idx = 1
        DO i = 1 , MNEloc
          ki = INEloc(1,i) ! + 1
          kj = INEloc(2,i) ! + 1
          kk = INEloc(3,i) ! + 1
          xi = REAL(NODESloc(ki)%X,JWRU)
          yi = REAL(NODESloc(ki)%Y,JWRU)
          xj = REAL(NODESloc(kj)%X,JWRU)
          yj = REAL(NODESloc(kj)%Y,JWRU)
          xk = REAL(NODESloc(kk)%X,JWRU)
          yk = REAL(NODESloc(kk)%Y,JWRU)
          IF ( i==1 ) THEN
            dx = MAX(ABS(xi-xj),ABS(xi-xk),ABS(xj-xk))
            dy = MAX(ABS(yi-yj),ABS(yi-yk),ABS(yj-yk))
          ELSE
            dx = MAX(dx,ABS(xi-xj),ABS(xi-xk),ABS(xj-xk))
            dy = MAX(dy,ABS(yi-yj),ABS(yi-yk),ABS(yj-yk))
          END IF
        END DO
      END IF
!     ------------------------------------------------------------------
!     TEST, OB DER PUNKT IM ZULETZT ANGESPROCHENEN ELEMENT LIEGT
!     ------------------------------------------------------------------
      IF ( i0==1 .AND. IE/=-1 ) THEN
         IF ( Yo8-ymin > THR ) THEN
            IF ( Yo8-ymax < -THR ) THEN
               IF ( Xo8-xmin > THR ) THEN
                  IF ( Xo8-xmax < -THR ) THEN
                     f = xi*(yj-Yo8) + xj*(Yo8-yi) + Xo8*(yi-yj)
                     IF ( f > THR ) THEN
                        f = xj*(yk-Yo8) + xk*(Yo8-yj) + Xo8*(yj-yk)
                        IF ( f > THR  ) THEN
                           f = xk*(yi-Yo8) + xi*(Yo8-yk) + Xo8*(yk-yi)
                           IF ( f > THR ) THEN
                              IE = ielem ! Element gefunden -->RETURN
                              RETURN
                           endif
                        endif
                     endif
                  endif
               endif
            endif
         endif
      endif
!     ------------------------------------------------------------------
!     Element suchen
!     ------------------------------------------------------------------
      i0 = 1
      i = ielem
      IF ( i<1 ) i = 1
      k = i
      l = i
      ijk = 0

100   DO
         ijk = ijk + 1
!.....   ABFRAGE AUF X-Richtung
         ki = INEloc(1,i)! + 1
         xi = REAL(NODESloc(ki)%X,JWRU)
         IF ( DABS(xi-Xo8)<=dx ) THEN
            kj = INEloc(2,i)! + 1
            kk = INEloc(3,i)! + 1
            xj = REAL(NODESloc(kj)%X,JWRU)
            xk = REAL(NODESloc(kk)%X,JWRU)
!.....    Punkt ausserhalb Element:
            xmin = MIN(xi,xj,xk)
            IF ( Xo8>=xmin ) THEN
               xmax = MAX(xi,xj,xk)
               IF ( Xo8<=xmax ) THEN
!.....        ABFRAGE AUF Y-Richtung
                  yi = REAL(NODESloc(ki)%Y,JWRU)
                  IF ( DABS(yi-Yo8)<=dy ) THEN
                     yj = REAL(NODESloc(kj)%Y,JWRU)
                     yk = REAL(NODESloc(kk)%Y,JWRU)
!.....          Punkt ausserhalb Element:
                     ymin = MIN(yi,yj,yk)
                     IF ( Yo8>=ymin ) THEN
                        ymax = MAX(yi,yj,yk)
                        IF ( Yo8<=ymax ) THEN
!.....              Bis jetzt liegt Punkt innerhalb des das Element
!                   umschlieszenden Rechtecks XMIN/XMAX, YMIN/YMAX
!                   Pruefen, ob Punkt wirklich innerhalb DREIECK-Element
!                   liegt: BERECHNUNG DER TEILFLAECHEN (ohne 0.5)
                           f = xi*(yj-Yo8) + xj*(Yo8-yi) + Xo8*(yi-yj)
                           IF ( f>=0.0_JWRU) THEN
                              f = xj*(yk-Yo8) + xk*(Yo8-yj)+Xo8*(yj-yk)
                              IF ( f>=0.0_JWRU ) THEN
                                 f = xk*(yi-Yo8)+xi*(Yo8-yk)+Xo8*(yk-yi)
                                 IF ( f>=0.0_JWRU ) THEN
                                    IE = i
                                    ielem = IE 
                                    RETURN
                                 endif
                              endif
                           endif
                        endif
                     endif
                  endif
               endif
            endif
         endif
!     SCHLEIFE UEBER ALLE ELEMENTE wird hier folgendermassen hochgezaehlt:
!     beginnend bei IEALT, im Wechsel nach vorn und rueckwaerts suchend
         IF ( k<MNEloc .AND. if1==0 ) THEN
            if0 = 0
            IF ( l>1 ) if1 = 1
            k = k + 1
            i = k
            IF ( ijk<=MNEloc ) CYCLE
         endif
         CONTINUE
         EXIT
      ENDDO

      IF ( l>1 .AND. if0==0 ) THEN
         if1 = 0
         IF ( k<MNEloc ) if0 = 1
         l = l - 1
         i = l
      endif

      IF ( ijk<=MNEloc ) GOTO 100

      IE = -1

      END SUBROUTINE FIND_ELE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INTELEMENT_IPOL(LLCLST,XYELE,SKALAR,Xo,Yo,IE,PMS8,VAL)

      USE YOWUNPOOL, ONLY : DEGRAD

      IMPLICIT NONE

!     Wrapper for linter interpolation ... 
!     however if LLCLST is true then the value of the closest point is taken.
!     XYELE contain the XY coordinates for the 3 points of the triangle IE
!     Xo, Yo is the point of interest
!     IE is the element number the contains this point
!     SKALAR are the 3 values of interest stored at the nodes of each triangle 
!     VAL is the answer of the interpolation 
!
! Dummy arguments
!
      INTEGER(KIND=JWIM), INTENT(IN)  :: IE

      REAL(KIND=JWRU),  INTENT(IN)  :: XYELE(2,3), SKALAR(3), Xo, Yo, PMS8 
      REAL(KIND=JWRU),  INTENT(OUT) :: VAL

      LOGICAL, INTENT(IN)  :: LLCLST
!
! Local variables

      INTEGER(KIND=JWIM) :: i
      REAL(KIND=JWRU) :: distmin, distmax 
      REAL(KIND=JWRU), DIMENSION(3) :: x, y, dist

      do i=1,3
        x(i) = XYELE(1,i)
        y(i) = XYELE(2,i)
      enddo

      IF (LLCLST) THEN
        do i=1,3
          call SPHERICAL_COORDINATE_DISTANCE(Xo, x(i), Yo, y(i), dist(i))
        enddo
        distmax=maxval(dist)
        do i=1,3
          if(SKALAR(i) == PMS8) dist(i)=distmax
        enddo
        distmin = minval(dist)
        val=SKALAR(1)
        if(dist(2) == distmin) then
          val=SKALAR(2)
        else if(dist(3) == distmin) then
          val=SKALAR(3)
        endif
      ELSE
        call linearInterpolationTriangle(x(1), y(1), SKALAR(1),        &
     &                                   x(2), y(2), SKALAR(2),        &
     &                                   x(3), y(3), SKALAR(3),        &
     &                                   Xo, Yo, val)
      ENDIF

      END SUBROUTINE INTELEMENT_IPOL
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE IPELEMENT_CLOSEST(Xo,Yo,IP,IE,Xn,Yn,distmin)

      USE YOWUNPOOL, ONLY : DEGRAD
      USE YOWPD,     ONLY : NE_GLOBAL, NP_GLOBAL, INE_GLOBAL, NODES_GLOBAL
      USE YOWPD,     ONLY : NODES=>nodes_global
      USE YOWMAP   , ONLY : IPER


      IMPLICIT NONE

!    Finds the GLOBAL node index that is closest to a point(Xo,Yo) contained in an element
!     Xo, Yo is the point of interest
!     IE is the global element number the contains this point
!     IP is the global node index
!     distmin is the distance on a sphere of radius 1 between that point and node IP
!
! Dummy arguments
!
      INTEGER(KIND=JWIM), INTENT(OUT)  :: IP, IE

      REAL(KIND=JWRU),  INTENT(IN)  :: Xo, Yo
      REAL(KIND=JWRU),  INTENT(OUT)  :: Xn, Yn, distmin

!
! Local variables

      INTEGER(KIND=JWIM) :: i, ic
      INTEGER(KIND=JWIM) :: NI(3)
      REAL(KIND=JWRU), DIMENSION(3) :: x, y, dist


!     FIND IF  POINT (Xo, Yo) is INSIDE AN ELEMENT
      CALL FIND_ELE(NE_GLOBAL, NP_GLOBAL, NP_GLOBAL, INE_GLOBAL, NODES_GLOBAL, Xo, Yo, IE)

!     IF INSIDE AN ELEMENT, FIND THE CLOSED NODE
      IF(IE.GT.-1) THEN
        NI = INE_GLOBAL(:,IE)
        X(:)=NODES(NI)%X
        Y(:)=NODES(NI)%Y

        do i=1,3
          CALL SPHERICAL_COORDINATE_DISTANCE(Xo, x(i), Yo, y(i), dist(i))
        enddo

        distmin = minval(dist)
        IP=NI(1)
        ic=1
        if(dist(2) == distmin) then
          IP=NI(2)
          ic=2
        else if(dist(3) == distmin) then
          IP=NI(3)
          ic=3
        endif

        Xn=X(ic)
        Yn=Y(ic)

      ELSE
        IP=0
        Xn=Xo
        Yn=Yo
        distmin=-1._JWRU
      ENDIF

      END SUBROUTINE IPELEMENT_CLOSEST
!**********************************************************************
!*                                                                    *
!**********************************************************************
      subroutine linearInterpolationTriangle(x1, y1, z1, x2, y2,       &
     &                             z2, x3, y3, z3, Xo, Yo, zout)
      !> calc a linear interpolation of a triangle.
      !> describes on this side:
      !> http://www.ems-i.com/smshelp/Data_Module/Interpolation/Linear_Interpolationsms.htm
      !> @param[in] x1, y1, z1 The first Point of the Triangle
      !> @param[in] x2, y2, z2 The second Point of the Triangle
      !> @param[in] x3, y3, z3 The third Point of the Triangle
      !> @param[in] Xo, Yo Point to interpolate
      !> @param[out] zout Interpolated Z-value
        implicit none 
        real(KIND=JWRU), intent(in) :: x1,y1,z1,x2,y2,z2,x3,y3,z3,Xo,Yo
        real(KIND=JWRU), intent(out) :: zout

        real(KIND=JWRU) :: A, B, C, D

        A = y1*(z2 - z3) + y2*(z3 - z1) + y3*(z1 - z2)
        B = z1*(x2 - x3) + z2*(x3 - x1) + z3*(x1 - x2)
        C = x1*(y2 - y3) + x2*(y3 - y1) + x3*(y1 - y2)
        D = -A*x1 - B*y1 - C*z1

        zout = - A/C*Xo - B/C*Yo - D/C
      end subroutine linearInterpolationTriangle
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INTELEMENT_COEFF(X, Y, Xp, Yp, WI)
#ifdef DEBUG
      USE YOW_RANK_GLOLOC, ONLY : MyRankGlobal
#endif
      IMPLICIT NONE
      REAL(KIND=JWRU),    INTENT(IN)  :: X(3), Y(3)
      REAL(KIND=JWRU),    INTENT(IN)  :: Xp, Yp
      REAL(KIND=JWRU),  INTENT(OUT) :: WI(3)
      REAL(KIND=JWRU) :: y1,y2,y3,x1,x2,x3
      REAL(KIND=JWRU) :: n1, d1, n2, d2, n3, d3
#ifdef DEBUG
      REAL(KIND=JWRU) :: ErrX, ErrY
#endif
      x1 = X(1)
      x2 = X(2)
      x3 = X(3)
      y1 = Y(1)
      y2 = Y(2)
      y3 = Y(3)
      n1=(Xp-x2)*(y3-y2) - (Yp-y2)*(x3-x2)
      d1=(x1-x2)*(y3-y2) - (y1-y2)*(x3-x2)
      n2=(Xp-x1)*(y3-y1) - (Yp-y1)*(x3-x1)
      d2=(x2-x1)*(y3-y1) - (y2-y1)*(x3-x1)
      n3=(Xp-x1)*(y2-y1) - (Yp-y1)*(x2-x1)
      d3=(x3-x1)*(y2-y1) - (y3-y1)*(x2-x1)
      Wi(1)=n1/d1
      Wi(2)=n2/d2
      Wi(3)=n3/d3
#ifdef DEBUG
      ErrX=Xp - x1*Wi(1) - x2*Wi(2) - x3*Wi(3)
      ErrY=Yp - y1*Wi(1) - y2*Wi(2) - y3*Wi(3)
!      WRITE(740+MyRankGlobal,*) 'ErrX/Y=', ErrX, ErrY
#endif
      END SUBROUTINE INTELEMENT_COEFF
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INTELEMENT(X,Y,Z,XP,YP,Wi,Zi,LSAME)

!     Purpose.: Advects the spectra using non-vectorized RD-schemes
!     --------


!        Explicit arguments :  
!        --------------------   


!        Implicit arguments :     N1.d0
!        --------------------

!     Method.
!     -------
!     Estimate the max. integration time step and amount of iterations ...

!     Externals.  Interpolate any point within 1.d0 element based on linear shape functions ...
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Aron Roland 2011  *ECMWF* based on Aron Roland & Mathieu Dutour 

!     Modifications.
!     --------------
!        Original : Aron Roland, 11,2011
!     --------------------------------------------------------------

      IMPLICIT NONE

      REAL(KIND=JWRU), INTENT(IN)  :: X(3), Y(3), Z(3)
      REAL(KIND=JWRU), INTENT(IN)  :: XP, YP
      REAL(KIND=JWRU), INTENT(OUT) :: Zi
      REAL(KIND=JWRU), INTENT(OUT) :: WI(3)

      REAL(KIND=JWRU)               :: y1,y2,y3,x1,x2,x3,z1,z2,z3
      REAL(KIND=JWRU), SAVE         :: A,B,C,D
      REAL(KIND=JWRU), PARAMETER    :: THR = TINY(1.0_JWRU)

      LOGICAL, INTENT(IN)  :: LSAME

      IF (.NOT. LSAME) THEN
        x1 = X(1); x2 = X(2); x3 = X(3)
        y1 = Y(1); y2 = Y(2); y3 = Y(3)
        z1 = Z(1); z2 = Z(2); z3 = Z(3)
        A = y1*(z2 - z3)  +  y2*(z3 - z1) +  y3*(z1 - z2)
        B = z1*(x2 - x3)  +  z2*(x3 - x1) +  z3*(x1 - x2)
        C = x1*(y2 - y3)  +  x2*(y3 - y1) +  x3*(y1 - y2)
        D = -A*x1 - B*y1 - C*z1
        IF (ABS(C) .GT. THR ) THEN
          WI(1) = -A/C
          WI(2) = -B/C
          WI(3) = -D/C
        ELSE
          WI    = 0.0_JWRU
        END IF 
      END IF
      Zi = REAL(WI(1) * REAL(XP,JWRU) + WI(2) * REAL(YP,JWRU) + WI(3))

      END SUBROUTINE INTELEMENT
!**********************************************************************
!*                                                                    *
!**********************************************************************
END MODULE
