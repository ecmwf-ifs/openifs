! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

      SUBROUTINE OUTGRID(BOUT)
! ----------------------------------------------------------------------

!*    PURPOSE.
!     --------

!    BOUT (output of integrated parameters per block)  => GOUT (global ready for output) 

!**   INTERFACE.
!     ----------

!        *CALL* *OUTGRID(BOUT)
!          *BOUT*    - OUTPUT PARAMETERS BUFFER

! ----------------------------------------------------------------------
      
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUT  , ONLY : JPPFLAG  ,IPFGTBL ,NIPRMOUT ,ITOBOUT
      USE YOWGRID  , ONLY : IJSLOC   ,IJLLOC  ,NPROMA_WAM, NCHNK, KIJL4CHNK
      USE YOWINTP  , ONLY : GOUT
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,MPMAXLENGTH ,KTAG
      USE YOWPARAM , ONLY : NIBLO    ,NGX      ,NGY
      USE YOWPCONS , ONLY : ZMISS
      USE YOWSPEC  , ONLY : NBLKS    ,NBLKE
      USE YOWTEST  , ONLY : IU06
      USE OUTPUT_STRUCT, ONLY : ARR_OUT_RECV, LocalPosICT
      USE OUTPUT_STRUCT, ONLY : INITIAL_OUTPUT_INITS_NEXTGEN
      USE OUTPUT_STRUCT, ONLY : SET_UP_ARR_OUT_RECV
      USE YOWUNPOOL, ONLY : NIBLO_FD, NIBLO_OUT, OUT_METHOD, LLUNSTR

      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
      USE MPL_MODULE, ONLY : MPL_RECV, MPL_SEND, MPL_WAIT, &
                           & JP_NON_BLOCKING_STANDARD

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "abort1.intfb.h"
#include "makegrid.intfb.h"

      REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, NIPRMOUT, NCHNK), INTENT(IN) :: BOUT


      INTEGER(KIND=JWIM) :: IJ, ITT, ICT, ITG, IFLD, IR, IPR, ICOUNT
      INTEGER(KIND=JWIM) :: ICHNK, IPRM
      INTEGER(KIND=JWIM) :: NFLDTOT, NFLDPPEMAX
      INTEGER(KIND=JWIM), DIMENSION(NPROC) :: ICNT, NFLDPPE
      INTEGER(KIND=JWIM), DIMENSION(NPROC+JPPFLAG) :: IREQ
      INTEGER(KIND=JWIM), DIMENSION(NPROC) :: ISENDCOUNTS,IRECVCOUNTS

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB) :: ZDUM(2)
      REAL(KIND=JWRB), ALLOCATABLE, DIMENSION(:) :: GTEMP
      REAL(KIND=JWRB), ALLOCATABLE, DIMENSION(:,:) :: ZSENDBUF, ZRECVBUF

      LOGICAL, SAVE :: HaveMPI_arrays = .FALSE.
!----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('OUTGRID',0,ZHOOK_HANDLE)

!     FIND THE NUMBER OF FIELDS EACH PE HAS TO DEAL WITH
      NFLDPPE(:)=0
      NFLDTOT=0
      NFLDPPEMAX=1
      DO ICT=1,JPPFLAG
        IF (IPFGTBL(ICT) > 0) THEN 
          NFLDPPE(IPFGTBL(ICT))=NFLDPPE(IPFGTBL(ICT))+1
          NFLDTOT=NFLDTOT+1
          NFLDPPEMAX=MAX(NFLDPPEMAX,NFLDPPE(IPFGTBL(ICT)))
        ENDIF
      ENDDO
      IF ( .NOT. HaveMPI_arrays ) THEN
        IF (LLUNSTR .AND. (OUT_METHOD == 2)) THEN
           CALL INITIAL_OUTPUT_INITS_NEXTGEN
        END IF
      END IF
      HaveMPI_arrays=.TRUE.
      IF (NFLDTOT == 0) THEN
        IF (LHOOK) CALL DR_HOOK('OUTGRID',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
 

!     SENDING TO RELEVANT PE'S
!     ------------------------

      DO IPR=1,NPROC
        ISENDCOUNTS(IPR)=NFLDPPE(IPR)*(NBLKE(IRANK)-NBLKS(IRANK)+1)
      ENDDO

      DO IPR=1,NPROC
        IRECVCOUNTS(IPR)=NFLDPPE(IRANK)*(NBLKE(IPR)-NBLKS(IPR)+1)
      ENDDO

      IF (.NOT.LLUNSTR) THEN
      
!     LOADING THE COMMUNICATION BUFFER
        ALLOCATE(ZSENDBUF(NFLDPPEMAX * MPMAXLENGTH,NPROC))

        ICNT(:)=0

        DO ICT = 1, JPPFLAG
          IPR = IPFGTBL(ICT)
          IF (IPR > 0) THEN 
            DO ICHNK = 1, NCHNK
              DO IPRM = 1, KIJL4CHNK(ICHNK)
                ICNT(IPR) = ICNT(IPR) + 1
                ZSENDBUF(ICNT(IPR), IPR) = BOUT(IPRM, ITOBOUT(ICT), ICHNK)
              ENDDO
            ENDDO
          ENDIF
        ENDDO

      ELSEIF (OUT_METHOD == 1) THEN
      
!!!!! this will need to be adapted to use BOUT
        WRITE(0,*) '!!! ********************************* !!'
        WRITE(0,*) '!!! in outgrid. Not yet ready !!!' 
        WRITE(0,*) '!!! ********************************* !!'
        CALL ABORT1

!     LOADING THE COMMUNICATION BUFFER
        ALLOCATE(ZSENDBUF(NFLDPPEMAX * MPMAXLENGTH,NPROC))

!        ICNT(:)=0
!        DO ICT=1,JPPFLAG
!          IPR=IPFGTBL(ICT)
!          IF (IPR > 0) THEN 
!            DO IJ=IJSLOC,IJLLOC
!              ICNT(IPR) = ICNT(IPR) + 1
!              ZSENDBUF(ICNT(IPR),IPR) = BOUT(IJ,ITOBOUT(ICT)) !!!!!
!            ENDDO
!          ENDIF
!        ENDDO
      ELSE

!!!!! this will need to be adapted to use BOUT
        WRITE(0,*) '!!! ********************************* !!'
        WRITE(0,*) '!!! in outgrid. Not yet ready !!!' 
        WRITE(0,*) '!!! ********************************* !!'
        CALL ABORT1

!!!        CALL SET_UP_ARR_OUT_RECV(IJSLOC, IJLLOC, BOUT(IJSLOC:IJLLOC,:), NFLDPPE)
      END IF


!     GLOBAL EXCHANGE

      IF ((.NOT.LLUNSTR) .OR. (OUT_METHOD == 1)) THEN
        CALL GSTATS(693,0)
        IR=0
        IF (NFLDPPE(IRANK) > 0) THEN
          ALLOCATE(ZRECVBUF(NFLDPPE(IRANK)*MPMAXLENGTH,NPROC))
          DO IPR=1,NPROC
            IR=IR+1
            CALL MPL_RECV(ZRECVBUF(1:IRECVCOUNTS(IPR),IPR),KSOURCE=IPR, &
     &        KTAG=KTAG,                                                &
     &        KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=IREQ(IR),      &
     &        CDSTRING='OUTGRID:')
          ENDDO
        ENDIF

        DO IPR=1,NPROC
          IF (NFLDPPE(IPR) > 0) THEN
            IR=IR+1
            CALL MPL_SEND(ZSENDBUF(1:ICNT(IPR),IPR),KDEST=IPR,          &
     &        KTAG=KTAG,                                                &
     &        KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=IREQ(IR),      &
     &        CDSTRING='OUTGRID:')
          ENDIF
        ENDDO
!       NOW WAIT FOR ALL TO COMPLETE
        CALL MPL_WAIT(KREQUEST=IREQ(1:IR),CDSTRING='OUTGRID:')
        CALL GSTATS(693,1)
      END IF


!     RETRIEVE THE INFORMATION
!     ------------------------

      IF (LLUNSTR .AND. (OUT_METHOD == 2)) THEN
        NIBLO_OUT = NIBLO_FD
      ELSE
        NIBLO_OUT = NIBLO
      END IF
      ALLOCATE(GTEMP(NIBLO_OUT))


!     RECEIVING AND TRANSFERING FROM BLOCK TO GRID
!     --------------------------------------------

      IF (NFLDPPE(IRANK) > 0) THEN
        IF (ALLOCATED(GOUT)) DEALLOCATE(GOUT)
        ALLOCATE(GOUT(NFLDPPE(IRANK),NGX,NGY))
      ENDIF

      ICNT(:)=0
      ICT=1
      IFLD=1
      DO WHILE ( IFLD <= NFLDPPE(IRANK) .AND. ICT <= JPPFLAG )
        IF (IPFGTBL(ICT) == IRANK ) THEN
          IF (.NOT.(LLUNSTR) .OR. (OUT_METHOD == 1)) THEN
            DO IPR=1,NPROC
              ICOUNT=ICNT(IPR)
              DO IJ=NBLKS(IPR),NBLKE(IPR)
                ICOUNT=ICOUNT+1
                GTEMP(IJ)=ZRECVBUF(ICOUNT,IPR)
              ENDDO
              ICNT(IPR)=ICOUNT
            ENDDO
          ELSE
            GTEMP=ARR_OUT_RECV(LocalPosICT(ICT),:)
          END IF

          CALL MAKEGRID (GTEMP, GOUT(IFLD,:,:), ZMISS)
          IFLD=IFLD+1

        ENDIF ! (IPFGTBL) 

        ICT=ICT+1
      ENDDO  ! WHILE

      DEALLOCATE(GTEMP)
      IF (ALLOCATED(ZRECVBUF)) DEALLOCATE(ZRECVBUF)
      IF (ALLOCATED(ZSENDBUF)) DEALLOCATE(ZSENDBUF)

      IF (LHOOK) CALL DR_HOOK('OUTGRID',1,ZHOOK_HANDLE)

      END SUBROUTINE OUTGRID
