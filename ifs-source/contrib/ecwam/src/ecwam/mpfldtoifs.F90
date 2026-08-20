! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE MPFLDTOIFS(IJS, IJL, BLK2GLO, LLINIT_GRID, NWVFIELDS, BLOCK,       &
 &                    GRID, DEFVAL, MASK_OUT, LLGLOBAL)

!****  *MPFLDTOIFS* - TRANSFORMS BLOCK DATA TO GRID DATA FOR
!****                 FIELDS THAT WILL BE RETURNED TO IFS.                

!     J. BIDLOT    ECMWF AUGUST 2008

!     PURPOSE.
!     --------

!     TRANSFORMS BLOCK DATA TO GRID DATA.

!*    INTERFACE.
!     ----------

!     CALL *MPFLDTOIFS(IJS, IJL, BLK2GLO, LLINIT_GRID, NWVFIELDS, BLOCK,
!    &                 GRID, DEFVAL, MASK_OUT )*
!         *IJS*        - BLOCK INDEX OF FIRST GRIDPOINT.
!         *IJL*        - BLOCK INDEX OF LAST GRIDPOINT.
!         *BLK2GLO*    - BLOCK TO GRID TRANSFORMATION
!         *LLINIT_GRID*- IF TRUE GROBAL ARRAY GRID IS INITIALISED GLOBALLY. !!! GRID IS INTENT(INOUT) !!!!
!         *NWVFIELDS*  - TOTAl NUMBER OF FIELDS RETURNED TO IFS
!         *BLOCK*      - FIELDS IN BLOCK FORM (INPUT)
!                        ONLY DEFINED LOCALLY.
!         *GRID*       - FIELDS IN GRID FORM (OUTPUT)
!                        DEFINED  GLOBALLY OR JUST ON THE
!                        SUBAREA NEEDED (SEE MASK_OUT) DEPENDING
!                        ON LWCOUNORMS OR LMASK_OUT_NOT_SET
!         *DEFVAL*     - DEFAULT VALUE TO ASSIGN EACH FIELD WHEN
!                        NOT DEFINED BY WAVE MODEL.
!         *MASK_OUT*   - MASK POINTING TO VALUES OF FIELDS THAT ARE
!                        NEEDED BY IFS ON CURRENT TASK.
!         *LLGLOBAL*   - TRUE IF GRID HAS BEEN DEFINED GLOBALLY (OUTPUT).


!     METHOD.
!     -------
!     MPL ALLGATHERV OF ARRAY FIELD TO ALL PE's FOR FIRST CALL
!     FOR SUBSEQUENT CALLS COMMUNICATION TAKES PLACE WITH ONLY THE
!     TASKS THAT WE NEED TO SEND DATA TO.

!     REFERENCES.
!     -----------
!         NONE
! -------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWDRVTYPE  , ONLY : WVGRIDGLO

      USE YOWCOUP  , ONLY : LWCOUNORMS,LMASK_OUT_NOT_SET,               &
     &            LMASK_TASK_STR,                                       &
     &            LFROMTASK,IJFROMTASK,NFROMTASKS,ISTFROMTASK,          &
     &            LTOTASK  ,IJTOTASK  ,NTOTASKS  ,ISTTOTASK
      USE YOWMAP   , ONLY : NGX      ,NGY      ,NIBLO
      USE YOWMPP   , ONLY : IRANK    ,NPROC
      USE YOWSPEC  , ONLY : NSTART   ,NEND
      USE YOWTEST  , ONLY : IU06

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
      USE MPL_MODULE, ONLY : MPL_ALLGATHERV, MPL_SEND, MPL_RECV, MPL_WAIT, &
                           & JP_NON_BLOCKING_STANDARD, JP_BLOCKING_STANDARD
!----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL
      TYPE(WVGRIDGLO), INTENT(IN) :: BLK2GLO
      LOGICAL, INTENT(IN) :: LLINIT_GRID
      INTEGER(KIND=JWIM), INTENT(IN) :: NWVFIELDS
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NWVFIELDS), INTENT(IN) :: BLOCK
      REAL(KIND=JWRB), DIMENSION(:,:,:), INTENT(INOUT) :: GRID
      REAL(KIND=JWRB), DIMENSION(NWVFIELDS), INTENT(IN) :: DEFVAL 
      INTEGER(KIND=JWIM), INTENT(IN) :: MASK_OUT(NGX,NGY)
      LOGICAL, INTENT(OUT) :: LLGLOBAL


      INTEGER(KIND=JWIM) :: IFLD, IP, IP1, IJ, I, J, IX, IY
      INTEGER(KIND=JWIM) :: ITAG, IREQ, ICOUNT, IC, IST, IEND
      INTEGER(KIND=JWIM), DIMENSION(NPROC) :: KRECVCOUNTS
      INTEGER(KIND=JWIM) :: ISENDREQ(NPROC)

      REAL(KIND=JWRB), ALLOCATABLE :: ZBUFS(:) ,ZBUFR(:)
      REAL(KIND=JWRB), ALLOCATABLE :: ZSENDBUF(:), ZRECVBUF(:,:) 
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE1

!----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('MPFLDTOIFS',0,ZHOOK_HANDLE)

      CALL GSTATS_BARRIER(734)

      LLGLOBAL=(LWCOUNORMS .OR. LMASK_OUT_NOT_SET)

      IF (LLINIT_GRID) THEN
        CALL GSTATS(1503,0)
!       DEFAULT VALUES
!$OMP   PARALLEL DO SCHEDULE(STATIC) PRIVATE(IFLD,J,I) COLLAPSE(2)
        DO IFLD=1,NWVFIELDS
          DO J = 1,NGY
            DO I = 1,NGX
              GRID(I,J,IFLD) = DEFVAL(IFLD) 
            ENDDO
          ENDDO
        ENDDO
!$OMP   END PARALLEL DO
        CALL GSTATS(1503,1)
      ENDIF


      CALL GSTATS(686,0)

!     LOOP OVER ALL INPUT FIELDS

      DO IFLD=1,NWVFIELDS

        IF (NPROC > 1) THEN

!         GLOBAL EXCHANGE

          IF (LLGLOBAL) THEN

            IF(.NOT. ALLOCATED(ZBUFS)) ALLOCATE(ZBUFS(NIBLO))
            IF(.NOT. ALLOCATED(ZBUFR)) ALLOCATE(ZBUFR(NIBLO))

            DO IJ=IJS,IJL
              ZBUFS(IJ)=BLOCK(IJ,IFLD)
            ENDDO

            DO IP=1,NPROC
              KRECVCOUNTS(IP)=NEND(IP) - NSTART(IP)+1
            ENDDO

            CALL MPL_ALLGATHERV(ZBUFS(NSTART(IRANK):NEND(IRANK)),       &
     &                          ZBUFR(1:NIBLO),KRECVCOUNTS,             &
     &                          CDSTRING='MPFLDTOIFS:')

!           TRANSFORM FROM BLOCK TO GLOBAL GRID
            DO IP=1,NPROC
              DO IJ = NSTART(IP), NEND(IP)
               IX = BLK2GLO%IXLG(IJ)
               IY = NGY- BLK2GLO%KXLT(IJ) +1
               GRID(IX,IY,IFLD) = ZBUFR(IJ)
              ENDDO
            ENDDO

            LMASK_OUT_NOT_SET=.FALSE.

          ELSE

!           INITIALISE THE TASK STRUCTURE FOR OPTIMAL DATA EXCHANGE
!           FOR FIELDS RETURNED TO IFS.

            IF (LMASK_TASK_STR) THEN

              LMASK_TASK_STR=.FALSE.

              IF (ALLOCATED(LFROMTASK)) DEALLOCATE(LFROMTASK)
              ALLOCATE(LFROMTASK(NPROC))

              DO IP=1,NPROC
                LFROMTASK(IP)=0
                DO IJ = NSTART(IP), NEND(IP)
                  IX = BLK2GLO%IXLG(IJ)
                  IY = NGY- BLK2GLO%KXLT(IJ) +1
                  IF (MASK_OUT(IX,IY) == 1) LFROMTASK(IP)=LFROMTASK(IP)+1
                ENDDO
              ENDDO

              NFROMTASKS=SUM(LFROMTASK)
              IF (ALLOCATED(IJFROMTASK)) DEALLOCATE(IJFROMTASK)
              ALLOCATE(IJFROMTASK(NFROMTASKS))

              IF (ALLOCATED(ISTFROMTASK)) DEALLOCATE(ISTFROMTASK)
              ALLOCATE(ISTFROMTASK(NPROC))
              DO IP=1,NPROC
                ISTFROMTASK(IP)=0
              ENDDO

              ICOUNT=0
              DO IP=1,NPROC
                IF (LFROMTASK(IP) > 0) THEN
                  ISTFROMTASK(IP)=ICOUNT+1
                  DO IJ = NSTART(IP), NEND(IP)
                    IX = BLK2GLO%IXLG(IJ)
                    IY = NGY- BLK2GLO%KXLT(IJ) +1
                    IF (MASK_OUT(IX,IY) == 1) THEN
                      ICOUNT=ICOUNT+1
                      IJFROMTASK(ICOUNT)=IJ
                    ENDIF
                  ENDDO
                ENDIF
              ENDDO

              IF (ALLOCATED(LTOTASK)) DEALLOCATE(LTOTASK)
              ALLOCATE(LTOTASK(NPROC))
              DO IP=1,NPROC
                LTOTASK(IP)=0
              ENDDO

!             SEND TO ALL OTHER TASKS TO LET THEM KNOW IF CONTRIBUTIONS
!             FROM THEM ARE NEEDED.
              ITAG=1
              DO IP=1,NPROC
                CALL MPL_SEND(LFROMTASK(IP),KDEST=IP,KTAG=ITAG,         &
     &                        KMP_TYPE=JP_NON_BLOCKING_STANDARD,        &
     &                        KREQUEST=ISENDREQ(IP),                    &
     &                        CDSTRING='MPFLDTOIFS: SEND COUNT ' )
              ENDDO
!             RECEIVE INFORMATION ON WHICH TASKS WILL NEED TO BE SENT
!             SOME CONTRIBUTIONS.
              DO IP=1,NPROC
                CALL MPL_RECV(LTOTASK(IP),KSOURCE=IP,KTAG=ITAG,         &
     &                        KMP_TYPE=JP_BLOCKING_STANDARD,            &
     &                        CDSTRING='MPFLDTOIFS: RECV COUNT ' )
              ENDDO

              CALL MPL_WAIT(KREQUEST=ISENDREQ(1:NPROC),                 &
     &                      CDSTRING='MPFLDTOIFS: WAIT SEND COUNT')


              NTOTASKS=SUM(LTOTASK)

              IF (ALLOCATED(IJTOTASK)) DEALLOCATE(IJTOTASK)
              ALLOCATE(IJTOTASK(NTOTASKS))

              IF (ALLOCATED(ISTTOTASK)) DEALLOCATE(ISTTOTASK)
              ALLOCATE(ISTTOTASK(NPROC))
              DO IP=1,NPROC
                ISTTOTASK(IP)=0
              ENDDO

!             SEND TO ALL OTHER TASKS TO LET THEM KNOW WHICH CONTRIBUTIONS
!             FROM THEM ARE NEEDED.
              ITAG=2
              IREQ=0

              DO IP=1,NPROC
                IF (LFROMTASK(IP) > 0) THEN
                  IREQ=IREQ+1

                  IST=ISTFROMTASK(IP)
                  IEND=IST+LFROMTASK(IP)-1
                  CALL MPL_SEND(IJFROMTASK(IST:IEND),                   &
     &                          KDEST=IP,KTAG=ITAG,                     &
     &                          KMP_TYPE=JP_NON_BLOCKING_STANDARD,      &
     &                          KREQUEST=ISENDREQ(IREQ),                &
     &                          CDSTRING='MPFLDTOIFS: SEND IJ ' )
                ENDIF
              ENDDO
!             RECEIVE INFORMATION ON WHAT CONTRIBUTIONS WILL NEED TO BE 
!             SENT.
              IST=1
              DO IP=1,NPROC
                IF (LTOTASK(IP) > 0) THEN
                  IEND=IST+LTOTASK(IP)-1
                  ISTTOTASK(IP)=IST
                  CALL MPL_RECV(IJTOTASK(IST:IEND),                     &
     &                          KSOURCE=IP,KTAG=ITAG,                   &
     &                          KMP_TYPE=JP_BLOCKING_STANDARD,          &
     &                          CDSTRING='MPFLDTOIFS: RECV IJ ' )
                  IST=IST+LTOTASK(IP)
                ENDIF
              ENDDO

              IF ( IREQ > 0 )THEN
                CALL MPL_WAIT(KREQUEST=ISENDREQ(1:IREQ),                &
     &                        CDSTRING='MPFLDTOIFS: WAIT SEND IJ')
              ENDIF

            ENDIF   !! LMASK_TASK_STR


!           EXCHANGE INFORMATION

!           SEND LOCAL TASK CONTRIBUTION TO TASKS THAT NEED IT

            IF(.NOT. ALLOCATED(ZSENDBUF)) ALLOCATE(ZSENDBUF(MAX(1,NTOTASKS)))
            IF(.NOT. ALLOCATED(ZRECVBUF)) ALLOCATE(ZRECVBUF(NFROMTASKS,NWVFIELDS))

            DO IP=1,NPROC
              IF (LTOTASK(IP) > 0) THEN
                IST=ISTTOTASK(IP)
                IEND=IST+LTOTASK(IP)-1
                DO IC=IST,IEND
                  ZSENDBUF(IC)=BLOCK(IJTOTASK(IC),IFLD)
                ENDDO
              ENDIF
            ENDDO

            ITAG=3
            IREQ=0

            DO IP=1,NPROC
              IF (LTOTASK(IP) > 0) THEN
                IREQ=IREQ+1
                IST=ISTTOTASK(IP)
                IEND=IST+LTOTASK(IP)-1
                CALL MPL_SEND(ZSENDBUF(IST:IEND),                       &
     &                        KDEST=IP,KTAG=ITAG,                       &
     &                        KMP_TYPE=JP_NON_BLOCKING_STANDARD,        &
     &                        KREQUEST=ISENDREQ(IREQ),                  &
     &                        CDSTRING='MPFLDTOIFS: SEND DATA')
              ENDIF
            ENDDO


!           RECEIVE INFORMATION FROM OTHER TASKS (IF NEEDED)
            DO IP=1,NPROC

              IF (LFROMTASK(IP) > 0) THEN
                IST=ISTFROMTASK(IP)
                IEND=IST+LFROMTASK(IP)-1
                CALL MPL_RECV(ZRECVBUF(IST:IEND,IFLD),                  &
     &                        KSOURCE=IP,KTAG=ITAG,                     &
     &                        KMP_TYPE=JP_BLOCKING_STANDARD,            &
     &                        CDSTRING='MPFLDTOIFS: RECV DATA ')

              ENDIF
            ENDDO

!           ENSURE ALL SENDS ARE FINISHED.

            IF (IREQ > 0) THEN
              CALL MPL_WAIT(KREQUEST=ISENDREQ(1:IREQ),                  &
     &                      CDSTRING='MPFLDTOIFS: WAIT SEND DATA ')
            ENDIF

!           THE TRANSFER OF ZRECVBUF TO GRID HAS BEEN MOVED OUTSIDE LOOP on IFLD

          ENDIF

        ELSE ! ONE PE

!         TRANSFORM FROM BLOCK TO GRID
          DO IP=1,NPROC
            DO IJ = NSTART(IP), NEND(IP)
             IX = BLK2GLO%IXLG(IJ)
             IY = NGY- BLK2GLO%KXLT(IJ) +1
             GRID(IX,IY,IFLD) = BLOCK(IJ,IFLD)
            ENDDO
          ENDDO

        ENDIF

      ENDDO !! IFLD

      IF(ALLOCATED(ZBUFS)) DEALLOCATE(ZBUFS)
      IF(ALLOCATED(ZBUFR)) DEALLOCATE(ZBUFR)

      IF(ALLOCATED(ZSENDBUF)) DEALLOCATE(ZSENDBUF)

      CALL GSTATS(686,1)

      IF(ALLOCATED(ZRECVBUF)) THEN
        CALL GSTATS(1503,0)
!       THE TRANSFER OF ZRECVBUF TO GRID HAS BEEN MOVED OUTSIDE THE MAIN LOOP on IFLD
!$OMP   PARALLEL DO SCHEDULE(STATIC) PRIVATE(IFLD,IP,IST,IEND,IC,IJ,IX,IY)
        DO IFLD=1,NWVFIELDS
          DO IP=1,NPROC
            IF (LFROMTASK(IP) > 0) THEN
              IST=ISTFROMTASK(IP)
              IEND=IST+LFROMTASK(IP)-1
!             TRANSFORM FROM BLOCK TO GRID
              DO IC = IST,IEND
                IJ=IJFROMTASK(IC)
                IX = BLK2GLO%IXLG(IJ)
                IY = NGY- BLK2GLO%KXLT(IJ) +1
                GRID(IX,IY,IFLD) = ZRECVBUF(IC,IFLD)
              ENDDO
            ENDIF
          ENDDO
        ENDDO
!$OMP   END PARALLEL DO

        DEALLOCATE(ZRECVBUF)
        CALL GSTATS(1503,1)
      ENDIF

IF (LHOOK) CALL DR_HOOK('MPFLDTOIFS',1,ZHOOK_HANDLE)

END SUBROUTINE MPFLDTOIFS
