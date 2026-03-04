! (C) Copyright 2000- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE GATH_SPEC_CONTROL_MOD
CONTAINS
SUBROUTINE GATH_SPEC_CONTROL(PSPECG,KFGATHG,KTO,KVSET,PSPEC,LDIM1_IS_FLD,&
 & KSMAX,KSPEC2,KSPEC2_G,KPOSSP,KDIM0G,LDZA0IP)

!**** *GATH_SPEC_CONTROL* - Gather global spectral array from processors

!     Purpose.
!     --------
!        Routine for gathering spectral array

!**   Interface.
!     ----------
!     CALL GATH_SPEC_CONTROL(...)

!     Explicit arguments :
!     --------------------
!     PSPECG(:,:) - Global spectral array
!     KFGATHG     - Global number of fields to be distributed
!     KTO(:)    - Processor responsible for distributing each field
!     KVSET(:)    - "B-Set" for each field
!     PSPEC(:,:)  - Local spectral array
!     LDZA0IP     - Set first coefficients (imaginary part) to zero

!     ------------------------------------------------------------------


USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE MPL_MODULE  ,ONLY : MPL_RECV, MPL_SEND, MPL_BARRIER, MPL_WAIT,     &
     &                  JP_BLOCKING_STANDARD, JP_NON_BLOCKING_STANDARD

!USE TPM_GEN
!USE TPM_DIM
USE TPM_DISTR       ,ONLY : MTAGDISTSP, NPRCIDS, NPRTRW,  &
     &                      MYSETV, MYSETW, MYPROC, NPROC
USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS

USE SET2PE_MOD      ,ONLY : SET2PE
!USE SUWAVEDI_MOD
!

IMPLICIT NONE

REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT)  :: PSPECG(:,:)
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KFGATHG
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KTO(:)
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KVSET(:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(IN)  :: PSPEC(:,:)
LOGICAL            ,OPTIONAL, INTENT(IN)  :: LDIM1_IS_FLD
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KSMAX
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KSPEC2
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KSPEC2_G
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KPOSSP(:)
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KDIM0G(0:)
LOGICAL            ,OPTIONAL, INTENT(IN)  :: LDZA0IP

REAL(KIND=JPRB)    :: ZFLD(KSPEC2,KFGATHG),ZDUM(KSPEC2)
REAL(KIND=JPRB),ALLOCATABLE :: ZRECV(:,:)
INTEGER(KIND=JPIM) :: JM,JN,II,IFLDR,IFLDS,JFLD,ITAG,IBSET,ILEN,JA,ISND
INTEGER(KIND=JPIM) :: IRCV,ISP,ILENR,ISTA,ISTP,ISENDREQ(KFGATHG),IPOS0,JNM
INTEGER(KIND=JPIM) :: IDIST(KSPEC2_G),IMYFIELDS
LOGICAL            :: LLZA0IP

!     ------------------------------------------------------------------

LLZA0IP=.TRUE.
IF (PRESENT (LDZA0IP)) LLZA0IP=LDZA0IP

!GATHER SPECTRAL ARRAY

IF( NPROC == 1 ) THEN
  CALL GSTATS(1644,0)
  IF(LDIM1_IS_FLD) THEN
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JM,JFLD)
    DO JM=1,KSPEC2_G
      DO JFLD=1,KFGATHG
        PSPECG(JFLD,JM) =PSPEC(JFLD,JM)
      ENDDO
    ENDDO
!$OMP END PARALLEL DO
  ELSE
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JM,JFLD)
    DO JFLD=1,KFGATHG
      DO JM=1,KSPEC2_G
        PSPECG(JM,JFLD) =PSPEC(JM,JFLD)
      ENDDO
    ENDDO
!$OMP END PARALLEL DO
  ENDIF
  CALL GSTATS(1644,1)
ELSE
  IMYFIELDS = 0
  DO JFLD=1,KFGATHG
    IF(KTO(JFLD) == MYPROC) THEN
      IMYFIELDS = IMYFIELDS+1
    ENDIF
  ENDDO
  IF(IMYFIELDS>0) THEN
    ALLOCATE(ZRECV(KSPEC2_G,IMYFIELDS))
    II = 0
    CALL GSTATS(1804,0)
    DO JM=0,KSMAX
      DO JN=JM,KSMAX
        IDIST(II+1) = KDIM0G(JM)+(JN-JM)*2
        IDIST(II+2) = KDIM0G(JM)+(JN-JM)*2+1
        II = II+2
      ENDDO
    ENDDO
    CALL GSTATS(1804,1)
  ENDIF

  CALL GSTATS_BARRIER(788)

  !Send
  CALL GSTATS(810,0)
  IFLDS = 0
  IF(KSPEC2 > 0 )THEN
    DO JFLD=1,KFGATHG

      IBSET = KVSET(JFLD)
      IF( IBSET == MYSETV )THEN

        IFLDS = IFLDS+1
        ISND  = KTO(JFLD)
        ITAG  = MTAGDISTSP+JFLD+17
        IF(LDIM1_IS_FLD) THEN
          ZFLD(1:KSPEC2,IFLDS)=PSPEC(IFLDS,1:KSPEC2)
          CALL MPL_SEND(ZFLD(1:KSPEC2,IFLDS),KDEST=NPRCIDS(ISND),KTAG=ITAG,&
           &KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=ISENDREQ(JFLD),&
           &CDSTRING='GATH_SPEC_CONTROL')
        ELSE
          CALL MPL_SEND(PSPEC(1:KSPEC2,IFLDS),KDEST=NPRCIDS(ISND),KTAG=ITAG,&
           &KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=ISENDREQ(JFLD),&
           &CDSTRING='GATH_SPEC_CONTROL')
        ENDIF
      ENDIF
    ENDDO
  ENDIF

  ! Recieve
  IFLDR = 0
  DO JFLD=1,KFGATHG
    IF(KTO(JFLD) == MYPROC) THEN
      IBSET = KVSET(JFLD)
      IFLDR = IFLDR+1
      DO JA=1,NPRTRW
        ILEN = KPOSSP(JA+1)-KPOSSP(JA)
        IF( ILEN > 0 )THEN
          CALL SET2PE(IRCV,0,0,JA,IBSET)
          ITAG = MTAGDISTSP+JFLD+17
          ISTA = KPOSSP(JA)
          ISTP = ISTA+ILEN-1
          CALL MPL_RECV(ZRECV(ISTA:ISTP,IFLDR),KSOURCE=NPRCIDS(IRCV),KTAG=ITAG,&
           &KMP_TYPE=JP_BLOCKING_STANDARD,KOUNT=ILENR, &
           &CDSTRING='GATH_SPEC_CONTROL')
          IF( ILENR /= ILEN )THEN
            WRITE(0,'("GATH_SPEC_CONTROL: JFLD=",I4," JA=",I4," ILEN=",I10," ILENR=",I10)')&
            &JFLD,JA,ILEN,ILENR
            CALL ABORT_TRANS('GATH_SPEC_CONTROL:INVALID RECEIVE MESSAGE LENGTH')
          ENDIF
        ENDIF
      ENDDO
    ENDIF
  ENDDO

  ! Check for completion of sends
  IF(KSPEC2 > 0 )THEN
    DO JFLD=1,KFGATHG
      IBSET = KVSET(JFLD)
      IF( IBSET == MYSETV )THEN
        CALL MPL_WAIT(KREQUEST=ISENDREQ(JFLD), &
         & CDSTRING='GATH_GRID_CTL: WAIT')
      ENDIF
    ENDDO
  ENDIF
  CALL GSTATS(810,1)
  CALL GSTATS_BARRIER2(788)

  CALL GSTATS(1644,0)
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JFLD,JNM,II,JN,ISP)
  DO JFLD=1,IMYFIELDS
    IF(LDIM1_IS_FLD) THEN
      DO JNM=1,KSPEC2_G
        PSPECG(JFLD,JNM) = ZRECV(IDIST(JNM),JFLD)
      ENDDO
      IF (LLZA0IP) THEN
        II = 0
        DO JN=0,KSMAX
          ISP = KDIM0G(0)+JN*2+1
          II = II+2
          PSPECG(JFLD,II) = 0.0_JPRB
        ENDDO
      ENDIF
    ELSE
      DO JNM=1,KSPEC2_G
        PSPECG(JNM,JFLD) = ZRECV(IDIST(JNM),JFLD)
      ENDDO
      IF (LLZA0IP) THEN
        II = 0
        DO JN=0,KSMAX
          ISP = KDIM0G(0)+JN*2+1
          II = II+2
          PSPECG(II,JFLD) = 0.0_JPRB
        ENDDO
      ENDIF
    ENDIF
  ENDDO
!$OMP END PARALLEL DO
  CALL GSTATS(1644,1)
  IF(ALLOCATED(ZRECV)) DEALLOCATE(ZRECV)

  !Synchronize processors
  CALL GSTATS(785,0)
  CALL MPL_BARRIER(CDSTRING='GATH_SPEC_CONTROL:')
  CALL GSTATS(785,1)
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE GATH_SPEC_CONTROL
END MODULE GATH_SPEC_CONTROL_MOD


