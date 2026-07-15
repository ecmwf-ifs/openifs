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
 & KSMAX,KSPEC2,KSPEC2MX,KSPEC2G,KPOSSP,KDIM0G,KUMPP,KALLMS,KPTRMS,KN,LDZA0IP)

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
!     KFGATHG     - Global number of fields to be gathered
!     KTO(:)    - Processor responsible for distributing each field
!     KVSET(:)    - "B-Set" for each field
!     PSPEC(:,:)  - Local spectral array
!     LDIM1_IS_FLD - .TRUE. if first dimension contains the fields
!     KSMAX        - Spectral truncation limit
!     KSPEC2       - Local number of spectral coefficients
!     KSPEC2MX     - Maximum local number of spectral coefficients
!     KSPEC2G      - Global number of spectral coefficients
!     KPOSSP       - Position of local waves for each task
!     KDIM0G       - Defines partitioning of global spectral fields among PEs
!     KUMPP        - Number of spectral waves on this a-set
!     KALLMS       - Wave numbers for all a-set concatenated together to give all wave numbers in a-set order
!     KPTRMS       - Pointer to the first wave number of a given a-set in kallms array.
!     KN           - Number of spectral coefficients for each m wave
!     LDZA0IP     - Set first coefficients (imaginary part) to zero (global model only)

!     Externals.  SET2PE - compute "A and B" set from PE
!     ----------  MPL..  - message passing routines

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 2000-04-01
!    R. El Khatib 02-Dec-2020 re-write for optimizations and merge with LAM counterpart
!     ------------------------------------------------------------------


USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE MPL_MODULE  ,ONLY : MPL_RECV, MPL_SEND, MPL_BARRIER, MPL_WAIT, &
     &                  JP_BLOCKING_STANDARD, JP_NON_BLOCKING_STANDARD

USE TPM_DISTR       ,ONLY : MTAGDISTSP, NPRCIDS, NPRTRW, MYSETV, MYSETW, MYPROC, NPROC, NPRTRV
USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS
USE SET2PE_MOD      ,ONLY : SET2PE
USE TPM_GEOMETRY    ,ONLY : G

IMPLICIT NONE

REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT)  :: PSPECG(:,:)
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KFGATHG
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KTO(:)
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KVSET(:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(IN)  :: PSPEC(:,:)
LOGICAL            ,OPTIONAL, INTENT(IN)  :: LDIM1_IS_FLD
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KSMAX
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KSPEC2
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KSPEC2MX
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KSPEC2G
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KPOSSP(:)
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KDIM0G(0:)
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KUMPP(NPRTRW)
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KALLMS(KSMAX+1)
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KPTRMS(NPRTRW)
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KN(0:KSMAX)
LOGICAL            ,OPTIONAL, INTENT(IN)  :: LDZA0IP

REAL(KIND=JPRB) :: ZBUFSEND(KSPEC2MX,COUNT(KVSET(1:KFGATHG) == MYSETV))
REAL(KIND=JPRB) :: ZRECV(KSPEC2MX,COUNT(KTO(1:KFGATHG) == MYPROC))
INTEGER(KIND=JPIM) :: IASM0G(0:KSMAX)
INTEGER(KIND=JPIM) :: JM,JN,II,IFLDR,IFLDS,JFLD,ITAG,IB,ILEN(NPRTRW),JA,JB,ISND,JMLOC
INTEGER(KIND=JPIM) :: IPE(NPRTRV,NPRTRW),ILENR,ISENDREQ(NPROC),IPOSSP,JNM,JROC
INTEGER(KIND=JPIM) :: IFLD,IFLDLOC(COUNT(KTO(1:KFGATHG) == MYPROC)),IOFFPROC
INTEGER(KIND=JPIM) :: ILOCFLD(COUNT(KVSET(1:KFGATHG) == MYSETV))
LOGICAL            :: LLZA0IP

!     ------------------------------------------------------------------

! Compute help array for distribution

DO JA=1,NPRTRW
  ILEN(JA) = KPOSSP(JA+1)-KPOSSP(JA)
ENDDO
DO JA=1,NPRTRW
  DO JB=1,NPRTRV
    CALL SET2PE(IPE(JB,JA),0,0,JA,JB)
  ENDDO
ENDDO
IASM0G(0)=1
DO JM=1,KSMAX
  IASM0G(JM)=IASM0G(JM-1)+KN(JM-1)
ENDDO

LLZA0IP=.NOT.G%LAM ! or it should have been coded in the original code, please :-(
IF (PRESENT (LDZA0IP)) LLZA0IP=LDZA0IP

!GATHER SPECTRAL ARRAY

!Send
ISND=0
IOFFPROC=0
IF (KSPEC2 > 0) THEN
  CALL GSTATS(810,0)
  DO JROC=1,NPROC
    IF (JROC /= MYPROC) THEN
      IFLD=0 ! counter of fields in PSPEC
      IFLDS=0 ! counter of fields in ZBUFSEND
      DO JFLD=1,KFGATHG
        IF (KVSET(JFLD) == MYSETV) THEN
          IFLD=IFLD+1
          IF (JROC==KTO(JFLD)) THEN
            IFLDS=IFLDS+1
            IF (LDIM1_IS_FLD) THEN
              ZBUFSEND(1:KSPEC2,IOFFPROC+IFLDS)=PSPEC(IFLD,1:KSPEC2)
            ELSE
              ZBUFSEND(1:KSPEC2,IOFFPROC+IFLDS)=PSPEC(1:KSPEC2,IFLD)
            ENDIF
          ENDIF
        ENDIF
      ENDDO
      IF (IFLDS > 0) THEN
        ITAG=MTAGDISTSP+MYPROC
        ISND=ISND+1
        CALL MPL_SEND(ZBUFSEND(:,IOFFPROC+1:IOFFPROC+IFLDS),KDEST=NPRCIDS(JROC),KTAG=ITAG,&
         & KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=ISENDREQ(ISND),&
         & CDSTRING='GATH_SPEC_CONTROL')
      ENDIF
      IOFFPROC=IOFFPROC+IFLDS
    ENDIF
  ENDDO
  CALL GSTATS(810,1)

! Myself :
  IFLD=0
  IFLDR=0
  DO JFLD=1,KFGATHG
    IF (KTO(JFLD) == MYPROC) THEN
      IFLD=IFLD+1
      IF (KVSET(JFLD)==MYSETV) THEN
        IFLDR = IFLDR+1
        IFLDLOC(IFLDR)=IFLD
      ENDIF 
    ENDIF 
  ENDDO
  IFLD=0
  IFLDR=0
  DO JFLD=1,KFGATHG
    IF (KVSET(JFLD)==MYSETV) THEN
      IFLD=IFLD+1
      IF (KTO(JFLD) == MYPROC) THEN
        IFLDR = IFLDR+1
        ILOCFLD(IFLDR)=IFLD
      ENDIF
    ENDIF 
  ENDDO
  IF (IFLDR > 0) THEN
    IF (LDIM1_IS_FLD) THEN
      CALL GSTATS(1644,0)
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JFLD,JMLOC,JM,IPOSSP,II,JN)
      DO JFLD=1,IFLDR
        DO JMLOC=1,KUMPP(MYSETW)
          JM=KALLMS(KPTRMS(MYSETW)+JMLOC-1)
          IPOSSP=KDIM0G(JM)-KPOSSP(MYSETW)+1
          PSPECG(IFLDLOC(JFLD),IASM0G(JM):IASM0G(JM)+KN(JM)-1)=PSPEC(ILOCFLD(JFLD),IPOSSP:IPOSSP+KN(JM)-1)
        ENDDO
        IF (LLZA0IP) THEN
          II = 0
          DO JN=0,KSMAX
            II = II+2
            PSPECG(IFLDLOC(JFLD),II) = 0.0_JPRB
          ENDDO
        ENDIF
      ENDDO
!$OMP END PARALLEL DO
      CALL GSTATS(1644,1)
    ELSE
      CALL GSTATS(1644,0)
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JFLD,JMLOC,JM,IPOSSP,II,JN)
      DO JFLD=1,IFLDR
        DO JMLOC=1,KUMPP(MYSETW)
          JM=KALLMS(KPTRMS(MYSETW)+JMLOC-1)
          IPOSSP=KDIM0G(JM)-KPOSSP(MYSETW)+1
          PSPECG(IASM0G(JM):IASM0G(JM)+KN(JM)-1,IFLDLOC(JFLD))=PSPEC(IPOSSP:IPOSSP+KN(JM)-1,ILOCFLD(JFLD))
        ENDDO
        IF (LLZA0IP) THEN
          II = 0
          DO JN=0,KSMAX
            II = II+2
            PSPECG(II,IFLDLOC(JFLD)) = 0.0_JPRB
          ENDDO
        ENDIF
      ENDDO
!$OMP END PARALLEL DO
      CALL GSTATS(1644,1)
    ENDIF
  ENDIF

ENDIF

! Receive
DO JA=1,NPRTRW
  IF (ILEN(JA) > 0) THEN
    DO JB=1,NPRTRV
      IF (IPE(JB,JA) /= MYPROC) THEN
        ! Locate received fields in source array :
        IFLD=0
        IFLDR=0
        DO JFLD=1,KFGATHG
          IF (KTO(JFLD) == MYPROC) THEN
            IFLD=IFLD+1
            IF (KVSET(JFLD)==JB) THEN
              IFLDR = IFLDR+1
              IFLDLOC(IFLDR)=IFLD
            ENDIF 
          ENDIF 
        ENDDO
        IF (IFLDR > 0) THEN
          ITAG=MTAGDISTSP+IPE(JB,JA)
          CALL GSTATS(810,0)
          CALL MPL_RECV(ZRECV(:,1:IFLDR),KSOURCE=NPRCIDS(IPE(JB,JA)),KTAG=ITAG,&
           & KMP_TYPE=JP_BLOCKING_STANDARD,KOUNT=ILENR, &
           & CDSTRING='GATH_SPEC_CONTROL')
          IF (ILENR /= KSPEC2MX*IFLDR) THEN
            CALL ABORT_TRANS('GATH_SPEC_CONTROL:INVALID RECEIVE MESSAGE LENGTH')
          ENDIF
          CALL GSTATS(810,1)
          CALL GSTATS(1644,0)
          IF (LDIM1_IS_FLD) THEN
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JFLD,JMLOC,JM,IPOSSP,II,JN)
            DO JFLD=1,IFLDR
              DO JMLOC=1,KUMPP(JA)
                JM=KALLMS(KPTRMS(JA)+JMLOC-1)
                IPOSSP=KDIM0G(JM)-KPOSSP(JA)+1
                PSPECG(IFLDLOC(JFLD),IASM0G(JM):IASM0G(JM)+KN(JM)-1)=ZRECV(IPOSSP:IPOSSP+KN(JM)-1,JFLD)
              ENDDO
              IF (LLZA0IP) THEN
                II = 0
                DO JN=0,KSMAX
                  II = II+2
                  PSPECG(IFLDLOC(JFLD),II) = 0.0_JPRB
                ENDDO
              ENDIF
            ENDDO
!$OMP END PARALLEL DO
          ELSE
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JFLD,JMLOC,JM,IPOSSP,II,JN)
            DO JFLD=1,IFLDR
              DO JMLOC=1,KUMPP(JA)
                JM=KALLMS(KPTRMS(JA)+JMLOC-1)
                IPOSSP=KDIM0G(JM)-KPOSSP(JA)+1
                PSPECG(IASM0G(JM):IASM0G(JM)+KN(JM)-1,IFLDLOC(JFLD))=ZRECV(IPOSSP:IPOSSP+KN(JM)-1,JFLD)
              ENDDO
              IF (LLZA0IP) THEN
                II = 0
                DO JN=0,KSMAX
                  II = II+2
                  PSPECG(II,IFLDLOC(JFLD)) = 0.0_JPRB
                ENDDO
              ENDIF
            ENDDO
!$OMP END PARALLEL DO
          ENDIF
          CALL GSTATS(1644,1)
        ENDIF
      ENDIF
    ENDDO
  ENDIF
ENDDO
CALL GSTATS_BARRIER2(788)

! Check for completion of sends
CALL GSTATS(810,0)
IF (ISND > 0) THEN
  CALL MPL_WAIT(ISENDREQ(1:ISND),CDSTRING='GATH_GRID_CTL: WAIT')
ENDIF
CALL GSTATS(810,1)

!Synchronize processors. Useful ??
CALL GSTATS(785,0)
!rekCALL MPL_BARRIER(CDSTRING='GATH_SPEC_CONTROL:')
CALL GSTATS(785,1)

CALL GSTATS_BARRIER(788)

!     ------------------------------------------------------------------

END SUBROUTINE GATH_SPEC_CONTROL
END MODULE GATH_SPEC_CONTROL_MOD
