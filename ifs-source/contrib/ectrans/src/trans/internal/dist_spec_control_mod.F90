! (C) Copyright 2000- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE DIST_SPEC_CONTROL_MOD
CONTAINS
SUBROUTINE DIST_SPEC_CONTROL(PSPECG,KFDISTG,KFROM,KVSET,PSPEC,LDIM1_IS_FLD,&
 & KSMAX,KSPEC2,KSPEC2MX,KSPEC2G,KPOSSP,KDIM0G,KUMPP,KALLMS,KPTRMS,KN,KSORT)

!**** *DIST_SPEC_CONTROL* - Distribute global spectral array among processors

!     Purpose.
!     --------
!        Routine for distributing spectral array

!**   Interface.
!     ----------
!     CALL DIST_SPEC_CONTROL(...)

!     Explicit arguments :
!     --------------------
!     PSPECG(:,:)  - Global spectral array
!     KFDISTG      - Global number of fields to be distributed
!     KFROM(:)     - Processor resposible for distributing each field
!     KVSET(:)     - "B-Set" for each field
!     PSPEC(:,:)   - Local spectral array
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
!     KSORT(:)     - Re-order fields on output

!     Externals.  SET2PE - compute "A and B" set from PE
!     ----------  MPL..  - message passing routines

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 2000-04-01
!    P.Marguinaud : 2014-10-10
!    R. El Khatib 25-Jul-2019 Optimization by vectorization, proper non-blocking comms 
!                             and overlapp send/recv with pack/unpack
!    R. El Khatib 02-Jun-2022 Optimization/Cleaning
!     ------------------------------------------------------------------

USE PARKIND1        ,ONLY : JPIM     ,JPRB
USE MPL_MODULE      ,ONLY : MPL_RECV, MPL_SEND, MPL_WAIT, MPL_WAITANY, JP_NON_BLOCKING_STANDARD
USE TPM_DISTR       ,ONLY : MTAGDISTSP, MYSETV, MYSETW, NPRCIDS, NPRTRW, MYPROC, NPROC, NPRTRV, D
USE SET2PE_MOD      ,ONLY : SET2PE
USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS

IMPLICIT NONE

REAL(KIND=JPRB)    ,OPTIONAL, INTENT(IN), CONTIGUOUS  :: PSPECG(:,:)
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KFDISTG
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KFROM(:)
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KVSET(:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT), CONTIGUOUS :: PSPEC(:,:)
LOGICAL                     , INTENT(IN)  :: LDIM1_IS_FLD
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KSMAX
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KSPEC2
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KSPEC2MX
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KSPEC2G
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KPOSSP(NPRTRW+1)
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KDIM0G(0:)
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KUMPP(NPRTRW)
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KALLMS(KSMAX+1)
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KPTRMS(NPRTRW)
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KN(0:KSMAX)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN), TARGET :: KSORT (:)
    
REAL(KIND=JPRB)    :: ZSPEC(KSPEC2MX,COUNT(KVSET(:)==MYSETV))
REAL(KIND=JPRB), ALLOCATABLE :: ZBUF(:,:,:)
INTEGER(KIND=JPIM) :: IASM0G(0:KSMAX)
INTEGER(KIND=JPIM) :: JM,JN,IFLDR,IFLD,JFLD,ITAG,JNM,ILEN(NPRTRW),JA,ISND(NPRTRV,NPRTRW), JB, IFLDOFF
INTEGER(KIND=JPIM) :: IRCV,ISENDREQ(NPROC), IREQRCV(NPROC), IPROC(NPROC), JMLOC, IFLDBUF, IFLDSPG, IPOSSP
INTEGER(KIND=JPIM) :: ISMAX, ISPEC2, IPOS0,ISENT, INR, IOFFPROC(NPROC+1), IFLDLOC(KFDISTG), IOFF, ILOCFLD(KFDISTG)
INTEGER(KIND=JPIM), POINTER :: ISORT (:)

!     ------------------------------------------------------------------


! Compute help array for distribution

IF (PRESENT (KSORT)) THEN
  ISORT => KSORT
ELSE
  ALLOCATE (ISORT (KFDISTG))
  DO JFLD = 1, KFDISTG
    ISORT (JFLD) = JFLD
  ENDDO
ENDIF

DO JA=1,NPRTRW
  ILEN(JA) = KPOSSP(JA+1)-KPOSSP(JA)
ENDDO
DO JA=1,NPRTRW
  DO JB=1,NPRTRV
    CALL SET2PE(ISND(JB,JA),0,0,JA,JB)
  ENDDO
ENDDO

! Post receive
CALL GSTATS_BARRIER(790)
CALL GSTATS(812,0)
IRCV=0
IOFFPROC(1)=0
IF (ILEN(MYSETW) > 0) THEN
  DO JA=1,NPRTRW
    DO JB=1,NPRTRV
      IF (ISND(JB,JA) /= MYPROC) THEN
        ! count number of fields to receive from each task:
        IFLDR=0
        DO JFLD=1,KFDISTG
          IF (KFROM(JFLD)==ISND(JB,JA)) THEN
            IF (KVSET(JFLD)==MYSETV) THEN
              IFLDR = IFLDR+1
            ENDIF
          ENDIF
        ENDDO
        IF (IFLDR > 0) THEN
          ITAG=MTAGDISTSP+ISND(JB,JA)
          IRCV=IRCV+1
          CALL MPL_RECV(ZSPEC(:,IOFFPROC(IRCV)+1:IOFFPROC(IRCV)+IFLDR),KSOURCE=NPRCIDS(ISND(JB,JA)),KTAG=ITAG,&
           & KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=IREQRCV(IRCV),&
           & CDSTRING='DIST_SPEC_CONTROL:')
          IPROC(IRCV)=ISND(JB,JA)
          IOFFPROC(IRCV+1)=IOFFPROC(IRCV)+IFLDR
        ENDIF
      ENDIF
    ENDDO
  ENDDO
ENDIF
CALL GSTATS(812,1)

!Distribute spectral array

CALL GSTATS(1804,0)

IASM0G(0)=1
DO JM=1,KSMAX
  IASM0G(JM)=IASM0G(JM-1)+KN(JM-1)
ENDDO

CALL GSTATS(1804,1)

ALLOCATE(ZBUF(KSPEC2MX,COUNT(KFROM(:)==MYPROC),NPRTRW))
! The next lines ensure the large array zbuf is allocated right here and not inside an omp loop below,
! where an extra omp synchro might be needed :
IF (SIZE(ZBUF) > 0) THEN
  ZBUF(LBOUND(ZBUF,DIM=1),LBOUND(ZBUF,DIM=2),LBOUND(ZBUF,DIM=3))=HUGE(1._JPRB)
ENDIF

IF (LDIM1_IS_FLD) THEN

  ISENT=0
  DO JA=1,NPRTRW
    IF (ILEN(JA) > 0) THEN
      IFLDOFF=0
      DO JB=1,NPRTRV
        IF (ISND(JB,JA) /= MYPROC) THEN
          ! Locate received fields in source array :
          IFLD=0
          IFLDR=0
          DO JFLD=1,KFDISTG
            IF (KFROM(JFLD)==MYPROC) THEN
              IFLD = IFLD+1
              IF (KVSET(JFLD)==JB) THEN
                IFLDR = IFLDR+1
                IFLDLOC(IFLDR)=IFLD
              ENDIF
            ENDIF
          ENDDO
          IF (IFLDR > 0) THEN
            CALL GSTATS(1644,0)
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JFLD,JMLOC,JM,IFLDBUF,IFLDSPG,IPOSSP)
            DO JFLD=1,IFLDR
              IFLDBUF=IFLDOFF+JFLD
              IFLDSPG=IFLDLOC(JFLD)
              DO JMLOC=1,KUMPP(JA)
                JM=KALLMS(KPTRMS(JA)+JMLOC-1)
                IPOSSP=KDIM0G(JM)-KPOSSP(JA)+1
                ZBUF(IPOSSP:IPOSSP+KN(JM)-1,IFLDBUF,JA) = PSPECG(IFLDSPG,IASM0G(JM):IASM0G(JM)+KN(JM)-1)
              ENDDO
            ENDDO
!$OMP END PARALLEL DO
            CALL GSTATS(1644,1)
            CALL GSTATS(812,0)
            ISENT = ISENT+1
            ITAG  = MTAGDISTSP+MYPROC
            CALL MPL_SEND(ZBUF(:,IFLDOFF+1:IFLDOFF+IFLDR,JA),KDEST=NPRCIDS(ISND(JB,JA)),KTAG=ITAG,&
             & KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=ISENDREQ(ISENT),&
             & CDSTRING='DIST_SPEC_CONTROL:')
            IFLDOFF=IFLDOFF+IFLDR
            CALL GSTATS(812,1)
          ENDIF
        ENDIF
      ENDDO
    ENDIF
  ENDDO

  ! Myself:
  IF (ILEN(MYSETW) > 0) THEN
    ! Locate received fields in target and source arrays:
    IFLD=0
    IFLDR=0
    DO JFLD=1,KFDISTG
      IF (KFROM(JFLD)==MYPROC) THEN
        IFLD = IFLD+1
        IF (KVSET(JFLD)==MYSETV) THEN
          IFLDR = IFLDR+1
          IFLDLOC(IFLDR)=IFLD
        ENDIF
      ENDIF
    ENDDO
    IFLD=0
    IFLDR=0
    DO JFLD=1,KFDISTG
      IF (KVSET(JFLD)==MYSETV) THEN
        IFLD = IFLD+1
        IF (KFROM(JFLD)==MYPROC) THEN
          IFLDR = IFLDR+1
          ILOCFLD(IFLDR)=IFLD
        ENDIF
      ENDIF
    ENDDO
    IF (IFLDR > 0) THEN
      CALL GSTATS(1644,0)
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JFLD,JMLOC,JM,IFLDBUF,IFLDSPG,IPOSSP)
      DO JFLD=1,IFLDR
        IFLDBUF=ISORT(ILOCFLD(JFLD))
        IFLDSPG=IFLDLOC(JFLD)
        DO JMLOC=1,KUMPP(MYSETW)
          JM=KALLMS(KPTRMS(MYSETW)+JMLOC-1)
          IPOSSP=KDIM0G(JM)-KPOSSP(MYSETW)+1
          PSPEC(IFLDBUF,IPOSSP:IPOSSP+KN(JM)-1) = PSPECG(IFLDSPG,IASM0G(JM):IASM0G(JM)+KN(JM)-1)
        ENDDO
      ENDDO
!$OMP END PARALLEL DO
      CALL GSTATS(1644,1)
    ENDIF
  ENDIF

  DO JA=1,IRCV
    CALL GSTATS(812,0)
    CALL MPL_WAITANY(KREQUEST=IREQRCV(1:IRCV),KINDEX=INR,CDSTRING='DIST_SPEC_CTL: WAIT FOR RECV')
    CALL GSTATS(812,1)
    ! Locate received fields in target array :
    IFLD=0
    IFLDR=0
    DO JFLD=1,KFDISTG
      IF (KVSET(JFLD)==MYSETV) THEN
        IFLD=IFLD+1
        IF (KFROM(JFLD)==IPROC(INR)) THEN
          IFLDR = IFLDR+1
          IFLDLOC(IFLDR)=IFLD
        ENDIF
      ENDIF
    ENDDO
    IF (IFLDR > 0) THEN
      CALL GSTATS(1644,0)
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JFLD)
      DO JFLD=1,IFLDR
        PSPEC(ISORT(IFLDLOC(JFLD)),1:KSPEC2) = ZSPEC(1:KSPEC2,IOFFPROC(INR)+JFLD)
      ENDDO
!$OMP END PARALLEL DO
      CALL GSTATS(1644,1)
    ENDIF
  ENDDO

ELSE

  ISENT=0
  DO JA=1,NPRTRW
    IF (ILEN(JA) > 0) THEN
      IFLDOFF=0
      DO JB=1,NPRTRV
        IF (ISND(JB,JA) /= MYPROC) THEN
          ! Locate received fields in source array :
          IFLD=0
          IFLDR=0
          DO JFLD=1,KFDISTG
            IF (KFROM(JFLD)==MYPROC) THEN
              IFLD = IFLD+1
              IF (KVSET(JFLD)==JB) THEN
                IFLDR = IFLDR+1
                IFLDLOC(IFLDR)=IFLD
              ENDIF
            ENDIF
          ENDDO
          IF (IFLDR > 0) THEN
            CALL GSTATS(1644,0)
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JFLD,JMLOC,JM,IFLDBUF,IFLDSPG,IPOSSP)
            DO JFLD=1,IFLDR
              IFLDBUF=IFLDOFF+JFLD
              IFLDSPG=IFLDLOC(JFLD)
              DO JMLOC=1,KUMPP(JA)
                JM=KALLMS(KPTRMS(JA)+JMLOC-1)
                IPOSSP=KDIM0G(JM)-KPOSSP(JA)+1
                ZBUF(IPOSSP:IPOSSP+KN(JM)-1,IFLDBUF,JA) = PSPECG(IASM0G(JM):IASM0G(JM)+KN(JM)-1,IFLDSPG)
              ENDDO
            ENDDO
!$OMP END PARALLEL DO
            CALL GSTATS(1644,1)
            CALL GSTATS(812,0)
            ISENT = ISENT+1
            ITAG  = MTAGDISTSP+MYPROC
            CALL MPL_SEND(ZBUF(:,IFLDOFF+1:IFLDOFF+IFLDR,JA),KDEST=NPRCIDS(ISND(JB,JA)),KTAG=ITAG,&
             & KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=ISENDREQ(ISENT),&
             & CDSTRING='DIST_SPEC_CONTROL:')
            IFLDOFF=IFLDOFF+IFLDR
            CALL GSTATS(812,1)
          ENDIF
        ENDIF
      ENDDO
    ENDIF
  ENDDO
  ! Myself:
  IF (ILEN(MYSETW) > 0) THEN
    ! Locate received fields in target and source arrays:
    IFLD=0
    IFLDR=0
    DO JFLD=1,KFDISTG
      IF (KFROM(JFLD)==MYPROC) THEN
        IFLD = IFLD+1
        IF (KVSET(JFLD)==MYSETV) THEN
          IFLDR = IFLDR+1
          IFLDLOC(IFLDR)=IFLD
        ENDIF
      ENDIF
    ENDDO
    IFLD=0
    IFLDR=0
    DO JFLD=1,KFDISTG
      IF (KVSET(JFLD)==MYSETV) THEN
        IFLD = IFLD+1
        IF (KFROM(JFLD)==MYPROC) THEN
          IFLDR = IFLDR+1
          ILOCFLD(IFLDR)=IFLD
        ENDIF
      ENDIF
    ENDDO
    IF (IFLDR > 0) THEN
      CALL GSTATS(1644,0)
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JFLD,JMLOC,JM,IFLDBUF,IFLDSPG,IPOSSP)
      DO JFLD=1,IFLDR
        IFLDBUF=ISORT(ILOCFLD(JFLD))
        IFLDSPG=IFLDLOC(JFLD)
        DO JMLOC=1,KUMPP(MYSETW)
          JM=KALLMS(KPTRMS(MYSETW)+JMLOC-1)
          IPOSSP=KDIM0G(JM)-KPOSSP(MYSETW)+1
          PSPEC(IPOSSP:IPOSSP+KN(JM)-1,IFLDBUF) = PSPECG(IASM0G(JM):IASM0G(JM)+KN(JM)-1,IFLDSPG)
        ENDDO
      ENDDO
!$OMP END PARALLEL DO
      CALL GSTATS(1644,1)
    ENDIF
  ENDIF
  

  DO JA=1,IRCV
    CALL GSTATS(812,0)
    CALL MPL_WAITANY(KREQUEST=IREQRCV(1:IRCV),KINDEX=INR,CDSTRING='DIST_SPEC_CTL: WAIT FOR RECV')
    CALL GSTATS(812,1)
    ! Locate received fields in target array :
    IFLD=0
    IFLDR=0
    DO JFLD=1,KFDISTG
      IF (KVSET(JFLD)==MYSETV) THEN
        IFLD=IFLD+1
        IF (KFROM(JFLD)==IPROC(INR)) THEN
          IFLDR = IFLDR+1
          IFLDLOC(IFLDR)=IFLD
        ENDIF
      ENDIF
    ENDDO
    IF (IFLDR > 0) THEN
      CALL GSTATS(1644,0)
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JFLD)
      DO JFLD=1,IFLDR
        PSPEC(1:KSPEC2,ISORT(IFLDLOC(JFLD))) = ZSPEC(1:KSPEC2,IOFFPROC(INR)+JFLD)
      ENDDO
!$OMP END PARALLEL DO
      CALL GSTATS(1644,1)
    ENDIF
  ENDDO

ENDIF

CALL GSTATS(812,0)
DO JA=1,ISENT
  CALL MPL_WAIT(KREQUEST=ISENDREQ(JA),CDSTRING='DIST_SPEC_CTL: WAIT FOR SEND')
ENDDO
CALL GSTATS(812,1)

CALL GSTATS_BARRIER2(790)

IF (.NOT. PRESENT (KSORT)) THEN
  DEALLOCATE (ISORT)
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE DIST_SPEC_CONTROL
END MODULE DIST_SPEC_CONTROL_MOD
