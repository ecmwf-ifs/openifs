! (C) Copyright 2000- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE DIST_GRID_CTL_MOD
CONTAINS
SUBROUTINE DIST_GRID_CTL(PGPG,KFDISTG,KPROMA,KFROM,PGP,KSORT)

!**** *DIST_GRID_CTL* - Distributing global gridpoint array to processors

!     Purpose.
!     --------
!        Routine for distributing gridpoint array

!**   Interface.
!     ----------
!     CALL DIST_GRID_CTL(...)

!     Explicit arguments :
!     --------------------
!     PGPG(:,:)   - Global gridpoint array
!     KFDISTG     - Global number of fields to be distributed
!     KPROMA      - required blocking factor for gridpoint output
!     KFROM(:)    - Processor responsible for distributing each field
!     PGP(:,:,:)  - Local spectral array
!     KSORT(:)    - Add KSORT

!     Externals.  SET2PE - compute "A and B" set from PE
!     ----------  MPL..  - message passing routines

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 2000-04-01
!    P.Marguinaud : 2014-10-10

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE MPL_MODULE  ,ONLY : MPL_RECV, MPL_SEND, MPL_BARRIER, MPL_WAIT,     &
     &                  JP_BLOCKING_STANDARD, JP_NON_BLOCKING_STANDARD

USE TPM_DISTR       ,ONLY : D, MTAGDISTGP, NPRCIDS, MYPROC, NPROC
USE TPM_GEOMETRY    ,ONLY : G

USE SET2PE_MOD      ,ONLY : SET2PE
USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS
USE EQ_REGIONS_MOD  ,ONLY : N_REGIONS, N_REGIONS_NS
!

IMPLICIT NONE

! Declaration of arguments

REAL(KIND=JPRB)    ,OPTIONAL, INTENT(IN)  :: PGPG(:,:)
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KFDISTG
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KPROMA
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KFROM(:)
REAL(KIND=JPRB)             , INTENT(OUT) :: PGP(:,:,:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN), TARGET :: KSORT (:)

! Declaration of local variables

! SS/2018: Removed stack hogs

!REAL(KIND=JPRB)    :: ZDUM(D%NGPTOTMX) -- not used
REAL(KIND=JPRB),ALLOCATABLE :: ZBUF(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE :: ZRCV(:,:) ! (D%NGPTOTMX,KFDISTG)
INTEGER(KIND=JPIM) :: JFLD,JB,JA,IGLOFF,IGL1,IGL2,IOFF,ILAST,ILOFF,ILENR
INTEGER(KIND=JPIM) :: JGL,JLON,ISND,ITAG,J,IRCV
INTEGER(KIND=JPIM) :: JKGLO,IEND,JROF,IBL,JROC
INTEGER(KIND=JPIM) :: ISENDREQ(NPROC,KFDISTG),ILEN(NPROC,KFDISTG), IRECVREQ(KFDISTG)
INTEGER(KIND=JPIM) :: IFROM,IMYFIELDS,IFLD
INTEGER(KIND=JPIM), POINTER :: ISORT (:)
LOGICAL :: LLSAME

!     ------------------------------------------------------------------

IF (PRESENT (KSORT)) THEN
  ISORT => KSORT
ELSE
  ALLOCATE (ISORT (KFDISTG))
  DO JFLD = 1, KFDISTG
    ISORT (JFLD) = JFLD
  ENDDO
ENDIF

! Copy for single PE

IF(NPROC == 1) THEN
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JKGLO,IEND,IOFF,IBL,JFLD,JROF)
  DO JKGLO=1,D%NGPTOT,KPROMA
    IEND = MIN(KPROMA,D%NGPTOT-JKGLO+1)
    IOFF = JKGLO-1
    IBL  = (JKGLO-1)/KPROMA+1
    DO JFLD=1,KFDISTG
      DO JROF=1,IEND
        PGP(JROF,ISORT(JFLD),IBL) = PGPG(IOFF+JROF,JFLD)
      ENDDO
    ENDDO
  ENDDO
!$OMP END PARALLEL DO

ELSEIF(KFDISTG>0) THEN

! test if values in KFROM are all the same
  LLSAME=.TRUE.
  IFROM=KFROM(1)
  DO JFLD=2,KFDISTG
    IF(KFROM(JFLD) /= IFROM) THEN
      LLSAME=.FALSE.
      EXIT
    ENDIF
  ENDDO

  IMYFIELDS = 0
  DO JFLD=1,KFDISTG
    IF(KFROM(JFLD) == MYPROC) THEN
      IMYFIELDS = IMYFIELDS+1
    ENDIF
  ENDDO

  CALL GSTATS(1663,0)
  IF(IMYFIELDS > 0) THEN
    ALLOCATE(ZBUF(D%NGPTOTMX,IMYFIELDS,NPROC))

!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1)&
!$OMP&PRIVATE(JFLD,JA,JB,ISND,IGLOFF,IGL1,IGL2,IOFF,ILAST,J,&
!$OMP&ILOFF,JGL,JLON)
    DO JFLD=1,IMYFIELDS
      DO JA=1,N_REGIONS_NS
        DO JB=1,N_REGIONS(JA)
          CALL SET2PE(ISND,JA,JB,0,0)

          IGLOFF = D%NPTRFRSTLAT(JA)
          IGL1   = D%NFRSTLAT(JA)
          IGL2   = D%NLSTLAT(JA)
          IOFF = 0
          IF(JA > 1) THEN
            IF( D%NLSTLAT(JA-1) == D%NFRSTLAT(JA) )THEN
              ILAST = D%NLSTLAT(JA-1)-1
            ELSE
              ILAST = D%NLSTLAT(JA-1)
            ENDIF
            DO J=D%NFRSTLAT(1),ILAST
              IOFF = IOFF+G%NLOEN(J)
            ENDDO
          ENDIF
          
          ILEN(ISND,JFLD) = 0
          ILOFF = 0
          DO JGL=IGL1,IGL2
            DO JLON=1,D%NONL(IGLOFF+JGL-IGL1,JB)
              ZBUF(ILEN(ISND,JFLD)+JLON,JFLD,ISND) = &
               & PGPG(IOFF+ILOFF+D%NSTA(IGLOFF+JGL-IGL1,JB)+JLON-1,JFLD)
            ENDDO
            ILEN(ISND,JFLD) = ILEN(ISND,JFLD) + D%NONL(IGLOFF+JGL-IGL1,JB)
            ILOFF = ILOFF + G%NLOEN(JGL)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
!$OMP END PARALLEL DO
  ENDIF
  CALL GSTATS(1663,1)
    
  ! Message passing
  CALL GSTATS_BARRIER(791)
  CALL GSTATS(811,0)
  ! Receive

  ALLOCATE(ZRCV(D%NGPTOTMX,KFDISTG))

  IF( LLSAME )THEN
    IRCV = KFROM(1)
    ITAG = MTAGDISTGP
    CALL MPL_RECV(ZRCV,KSOURCE=NPRCIDS(IRCV),KTAG=ITAG,&
     &KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=IRECVREQ(1),CDSTRING='DIST_GRID_CTL:')
  ELSE
    DO JFLD=1,KFDISTG
      IRCV = KFROM(JFLD)
      ITAG = MTAGDISTGP+JFLD
      CALL MPL_RECV(ZRCV(:,JFLD),KSOURCE=NPRCIDS(IRCV),KTAG=ITAG,&
       &KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=IRECVREQ(JFLD),CDSTRING='DIST_GRID_CTL:')
    ENDDO
  ENDIF


  ! Send
  IF( LLSAME )THEN
    IF(KFROM(1) == MYPROC) THEN
      ITAG = MTAGDISTGP
      DO JROC=1,NPROC
        CALL MPL_SEND(ZBUF(:,:,JROC),KDEST=NPRCIDS(JROC),KTAG=ITAG,&
         &KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=ISENDREQ(JROC,1),&
         &CDSTRING='DIST_GRID_CTL')
      ENDDO
    ENDIF
  ELSE
    IFLD = 0
    DO JFLD=1,KFDISTG
      IF(KFROM(JFLD) == MYPROC) THEN
        IFLD = IFLD+1
        ITAG = MTAGDISTGP+JFLD
        DO JROC=1,NPROC
          CALL MPL_SEND(ZBUF(1:ILEN(JROC,IFLD),IFLD,JROC),KDEST=NPRCIDS(JROC),KTAG=ITAG,&
           &KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=ISENDREQ(JROC,JFLD),&
           &CDSTRING='DIST_GRID_CTL')
        ENDDO
      ENDIF
    ENDDO
  ENDIF
      

! Wait for sends and receives to complete

  IF( LLSAME )THEN
    IF(KFROM(1) == MYPROC) THEN
      CALL MPL_WAIT(KREQUEST=ISENDREQ(:,1), &
       & CDSTRING='DIST_GRID_CTL: WAIT 1')
    ENDIF
    CALL MPL_WAIT(KREQUEST=IRECVREQ(1), &
       & CDSTRING='DIST_GRID_CTL: WAIT 2')
  ELSE
    DO JFLD=1,KFDISTG
      IF(KFROM(JFLD) == MYPROC) THEN
        CALL MPL_WAIT(KREQUEST=ISENDREQ(:,JFLD), &
         & CDSTRING='DIST_GRID_CTL: WAIT 3')
      ENDIF
      CALL MPL_WAIT(KREQUEST=IRECVREQ(JFLD), &
         & CDSTRING='DIST_GRID_CTL: WAIT 4')
    ENDDO
  ENDIF

  CALL GSTATS(811,1)
  CALL GSTATS_BARRIER2(791)

  CALL GSTATS(1663,0)
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JKGLO,IEND,IOFF,IBL,JFLD,JROF)
  DO JKGLO=1,D%NGPTOT,KPROMA
    IEND = MIN(KPROMA,D%NGPTOT-JKGLO+1)
    IOFF = JKGLO-1
    IBL  = (JKGLO-1)/KPROMA+1
    DO JFLD=1,KFDISTG
      DO JROF=1,IEND
        PGP(JROF,ISORT(JFLD),IBL) = ZRCV(IOFF+JROF,JFLD)
      ENDDO
    ENDDO
  ENDDO
!$OMP END PARALLEL DO
  CALL GSTATS(1663,1)
  DEALLOCATE(ZRCV)
  !Synchronize processors
  CALL GSTATS(786,0)
  CALL MPL_BARRIER(CDSTRING='DIST_GRID_CTL:')
  CALL GSTATS(786,1)
  IF(ALLOCATED(ZBUF)) DEALLOCATE(ZBUF)
ENDIF

IF (.NOT. PRESENT (KSORT)) THEN
  DEALLOCATE (ISORT)
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE DIST_GRID_CTL
END MODULE DIST_GRID_CTL_MOD




