! (C) Copyright 2000- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE DIST_GRID_32_CTL_MOD
CONTAINS
SUBROUTINE DIST_GRID_32_CTL(PGPG,KFDISTG,KPROMA,KFROM,PGP)

!**** *DIST_GRID_32_CTL* - Distributing global gridpoint array to processors

!     Purpose.
!     --------
!        Routine for distributing gridpoint array

!**   Interface.
!     ----------
!     CALL DIST_GRID_32_CTL(...)

!     Explicit arguments : 
!     -------------------- 
!     PGPG(:,:)   - Global gridpoint array
!     KFDISTG     - Global number of fields to be distributed
!     KPROMA      - required blocking factor for gridpoint output
!     KFROM(:)    - Processor responsible for distributing each field
!     PGP(:,:,:)  - Local spectral array

!     Externals.  SET2PE - compute "A and B" set from PE
!     ----------  MPL..  - message passing routines

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 2000-04-01

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB  ,JPRM
USE MPL_MODULE

USE TPM_DISTR
USE TPM_GEOMETRY

USE SET2PE_MOD
USE ABORT_TRANS_MOD
USE EQ_REGIONS_MOD

IMPLICIT NONE

! Declaration of arguments

REAL(KIND=JPRM)    ,OPTIONAL, INTENT(IN)  :: PGPG(:,:)
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KFDISTG
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KPROMA
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KFROM(:)
REAL(KIND=JPRM)             , INTENT(OUT) :: PGP(:,:,:)

! Declaration of local variables

REAL(KIND=JPRM)    :: ZDUM(D%NGPTOTMX)
REAL(KIND=JPRM),ALLOCATABLE :: ZBUF(:,:,:),ZRCV2(:,:)
REAL(KIND=JPRM)    :: ZRCV(D%NGPTOTMX,KFDISTG)
INTEGER(KIND=JPIM) :: JFLD,JB,JA,IGLOFF,IGL1,IGL2,IOFF,ILAST,ILOFF,ILENR
INTEGER(KIND=JPIM) :: JGL,JLON,ISND,ITAG,J,IRCV
INTEGER(KIND=JPIM) :: JKGLO,IEND,JROF,IBL,JROC
INTEGER(KIND=JPIM) :: ISENDREQ(NPROC,KFDISTG),ILEN(NPROC,KFDISTG)
INTEGER(KIND=JPIM) :: IFROM,IMYFIELDS,IFLD,IFLDSFROM(NPROC)
LOGICAL :: LLSAME

!     ------------------------------------------------------------------

! Copy for single PE

IF(NPROC == 1) THEN
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JKGLO,IEND,IOFF,IBL,JFLD,JROF)
  DO JKGLO=1,D%NGPTOT,KPROMA
    IEND = MIN(KPROMA,D%NGPTOT-JKGLO+1)
    IOFF = JKGLO-1
    IBL  = (JKGLO-1)/KPROMA+1
    DO JFLD=1,KFDISTG
      DO JROF=1,IEND
        PGP(JROF,JFLD,IBL) = PGPG(IOFF+JROF,JFLD) 
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
  ! Send
  IF( LLSAME )THEN
    IF(KFROM(1) == MYPROC) THEN
      ITAG = MTAGDISTGP
      DO JROC=1,NPROC
        CALL MPL_SEND(ZBUF(:,:,JROC),KDEST=NPRCIDS(JROC),KTAG=ITAG,&
         &KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=ISENDREQ(JROC,1),&
         &CDSTRING='DIST_GRID_32_CTL')
      ENDDO
    ENDIF
  ELSE
   IF(IMYFIELDS > 0) THEN
      ITAG = MTAGDISTGP
      DO JROC=1,NPROC
        CALL MPL_SEND(ZBUF(:,:,JROC),KDEST=NPRCIDS(JROC),KTAG=ITAG,&
         &KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=ISENDREQ(JROC,1),&
         &CDSTRING='DIST_GRID_32_CTL')
      ENDDO
    ENDIF
  ENDIF
      
    ! Receive

  IF( LLSAME )THEN
    IRCV = KFROM(1)
    ITAG = MTAGDISTGP
    CALL MPL_RECV(ZRCV,KSOURCE=NPRCIDS(IRCV),KTAG=ITAG,&
     &KMP_TYPE=JP_BLOCKING_STANDARD,KOUNT=ILENR,CDSTRING='DIST_GRID_32_CTL:')
    IF( ILENR /= D%NGPTOTMX*KFDISTG )THEN
      CALL ABORT_TRANS(' DIST_GRID_32_CTL: INVALID RECEIVE MESSAGE LENGTH 1')
    ENDIF
  ELSE
    IFLDSFROM(:)=0
    DO JFLD=1,KFDISTG
      IFLDSFROM(KFROM(JFLD)) = IFLDSFROM(KFROM(JFLD))+1
    ENDDO
    ITAG = MTAGDISTGP
    DO JROC=1,NPROC
      IF(IFLDSFROM(JROC) > 0 ) THEN
        IRCV = JROC
        ALLOCATE(ZRCV2(D%NGPTOTMX,IFLDSFROM(JROC)))
        CALL MPL_RECV(ZRCV2,KSOURCE=NPRCIDS(IRCV),KTAG=ITAG,&
         &KMP_TYPE=JP_BLOCKING_STANDARD,KOUNT=ILENR,CDSTRING='DIST_GRID_32_CTL:')
        IF( ILENR /= D%NGPTOTMX*IFLDSFROM(JROC) )THEN
          CALL ABORT_TRANS(' DIST_GRID_32_CTL: INVALID RECEIVE MESSAGE LENGTH 2')
        ENDIF
        IFLD = 0
        DO JFLD=1,KFDISTG
          IF(KFROM(JFLD) == JROC) THEN
            IFLD = IFLD+1
            ZRCV(1:D%NGPTOT,JFLD) = ZRCV2(1:D%NGPTOT,IFLD)
          ENDIF
        ENDDO
        DEALLOCATE(ZRCV2)
      ENDIF
    ENDDO
  ENDIF

! Wait for send to complete

  IF( LLSAME )THEN
    IF(KFROM(1) == MYPROC) THEN
      CALL MPL_WAIT(KREQUEST=ISENDREQ(:,1), &
       & CDSTRING='DIST_GRID_32_CTL: WAIT 1')
    ENDIF
  ELSEIF(IMYFIELDS > 0) THEN
    CALL MPL_WAIT(KREQUEST=ISENDREQ(:,1), &
     & CDSTRING='DIST_GRID_32_CTL: WAIT 2')
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
        PGP(JROF,JFLD,IBL) = ZRCV(IOFF+JROF,JFLD) 
      ENDDO
    ENDDO
  ENDDO
!$OMP END PARALLEL DO
  CALL GSTATS(1663,1)
  !Synchronize processors
  CALL GSTATS(786,0)
  CALL MPL_BARRIER(CDSTRING='DIST_GRID_32_CTL:')
  CALL GSTATS(786,1)
  IF(ALLOCATED(ZBUF)) DEALLOCATE(ZBUF)
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE DIST_GRID_32_CTL
END MODULE DIST_GRID_32_CTL_MOD




