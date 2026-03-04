! (C) Copyright 1989- Meteo-France.

SUBROUTINE DRESDDH(YDDIMV,YDMDDH,YDPADDH,PRESULT,KTID,KMASK)

!**** *DRESDDH* Communicate intermediate results to master

!     Purpose.
!     --------

!**   Interface.
!     ----------
!        *CALL* *DRESDDH*

!        Explicit arguments : presult  - local result(in) on master
!        --------------------            global result(out)
!                             ktid     - flag to distinguish cases
!                                        global (0), bands (1) and
!                                        domaines/points (-1)
!                             KMASK    - domain number

!        Implicit arguments :
!        --------------------

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
!      MPP Group *ECMWF*
!      Original: 95-10-01 

!     Modifications.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      S. Serrar     07-Sep-2006   PRESULT and ZBUFFER extended to include intial
!                                  values of Free Style variables
!      G. Mozdzynski 11-May-2009 scalability optimisation for high task counts
!      T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!      P.Towers      Nov-2015  Remove barrier from end of routine at ECMWF
!      R. El Khatib 12-Aug-2016 optimization
!     ------------------------------------------------------------------

USE YOMDIMV  , ONLY : TDIMV
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMMP0   , ONLY : NPRCIDS, MYPROC, NPROC
USE YOMTAG   , ONLY : MTAGDDH4
USE YOMMDDH  , ONLY : TMDDH
USE YOMPADDH , ONLY : TPADDH
USE YOMCT0   , ONLY : LECMWF
USE MPL_MODULE

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TDIMV)       ,INTENT(IN)    :: YDDIMV
TYPE(TMDDH)       ,INTENT(INOUT) :: YDMDDH
TYPE(TPADDH)      ,INTENT(INOUT) :: YDPADDH
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PRESULT((YDDIMV%NFLEVG+1)*(YDMDDH%NDHCVSU+YDMDDH%NDHVV) &
 &                                          +YDMDDH%NDHCSSU+YDMDDH%NDHVS+YDMDDH%NDHVFS+1)
INTEGER(KIND=JPIM),INTENT(IN) :: KTID 
INTEGER(KIND=JPIM),INTENT(IN) :: KMASK 

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZRESULTS((YDDIMV%NFLEVG+1)*(YDMDDH%NDHCVSU+YDMDDH%NDHVV) &
 & +YDMDDH%NDHCSSU+YDMDDH%NDHVS+YDMDDH%NDHVFS+1)
REAL(KIND=JPRB) :: ZBUFFER((YDDIMV%NFLEVG+1)*(YDMDDH%NDHCVSU+YDMDDH%NDHVV) &
 & +YDMDDH%NDHCSSU+YDMDDH%NDHVS+YDMDDH%NDHVFS+1)

INTEGER(KIND=JPIM) :: I, IMSGLEN, IRCV, ISIZE,&
 & ISND, ITAG1, ITAG11, J, JROC, JI
INTEGER(KIND=JPIM) :: IGROUPSIZE,IMYGROUPBEG,IMYGROUPEND,IGROUPBEG,IGROUPEND

LOGICAL :: LLFLAG(NPROC), LLFLAG_GROUP(NPROC)
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('DRESDDH',0,ZHOOK_HANDLE)
ASSOCIATE(NFLEVG=>YDDIMV%NFLEVG, &
 & NDHCSSU=>YDMDDH%NDHCSSU, NDHCVSU=>YDMDDH%NDHCVSU, NDHVFS=>YDMDDH%NDHVFS, &
 & NDHVS=>YDMDDH%NDHVS, NDHVV=>YDMDDH%NDHVV, &
 & NGLALIST=>YDPADDH%NGLALIST, NGPUMASK=>YDPADDH%NGPUMASK)
!     ------------------------------------------------------------------

IF(NPROC == 1 .AND. LHOOK) CALL DR_HOOK('DRESDDH',1,ZHOOK_HANDLE)
IF(NPROC == 1) RETURN

IGROUPSIZE=SQRT(FLOAT(NPROC))
IMYGROUPBEG=((MYPROC-1)/IGROUPSIZE)*IGROUPSIZE+1
IMYGROUPEND=MIN(NPROC,IMYGROUPBEG+IGROUPSIZE-1)

ISIZE=SIZE(PRESULT)

!     domaines and points

IF( KTID == -1 )THEN

  ITAG1=MTAGDDH4+10*KMASK

  DO JROC=1,NPROC
    IF( NGPUMASK(KMASK,JROC) == 1 )THEN
      LLFLAG(JROC)=.TRUE.
    ELSE
      LLFLAG(JROC)=.FALSE.
    ENDIF
  ENDDO

ELSE

!     global domain or latitude bands

  ITAG1=MTAGDDH4

  DO JROC=1,NPROC
    LLFLAG(JROC)=.FALSE.
    IF(KTID == 0) THEN
      LLFLAG(JROC)=.TRUE.
    ELSE
      DO I=2,NGLALIST(1,JROC)+1
        IF(NGLALIST(I,JROC) == KTID) LLFLAG(JROC)=.TRUE.
      ENDDO
    ENDIF
  ENDDO

ENDIF

!     if group master initialize global result array

IF( MYPROC == IMYGROUPBEG)THEN
  CALL GSTATS(1140,0)
  ZRESULTS(:)=0.0_JPRB
  CALL GSTATS(1140,1)
ENDIF


!    send local results to group master 

IF( LLFLAG(MYPROC) )THEN
  CALL GSTATS(520,0)
  ITAG11=ITAG1+MYPROC
  ISND=NPRCIDS(IMYGROUPBEG)
  CALL MPL_SEND(PRESULT(1:ISIZE),KDEST=NPRCIDS(ISND),KTAG=ITAG11,&
   & CDSTRING='DRESDDH:')  
  CALL GSTATS(520,1)
ENDIF

!    if group master receive local result and derive group global results

IF (MYPROC == IMYGROUPBEG) THEN
  DO JROC = IMYGROUPBEG,IMYGROUPEND
    IF( LLFLAG(JROC) ) THEN
      ITAG11=ITAG1+JROC
      IRCV=NPRCIDS(JROC)
      CALL GSTATS(520,0)
      CALL MPL_RECV(ZBUFFER(1:ISIZE),KSOURCE=NPRCIDS(IRCV),KTAG=ITAG11,&
       & KOUNT=IMSGLEN,CDSTRING='DRESDDH:')  
      CALL GSTATS(520,1)
      IF( IMSGLEN /= ISIZE )THEN
        CALL ABOR1('DRESDDH: ERROR IN RECEIVED MESSAGE LENGTH')
      ENDIF
      CALL GSTATS(1140,0)
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JI)
      DO JI=1,ISIZE
        ZRESULTS(JI)=ZRESULTS(JI)+ZBUFFER(JI)
      ENDDO
!$OMP END PARALLEL DO
      CALL GSTATS(1140,1)
    ENDIF
  ENDDO
ENDIF

!    send group global results to master

LLFLAG_GROUP(:)=.FALSE.
DO J=1,NPROC,IGROUPSIZE
  IGROUPBEG=((J-1)/IGROUPSIZE)*IGROUPSIZE+1
  IGROUPEND=MIN(NPROC,IGROUPBEG+IGROUPSIZE-1)
  DO JROC=IGROUPBEG,IGROUPEND
    IF( LLFLAG(JROC) )THEN
      LLFLAG_GROUP(((JROC-1)/IGROUPSIZE)*IGROUPSIZE+1)=.TRUE.
    ENDIF
  ENDDO
ENDDO

IF (MYPROC == IMYGROUPBEG .AND. LLFLAG_GROUP(MYPROC) ) THEN
  CALL GSTATS(520,0)
  ITAG11=ITAG1+MYPROC
  ISND=NPRCIDS(1)
  CALL MPL_SEND(ZRESULTS(1:ISIZE),KDEST=NPRCIDS(ISND),KTAG=ITAG11,&
   & CDSTRING='DRESDDH:')  
  CALL GSTATS(520,1)
ENDIF

!    if master receive group global results and derive full global results

IF (MYPROC == 1) THEN
  CALL GSTATS(1140,0)
  ZRESULTS(:)=0.0_JPRB
  CALL GSTATS(1140,1)
  DO JROC = 1, NPROC, IGROUPSIZE
    IF( LLFLAG_GROUP(JROC) )THEN
      ITAG11=ITAG1+JROC
      IRCV=NPRCIDS(JROC)
      CALL GSTATS(520,0)
      CALL MPL_RECV(ZBUFFER(1:ISIZE),KSOURCE=NPRCIDS(IRCV),KTAG=ITAG11,&
       & KOUNT=IMSGLEN,CDSTRING='DRESDDH:')  
      CALL GSTATS(520,1)
      IF( IMSGLEN /= ISIZE )THEN
        CALL ABOR1('DRESDDH: ERROR IN RECEIVED MESSAGE LENGTH')
      ENDIF
      CALL GSTATS(1140,0)
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JI)
      DO JI=1,ISIZE
        ZRESULTS(JI)=ZRESULTS(JI)+ZBUFFER(JI)
      ENDDO
!$OMP END PARALLEL DO
      CALL GSTATS(1140,1)
    ENDIF
  ENDDO

!   if master put full global results back in original array

  CALL GSTATS(1140,0)
  PRESULT(:)=ZRESULTS(:)
  CALL GSTATS(1140,1)
ENDIF

IF(.NOT.LECMWF) THEN
  IF( NPROC > 1 )THEN
    CALL GSTATS(794,0)
    CALL MPL_BARRIER(CDSTRING='DRESDDH:')
    CALL GSTATS(794,1)
  ENDIF
ENDIF

!     ------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('DRESDDH',1,ZHOOK_HANDLE)
END SUBROUTINE DRESDDH
