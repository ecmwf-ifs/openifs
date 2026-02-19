! (C) Copyright 1989- Meteo-France.

SUBROUTINE DMADDH(YDMDDH,YDPADDH,KMALIST)

!**** *DMADDH* Communicate communication tables for
!              domain masks to master process

!     Purpose.
!     --------

!**   Interface.
!     ----------
!        *CALL* *dmaddh(...)

!        Explicit arguments : KMALIST    - Array of masks 0,1
!        -------------------- 

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
!      R. El Khatib 06-Jul-2009 Preserve bounds checking
!      T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMMP0   , ONLY : NPRCIDS, MYPROC, NPROC
USE YOMTAG   , ONLY : MTAGDDH3
USE YOMMDDH  , ONLY : TMDDH
USE YOMPADDH , ONLY : TPADDH
USE YOMLUN   , ONLY : NULOUT, NULERR
USE MPL_MODULE

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TMDDH)       ,INTENT(INOUT) :: YDMDDH
TYPE(TPADDH)      ,INTENT(INOUT) :: YDPADDH
INTEGER(KIND=JPIM),INTENT(IN)    :: KMALIST(YDMDDH%NDHNOM)

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: I, IBFLEN, IMSGLEN, ITAG1,&
 & ITYPE, JBOX, JROC  

LOGICAL :: LLFOUND
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('DMADDH',0,ZHOOK_HANDLE)
ASSOCIATE(FNODDH=>YDMDDH%FNODDH, NDHNOM=>YDMDDH%NDHNOM, &
 & NGPUMASK=>YDPADDH%NGPUMASK)
!     ------------------------------------------------------------------

IBFLEN = NDHNOM

WRITE(NULOUT,'(" DMADDH CALLED")')
ITAG1=MTAGDDH3
IF (MYPROC == 1) THEN

!       IF MASTER RECEIVE MASK AND FILL GLOBAL MASK

  DO I=1,NDHNOM
    NGPUMASK(I,1)=KMALIST(I)
  ENDDO

  DO JROC=2,NPROC
    CALL MPL_RECV(NGPUMASK(1:IBFLEN,JROC),KSOURCE=NPRCIDS(JROC),KTAG=ITAG1,&
     & KOUNT=IMSGLEN,CDSTRING='DMADDH:')  
    IF (IMSGLEN /= IBFLEN) THEN
      WRITE(NULERR,*)' MASTERRECV: IMSGLEN,IBFLEN ',IMSGLEN,IBFLEN
      CALL ABOR1('DMADDH: RECEIVED'//' MESSAGE LENGTH OF WRONG SIZE')
    ENDIF

  ENDDO

!        Update NGPUMASK for single points to eliminate duplicates
!        (can occur at poles or equator)

  DO JBOX=1,NDHNOM
    ITYPE=NINT(FNODDH(11,JBOX))
    IF( ITYPE == 4 )THEN
      LLFOUND=.FALSE.
      DO JROC=1,NPROC
        IF( NGPUMASK(JBOX,JROC) == 1 )THEN
          IF( LLFOUND )THEN
            NGPUMASK(JBOX,JROC)=0
          ELSE
            LLFOUND=.TRUE.
          ENDIF
        ENDIF
      ENDDO
    ENDIF
  ENDDO

  DO JROC=2,NPROC
    CALL MPL_SEND(NGPUMASK(1:IBFLEN,1:NPROC),KDEST=NPRCIDS(JROC),KTAG=ITAG1,&
     & CDSTRING='DMADDH:')  
  ENDDO

ELSE

!     IF NOT MASTER SEND LOCAL MASK TO MASTER AND RECEIVE GLOBAL MASK

  CALL MPL_SEND(KMALIST(1:IBFLEN),KDEST=NPRCIDS(1),KTAG=ITAG1,&
   & CDSTRING='DMADDH:')  
  CALL MPL_RECV(NGPUMASK(1:IBFLEN,1:NPROC),KSOURCE=NPRCIDS(1),KTAG=ITAG1,&
   & KOUNT=IMSGLEN,CDSTRING='DMADDH:')  
  IF (IMSGLEN /= IBFLEN*NPROC) THEN
    WRITE(NULERR,*)' MASTERRECV: IMSGLEN,IBFLEN ',IMSGLEN,IBFLEN
    CALL ABOR1('DMADDH: RECEIVED'//' MESSAGE LENGTH OF WRONG SIZE')
  ENDIF
ENDIF

!     ------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('DMADDH',1,ZHOOK_HANDLE)
END SUBROUTINE DMADDH
