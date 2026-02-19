! (C) Copyright 1989- Meteo-France.

SUBROUTINE DLADDH(YDMDDH,YDPADDH,KLALIST)

!**** *DLADDH* Communicate communication tables for
!              latitude bands to master process

!     Purpose.
!     --------

!**   Interface.
!     ----------
!        *CALL* *dladdh(..)

!        Explicit arguments : KLALIST    - Array of latitude band numbers
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
USE YOMLUN   , ONLY : NULERR
USE YOMMP0   , ONLY : NPRCIDS, MYPROC, NPROC
USE YOMTAG   , ONLY : MTAGDDH2
USE YOMMDDH  , ONLY : TMDDH
USE YOMPADDH , ONLY : TPADDH
USE MPL_MODULE

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TMDDH)       , INTENT(INOUT):: YDMDDH
TYPE(TPADDH)      , INTENT(INOUT):: YDPADDH
INTEGER(KIND=JPIM), INTENT(IN)   :: KLALIST(YDMDDH%NDHKD+1)
INTEGER(KIND=JPIM) :: I, IBFLEN, IMSGLEN, ITAG1, JROC
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('DLADDH',0,ZHOOK_HANDLE)
ASSOCIATE(NDHKD=>YDMDDH%NDHKD, &
 & NGLALIST=>YDPADDH%NGLALIST)
!     ------------------------------------------------------------------

IBFLEN=NDHKD+1

ITAG1=MTAGDDH2

!    master creates global communication table

IF (MYPROC == 1) THEN
  DO I=1,NDHKD+1
    NGLALIST(I,1)=KLALIST(I)
  ENDDO

!    master receives local tables

  DO JROC=2,NPROC
    CALL MPL_RECV(NGLALIST(1:IBFLEN,JROC),KSOURCE=NPRCIDS(JROC),KTAG=ITAG1,&
     & KOUNT=IMSGLEN,CDSTRING='DLADDH:')  
    IF (IMSGLEN /= IBFLEN) THEN
      WRITE(NULERR,*)' MASTERRECV: IMSGLEN,IBFLEN ',IMSGLEN,IBFLEN
      CALL ABOR1('DLADDH: RECEIVED'//' MESSAGE LENGTH OF WRONG SIZE')
    ENDIF

  ENDDO

  DO JROC=2,NPROC
    CALL MPL_SEND(NGLALIST(1:IBFLEN,1:NPROC),KDEST=NPRCIDS(JROC),KTAG=ITAG1,&
     & CDSTRING='DLADDH:')  
  ENDDO

ELSE

!     IF NOT MASTER SEND LOCAL TABLE TO MASTER AND RECEIVE GLOBAL TABLE

  CALL MPL_SEND(KLALIST(1:IBFLEN),KDEST=NPRCIDS(1),KTAG=ITAG1,&
   & CDSTRING='DLADDH:')  
  CALL MPL_RECV(NGLALIST(1:IBFLEN,1:NPROC),KSOURCE=NPRCIDS(1),KTAG=ITAG1,&
   & KOUNT=IMSGLEN,CDSTRING='DLADDH:')  
  IF (IMSGLEN /= IBFLEN*NPROC) THEN
    WRITE(NULERR,*)' MASTERRECV: IMSGLEN,IBFLEN ',IMSGLEN,IBFLEN
    CALL ABOR1('DLADDH: RECEIVED'//' MESSAGE LENGTH OF WRONG SIZE')
  ENDIF
ENDIF

!     ------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('DLADDH',1,ZHOOK_HANDLE)
END SUBROUTINE DLADDH

