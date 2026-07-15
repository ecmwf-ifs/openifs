SUBROUTINE GSTATS_QUERY(KNUM,PTIME)

!**** *GSTATS_QUERY*  - Get current value of gstats timer

!     PURPOSE.
!     --------
!       To query values of gstats timer for use in live output 


!**   INTERFACE.
!     ----------
!       *CALL* *GSTATS(KNUM,PTIME)

!        EXPLICIT ARGUMENTS
!        --------------------
!        KNUM - timing event number (for list of already defined events
!               see routine STATS_OUTPUT)
!        PTIME    - Output current value of timer

!        IMPLICIT ARGUMENTS
!        --------------------
!        Module YOMGSTATS

!     METHOD.
!     -------

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     AUTHOR.
!     -------
!        P. Gillies ECMWF

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 2021-03-03 

!     ------------------------------------------------------------------

USE EC_PARKIND  ,ONLY : JPRD, JPIM
USE YOMHOOK   ,ONLY : LHOOK, DR_HOOK, JPHOOK

USE YOMGSTATS  
USE OML_MOD

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KNUM
REAL(KIND=JPRD),INTENT(OUT) :: PTIME
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('GSTATS_QUERY',0,ZHOOK_HANDLE)

IF(LSTATS) THEN

! only process gstats calls for master thread

  IF(OML_MY_THREAD() <= 1) THEN

    ! Return current total value of specified timer
    IF(NCALLS(KNUM)>1) THEN
      PTIME=TIMESUM(KNUM)
    ELSE
      PTIME=0.0_JPRD
    ENDIF

  ENDIF

ELSE
  PTIME=0.0_JPRD
ENDIF

IF (LHOOK) CALL DR_HOOK('GSTATS_QUERY',1,ZHOOK_HANDLE)

END SUBROUTINE GSTATS_QUERY
