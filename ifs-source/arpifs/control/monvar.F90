! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction
! 
! (C) Copyright 1989- Meteo-France.
! 

SUBROUTINE MONVAR(YDRIP,KREFTS,KANATS,KGRATS)

!**** *MONVAR* - Management of the 3-D/4-D VAR  events

!     Purpose.
!     --------
!     Set up the 3-D/4-D VAR events control arrays

!**   Interface.
!     ----------
!        *CALL* *MONVAR(...)

!        Explicit arguments :
!        --------------------
!        KREFTS : ARRAY CONTAINING SIMULATED OBS. EVENTS STEPS
!        KANATS : ARRAY CONTAINING ANALYSIS WRITE-OUT EVENTS
!        KGRATS : ARRAY CONTAINING GRADIENT WRITE-OUT EVENTS

!        Implicit arguments :
!        --------------------
!        None

!     Method.
!     -------
!        See documentation (mon oeil !!)

!     Externals.
!     ----------

!     Reference.
!     ----------
!        Documentation ARPEGE  (mon oeil !!)

!     Author.
!     -------
!      Jean-Noel Thepaut  *ECMWF*
!      somewhat based on MONIO written by P.Courtier.
!      Original : 91-09-05

!     Modifications.
!     --------------
!      Modified by R. Brozkova   : 01-04-27 AVARC events added
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      K. Yessad (July 2014): Move some variables.
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMRIP   , ONLY : TRIP
USE YOMVAR   , ONLY : NREFTS   ,NANATS   ,NGRATS   ,&
 &                    NSIMU    ,NDIAG    ,NFRREF   ,NFRANA   ,NFRGRA   ,NSIM4DL  

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TRIP)        ,INTENT(INOUT) :: YDRIP
INTEGER(KIND=JPIM),INTENT(OUT)   :: KREFTS(0:YDRIP%NSTOP/NFRREF) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KANATS(0:NSIM4DL/NFRANA) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KGRATS(0:NSIM4DL/NFRGRA) 

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: ISS, JS
REAL(KIND=JPRB) :: ZLS, ZTSTEP, ZUNIT
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MONVAR',0,ZHOOK_HANDLE)
ASSOCIATE(NSTART=>YDRIP%NSTART, NSTOP=>YDRIP%NSTOP, TSTEP=>YDRIP%TSTEP)
!     ------------------------------------------------------------------

!      -----------------------------------------------------------

!*       1.    Prepare occurences of simulated obs. events.
!              --------------------------------------------

!  ZUNIT  : NUMBER OF SECONDS IN TIME UNIT
ZUNIT=3600._JPRB
IF(ABS(TSTEP) < TINY(ZLS))THEN
  ZTSTEP=TINY(ZLS)
ELSE
  ZTSTEP=TSTEP
ENDIF

DO JS=NSTART,NSTOP/NFRREF
  KREFTS(JS)=0
ENDDO

!*    1.5   SIMULATED OBS. EVENTS

IF(NREFTS(0) >= 1)THEN
  DO JS=1,NREFTS(0)
    ISS=NREFTS(JS)*NFRREF
    IF(ISS >= NSTART.AND.ISS <= NSTOP)THEN
      KREFTS(NREFTS(JS))=1
    ENDIF
  ENDDO
ELSE
  DO JS=NSTART,NSTOP
    IF(MOD(JS,NFRREF) == 0)THEN
      KREFTS(JS/NFRREF)=1
    ENDIF
  ENDDO
ENDIF

!      -----------------------------------------------------------------

!*       3.    Prepare occurences of ANALYSIS write-outs
!              -----------------------------------------

DO JS=0,NSIM4DL/NFRANA
  KANATS(JS)=0
ENDDO

!*    3.5   ANALYSIS WRITE OUT EVENTS

IF (NANATS(0) >= 1)THEN
  DO JS=1,NANATS(0)
    ISS=NANATS(JS)*NFRANA
    IF (ISS <= NSIMU) KANATS(NANATS(JS))=NDIAG
  ENDDO
ELSE
  DO JS=0,NSIMU
    IF (MOD(JS,NFRANA) == 0) KANATS(JS/NFRANA)=NDIAG
  ENDDO
ENDIF

!      -----------------------------------------------------------------

!*       4.    Prepare occurences of GRADIENT write-outs
!              -----------------------------------------

DO JS=0,NSIM4DL/NFRGRA
  KGRATS(JS)=0
ENDDO

!*    3.5   GRADIENT WRITE OUT EVENTS

IF (NGRATS(0) >= 1)THEN
  DO JS=1,NGRATS(0)
    ISS=NGRATS(JS)*NFRGRA
    IF (ISS <= NSIMU) KGRATS(NGRATS(JS))=NDIAG
  ENDDO
ELSE
  DO JS=0,NSIMU
    IF (MOD(JS,NFRGRA) == 0) KGRATS(JS/NFRGRA)=NDIAG
  ENDDO
ENDIF

!     ------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('MONVAR',1,ZHOOK_HANDLE)
END SUBROUTINE MONVAR

