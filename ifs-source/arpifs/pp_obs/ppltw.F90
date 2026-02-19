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

SUBROUTINE PPLTW(YDSTA,KPROMA,KST,KND,KFLEV,PGEO,PT,PXTEMP,PXGEO)

!**** *PPLTW* Compute the geopotential for a given temperature in a vertical profile

!      PURPOSE.
!      --------
!           COMPUTE THE PRESSURE OF A GIVEN WET-BULB TEMPERATURE

!**    INTERFACE.
!      ----------
!           *CALL* PPLTW( ... )

!           EXPLICITE ARGUMENTS.
!           --------------------
!           YDSTA    : replaces YOMSTA%YRSTA                  (INPUT)
!           KPROMA   : horizontal dimension.                  (INPUT)
!           KST      : start of work.                         (INPUT)
!           KND      : depth of work.                         (INPUT)
!           KFLEV    : number of input model levels.          (INPUT)
!           PGEO     : model full level geopotential          (INPUT)
!           PT       : temperature at full levels             (INPUT)
!           PXTEMP   : value of the given temperature         (INPUT)
!           PXGEO   : geopotential of the given temperature  (OUTPUT)

!     METHOD.
!     -------
!        SEE DOCUMENTATION
!        Vertical interpolation : linear in geopotential. 
!        Research is performed below the tropopause to avoid selecting
!        a level inthe statosphere.
!        If no level is found, the outut value is either :
!         - the top of the atmosphere if the temperature profile is always 
!           above PXTEMP,
!         - an extrapolated value (following the standard gradient) if the 
!           temperature profile is always below PXTEMP

!      EXTERNALS.
!      ----------
!           NONE

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS

!      AUTHOR.
!      -------
!           Ryad El Khatib

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 2000-12-04
!    R. El Khatib : 03-04-17 Extrapolation below surface.
!    M.Hamrud      01-Oct-2003 CY28 Cleaning
!    T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!    I. Etchevers : January 2016 : adapted from ppltemp to wet-bulb temperature
!      -----------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMSTA   , ONLY : TSTA,NLEXTRAP, RDTDZ1, RZTROP
USE YOMCST   , ONLY : RG

IMPLICIT NONE

TYPE(TSTA), INTENT(IN)           :: YDSTA
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KND 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEO(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PXTEMP 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PXGEO(KPROMA) 
INTEGER(KIND=JPIM) :: ILEVT(KPROMA,5), ILEVB(KPROMA,5), INLEV(KPROMA)

INTEGER(KIND=JPIM) :: JROF, JLEV, ISTLEV, IENDLEV, IINCLEV, ILEVTOP

REAL(KIND=JPRB) :: ZDT, ZDIFF, ZPROD, ZDP, ZINCT

LOGICAL :: LLDONE(KPROMA), LLEXTRA(KPROMA)
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

!*       1.    COMPUTE LEVELS SURROUNDING THE SEARCHED VALUE
!              ---------------------------------------------

IF (LHOOK) CALL DR_HOOK('PPLTW',0,ZHOOK_HANDLE)
ASSOCIATE(STZ=>YDSTA%STZ)
DO JLEV=2,KFLEV
  IF (STZ(JLEV) < RZTROP) THEN
    ILEVTOP=JLEV-1
    EXIT
  ENDIF
ENDDO

ISTLEV=KFLEV
IENDLEV=ILEVTOP+1
IINCLEV=-1
!WRITE (NULOUT,*) ISTLEV, IENDLEV

ILEVT(KST:KND,1:4)=KFLEV
ILEVB(KST:KND,1:4)=KFLEV-1
INLEV(KST:KND)=0
LLDONE(KST:KND)=.FALSE.
PXGEO(KST:KND)=0.0000001_JPRB

DO JLEV=ISTLEV,IENDLEV,IINCLEV
  DO JROF=KST,KND
    IF (.NOT.LLDONE(JROF)) THEN
      ZPROD=(PT(JROF,JLEV)-PXTEMP)*(PT(JROF,JLEV-1)-PXTEMP)
      ZDIFF=(PT(JROF,JLEV)-PXTEMP)
      IF ((ZDIFF >0.).AND.(ZPROD<0.)) THEN
        INLEV(JROF)=INLEV(JROF)+1
        ILEVT(JROF,INLEV(JROF))=JLEV-1
        ILEVB(JROF,INLEV(JROF))=JLEV
      ENDIF
      IF (INLEV(JROF)>=2) THEN
        LLDONE(JROF)=.TRUE.
      ENDIF
    ENDIF
  ENDDO
ENDDO


!*       2.    INTERPOLATIONS/EXTRAPOLATIONS
!              -----------------------------

DO JROF=KST,KND
  IF (INLEV(JROF)==0) THEN
        IF (MAXVAL(PT(JROF,ILEVTOP:KFLEV)) < PXTEMP) THEN
           PXGEO(JROF)=PGEO(JROF,NLEXTRAP)+ &
             & (PXTEMP-PT(JROF,NLEXTRAP))*RG/RDTDZ1
        ELSEIF  (MINVAL(PT(JROF,ILEVTOP:KFLEV)) > PXTEMP) THEN
           PXGEO(JROF)=PGEO(JROF,ILEVTOP)
        ENDIF
  ELSE
       ZDT=PT(JROF,ILEVT(JROF,INLEV(JROF)))-PT(JROF,ILEVB(JROF,INLEV(JROF)))
       ZDP=PGEO(JROF,ILEVT(JROF,INLEV(JROF)))-PGEO(JROF,ILEVB(JROF,INLEV(JROF)))
       ZINCT=PXTEMP-PT(JROF,ILEVB(JROF,INLEV(JROF)))
       PXGEO(JROF)=PGEO(JROF,ILEVB(JROF,INLEV(JROF)))+ZINCT*ZDP/ZDT
  ENDIF
  IF (PXGEO(JROF)<0.) THEN
     PXGEO(JROF)=0.
  ENDIF
ENDDO

!     -----------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('PPLTW',1,ZHOOK_HANDLE)
END SUBROUTINE PPLTW
