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

SUBROUTINE ERIEN(CDNAMC,KTYPTR,LDMAP,&
 & KTRONC,KDGL,KNXLON,KNLOPA,PSINLA,&
 & KFLEV,PREF,PVALH,PVBH,&
 & PEPS,KULOUT)  

!**** *ERIEN*  - Read Input ENvironment - LAM case.

!     Purpose.
!     --------
!           It extracts geometry information from ALADIN file.

!**   Interface.
!     ----------
!        *CALL* *ERIEN(...)

!        Explicit arguments :
!        --------------------
!        Input - Output :
!        ----------------------------

!        CDNAMC        ...   Name of the cadre

!        Determination of reference geometry:

!        KTYPTR        ...   Truncation NMSMAX
!        LDMAP         ...   .TRUE. : Map projection calculated by EGGPACK
!                      ...   .FALSE.: Biperiodic experiment file, EGGPACK
!                                     routine not called
!        KTRONC        ...   Truncation NSMAX
!        KDGL          ...   Number of latitudes without poles
!        KNXLON        ...   Max. number of longitudes at a parallel
!        KNLOPA        ...   Limited Area characteristics
!        PSINLA        ...   Horizontal geometry characteristics
!        KFLEV         ...   Number of vertical levels
!        PREF          ...   Reference pressure
!        PVALH         ...   "A" coefficients of vertical system
!        PVBH          ...   "B" coefficients of vertical system

!   -----------------------------------------------------------------
!        Input :
!        -------

!        PEPS          ...   Precision of the tests on real variables
!        KULOUT        ...   Output file unit
!   -----------------------------------------------------------------

!        Implicit arguments :
!        --------------------
!        YOMCST

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!        EGGX_N, some FA.. routines.

!     Reference.
!     ----------
!        ARPEGE/ALADIN Documentation.
!        Document 'Control of coherence between namelist and Arpege File'
!        by R. El Khatib

!     Author.
!     -------
!        R. El Khatib *METEO-FRANCE*

!     Modifications.
!     --------------
!        Original : 91-12-10
!        R. El Khatib 24-Mar-2017 Extracted the KINF==1 case from ECHIEN, to clean up GEOMETRY intents
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMCST_IFSAUX   , ONLY : XRPI , XRA

IMPLICIT NONE

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: JPXGEO=18
INTEGER(KIND=JPIM), PARAMETER :: JPXPAH=8

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM),INTENT(INOUT) :: KFLEV 
CHARACTER(LEN=16),INTENT(IN)     :: CDNAMC
INTEGER(KIND=JPIM),INTENT(INOUT) :: KTYPTR 
LOGICAL           ,INTENT(INOUT) :: LDMAP 
INTEGER(KIND=JPIM),INTENT(INOUT) :: KTRONC 
INTEGER(KIND=JPIM),INTENT(INOUT) :: KDGL 
INTEGER(KIND=JPIM),INTENT(INOUT) :: KNXLON 
INTEGER(KIND=JPIM),INTENT(INOUT) :: KNLOPA(JPXPAH) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSINLA(JPXGEO) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PREF 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVALH(0:KFLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVBH(0:KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEPS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KULOUT

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM), ALLOCATABLE :: INLOPA(:),INOZPA(:)
REAL(KIND=JPRB), ALLOCATABLE :: ZVALH(:),ZVBH(:)
REAL(KIND=JPRB), ALLOCATABLE :: ZSINLA(:)

! Work files for EGGX only
REAL(KIND=JPRB), ALLOCATABLE :: ZGELAM(:,:), ZGELAT(:,:), ZGM(:,:),&
 & ZGENORX(:,:),ZGENORY(:,:)  

INTEGER(KIND=JPIM) :: IERR, IERRA, II, INIVER, INLATI, INXLON, ITRONC, &
 & ITYPTR, JFLEV, JL, JLEV, IROTEQ, ISOTRP, IGIVO, IMAXLEV, IMAXGL, & 
 & IMAXLON, IMAXTRUNC , IBWX, IBWY

LOGICAL :: LLMAP, LLGARD

REAL(KIND=JPRB) :: Z2PI, ZCLOPO, ZCODIL, ZDIFF, ZREF, ZSLAPO, ZSLOPO, ZEPS
REAL(KIND=JPRB) :: ZRPK, ZLON0, ZLAT0, ZLONC, ZLATC, ZDELX, ZDELY, ZELX, ZELY
REAL(KIND=JPRB) :: ZEXWN,ZEYWN, ZLON1, ZLAT1, ZLON2, ZLAT2
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "eggx_n.h"

#include "abor1.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('ERIEN',0,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

!*       0.    Get software limits
!              -------------------

CALL FALIMU(IMAXLEV,IMAXTRUNC,IMAXGL,IMAXLON)
ALLOCATE(INLOPA(IMAXGL))
ALLOCATE(INOZPA(IMAXGL))
ALLOCATE(ZSINLA(IMAXGL))
ALLOCATE(ZVALH(0:IMAXLEV))
ALLOCATE(ZVBH(0:IMAXLEV))

!*       1.    Read file characteristics
!              -------------------------

WRITE(KULOUT,*) ' HAF, HAF : CADRE : ',CDNAMC
LLGARD=.FALSE.
CALL FACIES(CDNAMC,ITYPTR,ZSLAPO,ZCLOPO,ZSLOPO,ZCODIL,ITRONC,&
 & INLATI,INXLON,INLOPA,INOZPA,ZSINLA,INIVER,ZREF,ZVALH,&
 & ZVBH,LLGARD)  

IF (INIVER > KFLEV) THEN
  CALL ABOR1('ERIEN : MAX. NUMBER OF LEVEL IN MODEL TOO SMALL !')
ENDIF

!     ------------------------------------------------------------------

!*       2.    Preliminary test
!              ----------------

IF(ITYPTR > 0) THEN
  WRITE(KULOUT,*) 'YOU ARE USING A FILE ARPEGE ',&
   & 'WHILE THE MODEL EXPECTS A FILE ALADIN'  
  CALL ABOR1('ERIEN: ABOR1 CALLED 2')
ELSE
  ITYPTR = - ITYPTR
  LLMAP=ZCODIL >= 0.0_JPRB
ENDIF

Z2PI = 2.0_JPRB*XRPI

IF(ZSINLA(1) >= 0.0_JPRB) THEN
  ! Echien smells the old EGGX (i.e. the old format of the cadre)
  WRITE(KULOUT,*) ' the cadre >>',CDNAMC,'<< has the old EGGX format'
  WRITE(KULOUT,*) ' => consistency check of the geometry in the cadre&
   & will be more forgiving'

  ZEPS=PEPS*1000._JPRB

  ZRPK=ZSINLA(10)
  ZLON0=ZSINLA(8)
  ZLAT0=ZSINLA(9)
  ZLON1=ZSINLA(4)
  ZLAT1=ZSINLA(5)
  ZLON2=ZSINLA(6)
  ZLAT2=ZSINLA(7)
  ZELX=ZSINLA(13)
  ZELY=ZSINLA(14)
  ZDELX=ZSINLA(15)
  ZDELY=ZSINLA(16)
  ZEXWN=ZSINLA(17)
  ZEYWN=ZSINLA(18)
  ZLONC=ZSINLA(2)
  ZLATC=ZSINLA(3)

  IF (ZRPK < 0.0_JPRB) THEN
    ! latlon case :
    IF (ZLON1 <= ZLON2) THEN 
      ZLONC=MOD(0.5_JPRB*(ZLON1+ZLON2),Z2PI)
    ELSE
      ZLONC=MOD(0.5_JPRB*(ZLON1-Z2PI+ZLON2),Z2PI)
    ENDIF
    ZLATC=0.5_JPRB*(ZLAT1+ZLAT2)
    ZDELX=ZSINLA(15)
    ZDELY=ZSINLA(16)
  ELSEIF(LLMAP) THEN
    ! projection
    ALLOCATE(ZGELAM(INLOPA(3):INLOPA(4),INLOPA(5):INLOPA(6)))
    ALLOCATE(ZGELAT(INLOPA(3):INLOPA(4),INLOPA(5):INLOPA(6)))
    ALLOCATE(ZGM(INLOPA(3):INLOPA(4),INLOPA(5):INLOPA(6)))
    ALLOCATE(ZGENORX(INLOPA(3):INLOPA(4),INLOPA(5):INLOPA(6)))
    ALLOCATE(ZGENORY(INLOPA(3):INLOPA(4),INLOPA(5):INLOPA(6)))
    IROTEQ=INT(ZSINLA(1))
    ISOTRP=INT(ZSINLA(11))
    IGIVO=INT(ZSINLA(12))

    WRITE(KULOUT,*) 'Call EGGX_N by ERIEN'

    CALL EGGX_N(XRPI,XRA,IROTEQ,ZSINLA(2),ZSINLA(3),ZSLAPO,&
     & ZSINLA(4),ZSINLA(5),ZSINLA(6),ZSINLA(7),ZLON0,ZLAT0,&
     & ZSINLA(10),KULOUT,ISOTRP,IGIVO,&
     & ZGELAM,ZGELAT,ZGM,ZGENORX,ZGENORY,&
     & INLOPA(3),INLOPA(4),INLOPA(5),INLOPA(6),&
     & INLOPA(3),INLOPA(4),INLOPA(5),INLOPA(6),&
     & ZDELX,ZDELY,ZLONC,ZLATC)  
    DEALLOCATE(ZGELAM)
    DEALLOCATE(ZGELAT)
    DEALLOCATE(ZGM)
    DEALLOCATE(ZGENORX)
    DEALLOCATE(ZGENORY)
    ZSINLA(1)=REAL(IROTEQ,JPRB)
    ZSINLA(11)=REAL(ISOTRP,JPRB)
    ZSINLA(12)=REAL(IGIVO,JPRB)
  ENDIF

ELSE

  ZEPS=PEPS

  ZRPK =ZSINLA(2)
  ZLON0=ZSINLA(3)
  ZLAT0=ZSINLA(4)
  ZLONC=ZSINLA(5)
  ZLATC=ZSINLA(6)
  ZDELX=ZSINLA(7)
  ZDELY=ZSINLA(8)
  ZELX =ZSINLA(9)
  ZELY =ZSINLA(10)
  ZEXWN=ZSINLA(11)
  ZEYWN=ZSINLA(12)
  ZLON1=ZSINLA(13)
  ZLAT1=ZSINLA(14)
  ZLON2=ZSINLA(15)
  ZLAT2=ZSINLA(16)
  IBWX=INT(ZSINLA(17))
  IBWY=INT(ZSINLA(18))

ENDIF


!*       4.    Read information from file (extracted from ECHIEN)
!              --------------------------------------------------

!*       4.1  Truncation, number of levels, ref. pressure, coef. A, B

  KTYPTR = ITYPTR
  LDMAP  = ZCODIL >= 0.0_JPRB
  KTRONC = ITRONC
  KFLEV  = INIVER
  PREF   = ZREF
  DO JLEV = 0,KFLEV
    PVALH(JLEV) = ZVALH(JLEV)
    PVBH (JLEV) = ZVBH (JLEV)
  ENDDO

!*      4.2  Geometrical characteristics

  KNXLON = INXLON
  KDGL   = INLATI
  IF (ZSINLA(1) >= 0.0_JPRB) THEN
    PSINLA(1) = -1.0_JPRB
  ELSE
    PSINLA(1) = ZSINLA(1)
    PSINLA(17)= ZSINLA(17)
    PSINLA(18)= ZSINLA(18)
  ENDIF
  PSINLA(2) = ZRPK
  PSINLA(3) = ZLON0
  PSINLA(4) = ZLAT0
  PSINLA(5) = ZLONC
  PSINLA(6) = ZLATC
  PSINLA(7) = ZDELX
  PSINLA(8) = ZDELY
  PSINLA(9) = ZELX
  PSINLA(10)= ZELY
  PSINLA(11)= ZEXWN
  PSINLA(12)= ZEYWN
  PSINLA(13)= ZLON1
  PSINLA(14)= ZLAT1
  PSINLA(15)= ZLON2
  PSINLA(16)= ZLAT2
  DO JL= 1, JPXPAH
    KNLOPA(JL) = INLOPA(JL)
  ENDDO

DEALLOCATE(INLOPA)
DEALLOCATE(INOZPA)
DEALLOCATE(ZSINLA)
DEALLOCATE(ZVALH)
DEALLOCATE(ZVBH)

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('ERIEN',1,ZHOOK_HANDLE)
END SUBROUTINE ERIEN
