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

SUBROUTINE ECHIEN(CDNAMC,KTYPTR,LDMAP,&
 & KTRONC,KDGL,KNXLON,KNLOPA,PSINLA,&
 & KFLEV,PREF,PVALH,PVBH,KINF,&
 & PEPS,KULOUT)  

!**** *ECHIEN*  - CHeck Input ENvironment:   LAM case

!     Purpose.
!     --------
!           It controls coherence between defined geometry and ALADIN
!       file. In the case of inconsistency it calls ABOR1. This
!       routine could be also used in order to simply get full
!       information from the cadre.

!**   Interface.
!     ----------
!        *CALL* *ECHIEN(...)

!        Explicit arguments :
!        --------------------
!        Input (Output case is now done in *ERIEN*) :
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
!        KINF          ...   Key:
!                            -1 ==> Checks for climate file
!                            and call abort at "warning" mismatch.
!                            0 ==> Checks for all files 
!                            and call abort at "warning" mismatch.
!                            -2 ==> Checks for all files
!                            and call abort at "fatal" mismatch.
!                            1 ==> Simply gives back full information
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
!        Radmila Bubnova *GMAP/COMPAS - stage MICECO*

!     Remark.
!     -------
!******  This routine is a small christmas gift for our friend
!******  Ryad El Khatib **************************************
! (hopefully without too many bugs)

!     Modifications.
!     --------------
!        Original : 91-12-10
!        Modification : 92-02-07  R El Khatib
!        Modification : 92-06-01  R El Khatib (option KINF=-1)
!        Modification : 92-06-21  R Bubnova   (LAM: ECHIEN * )
!        Modification : 94-07-20  R El Khatib (No test on truncation
!           if at least one of the two geometry is fully gridpoint)
!        Modification : 96-04-03  R El Khatib (Test on truncation only
!                       when both geometries are spectral)
!        Modification : 97-07-17  R El Khatib (Remove test on NSOTRP since
!          all actual four corners are controlled)                    
!        Modification : 97-07-22  R El Khatib (Deep cleanup+KINF=-2/-3)
!        Modification : 97-09-17  R El Khatib + J.-F. Estrade (Bugfix on 
!                                 arrays  dimensionnings)
!        Modified 98-08-10 by K. YESSAD: removal of LRPOLE option.
!          - removal of some obsolete comments about poles.
!        Modified 01-04-09 by M. Janousek: New geographic parameters
!        Modified 03-02-27 by S. Petitcol: Correct ZLONC for latlon domains
!        Modified 12-10-2002 by J. Masek : Bugfix for 2D model (LMAP=.F.).
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        Modification : 17-Nov-2004 JD Gril (Mercator Rotated-Tilted)
!      R. El Khatib 27-Sep-2013 Boyd window in frame
!      R. El Khatib 01-Sep-2014 Ref point and Center of domain printed in
!        radians and degrees for an easier debugging of namelists
!      R. El Khatib 24-Mar-2017 Moved the KINF==1 case to a new routine (ERIEN) 
!                                and changed argument intents to IN wherever possible
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMCST_IFSAUX   , ONLY : XRPI , XRA

IMPLICIT NONE

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: JPXGEO=18
INTEGER(KIND=JPIM), PARAMETER :: JPXPAH=8

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV 
CHARACTER(LEN=16),INTENT(IN)     :: CDNAMC
INTEGER(KIND=JPIM),INTENT(IN)    :: KTYPTR 
LOGICAL           ,INTENT(IN)    :: LDMAP 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTRONC 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDGL 
INTEGER(KIND=JPIM),INTENT(IN)    :: KNXLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KNLOPA(JPXPAH) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSINLA(JPXGEO) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PREF 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVALH(0:KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVBH(0:KFLEV) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KINF 
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

IF (LHOOK) CALL DR_HOOK('ECHIEN',0,ZHOOK_HANDLE)

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

IF (KINF == 1) THEN
  IF (INIVER > KFLEV) THEN
    CALL ABOR1('ECHIEN : MAX. NUMBER OF LEVEL IN MODEL TOO SMALL !')
  ENDIF
ENDIF

!     ------------------------------------------------------------------

!*       2.    Preliminary test
!              ----------------

IF(ITYPTR > 0) THEN
  WRITE(KULOUT,*) 'YOU ARE USING A FILE ARPEGE ',&
   & 'WHILE THE MODEL EXPECTS A FILE ALADIN'  
  CALL ABOR1('ECHIEN: ABOR1 CALLED 2')
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

    WRITE(KULOUT,*) 'Call EGGX_N by ECHIEN'

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

IF((KINF == 0).OR.(KINF == -1).OR.(KINF == -2).OR.(KINF == -3)) THEN

!*       3.    Checklist
!              ---------

  IERR=0

!*      3.1  Spectral dimensions

  IF(INLOPA(2) == 1.AND.KNLOPA(2) == 1) THEN
    IF(ITRONC /= KTRONC) THEN
      WRITE(KULOUT,*) ' TRUNCATION NSMAX MISMATCH : '&
       & ,'FILE = ',ITRONC, ' ; ARGUMENT = ',KTRONC  
      IERR=1
    ENDIF
    IF(ITYPTR /= KTYPTR) THEN
      WRITE(KULOUT,*) ' TRUNCATION NMSMAX MISMATCH : '&
       & ,'FILE = ',ITYPTR, ' ; ARGUMENT = ',KTYPTR  
      IERR=1
    ENDIF
  ENDIF
  IF ((INLOPA(2) /= 0.AND.KNLOPA(2) /= 0).OR.&
     & (.NOT.LLMAP.AND..NOT.LDMAP)) THEN  
    IF(INXLON /= KNXLON) THEN
      WRITE(KULOUT,*) ' NUMBER OF LONGITUDES MISMATCH : '&
       & ,'FILE = ',INXLON, ' ; ARGUMENT = ',KNXLON  
      IERR=1
    ENDIF
    IF(INLATI /= KDGL) THEN
      WRITE(KULOUT,*) ' NUMBER OF LATITUDES MISMATCH : '&
       & ,'FILE = ',INLATI, ' ; ARGUMENT = ',KDGL  
      IERR=1
    ENDIF
  ENDIF

!*      3.3  Horizontal geometry

  IF((LDMAP.AND..NOT.LLMAP).OR.(LLMAP.AND..NOT.LDMAP)) THEN

    WRITE(KULOUT,*) ' HORIZONTAL REPRESENTATION LMAP MISMATCH : '&
     & ,'FILE = ',LLMAP, ' ; ARGUMENT = ',LDMAP  
    IERR=1

  ELSEIF(LLMAP.AND.LDMAP) THEN

    IF((ZRPK >= 0.0_JPRB .AND. PSINLA(2) < 0.0_JPRB) .OR.&
       & (ZRPK < 0.0_JPRB .AND. PSINLA(2) >= 0.0_JPRB)) THEN  
      WRITE(KULOUT,*) ' PROJECTION TYPE MISMATCH : '&
       & ,'FILE = ',ZRPK, ' ; ARGUMENT = ',PSINLA(2)  
      IERR=1
    ENDIF

    ZDIFF=ABS(MOD(ZLON0-PSINLA(3),Z2PI))
    IF(ZDIFF > ZEPS.AND.(Z2PI-ZDIFF) > ZEPS) THEN
      WRITE(KULOUT,*) ' REFERENCE LONGITUDE MISMATCH : '&
       & ,'FILE = ',ZLON0,' (',ZLON0*180._JPRB/XRPI,' DEGREES)', &
       & ' ; ARGUMENT = ',PSINLA(3),' (',PSINLA(3)*180._JPRB/XRPI,' DEGREES)'
      IERR=1
    ENDIF

    IF(ABS(ZLAT0-PSINLA(4)) > ZEPS) THEN
      WRITE(KULOUT,*) ' REFERENCE LATITUDE MISMATCH : '&
       & ,'FILE = ',ZLAT0,' (',ZLAT0*180._JPRB/XRPI,' DEGREES)', &
       & ' ; ARGUMENT = ',PSINLA(4),' (',PSINLA(4)*180._JPRB/XRPI,' DEGREES)'  
      IERR=1
    ENDIF

    ZDIFF=ABS(MOD(ZLONC-PSINLA(5),Z2PI))
    IF(ZDIFF > ZEPS.AND.(Z2PI-ZDIFF) > ZEPS) THEN
      WRITE(KULOUT,*) ' DOMAIN CENTRE LONGITUDE MISMATCH : '&
       & ,'FILE = ',ZLONC,' (',ZLONC*180._JPRB/XRPI,' DEGREES)', & 
       & ' ; ARGUMENT = ',PSINLA(5),' (',PSINLA(5)*180._JPRB/XRPI,' DEGREES)'  
      IERR=1
    ENDIF

    IF(ABS(ZLATC-PSINLA(6)) > ZEPS) THEN
      WRITE(KULOUT,*) ' DOMAIN CENTRE LATITUDE MISMATCH : '&
       & ,'FILE = ',ZLATC,' (',ZLATC*180._JPRB/XRPI,' DEGREES)', & 
       & ' ; ARGUMENT = ',PSINLA(6),' (',PSINLA(6)*180._JPRB/XRPI,' DEGREES)'  
      IERR=1
    ENDIF

    IF(ABS(ZDELX-PSINLA(7)) > ZEPS*10000.) THEN
      WRITE(KULOUT,*) ' RESOLUTION IN X MISMATCH : '&
       & ,'FILE = ',ZDELX, ' ; ARGUMENT = ',PSINLA(7)  
      IERR=1
    ENDIF

    IF(ABS(ZDELY-PSINLA(8)) > ZEPS*10000.) THEN
      WRITE(KULOUT,*) ' RESOLUTION IN Y MISMATCH : '&
       & ,'FILE = ',ZDELY, ' ; ARGUMENT = ',PSINLA(8)  
      IERR=1
    ENDIF

    IF(INLOPA(2) == 0) THEN
      IF(KNLOPA(2) /= 0) THEN
        ! Abort when extension zone in argument is NOT null
        IF ((KNLOPA(4)-KNLOPA(3)+1 /= KNXLON).OR.&
         & (KNLOPA(6)-KNLOPA(5)+1 /= KDGL)) THEN  
          IF(KINF == 0.OR.KINF == -1) THEN
            WRITE(KULOUT,*) 'HORIZONTAL DOMAIN INDICATOR (NDOM) ',&
             & 'MISMATCH : ',&
             & 'FILE = ',INLOPA(2), ' (C+I) ; ARGUMENT = ',KNLOPA(2),&
             & ' (C+I+E)'  
            WRITE(KULOUT,*) ' PROPER INITIALIZATION OF (E) '&
             & ,'IS EXPECTED IN THE CALLING SUBROUTINE'  
            IF(KINF == 0) THEN
              II=-2
            ELSE
              II=-3
            ENDIF
            WRITE(KULOUT,*) ' WHEN THIS IS OK, SET KINF=',II,&
             & ' IN THE CALLING SUBROUTINE TO ANIHILATE THIS ABORT'  
            IERR=1
          ENDIF
        ENDIF
      ENDIF
    ELSE
      IF(KNLOPA(2) == 0) THEN
        ! Warning when extension zone in file is NOT null
        IF ((INLOPA(4)-INLOPA(3)+1 /= INXLON).OR.&
         & (INLOPA(6)-INLOPA(5)+1 /= INLATI)) THEN  
          IF(KINF == 0.OR.KINF == -1) THEN
            WRITE(KULOUT,*) 'HORIZONTAL DOMAIN INDICATOR (NDOM) ',&
             & 'MISMATCH : ',&
             & 'FILE = ',INLOPA(2), ' (C+I+E) ; ARGUMENT = ',&
             & KNLOPA(2),' (C+I)'  
            WRITE(KULOUT,*) ' PROPER INITIALIZATION OF (E) '&
             & ,'IS EXPECTED IN THE CALLING SUBROUTINE'  
            IF(KINF == 0) THEN
              II=-2
            ELSE
              II=-3
            ENDIF
            WRITE(KULOUT,*) ' WHEN THIS IS OK, SET KINF=',II,&
             & ' IN THE CALLING SUBROUTINE TO ANIHILATE THIS ABORT'  
            IERR=1
          ENDIF
        ENDIF
      ELSEIF(INLOPA(2) == 1.AND.KNLOPA(2) == -1) THEN
        WRITE(KULOUT,*) ' CAUTION : FILE CONTAINS SPECTRALLY ','FITTED DATA'
      ELSEIF(INLOPA(2) == -1.AND.KNLOPA(2) == 1) THEN
        WRITE(KULOUT,*) ' CAUTION : FILE CONTAINS UNFITTED DATA'
      ENDIF
    ENDIF

    IF(INLOPA(3) /= KNLOPA(3)) THEN
      WRITE(KULOUT,*) ' START INDEX FOR C+I IN X DIRECTION '&
       & ,'(NDLUNG) MISMATCH : '&
       & ,' FILE = ',INLOPA(3), ' ; ARGUMENT = ',KNLOPA(3)  
      IERR=1
    ENDIF

    IF(INLOPA(4) /= KNLOPA(4)) THEN
      WRITE(KULOUT,*) ' END INDEX FOR C+I IN X DIRECTION '&
       & ,'(NDLUXG) MISMATCH : '&
       & ,' FILE = ',INLOPA(4), ' ; ARGUMENT = ',KNLOPA(4)  
      IERR=1
    ENDIF

    IF(INLOPA(5) /= KNLOPA(5)) THEN
      WRITE(KULOUT,*) ' START INDEX FOR C+I IN Y DIRECTION '&
       & ,'(NDGUNG) MISMATCH : '&
       & ,' FILE = ',INLOPA(5), ' ; ARGUMENT = ',KNLOPA(5)  
      IERR=1
    ENDIF

    IF(INLOPA(6) /= KNLOPA(6)) THEN
      WRITE(KULOUT,*) ' END INDEX FOR C+I IN Y DIRECTION '&
       & ,'(NDGUXG) MISMATCH : '&
       & ,' FILE = ',INLOPA(6), ' ; ARGUMENT = ',KNLOPA(6)  
      IERR=1
    ENDIF

    IF(INLOPA(7) /= KNLOPA(7)) THEN
      WRITE(KULOUT,*) 'CAUTION : LENGTH OF I ZONE IN X DIRECTION '&
       & ,'(NBZONL) MISMATCH : '&
       & ,' FILE = ',INLOPA(7), ' ; ARGUMENT = ',KNLOPA(7)  
    ENDIF

    IF(INLOPA(8) /= KNLOPA(8)) THEN
      WRITE(KULOUT,*) 'CAUTION : LENGTH OF I ZONE IN Y DIRECTION '&
       & ,'(NBZONG) MISMATCH : '&
       & ,' FILE = ',INLOPA(8), ' ; ARGUMENT = ',KNLOPA(8)  
    ENDIF

    IF (KINF == 0 .AND. ZSINLA(1) < 0.0_JPRB) THEN
      IF(IBWX < INT(PSINLA(17))) THEN
        WRITE(KULOUT,*) ' PORTION OF SCIENTIFIC E-ZONE LYING INSIDE C+I (X AXIS) TOO BIG : '&
         & ,' FILE = ',IBWX, ' ; ARGUMENT = ',REAL(PSINLA(17),KIND=JPRB)
        IERR=1
      ENDIF
      IF(IBWY < INT(PSINLA(18))) THEN
        WRITE(KULOUT,*) ' PORTION OF SCIENTIFIC E-ZONE LYING INSIDE C+I (Y AXIS) TOO BIG : '&
         & ,' FILE = ',IBWY, ' ; ARGUMENT = ',REAL(PSINLA(18),KIND=JPRB)
        IERR=1
      ENDIF
    ENDIF

  ELSE

    IF(ABS(ZELX-PSINLA(9)) > ZEPS) THEN
      WRITE(KULOUT,*) ' WAVE LENGTH IN X DIRECTION  MISMATCH : '&
       & ,'FILE = ',ZELX, ' ; ARGUMENT = ',PSINLA(9)  
      IERR=1
    ENDIF

    IF(ABS(ZELY-PSINLA(10)) > ZEPS) THEN
      WRITE(KULOUT,*) ' WAVE LENGTH IN Y DIRECTION  MISMATCH : '&
       & ,'FILE = ',ZELY, ' ; ARGUMENT = ',PSINLA(10)  
      IERR=1
    ENDIF

  ENDIF

!*      3.4  Vertical levels

  IF (KINF == 0.OR.KINF == -2) THEN
    IF(INIVER /= KFLEV) THEN
      WRITE(KULOUT,*) ' NUMBER OF LEVELS MISMATCH : '&
       & ,'FILE = ',INIVER, ' ; ARGUMENT = ',KFLEV  
      IERR=1
    ELSE
      IERRA=0
      DO JFLEV = 0,KFLEV
        IF(ABS(ZVALH(JFLEV)*ZREF-PVALH(JFLEV)*PREF) > PEPS) THEN
          WRITE(KULOUT,*) ' VERTICAL FUNCTION *A* MISMATCH ON ',&
           & 'LEVEL ',JFLEV,' : ',&
           & 'FILE = ',ZVALH(JFLEV), ' ; ARGUMENT = ',PVALH(JFLEV)  
          IERRA=1
          IERR=1
        ENDIF
        IF(ABS(ZVBH(JFLEV)-PVBH(JFLEV)) > PEPS) THEN
          WRITE(KULOUT,*) ' VERTICAL FUNCTION *B* MISMATCH ON ',&
           & 'LEVEL ',JFLEV,' : ',&
           & 'FILE = ',ZVBH(JFLEV), ' ; ARGUMENT = ',PVBH(JFLEV)  
          IERR=1
        ENDIF
      ENDDO
      IF (IERRA /= 0) THEN
        WRITE(KULOUT,*) ' REFERENCE PRESSURE : ',&
         & 'FILE = ',ZREF, ' ; ARGUMENT = ',PREF  
      ENDIF
    ENDIF
  ENDIF

!*      3.5  Packing characteristics (fatal ???)

  IF(INLOPA(1) /= KNLOPA(1) .AND. (LLMAP .OR. LDMAP) ) THEN
    WRITE(KULOUT,*) ' PACKING PARAMETER MISMATCH : '&
     & ,'FILE = ',INLOPA(1), ' ; ARGUMENT = ',KNLOPA(1)  
    IERR=1
  ENDIF

  IF(IERR /= 0) THEN
    CALL ABOR1('ECHIEN: ABOR1 CALLED 3.5')
  ENDIF

!     ------------------------------------------------------------------

!*       4.    Bring back information on file
!              ------------------------------

ELSEIF(KINF == 1) THEN

  WRITE(KULOUT,*) 'ECHIEN ERROR : ECHIEN(..., KINF=1,...) HAS BEEN REPLACED BY ERIEN(...)'
  CALL ABOR1('ECHIEN: ABOR1 CALLED 4.1')

ELSE
  WRITE(KULOUT,*) 'INTERNAL ERROR : KINF = ',KINF
  CALL ABOR1('ECHIEN: ABOR1 CALLED 4.2')
ENDIF

DEALLOCATE(INLOPA)
DEALLOCATE(INOZPA)
DEALLOCATE(ZSINLA)
DEALLOCATE(ZVALH)
DEALLOCATE(ZVBH)

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('ECHIEN',1,ZHOOK_HANDLE)
END SUBROUTINE ECHIEN
