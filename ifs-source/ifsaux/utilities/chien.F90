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

SUBROUTINE CHIEN(CDNAMC,KTYPTR,PSLAPO,PLOCEN,&
 & PCODIL,KTRONC,KDGL,KNXLON,KNLOPA,KNOZPA,&
 & KHTYP,KFLEV,PREF,PVALH,PVBH,KQUAD,KINF,&
 & KDGSA,KDGEN,PEPS,LDFICP,KULOUT)  

!**** *CHIEN*  - CHeck Input ENvironment

!     Purpose.
!     --------
!           It controls coherence between defined geometry and ARPEGE
!       file. In the case of inconsistency it calls ABORT. This
!       routine could be also used in order to simply get full
!       information from the cadre.

!**   Interface.
!     ----------
!        *CALL* *CHIEN(CDNAMC,KTYPTR,PSLAPO,PLOCEN,
!    &           PCODIL,KTRONC,KDGL,KNXLON,KNLOPA,KNOZPA,
!    &           KHTYP,KFLEV,PREF,PVALH,PVBH,KQUAD,KINF,
!    &           KDGSA,KDGEN,PEPS,LDFICP,KULOUT)

!        Explicit arguments :
!        --------------------

!        Input (Output case is now done in *RIEN*) :
!        ----------------------------

!        CDNAMC        ...   Name of the cadre

!        Determination of reference geometry:

!        KTYPTR        ...   Type of Schmidt transform
!                            1  ===>  Pole is at geog. North Pole
!                            and stretching is equal to 1
!                            2  ===>  General case
!        PSLAPO        ...   Sinus latitude of pole of dilatation
!        PLOCEN        ...   Longitude of pole of dilatation
!        PCODIL        ...   Stretching factor
!        KTRONC        ...   Truncation
!        KDGL          ...   Number of latitudes without poles
!        KNXLON        ...   Max. number of longitudes at a parallel
!        KNLOPA        ...   Number of longitudes at a parallel
!        KNOZPA        ...   Max. wave number at a parallel
!        KHTYP         ...   Type of collocation grid
!                            0 ==>   regular grid
!                            2 ==>   reduced grid towards the poles
!        KFLEV         ...   Number of vertical levels
!        PREF          ...   Reference pressure
!        PVALH         ...   "A" coefficients of vertical system
!        PVBH          ...   "B" coefficients of vertical system
!        KQUAD         ...   Quadrature ( 1 : Gauss ; 2 : Lobatto)
!        LDFICP        ...   .TRUE.  if file contains the poles
!   -----------------------------------------------------------------
!        Input :
!        -------

!        KDGSA         ...   First row of arrays KNLOPA and KNOZPA
!        KDGEN         ...   Last row of arrays KNLOPA and KNOZPA
!        PEPS          ...   Precision of the tests on real variables
!        KINF          ...   Key:
!                            -1 ==> Minimum checks for climate file
!                            and call abort; if O.K. it
!                            gives back LDFICP
!                            0 ==> Check and call abort; if O.K. it
!                            gives back LDFICP
!                            1 ==> Simply gives back full information
!        KULOUT        ...   Output file unit
!   -----------------------------------------------------------------

!   -----------------------------------------------------------------

!        Implicit arguments :
!        --------------------
!        None.

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

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
! (hopefully without too much bugs)

!     Modifications.
!     --------------
!        Original : 91-12-10
!        R El Khatib : 92-02-07
!        R El Khatib : 92-06-01 (option KINF=-1)
!        M Hamrud    : 92-10-01 (NHTYP=2)
!        R El Khatib : 93-03-03 (NHTYP=2 recoded)
!        R El Khatib : 93-05-04 (KNOZPA NOT tested when KINF=-1)
!        R El Khatib : 97-07-22 (Deep cleanup)
!        K. YESSAD   : 98-08-10 removal of LRPOLE option.
!         -> LDPOLE, LLPOLE become .false. and disappear.
!        R El Khatib : 99-09-02  (KNOZPA NOW tested again when KINF=-1)
!        M.Hamrud    : 01-Oct-2003 CY28 Cleaning
!        R El Khatib : 05-03-01  Cleanups
!        O. Marsden  : May 2016  Moved the KINF==1 case to a new routine (RIEN) 
!                                and changed argument intents to IN wherever possible
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)  :: KFLEV 
INTEGER(KIND=JPIM),INTENT(IN)  :: KDGSA 
INTEGER(KIND=JPIM),INTENT(IN)  :: KDGEN 
CHARACTER(LEN=16) ,INTENT(IN)  :: CDNAMC
INTEGER(KIND=JPIM),INTENT(IN)  :: KTYPTR 
REAL(KIND=JPRB)   ,INTENT(IN)  :: PSLAPO 
REAL(KIND=JPRB)   ,INTENT(IN)  :: PLOCEN 
REAL(KIND=JPRB)   ,INTENT(IN)  :: PCODIL 
INTEGER(KIND=JPIM),INTENT(IN)  :: KTRONC 
INTEGER(KIND=JPIM),INTENT(IN)  :: KDGL 
INTEGER(KIND=JPIM),INTENT(IN)  :: KNXLON 
INTEGER(KIND=JPIM),INTENT(IN)  :: KNLOPA(KDGSA:KDGEN) 
INTEGER(KIND=JPIM),INTENT(IN)  :: KNOZPA(KDGSA:KDGEN) 
INTEGER(KIND=JPIM),INTENT(IN)  :: KHTYP 
REAL(KIND=JPRB)   ,INTENT(IN)  :: PREF 
REAL(KIND=JPRB)   ,INTENT(IN)  :: PVALH(0:KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)  :: PVBH(0:KFLEV) 
INTEGER(KIND=JPIM),INTENT(IN)  :: KQUAD 
INTEGER(KIND=JPIM),INTENT(IN)  :: KINF 
REAL(KIND=JPRB)   ,INTENT(IN)  :: PEPS 
LOGICAL           ,INTENT(OUT) :: LDFICP 
INTEGER(KIND=JPIM),INTENT(IN)  :: KULOUT 

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM), ALLOCATABLE :: INLOPA(:),INOZPA(:)
REAL(KIND=JPRB), ALLOCATABLE :: ZVALH(:),ZVBH(:)
REAL(KIND=JPRB), ALLOCATABLE :: ZSINLA(:)
LOGICAL :: LLGARD

INTEGER(KIND=JPIM) :: IDGL, IDGNH, IERR, IERRA, IHTYP, INIVER, INLATI, &
 & INXLON, IQUADF, ISTROW, ITRONC, ITYPTR, JFLEV, JL, JLAT, JLEV, IMAXLEV, &
 & IMAXGL, IMAXLON, IMAXTRUNC 

REAL(KIND=JPRB) :: ZCLOPO, ZCODIL, ZEPS, ZMUNPOL, ZREF, ZSLAPO, ZSLOPO, ZX1
REAL(KIND=JPRB) :: ZX2
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CHIEN',0,ZHOOK_HANDLE)

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
CALL FACIES(CDNAMC,ITYPTR,ZSLAPO,ZCLOPO,ZSLOPO,ZCODIL,ITRONC,&
 & INLATI,INXLON,INLOPA,INOZPA,ZSINLA,INIVER,ZREF,ZVALH,ZVBH,LLGARD)  

IF (KINF == 1) THEN
  IF (INLATI > KDGEN-KDGSA+1) THEN
    CALL ABOR1('CHIEN : MAX. NUMBER OF LATITUDE ROWS IN MODEL TOO SMALL !')
  ENDIF
  IF (INIVER > KFLEV) THEN
    CALL ABOR1('CHIEN : MAX. NUMBER OF LEVEL IN MODEL TOO SMALL !')
  ENDIF
ENDIF

!*       2.    Preliminary tests
!              -----------------

!              Test - type of file

IF(ITYPTR < 0) THEN
  WRITE(KULOUT,*) 'YOU ARE USING A FILE ALADIN ',&
   & 'WHILE THE MODEL EXPECTS A FILE ARPEGE'  
  CALL ABOR1('CHIEN: ABOR1 CALLED 2a')
ENDIF

!              Test - type of collocation grid

IF(INLOPA(1) == INLOPA(INT(INLATI/2))) THEN
  IHTYP = 0
  WRITE(KULOUT,*) 'FILE HAS REGULAR GRID '
ELSE
  IHTYP = 2
  WRITE(KULOUT,*) 'FILE HAS REDUCED GRID '
ENDIF

!              Poles story

ZMUNPOL = 0.9999999999_JPRB
IF(ZSINLA(1) >= ZMUNPOL) THEN
  WRITE(KULOUT,*) ' FILE CONTAINS THE POLES '
  LDFICP = .TRUE.
  IDGL = INLATI - 2
  IDGNH = (IDGL+1)/2
!       The following test ensures that the fields of the file
!       will be read properly :
  IF(INLOPA(1) /= INLOPA(2)) THEN
    WRITE(KULOUT,*) ' FILE ROWS #1 AND #2 DO NOT HAVE THE ',&
     & 'SAME NUMBER OF LONGITUDES'  
    WRITE(KULOUT,*) ' THIS MAKES THE MODEL UNABLE TO READ THE '&
     & ,'FILE PROPERLY'  
    CALL ABOR1('CHIEN: ABOR1 CALLED 2b')
  ENDIF
ELSE
  WRITE(KULOUT,*) ' FILE DOES NOT CONTAINS THE POLES '
  LDFICP = .FALSE.
  IDGL = INLATI
  IDGNH = (IDGL+1)/2
ENDIF

!              Test - Gaussian or Lobatto truncation

IF(LDFICP) THEN
  ZX1 = 1.0_JPRB - ZSINLA(2)
  ZX2 = ZSINLA(2) - ZSINLA(3)
  IF (ZX1 > ZX2) THEN
    IQUADF = 2
    WRITE(KULOUT,*) 'FILE HAS LOBATTO QUADRATURE'
    IF(MOD(INLATI,2) == 0) THEN
      WRITE(KULOUT,*) ' WARNING ! INLATI IS EVEN !'
    ENDIF
  ELSE
    IQUADF = 1
    WRITE(KULOUT,*) 'FILE HAS GAUSSIAN QUADRATURE'
    IF(MOD(INLATI,2) == 1) THEN
      WRITE(KULOUT,*) ' WARNING ! INLATI IS ODD !'
    ENDIF
  ENDIF
ELSE
  IQUADF = 1
  WRITE(KULOUT,*) 'FILE HAS GAUSSIAN QUADRATURE'
  IF(MOD(INLATI,2) == 1) THEN
    WRITE(KULOUT,*) ' WARNING ! INLATI IS ODD !'
  ENDIF
ENDIF

IF ((KINF == 0).OR.(KINF == -1)) THEN

!*       3.    Checklist
!              ---------

  IERR=0

!*      3.1  Spectral dimensions

  IF(ITRONC /= KTRONC) THEN
    WRITE(KULOUT,*) ' TRUNCATION MISMATCH : '&
     & ,'FILE = ',ITRONC, ' ; ARGUMENT = ',KTRONC  
    IERR=1
  ENDIF

!*      3.2  Spectral-related dimensions

  IF(INXLON /= KNXLON) THEN
    WRITE(KULOUT,*) ' MAX. NUMBER OF LONGITUDES MISMATCH : '&
     & ,'FILE = ',INXLON, ' ; ARGUMENT = ',KNXLON  
    IERR=1
  ENDIF
  IF (LDFICP) THEN
    IF(INLATI /= (KDGL+2)) THEN
      WRITE(KULOUT,*) 'MAX. NUMBER OF LATITUDES MISMATCH : '&
       & ,'FILE = ',INLATI, ' INCLUDING POLES ; ARGUMENT = ',KDGL  
      IERR=1
    ELSE
      ISTROW=1
      DO JLAT = ISTROW, (INLATI-1)/2
        IF(KNLOPA(JLAT) /= INLOPA(JLAT+1)) THEN
          WRITE(KULOUT,*) ' NUMBER OF LONGITUDES MISMATCH ON ',&
           & 'ROW ',JLAT,' : ', &
           & 'FILE = ',INLOPA(JLAT+1), ' ; ARGUMENT = ',KNLOPA(JLAT)  
          IERR=1
        ENDIF
        IF(KNOZPA(JLAT) /= INOZPA(JLAT+1)) THEN
          WRITE(KULOUT,*) ' WAVES NUMBER MISMATCH ON ',&
           & 'ROW ',JLAT,' : ', &
           & 'FILE = ',INOZPA(JLAT+1), ' ; ARGUMENT = ',KNOZPA(JLAT)  
          IERR=1
        ENDIF
      ENDDO
    ENDIF
  ELSE
    IF(INLATI /= (KDGL)) THEN
      WRITE(KULOUT,*) 'NUMBER OF LATITUDES MISMATCH : '&
       & ,'FILE = ',INLATI, ' (NO POLES) ; ARGUMENT = ',KDGL  
      IERR=1
    ELSE
      DO JLAT = 1, (INLATI+1)/2
        IF(KNLOPA(JLAT) /= INLOPA(JLAT)) THEN
          WRITE(KULOUT,*) ' NUMBER OF LONGITUDES MISMATCH ON ',&
           & 'ROW ',JLAT,' : ', &
           & 'FILE = ',INLOPA(JLAT), ' ; ARGUMENT = ',KNLOPA(JLAT)  
          IERR=1
        ENDIF
        IF(KNOZPA(JLAT) /= INOZPA(JLAT)) THEN
          WRITE(KULOUT,*) ' WAVES NUMBER MISMATCH ON ',&
           & 'ROW ',JLAT,' : ', &
           & 'FILE = ',INOZPA(JLAT), ' ; ARGUMENT = ',KNOZPA(JLAT)  
          IERR=1
        ENDIF
      ENDDO
    ENDIF
  ENDIF

!*      3.3  Horizontal geometry

  IF(IHTYP /= KHTYP) THEN
    WRITE(KULOUT,*) ' HORIZONTAL GRID MISMATCH : '&
     & ,'FILE = ',IHTYP, ' ; ARGUMENT = ',KHTYP  
    IERR=1
  ENDIF
  IF(ITYPTR /= KTYPTR) THEN
    WRITE(KULOUT,*) ' TRANSFORMATION MISMATCH : '&
     & ,'FILE = ',ITYPTR, ' ; ARGUMENT = ',KTYPTR  
    IERR=1
  ENDIF
  IF(ABS(PSLAPO-ZSLAPO) > PEPS) THEN
    WRITE(KULOUT,*) ' SINE OF LATITUDE OF POLE MISMATCH  : '&
     & ,'FILE = ',ZSLAPO, ' ; ARGUMENT = ',PSLAPO  
    IERR=1
  ENDIF
  IF(ABS(COS(PLOCEN)-ZCLOPO) > PEPS) THEN
    WRITE(KULOUT,*) ' COSINE OF LONGITUDE OF POLE MISMATCH : '&
     & ,'FILE = ',ZCLOPO, ' ; ARGUMENT = ',COS(PLOCEN)  
    IERR=1
  ENDIF
  IF(ABS(SIN(PLOCEN)-ZSLOPO) > PEPS) THEN
    WRITE(KULOUT,*) ' SINE OF LONGITUDE OF POLE MISMATCH : '&
     & ,'FILE = ',ZSLOPO, ' ; ARGUMENT = ',SIN(PLOCEN)  
    IERR=1
  ENDIF
  IF(ABS(ZCODIL-PCODIL) > PEPS) THEN
    WRITE(KULOUT,*) ' STRETCHING MISMATCH : '&
     & ,'FILE = ',ZCODIL, ' ; ARGUMENT = ',PCODIL  
    IERR=1
  ENDIF
  IF(KQUAD /= IQUADF) THEN
    WRITE(KULOUT,*) ' QUADRATURE MISMATCH : ',&
     & 'FILE = ',IQUADF, ' ; ARGUMENT = ',KQUAD  
    IERR=1
  ENDIF

!*      3.4  Vertical levels

  IF (KINF == 0) THEN
    IF(INIVER /= KFLEV) THEN
      WRITE(KULOUT,*) ' NUMBER OF LEVELS MISMATCH : '&
       & ,'FILE = ',INIVER, ' ; ARGUMENT = ',KFLEV  
      IERR=1
    ELSE
      ZEPS=PEPS*10._JPRB*MAX(ZREF,PREF)
      IERRA=0
      DO JFLEV = 0,KFLEV
        IF(ABS(ZVALH(JFLEV)*ZREF-PVALH(JFLEV)*PREF) > ZEPS) THEN
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

  IF(IERR /= 0) THEN
    CALL ABOR1('CHIEN: ABOR1 CALLED 3.4')
  ENDIF

ELSEIF(KINF == 1) THEN

  WRITE(KULOUT,*) 'CHIEN ERROR : CHIEN(..., KINF=1,...) HAS BEEN REPLACED BY RIEN(...)'
  CALL ABOR1('CHIEN: ABOR1 CALLED 4.1')

ELSE
  WRITE(KULOUT,*) 'INTERNAL ERROR : KINF=',KINF
  CALL ABOR1('CHIEN: ABOR1 CALLED 4.2')
ENDIF

DEALLOCATE(INLOPA)
DEALLOCATE(INOZPA)
DEALLOCATE(ZSINLA)
DEALLOCATE(ZVALH)
DEALLOCATE(ZVBH)

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CHIEN',1,ZHOOK_HANDLE)
END SUBROUTINE CHIEN

