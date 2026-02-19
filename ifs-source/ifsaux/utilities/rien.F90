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

SUBROUTINE RIEN(CDNAMC,KTYPTR,PSLAPO,PLOCEN,&
 & PCODIL,KTRONC,KDGL,KNXLON,KNLOPA,KNOZPA,PSINLA,&
 & KHTYP,KFLEV,PREF,PVALH,PVBH,KQUAD,&
 & KDGSA,KDGEN,PEPS,LDFICP,KULOUT)  

!**** *RIEN*  - Read Input ENvironment

!     Purpose.
!     --------
!           It extracts geometry information from ARPEGE file. 


!**   Interface.
!     ----------
!        *CALL* *RIEN(CDNAMC,KTYPTR,PSLAPO,PLOCEN,
!    &           PCODIL,KTRONC,KDGL,KNXLON,KNLOPA,KNOZPA,PSINLA,
!    &           KHTYP,KFLEV,PREF,PVALH,PVBH,KQUAD,
!    &           KDGSA,KDGEN,PEPS,LDFICP,KULOUT)

!        Explicit arguments :
!        --------------------

!        Input-Output:
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
!        KULOUT        ...   Output file unit
!   -----------------------------------------------------------------

!        Output:
!        ----------------------------

!        PSINLA        ...   Sinus of latitudes

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

!     Original CHIEN Author
!     -------
!        Radmila Bubnova *GMAP/COMPAS - stage MICECO*


!     Modifications.
!     --------------
!        Original : 91-12-10
!        O. Marsden : May 2016 Extracted the KINF==1 case from CHIEN, to clean up GEOMETRY intents
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(INOUT) :: KFLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDGSA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDGEN 
CHARACTER(LEN=16) ,INTENT(IN)    :: CDNAMC
INTEGER(KIND=JPIM),INTENT(INOUT) :: KTYPTR 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSLAPO 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLOCEN 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCODIL 
INTEGER(KIND=JPIM),INTENT(INOUT) :: KTRONC 
INTEGER(KIND=JPIM),INTENT(INOUT) :: KDGL 
INTEGER(KIND=JPIM),INTENT(INOUT) :: KNXLON 
INTEGER(KIND=JPIM),INTENT(INOUT) :: KNLOPA(KDGSA:KDGEN) 
INTEGER(KIND=JPIM),INTENT(INOUT) :: KNOZPA(KDGSA:KDGEN) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSINLA(KDGSA:KDGEN) 
INTEGER(KIND=JPIM),INTENT(INOUT) :: KHTYP 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PREF 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVALH(0:KFLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVBH(0:KFLEV) 
INTEGER(KIND=JPIM),INTENT(INOUT) :: KQUAD 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEPS 
LOGICAL           ,INTENT(OUT)   :: LDFICP 
INTEGER(KIND=JPIM),INTENT(IN)    :: KULOUT 

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

IF (LHOOK) CALL DR_HOOK('RIEN',0,ZHOOK_HANDLE)

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

IF (INLATI > KDGEN-KDGSA+1) THEN
  CALL ABOR1('RIEN : MAX. NUMBER OF LATITUDE ROWS IN MODEL TOO SMALL !')
ENDIF
IF (INIVER > KFLEV) THEN
  CALL ABOR1('RIEN : MAX. NUMBER OF LEVEL IN MODEL TOO SMALL !')
ENDIF

!*       2.    Preliminary tests
!              -----------------

!              Test - type of file

IF(ITYPTR < 0) THEN
  WRITE(KULOUT,*) 'YOU ARE USING A FILE ALADIN ',&
   & 'WHILE THE MODEL EXPECTS A FILE ARPEGE'  
  CALL ABOR1('RIEN: ABOR1 CALLED 2a')
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
    CALL ABOR1('RIEN: ABOR1 CALLED 2b')
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



!*       4.    Read information from file (extracted from CHIEN) 

!*      4.1  Pole of dilatation, stretching, truncation, coef. A, B

  KHTYP  = IHTYP
  KQUAD  = IQUADF
  KTYPTR = ITYPTR
  PCODIL = ZCODIL
  PSLAPO = ZSLAPO
  PLOCEN = SIGN(1.0_JPRB,ZSLOPO)*ACOS(ZCLOPO)
  KTRONC = ITRONC
  KFLEV  = INIVER
  PREF   = ZREF
  DO JLEV = 0,KFLEV
    PVALH(JLEV) = ZVALH(JLEV)
    PVBH(JLEV)  = ZVBH(JLEV)
  ENDDO

!*      4.2  Latitudes and longitudes.

  KNXLON = INXLON
  IF(LDFICP) THEN
    KDGL = IDGL
    DO JL= 0, IDGNH
      PSINLA(JL) = ZSINLA(JL+1)
      KNLOPA(JL) = INLOPA(JL+1)
      KNOZPA(JL) = INOZPA(JL+1)
    ENDDO
    DO JL=0, IDGNH
      PSINLA(KDGL+1 - JL) =  - ZSINLA(JL+1)
      KNLOPA(KDGL+1 - JL) = INLOPA(JL+1)
      KNOZPA(KDGL+1 - JL) = INOZPA(JL+1)
    ENDDO
  ELSE
    KDGL = IDGL
    DO JL= 1, IDGNH
      PSINLA(JL) = ZSINLA(JL)
      KNLOPA(JL) = INLOPA(JL)
      KNOZPA(JL) = INOZPA(JL)
    ENDDO
    DO JL=1, IDGNH
      PSINLA(KDGL - JL + 1) =  - ZSINLA(JL)
      KNLOPA(KDGL - JL + 1) = INLOPA(JL)
      KNOZPA(KDGL - JL + 1) = INOZPA(JL)
    ENDDO
    IF (KDGSA < 1) THEN
      PSINLA(0) = 1.0_JPRB
      KNLOPA(0) = INLOPA(1)
      KNOZPA(0) = INOZPA(1)
    ENDIF
    IF (KDGEN > KDGL) THEN
      PSINLA(KDGL+1) =  - 1.0_JPRB
      KNLOPA(KDGL+1) = INLOPA(1)
      KNOZPA(KDGL+1) = INOZPA(1)
    ENDIF
  ENDIF



DEALLOCATE(INLOPA)
DEALLOCATE(INOZPA)
DEALLOCATE(ZSINLA)
DEALLOCATE(ZVALH)
DEALLOCATE(ZVBH)

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('RIEN',1,ZHOOK_HANDLE)
END SUBROUTINE RIEN

