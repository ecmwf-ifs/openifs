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

SUBROUTINE SUFRAME(CDMCA,PMUCEN,PLOCEN,PSTRET,KSTTYP,KSMAX,KDGL,KDLON,KLOEN,&
 & KMEN,PMU,LDGARD,KFLEV,PVP00,PVALH,PVBH)  

!**** *SUFRAME* - Set up the frame of ARPEGE file

!     Purpose.   Initialize the frame of ARPEGE file
!     --------

!**   Interface.
!     ----------
!        *CALL* *SUFRAME*

!        Explicit arguments :
!        --------------------
!     CDMCA : name of the frame
!     PMUCEN: Sine of the latitude of the pole of interest
!     PLOCEN: longitude of the pole of interest
!     PSTRET: stretching factor
!     KSTTYP: type of transformation
!     KSMAX : truncation
!     KDGL  : number of latitudes
!     KDLON : number of longitude
!     KLOEN : number of points on each latitude row
!     KMEN  : max. zonal wave number for each latitude row
!     PMU   : Sine of latitudes of the Gaussian grid
!     LDGARD: .TRUE. if the frame should be kept after closing the last file
!     KFLEV : number of vertical levels
!     PVP00 : ref. pressure
!     PVALH : vertical function A
!     PVBH  : vertical function B

!        Implicit arguments : None.
!        --------------------

!     Method.
!     -------

!     Externals : FACADE
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      RYAD EL KHATIB *METEO-FRANCE*
!      Original : 94-04-29 from SUOPH

!     Modifications.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      Modified : 05-03-01 Ryad El Khatib : Cleanups
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KDGL
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV
CHARACTER(LEN=16) ,INTENT(IN)    :: CDMCA
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMUCEN 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLOCEN 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSTRET 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTTYP 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSMAX 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLOEN((KDGL+1)/2) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KMEN((KDGL+1)/2) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMU((KDGL+1)/2) 
LOGICAL           ,INTENT(IN)    :: LDGARD
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVP00
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVALH(0:KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVBH(0:KFLEV)

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZSINLA((KDGL+1)/2)
INTEGER(KIND=JPIM) :: INLOPA((KDGL+1)/2),INOZPA((KDGL+1)/2)
REAL(KIND=JPRB) :: ZVALH(2),ZVBH(2)

INTEGER(KIND=JPIM) :: IDGNH, INLATI, INXLON, ITRONC, ITYPTR, JGL
INTEGER(KIND=JPIM) :: IMAXLEV, IMAXGL, IMAXLON, IMAXTRUNC

REAL(KIND=JPRB) :: ZCLOPO, ZCODIL, ZSLAPO, ZSLOPO
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUFRAME',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

CALL FALIMU(IMAXLEV,IMAXTRUNC,IMAXGL,IMAXLON)
IF (IMAXLEV < KFLEV) THEN
  CALL ABOR1('SUFRAME : MAX. NUMBER OF LEVELS IN *FA* TOO SMALL !')
ELSEIF (IMAXGL < KDGL) THEN
  CALL ABOR1('SUFRAME : MAX. NUMBER OF LATITUDE ROWS IN *FA* TOO SMALL !')
ELSEIF (IMAXLON < KDLON) THEN
  CALL ABOR1('SUFRAME : MAX. NUMBER OF LONGITUDE ROWS IN *FA* TOO SMALL !')
ELSEIF (IMAXTRUNC < KSMAX) THEN
  CALL ABOR1('SUFRAME : MAX. TRUNCATION IN *FA* TOO SMALL !')
ENDIF

ZSLAPO = PMUCEN
ZCLOPO = COS(PLOCEN)
ZSLOPO = SIN(PLOCEN)
ZCODIL = PSTRET
ITYPTR = KSTTYP
ITRONC = KSMAX
INLATI = KDGL
INXLON = KDLON
IDGNH  = (KDGL+1)/2
DO JGL=1,IDGNH
  INLOPA(JGL) = KLOEN(JGL)
  INOZPA(JGL) = KMEN(JGL)
  ZSINLA(JGL) = PMU(JGL)
ENDDO
IF (KFLEV > 1) THEN
  CALL FACADE(CDMCA,ITYPTR,ZSLAPO,ZCLOPO,ZSLOPO,ZCODIL,ITRONC,&
   & INLATI,INXLON,INLOPA,INOZPA,ZSINLA,KFLEV,PVP00,&
   & PVALH,PVBH,LDGARD)  
ELSE
  ZVALH(1) = 0._JPRB
  ZVALH(2) = 0._JPRB
  ZVBH(1) = 0._JPRB
  ZVBH(2) = 1._JPRB
  CALL FACADE(CDMCA,ITYPTR,ZSLAPO,ZCLOPO,ZSLOPO,ZCODIL,ITRONC,&
   & INLATI,INXLON,INLOPA,INOZPA,ZSINLA,KFLEV,PVP00,&
   & ZVALH,ZVBH,LDGARD)  
ENDIF

! -------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUFRAME',1,ZHOOK_HANDLE)
END SUBROUTINE SUFRAME

