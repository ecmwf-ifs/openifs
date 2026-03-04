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

SUBROUTINE SUEFRAME(CDMCA,LDMAP,KMSMAX,KSMAX,KDGL,KDLON,&
 & LDGARD,LDMRT,PELON0,PELAT0,PELONC,PELATC,PEDELX,PEDELY,&
 & PELON1,PELAT1,PELON2,PELAT2,&
 & PERPK,PELX,PELY,PEXWN,PEYWN,&
 & KDLUN,KDLUX,KDGUN,KDGUX,KBZONL,KBZONG,KFLEV,&
 & PVP00,PVALH,PVBH,KDOM,KS,KBWX,KBWY)

!**** *SUEFRAME* - Set up the frame of ARPEGE/ALADIN file

!     Purpose.   Initialize the frame of ARPEGE/ALADIN file
!     --------

!**   Interface.
!     ----------
!        *CALL* *SUEFRAME*

!        Explicit arguments :
!        --------------------
!     CDMCA : name of the frame
!     LDMAP : .T./.F. if coordinates/wavelenghts define the domain
!     KMSMAX: truncation in y direction
!     KSMAX : truncation in x direction
!     KDGL  : number of latitude rows in C+I+E
!     KDLON : number of longitude rows in C+I+E
!     LDGARD: .TRUE. if the frame should be kept after closing the last file
!     LDMRT : .TRUE. if Mercator rotated/tilted used
!     PELON0: LA0 geographic longitude of reference for the projection
!     PELAT0: FI0 geographic latitude of reference for the projection
!     PELONC: LAR geographic longitude of the centre of domain
!     PELATC: FIR geographic latitude of the centre of domain
!     PEDELX: grid size along x
!     PEDELY: grid size along y
!     PELON1: LA1G geographic longitude of the SW corner of useful domain
!     PELAT1: FI1G geographic latitude of the SW corner of useful domain
!     PELON2: LA1G geographic longitude of the NE corner of useful domain
!     PELAT2: FI1G geographic latitude of the NE corner of useful domain
!     PERPK : K projection parameter
!     PELX  : wavelenght of the domain in x direction
!     PELY  : wavelenght of the domain in y direction
!     PEXWN : wavenumber in x direction
!     PEYWN : wavenumber in y direction
!     KDLUN : lower first dimension of the domain of interest
!     KDLUX : upper first dimension of the domain of interest
!     KDGUN : lower second dimension of the domain of interest
!     KDGUX : upper second dimension of the domain of interest
!     KBZONL: half-size of I zone in x
!     KBZONG: half-size of I zone in y
!     KFLEV : number of vertical levels
!     PVP00 : ref. pressure
!     PVALH : vertical function A
!     PVBH  : vertical function B
!     KDOM  : how the area C+I and E are computed
!     KS    : unpacked undertruncation
!     KBWX  : portion of scientific E-zone lying inside C+I, x axis
!     KBWY  : portion of scientific E-zone lying inside C+I, y axis

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
!        RYAD EL KHATIB *METEO-FRANCE*

!     Modifications.
!     --------------
!        Original : 94-04-29 from SUEOPH
!        R. El Khatib : 01-08-07 Pruning options
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        Modified : 05-03-01 Ryad El Khatib : Cleanups
!        JD. Gril : 10-02-2005 Mercator Rotated/tilted Arg
!        JD. Gril : 01-02-2006 Test for LMRT=.T. and NCADFORM=0 => abort
!      R. El Khatib 27-Sep-2013 Boyd window in frame
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMOPH0  , ONLY : NCADFORM

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV
CHARACTER(LEN=16) ,INTENT(IN)    :: CDMCA
LOGICAL           ,INTENT(IN)    :: LDMAP 
INTEGER(KIND=JPIM),INTENT(IN)    :: KMSMAX 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSMAX 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDGL 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDLON 
LOGICAL           ,INTENT(IN)    :: LDGARD
LOGICAL           ,INTENT(IN)    :: LDMRT
REAL(KIND=JPRB)   ,INTENT(IN)    :: PELON0 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PELAT0 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PELONC 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PELATC 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEDELX 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEDELY 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PELON1 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PELAT1 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PELON2 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PELAT2 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PERPK 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PELX 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PELY 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEXWN 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEYWN 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDLUN 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDLUX 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDGUN 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDGUX 
INTEGER(KIND=JPIM),INTENT(IN)    :: KBZONL 
INTEGER(KIND=JPIM),INTENT(IN)    :: KBZONG 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVP00
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVALH(0:KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVBH(0:KFLEV)
INTEGER(KIND=JPIM),INTENT(IN)    :: KDOM 
INTEGER(KIND=JPIM),INTENT(IN)    :: KS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KBWX
INTEGER(KIND=JPIM),INTENT(IN)    :: KBWY

REAL(KIND=JPRB) :: ZGEOM(18)
INTEGER(KIND=JPIM) :: IDOM(8),ISUBTR(1)
REAL(KIND=JPRB) :: ZVALH(2),ZVBH(2) 

INTEGER(KIND=JPIM) :: INLATI, INXLON, ITRONC, ITYPTR
INTEGER(KIND=JPIM) :: IMAXLEV, IMAXGL, IMAXLON, IMAXTRUNC

REAL(KIND=JPRB) :: ZCLOPO, ZCODIL, ZSLAPO, ZSLOPO
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

!*       1.   DEFINE THE FRAME

IF (LHOOK) CALL DR_HOOK('SUEFRAME',0,ZHOOK_HANDLE)

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

ZSLAPO = 0.0_JPRB
ZCLOPO = 0._JPRB
ZSLOPO = 0._JPRB
IF (LDMAP) THEN
  ZCODIL = 0._JPRB
ELSE
  ZCODIL = -1.0_JPRB
ENDIF
ITYPTR = -KMSMAX
ITRONC = KSMAX
INLATI = KDGL
INXLON = KDLON

IF (NCADFORM == 0) THEN
  IF (LDMRT) THEN
    CALL ABOR1('SUFRAME : Mercator Rot/Til and OLD CADRE not allowed !')
  ENDIF 
  ! old format of file frame (CADRE)
  ZGEOM (1)  = 0.0_JPRB             ! ex-NROTEQ
  ZGEOM (2)  = 0.0_JPRB             ! ex-ELONR
  ZGEOM (3)  = 0.0_JPRB             ! ex-ELATR
  ZGEOM (4)  = PELON1
  ZGEOM (5)  = PELAT1
  ZGEOM (6)  = PELON2
  ZGEOM (7)  = PELAT2
  ZGEOM (8)  = PELON0
  ZGEOM (9)  = PELAT0
  ZGEOM (10) = PERPK
  ZGEOM (11) = 0.0_JPRB             ! ex-NSOTRP
  ZGEOM (12) = 0.0_JPRB             ! ex-NGIV0
  ZGEOM (13) = PELX
  ZGEOM (14) = PELY
  ZGEOM (15) = PEDELX
  ZGEOM (16) = PEDELY
  ZGEOM (17) = PEXWN
  ZGEOM (18) = PEYWN
ELSE
  ! new format of file frame (CADRE)
  IF (LDMRT) THEN
    ZGEOM (1)  = -2                 ! Mercator Rot/Til
  ELSE
    ZGEOM (1)  = -1                 ! 
  ENDIF
  ZGEOM (2)  = PERPK
  ZGEOM (3)  = PELON0
  ZGEOM (4)  = PELAT0
  ZGEOM (5)  = PELONC
  ZGEOM (6)  = PELATC
  ZGEOM (7)  = PEDELX
  ZGEOM (8)  = PEDELY
  ZGEOM (9)  = PELX
  ZGEOM (10) = PELY
  ZGEOM (11) = PEXWN
  ZGEOM (12) = PEYWN
  ZGEOM (13) = PELON1
  ZGEOM (14) = PELAT1
  ZGEOM (15) = PELON2
  ZGEOM (16) = PELAT2
  ZGEOM (17) = REAL(KBWX,KIND=JPRB)
  ZGEOM (18) = REAL(KBWY,KIND=JPRB)
ENDIF

IDOM  (1)  = KS
IDOM  (2)  = KDOM
IDOM  (3)  = KDLUN
IDOM  (4)  = KDLUX
IDOM  (5)  = KDGUN
IDOM  (6)  = KDGUX
IDOM  (7)  = KBZONL
IDOM  (8)  = KBZONG

IF (KFLEV > 1) THEN
  CALL FACADE(CDMCA,ITYPTR,ZSLAPO,ZCLOPO,ZSLOPO,ZCODIL,ITRONC,&
   & INLATI,INXLON,IDOM,ISUBTR,ZGEOM,KFLEV,PVP00,PVALH,PVBH,LDGARD)  
ELSE
  ZVALH(1) = 0._JPRB
  ZVALH(2) = 0._JPRB
  ZVBH(1) = 0._JPRB
  ZVBH(2) = 1._JPRB
  CALL FACADE(CDMCA,ITYPTR,ZSLAPO,ZCLOPO,ZSLOPO,ZCODIL,ITRONC,&
   & INLATI,INXLON,IDOM,ISUBTR,ZGEOM,KFLEV,PVP00,ZVALH,ZVBH,LDGARD)  
ENDIF

! -------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUEFRAME',1,ZHOOK_HANDLE)
END SUBROUTINE SUEFRAME
