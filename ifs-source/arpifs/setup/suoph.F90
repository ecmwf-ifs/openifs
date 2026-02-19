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

SUBROUTINE SUOPH(YDGEOMETRY)

!**** *SUOPH* - Routine to prepare input/output handling

!     Purpose.   To set up common block YOMOPH which contains file-
!     --------   handling parameters.

!**   Interface.
!     ----------
!        *CALL* *SUOPH*

!        Explicit arguments :  None.
!        --------------------

!        Implicit arguments :
!        --------------------

!     Method.
!     -------

!     Externals :
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      Mats Hamrud and Philippe Courtier  *ECMWF*
!      Original : 88-03-01

!     Modifications.
!     --------------
!      M. Janousek: 01-11-26 : NCADFORM - format of FA file frame
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      R.Buizza      23-May-2005 Added CFNVAREPS (input file with VAREPS accum fields)
!      G. Desroziers and K. Yessad (sept 2005):
!       - split option LTRAJHR into LTRAJHR_ALTI and LTRAJHR_SURF.
!       - adapt option LTRAJHR to METEO-FRANCE configurations.
!      JD Gril : 03-feb-2006 : default NCADFORM=1
!      R. El Khatib : 23-Oct-2008 : CFPATH now a directory ending with a '/'
!        and used also for CFNHWF + merge with sueoph + cleanings + CEFNLSH in namelist
!      E.Holm        13-Nov-2008 Added CFNBGHRSH/GG, background 
!                                at outer loop resolution.
!      K. Yessad (Jan 2010): remove useless variables.
!      G.Mozdzynski (Feb 2011): OOPS cleaning, use of derived type TCSGLEG
!       R. El Khatib 10-Aug-2011 NIOLEVG management
!      P.Marguinaud : 26-Apr-2012 : Handle new parameter NTIMEFMT
!      P. Bechtold 14/05/2012 replace 86400 by INT(RDAY)
!      R. El Khatib 27-Sep-2013 Boyd window in frame
!      T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!      K. Yessad (July 2014): Move some variables.
!      M.Hamrud  (Jan 2015) ! Remove FDB stuff from here
!      R. El Khatib 08-Dec-2015 Interoperability GRIB2 vs FA
!      R. El Khatib 17-Aug-2016 Cleaning
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCT0   , ONLY : LARPEGEF, LELAM
USE YOMARG   , ONLY : NGRIBFILE
USE YOMOPH0  , ONLY : CNMCA, YMDLOPH
USE YOMVERT  , ONLY : VP00
USE FA_MOD, ONLY : JPPRCM

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY), INTENT(INOUT) :: YDGEOMETRY
INTEGER(KIND=JPIM) :: IDMOPL, INBCSP, INBPDG, INGRIB, IPUILA, ISTRON, IREP

LOGICAL :: LLGARD=.TRUE.
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
#ifdef __INTEL_COMPILER
INTRINSIC :: SHAPE ! Fails with Intel compiler as it thinks we use ATLAS shape function
#endif


!     ------------------------------------------------------------------

#include "suframe.intfb.h"
#include "sueframe.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUOPH',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, &
 & YDCSGLEG=>YDGEOMETRY%YRCSGLEG, YDEDIM=>YDGEOMETRY%YREDIM, YDEGEO=>YDGEOMETRY%YREGEO,&
 & YVABIO=>YDGEOMETRY%YVABIO)
ASSOCIATE(NBZONG=>YDEDIM%NBZONG, NEDOM=>YDEDIM%NEDOM, NBZONL=>YDEDIM%NBZONL, &
 & NBIPINCIY=>YDEDIM%NBIPINCIY, NBIPINCIX=>YDEDIM%NBIPINCIX, &
 & NDGLG=>YDDIM%NDGLG, NSMAX=>YDDIM%NSMAX, NDLUXG=>YDDIM%NDLUXG, &
 & NDGUNG=>YDDIM%NDGUNG, NMSMAX=>YDDIM%NMSMAX, NDGUXG=>YDDIM%NDGUXG, &
 & NDLUNG=>YDDIM%NDLUNG, NDLON=>YDDIM%NDLON, NSPEC2G=>YDDIM%NSPEC2G, &
 & RMUCEN=>YDGEM%RMUCEN, RLOCEN=>YDGEM%RLOCEN, NLOENG=>YDGEM%NLOENG, &
 & NSTTYP=>YDGEM%NSTTYP, RSTRET=>YDGEM%RSTRET, NMENG=>YDGEM%NMENG, NGPTOTG=>YDGEM%NGPTOTG, &
 & ELON1=>YDEGEO%ELON1, ELON0=>YDEGEO%ELON0, ELATC=>YDEGEO%ELATC, &
 & ELON2=>YDEGEO%ELON2, EDELY=>YDEGEO%EDELY, LMRT=>YDEGEO%LMRT, &
 & ELONC=>YDEGEO%ELONC, EDELX=>YDEGEO%EDELX, ELAT1=>YDEGEO%ELAT1, &
 & ELAT0=>YDEGEO%ELAT0, ELX=>YDEGEO%ELX, ELAT2=>YDEGEO%ELAT2, EYWN=>YDEGEO%EYWN, &
 & LMAP=>YDEGEO%LMAP, EXWN=>YDEGEO%EXWN, ERPK=>YDEGEO%ERPK, ELY=>YDEGEO%ELY, &
 & NIOLEVG=>YDDIMV%NIOLEVG)

!     ------------------------------------------------------------------

!*       1.    INITIALISE STANDARD FRAME OF *FA*
!              ---------------------------------

IF (NGRIBFILE==0.OR.LARPEGEF) THEN

  IF (LELAM) THEN

    CNMCA='CADRE.STANDARD.E'
    ISTRON=MAX(0,MIN(11,NMSMAX,NSMAX)-1)

    CALL SUEFRAME(CNMCA,LMAP,NMSMAX,NSMAX,NDGLG,NDLON,LLGARD,LMRT,ELON0,ELAT0,&
     & ELONC,ELATC,EDELX,EDELY,ELON1,ELAT1,ELON2,ELAT2,ERPK,ELX,ELY,EXWN,EYWN,&
     & NDLUNG,NDLUXG,NDGUNG,NDGUXG,NBZONL,NBZONG,NIOLEVG,VP00,YVABIO%VALH,YVABIO%VBH, &
     & NEDOM,ISTRON,NBIPINCIX,NBIPINCIY)

  ELSE

!*         3.1  MODIFICATION OF IMPLICITS OPTIONS

    IF (NSMAX <= 11) THEN
      INGRIB=2
      INBPDG=24
      INBCSP=24
      ISTRON=NSMAX-1
      IPUILA=0
      IDMOPL=10
      CALL FAGIOT(INGRIB,INBPDG,INBCSP,ISTRON,IPUILA,IDMOPL)
    ENDIF

!*       3.2  DEFINE STANDARD FRAME

    CNMCA = 'CADRE.STANDARD  '
    CALL SUFRAME(CNMCA,RMUCEN,RLOCEN,RSTRET,NSTTYP,NSMAX,NDGLG,NDLON,&
     & NLOENG(1:),NMENG(1:),YDCSGLEG%RMU(1:),LLGARD,NIOLEVG,VP00,YVABIO%VALH,YVABIO%VBH)  

  ENDIF

  CALL FASGRA (IREP, CNMCA , YMDLOPH(1)%NHEADMAX)
  YMDLOPH(1)%CFPCA = CNMCA
  YMDLOPH(1)%NHEADMAX = JPPRCM*YMDLOPH(1)%NHEADMAX
  YMDLOPH(1)%NGPSIZPK = NGPTOTG  + YMDLOPH(1)%NHEADMAX ! See FACOND for justification
  YMDLOPH(1)%NSPSIZPK = NSPEC2G  + YMDLOPH(1)%NHEADMAX ! See FACOND for justification

ENDIF

! -------------------------------------------------------------------


END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUOPH',1,ZHOOK_HANDLE)
END SUBROUTINE SUOPH
