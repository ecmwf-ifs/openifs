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

SUBROUTINE SUGRIB(YDDIM,YDEPHY,YDDPHY,YDPHY)

!**** *SUGRIB* - Routine to intitialize parameters for GRIB coding

!     Purpose.
!     --------
!         Set up parameters for GRIB coding including GRIB codes
!     for all parameters.

!**   Interface.
!     ----------
!        *CALL* *SUGRIB*

!        Explicit arguments :  None.
!        --------------------

!     Method.  
!     -------  

!     Externals.
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
!      Modified : 01-05-31 D.Richardson - ECMWF local use 18
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      Modified : 09-03-2005 R. Buizza - Change GRIB header for VAREPS
!      K. Yessad (Jan 2010): remove useless variables.
!      G.Mozdzynski (Feb 2011): OOPS cleaning, use of derived type TCSGLEG
!      M.Leutbecher  28-Mar-2011 correction for VAREPS
!      G. Radnoti    20-Nov-2012 NWINOFF window offset for long window 4DVAR
!      R. El Khatib 27-Apr-2015 systematic inquiry of the surface
!      P. Lopez     27-Apr-2017: added lightning parameters.
!     ------------------------------------------------------------------

USE YOMDIM            , ONLY : TDIM
USE PARKIND1          , ONLY : JPIM, JPRB
USE YOMCT0            , ONLY : NCONF
USE YOMHOOK           , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN            , ONLY : NULNAM, NULOUT
USE YOMPHY            , ONLY : TPHY
USE YOEPHY            , ONLY : TEPHY
USE YOMDPHY           , ONLY : TDPHY
USE YOM_GRIB_CODES    , ONLY : NGRBMX2T, NGRBMN2T, NGRB10FG, NGRBMXTPR, NGRBMNTPR,&
 &                             NGRBSTL1, NGRBSTL2, NGRBSTL3, NGRBSTL4, NGRBSWVL1, NGRBSWVL2, NGRBSWVL3,&
 &                             NGRBSWVL4, NGRBISTL1, NGRBISTL2, NGRBISTL3, NGRBISTL4
USE YOMGRIB           , ONLY : NSTEPLPP, NLOCGRB, NTOTENS, NENSFNB, NWINOFF, NWINOFF_4V,&
 &                             NBITSSH, NBITSSHLNSP, NJDIAG, NJDOMAI, NJITER, NSTREAM, NSYSTEM, NMETHOD,&
 &                             NREFERENCE, NCONSENSUS, NDWD, NMFR, NNCEP, NUKM, NLEG, NSFLEVS,&
 &                             NCYCLE, CFCLASS, CTYPE, NBITSEXPR
USE FITSPECTRUM_MOD   , ONLY : LFORCEZ, LFORCELNSP

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TDIM) , INTENT(IN) :: YDDIM
TYPE(TDPHY) ,INTENT(INOUT):: YDDPHY
TYPE(TEPHY) ,INTENT(INOUT):: YDEPHY
TYPE(TPHY)  ,INTENT(INOUT):: YDPHY

INTEGER(KIND=JPIM) :: JC,JV,II,ITOP,IBOT
INTEGER(KIND=JPIM) :: IGRBCODES(YDDPHY%NCSS), IGRBTMP (4)
REAL(KIND=JPRB) :: ZRDA(YDDPHY%NCSS)
REAL(KIND=JPRB) :: ZBOT,ZTOP
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
#include "surf_inq.h"

#include "namgrib.nam.h"
#include "posnam.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUGRIB',0,ZHOOK_HANDLE)
ASSOCIATE(NSMAX=>YDDIM%NSMAX, &
 & NCSS=>YDDPHY%NCSS, &
 & LRAY=>YDPHY%LRAY, LRAYFM15=>YDPHY%LRAYFM15)
!     ------------------------------------------------------------------

!*       0.    READ NAMELIST
!              -------------

NLOCGRB=1
NENSFNB=0
NTOTENS=0
NWINOFF=0
NWINOFF_4V=0
NSTREAM=0
NSYSTEM=0
NMETHOD=1
NREFERENCE=20000101
NCONSENSUS=0
NDWD=0
NMFR=0
NNCEP=0
NUKM=0
NBITSSH=16
NBITSSHLNSP=16
NBITSEXPR=-1
NJDIAG =99
NJDOMAI=99
NJITER=0
NLEG=1
CFCLASS='rd'
IF (ABS(NCONF)/100 == 8) THEN
  CTYPE='sg'
ELSE
  CTYPE='fc'
ENDIF
LFORCEZ=.TRUE.
LFORCELNSP=.TRUE.

CALL POSNAM(NULNAM,'NAMGRIB')
READ(NULNAM,NAMGRIB)

!*       0.1  CHECK STREAM IS RECOGNIZED
!             --------------------------

WRITE(NULOUT,'(''SUGRIB: NLOCGRB = '',I0,'' NSTREAM = '',I0)') NLOCGRB,NSTREAM

IF(NSTREAM == 0) THEN
  IF(NLOCGRB == 7) THEN !Sensitivity data
    NSTREAM = 1036 ! Sensitivity forecast
  ELSEIF(NLOCGRB == 9) THEN !Singular vectors and ensemble perturbations
    NSTREAM = 1035 ! Ensemble prediction system
  ELSEIF ( NLOCGRB == 15) THEN ! Seasonal forecast data
    NSTREAM = 1090 ! Seasonal forecast
  ELSEIF ( NLOCGRB == 23) THEN !Coupled atmospheric, wave and ocean means (HC)
    NSTREAM = 1200 ! Real-time Monthly forecast
  ELSEIF ( NLOCGRB == 26) THEN !MARS labelling or ensemble forecast data (HC)
    IF(NTOTENS == 0) THEN
      NSTREAM=1024 ! Daily archive hindcast
    ELSE
      NSTREAM=1039 ! Ensemble forecast hindcasts
    ENDIF
  ELSEIF ( NLOCGRB == 27 .OR. NLOCGRB == 30 ) THEN !VarEPS
    IF(NREFERENCE > 0 ) THEN
       NSTREAM=1033 ! Ensemble prediction system hindcast
    ELSE
       NSTREAM=1035 ! Ensemble prediction system
    ENDIF
  ELSEIF ( NLOCGRB == 18) THEN ! Multi-analysis ensemble data
    NSTREAM = 1037 ! Multianalysis ensemble data
  ELSE
    IF(NTOTENS == 0) THEN
      NSTREAM=1025 ! Atmospheric model
    ELSE
      NSTREAM=1035 ! Ensemble prediction system
    ENDIF
  ENDIF
ENDIF
WRITE(NULOUT,'(''SUGRIB: OUTPUT STREAM SET TO '',I6)') NSTREAM

!     ------------------------------------------------------------------


II=0

DO JV = 1, 3
  SELECT CASE(JV)
    CASE(1)
      IGRBTMP = (/NGRBSTL1,NGRBSTL2,NGRBSTL3,NGRBSTL4/)
      IGRBCODES(1:MIN (NCSS, 4)) = IGRBTMP (1:MIN (NCSS, 4))
      CALL SURF_INQ(YDEPHY%YSURF,PRDAT=ZRDA)
    CASE(2)
      IGRBTMP = (/NGRBSWVL1,NGRBSWVL2,NGRBSWVL3,NGRBSWVL4/)
      IGRBCODES(1:MIN (NCSS, 4)) = IGRBTMP (1:MIN (NCSS, 4))
      CALL SURF_INQ(YDEPHY%YSURF,PRDAW=ZRDA)
    CASE(3)
      IGRBTMP = (/NGRBISTL1,NGRBISTL2,NGRBISTL3,NGRBISTL4/)
      IGRBCODES(1:MIN (NCSS, 4)) = IGRBTMP (1:MIN (NCSS, 4))
! Should be:     CALL SURF_INQ(YDEPHY%YSURF,PRDAI=ZRDA)
      CALL SURF_INQ(YDEPHY%YSURF,PRDAW=ZRDA)
    END SELECT
    
    DO JC=1,NCSS
      ZBOT = SUM(ZRDA(1:JC))
      ZTOP = ZBOT-ZRDA(JC)
      ITOP = MIN( ZTOP*100.0_JPRB+0.1_JPRB, 255._JPRB)
      IBOT = MIN( ZBOT*100.0_JPRB+0.1_JPRB, 255._JPRB)
      II=II+1
      NSFLEVS(II,1) = IGRBCODES(JC)
      NSFLEVS(II,2) = IBOT
      NSFLEVS(II,3) = ITOP
  ENDDO
ENDDO



!     ------------------------------------------------------------------


!*       3.    SET TIME OF PREVIOUS P.P. TO 0 FOR GRIB CODES.
!              --------------------------------------------------

NSTEPLPP(:,2) = 0
NSTEPLPP(1,1) = NGRBMX2T
NSTEPLPP(2,1) = NGRBMN2T
NSTEPLPP(3,1) = NGRB10FG
NSTEPLPP(4,1) = NGRBMXTPR
NSTEPLPP(5,1) = NGRBMNTPR

WRITE(UNIT=NULOUT,FMT='('' NCYCLE = '',I6 )') NCYCLE
!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUGRIB',1,ZHOOK_HANDLE)
END SUBROUTINE SUGRIB
