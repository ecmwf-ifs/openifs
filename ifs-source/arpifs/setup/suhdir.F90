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

SUBROUTINE SUHDIR(YDGEOMETRY,YDDYN,YDDYNA)

!**** *SUHDIR*  - horizontal diffusion coefficients
!                 (computation from NAMELIST global parameters)

!     Purpose.
!     --------
!           Computes value of horizontal diffusion
!           coefficients (HDIR*) from NAMELIST
!           parameters (LNEWHD,RRDXTAU,RDAMPDIV) and
!           configuration parameters (NSMAX, NDLON,...)
!           Remark for LELAM: the mapping factor is assumed to be
!            very close to 1 in the projection.

!**   Interface.
!     ----------
!        *CALL* *SUHDIR( ... )

!        Explicit arguments :
!        --------------------

!        Implicit arguments :
!        --------------------
!        none.

!     Method.
!     -------
!        See documentation about horizontal diffusion scheme.

!     Externals.
!     ----------
!        none.
!        Called by SUDYN.

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        K. YESSAD (METEO-FRANCE/CNRM/GMAP).
!        Original: APRIL 1996.

!     Modifications.
!     --------------
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        F. Vana and K. YESSAD (Dec 2003): rewriting the SLHD diffusion.
!        F. Vana (29-Nov-2004): more flexibility in new HD set-up for M-F
!        R. El Khatib 05-02-15 Fix on LREPHD
!        F. Vana (15-Feb-2005): Split of the global LSLHD key.
!        F. Vana (September 2005): more freedom to SLHD sup. diff.
!        K. Yessad (07-Mar-2007): remove LREPHD.
!        O.Spaniel, F.Vana (04-Nov-2007): bug fixed in SLHD.
!        F. Vana 13-Jan-2009 : No supp. diffusion when LSLHD_STATIC
!        N.Semane+P.Bechtold   04-10-2012 add RPLRADI for small planet
!        F. Vana 15-Mar-2013 : Simplified sup. diff. setup when LECMWF=.T.
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPRB, JPIB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCT0   , ONLY : LELAM, LECMWF
USE YOMDYNA  , ONLY : TDYNA
USE YOMDYN   , ONLY : TDYN
USE YOMCST   , ONLY : RA, RPI
USE YOMDYNCORE,ONLY : RPLRADI

!     ------------------------------------------------------------------

IMPLICIT NONE

!     ------------------------------------------------------------------

TYPE(GEOMETRY), INTENT(INOUT) :: YDGEOMETRY
TYPE(TDYN)     ,INTENT(INOUT) :: YDDYN
TYPE(TDYNA)    ,INTENT(INOUT) :: YDDYNA
REAL(KIND=JPRB) :: ZHDIRDIV, ZHDIRVOR, ZDELX, ZHDIR, ZHDSR, ZRHDSR,&
  &  ZRHDIR

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"

!     ------------------------------------------------------------------


IF (LHOOK) CALL DR_HOOK('SUHDIR',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, &
 & YDEGEO=>YDGEOMETRY%YREGEO)
ASSOCIATE(NDLON=>YDDIM%NDLON, NSMAX=>YDDIM%NSMAX, &
 & HDIRDIV=>YDDYN%HDIRDIV, HDIRO3=>YDDYN%HDIRO3, HDIRPD=>YDDYN%HDIRPD, &
 & HDIRQ=>YDDYN%HDIRQ, HDIRSP=>YDDYN%HDIRSP, HDIRT=>YDDYN%HDIRT, &
 & HDIRVD=>YDDYN%HDIRVD, HDIRVOR=>YDDYN%HDIRVOR, HDSRDIV=>YDDYN%HDSRDIV, &
 & HDSRVD=>YDDYN%HDSRVD, HDSRVOR=>YDDYN%HDSRVOR, HRDIRDIV=>YDDYN%HRDIRDIV, &
 & HRDIRO3=>YDDYN%HRDIRO3, HRDIRPD=>YDDYN%HRDIRPD, HRDIRQ=>YDDYN%HRDIRQ, &
 & HRDIRSP=>YDDYN%HRDIRSP, HRDIRT=>YDDYN%HRDIRT, HRDIRVD=>YDDYN%HRDIRVD, &
 & HRDIRVOR=>YDDYN%HRDIRVOR, HRDSRDIV=>YDDYN%HRDSRDIV, HRDSRVD=>YDDYN%HRDSRVD, &
 & HRDSRVOR=>YDDYN%HRDSRVOR, LNEWHD=>YDDYN%LNEWHD, RDAMPDIV=>YDDYN%RDAMPDIV, &
 & RDAMPDIVS=>YDDYN%RDAMPDIVS, RDAMPHDS=>YDDYN%RDAMPHDS, RDAMPO3=>YDDYN%RDAMPO3, &
 & RDAMPPD=>YDDYN%RDAMPPD, RDAMPQ=>YDDYN%RDAMPQ, RDAMPSP=>YDDYN%RDAMPSP, &
 & RDAMPT=>YDDYN%RDAMPT, RDAMPVD=>YDDYN%RDAMPVD, RDAMPVDS=>YDDYN%RDAMPVDS, &
 & RDAMPVOR=>YDDYN%RDAMPVOR, RDAMPVORS=>YDDYN%RDAMPVORS, RRDXTAU=>YDDYN%RRDXTAU, &
 & LSLHD_W=>YDDYNA%LSLHD_W,LSLHD_SVD=>YDDYNA%LSLHD_SVD, &
 & LSLHD_STATIC=>YDDYNA%LSLHD_STATIC, &
 & EDELX=>YDEGEO%EDELX, EDELY=>YDEGEO%EDELY, &
 & RNLGINC=>YDGEM%RNLGINC)
!     ------------------------------------------------------------------

IF (LECMWF .AND. (.NOT. LNEWHD)) THEN


!        1. FOR ECMWF (WHEN THE "OLD STYLE" IS KEPT):
!           COMPUTATION OF HORIZONTAL DIFFUSION COEFFICIENTS.
!        ----------------------------------------------------
!    (this computation is performed before namelist values are readed)

!        HD coefficients (unstretched geometry).

  IF(NSMAX >= 3999) THEN
    ZHDIRVOR=240._JPRB
  ELSEIF(NSMAX >= 2047) THEN
    ZHDIRVOR=450._JPRB
  ELSEIF(NSMAX >= 1279) THEN
    ZHDIRVOR=600._JPRB
  ELSEIF(NSMAX >= 799) THEN
    ZHDIRVOR=1200._JPRB
  ELSEIF(NSMAX >= 511) THEN
    ZHDIRVOR=1800._JPRB
  ELSEIF(NSMAX >= 319) THEN
    ZHDIRVOR=2100._JPRB
  ELSEIF(NSMAX >= 255) THEN
    ZHDIRVOR=2500._JPRB
  ELSEIF(NSMAX >= 213) THEN
    ZHDIRVOR=2700._JPRB
  ELSEIF(NSMAX >= 159) THEN
    ZHDIRVOR=5400._JPRB
  ELSEIF(NSMAX >= 106) THEN
    ZHDIRVOR=12570._JPRB
  ELSEIF(NSMAX >= 95) THEN
    ZHDIRVOR=15000._JPRB
  ELSEIF(NSMAX >= 63) THEN
    ZHDIRVOR=27000._JPRB
  ELSEIF(NSMAX >= 42) THEN
    ZHDIRVOR=54000._JPRB
  ELSEIF(NSMAX >= 21) THEN
    ZHDIRVOR=108000._JPRB
  ELSE
    CALL ABOR1(' WRONG HORIZONTAL RESOLUTION')
  ENDIF

  ! This is just computation from defaults, so RDAMPOR should
  ! not be zero, but anyway...
  IF (RDAMPVOR /= 0._JPRB) THEN
    ZHDIRDIV = ZHDIRVOR / RDAMPVOR
  ELSE
    ZHDIRDIV = 0._JPRB
  ENDIF

  ZHDIRVOR = REAL(NINT(ZHDIRVOR,JPIB),JPRB)
  ZHDIRDIV = REAL(NINT(ZHDIRDIV,JPIB),JPRB)

!     transfer to module HD coefficients.

  HDIRDIV = ZHDIRDIV*RPLRADI
  HDIRVOR = ZHDIRVOR*RPLRADI
  HDIRT   = ZHDIRVOR*RPLRADI
  HDIRQ   = ZHDIRVOR*RPLRADI

  HDIRSP  = 0.0_JPRB
  HDIRO3  = 0.0_JPRB
  HDIRPD  = 0.0_JPRB
  HDIRVD  = 0.0_JPRB

  ! Supporting diffusion (not really of any use in IFS)
  HDSRVOR  = 0.0_JPRB
  HDSRDIV  = 0.0_JPRB
  HDSRVD   = 0.0_JPRB
  HRDSRVOR = 0.0_JPRB
  HRDSRDIV = 0.0_JPRB
  HRDSRVD  = 0.0_JPRB
 
ELSE

!     ------------------------------------------------------------------

!        2. FOR METEO-FRANCE AND LNEWHD OF ECMWF
!           COMPUTATION OF HORIZONTAL DIFFUSION COEFFICIENTS.
!        ----------------------------------------------------

!        2.1: calculation of ZDELX: horizontal grid-point mesh size.

  IF (LELAM) THEN
    ZDELX=SQRT(0.5_JPRB*(EDELX*EDELX+EDELY*EDELY))
  ELSE
    ZDELX=2.0_JPRB*RPI*RA/REAL(NDLON,JPRB)
  ENDIF

!        2.2: horizontal diffusion coefficients.
!             (For ECMWF the HDIR[x] are used,
!              for M-F the reverse HRDIR[x] are used instead)

  IF (RRDXTAU /= 0.0_JPRB) THEN
    ZHDIR = ZDELX*((1.0_JPRB+0.5_JPRB*RNLGINC)**2.5_JPRB)/RRDXTAU
    ZRHDIR=RRDXTAU/(ZDELX*((1.0_JPRB+0.5_JPRB*RNLGINC)**2.5_JPRB))

!    transfer to module HD coefficients

    ! divergence
    IF(RDAMPDIV > 0._JPRB) THEN
      HDIRDIV = ZHDIR * RDAMPDIV
      HRDIRDIV=ZRHDIR/RDAMPDIV
    ELSE
      HDIRDIV =0._JPRB
      HRDIRDIV=0._JPRB
    ENDIF

    ! vorticity
    IF(RDAMPVOR > 0._JPRB) THEN
      HDIRVOR = ZHDIR * RDAMPVOR
      HRDIRVOR=ZRHDIR/RDAMPVOR
    ELSE
      HDIRVOR =0._JPRB
      HRDIRVOR=0._JPRB
    ENDIF

    ! temperature
    IF(RDAMPT > 0._JPRB) THEN
      HDIRT   = ZHDIR * RDAMPT
      HRDIRT  =ZRHDIR/RDAMPT
    ELSE
      HDIRT   =0._JPRB
      HRDIRT  =0._JPRB
    ENDIF

    !  humidity
    IF(RDAMPQ > 0._JPRB) THEN
      HDIRQ   = ZHDIR * RDAMPQ
      HRDIRQ  =ZRHDIR/RDAMPQ
    ELSE
     HDIRQ   =0._JPRB
     HRDIRQ  =0._JPRB
    ENDIF

    ! surface pressure
    IF(RDAMPSP > 0._JPRB) THEN
      HDIRSP  = ZHDIR * RDAMPSP
      HRDIRSP =ZRHDIR/RDAMPSP
    ELSE
      HDIRSP  =0._JPRB
      HRDIRSP =0._JPRB
    ENDIF

    ! ozone
    IF(RDAMPO3 > 0._JPRB) THEN
      HDIRO3  = ZHDIR * RDAMPO3
      HRDIRO3 =ZRHDIR/RDAMPO3
    ELSE
      HDIRO3  =0._JPRB
      HRDIRO3 =0._JPRB
    ENDIF

    ! pressure departure (NH)
    IF(RDAMPPD > 0._JPRB) THEN
      HDIRPD  = ZHDIR * RDAMPPD
      HRDIRPD =ZRHDIR/RDAMPPD
    ELSE
      HDIRPD  =0._JPRB
      HRDIRPD =0._JPRB
    ENDIF

    ! vertical divergence (NH)
    IF(RDAMPVD > 0._JPRB) THEN
      HDIRVD  = ZHDIR * RDAMPVD
      HRDIRVD =ZRHDIR/RDAMPVD
    ELSE
      HDIRVD  =0._JPRB
      HRDIRVD =0._JPRB
    ENDIF

  ELSE
    HDIRVOR = 0.0_JPRB
    HRDIRVOR= 0.0_JPRB
    HDIRDIV = 0.0_JPRB
    HRDIRDIV= 0.0_JPRB
    HDIRT   = 0.0_JPRB
    HRDIRT  = 0.0_JPRB
    HDIRQ   = 0.0_JPRB
    HRDIRQ  = 0.0_JPRB
    HDIRO3  = 0.0_JPRB
    HRDIRO3 = 0.0_JPRB
    HDIRPD  = 0.0_JPRB
    HRDIRPD = 0.0_JPRB
    HDIRVD  = 0.0_JPRB
    HRDIRVD = 0.0_JPRB
    HDIRSP  = 0.0_JPRB
    HRDIRSP = 0.0_JPRB
  ENDIF

  

  !     ------------------------------------------------------------------

  !        3. COMPUTATION OF SLHD SUPPORTING DIFFUSION
  !        ----------------------------------------------------
  !        This diffusion is available just in the modern
  !        spectral diffusion set-up form

  IF ((LSLHD_W .OR. LSLHD_SVD).AND.(RRDXTAU /= 0.0_JPRB)&
   &  .AND.(.NOT.LSLHD_STATIC)) THEN
    ZRHDSR = RDAMPHDS * ZRHDIR
    ZHDSR = 1._JPRB/ZRHDSR

    ! vorticity
    IF (LSLHD_W .AND. (RDAMPVORS /= 0._JPRB)) THEN
      HDSRVOR  = ZHDSR  * RDAMPVORS
      HRDSRVOR = ZRHDSR / RDAMPVORS
    ELSE
      HDSRVOR  = 0._JPRB
      HRDSRVOR = 0._JPRB
    ENDIF

    ! divergence
    IF (LSLHD_W .AND. (RDAMPDIVS /= 0._JPRB)) THEN
      HDSRDIV = ZHDSR  *  RDAMPDIVS
      HRDSRDIV = ZRHDSR / RDAMPDIVS
    ELSE
      HDSRDIV  = 0._JPRB
      HRDSRDIV = 0._JPRB
    ENDIF

    ! vertical divergence (NH)
    IF (LSLHD_SVD .AND. (RDAMPVDS /= 0._JPRB)) THEN
      HDSRVD = ZHDSR * RDAMPVDS
      HRDSRVD = ZRHDSR / RDAMPVDS
    ELSE
      HDSRVD  = 0._JPRB
      HRDSRVD = 0._JPRB
    ENDIF

  ELSE
    HDSRVOR  = 0.0_JPRB
    HDSRDIV  = 0.0_JPRB
    HDSRVD   = 0.0_JPRB
    HRDSRVOR = 0.0_JPRB
    HRDSRDIV = 0.0_JPRB
    HRDSRVD  = 0.0_JPRB
  ENDIF

ENDIF
!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUHDIR',1,ZHOOK_HANDLE)
END SUBROUTINE SUHDIR
