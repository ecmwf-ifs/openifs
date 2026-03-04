! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SUINIF21C_NC(YDDIMV,PDUG0,PDVG0,PDVVEL0,PDUADV,PDVADV,&
                        &PDTADV,PDQADV, PDETADOT, PDEXTSHF, PDEXTLHF)

!**** *SUINIF21C_NC* - Read time-varying large-scale forcing.

!     Purpose.
!     --------
!           Read time-varying large-scale dynamical
!           forcing of the single column model.

!**   Interface.
!     ----------
!        *CALL* *SUINIF21C_NC*

!        Explicit arguments :
!        --------------------
!         PD..(NFLEVG) - tendencies of the large-scale forcing 

!        Implicit arguments :
!        --------------------
!         NONE

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!        Called by CNT41C

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the SCM


!     Author.
!     -------
!        Joao Teixeira   *ECMWF*

!     Modifications.
!     --------------
!        Original :        94-04-28
!                          98-02     From GF arrays to UG0,... arrays
!        Martin Koehler: 2000-12-04  modification to NetCDF 
!     ------------------------------------------------------------------

USE YOMDIMV  , ONLY : TDIMV
USE PARKIND1 , ONLY : JPIM     ,JPRB     ,JPRM
USE YOMHOOK  , ONLY : LHOOK    ,DR_HOOK, JPHOOK
USE YOMLUN   , ONLY : NULOUT
USE YOMCT0   , ONLY : REXTSHF  ,REXTLHF
USE YOMCT3   , ONLY : NSTEP
USE YOMGF1C  , ONLY : UG0      ,VG0      ,VVEL0    ,UADV     ,&
                     &VADV     ,TADV     ,QADV     ,ETADOTDPDETA
USE YOMLOG1C , ONLY : LDYNFOR  ,LUGVG    ,LVERVEL  ,LETADOT  ,&
                     &LTADV    ,LQADV    ,LUVADV   ,NSTRTINI ,&
                     &NFRFOR

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TDIMV), INTENT(INOUT) :: YDDIMV
REAL(KIND=JPRB)    :: PDUG0(YDDIMV%NFLEVG)  , PDVG0(YDDIMV%NFLEVG)  , PDVVEL0(YDDIMV%NFLEVG)
REAL(KIND=JPRB)    :: PDUADV(YDDIMV%NFLEVG) , PDVADV(YDDIMV%NFLEVG)
REAL(KIND=JPRB)    :: PDTADV(YDDIMV%NFLEVG) , PDQADV(YDDIMV%NFLEVG) , PDETADOT(0:YDDIMV%NFLEVG)
REAL(KIND=JPRB)    :: PDEXTSHF       , PDEXTLHF
REAL(KIND=JPRM)    :: TEMP0, TEMP1A(YDDIMV%NFLEVG) , TEMP2A(0:YDDIMV%NFLEVG)
REAL(KIND=JPRB)    :: ZUG0(2,YDDIMV%NFLEVG) , ZVG0(2,YDDIMV%NFLEVG) , ZVVEL0(2,YDDIMV%NFLEVG)
REAL(KIND=JPRB)    :: ZUADV(2,YDDIMV%NFLEVG), ZVADV(2,YDDIMV%NFLEVG)
REAL(KIND=JPRB)    :: ZTADV(2,YDDIMV%NFLEVG), ZQADV(2,YDDIMV%NFLEVG), ZETADOT(2,0:YDDIMV%NFLEVG)
REAL(KIND=JPRB)    :: ZEXTSHF(2)     , ZEXTLHF(2)

INTEGER(KIND=JPIM) :: IINT, JALEV, JTIME
INTEGER(KIND=JPIM) :: INCID, VARID, START1, COUNT1, START2(2), COUNT2(2), ISTATUS, NT

REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

#include "netcdf.inc"

!     ------------------------------------------------------------------
#include "handle_err_nc.intfb.h"
!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUINIF21C_NC',0,ZHOOK_HANDLE)
ASSOCIATE(NFLEVG=>YDDIMV%NFLEVG)
!     ------------------------------------------------------------------

!*       1.    OPEN INPUT FILE.
!              ----------------

ISTATUS = NF_OPEN ('scm_in.nc', NF_NOWRITE, INCID)
CALL HANDLE_ERR_NC(ISTATUS)


!*       2.    SETUP FOR INTERPOLATION.
!              ------------------------

IF (NSTEP == 0.0_JPRB) THEN
  IINT=1
ELSE
  IINT=2

!     RELATED WITH SWAPPING AND LINEAR INTERPOLATIONS IN TIME.
!     ...time level 1  =  previous model time step variables

  ZUG0   (1,1:NFLEVG) = UG0         (1:NFLEVG) + PDUG0   (1:NFLEVG)
  ZVG0   (1,1:NFLEVG) = VG0         (1:NFLEVG) + PDVG0   (1:NFLEVG)
  ZVVEL0 (1,1:NFLEVG) = VVEL0       (1,1:NFLEVG) + PDVVEL0 (1:NFLEVG)

  ZUADV  (1,1:NFLEVG) = UADV        (1:NFLEVG) + PDUADV  (1:NFLEVG)
  ZVADV  (1,1:NFLEVG) = VADV        (1:NFLEVG) + PDVADV  (1:NFLEVG)
  ZTADV  (1,1:NFLEVG) = TADV        (1:NFLEVG) + PDTADV  (1:NFLEVG)
  ZQADV  (1,1:NFLEVG) = QADV        (1:NFLEVG) + PDQADV  (1:NFLEVG)

  ZETADOT(1,0:NFLEVG) = ETADOTDPDETA(0:NFLEVG) + PDETADOT(0:NFLEVG)
 
  ZEXTSHF(1)          = REXTSHF                + PDEXTSHF
  ZEXTLHF(1)          = REXTLHF                + PDEXTLHF

ENDIF

DO JTIME=IINT,2


!        3.    READ GEOSTROPHIC WIND AND VERTICAL VELOCITY.
!              --------------------------------------------

! ... data step selection:
! ... (nstep=model step (0...n), nfrfor=freq. forcing)
  IF ( JTIME == 1 ) THEN
    NT = NSTRTINI
  ELSE
    NT = NSTRTINI + ( NSTEP / NFRFOR ) + 1  !apparently needs to be one ahead? 
  ENDIF
  START2 = (/ 1     , NT /)
  COUNT2 = (/ NFLEVG, 1  /)

  IF ( LDYNFOR ) THEN

    ISTATUS = NF_INQ_VARID     (INCID, 'ug', VARID)
    CALL HANDLE_ERR_NC(ISTATUS)
    ISTATUS = NF_GET_VARA_REAL (INCID, VARID, START2, COUNT2, TEMP1A)
    CALL HANDLE_ERR_NC(ISTATUS)
    UG0(1:NFLEVG) = TEMP1A

    ISTATUS = NF_INQ_VARID     (INCID, 'vg', VARID)
    CALL HANDLE_ERR_NC(ISTATUS)
    ISTATUS = NF_GET_VARA_REAL (INCID, VARID, START2, COUNT2, TEMP1A)
    CALL HANDLE_ERR_NC(ISTATUS)
    VG0(1:NFLEVG) = TEMP1A

    ISTATUS = NF_INQ_VARID     (INCID, 'omega', VARID)
    CALL HANDLE_ERR_NC(ISTATUS)
    ISTATUS = NF_GET_VARA_REAL (INCID, VARID, START2, COUNT2, TEMP1A)
    CALL HANDLE_ERR_NC(ISTATUS)
    VVEL0(1,1:NFLEVG) = TEMP1A

    COUNT2 = (/ NFLEVG+1, 1  /)
    ISTATUS = NF_INQ_VARID     (INCID, 'etadotdpdeta', VARID)
    CALL HANDLE_ERR_NC(ISTATUS)
    ISTATUS = NF_GET_VARA_REAL (INCID, VARID, START2, COUNT2, TEMP2A)
    CALL HANDLE_ERR_NC(ISTATUS)
    ETADOTDPDETA(0:NFLEVG) = TEMP2A

    IF (.NOT.LUGVG) THEN
      UG0   =0.0_JPRB
      VG0   =0.0_JPRB
    ENDIF
    IF (.NOT.LVERVEL) THEN
      VVEL0 =0.0_JPRB
    ENDIF
    IF (.NOT.LETADOT) THEN
      ETADOTDPDETA=0.0_JPRB
    ENDIF


!        4.    READ HORIZONTAL ADVECTION.
!              --------------------------

    COUNT2 = (/ NFLEVG, 1  /)
    ISTATUS = NF_INQ_VARID     (INCID, 'uadv', VARID)
    CALL HANDLE_ERR_NC(ISTATUS)
    ISTATUS = NF_GET_VARA_REAL (INCID, VARID, START2, COUNT2, TEMP1A)
    CALL HANDLE_ERR_NC(ISTATUS)
    UADV(1:NFLEVG) = TEMP1A

    ISTATUS = NF_INQ_VARID     (INCID, 'vadv', VARID)
    CALL HANDLE_ERR_NC(ISTATUS)
    ISTATUS = NF_GET_VARA_REAL (INCID, VARID, START2, COUNT2, TEMP1A)
    CALL HANDLE_ERR_NC(ISTATUS)
    VADV(1:NFLEVG) = TEMP1A

    ISTATUS = NF_INQ_VARID     (INCID, 'tadv', VARID)
    CALL HANDLE_ERR_NC(ISTATUS)
    ISTATUS = NF_GET_VARA_REAL (INCID, VARID, START2, COUNT2, TEMP1A)
    CALL HANDLE_ERR_NC(ISTATUS)
    TADV(1:NFLEVG) = TEMP1A

    ISTATUS = NF_INQ_VARID     (INCID, 'qadv', VARID)
    CALL HANDLE_ERR_NC(ISTATUS)
    ISTATUS = NF_GET_VARA_REAL (INCID, VARID, START2, COUNT2, TEMP1A)
    CALL HANDLE_ERR_NC(ISTATUS)
    QADV(1:NFLEVG) = TEMP1A

    IF (.NOT.LUVADV) THEN
      UADV=0.0_JPRB
      VADV=0.0_JPRB
    ENDIF
    IF (.NOT.LTADV) THEN
      TADV=0.0_JPRB
    ENDIF
    IF (.NOT.LQADV) THEN
      QADV=0.0_JPRB
    ENDIF

  ELSE

    UG0  =0.0_JPRB
    VG0  =0.0_JPRB
    VVEL0=0.0_JPRB
    UADV =0.0_JPRB
    VADV =0.0_JPRB
    TADV =0.0_JPRB
    QADV =0.0_JPRB
    ETADOTDPDETA =0.0_JPRB

  ENDIF


!        5.    READ SURFACE FLUXES.
!              --------------------
 
  START1 = NT
  COUNT1 = 1

  ISTATUS = NF_INQ_VARID     (INCID, 'sfc_sens_flx', VARID)
  CALL HANDLE_ERR_NC(ISTATUS)
  ISTATUS = NF_GET_VARA_REAL (INCID, VARID, START1, COUNT1, TEMP0)
  CALL HANDLE_ERR_NC(ISTATUS)
  REXTSHF=TEMP0

  ISTATUS = NF_INQ_VARID     (INCID, 'sfc_lat_flx', VARID)
  CALL HANDLE_ERR_NC(ISTATUS)
  ISTATUS = NF_GET_VARA_REAL (INCID, VARID, START1, COUNT1, TEMP0)
  CALL HANDLE_ERR_NC(ISTATUS)
  REXTLHF=TEMP0


!        6.    STORING.
!              --------

  ZUG0   (JTIME,1:NFLEVG) = UG0  (1:NFLEVG)
  ZVG0   (JTIME,1:NFLEVG) = VG0  (1:NFLEVG)
  ZVVEL0 (JTIME,1:NFLEVG) = VVEL0(1,1:NFLEVG)

  ZUADV  (JTIME,1:NFLEVG) = UADV (1:NFLEVG)
  ZVADV  (JTIME,1:NFLEVG) = VADV (1:NFLEVG)
  ZTADV  (JTIME,1:NFLEVG) = TADV (1:NFLEVG)
  ZQADV  (JTIME,1:NFLEVG) = QADV (1:NFLEVG)

  ZETADOT(JTIME,0:NFLEVG) = ETADOTDPDETA(0:NFLEVG)

  ZEXTSHF(JTIME)         = REXTSHF
  ZEXTLHF(JTIME)         = REXTLHF

ENDDO


!        7.    LINEAR INTERPOLATION IN TIME.
!              -----------------------------

PDUG0   (1:NFLEVG) = (ZUG0   (2,1:NFLEVG) - ZUG0   (1,1:NFLEVG)) / NFRFOR
PDVG0   (1:NFLEVG) = (ZVG0   (2,1:NFLEVG) - ZVG0   (1,1:NFLEVG)) / NFRFOR
PDVVEL0 (1:NFLEVG) = (ZVVEL0 (2,1:NFLEVG) - ZVVEL0 (1,1:NFLEVG)) / NFRFOR

PDUADV  (1:NFLEVG) = (ZUADV  (2,1:NFLEVG) - ZUADV  (1,1:NFLEVG)) / NFRFOR
PDVADV  (1:NFLEVG) = (ZVADV  (2,1:NFLEVG) - ZVADV  (1,1:NFLEVG)) / NFRFOR
PDTADV  (1:NFLEVG) = (ZTADV  (2,1:NFLEVG) - ZTADV  (1,1:NFLEVG)) / NFRFOR
PDQADV  (1:NFLEVG) = (ZQADV  (2,1:NFLEVG) - ZQADV  (1,1:NFLEVG)) / NFRFOR

PDETADOT(0:NFLEVG) = (ZETADOT(2,0:NFLEVG) - ZETADOT(1,0:NFLEVG)) / NFRFOR

PDEXTSHF           = (ZEXTSHF(2)          - ZEXTSHF(1))         / NFRFOR
PDEXTLHF           = (ZEXTLHF(2)          - ZEXTLHF(1))         / NFRFOR

!        8.    SWAPPING (previous model time step .
!              ---------

UG0         (1:NFLEVG) = ZUG0   (1,1:NFLEVG)
VG0         (1:NFLEVG) = ZVG0   (1,1:NFLEVG)
VVEL0       (1,1:NFLEVG) = ZVVEL0 (1,1:NFLEVG)

UADV        (1:NFLEVG) = ZUADV  (1,1:NFLEVG)
VADV        (1:NFLEVG) = ZVADV  (1,1:NFLEVG)
TADV        (1:NFLEVG) = ZTADV  (1,1:NFLEVG)
QADV        (1:NFLEVG) = ZQADV  (1,1:NFLEVG)

ETADOTDPDETA(0:NFLEVG) = ZETADOT(1,0:NFLEVG)

REXTSHF                = ZEXTSHF(1)
REXTLHF                = ZEXTLHF(1)

!        9.    CLOSE INPUT FILE.
!              -----------------

ISTATUS = NF_CLOSE (INCID)
CALL HANDLE_ERR_NC(ISTATUS)


!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUINIF21C_NC',1,ZHOOK_HANDLE)
END SUBROUTINE SUINIF21C_NC


