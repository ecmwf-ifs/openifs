! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE CNT41C(YDGEOMETRY,YDMODEL,YDFIELDS)
!SUBROUTINE CNT41C(YDGEOMETRY,YDSURF)

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE TYPE_MODEL   , ONLY : MODEL
USE MODEL_PHYSICS_RADIATION_MOD , ONLY : MODEL_PHYSICS_RADIATION_TYPE
USE FIELDS_MOD , ONLY : FIELDS
USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK    ,DR_HOOK, JPHOOK
USE PARDIM1C
USE YOMLUN   , ONLY : NULOUT
USE YOMCT0   , ONLY : REXTSHF  ,REXTLHF
USE YOMCT2   , ONLY : NSTAR2   ,NSTOP2
USE YOMCT3   , ONLY : NSTEP
USE YOMLOG1C
USE YOMGF1C  , ONLY : UG0      ,VG0      ,VVEL0    ,UADV     ,&
                     &VADV     ,TADV     ,QADV     ,ETADOTDPDETA
USE YOMGP1C0
USE YOMCST
USE YOMGT1C0

!**** *CNT41C*  - Controls integration job at level 4

!     Purpose.
!     --------
!     CONTROLS THE INTEGRATION

!**   Interface.
!     ----------
!        *CALL* *CNT41C

!        Explicit arguments :
!        --------------------
!        None

!        Implicit arguments :
!        --------------------
!        None

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!        Called by CNT1C.

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the single column model

!     Author.
!     -------
!        Joao Teixeira  *ECMWF*

!     Modifications.
!     --------------
!        Original : 93-12-31
!        J.Teixeira, Jan-95: new output files.
!        J.Teixeira, May-95: two time level scheme.
!        J.Teixeira, Feb-98: lat, lon, sst, ps varying in time.
!                   + relaxation to obs: u,v,T,q.
!        M. Ko"hler  6-6-2006 Single Column Model integration within IFS 
!----------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY), INTENT(INOUT) :: YDGEOMETRY
TYPE(MODEL),    INTENT(INOUT) :: YDMODEL
TYPE(FIELDS),   INTENT(INOUT) :: YDFIELDS

REAL(KIND=JPRB) :: ZDUG0(YDGEOMETRY%YRDIMV%NFLEVG) ,ZDVG0(YDGEOMETRY%YRDIMV%NFLEVG),ZDVVEL0(YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZDUADV(YDGEOMETRY%YRDIMV%NFLEVG),ZDVADV(YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZDTADV(YDGEOMETRY%YRDIMV%NFLEVG),ZDQADV(YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZDETADOT(0:YDGEOMETRY%YRDIMV%NFLEVG), ZDEXTSHF, ZDEXTLHF

REAL(KIND=JPRB) :: ZLA(5), ZLO(5), ZST(5), ZSP(5)

REAL(KIND=JPRB) :: ZUOBS(5,YDGEOMETRY%YRDIMV%NFLEVG),ZVOBS(5,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZTOBS(5,YDGEOMETRY%YRDIMV%NFLEVG),ZQOBS(5,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZDUOBS(YDGEOMETRY%YRDIMV%NFLEVG) ,ZDVOBS(YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZDTOBS(YDGEOMETRY%YRDIMV%NFLEVG) ,ZDQOBS(YDGEOMETRY%YRDIMV%NFLEVG)

!     LOCAL INTEGER SCALARS
INTEGER(KIND=JPIM) :: IA, IC, IINT, ISTOP, JALEV, JL, JSTEP, JTIME, NT

!     LOCAL REAL SCALARS
REAL(KIND=JPRB) :: ZDLAT, ZDLON, ZDSPL, ZDSST, ZLAT, ZLON, ZSPL, ZSST

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!      -----------------------------------------------------------
#include "suinif21c_nc.intfb.h"
#include "ecradfr.intfb.h"
#include "stepo1c.intfb.h"
#include "updtim.intfb.h"
#include "read_obs.intfb.h"
#include "read_sst.intfb.h"
!      -----------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CNT41C',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, &
 & YDGSGEOM=>YDGEOMETRY%YRGSGEOM,YDGSGEOM_NB=>YDGEOMETRY%YRGSGEOM_NB,  &
 & YDRIP=>YDMODEL%YRML_GCONF%YRRIP)
ASSOCIATE(NFLEVG=>YDDIMV%NFLEVG, &
 & LTWOTL=>YDMODEL%YRML_DYN%YRDYNA%LTWOTL, &
 & NSTOP=>YDRIP%NSTOP, TDT=>YDRIP%TDT, &
 & TSTEP=>YDRIP%TSTEP, &
 & YSD_VF=>YDFIELDS%YRSURF%YSD_VF, YDML_PHY_RAD=>YDMODEL%YRML_PHY_RAD)
!      -----------------------------------------------------------


!*       1.    PREPARING MAIN LOOP.
!              --------------------


IF(NSTOP2 == NSTOP)THEN
  ISTOP=NSTOP
ELSE
  ISTOP=NSTOP2-1
ENDIF


!*       2.    OPEN FILES.
!              -----------   

!*       2.1   OPEN OUTPUT FILES.
!              -----------------

OPEN(NPOSASC,FILE='onecol.r')

IF ( OUTFORM /= 'netcdf' ) THEN

  OPEN(51,FILE='upperair.r')
  OPEN(53,FILE='diagvar.r')
  OPEN(54,FILE='uatend.r')
  OPEN(55,FILE='uaflux.r')
  OPEN(56,FILE='surfsoil.r')

ENDIF


!*       3.    MAIN TIME LOOP.
!              ---------------

DO JSTEP=NSTAR2,ISTOP

!*       3.1   CURRENT VALUE OF TIME STEP LENGTH.
!              ----------------------------------

  WRITE(NULOUT,*) 'STEP=',JSTEP

  IF(JSTEP == 0.OR.LTWOTL)THEN
    TDT=TSTEP
  ELSE
    TDT=2.0_JPRB*TSTEP
  ENDIF


!*       3.2   RESET OF TIME DEPENDENT CONSTANTS.
!              ----------------------------------

  CALL UPDTIM(YDGEOMETRY,YDFIELDS%YRSURF,YDMODEL,JSTEP,TDT,TSTEP,.FALSE.)

!*       3.3   FREQUENCY OF FULL-RADIATION COMPUTATIONS.
!              -----------------------------------------

  NSTEP=JSTEP

  CALL ECRADFR(YDDIM,YDML_PHY_RAD,YDMODEL%YRML_DYN%YRDYNA,YDRIP)

!*       3.4.1 READ CHANGING LAT,LON,SST.
!              --------------------------

  IF(LVARSST) THEN

    IF ( MOD(NSTEP,NFRSST) == 0 .AND. NSTEP /= ISTOP ) THEN

      WRITE(NULOUT,*) 'TIME-VARYING LAT/LON'

      IF (NSTEP == 0.0_JPRB) THEN
        IINT=1
      ELSE
        IINT=2
      ENDIF

      DO JTIME=IINT,2

        IF ( JTIME == 1 ) THEN
          NT = NSTRTINI
        ELSE
          NT = NSTRTINI + ( NSTEP / NFRSST ) + 1  !apparently needs to be one ahead? 
        ENDIF
        CALL READ_SST (NT, ZLAT, ZLON, ZSST, ZSPL)

        ZLA(JTIME) = ZLAT
        ZLO(JTIME) = ZLON
        ZST(JTIME) = ZSST
        ZSP(JTIME) = ZSPL

      ENDDO

      ZDLAT = (ZLA(2) - ZLA(1)) / NFRSST
      ZDLON = (ZLO(2) - ZLO(1)) / NFRSST
      ZDSST = (ZST(2) - ZST(1)) / NFRSST
      ZDSPL = (ZSP(2) - ZSP(1)) / NFRSST

      ZLAT = ZLA(1)
      ZLON = ZLO(1)
      ZSST = ZST(1)
      ZSPL = ZSP(1)

      ZLA(1) = ZLA(2)
      ZLO(1) = ZLO(2)
      ZST(1) = ZST(2)
      ZSP(1) = ZSP(2)

    ELSE

      WRITE(NULOUT,*) 'INTERPOLATED LAT/LON'
      ZLAT = ZLAT + ZDLAT
      ZLON = ZLON + ZDLON
      ZSST = ZSST + ZDSST
      ZSPL = ZSPL + ZDSPL

    ENDIF

    YDGSGEOM(1)%GELAT = ZLAT * RPI / 180._JPRB
    YDGSGEOM(1)%GELAM = ZLON * RPI / 180._JPRB
    TL0   = ZSST
    SPT0  = ZSPL
    YDFIELDS%YRSURF%SD_VF(1,YSD_VF%YSST%MP,1)  = ZSST

    IF (YDGSGEOM(1)%GELAM(1) < 0.0_JPRB) THEN

      YDGSGEOM(1)%GELAM = YDGSGEOM(1)%GELAM + 2.0_JPRB * RPI

    ENDIF

    YDGSGEOM(1)%GEMU  = SIN(YDGSGEOM(1)%GELAT)
    YDGSGEOM(1)%GSQM2 = COS(YDGSGEOM(1)%GELAT)

    YDGSGEOM(1)%GESLO = SIN(YDGSGEOM(1)%GELAM)
    YDGSGEOM(1)%GECLO = COS(YDGSGEOM(1)%GELAM)

    YDGSGEOM(1)%RCORI = 2.0_JPRB*ROMEGA*YDGSGEOM(1)%GEMU

  ENDIF


!*       3.4.2 READ TIME-VARYING LARGE-SCALE FORCING.
!              --------------------------------------

  IF (LVARFOR) THEN

    IF ( MOD(NSTEP,NFRFOR) == 0 .AND. NSTEP /= ISTOP ) THEN

      WRITE(NULOUT,*) 'TIME-VARYING LARGE-SCALE FORCING'
      CALL SUINIF21C_NC(YDGEOMETRY%YRDIMV,ZDUG0,ZDVG0,ZDVVEL0,ZDUADV,ZDVADV,&
        &ZDTADV,ZDQADV, ZDETADOT, ZDEXTSHF, ZDEXTLHF)

    ELSE

      WRITE(NULOUT,*) 'INTERPOLATION OR LAST TIME STEP'
      UG0  (1:NFLEVG)        = UG0  (1:NFLEVG)        + ZDUG0   (1:NFLEVG)
      VG0  (1:NFLEVG)        = VG0  (1:NFLEVG)        + ZDVG0   (1:NFLEVG)
      VVEL0(1,1:NFLEVG)      = VVEL0(1,1:NFLEVG)      + ZDVVEL0 (1:NFLEVG)

      UADV (1:NFLEVG)        = UADV (1:NFLEVG)        + ZDUADV  (1:NFLEVG)
      VADV (1:NFLEVG)        = VADV (1:NFLEVG)        + ZDVADV  (1:NFLEVG)
      TADV (1:NFLEVG)        = TADV (1:NFLEVG)        + ZDTADV  (1:NFLEVG)
      QADV (1:NFLEVG)        = QADV (1:NFLEVG)        + ZDQADV  (1:NFLEVG)

      ETADOTDPDETA(0:NFLEVG) = ETADOTDPDETA(0:NFLEVG) + ZDETADOT(0:NFLEVG)

      REXTSHF                = REXTSHF                + ZDEXTSHF
      REXTLHF                = REXTLHF                + ZDEXTLHF

    ENDIF

  ENDIF

!        3.4.3. READ RELAXATION FIELDS.         
!               -----------------------

  IF(LRELAX) THEN

    IF ( MOD(NSTEP,NFROBS) == 0 .AND. NSTEP /= ISTOP ) THEN

      IF (NSTEP == 0.0_JPRB) THEN
        IINT=1
      ELSE
        IINT=2
      ENDIF

      DO JTIME=IINT,2

        DO JL=1,NFLEVG
          IF ( JTIME == 1 ) THEN
            NT = NSTRTINI
          ELSE
            NT = NSTRTINI + ( NSTEP / NFROBS ) + 1  !apparently needs to be one ahead? 
          ENDIF
          CALL READ_OBS(YDGEOMETRY%YRDIMV,NT, UOBS(1:YDDIMV%NFLEVG), VOBS(1:YDDIMV%NFLEVG), &
 &                      TOBS(1:YDDIMV%NFLEVG), QOBS(1:YDDIMV%NFLEVG))

          ZUOBS(JTIME,JL) = UOBS(JL)
          ZVOBS(JTIME,JL) = VOBS(JL)
          ZTOBS(JTIME,JL) = TOBS(JL)
          ZQOBS(JTIME,JL) = QOBS(JL)
        ENDDO

      ENDDO

      DO JL=1,NFLEVG
        ZDUOBS(JL) = (ZUOBS(2,JL) - ZUOBS(1,JL)) / NFROBS
        ZDVOBS(JL) = (ZVOBS(2,JL) - ZVOBS(1,JL)) / NFROBS
        ZDTOBS(JL) = (ZTOBS(2,JL) - ZTOBS(1,JL)) / NFROBS
        ZDQOBS(JL) = (ZQOBS(2,JL) - ZQOBS(1,JL)) / NFROBS
      ENDDO

      DO JL=1,NFLEVG
        UOBS(JL) = ZUOBS(1,JL)
        VOBS(JL) = ZVOBS(1,JL)
        TOBS(JL) = ZTOBS(1,JL)
        QOBS(JL) = ZQOBS(1,JL)
      ENDDO

      DO JL=1,NFLEVG
        ZUOBS(1,JL) = ZUOBS(2,JL)
        ZVOBS(1,JL) = ZVOBS(2,JL)
        ZTOBS(1,JL) = ZTOBS(2,JL)
        ZQOBS(1,JL) = ZQOBS(2,JL)
      ENDDO

      WRITE(NULOUT,*) 'TIME-VARYING RELAXATION'

    ELSE

      DO JL=1,NFLEVG
        UOBS(JL) = UOBS(JL) + ZDUOBS(JL)
        VOBS(JL) = VOBS(JL) + ZDVOBS(JL)
        TOBS(JL) = TOBS(JL) + ZDTOBS(JL)
        QOBS(JL) = QOBS(JL) + ZDQOBS(JL)
      ENDDO

      WRITE(NULOUT,*) 'INTERPOLATED RELAXATION'

    ENDIF

  ENDIF

!*       3.5   COMPUTATION OF THE ENTIRE TIME STEP.
!              ------------------------------------

  CALL STEPO1C(YDGEOMETRY,YDMODEL,YDFIELDS)

ENDDO


!*       4.    CLOSE FILES.
!              ------------

!*       4.1   CLOSE OUTPUT FILES.
!              ------------------

CLOSE(NPOSASC)
CLOSE(51)
CLOSE(53)
CLOSE(54)
CLOSE(55)
CLOSE(56)


!      ----------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CNT41C',1,ZHOOK_HANDLE)
END SUBROUTINE CNT41C
