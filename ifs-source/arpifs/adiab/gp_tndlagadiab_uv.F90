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

!OCL  NOEVAL
SUBROUTINE GP_TNDLAGADIAB_UV(YDGEOMETRY,YDGMV,YDEPHY,YDDYN,KST,KPROF,PRCORI,PGNORDL,PGNORDM,&
 & PSGRTL,PSGRTM,PU,PV,PTNDU,PTNDV,PTNDU_NOC,PTNDV_NOC)

!**** *GP_TNDLAGADIAB_UV*   Compute adiabatic Lagrangian tendency of horizontal wind.

!     Purpose.
!     --------
!          Compute adiabatic Lagrangian tendency of horizontal wind, with the following assumptions:
!          - explicit representation of Coriolis term.
!          - no curvature term.
!          - Rayleigh friction taken into account.
!          - pressure gradient term taken into account.

!**   Interface.
!     ----------
!        *CALL* *GP_TNDLAGADIAB_UV(..)

!        Explicit arguments :
!        --------------------

!        INPUT:
!          KST       - first element of work.
!          KPROF     - depth of work.
!          PRCORI    - Coriolis parameter "f = 2 Omega sin(theta)".
!          PGNORDL   - zonal component ("Gnordl") of the unit vector
!                      directed towards the true North pole.
!          PGNORDM   - meridian component ("Gnordm") of the unit vector
!                      directed towards the true North pole.
!          PSGRTL    - zonal component of the pressure force grad.
!          PSGRTM    - merid component of the pressure force grad.
!          PGMV      - GMV variables at t-dt and t.

!        OUTPUT:
!          PTNDU     - Tendency for U-wind.
!          PTNDV     - Tendency for V-wind.
!          PTNDU_NOC - Tendency for U-wind without Coriolis term.
!          PTNDV_NOC - Tendency for V-wind without Coriolis term.

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!           none

!     Reference.
!     ----------
!             Arpege documentation about model equations.

!     Author.
!     -------
!        K. YESSAD (METEO FRANCE/CNRM/GMAP)
!         after some code present in LAVENT, CPEULDYN, LATTEX.
!        Original : SEPT 2010.

! Modifications
! -------------
!  K. Yessad (Nov 2011): more flexible Rayleigh friction, and adaptations.
!  K. Yessad (July 2014): Move some variables.
!  K. Yessad (Feb 2018): remove deep-layer formulations.
!  H. Petithomme (Dec 2020): test reordering and optimisation
!------------------------------------------------------------------------------
! End Modifications

USE PARKIND1,ONLY: JPIM,JPRB
USE GEOMETRY_MOD,ONLY: GEOMETRY
USE YOMGMV,ONLY: TGMV
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCT0,ONLY: LRPLANE
USE YOMDYN,ONLY: TDYN
USE YOEPHY,ONLY: TEPHY

IMPLICIT NONE

TYPE(GEOMETRY),INTENT(IN) :: YDGEOMETRY
TYPE(TGMV),INTENT(INOUT) :: YDGMV
TYPE(TDYN),INTENT(IN) :: YDDYN
TYPE(TEPHY),INTENT(IN) :: YDEPHY
INTEGER(KIND=JPIM),INTENT(IN) :: KST,KPROF
REAL(KIND=JPRB),INTENT(IN) :: PRCORI(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB),INTENT(IN) :: PGNORDL(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB),INTENT(IN) :: PGNORDM(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB),INTENT(IN) :: PSGRTL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB),INTENT(IN) :: PSGRTM(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB),INTENT(IN) :: PU(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB),INTENT(IN) :: PV(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB),INTENT(OUT) :: PTNDU(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB),INTENT(OUT) :: PTNDV(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB),INTENT(OUT) :: PTNDU_NOC(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB),INTENT(OUT) :: PTNDV_NOC(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

INTEGER(KIND=JPIM) :: JROF,JLEV
REAL(KIND=JPRB) :: ZGNORDLM(YDGEOMETRY%YRDIM%NPROMA),ZGNORDL2(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB) :: ZGNORDM2(YDGEOMETRY%YRDIM%NPROMA)

IF (LHOOK) CALL DR_HOOK('GP_TNDLAGADIAB_UV',0,ZHOOK_HANDLE)

ASSOCIATE(DIMV=>YDGEOMETRY%YRDIMV,GEM=>YDGEOMETRY%YRGEM)
ASSOCIATE(NFLEVG=>DIMV%NFLEVG,RKRF=>YDDYN%RKRF,LEGWWMS=>YDEPHY%LEGWWMS,&
  NSTTYP=>GEM%NSTTYP,YT0=>YDGMV%YT0)

IF (LEGWWMS.OR..NOT.YDDYN%LRFRIC) THEN
  ! no Rayleigh friction:
  DO JROF=KST,KPROF
    PTNDU_NOC(JROF,1:NFLEVG) = -PSGRTL(JROF,1:NFLEVG)
    PTNDV_NOC(JROF,1:NFLEVG) = -PSGRTM(JROF,1:NFLEVG)
  ENDDO
ELSE IF (YDDYN%LRFRICISOTR) THEN
  ! isotropic Rayleigh friction:
  DO JLEV=1,NFLEVG
    DO JROF=KST,KPROF
      PTNDU_NOC(JROF,JLEV) = -RKRF(JLEV)*PU(JROF,JLEV)-PSGRTL(JROF,JLEV)
      PTNDV_NOC(JROF,JLEV) = -RKRF(JLEV)*PV(JROF,JLEV)-PSGRTM(JROF,JLEV)
    ENDDO
  ENDDO
ELSE IF (.NOT.LRPLANE.AND.NSTTYP == 1) THEN
  ! non-isotropic Rayleigh friction in untilted spherical geometry:
  DO JLEV=1,NFLEVG
    DO JROF=KST,KPROF
      PTNDU_NOC(JROF,JLEV) = -RKRF(JLEV)*PU(JROF,JLEV)-PSGRTL(JROF,JLEV)
      PTNDV_NOC(JROF,JLEV) = -PSGRTM(JROF,JLEV)
    ENDDO
  ENDDO
ELSE
  ! non-isotropic Rayleigh friction, other cases:
  DO JROF=KST,KPROF
    ZGNORDLM(JROF) = PGNORDM(JROF)*PGNORDL(JROF)
    ZGNORDL2(JROF) = PGNORDL(JROF)**2
    ZGNORDM2(JROF) = PGNORDM(JROF)**2
  ENDDO

  ! simplified version from (with rkrfu=rkrf, rkrfv=0):
  ! ptndu=(rkrfu-rkrfv)*nordm*nordl*v-(rkrfu*nordm**2+rkrfv*nordl**2)*u
  ! ptndv=(rkrfu-rkrfv)*nordm*nordl*u-(rkrfu*nordl**2+rkrfv*nordm**2)*v
  DO JLEV=1,NFLEVG
    DO JROF=KST,KPROF
      PTNDU_NOC(JROF,JLEV) = RKRF(JLEV)*(ZGNORDLM(JROF)*PV(JROF,JLEV)-&
        ZGNORDM2(JROF)*PU(JROF,JLEV))-PSGRTL(JROF,JLEV)
      PTNDV_NOC(JROF,JLEV) = RKRF(JLEV)*(ZGNORDLM(JROF)*PU(JROF,JLEV)-&
        ZGNORDL2(JROF)*PV(JROF,JLEV))-PSGRTM(JROF,JLEV)
    ENDDO
  ENDDO
ENDIF

DO JLEV=1,NFLEVG
  DO JROF=KST,KPROF
    PTNDU(JROF,JLEV) = PTNDU_NOC(JROF,JLEV)+PRCORI(JROF)*PV(JROF,JLEV)
    PTNDV(JROF,JLEV) = PTNDV_NOC(JROF,JLEV)-PRCORI(JROF)*PU(JROF,JLEV)
  ENDDO
ENDDO
END ASSOCIATE
END ASSOCIATE

IF (LHOOK) CALL DR_HOOK('GP_TNDLAGADIAB_UV',1,ZHOOK_HANDLE)
END SUBROUTINE GP_TNDLAGADIAB_UV
