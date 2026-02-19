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
 & PSGRTL,PSGRTM,&
 & PGMV,PTNDU,PTNDV,PTNDU_NOC,PTNDV_NOC)

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
! End Modifications
!------------------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMGMV       , ONLY : TGMV
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCT0       , ONLY : LRPLANE
USE YOMDYN       , ONLY : TDYN
USE YOEPHY       , ONLY : TEPHY

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TGMV)        ,INTENT(INOUT) :: YDGMV
TYPE(TDYN)        ,INTENT(IN)    :: YDDYN
TYPE(TEPHY)       ,INTENT(IN)    :: YDEPHY
INTEGER(KIND=JPIM),INTENT(IN)    :: KST
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRCORI(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGNORDL(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGNORDM(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSGRTL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSGRTM(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGMV(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDGMV%NDIMGMV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTNDU(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTNDV(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTNDU_NOC(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTNDV_NOC(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
!     ------------------------------------------------------------------
INTEGER(KIND=JPIM) ::JROF, JLEV
REAL(KIND=JPRB) :: ZKRFU(YDGEOMETRY%YRDIMV%NFLEVG),ZKRFV(YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZTNDU_RF(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG), &
 & ZTNDV_RF(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GP_TNDLAGADIAB_UV',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NPROMA=>YDDIM%NPROMA, &
 & NFLEVG=>YDDIMV%NFLEVG, &
 & RKRF=>YDDYN%RKRF, &
 & LEGWWMS=>YDEPHY%LEGWWMS, &
 & NSTTYP=>YDGEM%NSTTYP, &
 & NDIMGMV=>YDGMV%NDIMGMV, YT0=>YDGMV%YT0)
!     ------------------------------------------------------------------

! * tendency due to Rayleigh friction:
IF ((YDDYN%LRFRIC .AND. .NOT.LEGWWMS).AND.YDDYN%LRFRICISOTR) THEN
  ! isotropic Rayleigh friction:
  DO JLEV=1,NFLEVG
    DO JROF=KST,KPROF
      ZTNDU_RF(JROF,JLEV)=-RKRF(JLEV)*PGMV(JROF,JLEV,YT0%MU)
      ZTNDV_RF(JROF,JLEV)=-RKRF(JLEV)*PGMV(JROF,JLEV,YT0%MV)
    ENDDO
  ENDDO
ELSEIF ((YDDYN%LRFRIC .AND. .NOT.LEGWWMS).AND..NOT.YDDYN%LRFRICISOTR) THEN
  ! non-isotropic Rayleigh friction:
  IF (.NOT.LRPLANE .AND. NSTTYP==1) THEN
    ! not tilted spherical geometry (gnordl=0, gnordm=1 everywhere):
    DO JLEV=1,NFLEVG
      DO JROF=KST,KPROF
        ZTNDU_RF(JROF,JLEV)=-RKRF(JLEV)*PGMV(JROF,JLEV,YT0%MU)
        ZTNDV_RF(JROF,JLEV)=0._JPRB
      ENDDO
    ENDDO
  ELSE
    ! tilted spherical geometry or plane projection:
    ZKRFU(1:NFLEVG)=RKRF(1:NFLEVG)
    ZKRFV(1:NFLEVG)=0._JPRB
    ! what follows is also valid with any value of ZKRFV different from ZKRFU:
    DO JLEV=1,NFLEVG
      DO JROF=KST,KPROF
        ZTNDU_RF(JROF,JLEV)=&
         & -(ZKRFU(JLEV)*PGNORDM(JROF)*PGNORDM(JROF)+ZKRFV(JLEV)*PGNORDL(JROF)*PGNORDL(JROF))*PGMV(JROF,JLEV,YT0%MU)&
         & + (ZKRFU(JLEV)-ZKRFV(JLEV))*PGNORDM(JROF)*PGNORDL(JROF)*PGMV(JROF,JLEV,YT0%MV)
        ZTNDV_RF(JROF,JLEV)=&
         & (ZKRFU(JLEV)-ZKRFV(JLEV))*PGNORDM(JROF)*PGNORDL(JROF)*PGMV(JROF,JLEV,YT0%MU)&
         & -(ZKRFU(JLEV)*PGNORDL(JROF)*PGNORDL(JROF)+ZKRFV(JLEV)*PGNORDM(JROF)*PGNORDM(JROF))*PGMV(JROF,JLEV,YT0%MV)
      ENDDO
    ENDDO
  ENDIF
ELSE
  ! no Rayleigh friction:
  ZTNDU_RF(KST:KPROF,1:NFLEVG)=0._JPRB
  ZTNDV_RF(KST:KPROF,1:NFLEVG)=0._JPRB
ENDIF

! * total tendency:
DO JLEV=1,NFLEVG
  DO JROF=KST,KPROF
    PTNDU_NOC(JROF,JLEV)=-PSGRTL(JROF,JLEV)+ZTNDU_RF(JROF,JLEV)
    PTNDU(JROF,JLEV)=PTNDU_NOC(JROF,JLEV)+PRCORI(JROF)*PGMV(JROF,JLEV,YT0%MV)
    PTNDV_NOC(JROF,JLEV)=-PSGRTM(JROF,JLEV)+ZTNDV_RF(JROF,JLEV)
    PTNDV(JROF,JLEV)=PTNDV_NOC(JROF,JLEV)-PRCORI(JROF)*PGMV(JROF,JLEV,YT0%MU)
  ENDDO
ENDDO

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GP_TNDLAGADIAB_UV',1,ZHOOK_HANDLE)
END SUBROUTINE GP_TNDLAGADIAB_UV

