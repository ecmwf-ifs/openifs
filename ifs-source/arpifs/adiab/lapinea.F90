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

#ifdef RS6K
@PROCESS NOCHECK
#endif
SUBROUTINE LAPINEA(&
 ! --- INPUT ------------------------------------------------------------------
 & YDGEOMETRY, YDML_GCONF,YDML_DYN,KST,KPROF,YDSL,KIBL,PB1,PB2,PWRL9,&
 ! --- INPUT/OUTPUT -----------------------------------------------------------
 & KVSEPC,KVSEPL,&
 ! --- OUTPUT -----------------------------------------------------------------
 & PSAVEDP,PCCO,PUF,PVF,KL0,KLH0,PLSCAW,PRSCAW,KL0H,PLSCAWH,PRSCAWH,&
 & PSCO,PGFLT1,KNOWENO)

!**** *LAPINEA* - semi-LAgrangian scheme(Trajectory):
!                Interface subroutine for interpolations and lagrangian
!                trends. (Programme INterface d'Ensemble).

!     Purpose.
!     --------
!           Grid point calculations in dynamics.

!           Computes the semi-Lagrangian trajectory. 
!           Computes the weights ready for interpolations at origin point.

!**   Interface.
!     ----------
!        *CALL* *LAPINEA(.....)

!        Explicit arguments :
!        --------------------
!        INPUT:
!          KST       - first element of arrays where computations are performed.
!          KPROF     - depth of work.
!          YDSL      - SL_STRUCT definition.
!          KIBL      - index into YRGSGEOM/YRCSGEOM instances in YDGEOMETRY
!          PB1       - SLBUF1-buffer for interpolations.
!          PB2       - SLBUF2-buf to communicate info from non lag. to lag. dyn.

!        INPUT AND OUTPUT:
!          KVSEPC    - vertical separation (used in S/L adjoint of cubic interp.)
!          KVSEPL    - vertical separation (used in S/L adjoint, linear interp.)

!        OUTPUT:

!          PCCO      - information about comput. space position of interpol. point.
!          PUF,PVF   - U-comp and V-comp of "(a/rs)*wind" necessary to
!                      find the position of the origin point,
!                      in a local repere linked to computational sphere.
!          KL0       - index of the four western points
!                      of the 16 points interpolation grid.
!          KLH0      - second value of index of the four western points
!                      of the 16 points interpolation grid if needed.
!          PLSCAW    - linear weights (distances) for interpolations.
!          PRSCAW    - non-linear weights for interpolations.
!          KL0H      - half-level KL0
!          PLSCAWH   - linear weights (distances) for interpolations at half-levels.
!          PRSCAWH   - non-linear weights for interpolations at half-levels.
!          PSCO      - information about geographic position of interpol. point.
!          PGFLT1    - t+dt unified_treatment grid-point fields (for SL diagn.)
!          KNOWENO   - quintic interpolation stencil indicator

!        Implicit arguments :
!        --------------------

!     Method.
!     -------

!     Externals.
!     ----------
!       See includes below.
!       Called by CALL_SL.

!     Reference.
!     ----------
!        ARPEGE documentation about semi-Lagrangian scheme.

!     Author.
!     -------
!        K. YESSAD after the subroutine CPLGDY2 written by
!        Maurice IMBARD   METEO FRANCE/EERM/CRMD
!        Original : JUNE 1991.

!     Modifications.
!     --------------
!   F.Vana  09-Jan-2007 new argument to LARCINA
!   K. Yessad 07-03-2007: Remove useless (gw)_surf interpol. in NH+LGWADV.
!   F. Vana   28-Aug-2007: removing 3D spline interpolation
!   N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!   30-Jun-2008 J. Masek    New dataflow for SLHD scheme.
!   K. Yessad Nov 2008: rationalisation of dummy argument interfaces
!   K. Yessad (Aug 2009): always use root (QX,QY) for (p,q) variables names
!   K. Yessad (Jan 2011): introduce INTDYN_MOD structures.
!   F. Vana 21-Feb-2011: Extended dimensions of weights for hor. turbulence
!   G. Mozdzynski (Jan 2011): OOPS cleaning, use of derived type SL_STRUCT
!   G.Mozdzynski (Feb 2011): OOPS cleaning, use of derived types TGSGEOM and TCSGEOM
!   G. Mozdzynski (May 2012): further cleaning
!   F. Vana 13-Feb-2014 KAPPAT for heat quantities.
!   T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!   S. Malardel (Nov 2013): COMAD weights for SL interpolations
!   K. Yessad (July 2014): Move some variables.
!   K. Yessad (March 2017): simplify level numbering in interpolator.
!   K. Yessad (June 2017): Introduce NHQE model.
!   F. Vana    21-Nov-2017: Options LSLDP_CURV and LHOISLT
!   F. Vana July 2018: RK4 scheme for trajectory research.
!   F. Vana October 2018: Extended LSLDP_CURV.
!   F. Vana 20-Fev-2019: Quintic vertical interpolation
!   R. El Khatib 27-02-2019 Use pointer function SC2PRG to avoid bounds violation
!   F. Vana 11-Jul-2019: Option LRHS_CURV
!   M. Diamantakis Feb 2021: Add cartesian geometry option LSLDP_XYZ for departure point iterations
! End Modifications
!-------------------------------------------------------------------------------

USE MODEL_DYNAMICS_MOD     , ONLY : MODEL_DYNAMICS_TYPE
USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE GEOMETRY_MOD           , ONLY : GEOMETRY
USE PARKIND1               , ONLY : JPIM, JPRB
USE YOMHOOK                , ONLY : LHOOK, DR_HOOK, JPHOOK
USE SC2PRG_MOD             , ONLY : SC2PRG
USE YOMCT0                 , ONLY : LRPLANE
USE YOMCST                 , ONLY : RPI
USE EINT_MOD               , ONLY : SL_STRUCT

!    -------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(MODEL_DYNAMICS_TYPE),INTENT(IN):: YDML_DYN
TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(IN):: YDML_GCONF
INTEGER(KIND=JPIM),INTENT(IN)    :: KST
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF
TYPE(SL_STRUCT)   ,INTENT(INOUT) :: YDSL
INTEGER(KIND=JPIM),INTENT(IN)    :: KIBL
REAL(KIND=JPRB)   ,INTENT(IN)    :: PB1(YDSL%NASLB1,YDML_DYN%YRPTRSLB1%NFLDSLB1)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PB2(YDGEOMETRY%YRDIM%NPROMA,YDML_DYN%YRPTRSLB2%NFLDSLB2)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWRL9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
INTEGER(KIND=JPIM),INTENT(INOUT) :: KVSEPC
INTEGER(KIND=JPIM),INTENT(INOUT) :: KVSEPL
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSAVEDP(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,3)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCCO(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YYTCCO%NDIM)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PUF(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVF(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
INTEGER(KIND=JPIM),INTENT(OUT)   :: KL0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,0:3)
INTEGER(KIND=JPIM),INTENT(OUT)   :: KLH0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,0:3)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PLSCAW(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YYTLSCAW%NDIM)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRSCAW(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YYTRSCAW%NDIM)
INTEGER(KIND=JPIM),INTENT(OUT)   :: KL0H(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,0:3)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PLSCAWH(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YYTLSCAWH%NDIM)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRSCAWH(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YYTRSCAWH%NDIM)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSCO(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YYTSCO%NDIM)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PGFLT1(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_GCONF%YGFL%NDIM1)
INTEGER(KIND=JPIM),INTENT(OUT)   :: KNOWENO(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: ILEV(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
INTEGER(KIND=JPIM) :: IGLGLO
INTEGER(KIND=JPIM) :: IHVI
INTEGER(KIND=JPIM) :: IROT
INTEGER(KIND=JPIM) :: ISTALAT(YDSL%NDGSAH:YDSL%NDGENH)
INTEGER(KIND=JPIM) :: ITIP
INTEGER(KIND=JPIM) :: JGL
INTEGER(KIND=JPIM) :: JROF
INTEGER(KIND=JPIM) :: JLEV
INTEGER(KIND=JPIM) :: JGFL
INTEGER(KIND=JPIM) :: IDEP(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)

LOGICAL :: LLINTV

REAL(KIND=JPRB) :: ZUF0(1), ZURL0(1), ZVF0(1), ZZRL0(1)
REAL(KIND=JPRB) :: ZVRL0(1), ZZF0(1), ZWF0(1), ZWFSM(1), ZWRL0(1)

! SL diagnostics
REAL(KIND=JPRB) :: ZWF(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)

REAL(KIND=JPRB) :: ZDEPI                     ! Now only used by ELARMES
REAL(KIND=JPRB) :: ZGMDTX(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB) :: ZGMDTY(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB) :: ZLEV(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZLSDEPI(YDSL%NDGSAH:YDSL%NDGENH)

REAL(KIND=JPRB), POINTER :: ZSLB2KAPPA(:), ZSLB2KAPPAH(:), ZSLB2KAPPAM(:), ZSLB2KAPPAT(:)
REAL(KIND=JPRB), POINTER :: ZSLB2STDDISU(:), ZSLB2STDDISV(:), ZSLB2STDDISW(:)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "elarmes.intfb.h"
#include "larcina.intfb.h"
#include "larcinha.intfb.h"
#include "larmes.intfb.h"
#include "larmes_rk.intfb.h"
#include "larmes_xyz.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LAPINEA',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, &
 &  YDCSGEOM=>YDGEOMETRY%YRCSGEOM(KIBL), &
 & YDGSGEOM=>YDGEOMETRY%YRGSGEOM(KIBL), YDEGEO=>YDGEOMETRY%YREGEO, YDDYN=>YDML_DYN%YRDYN, &
 & YDDYNA=>YDML_DYN%YRDYNA, &
 & YDPTRSLB1=>YDML_DYN%YRPTRSLB1, YDPTRSLB2=>YDML_DYN%YRPTRSLB2,YDRIP=>YDML_GCONF%YRRIP, &
 & YGFL=>YDML_GCONF%YGFL, YDTCCO=>YDML_DYN%YYTCCO, YDTSCO=>YDML_DYN%YYTSCO)

ASSOCIATE(NDIM1=>YGFL%NDIM1, NSLDIA=>YGFL%NSLDIA, NSLDIAGP=>YGFL%NSLDIAGP, &
 & NUMFLDS=>YGFL%NUMFLDS, YCOMP=>YGFL%YCOMP, YSLDIA=>YGFL%YSLDIA, &
 & NPROMA=>YDDIM%NPROMA, &
 & NFLEN=>YDDIMV%NFLEN, NFLEVG=>YDDIMV%NFLEVG, NFLSA=>YDDIMV%NFLSA, &
 & LADVF=>YDDYN%LADVF, LFINDVSEP=>YDDYN%LFINDVSEP, & 
 & LSLDP_RK=>YDDYN%LSLDP_RK, LSLDP_XYZ=>YDDYN%LSLDP_XYZ, LSLDP_SAVE=>YDDYN%LSLDP_SAVE, & 
 & LRHS_CURV=>YDDYN%LRHS_CURV, NCURRENT_ITER=>YDDYN%NCURRENT_ITER, NSITER=>YDDYN%NSITER, &
 & EDELX=>YDEGEO%EDELX, EDELY=>YDEGEO%EDELY, &
 & NFLDSLB1=>YDPTRSLB1%NFLDSLB1, &
 & MSLB2KAPPA=>YDPTRSLB2%MSLB2KAPPA, MSLB2KAPPAH=>YDPTRSLB2%MSLB2KAPPAH, &
 & MSLB2KAPPAM=>YDPTRSLB2%MSLB2KAPPAM, MSLB2KAPPAT=>YDPTRSLB2%MSLB2KAPPAT, &
 & MSLB2STDDISU=>YDPTRSLB2%MSLB2STDDISU, MSLB2STDDISV=>YDPTRSLB2%MSLB2STDDISV, &
 & MSLB2STDDISW=>YDPTRSLB2%MSLB2STDDISW, NFLDSLB2=>YDPTRSLB2%NFLDSLB2, &
 & RTDT=>YDRIP%RTDT)

CALL SC2PRG(MSLB2KAPPA,PB2      ,ZSLB2KAPPA)
CALL SC2PRG(MSLB2KAPPAH,PB2      ,ZSLB2KAPPAH)
CALL SC2PRG(MSLB2KAPPAM,PB2      ,ZSLB2KAPPAM)
CALL SC2PRG(MSLB2KAPPAT,PB2      ,ZSLB2KAPPAT)
CALL SC2PRG(MSLB2STDDISU,PB2      ,ZSLB2STDDISU)
CALL SC2PRG(MSLB2STDDISV,PB2      ,ZSLB2STDDISV)
CALL SC2PRG(MSLB2STDDISW,PB2      ,ZSLB2STDDISW)

!     ------------------------------------------------------------------

!*       1.    PRELIMINARY INITIALISATIONS:
!              ---------------------------

! The following are now only used by ELARMES:
ZDEPI=2.0_JPRB*RPI

DO JGL=MAX(YDSL%NDGSAG,YDSL%NDGSAL+YDSL%NFRSTLOFF-YDSL%NSLWIDE)-YDSL%NFRSTLOFF,&
   & MIN(YDSL%NDGENG,YDSL%NDGENL+YDSL%NFRSTLOFF+YDSL%NSLWIDE)-YDSL%NFRSTLOFF  
  IGLGLO=JGL+YDSL%NFRSTLOFF
  ZLSDEPI(JGL)=REAL(YDSL%NLOENG(IGLGLO),JPRB)/ZDEPI
ENDDO

IHVI=0
DO JGFL=1,NUMFLDS
  IF(YCOMP(JGFL)%CSLINT == 'LAIHVT      ' .OR.&
     & YCOMP(JGFL)%CSLINT == 'LAIHVTQM    ' .OR.&
     & YCOMP(JGFL)%CSLINT == 'LAIHVTQMH   ' ) THEN  
    IHVI=1
  ENDIF
ENDDO

IF(LRPLANE) THEN
  DO JROF = KST, KPROF
    ZGMDTX (JROF)=0.5_JPRB*YDGSGEOM%GM(JROF)*RTDT/EDELX
    ZGMDTY (JROF)=0.5_JPRB*YDGSGEOM%GM(JROF)*RTDT/EDELY
  ENDDO
ENDIF

DO JGL=MAX(YDSL%NDGSAH,LBOUND(YDSL%NSLOFF,1)),MIN(YDSL%NDGENH,UBOUND(YDSL%NSLOFF,1))
  ISTALAT(JGL)=YDSL%NSLOFF(JGL)
ENDDO

!     ------------------------------------------------------------------

!*       2.    COMPUTATION OF THE SL TRAJECTORY.
!              ---------------------------------

IF (LRPLANE) THEN

  IF (LSLDP_RK) THEN
    CALL ABOR1(" LAPINEA : RK4 scheme is not yet available for LAM. ")
    !CALL ELARMES_RK(...)
  ELSEIF (LSLDP_XYZ) THEN
    CALL ABOR1(" LAPINEA : DP scheme in XYZ coordinates is not yet available for LAM. ")
    !CALL ELARMES_XYZ(...)
  ELSE
  CALL ELARMES(YDGEOMETRY,YDML_DYN,YDRIP,KST,KPROF,YDSL,ISTALAT,ZGMDTX,ZGMDTY,PB1,PB2,&
   & ZLSDEPI,KIBL,&
   & KVSEPC,KVSEPL,PSCO,ZLEV,PCCO,PUF,PVF,&
   & KL0,KLH0,ILEV,PLSCAW,PRSCAW)  
  ENDIF

ELSE

  IF (LSLDP_RK) THEN
    CALL LARMES_RK(YDGEOMETRY,YDRIP,YDML_DYN,KST,KPROF,YDSL,ISTALAT,PB1,PB2,&
     & ZLSDEPI,KIBL,&
     & KVSEPC,KVSEPL,PSCO,ZLEV,PCCO,PUF,PVF,ZWF,&
     & KL0,KLH0,ILEV,PLSCAW,PRSCAW)
  ELSEIF (LSLDP_XYZ) THEN
    CALL LARMES_XYZ(YDGEOMETRY,YDRIP,YDML_DYN,KST,KPROF,&
      & YDSL,ISTALAT,PB1,PB2,PWRL9,ZLSDEPI,KIBL,&
      & KVSEPC,KVSEPL,PSAVEDP,PSCO,ZLEV,PCCO,PUF,PVF,ZWF,&
      & KL0,KLH0,ILEV)
  ELSE
    CALL LARMES(YDGEOMETRY,YDRIP,YDML_DYN,KST,KPROF,YDSL,ISTALAT,PB1,PB2,&
     & ZLSDEPI,KIBL,&
     & KVSEPC,KVSEPL,PSCO,ZLEV,PCCO,PUF,PVF,ZWF,&
     & KL0,KLH0,ILEV,PLSCAW,PRSCAW)
  ENDIF

  ! diagnostics for semi-Lagrangian scheme
  IF( YDDYNA%LSLDIA .AND. ( YDDYNA%LPC_CHEAP.OR.(NCURRENT_ITER == NSITER) )) THEN
    DO JGFL=NSLDIAGP+1,NSLDIA
      IF (TRIM(YSLDIA(JGFL)%CNAME) == 'UTRAJ           ') THEN
        PGFLT1(KST:KPROF,1:NFLEVG,YSLDIA(JGFL)%MP1)=PUF(KST:KPROF,1:NFLEVG)
      ENDIF
      IF (TRIM(YSLDIA(JGFL)%CNAME) == 'VTRAJ           ') THEN     
        PGFLT1(KST:KPROF,1:NFLEVG,YSLDIA(JGFL)%MP1)=PVF(KST:KPROF,1:NFLEVG)
      ENDIF
      IF (TRIM(YSLDIA(JGFL)%CNAME) == 'WTRAJ           ') THEN     
        PGFLT1(KST:KPROF,1:NFLEVG,YSLDIA(JGFL)%MP1)=ZWF(KST:KPROF,1:NFLEVG)
      ENDIF
      IF (TRIM(YSLDIA(JGFL)%CNAME) == 'LAMBDA0         ') THEN
        DO JLEV=1,NFLEVG
          DO JROF=KST,KPROF     
            PGFLT1(JROF,JLEV,YSLDIA(JGFL)%MP1)=SIN(PCCO(JROF,JLEV,YDTCCO%M_RLON))
          ENDDO
        ENDDO
      ENDIF
      IF (TRIM(YSLDIA(JGFL)%CNAME) == 'THETA0          ') THEN     
        DO JLEV=1,NFLEVG
          DO JROF=KST,KPROF     
            PGFLT1(JROF,JLEV,YSLDIA(JGFL)%MP1)=SIN(PCCO(JROF,JLEV,YDTCCO%M_RLAT))
          ENDDO
        ENDDO
      ENDIF
      IF (TRIM(YSLDIA(JGFL)%CNAME) == 'ETA0            ') THEN  
        PGFLT1(KST:KPROF,1:NFLEVG,YSLDIA(JGFL)%MP1)=ZLEV(KST:KPROF,1:NFLEVG)
      ENDIF
    ENDDO
  ENDIF

ENDIF

!     ------------------------------------------------------------------

!*       3.    COMPUTATION OF WEIGHTS AND INDICES ARRAYS FOR INTERPOLATIONS.
!              -------------------------------------------------------------

!*       3.1   Full-level Quantities.
!              ---------------------

IF (LRHS_CURV) THEN
  IROT=0
ELSE
  IROT=1
ENDIF
ITIP=2

LLINTV=.FALSE.

!  Warning: ZURL0,ZVRL0,ZZRL0,ZWRL0,ZUF0,ZVF0,ZZF0,ZWF0,ZWFSM are not
!           used in this call to LARCINA.

CALL LARCINA(YDGEOMETRY,YDML_DYN,KST,KPROF,YDSL,IHVI,ISTALAT,LFINDVSEP,YDDYNA%LSLHD,YDDYNA%LSLHDQUAD,LLINTV,&
 & ITIP,IROT,LRPLANE,ZLSDEPI,&
 & KIBL,&
 & PSCO,ZLEV,&
 & ZSLB2KAPPA,ZSLB2KAPPAT,ZSLB2KAPPAM,ZSLB2KAPPAH,&
 & ZSLB2STDDISU,ZSLB2STDDISV,ZSLB2STDDISW,&
 & ZURL0,ZVRL0,ZZRL0,ZWRL0,&
 & KVSEPC,KVSEPL,PCCO,&
 & ZUF0,ZVF0,ZZF0,ZWF0,ZWFSM,&
 & KL0,KLH0,ILEV,PLSCAW,PRSCAW,IDEP,KNOWENO)

IF (LRPLANE) THEN
! * Use PSCO (SINLA and COPHI) to transfer origin point coordinates to LAPINEB
!   for "2 Omega vectorial a k" recalculation.
  IF(LADVF) THEN
    PSCO(KST:KPROF,1:NFLEVG,YDTSCO%M_SINLA)=PCCO(KST:KPROF,1:NFLEVG,YDTCCO%M_RLON)
    PSCO(KST:KPROF,1:NFLEVG,YDTSCO%M_COPHI)=PCCO(KST:KPROF,1:NFLEVG,YDTCCO%M_RLAT)
  ENDIF
ENDIF

!*       3.2   Half-level Quantities.
!              ---------------------

IF (YDDYNA%LNHDYN.AND.YDDYNA%LGWADV.AND.(.NOT.YDGEOMETRY%YRCVER%LVFE_GW)) THEN

  IROT=0
  ITIP=2

  CALL LARCINHA(YDGEOMETRY,YDML_DYN,KST,KPROF,YDSL,IHVI,ISTALAT,ITIP,IROT,ZLSDEPI,&
   & KIBL,&
   & PSCO,PCCO,ZLEV,ZSLB2KAPPA,ZSLB2KAPPAT,ZSLB2KAPPAM,ZSLB2KAPPAH,&
   & ZSLB2STDDISU,ZSLB2STDDISV,ZSLB2STDDISW,&
   & KL0H,PLSCAWH,PRSCAWH)

ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('LAPINEA',1,ZHOOK_HANDLE)
END SUBROUTINE LAPINEA
