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

SUBROUTINE LASSIE(YDGEOMETRY,YDGMV,YGFL,YDDYN,YDDYNA,KST,KPROF,PRCORI,PGMV,PGMVS,PGFL,&
 & PSDIV0,PSDIV9,PTOD0,PTOD9,PGAGT0L,PGAGT0M,PGAGT9L,PGAGT9M)  

!**** *LASSIE*   Semi-Lagrangian scheme.
!                Computation of linear terms used in the semi-implicit scheme. 

!     Purpose.
!     --------
!        Computation of linear terms used in the semi-implicit scheme:
!        Nabla(Gamma*T+Mu*Pi), (Tau*D) and (Nu*D).

!**   Interface.
!     ----------
!        *CALL* *LASSIE(..)

!        Explicit arguments :
!        --------------------

!        INPUT:
!          KST      - first element of work.
!          KPROF    - depth of work.
!          PRCORI   - Coriolis parameter 2.OMEGA .
!          PGMV     - GMV variables at t-dt and t.
!          PGMVS    - GMVS variables at t-dt and t.
!          PGFL     - unified_treatment grid-point (GFL) fields.

!        OUTPUT:
!          PSDIV0   - semi-implicit term at time t for continuity equation
!                     (Nu*D).
!          PSDIV9   - semi-implicit term at time t-dt for continuity equation
!                     (Nu*D).
!          PTOD0    - semi-implicit term at time t for temperature equation
!                     (Tau*D).
!          PTOD9    - semi-implicit term at time t-dt for temperature equation
!                     (Tau*D).
!          PGAGT0L  - semi-implicit term at time t for U-wind equation
!                     (zonal component of Nabla(Gamma*T+Mu*Pi)).
!          PGAGT0M  - semi-implicit term at time t for V-wind equation
!                     (meridian component of Nabla(Gamma*T+Mu*Pi)).
!          PGAGT9L  - semi-implicit term at time t-dt for U-wind equation
!                     (zonal component of Nabla(Gamma*T+Mu*Pi)).
!          PGAGT9M  - semi-implicit term at time t-dt for V-wind equation
!                     (meridian component of Nabla(Gamma*T+Mu*Pi)).

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!           Called by LACDYN.

!     Reference.
!     ----------
!             Arpege documentation about semi-lagrangian scheme.

!     Author.
!     -------
!      K. YESSAD (METEO FRANCE/CNRM/GMAP) after part 3.1 of LACDYN.
!      Loops have been recoded according to F90 norms.
!      Original : AUGUST 1995.

!     Modifications.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      Modified 04-02-06: Y. Seity : new arguments to GPRCP (rain, snow and graupel)
!      M.Hamrud  15-Jan-2006  Revised GPRCP
!      K. Yessad Aug 2008: rationalisation of dummy argument interfaces
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMGMV       , ONLY : TGMV
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST       , ONLY : RD
USE YOMDYN       , ONLY : TDYN
USE YOMDYNA      , ONLY : TDYNA
USE YOM_YGFL     , ONLY : TYPE_GFLD

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TGMV)        ,INTENT(INOUT) :: YDGMV
TYPE(TDYN)        ,INTENT(IN)    :: YDDYN
TYPE(TDYNA)       ,INTENT(IN)    :: YDDYNA
TYPE(TYPE_GFLD)   ,INTENT(IN)    :: YGFL
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRCORI(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGMV(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDGMV%NDIMGMV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGMVS(YDGEOMETRY%YRDIM%NPROMA,YDGMV%NDIMGMVS)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGFL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YGFL%NDIM)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSDIV0(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSDIV9(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTOD0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTOD9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGAGT0L(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGAGT0M(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGAGT9L(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGAGT9M(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
!     ------------------------------------------------------------------
REAL(KIND=JPRB) :: ZR9  (YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)

INTEGER(KIND=JPIM) :: IPROFS, JLEV, JROF
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "gprcp.intfb.h"
#include "sigam.intfb.h"
#include "sitnu.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LASSIE',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP,  &
 & YDVAB=>YDGEOMETRY%YRVAB, YDVETA=>YDGEOMETRY%YRVETA, YDVFE=>YDGEOMETRY%YRVFE, YDSTA=>YDGEOMETRY%YRSTA,  &
 & YDLAP=>YDGEOMETRY%YRLAP, YDCSGLEG=>YDGEOMETRY%YRCSGLEG,  &
 & YDCSGEOM=>YDGEOMETRY%YRCSGEOM, &
  & YDCSGEOM_NB=>YDGEOMETRY%YRCSGEOM_NB, YDGSGEOM=>YDGEOMETRY%YRGSGEOM, YDGSGEOM_NB=>YDGEOMETRY%YRGSGEOM_NB,  &
  & YDSPGEOM=>YDGEOMETRY%YSPGEOM)
ASSOCIATE(NDIM=>YGFL%NDIM, &
 & NPROMA=>YDDIM%NPROMA, &
 & NFLEVG=>YDDIMV%NFLEVG, &
 & LIMPF=>YDDYN%LIMPF, &
 & NDIMGMV=>YDGMV%NDIMGMV, NDIMGMVS=>YDGMV%NDIMGMVS, YT0=>YDGMV%YT0, &
 & YT9=>YDGMV%YT9)
!     ------------------------------------------------------------------

IPROFS=KPROF-KST+1

! * Variables at time 0 (SL3TL and SL2TL).

!   - Computation of Nu*D (SI term for continuity equation)
!     and Tau*D (SI term for temperature equation).
CALL SITNU(YDGEOMETRY,YDDYN,NPROMA,1,PGMV(KST,1,YT0%MDIV),PTOD0(KST,1),PSDIV0(KST),IPROFS)
!   - Computation of Nabla(Gamma*T+Mu*Pi) (SI term for momentum equation).
CALL SIGAM(YDGEOMETRY,YDDYN,NPROMA,1,PGAGT0L(KST,1),PGMV(KST,1,YT0%MTL),PGMVS(KST,YT0%MSPL),&
 & IPROFS,NFLEVG)
CALL SIGAM(YDGEOMETRY,YDDYN,NPROMA,1,PGAGT0M(KST,1),PGMV(KST,1,YT0%MTM),PGMVS(KST,YT0%MSPM),&
 & IPROFS,NFLEVG)

! * Variables at time 9 (SL3TL only).

IF (.NOT.YDDYNA%LTWOTL) THEN
  ! - Computation of Nu*D (SI term for continuity equation)
  !   and Tau*D (SI term for temperature equation).
  CALL SITNU(YDGEOMETRY,YDDYN,NPROMA,1,PGMV(KST,1,YT9%MDIV),PTOD9(KST,1),PSDIV9(KST),IPROFS)
  ! - Computation of Nabla(Gamma*T+Mu*Pi) (SI term for momentum equation).
  CALL SIGAM(YDGEOMETRY,YDDYN,NPROMA,1,PGAGT9L(KST,1),PGMV(KST,1,YT9%MTL),PGMVS(KST,YT9%MSPL),&
   & IPROFS,NFLEVG)
  CALL SIGAM(YDGEOMETRY,YDDYN,NPROMA,1,PGAGT9M(KST,1),PGMV(KST,1,YT9%MTM),PGMVS(KST,YT9%MSPM),&
   & IPROFS,NFLEVG)
ENDIF

! * For "spectral RT" option, adjust semi-implicit term (Tau*D) in
!   T-equation to compensate for later multiplication by R/Rd

IF (YDDYNA%LSPRT) THEN
  IF (YDDYNA%LTWOTL) THEN
    ! remark for lpc_full:
    !  predictor: treatment of "t" data.
    !  corrector: treatment of provisional "t+dt" data.
    ! So in this case this is always the Y[X]%MP data which are used.
    CALL GPRCP(NPROMA,KST,KPROF,NFLEVG,PGFL=PGFL,PR=ZR9)  
    DO JLEV=1,NFLEVG
      DO JROF=KST,KPROF
        PTOD0(JROF,JLEV)=RD*PTOD0(JROF,JLEV)/ZR9(JROF,JLEV)
      ENDDO
    ENDDO
  ELSE
    CALL GPRCP(NPROMA,KST,KPROF,NFLEVG,PGFL=PGFL,KGFLTYP=9,PR=ZR9)  
    DO JLEV=1,NFLEVG
      DO JROF=KST,KPROF
        PTOD0(JROF,JLEV)=RD*PTOD0(JROF,JLEV)/ZR9(JROF,JLEV)
        PTOD9(JROF,JLEV)=RD*PTOD9(JROF,JLEV)/ZR9(JROF,JLEV)
      ENDDO
    ENDDO
  ENDIF
ENDIF

! * Add semi-implicit Coriolis terms to Nabla(Gamma*T+Mu*Pi)
!   if required (LIMPF=.T.).

IF (LIMPF) THEN
  DO JLEV=1,NFLEVG
    DO JROF=KST,KPROF
      PGAGT0L(JROF,JLEV)=PGAGT0L(JROF,JLEV)-PRCORI(JROF)*PGMV(JROF,JLEV,YT0%MV)
      PGAGT0M(JROF,JLEV)=PGAGT0M(JROF,JLEV)+PRCORI(JROF)*PGMV(JROF,JLEV,YT0%MU)
    ENDDO
  ENDDO
  IF (.NOT.YDDYNA%LTWOTL) THEN
    DO JLEV=1,NFLEVG
      DO JROF=KST,KPROF
        PGAGT9L(JROF,JLEV)=PGAGT9L(JROF,JLEV)&
         & -PRCORI(JROF)*PGMV(JROF,JLEV,YT9%MV)
        PGAGT9M(JROF,JLEV)=PGAGT9M(JROF,JLEV)&
         & +PRCORI(JROF)*PGMV(JROF,JLEV,YT9%MU)
      ENDDO
    ENDDO
  ENDIF
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('LASSIE',1,ZHOOK_HANDLE)
END SUBROUTINE LASSIE

