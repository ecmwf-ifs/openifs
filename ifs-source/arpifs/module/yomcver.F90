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

MODULE YOMCVER

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK    ,DR_HOOK,  JPHOOK

USE YOMLUN   , ONLY : NULOUT   ,NULNAM   
USE YOMCT0   , ONLY : LR2D     ,LECMWF

IMPLICIT NONE

SAVE

! =============================================================================

TYPE TCVER
! ------ Vertical discretisation --------------------------------------------

! NDLNPR  : NDLNPR=0: conventional formulation of delta, i.e. ln(P(l)/P(l-1)).
!           NDLNPR=1: formulation of delta used in non hydrostatic model,
!                     i.e. (P(l)-P(l-1))/SQRT(P(l)*P(l-1)).
! RHYDR0  : value given to "alpha(1) = depth of log(Pi) between top and full level nr 1"
!           in case where general formula to compute "alpha" gives an infinite value
!           (used only if LVERTFE=F, NDLNPR=0).
!           This quantity is never used in the following cases:
!            LVERTFE=T.
!            LVERTFE=F with NDLNPR=1.
! LAPRXPK : way of computing full-levels pressures in primitive equation
!           hydrostatic model.
!           .T.: full levels are computed by PK=(PK+1/2 + PK-1/2)*0.5
!           .F.: full levels are computed by a more complicated formula
!                consistent with "alpha" in geopotential formula.

LOGICAL :: LAPRXPK
INTEGER(KIND=JPIM) :: NDLNPR
REAL(KIND=JPRB) :: RHYDR0

! ----- vertical discretisation, vertical boundaries:
! LREGETA   : .T.: for the interlayer L, ETA(L)=L/NFLEVG
!             .F.: for the interlayer L, ETA(L)=A(L)/P0+B(L)
! LVFE_REGETA: cf. LREGETA for "eta" used in VFE operators.
LOGICAL :: LREGETA
LOGICAL :: LVFE_REGETA


! * Variables related to vertical discretisation in finite elements:

! NVSCH         : type of basis if the finite element vertical discretisation is used.
!               (1=>linear functions, 3=>Hermite cubic functions)
! NVFE_TYPE     : Type of spline basis used for finite element vertical discretisation.
!               (1 = linear, 3 = cubic)
! NVFE_ORDER    : Order of spline used in VFE; NVFE_ORDER=NVFE_TYPE+1
! NVFE_INTBC/DERBC : Boundary conditions used for integrals/derivatives
!               0 - no bc applied (RINTE/RDERI/RDDERI used)
!               1 - old LVFE_INTB/DERIB (as in cy45)
!               2 - implicit definition through basis
!               3 - explicit definition
! NVFE_INTERNALS: number of internals knots
! NVFE_BC       : integer that determines the way the boundary knots are defined
!               0 - full levels used
!               1 - regular distribution between last internal knot and boundary

! LVERTFE       : .T./.F. Finite element/conventional vertical discretisation.
! LVFE_LAPL     : VFE for vertical Laplacian term (NH model)
! LVFE_LAPL_BC  : VFE for boundary cond. in vert. Laplacian term (NH model)
! LVFE_LAPL_TBC/BBC: VFE for top/bottom boundary cond. in vert. Laplacian term (NH model)
!               : if inner domain is purely in VFE manner
! LVFE_LAPL2PI  : simpler formula for vertical Laplacian in VFE
! RLAPL2PI      : parameter used for simpler Laplacian in VFE
! LVFE_X_TERM   : VFE X-term (NH model)
! LVFE_Z_TERM   : VFE Z-term (w on full levels in NH model)
! LVFE_GW       : VFE for vertical velocity (NH model); in this case
!                 vertical velocity is at full levels.
! LVFE_DELNHPRE : VFE to compute [Delta pre] at full levels.
! LVFE_GWMPA    : VFE for AROME physics vertical velocity
!                 (NH model with AROME physics)
! LVFE_CHEB     : chebyshev nodes (dense distribution of levels near BCs in eta space)
! LVFE_CENTRI   : Centripetal method for full levels eta_vfe calculation.
! RVFE_CENTRI   : Exponent in function computing eta_vfe if LVFE_CENTRI=T
! RVFE_ALPHA    : Exponent that constrols definition of eta. 
!                 RVFE_ALPHA =  0   -  gives regular (the same like LVFE_REGETA)
!                 RVFE_ALPHA =  1   -  gives classic sigma (eta = sigma for pressure VP00)
! RVFE_BETA     : Exponent that constrols density of levels close to boundaries.
!                 RVFE_BETA  =  0   -  there is density transformation
!                 RVFE_BETA  =  1   -  chanyshev definition (when combined with RVFE_ALPHA=0.0)
! LVFE_APPROX   : Approximation (or interpolation) used to represent a function.
! RVFE_KNOT_STRETCH   : stretching of knots 
! LVFE_ECMWF    : T if original ECMWF way to compute vertical integral and derivative
! LVFE_LAPL_HALF: Vertical Laplacian uses derivative operators full->half->full
! LVFE_FIX_ORDER: T/F - VFE operators defined with fixed order splines/fixed knot sequence
!                 (BCs are included by changing order of splines)
! LVFE_GW_HALF  : T - GW on HALF levels under key LGWADV 
! LVFE_MAXIMAS  : T/F - full levels at maximas of spline basis functions/full
!                 levels - Greville abscissa - Variation diminishing approach
! LVFE_VERBOSE  : print several diagnostics or not
! RMINDETA      : minimum distance between knots; for smaller intervals knots
!                 are considered to be equal (multiple knots)
! LDYN_ANALYSIS_STABILITY : this key is more general than VFE itself
!                 It turn on analysis of stability of linear operator under SUSI.
!                 We compute eigenvalues of mastrix "M = (I - tau L)^-1 (I + tau L)"
!                 M has dimension (2*NFLEVG + 1)*(2*NFLEVG + 1) and 
!                 correctly design L operator must have all eigenvalues Abs(eval) <= 1.0.
!                 (please move this key into NAMDYN and yomdyn. I did not done it to save compilation time .. lazyness:-))
! ----------------------------------------------------------------------

INTEGER(KIND=JPIM) :: NVSCH
INTEGER(KIND=JPIM) :: NVFE_TYPE
INTEGER(KIND=JPIM) :: NVFE_ORDER
INTEGER(KIND=JPIM) :: NVFE_INTBC
INTEGER(KIND=JPIM) :: NVFE_DERBC
INTEGER(KIND=JPIM) :: NVFE_INTERNALS
INTEGER(KIND=JPIM) :: NVFE_BC
LOGICAL :: LVERTFE
LOGICAL :: LVFE_LAPL
LOGICAL :: LVFE_LAPL_BC
LOGICAL :: LVFE_LAPL_TBC
LOGICAL :: LVFE_LAPL_BBC
LOGICAL :: LVFE_LAPL2PI
REAL(KIND=JPRB) :: RLAPL2PI = 0._JPRB
LOGICAL :: LVFE_X_TERM
LOGICAL :: LVFE_Z_TERM
LOGICAL :: LVFE_GW
LOGICAL :: LVFE_DELNHPRE
LOGICAL :: LVFE_GWMPA
LOGICAL :: LVFE_CENTRI
LOGICAL :: LVFE_CHEB
REAL(KIND=JPRB) :: RVFE_CENTRI
REAL(KIND=JPRB) :: RVFE_ALPHA
REAL(KIND=JPRB) :: RVFE_BETA
REAL(KIND=JPRB) :: RVFE_KNOT_STRETCH
LOGICAL :: LVFE_APPROX
LOGICAL :: LVFE_ECMWF
LOGICAL :: LVFE_LAPL_HALF
LOGICAL :: LVFE_FIX_ORDER
LOGICAL :: LVFE_GW_HALF
LOGICAL :: LVFE_MAXIMAS
LOGICAL :: LVFE_VERBOSE
LOGICAL :: LVFE_NORMALIZE
LOGICAL :: LDYN_ANALYSIS_STABILITY
REAL(KIND=JPRB) :: RMINDETA, RFAC1, RFAC2
END TYPE TCVER
! =============================================================================


CONTAINS

! =============================================================================
SUBROUTINE SUCVER_GEOM(YDCVER,LDNHDYN_GEOM)

!**** *SUCVER*   - Set-up for some keys
!                  used in the vertical finite elements discretisation

!     Purpose.
!     --------
!      sets-up YOMCVER

!**   Interface.
!     ----------
!        *CALL* *SUCVER(...)

!        Explicit arguments :
!        --------------------

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      K. Yessad (from some SUCT0 and SUDYN code).
!      Original : May 2012

! Modifications
! -------------
!      F. Vana  26-Sep-2019   Defaults changed to the new VFE scheme (LVFE_ECMWF=F)

! End Modifications
!      ----------------------------------------------------------------

IMPLICIT NONE

TYPE(TCVER),TARGET,INTENT(INOUT) :: YDCVER
LOGICAL,INTENT(IN)        :: LDNHDYN_GEOM
 
LOGICAL           ,POINTER :: LAPRXPK
REAL(KIND=JPRB)   ,POINTER :: RHYDR0
INTEGER(KIND=JPIM),POINTER :: NDLNPR
LOGICAL           ,POINTER :: LREGETA
LOGICAL           ,POINTER :: LVFE_REGETA
INTEGER(KIND=JPIM),POINTER :: NVSCH
INTEGER(KIND=JPIM),POINTER :: NVFE_TYPE
INTEGER(KIND=JPIM),POINTER :: NVFE_ORDER
INTEGER(KIND=JPIM),POINTER :: NVFE_INTBC
INTEGER(KIND=JPIM),POINTER :: NVFE_DERBC
INTEGER(KIND=JPIM),POINTER :: NVFE_INTERNALS
INTEGER(KIND=JPIM),POINTER :: NVFE_BC
LOGICAL,POINTER :: LVERTFE
LOGICAL,POINTER :: LVFE_LAPL
LOGICAL,POINTER :: LVFE_LAPL_BC
LOGICAL,POINTER :: LVFE_LAPL_TBC
LOGICAL,POINTER :: LVFE_LAPL_BBC
LOGICAL,POINTER :: LVFE_LAPL2PI
REAL(KIND=JPRB),POINTER :: RLAPL2PI
LOGICAL,POINTER :: LVFE_X_TERM
LOGICAL,POINTER :: LVFE_Z_TERM
LOGICAL,POINTER :: LVFE_GW
LOGICAL,POINTER :: LVFE_DELNHPRE
LOGICAL,POINTER :: LVFE_GWMPA
LOGICAL,POINTER :: LVFE_CENTRI
LOGICAL,POINTER :: LVFE_CHEB
REAL(KIND=JPRB),POINTER :: RVFE_CENTRI
REAL(KIND=JPRB),POINTER :: RVFE_ALPHA
REAL(KIND=JPRB),POINTER :: RVFE_BETA
REAL(KIND=JPRB),POINTER :: RVFE_KNOT_STRETCH
LOGICAL,POINTER :: LVFE_APPROX
LOGICAL,POINTER :: LVFE_ECMWF
LOGICAL,POINTER :: LVFE_LAPL_HALF
LOGICAL,POINTER :: LVFE_FIX_ORDER
LOGICAL,POINTER :: LVFE_GW_HALF
LOGICAL,POINTER :: LVFE_MAXIMAS
LOGICAL,POINTER :: LVFE_VERBOSE
LOGICAL,POINTER :: LVFE_NORMALIZE
LOGICAL,POINTER :: LDYN_ANALYSIS_STABILITY
REAL(KIND=JPRB),POINTER :: RMINDETA, RFAC1, RFAC2
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "namcver.nam.h"

! =============================================================================

#include "abor1.intfb.h"
#include "posnam.intfb.h"

!      ----------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('YOMCVER:SUCVER',0,ZHOOK_HANDLE)
!      ----------------------------------------------------------------

!*       1.    SET DEFAULT VALUES.
!    
!              -------------------
LAPRXPK=>YDCVER%LAPRXPK
RHYDR0=>YDCVER%RHYDR0
NDLNPR=>YDCVER%NDLNPR
LREGETA=>YDCVER%LREGETA
LVFE_REGETA=>YDCVER%LVFE_REGETA
NVSCH=>YDCVER%NVSCH
NVFE_TYPE=>YDCVER%NVFE_TYPE
NVFE_ORDER=>YDCVER%NVFE_ORDER
NVFE_INTBC=>YDCVER%NVFE_INTBC
NVFE_DERBC=>YDCVER%NVFE_DERBC
NVFE_INTERNALS=>YDCVER%NVFE_INTERNALS
NVFE_BC=>YDCVER%NVFE_BC
LVERTFE=>YDCVER%LVERTFE
LVFE_LAPL=>YDCVER%LVFE_LAPL
LVFE_LAPL_BC=>YDCVER%LVFE_LAPL_BC
LVFE_LAPL_TBC=>YDCVER%LVFE_LAPL_TBC
LVFE_LAPL_BBC=>YDCVER%LVFE_LAPL_BBC
LVFE_LAPL2PI=>YDCVER%LVFE_LAPL2PI
RLAPL2PI=>YDCVER%RLAPL2PI
LVFE_X_TERM=>YDCVER%LVFE_X_TERM
LVFE_Z_TERM=>YDCVER%LVFE_Z_TERM
LVFE_GW=>YDCVER%LVFE_GW
LVFE_DELNHPRE=>YDCVER%LVFE_DELNHPRE
LVFE_GWMPA=>YDCVER%LVFE_GWMPA
LVFE_CENTRI=>YDCVER%LVFE_CENTRI
LVFE_CHEB=>YDCVER%LVFE_CHEB
RVFE_CENTRI=>YDCVER%RVFE_CENTRI
RVFE_ALPHA=>YDCVER%RVFE_ALPHA
RVFE_BETA=>YDCVER%RVFE_BETA
RVFE_KNOT_STRETCH=>YDCVER%RVFE_KNOT_STRETCH
LVFE_APPROX=>YDCVER%LVFE_APPROX
LVFE_ECMWF=>YDCVER%LVFE_ECMWF
LVFE_LAPL_HALF=>YDCVER%LVFE_LAPL_HALF
LVFE_FIX_ORDER=>YDCVER%LVFE_FIX_ORDER
LVFE_GW_HALF=>YDCVER%LVFE_GW_HALF
LVFE_MAXIMAS=>YDCVER%LVFE_MAXIMAS
LVFE_VERBOSE=>YDCVER%LVFE_VERBOSE
LVFE_NORMALIZE=>YDCVER%LVFE_NORMALIZE
LDYN_ANALYSIS_STABILITY=>YDCVER%LDYN_ANALYSIS_STABILITY
RMINDETA=>YDCVER%RMINDETA
RFAC1=>YDCVER%RFAC1
RFAC2=>YDCVER%RFAC2


NDLNPR=0
LREGETA=.FALSE.
LVFE_REGETA=.FALSE.
IF (LECMWF) THEN
  LAPRXPK=.TRUE.
  RHYDR0=LOG(2._JPRB)
  IF(.NOT.LDNHDYN_GEOM .AND. .NOT.LR2D) THEN
    ! Assume semi-lagrangien to avoid dependency on model/MH
!    IF(LSLAG) THEN
      LVERTFE=.TRUE.
      NVFE_TYPE=3
!    ELSE
!      LVERTFE=.FALSE.
!      NVFE_TYPE=0
!    ENDIF
  ELSE
    LVERTFE=.FALSE.
    NVFE_TYPE=0
  ENDIF
ELSE
  LAPRXPK=.FALSE.
  RHYDR0=1._JPRB
  LVERTFE=.FALSE.
  NVFE_TYPE=0
ENDIF
NVFE_ORDER=NVFE_TYPE+1

NVSCH=0
NVFE_INTBC=3    ! 3 new VFE
NVFE_DERBC=3    ! 3 new VFE
NVFE_INTERNALS=0
NVFE_BC=0

LDYN_ANALYSIS_STABILITY = .FALSE.
LVFE_NORMALIZE = .FALSE.
LVFE_LAPL=LVERTFE.AND.LDNHDYN_GEOM
LVFE_LAPL_BC=LVERTFE.AND.LDNHDYN_GEOM
LVFE_LAPL_TBC=LVERTFE.AND.LDNHDYN_GEOM
LVFE_LAPL_BBC=LVERTFE.AND.LDNHDYN_GEOM
LVFE_LAPL2PI=LVERTFE.AND.LDNHDYN_GEOM
LVFE_X_TERM=LVERTFE.AND.LDNHDYN_GEOM
LVFE_GW=LVERTFE.AND.LDNHDYN_GEOM
LVFE_Z_TERM=LVERTFE.AND.LDNHDYN_GEOM
LVFE_DELNHPRE=LVERTFE.AND.LDNHDYN_GEOM
LVFE_GWMPA=.FALSE.
LVFE_CENTRI=.TRUE.  ! .TRUE. new VFE
LVFE_CHEB=.FALSE.
RVFE_ALPHA=0.0_JPRB
RVFE_BETA =0.5_JPRB
RVFE_KNOT_STRETCH=1.0_JPRB
LVFE_APPROX=.FALSE.
!!later LVFE_ECMWF=.NOT.LDNHDYN_GEOM
LVFE_ECMWF=.FALSE. ! .false. new VFE
LVFE_LAPL_HALF=.FALSE.
LVFE_FIX_ORDER=.TRUE.
LVFE_GW_HALF=.FALSE.
LVFE_MAXIMAS=.FALSE.
LVFE_VERBOSE=.FALSE.
RMINDETA=0.0_JPRB
RFAC1=0.0_JPRB
RFAC2=0.0_JPRB

!      ----------------------------------------------------------------
!*       2.    Modifies default values.
!              -----------------------

CALL POSNAM(NULNAM,'NAMCVER')
READ(NULNAM,NAMCVER)

!     ------------------------------------------------------------------

!*       3.    Reset variables and test.
!              -------------------------

! * (LAPRXPK,NDLNPR) reset to (T,0) if LVERTFE=T
IF (LVERTFE) THEN
  LAPRXPK=.TRUE.
  WRITE(UNIT=NULOUT,FMT='('' SUCVER_GEOM: VFE => LAPRXPK reset to TRUE '')')
  NDLNPR=0
  WRITE(UNIT=NULOUT,FMT='('' SUCVER_GEOM: VFE => NDLNPR reset to 0 '')')
ENDIF

IF (NVSCH/=0) THEN
  NVFE_TYPE=NVSCH
ENDIF
IF(.NOT.LVERTFE) THEN
  NVFE_TYPE=0
  NVFE_INTBC = 0
  NVFE_DERBC = 0
ENDIF
IF (LVFE_ECMWF) THEN
  NVFE_INTBC = 0
  NVFE_DERBC = 0
ENDIF
NVFE_ORDER=NVFE_TYPE+1

! * Reset the LVFE_... keys to F if LVERTFE=F.
LVFE_LAPL=(LVERTFE.AND.LDNHDYN_GEOM).AND.LVFE_LAPL
LVFE_LAPL_BC=(LVERTFE.AND.LDNHDYN_GEOM).AND.LVFE_LAPL_BC
LVFE_LAPL_TBC=(LVERTFE.AND.LDNHDYN_GEOM).AND.LVFE_LAPL_TBC
LVFE_LAPL_BBC=(LVERTFE.AND.LDNHDYN_GEOM).AND.LVFE_LAPL_BBC
LVFE_LAPL2PI=(LVERTFE.AND.LDNHDYN_GEOM).AND.LVFE_LAPL2PI
LVFE_X_TERM=LVERTFE.AND.LVFE_X_TERM
LVFE_Z_TERM=(LVERTFE.AND.LDNHDYN_GEOM).AND.LVFE_Z_TERM
LVFE_GW=LVERTFE.AND.LVFE_GW
LVFE_DELNHPRE=(LVERTFE.AND.LDNHDYN_GEOM).AND.LVFE_DELNHPRE
LVFE_GWMPA=LVERTFE.AND.LVFE_GWMPA
LVFE_CENTRI=LVERTFE.AND.LVFE_CENTRI
LVFE_CHEB=LVERTFE.AND.LVFE_CHEB
LVFE_APPROX=(LVERTFE.AND. .NOT.LVFE_ECMWF).AND.LVFE_APPROX
LVFE_LAPL_HALF=LVERTFE.AND.LVFE_LAPL_HALF
LVFE_FIX_ORDER=LVERTFE.AND.LVFE_FIX_ORDER
LVFE_GW_HALF=LVERTFE.AND.LVFE_GW_HALF
LVFE_MAXIMAS=LVERTFE.AND.LVFE_MAXIMAS
LVFE_VERBOSE=LVERTFE.AND.LVFE_VERBOSE

! * Reset LVFE_ECMWF to .T. if LVERTFE.AND.NVFE_TYPE==1.
IF (LVERTFE.AND.NVFE_TYPE==1) THEN
  LVFE_ECMWF=.TRUE.
  WRITE(NULOUT,*) ' SUCVER_GEOM: LVERTFE.AND.NVFE_TYPE=1 => LVFE_ECMWF set to TRUE'
ENDIF

! * Reset LVFE_ECMWF to .T. if .NOT.LDNHDYN_GEOM
IF (LVERTFE.AND..NOT.LDNHDYN_GEOM) THEN

  IF(LVFE_ECMWF)THEN
    WRITE(NULOUT,*) ' SUCVER_GEOM: LVERTFE in hydrostatic model => ECMWF operators'
  ELSE
    WRITE(NULOUT,*) ' SUCVER_GEOM: LVERTFE in hydrostatic model => LACE operators'
  ENDIF

ENDIF

! * Reset LVFE_CENTRI to .F. if LVERTFE.AND.LVFE_APPROX.
IF( LVERTFE.AND.LVFE_APPROX.AND.LVFE_CENTRI ) THEN
  ! * Approximation in VFE and centripetal definition of VFE_ETAH?
  LVFE_CENTRI=.FALSE.
  WRITE(NULOUT,*) ' SUCVER_GEOM: LVERTFE.AND.LVFE_APPROX => LVFE_CENTRI set to FALSE'
ENDIF

! * Bound RVFE_ALPHA between <-1, 1>
IF(LVFE_CENTRI)THEN
  RVFE_ALPHA=MAX(0._JPRB,MIN(1._JPRB,RVFE_ALPHA))
  RVFE_BETA =MAX(0._JPRB,MIN(1._JPRB,RVFE_BETA))
ENDIF

RVFE_KNOT_STRETCH=MAX(1._JPRB,RVFE_KNOT_STRETCH)

! * NVFE_TYPE must be >=1.
IF( LVERTFE .AND. NVFE_TYPE < 1 ) THEN
  CALL ABOR1(' SUCVER_GEOM: for VFE NVFE_TYPE must be >= 1')
ENDIF

! * NVFE_INTBC/DERBC must be <=3.
IF( LVERTFE .AND. NVFE_INTBC > 3 ) THEN
  CALL ABOR1(' SUCVER_GEOM: for VFE NVFE_INTBC must be <= 3')
ENDIF
IF( LVERTFE .AND. NVFE_DERBC > 3 ) THEN
  CALL ABOR1(' SUCVER_GEOM: for VFE NVFE_DERBC must be <= 3')
ENDIF

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('YOMCVER:SUCVER',1,ZHOOK_HANDLE)
END SUBROUTINE SUCVER_GEOM

SUBROUTINE PRT_CVER_GEOM(YDCVER)

!**** *PRT_CVER*   - Prints YOMCVER/NAMCVER keys 
!      ----------------------------------------------------------------
IMPLICIT NONE
TYPE(TCVER),INTENT(INOUT) :: YDCVER
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!      ----------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('YOMCVER:PRT_CVER_GEOM',0,ZHOOK_HANDLE)
!      ----------------------------------------------------------------

WRITE(UNIT=NULOUT,FMT='('' '')')
WRITE(UNIT=NULOUT,FMT='('' Printings of YOMCVER/NAMCVER variables '')')
WRITE(UNIT=NULOUT,FMT='('' LAPRXPK = '',L2,'' NDLNPR = '',I2,'' RHYDR0 = '',E10.4)') &
 & YDCVER%LAPRXPK,YDCVER%NDLNPR,YDCVER%RHYDR0
WRITE(UNIT=NULOUT,FMT='('' LREGETA = '',L2)') YDCVER%LREGETA
WRITE(UNIT=NULOUT,FMT='('' LVFE_REGETA = '',L2)') YDCVER%LVFE_REGETA
WRITE(UNIT=NULOUT,FMT='('' LVERTFE= '',L2,'' NVFE_TYPE= '',I2)') YDCVER%LVERTFE,YDCVER%NVFE_TYPE
WRITE(UNIT=NULOUT,FMT='('' NVFE_ORDER= '',I2)') YDCVER%NVFE_ORDER
WRITE(UNIT=NULOUT,FMT='('' NVFE_INTBC= '',I2)') YDCVER%NVFE_INTBC
WRITE(UNIT=NULOUT,FMT='('' NVFE_DERBC= '',I2)') YDCVER%NVFE_DERBC
WRITE(UNIT=NULOUT,FMT='('' LVFE_LAPL = '',L2)') YDCVER%LVFE_LAPL
WRITE(UNIT=NULOUT,FMT='('' LVFE_LAPL_BC = '',L2)') YDCVER%LVFE_LAPL_BC
WRITE(UNIT=NULOUT,FMT='('' LVFE_LAPL_TBC = '',L2)') YDCVER%LVFE_LAPL_TBC
WRITE(UNIT=NULOUT,FMT='('' LVFE_LAPL_BBC = '',L2)') YDCVER%LVFE_LAPL_BBC
WRITE(UNIT=NULOUT,FMT='('' LVFE_LAPL2PI = '',L2,'' RLAPL2PI ='',F20.14)') YDCVER%LVFE_LAPL2PI,YDCVER%RLAPL2PI
WRITE(UNIT=NULOUT,FMT='('' LVFE_X_TERM = '',L2)') YDCVER%LVFE_X_TERM
WRITE(UNIT=NULOUT,FMT='('' LVFE_Z_TERM = '',L2)') YDCVER%LVFE_Z_TERM
WRITE(UNIT=NULOUT,FMT='('' LVFE_GW = '',L2)') YDCVER%LVFE_GW
WRITE(UNIT=NULOUT,FMT='('' LVFE_DELNHPRE = '',L2)') YDCVER%LVFE_DELNHPRE
WRITE(UNIT=NULOUT,FMT='('' LVFE_GWMPA = '',L2)') YDCVER%LVFE_GWMPA
WRITE(UNIT=NULOUT,FMT='('' LVFE_CENTRI = '',L2,'' RVFE_ALPHA= '',F20.14,'' RVFE_BETA= '',F20.14)') &
 & YDCVER%LVFE_CENTRI, YDCVER%RVFE_ALPHA, YDCVER%RVFE_BETA
WRITE(UNIT=NULOUT,FMT='('' LVFE_CHEB = '',L2)') YDCVER%LVFE_CHEB
WRITE(UNIT=NULOUT,FMT='('' LVFE_APPROX = '',L2)') YDCVER%LVFE_APPROX
WRITE(UNIT=NULOUT,FMT='('' LVFE_ECMWF= '',L2)') YDCVER%LVFE_ECMWF
WRITE(UNIT=NULOUT,FMT='('' LVFE_LAPL_HALF = '',L2)') YDCVER%LVFE_LAPL_HALF
WRITE(UNIT=NULOUT,FMT='('' LVFE_FIX_ORDER = '',L2)') YDCVER%LVFE_FIX_ORDER
WRITE(UNIT=NULOUT,FMT='('' LVFE_GW_HALF = '',L2)') YDCVER%LVFE_GW_HALF
WRITE(UNIT=NULOUT,FMT='('' LVFE_MAXIMAS = '',L2)') YDCVER%LVFE_MAXIMAS
WRITE(UNIT=NULOUT,FMT='('' LVFE_VERBOSE = '',L2)') YDCVER%LVFE_VERBOSE
WRITE(UNIT=NULOUT,FMT='('' RMINDETA = '',F20.14)') YDCVER%RMINDETA
WRITE(UNIT=NULOUT,FMT='('' RVFE_KNOT_STRETCH = '',F20.14)') YDCVER%RVFE_KNOT_STRETCH
WRITE(UNIT=NULOUT,FMT='('' LVFE_NORMALIZE = '',L2)') YDCVER%LVFE_NORMALIZE
WRITE(UNIT=NULOUT,FMT='('' LDYN_ANALYSIS_STABILITY = '',L2)') YDCVER%LDYN_ANALYSIS_STABILITY
WRITE(UNIT=NULOUT,FMT='('' RFAC1 = '',F20.14)') YDCVER%RFAC1
WRITE(UNIT=NULOUT,FMT='('' RFAC2 = '',F20.14)') YDCVER%RFAC2

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('YOMCVER:PRT_CVER_GEOM',1,ZHOOK_HANDLE)
END SUBROUTINE PRT_CVER_GEOM

! =============================================================================
END MODULE YOMCVER
