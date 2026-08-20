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
USE YOMCT0   , ONLY : LR2D     ,LECMWF   , LARPEGEF

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

! NVFE_TYPE     : Type of spline basis used for finite element vertical discretisation.
!               (1 = linear, 3 = cubic)
! NVFE_ORDER    : Order of spline used in VFE; NVFE_ORDER=NVFE_TYPE+1
! NVFE_INTERNALS: number of internals knots

! LVERTFE       : .T./.F. Finite element/conventional vertical discretisation.
! LVFE_LAPL_BC  : VFE for boundary cond. in vert. Laplacian term (NH model)

! VFE for vertical velocity (NH model):
! LVFE_GW       : T - invertible RINTGW/RDERGW used under key LGWADV, full levels gw.
! LVFE_GW_HALF  : T - invertible RINTGW/RDERGW used under key LGWADV, half levels gw.
! LVFE_GWMPA    : T - VFE for AROME physics vertical velocity

! LVFE_CHEB     : chebyshev nodes (dense distribution of levels near BCs in eta space)
! LVFE_ECMWF    : T if original ECMWF way to compute vertical integral and derivative
! LVFE_LAPL_HALF: Vertical Laplacian uses derivative operators full->half->full
!               : this is consistent with algorithm applied infull nonlinear model
! LVFE_NOBC     : T no boundary conditions applied for vert.derivative (RDERI)
!                 F boundary conditions applied for vert.derivative (RDERB)
! LPERCENTILS   : used in SUVFE_KNOT to determine method for computing knots for basis functions
! LVFE_VERBOSE  : print several diagnostics or not
! CVFE_ETAH     : half levels eta definition
!           REGETA - regular distribution
!           MCONST - general definition
!           MCNTRL - explicitly prescribed m = dpi/eta at boundaries; density control
!           CHORDAL - using A and Bs
! LDYN_ANALYSIS_STABILITY : this key is more general than VFE itself
!                 It turn on analysis of stability of linear operator under SUSI.
!                 We compute eigenvalues of mastrix "M = (I - tau L)^-1 (I + tau L)"
!                 M has dimension (2*NFLEVG + 1)*(2*NFLEVG + 1) and 
!                 correctly design L operator must have all eigenvalues Abs(eval) <= 1.0.
!                 (please move this key into NAMDYN and yomdyn. I did not done it to save compilation time .. lazyness:-))
! RVFE_ALPHA    : Exponent that constrols definition of eta. 
!                 RVFE_ALPHA =  0   -  gives regular (the same like LVFE_REGETA)
!                 RVFE_ALPHA =  1   -  gives classic sigma (eta = sigma for pressure VP00)
! RVFE_BETA     : Exponent that constrols density of levels close to boundaries.
!                 RVFE_BETA  =  0   -  there is density transformation
!                 RVFE_BETA  =  1   -  chanyshev definition (when combined with RVFE_ALPHA=0.0)
! RVFE_ALPHA/RVFE_BETA : control definition of eta close to boundaries
!           RVFE_ALPHA/BETA =  0   -  gives regular (the same like LVFE_REGETA)
!           RVFE_ALPHA/BETA >  0   -  denser close to boundaries
!           RVFE_ALPHA/BETA <  0   -  denser close to inner domain (MCNTRL only)
! RVFE_KNOT_STRETCH   : stretching of knots 
! ----------------------------------------------------------------------

INTEGER(KIND=JPIM) :: NVFE_TYPE
INTEGER(KIND=JPIM) :: NVFE_ORDER
INTEGER(KIND=JPIM) :: NVFE_INTERNALS
LOGICAL :: LVERTFE
LOGICAL :: LVFE_LAPL_BC
LOGICAL :: LVFE_GW
LOGICAL :: LVFE_GW_HALF
LOGICAL :: LVFE_GWMPA
LOGICAL :: LVFE_CHEB
LOGICAL :: LVFE_ECMWF
LOGICAL :: LVFE_LAPL_HALF
LOGICAL :: LVFE_NOBC
LOGICAL :: LVFE_VERBOSE
LOGICAL :: LVFE_NORMALIZE
LOGICAL :: LDYN_ANALYSIS_STABILITY
LOGICAL :: LPERCENTILS
REAL(KIND=JPRB) :: RVFE_ALPHA
REAL(KIND=JPRB) :: RVFE_BETA
REAL(KIND=JPRB) :: RVFE_KNOT_STRETCH
REAL(KIND=JPRB) :: RFAC1, RFAC2
CHARACTER(LEN=8) :: CVFE_ETAH
REAL(KIND=JPRB)  :: REXP_VRAT
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
!      P. Smolikova (Sep 2020): VFE pruning.

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
INTEGER(KIND=JPIM),POINTER :: NVFE_TYPE
INTEGER(KIND=JPIM),POINTER :: NVFE_ORDER
INTEGER(KIND=JPIM),POINTER :: NVFE_INTERNALS
LOGICAL,POINTER :: LVERTFE
LOGICAL,POINTER :: LVFE_LAPL_BC
LOGICAL,POINTER :: LVFE_GW
LOGICAL,POINTER :: LVFE_GW_HALF
LOGICAL,POINTER :: LVFE_GWMPA
LOGICAL,POINTER :: LVFE_CHEB
LOGICAL,POINTER :: LPERCENTILS
REAL(KIND=JPRB),POINTER :: RVFE_ALPHA
REAL(KIND=JPRB),POINTER :: RVFE_BETA
REAL(KIND=JPRB),POINTER :: RVFE_KNOT_STRETCH
LOGICAL,POINTER :: LVFE_ECMWF
LOGICAL,POINTER :: LVFE_LAPL_HALF
LOGICAL,POINTER :: LVFE_NOBC
LOGICAL,POINTER :: LVFE_VERBOSE
LOGICAL,POINTER :: LVFE_NORMALIZE
LOGICAL,POINTER :: LDYN_ANALYSIS_STABILITY
REAL(KIND=JPRB),POINTER :: RFAC1, RFAC2
CHARACTER(LEN=8),POINTER :: CVFE_ETAH
REAL(KIND=JPRB),POINTER  :: REXP_VRAT

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "namcver.nam.h"

! =============================================================================

#include "abor1.intfb.h"
#include "posnam.intfb.h"

!      ----------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('YOMCVER:SUCVER_GEOM',0,ZHOOK_HANDLE)
!      ----------------------------------------------------------------

!*       1.    SET DEFAULT VALUES.
!
!              -------------------
LAPRXPK=>YDCVER%LAPRXPK
RHYDR0=>YDCVER%RHYDR0
NDLNPR=>YDCVER%NDLNPR
LREGETA=>YDCVER%LREGETA
LVFE_REGETA=>YDCVER%LVFE_REGETA
NVFE_TYPE=>YDCVER%NVFE_TYPE
NVFE_ORDER=>YDCVER%NVFE_ORDER
NVFE_INTERNALS=>YDCVER%NVFE_INTERNALS
LVERTFE=>YDCVER%LVERTFE
LVFE_LAPL_BC=>YDCVER%LVFE_LAPL_BC
LVFE_GW=>YDCVER%LVFE_GW
LVFE_GW_HALF=>YDCVER%LVFE_GW_HALF
LVFE_GWMPA=>YDCVER%LVFE_GWMPA
LVFE_CHEB=>YDCVER%LVFE_CHEB
LPERCENTILS=>YDCVER%LPERCENTILS
RVFE_ALPHA=>YDCVER%RVFE_ALPHA
RVFE_BETA=>YDCVER%RVFE_BETA
RVFE_KNOT_STRETCH=>YDCVER%RVFE_KNOT_STRETCH
LVFE_ECMWF=>YDCVER%LVFE_ECMWF
LVFE_LAPL_HALF=>YDCVER%LVFE_LAPL_HALF
LVFE_NOBC=>YDCVER%LVFE_NOBC
LVFE_VERBOSE=>YDCVER%LVFE_VERBOSE
LVFE_NORMALIZE=>YDCVER%LVFE_NORMALIZE
LDYN_ANALYSIS_STABILITY=>YDCVER%LDYN_ANALYSIS_STABILITY
RFAC1=>YDCVER%RFAC1
RFAC2=>YDCVER%RFAC2
CVFE_ETAH=>YDCVER%CVFE_ETAH
REXP_VRAT=>YDCVER%REXP_VRAT

NDLNPR=0
LREGETA=.FALSE.
LVFE_REGETA=.FALSE.

IF (LECMWF .AND. .NOT. LARPEGEF) THEN
  LAPRXPK=.TRUE.
  RHYDR0=LOG(2._JPRB)
  IF(.NOT.LDNHDYN_GEOM .AND. .NOT.LR2D) THEN
    ! Assume semi-lagrangian to avoid dependency on model/MH
      LVERTFE=.TRUE.
      NVFE_TYPE=3
  ELSE
    LVERTFE=.TRUE.
    NVFE_TYPE=3
  ENDIF
  LPERCENTILS=.TRUE.
  RVFE_ALPHA=0.0_JPRB
  RVFE_BETA =0.5_JPRB
  CVFE_ETAH='MCONST'
ELSE
  LAPRXPK=.FALSE.
  RHYDR0=1._JPRB
  LVERTFE=.FALSE.
  NVFE_TYPE=0
  LPERCENTILS=.FALSE.
  RVFE_ALPHA=0.0_JPRB
  RVFE_BETA =0.0_JPRB
  CVFE_ETAH='CHORDAL'
ENDIF

NVFE_ORDER=NVFE_TYPE+1

NVFE_INTERNALS=0
LDYN_ANALYSIS_STABILITY = .FALSE.
LVFE_NORMALIZE = .FALSE.
LVFE_LAPL_BC=.FALSE.
LVFE_GW=.FALSE.

LVFE_CHEB=.FALSE.
RVFE_KNOT_STRETCH=1.0_JPRB
LVFE_ECMWF=.NOT.(LECMWF.OR.LDNHDYN_GEOM)
LVFE_LAPL_HALF=.TRUE.
LVFE_GW_HALF=.FALSE.
LVFE_GWMPA=.FALSE.
LVFE_NOBC=.FALSE.
LVFE_VERBOSE=.FALSE.
RFAC1=0.0_JPRB
RFAC2=0.0_JPRB
REXP_VRAT=1.0_JPRB

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

IF(.NOT.LVERTFE) THEN
  NVFE_TYPE=0
ENDIF
NVFE_ORDER=NVFE_TYPE+1

! * Reset the LVFE_... keys to F if LVERTFE=F.
LVFE_LAPL_BC=(LVERTFE.AND.LDNHDYN_GEOM).AND.LVFE_LAPL_BC
LVFE_GW=LVERTFE.AND.LVFE_GW
LVFE_CHEB=LVERTFE.AND.LVFE_CHEB
LVFE_GW_HALF=LVERTFE.AND.LVFE_GW_HALF
LVFE_GWMPA=LVERTFE.AND.LVFE_GWMPA
LVFE_NOBC=LVERTFE.AND.LVFE_NOBC
LVFE_VERBOSE=LVERTFE.AND.LVFE_VERBOSE
LVFE_LAPL_HALF=LVERTFE.AND.LVFE_LAPL_HALF

IF(LVFE_LAPL_HALF)THEN
  LVFE_LAPL_BC = .TRUE.
  WRITE(NULOUT,*) 'SUCVER_GEOM: LVFE_LAPL_BC set to TRUE to be consistent with LVFE_LAPL_HALF'
ENDIF

! * Reset LVFE_ECMWF to .F. if LNHDYN
IF (LVERTFE.AND.LDNHDYN_GEOM) THEN
  IF (NVFE_TYPE==1) THEN
    CALL ABOR1(' SUCVER_GEOM: LVERTFE.AND.NVFE_TYPE=1 in NH model not possible.')
  ENDIF
  IF (LVFE_ECMWF) THEN
    WRITE(NULOUT,*) ' SUCVER_GEOM: LVERTFE in NH model => LACE operators; LVFE_ECMWF set to FALSE.'
    LVFE_ECMWF=.FALSE.
  ELSE
    WRITE(NULOUT,*) ' SUCVER_GEOM: LVERTFE in NH model => LACE operators.'
  ENDIF
ENDIF

! * Reset LVFE_ECMWF to .T. if .NOT.LDNHDYN_GEOM
IF (LVERTFE.AND..NOT.LDNHDYN_GEOM) THEN
  IF(LVFE_ECMWF)THEN
    WRITE(NULOUT,*) ' SUCVER_GEOM: LVERTFE in hydrostatic model => ECMWF operators'
  ELSE
    WRITE(NULOUT,*) ' SUCVER_GEOM: LVERTFE in hydrostatic model => LACE operators'
  ENDIF
ENDIF

! Only ECMWF operators available without boundary conditions.
IF (LVFE_NOBC.AND..NOT.LVFE_ECMWF) THEN
  LVFE_NOBC=.FALSE.
  WRITE(NULOUT,*) ' SUCVER_GEOM: LACE operators in VFE => LVFE_NOBC set to FALSE.'
ENDIF

! * Only one from LVFE_GW and LVFE_GW_HALF may be true.
IF( LVERTFE .AND. LVFE_GW .AND. LVFE_GW_HALF ) THEN
  CALL ABOR1(' SUCVER_GEOM: LVFE_GW and LVFE_GW_HALF can not be both set to TRUE')
ENDIF

! * Bound RVFE_ALPHA between <0, 1>
IF (CVFE_ETAH=='MCONST') THEN
  RVFE_ALPHA=MAX(0._JPRB,MIN(1._JPRB,RVFE_ALPHA))
  RVFE_BETA =MAX(0._JPRB,MIN(1._JPRB,RVFE_BETA))
  WRITE(NULOUT,*) ' SUCVER_GEOM: CVFE_ETAH=MCONST, 0< RVFE_ALPHA/BETA <1.'
ELSEIF (.NOT.(CVFE_ETAH == "CHORDAL".OR.CVFE_ETAH == "REGETA".OR.CVFE_ETAH == "MCNTRL")) THEN
  WRITE(NULOUT,*) " SUCVER_GEOM: CVFE_ETAH '",CVFE_ETAH,"'"
  CALL ABOR1(" SUCVER: unknown option for CVFE_ETAH")
ENDIF

RVFE_KNOT_STRETCH=MAX(1._JPRB,RVFE_KNOT_STRETCH)

! * NVFE_TYPE must be >=1.
IF( LVERTFE .AND. NVFE_TYPE < 1 ) THEN
  CALL ABOR1(' SUCVER_GEOM: for VFE NVFE_TYPE must be >= 1')
ENDIF

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('YOMCVER:SUCVER_GEOM',1,ZHOOK_HANDLE)
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
WRITE(UNIT=NULOUT,FMT='('' CVFE_ETAH= '',A)') YDCVER%CVFE_ETAH
WRITE(UNIT=NULOUT,FMT='('' RVFE_ALPHA= '',F20.14,'' RVFE_BETA= '',F20.14)') YDCVER%RVFE_ALPHA, YDCVER%RVFE_BETA
WRITE(UNIT=NULOUT,FMT='('' LPERCENTILS= '',L2)') YDCVER%LPERCENTILS
WRITE(UNIT=NULOUT,FMT='('' NVFE_ORDER= '',I2)') YDCVER%NVFE_ORDER
WRITE(UNIT=NULOUT,FMT='('' LVFE_LAPL_BC = '',L2)') YDCVER%LVFE_LAPL_BC
WRITE(UNIT=NULOUT,FMT='('' LVFE_GW = '',L2)') YDCVER%LVFE_GW
WRITE(UNIT=NULOUT,FMT='('' LVFE_GW_HALF = '',L2)') YDCVER%LVFE_GW_HALF
WRITE(UNIT=NULOUT,FMT='('' LVFE_GWMPA = '',L2)') YDCVER%LVFE_GWMPA
WRITE(UNIT=NULOUT,FMT='('' LVFE_ECMWF = '',L2)') YDCVER%LVFE_ECMWF
WRITE(UNIT=NULOUT,FMT='('' LVFE_LAPL_HALF = '',L2)') YDCVER%LVFE_LAPL_HALF
WRITE(UNIT=NULOUT,FMT='('' LVFE_NOBC= '',L2)') YDCVER%LVFE_NOBC
WRITE(UNIT=NULOUT,FMT='('' LVFE_CHEB = '',L2)') YDCVER%LVFE_CHEB
WRITE(UNIT=NULOUT,FMT='('' LVFE_VERBOSE = '',L2)') YDCVER%LVFE_VERBOSE
WRITE(UNIT=NULOUT,FMT='('' RVFE_KNOT_STRETCH = '',F20.14)') YDCVER%RVFE_KNOT_STRETCH
WRITE(UNIT=NULOUT,FMT='('' LVFE_NORMALIZE = '',L2)') YDCVER%LVFE_NORMALIZE
WRITE(UNIT=NULOUT,FMT='('' LDYN_ANALYSIS_STABILITY = '',L2)') YDCVER%LDYN_ANALYSIS_STABILITY
WRITE(UNIT=NULOUT,FMT='('' RFAC1 = '',F20.14)') YDCVER%RFAC1
WRITE(UNIT=NULOUT,FMT='('' RFAC2 = '',F20.14)') YDCVER%RFAC2
WRITE(UNIT=NULOUT,FMT='('' REXP_VRAT = '',F20.14)') YDCVER%REXP_VRAT

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('YOMCVER:PRT_CVER_GEOM',1,ZHOOK_HANDLE)
END SUBROUTINE PRT_CVER_GEOM

! =============================================================================
END MODULE YOMCVER
