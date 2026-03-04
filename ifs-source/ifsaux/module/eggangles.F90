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

MODULE EGGANGLES

! Version 2009.0317 by JD GRIL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! DOC !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  All these functions make a package tool for angle.
!  In functions where appears DOM and UNIT indicate the domain of validity :
!         DOM     UNIT           Longitudes           Latitudes
!         "-+"    "D"          [-180.0,180.0[       [-90.0,90.0]
!         "0+"    "D"             [0.0,360.0[       [-90.0,90.0]
!         "-+"    "R"             [-pi,pi[        [-pi/2.0,pi/2.0] 
!         "0+"    "R"              [0,pi[         [-pi/2.0,pi/2.0]
!  (defaults values are DOM = "-+" and UNIT = "D").

!  All functions work for scalar or one dimensional array in input.

!  -1- ANGLE_DOMAIN function

!->function ANGLE_DOMAIN(ALPHA,PI,DOM,UNIT)

!     Converts longitudes in UNIT values under choisen DOMain.
!     The input (ALPHA) is a longitude (REAL) or a LOLA type structure ( or
!     array of them). The output has the same type than the input.       

!  -2- VAL_ functions

!->integer function VAL_LAT(LAT,NUM_ERR,PI,UNIT)

!     Test validity of LAT [-90.0,90.0] 
!     Return -1 or NUM_ERR if it's present in error case, 1 if it's ok.

!->integer function VAL_LON(LON,NUM_ERR,PI,DOM,UNIT)

!     Test validity of LON [-180.0,180.0[ or [0.0,360.0[ 
!     Return -1 or NUM_ERR if it's present in error case, 1 if it's ok.

!->integer function VAL_COORD(PT_COORD,NUM_ERR,PI,DOM,UNIT) 

!     Test validity of LAT [-90.0,90.0] and LON [-180.0,180.0[ or [0.0,360.0[
!     (depends the value of DOM) of a PT_COORD structure of type LOLA (in UNIT).
!     Return -1 or NUM_ERR if it's present in error case, 1 if it's ok.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Author : Jean-Daniel GRIL , CNRM/GMAP/COOPE , Februry 08 2000

! Modified:
! In April 2001 by M. Janousek (A few modifs to port the deck to the model code)
! In November 2004 by JD Gril : more routines to manage angles
!                             : debug VAL_COORD_x
! 2005 by JD Gril : more functions for Mercator RT
! In June 2006 by JD Gril     : line too long (L607 > 132 col.)
! July 2008 by JD Gril        : add 2 new functions to compute distance or size
!                             : on longitude : DIST_2REF and SIZE_W2E
! March 2009 by JD Gril       : add Vector routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! ******************* Definition of parameters **********************************

! Include Kinds
! -------------

!* kindef: define default KIND macros
! --------------------------------------
USE PARKIND1  ,ONLY : JPIM,    JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
! --------------------------------------

IMPLICIT NONE

! ******************* Definition of type ****************************************

TYPE LOLA
  SEQUENCE
  REAL(KIND=JPRB) :: LON, LAT
END TYPE LOLA

! ******************* Definition of Interface ***********************************

INTERFACE ANGLE_DOMAIN
  MODULE PROCEDURE ANGLE_DOMAIN_RS, ANGLE_DOMAIN_LOLAS, ANGLE_DOMAIN_RV, ANGLE_DOMAIN_LOLAV
END INTERFACE
INTERFACE VAL_LAT
  MODULE PROCEDURE VAL_LAT_S, VAL_LAT_V
END INTERFACE
INTERFACE VAL_LON
  MODULE PROCEDURE VAL_LON_S, VAL_LON_V
END INTERFACE
INTERFACE VAL_COORD
  MODULE PROCEDURE VAL_COORD_S, VAL_COORD_V
END INTERFACE
INTERFACE LOLAD
  MODULE PROCEDURE LOLAD_S, LOLAD_V
END INTERFACE
INTERFACE LOLAR
  MODULE PROCEDURE LOLAR_S, LOLAR_V
END INTERFACE
INTERFACE MINIMAX
  MODULE PROCEDURE  MINIMAX_S, MINIMAX_V
END INTERFACE
INTERFACE COSIN_TO_ANGLE
  MODULE PROCEDURE COSIN_TO_ANGLE_S, COSIN_TO_ANGLE_V
END INTERFACE
INTERFACE P_ASIN
  MODULE PROCEDURE P_ASIN_S, P_ASIN_V
END INTERFACE
INTERFACE P_ACOS
  MODULE PROCEDURE P_ACOS_S, P_ACOS_V
END INTERFACE
INTERFACE DIST_2REF
  MODULE PROCEDURE DIST_2REF_S, DIST_2REF_V, DIST_2REF_L
END INTERFACE
INTERFACE SIZE_W2E
  MODULE PROCEDURE SIZE_W2E_S, SIZE_W2E_L
END INTERFACE
CONTAINS

! =================== FUNCTIONS =================================================

! ******************* Independants functions ************************************

! -------------------------------------------------------------------------------
REAL(KIND=JPRB) FUNCTION ANGLE_DOMAIN_RS(ALPHA,PI,DOM,UNIT) RESULT (BETA)
REAL(KIND=JPRB), INTENT(IN)                           :: ALPHA
CHARACTER (LEN=2), INTENT(IN), OPTIONAL               :: DOM
CHARACTER (LEN=1), INTENT(IN), OPTIONAL               :: UNIT
REAL(KIND=JPRB), INTENT(IN), OPTIONAL                 :: PI


REAL(KIND=JPRB)   :: CVT, TPI, M
CHARACTER (LEN=2) :: TDOM
CHARACTER (LEN=1) :: TUNIT
REAL(KIND=JPHOOK)   :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('EGGANGLES:ANGLE_DOMAIN_RS',0,ZHOOK_HANDLE)
IF (PRESENT(PI)) THEN
  TPI = PI
ELSE
  TPI = ASIN(1.0_JPRB)*2.0_JPRB
ENDIF
IF (PRESENT(DOM)) THEN
  IF ((DOM=='0+').OR.(DOM=='-+')) THEN
    TDOM = DOM
  ELSE
    TDOM = "-+"
  ENDIF
ELSE
  TDOM = "-+"
ENDIF
IF (PRESENT(UNIT)) THEN
  IF ((UNIT=='R').OR.(UNIT=='D')) THEN
    TUNIT = UNIT
  ELSE
    TUNIT = "D"
  ENDIF
ELSE
  TUNIT = "D"
ENDIF

IF (TUNIT=='R') THEN
  CVT = TPI
ELSE
  CVT = 180.0_JPRB
ENDIF

IF (TDOM=='-+') THEN
  M = MOD(ALPHA,CVT)
  BETA = (M-CVT*MOD(REAL(INT(ALPHA/CVT),KIND=JPRB),2.0_JPRB))*SIGN(1.0_JPRB,ALPHA)*SIGN(1.0_JPRB,M)
ELSE
  M = MOD(ALPHA,2.0_JPRB*CVT)
  BETA = M-2.0_JPRB*CVT*(SIGN(0.5_JPRB,ALPHA)-0.5_JPRB)
ENDIF
IF (LHOOK) CALL DR_HOOK('EGGANGLES:ANGLE_DOMAIN_RS',1,ZHOOK_HANDLE)
END FUNCTION ANGLE_DOMAIN_RS
! -------------------------------------------------------------------------------
TYPE (LOLA) FUNCTION ANGLE_DOMAIN_LOLAS(ALPHA,PI,DOM,UNIT) RESULT (BETA)
TYPE (LOLA), INTENT(IN)                                :: ALPHA
CHARACTER (LEN=2), INTENT(IN), OPTIONAL                :: DOM
CHARACTER (LEN=1), INTENT(IN), OPTIONAL                :: UNIT
REAL(KIND=JPRB), INTENT(IN), OPTIONAL                  :: PI

REAL(KIND=JPRB)   :: TPI
CHARACTER (LEN=2) :: TDOM
CHARACTER (LEN=1) :: TUNIT
REAL(KIND=JPHOOK)   :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('EGGANGLES:ANGLE_DOMAIN_LOLAS',0,ZHOOK_HANDLE)
IF (PRESENT(PI)) THEN
  TPI = PI
ELSE
  TPI = ASIN(1.0_JPRB)*2.0_JPRB
ENDIF
IF (PRESENT(DOM)) THEN
  IF ((DOM=='0+').OR.(DOM=='-+')) THEN
    TDOM = DOM
  ELSE
    TDOM = "-+"
  ENDIF
ELSE
  TDOM = "-+"
ENDIF
IF (PRESENT(UNIT)) THEN
  IF ((UNIT=='R').OR.(UNIT=='D')) THEN
    TUNIT = UNIT
  ELSE
    TUNIT = "D"
  ENDIF
ELSE
  TUNIT = "D"
ENDIF

BETA%LON = ANGLE_DOMAIN(ALPHA%LON,TPI,TDOM,TUNIT)
BETA%LAT = ALPHA%LAT
IF (LHOOK) CALL DR_HOOK('EGGANGLES:ANGLE_DOMAIN_LOLAS',1,ZHOOK_HANDLE)
END FUNCTION ANGLE_DOMAIN_LOLAS
! -------------------------------------------------------------------------------
FUNCTION ANGLE_DOMAIN_RV(ALPHA,PI,DOM,UNIT) RESULT (BETA)
REAL(KIND=JPRB), DIMENSION(:), INTENT(IN)             :: ALPHA
CHARACTER (LEN=2), INTENT(IN), OPTIONAL               :: DOM
CHARACTER (LEN=1), INTENT(IN), OPTIONAL               :: UNIT
REAL(KIND=JPRB), INTENT(IN), OPTIONAL                 :: PI
REAL(KIND=JPRB), DIMENSION(SIZE(ALPHA)) :: BETA

REAL(KIND=JPRB)                         :: CVT, TPI
REAL(KIND=JPRB), DIMENSION(SIZE(ALPHA)) :: Z_M
CHARACTER (LEN=2)                       :: TDOM
CHARACTER (LEN=1)                       :: TUNIT
REAL(KIND=JPHOOK)                         :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('EGGANGLES:ANGLE_DOMAIN_RV',0,ZHOOK_HANDLE)
IF (PRESENT(PI)) THEN
  TPI = PI
ELSE
  TPI = ASIN(1.0_JPRB)*2.0_JPRB
ENDIF
IF (PRESENT(DOM)) THEN
  IF ((DOM=='0+').OR.(DOM=='-+')) THEN
    TDOM = DOM
  ELSE
    TDOM = "-+"
  ENDIF
ELSE
  TDOM = "-+"
ENDIF
IF (PRESENT(UNIT)) THEN
  IF ((UNIT=='R').OR.(UNIT=='D')) THEN
    TUNIT = UNIT
  ELSE
    TUNIT = "D"
  ENDIF
ELSE
  TUNIT = "D"
ENDIF

IF (TUNIT=='R') THEN
  CVT = TPI
ELSE
  CVT = 180.0_JPRB
ENDIF

IF (TDOM=='-+') THEN
  Z_M(:) = MOD(ALPHA(:),CVT)
  BETA = (Z_M(:)-CVT*MOD(REAL(INT(ALPHA(:)/CVT),KIND=JPRB),2.0_JPRB))*SIGN(1.0_JPRB,ALPHA(:))*SIGN(1.0_JPRB,Z_M(:))
ELSE
  Z_M(:) = MOD(ALPHA(:),2.0_JPRB*CVT)
  BETA = Z_M(:)-2.0_JPRB*CVT*(SIGN(0.5_JPRB,ALPHA(:))-0.5_JPRB)
ENDIF
IF (LHOOK) CALL DR_HOOK('EGGANGLES:ANGLE_DOMAIN_RV',1,ZHOOK_HANDLE)
END FUNCTION ANGLE_DOMAIN_RV
! -------------------------------------------------------------------------------
FUNCTION ANGLE_DOMAIN_LOLAV(YL_ALPHA,PI,DOM,UNIT) RESULT (YD_BETA)
TYPE (LOLA), DIMENSION(:), INTENT(IN)               :: YL_ALPHA
CHARACTER (LEN=2), INTENT(IN), OPTIONAL             :: DOM
CHARACTER (LEN=1), INTENT(IN), OPTIONAL             :: UNIT
REAL(KIND=JPRB), INTENT(IN), OPTIONAL               :: PI
TYPE (LOLA), DIMENSION(SIZE(YL_ALPHA)) :: YD_BETA

REAL(KIND=JPRB)   :: TPI
CHARACTER (LEN=2) :: TDOM
CHARACTER (LEN=1) :: TUNIT
REAL(KIND=JPHOOK)   :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('EGGANGLES:ANGLE_DOMAIN_LOLAV',0,ZHOOK_HANDLE)
IF (PRESENT(PI)) THEN
  TPI = PI
ELSE
  TPI = ASIN(1.0_JPRB)*2.0_JPRB
ENDIF
IF (PRESENT(DOM)) THEN
  IF ((DOM=='0+').OR.(DOM=='-+')) THEN
    TDOM = DOM
  ELSE
    TDOM = "-+"
  ENDIF
ELSE
  TDOM = "-+"
ENDIF
IF (PRESENT(UNIT)) THEN
  IF ((UNIT=='R').OR.(UNIT=='D')) THEN
    TUNIT = UNIT
  ELSE
    TUNIT = "D"
  ENDIF
ELSE
  TUNIT = "D"
ENDIF

YD_BETA(:)%LON = ANGLE_DOMAIN(YL_ALPHA(:)%LON,TPI,TDOM,TUNIT)
YD_BETA(:)%LAT = YL_ALPHA(:)%LAT
IF (LHOOK) CALL DR_HOOK('EGGANGLES:ANGLE_DOMAIN_LOLAV',1,ZHOOK_HANDLE)
END FUNCTION ANGLE_DOMAIN_LOLAV
! -------------------------------------------------------------------------------
INTEGER(KIND=JPIM) FUNCTION VAL_LAT_S(LAT,NUM_ERR,PI,UNIT) RESULT(ETAT)
REAL(KIND=JPRB), INTENT(IN)                          :: LAT
CHARACTER (LEN=1), INTENT(IN), OPTIONAL              :: UNIT
REAL(KIND=JPRB), INTENT(IN), OPTIONAL                :: PI
INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL             :: NUM_ERR

INTEGER(KIND=JPIM) :: TNE
REAL(KIND=JPRB)    :: TPI, LATMXABS
CHARACTER (LEN=1)  :: TUNIT
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('EGGANGLES:VAL_LAT_S',0,ZHOOK_HANDLE)
IF (PRESENT(NUM_ERR))THEN
  TNE = NUM_ERR
ELSE
  TNE = -1_JPIM
ENDIF
IF (PRESENT(PI)) THEN
  TPI = PI
ELSE
  TPI = ASIN(1.0_JPRB)*2.0_JPRB
ENDIF
IF (PRESENT(UNIT)) THEN
  IF ((UNIT=='R').OR.(UNIT=='D')) THEN
    TUNIT = UNIT
  ELSE
    TUNIT = "D"
  ENDIF
ELSE
  TUNIT = "D"
ENDIF

IF (TUNIT=='R') THEN
  LATMXABS = TPI/2.0_JPRB
ELSE
  LATMXABS = 90.0_JPRB
ENDIF

IF (ABS(LAT) > LATMXABS) THEN
  ETAT = TNE
ELSE
  ETAT = 1_JPIM
ENDIF
IF (LHOOK) CALL DR_HOOK('EGGANGLES:VAL_LAT_S',1,ZHOOK_HANDLE)
END FUNCTION VAL_LAT_S
! -------------------------------------------------------------------------------
INTEGER(KIND=JPIM) FUNCTION VAL_LAT_V(P_LAT,NUM_ERR,PI,UNIT) RESULT(ETAT)
REAL(KIND=JPRB), DIMENSION(:), INTENT(IN)                 :: P_LAT
CHARACTER (LEN=1), INTENT(IN), OPTIONAL                   :: UNIT
REAL(KIND=JPRB), INTENT(IN), OPTIONAL                     :: PI
INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL                  :: NUM_ERR

INTEGER(KIND=JPIM) :: TNE
REAL(KIND=JPRB)    :: TPI, Z_LATMXABS
CHARACTER (LEN=1)  :: TUNIT
REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('EGGANGLES:VAL_LAT_V',0,ZHOOK_HANDLE)
IF (PRESENT(NUM_ERR))THEN
  TNE = NUM_ERR
ELSE
  TNE = -1_JPIM
ENDIF
IF (PRESENT(PI)) THEN
  TPI = PI
ELSE
  TPI = ASIN(1.0_JPRB)*2.0_JPRB
ENDIF
IF (PRESENT(UNIT)) THEN
  IF ((UNIT=='R').OR.(UNIT=='D')) THEN
    TUNIT = UNIT
  ELSE
    TUNIT = "D"
  ENDIF
ELSE
  TUNIT = "D"
ENDIF

IF (TUNIT=='R') THEN
  Z_LATMXABS = TPI/2.0_JPRB
ELSE
  Z_LATMXABS = 90.0_JPRB
ENDIF

IF (ANY(ABS(P_LAT(:)) > Z_LATMXABS)) THEN
  ETAT = TNE
ELSE
  ETAT = 1_JPIM
ENDIF
IF (LHOOK) CALL DR_HOOK('EGGANGLES:VAL_LAT_V',1,ZHOOK_HANDLE)
END FUNCTION VAL_LAT_V
! -------------------------------------------------------------------------------
INTEGER(KIND=JPIM) FUNCTION VAL_LON_S(LON,NUM_ERR,PI,DOM,UNIT) RESULT(ETAT)
REAL(KIND=JPRB), INTENT(IN)                                :: LON
CHARACTER (LEN=2), INTENT(IN), OPTIONAL                    :: DOM
CHARACTER (LEN=1), INTENT(IN), OPTIONAL                    :: UNIT
REAL(KIND=JPRB), INTENT(IN), OPTIONAL                      :: PI
INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL                   :: NUM_ERR

INTEGER(KIND=JPIM) :: TNE
REAL(KIND=JPRB)    :: TPI, CVT, S, LONMIN, LONMAX
CHARACTER (LEN=2)  :: TDOM
CHARACTER (LEN=1)  :: TUNIT
REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('EGGANGLES:VAL_LON_S',0,ZHOOK_HANDLE)
IF (PRESENT(NUM_ERR))THEN
  TNE = NUM_ERR
ELSE
  TNE = -1_JPIM
ENDIF
IF (PRESENT(PI)) THEN
  TPI = PI
ELSE
  TPI = ASIN(1.0_JPRB)*2.0_JPRB
ENDIF
IF (PRESENT(DOM)) THEN
  IF ((DOM=='0+').OR.(DOM=='-+')) THEN
    TDOM = DOM
  ELSE
    TDOM = "-+"
  ENDIF
ELSE
  TDOM = "-+"
ENDIF
IF (PRESENT(UNIT)) THEN
  IF ((UNIT=='R').OR.(UNIT=='D')) THEN
    TUNIT = UNIT
  ELSE
    TUNIT = "D"
  ENDIF
ELSE
  TUNIT = "D"
ENDIF

IF (TUNIT=='R') THEN
  CVT = TPI
ELSE
  CVT = 180.0_JPRB
ENDIF
IF (TDOM=='-+') THEN
  S = -1.0_JPRB
ELSE
  S = 0.0_JPRB
ENDIF
LONMIN = S*CVT
LONMAX =(2.0_JPRB +S)*CVT

IF ((LON < LONMIN).OR.(LON >= LONMAX)) THEN
  ETAT = TNE
ELSE
  ETAT = 1_JPIM
ENDIF
IF (LHOOK) CALL DR_HOOK('EGGANGLES:VAL_LON_S',1,ZHOOK_HANDLE)
END FUNCTION VAL_LON_S
! -------------------------------------------------------------------------------
INTEGER(KIND=JPIM) FUNCTION VAL_LON_V(LON,NUM_ERR,PI,DOM,UNIT) RESULT(ETAT)
REAL(KIND=JPRB), DIMENSION(:), INTENT(IN)                       :: LON
CHARACTER (LEN=2), INTENT(IN), OPTIONAL                         :: DOM
CHARACTER (LEN=1), INTENT(IN), OPTIONAL                         :: UNIT
REAL(KIND=JPRB), INTENT(IN), OPTIONAL                           :: PI
INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL                        :: NUM_ERR

INTEGER(KIND=JPIM) :: TNE
REAL(KIND=JPRB)    :: TPI, Z_CVT, Z_S, Z_LONMIN, Z_LONMAX
CHARACTER (LEN=2)  :: TDOM
CHARACTER (LEN=1)  :: TUNIT
REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('EGGANGLES:VAL_LON_V',0,ZHOOK_HANDLE)
IF (PRESENT(NUM_ERR))THEN
  TNE = NUM_ERR
ELSE
  TNE = -1_JPIM
ENDIF
IF (PRESENT(PI)) THEN
  TPI = PI
ELSE
  TPI = ASIN(1.0_JPRB)*2.0_JPRB
ENDIF
IF (PRESENT(DOM)) THEN
  IF ((DOM=='0+').OR.(DOM=='-+')) THEN
    TDOM = DOM
  ELSE
    TDOM = "-+"
  ENDIF
ELSE
  TDOM = "-+"
ENDIF
IF (PRESENT(UNIT)) THEN
  IF ((UNIT=='R').OR.(UNIT=='D')) THEN
    TUNIT = UNIT
  ELSE
    TUNIT = "D"
  ENDIF
ELSE
  TUNIT = "D"
ENDIF

IF (TUNIT=='R') THEN
  Z_CVT = TPI
ELSE
  Z_CVT = 180.0_JPRB
ENDIF
IF (TDOM=='-+') THEN
  Z_S = -1.0_JPRB
ELSE
  Z_S = 0.0_JPRB
ENDIF
Z_LONMIN = Z_S*Z_CVT
Z_LONMAX =(2.0_JPRB +Z_S)*Z_CVT

IF ((ANY(LON(:) < Z_LONMIN)).OR.(ANY(LON(:) >= Z_LONMAX))) THEN
  ETAT = TNE
ELSE
  ETAT = 1_JPIM
ENDIF
IF (LHOOK) CALL DR_HOOK('EGGANGLES:VAL_LON_V',1,ZHOOK_HANDLE)
END FUNCTION VAL_LON_V
! -------------------------------------------------------------------------------
INTEGER(KIND=JPIM) FUNCTION VAL_COORD_S(PT_COORD,NUM_ERR,PI,DOM,UNIT) RESULT(ETAT)
TYPE (LOLA), INTENT(IN)                               :: PT_COORD
CHARACTER (LEN=2), INTENT(IN), OPTIONAL               :: DOM
CHARACTER (LEN=1), INTENT(IN), OPTIONAL               :: UNIT
REAL(KIND=JPRB), INTENT(IN), OPTIONAL                 :: PI
INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL              :: NUM_ERR

INTEGER(KIND=JPIM) :: TNE
REAL(KIND=JPRB)    :: TPI
CHARACTER (LEN=2)  :: TDOM
CHARACTER (LEN=1)  :: TUNIT
REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('EGGANGLES:VAL_COORD_S',0,ZHOOK_HANDLE)
IF (PRESENT(NUM_ERR))THEN
  TNE = NUM_ERR
ELSE
  TNE = -1_JPIM
ENDIF
IF (PRESENT(PI)) THEN
  TPI = PI
ELSE
  TPI = ASIN(1.0_JPRB)*2.0_JPRB
ENDIF
IF (PRESENT(DOM)) THEN
  IF ((DOM=='0+').OR.(DOM=='-+')) THEN
    TDOM = DOM
  ELSE
    TDOM = "-+"
  ENDIF
ELSE
  TDOM = "-+"
ENDIF
IF (PRESENT(UNIT)) THEN
  IF ((UNIT=='R').OR.(UNIT=='D')) THEN
    TUNIT = UNIT
  ELSE
    TUNIT = "D"
  ENDIF
ELSE
  TUNIT = "D"
ENDIF

IF ((VAL_LON(PT_COORD%LON,TNE,TPI,TDOM,TUNIT) == 1_JPIM).AND.(VAL_LAT(PT_COORD%LAT,TNE,TPI,TUNIT) == 1_JPIM)) THEN
  ETAT = 1_JPIM
ELSE
  ETAT = TNE
ENDIF
IF (LHOOK) CALL DR_HOOK('EGGANGLES:VAL_COORD_S',1,ZHOOK_HANDLE)
END FUNCTION VAL_COORD_S
! -------------------------------------------------------------------------------
INTEGER(KIND=JPIM) FUNCTION VAL_COORD_V(YD_PT_COORD,K_NUM_ERR,PI,CD_DOM,CD_UNIT) RESULT(ETAT)
TYPE (LOLA), DIMENSION(:), INTENT(IN)                   :: YD_PT_COORD
CHARACTER (LEN=2), INTENT(IN), OPTIONAL                 :: CD_DOM
CHARACTER (LEN=1), INTENT(IN), OPTIONAL                 :: CD_UNIT
REAL(KIND=JPRB), INTENT(IN), OPTIONAL                   :: PI
INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL                :: K_NUM_ERR

INTEGER(KIND=JPIM) :: I_TNE
CHARACTER (LEN=2)  :: CL_TDOM
REAL(KIND=JPRB)    :: Z_TPI
CHARACTER (LEN=1)  :: CL_TUNIT
REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('EGGANGLES:VAL_COORD_V',0,ZHOOK_HANDLE)
IF (PRESENT(K_NUM_ERR))THEN
  I_TNE = K_NUM_ERR
ELSE
  I_TNE = -1_JPIM
ENDIF
IF (PRESENT(PI)) THEN
  Z_TPI = PI
ELSE
  Z_TPI = ASIN(1.0_JPRB)*2.0_JPRB
ENDIF
IF (PRESENT(CD_DOM)) THEN
  IF ((CD_DOM=='0+').OR.(CD_DOM=='-+')) THEN
    CL_TDOM = CD_DOM
  ELSE
    CL_TDOM = "-+"
  ENDIF
ELSE
  CL_TDOM = "-+"
ENDIF
IF (PRESENT(CD_UNIT)) THEN
  IF ((CD_UNIT=='R').OR.(CD_UNIT=='D')) THEN
    CL_TUNIT = CD_UNIT
  ELSE
    CL_TUNIT = "D"
  ENDIF
ELSE
  CL_TUNIT = "D"
ENDIF

IF ((VAL_LON(YD_PT_COORD(:)%LON,I_TNE,Z_TPI,CL_TDOM,CL_TUNIT) == 1_JPIM).AND. &
 & (VAL_LAT(YD_PT_COORD(:)%LAT,I_TNE,Z_TPI,CL_TUNIT) == 1_JPIM)) THEN
  ETAT = 1_JPIM
ELSE
  ETAT = I_TNE
ENDIF
IF (LHOOK) CALL DR_HOOK('EGGANGLES:VAL_COORD_V',1,ZHOOK_HANDLE)
END FUNCTION VAL_COORD_V
! -------------------------------------------------------------------------------
TYPE(LOLA) FUNCTION LOLAR_S (COORD_DEG) RESULT (COORD_RAD)
! DEG => RAD for lola type
TYPE(LOLA), INTENT(IN)                      :: COORD_DEG

REAL(KIND=JPRB) :: TPI,DTR
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('EGGANGLES:LOLAR_S',0,ZHOOK_HANDLE)
TPI = ASIN(1.0_JPRB)*2.0_JPRB
DTR = TPI/180.0_JPRB
COORD_RAD%LON = COORD_DEG%LON*DTR
COORD_RAD%LAT = COORD_DEG%LAT*DTR
IF (LHOOK) CALL DR_HOOK('EGGANGLES:LOLAR_S',1,ZHOOK_HANDLE)
END FUNCTION LOLAR_S

FUNCTION LOLAR_V (COORD_DEG) RESULT (COORD_RAD)
! DEG => RAD for lola type
TYPE(LOLA), DIMENSION(:), INTENT(IN)                :: COORD_DEG
TYPE(LOLA), DIMENSION(SIZE(COORD_DEG)) :: COORD_RAD

REAL(KIND=JPRB) :: TPI,DTR
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('EGGANGLES:LOLAR_V',0,ZHOOK_HANDLE)
TPI = ASIN(1.0_JPRB)*2.0_JPRB
DTR = TPI/180.0_JPRB
COORD_RAD(:)%LON = COORD_DEG(:)%LON*DTR
COORD_RAD(:)%LAT = COORD_DEG(:)%LAT*DTR
IF (LHOOK) CALL DR_HOOK('EGGANGLES:LOLAR_V',1,ZHOOK_HANDLE)
END FUNCTION LOLAR_V
! -------------------------------------------------------------------------------
TYPE(LOLA) FUNCTION LOLAD_S (COORD_RAD) RESULT (COORD_DEG)
! RAD => DEG for lola type
TYPE(LOLA), INTENT(IN)                      :: COORD_RAD

REAL(KIND=JPRB) :: TPI,RTD
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('EGGANGLES:LOLAD_S',0,ZHOOK_HANDLE)
TPI = ASIN(1.0_JPRB)*2.0_JPRB
RTD = 180.0_JPRB/TPI
COORD_DEG%LON = COORD_RAD%LON*RTD
COORD_DEG%LAT = COORD_RAD%LAT*RTD
IF (LHOOK) CALL DR_HOOK('EGGANGLES:LOLAD_S',1,ZHOOK_HANDLE)
END FUNCTION LOLAD_S

FUNCTION LOLAD_V (COORD_RAD) RESULT (COORD_DEG)
! RAD => DEG for lola type
TYPE(LOLA), DIMENSION(:), INTENT(IN)                :: COORD_RAD
TYPE(LOLA), DIMENSION(SIZE(COORD_RAD)) :: COORD_DEG

REAL(KIND=JPRB) :: TPI,RTD
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('EGGANGLES:LOLAD_V',0,ZHOOK_HANDLE)
TPI = ASIN(1.0_JPRB)*2.0_JPRB
RTD = 180.0_JPRB/TPI
COORD_DEG(:)%LON = COORD_RAD(:)%LON*RTD
COORD_DEG(:)%LAT = COORD_RAD(:)%LAT*RTD
IF (LHOOK) CALL DR_HOOK('EGGANGLES:LOLAD_V',1,ZHOOK_HANDLE)
END FUNCTION LOLAD_V
! -------------------------------------------------------------------------------
! Function to compute Cosine,Sine to Angle
! -------------------------------------------------------------------------------
REAL(KIND=JPRB) FUNCTION COSIN_TO_ANGLE_S(COSINUS,SINUS) RESULT (ANGLE)
! (Cosinus,Sinus) => Angle
REAL(KIND=JPRB), INTENT(IN)                  :: COSINUS,SINUS

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('EGGANGLES:COSIN_TO_ANGLE_S',0,ZHOOK_HANDLE)
ANGLE = P_ACOS(COSINUS)*SIGN(1.0_JPRB,SINUS)
IF (LHOOK) CALL DR_HOOK('EGGANGLES:COSIN_TO_ANGLE_S',1,ZHOOK_HANDLE)
END FUNCTION COSIN_TO_ANGLE_S

FUNCTION COSIN_TO_ANGLE_V(COSINUS,SINUS) RESULT (ANGLE)
! (Cosinus,Sinus) => Angle
REAL(KIND=JPRB), DIMENSION(:), INTENT(IN)              :: COSINUS,SINUS
REAL(KIND=JPRB), DIMENSION(SIZE(COSINUS)) :: ANGLE

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('EGGANGLES:COSIN_TO_ANGLE_V',0,ZHOOK_HANDLE)
ANGLE(:) = P_ACOS(COSINUS(:))*SIGN(1.0_JPRB,SINUS(:))
IF (LHOOK) CALL DR_HOOK('EGGANGLES:COSIN_TO_ANGLE_V',1,ZHOOK_HANDLE)
END FUNCTION COSIN_TO_ANGLE_V
! -------------------------------------------------------------------------------
! -------------------------------------------------------------------------------
! Function to compute Acos without error
! -------------------------------------------------------------------------------
REAL(KIND=JPRB) FUNCTION P_ACOS_S(COSINUS) RESULT (ANGLE)
! Protected ACOS
REAL(KIND=JPRB), INTENT(IN)                  :: COSINUS

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('EGGANGLES:P_ACOS_S',0,ZHOOK_HANDLE)
ANGLE = ACOS(MINIMAX(COSINUS))
IF (LHOOK) CALL DR_HOOK('EGGANGLES:P_ACOS_S',1,ZHOOK_HANDLE)
END FUNCTION P_ACOS_S

FUNCTION P_ACOS_V(COSINUS) RESULT (ANGLE)
! Protected ACOS
REAL(KIND=JPRB), DIMENSION(:), INTENT(IN)              :: COSINUS
REAL(KIND=JPRB), DIMENSION(SIZE(COSINUS)) :: ANGLE

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('EGGANGLES:P_ACOS_V',0,ZHOOK_HANDLE)
ANGLE(:) = ACOS(MINIMAX(COSINUS(:)))
IF (LHOOK) CALL DR_HOOK('EGGANGLES:P_ACOS_V',1,ZHOOK_HANDLE)
END FUNCTION P_ACOS_V
! -------------------------------------------------------------------------------
! -------------------------------------------------------------------------------
! Function to compute Asin without error
! -------------------------------------------------------------------------------
REAL(KIND=JPRB) FUNCTION P_ASIN_S(SINUS) RESULT (ANGLE)
! Protected ASIN
REAL(KIND=JPRB), INTENT(IN)                  :: SINUS

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('EGGANGLES:P_ASIN_S',0,ZHOOK_HANDLE)
ANGLE = ASIN(MINIMAX(SINUS))
IF (LHOOK) CALL DR_HOOK('EGGANGLES:P_ASIN_S',1,ZHOOK_HANDLE)
END FUNCTION P_ASIN_S

FUNCTION P_ASIN_V(SINUS) RESULT (ANGLE)
! Protected ASIN
REAL(KIND=JPRB), DIMENSION(:), INTENT(IN)            :: SINUS
REAL(KIND=JPRB), DIMENSION(SIZE(SINUS)) :: ANGLE

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('EGGANGLES:P_ASIN_V',0,ZHOOK_HANDLE)
ANGLE(:) = ASIN(MINIMAX(SINUS(:)))
IF (LHOOK) CALL DR_HOOK('EGGANGLES:P_ASIN_V',1,ZHOOK_HANDLE)
END FUNCTION P_ASIN_V
! -------------------------------------------------------------------------------
! -------------------------------------------------------------------------------
! Function MinMax
! -------------------------------------------------------------------------------
REAL(KIND=JPRB) FUNCTION MINIMAX_S(VAL,LIM) RESULT (VALO)
! Return Value in [-LIM,LIM]
REAL(KIND=JPRB), INTENT(IN)                      :: VAL
REAL(KIND=JPRB), INTENT(IN), OPTIONAL            :: LIM

REAL(KIND=JPRB) :: TLIM
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('EGGANGLES:MINIMAX_S',0,ZHOOK_HANDLE)
IF (PRESENT(LIM)) THEN
  TLIM = LIM
ELSE
  TLIM = 1.0_JPRB
ENDIF
VALO = MIN(TLIM,MAX(-1.0_JPRB*TLIM,VAL))
IF (LHOOK) CALL DR_HOOK('EGGANGLES:MINIMAX_S',1,ZHOOK_HANDLE)
END FUNCTION MINIMAX_S

FUNCTION MINIMAX_V(VAL,LIM) RESULT (VALO)
! Return Value in [-LIM,LIM]
REAL(KIND=JPRB), DIMENSION(:), INTENT(IN)          :: VAL
REAL(KIND=JPRB), INTENT(IN), OPTIONAL              :: LIM
REAL(KIND=JPRB), DIMENSION(SIZE(VAL)) :: VALO

REAL(KIND=JPRB) :: TLIM
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('EGGANGLES:MINIMAX_V',0,ZHOOK_HANDLE)
IF (PRESENT(LIM)) THEN
  TLIM = LIM
ELSE
  TLIM = 1.0_JPRB
ENDIF
VALO(:) = MIN(TLIM,MAX(-1.0_JPRB*TLIM,VAL(:)))
IF (LHOOK) CALL DR_HOOK('EGGANGLES:MINIMAX_V',1,ZHOOK_HANDLE)
END FUNCTION MINIMAX_V
! -------------------------------------------------------------------------------
! -------------------------------------------------------------------------------
! Functions Longitude size/distance
! -------------------------------------------------------------------------------
REAL(KIND=JPRB) FUNCTION DIST_2REF_L(COORD_LON,REF_LON,PI) RESULT(DIST)
! COORD_LON, REF_LON in -+Radians
! DIST in -+Radians

! Calcule la distance orientee DIST (abscisse dans l'intervale [-pi,pi[ et d'origine le meridien de
! reference REF_LON) de COORD_LON a REF_LON (coordonnees en radians dans l'intervale
! [-pi,pi[ et avec comme origine le meridien de GreenWiTch). Les valeurs negatives sont vers l'Ouest
! des origines.

! Computes oriented distance DIST (as an absciss in [-pi,pi[ with origin at REF_LON meridian)
! from COORD_LON to REF_LON (coordinates in rad [-pi,pi[ with origin at GreenWiTch meridian)
! Negatives values are on West of origins.

REAL(KIND=JPRB), INTENT(IN)                   :: COORD_LON, REF_LON
REAL(KIND=JPRB), INTENT(IN), OPTIONAL         :: PI 

REAL(KIND=JPRB)                            :: Z  
REAL(KIND=JPRB)                            :: TPI
REAL(KIND=JPHOOK)                            :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('EGGANGLES:DIST2REF_L',0,ZHOOK_HANDLE)    

IF (PRESENT(PI)) THEN
  TPI = PI
ELSE
  TPI = ASIN(1.0_JPRB)*2.0_JPRB
ENDIF
Z = (COORD_LON-REF_LON)
Z = Z-SIGN(TPI,Z)*(1.0_JPRB+SIGN(1.0_JPRB,ABS(Z)-TPI))
DIST = -Z*SIGN(1.0_JPRB,Z-TPI) ! because [-pi,pi[ : if pi then -pi 

IF (LHOOK) CALL DR_HOOK('EGGANGLES:DIST2REF_L',1,ZHOOK_HANDLE)
END FUNCTION DIST_2REF_L
! -------------------------------------------------------------------------------
REAL(KIND=JPRB) FUNCTION DIST_2REF_S(PT_COORD,REF_COORD,PI) RESULT(DIST)
! PT_COORD, REF_COORD in -+Radians
! DIST in -+Radians

! Calcule la distance orientee DIST (abscisse dans l'intervale [-pi,pi[ et d'origine le meridien de
! reference REF_COORD%LON) de PT_COORD%LON a REF_COORD%LON (coordonnees en radians dans l'intervale
! [-pi,pi[ et avec comme origine le meridien de GreenWiTch). Les valeurs negatives sont vers l'Ouest
! des origines.

! Computes oriented distance DIST (as an absciss in [-pi,pi[ with origin at REF_COORD%LON meridian)
! from PT_COORD%LON to REF_COORD%LON (coordinates in rad [-pi,pi[ with origin at GreenWiTch meridian)
! Negatives values are on West of origins.

TYPE (LOLA), INTENT(IN)                       :: PT_COORD, REF_COORD
REAL(KIND=JPRB), INTENT(IN), OPTIONAL         :: PI 

REAL(KIND=JPRB)                            :: TPI
REAL(KIND=JPHOOK)                            :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('EGGANGLES:DIST2REF_S',0,ZHOOK_HANDLE)    

IF (PRESENT(PI)) THEN
  TPI = PI
ELSE
  TPI = ASIN(1.0_JPRB)*2.0_JPRB
ENDIF
DIST = DIST_2REF(PT_COORD%LON,REF_COORD%LON,TPI)

IF (LHOOK) CALL DR_HOOK('EGGANGLES:DIST2REF_S',1,ZHOOK_HANDLE)
END FUNCTION DIST_2REF_S
! -------------------------------------------------------------------------------
FUNCTION DIST_2REF_V(PT_COORD,REF_COORD,PI) RESULT(DIST)
! PT_COORD, REF_COORD in -+Radians
! DIST in -+Radians

! Calcule la distance orientee DIST (abscisse dans l'intervale [-pi,pi[ et d'origine le meridien de
! reference REF_COORD%LON) de PT_COORD%LON a REF_COORD%LON (coordonnees en radians dans l'intervale
! [-pi,pi[ et avec comme origine le meridien de GreenWiTch). Les valeurs negatives sont vers l'Ouest
! des origines.

! Computes oriented distance DIST (as an absciss in [-pi,pi[ with origin at REF_COORD%LON meridian)
! from PT_COORD%LON to REF_COORD%LON (coordinates in rad [-pi,pi[ with origin at GreenWiTch meridian)
! Negatives values are on West of origins.

TYPE (LOLA), DIMENSION(:), INTENT(IN)                   :: PT_COORD
TYPE (LOLA), INTENT(IN)                                 :: REF_COORD
REAL(KIND=JPRB), INTENT(IN), OPTIONAL                   :: PI 
REAL(KIND=JPRB), DIMENSION(SIZE(PT_COORD)) :: DIST

REAL(KIND=JPRB), DIMENSION(SIZE(PT_COORD)) :: Z
REAL(KIND=JPRB)                            :: TPI
REAL(KIND=JPHOOK)                            :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('EGGANGLES:DIST2REF_V',0,ZHOOK_HANDLE)    

IF (PRESENT(PI)) THEN
  TPI = PI
ELSE
  TPI = ASIN(1.0_JPRB)*2.0_JPRB
ENDIF
Z(:) = PT_COORD(:)%LON-REF_COORD%LON
Z(:) = Z(:)-SIGN(TPI,Z(:))*(1.0_JPRB+SIGN(1.0_JPRB,ABS(Z(:))-TPI))
DIST(:) = -Z(:)*SIGN(1.0_JPRB,Z(:)-TPI) ! because [-pi,pi[ : if pi then -pi 

IF (LHOOK) CALL DR_HOOK('EGGANGLES:DIST2REF_V',1,ZHOOK_HANDLE)
END FUNCTION DIST_2REF_V
! -------------------------------------------------------------------------------
REAL(KIND=JPRB) FUNCTION SIZE_W2E_L(WEST_LON,EAST_LON,PI) RESULT(TAILLE)
! WEST_LON, EAST_LON in -+Radians
! SIZE in 0+Radians 

! Calcule la distance ou taille (norme entre ]0,2pi]) entre WEST_LON et EAST_LON
! dans le sens des aiguilles d'une montre en regardant du Pole Sud vers le Pole Nord ( ceci
! pour ne pas calculer le complementaire a 2pi ).

! Computes distance or length (norm in ]0,2pi]) between WEST_LON and EAST_LON in
! clockwise seeing from South Pole to North Pole.

REAL(KIND=JPRB), INTENT(IN)                   :: WEST_LON, EAST_LON
REAL(KIND=JPRB), INTENT(IN), OPTIONAL         :: PI 

REAL(KIND=JPRB)                            :: Z  
REAL(KIND=JPRB)                            :: TPI
REAL(KIND=JPHOOK)                            :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('EGGANGLES:SIZE_W2E_L',0,ZHOOK_HANDLE)    

IF (PRESENT(PI)) THEN
  TPI = PI
ELSE
  TPI = ASIN(1.0_JPRB)*2.0_JPRB
ENDIF
Z =  DIST_2REF(WEST_LON,EAST_LON,TPI)
TAILLE = TPI*(1.0_JPRB+SIGN(1.0_JPRB,Z))-Z

IF (LHOOK) CALL DR_HOOK('EGGANGLES:SIZE_W2E_L',1,ZHOOK_HANDLE)
END FUNCTION SIZE_W2E_L
! -------------------------------------------------------------------------------
REAL(KIND=JPRB) FUNCTION SIZE_W2E_S(WEST_COORD,EAST_COORD,PI) RESULT(TAILLE)
! WEST_COORD, EAST_COORD in -+Radians
! SIZE in 0+Radians 

! Calcule la distance ou taille (norme entre ]0,2pi]) entre WEST_COORD%LON et EAST_COORD%LON
! dans le sens des aiguilles d'une montre en regardant du Pole Sud vers le Pole Nord ( ceci
! pour ne pas calculer le complementaire a 2pi ).

! Computes distance or length (norm in ]0,2pi]) between WEST_COORD%LON and EAST_COORD%LON in
! clockwise seeing from South Pole to North Pole.

TYPE (LOLA), INTENT(IN)                       :: WEST_COORD, EAST_COORD
REAL(KIND=JPRB), INTENT(IN), OPTIONAL         :: PI 

REAL(KIND=JPRB)                            :: TPI
REAL(KIND=JPHOOK)                            :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('EGGANGLES:SIZE_W2E_S',0,ZHOOK_HANDLE)    

IF (PRESENT(PI)) THEN
  TPI = PI
ELSE
  TPI = ASIN(1.0_JPRB)*2.0_JPRB
ENDIF
TAILLE = SIZE_W2E(WEST_COORD%LON,EAST_COORD%LON,TPI)

IF (LHOOK) CALL DR_HOOK('EGGANGLES:SIZE_W2E_S',1,ZHOOK_HANDLE)
END FUNCTION SIZE_W2E_S
! -------------------------------------------------------------------------------
END MODULE EGGANGLES
