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

MODULE EGGPACK

! Version 2009.0317 by JD GRIL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! DOC !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! This package contains many usefull functions about Projection in "Tangent 
!Case" for Lambert Conic Comformal, Polar Stereographic and Mercator. A Domain
!Maker subroutine for grid points domain defined by its center is present.
! * The functions are classed in three types :
!  - Independants functions
!  - Specifics functions -| theses functions work for single point data structure
!  - Generics functions  -| or for array of single point data structure (rank of
!                         | array is one)
! * One subroutine to make domain

! -1- Decription of Functions :

!  -1.1- Independants functions : 

!    Theses functions have not direct relations with projection but are usefull 
!    for this package.

!->logical function RETURN_PRINT(CODE_ERR,NUM_TEST,AUTO_STOP) 

!     Return a logical to know if in case of error the the program will stop
!     or return an error code to the caller depending AUTO_STOP logical flag.
!     In every cases it prints informations (ERROR,WARNING,INFO,OK) contained in
!     CODE_ERR structure of type ERROR and NUM_TEST. It's used in subroutine
!     MAKDO. See how inside. 

!->function TYPE_PROJ(REF_COORD)

!     Return 1 character (M,L,S) depending the coordinates of reference point
!     REF_COORD structure of type LOLA (in Degrees).

!->real function POLE_IS (REF_COORD) 

!     Return pole location (South is -1.0, North is 1.0) for Polar Stereographic or
!     Lambert projections (return 0.0 in Mercator : no sense) depending the coor-
!     dinates of reference point REF_COORD structure of type LOLA (in Degrees).

!  -1.2- Specifics functions : 

!    These functions are only used in Polar Stereographic or Lambert projection
!    type. Depending the type of arguments (arrays or not) the generic caller-
!    name defined by interface selects the good one (suffixed by _V for arrays
!     or _S for scalars).
!    All of these functions are prefixed by STPL_; they use P_PJ structure of
!    type PARAM_PROJ (the most important of this package) that defines the 
!    parameters of the projection depending the reference point, see below.

!->type (RTETA) function STPL_LATLON_TO_RTETA(PT_COORD,P_PJ,PI) result (PT_RTETA)

!->type (LOLA)  function STLP_RTETA_TO_LATLON(PT_RTETA,P_PJ,PI) result (PT_COORD)

!->type (XY)    function STLP_RTETA_TO_XY(PT_RTETA,P_PJ)        result (PT_XY)

!->type (RTETA) function STLP_XY_TO_RTETA(PT_XY,P_PJ,PI)        result (PT_RTETA)

!     These functions convert coordinates from LAT-LON (in Radians) to R-THETA
!     (in meters, radians) and from R-THETA to XY (in meters) and reverse. PT_???
!     are structures or arrays of structures. PT_XY coordinates are referenced by
!     STD (standard) origin that is projection of pole, the X axis is West to
!     East, the Y axis is South to North. PT_RTETA coordinates are referenced by
!     STD (standard) origin that is projection of pole, R is the distance to it,
!     TETA is angle with reference point longitude (negative in clockwise).
!            Remark -> PT_XY and PT_RTETA are on projection plane !

!  -1.3- Generics functions :

!    These functions are standards for all types of projections (M,L,S). the most
!    important is :

!->type (PARAM_PROJ) function REF_DATAS (REF_COORD,RA,TOZERO_COORD,LRT) 
!  result (P_P) 

!      It defines a P_P structure of type PARAM_PROJ that creates all the
!    parameters needs by all specifics and generics functions : this structure
!    defines the type of projection depending the reference point coordinates
!    (in Degrees, inside you have the reference point coordinates in Radians).
!      For all others functions, the generic caller-name selects the good one
!    (suffixed by _V for arrays or _S for scalars) depending the type of 
!    arguments (arrays or not). they are three classes of functions :

!    a) changing origin of XY type of coordinates :

!     In Mercator projection the YX origin is the reference point, because 
!     projection of pole is non sense. 

!->type (XY) function XY_NEW_TO_STD_ORIGIN(NEW_ORIGIN_COORD,PT_XY_IN_NEW_ORIGIN,
!  P_PJ,PI) result (PT_XY_IN_STD_ORIGIN) 

!     This function changes coordinates PT_XY_IN_NEW_ORIGIN (may be an array)
!     relative to NEW_ORIGIN_COORD (in Degrees) to STD origin (depending of 
!     projection).

!->type (XY) function XY_STD_TO_NEW_ORIGIN(NEW_ORIGIN_COORD,PT_XY_IN_STD_ORIGIN,
!  P_PJ,PI) result (PT_XY_IN_NEW_ORIGIN)

!     This function changes coordinates PT_XY_IN_STD_ORIGIN (may be an array)
!     relative to STD origin (depending of projection) to NEW_ORIGIN_COORD 
!     (in Degrees).

!            Remark -> To pass PT_XY_IN_OLD_ORIGIN to PT_XY_IN_NEW_ORIGIN, you
!                      can do (or create your own function) :
!       PT_XY_IN_NEW_ORIGIN =                                                   &
!           XY_STD_TO_NEW_ORIGIN(NEW_ORIGIN_COORD,                              &
!           XY_NEW_TO_STD_ORIGIN(OLD_ORIGIN_COORD,PT_XY_IN_OLD_ORIGIN,P_PJ,PI), &
!                                                 P_PJ,PI)  

!    b) converting LAT-LON coordinates in XY and reverse :

!->type (XY)   function LATLON_TO_XY(PT_COORD,P_PJ,PI) result (PT_XY)

!->type (LOLA) function XY_TO_LATLON(PT_XY,P_PJ,PI)    result (PT_COORD)

!     These functions convert coordinates from LAT-LON (in Radians) to XY (in
!     meters) and reverse. PT_??? are structures or arrays of structures. PT_XY
!     coordinates are referenced by STD (standard) origin that is projection of
!     pole, for Lambert or Polar Stereographic projection or reference point for
!     Mercator; the X axis is West to East, the Y axis is South to North.
!            Remark -> PT_XY is on projection plane !

!    c) Map Factor,  projection of Geographical North : 

!->real function MAP_FACTOR(PT_COORD,P_PJ,PI,RA)

!     Compute Map Factor for PT_COORD (in Radians) (array of) structure of type
!     LOLA depending type of projection.

!->type (PGN) function GN(PT_COORD,P_PJ)  

!     Compute Projection of Geographical North unit vector on XY for PT_COORD (in
!     Radians) (array of) structure of type LOLA depending type of projection.

! -2- Decription of Subroutine :       

!  Subroutine MAKDO (make domain) use CENTER point of the Domain as the Origin 
!  of XY points.

!->subroutine MAKDO(REF_COORD,CENTER_COORD,PDEL,NB_PTS,GRID_COORD,GRID_MF,
!  GRID_PGN,GRID_INFO,ERR_CODE,AUTO_STOP,PI,RA,LMRT)

!   This subroutine creates a grid domain usefull datas depending projection.
!  It makes validations tests, using functions of this package computes
!  outputs.
!   The X and Y points are computed from the center in a 2 dimensional array
!  GRID_XY_C structure of type XY
!   To pass 2 dimensional arrays to one dimensional in functions we must
!  use the pack F90 function and unpack for the reverse because all the
!  functions are defined in interfaces.

!  - Inputs-Outputs :
!    REF_COORD, CENTER_COORD : in Degrees in input, in Radians in output
!                              type LOLA
!                              defining Reference coordinates for the
!                              projection and center coordinates of domain.
!  - Inputs :
!    PDEL                    : resolution in meters of the grid in projection
!                              plane.
!                              type DELTA
!    NB_PTS                  : number of points on X and Y
!                              type NBPTS
!  - Optionals Inputs :
!    PI, RA                  : Pi and Raduis of Earth
!                              type real
!    AUTO_STOP               : specify if ERROR stops the program (default)
!                              or return ERR_CODE to caller
!                              type logical
!  - Outputs :
!    ERR_CODE                : error code returned
!                              type ERROR
!    GRID_INFO               : informations on grid domain (projections,
!                              reference point, center, corners, map factors)
!                              type DOMI
!    GRID_COORD              : grid coordinates LAT [-90.0,90.0] in Radian
!                                               LON [0.0,360.0[ in Radians
!                              array dim = 2 : NB_PTS%ONX,NB_PTS%ONY
!                              type LOLA
!    GRID_MF                 : grid map factor
!                              array dim = 2 : NB_PTS%ONX,NB_PTS%ONY
!                              type real
!    GRID_PGN                : projection of Geographical North unit vector
!                              on X and Y
!                              array dim = 2 : NB_PTS%ONX,NB_PTS%ONY
!                              type PGN

!            Remark -> the EBETA parameter in old EGGX is fully automatic
!                      for Lambert and Polar Stereographic : it depends only
!                      of reference longitude. In mercator, it has non sense
!                      because alls meridians are parallels to reference
!                      longitude.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Author : Jean-Daniel GRIL , CNRM/GMAP/EXT , March 28 2002
! Modified 04-08-19 R. El Khatib & J.-D. Gril : cleanups
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modifications :
!  _______________

!  November 09, 2004 : JD GRIL
!  --------------------------
!      - supress of many "#define" values
!      - add 2 news optionals arguments in REF_DATAS :
!            - TOZERO_COORD coordinates of domain's center that is placed in (0,0)
!              coord. if Rotated-Tilted Mercator Projection is used
!            - LRT flag to .TRUE. if you use Rotated-Tilted Mercator Projection
!            The use of this flag defines a "W" value for P_P%TYPE_PJ (projection
!            type) and the use of EGGMRT module. The functions and MAKDO routine 
!            was modified according this new feature. 
!      - MAKDO has a new optional parameter to define Rotated-Tilted Mercator
!        Projection :
!            - LD_LMRT 

!      To use Rotated-Tilted Mercator projection, remember that reference coordinates
!      YD_REF_COORD%LAT = 0 (mercator) and YD_REF_COORD%LON is the tilting angle
!      of the domain (at the point of coordinate (0,0)),
!      after the rotation of the center of the domain (YD_CENTER_COORD coordinates)
!      to the point of coordinates (0,0). 
!      All use of functions is unchanged (the new description of projection is
!      done by the P_P argument).

!      - add some protections in Lambert (STLP_XY_TO_RTETA functions) to prevent
!        the nosense of XY values in the cutting part of Lambert projected plane.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modified 11-07-2006 : JD Gril : use of LOLAD & LOLAR functions vs /DTR & *DTR
!  Modified 16-10-2006 : JD Gril : replace pack/unpack by reshape intrinsic func.
!  Modified 12-01-2007 : JD Gril : cleaning
!  Modified 04-07-2008 : JD Gril : latlon_to_xy and stpl_latlon_to_rteta edge
!                                : effect for huge domain : pb diff longitudes.
!                                : replace reshape by do loop
!  Modified 13-03-2009 : JD Gril : Optimization et cleanup
!  R. El Khatib 22-Mar-2017 disable verbosity if TKOUT < 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

! ******************* Definition of parameters **********************************

!     Include Kinds
!     -------------

#ifndef DEBUG
! force no print info in makdo by default
#define _DEFP_ .FALSE.
#else
! force print info in makdo by default
#define _DEFP_ .TRUE.
#endif

USE PARKIND1  ,ONLY : JPIM,    JPRB, JPRD
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

! ******************* Loading module ********************************************

USE EGGANGLES  ,ONLY : LOLA, VAL_COORD, LOLAR, LOLAD, ANGLE_DOMAIN, DIST_2REF
USE EGGMRT  ,ONLY : MEROTIL, METILROT 

IMPLICIT NONE

! Physics constant
! ----------------

REAL(KIND=JPRB), PARAMETER, PUBLIC :: R_EARTH = 6371229._JPRB

! ******************* Definition of type ****************************************

TYPE ERROR
  INTEGER(KIND=JPIM) :: NUM
  CHARACTER(LEN=100) :: TXT
END TYPE ERROR

TYPE PGN
  REAL(KIND=JPRB) :: ONX, ONY
END TYPE PGN

TYPE NBPTS
  INTEGER(KIND=JPIM) :: ONX, ONY
END TYPE NBPTS

TYPE DELTA
  REAL(KIND=JPRB) :: ONX, ONY
END TYPE DELTA

TYPE RTETA
  REAL(KIND=JPRD) :: R, TETA
END TYPE RTETA

TYPE XY
  REAL(KIND=JPRB) :: X, Y
END TYPE XY

TYPE PARAM_PROJ
  SEQUENCE
  TYPE (LOLA) :: REF_PT,TZO_PT
  REAL(KIND=JPRD) :: KL
  REAL(KIND=JPRD) :: R_EQUATEUR
  REAL(KIND=JPRB) :: POLE
  CHARACTER(LEN=1) :: TYPE_PJ
END TYPE PARAM_PROJ

TYPE DOMI
  TYPE (NBPTS) :: G_SIZE
  TYPE (LOLA) :: CT_COORD, RF_COORD, SW_COORD, SE_COORD, NE_COORD, NW_COORD
  REAL(KIND=JPRB) :: MF_CT, MF_RF, MF_SW, MF_SE, MF_NE, MF_NW
  TYPE (PARAM_PROJ) :: INFO_PROJ
END TYPE DOMI

! ******************* Definition of Interface ***********************************

INTERFACE INFO_PRINT
  MODULE PROCEDURE INFO_PP_PRINT, INFO_DOMI_PRINT
END INTERFACE

INTERFACE XY_NEW_TO_STD_ORIGIN
  MODULE PROCEDURE XY_NEW_TO_STD_ORIGIN_V, XY_NEW_TO_STD_ORIGIN_S
END INTERFACE

INTERFACE XY_STD_TO_NEW_ORIGIN
  MODULE PROCEDURE XY_STD_TO_NEW_ORIGIN_V, XY_STD_TO_NEW_ORIGIN_S
END INTERFACE

INTERFACE STPL_LATLON_TO_RTETA
  MODULE PROCEDURE STPL_LATLON_TO_RTETA_V, STPL_LATLON_TO_RTETA_S
END INTERFACE

INTERFACE STLP_RTETA_TO_XY
  MODULE PROCEDURE STLP_RTETA_TO_XY_V, STLP_RTETA_TO_XY_S
END INTERFACE

INTERFACE STLP_XY_TO_RTETA
  MODULE PROCEDURE STLP_XY_TO_RTETA_V, STLP_XY_TO_RTETA_S
END INTERFACE

INTERFACE STLP_RTETA_TO_LATLON
  MODULE PROCEDURE STLP_RTETA_TO_LATLON_V, STLP_RTETA_TO_LATLON_S
END INTERFACE

INTERFACE LATLON_TO_XY
  MODULE PROCEDURE LATLON_TO_XY_V, LATLON_TO_XY_S
END INTERFACE

INTERFACE XY_TO_LATLON
  MODULE PROCEDURE XY_TO_LATLON_V, XY_TO_LATLON_S
END INTERFACE

INTERFACE MAP_FACTOR
  MODULE PROCEDURE MAP_FACTOR_V, MAP_FACTOR_S
END INTERFACE

INTERFACE GN
  MODULE PROCEDURE GN_V, GN_S
END INTERFACE

#include "abor1.intfb.h"

CONTAINS

! =================== FUNCTIONS =================================================

! ******************* Independants functions ************************************

! -------------------------------------------------------------------------------
SUBROUTINE INFO_DOMI_PRINT(YD_G_INFO,KOUT,PI)
TYPE (DOMI), INTENT(IN)                  :: YD_G_INFO
INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: KOUT
REAL(KIND=JPRB), INTENT(IN), OPTIONAL    :: PI

REAL(KIND=JPRB)      :: TPI
INTEGER(KIND=JPIM)   :: TKOUT
REAL(KIND=JPHOOK)      :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('EGGPACK:INFO_DOMI_PRINT',0,ZHOOK_HANDLE)
IF (PRESENT(KOUT))THEN
  TKOUT = KOUT
ELSE
  TKOUT = 6_JPIM
ENDIF
IF (PRESENT(PI)) THEN
  TPI = PI
ELSE
  TPI = ASIN(1.0_JPRB)*2.0_JPRB
ENDIF
CALL INFO_PRINT(YD_G_INFO%INFO_PROJ,TKOUT,TPI)
IF (TKOUT >= 0) WRITE(TKOUT,*) "============================================================="
IF (TKOUT >= 0) WRITE(TKOUT,*) "===   Informations about Domain Information Structure    ===="
IF (TKOUT >= 0) WRITE(TKOUT,*) "============================================================="
IF (TKOUT >= 0) WRITE(TKOUT,*)
IF (TKOUT >= 0) WRITE(TKOUT,*) "  -Size of Domain (in points) :"
IF (TKOUT >= 0) WRITE(TKOUT,'(13X,A7,10X,A7)') " On X  "," On Y  "
IF (TKOUT >= 0) WRITE(TKOUT,'(13X,I7,10X,I7)') YD_G_INFO%G_SIZE%ONX,YD_G_INFO%G_SIZE%ONY
IF (TKOUT >= 0) WRITE(TKOUT,*) "  -Most important points informations :"
IF (TKOUT >= 0) WRITE(TKOUT,'(1X,A7,1X,"|",1X,A9,1X,"|",1X,A8,1X,"|",1X,A16)') &
 & " Points ","Longitude","Latitude","   Map Factor   "
IF (TKOUT >= 0) WRITE(TKOUT,'(A47)') "------------------------------------------------"
IF (TKOUT >= 0) WRITE(TKOUT,'(1X,A7,1X,"|",2X,F7.2,2X,"|",2X,F7.2,1X,"|",1X,G16.10)') &
 & "Center ",YD_G_INFO%CT_COORD%LON,YD_G_INFO%CT_COORD%LAT,YD_G_INFO%MF_CT
IF (TKOUT >= 0) WRITE(TKOUT,'(1X,A7,1X,"|",2X,F7.2,2X,"|",2X,F7.2,1X,"|",1X,G16.10)') &
 & "Refer. ",YD_G_INFO%RF_COORD%LON,YD_G_INFO%RF_COORD%LAT,YD_G_INFO%MF_RF
IF (TKOUT >= 0) WRITE(TKOUT,'(1X,A7,1X,"|",2X,F7.2,2X,"|",2X,F7.2,1X,"|",1X,G16.10)') &
 & "S.West ",YD_G_INFO%SW_COORD%LON,YD_G_INFO%SW_COORD%LAT,YD_G_INFO%MF_SW
IF (TKOUT >= 0) WRITE(TKOUT,'(1X,A7,1X,"|",2X,F7.2,2X,"|",2X,F7.2,1X,"|",1X,G16.10)') &
 & "S.East ",YD_G_INFO%SE_COORD%LON,YD_G_INFO%SE_COORD%LAT,YD_G_INFO%MF_SE
IF (TKOUT >= 0) WRITE(TKOUT,'(1X,A7,1X,"|",2X,F7.2,2X,"|",2X,F7.2,1X,"|",1X,G16.10)') &
 & "N.East ",YD_G_INFO%NE_COORD%LON,YD_G_INFO%NE_COORD%LAT,YD_G_INFO%MF_NE
IF (TKOUT >= 0) WRITE(TKOUT,'(1X,A7,1X,"|",2X,F7.2,2X,"|",2X,F7.2,1X,"|",1X,G16.10)') &
 & "N.West ",YD_G_INFO%NW_COORD%LON,YD_G_INFO%NW_COORD%LAT,YD_G_INFO%MF_NW
IF (TKOUT >= 0) WRITE(TKOUT,*) "============================================================="
IF (LHOOK) CALL DR_HOOK('EGGPACK:INFO_DOMI_PRINT',1,ZHOOK_HANDLE)
END SUBROUTINE INFO_DOMI_PRINT
! -------------------------------------------------------------------------------
SUBROUTINE INFO_PP_PRINT(P_P,KOUT,PI)
TYPE (PARAM_PROJ), INTENT(IN)                :: P_P
INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL     :: KOUT
REAL(KIND=JPRB), INTENT(IN), OPTIONAL        :: PI

REAL(KIND=JPRB)    :: TPI, DTR
INTEGER(KIND=JPIM) :: TKOUT
REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('EGGPACK:INFO_PP_PRINT',0,ZHOOK_HANDLE)
IF (PRESENT(KOUT))THEN
  TKOUT = KOUT
ELSE
  TKOUT = 6_JPIM
ENDIF
IF (PRESENT(PI)) THEN
  TPI = PI
ELSE
  TPI = ASIN(1.0_JPRB)*2.0_JPRB
ENDIF
DTR = TPI/180.0_JPRB
IF (TKOUT >= 0) WRITE(TKOUT,*) "============================================================="
IF (TKOUT >= 0) WRITE(TKOUT,*) "===  Informations about Parameters Projection Structure  ===="
IF (TKOUT >= 0) WRITE(TKOUT,*) "============================================================="
IF (TKOUT >= 0) WRITE(TKOUT,*)
IF (TKOUT >= 0) WRITE(TKOUT,*) "  -Reference Point Coordinates :"
IF (TKOUT >= 0) WRITE(TKOUT,'(13X,A7,10X,A7)') "Degrees","Radians"
IF (TKOUT >= 0) WRITE(TKOUT,'(1X,"Longitude : ",F7.2,5X,G16.10)') P_P%REF_PT%LON/DTR,P_P%REF_PT%LON
IF (TKOUT >= 0) WRITE(TKOUT,'(1X,"Latitude  : ",F7.2,5X,G16.10)') P_P%REF_PT%LAT/DTR,P_P%REF_PT%LAT
IF (TKOUT >= 0) WRITE(TKOUT,*) "  -Projection Characteristics  :"
IF (TKOUT >= 0) WRITE(TKOUT,'(13X,A16,5X,A18)') "   ERPK or KL   ","Type of Projection"
IF (TKOUT >= 0) WRITE(TKOUT,'(13X,G16.10,13X,A1)') P_P%KL,P_P%TYPE_PJ
IF ((P_P%TYPE_PJ == "M").OR.(P_P%TYPE_PJ == "W")) THEN
  IF (TKOUT >= 0) WRITE(TKOUT,*) "  -Rayon of Earth (in meters) :"
ELSE
  IF (TKOUT >= 0) WRITE(TKOUT,*) "  -Distance between Equator and Pole of projection on"
  IF (TKOUT >= 0) WRITE(TKOUT,*) "   projection plane (in meters) :"
ENDIF
IF (TKOUT >= 0) WRITE(TKOUT,'(13X,G16.10)') P_P%R_EQUATEUR
IF (TKOUT >= 0) WRITE(TKOUT,*) "  -Pole of projection (-1.0 for South, 1.0 for North, 0.0 for"
IF (TKOUT >= 0) WRITE(TKOUT,*) "   Mercator projection : no sense) :"
IF (TKOUT >= 0) WRITE(TKOUT,'(13X,F4.1)') P_P%POLE
IF (TKOUT >= 0) WRITE(TKOUT,*) "============================================================="
IF (LHOOK) CALL DR_HOOK('EGGPACK:INFO_PP_PRINT',1,ZHOOK_HANDLE)
END SUBROUTINE INFO_PP_PRINT
! -------------------------------------------------------------------------------
LOGICAL FUNCTION RETURN_PRINT(YD_CODE_ERR,K_NUM_TEST,AUTO_STOP,KOUT)
CHARACTER(LEN=36), PARAMETER                 :: FMT = '(A25,I4," ; test number = ",I3,A100)'
TYPE (ERROR), INTENT(IN)                     :: YD_CODE_ERR
INTEGER(KIND=JPIM), INTENT(IN)               :: K_NUM_TEST
LOGICAL, INTENT(IN), OPTIONAL                :: AUTO_STOP
INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL     :: KOUT

CHARACTER(LEN=25)  :: CL_ADD_TEXT
LOGICAL            :: TAS
INTEGER(KIND=JPIM) :: TKOUT
REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('EGGPACK:RETURN_PRINT',0,ZHOOK_HANDLE)
IF (PRESENT(KOUT))THEN
  TKOUT = KOUT
ELSE
  TKOUT = 6_JPIM
ENDIF
IF (PRESENT(AUTO_STOP))THEN
  TAS = AUTO_STOP
ELSE
  TAS = .TRUE.
ENDIF
SELECT CASE (YD_CODE_ERR%NUM)
CASE(:-1_JPIM) ; CL_ADD_TEXT = 'ERROR   : return value = '
CASE(0_JPIM)   ; CL_ADD_TEXT = 'OK      : return value = '
CASE(1_JPIM)   ; CL_ADD_TEXT = 'INFO    : return value = '
CASE(2_JPIM:)  ; CL_ADD_TEXT = 'WARNING : return value = '
END SELECT
IF (YD_CODE_ERR%NUM < 0_JPIM) WRITE (ABS(TKOUT),*) "Subroutine last status when aborted : "
IF (YD_CODE_ERR%NUM < 0_JPIM) THEN
  WRITE (ABS(TKOUT),FMT) CL_ADD_TEXT,YD_CODE_ERR%NUM,K_NUM_TEST,YD_CODE_ERR%TXT
  IF (TAS) THEN
    CALL ABOR1("Abort by EGGPACK:RETURN_PRINT") 
  ELSE
    RETURN_PRINT = .TRUE.
  ENDIF
ELSE
  IF (TKOUT >= 0) WRITE (TKOUT,FMT) CL_ADD_TEXT,YD_CODE_ERR%NUM,K_NUM_TEST,YD_CODE_ERR%TXT
  RETURN_PRINT = .FALSE.
ENDIF
IF (LHOOK) CALL DR_HOOK('EGGPACK:RETURN_PRINT',1,ZHOOK_HANDLE)
END FUNCTION RETURN_PRINT
! -------------------------------------------------------------------------------
FUNCTION TYPE_PROJ(REF_COORD) RESULT (TY_PJ)
! REF_COORD in Degrees
TYPE (LOLA), INTENT(IN)       :: REF_COORD
CHARACTER(LEN=1) :: TY_PJ
 
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('EGGPACK:TYPE_PROJ',0,ZHOOK_HANDLE)    
IF (REF_COORD%LAT == 0.0_JPRB) THEN
  TY_PJ = "M"
ELSEIF (ABS(REF_COORD%LAT) == 90.0_JPRB) THEN
  TY_PJ = "S"
ELSE
  TY_PJ = "L"
ENDIF
IF (LHOOK) CALL DR_HOOK('EGGPACK:TYPE_PROJ',1,ZHOOK_HANDLE)
END FUNCTION TYPE_PROJ
! -------------------------------------------------------------------------------
REAL(KIND=JPRB) FUNCTION POLE_IS (REF_COORD) RESULT (POL)
! REF_COORD in Degrees
TYPE (LOLA), INTENT(IN)         :: REF_COORD

REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('EGGPACK:POLE_IS',0,ZHOOK_HANDLE)
IF (REF_COORD%LAT == 0.0_JPRB) THEN
  POL = 0.0_JPRB
ELSE
  POL = SIGN(1.0_JPRB,REF_COORD%LAT)
ENDIF
IF (LHOOK) CALL DR_HOOK('EGGPACK:POLE_IS',1,ZHOOK_HANDLE)
END FUNCTION POLE_IS
! -------------------------------------------------------------------------------

! ******************* Specifics functions ***************************************

! STPL : STEREOGRAPHIQUE POLAIRE / LAMBERT MODE FONCTIONS 
! Transform LAT-LON to R-TETA, R-TETA to XY in STD Origin and reverse (-+R)
! -------------------------------------------------------------------------------
TYPE (RTETA) FUNCTION STPL_LATLON_TO_RTETA_S(PT_COORD,P_PJ,PI) RESULT (PT_RTETA)
! PT_COORD in Radians
TYPE (LOLA), INTENT(IN)                             :: PT_COORD
TYPE (PARAM_PROJ), INTENT(IN)                       :: P_PJ
REAL(KIND=JPRB), INTENT(IN), OPTIONAL               :: PI

REAL(KIND=JPRB) :: TPI
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('EGGPACK:STPL_LATLON_TO_RTETA_S',0,ZHOOK_HANDLE)
IF (PRESENT(PI)) THEN
  TPI = PI
ELSE
  TPI = ASIN(1.0_JPRB)*2.0_JPRB
ENDIF
PT_RTETA%R = P_PJ%R_EQUATEUR*((TAN((TPI/4.0_JPRB)-(P_PJ%POLE*PT_COORD%LAT/2.0_JPRB)))**(P_PJ%KL))
PT_RTETA%TETA = P_PJ%KL*DIST_2REF(PT_COORD,P_PJ%REF_PT,TPI)
IF (LHOOK) CALL DR_HOOK('EGGPACK:STPL_LATLON_TO_RTETA_S',1,ZHOOK_HANDLE)
END FUNCTION STPL_LATLON_TO_RTETA_S
! -------------------------------------------------------------------------------
TYPE (XY) FUNCTION STLP_RTETA_TO_XY_S(PT_RTETA,P_PJ) RESULT (PT_XY)
! PT_XY in STD Origin that is pole ( Y positive to North, X positive to East )
TYPE (RTETA), INTENT(IN)            :: PT_RTETA
TYPE (PARAM_PROJ), INTENT(IN)       :: P_PJ

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('EGGPACK:STLP_RTETA_TO_XY_S',0,ZHOOK_HANDLE)
PT_XY%X = PT_RTETA%R*SIN(PT_RTETA%TETA)
PT_XY%Y = -P_PJ%POLE*PT_RTETA%R*COS(PT_RTETA%TETA)
IF (LHOOK) CALL DR_HOOK('EGGPACK:STLP_RTETA_TO_XY_S',1,ZHOOK_HANDLE)
END FUNCTION STLP_RTETA_TO_XY_S
! -------------------------------------------------------------------------------
TYPE (RTETA) FUNCTION STLP_XY_TO_RTETA_S(PT_XY,P_PJ,PI) RESULT (PT_RTETA)
! PT_XY in STD Origin that is pole ( Y positive to North, X positive to East )
TYPE (XY), INTENT(IN)                     :: PT_XY
TYPE (PARAM_PROJ), INTENT(IN)             :: P_PJ
REAL(KIND=JPRB), INTENT(IN), OPTIONAL     :: PI

REAL(KIND=JPRB)     :: TPI
REAL(KIND=JPRD)     :: TATNG
REAL(KIND=JPHOOK)     :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('EGGPACK:STLP_XY_TO_RTETA_S',0,ZHOOK_HANDLE)
IF (PRESENT(PI)) THEN
  TPI = PI
ELSE
  TPI = ASIN(1.0_JPRB)*2.0_JPRB
ENDIF
PT_RTETA%R = SQRT((REAL(PT_XY%X,KIND=JPRD)*REAL(PT_XY%X,KIND=JPRD))+(REAL(PT_XY%Y,KIND=JPRD)*REAL(PT_XY%Y,KIND=JPRD)))
IF (PT_XY%Y == 0.0_JPRB) THEN
  IF (PT_XY%X == 0.0_JPRB) THEN
    TATNG = TPI
  ELSE
    TATNG = SIGN(TPI/2.0_JPRB,-P_PJ%POLE*PT_XY%X)
  ENDIF
ELSE
  TATNG = ATAN(-P_PJ%POLE*(REAL(PT_XY%X,KIND=JPRD)/REAL(PT_XY%Y,KIND=JPRD)))
ENDIF
PT_RTETA%TETA = TPI*SIGN(1.0_JPRB,PT_XY%X)*(SIGN(0.5_JPRB,P_PJ%POLE*PT_XY%Y)+0.5_JPRB)+TATNG
! This term : <-------------------------------------------------------------> modifies the value of atan() according
! the square defined by : y<0 => 0 ; x>0 and y>0 => pi ; x<0 and y>0 => -pi
!========>if PT_RTETA%TETA > TPI*P_PJ%KL then error (bad section of planed cone)
IF (ABS(PT_RTETA%TETA) > TPI*P_PJ%KL) THEN
  PRINT *,"Point at x = ",PT_XY%X," , y = ",PT_XY%Y
  PRINT *,"Is out of the planed cone section ! Abort !!!"
  CALL ABOR1("Abort by EGGPACK:STLP_XY_TO_RTETA_S")
ENDIF
IF (LHOOK) CALL DR_HOOK('EGGPACK:STLP_XY_TO_RTETA_S',1,ZHOOK_HANDLE)
END FUNCTION STLP_XY_TO_RTETA_S
! -------------------------------------------------------------------------------
TYPE (LOLA) FUNCTION STLP_RTETA_TO_LATLON_S(PT_RTETA,P_PJ,PI) RESULT (PT_COORD)
! PT_COORD in Radians
TYPE (RTETA), INTENT(IN)                    :: PT_RTETA
TYPE (PARAM_PROJ), INTENT(IN)               :: P_PJ
REAL(KIND=JPRB), INTENT(IN), OPTIONAL       :: PI

REAL(KIND=JPRB) :: TPI
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('EGGPACK:STLP_RTETA_TO_LATLON_S',0,ZHOOK_HANDLE)
IF (PRESENT(PI)) THEN
  TPI = PI
ELSE
  TPI = ASIN(1.0_JPRB)*2.0_JPRB
ENDIF
PT_COORD%LON = P_PJ%REF_PT%LON + PT_RTETA%TETA/P_PJ%KL
PT_COORD%LAT = P_PJ%POLE*((TPI/2.0_JPRB)-2.0_JPRB*ATAN((PT_RTETA%R/P_PJ%R_EQUATEUR)**(1.0_JPRD/P_PJ%KL)))
IF (LHOOK) CALL DR_HOOK('EGGPACK:STLP_RTETA_TO_LATLON_S',1,ZHOOK_HANDLE)
END FUNCTION STLP_RTETA_TO_LATLON_S
! -------------------------------------------------------------------------------
FUNCTION STPL_LATLON_TO_RTETA_V(PT_COORD,P_PJ,PI) RESULT (PT_RTETA)
! PT_COORD in Radians
TYPE (LOLA), DIMENSION(:), INTENT(IN)                :: PT_COORD
TYPE (PARAM_PROJ), INTENT(IN)                        :: P_PJ
REAL(KIND=JPRB), INTENT(IN), OPTIONAL                :: PI
TYPE (RTETA), DIMENSION(SIZE(PT_COORD)) :: PT_RTETA

REAL(KIND=JPRB) :: TPI
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('EGGPACK:STPL_LATLON_TO_RTETA_V',0,ZHOOK_HANDLE)
IF (PRESENT(PI)) THEN
  TPI = PI
ELSE
  TPI = ASIN(1.0_JPRB)*2.0_JPRB
ENDIF
PT_RTETA(:)%R = P_PJ%R_EQUATEUR*((TAN((TPI/4.0_JPRB)-(P_PJ%POLE*PT_COORD(:)%LAT/2.0_JPRD)))**P_PJ%KL)
PT_RTETA(:)%TETA = P_PJ%KL*DIST_2REF(PT_COORD(:),P_PJ%REF_PT,TPI)
IF (LHOOK) CALL DR_HOOK('EGGPACK:STPL_LATLON_TO_RTETA_V',1,ZHOOK_HANDLE)
END FUNCTION STPL_LATLON_TO_RTETA_V
! -------------------------------------------------------------------------------
FUNCTION STLP_RTETA_TO_XY_V(PT_RTETA,P_PJ) RESULT (PT_XY)
! PT_XY in STD Origin that is pole ( Y positive to North, X positive to East )
TYPE (RTETA), DIMENSION(:), INTENT(IN)            :: PT_RTETA
TYPE (PARAM_PROJ), INTENT(IN)                     :: P_PJ
TYPE (XY), DIMENSION(SIZE(PT_RTETA)) :: PT_XY

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('EGGPACK:STLP_RTETA_TO_XY_V',0,ZHOOK_HANDLE)
PT_XY(:)%X = PT_RTETA(:)%R*SIN(PT_RTETA(:)%TETA)
PT_XY(:)%Y = -P_PJ%POLE*PT_RTETA(:)%R*COS(PT_RTETA(:)%TETA)
IF (LHOOK) CALL DR_HOOK('EGGPACK:STLP_RTETA_TO_XY_V',1,ZHOOK_HANDLE)
END FUNCTION STLP_RTETA_TO_XY_V
! -------------------------------------------------------------------------------
FUNCTION STLP_XY_TO_RTETA_V(PT_XY,P_PJ,PI) RESULT (PT_RTETA)
! PT_XY in STD Origin that is pole ( Y positive to North, X positive to East )
TYPE (XY), DIMENSION(:), INTENT(IN)               :: PT_XY
TYPE (PARAM_PROJ), INTENT(IN)                     :: P_PJ
REAL(KIND=JPRB), INTENT(IN), OPTIONAL             :: PI
TYPE (RTETA), DIMENSION(SIZE(PT_XY)) :: PT_RTETA

REAL(KIND=JPRB)                         :: TPI
REAL(KIND=JPRD), DIMENSION(SIZE(PT_XY)) :: TATNG
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('EGGPACK:STLP_XY_TO_RTETA_V',0,ZHOOK_HANDLE)
IF (PRESENT(PI)) THEN
  TPI = PI
ELSE
  TPI = ASIN(1.0_JPRB)*2.0_JPRB
ENDIF
PT_RTETA(:)%R = SQRT((REAL(PT_XY(:)%X,KIND=JPRD)*REAL(PT_XY(:)%X,KIND=JPRD))+(REAL(PT_XY(:)%Y,KIND=JPRD)*REAL(PT_XY(:)%Y,KIND=JPRD)))
WHERE (PT_XY(:)%Y == 0.0_JPRB)
  WHERE (PT_XY(:)%X == 0.0_JPRB)
    TATNG = TPI
  ELSEWHERE
    TATNG = SIGN(TPI/2.0_JPRB,-P_PJ%POLE*PT_XY(:)%X)
  ENDWHERE
ELSEWHERE
  TATNG = ATAN(-P_PJ%POLE*(REAL(PT_XY(:)%X,KIND=JPRD)/REAL(PT_XY(:)%Y,KIND=JPRD)))
ENDWHERE
PT_RTETA(:)%TETA = TPI*SIGN(1.0_JPRB,PT_XY(:)%X)*(SIGN(0.5_JPRB,P_PJ%POLE*PT_XY(:)%Y)+0.5_JPRB)+TATNG(:)
! This term : <-------------------------------------------------------------------> modifies the value of atan()
! according the square defined by : y<0 => 0 ; x>0 and y>0 => pi ; x<0 and y>0 => -pi
!========>if one of PT_RTETA(:)%TETA > TPI*P_PJ%KL then error (bad section of planed cone)
IF (ANY(ABS(PT_RTETA(:)%TETA) > TPI*P_PJ%KL)) THEN
  PRINT *,"Some points are out of planed cone section ! Abort !!!"
  CALL ABOR1("Abort by EGGPACK:STLP_XY_TO_RTETA_V")       
ENDIF
IF (LHOOK) CALL DR_HOOK('EGGPACK:STLP_XY_TO_RTETA_V',1,ZHOOK_HANDLE)
END FUNCTION STLP_XY_TO_RTETA_V
! -------------------------------------------------------------------------------
FUNCTION STLP_RTETA_TO_LATLON_V(PT_RTETA,P_PJ,PI) RESULT (PT_COORD)
! PT_COORD in Radians
TYPE (RTETA), DIMENSION(:), INTENT(IN)              :: PT_RTETA
TYPE (PARAM_PROJ), INTENT(IN)                       :: P_PJ
REAL(KIND=JPRB), INTENT(IN), OPTIONAL               :: PI
TYPE (LOLA), DIMENSION(SIZE(PT_RTETA)) :: PT_COORD

REAL(KIND=JPRB) :: TPI
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('EGGPACK:STLP_RTETA_TO_LATLON_V',0,ZHOOK_HANDLE)
IF (PRESENT(PI)) THEN
  TPI = PI
ELSE
  TPI = ASIN(1.0_JPRB)*2.0_JPRB
ENDIF
PT_COORD(:)%LON = P_PJ%REF_PT%LON + (PT_RTETA(:)%TETA/P_PJ%KL)
PT_COORD(:)%LAT = P_PJ%POLE*((TPI/2.0_JPRB)-2.0_JPRB*ATAN((PT_RTETA(:)%R/P_PJ%R_EQUATEUR)**(1.0_JPRD/P_PJ%KL)))
IF (LHOOK) CALL DR_HOOK('EGGPACK:STLP_RTETA_TO_LATLON_V',1,ZHOOK_HANDLE)
END FUNCTION STLP_RTETA_TO_LATLON_V
! -------------------------------------------------------------------------------

! ******************* Generics functions ****************************************

! Creation of Structure Definition projection type (-+D)
! -------------------------------------------------------------------------------
TYPE (PARAM_PROJ) FUNCTION REF_DATAS (REF_COORD,RA,TOZERO_COORD,LRT) RESULT (P_P)
! REF_COORD,TOZERO_COORD in -+Degrees
TYPE (LOLA), INTENT(IN)                    :: REF_COORD
REAL(KIND=JPRB), INTENT(IN), OPTIONAL      :: RA
TYPE (LOLA), INTENT(IN), OPTIONAL          :: TOZERO_COORD
LOGICAL, INTENT(IN), OPTIONAL              :: LRT

REAL(KIND=JPRB) :: RT
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('EGGPACK:REF_DATAS',0,ZHOOK_HANDLE)
IF (PRESENT(RA))THEN
  RT = RA
ELSE
  RT = R_EARTH
ENDIF
! P_P%REF in -+Radians
P_P%REF_PT = LOLAR(ANGLE_DOMAIN(REF_COORD)) ! expect REF_COORD in Dg & put it in -+
IF (PRESENT(LRT)) THEN ! flag rot-tilt present
  IF ((LRT).AND.(P_P%REF_PT%LAT==0.0_JPRB)) THEN ! flag rot-tilt true and mercator
    P_P%TYPE_PJ = "W"
    IF (PRESENT(TOZERO_COORD)) THEN ! value of TOZERO_COORD
      P_P%TZO_PT = LOLAR(ANGLE_DOMAIN(TOZERO_COORD)) ! expect TOZERO_COORD in Dg & put it in -+
    ELSE ! default TOZERO_COORD value is (0,0)
      P_P%TZO_PT%LAT = 0.0_JPRB
      P_P%TZO_PT%LON = 0.0_JPRB
    ENDIF
  ELSE ! either flag false either not mercator so
    P_P%TYPE_PJ = TYPE_PROJ(REF_COORD) ! doing standard
    P_P%TZO_PT%LAT = -999.999_JPRB ! no need of TOZERO_COORD
    P_P%TZO_PT%LON = -999.999_JPRB
  ENDIF
ELSE ! not rot-tilt mode so
  P_P%TYPE_PJ = TYPE_PROJ(REF_COORD) ! doing standard
  P_P%TZO_PT%LAT = -999.999_JPRB ! no need of TOZERO_COORD
  P_P%TZO_PT%LON = -999.999_JPRB
ENDIF
P_P%POLE = POLE_IS(REF_COORD)
P_P%KL = P_P%POLE*SIN(P_P%REF_PT%LAT)
IF (P_P%KL /= 0.0_JPRB) THEN
  ! Rho of projection of equator on polar plane (Rho-Theta)
  P_P%R_EQUATEUR = RT*((COS(P_P%REF_PT%LAT))**(1.0_JPRB -P_P%KL))*(((1.0_JPRB +P_P%KL)**P_P%KL)/P_P%KL)
ELSE
  P_P%R_EQUATEUR = RT
ENDIF
IF (LHOOK) CALL DR_HOOK('EGGPACK:REF_DATAS',1,ZHOOK_HANDLE)
END FUNCTION REF_DATAS
! -------------------------------------------------------------------------------
! Change of origine between STD and NEW and reverse (-+D)
! -------------------------------------------------------------------------------
TYPE (XY) FUNCTION XY_NEW_TO_STD_ORIGIN_S(NEW_ORIGIN_COORD,PT_XY_IN_NEW_ORIGIN,P_PJ,PI) RESULT (PT_XY_IN_STD_ORIGIN)
! NEW_ORIGIN_COORD in -+Degrees
TYPE (LOLA), INTENT(IN)                               :: NEW_ORIGIN_COORD
TYPE (XY), INTENT(IN)                                 :: PT_XY_IN_NEW_ORIGIN
TYPE (PARAM_PROJ), INTENT(IN)                         :: P_PJ
REAL(KIND=JPRB), INTENT(IN), OPTIONAL                 :: PI

REAL(KIND=JPRB) :: TPI
TYPE (XY)       :: N_O_PT_XY
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('EGGPACK:XY_NEW_TO_STD_ORIGIN_S',0,ZHOOK_HANDLE)
IF (PRESENT(PI)) THEN
  TPI = PI
ELSE
  TPI = ASIN(1.0_JPRB)*2.0_JPRB
ENDIF
N_O_PT_XY = LATLON_TO_XY(LOLAR(NEW_ORIGIN_COORD),P_PJ,TPI)
PT_XY_IN_STD_ORIGIN%X = PT_XY_IN_NEW_ORIGIN%X+N_O_PT_XY%X
PT_XY_IN_STD_ORIGIN%Y = PT_XY_IN_NEW_ORIGIN%Y+N_O_PT_XY%Y
IF (LHOOK) CALL DR_HOOK('EGGPACK:XY_NEW_TO_STD_ORIGIN_S',1,ZHOOK_HANDLE)
END FUNCTION XY_NEW_TO_STD_ORIGIN_S
! -------------------------------------------------------------------------------
TYPE (XY) FUNCTION XY_STD_TO_NEW_ORIGIN_S(NEW_ORIGIN_COORD,PT_XY_IN_STD_ORIGIN,P_PJ,PI) RESULT (PT_XY_IN_NEW_ORIGIN)
! NEW_ORIGIN_COORD in -+Degrees
TYPE (LOLA), INTENT(IN)                             :: NEW_ORIGIN_COORD
TYPE (XY), INTENT(IN)                               :: PT_XY_IN_STD_ORIGIN
TYPE (PARAM_PROJ), INTENT(IN)                       :: P_PJ
REAL(KIND=JPRB), INTENT(IN), OPTIONAL               :: PI

REAL(KIND=JPRB) :: TPI
TYPE (XY)       :: N_O_PT_XY
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('EGGPACK:XY_STD_TO_NEW_ORIGIN_S',0,ZHOOK_HANDLE)
IF (PRESENT(PI)) THEN
  TPI = PI
ELSE
  TPI = ASIN(1.0_JPRB)*2.0_JPRB
ENDIF
N_O_PT_XY = LATLON_TO_XY(LOLAR(NEW_ORIGIN_COORD),P_PJ,TPI)
PT_XY_IN_NEW_ORIGIN%X = PT_XY_IN_STD_ORIGIN%X-N_O_PT_XY%X
PT_XY_IN_NEW_ORIGIN%Y = PT_XY_IN_STD_ORIGIN%Y-N_O_PT_XY%Y
IF (LHOOK) CALL DR_HOOK('EGGPACK:XY_STD_TO_NEW_ORIGIN_S',1,ZHOOK_HANDLE)
END FUNCTION XY_STD_TO_NEW_ORIGIN_S
! -------------------------------------------------------------------------------
FUNCTION XY_NEW_TO_STD_ORIGIN_V(NEW_ORIGIN_COORD,PT_XY_IN_NEW_ORIGIN,P_PJ,PI) RESULT (PT_XY_IN_STD_ORIGIN)
! NEW_ORIGIN_COORD in -+Degrees
TYPE (LOLA), INTENT(IN)                                      :: NEW_ORIGIN_COORD
TYPE (XY), DIMENSION(:), INTENT(IN)                          :: PT_XY_IN_NEW_ORIGIN
TYPE (PARAM_PROJ), INTENT(IN)                                :: P_PJ
REAL(KIND=JPRB), INTENT(IN), OPTIONAL                        :: PI
TYPE (XY), DIMENSION(SIZE(PT_XY_IN_NEW_ORIGIN)) :: PT_XY_IN_STD_ORIGIN

REAL(KIND=JPRB) :: TPI
TYPE (XY)       :: N_O_PT_XY
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('EGGPACK:XY_NEW_TO_STD_ORIGIN_V',0,ZHOOK_HANDLE)
IF (PRESENT(PI)) THEN
  TPI = PI
ELSE
  TPI = ASIN(1.0_JPRB)*2.0_JPRB
ENDIF
N_O_PT_XY = LATLON_TO_XY(LOLAR(NEW_ORIGIN_COORD),P_PJ,TPI)
PT_XY_IN_STD_ORIGIN(:)%X = PT_XY_IN_NEW_ORIGIN(:)%X+N_O_PT_XY%X
PT_XY_IN_STD_ORIGIN(:)%Y = PT_XY_IN_NEW_ORIGIN(:)%Y+N_O_PT_XY%Y
IF (LHOOK) CALL DR_HOOK('EGGPACK:XY_NEW_TO_STD_ORIGIN_V',1,ZHOOK_HANDLE)
END FUNCTION XY_NEW_TO_STD_ORIGIN_V
! -------------------------------------------------------------------------------
FUNCTION XY_STD_TO_NEW_ORIGIN_V(YL_NEW_ORIGIN_COORD,YL_PT_XY_IN_STD_ORIGIN,P_PJ,PI) RESULT (YD_PT_XY_IN_NEW_ORIGIN)
! YL_NEW_ORIGIN_COORD in -+Degrees
TYPE (LOLA), INTENT(IN)                                         :: YL_NEW_ORIGIN_COORD
TYPE (XY), DIMENSION(:), INTENT(IN)                             :: YL_PT_XY_IN_STD_ORIGIN
TYPE (PARAM_PROJ), INTENT(IN)                                   :: P_PJ
REAL(KIND=JPRB), INTENT(IN), OPTIONAL                           :: PI
TYPE (XY), DIMENSION(SIZE(YL_PT_XY_IN_STD_ORIGIN)) :: YD_PT_XY_IN_NEW_ORIGIN

REAL(KIND=JPRB) :: TPI
TYPE (XY)       :: YL_N_O_PT_XY
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('EGGPACK:XY_STD_TO_NEW_ORIGIN_V',0,ZHOOK_HANDLE)
IF (PRESENT(PI)) THEN
  TPI = PI
ELSE
  TPI = ASIN(1.0_JPRB)*2.0_JPRB
ENDIF
YL_N_O_PT_XY = LATLON_TO_XY(LOLAR(YL_NEW_ORIGIN_COORD),P_PJ,TPI)
YD_PT_XY_IN_NEW_ORIGIN(:)%X = YL_PT_XY_IN_STD_ORIGIN(:)%X-YL_N_O_PT_XY%X
YD_PT_XY_IN_NEW_ORIGIN(:)%Y = YL_PT_XY_IN_STD_ORIGIN(:)%Y-YL_N_O_PT_XY%Y
IF (LHOOK) CALL DR_HOOK('EGGPACK:XY_STD_TO_NEW_ORIGIN_V',1,ZHOOK_HANDLE)
END FUNCTION XY_STD_TO_NEW_ORIGIN_V
! -------------------------------------------------------------------------------
! Coordinates transforms between XY in STD Origin and LAT-LON (-+R)
! -------------------------------------------------------------------------------
TYPE (XY) FUNCTION LATLON_TO_XY_S(PT_COORD,P_PJ,PI) RESULT (PT_XY)
! PT_COORD in -+Radians
! PT_XY in STD Origin
TYPE (LOLA), INTENT(IN)                           :: PT_COORD
TYPE (PARAM_PROJ), INTENT(IN)                     :: P_PJ
REAL(KIND=JPRB), INTENT(IN), OPTIONAL             :: PI

REAL(KIND=JPRB) :: TPI
TYPE (LOLA)     :: PT_COORD2
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('EGGPACK:LATLON_TO_XY_S',0,ZHOOK_HANDLE)
IF (PRESENT(PI)) THEN
  TPI = PI
ELSE
  TPI = ASIN(1.0_JPRB)*2.0_JPRB
ENDIF
IF ((P_PJ%TYPE_PJ == "S").OR.(P_PJ%TYPE_PJ == "L")) THEN
  PT_XY = STLP_RTETA_TO_XY(STPL_LATLON_TO_RTETA(PT_COORD,P_PJ,TPI),P_PJ)
ELSE
  PT_COORD2 = PT_COORD
  IF (P_PJ%TYPE_PJ == "W") THEN
    PT_COORD2 = METILROT(P_PJ%REF_PT,P_PJ%TZO_PT,PT_COORD)
    PT_XY%X = PT_COORD2%LON*P_PJ%R_EQUATEUR
  ELSE
    PT_XY%X = P_PJ%R_EQUATEUR*DIST_2REF(PT_COORD,P_PJ%REF_PT,TPI)
  ENDIF
  PT_XY%Y = -P_PJ%R_EQUATEUR*LOG(TAN((TPI/4.0_JPRB)-(PT_COORD2%LAT/2.0_JPRB)))
ENDIF
IF (LHOOK) CALL DR_HOOK('EGGPACK:LATLON_TO_XY_S',1,ZHOOK_HANDLE)
END FUNCTION LATLON_TO_XY_S
! -------------------------------------------------------------------------------
TYPE (LOLA) FUNCTION XY_TO_LATLON_S(PT_XY,P_PJ,PI) RESULT (PT_COORD)
! PT_XY in STD Origin
! PT_COORD in -+Radians??
TYPE (XY), INTENT(IN)                           :: PT_XY
TYPE (PARAM_PROJ), INTENT(IN)                   :: P_PJ
REAL(KIND=JPRB), INTENT(IN), OPTIONAL           :: PI

REAL(KIND=JPRB) :: TPI
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('EGGPACK:XY_TO_LATLON_S',0,ZHOOK_HANDLE)
IF (PRESENT(PI)) THEN
  TPI = PI
ELSE
  TPI = ASIN(1.0_JPRB)*2.0_JPRB
ENDIF
IF ((P_PJ%TYPE_PJ == "S").OR.(P_PJ%TYPE_PJ == "L")) THEN
  PT_COORD = STLP_RTETA_TO_LATLON(STLP_XY_TO_RTETA(PT_XY,P_PJ,TPI),P_PJ,TPI)
ELSE
  PT_COORD%LON = (PT_XY%X/P_PJ%R_EQUATEUR)
  PT_COORD%LAT = (TPI/2.0_JPRB)-2.0_JPRB*ATAN(EXP(-(PT_XY%Y/P_PJ%R_EQUATEUR)))
  IF (P_PJ%TYPE_PJ == "W") THEN
    PT_COORD=MEROTIL(P_PJ%REF_PT,P_PJ%TZO_PT,PT_COORD)
  ELSE
    PT_COORD%LON = P_PJ%REF_PT%LON+PT_COORD%LON
  ENDIF
ENDIF
IF (LHOOK) CALL DR_HOOK('EGGPACK:XY_TO_LATLON_S',1,ZHOOK_HANDLE)
END FUNCTION XY_TO_LATLON_S
! -------------------------------------------------------------------------------
FUNCTION LATLON_TO_XY_V(PT_COORD,P_PJ,PI) RESULT (PT_XY)
! PT_COORD in -+Radians
! PT_XY in STD Origin
TYPE (LOLA), DIMENSION(:), INTENT(IN)             :: PT_COORD
TYPE (PARAM_PROJ), INTENT(IN)                     :: P_PJ
REAL(KIND=JPRB), INTENT(IN), OPTIONAL             :: PI
TYPE (XY), DIMENSION(SIZE(PT_COORD)) :: PT_XY

REAL(KIND=JPRB)                        :: TPI
TYPE (LOLA), DIMENSION(SIZE(PT_COORD)) :: PT_COORD2
REAL(KIND=JPHOOK)                        :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('EGGPACK:LATLON_TO_XY_V',0,ZHOOK_HANDLE)
IF (PRESENT(PI)) THEN
  TPI = PI
ELSE
  TPI = ASIN(1.0_JPRB)*2.0_JPRB
ENDIF
IF ((P_PJ%TYPE_PJ == "S").OR.(P_PJ%TYPE_PJ == "L")) THEN
  PT_XY = STLP_RTETA_TO_XY(STPL_LATLON_TO_RTETA(PT_COORD,P_PJ,TPI),P_PJ)
ELSE
  PT_COORD2(:) = PT_COORD(:)
  IF (P_PJ%TYPE_PJ == "W") THEN
    PT_COORD2(:) = METILROT(P_PJ%REF_PT,P_PJ%TZO_PT,PT_COORD(:))
    PT_XY(:)%X = PT_COORD2(:)%LON*P_PJ%R_EQUATEUR
  ELSE
    PT_XY(:)%X = P_PJ%R_EQUATEUR*DIST_2REF(PT_COORD(:),P_PJ%REF_PT,TPI)
  ENDIF
  PT_XY(:)%Y = -P_PJ%R_EQUATEUR*LOG(TAN((TPI/4.0_JPRB)-(PT_COORD2(:)%LAT/2.0_JPRB)))
ENDIF
IF (LHOOK) CALL DR_HOOK('EGGPACK:LATLON_TO_XY_V',1,ZHOOK_HANDLE)
END FUNCTION LATLON_TO_XY_V
! -------------------------------------------------------------------------------
FUNCTION XY_TO_LATLON_V(YL_PT_XY,P_PJ,PI) RESULT (PT_COORD)
! YL_PT_XY in STD Origin
! PT_COORD in -+Radians??
TYPE (XY), DIMENSION(:), INTENT(IN)                 :: YL_PT_XY
TYPE (PARAM_PROJ), INTENT(IN)                       :: P_PJ
REAL(KIND=JPRB), INTENT(IN), OPTIONAL               :: PI
TYPE (LOLA), DIMENSION(SIZE(YL_PT_XY)) :: PT_COORD

REAL(KIND=JPRB) :: TPI
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('EGGPACK:XY_TO_LATLON_V',0,ZHOOK_HANDLE)
IF (PRESENT(PI)) THEN
  TPI = PI
ELSE
  TPI = ASIN(1.0_JPRB)*2.0_JPRB
ENDIF
IF ((P_PJ%TYPE_PJ == "S").OR.(P_PJ%TYPE_PJ == "L")) THEN
  PT_COORD = STLP_RTETA_TO_LATLON(STLP_XY_TO_RTETA(YL_PT_XY,P_PJ,TPI),P_PJ,TPI)
ELSE
  PT_COORD(:)%LON = (YL_PT_XY(:)%X/P_PJ%R_EQUATEUR)
  PT_COORD(:)%LAT = (TPI/2.0_JPRB)-2.0_JPRB*ATAN(EXP(-(YL_PT_XY(:)%Y/P_PJ%R_EQUATEUR)))
  IF (P_PJ%TYPE_PJ == "W") THEN
    PT_COORD(:)=MEROTIL(P_PJ%REF_PT,P_PJ%TZO_PT,PT_COORD(:))
  ELSE
    PT_COORD(:)%LON = P_PJ%REF_PT%LON+PT_COORD(:)%LON
  ENDIF
ENDIF
IF (LHOOK) CALL DR_HOOK('EGGPACK:XY_TO_LATLON_V',1,ZHOOK_HANDLE)
END FUNCTION XY_TO_LATLON_V
! -------------------------------------------------------------------------------
! Functions MAP_FACTOR and GN
! -------------------------------------------------------------------------------
REAL(KIND=JPRB) FUNCTION MAP_FACTOR_S(PT_COORD,P_PJ,PI,RA)
! PT_COORD in -+Radians??
TYPE (LOLA), INTENT(IN)                    :: PT_COORD
TYPE (PARAM_PROJ), INTENT(IN)              :: P_PJ
REAL(KIND=JPRB), INTENT(IN), OPTIONAL      :: RA, PI

REAL(KIND=JPRB) :: RT, TPI
TYPE (RTETA)    :: PT_RTETA
TYPE (LOLA)     :: PT_COORD2
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('EGGPACK:MAP_FACTOR_S',0,ZHOOK_HANDLE)
IF (PRESENT(RA))THEN
  RT = RA
ELSE
  RT = R_EARTH
ENDIF
IF (PRESENT(PI)) THEN
  TPI = PI
ELSE
  TPI = ASIN(1.0_JPRB)*2.0_JPRB
ENDIF
SELECT CASE(P_PJ%TYPE_PJ)
CASE('W') ; PT_COORD2 = METILROT(P_PJ%REF_PT,P_PJ%TZO_PT,PT_COORD)
  MAP_FACTOR_S = 1.0_JPRB/(COS(PT_COORD2%LAT))
CASE('M') ; MAP_FACTOR_S = 1.0_JPRB/(COS(PT_COORD%LAT))
CASE('S') ; MAP_FACTOR_S = 2.0_JPRB/(1.0_JPRB +P_PJ%POLE*SIN(PT_COORD%LAT))
CASE('L') ; PT_RTETA = STPL_LATLON_TO_RTETA(PT_COORD,P_PJ,TPI)
  MAP_FACTOR_S = (P_PJ%KL*PT_RTETA%R)/(RT*COS(PT_COORD%LAT))
END SELECT
IF (LHOOK) CALL DR_HOOK('EGGPACK:MAP_FACTOR_S',1,ZHOOK_HANDLE)
END FUNCTION MAP_FACTOR_S
! -------------------------------------------------------------------------------
FUNCTION MAP_FACTOR_V(PT_COORD,P_PJ,PI,RA)
! PT_COORD in -+Radians??
TYPE (LOLA), DIMENSION(:), INTENT(IN)                    :: PT_COORD
TYPE (PARAM_PROJ), INTENT(IN)                            :: P_PJ
REAL(KIND=JPRB), INTENT(IN), OPTIONAL                    :: RA, PI
REAL(KIND=JPRB), DIMENSION(SIZE(PT_COORD))  :: MAP_FACTOR_V

REAL(KIND=JPRB)                         :: RT, TPI
TYPE (RTETA), DIMENSION(SIZE(PT_COORD)) :: YL_PT_RTETA
TYPE (LOLA), DIMENSION(SIZE(PT_COORD))  :: PT_COORD2
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('EGGPACK:MAP_FACTOR_V',0,ZHOOK_HANDLE)
IF (PRESENT(RA))THEN
  RT = RA
ELSE
  RT = R_EARTH
ENDIF
IF (PRESENT(PI)) THEN
  TPI = PI
ELSE
  TPI = ASIN(1.0_JPRB)*2.0_JPRB
ENDIF
SELECT CASE(P_PJ%TYPE_PJ)
CASE('W') ; PT_COORD2(:) = METILROT(P_PJ%REF_PT,P_PJ%TZO_PT,PT_COORD(:))
  MAP_FACTOR_V(:) = 1.0_JPRB/(COS(PT_COORD2(:)%LAT))
CASE('M') ; MAP_FACTOR_V(:) = 1.0_JPRB/(COS(PT_COORD(:)%LAT))
CASE('S') ; MAP_FACTOR_V(:) = 2.0_JPRB/(1.0_JPRB +P_PJ%POLE*SIN(PT_COORD(:)%LAT))
CASE('L') ; YL_PT_RTETA     = STPL_LATLON_TO_RTETA(PT_COORD,P_PJ,TPI)
  MAP_FACTOR_V(:) = (P_PJ%KL*YL_PT_RTETA(:)%R)/(RT*COS(PT_COORD(:)%LAT))
END SELECT
IF (LHOOK) CALL DR_HOOK('EGGPACK:MAP_FACTOR_V',1,ZHOOK_HANDLE)
END FUNCTION MAP_FACTOR_V
! -------------------------------------------------------------------------------
TYPE (PGN) FUNCTION GN_S(PT_COORD,P_PJ)
! PT_COORD in -+Radians??
TYPE (LOLA), INTENT(IN)                           :: PT_COORD
TYPE (PARAM_PROJ), INTENT(IN)                     :: P_PJ

TYPE (LOLA)     :: PT_COORD2
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('EGGPACK:GN_S',0,ZHOOK_HANDLE)
SELECT CASE(P_PJ%TYPE_PJ)
CASE('W') ; PT_COORD2 = METILROT(P_PJ%REF_PT,P_PJ%TZO_PT,PT_COORD)
  GN_S%ONY = (((COS(P_PJ%REF_PT%LON)*((COS(P_PJ%TZO_PT%LAT)*COS(PT_COORD%LAT))+ &
   & (SIN(P_PJ%TZO_PT%LAT)*SIN(PT_COORD%LAT)*COS(PT_COORD%LON-P_PJ%TZO_PT%LON)))) &
   & -(SIN(P_PJ%REF_PT%LON)*SIN(PT_COORD%LAT)*SIN(PT_COORD%LON-P_PJ%TZO_PT%LON))) &
   & /COS(PT_COORD2%LAT))
  GN_S%ONX = -((COS(P_PJ%REF_PT%LON)*SIN(P_PJ%TZO_PT%LAT)*SIN(PT_COORD%LON-P_PJ%TZO_PT%LON)) &
   & +(SIN(P_PJ%REF_PT%LON)*COS(PT_COORD%LON-P_PJ%TZO_PT%LON)))/COS(PT_COORD2%LAT)
CASE('M') ; GN_S%ONX = 0.0_JPRB
  GN_S%ONY = 1.0_JPRB
CASE DEFAULT ; GN_S%ONX = -P_PJ%POLE*SIN(P_PJ%KL*(PT_COORD%LON-P_PJ%REF_PT%LON))
  GN_S%ONY = COS(P_PJ%KL*(PT_COORD%LON-P_PJ%REF_PT%LON))
END SELECT
IF (LHOOK) CALL DR_HOOK('EGGPACK:GN_S',1,ZHOOK_HANDLE)
END FUNCTION GN_S
! -------------------------------------------------------------------------------
FUNCTION GN_V(YD_PT_COORD,YD_P_PJ)
! YD_PT_COORD in -+Radians??
TYPE (LOLA), DIMENSION(:), INTENT(IN)                      :: YD_PT_COORD
TYPE (PARAM_PROJ), INTENT(IN)                              :: YD_P_PJ
TYPE (PGN), DIMENSION(SIZE(YD_PT_COORD))      :: GN_V

TYPE (LOLA), DIMENSION(SIZE(YD_PT_COORD)) :: YD_PT_COORD2
REAL(KIND=JPHOOK)                           :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('EGGPACK:GN_V',0,ZHOOK_HANDLE)
SELECT CASE(YD_P_PJ%TYPE_PJ)
CASE('W') ; YD_PT_COORD2(:) = METILROT(YD_P_PJ%REF_PT,YD_P_PJ%TZO_PT,YD_PT_COORD(:))
  GN_V(:)%ONY = (((COS(YD_P_PJ%REF_PT%LON)*((COS(YD_P_PJ%TZO_PT%LAT)*COS(YD_PT_COORD(:)%LAT))+ &
   & (SIN(YD_P_PJ%TZO_PT%LAT)*SIN(YD_PT_COORD(:)%LAT)*COS(YD_PT_COORD(:)%LON-YD_P_PJ%TZO_PT%LON)))) &
   & -(SIN(YD_P_PJ%REF_PT%LON)*SIN(YD_PT_COORD(:)%LAT)*SIN(YD_PT_COORD(:)%LON-YD_P_PJ%TZO_PT%LON))) &
   & /COS(YD_PT_COORD2(:)%LAT))
  GN_V(:)%ONX = -((COS(YD_P_PJ%REF_PT%LON)*SIN(YD_P_PJ%TZO_PT%LAT)*SIN(YD_PT_COORD(:)%LON-YD_P_PJ%TZO_PT%LON)) &
   & +(SIN(YD_P_PJ%REF_PT%LON)*COS(YD_PT_COORD(:)%LON-YD_P_PJ%TZO_PT%LON)))/COS(YD_PT_COORD2(:)%LAT)
CASE('M') ; GN_V(:)%ONX = 0.0_JPRB
  GN_V(:)%ONY = 1.0_JPRB
CASE DEFAULT ; GN_V(:)%ONX = -YD_P_PJ%POLE*SIN(YD_P_PJ%KL*(YD_PT_COORD(:)%LON-YD_P_PJ%REF_PT%LON))
  GN_V(:)%ONY = COS(YD_P_PJ%KL*(YD_PT_COORD(:)%LON-YD_P_PJ%REF_PT%LON))
END SELECT
IF (LHOOK) CALL DR_HOOK('EGGPACK:GN_V',1,ZHOOK_HANDLE)
END FUNCTION GN_V
! -------------------------------------------------------------------------------

! =================== SUBROUTINE ================================================

! ******************* subroutine to make grid domain ****************************

! -------------------------------------------------------------------------------
SUBROUTINE MAKDO(YD_REF_COORD,YD_CENTER_COORD,YD_PDEL,YD_NB_PTS,YD_GRID_COORD,P_GRID_MF, &
 & YD_GRID_PGN,YD_GRID_INFO,YD_ERR_CODE,LD_LIP,LD_AUTO_STOP,PI,P_RA,KOUT,LD_LMRT)
! Input Coordinates are in -+Degrees, output in 0+Radians excepted in YD_GRID_INFO : 0+D sauf Ref & Cen -+D
TYPE (LOLA), INTENT(INOUT)                                            :: YD_REF_COORD, YD_CENTER_COORD
TYPE (DELTA), INTENT(IN)                                              :: YD_PDEL
TYPE (NBPTS), INTENT(IN)                                              :: YD_NB_PTS
REAL(KIND=JPRB), INTENT(IN), OPTIONAL                                 :: PI, P_RA
LOGICAL, INTENT(IN), OPTIONAL                                         :: LD_AUTO_STOP, LD_LIP, LD_LMRT
TYPE (LOLA), DIMENSION(YD_NB_PTS%ONX,YD_NB_PTS%ONY), INTENT(OUT)      :: YD_GRID_COORD
REAL(KIND=JPRB), DIMENSION(YD_NB_PTS%ONX,YD_NB_PTS%ONY), INTENT(OUT)  :: P_GRID_MF
TYPE (PGN), DIMENSION(YD_NB_PTS%ONX,YD_NB_PTS%ONY), INTENT(OUT)       :: YD_GRID_PGN
TYPE (DOMI), INTENT(OUT)                                              :: YD_GRID_INFO
TYPE (ERROR), INTENT(OUT)                                             :: YD_ERR_CODE
INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL                              :: KOUT
TYPE (ERROR), DIMENSION(-6:4)                     :: YL_TAB_ERR_DEF = (/                                            &
 & ERROR (-6_JPIM, " : In Mercator, Deformations are too big out of [-85.0,85.0] (Rotated/Tilted mode)!"),           &
 & ERROR (-5_JPIM, " : In Lambert, the pole must be out of the domain !"),                                           &
 & ERROR (-4_JPIM, " : In Mercator, Deformations are too big out of [-85.0,85.0] !"),                                &
 & ERROR (-3_JPIM, " : Center point degrees coordinates out of bounds !"),                                           &
 & ERROR (-2_JPIM, " : Reference point degrees coordinates out of bounds !"),                                        &
 & ERROR (-1_JPIM, " : Subroutine aborted, sorry ..."),                                                              & 
 & ERROR ( 0_JPIM, " : Subroutine finished successly."),                                                             &
 & ERROR ( 1_JPIM, " : Test OK, subroutine go ahead !"),                                                             &
 & ERROR ( 2_JPIM, " : In Mercator, It's better to use Lambert or St.Pol.  if ABS(CENTER_LAT) > 20.0 !"),            &
 & ERROR ( 3_JPIM, " : In St.Pol. , It's better to use Lambert or Mercator if ABS(CENTER_LAT) < 70.0 !"),            &
 & ERROR ( 4_JPIM, " : In Lambert , It's better to use St.Pol. or Mercator if ABS(REF_LAT) is out of [20.0,70.0] !") &
 & /)
REAL(KIND=JPRB)                                   :: Z_RT, Z_TPI
LOGICAL                                           :: LL_TAS, LL_TLIP, LL_TLMRT
INTEGER(KIND=JPIM)                                :: I, J, I_TKOUT
TYPE (PARAM_PROJ)                                 :: YL_P_P
TYPE (XY), DIMENSION(YD_NB_PTS%ONX,YD_NB_PTS%ONY) :: YL_GRID_XY_C, YL_GRID_XY_P
REAL(KIND=JPHOOK)                                   :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('EGGPACK:MAKDO',0,ZHOOK_HANDLE)
IF (PRESENT(LD_LIP))THEN
  LL_TLIP = LD_LIP
ELSE
  LL_TLIP = _DEFP_
ENDIF
IF (PRESENT(LD_LMRT))THEN
  LL_TLMRT = LD_LMRT
ELSE
  LL_TLMRT = .FALSE.
ENDIF
IF (PRESENT(KOUT))THEN
  I_TKOUT = KOUT
ELSE
  I_TKOUT = 6_JPIM
ENDIF
IF (PRESENT(LD_AUTO_STOP))THEN
  LL_TAS = LD_AUTO_STOP
ELSE
  LL_TAS = .TRUE.
ENDIF
IF (PRESENT(P_RA))THEN
  Z_RT = P_RA
ELSE
  Z_RT = R_EARTH
ENDIF
IF (PRESENT(PI)) THEN
  Z_TPI = PI
ELSE
  Z_TPI = ASIN(1.0_JPRB)*2.0_JPRB
ENDIF
IF (I_TKOUT >=0 ) WRITE(I_TKOUT,*) "Begining of Makin Domain MAKDO subroutine"
! Validation of inputs LAT_LON coordinates or init YD_ERR_CODE to no error
YD_ERR_CODE = YL_TAB_ERR_DEF(VAL_COORD(YD_REF_COORD,-2_JPIM,Z_TPI))
TEST_1: IF (RETURN_PRINT(YD_ERR_CODE,1_JPIM,LL_TAS,I_TKOUT)) THEN
  IF (LHOOK) CALL DR_HOOK('EGGPACK:MAKDO',1,ZHOOK_HANDLE) 
  RETURN
ENDIF TEST_1
YD_ERR_CODE = YL_TAB_ERR_DEF(VAL_COORD(YD_CENTER_COORD,-3_JPIM,Z_TPI))
TEST_2: IF (RETURN_PRINT(YD_ERR_CODE,2_JPIM,LL_TAS,I_TKOUT)) THEN
  IF (LHOOK) CALL DR_HOOK('EGGPACK:MAKDO',1,ZHOOK_HANDLE) 
  RETURN
ENDIF TEST_2
! Compute parameters of projection and init of domain structures
YD_GRID_INFO%CT_COORD  = YD_CENTER_COORD
YD_GRID_INFO%RF_COORD  = YD_REF_COORD
YL_P_P                 = REF_DATAS(YD_REF_COORD,Z_RT,YD_CENTER_COORD,LL_TLMRT)
YD_GRID_INFO%INFO_PROJ = YL_P_P
! Degrees to Radians
YD_CENTER_COORD        = LOLAR(YD_CENTER_COORD)
YD_REF_COORD           = YL_P_P%REF_PT
! Compute XY grid points under CENTER origin
DO I=1_JPIM,YD_NB_PTS%ONX
  DO J=1_JPIM,YD_NB_PTS%ONY
    YL_GRID_XY_C(I,J)%X = (REAL(I,KIND=JPRB)-(REAL(YD_NB_PTS%ONX+1_JPIM,KIND=JPRB)/2.0_JPRB))*YD_PDEL%ONX
    YL_GRID_XY_C(I,J)%Y = (REAL(J,KIND=JPRB)-(REAL(YD_NB_PTS%ONY+1_JPIM,KIND=JPRB)/2.0_JPRB))*YD_PDEL%ONY
  ENDDO
ENDDO
IF (YL_P_P%TYPE_PJ/='W') THEN
  ! Change XY coordinates in CENTER Origin in STD Origin
  DO J=1_JPIM,YD_NB_PTS%ONY
    YL_GRID_XY_P(1:YD_NB_PTS%ONX,J) = &
     & XY_NEW_TO_STD_ORIGIN(YD_GRID_INFO%CT_COORD,YL_GRID_XY_C(1:YD_NB_PTS%ONX,J),YL_P_P,Z_TPI)
  ENDDO
ELSE
  YL_GRID_XY_P = YL_GRID_XY_C
ENDIF
! Validation depending projection type
SELECT CASE (YL_P_P%TYPE_PJ)
CASE('M')
  ! Validation of Mercator projection
  IF (ABS(YD_GRID_INFO%CT_COORD%LAT) > 20.0_JPRB) YD_ERR_CODE = YL_TAB_ERR_DEF(2_JPIM)
  TEST_3: IF (RETURN_PRINT(YD_ERR_CODE,3_JPIM,LL_TAS,I_TKOUT)) THEN
    IF (LHOOK) CALL DR_HOOK('EGGPACK:MAKDO',1,ZHOOK_HANDLE) 
    RETURN
  ENDIF TEST_3
  IF (MAXVAL(ABS(YL_GRID_XY_P(:,:)%Y)) > 3.0_JPRB*Z_RT) YD_ERR_CODE = YL_TAB_ERR_DEF(-4_JPIM)
  TEST_4: IF (RETURN_PRINT(YD_ERR_CODE,4_JPIM,LL_TAS,I_TKOUT)) THEN
    IF (LHOOK) CALL DR_HOOK('EGGPACK:MAKDO',1,ZHOOK_HANDLE) 
    RETURN
  ENDIF TEST_4
CASE('S')
  ! Validation of Polar Stereographic projection
  IF (ABS(YD_GRID_INFO%CT_COORD%LAT) < 70.0_JPRB) YD_ERR_CODE = YL_TAB_ERR_DEF(3_JPIM)
  TEST_5: IF (RETURN_PRINT(YD_ERR_CODE,5_JPIM,LL_TAS,I_TKOUT)) THEN
    IF (LHOOK) CALL DR_HOOK('EGGPACK:MAKDO',1,ZHOOK_HANDLE) 
    RETURN
  ENDIF TEST_5
CASE('L')
  ! Validation of Lambert projection
  IF ((ABS(YD_GRID_INFO%RF_COORD%LAT) < 20.0_JPRB).OR.(ABS(YD_GRID_INFO%RF_COORD%LAT) > 70.0_JPRB)) &
   & YD_ERR_CODE = YL_TAB_ERR_DEF(4_JPIM)
  TEST_6: IF (RETURN_PRINT(YD_ERR_CODE,6_JPIM,LL_TAS,I_TKOUT)) THEN
    IF (LHOOK) CALL DR_HOOK('EGGPACK:MAKDO',1,ZHOOK_HANDLE) 
    RETURN
  ENDIF TEST_6
  IF (((YL_GRID_XY_P(1,1)%X*YL_GRID_XY_P(YD_NB_PTS%ONX,1)%X) < 0.0_JPRB).AND. &
   & ((YL_GRID_XY_P(1,1)%Y*YL_GRID_XY_P(1,YD_NB_PTS%ONY)%Y) < 0.0_JPRB)) YD_ERR_CODE = YL_TAB_ERR_DEF(-5_JPIM)
  TEST_7: IF (RETURN_PRINT(YD_ERR_CODE,7_JPIM,LL_TAS,I_TKOUT)) THEN
    IF (LHOOK) CALL DR_HOOK('EGGPACK:MAKDO',1,ZHOOK_HANDLE) 
    RETURN
  ENDIF TEST_7
CASE('W')
  ! Validation of Mercator projection in rotated/tilted case
  IF (MAXVAL(ABS(YL_GRID_XY_P(:,:)%Y)) > 3.0_JPRB*Z_RT) YD_ERR_CODE = YL_TAB_ERR_DEF(-6_JPIM)
  TEST_8: IF (RETURN_PRINT(YD_ERR_CODE,8_JPIM,LL_TAS,I_TKOUT)) THEN
    IF (LHOOK) CALL DR_HOOK('EGGPACK:MAKDO',1,ZHOOK_HANDLE) 
    RETURN
  ENDIF TEST_8
END SELECT
! Compute ouputs datas depending projection type
DO J=1_JPIM,YD_NB_PTS%ONY
  YD_GRID_COORD(1:YD_NB_PTS%ONX,J) = XY_TO_LATLON(YL_GRID_XY_P(1:YD_NB_PTS%ONX,J),YL_P_P,Z_TPI)
  P_GRID_MF(1:YD_NB_PTS%ONX,J)     = MAP_FACTOR(YD_GRID_COORD(1:YD_NB_PTS%ONX,J),YL_P_P,Z_TPI,Z_RT)
  YD_GRID_PGN(1:YD_NB_PTS%ONX,J)   = GN(YD_GRID_COORD(1:YD_NB_PTS%ONX,J),YL_P_P)
  YD_GRID_COORD(1:YD_NB_PTS%ONX,J) = ANGLE_DOMAIN(YD_GRID_COORD(1:YD_NB_PTS%ONX,J),Z_TPI,'0+','R')
ENDDO
! Fill info structure
YD_GRID_INFO%G_SIZE   = YD_NB_PTS
YD_GRID_INFO%SW_COORD = LOLAD(YD_GRID_COORD(1,1))
YD_GRID_INFO%MF_SW    = P_GRID_MF(1,1)
YD_GRID_INFO%SE_COORD = LOLAD(YD_GRID_COORD(YD_NB_PTS%ONX,1))
YD_GRID_INFO%MF_SE    = P_GRID_MF(YD_NB_PTS%ONX,1)
YD_GRID_INFO%NE_COORD = LOLAD(YD_GRID_COORD(YD_NB_PTS%ONX,YD_NB_PTS%ONY))
YD_GRID_INFO%MF_NE    = P_GRID_MF(YD_NB_PTS%ONX,YD_NB_PTS%ONY)
YD_GRID_INFO%NW_COORD = LOLAD(YD_GRID_COORD(1,YD_NB_PTS%ONY))
YD_GRID_INFO%MF_NW    = P_GRID_MF(1,YD_NB_PTS%ONY)
YD_GRID_INFO%MF_RF    = MAP_FACTOR(YD_REF_COORD,YL_P_P,Z_TPI,Z_RT)
YD_GRID_INFO%MF_CT    = MAP_FACTOR(YD_CENTER_COORD,YL_P_P,Z_TPI,Z_RT)
YD_REF_COORD          = ANGLE_DOMAIN(YD_REF_COORD,Z_TPI,'0+','R')
YD_CENTER_COORD       = ANGLE_DOMAIN(YD_CENTER_COORD,Z_TPI,'0+','R')
IF (LL_TLIP) CALL INFO_PRINT(YD_GRID_INFO,I_TKOUT,Z_TPI)
IF (I_TKOUT >=0) WRITE (I_TKOUT,*) "Subroutine last status when finished : "
IF (YD_ERR_CODE%NUM == 1_JPIM) YD_ERR_CODE = YL_TAB_ERR_DEF(0_JPIM)
TEST_9: IF (RETURN_PRINT(YD_ERR_CODE,9_JPIM,LL_TAS,I_TKOUT)) THEN
  IF (LHOOK) CALL DR_HOOK('EGGPACK:MAKDO',1,ZHOOK_HANDLE) 
  RETURN
ENDIF TEST_9
IF (LHOOK) CALL DR_HOOK('EGGPACK:MAKDO',1,ZHOOK_HANDLE)
END SUBROUTINE MAKDO

#undef _DEFP_
! -------------------------------------------------------------------------------
END MODULE EGGPACK
