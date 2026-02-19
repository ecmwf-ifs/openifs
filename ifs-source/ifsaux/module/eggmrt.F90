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

MODULE EGGMRT

! Version 2009.0317 by JD GRIL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! DOC !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Tilt and Rotate routines and reverse (see Pierre Benard Document.)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Author : Jean-Daniel GRIL , CNRM/GMAP/COOPE , October 19, 2004
! Modifs : JD Gril : March 2009 : Optimisation of vectirization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! ******************* Definition of parameters **********************************

! Include Kinds
! -------------
USE PARKIND1  ,ONLY : JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

! ******************* Loading module ********************************************

USE EGGANGLES, ONLY : LOLA, COSIN_TO_ANGLE, P_ASIN

IMPLICIT NONE

! ******************* Definition of type ****************************************

! ******************* Definition of Interface ***********************************

INTERFACE TILT
  MODULE PROCEDURE TILT_V, TILT_S
END INTERFACE

INTERFACE ROTATE
  MODULE PROCEDURE ROTATE_V, ROTATE_S
END INTERFACE

INTERFACE ANTI_TILT
  MODULE PROCEDURE ANTI_TILT_V, ANTI_TILT_S
END INTERFACE

INTERFACE ANTI_ROTATE
  MODULE PROCEDURE ANTI_ROTATE_V, ANTI_ROTATE_S
END INTERFACE

INTERFACE MEROTIL
  MODULE PROCEDURE MEROTIL_V, MEROTIL_S
END INTERFACE

INTERFACE METILROT
  MODULE PROCEDURE METILROT_V, METILROT_S
END INTERFACE

CONTAINS

! =================== FUNCTIONS =================================================

! ******************* Specifics functions ***************************************
! -------------------------------------------------------------------------------
! Function to TILT Scalar
! -------------------------------------------------------------------------------
TYPE(LOLA) FUNCTION TILT_S(REF_COORD,PT_COORD2) RESULT (PT_COORD1)
TYPE (LOLA), INTENT(IN)          :: REF_COORD, PT_COORD2

REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('EGGMRT:TILT_S',0,ZHOOK_HANDLE)    
PT_COORD1%LAT = P_ASIN(COS(REF_COORD%LON)*SIN(PT_COORD2%LAT)-SIN(REF_COORD%LON)* &
 & COS(PT_COORD2%LAT)*SIN(PT_COORD2%LON))
IF (COS(PT_COORD1%LAT) /= 0.0_JPRB) THEN
  PT_COORD1%LON = COSIN_TO_ANGLE((COS(PT_COORD2%LAT)*COS(PT_COORD2%LON))/COS(PT_COORD1%LAT), &
   & (SIN(REF_COORD%LON)*SIN(PT_COORD2%LAT)+COS(REF_COORD%LON)* &
   & COS(PT_COORD2%LAT)*SIN(PT_COORD2%LON))/COS(PT_COORD1%LAT))
ELSE
  PT_COORD1%LON = 0.0_JPRB
ENDIF
IF (LHOOK) CALL DR_HOOK('EGGMRT:TILT_S',1,ZHOOK_HANDLE) 
END FUNCTION TILT_S
! -------------------------------------------------------------------------------
! Function to TILT Vector
! -------------------------------------------------------------------------------
FUNCTION TILT_V(REF_COORD,PT_COORD2) RESULT (PT_COORD1)
TYPE (LOLA), INTENT(IN)                              :: REF_COORD
TYPE (LOLA), DIMENSION(:),INTENT(IN)                 :: PT_COORD2
TYPE (LOLA), DIMENSION(SIZE(PT_COORD2)) :: PT_COORD1

REAL(KIND=JPHOOK)                         :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('EGGMRT:TILT_V',0,ZHOOK_HANDLE)

PT_COORD1(:)%LAT = P_ASIN(COS(REF_COORD%LON)*SIN(PT_COORD2(:)%LAT)-SIN(REF_COORD%LON)* &
 & COS(PT_COORD2(:)%LAT)*SIN(PT_COORD2(:)%LON))
WHERE (COS(PT_COORD1(:)%LAT) /= 0.0_JPRB)
  PT_COORD1(:)%LON = COSIN_TO_ANGLE((COS(PT_COORD2(:)%LAT)*COS(PT_COORD2(:)%LON))/COS(PT_COORD1(:)%LAT), &
   & (SIN(REF_COORD%LON)*SIN(PT_COORD2(:)%LAT)+COS(REF_COORD%LON)* &
   & COS(PT_COORD2(:)%LAT)*SIN(PT_COORD2(:)%LON))/COS(PT_COORD1(:)%LAT))
ELSEWHERE
  PT_COORD1(:)%LON = 0.0_JPRB
ENDWHERE

IF (LHOOK) CALL DR_HOOK('EGGMRT:TILT_V',1,ZHOOK_HANDLE)
END FUNCTION TILT_V

! -------------------------------------------------------------------------------
! Function to ROTATE Scalar
! -------------------------------------------------------------------------------
TYPE(LOLA) FUNCTION ROTATE_S(CENTER_COORD,PT_COORD1) RESULT (PT_COORD)
TYPE (LOLA), INTENT(IN)                           :: CENTER_COORD, PT_COORD1

REAL(KIND=JPHOOK)                         :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('EGGMRT:ROTATE_S',0,ZHOOK_HANDLE)
PT_COORD%LAT = P_ASIN(COS(CENTER_COORD%LAT)*SIN(PT_COORD1%LAT)+SIN(CENTER_COORD%LAT)* &
 & COS(PT_COORD1%LAT)*COS(PT_COORD1%LON))

IF (COS(PT_COORD%LAT) /= 0.0_JPRB) THEN
  PT_COORD%LON = COSIN_TO_ANGLE((-SIN(CENTER_COORD%LAT)*SIN(PT_COORD1%LAT)+ &
   & COS(CENTER_COORD%LAT)*COS(PT_COORD1%LAT)*COS(PT_COORD1%LON))/COS(PT_COORD%LAT), &
   & (COS(PT_COORD1%LAT)*SIN(PT_COORD1%LON))/COS(PT_COORD%LAT))+CENTER_COORD%LON
ELSE
  PT_COORD%LON = 0.0_JPRB
ENDIF
IF (LHOOK) CALL DR_HOOK('EGGMRT:ROTATE_S',1,ZHOOK_HANDLE)
END FUNCTION ROTATE_S
! -------------------------------------------------------------------------------
! Function to ROTATE Vector
! -------------------------------------------------------------------------------
FUNCTION ROTATE_V(CENTER_COORD,PT_COORD1) RESULT (PT_COORD)
TYPE (LOLA), INTENT(IN)                               :: CENTER_COORD
TYPE (LOLA), DIMENSION(:),INTENT(IN)                  :: PT_COORD1
TYPE (LOLA), DIMENSION(SIZE(PT_COORD1)) :: PT_COORD

REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('EGGMRT:ROTATE_V',0,ZHOOK_HANDLE)
PT_COORD(:)%LAT = P_ASIN(COS(CENTER_COORD%LAT)*SIN(PT_COORD1(:)%LAT)+SIN(CENTER_COORD%LAT)* &
 & COS(PT_COORD1(:)%LAT)*COS(PT_COORD1(:)%LON))

WHERE (COS(PT_COORD(:)%LAT) /= 0.0_JPRB)
  PT_COORD(:)%LON = COSIN_TO_ANGLE((-SIN(CENTER_COORD%LAT)*SIN(PT_COORD1(:)%LAT)+ &
   & COS(CENTER_COORD%LAT)*COS(PT_COORD1(:)%LAT)*COS(PT_COORD1(:)%LON))/COS(PT_COORD(:)%LAT), &
   & (COS(PT_COORD1(:)%LAT)*SIN(PT_COORD1(:)%LON))/COS(PT_COORD(:)%LAT))+CENTER_COORD%LON
ELSEWHERE
  PT_COORD(:)%LON = 0.0_JPRB
ENDWHERE
IF (LHOOK) CALL DR_HOOK('EGGMRT:ROTATE_V',1,ZHOOK_HANDLE)
END FUNCTION ROTATE_V
! -------------------------------------------------------------------------------
! Function to ROTATE & TILT Scalar
! -------------------------------------------------------------------------------
TYPE(LOLA) FUNCTION MEROTIL_S(REF_COORD,CENTER_COORD,PT_COORD2) RESULT (PT_COORD)
TYPE (LOLA), INTENT(IN)                  :: REF_COORD, CENTER_COORD, PT_COORD2

REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('EGGMRT:MEROTIL_S',0,ZHOOK_HANDLE)
PT_COORD=ROTATE(CENTER_COORD,TILT(REF_COORD,PT_COORD2))
IF (LHOOK) CALL DR_HOOK('EGGMRT:MEROTIL_S',1,ZHOOK_HANDLE)
END FUNCTION MEROTIL_S
! -------------------------------------------------------------------------------
! Function to ROTATE & TILT Vector
! -------------------------------------------------------------------------------
FUNCTION MEROTIL_V(REF_COORD,CENTER_COORD,PT_COORD2) RESULT (PT_COORD)
TYPE (LOLA), INTENT(IN)                              :: REF_COORD, CENTER_COORD
TYPE (LOLA), DIMENSION(:),INTENT(IN)                 :: PT_COORD2
TYPE (LOLA), DIMENSION(SIZE(PT_COORD2)) :: PT_COORD

REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('EGGMRT:MEROTIL_V',0,ZHOOK_HANDLE)
PT_COORD(:)=ROTATE(CENTER_COORD,TILT(REF_COORD,PT_COORD2(:)))
IF (LHOOK) CALL DR_HOOK('EGGMRT:MEROTIL_V',1,ZHOOK_HANDLE)
END FUNCTION MEROTIL_V
! -------------------------------------------------------------------------------
! -------------------------------------------------------------------------------
! Function to ANTI_TILT Scalar
! -------------------------------------------------------------------------------
TYPE(LOLA) FUNCTION ANTI_TILT_S(REF_COORD,PT_COORD1) RESULT (PT_COORD2)
TYPE (LOLA), INTENT(IN)                       :: REF_COORD, PT_COORD1

TYPE (LOLA)        :: ANTI_REF_COORD
REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('EGGMRT:ANTI_TILT_S',0,ZHOOK_HANDLE)
ANTI_REF_COORD%LAT = REF_COORD%LAT
ANTI_REF_COORD%LON = -(1.0_JPRB)*REF_COORD%LON
PT_COORD2 = TILT(ANTI_REF_COORD,PT_COORD1)
IF (LHOOK) CALL DR_HOOK('EGGMRT:ANTI_TILT_S',1,ZHOOK_HANDLE)
END FUNCTION ANTI_TILT_S
! -------------------------------------------------------------------------------
! Function to ANTI_TILT Vector
! -------------------------------------------------------------------------------
FUNCTION ANTI_TILT_V(REF_COORD,PT_COORD1) RESULT (PT_COORD2)
TYPE (LOLA), INTENT(IN)                                     :: REF_COORD
TYPE (LOLA), DIMENSION(:),INTENT(IN)                        :: PT_COORD1
TYPE (LOLA), DIMENSION(SIZE(PT_COORD1))        :: PT_COORD2

TYPE (LOLA)        :: ANTI_REF_COORD
REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('EGGMRT:ANTI_TILT_V',0,ZHOOK_HANDLE)
ANTI_REF_COORD%LAT = REF_COORD%LAT
ANTI_REF_COORD%LON = -(1.0_JPRB)*REF_COORD%LON
PT_COORD2(:) = TILT(ANTI_REF_COORD,PT_COORD1(:))
IF (LHOOK) CALL DR_HOOK('EGGMRT:ANTI_TILT_V',1,ZHOOK_HANDLE)
END FUNCTION ANTI_TILT_V
! -------------------------------------------------------------------------------
! Function to ANTI_ROTATE Scalar
! -------------------------------------------------------------------------------
TYPE(LOLA) FUNCTION ANTI_ROTATE_S(CENTER_COORD,PT_COORD) RESULT (PT_COORD1)
TYPE (LOLA), INTENT(IN)                       :: CENTER_COORD, PT_COORD

TYPE (LOLA)        :: ANTI_CENTER_COORD, ANTI_PT_COORD
REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('EGGMRT:ANTI_ROTATE_S',0,ZHOOK_HANDLE)
ANTI_CENTER_COORD%LAT = -(1.0_JPRB)*CENTER_COORD%LAT
ANTI_CENTER_COORD%LON = CENTER_COORD%LON
ANTI_PT_COORD%LAT = PT_COORD%LAT
ANTI_PT_COORD%LON = PT_COORD%LON-CENTER_COORD%LON
PT_COORD1 = ROTATE(ANTI_CENTER_COORD,ANTI_PT_COORD)
PT_COORD1%LON = PT_COORD1%LON-CENTER_COORD%LON
IF (LHOOK) CALL DR_HOOK('EGGMRT:ANTI_ROTATE_S',1,ZHOOK_HANDLE)
END FUNCTION ANTI_ROTATE_S
! -------------------------------------------------------------------------------
! Function to ANTI_ROTATE Vector
! -------------------------------------------------------------------------------
FUNCTION ANTI_ROTATE_V(CENTER_COORD,PT_COORD) RESULT (PT_COORD1)
TYPE (LOLA), INTENT(IN)                                    :: CENTER_COORD
TYPE (LOLA), DIMENSION(:),INTENT(IN)                       :: PT_COORD
TYPE (LOLA), DIMENSION(SIZE(PT_COORD))        :: PT_COORD1

TYPE (LOLA)                            :: ANTI_CENTER_COORD
TYPE (LOLA), DIMENSION(SIZE(PT_COORD)) :: ANTI_PT_COORD
REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('EGGMRT:ANTI_ROTATE_V',0,ZHOOK_HANDLE)
ANTI_CENTER_COORD%LAT = -(1.0_JPRB)*CENTER_COORD%LAT
ANTI_CENTER_COORD%LON = CENTER_COORD%LON
ANTI_PT_COORD(:)%LAT = PT_COORD(:)%LAT
ANTI_PT_COORD(:)%LON = PT_COORD(:)%LON-CENTER_COORD%LON
PT_COORD1(:) = ROTATE(ANTI_CENTER_COORD,ANTI_PT_COORD(:))
PT_COORD1(:)%LON = PT_COORD1(:)%LON-CENTER_COORD%LON
IF (LHOOK) CALL DR_HOOK('EGGMRT:ANTI_ROTATE_V',1,ZHOOK_HANDLE)
END FUNCTION ANTI_ROTATE_V
! -------------------------------------------------------------------------------
! Function to ANTI_ROTATE & ANTI_TILT Scalar
! -------------------------------------------------------------------------------
TYPE(LOLA) FUNCTION METILROT_S(REF_COORD,CENTER_COORD,PT_COORD) RESULT (PT_COORD2)
TYPE (LOLA), INTENT(IN)                   :: REF_COORD, CENTER_COORD, PT_COORD

REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('EGGMRT:METILROT_S',0,ZHOOK_HANDLE)
PT_COORD2=ANTI_TILT(REF_COORD,ANTI_ROTATE(CENTER_COORD,PT_COORD))
IF (LHOOK) CALL DR_HOOK('EGGMRT:METILROT_S',1,ZHOOK_HANDLE)
END FUNCTION METILROT_S
! -------------------------------------------------------------------------------
! Function to ANTI_ROTATE & ANTI_TILT Vector
! -------------------------------------------------------------------------------
FUNCTION METILROT_V(REF_COORD,CENTER_COORD,PT_COORD) RESULT (PT_COORD2)
TYPE (LOLA), INTENT(IN)                             :: REF_COORD, CENTER_COORD
TYPE (LOLA), DIMENSION(:),INTENT(IN)                :: PT_COORD
TYPE (LOLA), DIMENSION(SIZE(PT_COORD)) :: PT_COORD2

REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('EGGMRT:METILROT_V',0,ZHOOK_HANDLE)
PT_COORD2(:)=ANTI_TILT(REF_COORD,ANTI_ROTATE(CENTER_COORD,PT_COORD(:)))
IF (LHOOK) CALL DR_HOOK('EGGMRT:METILROT_V',1,ZHOOK_HANDLE)
END FUNCTION METILROT_V
! -------------------------------------------------------------------------------
END MODULE EGGMRT
