! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

!OPTIONS XOPT(HSFUN)
SUBROUTINE VDFDPBL(YDEPHY,KIDIA,KFDIA,KLON,KLEV,&
 & PUM1,PVM1,PTM1,PQM1,PGEOM1,&
 & PKMFL,PKHFL,PKQFL,PDHPBL)  
!     ------------------------------------------------------------------

!**   *VDFDPBL* - VDFDPBL (Diagnostic PBL height) determines  
!                 PBL height for diagnostic purposes only

!     A.C.M. BELJAARS       E.C.M.W.F.    17/02/1998.

!     PURPOSE
!     -------

!     Determine PBL height

!     INTERFACE
!     ---------

!     *VDFDPBL* IS CALLED BY *VDFMAIN*

!     INPUT PARAMETERS (INTEGER):

!     *KIDIA*        START POINT
!     *KFDIA*        END POINT
!     *KLEV*         NUMBER OF LEVELS
!     *KLON*         NUMBER OF GRID POINTS PER PACKET

!     INPUT PARAMETERS (REAL):

!     *PUM1*         X-VELOCITY COMPONENT AT T-1
!     *PVM1*         Y-VELOCITY COMPONENT AT T-1
!     *PTM1*         TEMPERATURE AT T-1
!     *PQM1*         SPECIFIC HUMUDITY AT T-1
!     *PGEOM1*       GEOPOTENTIAL AT T-1
!     *PKMFL*        KINEMATIC MOMENTUM FLUX                
!     *PKHFL*        KINEMATIC HEAT FLUX                    
!     *PKQFL*        KINEMATIC MOISTURE FLUX           

!     OUTPUT PARAMETERS (REAL):

!     *PDHPBL*        Boundary layer height                  m

!     METHOD
!     ------

!     Troen and Mahrt method using bulk Richardson criterion
!     (see documentation)
!
!     MODIFICATIONS.
!     --------------
!      N. Semane+P.Bechtold     04-10-2012 Add RPLRG/RPLDARE factor for small planet
!     I.Sandu 15/03/2013 new algorithm based on Seidel et al. 2012
!     ------------------------------------------------------------------

USE PARKIND1  , ONLY : JPIM, JPRB
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST    , ONLY : RG, RCPD, RETV
USE PARPHY    , ONLY : REPDU2
USE YOEPHY    , ONLY : TEPHY
USE YOMDYNCORE, ONLY : RPLRG

IMPLICIT NONE

TYPE(TEPHY)       ,INTENT(INOUT) :: YDEPHY
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKMFL(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKHFL(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKQFL(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDHPBL(KLON) 

!*            LOCAL STORAGE
!             ----- -------

LOGICAL :: LLDONE(KLON)
REAL(KIND=JPRB) :: ZRI(KLON),ZDU2(KLON),&
 & ZSVBOT(KLON),ZSVBOTP(KLON)


INTEGER(KIND=JPIM) :: ILEVM1, ITOT, JIT, JK, JL

REAL(KIND=JPRB) :: ZBUST, ZCONS13, ZCONS14, ZCONS15, ZDRORO,&
 & ZEPS, ZPAR, ZPAR1, ZPARZI, ZRICRI, ZRILEV, ZSV, ZREPUST
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "surf_inq.h"

!     ------------------------------------------------------------------

!*       1.     INITIALIZE CONSTANTS
!               ---------- ----------

IF (LHOOK) CALL DR_HOOK('VDFDPBL',0,ZHOOK_HANDLE)

ZPAR   = 8.5_JPRB
ZPAR1  = 0.6_JPRB
ZCONS14= RCPD*ZPAR
ZCONS15= ZPAR1*RG
ZEPS   = 1.E-10_JPRB
ZBUST  = 100._JPRB
ZRICRI = 0.25_JPRB
ZPARZI = 1000._JPRB/RPLRG

CALL SURF_INQ(YDEPHY%YSURF,PREPUST=ZREPUST)

ZCONS13=1.0_JPRB/3._JPRB

ILEVM1 = KLEV-1

!        1.    PREPARE SURFACE PARAMETERS
!              ------- ------- ----------

DO JL=KIDIA,KFDIA
  ZSVBOT(JL)=RCPD*PTM1(JL,KLEV)*(1.0_JPRB+RETV*PQM1(JL,KLEV))+PGEOM1(JL,KLEV)
  PDHPBL(JL)=ZPARZI
ENDDO

!        2.    DO 1 OR 2 ITERATIONS ON PBL HEIGHT
!              -- - -- - ---------- -- --- ------

!*************************************
JIT=1
!*************************************

!        2.1   UPDATE VELOCITY SCALE AND EXCESS TEMPERATURE
!              ------ -------- ----- --- ------ -----------

  DO JL=KIDIA,KFDIA
    LLDONE(JL)=.FALSE.
    ZRI(JL)=0.0_JPRB
    ZSVBOTP(JL)=ZSVBOT(JL)+ZEPS
  ENDDO

!        2.2    VERTICAL SCAN TO DETERMINE MIXED LAYER DEPTH
!               -------- ---- -- --------- ----- ----- -----

!***
  DO JK=ILEVM1,1,-1
!***

    ITOT=KFDIA-KIDIA+1
    DO JL=KIDIA,KFDIA
      IF (.NOT. LLDONE(JL)) THEN
        ZSV=RCPD*PTM1(JL,JK)*(1.0_JPRB+RETV*PQM1(JL,JK))+PGEOM1(JL,JK)
        ZDU2(JL)=MAX(REPDU2, PUM1(JL,JK)**2+PVM1(JL,JK)**2) 
        ZDRORO=(ZSV-ZSVBOTP(JL)) &
         & /(ZSV-PGEOM1(JL,JK)) 
!
        ZRILEV=(PGEOM1(JL,JK)-PGEOM1(JL,KLEV))*ZDRORO/ZDU2(JL)
        IF (ZRILEV  >  ZRICRI) THEN
          PDHPBL(JL)=( (ZRILEV-ZRICRI)*PGEOM1(JL,JK+1)&
           & +(ZRICRI-ZRI(JL))*PGEOM1(JL,JK) )/&
           & ((ZRILEV-ZRI(JL))*RG)  
          LLDONE(JL)=.TRUE.
          ITOT=ITOT-1
        ELSE
          ZRI(JL)=ZRILEV
        ENDIF
      ELSE
        ITOT=ITOT-1
      ENDIF
    ENDDO
    IF (ITOT  <=  0) EXIT
!***
  ENDDO
!***

!*************************************
!*************************************

IF (LHOOK) CALL DR_HOOK('VDFDPBL',1,ZHOOK_HANDLE)
END SUBROUTINE VDFDPBL
