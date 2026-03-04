! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE VDFVINT(KIDIA,KFDIA,KLON,KLEV,&
 & PUM1  ,PVM1  ,PGEOM1, PZZ,&
 ! OUTPUTS
 & PUZ, PVZ )

!     ------------------------------------------------------------------

!**   *VDFVINT* - COMPUTES THE WIND SPEED AT A SPECIFIED HEIGHT

!     A. Beljaars        E.C.M.W.F.    19/09/2009. (BASED ON VDFFBLEND)
!     N. Semane+P.Bechtold     04-10-2012 Add RPLRG factor for small planet

!     PURPOSE
!     -------

!     COMPUTE WIND SPEED AT SPECIFIED HEIGHT (100m) for post-processing
!

!     INTERFACE
!     ---------

!     *VDFFBLEND* IS CALLED BY *VDFMAIN*

!     INPUT PARAMETERS (INTEGER):

!     *KIDIA*        START POINT
!     *KFDIA*        END POINT
!     *KLEV*         NUMBER OF LEVELS
!     *KLON*         NUMBER OF GRID POINTS PER PACKET

!     INPUT PARAMETERS (REAL):

!     *PUM1*         U-COMPONENT WIND AT T-1
!     *PVM1*         V-COMPONENT WIND AT T-1
!     *PGEOM1*       GEOPOTENTIAL AT T-1
!     *PZZ*          HEIGHT (m)

!     OUTPUT PARAMETERS (REAL):

!     *PUZ*          U-COMPONENT OF WIND AT HEIGHT PZZ
!     *PVZ*          V-COMPONENT OF WIND AT HEIGHT PZZ

!     METHOD
!     ------

!     LINEAR INTERPOLATION IN HEIGHT FROM THE MODEL LEVEL WIND SPEED TO THE
!     REQUESTED HEIGHT

!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST   , ONLY : RG
USE YOMDYNCORE,ONLY : RPLRG

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PZZ 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PUZ(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVZ(KLON) 

INTEGER(KIND=JPIM) :: JK, JL

REAL(KIND=JPRB) :: ZZR
 
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('VDFVINT',0,ZHOOK_HANDLE)

ZZR=PZZ*(RG/RPLRG)

DO JL=KIDIA,KFDIA
  DO JK=KLEV,2,-1
    IF (ZZR  <  PGEOM1(JL,JK-1) .AND.&
       & ZZR  >=  PGEOM1(JL,JK)) THEN  
      PUZ(JL)=( PUM1(JL,JK-1)*(ZZR-PGEOM1(JL,JK))&
        & +PUM1(JL,JK)*(PGEOM1(JL,JK-1)-ZZR)&
        & )/(PGEOM1(JL,JK-1)-PGEOM1(JL,JK))  
      PVZ(JL)=( PVM1(JL,JK-1)*(ZZR-PGEOM1(JL,JK))&
        & +PVM1(JL,JK)*(PGEOM1(JL,JK-1)-ZZR)&
        & )/(PGEOM1(JL,JK-1)-PGEOM1(JL,JK))  
      EXIT 
    ENDIF
  ENDDO
ENDDO

IF (LHOOK) CALL DR_HOOK('VDFVINT',1,ZHOOK_HANDLE)
END SUBROUTINE VDFVINT
