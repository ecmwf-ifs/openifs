MODULE SRFWL_MOD
CONTAINS
SUBROUTINE SRFWL(KIDIA,KFDIA,KTILES,&
 & PTMST,PWLM1M,PCVL,PCVH,PWLMX,&
 & PFRTI,PEVAPTI,PRSFC,PRSFL,PEVAPSNW,&
 & PWL,PFWEL1,PTSFC,PTSFL,PEINTTI,&
 & LDLAND,&
 & YDSOIL,YDVEG,&
 & PDHIIS)  
 
USE PARKIND1  , ONLY : JPIM, JPRB
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOS_THF   , ONLY : RHOH2O
USE YOS_SOIL  , ONLY : TSOIL
USE YOS_VEG   , ONLY : TVEG
     
! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!**** *SRFWL* - COMPUTES CHANGES IN THE SKIN RESERVOIR.
!     PURPOSE.
!     --------
!          THIS ROUTINE COMPUTES THE CHANGES IN THE SKIN RESERVOIR AND
!     THE RUN-OFF BEFORE THE SNOW MELTS.

!**   INTERFACE.
!     ----------
!          *SRFWL* IS CALLED FROM *SURF*.

!     PARAMETER   DESCRIPTION                                    UNITS
!     ---------   -----------                                    -----

!     INPUT PARAMETERS (INTEGER):
!    *KIDIA*      START POINT
!    *KFDIA*      END POINT
!    *KLON*       NUMBER OF GRID POINTS PER PACKET
!    *KTILES*     NUMBER OF TILES (I.E. SUBGRID AREAS WITH DIFFERENT 
!                 SURFACE BOUNDARY CONDITION)
!    *KDHVIIS*    Number of variables for interception layer water budget
!    *KDHFIIS*    Number of fluxes for interception layer water budget

!     INPUT PARAMETERS (REAL):
!    *PTMST*      TIME STEP                                      S

!     INPUT PARAMETERS (LOGICAL):
!    *LDLAND*     LAND/SEA MASK (TRUE/FALSE)

!     INPUT PARAMETERS AT T-1 OR CONSTANT IN TIME (REAL):
!    *PWLM1M*     SKIN RESERVOIR WATER CONTENT                   kg/m**2
!    *PCVL*       LOW VEGETATION COVER  (CORRECTED)              (0-1)
!    *PCVH*       HIGH VEGETATION COVER (CORRECTED)              (0-1)
!    *PWLMX*      MAXIMUM SKIN RESERVOIR CAPACITY                kg/m**2
!    *PFRTI*      TILE FRACTIONS                                 (0-1)
!            1 : WATER                  5 : SNOW ON LOW-VEG+BARE-SOIL
!            2 : ICE                    6 : DRY SNOW-FREE HIGH-VEG
!            3 : WET SKIN               7 : SNOW UNDER HIGH-VEG
!            4 : DRY SNOW-FREE LOW-VEG  8 : BARE SOIL
!            9 : LAKE                  10 : URBAN
!    *PEVAPTI*      SURFACE MOISTURE FLUX, FOR EACH TILE       KG/M2/S
!    *PRSFC*      CONVECTIVE RAIN FLUX AT THE SURFACE          KG/M**2/S
!    *PRSFL*      LARGE SCALE RAIN FLUX AT THE SURFACE         KG/M**2/S
!    *PEVAPSNW*   EVAPORATION FROM SNOW UNDER FOREST           KG/M**2/S

!     OUTPUT PARAMETERS AT T+1 (UNFILTERED,REAL):
!    *PWL*        SKIN RESERVOIR WATER CONTENT                   kg/m**2
!    *PFWEL1*     SKIN CONTRIBUTION TO THE EVAPORATION FROM
!                        SKIN AND TOP LAYER                    KG/M**2/S
!    *PTSFC*      CONVECTIVE THROUGHFALL AT THE SURFACE        KG/M**2/S
!    *PTSFL*      LARGE SCALE THROUGHFALL AT THE SURFACE       KG/M**2/S
!                  (NB: THROUGHFALL=RAINFALL-INTERCEPTION)
!    *PEINTTI*    TILE EVAPORATION SEEN BY THE INTERCEPTION
!                 LAYER (INCLUDES NUMERICAL EVAPORATION
!                 MISMATCHES, FOR TILE 3, AND DEW DEPOSITION
!                 FOR TILES 3,4,6,7,8)                         KG/M**2/S

!     OUTPUT PARAMETERS (DIAGNOSTIC):
!    *PDHIIS*     Diagnostic array for interception layer (see module yomcdh)

!     METHOD.
!     -------
!          STRAIGHTFORWARD ONCE THE DEFINITION OF THE CONSTANTS IS
!     UNDERSTOOD. FOR THIS REFER TO DOCUMENTATION.

!     EXTERNALS.
!     ----------
!          NONE.

!     REFERENCE.
!     ----------
!          SEE SOIL PROCESSES' PART OF THE MODEL'S DOCUMENTATION FOR
!     DETAILS ABOUT THE MATHEMATICS OF THIS ROUTINE.

!     Modifications:
!     Original   P.VITERBO      E.C.M.W.F.     16/01/89
!     Modified   P.VITERBO    99-03-26   Tiling of the land surface
!     Modified   P.VITERBO    17-05-2000 Surface DDH for TILES
!     Modified   J.F. Estrade *ECMWF* 03-10-01 move in surf vob
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!     P. Viterbo    24-05-2004      Change surface units
!     E. Dutra      07-07-2008      clean number of tiles dependence 
!     ------------------------------------------------------------------

IMPLICIT NONE

! Declaration of arguments

INTEGER(KIND=JPIM), INTENT(IN)   :: KIDIA
INTEGER(KIND=JPIM), INTENT(IN)   :: KFDIA
INTEGER(KIND=JPIM), INTENT(IN)   :: KTILES

REAL(KIND=JPRB),    INTENT(IN)   :: PTMST
REAL(KIND=JPRB),    INTENT(IN)   :: PWLM1M(:)
REAL(KIND=JPRB),    INTENT(IN)   :: PCVL(:)
REAL(KIND=JPRB),    INTENT(IN)   :: PCVH(:)
REAL(KIND=JPRB),    INTENT(IN)   :: PWLMX(:)
REAL(KIND=JPRB),    INTENT(IN)   :: PFRTI(:,:)
REAL(KIND=JPRB),    INTENT(IN)   :: PEVAPTI(:,:)
REAL(KIND=JPRB),    INTENT(IN)   :: PRSFC(:)
REAL(KIND=JPRB),    INTENT(IN)   :: PRSFL(:)
REAL(KIND=JPRB),    INTENT(IN)   :: PEVAPSNW(:)
LOGICAL,            INTENT(IN)   :: LDLAND(:)
TYPE(TSOIL),        INTENT(IN)   :: YDSOIL
TYPE(TVEG),         INTENT(IN)   :: YDVEG

REAL(KIND=JPRB),    INTENT(OUT)  :: PWL(:)
REAL(KIND=JPRB),    INTENT(OUT)  :: PFWEL1(:)
REAL(KIND=JPRB),    INTENT(OUT)  :: PTSFC(:)
REAL(KIND=JPRB),    INTENT(OUT)  :: PTSFL(:)
REAL(KIND=JPRB),    INTENT(OUT)  :: PEINTTI(:,:)
REAL(KIND=JPRB),    INTENT(OUT)  :: PDHIIS(:,:)

INTEGER(KIND=JPIM) :: JL, JTILE

REAL(KIND=JPRB) :: ZCONS1, ZELIQUID, ZEPFR, ZEPPRCP,&
 & ZEPTINY, ZIPRCP, ZMPRCP, ZPSFR, ZQHFLW, ZTMST, &
 & ZTPRCP, ZVINTER, ZWL
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
!*         1.    SET UP SOME CONSTANTS.
!                --- -- ---- ----------
!*               PHYSICAL CONSTANTS.
IF (LHOOK) CALL DR_HOOK('SRFWL_MOD:SRFWL',0,ZHOOK_HANDLE)
ASSOCIATE(RPSFR=>YDSOIL%RPSFR, &
 & RVINTER=>YDVEG%RVINTER)

ZVINTER=RVINTER
ZPSFR=1.0_JPRB/RPSFR

!*    SECURITY PARAMETERS
ZEPTINY=10._JPRB*TINY(RHOH2O)
ZEPPRCP=ZEPTINY
ZEPFR=10._JPRB*ZEPTINY

!*    COMPUTATIONAL CONSTANTS.
ZTMST=1.0_JPRB/PTMST
ZCONS1=ZTMST/RHOH2O

DO JTILE=1,KTILES
  DO JL=KIDIA,KFDIA
    PEINTTI(JL,JTILE)=0.0_JPRB
  ENDDO
ENDDO

!     ------------------------------------------------------------------
!*          2.  CORRECT FOR 0 < PWLM1M < WLMX.
!               ------- --- - - ------ - -----
DO JL=KIDIA,KFDIA
  IF (LDLAND(JL)) THEN
    PWL(JL)=MAX(0.0_JPRB,MIN(PWLMX(JL)-ZEPTINY,PWLM1M(JL)))
! Labelled evaporation seen by the interception reservoir
    PEINTTI(JL,3)=PEINTTI(JL,3)+(PWL(JL)-PWLM1M(JL))*ZTMST
  ELSE
!    Sea points
    PWL(JL)=0.0_JPRB
  ENDIF
ENDDO

!     ------------------------------------------------------------------
!*          3.  UPWARDS EVAPORATION.
!               ------- ------------
DO JL=KIDIA,KFDIA
  IF (LDLAND(JL)) THEN

!           INITIALISE PWL (TO MAKE THE CODE SIMPLER).
  IF ( .NOT. YDSOIL%LEWBSOILFIX ) THEN  
    PWL(JL)=PWLM1M(JL)
  ENDIF
    
!           EVAPORATION OF THE SKIN RESERVOIR (EL < 0).
    IF (PEVAPTI(JL,3) < 0.0_JPRB) THEN
      IF ( YDSOIL%LEWBSOILFIX ) THEN
        ZWL=PWL(JL)
      ELSE
        ZWL=PWLM1M(JL)
      ENDIF
      
      ZQHFLW=PTMST*PFRTI(JL,3)*PEVAPTI(JL,3)
      PWL(JL)=ZWL+ZQHFLW
      PWL(JL)=MAX(0.0_JPRB,PWL(JL))
! Evaporation seen by the interception reservoir
      PEINTTI(JL,3)=PEINTTI(JL,3)+(PWL(JL)-ZWL)*ZTMST
    ENDIF
  ENDIF
ENDDO

!           4.  COLLECTION OF DEW BY THE SKIN RESERVOIR.
!               ---------- -- --- -- --- ---- ----------
!              (Ci*Ei > 0.), for i=3,4,6,7,8
DO JTILE=1,KTILES
  IF ( JTILE /= 3 .AND. JTILE /= 4 .AND. JTILE /= 6 .AND. JTILE /= 7 .AND. JTILE /= 8 .AND. JTILE /= 10 ) CYCLE  
  DO JL=KIDIA,KFDIA
    IF (LDLAND(JL)) THEN
      IF (JTILE == 7) THEN
        ZELIQUID=PFRTI(JL,JTILE)*(PEVAPTI(JL,JTILE)-PEVAPSNW(JL))
      ELSE
        ZELIQUID=PFRTI(JL,JTILE)*PEVAPTI(JL,JTILE)
      ENDIF
      IF (ZELIQUID > 0.0_JPRB) THEN
        ZWL=PWL(JL)
        ZQHFLW=PTMST*ZELIQUID
        PWL(JL)=ZWL+MIN(PWLMX(JL)-ZEPTINY-ZWL,ZQHFLW)
        PEINTTI(JL,JTILE)=PEINTTI(JL,JTILE)+(PWL(JL)-ZWL)*ZTMST
      ENDIF
    ENDIF
  ENDDO
ENDDO

!           5.  BUDGETS.
!               --------
DO JL=KIDIA,KFDIA
  IF (LDLAND(JL)) THEN
    PFWEL1(JL)=PFWEL1(JL)+(PWL(JL)-PWLM1M(JL))*ZTMST
  ELSE
    PFWEL1(JL)=0.0_JPRB
  ENDIF
ENDDO

!           6.  INTERCEPTION OF PRECIPITATION BY THE VEGETATION.
!               ------------ -- ------------- -- --- -----------
DO JL=KIDIA,KFDIA
  IF (LDLAND(JL)) THEN
!          LARGE SCALE PRECIPITATION.
    IF (PRSFL(JL) > ZEPPRCP) THEN
      ZTPRCP=PRSFL(JL)
      ZIPRCP=ZTPRCP*(PCVL(JL)+PCVH(JL))*ZVINTER
      ZMPRCP=MIN(PWLMX(JL)-ZEPTINY-PWL(JL),PTMST*ZIPRCP)
      PWL(JL)=PWL(JL)+ZMPRCP
      PTSFL(JL)=PRSFL(JL)-ZMPRCP*ZTMST
    ELSE
      PTSFL(JL)=0.0_JPRB
    ENDIF

!          CONVECTIVE PRECIPITATION.
    IF (PRSFC(JL) > ZEPPRCP) THEN
      ZTPRCP=PRSFC(JL)*RPSFR
      ZIPRCP=ZTPRCP*(PCVL(JL)+PCVH(JL))*ZVINTER
      ZMPRCP=MIN(PWLMX(JL)-ZEPTINY-PWL(JL),PTMST*ZIPRCP)*ZPSFR
      PWL(JL)=PWL(JL)+ZMPRCP
      PTSFC(JL)=PRSFC(JL)-ZMPRCP*ZTMST
    ELSE
      PTSFC(JL)=0.0_JPRB
    ENDIF

!          SEA POINTS.
  ELSE
    PTSFC(JL)=PRSFC(JL)
    PTSFL(JL)=PRSFL(JL)
  ENDIF
ENDDO

!           7.  DDH Diagnostics.
!               ----------------
IF (SIZE(PDHIIS) > 0) THEN
! Interception layer water
  DO JL=KIDIA,KFDIA
    PDHIIS(JL,1)=PWLM1M(JL)
! Large-scale interception
    PDHIIS(JL,2)=PRSFL(JL)-PTSFL(JL)
! Convective interception
    PDHIIS(JL,3)=PRSFC(JL)-PTSFC(JL)
! Evaporation of intercepted water
    PDHIIS(JL,4)=PEINTTI(JL,3)+PEINTTI(JL,4)+&
      & PEINTTI(JL,6)+PEINTTI(JL,7)+&
      & PEINTTI(JL,8)  
  ENDDO
ENDIF

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SRFWL_MOD:SRFWL',1,ZHOOK_HANDLE)

END SUBROUTINE SRFWL
END MODULE SRFWL_MOD
