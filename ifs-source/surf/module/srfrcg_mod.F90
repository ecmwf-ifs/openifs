MODULE SRFRCG_MOD
CONTAINS
SUBROUTINE SRFRCG(KIDIA  , KFDIA  , KLON , KTILES, KLEVS ,&
 & LDLAND , LDSICE ,&
 & PTSAM1M, KSOTY  , PCVL   , PCVH  , PCUR,&
 & YDCST, YDSOIL   , YDURB  ,&
 & PCTSA)  
 
USE PARKIND1  , ONLY : JPIM, JPRB
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOS_THF   , ONLY : RHOH2O
USE YOS_CST   , ONLY : TCST
USE YOS_SOIL  , ONLY : TSOIL
USE YOS_URB   , ONLY : TURB

! (C) Copyright 1993- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!**** *SRFT* - COMPUTES SOIL VOLUMETRIC HEAT CAPACITY.
!     PURPOSE.
!     --------
!          THIS ROUTINE COMPUTES THE APARENT VOLUMETRIC HEAT CAPACITY
!          IN THE SOIL, TAKING INTO ACCOUNT SNOW. APPARENT STANDS FOR
!          THE FACT THAT THE EFFECTS OF FREEZING AND MELTING OF WATER
!          IN THE SOIL ARE TAKEN INTO ACCOUNT.

!**   INTERFACE.
!     ----------
!          *SRFRCG* IS CALLED FROM *SRFT* AND DIAGNOSTIC (DDH) ROUTINES.

!     PARAMETER   DESCRIPTION                                    UNITS
!     ---------   -----------                                    -----
!     INPUT PARAMETERS (INTEGER):
!    *KIDIA*      START POINT
!    *KFDIA*      END POINT
!    *KLON*       NUMBER OF GRID POINTS PER PACKET
!    *KTILES*     NUMBER OF TILES
!    *KLEVS*      NUMBER OF SOIL LAYERS

!     INPUT PARAMETERS (LOGICAL):
!    *LDLAND*     LAND/SEA MASK (TRUE/FALSE)
!    *LDSICE*     SEA ICE MASK (.T. OVER SEA ICE)

!     INPUT PARAMETERS AT T-1 OR CONSTANT IN TIME (REAL):
!    *PTSAM1M*    SOIL TEMPERATURE                               K
!    *PCVL*       LOW VEGETATION COVER  (CORRECTED)              (0-1)
!    *PCVH*       HIGH VEGETATION COVER (CORRECTED)              (0-1)
!    *PCUR*       URBAN COVER                                    (0-1)

!     OUTPUT PARAMETERS:
!    *PCTSA*      VOLUMETRIC HEAT CAPACITY                      J/K/M**3

!     METHOD.
!     -------
!          STRAIGHTFORWARD ONCE THE DEFINITION OF THE CONSTANTS IS
!     UNDERSTOOD. FOR THIS REFER TO DOCUMENTATION.

!     EXTERNALS.
!     ----------
!          NONE.

!     REFERENCE.
!     ----------
!     Original    P.VITERBO      E.C.M.W.F.     14/10/93
!     Modified    P.VITERBO  99-03-26   Tiling of the land surface
!     Modified    J.F. Estrade *ECMWF* 03-10-01 move in surf vob
!     Modified    P. Viterbo   24/05/2004   Change surface units
!     Modified    G. Balsamo   10/01/2006   Include Van Genuchten Hydro.
!     Modified    G. Balsamo   03/07/2006   Add soil type
!     ------------------------------------------------------------------

IMPLICIT NONE

! Declaration of arguments

INTEGER(KIND=JPIM), INTENT(IN)   :: KIDIA
INTEGER(KIND=JPIM), INTENT(IN)   :: KFDIA
INTEGER(KIND=JPIM), INTENT(IN)   :: KLON
INTEGER(KIND=JPIM), INTENT(IN)   :: KTILES
INTEGER(KIND=JPIM), INTENT(IN)   :: KLEVS
INTEGER(KIND=JPIM), INTENT(IN)   :: KSOTY(:)

REAL(KIND=JPRB),    INTENT(IN)   :: PTSAM1M(:,:)
REAL(KIND=JPRB),    INTENT(IN)   :: PCVL(:)
REAL(KIND=JPRB),    INTENT(IN)   :: PCVH(:)
REAL(KIND=JPRB),    INTENT(IN)   :: PCUR(:)

LOGICAL,   INTENT(IN)   :: LDLAND(:)
LOGICAL,   INTENT(IN)   :: LDSICE(:)

TYPE(TCST),         INTENT(IN)   :: YDCST
TYPE(TSOIL),        INTENT(IN)   :: YDSOIL
TYPE(TURB),         INTENT(IN)   :: YDURB

REAL(KIND=JPRB),    INTENT(OUT)  :: PCTSA(:,:)

INTEGER(KIND=JPIM) :: JK, JL, JS

REAL(KIND=JPRB) :: ZGICE, ZRCSICE, ZSNOWI, ZWA
REAL(KIND=JPRB) :: ZRCSOIL, ZWCAP
REAL(KIND=JPRB) :: ZDFDT(KLON,KLEVS)
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
!*         1.    SET UP SOME CONSTANTS.
!                --- -- ---- ----------
!*               PHYSICAL CONSTANTS.
!                -------- ----------

IF (LHOOK) CALL DR_HOOK('SRFRCG_MOD:SRFRCG',0,ZHOOK_HANDLE)
ASSOCIATE(RLMLT=>YDCST%RLMLT, &
 & LEVGEN=>YDSOIL%LEVGEN, RDAT=>YDSOIL%RDAT, RGH2O=>YDSOIL%RGH2O, &
 & RRCSICE=>YDSOIL%RRCSICE, RRCSOIL=>YDSOIL%RRCSOIL, RRCSOILM=>YDSOIL%RRCSOILM, &
 & RTF1=>YDSOIL%RTF1, RTF2=>YDSOIL%RTF2, RTF3=>YDSOIL%RTF3, RTF4=>YDSOIL%RTF4, &
 & RWCAP=>YDSOIL%RWCAP, RWCAPM=>YDSOIL%RWCAPM, RURBVHC=>YDURB%RURBVHC)

ZGICE=0.5_JPRB*RGH2O
ZSNOWI=1.0_JPRB/RDAT(1)
ZRCSICE=RRCSICE

!     ------------------------------------------------------------------
!*         2. CONTRIBUTION TO APPARENT HEAT CAPACITY.
!             ---------------------------------------
!          CONTRIBUTION TO APPARENT HEAT CAPACITY, TAKING INTO ACCOUNT
!          FREEZING/MELTING OF SOIL WATER.

DO JK=1,KLEVS
  DO JL=KIDIA,KFDIA
    IF (LDLAND(JL)) THEN

!     NOTE: ZDFDT STANDS FOR D/DT OF THE FUNCTION F(T) IN THE
!           ROUTINE SRFENE

      IF(PTSAM1M(JL,JK) < RTF1.AND.PTSAM1M(JL,JK) > RTF2) THEN
        ZDFDT(JL,JK)=-0.5_JPRB*RTF4*COS(RTF4*(PTSAM1M(JL,JK)-RTF3))
      ELSE
        ZDFDT(JL,JK)=0.0_JPRB
      ENDIF
    ENDIF
  ENDDO
ENDDO

!     ------------------------------------------------------------------
!*         3. COMPUTE HEAT CAPACITIES.
!             ------------------------

DO JL=KIDIA,KFDIA

  IF (LDLAND(JL)) THEN

!          SOIL THERMAL COEFFICIENTS MODIFIED WHEN SNOW COVERS
!          THE GROUND AND IS PARTIALLY MASKED BY THE VEGETATION.
    IF(LEVGEN)THEN
      JS=KSOTY(JL)
      ZWCAP=RWCAPM(JS)
      ZRCSOIL=RRCSOILM(JS)
    ELSE
      ZWCAP=RWCAP
      ZRCSOIL=RRCSOIL
    ENDIF
    ZWA=(PCVL(JL)+PCVH(JL))*ZWCAP
    DO JK=1,KLEVS
      PCTSA(JL,JK)=ZRCSOIL-RLMLT*RHOH2O*ZWA*ZDFDT(JL,JK)
    ENDDO

!          URBAN TOP-LAYER
    IF (KTILES .GT. 9) THEN
     IF(LEVGEN)THEN
       ZRCSOIL=RRCSOILM(JS)*(1-PCUR(JL))+RURBVHC*PCUR(JL)
     ELSE
       ZRCSOIL=RRCSOIL*(1-PCUR(JL))+RURBVHC*PCUR(JL)
     ENDIF
     PCTSA(JL,1)=ZRCSOIL-RLMLT*RHOH2O*ZWA*ZDFDT(JL,1)
    ENDIF

!          SEA ICE POINTS

  ELSEIF (LDSICE(JL)) THEN
    PCTSA(JL,1)=ZRCSICE
    DO JK=2,KLEVS
      PCTSA(JL,JK)=0.0_JPRB
    ENDDO

!          SEA POINTS

  ELSE
    DO JK=1,KLEVS
      PCTSA(JL,JK)=0.0_JPRB
    ENDDO
  ENDIF
ENDDO
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SRFRCG_MOD:SRFRCG',1,ZHOOK_HANDLE)

END SUBROUTINE SRFRCG
END MODULE SRFRCG_MOD
