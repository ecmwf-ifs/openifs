MODULE SRFENE_MOD
CONTAINS
SUBROUTINE SRFENE(&
 & KIDIA  , KFDIA  , KLON  , KLEVS,&
 & LDLAND , LDSICE ,&
 & PTSAM1M, KSOTY, PCVL , PCVH ,&
 & YDCST  , YDSOIL ,&
 & PENES)  

! (C) Copyright 1996- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!**** *SRFENE* - COMPUTES SOIL ENERGY FOR EACH LAYER.

!     Original  P.VITERBO      E.C.M.W.F.     26/03/96
!     Modified  P.VITERBO  99-03-26   Tiling of the land surface
!               P.VITERBO  2004-05-24 Move to surf library
!               G.BALSAMO  2006-07-03 Add soil type 

!     PURPOSE.
!     --------

!          THIS ROUTINE COMPUTES THE SOIL ENERGY
!          IN THE SOIL, AVOIDING DELICATE SNOW SITUATIONS. APPARENT
!          STANDS FOR THE FACT THAT THE EFFECTS OF FREEZING AND MELTING
!          OF WATER IN THE SOIL ARE TAKEN INTO ACCOUNT.

!**   INTERFACE.
!     ----------

!          *SRFENE* IS CALLED FROM DIAGNOSTIC (DDH) ROUTINES.

!     PARAMETER   DESCRIPTION                                    UNITS
!     ---------   -----------                                    -----

!     INPUT PARAMETERS (INTEGER):

!    *KIDIA*      START POINT
!    *KFDIA*      END POINT
!    *KLON*       NUMBER OF GRID POINTS PER PACKET
!    *KLEVS*      NUMBER OF SOIL LAYERS
!    *KSOTY*      SOIL TYPE                                        (1-7)

!     INPUT PARAMETERS (LOGICAL):

!    *LDLAND*     LAND/SEA MASK (TRUE/FALSE)
!    *LDSICE*     SEA ICE MASK (.T. OVER SEA ICE)

!     INPUT PARAMETERS AT T-1 OR CONSTANT IN TIME (REAL):

!    *PTSAM1M*    SOIL TEMPERATURE                                  K
!    *PCVL*       LOW VEGETATION COVER  (CORRECTED)                (0-1)
!    *PCVH*       HIGH VEGETATION COVER (CORRECTED)                (0-1)

!     OUTPUT PARAMETERS:

!    *PENES*      SOIL ENERGY per unit area                        J/M**2

!     METHOD.
!     -------

!          STRAIGHTFORWARD ONCE THE DEFINITION OF THE CONSTANTS IS
!     UNDERSTOOD. FOR THIS REFER TO DOCUMENTATION.

!     EXTERNALS.
!     ----------

!          NONE.

!     REFERENCE.
!     ----------

!     ------------------------------------------------------------------

USE PARKIND1  , ONLY : JPIM, JPRB
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOS_THF   , ONLY : RHOH2O
USE YOS_CST   , ONLY : TCST
USE YOS_SOIL  , ONLY : TSOIL

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEVS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSOTY(:) 
LOGICAL           ,INTENT(IN)    :: LDLAND(:) 
LOGICAL           ,INTENT(IN)    :: LDSICE(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSAM1M(:,:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCVL(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCVH(:)
TYPE(TCST)        ,INTENT(IN)    :: YDCST
TYPE(TSOIL)       ,INTENT(IN)    :: YDSOIL
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PENES(:,:) 
!*         0.1    DECLARATION OF GLOBAL VARIABLES.
!                 ----------- -- ------ ----------

REAL(KIND=JPRB) :: ZF(KLON,KLEVS)

INTEGER(KIND=JPIM) :: JK, JL, JS

REAL(KIND=JPRB) :: ZWA,ZRCSOIL,ZWCAP
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

!*         1.    SET UP SOME CONSTANTS.
!                --- -- ---- ----------

IF (LHOOK) CALL DR_HOOK('SRFENE_MOD:SRFENE',0,ZHOOK_HANDLE)
ASSOCIATE(RLMLT=>YDCST%RLMLT, &
 & LEVGEN=>YDSOIL%LEVGEN, RDAI=>YDSOIL%RDAI, RDAT=>YDSOIL%RDAT, &
 & RRCSICE=>YDSOIL%RRCSICE, RRCSOIL=>YDSOIL%RRCSOIL, RRCSOILM=>YDSOIL%RRCSOILM, &
 & RTF1=>YDSOIL%RTF1, RTF2=>YDSOIL%RTF2, RTF3=>YDSOIL%RTF3, RTF4=>YDSOIL%RTF4, &
 & RWCAP=>YDSOIL%RWCAP, RWCAPM=>YDSOIL%RWCAPM)

!     ------------------------------------------------------------------

!*         2. CONTRIBUTION TO APPARENT ENERGY.
!             --------------------------------

!          CONTRIBUTION TO APPARENT ENERGY, TAKING INTO ACCOUNT
!          FREEZING/MELTING OF SOIL WATER.

DO JK=1,KLEVS
  DO JL=KIDIA,KFDIA
    IF (LDLAND(JL)) THEN

!     NOTE: FUNCTION F(T). THE FUNCTION DFDT IN ROUTINE SRFRCG HAS
!           TO BE D/DT (F(T))

      IF(PTSAM1M(JL,JK) < RTF1.AND.PTSAM1M(JL,JK) > RTF2) THEN
        ZF(JL,JK)=0.5_JPRB*(1.0_JPRB-SIN(RTF4*(PTSAM1M(JL,JK)-RTF3)))
      ELSEIF (PTSAM1M(JL,JK) <= RTF2) THEN
        ZF(JL,JK)=1.0_JPRB
      ELSE
        ZF(JL,JK)=0.0_JPRB
      ENDIF
    ENDIF
  ENDDO
ENDDO

!     ------------------------------------------------------------------

!*         3. COMPUTE ENERGY.
!             ---------------

!          APPARENT SOIL ENERGY, TAKING INTO ACCOUNT FREEZING/
!          MELTING OF SOIL WATER. NOTHING IS DONE TO THE FIRST LAYER
!          IN THE PRESENCE OF A SIGNIFICANT AMOUNT OF SNOW, BECAUSE OF
!          THE MIXED THERMAL NATURE OF THE SOIL. 

DO JK=1,KLEVS
  DO JL=KIDIA,KFDIA

!          SOIL THERMAL COEFFICIENTS MODIFIED WHEN SNOW COVERS
!          THE GROUND AND IS PARTIALLY MASKED BY THE VEGETATION.

    IF (LDLAND(JL)) THEN
      IF (LEVGEN) THEN
         JS=KSOTY(JL)
         ZWCAP=RWCAPM(JS)
         ZRCSOIL=RRCSOILM(JS)
      ELSE
         ZWCAP=RWCAP
         ZRCSOIL=RRCSOIL
      ENDIF
      ZWA=(PCVL(JL)+PCVH(JL))*ZWCAP
      PENES(JL,JK)=ZRCSOIL*PTSAM1M(JL,JK)
      PENES(JL,JK)=(PENES(JL,JK)-RLMLT*RHOH2O*ZWA*ZF(JL,JK))*RDAT(JK)

!          SEA ICE POINTS

    ELSEIF (LDSICE(JL)) THEN
      PENES(JL,JK)=RRCSICE*PTSAM1M(JL,JK)*RDAI(JK)

!          SEA POINTS

    ELSE
      PENES(JL,JK)=0.0_JPRB
    ENDIF
  ENDDO
ENDDO

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SRFENE_MOD:SRFENE',1,ZHOOK_HANDLE)
END SUBROUTINE SRFENE
END MODULE SRFENE_MOD
