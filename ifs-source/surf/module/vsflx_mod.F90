MODULE VSFLX_MOD
CONTAINS
SUBROUTINE VSFLX(KIDIA,KFDIA,KLON,PTMST,PRVDIFTS,&
 & PUMLEV,PVMLEV,PTMLEV,PQMLEV,PAPHMS,&
 & PCPTGZLEV,PCPTS,PQS,PCFM,PCFH,PCFQ,PCAIR,PCSAT,&
 & PUCURR,PVCURR,&
 & YDCST,&
 & PKMFL,PKHFL,PKQFL)  

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOS_THF  , ONLY : RVTMP2
USE YOS_CST  , ONLY : TCST

! (C) Copyright 1991- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!     ------------------------------------------------------------------

!**   *VSFLX* - COMPUTES SURFACE FLUXES

!     Original   A.C.M. BELJAARS       E.C.M.W.F.    24-02-91
!     Modified   A.C.M. BELJAARS  26-03-99  Tiling of the land surface
!     Modified   P. Viterbo       15-05-2005  move to SURF library
!                                               (based on VDFSFLX)

!     PURPOSE
!     -------

!     COMPUTE SURFACE FLUXES FOR LATER USE IN SIMILARITY FUNCTIONS

!     INTERFACE
!     ---------

!     *VSFLX* IS CALLED BY *SURFEXCDRIVER*

!     INPUT PARAMETERS (INTEGER):

!     *KIDIA*        START POINT
!     *KFDIA*        END POINT
!     *KLON*         NUMBER OF POINTS IN PACKET

!     INPUT PARAMETERS (REAL):

!     *PTMST*        TIME STEP
!     *PRVDIFTS*     Semi-implicit factor for vertical diffusion discretization
!     *PUMLEV*       X-VELOCITY COMPONENT AT T-1, lowest model level
!     *PVMLEV*       Y-VELOCITY COMPONENT AT T-1, lowest model level
!     *PTMLEV*       TEMPERATURE AT T-1, lowest model level
!     *PQMLEV*       SPECIFIC HUMIDITY AT T-1, lowest model level
!     *PAPHMS*       PRESSURE AT T-1
!     *PCPTGZLEV*    DRY STATIC ENERGY AT T-1, lowest model level
!     *PCPTS*        DRY STATIC ENERGY AT SURFACE
!     *PQS*          SPECIFIC HUMIDITY AT SURFACE
!     *PCFM*         PROP. TO EXCH. COEFF. FOR MOMENTUM (C-STAR IN DOC.)
!                    (SURFACE LAYER ONLY)
!     *PCFH*         PROP. TO EXCH. COEFF. FOR HEAT     (C-STAR IN DOC.)
!                    (SURFACE LAYER ONLY)
!     *PCFQ*         PROP. TO EXCH. COEFF. FOR MOISTURE (C-STAR IN DOC.)
!                    (SURFACE LAYER ONLY)
!     *PCAIR*        MULTIPLICATION FACTOR FOR Q AT LOWEST MODEL LEVEL
!                    FOR SURFACE FLUX COMPUTATION
!     *PCSAT*        MULTIPLICATION FACTOR FOR QS AT SURFACE
!                    FOR SURFACE FLUX COMPUTATION
!     *PUCURR*       OCEAN CURRENT X-COMPONENT
!     *PVCURR*       OCEAN CURRENT Y-COMPONENT

!     OUTPUT PARAMETERS (REAL):

!     *PKMFL*        KINEMATIC MOMENTUM FLUX (DOWN=POS.)
!     *PKHFL*        KINEMATIC HEAT FLUX     (DOWN=POS.)
!     *PKQFL*        KINEMATIC MOISTURE FLUX (DOWN=POS.)

!     METHOD
!     ------

!     MULTIPLY DIFFERENCES BETWEEN SURFACE AND MODEL VARIABLES AT
!     THE LOWEST MODEL LEVEL WITH THE RESPECTIVE EXCHANGE COEFFICIENTS.

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTMST 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRVDIFTS
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUMLEV(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVMLEV(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTMLEV(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQMLEV(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHMS(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCPTGZLEV(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCPTS(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQS(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCFM(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCFH(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCFQ(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCAIR(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCSAT(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUCURR(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVCURR(:)
TYPE(TCST)        ,INTENT(IN)    :: YDCST
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PKMFL(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PKHFL(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PKQFL(:) 
INTEGER(KIND=JPIM) :: JL

REAL(KIND=JPRB) :: ZCDRO, ZCONS1, ZCONS16, ZKHFL
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!*    LOCAL STORAGE
!     ----- -------

!     ------------------------------------------------------------------

!*       1.     INITIALIZE CONSTANTS
!               ---------- ----------

IF (LHOOK) CALL DR_HOOK('VSFLX_MOD:VSFLX',0,ZHOOK_HANDLE)
ASSOCIATE(RCPD=>YDCST%RCPD, RD=>YDCST%RD, RETV=>YDCST%RETV, RG=>YDCST%RG)

ZCONS1 =RD/(RG*PTMST*PRVDIFTS)
ZCONS16=RCPD*RVTMP2

!     ------------------------------------------------------------------

!        2. COMPUTE DENSITY AND SURFACE FLUXES AT T-1
!           ------- ------- --- ------- ------ -- ---

DO JL=KIDIA,KFDIA
  ZCDRO    =ZCONS1*PTMLEV(JL)*(1.0_JPRB+RETV*PQMLEV(JL))/PAPHMS(JL)
  PKMFL(JL)=ZCDRO*PCFM(JL)*SQRT((PUMLEV(JL)-PUCURR(JL))**2&
   & +(PVMLEV(JL)-PVCURR(JL))**2)
  ZKHFL    =ZCDRO*PCFH(JL)*(PCPTGZLEV(JL)-PCPTS(JL))
  PKQFL(JL)=ZCDRO*PCFQ(JL)*(PCAIR(JL)*PQMLEV(JL)-PCSAT(JL)*PQS(JL))

  PKHFL(JL)=(ZKHFL-ZCONS16*PTMLEV(JL)*PKQFL(JL))&
   & /(RCPD*(1.0_JPRB+RVTMP2*PQMLEV(JL)))  

ENDDO

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('VSFLX_MOD:VSFLX',1,ZHOOK_HANDLE)
END SUBROUTINE VSFLX
END MODULE VSFLX_MOD
