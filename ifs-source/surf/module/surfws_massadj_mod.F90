MODULE SURFWS_MASSADJ_MOD
CONTAINS

SUBROUTINE SURFWS_MASSADJ(KIDIA, KFDIA, KLON, KLEVSN,    &
                       &   KLEVSNA, ZTHRESWS,             &
                       &  PDSN,PDSNREAL,PSNDEPTH,PSNDEPTHR,&
                       &  PRSN, PSSN, PRSNMAX, PDSNTOT,   &
                       &  PRCONSTAVG, PRMINCL,            &
                       &  ZSNPERT,                        &
                       &  PTSNWS,                         &
                       &  PRSNWS, PSSNWS,PWSNWS, PRSNTOP, &
                       &  YDCST,YDSOIL )

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOS_CST  , ONLY : TCST
USE YOS_SOIL , ONLY : TSOIL

USE ABORT_SURF_MOD

! (C) Copyright 2017- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!**** *SURFWS_MASSADJ* - Snow warm start multi-layer 
!     PURPOSE.
!     --------
!          THIS ROUTINE ADJUST SNOW DENSITY PROFILES 
!          TO CONSERVE TOTAL MASS OF SNOW LAYER

!**   INTERFACE.
!     ----------
!          *SURFWS_MASSADJ* IS CALLED FROM *SURFWS_CTL*.

!     PARAMETER   DESCRIPTION                                    UNITS
!     ---------   -----------                                    -----

!     INPUT PARAMETERS (INTEGER):
!    *KIDIA*      START POINT
!    *KFDIA*      END POINT
!    *KLON*       Length of arrays
!    *KLEVSN*     Snow vertical levels
!    *KLEVSNA*    Snow vertical Active levels
!    *PRMINCL*    Cluster index for density exp function

!     INPUT PARAMETERS for WARM START (REAL):
!    *ZTHRESWS*      snow depth threhsold for warm start
!    *ZSNPERT*       snow depth threshold for glaciers
!    *PDSNTOT*    Total snow depth                            (m)
!    *PSNDEPTH*       Snow depth per layer wrt 0              (m)
!    *PSNDEPTHR*   Snow depth per layer (full) wrt 0          (m)
!    *PDSN*       Snow depth per layer                        (m)
!    *PDSNREAL*   Snow depth per layer (full)                 (m)
!    *PRSNMAX*    Snow density MAX allowed                    (kg m-3)
!    *PRCONSTAVG* Constants for density exp function
!    *PRSNTOP*    Snow density top layer                      (kg m-3)


!     INPUT PARAMETERS PROGNOSTIC FROM SINGLE-LAYER (REAL):
!    *PSSN*       SNOW MASS        single layer               (kg m-2)
!    *PRSN*       SNOW DENSITY     single layer               (kg m-3)


!     OUTPUT PARAMETERS (REAL):
!    *PRSNWS*        Snow density    warm started
!    *PSSNWS*        Snow mass       warm started
!    *PWSNWS*        Snow liq water  warm started
!    *PTSNWS*        Snow temperature warm started



!     METHOD.
!     -------
!     Mass is adjusted in each layer by iteratively increase/decrease 
!     the snow density in order to conserve the total mass (input).
!     Snow liquid water is initialised here using the diagnostic
!     formulation assuming a temperature dependency of the liq water content
!     in each layer.

!     EXTERNALS.
!     ----------
!          NONE.

!     REFERENCE.
!     ----------
!          Arduini et al. (2019)

!     Modifications:
!     Original   G. Arduini      ECMWF     28/07/2017

!     ------------------------------------------------------------------

IMPLICIT NONE

! Input variables:
INTEGER(KIND=JPIM), INTENT(IN)  :: KLON, KIDIA, KFDIA
INTEGER(KIND=JPIM), INTENT(IN)  :: KLEVSN
INTEGER(KIND=JPIM), INTENT(IN)  :: KLEVSNA(:)

REAL(KIND=JPRB),    INTENT(IN)  :: PDSN(:,:), PDSNREAL(:,:), PSNDEPTH(:, :), PSNDEPTHR(:, :)
REAL(KIND=JPRB),    INTENT(IN)  :: PRSN(:), PSSN(:)
REAL(KIND=JPRB),    INTENT(IN)  :: PRSNMAX(:)
REAL(KIND=JPRB),    INTENT(IN)  :: PDSNTOT(:)
REAL(KIND=JPRB),    INTENT(IN)  :: PRCONSTAVG(:,:)
INTEGER(KIND=JPIM),    INTENT(IN)  :: PRMINCL(:)

REAL(KIND=JPRB),    INTENT(IN)  :: ZSNPERT  ! move to sussoil
REAL(KIND=JPRB),    INTENT(IN)  :: ZTHRESWS ! move to sussoil
REAL(KIND=JPRB),    INTENT(IN)  :: PTSNWS(:,:)

TYPE(TCST)     , INTENT(IN) :: YDCST
TYPE(TSOIL)    , INTENT(IN) :: YDSOIL


! Output
REAL(KIND=JPRB),    INTENT(INOUT) :: PRSNWS(:,:), PSSNWS(:,:),PWSNWS(:,:)
REAL(KIND=JPRB),    INTENT(INOUT) :: PRSNTOP(:)

! Local variables
REAL(KIND=JPRB)                   :: ZRSNWSTST(KLON,KLEVSN)
REAL(KIND=JPRB)                   :: ZERRORTST, ZERROR
REAL(KIND=JPRB)                   :: ZRBOTTOM
REAL(KIND=JPRB)                   :: ZRBOTTOMTST
REAL(KIND=JPRB)                   :: ZRBOTTOMSTORE
REAL(KIND=JPRB)                   :: ZRSNTOPTST
REAL(KIND=JPRB)                   :: ZRSNTOPSTORE
REAL(KIND=JPRB)                   :: ZDELTAC
REAL(KIND=JPRB)                   :: ZDRSN
REAL(KIND=JPRB)                   :: ZBR, ZAR, ZCR

REAL(KIND=JPRB)                   :: ZLWC , ZFUNC
INTEGER(KIND=JPIM) :: KFILLING
INTEGER(KIND=JPIM) :: ICOUNTMAX, ICOUNT
INTEGER(KIND=JPIM) :: JL,JK
REAL(KIND=JPRB)                   :: ZEPSILON

REAL(KIND=JPHOOK)                  :: ZHOOK_HANDLE

! INCLUDE FUNCTIONS
#include "fcsurf.h"

!    -----------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SURFWS_MASSADJ_MOD:SURFWS_MASSADJ',0,ZHOOK_HANDLE)

!    -----------------------------------------------------------------

ASSOCIATE(RTT=>YDCST%RTT,RPI=>YDCST%RPI, &
        & RHOMINSND=>YDSOIL%RHOMINSND,RTEMPAMP=>YDSOIL%RTEMPAMP )

ZEPSILON  = 10E4*EPSILON(ZEPSILON)
KFILLING  = 1._JPIM
! 2.3.1 snow density: recompute profiles spanning bottom snow density value:
!       Looking for value which minimizes the error in mass:
DO JL=KIDIA, KFDIA
  IF (PDSNTOT(JL) >= ZTHRESWS) THEN
    ZRSNWSTST(JL,1:KLEVSN)=PRSNWS(JL,1:KLEVSN)
    ZERRORTST = 0._JPRB
    ZRBOTTOMTST=PRSNWS(JL, KLEVSNA(JL))
    ZRBOTTOMSTORE=PRSNWS(JL, KLEVSNA(JL))
    ZDELTAC   = 0._JPRB
    ZDRSN     = 50.0_JPRB 
    ICOUNTMAX = 50_JPIM*2_JPIM+1_JPIM 
    ICOUNT    = 1_JPIM

! HERE I AM NOT USING PDSNREAL... to be re-evaluated...
    !ZERROR    = ( SUM( ( PRSNWS(JL, 1:KLEVSNA(JL)) * PDSN(JL, 1:KLEVSNA(JL)) ), DIM=1 ) - PSSN(JL) ) / PSSN(JL)
    ZERROR    = ( SUM( ( PRSNWS(JL, 1:KLEVSNA(JL)) * PDSNREAL(JL, 1:KLEVSNA(JL)) ), DIM=1 ) - PSSN(JL) ) / PSSN(JL)
  
    ZCR    = PRCONSTAVG(JL, PRMINCL(JL) )
    DO WHILE ( (ABS(ZERROR) > ZEPSILON) .AND. (ICOUNT < ICOUNTMAX) )
      ZDELTAC =ZDRSN - (ICOUNT)*ZDRSN/(0.5_JPRB*ICOUNTMAX)
      ZRBOTTOMTST=PRSNWS(JL, KLEVSNA(JL)) - ZDELTAC
  
      ZBR=(ZRBOTTOMTST-PRSNTOP(JL))/SIGN(MAX(ZEPSILON,             &
                                &ABS(MAX(ZEPSILON,EXP(-PDSNTOT(JL)*ZCR))-1)),&
                                &(MAX(ZEPSILON,EXP(-PDSNTOT(JL)*ZCR))-1)     )
     !ZAR=PRSNTOP(JL) - ZBR*EXP(-PSNDEPTH(JL,1)*ZCR)
      ZAR=PRSNTOP(JL) - ZBR*EXP(-PSNDEPTHR(JL,1)*ZCR)
      DO JK=1, KLEVSNA(JL)
        !ZRSNWSTST(JL, JK) = ZAR + ZBR *  MAX(ZEPSILON,EXP( - PSNDEPTH(JL, JK) * ZCR )) 
        ZRSNWSTST(JL, JK) = ZAR + ZBR *  MAX(ZEPSILON,EXP( - PSNDEPTHR(JL, JK) * ZCR )) 
        IF (ZRSNWSTST(JL, JK) < RHOMINSND) ZRSNWSTST(JL, JK) = RHOMINSND
        IF (ZRSNWSTST(JL, JK) > PRSNMAX(JL) ) ZRSNWSTST(JL, JK) = PRSNMAX(JL)
      ENDDO
      IF (KLEVSNA(JL) < KLEVSN) THEN
        ZRSNWSTST(JL, KLEVSNA(JL)+1:KLEVSN)=100._JPRB
      ENDIF
! 
      ZERRORTST=( SUM( ( ZRSNWSTST(JL, 1:KLEVSNA(JL)) * PDSNREAL(JL, 1:KLEVSNA(JL)) ), DIM=1 ) - PSSN(JL) ) / PSSN(JL)
      IF (ABS(ZERRORTST) < ABS(ZERROR)) THEN
        PRSNWS(JL,1:KLEVSNA(JL))=ZRSNWSTST(JL,1:KLEVSNA(JL))
        ZERROR=ZERRORTST
        ZRBOTTOMSTORE=ZRBOTTOMTST
      ENDIF
      ICOUNT = ICOUNT + 1_JPIM
    ENDDO
    ZRBOTTOM=ZRBOTTOMSTORE
  
! 2.3.2 snow density: recompute profiles spanning top snow density value:
!       Looking for value which minimizes the error in mass:
    IF ( PSSN(JL) < ZSNPERT ) THEN
      ZRSNWSTST(JL,:) = PRSNWS(JL,:)
      ZRSNTOPSTORE    = PRSNTOP(JL)
      ZRSNTOPTST      = PRSNTOP(JL)
      ZDELTAC         = 0._JPRB
      ZERRORTST       = 0._JPRB
    
      ZERROR    = ( SUM( ( PRSNWS(JL, 1:KLEVSNA(JL)) * PDSNREAL(JL, 1:KLEVSNA(JL)) ), DIM=1 ) - PSSN(JL) ) / PSSN(JL)
      ZDRSN     = 20.0_JPRB 
      ICOUNTMAX = INT(ZDRSN)*2_JPIM+1_JPIM 
      ICOUNT    = 1_JPIM
      DO WHILE ( (ABS(ZERROR) > ZEPSILON) .AND. (ICOUNT < ICOUNTMAX) )
        ZDELTAC   = ZDRSN - (ICOUNT)*ZDRSN/(0.5_JPRB*ICOUNTMAX)
        ZRSNTOPTST= PRSNTOP(JL) - ZDELTAC
    
        ZBR=(ZRBOTTOM-ZRSNTOPTST)/SIGN(MAX(ZEPSILON,            &
                                 &ABS(MAX(ZEPSILON,EXP(-PDSNTOT(JL)*ZCR))-1)),&
                                 &(MAX(ZEPSILON,EXP(-PDSNTOT(JL)*ZCR))-1)     )
        !ZAR=ZRSNTOPTST - ZBR*MAX(ZEPSILON,EXP(-PSNDEPTH(JL,1)*ZCR))
        ZAR=ZRSNTOPTST - ZBR*MAX(ZEPSILON,EXP(-PSNDEPTHR(JL,1)*ZCR))
        DO JK=1, KLEVSNA(JL)
          !ZRSNWSTST(JL, JK) = ZAR + ZBR *  exp( - PSNDEPTH(JL, JK) * ZCR ) 
          ZRSNWSTST(JL, JK) = ZAR + ZBR *  exp( - PSNDEPTHR(JL, JK) * ZCR ) 
          IF (ZRSNWSTST(JL, JK) < RHOMINSND) THEN
            ZRSNWSTST(JL, JK) = RHOMINSND
          ENDIF
          IF (ZRSNWSTST(JL, JK) > PRSNMAX(JL) ) THEN
            ZRSNWSTST(JL, JK) = PRSNMAX(JL)
          ENDIF
        ENDDO
        IF (KLEVSNA(JL) < KLEVSN) THEN
          ZRSNWSTST(JL, KLEVSNA(JL)+1:KLEVSN)=100._JPRB
        ENDIF
        ZERRORTST=( SUM( ( ZRSNWSTST(JL, 1:KLEVSNA(JL)) * PDSNREAL(JL, 1:KLEVSNA(JL)) ), DIM=1 ) - PSSN(JL) ) / PSSN(JL)
        IF (ABS(ZERRORTST) < ABS(ZERROR)) THEN
          PRSNWS(JL,:)=ZRSNWSTST(JL,:)
          ZERROR=ZERRORTST
          ZRSNTOPSTORE=ZRSNTOPTST
        ENDIF
        ICOUNT = ICOUNT + 1_JPIM
      ENDDO
      PRSNTOP(JL)=ZRSNTOPSTORE

! 2.3.2 snow density: filling 
      ICOUNTMAX = 1000_JPIM 
      ZERROR   = ( SUM( ( PRSNWS(JL, 1:KLEVSNA(JL)) * PDSNREAL(JL, 1:KLEVSNA(JL)) ), DIM=1 ) - PSSN(JL) ) / PSSN(JL)
  
      ICOUNT = 1_JPIM
      DO WHILE ( (ABS(ZERROR) > ZEPSILON) .AND. (ICOUNT < ICOUNTMAX) )
        IF (PDSNTOT(JL) >= 0.15_JPRB) THEN
          KFILLING=MAX(1,KLEVSNA(JL)-1_JPIM)
        ELSEIF (PDSNTOT(JL) < 0.15_JPRB) THEN
          KFILLING=MIN(2_JPIM, KLEVSNA(JL))
        ENDIF
        IF ( ZERROR < 0._JPRB ) THEN ! LESS MASS
          IF (KFILLING > 1) THEN
            PRSNWS(JL, KFILLING:KLEVSNA(JL)) = PRSNWS(JL, KFILLING:KLEVSNA(JL)) + 0.1_JPRB
            IF (KFILLING > 2) THEN
              PRSNWS(JL, 2:KFILLING-1)         = PRSNWS(JL, 2:KFILLING-1) + 0.005_JPRB
            ENDIF
          ELSE
            PRSNWS(JL, KFILLING)             = PRSNWS(JL, KFILLING) + 0.1_JPRB
          ENDIF
  
        ELSEIF ( ZERROR > 0._JPRB ) THEN
          IF (KFILLING > 1) THEN
            IF (PRSNWS(JL, KFILLING) > PRSNWS(JL, KFILLING-1) ) THEN
              PRSNWS(JL, KFILLING:KLEVSNA(JL) ) = PRSNWS(JL, KFILLING:KLEVSNA(JL)) - 0.1_JPRB
              PRSNWS(JL, 1:KFILLING-1)   = PRSNWS(JL, 1:KFILLING-1) - 0.01_JPRB
            ELSE
              PRSNWS(JL, 1:KLEVSNA(JL) ) = PRSNWS(JL, 1:KLEVSNA(JL)) - 0.01_JPRB
            ENDIF
          ELSE
              PRSNWS(JL, KFILLING ) = PRSNWS(JL, KFILLING) - 0.1_JPRB
          ENDIF
        ENDIF
      
        ZERROR   = ( SUM( ( PRSNWS(JL, 1:KLEVSNA(JL)) * PDSNREAL(JL, 1:KLEVSNA(JL))), DIM=1 ) - PSSN(JL) ) / PSSN(JL)
        ICOUNT = ICOUNT + 1_JPIM
      END DO
  
    ELSE ! Glacier filling:
         ! We fill/remove density (mass) from the first three layers, the active
         ! ones here so why 1:KLEVSNA(JL)-2 over glacier. The remaining mass is
         ! put in the accumulation layer KLEVSNA(JL)-1
      ICOUNTMAX = 1000_JPIM 
      ZERROR   = ( SUM( ( PRSNWS(JL, 1:KLEVSNA(JL)) * PDSNREAL(JL, 1:KLEVSNA(JL)) ), DIM=1 ) - PSSN(JL) ) / PSSN(JL)
      ICOUNT = 1_JPIM
      DO WHILE ( (ABS(ZERROR) > ZEPSILON) .AND. (ICOUNT < ICOUNTMAX) )
        KFILLING=1
        IF ( ZERROR < 0._JPRB ) THEN ! LESS MASS
          PRSNWS(JL, 1:KLEVSNA(JL)-2)   = PRSNWS(JL, 1:KLEVSNA(JL)-2) + 0.01_JPRB
        ELSEIF ( ZERROR > 0._JPRB ) THEN    
          PRSNWS(JL, 1:KLEVSNA(JL)-2)   = PRSNWS(JL, 1:KLEVSNA(JL)-2) - 0.01_JPRB
        ENDIF
        ZERROR   = ( SUM( ( PRSNWS(JL, 1:KLEVSNA(JL)) * PDSNREAL(JL, 1:KLEVSNA(JL))), DIM=1 ) - PSSN(JL) ) / PSSN(JL)
        ICOUNT = ICOUNT + 1_JPIM
      ENDDO
    ENDIF

!**********************************************************************
! 2.3.2 Compute snow mass in each layer
!       Update now the snow mass in each layer accordingly to its depth
!       For glaciers put in the accumulation layer KLEVSNA(JL)-1
    DO JK=1, KLEVSNA(JL)
      PSSNWS(JL, JK)    = ( PRSNWS(JL, JK) * PDSNREAL(JL, JK) )
    ENDDO
    IF  (SUM(PSSNWS(JL, 1:KLEVSNA(JL)), DIM=1) .NE. PSSN(JL) ) THEN
      ZERROR = (SUM(PSSNWS(JL, 1:KLEVSNA(JL)), DIM=1)-PSSN(JL))/PSSN(JL)
      DO JK=1, KLEVSNA(JL)-2
        PSSNWS(JL, JK) = PSSNWS(JL, JK) - ( ( ZERROR * PSSN(JL) ) * PDSNREAL(JL, JK) / SUM(PDSNREAL(JL,1:KLEVSNA(JL)),DIM=1) )
        PRSNWS(JL, JK) = PSSNWS(JL, JK) / PDSNREAL(JL, JK)
        IF (PRSNWS(JL, JK)  < RHOMINSND)  THEN
          PRSNWS(JL, JK) = RHOMINSND
          PSSNWS(JL, JK) = ( PRSNWS(JL, JK) * PDSNREAL(JL, JK) )
        ENDIF
        IF (PRSNWS(JL, JK) > PRSNMAX(JL) ) THEN
          PRSNWS(JL, JK) = PRSNMAX(JL)
          PSSNWS(JL, JK) = ( PRSNWS(JL, JK) * PDSNREAL(JL, JK) )
        ENDIF
      ENDDO
      ! double-ckeck for glaciers, adding to the accumulation layer
      IF  (SUM(PSSNWS(JL, 1:KLEVSNA(JL)), DIM=1) .NE. PSSN(JL) ) THEN
        ZERROR = (SUM(PSSNWS(JL, 1:KLEVSNA(JL)), DIM=1)-PSSN(JL))/PSSN(JL)
        PSSNWS(JL, KLEVSNA(JL)-1) = PSSNWS(JL, KLEVSNA(JL)-1) - ( ZERROR * PSSN(JL) )
      ENDIF
    ENDIF

!**********************************************************************
! 2.3.3 Update the liquid water part.
! We diagnose slw using the diagnostic formulation used in srfsn_lwimp for each
! layer
    IF ( PSSN(JL) < ZSNPERT ) THEN
      DO JK=1, KLEVSNA(JL)
      ! SNOW LIQUID WATER CAPACITY
        ZLWC     = FLWC( PSSNWS(JL, JK), PRSNWS(JL, JK) )
  
      ! ANALYTICAL FUNCTIONS - SNOW LIQUID WATER
        IF ( PTSNWS(JL, JK) < RTT-0.5_JPRB*RTEMPAMP ) THEN
          ZFUNC  = 0._JPRB
        ELSEIF (PTSNWS(JL, JK) < RTT ) THEN  
          ZFUNC  = (1._JPRB + SIN(RPI*(PTSNWS(JL,JK)-RTT)/(RTEMPAMP)))
        ELSE
          ZFUNC  = 0._JPRB
        ENDIF
  
        PWSNWS(JL, JK)  = ZFUNC*ZLWC
      ! PWSNWS(JL, JK)  = ZLWC
      ENDDO
    ELSE 
      PWSNWS(JL, :)  = 0._JPRB
    ENDIF

  ENDIF ! END IF WARM START THR
END DO  ! END DO IN JL


END ASSOCIATE

!    -----------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SURFWS_MASSADJ_MOD:SURFWS_MASSADJ',1,ZHOOK_HANDLE)

END SUBROUTINE SURFWS_MASSADJ
END MODULE SURFWS_MASSADJ_MOD


