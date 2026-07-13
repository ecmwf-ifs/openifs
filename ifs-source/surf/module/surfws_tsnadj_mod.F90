MODULE SURFWS_TSNADJ_MOD
CONTAINS

SUBROUTINE SURFWS_TSNADJ(KIDIA, KFDIA,KLON,KLEVSN,                  &
                       & KLEVSNA, KLEVMID, ZTHRESWS,                &
                       & PSNDEPTH, PTCONSTAVG, PTCONSTSTD, PTMINCL, &
                       & PRSN, PSSN, PTSN, PHSN,                    &
                       & PTSKIN, PDSNTOT, PACTDEPTH, PSADEPTH,      &
                       & PTSNBOTTOM, PTSNTOP,PTSNMIDDLE,            &
                       & ZSNPERT,                                   &
                       & PTSNWS, PSSNWS,PWSNWS,                     &
                       & YDCST, YDSOIL                              )


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

!**** *SURFWS_TSNADJ* - Snow warm start multi-layer 
!     PURPOSE.
!     --------
!          THIS ROUTINE CONTROLS THE WARM START OF MULTI-LAYER SCHEME IN
!          FC MODE

!**   INTERFACE.
!     ----------
!          *SURFWS_CTL* IS CALLED FROM *SRFSN_DRIVER*.

!     PARAMETER   DESCRIPTION                                    UNITS
!     ---------   -----------                                    -----

!     INPUT PARAMETERS (INTEGER):
!    *KIDIA*      START POINT


!     INPUT PARAMETERS (REAL):
!    *PTMST*      TIME STEP                                      S

!     INPUT PARAMETERS (LOGICAL):
!    *LDLAND*     LAND/SEA MASK (TRUE/FALSE)

!     INPUT PARAMETERS AT T-1 OR CONSTANT IN TIME (REAL):
!    *PWLM1M*     SKIN RESERVOIR WATER CONTENT                   kg/m**2


!     OUTPUT PARAMETERS AT T+1 (UNFILTERED,REAL):
!    *PWL*        SKIN RESERVOIR WATER CONTENT                   kg/m**2


!     OUTPUT PARAMETERS (DIAGNOSTIC):

!    *PDHIIS*     Diagnostic array for interception layer (see module yomcdh)

!     METHOD.
!     -------
!          

!     EXTERNALS.
!     ----------
!          NONE.

!     REFERENCE.
!     ----------
!          

!     Modifications:
!     Original   G. Arduini      ECMWF     28/07/2017

!     ------------------------------------------------------------------

IMPLICIT NONE

! Input variables:
INTEGER(KIND=JPIM), INTENT(IN)  :: KIDIA, KFDIA,KLON
INTEGER(KIND=JPIM), INTENT(IN)  :: KLEVSN
INTEGER(KIND=JPIM), INTENT(IN)  :: KLEVSNA(:),KLEVMID(:)

REAL(KIND=JPRB),    INTENT(IN)  :: PSNDEPTH(:, :)
REAL(KIND=JPRB),    INTENT(IN)  :: PTCONSTAVG(:,:),PTCONSTSTD(:,:)
INTEGER(KIND=JPIM), INTENT(IN)  :: PTMINCL(:)
REAL(KIND=JPRB),    INTENT(IN)  :: PRSN(:), PSSN(:), PTSN(:), PHSN(:)  
REAL(KIND=JPRB),    INTENT(IN)  :: PTSKIN(:)
REAL(KIND=JPRB),    INTENT(IN)  :: PDSNTOT(:), PACTDEPTH(:), PSADEPTH(:)
REAL(KIND=JPRB),    INTENT(IN)  :: PTSNBOTTOM(:), PTSNTOP(:),PTSNMIDDLE(:)
REAL(KIND=JPRB),    INTENT(IN)  :: ZSNPERT   ! move to sussoil
REAL(KIND=JPRB),    INTENT(IN)  :: ZTHRESWS  ! move to sussoil
                                    
TYPE(TCST)     , INTENT(IN) :: YDCST
TYPE(TSOIL)    , INTENT(IN) :: YDSOIL

! Input/Output
REAL(KIND=JPRB),    INTENT(INOUT) :: PTSNWS(:,:), PWSNWS(:,:),PSSNWS(:,:)


! Local variables
!REAL(KIND=JPRB)    :: ZTSNWSTST(KLON,0:KLEVSN+1)
REAL(KIND=JPRB)    :: ZTSNWSTST(KLON,KLEVSN)
REAL(KIND=JPRB)    :: ZWSNWSTST(KLON,KLEVSN)
REAL(KIND=JPRB)    :: ZHSNWS(KLON)
REAL(KIND=JPRB)    :: ZERRORTST, ZERROR
REAL(KIND=JPRB)    :: ZEPSILON
REAL(KIND=JPRB)    :: ZEPSILONL
REAL(KIND=JPRB)    :: ZDELTAC
REAL(KIND=JPRB)    :: ZAT, ZBT
REAL(KIND=JPRB)    :: ZCT,ZCTTST
REAL(KIND=JPRB)    :: ZHSNWSTST
REAL(KIND=JPRB)    :: ZTINC
REAL(KIND=JPRB)    :: ZTINCTOP
REAL(KIND=JPRB)    :: ZWINC
REAL(KIND=JPRB)    :: ZIHCAP

INTEGER(KIND=JPIM) :: JL, JK
INTEGER(KIND=JPIM) :: ICOUNTMAX, ICOUNT
INTEGER(KIND=JPIM) :: KFILLING, KBOTTOM

REAL(KIND=JPHOOK)                  :: ZHOOK_HANDLE

! INCLUDE FUNCTIONS
#include "fcsurf.h"

!    -----------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SURFWS_TSNADJ_MOD:SURFWS_TSNADJ',0,ZHOOK_HANDLE)

!    -----------------------------------------------------------------

ASSOCIATE(RTT=>YDCST%RTT, RLMLT=>YDCST%RLMLT  )

ZIHCAP    = YDSOIL%RHOCI / YDSOIL%RHOICE
ZEPSILON  = 10E4*EPSILON(ZEPSILON)
ZEPSILONL = 1E-5


!**************************************************************
! 2.3.4 snow temperature: recompute profiles spanning C-const snow temperature:
!       Looking for value which minimizes the error in energy:
DO JL=KIDIA, KFDIA

  IF (PDSNTOT(JL) >= ZTHRESWS) THEN
    IF ( PSSN(JL) < ZSNPERT ) THEN ! not a glacier point
      ZTSNWSTST(JL,:)=PTSNWS(JL,:)
      ZWSNWSTST(JL,:)=PWSNWS(JL,:)
  
      ZCT    = PTCONSTAVG(JL,PTMINCL(JL))
      ZCTTST = ZCT
  
      ZHSNWS(JL)  = SUM( ZIHCAP*PSSNWS(JL, 1:KLEVSNA(JL)) * (PTSNWS(JL, 1:KLEVSNA(JL)) - RTT) &
                   &- RLMLT * (PSSNWS(JL,1:KLEVSNA(JL))-PWSNWS(JL, 1:KLEVSNA(JL))), DIM=1 )
  
      ZERROR=(ZHSNWS(JL)-PHSN(JL))/(ABS(PHSN(JL)))
      ICOUNTMAX = 50_JPIM*2_JPIM+1_JPIM !INT(2*ZRCONSTSTD(PTMINCL(JL))/20.0, JPIM)+1
      ICOUNT    = 1_JPIM
      DO WHILE ( (ABS(ZERROR) > ZEPSILONL) .AND. (ICOUNT < ICOUNTMAX) )
        ZDELTAC = PTCONSTSTD(JL,PTMINCL(JL)) - (ICOUNT)*PTCONSTSTD(JL,PTMINCL(JL)) / (0.5_JPRB*ICOUNTMAX)
        ZCTTST  = PTCONSTAVG(JL,PTMINCL(JL)) - ZDELTAC
! Recompute profiles:
        IF ( PDSNTOT(JL) < 0.20_JPRB ) THEN
          ZBT=(PTSNBOTTOM(JL)-PTSNTOP(JL)) / SIGN(MAX(ZEPSILON,                 &
                                            &ABS(MAX(ZEPSILON,EXP(-PDSNTOT(JL)*ZCTTST ))-1.0_JPRB) ),&
                                            &(MAX(ZEPSILON,EXP(-PDSNTOT(JL)*ZCTTST ))-1.0_JPRB)      )
          ZAT=PTSNTOP(JL) - ZBT*MAX(ZEPSILON,EXP( -PSNDEPTH(JL, 1)*ZCTTST ))
        ELSE
        !!  IF (PTSNTOP(JL)>=PTSN(JL) .OR. (PSSN(JL) > ZSNPERT)) THEN
          IF (PTSNTOP(JL)>=PTSNMIDDLE(JL) .OR. (PSSN(JL) > ZSNPERT)) THEN
            ZBT=(PTSNBOTTOM(JL) - PTSNTOP(JL)) / SIGN(MAX(ZEPSILON,                                        &
                                                &ABS(MAX(ZEPSILON,EXP(-PACTDEPTH(JL) * ZCTTST))-1.0_JPRB)),&
                                                &(MAX(ZEPSILON,EXP(-PACTDEPTH(JL) * ZCTTST ))-1.0_JPRB)    )
            IF (PTSNTOP(JL)<=RTT) THEN
              ZAT=PTSNTOP(JL)-ZBT
            ELSE 
              ZAT=RTT-ZBT
            ENDIF
          ELSE
            ZBT=(PTSNBOTTOM(JL)-PTSNTOP(JL)) / SIGN(MAX(ZEPSILON,                                  &
                                              &ABS(MAX(ZEPSILON,EXP(-(PDSNTOT(JL)+PSADEPTH(JL)) * ZCTTST))-1.0_JPRB)),&
                                              &(MAX(ZEPSILON,EXP(-(PDSNTOT(JL)+PSADEPTH(JL)) * ZCTTST)) - 1.0_JPRB) )
            ZAT=PTSNTOP(JL) - ZBT
          ENDIF
        ENDIF
        DO JK=1, KLEVSNA(JL)
          ! -- Temperature profile parametrisation
          ZTSNWSTST(JL, JK) = ZAT + ZBT * MAX(ZEPSILON,EXP( -PSNDEPTH(JL, JK) * ZCTTST )) 
         ! Safety checks:
          IF (ZTSNWSTST(JL, JK) > RTT ) ZTSNWSTST(JL, JK) = RTT
          IF (PDSNTOT(JL) >= 0.5_JPRB) THEN
            ZTSNWSTST(JL, KLEVSNA(JL)-2:KLEVSNA(JL))=PTSNMIDDLE(JL) !PTSN(JL)
          ELSEIF (PDSNTOT(JL) >= 0.20_JPRB ) THEN
            ZTSNWSTST(JL, KLEVMID(JL))=PTSNMIDDLE(JL) !PTSN(JL)
          ELSEIF (PDSNTOT(JL) < 0.20_JPRB ) THEN
            ZTSNWSTST(JL, KLEVMID(JL))=PTSNMIDDLE(JL) !PTSN(JL)
          ENDIF
          IF ( PSSN(JL) >= ZSNPERT ) ZTSNWSTST(JL, KLEVSNA(JL))=PTSNBOTTOM(JL)
   
        ENDDO
   
        IF (KLEVSNA(JL) < KLEVSN) THEN
          ZTSNWSTST(JL, KLEVSNA(JL)+1:KLEVSN)=ZTSNWSTST(JL, KLEVSNA(JL))
        ENDIF
        ZHSNWSTST  = SUM( ZIHCAP*PSSNWS(JL, 1:KLEVSNA(JL)) * (ZTSNWSTST(JL, 1:KLEVSNA(JL)) - RTT) &
                     &- RLMLT * (PSSNWS(JL,1:KLEVSNA(JL))-ZWSNWSTST(JL, 1:KLEVSNA(JL))), DIM=1 )
   
        ZERRORTST=(ZHSNWSTST-PHSN(JL))/(ABS(PHSN(JL)))
   
        IF (ABS(ZERRORTST) < ABS(ZERROR)) THEN
          PTSNWS(JL,:) = ZTSNWSTST(JL,:)
          PWSNWS(JL,:) = ZWSNWSTST(JL,:)
          ZERROR       = ZERRORTST
          ZCT          = ZCTTST
        ENDIF
        ICOUNT = ICOUNT + 1_JPIM
      ENDDO
    ENDIF

!***********************************************************
! 2.3.4 Check that heat (energy) content is conserved 
    ZTSNWSTST(JL,:) = PTSNWS(JL,:)
    ZWSNWSTST(JL,:) = PWSNWS(JL,:)
    ZHSNWS(JL)  = SUM( ZIHCAP*PSSNWS(JL, 1:KLEVSNA(JL)) * (PTSNWS(JL, 1:KLEVSNA(JL)) - RTT) -&
                &     RLMLT * (PSSNWS(JL,1:KLEVSNA(JL))-PWSNWS(JL, 1:KLEVSNA(JL))), DIM=1 )
  
    ZERROR=(ZHSNWS(JL)-PHSN(JL))/(ABS(PHSN(JL)))
    !PRINT*,'----JL,',JL,' HEAT CONTENT BEFORE Heat cont adjustment----'
    !PRINT*,'ZHSNWS=',ZHSNWS(JL),' ZHSN IN=',PHSN(JL), 'rel error=',ZERROR
    ICOUNT = 1_JPIM
    ICOUNTMAX=1000_JPIM

    !!-- Final adjustment temperature and water profiles: 
    !!   Heat filling or removing
    !!   - moved ICOUNMAX=100
    ZTINC = 0.005_JPRB
    ZTINCTOP=ZTINC
    ZWINC = 0.01_JPRB
    DO WHILE ( ( ABS(ZERROR) .GT. ZEPSILONL ) .AND. (ICOUNT .LT. ICOUNTMAX) )
      IF ( PSSN(JL) < ZSNPERT ) THEN
        IF (PDSNTOT(JL) >= 0.15_JPRB) THEN
          KFILLING=KLEVSNA(JL)-1_JPIM
          KBOTTOM=1_JPIM
        ELSEIF (PDSNTOT(JL) < 0.15_JPRB) THEN
          KFILLING=2_JPIM
          KBOTTOM=2_JPIM
        ENDIF
        IF ( PTSNWS(JL, 1) < (PTSKIN(JL)-1._JPRB) ) THEN
          KBOTTOM=2_JPIM
        ENDIF
        IF ( ZERROR > 0._JPRB ) THEN
          ! Temperature::
          IF (KFILLING > 1) THEN
            ZTSNWSTST(JL, KFILLING:KLEVSNA(JL)) = ZTSNWSTST(JL, KFILLING:KLEVSNA(JL)) - ZTINC 
            ZTSNWSTST(JL, KBOTTOM:KFILLING-1)   = ZTSNWSTST(JL, KBOTTOM:KFILLING-1) - ZTINCTOP
            IF (ZTSNWSTST(JL, KLEVSNA(JL)) < (ZTSNWSTST(JL, KLEVSNA(JL)-1))) THEN
              ZTSNWSTST(JL, KLEVSNA(JL)) = ZTSNWSTST(JL, KLEVSNA(JL)-1)
            ENDIF
          ELSE
            ZTSNWSTST(JL, KFILLING) = ZTSNWSTST(JL, KFILLING) - ZTINCTOP
          ENDIF
          ! Water::
          WHERE( ZWSNWSTST(JL, :) > 0._JPRB )
            ZWSNWSTST(JL, 1:KLEVSNA(JL) ) = ZWSNWSTST(JL, 1:KLEVSNA(JL)) - ZWINC 
          ELSEWHERE
            ZWSNWSTST(JL, 1:KLEVSNA(JL) ) = 0._JPRB
          ENDWHERE
          
        ELSEIF ( ZERROR < 0._JPRB ) THEN
          IF (KFILLING > 1) THEN
            ZTSNWSTST(JL, KFILLING:KLEVSNA(JL)) = ZTSNWSTST(JL, KFILLING:KLEVSNA(JL)) + ZTINC 
            ZTSNWSTST(JL, KBOTTOM:KFILLING-1)   = ZTSNWSTST(JL, KBOTTOM:KFILLING-1)   + ZTINCTOP
          ELSE
            ZTSNWSTST(JL, KFILLING) = ZTSNWSTST(JL, KFILLING) + ZTINCTOP
          ENDIF
          WHERE ( ZTSNWSTST(JL, :) > RTT )
            ZTSNWSTST(JL, :) = RTT
          ENDWHERE
          IF (KLEVSNA(JL) < KLEVSN) THEN
            ZTSNWSTST(JL, KLEVSNA(JL)+1:KLEVSN) = ZTSNWSTST(JL, KLEVSNA(JL))
          ENDIF
        ENDIF
  
      ELSE ! Glacier heat content ''filling'':
        KFILLING=1
        IF ( ZERROR > 0._JPRB ) THEN
          ZTSNWSTST(JL, KFILLING:KLEVSNA(JL)-2) = ZTSNWSTST(JL, KFILLING:KLEVSNA(JL)-2) - ZTINC 
        ELSEIF ( ZERROR < 0._JPRB ) THEN
          ZTSNWSTST(JL, KFILLING:KLEVSNA(JL)-2) = ZTSNWSTST(JL, KFILLING:KLEVSNA(JL)-2) + ZTINC 
        ENDIF 
      ENDIF
  
      ZHSNWSTST  = SUM( ZIHCAP*PSSNWS(JL, 1:KLEVSNA(JL)) * (ZTSNWSTST(JL, 1:KLEVSNA(JL)) - RTT) - RLMLT * (PSSNWS(JL,1:KLEVSNA(JL))-ZWSNWSTST(JL, 1:KLEVSNA(JL))), DIM=1 )
      ZERRORTST=(ZHSNWSTST-PHSN(JL))/(ABS(PHSN(JL)))
      IF (ABS(ZERRORTST) <= ABS(ZERROR)) THEN
        PTSNWS(JL,:) = ZTSNWSTST(JL,:)
        PWSNWS(JL,:) = ZWSNWSTST(JL,:)
        ZERROR=ZERRORTST
        ICOUNT = ICOUNT + 1_JPIM
      ! PRINT*,'----HEAT CONT FILLING----'
      ! PRINT*,PTSNWS(JL,:)
      ELSE
        ZTINCTOP  = - 0.5_JPRB*ZTINCTOP 
        ZTINC     = - 0.5_JPRB*ZTINC
        ICOUNT = ICOUNT + 1_JPIM
      ENDIF
    END DO

  END IF ! END IF IN SURFWS THR
END DO ! END LOOP IN JL

END ASSOCIATE

!    -----------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SURFWS_TSN_ADJ_MOD:SURFWS_TSNADJ',1,ZHOOK_HANDLE)

END SUBROUTINE SURFWS_TSNADJ
END MODULE SURFWS_TSNADJ_MOD
