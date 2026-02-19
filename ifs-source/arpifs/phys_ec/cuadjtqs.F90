! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE CUADJTQS &
 &(KIDIA,    KFDIA,    KLON,    KLEV,&
 & KK,&
 & PSP,      PT,       PQ,       LDFLAG,   KCALL)  

!**   *CUADJTQS* - SIMPLIFIED VERSION OF MOIST ADJUSTMENT

!     PURPOSE.
!     --------
!     TO PRODUCE T,Q AND L VALUES FOR CLOUD ASCENT

!     INTERFACE
!     ---------
!     THIS ROUTINE IS CALLED FROM SUBROUTINES:

!       *COND*       
!       *CUBMADJ*    
!       *CUBMD*      
!       *CONDAD*     
!       *CUBMADJAD*  
!       *CUBMDAD*    

!     INPUT ARE UNADJUSTED T AND Q VALUES,
!     IT RETURNS ADJUSTED VALUES OF T AND Q

!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):

!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KLEV*         NUMBER OF LEVELS
!    *KK*           LEVEL
!    *KCALL*        DEFINES CALCULATION AS
!                      KCALL=0  ENV. T AND QS IN*CUINI*
!                      KCALL=1  CONDENSATION IN UPDRAFTS  (E.G. CUBASE, CUASC)
!                      KCALL=2  EVAPORATION IN DOWNDRAFTS (E.G. CUDLFS,CUDDRAF)

!     INPUT PARAMETERS (LOGICAL):

!    *LDLAND*       LAND-SEA MASK (.TRUE. FOR LAND POINTS)

!     INPUT PARAMETERS (REAL):

!    *PSP*          PRESSURE                                        PA

!     UPDATED PARAMETERS (REAL):

!    *PT*           TEMPERATURE                                     K
!    *PQ*           SPECIFIC HUMIDITY                             KG/KG

!     AUTHOR.
!     -------
!      J.F. MAHFOUF      ECMWF         

!     MODIFICATIONS.
!     --------------
!      M.Hamrud     01-Oct-2003 CY28 Cleaning  
!      20180303 : Gabor: Just a comment line to force recompilation due to
!                        compiler wrapper optimation exception liat change

!----------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMCST   , ONLY : RETV     ,RLVTT    ,RLSTT    ,RTT
USE YOETHF   , ONLY : R2ES     ,R3LES    ,R3IES    ,R4LES    ,&
 & R4IES    ,R5LES    ,R5IES    ,R5ALVCP  ,R5ALSCP  ,&
 & RALVDCP  ,RALSDCP  ,RTWAT    ,RTICE    ,RTICECU  ,&
 & RTWAT_RTICE_R      ,RTWAT_RTICECU_R  

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KK 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSP(KLON) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PT(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PQ(KLON,KLEV) 
LOGICAL           ,INTENT(IN)    :: LDFLAG(KLON) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KCALL 
REAL(KIND=JPRB) ::     Z3ES(KLON),             Z4ES(KLON),&
 & Z5ALCP(KLON),           ZALDCP(KLON)  

INTEGER(KIND=JPIM) :: JL

REAL(KIND=JPRB) :: ZQMAX, ZQP, ZCOND, ZCOND1, ZTARG, ZCOR, ZQSAT, ZFOEEW, Z2S
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!DIR$ VFUNCTION EXPHF
#include "fcttre.func.h"
!----------------------------------------------------------------------

!     1.           DEFINE CONSTANTS
!                  ----------------

IF (LHOOK) CALL DR_HOOK('CUADJTQS',0,ZHOOK_HANDLE)
ZQMAX=0.5_JPRB

!     2.           CALCULATE CONDENSATION AND ADJUST T AND Q ACCORDINGLY
!                  -----------------------------------------------------

!*    ICE-WATER THERMODYNAMICAL FUNCTIONS

DO JL=KIDIA,KFDIA
  IF (PT(JL,KK) > RTT) THEN
    Z3ES(JL)=R3LES
    Z4ES(JL)=R4LES
    Z5ALCP(JL)=R5ALVCP
    ZALDCP(JL)=RALVDCP
  ELSE
    Z3ES(JL)=R3IES
    Z4ES(JL)=R4IES
    Z5ALCP(JL)=R5ALSCP
    ZALDCP(JL)=RALSDCP
  ENDIF
ENDDO

IF (KCALL == 1 ) THEN

!DIR$    IVDEP
!OCL NOVREC
  DO JL=KIDIA,KFDIA
    IF(LDFLAG(JL)) THEN
      ZQP    =1.0_JPRB/PSP(JL)
      ZTARG    =PT(JL,KK)
      ZFOEEW    =R2ES*EXP(Z3ES(JL)*(ZTARG    -RTT)/(ZTARG    -Z4ES(JL)))
      ZQSAT    =ZQP    *ZFOEEW    
      IF (ZQSAT     > ZQMAX) THEN
        ZQSAT    =ZQMAX
      ENDIF
      ZCOR    =1.0_JPRB/(1.0_JPRB-RETV*ZQSAT    )
      ZQSAT    =ZQSAT    *ZCOR    
      Z2S    =Z5ALCP(JL)/(ZTARG    -Z4ES(JL))**2
      ZCOND    =(PQ(JL,KK)-ZQSAT    )/(1.0_JPRB+ZQSAT    *ZCOR    *Z2S    )
      ZCOND    =MAX(ZCOND    ,0.0_JPRB)
!     IF(ZCOND /= _ZERO_) THEN
      PT(JL,KK)=PT(JL,KK)+ZALDCP(JL)*ZCOND    
      PQ(JL,KK)=PQ(JL,KK)-ZCOND    
      ZTARG    =PT(JL,KK)
      ZFOEEW    =R2ES*EXP(Z3ES(JL)*(ZTARG    -RTT)/(ZTARG    -Z4ES(JL)))
      ZQSAT    =ZQP    *ZFOEEW    
      IF (ZQSAT     > ZQMAX) THEN
        ZQSAT    =ZQMAX
      ENDIF
      ZCOR    =1.0_JPRB/(1.0_JPRB-RETV*ZQSAT    )
      ZQSAT    =ZQSAT    *ZCOR    
      Z2S    =Z5ALCP(JL)/(ZTARG    -Z4ES(JL))**2
      ZCOND1    =(PQ(JL,KK)-ZQSAT    )/(1.0_JPRB+ZQSAT    *ZCOR    *Z2S    )
      IF(ZCOND ==  0.0_JPRB)ZCOND1=0.0_JPRB
      PT(JL,KK)=PT(JL,KK)+ZALDCP(JL)*ZCOND1    
      PQ(JL,KK)=PQ(JL,KK)-ZCOND1    
!     ENDIF
    ENDIF
  ENDDO

ENDIF

IF(KCALL == 2) THEN

!DIR$    IVDEP
!OCL NOVREC
  DO JL=KIDIA,KFDIA
    IF(LDFLAG(JL)) THEN
      ZQP    =1.0_JPRB/PSP(JL)
      ZTARG    =PT(JL,KK)
      ZFOEEW    =R2ES*EXP(Z3ES(JL)*(ZTARG    -RTT)/(ZTARG    -Z4ES(JL)))
      ZQSAT    =ZQP    *ZFOEEW    
      IF (ZQSAT     > ZQMAX) THEN
        ZQSAT    =ZQMAX
      ENDIF
      ZCOR    =1.0_JPRB/(1.0_JPRB-RETV  *ZQSAT    )
      ZQSAT    =ZQSAT    *ZCOR    
      Z2S    =Z5ALCP(JL)/(ZTARG    -Z4ES(JL))**2
      ZCOND    =(PQ(JL,KK)-ZQSAT    )/(1.0_JPRB+ZQSAT    *ZCOR    *Z2S    )
      ZCOND    =MIN(ZCOND    ,0.0_JPRB)
!     IF(ZCOND /= _ZERO_) THEN
      PT(JL,KK)=PT(JL,KK)+ZALDCP(JL)*ZCOND    
      PQ(JL,KK)=PQ(JL,KK)-ZCOND    
      ZTARG    =PT(JL,KK)
      ZFOEEW    =R2ES*EXP(Z3ES(JL)*(ZTARG    -RTT)/(ZTARG    -Z4ES(JL)))
      ZQSAT    =ZQP    *ZFOEEW    
      IF (ZQSAT     > ZQMAX) THEN
        ZQSAT    =ZQMAX
      ENDIF
      ZCOR    =1.0_JPRB/(1.0_JPRB-RETV  *ZQSAT    )
      ZQSAT    =ZQSAT    *ZCOR    
      Z2S    =Z5ALCP(JL)/(ZTARG    -Z4ES(JL))**2
      ZCOND1    =(PQ(JL,KK)-ZQSAT    )/(1.0_JPRB+ZQSAT    *ZCOR    *Z2S    )
      IF(ZCOND ==  0.0_JPRB)ZCOND1=0.0_JPRB
      PT(JL,KK)=PT(JL,KK)+ZALDCP(JL)*ZCOND1    
      PQ(JL,KK)=PQ(JL,KK)-ZCOND1    
!     ENDIF
    ENDIF
  ENDDO

ENDIF

IF(KCALL == 0) THEN

!DIR$    IVDEP
!OCL NOVREC
  DO JL=KIDIA,KFDIA
    ZQP    =1.0_JPRB/PSP(JL)
    ZTARG    =PT(JL,KK)
    ZFOEEW    =R2ES*EXP(Z3ES(JL)*(ZTARG    -RTT)/(ZTARG    -Z4ES(JL)))
    ZQSAT    =ZQP    *ZFOEEW    
    IF (ZQSAT     > ZQMAX) THEN
      ZQSAT    =ZQMAX
    ENDIF
    ZCOR    =1.0_JPRB/(1.0_JPRB-RETV  *ZQSAT    )
    ZQSAT    =ZQSAT    *ZCOR    
    Z2S    =Z5ALCP(JL)/(ZTARG    -Z4ES(JL))**2
    ZCOND1    =(PQ(JL,KK)-ZQSAT    )/(1.0_JPRB+ZQSAT    *ZCOR    *Z2S    )
    PT(JL,KK)=PT(JL,KK)+ZALDCP(JL)*ZCOND1    
    PQ(JL,KK)=PQ(JL,KK)-ZCOND1    
    ZTARG    =PT(JL,KK)
    ZFOEEW    =R2ES*EXP(Z3ES(JL)*(ZTARG    -RTT)/(ZTARG    -Z4ES(JL)))
    ZQSAT    =ZQP    *ZFOEEW    
    IF (ZQSAT     > ZQMAX) THEN
      ZQSAT    =ZQMAX
    ENDIF
    ZCOR    =1.0_JPRB/(1.0_JPRB-RETV  *ZQSAT    )
    ZQSAT    =ZQSAT    *ZCOR    
    Z2S    =Z5ALCP(JL)/(ZTARG    -Z4ES(JL))**2
    ZCOND1    =(PQ(JL,KK)-ZQSAT    )/(1.0_JPRB+ZQSAT    *ZCOR    *Z2S    )
    PT(JL,KK)=PT(JL,KK)+ZALDCP(JL)*ZCOND1    
    PQ(JL,KK)=PQ(JL,KK)-ZCOND1    
  ENDDO

ENDIF

IF(KCALL == 4) THEN

!DIR$    IVDEP
!OCL NOVREC
  DO JL=KIDIA,KFDIA
    ZQP    =1.0_JPRB/PSP(JL)
    ZTARG    =PT(JL,KK)
    ZFOEEW    =R2ES*EXP(Z3ES(JL)*(ZTARG    -RTT)/(ZTARG    -Z4ES(JL)))
    ZQSAT    =ZQP    *ZFOEEW    
    IF (ZQSAT     > ZQMAX) THEN
      ZQSAT    =ZQMAX
    ENDIF
    ZCOR    =1.0_JPRB/(1.0_JPRB-RETV  *ZQSAT    )
    ZQSAT    =ZQSAT    *ZCOR    
    Z2S    =Z5ALCP(JL)/(ZTARG    -Z4ES(JL))**2
    ZCOND    =(PQ(JL,KK)-ZQSAT    )/(1.0_JPRB+ZQSAT    *ZCOR    *Z2S    )
    PT(JL,KK)=PT(JL,KK)+ZALDCP(JL)*ZCOND    
    PQ(JL,KK)=PQ(JL,KK)-ZCOND    
    ZTARG    =PT(JL,KK)
    ZFOEEW    =R2ES*EXP(Z3ES(JL)*(ZTARG    -RTT)/(ZTARG    -Z4ES(JL)))
    ZQSAT    =ZQP    *ZFOEEW    
    IF (ZQSAT     > ZQMAX) THEN
      ZQSAT    =ZQMAX
    ENDIF
    ZCOR    =1.0_JPRB/(1.0_JPRB-RETV  *ZQSAT    )
    ZQSAT    =ZQSAT    *ZCOR    
    Z2S    =Z5ALCP(JL)/(ZTARG    -Z4ES(JL))**2
    ZCOND1    =(PQ(JL,KK)-ZQSAT    )/(1.0_JPRB+ZQSAT    *ZCOR    *Z2S    )
    PT(JL,KK)=PT(JL,KK)+ZALDCP(JL)*ZCOND1    
    PQ(JL,KK)=PQ(JL,KK)-ZCOND1    
  ENDDO

ENDIF

IF (LHOOK) CALL DR_HOOK('CUADJTQS',1,ZHOOK_HANDLE)
END SUBROUTINE CUADJTQS
