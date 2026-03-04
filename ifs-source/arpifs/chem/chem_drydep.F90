! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE CHEM_DRYDEP &
 & ( YGFL,KIDIA , KFDIA  , KLON , KLEV , KCHEM , KCHEM_DRYDEP, PTSTEP, &
 &    PGEOH, PDELP , PCHEMDV, PCEN, PTENC0,  PTENC1, PWLOSS )

!*** * CHEM_DRYDEP* - dry deposition separately
! INPUTS:
! -------
! KIDIA :  Start of Array
! KFDIA :  End  of Array
! KLON  :  Length of Arrays
! KLEV  :  Number of Levels
! KCHEM :  Number of chemistry tracers 
! KCHEM_DRYDEP :  Number of chemistry tracers with dry depostion


! PTSEP                       :  Time step length in seconds
! PDP(KLON,KLEV)              :  PRESSURE DELTA in PRESSURE UNITES      (Pa)
! PGEOH(KLON,KLEV+1)          :  Geopotential over model surface m**2/s**2
! PCHEMDV(KLON,KCHEM_DRYDEP)  :  dry deposition Velocity m/s
! PCEN(KLON,KLEV,KCHEM)       :  CONCENTRATION OF TRACERS           (kg/kg)
! PTENC0(KLON,KLEV,KCHEM)     :  TOTAL TENDENCY OF CONCENTRATION OF TRACERS BEFORE(kg/kg s-1)
!
! -------
! PTENC1 (KLON,KLEV,KCHEM)     : TENDENCY OF CONCENTRATION OF TRACERS after (kg/kg s-1)
! PLOSS(KLON,KCHEM_DRYDEP)          : Total Mass Loss due to scavening     (kg/m2 s-1)!
!
!**   INTERFACE.
!     ----------
!          *CHEM_DRYDEP* IS CALLED FROM *CHEM_MAIN*.
!
!     AUTHOR.
!     -------
!        Johannes Flemming 
!        
!
!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 2010-11-09

!-----------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST   , ONLY : RG 
USE YOMLUN   , ONLY : NULOUT
USE YOM_YGFL , ONLY : TYPE_GFLD

IMPLICIT NONE

!-----------------------------------------------------------------------

!*       0.1  ARGUMENTS
!             ---------

TYPE(TYPE_GFLD)   ,INTENT(INOUT):: YGFL
INTEGER(KIND=JPIM),INTENT(IN) :: KIDIA, KFDIA, KLON, KLEV, KCHEM, KCHEM_DRYDEP
REAL(KIND=JPRB),INTENT(IN)    :: PGEOH(KLON,0:KLEV)
REAL(KIND=JPRB),INTENT(IN)    :: PDELP(KLON,KLEV)    
REAL(KIND=JPRB),INTENT(IN)    :: PTENC0(KLON,KLEV,KCHEM), PCEN(KLON,KLEV,KCHEM), PCHEMDV(KLON,KCHEM_DRYDEP)
REAL(KIND=JPRB),INTENT(IN)    :: PTSTEP

REAL(KIND=JPRB),INTENT(OUT)   :: PWLOSS(KLON,KCHEM_DRYDEP), PTENC1(KLON,KLEV,KCHEM)

!*       0.5   LOCAL VARIABLES
!              ---------------

INTEGER(KIND=JPIM) :: JL, JT, IDRYDEP

REAL(KIND=JPRB) :: ZALPHA, ZRRG, ZHGT, ZRTSTEP 
REAL(KIND=JPRB) :: ZCEN1(KLON)
CHARACTER(LEN=3) :: CLNUM 
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('CHEM_DRYDEP',0,ZHOOK_HANDLE)
ASSOCIATE(YCHEM=>YGFL%YCHEM)
! Euler forward
! CLNUM='FWD'
! Euler backward
!CLNUM='BWD'
! exponential
! CLNUM='EXP'
! centrred 
! CLNUM='CTR'

! use analytical solution
CLNUM='EXP'

!* update tendecies  
PTENC1(KIDIA:KFDIA,1:KLEV,1:KCHEM) =PTENC0(KIDIA:KFDIA,1:KLEV,1:KCHEM)
PWLOSS(:,:)=0.0_JPRB

IDRYDEP=0

ZRRG = 1.0 / RG
ZRTSTEP = 1.0/ PTSTEP 
DO JT=1,KCHEM
 IF ( YCHEM(JT)%IGRIBDV <=  0 ) CYCLE
    IDRYDEP = IDRYDEP+1

     DO JL=KIDIA,KFDIA
        ZHGT= PGEOH(JL,KLEV-1) * ZRRG 
         IF (PCHEMDV(JL,IDRYDEP) < 0.0_JPRB .OR. ZHGT < 0.0_JPRB ) THEN
          WRITE (NULOUT,*) JL, JT, IDRYDEP,PCHEMDV(JL,IDRYDEP),ZHGT 
          CALL ABOR1( " sign error drydep")
         ENDIF  
     ENDDO
    IF (CLNUM == 'FWD' ) THEN
      DO JL=KIDIA,KFDIA
         ZHGT= PGEOH(JL,KLEV-1) * ZRRG 
         ZALPHA=PTSTEP*PCHEMDV(JL,IDRYDEP)/ZHGT 
         ZCEN1(JL)=PCEN(JL,KLEV,JT) * ( 1.0_JPRB -  ZALPHA) 
       ENDDO
    ENDIF
    IF (CLNUM == 'BWD' ) THEN
      DO JL=KIDIA,KFDIA
         ZHGT = PGEOH(JL,KLEV-1) * ZRRG 
         ZALPHA=PTSTEP*PCHEMDV(JL,IDRYDEP)/ZHGT 
         ZCEN1(JL)=PCEN(JL,KLEV,JT) / ( 1.0_JPRB + ZALPHA) 
       ENDDO
    ENDIF
    IF (CLNUM == 'CTR' ) THEN
      DO JL=KIDIA,KFDIA
         ZHGT = PGEOH(JL,KLEV-1) * ZRRG 
         ZALPHA=PTSTEP*PCHEMDV(JL,IDRYDEP)/ZHGT 
         ZCEN1(JL)=PCEN(JL,KLEV,JT) * (( 1.0_JPRB - ZALPHA) / ( 1.0_JPRB + ZALPHA)) 
       ENDDO
    ENDIF
    IF (CLNUM == 'EXP' ) THEN
      DO JL=KIDIA,KFDIA
         ZHGT = PGEOH(JL,KLEV-1) * ZRRG 
         ZALPHA=PTSTEP*PCHEMDV(JL,IDRYDEP)/ZHGT 
         ZCEN1(JL)=PCEN(JL,KLEV,JT) * EXP( -1.0_JPRB *  ZALPHA)  
       ENDDO
    ENDIF 
    DO JL=KIDIA,KFDIA
      PTENC1(JL,KLEV,JT) = PTENC0(JL,KLEV,JT) + ( ZCEN1(JL)-PCEN(JL, KLEV, JT) ) * ZRTSTEP 
    ENDDO
     
    DO JL=KIDIA,KFDIA  
     PWLOSS(JL,IDRYDEP)=  (( PCEN(JL, KLEV, JT) -  ZCEN1(JL)) * ZRTSTEP ) * PDELP(JL,KLEV)*ZRRG 
    ENDDO
ENDDO

!!!!! debug undo anything  
!!!!! PTENC1(KIDIA:KFDIA,1:KLEV,1:KCHEM) =PTENC0(KIDIA:KFDIA,1:KLEV,1:KCHEM)
!-----------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CHEM_DRYDEP',1,ZHOOK_HANDLE)
END SUBROUTINE CHEM_DRYDEP 
