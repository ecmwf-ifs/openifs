! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE CHEM_INEXT &
 & ( KIDIA , KFDIA  , KLON , KLEV , KCHEM, KCEXTR, PDP, PTSTEP,  PTENC1, PTENC0, PEXTRA, KTOP )

!*** * CHEM_INEXTR* - ifill extra arrays with total column fields of tracer
!**** Levels are used to store different tracers 
! INPUTS:
! -------
! KIDIA :  Start of Array
! KFDIA :  End  of Array
! KLON  :  Length of Arrays
! KLEV  :  Number of Levels
! KCHEM :  Dimension of 
! KCEXTR : 2nd Dimension (levels) of extra fields    
! PTSTEP : time step 
! PDP(KLON,KLEV)              :  PRESSURE DELTA in PRESSURE UNITES      (Pa)
! PTENC0(KLON,KLEV,KCHEM)     :  TENDENCY OF CONCENTRATION TO BEFORE PROCESS (kg/kg s-1)
! PTENC1(KLON,KLEV,KCHEM)     :  TENDENCY OF CONCENTRATION TO AFTER PROCESS (kg/kg s-1)
! KTOP : Top level for TC
! OUTPUT:
! -------
! PEXTRA (KLON,KCEXTR)         : TC of PTENC1 - PTENC0 (kg/kg s-1)
!
!**   INTERFACE.
!     ----------
!          *CHEM_I* IS CALLED FROM *CHEM_MAIN*.
!
!     AUTHOR.
!     -------
!        Johannes Flemming 
!        
!
!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 2011-03-31

!-----------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOMCST , ONLY :  RG

IMPLICIT NONE

#include "abor1.intfb.h"
!-----------------------------------------------------------------------

!*       0.1  ARGUMENTS
!             ---------

INTEGER(KIND=JPIM),INTENT(IN) :: KIDIA, KFDIA, KLON, KLEV, KCHEM, KCEXTR 
INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: KTOP(KLON)
REAL(KIND=JPRB),INTENT(IN)    :: PDP(KLON,KLEV)    
REAL(KIND=JPRB),INTENT(IN)    :: PTSTEP   
REAL(KIND=JPRB),INTENT(IN)    :: PTENC0(KLON,KLEV,KCHEM) 
REAL(KIND=JPRB),INTENT(IN)    :: PTENC1(KLON,KLEV,KCHEM) 
REAL(KIND=JPRB),INTENT(INOUT)   :: PEXTRA(KLON,KCEXTR)

!*       0.5   LOCAL VARIABLES
!              ---------------


INTEGER(KIND=JPIM) :: JK, JL, JT, ITOP(KLON)
REAL(KIND=JPRB) :: ZTC(KLON) 
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('CHEM_INEXT',0,ZHOOK_HANDLE)

ITOP(:)=1
IF (PRESENT(KTOP)) ITOP=KTOP
IF (KCHEM > KCEXTR )  CALL ABOR1(" NOT ENOUGH LEVELS IN EXTRA FIELDS IN CHEM_INEXT ")

DO JT=1,KCHEM
  ZTC(:) = 0.0_JPRB 
  DO JL=KIDIA,KFDIA
     DO JK=ITOP(JL),KLEV
       ZTC (JL) =  ZTC(JL) + ( PTENC1(JL,JK,JT) -  PTENC0(JL,JK,JT) )* PDP(JL,JK) / RG
!       IF (ABS( PTENC1(JL,JK,JT) -  PTENC0(JL,JK,JT) ) > 0.1 ) THEN
!         print*, '  STRANGE TEND  ', JL, JK, JT
!       ENDIF
     ENDDO
  ENDDO
  PEXTRA(KIDIA:KFDIA,JT) = PEXTRA(KIDIA:KFDIA,JT) + ZTC(KIDIA:KFDIA) * PTSTEP
ENDDO

!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('CHEM_INEXT',1,ZHOOK_HANDLE)
END SUBROUTINE CHEM_INEXT 
