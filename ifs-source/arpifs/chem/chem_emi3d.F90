! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE CHEM_EMI3D &
 & ( KIDIA , KFDIA  , KLON , KLEV ,  PDP , PEMI3D,  PTENC0, PTENC1 )

!*** * CHEM_EMI3D* - ADD NO Air craft emissions 
! INPUTS:
! -------
! KIDIA :  Start of Array
! KFDIA :  End  of Array
! KLON  :  Length of Arrays
! KLEV  :  Number of Levels


! PDP(KLON,KLEV)              :  PRESSURE DELTA in PRESSURE UNITES      (Pa)
! PEMI3D(KLON,KLEV)      :  NO Aircraft emissions in kg/ms per level 
! PTENC0(KLON,KLEV)     :  TOTAL TENDENCY OF CONCENTRATION OF TRACERS BEFORE(kg/kg s-1)
!
! NB: PCLWAT is the in-cloud water mixing ratio
! OUTPUTS:
! -------
! PTENC1 (KLON,KLEV)     : TENDENCY OF CONCENTRATION OF TRACERS after (kg/kg s-1)
!
!**   INTERFACE.
!     ----------
!          *CHEM_EMI3D* IS CALLED FROM *CHEM_MAIN*.
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

!-----------------------------------------------------------------------

!*       0.1  ARGUMENTS
!             ---------

INTEGER(KIND=JPIM),INTENT(IN) :: KIDIA, KFDIA, KLON, KLEV 

REAL(KIND=JPRB),INTENT(IN)    :: PDP(KLON,KLEV)    
REAL(KIND=JPRB),INTENT(IN)    :: PEMI3D(KLON,KLEV)
REAL(KIND=JPRB),INTENT(IN)    :: PTENC0(KLON,KLEV) 

REAL(KIND=JPRB),INTENT(OUT)   :: PTENC1(KLON,KLEV)

!*       0.5   LOCAL VARIABLES
!              ---------------

INTEGER(KIND=JPIM) :: JK, JL

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('CHEM_EMI3D',0,ZHOOK_HANDLE)

PTENC1(KIDIA:KFDIA,1:KLEV) =PTENC0(KIDIA:KFDIA,1:KLEV)


!* LOOP OVER LAYERS T 
    DO JK=1,KLEV
      DO JL=KIDIA,KFDIA
      PTENC1(JL,JK) =  PTENC0(JL,JK) + PEMI3D(JL,JK) * ( RG / PDP(JL,JK) )
     ENDDO
  ENDDO

!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('CHEM_EMI3D',1,ZHOOK_HANDLE)
END SUBROUTINE CHEM_EMI3D 
