! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

 SUBROUTINE CHEM_DECAY &
 &    (YGFL,KIDIA  , KFDIA , KLON, KLEV, PTSTEP, PRSF1, PGELAT, &
 &      PCEN ,PTENC) 

!**   DESCRIPTION 
!     ----------
!
!   DECAY routine for IFS chemistry 
!
!
!
!**   INTERFACE.
!     ----------
!          *CHEM_DECAY* IS CALLED FROM *CHEM_MASTER*.

! INPUTS:
! -------
! KIDIA :  Start of Array  
! KFDIA :  End  of Array 
! KLON  :  Length of Arrays 
! KLEV  :  Number of Levels
! PTSTEP:  Time step in seconds 
! PCEN(KLON,KLEV,NCHEM)        : CONCENTRATION OF TRACERS           (kg/kg)
! PRSF1(KLON,KLEV)            : FULL-LEVEL PRESSURE           (Pa)
! PGELAT(KLON)                : LATITUDE (RADIANS)

! OUTPUTS:
! -------
! PTENC  (KLON,KLEV,NCHEM)     : TENDENCY OF CONCENTRATION OF TRACERS BCEUASE OF DEACY(kg/kg s-1)
!
! LOCAL:
! -------
!
! ZCEN(KLON,KLEV,NCHEMNCHEM)        : CONCENTRATION OF TRACERS           (kg/kg)
! ZTENC(KLON,KLEV,NCHEM)       : TOTAL TENDENCY OF CONCENTRATION OF TRACERS BEFORE(kg/kg s-1)

!
!     AUTHOR.
!     -------
!        JOHANNES FLEMMING  *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 2009-07-22

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOM_YGFL , ONLY : TYPE_GFLD

IMPLICIT NONE

!-----------------------------------------------------------------------
!*       0.1  ARGUMENTS
!             ---------

TYPE(TYPE_GFLD)   ,INTENT(INOUT):: YGFL
INTEGER(KIND=JPIM),INTENT(IN) :: KIDIA , KFDIA , KLON , KLEV
REAL(KIND=JPRB),INTENT(IN)    :: PTSTEP
REAL(KIND=JPRB),INTENT(INOUT)   :: PTENC(KLON,KLEV, YGFL%NCHEM)
REAL(KIND=JPRB),INTENT(IN)    :: PCEN(KLON,KLEV, YGFL%NCHEM) 
REAL(KIND=JPRB),INTENT(IN)    :: PRSF1(KLON,KLEV) 
REAL(KIND=JPRB),INTENT(IN)    :: PGELAT(KLON) 


!
!*       0.5   LOCAL VARIABLES
!              ---------------

INTEGER(KIND=JPIM) :: JK, JL, JT, IO3ST
REAL(KIND=JPRB)    :: ZTDECAY
REAL(KIND=JPHOOK)  :: ZHOOK_HANDLE 


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('CHEM_DECAY',0,ZHOOK_HANDLE)
ASSOCIATE(NCHEM=>YGFL%NCHEM, YCHEM=>YGFL%YCHEM)
IO3ST = -1

!*       1.     CALCULATE TENDECIES  

! Loop over Chemical species 
DO JT=1,NCHEM 
  IF ( YCHEM(JT)%REFOLD > 0.0_JPRB) THEN 
! loss is ZTDECAY * Concentration 
   ZTDECAY=(EXP(-1.0_JPRB*PTSTEP/(YCHEM(JT)%REFOLD*86400.0_JPRB))-1.0_JPRB)/PTSTEP
! loop over levels
    DO JK=1,KLEV
! loop over points  
      DO JL=KIDIA,KFDIA
! calculate tendency  
        PTENC(JL,JK,JT)  =  PTENC(JL,JK,JT) + PCEN(JL,JK,JT)*ZTDECAY  
      ENDDO
    ENDDO
  ENDIF
ENDDO

!
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CHEM_DECAY',1,ZHOOK_HANDLE )
END SUBROUTINE CHEM_DECAY 

