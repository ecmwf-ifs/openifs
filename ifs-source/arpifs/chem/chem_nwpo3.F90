! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

 SUBROUTINE CHEM_NWPO3 &
 &    (YDMODEL, KIDIA, KFDIA , KLON , KLEV, KTRAC, KVCLIS, &
 &     PTSTEP, PKOZO, PCSZA, PGEMU, PDELP,  PRS1, PRSF1, PTP, PCEN, &
 &     PTENC) 

!**   DESCRIPTION 
!     ----------
!
!   Call Cariolle chemistry in O3CHEM 
!
!
!
!**   INTERFACE.
!     ----------
!          *CHEM_NWPO3* IS CALLED FROM *CHEM_MAIN*.

! INPUTS:
! -------
!
! - DIMENSIONS ETC.
!
! KIDIA   :  Start of Array  
! KFDIA   :  End  of Array 
! KLON    :  Length of Arrays 
! KLEV    :  Number of Levels
! KVCLIS                      : Number Cariolle chemistry coefficinets 
! PTSTEP  :  Time step in seconds 
!
! - 2D and 3D
! PDELP(KLON,KLEV)            : PRESSURE DELTA in PRESSURE UNITES      (Pa)
! PRS1(KLON,0:KLEV)           : HALF-LEVEL PRESSURE           (Pa)
! PRSF1(KLON,KLEV)            : FULL-LEVEL PRESSURE           (Pa)
! PTP(KLON,KLEV)              : FULL-LEVEL TEMPERATURE (W. DYN.TEND.) (K)
! PKOZO(KLON,KLEV,KVCLIS)     : PHOTOCHEMICAL COEFFICIENTS COMPUTED FROM A 2D PHOTOCHEMICAL MODEL (KVCLIS=8)!
! PCSZA(KLON)                  : COS of Solar Zenit Angle
! PGEMU(KLON)                 : SINE OF LATITUDE
! PCEN(KLON,KLEV,1)            : CONCENTRATION OF TRACERS           (kg/kg)
!
! OUTPUTS:
! -------
!
! PTENC  (KLON,KLEV,KTRAC)          : TENDENCY OF CONCENTRATION OF TRACERS including chemistry(kg/kg s-1)
!
! LOCAL:
! -------
!
! ZKCO(KLON,KLEV,NCHEM_LCOCOEF) : COEFICIENTS OF THE LINEAR SCHEME (volume mixing ratio)
! ZCEN(KLON,KLEV,1)             : CONCENTRATION OF TRACERS           (kg/kg)
!
!     AUTHORS.
!     -------
!        JOHANNES FLEMMING  *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 2019-07-22

USE TYPE_MODEL , ONLY : MODEL
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE TM5_CHEM_MODULE, ONLY: IFLO3
 

IMPLICIT NONE

!-----------------------------------------------------------------------
!*       0.1  ARGUMENTS
!             ---------
TYPE(MODEL)       ,INTENT(INOUT):: YDMODEL
INTEGER(KIND=JPIM),INTENT(IN) :: KIDIA , KFDIA , KLON , KLEV, KVCLIS, KTRAC
REAL(KIND=JPRB),INTENT(IN)    :: PTSTEP
REAL(KIND=JPRB),INTENT(IN)    :: PTP(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN) :: PRSF1(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN) :: PDELP(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN) :: PRS1(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(IN) :: PKOZO(KLON,KLEV,KVCLIS)
REAL(KIND=JPRB)   ,INTENT(IN) :: PCSZA(KLON)                  
REAL(KIND=JPRB)   ,INTENT(IN) :: PGEMU(KLON)                  
REAL(KIND=JPRB),INTENT(IN)    :: PCEN(KLON,KLEV,KTRAC) 
REAL(KIND=JPRB),INTENT(OUT)   :: PTENC(KLON,KLEV,KTRAC)

! * LOCAL 

REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
#include "o3chem.intfb.h"

IF (LHOOK) CALL DR_HOOK('CHEM_NWPO3',0,ZHOOK_HANDLE)

 CALL O3CHEM (YDMODEL%YRML_GCONF%YRRIP,YDMODEL%YRML_PHY_EC%YREPHY,YDMODEL%YRML_CHEM%YROZO,YDMODEL%YRML_PHY_MF%YRPHY2, &
    & KIDIA, KFDIA, KLON, 1, KLEV, KVCLIS, PGEMU, PCSZA,&
    & PRS1, PRSF1, PKOZO, PDELP, PTP, PCEN(:,:,IFLO3),&
    & PTENC(:,:,IFLO3)) 

IF (LHOOK) CALL DR_HOOK('CHEM_NWPO3',1,ZHOOK_HANDLE)
END SUBROUTINE CHEM_NWPO3
