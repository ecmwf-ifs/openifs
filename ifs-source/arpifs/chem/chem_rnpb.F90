! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

 SUBROUTINE CHEM_RNPB &
 &    (YGFL, KIDIA  , KFDIA , KLON , KLEV , KTRAC, &
 &     PTSTEP , PRSF1, PTP, PCEN ,&
 &     PTENC) 

!**   DESCRIPTION 
!     ----------
!
!   Radon (Rn) to Lead decay (Pb)  
!
!
!
!**   INTERFACE.
!     ----------
!          *CHEM_RNPB* IS CALLED FROM *CHEM_MAIN*.

! INPUTS:
! -------
!
! - DIMENSIONS ETC.
!
! KIDIA   :  Start of Array  
! KFDIA   :  End  of Array 
! KLON    :  Length of Arrays 
! KLEV    :  Number of Levels
! KTRAC :  Number tracers 
! PTSTEP  :  Time step in seconds 
!
! - 2D and 3D
!  PRSF1(KLON,KLEV)             : FULL-LEVEL PRESSURE           (Pa)
! PTP(KLON,KLEV)                : FULL-LEVEL TEMPERATURE (W. DYN.TEND.) (K)
! PCEN(KLON,KLEV,KTRAC)         : CONCENTRATION OF TRACERS           (kg/kg)
!
! OUTPUTS:
! -------
!
! PTENC  (KLON,KLEV,KTRAC)          : TENDENCY OF CONCENTRATION by chemistry(kg/kg s-1)
!
! LOCAL:
! -------
!
!     AUTHORS.
!     -------
!        JOHANNES FLEMMING  *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 2019-01-22


USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOM_YGFL  ,ONLY : TYPE_GFLD
USE YOMCST    ,ONLY : RMD
USE TM5_CHEM_MODULE, ONLY: IFLRN222, IFLPB210


IMPLICIT NONE

!-----------------------------------------------------------------------
!*       0.1  ARGUMENTS
!             ---------

TYPE(TYPE_GFLD)   ,INTENT(INOUT):: YGFL
INTEGER(KIND=JPIM),INTENT(IN) :: KIDIA , KFDIA , KLON , KLEV, KTRAC
REAL(KIND=JPRB),INTENT(IN)    :: PTSTEP
REAL(KIND=JPRB),INTENT(IN)    :: PTP(KLON,KLEV)   
REAL(KIND=JPRB),INTENT(IN)    :: PRSF1(KLON,KLEV)   
REAL(KIND=JPRB),INTENT(IN)    :: PCEN(KLON,KLEV,KTRAC) 
REAL(KIND=JPRB),INTENT(OUT)   :: PTENC(KLON,KLEV,KTRAC)

! * LOCAL 
INTEGER(KIND=JPIM) :: JL, JK
!
LOGICAL       :: LLPB 

REAL(KIND=JPRB)    :: ZRR, ZAIRDM, ZCEN1, ZCVM0RN222 

REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CHEM_RNPB',0,ZHOOK_HANDLE)

ASSOCIATE(YCHEM=>YGFL%YCHEM)

LLPB = ( IFLPB210 > 0_JPIM ) 

! see tm5_calcrates
ZRR=2.10E-6_JPRB

! Loop over the levels and the grid points
!
DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
! convert to mol/cm3
    ZAIRDM = (7.24291E16_JPRB*PRSF1(JL,JK)/PTP(JL,JK)) * RMD
    ZCVM0RN222  =  PCEN(JL,JK,IFLRN222) / YCHEM(IFLRN222)%RMOLMASS *ZAIRDM
    ZCEN1 = ( ZCVM0RN222 /(1._JPRB+ZRR*PTSTEP) - ZCVM0RN222 ) / PTSTEP
    PTENC(JL,JK,IFLRN222) =             ZCEN1 * YCHEM(IFLRN222)%RMOLMASS / ZAIRDM
    IF (LLPB) PTENC(JL,JK,IFLPB210) = - ZCEN1 * YCHEM(IFLPB210)%RMOLMASS / ZAIRDM
  ENDDO
ENDDO

END ASSOCIATE

IF (LHOOK) CALL DR_HOOK('CHEM_RNPB',1,ZHOOK_HANDLE)
END SUBROUTINE CHEM_RNPB
