! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction
! 
! (C) Copyright 1989- Meteo-France.
! 

SUBROUTINE LAIDDI_RAD(KPROMA,KFIELDS,PCLA,PDLO,PCLO,KL0,PRE,POST)  

!**** *LAIDDI_RAD-  Observation interpolator:
!                  BI-dimensional 12 points interpolations (adjoint)

!     Purpose.
!     --------
!       Performs bi-dimensional 12 points interpolations.
!       Version designed for observation interpolator.
!       Differences with LAIDDI are the following ones:
!        - the position of the interpolation point and the weights
!          are considered as constants.

!     Author.
!     -------
!       Mats Hamrud - copy of laiddiobs before it changed

!     Modifications.
!     --------------
!        Original : May 2018.

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA    ! Number of grid (model) points
INTEGER(KIND=JPIM),INTENT(IN)    :: KFIELDS   ! Number fields to be interpolated
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCLA(3)   ! Latitude interpolation weights
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDLO(0:3) ! Longitude interpolation weights
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCLO(3,2) ! Longitude interpolation weights
INTEGER(KIND=JPIM),INTENT(IN)    :: KL0(0:3,0:3)  ! Indices of surrounding grid points
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRE(KPROMA,KFIELDS) ! Fields on grid points
REAL(KIND=JPRB)   ,INTENT(OUT)   :: POST(KFIELDS)       ! Fields in osbervation space

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: IP01, IP02, IP10, IP11, IP12, IP13, IP20,&
 & IP21, IP22, IP23, IP31, IP32 

REAL(KIND=JPRB), DIMENSION(KFIELDS) :: Z0, Z1, Z2, Z3

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LAIDDI_RAD',0,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

!*       1.    INTERPOLATIONS.
!              ---------------
IP31 = KL0(3,1)
IP32 = KL0(3,2)
IP20 = KL0(2,0)
IP21 = KL0(2,1)
IP22 = KL0(2,2)
IP23 = KL0(2,3)
IP10 = KL0(1,0)
IP11 = KL0(1,1)
IP12 = KL0(1,2)
IP13 = KL0(1,3)
IP01 = KL0(0,1)
IP02 = KL0(0,2)

 ! interpolations in longitude
Z0(:)=PRE(IP01,:) + PDLO(0) * ( PRE(IP02,:)-PRE(IP01,:) )
Z1(:)=PRE(IP10,:)&
  & +PCLO(1,1)*( PRE(IP11,:)-PRE(IP10,:) )&
  & +PCLO(2,1)*( PRE(IP12,:)-PRE(IP10,:) )&
  & +PCLO(3,1)*( PRE(IP13,:)-PRE(IP10,:) )  
Z2(:)=PRE(IP20,:)&
  & +PCLO(1,2)*( PRE(IP21,:)-PRE(IP20,:) )&
  & +PCLO(2,2)*( PRE(IP22,:)-PRE(IP20,:) )&
  & +PCLO(3,2)*( PRE(IP23,:)-PRE(IP20,:) )  
Z3(:)=PRE(IP31,:) + PDLO(3) * ( PRE(IP32,:)-PRE(IP31,:) )

! final interpolation in latitude
POST(:)= Z0 +PCLA(1)*(Z1-Z0) +PCLA(2)*(Z2-Z0) +PCLA(3)*(Z3-Z0)

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LAIDDI_RAD',1,ZHOOK_HANDLE)
END SUBROUTINE LAIDDI_RAD
