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

SUBROUTINE LAIDLI_RAD(KPROMA,KFIELDS,PDLAT,PDLO,KL0,PRE,POST)

!**** *LAIDLI_RAD- Observation interpolator: 
!                 Bilinear horizontal interpolations for one variable.

!     Purpose.
!     --------
!       Performs bilinear horizontal interpolations for one variable.
!       Version designed for observation interpolator.
!       Differences with LAIDLI are the following ones:
!        - the position of the interpolation point and the weights
!          are considered as constants.

!     Author.
!     -------
!       Mats Hamrud - copy of laidliobs before it changed

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
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDLAT     ! Latitude weights
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDLO(1:2) ! Longitude weights
INTEGER(KIND=JPIM),INTENT(IN)    :: KL0(2,2)  ! Index of four surrounding grid points

REAL(KIND=JPRB)   ,INTENT(IN)    :: PRE(KPROMA,KFIELDS) ! Fields on grid points
REAL(KIND=JPRB)   ,INTENT(OUT)   :: POST(KFIELDS)       ! Fields in osbervation space

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: IP11, IP12, IP21, IP22, JF

REAL(KIND=JPRB) :: ZDLAT, ZDLO1, ZDLO2, ZVALLO1, ZVALLO2
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LAIDLI_RAD',0,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

!*       1.    INTERPOLATIONS.
!              ---------------

  !     Computation of coordinates and distances.

ZDLAT =PDLAT

ZDLO1 =PDLO(1)
ZDLO2 =PDLO(2)

IP11 = KL0(1,1)
IP12 = KL0(1,2)
IP21 = KL0(2,1)
IP22 = KL0(2,2)

  !     Interpolation.
  ! Yes, there is a stride in this copy, but that needs to happen somewhere
  ! Here's as good as any.
DO JF=1,KFIELDS
  ZVALLO1=PRE(IP11,JF) + ZDLO1*( PRE(IP12,JF)-PRE(IP11,JF) )
  ZVALLO2=PRE(IP21,JF) + ZDLO2*( PRE(IP22,JF)-PRE(IP21,JF) )

  POST(JF) = ZVALLO1 + ZDLAT*(ZVALLO2-ZVALLO1)
ENDDO

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('LAIDLI_RAD',1,ZHOOK_HANDLE)
END SUBROUTINE LAIDLI_RAD
