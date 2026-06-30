MODULE CANALB_MOD
CONTAINS
SUBROUTINE CANALB(KIDIA,KFDIA,YDURB,PMU0,PCANALB)

    ! (C) Copyright 2010- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!     PURPOSE
!     -------
!     THIS ROUTINE CALCULATES LOCAL URBAN ALBEDO AS A FUNCTION OF SOLAR ZENITH ANGLE

!     INTERFACE.
!     ----------
!     CALLLED FROM *Surfrad*

!     METHOD.
!     canyon albedo is calculated at each timestep to account for 
!     solar zenith angle

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     Albedo scheme based Masson 2000 shadowing, Harman 2004 matrix exchanges 
!     and Porson 2010 formulation

!     MODIFICATIONS
!     -------------
!==============================================================================

! Modules used : 
USE PARKIND1  , ONLY : JPIM, JPRB
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOS_URB   , ONLY : TURB

IMPLICIT NONE

REAL(KIND=JPRB)   , INTENT(IN)   :: PMU0(:)  ! Cosine of solar zenith angle
TYPE(TURB)        , INTENT(IN)   :: YDURB   
REAL(KIND=JPRB)   , INTENT(OUT)  :: PCANALB(:) ! Canyon Albedo (Not urban albedo)

! Local Arrays


INTEGER(KIND=JPIM),INTENT(IN)   :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN)   :: KFDIA
INTEGER(KIND=JPIM)              :: JL, i, j

REAL(KIND=JPRB) :: x_f      ! Fraction SW scattered by sky
REAL(KIND=JPRB) :: x_r      ! Fraction SW scattered by road
REAL(KIND=JPRB) :: x_w      ! Fraction SW scattered by wall
REAL(KIND=JPRB) :: omega0   ! F(solar zenith angle, aspect ratio)
REAL(KIND=JPRB) :: pi
REAL(KIND=JPRB) :: var1     ! Temporary variable
REAL(KIND=JPRB) :: aa       ! Temporary variable
REAL(KIND=JPRB) :: bb       ! Temporary variable
REAL(KIND=JPRB) :: cc       ! Temporary variable
REAL(KIND=JPRB) :: dd       ! Temporary variable

INTEGER, PARAMETER:: matdim = 3
REAL (KIND = JPRB):: shf(matdim,matdim)     ! shape factor
REAL (KIND = JPRB):: delta(matdim,matdim)   ! diagonal
REAL (KIND = JPRB):: rfl(matdim,matdim)     ! reflection matrix
REAL (KIND = JPRB):: rinv(matdim,matdim)    ! inverse matrix
REAL (KIND = JPRB):: detinv
REAL (KIND = JPRB):: ref(matdim)            ! reference fluxes

REAL (KIND = JPRB):: out(matdim)            ! outgoing radiation
REAL (KIND = JPRB):: inc(matdim)            ! incoming radiation

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('CANALB_MOD:CANALB',0,ZHOOK_HANDLE)
ASSOCIATE(RWALALB=>YDURB%RWALALB,RROAALB=>YDURB%RROAALB,RHWR=>YDURB%RHWR)

DO JL=KIDIA,KFDIA


 x_f = 0.3_JPRB
 pi  = 4.0_JPRB*ATAN(1.0_JPRB)
 
! Masson 2000 equation 13 & 14
 var1   = tan(acos(PMU0(JL)))

! limit check
 IF(var1 == 0.0_JPRB) THEN
 var1 = 0.0001_JPRB
 ENDIF
 omega0 = 1.0_JPRB/(var1*RHWR)
 var1 = min(omega0,1.0_JPRB)
 omega0 = asin(var1)
 x_r = (2.0_JPRB/pi)*(omega0 - RHWR*var1*(1.0_JPRB-cos(omega0)))
 x_w = (1.0_JPRB-x_r)/(2.0_JPRB*RHWR)

! Harman 2004 equations 2-5
! SET UP ARRAY
 aa = (1.0_JPRB+RHWR**2.0_JPRB)**(0.5_JPRB) - RHWR
 bb = (1.0_JPRB+(1.0_JPRB/RHWR)**2.0_JPRB)**(0.5_JPRB) - (1.0_JPRB/RHWR)
 cc = (1.0_JPRB-aa)
 dd = (1.0_JPRB-bb)/2.0_JPRB
! SHAPE FACTORS Porson equations 12-14
 shf(1,1) = 0.0_JPRB
 shf(1,2) = cc
 shf(1,3) = aa
 shf(2,1) = dd
 shf(2,2) = bb
 shf(2,3) = dd
 shf(3,1) = aa
 shf(3,2) = cc
 shf(3,3) = 0.0_JPRB

! REFLECTIONS
 rfl(1,1) = 1.0_JPRB
 rfl(1,2) = -RROAALB*cc
 rfl(1,3) = -RROAALB*aa
 rfl(2,1) = -RWALALB*dd
 rfl(2,2) = 1.0_JPRB-RWALALB*bb
 rfl(2,3) = -RWALALB*dd
 rfl(3,1) = 0.0_JPRB
 rfl(3,2) = 0.0_JPRB
 rfl(3,3) = 1.0_JPRB

    
! INVERT REFLECTION MATRIX
 detinv = 1/(rfl(1,1)*rfl(2,2)*rfl(3,3) - rfl(1,1)*rfl(2,3)*rfl(3,2)&
          - rfl(1,2)*rfl(2,1)*rfl(3,3) + rfl(1,2)*rfl(2,3)*rfl(3,1)&
          + rfl(1,3)*rfl(2,1)*rfl(3,2) - rfl(1,3)*rfl(2,2)*rfl(3,1))

! Calculate the inverse of the matrix
 rinv(1,1) = +detinv * (rfl(2,2)*rfl(3,3) - rfl(2,3)*rfl(3,2))
 rinv(2,1) = -detinv * (rfl(2,1)*rfl(3,3) - rfl(2,3)*rfl(3,1))
 rinv(3,1) = +detinv * (rfl(2,1)*rfl(3,2) - rfl(2,2)*rfl(3,1))
 rinv(1,2) = -detinv * (rfl(1,2)*rfl(3,3) - rfl(1,3)*rfl(3,2))
 rinv(2,2) = +detinv * (rfl(1,1)*rfl(3,3) - rfl(1,3)*rfl(3,1))
 rinv(3,2) = -detinv * (rfl(1,1)*rfl(3,2) - rfl(1,2)*rfl(3,1))
 rinv(1,3) = +detinv * (rfl(1,2)*rfl(2,3) - rfl(1,3)*rfl(2,2))
 rinv(2,3) = -detinv * (rfl(1,1)*rfl(2,3) - rfl(1,3)*rfl(2,1))
 rinv(3,3) = +detinv * (rfl(1,1)*rfl(2,2) - rfl(1,2)*rfl(2,1))

! Porson 2010 equation 10
    ref(1) = RROAALB*(1.0_JPRB-x_f)*x_r
    ref(2) = RWALALB*(1.0_JPRB-x_f)*x_w
    ref(3) = x_f

! TOT_OUT
 DO i=1,3
  OUT(i) = 0.0_JPRB
  DO j=1,3
   OUT(i) = OUT(i) + rinv(i,j)*ref(j)
  ENDDO
 ENDDO
! TOT_IN
  DO i=1,3
  INC(i) = 0.0_JPRB
  DO j=1,3
   INC(i) = INC(i) + shf(i,j)*OUT(j)
  ENDDO
 ENDDO

    inc(1)      = inc(1) + (1.0_JPRB-x_f)*x_r - OUT(1)
    inc(2)      = inc(2) + (1.0_JPRB-x_f)*x_w - OUT(2)
    PCANALB(JL) = inc(1) + 2.0_JPRB*RHWR*inc(2)

 PCANALB(JL) = 1.0_JPRB - PCANALB(JL) 

 ! limit check
 IF(PCANALB(JL) > 1.0_JPRB) THEN
 PCANALB(JL) = 1.0_JPRB
 ENDIF
 
 IF(PCANALB(JL) < 0.0_JPRB) THEN
 PCANALB(JL) = 0.0_JPRB
 ENDIF
 
ENDDO

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CANALB_MOD:CANALB',1,ZHOOK_HANDLE)

END SUBROUTINE CANALB
END MODULE CANALB_MOD
