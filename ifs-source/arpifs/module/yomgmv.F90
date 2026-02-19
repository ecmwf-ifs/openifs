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

MODULE YOMGMV

! Module containing t,t-dt,t+dt gridpoint arrays (apart from GFL) for dynamics
! and all "pointers" for accesssing elements in the arrays
! GMV and GMVS are permanently allocated, GMVT1 and GMVT1S temporary

!-------------------------------------------------------------------------

USE PARKIND1,  ONLY : JPIM, JPRB
USE TYPE_GMVS, ONLY : TYPE_T0, TYPE_T9, TYPE_T1, TYPE_PH9, TYPE_GP
USE YOMLUN,    ONLY : NULERR
IMPLICIT NONE
SAVE

TYPE :: TGMV

!-------------------------------------------------------------------------

REAL(KIND=JPRB), ALLOCATABLE :: GMV   (:,:,:,:) ! Multilevel fields at t and t-dt
REAL(KIND=JPRB), ALLOCATABLE :: GMVS  (:,:,:)   ! Single level fields at t snd t-dt
REAL(KIND=JPRB), ALLOCATABLE :: GMVT1 (:,:,:,:) ! Multilevel fields at t+dt
REAL(KIND=JPRB), ALLOCATABLE :: GMVT1S(:,:,:)   ! Single level fields at t+dt 

REAL(KIND=JPRB), ALLOCATABLE :: GMV5  (:,:,:,:)       ! Multilevel fields trajectory
REAL(KIND=JPRB), ALLOCATABLE :: GMV5S (:,:,:)         ! Single level fields trajectory 
! The following two fields are for 3-D FGAT (LIDMODEL)
REAL(KIND=JPRB), ALLOCATABLE :: GMV_DEPART  (:,:,:,:) ! Multilevel fields departure
REAL(KIND=JPRB), ALLOCATABLE :: GMVS_DEPART (:,:,:)   ! Single level fields departure

INTEGER(KIND=JPIM) :: NDIMGMV ! Third dim. of GMV "(NPROMA,NFLEVG,NDIMGMV,NGPBLKS)"
INTEGER(KIND=JPIM) :: NDIMGMVS ! Second dim. GMVS "(NPROMA,NDIMGMVS,NGPBLKS)"

!-------------------------------------------------------------------------

TYPE(TYPE_T0)  :: YT0  ! Pointers to time t quantities
TYPE(TYPE_T9)  :: YT9  ! Pointers to time t-dt quantities
TYPE(TYPE_T1)  :: YT1  ! Pointers to time t+dt quantities
TYPE(TYPE_PH9) :: YPH9 ! Pointers to physics time t-dt quantities
TYPE(TYPE_GP)  :: YGP  ! Pointers (offsets) for individual fields in GMV grid-point arrays

TYPE(TYPE_T0)  :: YT5  ! Pointers to trajectory quantities

!-------------------------------------------------------------------------

END TYPE TGMV

!-------------------------------------------------------------------------

CONTAINS

SUBROUTINE ZERO_YOMGMV(SELF)

USE YOMHOOK,   ONLY : LHOOK, DR_HOOK, JPHOOK
IMPLICIT NONE

TYPE(TGMV), INTENT(INOUT) :: SELF

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('YOMGMV:ZERO_YOMGMV',0,ZHOOK_HANDLE)

IF (ALLOCATED(SELF%GMV))  SELF%GMV  = 0.0_JPRB
IF (ALLOCATED(SELF%GMVS)) SELF%GMVS = 0.0_JPRB

IF (LHOOK) CALL DR_HOOK('YOMGMV:ZERO_YOMGMV',1,ZHOOK_HANDLE)

END SUBROUTINE ZERO_YOMGMV

!-------------------------------------------------------------------------

SUBROUTINE COPY_YOMGMV(SELF,RHS)

USE YOMHOOK,   ONLY : LHOOK, DR_HOOK, JPHOOK
IMPLICIT NONE

TYPE(TGMV), INTENT(INOUT) :: SELF
TYPE(TGMV), INTENT(IN)    :: RHS

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

IF (LHOOK) CALL DR_HOOK('YOMGMV:COPY_YOMGMV',0,ZHOOK_HANDLE)

IF (ALLOCATED(SELF%GMV) .AND. ALLOCATED(RHS%GMV)) THEN
  IF (ALL(SHAPE(SELF%GMV) == SHAPE(RHS%GMV))) THEN
    SELF%GMV = RHS%GMV
  ELSE
    WRITE(NULERR,*) 'SHAPE(SELF%GMV) ',SHAPE(SELF%GMV)
    WRITE(NULERR,*) 'SHAPE(RHS%GMV) ' ,SHAPE(RHS%GMV)
    CALL ABOR1("YOMGMV:COPY GMV different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%GMVS) .AND. ALLOCATED(RHS%GMVS)) THEN
  IF (ALL(SHAPE(SELF%GMVS) == SHAPE(RHS%GMVS))) THEN
    SELF%GMVS = RHS%GMVS
  ELSE
    CALL ABOR1("YOMGMV:COPY GMVS different shapes")
  ENDIF
ENDIF

IF (LHOOK) CALL DR_HOOK('YOMGMV:COPY_YOMGMV',1,ZHOOK_HANDLE)

END SUBROUTINE COPY_YOMGMV

!-------------------------------------------------------------------------

SUBROUTINE MUL_YOMGMV(SELF,PZ)

USE YOMHOOK,   ONLY : LHOOK, DR_HOOK, JPHOOK
IMPLICIT NONE

TYPE(TGMV),      INTENT(INOUT) :: SELF
REAL(KIND=JPRB), INTENT(IN)    :: PZ

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('YOMGMV:MUL_YOMGMV',0,ZHOOK_HANDLE)

IF (ALLOCATED(SELF%GMV))  SELF%GMV  = PZ*SELF%GMV
IF (ALLOCATED(SELF%GMVS)) SELF%GMVS = PZ*SELF%GMVS

IF (LHOOK) CALL DR_HOOK('YOMGMV:MUL_YOMGMV',1,ZHOOK_HANDLE)

END SUBROUTINE MUL_YOMGMV

!-------------------------------------------------------------------------

SUBROUTINE AXPBY_YOMGMV(SELF,PA,RHS,PB)

USE YOMHOOK,   ONLY : LHOOK, DR_HOOK, JPHOOK
IMPLICIT NONE

TYPE(TGMV),      INTENT(INOUT) :: SELF
REAL(KIND=JPRB), INTENT(IN)    :: PA
TYPE(TGMV),      INTENT(IN)    :: RHS
REAL(KIND=JPRB), INTENT(IN)    :: PB

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

IF (LHOOK) CALL DR_HOOK('YOMGMV:AXPBY_YOMGMV',0,ZHOOK_HANDLE)

IF (ALLOCATED(SELF%GMV) .AND. ALLOCATED(RHS%GMV)) THEN
  IF (ALL(SHAPE(SELF%GMV) == SHAPE(RHS%GMV))) THEN
    SELF%GMV = PA*SELF%GMV + PB*RHS%GMV
  ELSE
    CALL ABOR1("YOMGMV:AXPBY GMV different shapes")
  ENDIF
ENDIF
  
IF (ALLOCATED(SELF%GMVS) .AND. ALLOCATED(RHS%GMVS)) THEN
  IF (ALL(SHAPE(SELF%GMVS) == SHAPE(RHS%GMVS))) THEN
    SELF%GMVS = PA*SELF%GMVS + PB*RHS%GMVS
  ELSE
    CALL ABOR1("YOMGMV:AXPBY GMVS different shapes")
  ENDIF
ENDIF

IF (LHOOK) CALL DR_HOOK('YOMGMV:AXPBY_YOMGMV',1,ZHOOK_HANDLE)

END SUBROUTINE AXPBY_YOMGMV

!-------------------------------------------------------------------------

SUBROUTINE DIFF_YOMGMV(SELF,RHS)

USE YOMHOOK,   ONLY : LHOOK, DR_HOOK, JPHOOK
IMPLICIT NONE

TYPE(TGMV),      INTENT(INOUT) :: SELF
TYPE(TGMV),      INTENT(IN)    :: RHS

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

IF (LHOOK) CALL DR_HOOK('YOMGMV:DIFF_YOMGMV',0,ZHOOK_HANDLE)

IF (ALLOCATED(SELF%GMV) .AND. ALLOCATED(RHS%GMV)) THEN
  IF (ALL(SHAPE(SELF%GMV) == SHAPE(RHS%GMV))) THEN
    SELF%GMV = SELF%GMV - RHS%GMV
  ELSE
    CALL ABOR1("YOMGMV:DIFF GMV different shapes")
  ENDIF
ENDIF

IF (ALLOCATED(SELF%GMVS) .AND. ALLOCATED(RHS%GMVS)) THEN
  IF (ALL(SHAPE(SELF%GMVS) == SHAPE(RHS%GMVS))) THEN
    SELF%GMVS = SELF%GMVS - RHS%GMVS
  ELSE
    CALL ABOR1("YOMGMV:DIFF GMVS different shapes")
  ENDIF
ENDIF

IF (LHOOK) CALL DR_HOOK('YOMGMV:DIFF_YOMGMV',1,ZHOOK_HANDLE)

END SUBROUTINE DIFF_YOMGMV


!-------------------------------------------------------------------------

SUBROUTINE DOT_PROD_YOMGMV(YDGEOMETRY,FLD1,FLD2,PPROD)

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMHOOK,   ONLY : LHOOK, DR_HOOK, JPHOOK
IMPLICIT NONE

TYPE(GEOMETRY),  INTENT(IN)  :: YDGEOMETRY
TYPE(TGMV),      INTENT(IN)  :: FLD1
TYPE(TGMV),      INTENT(IN)  :: FLD2
REAL(KIND=JPRB), INTENT(OUT) :: PPROD

INTEGER(KIND=JPIM) :: JF
REAL(KIND=JPRB) :: ZTMP
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "dotprod2.intfb.h"
#include "dotprod3.intfb.h"

IF (LHOOK) CALL DR_HOOK('YOMGMV:DOT_PROD_YOMGMV',0,ZHOOK_HANDLE)

PPROD = 0.0_JPRB

IF (ALLOCATED(FLD1%GMV) .AND. ALLOCATED(FLD2%GMV)) THEN
  DO JF = 1, SIZE(FLD1%GMV, DIM=3)
    CALL DOTPROD3(YDGEOMETRY,JF,FLD1%GMV,FLD2%GMV,ZTMP)
    PPROD = PPROD + ZTMP
  ENDDO
ENDIF

IF (ALLOCATED(FLD1%GMVS) .AND. ALLOCATED(FLD2%GMVS)) THEN
  DO JF = 1, SIZE(FLD1%GMVS, DIM=2)
    CALL DOTPROD2(YDGEOMETRY,JF,FLD1%GMVS,FLD2%GMVS,ZTMP)
    PPROD = PPROD + ZTMP
  ENDDO
ENDIF

IF (LHOOK) CALL DR_HOOK('YOMGMV:DOT_PROD_YOMGMV',1,ZHOOK_HANDLE)

END SUBROUTINE DOT_PROD_YOMGMV

!-------------------------------------------------------------------------

SUBROUTINE RANDOM_YOMGMV(SELF)

USE YOMHOOK,   ONLY : LHOOK, DR_HOOK, JPHOOK
IMPLICIT NONE

TYPE(TGMV), INTENT(INOUT) :: SELF

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('YOMGMV:RANDOM_YOMGMV',0,ZHOOK_HANDLE)

! Just a plain RANOM_NUMBER for now. 

IF (ALLOCATED(SELF%GMV))  THEN
  CALL RANDOM_NUMBER(SELF%GMV)
  SELF%GMV = 2.0_JPRB*SELF%GMV - 1.0_JPRB
ENDIF
IF (ALLOCATED(SELF%GMVS))  THEN
  CALL RANDOM_NUMBER(SELF%GMVS)
  SELF%GMVS = 2.0_JPRB*SELF%GMVS - 1.0_JPRB
ENDIF

IF (LHOOK) CALL DR_HOOK('YOMGMV:RANDOM_YOMGMV',1,ZHOOK_HANDLE)

END SUBROUTINE RANDOM_YOMGMV

END MODULE YOMGMV

