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

MODULE TYPE_ECV

!     Purpose.
!     --------
!       Data and controls for extended control variable

!     Author.
!     -------
!       S. Massart

!     Modifications.
!     --------------
!       Original    April-2018
!       S. Massart 19-Feb-2019 : Addition of more types
! ------------------------------------------------------------------


USE PARKIND1, ONLY: JPIM, JPRB

IMPLICIT NONE

SAVE

! ------------------------------------------------------------------
INTERFACE ASSIGNMENT (=)
MODULE PROCEDURE ASSIGN_ECV_ECV
END INTERFACE

INTERFACE OPERATOR (.EQV.)
MODULE PROCEDURE EQUIV_ECV
END INTERFACE

! ------------------------------------------------------------------

TYPE, PUBLIC :: TECVDIM
  !   ------------------------------------------------------------------
  !   LECV_1D     True to add  1D fields in the control variable
  !   LECV_2D     True to add  2D fields in the control variable
  !   LECV_3D     True to add  3D fields in the control variable
  !
  !   NECV_1D     Number of 1D fields to add 
  !   NECVSP_2D   Number of 2D spectral fields 
  !   NECVGP_2D   Number of 2D gird point fields 
  !   NECV_2D     Number of 2D fields (spectral + grib points)
  !   NECVSP_3D   Number of 3D spectral fields 
  !   NECVGP_3D   Number of 3D gird point fields 
  !   NECV_3D     Number of 3D fields (spectral + grib points)
  !   ------------------------------------------------------------------
  LOGICAL              :: LECV_1D
  LOGICAL              :: LECV_2D
  LOGICAL              :: LECV_3D
  INTEGER(KIND=JPIM)   :: NECV_1D
  INTEGER(KIND=JPIM)   :: NECVSP_2D
  INTEGER(KIND=JPIM)   :: NECVGP_2D
  INTEGER(KIND=JPIM)   :: NECV_2D
  INTEGER(KIND=JPIM)   :: NECVSP_3D
  INTEGER(KIND=JPIM)   :: NECVGP_3D
  INTEGER(KIND=JPIM)   :: NECV_3D
  INTEGER(KIND=JPIM)   :: NECV_ALL
END TYPE TECVDIM

TYPE TECV_CONFIG_BASE
  !   ------------------------------------------------------------------
  !   L_IN_1D     True if in 1D 
  !   L_IN_2D_SP  True if in spectral 2D
  !   L_IN_2D_GP  True if in grid point 2D
  !   L_IN_3D_SP  True if in spectral 3D
  !   L_IN_3D_GP  True if in grid point 3D
  !   CSETDESC    Description of each variable
  !   CWNAME      Name in the wavelet file
  !
  !   ------------------------------------------------------------------
  LOGICAL              :: L_IN_1D
  LOGICAL              :: L_IN_2D_SP
  LOGICAL              :: L_IN_2D_GP
  LOGICAL              :: L_IN_3D_SP
  LOGICAL              :: L_IN_3D_GP
  CHARACTER(LEN =20)   :: CSETDESC
  CHARACTER(LEN =20)   :: CWNAME
END TYPE TECV_CONFIG_BASE

TYPE TECV_CONFIG
  !   ------------------------------------------------------------------
  !   NECV        Total number of ECV 
  !   YRCONFIG    Configuration array
  !
  !   ------------------------------------------------------------------
  INTEGER(KIND=JPIM)   :: NECV
  TYPE(TECV_CONFIG_BASE), ALLOCATABLE :: YRCONFIG(:)
END TYPE TECV_CONFIG

TYPE, PUBLIC :: TECVGRIB
  !   ------------------------------------------------------------------
  !   MECVGRB_1D   Array with the grib numbers of the 1D fields
  !   MECVGRB_2D   Array with the grib numbers of the 2D fields (first SP and then GP)
  !   MECVGRB_2D   Array with the grib numbers of the 3D fields (first SP and then GP)
  !   MECVGRB_ALL  Array with the grib numbers of the all fields (1D + 2D + 3D)
  !   ------------------------------------------------------------------
  INTEGER(KIND=JPIM), ALLOCATABLE :: MECVGRB_1D(:)
  INTEGER(KIND=JPIM), ALLOCATABLE :: MECVGRB_2D(:)
  INTEGER(KIND=JPIM), ALLOCATABLE :: MECVGRB_3D(:)
  INTEGER(KIND=JPIM), ALLOCATABLE :: MECVGRB_ALL(:)
END TYPE TECVGRIB

TYPE, PUBLIC :: ECV_CONTAINER
  !     ------------------------------------------------------------------
  !   RECV1D      ECV 1D fields
  !   RSPECV2D    ECV 2D spectral fields
  !   RGPECV2D    ECV 2D grid point fields
  !   RSPECV3D    ECV 3D spectral fields
  !   RGPECV3D    ECV 3D grid point fields
  !     ------------------------------------------------------------------
  REAL(KIND=JPRB), ALLOCATABLE :: RECV1D(:)
  REAL(KIND=JPRB), ALLOCATABLE :: RSPECV2D(:,:)
  REAL(KIND=JPRB), ALLOCATABLE :: RGPECV2D(:,:,:)
  REAL(KIND=JPRB), ALLOCATABLE :: RSPECV3D(:,:,:)
  REAL(KIND=JPRB), ALLOCATABLE :: RGPECV3D(:,:,:,:)
  REAL(KIND=JPRB), ALLOCATABLE :: RLSM(:,:)
END TYPE ECV_CONTAINER

!-----------------------------------------------------------------------
CONTAINS
!-----------------------------------------------------------------------

SUBROUTINE INITIALIZE_ECVDIM(YDDIMECV)

USE YOMHOOK       , ONLY : LHOOK, DR_HOOK, JPHOOK
IMPLICIT NONE

TYPE(TECVDIM)  ,INTENT(INOUT)    :: YDDIMECV
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('TYPE_ECV:INITIALIZE_ECVDIM',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

YDDIMECV%LECV_1D = .FALSE.
YDDIMECV%LECV_2D = .FALSE.
YDDIMECV%LECV_3D = .FALSE.
YDDIMECV%NECV_1D   = 0_JPIM
YDDIMECV%NECVSP_2D = 0_JPIM
YDDIMECV%NECVGP_2D = 0_JPIM
YDDIMECV%NECV_2D   = 0_JPIM
YDDIMECV%NECVSP_3D = 0_JPIM
YDDIMECV%NECVGP_3D = 0_JPIM
YDDIMECV%NECV_3D   = 0_JPIM
YDDIMECV%NECV_ALL  = 0_JPIM

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('TYPE_ECV:INITIALIZE_ECVDIM',1,ZHOOK_HANDLE)

END SUBROUTINE INITIALIZE_ECVDIM

!-----------------------------------------------------------------------

SUBROUTINE ASSIGN_ECV_ECV(YDECV1,YDECV2)

USE YOMHOOK       , ONLY : LHOOK, DR_HOOK, JPHOOK
IMPLICIT NONE

TYPE(ECV_CONTAINER)  ,INTENT(INOUT) :: YDECV1
TYPE(ECV_CONTAINER)  ,INTENT(IN)    :: YDECV2
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('TYPE_ECV:ASSIGN_ECV_ECV',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

IF (ALLOCATED(YDECV1%RECV1D) .AND. ALLOCATED(YDECV2%RECV1D)) &
  & YDECV1%RECV1D(:)       = YDECV2%RECV1D(:)
IF (ALLOCATED(YDECV1%RSPECV2D) .AND. ALLOCATED(YDECV2%RSPECV2D)) &
  & YDECV1%RSPECV2D(:,:)   = YDECV2%RSPECV2D(:,:)
IF (ALLOCATED(YDECV1%RGPECV2D) .AND. ALLOCATED(YDECV2%RGPECV2D)) &
  & YDECV1%RGPECV2D(:,:,:)   = YDECV2%RGPECV2D(:,:,:)
IF (ALLOCATED(YDECV1%RSPECV3D) .AND. ALLOCATED(YDECV2%RSPECV3D)) &
  & YDECV1%RSPECV3D(:,:,:) = YDECV2%RSPECV3D(:,:,:)
IF (ALLOCATED(YDECV1%RGPECV3D) .AND. ALLOCATED(YDECV2%RGPECV3D)) &
  & YDECV1%RGPECV3D(:,:,:,:) = YDECV2%RGPECV3D(:,:,:,:)


!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('TYPE_ECV:ASSIGN_ECV_ECV',1,ZHOOK_HANDLE)

END SUBROUTINE ASSIGN_ECV_ECV

! ------------------------------------------------------------------

LOGICAL FUNCTION EQUIV_ECV(YDECV1,YDECV2)

USE YOMHOOK       , ONLY : LHOOK, DR_HOOK, JPHOOK
IMPLICIT NONE

TYPE(ECV_CONTAINER), INTENT(IN) :: YDECV1
TYPE(ECV_CONTAINER), INTENT(IN) :: YDECV2
LOGICAL :: LL

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('TYPE_ECV:EQUIV_ECV',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

LL = .TRUE.

LL = LL .AND. (ALLOCATED(YDECV1%RECV1D) .AND. ALLOCATED(YDECV2%RECV1D))
LL = LL .AND. (ALLOCATED(YDECV1%RSPECV2D) .AND. ALLOCATED(YDECV2%RSPECV2D)) 
LL = LL .AND. (ALLOCATED(YDECV1%RGPECV2D) .AND. ALLOCATED(YDECV2%RGPECV2D)) 
LL = LL .AND. (ALLOCATED(YDECV1%RSPECV3D) .AND. ALLOCATED(YDECV2%RSPECV3D)) 
LL = LL .AND. (ALLOCATED(YDECV1%RGPECV3D) .AND. ALLOCATED(YDECV2%RGPECV3D)) 

IF (ALLOCATED(YDECV1%RECV1D) .AND. ALLOCATED(YDECV2%RECV1D)) &
  &  LL = LL .AND. ALL(YDECV1%RECV1D ==YDECV2%RECV1D)
IF (ALLOCATED(YDECV1%RSPECV2D) .AND. ALLOCATED(YDECV2%RSPECV2D)) &
  &  LL = LL .AND. ALL(YDECV1%RSPECV2D ==YDECV2%RSPECV2D)
IF (ALLOCATED(YDECV1%RGPECV2D) .AND. ALLOCATED(YDECV2%RGPECV2D)) &
  &  LL = LL .AND. ALL(YDECV1%RGPECV2D ==YDECV2%RGPECV2D)
IF (ALLOCATED(YDECV1%RSPECV3D) .AND. ALLOCATED(YDECV2%RSPECV3D)) &
  &  LL = LL .AND. ALL(YDECV1%RSPECV3D ==YDECV2%RSPECV3D)
IF (ALLOCATED(YDECV1%RGPECV3D) .AND. ALLOCATED(YDECV2%RGPECV3D)) &
  &  LL = LL .AND. ALL(YDECV1%RGPECV3D ==YDECV2%RGPECV3D)

EQUIV_ECV=LL

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('TYPE_ECV:EQUIV_ECV',1,ZHOOK_HANDLE)
END FUNCTION EQUIV_ECV

END MODULE TYPE_ECV
