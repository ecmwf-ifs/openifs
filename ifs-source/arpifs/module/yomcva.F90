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

MODULE YOMCVA

USE PARKIND1  ,ONLY : JPIM

USE CONTROL_VECTORS_BASE_MIX, ONLY : CONTROL_VECTOR

IMPLICIT NONE
SAVE

!     ------------------------------------------------------------------
!      YOMCVA
!
! Modifications:
!   15-Feb-2013 M. Fisher : Move globals into structures
!   19-Feb-2019 S. Massart: Paramater optimisation
!     ------------------------------------------------------------------

! NVA3SP        : Number of 3D fields in the control vector which
!                 are spectral fields in the model (ie in SPA3).
! NVA3D         : Total number of 3D fields in the control vector.
! NVA2D         : Total number of 2D fields in the control vector.
! NVA1D         : Dimension of LAM Aladin mean wind components in ctl vec
! NVATOV        : DIMENSION OF TOVS CONTOL VARIABLE (part of NVADIM)
! NVATOVV       : NUMBER OF TOVS CONTOL VARIABLES 
! NVAPARAM      : Number of scalar parameters in control variable
! NVPARECV      : Number of ECV parameters in control variable
! NGRBVAR       : Grib codes of fields in control variable

TYPE CVA_DATA_TYPE
  INTEGER(KIND=JPIM) :: NVA3D
  INTEGER(KIND=JPIM) :: NVA2D
  INTEGER(KIND=JPIM) :: NVA1D
  INTEGER(KIND=JPIM) :: NVATOV
  INTEGER(KIND=JPIM) :: NVATOVV
  INTEGER(KIND=JPIM) :: NVAPARAM
  INTEGER(KIND=JPIM) :: NVPARECV
  INTEGER(KIND=JPIM), POINTER :: NGRBVAR(:) => NULL()
END TYPE CVA_DATA_TYPE

TYPE CVA_STRUCT_TYPE
  TYPE(CONTROL_VECTOR) :: YVAZX ! control variable
  TYPE(CONTROL_VECTOR) :: YVAZG ! gradient of control variable
  TYPE(CONTROL_VECTOR) :: YVAZ0
END TYPE CVA_STRUCT_TYPE

TYPE SCALP_STRUCT_TYPE
  TYPE(CONTROL_VECTOR) :: YSCALP      ! scalar product array
  TYPE(CONTROL_VECTOR) :: YRSCALP     ! inverse of YSCALP
  TYPE(CONTROL_VECTOR) :: YSCALPSQRT  ! square root of YSCALP
  TYPE(CONTROL_VECTOR) :: YRSCALPSQRT ! inverse of YSCALPSQRT
END TYPE SCALP_STRUCT_TYPE

TYPE(CVA_DATA_TYPE), POINTER :: CVA_DATA => NULL()
TYPE(CVA_STRUCT_TYPE), POINTER :: CVA_STRUCT => NULL()
TYPE(SCALP_STRUCT_TYPE), POINTER :: SCALP_STRUCT => NULL()

!     ------------------------------------------------------------------
END MODULE YOMCVA
