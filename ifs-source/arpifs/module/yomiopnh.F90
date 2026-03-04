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

MODULE YOMIOPNH

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------

! -------- TRAJECTORY FOR ADJOINT OF NONHYDROSTATIC DYN. --------

! NG3NH95   : NUMBER OF 3D FIELDS - (NFLEVG)
! LTRAJNH   : .T. IF NHS TRAJECTORY IS WRITTEN INTO BUFFER
!             NEEDED FOR TL AND AD MODELS (TRAJNH AND TRAJIT)

INTEGER(KIND=JPIM) :: NG3NH95

LOGICAL :: LTRAJNH

!     ------------------------------------------------------------------
END MODULE YOMIOPNH
