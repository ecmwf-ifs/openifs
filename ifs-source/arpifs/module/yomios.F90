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

MODULE YOMIOS

IMPLICIT NONE

SAVE

! -------- RESTART FILES -------------------
! CFRCF    : pathname for restart control file
! CIOSPRF  : PATHNAME FOR SPECTRAL RESTART FILE
! LRCFTIME : ALSO WRITE A RESTART CONTROL FILE WITH TIME

CHARACTER (LEN = 120) ::  CIOSPRF
CHARACTER (LEN = 120) ::  CFRCF
LOGICAL :: LRCFTIME
!     ------------------------------------------------------------------
END MODULE YOMIOS
