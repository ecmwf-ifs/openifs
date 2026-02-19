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

MODULE SPECTRAL_VARIABLES_MOD

USE PARKIND1, ONLY: JPIM

IMPLICIT NONE

INTEGER(KIND=JPIM), PARAMETER :: JP_MAX_N_SPJB_VARS = 100

TYPE SPECTRAL_VARIABLES
  INTEGER(KIND=JPIM) :: NS3D
  INTEGER(KIND=JPIM) :: NS2D
  INTEGER(KIND=JPIM) :: NS1D
  INTEGER(KIND=JPIM) :: NGRBVAR(JP_MAX_N_SPJB_VARS)
END TYPE SPECTRAL_VARIABLES

END MODULE SPECTRAL_VARIABLES_MOD
