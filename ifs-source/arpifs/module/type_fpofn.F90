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

MODULE TYPE_FPOFN

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------

! CSH : output file name containing spectral (hybrid) data
! CGG : output file name containing gridpoint (surface) data
! CSX : output file name containing gridpoint SURFEX data
! CUA : output file name containing upper air gridpoint data

! CLI : input file name containing surface climatology data on target geometry
! CSU : input file name containing SURFEX climatology data on target geometry

!====== FULLPOS INPUT/OUTPUT FILES  ======

TYPE TFPOFN

CHARACTER(LEN=256)  :: CSH
CHARACTER(LEN=256)  :: CGG
CHARACTER(LEN=256)  :: CSX
CHARACTER(LEN=256)  :: CUA
CHARACTER(LEN=256)  :: CLI
CHARACTER(LEN=256)  :: CSU

END TYPE TFPOFN

!     ------------------------------------------------------------------
END MODULE TYPE_FPOFN
