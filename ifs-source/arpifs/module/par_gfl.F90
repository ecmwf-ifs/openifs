! (C) Copyright 2017- ECMWF.
! (C) Copyright 2017- Meteo-France.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE PAR_GFL

USE PARKIND1 , ONLY : JPIM

IMPLICIT NONE

! Parametars used by GFL (and elsewhere)

! JPGFL : Max number of GFL fields
! JPNAMED_GFL : Number of currently pre-defined components of GFL
! JPGHG : Number of greenhouse gas fields
! JPCHEM : Number of chemical species
! JPAERO : Number of active aerosol fields
! JPAEROUT: Number of output aerosol fields
! JPAEROCLIM: Number of 3D climatological aerosol fields
! JPUVP : Number of output from UV processor
! JPERA40 : Number of ERA40 diagnostic fields
! JPCH4S  : Number of added fields related to methane
! JPNOGW  : Number of diagnostic fields for NORO GWD SCHEME
! JPSLDIA : Number of SL dynamics diagnostic fields
! JPCHEM_ASSIM : Maximum number of assimilated of chemical species
! JPGHG_ASSIM  : Maximum number of assimilated of GHG species
! JPFSD   : Number of cloud heterogeneity fields
! JPEDRP  : Maximum number of EDR Parameters diagnostic turbulence fields
!-------------------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: JPGFL=2163
INTEGER(KIND=JPIM), PARAMETER :: JPNAMED_GFL=27
INTEGER(KIND=JPIM), PARAMETER :: JPGHG=3
INTEGER(KIND=JPIM), PARAMETER :: JPCHEM=220
INTEGER(KIND=JPIM), PARAMETER :: JPGHG_ASSIM=2
INTEGER(KIND=JPIM), PARAMETER :: JPCHEM_ASSIM=6
INTEGER(KIND=JPIM), PARAMETER :: JPAERO=42
INTEGER(KIND=JPIM), PARAMETER :: JPFORC=1100
INTEGER(KIND=JPIM), PARAMETER :: JPERA40=14
INTEGER(KIND=JPIM), PARAMETER :: JPSLDIA=7
INTEGER(KIND=JPIM), PARAMETER :: JPEZDIAG=50
INTEGER(KIND=JPIM), PARAMETER :: JPCH4S=1
INTEGER(KIND=JPIM), PARAMETER :: JPNOGW=2
INTEGER(KIND=JPIM), PARAMETER :: JPAEROUT=38
INTEGER(KIND=JPIM), PARAMETER :: JPAEROCLIM=3
INTEGER(KIND=JPIM), PARAMETER :: JPUVP=2
INTEGER(KIND=JPIM), PARAMETER :: JPPHYS=9
INTEGER(KIND=JPIM), PARAMETER :: JPPHYCTY=1
INTEGER(KIND=JPIM), PARAMETER :: JPFSD=1
INTEGER(KIND=JPIM), PARAMETER :: JPEDRP=2
INTEGER(KIND=JPIM), PARAMETER :: JPLIMA=50
INTEGER(KIND=JPIM), PARAMETER :: GRIB_CODE_GFL_PHYS=81  ! AJGDB hopefully harmless

END MODULE PAR_GFL
