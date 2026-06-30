MODULE YOMCDH1S
! (C) Copyright 2005- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE
SAVE

!*    sizes of diagnostic arrays from the 28r3 cycle, used to store
!     output from callpar1s
 
!    *NLEVI*        Number of sea ice layers (diagnostics)
!    *NDHVTLS*      Number of variables for individual tiles
!    *NDHFTLS*      Number of fluxes for individual tiles
!    *NDHVTSS*      Number of variables for snow energy budget
!    *NDHFTSS*      Number of fluxes for snow energy budget
!    *NDHVTTS*      Number of variables for soil energy budget
!    *NDHFTTS*      Number of fluxes for soil energy budget
!    *NDHVTIS*      Number of variables for sea ice energy budget
!    *NDHFTIS*      Number of fluxes for sea ice energy budget
!    *NDHVIIS*      Number of variables for interception energy budget
!    *NDHFIIS*      Number of fluxes for interception energy budget
!    *NDHVSSS*      Number of variables for snow mass energy budget
!    *NDHFSSS*      Number of fluxes for snow mass energy budget
!    *NDHVWLS*      Number of variables for soil water energy budget
!    *NDHFWLS*      Number of fluxes for soil water energy budget
!    *NDHVRESS*     Number of variables for resistances
!    *NDHFRESS*     Number of fluxes for resistances

!    *NDHVCO2S*     Number of variables for CO2
!    *NDHFCO2S*     Number of fluxes for CO2
!    *NDHVBIOS*     Number of variables for biomass
!    *NDHFBIOS*     Number of fluxes for biomass
!    *NDHVVEGS*     Number of variables for vegetation
!    *NDHFVEGS*     Number of fluxes for vegetation

INTEGER(KIND=JPIM) :: NLEVI
INTEGER(KIND=JPIM) :: NDHVTLS
INTEGER(KIND=JPIM) :: NDHFTLS
INTEGER(KIND=JPIM) :: NDHVTSS
INTEGER(KIND=JPIM) :: NDHFTSS
INTEGER(KIND=JPIM) :: NDHVTTS
INTEGER(KIND=JPIM) :: NDHFTTS
INTEGER(KIND=JPIM) :: NDHVTIS
INTEGER(KIND=JPIM) :: NDHFTIS
INTEGER(KIND=JPIM) :: NDHVIIS
INTEGER(KIND=JPIM) :: NDHFIIS
INTEGER(KIND=JPIM) :: NDHVSSS
INTEGER(KIND=JPIM) :: NDHFSSS
INTEGER(KIND=JPIM) :: NDHVWLS
INTEGER(KIND=JPIM) :: NDHFWLS
INTEGER(KIND=JPIM) :: NDHVRESS
INTEGER(KIND=JPIM) :: NDHFRESS

INTEGER(KIND=JPIM) :: NDHVCO2S
INTEGER(KIND=JPIM) :: NDHFCO2S
INTEGER(KIND=JPIM) :: NDHVBIOS
INTEGER(KIND=JPIM) :: NDHFBIOS
INTEGER(KIND=JPIM) :: NDHVVEGS
INTEGER(KIND=JPIM) :: NDHFVEGS
END MODULE YOMCDH1S
