MODULE YOMCST
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

!     ------------------------------------------------------------------

!*    Common of physical constants
!     You will find the meanings in the annex 1 of the documentation

! A1.0 Fundamental constants
REAL(KIND=JPRB) :: RPI
REAL(KIND=JPRB) :: RCLUM
REAL(KIND=JPRB) :: RHPLA
REAL(KIND=JPRB) :: RKBOL
REAL(KIND=JPRB) :: RNAVO
! A1.1 Astronomical constants
REAL(KIND=JPRB) :: RDAY
REAL(KIND=JPRB) :: REA
REAL(KIND=JPRB) :: REPSM
REAL(KIND=JPRB) :: RSIYEA
REAL(KIND=JPRB) :: RSIDAY
REAL(KIND=JPRB) :: ROMEGA
! A1.2 Geoide
REAL(KIND=JPRB) :: RA
REAL(KIND=JPRB) :: RG
REAL(KIND=JPRB) :: R1SA
! A1.3 Radiation
REAL(KIND=JPRB) :: RSIGMA
REAL(KIND=JPRB) :: RI0
! A1.4 Thermodynamic gas phase
REAL(KIND=JPRB) :: R
REAL(KIND=JPRB) :: RMD
REAL(KIND=JPRB) :: RMCO2
REAL(KIND=JPRB) :: RMV
REAL(KIND=JPRB) :: RMO3
REAL(KIND=JPRB) :: RD
REAL(KIND=JPRB) :: RV
REAL(KIND=JPRB) :: RCPD
REAL(KIND=JPRB) :: RCPV
REAL(KIND=JPRB) :: RCVD
REAL(KIND=JPRB) :: RCVV
REAL(KIND=JPRB) :: RKAPPA
REAL(KIND=JPRB) :: RETV
! A1.5,6 Thermodynamic liquid,solid phases
REAL(KIND=JPRB) :: RCW
REAL(KIND=JPRB) :: RCS
! A1.7 Thermodynamic transition of phase
REAL(KIND=JPRB) :: RLVTT
REAL(KIND=JPRB) :: RLSTT
REAL(KIND=JPRB) :: RLVZER
REAL(KIND=JPRB) :: RLSZER
REAL(KIND=JPRB) :: RLMLT
REAL(KIND=JPRB) :: RTT
REAL(KIND=JPRB) :: RATM
REAL(KIND=JPRB) :: RDT
! A1.8 Curve of saturation
REAL(KIND=JPRB) :: RESTT
REAL(KIND=JPRB) :: RALPW
REAL(KIND=JPRB) :: RBETW
REAL(KIND=JPRB) :: RGAMW
REAL(KIND=JPRB) :: RALPS
REAL(KIND=JPRB) :: RBETS
REAL(KIND=JPRB) :: RGAMS
REAL(KIND=JPRB) :: RALPD
REAL(KIND=JPRB) :: RBETD
REAL(KIND=JPRB) :: RGAMD
! CMPI6 ANNUAL GLOBAL MEAN CO2 (1950-2030)
INTEGER(KIND=JPIM):: RCO2REFYEAR
INTEGER(KIND=JPIM), PARAMETER :: RNCO2YEARS=81
REAL(KIND=JPRB) :: RANCO2(RNCO2YEARS)


!    ------------------------------------------------------------------
END MODULE YOMCST
