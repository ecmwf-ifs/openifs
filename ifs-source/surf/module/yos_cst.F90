MODULE YOS_CST

USE PARKIND1  ,ONLY : JPIM     ,JPRB

! (C) Copyright 2005- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
  
IMPLICIT NONE

SAVE

TYPE :: TCST 
REAL(KIND=JPRB) :: RTT     ! RT0=273.16= TRIPLE POINT TEMPERATURE
REAL(KIND=JPRB) :: RPI     ! PI=3.14...
REAL(KIND=JPRB) :: RDAY    ! Astronomical constants
REAL(KIND=JPRB) :: R       ! Thermodynamic gas phase
REAL(KIND=JPRB) :: RD      ! Thermodynamic gas phase
REAL(KIND=JPRB) :: RV      ! Thermodynamic gas phase
REAL(KIND=JPRB) :: RETV    ! Thermodynamic gas phase
REAL(KIND=JPRB) :: RLSTT   ! Thermodynamic transition of phase
REAL(KIND=JPRB) :: RLMLT   ! Thermodynamic transition of phase
REAL(KIND=JPRB) :: RLVTT   ! Thermodynamic transition of phase
REAL(KIND=JPRB) :: RCPD    ! Thermodynamic gas phase
REAL(KIND=JPRB) :: RCPV    ! Thermodynamic gas phase
REAL(KIND=JPRB) :: RSIGMA  ! Radiation
REAL(KIND=JPRB) :: RG      ! Geoide
REAL(KIND=JPRB) :: RATM    ! 
REAL(KIND=JPRB) :: ROMEGA
END TYPE TCST 

END MODULE YOS_CST
