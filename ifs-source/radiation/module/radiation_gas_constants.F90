! radiation_gas_constants.F90 - Molar mases and ID codes of the various gases
!
! (C) Copyright 2014- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
! License: see the COPYING file for details
!

module radiation_gas_constants

  use parkind1, only : jprb

  implicit none

  public

  ! Gas codes; these indices match those of RRTM-LW up to 7
  integer, parameter :: IGasNotPresent = 0
  integer, parameter :: IH2O   = 1
  integer, parameter :: ICO2   = 2
  integer, parameter :: IO3    = 3
  integer, parameter :: IN2O   = 4
  integer, parameter :: ICO    = 5
  integer, parameter :: ICH4   = 6
  integer, parameter :: IO2    = 7
  integer, parameter :: ICFC11 = 8
  integer, parameter :: ICFC12 = 9
  integer, parameter :: IHCFC22= 10
  integer, parameter :: ICCl4  = 11 
  integer, parameter :: INO2   = 12
  integer, parameter :: NMaxGases = 12

  ! Molar masses (g mol-1) of dry air and the various gases above
  real(jprb), parameter :: AirMolarMass = 28.970_jprb
  real(jprb), parameter, dimension(0:NMaxGases) :: GasMolarMass = (/ &
       & 0.0_jprb,        & ! Gas not present
       & 18.0152833_jprb, & ! H2O
       & 44.011_jprb,     & ! CO2
       & 47.9982_jprb,    & ! O3
       & 44.013_jprb,     & ! N2O
       & 28.0101_jprb,    & ! CO
       & 16.043_jprb,     & ! CH4
       & 31.9988_jprb,    & ! O2
       & 137.3686_jprb,   & ! CFC11
       & 120.914_jprb,    & ! CFC12
       & 86.469_jprb,     & ! HCFC22
       & 153.823_jprb,    & ! CCl4    
       & 46.0055_jprb /)    ! NO2

  ! The corresponding names of the gases in upper and lower case, used
  ! for reading variables from the input file
  character*6, dimension(NMaxGases), parameter :: GasName &
       &  = (/'H2O   ','CO2   ','O3    ','N2O   ','CO    ','CH4   ', &
       &      'O2    ','CFC11 ','CFC12 ','HCFC22','CCl4  ','NO2   '/)
  character*6, dimension(NMaxGases), parameter :: GasLowerCaseName &
       &  = (/'h2o   ','co2   ','o3    ','n2o   ','co    ','ch4   ', &
       &      'o2    ','cfc11 ','cfc12 ','hcfc22','ccl4  ','no2   '/)

end module radiation_gas_constants
