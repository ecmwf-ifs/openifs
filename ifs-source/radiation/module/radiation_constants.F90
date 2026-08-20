! radiation_constants.F90 - Constants used in radiation calculations
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
!

module radiation_constants

  use parkind1, only : jprb
!  use yomcst,   only : RPI, RSIGMA, RG, RD, RMV, RMO3, RHPLA

  implicit none
  public

  ! Rename some constants from their cryptic IFS names
  real(jprb), parameter :: Pi                 = 3.14159265358979323846_jprb
  real(jprb), parameter :: AccelDueToGravity  = 9.80665_jprb! m s-2
  real(jprb), parameter :: StefanBoltzmann    = 5.67037321e-8_jprb ! W m-2 K-4
  real(jprb), parameter :: DensityLiquidWater = 1000.0_jprb ! kg m-3
  real(jprb), parameter :: DensitySolidIce    = 916.7_jprb  ! kg m-3
  real(jprb), parameter :: GasConstantDryAir  = 287.058_jprb! J kg-1 K-1
  real(jprb), parameter :: PlanckConstant     = 6.6260695729e-34_jprb ! J s
  real(jprb), parameter :: BoltzmannConstant  = 1.380648813e-23_jprb ! J K-1
  real(jprb), parameter :: SpeedOfLight       = 299792458.0_jprb ! m s-1

end module radiation_constants
