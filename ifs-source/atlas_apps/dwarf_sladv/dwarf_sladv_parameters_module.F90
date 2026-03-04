! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

module dwarf_sladv_parameters_module

use dwarf_sladv_kind_module

implicit none

real(wp), parameter, public :: rpi=2.0_wp*asin(1.0_wp)
real(wp), parameter, public :: Deg2Rad=RPI/180._wp
real(wp), parameter, public :: Rad2Deg=180._wp/RPI
real(wp), parameter, public :: rndoff = epsilon(1.0_wp)
real(wp), parameter, public :: earth_radius = 6371229._wp

end module
