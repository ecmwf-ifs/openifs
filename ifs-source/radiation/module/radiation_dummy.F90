! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

! The IFS does not use the "PSRAD" capability of the radiation scheme;
! it is hidden behind #ifdef...#endif in radiation_interface.
! However, the interface generation script does not know about
! preprocessor directives, so believes that radiation_interface
! depends on a module radiation_psrad_rrtm.mod.  Therefore we need to
! create a dummy mod file that won't actually be used.

module radiation_psrad_rrtm
  integer :: this_will_never_be_used
end module radiation_psrad_rrtm
