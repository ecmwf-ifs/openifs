! (C) Copyright 2001- Aron Roland (Roland & Partner, Germany).
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

!> \file yowpd.F90
!> \brief This is an interface module. Simply call use yowpd in your fortran code to use pdlib
!> \author Thomas Huxhorn
!> \date 2013

module yowpd
  use yowError
  use yowDatapool
  use yowpdlibMain
  use yowNodepool
  use yowElementpool
  use yowSidepool
  use yowExchangeModule
  use yowRankModule
!  comment out because it cause cyclic dependency
!   use yowMpiModule
  use yowChecksModule
  implicit none
  public
end module yowpd
