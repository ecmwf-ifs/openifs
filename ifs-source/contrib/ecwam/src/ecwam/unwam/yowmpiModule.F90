! (C) Copyright 2001- Aron Roland (Roland & Partner, Germany).
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

!> \file yowmpimodule
!> \brief just a wrapper for different mpi modules

module yowMpiModule
#ifdef WAM_HAVE_MPI_F08
  use mpi_f08
#else
  use MPL_MPIF
#endif
end module
