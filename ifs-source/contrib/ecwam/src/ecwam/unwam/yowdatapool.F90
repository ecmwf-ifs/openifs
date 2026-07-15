! (C) Copyright 2001- Aron Roland (Roland & Partner, Germany).
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

!> Has fancy data
module yowDatapool
  USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
  use yowMpiModule
  implicit none


#ifdef USE_SINGLE
  !> single precision. Enable with compiler flag -DUSE_SINGLE
  integer(KIND=JWIM),parameter :: rkind = JWRB
#else
  !> double precision. Default real datatype
  integer(KIND=JWIM),parameter :: rkind = JWRU
#endif


  logical, parameter :: debugPrePartition = .false.
  logical, parameter :: debugPostPartition = .true.
  logical, parameter :: debugParmetis = .false.
  !> write partition information into file fort.600.
  !> one can display partition with katerfempresenter. just open the system.dat and
  !> click on the partitions icon
  logical, parameter :: debugPartition = .false.

  !> Number of threads
  integer(KIND=JWIM), save :: nTasks = 0

  !> The thread id.
  !> starts by 0. The first Thread has rank 0
  integer(KIND=JWIM), save :: myrank = 0

  !> MPI Communicator.
  !> Should be MPI_COMM_WORLD. If pdlib is run into a existing MPI enviroment, comm is set to a new communicator
#ifdef WAM_HAVE_MPI_F08
  type(mpi_comm),public,save :: comm = MPI_COMM_WORLD
#else
  integer(KIND=JWIM),public,save :: comm = MPI_COMM_WORLD
#endif
  !> MPI integer Type.
  !> Should be MPI_integer
#ifdef WAM_HAVE_MPI_F08
  type(mpi_datatype),parameter :: itype = MPI_integer
#else
  integer(KIND=JWIM),parameter :: itype = MPI_integer
#endif
  !> MPI Real Type
  !> Shpuld be MPI_REAL8
#ifdef USE_SINGLE
#ifdef WAM_HAVE_MPI_F08
  type(mpi_datatype), parameter :: rtype = MPI_REAL4
#else
  integer(KIND=JWIM), parameter :: rtype = MPI_REAL4
#endif
#else
#ifdef WAM_HAVE_MPI_F08
  type(mpi_datatype), parameter :: rtype = MPI_REAL8
#else
  integer(KIND=JWIM), parameter :: rtype = MPI_REAL8
#endif
#endif

end module yowDatapool
