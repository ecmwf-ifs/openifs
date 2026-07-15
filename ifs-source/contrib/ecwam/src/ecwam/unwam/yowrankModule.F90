! (C) Copyright 2001- Aron Roland (Roland & Partner, Germany).
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

!> \file yowrankModule.F90
!> \brief provides RANK array
!> \author Thomas Huxhorn
!> \date 2013

#include "yowincludes.h"

!> Provides access to some information of all threads e.g. iplg
module yowRankModule
  USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
  use yowError
  implicit none
  private
  public :: initRankModule, finalizeRankModule

  type, public :: t_rank
    !> number of local nodes of this rank
    integer(KIND=JWIM) :: np = 0

    !> nummer of ghost + resident nodes of this rank
    integer(KIND=JWIM) :: npa = 0

    !> Node local to gloabl mapping of this rank
    !> rank%npa long
    integer(KIND=JWIM), allocatable :: iplg(:)

    !> global start node number for every thread
    integer(KIND=JWIM):: IStart = 0
  end type

  !> Provides access to some information of all threads e.g. iplg
  !> \note range [1:nTasks]
  !> \note range[myrank] are filled with the values from this rank
  type(t_rank),public, allocatable :: rank(:)

  contains

  !> allocate and exchange
  subroutine initRankModule()
    use yowDatapool, only: nTasks
    implicit none
    integer(KIND=JWIM) :: stat

    if(allocated(rank)) deallocate(rank)
    allocate(rank(nTasks), stat=stat)
    if(stat/=0) ABORT('rank allocation failure')

    call exchangeIPLG()
    call calcISTART()
  end subroutine

  !> send iplg from this thread to every neighbor thread
  !> \internal
  subroutine exchangeIPLG()
    use yowNodepool, only: np, npa, iplg
    use yowDatapool, only: nTasks, myrank, comm, itype
    use yowMpiModule
    implicit none
    integer(KIND=JWIM) :: i, ierr, stat
#ifdef WAM_HAVE_MPI_F08
    type(mpi_request)  :: sendRqst(nTasks), recvRqst(nTasks)
    type(mpi_status)   :: recvStat(nTasks), sendStat(nTasks)
#else
    integer(KIND=JWIM) :: sendRqst(nTasks), recvRqst(nTasks)
    integer(KIND=JWIM) :: recvStat(MPI_STATUS_SIZE, nTasks), sendStat(MPI_STATUS_SIZE, nTasks)
#endif

    ! step1 exchange np
    ! step2 exchange npa
    ! step3 allocate rank%iplg
    ! step4 exchange iplg

    ! step1 exchange np
    ! post receives
    do i=1, nTasks
      if(i /= myrank+1) then
        call MPI_IRecv(rank(i)%np, &
                       1, &
                       itype, &
                       i-1, & ! rank
                       42, &
                       comm, &
                       recvRqst(i), &
                       ierr)
        if(ierr/=MPI_SUCCESS) then
          PARALLEL_ABORT("MPI_IRecv", ierr)
        endif
      else
        recvRqst(i) = MPI_REQUEST_NULL
      endif
    end do

    ! post sends
    do i=1, nTasks
      if(i /= myrank+1) then
        call MPI_ISend(np, &
                       1, &
                       itype, &
                       i-1, & ! rank
                       42, &
                       comm, &
                       sendRqst(i), &
                       ierr);
        if(ierr/=MPI_SUCCESS) then
          PARALLEL_ABORT("MPI_ISend", ierr)
        endif
      else
        sendRqst(i) = MPI_REQUEST_NULL
      endif
    end do

    rank(myrank+1)%np = np

    ! Wait for completion
    call mpi_waitall(nTasks, recvRqst, recvStat,ierr)
    if(ierr/=MPI_SUCCESS) PARALLEL_ABORT("waitall", ierr)
    call mpi_waitall(nTasks, sendRqst, sendStat,ierr)
    if(ierr/=MPI_SUCCESS) PARALLEL_ABORT("waitall", ierr)

    ! step2 exchange npa
    ! post receives
    do i=1, nTasks
      if(i /= myrank+1) then
        call MPI_IRecv(rank(i)%npa, &
                       1, &
                       itype, &
                       i-1, & ! rank
                       42, &
                       comm, &
                       recvRqst(i), &
                       ierr)
        if(ierr/=MPI_SUCCESS) then
          PARALLEL_ABORT("MPI_IRecv", ierr)
        endif
      else
        recvRqst(i) = MPI_REQUEST_NULL
      endif
    end do

    ! post sends
    do i=1, nTasks
      if(i /= myrank+1) then
        call MPI_ISend(npa, &
                       1, &
                       itype, &
                       i-1, & ! rank
                       42, &
                       comm, &
                       sendRqst(i), &
                       ierr);
        if(ierr/=MPI_SUCCESS) then
          PARALLEL_ABORT("MPI_ISend", ierr)
        endif
      else
        sendRqst(i) = MPI_REQUEST_NULL
      endif
    end do

    rank(myrank+1)%npa = npa

    ! Wait for completion
    call mpi_waitall(nTasks, recvRqst, recvStat,ierr)
    if(ierr/=MPI_SUCCESS) PARALLEL_ABORT("waitall", ierr)
    call mpi_waitall(nTasks, sendRqst, sendStat,ierr)
    if(ierr/=MPI_SUCCESS) PARALLEL_ABORT("waitall", ierr)

    ! step3 allocal rank%iplg
    do i=1, nTasks
      if(allocated(rank(i)%iplg)) deallocate(rank(i)%iplg)
      allocate(rank(i)%iplg(rank(i)%npa), stat=stat)
      if(stat/=0) ABORT('rank%iplg allocation failure')
      rank(i)%iplg = 0
    end do



    ! step4 exchange iplg
    ! post receives
    do i=1, nTasks
      if(i /= myrank+1) then
        call MPI_IRecv(rank(i)%iplg, &
                       rank(i)%npa, &
                       itype, &
                       i-1, & ! rank
                       42, &
                       comm, &
                       recvRqst(i), &
                       ierr)
        if(ierr/=MPI_SUCCESS) then
          PARALLEL_ABORT("MPI_IRecv", ierr)
        endif
      else
        recvRqst(i) = MPI_REQUEST_NULL
      endif
    end do

    ! post sends
    do i=1, nTasks
      if(i /= myrank+1) then
        call MPI_ISend(iplg, &
                       npa, &
                       itype, &
                       i-1, & ! rank
                       42, &
                       comm, &
                       sendRqst(i), &
                       ierr);
        if(ierr/=MPI_SUCCESS) then
          PARALLEL_ABORT("MPI_ISend", ierr)
        endif
      else
        sendRqst(i) = MPI_REQUEST_NULL
      endif
    end do

    rank(myrank+1)%iplg = iplg

    ! Wait for completion
    call mpi_waitall(nTasks, recvRqst, recvStat,ierr)
    if(ierr/=MPI_SUCCESS) PARALLEL_ABORT("waitall", ierr)
    call mpi_waitall(nTasks, sendRqst, sendStat,ierr)
    if(ierr/=MPI_SUCCESS) PARALLEL_ABORT("waitall", ierr)
  end subroutine

  !> \internal
  subroutine calcISTART()
    use yowDatapool, only: nTasks
    implicit none
    integer(KIND=JWIM) :: ir

    rank(1)%IStart = 1
    do ir=2, nTasks
      rank(ir)%IStart = rank(ir-1)%IStart + rank(ir-1)%np
    end do
  end subroutine

  subroutine finalizeRankModule()
    implicit none
    integer(KIND=JWIM) :: i
  
    if(allocated(rank)) then
      do i=1, size(rank)
        if(allocated(rank(i)%iplg)) deallocate(rank(i)%iplg)
      end do
      deallocate(rank)
    endif
  end subroutine
end module
