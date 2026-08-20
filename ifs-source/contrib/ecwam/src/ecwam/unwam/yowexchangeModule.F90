! (C) Copyright 2001- Aron Roland (Roland & Partner, Germany).
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

#include "yowincludes.h"

!> Has only the ghost nodes assign to a neighbor domain
module yowExchangeModule
  USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
  use yowMpiModule
  implicit none
  private
  public :: initNbrDomains, exchange, createMPITypes, setDimSize, finalizeExchangeModule

  !> Holds some data belong to a neighbor Domain
  type, public :: t_neighborDomain
    
    
    !> the domain ID
    !> The domain ID of the neighbor domain. Starts by 1
    integer(KIND=JWIM) :: domainID = 0

    !> number of ghosts nodes.
    !> holds the number of ghosts nodes the two domains share together
    !> (the neighbor domain has a copy of the ghosts. when you change some
    !> value in a ghost node, it dosen't change in the node on the other domain)
    integer(KIND=JWIM) :: numNodesToReceive = 0

    !> this are the ghosts that we.
    !> has in this neighbor domain. global node IDs
    integer(KIND=JWIM), allocatable :: nodesToReceive(:)

    !> number of nodes we have to send to this neighbor
    integer(KIND=JWIM) :: numNodesToSend = 0

    !> this are the ghosts from this neighbor.
    !> global node IDs to send
    integer(KIND=JWIM), allocatable :: nodesToSend(:)

#ifdef WAM_HAVE_MPI_F08
    !> MPI datatypes for 1D exchange
    type(mpi_datatype) :: p1DRsendType = MPI_DATATYPE_NULL
    type(mpi_datatype) :: p1DRrecvType = MPI_DATATYPE_NULL
    type(mpi_datatype) :: p1DIsendType = MPI_DATATYPE_NULL
    type(mpi_datatype) :: p1DIrecvType = MPI_DATATYPE_NULL
    !> MPI datatypes for 2D exchange
    type(mpi_datatype) :: p2DRsendType1 = MPI_DATATYPE_NULL
    type(mpi_datatype) :: p2DRrecvType1 = MPI_DATATYPE_NULL
    type(mpi_datatype) :: p2DRsendType2 = MPI_DATATYPE_NULL
    type(mpi_datatype) :: p2DRrecvType2 = MPI_DATATYPE_NULL
    !> MPI datatypes for 3D exchange
    type(mpi_datatype) :: p3DR8sendType = MPI_DATATYPE_NULL
    type(mpi_datatype) :: p3DR8recvType = MPI_DATATYPE_NULL
    type(mpi_datatype) :: p3DR4sendType = MPI_DATATYPE_NULL
    type(mpi_datatype) :: p3DR4recvType = MPI_DATATYPE_NULL
#else
    !> MPI datatypes for 1D exchange
    integer(KIND=JWIM) :: p1DRsendType = MPI_DATATYPE_NULL
    integer(KIND=JWIM) :: p1DRrecvType = MPI_DATATYPE_NULL
    integer(KIND=JWIM) :: p1DIsendType = MPI_DATATYPE_NULL
    integer(KIND=JWIM) :: p1DIrecvType = MPI_DATATYPE_NULL
    !> MPI datatypes for 2D exchange
    integer(KIND=JWIM) :: p2DRsendType1 = MPI_DATATYPE_NULL
    integer(KIND=JWIM) :: p2DRrecvType1 = MPI_DATATYPE_NULL
    integer(KIND=JWIM) :: p2DRsendType2 = MPI_DATATYPE_NULL
    integer(KIND=JWIM) :: p2DRrecvType2 = MPI_DATATYPE_NULL
    !> MPI datatypes for 3D exchange
    integer(KIND=JWIM) :: p3DR8sendType = MPI_DATATYPE_NULL
    integer(KIND=JWIM) :: p3DR8recvType = MPI_DATATYPE_NULL
    integer(KIND=JWIM) :: p3DR4sendType = MPI_DATATYPE_NULL
    integer(KIND=JWIM) :: p3DR4recvType = MPI_DATATYPE_NULL
#endif
    contains
!     procedure :: exchangeGhostIds
!     final :: finalizeNeighborDomain
    procedure :: finalize
    procedure :: createMPIType

  end type

  !> Knows for all domains neighbors, which node we must send or revc from neighbor domains
  !> from 1 to nConnDomains
  type(t_neighborDomain), public, allocatable :: neighborDomains(:)

  !> Number of neighbor domains
  integer(KIND=JWIM), public :: nConnDomains = 0

  !> number of the second dimension for exchange
  integer(KIND=JWIM), public :: n2ndDim = 1

  !> number fo the third dimension for exchange
  integer(KIND=JWIM), public :: n3ndDim = 1

  interface exchange
    module procedure exchange1Dreal, exchange1Dint, exchange2Dreal, exchange3Dreal4, exchange3Dreal8
  end interface

  interface nodeExchange
    module procedure nodeExchange1Dreal4
  end interface

  contains


  subroutine finalize(this)
    use yowError
    use yowMpiModule
    implicit none
    class(t_neighborDomain), intent(inout) :: this
    integer(KIND=JWIM) :: ierr

    if(allocated(this%nodesToSend))    deallocate(this%nodesToSend)
    if(allocated(this%nodesToReceive)) deallocate(this%nodesToReceive)

    call mpi_type_free(this%p1DRsendType, ierr)
    if(ierr /= MPI_SUCCESS) PARALLEL_ABORT("freeMPItype", ierr)
    call mpi_type_free(this%p1DRrecvType, ierr)
    if(ierr /= MPI_SUCCESS) PARALLEL_ABORT("freeMPItype", ierr)
    call mpi_type_free(this%p1DIsendType, ierr)
    if(ierr /= MPI_SUCCESS) PARALLEL_ABORT("freeMPItype", ierr)
    call mpi_type_free(this%p1DIrecvType, ierr)
    if(ierr /= MPI_SUCCESS) PARALLEL_ABORT("freeMPItype", ierr)
    call mpi_type_free(this%p2DRsendType1, ierr)
    if(ierr /= MPI_SUCCESS) PARALLEL_ABORT("freeMPItype", ierr)
    call mpi_type_free(this%p2DRrecvType1, ierr)
    if(ierr /= MPI_SUCCESS) PARALLEL_ABORT("freeMPItype", ierr)
    call mpi_type_free(this%p2DRsendType2, ierr)
    if(ierr /= MPI_SUCCESS) PARALLEL_ABORT("freeMPItype", ierr)
    call mpi_type_free(this%p2DRrecvType2, ierr)
    if(ierr /= MPI_SUCCESS) PARALLEL_ABORT("freeMPItype", ierr)
    call mpi_type_free(this%p3DR4sendType, ierr)
    if(ierr /= MPI_SUCCESS) PARALLEL_ABORT("freeMPItype", ierr)
    call mpi_type_free(this%p3DR4recvType, ierr)
    if(ierr /= MPI_SUCCESS) PARALLEL_ABORT("freeMPItype", ierr)
    call mpi_type_free(this%p3DR8sendType, ierr)
    if(ierr /= MPI_SUCCESS) PARALLEL_ABORT("freeMPItype", ierr)
    call mpi_type_free(this%p3DR8recvType, ierr)
    if(ierr /= MPI_SUCCESS) PARALLEL_ABORT("freeMPItype", ierr)
  end subroutine

  ! create MPI indexed datatype for this neighborDomain
  subroutine createMPIType(this)
    use yowError
    use yowMpiModule
    use yowNodepool, only: ghostgl, np, ipgl
    use yowDatapool, only: rtype, itype
    implicit none
    class(t_neighborDomain), intent(inout) :: this

    integer(KIND=JWIM) :: ierr
    integer(KIND=JWIM) :: dsplSend(this%numNodesToSend)
    integer(KIND=JWIM) :: dsplRecv(this%numNodesToReceive)

    
    dsplSend = ipgl(this%nodesToSend)-1
    dsplRecv = ghostgl(this%nodesToReceive) + np -1

    ! p1D real
    call mpi_type_create_indexed_block(this%numNodesToSend, 1, dsplSend, rtype, this%p1DRsendType,ierr)
    if(ierr /= MPI_SUCCESS) PARALLEL_ABORT("createMPIType", ierr)
    call mpi_type_commit(this%p1DRsendType,ierr)
    if(ierr /= MPI_SUCCESS) PARALLEL_ABORT("createMPIType", ierr)
    
    call mpi_type_create_indexed_block(this%numNodesToReceive, 1, dsplRecv, rtype, this%p1DRrecvType,ierr)
    if(ierr /= MPI_SUCCESS) PARALLEL_ABORT("createMPIType", ierr)
    call mpi_type_commit(this%p1DRrecvType,ierr)
    if(ierr /= MPI_SUCCESS) PARALLEL_ABORT("createMPIType", ierr)

    ! p1D integer
    call mpi_type_create_indexed_block(this%numNodesToSend, 1, dsplSend, itype, this%p1DIsendType,ierr)
    if(ierr /= MPI_SUCCESS) PARALLEL_ABORT("createMPIType", ierr)
    call mpi_type_commit(this%p1DIsendType,ierr)
    if(ierr /= MPI_SUCCESS) PARALLEL_ABORT("createMPIType", ierr)

    call mpi_type_create_indexed_block(this%numNodesToReceive, 1, dsplRecv, itype, this%p1DIrecvType,ierr)
    if(ierr /= MPI_SUCCESS) PARALLEL_ABORT("createMPIType", ierr)
    call mpi_type_commit(this%p1DIrecvType,ierr)
    if(ierr /= MPI_SUCCESS) PARALLEL_ABORT("createMPIType", ierr)

    ! p2D real
    dsplSend = (ipgl(this%nodesToSend)-1) * n2ndDim
    dsplRecv = (ghostgl(this%nodesToReceive) + np -1) * n2ndDim
    call mpi_type_create_indexed_block(this%numNodesToSend, n2ndDim, dsplSend, rtype, this%p2DRsendType1,ierr)
    if(ierr /= MPI_SUCCESS) PARALLEL_ABORT("createMPIType", ierr)
    call mpi_type_commit(this%p2DRsendType1,ierr)
    if(ierr /= MPI_SUCCESS) PARALLEL_ABORT("createMPIType", ierr)

    call mpi_type_create_indexed_block(this%numNodesToReceive, n2ndDim, dsplRecv, rtype, this%p2DRrecvType1,ierr)
    if(ierr /= MPI_SUCCESS) PARALLEL_ABORT("createMPIType", ierr)
    call mpi_type_commit(this%p2DRrecvType1,ierr)
    if(ierr /= MPI_SUCCESS) PARALLEL_ABORT("createMPIType", ierr)


    dsplSend = (ipgl(this%nodesToSend)-1) * n3ndDim
    dsplRecv = (ghostgl(this%nodesToReceive) + np -1) * n3ndDim
    call mpi_type_create_indexed_block(this%numNodesToSend, n3ndDim, dsplSend, rtype, this%p2DRsendType2,ierr)
    if(ierr /= MPI_SUCCESS) PARALLEL_ABORT("createMPIType", ierr)
    call mpi_type_commit(this%p2DRsendType2,ierr)
    if(ierr /= MPI_SUCCESS) PARALLEL_ABORT("createMPIType", ierr)

    call mpi_type_create_indexed_block(this%numNodesToReceive, n3ndDim, dsplRecv, rtype, this%p2DRrecvType2,ierr)
    if(ierr /= MPI_SUCCESS) PARALLEL_ABORT("createMPIType", ierr)
    call mpi_type_commit(this%p2DRrecvType2,ierr)
    if(ierr /= MPI_SUCCESS) PARALLEL_ABORT("createMPIType", ierr)

    ! p3D real4
    dsplSend = (ipgl(this%nodesToSend)-1) * n2ndDim*n3ndDim
    dsplRecv = (ghostgl(this%nodesToReceive) + np -1) * n2ndDim*n3ndDim
    call mpi_type_create_indexed_block(this%numNodesToSend, n2ndDim*n3ndDim, dsplSend, MPI_REAL4, this%p3DR4sendType,ierr)
    if(ierr /= MPI_SUCCESS) PARALLEL_ABORT("createMPIType", ierr)
    call mpi_type_commit(this%p3DR4sendType,ierr)
    if(ierr /= MPI_SUCCESS) PARALLEL_ABORT("createMPIType", ierr)

    call mpi_type_create_indexed_block(this%numNodesToReceive, n2ndDim*n3ndDim, dsplRecv, MPI_REAL4, this%p3DR4recvType,ierr)
    if(ierr /= MPI_SUCCESS) PARALLEL_ABORT("createMPIType", ierr)
    call mpi_type_commit(this%p3DR4recvType,ierr)
    if(ierr /= MPI_SUCCESS) PARALLEL_ABORT("createMPIType", ierr)

    ! p3D real8
    dsplSend = (ipgl(this%nodesToSend)-1) * n2ndDim*n3ndDim
    dsplRecv = (ghostgl(this%nodesToReceive) + np -1) * n2ndDim*n3ndDim
    call mpi_type_create_indexed_block(this%numNodesToSend, n2ndDim*n3ndDim, dsplSend, MPI_REAL8, this%p3DR8sendType,ierr)
    if(ierr /= MPI_SUCCESS) PARALLEL_ABORT("createMPIType", ierr)
    call mpi_type_commit(this%p3DR8sendType,ierr)
    if(ierr /= MPI_SUCCESS) PARALLEL_ABORT("createMPIType", ierr)

    call mpi_type_create_indexed_block(this%numNodesToReceive, n2ndDim*n3ndDim, dsplRecv, MPI_REAL8, this%p3DR8recvType,ierr)
    if(ierr /= MPI_SUCCESS) PARALLEL_ABORT("createMPIType", ierr)
    call mpi_type_commit(this%p3DR8recvType,ierr)
    if(ierr /= MPI_SUCCESS) PARALLEL_ABORT("createMPIType", ierr)
  end subroutine

  subroutine initNbrDomains(nConnD)
    use yowError
    implicit none
    integer(KIND=JWIM), intent(in) :: nConnD

    call finalizeExchangeModule()
    nConnDomains = nConnD
    allocate(neighborDomains(nConnDomains), stat=stat)
    if(stat/=0)  ABORT('neighborDomains allocation failure')
  end subroutine

  subroutine createMPITypes()    
    implicit none
    integer(KIND=JWIM) :: i

    do i=1, nConnDomains
      call neighborDomains(i)%createMPIType()
    end do
  end subroutine

  !> exchange values in U.
  !> \param[inout] U array with values to exchange. np+ng long.
  !> Send values from U(1:np) to other threads.
  !> Receive values from other threads and updates U(np+1:np+ng
  !> \note MPI recv tag: 10000 + MPI rank
  !> \note MPI send tag: 10000 + neighbor MPI rank
  subroutine exchange1Dreal(U)
    use yowDatapool, only: comm, myrank, rkind
    use yowNodepool, only: npa
    use yowError
    use yowMpiModule
    implicit none
    real(kind=rkind), intent(inout) :: U(:)

    integer(KIND=JWIM) :: i, ierr, tag
#ifdef WAM_HAVE_MPI_F08
    type(mpi_request) :: sendRqst(nConnDomains), recvRqst(nConnDomains)
    type(mpi_status)  :: recvStat(nConnDomains), sendStat(nConnDomains)
#else
    integer(KIND=JWIM) :: sendRqst(nConnDomains), recvRqst(nConnDomains)
    integer(KIND=JWIM) :: recvStat(MPI_STATUS_SIZE, nConnDomains), sendStat(MPI_STATUS_SIZE, nConnDomains)
#endif
    character(len=200) errstr

    if(size(U) /= npa) then
      write(errstr, *) "size(U) /= npa", size(U), "should be", npa
      ABORT(errstr)
    endif

    ! post receives
    do i=1, nConnDomains
      tag = 10000 + myrank
      call MPI_IRecv(U, &
                    1, &
                    neighborDomains(i)%p1DRrecvType, &
                    neighborDomains(i)%domainID-1, &
                    tag, &
                    comm, &
                    recvRqst(i), &
                    ierr)
      if(ierr/=MPI_SUCCESS) then
        PARALLEL_ABORT("MPI_IRecv", ierr)
      endif
    enddo

    ! post sends
    do i=1, nConnDomains
      tag = 10000 + (neighborDomains(i)%domainID-1)
      call MPI_ISend(U, &
                    1, &
                    neighborDomains(i)%p1DRsendType, &
                    neighborDomains(i)%domainID-1, &
                    tag, &
                    comm, &
                    sendRqst(i), &
                    ierr);
        if(ierr/=MPI_SUCCESS) then
          PARALLEL_ABORT("MPI_ISend", ierr)
        endif
    end do

    ! Wait for completion
    call mpi_waitall(nConnDomains, recvRqst, recvStat,ierr)
    if(ierr/=MPI_SUCCESS) PARALLEL_ABORT("waitall", ierr)
    call mpi_waitall(nConnDomains, sendRqst, sendStat,ierr)
    if(ierr/=MPI_SUCCESS) PARALLEL_ABORT("waitall", ierr)
  end subroutine exchange1Dreal

  !> \overload exchange1Dreal
  !> \note MPI recv tag: 20000 + MPI rank
  !> \note MPI send tag: 20000 + neighbor MPI rank
  subroutine exchange1Dint(U)
    use yowDatapool, only: comm, myrank
    use yowNodepool, only: npa
    use yowError
    use yowMpiModule
    implicit none
    integer(KIND=JWIM), intent(inout) :: U(:)

    integer(KIND=JWIM) :: i, ierr, tag
#ifdef WAM_HAVE_MPI_F08
    type(mpi_request) :: sendRqst(nConnDomains), recvRqst(nConnDomains)
    type(mpi_status)  :: recvStat(nConnDomains), sendStat(nConnDomains)
#else
    integer(KIND=JWIM) :: sendRqst(nConnDomains), recvRqst(nConnDomains)
    integer(KIND=JWIM) :: recvStat(MPI_STATUS_SIZE, nConnDomains), sendStat(MPI_STATUS_SIZE, nConnDomains)
#endif
    character(len=200) errstr

    if(size(U) /= npa) then
      write(errstr, *) "size(U) /= npa", size(U), "should be", npa
      ABORT(errstr)
    endif

    ! post receives
    do i=1, nConnDomains
      tag = 20000 + myrank
      call MPI_IRecv(U, &
                    1, &
                    neighborDomains(i)%p1DIrecvType, &
                    neighborDomains(i)%domainID-1, &
                    tag, &
                    comm, &
                    recvRqst(i), &
                    ierr)
      if(ierr/=MPI_SUCCESS) then
        PARALLEL_ABORT("MPI_IRecv", ierr)
      endif
    enddo

    ! post sends
    do i=1, nConnDomains
      tag = 20000 + (neighborDomains(i)%domainID-1)
      call MPI_ISend(U, &
                    1, &
                    neighborDomains(i)%p1DIsendType, &
                    neighborDomains(i)%domainID-1, &
                    tag, &
                    comm, &
                    sendRqst(i), &
                    ierr);
        if(ierr/=MPI_SUCCESS) then
          PARALLEL_ABORT("MPI_ISend", ierr)
        endif
    end do

    ! Wait for completion
    call mpi_waitall(nConnDomains, recvRqst, recvStat,ierr)
    if(ierr/=MPI_SUCCESS) PARALLEL_ABORT("waitall", ierr)
    call mpi_waitall(nConnDomains, sendRqst, sendStat,ierr)
    if(ierr/=MPI_SUCCESS) PARALLEL_ABORT("waitall", ierr)
  end subroutine exchange1Dint

  !> \overload exchange1Dreal
  !> 
  !> \note MPI recv tag: 30000 + MPI rank
  !> \note MPI send tag: 30000 + neighbor MPI rank
  subroutine exchange2Dreal(U)
    use yowDatapool, only: comm, myrank, rkind
    use yowNodepool, only: npa
    use yowError
    use yowMpiModule
    implicit none
    real(kind=rkind), intent(inout) :: U(:,:)

    integer(KIND=JWIM) :: i, ierr, tag
#ifdef WAM_HAVE_MPI_F08
    type(mpi_request) :: sendRqst(nConnDomains), recvRqst(nConnDomains)
    type(mpi_status)  :: recvStat(nConnDomains), sendStat(nConnDomains)
#else
    integer(KIND=JWIM) :: sendRqst(nConnDomains), recvRqst(nConnDomains)
    integer(KIND=JWIM) :: recvStat(MPI_STATUS_SIZE, nConnDomains), sendStat(MPI_STATUS_SIZE, nConnDomains)
#endif
    character(len=200) errstr

    if(size(U,2) /= npa) then
      write(errstr, *) "size(U,2) /= npa", size(U,2), "should be", npa
      ABORT(errstr)
    endif
 
    if((size(U,1) /= n2ndDim) ) then
      if((size(U,1) /= n3ndDim) ) then
        write(errstr, *) "size(U,1) /= n2ndDim or n3ndDim. size(U,1)=", size(U,1), " n2ndDim=", n2ndDim, " n3ndDim=", n3ndDim
        ABORT(errstr)
      endif
    endif

    !> \todo do that better
    if(size(U,1) == n2ndDim) then
      ! post receives
      do i=1, nConnDomains
        tag = 30000 + myrank
        call MPI_IRecv(U, &
                      1, &
                      neighborDomains(i)%p2DRrecvType1, &
                      neighborDomains(i)%domainID-1, &
                      tag, &
                      comm, &
                      recvRqst(i), &
                      ierr)
        if(ierr/=MPI_SUCCESS) then
          PARALLEL_ABORT("MPI_IRecv", ierr)
        endif
      enddo

      ! post sends
      do i=1, nConnDomains
        tag = 30000 + (neighborDomains(i)%domainID-1)
        call MPI_ISend(U, &
                      1, &
                      neighborDomains(i)%p2DRsendType1, &
                      neighborDomains(i)%domainID-1, &
                      tag, &
                      comm, &
                      sendRqst(i), &
                      ierr);
          if(ierr/=MPI_SUCCESS) then
            PARALLEL_ABORT("MPI_ISend", ierr)
          endif
      end do
    else if(size(U,1) == n3ndDim) then
      ! post receives
      do i=1, nConnDomains
        tag = 30000 + myrank
        call MPI_IRecv(U, &
                      1, &
                      neighborDomains(i)%p2DRrecvType2, &
                      neighborDomains(i)%domainID-1, &
                      tag, &
                      comm, &
                      recvRqst(i), &
                      ierr)
        if(ierr/=MPI_SUCCESS) then
          PARALLEL_ABORT("MPI_IRecv", ierr)
        endif
      enddo

      ! post sends
      do i=1, nConnDomains
        tag = 30000 + (neighborDomains(i)%domainID-1)
        call MPI_ISend(U, &
                      1, &
                      neighborDomains(i)%p2DRsendType2, &
                      neighborDomains(i)%domainID-1, &
                      tag, &
                      comm, &
                      sendRqst(i), &
                      ierr);
          if(ierr/=MPI_SUCCESS) then
            PARALLEL_ABORT("MPI_ISend", ierr)
          endif
      end do
    endif

    ! Wait for completion
    call mpi_waitall(nConnDomains, recvRqst, recvStat,ierr)
    if(ierr/=MPI_SUCCESS) PARALLEL_ABORT("waitall", ierr)
    call mpi_waitall(nConnDomains, sendRqst, sendStat,ierr)
    if(ierr/=MPI_SUCCESS) PARALLEL_ABORT("waitall", ierr)
  end subroutine exchange2Dreal

  !> \overload exchange1Dreal
  !> \note MPI recv tag: 40000 + MPI rank
  !> \note MPI send tag: 40000 + neighbor MPI rank
  subroutine exchange3Dreal8(U)
    use yowDatapool, only: comm, myrank, rkind
    use yowNodepool, only: npa
    use yowError
    use yowMpiModule
    implicit none
    real(kind=8), intent(inout) :: U(:,:,:)

    integer(KIND=JWIM) :: i, ierr, tag
#ifdef WAM_HAVE_MPI_F08
    type(mpi_request) :: sendRqst(nConnDomains), recvRqst(nConnDomains)
    type(mpi_status)  :: recvStat(nConnDomains), sendStat(nConnDomains)
#else
    integer(KIND=JWIM) :: sendRqst(nConnDomains), recvRqst(nConnDomains)
    integer(KIND=JWIM) :: recvStat(MPI_STATUS_SIZE, nConnDomains), sendStat(MPI_STATUS_SIZE, nConnDomains)
#endif
    character(len=200) errstr

    if(size(U,3) /= npa) then
      write(errstr, *) "size(U,3) /= npa", size(U), "should be", npa
      ABORT(errstr)
    endif

    if((size(U,2) /= n2ndDim) ) then
      write(errstr, *) "size(U,2) /= n2ndDim. size(U,2)=", size(U,2), " n2ndDim=", n2ndDim
      ABORT(errstr)
    endif

    if((size(U,1) /= n3ndDim) ) then
      write(errstr, *) "size(U,1) /= n3ndDim. size(U,1)=", size(U,1), " n2ndDim=", n3ndDim
      ABORT(errstr)
    endif

    ! post receives
    do i=1, nConnDomains
      tag = 40000 + myrank
      call MPI_IRecv(U, &
                    1, &
                    neighborDomains(i)%p3DR8recvType, &
                    neighborDomains(i)%domainID-1, &
                    tag, &
                    comm, &
                    recvRqst(i), &
                    ierr)
      if(ierr/=MPI_SUCCESS) then
        PARALLEL_ABORT("MPI_IRecv", ierr)
      endif
    enddo

    ! post sends
    do i=1, nConnDomains
      tag = 40000 + (neighborDomains(i)%domainID-1)
      call MPI_ISend(U, &
                    1, &
                    neighborDomains(i)%p3DR8sendType, &
                    neighborDomains(i)%domainID-1, &
                    tag, &
                    comm, &
                    sendRqst(i), &
                    ierr);
        if(ierr/=MPI_SUCCESS) then
          PARALLEL_ABORT("MPI_ISend", ierr)
        endif
    end do

    ! Wait for completion
    call mpi_waitall(nConnDomains, recvRqst, recvStat,ierr)
    if(ierr/=MPI_SUCCESS) PARALLEL_ABORT("waitall", ierr)
    call mpi_waitall(nConnDomains, sendRqst, sendStat,ierr)
    if(ierr/=MPI_SUCCESS) PARALLEL_ABORT("waitall", ierr)
  end subroutine exchange3Dreal8

  !> \overload exchange1Dreal
  !> \note MPI recv tag: 40000 + MPI rank
  !> \note MPI send tag: 40000 + neighbor MPI rank
  subroutine exchange3Dreal4(U)
    use yowDatapool, only: comm, myrank, rkind
    use yowNodepool, only: npa
    use yowError
    use yowMpiModule
    implicit none
    real(kind=4), intent(inout) :: U(:,:,:)

    integer(KIND=JWIM) :: i, ierr, tag
#ifdef WAM_HAVE_MPI_F08
    type(mpi_request) :: sendRqst(nConnDomains), recvRqst(nConnDomains)
    type(mpi_status)  :: recvStat(nConnDomains), sendStat(nConnDomains)
#else
    integer(KIND=JWIM) :: sendRqst(nConnDomains), recvRqst(nConnDomains)
    integer(KIND=JWIM) :: recvStat(MPI_STATUS_SIZE, nConnDomains), sendStat(MPI_STATUS_SIZE, nConnDomains)
#endif
    character(len=200) errstr

    if(size(U,3) /= npa) then
      write(errstr, *) "size(U,3) /= npa", size(U), "should be", npa
      ABORT(errstr)
    endif

    if((size(U,2) /= n2ndDim) ) then
      ABORT("size(U,2) /= n2ndDim")
    endif

    if((size(U,1) /= n3ndDim) ) then
      ABORT("size(U,1) /= n3ndDim")
    endif

    ! post receives
    do i=1, nConnDomains
      tag = 40000 + myrank
      call MPI_IRecv(U, &
                    1, &
                    neighborDomains(i)%p3DR4recvType, &
                    neighborDomains(i)%domainID-1, &
                    tag, &
                    comm, &
                    recvRqst(i), &
                    ierr)
      if(ierr/=MPI_SUCCESS) then
        PARALLEL_ABORT("MPI_IRecv", ierr)
      endif
    enddo

    ! post sends
    do i=1, nConnDomains
      tag = 40000 + (neighborDomains(i)%domainID-1)
      call MPI_ISend(U, &
                    1, &
                    neighborDomains(i)%p3DR4sendType, &
                    neighborDomains(i)%domainID-1, &
                    tag, &
                    comm, &
                    sendRqst(i), &
                    ierr);
        if(ierr/=MPI_SUCCESS) then
          PARALLEL_ABORT("MPI_ISend", ierr)
        endif
    end do

    ! Wait for completion
    call mpi_waitall(nConnDomains, recvRqst, recvStat,ierr)
    if(ierr/=MPI_SUCCESS) PARALLEL_ABORT("waitall", ierr)
    call mpi_waitall(nConnDomains, sendRqst, sendStat,ierr)
    if(ierr/=MPI_SUCCESS) PARALLEL_ABORT("waitall", ierr)
  end subroutine exchange3Dreal4

  
  !> set the size of the second and third dimension for exchange
  !> \note the size of the first dimension is npa
  !> \note call this before initPD()
  subroutine setDimSize(second, third)
    implicit none
    integer(KIND=JWIM), intent(in) :: second, third
    n2ndDim = second
    n3ndDim = third
  end subroutine setDimSize

  subroutine finalizeExchangeModule()
    implicit none
    integer(KIND=JWIM) :: i

    if(allocated(neighborDomains)) then
      do i=1, size(neighborDomains)
        call neighborDomains(i)%finalize()
      end do
      deallocate(neighborDomains)
    endif
  end subroutine

  !>
  !> The Array \a values is as long as array \a globalIDs.
  !> Return in array \a values the values of \a U from the other thread
  !> which are the owner of the given global node IDs in array \a globalIDs
  !>
  subroutine nodeExchange1Dreal4(U, globalIDs, values)
    use yowError
    use yowNodepool, only: npa, node2rank, np_global
    use yowDatapool, only: nTasks
    implicit none
    real(kind=4), intent(in) :: U(:)
    integer(KIND=JWIM), intent(in) :: globalIDs(:)
    real(kind=4), intent(inout) :: values(:)

    character(len=200) errstr
    integer(KIND=JWIM) :: i, rank
    integer(KIND=JWIM) :: nToGetPerRank(nTasks)

    if(size(U) /= npa) then
      write(errstr, *) "size(U) /= npa", size(U), "should be", npa
      ABORT(errstr)
    endif
    
    if(size(values) < size(globalIDs) ) then
      write(errstr, *) "size(values) < size(globalIds) need at least", size(globalIDs)
      ABORT(errstr)
    endif

    if(maxval(globalIDs) > np_global) then
      write(errstr, *) "globalIDs contains a value which is grather then np_global", maxval(globalIDs)
      ABORT(errstr)
    endif

    if(minval(globalIDs) < 1) then
      write(errstr, *) "globalIDs contains a value which is less than 1", minval(globalIDs)
      ABORT(errstr)
    endif

    ! detect the rank of every global ID
    ! send to every rand from which we want some values, to gloabl IDs
    ! receive the global IDs
    ! get the values from U
    ! send back

    ! what if one thread does not receive anything but calls mpi receice?
    ! aaaahhh, we send a "0" to every thread in the world. "nothing to receive" haha :D

    nToGetPerRank(:) = 0
    do i=1, size(globalIDs)
      rank = node2rank(globalIDs(i))
      nToGetPerRank(rank) = nToGetPerRank(rank)
    end do

    do i=1, nTasks
      ! receive. fill nToSendPerRank
    
      ! send task i-1 we want nToGetPerRank(i) nodes from it
    end do

    ! todo make some type equal to t_neighborDomain which store how many and which nodes we have to send
    ! to which rank
    
    
  end subroutine
end module yowExchangeModule
