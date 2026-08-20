! (C) Copyright 2001- Aron Roland (Roland & Partner, Germany).
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

!> \file yowchecks.F90
!> \brief Some very usefull checks to test the partition and MPI
!> \author Mathieu Dutour, Thomas Huxhorn
!> \data 2014

#include "yowincludes.h"

module yowChecksModule  
  USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
  implicit none

  contains

  !> rank 0 listen
  !> the other ranks send their ACloc data to rank 0
  !> rank 0 checks if no node is assigned to more than one thread
  function checkCoherency(ACin) result(sumError)
    use yowError, only: parallel_abort
    use yowdatapool, only: myrank, comm, nTasks, rkind, rtype
    use yownodepool, only: npa, np_global
    use yowrankModule, only: rank
    use yowMpiModule
    implicit none
    real(rkind), intent(in) :: ACin(npa)
    
    integer(KIND=JWIM) :: fhdl, oRank, IP, IPglobal
#ifdef WAM_HAVE_MPI_F08
    type(mpi_status) :: status
#else
    integer(KIND=JWIM) :: status(MPI_STATUS_SIZE)
#endif
    integer(KIND=JWIM) :: ierr
    integer(KIND=JWIM) :: nTimeSeen(np_global)
    real(rkind) :: newVal, sumError
    real(rkind) :: ACtotal(np_global)
    real(rkind), allocatable :: ACloc(:)
    real(rkind) :: eReal(1)
    fhdl = 740 + myrank
    
    
    ACtotal(:) = 0.0
    sumError = 0.0
    nTimeSeen(:) = 0
    
    if(myrank == 0) then

      do IP=1, npa
        ! note: this threads iplg does not map ghost nodes. have to fix this sometime
        IPglobal = rank(1)%iplg(IP)
        ACtotal(IPglobal) = ACin(IP)
        nTimeSeen(IPglobal) = 1
      end do
      
      do oRank = 2, nTasks
        allocate(ACloc(rank(orank)%npa))
        CALL MPI_RECV(ACloc, rank(orank)%npa ,rtype, oRank-1, 51, comm, status, ierr)
        do IP=1, rank(oRank)%npa
          IPglobal = rank(orank)%iplg(IP)
          newVal = ACloc(IP)
          
          if (nTimeSeen(IPglobal) .eq. 1) then
            sumError = sumError + abs(newVal - ACtotal(IPglobal))
          else
            nTimeSeen(IPglobal) = 1
            ACtotal(IPglobal) = newVal
          end if
        end do

        deallocate(ACloc)
      end do

      eReal(1)=SumError
      do oRank = 2, nTasks
        CALL MPI_SEND(eReal, 1, rtype, oRank-1, 50, comm, ierr)
      end do

    ! other ranks
    else
      CALL MPI_SEND(ACin, npa, rtype, 0, 51, comm, ierr)
      CALL MPI_RECV(eReal, 1, rtype,0, 50, comm, status, ierr)
      SumError=eReal(1)
    endif
  end function
  !
  ! Now computing overall integral
  !
  function ComputeSum(ACin) result(GlobalSum)
    use yowError, only: parallel_abort
    use yowdatapool, only: myrank, comm, nTasks, rkind, rtype
    use yownodepool, only: npa, np_global
    use yowrankModule, only: rank
    use yowMpiModule
    implicit none
    real(rkind), intent(in) :: ACin(npa)
    integer(KIND=JWIM) :: oRank, IP, IPglobal
#ifdef WAM_HAVE_MPI_F08
    type(mpi_status)   :: status
#else
    integer(KIND=JWIM) :: status(MPI_STATUS_SIZE)
#endif
    integer(KIND=JWIM) :: ierr
    integer(KIND=JWIM) :: nTimeSeen(np_global)
    real(rkind) :: newVal, GlobalSum
    real(rkind) :: ACtotal(np_global)
    real(rkind), allocatable :: ACloc(:)
    real(rkind) :: eReal(1)
    ACtotal(:) = 0.0
    GlobalSum = 0.0
    nTimeSeen(:) = 0
    if(myrank == 0) then
      do IP=1, npa
        ! note: this threads iplg does not map ghost nodes. have to fix this sometime
        IPglobal = rank(1)%iplg(IP)
        ACtotal(IPglobal) = ACin(IP)
        nTimeSeen(IPglobal) = 1
      end do
      do oRank = 2, nTasks
        allocate(ACloc(rank(orank)%npa))
        CALL MPI_RECV(ACloc, rank(orank)%npa ,rtype, oRank-1, 51, comm, status, ierr)
        do IP=1, rank(oRank)%npa
          IPglobal = rank(orank)%iplg(IP)
          newVal = ACloc(IP)
          ACtotal(IPglobal) = newVal
        end do
        deallocate(ACloc)
      end do
      GlobalSum=sum(ACtotal)
      eReal(1)=GlobalSum
      do oRank = 2, nTasks
        CALL MPI_SEND(eReal, 1, rtype, oRank-1, 50, comm, ierr)
      end do
    else
      CALL MPI_SEND(ACin, npa, rtype, 0, 51, comm, ierr)
      CALL MPI_RECV(eReal, 1, rtype,0, 50, comm, status, ierr)
      GlobalSum=eReal(1)
    endif

  end function


end module
