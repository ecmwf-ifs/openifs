! (C) Copyright 2001- Aron Roland (Roland & Partner, Germany).
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

#include "yowincludes.h"

!> Has some subroutine to make a nice error message
module yowError
  USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
  implicit none

  integer(KIND=JWIM), public :: stat = 0
  
  contains

  !> \note fro elfe_sgp.F90
  subroutine parallel_abort(string, error)
    use yowDatapool, only: comm
    use yowMpiModule
    implicit none

    character(*),optional,intent(in) :: string !string to print
    integer(KIND=JWIM),optional,intent(in) :: error       !mpi errorcode
    integer(KIND=JWIM) :: ierr,i
    logical :: lopen
    integer(KIND=JWIM) :: sl
    ! MPI_MAX_ERROR_STRING = 1024
    character(1024) :: errorstring
    integer(KIND=JWIM) :: myrank

    ! Get rank
    call mpi_comm_rank(comm, myrank,ierr)
    if(ierr/=MPI_SUCCESS) write(*,*) "parallel_abort: ierr=", ierr

    inquire(11,opened=lopen)

    if(present(string)) then
      write(*,'(i4,2a)') myrank,': ABORT: ', trim(string)
      if(lopen) write(11,'(i4,2a)') myrank,': ABORT: ',string
    endif

    if(present(error)) then
      if(error /= MPI_SUCCESS) then
        ! get errorstring associatet fuer errorcode err
        call mpi_error_string(error, errorstring, sl, ierr)
        if(ierr/=MPI_SUCCESS) write(*,*) "parallel_abort: ierr=", ierr
        write(*,'(i4,2a)') myrank,': MPI ERROR: ', errorstring(1:sl)
        if(lopen) write(11,'(i4,2a)') myrank,': MPI ERROR: ', errorstring
      endif
      do i=1,200; inquire(i,opened=lopen); if(lopen) close(i); enddo;
      call mpi_abort(comm, error,ierr)
      if(ierr/=MPI_SUCCESS) write(*,*) "parallel_abort: ierr=", ierr
    else
      do i=1,200; inquire(i,opened=lopen); if(lopen) close(i); enddo;
      call mpi_abort(comm, 0,ierr)
      if(ierr/=MPI_SUCCESS) write(*,*) "parallel_abort: ierr=", ierr
    endif

#if GCC_VERSION >= 40800
    call backtrace()
#endif
  end subroutine parallel_abort
  
  
  !> print various error strings and exit.
  !> Call this to print an error string and optional line number, file and MPI error string
  !> \param[in] string Errorstring
  !> \param[in] line Line number
  !> \param[in] file Filename
  !> \param[in] errno The MPI error number which is translated into an error string
  subroutine abort(string, line, file, errno)
    use yowDatapool, only: comm     
    use yowMpiModule
    implicit none
    ! Errorstring to print
    character(*), optional, intent(in) :: string
    ! Linenumber to print
    integer(KIND=JWIM),      optional, intent(in) :: line
    ! Filename to print
    character(*), optional, intent(in) :: file
    ! MPI error number to translate
    integer(KIND=JWIM),      optional, intent(in) :: errno
    ! Linenumber as string
    character(50) :: lineNumber
    ! MPI_MAX_ERROR_STRING = 1024
    ! MPI Errorstring
    character(MPI_MAX_ERROR_STRING) :: errorstring
    ! The rank of this thread
    integer(KIND=JWIM) :: myrank
    ! real MPI errorsting lengt
    integer(KIND=JWIM) :: stringLengh
    !
    integer(KIND=JWIM) :: ierr

    ! Get rank
    call mpi_comm_rank(comm, myrank,ierr)
!     if(ierr/=MPI_SUCCESS) write(*,*) "parallel_abort: ierr=", ierr

    ! Always print rank
    write(*, '(i2,a)', advance='no') myrank, " "

    ! Print a simple "ERROR" when no MPI error number was given because the MPI error string contain an "ERROR" allready
    if(.not. present(errno)) then
      write(*,'(a)', advance='no' ) " ERROR "
    endif
    
    ! print file and linenumber
    if(present(file)) then
      write(*,'(a)',advance='no' ) file
      
      if(present(line)) then
        Write(lineNumber, '(i10)') line
        write(*, '(2a)', advance='no') ":", trim(adjustl(lineNumber))
      endif
      
      write(*, '(a)', advance='no') " "
    endif
    
    ! if only linenumber is present, add an "Line:" string
    if(.not. present(file) .and. present(line)) then
        Write(lineNumber, '(i10)') line
        write(*, '(2a)', advance='no') "Line:", trim(adjustl(lineNumber))
        write(*, '(a)', advance='no') " "
    endif

    if(stat /= 0) then
      write(*,'(a,i4)', advance='no') "error id", stat
    endif
    
    ! print the errror string
    if(present(string)) then
      write(*,'(a)', advance='no') trim(string)
    endif
    
    ! translate and print the MPI error string
    if(present(errno)) then
      if(errno /= MPI_SUCCESS) then
        call mpi_error_string(errno, errorstring, stringLengh, ierr)
        write(*,'(2a)', advance='no') 'MPI ERROR: ', errorstring(1:stringLengh)
      endif
    endif


    
    write(*,*)
#if GCC_VERSION >= 40800
    call backtrace()
#endif
    stop
  end subroutine

  !> print warning
  !> Call this to print an warning string and optional line number, and file
  !> \param[in] string warnstring
  !> \param[in] line Line number
  !> \param[in] file Filename
  subroutine warn(string, line, file)
    use yowDatapool, only: comm
    use yowMpiModule
    implicit none
    ! Errorstring to print
    character(*), optional, intent(in) :: string
    ! Linenumber to print
    integer(KIND=JWIM),      optional, intent(in) :: line
    ! Filename to print
    character(*), optional, intent(in) :: file
    ! Linenumber as string
    character(50) :: lineNumber
    ! The rank of this thread
    integer(KIND=JWIM) :: myrank
    !
    integer(KIND=JWIM) :: ierr

    ! Get rank
    call mpi_comm_rank(comm, myrank,ierr)
!     if(ierr/=MPI_SUCCESS) write(*,*) "parallel_abort: ierr=", ierr

    ! Always print rank
    write(*, '(i2,a)', advance='no') myrank, " "

    write(*,'(a)', advance='no' ) " WARN "

    ! print file and linenumber
    if(present(file)) then
      write(*,'(a)',advance='no' ) file

      if(present(line)) then
        Write(lineNumber, '(i10)') line
        write(*, '(2a)', advance='no') ":", trim(adjustl(lineNumber))
      endif

      write(*, '(a)', advance='no') " "
    endif

    ! if only linenumber is present, add an "Line:" string
    if(.not. present(file) .and. present(line)) then
        Write(lineNumber, '(i10)') line
        write(*, '(2a)', advance='no') "Line:", trim(adjustl(lineNumber))
        write(*, '(a)', advance='no') " "
    endif

    ! print the errror string
    if(present(string)) then
      write(*,'(a)', advance='no') string
    endif

    write(*,*)
  end subroutine
end module yowError
