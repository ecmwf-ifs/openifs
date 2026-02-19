! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

module dwarf_module
implicit none
private
!------------------------------------------------------------------------------

type, public, abstract :: dwarf
  !! An abstract base type for a command type that holds only a method 
  !!    dwarf%execute( args ) with args a derived type of dwarf_args
contains
  procedure(dwarf_execute), deferred, public :: execute
  procedure(dwarf_final),   deferred, public :: final
end type

!------------------------------------------------------------------------------

type, public :: dwarf_args
  !! An abstract base type for argument to a dwarf%execute()
end type

!------------------------------------------------------------------------------

interface

  subroutine dwarf_execute( this, args )
    import dwarf
    import dwarf_args
    class(dwarf),      intent(inout) :: this
    class(dwarf_args), intent(inout) :: args
  end subroutine

  subroutine dwarf_final( this )
    import dwarf
    class(dwarf), intent(inout) :: this
  end subroutine
  
end interface

!------------------------------------------------------------------------------

end module

