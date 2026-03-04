! easy_netcdf.F90 - Module providing convenient NetCDF read/write capability
!
! Copyright (C) 2014-2019 ECMWF
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
! License: see the COPYING file for details
!
! Modifications
!   2017-04-28  R. Hogan   Fix "reshape" when writing 3D array
!   2017-10-23  A. Bozzo   Reading 4-D array
!   2018-03-14  R. Hogan   Fix "reshape" properly this time
!   2019-01-04  R. Hogan   Allow reading and writing a slice of a larger array
!   2019-01-07  R. Hogan   HDF5 writing support, allowing larger files, provided NC_NETCDF4 defined
!   2021-03-15  O. Marsden Add 'first-touch' option to the 4D array get method, to allow improved OpenMP access to arrays.

module easy_netcdf

  use netcdf
  use parkind1,      only : jprb
  use radiation_io,  only : nulout, nulerr, my_abort => radiation_abort

  implicit none

  !---------------------------------------------------------------------
  ! An object of this type provides convenient read or write access to
  ! a NetCDF file
  type netcdf_file
    integer :: ncid             ! NetCDF file ID
    integer :: iverbose         ! Verbosity: 0=errors,   1=warning,
                                !            2=info,     3=progress, 
                                !            4=detailed, 5=debug
    logical :: do_transpose_2d = .false.   ! Transpose 2D arrays on read/write?
    logical :: is_write_mode   = .false.   ! .false. for read, .true. for write
    logical :: is_define_mode  = .true.    ! .true. if in NetCDF define mode
    logical :: is_double_precision = .false. ! Write reals in double precision?
    logical :: do_permute_3d   = .false.   ! Permute 3D arrays on write?
    logical :: do_permute_4d   = .false.   ! Permute 3D arrays on write?
    integer :: i_permute_3d(3) = (/1,2,3/) ! New order of dimensions
    integer :: i_permute_4d(4) = (/1,2,3,4/) ! New order of dimensions
    character(len=511) :: file_name
  contains
    procedure :: open => open_netcdf_file
    procedure :: create => create_netcdf_file
    procedure :: close => close_netcdf_file
    procedure :: get_real_scalar
    procedure :: get_real_vector
    procedure :: get_real_matrix
    procedure :: get_real_array3
    procedure :: get_real_scalar_indexed
    procedure :: get_real_vector_indexed
    procedure :: get_real_matrix_indexed
    procedure :: get_real_array3_indexed
    procedure :: get_real_array4
    generic   :: get => get_real_scalar, get_real_vector, &
         &              get_real_matrix, get_real_array3, &
         &              get_real_array4, &
         &              get_real_scalar_indexed, get_real_vector_indexed, &
         &              get_real_matrix_indexed, get_real_array3_indexed
    procedure :: get_global_attribute

    procedure :: define_dimension
    procedure :: define_variable
    procedure :: put_global_attributes
    procedure :: put_global_attribute
    procedure :: put_real_scalar
    procedure :: put_real_vector
    procedure :: put_real_matrix
    procedure :: put_real_array3
    procedure :: put_real_scalar_indexed
    procedure :: put_real_vector_indexed
    procedure :: put_real_matrix_indexed
    generic   :: put => put_real_scalar, put_real_vector, &
         &              put_real_matrix, put_real_array3, &
         &              put_real_scalar_indexed, put_real_vector_indexed, &
         &              put_real_matrix_indexed
    procedure :: set_verbose
    procedure :: transpose_matrices
    procedure :: double_precision
    procedure :: permute_3d_arrays
    procedure :: get_rank
    procedure :: exists
    procedure :: get_outer_dimension
    procedure :: global_attribute_exists
    procedure, private :: get_array_dimensions
    procedure, private :: get_variable_id
    procedure, private :: end_define_mode
  end type netcdf_file

contains

  ! --- GENERIC SUBROUTINES ---

  !---------------------------------------------------------------------
  ! Open a NetCDF file with name "file_name", optionally specifying the
  ! verbosity level (0-5) and if the file is for writing (the default
  ! is read-only)
  subroutine open_netcdf_file(this, file_name, iverbose, is_write_mode, is_hdf5_file)
    class(netcdf_file)            :: this
    character(len=*), intent(in)  :: file_name
    integer, intent(in), optional :: iverbose
    logical, intent(in), optional :: is_write_mode
    logical, intent(in), optional :: is_hdf5_file ! Only for write mode

    integer                       :: istatus
    integer                       :: i_write_mode


    ! Store verbosity level in object
    if (present(iverbose)) then
      this%iverbose = iverbose
    else
      ! By default announce files being opened and closed, but not
      ! variables read/written
      this%iverbose = 2
    end if

    ! Store read/write mode in object
    if (present(is_write_mode)) then
      this%is_write_mode = is_write_mode
    else
      this%is_write_mode = .false.
    end if

    ! By default we don't transpose 2D arrays on read/write
    this%do_transpose_2d = .false.

    ! Store filename
    this%file_name = file_name

    ! Open file according to write mode
    if (.not. this%is_write_mode) then
      istatus = nf90_open(file_name, NF90_NOWRITE, this%ncid)
      if (this%iverbose >= 2) then
        write(nulout,'(a,a)') 'Reading NetCDF file ', file_name
        !write(nulout,'(a,a,a,i0,a)') 'Reading NetCDF file ', file_name, ' (ID=', this%ncid, ')'
      end if
      this%is_define_mode = .false.
    else
      i_write_mode = NF90_CLOBBER
      ! Check if HDF5 file is to be written (which can be larger)
      if (present(is_hdf5_file)) then
        if (is_hdf5_file) then
#ifdef NC_NETCDF4
          i_write_mode = ior(i_write_mode, NF90_HDF5)
#else
          if (this%iverbose >= 1) then
            write(nulout,'(a,a)') 'Warning: cannot use HDF5 format for writing ', file_name, &
                 &                ' unless compiled with NC_NETCDF4 defined'
          end if
#endif
        end if
      end if

      istatus = nf90_create(file_name, i_write_mode, this%ncid)
      if (this%iverbose >= 2) then
        write(nulout,'(a,a)') 'Writing NetCDF file ', file_name
      end if
      this%is_define_mode = .true.
    end if

    ! Check the file opened correctly
    if (istatus /= NF90_NOERR) then
      write(nulerr,'(a,a,a,a)') '*** Error opening NetCDF file ', &
           &       file_name, ': ', trim(nf90_strerror(istatus))
      call my_abort('Error opening NetCDF file')
    end if

  end subroutine open_netcdf_file


  !---------------------------------------------------------------------
  ! Open a NetCDF file for writing
  subroutine create_netcdf_file(this, file_name, iverbose, is_hdf5_file)
    class(netcdf_file)            :: this
    character(len=*), intent(in)  :: file_name
    integer, intent(in), optional :: iverbose
    logical, intent(in), optional :: is_hdf5_file

    integer                       :: istatus
    integer                       :: i_write_mode

    if (present(iverbose)) then
      this%iverbose = iverbose
    else
      this%iverbose = 2
    end if

    this%do_transpose_2d = .false.

    i_write_mode = NF90_CLOBBER
    ! Check if HDF5 file is to be written (which can be large)
    if (present(is_hdf5_file)) then
      if (is_hdf5_file) then
#ifdef NC_NETCDF4
        i_write_mode = ior(i_write_mode, NF90_HDF5)
#else
        if (this%iverbose >= 1) then
          write(nulout,'(a,a)') 'Warning: cannot use HDF5 format for writing ', file_name, &
               &                ' unless compiled with NC_NETCDF4 defined'
        end if
#endif
      end if
    end if

    istatus = nf90_create(file_name, i_write_mode, this%ncid)
    if (this%iverbose >= 2) then
      write(nulout,'(a,a)') 'Writing NetCDF file ', file_name
      !write(nulout,'(a,a,a,i0,a)') 'Writing NetCDF file ', file_name, ' (ID=', this%ncid, ')'
    end if
    this%is_define_mode = .true.

    if (istatus /= NF90_NOERR) then
      write(nulerr,'(a,a,a)') '*** Error opening NetCDF file ', file_name, &
           &                  ': ', trim(nf90_strerror(istatus))
      stop
    end if
    this%file_name = file_name

  end subroutine create_netcdf_file


  !---------------------------------------------------------------------
  ! Close the NetCDF file
  subroutine close_netcdf_file(this)
    class(netcdf_file) :: this
    integer            :: istatus

    if (this%iverbose >= 4) then
      write(nulout,'(a,a)') 'Closing NetCDF file ', trim(this%file_name)
    end if

    istatus = nf90_close(this%ncid)
    if (istatus /= NF90_NOERR) then
      write(nulerr,'(a,a,a,a)') '*** Error closing NetCDF file ', &
           & trim(this%file_name), ': ', trim(nf90_strerror(istatus))
      stop
    end if

  end subroutine close_netcdf_file


  !---------------------------------------------------------------------
  ! Set the verbosity level from 0 to 5, where the codes have the
  ! following meaning: 0=errors only, 1=warning, 2=info, 3=progress,
  ! 4=detailed, 5=debug
  subroutine set_verbose(this, ival)
    class(netcdf_file) :: this
    integer, optional  :: ival

    if (present(ival)) then
      this%iverbose = ival
    else
      this%iverbose = 2
    end if

  end subroutine set_verbose


  !---------------------------------------------------------------------
  ! Specify whether floating-point arrays should be written in double precision
  subroutine double_precision(this, is_double)
    class(netcdf_file) :: this
    logical, optional  :: is_double

    if (present(is_double)) then
      this%is_double_precision = is_double
    else
      this%is_double_precision = .true.
    end if

  end subroutine double_precision


  !---------------------------------------------------------------------
  ! Specify whether 2D arrays should be transposed on read/write
  subroutine transpose_matrices(this, do_transpose)
    class(netcdf_file) :: this
    logical, optional  :: do_transpose

    if (present(do_transpose)) then
      this%do_transpose_2d = do_transpose
    else
      this%do_transpose_2d = .true.
    end if

  end subroutine transpose_matrices


  !---------------------------------------------------------------------
  ! Specify that 3D arrays should be permuted on write, with the new
  ! dimension order in the input argument "ipermute" (e.g. 3,2,1)
  subroutine permute_3d_arrays(this, ipermute)
    class(netcdf_file)  :: this
    integer, intent(in) :: ipermute(3)

    this%do_permute_3d = .true.
    this%i_permute_3d  = ipermute

  end subroutine permute_3d_arrays


  !---------------------------------------------------------------------
  ! Specify that 4D arrays should be permuted on write, with the new
  ! dimension order in the input argument "ipermute" (e.g. 4,3,2,1)
  subroutine permute_4d_arrays(this, ipermute)
    class(netcdf_file)  :: this
    integer, intent(in) :: ipermute(4)

    this%do_permute_4d = .true.
    this%i_permute_4d  = ipermute

  end subroutine permute_4d_arrays


  ! --- PRIVATE SUBROUTINES ---

  !---------------------------------------------------------------------
  ! Return the NetCDF variable ID for variable "var_name", or abort if
  ! not present
  subroutine get_variable_id(this, var_name, ivarid)
    class(netcdf_file)           :: this
    character(len=*), intent(in) :: var_name
    integer, intent(out)         :: ivarid

    integer                      :: istatus

    istatus = nf90_inq_varid(this%ncid, var_name, ivarid)
    if (istatus /= NF90_NOERR) then
      write(nulerr,'(a,a,a,a)') '*** Error reading NetCDF variable ', &
           & var_name, ': ', trim(nf90_strerror(istatus))
      call my_abort('Error reading NetCDF file')
    end if

  end subroutine get_variable_id


  !---------------------------------------------------------------------
  ! Return the array dimensions of variable with specified ID, along
  ! with the number of dimensions and optionally the total number of
  ! elements, or abort if variable not present
  subroutine get_array_dimensions(this, ivarid, ndims, ndimlens, ntotal)
    class(netcdf_file)             :: this
    integer, intent(in)            :: ivarid
    integer, intent(out)           :: ndims
    integer, intent(out)           :: ndimlens(NF90_MAX_VAR_DIMS)
    integer, intent(out), optional :: ntotal

    integer                        :: j, istatus
    integer                        :: dimids(NF90_MAX_VAR_DIMS)

    istatus = nf90_inquire_variable(this%ncid, ivarid, &
         &                          ndims=ndims, dimids=dimids)
    if (istatus /= NF90_NOERR) then
      write(nulerr,'(a,i0,a,a)') '*** Error inquiring about NetCDF variable with id ', &
           & ivarid, ': ', trim(nf90_strerror(istatus))
      call my_abort('Error reading NetCDF file')
    end if

    ndimlens(:) = 0
    do j = 1,ndims
      istatus = nf90_inquire_dimension(this%ncid, dimids(j), len=ndimlens(j))
      if (istatus /= NF90_NOERR) then
        write(nulerr,'(a,i0,a,i0,a,a)') '*** Error reading length of dimension ', &
             & j, ' of NetCDF variable with id ', ivarid, ': ', &
             & trim(nf90_strerror(istatus))
        call my_abort('Error reading NetCDF file')
      end if
    end do

    if (present(ntotal)) then
      ntotal = 1
      do j = 1, ndims
        ntotal = ntotal * ndimlens(j)
      end do
    end if

  end subroutine get_array_dimensions


  !---------------------------------------------------------------------
  ! End define mode, if in define mode, and check the error code
  ! (errors are possible if variables are too large for the format,
  ! for example)
  subroutine end_define_mode(this)
    class(netcdf_file)             :: this
    integer                        :: istatus
    if (this%is_define_mode) then
      istatus = nf90_enddef(this%ncid)
      if (istatus /= NF90_NOERR) then
        write(nulerr,'(a,a,a,a)') '*** Error ending define mode when writing ', &
             &    trim(this%file_name), ': ', trim(nf90_strerror(istatus))
        call my_abort('Error writing NetCDF file')
      end if
      this%is_define_mode = .false.
    end if
  end subroutine end_define_mode


  ! --- READING SUBROUTINES ---

  !---------------------------------------------------------------------
  ! Return the number of dimensions of variable with name var_name, or
  ! -1 if the variable is not found
  function get_rank(this, var_name) result(ndims)
    class(netcdf_file)           :: this
    character(len=*), intent(in) :: var_name

    integer :: ndims
    integer :: ivarid
    integer :: ndimlens(NF90_MAX_VAR_DIMS)
    integer :: istatus

    istatus = nf90_inq_varid(this%ncid, var_name, ivarid)
    if (istatus /= NF90_NOERR) then
      if (istatus == NF90_ENOTVAR) then
        if (this%iverbose >= 1) then
          write(nulout,'(a,a,a)') '  Warning: variable ', var_name, ' not found'
        end if
      else
        write(nulerr,'(a,a,a)') '*** Error inquiring about variable ', &
             &                  var_name, ': ', trim(nf90_strerror(istatus))
        call my_abort('Error reading NetCDF file')
      end if
      ndims = -1
    else
      call this%get_array_dimensions(ivarid, ndims, ndimlens)
    end if

  end function get_rank


  !---------------------------------------------------------------------
  ! Return the length of the slowest-varying dimension of variable
  ! with name var_name, or -1 if the variable is not found
  function get_outer_dimension(this, var_name) result(n)
    class(netcdf_file)           :: this
    character(len=*), intent(in) :: var_name

    integer :: n
    integer :: ndims
    integer :: ivarid
    integer :: ndimlens(NF90_MAX_VAR_DIMS)
    integer :: istatus

    istatus = nf90_inq_varid(this%ncid, var_name, ivarid)
    if (istatus /= NF90_NOERR) then
      if (istatus == NF90_ENOTVAR) then
        if (this%iverbose >= 1) then
          write(nulout,'(a,a,a)') '  Warning: variable ', var_name, ' not found'
        end if
      else
        write(nulerr,'(a,a,a)') '*** Error inquiring about variable ', &
             &                  var_name, ': ', trim(nf90_strerror(istatus))
        call my_abort('Error reading NetCDF file')
      end if
      n = -1
    else
      call this%get_array_dimensions(ivarid, ndims, ndimlens)
      n = ndimlens(ndims)
    end if

  end function get_outer_dimension


  !---------------------------------------------------------------------
  ! Return true if the variable is present, false otherwise
  function exists(this, var_name) result(is_present)
    class(netcdf_file)           :: this
    character(len=*), intent(in) :: var_name

    logical :: is_present

    integer :: ivarid
    integer :: istatus

    istatus = nf90_inq_varid(this%ncid, var_name, ivarid)
    if (istatus /= NF90_NOERR) then
      is_present = .false.
    else
      is_present = .true.
    end if

  end function exists


  !---------------------------------------------------------------------
  ! Return true if the global attribute is present, false otherwise.
  ! If argument "len" is provided then return false if len is smaller
  ! than the length of the attribute.  This is useful if you have a
  ! fixed array size and want to check whether the attribute will fit
  ! into it.
  function global_attribute_exists(this, attr_name, len) result(is_present)
    class(netcdf_file)            :: this
    character(len=*), intent(in)  :: attr_name
    integer, optional, intent(in) :: len

    logical :: is_present
    integer :: i_attr_len
    integer :: istatus

    istatus = nf90_inquire_attribute(this%ncid, NF90_GLOBAL, attr_name, &
         &                           len=i_attr_len)
    if (istatus /= NF90_NOERR) then
      is_present = .false.
    else 
      is_present = .true.
      if (present(len)) then
        if (i_attr_len > len) then
          is_present = .false.
        end if
      end if
    end if

  end function global_attribute_exists


  !---------------------------------------------------------------------
  ! The method "get" will read either a scalar, vector or matrix
  ! depending on the rank of the output argument. This version reads a
  ! scalar.
  subroutine get_real_scalar(this, var_name, scalar)
    class(netcdf_file)           :: this
    character(len=*), intent(in) :: var_name
    real(jprb), intent(out)      :: scalar

    integer                      :: istatus
    integer                      :: ivarid, ndims
    integer                      :: ndimlens(NF90_MAX_VAR_DIMS)
    integer                      :: j, ntotal

    ! Inquire the ID, shape & size of the variable
    call this%get_variable_id(var_name, ivarid)
    call this%get_array_dimensions(ivarid, ndims, ndimlens)

    ! Compute number of elements of the variable in the file
    ntotal = 1
    do j = 1, ndims
      ntotal = ntotal * ndimlens(j)
    end do

    if (this%iverbose >= 4) then
      write(nulout,'(a,a)') '  Reading ', var_name
    end if

    ! Abort if the number of elements is anything other than 1
    if (ntotal /= 1) then
      write(nulerr,'(a,a,a,i0,a)') '*** Error reading NetCDF variable ', &
           &    var_name, ' with total length ', ntotal, ' as a scalar'
      call my_abort('Error reading NetCDF file')
    end if

    ! Read variable
    istatus = nf90_get_var(this%ncid, ivarid, scalar)
    if (istatus /= NF90_NOERR) then
      write(nulerr,'(a,a,a,a)') '*** Error reading NetCDF variable ', &
           &  var_name, ' as a scalar: ', trim(nf90_strerror(istatus))
      call my_abort('Error reading NetCDF file')
    end if

  end subroutine get_real_scalar


  !---------------------------------------------------------------------
  ! Read a scalar from a larger array, where "index" indexes the most
  ! slowly varying dimension
  subroutine get_real_scalar_indexed(this, var_name, index, scalar)
    class(netcdf_file)           :: this
    character(len=*), intent(in) :: var_name
    integer, intent(in)          :: index
    real(jprb), intent(out)      :: scalar

    integer                      :: istatus
    integer                      :: ivarid, ndims
    integer                      :: ndimlens(NF90_MAX_VAR_DIMS)
    integer                      :: vstart(NF90_MAX_VAR_DIMS)
    integer                      :: j, ntotal

    ! Inquire the ID, shape & size of the variable
    call this%get_variable_id(var_name, ivarid)
    call this%get_array_dimensions(ivarid, ndims, ndimlens)

    ! Compute number of elements of the slice in the file,
    ! i.e. excluding the slowest varying dimension, indexed by "index"
    ntotal = 1
    do j = 1, ndims-1
      ntotal = ntotal * ndimlens(j)
    end do

    if (this%iverbose >= 4) then
      write(nulout,'(a,a,i0,a,a)') '  Reading element ', index, ' of ', var_name
    end if

    ! Abort if the number of elements is anything other than 1
    if (ntotal /= 1) then
      write(nulerr,'(a,a,a,i0,a)') '*** Error reading NetCDF variable ', &
           &    var_name, ' with slice length ', ntotal, ' as a scalar'
      call my_abort('Error reading NetCDF file')
    end if

    if (index < 1 .or. index > ndimlens(ndims)) then
      write(nulerr,'(a,i0,a,a,a,i0)') '*** Error reading element ', index, &
           &  ' of NetCDF variable ', &
           &    var_name, ' with outer dimension ', ndimlens(ndims)
      call my_abort('Error reading NetCDF file')
    end if

    ! Read variable
    vstart(1:ndims-1) = 1
    vstart(ndims)     = index
    istatus = nf90_get_var(this%ncid, ivarid, scalar, start=vstart)
    if (istatus /= NF90_NOERR) then
      write(nulerr,'(a,a,a,a)') '*** Error reading NetCDF variable ', &
           &  var_name, ' as a scalar: ', trim(nf90_strerror(istatus))
      call my_abort('Error reading NetCDF file')
    end if

  end subroutine get_real_scalar_indexed


  !---------------------------------------------------------------------
  ! Read a 1D array into "vector", which must be allocatable and will
  ! be reallocated if necessary
  subroutine get_real_vector(this, var_name, vector)
    class(netcdf_file)           :: this
    character(len=*), intent(in) :: var_name
    real(jprb), allocatable, intent(out) :: vector(:)

    integer                      :: n  ! Length of vector
    integer                      :: istatus
    integer                      :: ivarid, ndims
    integer                      :: ndimlens(NF90_MAX_VAR_DIMS)
    integer                      :: j

    call this%get_variable_id(var_name, ivarid)
    call this%get_array_dimensions(ivarid, ndims, ndimlens)

    ! Ensure variable has only one dimension in the file
    n = 1
    do j = 1, ndims
      n = n * ndimlens(j)
      if (j > 1 .and. ndimlens(j) > 1) then
        write(nulerr,'(a,a,a)') '*** Error reading NetCDF variable ', &
             & var_name, &
             & ' as a vector: all dimensions above the first must be singletons'
        call my_abort('Error reading NetCDF file')
      end if
    end do

    ! Reallocate if necessary
    if (allocated(vector)) then
      if (size(vector) /= n) then
        if (this%iverbose >= 1) then
          write(nulout,'(a,a)') '  Warning: resizing vector to read ', var_name
        end if
        deallocate(vector)
        allocate(vector(n))
      end if
    else
      allocate(vector(n))
    end if

    if (this%iverbose >= 4) then
      write(nulout,'(a,a,a,i0,a)') '  Reading ', var_name, '(', n, ')'
    end if

    ! Read variable
    istatus = nf90_get_var(this%ncid, ivarid, vector)
    if (istatus /= NF90_NOERR) then
      write(nulerr,'(a,a,a,a)') '*** Error reading NetCDF variable ', &
           &  var_name, ' as a vector: ', trim(nf90_strerror(istatus))
      call my_abort('Error reading NetCDF file')
    end if

  end subroutine get_real_vector


  !---------------------------------------------------------------------
  ! Read a vector of data from a larger array; the vector must be
  ! allocatable and will be reallocated if necessary
  subroutine get_real_vector_indexed(this, var_name, index, vector)
    class(netcdf_file)           :: this
    character(len=*), intent(in) :: var_name
    integer, intent(in)          :: index
    real(jprb), allocatable, intent(out) :: vector(:)

    integer                      :: n  ! Length of vector
    integer                      :: istatus
    integer                      :: ivarid, ndims
    integer                      :: ndimlens(NF90_MAX_VAR_DIMS)
    integer                      :: vstart(NF90_MAX_VAR_DIMS)
    integer                      :: vcount(NF90_MAX_VAR_DIMS)
    integer                      :: j

    call this%get_variable_id(var_name, ivarid)
    call this%get_array_dimensions(ivarid, ndims, ndimlens)

    ! Ensure variable has only one dimension aside from the last one
    n = 1
    do j = 1, ndims-1
      n = n * ndimlens(j)
      if (j > 1 .and. ndimlens(j) > 1) then
        write(nulerr,'(a,a,a)') '*** Error reading 1D slice from NetCDF variable ', &
             & var_name, &
             & ': all dimensions except the first and last must be singletons'
        call my_abort('Error reading NetCDF file')
      end if
    end do

    ! Reallocate if necessary
    if (allocated(vector)) then
      if (size(vector) /= n) then
        if (this%iverbose >= 1) then
          write(nulout,'(a,i0,a,a)') '  Warning: resizing vector to length ', n, &
               &                ' to read slice of ', var_name
        end if
        deallocate(vector)
        allocate(vector(n))
      end if
    else
      allocate(vector(n))
    end if

    if (this%iverbose >= 4) then
      write(nulout,'(a,i0,a,a)') '  Reading column ', index, ' of ', var_name
    end if

    if (index < 1 .or. index > ndimlens(ndims)) then
      write(nulerr,'(a,i0,a,a,a,i0)') '*** Error reading element ', index, &
           &  ' of NetCDF variable ', &
           &    var_name, ' with outer dimension ', ndimlens(ndims)
      call my_abort('Error reading NetCDF file')
    end if

    ! Read variable
    vstart(1:ndims-1) = 1
    vstart(ndims)     = index
    vcount(1:ndims-1) = ndimlens(1:ndims-1)
    vcount(ndims)     = 1
    istatus = nf90_get_var(this%ncid, ivarid, vector, &
         &                 start=vstart, count=vcount)
    if (istatus /= NF90_NOERR) then
      write(nulerr,'(a,a,a,a)') '*** Error reading 1D slice of NetCDF variable ', &
           &  var_name, ': ', trim(nf90_strerror(istatus))
      call my_abort('Error reading NetCDF file')
    end if

  end subroutine get_real_vector_indexed


  !---------------------------------------------------------------------
  ! Read 2D array into "matrix", which must be allocatable and will be
  ! reallocated if necessary.  Whether to transpose is specifed by the
  ! final optional argument, but can also be specified by the
  ! do_transpose_2d class data member.
  subroutine get_real_matrix(this, var_name, matrix, do_transp)
    class(netcdf_file)           :: this
    character(len=*), intent(in) :: var_name
    real(jprb), allocatable, intent(out) :: matrix(:,:)
    logical, optional, intent(in):: do_transp ! Transpose data?

    real(jprb), allocatable      :: tmp_matrix(:,:)
    integer                      :: ndimlen1, ndimlen2
    integer                      :: istatus
    integer                      :: ivarid, ndims
    integer                      :: ndimlens(NF90_MAX_VAR_DIMS)
    integer                      :: j, ntotal
    logical                      :: do_transpose

    ! Decide whether to transpose the array
    if (present(do_transp)) then
      do_transpose = do_transp
    else
      do_transpose = this%do_transpose_2d
    end if

    call this%get_variable_id(var_name, ivarid)
    call this%get_array_dimensions(ivarid, ndims, ndimlens)

    ! Ensure the variable has no more than two non-singleton
    ! dimensions
    ntotal = 1
    do j = 1, ndims
      ntotal = ntotal * ndimlens(j)
      if (j > 2 .and. ndimlens(j) > 1) then
        write(nulerr,'(a,a,a)') '*** Error reading NetCDF variable ', &
           & var_name, &
           & ' as a matrix: all dimensions above the second must be singletons'
        call my_abort('Error reading NetCDF file')
      end if
    end do

    ! Work out dimension lengths
    if (ndims >= 2) then
      ndimlen1 = ndimlens(1)
      ndimlen2 = ntotal/ndimlen1
    else
      ndimlen1 = ntotal
      ndimlen2 = 1
    end if

    if (do_transpose) then
      ! Read and transpose
      allocate(tmp_matrix(ndimlen1, ndimlen2))

      ! Reallocate if necessary
      if (allocated(matrix)) then
        if (size(matrix,1) /= ndimlen2 .or. size(matrix,2) /= ndimlen1) then
          if (this%iverbose >= 1) then
            write(nulout,'(a,a)') '  Warning: resizing matrix to read ', var_name
          end if
          deallocate(matrix)
          allocate(matrix(ndimlen2, ndimlen1))
        end if
      else
        allocate(matrix(ndimlen2, ndimlen1))
      end if

      if (this%iverbose >= 4) then
        write(nulout,'(a,a,a,i0,a,i0,a)') '  Reading ', var_name, '(', &
             &                            ndimlen2, ',', ndimlen1, ')'
      end if

      istatus = nf90_get_var(this%ncid, ivarid, tmp_matrix)
      matrix = transpose(tmp_matrix)
      deallocate(tmp_matrix)
    else
      ! Read data without transposition

      ! Reallocate if necessary
      if (allocated(matrix)) then
        if (size(matrix,1) /= ndimlen1 .or. size(matrix,2) /= ndimlen2) then
          if (this%iverbose >= 1) then
            write(nulout,'(a,a)') '  Warning: resizing matrix to read ', var_name
          end if
          allocate(matrix(ndimlen1, ndimlen2))
        end if
      else
        allocate(matrix(ndimlen1, ndimlen2))
      end if

      if (this%iverbose >= 4) then
        write(nulout,'(a,a,a,i0,a,i0,a)') '  Reading ', var_name, '(', &
             &                            ndimlen1, ',', ndimlen2, ')'
      end if

      istatus = nf90_get_var(this%ncid, ivarid, matrix)
    end if

    if (istatus /= NF90_NOERR) then
      write(nulerr,'(a,a,a,a)') '*** Error reading NetCDF variable ', &
           &    var_name, ' as a matrix: ', trim(nf90_strerror(istatus))
      call my_abort('Error reading NetCDF file')
    end if

  end subroutine get_real_matrix


  !---------------------------------------------------------------------
  ! Read matrix of data from a larger array, which must be allocatable
  ! and will be reallocated if necessary.  Whether to transpose is
  ! specifed by the final optional argument, but can also be specified
  ! by the do_transpose_2d class data member.
  subroutine get_real_matrix_indexed(this, var_name, index, matrix, do_transp)
    class(netcdf_file)           :: this
    character(len=*), intent(in) :: var_name
    integer, intent(in)          :: index
    real(jprb), allocatable, intent(out) :: matrix(:,:)
    logical, optional, intent(in):: do_transp ! Transpose data?

    real(jprb), allocatable      :: tmp_matrix(:,:)
    integer                      :: ndimlen1, ndimlen2
    integer                      :: istatus
    integer                      :: ivarid, ndims
    integer                      :: ndimlens(NF90_MAX_VAR_DIMS)
    integer                      :: vstart(NF90_MAX_VAR_DIMS)
    integer                      :: vcount(NF90_MAX_VAR_DIMS)
    integer                      :: j, ntotal
    logical                      :: do_transpose

    ! Decide whether to transpose the array
    if (present(do_transp)) then
      do_transpose = do_transp
    else
      do_transpose = this%do_transpose_2d
    end if

    call this%get_variable_id(var_name, ivarid)
    call this%get_array_dimensions(ivarid, ndims, ndimlens)

    ! Ensure the variable has no more than two non-singleton
    ! dimensions aside from the last one
    ntotal = 1
    do j = 1, ndims-1
      ntotal = ntotal * ndimlens(j)
      if (j > 2 .and. ndimlens(j) > 1) then
        write(nulerr,'(a,a,a)') '*** Error reading 2D slice from NetCDF variable ', &
           & var_name, &
           & ': all dimensions except the first, second and last must be singletons'
        call my_abort('Error reading NetCDF file')
      end if
    end do

    if (index < 1 .or. index > ndimlens(ndims)) then
      write(nulerr,'(a,i0,a,a,a,i0)') '*** Error reading element ', index, &
           &  ' of NetCDF variable ', &
           &    var_name, ' with outer dimension ', ndimlens(ndims)
      call my_abort('Error reading NetCDF file')
    end if

    ! Work out dimension lengths
    if (ndims >= 2) then
      ndimlen1 = ndimlens(1)
      ndimlen2 = ntotal/ndimlen1
    else
      ndimlen1 = ntotal
      ndimlen2 = 1
    end if

    vstart(1:ndims-1) = 1
    vstart(ndims)     = index
    vcount(1:ndims-1) = ndimlens(1:ndims-1)
    vcount(ndims)     = 1

    if (do_transpose) then
      ! Read and transpose
      allocate(tmp_matrix(ndimlen1, ndimlen2))

      ! Reallocate if necessary
      if (allocated(matrix)) then
        if (size(matrix,1) /= ndimlen2 .or. size(matrix,2) /= ndimlen1) then
          if (this%iverbose >= 1) then
            write(nulout,'(a,a)') '  Warning: resizing matrix to read ', var_name
          end if
          allocate(matrix(ndimlen2, ndimlen1))
        end if
      else
        allocate(matrix(ndimlen2, ndimlen1))
      end if

      if (this%iverbose >= 4) then
        write(nulout,'(a,i0,a,a,a,i0,a,i0,a)') '  Reading slice ', index, &
             &  ' of ', var_name, ' as ', ndimlen2, 'x', ndimlen1, ' array'
      end if

      istatus = nf90_get_var(this%ncid, ivarid, tmp_matrix, &
           &                 start=vstart, count=vcount)
      matrix = transpose(tmp_matrix)
      deallocate(tmp_matrix)
    else
      ! Read data without transposition

      ! Reallocate if necessary
      if (allocated(matrix)) then
        if (size(matrix,1) /= ndimlen1 .or. size(matrix,2) /= ndimlen2) then
          if (this%iverbose >= 1) then
            write(nulout,'(a,a)') '  Warning: resizing matrix to read ', var_name
          end if
          allocate(matrix(ndimlen1, ndimlen2))
        end if
      else
        allocate(matrix(ndimlen1, ndimlen2))
      end if

      if (this%iverbose >= 4) then
        write(nulout,'(a,i0,a,a,a,i0,a,i0,a)') '  Reading slice ', index, &
             &  ' of ', var_name, ' as ', ndimlen1, 'x', ndimlen2, ' array'
      end if

      istatus = nf90_get_var(this%ncid, ivarid, matrix, &
           &                 start=vstart, count=vcount)
    end if

    if (istatus /= NF90_NOERR) then
      write(nulerr,'(a,a,a,a)') '*** Error reading 2D slice of NetCDF variable ', &
           &    var_name, ': ', trim(nf90_strerror(istatus))
      call my_abort('Error reading NetCDF file')
    end if

  end subroutine get_real_matrix_indexed


  !---------------------------------------------------------------------
  ! Read 3D array into "var", which must be allocatable and will be
  ! reallocated if necessary.  Whether to pemute is specifed by the
  ! final optional argument
  subroutine get_real_array3(this, var_name, var, ipermute)
    class(netcdf_file)                   :: this
    character(len=*), intent(in)         :: var_name
    real(jprb), allocatable, intent(out) :: var(:,:,:)
    integer, optional, intent(in)        :: ipermute(3)

    real(jprb), allocatable   :: var_permute(:,:,:)
    integer                   :: ndimlen1, ndimlen2, ndimlen3
    integer                   :: istatus
    integer                   :: ivarid, ndims
    integer                   :: ndimlens(NF90_MAX_VAR_DIMS)
    integer                   :: j, ntotal
    integer                   :: n_dimlens_permuted(3)
    integer                   :: i_permute_3d(3)
    logical                   :: do_permute

    ! Decide whether to permute
    if (present(ipermute)) then
      do_permute = .true.
      i_permute_3d = ipermute
    else
      do_permute = this%do_permute_3d
      i_permute_3d = this%i_permute_3d
    end if

    call this%get_variable_id(var_name, ivarid)
    call this%get_array_dimensions(ivarid, ndims, ndimlens)

    ! Ensure the variable has no more than three non-singleton
    ! dimensions
    ntotal = 1
    do j = 1, ndims
      ntotal = ntotal * ndimlens(j)
      if (j > 3 .and. ndimlens(j) > 1) then
        write(nulerr,'(a,a,a)') '*** Error reading NetCDF variable ', &
           & var_name, &
           & ' as a 3D array: all dimensions above the third must be singletons'
        call my_abort('Error reading NetCDF file')
      end if
    end do

    ! Work out dimension lengths
    if (ndims >= 3) then
      ndimlen1 = ndimlens(1)
      ndimlen2 = ndimlens(2)
      ndimlen3 = ntotal/(ndimlen1*ndimlen2)
    else if (ndims >= 2) then
      ndimlen1 = ndimlens(1)
      ndimlen2 = ntotal/ndimlen1
      ndimlen3 = 1
    else
      ndimlen1 = ntotal
      ndimlen2 = 1
      ndimlen3 = 1
    end if

    ! Deallocate if necessary
    if (allocated(var)) then
      deallocate(var)
    end if

    if (do_permute) then
      ! Read and permute
      allocate(var_permute(ndimlen1, ndimlen2, ndimlen3))
      n_dimlens_permuted(i_permute_3d) = ndimlens(1:3)

      ! Reallocate if necessary
      if (allocated(var)) then
        if (size(var,1) /= n_dimlens_permuted(1) &
             &  .or. size(var,2) /= n_dimlens_permuted(2) &
             &  .or. size(var,3) /= n_dimlens_permuted(3)) then
          if (this%iverbose >= 1) then
            write(nulout,'(a,a)') '  Warning: resizing array to read ', var_name
          end if
          deallocate(var)
          allocate(var(n_dimlens_permuted(1), n_dimlens_permuted(2), &
               &       n_dimlens_permuted(3)))
        end if
      else
        allocate(var(n_dimlens_permuted(1), n_dimlens_permuted(2), &
             &       n_dimlens_permuted(3)))
      end if

      if (this%iverbose >= 4) then
        write(nulout,'(a,a,a,i0,i0,i0,a)') '  Reading ', var_name, &
             & ' (permuted dimensions ', i_permute_3d, ')'
      end if

      istatus = nf90_get_var(this%ncid, ivarid, var_permute)
      var = reshape(var_permute, n_dimlens_permuted, order=i_permute_3d)
      deallocate(var_permute)

    else
      ! Read data without permutation

      ! Reallocate if necessary
      if (allocated(var)) then
        if (size(var,1) /= ndimlen1 &
             &  .or. size(var,2) /= ndimlen2 &
             &  .or. size(var,3) /= ndimlen3) then
          if (this%iverbose >= 1) then
            write(nulout,'(a,a)') '  Warning: resizing array to read ', var_name
          end if
          deallocate(var)
          allocate(var(ndimlen1, ndimlen2, ndimlen3))
        end if
      else
        allocate(var(ndimlen1, ndimlen2, ndimlen3))
      end if

      if (this%iverbose >= 4) then
        write(nulout,'(a,a,a,i0,a,i0,a,i0,a)') '  Reading ', var_name, &
             &            '(', ndimlen1, ',', ndimlen2, ',', ndimlen3, ')'
      end if

      istatus = nf90_get_var(this%ncid, ivarid, var)
    end if

    if (istatus /= NF90_NOERR) then
      write(nulerr,'(a,a,a,a)') '*** Error reading NetCDF variable ', &
           &  var_name, ' as a 3D array: ', trim(nf90_strerror(istatus))
      call my_abort('Error reading NetCDF file')
    end if

  end subroutine get_real_array3


  !---------------------------------------------------------------------
  ! Read 3D array of data from a larger-dimensional array, which must
  ! be allocatable and will be reallocated if necessary.  Whether to
  ! pemute is specifed by the final optional argument
  subroutine get_real_array3_indexed(this, var_name, index, var, ipermute)
    class(netcdf_file)                   :: this
    character(len=*), intent(in)         :: var_name
    integer, intent(in)                  :: index
    real(jprb), allocatable, intent(out) :: var(:,:,:)
    integer, optional, intent(in)        :: ipermute(3)

    real(jprb), allocatable   :: var_permute(:,:,:)
    integer                   :: ndimlen1, ndimlen2, ndimlen3
    integer                   :: istatus
    integer                   :: ivarid, ndims
    integer                   :: ndimlens(NF90_MAX_VAR_DIMS)
    integer                   :: vstart(NF90_MAX_VAR_DIMS)
    integer                   :: vcount(NF90_MAX_VAR_DIMS)
    integer                   :: j, ntotal
    integer                   :: n_dimlens_permuted(3)
    integer                   :: i_permute_3d(3)
    logical                   :: do_permute

    ! Decide whether to permute
    if (present(ipermute)) then
      do_permute = .true.
      i_permute_3d = ipermute
    else
      do_permute = this%do_permute_3d
      i_permute_3d = this%i_permute_3d
    end if

    call this%get_variable_id(var_name, ivarid)
    call this%get_array_dimensions(ivarid, ndims, ndimlens)

    ! Ensure the variable has no more than three non-singleton
    ! dimensions aside from the last one
    ntotal = 1
    do j = 1, ndims-1
      ntotal = ntotal * ndimlens(j)
      if (j > 3 .and. ndimlens(j) > 1) then
        write(nulerr,'(a,a,a)') '*** Error reading 3D slice from NetCDF variable ', &
           & var_name, &
           & ': all dimensions above the third must be singletons'
        call my_abort('Error reading NetCDF file')
      end if
    end do

    if (index < 1 .or. index > ndimlens(ndims)) then
      write(nulerr,'(a,i0,a,a,a,i0)') '*** Error reading element ', index, &
           &  ' of NetCDF variable ', &
           &    var_name, ' with outer dimension ', ndimlens(ndims)
      call my_abort('Error reading NetCDF file')
    end if

    ! Work out dimension lengths
    if (ndims >= 3) then
      ndimlen1 = ndimlens(1)
      ndimlen2 = ndimlens(2)
      ndimlen3 = ntotal/(ndimlen1*ndimlen2)
    else if (ndims >= 2) then
      ndimlen1 = ndimlens(1)
      ndimlen2 = ntotal/ndimlen1
      ndimlen3 = 1
    else
      ndimlen1 = ntotal
      ndimlen2 = 1
      ndimlen3 = 1
    end if

    vstart(1:ndims-1) = 1
    vstart(ndims)     = index
    vcount(1:ndims-1) = ndimlens(1:ndims-1)
    vcount(ndims)     = 1

    if (do_permute) then
      ! Read and permute
      allocate(var_permute(ndimlen1, ndimlen2, ndimlen3))
      n_dimlens_permuted(i_permute_3d) = ndimlens(1:3)

      ! Reallocate if necessary
      if (allocated(var)) then
        if (size(var,1) /= n_dimlens_permuted(1) &
             &  .or. size(var,2) /= n_dimlens_permuted(2) &
             &  .or. size(var,3) /= n_dimlens_permuted(3)) then
          if (this%iverbose >= 1) then
            write(nulout,'(a,a)') '  Warning: resizing array to read ', var_name
          end if
          deallocate(var)
          allocate(var(n_dimlens_permuted(1), n_dimlens_permuted(2), &
               &       n_dimlens_permuted(3)))
        end if
      else
        allocate(var(n_dimlens_permuted(1), n_dimlens_permuted(2), &
             &       n_dimlens_permuted(3)))
      end if

      if (this%iverbose >= 4) then
        write(nulout,'(a,i0,a,a,a,i0,i0,i0,a)') '  Reading slice ', index, &
             &  ' of ', var_name, &
             & ' (permuted dimensions ', i_permute_3d, ')'
      end if

      istatus = nf90_get_var(this%ncid, ivarid, var_permute, &
           &                 start=vstart, count=vcount)
      var = reshape(var_permute, n_dimlens_permuted, order=i_permute_3d)
      deallocate(var_permute)

    else
      ! Read data without permutation

      ! Reallocate if necessary
      if (allocated(var)) then
        if (size(var,1) /= ndimlen1 &
             &  .or. size(var,2) /= ndimlen2 &
             &  .or. size(var,3) /= ndimlen3) then
          if (this%iverbose >= 1) then
            write(nulout,'(a,a)') '  Warning: resizing array to read ', var_name
          end if
          deallocate(var)
          allocate(var(ndimlen1, ndimlen2, ndimlen3))
        end if
      else
        allocate(var(ndimlen1, ndimlen2, ndimlen3))
      end if

      if (this%iverbose >= 4) then
        write(nulout,'(a,i0,a,a,a,i0,a,i0,a,i0,a)') '  Reading slice ', index, &
             &  ' of ', var_name, ' as ', ndimlen1, 'x', ndimlen2, 'x', &
             &  ndimlen3, 'array'
      end if

      istatus = nf90_get_var(this%ncid, ivarid, var, &
           &                 start=vstart, count=vcount)
    end if

    if (istatus /= NF90_NOERR) then
      write(nulerr,'(a,a,a,a)') '*** Error reading 3D slice of NetCDF variable ', &
           &  var_name, ': ', trim(nf90_strerror(istatus))
      call my_abort('Error reading NetCDF file')
    end if

  end subroutine get_real_array3_indexed


  !---------------------------------------------------------------------
  ! Read 4D array into "var", which must be allocatable and will be
  ! reallocated if necessary.  Whether to pemute is specifed by the
  ! ipermute optional argument. For the non-permuted case, OpenMP 
  ! thread-optimized location of array is enabled by setting optional 
  ! argument ld_first_touch to true. This results in zero-ing of "var"
  ! inside an OpenMP loop before reading the array in from NetCDF.
 
  subroutine get_real_array4(this, var_name, var, ipermute, ld_first_touch)
    class(netcdf_file)                   :: this
    character(len=*), intent(in)         :: var_name
    real(jprb), allocatable, intent(out) :: var(:,:,:,:)
    integer, optional, intent(in)        :: ipermute(4)
    logical, optional, intent(in)        :: ld_first_touch

    real(jprb), allocatable   :: var_permute(:,:,:,:)
    integer                   :: ndimlen1, ndimlen2, ndimlen3, ndimlen4
    integer                   :: istatus
    integer                   :: ivarid, ndims
    integer                   :: ndimlens(NF90_MAX_VAR_DIMS)
    integer                   :: j, ntotal
    integer                   :: n_dimlens_permuted(4)
    integer                   :: i_permute_4d(4)
    logical                   :: do_permute

    logical                   :: ll_first_touch
    integer                   :: ii,jj,kk,mm

    ! Decide whether to permute
    if (present(ipermute)) then
      do_permute = .true.
      i_permute_4d = ipermute
    else
      do_permute = this%do_permute_4d
      i_permute_4d = this%i_permute_4d
    end if

    if (present(ld_first_touch)) then
      ll_first_touch = ld_first_touch
    else
      ll_first_touch = .false.
    endif

    call this%get_variable_id(var_name, ivarid)
    call this%get_array_dimensions(ivarid, ndims, ndimlens)

    ! Ensure the variable has no more than three non-singleton
    ! dimensions
    ntotal = 1
    do j = 1, ndims
      ntotal = ntotal * ndimlens(j)
      if (j > 4 .and. ndimlens(j) > 1) then
        write(nulerr,'(a,a,a)') '*** Error reading NetCDF variable ', &
           & var_name, &
           & ' as a 4D array: all dimensions above the third must be singletons'
        call my_abort('Error reading NetCDF file')
      end if
    end do

    ! Work out dimension lengths
    if (ndims >= 4) then
      ndimlen1 = ndimlens(1)
      ndimlen2 = ndimlens(2)
      ndimlen3 = ndimlens(3)
      ndimlen4 = ntotal/(ndimlen1*ndimlen2*ndimlen3)
    else if (ndims >= 3) then
      ndimlen1 = ndimlens(1)
      ndimlen2 = ndimlens(2)
      ndimlen3 = ntotal/(ndimlen1*ndimlen2)
    else if (ndims >= 2) then
      ndimlen1 = ndimlens(1)
      ndimlen2 = ntotal/ndimlen1
      ndimlen3 = 1
    else
      ndimlen1 = ntotal
      ndimlen2 = 1
      ndimlen3 = 1
    end if

    ! Deallocate if necessary
    if (allocated(var)) then
      deallocate(var)
    end if

    if (do_permute) then
      ! Read and permute - not tested
      allocate(var_permute(ndimlen1, ndimlen2, ndimlen3, ndimlen4))
      n_dimlens_permuted(i_permute_4d) = ndimlens(1:4)

      ! Reallocate if necessary
      if (allocated(var)) then
        if (size(var,1) /= n_dimlens_permuted(1) &
             &  .or. size(var,2) /= n_dimlens_permuted(2) &
             &  .or. size(var,3) /= n_dimlens_permuted(3) &
             &  .or. size(var,4) /= n_dimlens_permuted(4)) then
          if (this%iverbose >= 1) then
            write(nulout,'(a,a)') '  Warning: resizing array to read ', var_name
          end if
          deallocate(var)
          allocate(var(n_dimlens_permuted(1), n_dimlens_permuted(2), &
               &       n_dimlens_permuted(3), n_dimlens_permuted(4)))
        end if
      else
        allocate(var(n_dimlens_permuted(1), n_dimlens_permuted(2), &
             &       n_dimlens_permuted(3), n_dimlens_permuted(4)))
      end if

      if (this%iverbose >= 4) then
        write(nulout,'(a,a,a,i0,i0,i0,a)') '  Reading ', var_name, &
             & ' (permuted dimensions ', i_permute_4d, ')'
      end if

      istatus = nf90_get_var(this%ncid, ivarid, var_permute)
      var = reshape(var_permute, n_dimlens_permuted, order=i_permute_4d)
      deallocate(var_permute)

    else
      ! Read data without permutation

      ! Reallocate if necessary
      if (allocated(var)) then
        if (size(var,1) /= ndimlen1 &
             &  .or. size(var,2) /= ndimlen2 &
             &  .or. size(var,3) /= ndimlen3 &
             &  .or. size(var,4) /= ndimlen4) then
          if (this%iverbose >= 1) then
            write(nulout,'(a,a)') '  Warning: resizing array to read ', var_name
          end if
          deallocate(var)
          allocate(var(ndimlen1, ndimlen2, ndimlen3, ndimlen4))
        end if
      else
        allocate(var(ndimlen1, ndimlen2, ndimlen3, ndimlen4))
      end if

      if (this%iverbose >= 4) then
        write(nulout,'(a,a,a,i0,a,i0,a,i0,a,i0,a)') '  Reading ', var_name, &
             &            '(', ndimlen1, ',', ndimlen2, ',', ndimlen3,',', ndimlen4, ')'
      end if

      if (ll_first_touch) then
        !!first touch
        !$OMP parallel PRIVATE(ii,jj,kk,mm)
        do mm=1,ndimlen4
          !$OMP DO PRIVATE(ii,jj,kk) collapse(2) schedule(static) 
          do kk=1,ndimlen3
            do jj=1,ndimlen2
              !$omp simd
              do ii=1,ndimlen1
                var(ii,jj,kk,mm) = 0.0_jprb        
              enddo
            enddo
          enddo
          !$OMP END DO
        end do
        !$OMP END PARALLEL
      end if

      istatus = nf90_get_var(this%ncid, ivarid, var)
    end if

    if (istatus /= NF90_NOERR) then
      write(nulerr,'(a,a,a,a)') '*** Error reading NetCDF variable ', &
           &  var_name, ' as a 4D array: ', trim(nf90_strerror(istatus))
      call my_abort('Error reading NetCDF file')
    end if

  end subroutine get_real_array4


  !---------------------------------------------------------------------
  ! Get a global attribute as a character string
  subroutine get_global_attribute(this, attr_name, attr_str)
    class(netcdf_file) :: this

    character(len=*), intent(in)    :: attr_name
    character(len=*), intent(inout) :: attr_str

    integer :: i_attr_len
    integer :: istatus
    integer :: j

    istatus = nf90_inquire_attribute(this%ncid, NF90_GLOBAL, attr_name, &
         &                           len = i_attr_len)
    if (istatus /= NF90_NOERR) then
      write(nulerr,'(a,a,a,a)') '*** Error reading size of global attribute ', attr_name, &
           &                    ': ', trim(nf90_strerror(istatus))
      call my_abort('Error reading NetCDF file')
    end if

    ! Allocatable character strings not supported one enough compilers
    ! yet
    !    if (allocated(attr_str)) then
    !      deallocate(attr_str)
    !    end if
    !    allocate(character(len=i_attr_len) :: attr_str)
    if (len(attr_str) < i_attr_len) then
      write(nulerr,'(a,a)') '*** Not enough space to read global attribute ', attr_name
      call my_abort('Error reading NetCDF file')
    end if

    istatus = nf90_get_att(this%ncid, NF90_GLOBAL, attr_name, attr_str)
    if (istatus /= NF90_NOERR) then
      write(nulerr,'(a,a,a,a)') '*** Error reading global attribute ', attr_name, &
           &                    ': ', trim(nf90_strerror(istatus))
      call my_abort('Error reading NetCDF file')
    end if

    ! Pad with blanks since nf90_get_att does not do this
    do j = i_attr_len+1,len(attr_str)
      attr_str(j:j) = ' '
    end do

  end subroutine get_global_attribute


  ! --- WRITING SUBROUTINES ---

  !---------------------------------------------------------------------
  ! Define a dimension with name dim_name of length n (or 0 to
  ! indicate the unlimited dimension)
  subroutine define_dimension(this, dim_name, n)
    class(netcdf_file)           :: this
    character(len=*), intent(in) :: dim_name
    integer, intent(in)          :: n
    integer                      :: idimid, istatus

    istatus = nf90_def_dim(this%ncid, dim_name, n, idimid)
    if (istatus /= NF90_NOERR) then
      write(nulerr,'(a,a,a,a)') '*** Error defining ', dim_name, &
           &   ' as a dimension: ', trim(nf90_strerror(istatus))
      call my_abort('Error writing NetCDF file')
    end if

    if (this%iverbose >= 4) then
      write(nulout,'(a,a,a,i0)') '  Defining dimension ',trim(dim_name), &
           & ' of length ', n
    end if

  end subroutine define_dimension


  !---------------------------------------------------------------------
  ! Define a variable with name var_name, and specify its shape via
  ! the dim1_name, dim2_name and dim3_name optional arguments, which
  ! are strings referring to existing dimensions. If none are
  ! specified, the variable will be a scalar. The optional arguments
  ! long_name, units and comment write string attributes with these
  ! names.
  subroutine define_variable(this, var_name, dim1_name, dim2_name, dim3_name, &
       &                     long_name, units_str, comment_str, is_double)
    class(netcdf_file)                     :: this
    character(len=*), intent(in)           :: var_name
    character(len=*), intent(in), optional :: long_name, units_str, comment_str
    character(len=*), intent(in), optional :: dim1_name, dim2_name, dim3_name
    logical,          intent(in), optional :: is_double

    integer :: istatus, ndims, ivarid
    integer, dimension(NF90_MAX_VAR_DIMS) :: idimids
    logical :: is_double_precision

    if (present(dim1_name)) then
      ! Variable is at least one dimensional
      ndims = 1
      istatus = nf90_inq_dimid(this%ncid, dim1_name, idimids(1))
      if (istatus /= NF90_NOERR) then
        write(nulerr,'(a,a,a,a)') '*** Error inquiring ID of dimension ', &
             &             dim1_name, ': ', trim(nf90_strerror(istatus))
        call my_abort('Error writing NetCDF file')
      end if
      if (present(dim2_name)) then
        ! Variable is at least two dimensional
        ndims = 2
        istatus = nf90_inq_dimid(this%ncid, dim2_name, idimids(2))
        if (istatus /= NF90_NOERR) then
          write(nulerr,'(a,a,a)') '*** Error inquiring ID of dimension ', &
               &           dim2_name, ': ', trim(nf90_strerror(istatus))
          call my_abort('Error writing NetCDF file')
        end if
        if (present(dim3_name)) then
          ! Variable is at least three dimensional
          ndims = 3
          istatus = nf90_inq_dimid(this%ncid, dim3_name, idimids(3))
          if (istatus /= NF90_NOERR) then
            write(nulerr,'(a,a,a,a)') '*** Error inquiring ID of dimension ', &
                 &             dim3_name, ': ', trim(nf90_strerror(istatus))
            call my_abort('Error writing NetCDF file')
          end if
        end if
      end if
    else
      ! Variable is a scalar
      ndims = 0
    end if

    ! Read output precision from optional argument "is_double" if
    ! present, otherwise from default output precision for this file
    if (present(is_double)) then
      is_double_precision = is_double
    else
      is_double_precision = this%is_double_precision
    end if

    ! Define variable
    if (is_double_precision) then
      istatus = nf90_def_var(this%ncid, var_name, NF90_DOUBLE, idimids(1:ndims), &
           & ivarid)
    else
      istatus = nf90_def_var(this%ncid, var_name, NF90_FLOAT, idimids(1:ndims), &
           & ivarid)
    end if
    if (istatus /= NF90_NOERR) then
      write(nulerr,'(a,a,a,a)') '*** Error defining variable ', var_name, &
           &                    ': ', trim(nf90_strerror(istatus))
      call my_abort('Error writing NetCDF file')
    end if

    ! Add optional attributes
    if (present(long_name)) then
      istatus = nf90_put_att(this%ncid, ivarid, "long_name", long_name)
      if (this%iverbose >= 4) then
        write(nulout,'(a,a,a,a)') '  Defining ',trim(var_name),': ',long_name
      end if
    else
      if (this%iverbose >= 4) then
        write(nulout,'(a,a)') '  Defining ',trim(var_name)
      end if
    end if
    if (present(units_str)) then
      istatus = nf90_put_att(this%ncid, ivarid, "units", units_str)
    end if
    if (present(comment_str)) then
      istatus = nf90_put_att(this%ncid, ivarid, "comment", comment_str)
    end if

  end subroutine define_variable


  !---------------------------------------------------------------------
  ! Put CF-compliant global attributes into the file
  subroutine put_global_attributes(this, title_str, inst_str, source_str, &
       &  comment_str, ref_str)
    class(netcdf_file)                     :: this

    character(len=*), intent(in), optional :: title_str
    character(len=*), intent(in), optional :: inst_str
    character(len=*), intent(in), optional :: source_str
    character(len=*), intent(in), optional :: comment_str
    character(len=*), intent(in), optional :: ref_str

    character(len=32)   :: date_time_str
    character(len=1024) :: command_line_str
    character(len=1056) :: history_str

    integer :: time_vals(8)
    integer :: i ! status

    call date_and_time(values=time_vals)
    call get_command(command_line_str)

    write(date_time_str,"(i0.4,'-',i0.2,'-',i0.2,' ',i0.2,':',i0.2,':',i0.2)") &
         &   time_vals(1), time_vals(2), time_vals(3), time_vals(5), time_vals(6), time_vals(7)

    history_str = trim(date_time_str) // ': ' // trim(command_line_str)


    if (present(title_str))   i=nf90_put_att(this%ncid, NF90_GLOBAL, "title", title_str)
    if (present(inst_str))    i=nf90_put_att(this%ncid, NF90_GLOBAL, "institution", inst_str)
    if (present(source_str))  i=nf90_put_att(this%ncid, NF90_GLOBAL, "source", source_str)

    i = nf90_put_att(this%ncid, NF90_GLOBAL, "history", history_str)

    if (present(comment_str)) i=nf90_put_att(this%ncid, NF90_GLOBAL, "comment", comment_str)
    if (present(ref_str))     i=nf90_put_att(this%ncid, NF90_GLOBAL, "references", ref_str)

  end subroutine put_global_attributes


  !---------------------------------------------------------------------
  ! Put a non-standard global attribute into the file
  subroutine put_global_attribute(this, attr_name, attr_str)
    class(netcdf_file) :: this

    character(len=*), intent(in) :: attr_name, attr_str

    integer :: istatus

    istatus = nf90_put_att(this%ncid, NF90_GLOBAL, trim(attr_name), trim(attr_str))

    if (istatus /= NF90_NOERR) then
      write(nulerr,'(a,a,a,a)') '*** Error writing global attribute ', attr_name, &
           &                    ': ', trim(nf90_strerror(istatus))
      call my_abort('Error writing NetCDF file')
    end if

  end subroutine put_global_attribute


  !---------------------------------------------------------------------
  ! The "put" method saves a scalar, vector or matrix into the
  ! variable with name var_name, according to the rank of the var
  ! argument. This version saves a scalar.
  subroutine put_real_scalar(this, var_name, var)
    class(netcdf_file)             :: this
    character(len=*), intent(in)   :: var_name
    real(jprb), intent(in)         :: var

    integer :: ivarid, ndims, ntotal, istatus
    integer :: ndimlens(NF90_MAX_VAR_DIMS)

    ! If we are in define mode, exit define mode
    call this%end_define_mode()

    ! Check the variable is a scalar
    call this%get_variable_id(var_name, ivarid)
    call this%get_array_dimensions(ivarid, ndims, ndimlens, ntotal)
    if (ntotal /= 1) then
      write(nulerr,'(a,a,a,i0)') '*** Error: attempt to write scalar to ', &
           &                 var_name, ' which has total length ', ntotal
      call my_abort('Error writing NetCDF file')
    end if

    ! Save the scalar
    istatus = nf90_put_var(this%ncid, ivarid, var)
    if (istatus /= NF90_NOERR) then
      write(nulerr,'(a,a,a,a)') '*** Error writing scalar ', var_name, ': ', &
           &                    trim(nf90_strerror(istatus))
      call my_abort('Error writing NetCDF file')
    end if

  end subroutine put_real_scalar


  !---------------------------------------------------------------------
  ! Save a scalar.
  subroutine put_real_scalar_indexed(this, var_name, index, var)
    class(netcdf_file)             :: this
    character(len=*), intent(in)   :: var_name
    real(jprb), intent(in)         :: var
    integer, intent(in)            :: index

    integer :: ivarid, ndims, ntotal, istatus
    integer :: ndimlens(NF90_MAX_VAR_DIMS)
    integer :: vstart(NF90_MAX_VAR_DIMS)

    ! If we are in define mode, exit define mode
    call this%end_define_mode()

    ! Check the variable is a scalar
    call this%get_variable_id(var_name, ivarid)
    call this%get_array_dimensions(ivarid, ndims, ndimlens)
    if (index < 1 .or. index > ndimlens(ndims)) then
      write(nulerr,'(a,i0,a,a,a,i0)') '*** Error: attempt to write scalar to element ', &
           &  index, ' of ', var_name, ' which has outer dimension  ', ndimlens(ndims)
      call my_abort('Error writing NetCDF file')
    end if

    ! Save the scalar
    vstart(1:ndims-1) = 1
    vstart(ndims)     = index
    istatus = nf90_put_var(this%ncid, ivarid, var, start=vstart)
    if (istatus /= NF90_NOERR) then
      write(nulerr,'(a,a,a,a)') '*** Error writing scalar to ', var_name, ': ', &
           &                    trim(nf90_strerror(istatus))
      call my_abort('Error writing NetCDF file')
    end if

  end subroutine put_real_scalar_indexed


  !---------------------------------------------------------------------
  ! Save a vector with name var_name in the file
  subroutine put_real_vector(this, var_name, var)
    class(netcdf_file)             :: this
    character(len=*), intent(in)   :: var_name
    real(jprb), intent(in)         :: var(:)

    integer :: ivarid, ndims, ntotal, istatus
    integer :: ndimlens(NF90_MAX_VAR_DIMS)

    call this%end_define_mode()

    ! Check the vector is of the right length
    call this%get_variable_id(var_name, ivarid)
    call this%get_array_dimensions(ivarid, ndims, ndimlens, ntotal)
    if (ntotal /= size(var)) then
      write(nulerr,'(a,i0,a,a,a,i0)') '*** Error: attempt to write vector of length ', &
           & size(var), ' to ', var_name, ' which has total length ', ntotal
      call my_abort('Error writing NetCDF file')
    end if

    ! Save the vector
    istatus = nf90_put_var(this%ncid, ivarid, var)
    if (istatus /= NF90_NOERR) then
      write(nulerr,'(a,a,a,a)') '*** Error writing vector ', var_name, ': ', &
           &                    trim(nf90_strerror(istatus))
      call my_abort('Error writing NetCDF file')
    end if

  end subroutine put_real_vector


  !---------------------------------------------------------------------
  ! Save a vector slice with name var_name in the file
  subroutine put_real_vector_indexed(this, var_name, index, var)
    class(netcdf_file)             :: this
    character(len=*), intent(in)   :: var_name
    real(jprb), intent(in)         :: var(:)
    integer, intent(in)            :: index

    integer :: ivarid, ndims, ntotal, istatus
    integer :: ndimlens(NF90_MAX_VAR_DIMS)
    integer :: vstart(NF90_MAX_VAR_DIMS)
    integer :: vcount(NF90_MAX_VAR_DIMS)

    call this%end_define_mode()

    ! Check the vector is of the right length
    call this%get_variable_id(var_name, ivarid)
    call this%get_array_dimensions(ivarid, ndims, ndimlens, ntotal)
    ntotal = ntotal / ndimlens(ndims)
    if (ntotal /= size(var)) then
      write(nulerr,'(a,i0,a,a,a,i0)') '*** Error: attempt to write vector of length ', &
           & size(var), ' to slice of ', var_name, ' which has length ', ntotal
      call my_abort('Error writing NetCDF file')
    end if
    if (index < 1 .or. index > ndimlens(ndims)) then
      write(nulerr,'(a,i0,a,a,a,i0)') '*** Error: attempt to write vector to slice ', &
           &  index, ' of ', var_name, ' which has outer dimension  ', ndimlens(ndims)
      call my_abort('Error writing NetCDF file')
    end if

    ! Save the vector
    vstart(1:ndims-1) = 1
    vstart(ndims)     = index
    vcount(1:ndims-1) = ndimlens(1:ndims-1)
    vcount(ndims)     = 1
    istatus = nf90_put_var(this%ncid, ivarid, var, start=vstart, count=vcount)
    if (istatus /= NF90_NOERR) then
      write(nulerr,'(a,a,a,a)') '*** Error writing vector to ', var_name, ': ', &
           &                    trim(nf90_strerror(istatus))
      call my_abort('Error writing NetCDF file')
    end if

  end subroutine put_real_vector_indexed


  !---------------------------------------------------------------------
  ! Save a matrix with name var_name in the file, transposing its
  ! dimensions if either optional argument transp is .true., or the
  ! transpose_matrices method has already been called.
  subroutine put_real_matrix(this, var_name, var, do_transp)
    class(netcdf_file)             :: this
    character(len=*), intent(in)   :: var_name
    real(jprb), intent(in)         :: var(:,:)
    real(jprb), allocatable        :: var_transpose(:,:)
    logical, optional, intent(in):: do_transp

    integer :: ivarid, ndims, ntotal, nvarlen, istatus
    integer :: ndimlens(NF90_MAX_VAR_DIMS)

    logical :: do_transpose

    if (present(do_transp)) then
      do_transpose = do_transp
    else
      do_transpose = this%do_transpose_2d
    end if

    call this%end_define_mode()

    call this%get_variable_id(var_name, ivarid)
    call this%get_array_dimensions(ivarid, ndims, ndimlens, ntotal)

    nvarlen = size(var,1)*size(var,2)

    ! Check the total size of the variable to be stored (but receiving
    ! ntotal is zero then there must be an unlimited dimension)
    if (ntotal /= size(var) .and. ntotal /= 0) then
      write(nulerr,'(a,i0,a,a,a,i0)') '*** Error: attempt to write matrix of total size ', &
           & nvarlen, ' to ', var_name, ' which has total size ', ntotal
      call my_abort('Error writing NetCDF file')
    end if

    if (do_transpose) then
      ! Save the matrix with transposition
      if (this%iverbose >= 4) then
        write(nulout,'(a,a,a)') '  Writing ', var_name, &
             & ' (transposing dimensions)'
      end if
      allocate(var_transpose(size(var,2), size(var,1)))
      var_transpose = transpose(var)
      istatus = nf90_put_var(this%ncid, ivarid, var_transpose)
      deallocate(var_transpose)
    else
      ! Save the matrix without transposition
      if (this%iverbose >= 4) then
        write(nulout,'(a,a)') '  Writing ', var_name
      end if
      istatus = nf90_put_var(this%ncid, ivarid, var)
    end if

    if (istatus /= NF90_NOERR) then
      write(nulerr,'(a,a,a,a)') '*** Error writing matrix ', var_name, &
           &                    ': ', trim(nf90_strerror(istatus))
      call my_abort('Error writing NetCDF file')
    end if

  end subroutine put_real_matrix


  !---------------------------------------------------------------------
  ! Save a matrix slice with name var_name in the file, transposing its
  ! dimensions if either optional argument transp is .true., or the
  ! transpose_matrices method has already been called.
  subroutine put_real_matrix_indexed(this, var_name, index, var, do_transp)
    class(netcdf_file)             :: this
    character(len=*), intent(in)   :: var_name
    real(jprb), intent(in)         :: var(:,:)
    integer, intent(in)            :: index

    real(jprb), allocatable        :: var_transpose(:,:)
    logical, optional, intent(in):: do_transp

    integer :: ivarid, ndims, ntotal, nvarlen, istatus
    integer :: ndimlens(NF90_MAX_VAR_DIMS)
    integer :: vstart(NF90_MAX_VAR_DIMS)
    integer :: vcount(NF90_MAX_VAR_DIMS)

    logical :: do_transpose

    if (present(do_transp)) then
      do_transpose = do_transp
    else
      do_transpose = this%do_transpose_2d
    end if

    call this%end_define_mode()

    call this%get_variable_id(var_name, ivarid)
    call this%get_array_dimensions(ivarid, ndims, ndimlens, ntotal)

    nvarlen = size(var,1)*size(var,2)

    ! Check the total size of the variable to be stored (but receiving
    ! ntotal is zero then there must be an unlimited dimension)
    ntotal = ntotal / ndimlens(ndims)
    if (ntotal /= size(var) .and. ntotal /= 0) then
      write(nulerr,'(a,i0,a,a,a,i0)') '*** Error: attempt to write matrix of total size ', &
           & nvarlen, ' to ', var_name, ' which has total size ', ntotal
      call my_abort('Error writing NetCDF file')
    end if

    vstart(1:ndims-1) = 1
    vstart(ndims)     = index
    vcount(1:ndims-1) = ndimlens(1:ndims-1)
    vcount(ndims)     = 1

    if (do_transpose) then
      ! Save the matrix with transposition
      if (this%iverbose >= 4) then
        write(nulout,'(a,a,a)') '  Writing ', var_name, &
             & ' (transposing dimensions)'
      end if
      allocate(var_transpose(size(var,2), size(var,1)))
      var_transpose = transpose(var)
      istatus = nf90_put_var(this%ncid, ivarid, var_transpose, start=vstart, count=vcount)
      deallocate(var_transpose)
    else
      ! Save the matrix without transposition
      if (this%iverbose >= 4) then
        write(nulout,'(a,a)') '  Writing ', var_name
      end if
      istatus = nf90_put_var(this%ncid, ivarid, var, start=vstart, count=vcount)
    end if

    if (istatus /= NF90_NOERR) then
      write(nulerr,'(a,a,a,a)') '*** Error writing matrix ', var_name, &
           &                    ': ', trim(nf90_strerror(istatus))
      call my_abort('Error writing NetCDF file')
    end if

  end subroutine put_real_matrix_indexed


  !---------------------------------------------------------------------
  ! Save a 3D array with name var_name in the file.  The optional
  ! argument permute specifies that the dimensions should first be
  ! permuted according to the three integers therein (or if
  ! permute_3d_arrays has already been called). ipermute is
  ! interpretted such that if OLD and NEW are 3-element vectors
  ! containing the size of each dimension in memory and in the written
  ! file, respectively, then NEW=OLD(ipermute).
  subroutine put_real_array3(this, var_name, var, ipermute)
    class(netcdf_file)             :: this
    character(len=*), intent(in)   :: var_name
    real(jprb), intent(in)         :: var(:,:,:)
    real(jprb), allocatable        :: var_permute(:,:,:)
    integer, optional, intent(in)  :: ipermute(3)

    integer :: ivarid, ndims, ntotal, nvarlen, istatus
    integer :: ndimlens(NF90_MAX_VAR_DIMS)

    logical :: do_permute          ! Do we permute?
    integer :: i_permute_3d(3)
    integer :: n_dimlens_permuted(3)
    integer :: i_order(3)

    ! Decide whether to permute
    if (present(ipermute)) then
      do_permute   = .true.
      i_permute_3d = ipermute
    else
      do_permute   = this%do_permute_3d
      i_permute_3d = this%i_permute_3d
    end if

    call this%end_define_mode()

    ! Check total size
    call this%get_variable_id(var_name, ivarid)
    call this%get_array_dimensions(ivarid, ndims, ndimlens, ntotal)
    nvarlen = size(var,1)*size(var,2)*size(var,3)
    if (ntotal /= size(var)) then
      write(nulerr,'(a,i0,a,a,a,i0)') '*** Error: attempt to write array of total size ', &
           & nvarlen, ' to ', var_name, ' which has total size ', ntotal
      call my_abort('Error writing NetCDF file')
    end if

    if (do_permute) then
      ! Save array after permuting dimensions
      if (this%iverbose >= 4) then
        write(nulout,'(a,a,a,i0,i0,i0,a)') '  Writing ', var_name, &
             & ' (permuted dimensions: ', i_permute_3d, ')'
      end if
      n_dimlens_permuted = (/ size(var,i_permute_3d(1)), &
           &                  size(var,i_permute_3d(2)), &
           &                  size(var,i_permute_3d(3))  /)
      if (this%iverbose >= 5) then
        write(nulout,'(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a)') '    (', &
             &  n_dimlens_permuted(1), ',', n_dimlens_permuted(2), &
             &  ',', n_dimlens_permuted(3), ') -> (', ndimlens(1), &
             &  ',', ndimlens(2), ',', ndimlens(3), ')'
      end if
      allocate(var_permute(n_dimlens_permuted(1), &
           &   n_dimlens_permuted(2), n_dimlens_permuted(3)))
      ! Due to the odd way that ORDER works in Fortran RESHAPE, we
      ! need to do this:
      i_order(i_permute_3d) = (/ 1, 2, 3 /)
      var_permute = reshape(var, n_dimlens_permuted, order=i_order)
      istatus = nf90_put_var(this%ncid, ivarid, var_permute)
      deallocate(var_permute)
    else
      ! Save array without permuting dimensions
      if (this%iverbose >= 4) then
        write(nulout,'(a,a)') '  Writing ', var_name
      end if
      istatus = nf90_put_var(this%ncid, ivarid, var)
    end if

    if (istatus /= NF90_NOERR) then
      write(nulerr,'(a,a,a,a)') '*** Error writing array ', var_name, &
           &                    ': ', trim(nf90_strerror(istatus))
      call my_abort('Error writing NetCDF file')
    end if

  end subroutine put_real_array3

end module easy_netcdf
