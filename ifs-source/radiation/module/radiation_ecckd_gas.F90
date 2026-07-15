! radiation_ecckd_gas.F90 - type representing a single ecCKD gas
!
! (C) Copyright 2020- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
! License: see the COPYING file for details
!

module radiation_ecckd_gas

  use parkind1, only : jprb
  use radiation_gas_constants

  implicit none

  public

  ! Concentration dependence of individual gases
  enum, bind(c)
    enumerator :: IConcDependenceNone = 0, &
         &        IConcDependenceLinear, &
         &        IConcDependenceLUT, &
         &        IConcDependenceRelativeLinear
  end enum

  !---------------------------------------------------------------------
  ! This derived type describes a correlated k-distribution
  ! representation of an individual gas (including composite gases)
  type ckd_gas_type

    ! Code identifying the gas, from the codes in the
    ! radiation_gas_constants module
    integer :: i_gas_code = -1

    ! One of the IConcDependence* enumerators
    integer :: i_conc_dependence

    ! Molar absorption coefficient in m2 mol-1. If
    ! i_conc_dependence==IConcDependenceNone then it is the absorption
    ! cross section per mole of dry air.  If
    ! conc_dependence==IConcDependenceLinear|IConcDependenceRelativeLinear,
    ! it is the absorption cross section per mole of the gas in
    ! question. It is dimensioned (g_point,pressure,temperature).
    real(jprb), allocatable :: molar_abs(:,:,:)
    
    ! If i_conc_dependence==IConcDependenceLUT then we have an
    ! additional dimension for concentration. It is dimensioned
    ! (g_point,pressure,temperature,conc)
    real(jprb), allocatable :: molar_abs_conc(:,:,:,:)

    ! If i_conc_dependence==IConcDependenceRelativeLinear then the
    ! following reference concentration is subtracted from the actual
    ! concentration before the result is multiplied by the mass
    ! absorption coefficient
    real(jprb) :: reference_mole_frac = 0.0_jprb

    ! Mole fraction coordinate variable if
    ! i_conc_dependence==IConcDependenceLUT
    real(jprb) :: log_mole_frac1 = 0.0_jprb, d_log_mole_frac = 1.0_jprb
    integer    :: n_mole_frac = 0

  contains

    procedure :: read => read_ckd_gas
!    procedure :: deallocate => deallocate_ckd_gas

  end type ckd_gas_type

contains

  !---------------------------------------------------------------------
  ! Read information about the representation of a single gas from a
  ! NetCDF file, identifying it with code i_gas_code
  subroutine read_ckd_gas(this, file, gas_name, i_gas_code)

    use easy_netcdf_read_mpi, only : netcdf_file

    class(ckd_gas_type), intent(inout) :: this
    type(netcdf_file),   intent(inout) :: file
    character(len=*),    intent(in)    :: gas_name
    integer,             intent(in)    :: i_gas_code
    
    ! Local storage for mole fraction coordinate variable
    real(jprb), allocatable :: mole_fraction(:)

    this%i_gas_code = i_gas_code

    call file%get(gas_name // "_conc_dependence_code", this%i_conc_dependence)
    if (this%i_conc_dependence == IConcDependenceLut) then
      call file%get(gas_name // "_molar_absorption_coeff", &
           &        this%molar_abs_conc)
      call file%get(gas_name // "_mole_fraction", mole_fraction)
      this%log_mole_frac1  = log(mole_fraction(1))
      this%n_mole_frac     = size(mole_fraction)
      this%d_log_mole_frac = (log(mole_fraction(size(mole_fraction))) &
           &                  - this%log_mole_frac1) / (this%n_mole_frac-1)
      deallocate(mole_fraction)
    else
      call file%get(gas_name // "_molar_absorption_coeff", &
           &        this%molar_abs)
    end if

    if (this%i_conc_dependence == IConcDependenceRelativeLinear) then
      call file%get(gas_name // "_reference_mole_fraction", &
           &        this%reference_mole_frac)
    end if

  end subroutine read_ckd_gas

end module radiation_ecckd_gas
