! radiation_cloud_optics_data.f90 - Type to store cloud optical properties
!
! Copyright (C) 2014-2017 ECMWF
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
! License: see the COPYING file for details
!

module radiation_cloud_optics_data

  use parkind1, only : jprb

  implicit none

  !---------------------------------------------------------------------
  ! This type holds the configuration information to compute
  ! cloud optical properties
  type cloud_optics_type
     ! Band-specific coefficients are provided separately in the
     ! shortwave and longwave, and are dimensioned (nband,ncoeff),
     ! where ncoeff depends on the nature of the parameterization
     real(jprb), allocatable, dimension(:,:) :: &
          &  liq_coeff_lw, liq_coeff_sw, &
          &  ice_coeff_lw, ice_coeff_sw
     ! General coefficients are vectors of length ncoeffgen, which
     ! depends on the nature of the parameterization; note that most
     ! parameterizations use only band-specific coefficients
     real(jprb), allocatable, dimension(:) :: &
          &  liq_coeff_gen, ice_coeff_gen

   contains
     procedure :: setup => setup_cloud_optics

  end type cloud_optics_type

contains

  !---------------------------------------------------------------------
  ! Setup cloud optics coefficients by reading them from a file
  subroutine setup_cloud_optics(this, liq_file_name, ice_file_name, &
       &                        iverbose)

    use yomhook,  only           : lhook, dr_hook, jphook
    use easy_netcdf_read_mpi, only : netcdf_file

    class(cloud_optics_type), intent(inout) :: this
    character(len=*), intent(in)            :: liq_file_name, ice_file_name
    integer, intent(in), optional           :: iverbose

    ! The NetCDF file containing the coefficients
    type(netcdf_file)  :: file
    integer            :: iverb
    real(jphook)       :: hook_handle

    if (lhook) call dr_hook('radiation_cloud_optics_data:setup',0,hook_handle)

    if (present(iverbose)) then
      iverb = iverbose
    else
      iverb = 2
    end if

    ! Open the droplet scattering file and configure the way it is
    ! read
    call file%open(trim(liq_file_name), iverbose=iverb)
    call file%transpose_matrices()

    ! Read the band-specific coefficients
    call file%get('coeff_lw',this%liq_coeff_lw)
    call file%get('coeff_sw',this%liq_coeff_sw)

    ! Read the general  coefficients
    if (file%exists('coeff_gen')) then
      call file%get('coeff_gen',this%liq_coeff_gen)
    end if

    ! Close droplet scattering file
    call file%close()

    ! Open the ice scattering file and configure the way it is read
    call file%open(trim(ice_file_name), iverbose=iverb)
    call file%transpose_matrices()

    ! Read the band-specific  coefficients
    call file%get('coeff_lw',this%ice_coeff_lw)
    call file%get('coeff_sw',this%ice_coeff_sw)

    ! Read the general  coefficients
    if (file%exists('coeff_gen')) then
      call file%get('coeff_gen',this%ice_coeff_gen)
    end if

    ! Close ice scattering file
    call file%close()

    if (lhook) call dr_hook('radiation_cloud_optics_data:setup',1,hook_handle)

  end subroutine setup_cloud_optics

end module radiation_cloud_optics_data
