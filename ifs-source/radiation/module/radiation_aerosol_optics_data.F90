! radiation_aerosol_optics_data.F90 - Type to store aerosol optical properties
!
! Copyright (C) 2015-2018 ECMWF
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
! License: see the COPYING file for details
!
! Modifications
!   2017-10-23  R. Hogan  Renamed single-character variables
!   2018-04-20  A. Bozzo  Read optical properties at selected wavelengths


module radiation_aerosol_optics_data

  use parkind1,      only : jprb
  use radiation_io,  only : nulerr, radiation_abort

  implicit none

  private :: get_line

  ! The following provide possible values for
  ! aerosol_optics_type%iclass, which is used to map the user's type
  ! index to the hydrophobic or hydrophilic types loaded from the
  ! aerosol optics file. Initially iclass is equal to
  ! AerosolClassUndefined, which will throw an error if ever the user
  ! tries to use this aerosol type. The user may specify that an
  ! aerosol type is to be ignored in the radiation calculation, in
  ! which case iclass will be set equal to AerosolClassIgnored.
  enum, bind(c) 
     enumerator IAerosolClassUndefined,   IAerosolClassIgnored, &
          &     IAerosolClassHydrophobic, IAerosolClassHydrophilic
  end enum

  integer, parameter :: NMaxStringLength = 2000
  integer, parameter :: NMaxLineLength   = 200

  !---------------------------------------------------------------------
  ! This type holds the configuration information to compute
  ! aerosol optical properties
  type aerosol_optics_type
     ! A vector of length ntype, iclass maps user-defined types on to
     ! the hydrophilic or hydrophobic aerosol classes using the
     ! enumerators above
     integer, allocatable, dimension(:) :: iclass

     ! Also a vector of length ntype, itype maps user-defined types on
     ! to specific hydrophilic or hydrophobic aerosol types
     integer, allocatable, dimension(:) :: itype

     ! Scattering properties are provided separately in the shortwave
     ! and longwave for hydrophobic and hydrophilic aerosols.
     ! Hydrophobic aerosols are dimensioned (nband,n_type_phobic):
     real(jprb), allocatable, dimension(:,:) :: &
          &  mass_ext_sw_phobic, & ! Mass-extinction coefficient (m2 kg-1)
          &  ssa_sw_phobic,      & ! Single scattering albedo
          &  g_sw_phobic,        & ! Asymmetry factor
          &  mass_ext_lw_phobic, & ! Mass-extinction coefficient (m2 kg-1)
          &  ssa_lw_phobic,      & ! Single scattering albedo
          &  g_lw_phobic           ! Asymmetry factor

     ! Hydrophilic aerosols are dimensioned (nband, nrh, n_type_philic):
     real(jprb), allocatable, dimension(:,:,:) :: &
          &  mass_ext_sw_philic, & ! Mass-extinction coefficient (m2 kg-1)
          &  ssa_sw_philic,      & ! Single scattering albedo
          &  g_sw_philic,        & ! Asymmetry factor
          &  mass_ext_lw_philic, & ! Mass-extinction coefficient (m2 kg-1)
          &  ssa_lw_philic,      & ! Single scattering albedo
          &  g_lw_philic           ! Asymmetry factor

     ! Scattering properties at selected wavelengths
     ! (n_mono_wl,n_type_phobic/philic)
     real(jprb), allocatable, dimension(:,:) :: &
          &  mass_ext_mono_phobic, & ! Mass-extinction coefficient (m2 kg-1)
          &  ssa_mono_phobic,      & ! Single scattering albedo
          &  g_mono_phobic,        & ! Asymmetry factor
          &  lidar_ratio_mono_phobic ! Lidar Ratio
     real(jprb), allocatable, dimension(:,:,:) :: &
          &  mass_ext_mono_philic, & ! Mass-extinction coefficient (m2 kg-1)
          &  ssa_mono_philic,      & ! Single scattering albedo
          &  g_mono_philic,        & ! Asymmetry factor
          &  lidar_ratio_mono_philic ! Lidar Ratio

     ! For hydrophilic aerosols, the lower bounds of the relative
     ! humidity bins is a vector of length nrh:
     real(jprb), allocatable, dimension(:) :: &
          &  rh_lower    ! Dimensionless (1.0 = 100% humidity)

     ! Strings describing the aerosol types
     character(len=NMaxStringLength) :: description_phobic_str = ' '
     character(len=NMaxStringLength) :: description_philic_str = ' '

     ! The number of user-defined aerosol types
     integer :: ntype

     ! The number of hydrophobic and hydrophilic types read from the
     ! aerosol optics file
     integer :: n_type_phobic, n_type_philic

     ! Number of relative humidity bins
     integer :: nrh

     ! Number of longwave and shortwave bands of the data in the file,
     ! and monochromatic wavelengths
     integer :: n_bands_lw = 0, n_bands_sw = 0, n_mono_wl = 0

     ! Do we have any hydrophilic types?
     logical :: use_hydrophilic = .true.

     ! Do we have monochromatic optical properties
     logical :: use_monochromatic = .false.

   contains
     procedure :: setup => setup_aerosol_optics
     procedure :: set_hydrophobic_type
     procedure :: set_hydrophilic_type
     procedure :: set_empty_type
     procedure :: set_types
     procedure :: calc_rh_index
     procedure :: print_description

  end type aerosol_optics_type

contains


  !---------------------------------------------------------------------
  ! Setup aerosol optics coefficients by reading them from a file
  subroutine setup_aerosol_optics(this, file_name, ntype, iverbose)

    use yomhook,  only           : lhook, dr_hook, jphook
    use easy_netcdf_read_mpi, only : netcdf_file
    use radiation_io,         only : nulerr, radiation_abort

    class(aerosol_optics_type), intent(inout) :: this
    character(len=*), intent(in)              :: file_name
    integer, intent(in)                       :: ntype
    integer, intent(in), optional             :: iverbose

    ! The NetCDF file containing the aerosol optics data
    type(netcdf_file)  :: file
    integer            :: iverb
    real(jphook)       :: hook_handle

    if (lhook) call dr_hook('radiation_aerosol_optics_data:setup',0,hook_handle)

    if (present(iverbose)) then
       iverb = iverbose
    else
       iverb = 2
    end if

    ! Open the aerosol scattering file and configure the way it is
    ! read
    call file%open(trim(file_name), iverbose=iverb)

    if (file%exists('mass_ext_sw_hydrophilic')) then
      this%use_hydrophilic = .true.
    else
      this%use_hydrophilic = .false.
    end if

    ! Read the raw scattering data
    call file%get('mass_ext_sw_hydrophobic', this%mass_ext_sw_phobic)
    call file%get('ssa_sw_hydrophobic',      this%ssa_sw_phobic)
    call file%get('asymmetry_sw_hydrophobic',this%g_sw_phobic)
    call file%get('mass_ext_lw_hydrophobic', this%mass_ext_lw_phobic)
    call file%get('ssa_lw_hydrophobic',      this%ssa_lw_phobic)
    call file%get('asymmetry_lw_hydrophobic',this%g_lw_phobic)

    call file%get_global_attribute('description_hydrophobic', &
         &                         this%description_phobic_str)

    if (this%use_hydrophilic) then
      call file%get('mass_ext_sw_hydrophilic', this%mass_ext_sw_philic)
      call file%get('ssa_sw_hydrophilic',      this%ssa_sw_philic)
      call file%get('asymmetry_sw_hydrophilic',this%g_sw_philic)
      call file%get('mass_ext_lw_hydrophilic', this%mass_ext_lw_philic)
      call file%get('ssa_lw_hydrophilic',      this%ssa_lw_philic)
      call file%get('asymmetry_lw_hydrophilic',this%g_lw_philic)

      call file%get('relative_humidity1',      this%rh_lower)

      call file%get_global_attribute('description_hydrophilic', &
           &                         this%description_philic_str)
    end if

    ! Read the raw scattering data at selected wavelengths if
    ! available in the input file
    if (file%exists('mass_ext_mono_hydrophobic')) then
      this%use_monochromatic = .true.
      call file%get('mass_ext_mono_hydrophobic', this%mass_ext_mono_phobic)
      call file%get('ssa_mono_hydrophobic',      this%ssa_mono_phobic)
      call file%get('asymmetry_mono_hydrophobic',this%g_mono_phobic)
      call file%get('lidar_ratio_mono_hydrophobic',this%lidar_ratio_mono_phobic)
      if (this%use_hydrophilic) then
        call file%get('mass_ext_mono_hydrophilic', this%mass_ext_mono_philic)
        call file%get('ssa_mono_hydrophilic',      this%ssa_mono_philic)
        call file%get('asymmetry_mono_hydrophilic',this%g_mono_philic)
        call file%get('lidar_ratio_mono_hydrophilic',this%lidar_ratio_mono_philic)
      end if
    else
      this%use_monochromatic = .false.
    end if

    ! Close aerosol scattering file
    call file%close()

    ! Get array sizes
    this%n_bands_lw    = size(this%mass_ext_lw_phobic, 1)
    this%n_bands_sw    = size(this%mass_ext_sw_phobic, 1)
    if (this%use_monochromatic) then
      this%n_mono_wl   = size(this%mass_ext_mono_phobic, 1)
    else
      this%n_mono_wl   = 0
    end if
    this%n_type_phobic = size(this%mass_ext_lw_phobic, 2)

    if (this%use_hydrophilic) then
      this%n_type_philic = size(this%mass_ext_lw_philic, 3)
      this%nrh           = size(this%mass_ext_lw_philic, 2)

      ! Check agreement of dimensions
      if (size(this%mass_ext_lw_philic,1) /= this%n_bands_lw) then
        write(nulerr,'(a,a)') '*** Error: mass extinction for hydrophilic and hydrophobic ', &
             &                'aerosol have different numbers of longwave bands'
        call radiation_abort()
      end if
      if (size(this%mass_ext_sw_philic,1) /= this%n_bands_sw) then
        write(nulerr,'(a,a)') '*** Error: mass extinction for hydrophilic and hydrophobic ', &
             &                'aerosol have different numbers of shortwave bands'
        call radiation_abort()
      end if
      if (size(this%rh_lower) /= this%nrh) then
        write(nulerr,'(a)') '*** Error: size(relative_humidity1) /= size(mass_ext_sw_hydrophilic,2)'
        call radiation_abort()
      end if

    else
      this%n_type_philic = 0
      this%nrh           = 0
    end if

    ! Allocate memory for mapping arrays
    this%ntype = ntype
    allocate(this%iclass(ntype))
    allocate(this%itype(ntype))

    this%iclass = IAerosolClassUndefined
    this%itype  = 0

    if (lhook) call dr_hook('radiation_aerosol_optics_data:setup',1,hook_handle)

  end subroutine setup_aerosol_optics


  !---------------------------------------------------------------------
  ! Map user type "itype" onto stored hydrophobic type "i_type_phobic"
  subroutine set_hydrophobic_type(this, itype, i_type_phobic)

    use yomhook,  only           : lhook, dr_hook, jphook
    class(aerosol_optics_type), intent(inout) :: this
    integer, intent(in)                       :: itype, i_type_phobic
    real(jphook)                              :: hook_handle

    if (lhook) call dr_hook('radiation_aerosol_optics_data:set_hydrophobic_type',0,hook_handle)

    if (itype < 1 .or. itype > this%ntype) then
      write(nulerr,'(a,i0)') '*** Error: aerosol type must be in the range 1 to ', &
           &                  this%ntype
       call radiation_abort('Error setting up aerosols')
    end if
    if (i_type_phobic < 1 .or. i_type_phobic > this%n_type_phobic) then
      write(nulerr,'(a,i0)') '*** Error: hydrophobic type must be in the range 1 to ', &
           &                  this%n_type_phobic
      call radiation_abort('Error setting up aerosols')
    end if

    this%iclass(itype) = IAerosolClassHydrophobic
    this%itype (itype) = i_type_phobic

    if (lhook) call dr_hook('radiation_aerosol_optics_data:set_hydrophobic_type',1,hook_handle)

  end subroutine set_hydrophobic_type


  !---------------------------------------------------------------------
  ! Map user type "itype" onto stored hydrophilic type "i_type_philic"
  subroutine set_hydrophilic_type(this, itype, i_type_philic)

    use yomhook,  only           : lhook, dr_hook, jphook
    class(aerosol_optics_type), intent(inout) :: this
    integer, intent(in)                       :: itype, i_type_philic
    real(jphook)                              :: hook_handle

    if (lhook) call dr_hook('radiation_aerosol_optics_data:set_hydrophilic_type',0,hook_handle)

    if (.not. this%use_hydrophilic) then
      write(nulerr,'(a)') '*** Error: attempt to set hydrophilic aerosol type when no such types present'
      call radiation_abort('Error setting up aerosols')      
    end if

    if (itype < 1 .or. itype > this%ntype) then
      write(nulerr,'(a,i0)') '*** Error: aerosol type must be in the range 1 to ', &
            &          this%ntype
      call radiation_abort('Error setting up aerosols')
    end if
    if (i_type_philic < 1 .or. i_type_philic > this%n_type_philic) then
      write(nulerr,'(a,i0)') '*** Error: hydrophilic type must be in the range 1 to ', &
           &                 this%n_type_philic
      call radiation_abort('Error setting up aerosols')
    end if

    this%iclass(itype) = IAerosolClassHydrophilic
    this%itype (itype) = i_type_philic

    if (lhook) call dr_hook('radiation_aerosol_optics_data:set_hydrophilic_type',1,hook_handle)

  end subroutine set_hydrophilic_type


  !---------------------------------------------------------------------
  ! Set a user type "itype" to be ignored in the radiation scheme
  subroutine set_empty_type(this, itype)

    use yomhook,  only           : lhook, dr_hook, jphook
    class(aerosol_optics_type), intent(inout) :: this
    integer, intent(in)                       :: itype
    real(jphook)                              :: hook_handle

    if (lhook) call dr_hook('radiation_aerosol_optics_data:set_empty_type',0,hook_handle)

    if (itype < 1 .or. itype > this%ntype) then
      write(nulerr,'(a,i0)') '*** Error: aerosol type must be in the range 1 to ', &
           &                 this%ntype
      call radiation_abort('Error setting up aerosols')
    end if

    this%iclass(itype) = IAerosolClassIgnored

    if (lhook) call dr_hook('radiation_aerosol_optics_data:set_empty_type',1,hook_handle)

  end subroutine set_empty_type


  !---------------------------------------------------------------------
  ! Set user types "itypes" to map onto the stored hydrophobic and
  ! hydrophilic types according to its sign and value, with a value of
  ! 0 indicating that this type is to be ignored.  Thus if itypes=(/
  ! 3, 4, -6, 0 /) then user types 1 and 2 map on to hydrophobic types
  ! 3 and 4, user type 3 maps on to hydrophilic type 6 and user type 4
  ! is ignored.
  subroutine set_types(this, itypes)

    use yomhook,  only           : lhook, dr_hook, jphook
    class(aerosol_optics_type), intent(inout) :: this
    integer, dimension(:), intent(in)         :: itypes

    integer :: jtype
    integer :: istart, iend
    real(jphook)                              :: hook_handle

    if (lhook) call dr_hook('radiation_aerosol_optics_data:set_types',0,hook_handle)

    istart = lbound(itypes,1)
    iend   = ubound(itypes,1)

    do jtype = istart, iend
      if (itypes(jtype) == 0) then
        call this%set_empty_type(jtype)
      else if (itypes(jtype) > 0) then
        call this%set_hydrophobic_type(jtype, itypes(jtype))
      else
        call this%set_hydrophilic_type(jtype, -itypes(jtype))
      end if
    end do

    if (lhook) call dr_hook('radiation_aerosol_optics_data:set_types',1,hook_handle)

  end subroutine set_types


  !---------------------------------------------------------------------
  ! Return an index to the relative-humdity array, or zero if no
  ! hydrophilic types are present. This function does so little that
  ! it is best to remove the Dr Hook call.
  function calc_rh_index(this, rh)

    !use yomhook,  only           : lhook, dr_hook, jphook
    class(aerosol_optics_type), intent(inout) :: this
    real(jprb),                 intent(in)    :: rh
    integer                                   :: calc_rh_index
    !real(jphook)                              :: hook_handle

    !if (lhook) call dr_hook('radiation_aerosol_optics_data:calc_rh_index',0,hook_handle)

    if (.not. this%use_hydrophilic) then
      calc_rh_index = 0
    else if (rh > this%rh_lower(this%nrh)) then
      calc_rh_index = this%nrh
    else
      calc_rh_index = 1
      do while (rh > this%rh_lower(calc_rh_index + 1))
        calc_rh_index = calc_rh_index + 1
      end do
    end if

    !if (lhook) call dr_hook('radiation_aerosol_optics_data:calc_rh_index',1,hook_handle)

  end function calc_rh_index


  !---------------------------------------------------------------------
  ! Print a description of the aerosol types to nulout
  subroutine print_description(this, i_type_map)

    use radiation_io, only : nulout

    class(aerosol_optics_type), intent(in) :: this
    integer,                    intent(in) :: i_type_map(:)

    integer :: jtype

    if (size(i_type_map) > 0) then
      write(nulout,'(a)') 'Aerosol mapping:'
    else
      write(nulout,'(a)') 'No aerosol types in radiation scheme'
    end if

    do jtype = 1,size(i_type_map)
      if (i_type_map(jtype) > 0) then
        write(nulout,'(i4,a,a)') jtype, ' -> hydrophobic type ', &
             &  trim(get_line(this%description_phobic_str, i_type_map(jtype)))
      else if (i_type_map(jtype) < 0) then
        write(nulout,'(i4,a,a)') jtype, ' -> hydrophilic type ', &
             &  trim(get_line(this%description_philic_str, -i_type_map(jtype)))
      else
        write(nulout,'(i4,a)') jtype, ' is unused'
      end if
    end do
    
  end subroutine print_description


  !---------------------------------------------------------------------
  ! Private helper function for print_description
  pure function get_line(str,iline) result(line_str)
    character(len=*), intent(in)  :: str
    integer,          intent(in)  :: iline
    character(len=NMaxLineLength) :: line_str
    
    integer :: istart, iend, i_start_new, ioffset, ilength, i_line_current
    logical :: is_fail
    
    i_line_current = 1
    istart = 1
    iend = len(str)
    is_fail = .false.
    line_str = ' '

    ! Find index of first character
    do while (i_line_current < iline)
      i_start_new = scan(str(istart:iend), new_line(' '))
      if (i_start_new == 0) then
        is_fail = .true.
        cycle
      else
        istart = istart + i_start_new
      end if
      i_line_current = i_line_current + 1
    end do
    
    if (.not. is_fail) then
      ! Find index of last character
      ioffset = scan(str(istart:iend), new_line(' '))
      if (ioffset == 0) then
        ilength = len(trim(str(istart:iend)))
      else
        ilength = ioffset - 1
      end if
      
      if (ilength > NMaxLineLength) then
        ilength = NMaxLineLength
      end if
      iend = istart + ilength - 1
      
      line_str = str(istart:iend)
    else
      write(line_str,'(i0,a)') iline, ': <unknown>'
    end if
    
  end function get_line
  
end module radiation_aerosol_optics_data
