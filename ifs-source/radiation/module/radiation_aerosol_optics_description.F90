! radiation_aerosol_optics_description.F90 - Type to store aerosol optics metadata
!
! (C) Copyright 2022- ECMWF.
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
!

module radiation_aerosol_optics_description

  use parkind1,      only : jprb

  implicit none
  public

  !---------------------------------------------------------------------
  ! This type holds the metadata from an aerosol optical property
  ! file, enabling the user to request the index to the aerosol type
  ! with particular properties.  Note that string information is held
  ! in the form of 2D arrays of single characters, so comparison to
  ! character strings requires the to_string helper function at the
  ! end of this file.
  type aerosol_optics_description_type

    ! Two-character code describing the aerosol family, dimensioned
    ! (2,naer), e.g.
    !   SS: Sea salt
    !   OM: Organic matter
    !   SU: Sulfate
    !   OB: Secondary organic biogenic
    !   OA: Secondary organic anthropogenic
    !   AM: Fine-mode ammonium sulfate
    !   NI: Nitrate
    !   DD: Desert dust
    !   BC: Black carbon
    character(len=1), allocatable :: code_phobic(:,:)
    character(len=1), allocatable :: code_philic(:,:)

    ! Size bin, typically 1-2 or 1-3 in order fine to coarse, or zero
    ! if no division by size is used, dimensioned (naer)
    integer, allocatable :: bin_phobic(:)
    integer, allocatable :: bin_philic(:)

    ! Character string characterizing the optical model, e.g. OPAC,
    ! GACP, GLOMAP, Dubovik2002 etc.
    character(len=1), allocatable :: optical_model_phobic(:,:)
    character(len=1), allocatable :: optical_model_philic(:,:)

    ! The user can call preferred_optical_model to specify that a
    ! certain optical model for a certain aerosol family is to be
    ! preferred when get_index is called
    logical, allocatable :: is_preferred_phobic(:)
    logical, allocatable :: is_preferred_philic(:)

    ! Verbosity level
    integer :: iverbose
    
  contains
    procedure :: read
    procedure :: preferred_optical_model
    procedure :: get_index

  end type aerosol_optics_description_type

contains

  !---------------------------------------------------------------------
  ! Read optical property file file_name into an
  ! aerosol_optics_description_type object
  subroutine read(this, file_name, iverbose)

    use yomhook,              only : lhook, dr_hook, jphook
    use easy_netcdf,          only : netcdf_file

    class(aerosol_optics_description_type), intent(inout) :: this
    character(len=*), intent(in)              :: file_name
    integer, intent(in), optional             :: iverbose
    
    ! The NetCDF file containing the aerosol optics data
    type(netcdf_file)  :: file

    real(jphook)       :: hook_handle

    if (lhook) call dr_hook('radiation_aerosol_optics_description:load',0,hook_handle)

    ! Open the aerosol scattering file and configure the way it is
    ! read
    call file%open(trim(file_name), iverbose=iverbose)

    ! Read metadata variables
    call file%get('code_hydrophilic', this%code_philic)
    call file%get('code_hydrophobic', this%code_phobic)
    call file%get('bin_hydrophilic',  this%bin_philic)
    call file%get('bin_hydrophobic',  this%bin_phobic)
    call file%get('optical_model_hydrophilic', this%optical_model_philic)
    call file%get('optical_model_hydrophobic', this%optical_model_phobic)

    ! Allocate logical arrays of the appropriate size and set to FALSE
    allocate(this%is_preferred_philic(size(this%bin_philic)))
    allocate(this%is_preferred_phobic(size(this%bin_phobic)))
    this%is_preferred_philic = .false.
    this%is_preferred_phobic = .false.

    call file%close()

    if (present(iverbose)) then
      this%iverbose = iverbose
    else
      this%iverbose = 3
    end if
    
    if (lhook) call dr_hook('radiation_aerosol_optics_description:load',1,hook_handle)

  end subroutine read

  !---------------------------------------------------------------------
  ! Specify the preferred optical model for a particular aerosol
  ! family, e.g. "call
  ! aer_desc%preferred_optical_model('DD','Woodward2001')" would mean
  ! that subsequent calls to get_index in which the optical model is
  ! not specified would return the Woodward model rather than the
  ! first matching model in the file.  The check is only done on the
  ! first len(optical_model_str) characters, so "Woodward" and
  ! "Woodward2001" would both match the Woodward2001 model.
  subroutine preferred_optical_model(this, code_str, optical_model_str)

    use yomhook,              only : lhook, dr_hook, jphook
    use radiation_io,         only : nulout, nulerr, radiation_abort
    
    class(aerosol_optics_description_type), intent(inout) :: this
    character(len=2), intent(in) :: code_str
    character(len=*), intent(in) :: optical_model_str

    ! Aerosol loop counter
    integer :: ja

    logical :: is_found, is_philic, is_phobic

    real(jphook)         :: hook_handle

    if (lhook) call dr_hook('radiation_aerosol_optics_description:preferred_optical_model',0,hook_handle)

    ! Check for empty string
    if (optical_model_str == ' ') then
      if (lhook) call dr_hook('radiation_aerosol_optics_description:preferred_optical_model',1,hook_handle)
      return
    end if

    is_found  = .false.
    is_philic = .false.
    is_phobic = .false.
    
    ! Loop over hydrophilic types
    do ja = 1,size(this%bin_philic)
      ! Check if we have a match
      if (to_string(this%code_philic(:,ja)) == code_str &
           &  .and. trim(to_string(this%optical_model_philic(:,ja))) &
           &          == optical_model_str) then
        this%is_preferred_philic(ja) = .true.
        is_found  = .true.
        is_philic = .true.
      end if
    end do
    ! Repeat for the hydrophobic types
    do ja = 1,size(this%bin_phobic)
      if (to_string(this%code_phobic(:,ja)) == code_str &
           &  .and. trim(to_string(this%optical_model_phobic(:,ja))) &
           &          == optical_model_str) then
        this%is_preferred_phobic(ja) = .true.
        is_found  = .true.
        is_phobic = .true.
      end if
    end do

    if (.not. is_found) then
      write(nulerr,'(a,a2,a,a,a)') '*** Error: Preferred "', code_str ,'" aerosol optical model "', &
           &  trim(optical_model_str), '" not found in file'
      call radiation_abort()
    else if (this%iverbose > 2) then
      write(nulout,'(a,a2,a,a,a)',advance='no') 'Preferred "', code_str, '" aerosol optical model set to "', &
           &  trim(optical_model_str), '" ('
      if (is_phobic) then
        write(nulout,'(a)',advance='no') ' hydrophobic'
      end if
      if (is_philic) then
        write(nulout,'(a)',advance='no') ' hydrophilic'
      end if
      write(nulout,'(a)') ' )'
    end if
    
    if (lhook) call dr_hook('radiation_aerosol_optics_description:preferred_optical_model',1,hook_handle)

  end subroutine preferred_optical_model

  
  !---------------------------------------------------------------------
  ! Return the index to the aerosol optical properties corresponding
  ! to the aerosol family in code_str (e.g. SS, DD etc), whether or
  ! not the requested aerosol is hydrophilic in the logical
  ! lhydrophilic, and optionally the size bin ibin and optical model
  ! in optical_model_str. The return value may be used to populate the
  ! radiation_config%i_aerosol_map vector, where a positive number is
  ! a hydrophobic index, a negative number is a hydrophilic index and
  ! zero indicates that the aerosol type was not found in the file.
  ! This is a valid entry in i_aerosol_map meaning the aerosol is
  ! ignored, but the calling routine to get_index might wish to throw
  ! a warning or error. This routine works by assigning a score based
  ! on the closeness of the match.
  function get_index(this, code_str, lhydrophilic, ibin, optical_model_str)
    
    use yomhook,              only : lhook, dr_hook, jphook
    use easy_netcdf,          only : netcdf_file
    use radiation_io,         only : nulout

    class(aerosol_optics_description_type), intent(in) :: this
    character(len=2), intent(in) :: code_str
    logical, intent(in) :: lhydrophilic
    integer, intent(in), optional :: ibin
    character(len=*), intent(in), optional :: optical_model_str

    ! Score of the currently selected aerosol index, and the score of
    ! the current one under consideration
    integer :: score, current_score

    ! Loop index for aerosol type
    integer :: ja

    ! Return value
    integer :: get_index

    ! Issue a warning if there is more than one equal match
    logical :: is_ambiguous

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_aerosol_optics_description:get_index',0,hook_handle)

    ! Initial values
    get_index = 0
    score = 0
    is_ambiguous = .false.

    if (lhydrophilic) then
      ! Loop over hydrophilic aerosol types
      do ja = 1,size(this%bin_philic)
        current_score = 0
        if (to_string(this%code_philic(:,ja)) == code_str) then
          ! Aerosol code matches
          if (present(ibin) .and. this%bin_philic(ja) > 0) then
            if (ibin > 0) then
              if (ibin == this%bin_philic(ja)) then
                ! Requested bin number matches
                current_score = 4
              else
                ! Requested bin number does not match
                current_score = -1
              end if
            else
              ! Bin number is zero: no request
              current_score = 2
            end if
          else
            ! No bin number present
            current_score = 2
          end if
          if (present(optical_model_str)) then
            if (trim(to_string(this%optical_model_philic(:,ja))) &
                 &  == optical_model_str) then
              ! Requested optical model matches
              if (current_score >= 0) then
                current_score = current_score + 4
              end if
            else
              ! Requested optical model does not match
              current_score = -1
            end if
          else if (current_score >= 0) then
            ! No requested optical model
            current_score = current_score + 2
          end if
          if (current_score > 0 .and. this%is_preferred_philic(ja)) then
            current_score = current_score + 1
          end if
          if (current_score > score) then
            ! Better score than any existing aerosol type
            get_index = -ja
            score = current_score
            is_ambiguous = .false.
          else if (current_score > 0 .and. current_score == score) then
            is_ambiguous = .true.
          end if
        end if
      end do
    else
      ! Loop over hydrophobic aerosol types
      do ja = 1,size(this%bin_phobic)
        current_score = 0
        if (to_string(this%code_phobic(:,ja)) == code_str) then
          ! Aerosol code matches
          if (present(ibin) .and. this%bin_phobic(ja) > 0) then
            if (ibin > 0) then
              if (ibin == this%bin_phobic(ja)) then
                ! Requested bin number matches
                current_score = 4
              else
                ! Requested bin number does not match
                current_score = -1
              end if
            else
              ! Bin number is zero: no request
              current_score = 2
            end if
          else
            ! No bin number requested or present
            current_score = 2
          end if
          if (present(optical_model_str)) then
            if (trim(to_string(this%optical_model_phobic(:,ja))) &
                 &  == optical_model_str) then
              ! Requested optical model matches
              if (current_score >= 0) then
                current_score = current_score + 4
              end if
            else
              ! Requested optical model does not match
              current_score = -1
            end if
          else if (current_score >= 0) then
            ! No requested optical model
            current_score = current_score + 2
          end if
          if (current_score > 0 .and. this%is_preferred_phobic(ja)) then
            current_score = current_score + 1
          end if
          if (current_score > score) then
            ! Better score than any existing aerosol type
            get_index = ja
            score = current_score
            is_ambiguous = .false.
          else if (current_score > 0 .and. current_score == score) then
            is_ambiguous = .true.
          end if          
        end if
      end do
    end if

    if (is_ambiguous) then
      write(nulout,'(a,a2,a,l,a)') 'Warning: radiation_aerosol_optics_description:get_index("', &
           &  code_str, '",', lhydrophilic, &
           &  ',...) does not unambiguously identify an aerosol optical property index'
    end if
    
    if (lhook) call dr_hook('radiation_aerosol_optics_description:get_index',1,hook_handle)

  end function get_index

  !---------------------------------------------------------------------
  ! Utility function to convert an array of single characters to a
  ! character string (yes Fortran's string handling is a bit
  ! rubbish). We set NULL characters (ASCII code 0) returned from the
  ! NetCDF library to spaces, so that TRIM can remove them.
  pure function to_string(arr) result(str)
    character, intent(in)  :: arr(:)
    character(len=size(arr)) :: str
    integer :: jc
    do jc = 1,size(arr)
      if (ichar(arr(jc)) == 0) then
        ! Replace NULL character with a space
        str(jc:jc) = ' '
      else
        str(jc:jc) = arr(jc)
      end if
    end do
  end function to_string

end module radiation_aerosol_optics_description
