! radiation_gas.F90 - Derived type to store the gas mixing ratios
!
! (C) Copyright 2014- ECMWF.
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
! Modifications
!   2019-01-14  R. Hogan  Added out_of_physical_bounds routine

module radiation_gas

  use parkind1, only : jprb
  use radiation_gas_constants

  implicit none
  public

  ! Available units
  enum, bind(c)
    enumerator IMassMixingRatio, IVolumeMixingRatio
  end enum

  !---------------------------------------------------------------------
  ! This derived type describes the gaseous composition of the
  ! atmosphere; gases may be stored as part of a 3D array (if their
  ! variation with height/column is to be represented) or one 1D array
  ! (if they are to be assumed globally well-mixed).
  type gas_type
    ! Units of each stored gas (or 0 if not present)
    integer :: iunits(NMaxGases) = 0

    ! Scaling factor that should be applied to each stored gas to get
    ! a dimensionless result, e.g. if iunits=IVolumeMixingRatio then
    ! 1.0e-6 is used to indicate the units are actually PPMV: need to
    ! multiply by 1e-6 to get mol/mol.
    real(jprb) :: scale_factor(NMaxGases) = 1.0_jprb
    
    ! Mixing ratios of variable gases, dimensioned (ncol, nlev,
    ! NMaxGases)
    real(jprb), allocatable, dimension(:,:,:) :: mixing_ratio

    ! Flag to indicate whether a gas is present
    logical :: is_present(NMaxGases) = .false.

    ! Flag to indicate whether a gas is well mixed
    logical :: is_well_mixed(NMaxGases) = .false.

    integer :: ntype          = 0 ! Number of gas types described

    integer :: ncol           = 0 ! Number of columns in mixing_ratio
    integer :: nlev           = 0 ! Number of levels  in mixing_ratio

    ! A list of length ntype of gases whose volume mixing ratios have
    ! been provided
    integer :: icode(NMaxGases) = 0
    
   contains
     procedure :: allocate   => allocate_gas
     procedure :: deallocate => deallocate_gas
     procedure :: put        => put_gas
     procedure :: put_well_mixed => put_well_mixed_gas
     procedure :: scale      => scale_gas
     procedure :: set_units  => set_units_gas
     procedure :: assert_units => assert_units_gas
     procedure :: get        => get_gas
     procedure :: reverse    => reverse_gas
     procedure :: out_of_physical_bounds
  end type gas_type

contains


  !---------------------------------------------------------------------
  ! Allocate a derived type for holding gas mixing ratios given the
  ! number of columns and levels
  subroutine allocate_gas(this, ncol, nlev)

    use yomhook, only : lhook, dr_hook, jphook

    class(gas_type), intent(inout) :: this
    integer,         intent(in)    :: ncol, nlev

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_gas:allocate',0,hook_handle)

    call this%deallocate()

    allocate(this%mixing_ratio(ncol, nlev, NMaxGases))
    this%mixing_ratio = 0.0_jprb

    this%ncol = ncol
    this%nlev = nlev

    if (lhook) call dr_hook('radiation_gas:allocate',1,hook_handle)

  end subroutine allocate_gas


  !---------------------------------------------------------------------
  ! Deallocate memory and reset arrays
  subroutine deallocate_gas(this)

    use yomhook, only : lhook, dr_hook, jphook

    class(gas_type), intent(inout) :: this

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_gas:deallocate',0,hook_handle)

    if (allocated(this%mixing_ratio)) then
       deallocate(this%mixing_ratio)
    end if
    
    this%iunits = 0
    this%scale_factor = 0.0_jprb
    this%is_present = .false.
    this%is_well_mixed = .false.
    this%ntype = 0
    this%ncol = 0
    this%nlev = 0
    this%icode = 0

    if (lhook) call dr_hook('radiation_gas:deallocate',1,hook_handle)

  end subroutine deallocate_gas


  !---------------------------------------------------------------------
  ! Put gas mixing ratio corresponding to gas ID "igas" with units
  ! "iunits"
  subroutine put_gas(this, igas, iunits, mixing_ratio, scale_factor, &
       istartcol)

    use yomhook,        only : lhook, dr_hook, jphook
    use radiation_io,   only : nulerr, radiation_abort

    class(gas_type),      intent(inout) :: this
    integer,              intent(in)    :: igas
    integer,              intent(in)    :: iunits
    real(jprb),           intent(in)    :: mixing_ratio(:,:)
    real(jprb), optional, intent(in)    :: scale_factor
    integer,    optional, intent(in)    :: istartcol

    integer :: i1, i2, jc, jk


    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_gas:put',0,hook_handle)

    ! Check inputs
    if (igas <= IGasNotPresent .or. iunits > NMaxGases) then
      write(nulerr,'(a,i0,a,i0,a,i0)') '*** Error: provided gas ID (', &
           &   igas, ') must be in the range ', IGasNotPresent+1, ' to ', &
           &   NMaxGases
      call radiation_abort()
    end if
    if (iunits < IMassMixingRatio .or. iunits > IVolumeMixingRatio) then
      write(nulerr,'(a,i0,a,i0,a,i0)') '*** Error: provided gas units (', &
           &   iunits, ') must be in the range ', IMassMixingRatio, ' to ', &
           &   IVolumeMixingRatio
      call radiation_abort()
    end if

    if (.not. allocated(this%mixing_ratio)) then
      write(nulerr,'(a,i0,a,i0,a,i0)') '*** Error: attempt to put data to unallocated radiation_gas object'
      call radiation_abort()
    end if

    if (present(istartcol)) then
      i1 = istartcol
    else
      i1 = 1
    end if

    i2 = i1 + size(mixing_ratio,1) - 1

    if (i1 < 1 .or. i2 < 1 .or. i1 > this%ncol .or. i2 > this%ncol) then
      write(nulerr,'(a,i0,a,i0,a,i0)') '*** Error: attempt to put columns indexed ', &
           &   i1, ' to ', i2, ' to array indexed 1 to ', this%ncol
      call radiation_abort()
    end if

    if (size(mixing_ratio,2) /= this%nlev) then
      write(nulerr,'(a,i0,a)') &
           &  '*** Error: gas mixing ratio expected to have ', this%nlev, &
           &  ' levels'
      call radiation_abort()
    end if

    if (.not. this%is_present(igas)) then
      ! Gas not present until now
      this%ntype = this%ntype + 1
      this%icode(this%ntype) = igas
    end if
    this%is_present(igas) = .true.
    this%iunits(igas) = iunits
    this%is_well_mixed(igas) = .false.

    do jk = 1,this%nlev
      do jc = i1,i2
        this%mixing_ratio(jc,jk,igas) = mixing_ratio(jc-i1+1,jk)
      end do
    end do
    if (present(scale_factor)) then
      this%scale_factor(igas) = scale_factor
    else
      this%scale_factor(igas) = 1.0_jprb
    end if

    if (lhook) call dr_hook('radiation_gas:put',1,hook_handle)

  end subroutine put_gas


  !---------------------------------------------------------------------
  ! Put well-mixed gas mixing ratio corresponding to gas ID "igas"
  ! with units "iunits"
  subroutine put_well_mixed_gas(this, igas, iunits, mixing_ratio, &
       scale_factor, istartcol, iendcol)

    use yomhook,        only : lhook, dr_hook, jphook
    use radiation_io,   only : nulerr, radiation_abort

    class(gas_type),      intent(inout) :: this
    integer,              intent(in)    :: igas
    integer,              intent(in)    :: iunits
    real(jprb),           intent(in)    :: mixing_ratio
    real(jprb), optional, intent(in)    :: scale_factor
    integer,    optional, intent(in)    :: istartcol, iendcol

    real(jphook) :: hook_handle

    integer :: i1, i2, jc, jk

    if (lhook) call dr_hook('radiation_gas:put_well_mixed',0,hook_handle)

    ! Check inputs
    if (igas <= IGasNotPresent .or. igas > NMaxGases) then
      write(nulerr,'(a,i0,a,i0,a,i0)') '*** Error: provided gas ID (', &
           &   igas, ') must be in the range ', IGasNotPresent+1, ' to ', &
           &   NMaxGases
      call radiation_abort()
    end if
    if (iunits < IMassMixingRatio .or. iunits > IVolumeMixingRatio) then
      write(nulerr,'(a,i0,a,i0,a,i0)') '*** Error: provided gas units (', &
           &   iunits, ') must be in the range ', IMassMixingRatio, ' to ', &
           &   IVolumeMixingRatio
      call radiation_abort()
    end if

    if (.not. allocated(this%mixing_ratio)) then
      write(nulerr,'(a)') '*** Error: attempt to put well-mixed gas data to unallocated radiation_gas object'
      call radiation_abort()
    end if

    if (present(istartcol)) then
      i1 = istartcol
    else
      i1 = 1
    end if

    if (present(iendcol)) then
      i2 = iendcol
    else
      i2 = this%ncol
    end if

    if (i1 < 1 .or. i2 < 1 .or. i1 > this%ncol .or. i2 > this%ncol) then
      write(nulerr,'(a,i0,a,i0,a,i0)') '*** Error: attempt to put columns indexed ', &
           &   i1, ' to ', i2, ' to array indexed 1 to ', this%ncol
      call radiation_abort()
    end if

    if (.not. this%is_present(igas)) then
      ! Gas not present until now
      this%ntype = this%ntype + 1
      this%icode(this%ntype) = igas
    end if
    ! Map uses a negative value to indicate a well-mixed value
    this%is_present(igas)              = .true.
    this%iunits(igas)                  = iunits
    this%is_well_mixed(igas)           = .true.

    do jk = 1,this%nlev
      do jc = i1,i2
        this%mixing_ratio(jc,jk,igas) = mixing_ratio
      end do
    end do
    if (present(scale_factor)) then
      this%scale_factor(igas) = scale_factor
    else
      this%scale_factor(igas) = 1.0_jprb
    end if

    if (lhook) call dr_hook('radiation_gas:put_well_mixed',1,hook_handle)

  end subroutine put_well_mixed_gas


  !---------------------------------------------------------------------
  ! Scale gas concentrations, e.g. igas=ICO2 and set scale_factor=2 to
  ! double CO2.  Note that this does not perform the scaling
  ! immediately, but changes the scale factor for the specified gas,
  ! ready to be used in set_units_gas.
  subroutine scale_gas(this, igas, scale_factor, lverbose)

    use radiation_io, only : nulout

    class(gas_type),      intent(inout) :: this
    integer,              intent(in)    :: igas
    real(jprb),           intent(in)    :: scale_factor
    logical,    optional, intent(in)    :: lverbose

    if (scale_factor /= 1.0_jprb) then
      this%scale_factor(igas) = this%scale_factor(igas) * scale_factor
      if (present(lverbose)) then
        if (lverbose) then
          write(nulout,'(a,a,a,f0.3)') '  Scaling ', trim(GasName(igas)), &
               &  ' concentration by ', scale_factor
        end if
      end if
    end if

  end subroutine scale_gas


  !---------------------------------------------------------------------
  ! Scale the gas concentrations so that they have the units "iunits"
  ! and are therefore ready to be used by the gas optics model within
  ! ecRad with no further scaling.  The existing scale_factor for each
  ! gas is applied.  If "igas" is present then apply only to gas with
  ! ID "igas", otherwise to all gases. Optional argument scale_factor
  ! specifies scaling that any subsequent access would need to apply
  ! to get a dimensionless result (consistent with definition of
  ! gas_type). So say that your gas optics model requires gas
  ! concentrations in PPMV, specify iunits=IVolumeMixingRatio and
  ! scale_factor=1.0e-6. If the gas concentrations were currently
  ! dimensionless volume mixing ratios, then the values would be
  ! internally divided by 1.0e-6.
  recursive subroutine set_units_gas(this, iunits, igas, scale_factor)
    class(gas_type),      intent(inout) :: this
    integer,              intent(in)    :: iunits
    integer,    optional, intent(in)    :: igas
    real(jprb), optional, intent(in)    :: scale_factor    

    integer :: ig

    ! Scaling factor to convert from old to new
    real(jprb) :: sf

    ! New scaling factor to store inside the gas object
    real(jprb) :: new_sf

    if (present(scale_factor)) then
      ! "sf" is the scaling to be applied now to the numbers (and may
      ! be modified below), "new_sf" is the value to be stored along
      ! with the numbers, informing subsequent routines how much you
      ! would need to multiply the numbers by to get a dimensionless
      ! result.
      sf     = 1.0_jprb / scale_factor
      new_sf = scale_factor
    else
      sf     = 1.0_jprb
      new_sf = 1.0_jprb
    end if

    if (present(igas)) then
      if (this%is_present(igas)) then
        if (iunits == IMassMixingRatio &
             &   .and. this%iunits(igas) == IVolumeMixingRatio) then
          sf = sf * GasMolarMass(igas) / AirMolarMass
        else if (iunits == IVolumeMixingRatio &
             &   .and. this%iunits(igas) == IMassMixingRatio) then
          sf = sf * AirMolarMass / GasMolarMass(igas)
        end if
        sf = sf * this%scale_factor(igas)
        
        if (sf /= 1.0_jprb) then
          this%mixing_ratio(:,:,igas) = this%mixing_ratio(:,:,igas) * sf
        end if
        ! Store the new units and scale factor for this gas inside the
        ! gas object
        this%iunits(igas) = iunits
        this%scale_factor(igas) = new_sf
      end if
    else
      do ig = 1,this%ntype
        call this%set_units(iunits, igas=this%icode(ig), scale_factor=new_sf)
      end do
    end if

  end subroutine set_units_gas


  !---------------------------------------------------------------------
  ! Assert that gas mixing ratio units are "iunits", applying to gas
  ! with ID "igas" if present, otherwise to all gases. Otherwise the
  ! program will exit. Otional argument scale factor specifies any
  ! subsequent multiplication to apply; for PPMV one would use
  ! iunits=IVolumeMixingRatio and scale_factor=1.0e6.
  recursive subroutine assert_units_gas(this, iunits, igas, scale_factor)
    
    use radiation_io,   only : nulerr, radiation_abort    

    class(gas_type),      intent(in) :: this
    integer,              intent(in) :: iunits
    integer,    optional, intent(in) :: igas
    real(jprb), optional, intent(in) :: scale_factor    

    integer :: ig

    real(jprb) :: sf

    if (present(scale_factor)) then
      sf = scale_factor
    else
      sf = 1.0_jprb
    end if

    if (present(igas)) then
      if (this%is_present(igas)) then
        if (iunits /= this%iunits(igas)) then
          write(nulerr,'(a,a,a)') '*** Error: ', trim(GasName(igas)), &
               &  ' is not in the required units'
          call radiation_abort()
        else if (sf /= this%scale_factor(igas)) then
          write(nulerr,'(a,a,a,e12.4,a,e12.4)') '*** Error: ', GasName(igas), &
               &  ' scaling of ', this%scale_factor(igas), &
               &  ' does not match required ', sf
          call radiation_abort()
        end if
      end if
    else
      do ig = 1,this%ntype
        call this%assert_units(iunits, igas=this%icode(ig), scale_factor=sf)
      end do
    end if

  end subroutine assert_units_gas


  !---------------------------------------------------------------------
  ! Get gas mixing ratio corresponding to gas ID "igas" with units
  ! "iunits" and return as a 2D array of dimensions (ncol,nlev).  The
  ! array will contain zeros if the gas is not stored.
  subroutine get_gas(this, igas, iunits, mixing_ratio, scale_factor, &
       &   istartcol)

    use yomhook,        only : lhook, dr_hook, jphook
    use radiation_io,   only : nulerr, radiation_abort

    class(gas_type),      intent(in)  :: this
    integer,              intent(in)  :: igas
    integer,              intent(in)  :: iunits
    real(jprb),           intent(out) :: mixing_ratio(:,:)
    real(jprb), optional, intent(in)  :: scale_factor
    integer,    optional, intent(in)  :: istartcol

    real(jprb)                        :: sf
    integer                           :: i1, i2

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_gas:get',0,hook_handle)

    if (present(scale_factor)) then
      sf = scale_factor
    else
      sf = 1.0_jprb
    end if

    if (present(istartcol)) then
      i1 = istartcol
    else
      i1 = 1
    end if

    i2 = i1 + size(mixing_ratio,1) - 1

    if (i1 < 1 .or. i2 < 1 .or. i1 > this%ncol .or. i2 > this%ncol) then
      write(nulerr,'(a,i0,a,i0,a,i0)') '*** Error: attempt to get columns indexed ', &
           &   i1, ' to ', i2, ' from array indexed 1 to ', this%ncol
      call radiation_abort()
    end if

    if (size(mixing_ratio,2) /= this%nlev) then
      write(nulerr,'(a,i0,a)') &
           &  '*** Error: gas destination array expected to have ', this%nlev, &
           &  ' levels'
      call radiation_abort()
    end if

    if (.not. this%is_present(igas)) then
      mixing_ratio = 0.0_jprb
    else 
      if (iunits == IMassMixingRatio &
           &   .and. this%iunits(igas) == IVolumeMixingRatio) then
        sf = sf * GasMolarMass(igas) / AirMolarMass
      else if (iunits == IVolumeMixingRatio &
           &   .and. this%iunits(igas) == IMassMixingRatio) then
        sf = sf * AirMolarMass / GasMolarMass(igas)
      end if
      sf = sf * this%scale_factor(igas)
        
      if (sf /= 1.0_jprb) then
        mixing_ratio = this%mixing_ratio(i1:i2,:,igas) * sf
      else
        mixing_ratio = this%mixing_ratio(i1:i2,:,igas)
      end if
    end if

    if (lhook) call dr_hook('radiation_gas:get',1,hook_handle)

  end subroutine get_gas


  !---------------------------------------------------------------------
  ! Copy data to "gas_rev", reversing the height ordering of the gas
  ! data
  subroutine reverse_gas(this, istartcol, iendcol, gas_rev)

    class(gas_type)             :: this
    integer,        intent(in)  :: istartcol, iendcol
    type(gas_type), intent(out) :: gas_rev

    gas_rev%iunits = this%iunits
    gas_rev%scale_factor = this%scale_factor
    gas_rev%is_present = this%is_present
    gas_rev%is_well_mixed = this%is_well_mixed
    gas_rev%ntype = this%ntype
    gas_rev%ncol = this%ncol
    gas_rev%nlev = this%nlev
    gas_rev%icode = this%icode

    if (allocated(gas_rev%mixing_ratio)) deallocate(gas_rev%mixing_ratio)

    if (allocated(this%mixing_ratio)) then
      allocate(gas_rev%mixing_ratio(istartcol:iendcol,this%nlev,NMaxGases))
      gas_rev%mixing_ratio(istartcol:iendcol,:,:) &
           &  = this%mixing_ratio(istartcol:iendcol,this%nlev:1:-1,:)
    end if

  end subroutine reverse_gas

  !---------------------------------------------------------------------
  ! Return .true. if variables are out of a physically sensible range,
  ! optionally only considering columns between istartcol and iendcol
  function out_of_physical_bounds(this, istartcol, iendcol, do_fix) result(is_bad)

    use yomhook,          only : lhook, dr_hook, jphook
    use radiation_check,  only : out_of_bounds_3d

    class(gas_type),   intent(inout) :: this
    integer,  optional,intent(in) :: istartcol, iendcol
    logical,  optional,intent(in) :: do_fix
    logical                       :: is_bad

    logical    :: do_fix_local

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_gas:out_of_physical_bounds',0,hook_handle)

    if (present(do_fix)) then
      do_fix_local = do_fix
    else
      do_fix_local = .false.
    end if

    is_bad = out_of_bounds_3d(this%mixing_ratio, 'gas%mixing_ratio', &
         &                    0.0_jprb, 1.0_jprb, do_fix_local, i1=istartcol, i2=iendcol)

    if (lhook) call dr_hook('radiation_gas:out_of_physical_bounds',1,hook_handle)

  end function out_of_physical_bounds
  
end module radiation_gas
