! radiation_thermodynamics.F90 - Derived type for pressure & temperature
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
!   2017-05-11  R. Hogan  Fix startcol/endcol for get_layer_mass
!   2019-01-14  R. Hogan  Added out_of_physical_bounds routine
!   2019-01-14  R. Hogan  Capped h2o_sat_liq at 1

module radiation_thermodynamics

  use parkind1, only : jprb

  implicit none
  public

  !---------------------------------------------------------------------
  ! Derived type for storing pressure and temperature at half levels
  type thermodynamics_type
     real(jprb), allocatable, dimension(:,:) :: &
          &  pressure_hl, &   ! (ncol,nlev+1) pressure (Pa)
          &  temperature_hl   ! (ncol,nlev+1) temperature (K)

     ! The following is a function of pressure and temperature: you
     ! can calculate it according to your favourite formula, or the
     ! calc_saturation_wrt_liquid subroutine can be used to do this
     ! approximately
     real(jprb), allocatable, dimension(:,:) :: &
          &  h2o_sat_liq ! (ncol,nlev) specific humidity at liquid
                         ! saturation (kg/kg)
   contains
     procedure :: allocate   => allocate_thermodynamics_arrays
     procedure :: deallocate => deallocate_thermodynamics_arrays
     procedure :: calc_saturation_wrt_liquid
     procedure :: get_layer_mass
     procedure :: get_layer_mass_column
     procedure :: out_of_physical_bounds

  end type thermodynamics_type

contains


  !---------------------------------------------------------------------
  ! Allocate variables with specified dimensions
  subroutine allocate_thermodynamics_arrays(this, ncol, nlev, &
       &                                    use_h2o_sat)

    use yomhook,  only : lhook, dr_hook, jphook

    class(thermodynamics_type), intent(inout) :: this
    integer, intent(in)           :: ncol  ! Number of columns
    integer, intent(in)           :: nlev  ! Number of levels
    logical, intent(in), optional :: use_h2o_sat ! Allocate h2o_sat_liq?

    logical :: use_h2o_sat_local

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_thermodynamics:allocate',0,hook_handle)

    allocate(this%pressure_hl(ncol,nlev+1))
    allocate(this%temperature_hl(ncol,nlev+1))

    use_h2o_sat_local = .false.
    if (present(use_h2o_sat)) then
      use_h2o_sat_local = use_h2o_sat
    end if
    
    if (use_h2o_sat_local) then
      allocate(this%h2o_sat_liq(ncol,nlev))
    end if    

    if (lhook) call dr_hook('radiation_thermodynamics:allocate',1,hook_handle)

  end subroutine allocate_thermodynamics_arrays


  !---------------------------------------------------------------------
  ! Deallocate variables
  subroutine deallocate_thermodynamics_arrays(this)

    use yomhook,  only : lhook, dr_hook, jphook

    class(thermodynamics_type), intent(inout) :: this

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_thermodynamics:deallocate',0,hook_handle)

    if (allocated(this%pressure_hl)) then
      deallocate(this%pressure_hl)
    end if
    if (allocated(this%temperature_hl)) then
      deallocate(this%temperature_hl)
    end if
    if (allocated(this%h2o_sat_liq)) then
      deallocate(this%h2o_sat_liq)
    end if

    if (lhook) call dr_hook('radiation_thermodynamics:deallocate',1,hook_handle)
  
  end subroutine deallocate_thermodynamics_arrays


  !---------------------------------------------------------------------
  ! Calculate approximate saturation with respect to liquid
  subroutine calc_saturation_wrt_liquid(this,istartcol,iendcol)

    use yomhook,  only : lhook, dr_hook, jphook

    class(thermodynamics_type), intent(inout) :: this
    integer, intent(in)                       :: istartcol, iendcol

    ! Pressure and temperature at full levels
    real(jprb) :: pressure, temperature

    ! Vapour pressure (Pa)
    real(jprb) :: e_sat

    integer :: ncol, nlev ! Dimension sizes
    integer :: jcol, jlev ! Loop indices for column and level

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_thermodynamics:calc_saturation_wrt_liquid',0,hook_handle)

    ncol = size(this%pressure_hl,1)
    nlev = size(this%pressure_hl,2) - 1

    if (.not. allocated(this%h2o_sat_liq)) then
      allocate(this%h2o_sat_liq(ncol,nlev))
    end if

    do jlev = 1,nlev
       do jcol = istartcol,iendcol
          pressure = 0.5 * (this%pressure_hl(jcol,jlev)+this%pressure_hl(jcol,jlev+1))
          temperature = 0.5 * (this%temperature_hl(jcol,jlev)+this%temperature_hl(jcol,jlev+1))
          e_sat = 6.11e2_jprb * exp( 17.269_jprb * (temperature-273.16_jprb) / (temperature-35.86_jprb) )
          ! This formula can go above 1 at low pressure so needs to be
          ! capped
          this%h2o_sat_liq(jcol,jlev) = min(1.0_jprb, 0.622_jprb * e_sat / pressure)
       end do
    end do

    if (lhook) call dr_hook('radiation_thermodynamics:calc_saturation_wrt_liquid',1,hook_handle)

  end subroutine calc_saturation_wrt_liquid


  !---------------------------------------------------------------------
  ! Calculate the dry mass of each layer, neglecting humidity effects.
  ! The first version is for all columns.
  subroutine get_layer_mass(this,istartcol,iendcol,layer_mass)

    use yomhook,              only : lhook, dr_hook, jphook
    use radiation_constants,  only : AccelDueToGravity

    class(thermodynamics_type), intent(in)  :: this
    integer,                    intent(in)  :: istartcol, iendcol
    real(jprb),                 intent(out) :: layer_mass(:,:)

    integer    :: nlev
    real(jprb) :: inv_g

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_thermodynamics:get_layer_mass',0,hook_handle)

    nlev  = ubound(this%pressure_hl,2) - 1
    inv_g = 1.0_jprb / AccelDueToGravity

    layer_mass(istartcol:iendcol,1:nlev) &
         &  = ( this%pressure_hl(istartcol:iendcol,2:nlev+1) &
         &     -this%pressure_hl(istartcol:iendcol,1:nlev  )  ) &
         &  * inv_g 
    
    if (lhook) call dr_hook('radiation_thermodynamics:get_layer_mass',1,hook_handle)

  end subroutine get_layer_mass

  !---------------------------------------------------------------------
  ! Calculate the dry mass of each layer, neglecting humidity effects.
  ! The second version is for one column, the one numbered "icol".
  subroutine get_layer_mass_column(this, icol, layer_mass)

    use yomhook,              only : lhook, dr_hook, jphook
    use radiation_constants,  only : AccelDueToGravity

    class(thermodynamics_type), intent(in)  :: this
    integer,                    intent(in)  :: icol
    real(jprb),                 intent(out) :: layer_mass(:)

    integer    :: nlev
    real(jprb) :: inv_g

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_thermodynamics:get_layer_mass_column',0,hook_handle)

    nlev  = ubound(this%pressure_hl,2) - 1
    inv_g = 1.0_jprb / AccelDueToGravity

    layer_mass = ( this%pressure_hl(icol,2:nlev+1) &
             &    -this%pressure_hl(icol,1:nlev  )  ) &
             &   * inv_g
    
    if (lhook) call dr_hook('radiation_thermodynamics:get_layer_mass_column',1,hook_handle)

  end subroutine get_layer_mass_column


  !---------------------------------------------------------------------
  ! Estimate the separation between the mid-points of model layers
  ! given the half-level pressure and temperature.  This is not in
  ! terms of the "thermodynamics" type as it is useful for computing
  ! overlap decorrelation lengths and hence cloud cover outside the
  ! radiation scheme.
  subroutine get_layer_separation(pressure_hl, temperature_hl, layer_separation)

    use yomhook,              only : lhook, dr_hook, jphook
    use radiation_constants,  only : GasConstantDryAir, AccelDueToGravity

    ! Pressure (Pa) and temperature (K) at half-levels, dimensioned
    ! (ncol,nlev+1) where ncol is the number of columns and nlev is
    ! the number of model levels
    real(jprb), dimension(:,:), intent(in)  :: pressure_hl, temperature_hl

    ! Layer separation in metres, dimensioned (ncol,nlev-1)
    real(jprb), dimension(:,:), intent(out) :: layer_separation

    ! Ratio of gas constant for dry air to acceleration due to gravity
    real(jprb), parameter :: R_over_g = GasConstantDryAir / AccelDueToGravity

    ! Loop indices and array bounds
    integer    :: jlev
    integer    :: i1, i2, nlev

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_thermodynamics:get_layer_separation',0,hook_handle)

    i1   = lbound(pressure_hl,1)
    i2   = ubound(pressure_hl,1)
    nlev =   size(pressure_hl,2)-1

    if (pressure_hl(i1,2) > pressure_hl(i1,1)) then
      ! Pressure is increasing with index (order of layers is
      ! top-of-atmosphere to surface). In case pressure_hl(:,1)=0, we
      ! don't take the logarithm of the first pressure in each column.
      layer_separation(i1:i2,1) = R_over_g * temperature_hl(i1:i2,2) &
           &                    * log(pressure_hl(i1:i2,3)/pressure_hl(i1:i2,2))
      
      ! For other layers we take the separation between midpoints to
      ! be half the separation between the half-levels at the edge of
      ! the two adjacent layers
      do jlev = 2,nlev-1
        layer_separation(i1:i2,jlev) = (0.5_jprb * R_over_g) * temperature_hl(i1:i2,jlev+1) &
             &                    * log(pressure_hl(i1:i2,jlev+2)/pressure_hl(i1:i2,jlev))

      end do

    else
      ! Pressure is decreasing with index (order of layers is surface
      ! to top-of-atmosphere).  In case pressure_hl(:,nlev+1)=0, we
      ! don't take the logarithm of the last pressure in each column.

      do jlev = 1,nlev-2
        layer_separation(i1:i2,jlev) = (0.5_jprb * R_over_g) * temperature_hl(i1:i2,jlev+1) &
             &                    * log(pressure_hl(i1:i2,jlev)/pressure_hl(i1:i2,jlev+2))

      end do
      layer_separation(i1:i2,nlev-1) = R_over_g * temperature_hl(i1:i2,nlev) &
           &                    * log(pressure_hl(i1:i2,nlev-1)/pressure_hl(i1:i2,nlev))

    end if

    if (lhook) call dr_hook('radiation_thermodynamics:get_layer_separation',1,hook_handle)    

  end subroutine get_layer_separation


  !---------------------------------------------------------------------
  ! Return .true. if variables are out of a physically sensible range,
  ! optionally only considering columns between istartcol and iendcol
  function out_of_physical_bounds(this, istartcol, iendcol, do_fix) result(is_bad)

    use yomhook,          only : lhook, dr_hook, jphook
    use radiation_check,  only : out_of_bounds_2d

    class(thermodynamics_type), intent(inout) :: this
    integer,           optional,intent(in) :: istartcol, iendcol
    logical,           optional,intent(in) :: do_fix
    logical                                :: is_bad

    logical    :: do_fix_local

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_thermodynamics:out_of_physical_bounds',0,hook_handle)

    if (present(do_fix)) then
      do_fix_local = do_fix
    else
      do_fix_local = .false.
    end if

    ! Dangerous to cap pressure_hl as then the pressure difference across a layer could be zero
    is_bad =    out_of_bounds_2d(this%pressure_hl, 'pressure_hl', 0.0_jprb, 110000.0_jprb, &
         &                       .false., i1=istartcol, i2=iendcol) &
         & .or. out_of_bounds_2d(this%temperature_hl, 'temperature_hl', 100.0_jprb,  400.0_jprb, &
         &                       do_fix_local, i1=istartcol, i2=iendcol) &
         & .or. out_of_bounds_2d(this%h2o_sat_liq, 'h2o_sat_liq', 0.0_jprb, 1.0_jprb, &
         &                       do_fix_local, i1=istartcol, i2=iendcol)

    if (lhook) call dr_hook('radiation_thermodynamics:out_of_physical_bounds',1,hook_handle)

  end function out_of_physical_bounds
  
end module radiation_thermodynamics
