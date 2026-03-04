! radiation_flux.F90 - Derived type to store the output fluxes
!
! Copyright (C) 2014-2019 ECMWF
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
! License: see the COPYING file for details
!
! Modifications
!   2017-09-08  R. Hogan  Store g-point fluxes
!   2017-10-23  R. Hogan  Renamed single-character variables
!   2019-01-08  R. Hogan  Added "indexed_sum_profile"
!   2019-01-14  R. Hogan  out_of_physical_bounds calls routine in radiation_config

module radiation_flux

  use parkind1, only : jprb

  implicit none

  !---------------------------------------------------------------------
  ! This derived type contains the output from the radiation
  ! calculation.  Currently this is solely flux profiles, but in
  ! future surface fluxes in each band may be stored in order that the
  ! calling program can compute surface-radiation such as
  ! photosynthetically active radiation and UV index.
  type flux_type
     ! All the following are broad-band fluxes in W m-2 with
     ! dimensions (ncol,nlev+1).  Note that only those fluxes that are
     ! requested will be used, so clear-sky and direct-beam arrays may
     ! not be allocated
     real(jprb), allocatable, dimension(:,:) :: &
          &  lw_up, lw_dn, &   ! Upwelling and downwelling longwave
          &  sw_up, sw_dn, &   ! Upwelling and downwelling shortwave
          &  sw_dn_direct, &   ! Direct-beam shortwave into a horizontal plane
          &  lw_up_clear, lw_dn_clear, & ! Clear-sky quantities...
          &  sw_up_clear, sw_dn_clear, &
          &  sw_dn_direct_clear
     ! As above but fluxes in each spectral band in W m-2 with
     ! dimensions (nband,ncol,nlev+1).  These are only allocated if
     ! config%do_save_spectral_flux==.true.
     real(jprb), allocatable, dimension(:,:,:) :: &
          &  lw_up_band, lw_dn_band, &   ! Upwelling and downwelling longwave
          &  sw_up_band, sw_dn_band, &   ! Upwelling and downwelling shortwave
          &  sw_dn_direct_band, &        ! Direct-beam shortwave
          &  lw_up_clear_band, lw_dn_clear_band, & ! Clear-sky quantities...
          &  sw_up_clear_band, sw_dn_clear_band, &
          &  sw_dn_direct_clear_band
     ! Surface downwelling quantaties at each g point, dimensioned
     ! (ng,ncol), that are always saved by the solver, except for the
     ! clear-sky ones that are only produced if
     ! config%do_clear==.true.
     real(jprb), allocatable, dimension(:,:) :: &
          &  lw_dn_surf_g, lw_dn_surf_clear_g, &
          &  sw_dn_diffuse_surf_g, sw_dn_direct_surf_g, &
          &  sw_dn_diffuse_surf_clear_g, sw_dn_direct_surf_clear_g
     ! Shortwave downwelling spectral fluxes in W m-2 at the surface,
     ! from which quantities such as photosynthetically active and UV
     ! radiation can be computed. Only allocated in
     ! config%do_surface_sw_spectral_flux==.true.  Note that the
     ! clear-sky quantities are only computed if
     ! config%do_clear==.true., but direct fluxes are computed whether
     ! or not do_direct==.true.. The dimensions are (nband,ncol).
     real(jprb), allocatable, dimension(:,:) :: &
          &  sw_dn_surf_band, sw_dn_direct_surf_band, &
          &  sw_dn_surf_clear_band, sw_dn_direct_surf_clear_band
     ! Surface downwelling fluxes in W m-2 at the spectral resolution
     ! needed by any subsequent canopy radiative transfer.  If
     ! config%use_canopy_full_spectrum_[sw|lw] then these will be at
     ! g-point resolution; otherwise they will be at
     ! config%n_albedo_bands and config%n_emiss_bands resolution.
     real(jprb), allocatable, dimension(:,:) :: &
          &  lw_dn_surf_canopy, &
          &  sw_dn_diffuse_surf_canopy, sw_dn_direct_surf_canopy

     ! Diagnosed cloud cover from the short- and long-wave solvers
     real(jprb), allocatable, dimension(:) :: &
          &  cloud_cover_lw, cloud_cover_sw
     ! Longwave derivatives needed by Hogan and Bozzo (2015) method
     ! for approximate longwave updates in between the full radiation
     ! calls: rate of change of upwelling broad-band flux with respect
     ! to surface value, dimensioned (ncol,nlev+1)
     real(jprb), allocatable, dimension(:,:) :: &
          &  lw_derivatives

   contains
     procedure :: allocate   => allocate_flux_type
     procedure :: deallocate => deallocate_flux_type
     procedure :: calc_surface_spectral
     procedure :: out_of_physical_bounds
  end type flux_type

contains

  !---------------------------------------------------------------------
  ! Allocate arrays for flux profiles, using config to define which
  ! fluxes are needed.  The arrays are dimensioned for columns between
  ! istartcol, iendcol and levels from 1 to nlev+1
  subroutine allocate_flux_type(this, config, istartcol, iendcol, nlev)

    use yomhook,  only           : lhook, dr_hook, jphook
    use radiation_io,     only : nulerr, radiation_abort
    use radiation_config, only : config_type

    integer, intent(in)             :: istartcol, iendcol, nlev
    class(flux_type), intent(inout) :: this
    type(config_type), intent(in)   :: config

    real(jphook)                    :: hook_handle

    if (lhook) call dr_hook('radiation_flux:allocate',0,hook_handle)

    ! Allocate longwave arrays
    if (config%do_lw) then
      allocate(this%lw_up(istartcol:iendcol,nlev+1))
      allocate(this%lw_dn(istartcol:iendcol,nlev+1))
      if (config%do_clear) then
        allocate(this%lw_up_clear(istartcol:iendcol,nlev+1))
        allocate(this%lw_dn_clear(istartcol:iendcol,nlev+1))
      end if
      
      if (config%do_save_spectral_flux) then
        if (config%n_spec_lw == 0) then
          write(nulerr,'(a)') '*** Error: number of LW spectral points to save not yet defined ' &
               & // 'so cannot allocated spectral flux arrays'
          call radiation_abort()
        end if
        
        allocate(this%lw_up_band(config%n_spec_lw,istartcol:iendcol,nlev+1))
        allocate(this%lw_dn_band(config%n_spec_lw,istartcol:iendcol,nlev+1))
        if (config%do_clear) then
          allocate(this%lw_up_clear_band(config%n_spec_lw, &
               &                         istartcol:iendcol,nlev+1))
          allocate(this%lw_dn_clear_band(config%n_spec_lw, &
               &                         istartcol:iendcol,nlev+1))
        end if
      end if
      
      if (config%do_lw_derivatives) then
        allocate(this%lw_derivatives(istartcol:iendcol,nlev+1))
      end if

      ! Allocate g-point downwelling fluxes at surface passed from
      ! solver to surface_intermediate%partition
      allocate(this%lw_dn_surf_g(config%n_g_lw,istartcol:iendcol))
      if (config%do_clear) then
        allocate(this%lw_dn_surf_clear_g(config%n_g_lw,istartcol:iendcol))
      end if

      if (config%do_canopy_fluxes_lw) then
        ! Downward fluxes at top of canopy at the spectral resolution
        ! used in the canopy radiative transfer scheme
        allocate(this%lw_dn_surf_canopy(config%n_canopy_bands_lw,istartcol:iendcol))
      end if
    end if
    
    ! Allocate shortwave arrays
    if (config%do_sw) then
      allocate(this%sw_up(istartcol:iendcol,nlev+1))
      allocate(this%sw_dn(istartcol:iendcol,nlev+1))
      if (config%do_sw_direct) then
        allocate(this%sw_dn_direct(istartcol:iendcol,nlev+1))
      end if
      if (config%do_clear) then
        allocate(this%sw_up_clear(istartcol:iendcol,nlev+1))
        allocate(this%sw_dn_clear(istartcol:iendcol,nlev+1))
        if (config%do_sw_direct) then
          allocate(this%sw_dn_direct_clear(istartcol:iendcol,nlev+1))
        end if
      end if
      
      if (config%do_save_spectral_flux) then
        if (config%n_spec_sw == 0) then
          write(nulerr,'(a)') '*** Error: number of SW spectral points to save not yet defined ' &
               & // 'so cannot allocate spectral flux arrays'
          call radiation_abort()
        end if
        
        allocate(this%sw_up_band(config%n_spec_sw,istartcol:iendcol,nlev+1))
        allocate(this%sw_dn_band(config%n_spec_sw,istartcol:iendcol,nlev+1))
        
        if (config%do_sw_direct) then
          allocate(this%sw_dn_direct_band(config%n_spec_sw, &
               &                          istartcol:iendcol,nlev+1))
        end if
        if (config%do_clear) then
          allocate(this%sw_up_clear_band(config%n_spec_sw, &
               &                         istartcol:iendcol,nlev+1))
          allocate(this%sw_dn_clear_band(config%n_spec_sw, &
               &                         istartcol:iendcol,nlev+1))
          if (config%do_sw_direct) then
            allocate(this%sw_dn_direct_clear_band(config%n_spec_sw, &
                 &                                istartcol:iendcol, nlev+1))
          end if
        end if
      end if
      
      if (config%do_surface_sw_spectral_flux) then
        if (config%n_bands_sw == 0) then
          write(nulerr,'(a)') '*** Error: number of SW bands not yet defined ' &
               & // 'so cannot allocate surface spectral flux arrays'
          call radiation_abort()
        end if
        allocate(this%sw_dn_surf_band(config%n_bands_sw,istartcol:iendcol))
        allocate(this%sw_dn_direct_surf_band(config%n_bands_sw,istartcol:iendcol))
        if (config%do_clear) then
          allocate(this%sw_dn_surf_clear_band(config%n_bands_sw, &
               &                              istartcol:iendcol))
          allocate(this%sw_dn_direct_surf_clear_band(config%n_bands_sw, &
               &                                     istartcol:iendcol))
        end if
      end if

      ! Allocate g-point downwelling fluxes at surface passed from
      ! solver to surface_intermediate%partition
      allocate(this%sw_dn_diffuse_surf_g(config%n_g_sw,istartcol:iendcol))
      allocate(this%sw_dn_direct_surf_g (config%n_g_sw,istartcol:iendcol))
      if (config%do_clear) then
        allocate(this%sw_dn_diffuse_surf_clear_g(config%n_g_sw,istartcol:iendcol))
        allocate(this%sw_dn_direct_surf_clear_g (config%n_g_sw,istartcol:iendcol))
      end if

      if (config%do_canopy_fluxes_sw) then
        ! Downward fluxes at top of canopy at the spectral resolution
        ! used in the canopy radiative transfer scheme
        allocate(this%sw_dn_diffuse_surf_canopy(config%n_canopy_bands_sw,istartcol:iendcol))
        allocate(this%sw_dn_direct_surf_canopy (config%n_canopy_bands_sw,istartcol:iendcol))
      end if
    end if
    
    ! Allocate cloud cover arrays
    allocate(this%cloud_cover_lw(istartcol:iendcol))
    allocate(this%cloud_cover_sw(istartcol:iendcol))

    ! Some solvers may not write to cloud cover, so we initialize to
    ! an unphysical value
    this%cloud_cover_lw = -1.0_jprb
    this%cloud_cover_sw = -1.0_jprb

    if (lhook) call dr_hook('radiation_flux:allocate',1,hook_handle)
    
  end subroutine allocate_flux_type


  !---------------------------------------------------------------------
  ! Deallocate flux arrays
  subroutine deallocate_flux_type(this)

    use yomhook,  only           : lhook, dr_hook, jphook
    class(flux_type), intent(inout) :: this
    real(jphook)                    :: hook_handle

    if (lhook) call dr_hook('radiation_flux:deallocate',0,hook_handle)

    if (allocated(this%lw_up)) then
      deallocate(this%lw_up)
      if (allocated(this%lw_dn))       deallocate(this%lw_dn)
      if (allocated(this%lw_up_clear)) deallocate(this%lw_up_clear)
      if (allocated(this%lw_dn_clear)) deallocate(this%lw_dn_clear)
    end if

    if (allocated(this%sw_up)) then
      deallocate(this%sw_up)
      if (allocated(this%sw_dn))        deallocate(this%sw_dn)
      if (allocated(this%sw_up_clear))  deallocate(this%sw_up_clear)
      if (allocated(this%sw_dn_clear))  deallocate(this%sw_dn_clear)
      if (allocated(this%sw_dn_direct)) deallocate(this%sw_dn_direct)
      if (allocated(this%sw_dn_direct_clear)) &
           &   deallocate(this%sw_dn_direct_clear)
    end if

    if (allocated(this%lw_up_band)) then
      deallocate(this%lw_up_band)
      if (allocated(this%lw_dn_band))       deallocate(this%lw_dn_band)
      if (allocated(this%lw_up_clear_band)) deallocate(this%lw_up_clear_band)
      if (allocated(this%lw_dn_clear_band)) deallocate(this%lw_dn_clear_band)
    end if
    
    if (allocated(this%sw_up_band)) then
      deallocate(this%sw_up_band)
      if (allocated(this%sw_dn_band))        deallocate(this%sw_dn_band)
      if (allocated(this%sw_up_clear_band))  deallocate(this%sw_up_clear_band)
      if (allocated(this%sw_dn_clear_band))  deallocate(this%sw_dn_clear_band)
      if (allocated(this%sw_dn_direct_band)) deallocate(this%sw_dn_direct_band)
      if (allocated(this%sw_dn_direct_clear_band)) &
           &   deallocate(this%sw_dn_direct_clear_band)      
    end if

    if (allocated(this%sw_dn_surf_band)) then
      deallocate(this%sw_dn_surf_band)
      deallocate(this%sw_dn_direct_surf_band)
    end if
    if (allocated(this%sw_dn_surf_clear_band)) then
      deallocate(this%sw_dn_surf_clear_band)
      deallocate(this%sw_dn_direct_surf_clear_band)
    end if

    if (allocated(this%lw_dn_surf_canopy)) deallocate(this%lw_dn_surf_canopy)
    if (allocated(this%sw_dn_diffuse_surf_canopy)) deallocate(this%sw_dn_diffuse_surf_canopy)
    if (allocated(this%sw_dn_direct_surf_canopy)) deallocate(this%sw_dn_direct_surf_canopy)

    if (allocated(this%cloud_cover_sw)) then
      deallocate(this%cloud_cover_sw)
    end if
    if (allocated(this%cloud_cover_lw)) then
      deallocate(this%cloud_cover_lw)
    end if

    if (allocated(this%lw_derivatives)) then
      deallocate(this%lw_derivatives)
    end if

    if (lhook) call dr_hook('radiation_flux:deallocate',1,hook_handle)

  end subroutine deallocate_flux_type

  !---------------------------------------------------------------------
  ! Calculate surface downwelling fluxes in each band using the
  ! downwelling surface fluxes at each g point
  subroutine calc_surface_spectral(this, config, istartcol, iendcol)

    use yomhook,  only           : lhook, dr_hook, jphook
    use radiation_config, only : config_type

    class(flux_type),  intent(inout) :: this
    type(config_type), intent(in)    :: config
    integer,           intent(in)    :: istartcol, iendcol

    integer :: jcol, jband, jalbedoband, nalbedoband

    ! Longwave surface downwelling in each band needed to compute
    ! canopy fluxes
    real(jprb) :: lw_dn_surf_band(config%n_bands_lw,istartcol:iendcol)

    real(jphook)                     :: hook_handle

    if (lhook) call dr_hook('radiation_flux:calc_surface_spectral',0,hook_handle)

    if (config%do_sw .and. config%do_surface_sw_spectral_flux) then

      do jcol = istartcol,iendcol
        call indexed_sum(this%sw_dn_direct_surf_g(:,jcol), &
             &           config%i_band_from_reordered_g_sw, &
             &           this%sw_dn_direct_surf_band(:,jcol))
        call indexed_sum(this%sw_dn_diffuse_surf_g(:,jcol), &
             &           config%i_band_from_reordered_g_sw, &
             &           this%sw_dn_surf_band(:,jcol))
        this%sw_dn_surf_band(:,jcol) &
             &  = this%sw_dn_surf_band(:,jcol) &
             &  + this%sw_dn_direct_surf_band(:,jcol)
      end do

      if (config%do_clear) then
        do jcol = istartcol,iendcol
          call indexed_sum(this%sw_dn_direct_surf_clear_g(:,jcol), &
               &           config%i_band_from_reordered_g_sw, &
               &           this%sw_dn_direct_surf_clear_band(:,jcol))
          call indexed_sum(this%sw_dn_diffuse_surf_clear_g(:,jcol), &
               &           config%i_band_from_reordered_g_sw, &
               &           this%sw_dn_surf_clear_band(:,jcol))
          this%sw_dn_surf_clear_band(:,jcol) &
               &  = this%sw_dn_surf_clear_band(:,jcol) &
               &  + this%sw_dn_direct_surf_clear_band(:,jcol)
        end do
      end if

    end if ! do_surface_sw_spectral_flux

    ! Fluxes in bands required for canopy radiative transfer
    if (config%do_sw .and. config%do_canopy_fluxes_sw) then
      if (config%use_canopy_full_spectrum_sw) then
        this%sw_dn_diffuse_surf_canopy(:,istartcol:iendcol) = this%sw_dn_diffuse_surf_g(:,istartcol:iendcol)
        this%sw_dn_direct_surf_canopy (:,istartcol:iendcol) = this%sw_dn_direct_surf_g (:,istartcol:iendcol)
      else if (config%do_nearest_spectral_sw_albedo) then
        do jcol = istartcol,iendcol
          call indexed_sum(this%sw_dn_direct_surf_g(:,jcol), &
               &           config%i_albedo_from_band_sw(config%i_band_from_reordered_g_sw), &
               &           this%sw_dn_direct_surf_canopy(:,jcol))
          call indexed_sum(this%sw_dn_diffuse_surf_g(:,jcol), &
               &           config%i_albedo_from_band_sw(config%i_band_from_reordered_g_sw), &
               &           this%sw_dn_diffuse_surf_canopy(:,jcol))
        end do
      else
        ! More accurate calculations using weights, but requires
        ! this%sw_dn_[direct_]surf_band to be defined, i.e.
        ! config%do_surface_sw_spectral_flux == .true.
        nalbedoband = size(config%sw_albedo_weights,1)
        this%sw_dn_diffuse_surf_canopy(:,istartcol:iendcol) = 0.0_jprb
        this%sw_dn_direct_surf_canopy (:,istartcol:iendcol) = 0.0_jprb
        do jband = 1,config%n_bands_sw
          do jalbedoband = 1,nalbedoband
            if (config%sw_albedo_weights(jalbedoband,jband) /= 0.0_jprb) then
              ! Initially, "diffuse" is actually "total"
              this%sw_dn_diffuse_surf_canopy(jalbedoband,istartcol:iendcol) &
                   &  = this%sw_dn_diffuse_surf_canopy(jalbedoband,istartcol:iendcol) &
                   &  + config%sw_albedo_weights(jalbedoband,jband) &
                   &    * this%sw_dn_surf_band(jband,istartcol:iendcol)
              this%sw_dn_direct_surf_canopy(jalbedoband,istartcol:iendcol) &
                   &  = this%sw_dn_direct_surf_canopy(jalbedoband,istartcol:iendcol) &
                   &  + config%sw_albedo_weights(jalbedoband,jband) &
                   &    * this%sw_dn_direct_surf_band(jband,istartcol:iendcol)
            end if
          end do
        end do
        ! Subtract the direct from total to get diffuse
        this%sw_dn_diffuse_surf_canopy(:,istartcol:iendcol) &
             &  = this%sw_dn_diffuse_surf_canopy(:,istartcol:iendcol) &
             &  - this%sw_dn_direct_surf_canopy(:,istartcol:iendcol)
      end if

    end if ! do_canopy_fluxes_sw

    if (config%do_lw .and. config%do_canopy_fluxes_lw) then
      if (config%use_canopy_full_spectrum_lw) then
        this%lw_dn_surf_canopy(:,istartcol:iendcol) = this%lw_dn_surf_g(:,istartcol:iendcol)
      else if (config%do_nearest_spectral_lw_emiss) then
        do jcol = istartcol,iendcol
          call indexed_sum(this%lw_dn_surf_g(:,jcol), &
               &           config%i_emiss_from_band_lw(config%i_band_from_reordered_g_lw), &
               &           this%lw_dn_surf_canopy(:,jcol))
        end do
      else
        ! Compute fluxes in each longwave emissivity interval using
        ! weights; first sum over g points to get the values in bands
        do jcol = istartcol,iendcol
          call indexed_sum(this%lw_dn_surf_g(:,jcol), &
               &           config%i_band_from_reordered_g_lw, &
               &           lw_dn_surf_band(:,jcol))
        end do
        nalbedoband = size(config%lw_emiss_weights,1)
        this%lw_dn_surf_canopy(:,istartcol:iendcol) = 0.0_jprb
        do jband = 1,config%n_bands_lw
          do jalbedoband = 1,nalbedoband
            if (config%lw_emiss_weights(jalbedoband,jband) /= 0.0_jprb) then
              this%lw_dn_surf_canopy(jalbedoband,istartcol:iendcol) &
                   &  = this%lw_dn_surf_canopy(jalbedoband,istartcol:iendcol) &
                   &  + config%lw_emiss_weights(jalbedoband,jband) &
                   &    * lw_dn_surf_band(jband,istartcol:iendcol)
            end if
          end do
        end do
      end if
    end if

    if (lhook) call dr_hook('radiation_flux:calc_surface_spectral',1,hook_handle)

  end subroutine calc_surface_spectral


  !---------------------------------------------------------------------
  ! Return .true. if the most important flux variables are out of a
  ! physically sensible range, optionally only considering columns
  ! between istartcol and iendcol
  function out_of_physical_bounds(this, istartcol, iendcol) result(is_bad)

    use yomhook,  only           : lhook, dr_hook, jphook
    use radiation_config, only : out_of_bounds_2d

    class(flux_type), intent(inout) :: this
    integer, optional,intent(in) :: istartcol, iendcol
    logical                      :: is_bad

    real(jphook)                 :: hook_handle

    if (lhook) call dr_hook('radiation_flux:out_of_physical_bounds',0,hook_handle)

    is_bad =    out_of_bounds_2d(this%lw_up, 'lw_up', 10.0_jprb, 900.0_jprb, .false., i1=istartcol, i2=iendcol) &
         & .or. out_of_bounds_2d(this%lw_dn, 'lw_dn', 0.0_jprb,  800.0_jprb, .false., i1=istartcol, i2=iendcol) &
         & .or. out_of_bounds_2d(this%sw_up, 'sw_up', 0.0_jprb, 1500.0_jprb, .false., i1=istartcol, i2=iendcol) &
         & .or. out_of_bounds_2d(this%sw_dn, 'sw_dn', 0.0_jprb, 1500.0_jprb, .false., i1=istartcol, i2=iendcol) &
         & .or. out_of_bounds_2d(this%sw_dn_direct, 'sw_dn_direct', 0.0_jprb, 1500.0_jprb, .false., i1=istartcol, i2=iendcol) &
         & .or. out_of_bounds_2d(this%lw_derivatives, 'lw_derivatives', 0.0_jprb, 1.0_jprb, .false., i1=istartcol, i2=iendcol) &
         & .or. out_of_bounds_2d(this%sw_dn_surf_band, 'sw_dn_surf_band', 0.0_jprb, 1500.0_jprb, &
         &                       .false., j1=istartcol, j2=iendcol) &
         & .or. out_of_bounds_2d(this%sw_dn_surf_clear_band, 'sw_dn_surf_clear_band', 0.0_jprb, 1500.0_jprb, &
         &                       .false., j1=istartcol, j2=iendcol)

    if (lhook) call dr_hook('radiation_flux:out_of_physical_bounds',1,hook_handle)

  end function out_of_physical_bounds
  

  !---------------------------------------------------------------------
  ! Sum elements of "source" into "dest" according to index "ind".
  ! "source" and "ind" should have the same size and bounds, and no
  ! element of "ind" should refer outside the bounds of "dest".  This
  ! version increments existing contents of "dest".
  pure subroutine add_indexed_sum(source, ind, dest)

    real(jprb), intent(in)    :: source(:)
    integer,    intent(in)    :: ind(:)
    real(jprb), intent(inout) :: dest(:)

    integer :: ig, jg, istart, iend

    istart = lbound(source,1)
    iend   = ubound(source,1)

    do jg = istart, iend
      ig = ind(jg)
      dest(ig) = dest(ig) + source(jg)
    end do

  end subroutine add_indexed_sum


  !---------------------------------------------------------------------
  ! As "add_indexed_sum" but this version overwrites existing contents
  ! of "dest"
  pure subroutine indexed_sum(source, ind, dest)

    real(jprb), intent(in)  :: source(:)
    integer,    intent(in)  :: ind(:)
    real(jprb), intent(out) :: dest(:)

    integer :: ig, jg, istart, iend

    dest = 0.0

    istart = lbound(source,1)
    iend   = ubound(source,1)

    do jg = istart, iend
      ig = ind(jg)
      dest(ig) = dest(ig) + source(jg)
    end do

  end subroutine indexed_sum


  !---------------------------------------------------------------------
  ! As "add_indexed_sum" but a whole vertical profiles
  pure subroutine add_indexed_sum_profile(source, ind, dest)

    real(jprb), intent(in)  :: source(:,:)
    integer,    intent(in)  :: ind(:)
    real(jprb), intent(out) :: dest(:,:)

    integer :: ig, jg, istart, iend, jlev, nlev

    istart = lbound(source,1)
    iend   = ubound(source,1)
    nlev   = size(source,2)

    do jlev = 1,nlev
      do jg = istart, iend
        ig = ind(jg)
        dest(ig,jlev) = dest(ig,jlev) + source(jg,jlev)
      end do
    end do

  end subroutine add_indexed_sum_profile


  !---------------------------------------------------------------------
  ! As "indexed_sum" but a whole vertical profiles
  pure subroutine indexed_sum_profile(source, ind, dest)

    real(jprb), intent(in)  :: source(:,:)
    integer,    intent(in)  :: ind(:)
    real(jprb), intent(out) :: dest(:,:)

    integer :: ig, jg, istart, iend, jlev, nlev

    dest = 0.0

    istart = lbound(source,1)
    iend   = ubound(source,1)
    nlev   = size(source,2)

    do jlev = 1,nlev
      do jg = istart, iend
        ig = ind(jg)
        dest(ig,jlev) = dest(ig,jlev) + source(jg,jlev)
      end do
    end do

  end subroutine indexed_sum_profile
  
end module radiation_flux
