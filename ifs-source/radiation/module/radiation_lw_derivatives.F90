! radiation_lw_derivatives.F90 - Compute longwave derivatives for Hogan and Bozzo (2015) method
!
! (C) Copyright 2016- ECMWF.
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
! This module provides routines to compute the rate of change of
! broadband upwelling longwave flux at each half level with respect to
! the surface broadband upwelling flux.  This is done from the surface
! spectral fluxes and the spectral transmittance of each atmospheric
! layer, assuming no longwave scattering. The result may be used to
! perform approximate updates to the longwave flux profile in between
! calls to the full radiation scheme, accounting for the change in
! skin temperature, following the method of Hogan and Bozzo (JAMES
! 2015).  Separate routines are provided for each solver.
!
! Note that currently a more approximate calculation is performed from
! the exact one in Hogan and Bozzo (2015); here we assume that a
! change in temperature increases the spectral fluxes in proportion,
! when in reality there is a change in shape of the Planck function in
! addition to an overall increase in the total emission.
!
! Modifications
!   2017-10-23  R. Hogan  Renamed single-character variables
!   2022-11-22  P. Ukkonen / R. Hogan  Optimized calc_lw_derivatives_region

module radiation_lw_derivatives

  public

contains

  !---------------------------------------------------------------------
  ! Calculation for the Independent Column Approximation
  subroutine calc_lw_derivatives_ica(ng, nlev, icol, transmittance, flux_up_surf, lw_derivatives)

    use parkind1, only           : jprb
    use yomhook,  only           : lhook, dr_hook, jphook

    implicit none

    ! Inputs
    integer,    intent(in) :: ng   ! number of spectral intervals
    integer,    intent(in) :: nlev ! number of levels
    integer,    intent(in) :: icol ! Index of column for output
    real(jprb), intent(in) :: transmittance(ng,nlev)
    real(jprb), intent(in) :: flux_up_surf(ng) ! Upwelling surface spectral flux (W m-2)
    
    ! Output
    real(jprb), intent(out) :: lw_derivatives(:,:) ! dimensioned (ncol,nlev+1)

    ! Rate of change of spectral flux at a given height with respect
    ! to the surface value
    real(jprb) :: lw_derivatives_g(ng)

    integer    :: jlev

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_lw_derivatives:calc_lw_derivatives_ica',0,hook_handle)

    ! Initialize the derivatives at the surface
    lw_derivatives_g = flux_up_surf / sum(flux_up_surf)
    lw_derivatives(icol, nlev+1) = 1.0_jprb

    ! Move up through the atmosphere computing the derivatives at each
    ! half-level
    do jlev = nlev,1,-1
      lw_derivatives_g = lw_derivatives_g * transmittance(:,jlev)
      lw_derivatives(icol,jlev) = sum(lw_derivatives_g)
    end do

    if (lhook) call dr_hook('radiation_lw_derivatives:calc_lw_derivatives_ica',1,hook_handle)

  end subroutine calc_lw_derivatives_ica


  !---------------------------------------------------------------------
  ! Calculation for the Independent Column Approximation
  subroutine modify_lw_derivatives_ica(ng, nlev, icol, transmittance, &
       &                               flux_up_surf, weight, lw_derivatives)

    use parkind1, only           : jprb
    use yomhook,  only           : lhook, dr_hook, jphook

    implicit none

    ! Inputs
    integer,    intent(in) :: ng   ! number of spectral intervals
    integer,    intent(in) :: nlev ! number of levels
    integer,    intent(in) :: icol ! Index of column for output
    real(jprb), intent(in) :: transmittance(ng,nlev)
    real(jprb), intent(in) :: flux_up_surf(ng) ! Upwelling surface spectral flux (W m-2)
    real(jprb), intent(in) :: weight ! Weight new values against existing
    
    ! Output
    real(jprb), intent(inout) :: lw_derivatives(:,:) ! dimensioned (ncol,nlev+1)

    ! Rate of change of spectral flux at a given height with respect
    ! to the surface value
    real(jprb) :: lw_derivatives_g(ng)

    integer    :: jlev

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_lw_derivatives:modify_lw_derivatives_ica',0,hook_handle)

    ! Initialize the derivatives at the surface
    lw_derivatives_g = flux_up_surf / sum(flux_up_surf)
    ! This value must be 1 so no weighting applied
    lw_derivatives(icol, nlev+1) = 1.0_jprb

    ! Move up through the atmosphere computing the derivatives at each
    ! half-level
    do jlev = nlev,1,-1
      lw_derivatives_g = lw_derivatives_g * transmittance(:,jlev)
      lw_derivatives(icol,jlev) = (1.0_jprb - weight) * lw_derivatives(icol,jlev) &
           &                    + weight * sum(lw_derivatives_g)
    end do

    if (lhook) call dr_hook('radiation_lw_derivatives:modify_lw_derivatives_ica',1,hook_handle)

  end subroutine modify_lw_derivatives_ica



  !---------------------------------------------------------------------
  ! Calculation for solvers involving multiple regions and matrices
  subroutine calc_lw_derivatives_matrix(ng, nlev, nreg, icol, transmittance, &
       &                                u_matrix, flux_up_surf, lw_derivatives)

    use parkind1, only           : jprb
    use yomhook,  only           : lhook, dr_hook, jphook

    use radiation_matrix

    implicit none

    ! Inputs
    integer,    intent(in) :: ng   ! number of spectral intervals
    integer,    intent(in) :: nlev ! number of levels
    integer,    intent(in) :: nreg ! number of regions
    integer,    intent(in) :: icol ! Index of column for output
    real(jprb), intent(in) :: transmittance(ng,nreg,nreg,nlev)
    real(jprb), intent(in) :: u_matrix(nreg,nreg,nlev+1) ! Upward overlap matrix
    real(jprb), intent(in) :: flux_up_surf(ng) ! Upwelling surface spectral flux (W m-2)
    
    ! Output
    real(jprb), intent(out) :: lw_derivatives(:,:) ! dimensioned (ncol,nlev+1)

    ! Rate of change of spectral flux at a given height with respect
    ! to the surface value
    real(jprb) :: lw_derivatives_g_reg(ng,nreg)

    integer    :: jlev

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_lw_derivatives:calc_lw_derivatives_matrix',0,hook_handle)

    ! Initialize the derivatives at the surface; the surface is
    ! treated as a single clear-sky layer so we only need to put
    ! values in region 1.
    lw_derivatives_g_reg = 0.0_jprb
    lw_derivatives_g_reg(:,1) = flux_up_surf / sum(flux_up_surf)
    lw_derivatives(icol, nlev+1) = 1.0_jprb

    ! Move up through the atmosphere computing the derivatives at each
    ! half-level
    do jlev = nlev,1,-1
      ! Compute effect of overlap at half-level jlev+1, yielding
      ! derivatives just above that half-level
      lw_derivatives_g_reg = singlemat_x_vec(ng,ng,nreg,u_matrix(:,:,jlev+1),lw_derivatives_g_reg)

      ! Compute effect of transmittance of layer jlev, yielding
      ! derivatives just below the half-level above (jlev)
      lw_derivatives_g_reg = mat_x_vec(ng,ng,nreg,transmittance(:,:,:,jlev),lw_derivatives_g_reg)

      lw_derivatives(icol, jlev) = sum(lw_derivatives_g_reg)
    end do

    if (lhook) call dr_hook('radiation_lw_derivatives:calc_lw_derivatives_matrix',1,hook_handle)

  end subroutine calc_lw_derivatives_matrix


  !---------------------------------------------------------------------
  ! Calculation for solvers involving multiple regions but no 3D
  ! effects: the difference from calc_lw_derivatives_matrix is that transmittance
  ! has one fewer dimensions
  subroutine calc_lw_derivatives_region(ng, nlev, nreg, icol, transmittance, &
       &                                u_matrix, flux_up_surf, lw_derivatives)

    use parkind1, only           : jprb
    use yomhook,  only           : lhook, dr_hook, jphook

    use radiation_matrix

    implicit none

    ! Inputs
    integer,    intent(in) :: ng   ! number of spectral intervals
    integer,    intent(in) :: nlev ! number of levels
    integer,    intent(in) :: nreg ! number of regions
    integer,    intent(in) :: icol ! Index of column for output
    real(jprb), intent(in) :: transmittance(ng,nreg,nlev)
    real(jprb), intent(in) :: u_matrix(nreg,nreg,nlev+1) ! Upward overlap matrix
    real(jprb), intent(in) :: flux_up_surf(ng) ! Upwelling surface spectral flux (W m-2)
    
    ! Output
    real(jprb), intent(out) :: lw_derivatives(:,:) ! dimensioned (ncol,nlev+1)

    ! Rate of change of spectral flux at a given height with respect
    ! to the surface value
    real(jprb) :: lw_deriv(ng,nreg), lw_deriv_below(ng,nreg)
    real(jprb) :: partial_sum(ng)

    integer    :: jlev, jg

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_lw_derivatives:calc_lw_derivatives_region',0,hook_handle)

    ! Initialize the derivatives at the surface; the surface is
    ! treated as a single clear-sky layer so we only need to put
    ! values in region 1.
    lw_deriv = 0.0_jprb
    lw_deriv(:,1) = flux_up_surf / sum(flux_up_surf)
    lw_derivatives(icol, nlev+1) = 1.0_jprb

    if (nreg == 3) then 
      ! Optimize the most common case of 3 regions by removing the
      ! nested call to singlemat_x_vec and unrolling the matrix
      ! multiplication inline
      
      do jlev = nlev,1,-1
        ! Compute effect of overlap at half-level jlev+1, yielding
        ! derivatives just above that half-level
        lw_deriv_below = lw_deriv
        
        associate(A=>u_matrix(:,:,jlev+1), b=>lw_deriv_below)
          do jg = 1,ng   
            ! Both inner and outer loop of the matrix loops j1 and j2 unrolled
            ! inner loop:        j2=1             j2=2             j2=3 
            lw_deriv(jg,1) = A(1,1)*b(jg,1) + A(1,2)*b(jg,2) + A(1,3)*b(jg,3) 
            lw_deriv(jg,2) = A(2,1)*b(jg,1) + A(2,2)*b(jg,2) + A(2,3)*b(jg,3) 
            lw_deriv(jg,3) = A(3,1)*b(jg,1) + A(3,2)*b(jg,2) + A(3,3)*b(jg,3) 

            ! Compute effect of transmittance of layer jlev, yielding
            ! derivatives just below the half-level above (jlev)
            lw_deriv(jg,1) = lw_deriv(jg,1) * transmittance(jg,1,jlev)
            lw_deriv(jg,2) = lw_deriv(jg,2) * transmittance(jg,2,jlev)
            lw_deriv(jg,3) = lw_deriv(jg,3) * transmittance(jg,3,jlev)

            partial_sum(jg) = lw_deriv(jg,1) + lw_deriv(jg,2) + lw_deriv(jg,3)
          end do
        end associate

        lw_derivatives(icol, jlev) = sum(partial_sum)
      end do
    else
      ! General case when number of regions is not 3
      
      ! Move up through the atmosphere computing the derivatives at each
      ! half-level
      do jlev = nlev,1,-1
        ! Compute effect of overlap at half-level jlev+1, yielding
        ! derivatives just above that half-level
        lw_deriv = singlemat_x_vec(ng,ng,nreg,u_matrix(:,:,jlev+1),lw_deriv)
        
        ! Compute effect of transmittance of layer jlev, yielding
        ! derivatives just below the half-level above (jlev)
        lw_deriv = transmittance(:,:,jlev) * lw_deriv
        
        lw_derivatives(icol, jlev) = sum(lw_deriv)
      end do
    end if
    
    if (lhook) call dr_hook('radiation_lw_derivatives:calc_lw_derivatives_region',1,hook_handle)

  end subroutine calc_lw_derivatives_region


end module radiation_lw_derivatives
