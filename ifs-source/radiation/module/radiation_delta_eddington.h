! radiation_delta_eddington.h - Delta-Eddington scaling
!
! Copyright (C) 2015-2016 ECMWF
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
! License: see the COPYING file for details
!
! This file is intended to be included inside a module to ensure that
! these simple functions may be inlined

!---------------------------------------------------------------------
! Perform in-place delta-Eddington scaling of the phase function
elemental subroutine delta_eddington(od, ssa, g)

  use parkind1, only : jprb
  
  ! Total optical depth, single scattering albedo and asymmetry
  ! factor
  real(jprb), intent(inout) :: od, ssa, g
  
  ! Fraction of the phase function deemed to be in the forward lobe
  ! and therefore treated as if it is not scattered at all
  real(jprb) :: f
  
  f   = g*g
  od  = od * (1.0_jprb - ssa*f)
  ssa = ssa * (1.0_jprb - f) / (1.0_jprb - ssa*f)
  g   = g / (1.0_jprb + g)
  
end subroutine delta_eddington

subroutine delta_eddington_simd(od, ssa, g)
  use parkind1, only : jprb
#if defined(__PGI) || defined(__NEC__)
  !$omp declare simd(delta_eddington_simd) notinbranch
#else
#ifndef __GFORTRAN__
  !$omp declare simd(delta_eddington_simd) &
  !$omp         linear(ref(od)) linear(ref(ssa)) linear(ref(g)) notinbranch
#endif
#endif
  ! Total optical depth, single scattering albedo and asymmetry
  ! factor
  real(jprb), intent(inout) :: od, ssa, g
  
  ! Fraction of the phase function deemed to be in the forward lobe
  ! and therefore treated as if it is not scattered at all
  real(jprb) :: f
  
  f   = g*g
  od  = od * (1.0_jprb - ssa*f)
  ssa = ssa * (1.0_jprb - f) / (1.0_jprb - ssa*f)
  g   = g / (1.0_jprb + g)
  
end subroutine delta_eddington_simd

!---------------------------------------------------------------------
! Perform in-place delta-Eddington scaling of the phase function, but
! using extensive variables (i.e. the scattering optical depth,
! scat_od, rather than the single-scattering albedo, and the
! scattering-optical-depth-multiplied-by-asymmetry-factor, scat_od_g,
! rather than the asymmetry factor.
elemental subroutine delta_eddington_extensive(od, scat_od, scat_od_g)

  use parkind1, only : jprb

  ! Total optical depth, scattering optical depth and asymmetry factor
  ! multiplied by the scattering optical depth
  real(jprb), intent(inout) :: od, scat_od, scat_od_g

  ! Fraction of the phase function deemed to be in the forward lobe
  ! and therefore treated as if it is not scattered at all
  real(jprb) :: f, g

  if (scat_od > 0.0_jprb) then
    g = scat_od_g / scat_od
  else
    g = 0.0
  end if

  f         = g*g
  od        = od - scat_od * f
  scat_od   = scat_od * (1.0_jprb - f)
  scat_od_g = scat_od * g / (1.0_jprb + g)
  
end subroutine delta_eddington_extensive


!---------------------------------------------------------------------
! Perform in-place delta-Eddington scaling of the phase function,
! using the scattering optical depth rather than the single scattering
! albedo
elemental subroutine delta_eddington_scat_od(od, scat_od, g)

  use parkind1, only : jprb
  
  ! Total optical depth, scattering optical depth and asymmetry factor
  real(jprb), intent(inout) :: od, scat_od, g

  ! Fraction of the phase function deemed to be in the forward lobe
  ! and therefore treated as if it is not scattered at all
  real(jprb) :: f

  f       = g*g
  od      = od - scat_od * f
  scat_od = scat_od * (1.0_jprb - f)
  g       = g / (1.0_jprb + g)

end subroutine delta_eddington_scat_od

