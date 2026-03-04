! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

#include "dwarf_sladv.h"

module dwarf_sladv_compute_dp_iter_module

public

contains

subroutine compute_dp_iter(dt,dp_niter,geometry,wind_t0,wind_ext,dpmesh)
  use dwarf_sladv_use_module, only: atlas_real, atlas_Field, atlas_FieldSet, &
  &                        atlas_Trace, atlas_Config, atlas_Interpolation, rpi
  use dwarf_sladv_dp_rotmat_trans_module, only: dp_rotmat_trans
  !---------------------------------------------------------------------
  ! Aim: Compute the departure point of the semi-Lagrangian trajectory.
  !
  ! This version implements the rotation matrix method described
  ! by Temperton et al QJRMS 2001 
  !
  ! M. Diamantakis June 2018
  !---------------------------------------------------------------------
  use dwarf_sladv_kind_module
  use dwarf_sladv_geometry_module, only: Geometry_type
  implicit none
  real(wp),              intent(in)    :: dt
  integer(ip),           intent(in)    :: dp_niter
  type(Geometry_type),   intent(in)    :: geometry
  type(atlas_Field),     intent(in)    :: wind_t0, wind_ext
  type(atlas_FieldSet),  intent(inout) :: dpmesh
  !---------------------------------------------------------------------
  integer(ip) :: npts, jlev, nlev, ngp, iter, jnode
  real(wp) :: zver, dt2
  real(wp) :: px(geometry%nlev,geometry%ngp), qy(geometry%nlev,geometry%ngp)
  type(atlas_Field) :: lon_dp, lat_dp, z_dp
  type(atlas_Field) :: wind_dp
  real(wp), pointer :: gglon(:), gglat(:)
  real(wp), pointer :: lond(:,:), latd(:,:), zd(:,:), zco(:)
  real(wp), pointer :: uvw_t0(:,:,:), uvw_ext(:,:,:), uvw_dp(:,:,:)
  real(wp), pointer :: dp_view(:,:,:)
  real(wp) :: velmax
  real(wp) :: vcmin, vcmax
  type(atlas_Trace) :: trace

  trace = atlas_Trace(__FILENAME__,__LINE__,"compute_dp_iter")

  velmax = 250.0_wp
  !---------------------------------------------------------------------
  ! initialize
  !---------------------------------------------------------------------
  ngp=geometry%ngp
  nlev=geometry%nlev

  dt2=0.5_wp*dt
 
  wind_dp = geometry%fs_structuredcolumns%create_field( &
    name="wind(dp)", kind=atlas_real(wp), variables=3, type="vector" )

  lon_dp=dpmesh%field("dplonfield")
  lat_dp=dpmesh%field("dplatfield")
  z_dp=dpmesh%field("dpvertfield")

  call lon_dp%data(lond)
  call lat_dp%data(latd)
  call z_dp%data(zd)

  call geometry%glon%data(gglon)
  call geometry%glat%data(gglat)
  call geometry%zcoord%data(zco)
  call wind_ext%data(uvw_ext)
  call wind_t0%data(uvw_t0)
  call wind_dp%data(uvw_dp)

  call wind_ext%halo_exchange()

  vcmin=zco(1)
  vcmax=zco(nlev)
  
  !----------------------------------------------------------------------
  ! Compute initial d.p. estimate
  !----------------------------------------------------------------------
  iter=1
  call dp_rotmat_trans(iter,dt,dp_niter,geometry,uvw_t0,uvw_dp,px,qy,lond,latd)
  !! --> will have initialized px, qy

  !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,jlev,zver)
  do jnode=1, ngp
    do jlev=1, nlev
      zver           = zco(jlev) - dt*uvw_t0(3,jlev,jnode)
      zd(jlev,jnode) = min(vcmax,max(vcmin,zver))
    enddo
  enddo
  !$OMP END PARALLEL DO

  do iter=2, dp_niter

    call interpolate_wind_ext_to_dp()
    call dp_rotmat_trans(iter,dt,dp_niter,geometry,uvw_t0,uvw_dp,px,qy,lond,latd)

    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,jlev,zver)
    do jnode=1, ngp
      do jlev=1, nlev
        zver           = zco(jlev) - dt2*(uvw_dp(3,jlev,jnode) + uvw_t0(3,jlev,jnode))
        zd(jlev,jnode) = min(vcmax,max(vcmin,zver))
      enddo
    enddo
    !$OMP END PARALLEL DO

  enddo

  !-----------------------
  ! Finalize
  !-----------------------

  call lon_dp%final()
  call lat_dp%final()
  call z_dp%final()
  call wind_dp%final()
  call trace%final()

contains

  subroutine interpolate_wind_ext_to_dp()
    type(atlas_Config) :: config
    type(atlas_interpolation) :: interpolation
    type(atlas_Trace) :: trace
    trace = atlas_Trace(__FILENAME__,__LINE__,"interpolate_wind_ext_to_dp")

    config = atlas_Config()
    call config%set("type","linear3D")
    call config%set("matrix_free",.true.)
    interpolation = atlas_interpolation( config, geometry%fs_structuredcolumns, dpmesh )
    call interpolation%execute( wind_ext, wind_dp )
    call interpolation%final()

    call trace%final()
  end subroutine
end subroutine


end module
