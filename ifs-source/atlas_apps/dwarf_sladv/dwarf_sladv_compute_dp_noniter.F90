! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

module dwarf_sladv_compute_dp_noniter_module
public
contains
subroutine compute_dp_noniter(dt,geometry,wind_t0,wind_ext,nabla,dpmesh)
  use atlas_module
  use dwarf_sladv_dp_cart2latlon_module
  !---------------------------------------------------------------------
  ! Aim: Compute the departure point of the semi-Lagrangian trajectory.
  !
  ! This version is based on McGregor's 1993 MWR scheme.
  ! NOTE: gradient computation needs optimization as it is slow
  !
  ! M. Diamantakis June 2018
  !---------------------------------------------------------------------
  use dwarf_sladv_kind_module
  use dwarf_sladv_geometry_module
  implicit none
  real(wp),                              intent(in) :: dt
  type(Geometry_type),                   intent(in) :: geometry
  type(atlas_Field),                     intent(in) :: wind_t0, wind_ext
  type(atlas_nabla),                     intent(inout) :: nabla
  type(atlas_FieldSet),                  intent(inout) :: dpmesh
  !---------------------------------------------------------------------
  integer(ip) :: npts, i, jlev, nlev, ngp, jnode
  real(wp) :: dt2, z_tmp, dwdz
  type(atlas_Field) :: field_lon_dp, field_lat_dp, field_z_dp, field_etadot_midp
  real(wp) :: lon_tmp(geometry%nlev,geometry%ngp)
  real(wp) :: lat_tmp(geometry%nlev,geometry%ngp)
  real(wp), pointer :: lon_dp(:,:), lat_dp(:,:), z_dp(:,:), etadot_midp(:,:)
  real(wp), pointer :: zco(:)
  real(wp), pointer :: uvw_t0(:,:,:), uvw_ext(:,:,:)
  real(wp) :: vcmin, vcmax
  type(atlas_Trace) :: trace
  !---------------------------------------------------------------------
  ! initialize
  !---------------------------------------------------------------------
  trace = atlas_Trace( __FILENAME__, __LINE__, "compute_dp_noniter" )
  npts=geometry%ngptot
  ngp=geometry%ngp
  nlev=geometry%nlev
  dt2=0.5_wp*dt*dt

  field_etadot_midp = atlas_Field(name="etadot_midp",kind=atlas_real(wp),shape=[nlev,npts])

  field_lon_dp=dpmesh%field("dplonfield")
  field_lat_dp=dpmesh%field("dplatfield")
  field_z_dp=dpmesh%field("dpvertfield")

  call field_lon_dp%data(lon_dp)
  call field_lat_dp%data(lat_dp)
  call field_z_dp%data(z_dp)
  call geometry%zcoord%data(zco)
  call wind_ext%data(uvw_ext)
  call field_etadot_midp%data(etadot_midp)

  call wind_ext%halo_exchange()

  vcmin=zco(1)
  vcmax=zco(nlev)

  !----------------------------------------------------------------------
  ! Compute horizontal components of the d.p.
  !----------------------------------------------------------------------
  call dp_cart2latlon(dt,geometry,uvw_ext,nabla,lon_dp,lat_dp,lon_tmp,lat_tmp)

  ! w-component interpol

  call interpolate_atlas(uvw_ext,etadot_midp)

!------------------------------------------------------------------------------------------------
! add 2nd order term for the vertical separately
!------------------------------------------------------------------------------------------------
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,jlev,z_tmp,dwdz)
  do jnode=1,ngp

    do jlev=1, nlev
      lon_dp(jlev,jnode) = lon_tmp(jlev,jnode)! restore dp lon from temp buffer
      lat_dp(jlev,jnode) = lat_tmp(jlev,jnode)! restore dp lat from temp buffer
      z_dp(jlev,jnode)   = zco(jlev) - dt*etadot_midp(jlev,jnode)
    enddo

    do jlev=2, nlev-1
      dwdz           = (etadot_midp(jlev+1,jnode)-etadot_midp(jlev-1,jnode))/(zco(jlev+1)-zco(jlev-1))
      z_tmp          = z_dp(jlev,jnode) + dt2*etadot_midp(jlev,jnode)*dwdz
      z_dp(jlev,jnode) = min(vcmax,max(vcmin,z_tmp))
    enddo

    ! level 1
    dwdz           = (etadot_midp(2,jnode)-etadot_midp(1,jnode))/(zco(2)-zco(1))
    z_tmp          = z_dp(1,jnode) + dt2*etadot_midp(1,jnode)*dwdz
    z_dp(1,jnode)    = min(vcmax,max(vcmin,z_tmp))

    ! level nlev
    dwdz           = (etadot_midp(nlev,jnode)-etadot_midp(nlev-1,jnode))/(zco(nlev)-zco(nlev-1))
    z_tmp          = z_dp(nlev,jnode) + dt2*etadot_midp(nlev,jnode)*dwdz
    z_dp(nlev,jnode) = min(vcmax,max(vcmin,z_tmp))
  enddo
!$OMP END PARALLEL DO

  !-----------------------
  ! Finalize
  !-----------------------
  call field_etadot_midp%final()
  call field_z_dp%final()
  call field_lat_dp%final()
  call field_lon_dp%final()

  call trace%final()

contains

  subroutine interpolate_atlas(uvw_ext,edot_midp)
    real(wp), intent(in)    :: uvw_ext(:,:,:)
    real(wp), intent(inout) :: edot_midp(:,:)
    !--------------------------------------------------------
    type(atlas_interpolation) :: interpolation
    type(atlas_Config) :: config
    type(atlas_Field) :: wind_vert

    type(atlas_Field) :: londp_field, latdp_field
    real(wp), pointer :: londp(:,:), latdp(:,:)

    type(atlas_Field) :: field

    real(wp), pointer :: dplonlat(:,:)

    type(atlas_Field) :: etadot_midp_slice

    type(atlas_FieldSet) :: dplonlat_slice
    integer :: jlev

    londp_field = dpmesh%field("dplonfield")
    latdp_field = dpmesh%field("dplatfield")
    call londp_field%data( londp )
    call latdp_field%data( latdp )

    config = atlas_Config()
    call config%set("type","linear2D")
    call config%set("matrix_free",.true.)

    ! A lot of 2D horizontal interpolations (1 for each level).
    ! Each level has different departure points, so interpolators cannot be reused.
    do jlev=1,geometry%nlev
      dplonlat_slice = atlas_FieldSet()
      call dplonlat_slice%add( atlas_Field(name="dplon_level", data=londp(jlev,1:ngp) ) )
      call dplonlat_slice%add( atlas_Field(name="dplat_level", data=latdp(jlev,1:ngp) ) )
      field = dplonlat_slice%field("dplon_level")
      call set_metadata_radians( field )
      field = dplonlat_slice%field("dplat_level")
      call set_metadata_radians( field )
      wind_vert = atlas_Field( "wind_vert_level", uvw_ext(3,jlev,:) )
      etadot_midp_slice = atlas_Field("etadot_level", edot_midp(jlev,1:ngp) )
      interpolation = atlas_interpolation( config, geometry%fs_structuredcolumns, dplonlat_slice )
      call interpolation%execute( wind_vert, etadot_midp_slice )
      call interpolation%final()
      call wind_vert%final()
      call dplonlat_slice%final()
      call etadot_midp_slice%final()
    enddo

    call config%final()
    call londp_field%final()
    call latdp_field%final()
  end subroutine

  subroutine set_metadata_radians( field )
    type(atlas_Field) :: field
    type(atlas_Metadata) :: metadata
    metadata = field%metadata()
    call metadata%set("units","radians")
    call metadata%final()
  end subroutine

end subroutine
end module
