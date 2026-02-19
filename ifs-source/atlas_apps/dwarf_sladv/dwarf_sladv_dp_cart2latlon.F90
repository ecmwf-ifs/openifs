! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

module dwarf_sladv_dp_cart2latlon_module
public
contains
subroutine dp_cart2latlon(dt,geometry,uvw_ext,nabla,lonmidp,latmidp,lond,latd)
  !---------------------------------------------------------------------
  ! Aim: Compute d.p. in cartesian space and convert to lat/lon space.
  !      
  ! Based on McGregor's 1993 and Ritchie's 1987 MWR scheme.
  !
  ! Author:   M. Diamantakis June 2018
  ! Modified by G. Tumolo for bit reproducibility in September 2018
  !---------------------------------------------------------------------
  use dwarf_sladv_use_module, only: atlas_Field, atlas_real, atlas_Trace, wp, ip,    &
  &                                 atlas_Nabla, atlas_Functionspace, Geometry_type, &
  &                                 rpi
  implicit none
  real(wp), intent(in)           :: dt
  type(Geometry_type), intent(in)      :: geometry
  real(wp), intent(in)           :: uvw_ext(:,:,:)
  type(atlas_nabla), intent(inout) :: nabla
  real(wp), intent(out)          :: lond(:,:)
  real(wp), intent(out)          :: latd(:,:)
  real(wp), intent(out)          :: lonmidp(:,:) ! SL traj midp lon
  real(wp), intent(out)          :: latmidp(:,:) ! SL traj midp lat
  !---------------------------------------------------------------------
  integer(ip) :: jlev, jnode
  real(wp) :: lon2, pi2, dt2
  real(wp) :: lon_dist
  real(wp) :: X_dp, Y_dp, Z_dp, r_dp
  real(wp) :: sinlo, sinla, coslo, cosla
  type(atlas_Field) :: ufield, vfield, wfield
  type(atlas_Field) :: gradfieldu, gradfieldv, gradfieldw
  real(wp), pointer :: gglon(:), gglat(:)
  real(wp) :: Xcart,  Ycart,  Zcart
  real(wp), pointer :: Ufld(:,:), Vfld(:,:), Wfld(:,:) 
  real(wp), pointer :: gradu(:,:,:), gradv(:,:,:), gradw(:,:,:)
  real(wp) :: zdiff, max_dist

type(atlas_functionspace) :: function_space


type(atlas_Trace) :: trace, trace_convert, trace_nabla, trace_extrap

  trace = atlas_Trace(__FILENAME__,__LINE__,"dp_cart2latlon")
!---------------------------------------------------------------------
! initialize
!---------------------------------------------------------------------

  dt2=0.5_wp*dt*dt

  function_space = nabla%functionspace()

  ufield = function_space%create_field(name="ufld",kind=atlas_real(wp),levels=geometry%nlev)
  vfield = function_space%create_field(name="vfld",kind=atlas_real(wp),levels=geometry%nlev)
  wfield = function_space%create_field(name="wfld",kind=atlas_real(wp),levels=geometry%nlev)

  gradfieldu = function_space%create_field(name="gradu",kind=atlas_real(wp),levels=geometry%nlev,variables=2)
  gradfieldv = function_space%create_field(name="gradv",kind=atlas_real(wp),levels=geometry%nlev,variables=2)
  gradfieldw = function_space%create_field(name="gradw",kind=atlas_real(wp),levels=geometry%nlev,variables=2)

  call function_space%final()

  call ufield%data(Ufld)
  call vfield%data(Vfld)
  call wfield%data(Wfld)
  call gradfieldu%data(gradu)
  call gradfieldv%data(gradv)
  call gradfieldw%data(gradw)

  call geometry%glon%data(gglon)
  call geometry%glat%data(gglat)

  pi2=2.0_wp*rpi

  trace_convert = atlas_Trace(__FILENAME__,__LINE__,"convert_uvw_to_cart")
  !$OMP PARALLEL DO SCHEDULE(STATIC) &
  !$OMP& PRIVATE(jnode,jlev,sinlo,coslo,sinla,cosla)
  do jnode=1,geometry%ngp
    sinlo=sin(gglon(jnode))
    coslo=cos(gglon(jnode))
    sinla=sin(gglat(jnode))
    cosla=cos(gglat(jnode))
    do jlev=1,geometry%nlev
      Ufld(jlev,jnode) = -uvw_ext(1,jlev,jnode)*sinlo - uvw_ext(2,jlev,jnode)*coslo*sinla
      Vfld(jlev,jnode) =  uvw_ext(1,jlev,jnode)*coslo - uvw_ext(2,jlev,jnode)*sinlo*sinla
      Wfld(jlev,jnode) =  uvw_ext(2,jlev,jnode)*cosla
    enddo
  enddo
  !$OMP END PARALLEL DO
  call trace_convert%final()

  ! Halo exchange needed before gradients can be computed
  call ufield%halo_exchange()
  call vfield%halo_exchange()
  call wfield%halo_exchange()

  ! Compute gradients
  trace_nabla = atlas_Trace( __FILENAME__, __LINE__, "nabla%gradient" )
  call nabla%gradient(ufield,gradfieldu)
  call nabla%gradient(vfield,gradfieldv)
  call nabla%gradient(wfield,gradfieldw)
  call trace_nabla%final()

  ! These following halo exchanges are probably not needed
  call gradfieldu%halo_exchange()
  call gradfieldv%halo_exchange()
  call gradfieldw%halo_exchange()
  
  max_dist = 250.0_wp * dt ! a maximum velocity speed of 250.0 m/s is assumed

  trace_extrap = atlas_Trace( __FILENAME__, __LINE__, "extrapolate_dp" )
  
  !$OMP PARALLEL DO SCHEDULE(STATIC) &
  !$OMP& PRIVATE(jnode,jlev,X_dp,Y_dp,Z_dp,r_dp,lon2,Xcart,Ycart,Zcart,sinlo,coslo,sinla,cosla,lon_dist,zdiff)
  do jnode=1,geometry%ngp
    sinlo=sin(gglon(jnode))
    coslo=cos(gglon(jnode))
    sinla=sin(gglat(jnode))
    cosla=cos(gglat(jnode))

    Xcart=geometry%radius*cosla*coslo
    Ycart=geometry%radius*cosla*sinlo
    Zcart=geometry%radius*sinla

    do jlev=1,geometry%nlev

      X_dp = Xcart - dt*Ufld(jlev,jnode) + &
        &         dt2*(uvw_ext(1,jlev,jnode)*gradu(1,jlev,jnode) +    &
        &              uvw_ext(2,jlev,jnode)*gradu(2,jlev,jnode) )

      Y_dp = Ycart - dt*Vfld(jlev,jnode) + &
        &         dt2*(uvw_ext(1,jlev,jnode)*gradv(1,jlev,jnode) +    &
        &              uvw_ext(2,jlev,jnode)*gradv(2,jlev,jnode) )

      Z_dp = Zcart - dt*Wfld(jlev,jnode) + &
        &         dt2*(uvw_ext(1,jlev,jnode)*gradw(1,jlev,jnode) +    &
        &              uvw_ext(2,jlev,jnode)*gradw(2,jlev,jnode) )

      r_dp = sqrt(X_dp*X_dp+Y_dp*Y_dp+Z_dp*Z_dp)

      latd(jlev,jnode)=asin(Z_dp/r_dp)

      latmidp(jlev,jnode)=0.5_wp*(latd(jlev,jnode)+gglat(jnode))

      lon2 = atan2(Y_dp,X_dp)
      ! lond(jlev,jnode)=(0.5_wp-sign(0.5_wp,lon2))*pi2+lon2
      lond(jlev,jnode)=mod((0.5_wp-sign(0.5_wp,lon2))*pi2+lon2,pi2)
      !----------------------------------------------------------------------------------
      ! Assuming that we use a sufficiently wide halo, when the long distance between dp 
      ! and arrival pt exceeds 2*pi-length_of_halo then we adjust the longitude 
      ! to fall in the halo domain of the current mpi task. 
      ! This is a temp fix which will be addressed in future.
      !----------------------------------------------------------------------------------
      lon_dist=lond(jlev,jnode)-gglon(jnode)

      zdiff = geometry%radius*cos(latmidp(jlev,jnode))*abs(lon_dist)-max_dist

      lond(jlev,jnode)=lond(jlev,jnode)-(1.0_wp+sign(1.0_wp,zdiff))*sign(1.0_wp,lon_dist)*rpi
      
      lonmidp(jlev,jnode)=0.5_wp*(lond(jlev,jnode)+gglon(jnode))

    enddo
  enddo
  !$OMP END PARALLEL DO
  call trace_extrap%final()

  call ufield%final()
  call vfield%final()
  call wfield%final()
  call gradfieldu%final()
  call gradfieldv%final()
  call gradfieldw%final()
  call trace%final()

end subroutine

end module
