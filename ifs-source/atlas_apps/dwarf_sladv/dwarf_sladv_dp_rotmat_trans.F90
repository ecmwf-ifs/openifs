! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

#include "dwarf_sladv.h"

module dwarf_sladv_dp_rotmat_trans_module

public

contains

subroutine dp_rotmat_trans(iter,dt,dp_niter,geometry,uvw_t0,uvw_dp,px,qy,lond,latd)
  !---------------------------------------------------------------------
  ! Aim: Compute the departure point lon, lat, height of the 
  !      semi-Lagrangian trajectory applying the method by
  !      Temperton et al QJRMS 2001. 
  !      Performs iteration no "iter" of an iterative scheme (mid-point or SETTLS).
  !
  ! M. Diamantakis June 2018
  ! Modified by G. Tumolo for bit reproducibility in September 2018
  !---------------------------------------------------------------------
  use dwarf_sladv_use_module
  implicit none
  integer(ip),         intent(in)    :: iter
  type(Geometry_type), intent(in)    :: geometry
  real(wp),            intent(in)    :: dt
  integer(ip),         intent(in)    :: dp_niter
  real(wp),            intent(in)    :: uvw_t0(:,:,:)
  real(wp),            intent(in)    :: uvw_dp(:,:,:)
  real(wp),            intent(inout) :: px(:,:), qy(:,:)
  real(wp),            intent(inout) :: lond(:,:), latd(:,:)
  !---------------------------------------------------------------------
  integer(ip) :: jlev, nlev, ngp, jnode
  real(wp) :: ra, pi2, max_dist
  real(wp) :: rdtsa, zdtsa, rdtsa2, rdts62, zdts62, rdts22, zdts22
  real(wp) :: zpu, zpv
  real(wp), pointer :: gglon(:), gglat(:), zco(:)
  real(wp) :: sinlo, sinla, coslo, cosla
  type(atlas_Trace) :: trace
  !---------------------------------------------------------------------
  trace = atlas_Trace(__FILENAME__,__LINE__,"dp_rotmat_trans")
  !---------------------------------------------------------------------
  ! initialize
  !---------------------------------------------------------------------
  ngp=geometry%ngp
  nlev=geometry%nlev
  ra=geometry%radius

  call geometry%glon%data(gglon)
  call geometry%glat%data(gglat)

  rdtsa=0.5_wp*dt/RA
  zdtsa=2.0_wp*rdtsa
  rdtsa2=rdtsa**2
  rdts62=rdtsa2/6.0_wp
  zdts62=4.0_wp*rdts62
  rdts22=rdtsa2/2.0_wp
  zdts22=4.0_wp*rdts22
  pi2=2.0_wp*rpi
  max_dist = 250.0_wp * dt ! a maximum velocity speed of 250.0 m/s is assumed

  !------------------------------------------------------------------------------------
  ! Rotate winds at dp and compute d.p. at current iteration estimate of wind field on the spherical domain
  !------------------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  if (iter==1) then
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,jlev,sinlo,coslo,sinla,cosla,zpu,zpv)
    do jnode=1, ngp
      sinlo=sin(gglon(jnode))
      coslo=cos(gglon(jnode))
      sinla=sin(gglat(jnode))
      cosla=cos(gglat(jnode))
      do jlev=1, nlev
        zpu = uvw_t0(1,jlev,jnode)
        zpv = uvw_t0(2,jlev,jnode) 
        call compute_dp(jnode, jlev, sinlo, coslo, sinla, cosla, zpu, zpv)
      enddo
    enddo
    !$OMP END PARALLEL DO

  else ! SETTLS:
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,jlev,sinlo,coslo,sinla,cosla,zpu,zpv)
    do jnode=1, ngp
      sinlo=sin(gglon(jnode))
      coslo=cos(gglon(jnode))
      sinla=sin(gglat(jnode))
      cosla=cos(gglat(jnode))
      do jlev=1, nlev
        ! form 1/2(Vd+Vi) applying rotation matrix at Vd
        zpu = 0.5_wp*(px(jlev,jnode)*uvw_dp(1,jlev,jnode)+ & 
          &                        qy(jlev,jnode)*uvw_dp(2,jlev,jnode)+ &
          &                                       uvw_t0(1,jlev,jnode))  
        zpv = 0.5_wp*(-qy(jlev,jnode)*uvw_dp(1,jlev,jnode)+ & 
          &                         px(jlev,jnode)*uvw_dp(2,jlev,jnode)+ &
          &                                        uvw_t0(2,jlev,jnode))
        call compute_dp(jnode, jlev, sinlo, coslo, sinla, cosla, zpu, zpv)
      enddo
    enddo
    !$OMP END PARALLEL DO
  endif
  call trace%final()

contains

  subroutine compute_dp(jnode,jlev,sinlo,coslo,sinla,cosla,zpu,zpv)
    !------------------------------------------------------------------------------------
    ! compute d.p. at current iteration estimate of wind field on the spherical domain
    !------------------------------------------------------------------------------------
    integer(ip),intent(in) :: jnode, jlev
    real(wp),intent(in) :: sinlo,coslo,sinla,cosla
    real(wp),intent(in) :: zpu,zpv
    !------------------------------------------
    real(wp) :: znor2,zfac1,zfac2,cosphi,sinlat,coscos,sincos,zcosla,ztmp1,ztmp2
    real(wp) :: zsinla,zcoslo,latmidp,lon_dist,zdiff,zinv1,zinv2
    !------------------------------------------
    znor2=zpu*zpu+zpv*zpv
    zfac1  = zdtsa*(1.0_wp-zdts62*znor2)
    zfac2  = zpv*zfac1
    cosphi = 1.0_wp-zdts22*znor2
    sinlat = sinla*cosphi-zfac2*cosla
    coscos = cosla*cosphi+zfac2*sinla
    sincos = -zpu*zfac1
    ! compute lat/lon from trig products
    zcosla=max(rndoff,sqrt(coscos*coscos+sincos*sincos))
    ztmp1=coslo*coscos-sinlo*sincos
    ztmp2=sinlo*coscos+coslo*sincos
    zsinla=max(-1.0_wp,min(1.0_wp,sinlat))
    zcoslo=max(-1.0_wp,min(1.0_wp,ztmp1/zcosla))
    latd(jlev,jnode)=asin(zsinla)
    latmidp=0.5_wp*(latd(jlev,jnode)+gglat(jnode))
    lond(jlev,jnode)=rpi+sign(1.0_wp,ztmp2)*(acos(zcoslo)-rpi)
    !----------------------------------------------------------------------------------
    ! Assuming that we use a sufficiently wide halo, when the long distance between dp 
    ! and arrival pt exceeds 2*pi-length_of_halo then we adjust the longitude 
    ! to fall in the halo domain of the current mpi task. 
    ! This is a temp fix which will be addressed in future.
    !----------------------------------------------------------------------------------
    lon_dist=lond(jlev,jnode)-gglon(jnode)
    zdiff = ra*cos(latmidp)*abs(lon_dist)-max_dist
    lond(jlev,jnode)=lond(jlev,jnode)-(1.0_wp+sign(1.0_wp,zdiff))*sign(1.0_wp,lon_dist)*rpi
    ! update rotation matrix
    if (iter<dp_niter) then
      zinv1=1.0_wp/(1.0_wp+cosphi)
      zinv2=1.0_wp/zcosla
      px(jlev,jnode)=(cosla*zcosla &
        & +(1.0_wp+sinla*sinlat)*coscos*zinv2)*zinv1
      qy(jlev,jnode)=-(sinla+sinlat)*sincos*zinv2*zinv1
    endif
  end subroutine

end subroutine

end module
