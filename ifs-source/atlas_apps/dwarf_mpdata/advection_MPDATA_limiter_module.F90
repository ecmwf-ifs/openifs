! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

module advection_MPDATA_limiter_module

use fckit_module, only : fckit_log
use atlas_module
use advection_MPDATA_auxiliary_module
use advection_MPDATA_geometry_module, only : geometry_type

implicit none

public :: limiter_type

private

! Class Limiter
type :: limiter_type
  type(geometry_type), pointer :: geom
  real(wp) :: eps
  real(wp) :: mivbz

contains
  generic  , public  :: assignment(=) => copy
  procedure, private :: copy
  procedure, public  :: compute_min_max_xy => compute_scalar_max_and_min_xy
  procedure, public  :: compute_min_max_z  => compute_scalar_max_and_min_z 
  procedure, public  :: limit_flux         => limit_scalar_flux
  
end type

! Definition of the constructor
interface limiter_type
  module procedure limiter_constructor
end interface limiter_type


contains

! =============================================================================
! Constructor for class limiter
! =============================================================================
function limiter_constructor(geometry, config) &
& result(this)
type(limiter_type) :: this
type(geometry_type), target :: geometry
type(atlas_config)          :: config

this%geom   => geometry

if( .not. config%get("eps0" , this%eps   ) ) call fckit_log%error("eps0 not found in config")
if( .not. config%get("ivbz" , this%mivbz ) ) call fckit_log%error("ivbz not found in config")

end function

subroutine copy(obj_out,obj_in)
  class(limiter_type), intent(inout) :: obj_out
  class(limiter_type), intent(in) :: obj_in
  obj_out%geom         => obj_in%geom
  obj_out%eps          =  obj_in%eps
  obj_out%mivbz        =  obj_in%mivbz
end subroutine copy

! =============================================================================

subroutine compute_scalar_max_and_min_xy(this,pD,pDmax,pDmin,kflip)
class(limiter_type), intent(inout) :: this
real(wp), intent(in)  :: pD(:,:)
real(wp), intent(out) :: pDmax(:,:), pDmin(:,:)
integer, intent(in)   :: kflip
integer                               :: nb_levels
integer                               :: nb_nodes
integer                               :: nb_edges
real(wp) :: zD, zpbc
integer :: jnode, jlev, jedge, ip2, iedge
integer, pointer                   :: inode2edge  (:)
integer                            :: inode2edge_size

call fckit_log%debug('compute_scalar_max_and_min_xy')

nb_nodes  = this%geom%nb_nodes
nb_edges  = this%geom%nb_edges
nb_levels = this%geom%nb_levels

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,inode2edge, &
!$OMP & inode2edge_size,iedge,jlev,jedge,ip2,zpbc,zD)
do jnode  = 1,nb_nodes
  do jlev = 1,nb_levels
    pDmax(jlev,jnode) = pD(jlev,jnode)
    pDmin(jlev,jnode) = pD(jlev,jnode)
  enddo
  call this%geom%node2edge%row(jnode,inode2edge,inode2edge_size)
  do jedge = 1,inode2edge_size
    iedge = inode2edge(jedge)
    if ( this%geom%iedge2node(1,iedge) == jnode ) then
      ip2 = this%geom%iedge2node(2,iedge)
    else
      ip2 = this%geom%iedge2node(1,iedge)
    endif    
    zpbc  = (1-kflip)+kflip*this%geom%rpole_bc(iedge)
    do jlev = 1,nb_levels
      zD                = pD(jlev,ip2)
      pDmax(jlev,jnode) = max(pDmax(jlev,jnode),zpbc*zD)
      pDmin(jlev,jnode) = min(pDmin(jlev,jnode),zpbc*zD)
    enddo
  enddo
enddo
!$OMP END PARALLEL DO

end subroutine compute_scalar_max_and_min_xy

!###################################################################################

subroutine compute_scalar_max_and_min_z(this,pD,pDmax,pDmin)
class(limiter_type), intent(inout) :: this
real(wp), intent(in)  :: pD(:,:)
real(wp), intent(out) :: pDmax(:,:), pDmin(:,:)
integer                               :: nb_levels
integer                               :: nb_nodes
integer :: jnode, jlev

call fckit_log%debug('compute_scalar_max_and_min_z')

nb_nodes  = this%geom%nb_nodes
nb_levels = this%geom%nb_levels

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,jlev)
do jnode  = 1,nb_nodes
  do jlev = 2,nb_levels-1
    pDmax(jlev,jnode) = max(pD(jlev-1,jnode),pD(jlev,jnode),pD(jlev+1,jnode),pDmax(jlev,jnode))
    pDmin(jlev,jnode) = min(pD(jlev-1,jnode),pD(jlev,jnode),pD(jlev+1,jnode),pDmin(jlev,jnode))
  enddo
  pDmax(  1  ,jnode)     = max(pD(  1   ,jnode),pD(  2  ,jnode),pDmax(  1  ,jnode))
  pDmin(  1  ,jnode)     = min(pD(  1   ,jnode),pD(  2  ,jnode),pDmin(  1  ,jnode))
  pDmax(nb_levels,jnode) = max(pD(nb_levels-1,jnode),pD(nb_levels,jnode),pDmax(nb_levels,jnode))
  pDmin(nb_levels,jnode) = min(pD(nb_levels-1,jnode),pD(nb_levels,jnode),pDmin(nb_levels,jnode))
enddo
!$OMP END PARALLEL DO
end subroutine compute_scalar_max_and_min_z

!###################################################################################

subroutine limit_scalar_flux(this,pflux,pfluz,pD,pDmax,pDmin,plimit,pdt,pdz,pvol,prho)
class(limiter_type), intent(inout) :: this
real(wp), intent(inout) :: pflux(:,:),pfluz(:,:)
real(wp), intent(in)    :: pD(:,:), pDmax(:,:), pDmin(:,:), prho(:,:)
real(wp), intent(in)    :: pvol(:)
real(wp), intent(in)    :: plimit, pdt, pdz
type(atlas_HaloExchange)              :: halo_exchange
integer                               :: nb_levels
integer                               :: nb_nodes
integer                               :: nb_edges
real(wp) :: zsignp, zsignn, zadd, zpos, zneg
real(wp) :: zrhin(this%geom%nb_levels,this%geom%nb_nodes)
real(wp) :: zrhout(this%geom%nb_levels,this%geom%nb_nodes)
real(wp) :: cp(this%geom%nb_levels,this%geom%nb_nodes)
real(wp) :: cn(this%geom%nb_levels,this%geom%nb_nodes)
real(wp) :: zdzi
integer :: jnode, jedge, iedge, jlev, ip1, ip2
integer, pointer                   :: inode2edge  (:)
integer                            :: inode2edge_size

call fckit_log%debug('limit_scalar_flux')

! Assign required values
nb_nodes  = this%geom%nb_nodes
nb_edges  = this%geom%nb_edges
nb_levels = this%geom%nb_levels

zdzi = 1._wp/pdz

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(inode2edge,inode2edge_size, &
!$OMP & iedge,jedge,ip1,ip2,jlev,jnode,zadd,zsignp,zsignn,zpos,zneg)
do jnode  = 1,nb_nodes
  do jlev = 1,nb_levels
    zrhout(jlev,jnode) = 0._wp
    zrhin(jlev,jnode)  = 0._wp
  enddo
  call this%geom%node2edge%row(jnode,inode2edge,inode2edge_size)
  do jedge = 1,inode2edge_size
    iedge  = inode2edge(jedge)
    ip1    = this%geom%iedge2node(1,iedge)
    ip2    = this%geom%iedge2node(2,iedge)
    zadd   = real(this%geom%node2edge_sign(jedge,jnode),wp)
    zsignp = max(0._wp,zadd)
    zsignn = min(0._wp,zadd)
    do jlev = 1,nb_levels
      zpos               = max(0._wp,pflux(jlev,iedge))
      zneg               = min(0._wp,pflux(jlev,iedge))
      zrhin (jlev,jnode) = zrhin (jlev,jnode)-zsignp*zneg-zsignn*zpos
      zrhout(jlev,jnode) = zrhout(jlev,jnode)+zsignp*zpos+zsignn*zneg
    enddo
  enddo
  do jlev = 1,nb_levels
    zrhin(jlev,jnode)  = zrhin(jlev,jnode)/pvol(jnode)       &
                       & +(max(0._wp,pfluz( jlev ,jnode))-min(0._wp,pfluz(jlev+1,jnode)))*zdzi
    zrhout(jlev,jnode) = zrhout(jlev,jnode)/pvol(jnode)      &
                       & +(max(0._wp,pfluz(jlev+1,jnode))-min(0._wp,pfluz(jlev  ,jnode)))*zdzi
    cp(jlev,jnode)     = (pDmax(jlev,jnode)-pD(jlev,jnode))  &
                       & *prho(jlev,jnode)/(zrhin(jlev,jnode)*pdt+this%eps)
    cn(jlev,jnode)     = (pD(jlev,jnode)-pDmin(jlev,jnode))  &
                       & *prho(jlev,jnode)/( zrhout(jlev,jnode)* pdt + this%eps )
  enddo
enddo
!$OMP END PARALLEL DO

halo_exchange = this%geom%nodes%get_halo_exchange()
call halo_exchange%execute(cp)
call halo_exchange%execute(cn)

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(ip1,ip2,jlev,jedge)
do jedge  = 1,nb_edges
  ip1 = this%geom%iedge2node(1,jedge)
  ip2 = this%geom%iedge2node(2,jedge)
  do jlev = 1,nb_levels
    pflux(jlev,jedge) = max(0._wp,pflux(jlev,jedge))           &
                      & *min(plimit,cp(jlev,ip2),cn(jlev,ip1)) &
                      & +min(0._wp,pflux(jlev,jedge))          &
                      & *min(plimit,cn(jlev,ip2),cp(jlev,ip1))
  enddo
enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,jlev)
do jnode  = 1,nb_nodes
  do jlev = 2,nb_levels
    pfluz(jlev,jnode)      = max(0._wp,pfluz(jlev,jnode))                 &
                           & *min(plimit,cp(jlev,jnode),cn(jlev-1,jnode)) &
                           & +min(0._wp,pfluz(jlev,jnode))                &
                           & *min(plimit,cn(jlev,jnode),cp(jlev-1,jnode))
  enddo
  pfluz(  1        ,jnode) = this%mivbz*pfluz( 2       ,jnode)
  pfluz(nb_levels+1,jnode) = this%mivbz*pfluz(nb_levels,jnode)
enddo
!$OMP END PARALLEL DO

end subroutine limit_scalar_flux

!###################################################################################

end module advection_MPDATA_limiter_module
