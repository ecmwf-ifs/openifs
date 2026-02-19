! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

#define __FILENAME__ "advection_MPDATA_MPDATA_module.F90"
! =============================================================================
! MPDATA_module
! This module contains strictly algorithmic subroutines
! - MPDATA advection
! =============================================================================

module advection_MPDATA_MPDATA_module

use atlas_module
use fckit_module, only : fckit_log, fckit_exception
use advection_MPDATA_auxiliary_module
use advection_MPDATA_limiter_module
use advection_MPDATA_geometry_module, only : geometry_type

implicit none

private

public :: MPDATA_type

! Class MPDATA
type :: MPDATA_type
  type(geometry_type), pointer :: geom
  type(limiter_type)           :: limiter
  real(wp)                     :: rlimit
  real(wp)                     :: rIVBZ
  integer                      :: mporder
contains
  generic  , public  :: assignment(=) => copy
  procedure, private :: copy
  procedure, private :: MPDATA_scheme
  procedure, private :: MPDATA_density_scheme
  procedure, public :: execute
  procedure, public :: execute_density
  procedure, public :: advectors_in_edges
  procedure, public :: face_normals
end type

! Definition of the constructor
interface MPDATA_type
  module procedure MPDATA_constructor
end interface MPDATA_type

contains

! =============================================================================
! Constructor for class MPDATA
! =============================================================================
function MPDATA_constructor( geometry, config ) result(this)
type(MPDATA_type)           :: this
type(geometry_type), target :: geometry
type(atlas_config)          :: config

if( .not. config%get("limit",   this%rlimit) )  call fckit_exception%abort( "limit not found")
if( .not. config%get("ivbz",    this%rivbz) )   call fckit_exception%abort( "ivbz not found")
if( .not. config%get("mporder", this%mporder) ) call fckit_exception%abort( "mporder not found")

this%geom             => geometry

this%limiter = limiter_type(this%geom, config)
end function MPDATA_constructor
! =============================================================================

subroutine copy(obj_out,obj_in)
  class(MPDATA_type), intent(inout) :: obj_out
  class(MPDATA_type), intent(in) :: obj_in

  obj_out%geom            => obj_in%geom
  obj_out%limiter         =  obj_in%limiter
  
  obj_out%rlimit     = obj_in%rlimit
  obj_out%rivbz      = obj_in%rivbz
  obj_out%mporder    = obj_in%mporder
end subroutine copy

! =============================================================================
! Subroutine execute: driver to call core MPDATA routine
! =============================================================================
subroutine execute(this, dt, itrac, field_sol, field_rho, field_rhofac, field_rhodt, field_Vn, field_Wn)
class(MPDATA_type), intent(inout) :: this

real(wp), intent(in) :: dt

integer, intent(in)                        :: itrac 

type(atlas_Field), intent(inout)           :: field_sol
  !! field to advect

type(atlas_Field), intent(in)              :: field_rho
  !! density, includes geometry

type(atlas_Field), intent(in)              :: field_rhofac
  !! factor to multiply with density in second pass (value=1.0 hardcoded in this dwarf)

type(atlas_Field), intent(in)              :: field_rhodt
  !! For time-dependent geometry (value=0.0 hardcoded in this dwarf)

type(atlas_Field), intent(in)              :: field_Vn
  !! Horizontal components winds located at edges

type(atlas_Field), intent(in)              :: field_Wn
  !! Vertical component of wind located at half-levels, in nodes

real(wp), pointer :: sol   (:,:,:)
real(wp), pointer :: rho   (:,:)
real(wp), pointer :: rhofac(:,:)
real(wp), pointer :: rhodt (:,:)
real(wp), pointer :: Vn    (:,:)
real(wp), pointer :: Wn    (:,:)

call field_sol%data(sol)
call field_rho%data(rho)
call field_rhofac%data(rhofac)
call field_rhodt%data(rhodt)
call field_Vn%data(Vn)
call field_Wn%data(Wn)


call this%MPDATA_scheme(dt,sol(itrac,:,:),rho,rhofac,rhodt,Vn,Wn,kflip=0)

end subroutine execute

! =============================================================================

! =============================================================================
! Subroutine execute: driver to call core MPDATA routine
! =============================================================================
subroutine execute_density(this, dt, field_sol, field_rho, field_rhofac, field_rhodt, field_Vn, field_Wn)
class(MPDATA_type), intent(inout) :: this

real(wp), intent(in) :: dt

type(atlas_Field), intent(inout)           :: field_sol
  !! field to advect

type(atlas_Field), intent(in)              :: field_rho
  !! density, includes geometry

type(atlas_Field), intent(in)              :: field_rhofac
  !! factor to multiply with density in second pass (value=1.0 hardcoded in this dwarf)

type(atlas_Field), intent(in)              :: field_rhodt
  !! For time-dependent geometry (value=0.0 hardcoded in this dwarf)

type(atlas_Field), intent(inout)              :: field_Vn
  !! Horizontal components winds located at edges

type(atlas_Field), intent(inout)              :: field_Wn
  !! Vertical component of wind located at half-levels, in nodes

real(wp), pointer :: sol   (:,:)
real(wp), pointer :: rho   (:,:)
real(wp), pointer :: rhofac(:,:)
real(wp), pointer :: rhodt (:,:)
real(wp), pointer :: Vn    (:,:)
real(wp), pointer :: Wn    (:,:)

call field_sol%data(sol)
call field_rho%data(rho)
call field_rhofac%data(rhofac)
call field_rhodt%data(rhodt)
call field_Vn%data(Vn)
call field_Wn%data(Wn)

call this%MPDATA_density_scheme(dt,sol,rho,rhofac,rhodt,Vn,Wn,0)

end subroutine execute_density

! =============================================================================

subroutine MPDATA_density_scheme(this,dt,pD,prho,prhofac,prhodt,pVn,pWn,kflip)
class(MPDATA_type), intent(inout) :: this
real(wp), intent(in) :: dt
real(wp), intent(inout) :: pD(:,:), pVn(:,:), pWn(:,:)
real(wp), dimension(:,:), intent(in) :: prho, prhofac, prhodt
integer, intent(in) :: kflip
type(atlas_HaloExchange)              :: halo_exchange
integer                               :: jnode, jedge
integer                               :: jlev
real(wp)                              :: zDmin (ubound(pD,1),ubound(pD,2))
real(wp)                              :: zDmax (ubound(pD,1),ubound(pD,2))
real(wp)                              :: zDR   (ubound(pD,1),ubound(pD,2))
real(wp)                              :: zdivVD(ubound(pD,1),ubound(pD,2))
real(wp)                              :: zVn(this%geom%nb_levels,this%geom%nb_edges)
real(wp)                              :: zWn(this%geom%nb_levels+1,this%geom%nb_nodes)
real(wp)                              :: zflux(this%geom%nb_levels,this%geom%nb_edges)
real(wp)                              :: zfluz(this%geom%nb_levels+1,this%geom%nb_nodes)
real(wp)                              :: rlimit
real(wp)                              :: rIVBZ
integer                               :: mporder
integer                               :: nb_nodes
integer                               :: nb_edges
integer                               :: nb_levels
logical                               :: logical_test

call atlas_log%debug('mpdata_gauge_density')

halo_exchange = this%geom%nodes%get_halo_exchange()

nb_nodes  = this%geom%nb_nodes
nb_edges  = this%geom%nb_edges
nb_levels = this%geom%nb_levels

if( this%rlimit >= 0._wp .and. this%mporder >= 2 ) then
  call this%limiter%compute_min_max_xy(pD,zDmax,zDmin,kflip)
  call this%limiter%compute_min_max_z (pD,zDmax,zDmin)
endif

call compute_upwind_flux(this,zVn,pD,pVn)
call compute_upwind_fluz(this,zWn,pD,pWn,this%rIVBZ)

call compute_fluxzdiv(this,zdivVD,zVn,zWn,this%geom%dual_volumes)
call advance_solution(this,dt,pD,zdivVD,prho)

call halo_exchange%execute(pD)

if (this%mporder > 1) then

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode, jlev)
  do jnode  = 1,nb_nodes
    do jlev = 1,nb_levels
      zDR(jlev,jnode) = pD(jlev,jnode)*prhofac(jlev,jnode)
    enddo
  enddo
!$OMP END PARALLEL DO

  call compute_centred_flux(this,zflux,zDR,pVn)
  call compute_centred_fluz(this,zfluz,zDR,pWn)
  call compute_fluxzdiv(this,zdivVD,zflux,zfluz,this%geom%dual_volumes)
  call halo_exchange%execute(zdivVD)

  call compute_pseudovel_xy(this,dt,zflux,zDR,pVn,zdivVD,prho,prhodt)
  call compute_pseudovel_z(this,dt,zfluz,zDR,pWn,zdivVD,prho,prhodt)

  if (this%rlimit >= 0._wp) then
    call this%limiter%limit_flux(zflux,zfluz,zDR,zDmax,zDmin,this%rlimit,dt,this%geom%dz,this%geom%dual_volumes,prho)
  endif

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jedge, jlev)
  do jedge = 1,nb_edges
    do jlev = 1,nb_levels
      pVn(jlev,jedge) = zVn(jlev,jedge)+zflux(jlev,jedge)
    enddo
  enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode, jlev)
  do jnode  = 1,nb_nodes
    do jlev = 1,nb_levels+1
      pWn(jlev,jnode) = zWn(jlev,jnode)+zfluz(jlev,jnode)
    enddo
  enddo
!$OMP END PARALLEL DO

  call compute_fluxzdiv(this,zdivVD,zflux,zfluz,this%geom%dual_volumes)
  call advance_solution(this,dt,pD,zdivVD,prho)

  call halo_exchange%execute(pD)

else

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jedge, jlev)
  do jedge = 1,nb_edges
    do jlev = 1,nb_levels
      pVn(jlev,jedge) = zVn(jlev,jedge)
    enddo
  enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode, jlev)
  do jnode  = 1,nb_nodes
    do jlev = 1,nb_levels+1
      pWn(jlev,jnode) = zWn(jlev,jnode)
    enddo
  enddo
!$OMP END PARALLEL DO

  call halo_exchange%execute(pD)

endif

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode, jlev)
do jnode  = 1,nb_nodes
  do jlev = 1,nb_levels
    pD(jlev,jnode) = pD(jlev,jnode)*prhofac(jlev,jnode)
  enddo
enddo
!$OMP END PARALLEL DO

end subroutine MPDATA_density_scheme

! =============================================================================

subroutine MPDATA_scheme(this,dt,pD,prho,prhofac,prhodt,pVn,pWn,kflip)
class(MPDATA_type), intent(inout) :: this
real(wp), intent(in) :: dt
real(wp), intent(inout) :: pD(:,:)
real(wp), dimension(:,:), intent(in) :: prho, prhofac, prhodt, pVn, pWn
integer, intent(in) :: kflip
type(atlas_HaloExchange)              :: halo_exchange
integer                               :: jnode
integer                               :: jlev
real(wp)                              :: zDmin (ubound(pD,1),ubound(pD,2))
real(wp)                              :: zDmax (ubound(pD,1),ubound(pD,2))
real(wp)                              :: zDR   (ubound(pD,1),ubound(pD,2))
real(wp)                              :: zdivVD(ubound(pD,1),ubound(pD,2))
real(wp)                              :: zflux(this%geom%nb_levels,this%geom%nb_edges)
real(wp)                              :: zfluz(this%geom%nb_levels+1,this%geom%nb_nodes)
type(atlas_Trace)                     ::trace
trace = atlas_Trace(__FILENAME__,__LINE__,"MPDATA_scheme")
call fckit_log%debug('mpdata_gauge_scalar')

halo_exchange = this%geom%nodes%get_halo_exchange()

if( this%rlimit >= 0._wp .and. this%mporder >= 2 ) then
  call this%limiter%compute_min_max_xy(pD,zDmax,zDmin,kflip)
  call this%limiter%compute_min_max_z (pD,zDmax,zDmin)
endif

call compute_upwind_flux(this,zflux,pD,pVn)
call compute_upwind_fluz(this,zfluz,pD,pWn,this%rIVBZ)

call compute_fluxzdiv(this,zdivVD,zflux,zfluz,this%geom%dual_volumes)
call advance_solution(this,dt,pD,zdivVD,prho)

call halo_exchange%execute(pD)

MPCORR: if (this%mporder > 1) then

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode, jlev)
  do jnode  = 1,this%geom%nb_nodes
    do jlev = 1,this%geom%nb_levels
      zDR(jlev,jnode) = pD(jlev,jnode)*prhofac(jlev,jnode)
    enddo
  enddo
!$OMP END PARALLEL DO

  call compute_centred_flux(this,zflux,zDR,pVn)
  call compute_centred_fluz(this,zfluz,zDR,pWn)
  call compute_fluxzdiv(this,zdivVD,zflux,zfluz,this%geom%dual_volumes)
  call halo_exchange%execute(zdivVD)

  call compute_pseudovel_xy(this,dt,zflux,zDR,pVn,zdivVD,prho,prhodt)
  call compute_pseudovel_z (this,dt,zfluz,zDR,pWn,zdivVD,prho,prhodt)

  if (this%rlimit >= 0._wp) then
    call this%limiter%limit_flux(zflux,zfluz,zDR,zDmax,zDmin,this%rlimit,dt,this%geom%dz,this%geom%dual_volumes,prho)
  endif

  call compute_fluxzdiv(this,zdivVD,zflux,zfluz,this%geom%dual_volumes)
  call advance_solution(this,dt,pD,zdivVD,prho)

  call halo_exchange%execute(pD)

endif MPCORR

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode, jlev)
  do jnode  = 1,this%geom%nb_nodes
    do jlev = 1,this%geom%nb_levels
      pD(jlev,jnode) = pD(jlev,jnode)*prhofac(jlev,jnode)
    enddo
  enddo
!$OMP END PARALLEL DO
call trace%final()

end subroutine MPDATA_scheme


!###############################################################################


subroutine compute_fluxzdiv(this,pdivVD,pFx,pFz,pvol)
type(MPDATA_type), intent(inout) :: this
real(wp), intent(in)  :: pFx(:,:), pFz(:,:)
real(wp), intent(in)  :: pvol(:)
real(wp), intent(out) :: pdivVD(:,:)
real(wp) :: zadd
integer :: jnode, jlev, jedge, iedge
integer, pointer                   :: inode2edge  (:)
integer                            :: inode2edge_size
type(atlas_Trace) :: trace
trace=atlas_Trace(__FILENAME__,__LINE__,"compute_fluxzdiv")

call fckit_log%debug('compute_fluxzdiv')

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,jlev,jedge,iedge,inode2edge,inode2edge_size,zadd)
do jnode  = 1,this%geom%nb_nodes
  call this%geom%node2edge%row(jnode,inode2edge,inode2edge_size)
  do jlev = 1,this%geom%nb_levels
    pdivVD(jlev,jnode) = 0.0_wp
    do jedge = 1,inode2edge_size
      iedge  = inode2edge(jedge)
      zadd  = real(this%geom%node2edge_sign(jedge,jnode),wp)
      pdivVD(jlev,jnode) = pdivVD(jlev,jnode)+zadd*pFx(jlev,iedge)
    enddo
    pdivVD(jlev,jnode) = pdivVD(jlev,jnode)/pvol(jnode)              &
      & +(pFz(jlev+1,jnode)-pFz(jlev,jnode))/this%geom%dz
  enddo
enddo
!$OMP END PARALLEL DO
call trace%final()
end subroutine compute_fluxzdiv

!###############################################################################

subroutine advance_solution(this,dt,pD,pdivVD,prho)
  !! Update solution to next time step

type(MPDATA_type), intent(inout) :: this
real(wp), intent(in) :: dt
real(wp), intent(inout) :: pD(:,:)
real(wp), intent(in)  :: pdivVD(:,:), prho(:,:)
integer :: jnode, jlev
type(atlas_Trace) :: trace
trace=atlas_Trace(__FILENAME__,__LINE__,"advance_solution")

call fckit_log%debug('advance_solution')

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode, jlev)
do jnode  = 1,this%geom%nb_nodes
  do jlev = 1,this%geom%nb_levels
    pD(jlev,jnode) = pD(jlev,jnode)-dt*pdivVD(jlev,jnode)/prho(jlev,jnode)
  enddo
enddo
!$OMP END PARALLEL DO
call trace%final()
end subroutine advance_solution

!###############################################################################

subroutine compute_upwind_flux(this,pflux,pD,pVn)
type(MPDATA_type), intent(inout) :: this
real(wp), intent(out) :: pflux(:,:)
real(wp), intent(in)  :: pVn(:,:), pD(:,:)
real(wp) :: zpos, zneg
integer  :: jedge, jlev, ip1, ip2

type(atlas_Trace) :: trace
trace=atlas_Trace(__FILENAME__,__LINE__,"compute_upwind_flux")

call fckit_log%debug('compute_upwind_flux')

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jedge,jlev,ip1,ip2,zpos,zneg)
do jedge = 1,this%geom%nb_edges
  ip1 = this%geom%iedge2node(1,jedge)
  ip2 = this%geom%iedge2node(2,jedge)
  do jlev = 1,this%geom%nb_levels
    zpos               = max(0._wp,pVn(jlev,jedge))
    zneg               = min(0._wp,pVn(jlev,jedge))
    pflux(jlev,jedge)  = pD(jlev,ip1)*zpos+pD(jlev,ip2)*zneg
  enddo
enddo
!$OMP END PARALLEL DO
call trace%final()
end subroutine compute_upwind_flux

!###############################################################################

subroutine compute_upwind_fluz(this,pfluz,pD,pW,pivbz)
type(MPDATA_type), intent(inout) :: this
real(wp), intent(out) :: pfluz(:,:)
real(wp), intent(in)  :: pW(:,:), pD(:,:), pivbz
integer  :: jnode, jlev

type(atlas_Trace) :: trace
trace=atlas_Trace(__FILENAME__,__LINE__,"compute_upwind_fluz")

call fckit_log%debug('compute_upwind_fluz')

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode, jlev)
do jnode  = 1,this%geom%nb_nodes
  do jlev = 2,this%geom%nb_levels
    pfluz(jlev,jnode) =  max(0._wp,pW(jlev,jnode))*pD(jlev-1,jnode) &
                      & +min(0._wp,pW(jlev,jnode))*pD( jlev ,jnode) 
  enddo
  pfluz(   1 ,jnode) = pivbz*pfluz( 2 ,jnode)
  pfluz(this%geom%nb_levels+1,jnode) = pivbz*pfluz(this%geom%nb_levels,jnode)
enddo
!$OMP END PARALLEL DO
call trace%final()
end subroutine compute_upwind_fluz

!###############################################################################

subroutine compute_centred_flux(this,pflux,pD,pVn)
type(MPDATA_type), intent(inout) :: this
real(wp), intent(out) :: pflux(:,:)
real(wp), intent(in)  :: pVn(:,:), pD(:,:)
integer  :: jedge, jlev, ip1, ip2

type(atlas_Trace) :: trace
trace=atlas_Trace(__FILENAME__,__LINE__,"compute_centred_flux")

call fckit_log%debug('compute_centred_flux')


!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jedge,jlev,ip1,ip2)
do jedge  = 1,this%geom%nb_edges
  ip1 = this%geom%iedge2node(1,jedge)
  ip2 = this%geom%iedge2node(2,jedge)
  do jlev = 1,this%geom%nb_levels
    pflux(jlev,jedge) = 0.5_wp*pVn(jlev,jedge)*(pD(jlev,ip1)+pD(jlev,ip2))
  enddo
enddo
!$OMP END PARALLEL DO
call trace%final()
end subroutine compute_centred_flux

!###############################################################################

subroutine compute_centred_fluz(this,pfluz,pD,pW)
type(MPDATA_type), intent(inout) :: this
real(wp), intent(out) :: pfluz(:,:)
real(wp), intent(in)  :: pW(:,:), pD(:,:)
integer  :: jnode, jlev

call fckit_log%debug('compute_centred_fluz')

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode, jlev)
do jnode  = 1,this%geom%nb_nodes
  do jlev = 2,this%geom%nb_levels
    pfluz(jlev,jnode) = 0.5_wp*pW(jlev,jnode)*(pD(jlev,jnode)+pD(jlev-1,jnode))
  enddo
  pfluz(   1 ,jnode)= this%rIVBZ*pfluz( 2 ,jnode)
  pfluz(this%geom%nb_levels+1,jnode)= this%rIVBZ*pfluz(this%geom%nb_levels,jnode)
enddo
!$OMP END PARALLEL DO

end subroutine compute_centred_fluz

!###############################################################################

subroutine compute_pseudovel_xy(this,dt,pflux,pD,pV,pdivVD,prho,prhodt)
type(MPDATA_type), intent(inout) :: this
real(wp), intent(in) :: dt
real(wp), dimension(:,:), intent(out)    :: pflux
real(wp), dimension(:,:), intent(in)     :: pV, pD, pdivVD, prho, prhodt
real(wp) :: zpos, zneg, zdDdt, zrhoav, zdVDf, zdGdt, zflux
integer  :: jedge, jlev, ip1, ip2

call fckit_log%debug('compute_pseudovel_xy')

!!$!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(ip1,ip2,jlev,jedge)
!!$do jedge  = 1,this%geom%nb_edges
!!$  ip1    = this%geom%iedge2node(1,jedge)
!!$  ip2    = this%geom%iedge2node(2,jedge)
!!$  do jlev = 1,this%geom%nb_levels
!!$    pflux(jlev,jedge) = 0.5_wp*abs(pV(jlev,jedge))*(pD(jlev,ip2)-pD(jlev,ip1))
!!$  enddo
!!$enddo
!!$!$OMP END PARALLEL DO

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jedge,jlev,ip1,ip2,zrhoav,zdVDf,zdGdt,zflux)
 do jedge = 1,this%geom%nb_edges
   ip1    = this%geom%iedge2node(1,jedge)
   ip2    = this%geom%iedge2node(2,jedge)
   do jlev=1,this%geom%nb_levels
     zrhoav            = prho(jlev,ip1)+prho(jlev,ip2)
     zdVDf             = 0.5_wp*(pdivVD(jlev,ip1)+pdivVD(jlev,ip2))
     zdGdt             = 0.25_wp*(pD(jlev,ip2)+pD(jlev,ip1))*(prhodt(jlev,ip2)+prhodt(jlev,ip1))
     zflux             = 0.5_wp*abs(pV(jlev,jedge))*(pD(jlev,ip2)-pD(jlev,ip1))
     pflux(jlev,jedge) = zflux-dt*pV(jlev,jedge)*(zdVDf+zdGdt)/zrhoav
   enddo
 enddo
!$OMP END PARALLEL DO

end subroutine compute_pseudovel_xy

!###############################################################################

subroutine compute_pseudovel_z(this,dt,pfluz,pD,pW,pdivVD,prho,prhodt)
type(MPDATA_type), intent(inout) :: this
real(wp), intent(in) :: dt
real(wp), dimension(:,:), intent(out)    :: pfluz
real(wp), dimension(:,:), intent(in)     :: pW, pD, pdivVD, prho, prhodt
real(wp) :: zWz, zdDdt, zrhoav, zdVDf, zdGdt, zfluz
integer  :: jnode, jlev

call fckit_log%debug('compute_pseudovel_z')

!!$!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode, jlev)
!!$do jnode  = 1,this%geom%nb_nodes
!!$  do jlev = 2,this%geom%nb_levels
!!$    pfluz(jlev,jnode) = 0.5_wp*abs(pW(jlev,jnode))*(pD(jlev,jnode)-pD(jlev-1,jnode))
!!$  enddo
!!$enddo
!!$!$OMP END PARALLEL DO

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,jlev,zWz,zrhoav,zdVDf,zdGdt,zfluz)
do jnode  = 1,this%geom%nb_nodes
  do jlev = 2,this%geom%nb_levels
    zWz               = pW(jlev,jnode)
    zrhoav            = (prho(jlev,jnode)+prho(jlev-1,jnode))
    zdVDf             = 0.5_wp*(pdivVD(jlev,jnode)+pdivVD(jlev-1,jnode))
    zdGdt             = 0.25_wp*(pD(jlev,jnode)+pD(jlev-1,jnode))*(prhodt(jlev,jnode)+prhodt(jlev-1,jnode))
    zfluz             = 0.5_wp*abs(pW(jlev,jnode))*(pD(jlev,jnode)-pD(jlev-1,jnode))
    pfluz(jlev,jnode) = zfluz-dt*zWz*(zdVDf+zdGdt)/zrhoav
  enddo
enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode)
do jnode = 1,this%geom%nb_nodes
  pfluz(  1  ,jnode)       = this%rIVBZ*pfluz( 2 ,jnode)
  pfluz(this%geom%nb_levels+1,jnode) = this%rIVBZ*pfluz(this%geom%nb_levels,jnode)
enddo
!$OMP END PARALLEL DO

end subroutine compute_pseudovel_z

!###############################################################################

subroutine rhofac_correction(this,pD1,pD2,prf)
type(MPDATA_type), intent(inout) :: this
real(wp), intent(out) :: pD1(:,:)
real(wp), intent(inout)  :: pD2(:,:), prf(:,:)
integer  :: jnode, jlev

call fckit_log%debug('rhofac_correction')

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode, jlev)
do jnode  = 1,this%geom%nb_nodes
  do jlev = 1,this%geom%nb_levels
    pD1(jlev,jnode) = pD2(jlev,jnode)*prf(jlev,jnode)
  enddo
enddo
!$OMP END PARALLEL DO

end subroutine rhofac_correction

!###############################################################################

subroutine face_normals(this,field_V,field_Vn)
class(MPDATA_type), intent(inout) :: this
type(atlas_Field), intent(in)     :: field_V  !! input
type(atlas_Field), intent(inout)  :: field_Vn !! output
integer  :: jedge, jlev
real(wp), pointer :: pVn(:,:)
real(wp), pointer :: pV(:,:,:)
type(atlas_Trace) :: trace
trace=atlas_Trace(__FILENAME__,__LINE__,"advectors_in_edges")

call fckit_log%debug('compute_face_normals')

call field_Vn%data(pVn)
call field_V%data(pV)

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jedge, jlev)
do jedge  = 1,this%geom%nb_edges
  do jlev = 1,this%geom%nb_levels
    pVn(jlev,jedge) = dot_product(pV(:,jlev,jedge),this%geom%dual_normals(:,jedge))
  enddo
enddo
!$OMP END PARALLEL DO
 call trace%final()
end subroutine face_normals

!###############################################################################

subroutine advectors_in_edges(this,field_Vnodes,field_Vxy,field_Vz)
class(MPDATA_type), intent(inout) :: this
type(atlas_Field), intent(in)  :: field_Vnodes
type(atlas_Field), intent(inout) :: field_Vxy, field_Vz
integer :: jnode, jedge, iedge, jlev, ip1, ip2
real(wp), pointer :: Vnodes(:,:,:)
real(wp), pointer :: VXY(:,:,:), VZ(:,:)
type(atlas_Trace) :: trace
trace=atlas_Trace(__FILENAME__,__LINE__,"advectors_in_edges")
call fckit_log%debug('set_advectors_in_edges')

call field_Vnodes%data(Vnodes)
call field_Vxy%data(Vxy)
call field_Vz%data(Vz)

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jedge,jlev,ip1,ip2)
do jedge = 1,this%geom%nb_edges
  ip1  = this%geom%iedge2node(1,jedge)
  ip2  = this%geom%iedge2node(2,jedge)
  do jlev = 1,this%geom%nb_levels
    VXY(MXX,jlev,jedge) = (Vnodes(MXX,jlev,ip1)       &
      & +this%geom%rpole_bc(jedge)*Vnodes(MXX,jlev,ip2))*0.5_wp
    VXY(MYY,jlev,jedge) = (Vnodes(MYY,jlev,ip1)       &    
      & +this%geom%rpole_bc(jedge)*Vnodes(MYY,jlev,ip2))*0.5_wp
  enddo
enddo
!$OMP END PARALLEL DO

! Since the pole point lies outside the lon-lat domain, Vedges is wrongly
! calculated
! y_pole .ne. 0.5(y1+y2)

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jedge,iedge,jlev,ip1,ip2)
do jedge = 1,this%geom%nb_pole_edges
  iedge = this%geom%pole_edges(jedge)
  ip1  = this%geom%iedge2node(1,iedge)
  ip2  = this%geom%iedge2node(2,iedge)
  do jlev = 1,this%geom%nb_levels
    VXY(MYY,jlev,iedge) = 0._wp
  enddo
enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode, jlev)
do jnode  = 1,this%geom%nb_nodes
  do jlev = 2,this%geom%nb_levels
    VZ(jlev,jnode) = 0.5_wp*(Vnodes(MZZ,jlev,jnode)+Vnodes(MZZ,jlev-1,jnode))
  enddo
enddo
!$OMP END PARALLEL DO
call trace%final()
end subroutine advectors_in_edges

!###############################################################################

end module advection_MPDATA_MPDATA_module
