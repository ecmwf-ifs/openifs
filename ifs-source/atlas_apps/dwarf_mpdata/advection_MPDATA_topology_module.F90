! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

#define __FILENAME__ "advection_MPDATA_topology_module.F90"
module advection_MPDATA_topology_module

use atlas_module
use fckit_module, only : fckit_log
use advection_MPDATA_auxiliary_module
use advection_MPDATA_geometry_module, only : geometry_type

implicit none


private
public :: topology_type

! Class topology
type :: topology_type
  type(geometry_type), pointer, private :: geom
  integer , private :: vstretch
contains
  generic  , public          :: assignment(=) => copy
  procedure, private         :: copy
  procedure, public          :: execute
  procedure, private         :: compute_topology
  procedure, private, nopass :: compute_vertical_levels
end type

! Definition of the constructor
interface topology_type
  module procedure topology_constructor
end interface topology_type

contains


! =============================================================================
! Constructor for class topology
! =============================================================================
function topology_constructor(geometry, config) result(this)
type(topology_type) :: this

type(geometry_type), target    :: geometry
type(atlas_config), intent(in) :: config

logical  :: logical_test

this%geom => geometry

! Load parameters from config
if( .not. config%get("vstretch" , this%vstretch ) ) call fckit_log%error("vstretch not found in config")


end function topology_constructor
! =============================================================================

subroutine copy(obj_out,obj_in)
  class(topology_type), intent(inout) :: obj_out
  class(topology_type), intent(in) :: obj_in
  obj_out%geom          => obj_in%geom
  obj_out%vstretch      =  obj_in%vstretch
end subroutine copy

! =============================================================================
! Subroutine execute: driver to call core topology routine
! =============================================================================
subroutine execute(this, pheight_field, topology_fields)
class(topology_type), intent(inout) :: this
type(atlas_field),    intent(in)    :: pheight_field
type(atlas_fieldset), intent(inout) :: topology_fields

type(atlas_field)                   :: rzcr_field
type(atlas_Trace)                   :: trace
trace = atlas_Trace(__FILENAME__,__LINE__,"topology_type::execute")
if( .not. topology_fields%has_field("rzcr") ) then
  call topology_fields%add(                                               &
    &    this%geom%nodes%create_field( name="rzcr" , kind=atlas_real(wp), &
    &                                  levels=this%geom%nb_levels) )
endif
rzcr_field = topology_fields%field("rzcr")

call this%compute_topology(pheight_field,rzcr_field)

call rzcr_field%final()
call trace%final()

end subroutine execute
! =============================================================================



! =============================================================================
! SET TERRAIN TOPOLOGY
! =============================================================================
subroutine compute_topology(this, pheight_field, rzcr_field)
class(topology_type), intent(inout) :: this
type(atlas_Field), intent(in)  :: pheight_field
type(atlas_Field), intent(inout) :: rzcr_field
!-----------------------------------------------------
real(wp), pointer              :: rzcr(:,:)
real(wp), pointer              :: zheight(:,:)
!-----------------------------------------------------


call rzcr_field%data(rzcr)
call pheight_field%data(zheight)

call compute_vertical_levels(this%geom%nb_nodes, this%geom%nb_levels, this%geom%dz, zheight, rzcr, this%vstretch)

call this%geom%nodes%halo_exchange(rzcr_field)

end subroutine compute_topology
! =============================================================================



! =============================================================================
! SET THE VERTICAL LEVELS USING TERRAIN FOLLOWING COORDINATES
! =============================================================================
subroutine compute_vertical_levels(&
& nb_nodes, nb_levels, dz, zheight, rzcr, vstretch)

! Dummy parameters
integer , intent(in)  :: nb_nodes
integer , intent(in)  :: nb_levels
real(wp), intent(in)  :: dz
real(wp), intent(in)  :: zheight(:,:)
real(wp), intent(out) :: rzcr(:,:)
integer , intent(in)  :: vstretch

! Local variables
real(wp) :: zb (nb_levels)
real(wp) :: dzb(2:nb_levels)
real(wp) :: zbb(ubound(zb,1))
real(wp) :: ztrans
integer  :: jnode
integer  :: jlev

call fckit_log%debug('compute_vertical_levels')

!!$! Define levels based on vertical spacing dz
!!$!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jlev)
!!$do jlev = 1,nb_levels
!!$  zb(jlev) = real(jlev-1,wp) * dz
!!$enddo
!!$!$OMP END PARALLEL DO
!!$
!!$! If required, apply stretching in the vertical direction
!!$if (vstretch == 1) then
!!$  zb = compute_vertical_stretching(zb, zb(nb_levels), 2.2_wp, 'tanh')
!!$endif
!!$
!!$! Calculate the new vertical spacing dzb
!!$dzb(2:nb_levels) = zb(2:nb_levels)-zb(1:nb_levels-1)
!!$
!!$! Calculate the vertical levels using terrain following coordinates
!!$ztrans = zb(nb_levels)
!!$
!!$! Define terrain following coordinates
!!$where (zb <= ztrans)
!!$  zbb = 1.0_wp - zb / ztrans
!!$elsewhere
!!$  zbb = 0._wp
!!$end where

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,jlev)
do jnode = 1,nb_nodes
  do jlev = 1,nb_levels
    rzcr(jlev,jnode) = zheight(jlev,jnode) !zb(jlev)
  enddo
enddo
!$OMP END PARALLEL DO

end subroutine compute_vertical_levels
! =============================================================================



! =============================================================================
! STRETCH VERTICAL COORDINATES
! =============================================================================
elemental function compute_vertical_stretching(pz, ptop, pscale, cdmethod) &
& result(stretch)

real(wp) :: stretch

! Dummy parameters
real(wp)        , intent(in) :: pz
real(wp)        , intent(in) :: ptop
real(wp)        , intent(in) :: pscale
character(len=*), intent(in) :: cdmethod

select case (trim(cdmethod))
  case('exp')
    stretch = -pscale * log(1._wp - pz / ptop * (1._wp - exp(-ptop / pscale)))
  case('tanh')
    stretch = ptop*(1._wp + (tanh(pscale * (pz / ptop - 1._wp))) / &
                                                        & (tanh(pscale)))
  case default
    stretch = pz
end select

end function compute_vertical_stretching
! =============================================================================








end module advection_MPDATA_topology_module
