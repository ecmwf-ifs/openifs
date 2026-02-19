! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

! =============================================================================
! geometry_module
! This module contains geometry
! =============================================================================

module advection_MPDATA_geometry_module

use atlas_module
use fckit_module, only : fckit_log, fckit_mpi_comm, fckit_main, fckit_resource
use advection_MPDATA_auxiliary_module

implicit none

private

type, public :: geometry_type
  type(atlas_connectivity)             , public  :: edge2node
  type(atlas_connectivity)             , public  :: node2edge
  type(atlas_Mesh)                     , public  :: mesh
  type(atlas_functionspace_nodecolumns), public  :: nodes
  type(atlas_functionspace_edgecolumns), public  :: edges
  type(atlas_Field)                    , public  :: field_pole_bc
  type(atlas_Field)                    , public  :: field_lonlat
  type(atlas_Field)                    , public  :: field_dual_volumes
  type(atlas_Field)                    , public  :: field_dual_normals
  type(atlas_FieldSet)                 , public  :: topology_fields

  real(wp), pointer                , public  :: node2edge_sign(:,:) ! only computed data local to type !
  real(wp), pointer                , public  :: rpole_bc(:)
  integer,  pointer                , public  :: is_pole_edge(:)
  integer,  pointer                , public  :: pole_edges(:)
  integer,  pointer                , public  :: iedge2node(:,:)
  real(wp), pointer                , public  :: dual_volumes(:)
  real(wp), pointer                , public  :: dual_normals(:,:)
  real(wp), pointer                , public  :: lonlat(:,:)
  
  integer, public  :: nb_edges
  integer, public  :: nb_nodes
  integer, public  :: nb_levels
  integer, public  :: nb_pole_edges
  real(wp), public :: dz

contains
  generic  , public  :: assignment(=) => copy
  procedure, private :: copy
  procedure, public :: final
end type

interface geometry_type
  module procedure :: geometry_constructor
end interface

contains

! =============================================================================

function geometry_constructor(fvm, mesh, config) result(this)
  type(geometry_type) :: this
  type(atlas_fvm_Method), intent(in)    :: fvm
  type(atlas_Mesh),       intent(in)    :: mesh
  type(atlas_config),     intent(in)    :: config
  
  integer                               :: nb_levels

  ! Additional parameters
  integer                               :: jnode
  integer                               :: jlev
  integer                               :: nb_nodes
  integer                               :: nb_edges

  type(atlas_mesh_nodes)                :: mesh_nodes
  type(atlas_mesh_edges)                :: mesh_edges
  type(atlas_Field)                     :: field
  real(dp), dimension(:,:), pointer     :: lonlat
  type(atlas_field)                     :: field_lonlat
  integer, allocatable                  :: tmp(:)
  integer                               :: jedge
  integer                               :: iedge
  integer                               :: ip1
  integer                               :: ip2
  integer                               :: max_nb_neighbours
  real(dp), dimension(:)  , pointer     :: dual_volumes
  real(dp), dimension(:,:), pointer     :: dual_normals
  integer, pointer                      :: inode2edge  (:)
  integer                               :: inode2edge_size
  real(dp)                              :: radius
  logical :: found
  
  radius = 6371229._dp ! default value
  found = config%get("radius", radius)
  if( .not. config%get("nb_levels", this%nb_levels) ) call fckit_log%error("nb_levels not found in config")
  if( .not. config%get("dz",        this%dz)        ) call fckit_log%error("dz not found in config")

  this%mesh  = mesh

  this%nodes = fvm%node_columns()
  this%edges = fvm%edge_columns()

  mesh_nodes = this%nodes%nodes()
  mesh_edges = this%edges%edges()

  ! Create copy of Atlas fields in working precision
  this%field_lonlat = this%nodes%create_field(name="lonlat", kind=atlas_real(wp), variables=2)
  this%field_dual_volumes = this%nodes%create_field(name="dual_volumes"  , kind=atlas_real(wp))
  this%field_dual_normals = this%edges%create_field(name="dual_normals"  , kind=atlas_real(wp), variables=2)
  call this%field_lonlat%data(this%lonlat)
  call this%field_dual_volumes%data(this%dual_volumes)
  call this%field_dual_normals%data(this%dual_normals)

  ! Number of nodes and edges
  nb_nodes = mesh_nodes%size()
  nb_edges = mesh_edges%size()

  ! Edge to node connectivity
  this%edge2node = mesh_edges%node_connectivity()
  call this%edge2node%data(this%iedge2node)

  ! Node to edge connectivity
  this%node2edge = mesh_nodes%edge_connectivity()

  ! Node to node connectivity support
  max_nb_neighbours = this%node2edge%maxcols()
  allocate(this%node2edge_sign(max_nb_neighbours, nb_nodes))

  ! Compute node to node connectivity
  !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,jedge,iedge,ip1,ip2,inode2edge,inode2edge_size)
  do jnode=1,nb_nodes
    call this%node2edge%row(jnode,inode2edge,inode2edge_size)
    do jedge=1,inode2edge_size
      iedge = inode2edge(jedge)
      ip1 = this%iedge2node(1,iedge)
      ip2 = this%iedge2node(2,iedge)
      if (jnode == ip1) then
        this%node2edge_sign(jedge,jnode) = 1._wp
      else
        this%node2edge_sign(jedge,jnode) = -1._wp
      endif
    enddo
  enddo
  !$OMP END PARALLEL DO

  ! Boolean to identify polar edges
  field = mesh_edges%field("is_pole_edge")
  call field%data(this%is_pole_edge)
  call field%final()

  ! Compute polar edges
  this%nb_pole_edges = 0
  allocate(tmp(nb_edges))
  do jedge=1,nb_edges
    if (this%is_pole_edge(jedge) == 1) then
      this%nb_pole_edges = this%nb_pole_edges + 1
      tmp(this%nb_pole_edges) = jedge
    endif
  enddo
  allocate(this%pole_edges(this%nb_pole_edges))
  this%pole_edges(:) = tmp(1:this%nb_pole_edges)
  deallocate(tmp)

  ! Dual normals
  field = mesh_edges%field('dual_normals')
  call field%data(dual_normals)
  !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jedge)
  do jedge = 1,nb_edges
    this%dual_normals(:,jedge) = dual_normals(:,jedge) * radius * deg2rad
  enddo
  !$OMP END PARALLEL DO
  call field%final()

  ! Retrieve lonlat coordinates
  field_lonlat = mesh_nodes%lonlat()
  call field_lonlat%data(lonlat)
  !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode)
  do jnode=1,nb_nodes
    this%lonlat(:,jnode) = lonlat(:,jnode) * deg2rad
  enddo

  ! Define dual volumes with correct dimensions
  field = mesh_nodes%field('dual_volumes')
  call field%data(dual_volumes)
  call field%final()
  !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode)
  do jnode  = 1,nb_nodes 
    this%dual_volumes(jnode) = dual_volumes(jnode) * radius * radius * &
                        & deg2rad * deg2rad
  enddo
  !$OMP END PARALLEL DO

  ! Pole BCs
  this%field_pole_bc = this%edges%create_field(name="pole_bc", kind=atlas_real(wp))
  call this%field_pole_bc%data(this%rpole_bc)
  this%rpole_bc(:) = 1._wp
  do jedge = 1,this%nb_pole_edges
    iedge = this%pole_edges(jedge)
    this%rpole_bc(iedge) = -1._wp
  enddo

  ! Destroy object lonlat
  call field_lonlat%final()

  this%nb_nodes = nb_nodes
  this%nb_edges = nb_edges

  this%topology_fields = atlas_FieldSet()
  
!  call atlas_write_load_balance_report(this%mesh, "balance_report")

  
end function


subroutine copy(obj_out,obj_in)
  class(geometry_type), intent(inout) :: obj_out
  class(geometry_type), intent(in) :: obj_in
  
  obj_out%edge2node          =  obj_in%edge2node
  obj_out%node2edge          =  obj_in%node2edge
  obj_out%mesh               =  obj_in%mesh
  obj_out%nodes              =  obj_in%nodes
  obj_out%edges              =  obj_in%edges
  obj_out%field_pole_bc      =  obj_in%field_pole_bc
  obj_out%field_lonlat       =  obj_in%field_lonlat
  obj_out%field_dual_volumes =  obj_in%field_dual_volumes
  obj_out%field_dual_normals =  obj_in%field_dual_normals
  obj_out%topology_fields    =  obj_in%topology_fields
  obj_out%rpole_bc           => obj_in%rpole_bc
  obj_out%node2edge_sign     => obj_in%node2edge_sign
  obj_out%is_pole_edge       => obj_in%is_pole_edge
  obj_out%pole_edges         => obj_in%pole_edges
  obj_out%iedge2node         => obj_in%iedge2node
  obj_out%dual_volumes       => obj_in%dual_volumes
  obj_out%dual_normals       => obj_in%dual_normals
  obj_out%lonlat             => obj_in%lonlat
  obj_out%nb_edges           =  obj_in%nb_edges
  obj_out%nb_nodes           =  obj_in%nb_nodes
  obj_out%nb_levels          =  obj_in%nb_levels
  obj_out%nb_pole_edges      =  obj_in%nb_pole_edges
  obj_out%dz                 =  obj_in%dz

end subroutine copy

subroutine final(this)
  class(geometry_type) :: this
  
  call this%topology_fields%final()
  call this%mesh%final()
  call this%nodes%final()
  call this%edges%final()
  call this%field_lonlat%final()
  call this%field_pole_bc%final()
  call this%field_dual_normals%final()
  call this%field_dual_volumes%final()
  call this%node2edge%final()
  call this%edge2node%final()
  ! deallocate(this%node2edge_sign) ! memory leak for now !!!!!!!
end subroutine

end module
