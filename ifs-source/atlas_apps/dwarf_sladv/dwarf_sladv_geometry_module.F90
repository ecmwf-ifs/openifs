! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

#include "dwarf_sladv.h"

module dwarf_sladv_geometry_module
  !------------------------------------------------------------------------------
  ! Define mesh structure holding geometry and various metric terms used in
  ! calculations.
  !
  ! The horizontal data structure is 1-dimensional
  ! There is a 1-1 mapping between a node (gridpoint) value and a 2-dimensional
  ! lat-lon matrix:
  !    ij2node(ilat,jlon)=jnode
  !    ilat(jnode): row index
  !    jlon(jnode): column index for jnode
  !
  ! Author: M. Diamantakis & W. Deconinck June 2018
  ! Modified by G. Tumolo, August 2018
  !------------------------------------------------------------------------------
  use atlas_module, only : atlas_Field, atlas_functionspace_StructuredColumns, &
  &                        atlas_Vertical
  use dwarf_sladv_kind_module, only : wp, ip

  implicit none

  private :: wp, ip
  private :: atlas_Field, atlas_functionspace_StructuredColumns

  private

  type, public :: Geometry_type
    real(wp)    :: radius               !!< Radius of planet
    integer(ip) :: ngp                  !!< Internal points
    integer(ip) :: ngptot               !!< All points including halo
    integer(ip) :: nlev                 !!< Number of vertical levels

    type(atlas_Field) :: glon
    type(atlas_Field) :: glat
    type(atlas_Field) :: zcoord
    type(atlas_functionspace_StructuredColumns) :: fs_structuredcolumns

  contains

    procedure, public :: setup
    procedure, public :: final
    final :: final_auto

  end type

contains


  subroutine setup(this,config,functionspace)
    use iso_c_binding
    use fckit_module
    use atlas_module, only : atlas_Field, atlas_functionspace_StructuredColumns, &
    &                        atlas_Config, atlas_real
    use dwarf_sladv_parameters_module
    !---------------------------------------------------------------------
    implicit none
    !---------------------------------------------------------------------
    class(Geometry_type), intent(inout) :: this
    type(atlas_Config), intent(in) :: config
    type(atlas_functionspace_StructuredColumns), intent(inout) :: functionspace
    !---------------------------------------------------------------------
    integer(ip) :: jnode, jlev
    real(wp), pointer :: gglon(:), gglat(:)
    character(len=1024) :: logmsg
    real(wp), pointer :: lonlat(:,:)
    logical :: found
    type(atlas_Field) :: lonlat_field

    !---------------------------------------------------------------------

    this%fs_structuredcolumns = functionspace

    this%radius = earth_radius! default value
    found = config%get("radius", this%radius)

    !----number of computational and halo points is computed---:
    this%ngp=this%fs_structuredcolumns%size_owned()     ! my computational points
    this%ngptot=this%fs_structuredcolumns%size() ! my computational points + my halo (interprocess communication halo + intraprocess polar halo) points

    ! vertical levels settings
    this%nlev=this%fs_structuredcolumns%levels()

    ! numerical parameter setting
    this%glon    = atlas_Field( "glon", atlas_real(wp), [this%ngptot] )
    this%glat    = atlas_Field( "glat", atlas_real(wp), [this%ngptot] )
    this%zcoord  = this%fs_structuredcolumns%z()

    call this%glon%data(gglon)
    call this%glat%data(gglat)

    !-----------------------------------------------------------------
    ! Set grid coordinates.
    !-----------------------------------------------------------------
    lonlat_field = this%fs_structuredcolumns%xy() ! these nodes include halos
    call lonlat_field%data(lonlat)

    do jnode=1,this%fs_structuredcolumns%size()
      gglat(jnode)=lonlat(2,jnode)*Deg2Rad
      gglon(jnode)=lonlat(1,jnode)*Deg2Rad
    enddo
    
    call lonlat_field%final()

  end subroutine

  subroutine final(this)
    class(Geometry_type) :: this
    call this%glon%final()
    call this%glat%final()
    call this%zcoord%final()
    call this%fs_structuredcolumns%final()
  end subroutine

  subroutine final_auto(this)
    type(Geometry_type) :: this
    call this%final()
  end subroutine

end module dwarf_sladv_geometry_module
