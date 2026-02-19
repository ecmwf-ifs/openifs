! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

#define __FILENAME__ "advection_MPDATA_mappings_module.F90"
! =============================================================================
! mappings_module
! This module contains strictly the mappings required by MPDATA
! - compute_metrics     : compute the metric terms for sphere
! - compute_coefficients: compute the coefficient matrix C for calculating
!                         the gradient of a potential (V - C\nabla\phi)
! - transform_vector    : transform a generic vector backward and forward
! =============================================================================
module advection_MPDATA_mappings_module

use fckit_module, only : fckit_log
use atlas_module
use advection_MPDATA_auxiliary_module
use advection_MPDATA_nabla_module
use advection_MPDATA_geometry_module, only : geometry_type

implicit none



private
public :: mappings_type

! Class mappings
type :: mappings_type
  type(geometry_type), pointer , private :: geom
  type(nabla_type),    private :: nabla
    
  type(atlas_field), public :: jacobian_determinant
  type(atlas_field), public :: jacobian_matrix_field
  type(atlas_field), public :: coefficients_matrix_field

contains

  procedure, public  :: execute
  procedure, public  :: transform_vector
  procedure, private :: compute_metrics
  procedure, private :: compute_coefficients
  procedure, private :: create_fields

  generic  , public  :: assignment(=) => copy
  procedure, private :: copy
  procedure, public  :: final

end type

! Definition of the constructor
interface mappings_type
  module procedure mappings_constructor
end interface mappings_type

#define INTERNAL_INDEX(i,j) 3*(j-1)+i
integer, parameter :: idx_11 = INTERNAL_INDEX(1,1)
integer, parameter :: idx_12 = INTERNAL_INDEX(1,2)
integer, parameter :: idx_13 = INTERNAL_INDEX(1,3)
integer, parameter :: idx_21 = INTERNAL_INDEX(2,1)
integer, parameter :: idx_22 = INTERNAL_INDEX(2,2)
integer, parameter :: idx_23 = INTERNAL_INDEX(2,3)
integer, parameter :: idx_31 = INTERNAL_INDEX(3,1)
integer, parameter :: idx_32 = INTERNAL_INDEX(3,2)
integer, parameter :: idx_33 = INTERNAL_INDEX(3,3)


contains


! =============================================================================
! Constructor for class mappings
! =============================================================================
function mappings_constructor(geometry, nabla, config) &
& result(this)
type(mappings_type) :: this

type(geometry_type), target    :: geometry
type(nabla_type),   intent(in) :: nabla
type(atlas_config), intent(in) :: config

this%geom   => geometry
this%nabla  =  nabla

call this%create_fields()

end function mappings_constructor
! =============================================================================

subroutine copy(obj_out,obj_in)
  class(mappings_type), intent(inout) :: obj_out
  class(mappings_type), intent(in) :: obj_in
  obj_out%geom                       => obj_in%geom
  obj_out%nabla                      =  obj_in%nabla
  obj_out%jacobian_determinant       =  obj_in%jacobian_determinant
  obj_out%jacobian_matrix_field      =  obj_in%jacobian_matrix_field
  obj_out%coefficients_matrix_field  =  obj_in%coefficients_matrix_field
end subroutine copy

subroutine create_fields(this)
class(mappings_type), intent(inout) :: this
integer                             :: nb_levels

nb_levels = this%geom%nb_levels

this%jacobian_determinant =                                         &
  &       this%geom%nodes%create_field(                             &
  &                                    name="jacobian_determinant", &
  &                                    kind=atlas_real(wp)        , &
  &                                    levels=nb_levels)

this%jacobian_matrix_field =                                        &
  &       this%geom%nodes%create_field(                             &
  &                                    name="jacobian_matrix",      &
  &                                    kind=atlas_real(wp)   ,      &
  &                                    levels=nb_levels      ,      &
  &                                    variables=9 )

this%coefficients_matrix_field =                                    &
  &       this%geom%nodes%create_field(                             &
  &                                    name="coefficients_matrix",  &
  &                                    kind=atlas_real(wp)       ,  &
  &                                    levels=nb_levels          ,  &
  &                                    variables=9 )
end subroutine


! =============================================================================
! Subroutine execute: driver to call core mappings routine
! =============================================================================
subroutine execute(this)
class(mappings_type), intent(inout) :: this

! Metric and coefficient matrices and scalar jacobian
real(wp), pointer :: Jdet(:,:)
real(wp), pointer :: jacobian_matrix(:,:,:)
real(wp), pointer :: coefficients_matrix(:,:,:)

type(atlas_Trace) :: trace
trace = atlas_Trace(__FILENAME__,__LINE__,"mappings_type::execute")

call this%jacobian_determinant      % data(Jdet)
call this%jacobian_matrix_field     % data(jacobian_matrix)
call this%coefficients_matrix_field % data(coefficients_matrix)

call this%compute_metrics(Jdet, jacobian_matrix)
call this%compute_coefficients(jacobian_matrix, coefficients_matrix)

call trace%final()
end subroutine
! =============================================================================

! =============================================================================
! COMPUTE METRICS
! =============================================================================
subroutine compute_metrics(this, Jdet, jacobian_matrix)
class(mappings_type), intent(inout) :: this
real(wp), intent(out)               :: Jdet(:,:)
real(wp), intent(out)               :: jacobian_matrix(:,:,:)

! Local variables
integer                             :: jnode
integer                             :: jlev
type(atlas_field)                   :: rzcr_field
real(wp), pointer                   :: rzcr(:,:)

real(wp)                            :: rczx(this%geom%nb_levels,this%geom%nb_nodes)
real(wp)                            :: rczy(this%geom%nb_levels,this%geom%nb_nodes)
real(wp)                            :: rczz(this%geom%nb_levels,this%geom%nb_nodes)
real(wp), allocatable               :: zgrzcr(:,:,:)
real(wp)                            :: zgo11
real(wp)                            :: zgo22
real(wp)                            :: zgo33

! Calculate metric terms
call fckit_log%debug('compute_metrics')

allocate(zgrzcr(3,this%geom%nb_levels,this%geom%nb_nodes))

rzcr_field = this%geom%topology_fields%field("rzcr")
call rzcr_field%data(rzcr)

! Calculate gradient of vertical coordinate
call this%nabla%gradient_3d(rzcr, zgrzcr)

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,jlev)
do jnode = 1,this%geom%nb_nodes
  do jlev = 1,this%geom%nb_levels
    Jdet(jlev,jnode) = zgrzcr(MZZ,jlev,jnode)

    rczx(jlev,jnode) = -zgrzcr(MXX,jlev,jnode) / Jdet(jlev,jnode)
    rczy(jlev,jnode) = -zgrzcr(MYY,jlev,jnode) / Jdet(jlev,jnode)
    rczz(jlev,jnode) = 1._wp / zgrzcr(MZZ,jlev,jnode)
  enddo
enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,jlev,zgo11,zgo22,zgo33)
do jnode = 1,this%geom%nb_nodes
  do jlev = 1,this%geom%nb_levels

    ! Auxiliary terms
    zgo11 = 1._wp / cos(this%geom%lonlat(MYY,jnode))
    zgo22 = 1._wp
    zgo33 = 1._wp

    ! Metric terms (terms of the Jacobian matrix)
    jacobian_matrix(idx_11,jlev,jnode) = zgo11
    jacobian_matrix(idx_12,jlev,jnode) = 0._wp
    jacobian_matrix(idx_21,jlev,jnode) = 0._wp
    jacobian_matrix(idx_22,jlev,jnode) = zgo22
    jacobian_matrix(idx_13,jlev,jnode) = rczx(jlev,jnode) * zgo11
    jacobian_matrix(idx_23,jlev,jnode) = rczy(jlev,jnode) * zgo22
    jacobian_matrix(idx_33,jlev,jnode) = rczz(jlev,jnode) * zgo33

    ! Determinant of the Jacobian matrix
    Jdet(jlev,jnode) = Jdet(jlev,jnode) * cos(this%geom%lonlat(MYY,jnode))
  enddo
enddo
!$OMP END PARALLEL DO

deallocate(zgrzcr)
call rzcr_field%final()

end subroutine compute_metrics
! =============================================================================



! =============================================================================
! COMPUTE COEFFICIENTS
! =============================================================================
subroutine compute_coefficients(this, jacobian_matrix, coefficient_matrix)
class(mappings_type), intent(inout) :: this
real(wp), intent(in)                :: jacobian_matrix     (:,:,:)
real(wp), intent(out)               :: coefficient_matrix  (:,:,:)

real(wp) :: za11, za12, za13
real(wp) :: za21, za22, za23
real(wp) :: za31, za32, za33

integer :: jnode
integer :: jlev

call fckit_log%debug('compute_coefficients')

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,jlev,&
!$OMP & za11,za12,za13,za21,za22,za23,za31,za32,za33)
do jnode = 1,this%geom%nb_nodes
  do jlev = 1,this%geom%nb_levels

    za11 = jacobian_matrix(idx_11,jlev,jnode)
    za12 = jacobian_matrix(idx_12,jlev,jnode)
    za13 = jacobian_matrix(idx_13,jlev,jnode)
    za21 = jacobian_matrix(idx_21,jlev,jnode)
    za22 = jacobian_matrix(idx_22,jlev,jnode)
    za23 = jacobian_matrix(idx_23,jlev,jnode)
    za31 = 0._wp
    za32 = 0._wp
    za33 = jacobian_matrix(idx_33,jlev,jnode)

    coefficient_matrix(idx_11,jlev,jnode) = &
                            & jacobian_matrix(idx_11,jlev,jnode) * za11 + &
                            & jacobian_matrix(idx_21,jlev,jnode) * za21

    coefficient_matrix(idx_12,jlev,jnode) = &
                            & jacobian_matrix(idx_12,jlev,jnode) * za12 + &
                            & jacobian_matrix(idx_21,jlev,jnode) * za22

    coefficient_matrix(idx_13,jlev,jnode) = &
                            & jacobian_matrix(idx_11,jlev,jnode) * za13 + &
                            & jacobian_matrix(idx_21,jlev,jnode) * za23

    coefficient_matrix(idx_21,jlev,jnode) = &
                            & jacobian_matrix(idx_12,jlev,jnode) * za11 + &
                            & jacobian_matrix(idx_22,jlev,jnode) * za21

    coefficient_matrix(idx_22,jlev,jnode) = &
                            & jacobian_matrix(idx_12,jlev,jnode) * za12 + &
                            & jacobian_matrix(idx_22,jlev,jnode) * za22

    coefficient_matrix(idx_23,jlev,jnode) = &
                            & jacobian_matrix(idx_12,jlev,jnode) * za13 + &
                            & jacobian_matrix(idx_22,jlev,jnode) * za23

    coefficient_matrix(idx_31,jlev,jnode) = &
                            & jacobian_matrix(idx_13,jlev,jnode) * za11 + &
                            & jacobian_matrix(idx_23,jlev,jnode) * za21 + &
                            & jacobian_matrix(idx_33,jlev,jnode) * za31

    coefficient_matrix(idx_32,jlev,jnode) = &
                            & jacobian_matrix(idx_13,jlev,jnode) * za12 + &
                            & jacobian_matrix(idx_23,jlev,jnode) * za22 + &
                            & jacobian_matrix(idx_33,jlev,jnode) * za32

    coefficient_matrix(idx_33,jlev,jnode) = &
                            & jacobian_matrix(idx_13,jlev,jnode) * za13 + &
                            & jacobian_matrix(idx_23,jlev,jnode) * za23 + &
                            & jacobian_matrix(idx_33,jlev,jnode) * za33
  enddo
enddo
!$OMP END PARALLEL DO

end subroutine compute_coefficients
! =============================================================================




! =============================================================================
! VECTOR TRANSFORMATION
! =============================================================================
subroutine transform_vector(this, field_v, field_vt, mode)
class(mappings_type), intent(inout) :: this
type(atlas_Field),    intent(inout) :: field_v, field_vt
integer,              intent(in)    :: mode

real(wp)          :: zgxyi
integer           :: jnode
integer           :: jlev
real(wp), pointer :: jacobian_matrix(:,:,:)
real(wp), pointer :: zv(:,:,:)
real(wp), pointer :: zvt(:,:,:)
type(atlas_Trace) :: trace

trace=atlas_Trace(__FILENAME__,__LINE__,"mappings_type::transform_vector")
call fckit_log%debug('transform_vector')


! Retrieve jacobian_matrix from fieldset
call this%jacobian_matrix_field%data(jacobian_matrix)
call field_v %data(zv)
call field_vt%data(zvt)

! Forward transformation
if (mode == 1) then

  !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,jlev)
  do jnode = 1,this%geom%nb_nodes
    do jlev = 1,this%geom%nb_levels
      zvt(MXX,jlev,jnode) = &
            &   jacobian_matrix(idx_11,jlev,jnode) * zv(MXX,jlev,jnode) &
            & + jacobian_matrix(idx_21,jlev,jnode) * zv(MYY,jlev,jnode)

      zvt(MYY,jlev,jnode) = &
            &   jacobian_matrix(idx_12,jlev,jnode) * zv(MXX,jlev,jnode) &
            & + jacobian_matrix(idx_22,jlev,jnode) * zv(MYY,jlev,jnode)

      zvt(MZZ,jlev,jnode) = &
            &   jacobian_matrix(idx_13,jlev,jnode) * zv(MXX,jlev,jnode) &
            & + jacobian_matrix(idx_23,jlev,jnode) * zv(MYY,jlev,jnode) &
            & + jacobian_matrix(idx_33,jlev,jnode) * zv(MZZ,jlev,jnode)
    enddo
  enddo
  !$OMP END PARALLEL DO

endif

! Backward transformation
if (mode == -1) then

  !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,jlev,zgxyi)
  do jnode = 1,this%geom%nb_nodes
    do jlev = 1,this%geom%nb_levels
      zgxyi = 1.0_wp / (jacobian_matrix(idx_11,jlev,jnode) * &
                     &  jacobian_matrix(idx_22,jlev,jnode)   &
                     & -jacobian_matrix(idx_12,jlev,jnode) * &
                     &  jacobian_matrix(idx_21,jlev,jnode))

      zv(MXX,jlev,jnode) = &
        & zgxyi*(jacobian_matrix(idx_22,jlev,jnode) * zvt(MXX,jlev,jnode) &
        &       -jacobian_matrix(idx_21,jlev,jnode) * zvt(MYY,jlev,jnode))

      zv(MYY,jlev,jnode) = &
        & zgxyi*( jacobian_matrix(idx_11,jlev,jnode) * zvt(MYY,jlev,jnode) &
        &        -jacobian_matrix(idx_12,jlev,jnode) * zvt(MXX,jlev,jnode))

      zv(MZZ,jlev,jnode) = (zvt(MZZ,jlev,jnode) &
              &  -jacobian_matrix(idx_13,jlev,jnode) * zv(MXX,jlev,jnode)  &
              &  -jacobian_matrix(idx_23,jlev,jnode) * zv(MYY,jlev,jnode)) &
              & / jacobian_matrix(idx_33,jlev,jnode)
    enddo
  enddo
  !$OMP END PARALLEL DO
endif
call trace%final()

end subroutine transform_vector

! =============================================================================

subroutine final(this)
  class(mappings_type), intent(inout) :: this
  call this%jacobian_determinant%final()
  call this%jacobian_matrix_field%final()
  call this%coefficients_matrix_field%final()
end subroutine

! =============================================================================

end module 
