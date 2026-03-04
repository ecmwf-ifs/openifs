! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

module advection_MPDATA_nabla_module

use atlas_module
use fckit_module, only : fckit_log
use advection_MPDATA_auxiliary_module
use advection_MPDATA_geometry_module, only : geometry_type


implicit none

private
public :: nabla_type

! Class Preconditioner
type :: nabla_type
  type(geometry_type), pointer :: geom
contains
  generic  , public :: assignment(=) => copy
  procedure, private:: copy
  procedure, public :: gradient_2d
  procedure, public :: divergence_2d
  procedure, public :: gradient_3d
  procedure, public :: divergence_3d
  procedure, public :: gradient_potential
end type

! Definition of the constructor
interface nabla_type
  module procedure nabla_constructor
end interface nabla_type


contains


! =============================================================================
! Constructor for class preconditioner
! =============================================================================
function nabla_constructor(geometry, config) result(this)
type(nabla_type) :: this

type(geometry_type), target :: geometry
type(atlas_config)     :: config

this%geom       => geometry

end function nabla_constructor
! =============================================================================

subroutine copy(obj_out,obj_in)
  class(nabla_type), intent(inout) :: obj_out
  class(nabla_type), intent(in)    :: obj_in
  obj_out%geom            => obj_in%geom
end subroutine copy


! =============================================================================
! Calculate the 3d gradient of 'inarray' without vertical scaling
! using GAUSS-GREEN strategy
! =============================================================================
subroutine gradient_2d(this, inarray, outarray, kflp)
class(nabla_type), intent(inout) :: this

! Dummy arguments
real(wp)              , intent(in)    :: inarray (:,:)
real(wp)              , intent(out)   :: outarray(:,:,:)
integer, optional     , intent(in)    :: kflp

! Local variables
type(atlas_HaloExchange)              :: halo_exchange
real(wp)                              :: zavgQS(2,this%geom%nb_levels,this%geom%nb_edges)
real(wp)                              :: zadd
real(wp)                              :: zavgQ
real(wp)                              :: zbc
integer                               :: nb_nodes
integer                               :: nb_levels
integer                               :: nb_edges
integer                               :: jnode
integer                               :: jlev
integer                               :: jedge
integer                               :: iedge
integer                               :: ip1
integer                               :: ip2
integer                               :: kflip
logical                               :: logical_test
integer, pointer                      :: inode2edge  (:)
integer                               :: inode2edge_size

! Calculate gradient using greenGauss strategy
call fckit_log%debug('nabla_gradient_2d')

! Check whether we need to apply special treatment at the poles
if (.not.present(kflp)) then
  kflip = 0
else
  kflip = kflp
endif

! Assign required values
nb_nodes     = this%geom%nb_nodes
nb_edges     = this%geom%nb_edges
nb_levels    = this%geom%nb_levels

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jedge,ip1,ip2,zbc,jlev,zavgQ)
do jedge = 1,nb_edges
  ip1 = this%geom%iedge2node(1,jedge)
  ip2 = this%geom%iedge2node(2,jedge)
  zbc = (1 - kflip) + kflip * this%geom%rpole_bc(jedge)
  do jlev = 1,nb_levels
    zavgQ = 0.5_wp * (inarray(jlev,ip1) + zbc * inarray(jlev,ip2))
    zavgQS(:,jlev,jedge) = this%geom%dual_normals(:,jedge) * zavgQ
  enddo
enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO SCHEDULE(STATIC) &
!$OMP PRIVATE(jnode,jedge,iedge,zadd,inode2edge,inode2edge_size,jlev)
do jnode = 1,nb_nodes
  outarray(:,:,jnode) = 0._wp
  call this%geom%node2edge%row(jnode, inode2edge, inode2edge_size)
  do jedge=1,inode2edge_size
    iedge = inode2edge(jedge)
    zadd  = real(this%geom%node2edge_sign(jedge,jnode),wp)
    do jlev = 1,nb_levels
      outarray(MXX,jlev,jnode) = outarray(MXX,jlev,jnode) + &
                                & zadd * zavgQS(MXX,jlev,iedge)
      outarray(MYY,jlev,jnode) = outarray(MYY,jlev,jnode) + &
                                & zadd * zavgQS(MYY,jlev,iedge)
    enddo
  enddo
enddo
!$OMP END PARALLEL DO

! Special treatment for the north & south pole cell faces
! Sx == 0 at pole, and Sy has same sign at both sides of pole
if (.not. present(kflp)) then
do jedge = 1,this%geom%nb_pole_edges
  iedge = this%geom%pole_edges(jedge)
  ip1 = this%geom%iedge2node(1,iedge)
  ip2 = this%geom%iedge2node(2,iedge)
  do jlev = 1,nb_levels
    ! Correct for wrong Y-derivatives in previous loop
    outarray(MYY,jlev,ip2) = outarray(MYY,jlev,ip2) + &
                           & 2._wp * zavgQS(MYY,jlev,iedge)
  enddo
enddo
else
  ! do nothing
endif

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,jlev)
do jnode = 1,nb_nodes
  do jlev = 1,nb_levels
    outarray(MXX,jlev,jnode) = outarray(MXX,jlev,jnode) / this%geom%dual_volumes(jnode)
    outarray(MYY,jlev,jnode) = outarray(MYY,jlev,jnode) / this%geom%dual_volumes(jnode)
  enddo
enddo
!$OMP END PARALLEL DO

! Get halo_exchange
halo_exchange = this%geom%nodes%get_halo_exchange()
call halo_exchange%execute(outarray)

end subroutine gradient_2d
! =============================================================================



! =============================================================================
! Calculate the divergence of 'rho V' without vertical scaling
! using GAUSS-GREEN strategy
! =============================================================================
subroutine divergence_2d(this, prho, pVstar, pdiv, ldhaloex)
class(nabla_type), intent(inout) :: this

! Dummy arguments
real(wp)              , intent(in)    :: pVstar(:,:,:)
real(wp)              , intent(in)    :: prho  (:,:)
real(wp)              , intent(out)   :: pdiv  (:,:)
logical, optional     , intent(in)    :: ldhaloex

! Local variables
type(atlas_HaloExchange)              :: halo_exchange
real(wp)                              :: zavs(2,this%geom%nb_levels,this%geom%nb_edges)
real(wp)                              :: zaun(this%geom%nb_levels,this%geom%nb_edges)
real(wp)                              :: zadd
integer                               :: nb_nodes
integer                               :: nb_levels
integer                               :: nb_edges
integer                               :: jnode
integer                               :: jlev
integer                               :: jedge
integer                               :: iedge
integer                               :: ip1
integer                               :: ip2
logical                               :: llhaloex
logical                               :: logical_test
integer, pointer                      :: inode2edge  (:)
integer                               :: inode2edge_size


call fckit_log%debug('nabla_div_green_gauss_xy')

! Check whether we need to apply halo_exchange at the end of this routine
if (present(ldhaloex)) then
  llhaloex = ldhaloex
else
  llhaloex = .true.
endif

! Assign required values
nb_levels           = this%geom%nb_levels
nb_nodes            = this%geom%nb_nodes
nb_edges            = this%geom%nb_edges

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jedge,ip1,ip2,jlev)
do jedge = 1,nb_edges
  ip1 = this%geom%iedge2node(1,jedge)
  ip2 = this%geom%iedge2node(2,jedge)
  do jlev = 1,nb_levels
    zavs(MXX,jlev,jedge) = (pVstar(MXX,jlev,ip1) * prho(jlev,ip1) + &
                         & this%geom%rpole_bc(jedge) * pVstar(MXX,jlev,ip2) * &
                         & prho(jlev,ip2)) * 0.5_wp

    zavs(MYY,jlev,jedge) = (pVstar(MYY,jlev,ip1) * prho(jlev,ip1) + &
                         & this%geom%rpole_bc(jedge) * pVstar(MYY,jlev,ip2) * &
                         & prho(jlev,ip2)) * 0.5_wp
  enddo
enddo
!$OMP END PARALLEL DO

do jedge = 1,this%geom%nb_pole_edges
  iedge = this%geom%pole_edges(jedge)
  do jlev = 1,nb_levels
    zavs(MYY,jlev,iedge) = 0._wp
  enddo
enddo

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jedge,jlev)
do jedge = 1,nb_edges
  do jlev = 1,nb_levels
    zaun(jlev,jedge) = zavs(MXX,jlev,jedge) * this%geom%dual_normals(MXX,jedge) + &
                     & zavs(MYY,jlev,jedge) * this%geom%dual_normals(MYY,jedge)
  enddo
enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,inode2edge,inode2edge_size,jlev,jedge,iedge,zadd)
do jnode = 1,nb_nodes
  call this%geom%node2edge%row(jnode, inode2edge, inode2edge_size)
  do jlev = 1,nb_levels
    pdiv(jlev,jnode) = 0._wp
    do jedge=1,inode2edge_size
      iedge = inode2edge(jedge)
      zadd  = real(this%geom%node2edge_sign(jedge,jnode),wp)
      pdiv(jlev,jnode) = pdiv(jlev,jnode) + zadd * zaun(jlev,iedge)
    enddo
    pdiv(jlev,jnode) = pdiv(jlev,jnode) / (this%geom%dual_volumes(jnode))
  enddo
enddo
!$OMP END PARALLEL DO

! Get halo_exchange
if (llhaloex) then
  halo_exchange = this%geom%nodes%get_halo_exchange()
  call halo_exchange%execute(pdiv)
endif
end subroutine divergence_2d
! =============================================================================



! =============================================================================
! Calculate the 3d gradient of 'inarray' using GAUSS-GREEN strategy
! =============================================================================
subroutine gradient_3d(this, inarray, outarray, kflp)
class(nabla_type), intent(inout) :: this

! Dummy arguments
real(wp)                       , intent(in)  :: inarray (:,:)
real(wp)                       , intent(out) :: outarray(:,:,:)
integer, optional              , intent(in)  :: kflp

! Local variables
type(atlas_HaloExchange)              :: halo_exchange
real(wp)                              :: zavgQS(2,this%geom%nb_levels,this%geom%nb_edges)
real(wp)                              :: zadd
real(wp)                              :: zavgQ
real(wp)                              :: zbc
real(wp)                              :: zdzi
real(wp)                              :: zdzi2
integer                               :: nb_nodes
integer                               :: nb_levels
integer                               :: nb_edges
integer                               :: jnode
integer                               :: jlev
integer                               :: jedge
integer                               :: iedge
integer                               :: ip1
integer                               :: ip2
integer                               :: kflip
logical                               :: logical_test
integer, pointer                   :: inode2edge  (:)
integer                            :: inode2edge_size

call fckit_log%debug('nabla_grad_green_gauss_xy_3d')

! Check whether we need to apply special treatment at the poles
if (.not.present(kflp)) then
  kflip = 0
else
  kflip = kflp
endif

! Assign required values
nb_nodes  = this%geom%nb_nodes
nb_levels = this%geom%nb_levels
nb_edges  = this%geom%nb_edges

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jedge,ip1,ip2,zbc,jlev,zavgQ)
do jedge = 1,nb_edges
  ip1 = this%geom%iedge2node(1,jedge)
  ip2 = this%geom%iedge2node(2,jedge)
  zbc = (1 - kflip) + kflip * this%geom%rpole_bc(jedge)
  do jlev = 1,nb_levels
    zavgQ = 0.5_wp * (inarray(jlev,ip1) + zbc * inarray(jlev,ip2))
    zavgQS(MXX,jlev,jedge) = this%geom%dual_normals(MXX,jedge) * zavgQ
    zavgQS(MYY,jlev,jedge) = this%geom%dual_normals(MYY,jedge) * zavgQ
  enddo
enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,inode2edge,inode2edge_size,jedge,iedge,zadd,jlev)
do jnode = 1,nb_nodes
  outarray(:,:,jnode) = 0._wp
  call this%geom%node2edge%row(jnode, inode2edge, inode2edge_size)
  do jedge=1,inode2edge_size
    iedge = inode2edge(jedge)
    zadd  = real(this%geom%node2edge_sign(jedge,jnode),wp)
    do jlev = 1,nb_levels
      outarray(MXX,jlev,jnode) = outarray(MXX,jlev,jnode) + &
                                & zadd * zavgQS(MXX,jlev,iedge)
      outarray(MYY,jlev,jnode) = outarray(MYY,jlev,jnode) + &
                                & zadd * zavgQS(MYY,jlev,iedge)
    enddo
  enddo
enddo
!$OMP END PARALLEL DO


! Special treatment for the north & south pole cell faces
! Sx == 0 at pole, and Sy has same sign at both sides of pole
if (.not. present(kflp)) then
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jedge,iedge,ip1,ip2,jlev)
  do jedge = 1,this%geom%nb_pole_edges
    iedge = this%geom%pole_edges(jedge)
    ip1 = this%geom%iedge2node(1,iedge)
    ip2 = this%geom%iedge2node(2,iedge)
    do jlev = 1,nb_levels
      ! Correct for wrong Y-derivatives in previous loop
      outarray(MYY,jlev,ip2) = outarray(MYY,jlev,ip2) + &
                             & 2._wp * zavgQS(MYY,jlev,iedge)
    enddo
  enddo
!$OMP END PARALLEL DO
else
  ! do nothing
endif

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,jlev)
do jnode = 1,nb_nodes
  do jlev = 1,nb_levels
    outarray(MXX,jlev,jnode) = outarray(MXX,jlev,jnode) / this%geom%dual_volumes(jnode)
    outarray(MYY,jlev,jnode) = outarray(MYY,jlev,jnode) / this%geom%dual_volumes(jnode)
  enddo
enddo
!$OMP END PARALLEL DO

! Get halo_exchange
halo_exchange = this%geom%nodes%get_halo_exchange()
call halo_exchange%execute(outarray)

! Vertical part for every node
zdzi  = 1._wp / this%geom%dz
zdzi2 = 0.5_wp * zdzi

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode)
do jnode = 1,nb_nodes
  outarray(MZZ,2:nb_levels-1,jnode) = (inarray(3:nb_levels,jnode) - &
                               & inarray(1:nb_levels-2,jnode)) * zdzi2

  outarray(MZZ,1:nb_levels:nb_levels-1,jnode) = &
  & (inarray(2:nb_levels:nb_levels-2  ,jnode) - &
  &  inarray(1:nb_levels-1:nb_levels-2,jnode)) * zdzi
enddo
!$OMP END PARALLEL DO

end subroutine gradient_3d
! =============================================================================



! =============================================================================
! Calculate the 3d divergence of 'rho V' using GAUSS-GREEN strategy
! =============================================================================
subroutine divergence_3d(&
& this, prho, ppVstar, pdiv, psign, ldhaloex)
class(nabla_type), intent(inout) :: this

! Dummy arguments
real(wp)                       , intent(in)    :: prho   (:,:)
real(wp)                       , intent(in)    :: ppVstar(:,:,:)
real(wp)                       , intent(out)   :: pdiv   (:,:)
real(wp)                       , intent(in)    :: psign
logical, optional              , intent(in)    :: ldhaloex

! Local variables
type(atlas_HaloExchange)              :: halo_exchange
real(wp)                              :: zavs(2,this%geom%nb_levels,this%geom%nb_edges)
real(wp)                              :: zaun(this%geom%nb_levels,this%geom%nb_edges)
real(wp)                              :: zadd
real(wp)                              :: za1
real(wp)                              :: za2
real(wp)                              :: zavgQ
real(wp)                              :: zbc
real(wp)                              :: zdzi
real(wp)                              :: zdzi2
integer                               :: nb_nodes
integer                               :: nb_levels
integer                               :: nb_edges
integer                               :: jnode
integer                               :: jlev
integer                               :: jedge
integer                               :: iedge
integer                               :: ip1
integer                               :: ip2
logical                               :: llhaloex
logical                               :: logical_test
integer, pointer                      :: inode2edge  (:)
integer, pointer                      :: iedge2node  (:,:)
integer                               :: inode2edge_size
integer                               :: nb_pole_edges


call fckit_log%debug('nabla_divergence_3d')

! Check whether we need to apply halo_exchange at the end of this routine
if (present(ldhaloex)) then
  llhaloex = ldhaloex
else
  llhaloex = .true.
endif

! Assign required values
nb_levels           = this%geom%nb_levels
nb_nodes            = this%geom%nb_nodes
nb_edges            = this%geom%nb_nodes

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jedge,ip1,ip2,jlev)
do jedge  = 1,nb_edges
  ip1 = iedge2node(1,jedge)
  ip2 = iedge2node(2,jedge)
  do jlev = 1,nb_levels
    zavs(MXX,jlev,jedge) = (ppVstar(MXX,jlev,ip1) * prho(jlev,ip1) + &
                         & this%geom%rpole_bc(jedge) * ppVstar(MXX,jlev,ip2) * &
                         & prho(jlev,ip2)) * 0.5_wp

    zavs(MYY,jlev,jedge) = (ppVstar(MYY,jlev,ip1) * prho(jlev,ip1) + &
                         & this%geom%rpole_bc(jedge) * ppVstar(MYY,jlev,ip2) * &
                         & prho(jlev,ip2)) * 0.5_wp
  enddo
enddo
!$OMP END PARALLEL DO

! Boundary conditions
do jedge = 1,this%geom%nb_pole_edges
  iedge = this%geom%pole_edges(jedge)
  do jlev = 1,nb_levels
    zavs(MYY,jlev,iedge) = 0._wp
  enddo
enddo

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jedge,jlev)
do jedge = 1,nb_edges
  do jlev = 1,nb_levels
    zaun(jlev,jedge) =  zavs(MXX,jlev,jedge)*this%geom%dual_normals(MXX,jedge) &
                     & +zavs(MYY,jlev,jedge)*this%geom%dual_normals(MYY,jedge)
  enddo
enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,inode2edge,inode2edge_size,jlev,jedge,iedge,zadd)
do jnode = 1,nb_nodes
  call this%geom%node2edge%row(jnode,inode2edge,inode2edge_size)
  do jlev = 1,nb_levels
    pdiv(jlev,jnode) = 0._wp
    do jedge=1,inode2edge_size
      iedge = inode2edge(jedge)
      zadd  = real(this%geom%node2edge_sign(jedge,jnode),wp)
      pdiv(jlev,jnode) = pdiv(jlev,jnode) + zadd * zaun(jlev,iedge)
    enddo
    pdiv(jlev,jnode) = pdiv(jlev,jnode) / (this%geom%dual_volumes(jnode))
  enddo
enddo
!$OMP END PARALLEL DO

! Get halo_exchange
if (llhaloex) then
  halo_exchange = this%geom%nodes%get_halo_exchange()
  call halo_exchange%execute(pdiv)
endif

! Vertical part for every node
zdzi   = 1._wp / this%geom%dz
zdzi2 = 0.5_wp * zdzi

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode)
do jnode = 1,nb_nodes
  pdiv(2:nb_levels-1,jnode) = &
    & pdiv(2:nb_levels-1,jnode) + (ppVstar(MZZ,3:nb_levels,  jnode) * &
    & prho(3:nb_levels  ,jnode) -  ppVstar(MZZ,1:nb_levels-2,jnode) * &
    & prho(1:nb_levels-2,jnode)) * zdzi2

  pdiv(1:nb_levels:nb_levels-1,jnode) = &
    &  pdiv   (    1:nb_levels:nb_levels-1  ,jnode)  + &
    & (ppVstar(MZZ,2:nb_levels:nb_levels-2  ,jnode)  * &
    &  prho   (    2:nb_levels:nb_levels-2  ,jnode)  - &
    &  ppVstar(MZZ,1:nb_levels-1:nb_levels-2,jnode)  * &
    &  prho   (    1:nb_levels-1:nb_levels-2,jnode)) * zdzi
enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,jlev)
do jnode = 1,nb_nodes
  do jlev = 1,nb_levels
    pdiv(jlev,jnode) = psign * pdiv(jlev,jnode) / prho(jlev,jnode)
  enddo
enddo
!$OMP END PARALLEL DO

end subroutine divergence_3d
! =============================================================================



! =============================================================================
! Compute (V - \nabla \phi)
! =============================================================================
subroutine gradient_potential(&
& this, mappings_fields, pp, pVhat, pVstar, pfd)
class(nabla_type), intent(inout) :: this

! Dummy arguments
type(atlas_fieldset), intent(in)    :: mappings_fields
real(wp)            , intent(inout) :: pp    (:,:)
real(wp)            , intent(in)    :: pVhat (:,:,:)
real(wp)            , intent(out)   :: pVstar(:,:,:)
real(wp)            , intent(in)    :: pfd   (:,:)

! Local variables
type(atlas_field)                     :: coefficients_matrix_field
real(wp), pointer                     :: coefficients_matrix(:,:,:,:)
real(wp)                              :: zgradpp(3,this%geom%nb_levels,this%geom%nb_nodes)

integer :: nb_nodes
integer :: nb_levels
integer :: jnode
integer :: jlev
integer :: inl
integer :: inl1
logical :: logical_test

call fckit_log%debug('gradient_potential')

! Assign required values
nb_levels           = this%geom%nb_levels
nb_nodes            = this%geom%nb_nodes

! Retrieve coefficients_matrix from fieldset
coefficients_matrix_field = mappings_fields%field("coefficients_matrix")
call coefficients_matrix_field%data(coefficients_matrix)

call this%gradient_3d(pp, zgradpp)

! Redefine number of vertical levels for loop
inl  = nb_levels
inl1 = inl-1

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,jlev)
do jnode  = 1,nb_nodes
  do jlev = 2,inl1
    pVstar(:,jlev,jnode) = pVhat(:,jlev,jnode) -           &
                         & matmul(coefficients_matrix(:,:,jlev,jnode), &
                         &        zgradpp(:,jlev,jnode))
  enddo
enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode)
do jnode = 1,nb_nodes
  pVstar(MZZ,1:inl:inl1,jnode) = 0._wp

  zgradpp(MZZ,1:inl:inl1,jnode) = (pVhat(MZZ,1:inl:inl1,jnode) -   &
                                & pVstar(MZZ,1:inl:inl1,jnode) -   &
                                & coefficients_matrix(3,1,1:inl:inl1,jnode) * &
                                & zgradpp(MXX,1:inl:inl1,jnode) -  &
                                & coefficients_matrix(3,2,1:inl:inl1,jnode) * &
                                & zgradpp(MYY,1:inl:inl1,jnode)) / &
                                & coefficients_matrix(3,3,1:inl:inl1,jnode)

  pVstar(MXX:MYY, 1,jnode) = pVhat(MXX:MYY,1,jnode) -         &
                           & matmul(coefficients_matrix(1:2,1:3,1,jnode), &
                           &        zgradpp(:,1,jnode))

  pVstar(MXX:MYY,inl,jnode) = pVhat(MXX:MYY,inl,jnode) -         &
                            & matmul(coefficients_matrix(1:2,1:3,inl,jnode), &
                            &        zgradpp(:,inl,jnode))
enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,jlev)
do jnode = 1,nb_nodes
  do jlev = 1,nb_levels
    pVstar(:,jlev,jnode) = pfd(jlev,jnode) * pVstar(:,jlev,jnode)
  enddo
enddo
!$OMP END PARALLEL DO

end subroutine gradient_potential
! =============================================================================

end module advection_MPDATA_nabla_module


