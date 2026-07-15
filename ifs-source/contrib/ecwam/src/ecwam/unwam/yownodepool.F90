! (C) Copyright 2001- Aron Roland (Roland & Partner, Germany).
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

!> Has data that belong to nodes
module yowNodepool
  USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
  use yowDatapool, only: rkind
  implicit none
  private
  public :: finalizeNodepool, nodes, ghosts


  !> Holds the nodes data.
  !> Such as x, y, z data or the number of connected nodes
  type, public :: t_Node

    !> the local node number
    integer(KIND=JWIM) :: id = 0

    !> the global node number
    integer(KIND=JWIM) :: id_global = 0

    !> X Coordiante
    real(rkind) :: x = 0.0
    !> Y Coordiante
    real(rkind) :: y = 0.0
    !> Z Coordiante
    real(rkind) :: z = 0.0

    !> number of connected nodes.
    !> holds the number of neighbors conntected to this node.
    !> to get the connected nodes, iterate over the connNodes() Array
    integer(KIND=JWIM) :: nConnNodes = 0

    !> The domain ID to which this node belongs.
    !> The first domain starts by 1. Fortran Stye
    integer(KIND=JWIM) :: domainID = 0

    contains
      !> Insert a node to the connected Nodes array. See Node_insertConnNode()
      !> Just a helper subroutine to make nicer code
      procedure :: insertConnNode

      !> return a pointer to the i-th node number conntected to this node
      !> Just a helper function to make nicer code
      procedure :: connNodes

      !> returns true if this node is a ghost node
      procedure :: isGhost
  end type

  !> coordinates of the local +  ghost nodes. range [1:npa]
  real(rkind), public, target, allocatable :: x(:), y(:), z(:)

  !> number of nodes, global
  integer(KIND=JWIM), public :: np_global = 0

  !> number of nodes, local
  integer(KIND=JWIM), public :: np  = 0

  !> number of ghost nodes this partition holds
  integer(KIND=JWIM), public :: ng = 0

  !> number of ghost + resident nodes this partition holds
  integer(KIND=JWIM), public :: npa = 0

  !> all nodes with their data.
  !> to iterate over all local nodes, use funtion node(local id) to get a pointer to t_node
  type(t_Node), public, allocatable, target :: nodes_global(:)


  !> max number of conntected nodes to a node
  integer(KIND=JWIM), public :: maxConnNodes = 0

  !> conntected Node Array for global nodes
  !> 2D Array. Holds the global node numbers conntected to each other.
  !> \param 1 global node number from wich you want the neighbors
  !> \param 2 from 1 to t_Node::nConnNodes
  integer(KIND=JWIM), public, allocatable :: connNodes_data(:,:)

  !> Number of connected nodes per node. [1:np]
  !> this is the same as node%nConnNodes
  integer(KIND=JWIM), public, allocatable :: NCONN(:)

  !> Holds the local node numbers connected to the given local node number
  !> \param 1 i-th neighbor of node IP. [1:NCONN(IP)]
  !> \param 2 IP local node number from which you want the neighbors [1:np] 
  !> \note the neighbor nodes are sorted. If you iterate over all connected nodes, their x/y coord. 
  !> builds a valid, simple, non-selfintersecting polygon
  integer(KIND=JWIM), public, allocatable :: CONN(:,:)



  !> Node local to global mapping.
  !> np long. give the gobal node id
  integer(KIND=JWIM), public, allocatable :: iplg(:)

  !> Node global to local mapping
  !> np_global long. give the local node id but only for this rank. local node id for other ranks are set to 0!
  integer(KIND=JWIM), public, allocatable :: ipgl(:)

  !> Ghost local to global mapping
  !> ng long. give the global node id of nodes, which
  !> belong to adjacent domains
  integer(KIND=JWIM), public, allocatable :: ghostlg(:)

  !> Ghost global to local mapping
  !> np_global long. give the local ghost node id. local ghost node ids for other ranks are set to 0!
  integer(KIND=JWIM), public, allocatable :: ghostgl(:)

  !> Numbers of Nodes pro Processor.
  !> Has the number of nodes each thread ows. Array is nTasks long
  integer(KIND=JWIM), public, allocatable :: np_perProc(:)

  !> Number of Nodes pro Processor totalize.
  !> Has the sum of nodes each thread owen. Array in nTasks+1 long
  !> Processor i stores np_perProcSum(i)::np_perProcSum(i+1)-1 nodes
  integer(KIND=JWIM), public, allocatable :: np_perProcSum(:)

  !> Node to domain mapping.
  !> np_global long. give the domain number for die global node number
  !> will be allocated and filled in subroutine runparmetis
  !> The first rank starts by 1.
  integer(KIND=JWIM), public, allocatable :: node2rank(:)

  contains

  !> return a pointer to the i-th node number conntected to this node.
  !> \param i
  !> \return pointer to the i-th node number
  function connNodes(this, i)
    implicit none
    class(t_Node) :: this
    integer(KIND=JWIM), intent(in) :: i
    type(t_Node), pointer :: connNodes
    connNodes => nodes_global(connNodes_data(this%id_global, i))
  end function

  !> return pointer to the (global) node from the local id.
  !> This is in effekt iplg(id_local)
  !> \param id_local the local node number
  !> \return poiner to the (global) node
  function nodes(id_local)
    implicit none
    integer(KIND=JWIM), intent(in) :: id_local
    type(t_Node), pointer :: nodes
    nodes => nodes_global(iplg(id_local))
  end function

  !> return pointer to the (global) (ghost) node
  !> Ghost nodes are nodes in the global node array, with the particularity
  !> that their local id is the id from another domain
  !> This is in effekt ghostlg(1:ng)
  !> \param id Counts from 1 to ng
  !> \return pointer to the (global) node
  function ghosts(id)
    implicit none
    integer(KIND=JWIM), intent(in) :: id
    type(t_Node), pointer :: ghosts
    ghosts => nodes_global(ghostlg(id))
  end function

  !> Insert a node number to the end of the conntected node array
  !> \param index optional - node number to insert. If it not present, just increas temporarily array lenght for later allocation
  subroutine insertConnNode(this, ind)
    implicit none
    class(t_Node) :: this
    integer(KIND=JWIM)    , intent(in), optional :: ind
    integer(KIND=JWIM) :: i
    type(t_Node), pointer :: node

    ! if index is present, Check if the node has allreay be insert. Then insert it
    ! if index is not present, just increas temporarily array lenght for later allocation
    if(present(ind)) then
      do i = 1, this%nConnNodes
        node => this%connNodes(i)
        if(node%id_global == ind) then
          return
        end if
      end do

      this%nConnNodes = this%nConnNodes +1
!       connNode => this%connNodes(this%nConnNodes)
      connNodes_data(this%id_global, this%nConnNodes) = ind
!       connNode = index
    else
      this%nConnNodes = this%nConnNodes +1
    end if
  end subroutine insertConnNode

  !> Returns true if this node is a ghost node
  function isGhost(this)
    implicit none
    class(t_node), intent(in) :: this
    logical :: isGhost

    if(this%id <= np) then
      isGhost = .false.
    else
      isGhost = .true.
    endif
  end function

  subroutine finalizeNodepool()
    implicit none
    
    if(allocated(x))              deallocate(x)
    if(allocated(y))              deallocate(y)
    if(allocated(z))              deallocate(z)
    if(allocated(nodes_global))   deallocate(nodes_global)
    if(allocated(connNodes_data)) deallocate(connNodes_data)
    if(allocated(iplg))           deallocate(iplg)
    if(allocated(ipgl))           deallocate(ipgl)
    if(allocated(ghostlg))        deallocate(ghostlg)
    if(allocated(ghostgl))        deallocate(ghostgl)
    if(allocated(np_perProc))     deallocate(np_perProc)
    if(allocated(np_perProcSum))  deallocate(np_perProcSum)
    if(allocated(node2rank))      deallocate(node2rank)
  end subroutine

  
end module yowNodepool
