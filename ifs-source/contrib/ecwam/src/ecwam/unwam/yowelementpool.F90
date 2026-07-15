! (C) Copyright 2001- Aron Roland (Roland & Partner, Germany).
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

module yowElementpool
  USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
  implicit none
  private
  public :: finalizeElementpool, eleBelongTo

  !> number of elements, global
  integer(KIND=JWIM), public :: ne_global = 0

  !> number of local elements
  integer(KIND=JWIM), public :: ne = 0

  !> local element array. it stores the local node IDs
  !> first index from 1 to 3.
  !> second index from 1 to ne.
  !> local node IDs in [1:np]. local ghost IDs in [np+1:np+ng]
  integer(KIND=JWIM), public, target, allocatable :: INE(:,:)

  !> global element array. it stored the global node IDs
  !> first index from 1 to 3.
  !> second index from 1 to ne_global
  integer(KIND=JWIM), public, allocatable :: INE_global(:,:)

  !> Element local to global mapping
  !> ne long. give the global element id
  integer(KIND=JWIM), public, target, allocatable :: ielg(:)

  !> Element global to local mapping
  !> ne_global long. give the local element id but only for this rank. local element id for other ranks are set to 0!
  integer(KIND=JWIM), public, allocatable :: iegl(:)

  contains

  !> Returns the domainID to which this element belongs
  !> conversione: if a element has two nodes from domain 1 and one node from domain 2, the element belongs to domain 1.
  !> If a element adjoint to three different domains, it belongs to the with the lowest domain ID
  !> The one node from domain 2 is, of course, a ghost node.
  !> @param[in] ele_global three global node IDs e.g. INE(:,IE)
  function eleGetDomainID(ele_global) result(domainID)
    use yowDatapool, only: nTasks
    use yowNodepool, only: t_Node, nodes_global
    implicit none
    integer(KIND=JWIM), intent(in) :: ele_global(3)
    integer(KIND=JWIM) :: domainID

    integer(KIND=JWIM) :: j, itemp, ranks
    type(t_Node) ::  nodes(3)

    domainID = -1

    ! check if this element adjoint to three different domains.
    nodes(:) = nodes_global(ele_global)
    if(nodes(1)%domainID /= nodes(2)%domainID .and. &
    & nodes(1)%domainID /= nodes(3)%domainID .and. &
    & nodes(2)%domainID /= nodes(3)%domainID) then
      domainID = minval(nodes(:)%domainID)

    ! check if this element has two nodes wich belongs to this domain
    else
      do ranks = 0, nTasks-1
        itemp = 0
        do j=1, 3
          if(nodes(j)%domainID == ranks+1) then
            itemp = itemp + 1
          endif
        end do
        ! yes, this element belongs to rank
        if(itemp >= 2) then
          domainID = ranks+1
          exit
        endif
      end do
    endif
  end function

  !> Returns true if the element belongs to rank.
  !> conversione: If a element is connected to domain 1,2 and 3. It belongs to 1,2 and 3.
  !> @param[in] ele_global three global node IDs e.g. INE(:,IE)
  !> @param[in] rank optional. If not given, datapool:myrank is used
  function eleBelongTo(ele_global, rank) result(belongTo)
    use yowDatapool, only: myrank
    use yowNodepool, only: t_Node, nodes_global
    implicit none
    integer(KIND=JWIM), intent(in) :: ele_global(3)
    integer(KIND=JWIM), intent(in), optional :: rank
    logical :: belongTo

    integer(KIND=JWIM) :: myDomainID
    type(t_Node) ::  nodes(3)

    belongTo = .false.

    if(present(rank) .eqv. .true.) then
      myDomainID = rank +1
    else
      myDomainID = myrank + 1
    endif

    ! check if this element adjoint to three different domains.
    nodes(:) = nodes_global(ele_global)
    if(nodes(1)%domainID == myDomainID .or. &
    & nodes(2)%domainID == myDomainID .or. &
    & nodes(3)%domainID == myDomainID) then
      belongTo = .true.
    endif
  end function

  subroutine finalizeElementpool()
    implicit none

    if(allocated(INE))        deallocate(INE)
    if(allocated(INE_global)) deallocate(INE_global)
    if(allocated(ielg))       deallocate(ielg)
    if(allocated(iegl))       deallocate(iegl)
  end subroutine
end module yowElementpool
