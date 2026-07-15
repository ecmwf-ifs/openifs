! (C) Copyright 2001- Aron Roland (Roland & Partner, Germany).
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

!> \file yowpdlibmain.F90
!> \brief initialization
!> \author Thomas Huxhorn
!> \date 2011-2012

#include "yowincludes.h"

module yowpdlibMain
  USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
  use yowError
  use yowDatapool, only: rkind
  implicit none
  private
  public :: initPD, finalizePD, initFromGridDim, clearGlobalMeshVars

  interface initPD
    module procedure initPD0, initPD1, initPD2, initFromGrid, initFromGridDim
  end interface

  contains

  !> only initialize MPI
  subroutine initPD0()
    use yowMpiModule
    implicit none

    call initMPI(MPI_COMM_WORLD)
  end subroutine

  
  !> initialize all
  !> init MPI, readin the mesh, fill the datastructure ...
  !> @param[in] filename mesh e.g. system.dat
  !> @note call exchangeModule::setDimSize() before initPD to set the dimensions size for the exchange routine
  !> @overload
  subroutine initPD1(filename)
    use yowMpiModule
    implicit none
    character(len=*), intent(in) :: filename

    call initPD(filename, MPI_COMM_WORLD)
  end subroutine initPD1
  
  
  !> @overload initPD1
  !> @param[in] MPIComm MPI communicator to use with pdlib
  subroutine initPD2(filename, MPIcomm)
    use yowMPIModule
    implicit none
    character(len=*), intent(in) :: filename
#ifdef WAM_HAVE_MPI_F08
    type(mpi_comm), intent(in) :: MPIcomm
#else
    integer(KIND=JWIM), intent(in) :: MPIcomm
#endif

    integer(KIND=JWIM) :: MNP, MNE
    integer(KIND=JWIM), allocatable :: INE(:,:)
    real(kind=rkind), allocatable :: XP(:), YP(:), DEP(:)

    call readMesh(filename, MNP, XP, YP, DEP, MNE, INE)
    call initFromGrid(MNP, XP, YP, DEP, MNE, INE, MPIcomm)

    if(allocated(XP)) deallocate(XP)
    if(allocated(YP)) deallocate(YP)
    if(allocated(DEP)) deallocate(DEP)
    if(allocated(INE)) deallocate(INE)

    
  end subroutine initPD2
  
  !> @param[in] MNP number of nodes global
  !> @param[in] XP node X value
  !> @param[in] XY node Y value
  !> @param[in] DEP node Z value
  !> @param[in] MNE number of element global
  !> @param[in] INE element array
  !> @param[in] MPIComm MPI communicator to use with pdlib
  !> @overload initPD1
  !> alter: np_global, nodes_global(), ne_global, elements(), INE_global
  subroutine initFromGrid(MNP, XP, YP, DEP, MNE, INE, MPIcomm)
    use yowExchangeModule, only : n2ndDim, n3ndDim
    use yowMpiModule
    implicit none
    integer(KIND=JWIM), intent(in) :: MNP, MNE
    integer(KIND=JWIM), intent(in) :: INE(3,MNE)
    real(kind=rkind), intent(in) :: XP(MNP), YP(MNP), DEP(MNP)
#ifdef WAM_HAVE_MPI_F08
    type(mpi_comm), intent(in) :: MPIcomm
#else
    integer(KIND=JWIM), intent(in) :: MPIcomm
#endif

    call initFromGridDim(MNP, XP, YP, DEP, MNE, INE, n2ndDim, n3ndDim, MPIcomm)
  end subroutine

  !> @param[in] MNP number of nodes global
  !> @param[in] XP node X value
  !> @param[in] XY node Y value
  !> @param[in] DEP node Z value
  !> @param[in] MNE number of element global
  !> @param[in] INE element array
  !> @param[in] secDim size of the second dimensions to exchange
  !> @param[in] thirdDim size of the third dimensions to exchange
  !> @param[in] MPIComm MPI communicator to use with pdlib
  !> @overload initPD1
  subroutine initFromGridDim(MNP, XP, YP, DEP, MNE, INE, secDim, thirdDim, MPIcomm)
      use yowDatapool,       only: myrank, debugPrePartition, debugPostPartition
      use yowNodepool,       only: np_global, np, np_perProcSum, ng
      use yowElementpool,    only: ne_global,ne
      use yowSidepool,       only: ns, ns_global
      use yowExchangeModule, only: nConnDomains, setDimSize
      use yowRankModule,     only: initRankModule
      use yowMpiModule
    implicit none
    integer(KIND=JWIM), intent(in) :: MNP, MNE
    integer(KIND=JWIM), intent(in) :: INE(3,MNE)
    real(kind=rkind), intent(in) :: XP(MNP), YP(MNP), DEP(MNP)
    integer(KIND=JWIM), intent(in) :: secDim, thirdDim
#ifdef WAM_HAVE_MPI_F08
    type(mpi_comm), intent(in) :: MPIcomm
#else
    integer(KIND=JWIM), intent(in) :: MPIcomm
#endif

    call setDimSize(secDim, thirdDim)
    call initMPI(MPIcomm)
    call assignMesh(MNP, XP, YP, DEP, MNE, INE)
!     pause
    call prePartition()
    call findConnNodes()

    if(debugPrePartition) then
      if(myrank == 0) then
        write(*,*) "pre-partition"
        write(*,*) "# Nodes: ", np_global
        write(*,*) "# Elements: ", ne_global
        write(*,*) "# Sides ", ns_global
        write(*,*) "np_perProcSum :", np_perProcSum
      end if
      write(*,*) "Thread", myrank, "# local nodes ", np
      write(*,*) "Thread", myrank, "# local sides" , ns
    endif

    call runParmetis()
    call postPartition()
    call findGhostNodes()
    call findConnDomains()
    call exchangeGhostIds()
    call postPartition2()
    call initRankModule()

    if(debugPostPartition) then
      if(myrank == 0) then
        write(*,*) "New data after partition"
        write(*,*) "# Nodes: ", np_global
        write(*,*) "# Elements: ", ne_global
        write(*,*) "# Sides ", ns_global
        write(*,*) "np_perProcSum :", np_perProcSum
      end if
      write(*,*) "Thread", myrank, "# local elements ", ne
      write(*,*) "Thread", myrank, "# local nodes ", np
      write(*,*) "Thread", myrank, "# local sides" , ns
      write(*,*) "Thread", myrank, "# of ghosts", ng
      write(*,*) "Thread", myrank, "# of neighbor domains", nConnDomains
    endif
  end subroutine initFromGridDim




  !-------------------------------------------------------------------------------
  ! Init MPI
  !-------------------------------------------------------------------------------

  !> initialize MPI.
  subroutine initMPI(MPIcomm)
    use yowDatapool, only: comm, nTasks, myrank
    use yowError
    use yowMpiModule
    implicit none
#ifdef WAM_HAVE_MPI_F08
    type(mpi_comm), intent(in) :: MPIcomm
#else
  integer(KIND=JWIM), intent(in) :: MPIcomm
#endif
    logical :: flag
    integer(KIND=JWIM) :: ierr

    if(MPIcomm == MPI_COMM_NULL) then
      ABORT("A null communicator is not allowed")
    endif

    comm = MPIcomm
    call MPI_Initialized(flag, ierr)
    if(ierr/=MPI_SUCCESS) call parallel_abort(error=ierr)

    if(flag .eqv. .false.) then
      call mpi_init(ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort(error=ierr)
    endif

    ! Get number of processors
    call mpi_comm_size(comm, nTasks,ierr)
    if(ierr/=MPI_SUCCESS) call parallel_abort(error=ierr)

    ! Get rank
    call mpi_comm_rank(comm, myrank,ierr)
    if(ierr/=MPI_SUCCESS) call parallel_abort(error=ierr)
  end subroutine initMPI

  !> read the mesh
  !> @internal
  !> @param[inout] MNP number of nodes global
  !> @param[inout] XP node X value
  !> @param[inout] XY node Y value
  !> @param[inout] DEP node Z value
  !> @param[inout] MNE number of element global
  !> @param[inout] INE element array
  !> @note the caller in a responsible to deallocate XP,YP,DEP, INE
  subroutine readMesh(filename, MNP, XP, YP, DEP, MNE, INE)
    use yowError,  only: parallel_abort
#ifdef ECMWF
    use yowunpool, only: grid
#endif
    implicit none
    character(len=*), intent(in) :: filename
    integer(KIND=JWIM), intent(inout) :: MNP, MNE
    integer(KIND=JWIM), allocatable, intent(inout) :: INE(:,:)
    real(kind=rkind), allocatable, intent(inout) :: XP(:), YP(:), DEP(:)

    integer(KIND=JWIM) :: itemp
    integer(KIND=JWIM) :: i, fd

#ifdef ECMWF
    fd = GRID%FHNDL!newunit()
#else
    fd = newunit()
#endif
    open(fd, file=filename, status='old')
!     if(stat/=0) call parallel_abort(' open system.dat failure')

    ! Dummy read
    read(fd,*);read(fd,*);
    ! # boundary Nodes
    read(fd,*) itemp
    read(fd,*)
    ! # resistend nodes
    read(fd,*) MNP

    MNP = MNP + itemp

    if(allocated(XP)) deallocate(XP)
    allocate(XP(MNP), stat=stat);
    if(stat/=0) ABORT('XP allocate failure')

    if(allocated(YP)) deallocate(YP)
    allocate(YP(MNP), stat=stat);
    if(stat/=0) ABORT('YP allocate failure')

    if(allocated(DEP)) deallocate(DEP)
    allocate(DEP(MNP), stat=stat);
    if(stat/=0) ABORT('DEP allocate failure')

    read(fd,*);read(fd,*);read(fd,*);read(fd,*);read(fd,*);read(fd,*);read(fd,*);

    ! Read all nodes
    do i=1, MNP
      read(fd,*) itemp, XP(i), YP(i), DEP(i)
    end do

    read(fd,*);read(fd,*);
    ! # elements
    read(fd,*) MNE

    if(allocated(INE)) deallocate(INE)
    allocate(INE(3, MNE), stat=stat);
    if(stat/=0) ABORT('INE allocate failure')

    read(fd,*);read(fd,*);read(fd,*);

    ! Read all elements
    do i = 1, MNE
      read(fd,*) INE(1,i), INE(2,i), INE(3,i), itemp, itemp
      ! fortran counts from 1
      INE(:,i) = INE(:,i) + 1
    end do

    close(fd)
  end subroutine readMesh


  !> @param[in] MNP number of nodes global
  !> @param[in] XP node X value
  !> @param[in] XY node Y value
  !> @param[in] DEP node Z value
  !> @param[in] MNE number of element global
  !> @param[in] INE element array
  !> alter: np_global, nodes_global(), ne_global, INE_global
  subroutine assignMesh(MNP, XP, YP, DEP, MNE, INE)
    use yowNodepool,    only: nodes_global, np_global, t_Node
    use yowElementpool, only: ne_global, INE_global
    implicit none
    integer(KIND=JWIM), intent(in) :: MNP, MNE
    integer(KIND=JWIM), intent(in) :: INE(3,MNE)
    real(kind=rkind), intent(in) :: XP(MNP), YP(MNP), DEP(MNP)

    integer(KIND=JWIM) :: i
    logical :: nodeIDexist(MNP)

    character(len=1000) :: errstr

    ne_global=MNE
    np_global=MNP
    
    ! various checks

    if(ne_global < 1) then
      ABORT("assignMesh() ne_global < 1")
    endif

    if(np_global < 1) then
      ABORT("assignMesh() np_global < 1")
    endif
    
    if(minval(INE(:,1:ne_global)) < 1) then
      ABORT("assignMesh() minval INE < 1")
    endif

    if(maxval(INE(:,1:ne_global)) > np_global) then
      ABORT("assignMesh() maxval INE > np_global")
    endif

    ! check if any global node ID exist at least one time in INE
    nodeIDexist(:) = .false.
    do i=1, ne_global
      nodeIDexist(INE(:,i)) = .true.
    end do

!!!debile
    do i=1, MNP 
      if( .not. nodeIDexist(i) ) then
        write(*,*) 'debile in yowpdlibmain, node not in any element : ',i
      endif
    end do


    if(any(nodeIDexist(:) .eqv. .false.) .eqv. .true.) then
      ABORT("assignMesh() INE does not contains all node IDs")
    endif

    
    if(allocated(nodes_global)) deallocate(nodes_global)
    allocate(nodes_global(np_global), stat=stat);
    if(stat/=0) then
      write(errstr,*) "nodes_global allocate failure", np_global*B2MB, "MB"
      ABORT(errstr)
    endif

    do i=1, np_global
      nodes_global(i)%x = XP(I)
      nodes_global(i)%y = YP(I)
      nodes_global(i)%z = DEP(I)
      nodes_global(i)%id_global = i
    end do    

    if(allocated(INE_global)) deallocate(INE_global)
    allocate(INE_global(3, ne_global), stat=stat);
    if(stat/=0) then
      write(errstr,*) "INE_global allocate failure", 4*3*ne_global*B2MB, "MB"
      ABORT(errstr)
    endif

    do i = 1, ne_global
      INE_global(:, i) = INE(:,i)
    end do
  end subroutine assignMesh


  !-------------------------------------------------------------------------------
  ! pre-partition: divide the mesh into nTasks parts.
  !-------------------------------------------------------------------------------

  !> pre-partition the mesh
  !> just divide the mesh into nTasks parts
  !> and create a premature iplg
  !> alter: np_perProc, np_perProcSum, np, iplg
  subroutine prePartition
    use yowDatapool, only: nTasks ,myrank
    use yowNodepool, only: np_global, np, np_perProc, np_perProcSum, iplg
    implicit none

    integer(KIND=JWIM) :: i

    ! determine equal number of nodes in each processor (except for the last one).
    ! and create a provisional node local to global mapping iplg

    ! start the arrays from 0, because the first thread id is 0
    if(allocated(np_perProc)) deallocate(np_perProc)
    allocate(np_perProc(0:nTasks-1), stat=stat)
    if(stat/=0) call parallel_abort('np_perProc allocation failure')

    if(allocated(np_perProcSum)) deallocate(np_perProcSum)
    allocate(np_perProcSum(0:nTasks), stat=stat)
    if(stat/=0) call parallel_abort('np_perProcSum allocation failure')

    np_perProcSum = 0
    np = np_global / nTasks

    do i = 0, nTasks-2
      np_perProc(i) = np
      np_perProcSum(i+1) = np_perProcSum(i) + np
    end do

    np_perProc(nTasks-1) = np_global - np_perProcSum(nTasks-1)
    np_perProcSum(nTasks) = np_perProcSum(nTasks-1) + np_perProc(nTasks-1)
    np = np_perProc(myrank)

    ! create a provisional node local to global mapping iplg
    ! this will override later with the data from parmetis
    if(allocated(iplg)) deallocate(iplg)
    allocate(iplg(np), stat=stat);
    if(stat/=0) call parallel_abort(' iplg allocate failure')
    do i = 1, np
      iplg(i) = i + np_perProcSum(myrank)
    end do
  end subroutine prePartition


  !-------------------------------------------------------------------------------
  !  Create the connected Nodes array
  !-------------------------------------------------------------------------------

  !> create the connected Nodes array
  !> loop over all elements and their nodes. get then the neighbor nodes
  !> finally calculate the number of sides
  !> alter: maxConnNodes, connNodes_data, ns, ns_global, node%nConnNodes
  subroutine findConnNodes
    use yowNodepool,    only: np, np_global, nodes_global, nodes, maxConnNodes, t_Node, connNodes_data
    use yowElementpool, only: ne_global, INE_global
    use yowSidepool,    only: ns, ns_global
    implicit none

    integer(KIND=JWIM) :: i, j
    type(t_Node), pointer :: node
    character(len=1000) :: errstr

    ! Loop over all nlements
    ! look at their nodes
    ! get the node 1, insert node 2 and 3 into the connected nodes array
    ! do that for node 2 and 3 again

    ! implementation is some different
    ! loop over alle elements to get the # of connected nodes
    ! allocate space
    ! loop a second time to insert the connected nodes

    ! first loop
    do i = 1, ne_global
      do j = 1, 3
        node => nodes_global(INE_global(j,i))
        call node%insertConnNode()
!> \bug why is maxConnNodes twice as long as it should be if we comment in the next line???
        call node%insertConnNode()
      end do
    end do

    maxConnNodes = maxval(nodes_global(:)%nConnNodes)
    nodes_global(:)%nConnNodes = 0

    ! allocate space
    !> \todo we allocate more than we really need
    if(allocated(connNodes_data)) deallocate(connNodes_data)
    allocate(connNodes_data(np_global, maxConnNodes), stat = stat)
    if(stat/=0) then
      write(errstr,*) "connNodes allocation failure", 4*np_global*maxConnNodes*B2MB, "MB"
      ABORT(errstr)
    endif

    ! second loop
    do i = 1, ne_global
  !     do j = 1, 3
  !       node = nodes(elements(i)%node(j))
  !       call node%insertConnNode( elements(i)%node( mod (j+1, 3)+1 ))
  !       call node%insertConnNode( elements(i)%node( mod (j+2, 3)+1 ))
  !     end do

      node => nodes_global(INE_global(1,i))
      call node%insertConnNode(INE_global(2,i))
      call node%insertConnNode(INE_global(3,i))

      node => nodes_global(INE_global(2,i))
      call node%insertConnNode(INE_global(3,i))
      call node%insertConnNode(INE_global(1,i))

      node => nodes_global(INE_global(3,i))
      call node%insertConnNode(INE_global(1,i))
      call node%insertConnNode(INE_global(2,i))

    end do

    ns = 0
    ! calc # sides local
    do i = 1, np
      node => nodes(i)
      ns = ns + node%nConnNodes
    end do

    do i = 1, np_global
      node => nodes_global(i)
      ns_global = ns_global + node%nConnNodes
    end do
  end subroutine


  !-------------------------------------------------------------------------------
  ! Collect all data for parmetis und partition the mesh
  !-------------------------------------------------------------------------------

  !> Collect all data for parmetis und partition the mesh
  !> after that, we knoe for every node the domain ID
  !> alter: t_Node::domainID
  subroutine runParmetis
    use yowDatapool, only: debugParmetis,debugPartition, nTasks, myrank, itype, comm
    use yowNodepool, only: np, np_global, nodes, nodes_global, t_Node, np_perProcSum
    use yowNodepool, only: np_perProc, node2rank
    use yowSidepool, only: ns
    use yowMpiModule
    implicit none

    ! Parmetis
    ! Node neighbor information
    integer(KIND=JWIM) :: wgtflag, numflag, ndims, nparts, edgecut, ncon
    integer(KIND=JWIM), allocatable :: xadj(:), part(:), vwgt(:), adjwgt(:), vtxdist(:), options(:), adjncy(:)
    ! parmetis need single precision
    real(4), allocatable :: xyz(:), tpwgts(:), ubvec(:)

    ! Mics
    integer(KIND=JWIM) :: i, j, ierr
    type(t_Node), pointer :: node, nodeNeighbor
    character(len=1000) :: errstr


    ! Create xadj and adjncy arrays. They holds the nodes neighbors in CSR Format
    ! Here, the adjacency structure of a graph is represented by two arrays,
    ! xadj[n+1] and adjncy[m]; n vertices and m edges. Every edge is listen twice
    allocate(adjncy(ns), stat=stat)
    if(stat/=0) call parallel_abort('adjncy allocation failure')
    allocate(xadj(np+1), stat=stat)
    if(stat/=0) call parallel_abort('xadj allocation failure')

    xadj = 0
    xadj(1) = 1
    adjncy = 0
    do i=1, np
      node => nodes(i)
      xadj(i+1) = xadj(i) + node%nConnNodes
      do j=1, node%nConnNodes
        nodeNeighbor => node%connNodes(j)
        adjncy(j + xadj(i) - 1) = nodeNeighbor%id_global
      end do
    end do

    ! Option for Parmetis
    allocate(options(3))
    options(1)=1   ! 0: default options; 1: user options
    if(debugParmetis) then
      options(2)=15
    else
      options(2)=0  ! Level of information returned: see defs.h in ParMETIS-Lib dir
    endif
    options(3)=15  ! Random number seed

    ! Fortran-style numbering that starts from 1
    numflag = 1
    ! # dimensions of the space in which the graph is embedded
    ndims = 2
    ! # parts the mesh is divide. Usually nTasks
    nparts = nTasks
    ! # weight per node
    ncon = 1


    ! Create weights
    ! keep it simple. ignore weights

    ! wgtflag: 0: none (vwgt and adjwgt are NULL); 1: edges (vwgt is NULL); 2: vertices (adjwgt is NULL); 3: both vertices & edges;
    wgtflag = 0
    allocate(vwgt(np*ncon), stat=stat)
    if(stat/=0) call parallel_abort('vwgt allocation failure')
    !> \todo
    vwgt = 1
    allocate(adjwgt(ns), stat=stat)
    if(stat/=0) call parallel_abort('adjwgt allocation failure')
    !> \todo
    adjwgt = 1

    ! Vertex weight fraction
    allocate(tpwgts(ncon*nTasks),stat=stat)
    if(stat/=0) call parallel_abort('partition: tpwgts allocation failure')
    tpwgts=1.0/real(nTasks)

    ! Imbalance tolerance
    allocate(ubvec(ncon),stat=stat)
    if(stat/=0) call parallel_abort('partition: ubvec allocation failure')
    ubvec=1.01

    ! Partition dual graph
    allocate(xyz(2*np),stat=stat)
    if(stat/=0) call parallel_abort('xyz: ubvec allocation failure')
    do i=1, np
      node => nodes(i)
      xyz(2*(i-1)+1) = REAL(node%x)
      xyz(2*(i-1)+2) = REAL(node%y)
    end do

    ! Partition array returned from ParMeTiS
    allocate(part(np),stat=stat)
    if(stat/=0) call parallel_abort('part: ubvec allocation failure')

    ! ParMeTiS vertex distribution array (starts at 1)
    allocate(vtxdist(nTasks+1),stat=stat)
    if(stat/=0) call parallel_abort('partition: vtxdist allocation failure')

    call mpi_allgather(np_perProcSum(myrank)+1, 1, itype, vtxdist, 1, itype, comm, ierr)
    if(ierr/=MPI_SUCCESS) call parallel_abort('partition: mpi_allgather',ierr)
    vtxdist(nTasks+1)=np_global+1

    ! check vtxdist
    ! myrank stats from 0
    if((vtxdist(myrank+2) - vtxdist(myrank+1)) < 1) then
      write(*,*) "Thread", myrank, "has no nodes"
      write(*,*) "vtxdist", vtxdist
      ABORT("Poor initial vertex distribution detected")
    endif



  !  My notes from manual:
  !  p: # of processors;
  !  n: total # of vertices (local) in graph sense;
  !  m: total # of neighboring vertices ("edges"); double counted between neighboring vertice u and v.
  !  ncon: # of weights for each vertex;
  !  int(in) vtxdist(p+1): Processor j stores vertices vtxdist(j):vtxdist(j+1)-1
  !  int (in) xadj(n+1), adjncy(m):
  !           locally, vertex j's neighboring vertices are adjncy(xadj(j):xadj(j+1)-1). adjncy points to global index;
  !  int(in) vwgt(ncon*n), adjwgt(m): weights at vertices and "edges". Format of adjwgt follows adjncy;
  !  int(in) wgtflag: 0: none (vwgt and adjwgt are NULL);
  !          1: edges (vwgt is NULL); 2: vertices (adjwgt is NULL); 3: both vertices & edges;
  !  int(in) numflag: 0: C-style numbering from 0; 1: FORTRAN style from 1;
  !  int(in) ndims: 2 or 3 (D);
  !  float(in) xyz(ndims*n): coordinate for vertex j is xyz(j*ndims:(j+1)*ndims-1);
  !  int(in)   nparts: # of desired sub-domains (usually nTasks);
  !  float(in) tpwgts(ncon*nparts): =1/nparts if sub-domains are to be of same size for each vertex weight;
  !  float(in) ubvec(ncon): imbalance tolerance for each weight;
  !  int(in)   options: additonal parameters for the routine (see above);
  !  int(out)  edgecut: # of edges that are cut by the partitioning;
  !  int(out)  part(): array size = # of local vertices. It stores indices of local vertices.



    if(debugParmetis .and. myrank == 0) write(*,*) "Run ParMETIS now..."

       write(*,*) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       write(*,*) !!!!!!!!!!!!!!!!  call to ParMETIS_V3_PartGeomKway was disabled !!!!
       write(*,*) !!!!!!!!!!!!!!!!  until work with unstructured grid is continued !!!!
       write(*,*) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!! unwam not available for use !!!!!
!!!    call ParMETIS_V3_PartGeomKway(vtxdist, xadj, adjncy, &
!!!                                  vwgt, & !vwgt - ignore weights
!!!                                  adjwgt, & ! adjwgt - ignore weights
!!!                                  wgtflag, &
!!!                              numflag,ndims,xyz,ncon,nparts,tpwgts,ubvec,options, &
!!!                              edgecut,part, comm)


  ! write(*,*) myrank, "edge cuted", edgecut

    ! Collect the parmetis data from all threads
    ! and create a global node to domain number mapping

    allocate(node2rank(np_global),stat=stat)
    if(stat/=0) then
      write(errstr,*) 'node2rank allocation failure', 4*np_global*B2MB, "MB"
      ABORT(errstr)
    endif
  !
    call mpi_allgatherv(part, np, itype, node2rank, np_perProc, np_perProcSum, itype, comm, ierr)
    if(ierr/=MPI_SUCCESS) call parallel_abort('mpi_allgatherv ',ierr)
  !
    do i = 1, np_global
        node => nodes_global(i)
        ! this node structure is obsolete and will be removed soon
        node%domainID = node2rank(node%id_global)
    end do


    ! write out partition info for katerfempresenter
    if(debugPartition) write(600,*) node2rank

    if(allocated(xadj))        deallocate(xadj)
    if(allocated(adjncy))      deallocate(adjncy)
    if(allocated(part))        deallocate(part)
    if(allocated(vwgt))        deallocate(vwgt)
    if(allocated(adjwgt))      deallocate(adjwgt)
    if(allocated(xyz))         deallocate(xyz)
    if(allocated(tpwgts))      deallocate(tpwgts)
    if(allocated(ubvec))       deallocate(ubvec)
    if(allocated(vtxdist))     deallocate(vtxdist)
  end subroutine runParmetis


  !-------------------------------------------------------------------------------
  ! with the new data from parmetis, recalculate some variables
  !-------------------------------------------------------------------------------

  !> recalculate some variables
  !> parmetis change the number of sides per domain. So recalculate some variables
  !> alter: np, ns, np_perProc, np_perProcSum, iplg, ipgl, nodes_global%id
  !> @note connNodes_data(:) has not changed after the call to parmetis.
  subroutine postPartition
    use yowDatapool, only: myrank, nTasks
    use yowNodepool, only: np_global, np, nodes_global, nodes, np_perProc, np_perProcSum, iplg, ipgl, t_Node
    use yowSidepool, only: ns
    implicit none

    integer(KIND=JWIM) :: i, j
    type(t_Node), pointer :: node

    ! determine how many nodes now belong to which thread
    ! and set the nodes local id
    np_perProc = 0
    do i = 1, np_global
      ! fortran counts from 1. np_perProc from 0
      np_perProc(nodes_global(i)%domainID-1) = np_perProc(nodes_global(i)%domainID-1)+1
      ! set the new local id
      nodes_global(i)%id = np_perProc(nodes_global(i)%domainID-1)
    end do

    np_perProcSum(0) = 0
    do i = 1, nTasks-1
      np_perProcSum(i) = np_perProcSum(i-1) + np_perProc(i-1)
    end do

    np = np_perProc(myrank)

    ! create the new node local to global mapping iplg. This is not the final one.
    ! create a iplg with ghost nodes in findGhostNodes
    if(allocated(iplg)) deallocate(iplg)
    allocate(iplg(np), stat=stat)
    if(stat/=0) call parallel_abort('iplg second allocation failure')
    iplg = 0

    if(allocated(ipgl)) deallocate(ipgl)
    allocate(ipgl(np_global), stat=stat)
    if(stat/=0) call parallel_abort('ipgl allocation failure')
    ipgl = 0

    j = 1
    do i = 1, np_global
      node => nodes_global(i)
      if(node%domainID == myrank+1) then
        iplg(j) = i
        ipgl(i) = j
        j = j + 1
      endif
    end do

    ! calc # sides local again, because the nodes now belongs to another domain
    ns = 0
    do i = 1, np
      node => nodes(i)
      ns = ns + node%nConnNodes
    end do
  end subroutine postPartition


  !-------------------------------------------------------------------------------
  ! find the ghost nodes of the local domain
  !-------------------------------------------------------------------------------

  !> find the ghost nodes of the local domain
  !> alter: ng, ghosts(), ghostlg, ghostgl, npa, iplg
  subroutine findGhostNodes
    use yowDatapool, only: myrank
    use yowNodepool, only: t_Node, np, nodes, ghosts, nodes_global, ng, ghostlg, ghostgl, npa, np_global, iplg
    implicit none
    integer(KIND=JWIM) :: i, j, k
    type(t_Node), pointer :: node, nodeNeighbor, nodeGhost
    !> temporary hold the ghost numbers
    integer(KIND=JWIM), save, allocatable :: ghostTemp(:)

    ! iterate over all local nodes and look at their neighbors
    ! has the neighbor another domain id, than it is a ghost

    ! implementation is some different
    ! loop over all nodes to get the # of ghost nodes
    ! allocate space
    ! loop a second time to insert the ghost nodes
    !> \todo make this faster. dont check all neighbors from all local nodes. mark if an node has already been checked.

    ! first loop. find out how many ghost nodes we have (with double entries)
    ng = 0
    do i = 1, np
      node => nodes(i)
      do j = 1, node%nConnNodes
        nodeNeighbor => node%connNodes(j)
        if(nodeNeighbor%domainID /= node%domainID) then
          ! yes, we found a ghost
          ng = ng + 1
        end if
      end do
    end do

    allocate(ghostTemp(ng), stat=stat)
    if(stat/=0) call parallel_abort('ghostTemp allocation failure')

    allocate(ghostgl(np_global), stat=stat)
    if(stat/=0) call parallel_abort('ghostgl allocation failure')
    ghostgl = 0

    ! second loop. fill ghostlg. ignore double entries
    ng = 0
    ! iterate over all local nodes
    do i = 1, np
      node => nodes(i)

      ! check their neighbors
      secondloop: do j = 1, node%nConnNodes
        nodeNeighbor => node%connNodes(j)

        if(nodeNeighbor%domainID /= node%domainID) then
          ! yes, we found a ghost
          ! check if this ghost is allready in the ghost list
          do k = 1, ng
            nodeGhost => nodes_global(ghostTemp(k))
            if(nodeNeighbor%id_global == nodeGhost%id_global) then
              ! yes, we allready know this ghost.
              ! check the next neighbor
              cycle secondloop
            end if
          end do

          ! no we don't know this ghost. insert it
          ng = ng + 1
          ghostTemp(ng) = nodeNeighbor%id_global
          ghostgl(nodeNeighbor%id_global) = ng
        end if
      end do secondloop
    end do

    ! reallocate the gosttemp array becouse it is longer then the new ng
    if(allocated(ghostlg)) deallocate(ghostlg)
    allocate(ghostlg(ng), stat=stat)
    if(stat/=0) call parallel_abort('ghostlg allocation failure')
    ghostlg = ghostTemp(1:ng)
    deallocate(ghostTemp)

    npa = np + ng

    ! check if ghostlg contains only uniqe values
    do i=1, ng
      do j=i+1, ng
        if(ghostlg(i) == ghostlg(j)) then
          write(*,*) "double global ghost id in ghostlg(i,j)", i, j
          stop "double global ghost id in ghostlg(i,j)"
        endif
      end do
    end do


    ! create the new node local to global mapping iplg with ghost. final one.
    if(allocated(iplg)) deallocate(iplg)
    allocate(iplg(npa), stat=stat)
    if(stat/=0) call parallel_abort('iplg second allocation failure')
    iplg = 0

    j = 1
    do i = 1, np_global
      node => nodes_global(i)
      if(node%domainID == myrank+1) then
        iplg(j) = i
        j = j + 1
      endif
    end do

    iplg(np+1: npa) = ghostlg(1:ng)
  end subroutine findGhostNodes


  !-------------------------------------------------------------------------------
  ! find the number of connected domains and their ghosts
  !-------------------------------------------------------------------------------

  !> find the number of connected domains and their ghosts
  !> 1) Iterate over all ghost nodes and look at their thread id to find # neighbor domains
  !> 2) assign the ghost nodes to their domains
  !> alter: neighborDomains(), nConnDomains
  subroutine findConnDomains
    use yowNodepool,       only: ghosts, ng, t_Node
    use yowDatapool,       only: nTasks
    use yowExchangeModule, only: neighborDomains, initNbrDomains
    implicit none

    integer(KIND=JWIM) :: i, itemp
    type(t_Node), pointer :: ghost

    ! # of ghost per neighbor domain
    integer(KIND=JWIM), allocatable :: numberGhostPerNeighborDomainTemp(:)
    ! look up table. domainID to neighbor (ID)
    integer(KIND=JWIM), allocatable :: domainID2NeighborTemp(:)

    ! Part 1) find # neighbor domains

    ! allocate this array with a fixed size of nTasks. even if we not have
    ! so many neighbor domains, nTasks will never be very large
    allocate(numberGhostPerNeighborDomainTemp(nTasks), stat=stat)
    if(stat/=0) call parallel_abort('numberGhostPerNeighborDomainTemp allocation failure')
    numberGhostPerNeighborDomainTemp = 0

    allocate(domainID2NeighborTemp(nTasks), stat=stat)
    if(stat/=0) call parallel_abort('domainID2NeighborTemp allocation failure')
    domainID2NeighborTemp = 0


    ! iterate over all ghost nodes an get their thread id
    itemp = 0
    do i = 1, ng
      ghost => ghosts(i)

      ! sum how many ghost belongs to the ghost domainID
      numberGhostPerNeighborDomainTemp(ghost%domainID) = numberGhostPerNeighborDomainTemp(ghost%domainID) + 1

      ! check if this ghost domainID is allready in the domains list
      if(domainID2NeighborTemp(ghost%domainID) /= 0) then
        ! yes we have allready insert this domain id
        cycle
      end if

      ! no we dont know this domain id. insert it
      itemp = itemp + 1
      domainID2NeighborTemp(ghost%domainID) = itemp
    end do

    ! Part 2)  assign the ghost nodes to their domains
    call initNbrDomains(itemp)

    do i = 1, nTasks
      if(numberGhostPerNeighborDomainTemp(i) /= 0) then
        neighborDomains(domainID2NeighborTemp(i) )%domainID = i

        allocate(neighborDomains(domainID2NeighborTemp(i))%nodesToReceive(numberGhostPerNeighborDomainTemp(i)), stat=stat)
        if(stat/=0) call parallel_abort('neighborDomains%ghosts allocation failure')
        neighborDomains(domainID2NeighborTemp(i))%nodesToReceive = 0
      end if
    end do

    do i = 1, ng
      ghost => ghosts(i)
      itemp = domainID2NeighborTemp(ghost%domainID)

      neighborDomains(itemp)%numNodesToReceive = neighborDomains(itemp)%numNodesToReceive + 1
      neighborDomains(itemp)%nodesToReceive(neighborDomains(itemp)%numNodesToReceive) = ghost%id_global
    end do

    if(allocated(numberGhostPerNeighborDomainTemp)) deallocate(numberGhostPerNeighborDomainTemp)
    if(allocated(domainID2NeighborTemp)) deallocate(domainID2NeighborTemp)
  end subroutine findConnDomains


  !-------------------------------------------------------------------------------
  ! exchange Ghost Ids so every thread knows which nodes he has to send to the
  ! other parition
  !-------------------------------------------------------------------------------

  !> exchange Ghost Ids so every thread knows which nodes he has to send to the
  !> other parition. every parition has a list of ghost nodes from the other
  !> partition. this data is stored in the neighborDomains variable
  !> \todo make a better  explanation
  !> 1) send to the neighbor domain which ghost nodes we want from him
  !> 2) receive from the neighbor domain which ghost we must send to him
  !> alter: neighborDomains()%{numNodesToSend, nodesToSend},
  subroutine  exchangeGhostIds
    use yowNodepool,       only: np, t_node, nodes
    use yowDatapool,       only: myrank, comm
    use yowExchangeModule, only: neighborDomains, nConnDomains, createMPITypes
    use yowMpiModule
    implicit none

    integer(KIND=JWIM) :: i, j, k
    integer(KIND=JWIM) :: ierr
    ! uniq tag that identify the sender and which information he sends
    integer(KIND=JWIM) :: tag
#ifdef WAM_HAVE_MPI_F08
    ! we use non-blocking send and recv subroutines
    ! store the send status
    type(mpi_request) :: sendRequest(nConnDomains)
    ! store the revc status
    type(mpi_request) :: recvRequest(nConnDomains)
    ! status to verify if one communication fails or not
    type(mpi_status) :: status(nConnDomains);
#else
    ! we use non-blocking send and recv subroutines
    ! store the send status
    integer(KIND=JWIM) :: sendRequest(nConnDomains)
    ! store the revc status
    integer(KIND=JWIM) :: recvRequest(nConnDomains)
    ! status to verify if one communication fails or not
    integer(KIND=JWIM) :: status(MPI_STATUS_SIZE, nConnDomains);
#endif

    type(t_node) , pointer :: node

    ! send to all domain neighbors how many ghosts nodes we want from him and which ones
    do i=1, nConnDomains
        ! create a uniq tag for this domain
        tag = neighborDomains(i)%domainID*10 + 1
        ! send to the neighbor how many ghost nodes we want from him
        call MPI_Isend(neighborDomains(i)%numNodesToReceive, &
                1, &
                MPI_INT, &
                neighborDomains(i)%domainID-1, &
                tag, &
                comm, &
                sendRequest(i), &
                ierr);
        if(ierr/=MPI_SUCCESS) then
          write(*,*) "mpi send failure"
        endif

        tag = neighborDomains(i)%domainID*10 + 2
        ! send to the neighbor which ghost nodes we want from him
        call MPI_Isend(neighborDomains(i)%nodesToReceive, &
                neighborDomains(i)%numNodesToReceive, &
                MPI_INT, &
                neighborDomains(i)%domainID-1, &
                tag, &
                comm, &
!> todo use a second sendRequest array here
                sendRequest(i), &
                ierr);
        if(ierr/=MPI_SUCCESS) then
          write(*,*) "mpi send failure"
        endif

        ! receive from neighbor how many ghost nodes we have to send him
        tag = (myrank+1)*10 + 1
        call MPI_Irecv(neighborDomains(i)%numNodesToSend, &
                1, &
                MPI_INT, &
                neighborDomains(i)%domainID-1, &
                tag, &
                comm, &
                recvRequest(i), &
                ierr)
        if(ierr/=MPI_SUCCESS) then
          write(*,*) "mpi recv failure"
        endif
    end do

    ! wait for communication end
    call MPI_Waitall(nConnDomains, recvRequest, status, ierr)

    ! test for all neighbor domains
    do i=1, nConnDomains
        ! test if the neighbor wants more ghost nodes than we have
        if(neighborDomains(i)%numNodesToSend > np) then
          write(*,'(i5, a, i5, a, i5,a, i5, a, /)', advance='no') myrank, " ERROR neighbordomain ", neighborDomains(i)%domainID, &
                    " wants ", neighborDomains(i)%numNodesToSend, &
                    " nodes, but we have only ", np, " nodes"
            ABORT("")
        end if
    end do

    ! receive from all neighbor domains which nodes we must send him
    do i=1, nConnDomains
        allocate(neighborDomains(i)%nodesToSend(neighborDomains(i)%numNodesToSend))
        neighborDomains(i)%nodesToSend = 0

        ! receive from neighbor which nodes we must send
        tag = (myrank+1)*10 + 2
        call MPI_Irecv(neighborDomains(i)%nodesToSend, &
                neighborDomains(i)%numNodesToSend, &
                MPI_INT, &
                neighborDomains(i)%domainID-1, &
                tag, &
                comm, &
                recvRequest(i), &
                ierr)
        if(ierr/=MPI_SUCCESS) then
          PARALLEL_ABORT("mpi recv failure", ierr)
        endif
    end do

    ! wait for communication end
    call MPI_Waitall(nConnDomains, recvRequest, status, ierr)

  ! test for all neighbor domains
    do i=1, nConnDomains
      ! test if the neighbor wants nodes that we don't own
      outerloop: do j=1, neighborDomains(i)%numNodesToSend

        ! compare with all local nodes
        do k=1, np
          node => nodes(k)
          if(node%id_global == neighborDomains(i)%nodesToSend(j)) then
            cycle outerloop
          end if
        end do

        write(*,*) myrank, "Neighbordomain", neighborDomains(i)%domainID, &
                  " want Node", neighborDomains(i)%nodesToSend(j), &
                  " but we don't own this node"
        stop
      end do outerloop
    end do

    call createMPITypes()
  !   do i=1, nConnDomains
  !     write(*,*)  myrank, "send to ", neighborDomains(i)%domainID-1,"this ", neighborDomains(i)%nodesToSend
  !   end do
  end subroutine exchangeGhostIds

  !> this collects all data which depends on ghost information
  !> alter: ne, INE, x, y, z, ielg, iegl
  subroutine postPartition2()
    use yowElementpool, only: ne, ne_global, INE, ielg, iegl, INE_global, eleBelongTo
    use yowDatapool,    only: myrank
    use yowNodepool,    only: np, nodes_global, t_Node, ghostlg, ng, npa, nodes, ghosts
    use yowNodepool,    only: x, y, z
    implicit none

    integer(KIND=JWIM) :: i, j, k
    type(t_Node), pointer :: node
    logical :: assigned

    ! to find the local elements, iterate over all global elements and check their node domain IDs

    ! step 1: calc the number of local elements
    ne = 0
    do i=1, ne_global
      if(eleBelongTo(INE_global(:,i)) .eqv. .true.) then
        ne = ne +1
!         write(*,*) myrank, "ele global", i-1
      endif
    end do

    ! step 2: fill the local element index array
    if(allocated(INE)) deallocate(INE)
    allocate(INE(3, ne), stat=stat)
    if(stat/=0) call parallel_abort('INE allocation failure')
    INE = 0

    ne = 0
    do i=1, ne_global
      ! yes, this element belongs to this domain
      if(eleBelongTo(INE_global(:,i)) .eqv. .true.) then
        ne = ne + 1
        do j=1, 3
          assigned = .false.
          node => nodes_global(INE_global(j,i))
          
          if(node%domainID == myrank+1) then
            INE(j, ne) = node%id
            assigned = .true.
          else
            ! the element have a ghost node
            !> \todo create some localnode to localghost mapping
            !> What number is this ghost
            do k=1, ng
              if(node%id_global == ghostlg(k)) then
                ! conversion: the ghost nodes are stored behind the local nodes.
                if(INE(j,ne) /= 0) then
                  write(*,*) "will write to INE(j, ne) but there is allready a value", j, ne, INE(j, ne) 
                endif
                INE(j, ne) = np + k
                node%id = np+k
                assigned = .true.
!                 write(*,*) myrank, "node to ele", node%id_global-1, i-1, np+k
                exit
              endif
            end do
          endif
          if(assigned .eqv. .false.) then
            write(*,*) "Can't assign global node to INE", node%id_global
          endif
        end do
      endif
    end do

    ! check if INE contains 0. 
    do i=1, ne
      if(MINVAL(ABS(INE(:,i))) == 0) then
        write(*,*) "0 in INE ne=", ne
        stop "0 in INE"
      endif
    end do

    if(MAXVAL(INE) /= npa) then
      write(*,*) "MAXVAL(INE) /= npa ERROR?"
    endif

    ! create element local to global mapping ielg
    ! and element global to local mapping iegl
    if(allocated(ielg)) deallocate(ielg)
    allocate(ielg(ne), stat=stat)
    if(stat/=0) call parallel_abort('ielg allocation failure')
    ielg = 0

    if(allocated(iegl)) deallocate(iegl)
    allocate(iegl(ne_global), stat=stat)
    if(stat/=0) call parallel_abort('iegl allocation failure')
    iegl = 0

    j = 0
    do i=1, ne_global
      if(eleBelongTo(INE_global(:,i)) .eqv. .true.) then
        j = j +1
        ielg(j) = i
        iegl(i) = j
      end if
    end do



    ! fill the local x,y arrays
    if(allocated(x)) deallocate(x)
    allocate(x(npa), stat=stat)
    if(stat/=0) call parallel_abort('x allocation failure')

    if(allocated(y)) deallocate(y)
    allocate(y(npa), stat=stat)
    if(stat/=0) call parallel_abort('y allocation failure')

    if(allocated(z)) deallocate(z)
    allocate(z(npa), stat=stat)
    if(stat/=0) call parallel_abort('z allocation failure')

    do i=1, np
      node => nodes(i)
      x(i) = node%x
      y(i) = node%y
      z(i) = node%z
    end do

    do i=1, ng
      node => ghosts(i)
      x(np+i) = node%x
      y(np+i) = node%y
      z(np+i) = node%z
    end do

    call createConnNodesLocal()
  end subroutine


  !> create connected nodes array for local nodes
  !> alter CONN, NCONN
  subroutine createConnNodesLocal()
    use yowDatapool, only: myrank
    use yowNodepool, only: np, ng, npa, np_global, ipgl, iplg, ghostlg, nodes, maxConnNodes, CONN, NCONN, t_Node
    use yowElementpool, only: INE_global, ne_global
    implicit none
    integer(KIND=JWIM) :: IP_global, IE_global, IP, i, iteration

    ! to create a sorted patch, we need the connected elements per node
    integer(KIND=JWIM) :: maxNCONE
    integer(KIND=JWIM), allocatable :: NCONE(:), CONE(:,:)
    
    ! the three node IDs of an element
    integer(KIND=JWIM) :: eleNodes(3), lut(npa)

    logical :: skipElement(maxConnNodes)
    integer(KIND=JWIM) :: curNode, nextNode
  
    ! remove this later ;)
    integer(KIND=JWIM), allocatable :: myipgl(:)
    type(t_node), pointer :: node

    ! create connected elements array

    ! We need the connected element information. But there are no ghost elements in pdlib.
    ! Iterate over the global element array and store the global element number for each local node

    allocate(NCONE(np))
    NCONE(:) = 0

    ! first detect the maxConnNodes
    do IE_global=1, ne_global
      do i=1, 3
        IP = ipgl(INE_global(i, IE_global))
        if(IP /= 0) then
          NCONE(IP) = NCONE(IP) + 1
        endif
      end do
    end do

    maxNCONE = maxval(NCONE)
    NCONE(:) = 0 ! reset this,because we use this variable as a counter how many elements are currently connected
    allocate(CONE(maxNCONE, np))
    CONE(:,:) = 0

    ! second, store for every local node the connected global element numbers
    do IE_global=1, ne_global
      do i=1, 3
        IP = ipgl(INE_global(i, IE_global))
        if(IP /= 0) then
          NCONE(IP) = NCONE(IP) + 1
          CONE(NCONE(IP), IP) = IE_global 
        endif
      end do
    end do

    do ip=1, np
      if(NCONE(ip) == 0) then
        write(*,*) myrank, "Global node is not connected to any element", iplg(ip)
        ABORT("NCONE ERROR")
     endif
   end do
    
    do ip=1, np
      if(any(CONE(1:NCONE(ip), ip) .eq. 0)) then
        write(*,*) myrank, "Global node is not connected to any element", iplg(ip)
        ABORT("CONE ERROR")
      endif
    end do

    ! myipgl is the same as ipgl, expect it maps the global->ghost nodes too
    allocate(myipgl(np_global))
    myipgl(:) = 0
    myipgl(:) = ipgl(:)
    do i=1, ng
      myipgl(ghostlg(i)) = np+i
    end do


    ! now build the connected nodes arrays

    if(allocated(NCONN)) deallocate(NCONN)
    allocate(NCONN(np))
    NCONN(:) = 0

    if(allocated(CONN)) deallocate(CONN)
    allocate(CONN(maxConnNodes, np))
    CONN(:,:) = 0

    ! stuff to create a sorted patch
    ! first we must distinguish between boundary nodes and domain nodes.
    ! For domain nodes, the patch is simply all connected nodes.
    ! For boundary nodes, the patch is not a closed cycle. So we have a begin and end node.

   
    do IP = 1, np
      IP_global = iplg(IP)
      skipElement(:) = .false.   

      curNode = 0
     
      ! attention: do not use NCONN here because we create this array first
      node => nodes(IP)

      ! IP is a domain node
      if(NCONE(IP) == node%nConnNodes) then
        IE_global = CONE(1, IP)
        eleNodes(:) = INE_global(:, IE_global)

        if(eleNodes(1) == IP_global) then
          curNode = eleNodes(3)
        else if(eleNodes(2) == IP_global) then
          curNode = eleNodes(3)
        else if(eleNodes(3) == IP_global) then
          curNode = eleNodes(2)
        endif

      ! IP is a boundary node
      else
        ! the patch is not closed if IP is a boundary node. so search for an end of the patch line
        lut(:) = 0
        do i=1, NCONE(IP)
          IE_global = CONE(i, IP)
          eleNodes(:) = INE_global(:, IE_global)
          lut(myipgl(eleNodes(:))) = lut(myipgl(eleNodes(:))) + 1                   
        end do

        do i=1, npa
          if(lut(i) == 1 .and. i /= IP) then
            NCONN(IP) = NCONN(IP) + 1
            CONN(NCONN(IP), IP) = i
              curNode = iplg(i)
            exit
          endif
        end do
  
      endif

      do iteration=1, NCONE(IP)
        do i=1, NCONE(IP)
          if(skipElement(i) .eqv. .true.) cycle

          IE_global = CONE(i, IP)
          eleNodes(:) = INE_global(:, IE_global)

          nextNode = 0
          if(eleNodes(1) == curNode) then
            if(eleNodes(2) == IP_global) then
              nextNode = eleNodes(3)
            else
              nextNode = eleNodes(2)
            endif
          else if(eleNodes(2) == curNode) then
            if(eleNodes(1) == IP_global) then
              nextNode = eleNodes(3)
            else
              nextNode = eleNodes(1)
            endif
          else if(eleNodes(3) == curNode) then
            if(eleNodes(1) == IP_global) then
              nextNode = eleNodes(2)
            else
              nextNode = eleNodes(1)
            endif
          endif   

          if(nextNode /= 0) then
            curNode = nextNode;
            skipElement(i) = .true.

            NCONN(IP) = NCONN(IP) + 1
            CONN(NCONN(IP), IP) = myipgl(curNode)
          endif
        end do
      end do
    end do

   
!     ! create an unsorted connected nodes array
!     do ip=1, np
!       node => nodes(ip)
!       
!       NCONN(IP) = node%nConnNodes
!       do i=1, NCONN(IP)
!         nodeNeighbor => node%connNodes(i)
!         CONN(i, IP) =  myipgl(nodeNeighbor%id_global)
!       end do
!     end do

    deallocate(myipgl)

  end subroutine

  !> Deallocates all global lenght mesh arrays to save memory: nodes_global, INE_global, connNodes_data
  !> All mappings and local lenght arrays are not touched
  subroutine clearGlobalMeshVars()
    use yowNodepool, only: nodes_global, connNodes_data
    use yowElementpool, only: INE_global
    implicit none

    if(allocated(nodes_global))   deallocate(nodes_global)
    if(allocated(connNodes_data)) deallocate(connNodes_data)
    if(allocated(INE_global))     deallocate(INE_global)

    ! don't clear node2rank. we need this for the nodeexchange call
  end subroutine


  subroutine finalizePD()
    use yowExchangeModule, only: finalizeExchangeModule
    use yowNodepool,    only: finalizeNodepool
    use yowElementpool, only: finalizeElementpool
    use yowRankModule,  only: finalizeRankModule
    implicit none

    call finalizeRankModule()
    call finalizeExchangeModule()   
    call finalizeElementpool()
    call finalizeNodepool()
  end subroutine

  !> This is a simple function to search for an available unit.
  !> @param[out] unit optional the new free unit number
  !> @return next free unit number. If no units are available, -1 is returned.
  integer function newunit(unit)
    implicit none
    integer(KIND=JWIM), intent(out), optional :: unit

  ! LUN_MIN and LUN_MAX define the range of possible LUNs to check.
    integer(KIND=JWIM), parameter :: LUN_MIN=10, LUN_MAX=1000
    logical :: opened
    integer(KIND=JWIM) :: lun

    newunit=-1
    do lun=LUN_MIN, LUN_MAX
      inquire(unit=lun,opened=opened)
      if (.not. opened) then
        newunit=lun
        exit
      end if
    end do
    if (present(unit)) unit=newunit
  end function newunit
end module yowpdlibMain
