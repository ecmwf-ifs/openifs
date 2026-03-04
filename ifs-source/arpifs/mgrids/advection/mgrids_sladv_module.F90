! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction
! 
! (C) Copyright 1989- Meteo-France.
! 

#define __FILENAME__ "mgrids_sladv_module.F90"
#define ABORT( what ) fckit_exception%abort(what,__FILENAME__,__LINE__)

module mgrids_sladv_module
#ifdef WITH_MGRIDS

use atlas_module
use dwarf_module
use fckit_module              , only : fckit_exception
use dwarf_sladv_module        , only : dwarf_sladv
use mgrids_advection_module   , only : mgrids_advection_args, mgrids_advection_config
use yomrip                    , only : trip
use yomdyn                    , only : tdyn
use geometry_mod              , only : geometry
use parkind1                  , only : jprb, jprd, jpim

implicit none
private

public ::  create_mgrids_sladv
!------------------------------------------------------------------------------

type, private :: mgrids_interpolation
  type(atlas_Interpolation)                   :: interpolation_to_IFS
  type(atlas_Interpolation)                   :: interpolation_from_IFS
  type(atlas_functionspace_StructuredColumns) :: fs_from_IFS
  type(atlas_functionspace_StructuredColumns) :: fs_from_adv
  type(atlas_functionspace_StructuredColumns) :: fs_to_IFS
  type(atlas_functionspace_StructuredColumns) :: fs_to_adv
  integer(jpim)                               :: ngp_IFS
  integer(jpim)                               :: ngp_adv
contains
  procedure :: to_IFS
  procedure :: from_IFS
end type

!------------------------------------------------------------------------------

type, public, extends( dwarf_sladv ) :: mgrids_sladv

  type(trip)     , pointer :: yrrip       => null()
  type(tdyn)     , pointer :: yrdyn       => null()
  integer                  :: nflevg
  integer                  :: ngptot
  integer                  :: nproma
  integer                  :: ngpblks

  real(kind=jprb), pointer :: ztracer(:,:,:) => null()
    ! shape=( ygfl%nfmg, ydgeometry%yrdimv%nflevg, ydgeometry%yrgem%ngptot )
  real(kind=jprb), pointer :: zuvw0(:,:,:)   => null()
    ! shape=( 3, ydgeometry%yrdimv%nflevg, ydgeometry%yrgem%ngptot )
  real(kind=jprb), pointer :: zuvw9(:,:,:)   => null()
    ! shape=( 3, ydgeometry%yrdimv%nflevg, ydgeometry%yrgem%ngptot )

  type(mgrids_interpolation) :: interpolation

  logical       :: lremesh
  type(atlas_Field) :: field_IFS_uvw0
  type(atlas_Field) :: field_IFS_uvw9
  type(atlas_Field) :: field_IFS_tracer
  type(atlas_Field) :: field_adv_tracer

contains
  procedure, public :: setup
  procedure, public :: preprocess
  procedure, public :: postprocess
end type

!------------------------------------------------------------------------------

interface mgrids_sladv
  module procedure constructor
end interface

!------------------------------------------------------------------------------

contains

!------------------------------------------------------------------------------

function constructor(  &
  config,              &
  ydrip,               &
  yddyn,               &
  ydgeometry ) result(this)
  !-----------------------------------------------------------
  class(mgrids_sladv), pointer :: this
  type(mgrids_advection_config)      :: config
  type(trip),     intent(in), target :: ydrip
  type(tdyn),     intent(in), target :: yddyn
  type(geometry), intent(in), target :: ydgeometry
  !-----------------------------------------------------------
  type(atlas_Trace) :: trace
  trace = atlas_Trace(__FILENAME__,__LINE__,"mgrids_sladv::constructor")

  allocate( mgrids_sladv :: this )
  call this%setup(                                           &
              config,                                        &
              ydrip,                                         &
              yddyn,                                         &
              ydgeometry )
  call trace%final()
end function

!------------------------------------------------------------------------------

subroutine create_mgrids_sladv(  &
  this,                          &
  config,                        &
  ydrip,                         &
  yddyn,                         &
  ydgeometry )
  !-----------------------------------------------------------
  class(dwarf), pointer, intent(out)        :: this
  type(mgrids_advection_config), intent(in) :: config
  type(trip),     intent(in), target :: ydrip
  type(tdyn),     intent(in), target :: yddyn
  type(geometry), intent(in), target :: ydgeometry
  !-----------------------------------------------------------
  this => mgrids_sladv( config,ydrip,yddyn,ydgeometry )
end subroutine

!------------------------------------------------------------------------------

subroutine setup_sladv_grid_and_distribution(this,grid, distribution, gridname, ydgeometry)
  use yomlun,        only : nulnam
  use fckit_module,  only : fckit_log
  !-----------------------------------------------------------
  class(mgrids_sladv),          intent(inout) :: this
  type(atlas_Grid),             intent(inout) :: grid
  type(atlas_GridDistribution), intent(inout) :: distribution
  character(len=*),             intent(in)    :: gridname
  type(geometry),               intent(in)    :: ydgeometry
  !-----------------------------------------------------------
  type(atlas_MeshGenerator) :: meshgenerator
  type(atlas_Partitioner)   :: partitioner
  type(atlas_Config)        :: config
  character(len=1024)       :: string
  !-----------------------------------------------------------

  if( trim(gridname) == "same" ) then

    call fckit_log%info( "    Setup mgrids SLADV geometry using same grid as IFS" )

    grid         = ydgeometry%yratlas%grid
    distribution = ydgeometry%yratlas%griddistribution

  else

    write(string,'(A,A,A)') "    Setup mgrids SLADV geometry using grid ", trim(gridname)
    call fckit_log%info( string )

    grid         = atlas_Grid(trim(gridname))
    partitioner  = atlas_MatchingMeshPartitioner(ydgeometry%yratlas%mesh)
    distribution = partitioner%partition(grid)

    call partitioner%final()
  endif
end subroutine

!------------------------------------------------------------------------------

subroutine setup(                         &
    this,                                 &
    mgrids_config,                        &
    ydrip,                                &
    yddyn,                                &
    ydgeometry                            &
  )
  !-----------------------------------------------------------
  use fckit_module,      only : fckit_log
  use yom_ygfl,          only : ygfl
  use yomcst,            only : ra
  !-----------------------------------------------------------
  class(mgrids_sladv), intent(inout) :: this
  type(mgrids_advection_config)      :: mgrids_config
  type(trip),     intent(in), target :: ydrip
  type(tdyn),     intent(in), target :: yddyn
  type(geometry), intent(in), target :: ydgeometry
  !-----------------------------------------------------------
  integer(jpim) :: ndp_iter
  integer(jpim) :: dp_extrap
  type(atlas_config)           :: config
  type(atlas_Grid)             :: grid
  type(atlas_griddistribution) :: distribution
  type(atlas_Nabla)            :: nabla
  type(atlas_Vertical)         :: vertical
  integer :: jlev, nflevg
  character(len=:), allocatable :: gridID
  character(len=2048) :: msg
  real(kind=jprd), allocatable :: vetaf(:)
  !-----------------------------------------------------------
  call fckit_log%info( "Setup mgrids advection using semi-lagrangian scheme" )

  this%yrrip   => ydrip
  this%yrdyn   => yddyn
  this%nflevg  = ydgeometry%yrdimv%nflevg
  this%ngptot  = ydgeometry%yrgem%ngptot
  this%nproma  = ydgeometry%yrdim%nproma
  this%ngpblks = ydgeometry%yrdim%ngpblks

  nflevg = this%nflevg

  this%lremesh = mgrids_config%lremesh

  config = atlas_config()

  call config%set("radius",       ra )

  call config%set("halo",         mgrids_config%sladv_nhalo)

  call config%set("departurepoint_method", mgrids_config%sladv_departurepoint_method )
  call config%set("departurepoint_iterations", this%yrdyn%nitmp )
  if ( mgrids_config%sladv_departurepoint_method == "SETTLSVF" ) then
    call config%set("SETTLSVF_levels", this%yrdyn%nflevsf )
    call config%set("SETTLSVF_scale", this%yrdyn%rscale )
    call config%set("SETTLSVF_scale_off", this%yrdyn%rscaleoff )
  endif

  call config%set("interpolation_method",  mgrids_config%sladv_interpolation_method )
  call config%set("interpolation_limiter", mgrids_config%sladv_interpolation_limiter )

  call config%set("grid",         mgrids_config%cgrid)
  call config%set("dt", this%yrrip%tdt )
  call config%set("tracers", ygfl%nfmg )

  call fckit_log%info( "    sladv configuration = "//config%json() )

  if( .not. config%get("grid",gridID) )       &
    & call ABORT("grid parameter missing")
  if( .not. config%get("dt",this%dt) )         &
    & call ABORT("dt parameter missing")
  if( .not. config%get("tracers",this%ntrac) )   &
    & call ABORT("tracers parameter missing")
  if( .not. config%get("interpolation_method",this%interpolation_method) ) &
    & call ABORT("interpolation_method parameter missing")
  if( .not. config%get("interpolation_limiter",this%interpolation_limiter) ) &
    & call ABORT("interpolation_limiter parameter missing")

  call setup_sladv_grid_and_distribution( this, grid, distribution, gridID, ydgeometry )

  if( this%needs_nabla(config) ) then !!! experimental !!!
     nabla = create_nabla()
  endif
  ! vetaf contains full levels but also extra 0. at vetaf(0) and 1. at vetaf(nflevg+1)
  ! The valid range is vetaf(1:nflevg)
  ! We need to strip the 0. and 1. and provide these limits as interval instead
  allocate(vetaf(nflevg))
  do jlev=1,nflevg
    vetaf(jlev) = ydgeometry%yrveta%vetaf(jlev)
  enddo
  vertical = atlas_Vertical( levels=vetaf(1:nflevg), interval=[0._jprd,1._jprd] )
  call this%dwarf_sladv_setup(             &
    config,                                &
    atlas_functionspace_StructuredColumns( &
      grid=grid,                           &
      distribution=distribution,           &
      halo=mgrids_config%sladv_nhalo,      &
      vertical=vertical                    &
    ),                                     &
    nabla                                  &
  )
  call vertical%final()
  call distribution%final()
  call grid%final()
  call config%final()
  deallocate(vetaf)

  if( this%lremesh ) then
    call fckit_log%info( "    Setup mgrids interpolators..." )
    this%interpolation%fs_from_IFS    = ydgeometry%yratlas%fs_structuredcolumns
    this%interpolation%fs_from_adv    = this%geometry%fs_structuredcolumns ! setup during "dwarf_sladv_setup()"
    this%interpolation%fs_to_IFS      = this%interpolation%fs_from_IFS
    this%interpolation%fs_to_adv      = this%interpolation%fs_from_adv
    this%field_adv_tracer = this%geometry%fs_structuredcolumns%create_field(name="tracer",&
      levels=nflevg,variables=this%ntrac,kind=atlas_real(jprb))
    config = atlas_Config()
    call config%set("type","quasicubic2D")
    this%interpolation%interpolation_from_IFS = atlas_Interpolation(config, &
      & this%interpolation%fs_from_IFS, this%interpolation%fs_to_adv)
    this%interpolation%interpolation_to_IFS = atlas_Interpolation(config, &
        & this%interpolation%fs_from_adv, this%interpolation%fs_to_IFS)
    call config%final()
    this%interpolation%ngp_IFS = ydgeometry%yrgem%ngptot
    this%interpolation%ngp_adv = this%geometry%ngp

    call fckit_log%info( "    Setup mgrids interpolators... done" )
  else
    call fckit_log%info( "    No mgrids interpolators required ( lremesh = false )" )
  endif

  if( this%lremesh ) then
    this%field_IFS_uvw0 = this%interpolation%fs_from_IFS%create_field( kind=atlas_real(jprb), levels=nflevg, &
      & name="uvw0", variables=3, type="vector" )

    this%field_IFS_uvw9 = this%interpolation%fs_from_IFS%create_field( kind=atlas_real(jprb), levels=nflevg, &
      & name="uvw9", variables=3, type="vector" )

    this%field_IFS_tracer = this%interpolation%fs_from_IFS%create_field( kind=atlas_real(jprb), levels=nflevg, &
      & name="tracer", variables=this%ntrac, type="scalar" )

  else
    this%field_IFS_uvw0   = this%wind_t0
    this%field_IFS_uvw9   = this%wind_t9
    this%field_IFS_tracer = this%tracer_field
  endif
  call this%field_IFS_uvw0%data( this%zuvw0 )
  call this%field_IFS_uvw9%data( this%zuvw9 )
  call this%field_IFS_tracer%data( this%ztracer )

contains

  function create_nabla() result(nabla)
    !! nabla is required to support the non-iterative departure-point calculations (experimental)
    type(atlas_Nabla) :: nabla
    type(atlas_MeshGenerator) :: meshgenerator
    type(atlas_Mesh) :: mesh
    type(atlas_Config) :: fvm_config
    type(atlas_fvm_method) :: fvm
    real(jprb) :: radius
    logical :: present
    call fckit_log%info("Create default nabla")

    radius = ra
    present = config%get("radius",radius)
    
    meshgenerator = atlas_MeshGenerator()
    mesh = meshgenerator%generate(grid)
    fvm_config = atlas_Config()
    call fvm_config%set('halo',1)
    call fvm_config%set("radius", radius )
    fvm = atlas_fvm_method(mesh,fvm_config)
    nabla = atlas_nabla(fvm)

    call meshgenerator%final()
    call mesh%final()
    call fvm%final()
    call fvm_config%final()
  end function

end subroutine

!------------------------------------------------------------------------------

subroutine preprocess(this,args)
  use yom_ygfl,      only : ygfl
  use yomct3,        only : nstep
  !-----------------------------------------------------------
  class(mgrids_sladv), intent(inout) :: this
  class(dwarf_args),    intent(inout) :: args
  !-----------------------------------------------------------
  integer(jpim) :: itrac, jgfl, jkglo, icend, ibl, jlev, jrof
  type(mgrids_advection_args) :: arguments
  type(atlas_Trace) :: trace, trace_transposition, trace_interpolation
  !-----------------------------------------------------------
  trace = atlas_Trace(__FILENAME__,__LINE__,"mgrids_sladv::preprocess")
  arguments = mgrids_advection_args(args)
  this%istep = nstep

  associate(                                   &
    nfmg    => ygfl%nfmg,                      &
    nmgflds => ygfl%nmgflds,                   &
    ycomp   => ygfl%ycomp,                     &
    nflevg  => this%nflevg,                    &
    ngptot  => this%ngptot,                    &
    nproma  => this%nproma,                    &
    ngpblks => this%ngpblks,                   &
    yt0     => arguments%yrgmv%yt0,            &
    yt9     => arguments%yrgmv%yt9,            &
    pwrl9   => arguments%zwrl9,                &
    pgmv    => arguments%zgmv                  &
  )

  trace_transposition = atlas_Trace(__FILENAME__,__LINE__,"transposition")

  do itrac=1,nfmg
    jgfl=nmgflds(itrac)
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jkglo,icend,ibl,jrof,jlev)
    do jkglo=1,ngptot,nproma
      icend=min(nproma,ngptot-jkglo+1)
      ibl=(jkglo-1)/nproma+1
      do jlev=1,nflevg
        do jrof=1,icend
          this%ztracer(itrac,jlev,jkglo+jrof-1) = arguments%zgfl(jrof,jlev,ycomp(jgfl)%mp,ibl)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
  enddo

  !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jkglo,icend,ibl,jrof,jlev)
  do jkglo=1,ngptot,nproma
    icend=min(nproma,ngptot-jkglo+1)
    ibl=(jkglo-1)/nproma+1
    do jlev=1,nflevg
      do jrof=1,icend
        this%zuvw0(1,jlev,jkglo+jrof-1) = pgmv(jrof,jlev,yt0%mu,ibl)
        this%zuvw0(2,jlev,jkglo+jrof-1) = pgmv(jrof,jlev,yt0%mv,ibl)
        this%zuvw0(3,jlev,jkglo+jrof-1) = pgmv(jrof,jlev,yt9%medot,ibl) ! note this is actually holding yt0 value!
        this%zuvw9(1,jlev,jkglo+jrof-1) = pgmv(jrof,jlev,yt9%mu,ibl)
        this%zuvw9(2,jlev,jkglo+jrof-1) = pgmv(jrof,jlev,yt9%mv,ibl)
        this%zuvw9(3,jlev,jkglo+jrof-1) = pwrl9(jrof,jlev,ibl)
      enddo
    enddo
  enddo
  !$OMP END PARALLEL DO

  call this%field_IFS_tracer%set_dirty()
  call this%field_IFS_uvw0%set_dirty()
  call this%field_IFS_uvw9%set_dirty()


  call trace_transposition%final()

  if (this%lremesh) then
    trace_interpolation = atlas_Trace(__FILENAME__,__LINE__,"interpolation")
!!    if( this%istep==0 ) then
    !-----------------------------------------------------------------------
    ! In a real app ptracer must be always interp as it is not only 
    ! advected but modified by other parts of ifs (phys, chemistry)
    ! and the changes must be put in trac
    !-----------------------------------------------------------------------
     call this%interpolation%from_IFS( this%field_IFS_tracer, this%tracer_field )
!!    endif
    call this%interpolation%from_IFS( this%field_IFS_uvw0, this%wind_t0 )
    call this%interpolation%from_IFS( this%field_IFS_uvw9, this%wind_t9 )
    call trace_interpolation%final()
  endif

  end associate

  call trace%final()

end subroutine

!------------------------------------------------------------------------------

subroutine postprocess(this,args)
  use yom_ygfl, only : ygfl
  !-----------------------------------------------------------
  class(mgrids_sladv), intent(inout) :: this
  class(dwarf_args),    intent(inout) :: args
  !-----------------------------------------------------------
  integer(jpim) :: itrac, jgfl, jkglo, icend, ibl, jlev, jrof, jnode
  real(jprb), pointer :: tracer(:,:,:), adv_tracer(:,:,:)
  type(mgrids_advection_args) :: arguments
  type(atlas_Trace) :: trace, trace_transposition, trace_interpolation
  !-----------------------------------------------------------
  trace = atlas_Trace(__FILENAME__,__LINE__,"mgrids_sladv::postrocess")
  arguments = mgrids_advection_args(args)

  associate(                                   &
    nfmg    => ygfl%nfmg,                      &
    nmgflds => ygfl%nmgflds,                   &
    ycomp   => ygfl%ycomp,                     &
    nflevg  => this%nflevg,                    &
    ngptot  => this%ngptot,                    &
    nproma  => this%nproma,                    &
    pgflt1  => arguments%zgflt1                &
  )

  if (this%lremesh) then
    trace_interpolation = atlas_Trace(__FILENAME__,__LINE__,"interpolation")
    call this%tracer_field%data( tracer )
    call this%field_adv_tracer%data( adv_tracer )
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode)
    do jnode=1,this%interpolation%ngp_adv
      adv_tracer(:,:,jnode) = tracer(:,:,jnode)
    enddo
    !$OMP END PARALLEL DO 
    call this%field_adv_tracer%set_dirty()
    call this%interpolation%to_IFS( this%field_adv_tracer, this%field_IFS_tracer )
    call trace_interpolation%final()
  endif

  trace_transposition = atlas_Trace(__FILENAME__,__LINE__,"transposition")
  do itrac=1,nfmg
    jgfl=nmgflds(itrac)
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jkglo,icend,ibl,jrof,jlev)
    do jkglo=1,ngptot,nproma
      icend=min(nproma,ngptot-jkglo+1)
      ibl=(jkglo-1)/nproma+1
      do jlev=1,nflevg
        do jrof=1,icend
          pgflt1(jrof,jlev,ycomp(jgfl)%mp1,ibl)=this%ztracer(itrac,jlev,jkglo+jrof-1)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
  enddo
  call trace_transposition%final()
  end associate
  call trace%final()
end subroutine

!------------------------------------------------------------------------------

subroutine from_IFS(this,field_IFS,field_adv)
  class(mgrids_interpolation), intent(inout) :: this
  type(atlas_Field), intent(in)    :: field_IFS
  type(atlas_Field), intent(inout) :: field_adv
  call field_IFS%halo_exchange()
  call this%interpolation_from_IFS%execute( field_IFS, field_adv )
end subroutine

!------------------------------------------------------------------------------

subroutine to_IFS(this,field_adv,field_IFS)
  class(mgrids_interpolation), intent(inout) :: this
  type(atlas_Field), intent(in)    :: field_adv
  type(atlas_Field), intent(inout) :: field_IFS
  call field_adv%halo_exchange()
  call this%interpolation_to_IFS%execute( field_adv, field_IFS )
end subroutine

!------------------------------------------------------------------------------

#endif
end module
