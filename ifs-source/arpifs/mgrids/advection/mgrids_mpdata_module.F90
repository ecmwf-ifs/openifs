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

#define __FILENAME__ "mgrids_mpdata_module.F90"

!------------------------------------------------------------------------------------------------------
! Call MPDATA advection scheme on new Atlas mesh.
! NOTE: This is a preliminary unoptimised version is working on height levels in the vertical
!       and restcited by the high vertical courant number of IFS.
!       
! Authors: M.Diamantakis & W. Deconinck June 2018
!
!------------------------------------------------------------------------------------------------------
module mgrids_mpdata_module
#ifdef WITH_MGRIDS

use atlas_module
use dwarf_module
use dwarf_mpdata_module,               only : dwarf_mpdata
use mgrids_advection_module,           only : mgrids_advection_args, mgrids_advection_config
use yomrip,                            only : trip
use yomdyn,                            only : tdyn
use geometry_mod,                      only : geometry
use parkind1,                          only : jprb, jpim, jpib
use advection_MPDATA_MPDATA_module,    only : mpdata_type
use advection_MPDATA_mappings_module,  only : mappings_type
use advection_MPDATA_topology_module,  only : topology_type
use advection_MPDATA_geometry_module,  only : geometry_type
use advection_MPDATA_nabla_module,     only : nabla_type

implicit none
private

public :: create_mgrids_mpdata

!------------------------------------------------------------------------------

type, private :: mgrids_interpolation
  type(atlas_Interpolation)                   :: interpolation_to_IFS
  type(atlas_Interpolation)                   :: interpolation_from_IFS
  type(atlas_functionspace_StructuredColumns) :: fs_from_IFS
  type(atlas_functionspace_StructuredColumns) :: fs_from_adv
  type(atlas_functionspace_StructuredColumns) :: fs_to_IFS
  type(atlas_functionspace_NodeColumns)       :: fs_to_adv
  integer(jpim)                               :: ngp_IFS
  integer(jpim)                               :: ngp_adv
contains
  procedure :: to_IFS
  procedure :: from_IFS
end type

!------------------------------------------------------------------------------

type, public, extends( dwarf_mpdata ) :: mgrids_mpdata

  type(trip)     , pointer :: yrrip       => null()
  type(tdyn)     , pointer :: yrdyn       => null()
  integer :: nflevg 
  integer :: ngptot 
  integer :: nproma 
  integer :: ngpblks


  real(kind=jprb), pointer :: zuvw(:,:,:)    => null()
    ! shape=( 3, ydgeometry%yrdimv%nflevg, ydgeometry%yrgem%ngptot )
  real(kind=jprb), pointer :: zrho(:,:)      => null()
    ! shape=( ydgeometry%yrdimv%nflevg, ydgeometry%yrgem%ngptot )
  real(kind=jprb), pointer :: ztracer(:,:,:) => null()
    ! shape=( ygfl%nfmg, ydgeometry%yrdimv%nflevg, ydgeometry%yrgem%ngptot )
  real(kind=jprb), pointer :: zlevh_ifs(:,:) => null()
    ! shape=( ydgeometry%yrdimv%nflevg, ydgeometry%yrgem%ngptot )

  logical       :: lremesh
  logical       :: lusecfl  
  real(jprb)    :: config_zcflmax
  integer(jpim) :: config_nsubsteps

  type(mgrids_interpolation) :: interpolation
  type(atlas_functionspace_StructuredColumns) :: fs_adv


  type(geometry_type) :: geometry
  type(mappings_type) :: mappings
  type(topology_type) :: topology
  type(nabla_type)    :: nabla

  type(atlas_Field)   :: zlevh_field
  type(atlas_Config)  :: config

  type(atlas_Field)   :: field_IFS_rho
  type(atlas_Field)   :: field_IFS_vel
  type(atlas_Field)   :: field_IFS_tracer
  type(atlas_Field)   :: field_IFS_zlevh
  type(atlas_Field)   :: field_adv_tracer

contains
  procedure, public :: setup
  procedure, public :: preprocess
  procedure, public :: postprocess
  procedure, public :: final
end type

!------------------------------------------------------------------------------

interface mgrids_mpdata
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
  class(mgrids_mpdata), pointer :: this
  type(mgrids_advection_config)      :: config
  type(trip),     intent(in), target :: ydrip
  type(tdyn),     intent(in), target :: yddyn
  type(geometry), intent(in), target :: ydgeometry
  !-----------------------------------------------------------
  type(atlas_Trace) :: trace
  trace = atlas_Trace(__FILENAME__,__LINE__,"mgrids_mpdata::constructor")
  allocate( mgrids_mpdata :: this )
  call this%setup(                                           &
              config,                                        &
              ydrip,                                         &
              yddyn,                                         &
              ydgeometry )
  call trace%final()
end function

!------------------------------------------------------------------------------

subroutine create_mgrids_mpdata(  &
  this,                           &
  config,                         &
  ydrip,                          &
  yddyn,                          &
  ydgeometry )
  !-----------------------------------------------------------
  class(dwarf), pointer, intent(out)        :: this
  type(mgrids_advection_config), intent(in) :: config
  type(trip),     intent(in), target :: ydrip
  type(tdyn),     intent(in), target :: yddyn
  type(geometry), intent(in), target :: ydgeometry
  !-----------------------------------------------------------
  this => mgrids_mpdata( config,ydrip,yddyn,ydgeometry )
end subroutine


!------------------------------------------------------------------------------

subroutine setup_mpdata_geometry( this, gridname, config, ydgeometry)
  use fckit_module,  only : fckit_log, fckit_exception
!  use atlas_module
  !-----------------------------------------------------------
  class(mgrids_mpdata), intent(inout) :: this
  character(len=*),    intent(in)    :: gridname
  type(atlas_Config),  intent(in)    :: config
  type(geometry),      intent(in)    :: ydgeometry
  !-----------------------------------------------------------
  type(atlas_Grid)             :: grid
  type(atlas_MeshGenerator)    :: meshgenerator
  type(atlas_Partitioner)      :: partitioner
  type(atlas_GridDistribution) :: distribution
  type(atlas_Mesh)             :: mesh
  type(atlas_fvm_Method)       :: fvm
  character(len=1024)          :: string
  !-----------------------------------------------------------

  if( trim(gridname) == "same" ) then

    call fckit_log%info( "    Setup MPDATA geometry using same mesh as IFS" )

    mesh = ydgeometry%yratlas%mesh
    fvm  = ydgeometry%yratlas%fvm
    if( this%lremesh ) then
      this%fs_adv = ydgeometry%yratlas%fs_structuredcolumns
    endif

  else

    write(string,'(A,A,A)') "    Setup MPDATA geometry using grid ", trim(gridname), " ... "
    call fckit_log%info( string )

    grid          = atlas_Grid(trim(gridname))
    partitioner   = atlas_MatchingMeshPartitioner( ydgeometry%yratlas%mesh )
    distribution  = partitioner%partition(grid)
    meshgenerator = atlas_MeshGenerator()
    mesh          = meshgenerator%generate( grid, distribution )
    fvm           = atlas_FVM_method( mesh, config )

    if( this%lremesh ) then
      this%fs_adv = atlas_functionspace_StructuredColumns( grid, distribution, halo=3 )
        ! halo=2 may be enough, 3 for certainty
    endif

    call meshgenerator%final()
    call distribution%final()
    call partitioner%final()
    call grid%final()

    write(string,'(A,A,A)') "    Setup MPDATA geometry using grid ", trim(gridname), " ... done"
    call fckit_log%info( string )

  endif

  this%geometry = geometry_type( fvm, mesh, config )

  call mesh%final()
  call fvm%final()
end subroutine setup_mpdata_geometry

!------------------------------------------------------------------------------
! Invoked only once from sudyn() -->  MGRIDS_ADVECTION() -->  
!                        mgrids_advection_create() --> setup_ifs_mpdata() -->
!                        setup_ifs_mpdata() ....
!------------------------------------------------------------------------------
subroutine setup(                         &
    this,                                 &
    mgrids_config,                        &
    ydrip,                                &
    yddyn,                                &
    ydgeometry                            &
  )
  !-----------------------------------------------------------
  use fckit_module, only: fckit_log
  use yom_ygfl,          only : ygfl
  use yomcst,            only : RA
  !-----------------------------------------------------------
  class(mgrids_mpdata), intent(inout) :: this
  type(mgrids_advection_config), intent(in) :: mgrids_config
  type(trip),     intent(in), target :: ydrip
  type(tdyn),     intent(in), target :: yddyn
  type(geometry), intent(in), target :: ydgeometry
  !-----------------------------------------------------------
  type(atlas_config)     :: config
  real(jprb) :: ztop = 100000.0_jprb
  integer :: nb_levels, ngptot, nproma, ngpblks
  !-----------------------------------------------------------
  call fckit_log%info( "Setup mgrids advection using MPDATA scheme" )

  this%yrrip      => ydrip
  this%yrdyn      => yddyn
  this%nflevg  = ydgeometry%yrdimv%nflevg
  this%ngptot  = ydgeometry%yrgem%ngptot
  this%nproma  = ydgeometry%yrdim%nproma
  this%ngpblks = ydgeometry%yrdim%ngpblks

  this%ntrac = ygfl%nfmg
  nb_levels  = ydgeometry%yrdimv%nflevg
  ngptot     = ydgeometry%yrgem%ngptot
  nproma     = ydgeometry%yrdim%nproma
  ngpblks    = ydgeometry%yrdim%ngpblks

  this%lusecfl = mgrids_config%mpdata_lusecfl
  this%config_nsubsteps = mgrids_config%mpdata_nsubsteps
  this%config_zcflmax = mgrids_config%mpdata_zcflmax

  this%config = atlas_config()
  call this%config%set("nb_levels" , nb_levels )
  call this%config%set("iout"      , this%config_nsubsteps )
  call this%config%set("dz"        , ztop/(nb_levels-1) )
  call this%config%set("vstretch"  , 0 )
  call this%config%set("eps0"      , 1.0e-9_jprb )
  call this%config%set("mporder"   , 2 )
  call this%config%set("ivbz"      , -1.0_jprb )
  call this%config%set("limit"     , 1.0_jprb )
  call this%config%set("radius"    , RA )
  call this%config%set("halo"      , 2 )

  this%lremesh = mgrids_config%lremesh

  call setup_mpdata_geometry( this, mgrids_config%cgrid, this%config, ydgeometry )

  ! Fields in Nodes
  this%tracer_field = this%geometry%nodes%create_field( name="trac",   kind=atlas_real(jprb), levels=nb_levels, variables=max(1,this%ntrac) )
  this%rho_field    = this%geometry%nodes%create_field( name="rhophy", kind=atlas_real(jprb), levels=nb_levels )
  this%rhojacfield  = this%geometry%nodes%create_field( name="rhojac", kind=atlas_real(jprb), levels=nb_levels )
  this%rhofacfield  = this%geometry%nodes%create_field( name="rhofac", kind=atlas_real(jprb), levels=nb_levels )
  this%rhodtfield   = this%geometry%nodes%create_field( name="rhodt",  kind=atlas_real(jprb), levels=nb_levels )
  this%wnfield      = this%geometry%nodes%create_field( name="wn",     kind=atlas_real(jprb), levels=nb_levels+1 )
  this%velfield     = this%geometry%nodes%create_field( name="vel",    kind=atlas_real(jprb), levels=nb_levels, variables=3 )
  this%vhatfield    = this%geometry%nodes%create_field( name="vhat",   kind=atlas_real(jprb), levels=nb_levels, variables=3 )
  this%zlevh_field  = this%geometry%nodes%create_field( name="zlevh",  kind=atlas_real(jprb), levels=nb_levels )

  ! Fields in Edges
  this%vxyfield     = this%geometry%edges%create_field( name="vxy",    kind=atlas_real(jprb), levels=nb_levels, variables=2 )
  this%vnfield      = this%geometry%edges%create_field( name="vn",     kind=atlas_real(jprb), levels=nb_levels )

  if( this%lremesh ) then
    call fckit_log%info( "    Setup mgrids interpolators..." )
    this%interpolation%fs_from_IFS    = ydgeometry%yratlas%fs_structuredcolumns
    this%interpolation%fs_from_adv    = this%fs_adv ! setup during "setup_mpdata_geometry()"
    this%interpolation%fs_to_IFS      = this%interpolation%fs_from_IFS
    this%interpolation%fs_to_adv      = this%geometry%nodes ! setup during "setup_mpdata_geometry()"
    this%field_adv_tracer = this%fs_adv%create_field(name="tracer",levels=nb_levels,variables=this%ntrac,kind=atlas_real(jprb))
    config = atlas_Config()
    call config%set("type","quasicubic2D")
    this%interpolation%interpolation_from_IFS = atlas_Interpolation(config, &
      & this%interpolation%fs_from_IFS, this%interpolation%fs_to_adv)
    this%interpolation%interpolation_to_IFS = atlas_Interpolation(config, &
        & this%interpolation%fs_from_adv, this%interpolation%fs_to_IFS)
    call config%final()
    this%interpolation%ngp_IFS = ngptot
    this%interpolation%ngp_adv = this%geometry%nb_nodes
    call fckit_log%info( "    Setup mgrids interpolators... done" )
  else
    call fckit_log%info( "    No mgrids interpolators required ( lremesh = false )" )
  endif

  if( this%lremesh ) then
    this%field_IFS_vel = this%interpolation%fs_from_IFS%create_field( kind=atlas_real(jprb), levels=nb_levels, &
      name="vel", variables=3, type="vector" )

    this%field_IFS_rho = this%interpolation%fs_from_IFS%create_field( kind=atlas_real(jprb), levels=nb_levels, &
      name="rho", type="scalar" )

    this%field_IFS_tracer = this%interpolation%fs_from_IFS%create_field( kind=atlas_real(jprb), levels=nb_levels, &
      name="tracer", variables=this%ntrac, type="scalar" )

    this%field_IFS_zlevh = this%interpolation%fs_from_IFS%create_field( kind=atlas_real(jprb), levels=nb_levels, &
      name="zlevh", type="scalar" )
  else
    this%field_IFS_vel    = this%velfield
    this%field_IFS_rho    = this%rho_field
    this%field_IFS_tracer = this%tracer_field
    this%field_IFS_zlevh  = this%zlevh_field
  endif
  call this%field_IFS_vel%data( this%zuvw )
  call this%field_IFS_rho%data( this%zrho )
  call this%field_IFS_tracer%data( this%ztracer )
  call this%field_IFS_zlevh%data( this%zlevh_ifs )

  ! Build topology
  this%topology = topology_type( this%geometry, this%config )

  ! Set nabla operators and build mappings
  this%nabla    = nabla_type( this%geometry, this%config )
  this%mappings = mappings_type( this%geometry, this%nabla, this%config )
  this%jacobian_determinant = this%mappings%jacobian_determinant

  ! Create MPDATA structure  
  this%mpdata = mpdata_type( this%geometry, this%config )

  call this%dwarf_mpdata_setup()

end subroutine

!------------------------------------------------------------------------------
! Invoked at all tsteps to interpolate winds, tracers on mesh B and specify
! mpdata advection levels
!------------------------------------------------------------------------------
subroutine preprocess(this,args)
  use yom_ygfl,      only : ygfl
  use yomct3,        only : nstep
  use yomcst,        only : RG
  use intdyn_mod,    only : yytcty0, yytrcp0
  use yom_grib_codes,only : ngrbq
  !-----------------------------------------------------------
  class(mgrids_mpdata),  intent(inout) :: this
  class(dwarf_args),     intent(inout) :: args
  !-----------------------------------------------------------
  type(mgrids_advection_args) :: arguments
  integer(jpim) :: itrac, jgfl, jgflq, jkglo, icend, ibl
  integer(jpim) :: jlev, jrof, jjlev, jjrof, jnode
  real(jprb), pointer :: vel(:,:,:)
  real(jprb), pointer :: vhat(:,:,:)
  type(atlas_Trace) :: trace, trace_transposition, trace_interpolation
  !-----------------------------------------------------------
  trace = atlas_Trace(__FILENAME__,__LINE__,"mgrids_mpdata::preprocess")
  arguments = mgrids_advection_args(args)
  this%istep = nstep

  associate(                                   &
    numflds => ygfl%numflds,                   &
    nmgflds => ygfl%nmgflds,                   &
    ycomp   => ygfl%ycomp,                     &
    nflevg  => this%nflevg,  &
    ngptot  => this%ngptot,   &
    nproma  => this%nproma,   &
    ngpblks => this%ngpblks,  &
    yt0     => arguments%yrgmv%yt0,            &
    yt9     => arguments%yrgmv%yt9,            &
    pgeo    => arguments%zgeo0,                &
    pgmv    => arguments%zgmv,                 &
    pgfl    => arguments%zgfl,                 &
    pre0f   => arguments%zpre0f,               &
    pcty0   => arguments%zcty0,                &
    prcp0   => arguments%zrcp0,                &
    dz      => this%geometry%dz,               &
    nb_levels => this%geometry%nb_levels,      &
    nb_nodes  => this%geometry%nb_nodes        &
  )

! find specific humidity array index
  jgflq=1
  do while ( jgflq<=numflds .and. ycomp(jgflq)%igrbcode /= ngrbq )
    jgflq = jgflq +1 
  end do

!----------------------------------------------------------------------------------------
! Transfer ifs arrays to atlas arrays on mesh A and then interpolate to mesh B.
!----------------------------------------------------------------------------------------
  trace_transposition = atlas_Trace(__FILENAME__,__LINE__,"transposition")
  !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jkglo,icend,ibl,jrof,jlev,jjrof,jjlev,itrac)
  do jkglo=1,ngptot,nproma
    icend=min(nproma,ngptot-jkglo+1)
    ibl=(jkglo-1)/nproma+1
    do jlev=1,nflevg
      do jrof=1,icend
        jjrof = jkglo+jrof-1
        jjlev = nflevg-jlev+1
        this%zlevh_ifs(jjlev,jjrof) = pgeo(jrof,jlev,ibl)/RG
        this%zrho(jjlev,jjrof)      = pre0f(jrof,jlev,ibl)*            & 
        &                 (1.0_jprb-pgfl(jrof,jlev,ycomp(jgflq)%mp,ibl) ) /           & 
        &                 (prcp0(jrof,jlev,yytrcp0%m_r,ibl)*pgmv(jrof,jlev,yt0%mt,ibl))
        this%zuvw(1,jjlev,jjrof)    = pgmv(jrof,jlev,yt0%mu,ibl)
        this%zuvw(2,jjlev,jjrof)    = pgmv(jrof,jlev,yt0%mv,ibl)
        this%zuvw(3,jjlev,jjrof)    = -pre0f(jrof,jlev,ibl)* & 
        &    pcty0(jrof,jlev,yytcty0%m_vvel,ibl)/(RG*this%zrho(jjlev,jjrof))
        do itrac=1,this%ntrac
          this%ztracer(itrac,jjlev,jjrof) = pgfl(jrof,jlev,ycomp(nmgflds(itrac))%mp,ibl)
        enddo
      enddo
    enddo
  enddo
  !$OMP END PARALLEL DO


  call this%field_IFS_zlevh%set_dirty()
  call this%field_IFS_rho%set_dirty()
  call this%field_IFS_vel%set_dirty()
  call this%field_IFS_tracer%set_dirty()

  call trace_transposition%final()

  !-----------------------------------------------------------------------------
  ! Interpolate fields from IFS grid to Atlas mesh grid.
  !
  !   A new set of height vertical levels is created based on latest avail
  !   geopotential. The latest avail density (from IFS press) is also used. 
  !   Then the velocity and density are held at these new equivalent 
  !   height levels there is no need to do any vertical interpolation.
  !-----------------------------------------------------------------------------

  trace_interpolation = atlas_Trace(__FILENAME__,__LINE__,"interpolation")
  if (this%lremesh) then
    call this%interpolation%from_IFS(   this%field_IFS_zlevh, this%zlevh_field )
    call this%interpolation%from_IFS(   this%field_IFS_rho,   this%rho_field )
    call this%interpolation%from_IFS(   this%field_IFS_vel,   this%velfield )
    if( this%ntrac > 0 ) then
      call this%interpolation%from_IFS( this%field_IFS_tracer, this%tracer_field )
    endif
  endif
  call trace_interpolation%final()

  call this%zlevh_field%halo_exchange()
  call this%velfield%halo_exchange()
  call this%rho_field%halo_exchange()
  if( this%ntrac > 0 ) then
    call this%tracer_field%halo_exchange()
  endif

  ! -----------------------------------------------------------------------------

  ! Update topology
  call this%topology%execute( this%zlevh_field, this%geometry%topology_fields )
  
  ! Update mappings based on new topology
  call this%mappings%execute()
  call this%mappings%transform_vector( this%velfield, this%vhatfield, 1)

  ! -----------------------------------------------------------------------------

  ! Compute timestep for mpdata based on cfl
  call compute_tstep(this%config_nsubsteps, this%yrrip%tstep, dz, this%nsteps, this%dt )


  end associate
  call trace%final()
contains

  subroutine compute_tstep( nsubdt_mpdata, dt_ifs, dz, nsteps, dt_mpdata )
    use fckit_module, only : fckit_log
    integer(jpim), intent(in)  :: nsubdt_mpdata
    real(jprb),    intent(in)  :: dt_ifs
    real(jprb),    intent(in)  :: dz
    integer(jpim), intent(out) :: nsteps
    real(jprb),  intent(out) :: dt_mpdata
    !-----------------------------------------------------------------------------
    type(atlas_Field) :: cflw_fld, cflwmax_fld, loc_fld
    integer(jpim), pointer :: loc(:)
    real(jprb), pointer :: cflw(:,:)
    real(jprb), pointer :: cflwmax(:)
    real(jprb) :: cflwmax_all
    character(len=1024) :: string
    integer(jpim) :: jmaxnode, jlmax
    integer(jpim) :: idiag=0
    real(jprb) :: zdtrdz
    type(atlas_Trace) :: trace
    !-----------------------------------------------------------------------------
    trace = atlas_Trace(__FILENAME__,__LINE__,"compute_tstep")
    nsteps = nsubdt_mpdata

    if( this%lusecfl ) then

    associate(                                   &
      nb_levels => this%geometry%nb_levels,      &
      nb_nodes  => this%geometry%nb_nodes        &
    )

    if (idiag>0) then 
      cflwmax_fld = atlas_field(name="cflwmax",kind=atlas_real(jprb),shape=[nb_levels])
      loc_fld = atlas_field(name="loc",kind=atlas_integer(jpib),shape=[nb_levels])
      call cflwmax_fld%data(cflwmax) 
      call loc_fld%data(loc)
    endif
    cflw_fld = this%geometry%nodes%create_field(name="cflw",kind=atlas_real(jprb),levels=nb_levels)
    call cflw_fld%data(cflw)

    call this%vhatfield%data(vhat)

    zdtrdz=dt_ifs/dz
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JNODE,JLEV)
    do jnode=1, nb_nodes
      do jlev=1, nb_levels
        cflw(jlev,jnode)=abs(vhat(3,jlev,jnode))*zdtrdz
      enddo
    enddo
    !$OMP END PARALLEL DO

    if (idiag>0) then
      call this%geometry%nodes%maximum_and_location_per_level(cflw_fld,cflwmax_fld,loc_fld)
      cflwmax_all=0.0_jprb
      do jlev=1, nb_levels
        write(string,*) '#cflw max at lev,',jlev,'is:',cflwmax(jlev),loc(jlev)
        call fckit_log%info(string)
        cflwmax_all=max(cflwmax_all,cflwmax(jlev))
      enddo
    else
      call this%geometry%nodes%maximum(cflw_fld,cflwmax_all)
    endif

    nsteps=max(nsubdt_mpdata,ceiling(cflwmax_all/this%config_zcflmax))
    dt_mpdata=dt_ifs/nsteps
    
    write(string,*) '#MPDATA substeps, cflw max is:',nsteps,cflwmax_all/nsteps
    call fckit_log%info(string)

    call cflw_fld%final()
    if (idiag>0) THEN
      call cflwmax_fld%final()
      call loc_fld%final()
    endif

    end associate

    endif
    call trace%final()
  end subroutine compute_tstep
end subroutine preprocess

!------------------------------------------------------------------------------

subroutine postprocess(this,args)
  use yom_ygfl, only : ygfl
  !-----------------------------------------------------------
  class(mgrids_mpdata), intent(inout) :: this
  class(dwarf_args),    intent(inout) :: args
  !-----------------------------------------------------------
  integer(jpim) :: itrac, jgfl, jkglo, icend, ibl, jlev, jrof, jnode
  real(jprb), pointer :: tracer(:,:,:), adv_tracer(:,:,:)
  type(mgrids_advection_args) :: arguments
  type(atlas_Trace) :: trace, trace_interpolation, trace_transposition
  !-----------------------------------------------------------
  trace = atlas_Trace(__FILENAME__,__LINE__,"mgrids_mpdata::postprocess")
  arguments = mgrids_advection_args(args)

  associate(                                     &
    nmgflds   => ygfl%nmgflds,                   &
    ycomp     => ygfl%ycomp,                     &
    nflevg    => this%nflevg,                    &
    ngptot    => this%ngptot,                    &
    nproma    => this%nproma,                    &
    nb_levels => this%geometry%nb_levels,        &
    nb_nodes  => this%geometry%nb_nodes,         &
    pgflt1    => arguments%zgflt1                &
  )

  if (this%lremesh) then
    trace_interpolation = atlas_Trace(__FILENAME__,__LINE__,"interpolation")
    if( this%ntrac > 0 ) then
      call this%tracer_field%data( tracer )
      call this%field_adv_tracer%data( adv_tracer )
      !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode)
      do jnode=1,this%interpolation%ngp_adv
        adv_tracer(:,:,jnode) = tracer(:,:,jnode)
      enddo
      !$OMP END PARALLEL DO 
      call this%field_adv_tracer%set_dirty()
      call this%interpolation%to_IFS( this%field_adv_tracer, this%field_IFS_tracer )
    endif
    call trace_interpolation%final()
  else
    ! this%field_IFS_tracer is the same as this%tracer_field
    !   with its data already pointed to by this%ztracer
  endif

!-----------------------------------------------------------------------------
! interpolate the tracer at the newly specified levels.
!-----------------------------------------------------------------------------
  trace_transposition = atlas_Trace(__FILENAME__,__LINE__,"transposition")

  !$OMP PARALLEL DO SCHEDULE(STATIC) &
  !$OMP& PRIVATE(jkglo,icend,ibl,jrof,jlev,itrac)
  do jkglo=1,ngptot,nproma
    icend=min(nproma,ngptot-jkglo+1)
    ibl=(jkglo-1)/nproma+1
    do jlev=1,nflevg
      do jrof=1,icend
        do itrac=1,this%ntrac
          pgflt1(jrof,jlev,ycomp(nmgflds(itrac))%mp1,ibl) = this%ztracer(itrac,nflevg-jlev+1,jkglo+jrof-1)
        enddo
      enddo
    enddo
  enddo
  !$OMP END PARALLEL DO
  call trace_transposition%final()
  end associate

  call trace%final()
end subroutine postprocess

!------------------------------------------------------------------------------

subroutine final(this)
  class(mgrids_mpdata), intent(inout) :: this

  if( this%lremesh ) then

    call this%interpolation%interpolation_from_IFS%final()
    call this%interpolation%interpolation_to_IFS%final()
    call this%interpolation%fs_to_IFS%final()
    call this%interpolation%fs_to_adv%final()
    call this%interpolation%fs_from_IFS%final()
    call this%interpolation%fs_from_adv%final()

  endif

  call this%field_IFS_rho%final()
  call this%field_IFS_vel%final()
  call this%field_IFS_tracer%final()
  call this%field_IFS_zlevh%final()
  call this%field_adv_tracer%final()

  call this%config%final()
  call this%vhatfield%final()
  call this%zlevh_field%final()
  ! TODO call this%nabla%final()
  call this%mappings%final()
  ! TODO call this%topology%final()
  call this%geometry%final()

  call this%dwarf_mpdata_final()

end subroutine final

!------------------------------------------------------------------------------

subroutine from_IFS(this,field_IFS,field_adv)
  class(mgrids_interpolation), intent(inout) :: this
  type(atlas_Field), intent(in)    :: field_IFS
  type(atlas_Field), intent(inout) :: field_adv
  call field_IFS%halo_exchange()
  call this%interpolation_from_IFS%execute( field_IFS, field_adv )
  call field_adv%halo_exchange()
end subroutine

!------------------------------------------------------------------------------

subroutine to_IFS(this,field_adv,field_IFS)
  class(mgrids_interpolation), intent(inout) :: this
  type(atlas_Field), intent(in)    :: field_adv
  type(atlas_Field), intent(inout) :: field_IFS
  call field_adv%halo_exchange()
  call this%interpolation_to_IFS%execute( field_adv, field_IFS )
  call field_IFS%halo_exchange()
end subroutine

!------------------------------------------------------------------------------

#endif
end module
