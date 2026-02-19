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

module spp_gen_mod

USE PARKIND1 , ONLY : JPIM, JPRB

implicit none

integer(kind=jpim), parameter :: jp_lab_len=16      ! string length for perturbation labels
integer(kind=jpim), parameter :: jp_ver_len=32      ! string length for version
integer(kind=jpim), parameter :: jp_mag  =4         ! max size of xmag
integer(kind=jpim), parameter :: jp_off  =2         ! max size of nseed_off
integer(kind=jpim), parameter :: jpmaxperts=128     ! max number of perturbations in the context of namelist inputs
character(len=jp_lab_len), parameter :: cp_undefined='????????????????'

type spp_pert
  !
  !   define variables for an individual perturbation of 1 variable or 1 parameter
  !
  character(len=jp_lab_len) :: label=cp_undefined
  logical                   :: on=.false.       ! on/off
  integer(kind=jpim)        :: idistr           ! Determines the distribution; 0: normal, 1: log-normal, 2: ...
  real(kind=jprb), dimension(jp_mag)   :: xmag=0._jprb        ! amplitude
  real(kind=jprb), dimension(jp_mag)   :: mu=0._jprb          ! mean of Gaussian
  integer(kind=jpim)        :: nmag=1           !   number of used magnitudes
  logical                   :: ln1=.true.       ! T implies the mean of the log-normal distribution is 1
  integer(kind=jpim),dimension(jp_off) :: nseed_off=0       ! random number seed offsets
  integer(kind=jpim)        :: noff=1           !   number of used seed offsets

  integer(kind=jpim)        :: nrf=1            ! number of random fields for this perturbation
                                                !   [e.b. would be 2 for CUDUDV but 1 for RTAU and most others]
  integer(kind=jpim)        :: nrf_radgrid=0    ! number of random fields on the radiation grid
  logical                   :: radgrid=.false.  ! perturbation requires random field on radiation grid
  character(len=jp_lab_len) :: rf_pert_label=cp_undefined ! point to random field of perturbation rf_pert_label
  integer(kind=jpim)        :: mp               ! pointer to first random field
                                                !   [correspond to MP* IN TSPP_DATA, e.g. MPCFM ... MPZHSVDAERO]
  integer(kind=jpim)        :: mp_radgrid       ! pointer to first random field on radiation grid

  real(kind=jprb), dimension(2)        :: xclip
  real(kind=jprb)                      :: tau, sdev, xlcor ! random field characteristics
end type spp_pert

type spp_model
  !
  !    define a model of perturbations
  !
  character(len=jp_ver_len)  :: version
  integer(kind=jpim)         :: nmax  ! max. number of perturbations for array dimensioning
  character(len=jp_lab_len), dimension(:), allocatable :: defined_perts
  character(len=jp_lab_len), dimension(:), allocatable :: active_perts
  integer(kind=jpim) :: ndef=0 ! number of defined perturbations
  integer(kind=jpim) :: nact ! number of active perturbations
  integer(kind=jpim) :: nrftotal=0 ! total number of random fields
  integer(kind=jpim) :: nrftotal_radgrid=0 ! total number of random fields on the radiation grid
  integer(kind=jpim), dimension(:), allocatable :: nseed_off ! seed offsets for random fields 1..nrftotal
  type(spp_pert), dimension(:), allocatable :: pndef  ! configurations of defined individual perturbations
  type(spp_pert), dimension(:), allocatable :: pn   ! configurations of active individual perturbations
  real(kind=jprb)    :: tau    ! global default decorrelation time
  real(kind=jprb)    :: xlcor  ! global default correlation length scale
  real(kind=jprb)    :: sdev   ! global default standard deviation
  integer(kind=jpim) :: kseed_off ! random number seed offset for current perturbation
end type spp_model


contains

subroutine allocate_spp_model( sm, knmax)
  !
  !   allocate all dynamic arrays in the spp_model variable ydspp
  !
  type(spp_model), intent(inout) :: sm
  integer(kind=jpim), intent(in) :: knmax

  sm%nmax = knmax
  allocate(sm%defined_perts( knmax ) )
  allocate(sm%active_perts(  knmax ) )
  allocate(sm%pndef(         knmax ) )
  allocate(sm%pn(            knmax ) )
  allocate(sm%nseed_off(     knmax * jp_off ) )

end subroutine allocate_spp_model

subroutine deallocate_spp_model( sm )
  !
  !   deallocate all arrays in the spp_model variable ydspp
  !
  type(spp_model), intent(inout) :: sm

  deallocate(sm%defined_perts )
  deallocate(sm%active_perts )
  deallocate(sm%pndef )
  deallocate(sm%pn    )
  deallocate(sm%nseed_off)

end subroutine deallocate_spp_model

subroutine update_mu_spp_pertn( pn )
  !
  !   if distribution is log-normal and ln1==T, compute mean mu(.) of Gaussian distribution that
  !     yields mean of log-normal distribution equal to one
  !
  type(spp_pert), intent(inout) :: pn

  integer(kind=jpim) :: jmag

  if (pn%ln1.and.(pn%idistr==1)) then
    do jmag=1, pn%nmag
      pn%mu( jmag ) = -0.5_jprb * ( pn%xmag( jmag ) * pn%sdev )**2
    enddo
  else
    pn%mu= 0._jprb
  endif

end subroutine update_mu_spp_pertn

subroutine implement_spp_pertn( sm, cdlabel, pxmag, kidistr, knseed_off, knrf, ln1, radgrid, cd_rf_pert_label  )
  !
  !   configures a new defined perturbation in spp_model sm component.
  !    The subroutine increments sm%ndef and stores the configuration of the perturbation in sm%pndef( sm%ndef )
  !
  type(spp_model),           intent(inout) :: sm
  character(len=*),          intent(in)    :: cdlabel
  real(kind=jprb), dimension(:)        , intent(in)    :: pxmag
  ! optional arguments
  integer(kind=jpim), optional        , intent(in)    :: kidistr
  integer(kind=jpim), dimension(:), optional        , intent(in)    :: knseed_off
  integer(kind=jpim), optional        , intent(in)    :: knrf
  logical,            optional        , intent(in)    :: ln1
  logical,            optional        , intent(in)    :: radgrid
  character(len=*),   optional        , intent(in)    :: cd_rf_pert_label

  type(spp_pert) :: pn

  pn%label      = cdlabel    ! who I am
  pn%on         = .false.    ! default is off

  !
  !   inherit global settings for random field
  !
  pn%tau    = sm%tau
  pn%sdev   = sm%sdev
  pn%xlcor  = sm%xlcor

  pn%nmag = size(pxmag)
  if (pn%nmag > jp_mag) then
    call abor1('It is not possible to request more than jp_mag different magnitudes. Change jp_mag in spp_test_mod. ')
  endif
  pn%xmag(1:pn%nmag)= pxmag(:)

  if (present(kidistr)) then
    pn%idistr = kidistr
  else
    pn%idistr = 1
  endif

  if (present(ln1)) then
    pn%ln1 = ln1  ! default set in type declaration
  endif
  call update_mu_spp_pertn(pn)

  !
  !  random field required on radiation grid?
  !
  if (present(radgrid)) then
    pn%radgrid= radgrid
  endif
  !
  !  ML: The number of random fields nrf could be determined from size(nseed_off)
  !
  if (present(knrf)) then
    pn%nrf = knrf
  else
    pn%nrf = 1
  endif
  if (pn%radgrid) then
    pn%nrf_radgrid = pn%nrf
  else
    pn%nrf_radgrid = 0
  endif
  !
  !
  !
  if (present(cd_rf_pert_label)) then
    if (present(knrf).or.present(knseed_off)) then
      call abor1('Neither knrf nor knseed_off are meaningful arguments when cd_rf_pert_label is present.')
    endif
    pn%rf_pert_label = cd_rf_pert_label
    pn%nrf         = 0
    pn%nrf_radgrid = 0
  endif

  !
  !   potentially use hash function like CRC-32 to generate seed modifier
  !
  if (present(knseed_off)) then
    pn%noff=size(knseed_off)
    if (pn%noff > jp_off) then
      call abor1('Too many seed offsets; increase jp_off in spp_test_mod.')
    endif
    if (pn%noff /= pn%nrf ) then
      call abor1('The number of random fields must be equal to the number of random seeds.')
    endif
    pn%nseed_off(1:pn%noff) = knseed_off
  else if (present(cd_rf_pert_label)) then
    pn%nseed_off = 0
  else
    pn%nseed_off = sm%kseed_off
  endif
  !
  !   update SPP model
  !
  if ( (sm%ndef +1) > sm%nmax) then
    call abor1( 'nmax is too small.')
  endif
  sm%ndef = sm%ndef + 1
  sm%pndef( sm%ndef ) = pn
  sm%defined_perts(sm%ndef) = cdlabel
  !
  !   seed offsets for next perturbation if called without nseed_off arguments
  !
  if (sm%kseed_off < 2_jpim **30) then
    sm%kseed_off = 2*sm%kseed_off
  else
    sm%kseed_off = 2_jpim **8 + 2_jpim **9
  endif
end subroutine implement_spp_pertn





subroutine get_active_spp_perts( sm, ku_nam, ku_out )
  !
  !   read namelist input from unit kunit to determine active perturbations and make adjustments to their configuration
  !   create list of seeds for all random fields
  !
  type(spp_model), intent(inout) :: sm
  integer(kind=jpim), intent(in) :: ku_nam
  integer(kind=jpim), intent(in) :: ku_out

  integer(kind=jpim) :: inreq ! number of requested perturbations
  integer(kind=jpim) :: jjpert, kf0, kf1

  character(len=jp_lab_len), dimension(jpmaxperts) :: requested_perts

  namelist /nam_spp_active/ requested_perts

  if (sm%nmax>jpmaxperts) then
    write(ku_out,'(A,I5)') 'spp_gen_mod:get_active_spp_perts: increase jpmaxperts to at least ',sm%nmax
    call abor1('spp_gen_mod:get_active_spp_perts: array requested_perts is too small.')
  endif

  requested_perts= cp_undefined
  read(ku_nam,nml=nam_spp_active)

  inreq = count(requested_perts /= cp_undefined )
  write(ku_out,*) inreq," perturbations have been requested"
  sm%nact = inreq

  !
  !   check that requested perturbations are actually implemented
  !
  do jjpert=1, sm%nact
    if ( .not.any(sm%defined_perts == requested_perts(jjpert) ) ) then
      ! abort with suitable error message
      call abor1( "the requested perturbation "//trim(requested_perts(jjpert))//" is not known in "//trim( sm%version ))
    else
      sm%active_perts(jjpert)=requested_perts(jjpert)
      sm%pn(jjpert) = sm%pndef( get_spp_idx_def(requested_perts(jjpert),sm) )
      sm%pn(jjpert)%on  =.true.
    endif
  enddo
end subroutine get_active_spp_perts


subroutine map_indices_spp_nml( sm, ku_nam, ku_out, kidx_map )
  !
  !   need to map perturbations in nam_spp_modif to active perturbations in sm
  !
  type(spp_model), intent(inout) :: sm
  integer(kind=jpim), intent(in) :: ku_nam
  integer(kind=jpim), intent(in) :: ku_out
  integer(kind=jpim), dimension(jpmaxperts), intent(out):: kidx_map

  integer(kind=jpim) :: jjpert
  type(spp_pert), dimension(jpmaxperts) :: pnx   ! configurations of active individual perturbations

  namelist /nam_spp_modif/ pnx

  kidx_map = -1
  if (sm%nact>jpmaxperts) then
    write(ku_out,'(A,I5)') 'spp_gen_mod:modify_spp_perts: increase jpmaxperts to at least ',sm%nact
    call abor1('spp_gen_mod:map_indices_spp_nml: array pnx is too small.')
  endif

  pnx(:)%label = cp_undefined
  read(ku_nam,nml=nam_spp_modif)

  do jjpert=1, sm%nact
    if (pnx(jjpert)%label /= cp_undefined) then
      kidx_map(jjpert) = get_spp_idx(  pnx(jjpert)%label, sm )
    endif
  enddo

end subroutine map_indices_spp_nml


subroutine modify_spp_perts( sm, kidx_map, ku_nam, ku_out )
  !
  !   Read namelist input from unit kunit to determine active perturbations and make adjustments to their configuration.
  !
  !   kidx_map is obtained from initial pass through namelist nam_spp_modif and links the perturbations pnx with the perturbations sm%pn
  !
  !   create list of seeds for all random fields
  !
  type(spp_model), intent(inout) :: sm
  integer(kind=jpim), dimension(jpmaxperts), intent(in):: kidx_map
  integer(kind=jpim), intent(in) :: ku_nam
  integer(kind=jpim), intent(in) :: ku_out

  integer(kind=jpim) :: jjpert, kf0, kf1

  integer(kind=jpim) :: i_rf_donor

  type(spp_pert), dimension(jpmaxperts) :: pnx   ! configurations of active individual perturbations

  namelist /nam_spp_modif/ pnx

  if (sm%nact>jpmaxperts) then
    write(ku_out,'(A,I5)') 'spp_gen_mod:modify_spp_perts: increase jpmaxperts to at least ',sm%nact
    call abor1('spp_gen_mod:modify_spp_perts: array pnx is too small.')
  endif

  pnx(:)%label=cp_undefined

  do jjpert=1, jpmaxperts
    if (kidx_map(jjpert)>0) then
      pnx(jjpert) = sm%pn( kidx_map(jjpert) )
    endif
  enddo

  !
  !   make requested changes to settings
  !
  read(ku_nam,nml=nam_spp_modif)

  do jjpert=1, jpmaxperts
    if (kidx_map(jjpert)>0) then
      sm%pn( kidx_map(jjpert) ) = pnx( jjpert)
      call update_mu_spp_pertn(sm%pn( kidx_map(jjpert) ))
    endif
  enddo
  !
  !   calculate total number of random fields and
  !     pointers mp to random fields for individual perturbations
  !
  sm%nrftotal=sum( sm%pn(1:sm%nact)%nrf )
  sm%nrftotal_radgrid=sum( sm%pn(1:sm%nact)%nrf_radgrid )
  do jjpert=1, sm%nact
    if (sm%pn(jjpert)%nrf > 0) then
      sm%pn(jjpert)%mp         = 1 + sum( sm%pn(1:jjpert-1)%nrf )
      sm%pn(jjpert)%mp_radgrid = 1 + sum( sm%pn(1:jjpert-1)%nrf_radgrid )
    else
      sm%pn(jjpert)%mp         = -1
      sm%pn(jjpert)%mp_radgrid = -1
    endif
  enddo
  !
  !   pointers to random fields that are taken from another perturbation
  !
  do jjpert=1, sm%nact
    if (sm%pn(jjpert)%nrf == 0) then
      i_rf_donor = get_spp_idx(  sm%pn(jjpert)%rf_pert_label, sm )
      if (i_rf_donor < 0) then
        write(ku_out, *) 'jjpert = ',jjpert,', sm%pn(jjpert)%label = ',sm%pn(jjpert)%label
        write(ku_out, *) 'sm%pn(jjpert)%rf_pert_label = ',sm%pn(jjpert)%rf_pert_label
        call abor1('Cannot find sm%pn(jjpert)%rf_pert_label in active perturbations.')
      endif
      sm%pn(jjpert)%mp         = sm%pn(i_rf_donor)%mp
      sm%pn(jjpert)%mp_radgrid = sm%pn(i_rf_donor)%mp_radgrid
    endif
  enddo
  !
  !   create list of seed offsets for the random fields
  !
  kf1=0
  do jjpert=1, sm%nact
    if (sm%pn(jjpert)%nrf > 0) then
      kf0=kf1+1
      kf1=kf0+sm%pn(jjpert)%noff-1 !ML: the number of seeds does not need to match the number of random fields. This could cause problems. How to resolve this?
      sm%nseed_off(kf0:kf1) = sm%pn(jjpert)%nseed_off(1:sm%pn(jjpert)%noff )
    endif
  enddo
end subroutine modify_spp_perts

integer(kind=jpim) function get_spp_idx( cd_label, sm )
  !
  !   get index from active perturbations for perturbation named cd_label in spp_model sm
  !
  character(len=*), intent(in)      :: cd_label
  type(spp_model)          , intent(in)      :: sm

  integer(kind=jpim) :: ipn
  integer(kind=jpim) :: jj
  ipn=-1
  do jj=1, sm%nact
    if (cd_label == sm%pn(jj)%label) then
      ipn=jj
      exit  ! no further search needed
    endif
  enddo
  if (ipn == -1) then
    !
    !   this should never happen, abort immediately
    !
    call abor1( "get_spp_idx: there is no perturbation " //trim(cd_label))
  endif
  get_spp_idx=ipn
end function get_spp_idx

integer(kind=jpim) function get_spp_idx_def( cd_label, sm )
  !
  !   get index from defined perturbations for perturbation named cd_label
  !
  character(len=*), intent(in)      :: cd_label
  type(spp_model)          , intent(in)      :: sm

  integer(kind=jpim) :: ipn
  integer(kind=jpim) :: jj
  ipn=-1
  do jj=1, sm%ndef
    if (cd_label == sm%pndef(jj)%label) then
      ipn=jj
      exit  ! no further search needed
    endif
  enddo
  if (ipn == -1) then
    call abor1( "get_spp_idx_def: there is no perturbation "//trim(cd_label))
  endif
  get_spp_idx_def=ipn
end function get_spp_idx_def


logical function active_spp_pn(cd_label, sm )
  !
  !   return true if cd_label refers to an active_perturbation
  !
  character(len=*), intent(in)      :: cd_label
  type(spp_model)          , intent(in)      :: sm

  active_spp_pn = any( cd_label==sm%active_perts )
end function active_spp_pn


subroutine write_spp_pn_infoline( p, ku)
  !
  ! write a line with information about perturbation p to unit u
  !
  type(spp_pert),               intent(in) :: p
  integer(kind=jpim), optional, intent(in) :: ku

  integer(kind=jpim) :: iunit
  if (present(ku)) then
    iunit=ku
  else
    iunit=6
  endif
  write( iunit, '(A16,X,4F10.5,X,L3,X,I3,X,I3,X,I3,X,2I11)' ) p%label, p%xmag(:), p%ln1, &
       & p%mp, p%mp_radgrid, p%nrf, p%nseed_off(1), p%nseed_off(2)

end subroutine write_spp_pn_infoline

subroutine write_spp_model_table( sm, ku, ldefined )
  !
  !   write a table with the active perturbations (default behaviour) or
  !   the defined perturbations if ldefined=T
  !
  type(spp_model), intent(in) :: sm
  integer(kind=jpim), optional, intent(in) :: ku
  logical, optional, intent(in) :: ldefined

  logical :: lldef

  integer(kind=jpim) :: jjpert, nperts, iunit
  if (present(ku)) then
    iunit=ku
  else
    iunit=6
  endif

  if ( present(ldefined)) then
    lldef=ldefined
  else
    lldef=.false.
  endif

  if (lldef) then
    nperts = sm%ndef
  else
    nperts = sm%nact
  endif
  write( iunit, '("--------------------------------------------------------------------------------------------")' )
  write( iunit, '(A16,X,4A10,X,A3,X,A3,X,A3,X,A3,X,2A11)' ) &
       & "Perturbation    ", "Ampl 1", "Ampl 2", "Ampl 3", "Ampl 4", "ln1","Ptr","Prg", "NRF", "Seed 1", "Seed 2"

  write( iunit, '("--------------------------------------------------------------------------------------------")' )
  do jjpert=1, nperts
    if (lldef) then
      call write_spp_pn_infoline( sm%pndef(jjpert), ku=iunit)
    else
      call write_spp_pn_infoline( sm%pn(jjpert), ku=iunit)
    endif
  enddo
  write( iunit, '("--------------------------------------------------------------------------------------------")' )

end subroutine write_spp_model_table



end module spp_gen_mod
