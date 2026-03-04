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

module nvar_class
!--------------------------------------------------------------------------------------------
!   NVAR is an index into the Jo table and into VarQC settings: most VARNOs map onto an NVAR

!      Only if there is a mapping from VARNO to NVAR can an observation be processed in 4D-Var.

!      Previously NVAR_ was a static and hardcoded set of indices, prone to
!      terrible merge problems. It is now dynamically allocated at runtime, but the 
!      resulting NVAR numbering will be static as long as this code does not change.

!      NVAR is used to index the Jo-table and the VarQC settings (including Huber); at
!      Meteo-France it is the numbering scheme used for the NOTVAR mechanism to remove
!      observations from the processing (a sort of instant blacklist, see hop.F90)

!      The NVAR-associated character code is used in only two places: (i) in the Jo table 
!      printouts; (ii) to select the appropriate obs operator in ppnew.F90 and similar 
!      conventional operators.

!      To add a new one:
!         (**1) Add the NVAR definition of your choice 
!         (**2) Add a call to register_nvar

!      Why do we not use VARNO directly?

!        (a) There is no 1-to-1 mapping: for example U and V wind components both map
!            onto NVAR%U
!        (b) NVAR is a more compact range of numbers and reduces the size of arrays 
!            in Jo (JPXVAR << NVNUMAX)

!      There is a strong case for removing NVAR altogether (and using VARNOs directly), but 
!      there are the following issues:

!        (i) Meteo France use the NOTVAR function routinely in operations, assuming a
!            fixed set of numerical values (i.e. hardcoded) for each NVAR. Any replacement
!            needs significant consultation or adaptation.
!       (ii) Major improvements are needed in setting up the Huber norm (see defrun.F90) which
!            is indexed on NVAR.
!      (iii) A transparent alternative way of collapsing U and V onto one Jo variable 
!            is needed.
!       (iv) Use of the NVAR direct pointers (e.g. NVAR%V10) needs to be eliminated

!      Possible issues in 48r1:

!         - The CAMS variable GHG3 (N2O) was previously unmapped and probably had no effect
!           in data assimilation (not procesed in hjo.F90)
!         - NVARs CLW, S0 and X were not mapped and hence probably unusable in 4D-Var. There is
!           no point mapping them to a VARNO unless they are needed.
!         - After Q2/TD2, which previously both mapped to the same numerical NVAR (50) due to a
!           merge bug in CY47R1, NVAR numbering between ECMWF and Meteo-France has probably been
!           inconsistent. This does not affect ECMWF, which does not use hardcoded NVARs, but
!           may affect Meteo-France through their use of NOTVAR.
        
!   Author.
!   -------
!   Alan Geer 26-Feb-2021  Replaced code in yomcosjo.F90, map_varno_to_nvar.F90, fgchk.F90 etc.

!   Modifications.
!   --------------

!--------------------------------------------------------------------------------------------

use parkind1,     only : jpim, jplm, jprb
use yomhook,      only : lhook, dr_hook, jphook
use yomancs,      only : nmdi   
use yomlun,       only : nulout, nulerr
use varno_module, only : varno

implicit none

integer(kind=jpim), private, parameter :: max_varnos = 180

type class_nvar

  ! These variables should all be protected (it's part of F2003 standard) but 
  ! that's not supported on ECMWF compilers

  ! protected

  ! (**1) Define all NVARS
  integer(kind=jpim) :: u   = -1
  integer(kind=jpim) :: u10 = -1
  integer(kind=jpim) :: dd  = -1
  integer(kind=jpim) :: ff  = -1
  integer(kind=jpim) :: h   = -1
  integer(kind=jpim) :: h2  = -1
  integer(kind=jpim) :: t   = -1
  integer(kind=jpim) :: z   = -1
  integer(kind=jpim) :: dz  = -1
  integer(kind=jpim) :: lh  = -1
  integer(kind=jpim) :: t2  = -1
  integer(kind=jpim) :: ts  = -1
  integer(kind=jpim) :: rad = -1
  integer(kind=jpim) :: sn  = -1
  integer(kind=jpim) :: rr  = -1
  integer(kind=jpim) :: ps  = -1
  integer(kind=jpim) :: cc  = -1
  !integer(kind=jpim) :: clw = -1 ! not yet mapped to a varno
  integer(kind=jpim) :: q   = -1
  integer(kind=jpim) :: ffs = -1
  !integer(kind=jpim) :: s0 = -1 ! not yet mapped to a varno
  !integer(kind=jpim) :: x = -1  ! not yet mapped to a varno
  integer(kind=jpim) :: pwc = -1
  integer(kind=jpim) :: to3 = -1
  integer(kind=jpim) :: tcw = -1 ! Cloud water of a layer, or the total column
  integer(kind=jpim) :: rfl = -1
  integer(kind=jpim) :: apd = -1
  integer(kind=jpim) :: ro  = -1
  integer(kind=jpim) :: hls = -1
  integer(kind=jpim) :: aod = -1 ! aerosol optical depth at 0.55 microns
  integer(kind=jpim) :: lra = -1
  integer(kind=jpim), allocatable :: chem(:)
  integer(kind=jpim), allocatable :: ghg(:)
  integer(kind=jpim) :: cod = -1   ! cloud optical depth
  integer(kind=jpim) :: rao = -1   ! ratio of fine-mode aerosol optical depth to total at 0.55 mic
  integer(kind=jpim) :: rfa = -1   ! aerosol reflectance
  integer(kind=jpim) :: oda = -1   ! aerosol optical depth (multi-channel)
  integer(kind=jpim) :: dow = -1 
  integer(kind=jpim) :: sm  = -1   ! soil moisture
  integer(kind=jpim) :: smn = -1   ! normalised soil moisture
  integer(kind=jpim) :: sfr = -1   ! l-band brightness temperature
  integer(kind=jpim) :: lrr = -1   ! log transform of precipitation
  integer(kind=jpim) :: lbs = -1   ! aerosol lidar backscatter
  integer(kind=jpim) :: q2  = -1   ! 2m specific humidity
  integer(kind=jpim) :: td2 = -1   ! 2m td
  integer(kind=jpim) :: n   = -1   ! n 
  integer(kind=jpim) :: pmsl  = -1 ! mean sea-level pressure  
  integer(kind=jpim) :: ptend = -1 ! pressure tendency 
  integer(kind=jpim) :: du = -1    ! wind shear 
  integer(kind=jpim) :: dv = -1    ! wind shear 
  integer(kind=jpim) :: vt = -1    ! virtual temperature  
  integer(kind=jpim) :: td = -1    ! upper air dew point  
  integer(kind=jpim) :: w  = -1    ! past weather  
  integer(kind=jpim) :: ww = -1    ! present weather  
  integer(kind=jpim) :: vsp  = -1  ! verticalspeed  
  integer(kind=jpim) :: tsts = -1  ! sea water temperature   
  integer(kind=jpim) :: cbsc = -1  ! cloud lidar backscatter
  integer(kind=jpim) :: crfl = -1  ! cloud radar reflectivity
  integer(kind=jpim) :: lab  = -1  ! aerosol lidar attenuated backscatter

  integer(kind=jpim) :: jpxvar 

  ! Private variables for mappings between varno and nvar
  integer(kind=jpim), private :: nvar(max_varnos)      ! 1:NVARNO_MAPPINGS
  integer(kind=jpim), private :: varnum(max_varnos)    ! 1:NVARNO_MAPPINGS
  character(len=10),  private :: cvar_name(max_varnos) ! 1:JPXVAR
  integer(kind=jpim), private :: nvarno_mappings

contains

  ! Constructor
  procedure :: create         

  procedure :: map_varno
  procedure :: get_name

  procedure, private :: register_nvar_0d
  procedure, private :: register_nvar_1d

  procedure :: finish    ! manual destructor
#if !defined (__GFORTRAN__)
  final :: auto_finish   ! automatic destructor for the class, not available on ancient GFORTRAN versions
#endif

end type

#include "abor1.intfb.h"

contains

!------------------------------------------------------------------------
!
! Defines the NVAR-VARNO mappings
!
subroutine create(this,kpchem_assim,kpghg_assim)
  class(class_nvar), intent(inout) :: this
  integer(kind=jpim), intent(in) :: kpchem_assim, kpghg_assim

  integer(kind=jpim) :: i

  real(kind=jphook) :: zhook_handle
  if (lhook) call dr_hook('class_nvar%create',0,zhook_handle)

  this%nvarno_mappings = 0_jpim
  this%jpxvar          = 0_jpim

  ! (**2) Define mappings between NVARS and VARNOS
  ! ----------------------------------------------
  ! The order of these calls defines the position in the Jo table and the
  ! numerical value of the NVAR (which is unfortunately used in the NOTVAR table)

  ! This maps both u and v varnos onto the u position in the Jo table
  call this%register_nvar_0d(this%u,     'U         ', varno%u)
  call this%register_nvar_0d(this%u,     'V         ', varno%v,     ld_permit_duplicate=.true.)
  call this%register_nvar_0d(this%u10,   'U10       ', varno%u10m)
  call this%register_nvar_0d(this%u10,   'V10       ', varno%v10m,  ld_permit_duplicate=.true.)
  call this%register_nvar_0d(this%u10,   'U10       ', varno%scatu, ld_permit_duplicate=.true.)
  call this%register_nvar_0d(this%u10,   'V10       ', varno%scatv, ld_permit_duplicate=.true.)
  call this%register_nvar_0d(this%dd,    'DD        ', varno%dd)
  call this%register_nvar_0d(this%ff,    'FF        ', varno%ff)
  call this%register_nvar_0d(this%h,     'H         ', varno%rh)    ! 5
  call this%register_nvar_0d(this%h2,    'H2        ', varno%rh2m)
  call this%register_nvar_0d(this%t,     'T         ', varno%t)
  call this%register_nvar_0d(this%z,     'Z         ', varno%z)
  call this%register_nvar_0d(this%dz,    'DZ        ', varno%dz)
  call this%register_nvar_0d(this%lh,    'LH        ', varno%rhlay) ! 10
  call this%register_nvar_0d(this%t2,    'T2        ', varno%t2m)
  call this%register_nvar_0d(this%ts,    'TS        ', varno%ts)
  call this%register_nvar_0d(this%rad,   'RAD       ', varno%rawbt) 
  call this%register_nvar_0d(this%sn,    'SN        ', varno%sfall)
  call this%register_nvar_0d(this%rr,    'RR        ', varno%rr)    ! 15
  call this%register_nvar_0d(this%ps,    'PS        ', varno%ps)
  call this%register_nvar_0d(this%cc,    'CC        ', varno%satcl)
  !NVAR_CLW=18 was not used 
  call this%register_nvar_0d(this%q,     'Q         ', varno%q)
  call this%register_nvar_0d(this%ffs,   'FFS       ', varno%scatws) ! 20 in old money (see CLW)
  !NVAR_S0 =21 was not used
  !NVAR_X  =22 was not used
  call this%register_nvar_0d(this%pwc,   'PWC       ', varno%pwc) 
  call this%register_nvar_0d(this%to3,   'TO3       ', varno%o3lay) 
  call this%register_nvar_0d(this%tcw,   'TCW       ', varno%cllqw)  ! 25
  call this%register_nvar_0d(this%rfl,   'RFL       ', varno%refl) 
  call this%register_nvar_0d(this%apd,   'APD       ', varno%apdss) 
  call this%register_nvar_0d(this%ro,    'RO        ', varno%bend_angle)  
  call this%register_nvar_0d(this%hls,   'HLS       ', varno%los) 
  call this%register_nvar_0d(this%aod,   'AOD       ', varno%aerod)  ! 30
  call this%register_nvar_0d(this%lra,   'LRA       ', varno%limb_radiance) 

  call this%register_nvar_1d(this%chem, &
    & (/'NO2       ', 'SO2       ', 'CO        ', 'CH2O      ', 'O3        ', 'SO2VLC    '/), &
    & (/varno%chem1,varno%chem2,varno%chem3,varno%chem4,varno%chem5,varno%chem6/), &
    & kpchem_assim)
  !call this%register_nvar_1d(this%ghg, &
  !  & (/'CO2       ', 'CH4       ', 'N2O       '/), &
  !  & (/varno%ghg1,varno%ghg2,varno%ghg3/), &  ! AJGDB GHG3 added for first time - will it affect CAMS?
  !  & kpghg)
  ! jpghg_assim=2, jpghg=3 (previously nvar used jpghg=3 but ignored N2O in map_nvar...
  call this%register_nvar_1d(this%ghg, &
    & (/'CO2       ', 'CH4       '/), &
    & (/varno%ghg1,varno%ghg2/), &  
    & kpghg_assim)

  call this%register_nvar_0d(this%cod,   'COD         ', varno%cod)   ! 40
  call this%register_nvar_0d(this%rao,   'RAO         ', varno%rao) 
  call this%register_nvar_0d(this%rfa,   'RFA         ', varno%rfltnc) 
  call this%register_nvar_0d(this%oda,   'ODA         ', varno%od) 
  call this%register_nvar_0d(this%dow,   'DOW         ', varno%dopp) 
  call this%register_nvar_0d(this%sm,    'SM          ', varno%soilm) ! 45   

  ! AJGDB this next code preserves bug or intended feature inherited from old suvnmb.F90 (pre 2015?)
  ! This should probably be mapped to VARNO%NSOILM instead.
  call this%register_nvar_0d(this%smn,   'SMN         ', varno%flgt_phase) 

  call this%register_nvar_0d(this%sfr,   'SFR         ', varno%bt_real) 
  call this%register_nvar_0d(this%lrr,   'LRR         ', varno%lnprc) 
  call this%register_nvar_0d(this%lbs,   'LBS         ', varno%libksc) 
  call this%register_nvar_0d(this%q2,    'Q2          ', varno%q2m)      ! First previous 50
  call this%register_nvar_0d(this%td2,   'TD2         ', varno%td2m)     ! Second previous 50 (here things went wrong)
  call this%register_nvar_0d(this%n,     'N-CLOUD     ', varno%n) 
  call this%register_nvar_0d(this%pmsl,  'PMSL        ', varno%pmsl) 
  call this%register_nvar_0d(this%ptend, 'PTEND       ', varno%ptend) 
  call this%register_nvar_0d(this%du,    'U-SHEAR     ', varno%du)
  call this%register_nvar_0d(this%dv,    'V-SHEAR     ', varno%dv)
  call this%register_nvar_0d(this%vt,    'VIRTUAL-T   ', varno%vt)
  call this%register_nvar_0d(this%td,    'DEW-POINT   ', varno%td)
  call this%register_nvar_0d(this%w,     'PAST-WX     ', varno%w)
  call this%register_nvar_0d(this%ww,    'PRES-WX     ', varno%ww)
  call this%register_nvar_0d(this%vsp,   'VSP         ', varno%vsp)
  call this%register_nvar_0d(this%tsts,  'TSTS        ', varno%tsts)
  call this%register_nvar_0d(this%cbsc,  'CBSC        ', varno%lidar_cloud_backscatter)
  call this%register_nvar_0d(this%crfl,  'CRFL        ', varno%cloud_radar_reflectivity)
  call this%register_nvar_0d(this%lab,   'LAB         ', varno%liatbk)

  ! Print out current mappings (needed for setting NOTVAR correctly, for example)
  write(nulout,*) ' -------------------------------------'
  write(nulout,*) ' Mappings between NVAR and VARNO codes'
  write(nulout,*) ' -------------------------------------'
  write(nulout,*) '   NVAR NAME          VARNO'
  do i=1,this%nvarno_mappings
    write(nulout,'(I8,X,A10,X,I8)') this%nvar(i), this%cvar_name(this%nvar(i)), this%varnum(i) 
  enddo
  write(nulout,*) ' -------------------------------------'

  if (lhook) call dr_hook('class_nvar%create',1,zhook_handle)
end subroutine create

!------------------------------------------------------------------------
!
! Create internal mappings between VARNO and NVAR
!
subroutine register_nvar_0d(this,knvar,cdlv,kvarno,ld_permit_duplicate)
  class(class_nvar),  intent(inout)        :: this
  integer(kind=jpim), intent(inout)        :: knvar
  character(len=10),  intent(in)           :: cdlv
  integer(kind=jpim), intent(in)           :: kvarno
  logical(kind=jplm), intent(in), optional :: ld_permit_duplicate

  logical(kind=jplm) :: lfail_duplicate

  real(kind=jphook) :: zhook_handle
  if (lhook) call dr_hook('class_nvar%register_nvar_0d',0,zhook_handle)

  if ((this%jpxvar          >= max_varnos) .or. &
    & (this%nvarno_mappings >= max_varnos)) then
    call abor1('Increase MAX_VARNOS in nvar_class.F90')
  endif

  if (knvar == -1) then 

    ! New NVAR number assigned
    this%jpxvar = this%jpxvar + 1
    knvar       = this%jpxvar

    this%cvar_name(this%jpxvar) = cdlv

  else  

    lfail_duplicate = .true.
    if (present(ld_permit_duplicate)) then
      if (ld_permit_duplicate) then
        lfail_duplicate = .false.
      endif
    endif

    if (lfail_duplicate) then
      write(nulerr, *) 'CLASS_NVAR: duplicate NVAR mapping attempted from varno : ',kvarno
      call abor1('Duplicate varno mapping in CLASS_NVAR')
    endif

  endif

  if (this%nvarno_mappings >= 1) then
    if (any(this%varnum(1:this%nvarno_mappings) == kvarno)) then
      write(nulerr, *) 'CLASS_NVAR: attempt to re-register a VARNO : ',kvarno
      call abor1('Atempt to re-register a varno in CLASS_NVAR')
    endif
  endif

  ! New VARNO registered
  this%nvarno_mappings = this%nvarno_mappings + 1
  this%varnum(this%nvarno_mappings) = kvarno
  this%nvar  (this%nvarno_mappings) = this%jpxvar
  
  if (lhook) call dr_hook('class_nvar%register_nvar_0d',1,zhook_handle)
end subroutine register_nvar_0d

!------------------------------------------------------------------------
!
! Create internal mappings between VARNO and NVAR - CAMS variable arrays
!
subroutine register_nvar_1d(this,knvar,cdlv,kvarno,knvarno)
  class(class_nvar),           intent(inout) :: this
  integer(kind=jpim), allocatable, intent(inout) :: knvar(:)
  character(len=10),           intent(in)    :: cdlv(knvarno)
  integer(kind=jpim),          intent(in)    :: kvarno(knvarno)
  integer(kind=jpim),          intent(in)    :: knvarno  

  integer(kind=jpim) :: i

  real(kind=jphook) :: zhook_handle
  if (lhook) call dr_hook('class_nvar%register_nvar_1d',0,zhook_handle)

  ! To make sure that CAMS is using consistent CHEM/GHG array lengths 
  !if ((size(cdlv) /= knvarno) .or. (size(kvarno) /= knvarno)) then
  !  call abor1('CLASS_NVAR: inconsistent number of nvar and varno mappings.')
  !endif

  if (.not.allocated(knvar)) then
    allocate(knvar(knvarno))   
    knvar=-1
  else
    call abor1('CLASS_NVAR: re-initialisation is not allowed')
  endif

  do i=1,knvarno
    call this%register_nvar_0d(knvar(i), cdlv(i), kvarno(i))
  enddo

  if (lhook) call dr_hook('class_nvar%register_nvar_1d',1,zhook_handle)
end subroutine register_nvar_1d

!------------------------------------------------------------------------
!
! Given varno KVARNO, return nvar number KNVAR. KSTATUS is 0 if the
! mapping is successful, and 1 if not. Optionally return the character
! code associated with the nvar in CDLV
!
subroutine map_varno(this,kvarno,knvar,kstatus,cdlv)
  class(class_nvar),  intent(in)    :: this
  integer(kind=jpim), intent(in)    :: kvarno
  integer(kind=jpim), intent(out)   :: knvar, kstatus
  character(len=10) , intent(out), optional :: cdlv

  logical(kind=jplm) :: llfound
  integer(kind=jpim) :: i

  real(kind=jphook) :: zhook_handle
  if (lhook) call dr_hook('class_nvar%map_varno',0,zhook_handle)

  llfound=.false.
  do i=1,this%nvarno_mappings
    if(kvarno == this%varnum(i)) then
      llfound=.true.
      exit
    endif
  enddo

  if(llfound) then
    knvar   = this%nvar(i)
    kstatus = 0_jpim
    if (present(cdlv)) cdlv = this%cvar_name(knvar)
  else
    knvar   = nmdi
    kstatus = 1_jpim
    if (present(cdlv)) cdlv = 'xxxxxxxxxx'
  endif

  if (lhook) call dr_hook('class_nvar%map_varno',1,zhook_handle)
end subroutine map_varno

!------------------------------------------------------------------------
!
! Given an NVAR, return the corresponding character code
!
subroutine get_name(this,knvar,cdlv)
  class(class_nvar),  intent(in)    :: this
  integer(kind=jpim), intent(in)    :: knvar
  character(len=10) , intent(out)   :: cdlv

  real(kind=jphook) :: zhook_handle
  if (lhook) call dr_hook('class_nvar%get_name',0,zhook_handle)

  if ((knvar < 1) .or. (knvar > this%jpxvar)) call abor1('class_nvar%get_name: broken nvar indexing')
  cdlv=this%cvar_name(knvar)
  
  if (lhook) call dr_hook('class_nvar%get_name',1,zhook_handle)
end subroutine get_name

!------------------------------------------------------------------------
!
! Destructor
!
subroutine finish(this)
  class(class_nvar),  intent(inout) :: this

  real(kind=jphook) :: zhook_handle
  if (lhook) call dr_hook('class_nvar%finish',0,zhook_handle)

  if (allocated(this%chem)) deallocate(this%chem)
  if (allocated(this%ghg))  deallocate(this%ghg)
  
  if (lhook) call dr_hook('class_nvar%finish',1,zhook_handle)
end subroutine finish

!------------------------------------------------------------------------
!
! Automatic destructor
!
subroutine auto_finish(this)
  type(class_nvar),  intent(inout) :: this

  real(kind=jphook) :: zhook_handle
  if (lhook) call dr_hook('class_nvar%auto_finish',0,zhook_handle)

  call finish(this)
  
  if (lhook) call dr_hook('class_nvar%auto_finish',1,zhook_handle)
end subroutine auto_finish

end module nvar_class
 
