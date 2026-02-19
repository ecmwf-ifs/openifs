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

MODULE BASCOE_TUV_MODULE

USE PARKIND1, ONLY : JPIM, JPRB
USE BASCOE_J_MODULE, only : jprocmax => ndiss

IMPLICIT NONE

SAVE

!**   DESCRIPTION 
!     ----------
!
!   Module with variables and parameters for routines implementing an adapted version of
!       TUV radiation transfer
!
!     MODIFICATIONS.
!     --------------
!     2017-12-07: YVES CHRISTOPHE (YC) *BIRA*
!       imported and adapted from BASCOE tuv_midatm.com.f90
!     2019-12-18:    YVES CHRISTOPHE    *BIRA* !YC
!         allow date-dependent solar flux
!----------------

! from module TUV_MIDATM_GRID
  INTEGER(KIND=JPIM), parameter :: mxwvn = 171      ! Number of WaveLength intervals
  INTEGER(KIND=JPIM), parameter :: iv_lya = 8       ! Index of Lyman-a wl. wlmid(iv_lya)=121.6nm
  INTEGER(KIND=JPIM), parameter :: iv_srb0 = 46     ! Indexes of wavelengths where Schuman-Runge bands [SRB] begin and end:
  INTEGER(KIND=JPIM), parameter :: iv_srb1 = 61     !   wlmid(iv_srb0)=173nm, wlmid(iv_srb1)=204nm
  INTEGER(KIND=JPIM), parameter :: nwvn_srb = 16    !   nb of wavelength intervals in SRB, for parameterization (Kockarts, 1994)

  REAL(KIND=JPRB), dimension(mxwvn)   :: wlmid      ! Wavelengths at middle of wl intervals (nm)     
  REAL(KIND=JPRB), dimension(mxwvn+1) :: wlbnd      ! Wavelengths at boundaries of wl intervals (nm) 
                                                    !   *WARNING* this construct for wavelength boundaries implicitly
                                                    !   assume that the intervals are adjacent, which *IS NOT* the case
                                                    !   for the our grid with a discontinuity at lyman-alpha !: *USE WITH CAUTION*

! from  module TUV_MIDATM_PARMS
  LOGICAL, parameter :: sw_lya = .true.
  INTEGER(KIND=JPIM), parameter :: nabspec = 6       ! number of species absorbing radiation
  INTEGER(KIND=JPIM), parameter :: nheatspec = 2     ! nb of specs releasing heat
  REAL(KIND=JPRB), parameter ::  sza_first_light = 96. ! all J = 0 for sza > sza_first_light
  REAL(KIND=JPRB), parameter ::  sza_chap = 65.      ! threshold to call Chapman function

! from  module TUV_MIDATM_VARS
  LOGICAL :: init_done = .false.
  INTEGER(KIND=JPIM), dimension(jprocmax) ::  ivbegin, ivend

  REAL(KIND=JPRB), dimension(jprocmax,mxwvn) :: crs   ! T-indep cross-sections (cm2)
  REAL(KIND=JPRB), dimension(jprocmax,mxwvn) :: qy    ! quantum yields

  REAL(KIND=JPRB), dimension(mxwvn) :: fbeamr
  LOGICAL :: daily_solflux = .false.
  INTEGER(KIND=JPIM), dimension(:), allocatable    :: fbeam_dates      ! dates YYYYmmdd in solflux file
  REAL(KIND=JPRB), dimension(:,:), allocatable  :: fbeamr2d         ! date-dependent solflux data

  REAL(KIND=JPRB), dimension(mxwvn) :: raycrs
  REAL(KIND=JPRB), dimension(mxwvn) :: tb_o3, tc_o3, ta_no2, tb_hno3,  &! specific parameters for
       &                               ta1_clono2, ta2_clono2,         &!   computation of
       &                               ta_ch2o, tb_ch2o, tb_pan,       &!   temperature dependent
       &                               crs200co2, crs370co2,           &!   cross-sections and quantum yields
       &                               gamma_ch2o,                     &!   see CRS_TDEP
       &                               ta_n2o5, tb_n2o5,               &!   ...
       &                               crs295ocs
  REAL(KIND=JPRB), dimension(3,mxwvn) :: qy_jno3_o, qy_jno3_o2
  REAL(KIND=JPRB), dimension(12,nwvn_srb)      :: ako, bko       ! params to calc rm, ro2

! from  module TUV_MIDATM_SOLHEAT
  REAL(KIND=JPRB), dimension(mxwvn) :: hv         ! hPl*c/lambda (J)
  REAL(KIND=JPRB), dimension(mxwvn) :: denerg_O2
  REAL(KIND=JPRB), dimension(mxwvn) :: denerg_O3

! from  module TUV_MIDATM_QYO3
  REAL(KIND=JPRB), dimension(mxwvn) :: coeff1, coeff2, coeff3

! from  module TUV_MIDATM_NAMES
  INTEGER(KIND=JPIM), parameter :: air=1, O2abs=2, O3abs=3, &
       &                            NOabs=4, CO2abs=5, NO2abs=6
  REAL(KIND=JPRB), parameter, dimension(nabspec) :: abs_molmass = &
       &  (/  28.9644d0,      & ! average molar mass of dry air (g/mole),valid under ~ 100km
       &      32.d0,          & ! O2
       &      48.d0,          & ! O3
       &      30.d0,          & ! NO
       &      44.d0,          & ! CO2
       &      46.d0         /)  ! NO2
  CHARACTER(len=3), dimension(nabspec) :: abs_name = &
     &  (/ 'air', 'O2 ', 'O3 ', 'NO ', 'CO2', 'NO2' /)


END MODULE BASCOE_TUV_MODULE
