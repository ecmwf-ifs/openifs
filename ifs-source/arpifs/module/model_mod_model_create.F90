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

SUBMODULE (MODEL_MOD) MODEL_MOD_MODEL_CREATE
IMPLICIT NONE

CONTAINS

MODULE SUBROUTINE MODEL_CREATE(SELF, YDGEOMETRY, PTSTEP, CDNAMELIST, LDLINEAR, KGFLCONF,K_SELF)
USE PARKIND1,     ONLY : JPRB, JPIM
USE YOMHOOK,      ONLY : LHOOK, DR_HOOK, JPHOOK
USE GEOMETRY_MOD, ONLY : GEOMETRY
USE YOMLUN,       ONLY : NULNAM, NULOUT
USE YOE_CUCONVCA, ONLY : INI_CUCONVCA
USE YOMTRAJ,      ONLY : LTRAJALLOC
USE SPNG_MOD,     ONLY : SUSPNG
USE YOMCT0,       ONLY : NCONF, CNMEXP  ! TOWIL NAMARG HACK
USE YOMARG,       ONLY : LELAM, LECMWF  ! TOWIL NAMARG HACK
USE GFL_SUBS_MOD, ONLY : DEACT_CLOUD_GFL
USE YOMVAR,       ONLY : LMODERR
USE YOMFPC,       ONLY : LOCEDELAY
USE YEMLBC_MODEL    , ONLY : SUELBC_INIT

IMPLICIT NONE

TYPE(MODEL),            INTENT(INOUT) :: SELF
TYPE(GEOMETRY), TARGET, INTENT(INOUT) :: YDGEOMETRY
REAL(KIND=JPRB),        INTENT(IN)    :: PTSTEP
CHARACTER(LEN=*),       INTENT(IN)    :: CDNAMELIST
LOGICAL,                INTENT(IN)    :: LDLINEAR
INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: KGFLCONF
INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: K_SELF

INTEGER(KIND=JPIM) :: IOS,IMINUT
CHARACTER (LEN=35) :: CLINE = '----------------------------------'
LOGICAL :: LLOPENED

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"
#include "allocate_empty_trajectory.intfb.h"
#include "chem_init.intfb.h"




! #include "slcset.intfb.h"
#include "su0phy.intfb.h"
#include "suallo.intfb.h"
!!#include "sucfu.intfb.h"
#include "sudimf1.intfb.h"
#include "sudimf2.intfb.h"
#include "sudyna.intfb.h"
#include "sudyn.intfb.h"
#include "sugfl.intfb.h"
#include "sugrib.intfb.h"
#include "suinterpolator.intfb.h"
#include "sumcc.intfb.h"
#include "sumcclag.intfb.h"
#include "suphy.intfb.h"
#include "surand1.intfb.h"
#include "surip.intfb.h"
#include "suoph.intfb.h"
#include "sufa.intfb.h"
#include "susc2b.intfb.h"
#include "suspvariables.intfb.h"
#include "suvareps.intfb.h"
! #include "suxfu.intfb.h"
#include "suspsdt.intfb.h"
#include "get_spp_conf.intfb.h"
#include "ini_spp.intfb.h"
#include "suclopn.intfb.h"
#include "suswn.intfb.h"
#include "suaersn.intfb.h"
#include "sumoderrmod.intfb.h"
#include "sumddh.intfb.h"
#include "sunddh.intfb.h"
#include "sualdyn_ddh.intfb.h"
#include "sualmdh.intfb.h"
#include "sualtdh.intfb.h"
!#include "updecaec.intfb.h"

! TOWIL NAMARG HACK
!!!!!#include "namarg.nam.h"

IF (LHOOK) CALL DR_HOOK('MODEL_MOD:MODEL_CREATE',0,ZHOOK_HANDLE)

SELF%YRML_GCONF%GEOM => YDGEOMETRY
CALL MODEL_SET(SELF)
IF(PRESENT(K_SELF)) THEN
  SELF%MOBJECT_ID = K_SELF
  WRITE(CLINE,'(A13,I4.4,A18)') '  MOBJECT_ID=', SELF%MOBJECT_ID,' -----------------'
ENDIF
WRITE(SELF%COBJECT_ID,'(A,I4.4,A)') ' MOBJECT_ID=',SELF%MOBJECT_ID,' '

SELF%LINEAR_MODEL = LDLINEAR
WRITE(NULOUT,'(A,L2,A,I4.4,A)') '---- MODEL_CREATE ----------- LINEAR MODEL = ',SELF%LINEAR_MODEL, &
 & ' RESOL=',YDGEOMETRY%YRDIM%NSMAX,CLINE
!* Open namelist
INQUIRE(NULNAM,OPENED=LLOPENED)
IF (LLOPENED) CLOSE(NULNAM)
OPEN(NULNAM,FILE=CDNAMELIST,ACTION='READ',IOSTAT=IOS)
IF (IOS /= 0)&
  & CALL ABOR1("MODEL_MOD:MODEL_SETUP failed to open namelist "//TRIM(CDNAMELIST))
!*    Initialize dynamics: part A
WRITE(NULOUT,*) '--- Set up dynamics part A ---------',CLINE
CALL SUDYNA(YDGEOMETRY%YRDIM,SELF%YRML_DYN%YRDYNA,YDGEOMETRY%YRCVER%LVERTFE, &
 & YDGEOMETRY%YRCVER%NDLNPR,YDGEOMETRY%LNONHYD_GEOM,NULOUT)
!*    Initialize special keys for the climate version
WRITE(NULOUT,*) '---- Set up MCC climate model keys --',CLINE
CALL SUMCC(SELF%YRML_AOC%YRMCC,SELF%YRML_PHY_EC%YREPHY,NULOUT)

!*    Initialize YOMRIP variables
WRITE(NULOUT,*) '------ Set up YOMRIP variables ',CLINE
CALL SURIP(YDGEOMETRY%YRDIM,SELF%YRML_DYN%YRDYNA, &
 & SELF%YRML_GCONF%YRRIP,PTSTEP=PTSTEP)

!*    Initialize control of physical parameterizations
WRITE(NULOUT,*) '-- Set up physical parameterizations ',CLINE
CALL SU0PHY(SELF,NULOUT)

!! ky: maybe that too (formerly in SUDYN)?
!! !*    Initialize some dimensions for trajectory and background.
!! WRITE(NULOUT,*) '------ Set up some dimensions for trajectory and background ------',CLINE
!! CALL SUDIM_TRAJ

!*    Initialize file handling
WRITE(NULOUT,*) '---- Set up files handling, FA --',CLINE
CALL SUOPH(YDGEOMETRY)

!*    Initialize Arpege field names and GRIB packing options
WRITE(NULOUT,*) '--- Set up GRIB packing options ---',CLINE
CALL SUFA

!*    Initialize number of fields dimensions (part 2)
WRITE(NULOUT,*) '------ Set up number of fields dimensions, part 1  ------',CLINE
CALL SUDIMF1(SELF)

!*    Set up unified_treatment grid-point fields, from su0yoma:284
WRITE(NULOUT,*) '---- Set up unified_treatment grid-point fields -',CLINE
CALL SUGFL(YDGEOMETRY%YRDIMV,SELF,KGFLCONF)

!*    Initialize number of fields dimensions (part 2)
WRITE(NULOUT,*) '------ Set up number of fields dimensions, part 2  ------',CLINE
CALL SUDIMF2(SELF%YRML_GCONF,SELF%YRML_DYN%YRDYNA)

CALL SUSPVARIABLES(SELF%YRML_GCONF,SELF%YRML_DYN%YRDYNA%LNHX,&
 & SELF%YRML_DYN%YRDYNA%LNHDYN)

!*    Allocate grid point and spectral arrays, from su0yoma:347
WRITE(NULOUT,*) '------ Set up : array allocations ------',CLINE
CALL SUALLO(YDGEOMETRY,SELF)

!*    Initialize extended control variable options
IF ((NCONF == 1).OR.(NCONF == 131)) THEN
  WRITE(NULOUT,*) '--- Set up ECV geometry ---------',CLINE
  CALL SUECV(YDGEOMETRY%YRGEM,YDGEOMETRY%YRDIM,YDGEOMETRY%YRDIMV,SELF%YRML_GCONF%YRRIP,SELF%YRML_GCONF%YRDIMECV)
ENDIF

!*    Initialize DDH (Horizontal domains diagnostics)
!!$WRITE(NULOUT,*) '------ Turn off DDH diagnostics ----	----',CLINE
!!$CALL DDHOFF(YDGEOMETRY,SELF%YRML_GCONF,SELF%YRML_DIAG)

!! ky: maybe that too (formerly in SUDYN)?
!! !*    Setup TESTVAR.
!! WRITE(NULOUT,*) '---- Set up TESTVAR ----------',CLINE
!! CALL SETUP_TESTVAR

!*    Initialize Dynamics, from su0yomb:569 (necessary for susc2b)
WRITE(NULOUT,*) '---- Set up model dynamics ----------',CLINE
CALL SUDYN(YDGEOMETRY,SELF,NULOUT)

!*    Initialize vertical interpolator
WRITE(NULOUT,*) '---- Set up vertical interpolator -------',CLINE
CALL SUINTERPOLATOR(YDGEOMETRY,SELF%YRML_DYN%YRDYNA,SELF%YRML_DYN%YRSLINT)

!*    Initialize new sponge
WRITE(NULOUT,*) '---- Set up new sponge ----------',CLINE
CALL SUSPNG(SELF%YRML_DYN%YRSPNG,SELF%YRML_GCONF%YRRIP,SELF%YRML_DYN%YRDYNA, &
 & YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRSTA%STZ)

!*    Initialize Physics
WRITE(NULOUT,*) '---- Set up model physics -----------',CLINE
CALL SUPHY(YDGEOMETRY,SELF,NULOUT)
!*    Initialize special keys for the climate version 2nd part

WRITE(NULOUT,*) '---- Set up MCC climate model keys (lagged part) --',CLINE
CALL SUMCCLAG(YDGEOMETRY%YRGEM,SELF%YRML_GCONF,SELF%YRML_AOC,SELF%YRML_CHEM%YRCOMPO, &
     &  SELF%YRML_CHEM%YRCHEM, SELF%YRML_PHY_AER%YREAERSRC, SELF%YRML_PHY_EC%YREPHY, NULOUT)

!*    Initialize cumulated fluxes requests
!!WRITE(NULOUT,*) '------ Set up cumulated fluxes diags ---',CLINE
!!CALL SUCFU(YDGEOMETRY,SELF%YRML_DIAG%YRCFU,SELF%YRML_GCONF%YRRIP,SELF%YRML_PHY_RAD%YRERAD,SELF%YRML_PHY_MF%YRPHY,NULOUT)

!*    Initialize instantaneous fluxes requests
!!WRITE(NULOUT,*) '------ Set up instantaneous fluxes diags ',CLINE
!!CALL SUXFU(YDGEOMETRY,SELF%YRML_DIAG%YRXFU,SELF%YRML_GCONF%YRRIP,SELF%YRML_PHY_MF%YRPHY,NULOUT)

CALL SUNDDH(YDGEOMETRY,SELF,LDONLY_SWITCHES=.TRUE.)
!*    Initialize buffers for gridpoint scanning, part2, from su0yomb:624
WRITE(NULOUT,*) '---- Set up gridpoint scanning, part B ----',CLINE
IF (NCONF /= 901) CALL SUSC2B(YDGEOMETRY,SELF)

!*    Initialize forcing by coarser model: part A
WRITE(NULOUT,*) '--- Set up forcing by coarser model part A ---------',CLINE
CALL SUELBC_INIT(SELF%YRML_DYN%YRDYNA,SELF%YRML_LBC)

!!** OLIVIER addition, to fill in NGRIB_HANDLE_XX variables
CALL SUVAREPS(SELF%YRML_GCONF%YRRIP)
CALL SUGRIB(YDGEOMETRY%YRDIM,SELF%YRML_PHY_EC%YREPHY,SELF%YRML_PHY_G%YRDPHY,SELF%YRML_PHY_MF%YRPHY)

!    If required allocate empty skeleton of trajectory structure
IF (.NOT. LTRAJALLOC .AND. (NCONF /= 401) .AND. (NCONF /= 501) .AND.&
                         & (NCONF /= 601) .AND. (NCONF /= 801)) THEN
  CALL ALLOCATE_EMPTY_TRAJECTORY(YDGEOMETRY%YRDIM,SELF%YRML_GCONF%YRRIP)
ENDIF

!*    Set up Stochastic Physics
WRITE(NULOUT,*) '--- Set up stochastic physics, SPBS, CABS ',CLINE
CALL INI_CUCONVCA(YDGEOMETRY,SELF%YRML_DYN%YRDYNA,&
 & SELF%YRML_PHY_EC%YRECUCONVCA, SELF%YRML_DYN%YRSL)
CALL SURAND1(YDGEOMETRY,SELF%YRML_PHY_STOCH,SELF%YRML_DYN%YRDYN,SELF%YRML_GCONF%YRRIP,SELF%YRML_PHY_EC%YRECUCONVCA)

!     Set up spectral stochastic diabatic tendencies
WRITE(NULOUT,'(A72)') '--- Set up stochastically perturbed parametrization tendencies '//CLINE
WRITE(NULOUT,*)       '      SPPT a.k.a. stochastic physics with spectral pattern'
CALL SUSPSDT(YDGEOMETRY, SELF%YRML_GCONF%YRRIP, SELF%YRML_GCONF%YRSPPT_CONFIG, SELF%YRML_SPPT)

!     Set up stochastically perturbed parameterisation scheme
WRITE(NULOUT,*) '--- Set up stochastically perturbed parametrization scheme (SPP)  ',CLINE
CALL GET_SPP_CONF(SELF%YRML_GCONF%YRRIP, SELF%YRML_GCONF%YRSPP_CONFIG)
CALL INI_SPP(YDGEOMETRY, SELF%YRML_GCONF%YRRIP, SELF%YRML_GCONF%YRSPP_CONFIG, SELF%YRML_SPP)

!antje Chemistry not called in minimization

IF (SELF%YRML_GCONF%YGFL%NCHEM > 0 ) THEN
  IF(LDLINEAR) THEN
    IF(SELF%YRML_CHEM%YRCHEM%LCHEM_TL) THEN
      WRITE(NULOUT,*) '--- Set up linear chemistry  ',CLINE
      CALL CHEM_INIT_TLAD(YDGEOMETRY,SELF%YRML_GCONF,SELF%YRML_CHEM)
    ELSE
      WRITE(NULOUT,*) '--- No linear chemistry  ',CLINE
    ENDIF
  ELSE
    WRITE(NULOUT,*) '--- Set up chemistry  ',CLINE
    CALL CHEM_INIT(YDGEOMETRY,SELF%YRML_GCONF,SELF%YRML_DYN%YRDYNA,SELF%YRML_CHEM)
  ENDIF
ENDIF

IF (LMODERR) THEN
  !*    Initialize model error control
  WRITE(NULOUT,*) '------ Set up YOMODERRMOD variables ',CLINE
  CALL SUMODERRMOD(SELF%YRML_GCONF%YRMODERR,SELF%YRML_GCONF%YRRIP,CDNAMELIST=CDNAMELIST)
ENDIF

!*    Initialize DDH (Horizontal domains diagnostics)
WRITE(NULOUT,*) '------ Set up DDH diagnostics --------',CLINE
CALL SUNDDH(YDGEOMETRY,SELF)
CALL SUALMDH(YDGEOMETRY%YRGEM,SELF%YRML_DIAG)
IF(SELF%YRML_DIAG%YRLDDH%LSDDH) THEN
  WRITE(NULOUT,*) '---- Set up DDH diagnostic domains ',CLINE
  CALL SUMDDH(YDGEOMETRY,SELF%YRML_DIAG)
ENDIF
!*    Memory allocation for cumulated DDH arrays (horizontal domains diags)
WRITE(NULOUT,*) '---- Set up DDH diags allocation --',CLINE
CALL SUALTDH(YDGEOMETRY%YRDIMV,SELF%YRML_DIAG,SELF%YRML_PHY_MF%YRARPHY,SELF%YRML_PHY_MF%YRPHY)

!*    Memory allocation for dynamical DDH tendencies arrays
WRITE(NULOUT,*) '---- Set up dynamical DDH arrays allocation --',CLINE
CALL SUALDYN_DDH(YDGEOMETRY,SELF%YRML_DIAG,SELF%YRML_GCONF,&
 & SELF%YRML_DYN%YRDYNA%LNHDYN,SELF%YRML_DYN%YRDYNA%LNHX)


CLOSE(NULNAM)
! This should be a temporary fudge
OPEN(NULNAM,FILE='fort.4',ACTION='READ',IOSTAT=IOS)

!Modify some settings if it is an OOPS linear model being set up
IF(SELF%LINEAR_MODEL) THEN
  WRITE(NULOUT,*) '---- Mode for linear model -----------',CLINE
  IF(SELF%YRML_PHY_SLIN%YRPHNC%LENCLD2.OR.SELF%YRML_PHY_SLIN%YRPHNC%LEPCLD2) &
   &  CALL DEACT_CLOUD_GFL(YDGEOMETRY%YRDIMV,YDGEOMETRY%YRCVER,SELF)
  SELF%YRML_PHY_SLIN%YREPHLI%LPHYLIN = .TRUE.
  WRITE(NULOUT,*) 'NSW = ',SELF%YRML_PHY_RAD%YRERAD%NSW,SELF%YRML_PHY_SLIN%YRPHNC%LERADSW2
  IF(SELF%YRML_PHY_SLIN%YRPHNC%LERADSW2 .AND. SELF%YRML_PHY_RAD%YRERAD%NSW > 2) THEN
    SELF%YRML_PHY_RAD%YRERAD%NSW=2
    CALL SUSWN(SELF%YRML_PHY_RAD%YRESWRT,SELF%YRML_PHY_RAD%YRERAD,SELF%YRML_PHY_RAD%YRERAD%NTSW,SELF%YRML_PHY_RAD%YRERAD%NSW)
    CALL SUCLOPN(SELF%YRML_PHY_RAD%YRERAD,SELF%YRML_PHY_RAD%YRESWRT,SELF%YRML_PHY_RAD%YRERAD%NTSW,&
     & SELF%YRML_PHY_RAD%YRERAD%NSW,YDGEOMETRY%YRDIMV%NFLEVG)
    CALL SUAERSN (SELF%YRML_PHY_RAD%YRESWRT,SELF%YRML_PHY_RAD%YRERAD%NTSW, SELF%YRML_PHY_RAD%YRERAD%NSW)
  ENDIF
ENDIF

!*    Set up for nemo instead of SUFPC
LOCEDELAY = LECMWF

!CLOSE(NULNAM)

IF (LHOOK) CALL DR_HOOK('MODEL_MOD:MODEL_CREATE',1,ZHOOK_HANDLE)

END SUBROUTINE MODEL_CREATE

END SUBMODULE
