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

#ifdef VPP
!OCL SCALAR
#endif
SUBROUTINE SUDYNA(YDDIM,YDDYNA,LDVERTFE,KDLNPR,LDNONHYD_GEOM,KULOUT)

!**** *SUDYNA*   - Initialize constants and control for the dynamics: "IFS_INIT" part

!     Purpose.
!     --------
!           Initialize YOMDYNA, a "IFS_INIT" module that controls the dynamics of the model.

!**   Interface.
!     ----------
!        *CALL* *SUDYNA(KULOUT)

!        Explicit arguments :
!        --------------------
!        KULOUT : Logical unit for the output

!        Implicit arguments :
!        --------------------
!        COMMON YOMDYNA

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!     none

!     Called by SU0YOMA.

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        K. YESSAD (CNRM/GMAP)
!        Original : September 2002

! Modifications
! -------------
!   K. Yessad (Nov 2009): prune lpc_old.
!   K. Yessad (Dec 2009): LRWSDLW,LRWSDLR,LRWSDLG=T,T,T in NH model for LGWADV=F.
!   F. Vana  22-Feb-2011: option L3DTURB
!   K. Yessad (Nov 2011): introduce LRALTVDISP,LRPRSLTRJ,LCURVW,LRFRICISOTR
!   N. Wedi (Nov 2011): LGRADSP
!   K. Yessad (Sep 2013): final value of LSLHD, and LSLHDQUAD set-up there.
!   P. Smolikova and J. Vivoda (Oct 2013): new options for VFE-NH
!   S. Malardel (Nov 2013): setup COMAD switches and tests
!   K. Yessad (July 2014): Various setup and module refactoring.
!   F. Vana 13-Feb-2014  LSLHDVER
!   K. Yessad (Dec 2016): Prune obsolete options.
!   K. Yessad (June 2017): Introduce NHQE model.
!   F. Vana (21-Nov-2017): Option LHOISLT
!   K. Yessad (Feb 2018): remove deep-layer formulations.
!   K. Yessad (Apr 2018): introduce key L_RDRY_VD (ensure consistent definition of "dver" everywhere).
!   K. Yessad (June 2018): Alternate NHEE SI scheme elimination.
!   K. Yessad (July 2018): Add LNHEE_SVDLAPL_FIRST=T option.
!   J. Vivoda (July 2018): mixed NESC/SETTLS scheme.
!   F. Vana (20-Feb-2019): Option LRHSWENO & better defaults
!   F. Vana  18-Jul-2019: SLVF & removed LSLHDVER
!   I. Polichtchouk (June 2020): LPGFSF
!   F. Vana  21-Sep-2020: LPGREUSE
! End Modifications
!------------------------------------------------------------------

USE YOMDIM    , ONLY : TDIM
USE YOMDYNA   , ONLY : TDYNA
USE PARKIND1  , ONLY : JPIM, JPRB, JPRD
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMMP0    , ONLY : NPRINTLEV,LSCMEC
USE YOMCT0    , ONLY : LRPLANE, LELAM, LECMWF, LR2D, NCONF
USE YOMLUN    , ONLY : NULNAM
USE INTDYN_MOD, ONLY : SUPTR_TTND,SUPTR_TGMVT,SUPTR_TGFLT

!------------------------------------------------------------------

IMPLICIT NONE


TYPE(TDIM), INTENT(IN)            :: YDDIM
TYPE(TDYNA),TARGET, INTENT(INOUT) :: YDDYNA
LOGICAL, INTENT(IN)               :: LDVERTFE
INTEGER(KIND=JPIM), INTENT(IN)    :: KDLNPR
LOGICAL, INTENT(IN)               :: LDNONHYD_GEOM
INTEGER(KIND=JPIM), INTENT(IN)    :: KULOUT

!     ------------------------------------------------------------------
LOGICAL,POINTER :: LNHDYN,LNHEE,LNHQE
LOGICAL,POINTER :: LSPRT,LTWOTL,LSLAG,LSACC

INTEGER(KIND=JPIM),POINTER :: NPDVAR,NVDVAR,ND4SYS,NGWADVSI
LOGICAL,POINTER :: LNHX,LNHXDER,LGWADV,LRDBBC,LSI_NHEE,LNHQE_C2,LNHQE_SIHYD,LNHEE_SVDLAPL_FIRST

LOGICAL,POINTER :: LSLHD,LSLHD_W,LSLHD_T,LSLHD_SPD,LSLHD_SVD,LSLHD_GFL,LSLHD_OLD,LSLHD_STATIC,LSLHDQUAD
REAL(KIND=JPRB),POINTER :: SLHDKMIN,SLHDKMAX,SLHDKREF,SLHDEPSH,SLHDEPSV

LOGICAL,POINTER :: LSLVF
REAL(KIND=JPRB),POINTER :: SLVFEPS,SLVFMIN,SLVFTR,P_SLVF_BOT,P_SLVF_TR
LOGICAL,POINTER :: LGRADSP,LPGFSF,LPGREUSE,L3DTURB,LSLDIA,LRPRSLTRJ

LOGICAL,POINTER :: LRALTVDISP,LVSPLIP,LRHSVWENO,LCOMAD,LCOMADH,LCOMADV,LCOMAD_W,LCOMAD_T,&
 & LCOMAD_SPD,LCOMAD_SVD,LCOMAD_SP,LCOMAD_GFL,LNESCT,LNESCV,LNESC,LSETTLST,LSETTLSV,LSETTLS, &
 & LELTRA,LSLINLC2,LSLINL,LHOISLT,LSLTVWENO
REAL(KIND=JPRB),POINTER :: HOISLTV,HOISLTH

LOGICAL,POINTER :: LRUBC,LPC_FULL,LPC_CHEAP,LPC_CHEAP2,LMIXETTLS,LMIXETTLS_PRINT,L_RDRY_VD
REAL(KIND=JPRB),POINTER :: RMIXNL_TRH
LOGICAL :: LLGRID
INTEGER(KIND=JPIM) :: IERR
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "posnam.intfb.h"

#include "namdyna.nam.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUDYNA',0,ZHOOK_HANDLE)
! Associate pointers for variables in namelist
!     ------------------------------------------------------------------
LNHDYN => YDDYNA%LNHDYN
LNHEE => YDDYNA%LNHEE
LNHQE => YDDYNA%LNHQE
LSPRT => YDDYNA%LSPRT
LTWOTL => YDDYNA%LTWOTL
LSLAG  => YDDYNA%LSLAG
LSACC  => YDDYNA%LSACC
NPDVAR => YDDYNA%NPDVAR
NVDVAR => YDDYNA%NVDVAR
ND4SYS => YDDYNA%ND4SYS
LNHX => YDDYNA%LNHX
LNHXDER => YDDYNA%LNHXDER
LGWADV => YDDYNA%LGWADV
NGWADVSI => YDDYNA%NGWADVSI
LRDBBC => YDDYNA%LRDBBC
LSI_NHEE => YDDYNA%LSI_NHEE
LNHQE_C2 => YDDYNA%LNHQE_C2
LNHQE_SIHYD => YDDYNA%LNHQE_SIHYD
LNHEE_SVDLAPL_FIRST => YDDYNA%LNHEE_SVDLAPL_FIRST
LSLHD => YDDYNA%LSLHD
LSLHD_W => YDDYNA%LSLHD_W
LSLHD_T => YDDYNA%LSLHD_T
LSLHD_SPD => YDDYNA%LSLHD_SPD
LSLHD_SVD => YDDYNA%LSLHD_SVD
LSLHD_GFL => YDDYNA%LSLHD_GFL
LSLHD_OLD => YDDYNA%LSLHD_OLD
LSLHD_STATIC => YDDYNA%LSLHD_STATIC
LSLHDQUAD => YDDYNA%LSLHDQUAD
SLHDKMIN => YDDYNA%SLHDKMIN
SLHDKMAX => YDDYNA%SLHDKMAX
SLHDKREF => YDDYNA%SLHDKREF
SLHDEPSH => YDDYNA%SLHDEPSH
SLHDEPSV => YDDYNA%SLHDEPSV
LSLVF => YDDYNA%LSLVF
SLVFEPS => YDDYNA%SLVFEPS
SLVFMIN => YDDYNA%SLVFMIN
SLVFTR => YDDYNA%SLVFTR
P_SLVF_BOT => YDDYNA%P_SLVF_BOT
P_SLVF_TR => YDDYNA%P_SLVF_TR
LGRADSP => YDDYNA%LGRADSP
LPGFSF => YDDYNA%LPGFSF
LPGREUSE => YDDYNA%LPGREUSE
L3DTURB => YDDYNA%L3DTURB
LSLDIA => YDDYNA%LSLDIA
LRPRSLTRJ => YDDYNA%LRPRSLTRJ 
LRALTVDISP => YDDYNA%LRALTVDISP
LVSPLIP => YDDYNA%LVSPLIP
LRHSVWENO => YDDYNA%LRHSVWENO
LCOMAD => YDDYNA%LCOMAD
LCOMADH => YDDYNA%LCOMADH
LCOMADV => YDDYNA%LCOMADV
LCOMAD_W => YDDYNA%LCOMAD_W
LCOMAD_T => YDDYNA%LCOMAD_T
LCOMAD_SPD => YDDYNA%LCOMAD_SPD
LCOMAD_SVD => YDDYNA%LCOMAD_SVD
LCOMAD_SP => YDDYNA%LCOMAD_SP
LCOMAD_GFL => YDDYNA%LCOMAD_GFL
LNESCT => YDDYNA%LNESCT
LNESCV => YDDYNA%LNESCV
LNESC => YDDYNA%LNESC
LSETTLST => YDDYNA%LSETTLST
LSETTLSV => YDDYNA%LSETTLSV
LSETTLS => YDDYNA%LSETTLS
LELTRA => YDDYNA%LELTRA
LSLINLC2 => YDDYNA%LSLINLC2
LSLINL => YDDYNA%LSLINL
LHOISLT => YDDYNA%LHOISLT
HOISLTV => YDDYNA%HOISLTV
HOISLTH => YDDYNA%HOISLTH
LSLTVWENO => YDDYNA%LSLTVWENO
LRUBC => YDDYNA%LRUBC
LPC_FULL => YDDYNA%LPC_FULL
LPC_CHEAP => YDDYNA%LPC_CHEAP
LPC_CHEAP2 => YDDYNA%LPC_CHEAP2
LMIXETTLS => YDDYNA%LMIXETTLS
LMIXETTLS_PRINT => YDDYNA%LMIXETTLS_PRINT
RMIXNL_TRH => YDDYNA%RMIXNL_TRH
L_RDRY_VD => YDDYNA%L_RDRY_VD

!*       1.    Set default values.
!              -------------------
LNHEE=.FALSE.
LNHQE=.FALSE.
LNHDYN=(LNHEE.OR.LNHQE)

YDDYNA%LSACC=.FALSE.
LSLAG=.TRUE.
IF (LECMWF) THEN
  LSPRT =.TRUE.
  LTWOTL=.TRUE.
ELSE
  LSPRT =.FALSE.
  LTWOTL=.TRUE.
ENDIF

LVSPLIP=.FALSE. ! Provisional value, may be modified via NAMGFL.

LGWADV=.FALSE.
NGWADVSI=1
NPDVAR=2
IF (LDNONHYD_GEOM) THEN
  NVDVAR=4
  LRDBBC=LSLAG
ELSE
  NVDVAR=3
  LRDBBC=.FALSE.
ENDIF
LSLHD_SPD =.FALSE.
LSLHD_SVD =.FALSE.
LSLHD_GFL =.FALSE.
IF (LECMWF) THEN
  ! IFS specific defaults
  LSLHD_W   =.FALSE.
  LSLHD_T   =.FALSE.
  LSLHD_OLD =.FALSE.
  SLHDKMAX=6.0_JPRB
ELSE 
  LSLHD_W   =.FALSE.
  LSLHD_T   =.FALSE.
  LSLHD_OLD =.TRUE.
  SLHDKMAX=1.0_JPRB
ENDIF
LSLHD_STATIC=.FALSE.
SLHDKMIN=0.0_JPRB
SLHDKREF=-999.0_JPRB
SLHDEPSH=0.0_JPRB
SLHDEPSV=0.0_JPRB
! LSLHD defaults for GFLs are set in SUDEFO_GFLATTR, before reading NAMGFL!!!
LCOMAD    =.FALSE.
LCOMADH   =.FALSE.
LCOMADV   =.FALSE.
LCOMAD_W  =.FALSE.
LCOMAD_T  =.FALSE.
LCOMAD_SPD=.FALSE.
LCOMAD_SVD=.FALSE.
LCOMAD_SP =.FALSE.
LCOMAD_GFL =.FALSE.

LGRADSP=.FALSE.
LPGFSF=.FALSE.
LPGREUSE=.FALSE.
IF (LECMWF) THEN
  IF (LSLAG) THEN
    LGRADSP=.TRUE.
  ENDIF
ENDIF
ND4SYS=2

LSLDIA=.FALSE.

L3DTURB=.FALSE.


LRUBC=.FALSE.

LRALTVDISP=.FALSE.
IF (LECMWF) THEN
  LRPRSLTRJ=.TRUE.
ELSE
  LRPRSLTRJ=(NPRINTLEV > 0)
ENDIF

LPC_FULL=.FALSE.
LPC_CHEAP=.FALSE.
LPC_CHEAP2=.FALSE.

LNESCT=.FALSE.
LNESCV=.FALSE.
LNESC=.FALSE.

LELTRA=.FALSE.
LSETTLST=.TRUE.
LSETTLSV=.TRUE.
LSETTLS=.TRUE.

LMIXETTLS =.FALSE.
LMIXETTLS_PRINT =.FALSE.
RMIXNL_TRH = 0.5

L_RDRY_VD=.NOT.LNHQE

LNHQE_SIHYD=.FALSE.
LSI_NHEE=.FALSE.
LNHEE_SVDLAPL_FIRST=.FALSE.

LHOISLT=.false.
HOISLTV=0._JPRB
HOISLTH=0._JPRB
LSLTVWENO=.false.
IF (LECMWF) THEN
  LRHSVWENO=.true.
  LSLVF=.true.
ELSE
  LRHSVWENO=.false.
  LSLVF=.false.
ENDIF

SLVFEPS=0.005_JPRB
P_SLVF_BOT  = 0._JPRB ! SLVF applied from p=85hPa upwards
SLVFMIN=0.1_JPRB          
SLVFTR=0.0001_JPRB           
P_SLVF_TR= 300._JPRB    ! pressure where SLVF strength goes from SLVFMINx to 1x SLVFEPS
!     ------------------------------------------------------------------

!*       2.    Modify default values.
!              ----------------------

!        2.1   Read namelist.

CALL POSNAM(NULNAM,'NAMDYNA')
READ(NULNAM,NAMDYNA)

!*       2.2   Reset variables and test.

LNHDYN=(LNHEE.OR.LNHQE)

IF (NCONF == 701) LSPRT = .FALSE.
! * Keys SETTLS.. can be T only for SL2TL advection.
IF ( .NOT.LTWOTL ) THEN
  LSETTLST = .FALSE.
  LSETTLSV = .FALSE.
  LSETTLS  = .FALSE.
  WRITE(KULOUT,*) ' LTWOTL=F requires LSETTLST=LSETTLSV=LSETTLS=F, reset to F'
ENDIF
IF (.NOT.LSLAG) THEN
  IF (.NOT.LSCMEC) THEN
    LTWOTL=.FALSE.
  ENDIF
ENDIF
IF(LNHDYN .AND. .NOT. LDNONHYD_GEOM) THEN
  IERR=IERR+1
  WRITE(KULOUT,'(1X,A,I3,A)') ' ! SUDYNA: ERROR NR ',IERR,' !!!'
  WRITE (KULOUT,*) ' LNHDYN true REQUIRES LNONHYD_GEOM true'
ENDIF
! * Ask for TL or AD code of NH model?
IF (LNHDYN.AND. (ANY(SPREAD(NCONF,1,4) == (/401,501,601,801/)))) THEN  
  IERR=IERR+1
  WRITE(KULOUT,'(1X,A,I3,A)') ' ! SUDYNA: ERROR NR ',IERR,' !!!'
  WRITE (KULOUT,*) ' LNHDYN TL and AD currently not coded'
ENDIF

! * Save T9 arrays in I/O buffer is mandatory in NHS tl ad models
IF (LNHEE.AND.LNHQE) THEN
  IERR=IERR+1
  WRITE(KULOUT,'(1X,A,I3,A)') ' ! SUDYNA: ERROR NR ',IERR,' !!!'
  WRITE (KULOUT,*) ' LNHEE=T or LNHQE=T but both cannot be T at the same time'
ENDIF
! * Ask for SL scheme for configurations coded only for Eulerian scheme?
IF(LSLAG .AND. (NCONF == 202 .OR. NCONF == 422 .OR. NCONF == 522)) THEN
  IERR=IERR+1
  WRITE(KULOUT,'(1X,A,I3,A)') ' ! SUARG: ERROR NR ',IERR,' !!!'
  WRITE (KULOUT,*)' SEMI-LAGRANGIAN NOT CODED FOR THIS CONFIGURATION!'
ENDIF
! * NO LSLAG FOR CREATION OF CLIM FILES
IF (LSLAG .AND. NCONF == 923) LSLAG=.FALSE. 

! * NPDVAR, NVDVAR.
IERR=0
IF (LNHEE) THEN
  IF (.NOT.(NPDVAR == 2 .AND. (NVDVAR==3 .OR. NVDVAR==4))) THEN
    ! such configuration does not exist.
    IERR=IERR+1
    WRITE(KULOUT,'(1X,A,I3,A)') ' ! SUDYNA: ERROR NR ',IERR,' !!!'
    WRITE(KULOUT,*) ' Not existing configuration for (NPDVAR,NVDVAR)!'
  ENDIF
ENDIF
IF (LNHQE) THEN
  IF (.NOT.(NPDVAR == 2 .AND. NVDVAR==4)) THEN
    ! such configuration does not exist.
    IERR=IERR+1
    WRITE(KULOUT,'(1X,A,I3,A)') ' ! SUDYNA: ERROR NR ',IERR,' !!!'
    WRITE(KULOUT,*) ' Not existing configuration for (NPDVAR,NVDVAR)!'
  ENDIF
ENDIF
IF (.NOT.LNHDYN) THEN
  ! hydrostatic model.
  ! and in this case,
  !  * (LGWADV,NGWADVSI,LRDBBC) is always (F,1,F,1).
  !  * (NPDVAR,NVDVAR) is always (2,3) for safety.
  !  * (LSLHD_SPD,LSLHD_SVD) is always (F,F).
  LGWADV=.FALSE.
  NGWADVSI=1
  LRDBBC=.FALSE.
  NPDVAR=2
  NVDVAR=3
  LSLHD_SPD = .FALSE.
  LSLHD_SVD = .FALSE.
ENDIF

IF (LPGFSF.AND.(LNHQE.OR.LNHEE)) THEN
    IERR=IERR+1
    WRITE(KULOUT,'(1X,A,I3,A)') ' ! SUDYNA: ERROR NR ',IERR,' !!!'
    WRITE(KULOUT,*) ' LPGFSF not yet available for NH dynamics!'
ENDIF
! Wnen no precomputation of PG terms no point to allocate memory 
IF (.NOT.LPGFSF) LPGREUSE=.FALSE.

! L_RDRY_VD=T not coded in NHQE model (dver always defined with R_moist).
IF (LNHQE) L_RDRY_VD=.FALSE.

!*       2.3   Final set up of model internal keys LSLHD and COMAD

! LSLHD/LCOMAD SHOULD NOT be modified further, but additional checkings will be required in SUCTRL_GFLATTR.
LSLHD =LSLHD_W .OR. LSLHD_T .OR. LSLHD_SPD .OR. LSLHD_SVD .OR. LSLHD_GFL
LCOMAD=LCOMAD_W .OR. LCOMAD_T .OR. LCOMAD_SPD .OR. LCOMAD_SVD .OR. LCOMAD_SP .OR. LCOMAD_GFL

!*       2.4   Adjusting of SLHDKMIN and SLHDKREF, set-up LSLHDQUAD

IF ( LSLHD_OLD.AND.(SLHDKMIN /= 0.0_JPRB) ) THEN
  SLHDKMIN=0.0_JPRB
  WRITE(KULOUT,*) 'SUDYNA: Warning!'
  WRITE(KULOUT,*) '  SLHDKMIN reset to zero for LSLHD_OLD.'
ENDIF
IF(SLHDKREF == -999.0_JPRB) SLHDKREF = SLHDKMIN

! WARNING: Using non-zero SLHDKMIN redefines weights for high order
! semi-Lagrangian interpolators (subroutines LAIDDI, LAITRI) for
! all advected fields, even if they are not subject to SLHD scheme.
! Computation of modified interpolation weights uses same auxiliary
! quantities as the new SLHD interpolators. Their initialization
! is controlled by internal key LSLHDQUAD.

LSLHDQUAD=(SLHDKMIN /= 0.0_JPRB) .OR. &
 & (LSLHD.AND.(.NOT.LSLHD_OLD).AND.(SLHDKMAX /= 0.0_JPRB))

!*       2.5   Additional settings

LNHX=LNHDYN.AND.(NVDVAR==4)
LNHXDER=LNHX

!*       2.6   Some tests

!* --- Tests on LGWADV ---

IF (LGWADV.AND.(.NOT.LSLAG)) THEN
  ! * Use Eulerian scheme when LGWADV=true?
  IERR=IERR+1
  WRITE(KULOUT,'(1X,A,I3,A)') ' ! SUDYNA: ERROR NR ',IERR,' !!!'
  WRITE(KULOUT,*) ' LGWADV=true is not coded in the Eulerian dynamics.'
ENDIF
IF (LGWADV.AND.LSLAG.AND.(.NOT.LTWOTL)) THEN
  ! * Use SL3TL scheme when LGWADV=true?
  IERR=IERR+1
  WRITE(KULOUT,'(1X,A,I3,A)') ' ! SUDYNA: ERROR NR ',IERR,' !!!'
  WRITE(KULOUT,*) ' LGWADV=true is not coded in the SL3TL dynamics.'
  ! ky: in this case additional interpolations are required
  !     for the linear terms of "dver" equation, and they are not coded.
ENDIF

IF(LGWADV .AND. LRUBC) THEN
  ! * these two options actually can run together on an informatic
  !   point of view, but the current assumptions done with LGWADV=T
  !   are that a particle which is on the top of the atmosphere remains
  !   on the top of the atmosphere;
  !   this is equivalent to assume that "etadot_top" is always zero.
  !   This condition is satisfied when LRUBC=F
  !   but not when LRUBC=T which gives a non-zero value to "etadot_top".
  !   So one forbids combination "LGWADV .AND. LRUBC" which leads
  !   to inconsistencies in the model at the top of the atmosphere.
  IERR=IERR+1
  WRITE(KULOUT,'(1X,A,I3,A)') ' ! SUDYNA: ERROR NR ',IERR,' !!!'
  WRITE(KULOUT,*) ' LRUBC=true and LGWADV=true are inconsistent.'
ENDIF

!* --- Tests on NGWADVSI ---
IF (NGWADVSI<1 .OR. NGWADVSI>2) THEN
  IERR=IERR+1
  WRITE(KULOUT,'(1X,A,I3,A)') ' ! SUDYNA: ERROR NR ',IERR,' !!!'
  WRITE(KULOUT,*) ' NGWADVSI should be 1 or 2.'
ENDIF

!* --- Tests on LRDBBC ---
IF ( LRDBBC.AND.LGWADV ) THEN
  ! * Use of lrdbbc=true with lgwadv=true makes no sense!'
  IERR=IERR+1
  WRITE(KULOUT,'(1X,A,I3,A)') ' ! SUDYNA: ERROR NR ',IERR,' !!!'
  WRITE(KULOUT,*) ' Options LGWADV=T and LRDBBC=T make no sense together!'
ENDIF
IF (LRDBBC.AND.(.NOT.LSLAG)) THEN
  ! * Use Eulerian scheme when LRDBBC=true?
  IERR=IERR+1
  WRITE(KULOUT,'(1X,A,I3,A)') ' ! SUDYNA: ERROR NR ',IERR,' !!!'
  WRITE(KULOUT,*) ' LRDBBC=true is the S.-L. discretization of BBC'
ENDIF

IF( .NOT.LELAM )THEN
  LLGRID=YDDIM%NSMAX > (YDDIM%NDLON+3)/3
ELSE
  LLGRID=YDDIM%NSMAX > (YDDIM%NDGLG -1)/3 .OR. YDDIM%NMSMAX > (YDDIM%NDLON -1)/3
ENDIF
IF ((NCONF == 1.OR.NCONF == 302).AND.LLGRID.AND..NOT.LSLAG) THEN
  WRITE(KULOUT,*) 'WARNING IN SUDYNA: Linear grid requires SL advection to ensure stability'
ENDIF

!* --- Tests on SLHD variables ---

! Use SLHD scheme with Eulerian advection?
IF (LSLHD.AND.(.NOT.LSLAG)) THEN
  WRITE(KULOUT,*) ' SUDYNA: Warning!'
  WRITE(KULOUT,*) ' SLHD scheme acts only with SL advection!'
ENDIF

! Use SLHD in Arpege?
IF (LSLHD.AND.(.NOT.(LELAM.OR.LECMWF))) THEN
  WRITE(KULOUT,*) ' SUDYNA: Warning!'
  WRITE(KULOUT,*) '   The setup of SLHD was never adapted to suite the setup of stretched geometry.'
ENDIF

! Use nonzero SLHDKMIN outside of SLHD scheme?
IF ((.NOT.LSLHD).AND.(SLHDKMIN /= 0)) THEN
  WRITE(KULOUT,*) ' SUDYNA: Warning!'
  WRITE(KULOUT,*) '   Nonzero SLHDKMIN redefines high order SL interpolator'
  WRITE(KULOUT,*) '   even if SLHD scheme is off.'
ENDIF

! Use SLHD with SLHDKMIN > SLHDKMAX?
IF ((LSLHD.AND.(.NOT.LSLHD_STATIC)).AND.(SLHDKMIN > SLHDKMAX)) THEN
  IERR=IERR+1
  WRITE(KULOUT,'(1X,A,I3,A)') ' ! SUDYNA: ERROR NR ',IERR,' !!!'
  WRITE(KULOUT,*) ' SLHDKMIN must not be greater than SLHDKMAX in diffusion mode!'
ENDIF

! Is SLHDKMIN out of bounds?
IF ((SLHDKMIN < -2.0_JPRB ).OR.(SLHDKMIN > 12.0_JPRB)) THEN
  IERR=IERR+1
  WRITE(KULOUT,'(1X,A,I3,A)') ' ! SUDYNA: ERROR NR ',IERR,' !!!'
  WRITE(KULOUT,*) ' SLHDKMIN is restricted to interval <-2,12>.'
ENDIF

! Is SLHDKMAX, SLHDKREF out of bounds?
IF (LSLHD_OLD) THEN
  IF ((SLHDKMAX < 0.0_JPRB ).OR.(SLHDKMAX > 5.0_JPRB)) THEN
    IERR=IERR+1
    WRITE(KULOUT,'(1X,A,I3,A)') ' ! SUDYNA: ERROR NR ',IERR,' !!!'
    WRITE(KULOUT,*) ' SLHDKMAX is restricted to interval <0,5>.'
  ENDIF
  IF ((SLHDKREF < 0.0_JPRB ).OR.(SLHDKREF > 5.0_JPRB)) THEN
    IERR=IERR+1
    WRITE(KULOUT,'(1X,A,I3,A)') ' ! SUDYNA: ERROR NR ',IERR,' !!!'
    WRITE(KULOUT,*) ' SLHDKREF is restricted to interval <0,5>.'
  ENDIF
ELSE
  IF ((SLHDKMAX < -2.0_JPRB ).OR.(SLHDKMAX > 12.0_JPRB)) THEN
    IERR=IERR+1
    WRITE(KULOUT,'(1X,A,I3,A)') ' ! SUDYNA: ERROR NR ',IERR,' !!!'
    WRITE(KULOUT,*) ' SLHDKMAX is restricted to interval <-2,12>.'
  ENDIF
  IF (SLHDKMAX > 6.0_JPRB) THEN
    WRITE(KULOUT,'(1X,A)') ' ! SUDYNA: WARNING !!!'
    WRITE(KULOUT,*) ' SLHDKMAX greater than 6. may influence your model stability.'
  ENDIF
  IF ((SLHDKREF < -2.0_JPRB ).OR.(SLHDKREF > 12.0_JPRB)) THEN
    IERR=IERR+1
    WRITE(KULOUT,'(1X,A,I3,A)') ' ! SUDYNA: ERROR NR ',IERR,' !!!'
    WRITE(KULOUT,*) ' SLHDKREF is restricted to interval <-2,12>.'
  ENDIF
ENDIF

! Is SLHDEPSH out of bounds?
IF ((SLHDEPSH < 0.0_JPRB ).OR.(SLHDEPSH > 0.5_JPRB )) THEN
  IERR=IERR+1
  WRITE(KULOUT,'(1X,A,I3,A)') ' ! SUDYNA: ERROR NR ',IERR,' !!!'
  WRITE(KULOUT,*) &
   & ' SLHDEPSH must be within interval <0,0.5>.'
ENDIF

! Is SLHDEPSV out of bounds?
IF ((SLHDEPSV < 0.0_JPRB ).OR.(SLHDEPSV > 0.5_JPRB )) THEN
  IERR=IERR+1
  WRITE(KULOUT,'(1X,A,I3,A)') ' ! SUDYNA: ERROR NR ',IERR,' !!!'
  WRITE(KULOUT,*) &
   & ' SLHDEPSV must be within interval <0,0.5>.'
ENDIF

! Is SLVFEPS out of bounds?
IF ((SLVFEPS < 0.0_JPRB ).OR.(SLVFEPS > 0.5_JPRB )) THEN
  IERR=IERR+1
  WRITE(KULOUT,'(1X,A,I3,A)') ' ! SUDYNA: ERROR NR ',IERR,' !!!'
  WRITE(KULOUT,*) &
   & ' SLVFEPS must be within interval <0,0.5>.'
ENDIF

! Is P_SLVF_BOT out of bounds?
IF (P_SLVF_BOT < 0.0_JPRB ) THEN
  IERR=IERR+1
  WRITE(KULOUT,'(1X,A,I3,A)') ' ! SUDYNA: ERROR NR ',IERR,' !!!'
  WRITE(KULOUT,*) &
   & ' P_SLVF_BOT cannot be negative.'
ENDIF

!* --- Tests on COMAD variables ---
! Use LCOMAD with Eulerian advection?
IF (LCOMAD.AND.(.NOT.LSLAG)) THEN
  IERR=IERR+1
  WRITE(KULOUT,'(1X,A,I3,A)') ' ! SUDYNA: ERROR NR ',IERR,' !!!'
  WRITE(KULOUT,*) ' LCOMAD acts only with SL advection!'
ENDIF
! LCOMADH or LCOMADV has to be .T. if LCOMAD=.T.
IF (LCOMAD .AND.((.NOT. LCOMADH) .AND. (.NOT. LCOMADV)) ) THEN
  IERR=IERR+1
  WRITE(KULOUT,'(1X,A,I3,A)') ' ! SUDYNA: ERROR NR ',IERR,' !!!'
  WRITE(KULOUT,*) 'LCOMADH or LCOMADV has to be .T. if LCOMAD=.T.'
ENDIF
! LCOMADH and LCOMADV must to be .F. if LCOMAD=.F.
IF ((.NOT. LCOMAD) .AND.(LCOMADH .OR. LCOMADV) ) THEN
  IERR=IERR+1
  WRITE(KULOUT,'(1X,A,I3,A)') ' ! SUDYNA: ERROR NR ',IERR,' !!!'
  WRITE(KULOUT,*) 'LCOMADH and LCOMADV have to be .F. if LCOMAD=.F.'
ENDIF
! must not have COMAD_x and SLHD_x together
IF ((LCOMAD_W .AND. LSLHD_W).OR.(LCOMAD_T .AND. LSLHD_T) &
 & .OR.(LCOMAD_SPD .AND. LSLHD_SPD).OR.(LCOMAD_SVD .AND. LSLHD_SVD) ) THEN
  IERR=IERR+1
  WRITE(KULOUT,'(1X,A,I3,A)') ' ! SUDYNA: ERROR NR ',IERR,' !!!'
  WRITE(KULOUT,*) 'COMAD and SLHD can not be used at the same time for the same variable'
ENDIF
! must not have COMAD_t and weno together
IF (LCOMAD_T .AND. LRHSVWENO) THEN
  IERR=IERR+1
  WRITE(KULOUT,'(1X,A,I3,A)') ' ! SUDYNA: ERROR NR ',IERR,' !!!'
  WRITE(KULOUT,*) 'COMAD and WENO can not be used at the same time for T'
ENDIF


!* --- Test on WENO ---

IF (LSLTVWENO .AND. (.NOT.LHOISLT)) THEN
  IERR=IERR+1
  WRITE(KULOUT,'(1X,A,I3,A)') ' ! SUDYNA: ERROR NR ',IERR,' !!!'
  WRITE(KULOUT,*) 'LSLTVWENO=.T. MAKES ONLY SENSE TO BE USED WITH LHOISLT=T.'
ENDIF
IF (LSLTVWENO .AND. (HOISLTV /= 0._JPRB)) THEN
  IERR=IERR+1
  WRITE(KULOUT,'(1X,A,I3,A)') ' ! SUDYNA: ERROR NR ',IERR,' !!!'
  WRITE(KULOUT,*) 'WENO assumes the basic interpolator to be at least of 3rd order accuracy.'
  WRITE(KULOUT,*) 'With LSLTVWENO=.T. only HOISLTV=0. ensures then good results.'
ENDIF
IF ((LSLTVWENO .OR. LHOISLT).AND.LRPLANE) THEN
  ! Not yet coded case
  IERR=IERR+1
  WRITE(KULOUT,'(1X,A,I3,A)') ' ! SUDYNA: ERROR NR ',IERR,' !!!'
  WRITE(KULOUT,*) 'LSLTVWENO OR LHOISLT not yet coded for LRPLANE=.T.'
ENDIF




!* --- Tests on ND4SYS ---
IF ((LSLAG.AND.LNHDYN.AND.(NVDVAR == 4)) .AND. .NOT.(ND4SYS == 1 .OR. ND4SYS == 2)) THEN
  IERR=IERR+1
  WRITE(KULOUT,'(1X,A,I3,A)') ' ! SUDYNA: ERROR NR ',IERR,' !!!'
  WRITE(KULOUT,*) ' Allowed values for ND4SYS are 1 or 2.'
ENDIF

!* --- Tests on LSLDIA ---

IF( LSLDIA.AND.LRPLANE ) THEN
  ! Not yet coded case.
  IERR=IERR+1
  WRITE(KULOUT,'(1X,A,I3,A)') ' ! SUDYNA: ERROR NR ',IERR,' !!!'
  WRITE(KULOUT,*) ' LSLDIA is not coded for LRPLANE=T.'
ENDIF

!* --- Tests on LRALTVDISP ---
IF( LRALTVDISP.AND.(.NOT.LSLAG) ) THEN
  IERR=IERR+1
  WRITE(KULOUT,'(1X,A,I3,A)') ' ! SUDYNA: ERROR NR ',IERR,' !!!'
  WRITE(KULOUT,*) ' LRALTVDISP must be F for Eulerian advection scheme.'
ENDIF

!* --- Tests on SETTLS and NESC keys ---


! * Keys NESC.. can be T only for SL2TL advection.
IF ( (LNESCT.OR.LNESCV.OR.LNESC).AND..NOT.LTWOTL ) THEN
  IERR=IERR+1
  WRITE(KULOUT,'(1X,A,I3,A)') ' ! SUDYNA: ERROR NR ',IERR,' !!!'
  WRITE(KULOUT,*) ' LTWOTL=F requires LNESCT=LNESCV=LNESC=F '
ENDIF

! * Keys NESC and SETTLS are exclusive.
IF( LNESCT.AND.LSETTLST )THEN
  IERR=IERR+1
  WRITE(KULOUT,'(1X,A,I3,A)') ' ! SUDYNA: ERROR NR ',IERR,' !!!'
  WRITE(KULOUT,*) ' LNESCT and LSETTLST cannot be both TRUE'
ENDIF
IF( LNESCV.AND.LSETTLSV )THEN
  IERR=IERR+1
  WRITE(KULOUT,'(1X,A,I3,A)') ' ! SUDYNA: ERROR NR ',IERR,' !!!'
  WRITE(KULOUT,*) ' LNESCV and LSETTLSV cannot be both TRUE'
ENDIF
IF( LNESC.AND.LSETTLS )THEN
  IERR=IERR+1
  WRITE(KULOUT,'(1X,A,I3,A)') ' ! SUDYNA: ERROR NR ',IERR,' !!!'
  WRITE(KULOUT,*) ' LNESC and LSETTLS cannot be both TRUE'
ENDIF
  
! * Additional testings.
IF (.NOT.LR2D) THEN
  ! * Use .NOT.(LNESCT.OR.LSETTLST)=T (not coded any more)?
  IF (LSLAG.AND.LTWOTL.AND.(.NOT.(LNESCT.OR.LSETTLST))) THEN
    ! * Use SL2TL scheme with LNESCT=F and LSETTLST=F?
    IERR=IERR+1
    WRITE(KULOUT,'(1X,A,I3,A)') ' ! SUDYNA: ERROR NR ',IERR,' !!!'
    WRITE(KULOUT,*) ' LNESCT or LSETTLST must be T in the SL2TL.'
  ENDIF
  ! * Use .NOT.(LNESCV.OR.LSETTLSV)=T (not coded any more)?
  IF (LSLAG.AND.LTWOTL.AND.(.NOT.(LNESCV.OR.LSETTLSV))) THEN
    ! * Use SL2TL scheme with LNESCV=F and LSETTLSV=F?
    IERR=IERR+1
    WRITE(KULOUT,'(1X,A,I3,A)') ' ! SUDYNA: ERROR NR ',IERR,' !!!'
    WRITE(KULOUT,*) ' LNESCV or LSETTLSV must be T in the SL2TL.'
  ENDIF
  ! * Use .NOT.(LNESC.OR.LSETTLS)=T when not coded?
  IF (LGWADV.AND.LSLAG.AND.LTWOTL.AND.(.NOT.(LNESC.OR.LSETTLS))) THEN
    ! * Use SL2TL scheme with LNESC=F and LSETTLS=F when LGWADV=true?
    IERR=IERR+1
    WRITE(KULOUT,'(1X,A,I3,A)') ' ! SUDYNA: ERROR NR ',IERR,' !!!'
    WRITE(KULOUT,*) ' LGWADV=T requires LNESC=T or LSETTLS=T in the SL2TL.'
  ENDIF
ENDIF 

!* --- Tests on LPC_FULL and LPC_CHEAP keys ---
  
! *  Use LPC_CHEAP when not available?
IF ( (LPC_CHEAP.OR.LPC_CHEAP2).AND..NOT.LPC_FULL ) THEN
  IERR=IERR+1
  WRITE(KULOUT,'(1X,A,I3,A)') ' ! SUDYNA: ERROR NR ',IERR,' !!!'
  WRITE(KULOUT,*) ' LPC_CHEAP (cheap PC_FULL scheme) HAS MEANING ONLY WITH PC SCHEME LPC_FULL'
ENDIF
IF ( (LPC_CHEAP.OR.LPC_CHEAP2).AND..NOT.LSLAG ) THEN
  IERR=IERR+1
  WRITE(KULOUT,'(1X,A,I3,A)') ' ! SUDYNA: ERROR NR ',IERR,' !!!'
  WRITE(KULOUT,*) ' LPC_CHEAP IS CODED ONLY IN THE SL SCHEME'
ENDIF
IF ( LPC_CHEAP.AND.LPC_CHEAP2 ) THEN
  IERR=IERR+1
  WRITE(KULOUT,'(1X,A,I3,A)') ' ! SUDYNA: ERROR NR ',IERR,' !!!'
  WRITE(KULOUT,*) ' LPC_CHEAP and LPC_CHEAP2 cannot be both T '
ENDIF

! * Use LPC_FULL for SL3TL (not coded)?
IF ( LPC_FULL.AND.(.NOT.LTWOTL.AND.LSLAG) ) THEN
  IERR=IERR+1
  WRITE(KULOUT,'(1X,A,I3,A)') ' ! SUDYNA: ERROR NR ',IERR,' !!!'
  WRITE(KULOUT,*) ' LPC_FULL is not yet coded for SL3TL scheme'
ENDIF

! * Use LRUBC=T without LPC_FULL (not allowed)?
IF ( LRUBC .AND. .NOT.LPC_FULL ) THEN
  IERR=IERR+1
  WRITE(KULOUT,'(1X,A,I3,A)') ' ! SUDYNA: ERROR NR ',IERR,' !!!'
  WRITE(KULOUT,*) ' When LRUBC=true, LPC_FULL must be T'
ENDIF

!* --- Tests on LRUBC ---

! * Use LRUBC in spherical geometry (not allowed)?
IF (LRUBC .AND. .NOT.LELAM) THEN
  IERR=IERR+1
  WRITE(KULOUT,'(1X,A,I3,A)') ' ! SUDYNA: ERROR NR ',IERR,' !!!'
  WRITE(KULOUT,*) ' LRUBC=T can be used with LAM models only'
ENDIF

! * Use LRUBC with improper options for vertical discretisations?
IF ( LRUBC .AND. (LDVERTFE .OR. KDLNPR /= 1) ) THEN
  IERR=IERR+1
  WRITE(KULOUT,'(1X,A,I3,A)') ' ! SUDYNA: ERROR NR ',IERR,' !!!'
  WRITE(KULOUT,*) ' LRUBC=T requires FD vertical discretisation and NDLNPR=1 '
ENDIF

!* --- Tests on L3DTURB ---

IF (L3DTURB.AND.(.NOT.LSLAG)) THEN
  IERR=IERR+1
  WRITE(KULOUT,'(1X,A,I3,A)') ' ! SUDYNA: ERROR NR ',IERR,' !!!'
  WRITE(KULOUT,*) ' L3DTURB=.T. works only with SL advection!'
ENDIF
IF (L3DTURB.AND.(NCONF /= 1)) THEN
  IERR=IERR+1
  WRITE(KULOUT,'(1X,A,I3,A)') ' ! SUDYNA: ERROR NR ',IERR,' !!!'
  WRITE(KULOUT,*) ' L3DTURB=.T. works only with configuration 001'
ENDIF


IF (IERR >= 1) THEN
  CALL FLUSH(KULOUT)
  CALL ABOR1(' SUDYNA: ABOR1 CALLED')
ENDIF




!*      2.9   LSLINL.. keys

! Separate linear from other terms in continuity equation
LSLINLC2=LNHDYN.AND.(LTWOTL.AND.LGWADV.AND.LSETTLS)
! Separate linear from other terms in other equations
LSLINL=LNHDYN.AND.(LTWOTL.AND.LGWADV.AND.LSETTLS)

!*      2.10   LMIXETTLS for restricted choices

IF (LMIXETTLS) THEN
  IF (.NOT.LELAM) THEN
    CALL ABOR1(' SUDYNA: LMIXETTLS CODED ONLY FOR LELAM = .TRUE.')
  ENDIF
  IF (.NOT.LNHEE) THEN
    CALL ABOR1(' SUDYNA: LMIXETTLS CODED ONLY FOR LNHEE = .TRUE.')
  ENDIF
  IF (ND4SYS /= 2) THEN
    CALL ABOR1(' SUDYNA: LMIXETTLS CODED ONLY FOR ND4SYS = 2')
  ENDIF
  IF (.NOT.LGWADV) THEN
    CALL ABOR1(' SUDYNA: LMIXETTLS CODED ONLY FOR LGWADV = .TRUE.')
  ENDIF
  IF (.NOT.LSETTLS) THEN
    CALL ABOR1(' SUDYNA: LMIXETTLS REQUIRES LSETTLS = .TRUE.')
  ENDIF
ENDIF

!*      2.11   LNHQE_C2

! Is constraint C2 matched for NHQE model?
LNHQE_C2=LNHQE.AND.(.NOT.LDVERTFE).AND.(KDLNPR==2)

!*      2.12   LNHQE_SIHYD

IF (.NOT.LNHQE) LNHQE_SIHYD=.FALSE.

! Set up model dependent INTDYN structures
CALL SUPTR_TTND(YDDYNA%YYTTND,LNHEE,LNHQE)
CALL SUPTR_TGMVT(YDDYNA%YYTGMVT95,LNHDYN)
CALL SUPTR_TGFLT(YDDYNA%YYTGFLT95)
!     ------------------------------------------------------------------

!*       3.    Print final values.
!              -------------------
IF(JPRB==JPRD)THEN
   WRITE(UNIT=KULOUT,FMT='('' The model is using double precision. '')') 
ELSE
   WRITE(UNIT=KULOUT,FMT='('' The model is using single precision. '')') 
ENDIF
WRITE(UNIT=KULOUT,FMT='('' LNHDYN = '',L2)') LNHDYN
WRITE(UNIT=KULOUT,FMT='('' LNHEE = '',L2,'' LNHQE = '',L2)') LNHEE,LNHQE
WRITE(UNIT=KULOUT,FMT='('' Printings of YOMDYNA/NAMDYNA variables '')')
WRITE(UNIT=KULOUT,FMT='('' LSPRT = '',L2,'' LTWOTL = '',L2,'' LSLAG = '',L2)') LSPRT,LTWOTL,LSLAG
WRITE(UNIT=KULOUT,FMT='('' LSACC = '',L2)') LSACC
WRITE(UNIT=KULOUT,FMT='('' LGWADV = '',L2)') LGWADV
WRITE(UNIT=KULOUT,FMT='('' NGWADVSI = '',I2)') NGWADVSI
WRITE(UNIT=KULOUT,FMT='('' LRDBBC = '',L2)') LRDBBC
WRITE(UNIT=KULOUT,FMT='('' NPDVAR = '',I2)') NPDVAR
WRITE(UNIT=KULOUT,FMT='('' NVDVAR = '',I2)') NVDVAR
WRITE(UNIT=KULOUT,FMT='('' LSLHD_W       = '',L2)') LSLHD_W
WRITE(UNIT=KULOUT,FMT='('' LSLHD_T       = '',L2)') LSLHD_T
WRITE(UNIT=KULOUT,FMT='('' LSLHD_SPD     = '',L2)') LSLHD_SPD
WRITE(UNIT=KULOUT,FMT='('' LSLHD_SVD     = '',L2)') LSLHD_SVD
WRITE(UNIT=KULOUT,FMT='('' LSLHD_GFL     = '',L2)') LSLHD_GFL
WRITE(UNIT=KULOUT,FMT='('' LSLHD         = '',L2)') LSLHD
WRITE(UNIT=KULOUT,FMT='('' LSLHDQUAD     = '',L2)') LSLHDQUAD
WRITE(UNIT=KULOUT,FMT='('' LSLHD_OLD     = '',L2)') LSLHD_OLD
WRITE(UNIT=KULOUT,FMT='('' LSLHD_STATIC  = '',L2)') LSLHD_STATIC
WRITE(UNIT=KULOUT,FMT='('' SLHDKMIN      = '',E10.4)') SLHDKMIN
WRITE(UNIT=KULOUT,FMT='('' SLHDKMAX      = '',E10.4)') SLHDKMAX
WRITE(UNIT=KULOUT,FMT='('' SLHDKREF      = '',E10.4)') SLHDKREF
WRITE(UNIT=KULOUT,FMT='('' SLHDEPSH      = '',E10.4)') SLHDEPSH
WRITE(UNIT=KULOUT,FMT='('' SLHDEPSV      = '',E10.4)') SLHDEPSV
WRITE(UNIT=KULOUT,FMT='('' LSLVF         = '',L2)') LSLVF
WRITE(UNIT=KULOUT,FMT='('' P_SLVF_BOT    = '',E14.8)') P_SLVF_BOT
WRITE(UNIT=KULOUT,FMT='('' SLVFEPS       = '',E10.4)') SLVFEPS
WRITE(UNIT=KULOUT,FMT='('' SLVFMIN       = '',E14.8)') SLVFMIN 
WRITE(UNIT=KULOUT,FMT='('' SLVFTR        = '',E14.8)') SLVFTR
WRITE(UNIT=KULOUT,FMT='('' P_SLVF_TR     = '',E14.8)') P_SLVF_TR
WRITE(UNIT=KULOUT,FMT='('' LCOMAD        = '',L2)') LCOMAD
WRITE(UNIT=KULOUT,FMT='('' LCOMADH       = '',L2)') LCOMADH
WRITE(UNIT=KULOUT,FMT='('' LCOMADV       = '',L2)') LCOMADV
WRITE(UNIT=KULOUT,FMT='('' LCOMAD_W      = '',L2)') LCOMAD_W
WRITE(UNIT=KULOUT,FMT='('' LCOMAD_T      = '',L2)') LCOMAD_T
WRITE(UNIT=KULOUT,FMT='('' LCOMAD_SPD    = '',L2)') LCOMAD_SPD
WRITE(UNIT=KULOUT,FMT='('' LCOMAD_SVD    = '',L2)') LCOMAD_SVD
WRITE(UNIT=KULOUT,FMT='('' LCOMAD_SP     = '',L2)') LCOMAD_SP 
WRITE(UNIT=KULOUT,FMT='('' LCOMAD_GFL    = '',L2)') LCOMAD_GFL
WRITE(UNIT=KULOUT,FMT='('' ND4SYS = '',I1)') ND4SYS
WRITE(UNIT=KULOUT,FMT='('' LNHX = '',L2)') LNHX
WRITE(UNIT=KULOUT,FMT='('' LNHXDER = '',L2)') LNHXDER
WRITE(UNIT=KULOUT,FMT='('' LSLDIA = '',L2)') LSLDIA
WRITE(UNIT=KULOUT,FMT='('' LGRADSP = '',L2)') LGRADSP
WRITE(UNIT=KULOUT,FMT='('' LPGFSF = '',L2)') LPGFSF
WRITE(UNIT=KULOUT,FMT='('' LPGREUSE = '',L2)') LPGREUSE
WRITE(UNIT=KULOUT,FMT='('' L3DTURB = '',L2)') L3DTURB
WRITE(UNIT=KULOUT,FMT='('' LHOISLT = '',L2)') LHOISLT
WRITE(UNIT=KULOUT,FMT='('' HOISLTH = '',E10.4,'' HOISLTV = '',E10.4)') &
 & HOISLTH,HOISLTV
WRITE(UNIT=KULOUT,FMT='('' LSLTVWENO = '',L2)') LSLTVWENO
WRITE(UNIT=KULOUT,FMT='('' LRHSVWENO = '',L2)') LRHSVWENO
WRITE(UNIT=KULOUT,FMT='('' LRUBC = '',L2)') LRUBC
WRITE(UNIT=KULOUT,FMT='('' LRALTVDISP = '',L2)') LRALTVDISP
WRITE(UNIT=KULOUT,FMT='('' LRPRSLTRJ = '',L2)') LRPRSLTRJ
WRITE(UNIT=KULOUT,FMT='('' LPC_FULL = '',L2)') LPC_FULL
WRITE(UNIT=KULOUT,FMT='('' LPC_CHEAP = '',L2)') LPC_CHEAP
WRITE(UNIT=KULOUT,FMT='('' LPC_CHEAP2 = '',L2)') LPC_CHEAP2
WRITE(UNIT=KULOUT,FMT='('' LNESCT = '',L2)') LNESCT
WRITE(UNIT=KULOUT,FMT='('' LNESCV = '',L2)') LNESCV
WRITE(UNIT=KULOUT,FMT='('' LNESC = '',L2)') LNESC
WRITE(UNIT=KULOUT,FMT='('' LSETTLST = '',L2)') LSETTLST
WRITE(UNIT=KULOUT,FMT='('' LSETTLSV = '',L2)') LSETTLSV
WRITE(UNIT=KULOUT,FMT='('' LSETTLS = '',L2)') LSETTLS
WRITE(UNIT=KULOUT,FMT='('' LELTRA = '',L2)') LELTRA
WRITE(UNIT=KULOUT,FMT='('' LSLINLC2 = '',L2)') LSLINLC2
WRITE(UNIT=KULOUT,FMT='('' LSLINL = '',L2)') LSLINL
WRITE(UNIT=KULOUT,FMT='('' L_RDRY_VD= '',L2)') L_RDRY_VD
WRITE(UNIT=KULOUT,FMT='('' LNHQE_SIHYD = '',L2)') LNHQE_SIHYD
WRITE(UNIT=KULOUT,FMT='('' LSI_NHEE = '',L2)') LSI_NHEE
WRITE(UNIT=KULOUT,FMT='('' LNHEE_SVDLAPL_FIRST = '',L2)') LNHEE_SVDLAPL_FIRST
WRITE(UNIT=KULOUT,FMT='('' LMIXETTLS    = '',L2)') LMIXETTLS
WRITE(UNIT=KULOUT,FMT='('' LMIXETTLS_PRINT   = '',L2)') LMIXETTLS_PRINT
WRITE(UNIT=KULOUT,FMT='('' RMIXNL_TRH  = '',E10.4)') RMIXNL_TRH

!     ------------------------------------------------------------------

!9990 FORMAT(1X,'ARRAY ',A10,' ALLOCATED ',8I8)

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUDYNA',1,ZHOOK_HANDLE)
END SUBROUTINE SUDYNA
