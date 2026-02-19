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

SUBROUTINE SURAND1(YDGEOMETRY,YDML_PHY_STOCH,YDDYN,YDRIP,YDECUCONVCA)

!**** *SURAND1*  - Initialize stochastic physics parameters - part 1

!     Purpose.
!     --------
!           Initialize stochastic physics

!**   Interface.
!     ----------
!        *CALL* *SURAND1

!        Explicit arguments :
!        --------------------
!        None

!        Implicit arguments :
!        --------------------
!        COMMON STOPH_MIX

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!        None

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      R Buizza and L Isaksen - ECMWF - 
!      Original : 98-01-01

!     Modifications.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      J.Berner 06-08-11 added stoch phys options for cellular automaton 
!                        backscatter and  spectral AR backscatter 
!      G. Shutts 09-03-4 spectral backscatter update
!      M. Steinheimer 09-03-05 added option RVP
!      M. Steinheimer 09-10-27 spectral backscatter update (changes to RVP,
!                        added option LSTOPH_SPBS_FAST)
!      P. Bechtold    14/05/2012 replace 86400 by RDAY
!      M. Leutbecher  13-Mar-2013 introduced LSPBS_DISSGW
!      T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!      K. Yessad (July 2014): Move some variables.
!      M. Leutbecher  16-Oct-2014 cubic grid diffusion for SKEB (LSTOPH_SPBS=T) and LSPECVIS
!      SJ Lock:          Jan-2016 Removed LSTOPH option
!     ------------------------------------------------------------------

USE MODEL_PHYSICS_STOCHAST_MOD , ONLY : MODEL_PHYSICS_STOCHAST_TYPE
USE GEOMETRY_MOD               , ONLY : GEOMETRY
USE PARKIND1                   , ONLY : JPIM, JPRB, JPRD
USE YOMHOOK                    , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCT0                     , ONLY : LALLOPR, LSPBSBAL
USE YOMCST                     , ONLY : RA, RDAY
USE YOMDYN                     , ONLY : TDYN
USE YOMVERT                    , ONLY : VP00
USE YOMGRIB                    , ONLY : NENSFNB 
USE YOMLUN                     , ONLY : NULOUT, NULNAM, NULERR, RESERVE_LUN, FREE_LUN
USE YOMMP0                     , ONLY : NPROC, MYPROC, NPRINTLEV
USE YOMMSC                     , ONLY : N_DEFAULT_REAL_KIND, N_DOUBLE_KIND
USE YOMTAG                     , ONLY : MTAGFCE 
USE YOMRIP                     , ONLY : TRIP
USE MPL_MODULE                 , ONLY : MPL_BROADCAST
USE STOPH_MIX                  , ONLY : SIGTOFACT,INITIALIZE_CELLS
USE YOE_CUCONVCA               , ONLY : TECUCONVCA
!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)                   , INTENT(INOUT) :: YDGEOMETRY
TYPE(MODEL_PHYSICS_STOCHAST_TYPE), INTENT(INOUT), TARGET :: YDML_PHY_STOCH
TYPE(TDYN)                       , INTENT(INOUT) :: YDDYN
TYPE(TRIP)                       , INTENT(INOUT) :: YDRIP
TYPE(TECUCONVCA)                 , INTENT(INOUT) :: YDECUCONVCA
REAL(KIND=JPRB) ,DIMENSION(0:YDGEOMETRY%YRDIM%NSMAX) :: ZCHI

INTEGER(KIND=JPIM) :: IN,IFILT,ICUTOFF,IIP,JJP
INTEGER(KIND=JPIM) :: JN, I, IU, ISPINUP, ISTEP
INTEGER(KIND=JPIM) :: INLEVS,  JX, JY, ITAG, INFO, IOMASTER, J1, J2, JLEV
INTEGER(KIND=JPIM) :: IULTMP

LOGICAL :: LLP
CHARACTER(LEN=10) :: CLID
CHARACTER(LEN=40) :: CLFILE

REAL(KIND=JPRB) :: ZCORRMAT(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZSQRTCORR2(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZGAMMAN,ZCONSTF0
REAL(KIND=JPRB) :: ZDSIGMA
REAL(KIND=JPRB) :: ZSIGMA
REAL(KIND=JPRB) :: ZSIGMA_SMOOTH
REAL(KIND=JPRB) :: ZAFILT
REAL(KIND=JPRB) :: ZP
REAL(KIND=JPRB) :: ZX
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

REAL(KIND=JPRB) :: ZSIGMA_NFRSPBS
REAL(KIND=JPRB) :: ZRVP_MULFACT
INTEGER(KIND=JPIM) :: IRVP_MULNSMAX
INTEGER(KIND=JPIM) :: INSMAXEFF
!     ------------------------------------------------------------------

LOGICAL, POINTER :: LSTOPH_SPBS, LSTOPH_CASBS, LVORTCON, LEXTRAFIELDS,&
 & LSTOPH_JBCOR, LSTOPH_UNCORR, LSTOPH_UNIFORM, LFORCENL, LSTOPH_RVP,&
 & LSTOPH_TAPER, LSTOPH_INI, LSTOPH_VARALPHA, LSTOPH_SPBS_FAST,&
 & LSTOPH_SPBS_VORT, LSTOPH_SPBS_T, LSPBS_DISSGW, LSTOPH_GAUSS, LSPBSNORM,&
 & LSPBS_DISSNUM, LSPBS_DISSCU, LSPBS_DISSNUM_CT,&
 & LSTOPH_RVPOLD, LSPBSDISS

INTEGER(KIND=JPIM), POINTER :: NSTOCHOPT, NFORCESTART, NFORCEEND,&
 & NFRSTOPH_SPBS, NFRSTOPH_VC, NSMAXSPBS

REAL(KIND=JPRB), POINTER :: RATIO_BACKSCAT, RATIO_BACKSCAT_CON2NUM, AMAGSTOPH,&
 & ADLATSTOPH, ADLONSTOPH, AMAGSTOPH_CASBS, ADLATSTOPH_CA, ADLONSTOPH_CA,&
 & ALPHA_DEEP_CONV, ALPHA_SHAL_CONV, SLDISSFAC, REXPONENT, RFLUX_DET_CLIP,&
 & VC_CON, RVP_MULMIN, RVP_MULMAX, RVP_MULEXP, RVP_MULNSMAX, RVP_MUL_A,&
 & RVP_MUL_B, RVP_MUL_C, RVP_MUL_D, RVP_MUL_1, RVP_MUL_2, RVP_MUL_1_T,&
 & RVP_MUL_2_T, REXPONENT_T, RATIO_APE2KE, RVP_MULMIN_T, RVP_MULMAX_T,&
 & RVP_MULEXP_T, RVP_MULNSMAX_T, RVP_MUL_A_T, RVP_MUL_B_T, RVP_MUL_C_T,&
 & RVP_MUL_D_T, RSMOOTHSCALE, RSPBS_TAU

#include "namstoph.nam.h"

#include "fcttim.func.h"

#include "abor1.intfb.h"
#include "posnam.intfb.h"
#include "setran.intfb.h"
#include "surand2.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SURAND1',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM, &
  & YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, &
  & YDVAB=>YDGEOMETRY%YRVAB, YDVETA=>YDGEOMETRY%YRVETA, YDVFE=>YDGEOMETRY%YRVFE, YDSTOPH=>YDML_PHY_STOCH%YRSTOPH)

ASSOCIATE(LCUCONV_CA=>YDECUCONVCA%LCUCONV_CA, &
 & NDGLG =>YDDIM%NDGLG,  NGPBLKS=>YDDIM%NGPBLKS, NPROMA=>YDDIM%NPROMA, &
 & NSMAX=>YDDIM%NSMAX, NSPEC2=>YDDIM%NSPEC2, &
 & NFLSUR=>YDDIMV%NFLSUR, NFLEVG=>YDDIMV%NFLEVG, &
 & HDIRVOR=>YDDYN%HDIRVOR, &
 & NGPTOT=>YDGEM%NGPTOT, &
 & LHDIFFM =>YDDYN%LHDIFFM, LSPECVIS=>YDDYN%LSPECVIS, NDIFFACT=>YDDYN%NDIFFACT, &
 & TDT=>YDRIP%TDT, &
 & BIHARM=>YDSTOPH%BIHARM, NFRSTOPH_SPBS_PAT=>YDSTOPH%NFRSTOPH_SPBS_PAT, &
 & RSIGMA2_EPS=>YDSTOPH%RSIGMA2_EPS, TAPER_SIGMABOT=>YDSTOPH%TAPER_SIGMABOT, &
 & TAPER_SIGMATOP=>YDSTOPH%TAPER_SIGMATOP, TAPER0=>YDSTOPH%TAPER0, &
 & TAPER1=>YDSTOPH%TAPER1, TAPER2=>YDSTOPH%TAPER2, TAPER3=>YDSTOPH%TAPER3) 
! Associate pointers for variables in namelist
LSTOPH_SPBS            => YDSTOPH%LSTOPH_SPBS
LSTOPH_CASBS           => YDSTOPH%LSTOPH_CASBS
LVORTCON               => YDSTOPH%LVORTCON
LEXTRAFIELDS           => YDSTOPH%LEXTRAFIELDS
LSTOPH_JBCOR           => YDSTOPH%LSTOPH_JBCOR
LSTOPH_UNCORR          => YDSTOPH%LSTOPH_UNCORR
LSTOPH_UNIFORM         => YDSTOPH%LSTOPH_UNIFORM
LFORCENL               => YDSTOPH%LFORCENL
NSTOCHOPT              => YDSTOPH%NSTOCHOPT
NFORCESTART            => YDSTOPH%NFORCESTART
NFORCEEND              => YDSTOPH%NFORCEEND
NFRSTOPH_SPBS          => YDSTOPH%NFRSTOPH_SPBS
NFRSTOPH_VC            => YDSTOPH%NFRSTOPH_VC
RATIO_BACKSCAT         => YDSTOPH%RATIO_BACKSCAT
RATIO_BACKSCAT_CON2NUM => YDSTOPH%RATIO_BACKSCAT_CON2NUM
AMAGSTOPH_CASBS        => YDSTOPH%AMAGSTOPH_CASBS
ADLATSTOPH_CA          => YDSTOPH%ADLATSTOPH_CA
ADLONSTOPH_CA          => YDSTOPH%ADLONSTOPH_CA
ALPHA_DEEP_CONV        => YDSTOPH%ALPHA_DEEP_CONV
ALPHA_SHAL_CONV        => YDSTOPH%ALPHA_SHAL_CONV
SLDISSFAC              => YDSTOPH%SLDISSFAC
REXPONENT              => YDSTOPH%REXPONENT
RFLUX_DET_CLIP         => YDSTOPH%RFLUX_DET_CLIP
VC_CON                 => YDSTOPH%VC_CON
LSTOPH_RVP             => YDSTOPH%LSTOPH_RVP
LSTOPH_TAPER           => YDSTOPH%LSTOPH_TAPER
LSTOPH_INI             => YDSTOPH%LSTOPH_INI
RVP_MULMIN             => YDSTOPH%RVP_MULMIN
RVP_MULMAX             => YDSTOPH%RVP_MULMAX
RVP_MULEXP             => YDSTOPH%RVP_MULEXP
RVP_MULNSMAX           => YDSTOPH%RVP_MULNSMAX
RVP_MUL_A              => YDSTOPH%RVP_MUL_A
RVP_MUL_B              => YDSTOPH%RVP_MUL_B
RVP_MUL_C              => YDSTOPH%RVP_MUL_C
RVP_MUL_D              => YDSTOPH%RVP_MUL_D
LSTOPH_VARALPHA        => YDSTOPH%LSTOPH_VARALPHA
RVP_MUL_1              => YDSTOPH%RVP_MUL_1
RVP_MUL_2              => YDSTOPH%RVP_MUL_2
RVP_MUL_1_T            => YDSTOPH%RVP_MUL_1_T
RVP_MUL_2_T            => YDSTOPH%RVP_MUL_2_T
LSTOPH_SPBS_FAST       => YDSTOPH%LSTOPH_SPBS_FAST
LSTOPH_SPBS_VORT       => YDSTOPH%LSTOPH_SPBS_VORT
LSTOPH_SPBS_T          => YDSTOPH%LSTOPH_SPBS_T
LSPBS_DISSGW           => YDSTOPH%LSPBS_DISSGW
LSPBS_DISSNUM          => YDSTOPH%LSPBS_DISSNUM
LSPBS_DISSCU           => YDSTOPH%LSPBS_DISSCU
LSPBS_DISSNUM_CT       => YDSTOPH%LSPBS_DISSNUM_CT
REXPONENT_T            => YDSTOPH%REXPONENT_T
RATIO_APE2KE           => YDSTOPH%RATIO_APE2KE
RVP_MULMIN_T           => YDSTOPH%RVP_MULMIN_T
RVP_MULMAX_T           => YDSTOPH%RVP_MULMAX_T
RVP_MULEXP_T           => YDSTOPH%RVP_MULEXP_T
RVP_MULNSMAX_T         => YDSTOPH%RVP_MULNSMAX_T
RVP_MUL_A_T            => YDSTOPH%RVP_MUL_A_T
RVP_MUL_B_T            => YDSTOPH%RVP_MUL_B_T
RVP_MUL_C_T            => YDSTOPH%RVP_MUL_C_T
RVP_MUL_D_T            => YDSTOPH%RVP_MUL_D_T
LSTOPH_GAUSS           => YDSTOPH%LSTOPH_GAUSS
RSMOOTHSCALE           => YDSTOPH%RSMOOTHSCALE
LSPBSNORM              => YDSTOPH%LSPBSNORM
RSPBS_TAU              => YDSTOPH%RSPBS_TAU
LSTOPH_RVPOLD          => YDSTOPH%LSTOPH_RVPOLD
NSMAXSPBS              => YDSTOPH%NSMAXSPBS
LSPBSDISS              => YDSTOPH%LSPBSDISS

!     ------------------------------------------------------------------

!*    1. Set default values.
!     ----------------------

!     NFRSTOPH_SPBS : frequency (number of time steps) between calls to 
!                    refresh the streamfunction forcing in SPBS (must be >= 2)
!                    if CASBS run it must be set to 1
!     NFRSTOPH_VC   :  frequency (number of time steps) between calls to vorticity confinement

!     1.2 Set implicit default values for CASBS and/or SPBS
!     ALPHA_STO       : Autoregressive parameter (function of n) 
!     ALPHA_DEEP_CONV : Entrainment cloud fraction for deep convection
!     ALPHA_SHAL_CONV : Entrainment cloud ascent fraction for shallow and mid-level convection
!     AMAGSTOPH_CASBS : Magnitude of forcing 
!     ADLATSTOPH_CA   : Gridsize of cellular automaton in zonal direction
!     ADLONSTOPH_CA   : Gridsize of cellular automaton in meridional direction
!     RATIO_BACKSCAT  : Backscatter ratio
!     RATIO_BACKSCAT_CON2NUM  : Convective Backscatter ratio/Numerical backscatter ratio
!     LSTOPH_CASBS    : .T. to add cellular automaton streamfunction perturbations
!     LSTOPH_SPBS     : .T. to add spectral backscatter streamfunction perturbations
!     LVORTCON        : .T. to switch on vorticity confinement
!     LEXTRAFIELDS    : .T. to print extrafields
!     LSTOPH_JBCOR    : .T. if perturbations are vertically correlated like Jb vorticity
!     LSTOPH_UNCORR   : .T. to make SPBS perturbations vertically uncorrelated
!     LSTOPH_UNIFORM  : .T. to make SPBS perturbations uniformly distributed, instead
!                        of Gaussian distributed
!     LSTOPH_RVP      : .T. to use random vertical profiles in SPBS
!     LSTOPH_SPBS_FAST: .T. pattern is only updated every NFRSTOPH_SPBS timesteps (statistically equivalent)
!     LSTOPH_VARALPHA : .T. to use scale-dependen decorrelation time in SPBS

!     -----------------------------------------------------------------
!*    1.1 Set default values for forced S.V.s
!     -----------------------------------------------------------------
LFORCENL=.FALSE.
NFORCESTART=0
NFORCEEND=240

!     -----------------------------------------------------------------
!*    1.2 Set default values for spectral backscatter
!     -----------------------------------------------------------------
LEXTRAFIELDS=.FALSE.
LSTOPH_SPBS=.FALSE.
LSTOPH_JBCOR=.FALSE.
LSTOPH_UNCORR=.FALSE.
LSTOPH_UNIFORM=.FALSE.
LSTOPH_RVP=.TRUE.
LSTOPH_RVPOLD=.FALSE.
LSTOPH_TAPER=.TRUE.
LSTOPH_INI=.TRUE.
LSTOPH_VARALPHA=.FALSE.
LSTOPH_SPBS_FAST=.TRUE.
LSTOPH_SPBS_VORT=.FALSE.
LSTOPH_GAUSS=.TRUE.
LSPBS_DISSNUM=.TRUE.
LSPBS_DISSNUM_CT=.FALSE.
LSPBS_DISSCU=.TRUE.
LSPBS_DISSGW=.FALSE.
LSPBSNORM=.TRUE.
LSPBSDISS=.TRUE.

NSMAXSPBS=MIN(159,NSMAX)  !note default setting in scripts/modeleps might be different!

!   REXPONENT produces (3-2*REXPONENT)-kinetic energy forcing spectrum 
REXPONENT=-1.27_JPRB 
RATIO_BACKSCAT= 0.095_JPRB
RATIO_BACKSCAT_CON2NUM= 1.0_JPRB
RFLUX_DET_CLIP= 0.005_JPRB  !  upper limit for convective mass flux detrainment rate
RSPBS_TAU    =  2.5E04_JPRB 
RSIGMA2_EPS=2._JPRB/12.0_JPRB
ISPINUP=1    

! set defaults for temperature backscatter

LSTOPH_SPBS_T=.FALSE.
REXPONENT_T= REXPONENT
RATIO_APE2KE= 0.2

!     -----------------------------------------------------------------
!*    1.2 Set default values for random vertical profiles in spectral backscatter
!     -----------------------------------------------------------------
!for N dependency
RVP_MULMIN=0._JPRB
RVP_MULMAX=1.00_JPRB
RVP_MULEXP=0.10_JPRB
RVP_MULNSMAX=255._JPRB
!for p dependency
RVP_MUL_A=1.6857964_JPRB
RVP_MUL_B=-0.0059995210E-02_JPRB
RVP_MUL_C=1.1070694E-09
RVP_MUL_D=-6.5077043E-15
RVP_MUL_1=18.0192_JPRB
RVP_MUL_2=-0.2159_JPRB

!for T-backscatter
!for N dependency
RVP_MULMIN_T=0._JPRB
RVP_MULMAX_T=1.00_JPRB
RVP_MULEXP_T=0.10_JPRB
RVP_MULNSMAX_T=255._JPRB
!for p dependency
RVP_MUL_A_T=1.6857964_JPRB
RVP_MUL_B_T=-0.0059995210E-02_JPRB
RVP_MUL_C_T=1.1070694E-09
RVP_MUL_D_T=-6.5077043E-15
RVP_MUL_1_T=10.5608_JPRB
RVP_MUL_2_T=-0.2959_JPRB

!     -----------------------------------------------------------------
!*    1.3 Set default values for cellular automaton stochastic backscatter
!     -----------------------------------------------------------------
LSTOPH_CASBS=.FALSE.
AMAGSTOPH_CASBS=2.0_JPRB
ADLATSTOPH_CA= 180._JPRB/NINT(90._JPRB*NSMAX/159.0_JPRB)
ADLONSTOPH_CA= ADLATSTOPH_CA
IIP= NINT(360._JPRB/ADLONSTOPH_CA)
JJP= NINT(180._JPRB/ADLATSTOPH_CA)

!     -----------------------------------------------------------------
!*    1.4 Set default values common to spectral backscatter and CASBS
!     -----------------------------------------------------------------
NFRSTOPH_SPBS=4     ! must be > 2
NFRSTOPH_SPBS_PAT=NFRSTOPH_SPBS
NSTOCHOPT=2 
SLDISSFAC= 3._JPRB
RSMOOTHSCALE=350000._JPRB

!     -----------------------------------------------------------------
!*    1.5 Set default value for vorticity confinement constant
!     -----------------------------------------------------------------
LVORTCON=.FALSE.
NFRSTOPH_VC= 2
VC_CON= 0.3_JPRB
!     -----------------------------

IU = NULOUT

!*    2. Read Namelist
!     ----------------

CALL POSNAM(NULNAM,'NAMSTOPH')
READ(NULNAM,NAMSTOPH)

IF (LSTOPH_CASBS) NFRSTOPH_SPBS=1  ! variable frequency calls to CASBS not implemented yet

IF (NSTOCHOPT==1) THEN
  ALPHA_SHAL_CONV=2.6E-02_JPRB  
  ALPHA_DEEP_CONV=2.1E-03_JPRB 
ELSEIF (NSTOCHOPT==2) THEN
  ALPHA_SHAL_CONV=2.0_JPRB !not used at the moment (see callpar.F90)
  ALPHA_DEEP_CONV=2.0_JPRB 
ENDIF

!     -----------------------------------------------------------------
!        Set diffusion constant used in SKEB (LSTOPH_SPBS=T)
!          N.B. Consistency with settings for RDIVOR in suhdf_ec may be
!               desirable.
!     -----------------------------------------------------------------
IF( NDGLG > NSMAX+1 ) THEN
  ! cubic grid test, reset diffusion to linear grid
  INSMAXEFF=NDGLG-1
ELSE
  INSMAXEFF=NSMAX
ENDIF

IF (LSPECVIS) THEN
  ZX=0._JPRB
  LSPBS_DISSNUM=.FALSE.
ELSE
  IF(LHDIFFM.AND.LSPBS_DISSNUM_CT) THEN
    ZX=1._JPRB/REAL(NDIFFACT,JPRB)
  ELSE
    ZX=TDT/MAX(0.000001_JPRB,HDIRVOR)
  ENDIF
ENDIF


BIHARM=SLDISSFAC*ZX*(RA**2/(INSMAXEFF*(INSMAXEFF+1)))**2

WRITE(NULOUT,'('' >> Sub SURAND1 ,,'')')


IF (LFORCENL) WRITE(NULOUT,'('' >> add perturbation to model tendencies '')')

!     -----------------------------------------------------------------
!     2.2  If LSTOPH_CASBS=.T. stochastic physics, option:  CASBS
!     -----------------------------------------------------------------
IF (LSTOPH_CASBS) THEN
  WRITE(NULOUT,'('' Logical to control stochastic physics, LSTOPH_CASBS='',L1)') LSTOPH_CASBS  
  WRITE(NULOUT,'('' >> initializing stochastic physics, option:  CASBS ,,'')')

  WRITE(NULOUT,'('' >> Setup for stochastic physics '')')
  WRITE(NULOUT,'('' Variables in common STOPH_MIX'')')
  WRITE(NULOUT,'(&
   & '' Magnitude of stochastic perturbation AMAGSTOPH_CASBS ='',E12.5)')&
   & AMAGSTOPH_CASBS
  WRITE(NULOUT,'('' Box size in lat (deg) ADLATSTOPH_CA ='',E12.5)') ADLATSTOPH_CA
  WRITE(NULOUT,'('' Box size in long (deg) ADLONSTOPH_CA='',E12.5)') ADLONSTOPH_CA

!*    2.2.1  Allocate 
!     ---------------
  ALLOCATE( YDSTOPH%SPSTREAM_FORC(NFLSUR, NSPEC2) ) 

!*    2.2.2  Allocate CASBS only
!     --------------------------
  ALLOCATE(YDSTOPH%RSTOPHCA(NGPTOT))
  YDSTOPH%RSTOPHCA=0.0
  ALLOCATE(YDSTOPH%MCELL(IIP*4,JJP*4))
  ALLOCATE(YDSTOPH%RWGHT(IIP,JJP))

!*   2.2.3  Seed for random number generation
!       Generate a unique number from the date and the ensemble member
!       number (NENSFNB from yomgrb.h) and the processor element PE
!       initialise the cellular automaton by randomly seeding with living cell clumps
!     ----------------------------------------

  CALL SETRAN (NENSFNB,YDML_PHY_STOCH%YR_RANDOM_STREAMS%STOCHPHYS_CABS)

  IF (.NOT. LCUCONV_CA) CALL INITIALIZE_CELLS(YDML_PHY_STOCH%YR_RANDOM_STREAMS,YDSTOPH%MCELL,YDSTOPH%RWGHT,IIP,JJP)

ENDIF

!     -----------------------------------------------------------------
!     2.3  If LSTOPH_SPBS=.T. stochastic physics, option:  SPBS
!     -----------------------------------------------------------------
IF (LSTOPH_SPBS) THEN
  IF (NSMAXSPBS > NSMAX .OR. NSMAXSPBS < 1) THEN
     NSMAXSPBS=NSMAX
  ENDIF

  WRITE(NULOUT,'('' Logical to control new stochastic physics, LSTOPH_SPBS='',L1)') LSTOPH_SPBS
  WRITE(NULOUT,'('' >> initializing stochastic physics, option:  SPBS ,,'')')
  WRITE(NULOUT,'('' >> Setup for stochastic physics SPBS '')')
  WRITE(NULOUT,'('' Variables in module STOPH_MIX'')')
  WRITE(NULOUT,'('' Logical to control vorticity ansatz, LSTOPH_SPBS_VORT='',L1)') LSTOPH_SPBS_VORT
  WRITE(NULOUT,'('' Forcing only wavenumbers 1 to '',I5)')  NSMAXSPBS
  WRITE(NULOUT,'('' Backscatter ratio, RATIO_BACKSCAT '',E12.5)') RATIO_BACKSCAT
  WRITE(NULOUT,'('' Backscatter ratio, RATIO_BACKSCAT_CON2NUM '',E12.5)') RATIO_BACKSCAT_CON2NUM
  WRITE(NULOUT,'('' Exponent for energy spectra, REXPONENT ='',E12.5)') REXPONENT
  WRITE(NULOUT,'('' Decorrelation time of noise, RSPBS_TAU ='',E12.5)') RSPBS_TAU
  WRITE(NULOUT,'('' Variance of noise, RSIGMA2_EPS  ='',E12.5)') RSIGMA2_EPS
  WRITE(NULOUT,'('' Frequency of random generation NFRSTOPH_SPBS ='',I10)')  NFRSTOPH_SPBS
  WRITE(NULOUT,'('' Calling pattern update only every NFRSTOPH_SPBS timestep ='',L1)')  LSTOPH_SPBS_FAST
  WRITE(NULOUT,'('' Bi-harmonic diffusion coefficient (BIHARM) x dt x SLDISSFAC ='',E13.5)')  BIHARM
  WRITE(NULOUT,'('' Force balanced stochastic perturbations ='',L1)')  LSPBSBAL
  WRITE(NULOUT,'('' If true, random numbers drawn from uniform distribution,'',&
   & ''otherwise from a Gaussian distribution='',L1)')  LSTOPH_UNIFORM
  WRITE(NULOUT,'('' Total dissipion rate computed using: '')') 
  WRITE(NULOUT,'('' Convective Mass Flux clipped at '',E12.5)') RFLUX_DET_CLIP
  IF (NSTOCHOPT==1) THEN
    WRITE(NULOUT,'(''   Massflux-scheme '')') 
  ELSEIF (NSTOCHOPT==2.OR.NSTOCHOPT==3) THEN
    WRITE(NULOUT,'(''   Updraft scheme, NSTOCHOPT= '',I1)') NSTOCHOPT 
  ENDIF
  WRITE(NULOUT,'('' LSPBS_DISSNUM='',L1,'' LSPBS_DISSCU='',L1,'' LSPBS_DISSGW='',L1)') LSPBS_DISSNUM, LSPBS_DISSCU, LSPBS_DISSGW
  WRITE(NULOUT,'('' LSPBS_DISSNUM_CT='',L1)') LSPBS_DISSNUM_CT
  WRITE(NULOUT,'('' PBL tapering, LSTOPH_TAPER '',L1)')  LSTOPH_TAPER
  WRITE(NULOUT,'('' Markov chain initialization, LSTOPH_INI '',L1)')  LSTOPH_INI
  WRITE(NULOUT,'('' Gaussian smoothing of dissipation rate, LSTOPH_GAUSS '',L1)')  LSTOPH_GAUSS
  WRITE(NULOUT,'(''    smoothing decorrelation scale, RSMOOTHSCALE '',E12.5)')  RSMOOTHSCALE
  WRITE(NULOUT,'('' Pattern correlation for SPBS, LSTOPH_JBCOR, LSTOPH_UNCORR, LSTOPH_RVP '',3L1)')&
                & LSTOPH_JBCOR, LSTOPH_UNCORR, LSTOPH_RVP 
  IF (LSTOPH_JBCOR) THEN
    WRITE(NULOUT,'(''   Vorticity correlation matrix used for SPBS random perturbation correlation '')')
  ELSEIF (LSTOPH_UNCORR) THEN
    WRITE(NULOUT,'(''   Each model level is perturbed with uncorrelated random variables '')')
  ELSEIF (LSTOPH_RVP) THEN
    WRITE(NULOUT,'(''   Random profiles are used for vertical correlation in SPBS '')')
    WRITE(NULOUT,'(''      Parameters used for RVP: '')')
    WRITE(NULOUT,'(''         RVP_MULMIN= '',E12.5)') RVP_MULMIN
    WRITE(NULOUT,'(''         RVP_MULMAX= '',E12.5)') RVP_MULMAX
    WRITE(NULOUT,'(''         RVP_MULEXP= '',E12.5)') RVP_MULEXP
    WRITE(NULOUT,'(''         RVP_MULNSMAX= '',E12.5)') RVP_MULNSMAX
    IF (LSTOPH_RVPOLD) THEN
      WRITE(NULOUT,'(''      Old p-dependency is used '')')
      WRITE(NULOUT,'(''         RVP_MUL_A= '',E12.5)') RVP_MUL_A
      WRITE(NULOUT,'(''         RVP_MUL_B= '',E12.5)') RVP_MUL_B
      WRITE(NULOUT,'(''         RVP_MUL_C= '',E12.5)') RVP_MUL_C
      WRITE(NULOUT,'(''         RVP_MUL_D= '',E12.5)') RVP_MUL_D
    ELSE
      WRITE(NULOUT,'(''         RVP_MUL_1= '',E12.5)') RVP_MUL_1
      WRITE(NULOUT,'(''         RVP_MUL_2= '',E12.5)') RVP_MUL_2
    ENDIF
  ELSEIF ((.NOT.LSTOPH_UNCORR).AND.(.NOT.LSTOPH_JBCOR).AND.(.NOT.LSTOPH_RVP)) THEN
    WRITE(NULOUT,'(''   Each model level is perturbed with perfectly correlated random variables '')')
  ENDIF
  IF ((LSTOPH_UNCORR .AND. LSTOPH_JBCOR).OR.(LSTOPH_UNCORR .AND. LSTOPH_RVP).OR.(LSTOPH_JBCOR .AND. LSTOPH_RVP)) THEN
    CALL ABOR1('SURAND1: Pattern correlation settings in surand1.F90 are inconsistent:&
    & Pattern can only be one option of uncorrelated, jb-correlated&
    & or random vertical profile')
  ENDIF
  IF (LSTOPH_VARALPHA) THEN
    WRITE(NULOUT,'(''   Scale-dependent decorrelation time used for SPBS '')')
  ENDIF
  WRITE(NULOUT,'('' ------------------- vorticity confinement -------------------'')')
  WRITE(NULOUT,'('' Logical to control vorticity confinement, LVORTCON='',L1)') LVORTCON
  WRITE(NULOUT,'('' ------------------- temperature backscatter -------------------'')')
  WRITE(NULOUT,'('' Logical to control T-backscatter, LSTOPH_SPBS_T='',L1)') LSTOPH_SPBS_T
  IF (LSTOPH_SPBS_T) THEN
    WRITE(NULOUT,'('' settings for T-backscatter:'')')
    WRITE(NULOUT,'(''   Exponent for T spectra, REXPONENT_T ='',E12.5)') REXPONENT_T
    WRITE(NULOUT,'(''   RATIO_APE2KE ='',E12.5)') RATIO_APE2KE
    IF (LSTOPH_RVP) THEN
      WRITE(NULOUT,'(''   Random profiles are used for vertical correlation in T-backscatter '')')
      WRITE(NULOUT,'(''      Parameters used for RVP: '')')
      WRITE(NULOUT,'(''         RVP_MULMIN_T= '',E12.5)') RVP_MULMIN_T
      WRITE(NULOUT,'(''         RVP_MULMAX_T= '',E12.5)') RVP_MULMAX_T
      WRITE(NULOUT,'(''         RVP_MULEXP_T= '',E12.5)') RVP_MULEXP_T
      WRITE(NULOUT,'(''         RVP_MULNSMAX_T= '',E12.5)') RVP_MULNSMAX_T
      IF (LSTOPH_RVPOLD) THEN
        WRITE(NULOUT,'(''      Old p-dependency is used '')')
        WRITE(NULOUT,'(''         RVP_MUL_A_T= '',E12.5)') RVP_MUL_A_T
        WRITE(NULOUT,'(''         RVP_MUL_B_T= '',E12.5)') RVP_MUL_B_T
        WRITE(NULOUT,'(''         RVP_MUL_C_T= '',E12.5)') RVP_MUL_C_T
        WRITE(NULOUT,'(''         RVP_MUL_D_T= '',E12.5)') RVP_MUL_D_T
      ELSE
        WRITE(NULOUT,'(''         RVP_MUL_1_T= '',E12.5)') RVP_MUL_1_T
        WRITE(NULOUT,'(''         RVP_MUL_2_T= '',E12.5)') RVP_MUL_2_T
      ENDIF
    ENDIF
  ENDIF
!*    2.3.1  Allocate SPBS only
!     -------------------------
  ALLOCATE( YDSTOPH%SPSTREAM(NFLSUR, NSPEC2) )
  ALLOCATE( YDSTOPH%SPSTREAM_FORC(NFLSUR, NSPEC2) )   !  and CASBS
  ALLOCATE( YDSTOPH%SPVELPOT(NFLSUR, NSPEC2) ) 
  ALLOCATE( YDSTOPH%SPVELPOT_FORC(NFLSUR, NSPEC2) )
  ALLOCATE( YDSTOPH%SPG_AMP(0:NSMAX) )
  ALLOCATE(YDSTOPH%RSTOPHCA(NGPTOT))
  ALLOCATE( YDSTOPH%ONEMINALPHA_NFRSPBS(0:NSMAX) )
  ALLOCATE( YDSTOPH%ALPHA_STO(0:NSMAX) )

  YDSTOPH%SPSTREAM=  0._JPRB 
  YDSTOPH%SPSTREAM_FORC=  0._JPRB
  YDSTOPH%SPVELPOT=  0._JPRB 
  YDSTOPH%SPVELPOT_FORC=  0._JPRB
  YDSTOPH%ALPHA_STO=0._JPRB
  YDSTOPH%ONEMINALPHA_NFRSPBS=0._JPRB

  IF (LSTOPH_SPBS_T) THEN
    ALLOCATE( YDSTOPH%SPTEMP(NFLSUR, NSPEC2) )
    ALLOCATE( YDSTOPH%SPTEMP_FORC(NFLSUR, NSPEC2) )
    ALLOCATE( YDSTOPH%SPG_AMP_T(0:NSMAX) )
    ALLOCATE( YDSTOPH%ONEMINALPHA_NFRSPBS_T(0:NSMAX) )
    ALLOCATE( YDSTOPH%ALPHA_STO_T(0:NSMAX) )

    YDSTOPH%SPTEMP=0._JPRB 
    YDSTOPH%SPTEMP_FORC=0._JPRB
    YDSTOPH%ALPHA_STO_T=0._JPRB
    YDSTOPH%ONEMINALPHA_NFRSPBS_T=0._JPRB
  ENDIF


  IF (LSTOPH_SPBS_FAST) THEN
    NFRSTOPH_SPBS_PAT=NFRSTOPH_SPBS
  ELSE
    NFRSTOPH_SPBS_PAT=1
  ENDIF



  DO IN=0, NSMAXSPBS
!    ALPHA_STO(IN)   =  (1._JPRB- EXP(-TDT/RSPBS_TAU))*TANH(((IN+1)/40._JPRB))/TANH(1._JPRB)  
    IF (LSTOPH_VARALPHA) THEN  
      YDSTOPH%ALPHA_STO(IN)   = (1+IN)*TDT/(30*RDAY)  !  4 hrs < tau < 30 days
    ELSE
      YDSTOPH%ALPHA_STO(IN)   = 1._JPRB- EXP(-TDT/RSPBS_TAU) 
    ENDIF

   YDSTOPH%ONEMINALPHA_NFRSPBS(IN)=(1-YDSTOPH%ALPHA_STO(IN))**NFRSTOPH_SPBS_PAT
    !WRITE(NULOUT,*) 'n=',IN,'  alpha_sto(n)=', ALPHA_STO(IN)

    IF (LSTOPH_SPBS_T) THEN
      IF (LSTOPH_VARALPHA) THEN  
        YDSTOPH%ALPHA_STO_T(IN)   = (1+IN)*TDT/(30*RDAY)  !  4 hrs < tau < 30 days
      ELSE
        YDSTOPH%ALPHA_STO_T(IN)   = 1._JPRB- EXP(-TDT/RSPBS_TAU) 
      ENDIF
      YDSTOPH%ONEMINALPHA_NFRSPBS_T(IN)=(1-YDSTOPH%ALPHA_STO_T(IN))**NFRSTOPH_SPBS_PAT
    ENDIF

  ENDDO


!*    2.3.2.  More setup for stochastic backscatter: Set up forcing spectrum
!     ---------------------------------------------------------------------------
  ZCHI    =  0._JPRB
 
  DO IN=1,NSMAXSPBS
    ZCHI(IN)=(IN)**REXPONENT
  ENDDO
  
  ZGAMMAN  =  0._JPRB
  DO IN=0,NSMAXSPBS
    ZGAMMAN= ZGAMMAN + IN*(IN+1)*(2*IN+1)*ZCHI(IN)**2/YDSTOPH%ALPHA_STO(IN)
  ENDDO

  ZCONSTF0=RA*SQRT( 2._JPRB/(TDT*RSIGMA2_EPS*ZGAMMAN) ) 


  YDSTOPH%SPG_AMP= 0._JPRB
  DO JN=1,NSMAXSPBS
      ZSIGMA_NFRSPBS=1.0
      DO IN=1,NFRSTOPH_SPBS_PAT-1
        ZSIGMA_NFRSPBS=ZSIGMA_NFRSPBS+(1-YDSTOPH%ALPHA_STO(JN))**(2*IN)
      ENDDO
      ZSIGMA_NFRSPBS=SQRT(ZSIGMA_NFRSPBS)
      YDSTOPH%SPG_AMP(JN) = ZSIGMA_NFRSPBS*SQRT(YDSTOPH%ALPHA_STO(JN))*ZCONSTF0*ZCHI(JN)
  ENDDO

  IF (LSTOPH_SPBS_T) THEN
    ZCHI    =  0._JPRB
 
    DO IN=1,NSMAXSPBS
      ZCHI(IN)=(IN)**REXPONENT_T
    ENDDO
  
    ZGAMMAN  =  0._JPRB
    DO IN=0,NSMAXSPBS
      ZGAMMAN= ZGAMMAN + IN*(IN+1)*(2*IN+1)*ZCHI(IN)**2/YDSTOPH%ALPHA_STO_T(IN)
    ENDDO

    ZCONSTF0=RA*SQRT( 2._JPRB/(TDT*RSIGMA2_EPS*ZGAMMAN) ) 


    YDSTOPH%SPG_AMP_T = 0._JPRB
    DO JN=1,NSMAXSPBS
      ZSIGMA_NFRSPBS=1.0
      DO IN=1,NFRSTOPH_SPBS_PAT-1
        ZSIGMA_NFRSPBS=ZSIGMA_NFRSPBS+(1-YDSTOPH%ALPHA_STO_T(JN))**(2*IN)
      ENDDO
      ZSIGMA_NFRSPBS=SQRT(ZSIGMA_NFRSPBS)
      YDSTOPH%SPG_AMP_T(JN) = ZSIGMA_NFRSPBS*SQRT(YDSTOPH%ALPHA_STO_T(JN))*ZCONSTF0*ZCHI(JN)
    ENDDO
  ENDIF

!*    2.3.3   Setup index in a global spectral array for random numbers (m, n=m)
!        For random number arrays wavnumbers m must be in increasing order
!        ([m=0,n=0], [0,1], [0,2],...,[1,1], [1,2],...,[2,2], [2,3],..,[NSMAX,NSMAX])
!        to insure reproducibility of random numbers for different numbers of processors

  ALLOCATE(YDSTOPH%NIMRAN(0:NSMAX))

  IN=0
  DO JN=0,NSMAX
      YDSTOPH%NIMRAN(JN)=IN+1
      IN=IN+(NSMAX-JN+1)*2
  ENDDO



!*   2.3.4  Seed for random number generation
!       Generate a unique number from the date and the ensemble member
!       number (NENSFNB from yomgrb.h)
!     ----------------------------------------

  CALL SETRAN (NENSFNB,YDML_PHY_STOCH%YR_RANDOM_STREAMS%STOCHPHYS_SPBS) 
  IF (LSTOPH_SPBS_T)  CALL SETRAN ((NENSFNB+24)*420,YDML_PHY_STOCH%YR_RANDOM_STREAMS%STOPH_SPBS_T) 

!*    2.3.5  For random vertical profiles of SPBS perturbations: 
!       - Setup of parameters
!       - Generate a unique random number seed from the date, the ensemble member
!       number (NENSFNB from yomgrb.h) and the wavenumber.
!     ----------------------------------------

  IF (LSTOPH_RVP) THEN
    ALLOCATE(YDSTOPH%RVP_MULFACT(0:NSMAX))
    YDSTOPH%RVP_MULFACT=RVP_MULMAX
    ZRVP_MULFACT=(RVP_MULMAX-RVP_MULMIN)/((1.*RVP_MULNSMAX)**RVP_MULEXP)
    IRVP_MULNSMAX=MIN(NSMAX,INT(RVP_MULNSMAX))
!DEC$ IVDEP
    DO JN=0,IRVP_MULNSMAX
      YDSTOPH%RVP_MULFACT(JN)=ZRVP_MULFACT*(1.*JN)**RVP_MULEXP+RVP_MULMIN
      !WRITE(NULOUT,*) 'n=',JN,'  RVP_MULFACT(n)=', RVP_MULFACT(JN)
    ENDDO

    ALLOCATE(YDSTOPH%RVP_MUL(NFLEVG))

    IF(LSTOPH_RVPOLD) THEN
!DEC$ IVDEP
      DO JN=1,NFLEVG
          !calculate levelpressure
          ZP=YDVAB%VAF(JN)+YDVAB%VBF(JN)*VP00
          !WRITE (NULOUT,*) '  P(',JN,')=',ZP
          !pressure depentent random number scaling
          YDSTOPH%RVP_MUL(JN)=RVP_MUL_A+RVP_MUL_B*ZP+RVP_MUL_C*ZP**2+RVP_MUL_D*ZP**3
      ENDDO
    ELSE
      DO JN=2,NFLEVG
          !calculate levelpressure
          ZP=((YDVAB%VAF(JN)+YDVAB%VBF(JN)*VP00)*&
        & (YDVAB%VAF(JN-1)+YDVAB%VBF(JN-1)*VP00))**0.5_JPRB 
!       middle (in terms of log(p/ps)) between full levels

          !pressure depentent random number scaling
          YDSTOPH%RVP_MUL(JN)=SQRT((RVP_MUL_1*(-LOG(ZP/VP00))**RVP_MUL_2)*&
            &  LOG((YDVAB%VAF(JN)+YDVAB%VBF(JN)*VP00)/(YDVAB%VAF(JN-1)+YDVAB%VBF(JN-1)*VP00)))
          !WRITE(NULOUT,*) 'lev=',JN,'  RVP_MUL(lev)=',RVP_MUL(JN),'  VP00=',VP00
          YDSTOPH%RVP_MUL(JN)=SIGTOFACT(YDSTOPH%RVP_MUL(JN))
          !WRITE(NULOUT,*) 'lev=',JN,'  RVP_MUL(lev)=',RVP_MUL(JN)
      ENDDO
    ENDIF

    CALL SETRAN (NENSFNB+100,YDML_PHY_STOCH%YR_RANDOM_STREAMS%STOCHPHYS_RVP) 

    IF (LSTOPH_SPBS_T) THEN
      YDSTOPH%RVP_MULFACT_T=(RVP_MULMAX_T-RVP_MULMIN_T)/((1.*RVP_MULNSMAX_T)**RVP_MULEXP_T)
      ALLOCATE(YDSTOPH%RVP_MULFACT_T(0:NSMAX))
      YDSTOPH%RVP_MULFACT_T=RVP_MULMAX_T
      ZRVP_MULFACT=(RVP_MULMAX_T-RVP_MULMIN_T)/((1.*RVP_MULNSMAX_T)**RVP_MULEXP_T)
      IRVP_MULNSMAX=MIN(NSMAX,INT(RVP_MULNSMAX_T))
!DEC$ IVDEP
      DO JN=0,IRVP_MULNSMAX
        YDSTOPH%RVP_MULFACT_T(JN)=ZRVP_MULFACT*(1.*JN)**RVP_MULEXP_T+RVP_MULMIN_T
        !WRITE(NULOUT,*) 'n=',JN,'  RVP_MULFACT_T(n)=', RVP_MULFACT_T(JN)
      ENDDO

      ALLOCATE(YDSTOPH%RVP_MUL_T(NFLEVG))

      IF(LSTOPH_RVPOLD) THEN
!DEC$ IVDEP
        DO JN=1,NFLEVG
          !calculate levelpressure
          ZP=YDVAB%VAF(JN)+YDVAB%VBF(JN)*VP00
          !pressure depentent random number scaling
          YDSTOPH%RVP_MUL_T(JN)=RVP_MUL_A_T+RVP_MUL_B_T*ZP+RVP_MUL_C_T*ZP**2+RVP_MUL_D_T*ZP**3
        ENDDO
      ELSE
        DO JN=2,NFLEVG
          !calculate levelpressure
          ZP=((YDVAB%VAF(JN)+YDVAB%VBF(JN)*VP00)*&
        & (YDVAB%VAF(JN-1)+YDVAB%VBF(JN-1)*VP00))**0.5_JPRB 
!       middle (in terms of log(p/ps)) between full levels

          !pressure depentent random number scaling
          YDSTOPH%RVP_MUL_T(JN)=SQRT((RVP_MUL_1_T*(-LOG(ZP/VP00))**RVP_MUL_2_T)*&
            & LOG((YDVAB%VAF(JN)+YDVAB%VBF(JN)*VP00)/(YDVAB%VAF(JN-1)+YDVAB%VBF(JN-1)*VP00)))
          !WRITE(NULOUT,*) 'lev=',JN,'  RVP_MUL_T(lev)=',RVP_MUL_T(JN),'  VP00=',VP00
          YDSTOPH%RVP_MUL_T(JN)=SIGTOFACT(YDSTOPH%RVP_MUL_T(JN))
          !WRITE(NULOUT,*) 'lev=',JN,'  RVP_MUL_T(lev)=',RVP_MUL_T(JN)
        ENDDO
      ENDIF

      CALL SETRAN (NENSFNB+420,YDML_PHY_STOCH%YR_RANDOM_STREAMS%STOPH_RVP_T) 
    ENDIF
  ENDIF


!*    2.3.6  For JB type vertical correlations of SPBS perturbations read a vorticity correlation matrix
!            Calculate sqrt of matrix, distribute result and use in surand2.F90
!             to calculate vertically correlated patterns
!     ----------------------------------------

  IF (LSTOPH_JBCOR) THEN

    !       Open the correlation matrix file and skip past the header

    IOMASTER=1
    ITAG = MTAGFCE+1

    IF (ALLOCATED(YDSTOPH%SQRTCORR)) DEALLOCATE(YDSTOPH%SQRTCORR)
    ALLOCATE(YDSTOPH%SQRTCORR(NFLEVG,NFLEVG))

    IF (MYPROC==IOMASTER) THEN
      CLFILE='vorcorrmatrix.cv'
      WRITE (NULOUT,*) '  Opening input file ',CLFILE
      IULTMP = RESERVE_LUN()
      OPEN  (IULTMP,FILE=CLFILE,FORM='unformatted')
      READ  (IULTMP) CLID
      IF (CLID(1:7)/='VORCORR') CALL ABOR1('Bad Vorticity correlation matrix for SPBS')
      READ  (IULTMP)
      READ  (IULTMP)

      READ (IULTMP) INLEVS
      READ (IULTMP)

      IF (INLEVS /= NFLEVG) THEN
        WRITE (NULOUT,*) 'INLEVS=',INLEVS,'   NFLEVG=',NFLEVG
        CALL ABOR1('Vorticity correlation matrix size not NFLEVGxNFLEVG')
      ENDIF
      DO JY=1,NFLEVG
        READ (IULTMP) (ZCORRMAT(JX,JY),JX=1,NFLEVG)
      ENDDO

      CLOSE(IULTMP)
      CALL FREE_LUN(IULTMP)

      WRITE(NULOUT,'(/35X,'' Vorticity correlations (*100) '')')
      WRITE(NULOUT,832) (JLEV, JLEV=1,NFLEVG)
      DO J1=1,NFLEVG
        WRITE(NULOUT,831) (NINT(100._JPRB*ZCORRMAT(J2,J1)),J2=1,NFLEVG)
      ENDDO
      831 FORMAT(1X,100I3)
      832 FORMAT(/1X,100I3)

      IF (JPRB == N_DEFAULT_REAL_KIND) THEN
        CALL SPOTRF ('U',NFLEVG,ZCORRMAT,NFLEVG,INFO)
      ELSEIF (JPRB == N_DOUBLE_KIND) THEN
        CALL DPOTRF ('U',NFLEVG,ZCORRMAT,NFLEVG,INFO)
      ELSE
       CALL ABOR1 ('REAL(KIND=JPRB) is neither default real nor double precision')
      ENDIF

      IF (INFO > 0) THEN
        WRITE (NULERR,*)'surand1: error computing Cholesky decomposition'
        WRITE (NULERR,*) 'SPOTRF/DPOTRF returns info=',INFO
        WRITE (NULERR,*) 'The matrix is is not positive definite'
        CALL ABOR1('surand1: error computing Cholesky decomposition')
      ELSEIF (INFO < 0) THEN
        WRITE (NULERR,*)'surand1: error computing Cholesky decomposition'
        WRITE (NULERR,*) 'SPOTRF/DPOTRF returns info=',INFO
        WRITE (NULERR,*) 'Illegal value found in matrix'
        CALL ABOR1('surand1: error computing Cholesky decomposition')
      ENDIF

      YDSTOPH%SQRTCORR(:,:)=ZCORRMAT(:,:)


      DO J1=2,NFLEVG
        DO J2=1,J1-1
          YDSTOPH%SQRTCORR(J1,J2)=0.
        ENDDO
      ENDDO


      DO J2=1,NFLEVG
        DO J1=1,NFLEVG
          ZSQRTCORR2(J1,J2)=YDSTOPH%SQRTCORR(J1,J2)**2
        ENDDO
      ENDDO

      IF (LLP)  THEN
        WRITE(NULOUT,'(/35X,'' SQRT Vorticity correlations (*100) '')')
        WRITE(NULOUT,832) (JLEV, JLEV=1,NFLEVG)
        DO J1=1,NFLEVG
          WRITE(NULOUT,831) (NINT(100._JPRB*YDSTOPH%SQRTCORR(J2,J1)),J2=1,NFLEVG)
        ENDDO

        WRITE(NULOUT,'(/35X,'' SQRT Vorticity variance estimate for all levels '')')
        WRITE(NULOUT,833) (JLEV, JLEV=1,NFLEVG)
        WRITE(NULOUT,834) (SUM(ZSQRTCORR2(:,J2)),J2=1,NFLEVG)
        833 FORMAT(1X,100(I4,2X))
        834 FORMAT(/1X,100(F5.3,1X))
      ENDIF

    ENDIF

    ! Distribute the sqrt of the correlation matrix to all PEs

    IF (NPROC>1) THEN
      CALL MPL_BROADCAST (YDSTOPH%SQRTCORR,ITAG,IOMASTER,CDSTRING='SURAND1:')
    ENDIF

  ENDIF

!*    2.3.7  Setup for PBL tapering.
!     ----------------------------------------

  IF (LSTOPH_TAPER) THEN
    TAPER_SIGMATOP=0.87_JPRB
    TAPER_SIGMABOT=0.97_JPRB

    ZDSIGMA = 1._JPRB/(TAPER_SIGMABOT-TAPER_SIGMATOP)**3

   TAPER3 =  2._JPRB                      *ZDSIGMA
   TAPER2 = -3._JPRB*(TAPER_SIGMATOP+TAPER_SIGMABOT)*ZDSIGMA
   TAPER1 =  6._JPRB* TAPER_SIGMATOP*TAPER_SIGMABOT *ZDSIGMA
   TAPER0 =  TAPER_SIGMABOT**2 *&
       &  (TAPER_SIGMABOT - 3._JPRB*TAPER_SIGMATOP) *ZDSIGMA

   ALLOCATE(YDSTOPH%TAPER_FACT(NFLEVG))

    DO JN=1,NFLEVG
       !calculate levelsigma
       ZSIGMA=YDVAB%VAF(JN)/VP00+YDVAB%VBF(JN)
       IF (ZSIGMA > TAPER_SIGMABOT) THEN
          !
          !   below transition layer
          !
              YDSTOPH%TAPER_FACT(JN) =0._JPRB
       ELSE
          IF (ZSIGMA > TAPER_SIGMATOP) THEN
             !
             !   in transition layer
             !
             YDSTOPH%TAPER_FACT(JN)=((TAPER3*ZSIGMA+TAPER2)*ZSIGMA+TAPER1)*ZSIGMA+TAPER0
          ELSE
             YDSTOPH%TAPER_FACT(JN) =1._JPRB
          ENDIF
       ENDIF
    ENDDO
  ENDIF

ENDIF

!--------------------------------------------------------------------------

IF (LSTOPH_SPBS.OR.LSTOPH_CASBS.OR.LVORTCON) THEN
!     -----------------------------------------------------------------
!*    2.x.1  Allocate  fields common to SPBS and CASBS
!     -----------------------------------------------------------------
  WRITE(NULOUT,'('' Entrainment cloud fraction for deep convection, ALPHA_DEEP_CONV  ='',E12.5)')&
    &  ALPHA_DEEP_CONV
  WRITE(NULOUT,'('' Entrainment cloud fraction for shallow convection, ALPHA_DEEP_CONV  ='',E12.5)')&
    &  ALPHA_SHAL_CONV
  WRITE(NULOUT,'('' Factor for numerical dissipation from semi-Lagrangian scheme, SLDISSFAC  ='',E12.5)')&
    &  SLDISSFAC

 

  ALLOCATE(YDSTOPH%GPVORTGRAD(NPROMA,3*NFLEVG,NGPBLKS)) 
  ALLOCATE( YDSTOPH%GPTOTDISS(NPROMA,NFLEVG,NGPBLKS) )
  ALLOCATE( YDSTOPH%GPTOTDISS_SMOOTH(NPROMA,NFLEVG,NGPBLKS) )
  ALLOCATE( YDSTOPH%GPSTREAM(NPROMA,NFLEVG,NGPBLKS) )
  ALLOCATE( YDSTOPH%GPVELPOT(NPROMA,NFLEVG,NGPBLKS) )
  ALLOCATE( YDSTOPH%RSMOOTH(0:NSMAX) )
  

  YDSTOPH%GPVORTGRAD=  0._JPRB
  YDSTOPH%GPTOTDISS=  0._JPRB
  YDSTOPH%GPTOTDISS_SMOOTH = 0._JPRB
  YDSTOPH%GPSTREAM =  0._JPRB
  YDSTOPH%GPTOTDISS_SMOOTH = 0._JPRB
  YDSTOPH%GPVELPOT = 0._JPRB

  IF (LSTOPH_SPBS_T) THEN
    ALLOCATE(YDSTOPH%GPTEMP(NPROMA,NFLEVG,NGPBLKS)) ! ML: tempor. fix to avoid segmentation fault in backscatter_layer when running on cray
    YDSTOPH%GPTEMP =  0._JPRB
  ENDIF

  ALLOCATE(YDSTOPH%SPDP(NFLEVG))
  DO  JN=1,NFLEVG
    YDSTOPH%SPDP(JN)=YDVAB%VDELA(JN)+YDVAB%VDELB(JN)*101325._JPRB
  ENDDO

!*    2.x.2  More setup: Spectral filter for dissipation rates 
!     ---------------------------------------------------------


  IF (LSTOPH_GAUSS) THEN
    ZSIGMA_SMOOTH=0.25_JPRB*(RSMOOTHSCALE/RA)**2
    DO IN=0,NSMAX
      YDSTOPH%RSMOOTH(IN)= EXP(-ZSIGMA_SMOOTH*IN*(IN+1))
    ENDDO
  ELSE
    IFILT= 30
    ICUTOFF= IFILT + 10
    ZAFILT= 3._JPRB/(ICUTOFF*(ICUTOFF+1) - IFILT*(IFILT+1))
    WRITE(NULOUT,*) 'IFILT=',IFILT,'   ICUTOFF=',ICUTOFF

    DO IN=0,NSMAX
      IF(IN>=IFILT.AND.IN<=ICUTOFF) THEN
        YDSTOPH%RSMOOTH(IN)= EXP(-ZAFILT*(IN*(IN+1)-IFILT*(IFILT+1)))
      ELSEIF (IN<IFILT) THEN
        YDSTOPH%RSMOOTH(IN)= 1._JPRB
      ELSE
        YDSTOPH%RSMOOTH(IN)= 0._JPRB
      ENDIF
    ENDDO
  ENDIF
ENDIF

!SPIN UP
ISTEP=0_JPIM
IF (LSTOPH_SPBS) THEN
   CALL SURAND2(YDGEOMETRY,YDECUCONVCA,YDML_PHY_STOCH,ISTEP)
ENDIF
IF (LSTOPH_CASBS .AND. (.NOT. LCUCONV_CA)) THEN
  DO I=1,2*ISPINUP
   CALL SURAND2(YDGEOMETRY,YDECUCONVCA,YDML_PHY_STOCH,ISTEP)
  ENDDO
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SURAND1',1,ZHOOK_HANDLE)
END SUBROUTINE SURAND1
