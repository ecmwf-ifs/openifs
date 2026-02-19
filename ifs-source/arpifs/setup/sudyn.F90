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
SUBROUTINE SUDYN(YDGEOMETRY,YDMODEL,KULOUT)

!------------------------------------------------------------------------------
!**** *SUDYN*   - Initialize constants and control for the dynamics (MODEL object part).

!     Purpose.
!     --------
!           Initialize YOMDYN, the module that controls the
!           dynamics of the model (MODEL object part). Initialize semi-implicit
!           including computation of vertical structure matrix.
!           Also computes vertical modes .
!**   Interface.
!     ----------
!        *CALL* *SUDYN(KULOUT)

!        Explicit arguments :
!        --------------------
!        YDGEOMETRY,YDMODEL
!        KULOUT : Logical unit for the output

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      Mats Hamrud and Philippe Courtier  *ECMWF*
!      Original : 87-10-15

! Modifications
! -------------
!   K. Yessad (Aug 2009): prune some unused obsolete options.
!   F. Vana  16-Oct-2009: option NSPLTHOI + new SLHD triggering by D2
!   K. Yessad (Nov 2009): prune lpc_old.
!   I. Santos, I. Martinez (Mar 2010): Variable map factor for RTM.
!   K. Yessad (Jan 2011): introduce INTDYN_MOD structures; add RPRES_SVTSM.
!   F. Vana  22-Feb-2011: 3D turbulence, (MF) physical tendencies diff, N[X]LAG=4
!   G.Mozdzynski (Feb 2011): OOPS cleaning, use of derived type TGSGEOM
!   K. Yessad (Nov 2011): various contributions, code reordering, printings reordering.
!   M. Fisher   7-March-2012 Move SIBI out of Jb (to allow late Jb setup)
!   K. Yessad (Jul 2012): calls to SUSTA and SUSTADLR moved in SU0YOMB.
!   F. Vana    15-March-2013: Revisited diffusion setup for LECMWF=T
!   K. Yessad (Sep 2013): Some checkings and printings about YOMDYNA SLHD variables moved in SUDYNA.
!   J. Vivoda (Oct 2013): VFE-NH
!   2013-11, J. Masek: Conversion of ACRANEB2 thermal intermittent
!                      frequency NTHRAYFR from hours to timesteps.
!   F. Vana    13-Feb-2014  SLHD setup for IFS
!   M. Diamantakis (Feb 14) : code for SETTLSTF=T  (modifying SETTLST in upper atmosphere) +
!                             Bermejo Staniforth interpolation limiter with improved conservation
!   K. Yessad (July 2013): Allocation of FTHSOR moved in SUSC2B.
!   T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!   K. Yessad (Oct 2013): allow NESC without ICI-PC scheme.
!   K. Yessad (Oct 2013): use LSETTLSV instead of LSETTLST for vertical SL displacement.
!   K. Yessad (July 2014): some reorganisation in the set-up.
!   M. Diamantakis (Feb 2015): change NITMP defaults
!   F. Vana and M. Diamantakis (Aug 2016): simplification and regularization of LSETTLSVF=t
!   K. Yessad (Dec 2016): Prune obsolete options.
!   K. Yessad (June 2017): Introduce NHQE model.
!   2017-09, J. Masek: Safety for SW/LW intermittency with DFI.
!   F. Vana    21-Nov-2017: Options LSLDP_CURV and LHOISLT
!   K. Yessad (Feb 2018): remove deep-layer formulations.
!   K. Yessad (May 2018): use RATIO_HDI_TOP rather than NSREFDH.
!   K. Yessad (June 2018): Alternate NHEE SI scheme elimination.
!   F. Vana July 2018: RK4 scheme for trajectory research and LSLHDSPONGE.
!   R. El Khatib  03-Sep-2018 new configuration 904 which is a test program for change of resolution of an object FIELDS
!   F. Vana 20-Feb-2019: Vertical quintic interpolation for SL scheme
!   F. Vana 21-May-2019: Model internal parameter NLEV_SPONGE
!   F. Vana  24-Jun-2019  QM/QMH fixers available in TL/AD code.
!   F. Vana 11-Jul-2019: Option LRHS_CURV and better setup for SLHDSPONGE
!   F. Vana 14-Sep-2020: SLAVEPP setup from here and SLHD with 1st order extrapol
!   I. Polichtchouk (Nov 2020): added KLSPONGE, RBEGSPONGE, M0DAMP - sponge settings
!   M. Diamantakis (Nov 2021): Option LSLDP_XYZ, SAVE and control for OOPS
! End Modifications
!-------------------------------------------------------------------------------

USE TYPE_MODEL   , ONLY : MODEL
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB, JPRD
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN       , ONLY : NULNAM
USE YOMCT0       , ONLY : LR2D, LR3D, NCONF, LELAM, N2DINI, LECMWF
USE YOMCST       , ONLY : RD
USE YOMDFI       , ONLY : RTDFI
USE YOMVERT      , ONLY : VP00
USE YOMINI       , ONLY : LDFI
USE INTDYNSL_MOD , ONLY : SUINTDYNSL
#ifdef WITH_MGRIDS
USE MGRIDS_ADVECTION_MODULE, ONLY : MGRIDS_ADVECTION
#endif
!     ------------------------------------------------------------------
IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(INOUT) :: YDGEOMETRY
TYPE(MODEL),TARGET,INTENT(INOUT) :: YDMODEL
INTEGER(KIND=JPIM),INTENT(IN)    :: KULOUT 

!     ------------------------------------------------------------------

LOGICAL :: LL2TLFFX
LOGICAL :: LLSTRSG ! .T. => spherical stretched geometry.
LOGICAL :: LLCOMPUTE_CVGQ
LOGICAL :: LLNH_NOC1

INTEGER(KIND=JPIM) :: JLEV, I_NSREFDHX

REAL(KIND=JPRB) :: SLPHY_P_LOW,SLPHY_P_HIGH ! To be moved to namelist
REAL(KIND=JPRB) :: ZEPS, ZSLHDP1, ZSLHDP3, ZLXY, ZSTPHR, ZTSTEP, ZNOTUSED, ZMASK

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

REAL(KIND=JPRB), POINTER  :: REPS1, REPS2, REPSP1, REPSM1, REPSM2,&
 & HDIRVOR, HDIRDIV, HDIRT, HDIRQ, HDIRO3, HDIRPD, HDIRVD, HDIRSP,&
 & BETADT, XIDT, REFGEO, RW2TLFF, SIPR, SITR, SITRA, SITRUB, SIPRUB, &
 & RRDXTAU, RDAMPVOR, RDAMPDIV, RDAMPT, RDAMPQ, RDAMPO3, RDAMPPD, RDAMPVD,&
 & RDAMPSP, REXPDH, FRANDH, SLEVDH, SLEVDH1, SLEVDH2, SLEVDH3, VNORM,&
 & VMAX1, VMAX2, VESL, RMAX_D3, RCMSLP0, RTEMRB, VETAON, VETAOX, &
 & SLHDA0, SLHDB, SLHDDIV, SLHDRATDDIV, SLHDHOR, REXPDHS, SLEVDHS, SLEVDHS1,&
 & SLEVDHS2, RDAMPDIVS, SDRED, RDAMPVORS, RDAMPVDS, RRFZ1, RRFPLM, RRFTAU,&
 & RPRES_SVTSM, RPROFHDBT, RPROFHDTP, RPROFHDMX, RPROFHDEX, RCLSTRESS,RCLPOLE,&
 & SLHDA0T, SLHDBT, SLHDD00, SLHDD00T, RPRES_SETTLSVF, RSCALE, RSCALEOFF, RALPHA, &
 & RALPHA_TOP,RATIO_HDI_TOP,WENO_ALPHA_W, WENO_ALPHA_T, WENO_ALPHA_SP, &
 & WENO_ALPHA_SPD, WENO_ALPHA_SVD, KLSPONGE, RBEGSPONGE, SLHD_P_LOW, SLHD_P_HIGH, RINTOPT

LOGICAL,POINTER :: LSLHDHEAT, LSLHDSPONGE, LSETTLSVF, LSETFSTAT, LSETTLSVF_DIFF, &
 & LGPMASCOR, LSPECVIS, LTOP_VOR

INTEGER(KIND=JPIM),POINTER :: NGPMASCOR, NTOP_VOR_TRUNC, NTOP_VOR_BOT, NOPT_SITRA

INTEGER(KIND=JPIM), POINTER :: NVLAG, NITMP, NWLAG, NTLAG, NSPDLAG, NEDER, NLEV_ZALPHA,&
 & NSVDLAG, NSPLTHOI, NSITER, NSREFDH, NCOMP_CVGQ, NITERHELM, NPROFILEHD, NQMHOISLT,&
 & NDIFFACT, NVSEPC, NVSEPL, NFLEVSF

LOGICAL, POINTER :: LADVF, LIMPF, LNEWHD,&
 & LQMW, LQMHW, LQMT, LQMHT, LQMP, LQMHP, LQMPD, LQMHPD, LQMVD, LQMHVD, &
 & LRHDI_LASTITERPC, LSVTSM, LRDISPE_EC, LHDIFFM, LGPSTRESS, LM0DAMP, &
 & LMASDRY, LMASCOR, LBOUND_D3, LRFRICISOTR, LRFRIC, LWENOBC, LSLDP_CURV, LSLDP_RK, &
 & LSLDP_XYZ, LSLDP_SAVE, LRHS_CURV, LNODRYFLX, &
 & LVWENO_W, LVWENO_T, LVWENO_SP, LVWENO_SPD, LVWENO_SVD

#include "namdyn.nam.h"

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "posnam.intfb.h"
#include "sualdyn.intfb.h"
#include "sualdynb.intfb.h"
#include "suehdf.intfb.h"
#include "suedyn.intfb.h"
#include "sueldynb.intfb.h"
#include "suhdf.intfb.h"
#include "suhdf_ec.intfb.h"
#include "suhdf2.intfb.h"
#include "suhdir.intfb.h"
#include "suhdu.intfb.h"
#include "suheg.intfb.h"
#include "sunhsi.intfb.h"
#include "sunhsi_testconv.intfb.h"
#include "sunheesi.intfb.h"
#include "sunhqesi.intfb.h"
#include "surcordi.intfb.h"
#include "surayfric.intfb.h"
#include "susi.intfb.h"
#include "set_slhd_sponge.intfb.h"
#include "susavtend.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUDYN',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, &
 & YDSPGEOM=>YDGEOMETRY%YSPGEOM, YDCSGEOM=>YDGEOMETRY%YRCSGEOM, YDCSGEOM_NB=>YDGEOMETRY%YRCSGEOM_NB, &
 & YDGSGEOM_NB=>YDGEOMETRY%YRGSGEOM_NB,YDCVER =>YDGEOMETRY%YRCVER,&
 & YDLAP=>YDGEOMETRY%YRLAP, YDSTA=>YDGEOMETRY%YRSTA, &
 & YDVAB=>YDGEOMETRY%YRVAB, YDVETA=>YDGEOMETRY%YRVETA, YDVFE=>YDGEOMETRY%YRVFE,  &
 & YDDYN=>YDMODEL%YRML_DYN%YRDYN,YDPHY=>YDMODEL%YRML_PHY_MF%YRPHY,YDLDDH=>YDMODEL%YRML_DIAG%YRLDDH, &
 & YDSLPHY=>YDMODEL%YRML_PHY_G%YRSLPHY, &
 & YDRIP=>YDMODEL%YRML_GCONF%YRRIP,YDEPHY=>YDMODEL%YRML_PHY_EC%YREPHY, &
 & YDARPHY=>YDMODEL%YRML_PHY_MF%YRARPHY,YDEDYN=>YDMODEL%YRML_DYN%YREDYN,YDEWCOU=>YDMODEL%YREWCOU, &
 & YGFL=>YDMODEL%YRML_GCONF%YGFL,YDTLSCAW=>YDMODEL%YRML_DYN%YYTLSCAW,YDTLSCAWH=>YDMODEL%YRML_DYN%YYTLSCAWH, &
 & YDTRSCAW=>YDMODEL%YRML_DYN%YYTRSCAW,YDTRSCAWH=>YDMODEL%YRML_DYN%YYTRSCAWH, &
 & YDTSCO=>YDMODEL%YRML_DYN%YYTSCO,YDTCCO=>YDMODEL%YRML_DYN%YYTCCO,YDDYNA=>YDMODEL%YRML_DYN%YRDYNA)

ASSOCIATE(LMPA=>YDARPHY%LMPA, &
 & NSMAX=>YDDIM%NSMAX, NMSMAX=>YDDIM%NMSMAX, NDGLG=>YDDIM%NDGLG, &
 & NFLEVG=>YDDIMV%NFLEVG, &
 & HDSRDIV=>YDDYN%HDSRDIV, HDSRVD=>YDDYN%HDSRVD, HDTIME_STRHD=>YDDYN%HDTIME_STRHD, &
 & HRDIRO3=>YDDYN%HRDIRO3, HRDIRPD=>YDDYN%HRDIRPD, HRDIRQ=>YDDYN%HRDIRQ, &
 & HRDIRSP=>YDDYN%HRDIRSP, HRDIRT=>YDDYN%HRDIRT, HDSRVOR=>YDDYN%HDSRVOR, &
 & HRDIRDIV=>YDDYN%HRDIRDIV, HRDIRVD=>YDDYN%HRDIRVD, HRDIRVOR=>YDDYN%HRDIRVOR, &
 & HRDSRDIV=>YDDYN%HRDSRDIV, HRDSRVD=>YDDYN%HRDSRVD, HRDSRVOR=>YDDYN%HRDSRVOR, &
 & LFINDVSEP=>YDDYN%LFINDVSEP, LSIDG=>YDDYN%LSIDG,  L2TLFF=>YDDYN%L2TLFF, &
 & LSTRHD=>YDDYN%LSTRHD, GMASSI=>YDDYN%GMASSI, GPMASSI=>YDDYN%GPMASSI, &
 & NCURRENT_ITER=>YDDYN%NCURRENT_ITER, NRUBC=>YDDYN%NRUBC, NSLDIMK=>YDDYN%NSLDIMK, &
 & RBT=>YDDYN%RBT, RBTS2=>YDDYN%RBTS2, RDAMPHDS=>YDDYN%RDAMPHDS, &
 & RRINTOPT=>YDDYN%RRINTOPT, &
 & SIRPRG=>YDDYN%SIRPRG, SIRPRN=>YDDYN%SIRPRN, SITIME=>YDDYN%SITIME, SIVP=>YDDYN%SIVP, &
 & NLEV_SPONGE=>YDDYN%NLEV_SPONGE,LESIDG=>YDEDYN%LESIDG, &
 & LSLHD=>YDDYNA%LSLHD,LSLHD_STATIC=>YDDYNA%LSLHD_STATIC,LSLHD_W=>YDDYNA%LSLHD_W,&
 & LSLHD_SVD=>YDDYNA%LSLHD_SVD,LPC_FULL=>YDDYNA%LPC_FULL,LPC_CHEAP=>YDDYNA%LPC_CHEAP, &
 & LRUBC=>YDDYNA%LRUBC,LHOISLT=>YDDYNA%LHOISLT,LSLVF=>YDDYNA%LSLVF, &
 & SLHDKMIN=>YDDYNA%SLHDKMIN,SLHDKMAX=>YDDYNA%SLHDKMAX,LSLHD_OLD=>YDDYNA%LSLHD_OLD, &
 & SLHDEPSH=>YDDYNA%SLHDEPSH,L3DTURB=>YDDYNA%L3DTURB,LSLHD_T=>YDDYNA%LSLHD_T, &
 & LMIXETTLS=>YDDYNA%LMIXETTLS,LSETTLST=>YDDYNA%LSETTLST,LSETTLSV=>YDDYNA%LSETTLSV, &
 & LSETTLS=>YDDYNA%LSETTLS,LSLTVWENO=>YDDYNA%LSLTVWENO,LRHSVWENO=>YDDYNA%LRHSVWENO, &
 & LNHQE_SIHYD=>YDDYNA%LNHQE_SIHYD,LSI_NHEE=>YDDYNA%LSI_NHEE,SLHDKREF=>YDDYNA%SLHDKREF, &
 & LWCOU=>YDEWCOU%LWCOU, &
 & NGPTOT=>YDGEM%NGPTOT, NSTTYP=>YDGEM%NSTTYP, RNLGINC=>YDGEM%RNLGINC, &
 & RSTRET=>YDGEM%RSTRET, SLHDP=>YDGEM%SLHDP, &
 & YQ=>YGFL%YQ, YCVGQ=>YGFL%YCVGQ, YQ_NL=>YGFL%YQ_NL, &
 & LSDDH=>YDLDDH%LSDDH, LRSLDDH=>YDLDDH%LRSLDDH, &
 & LMPHYS=>YDPHY%LMPHYS, NSORAYFR=>YDPHY%NSORAYFR, NRAY=>YDPHY%NRAY, &
 & NRAUTOEV=>YDPHY%NRAUTOEV, &
 & NTHRAYFR=>YDPHY%NTHRAYFR, LSLPHY=>YDEPHY%LSLPHY,&
 & NSTOP=>YDRIP%NSTOP, TSTEP=>YDRIP%TSTEP, &
 & STPREH=>YDSTA%STPREH,STPRE=>YDSTA%STPRE,STTEM=>YDSTA%STTEM)

! Associate pointers for variables in namelist
REPS1            => YDDYN%REPS1
REPS2            => YDDYN%REPS2
REPSP1           => YDDYN%REPSP1
REPSM1           => YDDYN%REPSM1
REPSM2           => YDDYN%REPSM2
LTOP_VOR         => YDDYN%LTOP_VOR
NTOP_VOR_TRUNC   => YDDYN%NTOP_VOR_TRUNC
NTOP_VOR_BOT     => YDDYN%NTOP_VOR_BOT
HDIRVOR          => YDDYN%HDIRVOR
HDIRDIV          => YDDYN%HDIRDIV
HDIRT            => YDDYN%HDIRT
HDIRQ            => YDDYN%HDIRQ
HDIRO3           => YDDYN%HDIRO3
HDIRPD           => YDDYN%HDIRPD
HDIRVD           => YDDYN%HDIRVD
HDIRSP           => YDDYN%HDIRSP
BETADT           => YDDYN%BETADT
XIDT             => YDDYN%XIDT
REFGEO           => YDDYN%REFGEO
RW2TLFF          => YDDYN%RW2TLFF
RALPHA           => YDDYN%RALPHA
RALPHA_TOP       => YDDYN%RALPHA_TOP
NLEV_ZALPHA      => YDDYN%NLEV_ZALPHA
KLSPONGE         => YDDYN%KLSPONGE
LM0DAMP          => YDDYN%LM0DAMP
RBEGSPONGE       => YDDYN%RBEGSPONGE
NEDER            => YDDYN%NEDER
LWENOBC          => YDDYN%LWENOBC
NQMHOISLT        => YDDYN%NQMHOISLT
LSLDP_RK         => YDDYN%LSLDP_RK
LSLDP_CURV       => YDDYN%LSLDP_CURV
LRHS_CURV        => YDDYN%LRHS_CURV
LSLDP_XYZ        => YDDYN%LSLDP_XYZ
LSLDP_SAVE       => YDDYN%LSLDP_SAVE
SIPR             => YDDYN%SIPR
SITR             => YDDYN%SITR
SITRA            => YDDYN%SITRA
SITRUB           => YDDYN%SITRUB
SIPRUB           => YDDYN%SIPRUB
NOPT_SITRA       => YDDYN%NOPT_SITRA
RRDXTAU          => YDDYN%RRDXTAU
RDAMPVOR         => YDDYN%RDAMPVOR
RDAMPDIV         => YDDYN%RDAMPDIV
RDAMPT           => YDDYN%RDAMPT
RDAMPQ           => YDDYN%RDAMPQ
RDAMPO3          => YDDYN%RDAMPO3
RDAMPPD          => YDDYN%RDAMPPD
RDAMPVD          => YDDYN%RDAMPVD
RDAMPSP          => YDDYN%RDAMPSP
REXPDH           => YDDYN%REXPDH
FRANDH           => YDDYN%FRANDH
SLEVDH           => YDDYN%SLEVDH
SLEVDH1          => YDDYN%SLEVDH1
SLEVDH2          => YDDYN%SLEVDH2
SLEVDH3          => YDDYN%SLEVDH3
VNORM            => YDDYN%VNORM
VMAX1            => YDDYN%VMAX1
VMAX2            => YDDYN%VMAX2
VESL             => YDDYN%VESL
RMAX_D3          => YDDYN%RMAX_D3
LBOUND_D3        => YDDYN%LBOUND_D3
RCMSLP0          => YDDYN%RCMSLP0
RTEMRB           => YDDYN%RTEMRB
VETAON           => YDDYN%VETAON
VETAOX           => YDDYN%VETAOX
SLHDA0           => YDDYN%SLHDA0
SLHDA0T          => YDDYN%SLHDA0T
SLHDB            => YDDYN%SLHDB
SLHDBT           => YDDYN%SLHDBT
SLHDD00          => YDDYN%SLHDD00
SLHDD00T         => YDDYN%SLHDD00T
SLHDDIV          => YDDYN%SLHDDIV
SLHDRATDDIV      => YDDYN%SLHDRATDDIV
SLHDHOR          => YDDYN%SLHDHOR
SLHD_P_LOW       => YDDYN%SLHD_P_LOW
SLHD_P_HIGH      => YDDYN%SLHD_P_HIGH
LSLHDHEAT        => YDDYN%LSLHDHEAT
LSLHDSPONGE      => YDDYN%LSLHDSPONGE
REXPDHS          => YDDYN%REXPDHS
SLEVDHS          => YDDYN%SLEVDHS
SLEVDHS1         => YDDYN%SLEVDHS1
SLEVDHS2         => YDDYN%SLEVDHS2
RDAMPDIVS        => YDDYN%RDAMPDIVS
SDRED            => YDDYN%SDRED
RDAMPVORS        => YDDYN%RDAMPVORS
RDAMPVDS         => YDDYN%RDAMPVDS
LRFRIC           => YDDYN%LRFRIC
LRFRICISOTR      => YDDYN%LRFRICISOTR
RRFZ1            => YDDYN%RRFZ1
RRFPLM           => YDDYN%RRFPLM
RRFTAU           => YDDYN%RRFTAU
RPROFHDBT        => YDDYN%RPROFHDBT
RPROFHDTP        => YDDYN%RPROFHDTP
RPROFHDMX        => YDDYN%RPROFHDMX
RPROFHDEX        => YDDYN%RPROFHDEX
RCLSTRESS        => YDDYN%RCLSTRESS
RCLPOLE          => YDDYN%RCLPOLE
NVLAG            => YDDYN%NVLAG
NITMP            => YDDYN%NITMP
NWLAG            => YDDYN%NWLAG
NTLAG            => YDDYN%NTLAG
NSPDLAG          => YDDYN%NSPDLAG
NSVDLAG          => YDDYN%NSVDLAG
NSPLTHOI         => YDDYN%NSPLTHOI
NSITER           => YDDYN%NSITER
NSREFDH          => YDDYN%NSREFDH
RATIO_HDI_TOP    => YDDYN%RATIO_HDI_TOP
NCOMP_CVGQ       => YDDYN%NCOMP_CVGQ
NITERHELM        => YDDYN%NITERHELM
NPROFILEHD       => YDDYN%NPROFILEHD
NDIFFACT         => YDDYN%NDIFFACT
NVSEPC           => YDDYN%NVSEPC
NVSEPL           => YDDYN%NVSEPL
NFLEVSF          => YDDYN%NFLEVSF
RSCALE           => YDDYN%RSCALE
RSCALEOFF        => YDDYN%RSCALEOFF
RINTOPT          => YDDYN%RINTOPT
LIMPF            => YDDYN%LIMPF
LADVF            => YDDYN%LADVF
LNEWHD           => YDDYN%LNEWHD
LQMW             => YDDYN%LQMW
LQMHW            => YDDYN%LQMHW
LQMT             => YDDYN%LQMT
LQMHT            => YDDYN%LQMHT
LQMP             => YDDYN%LQMP
LQMHP            => YDDYN%LQMHP
LQMPD            => YDDYN%LQMPD
LQMHPD           => YDDYN%LQMHPD
LQMVD            => YDDYN%LQMVD
LQMHVD           => YDDYN%LQMHVD
LRHDI_LASTITERPC => YDDYN%LRHDI_LASTITERPC
LSVTSM           => YDDYN%LSVTSM
LRDISPE_EC       => YDDYN%LRDISPE_EC
LSPECVIS         => YDDYN%LSPECVIS
LHDIFFM          => YDDYN%LHDIFFM
LGPSTRESS        => YDDYN%LGPSTRESS
LMASDRY          => YDDYN%LMASDRY
LMASCOR          => YDDYN%LMASCOR
LSETTLSVF        => YDDYN%LSETTLSVF
RPRES_SVTSM      => YDDYN%RPRES_SVTSM
RPRES_SETTLSVF   => YDDYN%RPRES_SETTLSVF
LSETFSTAT        => YDDYN%LSETFSTAT
LSETTLSVF_DIFF   => YDDYN%LSETTLSVF_DIFF
LGPMASCOR        => YDDYN%LGPMASCOR
NGPMASCOR        => YDDYN%NGPMASCOR
LNODRYFLX        => YDDYN%LNODRYFLX
LVWENO_W         => YDDYN%LVWENO_W
LVWENO_T         => YDDYN%LVWENO_T
LVWENO_SP        => YDDYN%LVWENO_SP
LVWENO_SPD       => YDDYN%LVWENO_SPD
LVWENO_SVD       => YDDYN%LVWENO_SVD
WENO_ALPHA_W     => YDDYN%WENO_ALPHA_W
WENO_ALPHA_T     => YDDYN%WENO_ALPHA_T
WENO_ALPHA_SP    => YDDYN%WENO_ALPHA_SP
WENO_ALPHA_SPD   => YDDYN%WENO_ALPHA_SPD
WENO_ALPHA_SVD   => YDDYN%WENO_ALPHA_SVD

!     ------------------------------------------------------------------

!*       1.    Set default values.
!              -------------------

ZNOTUSED=-9999._JPRB
IF (LELAM) THEN
  I_NSREFDHX=MAX(NSMAX,NMSMAX)
ELSE
  I_NSREFDHX=NSMAX
ENDIF

!*       1.03  Temporal filter.

IF(LECMWF) THEN
  REPS1=.1_JPRB
  REPS2=.1_JPRB
  REPSP1=.1_JPRB
  REPSM1=.1_JPRB
  REPSM2=.1_JPRB
ELSE
  IF (YDDYNA%LSLAG.AND.YDDYNA%LTWOTL) THEN
    REPS1=0._JPRB
    REPS2=0._JPRB
    REPSP1=0._JPRB
    REPSM1=0._JPRB
    REPSM2=0._JPRB
  ELSE
    REPS1=.05_JPRB
    REPS2=.05_JPRB
    REPSP1=.05_JPRB
    REPSM1=.05_JPRB
    REPSM2=.05_JPRB
  ENDIF
ENDIF

!*       1.04a Spectral horizontal diffusion.

IF(LECMWF) THEN
  LNEWHD=.FALSE.
  REXPDH=4._JPRB
  FRANDH=0._JPRB
  RRDXTAU=57.32_JPRB
  RDAMPDIV= 1._JPRB
  RDAMPVOR= 1._JPRB
  RDAMPT=1._JPRB
  RDAMPQ=1._JPRB
  RDAMPSP=0._JPRB
  RDAMPO3=0._JPRB
  RDAMPVD=1._JPRB
  RDAMPPD=1._JPRB
  IF(.NOT.YDDYNA%LSLAG) THEN
    RRDXTAU=143.3_JPRB
    RDAMPVOR=2.5_JPRB
    RDAMPT  =2.5_JPRB
    RDAMPQ  =2.5_JPRB
  ENDIF
  NPROFILEHD=1
  LRDISPE_EC=.TRUE.
ELSE
  LNEWHD=.FALSE.
  REXPDH=4._JPRB
  FRANDH=0._JPRB
  RRDXTAU=123._JPRB
  RDAMPDIV= 1._JPRB
  RDAMPVOR= 5._JPRB
  RDAMPT=5._JPRB
  RDAMPQ=5._JPRB
  RDAMPSP=0._JPRB
  RDAMPO3=0._JPRB
  IF(YDDYNA%LNHEE) THEN
    RDAMPVD=1._JPRB
    RDAMPPD=5._JPRB
  ELSEIF(YDDYNA%LNHQE) THEN
    RDAMPVD=1._JPRB
    RDAMPPD=0._JPRB
  ELSE
    RDAMPVD=0._JPRB
    RDAMPPD=0._JPRB
  ENDIF
  NPROFILEHD=3
  LRDISPE_EC=.FALSE.
ENDIF

! Security for new HD setting through namelist 
HDIRVOR= ZNOTUSED
HDIRDIV= ZNOTUSED
HDIRT  = ZNOTUSED
HDIRQ  = ZNOTUSED
HDIRO3 = ZNOTUSED
HDIRSP = ZNOTUSED
HDIRPD = ZNOTUSED
HDIRVD = ZNOTUSED

SLEVDH=1._JPRB
SLEVDH1=0.05_JPRB
SLEVDH2=0.005_JPRB
SLEVDH3=100._JPRB/VP00 ! to be tuned later.
NSREFDH=I_NSREFDHX
RATIO_HDI_TOP=1._JPRB

! RPROFHDBT to RPROFHDEX are used only if NPROFILEHD=2
RPROFHDBT=10000._JPRB
RPROFHDTP=100._JPRB
RPROFHDMX=10000._JPRB
RPROFHDEX=1.0_JPRB

!*       1.04b Semi Lagrangian diffusion (SLHD)

! * RDAMP[X]S:
RDAMPDIVS= 0._JPRB
RDAMPVORS= 0._JPRB
RDAMPVDS = 0._JPRB
RDAMPHDS = 1._JPRB
IF(LSLHD_W .OR. LSLHD_SVD) RDAMPHDS= 1.5_JPRB - RNLGINC * 0.5_JPRB
IF(LSLHD_W) THEN
  RDAMPDIVS= 1._JPRB
  RDAMPVORS= 5._JPRB
ENDIF
IF(LSLHD_SVD) RDAMPVDS = 15._JPRB

! * Reset RDAMPQ to zero if SLHD is applied on moisture.
IF(YQ_NL%LSLHD) RDAMPQ= 0._JPRB

! * Other quantities:
SLHDA0  =  0.25_JPRB
SLHDA0T =  0.04_JPRB
SLHDB   =  4.0_JPRB
SLHDBT  =  4.0_JPRB
SLHDD00 = 0.000065_JPRB
SLHDD00T= 0.000081_JPRB
SLHDDIV =  0.0_JPRB
SLHDRATDDIV = 1.0_JPRB
SLHDHOR     = 0._JPRB
SLHD_P_LOW  = 2000._JPRB ! 20 hPa
SLHD_P_HIGH = 1000._JPRB ! 10 hPa
REXPDHS =  6.0_JPRB
SLEVDHS=0.25_JPRB
SLEVDHS1=0.05_JPRB
SLEVDHS2=0.005_JPRB
SDRED=0._JPRB
IF (LECMWF) THEN
  LSLHDHEAT=.TRUE.
  LSLHDSPONGE=.TRUE.
  RDAMPDIVS= 0._JPRB
  RDAMPVORS= 0._JPRB
ELSE
  LSLHDHEAT=.FALSE.
  LSLHDSPONGE=.FALSE.
ENDIF
! Default for SLHD only set if GMV variables are involved
IF ((LSLHD_W.OR.LSLHD_T.OR.LSLHD_SVD).AND.(.NOT.(LSLHD_STATIC.OR.LSLHDSPONGE))) SDRED=1._JPRB

! The following variables are local in SUDYN and are also present in NAMDYN.
ZSLHDP1 = 1.7_JPRB
ZSLHDP3 = 0.6_JPRB

!*       1.04c Set default horizontal diffusion coefficients for ECMWF
!              (If LNEWHD=T is set by namelist, the setup will be
!               relaunched second time. At this level LNEWHD is always F.)
IF (LECMWF) CALL SUHDIR(YDGEOMETRY,YDDYN,YDDYNA)

!*       1.04d Diffusion

IF(LECMWF) THEN
  LHDIFFM=.TRUE.
  LSPECVIS=.FALSE.
  NDIFFACT=6
  LGPSTRESS=.FALSE.
  RCLSTRESS=1.0_JPRB
  RCLPOLE=0.99_JPRB

  NTOP_VOR_TRUNC=255
  NTOP_VOR_BOT=3
  LM0DAMP   =.FALSE.
  RBEGSPONGE =9._JPRB

  IF (NFLEVG < 60) THEN
    LTOP_VOR=.FALSE.
    KLSPONGE=16.0_JPRB
  ELSEIF (NFLEVG < 91) THEN
    LTOP_VOR=.FALSE.
    KLSPONGE=8.0_JPRB
  ELSEIF (NFLEVG < 137) THEN
    LTOP_VOR=.FALSE.
    KLSPONGE=10.0_JPRB
  ELSE
    LTOP_VOR=.TRUE.
    KLSPONGE=14.0_JPRB
  ENDIF

  IF( (NDGLG > (NSMAX+1)).OR.(NFLEVG==1) ) THEN
  ! apply for quadratic and cubic grid and shallow water
    LHDIFFM=.FALSE.
    LSPECVIS=.TRUE.
  ENDIF
ELSE
  ! Not used for the time being at MF
  LSPECVIS=.FALSE.
  LHDIFFM=.FALSE.
  NDIFFACT=6
  LGPSTRESS=.FALSE.
  RCLSTRESS=1.0_JPRB
  RCLPOLE=0.99_JPRB
  LTOP_VOR=.FALSE.
  NTOP_VOR_TRUNC=255
  NTOP_VOR_BOT=3
  LM0DAMP   =.TRUE.
  RBEGSPONGE=1._JPRB

  IF (NFLEVG < 91) THEN
    KLSPONGE=9.0_JPRB
  ELSEIF (NFLEVG < 137) THEN
    KLSPONGE=11.0_JPRB
  ELSE
    KLSPONGE=16.0_JPRB
  ENDIF
ENDIF

!*       1.07  Dynamical quantities:
!         - SI scheme: BETADT,REFGEO,SIPR,SITR,SITRA,NOPT_SITRA,NITERHELM
!         - Coriolis term: LADVF,LIMPF
!         - Advec scheme, SL traj: RW2TLFF,NITMP,RCMSLP0,LSVTSM,LSLDP_CURV,LRHS_CURV,LSLDP_RK,NQMHOISLT
!         - WENO: RALPHA, RALPHA_TOP, NLEV_ZALPHA, NEDER, LWENOBC
!         - Uncentering: XIDT,VESL
!         - Discretisation of eqns in SL scheme: N[X]LAG

REFGEO=100000._JPRB

BETADT=1._JPRB

IF (LECMWF) THEN
  IF (LR2D) THEN
    XIDT=0.0_JPRB
    LADVF=.FALSE.
    RW2TLFF=0.0_JPRB
    LIMPF=.TRUE.
    NITMP=2
  ELSE
    XIDT=0.05_JPRB
    RW2TLFF=1.0_JPRB
    IF(YDDYNA%LSLAG.AND.YDDYNA%LTWOTL) THEN
      LIMPF=.TRUE.
      LADVF=.FALSE.
    ELSE
      LIMPF=.FALSE.
      LADVF=.TRUE.
    ENDIF
    NITMP=2
    IF (LSETTLST.OR.LSETTLSV) THEN
      IF ((NCONF == 1).OR.(NCONF == 302)) NITMP=3
    ENDIF 
  ENDIF
  VESL=0._JPRB
  IF (YDDYNA%LSLAG) THEN
    RCMSLP0=1._JPRB
  ELSE
    RCMSLP0=0._JPRB
  ENDIF
  IF (YDDYNA%LNHDYN) THEN
    IF (YDDYNA%LSLAG.AND.YDDYNA%LTWOTL) THEN
      ! * for models with top in the mesosphere (above 1hPa)
      IF (0.5_JPRB*YDVAB%VAH(1) < 100.0_JPRB.OR.YDCVER%LVERTFE) THEN
        SITR=350._JPRB
      ELSE
        SITR=300._JPRB
      ENDIF
    ELSE
      SITR=300._JPRB
    ENDIF
    SIPR=85000._JPRB
    NOPT_SITRA=0
    SITRA=70._JPRB
    XIDT=0._JPRB
  ELSE
    IF (YDDYNA%LSLAG.AND.YDDYNA%LTWOTL) THEN
      ! * for models with top in the mesosphere (above 1hPa)
      IF (0.5_JPRB*YDVAB%VAH(1) < 100.0_JPRB.OR.YDCVER%LVERTFE) THEN
        SITR=350._JPRB
        SIPR=100000._JPRB
        XIDT=0._JPRB
      ELSE
        SITR=300._JPRB
        SIPR=100000._JPRB
      ENDIF
    ELSE
      SITR=300._JPRB
      SIPR=80000._JPRB
    ENDIF
    ! ky: SITRA is not used in this case.
    NOPT_SITRA=0
    SITRA=SITR
  ENDIF
ELSE
  LIMPF=.FALSE.
  IF (YDDYNA%LSLAG.AND.YDDYNA%LTWOTL.AND.(.NOT.LR2D)) THEN
    ! * SL2TL scheme, 3D model.
    LADVF=.TRUE.
    RW2TLFF=1._JPRB
    NITMP=3
    RCMSLP0=1._JPRB
    XIDT=0.05_JPRB
    VESL=0._JPRB
    SITR=300._JPRB
    NOPT_SITRA=0
    SITRA=50._JPRB
    SIPR=100000._JPRB
  ELSEIF (YDDYNA%LSLAG.AND.YDDYNA%LTWOTL.AND.LR2D) THEN
    ! * SL2TL scheme, 2D model.
    LADVF=.TRUE.
    RW2TLFF=0._JPRB
    NITMP=3
    RCMSLP0=1._JPRB
    XIDT=0._JPRB
    VESL=0.1_JPRB
    SITR=350._JPRB
    ! ky: SITRA is not used in this case.
    NOPT_SITRA=0
    SITRA=SITR
    SIPR=100000._JPRB
  ELSEIF (YDDYNA%LSLAG.AND.(.NOT.YDDYNA%LTWOTL)) THEN
    ! * SL3TL scheme.
    LADVF=.FALSE.
    RW2TLFF=0._JPRB
    NITMP=3
    RCMSLP0=1._JPRB
    XIDT=0._JPRB
    VESL=0.1_JPRB
    SITR=350._JPRB
    NOPT_SITRA=0
    SITRA=100._JPRB
    SIPR=100000._JPRB
  ELSE
    ! * Eulerian scheme.
    LADVF=.FALSE.
    RW2TLFF=0._JPRB
    NITMP=0
    RCMSLP0=0._JPRB
    XIDT=0._JPRB
    VESL=0._JPRB
    SITR=300._JPRB
    NOPT_SITRA=0
    SITRA=100._JPRB
    SIPR=80000._JPRB
  ENDIF
ENDIF

NVLAG=3
NTLAG=3
NWLAG=3
NSPDLAG=3
NSVDLAG=3

NSPLTHOI=0

RALPHA=2._JPRB 
RALPHA_TOP=2._JPRB 
! Repeats what happens in SUHDF_EC to comply with ECMWF sponge
NLEV_ZALPHA=NINT(KLSPONGE)
NEDER=5
LWENOBC=.false.
NQMHOISLT=0

LRHS_CURV=.false.
IF (LECMWF.AND.LR3D) THEN
  LSLDP_CURV=.true.
  LSLDP_XYZ=.true.
  IF ((NCONF == 1).OR.(NCONF == 302)) THEN
    LSLDP_SAVE=.true.
  ELSE
    LSLDP_SAVE=.false.
  ENDIF
ELSE
  LSLDP_CURV=.false.
  LSLDP_XYZ=.false.
  LSLDP_SAVE=.false.
ENDIF

LSLDP_RK=.false.

! * Smoothing of the vertical velocity in the trajectory computation 
!   in pure pressure levels.

LSVTSM=.FALSE.

RPRES_SVTSM=1.01325_JPRB


! Initialise parameters for trajectory scheme when LSETTLSVF=T

LSETTLSVF = LECMWF .AND. YDDYNA%LSLAG .AND. NSTOP > 0 .AND. ANY(SPREAD(NCONF,1,6) == (/1,131,302,401,501,801/))

LSETFSTAT=.FALSE.
RPRES_SETTLSVF=6000.0_JPRB
RSCALE=1.E14_JPRB
RSCALEOFF=1.E-14_JPRB
LSETTLSVF_DIFF=.FALSE. ! LSETTLSVF.AND.LSLHD

IF (YDDYNA%LNHQE .OR. (YDDYNA%LNHEE .AND. YDCVER%LVERTFE)) THEN
  NITERHELM=1
ELSE
  NITERHELM=0
ENDIF

!*       1.08  Semi-Lagrangian interpolators for GMV and GMVS variables.
!              - Quasi-monotonic interpolations: LQM[X], LQMH[X].

IF (YDDYNA%LSLAG.AND.(NCONF == 1 .OR. NCONF == 302 .OR. NCONF == 904)) THEN
  LQMHW = .TRUE.
  LQMHT = .TRUE.
  LQMP  = .TRUE.
ELSEIF ((.NOT.LECMWF).AND.YDDYNA%LSLAG.AND.(NCONF == 131)) THEN
  LQMHW = .TRUE.
  LQMHT = .TRUE.
  LQMP  = .TRUE.
ELSE
  LQMHW = .FALSE.
  LQMHT = .FALSE.
  LQMP  = .FALSE.
ENDIF

LQMW  = .FALSE.
LQMT  = .FALSE.
LQMHP = .FALSE.
LQMPD = .FALSE.
LQMHPD= .FALSE.
LQMVD = .FALSE.
LQMHVD= .FALSE.

LVWENO_W   = .FALSE.
LVWENO_T   = LECMWF.AND.LRHSVWENO
LVWENO_SP  = .FALSE.
LVWENO_SPD = .FALSE.
LVWENO_SVD = .FALSE.
WENO_ALPHA_W   = 2.0_JPRB
WENO_ALPHA_T   = 0.0_JPRB
WENO_ALPHA_SP  = 0.0_JPRB
WENO_ALPHA_SPD = 0.5_JPRB
WENO_ALPHA_SVD = 2.0_JPRB

!*       Defining dyn. range for GFL interpolation
IF (NCONF==401.OR.NCONF==131.OR.NCONF==601.OR.NCONF==801) THEN
  ! Useful when adjoint is involved
  RINTOPT=1._JPRB  ! no need so far with DP
ELSE
  RINTOPT=1._JPRB
ENDIF

!*       1.09  ICI (PC) schemes.

NCURRENT_ITER=0
NSITER=0
LRHDI_LASTITERPC=.TRUE.

!*       1.10  Control the length of the SL trajectory.

IF(NFLEVG >= 91) THEN
  VMAX1=250._JPRB
  VMAX2=400._JPRB
ELSE
  VMAX1=220._JPRB
  VMAX2=320._JPRB
ENDIF

VETAON=1.0_JPRB
VETAOX=1.0_JPRB

!*       1.11  NH model specific options.

NRUBC=0
IF (LRUBC) THEN
  RTEMRB=280._JPRB
  SITRUB=0._JPRB
  SIPRUB=0._JPRB
  NRUBC=1
ENDIF

IF (ABS(TSTEP) > EPSILON(1._JPRB)) THEN
  RMAX_D3=0.5_JPRB/TSTEP
ELSE
  RMAX_D3=HUGE(1._JPRB)
ENDIF
LBOUND_D3=.FALSE.

!*       1.14  Initial values for safe vertical separation for updates
!              (vectorization of adjoint semi-Lagrangian scheme).

NVSEPC=4
NVSEPL=2
LFINDVSEP=.FALSE.

!*       1.15  Computation of moisture convergence CVGQ.

! Initialisation of NCOMP_CVGQ : switch for computation of Moisture Convergence 
! for french deep convection scheme
!   0 ==> YQ must be spectral and derivative are used 
!   1 ==> YQ unspecified and YCVGQ spectral 
!   2 ==> YQ unspecified and YCVGQ grid point, SL computation are performed
NCOMP_CVGQ=0

!*       1.16  Mass corrector.

LMASCOR = .FALSE.
LMASDRY = .FALSE.
LGPMASCOR = .FALSE.
NGPMASCOR = 0
GPMASSI = 0.0_JPRB
GMASSI  = 0.0_JPRB

!*       1.17  Miscellaneous quantities.

IF (LECMWF) THEN
  VNORM = 0.04_JPRB
ELSE
  VNORM = 1._JPRB  
ENDIF

!*       1.18  Rayleigh friction (if active)

!*  Set default behaviour for Rayleigh friction : true if LECMWF and NFLEVG<100, false otherwise
!*  LRFRICISOTR will only be true if specified so in the namdyn namelist
!*  -------------------------------------------------------------------------------------------
LRFRIC = .FALSE. 
IF (LECMWF .AND. NFLEVG<=91) LRFRIC = .TRUE.
LRFRICISOTR=.FALSE.

!*       1.19 Moist Continuity Equation

LNODRYFLX = .FALSE.

!  Defaults still same as 36r1, to be revised later.
!  Suggested values: RRFZ1=68, RRFPLM=100

RRFZ1=61._JPRB
RRFPLM=990._JPRB
RRFTAU=3._JPRB*86400._JPRB

!     ------------------------------------------------------------------

!*       2.    Modify default values.
!              ----------------------

!        2.1   Read namelist.

CALL POSNAM(NULNAM,'NAMDYN')
READ(NULNAM,NAMDYN)

!        2.2   Overwrite namelist with command line arguments

IF (LELAM) THEN
  CALL SUEDYN(YDEDYN,YDGEOMETRY%YREGEO,YDGEOMETRY%YRGEM)
ELSE
  LESIDG=.FALSE.
ENDIF

!     ------------------------------------------------------------------

!*       3.    Reset variables and test.
!              -------------------------
!        3.00  Reset WENO value
! Repeats what happens in SUHDF_EC to comply with ECMWF sponge

NLEV_ZALPHA=NINT(KLSPONGE)
!*       3.01  Bound VETAON and VETAOX between 1.E-6 and 1.

ZEPS=0.000001_JPRB
VETAON=MAX(ZEPS,MIN(1._JPRB,VETAON))
VETAOX=MAX(ZEPS,MIN(1._JPRB,VETAOX))

!*       3.02  Stretched geometry specific features.

ZEPS=100.0_JPRB*TINY(1.0_JPRB)
LLSTRSG=(.NOT.LELAM).AND.(ABS(RSTRET-1.0_JPRB)>ZEPS)
LSTRHD=LLSTRSG.AND.(NCONF<900.OR.NCONF==904)
LSIDG=LLSTRSG.AND.(NCONF<900.OR.NCONF==904).AND.(NSTOP > 0.OR.LDFI)

!*       3.04  Semi-implicit scheme reference profiles.

SIRPRG=RD*SITR
SIRPRN=1.0_JPRB/SIPR

! * Use negative value or zero value for REFGEO?
IF (REFGEO <= 0.0_JPRB) THEN
  CALL ABOR1(' SUDYN : REFGEO MUST BE GREATER THAN 0.')
ENDIF

!*       3.05  Horizontal diffusion scheme.

! * Check NSREFDH, RATIO_HDI_TOP, SLEVDH.. quantities, and reset SLEVDH2 and SLEVDHS2 when relevant:
ZEPS=100.0_JPRB*TINY(1.0_JPRB)
IF ((NSREFDH >= I_NSREFDHX) .AND. (RATIO_HDI_TOP-1._JPRB < ZEPS)) THEN
  ! no additional diffusion near the top
  NSREFDH=I_NSREFDHX
  RATIO_HDI_TOP=1._JPRB
  SLEVDH1=SLEVDH
  SLEVDH2=SLEVDH
  SLEVDHS1=SLEVDHS
  SLEVDHS2=SLEVDHS
ELSEIF (RATIO_HDI_TOP-1._JPRB >= ZEPS) THEN
  ! RATIO_HDI_TOP has been modified in NAMDYN, and additional diffusion is required near the top:
  ! ignore namelist value of NSREFDH and recompute it
  NSREFDH=NINT(REAL(I_NSREFDHX,JPRB)/(RATIO_HDI_TOP**(1._JPRB/REXPDH)))
ENDIF
IF ( (SLEVDH2 > SLEVDH1) .OR. (SLEVDHS2 > SLEVDHS1) ) THEN
  ! * Use (SLEVDH2 > SLEVDH1) or (SLEVDHS2 > SLEVDHS1) (not coded)?
  WRITE(KULOUT,*) 'SLEVDH2 must be <= SLEVDH1.'
  WRITE(KULOUT,*) 'SLEVDHS2 must be <= SLEVDHS1.'
  CALL ABOR1('SUDYN: ABOR1 CALLED')
ENDIF
IF (NPROFILEHD == 4) THEN
  SLEVDH1=SLEVDH
  SLEVDHS1=SLEVDHS
ENDIF

! * Check LRDISPE_EC:
IF (LRDISPE_EC .AND. LELAM) THEN
  ! * Use LRDISPE_EC in LAM models (not coded)?
  WRITE(KULOUT,*) 'If LELAM=T, LRDISPE_EC must be F.'
  CALL ABOR1('SUDYN: ABOR1 CALLED')
ENDIF

! * HDIR[X], HRDIR[X], HDSR[X], HRDSR[X].
! Check for MF setup 
! (It will always complain for LECMWF=true, as there
!  the diffusion is initialized in the old style: LNEWHD=F.)
IF (.NOT. LECMWF) THEN
  IF((HDIRVOR /= ZNOTUSED).OR.(HDIRDIV /= ZNOTUSED) .OR.&
   & (HDIRT   /= ZNOTUSED).OR.(HDIRQ   /= ZNOTUSED) .OR.&
   & (HDIRO3  /= ZNOTUSED).OR.(HDIRSP  /= ZNOTUSED) .OR.&
   & (HDIRVD  /= ZNOTUSED).OR.(HDIRPD  /= ZNOTUSED)) THEN
    WRITE(KULOUT,*)&
     & 'SUDYN WARNING: The spectral diffusion in this case ',&
     & 'is not any more driven by namelist parameters HDIR[x]!!! ',&
     & 'Please remove them from your namelist and use new set-up.'
    CALL ABOR1('SUDYN: ABOR1 CALLED')
  ENDIF
ENDIF

! Setup for the 'modern style' diffusion
! (i.e. in case of ECMWF the setup is launched second time with the correct settings)
IF ((LECMWF .AND. LNEWHD) .OR. (.NOT. LECMWF)) THEN
  CALL SUHDIR(YDGEOMETRY,YDDYN,YDDYNA)
ELSE
  ! Here we do some minimal promotion of the old setup to be a way
  ! compatible with the new one.
  IF(HDIRVOR /= 0._JPRB) THEN
    HRDIRVOR = 1._JPRB/HDIRVOR
  ELSE
    HRDIRVOR = 0._JPRB
  ENDIF
  IF(HDIRDIV /= 0._JPRB) THEN
    HRDIRDIV = 1._JPRB/HDIRDIV
  ELSE
    HRDIRDIV = 0._JPRB
  ENDIF
  IF(HDIRT   /= 0._JPRB) THEN
    HRDIRT   = 1._JPRB/HDIRT  
  ELSE
    HRDIRT   = 0._JPRB
  ENDIF
  IF(HDIRQ   /= 0._JPRB) THEN
    HRDIRQ   = 1._JPRB/HDIRQ  
  ELSE
    HRDIRQ   = 0._JPRB
  ENDIF
  IF(HDIRO3  /= 0._JPRB) THEN
    HRDIRO3  = 1._JPRB/HDIRO3 
  ELSE
    HRDIRO3  = 0._JPRB
  ENDIF
  IF(HDIRSP  /= 0._JPRB) THEN
    HRDIRSP  = 1._JPRB/HDIRSP 
  ELSE
    HRDIRSP  = 0._JPRB
  ENDIF
  IF(HDIRVD  /= 0._JPRB) THEN
    HRDIRVD  = 1._JPRB/HDIRVD 
  ELSE
    HRDIRVD  = 0._JPRB
  ENDIF
  IF(HDIRPD  /= 0._JPRB) THEN
    HRDIRPD  = 1._JPRB/HDIRPD 
  ELSE
    HRDIRPD  = 0._JPRB
  ENDIF
  ! Supporting diffusion is driven fully by HDIR[x] coef. in the IFS
  HDSRDIV=0._JPRB
  HDSRVOR=0._JPRB
  HDSRVD =0._JPRB
  HRDSRDIV=0._JPRB
  HRDSRVOR=0._JPRB
  HRDSRVD =0._JPRB
ENDIF

! * Warning: with ECMWF, LNEWHD=.F., nam. param. RRDXTAU is ignored
IF (LECMWF.AND..NOT.LNEWHD) THEN
  WRITE(KULOUT,*) 'WARNING: LNEWHD=.F. => RRDXTAU IGNORED (IF SET)'
ENDIF

! * Check value of NPROFILEHD
IF (NPROFILEHD < 1 .OR. NPROFILEHD > 4) THEN
  WRITE(KULOUT,*) ' SUDYN: NPROFILEHD must be between 1 and 4!'
  CALL ABOR1(' SUDYN: ABOR1 CALLED')
ENDIF
IF (LRDISPE_EC.AND.NPROFILEHD == 4) THEN
  WRITE(KULOUT,*) ' SUDYN: NPROFILEHD must be between 1 and 3 if LRDISPE_EC=T!'
  CALL ABOR1(' SUDYN: ABOR1 CALLED')
ENDIF
IF (LELAM.AND.NPROFILEHD == 4) THEN
  WRITE(KULOUT,*) ' SUDYN: NPROFILEHD must be between 1 and 3 if LELAM=T!'
  CALL ABOR1(' SUDYN: ABOR1 CALLED')
ENDIF
IF (LECMWF.AND.(.NOT.LNEWHD).AND.(.NOT.LRDISPE_EC)) THEN
  WRITE(KULOUT,*) ' SUDYN: if LECMWF=T and LRDISPE_EC=F, LNEWHD must be T!'
  CALL ABOR1(' SUDYN: ABOR1 CALLED')
ENDIF

! * Reset values of RPROFHDBT and RPROFHDTP between YDVAB%VAH(0) and VP00
RPROFHDBT=MIN(VP00,MAX(RPROFHDBT,YDVAB%VAH(0)))
RPROFHDTP=MIN(VP00,MAX(RPROFHDTP,YDVAB%VAH(0)))

! * Check that RPROFHDBT >= RPROFHDTP
IF (RPROFHDBT < RPROFHDTP) THEN
  WRITE(KULOUT,*) ' SUDYN: RPROFHDBT must be >= RPROFHDTP!'
  CALL ABOR1(' SUDYN: ABOR1 CALLED')
ENDIF

! * Check that RPROFHDEX >= 0
IF (RPROFHDEX < 0.0_JPRB) THEN
  WRITE(KULOUT,*) ' SUDYN: RPROFHDEX must be >= 0!'
  CALL ABOR1(' SUDYN: ABOR1 CALLED')
ENDIF

! * Check that RPROFHDMX >= 1
IF (RPROFHDMX < 1.0_JPRB) THEN
  WRITE(KULOUT,*) ' SUDYN: RPROFHDMX must be >= 1!'
  CALL ABOR1(' SUDYN: ABOR1 CALLED')
ENDIF

! * Check SLHD:
IF (LSLHD) THEN
  ! * Use SDRED out of interval [0;1]?
  IF ((SDRED < 0.0_JPRB ).OR.(SDRED > 1.0_JPRB )) THEN
    WRITE(KULOUT,*)&
     & 'SUDYN WARNING: SDRED should be within interval <0,1> !'  
    CALL ABOR1('SUDYN: ABOR1 CALLED')
  ENDIF
  IF (.NOT.LELAM .AND. LSLHD_SVD .AND. (RDAMPVDS /= 0._JPRB)) THEN
    WRITE(KULOUT,*)&
     & 'SUDYN WARNING!: When LNHDYN and SLHD are TRUE on sphere',' the supporting diffusion of SVD is not coded yet!'
  ENDIF
  IF (LSLHDHEAT .AND. LELAM) THEN
    WRITE(KULOUT,*)&
     & 'SUDYN WARNING: Different triggering for heat and momentum variables',&
     & '   is currently not coded.'
    CALL ABOR1('SUDYN: ABOR1 CALLED')
  ENDIF
  IF (SLHD_P_LOW <= SLHD_P_HIGH) THEN
    WRITE(KULOUT,*)&
     & 'SUDYN WARNING: The pressure of the lowest level where SLHD acts has to be higher ',&
     & '   than the pressure of full SLHD. Please make sure SLHD_P_LOW > SLHD_P_HIGH. '
    CALL ABOR1('SUDYN: ABOR1 CALLED')
  ENDIF
ENDIF

! * Quasi-horizontal triggering in TL/AD model
!   (Note that the LSLHD can be still redefined in SUCTRL_GFLATTR,
!    thus some bizzare settings might pass through this security.)
IF (LSLHD.AND.(SLHDHOR /= 0._JPRB).AND.(NCONF == 131.OR.NCONF == 401&
 &  .OR.NCONF == 501.OR.NCONF == 601.OR.NCONF == 801).AND.&
 &  (.NOT.LSLHD_STATIC)) THEN
  WRITE(KULOUT,*) ' SUDYN: TL/AD versions of SLHD are coded only for SLHDHOR=0.'
  CALL ABOR1(' SUDYN: ABOR1 CALLED')
ENDIF

!*       3.07  Advection, SL scheme, Coriolis term.
!         - Coriolis term: LADVF,LIMPF
!         - Advec scheme, SL traj: RW2TLFF,NITMP,RCMSLP0,LSVTSM
!         - Discretisation of eqns in SL scheme: N[X]LAG
!         - Option NSPLTHOI
!         - RINTOPT/RRINTOPT

IF (LRSLDDH .AND. (NSPLTHOI==0) .AND. .NOT.LECMWF)  THEN
  NSPLTHOI=-1
ENDIF

IF (LR2D) THEN

  ! * N[X]LAG out of bounds?
  IF ((NVLAG /= 2).AND.(NVLAG /= 3).AND.(NVLAG /= -2)) THEN
    WRITE(KULOUT,*) '2D MODEL. NVLAG MUST BE -2 OR 2 OR 3.'
    CALL ABOR1('SUDYN: ABOR1 CALLED')
  ENDIF
  IF ((NVLAG == -2).AND.YDDYNA%LTWOTL) THEN
    WRITE(KULOUT,*) '2D MODEL. 2TLSL NOT CODED FOR NVLAG=-2.'
    CALL ABOR1('SUDYN: ABOR1 CALLED')
  ENDIF
  IF (NWLAG /= 2.AND.NWLAG /= 3) THEN
    WRITE(KULOUT,*) '2D MODEL. NWLAG MUST BE 2 OR 3.'
    CALL ABOR1('SUDYN: ABOR1 CALLED')
  ENDIF
  ! * Use LSETTLS=T when not coded?
  IF (LSETTLS.AND..NOT.(NWLAG == 3.AND.NVLAG == 3.AND.YDDYNA%LTWOTL)) THEN
    WRITE(KULOUT,*) 'LSETTLS ONLY WORKS WITH LTWOTL and NWLAG=NVLAG=3'
    CALL ABOR1('SUDYN: ABOR1 CALLED')
  ENDIF
  ! *  Proper use of LADVF/LIMPF?
  IF (LADVF.AND.LIMPF) THEN
    WRITE(KULOUT,*) 'LADVF and LIMPF CANNOT BOTH BE TRUE.'
    CALL ABOR1('SUDYN: ABOR1 CALLED')
  ENDIF
  IF (LIMPF.AND.(NSTTYP /= 1.OR.LLSTRSG)) THEN
    WRITE(KULOUT,*) 'LIMPF CANNOT BE USED WITH TILTING OR',' STRETCHING.'
    CALL ABOR1('SUDYN: ABOR1 CALLED')
  ENDIF
  ! * Use RCMSLP0 /= 0 when not available?
  IF ((RCMSLP0 /= 0.0_JPRB).AND.(.NOT.YDDYNA%LSLAG)) THEN
    WRITE(KULOUT,*) 'RCMSLP0 /= ZERO CAN BE USED ONLY FOR SL SCHEME.'
    CALL ABOR1('SUDYN: ABOR1 CALLED')
  ENDIF
  IF ((RCMSLP0 /= 0.0_JPRB).AND.YDDYNA%LSLAG.AND.(NVLAG < 0)) THEN
    WRITE(KULOUT,*) 'RCMSLP0 /= ZERO CANNOT BE USED FOR LAGRANGIAN',&
     & ' FORMULATION OF CONTINUITY EQUATION.'  
    CALL ABOR1('SUDYN: ABOR1 CALLED')
  ENDIF
  IF (NSPLTHOI /= 0) THEN
    WRITE(KULOUT,*) 'NSPLTHOI /= 0 IS NOT CODED FOR 2D MODEL'
    CALL ABOR1('SUDYN: ABOR1 CALLED')
  ENDIF

ELSE

  ! * N[X]LAG out of bounds?
  IF (NVLAG < 2.OR.NVLAG > 3) THEN
    WRITE(KULOUT,*) '3D MODEL. NVLAG MUST BE BETWEEN 2 AND 3.'
    CALL ABOR1('SUDYN: ABOR1 CALLED')
  ENDIF
  IF (NWLAG < 2.OR.NWLAG > 4) THEN
    WRITE(KULOUT,*) '3D MODEL. NWLAG MUST BE BETWEEN 2 AND 4.'
    CALL ABOR1('SUDYN: ABOR1 CALLED')
  ENDIF
  IF (NTLAG < 2.OR.NTLAG > 4) THEN
    WRITE(KULOUT,*) '3D MODEL. NTLAG MUST BE BETWEEN 2 AND 4.'
    CALL ABOR1('SUDYN: ABOR1 CALLED')
  ENDIF
  IF (YDDYNA%LNHEE.AND.(NSPDLAG < 2.OR.NSPDLAG > 3)) THEN
    WRITE(KULOUT,*) '3D MODEL. NSPDLAG MUST BE BETWEEN 2 AND 3.'
    CALL ABOR1('SUDYN: ABOR1 CALLED')
  ENDIF
  IF (YDDYNA%LNHDYN.AND.(NSVDLAG < 2.OR.NSVDLAG > 4)) THEN
    WRITE(KULOUT,*) '3D MODEL. NSVDLAG MUST BE BETWEEN 2 AND 4.'
    CALL ABOR1('SUDYN: ABOR1 CALLED')
  ENDIF
  IF ((.NOT.LMPHYS).AND.(NWLAG==4 .OR. NTLAG==4 .OR. NSVDLAG==4)) THEN
    WRITE(KULOUT,*) '3D MODEL. N[x]LAG=4 HAS ONLY MEANING WITH M-F PHYSICS.'
    CALL ABOR1('SUDYN: ABOR1 CALLED')
  ENDIF
  IF (LSETTLS.AND.LPC_CHEAP.AND.(NWLAG/=3 .OR. NTLAG/=3 .OR. NSPDLAG/=3 .OR. NSVDLAG/=3)) THEN
    WRITE(KULOUT,*) 'LSETTLS AND LPC_CHEAP ONLY FOR N[X]LAG=3.'
    CALL ABOR1('SUDYN: ABOR1 CALLED')
  ENDIF
  ! * Check of NSPLTHOI
  IF ((NSPLTHOI /= -1).AND.(NSPLTHOI /= 0).AND.(NSPLTHOI /= 1)) THEN
    WRITE(KULOUT,*) 'WRONG VALUE OF NSPLTHOI: ',NSPLTHOI
    WRITE(KULOUT,*) 'IT CAN ONLY BE SET TO -1,0,1'
    CALL ABOR1('SUDYN: ABOR1 CALLED')
  ENDIF
  IF ((NSPLTHOI /= 0).AND.(NCONF /= 1)) THEN
    WRITE(KULOUT,*) 'NSPLTHOI /= 0 IS ONLY CODED FOR CONF 001.'
    CALL ABOR1('SUDYN: ABOR1 CALLED')
  ENDIF
  IF ((NSPLTHOI /= 0).AND.LECMWF) THEN
    WRITE(KULOUT,*) 'NSPLTHOI /= 0 IS ONLY CODED FOR THE MF PHYSICS DATAFLOW.'
    CALL ABOR1('SUDYN: ABOR1 CALLED')
  ENDIF
  IF ((NSPLTHOI /= 0).AND.LSLPHY) THEN
    WRITE(KULOUT,*) 'NSPLTHOI /= 0 IS ONLY CODED FOR THE LSLPHY=.F. CASE.'
    CALL ABOR1('SUDYN: ABOR1 CALLED')
  ENDIF
  IF ((NSPLTHOI == 1).AND.(SLHDEPSH == 0.0_JPRB).AND.&
    & ((.NOT.LSLHD_OLD).AND.(SLHDKMAX == SLHDKMIN))) THEN
    WRITE(KULOUT,*) 'YOUR NSPLTHOI=1 SETTING WILL JUST CONSUME CPUs BY PERFORMING'
    WRITE(KULOUT,*) 'TWICE THE SAME KIND OF INTERPOLATION: Set either NSPLTHOI=0'
    WRITE(KULOUT,*) 'or activate smoother by better values of SLHDEPSH or SLHDKMAX.'
    CALL ABOR1('SUDYN: ABOR1 CALLED')
  ENDIF
  ! * Use LSETTLS=T when not coded?
  IF (LSETTLS.AND..NOT.(NWLAG >= 3.AND.NTLAG >= 3.AND.NVLAG == 3.AND.YDDYNA%LTWOTL)) THEN  
    WRITE(KULOUT,*) 'LSETTLS ONLY WORKS WITH NWLAG=NTLAG=NVLAG=3,(4) LTWOTL=T'
    CALL ABOR1('SUDYN: ABOR1 CALLED')
  ENDIF
  IF (YDDYNA%LNHDYN.AND.LSETTLS.AND..NOT.(NSPDLAG == 3.AND.NSVDLAG >= 3.AND.YDDYNA%LTWOTL)) THEN  
    WRITE(KULOUT,*) 'LSETTLS ONLY WORKS WITH NSPDLAG=NSVDLAG=3,(4) LTWOTL=T'
    CALL ABOR1('SUDYN: ABOR1 CALLED')
  ENDIF
  ! *  Proper use of LADVF/LIMPF?
  IF (LADVF.AND.LIMPF) THEN
    WRITE(KULOUT,*) 'LADVF and LIMPF CANNOT BOTH BE TRUE.'
    CALL ABOR1('SUDYN: ABOR1 CALLED')
  ENDIF
  IF (LIMPF.AND.(NSTTYP /= 1.OR.LLSTRSG)) THEN
    WRITE(KULOUT,*) 'LIMPF CANNOT BE USED WITH TILTING OR',' STRETCHING.'
    CALL ABOR1('SUDYN: ABOR1 CALLED')
  ENDIF
  ! * Use RCMSLP0 /= 0 when not available?
  IF ((RCMSLP0 /= 0.0_JPRB).AND.(.NOT.YDDYNA%LSLAG)) THEN
    WRITE(KULOUT,*) 'RCMSLP0 /= ZERO CAN BE USED ONLY FOR SL SCHEME.'
    CALL ABOR1('SUDYN: ABOR1 CALLED')
  ENDIF

ENDIF

IF (YDDYNA%LSLAG) THEN
  ! Set up RRINTOPT
  IF (RINTOPT /= 0._JPRB) THEN
    RRINTOPT=1._JPRB/RINTOPT
  ELSE
    RRINTOPT=1._JPRB
  ENDIF
  ! Do the precision check
  IF ((RINTOPT*RRINTOPT < 1._JPRB-ZEPS) .OR. (RINTOPT*RRINTOPT > 1._JPRB+ZEPS)) THEN
     WRITE(KULOUT,*) 'Factor RINTOPT to rise accuracy of SL interpolation will degrade it.'
     CALL ABOR1('SUDYN: ABOR1 CALLED')
  ENDIF
ELSE
  RINTOPT=1._JPRB
  RRINTOPT=1._JPRB
ENDIF

!*       3.08  Check settings of BETADT and XIDT.

! Reset VESL to zero when XIDT is non zero.
IF (LR2D) THEN
  XIDT=0._JPRB
  IF(N2DINI == 11) BETADT=0.0_JPRB
ELSE
  IF (YDDYNA%LSLAG.AND.YDDYNA%LTWOTL) THEN
    IF (XIDT > 0.0_JPRB) THEN
      WRITE(KULOUT,*) ' SUDYN: NON-ZERO XIDT IS REQUESTED,',&
       & ' SO VESL IS RESET TO ZERO.'  
      VESL=0._JPRB
    ENDIF
    IF (XIDT < 0.0_JPRB) THEN
      WRITE(KULOUT,*) ' SUDYN: XIDT<0. NOT ALLOWED,',&
       & ' SO XIDT IS RESET TO ZERO.'  
      XIDT=0._JPRB
    ENDIF
  ELSE
    XIDT=0._JPRB
  ENDIF
ENDIF

!*       3.10  RW2TLFF.

!      Check RW2TLFF and reset if necessary; compute L2TLFF:
!      Refined computation of 2*Omega*Vec*r can be activated only if:
!      - Semi-Lagrangian scheme.
!      - LADVF=T.
!      - NWLAG>=3.

LL2TLFFX=(RW2TLFF > 0._JPRB).AND.LADVF
IF (YDDYNA%LSLAG.AND.(NWLAG >= 3).AND.LL2TLFFX) THEN
  ! * cases where l2tlff=.true. is available and asked for.
  IF (NCONF == 1.OR.NCONF == 302) THEN
    ! - case where l2tlff=.true. and all values of rw2tlff are available.
    L2TLFF=.TRUE.
  ELSE
    ! - case where (l2tlff=.true.; rw2tlff=1.) is available, but not
    !   rw2tlff different from one if non-zero => reset RW2TLFF to one.
    L2TLFF=.TRUE.
    RW2TLFF=1.0_JPRB
  ENDIF
ELSE
  ! * cases where either l2tlff=.true. is not available or not asked for.
  L2TLFF=.FALSE.
  RW2TLFF=0.0_JPRB
ENDIF

!*       3.11  VNORM.

! * Use negative value or zero value for VNORM?
IF (VNORM <= 0.0_JPRB) THEN
  WRITE(UNIT=KULOUT,FMT=&
   & '('' WRONG CHOICE FOR THE SCALING OF VERTICAL EIGENMODES '')')  
  WRITE(KULOUT,*) 'VNORM must be > 0'
  CALL ABOR1('SUDYN : ABOR1 CALLED')
ENDIF

!*       3.12  ICI (PC) schemes.

! * Consistency between NSITER and (LPC_FULL or LRUBC):
IF (NSITER < 0) THEN
  WRITE(KULOUT,*) 'Negative value of NSITER is not allowed'
  CALL ABOR1('SUDYN: ABOR1 CALLED')
ENDIF
IF( NSITER==0.AND.LPC_FULL )THEN
  WRITE(KULOUT,*) 'No iterations NSITER=0 but LPC_FULL is .T.'
  CALL ABOR1('SUDYN: ABOR1 CALLED')
ENDIF
IF( NSITER>0.AND..NOT.LPC_FULL)THEN
  WRITE(KULOUT,*) 'Dynamics with NSITER>0 needs LPC_FULL'          
  CALL ABOR1('SUDYN: ABOR1 CALLED')
ENDIF
IF (LSETFSTAT.AND.LMIXETTLS) THEN
  WRITE(KULOUT,*) 'LMIXETTLS does not work with LSETFSTAT.'
  CALL ABOR1(' SUDYN: ABOR1 CALLED')
ENDIF

!*       3.13  LRUBC.

! * Use LRUBC in spherical geometry or improper options for LRUBC?
IF (LRUBC .AND. NSITER /= 1) THEN
  WRITE(KULOUT,*) ' LRUBC=T requires NSITER=1 '          
  CALL ABOR1('SUDYN: ABOR1 CALLED')
ENDIF

! * Correct NRUBC if necessary.
IF (LRUBC) THEN
  NRUBC=MAX(1,MIN(2,NRUBC))
ELSE
  NRUBC=0
ENDIF

!*       3.15  Moisture convergence.

! * CVGQ for french deep convection scheme:

! ky: definition of LLCOMPUTE_CVGQ must be the same everywhere in the code.
!     (also in CPPHINP).
LLCOMPUTE_CVGQ=LMPHYS.AND.(.NOT.LMPA)

IF (LLCOMPUTE_CVGQ) THEN
  ! * Try to compute an Eulerian CVGQ with spectral YQ-GFL when YQ is
  !   available only in grid-point space?
  IF ((NCOMP_CVGQ==0) .AND. YQ%LGP .AND. (NSTOP>0.OR.LDFI)) THEN
    WRITE(KULOUT,*)&
     & 'SUDYN WARNING: when NCOMP_CVGQ = 0 YQ must be spectral'
    CALL ABOR1('SUDYN: NCOMP_CVGQ inconsistency ABOR1 CALLED')
  ENDIF
  ! * Try to compute an Eulerian CVGQ with spectral YCVGQ-GFL when
  !   the GFL variable YCVGQ is not available, or when this GFL variable
  !   is available in grid-point only?
  IF ((NCOMP_CVGQ == 1) .AND. ((.NOT. YCVGQ%LACTIVE).OR. (YCVGQ%LGP))) THEN
    WRITE(KULOUT,*)&
     & 'SUDYN WARNING: when NCOMP_CVGQ == 1 YCVGQ must be active and spectral'
    CALL ABOR1('SUDYN: NCOMP_CVGQ inconsistency ABOR1 CALLED')
  ENDIF
  ! * Try to compute a Semi-Lagrangian CVGQ with grid-point YCVGQ-GFL when
  !   the GFL variable YCVGQ is not available, or when this GFL variable
  !   is available in spectral space only?
  IF ((NCOMP_CVGQ == 2) .AND. ((.NOT. YCVGQ%LACTIVE).OR. (YCVGQ%LSP))) THEN
    WRITE(KULOUT,*)&
     & 'SUDYN WARNING: when NCOMP_CVGQ == 2 YCVGQ must be active and grid point'
    CALL ABOR1('SUDYN: NCOMP_CVGQ inconsistency ABOR1 CALLED')
  ENDIF
  ! * Try to use NCOMP_CVGQ>0 in an Eulerian model?
  IF ((NCOMP_CVGQ > 0) .AND. (.NOT.YDDYNA%LSLAG) .AND. (NSTOP>0.OR.LDFI)) THEN
    WRITE(KULOUT,*)&
     & 'SUDYN WARNING: in an Eulerian scheme only NCOMP_CVGQ=0 can be used'
    CALL ABOR1('SUDYN: NCOMP_CVGQ inconsistency ABOR1 CALLED')
  ENDIF
  ! * Try to use NCOMP_CVGQ=2 in a 3TLSL model (not yet coded)?
  IF ((NCOMP_CVGQ == 2) .AND. (.NOT.YDDYNA%LTWOTL)) THEN
    WRITE(KULOUT,*)&
     & 'SUDYN WARNING: NCOMP_CVGQ=2 is currently available in a SL2TL scheme only'
    CALL ABOR1('SUDYN: NCOMP_CVGQ inconsistency ABOR1 CALLED')
  ENDIF
ENDIF

!*       3.17  3D turbulence.

! ky: this checking must go in SU0PHY.
IF (L3DTURB.AND.(.NOT.LMPHYS)) THEN
  WRITE(KULOUT,*) 'L3DTURB=.T. works (for the moment) only with Meteo-France physics!'
  CALL ABOR1('SUDYN: ABOR1 CALLED')
ENDIF

!*       3.18  Convert frequency for update of gaseous transmissions
!              from hours to timesteps (this belongs to SU0PHY but TSTEP
!              is not known yet when SU0PHY is called)

! ky: in a OOPS frame, NTHRAYFR is model-dependent => its calculation,
!     and the piece of code below, must be moved in a routine called under SUPHY.
IF ( NRAY > 1 ) THEN
  ZTSTEP=MAX(TSTEP,1.0_JPRB)
  ZSTPHR=3600._JPRB/ZTSTEP
  IF ( NSORAYFR < 0 ) THEN
    NSORAYFR=-NSORAYFR*ZSTPHR+0.5_JPRB   
    WRITE(KULOUT,'('' NSORAYFR converted to timesteps = '',I5)') NSORAYFR
  ENDIF
  IF ( NTHRAYFR < 0 ) THEN
    NTHRAYFR=-NTHRAYFR*ZSTPHR+0.5_JPRB
    WRITE(KULOUT,'('' NTHRAYFR converted to timesteps = '',I5)') NTHRAYFR
  ENDIF
  IF ( (NSORAYFR /= 1.OR.NTHRAYFR /= 1.OR.NRAUTOEV > 1).AND. &
   & LDFI.AND.RTDFI /= TSTEP ) THEN
    CALL ABOR1('SUDYN: ACRANEB2 intermittency not allowed for RTDFI /= TSTEP!')
  ENDIF
ENDIF

!*       3.19  No implicit Coriolis for LELAM.AND.(.NOT.LMAP).

IF ( LELAM.AND.(.NOT.YDGEOMETRY%YREGEO%LMAP) ) THEN
  LADVF=.FALSE.
  RW2TLFF=0.0_JPRB
ENDIF

!*       3.20 Specific options for SL trajectory research or interpolation

IF (LSLTVWENO.OR.LRHSVWENO) THEN
  IF (NEDER < 1 .OR. NEDER > 5) THEN
    WRITE(KULOUT,*) ' Only available options for the smoothness indicators are NEDER=1/2/3/4/5 '
    CALL ABOR1('SUDYN: ABOR1 CALLED')
  ENDIF
  IF (LWENOBC .AND. (.NOT.YDCVER%LREGETA)) THEN
    WRITE(KULOUT,*) ' Smart boundary condition LWENOBC=true only works with LREGETA=true '
    CALL ABOR1('SUDYN: ABOR1 CALLED')
  ENDIF
ENDIF
IF (LELAM .AND. (LSLDP_CURV.OR.LRHS_CURV)) THEN
  WRITE(KULOUT,*) ' Options LSLDP_CURV and LRHS_CURV are essentially designed for spherical geometry!!! '
  WRITE(KULOUT,*) ' Being naturally in (x,y) coordinates there is little point to activate '
  WRITE(KULOUT,*) ' (lat,lon) -> (x,y,z) conversion for horizontal wind components interpolation'
  WRITE(KULOUT,*) ' when in LAM. This option however allows to skip all the expensive and '
  WRITE(KULOUT,*) ' precision sensitive rotational matrix computation which might be usefull also '
  WRITE(KULOUT,*) ' to plane geometry. So far however the relevant code does not exist. Sorry. '
  CALL ABOR1('SUDYN: ABOR1 CALLED')
ENDIF

IF (LLSTRSG .AND. LSLDP_CURV) THEN
  WRITE(KULOUT,*) ' Option LSLDP_CURV is designed for unstretched/untilted spherical geometry!!! '
  CALL ABOR1('SUDYN: ABOR1 CALLED')
ENDIF

IF (LSLDP_XYZ.AND..NOT.LSLDP_CURV) THEN
  WRITE(KULOUT,*) 'LSLDP_XYZ ONLY WORKS WITH LSLDP_CURV=T'
  CALL ABOR1('SUDYN: ABOR1 CALLED')
ENDIF

IF (LSLDP_SAVE.AND..NOT.LSLDP_XYZ) THEN
  WRITE(KULOUT,*) 'LSLDP_SAVE ONLY WORKS WITH LSLDP_XYZ=T'
  CALL ABOR1('SUDYN: ABOR1 CALLED')
ENDIF

IF (LSLHD .AND. LRHSVWENO) THEN
  IF (LSLHDSPONGE) THEN
    ! Bit tricky code for ECMWF sponge
    WRITE(KULOUT,*) 'SLHD loweres 4-point interpolation to second order only while WENO requires it'
    WRITE(KULOUT,*) 'to be of third order accurate. Combination LSLHD,LRHSVWENO=t,t effectively'
    WRITE(KULOUT,*) 'supress vertical SLHD part for all variables treated by WENO interpolation.'
    IF (SLHDKREF /= 0._JPRB) THEN
      WRITE(KULOUT,*) ' SLHDKREF /= 0. is not compatible with LSLHD,LRHSVWENO=t,t. '
      CALL ABOR1('SUDYN: ABOR1 CALLED')
    ENDIF
  ELSE
    ! Unless really necessary don't overcomplicate it
    WRITE(KULOUT,*) 'SLHD loweres 4-point interpolation to second order only while WENO requires it'
    WRITE(KULOUT,*) 'to be of third order accurate. Combination LSLHD,LRHSVWENO=t,t does not work.'
    CALL ABOR1('SUDYN: ABOR1 CALLED')
  ENDIF
ENDIF
IF ((LVWENO_W.OR.LVWENO_T.OR.LVWENO_SP.OR.LVWENO_SPD.OR.LVWENO_SVD).AND.(.NOT.LRHSVWENO)) THEN
  WRITE(KULOUT,*) 'By asking WENO interpolation for RHS the LRHSVWENO must be set to TRUE.'
  CALL ABOR1('SUDYN: ABOR1 CALLED')
ENDIF

IF (LSLDP_RK.AND.(.NOT.LSLDP_CURV)) THEN
  WRITE(KULOUT,*) 'Implicit RK4 scheme (LSLDP_RK=true) only works with LSLDP_CURV=true.'
  CALL ABOR1('SUDYN: ABOR1 CALLED')
ENDIF
IF (LSLDP_RK.AND.LSLDP_XYZ) THEN
  WRITE(KULOUT,*) 'Implicit RK4 scheme (LSLDP_RK=true) is not currently supported with LSLDP_XYZ=true.'
  CALL ABOR1('SUDYN: ABOR1 CALLED')
ENDIF

IF (NQMHOISLT > 3 .OR. NQMHOISLT < 0) THEN
  WRITE(KULOUT,*) 'Quasi-monotonic setting for wind interpolation during SL trajectory research'
  WRITE(KULOUT,*) 'only allows NQMHOISLT=0/1/2/3. '
  CALL ABOR1('SUDYN: ABOR1 CALLED')
ENDIF

IF ((LRHS_CURV .OR. LSLDP_CURV) .AND. (NSTTYP /= 1.OR.LLSTRSG)) THEN
  WRITE(KULOUT,*) 'LRHS_CURV/LSLDP_CURV ARE NOT CODED FOR TILTING OR STRETCHING.'
  CALL ABOR1('SUDYN: ABOR1 CALLED')
ENDIF

IF ( LRHS_CURV .AND. L2TLFF ) THEN
  WRITE(KULOUT,*) 'WARNING: LRHS_CURV/L2TLFF=T/T represent a hybrid solution when rotational matrix is still'
  WRITE(KULOUT,*) '         required for coriolis term while wind is being treated in Cartesian space.'
ENDIF

! Special security for SLVF in lowres TL test
IF (LSLVF.AND.(NCONF == 501).AND.(NFLEVG < 50)) THEN
  LSLVF=.false.
  WRITE(KULOUT,*) '=== WARNING === : LSLVF is being deactivated!'
  WRITE(KULOUT,*) ' It can not convert when low number of levels is requested.'
  WRITE(UNIT=KULOUT,FMT='('' LSLVF         = '',L2)') LSLVF
ENDIF

!     ------------------------------------------------------------------

!*       3B.  SET DIMENSION OF SL WEIGHTS 
!             ---------------------------

! Used in horizontal part of 3D turbulence and diffusion of phys. tendencies
IF (YDDYNA%LSLAG) THEN

  ! By default there's only one set of weights.
  NSLDIMK=1

  ! 3D turbulence adds two more sets for K_M and K_H
  IF (L3DTURB) NSLDIMK=NSLDIMK+2

  ! Diffusive interpolator for physics
  IF (NSPLTHOI == 1) NSLDIMK=NSLDIMK+1

ENDIF

!     ------------------------------------------------------------------

!*       4.    SET DYNAMICS STRUCTURES OF INTDYN_MOD
!              -------------------------------------

! Requires final calculation of NSPLTHOI and NSLDIMK.
! Must be called after calling SUGFL3.
! Useful if semi-Lagrangian scheme only.
IF (YDDYNA%LSLAG) CALL SUINTDYNSL(YDDYN,YDDYNA,YGFL,YDTLSCAW,YDTLSCAWH,YDTRSCAW,YDTRSCAWH,YDTSCO,YDTCCO)

!     ------------------------------------------------------------------

!*       5.    ALLOCATE SOME ARRAYS.
!              ---------------------

CALL SUALDYN(YDGEOMETRY%YRDIMV,YDDYNA,YDDYN)
IF (NSTOP > 0.OR.LDFI) THEN
  IF (LELAM) THEN
    CALL SUELDYNB(YDGEOMETRY,YDDYN,YDEDYN,YGFL)
  ELSE
    CALL SUALDYNB(YDGEOMETRY,YDMODEL%YRML_GCONF,YDDYN,YDDYNA)
  ENDIF
ENDIF

!-----------------------------------------------------------------------
! Compute model level to apply the stable composite scheme (LSETTLSVF=T)
! 
!-----------------------------------------------------------------------
NFLEVSF=0

IF (LSETTLSVF) THEN
  DO JLEV=1,NFLEVG
    IF (STPREH(JLEV)<= RPRES_SETTLSVF) NFLEVSF=NFLEVSF+1
  ENDDO
  NFLEVSF=MIN(NFLEVG,NFLEVSF)
ENDIF

!-----------------------------------------------------------------------
!  Initialize SLAVEPP related computation  
!-----------------------------------------------------------------------
! Initialize pointers to save ECMWF physical tendencies
CALL SUSAVTEND(NFLEVG,LSLPHY,YDSLPHY)
!  Initialize implicitness for ECMWF physics 
IF(LSLPHY)THEN
  SLPHY_P_LOW =15000._JPRB
  SLPHY_P_HIGH=10000._JPRB
  DO JLEV=1,NFLEVG
   !! mask
   !ZMASK=MIN(1._JPRB,MAX(0._JPRB,(SLPHY_P_LOW-STPRE(JLEV))/ &
   !  &                           (SLPHY_P_LOW-SLPHY_P_HIGH)))
   !YDSLPHY%RSLWX(JLEV)=0.5_JPRB*(ZMASK+1._JPRB)
    YDSLPHY%RSLWX(JLEV)=0.5_JPRB
  ENDDO
ENDIF

!-----------------------------------------------------------------------
! Setup Multiple grids advection
!-----------------------------------------------------------------------
#ifdef WITH_MGRIDS
YDMODEL%YRML_DYN%YR_MGRIDS_ADVECTION = MGRIDS_ADVECTION( YDRIP, YDDYN, YDGEOMETRY )
#endif

IF( YGFL%NFMG > 0 ) THEN
#ifdef WITH_MGRIDS
  IF( .NOT. YDMODEL%YRML_DYN%YR_MGRIDS_ADVECTION%ACTIVE() ) THEN
    CALL ABOR1('SUDYN: YGFL%NFMG>0 requires correct setup of MGRIDS_ADVECTION')
  ENDIF
#else
  CALL ABOR1('SUDYN: YGFL%NFMG>0 requires compilation with definition "WITH_MGRIDS"')
#endif
ENDIF 

!     ------------------------------------------------------------------

!*       6.    Print final values.
!              -------------------

WRITE(UNIT=KULOUT,FMT='('' Printings for module YOMDYN '')')

WRITE(KULOUT,*) '* SUDYN: Asselin temporal filter.'
WRITE(UNIT=KULOUT,FMT='('' REPS1  = '',E14.8,'' REPS2  = '',E14.8)') REPS1,REPS2
WRITE(UNIT=KULOUT,FMT='('' REPSM1  = '',E14.8,'' REPSM2  = '',E14.8&
 & ,'' REPSP1  = '',E14.8)') REPSM1,REPSM2,REPSP1

WRITE(KULOUT,*) '* SUDYN: Spectral horizontal diffusion.'
WRITE(UNIT=KULOUT,FMT='('' LNEWHD  = '',L2,'' RRDXTAU = '',E14.8)') LNEWHD,RRDXTAU
WRITE(UNIT=KULOUT,FMT='('' RDAMPVOR = '',E14.8,'' RDAMPDIV = '',E14.8)') RDAMPVOR,RDAMPDIV
WRITE(UNIT=KULOUT,FMT='('' RDAMPT   = '',E14.8,'' RDAMPQ   = '',E14.8)') RDAMPT,RDAMPQ
WRITE(UNIT=KULOUT,FMT='(A,E14.8)') ' RDAMPO3 = ', RDAMPO3
WRITE(UNIT=KULOUT,FMT='('' RDAMPPD  = '',E14.8,'' RDAMPVD  = '',E14.8)') RDAMPPD,RDAMPVD
WRITE(UNIT=KULOUT,FMT='('' RDAMPSP  = '',E14.8)') RDAMPSP
WRITE(UNIT=KULOUT,FMT='('' HDIRVOR = '',E14.8,'' HDIRDIV = '',E14.8)') HDIRVOR,HDIRDIV
WRITE(UNIT=KULOUT,FMT='('' HDIRT   = '',E14.8,'' HDIRQ   = '',E14.8)') HDIRT,HDIRQ
WRITE(UNIT=KULOUT,FMT='(A,E14.8)') ' HDIRO3 = ', HDIRO3
WRITE(UNIT=KULOUT,FMT='('' HDIRPD  = '',E14.8,'' HDIRVD  = '',E14.8)') HDIRPD,HDIRVD
WRITE(UNIT=KULOUT,FMT='('' HDIRSP  = '',E14.8)') HDIRSP
WRITE(UNIT=KULOUT,FMT='('' REXPDH = '',E14.8,'' FRANDH = '',E14.8)') REXPDH,FRANDH
WRITE(UNIT=KULOUT,FMT='('' SLEVDH = '',E14.8,'' SLEVDH3 = '',E14.8)') SLEVDH,SLEVDH3
WRITE(UNIT=KULOUT,FMT='('' SLEVDH1= '',E14.8,'' SLEVDH2 = '',E14.8,'' RATIO_HDI_TOP = '',E14.8,'' NSREFDH = '',I5)') &
 & SLEVDH1,SLEVDH2,RATIO_HDI_TOP,NSREFDH
WRITE(UNIT=KULOUT,FMT='('' NPROFILEHD= '',I2)') NPROFILEHD
WRITE(UNIT=KULOUT,FMT='('' RPROFHDBT = '',E14.8,'' RPROFHDTP = '',E14.8&
 & ,'' RPROFHDEX = '',E14.8,'' RPROFHDMX = '',E14.8)') RPROFHDBT,RPROFHDTP,RPROFHDMX,RPROFHDEX
WRITE(UNIT=KULOUT,FMT='('' LRDISPE_EC= '',L2)') LRDISPE_EC
WRITE(UNIT=KULOUT,FMT='('' LSTRHD = '',L2)') LSTRHD
WRITE(UNIT=KULOUT,FMT='('' LRHDI_LASTITERPC = '',L2)') LRHDI_LASTITERPC

WRITE(KULOUT,*) '* SUDYN: diffusion quantities.'
WRITE(UNIT=KULOUT,FMT='('' LSPECVIS = '',L2)') LSPECVIS
WRITE(UNIT=KULOUT,FMT='('' LHDIFFM = '',L2)') LHDIFFM
WRITE(UNIT=KULOUT,FMT='('' NDIFFACT =  '',I1)') NDIFFACT
WRITE(UNIT=KULOUT,FMT='('' LGPSTRESS = '',L2)') LGPSTRESS
WRITE(UNIT=KULOUT,FMT='('' RCLSTRESS = '',E10.4)') RCLSTRESS
WRITE(UNIT=KULOUT,FMT='('' RCLPOLE = '',E10.4)') RCLPOLE
WRITE(UNIT=KULOUT,FMT='('' LTOP_VOR = '',L2)') LTOP_VOR
WRITE(UNIT=KULOUT,FMT='('' NTOP_VOR_TRUNC =  '',I6)') NTOP_VOR_TRUNC
WRITE(UNIT=KULOUT,FMT='('' NTOP_VOR_BOT =  '',I3)') NTOP_VOR_BOT
IF (LSLHD) THEN
  WRITE(KULOUT,*) '* SUDYN: semi-Lagrangian diffusion (SLHD).'
  WRITE(UNIT=KULOUT,FMT='('' SLHDA0  = '',E14.8)') SLHDA0
  WRITE(UNIT=KULOUT,FMT='('' SLHDA0T = '',E10.4)') SLHDA0T
  WRITE(UNIT=KULOUT,FMT='('' SLHDB   = '',E14.8)') SLHDB
  WRITE(UNIT=KULOUT,FMT='('' SLHDBT  = '',E10.4)') SLHDBT
  WRITE(UNIT=KULOUT,FMT='('' SLHDD00 = '',E14.8)') SLHDD00
  WRITE(UNIT=KULOUT,FMT='('' SLHDD00T= '',E14.8)') SLHDD00T
  WRITE(UNIT=KULOUT,FMT='('' SLHDDIV = '',E14.8)') SLHDDIV
  WRITE(UNIT=KULOUT,FMT='('' SLHDRATDDIV = '',E14.8)') SLHDRATDDIV
  WRITE(UNIT=KULOUT,FMT='('' SLHDHOR = '',E14.8)') SLHDHOR
  WRITE(UNIT=KULOUT,FMT='('' SLHD_P_LOW = '',E14.8)') SLHD_P_LOW
  WRITE(UNIT=KULOUT,FMT='('' SLHD_P_HIGH = '',E14.8)') SLHD_P_HIGH
  WRITE(UNIT=KULOUT,FMT='('' ZSLHDP1 = '',E14.8,'' ZSLHDP3 = '',E14.8)') ZSLHDP1,ZSLHDP3
  WRITE(UNIT=KULOUT,FMT='('' REXPDHS = '',E14.8)') REXPDHS
  WRITE(UNIT=KULOUT,FMT='('' SLEVDHS = '',E14.8)') SLEVDHS
  WRITE(UNIT=KULOUT,FMT='('' SLEVDHS1= '',E14.8,'' SLEVDHS2 = '',E14.8)') SLEVDHS1,SLEVDHS2
  WRITE(UNIT=KULOUT,FMT='('' SDRED   = '',E14.8)') SDRED
  WRITE(UNIT=KULOUT,FMT='('' LSLHDHEAT= '',L2)') LSLHDHEAT
  WRITE(UNIT=KULOUT,FMT='('' LSLHDSPONGE= '',L2)') LSLHDSPONGE
  WRITE(UNIT=KULOUT,FMT='('' RDAMPDIVS= '',E14.8)') RDAMPDIVS
  WRITE(UNIT=KULOUT,FMT='('' RDAMPVORS= '',E14.8)') RDAMPVORS
  WRITE(UNIT=KULOUT,FMT='('' RDAMPVDS= '',E14.8)') RDAMPVDS
  WRITE(UNIT=KULOUT,FMT='('' RDAMPHDS = '',E14.8)') RDAMPHDS
  IF (HRDSRVOR /= 0.0_JPRB) THEN
    WRITE(UNIT=KULOUT,FMT='('' 1/HRDSRVOR = '',E14.8)') 1.0_JPRB/HRDSRVOR
  ENDIF
  IF (HRDSRDIV /= 0.0_JPRB) THEN
    WRITE(UNIT=KULOUT,FMT='('' 1/HRDSRDIV = '',E14.8)') 1.0_JPRB/HRDSRDIV
  ENDIF
  IF (YDDYNA%LNHDYN .AND. HRDSRVD /= 0.0_JPRB) THEN
    WRITE(UNIT=KULOUT,FMT='('' 1/HRDSRVD = '',E14.8)') 1.0_JPRB/HRDSRVD
  ENDIF
ENDIF

WRITE(KULOUT,*) '* SUDYN: Semi-implicit scheme.'
WRITE(UNIT=KULOUT,FMT='('' LSIDG = '',L2)') LSIDG
WRITE(UNIT=KULOUT,FMT='('' BETADT = '',E14.8)') BETADT
WRITE(UNIT=KULOUT,FMT='('' SITR   = '',E14.8,'' SITRA  = '',E14.8,'' SIPR   = '',E14.8)') SITR,SITRA,SIPR
WRITE(UNIT=KULOUT,FMT='('' NOPT_SITRA= '',I2)') NOPT_SITRA
WRITE(UNIT=KULOUT,FMT='('' REFGEO = '',E14.8)') REFGEO
WRITE(UNIT=KULOUT,FMT='('' NITERHELM= '',I2)') NITERHELM
WRITE(UNIT=KULOUT,FMT='('' LIMPF = '',L2)') LIMPF

WRITE(KULOUT,*) '* SUDYN: ICI (PC) scheme.'
WRITE(UNIT=KULOUT,FMT='('' NSITER = '',I2)') NSITER

WRITE(KULOUT,*) '* SUDYN: Semi-Lagrangian scheme.'
WRITE(UNIT=KULOUT,FMT='('' NITMP= '',I2)') NITMP
WRITE(UNIT=KULOUT,FMT='('' NVLAG  = '',I8,'' NWLAG  = '',I8,'' NTLAG = '',I8)') NVLAG,NWLAG,NTLAG
WRITE(UNIT=KULOUT,FMT='('' NSVDLAG = '',I8,'' NSPDLAG = '',I8)') NSVDLAG,NSPDLAG
WRITE(UNIT=KULOUT,FMT='('' NSPLTHOI = '',I8)') NSPLTHOI
WRITE(UNIT=KULOUT,FMT='('' LQMW = '',L2,'' LQMHW  = '',L2)') LQMW,LQMHW
WRITE(UNIT=KULOUT,FMT='('' LQMT = '',L2,'' LQMHT  = '',L2)') LQMT,LQMHT
WRITE(UNIT=KULOUT,FMT='('' LQMP = '',L2,'' LQMHP  = '',L2)') LQMP,LQMHP
WRITE(UNIT=KULOUT,FMT='('' LQMPD = '',L2,'' LQMHPD  = '',L2)') LQMPD,LQMHPD
WRITE(UNIT=KULOUT,FMT='('' LQMVD = '',L2,'' LQMHVD  = '',L2)') LQMVD,LQMHVD
WRITE(UNIT=KULOUT,FMT='('' LVWENO_W = '',L2,'' LVWENO_T  = '',L2)') LVWENO_W,LVWENO_T
WRITE(UNIT=KULOUT,FMT='('' LVWENO_SP = '',L2)') LVWENO_SP
WRITE(UNIT=KULOUT,FMT='('' LVWENO_SPD = '',L2,'' LVWENO_SVD  = '',L2)') LVWENO_SPD,LVWENO_SVD
WRITE(UNIT=KULOUT,FMT='('' WENO_ALPHA_W = '',F16.12,'' WENO_ALPHA_T  = '',F16.12)') WENO_ALPHA_W,WENO_ALPHA_T
WRITE(UNIT=KULOUT,FMT='('' WENO_ALPHA_SP = '',F16.12)') WENO_ALPHA_SP
WRITE(UNIT=KULOUT,FMT='('' WENO_ALPHA_SPD = '',F16.12,'' WENO_ALPHA_SVD  = '',F16.12)') &
  &                                       WENO_ALPHA_SPD,WENO_ALPHA_SVD
WRITE(UNIT=KULOUT,FMT='('' RINTOPT = '',F16.5,'' RRINTOPT  = '',F16.12)') RINTOPT,RRINTOPT
WRITE(UNIT=KULOUT,FMT='('' Provisional VMAX1 = '',E14.8,'' Provisional VMAX2  = '',E14.8)') VMAX1,VMAX2
WRITE(UNIT=KULOUT,FMT='('' N.B.: VMAX1, VMAX2 WILL BE RECOMPUTED IN SUSC2B'')')
WRITE(UNIT=KULOUT,FMT='('' LBOUND_D3 = '',L2)') LBOUND_D3
WRITE(UNIT=KULOUT,FMT='('' RMAX_D3 = '',E14.8)') RMAX_D3
WRITE(UNIT=KULOUT,FMT='('' VETAON = '',F16.12,'' VETAOX  = '',F16.12)') VETAON,VETAOX
WRITE(UNIT=KULOUT,FMT='('' VESL  = '',E14.8)') VESL
WRITE(UNIT=KULOUT,FMT='('' XIDT  = '',E14.8)') XIDT
WRITE(UNIT=KULOUT,FMT='('' RALPHA = '',E14.8)') RALPHA
WRITE(UNIT=KULOUT,FMT='('' RALPHA_TOP = '',E14.8)') RALPHA_TOP
WRITE(UNIT=KULOUT,FMT='('' NLEV_ZALPHA = '',I2)') NLEV_ZALPHA
WRITE(UNIT=KULOUT,FMT='('' KLSPONGE = '',E14.8)') KLSPONGE
WRITE(UNIT=KULOUT,FMT='('' RBEGSPONGE = '',E14.8)') RBEGSPONGE
WRITE(UNIT=KULOUT,FMT='('' LM0DAMP = '',L2)') LM0DAMP
WRITE(UNIT=KULOUT,FMT='('' NEDER = '',I2)') NEDER
WRITE(UNIT=KULOUT,FMT='('' LWENOBC = '',L2)') LWENOBC
IF (LHOISLT) WRITE(UNIT=KULOUT,FMT='('' NQMHOISLT = '',I2)') NQMHOISLT
WRITE(UNIT=KULOUT,FMT='('' LSLDP_RK = '',L2)') LSLDP_RK
WRITE(UNIT=KULOUT,FMT='('' LSLDP_CURV = '',L2)') LSLDP_CURV
WRITE(UNIT=KULOUT,FMT='('' LRHS_CURV = '',L2)') LRHS_CURV
WRITE(UNIT=KULOUT,FMT='('' LSLDP_XYZ = '',L2)') LSLDP_XYZ
WRITE(UNIT=KULOUT,FMT='('' LSLDP_SAVE = '',L2)') LSLDP_SAVE
WRITE(UNIT=KULOUT,FMT='('' RW2TLFF = '',E14.8)') RW2TLFF
WRITE(UNIT=KULOUT,FMT='('' L2TLFF = '',L2)') L2TLFF
WRITE(UNIT=KULOUT,FMT='('' LSVTSM = '',L2)') LSVTSM
WRITE(UNIT=KULOUT,FMT='('' RPRES_SVTSM ='',E14.8)') RPRES_SVTSM
WRITE(UNIT=KULOUT,FMT='('' LSETTLSVF = '',L2)') LSETTLSVF
IF (LSETTLSVF) THEN
  WRITE(UNIT=KULOUT,FMT='('' RPRES_SETTLSVF ='',E14.8)') RPRES_SETTLSVF
  WRITE(UNIT=KULOUT,FMT='('' SETTLSTF SCHEME APPLIED AT LEVS 1 TO:'',I3)') NFLEVSF
  WRITE(UNIT=KULOUT,FMT='('' RSCALE ='',E14.8)') RSCALE
  WRITE(UNIT=KULOUT,FMT='('' RSCALEOFF ='',E14.8)') RSCALEOFF
  WRITE(UNIT=KULOUT,FMT='('' LSETTLSVF_DIFF = '',L2)') LSETTLSVF_DIFF
ENDIF
IF (LMASCOR) THEN
  WRITE(UNIT=KULOUT,FMT='('' LMASCOR = '',L2)') LMASCOR
  WRITE(UNIT=KULOUT,FMT='('' LGPMASCOR = '',L2)') LGPMASCOR
  WRITE(UNIT=KULOUT,FMT='('' LMASDRY = '',L2)') LMASDRY
ENDIF
WRITE(UNIT=KULOUT,FMT='('' NVSEPC = '',I8,'' NVSEPL = '',I8)') NVSEPC, NVSEPL
IF (YDDYNA%LSLAG) THEN
  WRITE(UNIT=KULOUT,FMT='('' NSLDIMK = '',I2)') NSLDIMK
ENDIF

WRITE(KULOUT,*) '* SUDYN: change of variable in dynamics.'
WRITE(UNIT=KULOUT,FMT='('' LADVF = '',L2)') LADVF
WRITE(UNIT=KULOUT,FMT='('' RCMSLP0 ='',E14.8)') RCMSLP0

WRITE(KULOUT,*) '* SUDYN: option LRUBC.'
WRITE(UNIT=KULOUT,FMT='('' NRUBC = '',I1)') NRUBC
IF (LRUBC) THEN
  WRITE(UNIT=KULOUT,FMT='('' RTEMRB= '',E14.8)') RTEMRB
  WRITE(UNIT=KULOUT,FMT='('' SITRUB= '',E14.8)') SITRUB
  WRITE(UNIT=KULOUT,FMT='('' SIPRUB= '',E14.8)') SIPRUB
ENDIF

WRITE(KULOUT,*) '* SUDYN: miscellaneous.'
WRITE(UNIT=KULOUT,FMT='('' VNORM  = '',F6.3)') VNORM
WRITE(UNIT=KULOUT,FMT='('' NCOMP_CVGQ = '',I1)') NCOMP_CVGQ
WRITE(UNIT=KULOUT,FMT='('' LRFRIC = '',L2)') LRFRIC
WRITE(UNIT=KULOUT,FMT='('' LRFRICISOTR = '',L2)') LRFRICISOTR

WRITE(KULOUT,*) '* SUDYN: moist/dry air conservation'
WRITE(UNIT=KULOUT,FMT='('' LNODRYFLX= '',L2)') LNODRYFLX

WRITE(UNIT=KULOUT,FMT='('' '')')

!     ------------------------------------------------------------------

!*       7.    INITIALIZE SEMI-IMPLICIT AND VERTICAL NORMAL MODES.
!              ---------------------------------------------------

! ky: In the 3D model, when NSTOP=0, the setup of the SI scheme should be
!     called each time we need to call CPG at instant 0 (for example to do
!     diagnostics). This is the case for example if LSDDH=T.
!     The test should be consistent with CNT4 code.

IF (YDDYNA%LSLAG) THEN
  IF (LR2D) THEN
    ! XIDT>0 not coded in 2D models.
    RBT=(1.0_JPRB+VESL)*BETADT
  ELSE
    IF (YDDYNA%LTWOTL) THEN
      RBT=(1.0_JPRB+VESL)*BETADT*(1.0_JPRB+XIDT)
    ELSE
      RBT=(1.0_JPRB+VESL)*BETADT
    ENDIF
  ENDIF
ELSE
  RBT=BETADT
ENDIF
RBTS2=0.5_JPRB*RBT

IF (NOPT_SITRA == 1) THEN
  DO JLEV=1,NFLEVG
    YDDYN%SITRAM(JLEV)=YDDYN%SITRA*(STTEM(JLEV)/STTEM(NFLEVG))
  ENDDO
ELSE
  DO JLEV=1,NFLEVG
    YDDYN%SITRAM(JLEV)=YDDYN%SITRA
  ENDDO
ENDIF

IF (.NOT.LELAM) THEN
  IF (.NOT.LR2D) THEN
    ! Setup of the SI scheme for 3D models.
    IF ((NSTOP > 0).OR.LDFI.OR.(NSTOP==0.AND.LSDDH).OR.(NSTOP==0.AND.LWCOU)) THEN
      IF (YDDYNA%LNHEE.AND.(.NOT.LSI_NHEE)) THEN
        SITIME=-1.0_JPRB
        CALL SUNHSI(YDGEOMETRY,YDRIP,YDDYN,YDEDYN,YDDYNA%LNESC,KULOUT,.TRUE.)

        ! Test of convergence of the iterative Helmholtz solver if needed
        LLNH_NOC1=YDDYNA%LNHEE.AND.(.NOT.(YDCVER%NDLNPR == 1))
        IF (LLNH_NOC1 .AND. NFLEVG>1) THEN
          CALL SUNHSI_TESTCONV(YDGEOMETRY,YDRIP,YDDYN,YDDYNA)
        ENDIF
      ELSEIF (YDDYNA%LNHEE.AND.LSI_NHEE) THEN
        SITIME=-1.0_JPRB
        CALL SUNHEESI(YDGEOMETRY,YDRIP,YDDYN,YDDYNA,KULOUT,.TRUE.)
      ELSEIF (YDDYNA%LNHQE.AND.(.NOT.LNHQE_SIHYD)) THEN
        SITIME=-1.0_JPRB
        CALL SUNHQESI(YDGEOMETRY,YDRIP,YDDYN,YDDYNA,KULOUT,.TRUE.)
      ELSE
        CALL SUSI(YDGEOMETRY,YDRIP,YDDYN,YDDYNA,YDEDYN,KULOUT)
      ENDIF
    ENDIF
  ELSE
    ! Setup of the SI scheme for 2D models.
    IF (NSTOP > 0.OR.LDFI) THEN
      YDDYN%SIVP(1)=REFGEO
      IF (LSIDG) THEN
        SITIME=-1.0_JPRB
        CALL SUHEG(YDGEOMETRY,YDRIP,YDDYN,YDDYNA)
      ENDIF
      YDDYN%SIBI(1,1)=1.0_JPRB/YDDYN%SIVP(1)
    ENDIF
  ENDIF ! LR2D
ELSE
  IF (.NOT.LR2D) THEN
    ! Setup of the SI scheme for 3D models.
    IF ((NSTOP > 0).OR.LDFI.OR.(NSTOP==0.AND.LSDDH).OR.(NSTOP==0.AND.LWCOU)) THEN
      IF (YDDYNA%LNHEE.AND.(.NOT.LSI_NHEE)) THEN
        SITIME=-1.0_JPRB
        CALL SUNHSI(YDGEOMETRY,YDRIP,YDDYN,YDEDYN,YDDYNA%LNESC,KULOUT,.TRUE.)

        ! Test of convergence of the iterative Helmholtz solver if needed
        LLNH_NOC1=YDDYNA%LNHEE.AND.(.NOT.(YDCVER%NDLNPR == 1))
        IF (LLNH_NOC1 .AND. NFLEVG>1) THEN
          CALL SUNHSI_TESTCONV(YDGEOMETRY,YDRIP,YDDYN,YDDYNA)
        ENDIF
      ELSEIF (YDDYNA%LNHEE.AND.LSI_NHEE) THEN
        SITIME=-1.0_JPRB
        CALL SUNHEESI(YDGEOMETRY,YDRIP,YDDYN,YDDYNA,KULOUT,.TRUE.)
      ELSEIF (YDDYNA%LNHQE.AND.(.NOT.LNHQE_SIHYD)) THEN
        SITIME=-1.0_JPRB
        CALL SUNHQESI(YDGEOMETRY,YDRIP,YDDYN,YDDYNA,KULOUT,.TRUE.)
      ELSE
        CALL SUSI(YDGEOMETRY,YDRIP,YDDYN,YDDYNA,YDEDYN,KULOUT)
      ENDIF
    ENDIF
  ELSE
    ! LELAM shallow-water model not implemented.
    CALL ABOR1('SUDYN: ABOR1 CALLED')
  ENDIF ! LR2D
ENDIF ! .NOT.LELAM

!     ------------------------------------------------------------------

!*       8.   INITIALIZE HORIZONTAL DIFFUSION.
!             --------------------------------

IF (NSTOP > 0.OR.LDFI) THEN

  !*     8.1. resolution correction for SLHD

  IF (LSLHD.AND.(.NOT.LSLHD_STATIC)) THEN

    ! rescaling of SLHDA with respect to the generalised (horizontal)
    !  smallest wavenumber times sqrt(2.)
    ! and computing of SLHDD0 with respect of mesh-size (after Leith(1969))

    ZLXY=52349.828_JPRB  ! reference mesh-size at a location where M=1

    ! * NOTE: The proposed setting is valid only for resolutions finer
    !   than 40-50 km of horizontal mesh. When a coarser resolution
    !   is used, the linear formula for SLHDD0 (and SLHDA) has to be
    !   generalized by a saturation above the mentioned limit.
    !   Since it is quite easy to meet this limit when LLSTRSG=.T., 
    !   a temporal fix has been implemented for this option. It has
    !   been experimantally proven as fine for the operational ARPEGE
    !   model configuration at MF (June, 2006). For the other geometries
    !   with LLSTRSG=.T. the SLHD performance is not garranted.
    IF (LLSTRSG) THEN
      YDDYN%SLHDA(1:NGPTOT) =SLHDA0 * ( ZLXY * YDGSGEOM_NB%GM(1:NGPTOT) / SLHDP )**(ZSLHDP1)
      YDDYN%SLHDD0(1:NGPTOT)=SLHDD00 *( ZLXY * YDGSGEOM_NB%GM(1:NGPTOT) / SLHDP )**(ZSLHDP3)&
       & * SQRT(RSTRET/YDGSGEOM_NB%GM(1:NGPTOT))
    ELSE
      ! LAM model or unstretched ARP/IFS.
      YDDYN%SLHDA(1:NGPTOT) =SLHDA0 * ( ZLXY * YDGSGEOM_NB%GM(1:NGPTOT) / SLHDP )**(ZSLHDP1)
      YDDYN%SLHDD0(1:NGPTOT)=SLHDD00 *( ZLXY * YDGSGEOM_NB%GM(1:NGPTOT) / SLHDP )**(ZSLHDP3)
    ENDIF

  ENDIF

  ! masking functions for SLHD
  IF (LSLHD) THEN
    ALLOCATE (YDDYN%SLHD_MASK_U(NFLEVG))
    IF (LSLHDHEAT) ALLOCATE (YDDYN%SLHD_MASK_T(NFLEVG))
    NLEV_SPONGE=NFLEVG

    ! NOTE: The allocate statement effectively prevents any associate name to be used for 
    !       those structures. To prevent a run time error it is necessary to use exclusively
    !       full names of those variables.
    IF (LSLHDSPONGE) THEN
      ! SLHD is only used in high atmosphere to prevent model from grid-point storms there
      CALL SET_SLHD_SPONGE(YDGEOMETRY,YDDYN)
    ELSE
      YDDYN%SLHD_MASK_U(1:NFLEVG)=1._JPRB
      IF (LSLHDHEAT) YDDYN%SLHD_MASK_T(1:NFLEVG)=1._JPRB
    ENDIF
    WRITE(UNIT=KULOUT,FMT='('' NLEV_SPONGE = '',I3)') NLEV_SPONGE
  ENDIF

  !*     8.2. spectral HD

  IF (LELAM) THEN
    IF (LR2D) THEN
      ! LELAM shallow-water model not implemented.
      CALL ABOR1('SUDYN: ABOR1 CALLED')
    ELSE
      CALL SUEHDF(YDGEOMETRY,YDDYN,YDEDYN,YDMODEL%YRML_GCONF)
    ENDIF
  ELSE
    IF (LECMWF) THEN
      CALL SUHDF_EC(YDGEOMETRY,YDMODEL%YRML_GCONF,YDDYN,YDDYNA)
    ELSE
      IF (LR2D) THEN
        CALL SUHDF2(YDGEOMETRY%YRDIM,YDDYN)
      ELSE
        CALL SUHDF(YDGEOMETRY,YDMODEL%YRML_GCONF,YDDYN,YDDYNA)
      ENDIF
    ENDIF
  ENDIF

  IF (LSTRHD) THEN
    HDTIME_STRHD=-1.0_JPRB
    CALL SUHDU(YDGEOMETRY,YDMODEL%YRML_GCONF,YDDYN,YDDYNA)
  ENDIF

ENDIF

!     ------------------------------------------------------------------

!*       9.   INITIALIZE RCORDI... ARRAYS.
!             ----------------------------

CALL SURCORDI(YDGEOMETRY,YDDYN)

!     ------------------------------------------------------------------

!*       10.   INITIALIZE RKRF (RAYLEIGH FRICTION).
!              ------------------------------------

CALL SURAYFRIC(YDVAB,YDGEOMETRY%YRDIMV,YDDYN)

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUDYN',1,ZHOOK_HANDLE)
END SUBROUTINE SUDYN
