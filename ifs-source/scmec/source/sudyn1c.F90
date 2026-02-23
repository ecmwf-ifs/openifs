! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SUDYN1C(YDGEOMETRY,YDMODEL)

!**** *SUDYN1C* - Initialize the dynamics

!     Purpose.
!     --------
!           Initialize YOMDYN

!**   Interface.
!     ----------
!          *CALL* *SUDYN1C

!        Explicit arguments :
!        --------------------

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Called by SU0YOM1C

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the SCM

!     Author.
!     -------
!        Joao Teixeira   *ECMWF*

!     Modifications.
!     --------------
!        Original      94-01-21
!        J.Teixeira    Jun.-95   quasi-monotone interpolation and 
!                                cloud variables logicals
!                                for semi-lag. vert. advection.
!        M. Ko"hler    6-6-2006  Single Column Model integration within IFS
!        G. Carver     Nov 2012  Updated SCM version for 38r2.
!------------------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE TYPE_MODEL   , ONLY : MODEL
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN   , ONLY : NULNAM, NULOUT
!USE YOMCVER  , ONLY : &
! & LVFE_GW, LVFE_Z_TERM, LVFE_GWMPA, LVFE_DELNHPRE

IMPLICIT NONE

TYPE(GEOMETRY), INTENT(INOUT) :: YDGEOMETRY
TYPE(MODEL),TARGET,INTENT(INOUT) :: YDMODEL


!!LOGICAL :: LLOLDALPHA
LOGICAL :: LOLDALPHA
LOGICAL :: LL2TLFFX
LOGICAL :: LLSTRSG ! .T. => spherical stretched geometry.
LOGICAL :: LLCOMPUTE_CVGQ
LOGICAL :: LLNH_NOC1

INTEGER(KIND=JPIM) :: JLEV

REAL(KIND=JPRB) :: ZEPS, ZSLHDP1, ZSLHDP3, ZLXY, ZSTPHR, ZTSTEP, ZNOTUSED

REAL(KIND=JPHOOK) ::    ZHOOK_HANDLE

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
 & SLHDA0T, SLHDBT, SLHDD00, SLHDD00T, RPRES_SETTLSVF, RSCALE, RSCALEOFF, RALPHA,&
 & RALPHA_TOP, RATIO_HDI_TOP, WENO_ALPHA_SP, WENO_ALPHA_SPD, WENO_ALPHA_SVD, WENO_ALPHA_T, WENO_ALPHA_W, &
 & KLSPONGE, RBEGSPONGE, SLHD_P_LOW, SLHD_P_HIGH, RINTOPT, &
 & RSLOPE_MAX

LOGICAL,POINTER :: LSLHDHEAT, LSLHDSPONGE, LSETTLSVF, LSETFSTAT, LSETTLSVF_DIFF, &
   LGPMASCOR, LSPECVIS, LTOP_VOR

INTEGER(KIND=JPIM),POINTER :: NGPMASCOR, NTOP_VOR_TRUNC, NTOP_VOR_BOT, NOPT_SITRA, NLEV_ZALPHA, NQMHOISLT

INTEGER(KIND=JPIM), POINTER :: NVLAG, NITMP, NWLAG, NTLAG, NSPDLAG, NEDER,&
 & NSVDLAG, NSPLTHOI, NSITER, NSREFDH, NCOMP_CVGQ, NITERHELM, NPROFILEHD,&
 & NDIFFACT, NVSEPC, NVSEPL, NFLEVSF

LOGICAL, POINTER :: LADVF, LDYN_STABAN, LIMPF, LNEWHD,&
 & LQMW, LQMHW, LQMT, LQMHT, LQMP, LQMHP, LQMPD, LQMHPD, LQMVD, LQMHVD, &
 & LRHDI_LASTITERPC, LSVTSM, LRDISPE_EC, LHDIFFM, LGPSTRESS, LM0DAMP,&
 & LMASDRY, LMASCOR, LBOUND_D3, LRFRICISOTR, LRFRIC, LWENOBC, LSLDP_CURV, &
 & LSLDP_CURV_FIX, LSLDP_RK, &
 & LSLDP_XYZ, LSLDP_SAVE, LRHS_CURV, LRHS_CURV_FIX, LNODRYFLX, &
 & LVWENO_W, LVWENO_T, LVWENO_SP, LVWENO_SPD, LVWENO_SVD

#include "namdyn.nam.h"

!     ------------------------------------------------------------------
#include "posnam.intfb.h"
!#include "susta.intfb.h"
#include "surayfric.intfb.h"
!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUDYN1C',0,ZHOOK_HANDLE)
ASSOCIATE(YDDYN=>YDMODEL%YRML_DYN%YRDYN, &
& LSLPHY=>YDMODEL%YRML_PHY_EC%YREPHY%LSLPHY, &
& LSLAG=>YDMODEL%YRML_DYN%YRDYNA%LSLAG, &
& YRSLPHY=>YDMODEL%YRML_PHY_G%YRSLPHY)

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
LSLDP_RK         => YDDYN%LSLDP_RK
LSLDP_CURV       => YDDYN%LSLDP_CURV
!LSLDP_CURV_FIX   => YDDYN%LSLDP_CURV_FIX
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
NCOMP_CVGQ       => YDDYN%NCOMP_CVGQ
NITERHELM        => YDDYN%NITERHELM
NPROFILEHD       => YDDYN%NPROFILEHD
NDIFFACT         => YDDYN%NDIFFACT
NVSEPC           => YDDYN%NVSEPC
NVSEPL           => YDDYN%NVSEPL
NFLEVSF          => YDDYN%NFLEVSF
RSCALE           => YDDYN%RSCALE
RSCALEOFF        => YDDYN%RSCALEOFF
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
LGPMASCOR        => YDDYN%LGPMASCOR
NGPMASCOR        => YDDYN%NGPMASCOR
LNODRYFLX        => YDDYN%LNODRYFLX

LRHS_CURV        => YDDYN%LRHS_CURV
!LRHS_CURV_FIX    => YDDYN%LRHS_CURV_FIX
NQMHOISLT        => YDDYN%NQMHOISLT
SLHD_P_LOW       => YDDYN%SLHD_P_LOW
SLHD_P_HIGH      => YDDYN%SLHD_P_HIGH

!     ------------------------------------------------------------------

!*       1.    Set default values.
!              -------------------

REPS1   = .1_JPRB
REPS2   = .1_JPRB
REPSP1  = .1_JPRB
REPSM1  = .1_JPRB
REPSM2  = .1_JPRB

LQMW    = .FALSE.
LQMT    = .FALSE.
!LQMQ    = .TRUE.  ! No longer in yomdyn. Set directly in LARCIN1C
!LQMV    = .FALSE.

ZNOTUSED=-9999._JPRB
HDIRVOR= ZNOTUSED
HDIRDIV= ZNOTUSED
HDIRT  = ZNOTUSED
HDIRQ  = ZNOTUSED
HDIRO3 = ZNOTUSED
HDIRSP = ZNOTUSED
HDIRPD = ZNOTUSED
HDIRVD = ZNOTUSED

!*       1.18  Rayleigh friction (if active)

!  Defaults still same as 36r1, to be revised later.
!  Suggested values: RRFZ1=68, RRFPLM=100

RRFZ1=61._JPRB
RRFPLM=990._JPRB

YDGEOMETRY%YRCVER%LVERTFE=.FALSE.

! Suggested setting for SL advection
NITMP = 3
YDMODEL%YRML_DYN%YRDYNA%LSETTLS=.TRUE.  ! is YOMDYNA now

!     ------------------------------------------------------------------

!*       2.    READ NAMELIST.
!              --------------

CALL POSNAM(NULNAM,'NAMDYN')
READ(NULNAM,NAMDYN)


WRITE (NULOUT,*) ' NITMP= ',NITMP,'   LSETTLS= ',YDMODEL%YRML_DYN%YRDYNA%LSETTLS
!     ------------------------------------------------------------------

!  SLPHYS
IF (LSLPHY.AND.LSLAG) THEN
  ALLOCATE(YRSLPHY%SAVTEND(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YRSLPHY%NVTEND,YDGEOMETRY%YRDIM%NGPBLKS))
  WRITE(NULOUT,*) 'SAVTEND  ',SIZE(YRSLPHY%SAVTEND  ),SHAPE(YRSLPHY%SAVTEND  )
ELSE
  ! Minimal allocation
  ALLOCATE(YRSLPHY%SAVTEND(1,1,1,1))
ENDIF
YRSLPHY%SAVTEND(:,:,:,:) = 0.0_JPRB


!*       3.    INITIALIZE STANDARD ATMOSPHERE.
!              -------------------------------

!CALL SUSTA(YDGEOMETRY)


!*       4.    INITIALIZE RKRF (RAYLEIGH FRICTION).
!              ------------------------------------
!*       1.18  Rayleigh friction (if active)


!  Defaults still same as 36r1, to be revised later.
!  Suggested values: RRFZ1=68, RRFPLM=100

RRFZ1=61._JPRB
RRFPLM=990._JPRB
RRFTAU=3._JPRB*86400._JPRB

CALL SURAYFRIC(YDGEOMETRY%YRVAB,YDGEOMETRY%YRDIMV,YDDYN)

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUDYN1C',1,ZHOOK_HANDLE)
END SUBROUTINE SUDYN1C
