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

SUBROUTINE SPCHOR(YDGEOMETRY,YDMODEL,KSTA,KEND,&
 & PSPVOR,PSPDIV,PSPT,PSPSPD,PSPSVD,PSPSNHX,PSPGFL,PSPSP,&
 & PSPTNDHD_VOR,PSPTNDHD_DIV,PSPTNDHD_T,&
 & PSPTNDHD_SPD,PSPTNDHD_SVD,PSPTNDHD_SNHX,PSPTNDHD_GFL)

!     ------------------------------------------------------------------

!**** *SPCHOR* - HORIZONTAL SPECTRAL SPACE COMPUTATIONS

!     Purpose.
!     --------

!**   Interface.
!     ----------
!        *CALL* *SPCHOR(..)

!        Explicit arguments :
!        --------------------
!         INPUT:
!          KSTA         : First column processed
!          KEND         : Last columns processed

!         INOUT:
!          PSPVOR       : Vorticity columns
!          PSPDIV       : Divergence columns
!          PSPT         : Temperature columns
!          PSPSPD       : NH pressure departure var.
!          PSPSVD       : NH vertical divergence var.
!          PSPSNHX      : NH "X" part divergence var.
!          PSPGFL       : GFL
!          PSPSP        : Surface Pressure
!          PSPTNDHD_[X] : HD tendencies of variable [X] (for DDH)

!     Method.
!     -------
!        SPCHOR  works in the "normal" spectral domain where a
!                given processor owns NFLEVL levels and NSPEC2
!                spectral coefficients. NSPEC2 depends on MYSETW only,
!                and NFLEVL depends on MYSETV only

!     Externals.
!     ----------
!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      Mats Hamrud and Philippe Courtier  *ECMWF*
!      Original   : 87-11-24 (before 1997 spchor.F was part of spc.F)

! Modifications
! -------------
!   M. Tolstykh 02-02-12: No diffusion on the work sphere if
!                         Fourier diffusion is chosen
!   K. YESSAD 03-04-16: bug correction in call MXTURHD.
!   01-Oct-2003 M. Hamrud           CY28 Cleaning
!   Dec-2003 K. Yessad           cleaning in horizontal diffusion
!   F. Vana, K. Yessad  modify diffusion when LSLHD=T
!   10-Jun-2004 J. Masek            NH cleaning (LFULLIMP)
!   K. Yessad 04-11-15: update the enhanced diffusion
!   F. Vana   05-02-11: split of global LSLHD key + cleaning
!   K. Yessad 08-Feb-2005: horizontal diffusion for NH variables (case c=1)
!   K. Yessad 12-Dec-2005: horizontal diffusion for NH variables (case c>1)
!   K. Yessad 12-Dec-2005: bug corrections for SLHD applied to SVD,NHX.
!   J. Berner 06-09-2006 : added streamfunction perturbations and
!                          computation for dissipations rate for stochstic phys
!   A. Alias  06-03-09: Nudging
!   M. Deque : Nudging in correction mode
!   M. Deque : Nudging with time-variable coefficients
!   P. Marquet : XWNUDG in YOMSNU
!   M. Deque : Nudging with vert-variable coefficients
!   K. Yessad 15-May-2006: memory optimisations for stretched geometry
!   A. Alias 05-11-20  M. Deque : Nudging with wave-variable coefficients
!   K. Yessad (Sep 2008): prune enhanced diffusion (lfrein).
!   F. Vana 13-Jan-2009: recognized case LSLHD_STATIC
!   K. Yessad (Feb 2012): new sponge; tests on CDCONF and NCONF in the caller.
!   T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!   K. Yessad (July 2014): Move some variables.
!   O. Marsden (May 2016): Removed redundant geometry arguments
!   I. Polichtchouk (Nov 2020): Added option to not damp m=0 mode in the sponge
! End Modifications
!------------------------------------------------------------------------------

USE TYPE_MODEL   , ONLY : MODEL
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST       , ONLY : RA
USE YOMCT0       , ONLY : LSPBSBAL
USE YOMCT3       , ONLY : NSTEP
USE YOMJG        , ONLY : JB_STRUCT
USE STOPH_MIX    , ONLY : SPNORMBS,GMTOTDISS
USE YOMMP0       , ONLY : MYSETV
USE YOMNUD       , ONLY : NFNUDG, NTNUDG, LNUDG, LNUDDI, LNUDLP, LNUDSH, LNUDTE,&
 &                        LNUDVO, XNUDDI, XNUDLP,  XNUDSH, XNUDTE, XNUDVO, NSPNU1, NSPNU2
USE YOMSNU       , ONLY : XPNUDG, TNUDTE, TNUDSH, TNUDDI, TNUDVO, TNUDLP, XWNUDG

!------------------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(MODEL)       ,INTENT(INOUT) :: YDMODEL
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTA
INTEGER(KIND=JPIM),INTENT(IN)    :: KEND
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPVOR(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPDIV(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPT(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPSPD(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPSVD(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPSNHX(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPGFL(:,:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPSP(:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPTNDHD_VOR(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPTNDHD_DIV(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPTNDHD_T(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPTNDHD_SPD(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPTNDHD_SVD(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPTNDHD_SNHX(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPTNDHD_GFL(:,:,:)

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZSPX(YDGEOMETRY%YRDIMV%NFLEVL*YDMODEL%YRML_GCONF%YRDIMF%NS3D+1,0:YDGEOMETRY%YRDIM%NSMAX,2)
REAL(KIND=JPRB) :: ZSPY(YDGEOMETRY%YRDIMV%NFLEVL*YDMODEL%YRML_GCONF%YRDIMF%NS3D+1,0:YDGEOMETRY%YRDIM%NSMAX,2)
REAL(KIND=JPRB) :: ZUSNUD2(KSTA:KEND),ZUSNUD1(KSTA:KEND)
REAL(KIND=JPRB) :: ZCOEF(KSTA:KEND),ZCOEFVOR(KSTA:KEND)
REAL(KIND=JPRB),ALLOCATABLE :: ZSPX_HDS(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE :: ZSPY_HDS(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE :: ZGP7(:,:,:)
REAL(KIND=JPRB) :: ZREF(KSTA:KEND)

INTEGER(KIND=JPIM) :: IJSE, ILEV, ILEV0, ILEVG, IL, ILL, IM, IN, INUM_GFL0,&
 & INUM_THE0, IS0, IVTHS, JLEV, JMLOC, JN, JSP, JSTEP, JGFL
INTEGER(KIND=JPIM) :: IOFF_VOR, IOFF_DIV, IOFF_U, IOFF_V, IOFF_NS_T, &
                    & IOFF_EW_U, IOFF_EW_V, IOFF_EW_T, IN_MIN
REAL(KIND=JPRB) :: ZDT, ZLAPDI, ZLAPIN
REAL(KIND=JPRB) :: ZCON1
LOGICAL :: LLDO
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     local variables for stochastic backscatter------------------------
REAL(KIND=JPRB) :: ZSPTOTDISS(YDGEOMETRY%YRDIMV%NFLEVL,KSTA:KEND)  ! spectral coeffs of total dissipation
REAL(KIND=JPRB) :: ZZP(YDGEOMETRY%YRDIMV%NFLEVL,KSTA:KEND)
REAL(KIND=JPRB) :: ZSPVORPERT(YDGEOMETRY%YRDIMV%NFLEVL,KSTA:KEND)   ! spectral vorticity  increment due to backscatter
!REAL(KIND=JPRB) :: ZSPDIVPERT(NFLEVL,KSTA:KEND)   ! spectral divergence increment due to backscatter
REAL(KIND=JPRB) :: ZSPVORT(YDGEOMETRY%YRDIMV%NFLEVL,KSTA:KEND)      ! spectral representation of vorticity forcing pattern
REAL(KIND=JPRB) :: ZSPDIVBALPERT(YDGEOMETRY%YRDIMV%NFLEVL,KSTA:KEND)
REAL(KIND=JPRB) :: ZSPDIVOMEGAPERT(YDGEOMETRY%YRDIMV%NFLEVL,KSTA:KEND)
REAL(KIND=JPRB) :: ZSPTPERT(YDGEOMETRY%YRDIMV%NFLEVL,KSTA:KEND)
REAL(KIND=JPRB) :: ZSPSPPERT(KSTA:KEND)

!     ------------------------------------------------------------------

#include "inv_trans.h"
#include "dir_trans.h"

#include "mxturhd.h"
#include "balvert.intfb.h"





!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SPCHOR',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, &
 &  YDLAP=>YDGEOMETRY%YRLAP, YDDYN=>YDMODEL%YRML_DYN%YRDYN,YDDYNA=>YDMODEL%YRML_DYN%YRDYNA, &
 & YDSTOPH=>YDMODEL%YRML_PHY_STOCH%YRSTOPH, &
 & YDSPNG=>YDMODEL%YRML_DYN%YRSPNG, &
 & YDMDDH=>YDMODEL%YRML_DIAG%YRMDDH,YDRIP=>YDMODEL%YRML_GCONF%YRRIP,YDLDDH=>YDMODEL%YRML_DIAG%YRLDDH, &
 & YGFL=>YDMODEL%YRML_GCONF%YGFL, &
 & YDDIMF=>YDMODEL%YRML_GCONF%YRDIMF)

ASSOCIATE(NUMFLDS=>YGFL%NUMFLDS, YA=>YGFL%YA, YCOMP=>YGFL%YCOMP, YI=>YGFL%YI, &
 & YL=>YGFL%YL, YO3=>YGFL%YO3, YQ=>YGFL%YQ, &
 & NPROMA=>YDDIM%NPROMA, NRESOL=>YDDIM%NRESOL, NSMAX=>YDDIM%NSMAX, NGPTOT=>YDGEM%NGPTOT, &
 & NUMP=>YDDIM%NUMP, &
 & LSPT=>YDDIMF%LSPT, NFTHER=>YDDIMF%NFTHER, NS3D=>YDDIMF%NS3D, &
 & NFLEVG=>YDDIMV%NFLEVG, NFLEVL=>YDDIMV%NFLEVL, &
 & LSTRHD=>YDDYN%LSTRHD, RCORDIT=>YDDYN%RCORDIT, RDHI=>YDDYN%RDHI, &
 & RDHS=>YDDYN%RDHS, RDIDIV=>YDDYN%RDIDIV, RDIGFL=>YDDYN%RDIGFL, &
 & RDIPD=>YDDYN%RDIPD, RDISP=>YDDYN%RDISP, RDITG=>YDDYN%RDITG, &
 & RDIVD=>YDDYN%RDIVD, RDIVOR=>YDDYN%RDIVOR, RDSDIV=>YDDYN%RDSDIV, &
 & RDSVD=>YDDYN%RDSVD, RDSVOR=>YDDYN%RDSVOR, &
 & LRHDDDH=>YDLDDH%LRHDDDH, &
 & MHDDDH_Q=>YDMDDH%MHDDDH_Q, &
 & MYLEVS=>YDMP%MYLEVS, NBSETLEV=>YDMP%NBSETLEV, NPSP=>YDMP%NPSP, &
 & NPTRLL=>YDMP%NPTRLL, &
 & TDT=>YDRIP%TDT, LM0DAMP=>YDDYN%LM0DAMP, &
 & LNSPONGE=>YDSPNG%LNSPONGE, RSPONGF=>YDSPNG%RSPONGF, &
 & GPSTREAM=>YDSTOPH%GPSTREAM, GPTEMP=>YDSTOPH%GPTEMP, &
 & GPTOTDISS=>YDSTOPH%GPTOTDISS, GPTOTDISS_SMOOTH=>YDSTOPH%GPTOTDISS_SMOOTH, &
 & GPVORTGRAD=>YDSTOPH%GPVORTGRAD, LSPBSDISS=>YDSTOPH%LSPBSDISS, &
 & LSPBSNORM=>YDSTOPH%LSPBSNORM, LSTOPH_CASBS=>YDSTOPH%LSTOPH_CASBS, &
 & LSTOPH_SPBS=>YDSTOPH%LSTOPH_SPBS, LSTOPH_SPBS_T=>YDSTOPH%LSTOPH_SPBS_T, &
 & LSTOPH_SPBS_VORT=>YDSTOPH%LSTOPH_SPBS_VORT, LVORTCON=>YDSTOPH%LVORTCON, &
 & NFRSTOPH_SPBS=>YDSTOPH%NFRSTOPH_SPBS, NFRSTOPH_VC=>YDSTOPH%NFRSTOPH_VC, &
 & RSMOOTH=>YDSTOPH%RSMOOTH, SPSTREAM=>YDSTOPH%SPSTREAM, &
 & SPSTREAM_FORC=>YDSTOPH%SPSTREAM_FORC, SPTEMP=>YDSTOPH%SPTEMP, &
 & SPTEMP_FORC=>YDSTOPH%SPTEMP_FORC)
!     ------------------------------------------------------------------

!  ILEV0 is the ordinal of the first level for this processor in the interval 1:NFLEVG
ILEV0=NPTRLL(MYSETV) - 1

!     ------------------------------------------------------------------

!*        2.    MAIN HORIZONTAL DIFFUSION.
!              ---------------------------

! Memory transfert for spectral HD tendencies
IF (LRHDDDH) THEN
  PSPTNDHD_VOR(1:NFLEVG,KSTA:KEND)=-PSPVOR(1:NFLEVG,KSTA:KEND)
  PSPTNDHD_DIV(1:NFLEVG,KSTA:KEND)=-PSPDIV(1:NFLEVG,KSTA:KEND)
  PSPTNDHD_T(1:NFLEVG,KSTA:KEND)=-PSPT(1:NFLEVG,KSTA:KEND)
  IF (YDDYNA%LNHDYN) THEN
    PSPTNDHD_SPD(1:NFLEVG,KSTA:KEND)=-PSPSPD(1:NFLEVG,KSTA:KEND)
    PSPTNDHD_SVD(1:NFLEVG,KSTA:KEND)=-PSPSVD(1:NFLEVG,KSTA:KEND)
    IF (YDDYNA%LNHX) PSPTNDHD_SNHX(1:NFLEVG,KSTA:KEND)=-PSPSNHX(1:NFLEVG,KSTA:KEND)
  ENDIF
  IF (YQ%LSP) THEN
    DO JSP=KSTA,KEND
      DO JLEV=1,NFLEVG
        PSPTNDHD_GFL(JLEV,JSP,MHDDDH_Q)=-PSPGFL(JLEV,JSP,YQ%MPSP)
      ENDDO
    ENDDO
  ENDIF
ENDIF ! LRHDDDH

IF (LSTRHD) THEN

!           2.1  Diffusion when stretching.

!       * The order in ZSPX, ZSPY and also in RDHI is assumed to be the
!         following one:
!         - VOR,DIV.
!         - Thermodynamic variables.
!         - NH variables (SPD,SVD,NHX).
!         - GFL variables.
!         - GMVS variables.
  INUM_GFL0=2+NFTHER

  ZSPX(:,:,:)=0.0_JPRB
  ZSPY(:,:,:)=0.0_JPRB

  DO JMLOC=1,NUMP

    IM=YDLAP%MYMS(JMLOC)
    DO JN=IM,NSMAX
      IJSE=YDLAP%NASM0(IM)+2*(JN-IM)
      INUM_THE0=0
      ! * Vor,Div.
      ZLAPIN=YDLAP%RLAPIN(JN)
      DO JLEV=1,NFLEVL
        ZSPX(JLEV         ,JN,1)=PSPVOR(JLEV,IJSE  )*ZLAPIN
        ZSPX(JLEV         ,JN,2)=PSPVOR(JLEV,IJSE+1)*ZLAPIN
        ZSPX(JLEV+  NFLEVL,JN,1)=PSPDIV(JLEV,IJSE  )*ZLAPIN
        ZSPX(JLEV+  NFLEVL,JN,2)=PSPDIV(JLEV,IJSE+1)*ZLAPIN
      ENDDO
      INUM_THE0=INUM_THE0+2
      ! * Temperature.
      IF (LSPT) THEN
        DO JLEV=1,NFLEVL
          ILEVG=ILEV0+JLEV
          ZSPX(JLEV+INUM_THE0*NFLEVL,JN,1)=PSPT(JLEV,IJSE  )&
           & -RCORDIT(ILEVG)*PSPSP(IJSE  )
          ZSPX(JLEV+INUM_THE0*NFLEVL,JN,2)=PSPT(JLEV,IJSE+1)&
           & -RCORDIT(ILEVG)*PSPSP(IJSE+1)
        ENDDO
        INUM_THE0=INUM_THE0+1
      ENDIF
      ! * Pressure departure.
      IF (YDDYNA%LNHDYN) THEN
        DO JLEV=1,NFLEVL
          ZSPX(JLEV+INUM_THE0*NFLEVL,JN,1)=PSPSPD(JLEV,IJSE  )
          ZSPX(JLEV+INUM_THE0*NFLEVL,JN,2)=PSPSPD(JLEV,IJSE+1)
        ENDDO
        INUM_THE0=INUM_THE0+1
      ENDIF
      ! * Vertical divergence.
      IF (YDDYNA%LNHDYN) THEN
        DO JLEV=1,NFLEVL
          ZSPX(JLEV+INUM_THE0*NFLEVL,JN,1)=PSPSVD(JLEV,IJSE  )
          ZSPX(JLEV+INUM_THE0*NFLEVL,JN,2)=PSPSVD(JLEV,IJSE+1)
        ENDDO
        INUM_THE0=INUM_THE0+1
      ENDIF
      ! * NHX.
      IF (YDDYNA%LNHX) THEN
        DO JLEV=1,NFLEVL
          ZSPX(JLEV+INUM_THE0*NFLEVL,JN,1)=PSPSNHX(JLEV,IJSE  )
          ZSPX(JLEV+INUM_THE0*NFLEVL,JN,2)=PSPSNHX(JLEV,IJSE+1)
        ENDDO
        INUM_THE0=INUM_THE0+1
      ENDIF
      ! * GFL.
      DO JGFL=1,NUMFLDS
        IF (YCOMP(JGFL)%LSP) THEN
          DO JLEV=1,NFLEVL
            ZSPX(JLEV+(INUM_GFL0+YCOMP(JGFL)%MPSP-1)*NFLEVL,JN,1)=&
             & PSPGFL(JLEV,IJSE  ,YCOMP(JGFL)%MPSP)
            ZSPX(JLEV+(INUM_GFL0+YCOMP(JGFL)%MPSP-1)*NFLEVL,JN,2)=&
             & PSPGFL(JLEV,IJSE+1,YCOMP(JGFL)%MPSP)
          ENDDO
        ENDIF
      ENDDO
    ENDDO

    ! * ln(prehyds).
    ILL=NFLEVL*NS3D+1
    DO JN=IM,NSMAX
      IJSE=YDLAP%NASM0(IM)+2*(JN-IM)
      ZSPX(ILL,JN,1)=PSPSP(IJSE  )
      ZSPX(ILL,JN,2)=PSPSP(IJSE+1)
    ENDDO

    IS0=YDLAP%NSE0L(JMLOC)
    CALL MXTURHD(NSMAX+1-IM,ILL,ILL,-2,.TRUE.,&
     & RDHI(1,IS0+1,1),RDHI(1,IS0+1,2),ZSPX(1,IM,1),ZSPY(1,IM,1))
    CALL MXTURHD(NSMAX+1-IM,ILL,ILL,3,.FALSE.,&
     & RDHI(1,IS0+1,1),RDHI(1,IS0+1,3),ZSPX(1,IM,1),ZSPY(1,IM,1))

    IF (IM > 0) THEN
      CALL MXTURHD(NSMAX+1-IM,ILL,ILL,-2,.TRUE.,&
       & RDHI(1,IS0+1,1),RDHI(1,IS0+1,2),ZSPX(1,IM,2),ZSPY(1,IM,2))
      CALL MXTURHD(NSMAX+1-IM,ILL,ILL,3,.FALSE.,&
       & RDHI(1,IS0+1,1),RDHI(1,IS0+1,3),ZSPX(1,IM,2),ZSPY(1,IM,2))
    ENDIF

    DO JN=IM,NSMAX
      IJSE=YDLAP%NASM0(IM)+2*(JN-IM)
      INUM_THE0=0
      ! * Vor,Div.
      ZLAPDI=YDLAP%RLAPDI(JN)
      DO JLEV=1,NFLEVL
        ILEVG=ILEV0+JLEV
        PSPVOR(JLEV,IJSE  )=ZSPY(JLEV       ,JN,1)*ZLAPDI
        PSPVOR(JLEV,IJSE+1)=ZSPY(JLEV       ,JN,2)*ZLAPDI
        PSPDIV(JLEV,IJSE  )=ZSPY(JLEV+NFLEVL,JN,1)*ZLAPDI
        PSPDIV(JLEV,IJSE+1)=ZSPY(JLEV+NFLEVL,JN,2)*ZLAPDI
      ENDDO
      INUM_THE0=INUM_THE0+2
      ! * Temperature.
      IF (LSPT) THEN
        DO JLEV=1,NFLEVL
          ILEVG=ILEV0+JLEV
          PSPT(JLEV,IJSE  )=ZSPY(JLEV+INUM_THE0*NFLEVL,JN,1)&
           & +RCORDIT(ILEVG)*ZSPY(ILL,JN,1)
          PSPT(JLEV,IJSE+1)=ZSPY(JLEV+INUM_THE0*NFLEVL,JN,2)&
           & +RCORDIT(ILEVG)*ZSPY(ILL,JN,2)
        ENDDO
        INUM_THE0=INUM_THE0+1
      ENDIF
      ! * Pressure departure.
      IF (YDDYNA%LNHDYN) THEN
        DO JLEV=1,NFLEVL
          PSPSPD(JLEV,IJSE  )=ZSPY(JLEV+INUM_THE0*NFLEVL,JN,1)
          PSPSPD(JLEV,IJSE+1)=ZSPY(JLEV+INUM_THE0*NFLEVL,JN,2)
        ENDDO
        INUM_THE0=INUM_THE0+1
      ENDIF
      ! * Vertical divergence.
      IF (YDDYNA%LNHDYN) THEN
        DO JLEV=1,NFLEVL
          PSPSVD(JLEV,IJSE  )=ZSPY(JLEV+INUM_THE0*NFLEVL,JN,1)
          PSPSVD(JLEV,IJSE+1)=ZSPY(JLEV+INUM_THE0*NFLEVL,JN,2)
        ENDDO
        INUM_THE0=INUM_THE0+1
      ENDIF
      ! * NHX.
      IF (YDDYNA%LNHX) THEN
        DO JLEV=1,NFLEVL
          PSPSNHX(JLEV,IJSE  )=ZSPY(JLEV+INUM_THE0*NFLEVL,JN,1)
          PSPSNHX(JLEV,IJSE+1)=ZSPY(JLEV+INUM_THE0*NFLEVL,JN,2)
        ENDDO
        INUM_THE0=INUM_THE0+1
      ENDIF
      ! * GFL.
      DO JGFL=1,NUMFLDS
        IF (YCOMP(JGFL)%LSP) THEN
          DO JLEV=1,NFLEVL
            PSPGFL(JLEV,IJSE  ,YCOMP(JGFL)%MPSP)=&
             & ZSPY(JLEV+(INUM_GFL0+YCOMP(JGFL)%MPSP-1)*NFLEVL,JN,1)
            PSPGFL(JLEV,IJSE+1,YCOMP(JGFL)%MPSP)=&
             & ZSPY(JLEV+(INUM_GFL0+YCOMP(JGFL)%MPSP-1)*NFLEVL,JN,2)
          ENDDO
        ENDIF
      ENDDO
    ENDDO

    DO JN=IM,NSMAX
      IJSE=YDLAP%NASM0(IM)+2*(JN-IM)
      PSPSP(IJSE  )=ZSPY(ILL,JN,1)
      PSPSP(IJSE+1)=ZSPY(ILL,JN,2)
    ENDDO

  ENDDO

  IF ((YDDYNA%LSLHD_W.OR.YDDYNA%LSLHD_SVD).AND.(.NOT.YDDYNA%LSLHD_STATIC)) THEN

    ! * Additional diffusion done in the direct integration only.
    !   The order in ZSPX_HDS, ZSPY_HDS and also in RDHS is assumed to be
    !   the following one: VOR,DIV,SVD,NHX.

    ! * First compute the number of prognostic variables on which the HDS
    !   diffusion is applied (currently IVTHS=(IL/NFLEVL) is between 0 and 4).
    !   This piece of code must remain consistent with the equivalent one done
    !   in SUHDU and SUALDYNB.
    IVTHS=0
    IF (YDDYNA%LSLHD_W) IVTHS=IVTHS+2
    IF (YDDYNA%LNHDYN.AND.YDDYNA%LSLHD_SVD) IVTHS=IVTHS+1
    IF (YDDYNA%LNHX.AND.YDDYNA%LSLHD_SVD) IVTHS=IVTHS+1
    IL=IVTHS*NFLEVL

    ! * Allocate ZSPX_HDS,ZSPY_HDS, then fill them with zeros.
    ALLOCATE(ZSPX_HDS(IL,0:NSMAX,2))
    ALLOCATE(ZSPY_HDS(IL,0:NSMAX,2))
    ZSPX_HDS(:,:,:)=0.0_JPRB
    ZSPY_HDS(:,:,:)=0.0_JPRB

    DO JMLOC=1,NUMP

      IM=YDLAP%MYMS(JMLOC)
      DO JN=IM,NSMAX
        IJSE=YDLAP%NASM0(IM)+2*(JN-IM)
        INUM_THE0=0
        ! * Vor,Div.
        IF (YDDYNA%LSLHD_W) THEN
          ZLAPIN=YDLAP%RLAPIN(JN)
          DO JLEV=1,NFLEVL
            ZSPX_HDS(JLEV         ,JN,1)=PSPVOR(JLEV,IJSE  )*ZLAPIN
            ZSPX_HDS(JLEV         ,JN,2)=PSPVOR(JLEV,IJSE+1)*ZLAPIN
            ZSPX_HDS(JLEV+  NFLEVL,JN,1)=PSPDIV(JLEV,IJSE  )*ZLAPIN
            ZSPX_HDS(JLEV+  NFLEVL,JN,2)=PSPDIV(JLEV,IJSE+1)*ZLAPIN
          ENDDO
          INUM_THE0=INUM_THE0+2
        ENDIF
        ! * Vertical divergence.
        IF (YDDYNA%LNHDYN.AND.YDDYNA%LSLHD_SVD) THEN
          DO JLEV=1,NFLEVL
            ZSPX_HDS(JLEV+INUM_THE0*NFLEVL,JN,1)=PSPSVD(JLEV,IJSE  )
            ZSPX_HDS(JLEV+INUM_THE0*NFLEVL,JN,2)=PSPSVD(JLEV,IJSE+1)
          ENDDO
          INUM_THE0=INUM_THE0+1
        ENDIF
        ! * NHX.
        IF (YDDYNA%LNHX.AND.YDDYNA%LSLHD_SVD) THEN
          DO JLEV=1,NFLEVL
            ZSPX_HDS(JLEV+INUM_THE0*NFLEVL,JN,1)=PSPSNHX(JLEV,IJSE  )
            ZSPX_HDS(JLEV+INUM_THE0*NFLEVL,JN,2)=PSPSNHX(JLEV,IJSE+1)
          ENDDO
          INUM_THE0=INUM_THE0+1
        ENDIF
      ENDDO
      ILL=IL
      IS0=YDLAP%NSE0L(JMLOC)
      CALL MXTURHD(NSMAX+1-IM,ILL,ILL,-2,.TRUE.,&
       & RDHS(1,IS0+1,1),RDHS(1,IS0+1,2),&
       & ZSPX_HDS(1,IM,1),ZSPY_HDS(1,IM,1))
      CALL MXTURHD(NSMAX+1-IM,ILL,ILL,3,.FALSE.,&
       & RDHS(1,IS0+1,1),RDHS(1,IS0+1,3),&
       & ZSPX_HDS(1,IM,1),ZSPY_HDS(1,IM,1))

      IF (IM > 0) THEN
        CALL MXTURHD(NSMAX+1-IM,ILL,ILL,-2,.TRUE.,&
         & RDHS(1,IS0+1,1),RDHS(1,IS0+1,2),&
         & ZSPX_HDS(1,IM,2),ZSPY_HDS(1,IM,2))
        CALL MXTURHD(NSMAX+1-IM,ILL,ILL,3,.FALSE.,&
         & RDHS(1,IS0+1,1),RDHS(1,IS0+1,3),&
         & ZSPX_HDS(1,IM,2),ZSPY_HDS(1,IM,2))
      ENDIF

      DO JN=IM,NSMAX
        IJSE=YDLAP%NASM0(IM)+2*(JN-IM)
        INUM_THE0=0
        ! * Vor,Div.
        IF (YDDYNA%LSLHD_W) THEN
          ZLAPDI=YDLAP%RLAPDI(JN)
          DO JLEV=1,NFLEVL
            PSPVOR(JLEV,IJSE  )=ZSPY_HDS(JLEV       ,JN,1)*ZLAPDI
            PSPVOR(JLEV,IJSE+1)=ZSPY_HDS(JLEV       ,JN,2)*ZLAPDI
            PSPDIV(JLEV,IJSE  )=ZSPY_HDS(JLEV+NFLEVL,JN,1)*ZLAPDI
            PSPDIV(JLEV,IJSE+1)=ZSPY_HDS(JLEV+NFLEVL,JN,2)*ZLAPDI
          ENDDO
          INUM_THE0=INUM_THE0+2
        ENDIF
        ! * Vertical divergence.
        IF (YDDYNA%LNHDYN.AND.YDDYNA%LSLHD_SVD) THEN
          DO JLEV=1,NFLEVL
            PSPSVD(JLEV,IJSE  )=ZSPY_HDS(JLEV+INUM_THE0*NFLEVL,JN,1)
            PSPSVD(JLEV,IJSE+1)=ZSPY_HDS(JLEV+INUM_THE0*NFLEVL,JN,2)
          ENDDO
          INUM_THE0=INUM_THE0+1
        ENDIF
        ! * NHX.
        IF (YDDYNA%LNHX.AND.YDDYNA%LSLHD_SVD) THEN
          DO JLEV=1,NFLEVL
            PSPSNHX(JLEV,IJSE  )=ZSPY_HDS(JLEV+INUM_THE0*NFLEVL,JN,1)
            PSPSNHX(JLEV,IJSE+1)=ZSPY_HDS(JLEV+INUM_THE0*NFLEVL,JN,2)
          ENDDO
          INUM_THE0=INUM_THE0+1
        ENDIF
      ENDDO

    ENDDO

    DEALLOCATE(ZSPX_HDS)
    DEALLOCATE(ZSPY_HDS)

  ENDIF

ELSE

!           2.2  Diffusion when no stretching.

  ZDT=ABS(TDT)
  CALL GSTATS(1037,0)
  IF (.NOT.LM0DAMP) THEN  !no damping on m=0 mode
     DO JMLOC=1,NUMP
      IM=YDLAP%MYMS(JMLOC)
      IF (IM>0) THEN
        IN_MIN=IM
        DO JN=IN_MIN,NSMAX
          IJSE=YDLAP%NASM0(IM)+2*(JN-IM)
          INUM_THE0=0
          ! * Vor,Div.
          DO JLEV=1,NFLEVL
              PSPVOR(JLEV,IJSE  )=PSPVOR(JLEV,IJSE  )/(1.0_JPRB+ZDT*RDIVOR(JLEV,JN))
              PSPVOR(JLEV,IJSE+1)=PSPVOR(JLEV,IJSE+1)/(1.0_JPRB+ZDT*RDIVOR(JLEV,JN))
              PSPDIV(JLEV,IJSE  )=PSPDIV(JLEV,IJSE  )/(1.0_JPRB+ZDT*RDIDIV(JLEV,JN))
              PSPDIV(JLEV,IJSE+1)=PSPDIV(JLEV,IJSE+1)/(1.0_JPRB+ZDT*RDIDIV(JLEV,JN))
          ENDDO
          ! * Temperature.
          IF (LSPT) THEN
            DO JLEV=1,NFLEVL
             ILEVG=ILEV0+JLEV
               PSPT(JLEV,IJSE  )= (PSPT  (JLEV,IJSE  )&
               & +RCORDIT(ILEVG)*PSPSP(IJSE  )*ZDT*RDITG(ILEVG,JN))&
               & /(1.0_JPRB+ZDT*RDITG(ILEVG,JN))
               PSPT(JLEV,IJSE+1)= (PSPT  (JLEV,IJSE+1)&
               & +RCORDIT(ILEVG)*PSPSP(IJSE+1)*ZDT*RDITG(ILEVG,JN))&
               & /(1.0_JPRB+ZDT*RDITG(ILEVG,JN))
            ENDDO
          ENDIF
        ENDDO
      ENDIF ! IM>0
    ENDDO
  ELSE ! * sponge on m=0 mode
!$OMP  PARALLEL&
!$OMP& PRIVATE(JLEV,ILEVG,JSP,IN,JGFL)

!$OMP DO SCHEDULE(STATIC)
  DO JSP=KSTA,KEND
    IN=YDLAP%NVALUE(JSP)

    ! * 2.2.1: GMV variables
!OCL NOVREC
    DO JLEV=1,NFLEVL
      PSPVOR(JLEV,JSP)=PSPVOR(JLEV,JSP)/(1.0_JPRB+ZDT*RDIVOR(JLEV,IN))
      PSPDIV(JLEV,JSP)=PSPDIV(JLEV,JSP)/(1.0_JPRB+ZDT*RDIDIV(JLEV,IN))
    ENDDO
    IF (YDDYNA%LSLHD_W.AND.(.NOT.YDDYNA%LSLHD_STATIC)) THEN
      ! * Additional diffusion
!OCL NOVREC
      DO JLEV=1,NFLEVL
        PSPVOR(JLEV,JSP)=PSPVOR(JLEV,JSP)/(1.0_JPRB+ZDT*RDSVOR(JLEV,IN))
        PSPDIV(JLEV,JSP)=PSPDIV(JLEV,JSP)/(1.0_JPRB+ZDT*RDSDIV(JLEV,IN))
      ENDDO
    ENDIF
    IF (LSPT) THEN
!OCL NOVREC
      DO JLEV=1,NFLEVL
        ILEVG=ILEV0+JLEV
        PSPT  (JLEV,JSP)=(PSPT  (JLEV,JSP)&
         & +RCORDIT(ILEVG)*PSPSP(JSP)*ZDT*RDITG(ILEVG,IN))&
         & /(1.0_JPRB+ZDT*RDITG(ILEVG,IN))
      ENDDO
    ENDIF

    IF (YDDYNA%LNHDYN) THEN
!OCL NOVREC
      DO JLEV=1,NFLEVL
        PSPSPD(JLEV,JSP)=PSPSPD(JLEV,JSP)/(1.0_JPRB+ZDT*RDIPD(JLEV,IN))
        PSPSVD(JLEV,JSP)=PSPSVD(JLEV,JSP)/(1.0_JPRB+ZDT*RDIVD(JLEV,IN))
      ENDDO
      IF (YDDYNA%LNHX) THEN
        DO JLEV=1,NFLEVL
          PSPSNHX(JLEV,JSP)=PSPSNHX(JLEV,JSP)/(1.0_JPRB+ZDT*RDIVD(JLEV,IN))
        ENDDO
      ENDIF
      IF (YDDYNA%LSLHD_SVD.AND.(.NOT.YDDYNA%LSLHD_STATIC)) THEN
        ! * Additional diffusion
!OCL NOVREC
        DO JLEV=1,NFLEVL
          PSPSVD(JLEV,JSP)=PSPSVD(JLEV,JSP)/(1.0_JPRB+ZDT*RDSVD(JLEV,IN))
        ENDDO
        IF (YDDYNA%LNHX) THEN
          DO JLEV=1,NFLEVL
            PSPSNHX(JLEV,JSP)=PSPSNHX(JLEV,JSP)/(1.0_JPRB+ZDT*RDSVD(JLEV,IN))
          ENDDO
        ENDIF
      ENDIF
    ENDIF

    ! * 2.2.2: GFL variables
    DO JGFL=1,NUMFLDS
      IF(YCOMP(JGFL)%LSP) THEN
!OCL NOVREC
        DO JLEV=1,NFLEVL
          PSPGFL(JLEV,JSP,YCOMP(JGFL)%MPSP) =&
           & PSPGFL(JLEV,JSP,YCOMP(JGFL)%MPSP)/&
           & (1.0_JPRB+ZDT*RDIGFL(JLEV,IN,YCOMP(JGFL)%MPSP))
        ENDDO
      ENDIF
    ENDDO
  ENDDO
!$OMP END DO

  ! * 2.2.3: GMVS variables
!OCL NOVREC
!$OMP DO SCHEDULE(STATIC)
  DO JSP=KSTA,KEND
    PSPSP(JSP)=PSPSP(JSP)/(1.0_JPRB+ZDT*RDISP(YDLAP%NVALUE(JSP)))
  ENDDO
!$OMP END DO

!$OMP END PARALLEL
  ENDIF ! LM0DAMP
  CALL GSTATS(1037,1)

ENDIF ! LSTRHD

!           2.3  Sponge (= additional diffusion near the model top).

IF (LNSPONGE) THEN

  ! * 2.3.1: ZCOEF=0. for (m,n)=(0,0), 1. otherwise.
  DO JMLOC=1,NUMP
    IM=YDLAP%MYMS(JMLOC)
    DO JN=IM,NSMAX
      IJSE=YDLAP%NASM0(IM)+2*(JN-IM)
      ZCOEF(IJSE)=MAX(MIN(REAL(JN,JPRB)+REAL(IM,JPRB),1._JPRB),0._JPRB)
      ZCOEF(IJSE+1)=MAX(MIN(REAL(JN,JPRB)+REAL(IM,JPRB),1._JPRB),0._JPRB)
    ENDDO
  ENDDO

  ! * 2.3.2: ZCOEFVOR=0. for (m,n)=(0,0),(0,1) and (1,1), 1. otherwise.
  DO JMLOC=1,NUMP
    IM=YDLAP%MYMS(JMLOC)
    IF (IM <= 1) THEN
      DO JN=IM,1
        IJSE=YDLAP%NASM0(IM)+2*(JN-IM)
        ZCOEFVOR(IJSE)=0._JPRB
        ZCOEFVOR(IJSE+1)=0._JPRB
      ENDDO
      DO JN=2,NSMAX
        IJSE=YDLAP%NASM0(IM)+2*(JN-IM)
        ZCOEFVOR(IJSE)=1._JPRB
        ZCOEFVOR(IJSE+1)=1._JPRB
      ENDDO
    ELSE
      DO JN=IM,NSMAX
        IJSE=YDLAP%NASM0(IM)+2*(JN-IM)
        ZCOEFVOR(IJSE)=1._JPRB
        ZCOEFVOR(IJSE+1)=1._JPRB
      ENDDO
    ENDIF
  ENDDO

  ! * 2.3.3: Sponge for GMV.
  DO JSP=KSTA,KEND
!OCL NOVREC
    DO JLEV=1,NFLEVL
      PSPVOR(JLEV,JSP)=((1.0_JPRB-ZCOEFVOR(JSP))+ZCOEFVOR(JSP)*RSPONGF(JLEV,1))*PSPVOR(JLEV,JSP)
      PSPDIV(JLEV,JSP)=((1.0_JPRB-ZCOEF(JSP))+ZCOEF(JSP)*RSPONGF(JLEV,1))*PSPDIV(JLEV,JSP)
    ENDDO
    IF (LSPT) THEN
!OCL NOVREC
      DO JLEV=1,NFLEVL
        PSPT(JLEV,JSP)=((1.0_JPRB-ZCOEF(JSP))+ZCOEF(JSP)*RSPONGF(JLEV,1))*PSPT(JLEV,JSP)
      ENDDO
    ENDIF
    IF (YDDYNA%LNHDYN) THEN
!OCL NOVREC
      DO JLEV=1,NFLEVL
        PSPSPD(JLEV,JSP)=((1.0_JPRB-ZCOEF(JSP))+ZCOEF(JSP)*RSPONGF(JLEV,1))*PSPSPD(JLEV,JSP)
        PSPSVD(JLEV,JSP)=((1.0_JPRB-ZCOEF(JSP))+ZCOEF(JSP)*RSPONGF(JLEV,1))*PSPSVD(JLEV,JSP)
      ENDDO
    ENDIF
    IF (YDDYNA%LNHDYN.AND.YDDYNA%LNHX) THEN
!OCL NOVREC
      DO JLEV=1,NFLEVL
        PSPSNHX(JLEV,JSP)=((1.0_JPRB-ZCOEF(JSP))+ZCOEF(JSP)*RSPONGF(JLEV,1))*PSPSNHX(JLEV,JSP)
      ENDDO
    ENDIF
  ENDDO

  ! * 2.3.3: Sponge for GFL.
  DO JGFL=1,NUMFLDS
    IF(YCOMP(JGFL)%LSP) THEN
      ! sponge is applied to q,ql,qi,qa,O3 only.
      LLDO=( (YCOMP(JGFL)%MPSP == YQ%MPSP) .OR.&
       & (YCOMP(JGFL)%MPSP == YL%MPSP) .OR. (YCOMP(JGFL)%MPSP == YI%MPSP) .OR.&
       & (YCOMP(JGFL)%MPSP == YA%MPSP) .OR. (YCOMP(JGFL)%MPSP == YO3%MPSP) )
      IF (LLDO) THEN
!OCL NOVREC
        DO JSP=KSTA,KEND
          DO JLEV=1,NFLEVL
            PSPGFL(JLEV,JSP,YCOMP(JGFL)%MPSP) =&
             & ((1.0_JPRB-ZCOEF(JSP))+ZCOEF(JSP)*RSPONGF(JLEV,2))&
             & * PSPGFL(JLEV,JSP,YCOMP(JGFL)%MPSP)
          ENDDO
        ENDDO
      ENDIF
    ENDIF
  ENDDO

ENDIF ! LNSPONGE

! Memory transfert for spectral HD tendencies
IF (LRHDDDH) THEN
  PSPTNDHD_VOR(1:NFLEVG,KSTA:KEND)=PSPTNDHD_VOR(1:NFLEVG,KSTA:KEND)&
   & + PSPVOR(1:NFLEVG,KSTA:KEND)
  PSPTNDHD_DIV(1:NFLEVG,KSTA:KEND)=PSPTNDHD_DIV(1:NFLEVG,KSTA:KEND)&
   & + PSPDIV(1:NFLEVG,KSTA:KEND)
  PSPTNDHD_T(1:NFLEVG,KSTA:KEND)=PSPTNDHD_T(1:NFLEVG,KSTA:KEND)&
   & + PSPT(1:NFLEVG,KSTA:KEND)
  IF (YDDYNA%LNHDYN) THEN
    PSPTNDHD_SPD(1:NFLEVG,KSTA:KEND)=PSPTNDHD_SPD(1:NFLEVG,KSTA:KEND)&
     & + PSPSPD(1:NFLEVG,KSTA:KEND)
    PSPTNDHD_SVD(1:NFLEVG,KSTA:KEND)=PSPTNDHD_SVD(1:NFLEVG,KSTA:KEND)&
     & + PSPSVD(1:NFLEVG,KSTA:KEND)
    IF (YDDYNA%LNHX) PSPTNDHD_SNHX(1:NFLEVG,KSTA:KEND)=PSPTNDHD_SNHX(1:NFLEVG,KSTA:KEND)&
     & + PSPSNHX(1:NFLEVG,KSTA:KEND)
  ENDIF
  IF (YQ%LSP) THEN
    DO JSP=KSTA,KEND
      DO JLEV=1,NFLEVG
        PSPTNDHD_GFL(JLEV,JSP,MHDDDH_Q)=PSPTNDHD_GFL(JLEV,JSP,1)&
         & + PSPGFL(JLEV,JSP,YQ%MPSP)
      ENDDO
    ENDDO
  ENDIF
ENDIF ! LRHDDDH

!     ------------------------------------------------------------------

!*        2.B   UPDATE STOCHASTIC BACKSCATTER VORTICITY INCREMENT
!               -------------------------------------------------

! * Update stochastic backscatter vorticity increment every NFSTOPH_SPBS steps

IF ((LSTOPH_SPBS.OR.LSTOPH_CASBS).AND.(MOD(NSTEP,NFRSTOPH_SPBS) == 0)) THEN
  CALL DIR_TRANS(PSPSCALAR=SPSTREAM_FORC,KRESOL=NRESOL,KPROMA=NPROMA,KVSETSC=NBSETLEV,PGP=GPSTREAM)
  ! CALL DIR_TRANS(PSPSCALAR=SPVELPOT_FORC,KRESOL=NRESOL,KPROMA=NPROMA,KVSETSC=NBSETLEV,PGP=GPVELPOT)
  IF (LSTOPH_SPBS_T) CALL DIR_TRANS(PSPSCALAR=SPTEMP_FORC,KRESOL=NRESOL,KPROMA=NPROMA,KVSETSC=NBSETLEV,PGP=GPTEMP)

  IF (LSPBSBAL) THEN
    ZCON1= RA**(-2)
    DO JSP=KSTA,KEND
      IN=YDLAP%NVALUE(JSP)
      DO JLEV=1,NFLEVL
        ZSPVORPERT(JLEV,JSP)   = -SPSTREAM_FORC(JLEV,JSP)*IN*(IN+1)*ZCON1*TDT
        ! ZSPDIVPERT(JLEV,JSP)   = -SPVELPOT_FORC(JLEV,JSP)*IN*(IN+1)*ZCON1*YRRIP%TDT
      ENDDO
    ENDDO

!         Transform background to grid space and compute its horizontal derivatives

    IF (JB_STRUCT%CONFIG%LJB_NONLINEAR_BALANCE .OR. JB_STRUCT%CONFIG%LJB_OMEGA_BALANCE) THEN
      ALLOCATE(ZGP7(NPROMA, 9*NFLEVG, (NGPTOT+NPROMA-1)/NPROMA))
      ZSPDIVOMEGAPERT(:,:) = 0.0_JPRB
      CALL INV_TRANS (PSPVOR=PSPVOR, PSPDIV=ZSPDIVOMEGAPERT, PSPSCALAR=PSPT,&
                    & LDUVDER=.TRUE., LDSCDERS=.TRUE., LDVORGP=.TRUE., LDDIVGP=.TRUE.,&
                    & KRESOL=NRESOL, KPROMA=NPROMA, KVSETUV=NBSETLEV, KVSETSC=NBSETLEV, PGP=ZGP7)
      IOFF_VOR  = 0
      IOFF_DIV  = NFLEVG
      IOFF_U    = 2*NFLEVG
      IOFF_V    = 3*NFLEVG

      IOFF_NS_T = 5*NFLEVG
      IOFF_EW_U = 6*NFLEVG
      IOFF_EW_V = 7*NFLEVG
      IOFF_EW_T = 8*NFLEVG
    ENDIF

    IF (JB_STRUCT%CONFIG%LJB_NONLINEAR_BALANCE) THEN
      !                 In         In     out
      CALL BALNONLINTL(YDGEOMETRY,ZSPVORPERT,ZZP,ZGP7,&
                     & IOFF_VOR, IOFF_DIV, IOFF_U, IOFF_V, IOFF_EW_U, IOFF_EW_V)
      !            In                  out                  out
      CALL BALVERT(YDGEOMETRY,YDMODEL%YRML_GCONF%YRDIMF%NS3D,JB_STRUCT,ZZP,LDJBBAL=.FALSE.,PSPDIV=ZSPDIVBALPERT,PSPT=ZSPTPERT,&
      !            out
       & PSPSP=ZSPSPPERT)
    ENDIF
    IF (JB_STRUCT%CONFIG%LJB_OMEGA_BALANCE) THEN
      IF (.NOT.JB_STRUCT%CONFIG%LJB_NONLINEAR_BALANCE) ZSPTPERT(:,:)=0.0_JPRB
      !                In         In       In     In   Out
      CALL BALOMEGATL(YDGEOMETRY,ZSPVORPERT,ZSPTPERT,ZSPDIVOMEGAPERT,ZGP7,&
                    & IOFF_VOR, IOFF_DIV, IOFF_U, IOFF_V, IOFF_NS_T,&
                    & IOFF_EW_U, IOFF_EW_V, IOFF_EW_T)
    ENDIF

    IF (JB_STRUCT%CONFIG%LJB_NONLINEAR_BALANCE) THEN
      PSPDIV(:,:) = PSPDIV(:,:) + ZSPDIVBALPERT(:,:)
      PSPT(:,:)   = PSPT(:,:)   + ZSPTPERT(:,:)
      PSPSP(:)    = PSPSP(:)    + ZSPSPPERT(:)
    ENDIF
    IF (JB_STRUCT%CONFIG%LJB_OMEGA_BALANCE) THEN
      PSPDIV(:,:) = PSPDIV(:,:) + ZSPDIVOMEGAPERT(:,:)
    ENDIF
    IF (JB_STRUCT%CONFIG%LJB_NONLINEAR_BALANCE .OR. JB_STRUCT%CONFIG%LJB_OMEGA_BALANCE) THEN
      DEALLOCATE(ZGP7)
    ENDIF
  ENDIF
ENDIF

! * Apply spectral stochastic vorticity increment at every step.

IF (LSTOPH_SPBS.OR.LSTOPH_CASBS) THEN
  ZCON1= RA**(-2)
  DO JSP=KSTA,KEND
    IN=YDLAP%NVALUE(JSP)
    DO JLEV=1,NFLEVL
      ! convert streamfunction and velocity potential forcing to vorticity and divergence forcing
      IF (LSTOPH_SPBS_VORT) THEN
        ZSPVORPERT(JLEV,JSP)   =  SPSTREAM_FORC(JLEV,JSP)*TDT  ! SPSTREAM holds vorticity forcing
      ELSE
        ZSPVORPERT(JLEV,JSP)   = -SPSTREAM_FORC(JLEV,JSP)*IN*(IN+1)*ZCON1*TDT
      ENDIF
      ! ZSPDIVPERT(JLEV,JSP)   = -SPVELPOT_FORC(JLEV,JSP)*IN*(IN+1)*ZCON1*YRRIP%TDT
    ENDDO
  ENDDO
  IF (LSPBSNORM) THEN
      CALL SPNORMBS(YDGEOMETRY,YDMODEL%YRML_PHY_STOCH%YRSTOPH,YDRIP,PSPVOR,PSPDIV*0._JPRB,ZSPVORPERT,ZSPVORPERT*0._JPRB,.FALSE.)
  ENDIF

  PSPVOR(:,:) = PSPVOR(:,:) + ZSPVORPERT(:,:)   !  add spectral stochastic vorticity increment
  ! PSPDIV(:,:) = PSPDIV(:,:) + ZSPDIVPERT(:,:)   !  add spectral stochastic divergence increment
  ! IF (LSTOPH_SPBS_T) PSPT(:,:) = PSPT(:,:) + SPTEMP_FORC(:,:)*YRRIP%TDT    ! 2nd attempt scheme

  IF (LSTOPH_SPBS_T) THEN
    ZCON1= RA**(-1)
    DO JSP=KSTA,KEND
      IN=YDLAP%NVALUE(JSP)
      DO JLEV=1,NFLEVL
        PSPT(JLEV,JSP)=  PSPT(JLEV,JSP) + SPTEMP_FORC(JLEV,JSP)*IN*ZCON1*TDT
      ENDDO
    ENDDO
  ENDIF

ENDIF

!     ------------------------------------------------------------------

!*       3     NUDGING.
!              --------

IF (LNUDG) THEN

  !     3.1 GMV VARIABLES:

  ! * Vorticity:
  IF(LNUDVO)THEN
    IF(XNUDVO <= 0) THEN
      DO JSP=KSTA,KEND
        ZUSNUD2(JSP)=XNUDVO
        ZUSNUD1(JSP)=1.0_JPRB
      ENDDO
    ELSE
      DO JSP=KSTA,KEND
        ZUSNUD2(JSP)=XNUDVO*XWNUDG*&
         &MIN(1.0_JPRB,MAX(0.0_JPRB,1.0_JPRB*(NSPNU2-YDLAP%NVALUE(JSP))/REAL(NSPNU2-NSPNU1)))
        ZUSNUD1(JSP)=1.0_JPRB/(1.0_JPRB+ZUSNUD2(JSP))
      ENDDO
    ENDIF
    DO JLEV=1,NFLEVL
      ILEV=MYLEVS(JLEV)
      IF (ILEV >= NTNUDG) THEN
        ZREF(:)=0.0_JPRB
        DO JSTEP=1,NFNUDG
          DO JSP=KSTA,KEND
            ZREF(JSP)=ZREF(JSP)+TNUDVO(JSP,JLEV,JSTEP)*XPNUDG(JSTEP)
          ENDDO
        ENDDO

        DO JSP=KSTA,KEND
          PSPVOR(JLEV,JSP)=(PSPVOR(JLEV,JSP)+ZUSNUD2(JSP)*ZREF(JSP))*ZUSNUD1(JSP)
        ENDDO
      ENDIF
    ENDDO
  ENDIF

  ! * Divergence:
  IF(LNUDDI)THEN
    IF(XNUDDI <= 0) THEN
      DO JSP=KSTA,KEND
        ZUSNUD2(JSP)=XNUDDI
        ZUSNUD1(JSP)=1.0_JPRB
      ENDDO
    ELSE
      DO JSP=KSTA,KEND
        ZUSNUD2(JSP)=XNUDDI*XWNUDG*&
         &MIN(1.0_JPRB,MAX(0.0_JPRB,1.0_JPRB*(NSPNU2-YDLAP%NVALUE(JSP))/REAL(NSPNU2-NSPNU1)))
        ZUSNUD1(JSP)=1.0_JPRB/(1.0_JPRB+ZUSNUD2(JSP))
      ENDDO
    ENDIF
    DO JLEV=1,NFLEVL
      ILEV=MYLEVS(JLEV)
      IF (ILEV >= NTNUDG) THEN
        ZREF(:)=0.0_JPRB
        DO JSTEP=1,NFNUDG
          DO JSP=KSTA,KEND
            ZREF(JSP)=ZREF(JSP)+TNUDDI(JSP,JLEV,JSTEP)*XPNUDG(JSTEP)
          ENDDO
        ENDDO

        DO JSP=KSTA,KEND
          PSPDIV(JLEV,JSP)=(PSPDIV(JLEV,JSP)+ZUSNUD2(JSP)*ZREF(JSP))*ZUSNUD1(JSP)
        ENDDO
      ENDIF
    ENDDO
  ENDIF

  ! * Temperature:
  IF(LSPT.AND.LNUDTE)THEN
    IF(XNUDTE <= 0) THEN
      DO JSP=KSTA,KEND
        ZUSNUD2(JSP)=XNUDTE
        ZUSNUD1(JSP)=1.0_JPRB
      ENDDO
    ELSE
      DO JSP=KSTA,KEND
        ZUSNUD2(JSP)=XNUDTE*XWNUDG*&
         &MIN(1.0_JPRB,MAX(0.0_JPRB,1.0_JPRB*(NSPNU2-YDLAP%NVALUE(JSP))/REAL(NSPNU2-NSPNU1)))
        ZUSNUD1(JSP)=1.0_JPRB/(1.0_JPRB+ZUSNUD2(JSP))
      ENDDO
    ENDIF
    DO JLEV=1,NFLEVL
      ILEV=MYLEVS(JLEV)
      IF (ILEV >= NTNUDG) THEN
        ZREF(:)=0.0_JPRB
        DO JSTEP=1,NFNUDG
          DO JSP=KSTA,KEND
            ZREF(JSP)=ZREF(JSP)+TNUDTE(JSP,JLEV,JSTEP)*XPNUDG(JSTEP)
          ENDDO
        ENDDO

        DO JSP=KSTA,KEND
          PSPT(JLEV,JSP)=(PSPT(JLEV,JSP)+ZUSNUD2(JSP)*ZREF(JSP))*ZUSNUD1(JSP)
        ENDDO
      ENDIF
    ENDDO
  ENDIF

  ! * Pressure departure variable:
  !   Currently not yet coded, no nudging is done on this variable.

  ! * Vertical divergence variable:
  !   Currently not yet coded, no nudging is done on this variable.


  !     3.2 GFL VARIABLES

  ! * In a final stage, this part of code must be valid for all
  !   GFL variables with a loop on JGFL; for that it needs to
  !   update the setup, i.e. define quantities TNUDGFL, LNUDGFL, XNUDGFL
  !   replacing the current (TNUDSH,TNUDSV), (LNUDSH,LNUDSV),
  !   (XNUDSH,XNUDSV).
  !   The provisional way of updating the code which has been adopted
  !   is to use LNUDSH,XNUDSH for all GFL variables.
  !   Recoding the setup has not yet been done for the moment.

  DO JGFL=1,NUMFLDS
    IF (YCOMP(JGFL)%LSP) THEN
      IF(LNUDSH)THEN
        IF(XNUDSH <= 0) THEN
          DO JSP=KSTA,KEND
            ZUSNUD2(JSP)=XNUDSH
            ZUSNUD1(JSP)=1.0_JPRB
          ENDDO
        ELSE
          DO JSP=KSTA,KEND
            ZUSNUD2(JSP)=XNUDSH*XWNUDG*&
             &MIN(1.0_JPRB,MAX(0.0_JPRB,1.0_JPRB*(NSPNU2-YDLAP%NVALUE(JSP))/REAL(NSPNU2-NSPNU1)))
            ZUSNUD1(JSP)=1.0_JPRB/(1.0_JPRB+ZUSNUD2(JSP))
          ENDDO
        ENDIF
        DO JLEV=1,NFLEVL
          ILEV=MYLEVS(JLEV)
          IF (ILEV >= NTNUDG) THEN
            ZREF(:)=0.0_JPRB
            DO JSTEP=1,NFNUDG
              DO JSP=KSTA,KEND
                ZREF(JSP)=ZREF(JSP)+TNUDSH(JSP,JLEV,JSTEP)*XPNUDG(JSTEP)
              ENDDO
            ENDDO
            DO JSP=KSTA,KEND
              PSPGFL(JLEV,JSP,YQ%MPSP)=(PSPGFL(JLEV,JSP,YQ%MPSP)+&
               & ZUSNUD2(JSP)*ZREF(JSP))*ZUSNUD1(JSP)
            ENDDO
          ENDIF
        ENDDO
      ENDIF
    ENDIF
  ENDDO


  !     3.3 GMVS VARIABLES

  ! * Logarithm of hydrostatic surface pressure:
  IF(LNUDLP.AND.NPSP == 1)THEN
    IF(XNUDLP <= 0) THEN
      DO JSP=KSTA,KEND
        ZUSNUD2(JSP)=XNUDLP
        ZUSNUD1(JSP)=1.0_JPRB
      ENDDO
    ELSE
      DO JSP=KSTA,KEND
        ZUSNUD2(JSP)=XNUDLP*XWNUDG*&
         &MIN(1.0_JPRB,MAX(0.0_JPRB,1.0_JPRB*(NSPNU2-YDLAP%NVALUE(JSP))/REAL(NSPNU2-NSPNU1)))
        ZUSNUD1(JSP)=1.0_JPRB/(1.0_JPRB+ZUSNUD2(JSP))
      ENDDO
    ENDIF
    ZREF(:)=0.0_JPRB
    DO JSTEP=1,NFNUDG
      DO JSP=KSTA,KEND
        ZREF(JSP)=ZREF(JSP)+TNUDLP(JSP,JSTEP)*XPNUDG(JSTEP)
      ENDDO
    ENDDO
    DO JSP=KSTA,KEND
      PSPSP(JSP)=(PSPSP(JSP)+ZUSNUD2(JSP)*ZREF(JSP))*ZUSNUD1(JSP)
    ENDDO
  ENDIF

ENDIF ! LNUDG

!     ------------------------------------------------------------------

!*        4.    INVERSE SPECTRAL TRANSFORMS REQUIRED FOR STOCHASTIC BACKSCATTER
!               ---------------------------------------------------------------

! k.y.: these inverse spectral transforms should be moved under
!       STEPO -> TRANSINVH -> (new) TRANSINV_STOPH in the future,
!       SPCHOR is not the right place to put spectral transforms!

IF ( ((LSTOPH_SPBS.OR.LSTOPH_CASBS).AND.(MOD(NSTEP+2,NFRSTOPH_SPBS) == 0))&
  & .OR. (LVORTCON .AND. MOD(NSTEP,NFRSTOPH_VC) == 0) ) THEN
  ! compute horizontal vorticity gradients for +ve definite numerical diffusion scheme
  CALL INV_TRANS(PSPSCALAR=PSPVOR,KRESOL=NRESOL,KPROMA=NPROMA,LDSCDERS=.TRUE.,&
   &KVSETSC=NBSETLEV,PGP=GPVORTGRAD)
ENDIF

! Transform spectral 'pattern' fields to gridpoint fields
IF (LSTOPH_SPBS.AND.(MOD(NSTEP+1,NFRSTOPH_SPBS) == 0)) THEN
  ! WRITE(NULOUT,*) 'spchor: convert psi forcing pattern to gridpoint space  at NSTEP=',NSTEP
  IF (LSTOPH_SPBS_VORT) THEN
  ! make spectral pattern a vorticity field
    ZCON1= RA**(-2)
    DO JSP=KSTA,KEND
      IN=YDLAP%NVALUE(JSP)
      DO JLEV=1,NFLEVL
        ZSPVORT(JLEV,JSP)   = -SPSTREAM(JLEV,JSP)*IN*(IN+1)*ZCON1
      ENDDO
    ENDDO

    CALL INV_TRANS(PSPSCALAR=ZSPVORT,KRESOL=NRESOL,KPROMA=NPROMA,KVSETSC=NBSETLEV,PGP=GPSTREAM)
  ELSE
    CALL INV_TRANS(PSPSCALAR=SPSTREAM,KRESOL=NRESOL,KPROMA=NPROMA,KVSETSC=NBSETLEV,PGP=GPSTREAM)
  ENDIF
  ! CALL INV_TRANS(PSPSCALAR=SPVELPOT,KRESOL=NRESOL,KPROMA=NPROMA,KVSETSC=NBSETLEV,PGP=GPVELPOT)
  IF (LSTOPH_SPBS_T) CALL INV_TRANS(PSPSCALAR=SPTEMP,KRESOL=NRESOL,KPROMA=NPROMA,KVSETSC=NBSETLEV,PGP=GPTEMP)
ENDIF

IF ((LSTOPH_SPBS.OR.LSTOPH_CASBS).AND.(MOD(NSTEP+1,NFRSTOPH_SPBS) == 0)) THEN

  ! convert unsmoothed dissipation rate to spectral space

  CALL DIR_TRANS(PSPSCALAR=ZSPTOTDISS,KRESOL=NRESOL,KPROMA=NPROMA,KVSETSC=NBSETLEV,PGP=GPTOTDISS)

  ! spectral  smoothing of the  total dissipation rate field

  DO JSP=KSTA,KEND
    IN=YDLAP%NVALUE(JSP)
    DO JLEV=1,NFLEVL
      ZSPTOTDISS(JLEV,JSP)= ZSPTOTDISS(JLEV,JSP)*RSMOOTH(IN)
    ENDDO
  ENDDO

  !calc global mean diss and write to standard out
  IF (LSPBSDISS) THEN
    CALL GMTOTDISS(YDGEOMETRY,YDMODEL%YRML_PHY_STOCH%YRSTOPH,YDRIP,ZSPTOTDISS,.FALSE.)
  ENDIF

  ! transform smoothed total dissipation field back to gridpoint space (for use in callpar.F90)

  CALL INV_TRANS(PSPSCALAR=ZSPTOTDISS,KRESOL=NRESOL,KPROMA=NPROMA,KVSETSC=NBSETLEV,PGP=GPTOTDISS_SMOOTH)
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SPCHOR',1,ZHOOK_HANDLE)
END SUBROUTINE SPCHOR
