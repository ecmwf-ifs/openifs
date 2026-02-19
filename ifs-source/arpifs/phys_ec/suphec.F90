! (C) Copyright 1995- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction
SUBROUTINE SUPHEC(YDGEOMETRY,YDMODEL,KULOUT)

!**** *SUPHEC - INITIALISES PHYSICAL CONSTANTS OF UNCERTAIN VALUE.
!               WITHIN THE E.C.M.W.F. PHYSICS PACKAGE

!     PURPOSE.
!     --------

!          THIS ROUTINE SETS THE VALUES FOR THE PHYSICAL CONSTANTS USED
!     IN THE PARAMETERIZATION ROUTINES WHENEVER THESE VALUES ARE NOT
!     KNOWN WELL ENOUGH TO FORBID ANY TUNING OR WHENEVER THEY ARE
!     SUBJECT TO AN ARBITRARY CHOICE OF THE MODELLER. THESE CONSTANTS
!     ARE DISTRIBUTED IN COMMON DECKS *YOEXXXX* WHERE XXXX CORRESPONDS
!     TO THE INDIVIDUAL PHYSICAL PARAMETRIZATION

!**   INTERFACE.
!     ----------

!          *SUPHEC* IS CALLED FROM *SUPHY*

!     METHOD.
!     -------

!          NONE.

!     EXTERNALS.
!     ----------

!          *SUECRAD*, *SUCUMF*, *SUCUMF2*,*SUVDFS*, *SUSURF*
!          *SUGWD*, *SUGWWMS*, *SUCLD*, *SUCOND*, *SUPHLI*, *SUMETHOX*
!          *SU_AERW* *SUOPTSURF*

!     REFERENCE.
!     ----------

!          SEE PHYSICAL ROUTINES FOR AN EXACT DEFINITION OF THE CONSTANTS.

!     AUTHOR.
!     -------
!      J.-J. MORCRETTE  E.C.M.W.F.    91/06/15  ADAPTATION TO I.F.S.

!     MODIFICATIONS
!     -------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      P.Viterbo     24-May-2004 surf library
!      P.Viterbo     03-Dec-2004 Include user-defined RTHRFRTI
!      M.Ko"hler     03-Dec-2004 cp,moist=cp,dry
!      P.Viterbo     10-Jun-2005 Externalise surf
!      R. El Khatib & J-F Estrade  20-Jan-2005 Default PRSUN for FMR15
!      D.Salmond     22-Nov-2005 Mods for coarser/finer physics
!      P. Lopez      21-Aug-2006 Added call to SUCUMF2 (new linearized convec)
!      JJMorcrette   20060525    MODIS albedo
!      G. Balsamo    20070115    New Hydrology (Van Genuchten + Orog. runoff)
!      G. Balsamo    20080407    Lake model (FLAKE)
!      P. Bechtold   04-Oct-2008 Added call to SUGWWMS
!      G. Balsamo    20081013    Switch for snow liquid water
!      S. Boussetta/G.Balsamo      May 2009   (Add switch for variable LAI: LELAIV)
!      R. El Khatib  05-Aug-2009 protect call to su_aerw
!      R. Forbes     Oct 2010   Set up constants for 550nm optical depth diag
!      S. Boussetta  20101112   Add switch for CTESSEL interactive LAI: RLAIINT
!      P. Bechtold   11-Mar-2011 Changed wrong definition not value of RHOH2O
!      G. Balsamo    June 2011  Add swith for LEAGS (modularity of CO2&Evap)
!      P. Bechtold   26-03-2012 Add RCORIOI RPLRG for small planet
!      J. Hague      21-08-2013 Remove YSP_SBD and Compute ICCS
!      K. Yessad (July 2014): Move some variables.
!      F. Vana & M. Kharoutdinov 09-Feb-2015: Consistency check for CRM sub-model top.
!      R. El Khatib 27-Apr-2015 Systematic call to the surface
!      R. El Khatib 09-Sep-2016 unique code call to the surface
!      E. Dutra 11-Oct-2016 : Add NCSNEC to SUSURF CALL (preparation for snow multi-layer)
!      R. Hogan 24-Jan-2019: Removed Cycle 15 radiation scheme
!      R. Hogan      07-03-2019 Add NALBEDOSCHEME, NEMISSSCHEME, PALFMINPSN to SUSURF call
!      A. Agusti-Panareda 09-04-2021 Include LEAIRCO2COUP switch (CAMSCO2&GPP modularity)
!      A. Agusti-Panareda 06-07-2021 Include LEFARQUHAR switch (CAMSCO2&GPP modularity)
!      R. Forbes     15-11-2020 Pass YDPHY2 into SUCLDP
!     ------------------------------------------------------------------

USE TYPE_MODEL   , ONLY : MODEL
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOELW        , ONLY : NSIL, TSTAND, XP
USE YOMCST       , ONLY : RD, RV, RCPD, RLVTT, RLSTT, RLMLT, RTT
USE YOETHF       , ONLY : R2ES, R3LES, R3IES, R4LES, R4IES, R5LES, R5IES,&
 &                        RVTMP2, RHOH2O, R5ALVCP, R5ALSCP, RALVDCP, RALSDCP, RALFDCP,&
 &                        RTWAT, RTBER, RTBERCU, RTICE, RTICECU, RTWAT_RTICE_R, RTWAT_RTICECU_R,&
 &                        RKOOP1, RKOOP2 
USE YOMMP0       , ONLY : LSCMEC
USE YOMCT0       , ONLY : LROUGH, REXTZ0M, REXTZ0H, NFPOS
USE YOMDYNCORE   , ONLY : RCORIOI, RPLRG
USE YOMVERT      , ONLY : VP00
USE YOEOPTSURF, ONLY : RVR0VT, RVCMAX25, RHUMREL, RA1, RB1, RG0, RGM25, RE_VCMAX, RE_JMAX


!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(MODEL)       ,INTENT(INOUT) :: YDMODEL
INTEGER(KIND=JPIM),INTENT(IN)    :: KULOUT 
!     ------------------------------------------------------------------

#include "susurf.h"

#include "gphpre.intfb.h"
#include "sucld.intfb.h"
#include "sucldp.intfb.h"
#include "suclopn.intfb.h"
#include "su_clop550.intfb.h"
#include "sucond.intfb.h"
#include "sucumf.intfb.h"
#include "sucumf2.intfb.h"
#include "suecrad.intfb.h"
#include "sugwd.intfb.h"
#include "sugwwms.intfb.h"
#include "sumethox.intfb.h"
#include "suphli.intfb.h"
#include "suvdf.intfb.h"
#include "suvdfs.intfb.h"
#include "suwcou.intfb.h"
#include "suwam.intfb.h"
#include "su_aerw.intfb.h"
#include "suoptsurf.intfb.h"
!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZPRES(0:YDGEOMETRY%YRDIMV%NFLEVG),ZPRESF(YDGEOMETRY%YRDIMV%NFLEVG),&
 & ZETA(YDGEOMETRY%YRDIMV%NFLEVG),ZETAH(0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB), ALLOCATABLE :: ZSUN(:)

INTEGER(KIND=JPIM) :: JK,IGRID, ISMAX
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUPHEC',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV, YDMP=>YDGEOMETRY%YRMP, &
 & YDSTA=>YDGEOMETRY%YRSTA, YDVAB=>YDGEOMETRY%YRVAB, YDVETA=>YDGEOMETRY%YRVETA, &
 & YDVFE=>YDGEOMETRY%YRVFE, YDLAP=>YDGEOMETRY%YRLAP, YDCSGLEG=>YDGEOMETRY%YRCSGLEG, &
 & YDCSGEOM=>YDGEOMETRY%YRCSGEOM, YDCSGEOM_NB=>YDGEOMETRY%YRCSGEOM_NB,  &
 & YDCVER=>YDGEOMETRY%YRCVER,YDSPGEOM=>YDGEOMETRY%YSPGEOM, &
 & YDPHY=>YDMODEL%YRML_PHY_MF%YRPHY, YDPHY2=>YDMODEL%YRML_PHY_MF%YRPHY2, &
 & YDMCC=>YDMODEL%YRML_AOC%YRMCC,YDDPHY=>YDMODEL%YRML_PHY_G%YRDPHY, &
 & YDERAD=>YDMODEL%YRML_PHY_RAD%YRERAD, YDESWRT=>YDMODEL%YRML_PHY_RAD%YRESWRT, &
 & YDEPHY=>YDMODEL%YRML_PHY_EC%YREPHY, &
 & LWCOU => YDMODEL%YREWCOU%LWCOU, LWCOU2W => YDMODEL%YREWCOU%LWCOU2W, LWCOUHMF => YDMODEL%YREWCOU%LWCOUHMF, &
 & YGFL=>YDMODEL%YRML_GCONF%YGFL)


ASSOCIATE(NDGLG=>YDDIM%NDGLG, NGPBLKS=>YDDIM%NGPBLKS, NPROMA=>YDDIM%NPROMA, &
 & NSMAX=>YDDIM%NSMAX, &
 & NFLEVG=>YDDIMV%NFLEVG, &
 & LTPROF=>YDDPHY%LTPROF, NCOM=>YDDPHY%NCOM, NCSS=>YDDPHY%NCSS, &
 & NTILES=>YDDPHY%NTILES, NCSNEC=>YDDPHY%NCSNEC,&
 & LEAGS=>YDEPHY%LEAGS, LECTESSEL=>YDEPHY%LECTESSEL, LEFLAKE=>YDEPHY%LEFLAKE, &
 & LEAIRCO2COUP=>YDEPHY%LEAIRCO2COUP, LEFARQUHAR=>YDEPHY%LEFARQUHAR, &
 & LELAIV=>YDEPHY%LELAIV, LEOCCO=>YDEPHY%LEOCCO, LEOCLA=>YDEPHY%LEOCLA, &
 & LEOCML=>YDEPHY%LEOCML, LEOCSA=>YDEPHY%LEOCSA, LEOCWA=>YDEPHY%LEOCWA, &
 & LEPHYS=>YDEPHY%LEPHYS, LESN09=>YDEPHY%LESN09, LESSRO=>YDEPHY%LESSRO, LESNML=>YDEPHY%LESNML, &
 & LEURBAN=>YDEPHY%LEURBAN, LEVGEN=>YDEPHY%LEVGEN, LOCMLTKE=>YDEPHY%LOCMLTKE, RLAIINT=>YDEPHY%RLAIINT, &
 & RTHRFRTI=>YDEPHY%RTHRFRTI, YSURF=>YDEPHY%YSURF, &
 & LCCNL=>YDERAD%LCCNL, LCCNO=>YDERAD%LCCNO, NALBEDOSCHEME=>YDEPHY%NALBEDOSCHEME, &
 & NEMISSSCHEME=>YDEPHY%NEMISSSCHEME, NLWEMISS=>YDERAD%NLWEMISS, &
 & RALFMINPSN=>YDEPHY%RALFMINPSN,LUSEPRE2017RAD=>YDERAD%LUSEPRE2017RAD, &
 & NSW=>YDERAD%NSW, NTSW=>YDERAD%NTSW, NAERMACC=>YDERAD%NAERMACC, RCCNLND=>YDERAD%RCCNLND, &
 & RCCNSEA=>YDERAD%RCCNSEA, &
 & RCIMIN=>YDEPHY%RCIMIN,LMCCDYNSEAICE=>YDMCC%LMCCDYNSEAICE,&
 & LNEMO1WAY=>YDMCC%LNEMO1WAY, &
 & LRAY=>YDPHY%LRAY, LRAYFM15=>YDPHY%LRAYFM15, &
 & NACTAERO=>YGFL%NACTAERO, NSNMLWS=>YDEPHY%NSNMLWS)
!     ------------------------------------------------------------------

!*         0.2    DEFINING DERIVED CONSTANTS FROM UNIVERSAL CONSTANTS
!                 ---------------------------------------------------

CALL GSTATS(1811,0)
!RVTMP2=RCPV/RCPD-1.0_JPRB   !use cp,moist
RVTMP2=0.0_JPRB              !neglect cp,moist
RHOH2O=1000.0_JPRB
R2ES=611.21_JPRB*RD/RV
R3LES=17.502_JPRB
R3IES=22.587_JPRB
R4LES=32.19_JPRB
R4IES=-0.7_JPRB
R5LES=R3LES*(RTT-R4LES)
R5IES=R3IES*(RTT-R4IES)
R5ALVCP=R5LES*RLVTT/RCPD
R5ALSCP=R5IES*RLSTT/RCPD
RALVDCP=RLVTT/RCPD
RALSDCP=RLSTT/RCPD
RALFDCP=RLMLT/RCPD
RTWAT=RTT
RTBER=RTT-5._JPRB
RTBERCU=RTT-5.0_JPRB
RTICE=RTT-23._JPRB
RTICECU=RTT-23._JPRB

RTWAT_RTICE_R=1.0_JPRB/(RTWAT-RTICE)
RTWAT_RTICECU_R=1.0_JPRB/(RTWAT-RTICECU)
ISMAX=NSMAX
IGRID=NDGLG

RKOOP1=2.583_JPRB
RKOOP2=0.48116E-2_JPRB

!     ------------------------------------------------------------------
!*         0.5    DEFINE STANDARD ATMOSPHERE VERTICAL CONFIGURATION
!                 -------------------------------------------------

ZPRES(NFLEVG)=VP00

CALL GPHPRE(1,NFLEVG,1,1,YDVAB,YDCVER,ZPRES,PRESF=ZPRESF)

DO JK=0,NFLEVG
  ZETAH(JK)= ZPRES(JK)/ZPRES(NFLEVG)
ENDDO
DO JK=1,NFLEVG
  ZETA(JK)= ZPRESF(JK)/ZPRES(NFLEVG)
ENDDO

!     ------------------------------------------------------------------

!*         1.     SETTING CONSTANTS FOR DIAGNOSTIC CLOUD SCHEME
!                 ---------------------------------------------

CALL SUCLD(YDMODEL%YRML_PHY_EC%YRECLD,NFLEVG,ZETA)

!     ------------------------------------------------------------------

!*         2.     SETTING CONSTANTS FOR LARGE-SCALE CONDENSATION SCHEME
!                 -----------------------------------------------------

CALL SUCOND(YDMODEL%YRML_PHY_EC%YRECND,KULOUT)

!     ------------------------------------------------------------------

!*         3.     SETTING CONSTANTS FOR CONVECTION SCHEME
!                 ---------------------------------------

CALL SUCUMF(YDSTA,YDDIMV,YDMODEL%YRML_PHY_EC%YRECUMF)

!     ------------------------------------------------------------------

!*         3.     SETTING CONSTANTS FOR NEW LINEARIZED CONVECTION SCHEME
!                 ------------------------------------------------------

CALL SUCUMF2(YDSTA,YDDIMV,YDMODEL%YRML_GCONF%YRRIP,YDMODEL%YRML_PHY_SLIN%YRCUMFS,&
           & YDMODEL%YRML_PHY_SLIN%YRECUMF2,IGRID)

!     ------------------------------------------------------------------
!*         4a.     SETTING CONSTANTS FOR GRAVITY WAVE DRAG SCHEME
!                 ----------------------------------------------

CALL SUGWD(YDSTA, YDMODEL%YRML_PHY_EC%YREGWD, YDEPHY, KULOUT, NFLEVG, YDVAB%VAH, YDVAB%VBH)

!*         4b.   SETTING CONSTANTS FOR WARNER-MCINTYRE-SCINOCCA NON-OROGRAPHIC GW SCHEME
!                -----------------------------------------------------------------------

CALL SUGWWMS(YDSTA,YDDIMV,YDMODEL%YRML_PHY_EC%YREGWWMS,YDMODEL%YRML_GCONF%YRRIP,IGRID)

!     ------------------------------------------------------------------

!*         5.     SETTING CONSTANTS FOR VERTICAL DIFFUSION
!                 ----------------------------------------

CALL SUVDFS

CALL SUVDF(YDMODEL%YRML_PHY_G%YRVDF)

!cccc CALL SUVDFD ( NABLPFR, ABLPLL ) cccccccccccccccccccccccccccccccccc

!     ------------------------------------------------------------------

!*         6.     SETTING CONSTANTS FOR RADIATION SCHEME
!                 --------------------------------------

CALL GSTATS(1811,2)
IF (LRAYFM15) THEN
  CALL ABOR1('RADIATION SCHEME FROM CYCLE 15 NO LONGER AVAILABLE')
ELSE
  CALL SUECRAD(YDGEOMETRY, YDMODEL, KULOUT, NFLEVG, ZETAH)
ENDIF
CALL GSTATS(1811,3)

!     ------------------------------------------------------------------
!*         7.     SETTING CONSTANTS FOR SURFACE SCHEME
!                 ------------------------------------

! Check minimum ice fraction and set appropiate values if not set.
IF (RCIMIN<0.0_JPRB) THEN
  IF (LMCCDYNSEAICE.AND.(.NOT.LNEMO1WAY)) THEN
    RCIMIN=0.005_JPRB
  ELSE 
    RCIMIN=0.02_JPRB
  ENDIF
ENDIF
WRITE(KULOUT,'(A,F7.4)')&
   & 'SUPHEC: Minimum ice fraction is ', RCIMIN

IF (LRAYFM15) THEN
  CALL ABOR1('RADIATION SCHEME FROM CYCLE 15 NO LONGER AVAILABLE')
ELSE
  ALLOCATE(ZSUN(SIZE(YDMODEL%YRML_PHY_RAD%YRESWRT%RSUN)))
  ZSUN(:)=YDMODEL%YRML_PHY_RAD%YRESWRT%RSUN(:)
ENDIF


CALL SUOPTSURF(YDEPHY,KULOUT)

CALL SUSURF(KSW=NSW,KCSS=NCSS,KCSNEC=NCSNEC,KSIL=NSIL,KCOM=NCOM,KTILES=NTILES,&
 & KTSW=NTSW,KLWEMISS=NLWEMISS,&
 & LD_LEFLAKE=LEFLAKE,LD_LEOCML=LEOCML,LD_LOCMLTKE=LOCMLTKE,&
 & LD_LWCOU=LWCOU,LD_LWCOU2W=LWCOU2W,LD_LWCOUHMF=LWCOUHMF, &
 & LD_LLCCNL=LCCNL,LD_LLCCNO=LCCNO,LD_LEVGEN=LEVGEN,LD_LESSRO=LESSRO,&
 & LD_LELAIV=LELAIV,LD_LECTESSEL=LECTESSEL,LD_LEAGS=LEAGS,&
 & LD_LEAIRCO2COUP=LEAIRCO2COUP, LD_LEFARQUHAR=LEFARQUHAR, &
 & LD_LESN09=LESN09,LD_LESNML=LESNML,&
 & LD_LEOCWA=LEOCWA,LD_LEOCCO=LEOCCO,LD_LEOCSA=LEOCSA,LD_LEOCLA=LEOCLA,&
 & KALBEDOSCHEME=NALBEDOSCHEME,KEMISSSCHEME=NEMISSSCHEME,&
 & LD_LSCMEC=LSCMEC,LD_LROUGH=LROUGH,PEXTZ0M=REXTZ0M,PEXTZ0H=REXTZ0H,&
 & PTHRFRTI=RTHRFRTI,PTSTAND=TSTAND,PXP=XP,PRCCNSEA=RCCNSEA,PRCCNLND=RCCNLND,&
 & PRLAIINT=RLAIINT,PRSUN=ZSUN,PRCORIOI=RCORIOI,PRPLRG=RPLRG,PNSNMLWS=NSNMLWS,&
 & PRVR0VT=RVR0VT,PRVCMAX25=RVCMAX25,PRHUMREL=RHUMREL,PRA1=RA1,PRB1=RB1,PRG0=RG0,PRGM25=RGM25,&
 & PRE_VCMAX=RE_VCMAX,PRE_JMAX=RE_JMAX,&
 & YDSURF=YSURF,PRALFMINPSN=RALFMINPSN,PRCIMIN=RCIMIN)
DEALLOCATE(ZSUN)




!          7.1    Allocate working arrays
! See SUTILEPROP
CALL GSTATS(1811,1)

!     ------------------------------------------------------------------

!*         8.     SETTING CONSTANTS FOR CLOUD OPTICAL PROPERTIES
!                 ----------------------------------------------

IF (LRAYFM15) THEN
  CALL ABOR1('RADIATION SCHEME FROM CYCLE 15 NO LONGER AVAILABLE')
ELSEIF (.NOT.LRAY) THEN
  CALL SUCLOPN(YDERAD,YDESWRT,NTSW,NSW,NFLEVG)
ENDIF

! Set up constants for 550nm optical depth diagnostic  
CALL SU_CLOP550

!     ------------------------------------------------------------------

!*         9.     SETTING CONSTANTS FOR PROGNOSTIC CLOUD SCHEME
!                 ----------------------------------------------

CALL SUCLDP(YDSTA,YDDIMV,YDPHY2,YDMODEL%YRML_PHY_EC%YRECLDP)

!     ------------------------------------------------------------------

!*        10.     SETTING CONSTANTS FOR WAVE COUPLING
!                 -----------------------------------

CALL SUWCOU(YDMODEL%YRML_GCONF%YRRIP,YDMODEL%YREWCOU)

!     ------------------------------------------------------------------
!*         11.   SETTING CONSTANTS FOR LINEARIZED PHYSICS
!                ----------------------------------------

CALL SUPHLI(YDSTA,YDDIMV,YDMODEL)

!     ------------------------------------------------------------------
!*         12.   SETTING CONSTANTS FOR METHANE OXIDATION
!                ---------------------------------------

CALL SUMETHOX

!     ------------------------------------------------------------------
!*         13.   SETTING CONSTANTS FOR PROGNOSTIC AEROSOLS
!                -----------------------------------------

! This will have already been called "early" if either prognostic aerosol
! or the CAMS climatology is being used, in which case do NOT call it again
! as this is unnecessary and will wipe the optical properties that have been
! loaded from file.
IF (LEPHYS .AND. NAERMACC /= 1 .AND. NACTAERO == 0) THEN
  ! Expensive subroutine fur Fullpos 927:
  CALL SU_AERW(YDMODEL)
ENDIF


!     ------------------------------------------------------------------
!*         14.   SETTING CONSTANTS (and stuff) FOR WAVE MODEL
!                -----------------------------------------

IF (LWCOU .AND. NFPOS/=2) CALL SUWAM(YDGEOMETRY,YDMODEL%YREWCOU)

!     ------------------------------------------------------------------

WRITE(UNIT=KULOUT,FMT='('' SUPHEC IS OVER '')')

!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUPHEC',1,ZHOOK_HANDLE)
END SUBROUTINE SUPHEC
