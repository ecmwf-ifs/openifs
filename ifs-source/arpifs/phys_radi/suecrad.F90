! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SUECRAD(YDGEOMETRY,YDMODEL,KULOUT,KLEV,PETAH)

!**** *SUECRAD*   - INITIALIZE COMMONS YOERxx CONTROLLING RADIATION

!     PURPOSE.
!     --------
!           INITIALIZE YOERAD, THE COMMON THAT CONTROLS THE
!           RADIATION OF THE MODEL, AND YOERDU THAT INCLUDES
!           ADJUSTABLE PARAMETERS FOR RADIATION COMPUTATIONS

!**   INTERFACE.
!     ----------
!        CALL *SUECRAD* FROM *SUPHEC*
!              -------        ------

!        EXPLICIT ARGUMENTS :
!        --------------------
!        NONE

!        IMPLICIT ARGUMENTS :
!        --------------------
!        COMMONS YOERAD, YOERDU

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------
!        SUAER, SUAERH, SUAERV, SULW, SUSW, SUOCST, SUSAT
!        SUAERL, SUAERSN, SUSRTAER, SRTM_INIT, SUSRTCOP

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     AUTHOR.
!     -------
!        JEAN-JACQUES MORCRETTE  *ECMWF*
!        ORIGINAL : 88-12-15

!     MODIFICATIONS.
!     --------------
!        R. El Khatib 01-02-02 proper initialization of NFRRC moved in SUCFU
!        010129 JJMorcrette clean-up LERAD1H, NLNGR1H
!        011105 GMozdzynski support new radiation grid
!        011005 JJMorcrette CCN --> Re Water clouds
!        R. El Khatib 01-02-02 LRRTM=lecmwf by default
!        020909 GMozdzynski support NRADRES to specify radiation grid
!        021001 GMozdzynski support on-demand radiation communications
!        030422 GMozdzynski automatic min-halo
!        030501 JJMorcrette new radiation grid on, new aerosols on (default)
!        030513 JJMorcrette progn. O3 / radiation interactions off (default)
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        050315 JJMorcrette prog.aerosols v1
!        041214 JJMorcrette SRTM
!        050111 JJMorcrette new cloud optical properties
!        050415 GMozdzynski Reduced halo support for radiation interpolation
!        051004 JJMorcrette UV surface radiation processor
!        051220 JJMorcrette SRTM112g+LWSCAT+UVprocessor+(bgfx:swclr, radaca)
!        060113 M.Janiskova constants for lw transmission functions
!        060510 JJMorcrette MODIS albedo (UVis, NIR)x(parallel+diffuse)
!        JJMorcrette 20060721 PP of clear-sky PAR and TOA incident solar radiation
!        060625 JJMorcrette MODIS albedo (UVis, NIR)x(parallel+diffuse)
!        060726 JJMorcrette McICA default operational configuration
!        JJMorcrette 20080318 O.Boucher sulphate obs: 1920-1990, A1B: 2000-2100
!        080423 JJMorcrette 3D climatologies for CO2, CH4, N2O
!        080520 M.Janiskova NPROMALW for AD longwave radiation
!        081002 Y.Bouteloup : Two read of namelist NAERAD (usefull when LSRTM=FALSE)
!        JJMorcrette 20090213 GEMS-derived CO2, CH4, O3 climatologies
!        090113 JJMorcrette/MJIacono Look-up table for exponentials in SW
!        JJMorcrette 20090804 Change ice cloud De as f(latitude)
!        JJMorcrette 20090805 Decorrelation for CF and CW as f(latitude)
!        091201 JJMorcrette Total and clear-sky direct SW radiation flux at surface 
!        010411 HHersbach Initialization of various new CMIP5 switches
!        G. Mozdzynski (Aug 2011): support higher order interpolatio
!        31 Oct 2011  M Ahlgrimm Surface downward clear-sky LW and SW fluxes
!        JJMorcrette 20111007 Possiblity of different frequency of full rad. for EPS
!        R. El Khatib 22-Mar-2012 Fix uninitialized variables
!        P Bechtold 14/05/2012 replace 86400 by RDAY
!        G. Mozdzynski (May 2012): further cleaning
!        N.Semane+P.Bechtold   04-10-2012 replace 3600s by RHOUR and add RPLRADI for small planet
!        T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!        2013-11, J. Masek: Ensure NTSW=NSW=NSWB_MNH=1 for LRAY=.T.
!         (compatibility of ACRANEB/ACRANEB2 with SURFEX).
!        K. Yessad (July 2014): Move some variables.
!        JJMorcrette/ABozzo 20130805/201410 MACC-derived aerosol climatology
!        R J Hogan (June 2014) Added approx LW and SW updates
!        R J Hogan (Nov  2014) Added centred-time solar zenith angle
!        A. Bozzo (Oct 2014) new switch for the revised version of the UV processor
!        R J Hogan (15 Apr 2015) Added LMannersSwUpdate
!        R J Hogan (24 Apr 2015) Added LAverageSZA
!        R. El Khatib 27-Apr-2015 dummy init. if NWS==1
!        2015-5, E.Gleeson, L. Rontu and K.P. Nielsen: Default settings + cleaning
!        G Mozdzynski (Aug 2015) Correction to radiation halo width calculation
!        R J Hogan (Sept 2015) Added setup of new radiation scheme
!        R J Hogan (17 Dec 2015) Default NHINCSOL=3: lower solar irradiance with a solar cycle
!        R J Hogan (Mar 2016)  Added NLwScattering, NSwSolver, NLwSolver, LUsePre2017Rad
!        2016-4, L.Rontu: Further streamlining and cleaning
!        R J Hogan (Jun 2016)  Added RCLOUD_FRAC_STD, LFU_LW_ICE_OPTICS_BUG
!        R J Hogan (Sep 2016)  Default now NLIQOPT=4 for newer radiation scheme
!        C Roberts/R Senan (Feb 2017) Support for CMIP6 forcings
!        R J Hogan (Mar 2016)  Pass YRADIATION to SETUP_RADIATION_SCHEME
!        A Bozzo (Jan 2017) new switches for CAMS 3D aerosol climatology
!        R J Hogan (Feb 2018)  Added NSWWVCONTINUUM for CAVIAR water vapour continuum
!        F Vana (Feb 2018) More consistent setup when LERADI=false.
!        R J Hogan (Mar 2018)  Added NDUMPBADINPUTS
!        R J Hogan (Jan 2019)  Added LINTERINCLOUDMEAN
!        R J Hogan (Mar 2019)  Setup greenhouse gas climatology and timeseries here
!        R J Hogan (Mar 2019)  Set default GHG scenario to CMIP6 SSP3-7.0
!        R J Hogan (Jan 2019)  Removed YDERAD%LEMODAL: it does nothing and a copy of YDEPHY%LE4ALB
!        R J Hogan (Jan 2019)  Added NCLOUDOVERLAP, NLWEMISS, NDUMPINPUTS, NLWOUT
!        R J Hogan (Feb 2019)  Added RCLOUD_SEPARATION_SCALE_[TOA|SURF]
!        R J Hogan (Mar 2019)  Default radiation frequency 1h at all resolutions
!        R Forbes  (Nov 2020)  Changed NCLOUDOVERLAP default to 3 (exp-ran)
!        R Forbes  (Nov 2020)  Changed default LINTERINCLOUDMEAN to F as radintg bug now fixed
!        R J Hogan (Jun 2021)  Turned on interactive ozone by default (LEPO3RA) if prognostic
!        R J Hogan (Jun 2021)    ozone is available, and made NGHGRAD consistent with LEPO3RA
!        R J Hogan (Oct 2021)  Updated ozone climatology (if used) in greenhouse_gas_climatology_48r1.nc
! End Modifications
!-------------------------------------------------------------------------------

USE TYPE_MODEL        , ONLY : MODEL
USE GEOMETRY_MOD      , ONLY : GEOMETRY
USE PARKIND1          , ONLY : JPIM, JPRB, JPRD
USE YOMHOOK           , ONLY : LHOOK, DR_HOOK, JPHOOK
USE PARDIM            , ONLY : JPMXGL
USE YOESRTM           , ONLY : JPGPT
USE YOMCT0            , ONLY : LECMWF, LALLOPR, LARPEGEF, L_OOPS
USE YOMLUN            , ONLY : NULNAM, NULRAD, NULOUT
USE YOMRIP0           , ONLY : NINDAT
USE YOMCST            , ONLY : RPI, RI0, RNAVO , RDAY, RHOUR
USE YOERDU            , ONLY : NUAER, NTRAER, R10E,&
 &                             REPLOG, REPSC, REPSCA, REPSCO, REPSCQ, REPSCT, REPSCW, DIFF  
USE YOMMP0            , ONLY : NPROC, MYPROC, NPRCIDS, LSPLIT, MY_REGION_NS, MY_REGION_EW,&
 &                             N_REGIONS_NS, N_REGIONS_EW, LOUTPUT, LOPT_SCALAR, NPRINTLEV
USE YOEAERCLI         , ONLY : NAERADJDU, RDUMULF, RWGHTDU1, RWGHTDU2, RWGHTDU3
USE OML_MOD           , ONLY : OML_MAX_THREADS
USE YOMTAG            , ONLY : MTAGRAD
USE MPL_MODULE        , ONLY : MPL_BROADCAST, MPL_SEND, MPL_RECV
USE YOMDYNCORE        , ONLY : LAPE, RPLRADI
USE YOELWCONST        , ONLY : RCH4A, RCH4B, RCN2OA, RCN2OB
USE YOEAERC           , ONLY : LVOLCDATA, CLISTSO4, CVOLCDATA
USE YOMTRANS          , ONLY : LFFTW
USE YOMCLIM           , ONLY : YGHGCLIM ! Climatology from NetCDF file
USE YOMGHGTIMESERIES  , ONLY : YGHGTIMESERIES ! GHG multi-annual timeseries
USE YOMSOLARIRRADIANCE, ONLY :YSOLARIRRADIANCE ! Total solar irradiance multi-annual timeseries
USE YOERADGHG         , ONLY : YRADGHG  ! Climatology interpolated to current model time
USE RADIATION_SETUP   , ONLY : SETUP_RADIATION_SCHEME

!-------------------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(MODEL),TARGET,INTENT(INOUT) :: YDMODEL
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KULOUT 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PETAH(KLEV+1) 

!-------------------------------------------------------------------------------

!     LOCAL ARRAYS FOR THE PURPOSE OF READING NAMRGRI (RADIATION GRID)
INTEGER(KIND=JPIM) :: NRGRI(JPMXGL)

!     LOCAL ARRAYS FOR THE PURPOSE OF READING NAERAD
CHARACTER (LEN=256) :: CRTABLEDIR
CHARACTER (LEN=32)  :: CRTABLEFIL

INTEGER(KIND=JPIM) :: IDGL,IDGLG2,INBLW,IRADFR,IST1HR,ISTNHR,IDIR,IFIL
INTEGER(KIND=JPIM) :: JLON,JGLAT,JGL,JGLSUR,IDLSUR,IOFF,ILAT,ISTLON,IENDLON
INTEGER(KIND=JPIM) :: ILBRLATI,IUBRLATI,IGLGLO,IDUM,IU,IGG
INTEGER(KIND=JPIM) :: I,J,JJ,JROC,IGPTOT
INTEGER(KIND=JPIM) :: IROWIDEMAXN,IROWIDEMAXS,IROWIDEMAXW,IROWIDEMAXE
INTEGER(KIND=JPIM) :: IRIWIDEMAXN,IRIWIDEMAXS,IRIWIDEMAXW,IRIWIDEMAXE
INTEGER(KIND=JPIM) :: IARIB1MAX,IAROB1MAX
INTEGER(KIND=JPIM) :: IWIDE(10)
INTEGER(KIND=JPIM) :: ILATS_DIFF_F,ILATS_DIFF_C, ILONS_DIFF_F,ILONS_DIFF_C
INTEGER(KIND=JPIM) :: ITHRMAX,IRPROMA,IRPROMABEG,IRPROMAEND,IGPBLKSMX,IMX,IPR
INTEGER(KIND=JPIM) :: I_MIN_HALO
INTEGER(KIND=JPIM) :: ISW,JUV,IDAYUV
INTEGER(KIND=JPIM) :: IDGSAH,IDGENH

LOGICAL :: LLINEAR_GRID, LLREDURAD
LOGICAL :: LLDEBUG,LLP,LLGRIDONLY

REAL(KIND=JPRD) :: ZSTPHR, ZTSTEP, ZGEMU, ZLON, ZD1, ZD2, ZD3, ZD4, ZD5, ZD6
REAL(KIND=JPRD) :: ZMINRADLAT,ZMAXRADLAT,ZMINRADLON,ZMAXRADLON
REAL(KIND=JPRD) :: ZMINMDLLAT,ZMAXMDLLAT,ZMINMDLLON,ZMAXMDLLON
REAL(KIND=JPRD) :: ZLAT

REAL(KIND=JPRD) :: ZCCO2, ZLAMBDA, Z1LMBD2, Z1LMBD4, ZLAMBDACM, ZLMBDCM4
REAL(KIND=JPRD) :: ZFN2, ZFO2, ZNFAIR, ZDFAIR, ZFAIR, ZRIM1_300, ZRIAIRM1
REAL(KIND=JPRD) :: ZRIAIR, ZRIAIR2

CHARACTER (LEN = 512) ::  CLF1
INTEGER(KIND=JPIM), PARAMETER :: JPIOMASTER=1

INTEGER(KIND=JPIM), ALLOCATABLE :: IGLOBALINDEX(:)

REAL(KIND=JPRD),ALLOCATABLE :: ZLATX(:)
REAL(KIND=JPRD),ALLOCATABLE :: ZLONX(:)
CHARACTER(LEN = 512) :: CLZZZ

REAL(KIND=JPRD), ALLOCATABLE :: ZMU(:)

LOGICAL         :: LACR1, LACR6
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
#ifdef __INTEL_COMPILER
INTRINSIC :: SHAPE ! Fails with Intel compiler as it thinks we use ATLAS shape function
#endif


!-------------------------------------------------------------------------------

#include "setup_trans.h"
#include "trans_inq.h"

#include "abor1.intfb.h"
#include "posnam.intfb.h"
#include "rrtm_init_140gp.intfb.h"

#include "slcset.intfb.h"
#include "suaerh.intfb.h"
#include "suaerl.intfb.h"
#include "suaersn.intfb.h"
#include "suaerv.intfb.h"
#include "su_aerv.intfb.h"
#include "suclopn.intfb.h"
#include "sulwn.intfb.h"
#include "sulwneur.intfb.h"
#include "suovlp.intfb.h"
#include "surdi.intfb.h"
#include "surrtab.intfb.h"
!  #include "surrtftr.intfb.h"
#include "surrtpk.intfb.h"
#include "surrtrf.intfb.h"

#include "suswn.intfb.h"
#include "susrtab.intfb.h"
#include "susrtaer.intfb.h"
#include "srtm_init.intfb.h"
#include "susrtcop.intfb.h"
#include "su_aerw.intfb.h"
#include "su_uvradi.intfb.h"
#include "su_mcica.intfb.h"
#include "su_clop550.intfb.h"
#include "suecozc.intfb.h"
#include "suecozcaqua.intfb.h"
#include "suecozv.intfb.h"

!-------------------------------------------------------------------------------

LOGICAL, POINTER :: LERAD1H, LEPO3RA, LONEWSW, LCCNL, LCCNO, LRAYL, LRRTM,&
 & LSRTM, LHVOLCA, LNEWAER, LDIFFC, LNOTROAER, LETRACGMS, LAERCLIM,&
 & LAERVISI, LAERADJDU, LAERADCLI, LPERPET, LVOLCSPEC, LVOLCDAMP, LECO2VAR,&
 & LHGHG, LUVPROC, LUVTDEP, LUVDBG, LEDBUG, LESO4HIS,&
 & LDIAGFORCING, LUVAERP, LAPPROXLWUPDATE, LAPPROXSWUPDATE,&
 & LCENTREDTIMESZA, LMANNERSSWUPDATE, LAVERAGESZA, LECSRAD, LECOMPGRID,&
 & LDUSEASON, LUSEPRE2017RAD, LFU_LW_ICE_OPTICS_BUG, LO3_CHEM_UV, &
 & LAER3D, LINTERPINCLOUDMEAN

INTEGER(KIND=JPIM), POINTER :: NICEOPT, NLIQOPT, NSWICEOPT, NLWICEOPT, NSWLIQOPT,&
 & NLWLIQOPT,  NMCICA, NRADIP, NRADLP,&
 & NREDGLW, NREDGSW, NAER, NMODE, NOZOCL, NINHOM, NLAYINH, NOVLP, NSW,&
 & NRADFR, NLNGR1H, NRADELG, NRADPFR, NRADPLA, NRPROMA, NRADINT, NRADRES,&
 & NGHGRAD, NDECOLAT, NMINICE, NVOLCVERT, NPERTAER, NPERTOZ, NHINCSOL, NSCEN,&
 & NAERMACC, NMCLAT, NMCLON, NUV, NUVTIM, NRADUV, KMODTS, NSWSOLVER, NLWSOLVER, &
 & NO3CMIP, NGHGCMIP, NRCP, NCMIPFIXYR, &
 & NLWSCATTERING, NSOLARSPECTRUM, NSWWVCONTINUUM, NDUMPBADINPUTS, &
 & NCLOUDOVERLAP, NDUMPINPUTS

REAL(KIND=JPRB), POINTER :: RVOLCSPEC(:), RCCNSEA, RCCNLND, RPERTOZ, RLWINHF,&
 & RSWINHF, RRE2DE, RMINICE, RUVLAM(:), RMUZUV, RCCO2, RCCH4, RCN2O, RCCFC11, RCCFC12,&
 & RAESHSS, RAESHDU, RAESHOM, RAESHBC, RAESHSU, TRBKG, STBKG, RCLOUD_FRAC_STD, &
 & RCLOUD_SEPARATION_SCALE_TOA, RCLOUD_SEPARATION_SCALE_SURF

CHARACTER(LEN=260), POINTER :: CO3DATADIR, GHGDATADIR, CSOLARDATA
CHARACTER (LEN=80), POINTER :: CO3DATAFIL
CHARACTER (LEN=256), POINTER :: CGHGCLIMFILE, CGHGTIMESERIESFILE, CSOLARIRRADIANCEFILE

#include "naerad.nam.h"
#include "naercli.nam.h"
#include "namrgri.nam.h"

!-------------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUECRAD',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, &
 & YDGSGEOM_NB=>YDGEOMETRY%YRGSGEOM_NB, YDCSGLEG=>YDGEOMETRY%YRCSGLEG,  YDSTA=>YDGEOMETRY%YRSTA, &
 & YDVAB=>YDGEOMETRY%YRVAB, YDVETA=>YDGEOMETRY%YRVETA, YDVFE=>YDGEOMETRY%YRVFE, &
 & YDLAP=>YDGEOMETRY%YRLAP, &
 & YDCSGEOM_NB=>YDGEOMETRY%YRCSGEOM_NB,  YDSPGEOM=>YDGEOMETRY%YSPGEOM,  &
 & YDPHY=>YDMODEL%YRML_PHY_MF%YRPHY,YDEUVRAD=>YDMODEL%YRML_PHY_RAD%YREUVRAD,YDPARAR=>YDMODEL%YRML_PHY_MF%YRPARAR, &
 & YDERDI=>YDMODEL%YRML_PHY_RAD%YRERDI, YDESWRT=>YDMODEL%YRML_PHY_RAD%YRESWRT, &
 & YDRIP=>YDMODEL%YRML_GCONF%YRRIP,YDERAD=>YDMODEL%YRML_PHY_RAD%YRERAD,YDELWRAD=>YDMODEL%YRML_PHY_RAD%YRELWRAD,  &
 & YGFL=>YDMODEL%YRML_GCONF%YGFL, &
 & YDEPHY=>YDMODEL%YRML_PHY_EC%YREPHY,YDEAERATM=>YDMODEL%YRML_PHY_RAD%YREAERATM, &
 & YDCOMPO=>YDMODEL%YRML_CHEM%YRCOMPO, &
 & YDEAERD=>YDMODEL%YRML_PHY_RAD%YREAERD,  RADGRID=>YDMODEL%YRML_PHY_RAD%RADGRID, &
 & YDPHNC=>YDMODEL%YRML_PHY_SLIN%YRPHNC, YDRADF=>YDMODEL%YRML_PHY_RAD%YRRADF, YDRI=>YDMODEL%YRML_PHY_RAD%YRRI, &
 & YDRO=>YDMODEL%YRML_PHY_RAD%YRRO,YDRADIATION=>YDMODEL%YRML_PHY_RAD%YRADIATION)

ASSOCIATE(YO3=>YGFL%YO3, NACTAERO=>YGFL%NACTAERO, &
 & NDGENG=>YDDIM%NDGENG, NDGENH=>YDDIM%NDGENH, NDGENL=>YDDIM%NDGENL, &
 & NDGLG=>YDDIM%NDGLG, NDGSAG=>YDDIM%NDGSAG, NDGSAH=>YDDIM%NDGSAH, &
 & NDGSAL=>YDDIM%NDGSAL, NDGSUR=>YDDIM%NDGSUR, NDLON=>YDDIM%NDLON, &
 & NDLSUR=>YDDIM%NDLSUR, NDSUR1=>YDDIM%NDSUR1, NGPBLKS=>YDDIM%NGPBLKS, &
 & NPROMA=>YDDIM%NPROMA, NSMAX=>YDDIM%NSMAX, NSTENCILWIDE=>YDDIM%NSTENCILWIDE, &
 & NFLEVG=>YDDIMV%NFLEVG, &
 & LAERUVP=>YDEAERATM%LAERUVP, &
 & CVDAED=>YDEAERD%CVDAED, CVDAEL=>YDEAERD%CVDAEL, CVDAES=>YDEAERD%CVDAES, &
 & CVDAEU=>YDEAERD%CVDAEU, RCAEADK=>YDEAERD%RCAEADK, RCAEADM=>YDEAERD%RCAEADM, &
 & RCAEOPD=>YDEAERD%RCAEOPD, RCAEOPL=>YDEAERD%RCAEOPL, RCAEOPS=>YDEAERD%RCAEOPS, &
 & RCAEOPU=>YDEAERD%RCAEOPU, RCAEROS=>YDEAERD%RCAEROS, RCSTBGA=>YDEAERD%RCSTBGA, &
 & RCTRBGA=>YDEAERD%RCTRBGA, RCTRPT=>YDEAERD%RCTRPT, RCVOBGA=>YDEAERD%RCVOBGA, &
 & LOPTLWPR=>YDELWRAD%LOPTLWPR, NLASTLW=>YDELWRAD%NLASTLW, &
 & NLOOPLW=>YDELWRAD%NLOOPLW, NLOOPLWO=>YDELWRAD%NLOOPLWO, &
 & NPROMALW=>YDELWRAD%NPROMALW, &
 & LEPHYS=>YDEPHY%LEPHYS, LERADI=>YDEPHY%LERADI, &
 & CVDAEBC=>YDERAD%CVDAEBC, CVDAEDU=>YDERAD%CVDAEDU, CVDAEOM=>YDERAD%CVDAEOM, &
 & CVDAESS=>YDERAD%CVDAESS, CVDAESU=>YDERAD%CVDAESU, &
 & LOPTRPROMA=>YDERAD%LOPTRPROMA, NMCLEV=>YDERAD%NMCLEV, NMCVAR=>YDERAD%NMCVAR, &
 & NRADE1H=>YDERAD%NRADE1H, NRADE3H=>YDERAD%NRADE3H, NRADNFR=>YDERAD%NRADNFR, &
 & NRADSFR=>YDERAD%NRADSFR, NSPMAPL=>YDERAD%NSPMAPL, NSPMAPS=>YDERAD%NSPMAPS, &
 & NSWNL=>YDERAD%NSWNL, NSWTL=>YDERAD%NSWTL, NTSW=>YDERAD%NTSW, &
 & RNS=>YDERAD%RNS, RSIGAIR=>YDERAD%RSIGAIR, &
 & RCCCL4=>YDERDI%RCCCL4, RCCFC22=>YDERDI%RCCFC22, RCNO2=>YDERDI%RCNO2, &
 & RSOLINC=>YDERDI%RSOLINC, &
 & JUVLAM=>YDEUVRAD%JUVLAM, &
 & NGPTOT=>YDGEM%NGPTOT, NGPTOTG=>YDGEM%NGPTOTG, NGPTOTMX=>YDGEM%NGPTOTMX, &
 & NLOENG=>YDGEM%NLOENG, &
 & MYFRSTACTLAT=>YDMP%MYFRSTACTLAT, MYLSTACTLAT=>YDMP%MYLSTACTLAT, &
 & NFRSTLAT=>YDMP%NFRSTLAT, NFRSTLOFF=>YDMP%NFRSTLOFF, &
 & NGLOBALINDEX=>YDMP%NGLOBALINDEX, NLSTLAT=>YDMP%NLSTLAT, NONL=>YDMP%NONL, &
 & NPTRFLOFF=>YDMP%NPTRFLOFF, NPTRFRSTLAT=>YDMP%NPTRFRSTLAT, NSTA=>YDMP%NSTA, &
 & LERADLW2=>YDPHNC%LERADLW2, LERADN2=>YDPHNC%LERADN2, &
 & TSTEP=>YDRIP%TSTEP, NSWB_MNH=>YDPARAR%NSWB_MNH, LMPHYS=>YDPHY%LMPHYS, &
 & LRAY=>YDPHY%LRAY, LRAYFM=>YDPHY%LRAYFM,LRAYFM15=>YDPHY%LRAYFM15,&
 & LPHYLIN=>YDMODEL%YRML_PHY_SLIN%YREPHLI%LPHYLIN, &
 & LHLRADUPD=>YDPHY%LHLRADUPD, RADGR=>YDPARAR%RADGR, RADSN=>YDPARAR%RADSN)
! Associate pointers for variables in namelist NAERAD
LERAD1H         => YDERAD%LERAD1H
LEPO3RA         => YDERAD%LEPO3RA
LONEWSW         => YDERAD%LONEWSW
LCCNL           => YDERAD%LCCNL
LCCNO           => YDERAD%LCCNO
LRAYL           => YDERAD%LRAYL
LRRTM           => YDERAD%LRRTM
LSRTM           => YDERAD%LSRTM
LHVOLCA         => YDERAD%LHVOLCA
LNEWAER         => YDERAD%LNEWAER
LDIFFC          => YDERAD%LDIFFC
LNOTROAER       => YDERAD%LNOTROAER
LETRACGMS       => YDERAD%LETRACGMS
LAERCLIM        => YDERAD%LAERCLIM
LAERVISI        => YDERAD%LAERVISI
LAERADJDU       => YDERAD%LAERADJDU
LAERADCLI       => YDERAD%LAERADCLI
LDUSEASON       => YDERAD%LDUSEASON
LAER3D          => YDERAD%LAER3D
LPERPET         => YDERAD%LPERPET
NICEOPT         => YDERAD%NICEOPT
NLIQOPT         => YDERAD%NLIQOPT
NSWICEOPT       => YDERAD%NSWICEOPT
NLWICEOPT       => YDERAD%NLWICEOPT
NSWLIQOPT       => YDERAD%NSWLIQOPT
NLWLIQOPT       => YDERAD%NLWLIQOPT
NMCICA          => YDERAD%NMCICA
NRADIP          => YDERAD%NRADIP
NRADLP          => YDERAD%NRADLP
NREDGLW         => YDERAD%NREDGLW
NREDGSW         => YDERAD%NREDGSW
NAERMACC        => YDERAD%NAERMACC
NMCLAT          => YDERAD%NMCLAT
NMCLON          => YDERAD%NMCLON
NAER            => YDERAD%NAER
NMODE           => YDERAD%NMODE
NOZOCL          => YDERAD%NOZOCL
NINHOM          => YDERAD%NINHOM
NLAYINH         => YDERAD%NLAYINH
NOVLP           => YDERAD%NOVLP
NSW             => YDERAD%NSW
NRADFR          => YDERAD%NRADFR
NLNGR1H         => YDERAD%NLNGR1H
NRADELG         => YDERAD%NRADELG
NRADPFR         => YDERAD%NRADPFR
NRADPLA         => YDERAD%NRADPLA
NRPROMA         => YDERAD%NRPROMA
NRADINT         => YDERAD%NRADINT
NRADRES         => YDERAD%NRADRES
NGHGRAD         => YDERAD%NGHGRAD
NDECOLAT        => YDERAD%NDECOLAT
NMINICE         => YDERAD%NMINICE
NVOLCVERT       => YDERAD%NVOLCVERT
LVOLCSPEC       => YDERAD%LVOLCSPEC
LVOLCDAMP       => YDERAD%LVOLCDAMP
RVOLCSPEC       => YDERAD%RVOLCSPEC
RCCNSEA         => YDERAD%RCCNSEA
RCCNLND         => YDERAD%RCCNLND
NPERTAER        => YDERAD%NPERTAER
NPERTOZ         => YDERAD%NPERTOZ
RPERTOZ         => YDERAD%RPERTOZ
RLWINHF         => YDERAD%RLWINHF
RSWINHF         => YDERAD%RSWINHF
RRE2DE          => YDERAD%RRE2DE
RMINICE         => YDERAD%RMINICE
RCCO2           => YDERDI%RCCO2
RCCH4           => YDERDI%RCCH4
RCN2O           => YDERDI%RCN2O
RCCFC11         => YDERDI%RCCFC11
RCCFC12         => YDERDI%RCCFC12
RAESHSS         => YDERAD%RAESHSS
RAESHDU         => YDERAD%RAESHDU
RAESHOM         => YDERAD%RAESHOM
RAESHBC         => YDERAD%RAESHBC
RAESHSU         => YDERAD% RAESHSU 
NHINCSOL        => YDERAD%NHINCSOL
LECO2VAR        => YDERAD%LECO2VAR
LHGHG           => YDERAD%LHGHG
NSCEN           => YDERAD%NSCEN
LUVPROC         => YDEUVRAD%LUVPROC
LUVTDEP         => YDEUVRAD%LUVTDEP
LUVDBG          => YDEUVRAD%LUVDBG
NUV             => YDERAD%NUV
LEDBUG          => YDERAD%LEDBUG
LESO4HIS        => YDERAD%LESO4HIS
NUVTIM          => YDEUVRAD%NUVTIM
NRADUV          => YDEUVRAD%NRADUV
RUVLAM          => YDEUVRAD%RUVLAM
RMUZUV          => YDEUVRAD%RMUZUV
LDIAGFORCING    => YDERAD%LDIAGFORCING
KMODTS          => YDERAD%KMODTS
LUVAERP         => YDEUVRAD%LUVAERP
LO3_CHEM_UV     => YDEUVRAD%LO3_CHEM_UV
LECOMPGRID      => YDERAD%LECOMPGRID
LECSRAD         => YDERAD%LECSRAD
LAPPROXLWUPDATE => YDERAD%LAPPROXLWUPDATE
LAPPROXSWUPDATE => YDERAD%LAPPROXSWUPDATE
LCENTREDTIMESZA => YDERAD%LCENTREDTIMESZA
LMANNERSSWUPDATE=> YDERAD%LMANNERSSWUPDATE
LAVERAGESZA     => YDERAD%LAVERAGESZA
LUSEPRE2017RAD  => YDERAD%LUSEPRE2017RAD
NLWSCATTERING   => YDERAD%NLWSCATTERING
NSWSOLVER       => YDERAD%NSWSOLVER
NLWSOLVER       => YDERAD%NLWSOLVER
NCLOUDOVERLAP   => YDERAD%NCLOUDOVERLAP
RCLOUD_FRAC_STD => YDERAD%RCLOUD_FRAC_STD
RCLOUD_SEPARATION_SCALE_TOA  => YDERAD%RCLOUD_SEPARATION_SCALE_TOA
RCLOUD_SEPARATION_SCALE_SURF => YDERAD%RCLOUD_SEPARATION_SCALE_SURF
LFU_LW_ICE_OPTICS_BUG => YDERAD%LFU_LW_ICE_OPTICS_BUG
NSOLARSPECTRUM  => YDERAD%NSOLARSPECTRUM
NSWWVCONTINUUM  => YDERAD%NSWWVCONTINUUM
NDUMPBADINPUTS  => YDERAD%NDUMPBADINPUTS
NDUMPINPUTS     => YDERAD%NDUMPINPUTS
LINTERPINCLOUDMEAN => YDERAD%LINTERPINCLOUDMEAN
CGHGCLIMFILE    => YDERAD%CGHGCLIMFILE
CGHGTIMESERIESFILE => YDERAD%CGHGTIMESERIESFILE
CSOLARIRRADIANCEFILE => YDERAD%CSOLARIRRADIANCEFILE

TRBKG           => YDERAD%TRBKG
STBKG           => YDERAD%STBKG
NSWICEOPT       => YDERAD%NSWICEOPT
NSWLIQOPT       => YDERAD%NSWLIQOPT
NLWICEOPT       => YDERAD%NLWICEOPT
NLWLIQOPT       => YDERAD%NLWLIQOPT

NRCP            => YDRIP%YRECMIP%NRCP
NO3CMIP         => YDRIP%YRECMIP%NO3CMIP
NGHGCMIP        => YDRIP%YRECMIP%NGHGCMIP
NCMIPFIXYR      => YDRIP%YRECMIP%NCMIPFIXYR

GHGDATADIR      => YDRIP%YRECMIP%GHGDATADIR
CSOLARDATA      => YDRIP%YRECMIP%CSOLARDATA
CO3DATADIR      => YDRIP%YRECMIP%CO3DATADIR
CO3DATAFIL      => YDRIP%YRECMIP%CO3DATAFIL

!-------------------------------------------------------------------------------

!*         0.       PREINITIALIZE DIMENSIONINGS USED LATER IN SL
!                   --------------------------------------------

YDRI%NSLWIDEN=0
YDRI%NSLWIDES=0
YDRI%NSLWIDEW=0
YDRI%NSLWIDEE=0
YDRO%NSLWIDEN=0
YDRO%NSLWIDES=0
YDRO%NSLWIDEW=0
YDRO%NSLWIDEE=0
I_MIN_HALO=5

!*         1.       INITIALIZE NEUROFLUX LONGWAVE RADIATION
!                   ---------------------------------------

CALL GSTATS(1818,0)

IF (LERADN2) THEN
  CALL SULWNEUR(YDMODEL%YRML_PHY_RAD%YRENEUR,YDPHNC,KLEV)
ENDIF

!*         2.       SET DEFAULT VALUES.
!                   -------------------

IDUM=0
!*         2.1      PRESET INDICES IN *YOERAD*
!                   --------------------------

LERAD1H=.FALSE.
NLNGR1H=12

LONEWSW=.TRUE.

! If TRUE, compute the partial derivative of upwelling longwave flux
! with respect to surface upwelling longwave flux, and then use this
! to update the net longwave flux each timestep/gridpoint using the
! local value of skin temperature.  This corrects for errors in
! surface temperature at coastal points when the surface net longwave
! flux from the most recent call to the radiation scheme at a coarser
! resolution at a sea point was used at an adjacent land point. By
! adjusting also the net longwave flux profile and consequently the
! heating rate profile, energy is conserved.
LAPPROXLWUPDATE = .TRUE.

! If TRUE, at each gridpoint/timestep correct the shortwave net flux
! profile according to the local value of surface albedo, which at
! coasts can be different from the value used in the radiation scheme
! at coarser resolution. 
LAPPROXSWUPDATE = .TRUE.

! Implement the Manners et al. (2009) correction for the change in
! direct-beam path length through the atmosphere due to change in
! solar position with time between radiation calls.
LMANNERSSWUPDATE = .TRUE.

! The solar zenith angle "PMU0" used to compute incoming solar
! radiation is computed at timestep n+1/2, to most accurately evolve
! temperatures between timesteps n and n+1.  By default, the solar
! zenith angle "PMU0M" used in the radiation scheme to compute the
! path length of the direct beam through the atmosphere is computed a
! further NRADFR/2 timesteps ahead (where NRADFR is the frequency of
! calls to the radiation scheme). The problem with this is that for
! radiation every timestep, NRADFR=1, PMU0M is actually computed at
! the next timestep and is not equal to PMU0. Setting the following to
! TRUE corrects this behaviour.
LCENTREDTIMESZA = .TRUE.

! The solar zenith angle is computed as an average over the time
! interval in which it is used.  Thus PMU0 is computed as an average
! over model timestep n to n+1, which minimizes the variations in
! incoming solar radiation that can occur if a discrete time
! (typically n+1/2) is used. Similarly, PMU0M is computed as an
! average over the radiation timestep, which minimizes the variations
! of atmospheric solar absorption and fluxes that can occur if a
! discrete time is used.  Note that if this flag is TRUE then
! LCentredTimeSZA must also be TRUE.
LAVERAGESZA = .TRUE.

! Shortwave water vapour continuum model: 0=default MT_CKD2.5,
! 1=CAVIAR
NSWWVCONTINUUM = 0

! Each node checks the physical range of output fluxes and warns the
! first time they go out of bounds.  If this is greater than 0 then
! also write out NetCDF file of inputs and outputs to
! /home/parr/ifs_dump for offline analysis, up to NDUMPBADINPUTS
! times per node.  If this number is negative then write
! inputs/outputs then abort the first time bad fluxes are encountered.
NDUMPBADINPUTS = 0

! For debugging we can output up to this many inputs and outputs
! regardless of whether the fluxes were out of reasonable bounds.
NDUMPINPUTS = 0

! When interpolating model fields to the radiation grid, do we
! interpolate the in-cloud mean cloud and precipitation water
! contents (TRUE)?  Better conservation achieved by instead interpolating
! gridbox-mean water contents (FALSE).
LINTERPINCLOUDMEAN = .FALSE.

!- Use variable ozone instead of climatology (for CMIP runs)
NO3CMIP=0
CO3DATADIR='./'
CO3DATAFIL='not set'
GHGDATADIR='./'

! NGHGCMIP specifies CMIP era for the time evolution of well-mixed
! greenhouse gases: 0=CMIP3, 5=CMIP5, 6=CMIP6.  For CMIP3, NSCEN
! controls the scenario, while NRCP specifies the Representative
! Concentration Pathway (RCP) in CMIP5 and the equivalent Shared
! Socio-economic Pathway (SSP) in CMIP6. Later on in this file you can
! see which scenario is loaded by each NGHGCMIP/NRCP combination.

! Default up to 46R1: CMIP3 A1B scenario
!NGHGCMIP=0
!NRCP=1
NSCEN=1

! Default from 47R1: CMIP6 SSP3-7.0 scenario
IF (LEPHYS.AND..NOT.LARPEGEF) THEN ! LARPEGEF for the files conversions IFS->Arpege.
  NGHGCMIP=6
ELSE
  NGHGCMIP=0
ENDIF
NRCP=3

!- Positive values fix CMIP forcings at specified year. 
NCMIPFIXYR=-1

!-- LDIAGFORCING: to enable extra variables to store in grib files the
!-- input climatologies to radiation (GHGs and aerosols)
!-- When LDIAGFORCING=True, remember also to add the following variable to the
!-- right namelists and check that there are no conflicts with the number
!-- of extra variables associated to radiation
!-- NAMDPHY: NVEXTRRAD=15
!-- NAMPHYDS: NVEXRADGB=82,83,84,85,86,87,88,89,90,91,92,93,94,95,96
!-- NAMFPC: NFP3DFS=15, MFP3Dfs=82,83,84,85,86,87,88,89,90,91,92,93,94,95,96
LDIAGFORCING=.FALSE.

! Do we use the pre-2017 radiation scheme "McICA" by default? The
! alternative is the newer modular scheme "ecRad" in the separate
! "radiation" project.
LUSEPRE2017RAD=.NOT.LECMWF.OR.LARPEGEF ! LARPEGEF for the files conversions IFS->Arpege. REK 

! Default settings for the ecRad radiation scheme, ignored for the
! pre-2015 scheme.  First longwave scattering: at a little extra
! expense a two-stream calculation is done in the longwave, where
! 0=none, 1=cloud scattering only, 2=cloud and aerosol scattering.
NLWSCATTERING=1
! Select solver in shortwave and longwave, where 0=McICA,
! 1=SPARTACUS-1D, 2=SPARTACUS-3D, 3=Tripleclouds
NLWSOLVER = 0
NSWSOLVER = 0

! Cloud overlap scheme in ecRad: 1=maximum-random,
! 2=exponential-exponential (as McRad behaviour), 3=exponential-random
NCLOUDOVERLAP = 3

! If SPARTACUS is used then cloud scale is defined by these
! numbers. The default values were derived from ARM data by Mark
! Fielding.
RCLOUD_SEPARATION_SCALE_TOA  = 14000.0_JPRB
RCLOUD_SEPARATION_SCALE_SURF =  2500.0_JPRB

!This is read from SU0PHY in NAEPHY and put in YOEPHY

!- default setting of cloud optical properties
!  liquid water cloud 0: Fouquart    (SW), Smith-Shi   (LW)
!                     1: Slingo      (SW), Savijarvi   (LW)
!                     2: Slingo      (SW), Lindner-Li  (LW)
!                     3: Nielsen     (SW), Smith-Shi   (LW) *RADLSW only*
!                     4: SOCRATES (SW/LW) *ecRad only*
!  ice water cloud    0: Ebert-Curry (SW), Smith-Shi   (LW)
!                     1: Ebert-Curry (SW), Ebert-Curry (LW)
!                     2: Fu-Liou'93  (SW), Fu-Liou'93  (LW)
!                     3: Fu'96       (SW), Fu et al'98 (LW)
!                     4: Modified Baran et al. (2014) *ecRad only*
IF (LUSEPRE2017RAD) THEN
  NLIQOPT=2         ! before 32R1 default=0    2
ELSE
  ! ecRad radiation scheme
  NLIQOPT=4
ENDIF
NICEOPT=3           ! before 32R1 default=1    3
! LW and SW options possible to separate, but for RADLSW only
NSWLIQOPT=NLIQOPT
NLWLIQOPT=NLIQOPT
NSWICEOPT=NICEOPT
NLWICEOPT=NICEOPT
IF (LHLRADUPD) THEN
   NLWLIQOPT=3
   NSWLIQOPT=3
ENDIF

!- default setting of cloud effective radius/diameter
!  liquid water cloud 0: f(P) 10 to 45
!                     1: 13: ocean; 10: land
!                     2: Martin et al. CCN 50 over ocean, 900 over land
!  ice water cloud    0: 40 microns
!                     1: f(T) 40 to 130 microns
!                     2: f(T) 30 to 60
!                     3: f(T,IWC) Sun'01: 22.5 to 175 microns
!  conversion factor between effective radius and particle size for ice
NRADIP=3            ! before 32R1 default=2     3
NRADLP=2            ! before 32R1 default=2     2
RRE2DE=0.64952_JPRD ! before 32R1 default=0.5_JPRD

!-- changed by PBe with CY35R3
!-- RMINICE=48._JPRD
NMINICE=1 ! 0=Constant, 1=varies with latitude (larger at Equator)
RMINICE=60._JPRD
NDECOLAT=2

!- RRTM as LW scheme
LRRTM  = LECMWF
LECSRAD=.FALSE.
!- SRTM as SW scheme
LSRTM = .TRUE.     ! before 32R1 default was .FALSE.   
!-- reduced number of g-points in both RRTM_LW and _SW 
!  =0               original RRTM_LW 256 and RRTM_SW 224 g-point configurations
!  =1               ECMWF high-resolution model configuration _LW: 140, _SW: 112
!  =2               ECMWF EPS configuration _LW: 70, _SW: 56
NREDGLW=1
NREDGSW=1

! -- McICA treatment of cloud-radiation interactions 
! - 1 is maximum-random, 2 is generalized cloud overlap (before 31R1 default=0 no McICA)
NMCICA = 2          !  2 for generalized overlap

LLREDURAD=.FALSE.
IF (LSRTM .AND. NMCICA /= 0) THEN
   LLREDURAD=.TRUE.
ENDIF

IF ((LSRTM.OR.LRAY).AND.(RADGR > 0 .OR. RADSN > 0)) THEN
   WRITE(UNIT=KULOUT,&
        & FMT='('' WITH LRAY OR LSRTM, AVOID RADGR, RADSN > 0'',&
        & '' CALL ABORT'')')
   CALL ABOR1(' ABOR1 CALLED SUECRAD')
ENDIF

! Default in-cloud water content fractional standard deviation
RCLOUD_FRAC_STD = 1.0_JPRD

! Maike's bug fix attempt - set missing default for this variable
LECSRAD=.FALSE.

! Original Fu ice optics in RADLSWR coded with a bug such that
! longwave single scattering albedo is one minus the correct
! value. Set to TRUE to leave this bug in.
LFU_LW_ICE_OPTICS_BUG = .FALSE.

!- Inhomogeneity factors in LW and SW (0=F, 1=0.7 in both, 2=Barker's, 3=Cairns)
NINHOM = 0          ! before 32R1 default=1
NLAYINH= 0
RLWINHF = 1.0_JPRD  ! before 32R1 default=0.7
RSWINHF = 1.0_JPRD  ! before 32R1 default=0.7  

!- Diffusivity correction a la Savijarvi
LDIFFC = .FALSE.    ! before 32R1 default=.FALSE. no change

!- monthly climatol. of tropospheric aerosols from Tegen et al. (1997)
LNEWAER=.TRUE.
LNOTROAER=.FALSE.
NPERTAER=0
LESO4HIS=.FALSE.
CLISTSO4=""
!- monthly climatologies of CO2, CH4 and O3 derived from GEMS-analyses
!- ABozzo: as from 2014 CY40R3 this switch enables the climatologies of
!- CO2,CH4,O3 derived from 2003-2011 MACC reanalyses. See in su_ghgclim.F90.
!- to use the old GEMS climatologies, uncomment the relevant bits in
!- su_ghgclim.F90
LETRACGMS=.TRUE.
!- monthly climatology of aerosols (SS, DU, OM, BC, SO4) from CAMS reanalyses (former MACC)
IF (LEPHYS.AND..NOT.LARPEGEF) THEN ! LARPEGEF for the files conversions IFS->Arpege. REK
  NAERMACC=1    ! =0 inactive
ELSE
  NAERMACC=0
ENDIF
LAER3D=.TRUE. !T=3d CAMS aerosol climatology


NMCLAT=1
NMCLON=1
NMCLEV=1
NMCVAR=1
LAERADCLI=.FALSE.
LAERADJDU=.FALSE.
NAERADJDU=0
RDUMULF =1.0_JPRD
RWGHTDU1=1.0_JPRD
RWGHTDU2=1.0_JPRD
RWGHTDU3=1.0_JPRD

!-- factors for scale height definition of MACC-derived aerosol climatology

! to have a monthly-varying scale height 
! for dust 
LDUSEASON=.TRUE.

RAESHSS=1000._JPRD
RAESHDU=3000._JPRD !redefined in suecaec when LDUSEASON=T.
RAESHOM=2000._JPRD
RAESHBC=1000._JPRD
RAESHSU=4000._JPRD
!-- factor to define the tropospheric background. 
!-- Default in the Tegen clim is 0.03
TRBKG=0.0534_JPRD
STBKG=0.0045_JPRD

!-- spectral map for aerosol properties
NSPMAPL = (/ 1, 1, 2, 2, 2, 3, 4, 3, 3, 5, 5, 5, 5, 5, 5, 5 /)
NSPMAPS = (/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14 /)

!- New Rayleigh formulation
LRAYL=.TRUE.

!- perpetual run for radiation
LPERPET=.FALSE.
IF( LAPE ) THEN
  LPERPET=.TRUE.
ENDIF

!- Number concentration of aerosols if specified
LCCNL=.TRUE.        ! before 32R1 default=.FALSE.     true
LCCNO=.TRUE.        ! before 32R1 default=.FALSE.     true
RCCNLND=900._JPRD   ! before 32R1 default=900. now irrelevant
RCCNSEA=50._JPRD    ! before 32R1 default=50.  now irrelevant
LAERCLIM=.FALSE.    ! to save the optical depths of the 6 types of climatological aerosols
LAERVISI=.TRUE.     ! to save the visibility using climatological or prognostic aerosols (if active)

! Constants used for longwave transmission functions related to CH4 and N2O
RCH4A = 4._JPRD/0.103_JPRD
RCH4B = 4._JPRD/0.113_JPRD
RCN2OA = 4._JPRD/0.416_JPRD
RCN2OB = 4._JPRD/0.197_JPRD

!- history of volcanic aerosols
LHVOLCA=.FALSE.
LVOLCDATA=.FALSE.
CVOLCDATA=""
NVOLCVERT=0
LVOLCSPEC=.FALSE.
LVOLCDAMP=.FALSE.
RVOLCSPEC(:)=(/STBKG,STBKG,STBKG/)

! Interaction radiation / prognostic O3 on by default, provided that
! IFS ozone is present as a GFL field; LEPO3RA can be overridden from
! the namelist.
IF (YO3%LGP) THEN
  LEPO3RA=.TRUE.
ELSE
  LEPO3RA=.FALSE.
ENDIF
RPERTOZ=0._JPRD
NPERTOZ=0

! NGHGRAD encodes which gases are climatological - see
! upd_ghgclim.F90, radghg.F90 and radintg.F90. Unfortunately it must
! be consistent with LEPO3RA (whether prognostic ozone is interactive
! with radiation): normally LEPO3RA=false should be accompanied by
! NGHGRAD=21, while LEPO3RA=true should be accompanied by
! NGHGRAD=20. Therefore here we set it to an impossible value, and if
! that is not overridden by the namelist, we set it to a value
! consistent with LEPO3RA lower down in this file.
NGHGRAD=-1

NAER=1
IF( LAPE ) THEN
  NAER=0
ENDIF
NMODE  =0
NOZOCL =1
NRADPFR=0
NRADPLA=15

!-- Frequency of calls to full radiation 
! Negative NRADFR is interpretted as hours, positive numbers as
! timesteps. NRADELG is the number of hours to use a frequency of
! radiation every hour, before switching to the frequency specified in
! NRADFR. This was useful when the default radiation frequency was 3
! hourly. Since operational cycle 46R1, HRES, ENS and data
! assimilation trajectories use 1 hour radiation, and it is also used
! in ERA5, as well as internal low-resolution "climate" testing used
! at ECMWF. Therefore the default NRADFR has been changed to one hour
! at all resolutions.
NRADFR =-1
NRADELG= 0
! Previously the default frequency was 3 hours except if resolution
! was above a certain amount.
!NRADFR=-3
!IF ( REAL(NSMAX,JPRD)/RPLRADI >= 511._JPRD ) NRADFR =-1


! AB dec 2014 -- UV diagnostic of surface fluxes over the 280-400 nm interval 
!    with up-to 24 values (5 nm wide spectral intervals)
LUVPROC=.FALSE.
LUVTDEP=.TRUE.
LUVDBG =.FALSE.
KMODTS = 0 !AB default =0 the original Morcrette&Arola processor
!          = 1 eddington (joseph et al., 1976)
!          = 2 pifm (zdunkowski et al., 1980)
!          = 3 discrete ordinates (liou, 1973)  
LUVAERP= LAERUVP    !AB prognostic aerosols in UV processor
!like that LUVAERP is not affected by LAERUVP because this is set
!after. I added LUVAERP to the naerad namelist (read later on here) 
!to control its value externally

!to enable O3 mmr from chemistry in the UV
!if FALSE the prognostic O3 comes from the 
!Cariolle parameterization
!if TRUE comes from the full chemistry (if active)
LO3_CHEM_UV= .FALSE. 

NRADUV =-1
NUVTIM = 0
NUV    = 72 !default value. It causes the crash of the whole IFS 
            !if the processor is invoked by accident. NUV can be only
            !24,120,600
RMUZUV = 1.E-03_JPRD
DO JUV=1,NUV
  RUVLAM(JUV)=280._JPRD+(JUV-1)*5._JPRD
ENDDO

!- radiation interpolation (George M's grid on by default)
LLDEBUG=.FALSE.
LEDBUG=.FALSE.
NRADINT=3
NRADRES=0

!LECOMPGRID=.TRUE.
LECOMPGRID=.FALSE.
CRTABLEDIR='./'
CRTABLEFIL='not set'
YDRI%LSLONDEM=.TRUE.
YDRO%LSLONDEM=.TRUE.
!GM Temporary as per trans/external/setup_trans.F90
IF( LLDEBUG )THEN
  WRITE(NULOUT,'("SUECRAD: NSMAX=",I6)')NSMAX
  WRITE(NULOUT,'("SUECRAD: NDLON=",I6)')NDLON
ENDIF

NUAER  = 24
NTRAER = 15
! Cloud overlap in pre-McRad scheme (1=maximum-random)
NOVLP  = 1
IF ( LRAY ) THEN
  NTSW = 1
  NSW  = 1
ELSE
  NTSW = 14
  NSW  = 6
ENDIF
IF ( LHLRADUPD ) THEN
  NTSW = 14
  NSW  = 6
ENDIF
NSWNL  = 6
NSWTL  = 2

IF(LOPT_SCALAR) THEN
  NRPROMA=-8
ELSE
  NRPROMA=511
ENDIF
IRPROMA=NRPROMA

!*         2.3      SET SECURITY PARAMETERS
!                   -----------------------

REPSC  = 1.E-04_JPRD
REPSCA = 1.E-10_JPRD
REPSCO = 1.E-12_JPRD
REPSCQ = 1.E-12_JPRD
REPSCT = 1.E-12_JPRD
REPSCW = 1.E-12_JPRD
REPLOG = 1.E-12_JPRD


!*          2.4     BACKGROUND GAS CONCENTRATIONS (IPCC/SACC, 1990)
!                   -----------------------------------------------

! Allow for time variation of greenhouse gas concentrations?
LHGHG   = LEPHYS .AND. .NOT. LARPEGEF  ! LARPEGEF for FullPos IFS->Arp

! Completely superseded by LHGHG: ignored by subsequent code and only
! here so that namelists containing it do not fail
LECO2VAR=.TRUE.

! Set the default treatment of total solar irradiance TSI (can be
! overridden by namelist value).  NHINCSOL=0 uses 1366 W m-2.
! NHINCSOL=4 uses the CMIP6 values from NetCDF file.  Override by
! using CSOLARIRRADIANCEFILE to specify your own NetCDF file.
IF (LEPHYS.AND..NOT.LARPEGEF) THEN ! LARPEGEF for the files conversions IFS->Arpege. REK
  NHINCSOL = 4
ELSE
  NHINCSOL = 0
ENDIF

! Set the treatment of the solar spectrum: 0 is to use the original
! from RRTM, which is Kurucz; 1 is to scale this to match the more
! recent measurements from the SORCE instrument
NSOLARSPECTRUM = 0

!- CMIP6 total solar irradiance specified via file when NHINCSOL=4.
CSOLARDATA=''

RSOLINC = RI0

! Reference gas concentrations from IPCC 1990, unlikely to be used in
! any configuration (see surdi.F90 and updrgas.F90)
RCCO2   = 353.E-06_JPRD
RCCH4   = 1.72E-06_JPRD
RCN2O   = 310.E-09_JPRD
RCNO2   = 500.E-13_JPRD
RCCFC11 = 280.E-12_JPRD
RCCFC12 = 484.E-12_JPRD
RCCFC22 =   0.E-12_JPRD
RCCCL4  =   0.E-12_JPRD
!RCCFC22 = 170.E-12_JPRD
!RCCCL4  = 883.E-13_JPRD

IF (LAPE) THEN
  ! Aquaplanet standard values
  RCCO2   = 348.E-06_JPRD
  RCCH4   = 1.65E-06_JPRD
  RCN2O   = 306.E-09_JPRD
ENDIF


!*          2.5     PREPARATORY WORK FOR VISIBILITY CALCULATIONS
!                   --------------------------------------------

ZCCO2  =RCCO2 * 100._JPRD
ZLAMBDA=0.55_JPRD
Z1LMBD2=1._JPRD/(ZLAMBDA*ZLAMBDA)
Z1LMBD4=Z1LMBD2*Z1LMBD2
ZLAMBDACM=ZLAMBDA*1.E-04_JPRD
ZLMBDCM4 =ZLAMBDACM**4

!-- depolarization factors due to various gases
ZFN2   =1.034_JPRD + 3.17E-04_JPRD * Z1LMBD2
ZFO2   =1.096_JPRD + 1.385E-03_JPRD * Z1LMBD2 + 1.448E-04_JPRD * Z1LMBD4
ZNFAIR =78.084_JPRD * ZFN2 + 20.946_JPRD * ZFO2 + 0.934_JPRD + ZCCO2 * 1.15_JPRD
ZDFAIR =78.084_JPRD        + 20.496_JPRD        + 0.934_JPRD + ZCCO2
ZFAIR  =ZNFAIR/ZDFAIR
!-- molecular density (molecules cm-3) at T=273.15
RNS    =RNAVO / (22.414_JPRD * 1000._JPRD)
!-- refractive index of air for 300 ppm of CO2
ZRIM1_300 = 1E-08_JPRD*(8060.51_JPRD + 2480990._JPRD/(132.274_JPRD - Z1LMBD2)&
 &            + 17455.7_JPRD/(39.32957_JPRD - Z1LMBD2))
!-- applying conversion factor for refractive index of air at zcCO2 /= 300 ppm
ZRIAIRM1 = ZRIM1_300 * (1._JPRD + 0.54_JPRD * (RCCO2 - 0.0003_JPRD))
!-- refractive index of air
ZRIAIR   = 1._JPRD + ZRIAIRM1
ZRIAIR2  = ZRIAIR*ZRIAIR
!-- invariable bits of the scattering cross section at given zlambda 
!  (final one to be in cm-2 molec-1)  
RSIGAIR = 24._JPRD*RPI**3 * (ZRIAIR2 -1._JPRD)**2 * ZFAIR&
  &     / (ZLMBDCM4 * (ZRIAIR2 +2._JPRD)**2)

!     ------------------------------------------------------------------

!*         3.       READ VALUES OF RADIATION CONFIGURATION
!                   --------------------------------------

CALL POSNAM(NULNAM,'NAERAD')
READ (NULNAM,NAERAD)

IF (NRPROMA==0) NRPROMA=IRPROMA  ! restore default nrproma

!- reset some defaults if SW6 is used - AROME default
!(revert to pre-CY3?R1 operational configuration)
IF (.NOT.LSRTM) THEN
  NMCICA = 0
  LCCNL  = .FALSE.
  LCCNO  = .FALSE.
  LDIFFC = .FALSE.
  IF (.NOT.LHLRADUPD) THEN   ! MF operational defaults
    NICEOPT= 1
    NLIQOPT= 0
    NSWLIQOPT=NLIQOPT
    NLWLIQOPT=NLIQOPT
    NSWICEOPT=NICEOPT
    NLWICEOPT=NICEOPT
    NRADIP = 2
    NRADLP = 2
    RRE2DE = 0.5_JPRD
    NINHOM = 1
    RLWINHF= 0.7_JPRD
    RSWINHF= 0.7_JPRD
  ENDIF
ENDIF

!- if MACC-derived aerosol climatology is used, potential adjustments
IF (LAERADJDU) THEN
  CALL POSNAM(NULNAM,'NAERCLI')
  READ(NULNAM,NAERCLI)
ENDIF

! New read of namelist to modify default values when LSRTM=.FALSE. 
! (Meteo-France operational model)
CALL POSNAM(NULNAM,'NAERAD')
READ (NULNAM,NAERAD)

IF (NGHGRAD < 0) THEN
  IF ((LEPHYS .AND. .NOT. LARPEGEF) .OR. LRAYFM) THEN
    ! NGHGRAD has not been set in the namelist so we set it consistent
    ! with LEPO3RA (=prognostic ozone interactive with radiation)
    IF (LEPO3RA) THEN
      ! Use 3D climatologies for CO2, CH4, N2O, but not O3
      NGHGRAD = 20
    ELSE
      ! Use 3D climatologies for CO2, CH4, N2O and O3
      NGHGRAD = 21
    ENDIF
  ELSE
    NGHGRAD = 0
  ENDIF
ENDIF

! The total solar irradiance file may be specified in the namelist; if
! it starts with a '/' or '.' then the path is used directly,
! otherwise the directory in the DATA environment variable is
! searched.  If not specified in the namelist then set to the default
! value.
IF (CSOLARIRRADIANCEFILE == ' ') THEN
  CSOLARIRRADIANCEFILE = 'total_solar_irradiance_CMIP6_47r1.nc';
ENDIF

! The greenhouse gas monthly climatology may be specified in the
! namelist; if it starts with a '/' or '.' then the path is used
! directly, otherwise the directory in the DATA environment variable
! is searched.  If the greenhouse gas climatology file has not been
! specified in the namelist, set to the default value.
IF (CGHGCLIMFILE == ' ') THEN
  ! 91-level climatology used until 47r3
  !CGHGCLIMFILE = 'greenhouse_gas_climatology_46r1.nc'
  ! Identical for all gases except ozone, which has been updated to
  ! match the CAMS Reanalysis, merged with CMIP6 data in the
  ! mesosphere.  Note that in 48r1 we use radiatively interactive
  ! prognostic ozone, so this new ozone climatology is only seen if
  ! interactive ozone is turned off (LEPO3RA=false).
  CGHGCLIMFILE = 'greenhouse_gas_climatology_48r1.nc'
ENDIF

! The climatology is scaled by the surface value appropriate for the
! year. This again may be specified directly in the namelist;
! otherwise the appropriate file is determined from other
! configuration variables.
IF (CGHGTIMESERIESFILE == ' ') THEN
  IF (NGHGCMIP == 0) THEN
    ! CMIP3 forcings, which were hard-coded in updrgas.F90 up to and
    ! including cycle 46r1. The CFC11 concentrations appear to be an
    ! "equivalent" value including additional minor species.  It is
    ! therefore possible that HCFC22 and CCl4 (constant with time) are
    ! double counted.
    IF (NSCEN == 1) THEN
      ! This option was the default up to and including operational
      ! IFS cycle 46r1.
      CGHGTIMESERIESFILE = 'greenhouse_gas_timeseries_CMIP3_A1B_46r1.nc'
    ELSEIF (NSCEN == 2) THEN
      CGHGTIMESERIESFILE = 'greenhouse_gas_timeseries_CMIP3_A2_46r1.nc'
    ELSEIF (NSCEN == 3) THEN
      CGHGTIMESERIESFILE = 'greenhouse_gas_timeseries_CMIP3_B1_46r1.nc'
    ELSE
      CALL ABOR1('SUECRAD: if NGHGCMIP==0 then NSCEN must be between 1 and 3')
    ENDIF
  ELSEIF (NGHGCMIP == 5) THEN
    ! CMIP5 forcings, which were available to be read in from text
    ! files up to and including cycle 46r1. These files replicate
    ! those values. Note that CFC11 is not an "equivalent" value, and
    ! HCFC22 and CCl4 are constant with time, so the radiative forcing
    ! associated with more minor species is not accounted for.
    IF (NRCP == 1) THEN
      CGHGTIMESERIESFILE = 'greenhouse_gas_timeseries_CMIP5_RCP3PD_46r1.nc'
    ELSEIF (NRCP == 2) THEN
      CGHGTIMESERIESFILE = 'greenhouse_gas_timeseries_CMIP5_RCP45_46r1.nc'
    ELSEIF (NRCP == 3) THEN
      CGHGTIMESERIESFILE = 'greenhouse_gas_timeseries_CMIP5_RCP6_46r1.nc'
    ELSEIF (NRCP == 4) THEN
      CGHGTIMESERIESFILE = 'greenhouse_gas_timeseries_CMIP5_RCP85_46r1.nc'
    ELSE
      CALL ABOR1('SUECRAD: if NGHGCMIP==5 then NRCP must be between 1 and 4')
    ENDIF
  ELSEIF (NGHGCMIP == 6) THEN
    ! CMIP6 forcings. These files use an equivalent CFC11
    ! concentration scaled to account for all other minor greenhouse
    ! gas species.  Therefore HCFC22 and CCl4 have their
    ! concentrations set to zero.
    IF (NRCP == 1) THEN
      CGHGTIMESERIESFILE = 'greenhouse_gas_timeseries_CMIP6_SSP126_CFC11equiv_47r1.nc'
    ELSEIF (NRCP == 2) THEN
      CGHGTIMESERIESFILE = 'greenhouse_gas_timeseries_CMIP6_SSP245_CFC11equiv_47r1.nc'
    ELSEIF (NRCP == 3) THEN
      CGHGTIMESERIESFILE = 'greenhouse_gas_timeseries_CMIP6_SSP370_CFC11equiv_47r1.nc'
    ELSEIF (NRCP == 4) THEN
      CGHGTIMESERIESFILE = 'greenhouse_gas_timeseries_CMIP6_SSP585_CFC11equiv_47r1.nc'
    ELSE
      CALL ABOR1('SUECRAD: if NGHGCMIP==6 then NRCP must be between 1 and 4')
    ENDIF
  ELSE
    CALL ABOR1('SUECRAD: NGHGCMIP must be 0, 5 or 6')
  ENDIF
ENDIF

IF (NRPROMA==0) NRPROMA=IRPROMA  ! restore default nrproma

IF (LHLRADUPD) THEN
   IF (NRADIP /= 3) RRE2DE = 0.5_JPRD
   ! spherical crystals for non-default NRADIP
   IF (NINHOM == 1 .AND. RSWINHF == 1.0_JPRD .AND. RLWINHF == 1.0_JPRD) THEN
   ! old default LW/SW inhomogenity values if only NINHOM=1 requested
      RLWINHF=0.7_JPRD
      RSWINHF=0.7_JPRD
   ENDIF
ENDIF
IF (NINHOM == 0) THEN
   RLWINHF=1.0_JPRD
   RSWINHF=1.0_JPRD
ENDIF

! Number of longwave surface emissivity intervals to use
IF (YDEPHY%NEMISSSCHEME == 1) THEN
  ! We do a more accurate mapping for emissivity if NEMISSSCHEME==1.
  ! See susrad_mod.F90 for values for different surface types.
  YDERAD%NLWEMISS = 6
  IF (YDERAD%LAPPROXLWUPDATE) THEN
    ! Pass the same number of longwave downwelling surface spectral
    ! fluxes from ecRad to RADHEATN so that longwave approximate
    ! update scheme can be as accurate as possible
    YDERAD%NLWOUT = 6
  ELSE
    YDERAD%NLWOUT = 1
  ENDIF
  ! Create a spectral Planck look-up table, used by RADHEATN.  Note
  ! that this routine makes use of the length of its third argument.
  ! The following wavelength bounds (metres) match the RRTM band
  ! boundaries.
  CALL YDERAD%YSPECTPLANCK%INIT(6, &
       &  [ 8.4746E-6_JPRB, 10.2041E-6_JPRB, 12.1951E-6_JPRB, 15.8730E-6_JPRB, 28.5714E-6_JPRB ], &
       &  [ 1,2,3,4,5,6 ])
ELSEIF (YDEPHY%NEMISSSCHEME == 0) THEN
  ! Traditional approach: one value of emissivty for parts of the
  ! spectrum on either side of the infrared atmospheric window
  ! (PEMIR), and one value for the window itself (PEMIW)
  YDERAD%NLWEMISS = 2
  ! ...and the longwave approximate update scheme uses a single
  ! broadband emissivity
  YDERAD%NLWOUT   = 1
  ! Create a spectral Planck look-up table, used by RADHEATN.  Note
  ! that this routine makes use of the length of its third argument.
  ! The wavelength bounds (metres) allow for the first emissivity to
  ! represent values outside the infrared atmospheric window, and the
  ! second emissivity to represent values within it.
  CALL YDERAD%YSPECTPLANCK%INIT(2, [ 8.0E-6_JPRB, 13.0E-6_JPRB ], [ 1,2,1 ])
ELSE
  CALL ABOR1('RADIATION_SETUP: NEMISSSCHEME must be 0 or 1')
ENDIF


!- for McICA computations, make sure these parameters are as follows ...
IF (NMCICA /= 0) THEN
   NINHOM = 0
   RLWINHF= 1.0_JPRD
   RSWINHF= 1.0_JPRD
   IF (LUSEPRE2017RAD) THEN
!-- read the XCW values for Raisanen-Cole-Barker cloud generator
     CALL SU_MCICA(YDMODEL%YRML_PHY_RAD%YREMCICA)
   ENDIF
ENDIF

IF( LLDEBUG )THEN
  WRITE(NULOUT,'("SUECRAD: NRADINT=",I2)')NRADINT
  WRITE(NULOUT,'("SUECRAD: NRADRES=",I4)')NRADRES
ENDIF

!     DETERMINE WHETHER NRPROMA IS NEGATIVE AND SET LOPTRPROMA

LOPTRPROMA=NRPROMA > 0
NRPROMA=ABS(NRPROMA)

IF( NRADINT > 0 .AND. NRADRES == NSMAX )THEN
  WRITE(NULOUT,'("SUECRAD: NRADINT > 0 .AND. NRADRES = NSMAX, NRADINT RESET TO 0")')
  NRADINT=0
ENDIF

IF( NRADINT > 0 .AND. LRAYFM .AND. NAER /= 0 .AND. .NOT.LHVOLCA )THEN
!   This combination is not supported as aerosol data would be
!   required to be interpolated (see radintg)
  WRITE(NULOUT,'("SUECRAD: NRADINT>0, LRAYFM=T NAER /= 0 .AND. LHVOLCA=F,",&
   & " NRADRES RESET TO NSMAX (NO INTERPOLATION)")')  
  NRADRES=NSMAX
ENDIF
CALL GSTATS(1818,1)

100 CONTINUE

! Nullify some pointer for a preper deallocation at the end of the model:
! NULLIFY(RADGRID%NRGRI)
! NULLIFY(RADGRID%NLOENG)
! NULLIFY(RADGRID%NSTA)
! NULLIFY(RADGRID%NONL)
! NULLIFY(RADGRID%NPTRFRSTLAT)
! NULLIFY(RADGRID%NFRSTLAT)
! NULLIFY(RADGRID%NLSTLAT)
! NULLIFY(RADGRID%RMU)
! NULLIFY(RADGRID%RSQM2)
! NULLIFY(RADGRID%RLATIG)
! NULLIFY(RADGRID%NASM0)
! NULLIFY(RADGRID%MYMS)
! NULLIFY(RADGRID%GELAM)
! NULLIFY(RADGRID%GELAT)
! NULLIFY(RADGRID%GESLO)
! NULLIFY(RADGRID%GECLO)
! NULLIFY(RADGRID%GEMU)
! NULLIFY(RADGRID%RLATI)
! NULLIFY(RADGRID%RIPI)


IF( LERADI )THEN   ! START OF LERADI BLOCK

  IF( NRADINT == 0 )THEN

    IF( NRADRES /= NSMAX )THEN
      WRITE(NULOUT,'("SUECRAD: NRADINT=0 REQUESTED, NRADRES RESET TO NSMAX")')
      NRADRES=NSMAX
    ENDIF
    RADGRID%NGPTOT=NGPTOT
    RADGRID%NGPTOTMX=NGPTOTMX

    YDRI%NASLB1=0
    YDRO%NASLB1=0

  ELSEIF( NRADINT >=1 .AND. NRADINT <= 3 )THEN

    LLGRIDONLY = .NOT.NRADINT==1
    YDRI%NASLB1=0
    YDRO%NASLB1=0

! set the default radiation grid resolution for the current model resolution
! if not already specified, always based on linear grid
    LLINEAR_GRID = .TRUE.
    IF( NRADRES == 0 ) THEN
    IGG=NDGLG/2
    IF( IGG == 16 .OR.  IGG == 32 )  THEN
      NRADRES = 21
      LLINEAR_GRID = .FALSE.
    ENDIF
    IF( IGG == 48 )  NRADRES =  63   ! 
    IF( IGG == 64 )  NRADRES =  63   ! 
    IF( IGG == 80 )  NRADRES=   63   ! 5.84
    IF( IGG == 96 )  NRADRES =  63   ! 
    IF( IGG == 128 ) NRADRES = 95   ! 6.69
    IF( IGG == 160 ) NRADRES=  159   ! 3.87
    IF( IGG == 200 ) NRADRES=  159   !
    IF( IGG ==  256 ) NRADRES=  159   !
    IF( IGG ==  320 ) NRADRES=  255   !
    IF( IGG ==  400 ) NRADRES=  319   !
    IF( IGG == 512 ) NRADRES=  399   !
    IF( IGG == 640 ) NRADRES=  511   ! 
    IF( IGG == 800) NRADRES=  799   !
    IF( IGG == 912) NRADRES=  799   ! 
    IF( IGG == 1024) NRADRES=  799   ! 
    IF( IGG == 1280) NRADRES=  799   ! 
    IF( IGG == 1600) NRADRES=  1023   ! 
    IF( IGG == 2000) NRADRES=  1023   ! 
    IF( IGG == 2560) NRADRES=  2047  !
    IF( IGG == 4000) NRADRES=  2047  ! 
    IF( IGG == 8000) NRADRES=  2047  ! 
    ENDIF

! test if radiation grid resolution has been set
    IF( NRADRES == 0 )THEN
      WRITE(NULOUT,'("SUECRAD: NRADRES NOT SET OR DEFAULT FOUND,NSMAX=",I4)')NSMAX
      CALL ABOR1('SUECRAD: NRADRES NOT SET OR DEFAULT FOUND')
    ENDIF


! test if no interpolation is required
    IF( NRADINT > 0 .AND. NRADRES == NSMAX )THEN
      WRITE(NULOUT,'("SUECRAD: NRADINT > 0 .AND. NRADRES = NSMAX, NRADINT RESET TO 0")')
      NRADINT=0
      GOTO 100
    ENDIF

    IF( LECOMPGRID ) THEN
! compute rad grid
      WRITE(NULOUT,'("COMPUTING OCTAHEDRAL RADIATION GRID")')
      RADGRID%NDGLG=NRADRES+1
      IDGLG2=(RADGRID%NDGLG+1)/2
      ALLOCATE(RADGRID%NRGRI(RADGRID%NDGLG))
      DO IDGL=1,IDGLG2
        RADGRID%NRGRI(IDGL)=20+(IDGL-1)*4
        RADGRID%NRGRI(RADGRID%NDGLG-IDGL+1)=20+(IDGL-1)*4
      ENDDO
      RADGRID%NSMAX=0
    ELSE

    CALL GSTATS(1818,0)
    IF( CRTABLEFIL == 'not set' )THEN
      IF( LLINEAR_GRID )THEN
        IF( NRADRES < 1000 )THEN
          WRITE(CRTABLEFIL,'("rtablel_2",I3.3)')NRADRES
        ELSE
          WRITE(CRTABLEFIL,'("rtablel_2",I4.4)')NRADRES
        ENDIF
      ELSE
        IF( NRADRES < 1000 )THEN
          WRITE(CRTABLEFIL,'("rtable_2" ,I3.3)')NRADRES
        ELSE
          WRITE(CRTABLEFIL,'("rtable_2" ,I4.4)')NRADRES
        ENDIF
      ENDIF
    ENDIF
    CALL GSTATS(1818,1)

    RADGRID%NSMAX=NRADRES

    IF( MYPROC == JPIOMASTER )THEN
      IDIR=LEN_TRIM(CRTABLEDIR)
      IFIL=LEN_TRIM(CRTABLEFIL)
      CLF1=CRTABLEDIR(1:IDIR)//CRTABLEFIL(1:IFIL)
      CALL GET_ENVIRONMENT_VARIABLE("DATA",CLZZZ)
      IF(CLZZZ /= " ".AND.CRTABLEDIR=='./') THEN
        WRITE(0,'(1x,a)')'Reading '//TRIM(CLF1)
        CLF1=TRIM(CLZZZ)//"/ifs/"//CRTABLEDIR(1:IDIR)//CRTABLEFIL(1:IFIL)
      ENDIF
      OPEN(NULRAD,FILE=CLF1,ACTION="READ",STATUS="OLD",ERR=999)
      GOTO 1000
      999 CONTINUE
      WRITE(NULOUT,'("SUECRAD: UNABLE TO OPEN FILE ",A)') trim(CLF1)
      CALL ABOR1('SUECRAD: UNABLE TO OPEN RADIATION GRID RTABLE FILE '//trim(CLF1))
      1000 CONTINUE
      NRGRI(:)=0
      CALL POSNAM(NULRAD,'NAMRGRI')
      READ (NULRAD,NAMRGRI)
      IDGL=1
      DO WHILE( NRGRI(IDGL)>0 )
        IF( LLDEBUG )THEN
          WRITE(NULOUT,'("SUECRAD: NRGRI(",I4,")=",I4)')IDGL,NRGRI(IDGL)
        ENDIF
        IDGL=IDGL+1
      ENDDO
      IDGL=IDGL-1
      RADGRID%NDGLG=IDGL
      IF( LLDEBUG )THEN
        WRITE(NULOUT,'("SUECRAD: RADGRID%NDGLG=",I4)')RADGRID%NDGLG
      ENDIF
      CLOSE(NULRAD)
    ENDIF
    CALL GSTATS(667,0)
    IF( NPROC > 1 )THEN
      CALL MPL_BROADCAST (RADGRID%NDGLG,MTAGRAD,JPIOMASTER,CDSTRING='SUECRAD:')
    ENDIF
    ALLOCATE(RADGRID%NRGRI(RADGRID%NDGLG))
    IF( MYPROC == JPIOMASTER )THEN
      RADGRID%NRGRI(1:RADGRID%NDGLG)=NRGRI(1:RADGRID%NDGLG)
    ENDIF
    IF( NPROC > 1 )THEN
      CALL MPL_BROADCAST (RADGRID%NRGRI(1:RADGRID%NDGLG),MTAGRAD,JPIOMASTER,CDSTRING='SUECRAD:')
    ENDIF
    CALL GSTATS(667,1)

    ENDIF

    CALL GSTATS(1818,0)
    IF    ( NRADINT == 1 )THEN
      WRITE(NULOUT,'("SUECRAD: INTERPOLATION METHOD - SPECTRAL TRANSFORM")')
      RADGRID%NDGSUR=0
      YDRI%NSLWIDEN=0
      YDRI%NSLWIDES=0
      YDRI%NSLWIDEW=0
      YDRI%NSLWIDEE=0
      YDRO%NSLWIDEN=0
      YDRO%NSLWIDES=0
      YDRO%NSLWIDEW=0
      YDRO%NSLWIDEE=0
    ELSEIF( NRADINT == 2 )THEN
      WRITE(NULOUT,'("SUECRAD: INTERPOLATION METHOD - 4 POINT")')
      RADGRID%NDGSUR=2
    ELSEIF( NRADINT == 3 )THEN
      WRITE(NULOUT,'("SUECRAD: INTERPOLATION METHOD - 12 POINT")')
      RADGRID%NDGSUR=2
    ENDIF
    WRITE(NULOUT,'("SUECRAD: RADGRID%NDGSUR       =",I8)')RADGRID%NDGSUR

    RADGRID%NDGSAG=1-RADGRID%NDGSUR
    RADGRID%NDGENG=RADGRID%NDGLG+RADGRID%NDGSUR
    RADGRID%NDLON=RADGRID%NRGRI(RADGRID%NDGLG/2)
    WRITE(NULOUT,'("SUECRAD: RADGRID%NDGSAG       =",I8)')RADGRID%NDGSAG
    WRITE(NULOUT,'("SUECRAD: RADGRID%NDGENG       =",I8)')RADGRID%NDGENG
    WRITE(NULOUT,'("SUECRAD: RADGRID%NDGLG        =",I8)')RADGRID%NDGLG
    WRITE(NULOUT,'("SUECRAD: RADGRID%NDLON        =",I8)')RADGRID%NDLON
    CALL FLUSH(NULOUT)

    ALLOCATE(RADGRID%NLOENG(RADGRID%NDGSAG:RADGRID%NDGENG))
    RADGRID%NLOENG(1:RADGRID%NDGLG)=RADGRID%NRGRI(1:RADGRID%NDGLG)
    IF(RADGRID%NDGSUR >= 1)THEN
      DO JGLSUR=1,RADGRID%NDGSUR
        RADGRID%NLOENG(1-JGLSUR)=RADGRID%NLOENG(JGLSUR)
      ENDDO
      DO JGLSUR=1,RADGRID%NDGSUR
        RADGRID%NLOENG(RADGRID%NDGLG+JGLSUR)=RADGRID%NLOENG(RADGRID%NDGLG+1-JGLSUR)
      ENDDO
    ENDIF
    CALL GSTATS(1818,1)

    WRITE(NULOUT,'("SUECRAD: NRADRES =",I8)')NRADRES
    WRITE(NULOUT,'("SUECRAD: NRADINT =",I8)')NRADINT
    WRITE(NULOUT,'("SUECRAD: RADGRID%NSMAX =",I8)')RADGRID%NSMAX
    WRITE(NULOUT,'("SUECRAD: LLINEAR_GRID=",L8)')LLINEAR_GRID
    WRITE(NULOUT,'("SUECRAD: LLGRIDONLY=",L8)')LLGRIDONLY
! Setup the transform package for the radiation grid
    CALL SETUP_TRANS (KSMAX=RADGRID%NSMAX,&
     & KDGL=RADGRID%NDGLG,&
     & KLOEN=RADGRID%NLOENG(1:RADGRID%NDGLG),&
     & LDSPLIT=LSPLIT,LDUSEFFTW=LFFTW,&
     & KRESOL=RADGRID%NRESOL_ID,&
     & LDGRIDONLY=LLGRIDONLY)

    ALLOCATE(RADGRID%NSTA(RADGRID%NDGSAG:RADGRID%NDGENG+N_REGIONS_NS-1,N_REGIONS_EW))
    ALLOCATE(RADGRID%NONL(RADGRID%NDGSAG:RADGRID%NDGENG+N_REGIONS_NS-1,N_REGIONS_EW))
    ALLOCATE(RADGRID%NPTRFRSTLAT(N_REGIONS_NS))
    ALLOCATE(RADGRID%NFRSTLAT(N_REGIONS_NS))
    ALLOCATE(RADGRID%NLSTLAT(N_REGIONS_NS))
    ALLOCATE(RADGRID%RMU(RADGRID%NDGSAG:RADGRID%NDGENG))
    ALLOCATE(ZMU(RADGRID%NDGSAG:RADGRID%NDGENG))
    ALLOCATE(RADGRID%RSQM2(RADGRID%NDGSAG:RADGRID%NDGENG))
    ALLOCATE(RADGRID%RLATIG(RADGRID%NDGSAG:RADGRID%NDGENG))

! Interrogate the transform package for the radiation grid
    CALL GSTATS(1818,0)
    CALL TRANS_INQ (KRESOL     =RADGRID%NRESOL_ID,&
     & KSPEC2     =RADGRID%NSPEC2,&
     & KNUMP      =RADGRID%NUMP,&
     & KGPTOT     =RADGRID%NGPTOT,&
     & KGPTOTG    =RADGRID%NGPTOTG,&
     & KGPTOTMX   =RADGRID%NGPTOTMX,&
     & KPTRFRSTLAT=RADGRID%NPTRFRSTLAT,&
     & KFRSTLAT   =RADGRID%NFRSTLAT,&
     & KLSTLAT    =RADGRID%NLSTLAT,&
     & KFRSTLOFF  =RADGRID%NFRSTLOFF,&
     & KSTA       =RADGRID%NSTA(1:RADGRID%NDGLG+N_REGIONS_NS-1,:),&
     & KONL       =RADGRID%NONL(1:RADGRID%NDGLG+N_REGIONS_NS-1,:),&
     & KPTRFLOFF  =RADGRID%NPTRFLOFF,&
     & PMU        =ZMU(1:RADGRID%NDGLG) )  

    RADGRID%RMU(1:RADGRID%NDGLG) = ZMU(1:RADGRID%NDGLG)

    IF( NRADINT == 2 .OR. NRADINT == 3 )THEN
      DO JGL=1,RADGRID%NDGLG
        RADGRID%RSQM2(JGL) = SQRT(1.0_JPRD - ZMU(JGL)*ZMU(JGL))
        RADGRID%RLATIG(JGL) = ASIN(ZMU(JGL))
!       WRITE(NULOUT,'("SUECRAD: JGL=",I6," RADGRID%RLATIG=",F10.3)')&
!        & JGL,RADGRID%RLATIG(JGL)
      ENDDO
      IF(RADGRID%NDGSUR >= 1)THEN
        DO JGLSUR=1,RADGRID%NDGSUR
          ZMU(1-JGLSUR)=ZMU(JGLSUR)
          RADGRID%RMU(1-JGLSUR)=RADGRID%RMU(JGLSUR)
          RADGRID%RSQM2(1-JGLSUR)=RADGRID%RSQM2(JGLSUR)
          RADGRID%RLATIG(1-JGLSUR)=RPI-RADGRID%RLATIG(JGLSUR)
        ENDDO
        DO JGLSUR=1,RADGRID%NDGSUR
          ZMU(RADGRID%NDGLG+JGLSUR)=ZMU(RADGRID%NDGLG+1-JGLSUR)
          RADGRID%RMU(RADGRID%NDGLG+JGLSUR)=RADGRID%RMU(RADGRID%NDGLG+1-JGLSUR)
          RADGRID%RSQM2(RADGRID%NDGLG+JGLSUR)=RADGRID%RSQM2(RADGRID%NDGLG+1-JGLSUR)
          RADGRID%RLATIG(RADGRID%NDGLG+JGLSUR)=-RPI-RADGRID%RLATIG(RADGRID%NDGLG+1-JGLSUR)
        ENDDO
      ENDIF
    ENDIF

    RADGRID%NDGSAL=1
    RADGRID%NDGENL=RADGRID%NLSTLAT(MY_REGION_NS)-RADGRID%NFRSTLOFF
    RADGRID%NDSUR1=NSTENCILWIDE*2-1-MOD(RADGRID%NDLON,2)
    IDLSUR=MAX(RADGRID%NDLON,2*RADGRID%NSMAX+1)
    RADGRID%NDLSUR=IDLSUR+RADGRID%NDSUR1
    RADGRID%MYFRSTACTLAT=RADGRID%NFRSTLAT(MY_REGION_NS)
    RADGRID%MYLSTACTLAT=RADGRID%NLSTLAT(MY_REGION_NS)

    WRITE(NULOUT,'("SUECRAD: RADGRID%NRESOL_ID    =",I8)')RADGRID%NRESOL_ID
    WRITE(NULOUT,'("SUECRAD: RADGRID%NSMAX        =",I8)')RADGRID%NSMAX
    IF( NRADINT == 1 )THEN
      WRITE(NULOUT,'("SUECRAD: RADGRID%NSPEC2       =",I8)')RADGRID%NSPEC2
    ENDIF
    WRITE(NULOUT,'("SUECRAD: RADGRID%NGPTOT       =",I8)')RADGRID%NGPTOT
    WRITE(NULOUT,'("SUECRAD: RADGRID%NGPTOTG      =",I8)')RADGRID%NGPTOTG
    WRITE(NULOUT,'("SUECRAD: RADGRID%NDGSAL       =",I8)')RADGRID%NDGSAL
    WRITE(NULOUT,'("SUECRAD: RADGRID%NDGENL       =",I8)')RADGRID%NDGENL
    WRITE(NULOUT,'("SUECRAD: RADGRID%NDSUR1       =",I8)')RADGRID%NDSUR1
    WRITE(NULOUT,'("SUECRAD: RADGRID%NDLSUR       =",I8)')RADGRID%NDLSUR
    WRITE(NULOUT,'("SUECRAD: RADGRID%MYFRSTACTLAT =",I8)')RADGRID%MYFRSTACTLAT
    WRITE(NULOUT,'("SUECRAD: RADGRID%MYLSTACTLAT  =",I8)')RADGRID%MYLSTACTLAT
    CALL FLUSH(NULOUT)

    IF( NRADINT == 1 )THEN
      ALLOCATE(RADGRID%NASM0(0:RADGRID%NSMAX))
      ALLOCATE(RADGRID%MYMS(RADGRID%NUMP))
      CALL TRANS_INQ (KRESOL     =RADGRID%NRESOL_ID,&
       & KASM0      =RADGRID%NASM0,&
       & KMYMS      =RADGRID%MYMS )  
    ENDIF

    ALLOCATE(RADGRID%GELAM(RADGRID%NGPTOT))
    ALLOCATE(RADGRID%GELAT(RADGRID%NGPTOT))
    ALLOCATE(RADGRID%GESLO(RADGRID%NGPTOT))
    ALLOCATE(RADGRID%GECLO(RADGRID%NGPTOT))
    ALLOCATE(RADGRID%GEMU (RADGRID%NGPTOT))

    IOFF=0
    ILAT=RADGRID%NPTRFLOFF
    DO JGLAT=RADGRID%NFRSTLAT(MY_REGION_NS),&
       & RADGRID%NLSTLAT(MY_REGION_NS)  
      ZGEMU=ZMU(JGLAT)
      ILAT=ILAT+1
      ISTLON  = RADGRID%NSTA(ILAT,MY_REGION_EW)
      IENDLON = ISTLON-1 + RADGRID%NONL(ILAT,MY_REGION_EW)

      DO JLON=ISTLON,IENDLON
        ZLON=  REAL(JLON-1,JPRD)*2.0_JPRD*RPI&
         & /REAL(RADGRID%NLOENG(JGLAT),JPRD)  
        IOFF=IOFF+1
        RADGRID%GELAM(IOFF) = ZLON
        RADGRID%GELAT(IOFF) = ASIN(ZGEMU)
        RADGRID%GESLO(IOFF) = SIN(ZLON)
        RADGRID%GECLO(IOFF) = COS(ZLON)
        RADGRID%GEMU (IOFF) = ZGEMU
      ENDDO
    ENDDO

    DEALLOCATE(ZMU)

    IF( NRADINT == 2 .OR. NRADINT == 3 )THEN

!   For grid point interpolations we need to calculate the halo size
!   required by each processor

      ALLOCATE(ZLATX(RADGRID%NGPTOTMX))
      ALLOCATE(ZLONX(RADGRID%NGPTOTMX))
      DO J=1,RADGRID%NGPTOT
        ZLATX(J)=RADGRID%GELAT(J)/RPI*2.0_JPRD*90.0_JPRD
        ZLONX(J)=(RADGRID%GELAM(J)-RPI)/RPI*180.0_JPRD+180.0_JPRD
      ENDDO
      ZMINRADLAT=MINVAL(ZLATX(1:RADGRID%NGPTOT))
      ZMAXRADLAT=MAXVAL(ZLATX(1:RADGRID%NGPTOT))
      ZMINRADLON=MINVAL(ZLONX(1:RADGRID%NGPTOT))
      ZMAXRADLON=MAXVAL(ZLONX(1:RADGRID%NGPTOT))
      IF( LLDEBUG )THEN
        WRITE(NULOUT,'("RADGRID,BEGIN")')
        IF( MYPROC /= 1 )THEN
          CALL MPL_SEND(RADGRID%NGPTOT,KDEST=NPRCIDS(1),KTAG=1,CDSTRING='SUECRAD.R')
          CALL MPL_SEND(ZLATX(1:RADGRID%NGPTOT),KDEST=NPRCIDS(1),KTAG=2,CDSTRING='SUECRAD.R')
          CALL MPL_SEND(ZLONX(1:RADGRID%NGPTOT),KDEST=NPRCIDS(1),KTAG=3,CDSTRING='SUECRAD.R')
        ENDIF
        IF( MYPROC == 1 )THEN
          DO JROC=1,NPROC
            IF( JROC == MYPROC )THEN
              DO J=1,RADGRID%NGPTOT
                WRITE(NULOUT,'(F7.2,2X,F7.2,2X,I6)')ZLATX(J),ZLONX(J),MYPROC
              ENDDO
            ELSE
              CALL MPL_RECV(IGPTOT,KSOURCE=NPRCIDS(JROC),KTAG=1,CDSTRING='SUECRAD.M')
              CALL MPL_RECV(ZLATX(1:IGPTOT),KSOURCE=NPRCIDS(JROC),KTAG=2,CDSTRING='SUECRAD.M')
              CALL MPL_RECV(ZLONX(1:IGPTOT),KSOURCE=NPRCIDS(JROC),KTAG=3,CDSTRING='SUECRAD.M')
              DO J=1,IGPTOT
                WRITE(NULOUT,'(F7.2,2X,F7.2,2X,I6)')ZLATX(J),ZLONX(J),JROC
              ENDDO
            ENDIF
          ENDDO
        ENDIF
        WRITE(NULOUT,'("RADGRID,END")')
      ENDIF
      DEALLOCATE(ZLATX)
      DEALLOCATE(ZLONX)
  
      ALLOCATE(ZLATX(NGPTOTMX))
      ALLOCATE(ZLONX(NGPTOTMX))
      DO J=1,NGPTOT
        ZLATX(J)=YDGSGEOM_NB%GELAT(J)/RPI*2.0_JPRD*90.0_JPRD
        ZLONX(J)=(YDGSGEOM_NB%GELAM(J)-RPI)/RPI*180.0_JPRD+180.0_JPRD
      ENDDO
      ZMINMDLLAT=MINVAL(ZLATX(1:NGPTOT))
      ZMAXMDLLAT=MAXVAL(ZLATX(1:NGPTOT))
      ZMINMDLLON=MINVAL(ZLONX(1:NGPTOT))
      ZMAXMDLLON=MAXVAL(ZLONX(1:NGPTOT))
      IF( LLDEBUG )THEN
        WRITE(NULOUT,'("MODELGRID,BEGIN")')
        IF( MYPROC /= 1 )THEN
          CALL MPL_SEND(NGPTOT,KDEST=NPRCIDS(1),KTAG=1,CDSTRING='SUECRAD')
          CALL MPL_SEND(ZLATX(1:NGPTOT),KDEST=NPRCIDS(1),KTAG=2,CDSTRING='SUECRAD')
          CALL MPL_SEND(ZLONX(1:NGPTOT),KDEST=NPRCIDS(1),KTAG=3,CDSTRING='SUECRAD')
          CALL MPL_SEND(NGLOBALINDEX(1:NGPTOT),KDEST=NPRCIDS(1),KTAG=4,CDSTRING='SUECRAD')
        ENDIF
        IF( MYPROC == 1 )THEN
          DO JROC=1,NPROC
            IF( JROC == MYPROC )THEN
              DO J=1,NGPTOT
                WRITE(NULOUT,'(F7.2,2X,F7.2,2X,I6,2X,I12)')ZLATX(J),ZLONX(J),MYPROC,NGLOBALINDEX(J)
              ENDDO
            ELSE
              CALL MPL_RECV(IGPTOT,KSOURCE=NPRCIDS(JROC),KTAG=1,CDSTRING='SUECRAD')
              CALL MPL_RECV(ZLATX(1:IGPTOT),KSOURCE=NPRCIDS(JROC),KTAG=2,CDSTRING='SUECRAD')
              CALL MPL_RECV(ZLONX(1:IGPTOT),KSOURCE=NPRCIDS(JROC),KTAG=3,CDSTRING='SUECRAD')
              ALLOCATE(IGLOBALINDEX(1:IGPTOT))
              CALL MPL_RECV(IGLOBALINDEX(1:IGPTOT),KSOURCE=NPRCIDS(JROC),KTAG=4,CDSTRING='SUECRAD')
              DO J=1,IGPTOT
                WRITE(NULOUT,'(F7.2,2X,F7.2,2X,I6,2X,I12)')ZLATX(J),ZLONX(J),JROC,IGLOBALINDEX(J)
              ENDDO
              DEALLOCATE(IGLOBALINDEX)
            ENDIF
          ENDDO
        ENDIF
        WRITE(NULOUT,'("MODELGRID,END")')
      ENDIF
      DEALLOCATE(ZLATX)
      DEALLOCATE(ZLONX)
  
      IF( LLDEBUG )THEN
        WRITE(NULOUT,'("ZMINRADLAT=",F10.2)')ZMINRADLAT
        WRITE(NULOUT,'("ZMINMDLLAT=",F10.2)')ZMINMDLLAT
        WRITE(NULOUT,'("ZMAXRADLAT=",F10.2)')ZMAXRADLAT
        WRITE(NULOUT,'("ZMAXMDLLAT=",F10.2)')ZMAXMDLLAT
        WRITE(NULOUT,'("ZMINRADLON=",F10.2)')ZMINRADLON
        WRITE(NULOUT,'("ZMINMDLLON=",F10.2)')ZMINMDLLON
        WRITE(NULOUT,'("ZMAXRADLON=",F10.2)')ZMAXRADLON
        WRITE(NULOUT,'("ZMAXMDLLON=",F10.2)')ZMAXMDLLON
      ENDIF

      ZLAT=REAL(NDGLG,JPRD)/180.0_JPRD
      ILATS_DIFF_C=CEILING(ABS(ZMINRADLAT-ZMINMDLLAT)*ZLAT)
      ILATS_DIFF_F=FLOOR  (ABS(ZMINRADLAT-ZMINMDLLAT)*ZLAT)
      IF( ZMINRADLAT < ZMINMDLLAT )THEN
        YDRI%NSLWIDES=I_MIN_HALO+ILATS_DIFF_C
      ELSE
        YDRI%NSLWIDES=MAX(0,I_MIN_HALO-ILATS_DIFF_F)
      ENDIF
      ILATS_DIFF_C=CEILING(ABS(ZMAXRADLAT-ZMAXMDLLAT)*ZLAT)
      ILATS_DIFF_F=FLOOR  (ABS(ZMAXRADLAT-ZMAXMDLLAT)*ZLAT)
      IF( ZMAXRADLAT < ZMAXMDLLAT )THEN
        YDRI%NSLWIDEN=MAX(0,I_MIN_HALO-ILATS_DIFF_F)
      ELSE
        YDRI%NSLWIDEN=I_MIN_HALO+ILATS_DIFF_C
      ENDIF
      ZLON=REAL(MAX(NLOENG(NFRSTLAT(MY_REGION_NS)),&
                  & NLOENG(NLSTLAT(MY_REGION_NS))),JPRD)/360.0_JPRD
      ILONS_DIFF_C=CEILING(ABS(ZMINRADLON-ZMINMDLLON)*ZLON)
      ILONS_DIFF_F=FLOOR  (ABS(ZMINRADLON-ZMINMDLLON)*ZLON)
      IF( ZMINRADLON < ZMINMDLLON )THEN
        YDRI%NSLWIDEW=I_MIN_HALO+ILONS_DIFF_C
      ELSE
        YDRI%NSLWIDEW=MAX(0,I_MIN_HALO-ILONS_DIFF_F)
      ENDIF
      ILONS_DIFF_C=CEILING(ABS(ZMAXRADLON-ZMAXMDLLON)*ZLON)
      ILONS_DIFF_F=FLOOR  (ABS(ZMAXRADLON-ZMAXMDLLON)*ZLON)
      IF( ZMAXRADLON < ZMAXMDLLON )THEN
        YDRI%NSLWIDEE=MAX(0,I_MIN_HALO-ILONS_DIFF_F)
      ELSE
        YDRI%NSLWIDEE=I_MIN_HALO+ILONS_DIFF_C
      ENDIF

      ZLAT=REAL(RADGRID%NDGLG,JPRD)/180.0_JPRD
      ILATS_DIFF_C=CEILING(ABS(ZMINRADLAT-ZMINMDLLAT)*ZLAT)
      ILATS_DIFF_F=FLOOR  (ABS(ZMINRADLAT-ZMINMDLLAT)*ZLAT)
      IF( ZMINMDLLAT < ZMINRADLAT )THEN
        YDRO%NSLWIDES=I_MIN_HALO+ILATS_DIFF_C
      ELSE
        YDRO%NSLWIDES=MAX(0,I_MIN_HALO-ILATS_DIFF_F)
      ENDIF
      ILATS_DIFF_C=CEILING(ABS(ZMAXRADLAT-ZMAXMDLLAT)*ZLAT)
      ILATS_DIFF_F=FLOOR  (ABS(ZMAXRADLAT-ZMAXMDLLAT)*ZLAT)
      IF( ZMAXMDLLAT < ZMAXRADLAT )THEN
        YDRO%NSLWIDEN=MAX(0,I_MIN_HALO-ILATS_DIFF_F)
      ELSE
        YDRO%NSLWIDEN=I_MIN_HALO+ILATS_DIFF_C
      ENDIF
      ZLON=REAL(MAX(RADGRID%NLOENG(RADGRID%NFRSTLAT(MY_REGION_NS)),&
                  & RADGRID%NLOENG(RADGRID%NLSTLAT(MY_REGION_NS))),JPRD)/360.0_JPRD
      ILONS_DIFF_C=CEILING(ABS(ZMINRADLON-ZMINMDLLON)*ZLON)
      ILONS_DIFF_F=FLOOR  (ABS(ZMINRADLON-ZMINMDLLON)*ZLON)
      IF( ZMINMDLLON < ZMINRADLON )THEN
        YDRO%NSLWIDEW=I_MIN_HALO+ILONS_DIFF_C
      ELSE
        YDRO%NSLWIDEW=MAX(0,I_MIN_HALO-ILONS_DIFF_F)
      ENDIF
      ILONS_DIFF_C=CEILING(ABS(ZMAXRADLON-ZMAXMDLLON)*ZLON)
      ILONS_DIFF_F=FLOOR  (ABS(ZMAXRADLON-ZMAXMDLLON)*ZLON)
      IF( ZMAXMDLLON < ZMAXRADLON )THEN
        YDRO%NSLWIDEE=MAX(0,I_MIN_HALO-ILONS_DIFF_F)
      ELSE
        YDRO%NSLWIDEE=I_MIN_HALO+ILONS_DIFF_C
      ENDIF

    ENDIF

    RADGRID%NDGSAH=MAX(RADGRID%NDGSAG,&
     & RADGRID%NDGSAL+RADGRID%NFRSTLOFF-YDRO%NSLWIDEN)-RADGRID%NFRSTLOFF  
    RADGRID%NDGENH=MIN(RADGRID%NDGENG,&
     & RADGRID%NDGENL+RADGRID%NFRSTLOFF+YDRO%NSLWIDES)-RADGRID%NFRSTLOFF  
    WRITE(NULOUT,'("SUECRAD: RADGRID%NDGSAH       =",I8)')RADGRID%NDGSAH
    WRITE(NULOUT,'("SUECRAD: RADGRID%NDGENH       =",I8)')RADGRID%NDGENH

    IF( NRADINT == 2 .OR. NRADINT == 3 )THEN

      ILBRLATI = MAX(RADGRID%NDGSAG,&
       & RADGRID%NDGSAL+RADGRID%NFRSTLOFF-YDRO%NSLWIDEN)-RADGRID%NFRSTLOFF  
      IUBRLATI = MIN(RADGRID%NDGENG,&
       & RADGRID%NDGENL+RADGRID%NFRSTLOFF+YDRO%NSLWIDES)-RADGRID%NFRSTLOFF  
      ALLOCATE(RADGRID%RLATI(ILBRLATI:IUBRLATI))
      ALLOCATE(RADGRID%RIPI(ILBRLATI:IUBRLATI,3))
  
      DO JGL= ILBRLATI,IUBRLATI
        IGLGLO=JGL+RADGRID%NFRSTLOFF
        IF(IGLGLO >= 0.AND.IGLGLO <= RADGRID%NDGLG) THEN
          ZD1=RADGRID%RLATIG(IGLGLO-1)-RADGRID%RLATIG(IGLGLO)
          ZD2=RADGRID%RLATIG(IGLGLO-1)-RADGRID%RLATIG(IGLGLO+1)
          ZD3=RADGRID%RLATIG(IGLGLO-1)-RADGRID%RLATIG(IGLGLO+2)
          ZD4=RADGRID%RLATIG(IGLGLO  )-RADGRID%RLATIG(IGLGLO+1)
          ZD5=RADGRID%RLATIG(IGLGLO  )-RADGRID%RLATIG(IGLGLO+2)
          ZD6=RADGRID%RLATIG(IGLGLO+1)-RADGRID%RLATIG(IGLGLO+2)
          RADGRID%RIPI(JGL,1)=-1.0_JPRD/(ZD1*ZD4*ZD5)
          RADGRID%RIPI(JGL,2)= 1.0_JPRD/(ZD2*ZD4*ZD6)
          RADGRID%RIPI(JGL,3)=-1.0_JPRD/(ZD3*ZD5*ZD6)
        ENDIF
        RADGRID%RLATI(JGL)=RADGRID%RLATIG(IGLGLO)
      ENDDO

      YDRI%CVER='RI'
      IDGSAH = MAX(NDGSAG,NDGSAL+NFRSTLOFF-YDRI%NSLWIDEN)-NFRSTLOFF
      IDGENH = MIN(NDGENG,NDGENL+NFRSTLOFF+YDRI%NSLWIDES)-NFRSTLOFF
      CALL SLCSET(YDGEOMETRY,YDRI,MY_REGION_NS,MY_REGION_EW,&
       & NDGLG,NDLON,NDGSAG,NDGENG,NDGSAL,NDGENL,IDGSAH,IDGENH,&
       & IDUM,IDUM,IDUM,IDUM,IDUM,&
       & NDSUR1,NDLSUR,NDGSUR,NGPTOT,&
       & NPTRFLOFF,NFRSTLOFF,MYFRSTACTLAT,MYLSTACTLAT,&
       & NSTA,NONL,NLOENG,NPTRFRSTLAT,NFRSTLAT,NLSTLAT,&
       & YDCSGLEG%RMU,YDCSGLEG%RSQM2,YDSLREP=YDMODEL%YRML_DYN%YRSLREP)
      WRITE(NULOUT,'("SUECRAD: YDRI%NASLB1=",I12)')YDRI%NASLB1

      YDRO%CVER='RO'
      CALL SLCSET(YDGEOMETRY,YDRO,MY_REGION_NS,MY_REGION_EW,&
       & RADGRID%NDGLG,RADGRID%NDLON,RADGRID%NDGSAG,&
       & RADGRID%NDGENG,RADGRID%NDGSAL,RADGRID%NDGENL,RADGRID%NDGSAH,RADGRID%NDGENH,&
       & IDUM,IDUM,IDUM,IDUM,IDUM,&
       & RADGRID%NDSUR1,RADGRID%NDLSUR,RADGRID%NDGSUR,RADGRID%NGPTOT,&
       & RADGRID%NPTRFLOFF,RADGRID%NFRSTLOFF,RADGRID%MYFRSTACTLAT,RADGRID%MYLSTACTLAT,&
       & RADGRID%NSTA,RADGRID%NONL,RADGRID%NLOENG,RADGRID%NPTRFRSTLAT,&
       & RADGRID%NFRSTLAT,RADGRID%NLSTLAT,&
       & RADGRID%RMU,RADGRID%RSQM2,YDSLREP=YDMODEL%YRML_DYN%YRSLREP)
      WRITE(NULOUT,'("SUECRAD: YDRO%NASLB1=",I12)')YDRO%NASLB1

      IF( LLDEBUG )THEN
        WRITE(NULOUT,'("")')
        IRIWIDEMAXN=0
        IRIWIDEMAXS=0
        IRIWIDEMAXW=0
        IRIWIDEMAXE=0
        IROWIDEMAXN=0
        IROWIDEMAXS=0
        IROWIDEMAXW=0
        IROWIDEMAXE=0
        IARIB1MAX=0
        IAROB1MAX=0
        IWIDE(1)=YDRI%NSLWIDEN
        IWIDE(2)=YDRI%NSLWIDES
        IWIDE(3)=YDRI%NSLWIDEW
        IWIDE(4)=YDRI%NSLWIDEE
        IWIDE(5)=YDRO%NSLWIDEN
        IWIDE(6)=YDRO%NSLWIDES
        IWIDE(7)=YDRO%NSLWIDEW
        IWIDE(8)=YDRO%NSLWIDEE
        IWIDE(9)=YDRI%NASLB1
        IWIDE(10)=YDRO%NASLB1
        IF( MYPROC /= 1 )THEN
          CALL MPL_SEND(IWIDE(1:10),KDEST=NPRCIDS(1),KTAG=1,CDSTRING='SUECRAD.W')
        ENDIF
        IF( MYPROC == 1 )THEN
          DO JROC=1,NPROC
            IF( JROC /= MYPROC )THEN
              CALL MPL_RECV(IWIDE(1:10),KSOURCE=NPRCIDS(JROC),KTAG=1,CDSTRING='SUECRAD.W')
            ENDIF
            WRITE(NULOUT,'("SUECRAD: PROC=",I5,2X,"YDRI%NSLWIDEN=",I3,2X,"YDRO%NSLWIDEN=",I3 )')&
             & JROC,IWIDE(1),IWIDE(5)  
            WRITE(NULOUT,'("SUECRAD: PROC=",I5,2X,"YDRI%NSLWIDES=",I3,2X,"YDRO%NSLWIDES=",I3 )')&
             & JROC,IWIDE(2),IWIDE(6)  
            WRITE(NULOUT,'("SUECRAD: PROC=",I5,2X,"YDRI%NSLWIDEW=",I3,2X,"YDRO%NSLWIDEW=",I3 )')&
             & JROC,IWIDE(3),IWIDE(7)  
            WRITE(NULOUT,'("SUECRAD: PROC=",I5,2X,"YDRI%NSLWIDEE=",I3,2X,"YDRO%NSLWIDEE=",I3 )')&
             & JROC,IWIDE(4),IWIDE(8)  
            WRITE(NULOUT,'("SUECRAD: PROC=",I5,2X,"YDRI%NASLB1=",I10,2X,"YDRO%NASLB1=",I10 )')&
             & JROC,IWIDE(9),IWIDE(10)
            WRITE(NULOUT,'("")')
            IF( IWIDE(1) > IRIWIDEMAXN ) IRIWIDEMAXN=IWIDE(1)
            IF( IWIDE(2) > IRIWIDEMAXS ) IRIWIDEMAXS=IWIDE(2)
            IF( IWIDE(3) > IRIWIDEMAXW ) IRIWIDEMAXW=IWIDE(3)
            IF( IWIDE(4) > IRIWIDEMAXE ) IRIWIDEMAXE=IWIDE(4)
            IF( IWIDE(5) > IROWIDEMAXN ) IROWIDEMAXN=IWIDE(5)
            IF( IWIDE(6) > IROWIDEMAXS ) IROWIDEMAXS=IWIDE(6)
            IF( IWIDE(7) > IROWIDEMAXW ) IROWIDEMAXW=IWIDE(7)
            IF( IWIDE(8) > IROWIDEMAXE ) IROWIDEMAXE=IWIDE(8)
            IF( IWIDE(9)  > IARIB1MAX  ) IARIB1MAX  =IWIDE(9)
            IF( IWIDE(10) > IAROB1MAX  ) IAROB1MAX  =IWIDE(10)
          ENDDO
          WRITE(NULOUT,'("")')
          WRITE(NULOUT,'("SUECRAD: YDRI%NSLWIDEN(MAX)  =",I8)')IRIWIDEMAXN
          WRITE(NULOUT,'("SUECRAD: YDRI%NSLWIDES(MAX)  =",I8)')IRIWIDEMAXS
          WRITE(NULOUT,'("SUECRAD: YDRI%NSLWIDEW(MAX)  =",I8)')IRIWIDEMAXW
          WRITE(NULOUT,'("SUECRAD: YDRI%NSLWIDEE(MAX)  =",I8)')IRIWIDEMAXE
          WRITE(NULOUT,'("SUECRAD: YDRO%NSLWIDEN(MAX)  =",I8)')IROWIDEMAXN
          WRITE(NULOUT,'("SUECRAD: YDRO%NSLWIDES(MAX)  =",I8)')IROWIDEMAXS
          WRITE(NULOUT,'("SUECRAD: YDRO%NSLWIDEW(MAX)  =",I8)')IROWIDEMAXW
          WRITE(NULOUT,'("SUECRAD: YDRO%NSLWIDEE(MAX)  =",I8)')IROWIDEMAXE
          WRITE(NULOUT,'("SUECRAD: YDRI%NASLB1(MAX)    =",I10)')IARIB1MAX
          WRITE(NULOUT,'("SUECRAD: YDRO%NASLB1(MAX)    =",I10)')IAROB1MAX
          WRITE(NULOUT,'("")')
        ENDIF
        CALL FLUSH(NULOUT)
      ENDIF

    ENDIF
    CALL GSTATS(1818,1)

  ELSE

    WRITE(NULOUT,'("SUECRAD: INVALID VALUE FOR NRADINT=",I6)')NRADINT
    CALL ABOR1('SUECRAD: NRADINT INVALID')

  ENDIF

!  Adjust NRPROMA to optimise for memory usage and vector length
!  Only do this if LOPTRPROMA is true

  IF (LOPTRPROMA) THEN
!   Allow nrproma to drift upwards by up to 5 percent if this results
!   in a reduction of the number of grid point blocks
    ITHRMAX=OML_MAX_THREADS()
    IF( ITHRMAX > 1 )THEN
!    Optimise NRPROMA for maximum number of threads
      IRPROMA=NRPROMA*1.05
      NRPROMA=(RADGRID%NGPTOTMX-1)/ITHRMAX+1
      I=1
      DO WHILE( NRPROMA > IRPROMA )
        I=I+1
        NRPROMA=(RADGRID%NGPTOTMX-1)/(ITHRMAX*I)+1
      ENDDO
! Check if we have enough blocks to switch to a more optimal nrproma
! i.e. more than 10 blocks per thread
      IGPBLKSMX=(RADGRID%NGPTOTMX-1)/NRPROMA+1
      IF( IGPBLKSMX/ITHRMAX > 10 )THEN
        NRPROMA=8
      ENDIF
    ELSE
      IRPROMABEG=NRPROMA+1
      IRPROMAEND=IRPROMABEG+NRPROMA*5/100
      IGPBLKSMX=(RADGRID%NGPTOTMX-1)/NRPROMA+1
      NRPROMA=(RADGRID%NGPTOTMX-1)/IGPBLKSMX+1
      DO JJ=IRPROMABEG,IRPROMAEND
        IMX=(RADGRID%NGPTOTMX-1)/JJ+1
        IPR=(RADGRID%NGPTOTMX-1)/IMX+1
        IF (IMX < IGPBLKSMX) THEN
          IGPBLKSMX=IMX
          NRPROMA=IPR
        ENDIF
      ENDDO
    ENDIF
  ENDIF

ENDIF              ! END OF LERADI BLOCK


! Default value of NPROMALW set to NPROMA

NPROMALW = NPROMA

! Set-up for the linearized (TL/AD) longwave radiation

IF (LERADLW2) THEN
  
!  Adjust NPROMALW to optimise for memory usage and vector length
!  Only do this if LOPTLWPR is true

!mj  IF (LOPTLWPR .AND. LH2OCO2) THEN
  IF (LOPTLWPR) THEN

! First estimate based on number of suggested loops inside of LWAD

    NPROMALW = NPROMA/NLOOPLW

! Adjustment for optimal length    

    IRPROMABEG=NPROMALW+1
    IRPROMAEND=IRPROMABEG+NPROMALW*5/100
    IGPBLKSMX=(NPROMA-1)/NPROMALW+1
    NPROMALW=(NPROMA-1)/IGPBLKSMX+1
    DO JJ=IRPROMABEG,IRPROMAEND
      IMX=(NPROMA-1)/JJ+1
      IPR=(NPROMA-1)/IMX+1
      IF (IMX < IGPBLKSMX) THEN
        IGPBLKSMX=IMX
        NPROMALW=IPR
      ENDIF
    ENDDO
  ENDIF

! Length of the last loop inside of LWAD after NPROMALW optimization

  NLASTLW = MOD(NPROMA,NPROMALW)

! Number of optimzed loops 

  NLOOPLWO = NPROMA/NPROMALW
  IF (NLASTLW /= 0) THEN
    NLOOPLWO = NLOOPLWO +1
  ENDIF
   
ELSE
  NLOOPLWO=0
  NLASTLW=0
ENDIF

!      ----------------------------------------------------------------

!*       4.    INITIALIZE RADIATION COEFFICIENTS.
!              ----------------------------------

DIFF   = 1.66_JPRD
R10E   = 0.4342945_JPRD

CALL GSTATS(1818,0)
CALL SURDI(YDERDI)

IF (NINHOM == 0) THEN
  RLWINHF=1._JPRD
  RSWINHF=1._JPRD
ENDIF

!      ----------------------------------------------------------------

!*       5.    INITIALIZE RADIATION ABSORPTION COEFFICIENTS
!              --------------------------------------------

!*       5.1.  Initialization routine for RRTM
!              -------------------------------

CALL SURRTAB
CALL SURRTPK
CALL SURRTRF
!!! CALL SURRTFTR 
CALL SU_CLOP550

IF (LRRTM) THEN
  IF (KLEV > NFLEVG) THEN
    WRITE(UNIT=KULOUT,&
     & FMT='('' RRTM MAXIMUM NUMBER OF LAYERS IS REACHED'',&
     & '' CALL ABORT'')')  
    CALL ABOR1(' ABOR1 CALLED SUECRAD')
  ENDIF
    
! Read the absorption coefficient data and reduce from 256 to 140 g-points

  CALL RRTM_INIT_140GP

  INBLW=16

ELSE
  INBLW=6

ENDIF

CALL SULWN

LACR1=LRAY.AND.(NSW == 1.AND.NSWB_MNH == 1)
LACR6=LRAY.AND.(NSW == 6.AND.NSWB_MNH == 6).AND.LHLRADUPD

IF(.NOT.LRAY) THEN
  CALL SUSWN(YDESWRT,YDERAD,NTSW,NSW)
  CALL SUCLOPN(YDERAD,YDESWRT,NTSW,NSW,KLEV)
ENDIF

!-- routines specific to SRTM
IF (LSRTM) THEN
  NTSW=14
  ISW =14
  CALL SRTM_INIT(NSWWVCONTINUUM)
  CALL SUSRTAB
  CALL SUSRTAER
  CALL SUSRTCOP
  WRITE(UNIT=KULOUT,FMT='(''SRTM Configuration'',L8,3I4)')LSRTM,NTSW,ISW,JPGPT

!-- other IFS radiation
ELSEIF (.NOT.LRAY) THEN
  IF (.NOT.LONEWSW .OR. ((NSW /= 2).AND.(NSW /= 4).AND.(NSW /= 6)) ) THEN
    WRITE(UNIT=KULOUT,FMT='(''Wrong SW Configuration'',L8,I3)')LONEWSW,NSW
  ENDIF
  CALL SUSWN(YDESWRT,YDERAD,NTSW,NSW)
  CALL SUAERSN (YDESWRT,NTSW,NSW)

!-- acraneb
ELSEIF (.NOT.(LACR1.OR.LACR6)) THEN
  WRITE(UNIT=KULOUT,FMT='(A)')&
   & 'For LRAY=.T. number of solar intervals NSW and NSWB_MNH must be 1 or 6!'
  CALL ABOR1('SUECRAD: ABORT')
ELSE
  CALL SUSWN(YDESWRT,YDERAD,NTSW,NSW)

ENDIF
WRITE(UNIT=KULOUT,FMT='('' INBLW,NTSW,NSW SET EQUAL TO:'',3I3)') INBLW,NTSW,NSW


!-- routine specific to the UV processor
IF (LUVPROC) THEN
  NUVTIM = NUVTIM * INT(RDAY)
  CALL SU_UVRADI(YDEUVRAD,NUV)
ENDIF

!      ----------------------------------------------------------------

!*       6.    INITIALIZE AEROSOL OPTICAL PARAMETERS AND DISTRIBUTION
!              ------------------------------------------------------

! -      6.1.  PARAMETERS FOR THE STANDARD TEGEN AEROSOL CLIMATOLOGY
  !              -----------------------------------------------------

!- LW optical properties
CALL SUAERL(YDMODEL%YRML_PHY_RAD%YRESWRT%RAER)
!- SW optical properties moved above
!CALL SUAERSN (NTSW,NSW)

!- horizontal distribution      
CALL SUAERH(YDEAERD)

!- vertical distribution
CALL SUAERV ( YDERAD, KLEV  , PETAH,&
 & CVDAES , CVDAEL , CVDAEU , CVDAED,&
 & RCTRBGA, RCVOBGA, RCSTBGA, RCAEOPS, RCAEOPL, RCAEOPU,&
 & RCAEOPD, RCTRPT , RCAEADK, RCAEADM, RCAEROS&
 & ) 
 
! -      6.2.  PARAMETERS FOR THE MACC-DERIVED AEROSOL CLIMATOLOGY
!              ---------------------------------------------------

IF (NAERMACC == 1) THEN
  LAERADCLI=.TRUE.

  CALL SU_AERV ( YDERAD, KLEV, PETAH )

! - define numbers of latitudes, longitudes, levels and variables in 
!   the CAMS(MACC) aerosol climatology
  NMCLAT=61
  NMCLON=120
  NMCVAR=12

  IF (LAER3D) THEN
    NMCLEV=60
  ELSE
    NMCLEV=1
  ENDIF

  IF (NACTAERO == 0) THEN
    ! This won't have already been called "early" if there's no prognostic aerosol,
    ! but must still be called prior to initialising the radiation scheme, not
    ! fully "late" in SUPHEC.
    CALL SU_AERW(YDMODEL)
  END IF
ENDIF

! -      6.3.  OTHER PARAMETERS
!              ----------------

!-- Overlap function (only used if NOVLP=4)
CALL SUOVLP(YDSTA,YDMODEL%YRML_PHY_RAD%YREOVLP,YDERAD,KLEV)

!-- parameters for prognostic aerosols: now called from phys_ec/suphec.F90
!IF (NACTAERO > 0) THEN
!  CALL SU_AERW
!ENDIF

!      ----------------------------------------------------------------

!*       7.    INITIALIZE SATELLITE GEOMETRICAL/RADIOMETRIC PARAMETERS
!              -------------------------------------------------------

!IF (LEPHYS .AND. NMODE > 1) THEN
!  CALL SUSAT
!ENDIF
CALL GSTATS(1818,1)

!      ----------------------------------------------------------------

!*       8.    INITIALIZE CLIMATOLOGICAL OZONE DISTRIBUTION
!              --------------------------------------------
!                  (not done here!!!  called from APLPAR as it depends
!                     on model pressure levels!)

!      ----------------------------------------------------------------

!*       9.    SET UP MODEL CONFIGURATION FOR TIME-SPACE INTERPOLATION
!              -------------------------------------------------------

ZTSTEP=MAX(TSTEP,1.0_JPRD)
ZSTPHR=RHOUR/ZTSTEP
IRADFR=NRADFR
IF(NRADFR < 0) THEN
  NRADFR=-NRADFR*ZSTPHR+0.5_JPRD
ENDIF
NRADPFR=NRADPFR*NRADFR
IF (MOD(NRADPLA,2) == 0.AND. NRADPLA /= 0) THEN
  NRADPLA=NRADPLA+1
ENDIF

IF(NRADUV < 0) THEN
  NRADUV=-NRADUV*ZSTPHR+0.5_JPRD
ENDIF

IST1HR=ZSTPHR+0.05_JPRD
ISTNHR=  NLNGR1H *ZSTPHR+0.05_JPRD
IF (MOD(RHOUR,ZTSTEP) > 0.1_JPRD) THEN
  801 CONTINUE
  IST1HR=IST1HR+1
  IF (MOD(ISTNHR,IST1HR) /= 0) GO TO 801
ENDIF
IF (NRADFR == 1) THEN
  NRADSFR=NRADFR
ELSE
  NRADSFR=IST1HR
ENDIF
NRADNFR=NRADFR
NRADE1H=IST1HR
NRADE3H=3*IST1HR

IF(LRAYFM) THEN
  NRPROMA=NDLON+6+(1-MOD(NDLON,2))
ENDIF

!      ----------------------------------------------------------------

!*       10.    ALLOCATE WORK ARRAYS
!               --------------------

IU = NULOUT
LLP = NPRINTLEV >= 1.OR. LALLOPR

IF (LEPHYS) THEN
  ALLOCATE(YDRADF%EMTD(NPROMA,NFLEVG+1,NGPBLKS))
  IF(LLP)WRITE(IU,9) 'EMTD     ',SIZE(YDRADF%EMTD     ),SHAPE(YDRADF%EMTD     )
  ALLOCATE(YDRADF%TRSW(NPROMA,NFLEVG+1,NGPBLKS))
  IF(LLP)WRITE(IU,9) 'TRSW     ',SIZE(YDRADF%TRSW     ),SHAPE(YDRADF%TRSW     )
  ALLOCATE(YDRADF%EMTC(NPROMA,NFLEVG+1,NGPBLKS))
  IF(LLP)WRITE(IU,9) 'EMTC     ',SIZE(YDRADF%EMTC     ),SHAPE(YDRADF%EMTC     )
  ALLOCATE(YDRADF%TRSC(NPROMA,NFLEVG+1,NGPBLKS))
  IF(LLP)WRITE(IU,9) 'TRSC     ',SIZE(YDRADF%TRSC     ),SHAPE(YDRADF%TRSC     )

!  IF (LEPAERO) THEN
    ALLOCATE(YDRADF%TAUAER(NPROMA,NFLEVG,6,NGPBLKS))
    IF(LLP)WRITE(IU,9) 'TAUAER   ',SIZE(YDRADF%TAUAER   ),SHAPE(YDRADF%TAUAER   )
    YDRADF%TAUAER = 0._JPRD
!  ENDIF

  ALLOCATE(YDRADF%SRSWD(NPROMA,NGPBLKS))
  IF(LLP)WRITE(IU,9) 'SRSWD    ',SIZE(YDRADF%SRSWD    ),SHAPE(YDRADF%SRSWD    )
  ALLOCATE(YDRADF%SRSWDC(NPROMA,NGPBLKS))                               
  IF(LLP)WRITE(IU,9) 'SRSWDC   ',SIZE(YDRADF%SRSWDC   ),SHAPE(YDRADF%SRSWDC   )
  ALLOCATE(YDRADF%SRLWD(NPROMA,YDERAD%NLWOUT,NGPBLKS))
  IF(LLP)WRITE(IU,9) 'SRLWD    ',SIZE(YDRADF%SRLWD    ),SHAPE(YDRADF%SRLWD    )
  ALLOCATE(YDRADF%SRLWDC(NPROMA,NGPBLKS))                               
  IF(LLP)WRITE(IU,9) 'SRLWDC   ',SIZE(YDRADF%SRLWDC   ),SHAPE(YDRADF%SRLWDC   )
  ALLOCATE(YDRADF%SRSWDCS(NPROMA,NGPBLKS))
  IF(LLP)WRITE(IU,9) 'SRSWDCS  ',SIZE(YDRADF%SRSWDCS  ),SHAPE(YDRADF%SRSWDCS  )
  ALLOCATE(YDRADF%SRLWDCS(NPROMA,NGPBLKS))
  IF(LLP)WRITE(IU,9) 'SRLWDCS  ',SIZE(YDRADF%SRLWDCS  ),SHAPE(YDRADF%SRLWDCS  )
  ALLOCATE(YDRADF%SRSWDV(NPROMA,NGPBLKS))
  IF(LLP)WRITE(IU,9) 'SRSWDV   ',SIZE(YDRADF%SRSWDV   ),SHAPE(YDRADF%SRSWDV   )
  ALLOCATE(YDRADF%SRSWDUV(NPROMA,NGPBLKS))
  IF(LLP)WRITE(IU,9) 'SRSWDUV  ',SIZE(YDRADF%SRSWDUV  ),SHAPE(YDRADF%SRSWDUV  )
  ALLOCATE(YDRADF%EDRO(NPROMA,NGPBLKS))
  IF(LLP)WRITE(IU,9) 'EDRO     ',SIZE(YDRADF%EDRO     ),SHAPE(YDRADF%EDRO     )
  ALLOCATE(YDRADF%SRSWPAR(NPROMA,NGPBLKS))
  IF(LLP)WRITE(IU,9) 'SRSWPAR  ',SIZE(YDRADF%SRSWPAR  ),SHAPE(YDRADF%SRSWPAR  )
  ALLOCATE(YDRADF%SRSWUVB(NPROMA,NGPBLKS))
  IF(LLP)WRITE(IU,9) 'SRSWUVB  ',SIZE(YDRADF%SRSWUVB  ),SHAPE(YDRADF%SRSWUVB  )

ELSEIF(LMPHYS .AND. (LRAYFM.OR.LRAYFM15)) THEN
  ALLOCATE(YDRADF%EMTD(NPROMA,NFLEVG+1,NGPBLKS))
  IF(LLP)WRITE(IU,9) 'EMTD     ',SIZE(YDRADF%EMTD     ),SHAPE(YDRADF%EMTD     )
  ALLOCATE(YDRADF%TRSW(NPROMA,NFLEVG+1,NGPBLKS))
  IF(LLP)WRITE(IU,9) 'TRSW     ',SIZE(YDRADF%TRSW     ),SHAPE(YDRADF%TRSW     )
  ALLOCATE(YDRADF%EMTU(NPROMA,NFLEVG+1,NGPBLKS))
  IF(LLP)WRITE(IU,9) 'EMTU     ',SIZE(YDRADF%EMTU     ),SHAPE(YDRADF%EMTU     )
  ALLOCATE(YDRADF%RMOON(NPROMA,NGPBLKS))
  IF(LLP)WRITE(IU,9) 'RMOON    ',SIZE(YDRADF%RMOON    ),SHAPE(YDRADF%RMOON    )
ENDIF
ALLOCATE(YDRADF%SRSWPARC(NPROMA,NGPBLKS))
IF(LLP)WRITE(IU,9) 'SRSWPARC ',SIZE(YDRADF%SRSWPARC ),SHAPE(YDRADF%SRSWPARC )
ALLOCATE(YDRADF%SRSWTINC(NPROMA,NGPBLKS))
IF(LLP)WRITE(IU,9) 'SRSWTINC ',SIZE(YDRADF%SRSWTINC ),SHAPE(YDRADF%SRSWTINC )
ALLOCATE(YDRADF%SRFDIR(NPROMA,NGPBLKS))
IF(LLP)WRITE(IU,9) 'SRFDIR   ',SIZE(YDRADF%SRFDIR   ),SHAPE(YDRADF%SRFDIR   )
ALLOCATE(YDRADF%SRCDIR(NPROMA,NGPBLKS))
IF(LLP)WRITE(IU,9) 'SRCDIR   ',SIZE(YDRADF%SRCDIR   ),SHAPE(YDRADF%SRCDIR   )

IF (LAPPROXLWUPDATE) THEN
  ALLOCATE(YDRADF%DERIVATIVELW(NPROMA,NFLEVG+1,NGPBLKS))
  IF(LLP)WRITE(IU,9) 'DerivativeLw ',SIZE(YDRADF%DERIVATIVELW),SHAPE(YDRADF%DERIVATIVELW)
ELSE
  ! Allocation of a minimal sized array to avoid problems when it is
  ! passed to subroutines, even if it is not used
  ALLOCATE(YDRADF%DERIVATIVELW(1,1,1))
ENDIF

9 FORMAT(1X,'ARRAY ',A10,' ALLOCATED ',8I8)

!      ----------------------------------------------------------------

!*       10.    PRINT FINAL VALUES.
!               -------------------

IF (LOUTPUT) THEN
  WRITE(UNIT=KULOUT,FMT='('' COMMON YOERAD '')')
  WRITE(UNIT=KULOUT,FMT='('' LERADI  = '',L5&
   & ,'' LERAD1H = '',L5,'' LHGHG = '',L5&
   & ,'' NLNGR1H = '',I2,'' NRADSFR = '',I2)')&
   & LERADI,LERAD1H,LHGHG,NLNGR1H,NRADSFR  
  WRITE(UNIT=KULOUT,FMT='('' LEPO3RA  = '',L5,'' YO3%LGP = '',L5 )') LEPO3RA,YO3%LGP
  WRITE(UNIT=KULOUT,FMT='('' NRADFR  = '',I2&
   & ,'' NRADPFR = '',I3&
   & ,'' NRADE1H = '',I3&
   & ,'' NRADE3H = '',I3&
   & ,'' NRADELG = '',I3&
   & ,'' NRADPLA = '',I2&
   & ,'' NRPROMA = '',I5&
   & )')&
   & NRADFR,NRADPFR,NRADE1H,NRADE3H,NRADELG,NRADPLA,NRPROMA
  WRITE(UNIT=KULOUT,FMT='('' NLOOPLWO= '',I5&
   & ,'' NPROMALW= '',I5&
   & ,'' NLASTLW = '',I5&
   & )')&
   & NLOOPLWO,NPROMALW,NLASTLW
  WRITE(UNIT=KULOUT,FMT='('' LRRTM = '',L5&
   & ,'' LSRTM = '',L5&
   & ,'' NMODE = '',I1&
   & ,'' NOZOCL= '',I1&
   & ,'' NAER  = '',I1&
   & ,'' NAERMACC='',I1&
   & ,'' LAER3D=  '',L5&
   & ,'' NHINCSOL='',I2&
   & ,'' NSOLARSPECTRUM='',I2&
   & ,'' NSWWVCONTINUUM='',I2&
   & ,'' NDUMPBADINPUTS='',I2&
   & ,'' NDUMPINPUTS='',I2&
   & ,'' NLWEMISS='',I2&
   & ,'' NLWOUT='',I2&
   & ,'' LINTERPINCLOUDMEAN='',L5&
   & ,'' TRBKG= '',F8.3&
   & ,'' STBKG= '',F8.3&
   & )')&
   & LRRTM,LSRTM,NMODE,NOZOCL,NAER,NAERMACC, LAER3D,NHINCSOL,NSOLARSPECTRUM,&
   & NSWWVCONTINUUM, NDUMPBADINPUTS, NDUMPINPUTS, &
   & YDERAD%NLWEMISS, YDERAD%NLWOUT, LINTERPINCLOUDMEAN, TRBKG,STBKG

  WRITE(UNIT=KULOUT,FMT='('' CGHGCLIMFILE =         '',A)') TRIM(CGHGCLIMFILE)
  WRITE(UNIT=KULOUT,FMT='('' CGHGTIMESERIESFILE =   '',A)') TRIM(CGHGTIMESERIESFILE)
  WRITE(UNIT=KULOUT,FMT='('' CSOLARIRRADIANCEFILE = '',A)') TRIM(CSOLARIRRADIANCEFILE)

  IF (.NOT.LHGHG) WRITE(UNIT=KULOUT,FMT='('' RCCO2= '',E10.3&
    &,'' RCCH4= '',E10.3,'' RCN2O= '',E10.3,'' RCCFC11= '',E10.3,'' RCFC12= '',E10.3&
    &,'' RCCFC22= '',E10.3,'' RCCCL4= '',E10.3)')&
    & RCCO2,RCCH4,RCN2O,RCCFC11,RCCFC12,RCCFC22,RCCCL4
  IF (NAERMACC /= 0) THEN
  WRITE(UNIT=KULOUT,FMT='('' NMCLAT = '',I4&
   & ,'' NMCLON = '',I4&
   & ,'' NMCLEV = '',I4&
   & ,'' NMCVAR = '',I4&
   & ,'' LAERADCLI= '',L1&
   & )')&
   & NMCLAT, NMCLON, NMCLEV, NMCVAR, LAERADCLI
  WRITE(UNIT=KULOUT,FMT='('' LAERADJDU = '',L1&
   & ,'' NAERADJDU = '',I1&
   & ,'' RDUMULF  = '',F5.1&
   & ,'' RWGHTDU1 = '',F5.1&
   & ,'' RWGHTDU2 = '',F5.1&
   & ,'' RWGHTDU3 = '',F5.1&
   & )')&
   & LAERADJDU, NAERADJDU, RDUMULF, RWGHTDU1, RWGHTDU2, RWGHTDU3
  WRITE(UNIT=KULOUT,FMT='('' LDUSEASON = '',L5&
   & ,'' RAESHSS = '',F7.1&
   & ,'' RAESHDU = '',F7.1&
   & ,'' RAESHOM = '',F7.1&
   & ,'' RAESHBC = '',F7.1&
   & ,'' RAESHSU = '',F7.1&
   & )')&
   & LDUSEASON,RAESHSS, RAESHDU, RAESHOM, RAESHBC, RAESHSU
  ENDIF
  WRITE(UNIT=KULOUT,FMT='('' NINHOM = '',I1&
   & ,'' NLAYINH='',I1&
   & ,'' RLWINHF='',F4.2&
   & ,'' RSWINHF='',F4.2&
   & )')&
   & NINHOM,NLAYINH,RLWINHF,RSWINHF  
  IF (NPERTAER /= 0 .OR. NPERTOZ /= 0) THEN
    WRITE(UNIT=KULOUT,FMT='('' NPERTAER= '',I2&
   & ,'' LNOTROAER='',L5&
   & ,'' NPERTOZ = '',I1&
   & ,'' RPERTOZ = '',F5.0&
   & )')&
   & NPERTAER,LNOTROAER,NPERTOZ,RPERTOZ
  ENDIF
  WRITE(UNIT=KULOUT,FMT='('' NRADINT = '',I2)')NRADINT
  WRITE(UNIT=KULOUT,FMT='('' NRADRES = '',I4)')NRADRES
  WRITE(UNIT=KULOUT,FMT='('' YDRI%LSLONDEM = '',L5)')YDRI%LSLONDEM
  WRITE(UNIT=KULOUT,FMT='('' YDRO%LSLONDEM = '',L5)')YDRO%LSLONDEM
  IF( NRADINT > 0 )THEN
    IDIR=LEN_TRIM(CRTABLEDIR)
    IFIL=LEN_TRIM(CRTABLEFIL)
    WRITE(UNIT=KULOUT,FMT='('' CRTABLEDIR = '',A,'' CRTABLEFIL = '',A)')&
     & CRTABLEDIR(1:IDIR),CRTABLEFIL(1:IFIL)  
  ENDIF
  WRITE(UNIT=KULOUT,FMT='('' NO3CMIP = '',I2, '' NGHGCMIP = '',I2, '' NRCP = '',I2 )')&
  &  NO3CMIP, NGHGCMIP,NRCP
  IF( NO3CMIP /= 0  )THEN
    WRITE(UNIT=KULOUT,FMT='('' CO3DATADIR = '',A)') TRIM(CO3DATADIR)
  ENDIF
  IF( NGHGCMIP /= 0 )THEN
    WRITE(UNIT=KULOUT,FMT='('' GHGDATADIR = '',A)') TRIM(GHGDATADIR)
  ENDIF

  IF( NCMIPFIXYR /= -1 )THEN
    WRITE(UNIT=KULOUT,FMT='('' NCMIPFIXYR = '',I8 )') NCMIPFIXYR
  ENDIF

  WRITE(UNIT=KULOUT,FMT='('' LCCNL = '',L5&
   & ,'' LCCNO = '',L5&
   & ,'' RCCNLND= '',F5.0&
   & ,'' RCCNSEA= '',F5.0&
   & ,'' LAERCLIM='',L2&
   & ,'' LAERVISI='',L2&
   & )')&
   & LCCNL,LCCNO,RCCNLND,RCCNSEA,LAERCLIM,LAERVISI
  IF (NGHGRAD /= 0) THEN
    WRITE(UNIT=KULOUT,FMT='('' NGHGRAD= '',I3,'' LETRACGMS= '',L2&
     &)')&
     & NGHGRAD,LETRACGMS
    IF (.NOT. LETRACGMS) THEN
      WRITE(NULOUT,*) 'Warning - old Cariole ozone climatology (LETRACGMS=FALSE) no longer available'
    ENDIF
  ENDIF  
  WRITE(UNIT=KULOUT,FMT='('' HISTORY OF SO4 AEROSOLS, LESO4HIS '',L5)')LESO4HIS
  IF (LESO4HIS) THEN
      WRITE(UNIT=KULOUT,FMT='(''FILE CONTAINING  LIST OF DECADE-WISE DATA: '',A)')TRIM(CLISTSO4)
  ENDIF
  IF (LHVOLCA) THEN
    WRITE(UNIT=KULOUT,FMT='('' HISTORY OF VOLCANIC AEROSOLS= '',L5)')LHVOLCA
    WRITE(UNIT=KULOUT,FMT='('' NVOLCVERT = '',I2)') NVOLCVERT
    IF(NVOLCVERT==3.AND..NOT.LEPO3RA) THEN
      WRITE(UNIT=KULOUT,FMT='('' INCOMPATIBLE WITH LEPO3RA = '',L5)')LEPO3RA
      CALL ABOR1('SUECRAD: ABORT')
    ENDIF
    IF (LVOLCSPEC) THEN
      WRITE(UNIT=KULOUT,FMT='('' SPECIFIED VOLCANIC AEROSOL= '',L5)')LVOLCSPEC
      WRITE(UNIT=KULOUT,FMT='('' DAMPING FROM INITIAL VALUE= '',L5)')LVOLCDAMP
      WRITE(UNIT=KULOUT,FMT='('' RVOLCSPEC = '',3F8.4)') RVOLCSPEC
    ENDIF
    IF (LVOLCDATA) THEN
      WRITE(UNIT=KULOUT,FMT='('' LVOLCDATA= '',L5)')LVOLCDATA
      WRITE(UNIT=KULOUT,FMT='('' READ STRATOSPHERIC AOD FROM CVOLCDATA: '',A)') TRIM(CVOLCDATA)
    ENDIF
  ENDIF
  WRITE(UNIT=KULOUT,FMT='('' LDIAGFORCING =    '',L5)') LDIAGFORCING
  WRITE(UNIT=KULOUT,FMT='('' LApproxLwUpdate = '',L5)') LAPPROXLWUPDATE
  WRITE(UNIT=KULOUT,FMT='('' LApproxSwUpdate = '',L5)') LAPPROXSWUPDATE
  WRITE(UNIT=KULOUT,FMT='('' LMannersSwUpdate= '',L5)') LMANNERSSWUPDATE
  WRITE(UNIT=KULOUT,FMT='('' LCentredTimeSZA = '',L5)') LCENTREDTIMESZA
  WRITE(UNIT=KULOUT,FMT='('' LAverageSZA     = '',L5)') LAVERAGESZA
  WRITE(UNIT=KULOUT,FMT='('' LUsePre2017Rad  = '',L5)') LUSEPRE2017RAD

  ! Only print options for post-2017 schem if it is active
  IF (.NOT. LUSEPRE2017RAD) THEN
    WRITE(UNIT=KULOUT,FMT='('' NLwScattering   = '',I1)') NLWSCATTERING
    WRITE(UNIT=KULOUT,FMT='('' NLwSolver       = '',I1)') NLWSOLVER
    WRITE(UNIT=KULOUT,FMT='('' NSwSolver       = '',I1)') NSWSOLVER
    WRITE(UNIT=KULOUT,FMT='('' NCloudOverlap   = '',I1)') NCLOUDOVERLAP
    WRITE(UNIT=KULOUT,FMT='('' RCLOUD_SEPARATION_SCALE_TOA  = '',F0.1)') &
         &  RCLOUD_SEPARATION_SCALE_TOA
    WRITE(UNIT=KULOUT,FMT='('' RCLOUD_SEPARATION_SCALE_SURF = '',F0.1)') &
         &  RCLOUD_SEPARATION_SCALE_SURF
  ENDIF

  WRITE(UNIT=KULOUT,FMT='('' RCLOUD_FRAC_STD = '',F5.2)') RCLOUD_FRAC_STD
  WRITE(UNIT=KULOUT,FMT='('' LFU_LW_ICE_OPTICS_BUG = '',L5)') LFU_LW_ICE_OPTICS_BUG

  WRITE(UNIT=KULOUT,FMT='('' LONEWSW= '',L5&
   & ,'' NRADIP = '',I1&
   & ,'' NRADLP = '',I1&
   & ,'' NICEOPT= '',I1&
   & ,'' NLIQOPT= '',I1&
   & ,'' LDIFFC = '',L5&
   & )')&
   & LONEWSW,NRADIP,NRADLP,NICEOPT,NLIQOPT,LDIFFC

  WRITE(UNIT=KULOUT,FMT='(&
   &   '' NLWICEOPT= '',I1&
   &   ,'' NLWLIQOPT= '',I1&
   &   ,'' NSWICEOPT= '',I1&
   &   ,'' NSWLIQOPT= '',I1&
   &   ,'' NSW      = '',I1&
   &   ,'' NSWB_MNH = '',I1&
   &   )')&
   & NLWICEOPT,NLWLIQOPT,NSWICEOPT,NSWLIQOPT,NSW,NSWB_MNH
  WRITE(UNIT=KULOUT,FMT='('' NMINICE= '',I1&
   & ,'' NDECOLAT='',I3&
   & )')&
   & NMINICE,NDECOLAT
  WRITE(UNIT=KULOUT,FMT='('' NOVLP   = '',I2)') NOVLP  
  IF (LUVPROC) THEN
    IDAYUV=NUVTIM/INT(RDAY)
    WRITE(UNIT=KULOUT,FMT='('' LUVPROC = '',L5&
   & ,'' LUVTDEP= '',L5&
   & ,'' NRADUV = '',I2&
   & ,'' NUV = '',I2&
   & ,'' NDAYUV = '',I5&
   & ,'' RMUZUV = '',E9.3&
   & ,'' LUVAERP = '',L5&
   & )')&
   & LUVPROC,LUVTDEP,NRADUV,NUV,IDAYUV,RMUZUV,LUVAERP
    WRITE(UNIT=KULOUT,FMT='('' RUVLAM = '',24F6.1)') (RUVLAM(JUV),JUV=1,NUV) 
    WRITE(UNIT=KULOUT,FMT='('' JUVLAM = '',24(3X,I1,2X))') (JUVLAM(JUV),JUV=1,NUV) 
  ENDIF
  WRITE(UNIT=KULOUT,FMT='('' NMCICA= '',I2&
   & ,'' NREDGLW= '',I2&
   & ,'' NREDGSW= '',I2&
   & )')&
   & NMCICA,NREDGLW,NREDGSW

  IF (ALLOCATED(YDESWRT%RSUN))&
   &  WRITE(UNIT=KULOUT,FMT='('' RSUN  FROM YOESW = '',6F10.5)')(YDESWRT%RSUN(I),I=1,NSW)
  IF (ALLOCATED(YDESWRT%RSUN2))&
   &  WRITE(UNIT=KULOUT,FMT='('' RSUN2 FROM YOESW = '',6F10.5)')(YDESWRT%RSUN2(I),I=1,NSW)

ENDIF

! Load total solar irradiance
IF (.NOT. ASSOCIATED(YSOLARIRRADIANCE)) THEN
  ALLOCATE(YSOLARIRRADIANCE)
ENDIF
IF (NHINCSOL /= 0) THEN
  CALL YSOLARIRRADIANCE%READ(TRIM(CSOLARIRRADIANCEFILE))
  IF (NHINCSOL /= 4) THEN
    WRITE(UNIT=KULOUT,FMT='(A)') 'Warning: pretending NHINCSOL == 4'
  ENDIF
ENDIF

! Load greenhouse-gas climatology 
IF (.NOT. ASSOCIATED(YGHGCLIM)) THEN
  ALLOCATE(YGHGCLIM)
ENDIF
IF (NGHGRAD /= 0) THEN
  CALL YGHGCLIM%READ(TRIM(CGHGCLIMFILE))
ENDIF

! Load greenhouse-gas timeseries
IF (.NOT. ASSOCIATED(YGHGTIMESERIES)) THEN
  ALLOCATE(YGHGTIMESERIES)
ENDIF
IF (LHGHG) THEN
  CALL YGHGTIMESERIES%READ(TRIM(CGHGTIMESERIESFILE))
ENDIF

! Initialize the pressure and latitude and allocate the
! time-interpolated gas fields; if they are already initialized, this
! routine does nothing
IF (.NOT. ASSOCIATED(YRADGHG)) THEN
  ALLOCATE(YRADGHG)
ENDIF
IF (NGHGRAD /= 0) THEN
  CALL YRADGHG%INIT(YGHGCLIM%LATITUDE, YGHGCLIM%PRESSURE)
ENDIF

IF (.NOT. LUSEPRE2017RAD .AND. LERADI) THEN
  ! We are not using the pre-2017 radiation scheme framework (although
  ! might be using the pre-2017 implementation of the RRTM-G gas
  ! optics) but rather the ecRad scheme
  CALL SETUP_RADIATION_SCHEME(YDERDI,YDEAERATM,YDCOMPO,YDEPHY,YDERAD,YDRADIATION,LOUTPUT)
ENDIF

! setup of YOEOZOC arrays moved out from the call in UPDTIM to an actual setup routine, called here
IF (NOZOCL==1) THEN
  IF (LAPE) THEN
    CALL SUECOZCAQUA(YDRIP%YREOZOC)
  ELSEIF (YDRIP%YRECMIP%NO3CMIP /= 0) THEN
    CALL SUECOZV(YDRIP%YRECMIP,NINDAT)
  ELSE 
    CALL SUECOZC(YDRIP%YREOZOC)
  ENDIF
ENDIF

! Aerosol climatology Tegen et al. / GISS Volcanic aerosol climatology 
! and  O.Boucher sulphate history (obs 1920-1990, A1B scenario 2000-2100)
IF (LHVOLCA .OR. LNEWAER) THEN
!!  CALL YDRIP%YREAERC_TEGEN%SETUP(LCALL_SUECAEC = (.NOT. YDMODEL%YRML_PHY_SLIN%YREPHLI%LPHYLIN) )
  CALL YDRIP%YREAERC_TEGEN%SETUP(LCALL_SUECAEC = (.NOT. YDMODEL%YRML_PHY_SLIN%YREPHLI%LPHYLIN) .OR. L_OOPS )
  WRITE(NULOUT,FMT='(" SUECRAD TEGEN AEROSOL CLIMATOLOGY IS LOADED")')    

  IF (NAERMACC == 1) THEN
    CALL YDRIP%YREAERC_MACC%SETUP(LAER3D)
  ENDIF
  
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUECRAD',1,ZHOOK_HANDLE)
END SUBROUTINE SUECRAD
