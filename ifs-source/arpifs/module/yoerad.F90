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
MODULE YOERAD

USE PARKIND1,            ONLY : JPIM, JPRB
USE YOE_SPECTRAL_PLANCK, ONLY : TSPECTRALPLANCK
USE YOMHOOK, ONLY : LHOOK, DR_HOOK, JPHOOK

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------
!*    ** *YOERAD* - CONTROL OPTIONS FOR RADIATION CONFIGURATION
!     ------------------------------------------------------------------

TYPE :: TERAD
INTEGER(KIND=JPIM) :: NAER
INTEGER(KIND=JPIM) :: NMODE
INTEGER(KIND=JPIM) :: NOZOCL
INTEGER(KIND=JPIM) :: NRADFR
INTEGER(KIND=JPIM) :: NRADPFR
INTEGER(KIND=JPIM) :: NRADPLA
INTEGER(KIND=JPIM) :: NRADINT
INTEGER(KIND=JPIM) :: NRADRES
INTEGER(KIND=JPIM) :: NRADNFR
INTEGER(KIND=JPIM) :: NRADSFR
INTEGER(KIND=JPIM) :: NRADE1H, NRADE3H
INTEGER(KIND=JPIM) :: NRADELG
INTEGER(KIND=JPIM) :: NOVLP
INTEGER(KIND=JPIM) :: NRPROMA
INTEGER(KIND=JPIM) :: NSW
INTEGER(KIND=JPIM) :: NSWNL
INTEGER(KIND=JPIM) :: NSWTL
INTEGER(KIND=JPIM) :: NLWEMISS
INTEGER(KIND=JPIM) :: NLWOUT
INTEGER(KIND=JPIM) :: NTSW
INTEGER(KIND=JPIM) :: NUV
INTEGER(KIND=JPIM) :: NCSRADF
INTEGER(KIND=JPIM) :: NICEOPT
INTEGER(KIND=JPIM) :: NLIQOPT
INTEGER(KIND=JPIM) :: NSWICEOPT
INTEGER(KIND=JPIM) :: NSWLIQOPT
INTEGER(KIND=JPIM) :: NLWICEOPT
INTEGER(KIND=JPIM) :: NLWLIQOPT
INTEGER(KIND=JPIM) :: NRADIP
INTEGER(KIND=JPIM) :: NRADLP
INTEGER(KIND=JPIM) :: NINHOM
INTEGER(KIND=JPIM) :: NLAYINH
INTEGER(KIND=JPIM) :: NLNGR1H
INTEGER(KIND=JPIM) :: NPERTAER
INTEGER(KIND=JPIM) :: NPERTOZ
INTEGER(KIND=JPIM) :: NSCEN
INTEGER(KIND=JPIM) :: NHINCSOL
INTEGER(KIND=JPIM) :: NMCICA
INTEGER(KIND=JPIM) :: NGHGRAD
INTEGER(KIND=JPIM) :: NDECOLAT
INTEGER(KIND=JPIM) :: NMINICE
INTEGER(KIND=JPIM) :: NVOLCVERT
INTEGER(KIND=JPIM) :: NREDGLW
INTEGER(KIND=JPIM) :: NREDGSW
INTEGER(KIND=JPIM) :: NAERMACC, NMCLAT, NMCLON, NMCLEV, NMCVAR
INTEGER(KIND=JPIM) :: NSPMAPL(16), NSPMAPS(14)
INTEGER(KIND=JPIM) :: NLWSCATTERING
INTEGER(KIND=JPIM) :: NSWSOLVER, NLWSOLVER
INTEGER(KIND=JPIM) :: KMODTS
INTEGER(KIND=JPIM) :: NSOLARSPECTRUM
INTEGER(KIND=JPIM) :: NSWWVCONTINUUM
INTEGER(KIND=JPIM) :: NDUMPBADINPUTS
INTEGER(KIND=JPIM) :: NDUMPINPUTS
INTEGER(KIND=JPIM) :: NCLOUDOVERLAP
REAL(KIND=JPRB)    :: RCLOUD_FRAC_STD
REAL(KIND=JPRB)    :: RCLOUD_SEPARATION_SCALE_TOA, RCLOUD_SEPARATION_SCALE_SURF
LOGICAL :: LFU_LW_ICE_OPTICS_BUG
LOGICAL :: LINTERPINCLOUDMEAN

LOGICAL :: LERAD1H
LOGICAL :: LEPO3RA
LOGICAL :: LONEWSW
LOGICAL :: LECSRAD
LOGICAL :: LRRTM
LOGICAL :: LSRTM
LOGICAL :: LDIFFC
LOGICAL :: LHVOLCA
LOGICAL :: LNEWAER
LOGICAL :: LNOTROAER
LOGICAL :: LRAYL
LOGICAL :: LOPTRPROMA
LOGICAL :: LECO2VAR
LOGICAL :: LHGHG
LOGICAL :: LESO4HIS
LOGICAL :: LETRACGMS
LOGICAL :: LAERCLIM, LAERVISI, LAER3D
LOGICAL :: LVOLCSPEC
LOGICAL :: LVOLCDAMP
LOGICAL :: LDIAGFORCING
LOGICAL :: LAERADCLI
LOGICAL :: LAERADJDU
LOGICAL :: LAPPROXLWUPDATE
LOGICAL :: LAPPROXSWUPDATE
LOGICAL :: LMANNERSSWUPDATE
LOGICAL :: LCENTREDTIMESZA
LOGICAL :: LAVERAGESZA
LOGICAL :: LECOMPGRID
LOGICAL :: LUSEPRE2017RAD
LOGICAL :: LDUSEASON

LOGICAL :: LCCNL
LOGICAL :: LCCNO
LOGICAL :: LPERPET

REAL(KIND=JPRB) :: RAOVLP , RBOVLP
REAL(KIND=JPRB) :: RCCNLND, RCCNSEA
LOGICAL :: LEDBUG
REAL(KIND=JPRB) :: RPERTOZ, RRE2DE
REAL(KIND=JPRB) :: RLWINHF, RSWINHF
REAL(KIND=JPRB) :: RMINICE
REAL(KIND=JPRB) :: RVOLCSPEC(3)
REAL(KIND=JPRB) :: RNS, RSIGAIR

REAL(KIND=JPRB) :: RAESHSS, RAESHDU, RAESHOM, RAESHBC, RAESHSU
REAL(KIND=JPRB) :: TRBKG
REAL(KIND=JPRB) :: STBKG
CHARACTER(LEN=256) :: CGHGCLIMFILE = ' '
CHARACTER(LEN=256) :: CGHGTIMESERIESFILE = ' '
CHARACTER(LEN=256) :: CSOLARIRRADIANCEFILE = ' '

REAL(KIND=JPRB),ALLOCATABLE:: CVDAESS(:)
REAL(KIND=JPRB),ALLOCATABLE:: CVDAEDU(:)
REAL(KIND=JPRB),ALLOCATABLE:: CVDAEOM(:)
REAL(KIND=JPRB),ALLOCATABLE:: CVDAEBC(:)
REAL(KIND=JPRB),ALLOCATABLE:: CVDAESU(:)

! Look-up table for Planck function in emissivity intervals
TYPE(TSPECTRALPLANCK) :: YSPECTPLANCK
!----------------------------------------------------------------------------
CONTAINS
  PROCEDURE, PASS :: PRINT => PRINT_CONFIGURATION 

END TYPE TERAD
!============================================================================

TYPE(TERAD), POINTER :: YRERAD => NULL()


!        * E.C.M.W.F. PHYSICS PACKAGE *

!     J.-J. MORCRETTE       E.C.M.W.F.      89/07/14
! Modifications
!    R J Hogan 20 May  2014: Added LApproxLwUpdate
!    R J Hogan 19 June 2014: Added LApproxSwUpdate
!    R J Hogan 19 Nov  2014: Added LCentredTimeSZA
!    R J Hogan 15 Apr  2015: Added LMannersSwUpdate
!    R J Hogan 24 Apr  2015: Added LAverageSZA
!    R J Hogan 18 Sept 2015: Added LUsePre2017Rad (was LUsePre2015Rad)
!    R J Hogan 1  Mar  2016: Added NLwSolver, NSwSolver, NLwScattering
!    A Bozzo      Feb  2017: Added logical to enable 3D aerosol climatology
!    R J Hogan 9  Mar  2018: Added NDUMPBADINPUTS
!    R J Hogan 22 Jan  2019: Added NLWEMISS, NCLOUDOVERLAP, NDUMPINPUTS
!    R J Hogan 4  Feb  2019: Added NLWOUT
!    R J Hogan 5  Feb  2019: Added YSPECTPLANCK
!    R J Hogan 11 Mar  2019: Added CGHG*FILE, CSOLARIRRADIANCEFILE

!  NAME     TYPE     PURPOSE
!  ----  :  ----   : ---------------------------------------------------
! LERAD1H: LOGICAL : .T. TO ALLOW MORE FREQUENT RADIATION CALCULATIONS
!                  : DURING FIRST N HOURS OF FORECAST
! NLNGR1H: INTEGER : NUMBER FORECAST HOURS DURING WHICH MORE FREQUENT
!                    RADIATION CALCULATIONS ARE REQUIRED
! LEPO3RA: LOGICAL : .T. IF PROGNOSTIC OZONE (EC) IS PASSED TO RADIATION
! NAER   : INTEGER : CONFIGURATION INDEX FOR AEROSOLS
! NMODE  : INTEGER : CONFIGURATION FOR RADIATION CODE: FLUX VS. RADIANCE
! NOZOCL : INTEGER : CHOICE OF OZONE CLIMATOLOGY (0 old, 1 new)
! NRADFR : INTEGER : FREQUENCY OF FULL RADIATION COMPUTATIONS
!                    IF(NRADFR.GT.0): RAD EVERY 'NRADFR' TIME-STEPS
!                    IF(NRADFR.LT.0): RAD EVERY '-NRADFR' HOURS
! NRADPFR: INTEGER : PRINT FREQUENCY FOR RAD.STATISTICS (in RAD.T.STEPS)
! NRADPLA: INTEGER : PRINT RAD.STATISTICS EVERY 'NRADPLA' ROWS
! NRADINT: INTEGER : RADIATION INTERPOLATION METHOD
!                  : 1 = SPECTRAL TRANSFORM INTERPOLATION
!                  : 2 =  4 POINT HORIZONTAL INTERPOLATION
!                  : 3 = 12 POINT HORIZONTAL INTERPOLATION
! NRADRES: INTEGER : RADIATION GRID SPECTRAL RESOLUTION
! NRADNFR: INTEGER : NORMAL   FREQUENCY OF RADIATION STEPS
! NRADSFR: INTEGER : START-UP FREQUENCY OF RADIATION STEPS
! NRADE1H: INTEGER : START-UP FREQUENCY OF RADIATION STEPS FOR EPS
! NRADE3H: INTEGER : SUBSEQUENT FREQUENCY OF RADIATION STEPS FOR EPS
! NRADELG: INTEGER : LENGTH IN HOURS DURING WHICH THE FREQUENCY OF RADIATION IS INCREASED FOR EPS
! NOVLP  : INTEGER : CLOUD OVERLAP CONFIGURATION IN PRE-MCRAD/ECRAD SCHEME
!                  : 1 = Max-rand (Geleyn & Hollingsworth)
!                  : 2 = Maximum
!                  : 3 = Random
! NRPROMA: INTEGER : VECTOR LENGTH FOR RADIATION CALCULATIONS
! NSW    : INTEGER : NUMBER OF SHORTWAVE SPECTRAL INTERVALS
! NSWNL  : INTEGER : NUMBER OF SHORTWAVE SPECTRAL INTERVALS IN NL MODEL
! NSWTL  : INTEGER : NUMBER OF SHORTWAVE SPECTRAL INTERVALS IN TL MODEL
! NTSW   : INTEGER : MAXIMUM POSSIBLE NUMBER OF SW SPECTRAL INTERVALS 
! NUV    : INTEGER : NUMBER OF UV SPECTRAL INTERVALS FOR THE UV PROCESSOR   
! LOPTRPROMA:LOGICAL: .T. NRPROMA will be optimised
!                   : .F. NRPROMA will not be optimised (forced
!                   :         by negative NRPROMA in namelist)

! NRADIP : INTEGER : INDEX FOR DIAGNOSIS OF ICE CLOUD EFFECTIVE RADIUS
!          0 = fixed at 40 microns
!          1 = Liou & Ou (1994) capped between 40-130 microns
!          2 = Liou & Ou but capped between 30 and 60 microns
!          3 = Sun & Rikus (1999) revised by Sun (2001)
! NRADLP : INTEGER : INDEX FOR DIAGNOSIS OF LIQ. CLOUD EFFECTIVE RADIUS
!          0 = ERA-15 function of pressure
!          1 = 10 microns over land, 13 microns over sea
!          2 = Martin_et_al (1994) in terms of land-sea number conc
!          3 = Linked to prognostic aerosols
! NICEOPT: INTEGER : INDEX FOR ICE CLOUD OPTICAL PROPERTIES
!          0 = SW Ebert-Curry, LW Smith & Shi (1992)
!          1 = SW Ebert-Curry, LW Ebert-Curry (1992)
!          2 = SW & LW Fu-Liou (1993)
!          3 = SW Fu (1996) LW Fu et al. (1998) + Chou et al. (1999) LW scatt approx
!   the following only available in newer modular radiation scheme:
!          4 = SW/LW Baran data fitted versus ice mixing ratio
! NLIQOPT: INTEGER : INDEX FOR LIQUID WATER CLOUD OPTICAL PROPERTIES
!          0 = SW Fouquart (1991) LW Smith-Shi (1992) YF/SmSh 
!          1 = SW Slingo (1989) LW Savijarvi (1997)
!          2 = SW Slingo (1989) LW Lindner-Li (2000)
!   the following only available in RADLSW, not RADLSWR:
!          3 = SW Nielsen       LW Smith-Shi
!   the following only available in newer modular radiation scheme:
!          4 = SW/LW SOCRATES scheme
!
! LONEWSW: LOGICAL : .T. IF NEW SW CODE IS ACTIVE
! LECSRAD: LOGICAL : .T. IF CLEAR-SKY RADIATION IS ARCHIVED AS PEXTR2
! NCSRADF: INTEGER : 1 IF ACCUMULATED, 2 IF INSTANTANEOUS
! LRRTM  : LOGICAL : .T. IF RRTM140MR IS USED FOR LW RADIATION TRANSFER

! LHVOLCA: LOGICAL : .T. IF USING HISTORICAL VOLCANIC AEROSOLS 
! LNEWAER: LOGICAL : .T. IF AEROSOL MONTHLY DISTRIBUTIONS ARE USED
! LNOTROAER:LOGICAL: .T. IF NO TROPOSPHERIC AEROSOLS
! CRTABLEDIR: CHAR : IF NRADINT > 0 SPECIFIES DIRECTORY PATH FOR RADIATION
!                  : GRID RTABLE NAMELIST
! CRTABLEFIL: CHAR : IF NRADINT > 0 SPECIFIES FILE NAME OF RADIATION 
!                  : GRID RTABLE NAMELIST
! LRAYL  : LOGICAL : .T. NEW RAYLEIGH FOR SW-6 VERSION

! RAOVLP : REAL    : COEFFICIENTS FOR ALPHA1 FACTOR IN HOGAN & 
! RBOVLP : REAL    : ILLINGWORTH's PARAMETRIZATION

! LCCNL  : LOGICAL : .T. IF CCN CONCENTRATION OVER LAND IS DIAGNOSED
! LCCNO  : LOGICAL : .T. IF CCN CONCENTRATION OVER OCEAN IS DIAGNOSED
! RCCNLND: REAL    : NUMBER CONCENTRATION (CM-3) OF CCNs OVER LAND
! RCCNSEA: REAL    : NUMBER CONCENTRATION (CM-3) OF CCNs OVER SEA

! LDIFFC : LOGICAL : .T. IF SAVIJARVI'S DIFFUSIVITY CORRECTION IS ON

! NINHOM : INTEGER : 0 IF NO INHOMOGENEITY SCALING EFFECT 
!                    1 IF SIMPLE 0.7 SCALING
!                    2 IF BARKER, 3 IF CAIRNS ET AL.
! RLWINHF: REAL    : INHOMOG. SCALING FACTOR FOR CLOUD LW OPTICAL THICKNESS
! RSWINHF: REAL    : INHOMOG. SCALING FACTOR FOR CLOUD SW OPTICAL THICKNESS

! NPERTAER : INTERGER : PERCENTAGE OF PERTURBATION FOR AEROSOL   
! NPERTOZONE : INTEGER : PERCENTAGE OF PERTURBATION FOR OZONE 
! NHINCSOL : INTEGER :  0: Total Solar Irradiance (TSI) fixed at 1366.0 W m-2
!                       1: Deprecated - use default
!                       2: Deprecated - use default
!                       3: Deprecated (was CMIP5) - use default
!                       4: TSI from CMIP6 NetCDF (default), or override with CSOLARIRRADIANCEFILE
! NSWWVCONTINUUM : INTEGER : 0 MT_CKD2.5 (SRTM default WV continuum)
!                            1 CAVIAR continuum (Shine et al. 2016)
! LECO2VAR: LOGICAL: .T. IF ERA-40/AMIP2 VARIABILITY OF GHG IS ON (ignored)
! LHGHG  : LOGICAL : .T. IF VARIABILITY OF GREENHOUSE GASES (INCLUDING CO2) IS ON
! N.B.: LHGHG supersedes LECO2VAR and allows using better specification of trace gases
! NSCEN  : INTEGER : 21st CENTURY SCENARIO FOR GHG (1=A1B, 2=A2, 3=B1)
! RRe2De : REAL    : CONVERSION FACTOR BETWWEN EFFECTIVE RADIUS AND PARTICLE SIZE
! RMINICE: REAL    : MINIMUM SIZE FOR ICE PARTICLES (um)
!                    FOR ICE
! NMINICE: INTEGER : 1-6 MINIMUM ICE PARTICLE SIZE DEPENDS ON LATITUDE, 0=INDEPENDENT OF LATITUDE
! NDECOLAT:INTEGER : DECORRELATION LENGTH FOR CF AND CW 
!                     0: SPECIFIED INDEPENDENT OF LATITUDE, 1: SHONK-HOGAN, 2: IMPROVED
! NMCICA : INTEGER :  0: NO McICA
!                     1: McICA w maximum-random in cloud generator
!                     2: McICA w generalized overlap in cloud generator
! LESO4HIS: LOGICAL:.T.: Use historical/projected SO4 data per decade and month
! NGHGRAD: INTEGER : configuration of 3D GHG climatologies accounted for in radiation
!                     0: global values
!                     1: CO2       2: CH4    3: N2O    4: NO2    5:CFC11   6:CFC12
!                    12: CO2+CH4  13: CO2+CH4+N2O     
!                    16: CO2+CH4+N2O+CFC11+CFC12
! LETRACGMS: LOGICAL : F=Cariolle climatol. T=GEMS-derived clim for CO2, CH4, O3
! LAERCLIM : LOGICAL : .T. for output of the climatological aerosol optical depth at 550 nm
! LAERVISI : LOGICAL : .T. for output of the visibility (from diagnsotic or prognostic aerosols)
! NVOLCVERT: INTEGER : Vertical distribution of volcanic aerosol
!                       0: original profile, diagnosed from T
!                       1: original profile, but upper boundary at 10hPa
!                       2: lower boundary diagnosed from ozone, upper boundary at 10hPa
! LVOLCSPEC: LOGICAL : T for specified volcanic aerosol
! LVOLCDAMP: LOGICAL : T for damping of specified volcanic aerosol from initial value
! RVOLCSPEC: REAL    : Specified volcanic aerosol (total optical depth) in NH/Tropics/SH
! RNs                : derived from Avogadro
! RSIGAIR: invariant terms in expression of Rayleigh scattering cross-section
! NREDGSW  : INTEGER : 0 full resolution for RRTM_SW (224)
!                      1 ECMWF High resolution model configuration (_SW: 112)
!                      2 ECMWF EPS configuration (_SW: 56)
! NREDGLW  : INTEGER : 0 full resolution for RRTM_LW (256)
!                      1 ECMWF High resolution model configuration (_LW: 140)
!                      2 ECMWF EPS configuration (_LW: 70)
! LDIAGFORCING : LOGICAL : T Write input ozone, ghg and aerosol forcing to 3D fields 
!                            To be used for diagnostics only; do not use in production runs
! NAERMACC : INTEGER : MACC-derived aerosol climatology on a NMCLAT x NMCLON grid
! RAESHxx  : REAL    : parameters related to scale height of MACC-derived aerosol climatology
! CVDAExx  : REAL    : scale heights of MACC-derived aerosol climatology
! LAERADJDU: LOGICAL : T adjust MACC-derived DU climatology
! LAERADCLI: LOGICAL : T if radiation uses the MACC-derived aerosol climatology
! LApproxLwUpdate : LOGICAL : Update the longwave upwelling flux every
!                             timestep/gridpoint using the stored rate
!                             of change of the fluxes with respect to
!                             the surface upwelling longwave flux
! LApproxSwUpdate : LOGICAL : Update the shortwave upwelling flux
!                             every gridpoint to account for the local
!                             value of surface albedo
! LMannersSwUpdate: LOGICAL : Update the shortwave flux every timestep
!                             using Manners et al. (2009) correction
!                             for solar zenith angle change
! LCentredTimeSZA : LOGICAL : Compute solar zenith angle in radiation
!                             scheme half way between calls to
!                             radiation scheme (rather than previous
!                             behaviour, which is half way between
!                             calls plus half a model timestep)
! LAverageSZA     : LOGICAL : Compute an averaged solar zenith angle
!                             across the time interval required
!                             (either a model timestep or a radiation
!                             timestep). Should be used with 
!                             LCentredTimeSZA=TRUE.
! LUsePre2017Rad  : LOGICAL : Use the pre-2017 radiation scheme, rather 
!                             than the modular scheme contained in the
!                             separate "radiation" library.  Note that
!                             the radiation library may make use of the
!                             pre-2017 RRTM-G gas optics.
! RCloud_Frac_Std : REAL    : Cloud water content horizontal fractional
!                             standard deviation in a gridbox
! LInterpInCloudMean : LOGICAL : When interpolating model fields to
!                             radiation grid, interpolate in-cloud
!                             mean water contents?  Better
!                             conservation achieved by interpolating
!                             gridbox-means.
! CGHGCLIMFILE : STRING :     Location of greenhouse gas climatology file,
!                             or empty if the default is to be used. If it
!                             starts with "." or "/" then a relative path
!                             is assumed, otherwise the default directory.
! CGHGTIMESERIESFILE:STRING : Location of greenhouse gas timeseries
!                             file, or empty if it is to be worked out
!                             from the NSCEN, YOECMIP%NGHGCMIP and
!                             YOECMIP%NRCP variables. If it starts
!                             with "." or "/" then a relative path is
!                             assumed, otherwise the default
!                             directory.
! CSOLARIRRADIANCEFILE:STRING:Location of Total Solar Irradiance file,
!                             or empty if the default is to be
!                             used. If it starts with "." or "/" then
!                             a relative path is assumed, otherwise
!                             the default directory.
! NLWEMISS      : INTEGER :   Number of emissivity spectral intervals, set 
!                             according to the value of NEMISSSCHEME; traditionally
!                             this has always been 2: outside the IR window and within
! NLWOUT        : INTEGER :   Number of spectral intervals to pass LW downwelling flux
!                             to RADHEATN; traditionally this was 1, but this led
!                             to errors with LAPPROXLWUPDATE=TRUE, which updated 
!                             fluxes using a single broadband emissivity. Now we can
!                             do approximate updates using full spectral emissivity.
! ------------------------------------------------------------------
! THE FOLLOWING ARE ONLY USED FOR THE ECRAD SCHEME (LUsePre2017Rad = .FALSE.)
! NLwScattering   : INTEGER : 0: No longwave scattering
!                             1: Longwave scattering by clouds only
!                             2: Longwave scattering by clouds and aerosols
! NSwSolver       : INTEGER :
! NLwSolver       : INTEGER : 0: McICA
!                             1: SPARTACUS 1D
!                             2: SPARTACUS 3D
!                             3: TripleClouds
! LFU_LW_ICE_OPTICS_BUG : LOGICAL : Continue to use bug in Fu LW ice
!                             optics whereby single scattering albedo is
!                             one minus what it should be
! NSOLARSPECTRUM : INTEGER :  0: Kurucz
!                             1: WHI reference solar minimum spectrum
!                             2: Mean of Coddington et al. (BAMS 2016) spectrum
! NDUMPBADINPUTS : INTEGER :  0: Warn only if fluxes out of physical bounds
!                             n: Write netcdf file of bad inputs up to n times per task
!                            -n: Abort if fluxes ever out of physical bounds
! NDUMPINPUTS    : INTEGER :  0: Do nothing
!                             n: Write netcdf file of all inputs up to n times per task
! NCLOUDOVERLAP  : INTEGER :  Cloud overlap scheme
!                             1: Maximum-random
!                             2: Exponential-exponential (the actual behaviour of McRad)
!                             3: Exponential-random (only option for Tripleclouds and SPARTACUS)
! RCLOUD_SEPARATION_SCALE_TOA, RCLOUD_SEPARATION_SCALE_SURF : REAL
!                 Cloud horizontal length scale, in metres, used to
!                 compute rate of horizontal exchange of radiation
!                 between clouds and clear skies in SPARTACUS solver
! ------------------------------------------------------------------
! KMODTS : INTEGER   : (A Bozzo) switch for different radiative transfer schemes for UV 
!                       = 0 Fouquart&Bonnel adapted by Morcrette and Arola
!                       = 1 eddington (joseph et al., 1976)
!                       = 2 pifm (zdunkowski et al., 1980)
!                       = 3 discrete ordinates (liou, 1973)
!     ------------------------------------------------------------------
! TRBKG : REAL tropospheric background OD@550nm for aerosol climatology.
!                  default for Tegen climatology was 0.03
! STBKG : REAL stratospheric background OD@550nm for aerosol climatology.
!     ------------------------------------------------------------------
! LDUSEASON : LOGICAL enables a monthly-varying scale height for the 
!                     dust aerosol climatology
! LAER3D : LOGICAL : to enable aerosol climatology in 3D


CONTAINS

SUBROUTINE PRINT_CONFIGURATION(SELF, KDEPTH, KOUTNO)

IMPLICIT NONE
CLASS(TERAD), INTENT(IN) :: SELF
INTEGER     , INTENT(IN) :: KDEPTH
INTEGER     , INTENT(IN) :: KOUTNO

INTEGER :: IDEPTHLOC
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('YOERAD:PRINT_CONFIGURATION',0,ZHOOK_HANDLE)
IDEPTHLOC = KDEPTH+2

WRITE(KOUTNO,*) REPEAT(' ',KDEPTH   ) // 'model%yrml_phy_rad%yrerad : '
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NAER = ', SELF%NAER
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NMODE = ', SELF%NMODE
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NOZOCL = ', SELF%NOZOCL
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NRADFR = ', SELF%NRADFR
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NRADPFR = ', SELF%NRADPFR
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NRADPLA = ', SELF%NRADPLA
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NRADINT = ', SELF%NRADINT
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NRADRES = ', SELF%NRADRES
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NRADNFR = ', SELF%NRADNFR
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NRADSFR = ', SELF%NRADSFR
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NRADE1H = ', SELF%NRADE1H
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NRADE3H = ', SELF%NRADE3H
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NRADELG = ', SELF%NRADELG
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NOVLP = ', SELF%NOVLP
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NRPROMA = ', SELF%NRPROMA
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NSW = ', SELF%NSW
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NSWNL = ', SELF%NSWNL
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NSWTL = ', SELF%NSWTL
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NTSW = ', SELF%NTSW
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NLWEMISS = ', SELF%NLWEMISS
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NLWOUT = ', SELF%NLWOUT
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NUV = ', SELF%NUV
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NCSRADF = ', SELF%NCSRADF
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NICEOPT = ', SELF%NICEOPT
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NLIQOPT = ', SELF%NLIQOPT
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NSWICEOPT = ', SELF%NSWICEOPT
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NSWLIQOPT = ', SELF%NSWLIQOPT
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NLWICEOPT = ', SELF%NLWICEOPT
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NLWLIQOPT = ', SELF%NLWLIQOPT
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NRADIP = ', SELF%NRADIP
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NRADLP = ', SELF%NRADLP
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NINHOM = ', SELF%NINHOM
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NLAYINH = ', SELF%NLAYINH
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NLNGR1H = ', SELF%NLNGR1H
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NPERTAER = ', SELF%NPERTAER
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NPERTOZ = ', SELF%NPERTOZ
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NSCEN = ', SELF%NSCEN
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NHINCSOL = ', SELF%NHINCSOL
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NMCICA = ', SELF%NMCICA
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NGHGRAD = ', SELF%NGHGRAD
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NDECOLAT = ', SELF%NDECOLAT
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NMINICE = ', SELF%NMINICE
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NVOLCVERT = ', SELF%NVOLCVERT
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NREDGLW = ', SELF%NREDGLW
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NREDGSW = ', SELF%NREDGSW
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NAERMACC = ', SELF%NAERMACC
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NMCLAT = ', SELF%NMCLAT
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NMCLON = ', SELF%NMCLON
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NMCLEV = ', SELF%NMCLEV
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NMCVAR = ', SELF%NMCVAR
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NSPMAPL sum = ',SUM(SELF%NSPMAPL)
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NSPMAPS sum =', SUM(SELF%NSPMAPS)
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NLWSCATTERING = ', SELF%NLWSCATTERING
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NSWSOLVER = ', SELF%NSWSOLVER
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NLWSOLVER = ', SELF%NLWSOLVER
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'KMODTS = ', SELF%KMODTS
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RCLOUD_FRAC_STD = ', SELF%RCLOUD_FRAC_STD
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LFU_LW_ICE_OPTICS_BUG = ', SELF%LFU_LW_ICE_OPTICS_BUG
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LERAD1H = ', SELF%LERAD1H
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LEPO3RA = ', SELF%LEPO3RA
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LONEWSW = ', SELF%LONEWSW
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LECSRAD = ', SELF%LECSRAD
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LRRTM = ', SELF%LRRTM
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LSRTM = ', SELF%LSRTM
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LDIFFC = ', SELF%LDIFFC
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LHVOLCA = ', SELF%LHVOLCA
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LNEWAER = ', SELF%LNEWAER
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LNOTROAER = ', SELF%LNOTROAER
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LRAYL = ', SELF%LRAYL
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LOPTRPROMA = ', SELF%LOPTRPROMA
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LECO2VAR = ', SELF%LECO2VAR
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LHGHG = ', SELF%LHGHG
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LESO4HIS = ', SELF%LESO4HIS
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LETRACGMS = ', SELF%LETRACGMS
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LAERCLIM = ', SELF%LAERCLIM
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LAERVISI = ', SELF%LAERVISI
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LVOLCSPEC = ', SELF%LVOLCSPEC
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LVOLCDAMP = ', SELF%LVOLCDAMP
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LDIAGFORCING = ', SELF%LDIAGFORCING
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LAERADCLI = ', SELF%LAERADCLI
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LAERADJDU = ', SELF%LAERADJDU
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LAPPROXLWUPDATE = ', SELF%LAPPROXLWUPDATE
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LAPPROXSWUPDATE = ', SELF%LAPPROXSWUPDATE
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LMANNERSSWUPDATE = ', SELF%LMANNERSSWUPDATE
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LCENTREDTIMESZA = ', SELF%LCENTREDTIMESZA
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LAVERAGESZA = ', SELF%LAVERAGESZA
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LECOMPGRID = ', SELF%LECOMPGRID
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LUSEPRE2017RAD = ', SELF%LUSEPRE2017RAD
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LDUSEASON = ', SELF%LDUSEASON
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LCCNL = ', SELF%LCCNL
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LCCNO = ', SELF%LCCNO
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LPERPET = ', SELF%LPERPET
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RAOVLP = ', SELF%RAOVLP
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RBOVLP = ', SELF%RBOVLP
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RCCNLND = ', SELF%RCCNLND
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RCCNSEA = ', SELF%RCCNSEA
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LEDBUG = ', SELF%LEDBUG
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RPERTOZ = ', SELF%RPERTOZ
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RRE2DE = ', SELF%RRE2DE
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RLWINHF = ', SELF%RLWINHF
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RSWINHF = ', SELF%RSWINHF
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RMINICE = ', SELF%RMINICE
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RVOLCSPEC sum = ', SUM(SELF%RVOLCSPEC)
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RNS = ', SELF%RNS
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RSIGAIR = ', SELF%RSIGAIR
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RAESHSS = ', SELF%RAESHSS
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RAESHDU = ', SELF%RAESHDU
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RAESHOM = ', SELF%RAESHOM
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RAESHBC = ', SELF%RAESHBC
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RAESHSU = ', SELF%RAESHSU
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'TRBKG = ', SELF%TRBKG
IF (ALLOCATED(SELF%CVDAESS)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'CVDAESS ALLOCATED OF SHAPE ', &
 &        SHAPE(SELF%CVDAESS), ' SUM ',SUM(SELF%CVDAESS)
IF (ALLOCATED(SELF%CVDAEDU)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'CVDAEDU ALLOCATED OF SHAPE ', &
 &        SHAPE(SELF%CVDAEDU), ' SUM ',SUM(SELF%CVDAEDU)
IF (ALLOCATED(SELF%CVDAEOM)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'CVDAEOM ALLOCATED OF SHAPE ', &
 &        SHAPE(SELF%CVDAEOM), ' SUM ',SUM(SELF%CVDAEOM)
IF (ALLOCATED(SELF%CVDAEBC)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'CVDAEBC ALLOCATED OF SHAPE ', &
 &        SHAPE(SELF%CVDAEBC), ' SUM ',SUM(SELF%CVDAEBC)
IF (ALLOCATED(SELF%CVDAESU)) WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'CVDAESU ALLOCATED OF SHAPE ', &
 &        SHAPE(SELF%CVDAESU), ' SUM ',SUM(SELF%CVDAESU)
IF (LHOOK) CALL DR_HOOK('YOERAD:PRINT_CONFIGURATION',1,ZHOOK_HANDLE)

END SUBROUTINE PRINT_CONFIGURATION

END MODULE YOERAD
