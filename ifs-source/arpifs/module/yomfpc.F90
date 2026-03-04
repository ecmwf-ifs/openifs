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

MODULE YOMFPC

USE PARKIND1  , ONLY : JPIM     ,JPRB
USE YOMCT0    , ONLY : JPNPST

USE PARFPOS, ONLY : JPOSDIR, JPOSLEN  ,JPOSDOM, JPOS2DF  ,JPOS3DF  ,JPOS3H   ,JPOS3P   ,JPOS3PV  , &
 & JPOS3S   ,JPOS3TH  ,JPOSCFU  ,JPOSSGP  ,JPOSXFU  ,JPOS3I   ,JPOS3F

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------

! === TECHNICAL CONSTANTS ===

!     LTRACEFP: trace for Full-POS computations
!     LALLOFP : trace for Full-POS allocations/deallocations
!     LFPMOIS: Month contolled against climatology :
!     L_READ_MODEL_DATE: if: .TRUE. read date from the model
!               .F. => month of the model (forecast)
!               .T. => month of the file
!     LFPPACKING : Enable/disable packing/unpacking of model fields before doing inline post-processing
!     NFPCHKDAT  : Check the consistency of the date in pp files (1) or not (0)
!     NFPSFXWRT  : Write Surfex PREP fields after interpolations (0) ; or write interpolated fields directly (1)
!     LWIDER_DOM  : .TRUE. if the output domain is allowed to be larger than the input one
!     RWIDOM     : Global missing value


! Parameters needed to control different objects simultaneously

!     CFPDIR : prefix for the output files
!     CFPFMT : format of the output files ('MODEL','GAUSS','LELAM' or 'LALON')
!     CFPDOM : names of the subdomains
!     CFPMONIPATH_IN  : pathname of control files in input (allows to support multiple objects)
!     CFPMONIPATH_OUT : pathname of control file in output (allows to support multiple objects)
!     NFPGRIB    : level of GRIB coding
!                0 : no packing at all
!                1 : standard GRIB-0 encoding
!                2 : modified GRIB-0 encoding
!                3 : standard GRIB-1 encoding
!                4 : hollow fields GRIB-1 encoding
!              181 : GRIB-2 encoding
!     NFPOSTS    : array containing postprocessing steps
!     NFPOSTSMIN : array containing postprocessing steps in minutes for sub-hour outputs
!     NFRFPOS    : frequency of post-processing events


!*    Scientific and technical variables

! === TECHNICAL VARIABLES ===

!     CFP3DF : names 3D dynamics fields to compute
!     CFP2DF : names 2D dynamics fields to compute
!     CFPPHY : names of physical fields to be post-processed
!     CFPCFU : names of cumulated fluxes fields to be post-processed
!     CFPXFU : names of instantaneous fluxes fields to be post-proc
!     Remark: C1FPXXX=CFPXXX(1)//CFPXXX(2)//...//CFPXXX(JPOSXXX)
!      for XXX=DOM,3DF,2DF,PHY,CFU,XFU.

!     Variables used for ECMWF

!     MFP3DFS: Gribcodes of  3D dynamics fields to compute on model levels
!     MFP3DFH: Gribcodes of  3D dynamics fields to compute on H levels
!     MFP3DFT: Gribcodes of  3D dynamics fields to compute on THETA levels
!     MFP3DFV: Gribcodes of  3D dynamics fields to compute on PV levels
!     MFP3DFP: Gribcodes of  3D dynamics fields to compute on P levels
!     MFP3DFI: Gribcodes of  3D dynamics fields to compute on T levels
!     MFP3DFF: Gribcodes of  3D dynamics fields to compute on Flight levels
!     MFP2DF : Gribcodes of  2D dynamics fields to compute
!     MFPPHY : Gribcodes of  physical fields to processed

!     RFP3P  : post-processing pressure levels
!     RFP3H  : post-processing height (above orography) levels
!     RFP3TH : post-processing potential temperature levels
!     RFP3PV : post-processing potential vorticity levels
!     RFP3I  : post-processing temperature levels
!     RFP3F  : post-processing height (above sea) levels
!     NRFP3S : post-processing eta levels (CONF. 927 only)

!     NFP3DFS: useful dimension of MFP3DFS
!     NFP3DFH: useful dimension of MFP3DFH
!     NFP3DFT: useful dimension of MFP3DFT
!     NFP3DFV: useful dimension of MFP3DFV
!     NFP3DFP: useful dimension of MFP3DFP
!     NFP3DFI: useful dimension of MFP3DFI
!     NFP3DFF: useful dimension of MFP3DFF
!     NFP2DF : useful dimension of CFP2DF
!     NFPPHY : useful dimension of CFPPHY

!     LFPLOSP : .TRUE. = Fill Ps array with Log(Ps)
!     LFPRH100 : .TRUE. to convert relative humidity in percent

! === SCIENTIFIC VARIABLES ===

!     NFIT*  : keys for spectral fit of post-processed fields
!              =0 if no spectral fit on this post-processing surface
!              =1 if spectral fit in the model geometry
!              =2 if spectral fit in the target geometries
!     NFITP  : P levels
!     NFITT  : THETA levels
!     NFITV  : PV levels
!     NFITI  : T levels
!     NFITF  : Flight levels
!     NFITH  : height levels
!     NFITS  : eta levels

!     NSPFIL*  : kind of generic spectral filter on post-processing levels
!              =1 : no filter
!              =2 : classic (gaussian filter)
!              =3 : for stretched geometry (applied one the spectrum of homogenous resolution
!     NSPFILP  : P levels
!     NSPFILT  : THETA levels
!     NSPFILV  : PV levels
!     NSPFILI  : T levels
!     NSPFILF  : Flight levels
!     NSPFILH  : height levels
!     NSPFILS  : eta levels
!     NSPFIL2  : 2D fields

!     LFITBETWEEN : Spectral fit between horizontal and lagged vertical interpolation
!                   which was implicitly done in the old 927 configuration.
!     LSATURCAP : Cap humididy to saturation between horizontal and lagged vertical interpolation
!                 which was done for ECMWF in the old 927 configuration.

!     NFPCLI : usage level for climatology
!              =0 no climatology
!              =1 orography and land-sea mask of output only
!              =2 all available climatological fields of the current month
!              =3 shifting mean from the climatological fields of the current
!                 month to the ones of the closest month

!     LCLIMALBEDOS : to read albedos in clim file (for compatibility with old clim files without them)
!     LCLIMAEROSOL : to read aerosols in clim file (for compatibility with old clim files without them)
!     LCLIMOZONE   : to read ozone coefficients in clim file (for compatibility with old clim files without them)
!     NFPINDYN : type of interpolations for dynamical fields
!     NFPINPHY : type of interpolations for  physical fields
!                4 for bilinear, 12 for quadratic

!     LFPQ   : =.TRUE. if specific humidity is interpolated
!            : =.FALSE. if relative humidity is interpolated

!     RFPCORR: Critical orography difference for correcting surface
!              temperature through standart profile.
!     RFPCSAB: Critical sand percentage difference for computing relative
!              soil moisture in ISBA
!     RFPCD2 : Critical soil depth difference for computing relative
!              soil moisture in ISBA
!     RFPVCAP: Minimum pressure of model level to provide an equatorial
!              cap in the computation of variables on constant PV surfaces
!     LFPISOPV : .T. => new diagnostic for computing an iso-PV level
!     NITERPV : Nb of vertical iter (1, 2 or 3) used in iso-PV level computing
!     NFPLAKE : To overwrite created lakes or islands by specific data :
!               0 => do not overwrite
!              -1 => overwrite with rhoughly interpolated data
!              +1 => overwrite with climatology
!     NFPCAPE : Kind of computation for CAPE & CIN :
!               1 => from bottom model layer
!               2 => from the most unstable layer
!               3 => from mto standart height (2 meters) as recomputed values
!               4 => from mto standart height (2 meters) out of fluxes
!                    (used for analysis)
!     LFPCAPEX : if true XFU fields used for CAPE&CIN computation (with NFPCAPE)
!     NFPSURFEX : Subcontract surface fields to SURFEX
!               0 => no subcontract
!               1 => transform native arp/ald surface fields to surfex fields
!                    and write out by surfex
!     LPUTZS    : if true "atmospheric" orography is imposed to surfex
!     NFPMASK   : number of masks for the interpolation of surface fields
!                 0 => no mask
!                 1 => land-sea mask
!                 2 => land mask, sea mask
!     RENTRA     : entrainement coefficient, if 0, (default) no entrainement is performed,
!                  allows to have an even more complicated choice of options for cape computation
!     FPRHMIN,FPRHMAX: min and max allowed values for relative humidity in pp.
!     LFPML_STD  : multi-linear interpolations if T for interpolations using weights WSTD...
!     LFPML_LAN  : multi-linear interpolations if T for interpolations using weights WLAN...
!     LFPML_SEA  : multi-linear interpolations if T for interpolations using weights WSEA...
!     RFPCR_STD  : coefficients for the ratio of mesh-sizes used if LFPML_STD=T
!     RFPCR_LAN  : coefficients for the ratio of mesh-sizes used if LFPML_LAN=T
!     RFPCR_SEA  : coefficients for the ratio of mesh-sizes used if LFPML_SEA=T
!     LCRITSNOWTEMP: Use a critical surface temperature (0 C) to remove existing snow in input field
!     LISOT_ABOVEG: Set K[TB][NNN]ISOT_ALTIT to missing value when below the
!                   orography
!
!
! Options for Boyd biperiodization (daand, 02/2012)
!     NFPBOYD    : periodization with Boyd's windowing method
!                  0 = splines
!                  1 = Boyd's windowing method
!     RFPBSCAL   : scalar parameter for Boyd periodization ("L" in Boyd's paper)
!
!     NFPSFXINT  : 0 = 4 point interpolation for PREP fields, 1 = boxed average


! Options for Interoperability IFS+TESSEL -> ARPEGE/AROME+ISBA
!     TSRESERV1  : Icing temperature for surface soil frost
!     TSRESERV2  : Icing temperature for deep    soil frost
!     TDELTA1    : Threshold of icing temperature for surface soil frost
!     TDELTA2    : Threshold of icing temperature for deep    soil frost
!     NFPSWI     : Soil water or frost interpolation method for interoperability (from TESSEL to ISBA)
!                  0 = Interpolate relative soil water content
!                      Interpolate frost from soil water, according to the following formula :
!                      p = 3*Teta**2 + 2*Teta**3
!                      where p denotes the proportion of frost inside the whole soil wetness
!                      Teta  = (Tlue-T0)/(2*delta) - 0.5
!                      Tlue  = temperature of the layer (surface or deep soil) 
!                      T0    = Icing temperature of the layer (resp. TSRESERV1 or TSRESERV2)
!                      delta = threshold of T0 (resp. TDELTA1 or TDELTA2)
!                  1 = Interpolate Soil Wetness Index then scale following the so-called "OI trick"
!                  2 = Interpolate Soil Wetness Index then scale with Leaf Area Index
!                      Interpolate frost from Soil Wetness Index, according to the following formula :
!                      p = (1-sin(Teta))/2
!                      where p denotes the proportion of frost inside the whole soil wetness
!                      Teta     = Pi * (Tlue-T0)/(2*delta)
!                      Tlue  = temperature of the layer (surface or deep soil) 
!                      T0    = Icing temperature of the layer (resp. TSRESERV1 or TSRESERV2)
!                      delta = threshold of T0 (resp. TDELTA1 or TDELTA2)
!
! Options for ocean post processing
!     LOCEDELAY  : Output ocean fields after WAM/NEMO has been called.
!

! Options for Interoperability SURFEX -> ISBA
!    RWPITPN : minimum deep soil temperature for melting of frozen water
!    RWPITPX : maximum deep soil temperature for melting of frozen water
!    RSNSTPN : minimum deep soil temperature for melting of snow
!    RSNSTPX : minimum deep soil temperature for melting of snow
!    RSNSMOD : Reference snow depth for normalization


LOGICAL :: LTRACEFP = .TRUE.
LOGICAL :: LALLOFP = .FALSE.
LOGICAL :: LFPMOIS = .FALSE.
LOGICAL :: L_READ_MODEL_DATE = .FALSE.
LOGICAL :: LFPPACKING = .FALSE.
LOGICAL :: LOCEDELAY = .FALSE.
LOGICAL :: LWIDER_DOM = .FALSE.
INTEGER(KIND=JPIM) :: NFPCHKDAT = 1
INTEGER(KIND=JPIM) :: NFPSFXWRT = 0
REAL(KIND=JPRB) :: RWIDOM


TYPE TNAMFPL

CHARACTER (LEN = 12) ::  CFP3DF(JPOS3DF) = ' '
CHARACTER (LEN = 16) ::  CFPPHY(JPOSSGP) = ' '

INTEGER(KIND=JPIM) :: MFP3DFS(JPOS3DF) = 0
INTEGER(KIND=JPIM) :: MFP3DFP(JPOS3DF) = 0
INTEGER(KIND=JPIM) :: MFP3DFH(JPOS3DF) = 0
INTEGER(KIND=JPIM) :: MFP3DFT(JPOS3DF) = 0
INTEGER(KIND=JPIM) :: MFP3DFV(JPOS3DF) = 0
INTEGER(KIND=JPIM) :: MFP3DFI(JPOS3DF) = 0
INTEGER(KIND=JPIM) :: MFP3DFF(JPOS3DF) = 0
INTEGER(KIND=JPIM) :: MFP2DF(JPOS2DF) = 0

REAL(KIND=JPRB) :: RFP3P(JPOS3P)   = -9._JPRB
REAL(KIND=JPRB) :: RFP3H(JPOS3H)   = -9._JPRB
REAL(KIND=JPRB) :: RFP3TH(JPOS3TH) = -9._JPRB
REAL(KIND=JPRB) :: RFP3PV(JPOS3PV) = 9999000000._JPRB
REAL(KIND=JPRB) :: RFP3I(JPOS3I)   = -9999._JPRB
REAL(KIND=JPRB) :: RFP3F(JPOS3F)   = -9._JPRB
INTEGER(KIND=JPIM) :: NRFP3S(JPOS3S) = -9

INTEGER(KIND=JPIM) :: NFP3DFS = 0
INTEGER(KIND=JPIM) :: NFP3DFP = 0
INTEGER(KIND=JPIM) :: NFP3DFH = 0
INTEGER(KIND=JPIM) :: NFP3DFT = 0
INTEGER(KIND=JPIM) :: NFP3DFV = 0
INTEGER(KIND=JPIM) :: NFP3DFI = 0
INTEGER(KIND=JPIM) :: NFP3DFF = 0
INTEGER(KIND=JPIM) :: NFP2DF = 0

INTEGER(KIND=JPIM) :: MFPPHY(JPOSSGP) = -9999
INTEGER(KIND=JPIM) :: NFPPHY = 0
CHARACTER (LEN = 16) ::  CFPCFU(JPOSCFU) = ' '
CHARACTER (LEN = 16) ::  CFPXFU(JPOSXFU) = ' '
CHARACTER (LEN = 16) ::  CFP2DF(JPOS2DF) = ' '

END TYPE TNAMFPL


! parameters driving the horizontal interpolations
! Special structure (to be re-worked)

TYPE TNAMFPINT

LOGICAL :: LFPML_STD = .FALSE.
LOGICAL :: LFPML_LAN = .FALSE.
LOGICAL :: LFPML_SEA = .FALSE.
INTEGER(KIND=JPIM) :: NFPINDYN = 12
INTEGER(KIND=JPIM) :: NFPINPHY = 12
INTEGER(KIND=JPIM) :: NFPSFXINT = 0
REAL(KIND=JPRB) :: RFPCR_STD(2) = (/1._JPRB, 0._JPRB/)
REAL(KIND=JPRB) :: RFPCR_LAN(2) = (/1._JPRB, 0._JPRB/)
REAL(KIND=JPRB) :: RFPCR_SEA(2) = (/1._JPRB, 0._JPRB/)

END TYPE TNAMFPINT

TYPE TNAMFPSCI

REAL(KIND=JPRB)    :: RFPVCAP
REAL(KIND=JPRB)    :: FPRHMAX
LOGICAL            :: LCRITSNOWTEMP
LOGICAL            :: LISOT_ABOVEG
LOGICAL            :: LFPRH100

INTEGER(KIND=JPIM) :: NFPSURFEX = 0
INTEGER(KIND=JPIM) :: NFPBOYD = 0
REAL(KIND=JPRB)    :: RFPBSCAL = 3._JPRB
LOGICAL            :: LFPCAPEX = .FALSE.
INTEGER(KIND=JPIM) :: NFPCLI = 0
LOGICAL            :: LCLIMALBEDOS = .TRUE.
LOGICAL            :: LCLIMAEROSOL = .TRUE.
LOGICAL            :: LCLIMOZONE = .TRUE.
INTEGER(KIND=JPIM) :: NFPLAKE = 0
INTEGER(KIND=JPIM) :: NFPMASK = 1
INTEGER(KIND=JPIM) :: NFPSWI = 0
INTEGER(KIND=JPIM) :: NFPCAPE = 2
LOGICAL            :: LFPISOPV = .TRUE.
INTEGER(KIND=JPIM) :: NITERPV = 3
REAL(KIND=JPRB)    :: RENTRA = 0.0_JPRB
REAL(KIND=JPRB)    :: FPRHMIN = 0.0_JPRB
LOGICAL            :: LSATURCAP = .FALSE.
LOGICAL            :: LFITBETWEEN = .FALSE.
LOGICAL            :: LFPQ = .FALSE.
LOGICAL            :: LFPLOSP = .FALSE.
LOGICAL            :: LPUTZS = .TRUE.
LOGICAL            :: LFPCLSTOGMV = .FALSE.
REAL(KIND=JPRB)    :: TSRESERV1 = 272.15_JPRB ! as used in operation (2016)
REAL(KIND=JPRB)    :: TSRESERV2 = 269.15_JPRB ! as used in operation (2016)
REAL(KIND=JPRB)    :: TDELTA1 = 3._JPRB
REAL(KIND=JPRB)    :: TDELTA2 = 7._JPRB
REAL(KIND=JPRB)    :: RWPITPN = 273._JPRB
REAL(KIND=JPRB)    :: RWPITPX = 273._JPRB
REAL(KIND=JPRB)    :: RSNSTPN = 272._JPRB
REAL(KIND=JPRB)    :: RSNSTPX = 272._JPRB
REAL(KIND=JPRB)    :: RSNSMOD = 80._JPRB
REAL(KIND=JPRB)    :: RFPCORR = 2943._JPRB ! roughly 300.*g
REAL(KIND=JPRB)    :: RFPCSAB = 0.01_JPRB
REAL(KIND=JPRB)    :: RFPCD2 = 0.001_JPRB

REAL(KIND=JPRB)    :: RSTRMMH=6000._JPRB
REAL(KIND=JPRB)    :: RSRHH=3000._JPRB

INTEGER(KIND=JPIM) :: NFITP
INTEGER(KIND=JPIM) :: NFITT
INTEGER(KIND=JPIM) :: NFITV
INTEGER(KIND=JPIM) :: NFITI
INTEGER(KIND=JPIM) :: NFITF
INTEGER(KIND=JPIM) :: NFITH
INTEGER(KIND=JPIM) :: NFITS
INTEGER(KIND=JPIM) :: NSPFILP
INTEGER(KIND=JPIM) :: NSPFILT
INTEGER(KIND=JPIM) :: NSPFILV
INTEGER(KIND=JPIM) :: NSPFILI
INTEGER(KIND=JPIM) :: NSPFILF
INTEGER(KIND=JPIM) :: NSPFILH
INTEGER(KIND=JPIM) :: NSPFILS
INTEGER(KIND=JPIM) :: NSPFIL2

END TYPE TNAMFPSCI

TYPE TNAMFPOBJ

CHARACTER(LEN=5) :: CFPFMT
CHARACTER(LEN=JPOSLEN) :: CFPDOM(JPOSDOM)
CHARACTER(LEN=JPOSDIR) :: CFPDIR = ' '
CHARACTER(LEN=180) :: CFPMONIPATH_IN = '.'
CHARACTER(LEN=180) :: CFPMONIPATH_OUT = '.'
INTEGER(KIND=JPIM) :: NFPGRIB
INTEGER(KIND=JPIM) :: NFPOSTS(0:JPNPST) = 0
INTEGER(KIND=JPIM) :: NFPOSTSMIN(0:JPNPST) = 0
INTEGER(KIND=JPIM) :: NFRFPOS = 1

END TYPE TNAMFPOBJ
!     ------------------------------------------------------------------
END MODULE YOMFPC

