! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE RADLSWR&
 &(YDDIMV, YDMODEL,KIDIA, KFDIA , KLON , KLEV , KAERO ,&
 &PRII0,&
 &PAER , PAERO , PALBD, PALBP, PAPH  , PAP  ,&
 &PCCNL, PCCNO , PGELAM,PGEMU,&
 &PCO2 , PCH4  , PN2O , PNO2 , PC11  , PC12 , PC22 , PCL4 ,&
 &PCLFR, PDP   , PEMIS, PEMIW, PLSM  , PMU0 , POZON,&
 &PQ   , PQIWP , PQSWP, PQLWP, PQRWP ,&
 &PTH  , PT    , PTS  ,&
 &PEMIT, PFCT  , PFLT , PFCS , PFLS  , PFRSOD, PFRTED, PFRSODC, PFRTEDC ,&
 &PSUDU, PUVDF , PPARF, PPARCF, PTINCF,&
 &PFDIR, PFDIF,PCDIR , PLWDERIVATIVE, PSWDIFFUSEBAND, &
 &PSWDIRECTBAND , PPERT)

!**** *RADLSWR* - INTERFACE TO ECMWF McICA LW AND SW RADIATION SCHEMES 

!     PURPOSE.
!     --------
!           CONTROLS RADIATION COMPUTATIONS

!**   INTERFACE.
!     ----------

!        EXPLICIT ARGUMENTS :
!        --------------------
! PAER   : (KLON,6,KLEV)     ; OPTICAL THICKNESS OF THE AEROSOLS
! PAERO  : (KLON,KLEV,NAERO) ; PROGNOSTIC AEROSOL MIXING RATIO
! PALBD  : (KLON,NSW)        ; SURF. SW ALBEDO FOR DIFFUSE RADIATION
! PALBP  : (KLON,NSW)        ; SURF. SW ALBEDO FOR PARALLEL RADIATION
! PAPH   : (KLON,KLEV+1)     ; HALF LEVEL PRESSURE
! PAP    : (KLON,KLEV)       ; FULL LEVEL PRESSURE
! PCCNL  : (KLON)            ; CCN CONCENTRATION OVER LAND
! PCCNO  : (KLON)            ; CCN CONCENTRATION OVER OCEAN
! PCO2   : (KLON,KLEV)       ; CONCENTRATION IN CO2 (KG/KG)
! PCH4   : (KLON,KLEV)       ; CONCENTRATION IN CH4 (KG/KG)
! PN2O   : (KLON,KLEV)       ; CONCENTRATION IN N2O (KG/KG)
! PNO2   : (KLON,KLEV)       ; CONCENTRATION IN NO2 (KG/KG)
! PC11   : (KLON,KLEV)       ; CONCENTRATION IN CFC11 (KG/KG)
! PC12   : (KLON,KLEV)       ; CONCENTRATION IN CFC12 (KG/KG)
! PC22   : (KLON,KLEV)       ; CONCENTRATION IN CFC22 (KG/KG)
! PCL4   : (KLON,KLEV)       ; CONCENTRATION IN CCL4  (KG/KG)
! PCLFR  : (KLON,KLEV)       ; CLOUD FRACTIONAL COVER
! PDP    : (KLON,KLEV)       ; LAYER PRESSURE THICKNESS
! PEMIS  : (KLON)            ; SURFACE LW EMISSIVITY
! PEMIW  : (KLON)            ; SURFACE LW WINDOW EMISSIVITY
! PGELAM : (KLON)            ; LONGITUDE (RADIANS)
! PGEMU  : (KLON)            ; SINE OF LATITUDE
! PLSM   : (KLON)            ; LAND-SEA MASK
! PMU0   : (KLON)            ; SOLAR ANGLE
! POZON  : (KLON,KLEV)       ; OZONE AMOUNT in LAYER (KG/KG*PA)
! PQ     : (KLON,KLEV)       ; SPECIFIC HUMIDITY KG/KG
! PQIWP  : (KLON,KLEV)       ; ICE (Small particles) KG/KG
! PQSWP  : (KLON,KLEV)       ; SNOW (Large particle) KG/KG
! PQLWP  : (KLON,KLEV)       ; LIQUID WATER KG/KG
! PQRWP  : (KLON,KLEV)       ; RAIN WATER KG/KG
! PTH    : (KLON,KLEV+1)     ; HALF LEVEL TEMPERATURE
! PT     : (KLON,KLEV)       ; FULL LEVEL TEMPERATURE
! PTS    : (KLON)            ; SURFACE TEMPERATURE
! PPERT  : (KLON,YDSPP_CONFIG%SM%NRFTOTAL_RADGRID); random patterns from SPP
!     ==== OUTPUTS ===
! PEMIT  : (KLON)            ; SURFACE TOTAL LONGWAVE EMISSIVITY (weighted average across the spectrum)
! PFCT   : (KLON,KLEV+1)     ; CLEAR-SKY LW NET FLUXES
! PFLT   : (KLON,KLEV+1)     ; TOTAL LW NET FLUXES
! PFCS   : (KLON,KLEV+1)     ; CLEAR-SKY SW NET FLUXES
! PFLS   : (KLON,KLEV+1)     ; TOTAL SW NET FLUXES
! PFRSOD : (KLON)            ; TOTAL-SKY SURFACE SW DOWNWARD FLUX
! PFRTED : (KLON)            ; TOTAL-SKY SURFACE LW DOWNWARD FLUX
! PFRSODC: (KLON)            ; CLEAR-SKY SURFACE SW DOWNWARD FLUX 
! PFRTEDC: (KLON)            ; CLEAR-SKY SURFACE LW DOWNWARD FLUX 
! PSUDU  : (KLON)            ; SOLAR RADIANCE IN SUN'S DIRECTION
! PFDIR  : (KLON)            ; TOTAL SKY DIRECT SW DOWNWARD FLUX
! PFDIF  : (KLON)            ; TOTAL SKY DIFUSE SW DOWNWARD FLUX
! PCDIR  : (KLON)            ; CLEAR-SKY DIRECT SW DOWNWARD FLUX
! PLwDerivative:(KLON,KLEV+1); see below for full definition
! PSwDiffuseBand: (KLON,NSW) ; see below for full definition
! PSwDirectBand:  (KLON,NSW) ; see below for full definition
!        IMPLICIT ARGUMENTS :   NONE
!        --------------------

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS

!     AUTHORS.
!     --------
!      J.-J. MORCRETTE         *ECMWF*
!      ORIGINAL : 88-02-04

!     MODIFICATIONS.
!     --------------
!      JJMorcrette : 010112 Sun-Rikus ice particle Diameter
!      JJMorcrette : 010301 cleaning liq/ice cloud optical properties
!      JJMorcrette : 011005 CCN --> Re liquid water clouds
!      JJMorcrette : 011108 Safety checks
!      JJMorcrette : 011108 Safety checks
!      DJSalmond   : 020211 Check before R-To-R
!      JJMorcrette : 030430 Interface to SRTM-224gp
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      JJMorcrette 20060721 PP of clear-sky and TOA incident solar radiation
!      JJMorcrette : 060110 McICA approach to RRTM and SRTM-112g
!      JJMorcrette : 20071015 3D fields for CO2, CH4, N2O and NO2
!      JJMorcrette 20070321 Prognostic aerosols in radiation
!      R.Forbes    20071015 Add RDECORR_CF,RDECORR_CW
!      JJMorcrette 20080422 Prognostic aerosols in radiation
!      JJMorcrette 20090804 Change in the partition IW/LW and ice cloud De
!      JJMorcrette 20090805 Decorrelation=f(latitude) 
!      JJMorcrette 20091201 Total and clear-sky direct SW radiation flux at surface 
!   JJMorcrette 20100322 Prognostic aerosols in radiation with consistent optical properties
!      M Ahlgrimm 31 Oct 2011 Surface Downward clear-sky LW and SW fluxes
!      PBechtold+NSemane     09-Jul-2012 gravity
!      K. Yessad (July 2014): Move some variables.
!      A Bozzo     20130509 bug in LW emissivity for NICEOPT=3 - not fixed yet
!      A Bozzo     20130509 bugfix in RFUETC (wrong initialized in suclopn.F90 for bands 14,15,16)
!      A Bozzo     201305   added boundaries for g-LW (NICEOPT=3)
!      JJMorcrette+ABozzo 20130805/201411 MACC-derived aerosol climatology
!      RJHogan  20 May 2014 Added partial derivatives and total-sky surface LW down outputs
!      RJHogan  29 June 2014 Pass through PSw*Band
!      RJHogan  12 February 2016 Call cloud_overlap_decorr_len
!      RJHogan  21 June 2016 Correct ice single-scattering albedo
!      M. Leutbecher & S.-J. Lock (Jan 2016) Introduced SPP scheme (LSPP)
!      RJHogan  30 June 2016 Cloud fractional standard deviation configurable
!      RJHogan  15 Sept 2016 Added Dr Hook profiling for aerosol and cloud optics
!      RJHogan  29 Sept 2016 SPP applied to ice effective radius when NRADIP==3 (default)
!      Y. Bouteloup 21 October 2016 : Add new output flux PFDIF as PFDIR : Solar Downward diffuse
!      M Leutbecher Oct 2020 SPP abstraction
!-----------------------------------------------------------------------

USE TYPE_MODEL , ONLY : MODEL
USE YOMDIMV    , ONLY : TDIMV
USE PARKIND1   , ONLY : JPIM, JPRB
USE YOMHOOK    , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCT3     , ONLY : NSTEP
USE YOMCST     , ONLY : RG       ,RD       ,RTT      ,RPI
USE YOELW      , ONLY : NTRA     ,NUA
USE YOESRTCOP  , ONLY : &
 &                      RSYFWA   ,RSYFWB   ,RSYFWC   ,RSYFWD   ,RSYFWE   ,RSYFWF&
 &                    , RSECIA   ,RSECIB   ,RSECIC   ,RSECID   ,RSECIE   ,RSECIF&
 &                    , RSASWA   ,RSASWB   ,RSASWC   ,RSASWD   ,RSASWE   ,RSASWF&
 &                    , RSFUA0   ,RSFUA1   ,RSFUB0   ,RSFUB1   ,RSFUB2   ,RSFUB3&
 &                    , RSFUC0   ,RSFUC1   ,RSFUC2   ,RSFUC3&
 &                    , RSFLA0   ,RSFLA1   ,RSFLB0   ,RSFLB1   ,RSFLB2   ,RSFLB3&
 &                    , RSFLC0   ,RSFLC1   ,RSFLC2   ,RSFLC3   ,RSFLD0   ,RSFLD1&
 &                    , RSFLD2   ,RSFLD3  
USE YOERDU     , ONLY : NUAER  ,NTRAER   ,REPLOG   ,REPSC    ,REPSCA  ,REPSCW   ,DIFF
USE YOETHF     , ONLY : RTICE
USE YOERRTFTR  , ONLY : NGB
USE YOERRTM    , ONLY : JPGLW 
USE YOESRTM    , ONLY : JPGSW, NGBSW
USE YOESRTWN   , ONLY : NMPSRTM
USE YOEAEROP   , ONLY : ALF_SU
USE YOMDYNCORE , ONLY : RPLRG
USE YOMLUN     , ONLY : NULOUT
USE SPP_GEN_MOD, ONLY:SPP_PERT

!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE(TDIMV)       ,INTENT(IN)    :: YDDIMV
TYPE(MODEL)       ,INTENT(INOUT) :: YDMODEL
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KAERO
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRII0
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAER(KLON,6,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAERO(KLON,KLEV,KAERO)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PALBD(KLON,YDMODEL%YRML_PHY_RAD%YRERAD%NSW) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PALBP(KLON,YDMODEL%YRML_PHY_RAD%YRERAD%NSW) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPH(KLON,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAP(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCCNL(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCCNO(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCO2(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCH4(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PN2O(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PNO2(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PC11(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PC12(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PC22(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCL4(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCLFR(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDP(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEMIS(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEMIW(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGELAM(KLON), PGEMU(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLSM(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMU0(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: POZON(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQ(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQIWP(KLON,KLEV) ! Ice 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQSWP(KLON,KLEV) ! Snow
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQLWP(KLON,KLEV) ! Liquid
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQRWP(KLON,KLEV) ! Rain
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTH(KLON,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTS(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PEMIT(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFCT(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFLT(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFCS(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFLS(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFRSOD(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFRTED(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFRSODC(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFRTEDC(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSUDU(KLON)  
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PUVDF(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PPARF(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PPARCF(KLON), PTINCF(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFDIR(KLON) , PFDIF(KLON), PCDIR(KLON)

! Partial derivative of total-sky longwave upward flux at each level
! with respect to upward flux at surface, used to correct heating
! rates at gridpoints/timesteps between calls to the full radiation
! scheme.  Note that this version uses the convention of level index
! increasing downwards, unlike the local variable ZLwDerivative
! that is returned from the LW radiation scheme.
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PLWDERIVATIVE(KLON,KLEV+1) 

! Surface diffuse and direct downwelling shortwave flux in each
! shortwave albedo band, used in RADINTG to update the surface fluxes
! accounting for high-resolution albedo information
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSWDIFFUSEBAND(KLON,YDMODEL%YRML_PHY_RAD%YRERAD%NSW)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSWDIRECTBAND(KLON,YDMODEL%YRML_PHY_RAD%YRERAD%NSW)
REAL(KIND=JPRB)   ,INTENT(IN),OPTIONAL :: PPERT(KLON, YDMODEL%YRML_GCONF%YRSPP_CONFIG%SM%NRFTOTAL_RADGRID)

!     -----------------------------------------------------------------

REAL(KIND=JPRB) ::&
 & ZTAUCLD(KLON,KLEV,16), ZTCLEAR(KLON)  
REAL(KIND=JPRB) :: ZAER(KLON,6,KLEV)&
 & , ZCLDLW(KLON,KLEV)  , ZCLDSW(KLON,KLEV)&
 & , ZDT0(KLON)&
 & , ZEMIS(KLON)        , ZEMIW(KLON)&
 & , ZFDIR(KLON,KLEV+1) , ZCDIR(KLON,KLEV+1), ZFDIF(KLON,KLEV+1), ZCDIF(KLON,KLEV+1)&
 & , ZFLUX (KLON,2,KLEV+1)                  , ZFLUC(KLON,2,KLEV+1)&
 & , ZFIWP(KLON,KLEV)   , ZFSWP(KLON,KLEV) , ZFLWP(KLON,KLEV), ZFRWP(KLON,KLEV)&
 & , ZIWC(KLON,KLEV)    , ZSWC(KLON,KLEV)  , ZLWC(KLON,KLEV) , ZRWC(KLON,KLEV)&
 & , ZMU0(KLON)         , ZOZ(KLON,KLEV)   , ZOZN(KLON,KLEV)&
 & , ZPMB(KLON,KLEV+1)  , ZPSOL(KLON)&
 & , ZTAVE (KLON,KLEV)  , ZTL(KLON,KLEV+1)&
 & , ZVIEW(KLON)  
REAL(KIND=JPRB) ::&
 & ZFCDWN(KLON,KLEV+1), ZFCUP(KLON,KLEV+1)&
 & , ZFSDWN(KLON,KLEV+1), ZFSUP(KLON,KLEV+1)&
 & , ZFSUPN(KLON)       , ZFSUPV(KLON)&
 & , ZFCUPN(KLON)       , ZFCUPV(KLON)&
 & , ZFSDNN(KLON)       , ZFSDNV(KLON)&
 & , ZFCDNN(KLON)       , ZFCDNV(KLON)  
REAL(KIND=JPRB) ::&
 & ZDESR(KLON), ZRADIP(KLON)       , ZRADLP(KLON)    ,&
 !cc           , ZRADRD(KLON)
 & ZRAINT(KLON)       , ZEMIT(KLON)     
REAL(KIND=JPRB) :: ZSUDU(KLON), ZFUVF(KLON), ZFUVC(KLON), ZPARF(KLON), ZPARCF(KLON)

INTEGER(KIND=JPIM) :: ICLIM
! Partial derivative of total-sky longwave upward flux at each level
! with respect to upward flux at surface, used to correct heating
! rates at gridpoints/timesteps between calls to the full radiation
! scheme. Note that this is the local version that uses the radiation
! convention of level index increasing upwards; PLwDerivative is
! the version passed back to the caller of this function that uses the
! model convention of level index increasing downwards.
REAL(KIND=JPRB) :: ZLWDERIVATIVE(KLON,KLEV+1) 

INTEGER(KIND=JPIM) :: IFLAG, IK, ILW16, ISW, ISW14
INTEGER(KIND=JPIM) :: JK, JKL, JKLP1, JKP1, JL, JRTM, JSW
INTEGER(KIND=JPIM) :: JAER

REAL(KIND=JPRB) :: ZASYMX, ZDIFFD, ZGI, ZGL, ZGR 
REAL(KIND=JPRB) :: ZIWGKG, ZSWGKG, ZLWGKG, ZRWGKG,&
 & ZMULTI, ZMULTL, ZOI, ZOL,&
 & ZOMGMX, ZOR, ZRMUZ, ZTAUD, ZTAUMX, ZTEMPC,&
 & ZTOI, ZTOL, ZTOR, ZDPOG, ZPODT  
REAL(KIND=JPRB) :: ZALND, ZASEA, ZDEN, ZNTOT, ZNUM, Z1RADI, Z1RADL,&
 & ZBETAI, ZOMGI, ZOMGP, ZFDEL, ZTCELS, ZFSR, ZAIWC, ZBIWC, ZEXPN, ZRELVNT,&
 & ZDEFRE, ZREFDE, ZEXTCF, Z1MOMG, ZRSAIA, ZRSAID, ZRSAIE, ZRSAIF, ZRSAIG, ZRSALD, ZRSAIO

REAL(KIND=JPRB) :: ZSALBD(KLON,14)     , ZSALBP(KLON,14)
REAL(KIND=JPRB) :: ZSTAUC(KLON,14,KLEV), ZSASYC(KLON,14,KLEV), ZSOMGC(KLON,14,KLEV)
REAL(KIND=JPRB) :: ZFSUC(KLON,2,KLEV+1), ZFSUX(KLON,2,KLEV+1)

!!!!!INTEGER(KIND=JPIM), PARAMETER :: ICOLSLW=JPGLW, ICOLSSW=JPGSW
INTEGER(KIND=JPIM) :: ICOLS, IJK, IPRINT(KLON,KLEV), JCOLS, ICLW, ICSW, JCLW, JCSW

REAL(KIND=JPRB)    :: ZCLFR(KLON,KLEV) , ZQIWP(KLON,KLEV) , ZQSWP(KLON,KLEV) 
REAL(KIND=JPRB)    :: ZQLWP(KLON,KLEV) , ZQRWP(KLON,KLEV) , ZETA(KLON,KLEV)
REAL(KIND=JPRB)    :: ZGLAT(KLON)      , ZGLON(KLON)
REAL(KIND=JPRB)    :: ZLC_CF(KLON,KLEV), ZLC_CW(KLON,KLEV), ZSIGQCW(KLON,KLEV)
REAL(KIND=JPRB)    :: ZQI_SL(KLON,KLEV,JPGLW), ZQL_SL(KLON,KLEV,JPGLW)
REAL(KIND=JPRB)    :: ZQI_SS(KLON,KLEV,JPGSW), ZQL_SS(KLON,KLEV,JPGSW)
REAL(KIND=JPRB)    :: ZQS_SL(KLON,KLEV,JPGLW), ZQS_SS(KLON,KLEV,JPGSW)
REAL(KIND=JPRB)    :: ZQILW(KLON,JPGLW,KLEV) , ZQLLW(KLON,JPGLW,KLEV)
REAL(KIND=JPRB)    :: ZQISW(KLON,JPGLW,KLEV) , ZQLSW(KLON,JPGLW,KLEV)
REAL(KIND=JPRB)    :: ZQSLW(KLON,JPGLW,KLEV) , ZQSSW(KLON,JPGLW,KLEV)
INTEGER(KIND=JPIM) :: ICLDLW(KLON), ICLDSW(KLON)

REAL(KIND=JPRB)    :: ZRADLI(KLON,KLEV), ZRADIC(KLON,KLEV), ZDESIC(KLON,KLEV)

REAL(KIND=JPRB)    :: ZFCLW(KLON,JPGLW,KLEV), ZTOLW(KLON,KLEV,JPGLW)
REAL(KIND=JPRB)    :: ZFCSW(KLON,JPGSW,KLEV) , ZTOSW(KLON,JPGSW,KLEV)
REAL(KIND=JPRB)    :: ZASSW(KLON,JPGSW,KLEV) , ZOMSW(KLON,JPGSW,KLEV)
REAL(KIND=JPRB)    :: ZSWFC(KLON,JPGSW,KLEV) , ZLWFC(KLON,JPGLW,KLEV) 
INTEGER(KIND=JPIM) :: IFCLW(KLON,JPGLW,KLEV) , IFCSW(KLON,JPGSW,KLEV)

REAL(KIND=JPRB)    :: ZRG, ZTEMP, ZRPI

! Effective radius for water cloud
REAL(KIND=JPRB)    :: ZK  ! constant relating rvol to reff for cloud droplets
REAL(KIND=JPRB)    :: ZKL ! constant relating rvol to reff for drizzle/rain droplets
REAL(KIND=JPRB)    :: ZKRATIO   ! ratio of ZKL to ZK
REAL(KIND=JPRB)    :: ZRLRATIO  ! ratio of rain water content to liquid water content
REAL(KIND=JPRB)    :: ZREFACDEN ! Denominator for eff radius increase factor (Wood, 2000)
REAL(KIND=JPRB)    :: ZREFAC    ! Effective radius increase factor (Wood, 2000) 

!* Indices and Arrays linked to prognostic aerosols
INTEGER(KIND=JPIM) :: IACTAERO, IRADLP

INTEGER(KIND=JPIM) :: JAERSS, JAERDD , JAEROM, JAERBC, JAERSU
REAL(KIND=JPRB)    :: ZINVFC0, ZINVFCT
REAL(KIND=JPRB)    :: ZCCN0(KLON), ZCCN(KLON,KLEV), ZMAER(KLON,KLEV,5), ZMAERMN(5),&
  & ZAERO(KLON,KLEV,KAERO), ZQSAT(KLON,KLEV)    , ZRE_LIQ(KLON,KLEV)  , ZRHO(KLON,KLEV),&
  & ZRHCL(KLON,KLEV)      , ZAERSTR(KLON,KLEV)

REAL(KIND=JPRB)    :: ZAERTAUL(KLON,KLEV,16), ZAEROMGL(KLON,KLEV,16), ZAERASYL(KLON,KLEV,16),&
  & ZAERTAULJ(KLON,KLEV,12),ZAERTAUSJ(KLON,KLEV,12)
REAL(KIND=JPRB)    :: ZAERTAUS(KLON,KLEV,14), ZAEROMGS(KLON,KLEV,14), ZAERASYS(KLON,KLEV,14)

INTEGER(KIND=JPIM) :: IHWMCH

!* Cloud overlap decorrelation lengths, in km, as a function of
!* latitude
REAL(KIND=JPRB)    :: ZDECORR_CF(KLON), ZDECORR_CW(KLON)

!* Cosine of latitude
REAL(KIND=JPRB)    :: ZCLAT(KLON)

!* Arrays linked to minimum ice particle size
REAL(KIND=JPRB)    :: ZMINICE(KLON)

! Robert's Cloud generator
!  type(randomNumberStream) :: stream

! Raisanen and Barker's generator
!  REAL(KIND=JPRB) :: RLC_CF(klev), RLC_CW(klev), ZDROPS(klev), ZNUEFF(klev), ZRELW0(klev), ZREIW0(klev)
!  REAL(KIND=JPRB) :: ZF(klon,klev), ZRHOF(klon,klev), ZSIGQCW(klon,klev), ZCLDF(klev)
!  REAL(KIND=JPRB) :: YCF(JPGSW,klev), YQLW(JPGSW,klev), YQIW(JPGSW,klev), YRELW(JPGSW,klev), YREIW(JPGSW,klev)

LOGICAL :: LLCONST_RE, LLMAXRAN, LLPPH, LLPRINT

!Local switches to activate SPP perturbations
LOGICAL :: LLZDECORR, LLZSIGQCW, LLZRADEFFL, LLZRADEFFI

!Local SPP parameters and array integers
REAL(KIND=JPRB)    :: ZFACT
REAL(KIND=JPRB)    :: ZFACT1(KLON)
INTEGER(KIND=JPIM) :: IPZDECORR, IPZSIGQCW, IPZRADEFFLP, IPZRADEFFIP
INTEGER(KIND=JPIM) :: IPN
TYPE(SPP_PERT)     :: PN1ZDECORR, PN1ZSIGQCW, PN1ZRADEFFL, PN1ZRADEFFI


INTEGER(KIND=JPIM) :: IDBUG

!-- end of McICA specific -----------------------

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE, ZHOOK_HANDLE_AEROSOL, ZHOOK_HANDLE_CLOUD

!     -----------------------------------------------------------------

!#include "lw.intfb.h"
#include "satur.intfb.h"
#include "aer_rrtm.intfb.h"
#include "rrtm_rrtm_140gp.intfb.h"
#include "srtm_srtm_224gp.intfb.h"
#include "cloud_overlap_decorr_len.intfb.h"

!-- McICA specific -----------------------
#include "mcica_cld_gen.intfb.h"
#include "rrtm_rrtm_140gp_mcica.intfb.h"
#include "srtm_srtm_224gp_mcica.intfb.h"
!-- McICA specific -----------------------

!     -----------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('RADLSWR',0,ZHOOK_HANDLE)
ASSOCIATE(YDEPHLI=>YDMODEL%YRML_PHY_SLIN%YREPHLI,YDECLDP=>YDMODEL%YRML_PHY_EC%YRECLDP,YDRIP=>YDMODEL%YRML_GCONF%YRRIP, &
 & YDERAD=>YDMODEL%YRML_PHY_RAD%YRERAD, YDESWRT=>YDMODEL%YRML_PHY_RAD%YRESWRT, &
 & YDECLD=>YDMODEL%YRML_PHY_EC%YRECLD,YGFL=>YDMODEL%YRML_GCONF%YGFL,YDEAERATM=>YDMODEL%YRML_PHY_RAD%YREAERATM, &
 & YDSPP_CONFIG=>YDMODEL%YRML_GCONF%YRSPP_CONFIG)
ASSOCIATE(NACTAERO=>YGFL%NACTAERO, &
 & LAERCCN=>YDEAERATM%LAERCCN, LAERCSTR=>YDEAERATM%LAERCSTR, &
 & LAERRRTM=>YDEAERATM%LAERRRTM, &
 & RDECORR_CF=>YDECLD%RDECORR_CF, RDECORR_CW=>YDECLD%RDECORR_CW, &
 & NAERCLD=>YDECLDP%NAERCLD, RCCNOM=>YDECLDP%RCCNOM, RCCNSS=>YDECLDP%RCCNSS, &
 & RCCNSU=>YDECLDP%RCCNSU, &
 & LPHYLIN=>YDEPHLI%LPHYLIN, &
 & LAERADCLI=>YDERAD%LAERADCLI, LAPPROXLWUPDATE=>YDERAD%LAPPROXLWUPDATE, &
 & LCCNL=>YDERAD%LCCNL, LCCNO=>YDERAD%LCCNO, LDIFFC=>YDERAD%LDIFFC, &
 & LRRTM=>YDERAD%LRRTM, NAERMACC=>YDERAD%NAERMACC, NDECOLAT=>YDERAD%NDECOLAT, &
 & NICEOPT=>YDERAD%NICEOPT, NLIQOPT=>YDERAD%NLIQOPT, NMCICA=>YDERAD%NMCICA, &
 & NMCVAR=>YDERAD%NMCVAR, NMINICE=>YDERAD%NMINICE, NOVLP=>YDERAD%NOVLP, &
 & NRADIP=>YDERAD%NRADIP, NRADLP=>YDERAD%NRADLP, NSW=>YDERAD%NSW, &
 & RCCNLND=>YDERAD%RCCNLND, RCCNSEA=>YDERAD%RCCNSEA, RMINICE=>YDERAD%RMINICE, &
 & RRE2DE=>YDERAD%RRE2DE, &
 & NSTART=>YDRIP%NSTART)
!     -----------------------------------------------------------------

LLPRINT=.FALSE.

IDBUG=0
IDBUG=IDBUG+1

IF (LHOOK) CALL DR_HOOK('RADLSWR:AEROSOL',0,ZHOOK_HANDLE_AEROSOL)

!*         1.     PREPARATORY CALCULATIONS RELATED TO AEROSOLS
!                 --------------------------------------------

!-- climatological aerosols
!   1 = land   = organic + sulphate
!   2 = sea    = sea salt
!   3 = desert = dust
!   4 = urban  = black carbon
!   5 = volcanic
!   6 = stratospheric background

IF (NAERMACC == 0) THEN
  IACTAERO=0
ENDIF
DO JAER=1,6
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      ZAER(JL,JAER,JK)=PAER(JL,JAER,JK)
    ENDDO
  ENDDO
ENDDO


!-- prognostic or MACC-derived climatological aerosols (in ZAERO)
!-  to be accounted in radiation and/or cloud effective radius

ZAERTAUL(:,:,:) = 0._JPRB
ZAEROMGL(:,:,:) = 0._JPRB
ZAERASYL(:,:,:) = 0._JPRB
ZAERTAUS(:,:,:) = 0._JPRB
ZAEROMGS(:,:,:) = 0._JPRB
ZAERASYS(:,:,:) = 0._JPRB

!-- GEMS/MACC configuration with NACTAERO prognostic erosols
IF (NACTAERO > 0 ) THEN
  IACTAERO=NACTAERO
  ICLIM=0
ENDIF

!-- MACC-derived aerosol mmr climatology
IF (NACTAERO == 0 .AND. NAERMACC == 1 ) THEN
  IACTAERO=NMCVAR
  ICLIM=1
ENDIF

IF (LAERCCN .OR. LAERRRTM .OR. NAERMACC == 1) THEN
  DO JAER=1,IACTAERO
    DO JK=1,KLEV
      DO JL=KIDIA,KFDIA
        ZAERO(JL,JK,JAER)=MAX(0._JPRB,PAERO(JL,JK,JAER))
      ENDDO
    ENDDO
  ENDDO
ENDIF

IF (LHOOK) CALL DR_HOOK('RADLSWR:AEROSOL',1,ZHOOK_HANDLE_AEROSOL)

! prepare SPP if on
! SPP: Assign mean of perturbation distributions, ZMU_* (see
! spp_mod.F90 / get_spp_conf.F90)
IF (YDSPP_CONFIG%LSPP) THEN
  
  IPN = YDSPP_CONFIG%PPTR%ZDECORR
  LLZDECORR= IPN > 0
  IF (LLZDECORR) THEN
    PN1ZDECORR = YDSPP_CONFIG%SM%PN(IPN)
    IPZDECORR=PN1ZDECORR%MP_RADGRID
  ENDIF

  IPN = YDSPP_CONFIG%PPTR%ZSIGQCW
  LLZSIGQCW= IPN > 0
  IF (LLZSIGQCW) THEN
    PN1ZSIGQCW = YDSPP_CONFIG%SM%PN(IPN)
    IPZSIGQCW = PN1ZSIGQCW%MP_RADGRID
  ENDIF

  IPN = YDSPP_CONFIG%PPTR%ZRADEFFL
  LLZRADEFFL= IPN > 0
  IF (LLZRADEFFL) THEN
    PN1ZRADEFFL = YDSPP_CONFIG%SM%PN(IPN)
    IPZRADEFFLP=PN1ZRADEFFL%MP_RADGRID
  ENDIF

  IPN = YDSPP_CONFIG%PPTR%ZRADEFFI
  LLZRADEFFI= IPN > 0
  IF (LLZRADEFFI) THEN
    PN1ZRADEFFI = YDSPP_CONFIG%SM%PN(IPN)
    IPZRADEFFIP=PN1ZRADEFFI%MP_RADGRID
  ENDIF

ELSE
  LLZDECORR  = .FALSE.
  LLZSIGQCW  = .FALSE.
  LLZRADEFFL = .FALSE.
  LLZRADEFFI = .FALSE.
ENDIF

!          1.1    INITIALISATION OF FLUXES AND DIAGNOSTIC OUTPUT FIELDS 
!                 -----------------------------------------------------

ZRG=1.0_JPRB*RPLRG/RG
ZFCUP(:,:) = 0.0_JPRB
ZFSUP(:,:) = 0.0_JPRB
ZFCDWN(:,:) = 0.0_JPRB
ZFSDWN(:,:) = 0.0_JPRB
ZFCDWN(:,KLEV+1)=REPLOG
ZFSDWN(:,KLEV+1)=REPLOG
ZFLUX(:,:,:) = 0.0_JPRB
ZFLUC(:,:,:) = 0.0_JPRB
ZFSUX(:,:,:) = 0.0_JPRB
ZFSUC(:,:,:) = 0.0_JPRB
ZFSDNN(:) = 0.0_JPRB
ZFSDNV(:) = 0.0_JPRB
ZFCDNN(:) = 0.0_JPRB
ZFCDNV(:) = 0.0_JPRB
ZFSUPN(:) = 0.0_JPRB
ZFSUPV(:) = 0.0_JPRB
ZFCUPN(:) = 0.0_JPRB
ZFCUPV(:) = 0.0_JPRB

ZPSOL(KIDIA:KFDIA) = PAPH(KIDIA:KFDIA,KLEV+1)
ZPMB(KIDIA:KFDIA,1) = ZPSOL(KIDIA:KFDIA) *0.01_JPRB
ZDT0(KIDIA:KFDIA) = PTS(KIDIA:KFDIA) - PTH(KIDIA:KFDIA,KLEV+1)
ZEMIT(:) = 0.0_JPRB
ZSUDU(:) = 0.0_JPRB
PSUDU(:) = 0.0_JPRB
PUVDF(:) = 0.0_JPRB
PPARF(:) = 0.0_JPRB
PPARCF(:)= 0.0_JPRB
ZEMIS(:) = PEMIS(:)
ZEMIW(:) = PEMIW(:)
ZMU0(:)  = PMU0(:)

!IDBUG=IDBUG+1
!IF (LLPRINT) WRITE(NULOUT,7001) IDBUG

!     -------------------------------------------------------------------------

!*         1.2    SET-UP INPUT QUANTITIES FOR ALL RADIATION CONFIGURATIONS
!                 --------------------------------------------------------

!-- in-cloud ice and water mixing ratios; cloud ice and liquid contents and paths

DO JK=1,KLEV
  IJK=KLEV+1-JK
  DO JL=KIDIA,KFDIA
    ZCLFR(JL,JK)   = PCLFR(JL,JK)
    ZCLDSW(JL,IJK) = PCLFR(JL,JK)
    ZCLDLW(JL,IJK) = PCLFR(JL,JK)
    ZETA(JL,JK) =PAP(JL,JK)/PAPH(JL,KLEV+1)
    IF (PCLFR(JL,JK) >= 0.001_JPRB) THEN
      ZTEMP=1.0_JPRB/PCLFR(JL,JK)
      ZQIWP(JL,JK)=MAX(0._JPRB, PQIWP(JL,JK)*ZTEMP)
      ZQSWP(JL,JK)=MAX(0._JPRB, PQSWP(JL,JK)*ZTEMP)
! RMF SENS 20091223 Set snow to zero
!      ZQSWP(JL,JK)=0._JPRB
      ZQLWP(JL,JK)=MAX(0._JPRB, PQLWP(JL,JK)*ZTEMP)
      ZQRWP(JL,JK)=MAX(0._JPRB, PQRWP(JL,JK)*ZTEMP)
! RMF SENS 20100224 Set rain to zero
!      ZQRWP(JL,JK)=0._JPRB
      ZIWGKG=ZQIWP(JL,JK)*1000._JPRB
      ZSWGKG=ZQSWP(JL,JK)*1000._JPRB
      ZLWGKG=ZQLWP(JL,JK)*1000._JPRB 
      ZRWGKG=ZQRWP(JL,JK)*1000._JPRB 
    ELSE
      ZCLFR(JL,JK)=0._JPRB
      ZQIWP(JL,JK)=0._JPRB
      ZQSWP(JL,JK)=0._JPRB
      ZQLWP(JL,JK)=0._JPRB
      ZQRWP(JL,JK)=0._JPRB
      ZIWGKG=0._JPRB
      ZSWGKG=0._JPRB
      ZLWGKG=0._JPRB
      ZRWGKG=0._JPRB
    ENDIF
    ZRAINT(JL)=0._JPRB
    ZDPOG=PDP(JL,JK)*ZRG
    ZFIWP(JL,JK) = ZIWGKG*ZDPOG
    ZFSWP(JL,JK) = ZSWGKG*ZDPOG
    ZFLWP(JL,JK) = ZLWGKG*ZDPOG
    ZFRWP(JL,JK) = ZRWGKG*ZDPOG
    ZPODT = PAP(JL,JK)/(RD*PT(JL,JK))
    ZRHO(JL,JK) = ZPODT 
    ZIWC(JL,JK) = ZIWGKG*ZPODT
    ZSWC(JL,JK) = ZSWGKG*ZPODT
    ZLWC(JL,JK) = ZLWGKG*ZPODT
    ZRWC(JL,JK) = ZRWGKG*ZPODT
  ENDDO
ENDDO

! Reciprocal of pi
ZRPI=1.0_JPRB/RPI


! Compute the cloud overlap decorrelation lengths, in km, for cloud
! boundaries and cloud water content according to the parameterization
! specified by NDECOLAT
CALL CLOUD_OVERLAP_DECORR_LEN(YDECLD, KIDIA, KFDIA, KLON, PGEMU, NDECOLAT, &
     &                         PDECORR_LEN_EDGES_KM=ZDECORR_CF, &
     &                         PDECORR_LEN_WATER_KM=ZDECORR_CW)
!SPP: perturb vertical decorrelation lengths if LLZDECORR=T
IF (LLZDECORR) THEN
  DO JL=KIDIA,KFDIA
    ZFACT=EXP(PN1ZDECORR%MU(1)+PN1ZDECORR%XMAG(1)*PPERT(JL,IPZDECORR))
    ZDECORR_CF(JL)=ZDECORR_CF(JL)*ZFACT
    ZDECORR_CW(JL)=ZDECORR_CW(JL)*ZFACT
  ENDDO
ENDIF
 
IRADLP=NRADLP

IF (LHOOK) CALL DR_HOOK('RADLSWR:AEROSOL',0,ZHOOK_HANDLE_AEROSOL)

!-- otherwise, use the prognostic aerosols for computing the effective size 
!   of cloud particles and/or in the computation of aerosol optical thickness

!-- prognostic aerosols
!  PAERO are entered in kg/kg
!  ZMAER are in ug / m3              ( micrograms per m**3 )
!  ZCNN  are in 1/cm3
!  Menon et al., 2002: JAS, 59, 692-713 Eqns 1a, 1b

!-- first aerosol indirect effect: the relevant prognostic aerosols 
!   (OM, SS, SU) are used as cloud condensation nuclei further used
!   to determine the effective radius of droplets in LIQUID water 
!   clouds
ZCCN(KIDIA:KFDIA,1:KLEV)=0._JPRB
ZRE_LIQ(KIDIA:KFDIA,1:KLEV)=0._JPRB
IF ( IACTAERO >= 12 .AND. NAERCLD > 0 ) THEN
  IRADLP=3 
  JAERSS=1
  JAEROM=2
  JAERBC=3
  JAERSU=4
  JAERDD=5
  ZMAERMN(JAERSS)=2.12E-10_JPRB*1.E9_JPRB 
  ZMAERMN(JAERDD)=1.01E-09_JPRB*1.E9_JPRB
  ZMAERMN(JAEROM)=3.05E-11_JPRB*1.E9_JPRB
  ZMAERMN(JAERBC)=3.05E-11_JPRB*1.E9_JPRB
  ZMAERMN(JAERSU)=1.02E-09_JPRB*1.E9_JPRB
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
!-- 1st bin of SS
      ZMAER(JL,JK,1) = ZAERO(JL,JK, 1)*ZRHO(JL,JK)*1.E9_JPRB
!-- hydrophilic OM
      ZMAER(JL,JK,2) = ZAERO(JL,JK, 7)*ZRHO(JL,JK)*1.E9_JPRB
!-- hydrophilic BC
      ZMAER(JL,JK,3) = ZAERO(JL,JK, 9)*ZRHO(JL,JK)*1.E9_JPRB
!-- SO4
      ZMAER(JL,JK,4) = ZAERO(JL,JK,11)*ZRHO(JL,JK)*1.E9_JPRB
!-- 1 st bin of DD
      ZMAER(JL,JK,5) = ZAERO(JL,JK, 4)*ZRHO(JL,JK)*1.E9_JPRB

      ZRELVNT=ZMAER(JL,JK,1)+ZMAER(JL,JK,2)+ZMAER(JL,JK,4)

! CCNxx factors: for org.matter 0.13, sea salt 0.05, and sulphate 0.50
      ZEXPN= RCCNOM*LOG10(MAX(ZMAER(JL,JK,JAEROM),REPSCA)) +&
           & RCCNSS*LOG10(MAX(ZMAER(JL,JK,JAERSS),REPSCA)) +&
           & RCCNSU*LOG10(MAX(ZMAER(JL,JK,JAERSU),REPSCA))
      ZCCN0(JL) = 10.0_JPRB**(2.41 + ZEXPN)

!-- NB: the prognosed CCN is bounded as the diagnostic one
      IF ( ZRELVNT > 0._JPRB ) THEN
        ZCCN(JL,JK) = MIN( 1743._JPRB, MAX( 32._JPRB, ZCCN0(JL) ))
!-- ZRE_LIQ for subsequent radiative computations (in um once multiplied by 1.E+06)
        ZRE_LIQ(JL,JK) = 1.E+06_JPRB*(2.387E-10_JPRB*ZRHO(JL,JK)*ZQLWP(JL,JK)/ZCCN(JL,JK))**0.333_JPRB
      ENDIF
    ENDDO
  ENDDO
ENDIF


IDBUG=IDBUG+1

!-- Prognostic aerosol optical properties for RRTM: LAERRRTM=.T.
!   MACC-derived aerosol climatology for RRTM: LAERADCLI=.T.

IF (IACTAERO >= 12 .AND. (LAERRRTM .OR. LAERADCLI)) THEN

  IFLAG=2
  CALL SATUR (KIDIA, KFDIA, KLON  , 1, KLEV, YDMODEL%YRML_PHY_SLIN%YREPHLI%LPHYLIN, &
    & PAP, PT   , ZQSAT, IFLAG )  

  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      ZRHCL(JL,JK)=PQ(JL,JK)/ZQSAT(JL,JK)
    ENDDO
  ENDDO

!-- In the absence of a prognostic stratospheric aerosols, LAERCSTR=T uses 
!   the operational climatological stratospheric aerosol optical thickness 
!   to derive an equivalent climatological stratospheric aerosol mass mixing 
!   ratio to be added to the prognostic sulphate aerosol.
!-- The H2SO4 stratospheric aerosols are assumed to be at an ambient RH = 50 %
!   Inversion is performed around 0.5 um.

!   For LAERCSTR=F, the climatological stratospheric optical thickness is 
!   added up in aer_rrtm with their own optical properties.

  IF (LAERCSTR) THEN 
    ZINVFC0=RG/(1.E+03_JPRB*ALF_SU(6,7))
    DO JK=1,KLEV
      DO JL=KIDIA,KFDIA
        ZINVFCT = ZINVFC0 / PDP(JL,JK) !AB fixed a bug here. 
                                       !There was another RG multiplying ZINVFC0
        ZAERSTR(JL,JK) = ZAER(JL,6,JK) * ZINVFCT
        ZAERO(JL,JK,11)=ZAERO(JL,JK,11)+ZAERSTR(JL,JK)
      ENDDO
    ENDDO
  ENDIF

!-- computes the LW and SW optical properties for the GEMS/MACC prognostic aerosols
!-- also adds the AOD for stratospheric and tropospheric aerosols, passed via PAER
  CALL AER_RRTM&
    & ( YDEAERATM,YDMODEL%YRML_PHY_AER,YDMODEL%YRML_GCONF, &
    & KIDIA   , KFDIA   , KLON    , KLEV    , IACTAERO, NSTART  , NSTEP, ICLIM,&
    & PAPH    , ZAERO   , ZAER, ZRHCL   ,&
    & ZAERTAUL, ZAERASYL, ZAEROMGL, ZAERTAUS, ZAERASYS, ZAEROMGS,&
    &   ZAERTAULJ,ZAERTAUSJ)

  IDBUG=IDBUG+1
  IF (NSTEP == NSTART) THEN
    WRITE(NULOUT,FMT='(" AERTAUL ",2F7.4,16E10.3)') PGEMU(KIDIA),PGELAM(KIDIA),(ZAERTAUL(KIDIA,KLEV,JAER),JAER=1,16)
    WRITE(NULOUT,FMT='(" AERTAUS ",2F7.4,16E10.3)') PGEMU(KIDIA),PGELAM(KIDIA),(ZAERTAUS(KIDIA,KLEV,JAER),JAER=1,14)
  ENDIF

ENDIF

IF (LHOOK) CALL DR_HOOK('RADLSWR:AEROSOL',1,ZHOOK_HANDLE_AEROSOL)

IF (LHOOK) CALL DR_HOOK('RADLSWR:CLOUD_EFFECTIVE_RADIUS',0,ZHOOK_HANDLE_CLOUD)

!     -------------------------------------------------------------------------

!*         1.3    EFFECTIVE RADII AND PARTICLE SIZE
!                 ---------------------------------

ZREFDE = RRE2DE 
ZDEFRE = 1.0_JPRB / ZREFDE

DO JL=KIDIA,KFDIA
  ZCLAT(JL) = COS(ASIN(PGEMU(JL)))
  IF (NMINICE == 0) THEN
    ! Constant ice effective radius
    ZMINICE(JL) = RMINICE
  ELSEIF (NMINICE >= 1) THEN
    ! Ice effective radius varies with latitude, smaller at poles
    ZMINICE(JL) = 20._JPRB+(RMINICE-20._JPRB)*ZCLAT(JL)
  ENDIF
ENDDO

DO JK=1,KLEV
  DO JL = KIDIA,KFDIA
! --- EFFECTIVE RADIUS FOR WATER, ICE AND RAIN PARTICLES

! very old parametrization as f(pressure)

    IF (IRADLP == 0) THEN
!-- very old parametrization as f(pressure) ERA-15
      ZRADLP(JL)=10.0_JPRB + (100000.0_JPRB-PAP(JL,JK))*3.5_JPRB

    ELSEIF (IRADLP == 1) THEN
! simple distinction between land (10) and ocean (13) Zhang and Rossow
      IF (PLSM(JL) < 0.5_JPRB) THEN
        ZRADLP(JL)=13.0_JPRB
      ELSE
        ZRADLP(JL)=10.0_JPRB
      ENDIF
      
    ELSEIF (IRADLP == 2) THEN
!--  based on Martin et al., 1994, JAS
      IF (PLSM(JL) < 0.5_JPRB) THEN
        IF (LCCNO) THEN
!          ZASEA=50.0_JPRB
          ZASEA=PCCNO(JL)
        ELSE  
          ZASEA=RCCNSEA
        ENDIF  

        ! Spectral dispersion of the droplet size distribution over sea 
        ! ZD = 0.33_JPRB
        ! Constant (k) relating volume radius to effective radius rv^3=k*re^3
        ! k = (1.0+ZD*ZD)^3)/((1.0+3.0*ZD*ZD)^2)
        ! ZK=(1.0_JPRB+ZD*ZD)**3/(1.0_JPRB+3.0_JPRB*ZD*ZD)**2
        ZK = 0.77_JPRB

        ! Cloud droplet concentration in cm-3 (activated CCN) over ocean
        ZNTOT=-1.15E-03_JPRB*ZASEA*ZASEA+0.963_JPRB*ZASEA+5.30_JPRB

      ELSE
        IF (LCCNL) THEN 
!          ZALND=900.0_JPRB
          ZALND=PCCNL(JL)
        ELSE  
          ZALND=RCCNLND
        ENDIF  

        ! Spectral dispersion (d) of the droplet size distribution over land 
        ! ZD=0.43_JPRB
        ! Constant (k) relating volume radius to effective radius rv^3=k*re^3
        ! k = (1.0+ZD*ZD)^3)/((1.0+3.0*ZD*ZD)^2)
        ! ZK=(1.0_JPRB+ZD*ZD)**3/(1.0_JPRB+3.0_JPRB*ZD*ZD)**2
        ZK = 0.69_JPRB

        ! Cloud droplet concentration in cm-3 (activated CCN) over land
        ZNTOT=-2.10E-04_JPRB*ZALND*ZALND+0.568_JPRB*ZALND-27.9_JPRB

      ENDIF

      !SPP: perturb effective radius (liquid phase) if LLZRADEFFL=T
      IF (LLZRADEFFL) THEN
        ZNTOT=ZNTOT*EXP(PN1ZRADEFFL%MU(1)+PN1ZRADEFFL%XMAG(1)*PPERT(JL,IPZRADEFFLP))
      ENDIF

      ZNUM  = 3.0_JPRB*(ZLWC(JL,JK)+ZRWC(JL,JK))
      ZDEN  = 4.0_JPRB*RPI*ZNTOT*ZK ! *density of water 1000 kg m-3
      ZTEMP = 1.0_JPRB/ZDEN
      ! So for ZTEMP in SI units
      ! g m-3 and cm-3 units cancel out with density of water 10^6/(1000*1000)
      ! Need a factor of 10^6 to convert to microns 
      ! and cubed root is factor of 100 which appears in ZRADLP equation below

      ! Adjustment to Martin_et_al(1994) effective radius 
      ! from Wood(2000) param (Eq. 19)
      IF (ZLWC(JL,JK) > REPSCW) THEN
        ZRLRATIO  = ZRWC(JL,JK)/ZLWC(JL,JK)
        ZKL       = 0.222_JPRB
        ZKRATIO   = (ZKL/ZK)**0.333_JPRB
        ZREFACDEN = 1.0_JPRB+0.2_JPRB*ZKRATIO*ZRLRATIO
        ZREFAC    = ((1.0_JPRB+ZRLRATIO)**0.666_JPRB)/ZREFACDEN
      ELSE
        ZREFAC    = 1.0_JPRB
      ENDIF
 
      IF((ZNUM*ZTEMP) > REPLOG)THEN

        ! R_eff=100*(ZNUM/ZDEN)^0.333, factor of 100 to convert from m to microns
        ZRADLP(JL)=100._JPRB*ZREFAC*EXP(0.333_JPRB*LOG(ZNUM*ZTEMP))

        ! Limit effective radius to within defined range
        ZRADLP(JL)=MAX(ZRADLP(JL), 4.0_JPRB)
        ZRADLP(JL)=MIN(ZRADLP(JL),30.0_JPRB)
      ELSE
        ZRADLP(JL)=4.0_JPRB
      ENDIF

    ELSEIF (IRADLP == 3) THEN
!- effective radius of droplets linked to prognostic aerosol 
!  using Memon et al's CCN (see above)
!      ZRADLP(JL) = ZRE_LIQ(JL,JK)
      ZRADLP(JL) = ZRE_LIQ(JL,JK)
      ZRADLP(JL)=MAX(ZRADLP(JL), 2.0_JPRB)
      ZRADLP(JL)=MIN(ZRADLP(JL),24.0_JPRB)
    ENDIF  

    ZRADLI(JL,JK)=ZRADLP(JL)    
  ENDDO

  DO JL = KIDIA,KFDIA

! diagnosing the ice particle effective radius/diameter

!- ice particle effective radius =f(T) from Liou and Ou (1994)
 
    IF (PT(JL,JK) < RTICE) THEN
      ZTEMPC=PT(JL,JK)-RTT
    ELSE
      ZTEMPC=RTICE-RTT
    ENDIF
    ZRADIP(JL)=326.3_JPRB+ZTEMPC*(12.42_JPRB + ZTEMPC*(0.197_JPRB + ZTEMPC*&
      & 0.0012_JPRB))    
    
    IF (NRADIP == 0) THEN
!-- fixed 40 micron effective radius
      ZRADIP(JL)= 40.0_JPRB
      ZDESR(JL)= ZDEFRE * ZRADIP(JL)
      
    ELSEIF (NRADIP == 1) THEN 
!-- old formulation based on Liou & Ou (1994) temperature (40-130microns)    
      ZRADIP(JL)=MAX(ZRADIP(JL), 40.0_JPRB)
      ZRADIP(JL)=MIN(ZRADIP(JL),130.0_JPRB)
      ZDESR(JL)= ZDEFRE * ZRADIP(JL)
      
    ELSEIF (NRADIP == 2) THEN  
!-- formulation following Jakob, Klein modifications to ice content    
      ZRADIP(JL)=MAX(ZRADIP(JL),30.0_JPRB)
      ZRADIP(JL)=MIN(ZRADIP(JL),60.0_JPRB)
      ZDESR(JL)= ZDEFRE * ZRADIP(JL)
 
    ELSEIF (NRADIP == 3  ) THEN
!- ice particle effective radius =f(T,IWC) from Sun and Rikus (1999)
! revised by Sun (2001)
      IF (ZIWC(JL,JK)+ZSWC(JL,JK) > 0.0_JPRB ) THEN
        ZTEMPC = PT(JL,JK)-83.15_JPRB
        ZTCELS = PT(JL,JK)-RTT
        ZFSR = 1.2351_JPRB +0.0105_JPRB * ZTCELS
! Sun, 2001 (corrected from Sun & Rikus, 1999)
        ZAIWC = 45.8966_JPRB * (ZIWC(JL,JK)+ZSWC(JL,JK))**0.2214_JPRB
        ZBIWC = 0.7957_JPRB * (ZIWC(JL,JK)+ZSWC(JL,JK))**0.2535_JPRB
        ZDESR(JL) = ZFSR * (ZAIWC + ZBIWC*ZTEMPC)

        !SPP: perturb effective radius (ice phase) if LLZRADEFFI=T
        IF (LLZRADEFFI) THEN
          ZDESR(JL)=ZDESR(JL)*EXP(PN1ZRADEFFI%MU(1)+PN1ZRADEFFI%XMAG(1)*PPERT(JL,IPZRADEFFIP))
        ENDIF

        ZDESR(JL) = MIN ( MAX( ZDESR(JL), ZMINICE(JL)), 155.0_JPRB)
        ZRADIP(JL)= ZREFDE * ZDESR(JL)
      ELSE
        ZDESR(JL) = 80.0_JPRB
        ZRADIP(JL)= ZREFDE * ZDESR(JL)
      ENDIF  
    ENDIF  

    ZDESR(JL) = ZRADIP(JL) / ZREFDE
    
    ZRADIC(JL,JK)=ZRADIP(JL)
    ZDESIC(JL,JK)=ZDESR(JL)
  ENDDO
ENDDO

IF (LHOOK) CALL DR_HOOK('RADLSWR:CLOUD_EFFECTIVE_RADIUS',1,ZHOOK_HANDLE_CLOUD)

!     -------------------------------------------------------------------------

!*         2.0    PREPARATION FOR McICA COMPUTATIONS
!                 ----------------------------------

!-- initialise all arrays to zero

! ZFCLW(:,:,:)=0._JPRB
! ZLWFC(:,:,:)=0._JPRB
! ZQILW(:,:,:)=0._JPRB
! ZQLLW(:,:,:)=0._JPRB
! ZTOLW(:,:,:)=0._JPRB

! ZFCSW(:,:,:)=0._JPRB
! ZSWFC(:,:,:)=0._JPRB
! ZQISW(:,:,:)=0._JPRB
! ZQLSW(:,:,:)=0._JPRB
! ZTOSW(:,:,:)=0._JPRB
! ZASSW(:,:,:)=0._JPRB
! ZOMSW(:,:,:)=0._JPRB

IF (LHOOK) CALL DR_HOOK('RADLSWR:CLOUD_GENERATOR',0,ZHOOK_HANDLE_CLOUD)

IF (NMCICA /= 0) THEN 
  DO JL=KIDIA,KFDIA  ! Already defined above
    ZGLON(JL)=PGELAM(JL)*180._JPRB*ZRPI
    ZGLAT(JL)= ASIN(PGEMU(JL)) * 180._JPRB*ZRPI
  ENDDO



!*         2.1    SPECIFYING THE McICA CONFIGURATION
!                 ----------------------------------

  IF (NMCICA == 1) THEN

!-- equivalent PPH_MAXRAN

    LLCONST_RE= .TRUE.
    LLMAXRAN  = .TRUE.
    LLPPH     = .TRUE.

    DO JK=1,KLEV
      DO JL=KIDIA,KFDIA  
        ZLC_CF(JL,JK) = ZDECORR_CF(JL) ! decorrelation depth for cloud fraction
        ZLC_CW(JL,JK) = ZDECORR_CW(JL) ! decorrelation depth for condensate
        ZSIGQCW(JL,JK)= 0._JPRB
!        ZDROPS(JK) = 100.0_JPRB    ! droplet concentration (cm^-3)
!        ZNUEFF(JK) = 0.1_JPRB      ! effective variance of drop size dist'n
!        ZRELW0(JK)  = 10.0_JPRB    ! droplet effective radius
!        ZREIW0(JK)  = 40.0_JPRB    ! ice crystal effective "radius"
      ENDDO
    ENDDO

  ELSEIF (NMCICA == 2) THEN

!-- equivalent GWTSA with Nu=1

    LLCONST_RE= .TRUE.
    LLMAXRAN  = .FALSE.
    LLPPH     = .FALSE.
    IF (LLZSIGQCW) THEN
      DO JL=KIDIA,KFDIA
        ZFACT1(JL)=EXP(PN1ZSIGQCW%MU(1)+PN1ZSIGQCW%XMAG(1)*PPERT(JL,IPZSIGQCW))
      ENDDO
    ENDIF

    DO JK=1,KLEV
      DO JL=KIDIA,KFDIA
        ZLC_CF(JL,JK) = ZDECORR_CF(JL) ! decorrelation depth for cloud fraction
        ZLC_CW(JL,JK) = ZDECORR_CW(JL) ! decorrelation depth for condensate
        ZSIGQCW(JL,JK)= YDERAD%RCLOUD_FRAC_STD ! Formerly was 1, now configurable

        !SPP: perturb fractional stdev of cloud water disctribution if
        !LLZSIGQCW=T
        IF (LLZSIGQCW) THEN
          ZSIGQCW(JL,JK)= ZFACT1(JL)*ZSIGQCW(JL,JK)
        ENDIF
!        ZDROPS(JK) = 0.0_JPRB! droplet concentration (cm^-3)
!        ZNUEFF(JK) = 0.0_JPRB! effective variance of drop size dist'n
!        ZRELW0(JK)  = 10.0_JPRB! droplet effective radius
!        ZREIW0(JK)  = 40.0_JPRB! ice crystal effective "radius"
      ENDDO
    ENDDO
  ENDIF

!     -----------------------------------------------------------------

!*         2.2    CLOUD PARAMETERS FOR McICA COMPUTATIONS
!                 ---------------------------------------

! Before starting, copy original cloud fields for RRTM_LW and RRTM_SW
! done on cloud fraction and condensed ice and liquid water

!-- interface for Raisanen, Cole, Barker's cloud generator
!   it works from top to bottom: 
!   inputs are ZCLFR, ZQLWP, ZQIWP, ZQSWP 
!   outputs are ZQL_SL, ZQI_SL, and ICLDLW an index 0 or 1 for each layer and each g-point

!-- call cloud generator to create ICOLS=JPGLW for LW computations
  CALL MCICA_CLD_GEN&
    &( YDMODEL%YRML_PHY_RAD%YREMCICA,YDRIP, &
    &ZCLFR  , ZQLWP, ZQIWP, ZQSWP, ZLC_CF, ZLC_CW, ZGLAT , ZGLON,&
    &ZSIGQCW, PT   , PAP , KLON  , KIDIA , KFDIA, KLEV, JPGLW,&
    &LLPPH   , LLMAXRAN,&
    & ICLDLW , ZQL_SL, ZQI_SL, ZQS_SL )

  DO JK=1,KLEV
    IJK=KLEV+1-JK
    DO JL=KIDIA,KFDIA
      IF (ICLDLW(JL) /= 0) THEN
        IPRINT(JL,JK)=0
      ENDIF
    ENDDO
    DO JCOLS=1,JPGLW
      DO JL=KIDIA,KFDIA
        IF (ICLDLW(JL) /= 0) THEN
          ZFCLW(JL,JCOLS,JK)=0._JPRB
          IF (ZQL_SL(JL,JK,JCOLS)+ZQI_SL(JL,JK,JCOLS)+ZQS_SL(JL,JK,JCOLS) > 0._JPRB) THEN
            IFCLW(JL,JCOLS,JK)=1
            IPRINT(JL,JK)=IPRINT(JL,JK)+IFCLW(JL,JCOLS,JK)
            ZFCLW(JL,JCOLS,JK)=1._JPRB
            ZLWFC(JL,JCOLS,IJK)=1._JPRB
            ZQILW(JL,JCOLS,JK)=ZQI_SL(JL,JK,JCOLS)
            ZQSLW(JL,JCOLS,JK)=ZQS_SL(JL,JK,JCOLS)
            ZQLLW(JL,JCOLS,JK)=ZQL_SL(JL,JK,JCOLS)
          ELSE
            ZLWFC(JL,JCOLS,IJK)=0._JPRB
            ZQILW(JL,JCOLS,JK)=0._JPRB
            ZQSLW(JL,JCOLS,JK)=0._JPRB
            ZQLLW(JL,JCOLS,JK)=0._JPRB
          ENDIF

        ELSE
          ZFCLW(JL,JCOLS,JK)=0._JPRB
          ZLWFC(JL,JCOLS,IJK)=0._JPRB
          ZQILW(JL,JCOLS,JK)=0._JPRB
          ZQSLW(JL,JCOLS,JK)=0._JPRB
          ZQLLW(JL,JCOLS,JK)=0._JPRB
        ENDIF
      ENDDO
    ENDDO
  ENDDO

!-- call cloud generator to create ICOLS=JPGSW for SW computations
  CALL MCICA_CLD_GEN&
    &( YDMODEL%YRML_PHY_RAD%YREMCICA,YDRIP, &
    &ZCLFR  , ZQLWP, ZQIWP, ZQSWP, ZLC_CF, ZLC_CW, ZGLAT , ZGLON,&
    &ZSIGQCW, PT   , PAP , KLON  , KIDIA , KFDIA, KLEV, JPGSW,&
    &LLPPH   , LLMAXRAN,&
    & ICLDSW , ZQL_SS, ZQI_SS, ZQS_SS )

  DO JK=1,KLEV
    IJK=KLEV+1-JK 
    DO JL=KIDIA,KFDIA
      IF (ICLDSW(JL) /= 0) THEN
        IPRINT(JL,JK)=0
      ENDIF
    ENDDO
    DO JCOLS=1,JPGSW
      DO JL=KIDIA,KFDIA
        IF (ICLDSW(JL) /= 0) THEN
!! original, maybe wrong 
!!        ZFCSW(JL,JCOLS,IJK)=0._JPRB
!- as for LW
          ZFCSW(JL,JCOLS,JK)=0._JPRB
          IF (ZQL_SS(JL,JK,JCOLS)+ZQI_SS(JL,JK,JCOLS)+ZQS_SS(JL,JK,JCOLS) > 0._JPRB) THEN
            IFCSW(JL,JCOLS,JK)=1
            IPRINT(JL,JK)=IPRINT(JL,JK)+IFCSW(JL,JCOLS,JK)
            ZFCSW(JL,JCOLS,JK)=1._JPRB
            ZSWFC(JL,JCOLS,IJK)=1._JPRB
            ZQISW(JL,JCOLS,JK)=ZQI_SS(JL,JK,JCOLS)
            ZQSSW(JL,JCOLS,JK)=ZQS_SS(JL,JK,JCOLS)
            ZQLSW(JL,JCOLS,JK)=ZQL_SS(JL,JK,JCOLS)
          ELSE
            ZSWFC(JL,JCOLS,IJK)=0._JPRB
            ZQISW(JL,JCOLS,JK)=0._JPRB
            ZQSSW(JL,JCOLS,JK)=0._JPRB
            ZQLSW(JL,JCOLS,JK)=0._JPRB
          ENDIF
        ELSE
          ZFCSW(JL,JCOLS,JK)=0._JPRB
          ZSWFC(JL,JCOLS,IJK)=0._JPRB
          ZQISW(JL,JCOLS,JK)=0._JPRB
          ZQSSW(JL,JCOLS,JK)=0._JPRB
          ZQLSW(JL,JCOLS,JK)=0._JPRB
        ENDIF
      ENDDO
    ENDDO
  ENDDO

  ICLW=JPGLW
  ICSW=JPGSW

ELSE

!-- normal non-McICA calculations
  ILW16=16
  DO JCLW=1,ILW16
    DO JK=1,KLEV
      IJK=KLEV+1-JK
      DO JL=KIDIA,KFDIA
!        ZCLDLW(JL,IJK)   =PCLFR(JL,JK)
!        ZFCLW(JL,JCLW,JK)=PCLFR(JL,JK)
!        ZQILW(JL,JCLW,JK)=ZQIWP(JL,JK)
!        ZQLLW(JL,JCLW,JK)=ZQLWP(JL,JK)
      ENDDO
    ENDDO
  ENDDO

  ISW14=14   
  DO JCSW=1,ISW14   
    DO JK=1,KLEV  
      IJK=KLEV+1-JK
      DO JL=KIDIA,KFDIA
 !       ZCLDSW(JL,IJK)   =PCLFR(JL,JK)
 !       ZFCSW(JL,JCSW,JK)=PCLFR(JL,JK)
 !       ZQISW(JL,JCSW,JK)=ZQIWP(JL,JK)
 !       ZQLSW(JL,JCSW,JK)=ZQLWP(JL,JK)
      ENDDO
    ENDDO
  ENDDO
ENDIF

IF (LHOOK) CALL DR_HOOK('RADLSWR:CLOUD_GENERATOR',1,ZHOOK_HANDLE_CLOUD)

LLPRINT=.FALSE.
!     -----------------------------------------------------------------

!*         3.0    INITIALIZE VARIOUS OTHER FIELDS
!                 -------------------------------

!*         3.1    DIFFUSIVITY FACTOR OR SATELLITE VIEWING ANGLE
!                 ---------------------------------------------

DO JL = KIDIA,KFDIA
  ZVIEW(JL) = DIFF
ENDDO

!*         3.2    SURFACE ALBEDO
!                 --------------

!-- mapping SW[1:6] surface albedo into albedo for SRTM[1:14]
ISW=14
DO JSW=1,ISW
  ISW14=NMPSRTM(JSW)
  DO JL = KIDIA,KFDIA
    ZSALBD(JL,JSW)=PALBD(JL,ISW14)
    ZSALBP(JL,JSW)=PALBP(JL,ISW14)
  ENDDO
ENDDO

!*        3.3     PUTTING FIELDS UPSIDE DOWN
!                 --------------------------

DO JK = 1 , KLEV
  JKP1 = JK + 1
  JKL = KLEV+ 1 - JK
  JKLP1 = JKL + 1
  DO JL = KIDIA,KFDIA
    ZPMB(JL,JK+1)=PAPH(JL,JKL)*0.01_JPRB

!-- ZOZ in cm.atm for SW scheme    
    ZOZ(JL,JK)   = POZON(JL,JKL) * 46.6968_JPRB *ZRG 
    ZTL(JL,JK)=PTH(JL,JKLP1)
    ZTAVE(JL,JK)=PT(JL,JKL)
  ENDDO
ENDDO

DO JL=KIDIA,KFDIA
  ZTL(JL,KLEV+1)= PTH(JL,1)
  ZPMB(JL,KLEV+1) = PAPH(JL,1)*0.01_JPRB
ENDDO


!     ------------------------------------------------------------------
!*         4.     CLOUD AND AEROSOL PARAMETERS FOR REGULAR VERSIONS OF RT CODES
!                 ---------------------------- --------------------------------
IF (LHOOK) CALL DR_HOOK('RADLSWR:CLOUD_OPTICS',0,ZHOOK_HANDLE_CLOUD)

IF (NMCICA == 0) THEN 

!          4.1    INITIALIZE OPTICAL PROPERTIES TO CLEAR SKY VALUES
!                 -------------------------------------------------

  ZTAUCLD(:,:,:)=0.0_JPRB
  ZSTAUC(:,:,:)= 0.0_JPRB
  ZSOMGC(:,:,:)= 1.0_JPRB
  ZSASYC(:,:,:)= 0.0_JPRB

  ISW14=14

  DO JK = 1 , KLEV
    IJK=KLEV+1-JK


!          4.2    CLOUD SHORTWAVE OPTICAL PROPERTIES
!                 ----------------------------------
!   -------------------------
! --+ SW OPTICAL PARAMETERS +  Water clouds after Fouquart (1987)
!   -------------------------  Ice clouds (Ebert, Curry, 1992)

    DO JSW=1,ISW14
      DO JL = KIDIA,KFDIA
        ZTOL=0.0_JPRB
        ZGL =0.0_JPRB
        ZOL =0.0_JPRB
        ZTOI=0.0_JPRB
        ZGI =0.0_JPRB
        ZOI =0.0_JPRB
        ZTOR=0.0_JPRB
        ZGR =0.0_JPRB
        ZOR =0.0_JPRB
  
        IF (ZFLWP(JL,JK) > REPSCW .OR. ZFIWP(JL,JK) > REPSCW&
     &    .OR. ZFSWP(JL,JK) > REPSCW ) THEN
          IF (ZFLWP(JL,JK) > REPSCW ) THEN
            IF (NLIQOPT /= 0 ) THEN
!-- SW: Slingo, 1989
              ZTOL = ZFLWP(JL,JK)*(RSASWA(JSW)+RSASWB(JSW)/ZRADLI(JL,JK))
              ZGL  = RSASWE(JSW)+RSASWF(JSW)*ZRADLI(JL,JK)
              ZOL  = 1. - RSASWC(JSW)-RSASWD(JSW)*ZRADLI(JL,JK)
            ELSE          
!-- SW: Fouquart, 1991
              ZTOL = ZFLWP(JL,JK)*(RSYFWA(JSW)+RSYFWB(JSW)/ZRADLI(JL,JK))
              ZGL  = RSYFWF(JSW)
              ZOL  = RSYFWC(JSW)-RSYFWD(JSW)*EXP(-RSYFWE(JSW)*ZTOL)
            ENDIF 
          ENDIF
  
          IF ((ZFIWP(JL,JK)+ZFSWP(JL,JK)) > REPSCW ) THEN
            IF (NICEOPT <= 1) THEN
!-- SW: Ebert-Curry          
              ZTOI = (ZFIWP(JL,JK)+ZFSWP(JL,JK))*(RSECIA(JSW)+RSECIB(JSW)/ZRADIC(JL,JK))
              ZGI  = RSECIE(JSW)+RSECIF(JSW)*ZRADIC(JL,JK)
              ZOI  = 1.0_JPRB - RSECIC(JSW)-RSECID(JSW)*ZRADIC(JL,JK)
            
            ELSEIF (NICEOPT == 2) THEN  
!-- SW: Fu-Liou 1993
              Z1RADI = 1.0_JPRB / ZDESIC(JL,JK)
              ZBETAI = RSFLA0(JSW)+Z1RADI* RSFLA1(JSW)
              ZTOI = (ZFIWP(JL,JK)+ZFSWP(JL,JK)) * ZBETAI
              ZOMGI= RSFLB0(JSW)+ZRADIC(JL,JK)*(RSFLB1(JSW) + ZRADIC(JL,JK)&
               & *(RSFLB2(JSW)+ZRADIC(JL,JK)* RSFLB3(JSW) ))              
              ZOI  = 1.0_JPRB - ZOMGI
              ZOMGP= RSFLC0(JSW)+ZRADIC(JL,JK)*(RSFLC1(JSW) + ZRADIC(JL,JK)&
               & *(RSFLC2(JSW)+ZRADIC(JL,JK)* RSFLC3(JSW) ))   
              ZFDEL= RSFLD0(JSW)+ZRADIC(JL,JK)*(RSFLD1(JSW) + ZRADIC(JL,JK)&
               & *(RSFLD2(JSW)+ZRADIC(JL,JK)* RSFLD3(JSW) ))   
              ZGI  = ((1.0_JPRB-ZFDEL)*ZOMGP + ZFDEL*3.0_JPRB) *0.33333_JPRB
            
            ELSEIF (NICEOPT == 3) THEN  
!-- SW: Fu 1996
              Z1RADI = 1.0_JPRB / ZDESIC(JL,JK)
              ZBETAI = RSFUA0(JSW)+Z1RADI* RSFUA1(JSW)
              ZTOI = (ZFIWP(JL,JK)+ZFSWP(JL,JK)) * ZBETAI
              ZOMGI= RSFUB0(JSW)+ZDESIC(JL,JK)*(RSFUB1(JSW) + ZDESIC(JL,JK)&
               & *(RSFUB2(JSW)+ZDESIC(JL,JK)* RSFUB3(JSW) ))              
              ZOI  = 1.0_JPRB - ZOMGI
              ZGI  = RSFUC0(JSW)+ZDESIC(JL,JK)*(RSFUC1(JSW) + ZDESIC(JL,JK)&
               & *(RSFUC2(JSW)+ZDESIC(JL,JK)* RSFUC3(JSW) ))   
              ZGI  = MIN(1.0_JPRB, ZGI) 
            ENDIF
          ENDIF

!  - MIX of WATER and ICE CLOUDS
          ZTAUMX= ZTOL + ZTOI + ZTOR
          ZOMGMX= ZTOL*ZOL + ZTOI*ZOI + ZTOR*ZOR
          ZASYMX= ZTOL*ZOL*ZGL + ZTOI*ZOI*ZGI + ZTOR*ZOR*ZGR
        
          ZASYMX= ZASYMX/ZOMGMX
          ZOMGMX= ZOMGMX/ZTAUMX

! --- SW FINAL CLOUD OPTICAL PARAMETERS

          ZSTAUC(JL,JSW,IJK) = ZTAUMX
          ZSOMGC(JL,JSW,IJK) = ZOMGMX
          ZSASYC(JL,JSW,IJK) = ZASYMX
        ENDIF
      ENDDO
    ENDDO


!          4.3    CLOUD LONGWAVE OPTICAL PROPERTIES FOR RRTM
!                 ------------------------------------------

!   -------------------------
! --+ LW OPTICAL PARAMETERS +  Water (and Ice) from Savijarvi (1998)
!   -------------------------  Ice clouds (Ebert, Curry, 1992)

! No need for a fixed diffusivity factor, accounted for spectrally below
! The detailed spectral structure does not require defining upward and
! downward effective optical properties

    ILW16=16

    DO JRTM=1,ILW16
      DO JL = KIDIA,KFDIA
        ZRSALD = 0.0_JPRB
        ZRSAID = 0.0_JPRB
        ZRSAIO = 0.0_JPRB
        
        IF (ZFLWP(JL,JK)+ZFIWP(JL,JK)+ZFSWP(JL,JK) > REPSCW) THEN
    
          IF (NLIQOPT == 0 .OR. NLIQOPT >= 3) THEN
! water cloud total emissivity a la Smith and Shi (1992)
            ZMULTL=1.2_JPRB-0.006_JPRB*ZRADLI(JL,JK)
            ZRSALD= 0.144_JPRB*ZMULTL * 0.60240_JPRB
            
          ELSEIF (NLIQOPT == 1) THEN
! water cloud spectral emissivity a la Savijarvi (1997)
            ZRSALD= YDESWRT%RHSAVI(JRTM,1) + ZRADLI(JL,JK)&
             & *(YDESWRT%RHSAVI(JRTM,2) + ZRADLI(JL,JK)*YDESWRT%RHSAVI(JRTM,3))  
             
          ELSEIF (NLIQOPT == 2) THEN
! water cloud spectral emissivity a la Lindner and Li (2000)
            Z1RADL = 1.0_JPRB / ZRADLI(JL,JK)
            ZEXTCF = YDESWRT%RLILIA(JRTM,1)+ZRADLI(JL,JK)*YDESWRT%RLILIA(JRTM,2)+ Z1RADL*&
             & (YDESWRT%RLILIA(JRTM,3) + Z1RADL*(YDESWRT%RLILIA(JRTM,4) + Z1RADL*&
             & YDESWRT%RLILIA(JRTM,5) ))  
            Z1MOMG = YDESWRT%RLILIB(JRTM,1) + Z1RADL*YDESWRT%RLILIB(JRTM,2)&
             & + ZRADLI(JL,JK) *(YDESWRT%RLILIB(JRTM,3) + ZRADLI(JL,JK)*YDESWRT%RLILIB(JRTM,4) )
            ZRSALD = Z1MOMG * ZEXTCF
          ENDIF  
         
          IF (NICEOPT == 0) THEN
! ice cloud spectral emissivity a la Smith & Shi (1992)
            ZMULTI=1.2_JPRB-0.006_JPRB*ZRADIC(JL,JK)
            ZRSAID= 0.103_JPRB*ZMULTI * 0.60240_JPRB
            
          ELSEIF (NICEOPT == 1) THEN
! ice cloud spectral emissivity a la Ebert-Curry (1992)
            ZRSAID= YDESWRT%REBCUH(JRTM)+YDESWRT%REBCUG(JRTM)/ZRADIC(JL,JK)
          
          ELSEIF (NICEOPT == 2) THEN
! ice cloud spectral emissivity a la Fu & Liou (1993)
            Z1RADI= 1.0_JPRB / ZDESIC(JL,JK)
            ZRSAID = YDESWRT%RFULIO(JRTM,1) + Z1RADI&
             & *(YDESWRT%RFULIO(JRTM,2) + Z1RADI*YDESWRT%RFULIO(JRTM,3))  
             
          ELSEIF (NICEOPT == 3) THEN
! ice cloud spectral emissivity a la Fu et al (1998) including 
! parametrisation for LW scattering effect a la Chou et al (1999)
            Z1RADI = 1.0_JPRB / ZDESIC(JL,JK)
            ZRSAIE = YDESWRT%RFUETA(JRTM,1) + Z1RADI&
             & *(YDESWRT%RFUETA(JRTM,2) + Z1RADI*YDESWRT%RFUETA(JRTM,3))   
            ZRSAIA = Z1RADI*(YDESWRT%RFUETB(JRTM,1) +ZDESIC(JL,JK)*( YDESWRT%RFUETB(JRTM,2) +ZDESIC(JL,JK)*&
             & ( YDESWRT%RFUETB(JRTM,3) +ZDESIC(JL,JK)* YDESWRT%RFUETB(JRTM,4))))
            ZRSAIG = YDESWRT%RFUETC(JRTM,1) +ZDESIC(JL,JK)*( YDESWRT%RFUETC(JRTM,2) +ZDESIC(JL,JK)*&
             & ( YDESWRT%RFUETC(JRTM,3) +ZDESIC(JL,JK)* YDESWRT%RFUETC(JRTM,4))) 
            ZRSAIG = MIN(MAX(ZRSAIG,-1.0_JPRB),1.0_JPRB)
            IF (YDERAD%LFU_LW_ICE_OPTICS_BUG) THEN
              ! Incorrect form
              ZRSAIO = ZRSAIA/ZRSAIE
            ELSE
              ZRSAIO = 1.0_JPRB - ZRSAIA/ZRSAIE
            ENDIF

            ZRSAIO = MIN(MAX(ZRSAIO,0.0_JPRB),1.0_JPRB)
            ZRSAIF = 0.5_JPRB + ZRSAIG*( 0.3738_JPRB + ZRSAIG*( 0.0076_JPRB + ZRSAIG*0.1186_JPRB ) )
            ZRSAID = (1.0_JPRB - ZRSAIO * ZRSAIF) * ZRSAIE

          ENDIF    
         
          ZTAUD = ZRSALD*ZFLWP(JL,JK)+ZRSAID*(ZFIWP(JL,JK)+ZFSWP(JL,JK))

! Diffusivity correction within clouds a la Savijarvi
          IF (LDIFFC) THEN
            ZDIFFD=MIN(MAX(1.517_JPRB-0.156_JPRB*LOG(ZTAUD) , 1.0_JPRB) , 2.0_JPRB)
          ELSE
            ZDIFFD=1.66_JPRB
          ENDIF
          ZTAUCLD(JL,IJK,JRTM) = ZTAUD*ZDIFFD
        ENDIF
        
      ENDDO
    ENDDO
  ENDDO

  ZTOLW(:,:,:) = 0.0_JPRB
  ZTOSW(:,:,:) = 0.0_JPRB
  ZASSW(:,:,:) = 0.0_JPRB
  ZOMSW(:,:,:) = 0.0_JPRB

ELSE


!          5.1    INITIALIZE OPTICAL PROPERTIES TO CLEAR SKY VALUES
!                 -------------------------------------------------

  ZTOLW(:,:,:) = 0.0_JPRB
! ZTOSW(:,:,:) = 0.0_JPRB
! ZASSW(:,:,:) = 0.0_JPRB
! ZOMSW(:,:,:) = 1.0_JPRB
  DO JK=1,KLEV
    DO JCSW=ICSW+1,JPGSW
      DO JL = KIDIA,KFDIA
        ZTOSW(JL,JCSW,JK) = 0.0_JPRB
        ZASSW(JL,JCSW,JK) = 0.0_JPRB
        ZOMSW(JL,JCSW,JK) = 1.0_JPRB
      ENDDO
    ENDDO
  ENDDO

!          5.2    SHORTWAVE OPTICAL PROPERTIES
!                 ----------------------------

  DO JCSW=1,ICSW
    JSW=NGBSW(JCSW)-15

    DO JK = 1 , KLEV
      IJK=KLEV+1-JK

!          5.3    CLOUD ICE AND LIQUID CONTENT AND PATH
!                 -------------------------------------

      DO JL = KIDIA,KFDIA

! --- LIQUID WATER CONTENT (g.m-3) AND LIQUID WATER PATH (g.m-2)
       IF (ZFCSW(JL,JCSW,JK) > REPSC ) THEN
         ZLWGKG=MAX(ZQLSW(JL,JCSW,JK)*1000.0_JPRB,0.0_JPRB)
         ZIWGKG=MAX(ZQISW(JL,JCSW,JK)*1000.0_JPRB,0.0_JPRB)
         ZSWGKG=MAX(ZQSSW(JL,JCSW,JK)*1000.0_JPRB,0.0_JPRB)
      ELSE  
         ZLWGKG=0.0_JPRB
         ZIWGKG=0.0_JPRB
         ZSWGKG=0.0_JPRB
       ENDIF
!        ZLWGKG=0.0_JPRB
!        ZIWGKG=0.0_JPRB
!        ZSWGKG=0.0_JPRB
!!        ZRWGKG=0.0_JPRB
!        IF (ZFCSW(JL,JCSW,JK) > REPSC ) THEN
!          ZLWGKG=500.0_JPRB*(ABS(ZQLSW(JL,JCSW,JK))+ZQLSW(JL,JCSW,JK))
!          ZIWGKG=500.0_JPRB*(ABS(ZQISW(JL,JCSW,JK))+ZQISW(JL,JCSW,JK))
!          ZSWGKG=500.0_JPRB*(ABS(ZQSSW(JL,JCSW,JK))+ZQSSW(JL,JCSW,JK))
!!          ZRWGKG=500.0_JPRB*(ABS(ZQRSW(JL,JCSW,JK))+ZQRSW(JL,JCSW,JK))
!        ENDIF
        ZRWGKG=0.0_JPRB
!       ZRAINT(JL)=0.0_JPRB  !  Not used in this loop

        ZDPOG=PDP(JL,JK)*ZRG
        ZFLWP(JL,JK)= ZLWGKG*ZDPOG
        ZFIWP(JL,JK)= ZIWGKG*ZDPOG
        ZFSWP(JL,JK)= ZSWGKG*ZDPOG
!       ZFRWP(JL,JK)= ZRWGKG*ZDPOG   ! Always zero
!       ZPODT=PAP(JL,JK)/(RD*PT(JL,JK)) ! Not used in this loop
!       ZLWC(JL,JK)=ZLWGKG*ZPODT    ! Not used in this loop
!       ZIWC(JL,JK)=ZIWGKG*ZPODT    ! Not used in this loop
!       ZSWC(JL,JK)=ZSWGKG*ZPODT    ! Not used in this loop
      ENDDO


!          5.4    CLOUD SHORTWAVE OPTICAL PROPERTIES
!                 ----------------------------------
!   -------------------------
! --+ SW OPTICAL PARAMETERS +  Water clouds after Fouquart (1987)
!   -------------------------  Ice clouds (Ebert, Curry, 1992)

      DO JL = KIDIA,KFDIA
        ZTOSW(JL,JCSW,IJK) = 0.0_JPRB
        ZOMSW(JL,JCSW,IJK) = 1.0_JPRB
        ZASSW(JL,JCSW,IJK) = 0.0_JPRB
      ENDDO
      DO JL = KIDIA,KFDIA
        ZTOL=0.0_JPRB
        ZGL =0.0_JPRB
        ZOL =0.0_JPRB
        ZTOI=0.0_JPRB
        ZGI =0.0_JPRB
        ZOI =0.0_JPRB
        ZTOR=0.0_JPRB
        ZGR =0.0_JPRB
        ZOR =0.0_JPRB

        IF (ZFLWP(JL,JK) > REPSCW .OR. ZFIWP(JL,JK) > REPSCW&
     &   .OR. ZFSWP(JL,JK) > REPSCW ) THEN
          IF (ZFLWP(JL,JK) > REPSCW ) THEN
            IF (NLIQOPT /= 0 ) THEN
!-- SW: Slingo, 1989
              ZTOL = ZFLWP(JL,JK)*(RSASWA(JSW)+RSASWB(JSW)/ZRADLI(JL,JK))
              ZGL  = RSASWE(JSW)+RSASWF(JSW)*ZRADLI(JL,JK)
              ZOL  = 1. - RSASWC(JSW)-RSASWD(JSW)*ZRADLI(JL,JK)
            ELSE          
!-- SW: Fouquart, 1991
              ZTOL = ZFLWP(JL,JK)*(RSYFWA(JSW)+RSYFWB(JSW)/ZRADLI(JL,JK))
              ZGL  = RSYFWF(JSW)
              ZOL  = RSYFWC(JSW)-RSYFWD(JSW)*EXP(-RSYFWE(JSW)*ZTOL)
            ENDIF 
          ENDIF
  
          IF ((ZFIWP(JL,JK)+ZFSWP(JL,JK)) > REPSCW ) THEN
            IF (NICEOPT <= 1) THEN
!-- SW: Ebert-Curry          
              ZTOI = (ZFIWP(JL,JK)+ZFSWP(JL,JK))*(RSECIA(JSW)+RSECIB(JSW)/ZRADIC(JL,JK))
              ZGI  = RSECIE(JSW)+RSECIF(JSW)*ZRADIC(JL,JK)
              ZOI  = 1.0_JPRB - RSECIC(JSW)-RSECID(JSW)*ZRADIC(JL,JK)
            
            ELSEIF (NICEOPT == 2) THEN  
!-- SW: Fu-Liou 1993
              Z1RADI = 1.0_JPRB / ZDESIC(JL,JK)
              ZBETAI = RSFLA0(JSW)+Z1RADI* RSFLA1(JSW)
              ZTOI = (ZFIWP(JL,JK)+ZFSWP(JL,JK)) * ZBETAI
              ZOMGI= RSFLB0(JSW)+ZRADIC(JL,JK)*(RSFLB1(JSW) + ZRADIC(JL,JK)&
               & *(RSFLB2(JSW)+ZRADIC(JL,JK)* RSFLB3(JSW) ))              
              ZOI  = 1.0_JPRB - ZOMGI
              ZOMGP= RSFLC0(JSW)+ZRADIC(JL,JK)*(RSFLC1(JSW) + ZRADIC(JL,JK)&
               & *(RSFLC2(JSW)+ZRADIC(JL,JK)* RSFLC3(JSW) ))   
              ZFDEL= RSFLD0(JSW)+ZRADIC(JL,JK)*(RSFLD1(JSW) + ZRADIC(JL,JK)&
               & *(RSFLD2(JSW)+ZRADIC(JL,JK)* RSFLD3(JSW) ))   
              ZGI  = ((1.0_JPRB-ZFDEL)*ZOMGP + ZFDEL*3.0_JPRB) *0.33333_JPRB
            
            ELSEIF (NICEOPT == 3) THEN  
!-- SW: Fu 1996
              Z1RADI = 1.0_JPRB / ZDESIC(JL,JK)
              ZBETAI = RSFUA0(JSW)+Z1RADI* RSFUA1(JSW)
              ZTOI = (ZFIWP(JL,JK)+ZFSWP(JL,JK)) * ZBETAI
              ZOMGI= RSFUB0(JSW)+ZDESIC(JL,JK)*(RSFUB1(JSW) + ZDESIC(JL,JK)&
               & *(RSFUB2(JSW)+ZDESIC(JL,JK)* RSFUB3(JSW) ))              
              ZOI  = 1.0_JPRB - ZOMGI
              ZGI  = RSFUC0(JSW)+ZDESIC(JL,JK)*(RSFUC1(JSW) + ZDESIC(JL,JK)&
               & *(RSFUC2(JSW)+ZDESIC(JL,JK)* RSFUC3(JSW) ))   
              ZGI  = MIN(1.0_JPRB, ZGI)
            ENDIF
          ENDIF

!  - MIX of WATER and ICE CLOUDS
          ZTAUMX= ZTOL + ZTOI + ZTOR
          ZOMGMX= ZTOL*ZOL + ZTOI*ZOI + ZTOR*ZOR
          ZASYMX= ZTOL*ZOL*ZGL + ZTOI*ZOI*ZGI + ZTOR*ZOR*ZGR
        
          ZASYMX= ZASYMX/ZOMGMX
          ZOMGMX= ZOMGMX/ZTAUMX

! --- SW FINAL CLOUD OPTICAL PARAMETERS

          ZTOSW(JL,JCSW,IJK) = ZTAUMX
          ZOMSW(JL,JCSW,IJK) = ZOMGMX
          ZASSW(JL,JCSW,IJK) = ZASYMX
        ENDIF
      ENDDO
    ENDDO
  ENDDO

!          5.3    CLOUD LONGWAVE OPTICAL PROPERTIES FOR RRTM
!                 ------------------------------------------

!   -------------------------
! --+ LW OPTICAL PARAMETERS +  Water (and Ice) from Savijarvi (1998)
!   -------------------------  Ice clouds (Ebert, Curry, 1992)

! No need for a fixed diffusivity factor, accounted for spectrally below
! The detailed spectral structure does not require defining upward and
! downward effective optical properties

  DO JCLW=1,ICLW
    JRTM=NGB(JCLW)
    DO JK=1,KLEV
      IJK=KLEV+1-JK

!          5.4    CLOUD ICE AND LIQUID CONTENT AND PATH
!                 -------------------------------------

      DO JL = KIDIA,KFDIA

! --- LIQUID WATER CONTENT (g.m-3) AND LIQUID WATER PATH (g.m-2)
        IF (ZFCLW(JL,JCLW,JK) > REPSC ) THEN
          ZLWGKG=MAX(ZQLLW(JL,JCLW,JK)*1000.0_JPRB,0.0_JPRB)
          ZIWGKG=MAX(ZQILW(JL,JCLW,JK)*1000.0_JPRB,0.0_JPRB)
          ZSWGKG=MAX(ZQSLW(JL,JCLW,JK)*1000.0_JPRB,0.0_JPRB)
        ELSE  
          ZLWGKG=0.0_JPRB
          ZIWGKG=0.0_JPRB
          ZSWGKG=0.0_JPRB
        ENDIF
        ZRWGKG=0.0_JPRB
!       ZRAINT(JL)=0.0_JPRB ! Not used in this loop

        ZDPOG=PDP(JL,JK)*ZRG
        ZFLWP(JL,JK)= ZLWGKG*ZDPOG
        ZFIWP(JL,JK)= ZIWGKG*ZDPOG
        ZFSWP(JL,JK)= ZSWGKG*ZDPOG
!       ZFRWP(JL,JK)= ZRWGKG*ZDPOG   ! Always zero
!       ZPODT=PAP(JL,JK)/(RD*PT(JL,JK)) ! Not used in this loop
!       ZLWC(JL,JK)=ZLWGKG*ZPODT  ! Not used in this loop
!       ZIWC(JL,JK)=ZIWGKG*ZPODT  ! Not used in this loop
!       ZSWC(JL,JK)=ZSWGKG*ZPODT  ! Not used in this loop
      ENDDO


!          5.5    CLOUD SHORTWAVE OPTICAL PROPERTIES
!                 ---------------------------------- 

      DO JL = KIDIA,KFDIA
        ZRSALD = 0.0_JPRB
        ZRSAID = 0.0_JPRB
        
        IF (ZFLWP(JL,JK)+ZFIWP(JL,JK)+ZFSWP(JL,JK) > REPSCW) THEN

          IF (NLIQOPT == 0 .OR. NLIQOPT >= 3) THEN
! water cloud total emissivity a la Smith and Shi (1992)
            ZMULTL=1.2_JPRB-0.006_JPRB*ZRADLI(JL,JK)
            ZRSALD= 0.144_JPRB*ZMULTL *0.60240_JPRB
            
          ELSEIF (NLIQOPT == 1) THEN
! water cloud spectral emissivity a la Savijarvi (1997)
            ZRSALD= YDESWRT%RHSAVI(JRTM,1) + ZRADLI(JL,JK)&
             & *(YDESWRT%RHSAVI(JRTM,2) + ZRADLI(JL,JK)*YDESWRT%RHSAVI(JRTM,3))  
             
          ELSEIF (NLIQOPT == 2) THEN
! water cloud spectral emissivity a la Lindner and Li (2000)
            Z1RADL = 1.0_JPRB / ZRADLI(JL,JK)
          ZEXTCF = YDESWRT%RLILIA(JRTM,1)+ZRADLI(JL,JK)*YDESWRT%RLILIA(JRTM,2)+ Z1RADL*&
             & (YDESWRT%RLILIA(JRTM,3) + Z1RADL*(YDESWRT%RLILIA(JRTM,4) + Z1RADL*&
             & YDESWRT%RLILIA(JRTM,5) ))  
            Z1MOMG = YDESWRT%RLILIB(JRTM,1) + Z1RADL*YDESWRT%RLILIB(JRTM,2)&
             & + ZRADLI(JL,JK) *(YDESWRT%RLILIB(JRTM,3) + ZRADLI(JL,JK)*YDESWRT%RLILIB(JRTM,4) )
            ZRSALD = Z1MOMG * ZEXTCF
          ENDIF  
         
          IF (NICEOPT == 0) THEN
! ice cloud spectral emissivity a la Smith & Shi (1992)
            ZMULTI=1.2_JPRB-0.006_JPRB*ZRADIC(JL,JK)
            ZRSAID= 0.103_JPRB*ZMULTI *0.60240_JPRB
            
          ELSEIF (NICEOPT == 1) THEN
! ice cloud spectral emissivity a la Ebert-Curry (1992)
            ZRSAID= YDESWRT%REBCUH(JRTM)+YDESWRT%REBCUG(JRTM)/ZRADIC(JL,JK)
          
          ELSEIF (NICEOPT == 2) THEN
! ice cloud spectral emissivity a la Fu & Liou (1993)
            Z1RADI= 1.0_JPRB / ZDESIC(JL,JK)
            ZRSAID = YDESWRT%RFULIO(JRTM,1) + Z1RADI&
             & *(YDESWRT%RFULIO(JRTM,2) + Z1RADI*YDESWRT%RFULIO(JRTM,3))  
             
          ELSEIF (NICEOPT == 3) THEN
! ice cloud spectral emissivity a la Fu et al (1998) including 
! parametrisation for LW scattering effect a la Chou et al (1999)
            Z1RADI = 1.0_JPRB / ZDESIC(JL,JK)
            ZRSAIE = YDESWRT%RFUETA(JRTM,1) + Z1RADI&
             & *(YDESWRT%RFUETA(JRTM,2) + Z1RADI*YDESWRT%RFUETA(JRTM,3))   
            ZRSAIA = Z1RADI*(YDESWRT%RFUETB(JRTM,1) +ZDESIC(JL,JK)*( YDESWRT%RFUETB(JRTM,2) +ZDESIC(JL,JK)*&
             & ( YDESWRT%RFUETB(JRTM,3) +ZDESIC(JL,JK)* YDESWRT%RFUETB(JRTM,4))))
            ZRSAIG = YDESWRT%RFUETC(JRTM,1) +ZDESIC(JL,JK)*( YDESWRT%RFUETC(JRTM,2) +ZDESIC(JL,JK)*&
             & ( YDESWRT%RFUETC(JRTM,3) +ZDESIC(JL,JK)* YDESWRT%RFUETC(JRTM,4))) 
            ZRSAIG = MIN(MAX(ZRSAIG,-1.0_JPRB),1.0_JPRB)
            IF (YDERAD%LFU_LW_ICE_OPTICS_BUG) THEN
              ! Incorrect form
              ZRSAIO = ZRSAIA/ZRSAIE
            ELSE
              ZRSAIO = 1.0_JPRB - ZRSAIA/ZRSAIE
            ENDIF

            ZRSAIO = MIN(MAX(ZRSAIO,0.0_JPRB),1.0_JPRB)
            ZRSAIF = 0.5_JPRB + ZRSAIG*( 0.3738_JPRB + ZRSAIG*( 0.0076_JPRB + ZRSAIG*0.1186_JPRB ) )
            ZRSAID = (1.0_JPRB - ZRSAIO * ZRSAIF) * ZRSAIE

          ENDIF    
         
          ZTAUD = ZRSALD*ZFLWP(JL,JK)+ZRSAID*(ZFIWP(JL,JK)+ZFSWP(JL,JK))

! Diffusivity correction within clouds a la Savijarvi
          IF (LDIFFC) THEN
            ZDIFFD=MIN(MAX(1.517_JPRB-0.156_JPRB*LOG(ZTAUD) , 1.0_JPRB) , 2.0_JPRB)
          ELSE
            ZDIFFD=1.66_JPRB
          ENDIF
          ZTOLW(JL,IJK,JCLW) = ZTAUD*ZDIFFD
        ENDIF
        
      ENDDO
    ENDDO
  ENDDO

ENDIF

IF (LHOOK) CALL DR_HOOK('RADLSWR:CLOUD_OPTICS',1,ZHOOK_HANDLE_CLOUD)

NUAER = NUA
NTRAER = NTRA

!     ------------------------------------------------------------------

!*         6.     CALL LONGWAVE RADIATION CODE
!                 ----------------------------

!*         6.1    FULL LONGWAVE RADIATION COMPUTATIONS
!                 ------------------------------------

IF (.NOT.LPHYLIN) THEN
  IF ( .NOT. LRRTM) THEN

    CALL ABOR1('OLD LW SCHEME NOT AVAILABLE IN McICA CONFIGURATION')

  ELSE

!*         6.2    FULL LONGWAVE RADIATION COMPUTATIONS - RRTM
!                 ------------------------------------   ----

!  i)  pass ZOZN (ozone mass mixing ratio) to RRTM; remove pressure
!      weighting applied to POZON in driverMC (below)
!  ii) pass ZEMIS and ZEMIW to RRTM; return ZEMIT from RRTM
!  iii)pass ZTAUCLD, cloud optical depths (water+ice) to RRTM, 
!      computed from equations above
!  iv) pass ECRT arrays to RRTM arrays in interface routine ECRTATM
!      in module rrtm_ecrt.f

    DO JL = KIDIA,KFDIA
      DO JK = 1, KLEV
        ZOZN(JL,JK) = POZON(JL,JK)/PDP(JL,JK)
      ENDDO
    ENDDO

    IF (NMCICA == 0) THEN

     CALL RRTM_RRTM_140GP&
      & (YDDIMV, YDERAD, KIDIA , KFDIA , KLON  , KLEV,&
      & YDESWRT%RAER, ZAER  , PAPH  , PAP,&
      & PTS   , PTH   , PT,&
      & ZEMIS , ZEMIW,&
      & PQ    , PCO2  , PCH4 , PN2O   , PNO2, PC11, PC12, PC22, PCL4, ZOZN, ZCLDSW, ZTAUCLD,&
      & ZEMIT , ZFLUX , ZFLUC, ZTCLEAR&
      & )  

    ELSE

      ICOLS=JPGLW

      CALL RRTM_RRTM_140GP_MCICA&
       &(YDDIMV, YDEAERATM,YDERAD,YGFL, KIDIA, KFDIA, KLON , KLEV, ICOLS, ICLDLW ,&
       &YDESWRT%RAER, ZAER , PAPH , PAP  , ZAERTAUL, ZAERASYL, ZAEROMGL,&
       &PTS  , PTH  , PT   ,&
       &ZEMIS, ZEMIW,&
       &PQ   , PCO2 , PCH4 , PN2O, PNO2 , PC11, PC12, PC22, PCL4,&
       &ZOZN , ZLWFC , ZTOLW , ZCLDLW,&
       &ZEMIT, ZFLUX, ZFLUC,&
       &  ZLWDERIVATIVE)  

    ENDIF

  ENDIF

! Reverse the ordering of levels in order to pass the arrays back to
! the rest of the model
  DO JK = 1 , KLEV+1
    IK=KLEV+2-JK
    DO JL=KIDIA,KFDIA
      PFLT(JL,JK) =-ZFLUX(JL,2,IK) - ZFLUX(JL,1,IK)
      PFCT(JL,JK) =-ZFLUC(JL,2,IK) - ZFLUC(JL,1,IK)
    ENDDO
! If the partial derivatives are not required then we should not write
! to PLwDerivative since it may not have been allocated in
! RADINTG
    IF (LAPPROXLWUPDATE) THEN
       DO JL=KIDIA,KFDIA
          PLWDERIVATIVE(JL,JK) = ZLWDERIVATIVE(JL,IK)
       ENDDO
    ENDIF
  ENDDO

ELSE
  ZEMIT (:)   = 0.0_JPRB
  ZFLUX(:,:,:)= 0.0_JPRB
  ZFLUC(:,:,:)= 0.0_JPRB
ENDIF

!     ------------------------------------------------------------------

!*         7.     CALL SHORTWAVE RADIATION CODE
!                 -----------------------------

ZRMUZ=0.0_JPRB
DO JL = KIDIA,KFDIA
  ZRMUZ = MAX (ZRMUZ, ZMU0(JL))
ENDDO

IF (NMCICA == 0) THEN

  CALL SRTM_SRTM_224GP&
    & (YDDIMV, YDERAD,YDMODEL%YRML_PHY_RAD%YRERDI,YDMODEL%YRML_PHY_MF%YRPHY3, &
    & KIDIA , KFDIA , KLON  , KLEV  , ISW , NOVLP,&
    & ZAER  , ZSALBD, ZSALBP, PAPH  , PAP ,&
    & PTS   , PTH   , PT    ,&
    & PQ    , PCO2  , PCH4  , PN2O  , PNO2, ZOZN , ZMU0  ,&
    & ZCLDSW, ZSTAUC, ZSASYC, ZSOMGC,&
    &   ZFSUX , ZFSUC , ZFUVF ,ZSUDU,&
    &   ZFDIR , ZCDIR , ZFDIF ,ZCDIF )

ELSE

  ICOLS=JPGSW

  IF (LLPRINT .AND. NSTEP == NSTART) THEN
    JL=INT((KIDIA+KFDIA)/2)
    WRITE(NULOUT,FMT='(" albedoD ",I4,14F7.4)') JL,(ZSALBD(JL,JSW),JSW=1,ISW)
    WRITE(NULOUT,FMT='(" albedoP ",I4,14F7.4)') JL,(ZSALBP(JL,JSW),JSW=1,ISW)
    DO JK=1,KLEV
      WRITE(NULOUT,FMT='(" oAER_RRTM ToL",2I4,16F7.4)') JL,JK,(ZAERTAUL(JL,JK,JRTM),JRTM=1,16)
    ENDDO       
    DO JK=1,KLEV
      WRITE(NULOUT,FMT='(" oAER_RRTM ToL",2I4,16F7.4)') JL,JK,(ZAEROMGL(JL,JK,JRTM),JRTM=1,16)
    ENDDO       
    DO JK=1,KLEV
      WRITE(NULOUT,FMT='(" oAER_RRTM ToL",2I4,16F7.4)') JL,JK,(ZAERASYL(JL,JK,JRTM),JRTM=1,16)
    ENDDO       
    DO JK=1,KLEV
      WRITE(NULOUT,FMT='(" oAER_RRTM ToS",2I4,16F7.4)') JL,JK,(ZAERTAUS(JL,JK,JSW),JSW=1,14)
    ENDDO       
    DO JK=1,KLEV
      WRITE(NULOUT,FMT='(" oAER_RRTM OmS",2I4,16F7.4)') JL,JK,(ZAEROMGS(JL,JK,JSW),JSW=1,14)
    ENDDO       
    DO JK=1,KLEV
      WRITE(NULOUT,FMT='(" oAER_RRTM AsS",2I4,16F7.4)') JL,JK,(ZAERASYS(JL,JK,JSW),JSW=1,14)
    ENDDO       
    DO JK=1,KLEV
      WRITE(NULOUT,FMT='(" TraceGases" ,2I4,8E12.4)') JL,JK,PT(JL,JK),PTH(JL,JK),PQ(JL,JK),PCO2(JL,JK),PCH4(JL,JK),PN2O(JL,JK),&
        & PNO2(JL,JK),ZOZN(JL,JK)
    ENDDO
    WRITE(NULOUT,FMT='("PTH ",20F7.1)') (PTH(JL,KLEV+1),JL=KIDIA,KFDIA)
    WRITE(NULOUT,FMT='("PTS ",20F7.1)') (PTS(JL),JL=KIDIA,KFDIA)
    WRITE(NULOUT,FMT='("MU0 ",20F7.1)') (ZMU0(JL),JL=KIDIA,KFDIA)
  ENDIF

  CALL SRTM_SRTM_224GP_MCICA&
   & (YDDIMV, YDMODEL%YRML_PHY_RAD,YGFL,YDMODEL%YRML_PHY_MF%YRPHY3, &
   & KIDIA , KFDIA , KLON  , KLEV , ISW , ICOLS , ICLDSW ,&
   & ZAER  , ZSALBD, ZSALBP, PAPH , PAP , ZAERTAUS, ZAERASYS, ZAEROMGS,&
   & PTS   , PTH   , PT    ,&
   & PQ    , PCO2  , PCH4  , PN2O , PNO2, ZOZN  , ZMU0  ,&
   & ZSWFC , ZTOSW , ZASSW , ZOMSW,&
   & ZFSUX , ZFSUC , ZFUVF , ZFUVC,ZPARF, ZPARCF, ZSUDU ,&
   &   ZFDIR , ZCDIR , ZFDIF , ZCDIF, PSWDIFFUSEBAND, PSWDIRECTBAND )  

ENDIF
 
IDBUG=IDBUG+1
IF (LLPRINT .AND. NSTEP == NSTART) THEN
  JL=INT((KIDIA+KFDIA)/2)
  DO JK=1,KLEV+1
    WRITE(NULOUT,FMT='(" LongW ",2I4,4F12.2)') JL,JK,ZFLUX(JL,1,JK),ZFLUX(JL,2,JK),ZFLUC(JL,1,JK),ZFLUC(JL,2,JK)
  ENDDO
  DO JK=1,KLEV+1
    WRITE(NULOUT,FMT='(" Solar ",2I4,4F12.2)') JL,JK,ZFSUX(JL,1,JK),ZFSUX(JL,2,JK),ZFSUC(JL,1,JK),ZFSUC(JL,2,JK)
  ENDDO
ENDIF
!     ------------------------------------------------------------------

!*         8.     FILL UP THE MODEL NET LW AND SW RADIATIVE FLUXES
!                 ------------------------------------------------
IHWMCH=MAX(1, (KFDIA-KIDIA+1)/10)
DO JK = 1 , KLEV+1
  IK=KLEV+2-JK
  DO JL = KIDIA,KFDIA
    PFLS(JL,JK) = ZFSUX(JL,2,JK) - ZFSUX(JL,1,JK)
    PFCS(JL,JK) = ZFSUC(JL,2,JK) - ZFSUC(JL,1,JK)
  ENDDO
ENDDO

DO JL = KIDIA,KFDIA
  PFRSOD(JL) = ZFSUX(JL,2,KLEV+1)
  PFRSODC(JL)= ZFSUC(JL,2,KLEV+1)
  PFRTED(JL) = -ZFLUX(JL,2,1)
  PFRTEDC(JL)= -ZFLUC(JL,2,1)
  PEMIT (JL) = ZEMIT(JL)
  PSUDU (JL) = ZSUDU(JL)
  PUVDF (JL) = ZFUVF(JL)
  PPARF (JL) = ZPARF(JL)
  PPARCF(JL) = ZPARCF(JL)
  PTINCF(JL) = PRII0 * ZMU0(JL)
  PFDIR (JL) = ZFDIR(JL,KLEV+1)
  PFDIF (JL) = ZFDIF(JL,KLEV+1)
  PCDIR (JL) = ZCDIR(JL,KLEV+1)
ENDDO
!IF (LLPRINT .AND. NSTEP <= 5) THEN
!  WRITE(NULOUT,9401) (PUVDF(JL),JL=KIDIA,KFDIA)
!9401 FORMAT(1X,'RADLSWR ',20E10.3)
!ENDIF

IDBUG=IDBUG+1
IF (LLPRINT) WRITE(NULOUT,7011) IDBUG
7011 FORMAT(1X,'RADLSWR part_9 after fluxes IDBUG=',I3,' going out of RADLSWR')

!     --------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('RADLSWR',1,ZHOOK_HANDLE)
END SUBROUTINE RADLSWR
