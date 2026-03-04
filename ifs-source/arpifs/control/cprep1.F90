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

SUBROUTINE CPREP1(YDGEOMETRY,YDSURF,YDEPHY,YDML_GCONF,YDPHY1,PTEFRCL)

! **** *CPREP1*  - Transform GRIB file to FA file (NCONF=901)

! Purpose.
! --------
! 901 : GRIB file format

! **   Interface.
! ----------
! *CALL* *CPREP1

! Explicit arguments :
! --------------------

! Implicit arguments :
! --------------------
! None

! Method.
! -------
! See documentation

! Externals.
! ----------

! Reference.
! ----------
! ECMWF Research Department documentation of the IFS
! Note de travail ARPEGE NR 17 et  xx

! Author.
! -------
! Mats Hamrud and Philippe Courtier  *ECMWF*
! Original : 87-10-15

! Modifications.
! --------------
! R. El Khatib : 01-08-07 Pruning options
! Modified 02-02-04 (Patrick SAEZ *CNRM/GMAP*)
!  Corrections on conversion GRIB DATE to ARPEGE DATE
! Modified 03-02-25 (Patrick SAEZ *CNRM/GMAP*)
!  Add LLOLDSWL,LLHUPDG,LLCONTROL options
! Modified 03-04-10 (Patrick SAEZ *CNRM/GMAP*)
!  Add computations of 'SURFALBEDO NEIGE','SURFDENSIT.NEIGE','SURFALBEDO.SOLNU',
! 'SURFALBEDO.VEG' and LLALBEDO2 switch.
! Modified 03-04-17 (Patrick SAEZ *CNRM/GMAP*)
!  compute estimated values of SURFRESERV.GLACE and PROFRESERV.GLACE
! R. El Khatib : 03-08-18 Roughness lengths not packed
! M.Hamrud      01-Oct-2003 CY28 Cleaning
! D. Paradis & R. El Khatib : 04-07-22 GRIBEX
! R. Brozkova : 05-10-31 Set forecast range in the file date.
! R. Brozkova: 28-07-06 allow the treatment of the upper-air only by
!                           relaxing the LSM test via LLCONTROL switch
! Modified 07-11-21 (Patrick SAEZ *CNRM/GMAP*)
!   Add aerosols and ozone climatological fields
!   Compute RESERV_EAU with new SLT MARS field (Soil Moisture Index)
! K. Yessad (Jul 2009): remove CDLOCK + some cleanings
! S.Martinez (2010-10-19): rehabilit ZWP_SAT(0)=PPWSAT_CEP
! P.saez (2010-12-02) : norm doctor for variables and end of GOTO statement
! Mate Mile & JM Audoin (2011-02-01) : use grib-api subroutines
! T.Wilhelmsson 1-July-2013 : Use GRIB_API_INTERFACE
! R. EL Khatib 16-Aug-2016 Bugfix when a parameter was unattended or had an incorrect resolution
! R. El Khatib 09-Sep-2016 LLMODGRIN ; cleanings.
! E.Dutra/G.Arduini (Jan 2018): SP_SG changes, snow multi-layer
! -----------------------------------------------------------------

USE GEOMETRY_MOD           , ONLY : GEOMETRY
USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE SURFACE_FIELDS_MIX     , ONLY : TSURF
USE PARKIND1               , ONLY : JPIM, JPRB
USE YOMHOOK                , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCT0                 , ONLY : CNMEXP, NCONF
USE YOMCST                 , ONLY : RG, R,RPI
USE YOMLUN                 , ONLY : NULOUT, NULNAM, NPOSSH    
USE YOM_GRIB_CODES         , ONLY : NGRBSWVL1, NGRBSD, NGRBSDSL, NGRBSDOR, NGRBSWVL2, NGRBSR,&
 &                                  NGRBSWVL3, NGRBSRC , NGRBSWVL4, NGRBCVL, NGRBCVH, NGRBTVL, NGRBTVH  
USE YOMOPH0                , ONLY : CNMCA
USE YOEPHY                 , ONLY : TEPHY
USE YOMMSC                 , ONLY : NINTLEN
USE YOMPHY1                , ONLY : TPHY1
USE YOMCLI                 , ONLY : SARGN, SSABN, SDEPX
USE GRIB_API_INTERFACE     , ONLY : IGRIB_OPEN_FILE, IGRIB_NEW_FROM_FILE, IGRIB_GET_VALUE,&
  &                                 IGRIB_RELEASE, IGRIB_CLOSE_FILE, JPGRIB_SUCCESS, JPGRIB_END_OF_FILE

IMPLICIT NONE

TYPE(GEOMETRY) ,INTENT(IN)    :: YDGEOMETRY
TYPE(TSURF)    ,INTENT(INOUT) :: YDSURF
TYPE(TEPHY)    ,INTENT(INOUT) :: YDEPHY
TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(INOUT):: YDML_GCONF
TYPE(TPHY1)    ,INTENT(INOUT) :: YDPHY1
REAL(KIND=JPRB),INTENT(IN)    :: PTEFRCL ! max number of input GRIB files to be read

INTEGER(KIND=JPIM),PARAMETER :: JPFILN=10 ! max number of input GRIB files to be read
INTEGER(KIND=JPIM),PARAMETER :: JPSTOIA=4 ! maximum number of intermediate grid-point arrays to be stored.
INTEGER(KIND=JPIM),PARAMETER :: JPNTVMER=1 !valeur sur mer dans PIVEG (genre de surface)
INTEGER(KIND=JPIM),PARAMETER :: JPNTVTER=3 !valeur sur terre dans PIVEG (genre de surface)
INTEGER(KIND=JPIM),PARAMETER :: JPOLDSWL1=140 !value of NGRBSWVL1 before june 2000
INTEGER(KIND=JPIM),PARAMETER :: JPOLDSWL2=171 !value of NGRBSWVL2 before june 2000
INTEGER(KIND=JPIM),PARAMETER :: JPOLDSWL3=184 !value of NGRBSWVL3 before june 2000
INTEGER(KIND=JPIM),PARAMETER :: JPOLDSWL4=237 !value of NGRBSWVL4 before june 2000

REAL(KIND=JPRB), PARAMETER   :: PPWSAT_CEP=0.472_JPRB ! Teneur en eau a la saturation au CEP
REAL(KIND=JPRB), PARAMETER   :: PPWWILT_CEP=0.171_JPRB ! Point de fletrissement au CEP
REAL(KIND=JPRB), PARAMETER   :: PPWCAP_CEP=0.323_JPRB ! Capacite au champ au CEP
REAL(KIND=JPRB), PARAMETER   :: PPPSWL1=0.07_JPRB ! Reservoir superficiel du CEP
REAL(KIND=JPRB), PARAMETER   :: PPPSOL=1.5_JPRB ! Prodondeur pour SWL2 si LLCLIM=.false.
! epaisseur des couches au CEP

REAL(KIND=JPRB), PARAMETER,DIMENSION(4)   :: PPZD = (/ 0.07, 0.21, 0.72, 1.89/)

REAL(KIND=JPRB),    ALLOCATABLE :: ZGRBDAT(:)

INTEGER(KIND=JPIM) :: IDATEF(11)
INTEGER(KIND=JPIM) :: ISEC1(60),ISEC2(22+YDGEOMETRY%YRDIM%NDGLG)
INTEGER(KIND=JPIM) :: IDATEC(11) ! Date of climatological file if LLCLIM
INTEGER(KIND=JPIM) :: IGRBS(JPSTOIA) ! grib codes of the stored fields.
INTEGER(KIND=JPIM) :: IUCLIM,IDEB2,IFIN,IINDEX,ISOIL
INTEGER(KIND=JPIM) :: ISWL1,ISWL2,ISWL3,ISWL4,INBCLIB,INBCGLA

REAL(KIND=JPRB) :: ZFASP(YDGEOMETRY%YRDIM%NSPEC2)    ! Spectral data to write on FA file.
REAL(KIND=JPRB) :: ZFASPS(YDGEOMETRY%YRDIM%NSPEC2)   ! Spectral data to send/receive to reespe/speree.
REAL(KIND=JPRB) :: ZFPDAT(YDGEOMETRY%YRDIM%NSEFRE)   ! Spectral working array.
REAL(KIND=JPRB) :: ZFAGG(YDGEOMETRY%YRGEM%NGPTOTG)   ! Grid-point data to write on FA file.
REAL(KIND=JPRB) :: ZFAGGA(YDGEOMETRY%YRGEM%NGPTOTG,JPSTOIA) ! Grid-point data to store.
REAL(KIND=JPRB) :: ZFASTL2(YDGEOMETRY%YRGEM%NGPTOTG) ! store STL2 to compute PROFTEMPERATURE
REAL(KIND=JPRB) :: ZFASTL3(YDGEOMETRY%YRGEM%NGPTOTG) ! store STL3 to compute PROFTEMPERATURE
REAL(KIND=JPRB) :: ZFASTL4(YDGEOMETRY%YRGEM%NGPTOTG) ! store STL4 to compute PROFTEMPERATURE
REAL(KIND=JPRB) :: ZFALSM(YDGEOMETRY%YRGEM%NGPTOTG)  ! store Land/Sea Mask
REAL(KIND=JPRB) :: ZFZ0G(YDGEOMETRY%YRGEM%NGPTOTG)   ! store SURFZ0.FOIS.G to compute SURFGZ0.THERM
REAL(KIND=JPRB) :: ZFCVL(YDGEOMETRY%YRGEM%NGPTOTG)   ! store CVL (Low Vegetation Cover)
REAL(KIND=JPRB) :: ZFCVH(YDGEOMETRY%YRGEM%NGPTOTG)   ! store CVH (High Vegetation Cover)
REAL(KIND=JPRB) :: ZFTVL(YDGEOMETRY%YRGEM%NGPTOTG)   ! store Type of Low Vegetation
REAL(KIND=JPRB) :: ZFTVH(YDGEOMETRY%YRGEM%NGPTOTG)   ! store Type of High Vegetation

! Array used to compute SURFRESERV.EAU and PROFRESERV.EAU in ISBA scheme
REAL(KIND=JPRB) :: ZFSWL1(YDGEOMETRY%YRGEM%NGPTOTG)  ! store SWL1
REAL(KIND=JPRB) :: ZFSWL2(YDGEOMETRY%YRGEM%NGPTOTG)  ! store SWL2
REAL(KIND=JPRB) :: ZFSWL3(YDGEOMETRY%YRGEM%NGPTOTG)  ! store SWL3
REAL(KIND=JPRB) :: ZFSWL4(YDGEOMETRY%YRGEM%NGPTOTG)  ! store SWL4
! Array used in argument to ACSOLW to find Wsat_isba (ZFALSM is also used)
REAL(KIND=JPRB) :: ZFPARG(YDGEOMETRY%YRGEM%NGPTOTG)  ! store ARGILE
REAL(KIND=JPRB) :: ZFPSAB(YDGEOMETRY%YRGEM%NGPTOTG)  ! store SABLE
REAL(KIND=JPRB) :: ZFPIVEG(YDGEOMETRY%YRGEM%NGPTOTG) ! store TYPE DE SURFACE
REAL(KIND=JPRB) :: ZFPD2(YDGEOMETRY%YRGEM%NGPTOTG)   ! store PD2 (SURFEPAIS.SOL)

REAL(KIND=JPRB) :: ZFPLAI(YDGEOMETRY%YRGEM%NGPTOTG)  ! store LAI (SURFIND.FOLIAIRE)
REAL(KIND=JPRB) :: ZFPRSMIN(YDGEOMETRY%YRGEM%NGPTOTG)  ! store RSMIN (SURFRESI.STO.MIN)
REAL(KIND=JPRB) :: ZFPVEG(YDGEOMETRY%YRGEM%NGPTOTG) ! store ISBA VEGETATION (SURFPROP.VEGETAT)

! Work Array for ACSOLW()
REAL(KIND=JPRB) :: ZFPWFC(YDGEOMETRY%YRGEM%NGPTOTG)  
REAL(KIND=JPRB) :: ZFPWPMX(YDGEOMETRY%YRGEM%NGPTOTG) 
REAL(KIND=JPRB) :: ZFPWSAT(YDGEOMETRY%YRGEM%NGPTOTG) 
REAL(KIND=JPRB) :: ZFPWSMX(YDGEOMETRY%YRGEM%NGPTOTG) 
REAL(KIND=JPRB) :: ZFPWWILT(YDGEOMETRY%YRGEM%NGPTOTG)

! Array used to store SURFTEMPERATURE and PROFTEMPERATURE
REAL(KIND=JPRB) :: ZFSURFT(YDGEOMETRY%YRGEM%NGPTOTG)
REAL(KIND=JPRB) :: ZFPROFT(YDGEOMETRY%YRGEM%NGPTOTG)
REAL(KIND=JPRB) :: ZFGLACE(YDGEOMETRY%YRGEM%NGPTOTG)

! NAMELIST variables for SURFRESERV.GLACE and PROFRESERV.GLACE
REAL(KIND=JPRB) :: TSRESERV1,TSRESERV2,TDELTA1,TDELTA2

REAL(KIND=JPRB),ALLOCATABLE :: ZRVCOV(:)
REAL(KIND=JPRB),ALLOCATABLE :: ZRVLAI(:)
REAL(KIND=JPRB) :: ZRVRSMIN(0:20)
REAL(KIND=JPRB),ALLOCATABLE :: ZRVROOTSA(:,:)

INTEGER(KIND=JPIM) :: IDATRT, IDATRTG, IDATRTS, IDONR,&
 & IERR, IGRBSN, ILENBYT, ILENWOR, ILEVEL,&
 & ILEVTY,  INBARI, INBARP,IVTYPES,&
 & INDEX, INIVEAU, INTIMES, IOK, IPARAM, IREP,&
 & IREPRM, IRET, JB, JC, JFILN,&
 & JGL, JGPTOTG, JINDEX, JM, JN, JNM,&
 & JSTOIA, JV, IIND, NBUF901

INTEGER(KIND=JPIM) :: INGRIG,INGRIB,INBPDG,INBCSP,ISTRON,IPUILA,IDMOPL

LOGICAL :: LLMODGRIN,LL3DMOD, LLFIRAR, LLGRB190,&
 & LLGRB191, LLGRB192, LLGRB193,&
 & LLINPFS, LLOUTFS, LLSCAL, LLWRIF,&
 & LLCLIM, LLLSM, LLCDST, LLSTL4, LLSTL2, LLZ0,&
 & LLSWL1, LLSWL2, LLSWL3, LLSWL4, LLSLT, LLCVL, LLCVH, LLTVL, LLTVH,&
 & LLOLDSWL, LLHUPDG, LLCONTROL, LLALBEDO2, LLSURFT, LLAEROSOL
LOGICAL :: LLFILE_TO_MODEL

LOGICAL :: LLLAISCAL, LLSWI

CHARACTER (LEN = 80) :: CLFAFILN
CHARACTER (LEN = 80) :: CLLIBELLE
CHARACTER :: CLPREF*4,CLSUFF*16
CHARACTER (LEN = 200) ::  CLFILN(JPFILN)

REAL(KIND=JPRB) :: ZGOL, ZMAX, ZMEA, ZMIN, ZSCALE, ZVAL, ZPSOL, ZPGLACE, ZDELTA, ZWSATSLT
REAL(KIND=JPRB) :: ZPGLACE_LAYERS(4), ZPLIQUID_LAYERS(4), ZPSOILT_LAYERS(4)
REAL(KIND=JPRB) :: ZPWWILT_CEP, ZPWCAP_CEP, ZSWI_CEP
REAL(KIND=JPRB) :: ZTETA_MEAN_H, ZTETA_MEAN_L, ZSWI_CEP_H, ZSWI_CEP_L, ZVEGH, ZVEGL
REAL(KIND=JPRB) :: ZRSMIN_LAI_IFS,  ZRSMIN_LAI_H, ZRSMIN_LAI_L ,ZRSMIN_LAI_ISBA, ZLAISCAL
REAL(KIND=JPRB) :: ZEVH, ZEVL ! pseudo evaporations vegetations high and low
REAL(KIND=JPRB) :: ZEVIFS ! pseudo total evaporation for IFS

! WP_SAT, WP_WILT, WP_CAP and Soil Moisture Index
REAL(KIND=JPRB), DIMENSION(0:7) :: ZWP_SAT
REAL(KIND=JPRB), DIMENSION(0:7) :: ZWP_WILT
REAL(KIND=JPRB), DIMENSION(0:7) :: ZWP_CAP
REAL(KIND=JPRB) :: ZFSLT(YDGEOMETRY%YRGEM%NGPTOTG)  ! store SLT

! ** Array of fields to be simply done if LLCLIM
CHARACTER(LEN=20), DIMENSION(20) :: CLLIB =(/&
 & 'SURFALBEDO          ',&
 & 'SURFEMISSIVITE      ','SURFEPAIS.SOL       ','SURFZ0.FOIS.G       ',&
 & 'SURFGZ0.THERM       ','SURFIND.FOLIAIRE    ','SURFIND.VEG.DOMI    ',&
 & 'SURFPROP.ARGILE     ','SURFPROP.SABLE      ','SURFPROP.VEGETAT    ',&
 & 'SURFRESI.STO.MIN    ','SURFALBEDO.SOLNU    ','SURFALBEDO.VEG      ',&
 & 'SURFA.OF.OZONE      ','SURFB.OF.OZONE      ','SURFC.OF.OZONE      ',&
 & 'SURFAEROS.SEA       ','SURFAEROS.LAND      ','SURFAEROS.SOOT      ',&
 & 'SURFAEROS.DESERT    '/)  

! ** Structure pour les champs (constants) qui n'existent pas au CEP
! (depend parfois du Land/Sea Mask)

TYPE CHAMPSCST
CHARACTER(LEN=20) :: CLHAMP    ! Nom du champ
REAL(KIND=JPRB)   :: ZVALTERRE ! Valeur sur terre (ou globale)
REAL(KIND=JPRB)   :: ZVALMER   ! Valeur sur mer
LOGICAL           :: LLISLSM   ! .TRUE. si depend du LS Mask
END TYPE CHAMPSCST

TYPE(CHAMPSCST),DIMENSION(12) :: YL_FIELDCST
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
INTEGER(KIND=JPIM) :: JGRIB, IFILE, IEDITION
CHARACTER(LEN=16) :: CLWORKAROUND ! workaround to a bug in Intel compiler v16.x.x

! -------------------------------------------------------------------

#include "surf_inq.h"
#include "abor1.intfb.h"
#include "acsolw.intfb.h"
#include "openfa.intfb.h"
#include "posnam.intfb.h"
#include "reespe.intfb.h"
#include "speree.intfb.h"
#include "val923.intfb.h"
#include "spreord.intfb.h"

! -------------------------------------------------------------------

#include "nammars.nam.h"

! -------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('CPREP1',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM, &
  & YDDIMV=>YDGEOMETRY%YRDIMV, &
  & YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, YDLAP=>YDGEOMETRY%YRLAP, YDDIMF=>YDML_GCONF%YRDIMF)

ASSOCIATE(RD1=>YDPHY1%RD1, RD2MER=>YDPHY1%RD2MER, RSMAX=>YDPHY1%RSMAX, &
 & NTVMER=>YDPHY1%NTVMER, GCONV=>YDPHY1%GCONV, &
 & NS3D=>YDDIMF%NS3D, NS2D=>YDDIMF%NS2D, &
 & NSEFRE=>YDDIM%NSEFRE, NDGLG=>YDDIM%NDGLG, NSMAX=>YDDIM%NSMAX, &
 & NDGNH=>YDDIM%NDGNH, NSPEC2=>YDDIM%NSPEC2, NDGSAG=>YDDIM%NDGSAG, &
 & NDGENG=>YDDIM%NDGENG, NSPEC=>YDDIM%NSPEC, NDLON=>YDDIM%NDLON, &
 & NFLEVG=>YDDIMV%NFLEVG, &
 & NGPTOTG=>YDGEM%NGPTOTG, NHTYP=>YDGEM%NHTYP, NLOENG=>YDGEM%NLOENG, &
 & YSD_WSD=>YDSURF%YSD_WSD, YSP_SG=>YDSURF%YSP_SG, YSP_SBD=>YDSURF%YSP_SBD, &
 & YSP_SB=>YDSURF%YSP_SB, YSP_RRD=>YDSURF%YSP_RRD, YSD_WS=>YDSURF%YSD_WS, &
 & YSD_WWD=>YDSURF%YSD_WWD, YSD_VPD=>YDSURF%YSD_VPD, YSP_SGD=>YDSURF%YSP_SGD, &
 & YSP_RR=>YDSURF%YSP_RR, YSP_X2D=>YDSURF%YSP_X2D, YSD_VP=>YDSURF%YSD_VP, &
 & YSD_VF=>YDSURF%YSD_VF, YSP_EP=>YDSURF%YSP_EP, YSP_EPD=>YDSURF%YSP_EPD, &
 & YSD_VFD=>YDSURF%YSD_VFD, YSD_WW=>YDSURF%YSD_WW, YSP_X2=>YDSURF%YSP_X2, &
 & YSURF=>YDEPHY%YSURF, NVALUE=>YDLAP%NVALUE, RLAPIN=>YDLAP%RLAPIN)
! -------------------------------------------------------------------

! ** Initialization.
YL_FIELDCST(1)=CHAMPSCST('SURFRESERV.GLACE    ',0._JPRB,   0._JPRB,     .FALSE.)
YL_FIELDCST(2)=CHAMPSCST('PROFRESERV.GLACE    ',0._JPRB,   0._JPRB,     .FALSE.)
YL_FIELDCST(3)=CHAMPSCST('SURFRES.EVAPOTRA    ',1._JPRB,   1._JPRB,     .FALSE.)
YL_FIELDCST(4)=CHAMPSCST('SURFRESERV.INTER    ',0._JPRB,   0._JPRB,     .FALSE.)
YL_FIELDCST(5)=CHAMPSCST('SURFEMISSIVITE      ',0.96_JPRB, 0._JPRB,     .FALSE.)
YL_FIELDCST(6)=CHAMPSCST('SURFEPAIS.SOL       ',PPPSOL,    REAL(RD2MER),.TRUE.)
YL_FIELDCST(7)=CHAMPSCST('SURFIND.FOLIAIRE    ',2._JPRB,   0._JPRB,     .TRUE.)
YL_FIELDCST(8)=CHAMPSCST('SURFIND.VEG.DOMI    ',3._JPRB,   REAL(NTVMER),.TRUE.)
YL_FIELDCST(9)=CHAMPSCST('SURFPROP.ARGILE     ',30._JPRB,  REAL(SARGN), .TRUE.)
YL_FIELDCST(10)=CHAMPSCST('SURFPROP.SABLE      ',45._JPRB, REAL(SSABN), .TRUE.)
YL_FIELDCST(11)=CHAMPSCST('SURFRESI.STO.MIN    ',150._JPRB,REAL(RSMAX), .TRUE.)
YL_FIELDCST(12)=CHAMPSCST('SURFPROP.VEGETAT    ',1._JPRB,  0._JPRB,     .TRUE.)

! values associated with SLT coeff.

CALL SURF_INQ(YSURF,PRWSATM=ZWP_SAT)
ZWP_SAT(0)=PPWSAT_CEP

!presently ZWP_SAT = (/0.0_JPRB, 0.403_JPRB, 0.439_JPRB, 0.430_JPRB, 0.520_JPRB, 0.614_JPRB, &
!          & 0.766_JPRB, 0.439_JPRB/)

CALL SURF_INQ(YSURF,PRWPWPM=ZWP_WILT)

!presently ZWP_WILT = (/0.0_JPRB, 0.059_JPRB, 0.151_JPRB, 0.133_JPRB, 0.279_JPRB, &
!          & 0.335_JPRB, 0.267_JPRB, 0.151_JPRB/)

CALL SURF_INQ(YSURF,PRWCAPM=ZWP_CAP)

!presently ZWP_CAP = (/0.0_JPRB, 0.242_JPRB, 0.346_JPRB, 0.382_JPRB, 0.448_JPRB, &
!          & 0.541_JPRB, 0.662_JPRB, 0.346_JPRB/)

LLSWI=.FALSE.
LLLAISCAL=.FALSE.

LLFIRAR=.TRUE. ! .true. if first grib article read.
LLGRB190=.FALSE. ! .true. if this GRIB code filed has been read from file.
LLGRB191=.FALSE.
LLGRB192=.FALSE.
LLGRB193=.FALSE.
LLMODGRIN=.FALSE. ! .true. if input grib is modified for TSL2-TSL3-TSL4 as in IFS
IGRBSN=0 ! number of intermediate grid-point arrays stored.
LLCLIM=.FALSE. ! .true. if climatological file is used
LLLSM=.FALSE.  ! .true. id Land/sea Mask is read
LLCDST=.FALSE. ! to compute PROFTEMPERATURE whith CDST
LLSTL4=.FALSE. ! to compute PROFTEMPERATURE whith STL4
LLSTL2=.FALSE. ! to compute PROFTEMPERATURE whith STL2
LLZ0=.FALSE.   ! to compute SURFGZ0.THERM whith SURFZ0.FOIS.G if read
LLSWL1=.FALSE.
LLSWL2=.FALSE.
LLSWL3=.FALSE.
LLSWL4=.FALSE.
LLSLT=.FALSE.
LLCVL=.FALSE.
LLCVH=.FALSE.
LLTVL=.FALSE.
LLTVH=.FALSE.
LLOLDSWL=.FALSE. ! used to compute SWLi before june 2000 if true.
LLHUPDG=.FALSE.  ! used to not write U% in a spectral way if true.
LLCONTROL=.TRUE. ! to not abort (if false) and force to write fields.
LLALBEDO2=.TRUE.  ! used to compute 'SURFALBEDO.SOLNU' and 'SURFALBEDO.VEG'
LLAEROSOL=.TRUE. ! use to add aerosols and ozone climatological fields
LLSURFT=.FALSE.
INBCLIB=20       ! default: LLAEROSOL=.TRUE.
INBCGLA=3

TSRESERV1=273.15_JPRB
TSRESERV2=273.15_JPRB
TDELTA1=3.0_JPRB
TDELTA2=7.0_JPRB

NBUF901=MAX(2*NSPEC,NGPTOTG)

! ** Read namelist.

DO JFILN=1,JPFILN
  CLFILN(JFILN)=' '
ENDDO
CALL POSNAM(NULNAM,'NAMMARS')
READ(NULNAM,NAMMARS)

IF (LLOLDSWL) THEN
  ISWL1=JPOLDSWL1 
  ISWL2=JPOLDSWL2 
  ISWL3=JPOLDSWL3 
  ISWL4=JPOLDSWL4 
ELSE
  ISWL1=NGRBSWVL1 
  ISWL2=NGRBSWVL2 
  ISWL3=NGRBSWVL3 
  ISWL4=NGRBSWVL4 
ENDIF

IF (.NOT.LLAEROSOL) INBCLIB=13
IF (.NOT.LLALBEDO2) INBCLIB=11

DO JFILN=1,JPFILN
  IF(CLFILN(JFILN) /= ' ') THEN

! ** Open grib file.

    WRITE(NULOUT,FMT='(2A)')&
     & 'INITIAL GRIB DATA TO BE READ FROM FILE ',CLFILN(JFILN)  

    CALL IGRIB_OPEN_FILE(IFILE, CLFILN(JFILN),'r')
!    ILENWOR=MAX(2*NSPEC,NGPTOTG) ! array dimension in words.
    ILENWOR=NBUF901 ! array dimension in words.
    ILENBYT=ILENWOR*NINTLEN ! array dimension in bytes.

! Control dimensioning.

    WRITE(NULOUT,*) '* CPREP1 outputs:'
    WRITE(NULOUT,*) 'ILENWOR=',ILENWOR
    WRITE(NULOUT,*) 'ILENBYT=',ILENBYT
    WRITE(NULOUT,*) 'NSPEC=',NSPEC
    WRITE(NULOUT,*) 'NSPEC2=',NSPEC2
    WRITE(NULOUT,*) 'NSEFRE=',NSEFRE
    WRITE(NULOUT,*) 'NDLON=',NDLON
    WRITE(NULOUT,*) 'NDGLG=',NDGLG
    WRITE(NULOUT,*) 'NFLEVG=',NFLEVG
    WRITE(NULOUT,*) 'NLOENG='
    WRITE(NULOUT,'(25(1X,I4))')(NLOENG(JGL),JGL=NDGSAG,NDGENG)
    WRITE(NULOUT,*) 'NGPTOTG=',NGPTOTG
    WRITE(NULOUT,*) 'LLCLIM=',LLCLIM
    WRITE(NULOUT,*) 'LLOLDSWL=',LLOLDSWL
    WRITE(NULOUT,*) 'LLHUPDG=',LLHUPDG
    WRITE(NULOUT,*) 'LLCONTROL=',LLCONTROL
    WRITE(NULOUT,*) 'LLALBEDO2=',LLALBEDO2
    WRITE(NULOUT,*) 'LLAEROSOL=',LLAEROSOL
    WRITE(NULOUT,*) 'LLSWI=',LLSWI
    WRITE(NULOUT,*) 'LLLAISCAL=',LLLAISCAL
    WRITE(NULOUT,*) 'TSRESERV1=',TSRESERV1
    WRITE(NULOUT,*) 'TSRESERV2=',TSRESERV2
    WRITE(NULOUT,*) 'TDELTA1=',TDELTA1
    WRITE(NULOUT,*) 'TDELTA2=',TDELTA2
    WRITE(NULOUT,*) 'NBUF901=',NBUF901
    WRITE(NULOUT,*) 'LLMODGRIN=',LLMODGRIN

!Allocate buffers to read and stock grib
    ALLOCATE(ZGRBDAT(NBUF901), STAT=IRET)
    IF (IRET /= 0) THEN
     WRITE(NULOUT,*) '* NBUF901=',NBUF901,' IRET=',IRET
     CALL ABOR1('CPREP1: PROBLEM IN ALLOCATE ZGRBDAT buffer')
    ENDIF


! ** Read grib article.
    WRITE(NULOUT,*)'---------------------------------------------'


    IF(NCONF == 901) THEN
     CALL IGRIB_NEW_FROM_FILE(IFILE,JGRIB,IRET)
    ELSE
     CALL ABOR1('ERROR CPREP1 NCONF!')
    ENDIF

    LOOP: DO WHILE (IRET /= JPGRIB_END_OF_FILE)

! ** Decode grib article.
!    First check grib edition (temporary surf grib still grib1 ed.)
      CALL IGRIB_GET_VALUE(JGRIB,'editionNumber',IEDITION,IERR)
      WRITE(NULOUT,*)'GRIB Edition=', IEDITION

      ZGOL=8.888888E20_JPRB
      DO JB=1,ILENWOR
       ZGRBDAT(JB)=ZGOL
      ENDDO

      CALL IGRIB_GET_VALUE(JGRIB,'paramId',ISEC1(6),IERR)
      WRITE(NULOUT,*)'paramid=', ISEC1(6)
      IF (IERR /= JPGRIB_SUCCESS) THEN
       WRITE(NULOUT,*) 'paramid not found'
      ENDIF

      CALL IGRIB_GET_VALUE(JGRIB,'indicatorOfTypeOfLevel',ISEC1(7),IERR)
      IF (IERR /= JPGRIB_SUCCESS) THEN
      WRITE(NULOUT,*) "indicatorOfTypeOfLevel ==> typeOfFirstFixedSurface"
      CALL IGRIB_GET_VALUE(JGRIB,'typeOfFirstFixedSurface',ISEC1(7),IERR)
      ENDIF

      WRITE(NULOUT,*)'indicatorOfTypeOfLevel=',ISEC1(7)

      CALL IGRIB_GET_VALUE(JGRIB,'level',ISEC1(8),IERR)
      WRITE(NULOUT,*)'level'       ,ISEC1(8)
      CALL IGRIB_GET_VALUE(JGRIB,'year',ISEC1(10),IERR)
      WRITE(NULOUT,*)'year'       ,ISEC1(10)
      CALL IGRIB_GET_VALUE(JGRIB,'month',ISEC1(11),IERR)
      WRITE(NULOUT,*)'month'      ,ISEC1(11)
      CALL IGRIB_GET_VALUE(JGRIB,'day',ISEC1(12),IERR)
      WRITE(NULOUT,*)'day'      ,ISEC1(12)
      CALL IGRIB_GET_VALUE(JGRIB,'hour',ISEC1(13),IERR)
      WRITE(NULOUT,*)'hour'      ,ISEC1(13)
      CALL IGRIB_GET_VALUE(JGRIB,'minute',ISEC1(14),IERR)
      WRITE(NULOUT,*)'minute'      ,ISEC1(14)
      CALL IGRIB_GET_VALUE(JGRIB,'indicatorOfUnitOfTimeRange',ISEC1(15),IERR)
      WRITE(NULOUT,*)'indicatorOfUnitOfTimeRange',ISEC1(15)
      CALL IGRIB_GET_VALUE(JGRIB,'startStep',ISEC1(16),IERR)
      WRITE(NULOUT,*)'startStep'      ,ISEC1(16)
      ISEC1(17)=0
      ISEC1(18)=0
      ISEC1(19)=0
      ISEC1(20)=0

     CALL IGRIB_GET_VALUE(JGRIB,'values',ZGRBDAT,IERR)
     CALL IGRIB_GET_VALUE(JGRIB,'gridDefinitionTemplateNumber',ISEC2(1),IERR)
!   il faut un tableau de correspondance ou rester en grib2 
     IF (ISEC2(1)==50) THEN
      WRITE(NULOUT,*)'ISEC2(1)=50 - Spectral'
     ELSEIF (ISEC2(1)==30) THEN
      ISEC2(1)=3
     ELSE
      ISEC2(1)=4
     ENDIF

!   un traitement particulier pour une grille de gauss ou une lambert
    WRITE(NULOUT,*)'gridDefinitionTemplateNumber'     ,ISEC2(1)

    IF (ISEC2(1)==4) THEN
      CALL IGRIB_GET_VALUE(JGRIB,'resolutionAndComponentFlags',ISEC2(6),IERR)
      WRITE(NULOUT,*)'resolutionAndComponentFlags',ISEC2(6)
      IF (IERR /= JPGRIB_SUCCESS) THEN
        WRITE(NULOUT,*) "resolutionAndComponentFlags Key/value not found ==>&
                        & ISEC2(6)=0"
        ISEC2(6)=0
      ENDIF 
      CALL IGRIB_GET_VALUE(JGRIB,'Nj',ISEC2(10),IERR)
      ISEC2(10)=ISEC2(10)/2
      WRITE(NULOUT,*)'numberOfParallelsBetweenAPoleAndTheEquator',ISEC2(10)
      IF (IERR /= JPGRIB_SUCCESS) THEN
        WRITE(NULOUT,*) "numberOfParallelsBetweenAPoleAndTheEquator Key/value&
                        &  not found ==> ISEC2(10)=0" 
        ISEC2(10)=0
      ENDIF 
      CALL IGRIB_GET_VALUE(JGRIB,'pl',ISEC2(23:),IERR)
      IF (IERR /= JPGRIB_SUCCESS) THEN
        WRITE(NULOUT,*) "pl Key/value not found ==> ISEC2(23:)=0"
        ISEC2(23:)=0
      ENDIF
    ELSEIF(ISEC2(1)==50) THEN
      WRITE(NULOUT,*) " spectral tout est a faire "
        ISEC2(10)=0
    ENDIF

! Check extrema of decoded field.

  IF (ISEC2(1)==4) THEN
    CALL IGRIB_GET_VALUE(JGRIB,'max',ZMAX,IERR)
    CALL IGRIB_GET_VALUE(JGRIB,'min',ZMIN,IERR)
    CALL IGRIB_GET_VALUE(JGRIB,'average',ZMEA,IERR)
    CALL IGRIB_GET_VALUE(JGRIB,'numberOfValues',IDONR,IERR)
!    IDONR=IDONR+IEZGL*(IEZGL+2*ISEC2(10))

    WRITE(NULOUT,'(A,i4,A,i8,3(A,E16.7))')&
     & ' Parameter ',ISEC1(6),' read: size '&
     & ,IDONR,' min ',ZMIN,' max ',ZMAX,' mea ',ZMEA  

  ELSEIF(ISEC2(1)==50) THEN
    CALL IGRIB_GET_VALUE(JGRIB,'numberOfValues',IDONR,IERR)
    WRITE(NULOUT,*) 'IDONR=', IDONR
  ENDIF

    IPARAM=ISEC1(6) ! grib code of the field.
    ILEVTY=ISEC1(7) ! type of level.
    ILEVEL=ISEC1(8) ! level.
    IDATRT=ISEC2(1) ! data representation type.
    IREPRM=ISEC2(6) ! resolution flag.

    IF(IDONR > NGPTOTG) THEN
      WRITE(NULOUT,*) 'WARNING: skip parameter with incorrect resolution IPARAM=',IPARAM,&
       & ' ILEVTY=',ILEVTY,' ILEVEL=',ILEVEL,' IDATRT=',IDATRT,' IREPRM=',IREPRM
      CALL ABOR1('ABOR1 INCORRECT RESOL.')
    ENDIF

  ! ** Open FA file.

      IF(LLFIRAR) THEN

  ! First grib article read.
  ! The FA file has not yet been opened.
  ! It is necessary to open FA file here
  ! in order to get date and time from grib file.

        LLFIRAR=.FALSE.

  ! Initialize date.

        IDATEF(1)=ISEC1(10)
        IDATEF(2)=ISEC1(11)
        IDATEF(3)=ISEC1(12)
        IDATEF(4)=ISEC1(13)
        IDATEF(5)=ISEC1(14)
        IDATEF(6)=ISEC1(15)
        IDATEF(7)=ISEC1(16)
        IDATEF(8)=ISEC1(17)
        IDATEF(9)=ISEC1(18)
        IDATEF(10)=ISEC1(19)
        IDATEF(11)=ISEC1(20)

  ! Write out date.

        WRITE(NULOUT,*) 'Grib file date =',IDATEF

  ! Initialize FA parameters.

        INTIMES=1
        INBARP=NS3D*NFLEVG + NS2D

  ! Open file.

        CLFAFILN='CN90x'//CNMEXP(1:4)//'INIT' ! ARPEGE file name.
        WRITE(NULOUT,*) 'OPEN ARPEGE FILE ',CLFAFILN
        CALL FAITOU(IREP,NPOSSH,.TRUE.,CLFAFILN,'UNKNOWN',&
         & .TRUE.,.TRUE.,INTIMES,INBARP,INBARI,CNMCA)  
        CALL LFIMST(IREP,NPOSSH,.FALSE.)
        CALL FANDAR(IREP,NPOSSH,IDATEF)

  ! ** Read climatologies and check month

        IF (LLCLIM) THEN
          WRITE(UNIT=NULOUT,FMT='(A)')&
           & ' GRIDPOINT CLIMATOLOGIC INPUT DATA TO BE READ FROM FILE  '  
  ! Les controles de coherence entre fichier clim et ARP sont faits par OPENFA()
          CALL OPENFA(YDGEOMETRY,YDML_GCONF%YRRIP,10,IUCLIM,KDEB2=IDEB2,&
                     &KFIN=IFIN,KINDEX=IINDEX,KDATE=IDATEC)

  ! Initialize Constantes from YOMCLI
          CALL VAL923(.TRUE.)
        ENDIF

      ELSE

  ! This grib article is not the first read.
  ! We check now if its date is compatible
  ! with previous ones.

        IF    ((ISEC1(10)) /= IDATEF(1)&
           & .OR.ISEC1(11) /= IDATEF(2)&
           & .OR.ISEC1(12) /= IDATEF(3)&
           & .OR.ISEC1(13) /= IDATEF(4)&
           & .OR.ISEC1(14) /= IDATEF(5)) THEN  
          WRITE(NULOUT,*) 'CPREP1: WARNING DATE INCONSISTENCY!'
          WRITE(NULOUT,*) 'isec1(10:14)=',ISEC1(10:14)
          WRITE(NULOUT,*) 'idatef=',IDATEF
        ENDIF
      ENDIF

  ! Grib coding conventions values.

      IDATRTG=4 ! data representation type: grid-point fields.
      IDATRTS=50 ! data representation type: spectral fields.

      IF(IDATRT == IDATRTS) THEN

  ! The field is spectral type.

        WRITE(NULOUT,*) 'Spectral field read.'
        LLINPFS=.TRUE. ! .true. if the input field is spectral.

  ! Initialize to zero.

        DO JINDEX=1,NSPEC2
          ZFASP(JINDEX)=0.0_JPRB
        ENDDO

  ! Initialize from grib file.

        WRITE(NULOUT,*) 'ireprm=',IREPRM
        INDEX=0
        DO JM=0,NSMAX
          DO JN=JM,NSMAX
            INDEX=INDEX+1
            ZFASP(INDEX)=ZGRBDAT(INDEX)
            INDEX=INDEX+1
            IF(JM /= 0) THEN
              ZFASP(INDEX)=ZGRBDAT(INDEX)
            ELSE
              ZFASP(INDEX)=0.0_JPRB
            ENDIF
          ENDDO
        ENDDO

      ELSEIF(IDATRT == IDATRTG) THEN

  ! The field is grid-point type.

        WRITE(NULOUT,*) 'Grid-point field read.'
        LLINPFS=.FALSE.

  ! Initialize from grib file.

        ZFAGG(1:NGPTOTG)=ZGRBDAT(1:NGPTOTG)

        IF(NGPTOTG /= IDONR) THEN
           WRITE(NULOUT,*) 'WARNING: skip Grid-point parameter with incorrect resolution IPARAM=',IPARAM
          CALL IGRIB_RELEASE(JGRIB)
          CALL IGRIB_NEW_FROM_FILE(IFILE,JGRIB,IRET)
           CYCLE ! so we go to the next row of the current GRIB file
        ENDIF

      ELSE

        WRITE(NULOUT,*)'CPREP1 parameter',IPARAM,': idatrt wrong 1!...'
        WRITE(NULOUT,*) 'ISEC2=',ISEC2
        CALL ABOR1('CPREP1: ABOR1 CALLED')

      ENDIF

  ! Check consistency between grib file grid and ARPEGE grid.

     IF(ILEVTY /= 1.AND.ILEVTY /= 105.AND.ILEVTY /= 106.AND.&
    &  ILEVTY /= 109.AND.ILEVTY /= 112) THEN
        WRITE(NULOUT,*) 'LEVEL TYPE ' ,ILEVTY,'!'
        CALL ABOR1('ERROR CPREP1 LEVEL!')
      ENDIF
      IF(.NOT.LLINPFS.AND.ISEC2(10) /= NDGNH) THEN
        WRITE(NULOUT,*) 'RESOLUTION OF MODEL ',NDGNH&
         & ,' ,OF INITIAL DATA ',ISEC2(10),'!'  
        CALL ABOR1('ERROR CPREP1 NDGNH')
      ENDIF
      IF(.NOT.LLINPFS.AND.NHTYP /= 0) THEN
        DO JGL=1,NDGLG
          IF(NLOENG(JGL) /= ISEC2(22+JGL)) THEN
            WRITE(NULOUT,*) 'INCONSISTENT REDUCED GRID'
            WRITE(NULOUT,*) 'IN MODEL ',(NLOENG(JB),JB=1,NDGLG)
            WRITE(NULOUT,*) 'IN FILE  ',(ISEC2(22+JB),JB=1,NDGLG),'!'
            CALL ABOR1('ERROR CPREP1 NLOENG!')
          ENDIF
        ENDDO
      ENDIF

  ! ** Determine what to do on current parameter.

  ! Actions to do:

      LLSCAL=.FALSE. ! .true. if the current parameter has to be scaled.
      LLOUTFS=.TRUE. ! .true. if the output field has to be spectral.
      LL3DMOD=.FALSE. ! .true. if the spectral field has to be scaled by rlapin.
      LLWRIF=.TRUE. ! .true. if the field has simply to be rewritten on the FA file.

      IF(IPARAM == 138) THEN

  ! Vorticity.

        CLPREF='S'
        INIVEAU=ILEVEL
        CLSUFF='FONC.COURANT'
        LL3DMOD=.TRUE.
      ELSEIF(IPARAM == 155) THEN

  ! Divergence.

        CLPREF='S'
        INIVEAU=ILEVEL
        CLSUFF='POT.VITESSE '
        LL3DMOD=.TRUE.
      ELSEIF(IPARAM == 130) THEN

  ! Temperature.

        CLPREF='S'
        INIVEAU=ILEVEL
        CLSUFF='TEMPERATURE '
      ELSEIF(IPARAM == 133) THEN

  ! Specific humidity.
        IF (LLHUPDG) LLOUTFS=.FALSE.
        CLPREF='S'
        INIVEAU=ILEVEL
        CLSUFF='HUMI.SPECIFI'
      ELSEIF(IPARAM == 139) THEN

  ! surface temperature.

        LLSURFT=.TRUE.
        ZFSURFT(1:NGPTOTG)=ZFAGG(1:NGPTOTG)
        LLOUTFS=.FALSE.
        CLPREF='SURF'
        INIVEAU=0
        CLSUFF='TEMPERATURE '
      ELSEIF(IPARAM == 183) THEN

  ! CDST.

        LLOUTFS=.FALSE.
        CLPREF='SURF'
        INIVEAU=0
        CLSUFF='CDST        '
  ! Store DST to compute some fields at the end
        LLCDST=.TRUE.
        ZFASTL3(1:NGPTOTG)=ZFAGG(1:NGPTOTG)

      ELSEIF(IPARAM == 236) THEN

  ! STL4.

        LLOUTFS=.FALSE.
        CLPREF='SURF'
        INIVEAU=0
        CLSUFF='STL4        '
  ! Store STL4 to compute some fields at the end
        LLSTL4=.TRUE.
        ZFASTL4(1:NGPTOTG)=ZFAGG(1:NGPTOTG)

      ELSEIF(IPARAM == 170) THEN

  ! STL2

        LLOUTFS=.FALSE.
        CLPREF='SURF'
        INIVEAU=0
        CLSUFF='STL2        '
  ! Store STL2 to compute some fields at the end
        LLSTL2=.TRUE.
        ZFASTL2(1:NGPTOTG)=ZFAGG(1:NGPTOTG)

      ELSEIF(IPARAM == 172) THEN

  ! Land/Sea Mask.
        LLOUTFS=.FALSE.
        CLPREF='SURF'
        INIVEAU=0
        CLSUFF='IND.TERREMER'
  ! Store LSM to compute some fields at the end
        LLLSM=.TRUE.
        ZFALSM(1:NGPTOTG)=ZFAGG(1:NGPTOTG)

      ELSEIF(IPARAM == 174) THEN

  ! surface albedo.
        IF (LLCLIM) THEN
          LLWRIF=.FALSE.
        ENDIF
        LLOUTFS=.FALSE.
        CLPREF='SURF'
        INIVEAU=0
        CLSUFF='ALBEDO      '
      ELSEIF(IPARAM == 32) THEN

  ! snow albedo.
        LLOUTFS=.FALSE.
        CLPREF='SURF'
        INIVEAU=0
        CLSUFF='ALBEDO NEIGE'
      ELSEIF(IPARAM == 33) THEN

  ! snow density.
        LLSCAL=.TRUE.
        ZSCALE=0.001_JPRB
        LLOUTFS=.FALSE.
        CLPREF='SURF'
        INIVEAU=0
        CLSUFF='DENSIT.NEIGE'
      ELSEIF(IPARAM == 199) THEN

  ! surface prop.vegetat.
        IF (LLCLIM) THEN
          LLWRIF=.FALSE.
        ENDIF
        LLOUTFS=.FALSE.
        CLPREF='SURF'
        INIVEAU=0
        CLSUFF='PROP.VEGETAT'
      ELSEIF(IPARAM == 173) THEN

  ! surface Z0.FOIS.G.
        IF (LLCLIM) THEN
          LLWRIF=.FALSE.
        ENDIF
        LLOUTFS=.FALSE.
        CLPREF='SURF'
        INIVEAU=0
        CLSUFF='Z0.FOIS.G'

      ELSEIF(IPARAM == 152.OR.IPARAM == 134) THEN

  ! Log(surface pressure) OR surface pressure.

        CLPREF='SURF'
        INIVEAU=0
        CLSUFF='PRESSION    '
      ELSEIF(IPARAM == 129) THEN

  ! Orography.

        CLPREF='SPEC'
        INIVEAU=0
        CLSUFF='SURFGEOPOTENTI'
      ELSEIF(IPARAM == 190) THEN

  ! Orography: EW component of sub-grid scale orographic variance.

        LLOUTFS=.FALSE.
        LLWRIF=.FALSE.
        LLGRB190=.TRUE.
        CLPREF=' '
        INIVEAU=0
        CLSUFF=' '
      ELSEIF(IPARAM == 191) THEN

  ! Orography: NS component of sub-grid scale orographic variance.

        LLOUTFS=.FALSE.
        LLWRIF=.FALSE.
        LLGRB191=.TRUE.
        CLPREF=' '
        INIVEAU=0
        CLSUFF=' '
      ELSEIF(IPARAM == 192) THEN

  ! Orography: NWSE component of sub-grid scale orographic variance.

        LLOUTFS=.FALSE.
        LLWRIF=.FALSE.
        LLGRB192=.TRUE.
        CLPREF=' '
        INIVEAU=0
        CLSUFF=' '
      ELSEIF(IPARAM == 193) THEN

  ! Orography: NESW component of sub-grid scale orographic variance.

        LLOUTFS=.FALSE.
        LLWRIF=.FALSE.
        LLGRB193=.TRUE.
        CLPREF=' '
        INIVEAU=0
        CLSUFF=' '
      ELSEIF(IPARAM == 160) THEN

  ! Orography: sigma.

        LLOUTFS=.FALSE.
        CLPREF='SURF'
        INIVEAU=0
        CLSUFF='ET.GEOPOTENT'
      ELSEIF(IPARAM == 161) THEN

  ! Orography: anisotropy.

        LLOUTFS=.FALSE.
        CLPREF='SURF'
        INIVEAU=0
        CLSUFF='VAR.GEOP.ANI'
      ELSEIF(IPARAM == 162) THEN

  ! Orography: angle of principal axis.

        LLOUTFS=.FALSE.
        CLPREF='SURF'
        INIVEAU=0
        CLSUFF='VAR.GEOP.DIR'
      ELSEIF(IPARAM == 200) THEN

  ! Orography: variance (old).

        DO JGPTOTG=1,NGPTOTG
          ZFAGG(JGPTOTG)=RG*SQRT(ZFAGG(JGPTOTG))
        ENDDO
        LLOUTFS=.FALSE.
        CLPREF='SURF'
        INIVEAU=0
        CLSUFF='ET.GEOPOTENT'
      ELSEIF (IPARAM==ISWL1) THEN
        LLWRIF=.FALSE.
      ELSEIF (IPARAM==ISWL2) THEN
        LLWRIF=.FALSE.
      ELSEIF (IPARAM==ISWL3) THEN
        LLWRIF=.FALSE.
      ELSEIF (IPARAM==ISWL4) THEN
        LLWRIF=.FALSE.
      ELSEIF (IPARAM==43) THEN
        LLWRIF=.FALSE.
      ELSEIF (IPARAM==NGRBCVL) THEN
        LLWRIF=.FALSE.
      ELSEIF (IPARAM==NGRBCVH) THEN
        LLWRIF=.FALSE.
      ELSEIF (IPARAM==NGRBTVL) THEN
        LLWRIF=.FALSE.
      ELSEIF (IPARAM==NGRBTVH) THEN
        LLWRIF=.FALSE.

      ELSE

  ! Look for other surface fields.

        IOK=0

  ! soilb.

        DO JC=1,YSP_SBD%NLEVS
          DO JV=1,YSP_SBD%NUMFLDS
            CLWORKAROUND = YSP_SB%YSB(JV)%CNAME(JC)
            !IF(YSP_SB%YSB(JV)%CNAME(JC) /= ' '.AND.YSP_SB%YSB(JV)%IGRBCODE(JC) == IPARAM) THEN
            IF(CLWORKAROUND /= ' '.AND.YSP_SB%YSB(JV)%IGRBCODE(JC) == IPARAM) THEN
              IOK=1
              LLOUTFS=.FALSE.
              CLPREF=CLWORKAROUND(1:4)
              INIVEAU=0
              CLSUFF=CLWORKAROUND(5:)
            ENDIF
          ENDDO
        ENDDO

  ! snowg.
        DO JC=1,YSP_SGD%NLEVS
          DO JV=1,YSP_SGD%NUMFLDS
            CLWORKAROUND = YSP_SG%YSG(JV)%CNAME(JC)
            !IF(YSP_SG%YSG(JV)%CNAME(JC) /= ' '.AND.YSP_SG%YSG(JV)%IGRBCODE(JC) == IPARAM) THEN
            IF(CLWORKAROUND /= ' '.AND.YSP_SG%YSG(JV)%IGRBCODE(JC) == IPARAM) THEN
              IOK=1
              LLOUTFS=.FALSE.
              CLPREF=CLWORKAROUND(1:4)
              INIVEAU=0
              CLSUFF=CLWORKAROUND(5:)
            ENDIF
          ENDDO
        ENDDO
      
  ! resvr.

        DO JV=1,YSP_RRD%NUMFLDS
          IF(YSP_RR%YRR(JV)%CNAME /= ' '.AND.YSP_RR%YRR(JV)%IGRBCODE == IPARAM) THEN
            IOK=1
            LLOUTFS=.FALSE.
            CLPREF=YSP_RR%YRR(JV)%CNAME(1:4)
            INIVEAU=0
            CLSUFF=YSP_RR%YRR(JV)%CNAME(5:)
          ENDIF
        ENDDO


  ! extrp.

        DO JC=1,YSP_EPD%NLEVS
          DO JV=1,YSP_EPD%NUMFLDS
            CLWORKAROUND = YSP_EP%YEP(JV)%CNAME(JC)
            !IF(YSP_EP%YEP(JV)%CNAME(JC) /= ' '.AND. YSP_EP%YEP(JV)%IGRBCODE(JC) == IPARAM) THEN
            IF(CLWORKAROUND /= ' '.AND. YSP_EP%YEP(JV)%IGRBCODE(JC) == IPARAM) THEN
              IOK=1
              LLOUTFS=.FALSE.
              CLPREF=CLWORKAROUND(1:4)
              INIVEAU=0
              CLSUFF=CLWORKAROUND(5:)
            ENDIF
          ENDDO
        ENDDO

  ! xtrp2.

        DO JV=1,YSP_X2D%NUMFLDS
          IF(YSP_X2%YX2(JV)%CNAME /= ' ' .AND. YSP_X2%YX2(JV)%IGRBCODE == IPARAM) THEN
            IOK=1
            LLOUTFS=.FALSE.
            CLPREF=YSP_X2%YX2(JV)%CNAME(1:4)
            INIVEAU=0
            CLSUFF=YSP_X2%YX2(JV)%CNAME(5:)
          ENDIF
        ENDDO

  ! varsf.

        DO JV=1,YSD_VFD%NUMFLDS
          IF(YSD_VF%YVF(JV)%CNAME /= ' '.AND. YSD_VF%YVF(JV)%IGRBCODE == IPARAM) THEN
            IOK=1
            LLOUTFS=.FALSE.
            CLPREF=YSD_VF%YVF(JV)%CNAME(1:4)
            INIVEAU=0
            CLSUFF=YSD_VF%YVF(JV)%CNAME(5:)
          ENDIF
        ENDDO

  ! vclip.

        DO JV=1,YSD_VPD%NUMFLDS
          IF(YSD_VP%YVP(JV)%CNAME /= ' '.AND. YSD_VP%YVP(JV)%IGRBCODE == IPARAM) THEN
            IOK=1
            LLOUTFS=.FALSE.
            CLPREF=YSD_VP%YVP(JV)%CNAME(1:4)
            INIVEAU=0
            CLSUFF=YSD_VP%YVP(JV)%CNAME(5:)
          ENDIF
        ENDDO

  ! waves.

        DO JV=1,YSD_WSD%NUMFLDS
          IF(YSD_WS%YWS(JV)%CNAME /= ' '.AND. YSD_WS%YWS(JV)%IGRBCODE == IPARAM) THEN
            IOK=1
            LLOUTFS=.FALSE.
            CLPREF=YSD_WS%YWS(JV)%CNAME(1:4)
            INIVEAU=0
            CLSUFF=YSD_WS%YWS(JV)%CNAME(5:)
          ENDIF
        ENDDO

        DO JV=1,YSD_WWD%NUMFLDS
          IF(YSD_WW%YWW(JV)%CNAME /= ' '.AND. YSD_WW%YWW(JV)%IGRBCODE == IPARAM) THEN
            IOK=1
            LLOUTFS=.FALSE.
            CLPREF=YSD_WW%YWW(JV)%CNAME(1:4)
            INIVEAU=0
            CLSUFF=YSD_WW%YWW(JV)%CNAME(5:)
          ENDIF
        ENDDO

        IF(IOK == 0) THEN

  ! Unattended parameter.

          WRITE(NULOUT,*)'Parameter ',IPARAM,' unattended!...'
          CALL IGRIB_RELEASE(JGRIB)
          CALL IGRIB_NEW_FROM_FILE(IFILE,JGRIB,IRET)
          CYCLE ! so we go to the next row of the current GRIB file
        ENDIF
      ENDIF
      ! WRITE(NULOUT,*) IPARAM,SUM(ZFAGG)/REAL(NGPTOTG,KIND=JPRB),minval(ZFAGG), maxval(ZFAGG)
      IF(LL3DMOD.AND.LLINPFS) THEN

  ! Modify 3D spectral fields spectrum.

        WRITE(NULOUT,*) 'Modify 3D spectral fields spectrum.'
        DO JNM=1,NSPEC2
          ZFASP(JNM)=ZFASP(JNM)*RLAPIN(NVALUE(JNM))
        ENDDO
      ENDIF

  ! ** Spectral direct/inverse transform.

      IF(LLINPFS.AND..NOT.LLOUTFS) THEN

  ! The input field is spectral and the output one
  ! has to be grid-point.

        WRITE(NULOUT,*) 'Call speree.'
        CALL SPEREE(YDGEOMETRY,1,1,ZFASP,ZFAGG)
      ELSEIF(.NOT.LLINPFS.AND.LLOUTFS) THEN

  ! The input field is grid-point and the output one
  ! has to be spectral.
        WRITE(NULOUT,*) 'Call reespe.'
        CALL REESPE(YDGEOMETRY,1,1,ZFASP,ZFAGG)
      ENDIF

  ! ** Odd scalings.

      IF(IPARAM == NGRBSDSL  ) THEN
        LLSCAL=.TRUE.
        ZSCALE=1000.0_JPRB
      ELSEIF(IPARAM == NGRBSR  ) THEN
        LLSCAL=.TRUE.
        ZSCALE=RG
      ELSEIF(IPARAM == NGRBSRC ) THEN
        LLSCAL=.TRUE.
        ZSCALE=1000.0_JPRB
      ELSEIF(IPARAM == NGRBSDOR) THEN
        LLSCAL=.TRUE.
        ZSCALE=RG
      ENDIF
      IF(LLSCAL.AND..NOT.LLOUTFS) THEN

  ! A scaling is required and the output field
  ! is grid-point type.

        WRITE(NULOUT,*) 'Scaling: ',ZSCALE
        ZFAGG(1:NGPTOTG)=ZFAGG(1:NGPTOTG)*ZSCALE
      ELSEIF(LLSCAL) THEN

  ! A scaling is required and the output field
  ! is spectral type.

        CALL ABOR1('CPREP1 scaling on spectral!...')
      ENDIF

      IF(LLWRIF) THEN

  ! ** Write field on FA file.

        IF(LLOUTFS) THEN

  ! The output field has to be spectral.
  ! One computes and prints
  ! its minima/maxima in grid-point space.

  ! One copies zfasp on zfasps because
  ! the input spectral data will be corrupted
  ! by speree.

          WRITE(NULOUT,*) 'Write spectral field'
          ZFASPS(1:NSPEC2)=ZFASP(1:NSPEC2)
          CALL SPEREE(YDGEOMETRY,1,1,ZFASPS,ZFAGG)
          CLLIBELLE='Extrema in grid-point space:'
          CALL CPREP1_DIAG(ZFAGG,CLLIBELLE)

  ! Write spectral field.

          CALL FAVEUR(IREP,NPOSSH,INGRIB,INBPDG,INBCSP,ISTRON,IPUILA,IDMOPL)
  !       data should be reordered only if packing type is not -1 and not 3
  !       Notice : If INGRIB=-1 or 3 the field will be a bit larger than if
  !       INGRIB=0, 1 or 2 because then there is no reordering : so the
  !       column 0 is doubled.
          IF (INGRIB==-1 .OR. INGRIB==3) THEN
            INGRIG=-1
          ELSE
            INGRIG=0
          ENDIF
          IF (INGRIG==0) THEN
            LLFILE_TO_MODEL=.FALSE.
            CALL SPREORD(YDGEOMETRY%YRDIM,YDGEOMETRY%YREDIM,YDGEOMETRY%YRELAP,1,ZFPDAT,ZFASP,LLFILE_TO_MODEL)
          ENDIF
          IF (TRIM(CLPREF)//CLSUFF == 'SPECSURFGEOPOTENTI') THEN
  !         Do not pack surface geopotential
            CALL FAGOTE(IREP,NPOSSH,INGRIG,INBPDG,INBCSP,ISTRON,IPUILA,IDMOPL)
          ENDIF
          IF (INGRIG==0) THEN
            CALL FAIENC(IREP,NPOSSH,CLPREF,INIVEAU,CLSUFF,ZFPDAT,LLOUTFS)
          ELSE
            CALL FAIENC(IREP,NPOSSH,CLPREF,INIVEAU,CLSUFF,ZFASP,LLOUTFS)
          ENDIF
          IF (TRIM(CLPREF)//CLSUFF == 'SPECSURFGEOPOTENTI') THEN
  !         Reset packing after the treatment of surface geopotential
            CALL FAGOTE(IREP,NPOSSH,INGRIB,INBPDG,INBCSP,ISTRON,IPUILA,IDMOPL)
          ENDIF
          WRITE(NULOUT,*) 'Write spectral field param=',&
           & IPARAM,' ',CLPREF(1:4),INIVEAU,CLSUFF(1:16)  
        ELSE

  ! Write grid-point field.

          WRITE(NULOUT,*) 'Write grid-point field param=',&
           & IPARAM,' ',CLPREF(1:4),INIVEAU,CLSUFF(1:16)  
          IF (TRIM(CLPREF)//CLSUFF == 'SURFZ0.FOIS.G' .OR.&
             & TRIM(CLPREF)//CLSUFF == 'SURFZ0REL.FOIS.G' .OR.&
             & TRIM(CLPREF)//CLSUFF == 'SURFGZ0.THERM') THEN  
  !         Do not pack roughness lengths
            CALL FAVEUR(IREP,NPOSSH,INGRIB,INBPDG,INBCSP,ISTRON,IPUILA,IDMOPL)
            INGRIG=0
            CALL FAGOTE(IREP,NPOSSH,INGRIG,INBPDG,INBCSP,ISTRON,IPUILA,IDMOPL)
          ENDIF
          CALL FAIENC(IREP,NPOSSH,CLPREF,INIVEAU,CLSUFF,ZFAGG,LLOUTFS)
          IF (TRIM(CLPREF)//CLSUFF == 'SURFZ0.FOIS.G' .OR.&
             & TRIM(CLPREF)//CLSUFF == 'SURFZ0REL.FOIS.G' .OR.&
             & TRIM(CLPREF)//CLSUFF == 'SURFGZ0.THERM') THEN  
  !         Reset packing after the treatment of roughness lengths
            CALL FAGOTE(IREP,NPOSSH,INGRIB,INBPDG,INBCSP,ISTRON,IPUILA,IDMOPL)
          ENDIF
          IF (IPARAM == 173.AND..NOT.LLCLIM) THEN
            ZFZ0G(1:NGPTOTG)=ZFAGG(1:NGPTOTG)
            LLZ0=.TRUE.
          ENDIF
        ENDIF

  !The field has not to be written on FA file.
      ELSEIF (IPARAM==ISWL1) THEN
        ZFSWL1(1:NGPTOTG)=ZFAGG(1:NGPTOTG)
        LLSWL1=.TRUE.
      ELSEIF (IPARAM==ISWL2) THEN
        ZFSWL2(1:NGPTOTG)=ZFAGG(1:NGPTOTG)
        LLSWL2=.TRUE.
      ELSEIF (IPARAM==ISWL3) THEN
        ZFSWL3(1:NGPTOTG)=ZFAGG(1:NGPTOTG)
        LLSWL3=.TRUE.
      ELSEIF (IPARAM==ISWL4) THEN
        ZFSWL4(1:NGPTOTG)=ZFAGG(1:NGPTOTG)
        LLSWL4=.TRUE.
      ELSEIF (IPARAM==43) THEN
        ZFSLT(1:NGPTOTG)=ZFAGG(1:NGPTOTG)
        LLSLT=.TRUE.
      ELSEIF (IPARAM==NGRBCVL) THEN
        ZFCVL(1:NGPTOTG)=ZFAGG(1:NGPTOTG)
        LLCVL=.TRUE.
      ELSEIF (IPARAM==NGRBCVH) THEN
        ZFCVH(1:NGPTOTG)=ZFAGG(1:NGPTOTG)
        LLCVH=.TRUE.
      ELSEIF (IPARAM==NGRBTVL) THEN
        ZFTVL(1:NGPTOTG)=ZFAGG(1:NGPTOTG)
        LLTVL=.TRUE.
      ELSEIF (IPARAM==NGRBTVH) THEN
        ZFTVH(1:NGPTOTG)=ZFAGG(1:NGPTOTG)
        LLTVH=.TRUE.
      ELSEIF (IPARAM /= 173.AND.IPARAM /= 174.AND.IPARAM /= 199) THEN
  ! In this case,  surfz0.fois.g, surfalbedo and surfprop.vegetat are read
  ! from climatological file and not stored. 

  ! The field has not to be written on FA file.
  ! It has to be stored in an intermediate array.

        IGRBSN=IGRBSN+1
        IGRBS(IGRBSN)=IPARAM
        INDEX=0
  !      CLLIBELLE='Field '// IGRBSN //' stored in an intermediate array.'
        WRITE(CLLIBELLE,FMT=*) 'Field ',IGRBSN,' stored in an intermediate array.'
        CALL CPREP1_DIAG(ZFAGG,CLLIBELLE)
        ZFAGGA(1:NGPTOTG,IGRBSN)=ZFAGG(1:NGPTOTG)
      ENDIF

! Go to next article in the GRIB file.

 CALL IGRIB_RELEASE(JGRIB)
 CALL IGRIB_NEW_FROM_FILE(IFILE,JGRIB,IRET)

 ENDDO LOOP

! Deallocate
IF (ALLOCATED(ZGRBDAT)) DEALLOCATE(ZGRBDAT)

! Close grib file.

    CALL IGRIB_CLOSE_FILE(IFILE)
  ENDIF

ENDDO   ! end of work on each GRIB file

IF (.NOT.LLLSM) THEN
  IF (.NOT. LLCONTROL) THEN
    WRITE(NULOUT,*) 'WARNING: NO LAND/SEA MASK IN GRIB FILE!!!'
    WRITE(NULOUT,*) 'ALL POINTS ARE SEA POINTS!!!'
    ZFALSM(1:NGPTOTG)=0._JPRB
  ELSE
    CALL ABOR1('CPREP1: NO LAND/SEA MASK IN GRIB FILE')
  ENDIF
ENDIF

IF (LLCLIM) THEN ! ** Working simply whith climatological fields
  WRITE(NULOUT,*) 'COPY CLIMATOLOGICAL FIELDS FROM FILE'
  INIVEAU=0
  DO IIND=1,INBCLIB
    CLPREF=CLLIB(IIND)(1:4)
    CLSUFF=CLLIB(IIND)(5:20)
    CALL FACOCH(IREP,IUCLIM,NPOSSH,CLPREF,INIVEAU,CLSUFF)
    WRITE(NULOUT,*) 'Write grid-point field '&
     & ,CLPREF(1:4),INIVEAU,CLSUFF(1:16)  
  ENDDO
ELSE
  LLOUTFS=.FALSE.
  INIVEAU=0
  DO IIND=5,12
    CLPREF=YL_FIELDCST(IIND)%CLHAMP(1:4)
    CLSUFF=YL_FIELDCST(IIND)%CLHAMP(5:20)
    IF (YL_FIELDCST(IIND)%LLISLSM) THEN ! link whith Land/Sea Mask
      DO JGPTOTG=1,NGPTOTG
        IF (ZFALSM(JGPTOTG) > 0.5) THEN
          ZFAGG(JGPTOTG)=YL_FIELDCST(IIND)%ZVALTERRE
        ELSE
          ZFAGG(JGPTOTG)=YL_FIELDCST(IIND)%ZVALMER
        ENDIF 
      ENDDO
      ZMIN=MIN(YL_FIELDCST(IIND)%ZVALTERRE,YL_FIELDCST(IIND)%ZVALMER)
      ZMAX=MAX(YL_FIELDCST(IIND)%ZVALTERRE,YL_FIELDCST(IIND)%ZVALMER)
    ELSE
      ZVAL=YL_FIELDCST(IIND)%ZVALTERRE
      ZFAGG(1:NGPTOTG)=ZVAL
      ZMIN=ZVAL
      ZMAX=ZVAL
    ENDIF
    WRITE(NULOUT,*)'---------------------------------------------'
    WRITE(NULOUT,*) YL_FIELDCST(IIND)%CLHAMP ,' field'
    WRITE(NULOUT,*) 'Min=',ZMIN,' Max=',ZMAX
    CALL FAIENC(IREP,NPOSSH,CLPREF,INIVEAU,CLSUFF,ZFAGG,LLOUTFS)
    WRITE(NULOUT,*) 'Write grid-point field '&
     & ,CLPREF(1:4),INIVEAU,CLSUFF(1:16)  

!store ARGILE, SABLE, SURFEPAIS.SOL for ACSOLW (need to compute RESERV.EAU)
    IF (YL_FIELDCST(IIND)%CLHAMP=='SURFPROP.ARGILE     ') THEN
      ZFPARG(1:NGPTOTG)=ZFAGG(1:NGPTOTG)
    ELSEIF (YL_FIELDCST(IIND)%CLHAMP=='SURFPROP.SABLE      ') THEN
      ZFPSAB(1:NGPTOTG)=ZFAGG(1:NGPTOTG)
    ELSEIF (YL_FIELDCST(IIND)%CLHAMP=='SURFEPAIS.SOL       ') THEN
      ZFPD2(1:NGPTOTG)=ZFAGG(1:NGPTOTG)
    ELSEIF (YL_FIELDCST(IIND)%CLHAMP=='SURFRESI.STO.MIN    ') THEN
      ZFPRSMIN(1:NGPTOTG)=ZFAGG(1:NGPTOTG)
    ELSEIF (YL_FIELDCST(IIND)%CLHAMP=='SURFIND.FOLIAIRE    ') THEN
      ZFPLAI(1:NGPTOTG)=ZFAGG(1:NGPTOTG)
    ELSEIF (YL_FIELDCST(IIND)%CLHAMP=='SURFPROP.VEGETAT    ') THEN
      ZFPVEG(1:NGPTOTG)=ZFAGG(1:NGPTOTG)
    ENDIF

  ENDDO
ENDIF  ! OF (LLCLIM)

! Les champs SURFRESERV.GLACE,PROFRESERV.GLACE,SURFRES.EVAPOTRA et
!SURFRESERV.INTER n'existent pas en CLIM et sont codes independamment de LLCLIM.

LLOUTFS=.FALSE.
INIVEAU=0
DO IIND=INBCGLA,4
  CLPREF=YL_FIELDCST(IIND)%CLHAMP(1:4)
  CLSUFF=YL_FIELDCST(IIND)%CLHAMP(5:20)
  ZVAL=YL_FIELDCST(IIND)%ZVALTERRE
  ZFAGG(1:NGPTOTG)=ZVAL
  ZMIN=ZVAL
  ZMAX=ZVAL
  WRITE(NULOUT,*)'---------------------------------------------'
  WRITE(NULOUT,*) YL_FIELDCST(IIND)%CLHAMP ,' field'
  WRITE(NULOUT,*) 'Min=',ZMIN,' Max=',ZMAX
  CALL FAIENC(IREP,NPOSSH,CLPREF,INIVEAU,CLSUFF,ZFAGG,LLOUTFS)
  WRITE(NULOUT,*) 'Write grid-point field '&
   & ,CLPREF(1:4),INIVEAU,CLSUFF(1:16)  
ENDDO

CALL SURF_INQ(YSURF,KNVTYPES=IVTYPES)
ALLOCATE(ZRVCOV(0:IVTYPES))
CALL SURF_INQ(YSURF,PRVCOV=ZRVCOV)

IF (LLCVL.AND.LLCVH.AND.LLTVL.AND.LLTVH) THEN
  ALLOCATE(ZRVROOTSA(4,0:IVTYPES))
  ALLOCATE(ZRVLAI(0:IVTYPES))
  CALL SURF_INQ(YSURF,PRVROOTSA=ZRVROOTSA)
  CALL SURF_INQ(YSURF,PRVLAI=ZRVLAI)
!    CALL SURF_INQ(YREPHY%YSURF,PRVRSMIN=ZRVRSMIN)  ! does not work yet
 
  ZRVRSMIN(1)=180._JPRB    ! Crops, Mixed Farming
  ZRVRSMIN(2)=110._JPRB    ! Short Grass
  ZRVRSMIN(3)=500._JPRB    ! Evergreen Needleleaf Trees
  ZRVRSMIN(4)=500._JPRB    ! Deciduous Needleleaf Trees
  ZRVRSMIN(5)=175._JPRB    ! Deciduous Broadleaf Trees
  ZRVRSMIN(6)=240._JPRB    ! Evergreen Broadleaf Trees
  ZRVRSMIN(7)=100._JPRB    ! Tall Grass
  ZRVRSMIN(8)=250._JPRB    ! Desert
  ZRVRSMIN(9)=80._JPRB     ! Tundra
  ZRVRSMIN(10)=180._JPRB   ! Irrigated Crops
  ZRVRSMIN(11)=150._JPRB   ! Semidesert
  ZRVRSMIN(12)=0.0_JPRB      ! Ice Caps and Glaciers
  ZRVRSMIN(13)=240._JPRB   ! Bogs and Marshes
  ZRVRSMIN(14)=0.0_JPRB      ! Inland Water
  ZRVRSMIN(15)=0.0_JPRB      ! Ocean
  ZRVRSMIN(16)=225._JPRB   ! Evergreen Shrubs
  ZRVRSMIN(17)=225._JPRB   ! Deciduous Shrubs
  ZRVRSMIN(18)=250._JPRB   ! Mixed Forest/woodland
  ZRVRSMIN(19)=175._JPRB   ! Interrupted Forest
  ZRVRSMIN(20)=150._JPRB   ! Water and Land Mixtures
  ZRVRSMIN(0)=ZRVRSMIN(8)
ENDIF

! ** Case of SURFGZ0.THERM and SURFPROP.VEGETAT
IF (.NOT.LLCLIM) THEN
  LLOUTFS=.FALSE.
  INIVEAU=0
  CLPREF='SURF'
  IF (LLZ0) THEN
    CLSUFF='GZ0.THERM    '
    DO JGPTOTG=1,NGPTOTG
      IF (ZFALSM(JGPTOTG) > 0.5) THEN
        ZFAGG(JGPTOTG)=ZFZ0G(JGPTOTG)/10.0_JPRB
      ELSE
        ZFAGG(JGPTOTG)=ZFZ0G(JGPTOTG)
      ENDIF 
    ENDDO
    WRITE(NULOUT,*)'---------------------------------------------'
    CLLIBELLE='SURFGZ0.THERM field'
    CALL CPREP1_DIAG(ZFAGG,CLLIBELLE)
!    Do not pack roughness lengths
    CALL FAVEUR(IREP,NPOSSH,INGRIB,INBPDG,INBCSP,ISTRON,IPUILA,IDMOPL)
    INGRIG=0
    CALL FAGOTE(IREP,NPOSSH,INGRIG,INBPDG,INBCSP,ISTRON,IPUILA,IDMOPL)
    CALL FAIENC(IREP,NPOSSH,CLPREF,INIVEAU,CLSUFF,ZFAGG,LLOUTFS)
!    Reset packing after the treatment of roughness lengths
    CALL FAGOTE(IREP,NPOSSH,INGRIB,INBPDG,INBCSP,ISTRON,IPUILA,IDMOPL)
    WRITE(NULOUT,*) 'Write grid-point field '&
     & ,CLPREF(1:4),INIVEAU,CLSUFF(1:16)  
  ELSE
    WRITE(NULOUT,*)&
     & 'SURFGZ0.THERM CANNOT BE COMPUTED BECAUSE SURFZ0.FOIS.G IS NOT'  
    IF (LLCONTROL) CALL ABOR1('CPREP1 : ABOR1 CALLED')
  ENDIF

  IF (LLCVL.AND.LLCVH.AND.LLTVL.AND.LLTVH) THEN
    CLSUFF='PROP.VEGETAT ' 
    DO JGPTOTG=1,NGPTOTG
      ZFAGG(JGPTOTG)=(ZFCVH(JGPTOTG)*ZRVCOV(NINT(ZFTVH(JGPTOTG)))) +&
       & (ZFCVL(JGPTOTG)*ZRVCOV(NINT(ZFTVL(JGPTOTG))))  
    ENDDO

    WRITE(NULOUT,*)'---------------------------------------------'
    CLLIBELLE='SURFPROP.VEGETAT field'
    CALL CPREP1_DIAG(ZFAGG,CLLIBELLE)
    CALL FAIENC(IREP,NPOSSH,CLPREF,INIVEAU,CLSUFF,ZFAGG,LLOUTFS)
    WRITE(NULOUT,*) 'Write grid-point field '&
     & ,CLPREF(1:4),INIVEAU,CLSUFF(1:16)  
  ELSEIF (.NOT.LLOLDSWL) THEN
    WRITE(NULOUT,*) 'CPREP1: WARNING&
     & ! ','SURFPROP.VEGETAT NOT COMPUTED WHITH NEW METHOD (MISSING PARAMETERS)'  
  ENDIF

ENDIF

IF (LLMODGRIN .AND. LLLSM) THEN
  DO JGPTOTG=1,NGPTOTG
    IF (ZFALSM(JGPTOTG) <= 0.5) THEN
      ZFASTL2(JGPTOTG)=ZFSURFT(JGPTOTG)
      ZFASTL3(JGPTOTG)=ZFSURFT(JGPTOTG)
      ZFASTL4(JGPTOTG)=ZFSURFT(JGPTOTG)
    ENDIF
  ENDDO
ENDIF

! ** Set ProftempeRATURE
IF ((LLCDST.AND.LLSTL2.AND.LLSWI).OR.(LLCDST.AND.LLSTL4.AND..NOT.LLSWI)) THEN ! PROFTEMPERATURE can be compute

 IF (LLSWI) THEN
   ZFAGG=(ZFASTL3 + ZFASTL2) * 0.5_JPRB
 ELSE
   ZFAGG=(ZFASTL3 + ZFASTL4) * 0.5_JPRB
 ENDIF

  ZFPROFT(1:NGPTOTG)=ZFAGG(1:NGPTOTG)
  LLOUTFS=.FALSE.
  CLPREF='PROF'
  INIVEAU=0
  CLSUFF='TEMPERATURE '
  CALL FAIENC(IREP,NPOSSH,CLPREF,INIVEAU,CLSUFF,ZFAGG,LLOUTFS)
  WRITE(NULOUT,*)'---------------------------------------------'
  CLLIBELLE='PROFTEMPERATURE field'
  CALL CPREP1_DIAG(ZFAGG,CLLIBELLE)

  WRITE(NULOUT,*) 'Write grid-point field '&
   & ,CLPREF(1:4),INIVEAU,CLSUFF(1:16)  
ELSE
  WRITE(NULOUT,*)&
   & 'PROFTEMPERATURE CANNOT BE COMPUTED BECAUSE CDST OR STL2 OR STL4 IS NOT'  
  IF (LLCONTROL) CALL ABOR1('CPREP1 : ABOR1 CALLED')
ENDIF

! Set ZFPIVEG
DO JGPTOTG=1,NGPTOTG
  IF (ZFALSM(JGPTOTG) > 0.5) THEN
    ZFPIVEG(JGPTOTG)=JPNTVTER
  ELSE
    ZFPIVEG(JGPTOTG)=JPNTVMER
  ENDIF
ENDDO

! ** Set SURFRESERV.EAU and PROFRESERV.EAU
IF (.NOT. LLOLDSWL) THEN
  ZPSOL=1.0_JPRB
ELSE
  ZPSOL=PPPSWL1
ENDIF

!formule utilisee (LLSWI=.FALSE.) pour estimer les reservoirs de glace a partir des reservoirs d'eau:
! T0=TSRESERV1 ou TSRESERV2 : temperature de base proche de 272
! delta=TDELTA1 ou TDELTA2  : intervalle autour de la temperature de base
! p= proportion de glace prise dans les reservoirs d'eau
! Tlue= temperature dans le sol profond ou superficiel
! p= 3 * (Teta*Teta) + 2 * (Teta*Teta*Teta)
! avec Teta= ((Tlue - T0)/(2 * delta)) - 0.5

!Compute SURFRESERV.EAU and SURFRESERV.GLACE
IF (LLSWL1) THEN

  IF (.NOT.LLSURFT) THEN
    WRITE(NULOUT,*)&
    & 'SURFRESERV.GLACE CANNOT BE COMPUTED BECAUSE SURFTEMPERATURE IS NOT READ'
    IF (LLCONTROL) CALL ABOR1('CPREP1 : ABOR1 CALLED')
  ENDIF

  IF (LLCLIM) THEN
! Read ARGILE, SABLE, SURFEPAIS.SOL from climatological FILE if LLCLIM
    CALL FACILE(IREP,IUCLIM,'SURF',0,'PROP.ARGILE',ZFPARG,.FALSE.)
    CALL FACILE(IREP,IUCLIM,'SURF',0,'PROP.SABLE',ZFPSAB,.FALSE.)
    CALL FACILE(IREP,IUCLIM,'SURF',0,'EPAIS.SOL',ZFPD2,.FALSE.)
    CALL FACILE(IREP,IUCLIM,'SURF',0,'IND.FOLIAIRE',ZFPLAI,.FALSE.)
    CALL FACILE(IREP,IUCLIM,'SURF',0,'RESI.STO.MIN',ZFPRSMIN,.FALSE.)
    CALL FACILE(IREP,IUCLIM,'SURF',0,'PROP.VEGETAT',ZFPVEG,.FALSE.)
  ENDIF

  CALL ACSOLW(YDPHY1,1,NGPTOTG,NGPTOTG,ZFPARG,ZFPD2,ZFALSM,ZFPIVEG,ZFPSAB,&
   & .FALSE.,ZFPWFC,ZFPWPMX,ZFPWSAT,ZFPWSMX,ZFPWWILT)

  ZSCALE=1000.0_JPRB
  IF (.NOT.LLSLT) THEN 
    WRITE(NULOUT,*) 'CPREP1 WARNING: SURFRESERV.EAU computed without LST !!'  
  ENDIF

  DO JGPTOTG=1,NGPTOTG
    IF (ZFALSM(JGPTOTG) > 0.5) THEN
     IF (LLSWI) THEN
      IF (LLSLT) THEN
        ZPWWILT_CEP=ZWP_WILT(NINT(ZFSLT(JGPTOTG)))
        ZPWCAP_CEP=ZWP_CAP(NINT(ZFSLT(JGPTOTG)))
      ELSE
        ZPWWILT_CEP=PPWWILT_CEP
        ZPWCAP_CEP=PPWCAP_CEP
      ENDIF

      IF (ZFSURFT(JGPTOTG) > (TSRESERV1+TDELTA1)) THEN
        ZPGLACE=0.0
      ELSEIF (ZFSURFT(JGPTOTG) < (TSRESERV1-TDELTA1)) THEN
        ZPGLACE=1.0
      ELSE
        ZDELTA = RPI*(ZFSURFT(JGPTOTG) - 0.5*(TSRESERV1+TDELTA1) - 0.5*(TSRESERV1-TDELTA1))/(2.0*TDELTA1)
        ZPGLACE = 0.5 * (1.0 - SIN(ZDELTA))
      ENDIF

      ZPLIQUID_LAYERS(1) = (1.0 - ZPGLACE) * ZFSWL1(JGPTOTG)
      ZPGLACE =  ZPGLACE * ZFSWL1(JGPTOTG)/ZPSOL
      ZPGLACE =  ZPGLACE * (ZSCALE*RD1)
      ZFGLACE(JGPTOTG)=MIN(ZPGLACE,150.0_JPRB)

      !CALCUL DU SWI CEP
      ZSWI_CEP= ((ZPLIQUID_LAYERS(1)/ZPSOL))/(ZPWCAP_CEP)

      !CALCUL Wp ISBA
      IF (ZSWI_CEP > 1.0_JPRB) THEN
        ZGOL = ZFPWFC(JGPTOTG)
      ELSE
        ZGOL = ZSWI_CEP * ZFPWFC(JGPTOTG)
      ENDIF

      !Calcul Wp ISBA en m
      ZGOL=ZGOL*(ZSCALE*RD1)
      ZFAGG(JGPTOTG)=ZGOL
     ELSE !if not LLSWI, old computation
      IF (LLSLT) THEN
        ZWSATSLT =ZWP_SAT(NINT(ZFSLT(JGPTOTG)))
      ELSE
        ZWSATSLT=PPWSAT_CEP
      ENDIF
      ZGOL=((ZFSWL1(JGPTOTG)/ZPSOL)*ZSCALE*RD1*ZFPWSAT(JGPTOTG))/&
       & ZWSATSLT  
        IF (.NOT.LLSURFT) THEN
          WRITE(NULOUT,*)&
           & 'SURFRESERV.GLACE CANNOT BE COMPUTED BECAUSE SURFTEMPERATURE IS NOT READ'  
          IF (LLCONTROL) CALL ABOR1('CPREP1 : ABOR1 CALLED')
        ENDIF

        IF (ZFSURFT(JGPTOTG) > (TSRESERV1+TDELTA1)) THEN
          ZPGLACE=0.0
        ELSEIF (ZFSURFT(JGPTOTG) < (TSRESERV1-TDELTA1)) THEN
          ZPGLACE=1.0
        ELSE
          ZDELTA=((ZFSURFT(JGPTOTG)-TSRESERV1)/(2.0*TDELTA1)) - 0.5
          ZPGLACE=MIN(((3.0_JPRB*ZDELTA*ZDELTA)&
           & + (2.0_JPRB*ZDELTA*ZDELTA*ZDELTA)),1.0_JPRB)  
        ENDIF
        ZFGLACE(JGPTOTG)=ZGOL*ZPGLACE
        ZFAGG(JGPTOTG)=ZGOL-ZFGLACE(JGPTOTG)
     ENDIF
    ELSE
      ZFAGG(JGPTOTG)=RD1*GCONV
      ZFGLACE(JGPTOTG)=0.0_JPRB
    ENDIF
  ENDDO

  CALL FAIENC(IREP,NPOSSH,'SURF',0,'RESERV.EAU  ',ZFAGG,.FALSE.)
  WRITE(NULOUT,*)'---------------------------------------------'
  CLLIBELLE='SURFRESERV.EAU field'
  CALL CPREP1_DIAG(ZFAGG,CLLIBELLE)
  WRITE(NULOUT,*) 'Write grid-point field '
  IF (LLSURFT) THEN
    CALL FAIENC(IREP,NPOSSH,'SURF',0,'RESERV.GLACE',ZFGLACE,.FALSE.)
    WRITE(NULOUT,*)'---------------------------------------------'
    CLLIBELLE='SURFRESERV.GLACE field'
    CALL CPREP1_DIAG(ZFGLACE,CLLIBELLE)
    WRITE(NULOUT,*) 'Write grid-point field '
  ENDIF
ELSE
  WRITE(NULOUT,*)&
   & 'SURFRESERV.EAU NOT COMPUTED BECAUSE MISSING PARAMETERS IN GRIB FILE'  
  IF (LLCONTROL) CALL ABOR1('CPREP1 : ABOR1 CALLED')
ENDIF

!Compute PROFRESERV.EAU and PROFRESERV.GLACE
IF (LLSWL1.AND.LLSWL2.AND.LLSWL3.AND.LLSWL4) THEN
  IF (.NOT.LLSLT) THEN 
    WRITE(NULOUT,*) 'CPREP1 WARNING: PROFRESERV.EAU computed without LST !!'  
  ENDIF
IF (((.NOT.LLCDST.OR..NOT.LLSTL4).AND.(.NOT.LLSWI)).OR.&
   &((.NOT.LLCDST.OR..NOT.LLSTL2).AND.(LLSWI))) THEN
  WRITE(NULOUT,*)&
  & 'PROFRESERV.GLACE CANNOT BE COMPUTED BECAUSE CDST OR STL4 OR STL2IS NOT READ'
  IF (LLCONTROL) CALL ABOR1('CPREP1 : ABOR1 CALLED')
ENDIF

  DO JGPTOTG=1,NGPTOTG
    IF (ZFALSM(JGPTOTG) > 0.5) THEN
     IF (LLSWI) THEN
       IF (LLSLT) THEN
         ZPWWILT_CEP=ZWP_WILT(NINT(ZFSLT(JGPTOTG)))
         ZPWCAP_CEP=ZWP_CAP(NINT(ZFSLT(JGPTOTG)))
       ELSE
         ZPWWILT_CEP=PPWWILT_CEP
         ZPWCAP_CEP=PPWCAP_CEP
       ENDIF

        ZPSOILT_LAYERS(1) = ZFSURFT(JGPTOTG)
        ZPSOILT_LAYERS(2) = ZFASTL2(JGPTOTG)
        ZPSOILT_LAYERS(3) = ZFASTL3(JGPTOTG)
        ZPSOILT_LAYERS(4) = ZFASTL4(JGPTOTG)

        DO ISOIL = 1,4
          IF (ZPSOILT_LAYERS(ISOIL) > (TSRESERV2+TDELTA2)) THEN
           ZPGLACE_LAYERS(ISOIL) = 0.0
          ELSEIF (ZPSOILT_LAYERS(ISOIL) < (TSRESERV2-TDELTA2)) THEN
           ZPGLACE_LAYERS(ISOIL)=1.0
          ELSE
            ZDELTA = RPI*(ZPSOILT_LAYERS(ISOIL) - 0.5*(TSRESERV2+TDELTA2) - 0.5*(TSRESERV2-TDELTA2))/(2.0*TDELTA2)
            ZPGLACE_LAYERS(ISOIL) = 0.5 * (1.0 - SIN(ZDELTA))
          ENDIF

         ENDDO

!       liquid water on each layer
         ZPLIQUID_LAYERS(1) = (1.0 - ZPGLACE_LAYERS(1)) * ZFSWL1(JGPTOTG)
         ZPLIQUID_LAYERS(2) = (1.0 - ZPGLACE_LAYERS(2)) * ZFSWL2(JGPTOTG)
         ZPLIQUID_LAYERS(3) = (1.0 - ZPGLACE_LAYERS(3)) * ZFSWL3(JGPTOTG)
         ZPLIQUID_LAYERS(4) = (1.0 - ZPGLACE_LAYERS(4)) * ZFSWL4(JGPTOTG)


         ZPGLACE = ZPGLACE_LAYERS(1) * ZFSWL1(JGPTOTG) * PPZD(1) + ZPGLACE_LAYERS(2) * ZFSWL2(JGPTOTG) * PPZD(2) +&
                 & ZPGLACE_LAYERS(3) * ZFSWL3(JGPTOTG) * PPZD(3) + ZPGLACE_LAYERS(4) * ZFSWL4(JGPTOTG) * PPZD(4)

         ZPGLACE = ZPGLACE / (PPZD(1) + PPZD(2) + PPZD(3) + PPZD(4))
         ZPGLACE = ZPGLACE * ZSCALE * ZFPD2(JGPTOTG)

!    Ice water over the four layers
         ZFGLACE(JGPTOTG)=MIN(ZPGLACE,150.0_JPRB)

         ZTETA_MEAN_H = MAX(ZPLIQUID_LAYERS(1)/ZPSOL,ZPWWILT_CEP)*ZRVROOTSA(1,NINT(ZFTVH(JGPTOTG)))&
 &            + MAX(ZPLIQUID_LAYERS(2)/ZPSOL,ZPWWILT_CEP)*ZRVROOTSA(2,NINT(ZFTVH(JGPTOTG)))&
 &            + MAX(ZPLIQUID_LAYERS(3)/ZPSOL,ZPWWILT_CEP)*ZRVROOTSA(3,NINT(ZFTVH(JGPTOTG)))&
 &            + MAX(ZPLIQUID_LAYERS(4)/ZPSOL,ZPWWILT_CEP)*ZRVROOTSA(4,NINT(ZFTVH(JGPTOTG))) 

         IF (ZTETA_MEAN_H > ZPWCAP_CEP) THEN
           ZSWI_CEP_H = 1.0_JPRB
         ELSEIF (ZTETA_MEAN_H > ZPWWILT_CEP) THEN
           ZSWI_CEP_H = (ZTETA_MEAN_H - ZPWWILT_CEP)/(ZPWCAP_CEP - ZPWWILT_CEP)
         ELSE
           ZSWI_CEP_H = 0.0_JPRB
         ENDIF

         ZTETA_MEAN_L = MAX(ZPLIQUID_LAYERS(1)/ZPSOL,ZPWWILT_CEP)*ZRVROOTSA(1,NINT(ZFTVL(JGPTOTG)))&
 &            + MAX(ZPLIQUID_LAYERS(2)/ZPSOL,ZPWWILT_CEP)*ZRVROOTSA(2,NINT(ZFTVL(JGPTOTG)))&
 &            + MAX(ZPLIQUID_LAYERS(3)/ZPSOL,ZPWWILT_CEP)*ZRVROOTSA(3,NINT(ZFTVL(JGPTOTG)))&
 &            + MAX(ZPLIQUID_LAYERS(4)/ZPSOL,ZPWWILT_CEP)*ZRVROOTSA(4,NINT(ZFTVL(JGPTOTG)))

         IF (ZTETA_MEAN_L > ZPWCAP_CEP) THEN
           ZSWI_CEP_L = 1.0_JPRB
         ELSEIF (ZTETA_MEAN_L > ZPWWILT_CEP) THEN
           ZSWI_CEP_L = (ZTETA_MEAN_L - ZPWWILT_CEP)/(ZPWCAP_CEP - ZPWWILT_CEP)
         ELSE
           ZSWI_CEP_L = 0.0_JPRB
         ENDIF

! pseudo evaporation vegetation high and low
         ZEVH = ZSWI_CEP_H
         ZEVL = ZSWI_CEP_L

         IF (ZPLIQUID_LAYERS(1)/ZPSOL > ZPWCAP_CEP) THEN
           ZSWI_CEP = 1.0_JPRB
         ELSEIF (ZPLIQUID_LAYERS(1)/ZPSOL > ZPWWILT_CEP) THEN
           ZSWI_CEP= ((ZPLIQUID_LAYERS(1)/ZPSOL)-ZPWWILT_CEP)/(ZPWCAP_CEP-ZPWWILT_CEP)
         ELSE
           ZSWI_CEP = 0.0_JPRB
         ENDIF

         ZVEGH = ZFCVH(JGPTOTG)*ZRVCOV(NINT(ZFTVH(JGPTOTG)))
         ZVEGL = ZFCVL(JGPTOTG)*ZRVCOV(NINT(ZFTVL(JGPTOTG)))

! total pseudo evaporation IFS
         ZEVIFS = (ZVEGH*ZEVH+ZVEGL*ZEVL+(ZFALSM(JGPTOTG)-ZVEGH-ZVEGL)*ZSWI_CEP) /&
               &ZFALSM(JGPTOTG)

         ZSWI_CEP = ZEVIFS

         IF (LLLAISCAL) THEN

!     does rsmin/lai scaling

           IF ((ZFPLAI(JGPTOTG) == 0.0).OR.( ZFPRSMIN(JGPTOTG) == 0.0 ).OR.&
           &(ZRVLAI(NINT(ZFTVH(JGPTOTG))) == 0.0).OR.(ZRVLAI(NINT(ZFTVL(JGPTOTG))) == 0.0).OR.&
           &(ZVEGH + ZVEGL == 0.0)) THEN
             ZLAISCAL = 1.0_JPRB
           ELSE
             ZRSMIN_LAI_H = ZRVRSMIN(NINT(ZFTVH(JGPTOTG))) / ZRVLAI(NINT(ZFTVH(JGPTOTG)))
             ZRSMIN_LAI_L = ZRVRSMIN(NINT(ZFTVL(JGPTOTG))) / ZRVLAI(NINT(ZFTVL(JGPTOTG)))
             ZRSMIN_LAI_IFS = (ZRSMIN_LAI_H * ZVEGH + ZRSMIN_LAI_L * ZVEGL ) / (ZVEGH + ZVEGL)
             ZRSMIN_LAI_ISBA = ZFPRSMIN(JGPTOTG) / ZFPLAI(JGPTOTG)
             ZLAISCAL = ZRSMIN_LAI_IFS / ZRSMIN_LAI_ISBA
           ENDIF

      !CALCUL w ISBA
           ZGOL = (ZEVIFS*(ZFPWFC(JGPTOTG)-ZFPWWILT(JGPTOTG)) + ZLAISCAL*ZFPVEG(JGPTOTG)*ZFPWWILT(JGPTOTG))*ZFPWFC(JGPTOTG) /&
             &( ZFPVEG(JGPTOTG)*ZFPWFC(JGPTOTG)*(ZLAISCAL-1.0) + ZFPWFC(JGPTOTG) + (ZFPVEG(JGPTOTG)-1.0)*ZFPWWILT(JGPTOTG)  )

         ELSE

!     does OI_Trick

           ZGOL = ZEVIFS*(ZFPWFC(JGPTOTG)-ZFPVEG(JGPTOTG)*ZFPWWILT(JGPTOTG)) + ZFPVEG(JGPTOTG)*ZFPWWILT(JGPTOTG)

         ENDIF

         IF (ZGOL > ZFPWFC(JGPTOTG)) ZGOL = ZFPWFC(JGPTOTG)
         IF (ZGOL < ZFPWWILT(JGPTOTG)) ZGOL = ZFPWWILT(JGPTOTG)

      !Calcul w ISBA en kg/m3
          ZGOL=ZGOL*ZSCALE
          ZGOL=ZGOL*ZFPD2(JGPTOTG)
          ZFAGG(JGPTOTG)=ZGOL

     ELSE !if not LLSWI, old computation

      IF (LLSLT) THEN
        ZWSATSLT =ZWP_SAT(NINT(ZFSLT(JGPTOTG)))
      ELSE
        ZWSATSLT=PPWSAT_CEP
      ENDIF
      ZVAL=((ZFSWL1(JGPTOTG)*PPZD(1)/ZPSOL) +&
       & (ZFSWL2(JGPTOTG)*PPZD(2)/ZPSOL) +&
       & (ZFSWL3(JGPTOTG)*PPZD(3)/ZPSOL) +&
       & (ZFSWL4(JGPTOTG)*PPZD(4)/ZPSOL)) /&
       & (PPZD(1)+PPZD(2)+PPZD(3)+PPZD(4))
      ZGOL=(ZVAL*ZSCALE*ZFPD2(JGPTOTG)*ZFPWSAT(JGPTOTG))/&
       & ZWSATSLT  
       IF (ZFPROFT(JGPTOTG) > (TSRESERV2+TDELTA2)) THEN
         ZPGLACE=0.0
       ELSEIF (ZFPROFT(JGPTOTG) < (TSRESERV2-TDELTA2)) THEN
         ZPGLACE=1.0
       ELSE
         ZDELTA=((ZFPROFT(JGPTOTG)-TSRESERV2)/(2.0*TDELTA2)) - 0.5
         ZPGLACE=MIN(((3.0_JPRB*ZDELTA*ZDELTA)&
          & + (2.0_JPRB*ZDELTA*ZDELTA*ZDELTA)),1.0_JPRB)  
       ENDIF
       ZFGLACE(JGPTOTG)=MIN(ZGOL*ZPGLACE,150.0_JPRB)
       ZFAGG(JGPTOTG)=ZGOL-ZFGLACE(JGPTOTG)
     ENDIF !OF IF LLSWI
    ELSE
      ZFAGG(JGPTOTG)=SDEPX*GCONV
      ZFGLACE(JGPTOTG)=0.0_JPRB
    ENDIF
  ENDDO

  CALL FAIENC(IREP,NPOSSH,'PROF',0,'RESERV.EAU  ',ZFAGG,.FALSE.)
  WRITE(NULOUT,*)'---------------------------------------------'
  CLLIBELLE='PROFRESERV.EAU field'
  CALL CPREP1_DIAG(ZFAGG,CLLIBELLE)
  WRITE(NULOUT,*) 'Write grid-point field '
  IF (LLCDST.AND.LLSTL4) THEN
    CALL FAIENC(IREP,NPOSSH,'PROF',0,'RESERV.GLACE',ZFGLACE,.FALSE.)
    WRITE(NULOUT,*)'---------------------------------------------'
    CLLIBELLE='PROFRESERV.GLACE field'
    CALL CPREP1_DIAG(ZFGLACE,CLLIBELLE)
    WRITE(NULOUT,*) 'Write grid-point field '
  ENDIF

ELSE
  WRITE(NULOUT,*)&
   & 'PROFRESERV.EAU NOT COMPUTED BECAUSE MISSING PARAMETERS IN GRIB FILE'  
  IF (LLCONTROL) CALL ABOR1('CPREP1 : ABOR1 CALLED')
ENDIF


IF(LLGRB190.AND.LLGRB191.AND.LLGRB192.AND.LLGRB193) THEN

! ** Set orography variance.

! One first sets zfagg to zero.
! Then one cumulates in zfagg the 190, 191, 192, 193
! grib code fields.
! The result is put on FA file.

  ZFAGG(1:NGPTOTG)=0.0_JPRB
  DO JSTOIA=1,IGRBSN
    IF(IGRBS(JSTOIA) == 190&
       & .OR.IGRBS(JSTOIA) == 191&
       & .OR.IGRBS(JSTOIA) == 192&
       & .OR.IGRBS(JSTOIA) == 193) THEN  

! The stored field has to be cumulated in the total
! variance.

      ZFAGG(1:NGPTOTG)=ZFAGG(1:NGPTOTG)+ZFAGGA(1:NGPTOTG,JSTOIA)
    ENDIF
  ENDDO
  DO JGPTOTG=1,NGPTOTG
    ZFAGG(JGPTOTG)=RG*SQRT(ZFAGG(JGPTOTG))
  ENDDO
  CLLIBELLE='Orography variance'
  CALL CPREP1_DIAG(ZFAGG,CLLIBELLE)

  LLOUTFS=.FALSE.
  CLPREF='SURF'
  INIVEAU=0
  CLSUFF='ET.GEOPOTENT'
  CALL FAIENC(IREP,NPOSSH,CLPREF,INIVEAU,CLSUFF,ZFAGG,LLOUTFS)
  WRITE(NULOUT,*) 'Write grid-point field ',CLPREF(1:4),INIVEAU,CLSUFF(1:16)
ENDIF


! ** Close FA file.

CALL LFILAF(IREP,NPOSSH,.FALSE.)
CALL FAIRME(IREP,NPOSSH,'UNKNOWN')
IF (LLCLIM) THEN
  CALL FAIRME(IREP,IUCLIM,'UNKNOWN')
ENDIF

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CPREP1',1,ZHOOK_HANDLE)

! -------------------------------------------------------------------

CONTAINS

SUBROUTINE CPREP1_DIAG(P_ZZFAGG,CDFIELD)
IMPLICIT NONE

REAL(KIND=JPRB), INTENT(IN), DIMENSION(YDGEOMETRY%YRGEM%NGPTOTG) :: P_ZZFAGG
CHARACTER(LEN=*), INTENT(IN) :: CDFIELD

REAL(KIND=JPRB) :: ZZMAX, ZZMIN
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('CPREP1:CPREP1_DIAG',0,ZHOOK_HANDLE)
ASSOCIATE(NGPTOTG=>YDGEOMETRY%YRGEM%NGPTOTG)
ZZMIN=MINVAL(P_ZZFAGG(1:NGPTOTG))
ZZMAX=MAXVAL(P_ZZFAGG(1:NGPTOTG))
!WRITE(NULOUT,*)'---------------------------------------------'
WRITE(NULOUT,*) TRIM(CDFIELD),' Ave=',SUM(P_ZZFAGG(1:NGPTOTG))/REAL(NGPTOTG,KIND=JPRB),'Min=',ZZMIN,' Max=',ZZMAX
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CPREP1:CPREP1_DIAG',1,ZHOOK_HANDLE)

END SUBROUTINE CPREP1_DIAG

END SUBROUTINE CPREP1

