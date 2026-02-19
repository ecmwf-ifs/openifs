! (C) Copyright 1989- Meteo-France.

SUBROUTINE HPOS(YDRQPHY,YDRQCLI,YDNAMFPSCI,YDAFN,LDFPOSHOR,YDGEOMETRY,YDSURF,YDMODEL,KEND,KSTGLO,&
 & PSP_SB, PSP_SG, PSP_SL, PSP_RR, PSP_CI, PSD_WS, PSD_VD, PSD_VX,&
 & PSD_VF, PSD_VV, PSD_VP, PSD_VA, PSD_VC, PSD_X2,PSD_SM, PSD_OC, PWSXI, PWDXI,&
 & PFPBUF1)

!**** *HPOS*  - HORIZONTAL POST-PROCESSING

!     PURPOSE.
!     --------
!        FILL THE "SEMI-LAGRANGIAN" BUFFERS WITH THE FIELDS TO TRANSFER
!        OR INTERPOLATE

!**   INTERFACE.
!     ----------
!       *CALL* *HPOS*

!        EXPLICIT ARGUMENTS
!        --------------------
!         INPUT:
!          LDFPOSHOR  : .true. if actual horizontal interpolations
!          KEND       : last  adress in post-processing buffers to read
!          KSTGLO     : global offset (used in message passing version)
!          PSP_SB     : soil prognostic quantities for the different reservoirs
!          PSP_SG     : surface snow prognostic quantities
!          PSP_SL     : lake prognostic quantities
!          PSP_RR     : surface prognostic quantities 
!          PSP_CI     : 2-d prognostic fields for CANARI
!          PSD_WS     : surface prognostic quantities over sea 
!          PSD_VD     : (ECMWF) diagnostic fields
!          PSD_VX     : auxilary climatological diagnostic fields
!          PSD_VV     : vegetation diagnostic fields
!          PSD_VP     : deep soil diagnostic fields
!          PSD_VA     : aerosol diagnostic fields
!          PSD_VC     : climatological ozone profiles diagnostic fields
!          PSD_X2     : extra 2-d diagnostic fields
!          PSD_SM     : simulated satellite images
!          PSD_VF     : climatological/geographical diagnostic fields
!          PSD_OC     ! ocean fields

!         INPUT/OUTPUT:
!          PFPBUF1    : interpolation buffer containing the fields to interpolate

!        IMPLICIT ARGUMENTS
!        --------------------
!        se #include below.

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------
!     CVLANISO,ACSOLW

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS
!        See documentation about FULL-POS.

!     AUTHOR.
!     -------
!        RYAD EL KHATIB *METEO-FRANCE*
!        ORIGINAL : 94-04-08

!     MODIFICATIONS.
!     --------------
!    R. El Khatib : 01-04-03 Derived type  - cleanings
!    R. El Khatib : 01-05-04 Duration of total precipitations,HUn,HUx
!    R. El Khatib : 01-08-07 Pruning options
!    R. El Khatib : 02-21-20 Fullpos B-level distribution + remove IO scheme
!    D.Dent       : 02-05-15 Improve handling of 2D extra fields
!    R. El Khatib : 03-04-17 Fullpos improvemnts
!    JJMorcrette  : 02-11-04 PAR, UV-B, CAPE
!    M.Hamrud     : 01-Oct-2003 CY28 Cleaning
!    F.Taillefer  : 23-Jun-2004 Add O3 and aerosols fields
!    Y.Bouteloup  : 2004-03-26 introduction of MTS radiance
!    M.Hamrud      10-Jan-2004 CY28R1 Cleaning
!    A. Tompkins  : 11-02-04 TCLW,TCIW
!    M.Hamrud      01-Jul-2006 Revised surface fields
!    S.Serrar     : 26-09-2005  total column for CO2 and SF6
!    S. Serrar    : 07-Sep-06 include tracers for diagnostics (GEMS)
!    R. Engelen   : 02-Jun_06 CO2 and SF6 replaced by generic GHG
!    JJMorcrette 20060721 PP of clear-sky PAR and TOA solar incident radiation
!    JJMorcrette 20060807 PP of VIMD
!    JJMorcrette  : 20060625 MODIS albedo
!    R. El Khatib : 02-21-20 Bugfix on the ordering of dynamic fields
!    JJMorcrette  : 20060925 DU, BC, OM, SO2, SOA climatological fields
!    G. Balsamo   : 20070115 Soil type
!    S. Serrar    : 20070717 methane surface fields
!    Y. Seity     : 07-01-15 Add graupel and hail precipitations 
!    H. Hersbach  : 20080606 Add ocean current
!    P. Lopez     : 26-Jun-2009 Added 100 m zonal and meridional wind fields
!    P. Lopez     : 22-Oct-2008 Added ducting diagnostic fields
!    K. Yessad (Jul 2009): remove CDLOCK + some cleanings
!    JJMorcrette  : 20091201 Total and clear-sky direct SW radiation flux at surface 
!    H. Hersbach  : 04-Dec-2009 Introduce 10m-neutral wind and friction velocity
!    S. Boussetta/G.Balsamo  05-June-2009   Added variable LAI fields
!    R. Forbes    : 01-Mar-2010: Add TCRW, TCSW diagnostics
!    JJMorcrette  : 20100212  PP of CBASE, 0DEGL and VISIH
!    T. Wilhelmsson : 25-Mar-10 Add 6 hourly min/max fields
!    T. Wilhelmsson : 20-Dec-10 Add 3 hourly min/max fields
!    G.Balsamo/S.Boussetta : 17-Apr-11 Add land carbon dioxide
!    R. El Khatib 20-Aug-2012 GAUXBUF removed and replaced by HFPBUF
!    M. Ahlgrimm 31 Oct 2011   Clear-sky downward radiation at surface
!    R.Engelen    (Oct 2011): Clean-up SF6
!    A. Inness    : 23-Mar-2012 Add total column CHEM fields
!    R. El Khatib  04-Dec-2012 Fix bounds checking issues created by recent "cleanings"
!    JJMorcrette 20130213 PP optical depths GEMS/MACC aerosols
!    G. Balsamo   : 14-Jan-2014 Add lake model fields
!    R. Forbes     01-Mar-2014 Added precip rates/type,TCSLW,I10FG,PEV
!    T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!    JJMorcrette 20130730 15-aerosol variable model + oceanic DMS
!    A. Agusti-Panareda (3 Nov 2013): Add GPP and REC flux adjustment coefficients
!    P.Marguinaud 20141001 Call HPOSSFX for PREP fields
!    P.Marguinaud 20141010 More fields
!    R. Forbes     10-Jan-2015 Added freezing rain FZRA
!    P. Lopez 16-Nov-2015 Added lightning density fields
!    P. Bechtold: 10-Nov-2015 Add CBASEA, CTOPC, ZTWETB
!    R. El Khatib 17-Aug-2016 Interoperability IFS vs ARPEGE
!    R. El Khatib 09-Sep-2016 More interoperability + cleanings
!    R. El Khatib 20-Sep-2016 Split HPOS too big
!    R. Forbes    10-Apr-2017 Added total precip rate
!    R. El Khatib 13-Sep-2017 Interoperability Surfex -> ISBA (following Bouyssel, Fischer & Degrauwe, 2017)
!    E.Dutra/G.Arduini Jan 2018: change of PSP_SG to 4 Dimensions, snow multi-layer 
!    R. El Khatib 12-Mar-2018 total albedo
!    R. El Khatib 14-May-2018 fix uninitialized variable
!    B. Ingleby   14-Jan-2019 add Y2SH
!    R. Hogan     15-Jan-2019 Added 6-component MODIS albedo
!    P. Bechtold  01-Jan-2020 Added M10DEGL more CAPE/CIN parameters
!    R. Forbes    May-2022    Added calculations of most severe and most frequent precip-types
!     ------------------------------------------------------------------

USE TYPE_MODEL         , ONLY : MODEL
USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE PARKIND1           , ONLY : JPIM, JPRB
USE YOMHOOK            , ONLY : LHOOK, DR_HOOK, JPHOOK
USE PARFPOS            , ONLY : JPOSFSU, JPOSVX2
USE YOMAFN             , ONLY : TAFN
USE YOMFPC             , ONLY : TNAMFPSCI
USE YOMCST             , ONLY : RPI, RG
USE YOMSTA             , ONLY : RDTDZ1
USE YOMFP4L            , ONLY : TRQFP, IFPSEARCH
USE YOM_YGFL           , ONLY : JPGHG
USE YOM_GRIB_CODES     , ONLY : NGRBGHG
USE YOMCLI             , ONLY : SMASK, NTPDES, NTPMER
USE YOE_AERODIAG       , ONLY : NPAERODIAG
USE YOECLDP            , ONLY : NPRECTYPES

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE (TRQFP)      ,INTENT(IN)    :: YDRQPHY
TYPE (TRQFP)      ,INTENT(IN)    :: YDRQCLI
TYPE (TNAMFPSCI)  ,INTENT(IN)    :: YDNAMFPSCI
TYPE (TAFN)       ,INTENT(IN)    :: YDAFN
LOGICAL           ,INTENT(IN)    :: LDFPOSHOR
TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TSURF)       ,INTENT(IN)    :: YDSURF
TYPE(MODEL)       ,INTENT(IN)    :: YDMODEL
INTEGER(KIND=JPIM),INTENT(IN)    :: KEND 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTGLO 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSP_SB(:,:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSP_SG(:,:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSP_SL(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSP_RR(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSP_CI(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSD_WS(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSD_VD(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSD_VX(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSD_VF(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSD_VV(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSD_VP(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSD_VA(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSD_VC(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSD_X2(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSD_SM(:,:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSD_OC(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWSXI
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWDXI
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFPBUF1(YDGEOMETRY%YRDIM%NPROMA,YDRQPHY%NFIELDG) 

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZTD(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB) :: ZWFC(YDGEOMETRY%YRDIM%NPROMA),ZWPMX(YDGEOMETRY%YRDIM%NPROMA),ZWSAT(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB) :: ZWSMX(YDGEOMETRY%YRDIM%NPROMA), ZIVG(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB) :: ZWWILT(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB) :: ZSSN(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB) :: ZTSN(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB) :: ZRSN(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB) :: ZDUMM(YDGEOMETRY%YRDIM%NPROMA),ZGAMMA(YDGEOMETRY%YRDIM%NPROMA), ZALFA(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB) :: ZGM(YDGEOMETRY%YRDIM%NPROMA), ZOROG(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB) :: ZECWSAT(YDGEOMETRY%YRDIM%NPROMA), ZECWCAP(YDGEOMETRY%YRDIM%NPROMA), ZECWWILT(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB) :: ZFSWL(YDMODEL%YRML_PHY_G%YRDPHY%NCSS), ZICEL(YDMODEL%YRML_PHY_G%YRDPHY%NCSS), &
 &  ZRDAT(YDMODEL%YRML_PHY_G%YRDPHY%NCSS),ZTETA_H(YDMODEL%YRML_PHY_G%YRDPHY%NCSS),ZTETA_L(YDMODEL%YRML_PHY_G%YRDPHY%NCSS)
REAL(KIND=JPRB), ALLOCATABLE :: ZRVCOV(:), ZRVROOTSA(:,:), ZRWSATM(:), ZRWPWPM(:), ZRWCAPM(:), ZRVLAI(:), ZRVRSMIN(:)

REAL(KIND=JPRB) :: ZRSMIN_LAI_H, ZRSMIN_LAI_L, ZLAISCAL, ZRSMIN_LAI_IFS, ZRSMIN_LAI_ISBA, ZVEG_ISBA

INTEGER(KIND=JPIM) :: IEXTR2, I1DIM, I3DIM, JVAR, JFLD, IC, JI, IANO, IFLD, IST, JK, JP
INTEGER(KIND=JPIM) ::  IBL, ICLBT, ICSBT, JAERO, JEMIS, ISDML,IRSNML,ITSNML,IWSNML
INTEGER(KIND=JPIM) :: JGHG, IGHG, JCHEM, JCHEM1, INDEX, IVTYPES, ISOILTYPE, JCSS, ITVH, ITVL, ISOTY
REAL(KIND=JPRB) :: ZU(YDGEOMETRY%YRDIM%NPROMA), ZV(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB) :: ZEPS, ZWPI(YDGEOMETRY%YRDIM%NPROMA), ZSNS(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB) :: ZWPL(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB) :: ZSCALWPI(YDGEOMETRY%YRDIM%NPROMA), ZSCALSNS(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB) :: ZPTYPESUM(YDGEOMETRY%YRDIM%NPROMA, NPRECTYPES)

LOGICAL :: LLHMT, LLU, LLV
LOGICAL :: LLSCALWPI, LLSCALSNS
LOGICAL :: LLPRECIPMASK(YDGEOMETRY%YRDIM%NPROMA, NPRECTYPES)

REAL(KIND=JPRB) :: Z1SWSX, ZDTSGDZ, ZPSOL
REAL(KIND=JPRB) :: ZVEGH, ZVEGL, ZTETA_MEAN_H, ZTETA_MEAN_L, ZEVH, ZEVL, ZEVIFS, ZECSWI, ZDELTA

INTEGER(KIND=JPIM)  :: JDIAG
REAL(KIND=JPRB) :: ZTMP,ZTMP1,ZRTT,ZRLMLT,ZIHCAP
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "acsolw.intfb.h"
#include "cvlaniso.intfb.h"

#include "surf_inq.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('HPOS',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM, YDDIMV=>YDGEOMETRY%YRDIMV, YDPHY=>YDMODEL%YRML_PHY_MF%YRPHY, &
 & YDDPHY=>YDMODEL%YRML_PHY_G%YRDPHY,YDPHY1=>YDMODEL%YRML_PHY_MF%YRPHY1, &
 & YGFL=>YDMODEL%YRML_GCONF%YGFL,YDEPHY=>YDMODEL%YRML_PHY_EC%YREPHY, &
 & YDEAERATM=>YDMODEL%YRML_PHY_RAD%YREAERATM,YDECLDP=>YDMODEL%YRML_PHY_EC%YRECLDP)
ASSOCIATE(NFIELDG=>YDRQPHY%NFIELDG, &
 & TSRESERV2=>YDNAMFPSCI%TSRESERV2, TDELTA2=>YDNAMFPSCI%TDELTA2, NFPSWI=>YDNAMFPSCI%NFPSWI, &
 & RWPITPN=>YDNAMFPSCI%RWPITPN, RWPITPX=>YDNAMFPSCI%RWPITPX, RSNSTPN=>YDNAMFPSCI%RSNSTPN, &
 & RSNSTPX=>YDNAMFPSCI%RSNSTPX, RSNSMOD=>YDNAMFPSCI%RSNSMOD, NFPCLI=>YDNAMFPSCI%NFPCLI, &
 & GFP=>YDAFN%GFP, GFP_PHYDS=>YDAFN%GFP_PHYDS, &
 & NCHEM=>YGFL%NCHEM, NCHEM_FLXO=>YGFL%NCHEM_FLXO, NGHG=>YGFL%NGHG, NACTAERO=>YGFL%NACTAERO, &
 & NCHEM_WDFLX=>YGFL%NCHEM_WDFLX,  NCHEM_DDFLX=>YGFL%NCHEM_DDFLX, &
 & NEMIS2D=>YGFL%NEMIS2D, NEMIS2DAUX=>YGFL%NEMIS2DAUX, &
 & YCHEM=>YGFL%YCHEM, YGHG=>YGFL%YGHG, &
 & YAERO_WVL_DIAG_NL=>YGFL%YAERO_WVL_DIAG_NL, &
 & NAERO_WVL_DIAG=>YGFL%NAERO_WVL_DIAG, &
 & NAERO_WVL_DIAG_TYPES=>YGFL%NAERO_WVL_DIAG_TYPES, &
 & YAERO_DESC=>YDEAERATM%YAERO_DESC, &
 & NPROMA=>YDDIM%NPROMA, NFLEVG=>YDDIMV%NFLEVG, &
 & LEPHYS=>YDEPHY%LEPHYS, YSURF=>YDEPHY%YSURF, NCSS=>YDDPHY%NCSS, &
 & YSD_SM=>YDSURF%YSD_SM, YSD_VA=>YDSURF%YSD_VA, YSD_VC=>YDSURF%YSD_VC, &
 & YSD_VD=>YDSURF%YSD_VD, YSD_VF=>YDSURF%YSD_VF, YSD_VP=>YDSURF%YSD_VP, &
 & YSD_VV=>YDSURF%YSD_VV, YSD_VX=>YDSURF%YSD_VX, YSD_WS=>YDSURF%YSD_WS, &
 & YSD_X2=>YDSURF%YSD_X2, YSP_CI=>YDSURF%YSP_CI, YSP_RR=>YDSURF%YSP_RR, &
 & YSP_SB=>YDSURF%YSP_SB, YSP_SG=>YDSURF%YSP_SG, YSP_SL=>YDSURF%YSP_SL, &
 & YSD_OC=>YDSURF%YSD_OC, &
 & LMPHYS=>YDPHY%LMPHYS, NCSNEC=>YDDPHY%NCSNEC, &
 & GCONV=>YDPHY1%GCONV, RD1=>YDPHY1%RD1, NTVGLA=>YDPHY1%NTVGLA, &
 & NPTYPE_SEV2WMO=>YDECLDP%NPTYPE_SEV2WMO, &
 & LNEMOCOUP=>YDMODEL%YRML_AOC%YRMCC%LNEMOCOUP, &
 & LNEMOGRIBMASK=>YDMODEL%YRML_AOC%YRMCC%LNEMOGRIBMASK, &
 & LNEMOGRIBFLDS=>YDMODEL%YRML_AOC%YRMCC%LNEMOGRIBFLDS )

!     ------------------------------------------------------------------

!*       1. PRELIMINARY INITIALISATIONS.
!           ----------------------------

!     IST     : first adress in post-processing buffers to read

IBL=(KSTGLO-1)/NPROMA+1
IST=1


!*       5. COPY PHYSICAL FIELDS FROM MODEL ARRAYS/WORKFILES TO BUFFER PFPBUF1.
!           ------------------------------------------------------------------

!*       5.1 various preinitializations

ZDTSGDZ=RDTDZ1/RG
ZGM(IST:KEND)=YDGEOMETRY%YRGSGEOM(IBL)%GM(IST:KEND)
ZOROG (1:KEND)=YDGEOMETRY%YROROG(IBL)%OROG (1:KEND)
ZDUMM(1:KEND)=1.0_JPRB
ICLBT=0
ICSBT=0
ISDML=0
IRSNML=0
ITSNML=0
IWSNML=0


IF (LEPHYS) THEN
  ! prepare to convert from TESSEL to ISBA
  DO JFLD=1,NFIELDG
    IFLD=JFLD
    IC=YDRQPHY%ICOD(IFLD)
    IF (IC == GFP%SSW%ICOD .OR. IC == GFP%DSW%ICOD .OR. IC == GFP%VEG%ICOD&
     & .OR. IC == GFP%FSSW%ICOD .OR. IC == GFP%FDSW%ICOD) THEN
      CALL SURF_INQ(YSURF,KNVTYPES=IVTYPES)
      ALLOCATE(ZRVCOV(0:IVTYPES))
      ALLOCATE(ZRVLAI(0:IVTYPES))
      ALLOCATE(ZRVRSMIN(0:IVTYPES))
      ALLOCATE(ZRVROOTSA(NCSS,0:IVTYPES))
      CALL SURF_INQ(YSURF,PRVCOV=ZRVCOV,PRVROOTSA=ZRVROOTSA,PRDAT=ZRDAT,PRVLAI=ZRVLAI,PRVRSMIN=ZRVRSMIN)
      IF (YSD_VF%YSOTY%LSET) THEN
        CALL SURF_INQ(YSURF,KNSOTY=ISOTY)
        ALLOCATE(ZRWPWPM(0:ISOTY))
        ALLOCATE(ZRWCAPM(0:ISOTY))
        ALLOCATE(ZRWSATM(0:ISOTY))
        CALL SURF_INQ(YSURF,PRWPWPM=ZRWPWPM,PRWCAPM=ZRWCAPM,PRWSATM=ZRWSATM)
        IF (NFPSWI > 0) THEN
          DO JI=IST,KEND
            ISOILTYPE=NINT(PSD_VF(JI,YSD_VF%YSOTY%MP))
            ZECWWILT(JI)=ZRWPWPM(ISOILTYPE)
            ZECWCAP(JI)=ZRWCAPM(ISOILTYPE)
          ENDDO
        ELSE
          DO JI=IST,KEND
            ISOILTYPE=NINT(PSD_VF(JI,YSD_VF%YSOTY%MP))
            ZECWWILT(JI)=ZRWPWPM(ISOILTYPE)
            ZECWSAT(JI) =ZRWSATM(ISOILTYPE)
          ENDDO
        ENDIF
        DEALLOCATE(ZRWCAPM)
        DEALLOCATE(ZRWSATM)
        DEALLOCATE(ZRWPWPM)
      ELSE
        ZECWWILT(IST:KEND)=0.171_JPRB ! conventional wilting point for TESSEL ?
        ZECWSAT(IST:KEND)=0.472_JPRB ! conventional field capacity for TESSEL ?
      ENDIF
      EXIT
    ENDIF
  ENDDO
  ! Notice : the value anterior to June 2000 = 1._JPRB/ZRDAT(1) is not supported here
  ! The variable below is used to keep the history of the code
  ZPSOL=1._JPRB
ENDIF

ZEPS=EPSILON(1._JPRB)*1000._JPRB
IF (LMPHYS) THEN
  ! prepare to convert from SURFEX to ISBA
  LLSCALWPI=(RWPITPX-RWPITPN > ZEPS)
  IF (LLSCALWPI) THEN
    IF (YSP_SB%YT%LSET.AND.NCSS>=1) THEN
      ZSCALWPI(IST:KEND)=(0.5_JPRB+0.5_JPRB*COS(RPI*(MAX(RWPITPN,MIN(PSP_SB(IST:KEND,1,YSP_SB%YT%MP0),RWPITPX))-RWPITPN) &
       & /(RWPITPX-RWPITPN)))
    ELSE
      CALL ABOR1('HPOS : Tp IS NEEDED FOR CORRECTION SURFEX -> ISBA !')
    ENDIF
  ENDIF
  LLSCALSNS=(RSNSTPX-RSNSTPN > ZEPS .AND. ABS(RSNSMOD-1._JPRB) > ZEPS)
  IF (LLSCALSNS) THEN
    IF (YSP_SB%YT%LSET.AND.NCSS>=1) THEN
      ZSCALSNS(IST:KEND)=(0.5_JPRB+0.5_JPRB*COS(RPI*(MAX(RSNSTPN,MIN(PSP_SB(IST:KEND,1,YSP_SB%YT%MP0),RSNSTPX))-RSNSTPN) &
       & /(RSNSTPX-RSNSTPN)))
    ELSE
      CALL ABOR1('HPOS : Tp IS NEEDED FOR CORRECTION SURFEX -> ISBA !')
    ENDIF
  ENDIF
ELSE
  LLSCALWPI=.FALSE.
  LLSCALSNS=.FALSE. 
ENDIF


!*       5.2 Fill buffer.

! preset to zero because fields can be unknown to the post-processing
PFPBUF1(:,:)=0.0_JPRB

IEXTR2=0

  !* snow single-layer
  !* pre-compute aggregater snow variables before loop
IF ( NCSNEC > 1 .AND. YSP_SG%YF%LSET .AND. YSP_SG%YR%LSET .AND. YSP_SG%YT%LSET) THEN
  ZRTT=273.16_JPRB
  ZRLMLT=2.8345E+6_JPRB-2.5008E+6_JPRB
  ZIHCAP=2.05E6_JPRB/920._JPRB
  DO JI=IST,KEND
    ZSSN(JI)=0._JPRB ! snow mass 
    ZTMP=0._JPRB     ! snow depth 
    ZTMP1=0._JPRB    ! snow internal energy
    DO JK=1,NCSNEC
      IF (PSP_SG(JI,JK,YSP_SG%YF%MP0) > 0._JPRB ) THEN
        ZSSN(JI) = ZSSN(JI) + PSP_SG(JI,JK,YSP_SG%YF%MP0)
        ZTMP = ZTMP + PSP_SG(JI,JK,YSP_SG%YF%MP0) / PSP_SG(JI,JK,YSP_SG%YR%MP0)
        ZTMP1 = ZTMP1 + ZIHCAP*PSP_SG(JI,JK,YSP_SG%YF%MP0)*&
         &              (PSP_SG(JI,JK,YSP_SG%YT%MP0)-ZRTT) -&
         &               ZRLMLT*(PSP_SG(JI,JK,YSP_SG%YF%MP0)-PSP_SG(JI,JK,YSP_SG%YW%MP0))
      ENDIF
    ENDDO
    IF (ZSSN(JI) > 0._JPRB .AND. ZTMP > 0._JPRB ) THEN
      ZRSN(JI) = ZSSN(JI) / ZTMP ! mean snow density
    ! ZTSN(JI)=MIN(ZRTT , ZRTT + (ZTMP1+ZRLMLT*ZSSN(JI))/&
    !                            (ZIHCAP*ZSSN(JI)))
      ZTSN(JI)=PSP_SG(JI,1,YSP_SG%YT%MP0)
    ELSE
      ZRSN(JI)=PSP_SG(JI,1,YSP_SG%YR%MP0)
      ZTSN(JI)=PSP_SG(JI,1,YSP_SG%YT%MP0)
    ENDIF
  ENDDO
ELSEIF (NCSNEC == 1) THEN
  IF (YSP_SG%YF%LSET) ZSSN(IST:KEND) = PSP_SG(IST:KEND,1,YSP_SG%YF%MP0)
  IF (YSP_SG%YR%LSET) ZRSN(IST:KEND) = PSP_SG(IST:KEND,1,YSP_SG%YR%MP0)
  IF (YSP_SG%YT%LSET) ZTSN(IST:KEND) = PSP_SG(IST:KEND,1,YSP_SG%YT%MP0)
ENDIF
               

DO JFLD=1,NFIELDG

  IFLD=JFLD
  IC=YDRQPHY%ICOD(IFLD)

  IF (LDFPOSHOR.OR.NFPCLI > 0) THEN
    IANO=GFP_PHYDS(IC)%IANO
  ELSE
    ! No interpolation nor change of physiography <=> only buffers transfers
    ! (unless TESSEL to ISBA conversion)
    IANO=0
  ENDIF

    IF   (IC == GFP%LSM%ICOD .OR.(IC == GFP%LAN%ICOD .AND. NFPCLI == 0)) THEN
      IF (YSD_VF%YLSM%LSET) THEN
      ! Land-sea mask
        PFPBUF1(IST:KEND,IFLD)=PSD_VF(IST:KEND,YSD_VF%YLSM%MP)
      ELSEIF (IFPSEARCH(YDRQCLI,GFP%LSM%ICOD) <= 0) THEN
        CALL ABOR1('HPOS : MISSING LAND-SEA MASK IN INPUT/OUTPUT DATA')
      ENDIF
    ELSEIF (IC == GFP%RDST%ICOD) THEN
      IF (LMPHYS.AND.YSP_RR%YT%LSET) THEN
        ! Interpolated surface temperature
        PFPBUF1(IST:KEND,IFLD)=PSP_RR(IST:KEND,YSP_RR%YT%MP0)
      ELSEIF (LEPHYS.AND.YSP_SB%YT%LSET.AND.NCSS>=1) THEN
        ! Interpolated TsLev1
        PFPBUF1(IST:KEND,IFLD)=PSP_SB(IST:KEND,1,YSP_SB%YT%MP0)
      ENDIF
    ELSEIF (IC == GFP%ST%ICOD .OR. (IC == GFP%SST%ICOD.AND..NOT.YSD_VF%YSST%LSET)) THEN
      ! Ts or pseudo-SST
      IF (LMPHYS.AND.YSP_RR%YT%LSET) THEN
        ! Surface temperature
        IF (IANO == 0.OR.NFPCLI <=2) THEN
          ! Interpolation of Ts
          PFPBUF1(IST:KEND,IFLD)=PSP_RR(IST:KEND,YSP_RR%YT%MP0)
        ELSEIF (YSD_VX%YTSC%LSET .AND. YSD_VX%YORO%LSET) THEN
          ! Interpolation of the anomaly Ts-Ts*
          DO JI=IST,KEND
            PFPBUF1(JI,IFLD)=PSP_RR(JI,YSP_RR%YT%MP0)-PSD_VX(JI,YSD_VX%YTSC%MP)&
             & -ZDTSGDZ*(ZOROG(JI)-PSD_VX(JI,YSD_VX%YORO%MP))  
          ENDDO
        ELSE
          CALL ABOR1('HPOS : MISSING FIELDS TO COMPUTE TS ANOMALY')
        ENDIF
      ELSEIF (LEPHYS.AND.YSP_SB%YT%LSET.AND.NCSS>=1) THEN
        ! TsLev1
        IF (IANO == 0.OR.NFPCLI <=2) THEN
          ! Interpolation of TsLev1
          PFPBUF1(IST:KEND,IFLD)=PSP_SB(IST:KEND,1,YSP_SB%YT%MP0)
        ELSEIF (YSD_VX%YTSC%LSET .AND. YSD_VX%YORO%LSET) THEN
          ! Interpolation of the anomaly Ts-Ts*
          DO JI=IST,KEND
            PFPBUF1(JI,IFLD)=PSP_SB(JI,1,YSP_SB%YT%MP0)-PSD_VX(JI,YSD_VX%YTSC%MP)&
             & -ZDTSGDZ*(ZOROG(JI)-PSD_VX(JI,YSD_VX%YORO%MP))  
          ENDDO
        ELSE
          CALL ABOR1('HPOS : MISSING FIELDS TO COMPUTE TSLEV1 ANOMALY')
        ENDIF
      ENDIF
    ELSEIF (IC == GFP%SST%ICOD .AND. YSD_VF%YSST%LSET) THEN
      ! True SST (ECMWF)
      IF (LNEMOCOUP.AND.LNEMOGRIBMASK) THEN
        PFPBUF1(IST:KEND,IFLD)=PSD_OC(IST:KEND,YSD_OC%YSSTM%MP)
      ELSE
        PFPBUF1(IST:KEND,IFLD)=PSD_VF(IST:KEND,YSD_VF%YSST%MP)
      ENDIF

      !        ECMWF extra fields (CNSKT to CNDPAT).

    ELSEIF (IC == GFP%SKT%ICOD .AND.YSP_RR%YT%LSET) THEN
      !  Interpolation of skin temp.
      PFPBUF1(IST:KEND,IFLD)=PSP_RR(IST:KEND,YSP_RR%YT%MP0)
    ELSEIF (IC == GFP%STL1%ICOD) THEN
      !  Interpolation of TsLev1
      PFPBUF1(IST:KEND,IFLD)=PSP_SB(IST:KEND,1,YSP_SB%YT%MP0)
    ELSEIF (IC == GFP%STL2%ICOD) THEN
      !  Interpolation of TsLev2
      PFPBUF1(IST:KEND,IFLD)=PSP_SB(IST:KEND,2,YSP_SB%YT%MP0)
    ELSEIF (IC == GFP%STL3%ICOD ) THEN
      !  Interpolation of TsLev3
      PFPBUF1(IST:KEND,IFLD)=PSP_SB(IST:KEND,3,YSP_SB%YT%MP0)
    ELSEIF (IC == GFP%STL4%ICOD ) THEN
      !  Interpolation of TsLev4
      PFPBUF1(IST:KEND,IFLD)=PSP_SB(IST:KEND,4,YSP_SB%YT%MP0)
    ELSEIF (IC == GFP%SRC%ICOD .AND.YSP_RR%YW%LSET ) THEN
      !  Interpolation of SRC
      PFPBUF1(IST:KEND,IFLD)=PSP_RR(IST:KEND,YSP_RR%YW%MP0)
    ELSEIF (IC == GFP%SWL1%ICOD ) THEN
      !  Interpolation of SWLev1
      PFPBUF1(IST:KEND,IFLD)=PSP_SB(IST:KEND,1,YSP_SB%YQ%MP0)
    ELSEIF (IC == GFP%SWL2%ICOD ) THEN
      !  Interpolation of SWLev2
      PFPBUF1(IST:KEND,IFLD)=PSP_SB(IST:KEND,2,YSP_SB%YQ%MP0)
    ELSEIF (IC == GFP%SWL3%ICOD ) THEN
      !  Interpolation of SWLev3
      PFPBUF1(IST:KEND,IFLD)=PSP_SB(IST:KEND,3,YSP_SB%YQ%MP0)
    ELSEIF (IC == GFP%SWL4%ICOD ) THEN
      !  Interpolation of SWLev4
      PFPBUF1(IST:KEND,IFLD)=PSP_SB(IST:KEND,4,YSP_SB%YQ%MP0)
    ELSEIF (IC == GFP%ISTL1%ICOD ) THEN
      !  Interpolation of ISTLev1
      PFPBUF1(IST:KEND,IFLD)=PSP_SB(IST:KEND,1,YSP_SB%YTL%MP0)
    ELSEIF (IC == GFP%ISTL2%ICOD ) THEN
      !  Interpolation of ISTLev2
      PFPBUF1(IST:KEND,IFLD)=PSP_SB(IST:KEND,2,YSP_SB%YTL%MP0)
    ELSEIF (IC == GFP%ISTL3%ICOD ) THEN
      !  Interpolation of ISTLev3
      PFPBUF1(IST:KEND,IFLD)=PSP_SB(IST:KEND,3,YSP_SB%YTL%MP0)
    ELSEIF (IC == GFP%ISTL4%ICOD ) THEN
      !  Interpolation of ISTLev4
      PFPBUF1(IST:KEND,IFLD)=PSP_SB(IST:KEND,4,YSP_SB%YTL%MP0)
    ELSEIF (IC == GFP%CI%ICOD ) THEN
      !  Interpolation of ice cover.                                
      IF (LNEMOCOUP.AND.LNEMOGRIBMASK) THEN
        PFPBUF1(IST:KEND,IFLD)=PSD_OC(IST:KEND,YSD_OC%YCIM%MP)
      ELSE
        PFPBUF1(IST:KEND,IFLD)=PSD_VF(IST:KEND,YSD_VF%YCI%MP)
      ENDIF
    ELSEIF (IC == GFP%CVL%ICOD ) THEN
      !  Interpolation of CVL.                                
      PFPBUF1(IST:KEND,IFLD)=PSD_VF(IST:KEND,YSD_VF%YCVL%MP)
    ELSEIF (IC == GFP%CVH%ICOD ) THEN
      !  Interpolation of CVH.                              
      PFPBUF1(IST:KEND,IFLD)=PSD_VF(IST:KEND,YSD_VF%YCVH%MP)  
    ELSEIF (IC == GFP%CO2TYP%ICOD ) THEN
      !  Interpolation of CO2TYP.                                
      PFPBUF1(IST:KEND,IFLD)=PSD_VF(IST:KEND,YSD_VF%YCO2TYP%MP)
    ELSEIF (IC == GFP%TVL%ICOD ) THEN
      !  Interpolation of TVL.                                
      PFPBUF1(IST:KEND,IFLD)=PSD_VF(IST:KEND,YSD_VF%YTVL%MP)
    ELSEIF (IC == GFP%TVH%ICOD ) THEN
      !  Interpolation of TVH.                                
      PFPBUF1(IST:KEND,IFLD)=PSD_VF(IST:KEND,YSD_VF%YTVH%MP)
    ELSEIF (IC == GFP%TSN%ICOD ) THEN
      !  Interpolation of snow temperature
      IF (ITSNML == 0) THEN
        PFPBUF1(IST:KEND,IFLD)=ZTSN(IST:KEND)
      ELSE
        IF (YSP_SG%YT%LSET .AND. ITSNML <= SIZE(PSP_SG,DIM=2)) THEN
          PFPBUF1(IST:KEND,IFLD)=PSP_SG(IST:KEND,ITSNML,YSP_SG%YT%MP0)
        ELSE
          CALL ABOR1('HPOS : ERROR IN SNOW TEMPERATURE')
        ENDIF
      ENDIF
      ITSNML = ITSNML + 1
    ELSEIF (IC == GFP%ASN%ICOD ) THEN
      !  Interpolation of snow albedo.                                
      PFPBUF1(IST:KEND,IFLD)=PSP_SG(IST:KEND,1,YSP_SG%YA%MP0)
    ELSEIF (IC == GFP%RSN%ICOD) THEN
      !  Interpolation of snow density.                                
      IF (LMPHYS) THEN
        ! Snow density (adimensional : in kg/kg)
        ! added this *1000. that wasn't here in the *_CY43R3_snowML       
        PFPBUF1(IST:KEND,IFLD)=ZRSN(IST:KEND)*1000._JPRB
      ELSEIF (LEPHYS) THEN
      ! Snow mass per unit of volume in kg.m-3 
      !  PFPBUF1(IST:KEND,IFLD)=PSP_SG(IST:KEND,YSP_SG%YR%MP0)
        IF (IRSNML == 0) THEN
           PFPBUF1(IST:KEND,IFLD)=ZRSN(IST:KEND)
        ELSE
          IF (YSP_SG%YR%LSET .AND. IRSNML <= SIZE(PSP_SG,DIM=2)) THEN
            PFPBUF1(IST:KEND,IFLD)=PSP_SG(IST:KEND,IRSNML,YSP_SG%YR%MP0)
          ELSE
            CALL ABOR1('HPOS : ERROR IN SNOW DENSITy')
          ENDIF
        ENDIF
        IRSNML = IRSNML + 1
      ENDIF
    ELSEIF (IC == GFP%WSN%ICOD) THEN
      IF (LEPHYS) THEN
        IWSNML = IWSNML + 1
        IF (YSP_SG%YW%LSET .AND. IWSNML <= SIZE(PSP_SG,DIM=2)) THEN
          PFPBUF1(IST:KEND,IFLD)=MAX(0._JPRB,PSP_SG(IST:KEND,IWSNML,YSP_SG%YW%MP0))
        ELSE
          CALL ABOR1('HPOS : ERROR IN SNOW LIQUID WATER CONT')
        ENDIF
      ENDIF
    ELSEIF (IC == GFP%LICT%ICOD ) THEN
      !  Interpolation of lake ice-layer T.                                
      PFPBUF1(IST:KEND,IFLD)=PSP_SL(IST:KEND,YSP_SL%YLICT%MP0)     
    ELSEIF (IC == GFP%LTLT%ICOD ) THEN
      !  Interpolation of lake total-layer T.                     
      PFPBUF1(IST:KEND,IFLD)=PSP_SL(IST:KEND,YSP_SL%YLTLT%MP0)
    ELSEIF (IC == GFP%LMLT%ICOD ) THEN
      !  Interpolation of lake mixed-layer T.                     
      PFPBUF1(IST:KEND,IFLD)=PSP_SL(IST:KEND,YSP_SL%YLMLT%MP0)
    ELSEIF (IC == GFP%LBLT%ICOD ) THEN
      !  Interpolation of lake bottom-layer T.                      
      PFPBUF1(IST:KEND,IFLD)=PSP_SL(IST:KEND,YSP_SL%YLBLT%MP0)
    ELSEIF (IC == GFP%LSHF%ICOD ) THEN
      !  Interpolation of lake shape factor                                
      PFPBUF1(IST:KEND,IFLD)=PSP_SL(IST:KEND,YSP_SL%YLSHF%MP0)
    ELSEIF (IC == GFP%LICD%ICOD ) THEN
      !  Interpolation of lake ice depth                                
      PFPBUF1(IST:KEND,IFLD)=PSP_SL(IST:KEND,YSP_SL%YLICD%MP0)
    ELSEIF (IC == GFP%LMLD%ICOD ) THEN
      !  Interpolation of lake mixed-layer depth
      PFPBUF1(IST:KEND,IFLD)=PSP_SL(IST:KEND,YSP_SL%YLMLD%MP0)
    ELSEIF (IC == GFP%DST%ICOD .AND. YSP_SB%YT%LSET.AND.NCSS>=1) THEN
      IF (LMPHYS) THEN
        ZTD(IST:KEND)=PSP_SB(IST:KEND,1,YSP_SB%YT%MP0)
        IF (IANO==0) THEN
          ! Interpolation of Tp
          PFPBUF1(IST:KEND,IFLD)=ZTD(IST:KEND)
        ELSEIF(NFPCLI <= 2 .AND. YSP_RR%YT%LSET) THEN
          ! Interpolation of Tp-Ts
          PFPBUF1(IST:KEND,IFLD)=ZTD(IST:KEND)-PSP_RR(IST:KEND,YSP_RR%YT%MP0)
        ELSEIF (YSD_VX%YTPC%LSET.AND.YSD_VX%YORO%LSET) THEN
          ! Interpolation of Tp-Tp*
          PFPBUF1(IST:KEND,IFLD)=ZTD(IST:KEND)-PSD_VX(IST:KEND,YSD_VX%YTPC%MP)&
           & -ZDTSGDZ*(ZOROG(IST:KEND)-PSD_VX(IST:KEND,YSD_VX%YORO%MP))  
        ELSE
          CALL ABOR1('HPOS : SOMETHING IS MISSING FOR THE INTERPOLATION OF TD')
        ENDIF
      ELSEIF (LEPHYS) THEN
        IF (NFPSWI==0.AND.NCSS>=4) THEN
          ! Mean of TSL3 and TSL4
          ZTD(IST:KEND)=(PSP_SB(IST:KEND,3,YSP_SB%YT%MP0)+PSP_SB(IST:KEND,4,YSP_SB%YT%MP0))*0.5_JPRB
        ELSEIF(NCSS>=3) THEN
          ! Mean of TSL2 and TSL3
          ZTD(IST:KEND)=(PSP_SB(IST:KEND,2,YSP_SB%YT%MP0)+PSP_SB(IST:KEND,3,YSP_SB%YT%MP0))*0.5_JPRB
        ELSE
          CALL ABOR1('HPOS : NCSS TOO SMALL')
        ENDIF
        IF (IANO==0) THEN
          ! Interpolation of Tp
          PFPBUF1(IST:KEND,IFLD)=ZTD(IST:KEND)
        ELSE
          ! Interpolation of Tp-Ts with Ts=Tslev1
          PFPBUF1(IST:KEND,IFLD)=ZTD(IST:KEND)-PSP_SB(IST:KEND,1,YSP_SB%YT%MP0)
        ENDIF
      ENDIF
    ELSEIF (IC == GFP%SSW%ICOD) THEN
      IF (LMPHYS.AND.YSP_RR%YW%LSET) THEN
        IF (IANO == 0) THEN
          ! Interpolation of Ws
          DO JI=IST,KEND
            PFPBUF1(JI,IFLD)=PSP_RR(JI,YSP_RR%YW%MP0)
          ENDDO
        ELSE
          IF (YSD_VV%YARG%LSET.AND.YSD_VV%YD2%LSET.AND.YSD_VV%YSAB%LSET) THEN
            ! Compute Wsmax 
            LLHMT=.FALSE.
            IF (YSD_VX%YLSM%LSET.AND.YSD_VX%YIVEG%LSET.AND.&
               & .NOT.(YSD_VF%YLSM%LSET.AND.YSD_VV%YIVEG%LSET)) THEN
              ! Use extra land-sea mask and vegetation index (SURFEX case)
              CALL ACSOLW(YDPHY1,IST,KEND,NPROMA,PSD_VV(:,YSD_VV%YARG%MP),PSD_VV(:,YSD_VV%YD2%MP),&
               & PSD_VX(:,YSD_VX%YLSM%MP),PSD_VX(:,YSD_VX%YIVEG%MP),PSD_VV(:,YSD_VV%YSAB%MP),&
               & LLHMT,ZWFC,ZWPMX,ZWSAT,ZWSMX,ZWWILT)
            ELSEIF (YSD_VF%YLSM%LSET.AND.YSD_VV%YIVEG%LSET) THEN
              ! Current ISBA case
              CALL ACSOLW(YDPHY1,IST,KEND,NPROMA,PSD_VV(:,YSD_VV%YARG%MP),PSD_VV(:,YSD_VV%YD2%MP),&
               & PSD_VF(:,YSD_VF%YLSM%MP),PSD_VV(:,YSD_VV%YIVEG%MP),PSD_VV(:,YSD_VV%YSAB%MP),&
               & LLHMT,ZWFC,ZWPMX,ZWSAT,ZWSMX,ZWWILT)
            ELSE
              ZWSMX(IST:KEND)=PWSXI
            ENDIF
          ELSE
            ZWSMX(IST:KEND)=PWSXI
          ENDIF
          IF (NFPCLI <= 2) THEN
            ! Interpolation of relative value Ws/Wsmax
            PFPBUF1(IST:KEND,IFLD)=PSP_RR(IST:KEND,YSP_RR%YW%MP0)/ZWSMX(IST:KEND)
          ELSEIF (YSD_VX%YPWS%LSET) THEN
            ! Interpolation of the anomaly of Ws/Wsmax 
            PFPBUF1(IST:KEND,IFLD)=PSP_RR(IST:KEND,YSP_RR%YW%MP0)/ZWSMX(IST:KEND)&
             & -PSD_VX(IST:KEND,YSD_VX%YPWS%MP)  
          ELSE
            CALL ABOR1('HPOS : MISSING CLIM. RELATIVE SOIL WATER CONTENT')
          ENDIF
        ENDIF
      ELSEIF(LEPHYS.AND.YSP_SB%YQ%LSET.AND.YSD_VF%YSOTY%LSET.AND.NCSS >= 1) THEN
        IF (NFPSWI == 0) THEN
          ! Relative soil water content
          DO JI=IST,KEND
            ISOILTYPE=NINT(PSD_VF(JI,YSD_VF%YSOTY%MP))
            IF (ISOILTYPE==0) THEN
              PFPBUF1(JI,IFLD)=ZPSOL
            ELSE
              PFPBUF1(JI,IFLD)=PSP_SB(JI,1,YSP_SB%YQ%MP0)*ZPSOL/ZECWSAT(JI)
            ENDIF
          ENDDO
        ELSE
          ! Soil Wetness Index
          DO JI=IST,KEND
            ISOILTYPE=NINT(PSD_VF(JI,YSD_VF%YSOTY%MP))
            IF (ISOILTYPE==0) THEN
              PFPBUF1(JI,IFLD)=ZPSOL
            ELSE
              PFPBUF1(JI,IFLD)=PSP_SB(JI,1,YSP_SB%YQ%MP0)*ZPSOL/ZECWCAP(JI)
            ENDIF
          ENDDO
        ENDIF
      ENDIF
    ELSEIF (IC == GFP%DSW%ICOD) THEN
      IF (LMPHYS.AND.YSP_SB%YQ%LSET.AND.NCSS>=1) THEN
        IF (LLSCALWPI.AND.YSP_SB%YTL%LSET) THEN
          ZWPL(IST:KEND)=PSP_SB(IST:KEND,1,YSP_SB%YQ%MP0)+(1._JPRB-ZSCALWPI(IST:KEND))*PSP_SB(IST:KEND,1,YSP_SB%YTL%MP0)
        ELSE
          ZWPL(IST:KEND)=PSP_SB(IST:KEND,1,YSP_SB%YQ%MP0)
        ENDIF
        IF (LLSCALSNS.AND.YSP_SG%YF%LSET) THEN
          ZSNS(IST:KEND)=ZSSN(IST:KEND)*(MIN(1._JPRB,ZSCALSNS(IST:KEND) &
           & /(1._JPRB+ZSSN(IST:KEND)/RSNSMOD)))
          ! No correction on sea ice
          IF (YSD_VV%YIVEG%LSET) THEN
            DO JI=IST,KEND
              IF (NINT(PSD_VV(JI,YSD_VV%YIVEG%MP)) == NTVGLA ) THEN
                ZSNS(JI)=ZSSN(JI)
              ENDIF
            ENDDO
          ELSEIF (YSD_VX%YIVEG%LSET) THEN
            DO JI=IST,KEND
              IF (NINT(PSD_VX(JI,YSD_VX%YIVEG%MP)) == NTVGLA ) THEN
                ZSNS(JI)=ZSSN(JI)
              ENDIF
            ENDDO
          ENDIF
          ZWPL(IST:KEND)=ZWPL(IST:KEND) + ZSSN(IST:KEND) - ZSNS(IST:KEND)
        ENDIF
        IF (IANO == 0) THEN
          ! Interpolation of Wd
          PFPBUF1(IST:KEND,IFLD)=ZWPL(IST:KEND)
        ELSE
          ! Compute Wdmax
          IF (YSD_VV%YARG%LSET.AND.YSD_VV%YD2%LSET.AND.YSD_VV%YSAB%LSET) THEN
            LLHMT=.FALSE.
            IF (YSD_VX%YLSM%LSET.AND.YSD_VX%YIVEG%LSET.AND.&
               & .NOT.(YSD_VF%YLSM%LSET.AND.YSD_VV%YIVEG%LSET)) THEN
              ! Use extra land-sea mask and vegetation index (SURFEX case)
              CALL ACSOLW(YDPHY1,IST,KEND,NPROMA,PSD_VV(:,YSD_VV%YARG%MP),PSD_VV(:,YSD_VV%YD2%MP),&
               & PSD_VX(:,YSD_VX%YLSM%MP),PSD_VX(:,YSD_VX%YIVEG%MP),PSD_VV(:,YSD_VV%YSAB%MP),&
               & LLHMT,ZWFC,ZWPMX,ZWSAT,ZWSMX,ZWWILT)
            ELSEIF (YSD_VF%YLSM%LSET.AND.YSD_VV%YIVEG%LSET) THEN
              ! Current ISBA case
              CALL ACSOLW(YDPHY1,IST,KEND,NPROMA,PSD_VV(:,YSD_VV%YARG%MP),PSD_VV(:,YSD_VV%YD2%MP),&
               & PSD_VF(:,YSD_VF%YLSM%MP),PSD_VV(:,YSD_VV%YIVEG%MP),PSD_VV(:,YSD_VV%YSAB%MP),&
               & LLHMT,ZWFC,ZWPMX,ZWSAT,ZWSMX,ZWWILT)
            ELSE
              ZWPMX(IST:KEND)=PWDXI
            ENDIF
          ELSE
            ZWPMX(IST:KEND)=PWDXI
          ENDIF
          IF (NFPCLI <= 2) THEN
            ! Interpolation of relative value Wd/Wdmax
            PFPBUF1(IST:KEND,IFLD)=ZWPL(IST:KEND)/ZWPMX(IST:KEND)
          ELSEIF (YSD_VX%YPWP%LSET) THEN
            ! Interpolation of the anomaly of Wd/Wdmax
            PFPBUF1(IST:KEND,IFLD)=ZWPL(IST:KEND)/ZWPMX(IST:KEND)&
             & -PSD_VX(IST:KEND,YSD_VX%YPWP%MP)
          ELSE
            CALL ABOR1('HPOS : MISSING CLIM. RELATIVE DEEP SOIL WATER CONTENT')
          ENDIF
        ENDIF
      ELSEIF(LEPHYS.AND.YSP_SB%YQ%LSET.AND.YSP_SB%YT%LSET.AND.YSD_VF%YSOTY%LSET.AND.NCSS>=1) THEN
        IF (NFPSWI==2) THEN
          IF (YSD_VX%YARG%LSET.AND.YSD_VX%YXD2%LSET.AND.YSD_VX%YSAB%LSET.AND.&
            & YSD_VX%YRSMIN%LSET.AND. YSD_VX%YLAI%LSET.AND.YSD_VX%YVEG%LSET.AND.&
            & YSD_VF%YLSM%LSET) THEN
            ! "binary" index to compute the soil water/ice content from TESSEL to ISBA
            DO JI=IST,KEND
              IF (PSD_VF(JI,YSD_VF%YLSM%MP) > SMASK) THEN
                ZIVG(JI)=NTPDES
              ELSE
                ZIVG(JI)=NTPMER
              ENDIF
            ENDDO
            LLHMT=.FALSE.
            CALL ACSOLW(YDPHY1,IST,KEND,NPROMA,PSD_VX(:,YSD_VX%YARG%MP),PSD_VX(:,YSD_VX%YXD2%MP),&
             & PSD_VF(:,YSD_VF%YLSM%MP),ZIVG(:),PSD_VX(:,YSD_VX%YSAB%MP),LLHMT,&
             & ZWFC,ZWPMX,ZWSAT,ZWSMX,ZWWILT)
          ELSE
            CALL ABOR1('HPOS : MISSING FIELDS FOR LAI SCALING')
          ENDIF
        ENDIF
        DO JI=IST,KEND
          IF (NFPSWI == 0) THEN
            ! Ponderated mean value of relative deep soil water content
            ISOILTYPE=NINT(PSD_VF(JI,YSD_VF%YSOTY%MP))
            IF (ISOILTYPE==0) THEN
              PFPBUF1(JI,IFLD)=ZPSOL
            ELSE
              DO JCSS = 1,NCSS
                ZFSWL(JCSS) = PSP_SB(JI,JCSS,YSP_SB%YQ%MP0)*ZRDAT(JCSS)*ZPSOL
              ENDDO
              PFPBUF1(JI,IFLD)=SUM(ZFSWL(:))/ZECWSAT(JI)/SUM(ZRDAT(1:NCSS))
            ENDIF
          ELSE
            IF (YSD_VF%YCVL%LSET.AND.YSD_VF%YCVH%LSET.AND.&
             &  YSD_VF%YTVL%LSET.AND.YSD_VF%YTVH%LSET.AND.YSD_VF%YLSM%LSET) THEN
              IF (PSD_VF(JI,YSD_VF%YLSM%MP) > SMASK) THEN
                DO JCSS = 1,NCSS
                  IF (PSP_SB(JI,JCSS,YSP_SB%YT%MP0) > (TSRESERV2+TDELTA2)) THEN
                    ZICEL(JCSS) = 0.0_JPRB
                  ELSEIF (PSP_SB(JI,JCSS,YSP_SB%YT%MP0) < (TSRESERV2-TDELTA2)) THEN
                    ZICEL(JCSS)=1.0_JPRB
                  ELSE
                    ZDELTA=RPI*(PSP_SB(JI,JCSS,YSP_SB%YT%MP0)-TSRESERV2)/(2.0_JPRB*TDELTA2)
                    ZICEL(JCSS) = 0.5_JPRB * (1.0_JPRB - SIN(ZDELTA))
                  ENDIF
                ENDDO
                ITVH=NINT(PSD_VF(JI,YSD_VF%YTVH%MP))
                ITVL=NINT(PSD_VF(JI,YSD_VF%YTVL%MP))
                DO JCSS = 1,NCSS
                  ! Liquid water on each layer
                  ZFSWL(JCSS) = (1.0_JPRB-ZICEL(JCSS))*PSP_SB(JI,JCSS,YSP_SB%YQ%MP0)*ZPSOL
                  ! Pseudo evaporation of vegetation (high and low)
                  ZTETA_H(JCSS)=MAX(ZFSWL(JCSS),ZECWWILT(JI))*ZRVROOTSA(JCSS,ITVH)
                  ZTETA_L(JCSS)=MAX(ZFSWL(JCSS),ZECWWILT(JI))*ZRVROOTSA(JCSS,ITVL)
                ENDDO
                ZTETA_MEAN_H=SUM(ZTETA_H(:))
                ZTETA_MEAN_L=SUM(ZTETA_L(:))
                IF (ZTETA_MEAN_H > ZECWCAP(JI)) THEN
                  ZEVH = 1.0_JPRB
                ELSEIF (ZTETA_MEAN_H > ZECWWILT(JI)) THEN
                  ZEVH = (ZTETA_MEAN_H - ZECWWILT(JI))/(ZECWCAP(JI) - ZECWWILT(JI))
                ELSE
                  ZEVH = 0.0_JPRB
                ENDIF
                IF (ZTETA_MEAN_L > ZECWCAP(JI)) THEN
                  ZEVL = 1.0_JPRB
                ELSEIF (ZTETA_MEAN_L > ZECWWILT(JI)) THEN
                  ZEVL = (ZTETA_MEAN_L - ZECWWILT(JI))/(ZECWCAP(JI) - ZECWWILT(JI))
                ELSE
                  ZEVL = 0.0_JPRB
                ENDIF
                IF (ZFSWL(1) > ZECWCAP(JI)) THEN
                  ZECSWI = 1.0_JPRB
                ELSEIF (ZFSWL(1) > ZECWWILT(JI)) THEN
                  ZECSWI= ((ZFSWL(1))-ZECWWILT(JI))/(ZECWCAP(JI)-ZECWWILT(JI))
                ELSE
                  ZECSWI = 0.0_JPRB
                ENDIF
                ! Total pseudo evaporation IFS
                ! Land-sea mask is considered in this formula as land proportion, used to take into account the bare soil
                ZVEGH=PSD_VF(JI,YSD_VF%YCVH%MP)*ZRVCOV(ITVH)
                ZVEGL=PSD_VF(JI,YSD_VF%YCVL%MP)*ZRVCOV(ITVL)
                ZEVIFS = (ZVEGH*ZEVH+ZVEGL*ZEVL+(PSD_VF(JI,YSD_VF%YLSM%MP)-ZVEGH-ZVEGL)*ZECSWI)/PSD_VF(JI,YSD_VF%YLSM%MP)
                IF (NFPSWI == 1) THEN
                  PFPBUF1(JI,IFLD)=ZEVIFS
                ELSE
                  IF ( PSD_VX(JI,YSD_VX%YLAI%MP)==0._JPRB .OR. PSD_VX(JI,YSD_VX%YRSMIN%MP)==0._JPRB .OR.&
                   &   ZRVLAI(ITVH)==0._JPRB  .OR. ZRVLAI(ITVL)==0._JPRB .OR.&
                   &   ZVEGH+ZVEGL == 0._JPRB) THEN
                     ZLAISCAL = 1.0_JPRB
                   ELSE
                     ZRSMIN_LAI_H = ZRVRSMIN(ITVH)/ZRVLAI(ITVH)
                     ZRSMIN_LAI_L = ZRVRSMIN(ITVL)/ZRVLAI(ITVL)
                     ZRSMIN_LAI_IFS = (ZRSMIN_LAI_H * ZVEGH + ZRSMIN_LAI_L * ZVEGL ) / (ZVEGH + ZVEGL)
                     ZRSMIN_LAI_ISBA = PSD_VX(JI,YSD_VX%YRSMIN%MP) / PSD_VX(JI,YSD_VX%YLAI%MP)
                     ZLAISCAL = ZRSMIN_LAI_IFS / ZRSMIN_LAI_ISBA
                   ENDIF
                   ZVEG_ISBA=PSD_VX(JI,YSD_VX%YVEG%MP)
                   PFPBUF1(JI,IFLD)= (ZEVIFS*(ZWFC(JI)-ZWWILT(JI)) + ZLAISCAL*ZVEG_ISBA*ZWWILT(JI))*ZWFC(JI) /&
                    &( ZVEG_ISBA*ZWFC(JI)*(ZLAISCAL-1._JPRB) + ZWFC(JI) + (ZVEG_ISBA-1._JPRB)*ZWWILT(JI)  )
                ENDIF
              ELSE
                PFPBUF1(JI,IFLD)=1._JPRB
              ENDIF
            ELSE
              CALL ABOR1('HPOS : MISSING FIELDS IF YCVL OR YCVH OR YTVL OR YTVH OR YLSM')
            ENDIF
          ENDIF
        ENDDO
      ENDIF
    ELSEIF (IC == GFP%FSSW%ICOD) THEN
      IF (LMPHYS .AND. YSP_RR%YIC%LSET) THEN
        IF (IANO == 0) THEN
        ! Interpolation of Wsf
          PFPBUF1(IST:KEND,IFLD)=PSP_RR(IST:KEND,YSP_RR%YIC%MP0)
        ELSE
          ! Interpolation of Wsf/Wsmax
          IF (YSD_VV%YARG%LSET.AND.YSD_VV%YD2%LSET.AND.YSD_VV%YSAB%LSET) THEN
            LLHMT=.FALSE.
            IF (YSD_VX%YLSM%LSET.AND.YSD_VX%YIVEG%LSET.AND.&
               & .NOT.(YSD_VF%YLSM%LSET.AND.YSD_VV%YIVEG%LSET)) THEN
              ! Use extra land-sea mask and vegetation index (SURFEX case)
              CALL ACSOLW(YDPHY1,IST,KEND,NPROMA,PSD_VV(:,YSD_VV%YARG%MP),PSD_VV(:,YSD_VV%YD2%MP),&
               & PSD_VX(:,YSD_VX%YLSM%MP),PSD_VX(:,YSD_VX%YIVEG%MP),PSD_VV(:,YSD_VV%YSAB%MP),&
               & LLHMT,ZWFC,ZWPMX,ZWSAT,ZWSMX,ZWWILT)
            ELSEIF (YSD_VF%YLSM%LSET.AND.YSD_VV%YIVEG%LSET) THEN
              ! Current ISBA case
              CALL ACSOLW(YDPHY1,IST,KEND,NPROMA,PSD_VV(:,YSD_VV%YARG%MP),PSD_VV(:,YSD_VV%YD2%MP),&
               & PSD_VF(:,YSD_VF%YLSM%MP),PSD_VV(:,YSD_VV%YIVEG%MP),PSD_VV(:,YSD_VV%YSAB%MP),&
               & LLHMT,ZWFC,ZWPMX,ZWSAT,ZWSMX,ZWWILT)
            ELSE
              ZWSMX(IST:KEND)=PWSXI
            ENDIF
          ELSE
            ZWSMX(IST:KEND)=PWSXI
          ENDIF
          PFPBUF1(IST:KEND,IFLD)=PSP_RR(IST:KEND,YSP_RR%YIC%MP0)/ZWSMX(IST:KEND)
        ENDIF
      ELSEIF(LEPHYS.AND.YSP_SB%YQ%LSET.AND.NCSS >= 1) THEN
        IF (NFPSWI == 0) THEN
          ! Computed from Relative soil water content
          DO JI=IST,KEND
            ISOILTYPE=NINT(PSD_VF(JI,YSD_VF%YSOTY%MP))
            IF (ISOILTYPE==0) THEN
              PFPBUF1(JI,IFLD)=ZPSOL
            ELSE
              PFPBUF1(JI,IFLD)=PSP_SB(JI,1,YSP_SB%YQ%MP0)*ZPSOL/ZECWSAT(JI)
            ENDIF
          ENDDO
        ELSE
          ! Soil Wetness Index
          DO JI=IST,KEND
            ISOILTYPE=NINT(PSD_VF(JI,YSD_VF%YSOTY%MP))
            IF (ISOILTYPE==0) THEN
              PFPBUF1(JI,IFLD)=ZPSOL
            ELSE
              PFPBUF1(JI,IFLD)=PSP_SB(JI,1,YSP_SB%YQ%MP0)*ZPSOL/ZECWCAP(JI)
            ENDIF
          ENDDO
        ENDIF
      ELSE
        ! No surface frost
        PFPBUF1(IST:KEND,IFLD)=0._JPRB
      ENDIF
    ELSEIF (IC == GFP%FDSW%ICOD) THEN
      IF (LMPHYS.AND.YSP_SB%YTL%LSET.AND.NCSS>=1) THEN
        IF (LLSCALWPI) THEN
          ZWPI(IST:KEND)=PSP_SB(IST:KEND,1,YSP_SB%YTL%MP0)*ZSCALWPI(IST:KEND)
        ELSE
          ZWPI(IST:KEND)=PSP_SB(IST:KEND,1,YSP_SB%YTL%MP0)
        ENDIF
        IF (IANO == 0) THEN
          ! Interpolation of Wdf
          PFPBUF1(IST:KEND,IFLD)=ZWPI(IST:KEND)
        ELSE
          ! Interpolation of Wdf/Wdmax
          IF (YSD_VV%YARG%LSET.AND.YSD_VV%YD2%LSET.AND.YSD_VV%YSAB%LSET) THEN
            LLHMT=.FALSE.
            IF (YSD_VX%YLSM%LSET.AND.YSD_VX%YIVEG%LSET.AND.&
               & .NOT.(YSD_VF%YLSM%LSET.AND.YSD_VV%YIVEG%LSET)) THEN
              ! Use extra land-sea mask and vegetation index (SURFEX case)
              CALL ACSOLW(YDPHY1,IST,KEND,NPROMA,PSD_VV(:,YSD_VV%YARG%MP),PSD_VV(:,YSD_VV%YD2%MP),&
               & PSD_VX(:,YSD_VX%YLSM%MP),PSD_VX(:,YSD_VX%YIVEG%MP),PSD_VV(:,YSD_VV%YSAB%MP),&
               & LLHMT,ZWFC,ZWPMX,ZWSAT,ZWSMX,ZWWILT)
            ELSEIF (YSD_VF%YLSM%LSET.AND.YSD_VV%YIVEG%LSET) THEN
              ! Current ISBA case
              CALL ACSOLW(YDPHY1,IST,KEND,NPROMA,PSD_VV(:,YSD_VV%YARG%MP),PSD_VV(:,YSD_VV%YD2%MP),&
               & PSD_VF(:,YSD_VF%YLSM%MP),PSD_VV(:,YSD_VV%YIVEG%MP),PSD_VV(:,YSD_VV%YSAB%MP),&
               & LLHMT,ZWFC,ZWPMX,ZWSAT,ZWSMX,ZWWILT)
            ELSE
              ZWPMX(IST:KEND)=PWDXI
            ENDIF
          ELSE
            ZWPMX(IST:KEND)=PWDXI
          ENDIF
          PFPBUF1(IST:KEND,IFLD)=ZWPI(IST:KEND)/ZWPMX(IST:KEND)
        ENDIF
      ELSEIF(LEPHYS.AND.YSP_SB%YQ%LSET.AND.YSD_VF%YSOTY%LSET.AND.NCSS>=1) THEN
        DO JI=IST,KEND
          IF (NFPSWI == 0) THEN
          ! Computed from Relative deep soil water content
            ISOILTYPE=NINT(PSD_VF(JI,YSD_VF%YSOTY%MP))
            IF (ISOILTYPE==0) THEN
              PFPBUF1(JI,IFLD)=ZPSOL
            ELSE
              DO JCSS = 1,NCSS
                ZFSWL(JCSS) = PSP_SB(JI,JCSS,YSP_SB%YQ%MP0)*ZRDAT(JCSS)*ZPSOL
              ENDDO
              PFPBUF1(JI,IFLD)=SUM(ZFSWL(:))/ZECWSAT(JI)/SUM(ZRDAT(1:NCSS))
            ENDIF
          ELSE
            IF (YSP_SB%YT%LSET) THEN
              DO JCSS = 1,NCSS
                IF (PSP_SB(JI,JCSS,YSP_SB%YT%MP0) > (TSRESERV2+TDELTA2)) THEN
                  ZICEL(JCSS) = 0.0_JPRB
                ELSEIF (PSP_SB(JI,JCSS,YSP_SB%YT%MP0) < (TSRESERV2-TDELTA2)) THEN
                  ZICEL(JCSS)=1.0_JPRB
                ELSE
                  ZDELTA=RPI*(PSP_SB(JI,JCSS,YSP_SB%YT%MP0)-TSRESERV2)/(2.0_JPRB*TDELTA2)
                  ZICEL(JCSS) = 0.5_JPRB * (1.0_JPRB - SIN(ZDELTA))
                ENDIF
              ENDDO
              ! Ponderated mean frost
              DO JCSS = 1,NCSS
                ZFSWL(JCSS)=ZICEL(JCSS)*PSP_SB(JI,JCSS,YSP_SB%YQ%MP0)*ZRDAT(JCSS)
              ENDDO
              PFPBUF1(JI,IFLD)=SUM(ZFSWL(:))/SUM(ZRDAT(1:NCSS))
            ELSE
              CALL ABOR1('HPOS : TS NEEDED')
            ENDIF
          ENDIF
        ENDDO
      ELSE
        ! No deep soil frost
        PFPBUF1(IST:KEND,IFLD)=0._JPRB
      ENDIF
    ELSEIF (IC == GFP%CSSW%ICOD) THEN
      IF (YSD_VX%YPWS%LSET) THEN
        PFPBUF1(IST:KEND,IFLD)=PSD_VX(IST:KEND,YSD_VX%YPWS%MP)
      ELSEIF(IFPSEARCH(YDRQCLI,GFP%CSSW%ICOD) <= 0) THEN
        CALL ABOR1('HPOS : MISSING SURFPROP.RMAX.EAU IN INPUT/OUTPUT DATA')
      ENDIF
    ELSEIF (IC == GFP%CDSW%ICOD) THEN
      IF (YSD_VX%YPWP%LSET) THEN
        PFPBUF1(IST:KEND,IFLD)=PSD_VX(IST:KEND,YSD_VX%YPWP%MP)
      ELSEIF(IFPSEARCH(YDRQCLI,GFP%CDSW%ICOD) <= 0) THEN
        CALL ABOR1('HPOS : MISSING PROFPROP.RMAX.EAU IN INPUT/OUTPUT DATA')
      ENDIF
    ELSEIF (IC == GFP%SD%ICOD.AND.YSP_SG%YF%LSET) THEN
        IF (LMPHYS) THEN
        ! Interpolation of absolute snow depth
          PFPBUF1(IST:KEND,IFLD)=ZSSN(IST:KEND)
        ELSEIF (LEPHYS) THEN
          ISDML = ISDML + 1
          IF (ISDML <= SIZE(PSP_SG,DIM=2)) THEN
            PFPBUF1(IST:KEND,IFLD)=MAX(0._JPRB,PSP_SG(IST:KEND,ISDML,YSP_SG%YF%MP0))
          ELSE
            CALL ABOR1('HPOS : ERROR IN SNOW DEPTH')
          ENDIF
        ENDIF
    ELSEIF (IC == GFP%SDSL%ICOD.AND.YSP_SG%YF%LSET) THEN
        IF (LLSCALSNS) THEN
          ZSNS(IST:KEND)=ZSSN(IST:KEND)*(MIN(1._JPRB,ZSCALSNS(IST:KEND)/(1._JPRB+ZSSN(IST:KEND)/RSNSMOD)))
      ! ZSNS(IST:KEND)=PSP_SG(IST:KEND,JK,YSP_SG%YF%MP0)*(MIN(1._JPRB,ZSCALSNS(IST:KEND)/(1._JPRB+ZSSN(IST:KEND)/RSNSMOD)))
        ! No correction on sea ice
          IF (YSD_VV%YIVEG%LSET) THEN
            DO JI=IST,KEND
              IF (NINT(PSD_VV(JI,YSD_VV%YIVEG%MP)) == NTVGLA ) THEN
                ZSNS(JI)=ZSSN(JI)
              ENDIF
            ENDDO
          ELSEIF (YSD_VX%YIVEG%LSET) THEN
            DO JI=IST,KEND
              IF (NINT(PSD_VX(JI,YSD_VX%YIVEG%MP)) == NTVGLA ) THEN
                ZSNS(JI)=ZSSN(JI)
              ENDIF
            ENDDO
          ENDIF
        ELSE
          ZSNS(IST:KEND)=ZSSN(IST:KEND)
        ENDIF
        IF (IANO == 0.OR.NFPCLI <=2) THEN
          ! Interpolation of absolute snow depth
            PFPBUF1(IST:KEND,IFLD)=ZSNS(IST:KEND)
        ELSEIF (YSD_VX%YSNO%LSET) THEN
          ! Interpolation of snow depth anomaly.. to be re eval with snow ml
          PFPBUF1(IST:KEND,IFLD)=ZSNS(IST:KEND)-PSD_VX(IST:KEND,YSD_VX%YSNO%MP)
        ELSE
          CALL ABOR1('HPOS : CLIM. SNOW DEPTH NEEDED')
        ENDIF

    ELSEIF (IC == GFP%SR%ICOD  ) THEN
      IF (YSD_VF%YZ0F%LSET) THEN
        PFPBUF1(IST:KEND,IFLD)=PSD_VF(IST:KEND,YSD_VF%YZ0F%MP)
      ELSEIF (IFPSEARCH(YDRQCLI,GFP%SR%ICOD) <= 0) THEN
        CALL ABOR1('HPOS : MISSING SURFZ0.FOIS.G IN INPUT/OUTPUT DATA')
      ENDIF
    ELSEIF (IC == GFP%IDZ0%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VF(IST:KEND,YSD_VF%YZ0F%MP)
    ELSEIF (IC == GFP%ITZ0%ICOD) THEN
      IF (YSD_VV%YZ0H%LSET) THEN
        PFPBUF1(IST:KEND,IFLD)=PSD_VV(IST:KEND,YSD_VV%YZ0H%MP)
      ELSE
        CALL ABOR1('HPOS : MISSING THERMAL ROUGHNESS LENGTH IN INPUT DATA')
      ENDIF
    ELSEIF (IC == GFP%Z0H%ICOD) THEN
      IF (YSD_VV%YZ0H%LSET) THEN
        PFPBUF1(IST:KEND,IFLD)=PSD_VV(IST:KEND,YSD_VV%YZ0H%MP)
      ELSEIF (IFPSEARCH(YDRQCLI,GFP%Z0H%ICOD) <= 0) THEN
        CALL ABOR1('HPOS : MISSING THERMAL ROUGHNESS LENGTH IN INPUT/OUTPUT DATA')
      ENDIF
    ELSEIF (IC == GFP%AL%ICOD  ) THEN
      IF (YSD_VF%YALBF%LSET) THEN
        PFPBUF1(IST:KEND,IFLD)=PSD_VF(IST:KEND,YSD_VF%YALBF%MP)
      ELSEIF (IFPSEARCH(YDRQCLI,GFP%AL%ICOD) <= 0) THEN
        CALL ABOR1('HPOS : MISSING SURFALBEDO IN INPUT/OUTPUT DATA')
      ENDIF

!-- 4-component MODIS albedo
    ELSEIF (IC == GFP%ALUVP%ICOD  ) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VF(IST:KEND,YSD_VF%YALUVP%MP)
    ELSEIF (IC == GFP%ALUVD%ICOD  ) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VF(IST:KEND,YSD_VF%YALUVD%MP)
    ELSEIF (IC == GFP%ALNIP%ICOD  ) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VF(IST:KEND,YSD_VF%YALNIP%MP)
    ELSEIF (IC == GFP%ALNID%ICOD  ) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VF(IST:KEND,YSD_VF%YALNID%MP)

!-- 6-component MODIS albedo
    ELSEIF (IC == GFP%ALUVI%ICOD  ) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VF(IST:KEND,YSD_VF%YALUVI%MP)
    ELSEIF (IC == GFP%ALUVV%ICOD  ) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VF(IST:KEND,YSD_VF%YALUVV%MP)
    ELSEIF (IC == GFP%ALUVG%ICOD  ) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VF(IST:KEND,YSD_VF%YALUVG%MP)
    ELSEIF (IC == GFP%ALNII%ICOD  ) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VF(IST:KEND,YSD_VF%YALNII%MP)
    ELSEIF (IC == GFP%ALNIV%ICOD  ) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VF(IST:KEND,YSD_VF%YALNIV%MP)
    ELSEIF (IC == GFP%ALNIG%ICOD  ) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VF(IST:KEND,YSD_VF%YALNIG%MP)

      ! aerosol climatological fields
    ELSEIF (IC == GFP%AERDEP%ICOD  ) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VF(IST:KEND,YSD_VF%YAERDEP%MP)
    ELSEIF (IC == GFP%AERLTS%ICOD  ) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VF(IST:KEND,YSD_VF%YAERLTS%MP)
    ELSEIF (IC == GFP%AERSCC%ICOD  ) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VF(IST:KEND,YSD_VF%YAERSCC%MP)
    ELSEIF (IC == GFP%FCA1%ICOD    ) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VF(IST:KEND,YSD_VF%YFCA1%MP)
    ELSEIF (IC == GFP%FCA2%ICOD    ) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VF(IST:KEND,YSD_VF%YFCA2%MP)
    ELSEIF (IC == GFP%DSF%ICOD     ) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VF(IST:KEND,YSD_VF%YDSF%MP)

    ELSEIF (IC == GFP%DMSO%ICOD  ) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VF(IST:KEND,YSD_VF%YDMSO%MP)
    ELSEIF (IC == GFP%SO2DD%ICOD  ) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VF(IST:KEND,YSD_VF%YSO2DD%MP)

    ELSEIF (IC == GFP%SDFOR%ICOD  ) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VF(IST:KEND,YSD_VF%YSDFOR%MP)
    ELSEIF (IC == GFP%EMIS%ICOD) THEN
      IF (YSD_VF%YEMISF%LSET) THEN
        PFPBUF1(IST:KEND,IFLD)=PSD_VF(IST:KEND,YSD_VF%YEMISF%MP)
      ELSEIF (IFPSEARCH(YDRQCLI,GFP%EMIS%ICOD) <= 0) THEN
        CALL ABOR1('HPOS : MISSING SURFEMISSIVITY IN INPUT/OUTPUT DATA')
      ENDIF
! GPP/REC flux adjustment coefficient
    ELSEIF (IC == GFP%CGPP%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VF(IST:KEND,YSD_VF%YCGPP%MP)
    ELSEIF (IC == GFP%CREC%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VF(IST:KEND,YSD_VF%YCREC%MP)

      ! LAI variable fields 
    ELSEIF (IC == GFP%LAIH%ICOD  ) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VF(IST:KEND,YSD_VF%YLAIH%MP)
    ELSEIF (IC == GFP%LAIL%ICOD  ) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VF(IST:KEND,YSD_VF%YLAIL%MP)
      ! LAKE ancillary
    ELSEIF (IC == GFP%CLK%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VF(IST:KEND,YSD_VF%YCLK%MP)
    ELSEIF (IC == GFP%DL%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VF(IST:KEND,YSD_VF%YDL%MP)
      ! LAKE prognostic     
    ELSEIF (IC == GFP%LMLT%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSP_SL(IST:KEND,YSP_SL%YLMLT%MP)
    ELSEIF (IC == GFP%LMLD%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSP_SL(IST:KEND,YSP_SL%YLMLD%MP)
    ELSEIF (IC == GFP%LBLT%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSP_SL(IST:KEND,YSP_SL%YLBLT%MP)
    ELSEIF (IC == GFP%LTLT%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSP_SL(IST:KEND,YSP_SL%YLTLT%MP)
    ELSEIF (IC == GFP%LSHF%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSP_SL(IST:KEND,YSP_SL%YLSHF%MP)      
    ELSEIF (IC == GFP%LICT%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSP_SL(IST:KEND,YSP_SL%YLICT%MP)
    ELSEIF (IC == GFP%LICD%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSP_SL(IST:KEND,YSP_SL%YLICD%MP)

    ELSEIF (IC == GFP%SOTY%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VF(IST:KEND,YSD_VF%YSOTY%MP)
    ELSEIF (IC == GFP%SDOG%ICOD) THEN
      IF (YSD_VF%YGETRL%LSET) THEN
        PFPBUF1(IST:KEND,IFLD)=PSD_VF(IST:KEND,YSD_VF%YGETRL%MP)
      ELSEIF (IFPSEARCH(YDRQCLI,GFP%SDOG%ICOD) <= 0) THEN
        CALL ABOR1('HPOS : MISSING SURFET.GEOPOTENT IN INPUT/OUTPUT DATA')
      ENDIF
    ELSEIF (IC == GFP%VEG%ICOD) THEN
      IF (LMPHYS.AND.YSD_VF%YVEG%LSET) THEN
        ! MF field
        PFPBUF1(IST:KEND,IFLD)=PSD_VF(IST:KEND,YSD_VF%YVEG%MP)
      ELSEIF (LEPHYS .AND. YSD_VF%YCVL%LSET .AND. YSD_VF%YCVH%LSET .AND. YSD_VF%YTVL%LSET .AND. YSD_VF%YTVH%LSET) THEN
        ! Conversion from EC fields
        DO JI=IST,KEND
          PFPBUF1(JI,IFLD)=PSD_VF(JI,YSD_VF%YCVH%MP)*ZRVCOV(NINT(PSD_VF(JI,YSD_VF%YTVH%MP))) +&
           & (PSD_VF(JI,YSD_VF%YCVL%MP)*ZRVCOV(NINT(PSD_VF(JI,YSD_VF%YTVL%MP))))
        ENDDO
      ELSEIF (LMPHYS.AND.IFPSEARCH(YDRQCLI,GFP%VEG%ICOD) <= 0) THEN
        CALL ABOR1('HPOS : MISSING SURFPROP.VEGETAT IN INPUT/OUTPUT DATA')
      ENDIF
    ELSEIF (IC == GFP%ISOR%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VF(IST:KEND,YSD_VF%YVRLAN%MP)
    ELSEIF (IC == GFP%ANOR%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VF(IST:KEND,YSD_VF%YVRLDI%MP)
    ELSEIF (IC == GFP%SLOR%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VF(IST:KEND,YSD_VF%YSIG%MP)
    ELSEIF (IC == GFP%LSRH%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VF(IST:KEND,YSD_VF%YLZ0H%MP)
    ELSEIF (IC == GFP%GFIS%ICOD .OR. IC == GFP%SFIS%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=ZOROG(IST:KEND)
    ELSEIF (IC == GFP%ACOT%ICOD) THEN
      IF (YSD_VF%YVRLAN%LSET) THEN
        PFPBUF1(IST:KEND,IFLD)=PSD_VF(IST:KEND,YSD_VF%YVRLAN%MP)
      ELSEIF (IFPSEARCH(YDRQCLI,GFP%ACOT%ICOD) <= 0) THEN
        CALL ABOR1('HPOS : MISSING SURFVAR.GEOP.ANI IN INPUT/OUTPUT DATA')
      ENDIF
    ELSEIF (IC == GFP%DPAT%ICOD) THEN
      IF (YSD_VF%YVRLDI%LSET) THEN
        PFPBUF1(IST:KEND,IFLD)=PSD_VF(IST:KEND,YSD_VF%YVRLDI%MP)
      ELSEIF (IFPSEARCH(YDRQCLI,GFP%DPAT%ICOD) <= 0) THEN
        CALL ABOR1('HPOS : MISSING SURFVAR.GEOP.DIR IN INPUT/OUTPUT DATA')
      ENDIF
    ELSEIF (IC == GFP%IVEG%ICOD) THEN
      IF (YSD_VV%YIVEG%LSET) THEN
        PFPBUF1(IST:KEND,IFLD)=PSD_VV(IST:KEND,YSD_VV%YIVEG%MP)
      ELSEIF (IFPSEARCH(YDRQCLI,GFP%IVEG%ICOD) <= 0) THEN
        CALL ABOR1('HPOS : MISSING VEGETATION INDEX IN INPUT/OUTPUT DATA')
      ENDIF
    ELSEIF (IC == GFP%RSMIN%ICOD) THEN
      IF (YSD_VV%YRSMIN%LSET) THEN
        PFPBUF1(IST:KEND,IFLD)=PSD_VV(IST:KEND,YSD_VV%YRSMIN%MP)
      ELSEIF (IFPSEARCH(YDRQCLI,GFP%RSMIN%ICOD) <= 0) THEN
        CALL ABOR1('HPOS : MISSING MIN. STOMATAL RESISTANCE IN INPUT/OUTPUT DATA')
      ENDIF
    ELSEIF (IC == GFP%ARG%ICOD) THEN
      IF (YSD_VV%YARG%LSET) THEN
        PFPBUF1(IST:KEND,IFLD)=PSD_VV(IST:KEND,YSD_VV%YARG%MP)
      ELSEIF (IFPSEARCH(YDRQCLI,GFP%ARG%ICOD) <= 0) THEN
        CALL ABOR1('HPOS : MISSING CLAY PROPORTION IN INPUT/OUTPUT DATA') 
      ENDIF
    ELSEIF (IC == GFP%SAB%ICOD) THEN
      IF (YSD_VV%YSAB%LSET) THEN
        PFPBUF1(IST:KEND,IFLD)=PSD_VV(IST:KEND,YSD_VV%YSAB%MP)
      ELSEIF (IFPSEARCH(YDRQCLI,GFP%SAB%ICOD) <= 0) THEN
        CALL ABOR1('HPOS : MISSING PROPORTION OF SAND IN INPUT/OUTPUT DATA')
      ENDIF
    ELSEIF (IC == GFP%D2%ICOD) THEN
      IF (YSD_VV%YD2%LSET) THEN
        PFPBUF1(IST:KEND,IFLD)=PSD_VV(IST:KEND,YSD_VV%YD2%MP)
      ELSEIF (IFPSEARCH(YDRQCLI,GFP%D2%ICOD) <= 0) THEN
        CALL ABOR1('HPOS : MISSING SOIL DEPTH IN INPUT/OUTPUT DATA')
      ENDIF
    ELSEIF (IC == GFP%LAI%ICOD) THEN
      IF (YSD_VV%YLAI%LSET) THEN
        PFPBUF1(IST:KEND,IFLD)=PSD_VV(IST:KEND,YSD_VV%YLAI%MP)
      ELSEIF (IFPSEARCH(YDRQCLI,GFP%LAI%ICOD) <= 0) THEN
        CALL ABOR1('HPOS : MISSING LEAF AREA INDEX IN INPUT/OUTPUT DATA')
      ENDIF
    ELSEIF (IC == GFP%HV%ICOD) THEN
      IF (YSD_VV%YHV%LSET) THEN
        PFPBUF1(IST:KEND,IFLD)=PSD_VV(IST:KEND,YSD_VV%YHV%MP)
      ELSE
        ! Default value had better be 1. like over sea
        PFPBUF1(IST:KEND,IFLD)=1._JPRB
      ENDIF
    ELSEIF (IC == GFP%ALS%ICOD) THEN
      IF (YSD_VV%YALS%LSET) THEN
        PFPBUF1(IST:KEND,IFLD)=PSD_VV(IST:KEND,YSD_VV%YALS%MP)
      ELSEIF (IFPSEARCH(YDRQCLI,GFP%ALS%ICOD) <= 0) THEN
        CALL ABOR1('HPOS : MISSING ALBEDO OF BARE GROUND IN INPUT/OUTPUT DATA')
      ENDIF
    ELSEIF (IC == GFP%ALV%ICOD) THEN
      IF (YSD_VV%YALV%LSET) THEN
        PFPBUF1(IST:KEND,IFLD)=PSD_VV(IST:KEND,YSD_VV%YALV%MP)
      ELSEIF (IFPSEARCH(YDRQCLI,GFP%ALV%ICOD) <= 0) THEN
        CALL ABOR1('HPOS : MISSING ALBEDO OF VEGETATION IN INPUT/OUTPUT DATA')
      ENDIF
    ELSEIF (IC == GFP%IC%ICOD) THEN
      IF (YSP_RR%YFC%LSET) THEN
        PFPBUF1(IST:KEND,IFLD)=PSP_RR(IST:KEND,YSP_RR%YFC%MP0)
      ELSE
        PFPBUF1(IST:KEND,IFLD)=0._JPRB
      ENDIF
    ELSEIF (IC == GFP%ALSN%ICOD) THEN
      IF (YSP_SG%YA%LSET) THEN
!              EDSNOWML: only 1 layer handled here, double check 
        !PFPBUF1(IST:KEND,IFLD)=PSP_SG(IST:KEND,YSP_SG%YA%MP0)
        !PFPBUF1(IST:KEND,IFLD)=PSP_SG(IST:KEND,1,YSP_SG%YA%MP0)
        PFPBUF1(IST:KEND,IFLD)=PSP_SG(IST:KEND,1,YSP_SG%YSG(2)%MP0)
      ELSE
        CALL ABOR1('HPOS : MISSING SURFALBEDO NEIGE IN INPUT DATA')
      ENDIF
    ELSEIF (IC == GFP%SNDE%ICOD) THEN
      ! Snow density (adimensional : in kg/kg)
      IF (LMPHYS.AND.YSP_SG%YR%LSET) THEN
        PFPBUF1(IST:KEND,IFLD)=ZRSN(IST:KEND)
      ELSEIF (LEPHYS.AND.YSP_SG%YR%LSET) THEN
      ! Snow mass per unit of volume in kg.m-3
        PFPBUF1(IST:KEND,IFLD)=ZRSN(IST:KEND)*0.001_JPRB
      ELSEIF(LMPHYS) THEN
        CALL ABOR1('HPOS : MISSING SURFDENSIT.NEIGE IN INPUT DATA')
      ENDIF
    ELSEIF (IC == GFP%BSR%ICOD) THEN
      IF (YSD_VF%YZ0RLF%LSET) THEN
        PFPBUF1(IST:KEND,IFLD)=PSD_VF(IST:KEND,YSD_VF%YZ0RLF%MP)
      ELSEIF (IFPSEARCH(YDRQCLI,GFP%BSR%ICOD) <= 0) THEN
        CALL ABOR1('HPOS : MISSING REL. ROUGHNESS LENGTH IN INPUT/OUTPUT DATA')
      ENDIF
    ELSEIF (IC == GFP%BAAL%ICOD) THEN
      IF (YSD_VF%YALBSF%LSET) THEN
        PFPBUF1(IST:KEND,IFLD)=PSD_VF(IST:KEND,YSD_VF%YALBSF%MP)
      ELSEIF (IFPSEARCH(YDRQCLI,GFP%BAAL%ICOD) <= 0) THEN
        CALL ABOR1('HPOS : MISSING SOIL SHORTWAVE ALBEDO IN INPUT/OUTPUT DATA')
      ENDIF
    ELSEIF (IC == GFP%ALBHIS%ICOD) THEN
      IF (YSP_SG%YT%LSET) THEN
        ! Warning : ambiguous pointer usage.
        ! Here we are supposed to point to historical albedo, not snow temerature
        ! (see su_surf_flds). REK
        PFPBUF1(IST:KEND,IFLD)=ZTSN(IST:KEND)
      ELSEIF (IFPSEARCH(YDRQCLI,GFP%ALBHIS%ICOD) <= 0) THEN
        CALL ABOR1('HPOS : MISSING TOTAL ALBEDO IN INPUT/OUTPUT DATA')
      ENDIF
    ELSEIF (IC == GFP%ASEA%ICOD) THEN
      IF (YSD_VA%YSEA%LSET) THEN
        PFPBUF1(IST:KEND,IFLD)=PSD_VA(IST:KEND,YSD_VA%YSEA%MP)
      ELSEIF(IFPSEARCH(YDRQCLI,GFP%ASEA%ICOD) <= 0) THEN
        CALL ABOR1('HPOS : MISSING SURFAEROS.SEA IN INPUT/OUTPUT DATA')
      ENDIF
    ELSEIF (IC == GFP%ALAN%ICOD) THEN
      IF (YSD_VA%YLAN%LSET) THEN
        PFPBUF1(IST:KEND,IFLD)=PSD_VA(IST:KEND,YSD_VA%YLAN%MP)
      ELSEIF(IFPSEARCH(YDRQCLI,GFP%ALAN%ICOD) <= 0) THEN
        CALL ABOR1('HPOS : MISSING SURFAEROS.LAND IN INPUT/OUTPUT DATA')
      ENDIF
    ELSEIF (IC == GFP%ASOO%ICOD) THEN
      IF (YSD_VA%YSOO%LSET) THEN
        PFPBUF1(IST:KEND,IFLD)=PSD_VA(IST:KEND,YSD_VA%YSOO%MP)
      ELSEIF(IFPSEARCH(YDRQCLI,GFP%ASOO%ICOD) <= 0) THEN
        CALL ABOR1('HPOS : MISSING SURFAEROS.SOOT IN INPUT/OUTPUT DATA')
      ENDIF
    ELSEIF (IC == GFP%ADES%ICOD) THEN
      IF (YSD_VA%YDES%LSET) THEN
        PFPBUF1(IST:KEND,IFLD)=PSD_VA(IST:KEND,YSD_VA%YDES%MP)
      ELSEIF(IFPSEARCH(YDRQCLI,GFP%ADES%ICOD) <= 0) THEN
        CALL ABOR1('HPOS : MISSING SURFAEROS.DESERT IN INPUT/OUTPUT DATA')
      ENDIF
    ELSEIF (IC == GFP%O3A%ICOD) THEN
      IF (YSD_VC%YA%LSET) THEN
        PFPBUF1(IST:KEND,IFLD)=PSD_VC(IST:KEND,YSD_VC%YA%MP)
      ELSEIF (IFPSEARCH(YDRQCLI,GFP%O3A%ICOD) <= 0) THEN
        CALL ABOR1('HPOS : MISSING SURFA.OF.OZONE IN INPUT/OUTPUT DATA')
      ENDIF
    ELSEIF (IC == GFP%O3B%ICOD) THEN
      IF (YSD_VC%YB%LSET) THEN
        PFPBUF1(IST:KEND,IFLD)=PSD_VC(IST:KEND,YSD_VC%YB%MP)
      ELSEIF (IFPSEARCH(YDRQCLI,GFP%O3B%ICOD) <= 0) THEN
        CALL ABOR1('HPOS : MISSING SURFB.OF.OZONE IN INPUT/OUTPUT DATA')
      ENDIF
    ELSEIF (IC == GFP%O3C%ICOD) THEN
      IF (YSD_VC%YC%LSET) THEN
        PFPBUF1(IST:KEND,IFLD)=PSD_VC(IST:KEND,YSD_VC%YC%MP)
      ELSEIF (IFPSEARCH(YDRQCLI,GFP%O3C%ICOD) <= 0) THEN
        CALL ABOR1('HPOS : MISSING SURFC.OF.OZONE IN INPUT/OUTPUT DATA')
      ENDIF
    ELSEIF (IC == GFP%PADOU%ICOD) THEN
      IF (YSD_VF%YVRLAN%LSET.AND.YSD_VF%YVRLDI%LSET.AND.YSD_VF%YGETRL%LSET) THEN
        !     first component of vector
        !        (sigma**2)*(1-gamma)/(1+gamma) * (cos(2*alpha),sin(2*alpha))
        !     where sigma = standard deviation of orography
        !           gamma = Anisotropy coefficient of topography
        !           alpha = Direction of the principal axis of the topography
        LLU=.TRUE.
        LLV=.FALSE.
        ZGAMMA(IST:KEND)=PSD_VF(IST:KEND,YSD_VF%YVRLAN%MP)
        ZALFA(IST:KEND)=PSD_VF(IST:KEND,YSD_VF%YVRLDI%MP)
        CALL CVLANISO(IST,KEND,NPROMA,LLU,LLV,ZGAMMA(:),ZALFA(:),&
         & PFPBUF1(1,JFLD),ZDUMM,PSD_VF(:,YSD_VF%YGETRL%MP))  
        IF (LDFPOSHOR) THEN
        !       Simulate a "reduced" field :
          DO JI=IST,KEND
            PFPBUF1(JI,IFLD)=PFPBUF1(JI,IFLD)/ZGM(JI)
          ENDDO
        ENDIF
      ELSEIF (IFPSEARCH(YDRQCLI,GFP%ACOT%ICOD) <= 0 .AND. IFPSEARCH(YDRQCLI,GFP%DPAT%ICOD) <= 0 &
         & .AND. IFPSEARCH(YDRQCLI,GFP%SDOG%ICOD) <= 0) THEN
        CALL ABOR1('HPOS : MISSING SURFET.GEOPOTENT/SURFVAR.GEOP.ANI/SURFVAR.GEOP.DIR IN INPUT/OUTPUT DATA')
      ENDIF
    ELSEIF (IC == GFP%PADOV%ICOD) THEN
      IF (YSD_VF%YVRLAN%LSET.AND.YSD_VF%YVRLDI%LSET.AND.YSD_VF%YGETRL%LSET) THEN
        !     Second component of vector
        !        (sigma**2)*(1-gamma)/(1+gamma) * (cos(2*alpha),sin(2*alpha))
        !     where sigma = standard deviation of orography
        !           gamma = Anisotropy coefficient of topography
        !           alpha = Direction of the principal axis of the topography
        LLU=.FALSE.
        LLV=.TRUE.
        ZGAMMA(IST:KEND)=PSD_VF(IST:KEND,YSD_VF%YVRLAN%MP)
        ZALFA(IST:KEND)=PSD_VF(IST:KEND,YSD_VF%YVRLDI%MP)
        CALL CVLANISO(IST,KEND,NPROMA,LLU,LLV,ZGAMMA(:),ZALFA(:),&
         & ZDUMM,PFPBUF1(1,JFLD),PSD_VF(:,YSD_VF%YGETRL%MP))  
        IF (LDFPOSHOR) THEN
        !       Simulate a "reduced" field :
          DO JI=IST,KEND
            PFPBUF1(JI,IFLD)=PFPBUF1(JI,IFLD)/ZGM(JI)
          ENDDO
        ENDIF
      ELSEIF (IFPSEARCH(YDRQCLI,GFP%ACOT%ICOD) <= 0 .AND. IFPSEARCH(YDRQCLI,GFP%DPAT%ICOD) <= 0 .AND. &
       & IFPSEARCH(YDRQCLI,GFP%SDOG%ICOD) <= 0) THEN
        CALL ABOR1('HPOS : MISSING SURFET.GEOPOTENT/SURFVAR.GEOP.ANI/SURFVAR.GEOP.DIR IN INPUT/OUTPUT DATA')
      ENDIF
    ELSEIF (IC == GFP%PSRHU%ICOD .AND.YSP_RR%YW%LSET) THEN
      !     Relative surface moisture
      IF (NFPCLI <= 2) THEN
        Z1SWSX=1.0_JPRB/PWSXI
        PFPBUF1(IST:KEND,JFLD)=0.5_JPRB*&
         & (1.0_JPRB-COS(RPI*PSP_RR(IST:KEND,YSP_RR%YW%MP0)*Z1SWSX))  
      ELSE
        IF (YSD_VV%YARG%LSET) THEN
          LLHMT=.TRUE. ! no need for vegetation index, d2, sand
          IF (YSD_VX%YLSM%LSET.AND..NOT.YSD_VF%YLSM%LSET) THEN
            ! Use extra land-sea mask (SURFEX case)
            CALL ACSOLW(YDPHY1,IST,KEND,NPROMA,PSD_VV(:,YSD_VV%YARG%MP),ZDUMM,&
             & PSD_VX(:,YSD_VX%YLSM%MP),ZDUMM,ZDUMM,LLHMT,ZWFC,ZWPMX,ZWSAT,ZWSMX,ZWWILT)
            DO JI=IST,KEND
              Z1SWSX=1.0_JPRB/(ZWFC(JI)*GCONV*RD1)
              PFPBUF1(JI,JFLD)=0.5_JPRB*&
               & (1.0_JPRB-COS(RPI*MIN(PSP_RR(JI,YSP_RR%YW%MP0)*Z1SWSX,1.0_JPRB)))
            ENDDO
          ELSEIF (YSD_VF%YLSM%LSET) THEN
            ! Current ISBA case
            CALL ACSOLW(YDPHY1,IST,KEND,NPROMA,PSD_VV(:,YSD_VV%YARG%MP),ZDUMM,&
             & PSD_VF(:,YSD_VF%YLSM%MP),ZDUMM,ZDUMM,LLHMT,ZWFC,ZWPMX,ZWSAT,ZWSMX,ZWWILT)
            DO JI=IST,KEND
              Z1SWSX=1.0_JPRB/(ZWFC(JI)*GCONV*RD1)
              PFPBUF1(JI,JFLD)=0.5_JPRB*&
               & (1.0_JPRB-COS(RPI*MIN(PSP_RR(JI,YSP_RR%YW%MP0)*Z1SWSX,1.0_JPRB)))
            ENDDO
          ELSE
            Z1SWSX=1.0_JPRB/PWSXI
            PFPBUF1(IST:KEND,JFLD)=0.5_JPRB*&
             & (1.0_JPRB-COS(RPI*PSP_RR(IST:KEND,YSP_RR%YW%MP0)*Z1SWSX))  
          ENDIF
        ELSE
          Z1SWSX=1.0_JPRB/PWSXI
          PFPBUF1(IST:KEND,JFLD)=0.5_JPRB*&
           & (1.0_JPRB-COS(RPI*PSP_RR(IST:KEND,YSP_RR%YW%MP0)*Z1SWSX))  
        ENDIF
      ENDIF
    ELSEIF (IC == GFP%PCAPG%ICOD) THEN
      IF (YSP_CI%YCI(1)%LSET) THEN
        PFPBUF1(IST:KEND,IFLD)=PSP_CI(IST:KEND,YSP_CI%YCI(1)%MP0)
      ELSE
        CALL ABOR1('HPOS : MISSING SURFETP.GEOPOTEN IN INPUT DATA')
      ENDIF
    ELSEIF (IC == GFP%PCAAG%ICOD) THEN
      IF (YSP_CI%YCI(2)%LSET) THEN
        PFPBUF1(IST:KEND,IFLD)=PSP_CI(IST:KEND,YSP_CI%YCI(2)%MP0)
      ELSE
        CALL ABOR1('HPOS : MISSING SURFETA.GEOPOTEN IN INPUT DATA')
      ENDIF

      !        add the ECMWF diagnostics/physics fields (CNLSP to CNTCO3).

    ELSEIF (IC == GFP%LSP%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YLSP%MP)
    ELSEIF (IC == GFP%CP%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YCP%MP)
      !           Total precipitation
    ELSEIF (IC == GFP%TP%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YLSP%MP)+PSD_VD(IST:KEND,YSD_VD%YCP%MP)
    ELSEIF (IC == GFP%SF%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YSF%MP)
    ELSEIF (IC == GFP%FZRA%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YFZRA%MP)
    ELSEIF (IC == GFP%BLD%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YBLD%MP)
    ELSEIF (IC == GFP%SSHF%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YSSHF%MP)
    ELSEIF (IC == GFP%SLHF%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YSLHF%MP)
    ELSEIF (IC == GFP%MSLD%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YMSL%MP)
    ELSEIF (IC == GFP%SP%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YSP%MP)
    ELSEIF (IC == GFP%TCC%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YTCC%MP)
    ELSEIF (IC == GFP%EC10U%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%Y10U%MP)
    ELSEIF (IC == GFP%EC10V%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%Y10V%MP)
    ELSEIF (IC == GFP%EC10SI%ICOD) THEN
      ZU(IST:KEND) = PSD_VD(IST:KEND,YSD_VD%Y10U%MP)
      ZV(IST:KEND) = PSD_VD(IST:KEND,YSD_VD%Y10V%MP)
      PFPBUF1(IST:KEND,IFLD)=SQRT(ZU(IST:KEND)*ZU(IST:KEND) + ZV(IST:KEND)*ZV(IST:KEND))
    ELSEIF (IC == GFP%EC2T%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%Y2T%MP)
    ELSEIF (IC == GFP%EC2D%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%Y2D%MP)
    ELSEIF (IC == GFP%EC2SH%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%Y2SH%MP)
    ELSEIF (IC == GFP%SSR%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YSSR%MP)
    ELSEIF (IC == GFP%STR%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YSTR%MP)
    ELSEIF (IC == GFP%TSR%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YTSR%MP)
    ELSEIF (IC == GFP%TTR%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YTTR%MP)
    ELSEIF (IC == GFP%EWSS%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YEWSS%MP)
    ELSEIF (IC == GFP%NSSS%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YNSSS%MP)
    ELSEIF (IC == GFP%E%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YE%MP)
    ELSEIF (IC == GFP%PEV%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YPEV%MP)
    ELSEIF (IC == GFP%CCC%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YCCC%MP)
    ELSEIF (IC == GFP%LCC%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YLCC%MP)
    ELSEIF (IC == GFP%MCC%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YMCC%MP)
    ELSEIF (IC == GFP%HCC%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YHCC%MP)
    ELSEIF (IC == GFP%LGWS%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YLGWS%MP)
    ELSEIF (IC == GFP%MGWS%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YMGWS%MP)
    ELSEIF (IC == GFP%GWD%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YGWD%MP)
    ELSEIF (IC == GFP%MX2T%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YMX2T%MP)
    ELSEIF (IC == GFP%MN2T%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YMN2T%MP)
    ELSEIF (IC == GFP%MX2T3%ICOD) THEN
      I3DIM = SIZE(YSD_VD%YMX2T6)/2
      PFPBUF1(IST:KEND,IFLD)=MAXVAL(PSD_VD(IST:KEND,YSD_VD%YMX2T6(1:I3DIM)%MP),DIM=2)
    ELSEIF (IC == GFP%MN2T3%ICOD) THEN
      I3DIM = SIZE(YSD_VD%YMN2T6)/2
      PFPBUF1(IST:KEND,IFLD)=MINVAL(PSD_VD(IST:KEND,YSD_VD%YMN2T6(1:I3DIM)%MP),DIM=2)
    ELSEIF (IC == GFP%MX2T6%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=MAXVAL(PSD_VD(IST:KEND,YSD_VD%YMX2T6(:)%MP),DIM=2)
    ELSEIF (IC == GFP%MN2T6%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=MINVAL(PSD_VD(IST:KEND,YSD_VD%YMN2T6(:)%MP),DIM=2)
    ELSEIF (IC == GFP%RO%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YRO%MP)
    ELSEIF (IC == GFP%SRO%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YSRO%MP)    
    ELSEIF (IC == GFP%SSRO%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YSSRO%MP)  
    ELSEIF (IC == GFP%NEE%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YNEE%MP)
    ELSEIF (IC == GFP%GPP%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YGPP%MP)
    ELSEIF (IC == GFP%REC%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YREC%MP)
    ELSEIF (IC == GFP%INEE%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YINEE%MP)
    ELSEIF (IC == GFP%IGPP%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YIGPP%MP)
    ELSEIF (IC == GFP%IREC%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YIREC%MP)
    ELSEIF (IC == GFP%ALB%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YALB%MP)
    ELSEIF (IC == GFP%IEWSS%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YIEWSS%MP)
    ELSEIF (IC == GFP%INSSS%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YINSSS%MP)
    ELSEIF (IC == GFP%ISSHF%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YISSHF%MP)
    ELSEIF (IC == GFP%IE%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YIE%MP)
    ELSEIF (IC == GFP%CSF%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YCSF%MP)
    ELSEIF (IC == GFP%LSF%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YLSSF%MP)
    ELSEIF (IC == GFP%MXTPR%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YMXTPR%MP)
    ELSEIF (IC == GFP%MNTPR%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YMNTPR%MP)
    ELSEIF (IC == GFP%MXTPR3%ICOD) THEN
      I3DIM = SIZE(YSD_VD%YMXTPR6)/2
      PFPBUF1(IST:KEND,IFLD)=MAXVAL(PSD_VD(IST:KEND,YSD_VD%YMXTPR6(1:I3DIM)%MP),DIM=2)
    ELSEIF (IC == GFP%MNTPR3%ICOD) THEN
      I3DIM = SIZE(YSD_VD%YMNTPR6)/2
      PFPBUF1(IST:KEND,IFLD)=MINVAL(PSD_VD(IST:KEND,YSD_VD%YMNTPR6(1:I3DIM)%MP),DIM=2)
    ELSEIF (IC == GFP%MXTPR6%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=MAXVAL(PSD_VD(IST:KEND,YSD_VD%YMXTPR6(:)%MP),DIM=2)
    ELSEIF (IC == GFP%MNTPR6%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=MINVAL(PSD_VD(IST:KEND,YSD_VD%YMNTPR6(:)%MP),DIM=2)
    ELSEIF (IC == GFP%TPR%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YTPR%MP)
    ELSEIF (IC == GFP%LSRR%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YLSRR%MP)
    ELSEIF (IC == GFP%CRR%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YCRR%MP)
    ELSEIF (IC == GFP%LSSFR%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YLSSFR%MP)
    ELSEIF (IC == GFP%CSFR%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YCSFR%MP)
    ELSEIF (IC == GFP%PTYPE%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YPTYPE%MP)
    ELSEIF (IC == GFP%ILSPF%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YILSPF%MP)
    ELSEIF (IC == GFP%Z0F%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YZ0F%MP)
    ELSEIF (IC == GFP%LZ0H%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YLZ0H%MP)
    ELSEIF (IC == GFP%VIWVE%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YVIWVE%MP)
    ELSEIF (IC == GFP%VIWVN%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YVIWVN%MP)
    ELSEIF (IC == GFP%TCW%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YTCW%MP)
    ELSEIF (IC == GFP%TCWV%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YTCWV%MP)
    ELSEIF (IC == GFP%TCLW%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YTCLW%MP)
    ELSEIF (IC == GFP%TCIW%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YTCIW%MP)
    ELSEIF (IC == GFP%TCRW%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YTCRW%MP)
    ELSEIF (IC == GFP%TCSW%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YTCSW%MP)
    ELSEIF (IC == GFP%TCSLW%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YTCSLW%MP)
    ELSEIF (IC == GFP%SSRD%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YSSRD%MP)
    ELSEIF (IC == GFP%STRD%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YSTRD%MP)
    ELSEIF (IC == GFP%SSRDC%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YSSRDC%MP)
    ELSEIF (IC == GFP%STRDC%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YSTRDC%MP)
    ELSEIF (IC == GFP%BLH%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YBLH%MP)
    ELSEIF (IC == GFP%SUND%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YSUND%MP)
    ELSEIF (IC == GFP%TSRC%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YTSRC%MP)
    ELSEIF (IC == GFP%TTRC%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YTTRC%MP)
    ELSEIF (IC == GFP%SSRC%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YSSRC%MP)
    ELSEIF (IC == GFP%STRC%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YSTRC%MP)
    ELSEIF (IC == GFP%ES%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YES%MP)
    ELSEIF (IC == GFP%SMLT%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YSMLT%MP)
    ELSEIF (IC == GFP%EC10FG%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%Y10FG%MP)
    ELSEIF (IC == GFP%EC10FG3%ICOD) THEN
      I3DIM = SIZE(YSD_VD%Y10FG6)/2
      PFPBUF1(IST:KEND,IFLD)=MAXVAL(PSD_VD(IST:KEND,YSD_VD%Y10FG6(1:I3DIM)%MP),DIM=2)
    ELSEIF (IC == GFP%EC10FG6%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=MAXVAL(PSD_VD(IST:KEND,YSD_VD%Y10FG6(:)%MP),DIM=2)
    ELSEIF (IC == GFP%I10FG%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YI10FG%MP)
    ELSEIF (IC == GFP%LSPF%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YLSPF%MP)
    ELSEIF (IC == GFP%VIMD%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YVIMD%MP)
    ELSEIF (IC == GFP%TCO3%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YTCO3%MP)

!- aerosol optical depths
    ELSEIF (IC == GFP%ODSS%ICOD  ) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YODSS%MP)
    ELSEIF (IC == GFP%ODDU%ICOD  ) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YODDU%MP)
    ELSEIF (IC == GFP%ODOM%ICOD  ) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YODOM%MP)
    ELSEIF (IC == GFP%ODBC%ICOD  ) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YODBC%MP)
    ELSEIF (IC == GFP%ODSU%ICOD  ) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YODSU%MP)
    ELSEIF (IC == GFP%ODSOA%ICOD  ) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YODSOA%MP)
    ELSEIF (IC == GFP%ODNI%ICOD  ) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YODNI%MP)
    ELSEIF (IC == GFP%ODAM%ICOD  ) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YODAM%MP)
    ELSEIF (IC == GFP%ODVFA%ICOD  ) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YODVFA%MP)
    ELSEIF (IC == GFP%ODVSU%ICOD  ) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YODVSU%MP)
    ELSEIF (IC == GFP%ODTOACC%ICOD  ) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YODTOACC%MP)
      !- aerosol (PM) particulate matter 
    ELSEIF (IC == GFP%AEPM1%ICOD ) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YAEPM1%MP)
    ELSEIF (IC == GFP%AEPM25%ICOD ) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YAEPM25%MP)
    ELSEIF (IC == GFP%AEPM10%ICOD ) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YAEPM10%MP)
    ELSEIF (IC == GFP%UVBED%ICOD ) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YUVBED%MP)
    ELSEIF (IC == GFP%UVBEDCS%ICOD ) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YUVBEDCS%MP)

    ELSEIF (ANY(IC == GFP%TCGHG(:)%ICOD)) THEN
       DO JGHG=1,JPGHG
         IF (IC == GFP%TCGHG(JGHG)%ICOD) THEN
           DO IGHG=1,NGHG
             IF(YGHG(IGHG)%IGRBCODE == NGRBGHG(JGHG)) INDEX=IGHG
           ENDDO
           PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YTCGHG(INDEX)%MP)
         ENDIF
       ENDDO
    ELSEIF (ANY(IC == GFP%TCCHEM(:)%ICOD)) THEN
      DO JCHEM=1,NCHEM
        IF (IC == GFP%TCCHEM(JCHEM)%ICOD) THEN
          IF(YCHEM(JCHEM)%IGRIBTC > 0._JPRB )  PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YTCCHEM(JCHEM)%MP)
        ENDIF
      ENDDO
    ELSEIF (ANY(IC == GFP%CHEMFLXO(:)%ICOD)) THEN
      JCHEM1=0_JPIM
      DO JCHEM=1,NCHEM
        IF (IC == GFP%CHEMFLXO(JCHEM)%ICOD) THEN
          DO JCHEM1=1,NCHEM_FLXO
            IF(YCHEM(JCHEM)%IGRBFLXO == YSD_VF%YCHEMFLXO(JCHEM1)%IGRBCODE ) THEN
              PFPBUF1(IST:KEND,IFLD)=PSD_VF(IST:KEND,YSD_VF%YCHEMFLXO(JCHEM1)%MP)
            ENDIF
          ENDDO
        ENDIF
      ENDDO
    ELSEIF (ANY(IC == GFP%CHEMWDFLX(:)%ICOD)) THEN
      JCHEM1=0_JPIM
      DO JCHEM=1,NCHEM
        IF (IC == GFP%CHEMWDFLX(JCHEM)%ICOD) THEN
          DO JCHEM1=1,NCHEM_WDFLX
            IF(YCHEM(JCHEM)%IGRBWDFLX == YSD_VF%YCHEMWDFLX(JCHEM1)%IGRBCODE ) THEN
              PFPBUF1(IST:KEND,IFLD)=PSD_VF(IST:KEND,YSD_VF%YCHEMWDFLX(JCHEM1)%MP)
            ENDIF
          ENDDO
        ENDIF
      ENDDO
    ELSEIF (ANY(IC == GFP%CHEMDDFLX(:)%ICOD)) THEN
      JCHEM1=0_JPIM
      DO JCHEM=1,NCHEM
        IF (IC == GFP%CHEMDDFLX(JCHEM)%ICOD) THEN
          DO JCHEM1=1,NCHEM_DDFLX
            IF(YCHEM(JCHEM)%IGRBDDFLX == YSD_VF%YCHEMDDFLX(JCHEM1)%IGRBCODE ) THEN
              PFPBUF1(IST:KEND,IFLD)=PSD_VF(IST:KEND,YSD_VF%YCHEMDDFLX(JCHEM1)%MP)
            ENDIF
          ENDDO
        ENDIF
      ENDDO
    ELSEIF (ANY(IC == GFP%AERODIAG(:,:)%ICOD)) THEN
      DO JDIAG=1,NPAERODIAG
        DO JAERO=1,NACTAERO
          IF (IC == GFP%AERODIAG(JAERO,JDIAG)%ICOD) THEN
            IF(YAERO_DESC(JAERO)%IGRIBDIAG(JDIAG) > 0._JPRB ) THEN
              PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YAERODIAG(JAERO,JDIAG)%MP)
            ENDIF
          ENDIF
        ENDDO
      ENDDO
    ELSEIF (ANY(IC == GFP%AERO_WVL_DIAG(:,:)%ICOD)) THEN
      DO JDIAG=1,NAERO_WVL_DIAG_TYPES
        DO JAERO=1,NAERO_WVL_DIAG
          IF (IC == GFP%AERO_WVL_DIAG(JAERO,JDIAG)%ICOD) THEN
            IF(YAERO_WVL_DIAG_NL(JAERO)%IGRIBDIAG(JDIAG) > 0._JPRB ) THEN
              PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YAERO_WVL_DIAG(JAERO,JDIAG)%MP)
            ENDIF
          ENDIF
        ENDDO
      ENDDO
    ELSEIF (ANY(IC == GFP%EMIS2D(:)%ICOD)) THEN
      DO JEMIS=1,NEMIS2D
        IF (IC == GFP%EMIS2D(JEMIS)%ICOD) THEN
          PFPBUF1(IST:KEND,IFLD)=PSD_VF(IST:KEND,YSD_VF%YEMIS2D(JEMIS)%MP)
        ENDIF
      ENDDO
    ELSEIF (ANY(IC == GFP%EMIS2DAUX(:)%ICOD)) THEN
      DO JEMIS=1,NEMIS2DAUX
        IF (IC == GFP%EMIS2DAUX(JEMIS)%ICOD) THEN
          PFPBUF1(IST:KEND,IFLD)=PSD_VF(IST:KEND,YSD_VF%YEMIS2DAUX(JEMIS)%MP)
        ENDIF
      ENDDO
    ELSEIF (IC == GFP%CHAR%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_WS(IST:KEND,YSD_WS%YCHAR%MP)
    ELSEIF (IC == GFP%SPAR%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YSPAR%MP)
    ELSEIF (IC == GFP%SUVB%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YSUVB%MP)
    ELSEIF (IC == GFP%CAPE%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YCAPE%MP)
    ELSEIF (IC == GFP%MUCAPE%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YMUCAPE%MP)
    ELSEIF (IC == GFP%MLCAPE50%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YMLCAPE50%MP)
    ELSEIF (IC == GFP%MLCAPE100%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YMLCAPE100%MP)
    ELSEIF (IC == GFP%PDEPL%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YPDEPL%MP)
    ELSEIF (IC == GFP%CAPES%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YCAPES%MP)
    ELSEIF (IC == GFP%SPARC%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YSPARC%MP)
    ELSEIF (IC == GFP%MXCAP6%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=MAXVAL(PSD_VD(IST:KEND,YSD_VD%YMXCAP6(:)%MP),DIM=2)
    ELSEIF (IC == GFP%MXCAPS6%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=MAXVAL(PSD_VD(IST:KEND,YSD_VD%YMXCAPS6(:)%MP),DIM=2)
    ELSEIF (IC == GFP%TROPOTP%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YTROPOTP%MP)
    ELSEIF (IC == GFP%STINC%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YSTINC%MP)
    ELSEIF (IC == GFP%SFDIR%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YSFDIR%MP)
    ELSEIF (IC == GFP%SCDIR%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YSCDIR%MP)
    ELSEIF (IC == GFP%CBASE%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YCBASE%MP)
    ELSEIF (IC == GFP%EC0DEGL%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%Y0DEGL%MP)
    ELSEIF (IC == GFP%ECM10DEGL%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YM10DEGL%MP)
    ELSEIF (IC == GFP%VISIH%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YVISIH%MP)
    ELSEIF (IC == GFP%CIN%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YCIN%MP)
    ELSEIF (IC == GFP%MLCIN50%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YMLCIN50%MP)
    ELSEIF (IC == GFP%MLCIN100%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YMLCIN100%MP)
    ELSEIF (IC == GFP%KINDEX%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YKINDEX%MP)
    ELSEIF (IC == GFP%TTINDEX%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YTTINDEX%MP)
    ELSEIF (IC == GFP%CBASEA%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YCBASEA%MP)
    ELSEIF (IC == GFP%CTOPC%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YCTOPC%MP)
    ELSEIF (IC == GFP%ZTWETB0%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YZTWETB0%MP)
    ELSEIF (IC == GFP%ZTWETB1%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YZTWETB1%MP)
    ELSEIF (IC == GFP%EC100U%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%Y100U%MP)
    ELSEIF (IC == GFP%EC100V%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%Y100V%MP)
    ELSEIF (IC == GFP%EC100SI%ICOD) THEN
      ZU(IST:KEND) = PSD_VD(IST:KEND,YSD_VD%Y100U%MP)
      ZV(IST:KEND) = PSD_VD(IST:KEND,YSD_VD%Y100V%MP)
      PFPBUF1(IST:KEND,IFLD)=SQRT(ZU(IST:KEND)*ZU(IST:KEND) + ZV(IST:KEND)*ZV(IST:KEND))
    ELSEIF (IC == GFP%EC200U%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%Y200U%MP)
    ELSEIF (IC == GFP%EC200V%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%Y200V%MP)
    ELSEIF (IC == GFP%EC200SI%ICOD) THEN
      ZU(IST:KEND) = PSD_VD(IST:KEND,YSD_VD%Y200U%MP)
      ZV(IST:KEND) = PSD_VD(IST:KEND,YSD_VD%Y200V%MP)
      PFPBUF1(IST:KEND,IFLD)=SQRT(ZU(IST:KEND)*ZU(IST:KEND) + ZV(IST:KEND)*ZV(IST:KEND))
    ELSEIF (IC == GFP%ZUST%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YZUST%MP)
    ELSEIF (IC == GFP%EC10NU%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%Y10NU%MP)
    ELSEIF (IC == GFP%EC10NV%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%Y10NV%MP)
    ELSEIF (IC == GFP%DNDZN%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YDNDZN%MP)
    ELSEIF (IC == GFP%DNDZA%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YDNDZA%MP)
    ELSEIF (IC == GFP%DCTB%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YDCTB%MP)
    ELSEIF (IC == GFP%TPLB%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YTPLB%MP)
    ELSEIF (IC == GFP%TPLT%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YTPLT%MP)
    ELSEIF (IC == GFP%CLBT%ICOD) THEN
!     Cloudy brightness temperatures
      ICLBT = ICLBT + 1
      IF (YSD_SM%YCLBT%LSET .AND. ICLBT <= SIZE(PSD_SM,DIM=2)) THEN
        PFPBUF1(IST:KEND,IFLD)=PSD_SM(IST:KEND,ICLBT,YSD_SM%YCLBT%MP)
      ELSE
        CALL ABOR1('HPOS : ERROR IN CLOUDY BRIGHTNESS TEMPERATURES')
      ENDIF
    ELSEIF (IC == GFP%CSBT%ICOD) THEN
!     Clear-sky brightness temperatures
      ICSBT = ICSBT + 1
      IF (YSD_SM%YCSBT%LSET .AND. ICSBT <= SIZE(PSD_SM,DIM=2)) THEN
        PFPBUF1(IST:KEND,IFLD)=PSD_SM(IST:KEND,ICSBT,YSD_SM%YCSBT%MP)
      ELSE
        CALL ABOR1('HPOS : ERROR IN CLEAR SKY BRIGHTNESS TEMPERATURES')
      ENDIF
    ELSEIF (IC == GFP%UCUR%ICOD) THEN
      !     u-component of ocean current.
      IF (LNEMOCOUP) THEN
        PFPBUF1(IST:KEND,IFLD)=PSD_OC(IST:KEND,YSD_OC%YUCURM%MP)
      ELSE
        PFPBUF1(IST:KEND,IFLD)=PSD_VF(IST:KEND,YSD_VF%YUCUR%MP)
      ENDIF
    ELSEIF (IC == GFP%VCUR%ICOD) THEN
      !     v-component of ocean current.
      IF (LNEMOCOUP) THEN
        PFPBUF1(IST:KEND,IFLD)=PSD_OC(IST:KEND,YSD_OC%YVCURM%MP)
      ELSE
        PFPBUF1(IST:KEND,IFLD)=PSD_VF(IST:KEND,YSD_VF%YVCUR%MP)
      ENDIF
    ELSEIF (IC == GFP%LITOTI%ICOD) THEN
!     Total lightning flash density (instantaneous)
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YLITOTI%MP)
    ELSEIF (IC == GFP%LITOTA1%ICOD) THEN
!     Total lightning flash density (averaged over last 1h)
      I1DIM = SIZE(YSD_VD%YLITOTA6)/6
      PFPBUF1(IST:KEND,IFLD)=SUM(PSD_VD(IST:KEND,YSD_VD%YLITOTA6(1:I1DIM)%MP),DIM=2)
    ELSEIF (IC == GFP%LITOTA3%ICOD) THEN
!     Total lightning flash density (averaged over last 3h)
!     Note: This is the second part of the accumulation process, in which hourly accumulations computed 
!           in CPEDIA are summed up to produce the final 3-hourly averages.
      I3DIM = SIZE(YSD_VD%YLITOTA6)/2
      PFPBUF1(IST:KEND,IFLD)=SUM(PSD_VD(IST:KEND,YSD_VD%YLITOTA6(1:I3DIM)%MP),DIM=2)/3._JPRB
    ELSEIF (IC == GFP%LITOTA6%ICOD) THEN
!     Total lightning flash density (averaged over last 6h)
!     Note: This is the second part of the accumulation process, in which hourly accumulations computed 
!           in CPEDIA are summed up to produce the final 6-hourly averages.
      PFPBUF1(IST:KEND,IFLD)=SUM(PSD_VD(IST:KEND,YSD_VD%YLITOTA6(:)%MP),DIM=2)/6._JPRB
    ELSEIF (IC == GFP%LICGI%ICOD) THEN
!     Cloud-to-ground lightning flash density (instantaneous)
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YLICGI%MP)
    ELSEIF (IC == GFP%LICGA1%ICOD) THEN
!     Cloud-to-ground lightning flash density (averaged over last 1h)
      I1DIM = SIZE(YSD_VD%YLICGA6)/6
      PFPBUF1(IST:KEND,IFLD)=SUM(PSD_VD(IST:KEND,YSD_VD%YLICGA6(1:I1DIM)%MP),DIM=2)
    ELSEIF (IC == GFP%LICGA3%ICOD) THEN
!     Cloud-to-ground lightning flash density (averaged over last 3h)
!     Note: This is the second part of the accumulation process, in which hourly accumulations computed 
!           in CPEDIA are summed up to produce the final 3-hourly averages.
      I3DIM = SIZE(YSD_VD%YLICGA6)/2
      PFPBUF1(IST:KEND,IFLD)=SUM(PSD_VD(IST:KEND,YSD_VD%YLICGA6(1:I3DIM)%MP),DIM=2)/3._JPRB
    ELSEIF (IC == GFP%LICGA6%ICOD) THEN
!     Cloud-to-ground lightning flash density (averaged over last 6h)
!     Note: This is the second part of the accumulation process, in which hourly accumulations computed 
!           in CPEDIA are summed up to produce the final 6-hourly averages.
      PFPBUF1(IST:KEND,IFLD)=SUM(PSD_VD(IST:KEND,YSD_VD%YLICGA6(:)%MP),DIM=2)/6._JPRB

    ELSEIF (IC == GFP%PTYPEMODE1%ICOD) THEN
!     Precipitation type most-frequent (mode) (over last 1h)
      I1DIM = SIZE(YSD_VD%YPTYPEOCC6,DIM=1)/6
      DO JP = 1,NPRECTYPES
        ZPTYPESUM(IST:KEND,JP)=SUM(PSD_VD(IST:KEND,YSD_VD%YPTYPEOCC6(1:I1DIM,JP)%MP),DIM=2)
      ENDDO      
      LLPRECIPMASK(IST:KEND,:) = ZPTYPESUM(IST:KEND,:) > 0.1_JPRB
      PFPBUF1(IST:KEND,IFLD) = REAL(MAXLOC(ZPTYPESUM(IST:KEND,:),DIM=2,MASK=LLPRECIPMASK(IST:KEND,:)))
    ELSEIF (IC == GFP%PTYPEMODE3%ICOD) THEN
!     Precipitation type most-frequent (mode) (over last 3h)
!     Note: This is the second part of the accumulation process, in which hourly frequency occurrence computed 
!           in CPEDIA are summed up to produce the final 3-hourly values.
      I3DIM = SIZE(YSD_VD%YPTYPEOCC6,DIM=1)/2
      DO JP = 1,NPRECTYPES
        ZPTYPESUM(IST:KEND,JP)=SUM(PSD_VD(IST:KEND,YSD_VD%YPTYPEOCC6(1:I3DIM,JP)%MP),DIM=2)
      ENDDO      
      LLPRECIPMASK(IST:KEND,:) = ZPTYPESUM(IST:KEND,:) > 0.1_JPRB
      PFPBUF1(IST:KEND,IFLD) = REAL(MAXLOC(ZPTYPESUM(IST:KEND,:),DIM=2,MASK=LLPRECIPMASK(IST:KEND,:)))
    ELSEIF (IC == GFP%PTYPEMODE6%ICOD) THEN
!     Precipitation type most-frequent (mode) (over last 6h)
!     Note: This is the second part of the accumulation process, in which hourly frequency occurrence computed 
!           in CPEDIA are summed up to produce the final 6-hourly values.
      I3DIM = SIZE(YSD_VD%YPTYPEOCC6,DIM=1)
      DO JP = 1,NPRECTYPES
        ZPTYPESUM(IST:KEND,JP)=SUM(PSD_VD(IST:KEND,YSD_VD%YPTYPEOCC6(:,JP)%MP),DIM=2)
      ENDDO      
      LLPRECIPMASK(IST:KEND,:) = ZPTYPESUM(IST:KEND,:) > 0.1_JPRB
      PFPBUF1(IST:KEND,IFLD) = REAL(MAXLOC(ZPTYPESUM(IST:KEND,:),DIM=2,MASK=LLPRECIPMASK(IST:KEND,:)))
      
    ELSEIF (IC == GFP%PTYPESEVR1%ICOD) THEN
!     Precipitation type most-severe (over last 1h)
      I1DIM = SIZE(YSD_VD%YPTYPEOCC6,DIM=1)/6
      DO JP = 1,NPRECTYPES
        ZPTYPESUM(IST:KEND,JP)=SUM(PSD_VD(IST:KEND,YSD_VD%YPTYPEOCC6(1:I1DIM,JP)%MP),DIM=2)
      ENDDO      
      ! Loop over increasing severity and find last (most severe) where occurence > 0
      DO JI=IST,KEND
        DO JP=1,NPRECTYPES
          IF (ZPTYPESUM(JI,NPTYPE_SEV2WMO(JP)) > 0) PFPBUF1(JI,IFLD) = REAL(NPTYPE_SEV2WMO(JP))
        ENDDO
      ENDDO
    ELSEIF (IC == GFP%PTYPESEVR3%ICOD) THEN
!     Precipitation type most-severe (over last 3h)
!     Note: This is the second part of the accumulation process, in which hourly severity occurrence computed 
!           in CPEDIA are combined to produce the final 3-hourly values.
      I3DIM = SIZE(YSD_VD%YPTYPEOCC6,DIM=1)/2
      DO JP = 1,NPRECTYPES
        ZPTYPESUM(IST:KEND,JP)=SUM(PSD_VD(IST:KEND,YSD_VD%YPTYPEOCC6(1:I3DIM,JP)%MP),DIM=2)
      ENDDO      
      ! Loop over increasing severity and find last (most severe) where occurence > 0
      DO JI=IST,KEND
        DO JP=1,NPRECTYPES
          IF (ZPTYPESUM(JI,NPTYPE_SEV2WMO(JP)) > 0._JPRB) PFPBUF1(JI,IFLD) = REAL(NPTYPE_SEV2WMO(JP))
        ENDDO
      ENDDO
    ELSEIF (IC == GFP%PTYPESEVR6%ICOD) THEN
!     Precipitation type most-severe (over last 6h)
!     Note: This is the second part of the accumulation process, in which hourly severity occurrence computed 
!           in CPEDIA are combined to produce the final 6-hourly values.
      I3DIM = SIZE(YSD_VD%YPTYPEOCC6,DIM=1)
      DO JP = 1,NPRECTYPES
        ZPTYPESUM(IST:KEND,JP)=SUM(PSD_VD(IST:KEND,YSD_VD%YPTYPEOCC6(:,JP)%MP),DIM=2)
      ENDDO      
      ! Loop over increasing severity and find last (most severe) where occurence > 0
      DO JI=IST,KEND
        DO JP=1,NPRECTYPES
          IF (ZPTYPESUM(JI,NPTYPE_SEV2WMO(JP)) > 0) PFPBUF1(JI,IFLD) = REAL(NPTYPE_SEV2WMO(JP))
        ENDDO
      ENDDO

    ELSEIF (IC == GFP%SDSRP%ICOD) THEN
      PFPBUF1(IST:KEND,IFLD)=PSD_VD(IST:KEND,YSD_VD%YSDSRP%MP)
    ELSEIF (IC == GFP%ICTH%ICOD) THEN
      IF (LNEMOCOUP.AND.LNEMOGRIBFLDS) THEN
        PFPBUF1(IST:KEND,IFLD)=PSD_OC(IST:KEND,YSD_OC%YICTH%MP)
      ELSE
        PFPBUF1(IST:KEND,IFLD)=0.0_JPRB
      ENDIF
    ELSEIF (IC == GFP%SSH%ICOD) THEN
      IF (LNEMOCOUP.AND.LNEMOGRIBFLDS) THEN
        PFPBUF1(IST:KEND,IFLD)=PSD_OC(IST:KEND,YSD_OC%YSSH%MP)
      ELSE
        PFPBUF1(IST:KEND,IFLD)=0.0_JPRB
      ENDIF
    ELSEIF (IC == GFP%EC20D%ICOD) THEN
      IF (LNEMOCOUP.AND.LNEMOGRIBFLDS) THEN
        PFPBUF1(IST:KEND,IFLD)=PSD_OC(IST:KEND,YSD_OC%Y20D%MP)
      ELSE
        PFPBUF1(IST:KEND,IFLD)=0.0_JPRB
      ENDIF
    ELSEIF (IC == GFP%MLD%ICOD) THEN
      IF (LNEMOCOUP.AND.LNEMOGRIBFLDS) THEN
        PFPBUF1(IST:KEND,IFLD)=PSD_OC(IST:KEND,YSD_OC%YMLD%MP)
      ELSE
        PFPBUF1(IST:KEND,IFLD)=0.0_JPRB
      ENDIF
    ELSEIF (IC == GFP%SSS%ICOD) THEN
      IF (LNEMOCOUP.AND.LNEMOGRIBFLDS) THEN
        PFPBUF1(IST:KEND,IFLD)=PSD_OC(IST:KEND,YSD_OC%YSSS%MP)
      ELSE
        PFPBUF1(IST:KEND,IFLD)=0.0_JPRB
      ENDIF
    ELSEIF (IC == GFP%TEM3%ICOD) THEN
      IF (LNEMOCOUP.AND.LNEMOGRIBFLDS) THEN
        PFPBUF1(IST:KEND,IFLD)=PSD_OC(IST:KEND,YSD_OC%YTEM3%MP)
      ELSE
        PFPBUF1(IST:KEND,IFLD)=0.0_JPRB
      ENDIF
    ELSEIF (IC == GFP%SAL3%ICOD) THEN
      IF (LNEMOCOUP.AND.LNEMOGRIBFLDS) THEN
        PFPBUF1(IST:KEND,IFLD)=PSD_OC(IST:KEND,YSD_OC%YSAL3%MP)
      ELSE
        PFPBUF1(IST:KEND,IFLD)=0.0_JPRB
      ENDIF
    ELSE

      DO JVAR=1,JPOSFSU
        IF (IC == GFP%FSU(JVAR)%ICOD) THEN
          !     Compute here !
          PFPBUF1(IST:KEND,IFLD)=0.0_JPRB
          EXIT
        ENDIF
      ENDDO
      IF (JPOSVX2 >= YDMODEL%YRML_PHY_G%YRDPHY%NVXTR2) THEN
        DO JVAR=1,JPOSVX2
          IF (IC == GFP%VX2(JVAR)%ICOD) THEN
            IEXTR2=IEXTR2+1
            !     Additional extra physics fields
            PFPBUF1(IST:KEND,IFLD)=PSD_X2(IST:KEND,YSD_X2%YX2(IEXTR2)%MP)
          ENDIF
        ENDDO
      ELSE
        CALL ABOR1('HPOS : JPOSVX2 SMALLER THAN NVXTR2 !')
      ENDIF

    ENDIF

  ENDDO

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('HPOS',1,ZHOOK_HANDLE)
END SUBROUTINE HPOS
