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

SUBROUTINE UPDTIM(YDGEOMETRY,YDSURF,YDMODEL,KSTEP,PTDT,PTSTEP,LDCLUPD,LDUPDECAEC,PI05)

!**** *UPDTIM*  - UPDATE TIME OF THE MODEL

!     Purpose.
!     --------
!     UPDATE TIME OF THE MODEL

!**   Interface.
!     ----------
!        *CALL* *UPDTIM(...)

!        Explicit arguments :
!        --------------------
!        KSTEP   : TIME STEP INDEX
!        PTDT    : TIME STEP LEAPFROG
!        PTSTEP  : TIME STEP
!        LDCLUPD : .TRUE. IF CALLED DURING A MODEL RUN DURING
!                    INITIALISATION OF TENDENCIES
!                    (ONLY RELEVANT FOR A CLIMATE FILE UPDATE)

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
!        Mats Hamrud and Philippe Courtier  *ECMWF*
!        Original : 87-10-15

!     Modifications.
!     --------------
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        M.Hamrud      20-Apr-2004 Save time by supressing repeated calls to SUECAEOR etc.
!        Modified by A. Alias: 05-08-30 GMGEC/EAC modif.
!                    M. Deque   -  nudging with model data
!                    JPh Piedelievre - Coupled modes (UPDCPL)
!                    M. Deque   - nudging with time-variable coefficients
!                    P. Marquet - LWNUDG added
!                    M. Deque   - Call to CORMASS2
!                    M. Deque   - grid point nudging
!        M.Janiskova 060111 bug-fix for use of new aerosols in TL/AD
!        S.Serrar  : 06-01-01 time updating of CO2 surface fields
!        M. Ko"hler  6-6-2006 Single Column Model option (LSCMEC)
!        M.Hamrud      01-Jul-2006 Revised surface fields
!        Modified by A. Alias: 06-03-10 clean up
!        JJMorcrette 20080318 O.Boucher SO4 climatol (Observ and A1B)
!        JJMorcrette 20080423 3D climatologies of CO2, CH4, N2O for radiation
!        Modified by A. Alias: 07-10-19 NTOTFNUDG  moved to YOMDIM
!        Y. Takaya 01-Feb-2009 call UPDCLIE for the ocean mixed layer model
!        Modified by A. Alias : 07-08-09 Frequency of nudging is now control by NFRNUDG
!        Modified by A. Voldoire : Dec 2010 Add RCODECN,RSIDECN,RCOVSRN,RSIVSRN
!                                  (next time step value for LMSE)
!        Modified by A. Alias : 03-2011 LASTRF to prevent any drift in insolation (A.Voldoire)
!                               Add RCODECN,RSIDECN,RCOVSRN,RSIVSRN
!                                  (next time step value for LMSE not for LMPA) (A.Voldoire)
!        D.Lindstedt 20101124: New updcli_mse for SURFEX and climate runs. Frequency
!                              of updating is the interval of LBC
!        JJMorcrette 20110726 volcanic eruptions in MACC configuration
!        S. YANG, H. Hersbach: 09-08-2010: CALL SUECOZV for varying ozone in AMIP runs
!        JJMorcrette 20111007  Possibility of different frequency of full rad. for EPS
!        Linus Magnusson 12-07-13: Modifications to the relaxation part
!        P Bechtold add invariant RDAYI for SW radiation
!        N.Semane+P.Bechtold 04-10-2012 replace 86400/RDAYI by RDAY and 3600 by RHOUR consistently with phys_ec/updtier.F90
!        T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!        K. Yessad (July 2014): Move some variables.
!        JJMorcrette 20130805  MACC-derived aerosol climatology
!        J. Flemming update of emission and dry dep fluxes for composition : call to updclie_compo
!        R Hogan (Oct 2014) Update solar time for radiation scheme only at radiation timesteps
!        R Hogan (Apr 2015) Solar time for radiation scheme controlled by LMannersSwUpdate
!        R Hogan (Apr 2015) Save solar time for radiation scheme
!        C. Roberts/R.Senan 26/01/2017 Support for CMIP6 forcings
!        R Hogan (Nov 2017) Call UPDRGAS every timestep for correct solar irradiance
!        A Bozzo (Jan 2018) Support for 3D aerosol climatology in input
!        R Hogan (Jan 2019) Removed stuff for old cycle 15 radiation scheme
!     ------------------------------------------------------------------

USE TYPE_MODEL         , ONLY : MODEL
USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE PARKIND1           , ONLY : JPIM, JPRB, JPRD
USE YOMHOOK            , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMMP0             , ONLY : LSCMEC
USE YOMCT0             , ONLY : NFRCO, L_OOPS
USE YOMLUN             , ONLY : NULOUT
USE YOMRIP0            , ONLY : NINDAT, NSSSSS, RTIMST, LASTRF
USE YOMCST             , ONLY : RPI, RDAY, RHOUR, REA, REPSM, RA, RI0, RV, RCPV, RETV,&
 &                              RCW, RCS, RLVTT, RLSTT, RTT, RALPW, RBETW, RGAMW, RALPS, RBETS, RGAMS,&
 &                              RALPD, RBETD, RGAMD
USE YOMNUD             , ONLY : NFNUDG, NFRNUDG, LNUDG, LWNUDG, NTOTFNUDG2, NTOTFNUDG3
USE YOMSNU             , ONLY : XPNUDG, XWNUDG
USE YOMRLX             , ONLY : NFRLXG, NFRLXU, LRLXG
USE YOMSRLX            , ONLY : XPRLXG
USE YOMDYNCORE         , ONLY : LAPE

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)      ,INTENT(IN)    :: YDGEOMETRY
TYPE(TSURF)         ,INTENT(INOUT) :: YDSURF
TYPE(MODEL)         ,INTENT(INOUT) :: YDMODEL
INTEGER(KIND=JPIM)  ,INTENT(IN)    :: KSTEP
REAL(KIND=JPRB)     ,INTENT(IN)    :: PTDT
REAL(KIND=JPRB)     ,INTENT(IN)    :: PTSTEP
LOGICAL             ,INTENT(IN)    :: LDCLUPD
LOGICAL             ,INTENT(IN), OPTIONAL :: LDUPDECAEC
REAL(KIND=JPRB)     ,INTENT(IN), OPTIONAL :: PI05

!     ------------------------------------------------------------------
INTEGER(KIND=JPIM) :: IBASE, IFRHIS, IGP, IGP_MSE, IMINUT, IPR, ISTADD,&
 & ISTASS, ISTASS0, ISTP1, ISTPF, ITIME, IZT, JSTEP, IZTN, ISS, IFRLXG
INTEGER(KIND=JPIM) :: ISEC, ICNT, ISTEP, IUPGHG

REAL(KIND=JPRB) :: ZANGOZC, ZDEASOM, ZDECLIM, ZDEL,&
 & ZEQTIMM, ZHGMT, ZI0, ZSOVRF, ZSTATI,&
 & ZTETA, ZTHETOZ, ZTIMTR, ZWSOVRF, ZSEASON,&
 & ZSTATIN, ZHGMTN, ZTIMTRN, ZTETAN, ZEQTIMN, ZSOVRN,&
 & ZWSOVRN, ZDECLIN, ZTI, ZTIN, ZTETAR
! The following variables are not used so have been commented out
! ZCOTHOZ, ZSITHOZ

INTEGER(KIND=JPIM) :: JUV
INTEGER(KIND=JPIM) :: IFRNUDG
INTEGER(KIND=JPIM) :: IJ0,IM0,IA0,IMOIS,IJOUR,IAN,ILMOIS(12)

INTEGER(KIND=JPIM) :: ISA
REAL(KIND=JPRB) :: ZHGMT0, ZTIA, ZTETAA, ZDECLIA, ZEQTIMA, ZSOVRA, ZWSOVRA
INTEGER(KIND=JPIM) :: IRLXI   ! number of intervals between two reference fields

REAL(KIND=JPRB) :: ZUNIT
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
LOGICAL :: LLFIRSTCALL = .TRUE.
LOGICAL :: LLUPDECAEC

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "updecaec.intfb.h"

#include "suecso4.intfb.h"
#include "updecozc.intfb.h"
#include "updecozcaqua.intfb.h"
#include "suecozo.intfb.h"
#include "updecozv.intfb.h"
#include "upd_ghgclim.intfb.h"
#include "updcli.intfb.h"
#include "updcli_mse.intfb.h"
#include "updclie.intfb.h"
#include "updclie_oasis.intfb.h"
#include "updcpl.intfb.h"
#include "updmoon.intfb.h"
#include "updnud.intfb.h"
#include "updrlxref.intfb.h"
#include "updo3ch.intfb.h"
#include "updrgas.intfb.h"
#include "su_aervole.intfb.h"
#include "updclie_compo.intfb.h"

#include "fctast.func.h"
#include "fcttim.func.h"
#include "fcttrm.func.h"
#include "updcal.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('UPDTIM',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, &
 &  YDPHY=>YDMODEL%YRML_PHY_MF%YRPHY,YDEAERVOL=>YDMODEL%YRML_PHY_AER%YREAERVOL,YDMCC=>YDMODEL%YRML_AOC%YRMCC, &
 & YDRIP=>YDMODEL%YRML_GCONF%YRRIP,YDDPHY=>YDMODEL%YRML_PHY_G%YRDPHY,YDERAD=>YDMODEL%YRML_PHY_RAD%YRERAD, &
 & YDARPHY=>YDMODEL%YRML_PHY_MF%YRARPHY,YDEUVRAD=>YDMODEL%YRML_PHY_RAD%YREUVRAD,YDERIP=>YDMODEL%YRML_PHY_RAD%YRERIP, &
 & YDEPHLI=>YDMODEL%YRML_PHY_SLIN%YREPHLI,YDPHY3=>YDMODEL%YRML_PHY_MF%YRPHY3,YDPHY2=>YDMODEL%YRML_PHY_MF%YRPHY2, &
  & YDEAERSRC=>YDMODEL%YRML_PHY_AER%YREAERSRC,YGFL=>YDMODEL%YRML_GCONF%YGFL,YDEPHY=>YDMODEL%YRML_PHY_EC%YREPHY, &
  & YDERDI=>YDMODEL%YRML_PHY_RAD%YRERDI, YDSIMPHL=>YDMODEL%YRML_PHY_MF%YRSIMPHL,YDDYNA=>YDMODEL%YRML_DYN%YRDYNA)
ASSOCIATE(NAERO=>YGFL%NAERO, &
 & NSMAX=>YDDIM%NSMAX, &
 & NFLEVG=>YDDIMV%NFLEVG, &
 & NTSSG=>YDDPHY%NTSSG, &
 & RCODECA=>YDEAERSRC%RCODECA, RCOVSRA=>YDEAERSRC%RCOVSRA, &
 & RSIDECA=>YDEAERSRC%RSIDECA, RSIVSRA=>YDEAERSRC%RSIVSRA, &
 & NVOLERUP=>YDEAERVOL%NVOLERUP, &
 & LPHYLIN=>YDEPHLI%LPHYLIN, &
 & LEO3CH=>YDEPHY%LEO3CH, LEOCML=>YDEPHY%LEOCML, LOCMLTKE=>YDEPHY%LOCMLTKE, &
 & LEPHYS=>YDEPHY%LEPHYS, LNEEONLINE=>YDEPHY%LNEEONLINE, &
 & LMANNERSSWUPDATE=>YDERAD%LMANNERSSWUPDATE, &
 & LCENTREDTIMESZA=>YDERAD%LCENTREDTIMESZA, LECO2VAR=>YDERAD%LECO2VAR, &
 & LERAD1H=>YDERAD%LERAD1H, LESO4HIS=>YDERAD%LESO4HIS, LHGHG=>YDERAD%LHGHG, &
 & LHVOLCA=>YDERAD%LHVOLCA, LNEWAER=>YDERAD%LNEWAER, LPERPET=>YDERAD%LPERPET, &
 & NAERMACC=>YDERAD%NAERMACC, NGHGRAD=>YDERAD%NGHGRAD, &
 & NHINCSOL=>YDERAD%NHINCSOL, NLNGR1H=>YDERAD%NLNGR1H, NOZOCL=>YDERAD%NOZOCL, &
 & NRADE1H=>YDERAD%NRADE1H, NRADE3H=>YDERAD%NRADE3H, NRADELG=>YDERAD%NRADELG, &
 & NRADFR=>YDERAD%NRADFR, NRADNFR=>YDERAD%NRADNFR, NRADSFR=>YDERAD%NRADSFR, &
 & NUV=>YDERAD%NUV, LAER3D=>YDERAD%LAER3D, &
 & RCARDI=>YDERDI%RCARDI, RSOLINC=>YDERDI%RSOLINC, &
 & LUVPROC=>YDEUVRAD%LUVPROC, RSUVB=>YDEUVRAD%RSUVB, RSUVB0=>YDEUVRAD%RSUVB0, &
 & NGPTOT=>YDGEM%NGPTOT, &
 & LMCC01=>YDMCC%LMCC01, LMCC01_MSE=>YDMCC%LMCC01_MSE, &
 & LMCC03=>YDMCC%LMCC03, LMCC04=>YDMCC%LMCC04, LMCC05=>YDMCC%LMCC05, &
 & LMCCEC=>YDMCC%LMCCEC, LMCC_COMPO=>YDMCC%LMCC_COMPO, &
 & LNEMO1WAY=>YDMCC%LNEMO1WAY, NFRCPL=>YDMCC%NFRCPL, &
 & RCODECM=>YDERIP%RCODECM, RSIDECM=>YDERIP%RSIDECM, RCOVSRM=>YDERIP%RCOVSRM, &
 & RSIVSRM=>YDERIP%RSIVSRM, &
 & RSOVRM=>YDERIP%RSOVRM, &
 & RWSOVRM=>YDERIP%RWSOVRM, &
 & NSTADD=>YDRIP%NSTADD, NSTASS=>YDRIP%NSTASS, RCODEC=>YDRIP%RCODEC, &
 & RCODECF=>YDRIP%RCODECF, RCODECLU=>YDRIP%RCODECLU, RCODECN=>YDRIP%RCODECN, &
 & RCOVSR=>YDRIP%RCOVSR, RCOVSRF=>YDRIP%RCOVSRF, RCOVSRLU=>YDRIP%RCOVSRLU, &
 & RCOVSRN=>YDRIP%RCOVSRN, RDEASO=>YDRIP%RDEASO, RDECLI=>YDRIP%RDECLI, &
 & RDECLU=>YDRIP%RDECLU, RDTS22=>YDRIP%RDTS22, RDTS62=>YDRIP%RDTS62, &
 & RDTSA=>YDRIP%RDTSA, RDTSA2=>YDRIP%RDTSA2, REQTIM=>YDRIP%REQTIM, &
 & RHGMT=>YDRIP%RHGMT, RIP0=>YDRIP%RIP0, RSIDEC=>YDRIP%RSIDEC, &
 & RSIDECF=>YDRIP%RSIDECF, RSIDECLU=>YDRIP%RSIDECLU, RSIDECN=>YDRIP%RSIDECN, &
 & RSIVSR=>YDRIP%RSIVSR, RSIVSRF=>YDRIP%RSIVSRF, RSIVSRLU=>YDRIP%RSIVSRLU, &
 & RSIVSRN=>YDRIP%RSIVSRN, RSOVR=>YDRIP%RSOVR, RSTATI=>YDRIP%RSTATI, &
 & RTDT=>YDRIP%RTDT, RTIMTR=>YDRIP%RTIMTR, RTMOLT=>YDRIP%RTMOLT, &
 & RWSOVR=>YDRIP%RWSOVR, TSTEP=>YDRIP%TSTEP, &
 & YSD_VAD=>YDSURF%YSD_VAD, YSD_VFD=>YDSURF%YSD_VFD, YSD_VPD=>YDSURF%YSD_VPD, &
 & YSD_VVD=>YDSURF%YSD_VVD, LMPA=>YDARPHY%LMPA, LMSE=>YDARPHY%LMSE, &
 & LMPHYS=>YDPHY%LMPHYS, LSIMPH=>YDSIMPHL%LSIMPH, LOZONE=>YDPHY%LOZONE, &
 & LRAYFM=>YDPHY%LRAYFM, LRAYFM15=>YDPHY%LRAYFM15, LRAYLU=>YDPHY%LRAYLU, &
 & LRMU0M=>YDPHY%LRMU0M, LSOLV=>YDPHY%LSOLV, TSPHY=>YDPHY2%TSPHY, &
 & RII0=>YDPHY3%RII0)
!     ------------------------------------------------------------------

LLUPDECAEC=.TRUE.
IF (PRESENT(LDUPDECAEC)) LLUPDECAEC=LDUPDECAEC

! Time-step length (seconds)
ITIME=NINT(PTSTEP)

IF (YDDYNA%LTWOTL) THEN
  ! In the two-level timestepping scheme, the solar zenith angle
  ! should be computed at a time half-way between the current and
  ! future timestep, hence the 0.5 here.  IZT is the number of seconds
  ! since the start of the forecast.
  IZT=NINT(PTSTEP*(REAL(KSTEP,JPRB)+0.5_JPRB))
ELSE
  IZT=ITIME*KSTEP
ENDIF

!--
IF (LPERPET) THEN
  ISEC=IZT/NINT(RDAY)
  IZT=IZT-ISEC*NINT(RDAY)
ENDIF
!--

RSTATI=REAL(IZT,JPRB)
NSTADD=IZT/NINT(RDAY)
NSTASS=MOD(IZT,NINT(RDAY))
ISTASS=MOD(IZT+NINT(RDAY)/2-NSSSSS,NINT(RDAY))
RTIMTR=RTIMST+RSTATI
IPR=0
IF(IPR == 1)THEN
  WRITE(UNIT=NULOUT,FMT='(1X,'' TIME OF THE MODEL '',E20.14,&
   & '' TIME SINCE START '',E20.14)') RTIMTR,RSTATI
ENDIF
RHGMT=REAL(MOD(NINT(RSTATI)+NSSSSS,NINT(RDAY)),JPRB)

IF (NAERO /= 0) THEN
  IJ0=NDD(NINDAT)
  IM0=NMM(NINDAT)
  IA0=NCCAA(NINDAT)
  CALL UPDCAL (IJ0,IM0,IA0,NSTADD,IJOUR,IMOIS,IAN,ILMOIS,NULOUT)
!-- for computing the mean insolation for the day, set up the clock
!   at 12 GMT for 0E-W
  ISA=0
!  ZHGMT0=0._JPRB
  ZHGMT0=43200._JPRB
  ZTIA=RTIME(IAN,IMOIS,IJOUR,ISA,RDAY)
  ZTETAA=RTETA(ZTIA)
  ZDECLIA=RDS(ZTETAA)
  ZEQTIMA=RET(ZTETAA)
  ZSOVRA =ZEQTIMA+ZHGMT0
  ZWSOVRA=ZSOVRA*2.0_JPRB*RPI/RDAY

  RCODECA=COS(ZDECLIA)
  RSIDECA=SIN(ZDECLIA)

  RCOVSRA=COS(ZWSOVRA)
  RSIVSRA=SIN(ZWSOVRA)
ENDIF

IF (LASTRF) THEN
  IJ0=NDD(NINDAT)
  IM0=NMM(NINDAT)
  IA0=NCCAA(NINDAT)
  CALL UPDCAL (IJ0,IM0,IA0,NSTADD,IJOUR,IMOIS,IAN,ILMOIS,NULOUT)
  ! Keep a year close to year 2000 (reference year)
  ! so that RET is correctly used
  IAN = MOD(IAN,4)+2000
  ISS = INT(RHGMT)
  ZTI=RTIME(IAN,IMOIS,IJOUR,ISS,RDAY)
  ZTETA=RTETA(ZTI)
ELSE
  ZTETA=RTETA(RTIMTR)
ENDIF
ZTETAR=ZTETA
IF( LAPE ) THEN
  RDEASO=RRSAQUA(ZTETA)
  RDECLI=RDSAQUA(ZTETA)
  REQTIM=RETAQUA(ZTETA)
ELSE
  RDEASO=RRS(ZTETA)
  RDECLI=RDS(ZTETA)
  REQTIM=RET(ZTETA)
ENDIF
RSOVR =REQTIM+RHGMT
RWSOVR=RSOVR*2.0_JPRB*RPI/RDAY
IF (PRESENT(PI05)) THEN
    ZI0=PI05
ELSE
  IF (NHINCSOL /= 0) THEN
    ZI0=RSOLINC
  ELSE
    ZI0=RI0
  ENDIF
ENDIF
ZSEASON= REA*REA/(RDEASO*RDEASO)
IF( LAPE ) THEN
  ZSEASON=1._JPRB
ENDIF
RIP0=ZI0*ZSEASON
RII0 = RIP0
IF (LUVPROC) THEN
  DO JUV=1,NUV
    RSUVB(JUV)=RSUVB0(JUV)*ZSEASON
  ENDDO
ENDIF

RCODEC=COS(RDECLI)
RSIDEC=SIN(RDECLI)

RCOVSR=COS(RWSOVR)
RSIVSR=SIN(RWSOVR)

IF (NAERO /= 0) THEN
  IF (LASTRF) THEN
    WRITE(NULOUT,FMT='("Std: ",8E12.5)') RCODEC,RSIDEC,RCOVSR,RSIVSR,ZTI,ZTETA,RDECLI,REQTIM
  ELSE
    ! ZTI will be uninitialised in this case
    WRITE(NULOUT,FMT='("Std: ",4E12.5,A12,3E12.5)') RCODEC,RSIDEC,RCOVSR,RSIVSR,'-----',ZTETA,RDECLI,REQTIM
  ENDIF
  WRITE(NULOUT,FMT='("Aer: ",8E12.5)') RCODECA,RSIDECA,RCOVSRA,RSIVSRA, ZTIA,ZTETAA,ZDECLIA,ZEQTIMA
ENDIF


IF(LMSE)THEN
! Calculation of RCODEC, RSIDEC, RCOVSR, RSIVSR at next time step for SURFEX
  IF (YDDYNA%LTWOTL) THEN
    IZTN=NINT(PTSTEP*(REAL(KSTEP+1,JPRB)+0.5_JPRB))
  ELSE
    IZTN=ITIME*(KSTEP+1)
  ENDIF
  ZSTATIN=REAL(IZTN,JPRB)
  ZHGMTN=REAL(MOD(NINT(ZSTATIN)+NSSSSS,NINT(RDAY)),JPRB)
  ZTIMTRN=RTIMST+KSTEP+ZSTATIN
  ZTETAN=RTETA(ZTIMTRN)
  ZEQTIMN=RET(ZTETAN)
  ZSOVRN=ZEQTIMN+ZHGMTN
  ZWSOVRN=ZSOVRN*2.0_JPRB*RPI/RDAY

  ZDECLIN=RDS(ZTETAN)

  RCODECN=COS(ZDECLIN)
  RSIDECN=SIN(ZDECLIN)

  RCOVSRN=COS(ZWSOVRN)
  RSIVSRN=SIN(ZWSOVRN)
ENDIF

ICNT=1
ISTEP=ICNT*INT(KSTEP/ICNT)
IF (ISTEP == KSTEP) THEN
  WRITE(NULOUT,FMT='(''ISTEP='',I8,'' RSOLINC='',F10.4,'' RI0='',F10.4,&
    & '' ZI0='',F10.4,'' RIP0='',F10.4)') ISTEP,RSOLINC,RI0,ZI0,RIP0
  WRITE(NULOUT,FMT='(''ISTEP='',I8,'' ZI0='',F10.4,'' ZSEASON='',F12.9,&
    & '' REA='',F15.2,'' RDEASO='',F15.2)') ISTEP,ZI0,ZSEASON,REA,RDEASO
ENDIF
! Not for AROME
IF(LMSE.AND..NOT.LMPA)THEN
! Calculation of RCODEC, RSIDEC, RCOVSR, RSIVSR at next time step for SURFEX
  IF (YDDYNA%LTWOTL) THEN
    IZTN=NINT(PTSTEP*(REAL(KSTEP+1,JPRB)+0.5_JPRB))
  ELSE
    IZTN=ITIME*(KSTEP+1)
  ENDIF
  ZSTATIN=REAL(IZTN,JPRB)
  ZHGMTN=REAL(MOD(NINT(ZSTATIN)+NSSSSS,NINT(RDAY)),JPRB)
  ZTIMTRN=RTIMST+KSTEP+ZSTATIN
  IF (LASTRF) THEN
    ISTADD = NSTADD + PTSTEP/NINT(RDAY)
    CALL UPDCAL (IJ0,IM0,IA0,ISTADD,IJOUR,IMOIS,IAN,ILMOIS,NULOUT)
    IAN = MOD(IAN,4)+2000
    ISS = INT(ZHGMTN)
    ZTIN=RTIME(IAN,IMOIS,IJOUR,ISS,RDAY)
    ZTETAN=RTETA(ZTIN)
  ELSE
    ZTETAN=RTETA(ZTIMTRN)
  ENDIF
  ZEQTIMN=RET(ZTETAN)
  ZSOVRN=ZEQTIMN+ZHGMTN
  ZWSOVRN=ZSOVRN*2.0_JPRB*RPI/RDAY

  ZDECLIN=RDS(ZTETAN)

  RCODECN=COS(ZDECLIN)
  RSIDECN=SIN(ZDECLIN)

  RCOVSRN=COS(ZWSOVRN)
  RSIVSRN=SIN(ZWSOVRN)
ENDIF

! MOON LOCATION COMPUTATIONS.

IF(LRAYLU) THEN

  CALL UPDMOON(YDRIP)

  RCODECLU=COS(RDECLU)
  RSIDECLU=SIN(RDECLU)

  RCOVSRLU=COS(RTMOLT)
  RSIVSRLU=SIN(RTMOLT)

ENDIF

RDTSA=0.5_JPRB*PTDT/RA
RDTSA2=RDTSA**2
RDTS62=RDTSA2/6._JPRB
RDTS22=RDTSA2/2.0_JPRB

RTDT=PTDT
TSPHY = MAX(PTDT,1.0_JPRB)

IF(LNUDG)THEN
  ! CALCULATION OF THE TIME WEIGHTS FOR NUDGING
  ! one assumes that the frequency of input restart
  ! is the same as the output restart, but at least once a day
  IFRHIS=MIN(NFRNUDG,NINT(RDAY)/NINT(TSTEP))
  ZUNIT=RHOUR
  IF ((TSTEP > 0.0_JPRB).AND.(NFRNUDG < 0)) THEN
    NFRNUDG=NINT((REAL(-NFRNUDG,JPRB)*ZUNIT)/TSTEP)
  ENDIF

  ! If LNUDG=.TRUE. make sure NFNUDG has the right value and NFRNUDG
  IF (.NOT.(NFNUDG == 2.OR.NFNUDG == 7)) THEN
    WRITE(NULOUT,'('' ERROR NFNUDG != '',I2)')NFNUDG
    CALL ABOR1('SUNUD: ABOR1 CALLED')
  ENDIF
  IF (NFNUDG == 7) THEN
    IFRNUDG=NINT(6*ZUNIT/TSTEP)
    IF (NFRNUDG /= IFRNUDG) THEN
      WRITE(NULOUT,'('' ERROR with NFNUDG=7 you must have NFRNUDG=6hours'')')
      CALL ABOR1('SUNUD: ABOR1 CALLED')
    ENDIF
  ELSE
    !  NFNUDG=2
    IFRNUDG=NINT(24*ZUNIT/TSTEP)
    IFRNUDG=NINT(24*ZUNIT/TSTEP)
    IF ((NFRNUDG == 0).OR.(NFRNUDG > IFRNUDG)) THEN
      WRITE(NULOUT,'('' ERROR : NFNUDG=2  Nudging at least once a day NFRNUDG= '',I2)')NFRNUDG
      CALL ABOR1('SUNUD: ABOR1 CALLED')
    ENDIF
  ENDIF

  DO JSTEP=1,NFNUDG
    XPNUDG(JSTEP)=0.0_JPRB
  ENDDO
! ISTPF: LENGTH BETWEEN TWO OBSERVATIONS (S)
  ISTPF=IFRHIS*ITIME
! ISTP1: TIME OF THE STEP (S)
  ISTP1=MOD(KSTEP*ITIME,NINT(RDAY))
  ZDEL=REAL(ISTP1,JPRB)/REAL(ISTPF,JPRB)-REAL(ISTP1/ISTPF,JPRB)
  IF(NFNUDG == 7) THEN
    IBASE=2+ISTP1/ISTPF
    XPNUDG(IBASE-1)=-0.5_JPRB*ZDEL*(1.0_JPRB-ZDEL)**2
    XPNUDG(IBASE)=1.0_JPRB-0.5_JPRB*ZDEL**2*(5._JPRB-3._JPRB*ZDEL)
    XPNUDG(IBASE+1)=1.0_JPRB-0.5_JPRB*(1.0_JPRB-ZDEL)**2*(2.0_JPRB+3._JPRB*ZDEL)
    XPNUDG(IBASE+2)=-0.5_JPRB*ZDEL**2*(1.0_JPRB-ZDEL)
  ELSEIF(NFNUDG == 2) THEN
    IBASE=1
    XPNUDG(IBASE)=1.0_JPRB-ZDEL
    XPNUDG(IBASE+1)=ZDEL
  ENDIF
  IF(LWNUDG)THEN
    XWNUDG=MAX(0.0_JPRB,COS(REAL(2*KSTEP,JPRB)/REAL(NFRNUDG,JPRB)*RPI))
  ELSE
    XWNUDG=1.0_JPRB
  ENDIF
ENDIF

!  CALCULATION OF THE TIME WEIGHTS FOR RELAXATION
!  one assumes that the frequency of input restart
!  is the same as the output restart, but at least once a day
IF(LRLXG)THEN
  IF (.NOT.ALLOCATED (XPRLXG)) ALLOCATE(XPRLXG(NFRLXG))

! IRLXI: Number of intervals between two reference fields
  IRLXI=NINT(REAL(NFRLXU)/(TSTEP/3600.0))

  IFRHIS=MIN(IRLXI,86400/NINT(TSTEP))
  IF(NFRLXG == 2 ) THEN
    IFRLXG=43200/(IFRHIS*NINT(TSTEP))
  ENDIF
  IF(NFRLXG /= IFRLXG)THEN
    WRITE(NULOUT,'('' ERROR NFRLXG != '',I2)')IFRLXG
    CALL ABOR1('UPDTIM: ABOR1 CALLED')
  ENDIF
  DO JSTEP=1,NFRLXG
    XPRLXG(JSTEP)=0.0_JPRB
  ENDDO
! ISTPF: LENGTH BETWEEN TWO OBSERVATIONS (S)
  ISTPF=IFRHIS*ITIME
! ISTP1: TIME OF THE STEP (S)
  ISTP1=MOD((KSTEP)*ITIME,NINT(RDAY))
  ZDEL=REAL(ISTP1,JPRB)/REAL(ISTPF,JPRB)-REAL(ISTP1/ISTPF,JPRB)
  IF (NFRLXG == 2) THEN
    IBASE=1
    XPRLXG(IBASE)=1.0_JPRB-ZDEL
    XPRLXG(IBASE+1)=ZDEL
  ENDIF
  WRITE(NULOUT,'(''TEST: XPRLXG='',2F9.3)')XPRLXG(1),XPRLXG(NFRLXG)
ENDIF


!          2.   PARAMETERS FOR ECMWF-STYLE INTERMITTENT RADIATION
!               -------------------------------------------------

IF (LEPHYS.OR.((LMPHYS.OR.LSIMPH).AND.LRAYFM)) THEN

  CALL GSTATS(1904,0)
  IF (NSMAX >= 63 .AND. LERAD1H) THEN
    ! For resolutions of T63 and higher, have the option of changing
    ! the frequency of radiation calls within the run
    IF (KSTEP*TSTEP  <  NLNGR1H*RHOUR) THEN
      NRADFR=NRADSFR
    ELSE
      NRADFR=NRADNFR
    ENDIF
  ENDIF

!-- For EPS (characterized by NRADELG > 0), possibility of having different
!   frequency for full radiation computations during the first part of the
!   forecast (of length NRADELG hours)
  IF (NRADELG > 0) THEN
    IF (KSTEP*TSTEP < NRADELG*RHOUR) THEN
      NRADFR=NRADE1H
    ELSE
      NRADFR=NRADE3H
    ENDIF
  ENDIF

  ITIME=NINT( TSTEP)
  IF (YDDYNA%LTWOTL) THEN
    IZT=NINT( TSTEP*(REAL(KSTEP,JPRB)+0.5_JPRB))
  ELSE
    IZT=ITIME*KSTEP
  ENDIF

  !--
  IF (LPERPET) THEN
    ISEC=IZT/NINT(RDAY)
    IZT=IZT-ISEC*NINT(RDAY)
  ENDIF
  !--

  ! Set the time used for computing solar zenith angle (SZA) in the
  ! shortwave radiation scheme to be half a radiation timestep into
  ! the future.  Note that the shortwave scheme produces fluxes
  ! normalized by the TOA incoming solar radiation, so SZA here is
  ! used to calculate the path length of the direct beam through the
  ! atmosphere only.
  IF (LCENTREDTIMESZA) THEN
    ! IZT is the number of seconds into the forecast and is already
    ! half a model timestep ahead.  For radiation every timestep
    ! (NRADFR=1), we don't want to modify this, hence the -1 below.
    ZSTATI=REAL(IZT,JPRB)+0.5_JPRB*(NRADFR-1)*ITIME
  ELSE
    ! The older scheme adds half a radiation timestep to a time that
    ! is already half a model timestep ahead.
    ZSTATI=REAL(IZT,JPRB)+0.5_JPRB*NRADFR*ITIME
  ENDIF
  ISTADD=IZT/NINT(RDAY)
  ISTASS=MOD(IZT,NINT(RDAY))
  ZTIMTR=RTIMST+ZSTATI
  ZHGMT=REAL(MOD(NINT(ZSTATI)+NSSSSS,NINT(RDAY)),JPRB)

  ! Updates concentrations for uniformly mixed gases and solar
  ! irradiance if required.
  IUPGHG=0

  ! In cycles 45R1 and earlier, UPDRGAS was only called every
  ! radiation timestep (via the following IF statement having also
  ! .AND.MOD(KSTEP,NRADFR)==0). However, the solar irradiance RSOLINC
  ! is used for heating-rate calculations every model timestep in
  ! RADHEATN, and a problem can arise in adjoint calculations if
  ! RSOLINC is not updated properly because the backward pass does not
  ! start on a radiation timestep. Therefore we now call it every
  ! model timestep. This leads to UPDRGAS doing a redundant update of
  ! the multipliers for trace gas concentrations (only used within the
  ! radiation scheme), but this is very cheap.  Better for the future
  ! would be to separate UPDRGAS into UPDATE_TRACE_GASES and
  ! UPDATE_SOLAR_IRRADIANCE.
  IF (LHGHG) THEN
    CALL UPDRGAS(YDDYNA,YDERAD,YDERDI,YDRIP,PI05)
    IF (NHINCSOL /= 0) THEN
      ZI0=RSOLINC
    ELSE
      IF (PRESENT(PI05)) THEN
        ZI0=PI05
      ELSE
        ZI0=RI0
      ENDIF
    ENDIF
    RIP0=ZI0*ZSEASON
    RII0 = RIP0
    IUPGHG=1
    ICNT=1
    ISTEP=ICNT*INT(KSTEP/ICNT)
    IF (ISTEP == KSTEP) THEN
      WRITE(NULOUT,FMT='(''ISTEP='',I8,'': Reset RSOLINC='',F10.4,&
        & '' RIP0='',F10.4)') ISTEP,RSOLINC,RIP0
    ENDIF
  ENDIF
  IF (LASTRF) THEN
    CALL UPDCAL (IJ0,IM0,IA0,ISTADD,IJOUR,IMOIS,IAN,ILMOIS,NULOUT)
    IAN = MOD(IAN,4)+2000
    ISS = INT(ZHGMT)
    ZTI=RTIME(IAN,IMOIS,IJOUR,ISS,RDAY)
    ZTETA=RTETA(ZTI)
  ELSE
    ZTETA=RTETA(ZTIMTR)
  ENDIF

  IF( LAPE ) THEN
    ZDEASOM=RRSAQUA(ZTETA)
    ZDECLIM=RDSAQUA(ZTETA)
    ZEQTIMM=RETAQUA(ZTETA)
  ELSE
    ZDEASOM=RRS(ZTETA)
    ZDECLIM=RDS(ZTETA)
    ZEQTIMM=RET(ZTETA)
  ENDIF
  ! Unused so commented out
  !  IF (NHINCSOL /= 0) THEN
  !    ZI0=RSOLINC
  !  ELSE
  !    ZI0=RI0
  !  ENDIF

  IF (.NOT. LPHYLIN) THEN
    ! The following variables are used to compute the solar zenith
    ! angle for use by the shortwave radiation scheme.  Since the
    ! transmissivities from the SW scheme are treated as constant for
    ! NRADFR timesteps, the time used to compute these variables is
    ! actually NRADFR/2 into the future. If LMannersSwUpdate is true
    ! then at every timestep, the subroutine RADHEATN makes an
    ! approximate correction for the change in path length with solar
    ! zenith angle.  To do this it needs to know the cosine of the
    ! solar zenith angle that was used by the radiation scheme, PMU0M.
    ! Since this variable is recomputed every timestep in EC_PHYS, it
    ! is important to modify the following variables only if we are at
    ! a radiation time step.
    IF ( (.NOT.LMANNERSSWUPDATE) .OR. MOD(KSTEP,NRADFR) == 0) THEN
      RSOVRM =ZEQTIMM+ZHGMT
      RWSOVRM=RSOVRM*2.0_JPRB*RPI/RDAY
      RCODECM=COS(ZDECLIM)
      RSIDECM=SIN(ZDECLIM)
      RCOVSRM=COS(RWSOVRM)
      RSIVSRM=SIN(RWSOVRM)
    ENDIF
  ELSE
    ! For the adjoint test to be passed, the simplified physics
    ! requires these variables to be updated every model timestep,
    ! even though that means they won't correctly indicate the values
    ! used in the most recent call to the radiation scheme
    RSOVRM =ZEQTIMM+ZHGMT
    RWSOVRM=RSOVRM*2.0_JPRB*RPI/RDAY
    RCODECM=COS(ZDECLIM)
    RSIDECM=SIN(ZDECLIM)
    RCOVSRM=COS(RWSOVRM)
    RSIVSRM=SIN(RWSOVRM)
  ENDIF

  ZTHETOZ=ZTETAR
  ZANGOZC=REL(ZTHETOZ)-1.7535_JPRB
  ! Unused so commented out
  !  ZCOTHOZ=COS(ZANGOZC)
  !  ZSITHOZ=SIN(ZANGOZC)

  CALL SUECOZO (ZANGOZC, YDRIP%YREOZOC)
! Fortuin-Langematz O3 CLIMATOLOGY:
  IF (NOZOCL == 1.AND. .NOT. LRAYFM15 .AND. .NOT.LPHYLIN) THEN
    IF(MOD(KSTEP,NRADFR) == 0) THEN
      IMINUT=NINT((ZSTATI + REAL(NSSSSS,JPRB))/60._JPRB)
      IF( LAPE ) THEN
        CALL UPDECOZCAQUA(YDERAD,YDERDI,YDRIP%YREOZOC,NINDAT,IMINUT)
      ELSEIF (YDRIP%YRECMIP%NO3CMIP /= 0) THEN
        CALL UPDECOZV(YDERDI,YDRIP%YRECMIP,NINDAT,IMINUT)
      ELSE
        !!CALL SUECOZC(YDERAD,YDERDI,YDRIP%YREOZOC,NINDAT,IMINUT)
        CALL UPDECOZC(YDERAD,YDERDI,YDRIP%YREOZOC,NINDAT,IMINUT)
      ENDIF
    ENDIF
  ENDIF

! Greenhouse gas climatologies  (CO2, O3, CH4, N2O, NO2, CFC11, CFC12, CFC22, CCl4)
  ICNT=100
  IF (NGHGRAD /= 0 .AND. .NOT.LRAYFM15 .AND. .NOT.LPHYLIN) THEN
    IF(MOD(KSTEP,NRADFR) == 0) THEN
      IMINUT=NINT((ZSTATI + REAL(NSSSSS,JPRB))/60._JPRB)
      CALL UPD_GHGCLIM(YDERAD,YDERDI,NINDAT,IMINUT,KSTEP,IUPGHG)
      ISTEP=INT(KSTEP/ICNT)*ICNT
      IF(ISTEP == KSTEP) THEN
        WRITE(NULOUT,FMT='(1X,''STEP='',I8,'' UPDTIM: IUPGHG='',I4,'' RCARDI='',E12.5)') KSTEP,IUPGHG, RCARDI
      ENDIF
    ENDIF
  ENDIF

! Aerosol climatology Tegen et al. / GISS Volcanic aerosol climatology
! and  O.Boucher sulphate history (obs 1920-1990, A1B scenario 2000-2100)
  IF (LHVOLCA .OR. LNEWAER) THEN
    IF ((L_OOPS .AND. LLUPDECAEC) .OR. .NOT.LPHYLIN) THEN
      IF(MOD(KSTEP,NRADFR) == 0) THEN
        IMINUT=NINT((ZSTATI + REAL(NSSSSS,JPRB))/60._JPRB)
        CALL UPDECAEC(YDGEOMETRY,YDERAD,YDRIP,NINDAT,IMINUT)
        IF (LESO4HIS) THEN
          CALL SUECSO4 ( YDRIP%RAERSO4, YDRIP%YRECMIP%NCMIPFIXYR,NINDAT, IMINUT )
        ENDIF
      ENDIF
    ENDIF
  ENDIF
  CALL GSTATS(1904,1)

  IF (LRMU0M) THEN
    ZSTATI=REAL(IZT,JPRB)+NRADFR*ITIME
    ZTIMTR=RTIMST+ZSTATI
    ZHGMT=REAL(MOD(NINT(ZSTATI)+NSSSSS,NINT(RDAY)),JPRB)

    ZTETA=RTETA(ZTIMTR)
    IF( LAPE ) THEN
      ZDEASOM=RRSAQUA(ZTETA)
      ZDECLIM=RDSAQUA(ZTETA)
      ZEQTIMM=RETAQUA(ZTETA)
    ELSE
      ZDEASOM=RRS(ZTETA)
      ZDECLIM=RDS(ZTETA)
      ZEQTIMM=RET(ZTETA)
    ENDIF
    ZSOVRF =ZEQTIMM+ZHGMT
    ZWSOVRF=ZSOVRF*2.0_JPRB*RPI/RDAY

    RCODECF=COS(ZDECLIM)
    RSIDECF=SIN(ZDECLIM)

    RCOVSRF=COS(ZWSOVRF)
    RSIVSRF=SIN(ZWSOVRF)
  ENDIF

ELSEIF (LRAYFM15) THEN
  CALL ABOR1('RADIATION SCHEME FROM CYCLE 15 NO LONGER AVAILABLE')
ENDIF


!          2.5    PARAMETERS FOR POTENTIAL VOLCANIC ERUPTIONS
!                 -------------------------------------------

IF (NAERO /=0 .AND. NVOLERUP == 2) THEN
!--NB: test for only the April 2010 Icelandic erupting volcano
  IF (NINDAT >= 20100414 .AND. NINDAT <= 20100524) THEN
    IMINUT=NINT((ZSTATI + REAL(NSSSSS,JPRB))/60._JPRB)

    CALL SU_AERVOLE(YDEAERVOL,YDRIP,NINDAT,IMINUT)
  ENDIF
ENDIF


!          3.     PARAMETERS USED FOR CLIMATE-TYPE RUNS
!                 -------------------------------------

! Update every time new LBC exits and from start
IF (LMCC01_MSE) THEN
  IF( MOD( RSTATI-ITIME/2,YDMODEL%YRML_LBC%TEFRCL) == 0 )THEN
!   only SST for the moment
    IGP_MSE=1
    CALL UPDCLI_MSE(YDGEOMETRY,YDRIP,IGP_MSE,YDMODEL%YRML_LBC%TEFRCL)
  ENDIF
ENDIF

IF(.NOT.LSCMEC) THEN

!        Updates climatology Meteo-France style

  IF(NSTASS < ITIME)THEN
    IF(LMCC01)THEN
!     basic fields for the climatology
      IGP=YSD_VFD%NUMFLDS+2*YSD_VPD%NUMFLDS+YSD_VAD%NUMFLDS
      IF(LMCC03)THEN
!                  coupled SST and Sea-Ice Mask and albedo
        IGP=IGP+3+NTSSG
      ELSE
!                  clim SST1+SST2 (2 consecutive months)
        IGP=IGP+2
      ENDIF
      IF(LSOLV)THEN
!                  dom.veget+min.sto.res+clay+sand+depth+LAI1+LAI2+(Z0 therm)
        IGP=IGP+YSD_VVD%NUMFLDS
      ENDIF
      IF(LOZONE)THEN
!                  7 basic fields
        IGP=MAX(IGP,14)
      ENDIF
      CALL UPDCLI(YDGEOMETRY,YDSURF,YDMODEL,IGP)
    ENDIF
  ENDIF

! Nudging of upper-air fields

  IF(LNUDG)THEN
    IF(NSTASS < ITIME .AND. NFNUDG == 7 )THEN
      IGP=NTOTFNUDG3*NFLEVG+NTOTFNUDG2
      CALL UPDNUD(YDGEOMETRY,YDRIP,IGP)
    ENDIF
    IF (MOD(KSTEP,NFRNUDG) == 0 .AND. NFNUDG == 2 ) THEN
      IGP=NTOTFNUDG3*NFLEVG+NTOTFNUDG2
      CALL UPDNUD(YDGEOMETRY,YDRIP,IGP)
    ENDIF
  ENDIF

  IF (LRLXG) THEN
    IRLXI=NINT(REAL(NFRLXU)/(TSTEP/3600.0))
    IF (MOD(KSTEP,IRLXI) == 0 .AND. NFRLXG == 2 ) THEN
      WRITE(NULOUT,*) 'Here we read data soon, KSTEP='',1X,I6)',KSTEP
      CALL FLUSH(NULOUT)
      CALL UPDRLXREF(YDGEOMETRY,YDRIP)
      WRITE(NULOUT,*) 'UPDRLXREF sucessfully finished'
      CALL FLUSH(NULOUT)
    ENDIF
  ENDIF

  IF (NFRCPL>0) THEN
    IF ((LMCC03.OR.LMCC05).AND.MOD(KSTEP,NFRCPL) == 0) THEN
      IF(LMCC03)THEN
        IGP=6+NTSSG
      ELSE
        IGP=6
      ENDIF
      CALL UPDCPL(YDGEOMETRY,YDSURF,YDMODEL%YRML_AOC,YDRIP,YDMODEL%YRML_PHY_MF%YRPHY1,IGP)
    ENDIF
  ENDIF

#ifdef WITH_OASIS

  CALL UPDCLIE_OASIS(YDGEOMETRY,YDSURF,YDMCC,YDRIP,YDERAD,YDDYNA,PTSTEP)

#else

  IF (LMCCEC.AND.LDCLUPD) THEN

    ! For coupled integrations, update SSTs at the frequency of
    ! the ocean-atmosphere coupling instead of once a day.
    IF (LMCC04.AND.(.NOT.LNEMO1WAY)) THEN
      ISTASS0=MOD(IZT,NFRCO*ITIME)
      IF(ISTASS0 < ITIME)THEN
        CALL UPDCLIE(YDGEOMETRY,YDDYNA,YDSURF,YDMODEL%YRML_AOC,YDERAD,YDEPHY,YDMODEL%YRML_GCONF,PTSTEP)
      ENDIF
    ENDIF

    !  For coupling to an ocean mixed layer model.
    !  UPDCLIE is called at every timestep to update SST.
    IF (LEOCML .OR. LOCMLTKE) THEN
      CALL UPDCLIE(YDGEOMETRY,YDDYNA,YDSURF,YDMODEL%YRML_AOC,YDERAD,YDEPHY,YDMODEL%YRML_GCONF,PTSTEP)
    ENDIF

  ENDIF

#endif

ENDIF ! .NOT. LSCMEC

IF(ISTASS < ITIME)THEN

!        Updates daily emissions and dry dep velocities
  IF(.NOT.LSCMEC) THEN
    IF ( LMCC_COMPO ) THEN
      CALL UPDCLIE_COMPO(YDGEOMETRY,YDDYNA,YDSURF,YDMODEL%YRML_CHEM%YRCOMPO,YDMCC,YDERAD,YDMODEL%YRML_GCONF,PTSTEP)
    ENDIF
  ENDIF

!        Updates climatology (SST) ECMWF style

  IF(.NOT.LSCMEC) THEN
!!!    IF (.NOT.(LMCC04.AND.(.NOT.LNEMO1WAY)) .AND. .NOT.LEOCML .AND. .NOT.LOCMLTKE) THEN
!!! for the Ocean TKE, I have assume that it will not represent the seasonal change of the foundation SST
!!! and therefore I decided to improse the seasonal change
    IF (.NOT.(LMCC04.AND.(.NOT.LNEMO1WAY)) .AND. .NOT.LEOCML) THEN
      IF (LMCCEC.AND.LDCLUPD) THEN
        CALL UPDCLIE(YDGEOMETRY,YDDYNA,YDSURF,YDMODEL%YRML_AOC,YDERAD,YDEPHY,YDMODEL%YRML_GCONF,PTSTEP)
      ENDIF
    ENDIF
  ENDIF

!        Updates ozone chemistry

  IF (.NOT. LSCMEC) THEN
    IF (LEO3CH.AND.KSTEP >= 0) THEN
      CALL UPDO3CH(YDGEOMETRY,YDDYNA,YDMODEL%YRML_CHEM%YROZO,YDEPHY,YDDPHY,YDRIP)
    ENDIF
  ENDIF

!        Updates concentrations for uniformly mixed gases

  IF (LHGHG.AND.KSTEP >= 0) THEN
    CALL UPDRGAS(YDDYNA,YDERAD,YDERDI,YDRIP,PI05)
  ENDIF

ENDIF

!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('UPDTIM',1,ZHOOK_HANDLE)
END SUBROUTINE UPDTIM
