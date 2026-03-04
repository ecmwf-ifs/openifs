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

SUBROUTINE CPEDIA(YDVAB,YDECUMF,YDML_GCONF,KPROMA,KIDIA,KFDIA,KFLEVG,KSATSIM,&
 & LDFSTEP,LDMLPP,LDRS6,LDRSACC,LDACCRS,&
 & PGELAM ,PGELAT ,PFPLCL ,PFPLCN ,PFPLSL ,PFPLSN ,&
 & PFRSO  ,PFRSOD ,PFRSODC,PFRTH  ,PFRTHD ,PFRTHDC,PUVDF ,PPARF ,&
 & PCAPE  ,PPARCF ,PTINCF ,PFDIR  ,PDSRP  ,PCDIR  ,&
 & PCBASE ,PPDEPL ,P0DEGL ,PM10DEGL,PVISIH,PCIN   ,PCONVIND,&
 & PLIGH_TOT ,PLIGH_CTG,&
 & PCBASEA,PCTOPC ,PTROPOTP,PZTWETB,&
 & PDIFTS ,PDIFTQ ,PSTRTU ,PSTRTV,&
 & PFCLL  ,PFCLN  ,PTCLS  ,PQCLS  ,PUCLS  ,PVCLS ,&
 & PTRENU ,PQLINU ,PLSM   ,PCI    ,PISUND,&
 & PVDIS  ,PVDISG ,PUSTRG,&
 & PVSTRG ,PRUISS ,PRUISP,&
 & PINEE  ,PIGPP  ,PIREC,&
 & PRESH  ,PRESF  ,PDELP,&
 & PU     ,PV     ,PT     ,PQ     ,PO3    ,PGHG   ,PCHEM ,PAERO,&
 & PA     ,PL     ,PR     ,PI     ,PS     ,PTENQ ,POROG,&
 & PDLSP  ,PDCP   ,PDSF   ,PDFZRA ,PDBLD  ,PSUND,&
 & PITTRC ,PITSRC ,PISTRC ,PISSRC ,PIES   ,PISMLT,PI10FG, PILSPF,&
 & PPRECTYPE      ,PFZRA  ,PIVIMD ,PEVAPMU,&
 & PVTTRC ,PVTSRC ,PVSTRC ,PVSSRC ,PVES, PVSMLT, PV10FG, PV10FG6,PVI10FG,&
 & PVLSPF ,PVVIMD ,&
 & PDSSHF ,PDSLHF ,PDSP   ,PDMSL  ,&
 & PDSSR  ,PDSTR  ,PDSSRD ,PDSTRD ,PDSSRDC,PDSTRDC,PDTSR  ,PDTTR,&
 & PDUVDF ,PDPARF ,PDCAPE ,PDCAPES,PDMXCAP6,PDMXCAPS6,&
 & PDPARCF,PDTINCF,PDFDIR,PDDSRP,PDCDIR,&
 & PDCBASE,PD0DEGL,PDM10DEGL,PDVISIH,PDCIN ,PDKINDEX,PDTTINDEX,&
 & PDPDEPL,PDMUCAPE,PDMLCAPE50,PDMLCAPE100,PDMLCIN50,PDMLCIN100,&
 & PDCBASEA,PDCTOPC,PDTROPOTP, PDZTWETB0, PDZTWETB1,&
 & PDEWSS ,PDNSSS ,PDE    ,PDPEV   ,PDLGWS,&
 & PDMGWS ,PDGWD  ,PDMX2T ,PDMN2T ,PDMX2T6,PDMN2T6,PDRO,PDSSRO,PDSRO,&
 & PDNEE  ,PDGPP  ,PDREC,&
 & PDCSF  ,PDLSSF ,PDMXTPR ,PDMNTPR ,PDMXTPR6,PDMNTPR6,&
 & PDTPR  ,PDLSRR ,PDCRR  ,PDLSSFR,PDCSFR,PDPTYPE,PVILSPF,&
 & PDIEWSS,PDINSSS,PDISSHF,PDIE   ,PDVIWVE,PDVIWVN,&
 & PDTCW  ,PDTCWV ,PDTCLW ,PDTCIW ,PDTCRW ,PDTCSW ,PDTCSLW, PDTCO3,&
 & PDTCGHG,PDTCCHEM,PDTCAERO,&
 & PDLITOTI,PDLITOTA6,PDLICGI,PDLICGA6,&
 & PDPTYPEOCC6, PDCLBT ,PDCSBT)

!**** *CPEDIA* - DIAGNOSTIC FIELDS FROM ECMWF PHYSICS

!     Purpose.
!     --------

!           COMPUTE THE DIAGNOSTIC FIELDS FOR FROM ECMWF PHYSICS.

!**   Interface.
!     ----------
!        *CALL* *CPEDIA(...)*

!        Explicit arguments :
!        --------------------
!        YDVAB    : vertical geometry
!        KPROMA   : dimension
!        KIDIA    : start of horizontal loop
!        KFDIA    : end of horizontal loop
!        KSATSIM  : number of simulated satellite images
!        LDFSTEP  : .true. if first step
!        LDMLPP   : .true. if first step after a write-out
!        LDRS6    : .true. to reset 6 hourly min/max bins
!        LDRSACC  : .true. to reset accumulated fluxes at PP time
!        LDACCRS  : .true. to reset accumulated fluxes at now
! Inputs for CPEDIA = outputs from APLPAR (same names as in APLPAR):

!        PGELAM   : longitude
!        PGELAT   : latitude
!        PFPLCL   : convective precipitation as rain.
!        PFPLCN   : convective precipitation as snow.
!        PFPLSL   : stratiform precipitation as rain.
!        PFPLSN   : stratiform precipitation as snow.
!        PFZRA    : stratiform precipitation as freezing rain.
!        PFRSO    : shortwave radiative flux.
!        PFRTH    : longwave radiative flux.
!        PFRSOD   : downward shortwave radiative flux at surface.
!        PFRSODC  : clear-sky downward showrtwave flux at surface
!        PFRTHD   : downward longwave radiative flux at surface.
!        PFRTHDC  : clear-sky downward longwave flux at surface
!        PUVDF    : downward UV flux at surface
!        PPARF    : downward PAR flux at surface
!        PCAPE    : CAPE array of different diagnostics
!        PPDEPL   : departure level for most unstable CAPE
!        PPARCF   : clear-sky downward PAR flux at surface
!        PTINCF   : TOA inclident solar radiation
!        PFDIR    : total sky surface direct SW radiation
!        PCDIR    : clear-sky surface direct SW radiation
!        PDSRP    : total sky surface direct beam SW radiation
!        PCBASE   : cloud base height
!        P0DEGL   : zero deg level
!        PM10DEGL : -10 deg level
!        PVISIH   : horizontal visibility
!        PCIN     : convective inhibition array
!        PCONVIND : array with convective indices
!        PLIGH_TOT: total lightning flash density
!        PLIGH_CTG: cloud-to-ground lightning flash density
!        PCBASEA  : cloud base height aviation
!        PCTOPC   : convective cloud top
!        PTROPOTP : pressure at thermal tropopause
!        PZTWETB  : height of 0 and 1 deg wet bulb T
!        PDIFTS   : sensible heat flux from vertical diffusion
!        PDIFTQ   : moisture flux from vertical diffusion
!        PPOTEVAP : surface potential evaporation
!        PSTRTU   : flux of u momentum from vertical diffusion
!        PSTRTV   : flux of v momentum from vertical diffusion
!        PFCLL    : latent heat flux corresponding to evaporation
!        PFCLN    : latent heat flux corresponding to sublimation
!        PTCLS    : temperature at 2 meters (diagnostic)
!        PQCLS    : humidity at 2 meters (diagnostic)
!        PUCLS    : U wind component at 10 meters (diagnostic)
!        PVCLS    : V wind component at 10 meters (diagnostic)
!        PTRENU   : Skin temperature
!        PQLINU   : Soil water content
!        PLSM     : Land-sea mask
!        PCI      : Sea ice
!        PRESH    : half level pressure
!        PRESF    : full level pressure
!        PDELP    : pressure difference across layers
!        PU       : zonal wind
!        PV       : meridional wind
!        PT       : temperature
!        PQ       : Q
!        PO3      : ozone mixing ratio
!        PGHG(NGHG): GHG mixing ratio
!        PCHEM(NCHEM_TC): CHEM mixing ratio
!        PAERO(NACTAERO): AERO mixing ratio
!        PL       : Cloud liquid water
!        PR       : Rain liquid water
!        PI       : Cloud ice
!        PS       : Cloud snow
!        PPRECTYPE: Precipitation type
!        PEVAPMU  : Potential evaporation
!        POROG    : orography
!        PITTRC   : top thermal radiation clear sky
!        PITSRC   : top solar radiation clear sky
!        PISTRC   : surface thermal radiation clear sky
!        PISSRC   : surface solar radiation clear sky
!        PIES     : snow evaporation
!        PISMLT   : snow melt
!        PI10FG   : gust a 10 m
!        PILSPF   : large scale precipitation fraction
!        PIVIMD   : vertically-integrated mass divergence
! Output
!        PDLSP    : large scale precipitation             (output)
!        PDCP     : convective precipitation              (output)
!        PDSF     : snow fall                             (output)
!        PDFZRA   : freezing rain accumulation            (output)
!        PDBLD    : boundary layer dissipitaion           (output)
!        PDSSHF   : surface sensible heat flux            (output)
!        PDSLHF   : surface latent heat flux              (output)
!        PDSP     : surface pressure                      (output)
!        PDMSL    : mean sea level pressure               (output)
!        PDSSR    : surface solar radiation               (output)
!        PDSTR    : surface thermal radiation             (output)
!        PDSSRD   : surface solar radiation downwards     (output)
!        PDSTRD   : surface thermal radiation downwards   (output)
!        PDSSRDC  : clear-sky surface solar radiation downwards     (output)
!        PDSTRDC  : clear-sky surface thermal radiation downwards   (output)
!        PDTSR    : top solar radiation                   (output)
!        PDTTR    : top thermal radiation                 (output)
!        PDUVDF   : surface UV                            (output)
!        PDPARF   : surface PAR                           (output)
!        PDCAPE   : CAPE                                  (output)
!        PDMUCAPE : most unstable CAPE                    (output)
!        PDMLCAPE50: CAPE from 50 hPa mixed layer         (output)
!        PDMLCAPE100:CAPE from 100 hPa mixed layer        (output)
!        PDCAPES  : CAPE-Shear                            (output)
!        PDMXCAP6 : max CAPE since last 6 hours           (output)
!        PDMXCAPS6: max CAPE-Shear since last 6 hours     (output)
!        PDPARCF  : surface clear-sky PAR                 (output)
!        PDTINCF  : TOA incident solar radiation          (output)
!        PDCBASE  : cloud base height                     (output)
!        PD0DEGL  : zero deg level                        (output)
!        PDVISIH  : horizontal visibility                 (output)
!        PDCIN    : convective inhibition                 (output)
!        PDMLCIN50: convective inhibition 50 hPa mixlayer (output)
!        PDMLCIN100:convective inhibition 100 hPa mixlayer(output)
!        PDKINDEX : convective K-Index                    (output)
!        PDTTINDEX: convective TT-Index                   (output)
!        PDCBASEA : ceiling cloud base heightn            (output)
!        PDCTOPC  : convective cloud top                  (output)
!        PDTROPOTP: pressure at thermal tropopause        (output)
!        PDZTWETB0: height of 0 deg wet bulb T            (output)
!        PDZTWETB1: height of 0 deg wet bulb T            (output)
!        PDFDIR   : total sky surface direct SW radiation (output)
!        PDCDIR   : clear-sky surface direct SW radiation (output)
!        PDDSRP   : total sky surface direct beam SW radiation (output)
!        PDEWSS   : U-sterss                              (output)
!        PDNSSS   : V-stress                              (output)
!        PDE      : evaporation                           (output)
!        PDPEV    : potential evaporation                 (output)
!        PDLGWS   : lat. comp. of gravity wave stress     (output)
!        PDMGWS   : mer. comp. of gravity wave stress     (output)
!        PDGWD    : gravity wave dissipation              (output)
!        PDMX2T   : max temp. at 2 m since prev. p.p.     (output)
!        PDMN2T   : min temp. at 2 m since prev. p.p.     (output)
!        PDMX2T6  : max temp. at 2 m since last 6 hours   (output)
!        PDMN2T6  : min temp. at 2 m since last 6 hours   (output)
!        PDRO     : runoff                                (output)
!        PDSRO    : Surface runoff                        (output)
!        PDSSRO   : Sub-surface runoff                    (output)
!        PDCSF    ; convective snow fall                  (output)
!        PDLSSF   : large scale snow fall                 (output)
!        PDMXTPR  : max precip rate since prev postproces (output)
!        PDMNTPR  : min precip rate since prev postproces (output)
!        PDMXTPR6 : max precip rate in last 6 hours       (output)
!        PDMNTPR6 : min precip rate in last 6 hours       (output)
!        PDTPR    : total precipiation rate               (output)
!        PDLSRR   : large scale rainfall rate             (output)
!        PDCRR    : convective rainfall rate              (output)
!        PDLSSFR  : large scale snowfall rate             (output)
!        PDCSFR   : convective snowfall rate              (output)
!        PDPTYPE  : precipitation type                    (output)
!        PVILSPF  : large scale precip fraction (inst)    (output)
!        PDIE     : instantaneous water vapour flux over  (output)
!                   liquid water (or wet soil)
!        PDISSHF  : instantaneous sensible heat flux at   (output)
!                   surface level.
!        PDIEWSS  : instantaneous surface u-stress        (output)
!        PDINSSS  : instantaneous surface v-stress        (output)
!        PDVIWVE  : eastward water vapour flux            (output)
!        PDVIWVN  : northward water vapour flux           (output)
!        PDTCW    : total column water (all phases)       (output)
!        PDTCWV   : total column water vapour             (output)
!        PDTCLW   : total column liquid water             (output)
!        PDTCIW   : total column ice water                (output)
!        PDTCRW   : total column rain water               (output)
!        PDTCSW   : total column snow water               (output)
!        PDTCSLW  : total column supercooled liquid water (output)
!        PDTCO3   : total column ozone                    (output)
!        PDTCGHG  : total column GHG fields               (output)
!        PDTCCHEM : total column CHEM fields              (output)
!        PDTCAERO : total column AERO fields              (output)
!        PDLITOTI : inst. total lightning flash density   (output)
!        PDLITOTA6: average total lightning flash density (output)
!        PDLICGI  : inst. CTG lightning flash density     (output)
!        PDLICGA6 : average CTG lightning flash density   (output)
!        PDPTYPEOCC6 : freq occurrence of each precip type(output)
!        PVTTRC   : top thermal radiation clear sky       (output)
!        PVTSRC   : top solar radiation clear sky         (output)
!        PVSTRC   : surface thermal radiation clear sky   (output)
!        PVSSRC   : surface solar radiation clear sky     (output)
!        PVES     : snow evaporation                      (output)
!        PVSMLT   : snow melt                             (output)
!        PV10FG   : gust a 10 m since prev. p.p.          (output)
!        PV10FG6  : gust a 10 m since last 6 hours        (output)
!        PVI10FG  : "instantaneous" gust at 10 m          (output)
!        PVLSPF   : large scale precipitation fraction (acc) (output)
!        PVVIMD   : vertically-integrated mass divergence (output)
!        PDCLBT   : Cloudy brightness temperature         (output)
!        PDSLBT   : Clear-sky brightness temperature      (output)
!        --------------------

!     Method.
!     -------

!     Externals.
!     ----------

!     Reference.
!     ----------

!     Author.
!     -------
!      M. Hamrud    ECMWF
!      Original     : 93-02-15

!     Modifications.
!     --------------
!      Modified     : 15-10-01 D.Salmond FULLIMP mods
!      Modified     : 02-09-01 JJMorcrette UV, PAR, CAPE
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      Modified     : 04-02-12 A.Tompkins TCLW, TCIW
!      Modified     : 26-Sep-05 S.Serrar Add total column CO2, SF6
!      Modified     : 02-Jun-06 R. Engelen CO2 replaced by generic GHG
!      Modified     : 07-Sep-06 S. Serrar tracers for diagnostics
!      JJMorcrette 20060721 PP of clear-sky PAR and TOA incident solar radiation
!      JJMorcrette 20060807 PP of VIMD
!      Modified     : 17-Jul-07 S.Serrar variable for methane atmospheric sink
!      Modified     : 20061002  JJMorcrette DDH for aerosol physics
!      JJMorcrette    20100212  PP for CBASE, 0DEGL, VISIH
!      JJMorcrette 20091201 Total and clear-sky direct SW radiation at the surface
!      Modified     : 01-Mar-10 R. Forbes TCRW, TCSW
!      Modified     : 25-Mar-10 T. Wilhelmsson Add 6 hourly min/max fields
!      K. Yessad (Jan 2011): remove useless overdimension.
!      Modified     : 17-Apr-11 G.Balsamo/S.Boussetta Added land carbon dioxide
!      P. Bechtold  : 09-Aug-11 add CIN and ConvIndices
!      Modified     : 31-Oct-2011 M. Ahlgrimm Clear-sky downward radiation at surface
!      R. El Khatib : 01-Mar-2012 LFPOS => LECFPOS
!      A. Inness    : 28-Mar-2012 Add total column CHEM fields
!      R. Forbes    : 01-Mar-2014 Add precip rates/type,I10FG,TCSLW,PDPEV
!      T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!      K. Yessad (July 2014): Move some variables.
!      R. Forbes    : 10-Jan-2015 Add freezing rain FZRA
!      P. Lopez     : Nov 2015 Added lightning fields (including time averaging between pp steps)
!      P. Bechtold  : 10-Nov-2015 Add CBASEA, CTOPC, ZTWETB0-1
!      R. Forbes    : 10-Apr-2017 Add total precipitation rate
!      P. Bechtold  : 29-Apr-2017 Add MXCAP6, MXCAPS6
!      P. Bechtold  : 26-Mar-2021 Add TROPOTP
!      R. Forbes    : May-2022    Add precip type occurrence frequency PTYPEOCC6
!     ------------------------------------------------------------------

USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE YOMSTA                 , ONLY : NLEXTRAP
USE YOMVERT                , ONLY : TVAB
USE PARKIND1               , ONLY : JPIM, JPRB
USE YOMHOOK                , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMPPC                 , ONLY : M2DGGP ,NO2DGG, MFPPHY ,NFPPHY
USE YOMCT0                 , ONLY : LECFPOS, NFRPOS
USE YOMCST                 , ONLY : RTT
USE YOECUMF                , ONLY : TECUMF
USE IOSTREAM_MIX           , ONLY : RMDI
USE YOECLDP                , ONLY : NPRECTYPES
USE YOM_GRIB_CODES         , ONLY : NGRBMX2T, NGRBMN2T, NGRB10FG, NGRBMXTPR, NGRBMNTPR, &
 &                                  NGRBCLBT, NGRBCSBT, NGRBCAPES, NGRBMXCAPS6
!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TVAB)        ,INTENT(IN)    :: YDVAB
TYPE(TECUMF)      ,INTENT(IN)    :: YDECUMF
TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(IN):: YDML_GCONF
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEVG
INTEGER(KIND=JPIM),INTENT(IN)    :: KSATSIM
LOGICAL           ,INTENT(IN)    :: LDFSTEP
LOGICAL           ,INTENT(IN)    :: LDMLPP
LOGICAL           ,INTENT(IN)    :: LDRSACC
LOGICAL           ,INTENT(IN)    :: LDACCRS
LOGICAL           ,INTENT(IN)    :: LDRS6
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGELAM(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGELAT(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFPLCL(KPROMA,0:KFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFPLCN(KPROMA,0:KFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFPLSL(KPROMA,0:KFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFPLSN(KPROMA,0:KFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFZRA(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFRSO(KPROMA,0:KFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFRSOD(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFRSODC(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFRTH(KPROMA,0:KFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFRTHD(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFRTHDC(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUVDF(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PPARF(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCAPE(KPROMA,4)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PPDEPL(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PPARCF(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTINCF(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCIN(KPROMA,3), PCONVIND(KPROMA,2)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLIGH_TOT(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCBASE(KPROMA), PVISIH(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: P0DEGL(KPROMA), PM10DEGL(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLIGH_CTG(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCBASEA(KPROMA), PCTOPC(KPROMA), PZTWETB(KPROMA,2)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTROPOTP(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFDIR(KPROMA), PDSRP(KPROMA), PCDIR(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDIFTS(KPROMA,0:KFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDIFTQ(KPROMA,0:KFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSTRTU(KPROMA,0:KFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSTRTV(KPROMA,0:KFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFCLL(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFCLN(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTCLS(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQCLS(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUCLS(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVCLS(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTRENU(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQLINU(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLSM(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCI(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PISUND(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVDIS(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVDISG(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUSTRG(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVSTRG(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRUISS(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRUISP(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PINEE(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PIGPP(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PIREC(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRESH(KPROMA,0:KFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRESF(KPROMA,KFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDELP(KPROMA,KFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PU(KPROMA,KFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PV(KPROMA,KFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KPROMA,KFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQ(KPROMA,KFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PO3(KPROMA,KFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGHG(KPROMA,KFLEVG,YDML_GCONF%YGFL%NGHG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCHEM(KPROMA,KFLEVG,YDML_GCONF%YGFL%NCHEM_TC)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAERO(KPROMA,KFLEVG,YDML_GCONF%YGFL%NACTAERO)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTENQ(KPROMA,KFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PA(KPROMA,KFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PL(KPROMA,KFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PR(KPROMA,KFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PI(KPROMA,KFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PS(KPROMA,KFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: POROG(KPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDLSP(KPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDCP(KPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDSF(KPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDFZRA(KPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDBLD(KPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSUND(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PITTRC(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PITSRC(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PISTRC(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PISSRC(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PIES(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PISMLT(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PI10FG(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PILSPF(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PPRECTYPE(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEVAPMU(KPROMA)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PIVIMD(KPROMA)

REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVTTRC(KPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVTSRC(KPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVSTRC(KPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVSSRC(KPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVES(KPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVSMLT(KPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PV10FG(KPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PV10FG6(KPROMA)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVI10FG(KPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVLSPF(KPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDSSHF(KPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDSLHF(KPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVVIMD(KPROMA)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDSP(KPROMA)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDMSL(KPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDSSR(KPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDSTR(KPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDSSRD(KPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDSTRD(KPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDSSRDC(KPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDSTRDC(KPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDTSR(KPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDTTR(KPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDUVDF(KPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDPARF(KPROMA)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDCAPE(KPROMA)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDPDEPL(KPROMA)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDMUCAPE(KPROMA)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDMLCAPE50(KPROMA)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDMLCAPE100(KPROMA)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDCAPES(KPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDMXCAP6(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDMXCAPS6(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDPARCF(KPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDTINCF(KPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDCIN(KPROMA), PDKINDEX(KPROMA), PDTTINDEX(KPROMA)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDMLCIN50(KPROMA)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDMLCIN100(KPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDCBASEA(KPROMA), PDCTOPC(KPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDCBASE(KPROMA), PDVISIH(KPROMA)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDTROPOTP(KPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PD0DEGL(KPROMA), PDM10DEGL(KPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDZTWETB0(KPROMA), PDZTWETB1(KPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDFDIR(KPROMA), PDDSRP(KPROMA), PDCDIR(KPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDEWSS(KPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDNSSS(KPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDE(KPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDPEV(KPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDLGWS(KPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDMGWS(KPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDGWD(KPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDMX2T(KPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDMN2T(KPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDMX2T6(KPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDMN2T6(KPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDRO(KPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDSRO(KPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDSSRO(KPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDNEE(KPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDGPP(KPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDREC(KPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDCSF(KPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDLSSF(KPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDMXTPR(KPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDMNTPR(KPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDMXTPR6(KPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDMNTPR6(KPROMA)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDTPR(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDLSRR(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDCRR(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDLSSFR(KPROMA)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDCSFR(KPROMA)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDPTYPE(KPROMA)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVILSPF(KPROMA)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDIEWSS(KPROMA)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDINSSS(KPROMA)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDISSHF(KPROMA)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDIE(KPROMA)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDVIWVE(KPROMA)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDVIWVN(KPROMA)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDTCW(KPROMA)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDTCWV(KPROMA)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDTCLW(KPROMA)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDTCIW(KPROMA)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDTCRW(KPROMA)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDTCSW(KPROMA)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDTCSLW(KPROMA)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDTCO3(KPROMA)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDTCGHG(KPROMA,YDML_GCONF%YGFL%NGHG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDTCCHEM(KPROMA,YDML_GCONF%YGFL%NCHEM_TC)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDTCAERO(KPROMA,YDML_GCONF%YGFL%NACTAERO)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDLITOTI(KPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDLITOTA6(KPROMA)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDLICGI(KPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDLICGA6(KPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDPTYPEOCC6(KPROMA,NPRECTYPES)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDCLBT(KPROMA,KSATSIM)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDCSBT(KPROMA,KSATSIM)

!     ------------------------------------------------------------------


INTEGER(KIND=JPIM), EXTERNAL :: ISRCHEQ
REAL(KIND=JPRB) :: ZTB(KPROMA), ZPRESBH(KPROMA), ZPRESBF(KPROMA)
REAL(KIND=JPRB) :: ZTSTAR(KPROMA), ZT0(KPROMA), ZCAPES(KPROMA)
REAL(KIND=JPRB) :: ZWCH(KPROMA,0:KFLEVG)
REAL(KIND=JPRB) :: ZWVLI(KPROMA,KFLEVG)
REAL(KIND=JPRB) :: ZCLD(KPROMA,KFLEVG)
REAL(KIND=JPRB) :: ZTOTPREC
INTEGER(KIND=JPIM) :: IMN2T, IMX2T, I10FG, IMNTPR, IMXTPR, ICAPES, IMXCAPS6
INTEGER(KIND=JPIM) :: ICLBT, ICSBT, IPRECTYPE
INTEGER(KIND=JPIM) :: JK, JROF, JGHG, JCHEM, JAERO

LOGICAL :: LLMN2TPP, LLMX2TPP, LL10FGPP, LLMNTPRPP, LLMXTPRPP, LLCAPES

REAL(KIND=JPRB) :: ZTSTEP,ZTSTEP_IN_HOURS
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "capeshear.intfb.h"
#include "ctstar.intfb.h"
#include "gppwc.intfb.h"
#include "gptco3.intfb.h"
#include "pppmer.intfb.h"
#include "satsim.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CPEDIA',0,ZHOOK_HANDLE)
ASSOCIATE(NACTAERO=>YDML_GCONF%YGFL%NACTAERO, NCHEM_TC=>YDML_GCONF%YGFL%NCHEM_TC, NGHG=>YDML_GCONF%YGFL%NGHG, &
 & YAERO=>YDML_GCONF%YGFL%YAERO, YCHEM=>YDML_GCONF%YGFL%YCHEM, YGHG=>YDML_GCONF%YGFL%YGHG, &
 & YI=>YDML_GCONF%YGFL%YI, YL=>YDML_GCONF%YGFL%YL, YO3=>YDML_GCONF%YGFL%YO3, YR=>YDML_GCONF%YGFL%YR, &
 & YS=>YDML_GCONF%YGFL%YS, &
 & RBASE0=>YDECUMF%RBASE0, RMINCIN=>YDECUMF%RMINCIN, &
 & TSTEP=>YDML_GCONF%YRRIP%TSTEP )
!     ------------------------------------------------------------------

!*       1.    COMPUTE DIAGNOSTIC QUANTITIES.
!              ------------------------------

!*       1.1  ACCUMULATED FIELDS

ZTSTEP=TSTEP
ZTSTEP_IN_HOURS=TSTEP / 3600._JPRB

! Vertically-integrated mass divergence

CALL GPTCO3(KPROMA,KIDIA,KFDIA,KFLEVG,PDELP,PTENQ,PIVIMD)

IF((LDRSACC.AND.(LDFSTEP.OR.LDMLPP)).OR.LDACCRS) THEN
!        reset accumulated fluxes
  DO JROF=KIDIA,KFDIA
    PDCP  (JROF)=(PFPLCL(JROF,KFLEVG)+PFPLCN(JROF,KFLEVG))*ZTSTEP
    PDLSP (JROF)=(PFPLSL(JROF,KFLEVG)+PFPLSN(JROF,KFLEVG))*ZTSTEP
    PDSF  (JROF)=(PFPLCN(JROF,KFLEVG)+PFPLSN(JROF,KFLEVG))*ZTSTEP
    PDFZRA(JROF)= PFZRA(JROF)*ZTSTEP
    PDCSF (JROF)= PFPLCN(JROF,KFLEVG)*ZTSTEP
    PDLSSF(JROF)= PFPLSN(JROF,KFLEVG)*ZTSTEP
    PDSLHF(JROF)=(PFCLL (JROF)+PFCLN(JROF))*ZTSTEP
    PDSSHF(JROF)= PDIFTS(JROF,KFLEVG)*ZTSTEP
    PDE   (JROF)= PDIFTQ(JROF,KFLEVG)*ZTSTEP
    PDPEV (JROF)= PEVAPMU(JROF)      *ZTSTEP
    PDBLD (JROF)= PVDIS (JROF)       *ZTSTEP
    PDEWSS(JROF)= PSTRTU(JROF,KFLEVG)*ZTSTEP
    PDNSSS(JROF)= PSTRTV(JROF,KFLEVG)*ZTSTEP
    PDGWD (JROF)= PVDISG(JROF)       *ZTSTEP
    PDLGWS(JROF)= PUSTRG(JROF)       *ZTSTEP
    PDMGWS(JROF)= PVSTRG(JROF)       *ZTSTEP
    PDRO  (JROF)=(PRUISS(JROF) + PRUISP(JROF)      )*ZTSTEP
    PDSSRO(JROF)=(PRUISP(JROF)       )*ZTSTEP
    PDSRO (JROF)=(PRUISS(JROF)       )*ZTSTEP
    PDNEE (JROF)=(PINEE (JROF)       )*ZTSTEP
    PDGPP (JROF)=(PIGPP (JROF)       )*ZTSTEP
    PDREC (JROF)=(PIREC (JROF)       )*ZTSTEP
    PDSSR (JROF)= PFRSO (JROF,KFLEVG)*ZTSTEP
    PDSSRD(JROF)= PFRSOD(JROF)       *ZTSTEP
    PDSSRDC(JROF)=PFRSODC(JROF)      *ZTSTEP
    PDSTR (JROF)= PFRTH (JROF,KFLEVG)*ZTSTEP
    PDSTRD(JROF)= PFRTHD(JROF)       *ZTSTEP
    PDSTRDC(JROF)=PFRTHDC(JROF)      *ZTSTEP
    PDTSR (JROF)= PFRSO (JROF,0     )*ZTSTEP
    PDTTR (JROF)= PFRTH (JROF,0     )*ZTSTEP
    PSUND (JROF)= PISUND(JROF)       *ZTSTEP
    PDUVDF(JROF)= PUVDF (JROF)       *ZTSTEP
    PDPARF(JROF)= PPARF (JROF)       *ZTSTEP
    PDPARCF(JROF)=PPARCF(JROF)       *ZTSTEP
    PDTINCF(JROF)=PTINCF(JROF)       *ZTSTEP
    PDFDIR(JROF)= PFDIR(JROF)        *ZTSTEP
    PDCDIR(JROF)= PCDIR(JROF)        *ZTSTEP
    PDDSRP(JROF)= PDSRP(JROF)        *ZTSTEP
    PVTTRC(JROF)= PITTRC(JROF)       *ZTSTEP
    PVTSRC(JROF)= PITSRC(JROF)       *ZTSTEP
    PVSTRC(JROF)= PISTRC(JROF)       *ZTSTEP
    PVSSRC(JROF)= PISSRC(JROF)       *ZTSTEP
    PVES  (JROF)= PIES  (JROF)       *ZTSTEP
    PVSMLT(JROF)= PISMLT(JROF)       *ZTSTEP
    PVLSPF(JROF)= PILSPF(JROF)       *ZTSTEP
    PVVIMD(JROF)= PIVIMD(JROF)       *ZTSTEP

  ENDDO
ELSE

  DO JROF=KIDIA,KFDIA
    PDCP  (JROF)=PDCP  (JROF)+(PFPLCL(JROF,KFLEVG)+PFPLCN(JROF,KFLEVG))*ZTSTEP
    PDLSP (JROF)=PDLSP (JROF)+(PFPLSL(JROF,KFLEVG)+PFPLSN(JROF,KFLEVG))*ZTSTEP
    PDSF  (JROF)=PDSF  (JROF)+(PFPLCN(JROF,KFLEVG)+PFPLSN(JROF,KFLEVG))*ZTSTEP
    PDFZRA(JROF)=PDFZRA(JROF)+ PFZRA(JROF)*ZTSTEP
    PDCSF (JROF)=PDCSF (JROF)+ PFPLCN(JROF,KFLEVG)*ZTSTEP
    PDLSSF(JROF)=PDLSSF(JROF)+ PFPLSN(JROF,KFLEVG)*ZTSTEP
    PDSLHF(JROF)=PDSLHF(JROF)+(PFCLL (JROF)+PFCLN(JROF))*ZTSTEP
    PDSSHF(JROF)=PDSSHF(JROF)+ PDIFTS(JROF,KFLEVG)*ZTSTEP
    PDE   (JROF)=PDE   (JROF)+ PDIFTQ(JROF,KFLEVG)*ZTSTEP
    PDPEV (JROF)=PDPEV (JROF)+ PEVAPMU(JROF)      *ZTSTEP
    PDBLD (JROF)=PDBLD (JROF)+ PVDIS (JROF)       *ZTSTEP
    PDEWSS(JROF)=PDEWSS(JROF)+ PSTRTU(JROF,KFLEVG)*ZTSTEP
    PDNSSS(JROF)=PDNSSS(JROF)+ PSTRTV(JROF,KFLEVG)*ZTSTEP
    PDGWD (JROF)=PDGWD (JROF)+ PVDISG(JROF)       *ZTSTEP
    PDLGWS(JROF)=PDLGWS(JROF)+ PUSTRG(JROF)       *ZTSTEP
    PDMGWS(JROF)=PDMGWS(JROF)+ PVSTRG(JROF)       *ZTSTEP
    PDRO  (JROF)=PDRO  (JROF)+(PRUISS(JROF)+ PRUISP(JROF)      )*ZTSTEP
    PDSSRO(JROF)=PDSSRO(JROF)+(PRUISP(JROF)      )*ZTSTEP
    PDSRO (JROF)=PDSRO (JROF)+(PRUISS(JROF)      )*ZTSTEP
    PDNEE (JROF)=PDNEE (JROF)+(PINEE (JROF)      )*ZTSTEP
    PDGPP (JROF)=PDGPP (JROF)+(PIGPP (JROF)      )*ZTSTEP
    PDREC (JROF)=PDREC (JROF)+(PIREC (JROF)      )*ZTSTEP
    PDSSR (JROF)=PDSSR (JROF)+ PFRSO (JROF,KFLEVG)*ZTSTEP
    PDSSRD(JROF)=PDSSRD(JROF)+ PFRSOD(JROF)       *ZTSTEP
    PDSSRDC(JROF)=PDSSRDC(JROF)+PFRSODC(JROF)     *ZTSTEP
    PDSTR (JROF)=PDSTR (JROF)+ PFRTH (JROF,KFLEVG)*ZTSTEP
    PDSTRD(JROF)=PDSTRD(JROF)+ PFRTHD(JROF)       *ZTSTEP
    PDSTRDC(JROF)=PDSTRDC(JROF)+PFRTHDC(JROF)     *ZTSTEP
    PDTSR (JROF)=PDTSR (JROF)+ PFRSO (JROF,0     )*ZTSTEP
    PDTTR (JROF)=PDTTR (JROF)+ PFRTH (JROF,0     )*ZTSTEP
    PSUND (JROF)=PSUND (JROF)+ PISUND(JROF)       *ZTSTEP
    PDUVDF(JROF)=PDUVDF(JROF)+ PUVDF (JROF)       *ZTSTEP
    PDPARF(JROF)=PDPARF(JROF)+ PPARF (JROF)       *ZTSTEP
    PDPARCF(JROF)=PDPARCF(JROF)+PPARCF(JROF)      *ZTSTEP
    PDTINCF(JROF)=PDTINCF(JROF)+PTINCF(JROF)      *ZTSTEP
    PDFDIR(JROF)=PDFDIR(JROF)+ PFDIR(JROF)        *ZTSTEP
    PDCDIR(JROF)=PDCDIR(JROF)+ PCDIR(JROF)        *ZTSTEP
    PDDSRP(JROF)=PDDSRP(JROF)+ PDSRP(JROF)        *ZTSTEP
    PVTTRC(JROF)=PVTTRC(JROF)+ PITTRC(JROF)       *ZTSTEP
    PVTSRC(JROF)=PVTSRC(JROF)+ PITSRC(JROF)       *ZTSTEP
    PVSTRC(JROF)=PVSTRC(JROF)+ PISTRC(JROF)       *ZTSTEP
    PVSSRC(JROF)=PVSSRC(JROF)+ PISSRC(JROF)       *ZTSTEP
    PVES  (JROF)=PVES  (JROF)+ PIES  (JROF)       *ZTSTEP
    PVSMLT(JROF)=PVSMLT(JROF)+ PISMLT(JROF)       *ZTSTEP
    PVLSPF(JROF)=PVLSPF(JROF)+ PILSPF(JROF)       *ZTSTEP
    PVVIMD(JROF)=PVVIMD(JROF)+ PIVIMD(JROF)       *ZTSTEP

  ENDDO

ENDIF

!*       1.2  INSTANTANEOUS FIELDS

LLMX2TPP =.FALSE.
LLMN2TPP =.FALSE.
LL10FGPP =.FALSE.
LLMXTPRPP =.FALSE.
LLMNTPRPP =.FALSE.
LLCAPES = .FALSE.
IF(LDMLPP ) THEN
  IF(LECFPOS) THEN
    IMX2T = ISRCHEQ(NFPPHY,MFPPHY,1,NGRBMX2T)
    IF(IMX2T <= NFPPHY) LLMX2TPP=.TRUE.
    IMN2T = ISRCHEQ(NFPPHY,MFPPHY,1,NGRBMN2T)
    IF(IMN2T <= NFPPHY) LLMN2TPP=.TRUE.
    I10FG = ISRCHEQ(NFPPHY,MFPPHY,1,NGRB10FG)
    IF(I10FG <= NFPPHY) LL10FGPP=.TRUE.
    IMXTPR = ISRCHEQ(NFPPHY,MFPPHY,1,NGRBMXTPR)
    IF(IMXTPR <= NFPPHY) LLMXTPRPP=.TRUE.
    IMNTPR = ISRCHEQ(NFPPHY,MFPPHY,1,NGRBMNTPR)
    IF(IMNTPR <= NFPPHY) LLMNTPRPP=.TRUE.
    ICAPES = ISRCHEQ(NFPPHY,MFPPHY,1,NGRBCAPES)
    IF(ICAPES <= NFPPHY) LLCAPES=.TRUE.
  ELSE
    IMX2T = ISRCHEQ(NO2DGG,M2DGGP,1,NGRBMX2T)
    IF(IMX2T <= NO2DGG) LLMX2TPP=.TRUE.
    IMN2T = ISRCHEQ(NO2DGG,M2DGGP,1,NGRBMN2T)
    IF(IMN2T <= NO2DGG) LLMN2TPP=.TRUE.
    I10FG = ISRCHEQ(NO2DGG,M2DGGP,1,NGRB10FG)
    IF(I10FG <= NO2DGG) LL10FGPP=.TRUE.
    IMXTPR = ISRCHEQ(NO2DGG,M2DGGP,1,NGRBMXTPR)
    IF(IMXTPR <= NO2DGG) LLMXTPRPP=.TRUE.
    IMNTPR = ISRCHEQ(NO2DGG,M2DGGP,1,NGRBMNTPR)
    IF(IMNTPR <= NO2DGG) LLMNTPRPP=.TRUE.
    ICAPES = ISRCHEQ(NO2DGG,M2DGGP,1,NGRBCAPES)
    IF(ICAPES <= NO2DGG) LLCAPES=.TRUE.
  ENDIF
ENDIF

IF(LECFPOS) THEN
   IMXCAPS6 = ISRCHEQ(NFPPHY,MFPPHY,1,NGRBMXCAPS6)
   IF(IMXCAPS6 <= NFPPHY) LLCAPES=.TRUE.
ELSE
   IMXCAPS6 = ISRCHEQ(NO2DGG,M2DGGP,1,NGRBMXCAPS6)
   IF(IMXCAPS6 <= NO2DGG) LLCAPES=.TRUE.
ENDIF

DO JROF=KIDIA,KFDIA

  ! Calculate total precipitation rate and check positive
  ZTOTPREC = MAX(PFPLSL(JROF,KFLEVG) + PFPLSN(JROF,KFLEVG)&
 &             + PFPLCL(JROF,KFLEVG) + PFPLCN(JROF,KFLEVG),0._JPRB)

  IF(LLMX2TPP.OR.LDFSTEP) THEN
    PDMX2T(JROF) = PTCLS(JROF)
  ELSE
    PDMX2T(JROF) = MAX(PTCLS(JROF),PDMX2T(JROF))
  ENDIF
  IF(LLMN2TPP.OR.LDFSTEP) THEN
    PDMN2T(JROF) = PTCLS(JROF)
  ELSE
    PDMN2T(JROF) = MIN(PTCLS(JROF),PDMN2T(JROF))
  ENDIF
  IF(LL10FGPP.OR.LDFSTEP) THEN
    PV10FG(JROF) = PI10FG(JROF)
  ELSE
    PV10FG(JROF) = MAX(PI10FG(JROF),PV10FG(JROF))
  ENDIF
  IF(LLMXTPRPP.OR.LDFSTEP) THEN
    PDMXTPR(JROF) = ZTOTPREC
  ELSE
    PDMXTPR(JROF) = MAX(ZTOTPREC,PDMXTPR(JROF))
  ENDIF
  IF(LLMNTPRPP.OR.LDFSTEP) THEN
    PDMNTPR(JROF) = ZTOTPREC
  ELSE
    PDMNTPR(JROF) = MIN(ZTOTPREC,PDMNTPR(JROF))
  ENDIF

  PDIE(JROF)    = PDIFTQ(JROF,KFLEVG)
  PDISSHF(JROF) = PDIFTS(JROF,KFLEVG)
  PDIEWSS(JROF) = PSTRTU(JROF,KFLEVG)
  PDINSSS(JROF) = PSTRTV(JROF,KFLEVG)
  PDSP(JROF)    =  PRESH(JROF,KFLEVG)
  PVI10FG(JROF) = PI10FG(JROF)

  ! Instantaneous precipitation variables
  PDTPR(JROF)   = PFPLSL(JROF,KFLEVG) + PFPLSN(JROF,KFLEVG) &
     &          + PFPLCL(JROF,KFLEVG) + PFPLCN(JROF,KFLEVG)
  PDLSRR(JROF)  = PFPLSL(JROF,KFLEVG)
  PDLSSFR(JROF) = PFPLSN(JROF,KFLEVG)
  PDCRR(JROF)   = PFPLCL(JROF,KFLEVG)
  PDCSFR(JROF)  = PFPLCN(JROF,KFLEVG)
  PVILSPF(JROF) = PILSPF(JROF)
  PDPTYPE(JROF) = PPRECTYPE(JROF)   

!*       1.3  ACCUMULATED FIELDS.

  ! Special fields for accumulations over past 6 hours (also used for 1h and 3h).
  ! Note: Here, we actually compute hourly accumulations, which will then be further 
  !       accumulated in subroutine HPOS to obtain the desired 1h, 3h and 6h fields.

  IF(LDRS6.OR.LDFSTEP) THEN
    PDMX2T6(JROF)   = PTCLS(JROF)
    PDMN2T6(JROF)   = PTCLS(JROF)
    PV10FG6(JROF)   = PI10FG(JROF)
    PDMXTPR6(JROF)  = ZTOTPREC
    PDMNTPR6(JROF)  = ZTOTPREC
    PDLITOTA6(JROF) = PLIGH_TOT(JROF)*ZTSTEP_IN_HOURS
    PDLICGA6 (JROF) = PLIGH_CTG(JROF)*ZTSTEP_IN_HOURS
    ! Reset precipitation type occurrence counter
    PDPTYPEOCC6(JROF,:) = 0
  ELSE
    PDMX2T6(JROF)   = MAX(PTCLS(JROF),PDMX2T6(JROF))
    PDMN2T6(JROF)   = MIN(PTCLS(JROF),PDMN2T6(JROF))
    PV10FG6(JROF)   = MAX(PI10FG(JROF),PV10FG6(JROF))
    PDMXTPR6(JROF)  = MAX(ZTOTPREC,PDMXTPR6(JROF))
    PDMNTPR6(JROF)  = MIN(ZTOTPREC,PDMNTPR6(JROF))
    PDLITOTA6(JROF) = PDLITOTA6(JROF) + PLIGH_TOT(JROF)*ZTSTEP_IN_HOURS
    PDLICGA6 (JROF) = PDLICGA6 (JROF) + PLIGH_CTG(JROF)*ZTSTEP_IN_HOURS
  ENDIF

  ! Store precipitation type occurrence frequency
  IPRECTYPE = INT(PPRECTYPE(JROF))
  IF (IPRECTYPE > 0) PDPTYPEOCC6(JROF,IPRECTYPE) = PDPTYPEOCC6(JROF,IPRECTYPE) + 1._JPRB
  
ENDDO

!*       1.4  MEAN SEA LEVEL PRESSURE

DO JROF=KIDIA,KFDIA
  ZTB    (JROF) = PT   (JROF,NLEXTRAP)
  ZPRESBH(JROF) = PRESH(JROF,KFLEVG)
  ZPRESBF(JROF) = PRESF(JROF,NLEXTRAP)
ENDDO
CALL CTSTAR(KPROMA,KIDIA,KFDIA,ZTB,ZPRESBH,ZPRESBF,POROG,ZTSTAR,ZT0)
CALL PPPMER(KPROMA,KIDIA,KFDIA,ZPRESBH,POROG,ZTSTAR,ZT0,PDMSL)

!*       1.5  TOTAL COLUMN WATER VAPOUR, CAPE, CLOUD BASE, 0DEGL,  WATER, ICE, LIQUID
!                                         CIN, KINDEX    ,TTINDEX

IF (LLCAPES) THEN
   CALL CAPESHEAR(YDVAB,KPROMA,KIDIA,KFDIA,KFLEVG,PRESH,PRESF,PU,PV,PCAPE(:,2),ZCAPES)
ELSE
   ZCAPES(:)=0._JPRB
ENDIF

DO JROF=KIDIA,KFDIA
   IF(LDRS6.OR.LDFSTEP) THEN
      PDMXCAP6(JROF) = PCAPE(JROF,2)
      IF (LLCAPES) PDMXCAPS6(JROF)= ZCAPES(JROF)
   ELSE
      PDMXCAP6(JROF) = MAX(PCAPE(JROF,2),PDMXCAP6(JROF))
      IF (LLCAPES) PDMXCAPS6(JROF)= MAX(ZCAPES(JROF),PDMXCAPS6(JROF))
   ENDIF
ENDDO

IF(LDMLPP ) THEN
  CALL GPPWC(KPROMA,KIDIA,KFDIA,KFLEVG,ZWCH,PQ,PRESH)
  DO JROF=KIDIA,KFDIA
    PDTCWV (JROF)=ZWCH(JROF,KFLEVG)
    PDCAPE (JROF)=PCAPE(JROF,1)
    PDMUCAPE(JROF)=PCAPE(JROF,2)
    PDMLCAPE50(JROF)=PCAPE(JROF,3)
    PDMLCAPE100(JROF)=PCAPE(JROF,4)
    PDPDEPL(JROF)=PPDEPL(JROF)
    PDCAPES(JROF)=ZCAPES(JROF)
    PDCBASE(JROF)=PCBASE(JROF)
    PD0DEGL(JROF)=P0DEGL(JROF)
    PDM10DEGL(JROF)=PM10DEGL(JROF)
    PDVISIH(JROF)=PVISIH(JROF)
    PDCIN(JROF)  =PCIN(JROF,1)
    PDMLCIN50(JROF)=PCIN(JROF,2)
    PDMLCIN100(JROF)=PCIN(JROF,3)
    PDTTINDEX(JROF)=PCONVIND(JROF,1)
    PDKINDEX(JROF)=PCONVIND(JROF,2)
    PDCBASEA(JROF)=PCBASEA(JROF)
    PDCTOPC(JROF)=PCTOPC(JROF)
    PDTROPOTP(JROF)=PTROPOTP(JROF)
    PDZTWETB0(JROF)=PZTWETB(JROF,1)
    PDZTWETB1(JROF)=PZTWETB(JROF,2)
  ENDDO

  ! Trigger bitmap in GRIB fields
  DO JROF=KIDIA,KFDIA
    IF(PCBASE(JROF) ==  RBASE0)  PDCBASE(JROF) = RMDI
    IF(PCBASEA(JROF)==  RBASE0)  PDCBASEA(JROF)= RMDI
    IF(PCTOPC(JROF) ==  RBASE0)  PDCTOPC(JROF) = RMDI
    IF(PCIN(JROF,1) == -RMINCIN) PDCIN(JROF)   = RMDI
    IF(PCIN(JROF,2) == -RMINCIN) PDMLCIN50(JROF)= RMDI
    IF(PCIN(JROF,3) == -RMINCIN) PDMLCIN100(JROF)= RMDI
  ENDDO

  ! Vertically integrated eastward water vapour flux
  DO JK=1,KFLEVG
    DO JROF=KIDIA,KFDIA
      ZCLD(JROF,JK)=PQ(JROF,JK)*PU(JROF,JK)
    ENDDO
  ENDDO
  CALL GPPWC(KPROMA,KIDIA,KFDIA,KFLEVG,ZWCH,ZCLD,PRESH)
  DO JROF=KIDIA,KFDIA
    PDVIWVE(JROF)=ZWCH(JROF,KFLEVG)
  ENDDO

  ! Vertically integrated northward water vapour flux
  DO JK=1,KFLEVG
    DO JROF=KIDIA,KFDIA
      ZCLD(JROF,JK)=PQ(JROF,JK)*PV(JROF,JK)
    ENDDO
  ENDDO
  CALL GPPWC(KPROMA,KIDIA,KFDIA,KFLEVG,ZWCH,ZCLD,PRESH)
  DO JROF=KIDIA,KFDIA
    PDVIWVN(JROF)=ZWCH(JROF,KFLEVG)
  ENDDO


  ! Total column liquid water content
  IF (YL%LACTIVE) THEN
    DO JK=1,KFLEVG
      DO JROF=KIDIA,KFDIA
        ZCLD(JROF,JK)=PL(JROF,JK)
        IF (ZCLD(JROF,JK)<1.E-12_JPRB) ZCLD(JROF,JK)=0.0_JPRB
      ENDDO
    ENDDO
    CALL GPPWC(KPROMA,KIDIA,KFDIA,KFLEVG,ZWCH,ZCLD,PRESH)
    DO JROF=KIDIA,KFDIA
      PDTCLW(JROF)=ZWCH(JROF,KFLEVG)
    ENDDO
  ENDIF

  ! Total column rain water content
  IF (YR%LACTIVE) THEN
    DO JK=1,KFLEVG
      DO JROF=KIDIA,KFDIA
        ZCLD(JROF,JK)=PR(JROF,JK)
        IF (ZCLD(JROF,JK)<1.E-12_JPRB) ZCLD(JROF,JK)=0.0_JPRB
      ENDDO
    ENDDO
    CALL GPPWC(KPROMA,KIDIA,KFDIA,KFLEVG,ZWCH,ZCLD,PRESH)
    DO JROF=KIDIA,KFDIA
      PDTCRW(JROF)=ZWCH(JROF,KFLEVG)
    ENDDO
  ENDIF

  ! Total column ice water content
  IF (YI%LACTIVE) THEN
    DO JK=1,KFLEVG
      DO JROF=KIDIA,KFDIA
        ZCLD(JROF,JK)=PI(JROF,JK)
        IF (ZCLD(JROF,JK)<1.E-12_JPRB) ZCLD(JROF,JK)=0.0_JPRB
      ENDDO
    ENDDO
    CALL GPPWC(KPROMA,KIDIA,KFDIA,KFLEVG,ZWCH,ZCLD,PRESH)
    DO JROF=KIDIA,KFDIA
      PDTCIW(JROF)=ZWCH(JROF,KFLEVG)
    ENDDO
  ENDIF

  ! Total column snow water content
  IF (YS%LACTIVE) THEN
    DO JK=1,KFLEVG
      DO JROF=KIDIA,KFDIA
        ZCLD(JROF,JK)=PS(JROF,JK)
        IF (ZCLD(JROF,JK)<1.E-12_JPRB) ZCLD(JROF,JK)=0.0_JPRB
      ENDDO
    ENDDO
    CALL GPPWC(KPROMA,KIDIA,KFDIA,KFLEVG,ZWCH,ZCLD,PRESH)
    DO JROF=KIDIA,KFDIA
      PDTCSW(JROF)=ZWCH(JROF,KFLEVG)
    ENDDO
  ENDIF

  ! Total column supercooled liquid water content
  IF (YL%LACTIVE) THEN
    DO JK=1,KFLEVG
      DO JROF=KIDIA,KFDIA
        ZCLD(JROF,JK)=PL(JROF,JK)
        IF (ZCLD(JROF,JK)<1.E-12_JPRB .OR. PT(JROF,JK)>RTT) ZCLD(JROF,JK)=0.0_JPRB
      ENDDO
    ENDDO
    CALL GPPWC(KPROMA,KIDIA,KFDIA,KFLEVG,ZWCH,ZCLD,PRESH)
    DO JROF=KIDIA,KFDIA
      PDTCSLW(JROF)=ZWCH(JROF,KFLEVG)
    ENDDO
  ENDIF

  ! Total column water content (vapour+liquid+ice+rain+snow)
  IF (YL%LACTIVE .AND. YI%LACTIVE .AND. .NOT. YS%LACTIVE .AND.&
   & .NOT. YR%LACTIVE) THEN
    DO JK=1,KFLEVG
      DO JROF=KIDIA,KFDIA
        ZWVLI(JROF,JK)=((PI(JROF,JK))+PL(JROF,JK))+PQ(JROF,JK)
      ENDDO
    ENDDO
    CALL GPPWC(KPROMA,KIDIA,KFDIA,KFLEVG,ZWCH,ZWVLI,PRESH)
  ENDIF

  IF (YL%LACTIVE .AND. YI%LACTIVE .AND. YS%LACTIVE .AND. .NOT. YR%LACTIVE) THEN
    DO JK=1,KFLEVG
      DO JROF=KIDIA,KFDIA
        ZWVLI(JROF,JK)=(((PI(JROF,JK))+PS(JROF,JK))+PL(JROF,JK))+PQ(JROF,JK)
      ENDDO
    ENDDO
    CALL GPPWC(KPROMA,KIDIA,KFDIA,KFLEVG,ZWCH,ZWVLI,PRESH)
  ENDIF

  IF (YL%LACTIVE .AND. YI%LACTIVE .AND. YS%LACTIVE .AND. YR%LACTIVE) THEN
    DO JK=1,KFLEVG
      DO JROF=KIDIA,KFDIA
        ZWVLI(JROF,JK)=((((PI(JROF,JK))+PS(JROF,JK))+PL(JROF,JK))+&
        &                  PR(JROF,JK))+PQ(JROF,JK)
      ENDDO
    ENDDO
    CALL GPPWC(KPROMA,KIDIA,KFDIA,KFLEVG,ZWCH,ZWVLI,PRESH)
  ENDIF

  DO JROF=KIDIA,KFDIA
    PDTCW(JROF)=ZWCH(JROF,KFLEVG)
  ENDDO

  DO JROF=KIDIA,KFDIA
    ! Instantaneous lightning flash densities (total and cloud-to-ground).
    PDLITOTI(JROF)=PLIGH_TOT(JROF)
    PDLICGI (JROF)=PLIGH_CTG(JROF)
  ENDDO

ENDIF

!*       1.6  TOTAL COLUMN OZONE

IF(LDMLPP .AND. YO3%LACTIVE) THEN
  CALL GPTCO3(KPROMA,KIDIA,KFDIA,KFLEVG,PDELP,PO3,PDTCO3)
ENDIF

!*       1.7  TOTAL COLUMN GHG fields

DO JGHG=1,NGHG
  IF(YGHG(JGHG)%LACTIVE) THEN
    CALL GPPWC(KPROMA,KIDIA,KFDIA,KFLEVG,ZWCH,PGHG(:,:,JGHG),PRESH)
  ENDIF
  DO JROF=KIDIA,KFDIA
    PDTCGHG(JROF,JGHG)=ZWCH(JROF,KFLEVG)
  ENDDO
ENDDO

!*       1.8  TOTAL COLUMN CHEM fields

DO JCHEM=1,NCHEM_TC
  IF(YCHEM(JCHEM)%LACTIVE) THEN
    CALL GPPWC(KPROMA,KIDIA,KFDIA,KFLEVG,ZWCH,PCHEM(:,:,JCHEM),PRESH)
  ENDIF
  DO JROF=KIDIA,KFDIA
    PDTCCHEM(JROF,JCHEM)=ZWCH(JROF,KFLEVG)
  ENDDO
ENDDO

!*       1.9  TOTAL COLUMN AERO fields

DO JAERO=1,NACTAERO
  IF(YAERO(JAERO)%LACTIVE) THEN
    CALL GPPWC(KPROMA,KIDIA,KFDIA,KFLEVG,ZWCH,PAERO(:,:,JAERO),PRESH)
  ENDIF
  DO JROF=KIDIA,KFDIA
    PDTCAERO(JROF,JAERO)=ZWCH(JROF,KFLEVG)
  ENDDO
ENDDO

!*       1.10  SIMULATED SATELLITE images

IF (LDMLPP) THEN
  ICLBT = ISRCHEQ(NFPPHY,MFPPHY,1,NGRBCLBT)
  ICSBT = ISRCHEQ(NFPPHY,MFPPHY,1,NGRBCSBT)

  IF (KSATSIM > 0 .AND. (ICLBT <= NFPPHY .OR. ICSBT <= NFPPHY)) THEN
    CALL SATSIM(YDML_GCONF%YGFL,KPROMA,KIDIA,KFDIA,KSATSIM,KFLEVG,&
      & PGELAT,PGELAM,POROG,PLSM,PCI,PTRENU,PQLINU,PCAPE(:,2),PTCLS,PQCLS,PUCLS,PVCLS,PDSP,&
      & PRESF,PT,PQ,PO3,PA,PL,PI,PR,PS,&
      & PDCLBT,PDCSBT)
  ENDIF
ENDIF


!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CPEDIA',1,ZHOOK_HANDLE)
END SUBROUTINE CPEDIA
