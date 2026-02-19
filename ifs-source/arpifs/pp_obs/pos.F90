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

SUBROUTINE POS(YDQTYPE,YDNAMFPSCI,YDTFP,LDHPOS,YDFPVAB,YDGEOMETRY,YDSURF,&
 & YDML_GCONF,YDDYN,YDDYNA,YDPHY,KPROMA,KST,KND,&
 & KOPLEV,KGFL,KFLDSPP,KFLDAUX,KFLEVG,KPPM,CDCONF,YDIN,YDGFL,&
 & PRCORI,PGM,PRATATH,POROG,POROGL,POROGM,&
 & PXLEV,PTSI,PEV,&
 & PU,PUL,PDIV,PV,PVL,PVOR,PT,PTL,PTM,PSPD,PSPDL,PSPDM,PSVD,PNHX,&
 & PSP,PSPL,PSPM,&
 & PGFL,PEXTRA,PEXRAD,PEXTDI,KCUFNR,PRMCUFGP,PGPP,PAUX)

!**** *POS* - post-processing

!     PURPOSE.
!     -------
!           Post-processing calculations.

!**   INTERFACE.
!     ---------
!        *CALL* *POS(...)

!        EXPLICIT ARGUMENTS :
!        ------------------
!        * INPUT:
!        LDHPOS    : .TRUE. if any horizontal gridpoint post-processing (core or E-zone)
!        KPROMA    : horizontal dimension.
!        KST       : first point to post-process
!        KND       : last point to post-process
!        KOPLEV    : number of post-processing levels
!        KGFL      : number of GFL fields in PGFL
!        KFLDSPP   : number of fields to be spectrally fitted
!        KFLDAUX   : number of fields to remain grid point
!        KFLEVG    : number of model vertical levels
!        KPPM      : number of interpolation methods in post-processing (NPPM)
!        CDCONF    : configuration of work
!                    'B' <=> Pressure                 level pp.
!                    'H' <=> Height (above orography) level pp.
!                    'F' <=> aeronautic flight        level pp.
!                    'K' <=> Temperature              level pp.
!                    'T' <=> Potential temperature    level pp.
!                    'V' <=> Potential vorticity      level pp.
!                    'M' <=> Model                    level pp.
!                    'S' <=> Eta                      level pp.
!        YDIN      : pointers allowing to retrieve individual fields in PGFL
!        YDGFL     : says if GFL is active or not
!        PRCORI    : coriolis parameter
!        PGM       : map factor
!        PRATATH   : 2*tan(theta_bne)/RA
!        POROG     : orography
!        POROGL    : zonal gradient of orography
!        POROGM    : meridian gradient of orography
!        PXLEV     : array containing the pp levels (in P,H,PV,TH or s)
!        PTSI      : input surface temperature
!        PEV       : <eta dot * dp/deta> on top of the atmosphere
!        PU        : U-wind
!        PUL       : zonal gradient of U-wind
!        PDIV      : divergence
!        PV        : V-wind
!        PVL       : zonal gradient of V-wind
!        PVOR      : vorticity
!        PT        : temperature
!        PTL       : zonal gradient of temperature
!        PTM       : meridian gradient of temperature
!        PSPD      : pressure departure variable
!        PSPDL     : zonal gradient of pressure departure var
!        PSPDM     : merid gradient of pressure departure var
!        PSVD      : vertical divergence variable
!        PNHX      : term NHX, used only if LNHX=T
!        PSP       : surface pressure variable
!        PSPL      : zonal gradient of surface pressure variable
!        PSPM      : meridian gradient of surface pressure variable
!        PGFL      : GFL fields
!        PEXTRA    : extra input fields
!        PEXRAD    : extra input fields from radiation code
!        PRMCUFGP  : CUF input fields

!        * OUTPUT (post-processed fields):
!        PGPP      : fields to be spectrally fitted
!        PAUX      : fields to remain gridpoint

!     EXTERNALS.
!     ---------

! Author
! ------
!   ORIGINAL : 88-02-04

! Modifications
! -------------
!   JJMorcrette : 20090217 PP of aerosol and UV processor outputs
!   K. Yessad (Aug 2009): rewrite vertical interpolator in post-processing.
!   Y. Bouteloup  : 14-07-2009 Add radiative cloud water and ice (YIRAD and YLRAD)
!   K. Yessad (Jan 2011): introduce INTDYN_MOD structures.
!   H. Hersbach   : 01-04-2011 auxiliary diagnostic radiation fields
!   K. Yessad (Dec 2011): various contributions.
!   R. El Khatib 22-Aug-2012 cfpfmt=model replaced by .not.lhpos
!      R. El Khatib 17-Jul-2013 FABEC post-processing
!   T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!   K. Yessad (July 2014): Move some variables.
!   A. Geer        Dec 2015  Remove spurious 0:KFLEVG dimension in PP routines
!   R. El Khatib 18-Feb-2015 Revert the change 0:KFLEVG for ZTT0 because PP2DINT
!   is still awaiting for 0:KFLEVG arrays
!   K. Yessad (Dec 2016): Prune obsolete options.
!   K. Yessad (June 2017): Introduce NHQE model.
!   K. Yessad (Feb 2018): remove deep-layer formulations.
!   K. Yessad (Apr 2018): introduce key L_RDRY_VD (ensure consistent definition of "dver" everywhere).
!   O. Jaron (Sept 2018): Correction of deformation terms
! End Modifications
!------------------------------------------------------------------------------

USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE GEOMETRY_MOD           , ONLY : GEOMETRY
USE SURFACE_FIELDS_MIX     , ONLY : TSURF
USE YOMDYN                 , ONLY : TDYN
USE YOMDYNA                , ONLY : TDYNA
USE PARKIND1               , ONLY : JPIM, JPRB
USE YOMHOOK                , ONLY : LHOOK, DR_HOOK, JPHOOK
USE PARFPOS                , ONLY : JPOSFSU, JPOSVX2
USE YOMCT0                 , ONLY : LRPLANE
USE YOMCST                 , ONLY : RG, RATM, RD, RV, RCPD, RKAPPA, RTT
USE YOMSTA                 , ONLY : NLEXTRAP
USE YOMFPC                 , ONLY : TNAMFPSCI
USE TYPE_FPRQDYNS, ONLY : TYPE_FPRQDYN
USE YOMAFN                 , ONLY : ALL_FULLPOS_TYPES
USE TYPE_GFLFLDS           , ONLY : TYPE_IGFLFLD,TYPE_IGFLFLDD,TYPE_LGFLFLDD
USE INTDYN_MOD             , ONLY : YYTHWPP, YYTCTYPP, YYTXYBDERPP, YYTXYBPP
USE YOMPHY                 , ONLY : TPHY
USE YOE_AERODIAG           , ONLY : NPAERAOT, NPAERLISI
USE YOMVERT                , ONLY : TVAB

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE (TYPE_FPRQDYN),  INTENT(IN) :: YDQTYPE
TYPE (TNAMFPSCI),  INTENT(IN) :: YDNAMFPSCI
TYPE(ALL_FULLPOS_TYPES), INTENT(IN) :: YDTFP
LOGICAL           ,INTENT(IN)    :: LDHPOS
TYPE(TVAB)        ,INTENT(IN)    :: YDFPVAB
TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TSURF)       ,INTENT(IN)    :: YDSURF
TYPE(TDYN)        ,INTENT(IN)    :: YDDYN
TYPE(TDYNA)       ,INTENT(IN)    :: YDDYNA
TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(IN) :: YDML_GCONF
TYPE(TPHY)        ,INTENT(IN)    :: YDPHY
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA
INTEGER(KIND=JPIM),INTENT(IN)    :: KST
INTEGER(KIND=JPIM),INTENT(IN)    :: KND
INTEGER(KIND=JPIM),INTENT(IN)    :: KOPLEV
INTEGER(KIND=JPIM),INTENT(IN)    :: KGFL
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDSPP
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDAUX
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEVG
INTEGER(KIND=JPIM),INTENT(IN)    :: KPPM
CHARACTER(LEN=1)  ,INTENT(IN)    :: CDCONF
TYPE(TYPE_IGFLFLDD),INTENT(IN)   :: YDIN
TYPE(TYPE_LGFLFLDD),INTENT(IN)   :: YDGFL
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRCORI(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGM(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRATATH(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: POROG(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: POROGL(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: POROGM(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PXLEV(KOPLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSI(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEV(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PU(KPROMA,KFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUL(KPROMA,KFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDIV(KPROMA,KFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PV(KPROMA,KFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVL(KPROMA,KFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVOR(KPROMA,KFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KPROMA,KFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTL(KPROMA,KFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTM(KPROMA,KFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPD(KPROMA,KFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPDL(KPROMA,KFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPDM(KPROMA,KFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSVD(KPROMA,KFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PNHX(KPROMA,KFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSP(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPL(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPM(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN), TARGET :: PGFL(KPROMA,KFLEVG,KGFL)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEXTRA(KPROMA,YDSURF%YSD_XAD%NLEVS,YDSURF%YSD_XAD%NDIM)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEXRAD(KPROMA,YDSURF%YSD_XRD%NLEVS,YDSURF%YSD_XRD%NDIM)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEXTDI(KPROMA,YDSURF%YSD_DID%NLEVS,YDSURF%YSD_DID%NDIM)
INTEGER(KIND=JPIM),INTENT(IN)    :: KCUFNR
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRMCUFGP(KPROMA,KCUFNR)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PGPP(KPROMA,KFLDSPP)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PAUX(KPROMA,KFLDAUX)
!     ------------------------------------------------------------------
! ==== REAL COMPUTED IN PARTS 1.2.1, 1.2.2, 1.2.3  ===================
REAL(KIND=JPRB) :: ZOROGL(KPROMA),ZOROGM(KPROMA),ZGMAPP(KPROMA)
REAL(KIND=JPRB), TARGET :: ZNULL(KPROMA,KFLEVG)
REAL(KIND=JPRB), POINTER :: ZQT0(:,:), ZLT0(:,:), ZIT0(:,:), ZRT0(:,:), ZST0(:,:), ZGT0(:,:), ZHT0(:,:), ZGFLT0(:,:)
REAL(KIND=JPRB), POINTER :: ZDPHYCTY(:,:) ! mass change rate from physics (s-1)
REAL(KIND=JPRB) :: ZQT0L(KPROMA,KFLEVG),ZQT0M(KPROMA,KFLEVG)
REAL(KIND=JPRB) :: ZGRHLT0(KPROMA,KFLEVG)
REAL(KIND=JPRB) :: ZUT0(KPROMA,KFLEVG),ZVT0(KPROMA,KFLEVG)
REAL(KIND=JPRB) :: ZUT0L(KPROMA,KFLEVG),ZVT0L(KPROMA,KFLEVG)
REAL(KIND=JPRB) :: ZDIVT0(KPROMA,0:KFLEVG),ZVORT0(KPROMA,0:KFLEVG)
REAL(KIND=JPRB) :: ZTT0(KPROMA,0:KFLEVG)
REAL(KIND=JPRB) :: ZTT0L(KPROMA,KFLEVG),ZTT0M(KPROMA,KFLEVG)
REAL(KIND=JPRB) :: ZPDT0(KPROMA,KFLEVG)
REAL(KIND=JPRB) :: ZPDT0L(KPROMA,KFLEVG),ZPDT0M(KPROMA,KFLEVG)
REAL(KIND=JPRB) :: ZSP0L(KPROMA),ZSP0M(KPROMA)
!!!
! ==== REAL COMPUTED IN PART 2.1  ====================================
REAL(KIND=JPRB) :: ZR(KPROMA,KFLEVG),ZCP(KPROMA,KFLEVG),ZKAP(KPROMA,KFLEVG)
REAL(KIND=JPRB) :: ZPRESH(KPROMA,0:KFLEVG),ZPRESF(KPROMA,KFLEVG)
REAL(KIND=JPRB) :: ZXYB(KPROMA,KFLEVG,YYTXYBPP%NDIM)
REAL(KIND=JPRB) :: ZRPRESF(KPROMA,KFLEVG)
REAL(KIND=JPRB) :: ZXYBDER(KPROMA,KFLEVG,YYTXYBDERPP%NDIM)
REAL(KIND=JPRB) :: ZRT(KPROMA,KFLEVG),ZRTL(KPROMA,KFLEVG),ZRTM(KPROMA,KFLEVG)
REAL(KIND=JPRB) :: ZRDT(KPROMA,KFLEVG)
REAL(KIND=JPRB) :: ZPDEP(KPROMA,KFLEVG),ZNH_DELP(KPROMA,KFLEVG)
REAL(KIND=JPRB) :: ZNHPREF(KPROMA,KFLEVG),ZNHPREH(KPROMA,0:KFLEVG)
REAL(KIND=JPRB) :: ZNHPPI(KPROMA,KFLEVG) ,ZRNHPPI(KPROMA,KFLEVG)
REAL(KIND=JPRB) :: ZRRED(KPROMA,KFLEVG)
REAL(KIND=JPRB) :: ZUSL_NHPREL(KPROMA,KFLEVG),ZUSL_NHPREM(KPROMA,KFLEVG)
REAL(KIND=JPRB) :: ZUSL_LNNHPREFL(KPROMA,KFLEVG),ZUSL_LNNHPREFM(KPROMA,KFLEVG)
REAL(KIND=JPRB) :: ZQCHAL(KPROMA,KFLEVG),ZQCHAM(KPROMA,KFLEVG)
REAL(KIND=JPRB) :: ZGEOPH(KPROMA,0:KFLEVG),ZGEOPF(KPROMA,0:KFLEVG)
REAL(KIND=JPRB) :: ZGPFL(KPROMA,KFLEVG),ZGPFM(KPROMA,KFLEVG)
REAL(KIND=JPRB) :: ZGPHL(KPROMA,0:KFLEVG),ZGPHM(KPROMA,0:KFLEVG)
REAL(KIND=JPRB) :: ZETADOT(KPROMA,KFLEVG)
REAL(KIND=JPRB) :: ZCTY(KPROMA,0:KFLEVG,YYTCTYPP%NDIM)
REAL(KIND=JPRB) :: ZUVH(KPROMA,0:KFLEVG,YYTHWPP%NDIM)
REAL(KIND=JPRB) :: ZUS(KPROMA),ZVS(KPROMA)
REAL(KIND=JPRB) :: ZNHX(KPROMA,KFLEVG)
REAL(KIND=JPRB) :: ZDVER(KPROMA,KFLEVG),ZWW(KPROMA,0:KFLEVG)
REAL(KIND=JPRB) :: ZMGDW(KPROMA,KFLEVG),ZGDW(KPROMA,KFLEVG)
REAL(KIND=JPRB) :: ZGW(KPROMA,0:KFLEVG),ZGWF(KPROMA,KFLEVG)
REAL(KIND=JPRB) :: ZR2(KPROMA,KFLEVG)

! ==== REAL COMPUTED IN PART 2.2  ====================================
REAL(KIND=JPRB) :: ZES(KPROMA,KFLEVG),ZRHF(KPROMA,KFLEVG)
REAL(KIND=JPRB) :: ZTOUP(KPROMA,0:KFLEVG),ZTETA(KPROMA,0:KFLEVG)
REAL(KIND=JPRB) :: ZSTDF(KPROMA,0:KFLEVG),ZSHDF(KPROMA,0:KFLEVG)
REAL(KIND=JPRB) :: ZEPTH(KPROMA,0:KFLEVG),ZIBET(KPROMA,0:KFLEVG)
REAL(KIND=JPRB) :: ZSIMR(KPROMA,0:KFLEVG),ZSIMRDB(KPROMA,0:KFLEVG)
REAL(KIND=JPRB) :: ZVIRT(KPROMA,0:KFLEVG)

! ==== REAL COMPUTED IN PART 3    ====================================
REAL(KIND=JPRB) :: ZTSTAR(KPROMA),ZT0(KPROMA)
REAL(KIND=JPRB) :: ZPRES(KPROMA,KOPLEV),ZLNPRES(KPROMA,KOPLEV)
REAL(KIND=JPRB) :: ZGEOP(KPROMA,KOPLEV)
REAL(KIND=JPRB) :: ZPXPD(KPROMA,0:KFLEVG,KPPM),ZPXP(KPROMA,0:KFLEVG,KPPM)
REAL(KIND=JPRB) :: ZSTTF(KPROMA,0:KFLEVG),ZSTZF(KPROMA,0:KFLEVG)

! ==== REAL COMPUTED IN PART 4    ====================================
REAL(KIND=JPRB),ALLOCATABLE :: ZEXTRA(:,:)
REAL(KIND=JPRB),ALLOCATABLE :: ZEXRAD(:,:)
REAL(KIND=JPRB) :: ZEDR(KPROMA,KFLEVG)
REAL(KIND=JPRB) :: ZTPW(KPROMA,0:KFLEVG)

! ==== LOCAL POST-PROCESSED FIELDS  ==================================
REAL(KIND=JPRB),ALLOCATABLE :: ZPPGFL(:,:,:)   ! GFL
REAL(KIND=JPRB) :: ZPPU(KPROMA,KOPLEV)         ! U
REAL(KIND=JPRB) :: ZPPV(KPROMA,KOPLEV)         ! V
REAL(KIND=JPRB) :: ZPPT(KPROMA,KOPLEV)         ! T
REAL(KIND=JPRB) :: ZPPUAP(KPROMA,KOPLEV)       ! upper air prehyd
REAL(KIND=JPRB) :: ZPPSP(KPROMA)               ! prehyds
REAL(KIND=JPRB) :: ZPPLNSP(KPROMA)             ! log(prehyds)
REAL(KIND=JPRB) :: ZPPPDEP(KPROMA,KOPLEV)      ! pressure departure
REAL(KIND=JPRB) :: ZPPWVEL(KPROMA,KOPLEV)      ! w
REAL(KIND=JPRB) :: ZPPMGDW(KPROMA,KOPLEV)      ! -G*Delta w
REAL(KIND=JPRB) :: ZPPPHI(KPROMA,KOPLEV)       ! geopotential height
REAL(KIND=JPRB) :: ZPPMG(KPROMA,KOPLEV)        ! Montgomery potential
REAL(KIND=JPRB) :: ZPPHRL(KPROMA,KOPLEV)       ! relative humidity
REAL(KIND=JPRB) :: ZPPVVE(KPROMA,KOPLEV)       ! omega=Dprehyd/Dt
REAL(KIND=JPRB) :: ZPPEVEL(KPROMA,KOPLEV)      ! etadot
REAL(KIND=JPRB) :: ZPPMSLP(KPROMA)             ! MSL pressure
REAL(KIND=JPRB) :: ZPPCUF(KPROMA,KCUFNR)       ! CUF
REAL(KIND=JPRB) :: ZPPPOTV(KPROMA,KOPLEV)      ! potential vorticity
REAL(KIND=JPRB) :: ZPPTETA(KPROMA,KOPLEV)      ! potential temperature
REAL(KIND=JPRB) :: ZPPEPTH(KPROMA,KOPLEV)      ! equivalent potential temp.
REAL(KIND=JPRB) :: ZPPSTDF(KPROMA,KOPLEV)      ! stretching deformation
REAL(KIND=JPRB) :: ZPPSHDF(KPROMA,KOPLEV)      ! shearing deformation
REAL(KIND=JPRB) :: ZPPTWV(KPROMA)              ! total water vapour
REAL(KIND=JPRB) :: ZPPFOLDIN(KPROMA)           ! tropopause folding indicator
REAL(KIND=JPRB) :: ZPPIBET(KPROMA,KOPLEV)      ! isobaric equivalent temperature
REAL(KIND=JPRB) :: ZPPSIMR(KPROMA,KOPLEV)      ! simulated reflectivities in mm/h
REAL(KIND=JPRB) :: ZPPSIMRDB(KPROMA,KOPLEV)    ! simulated reflectivities in dBZ
REAL(KIND=JPRB) :: ZPPECHOT(KPROMA)            ! Pressure of echotop
REAL(KIND=JPRB) :: ZPPVIRT(KPROMA,KOPLEV)      ! virtual potential temperature
REAL(KIND=JPRB) :: ZPPTHPW(KPROMA,KOPLEV)      ! theta'w
REAL(KIND=JPRB) :: ZPPTPW(KPROMA,KOPLEV)       ! t'w
REAL(KIND=JPRB) :: ZPPWIND(KPROMA,KOPLEV)      ! wind velocity
REAL(KIND=JPRB) :: ZPPDIV(KPROMA,KOPLEV)       ! horizontal divergence
REAL(KIND=JPRB) :: ZPPVOR(KPROMA,KOPLEV)       ! vorticity
REAL(KIND=JPRB) :: ZPPABS(KPROMA,KOPLEV)       ! absolute vorticity
REAL(KIND=JPRB) :: ZPPPTB(KPROMA,KOPLEV)       ! pressure of iso-T
REAL(KIND=JPRB) :: ZPPHTB(KPROMA,KOPLEV)       ! height at iso-T
REAL(KIND=JPRB) :: ZPPPICAO(KPROMA)            ! ICAO tropopause pressure
REAL(KIND=JPRB) :: ZPPTICAO(KPROMA)            ! ICAO tropopause temperature
REAL(KIND=JPRB) :: ZPPUJET(KPROMA)             ! max U at the jet
REAL(KIND=JPRB) :: ZPPVJET(KPROMA)             ! max V at the jet
REAL(KIND=JPRB) :: ZPPPJET(KPROMA)             ! pressure at the max wind level
REAL(KIND=JPRB) :: ZPPHTPW(KPROMA)             ! height at iso-T'w=0C
REAL(KIND=JPRB) :: ZPPHTPW1(KPROMA)            ! height at iso-T'w=1C
REAL(KIND=JPRB) :: ZPPHTPW2(KPROMA)            ! height at iso-T'w=1.5C
REAL(KIND=JPRB) :: ZPPMOCO(KPROMA)             ! moisture convergence
REAL(KIND=JPRB) :: ZPPRHO(KPROMA,KOPLEV)       ! humid density
REAL(KIND=JPRB) :: ZPPQNH(KPROMA)              ! QNH pressure
REAL(KIND=JPRB) :: ZPPMSLNH(KPROMA)            ! mean sea level NH pressure
REAL(KIND=JPRB) :: ZPPSPNH(KPROMA)             ! Surface NH pressure
REAL(KIND=JPRB),ALLOCATABLE :: ZPPFUA(:,:,:)   ! free upper air fields
REAL(KIND=JPRB) :: ZPPEDR(KPROMA,KOPLEV)       ! EDR (3D physical diagnotic)
REAL(KIND=JPRB),ALLOCATABLE :: ZPPFSU(:,:)     ! free surface fields
REAL(KIND=JPRB) :: ZPPIN(KPROMA,KOPLEV),ZPPIN1(KPROMA) ! work arrays
! Z(X)H: (X) post-processed at half-levels of pp.
!        Z(X)H(.,JLEV,0) : half level above full level JLEV
!        Z(X)H(.,JLEV,1) : half level below full level JLEV
REAL(KIND=JPRB) :: ZPPWVELH(KPROMA,KOPLEV,0:1)  ! w
REAL(KIND=JPRB) :: ZPH     (KPROMA,KOPLEV,0:1)  ! prehyd
REAL(KIND=JPRB) :: ZPPR(KPROMA,KOPLEV)
REAL(KIND=JPRB) :: ZPPRT(KPROMA,KOPLEV)

! ==== OTHER REAL VARIABLES ==========================================
REAL(KIND=JPRB) :: ZWORK(KPROMA),ZLNWORK(KPROMA)
REAL(KIND=JPRB) :: ZEPS, Z1SRG, ZTROPO, ZXTEMP, ZRHBNDS(2), ZVALUE
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ==== INTEGER VARIABLES =============================================
INTEGER(KIND=JPIM) :: IHILEV, ILOLEV, ILOC, ILOCU, ILOCV, IPTR, ILPT, I, INIT0
INTEGER(KIND=JPIM) :: JH, JL, JLEV, JLEVP, JROF, JVAR, JVEXTR, JCUFNR
INTEGER(KIND=JPIM) :: IPOL2 ! input to POAERO
INTEGER(KIND=JPIM) :: ISLCT1, ISLCT2 ! input to PP2DINT
INTEGER(KIND=JPIM) :: IPPGFL, INFUA, INFSU ! last dimension for ZPP... arrays
INTEGER(KIND=JPIM),ALLOCATABLE :: IGFLCOD(:), IDIN_GFL(:)
INTEGER(KIND=JPIM) :: ICUFCOD(KCUFNR)
INTEGER(KIND=JPIM) :: ILEVB(KPROMA,KOPLEV,KPPM)
INTEGER(KIND=JPIM) :: ILEVBH(KPROMA,KOPLEV,KPPM,0:1)
INTEGER(KIND=JPIM) :: IPTRLEV(KOPLEV) ! IPTRLEV : absolute indexes of levels
TYPE(TYPE_IGFLFLD) :: YLPPIC

! ==== LOGICAL LL[X] VARIABLES FOR INPUT MEMORY TRANSFER =============
LOGICAL :: LLUVDER, LLDIV0, LLVOR0, LLTDER, LLQDER, LLPVO
LOGICAL,ALLOCATABLE :: LLIN_GFL(:) ! Input GFL

! ==== LOGICAL LL[X] VARIABLES FOR PP ================================
LOGICAL,ALLOCATABLE :: LL_GFL(:) ! GFL
LOGICAL :: LL_PPUV    ! Wind vector (U,V)
LOGICAL :: LL_T       ! temperature
LOGICAL :: LL_P       ! Pressure
LOGICAL :: LL_SP      ! surface pressure
LOGICAL :: LL_LNSP    ! log. surface pressure
LOGICAL :: LL_PD      ! n.-h. pressure departure
LOGICAL :: LL_VW      ! true vertical velocity w
LOGICAL :: LL_WWS     ! surface vertical velocity
LOGICAL :: LL_VD      ! vertical divergence variable -G*(Delta w)
LOGICAL :: LL_Z       ! geopotential height
LOGICAL :: LL_MG      ! Montgomery potential
LOGICAL :: LL_HU      ! relative humidity
LOGICAL :: LL_VV      ! vertical velocity omega
LOGICAL :: LL_ETAD    ! vertical velocity etadot
LOGICAL :: LL_MSL     ! mean sea level pressure
LOGICAL :: LL_CUF(KCUFNR) ! Filtred ln(Ps) for monitoring Coupling Updates Frequency
LOGICAL :: LL_PV      ! Potential Vorticity
LOGICAL :: LL_TH      ! Potential temperature
LOGICAL :: LL_ETH     ! Equivalent pot. temperature
LOGICAL :: LL_STD     ! Stretching Deformation
LOGICAL :: LL_SHD     ! Shearing Deformation
LOGICAL :: LL_TWV     ! Total water vapour
LOGICAL :: LL_FOL     ! tropopause folding ind.
LOGICAL :: LL_IET     ! Isobaric equivalent temperature
LOGICAL :: LL_SRE     ! Simulated reflectivity in mm/h
LOGICAL :: LL_SREDB   ! Simulated reflectivity in dbZ
LOGICAL :: LL_ECHOTOP ! Pressure of top of Simul. reflectlectivity
LOGICAL :: LL_THV     ! virtual t or theta
LOGICAL :: LL_THPW    ! theta'w
LOGICAL :: LL_TPW     ! t'w
LOGICAL :: LL_WND     ! wind velocity
LOGICAL :: LL_DIV     ! Divergence
LOGICAL :: LL_VOR     ! Relative Vorticity
LOGICAL :: LL_ABS     ! Absolute Vorticity
LOGICAL :: LL_PCAO    ! ICAO tropopause pressure
LOGICAL :: LL_TCAO    ! ICAO tropopause temperature
LOGICAL :: LL_UJET    ! ICAO U component of the jet
LOGICAL :: LL_VJET    ! ICAO V component of the jet
LOGICAL :: LL_PJET    ! ICAO pressure at jet level
LOGICAL :: LL_PTB     ! Pressure of iso-T
LOGICAL :: LL_HTB     ! height of iso-T
LOGICAL :: LL_HTPW    ! height of iso-T'w=0 C, from bottom
LOGICAL :: LL_HTPW1   ! height of iso-T'w=1 C, from bottom
LOGICAL :: LL_HTPW2   ! height of iso-T'w=1.5 C, from bottom
LOGICAL :: LL_MOCO    ! Moisture convergence
LOGICAL :: LL_RHO     ! humid density
LOGICAL :: LL_QNH     ! qnh pressure
LOGICAL :: LL_FUA(JPOSVX2) ! Free upper air fields
LOGICAL :: LL_EDR     ! EDR (3D physical diagnostic)
LOGICAL :: LL_FSU(JPOSFSU) ! Free surface fields
LOGICAL :: LL_MSLNH     ! msl  nh.  pressure
LOGICAL :: LL_SPNH     ! surf nh. pressure

! ==== OTHER LOGICAL VARIABLES =======================================
LOGICAL :: LL,LLAERO,LLTOP,LLEXTR,LLDYN(4),LLPD,LLGEO,LLCTY,LLWW,LLPOSVER
LOGICAL :: LLGWF,LLGDWI,LLUNIT,LLVFE,LL_SR,LL_SRDB
LOGICAL :: LLBELO(KPROMA,KOPLEV), LLBLOW(KOPLEV)
LOGICAL :: LLBELS(KPROMA,KOPLEV), LLBLES(KOPLEV)
LOGICAL :: LLBELOH(KPROMA,KOPLEV,0:1), LLBLOWH(KOPLEV,0:1)
LOGICAL :: LLBELSH(KPROMA,KOPLEV,0:1), LLBLESH(KOPLEV,0:1)
LOGICAL :: LL_SREDB_ECHOT

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "ctstar.intfb.h"
#include "fpps.intfb.h"
#include "gnhgrpre.intfb.h"
#include "gnhpre.intfb.h"
#include "gnhpreh.intfb.h"
#include "gpcty.intfb.h"
#include "gpept.intfb.h"
#include "gpgeo.intfb.h"
#include "gpgrgeo.intfb.h"
#include "gpxx.intfb.h"
#include "gphlwi.intfb.h"
#include "gphluv.intfb.h"
#include "gpgrxyb.intfb.h"
#include "gpgw.intfb.h"
#include "gpiet.intfb.h"
#include "gphpre.intfb.h"
#include "gpprs0d.intfb.h"
#include "gppvo.intfb.h"
#include "gprcp.intfb.h"
#include "gprh.intfb.h"
#include "gprt.intfb.h"
#include "gpuvs.intfb.h"






#include "pos_prepgfl.intfb.h"
#include "pp2dint.intfb.h"
#include "ppflev.intfb.h"
#include "ppgeop.intfb.h"
#include "ppinit.intfb.h"
#include "ppleta.intfb.h"
#include "pplteta.intfb.h"
#include "ppltp.intfb.h"
#include "pppmer.intfb.h"
#include "ppq.intfb.h"
#include "ppt.intfb.h"
#include "ppthpw.intfb.h"
#include "ppuv.intfb.h"
#include "ppvvel.intfb.h"
#include "ppsta.intfb.h"
#include "ppltw.intfb.h"
#include "gppechot.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('POS',0,ZHOOK_HANDLE)
ASSOCIATE(YDVAB=>YDGEOMETRY%YRVAB,YDVETA=>YDGEOMETRY%YRVETA, &
 & YDVFE=>YDGEOMETRY%YRVFE,YDSTA=>YDGEOMETRY%YRSTA, YDLAP=>YDGEOMETRY%YRLAP, &
 & YDCSGLEG=>YDGEOMETRY%YRCSGLEG, &
 & YDDIMV=>YDGEOMETRY%YRDIMV,YDCVER=>YDGEOMETRY%YRCVER, &
 & YDCSGEOM=>YDGEOMETRY%YRCSGEOM,YDCSGEOM_NB=>YDGEOMETRY%YRCSGEOM_NB, &
 & YDGSGEOM=>YDGEOMETRY%YRGSGEOM, YDGSGEOM_NB=>YDGEOMETRY%YRGSGEOM_NB, YDSPGEOM=>YDGEOMETRY%YSPGEOM, &
 & YGFL=>YDML_GCONF%YGFL,YDDIMF=>YDML_GCONF%YRDIMF)
ASSOCIATE(LFPLOSP=>YDNAMFPSCI%LFPLOSP, LFPRH100=>YDNAMFPSCI%LFPRH100, FPRHMIN=>YDNAMFPSCI%FPRHMIN, &
 & FPRHMAX=>YDNAMFPSCI%FPRHMAX, RFPVCAP=>YDNAMFPSCI%RFPVCAP, NITERPV=>YDNAMFPSCI%NITERPV, &
 & LFPISOPV=>YDNAMFPSCI%LFPISOPV, LISOT_ABOVEG=>YDNAMFPSCI%LISOT_ABOVEG, &
 & LUVDER=>YDDIMF%LUVDER, &
 & NGHG=>YGFL%NGHG, NAEROUT=>YGFL%NAEROUT, NAERO=>YGFL%NAERO, &
 & NUVP=>YGFL%NUVP, NERA40=>YGFL%NERA40, NNOGW=>YGFL%NNOGW, NCHEM=>YGFL%NCHEM, &
 & NEDRP=>YGFL%NEDRP, NGFL_EXT=>YGFL%NGFL_EXT, NGFL_EZDIAG=>YGFL%NGFL_EZDIAG, &
 & YSD_XRD=>YDSURF%YSD_XRD, YSD_XA=>YDSURF%YSD_XA, YSD_XAD=>YDSURF%YSD_XAD, &
 & YSD_XR=>YDSURF%YSD_XR, YSD_DID=>YDSURF%YSD_DID,YSD_DI=>YDSURF%YSD_DI, &
 & LEDR=>YDPHY%LEDR, YLRCH4=>YGFL%YLRCH4, &
 & NEMIS3D=>YGFL%NEMIS3D, YEMIS3D=>YGFL%YEMIS3D, &
 & LAERAOT=>YGFL%LAERAOT, LAERLISI=>YGFL%LAERLISI)

!     ------------------------------------------------------------------

!*       1.    INITIALIZATIONS.
!              ---------------

!*    1.1  LOGICAL AND OTHER SCALAR INITIALIZATIONS

!     1.1.1  Define the "standard order" and pointers for post-processed GFL

! we compute here pointers YLPPIC%I[X] to locate GFL in LL_GFL and ZPPGFL,
! and the number of 2D fields (IPPGFL) we have to store in LL_GFL and ZPPGFL.
! Assumptions:
! - memory is conditionally allocated for multi-fields GFL and also for LRCH4.
! - for the other GFL variables allocation is always done for the time being,
!   even for not post-processed fields.
! These pointers will be re-used for the non-derivative fields stored in ZGFLT0.

IPPGFL=0
CALL SUPTRPPGFL_POS(1,IPPGFL,YLPPIC%IQ)
CALL SUPTRPPGFL_POS(1,IPPGFL,YLPPIC%IL)
CALL SUPTRPPGFL_POS(1,IPPGFL,YLPPIC%ILRAD)
CALL SUPTRPPGFL_POS(1,IPPGFL,YLPPIC%II)
CALL SUPTRPPGFL_POS(1,IPPGFL,YLPPIC%IIRAD)
CALL SUPTRPPGFL_POS(1,IPPGFL,YLPPIC%IA)
CALL SUPTRPPGFL_POS(1,IPPGFL,YLPPIC%IO3)
CALL SUPTRPPGFL_POS(1,IPPGFL,YLPPIC%ICPF)
CALL SUPTRPPGFL_POS(1,IPPGFL,YLPPIC%ISPF)
CALL SUPTRPPGFL_POS(1,IPPGFL,YLPPIC%IDAL)
CALL SUPTRPPGFL_POS(1,IPPGFL,YLPPIC%IDOM)
CALL SUPTRPPGFL_POS(1,IPPGFL,YLPPIC%IUAL)
CALL SUPTRPPGFL_POS(1,IPPGFL,YLPPIC%IUOM)
CALL SUPTRPPGFL_POS(1,IPPGFL,YLPPIC%ILCONV)
CALL SUPTRPPGFL_POS(1,IPPGFL,YLPPIC%IICONV)
CALL SUPTRPPGFL_POS(1,IPPGFL,YLPPIC%IRCONV)
CALL SUPTRPPGFL_POS(1,IPPGFL,YLPPIC%ISCONV)
CALL SUPTRPPGFL_POS(1,IPPGFL,YLPPIC%IUEN)
CALL SUPTRPPGFL_POS(1,IPPGFL,YLPPIC%IUNEBH)
CALL SUPTRPPGFL_POS(1,IPPGFL,YLPPIC%IS)
CALL SUPTRPPGFL_POS(1,IPPGFL,YLPPIC%IRR)
CALL SUPTRPPGFL_POS(1,IPPGFL,YLPPIC%IG)
CALL SUPTRPPGFL_POS(1,IPPGFL,YLPPIC%IH)
CALL SUPTRPPGFL_POS(1,IPPGFL,YLPPIC%ITKE)
CALL SUPTRPPGFL_POS(1,IPPGFL,YLPPIC%IPHYCTY)
IF(NGHG>0 .AND. YLRCH4%LACTIVE) CALL SUPTRPPGFL_POS(1,IPPGFL,YLPPIC%ILRCH4)
CALL SUPTRPPGFL_POS(NEMIS3D,IPPGFL,YLPPIC%IEMIS3D)
CALL SUPTRPPGFL_POS(NGHG,IPPGFL,YLPPIC%IGHG)
CALL SUPTRPPGFL_POS(NGFL_EZDIAG,IPPGFL,YLPPIC%IEZDIAG)
CALL SUPTRPPGFL_POS(NCHEM,IPPGFL,YLPPIC%ICHEM)
CALL SUPTRPPGFL_POS(NAERO,IPPGFL,YLPPIC%IAERO)
CALL SUPTRPPGFL_POS(NGFL_EXT,IPPGFL,YLPPIC%IEXT)
CALL SUPTRPPGFL_POS(NNOGW,IPPGFL,YLPPIC%INOGW)
CALL SUPTRPPGFL_POS(NEDRP,IPPGFL,YLPPIC%IEDRP)
IF (LAERAOT) CALL SUPTRPPGFL_POS(NPAERAOT,IPPGFL,YLPPIC%IAERAOT)
IF (LAERLISI) CALL SUPTRPPGFL_POS(NPAERLISI,IPPGFL,YLPPIC%IAERLISI)
CALL SUPTRPPGFL_POS(NAEROUT,IPPGFL,YLPPIC%IAEROUT)
CALL SUPTRPPGFL_POS(NUVP,IPPGFL,YLPPIC%IUVP)
CALL SUPTRPPGFL_POS(NERA40,IPPGFL,YLPPIC%IERA40)

!     1.1.2  LL_GFL,LLIN_GFL,IGFLCOD,IDIN_GFL

ALLOCATE(LL_GFL(IPPGFL))
ALLOCATE(IGFLCOD(IPPGFL))
ALLOCATE(IDIN_GFL(IPPGFL))
ALLOCATE(LLIN_GFL(IPPGFL))

CALL POS_PREPGFL(YDQTYPE,YDTFP,YGFL,IPPGFL,YDIN,YDGFL,YLPPIC,&
 & LL_GFL,LLIN_GFL,IGFLCOD,IDIN_GFL)

!     1.1.3  LL_[X] for other variables, and ICUFCOD

! provisional piece of code, but later it would be desirable to replace
! YDTFP%CUF1 to YDTFP%CUF5 by YDTFP%CUF.
IF (KCUFNR>=1) ICUFCOD(1)=YDTFP%CUF1%ICOD
IF (KCUFNR>=2) ICUFCOD(2)=YDTFP%CUF2%ICOD
IF (KCUFNR>=3) ICUFCOD(3)=YDTFP%CUF3%ICOD
IF (KCUFNR>=4) ICUFCOD(4)=YDTFP%CUF4%ICOD
IF (KCUFNR>=5) ICUFCOD(5)=YDTFP%CUF5%ICOD

LL_PPUV = YDQTYPE%LL(YDTFP%U%ICOD) .OR. YDQTYPE%LL(YDTFP%V%ICOD) .OR.&
 & (YDQTYPE%LL(YDTFP%VOR%ICOD).AND.(YDQTYPE%ISKP(YDTFP%VOR%ICOD) == 2)) .OR.&
 & (YDQTYPE%LL(YDTFP%DIV%ICOD).AND.(YDQTYPE%ISKP(YDTFP%DIV%ICOD) == 2)) .OR.&
 & (YDQTYPE%LL(YDTFP%ABS%ICOD).AND.(YDQTYPE%ISKP(YDTFP%ABS%ICOD) == 2)) .OR.&
 & YDQTYPE%LL(YDTFP%KHI%ICOD) .OR. YDQTYPE%LL(YDTFP%PSI%ICOD)
LL_T    = YDQTYPE%LL(YDTFP%T%ICOD)
LL_P    = YDQTYPE%LL(YDTFP%P%ICOD)
LL_SP   = YDQTYPE%LL(YDTFP%SP%ICOD)
LL_LNSP = YDQTYPE%LL(YDTFP%LNSP%ICOD)
LL_PD   = YDQTYPE%LL(YDTFP%PD%ICOD)
LL_VW   = YDQTYPE%LL(YDTFP%VW%ICOD)
LL_WWS  = YDQTYPE%LL(YDTFP%WWS%ICOD)
LL_VD   = YDQTYPE%LL(YDTFP%VD%ICOD)
LL_Z    = YDQTYPE%LL(YDTFP%Z%ICOD)
LL_MG   = YDQTYPE%LL(YDTFP%MG%ICOD)
LL_HU   = YDQTYPE%LL(YDTFP%HU%ICOD)
LL_VV   = YDQTYPE%LL(YDTFP%VV%ICOD)
LL_ETAD = YDQTYPE%LL(YDTFP%ETAD%ICOD)
LL_MSL  = YDQTYPE%LL(YDTFP%MSL%ICOD)
DO JCUFNR=1,KCUFNR
  LL_CUF(JCUFNR)=YDQTYPE%LL(ICUFCOD(JCUFNR))
ENDDO
LL_PV   = YDQTYPE%LL(YDTFP%PV%ICOD)
LL_TH   = YDQTYPE%LL(YDTFP%TH%ICOD)
LL_ETH  = YDQTYPE%LL(YDTFP%ETH%ICOD)
LL_STD  = YDQTYPE%LL(YDTFP%STD%ICOD)
LL_SHD  = YDQTYPE%LL(YDTFP%SHD%ICOD)
LL_TWV  = YDQTYPE%LL(YDTFP%TWV%ICOD)
LL_FOL  = YDQTYPE%LL(YDTFP%FOL%ICOD)
LL_IET  = YDQTYPE%LL(YDTFP%IET%ICOD)
LL_SRE  = YDQTYPE%LL(YDTFP%SRE%ICOD)
LL_SREDB =   YDQTYPE%LL(YDTFP%SREDB%ICOD)
LL_ECHOTOP = YDQTYPE%LL(YDTFP%TOPR%ICOD)
LL_THV  = YDQTYPE%LL(YDTFP%THV%ICOD)
LL_THPW = YDQTYPE%LL(YDTFP%THPW%ICOD)
LL_TPW  = YDQTYPE%LL(YDTFP%TPW%ICOD)
LL_WND  = YDQTYPE%LL(YDTFP%WND%ICOD)
LL_DIV  = YDQTYPE%LL(YDTFP%DIV%ICOD).AND.(YDQTYPE%ISKP(YDTFP%DIV%ICOD) == 1)
LL_VOR  = YDQTYPE%LL(YDTFP%VOR%ICOD).AND.(YDQTYPE%ISKP(YDTFP%VOR%ICOD) == 1) .OR.&
 & YDQTYPE%LL(YDTFP%ABS%ICOD).AND.(YDQTYPE%ISKP(YDTFP%ABS%ICOD) == 1)
LL_ABS  = YDQTYPE%LL(YDTFP%ABS%ICOD).AND.(YDQTYPE%ISKP(YDTFP%ABS%ICOD) == 1)
LL_PCAO = YDQTYPE%LL(YDTFP%PCAO%ICOD)
LL_TCAO = YDQTYPE%LL(YDTFP%TCAO%ICOD)
LL_UJET = YDQTYPE%LL(YDTFP%UJET%ICOD)
LL_VJET = YDQTYPE%LL(YDTFP%VJET%ICOD)
LL_PJET = YDQTYPE%LL(YDTFP%PJET%ICOD)
LL_PTB  = YDQTYPE%LL(YDTFP%PTB%ICOD)
LL_HTB  = YDQTYPE%LL(YDTFP%HTB%ICOD)
LL_HTPW = YDQTYPE%LL(YDTFP%HTPW%ICOD)
LL_HTPW1 = YDQTYPE%LL(YDTFP%HTPW1%ICOD)
LL_MOCO = YDQTYPE%LL(YDTFP%MOCO%ICOD)
LL_HTPW2 = YDQTYPE%LL(YDTFP%HTPW2%ICOD)
LL_RHO = YDQTYPE%LL(YDTFP%RHO%ICOD)
LL_QNH = YDQTYPE%LL(YDTFP%QNH%ICOD)
DO JVAR=1,JPOSVX2
  LL_FUA(JVAR) = YDQTYPE%LL(YDTFP%FUA(JVAR)%ICOD)
ENDDO
DO JVAR=1,JPOSFSU
  LL_FSU(JVAR) = YDQTYPE%LL(YDTFP%FSU(JVAR)%ICOD)
ENDDO
LL_EDR = YDQTYPE%LL(YDTFP%EDR%ICOD)
LL_MSLNH = YDQTYPE%LL(YDTFP%MSLNH%ICOD)
LL_SPNH  = YDQTYPE%LL(YDTFP%SPNH%ICOD)

!     1.1.4  LLDYN

! LLDYN says about the calculations of part 2.1 which must be done.
! * If CDCONF='B' or 'M' or 'V' or 'T' or 'S', LLDYN is set to T only
!   when a subset of the calculations under LLDYN (part 2.1) is required.
! * If CDCONF='H', LLDYN(3:4) is set to T only
!   when a subset of the calculations under LLDYN (part 2.1) is required;
!   LLDYN(1:2) is always set to T (calculation of geopotential height required).
! * In all the other cases: LLDYN=T.

! * LLDYN(1)=T => requires pressure departure calculation.
! * LLDYN(2)=T => requires geopotential height calculation.
! * LLDYN(3)=T => requires call to GPCTY
! * LLDYN(4)=T => requires calculation of "w" and -G*(Delta w)
! The following intermediate variables are defined:
! * LLAERO: call to POAERO is needed (uses geopot. as input)
! * LLPD: calculation of ZPDEP is needed
! * LLGEO: calculation of geopotential height is needed
! * LLCTY: call to GPCTY is needed
! * LLWW: calculation of "w" or -G*(Delta w) is needed
LLDYN(:)=.FALSE.
IF (CDCONF=='H'.OR.CDCONF == 'F') LLDYN(1:2)=.TRUE.
LLAERO=LL_UJET.OR.LL_VJET.OR.LL_PJET.OR.LL_PCAO.OR.LL_TCAO
LLPD=LL_PD.OR.LL_RHO.OR.LL_MSLNH.OR.LL_SPNH
LLGEO=LL_Z.OR.LL_MG.OR.LLAERO.OR.LL_PTB.OR.LL_HTB.OR.&
 & LL_HTPW.OR.LL_MOCO.OR.LL_HTPW1.OR.LL_HTPW2.OR.LL_TPW
LLCTY=LL_VV.OR.LL_ETAD
LLWW=LL_VW.OR.LL_WWS.OR.LL_VD
IF(LLPD.OR.LLGEO.OR.LLCTY.OR.LLWW) LLDYN(1)=.TRUE.
IF(LLGEO.OR.LLCTY.OR.LLWW) LLDYN(2)=.TRUE.
IF(LLCTY.OR.LLWW) LLDYN(3)=.TRUE.
IF(LLWW) LLDYN(4)=.TRUE.

!     1.1.5  Other variables.

! LLPVO: call to GPPVO is required
LLPVO=(CDCONF == 'V' .OR. CDCONF == 'T').OR.&
 & (LL_PV.OR.LL_FOL.OR.LL_ETH.OR.LL_THV.OR.LL_TH)
! LLTDER: horizontal derivatives of T are required
LLTDER=LLDYN(1).OR.LLPVO
! LLUVDER: zonal derivatives of (U,V) are required
LLUVDER=LL_STD.OR.LL_SHD
! LLDIV0: horizontal divergence is required
LLDIV0=LLDYN(3).OR.LL_DIV.OR.LL_STD.OR.LL_SHD.OR.LL_MOCO
! LLVOR0: vorticity is required
LLVOR0=LLPVO.OR.LL_VOR.OR.LL_ABS.OR.LL_STD.OR.LL_SHD
! LLQDER: horizontal derivatives of specific humidity are required
LLQDER=LLDYN(1).OR.LL_MOCO

! ----------------------------------------------------------------------------

ILOLEV=1
IHILEV=KOPLEV
ISLCT1=1
ISLCT2=2

ZRHBNDS(1)=FPRHMIN
ZRHBNDS(2)=FPRHMAX

!     1.1.6  tests.

IF (LLQDER .AND. .NOT.(YDGFL%LLQL.AND.YDGFL%LLQM) .AND. .NOT.YDDYNA%LSPRT) THEN
  CALL ABOR1 ('POS : NEEDS Q DERIVATIVES')
ENDIF

!*    1.2  ARRAYS INITIALIZATIONS

!*      1.2.1   Constants

DO JROF=KST,KND
  ZGMAPP(JROF)=PGM(JROF)*PGM(JROF)
  ZOROGL(JROF)=POROGL(JROF)*PGM(JROF)
  ZOROGM(JROF)=POROGM(JROF)*PGM(JROF)
ENDDO

ZNULL(:,:)=0._JPRB

!* Mass change rate from physics as RHS of the continuity equation.
! Specific mass tendencies (unit: s-1) are computed
! in the last physics call (previous time step in IFS).
! They are stored in GFL YPHYCTY.
CALL INIT_GFL(YGFL%YPHYCTY,PGFL,ZNULL,ZDPHYCTY)
! Hydrometeores:
CALL INIT_GFL(YGFL%YQ,PGFL,ZNULL,ZQT0,P0L=ZQT0L,P0M=ZQT0M)
CALL INIT_GFL(YGFL%YL,PGFL,ZNULL,ZLT0)
CALL INIT_GFL(YGFL%YI,PGFL,ZNULL,ZIT0)
CALL INIT_GFL(YGFL%YR,PGFL,ZNULL,ZRT0)
CALL INIT_GFL(YGFL%YS,PGFL,ZNULL,ZST0)
CALL INIT_GFL(YGFL%YG,PGFL,ZNULL,ZGT0)
CALL INIT_GFL(YGFL%YH,PGFL,ZNULL,ZHT0)
ZGRHLT0(KST:KND,1:KFLEVG)=ZGT0(KST:KND,1:KFLEVG)+ZHT0(KST:KND,1:KFLEVG)

!*      1.2.3   Z[X]T0 for GMV and GMVS

DO JLEV=1,KFLEVG
  ZUT0(KST:KND,JLEV)=PU(KST:KND,JLEV)*PGM(KST:KND)
  ZVT0(KST:KND,JLEV)=PV(KST:KND,JLEV)*PGM(KST:KND)
  IF (LLUVDER) THEN
    ZUT0L(KST:KND,JLEV)=PUL(KST:KND,JLEV)*ZGMAPP(KST:KND)
    ZVT0L(KST:KND,JLEV)=PVL(KST:KND,JLEV)*ZGMAPP(KST:KND)
  ENDIF
  IF (LLDIV0) THEN
    ZDIVT0(KST:KND,JLEV)=PDIV(KST:KND,JLEV)*ZGMAPP(KST:KND)
  ENDIF
  IF (LLVOR0) THEN
    ZVORT0(KST:KND,JLEV)=PVOR(KST:KND,JLEV)*ZGMAPP(KST:KND)
  ENDIF
  IF (LLTDER) THEN
    ZTT0L(KST:KND,JLEV)=PTL(KST:KND,JLEV)*PGM(KST:KND)
    ZTT0M(KST:KND,JLEV)=PTM(KST:KND,JLEV)*PGM(KST:KND)
  ENDIF
  IF (YDDYNA%LNHDYN) THEN
    ZPDT0(KST:KND,JLEV)=PSPD(KST:KND,JLEV)
    ZPDT0L(KST:KND,JLEV)=PSPDL(KST:KND,JLEV)*PGM(KST:KND)
    ZPDT0M(KST:KND,JLEV)=PSPDM(KST:KND,JLEV)*PGM(KST:KND)
  ELSE
    ZPDT0(KST:KND,JLEV)=0.0_JPRB
    ZPDT0L(KST:KND,JLEV)=0.0_JPRB
    ZPDT0M(KST:KND,JLEV)=0.0_JPRB
  ENDIF
ENDDO
DO JROF=KST,KND
  ZPRESH(JROF,KFLEVG)=EXP(PSP(JROF))
  ZSP0L(JROF)=ZPRESH(JROF,KFLEVG)*PSPL(JROF)*PGM(JROF)
  ZSP0M(JROF)=ZPRESH(JROF,KFLEVG)*PSPM(JROF)*PGM(JROF)
ENDDO

!*      1.2.4 INITIALIZATION OF THE LOCAL PP FIELD ARRAYS ZPP[X]

IF(YSD_XAD%NUMFLDS > 0) ALLOCATE(ZEXTRA(KPROMA,KFLEVG))
IF(YSD_XRD%NUMFLDS > 0) ALLOCATE(ZEXRAD(KPROMA,KFLEVG))

INFUA=COUNT(LL_FUA)
INFSU=COUNT(LL_FSU)
ALLOCATE(ZPPFUA(KPROMA,KOPLEV,INFUA))
ALLOCATE(ZPPFSU(KPROMA,INFSU))

ALLOCATE(ZPPGFL(KPROMA,KOPLEV,IPPGFL))
!  INIT0 =-1 : No initialization, for best performance
!  INIT0 = 0 : Initialize to HUGE, for debugging
!  INIT0 > 0 : Initialization to zero, "as before". Discouraged !!
INIT0=-1
IF (INIT0 == 0) THEN
  ZVALUE=HUGE(1._JPRB)
ELSE
  ZVALUE=0._JPRB
ENDIF
IF (INIT0 >= 0) THEN
  ZPPGFL(:,:,:)=ZVALUE
  ZPPU(:,:)   =ZVALUE
  ZPPV(:,:)   =ZVALUE
  ZPPT(:,:)   =ZVALUE
  ZPPUAP(:,:) =ZVALUE
  ZPPSP(:)    =ZVALUE
  ZPPLNSP(:)  =ZVALUE
  ZPPPDEP(:,:)=ZVALUE
  ZPPWVEL(:,:)=ZVALUE
  ZPPMGDW(:,:)=ZVALUE
  ZPPPHI(:,:) =ZVALUE
  ZPPMG(:,:)  =ZVALUE
  ZPPHRL(:,:) =ZVALUE
  ZPPVVE(:,:) =ZVALUE
  ZPPEVEL(:,:)=ZVALUE
  ZPPMSLP(:)  =ZVALUE
  ZPPCUF(:,:) =ZVALUE
  ZPPPOTV(:,:)=ZVALUE
  ZPPTETA(:,:)=ZVALUE
  ZPPEPTH(:,:)=ZVALUE
  ZPPSTDF(:,:)=ZVALUE
  ZPPSHDF(:,:)=ZVALUE
  ZPPFOLDIN(:)=ZVALUE
  ZPPIBET(:,:)=ZVALUE
  ZPPSIMR(:,:)=ZVALUE
  ZPPSIMRDB(:,:)=ZVALUE
  ZPPVIRT(:,:)=ZVALUE
  ZPPTHPW(:,:)=ZVALUE
  ZPPTPW(:,:) =ZVALUE
  ZPPWIND(:,:)=ZVALUE
  ZPPDIV(:,:) =ZVALUE
  ZPPVOR(:,:) =ZVALUE
  ZPPABS(:,:) =ZVALUE
  ZPPPICAO(:) =ZVALUE
  ZPPTICAO(:) =ZVALUE
  ZPPUJET(:)  =ZVALUE
  ZPPVJET(:)  =ZVALUE
  ZPPPJET(:)  =ZVALUE
  ZPPPTB(:,:) =ZVALUE
  ZPPHTB(:,:) =ZVALUE
  ZPPHTPW(:)  =ZVALUE
  ZPPHTPW1(:) =ZVALUE
  ZPPHTPW2(:) =ZVALUE
  ZPPMOCO(:)  =ZVALUE
  ZPPQNH(:)  =ZVALUE
  ZPPRHO(:,:)=ZVALUE
  ZPPR(:,:)  =ZVALUE
  ZPPFSU(:,:)=ZVALUE
  ZPPFUA(:,:,:)=ZVALUE
  ZPPEDR(:,:)=ZVALUE
  ZPPMSLNH(:)=ZVALUE
  ZPPSPNH(:) =ZVALUE
  ZPPECHOT(:)=ZVALUE
ENDIF
!     ------------------------------------------------------------------

!*       2.    DYNAMICS.
!              ---------

!*      2.1.    CALCULATIONS ALSO DONE IN CPG_GP

! Calculations done in this part must remain consistent with CPG_GP ones.

!*      2.1.1   COMPUTATION OF Cp, R AND R/Cp

CALL GPRCP(KPROMA,KST,KND,KFLEVG,PQ=ZQT0,&
 & PQI=ZIT0,PQL=ZLT0,&
 & PQR=ZRT0,PQS=ZST0,PQG=ZGRHLT0,&
 & PCP=ZCP,PR=ZR,PKAP=ZKAP)

IF (LL_RHO) THEN
  CALL GPRCP(KPROMA,KST,KND,KFLEVG,PQ=ZQT0,&
   & PQI=ZIT0,PQL=ZLT0,PR=ZR2)
ENDIF

!*      2.1.2   MODEL LEVEL PRESSURES

! ZPRESH, ZXYB and ZPRESF.
CALL GPHPRE(KPROMA,KFLEVG,KST,KND,YDVAB,YDCVER,ZPRESH,PXYB=ZXYB,PRESF=ZPRESF)

IF (LLDYN(1)) THEN
  ! 1/surface pressure and output of GPGRXYB.
  IF(YDCVER%LVERTFE) THEN
    DO JLEV=1,KFLEVG
      DO JROF=KST,KND
        ZRPRESF(JROF,JLEV)=1.0_JPRB/ZPRESF(JROF,JLEV)
      ENDDO
    ENDDO
  ENDIF
  CALL GPGRXYB(KPROMA,KST,KND,KFLEVG,.FALSE.,YDVAB,YDCVER,ZSP0L,ZSP0M,ZXYB,ZXYBDER)
ENDIF

!*      2.1.3   COMPUTATION OF RT and grad(RT)

IF (LLDYN(1)) THEN
  CALL GPRT(YDDYNA%LSPRT,KPROMA,KST,KND,KFLEVG,RD,RV,&
   & ZR,PT,ZTT0L,ZTT0M,ZQT0L,ZQT0M,ZRT,ZRTL,ZRTM)
  IF (YDDYNA%L_RDRY_VD) THEN
    ZRDT(KST:KND,1:KFLEVG)=RD*PT(KST:KND,1:KFLEVG)
  ELSE
    ZRDT(KST:KND,1:KFLEVG)=ZRT(KST:KND,1:KFLEVG)
  ENDIF
ENDIF

!*      2.1.4   NH PRESSURE DEPARTURE

IF (LLDYN(1)) THEN
  IF (YDDYNA%LNHDYN) THEN
    ! ZPDT0 contains pressure departure prognostic variable.
    CALL GNHPRE(YDGEOMETRY,YDDYNA%NPDVAR,KPROMA,KFLEVG,KST,KND,ZPDT0,ZPRESF,&
     & PNHPREF=ZNHPREF,PNHPPI=ZNHPPI,PRNHPPI=ZRNHPPI,PDEP=ZPDEP)
    ! total nh pressure at half level
    CALL GNHPREH(YDGEOMETRY,YDDYNA%NPDVAR,KPROMA,KFLEVG,KST,KND,ZPDT0,ZPRESH,&
     & ZXYB(1,1,YYTXYBPP%M_DELP),ZXYB(1,1,YYTXYBPP%M_LNPR),&
     & ZNHPPI,ZNHPREF,ZNHPREH,ZNH_DELP)
    DO JLEV=1,KFLEVG
      DO JROF=KST,KND
        ZRRED(JROF,JLEV)=ZR(JROF,JLEV)/ZNHPPI(JROF,JLEV)
      ENDDO
    ENDDO
    CALL GNHGRPRE(YDGEOMETRY,YDDYNA%NPDVAR,KPROMA,KFLEVG,KST,KND,ZXYB(1,1,YYTXYBPP%M_RTGR),ZSP0L,ZSP0M,&
     & ZNHPREF,ZPDT0L,ZPDT0M,&
     & ZUSL_NHPREL,ZUSL_NHPREM,ZUSL_LNNHPREFL,ZUSL_LNNHPREFM,ZQCHAL,ZQCHAM)
  ELSE
    ZPDEP(:,:)=0.0_JPRB
    ZNHPPI(KST:KND,1:KFLEVG)=1.0_JPRB
    ZRNHPPI(KST:KND,1:KFLEVG)=1.0_JPRB
    ZNHPREF(KST:KND,1:KFLEVG)=ZPRESF(KST:KND,1:KFLEVG)
    ZQCHAL(KST:KND,1:KFLEVG)=0.0_JPRB
    ZQCHAM(KST:KND,1:KFLEVG)=0.0_JPRB
    ZRRED(KST:KND,1:KFLEVG)=ZR(KST:KND,1:KFLEVG)
    ZNHPREH(KST:KND,0:KFLEVG)=ZPRESH(KST:KND,0:KFLEVG)
    ZNH_DELP(KST:KND,1:KFLEVG)=ZXYB(KST:KND,1:KFLEVG,YYTXYBPP%M_DELP)
  ENDIF
ENDIF

!*      2.1.5   MODEL LEVEL GEOPOTENTIALS AND THEIR GRADIENTS

IF (LLDYN(2)) THEN
  ! Geopotential height with surface orography=POROG
  DO JROF=KST,KND
    ZGEOPH(JROF,KFLEVG)=POROG(JROF)
  ENDDO
  !! If LSACC=.T. the geop calc. is not strictly accurate, as should be calculated via GPGEO_SACC
  !! but for dignostics purposes, makes sense to output the true hyd. GEOP
  CALL GPGEO(KPROMA,KST,KND,KFLEVG,ZGEOPH,ZGEOPF(1,1),&
      & PT,ZRRED,ZXYB(1,1,YYTXYBPP%M_LNPR),ZXYB(1,1,YYTXYBPP%M_ALPH),YDGEOMETRY%YRVERT_GEOM)
  CALL GPGRGEO(YDGEOMETRY,KPROMA,KST,KND,KFLEVG,&
   & ZRT,ZRTL,ZRTM,&
   & ZXYB(1,1,YYTXYBPP%M_LNPR),ZXYB(1,1,YYTXYBPP%M_ALPH),ZXYBDER,&
   & ZOROGL,ZOROGM,ZGPFL,ZGPFM,ZGPHL,ZGPHM,&
   & LDNHEE=YDDYNA%LNHDYN,PRNHPPI=ZRNHPPI,PQCHAL=ZQCHAL,PQCHAM=ZQCHAM)
ENDIF

!*      2.1.7   INTEGRATION OF HYDROSTATIC AND CALCULATION OF etadot,
!               etadot Dprehyd/Dt and omega/prehyd.

IF (LLDYN(3)) THEN
  CALL GPCTY(YDGEOMETRY%YRVFE,YDGEOMETRY%YRCVER,KPROMA,KST,KND,KFLEVG,.FALSE.,YDVAB,YDVETA,&
   & ZUT0(1,1),ZVT0(1,1),ZDIVT0(1,1),PEV,&
   & ZXYB,ZSP0L,ZSP0M,ZRPRESF,ZCTY,PDPHYCTY=ZDPHYCTY)
  ! etadot is always computed at full levels.
  DO JLEV=1,KFLEVG
    IF(YDCVER%LVERTFE) THEN
      ZETADOT(KST:KND,JLEV)=ZCTY(KST:KND,JLEV,YYTCTYPP%M_EVEL)&
       & *(YDVETA%VETAH(JLEV)-YDVETA%VETAH(JLEV-1))*ZXYB(KST:KND,JLEV,YYTXYBPP%M_RDELP)
    ELSE
      ZETADOT(KST:KND,JLEV)=&
       & 0.5_JPRB*(ZCTY(KST:KND,JLEV,YYTCTYPP%M_EVEL)+ZCTY(KST:KND,JLEV-1,YYTCTYPP%M_EVEL))&
       & *(YDVETA%VETAH(JLEV)-YDVETA%VETAH(JLEV-1))*ZXYB(KST:KND,JLEV,YYTXYBPP%M_RDELP)
    ENDIF
  ENDDO
  !nyc IF ( LSFORC .AND. LSW_FRC ) : call to GPCTY_FORC not yet coded here.
ENDIF

!*      2.1.8   COMPUTES HALF-LEVEL AND SURFACE WINDS, "NHX", "w" AND -G*Dw.

IF (LLDYN(4)) THEN

  ! Half-level and surface winds:
  CALL GPHLWI(YDDIMV,KPROMA,KST,KND,ZXYB(1,1,YYTXYBPP%M_LNPR),ZXYB(1,1,YYTXYBPP%M_ALPH),ZUVH(1,1,YYTHWPP%M_WWI))
  CALL GPHLUV(YDDIMV,KPROMA,KST,KND,ZUT0(1,1),ZVT0(1,1),ZUVH)
  CALL GPUVS(KFLEVG,KPROMA,KST,KND,.FALSE.,ZUT0(1,1),ZVT0(1,1),ZUS,ZVS)

  ! Term NHX:
  IF (YDDYNA%LNHX) THEN
    ZNHX(KST:KND,1:KFLEVG)=PNHX(KST:KND,1:KFLEVG)
  ELSE
    CALL GPXX(YDGEOMETRY,KFLEVG,KPROMA,KST,KND,ZGPHL,ZGPHM,ZGPFL,ZGPFM,ZXYB(1,1,YYTXYBPP%M_LNPR),&
     & ZRT,ZUT0(1,1),ZVT0(1,1),ZUVH(1,0,YYTHWPP%M_UH),ZUVH(1,0,YYTHWPP%M_VH),&
     & ZNHX,PNHPPI=ZNHPPI)
  ENDIF

  ! Vertical divergence:
  IF (YDDYNA%LNHDYN) THEN
    IF (YDDYNA%NVDVAR == 3) THEN
      ! PSVD contains 'dver'
      ZDVER(KST:KND,1:KFLEVG)=PSVD(KST:KND,1:KFLEVG)
    ELSEIF (YDDYNA%NVDVAR == 4) THEN
      ! PSVD contains 'dver + NHX'
      DO JLEV=1,KFLEVG
        DO JROF=KST,KND
          ZDVER(JROF,JLEV)=PSVD(JROF,JLEV)-ZNHX(JROF,JLEV)
        ENDDO
      ENDDO
    ENDIF
  ELSE
    IF (YDDYNA%L_RDRY_VD) THEN
      ! dver_hyd=-(R/Rd)*(D+X+(cv/cp)(omega/prehyd))
      DO JLEV=1,KFLEVG
        DO JROF=KST,KND
          ZDVER(JROF,JLEV)=(-ZR(JROF,JLEV)/RD)*(ZDIVT0(JROF,JLEV)&
           & +ZNHX(JROF,JLEV)+(1._JPRB-ZKAP(JROF,JLEV))*ZCTY(JROF,JLEV,YYTCTYPP%M_VVEL))
        ENDDO
      ENDDO
    ELSE
      ! dver_hyd=-(D+X+(cv/cp)(omega/prehyd))
      DO JLEV=1,KFLEVG
        DO JROF=KST,KND
          ZDVER(JROF,JLEV)=-(ZDIVT0(JROF,JLEV)&
           & +ZNHX(JROF,JLEV)+(1._JPRB-ZKAP(JROF,JLEV))*ZCTY(JROF,JLEV,YYTCTYPP%M_VVEL))
        ENDDO
      ENDDO
    ENDIF
  ENDIF

  ! Vertical velocity "w", and -G*(Delta w):
  Z1SRG=1.0_JPRB/RG
  ! ky: for the time being, LLVFE is set to .F. anytime because half level
  !     values of ZWW may be required in parts 4 and 5, and maybe some other reasons.
  !     Provisionally, a FD calculation of "w" at half levels is kept,
  !     even in a FULL-POS task running with LVERTFE.AND.LVFE_GW
  !! LLVFE=LVERTFE.AND.LVFE_GW
  LLVFE=.FALSE.
  ! Input ZDVER contains "dver"
  LLGWF=.FALSE.
  LLGDWI=.FALSE.
  ! ky: use array ZRDT there to be consistent with definition of "dver"
  CALL GPGW(YDGEOMETRY,YDDYNA%LNHDYN,KFLEVG,KPROMA,KST,KND,LLGWF,LLGDWI,&
   & ZOROGL,ZOROGM,ZXYB(1,1,YYTXYBPP%M_LNPR),ZXYB(1,1,YYTXYBPP%M_ALPH),ZUS,ZVS,ZRDT,ZDVER,ZGW,ZGWF,&
   & LDVFE=LLVFE,PRNHPPI=ZRNHPPI,PGDW=ZGDW)
  IF (LLVFE) THEN
    ZWW(KST:KND,1:KFLEVG)=ZGWF(KST:KND,1:KFLEVG)*Z1SRG
    ZMGDW(KST:KND,1:KFLEVG)=-ZGDW(KST:KND,1:KFLEVG)
  ELSE
    ZWW(KST:KND,0:KFLEVG)=ZGW(KST:KND,0:KFLEVG)*Z1SRG
    DO JLEV=1,KFLEVG
      DO JROF=KST,KND
        ZMGDW(JROF,JLEV)=ZGW(JROF,JLEV-1)-ZGW(JROF,JLEV)
      ENDDO
    ENDDO
  ENDIF

ENDIF

!*      2.2.    ADDITIONAL DYNAMICAL DIAGNOSTICS

! We find there all dynamical quantities which are not required in CPG_GP,
! like potential vorticity, potential temperature, reflectivities.

!*      2.2.2   RELATIVE HUMIDITY

IF (LL_HU) THEN
  CALL GPRH(.FALSE.,KPROMA,KST,KND,KFLEVG,ZRHBNDS(2),ZRHBNDS(1),&
   & ZQT0,PT,ZPRESF,ZES,ZRHF(1,1))
ENDIF

!*      2.2.3   MODEL LEVEL POTENTIAL VORTICITY AND POTENTIAL TEMPERATURE

IF (LLPVO) THEN
  CALL GPPVO(YDGEOMETRY%YRVAB,KPROMA,KST,KND,KFLEVG,&
   & ZPRESF,ZXYB(1,1,YYTXYBPP%M_RDELP),ZKAP,PRCORI,&
   & ZVORT0(1,1),ZUT0(1,1),ZVT0(1,1),PT,&
   & ZTT0M,ZTT0L,ZSP0M,ZSP0L,ZTOUP(1,1),ZTETA(1,1))
ENDIF

!*      2.2.4   MODEL LEVEL DEFORMATIONS

IF ( LLUVDER.AND. .NOT. LRPLANE ) THEN
  DO JLEV=1,KFLEVG
    DO JROF=KST,KND
      ZSTDF(JROF,JLEV)=ZUT0L(JROF,JLEV)-0.5_JPRB*ZDIVT0(JROF,JLEV) &
        & -ZVT0(JROF,JLEV)*PRATATH(JROF)*0.5_JPRB
      ZSHDF(JROF,JLEV)=ZVT0L(JROF,JLEV)-0.5_JPRB*ZVORT0(JROF,JLEV) &
        & +ZUT0(JROF,JLEV)*PRATATH(JROF)*0.5_JPRB
    ENDDO
  ENDDO
ELSEIF (LLUVDER) THEN
  DO JLEV=1,KFLEVG
    DO JROF=KST,KND
      ZSTDF(JROF,JLEV)=ZUT0L(JROF,JLEV)-0.5_JPRB*ZDIVT0(JROF,JLEV)
      ZSHDF(JROF,JLEV)=ZVT0L(JROF,JLEV)-0.5_JPRB*ZVORT0(JROF,JLEV)
    ENDDO
  ENDDO
ENDIF

!*      2.2.5   MODEL LEVEL EQUIVALENT POTENTIAL TEPERATURE

IF (LL_ETH) THEN
  CALL GPEPT(YDPHY,KPROMA,KST,KND,KFLEVG,ZTETA(1,1),PT,&
   & ZPRESF,ZEPTH(1,1))
ENDIF

!*      2.2.6   MODEL LEVEL ISOBARIC EQUIVALENT TEMPERATURE

IF (LL_IET) THEN
  CALL GPIET(YDPHY,KPROMA,KST,KND,KFLEVG,ZQT0,PT,ZCP,&
   & ZPRESF,ZIBET(1,1))
ENDIF

!*      2.2.7    MODEL LEVEL SIMULATED REFLECTIVITIES

IF (LL_SRE.OR.LL_SREDB.OR.LL_ECHOTOP) THEN
  LLUNIT=.TRUE.
  LL_SREDB_ECHOT = LL_SREDB.OR.LL_ECHOTOP
  CALL GPPRS0D(KPROMA,KST,KND,KFLEVG,PT,ZQT0,ZLT0,ZRT0,ZIT0,&
   & ZST0,ZGRHLT0,LLUNIT,LL_SRE,LL_SREDB_ECHOT,ZSIMR(1,1),ZSIMRDB(1,1))
ENDIF

!*      2.2.8    MODEL LEVEL thetav

IF (LL_THV) THEN






  ! OIFS/FC-only - dummy, code will abort
  call abor1("OIFS - ppcvirt should never be called")

ENDIF

! From now on no additional GP.. or GNH.. routine must be called in this routine.

!     ------------------------------------------------------------------
!*       3.    PREPARE FOR INTERPOLATIONS.
!              ---------------------------

!*      3.1      COMPUTE SURFACE TEMPERATURE AND LAPSE RATE.

IF (LL_MOCO.OR.LL_MG.OR.LL_T.OR.LL_THPW.OR.LL_Z.OR.LL_MSL.OR.&
  & LL_QNH.OR.LL_RHO.OR.LL_MSLNH.OR.(CDCONF == 'H').OR.(CDCONF == 'F')) THEN
  CALL CTSTAR(KPROMA,KST,KND,PT(1,NLEXTRAP),&
   & ZPRESH(1,KFLEVG),ZPRESF(1,NLEXTRAP),POROG,ZTSTAR,ZT0)
ENDIF

!*      3.2      COMPUTE POST-PROCESSING LEVELS

IF (CDCONF == 'B') THEN
  DO JLEVP=1,KOPLEV
    DO JROF=KST,KND
      ZPRES(JROF,JLEVP)=PXLEV(JLEVP)
      ZLNPRES(JROF,JLEVP)=LOG(PXLEV(JLEVP))
    ENDDO
  ENDDO
ELSEIF (CDCONF == 'H') THEN
  DO JLEVP=1,KOPLEV
    DO JROF=KST,KND
      ZGEOP(JROF,JLEVP)=RG*PXLEV(JLEVP)+POROG(JROF)
    ENDDO
  ENDDO
  CALL FPPS(KPROMA,KST,KND,KFLEVG,KOPLEV,ZGEOPH,ZGEOP,PT,&
   & ZR,ZPRESH,PTSI,ZPRES)
  DO JLEVP=1,KOPLEV
    DO JROF=KST,KND
      ZLNPRES(JROF,JLEVP)=LOG(ZPRES(JROF,JLEVP))
    ENDDO
  ENDDO
ELSEIF (CDCONF == 'V') THEN
  DO JLEVP = 1 , KOPLEV
     CALL PPLTP (RFPVCAP,NITERPV,LFPISOPV,KPROMA,KST,KND,KFLEVG,ZPRESF,&
      & ZTOUP(1,1),PXLEV(JLEVP),.FALSE.,PRCORI,&
      & ZPRES(1,JLEVP),ZLNPRES(1,JLEVP),&
      & ZWORK)
  ENDDO
ELSEIF (CDCONF == 'T') THEN
  DO JLEVP=1,KOPLEV
    CALL PPLTETA(KPROMA,KST,KND,KFLEVG,ZPRESF,ZTETA(1,1),&
     & PXLEV(JLEVP),ZPRES(1,JLEVP),ZLNPRES(1,JLEVP))
  ENDDO
ELSEIF (CDCONF == 'K') THEN
  DO JLEVP=1,KOPLEV
    ZXTEMP=ABS(PXLEV(JLEVP))
    LLTOP= (PXLEV(JLEVP) < 0.0_JPRB)
    CALL PPLTEMP(YDGEOMETRY%YRSTA,KPROMA,KST,KND,KFLEVG,ZGEOPF(1,1),PT,&
     & ZXTEMP,LLTOP,ZGEOP(1,JLEVP))
  ENDDO
  CALL FPPS(KPROMA,KST,KND,KFLEVG,KOPLEV,ZGEOPH,ZGEOP,PT,&
   & ZR,ZPRESH,PTSI,ZPRES)
  DO JLEVP=1,KOPLEV
    DO JROF=KST,KND
      ZLNPRES(JROF,JLEVP)=LOG(ZPRES(JROF,JLEVP))
    ENDDO
  ENDDO
ELSEIF (CDCONF == 'F') THEN
  DO JLEVP=1,KOPLEV
    DO JROF=KST,KND
      ZGEOP(JROF,JLEVP)=RG*PXLEV(JLEVP)
    ENDDO
  ENDDO
  CALL FPPS(KPROMA,KST,KND,KFLEVG,KOPLEV,ZGEOPH,ZGEOP,PT,&
   & ZR,ZPRESH,PTSI,ZPRES)
  DO JLEVP=1,KOPLEV
    DO JROF=KST,KND
      ZLNPRES(JROF,JLEVP)=LOG(ZPRES(JROF,JLEVP))
    ENDDO
  ENDDO
ELSEIF (CDCONF == 'S') THEN
  IF (UBOUND(YDFPVAB%VBH,DIM=1) == UBOUND(YDVAB%VBH,DIM=1)) THEN
    ZEPS = EPSILON(1.0_JPRB)*10000._JPRB
    LLPOSVER=(MAXVAL(ABS(YDFPVAB%VAH(:)-YDVAB%VAH(:))) < ZEPS) &
      & .AND.(MAXVAL(ABS(YDFPVAB%VBH(:)-YDVAB%VBH(:))) < ZEPS)
  ELSE
    LLPOSVER=.FALSE.
  ENDIF
  IF (LLPOSVER) THEN
    ! No actual vertical interpolations because the two vertical coordinates are equal
    ! cf. PPLETA but without (useless) call to GPHPRE.
    DO JLEVP=1,KOPLEV
      IPTRLEV(JLEVP)=INT(PXLEV(JLEVP),JPIM)
      ZPRES(KST:KND,JLEVP)=ZPRESF(KST:KND,IPTRLEV(JLEVP))
      ZPH(KST:KND,JLEVP,1)=ZPRESH(KST:KND,IPTRLEV(JLEVP))
      ZPH(KST:KND,JLEVP,0)=ZPRESH(KST:KND,IPTRLEV(JLEVP)-1)
    ENDDO
  ELSE
    CALL PPLETA(YDFPVAB,YDCVER,KPROMA,KST,KND,KOPLEV,PXLEV,ZPRESH(1,KFLEVG),ZPRES,ZPH(1,1,1),ZPH(1,1,0))
  ENDIF
  DO JLEVP=1,KOPLEV
    DO JROF=KST,KND
      ZLNPRES(JROF,JLEVP)=LOG(ZPRES(JROF,JLEVP))
    ENDDO
  ENDDO
ENDIF

!*      3.3     SET UP HELP-ARRAYS FOR VERTICAL INTERPOLATION

IF (CDCONF /= 'M') THEN
  CALL PPINIT(KPROMA,KST,KND,KFLEVG,KPPM,ZPRESH,ZPRESF,ZPXP,ZPXPD)
ENDIF

!*      3.4     FIND MODEL LEVEL UNDER SPECIFIED OUTPUT LEVELS

IF ( CDCONF /= 'M') THEN
  CALL PPFLEV(KPROMA,KST,KND,KFLEVG,KOPLEV,KPPM,ZPRES,ZPRESH,&
   & ZPRESF,ILEVB,LLBELO,LLBELS,LLBLOW,LLBLES)
  IF (LL_VD.AND.CDCONF == 'S') THEN
    DO JH=0,1
      CALL PPFLEV(KPROMA,KST,KND,KFLEVG,KOPLEV,KPPM,&
       & ZPH(1,1,JH),ZPRESH,ZPRESF,&
       & ILEVBH(1,1,1,JH),LLBELOH(1,1,JH),&
       & LLBELSH(1,1,JH),LLBLOWH(1,JH),LLBLESH(1,JH))
    ENDDO
  ENDIF
ENDIF

!*      3.5     COMPUTE ICAO ATMOSPHERE

IF (CDCONF /= 'M'.AND.CDCONF /= 'V'.AND.(LL_T.OR.LL_THPW.OR.LL_Z.OR.LL_MG.OR.LL_RHO)) THEN
  CALL PPSTA('PPREF',KPROMA,KST,KND,KFLEVG+1,1,ZPXP(1,0,2),ZPXP(1,0,4),ZSTTF,&
   & ZSTZF)
ENDIF

!     ------------------------------------------------------------------
!*       4.    POST-PROCESSING (INTERPOLATIONS).
!              ---------------------------------

!*    4.1   TRADITIONAL GFL FIELDS

DO JVAR=1,IPPGFL
  IF (LL_GFL(JVAR)) THEN
    IF (LLIN_GFL(JVAR)) THEN
      ZGFLT0(1:,1:) => PGFL(:,:,IDIN_GFL(JVAR))
    ELSE
      ZGFLT0(1:,1:) => ZNULL(:,:)
    ENDIF
    IF (CDCONF /= 'M') THEN
      CALL PPQ(KPROMA,KST,KND,KFLEVG,KOPLEV,ILOLEV,KPPM,&
       & ILEVB,ZPRES,LLBELO,LLBLOW,ZPXP,ZPXPD,ZGFLT0,ZPPGFL(:,:,JVAR))
    ELSE
      ZPPGFL(:,:,JVAR)=ZGFLT0(:,:)
    ENDIF
  ELSE
    ZPPGFL(:,:,JVAR)=0.0_JPRB
  ENDIF
ENDDO

!*    4.2   TRADITIONAL GMV+GMVS FIELDS

!*      4.2.1   WINDS :

IF (LL_PPUV.OR.LL_WND) THEN
  IF (CDCONF /= 'M') THEN
    CALL PPUV(KPROMA,KST,KND,KFLEVG,KOPLEV,ILOLEV,KPPM,ILEVB,LLBELO,LLBLOW,&
     & ZLNPRES,ZPXP,ZPXPD,ZUT0,ZVT0,ZPPU,ZPPV)
  ELSE
    DO JLEVP=1,KFLEVG
      ZPPU(KST:KND,JLEVP)=ZUT0(KST:KND,JLEVP)
      ZPPV(KST:KND,JLEVP)=ZVT0(KST:KND,JLEVP)
    ENDDO
  ENDIF
ENDIF

!*      4.2.2    TEMPERATURE :

IF (LL_T .OR. LL_THPW .OR. LL_RHO .AND. CDCONF /= 'K' ) THEN
  IF (CDCONF == 'V') THEN
    ! Do not use PPT because it contains an implicit definition of
    ! the tropopause
    ZTT0(:,1:KFLEVG)=PT(:,1:KFLEVG)
    CALL PP2DINT(KPROMA,KST,KND,KFLEVG,KPPM,KOPLEV,ILOLEV,&
     & IHILEV,ILEVB,ISLCT2,LLBELO,LLBLOW,ZPRES,ZPXP,ZPXPD,ZTT0,ZPPT)
  ELSEIF (CDCONF /= 'M') THEN
    LL=.TRUE.
    IF (LL_RHO) THEN
      CALL PPT(KPROMA,KST,KND,KFLEVG,KOPLEV,1,KPPM,ILEVB,&
       & ZPRES,ZLNPRES,LLBELO,LLBELS,LLBLOW,LLBLES,LL,&
       & POROG,ZPXP,ZPXPD,ZTSTAR,ZT0,ZR,PT,ZPPT,ZSTTF,&
       & PR2=ZR2,PRTPP=ZPPRT)
    ELSE
      CALL PPT(KPROMA,KST,KND,KFLEVG,KOPLEV,1,KPPM,ILEVB,&
       & ZPRES,ZLNPRES,LLBELO,LLBELS,LLBLOW,LLBLES,LL,&
       & POROG,ZPXP,ZPXPD,ZTSTAR,ZT0,ZR,PT,ZPPT,ZSTTF)
    ENDIF
  ELSE
    DO JLEVP=1,KFLEVG
      ZPPT(KST:KND,JLEVP)=PT(KST:KND,JLEVP)
    ENDDO
    IF (LL_RHO) THEN
      DO JLEVP=1,KFLEVG
        ZPPRT(KST:KND,JLEVP)=ZR2(KST:KND,JLEVP)*PT(KST:KND,JLEVP)
      ENDDO
    ENDIF
  ENDIF
ENDIF

!*      4.2.3    UPPER AIR PRESSURE :

IF (LL_P.AND.CDCONF /= 'B') THEN
  IF (CDCONF /= 'M') THEN
    ZPPUAP(KST:KND,1:KOPLEV)=ZPRES(KST:KND,1:KOPLEV)
  ELSE
    DO JLEVP=1,KFLEVG
      ZPPUAP(KST:KND,JLEVP)=ZPRESF(KST:KND,JLEVP)
    ENDDO
  ENDIF
ENDIF

!*      4.2.4    SURFACE PRESSURE :

IF ( LL_SP.OR.LL_LNSP ) THEN
  ZPPSP(KST:KND)=ZPRESH(KST:KND,KFLEVG)
ENDIF

!*      4.2.5   PRESSURE DEPARTURE :

IF (LL_PD.OR.LL_RHO) THEN
  IF (CDCONF /= 'M') THEN
    CALL PPQ(KPROMA,KST,KND,KFLEVG,KOPLEV,ILOLEV,KPPM,&
     & ILEVB,ZPRES,LLBELO,LLBLOW,ZPXP,ZPXPD,ZPDEP,ZPPPDEP)
  ELSE
    DO JLEVP=1,KFLEVG
      ZPPPDEP(KST:KND,JLEVP)=ZPDEP(KST:KND,JLEVP)
    ENDDO
  ENDIF
ENDIF

!*      4.2.6   TRUE VERTICAL VELOCITY W :

IF (LL_VW) THEN
  IF (CDCONF /= 'M') THEN
    CALL PP2DINT(KPROMA,KST,KND,KFLEVG,KPPM,KOPLEV,&
     & ILOLEV,IHILEV,ILEVB,ISLCT1,LLBELS,LLBLES,&
     & ZPRES,ZPXP,ZPXPD,ZWW,ZPPWVEL)
  ELSE
    ! ZPPWVEL is filled with upper half level values
    DO JLEVP=1,KFLEVG
      ZPPWVEL(KST:KND,JLEVP)=ZWW(KST:KND,JLEVP-1)
    ENDDO
  ENDIF
ENDIF

!*      4.2.7   VERTICAL DIVERGENCE VARIABLE -G*(Delta w):

IF (LL_VD) THEN
  IF (CDCONF == 'S') THEN
    DO JH=0,1
      CALL PP2DINT(KPROMA,KST,KND,KFLEVG,KPPM,KOPLEV,ILOLEV,&
       & IHILEV,ILEVBH(1,1,1,JH),ISLCT1,&
       & LLBELSH(1,1,JH),LLBLESH(1,JH),ZPH(1,1,JH),&
       & ZPXP,ZPXPD,ZWW,ZPPWVELH(1,1,JH))
    ENDDO
    DO JLEVP=1,KOPLEV
      DO JROF=KST,KND
        ZPPMGDW(JROF,JLEVP)=-RG*(ZPPWVELH(JROF,JLEVP,1)&
         & -ZPPWVELH(JROF,JLEVP,0))
      ENDDO
    ENDDO
  ELSEIF (CDCONF /= 'M') THEN
    CALL ABOR1('VERT. DIVERGENCE IS COMPUTABLE ONLY ON ETA LEVELS')
  ELSE
    DO JLEVP=1,KFLEVG
      ZPPMGDW(KST:KND,JLEVP)=ZMGDW(KST:KND,JLEVP)
    ENDDO
  ENDIF
ENDIF

!*      4.2.8   HUMID DENSITY :

IF (LL_RHO) THEN
  IF (CDCONF /= 'M') THEN
    DO JLEVP=1,KOPLEV
      ZPPRHO(KST:KND,JLEVP)=(ZPPPDEP(KST:KND,JLEVP)+ZPRES(KST:KND,JLEVP))/&
        & ZPPRT(KST:KND,JLEVP)
    ENDDO
  ELSEIF (CDCONF == 'V') THEN
    CALL ABOR1('HUMID DENSITY IS NOT COMPUTED ON PV LEVELS')
  ELSE
    DO JLEVP=1,KFLEVG
      ZPPRHO(KST:KND,JLEVP)=(ZPPPDEP(KST:KND,JLEVP)+ZPRES(KST:KND,JLEVP))/&
        & ZPPRT(KST:KND,JLEVP)
    ENDDO
  ENDIF
ENDIF

!*      4.2.9    SURFACE NH PRESSURE :

IF (LL_SPNH) THEN
  ZPPSPNH(KST:KND)=ZNHPREH(KST:KND,KFLEVG)
ENDIF


!*    4.3   OTHER TRADITIONAL FIELDS

!*      4.3.1   GEOPOTENTIAL & MONTGOMERY POTENTIAL

IF (LL_Z .OR. LL_MG) THEN
  IF (CDCONF == 'B'.OR.CDCONF == 'T'.OR.CDCONF == 'S') THEN
    LLEXTR=.TRUE.
    CALL PPGEOP(KPROMA,KST,KND,KFLEVG,KOPLEV,ILOLEV,KPPM,ILEVB,ZPRES,ZLNPRES,&
     & LLBELO,LLBELS,LLBLOW,LLBLES,LLEXTR,POROG,YDGEOMETRY%YRVERT_GEOM,&
     & ZPXP,ZPXPD,ZTSTAR,ZR,ZXYB(1,1,YYTXYBPP%M_LNPR),ZXYB(1,1,YYTXYBPP%M_ALPH),PT,&
     & ZPPPHI,PSTTF=ZSTTF,PRRED=ZRRED)
  ELSEIF (CDCONF == 'V' .OR. CDCONF == 'F') THEN
    ! Do not use PPGEOP because it contains an implicit definition of
    ! the tropopause
    CALL PP2DINT(KPROMA,KST,KND,KFLEVG,KPPM,KOPLEV,ILOLEV,&
     & IHILEV,ILEVB,ISLCT2,LLBELO,LLBLOW,ZPRES,ZPXP,ZPXPD,ZGEOPF,ZPPPHI)
  ELSEIF (CDCONF == 'M') THEN
    DO JLEVP=1,KFLEVG
      ZPPPHI(KST:KND,JLEVP)=ZGEOPF(KST:KND,JLEVP)
    ENDDO
  ELSEIF (CDCONF == 'K') THEN ! ???
    DO JLEVP=1,KOPLEV
      ZPPPHI(KST:KND,JLEVP)=ZGEOP(KST:KND,JLEVP)
    ENDDO
  ENDIF
  IF(LL_MG) THEN
    DO JLEVP=1,KOPLEV
      DO JROF=KST,KND
        ZPPMG(JROF,JLEVP)=ZPPPHI(JROF,JLEVP)+&
         & RCPD*PXLEV(JLEVP)*((ZPRES(JROF,JLEVP)/RATM)**RKAPPA)
      ENDDO
    ENDDO
  ENDIF
ENDIF

!*      4.3.2   RELATIVE HUMIDITY :

IF (LL_HU) THEN
  IF (CDCONF /= 'M') THEN
    CALL PPQ(KPROMA,KST,KND,KFLEVG,KOPLEV,ILOLEV,KPPM,ILEVB,&
     & ZPRES,LLBELO,LLBLOW,ZPXP,ZPXPD,ZRHF,ZPPHRL)
  ELSE
    DO JLEVP=1,KFLEVG
      ZPPHRL(KST:KND,JLEVP)=ZRHF(KST:KND,JLEVP)
    ENDDO
  ENDIF
  IF (LFPRH100) THEN
    ! convert relative humidity field to percent
    DO JLEVP=1,KOPLEV
      DO JROF=KST,KND
        ZPPHRL(JROF,JLEVP)=ZPPHRL(JROF,JLEVP)*100.0_JPRB
      ENDDO
    ENDDO
  ENDIF
ENDIF

!*      4.3.3   VERTICAL VELOCITY "OMEGA":

IF (LL_VV) THEN
  IF (CDCONF /= 'M') THEN
    CALL PPVVEL(YDGEOMETRY%YRVETA,YDGEOMETRY%YRCVER,KPROMA,KST,KND,KFLEVG,KOPLEV,ILOLEV,KPPM,ILEVB,&
     & ZPRES,LLBELO,LLBLOW,ZPXP,ZPXPD,&
     & ZXYB(1,1,YYTXYBPP%M_RDELP),ZCTY(1,0,YYTCTYPP%M_EVEL),ZCTY(1,1,YYTCTYPP%M_VVEL),ZPPVVE,.FALSE.)
  ELSE
    DO JLEVP=1,KFLEVG
      DO JROF=KST,KND
        ZPPVVE(JROF,JLEVP)=ZPXP(JROF,JLEVP,2)*ZCTY(JROF,JLEVP,YYTCTYPP%M_VVEL)
      ENDDO
    ENDDO
  ENDIF
ENDIF

!*      4.3.4   ETADOT

IF (LL_ETAD) THEN
  IF (CDCONF /= 'M') THEN
    CALL PPVVEL(YDGEOMETRY%YRVETA,YDGEOMETRY%YRCVER,KPROMA,KST,KND,KFLEVG,KOPLEV,ILOLEV,KPPM,ILEVB,&
     & ZPRES,LLBELO,LLBLOW,ZPXP,ZPXPD,&
     & ZXYB(1,1,YYTXYBPP%M_RDELP),ZCTY(1,0,YYTCTYPP%M_EVEL),ZCTY(1,1,YYTCTYPP%M_VVEL),ZPPEVEL,.TRUE.)
  ELSE
    DO JLEVP=1,KFLEVG
      ZPPEVEL(KST:KND,JLEVP)=ZETADOT(KST:KND,JLEVP)
    ENDDO
  ENDIF
ENDIF

!*      4.3.5   MSL PRESSURE :

IF (LL_MSL) THEN
  CALL PPPMER(KPROMA,KST,KND,ZPRESH(1,KFLEVG),POROG,ZTSTAR,ZT0,ZPPMSLP)
ENDIF

!*      4.3.6   MSL PRESSURE :

IF (LL_MSLNH) THEN
  CALL PPPMER(KPROMA,KST,KND,ZNHPREH(1,KFLEVG),POROG,ZTSTAR,ZT0,ZPPMSLNH)
ENDIF

!*     4.3.7   Filtred ln(Ps) for monitoring Couplinf Update Frequency :

DO JCUFNR=1,KCUFNR
  IF (LL_CUF(JCUFNR)) ZPPCUF(KST:KND,JCUFNR)=PRMCUFGP(KST:KND,JCUFNR)
ENDDO


!*    4.5   MODERN DYNAMICAL METEOROLOGY FIELDS

!*      4.5.1   POTENTIAL VORTICITY

IF (LL_PV.AND.CDCONF /= 'V') THEN
  IF (CDCONF /= 'M') THEN
    CALL PP2DINT(KPROMA,KST,KND,KFLEVG,KPPM,KOPLEV,ILOLEV,&
     & IHILEV,ILEVB,ISLCT2,LLBELO,LLBLOW,ZPRES,ZPXP,ZPXPD,ZTOUP,ZPPPOTV)
  ELSE
    DO JLEVP=1,KFLEVG
      ZPPPOTV(KST:KND,JLEVP)=ZTOUP(KST:KND,JLEVP)
    ENDDO
  ENDIF
ENDIF

!*      4.5.2   POTENTIAL TEMPERATURE

IF (LL_TH.AND.CDCONF /= 'T') THEN
  IF (CDCONF /= 'M') THEN
    CALL PP2DINT(KPROMA,KST,KND,KFLEVG,KPPM,KOPLEV,ILOLEV,&
     & IHILEV,ILEVB,ISLCT2,LLBELO,LLBLOW,ZPRES,ZPXP,ZPXPD,ZTETA,ZPPTETA)
  ELSE
    DO JLEVP=1,KFLEVG
      ZPPTETA(KST:KND,JLEVP)=ZTETA(KST:KND,JLEVP)
    ENDDO
  ENDIF
ENDIF

!*      4.5.3   EQUIVALENT POTENTIAL TEMPERATURE

IF (LL_ETH) THEN
  IF (CDCONF /= 'M') THEN
    CALL PP2DINT(KPROMA,KST,KND,KFLEVG,KPPM,KOPLEV,ILOLEV,&
     & IHILEV,ILEVB,ISLCT2,LLBELO,LLBLOW,ZPRES,ZPXP,ZPXPD,ZEPTH,ZPPEPTH)
  ELSE
    DO JLEVP=1,KFLEVG
      ZPPEPTH(KST:KND,JLEVP)=ZEPTH(KST:KND,JLEVP)
    ENDDO
  ENDIF
ENDIF

!*      4.5.4   STRETCHING DEFORMATION

IF (LL_STD) THEN
  IF (CDCONF /= 'M') THEN
    CALL PP2DINT(KPROMA,KST,KND,KFLEVG,KPPM,KOPLEV,ILOLEV,&
     & IHILEV,ILEVB,ISLCT2,LLBELO,LLBLOW,ZPRES,ZPXP,ZPXPD,ZSTDF,ZPPSTDF)
  ELSE
    DO JLEVP=1,KFLEVG
      ZPPSTDF(KST:KND,JLEVP)=ZSTDF(KST:KND,JLEVP)
    ENDDO
  ENDIF
ENDIF

!*      4.5.5   SHEARING DEFORMATION

IF (LL_SHD) THEN
  IF (CDCONF /= 'M') THEN
    CALL PP2DINT(KPROMA,KST,KND,KFLEVG,KPPM,KOPLEV,ILOLEV,&
     & IHILEV,ILEVB,ISLCT2,LLBELO,LLBLOW,ZPRES,ZPXP,ZPXPD,ZSHDF,ZPPSHDF)
  ELSE
    DO JLEVP=1,KFLEVG
      ZPPSHDF(KST:KND,JLEVP)=ZSHDF(KST:KND,JLEVP)
    ENDDO
  ENDIF
ENDIF

!*      4.5.6   TOTAL WATER VAPOUR

IF (LL_TWV) THEN
  Z1SRG=1.0_JPRB/RG
  ZPPTWV(:)=0._JPRB
  DO JLEV=1,KFLEVG
    DO JROF=KST,KND
      ZPPTWV(JROF)=ZPPTWV(JROF)+ZQT0(JROF,JLEV)*ZXYB(JROF,JLEV,YYTXYBPP%M_DELP)*Z1SRG
    ENDDO
  ENDDO
ENDIF

!*      4.5.7   TROPOPAUSE FOLDING INDICATOR

IF (LL_FOL) THEN
  ZTROPO=2.E-6_JPRB ! THIS IS ANOTHER INTERNAL DEFINITION OF THE TROPOPAUSE !
  CALL PPLTP(RFPVCAP,NITERPV,LFPISOPV,KPROMA,KST,KND,KFLEVG,ZPRESF,ZTOUP(1,1),ZTROPO,.TRUE.,PRCORI,&
   & ZWORK,ZLNWORK,ZPPFOLDIN)
ENDIF

!*      4.5.8   ISOBARIC EQUIVALENT TEMPERATURE

IF (LL_IET) THEN
  IF (CDCONF /= 'M') THEN
    CALL PP2DINT(KPROMA,KST,KND,KFLEVG,KPPM,KOPLEV,ILOLEV,&
     & IHILEV,ILEVB,ISLCT2,LLBELO,LLBLOW,ZPRES,ZPXP,ZPXPD,ZIBET,ZPPIBET)
  ELSE
    DO JLEVP=1,KFLEVG
      ZPPIBET(KST:KND,JLEVP)=ZIBET(KST:KND,JLEVP)
    ENDDO
  ENDIF
ENDIF

!*      4.5.9   SIMULATED REFLECTIVITIES

! in mm/h
IF (LL_SRE) THEN
  IF (CDCONF /= 'M') THEN
    CALL PP2DINT(KPROMA,KST,KND,KFLEVG,KPPM,KOPLEV,ILOLEV,&
     & IHILEV,ILEVB,ISLCT2,LLBELO,LLBLOW,ZPRES,ZPXP,ZPXPD,ZSIMR,ZPPSIMR)
  ELSE
    DO JLEVP=1,KFLEVG
      ZPPSIMR(KST:KND,JLEVP)=ZSIMR(KST:KND,JLEVP)
    ENDDO
  ENDIF
ENDIF

!in dBZ
IF (LL_SREDB) THEN
  IF (CDCONF /= 'M') THEN
    CALL PP2DINT(KPROMA,KST,KND,KFLEVG,KPPM,KOPLEV,ILOLEV,&
     & IHILEV,ILEVB,ISLCT2,LLBELO,LLBLOW,ZPRES,ZPXP,ZPXPD,ZSIMRDB,ZPPSIMRDB)
  ELSE
    DO JLEVP=1,KOPLEV
      ZPPSIMRDB(KST:KND,JLEVP)=ZSIMRDB(KST:KND,JLEVP)
    ENDDO
  ENDIF
ENDIF

! Echotop
IF (LL_ECHOTOP) THEN
  CALL GPPECHOT(KPROMA,KST,KND,KFLEVG,ZPRESF,ZSIMRDB(1,1),ZPPECHOT)
ENDIF
!*      4.5.10  THETAv

IF (LL_THV) THEN
  IF (CDCONF /= 'M') THEN
    CALL PP2DINT(KPROMA,KST,KND,KFLEVG,KPPM,KOPLEV,ILOLEV,&
     & IHILEV,ILEVB,ISLCT2,LLBELO,LLBLOW,ZPRES,ZPXP,ZPXPD,ZVIRT,ZPPVIRT)
  ELSE
    DO JLEVP=1,KFLEVG
      ZPPVIRT(KST:KND,JLEVP)=ZVIRT(KST:KND,JLEVP)
    ENDDO
  ENDIF
ENDIF

!*    4.6   OTHER FIELDS

!*      4.6.1   MOIST IRREVERSIBLE ADIABATIC POTENTIAL TEMPERATURE

IF (LL_THPW) THEN
  I=YDTFP%THPW%ICOD
  DO JL=1,YDQTYPE%ILEV(I)
    ILPT=YDQTYPE%ILVP(JL,I)
    IF (CDCONF /= 'M') THEN
      CALL PPTHPW(YDPHY,KST,KND,KPROMA,ZPRES(1,ILPT),ZPPT(1,ILPT),ZPPGFL(1,ILPT,YLPPIC%IQ),&
       & ZPPGFL(1,ILPT,YLPPIC%IL),ZPPGFL(1,ILPT,YLPPIC%II),ZPPTHPW(1,ILPT))
    ELSE
      CALL PPTHPW(YDPHY,KST,KND,KPROMA,ZPRESF(1,ILPT),ZPPT(1,ILPT),&
       & ZPPGFL(1,ILPT,YLPPIC%IQ),ZPPGFL(1,ILPT,YLPPIC%IL),&
       & ZPPGFL(1,ILPT,YLPPIC%II),ZPPTHPW(1,ILPT))
    ENDIF
  ENDDO
ENDIF

!*      4.6.2   WIND VELOCITY

IF (LL_WND) THEN
  DO JLEVP=1,KOPLEV
    DO JROF=KST,KND
      ZPPWIND(JROF,JLEVP)=SQRT(&
       & ZPPU(JROF,JLEVP)*ZPPU(JROF,JLEVP)+ZPPV(JROF,JLEVP)*ZPPV(JROF,JLEVP) )
    ENDDO
  ENDDO
ENDIF

!*      4.6.3   DIVERGENCE

IF (LL_DIV) THEN
  IF (CDCONF /= 'M') THEN
    CALL PP2DINT(KPROMA,KST,KND,KFLEVG,KPPM,KOPLEV,ILOLEV,&
     & IHILEV,ILEVB,ISLCT2,LLBELO,LLBLOW,ZPRES,ZPXP,ZPXPD,ZDIVT0,ZPPDIV)
  ELSE
    DO JLEVP=1,KFLEVG
      ZPPDIV(KST:KND,JLEVP)=ZDIVT0(KST:KND,JLEVP)
    ENDDO
  ENDIF
ENDIF

!*      4.6.4   VORTICITY

IF (LL_VOR.OR.LL_ABS) THEN
  IF (CDCONF /= 'M') THEN
    CALL PP2DINT(KPROMA,KST,KND,KFLEVG,KPPM,KOPLEV,ILOLEV,&
     & IHILEV,ILEVB,ISLCT2,LLBELO,LLBLOW,ZPRES,ZPXP,ZPXPD,ZVORT0,ZPPVOR)
  ELSE
    DO JLEVP=1,KFLEVG
      ZPPVOR(KST:KND,JLEVP)=ZVORT0(KST:KND,JLEVP)
    ENDDO
  ENDIF
  IF (LL_ABS) THEN
    IF (.NOT.LDHPOS) THEN
      DO JLEVP=1,KOPLEV
        ZPPABS(KST:KND,JLEVP)=ZPPVOR(KST:KND,JLEVP)+PRCORI(KST:KND)
      ENDDO
    ELSE
      ! Coriolis parameter will be added during horizontal interpolations
      DO JLEVP=1,KOPLEV
        ZPPABS(KST:KND,JLEVP)=ZPPVOR(KST:KND,JLEVP)
      ENDDO
    ENDIF
  ENDIF
ENDIF

!*      4.6.5   ICAO FIELDS

IF (LLAERO) THEN
  IPOL2=4
  CALL POAERO(ZGEOPF(1,1),PT,ZUT0(1,1),ZVT0(1,1),&
   & ZPPPICAO,ZPPTICAO,ZPPUJET,ZPPVJET,ZPPPJET,&
   & ZPRESF,KPROMA,KST,KND,KFLEVG,IPOL2)
ENDIF

!*      4.6.6 ISO-TEMPERATURE FIELDS

IF (LL_PTB .OR. LL_HTB .AND. CDCONF=='K') THEN
  IF (LL_PTB) THEN
    ZPPPTB(KST:KND,1:KOPLEV)=ZPRES(KST:KND,1:KOPLEV)
  ENDIF
  IF (LL_HTB) THEN
    ! Convert from geopotentials to meters
    Z1SRG=1.0_JPRB/RG
    ZPPHTB(KST:KND,1:KOPLEV)=ZGEOP(KST:KND,1:KOPLEV)*Z1SRG
    IF ( LISOT_ABOVEG ) THEN
     DO JLEVP=1,KOPLEV
      DO JROF=KST,KND
       IF ( ZPPHTB(JROF,JLEVP) < POROG(JROF)*Z1SRG ) &
        ZPPHTB(JROF,JLEVP)=-9999._JPRB
      ENDDO
     ENDDO
    ENDIF
  ENDIF
ENDIF

IF ((LL_HTPW).OR.(LL_HTPW1).OR.(LL_HTPW2).OR.(LL_TPW)) THEN
  CALL PPWETPOINT_LOOP(ZPRESF,PT,ZQT0,ZLT0,ZIT0,ZTPW)

  DO JL=KST,KND
      ZTPW(JL,0)=ZTPW(JL,1)
  ENDDO

  IF (LL_TPW) THEN
    IF (CDCONF /= 'M') THEN
      CALL PP2DINT(KPROMA,KST,KND,KFLEVG,KPPM,KOPLEV,ILOLEV,&
      & IHILEV,ILEVB,ISLCT2,LLBELO,LLBLOW,ZPRES,ZPXP,ZPXPD,ZTPW,ZPPTPW)
    ELSE
      DO JLEVP=1,KOPLEV
         ZPPTPW(KST:KND,JLEVP)=ZTPW(KST:KND,JLEVP)
      ENDDO
    ENDIF
  ENDIF

   IF (LL_HTPW) THEN
     ZXTEMP=RTT
     CALL PPLTW(YDGEOMETRY%YRSTA,KPROMA,KST,KND,KFLEVG,ZGEOPF(1,1),ZTPW,ZXTEMP,ZPPHTPW)
     ! Convert from geopotentials to meters
     Z1SRG=1.0_JPRB/RG
     ZPPHTPW(KST:KND)=ZPPHTPW(KST:KND)*Z1SRG
   ENDIF
   IF (LL_HTPW1) THEN
     ZXTEMP=RTT+1.0_JPRB
     CALL PPLTW(YDGEOMETRY%YRSTA,KPROMA,KST,KND,KFLEVG,ZGEOPF(1,1),ZTPW,ZXTEMP,ZPPHTPW1)
     ! Convert from geopotentials to meters
     Z1SRG=1.0_JPRB/RG
     ZPPHTPW1(KST:KND)=ZPPHTPW1(KST:KND)*Z1SRG
   ENDIF
   IF (LL_HTPW2) THEN
     ZXTEMP=RTT+1.5_JPRB
     CALL PPLTW(YDGEOMETRY%YRSTA,KPROMA,KST,KND,KFLEVG,ZGEOPF(1,1),ZTPW,ZXTEMP,ZPPHTPW2)
     ! Convert from geopotentials to meters
     Z1SRG=1.0_JPRB/RG
     ZPPHTPW2(KST:KND)=ZPPHTPW2(KST:KND)*Z1SRG
  ENDIF

ENDIF

!*      4.6.8   QNH PRESSURE :

IF (LL_QNH) THEN
  CALL PPPMER(KPROMA,KST,KND,ZPRESH(1,KFLEVG),POROG,ZTSTAR,ZT0,ZPPQNH,&
   & LDQNH=.TRUE.)
ENDIF

!*      4.6.9  EDR :

IF (LL_EDR) THEN
  IF (YSD_DID%NUMFLDS > 0) THEN
    IF(YSD_DID%NLEVS /= KFLEVG) THEN
      CALL ABOR1('POS: YSD_DID%NLEVS /= KFLEVG not supported')
    ENDIF
  ENDIF
  IF (CDCONF /= 'M') THEN
    IF (LEDR) THEN
      DO JLEV=1,KFLEVG
        ZEDR(KST:KND,JLEV)=PEXTDI(KST:KND,JLEV,YSD_DI%YXEDR%MP)
      ENDDO
    ELSE
      DO JLEV=1,KFLEVG
        ZEDR(KST:KND,JLEV)=0._JPRB
      ENDDO
    ENDIF
    CALL PPQ(KPROMA,KST,KND,KFLEVG,KOPLEV,ILOLEV,KPPM,&
     & ILEVB,ZPRES,LLBELO,LLBLOW,ZPXP,ZPXPD,ZEDR,ZPPEDR)
  ELSE
    IF (LEDR) THEN
      DO JLEVP=1,KFLEVG
        ZPPEDR(KST:KND,JLEVP)=PEXTDI(KST:KND,JLEVP,YSD_DI%YXEDR%MP)
      ENDDO
    ENDIF
  ENDIF
ENDIF

!*    4.7   FREE-USE FIELDS

IF (ANY(LL_FUA)) THEN

  ! provide vertical post-processing here in pressure with result ZPPFUA
  ! in general PPQ is used for interpolation in the vertical in pressure

  IF (YSD_XAD%NUMFLDS > 0) THEN
    IF(YSD_XAD%NLEVS /= KFLEVG) THEN
      CALL ABOR1('POS: YSD_XAD%NLEVS /= KFLEVG not allowed with LFPOS=.T.')
    ENDIF
    IF(YSD_XAD%NUMFLDS > JPOSVX2) THEN
      CALL ABOR1('POS: JPOSVX2 too small for YSD_XAD%NUMFLDS')
    ENDIF
  ENDIF
  IPTR=1
  DO JVAR=1,JPOSVX2
    IF (LL_FUA(JVAR)) THEN
      DO JVEXTR=1,YSD_XAD%NUMFLDS
        IF (YDTFP%FUA(JVAR)%IGRIB == YSD_XA%YXA(JVEXTR)%IGRBCODE(1)) THEN
          DO JLEV=1,YSD_XAD%NLEVS
            ZEXTRA(KST:KND,JLEV)=PEXTRA(KST:KND,JLEV,JVEXTR)
          ENDDO
          CALL PPQ(KPROMA,KST,KND,KFLEVG,KOPLEV,ILOLEV,KPPM,&
           & ILEVB,ZPRES,LLBELO,LLBLOW,ZPXP,ZPXPD,ZEXTRA,ZPPFUA(1,1,IPTR))
          IPTR=IPTR+1
          EXIT
        ENDIF
      ENDDO
    ENDIF
  ENDDO

  IF (YSD_XRD%NUMFLDS > 0) THEN
    IF(YSD_XRD%NLEVS /= KFLEVG) THEN
      CALL ABOR1 ('POS: YSD_XRD%NLEVS /= KFLEVG not allowed with LFPOS=.T.')
    ENDIF
    IF(YSD_XRD%NUMFLDS > JPOSVX2) THEN
      CALL ABOR1 ('POS: JPOSVX2 too small for YSD_XRD%NUMFLDS')
    ENDIF
  ENDIF
  IPTR=1
  DO JVAR=1,JPOSVX2
    IF (LL_FUA(JVAR)) THEN
      DO JVEXTR=1,YSD_XRD%NUMFLDS
        IF (YDTFP%FUA(JVAR)%IGRIB == YSD_XR%YXA(JVEXTR)%IGRBCODE(1)) THEN
          DO JLEV = 1 , YSD_XRD%NLEVS
            ZEXRAD(KST:KND,JLEV) = PEXRAD(KST:KND,JLEV,JVEXTR)
          ENDDO
          CALL PPQ(KPROMA,KST,KND,KFLEVG,KOPLEV,ILOLEV,KPPM,&
           & ILEVB,ZPRES,LLBELO,LLBLOW,ZPXP,ZPXPD,ZEXRAD,ZPPFUA(1,1,IPTR))
          IPTR=IPTR+1
          EXIT
        ENDIF
      ENDDO
    ENDIF
  ENDDO

ENDIF

IF (ANY(LL_FSU)) THEN
ENDIF

!     ------------------------------------------------------------------

!*       5.     STORE POST-PROCESSED FIELDS.
!               ----------------------------

!*    5.1   CONSTANT QUANTITIES (orography, mapping factor)

I=YDTFP%FIS%ICOD
IF (YDQTYPE%LL(I)) CALL FILL_PP2(KPROMA,KST,KND,I,POROG,PGPP,PAUX)

I=YDTFP%GM%ICOD
IF (YDQTYPE%LL(I)) CALL FILL_PP2(KPROMA,KST,KND,I,PGM,PGPP,PAUX)

!*    5.2   TRADITIONAL GFL FIELDS

DO JVAR=1,IPPGFL
  I=IGFLCOD(JVAR)
  IF (YDQTYPE%LL(I)) THEN
    ZPPIN(KST:KND,1:KOPLEV)=ZPPGFL(KST:KND,1:KOPLEV,JVAR)
    CALL FILL_PP3(KPROMA,KST,KND,KOPLEV,I,ZPPIN,PGPP,PAUX)
  ENDIF
ENDDO

!*    5.3   TRADITIONAL GMV+GMVS FIELDS

I=YDTFP%U%ICOD
IF (YDQTYPE%LL(I)) THEN
  DO JLEVP=1,KOPLEV
    ZPPIN(KST:KND,JLEVP)=ZPPU(KST:KND,JLEVP)/PGM(KST:KND)
  ENDDO
  CALL FILL_PP3(KPROMA,KST,KND,KOPLEV,I,ZPPIN,PGPP,PAUX)
ENDIF

I=YDTFP%V%ICOD
IF (YDQTYPE%LL(I)) THEN
  DO JLEVP=1,KOPLEV
    ZPPIN(KST:KND,JLEVP)=ZPPV(KST:KND,JLEVP)/PGM(KST:KND)
  ENDDO
  CALL FILL_PP3(KPROMA,KST,KND,KOPLEV,I,ZPPIN,PGPP,PAUX)
ENDIF

I=YDTFP%T%ICOD
IF (YDQTYPE%LL(I)) CALL FILL_PP3(KPROMA,KST,KND,KOPLEV,I,ZPPT,PGPP,PAUX)

I=YDTFP%P%ICOD
IF (YDQTYPE%LL(I)) CALL FILL_PP3(KPROMA,KST,KND,KOPLEV,I,ZPPUAP,PGPP,PAUX)

IF ((LFPLOSP.AND.YDQTYPE%LL(YDTFP%SP%ICOD)).OR.YDQTYPE%LL(YDTFP%LNSP%ICOD)) THEN
  DO JROF=KST,KND
    ZPPLNSP(JROF)=LOG(ZPPSP(JROF))
  ENDDO
ENDIF
I=YDTFP%SP%ICOD
IF (LFPLOSP) THEN
  IF (YDQTYPE%LL(I)) CALL FILL_PP2(KPROMA,KST,KND,I,ZPPLNSP,PGPP,PAUX)
ELSE
  IF (YDQTYPE%LL(I)) CALL FILL_PP2(KPROMA,KST,KND,I,ZPPSP,PGPP,PAUX)
ENDIF
I=YDTFP%LNSP%ICOD
IF (YDQTYPE%LL(I)) CALL FILL_PP2(KPROMA,KST,KND,I,ZPPLNSP,PGPP,PAUX)

I=YDTFP%SPNH%ICOD
IF (YDQTYPE%LL(I)) CALL FILL_PP2(KPROMA,KST,KND,I,ZPPSPNH,PGPP,PAUX)

I=YDTFP%PD%ICOD
IF (YDQTYPE%LL(I)) CALL FILL_PP3(KPROMA,KST,KND,KOPLEV,I,ZPPPDEP,PGPP,PAUX)

I=YDTFP%VW%ICOD
IF (YDQTYPE%LL(I)) CALL FILL_PP3(KPROMA,KST,KND,KOPLEV,I,ZPPWVEL,PGPP,PAUX)

I=YDTFP%WWS%ICOD
IF (YDQTYPE%LL(I)) CALL FILL_PP2(KPROMA,KST,KND,I,ZWW(1,KFLEVG),PGPP,PAUX)

I=YDTFP%VD%ICOD
IF (YDQTYPE%LL(I)) CALL FILL_PP3(KPROMA,KST,KND,KOPLEV,I,ZPPMGDW,PGPP,PAUX)

!*    5.4   OTHER TRADITIONAL FIELDS

I=YDTFP%Z%ICOD
IF (YDQTYPE%LL(I)) CALL FILL_PP3(KPROMA,KST,KND,KOPLEV,I,ZPPPHI,PGPP,PAUX)

I=YDTFP%MG%ICOD
IF (YDQTYPE%LL(I)) CALL FILL_PP3(KPROMA,KST,KND,KOPLEV,I,ZPPMG,PGPP,PAUX)

I=YDTFP%HU%ICOD
IF (YDQTYPE%LL(I)) CALL FILL_PP3(KPROMA,KST,KND,KOPLEV,I,ZPPHRL,PGPP,PAUX)

I=YDTFP%VV%ICOD
IF (YDQTYPE%LL(I)) CALL FILL_PP3(KPROMA,KST,KND,KOPLEV,I,ZPPVVE,PGPP,PAUX)

I=YDTFP%ETAD%ICOD
IF (YDQTYPE%LL(I)) CALL FILL_PP3(KPROMA,KST,KND,KOPLEV,I,ZPPEVEL,PGPP,PAUX)

I=YDTFP%MSL%ICOD
IF (YDQTYPE%LL(I)) CALL FILL_PP2(KPROMA,KST,KND,I,ZPPMSLP,PGPP,PAUX)

I=YDTFP%MSLNH%ICOD
IF (YDQTYPE%LL(I)) CALL FILL_PP2(KPROMA,KST,KND,I,ZPPMSLNH,PGPP,PAUX)

DO JCUFNR=1,KCUFNR
  I=ICUFCOD(JCUFNR)
  IF (YDQTYPE%LL(I)) CALL FILL_PP2(KPROMA,KST,KND,I,ZPPCUF(1,JCUFNR),PGPP,PAUX)
ENDDO

!*    5.5   MODERN DYNAMICAL METEOROLOGY FIELDS

I=YDTFP%PV%ICOD
IF (YDQTYPE%LL(I)) CALL FILL_PP3(KPROMA,KST,KND,KOPLEV,I,ZPPPOTV,PGPP,PAUX)

I=YDTFP%TH%ICOD
IF (YDQTYPE%LL(I)) CALL FILL_PP3(KPROMA,KST,KND,KOPLEV,I,ZPPTETA,PGPP,PAUX)

I=YDTFP%ETH%ICOD
IF (YDQTYPE%LL(I)) CALL FILL_PP3(KPROMA,KST,KND,KOPLEV,I,ZPPEPTH,PGPP,PAUX)

I=YDTFP%STD%ICOD
IF (YDQTYPE%LL(I)) THEN
  DO JLEVP=1,KOPLEV
    ZPPIN(KST:KND,JLEVP)=ZPPSTDF(KST:KND,JLEVP)/ZGMAPP(KST:KND)
  ENDDO
  CALL FILL_PP3(KPROMA,KST,KND,KOPLEV,I,ZPPIN,PGPP,PAUX)
ENDIF

I=YDTFP%SHD%ICOD
IF (YDQTYPE%LL(I)) THEN
  DO JLEVP=1,KOPLEV
    ZPPIN(KST:KND,JLEVP)=ZPPSHDF(KST:KND,JLEVP)/ZGMAPP(KST:KND)
  ENDDO
  CALL FILL_PP3(KPROMA,KST,KND,KOPLEV,I,ZPPIN,PGPP,PAUX)
ENDIF

I=YDTFP%TPW%ICOD
IF (YDQTYPE%LL(I)) CALL FILL_PP3(KPROMA,KST,KND,KOPLEV,I,ZPPTPW,PGPP,PAUX)

I=YDTFP%TWV%ICOD
IF (YDQTYPE%LL(I)) CALL FILL_PP2(KPROMA,KST,KND,I,ZPPTWV,PGPP,PAUX)

I=YDTFP%FOL%ICOD
IF (YDQTYPE%LL(I)) CALL FILL_PP2(KPROMA,KST,KND,I,ZPPFOLDIN,PGPP,PAUX)

I=YDTFP%IET%ICOD
IF (YDQTYPE%LL(I)) CALL FILL_PP3(KPROMA,KST,KND,KOPLEV,I,ZPPIBET,PGPP,PAUX)

I=YDTFP%SRE%ICOD
IF (YDQTYPE%LL(I)) CALL FILL_PP3(KPROMA,KST,KND,KOPLEV,I,ZPPSIMR,PGPP,PAUX)

I=YDTFP%SREDB%ICOD
IF (YDQTYPE%LL(I)) CALL FILL_PP3(KPROMA,KST,KND,KOPLEV,I,ZPPSIMRDB,PGPP,PAUX)

I=YDTFP%THV%ICOD
IF (YDQTYPE%LL(I)) CALL FILL_PP3(KPROMA,KST,KND,KOPLEV,I,ZPPVIRT,PGPP,PAUX)

!*    5.6   OTHER FIELDS

I=YDTFP%THPW%ICOD
IF (YDQTYPE%LL(I)) CALL FILL_PP3(KPROMA,KST,KND,KOPLEV,I,ZPPTHPW,PGPP,PAUX)

I=YDTFP%WND%ICOD
IF (YDQTYPE%LL(I)) CALL FILL_PP3(KPROMA,KST,KND,KOPLEV,I,ZPPWIND,PGPP,PAUX)

I=YDTFP%DIV%ICOD
IF (YDQTYPE%LL(I)) THEN
  IF (YDQTYPE%ISKP(I) == 2) THEN
    IF (YDQTYPE%ISF(I) >= 1) THEN
      DO JL=1,YDQTYPE%ILEV(I)
        ILPT=YDQTYPE%ILVP(JL,I)
        ILOCU=YDQTYPE%IGT1(I,1)+JL-1
        ILOCV=YDQTYPE%IGT1(I,2)+JL-1
        PGPP(KST:KND,ILOCU)=ZPPU(KST:KND,ILPT)/PGM(KST:KND)
        PGPP(KST:KND,ILOCV)=ZPPV(KST:KND,ILPT)/PGM(KST:KND)
      ENDDO
    ELSE
      CALL ABOR1('POS : INTERNAL ERROR ON DIV')
    ENDIF
  ELSE
    DO JLEVP=1,KOPLEV
      ZPPIN(KST:KND,JLEVP)=ZPPDIV(KST:KND,JLEVP)/ZGMAPP(KST:KND)
    ENDDO
    CALL FILL_PP3(KPROMA,KST,KND,KOPLEV,I,ZPPIN,PGPP,PAUX)
  ENDIF
ENDIF

I=YDTFP%VOR%ICOD
IF (YDQTYPE%LL(I)) THEN
  IF (YDQTYPE%ISKP(I) == 2) THEN
    IF (YDQTYPE%ISF(I) >= 1) THEN
      DO JL=1,YDQTYPE%ILEV(I)
        ILPT=YDQTYPE%ILVP(JL,I)
        ILOCU=YDQTYPE%IGT1(I,1)+JL-1
        ILOCV=YDQTYPE%IGT1(I,2)+JL-1
        PGPP(KST:KND,ILOCU)=ZPPU(KST:KND,ILPT)/PGM(KST:KND)
        PGPP(KST:KND,ILOCV)=ZPPV(KST:KND,ILPT)/PGM(KST:KND)
      ENDDO
    ELSE
      CALL ABOR1('POS : INTERNAL ERROR ON DIV')
    ENDIF
  ELSE
    DO JLEVP=1,KOPLEV
      ZPPIN(KST:KND,JLEVP)=ZPPVOR(KST:KND,JLEVP)/ZGMAPP(KST:KND)
    ENDDO
    CALL FILL_PP3(KPROMA,KST,KND,KOPLEV,I,ZPPIN,PGPP,PAUX)
  ENDIF
ENDIF

I=YDTFP%ABS%ICOD
IF (YDQTYPE%LL(I)) THEN
  IF (YDQTYPE%ISKP(I) == 2) THEN
    IF (YDQTYPE%ISF(I) >= 1) THEN
      DO JL=1,YDQTYPE%ILEV(I)
        ILPT=YDQTYPE%ILVP(JL,I)
        ILOCU=YDQTYPE%IGT1(I,1)+JL-1
        ILOCV=YDQTYPE%IGT1(I,2)+JL-1
        PGPP(KST:KND,ILOCU)=ZPPU(KST:KND,ILPT)/PGM(KST:KND)
        PGPP(KST:KND,ILOCV)=ZPPV(KST:KND,ILPT)/PGM(KST:KND)
      ENDDO
    ELSE
      CALL ABOR1('POS : INTERNAL ERROR ON DIV')
    ENDIF
  ELSE
    DO JLEVP=1,KOPLEV
      ZPPIN(KST:KND,JLEVP)=ZPPABS(KST:KND,JLEVP)/ZGMAPP(KST:KND)
    ENDDO
    CALL FILL_PP3(KPROMA,KST,KND,KOPLEV,I,ZPPIN,PGPP,PAUX)
  ENDIF
ENDIF

I=YDTFP%PSI%ICOD
IF (YDQTYPE%LL(I)) THEN
  DO JL=1,YDQTYPE%ILEV(I)
    ILPT=YDQTYPE%ILVP(JL,I)
    ILOC=YDQTYPE%IGT1(I,1)+JL-1
    PGPP(KST:KND,ILOC)=ZPPV(KST:KND,ILPT)/PGM(KST:KND)
  ENDDO
ENDIF

I=YDTFP%KHI%ICOD
IF (YDQTYPE%LL(I)) THEN
  DO JL=1,YDQTYPE%ILEV(I)
    ILPT=YDQTYPE%ILVP(JL,I)
    ILOC=YDQTYPE%IGT1(I,1)+JL-1
    PGPP(KST:KND,ILOC)=ZPPU(KST:KND,ILPT)/PGM(KST:KND)
  ENDDO
ENDIF

I=YDTFP%PTB%ICOD
IF (YDQTYPE%LL(I)) CALL FILL_PP3(KPROMA,KST,KND,KOPLEV,I,ZPPPTB,PGPP,PAUX)

I=YDTFP%HTB%ICOD
IF (YDQTYPE%LL(I)) CALL FILL_PP3(KPROMA,KST,KND,KOPLEV,I,ZPPHTB,PGPP,PAUX)

I=YDTFP%RHO%ICOD
IF (YDQTYPE%LL(I)) CALL FILL_PP3(KPROMA,KST,KND,KOPLEV,I,ZPPRHO,PGPP,PAUX)

!     Eddy diffusivity rate (3D Phys diagostic)

I=YDTFP%EDR%ICOD
IF (YDQTYPE%LL(I)) CALL FILL_PP3(KPROMA,KST,KND,KOPLEV,I,ZPPEDR,PGPP,PAUX)

!---------------------------------
! surface 2D fields

I=YDTFP%PCAO%ICOD
IF (YDQTYPE%LL(I)) CALL FILL_PP2(KPROMA,KST,KND,I,ZPPPICAO,PGPP,PAUX)

I=YDTFP%TCAO%ICOD
IF (YDQTYPE%LL(I)) CALL FILL_PP2(KPROMA,KST,KND,I,ZPPTICAO,PGPP,PAUX)

I=YDTFP%UJET%ICOD
IF (YDQTYPE%LL(I)) THEN
  ZPPIN1(KST:KND)=ZPPUJET(KST:KND)/PGM(KST:KND)
  CALL FILL_PP2(KPROMA,KST,KND,I,ZPPIN1,PGPP,PAUX)
ENDIF

I=YDTFP%VJET%ICOD
IF (YDQTYPE%LL(I)) THEN
  ZPPIN1(KST:KND)=ZPPVJET(KST:KND)/PGM(KST:KND)
  CALL FILL_PP2(KPROMA,KST,KND,I,ZPPIN1,PGPP,PAUX)
ENDIF

I=YDTFP%PJET%ICOD
IF (YDQTYPE%LL(I)) CALL FILL_PP2(KPROMA,KST,KND,I,ZPPPJET,PGPP,PAUX)

I=YDTFP%HTPW%ICOD
IF (YDQTYPE%LL(I)) CALL FILL_PP2(KPROMA,KST,KND,I,ZPPHTPW,PGPP,PAUX)

I=YDTFP%MOCO%ICOD
IF (YDQTYPE%LL(I)) CALL FILL_PP2(KPROMA,KST,KND,I,ZPPMOCO,PGPP,PAUX)

I=YDTFP%QNH%ICOD
IF (YDQTYPE%LL(I)) CALL FILL_PP2(KPROMA,KST,KND,I,ZPPQNH,PGPP,PAUX)

I=YDTFP%HTPW1%ICOD
IF (YDQTYPE%LL(I)) CALL FILL_PP2(KPROMA,KST,KND,I,ZPPHTPW1,PGPP,PAUX)

I=YDTFP%TOPR%ICOD
IF (YDQTYPE%LL(I)) CALL FILL_PP2(KPROMA,KST,KND,I,ZPPECHOT,PGPP,PAUX)


!*    5.7   FREE-USE FIELDS

IPTR=1
DO JVAR=1,JPOSVX2
  I=YDTFP%FUA(JVAR)%ICOD
  IF (YDQTYPE%LL(I)) THEN
    CALL FILL_PP3(KPROMA,KST,KND,KOPLEV,I,ZPPFUA(1,1,IPTR),PGPP,PAUX)
    IPTR=IPTR+1
  ENDIF
ENDDO

DO JVAR=1,JPOSFSU
  I=YDTFP%FSU(JVAR)%ICOD
  IF (YDQTYPE%LL(I)) CALL FILL_PP2(KPROMA,KST,KND,I,ZPPFSU(1,JVAR),PGPP,PAUX)
ENDDO

!*    5.8   DEALLOCATIONS

DEALLOCATE(ZPPFUA)
DEALLOCATE(ZPPFSU)
IF(YSD_XAD%NUMFLDS > 0) THEN
  DEALLOCATE(ZEXTRA)
ENDIF
IF(YSD_XRD%NUMFLDS > 0) THEN
  DEALLOCATE(ZEXRAD)
ENDIF
IF (ALLOCATED(LLIN_GFL)) DEALLOCATE(LLIN_GFL)
IF (ALLOCATED(LL_GFL)) DEALLOCATE(LL_GFL)
IF (ALLOCATED(IGFLCOD)) DEALLOCATE(IGFLCOD)
IF (ALLOCATED(IDIN_GFL)) DEALLOCATE(IDIN_GFL)
IF (ALLOCATED(ZPPGFL)) DEALLOCATE(ZPPGFL)

!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('POS',1,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

CONTAINS

SUBROUTINE FILL_PP2(KPROMA,KST,KND,KCOD,PIN,PGPP,PAUX)

INTEGER(KIND=JPIM), INTENT(IN)  :: KPROMA
INTEGER(KIND=JPIM), INTENT(IN)  :: KST
INTEGER(KIND=JPIM), INTENT(IN)  :: KND
INTEGER(KIND=JPIM), INTENT(IN)  :: KCOD
REAL(KIND=JPRB)   , INTENT(IN)  :: PIN(KPROMA)
REAL(KIND=JPRB)   , INTENT(OUT) :: PGPP(KPROMA,KFLDSPP)
REAL(KIND=JPRB)   , INTENT(OUT) :: PAUX(KPROMA,KFLDAUX)

INTEGER(KIND=JPIM) :: ILOC
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('POS:FILL_PP2',0,ZHOOK_HANDLE)

IF (YDQTYPE%ISF(KCOD) >= 1) THEN
  ILOC=YDQTYPE%IGT1(KCOD,1)
  PGPP(KST:KND,ILOC) = PIN(KST:KND)
ELSE
  ILOC=YDQTYPE%IGPX(KCOD)
  PAUX(KST:KND,ILOC) = PIN(KST:KND)
ENDIF

IF (LHOOK) CALL DR_HOOK('POS:FILL_PP2',1,ZHOOK_HANDLE)
END SUBROUTINE FILL_PP2

!     ------------------------------------------------------------------

SUBROUTINE FILL_PP3(KPROMA,KST,KND,KOPLEV,KCOD,PIN,PGPP,PAUX)

INTEGER(KIND=JPIM), INTENT(IN)  :: KPROMA
INTEGER(KIND=JPIM), INTENT(IN)  :: KST
INTEGER(KIND=JPIM), INTENT(IN)  :: KND
INTEGER(KIND=JPIM), INTENT(IN)  :: KOPLEV
INTEGER(KIND=JPIM), INTENT(IN)  :: KCOD
REAL(KIND=JPRB)   , INTENT(IN)  :: PIN(KPROMA,KOPLEV)
REAL(KIND=JPRB)   , INTENT(OUT) :: PGPP(KPROMA,KFLDSPP)
REAL(KIND=JPRB)   , INTENT(OUT) :: PAUX(KPROMA,KFLDAUX)

INTEGER(KIND=JPIM) :: ILOC,ILPT,JL
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('POS:FILL_PP3',0,ZHOOK_HANDLE)

DO JL=1,YDQTYPE%ILEV(KCOD)
  ILPT=YDQTYPE%ILVP(JL,KCOD)
  IF (YDQTYPE%ISF(KCOD) >= 1) THEN
    ILOC=YDQTYPE%IGT1(KCOD,1)+JL-1
    PGPP(KST:KND,ILOC) = PIN(KST:KND,ILPT)
  ELSE
    ILOC=YDQTYPE%IGPX(KCOD)+JL-1
    PAUX(KST:KND,ILOC) = PIN(KST:KND,ILPT)
  ENDIF
ENDDO

IF (LHOOK) CALL DR_HOOK('POS:FILL_PP3',1,ZHOOK_HANDLE)
END SUBROUTINE FILL_PP3

!     ------------------------------------------------------------------

SUBROUTINE SUPTRPPGFL_POS(KFL,KPPGFL,KMPP)
! computes pointers to store interpolated GFL variables.

INTEGER(KIND=JPIM), INTENT(IN)     :: KFL ! number of fields
INTEGER(KIND=JPIM), INTENT(INOUT)  :: KPPGFL ! cumulated number of fields
INTEGER(KIND=JPIM), INTENT(OUT)    :: KMPP ! pointer

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('POS:SUPTRPPGFL_POS',0,ZHOOK_HANDLE)

KMPP=KPPGFL+1
KPPGFL=KPPGFL+KFL

IF (LHOOK) CALL DR_HOOK('POS:SUPTRPPGFL_POS',1,ZHOOK_HANDLE)
END SUBROUTINE SUPTRPPGFL_POS

SUBROUTINE INIT_GFL(YDCOMP,PGFL,PDEFAULT,P0,P0L,P0M)

! Initialize an input GFL field according to its attribute
! Using pointers to avoid array copies

USE YOM_YGFL , ONLY : TYPE_GFL_COMP

TYPE(TYPE_GFL_COMP), INTENT(IN) :: YDCOMP
REAL(KIND=JPRB)    , INTENT(IN), TARGET :: PGFL(KPROMA,KFLEVG,KGFL)
REAL(KIND=JPRB)    , INTENT(IN), TARGET :: PDEFAULT(KPROMA,KFLEVG)
REAL(KIND=JPRB)    , INTENT(OUT), POINTER :: P0(:,:)
REAL(KIND=JPRB)    , INTENT(OUT), OPTIONAL :: P0L(KPROMA,KFLEVG)
REAL(KIND=JPRB)    , INTENT(OUT), OPTIONAL :: P0M(KPROMA,KFLEVG)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('POS:INIT_GFL',0,ZHOOK_HANDLE)

IF (YDCOMP%LACTIVE) THEN
  P0(1:,1:) => PGFL(:,:,YDCOMP%MP)
  IF (PRESENT(P0L).AND.PRESENT(P0M)) THEN
    IF (YDCOMP%LCDERS) THEN
      DO JLEV=1,KFLEVG
        DO JROF=KST,KND
          P0L(JROF,JLEV)=PGFL(JROF,JLEV,YDCOMP%MPL)*PGM(JROF)
          P0M(JROF,JLEV)=PGFL(JROF,JLEV,YDCOMP%MPM)*PGM(JROF)
        ENDDO
      ENDDO
    ELSE
      P0L(:,:)=PDEFAULT(:,:)
      P0M(:,:)=PDEFAULT(:,:)
    ENDIF
  ELSEIF(PRESENT(P0L).NEQV.PRESENT(P0M)) THEN
    ! Both derivatives should be present or absent
    CALL ABOR1('POS:INIT_GFL:BOTH DERIVATIVES SHOULD BE PRESENT OR ABSENT !')
  ENDIF
ELSE
  P0(1:,1:) => PDEFAULT(:,:)
ENDIF

IF (LHOOK) CALL DR_HOOK('POS:INIT_GFL',1,ZHOOK_HANDLE)
END SUBROUTINE INIT_GFL

SUBROUTINE PPWETPOINT_LOOP(PAPRS,PT,PQV,PQL,PQI,PWETPOINT)
! A wrapper to pass pointer arrays
REAL(KIND=JPRB) ,INTENT(IN)  :: PAPRS(KPROMA,KFLEVG)
REAL(KIND=JPRB) ,INTENT(IN)  :: PT(KPROMA,KFLEVG)
REAL(KIND=JPRB) ,INTENT(IN)  :: PQV(KPROMA,KFLEVG)
REAL(KIND=JPRB) ,INTENT(IN)  :: PQL(KPROMA,KFLEVG)
REAL(KIND=JPRB) ,INTENT(IN)  :: PQI(KPROMA,KFLEVG)
REAL(KIND=JPRB) ,INTENT(OUT) :: PWETPOINT(KPROMA,0:KFLEVG)
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('POS:PPWETPOINT_LOOP',0,ZHOOK_HANDLE)
DO JLEV=1,KFLEVG
  CALL PPWETPOINT(YDPHY,KST,KND,KPROMA,PAPRS(1,JLEV),PT(1,JLEV),PQV(1,JLEV),PQL(1,JLEV),PQI(1,JLEV),PWETPOINT(1,JLEV))
ENDDO
IF (LHOOK) CALL DR_HOOK('POS:PPWETPOINT_LOOP',1,ZHOOK_HANDLE)
END SUBROUTINE PPWETPOINT_LOOP
!     ------------------------------------------------------------------
END SUBROUTINE POS
