! (C) Copyright 1989- Meteo-France.

SUBROUTINE ENDPOS(YDQTYPE,YDNAMFPSCI,YDAFN,KLMOD,YDFPVAB,YDGEOMETRY,YDML_GCONF,YDDYNA,YDML_PHY_MF,LDZ0H,KPROMA,KST,KND,KAUX,KOPLEV,CDCONF, &
 & PXLEV,KGT1,YDIN_GFL,PABUF,KUT0,KVT0,KTT0,KDEP0,PWVEL0,KOROG,KFPCUFNR,KCUFCOD,KCUFT2,KEDR,PSPT0,PTS0,PGAUX,KITS0,KSD0,KTXI,KTNI, &
 & KSW0,KSF0,KHV0,KSR0,KTCLS,KRHCLS,KTR0,KHUXI,KHUNI,PCFIS,PCLSM,PCVEG,PCIVE,PCARG,PCD2,PCSAB,PGM,PGP,PFP,PMASK)

!**** *ENDPOS* - End of post-processing

!     PURPOSE.
!     --------
!       Post-processing calculations after horizontal interpolations of
!        model fields ; used to compute fields on surface-dependent levels
!        (height, eta levels or PBL physico-dynamic fields)

!**   INTERFACE.
!     ----------
!        *CALL* *ENDPOS(...)

!        EXPLICIT ARGUMENTS :
!        --------------------
!        * INPUT:
!          KLMOD : model/post-processing dynamic primitive variable indicator for each fullpos field :
!                  2 : upper air primitive dynamic variable
!                  1 : surface primitive dynamic variable
!                  0 : not a primitive variable
!        KPROMA    : horizontal dimension
!        KST       : starting index of work
!        KND       : depth of work
!        KAUX    : Number of fields in output array for field remaining gridpoint
!        KOPLEV    : number of post-processing levels
!        CDCONF    : Kind of vertical interpolations :
!                    'H' = height levels
!                    'S' = eta levels
!                    'K' = isothermic levels
!                    'F' = flight levels
!        PXLEV     : array containing the pp levels (heights or eta indexes)
!        YDIN_GFL  : pointers to locate individual GFL in PABUF.
!        KFPCUFNR  : Number of Coupling Updates Frequency fields
!        KCUFCOD   : Codes of Coupling Updates Frequency fields
!        KCUFT2    : Pointers of Coupling Updates Frequency fields
!        KEDR      : pointer to EDR diagnostics in PABUF (3D)
!   --  horizontally interpolated model variables ---
!        PABUF     : buffer containing GMV, GMVS and GFL variables
!        KUT0      : zonal wind
!        KVT0      : meridional wind
!        KTT0      : temperature
!        KDEP0     : n.-h. pressure departure
!        PWVEL0    : true vertical velocity w at half levels
!        KOROG     : surface geopotential
!        PSPT0     : surface pressure variable
!        PTS0      : output surface temperature
!        PGAUX     : Auxilary variables
!        KITS0     : interpolated surface temperature
!        KSD0      : output snow depth
!        KTXI      : increment (Tx-Tcls)
!        KTNI      : increment (Tn-Tcls)
!        KSW0      : output surface water content
!        KSF0      : output surface frost
!        KHV0      : output resistance to evapotranspiration
!        KSR0      : output model dynamical roughness length times g
!        KTCLS     : interpolated cls temperature
!        KRHCLS    : interpolated cls relative moisture
!        KTR0      : output model thermal roughness length times g
!        KHUXI     : increment (HUx-HUcls)
!        KHUNI     : increment (HUn-HUcls)
!        PCFIS     : output orography
!        PCLSM     : output land-sea mask
!        PCVEG     : climatology percentage of vegetation
!        PCIVE     : climatology index of vegetation
!        PCARG     : climatology percentage of silt
!        PCD2      : climatology soild depth
!        PCSAB     : climatology percentage of sand
!        PGM       : map factor
!        PMASK     : mask extraneous output domain part
!        * OUTPUT:
!        PGP       : output fields array for gridpoint fields
!        PFP       : output fields array for fitted fields

!        IMPLICIT ARGUMENTS
!        --------------------

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     AUTHOR.
!     -------
!      RYAD EL KHATIB *METEO-FRANCE*
!      ORIGINAL : 94-04-08

!     MODIFICATIONS.
!     --------------
!      N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!      R. El Khatib : 05-May-2008 Bugfixes on recent developments
!      F. Vana      : 27-Aug-2008 vectorization support for NEC
!      K. Yessad (Aug 2009): rewrite vertical interpolator in post-processing.
!      Y. Bouteloup  : 14-07-2009 Add radiative cloud water and ice (YIRAD and YLRAD)
!      F. Bouyssel  : 15-Jul-2009 LFPCAPEX
!      K. Yessad (Dec 2011): use GPHPRE.
!      R. El Khatib 13-Dec-2012 Fullpos buffers reshaping
!      R. El Khatib 15-Mar-2013 Enable CUF post-processing in NFPOS=2
!      R. El Khatib 17-Jul-2013 FABEC post-processing
!      T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!      F. Voitus    28-Jul-2014 Mask for too wide Domain
!      R. El Khatib 01-Sep-2014 Generalize option LFPCLSTOGMV
!      A. Geer      31-Dec-2015 Remove spurious 0:KFLEVG dimension from PP routines
!      R. El Khatib 28-Jul-2016 Recode LWIDER_DOM
!      Y. Bouteloup: 11-Jan-2017 Convective variables
!      K. Yessad (Feb 2018): remove deep-layer formulations.
!      09-2018 R. Brozkova, A. Bucanek: MOCON diagnostics in offline fullpos
!      09-2018 R. Brozkova: Dataflow for convective temperature
!      10-2019 Y. Seity : bf for thetav in H levels
!      I. Polichtchouk 2021: added geop. height calc for shallow-atm complete-Coriolis eqns.
!     ------------------------------------------------------------------

USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE MODEL_PHYSICS_MF_MOD   , ONLY : MODEL_PHYSICS_MF_TYPE
USE GEOMETRY_MOD           , ONLY : GEOMETRY
USE YOMVERT  , ONLY : TVAB
USE PARKIND1               , ONLY : JPIM     ,JPRB
USE YOMHOOK                , ONLY : LHOOK    ,DR_HOOK, JPHOOK
USE PARFPOS                , ONLY : JPOSVX2, JPOSDYN
USE YOMFPC                 , ONLY : TNAMFPSCI
USE YOMDYNA                , ONLY : TDYNA
USE YOMCST                 , ONLY : RV       ,RCPV     ,RETV     ,&
 &                                  RCW      ,RCS      ,RLVTT    ,RLSTT    ,RTT      ,&
 &                                  RALPW    ,RBETW    ,RGAMW    ,RALPS    ,RBETS    ,&
 &                                  RGAMS    ,RALPD    ,RBETD    ,RGAMD    ,RG
USE YOMSTA                 , ONLY : NLEXTRAP,RDTDZ1
USE TYPE_FPRQDYNS, ONLY : TYPE_FPRQDYN
USE YOMAFN, ONLY : TAFN
USE TYPE_GFLFLDS           , ONLY : TYPE_IGFLFLD
USE INTDYN_MOD             , ONLY : YYTXYB
USE PARDIM                 , ONLY : JPNPPM

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE (TYPE_FPRQDYN),  INTENT(IN) :: YDQTYPE
TYPE (TNAMFPSCI),  INTENT(IN) :: YDNAMFPSCI
TYPE(TAFN)        ,INTENT(IN)    :: YDAFN
INTEGER(KIND=JPIM),INTENT(IN)    :: KLMOD(JPOSDYN)
TYPE(TVAB)        ,INTENT(IN)    :: YDFPVAB
TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(IN) :: YDML_GCONF
TYPE(TDYNA)       ,INTENT(IN)    :: YDDYNA
TYPE(MODEL_PHYSICS_MF_TYPE)  ,INTENT(IN) :: YDML_PHY_MF
LOGICAL           ,INTENT(IN)    :: LDZ0H
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA
INTEGER(KIND=JPIM),INTENT(IN)    :: KST
INTEGER(KIND=JPIM),INTENT(IN)    :: KND
INTEGER(KIND=JPIM),INTENT(IN)    :: KAUX
INTEGER(KIND=JPIM),INTENT(IN)    :: KOPLEV
CHARACTER(LEN=1)  ,INTENT(IN)    :: CDCONF
REAL(KIND=JPRB)   ,INTENT(IN)    :: PXLEV(KOPLEV)
INTEGER(KIND=JPIM),INTENT(IN)    :: KGT1
TYPE(TYPE_IGFLFLD),INTENT(IN)    :: YDIN_GFL
REAL(KIND=JPRB)   ,INTENT(IN)    :: PABUF(:,:)
INTEGER(KIND=JPIM),INTENT(IN)    :: KUT0
INTEGER(KIND=JPIM),INTENT(IN)    :: KVT0
INTEGER(KIND=JPIM),INTENT(IN)    :: KTT0
INTEGER(KIND=JPIM),INTENT(IN)    :: KDEP0
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWVEL0(KPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)
INTEGER(KIND=JPIM),INTENT(IN)    :: KOROG
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPT0(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTS0(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGAUX(:,:)
INTEGER(KIND=JPIM),INTENT(IN)    :: KITS0
INTEGER(KIND=JPIM),INTENT(IN)    :: KSD0
INTEGER(KIND=JPIM),INTENT(IN)    :: KTXI
INTEGER(KIND=JPIM),INTENT(IN)    :: KTNI
INTEGER(KIND=JPIM),INTENT(IN)    :: KSW0
INTEGER(KIND=JPIM),INTENT(IN)    :: KSF0
INTEGER(KIND=JPIM),INTENT(IN)    :: KHV0
INTEGER(KIND=JPIM),INTENT(IN)    :: KSR0
INTEGER(KIND=JPIM),INTENT(IN)    :: KTCLS
INTEGER(KIND=JPIM),INTENT(IN)    :: KRHCLS
INTEGER(KIND=JPIM),INTENT(IN)    :: KTR0
INTEGER(KIND=JPIM),INTENT(IN)    :: KHUXI
INTEGER(KIND=JPIM),INTENT(IN)    :: KHUNI
INTEGER(KIND=JPIM),INTENT(IN)    :: KFPCUFNR
INTEGER(KIND=JPIM),INTENT(IN)    :: KCUFCOD(:)
INTEGER(KIND=JPIM),INTENT(IN)    :: KCUFT2(:)
INTEGER(KIND=JPIM),INTENT(IN)    :: KEDR
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCFIS(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCLSM(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCVEG(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCIVE(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCARG(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCD2(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCSAB(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGM(KPROMA)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PGP(KPROMA,KAUX)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFP(KPROMA,KGT1)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMASK(KPROMA)
!     ------------------------------------------------------------------
! === LOGICAL VARIABLES COMPUTED IN PART 1.1.4 =========================
LOGICAL :: LL_ETH   ! equivalent potential temperature
LOGICAL :: LL_IET   ! isobaric equivalent temperature
LOGICAL :: LL_THV   ! virtual theta
LOGICAL :: LL_PD    ! NH pressure departure
LOGICAL :: LL_HU    ! Relative humidity
LOGICAL :: LL_TH    ! Potential temperature
LOGICAL :: LL_THPW  ! Wet bulb potential temperature
LOGICAL :: LL_VD    ! Vertical divergence
LOGICAL :: LL_VW    ! True vertical velocity
LOGICAL :: LL_WND   ! wind speed
LOGICAL :: LL_WWS   ! surface true vertical velocity
LOGICAL :: LL_Z     ! Geopotential
LOGICAL :: LL_PTB   ! Pressure of iso-T
LOGICAL :: LL_HTB   ! height of iso-T
LOGICAL :: LL_HTPW  ! height of iso-T'w=0 Celsius, from bottom
LOGICAL :: LL_HTPW1 ! height of iso-T'w=1 Celsius, from bottom
LOGICAL :: LL_HTPW2 ! height of iso-T'w=1.5 Celsius, from bottom
LOGICAL :: LL_TPW   ! Wet bulb temperature
LOGICAL :: LL_UJET  ! ICAO U component of the jet
LOGICAL :: LL_VJET  ! ICAO V component of the jet
LOGICAL :: LL_PJET  ! ICAO pressure at jet level
LOGICAL :: LL_PCAO  ! ICAO tropopause pressure
LOGICAL :: LL_TCAO  ! ICAO tropopause temperature
LOGICAL :: LL_RHO   ! Humid density
LOGICAL :: LL_SRE   ! simulated reflectivities in mm/h
LOGICAL :: LL_SREDB ! simulated reflectivities in dBZ
LOGICAL :: LL_SRH   ! to activate computation of storm relative helicity and storm motion diagnostics
LOGICAL :: LL_SREX  ! max simu. reflectivities in mm/h
LOGICAL :: LL_SREDBX! max simu. reflectivities in dBZ
LOGICAL :: LL_SREDBC! compute simu. reflectivities in dBZ
LOGICAL :: LL_SREC  ! compute simu. reflectivities in mm/h

! === OTHER LOGICAL VARIABLES ==========================================
LOGICAL :: LLPBL ! to activate PBL computation
LOGICAL :: LLICV ! to activate the computation of CAPE & CIEN & TCVS
LOGICAL :: LL1,LL2
LOGICAL :: LLWFULL ! to compute w at full levels
LOGICAL :: LLWHALF ! to compute w at half levels
LOGICAL :: LLGEOP  ! to compute geopotential
LOGICAL :: LLNHPD  ! to compute pressure departure (NH)
LOGICAL :: LLSPCH  ! to interpolate Q rather than Hu
LOGICAL,ALLOCATABLE :: LL_GFL(:) ! GFL
LOGICAL,ALLOCATABLE :: LLS1_GFL(:),LLS2_GFL(:),LLIN_GFL(:)
LOGICAL :: LLTOP, LLAERO, LLUNIT, LLAPA(2)

! === REAL ARRAYS COMPUTED IN PART 1.2 =================================
REAL(KIND=JPRB),ALLOCATABLE :: ZGFLT0(:,:,:)
REAL(KIND=JPRB) :: ZUT0 (KPROMA,YDGEOMETRY%YRDIMV%NFLEVG), ZVT0 (KPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZTT0 (KPROMA,YDGEOMETRY%YRDIMV%NFLEVG),ZGRHLT0(KPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZDEP0(KPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZOROG(KPROMA)

! === REAL ARRAYS COMPUTED IN PART 1.3 =================================
REAL(KIND=JPRB) :: ZI_PRESH(KPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG), ZI_PRESF(KPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZI_XYB(KPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YYTXYB%NDIM)
REAL(KIND=JPRB) :: ZI_GEOPH(KPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG),ZI_GEOPF(KPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZI_R(KPROMA,YDGEOMETRY%YRDIMV%NFLEVG),ZI_RRED(KPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZI_CP(KPROMA,YDGEOMETRY%YRDIMV%NFLEVG),ZI_KAPA(KPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZI_RHF(KPROMA,YDGEOMETRY%YRDIMV%NFLEVG),ZUSL_ES(KPROMA,YDGEOMETRY%YRDIMV%NFLEVG)

REAL(KIND=JPRB) :: ZI_DELP(KPROMA,YDGEOMETRY%YRDIMV%NFLEVG)

! === REAL ARRAYS COMPUTED IN PARTS 1.4 AND 1.5 ========================
REAL(KIND=JPRB) :: ZTSI  (KPROMA) ! interpolated surface temperature without
                                  ! altitude correction
REAL(KIND=JPRB) :: ZTSA  (KPROMA) ! interpolated surface temperature WITH
                                  ! altitude correction
REAL(KIND=JPRB) :: ZTSO  (KPROMA) ! output surface temperature
REAL(KIND=JPRB) :: ZTSEA (KPROMA) ! sea surface temperature
REAL(KIND=JPRB) :: ZSPP(KPROMA)   ! surface pressure
REAL(KIND=JPRB) :: ZLNSP(KPROMA)  ! Ln of surface pressure

! === OUTPUT OF PART 2 =================================================
! Arrays dimensioned with MAX(NFLEVG,KOPLEV) are used in both parts 2.1 and 2.2;
! arrays dimensioned with NFLEVG (resp. KOPLEV) are used in part 2.1 (resp. 2.2)
REAL(KIND=JPRB), ALLOCATABLE :: ZWHPP2(:,:,:) ! w half levels
REAL(KIND=JPRB), ALLOCATABLE :: ZPREH2(:,:,:) ! half levels pressures
REAL(KIND=JPRB) :: ZXLEVG(YDGEOMETRY%YRDIMV%NFLEVG) ! level indexes
REAL(KIND=JPRB) :: ZPREPP(KPROMA,MAX(YDGEOMETRY%YRDIMV%NFLEVG,KOPLEV)) ! upper air pressure
REAL(KIND=JPRB) :: ZUPP(KPROMA,MAX(YDGEOMETRY%YRDIMV%NFLEVG,KOPLEV))   ! U
REAL(KIND=JPRB) :: ZVPP(KPROMA,MAX(YDGEOMETRY%YRDIMV%NFLEVG,KOPLEV))   ! V
REAL(KIND=JPRB) :: ZTPP(KPROMA,MAX(YDGEOMETRY%YRDIMV%NFLEVG,KOPLEV))   ! T
REAL(KIND=JPRB) :: ZPDEPP(KPROMA,MAX(YDGEOMETRY%YRDIMV%NFLEVG,KOPLEV)) ! NH pressure departure
REAL(KIND=JPRB) :: ZWFPP(KPROMA,MAX(YDGEOMETRY%YRDIMV%NFLEVG,KOPLEV))  ! w full levels
REAL(KIND=JPRB),ALLOCATABLE :: ZGFLPP(:,:,:)         ! GFL
REAL(KIND=JPRB) :: ZPHIPP(KPROMA,MAX(YDGEOMETRY%YRDIMV%NFLEVG,KOPLEV)) ! geopotential
REAL(KIND=JPRB) :: ZHUPP(KPROMA,MAX(YDGEOMETRY%YRDIMV%NFLEVG,KOPLEV))  ! Hu
REAL(KIND=JPRB) :: ZGRHLPP(KPROMA,MAX(YDGEOMETRY%YRDIMV%NFLEVG,KOPLEV))! Hail+Graupel
REAL(KIND=JPRB) :: ZWSPP(KPROMA)            ! surface true vertical velocity
REAL(KIND=JPRB) :: ZKAPPP(KPROMA,MAX(YDGEOMETRY%YRDIMV%NFLEVG,KOPLEV)) ! Kappa on output full levels
REAL(KIND=JPRB) :: ZRPP(KPROMA,MAX(YDGEOMETRY%YRDIMV%NFLEVG,KOPLEV))   ! R on output full levels
REAL(KIND=JPRB) :: ZCPPP(KPROMA,MAX(YDGEOMETRY%YRDIMV%NFLEVG,KOPLEV))  ! Cp on output full levels
REAL(KIND=JPRB) :: ZPRES_PPH(KPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG),ZPRES_PPF(KPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZXYB_PP(KPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YYTXYB%NDIM)
REAL(KIND=JPRB) :: ZRRED_PP(KPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZGEO_PPH(KPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG),ZGEO_PPF(KPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZTETAPP(KPROMA,KOPLEV)   ! Theta
REAL(KIND=JPRB) :: ZEPTHPP(KPROMA,KOPLEV)   ! equivalent potential temperature
REAL(KIND=JPRB) :: ZIBETPP(KPROMA,KOPLEV)   ! virtual potential temperature
REAL(KIND=JPRB) :: ZTHPWPP(KPROMA,KOPLEV)   ! wet bulb potential temperature
REAL(KIND=JPRB) :: ZWINDPP(KPROMA,KOPLEV)   ! horizontal wind velocity
REAL(KIND=JPRB) :: ZSIMRPP(KPROMA,KOPLEV)   ! simulated reflectivities in mm/h
REAL(KIND=JPRB) :: ZSIMRDBPP(KPROMA,KOPLEV) ! simulated reflectivities in dBZ
REAL(KIND=JPRB) :: ZVIRTPP(KPROMA,KOPLEV)   ! virtual theta
REAL(KIND=JPRB) :: ZMGDWPP(KPROMA,KOPLEV)   ! -G*(Delta w)
REAL(KIND=JPRB) :: ZRHOPP(KPROMA,KOPLEV)    ! Humid density
REAL(KIND=JPRB) :: ZTTPWPP(KPROMA,KOPLEV)   ! T'w,wet bulb temperature, 3D field
REAL(KIND=JPRB) :: ZTPWPP(KPROMA,MAX(YDGEOMETRY%YRDIMV%NFLEVG,KOPLEV)) ! T'w for iso
REAL(KIND=JPRB) :: ZPEORG(KPROMA,MAX(YDGEOMETRY%YRDIMV%NFLEVG,KOPLEV)) ! pressures of data origin
REAL(KIND=JPRB) :: ZR2PP(KPROMA,MAX(YDGEOMETRY%YRDIMV%NFLEVG,KOPLEV))  ! ?
REAL(KIND=JPRB) :: ZQCLS(KPROMA) ! Q of pbl
REAL(KIND=JPRB) :: ZRHCLS(KPROMA)! RH of pbl
REAL(KIND=JPRB) :: ZTCLS(KPROMA) ! T of pbl
REAL(KIND=JPRB) :: ZUCLS(KPROMA) ! U of pbl
REAL(KIND=JPRB) :: ZVCLS(KPROMA) ! V of pbl
REAL(KIND=JPRB) :: ZNUCLS(KPROMA)! neutral U of pbl
REAL(KIND=JPRB) :: ZNVCLS(KPROMA)! neutral V of pbl
REAL(KIND=JPRB) :: ZFCLS(KPROMA) ! wind velocity of pbl
REAL(KIND=JPRB) :: ZPCLS(KPROMA) ! pressure of pbl
REAL(KIND=JPRB) :: ZCAPE(KPROMA) ! CAPE
REAL(KIND=JPRB) :: ZCIN(KPROMA)  ! CIN
REAL(KIND=JPRB) :: ZTCVS(KPROMA) ! Convective Temperature
REAL(KIND=JPRB) :: ZUGST(KPROMA) ! U of gusts
REAL(KIND=JPRB) :: ZVGST(KPROMA) ! V of gusts
REAL(KIND=JPRB) :: ZFGST(KPROMA) ! gusts
INTEGER(KIND=JPIM) :: ILCL(KPROMA) ! Lifting condensation level
INTEGER(KIND=JPIM) :: IFCL(KPROMA) ! Free convection level
INTEGER(KIND=JPIM) :: IEL(KPROMA)  ! Equilibrium level
REAL(KIND=JPRB) :: ZLCL(KPROMA) !  Lifting condensation level (height)
REAL(KIND=JPRB) :: ZFCL(KPROMA) ! Free convection level (height)
REAL(KIND=JPRB) :: ZEL(KPROMA)  ! Equilibrium level (height)
REAL(KIND=JPRB) :: ZHTPW(KPROMA) ! height of iso-T'w=0 Celsius, from bottom
REAL(KIND=JPRB) :: ZHTPW1(KPROMA)! height of iso-T'w=1 Celsius, from bottom
REAL(KIND=JPRB) :: ZHTPW2(KPROMA)! height of iso-T'w=1.5 Celsius, from bottom
REAL(KIND=JPRB) :: ZPICAO(KPROMA)! ICAO tropopause - pressure
REAL(KIND=JPRB) :: ZTICAO(KPROMA)! ICAO tropopause - temperature
REAL(KIND=JPRB) :: ZUJET(KPROMA) ! ICAO jet - U momentum
REAL(KIND=JPRB) :: ZVJET(KPROMA) ! ICAO jet - V momentum
REAL(KIND=JPRB) :: ZPJET(KPROMA) ! ICAO jet - pressure
REAL(KIND=JPRB) :: ZTMIN(KPROMA) ! Tcls min
REAL(KIND=JPRB) :: ZTMAX(KPROMA) ! Tcls max
REAL(KIND=JPRB) :: ZHUMIN(KPROMA)! Hu cls min
REAL(KIND=JPRB) :: ZHUMAX(KPROMA)! Hu cls max
REAL(KIND=JPRB) :: ZSREMAX(KPROMA)! Simu reflectivity max in mm/h
REAL(KIND=JPRB) :: ZSREDBMAX(KPROMA)! Simu reflectivity max in dBZ
REAL(KIND=JPRB) :: ZGEOP(KPROMA,KOPLEV) ! geopotentials of the height levels
REAL(KIND=JPRB) :: ZPSIPP(KPROMA,KOPLEV)  ! velocity potential
REAL(KIND=JPRB) :: ZKHIPP(KPROMA,KOPLEV)  ! stream function
REAL(KIND=JPRB) :: ZPTB(KPROMA,KOPLEV)   ! Pressure of iso-T
REAL(KIND=JPRB) :: ZHTB(KPROMA,KOPLEV)   ! Height of iso-T
REAL(KIND=JPRB) :: ZQS(KPROMA)           ! surface pressure

REAL(KIND=JPRB) :: ZSTRMMU(KPROMA),ZSTRMMV(KPROMA) ! storm motion (u,v)
REAL(KIND=JPRB) :: ZSRH(KPROMA)  ! storm-relative helicity

! === OTHER LOCAL VARIABLES ============================================

REAL(KIND=JPRB), PARAMETER :: PPROOF=30000._JPRB

REAL(KIND=JPRB) :: ZPPIN(KPROMA,MAX(YDGEOMETRY%YRDIMV%NFLEVG,KOPLEV))
REAL(KIND=JPRB) :: Z1SGM(KPROMA) ! Reverse of map factor
REAL(KIND=JPRB) :: ZEW(KPROMA)
REAL(KIND=JPRB) :: ZRHBNDS(2)

REAL(KIND=JPRB) :: ZDTDPHI  ! dT/(g*dz1)
REAL(KIND=JPRB) :: Z1SG     ! 1/g
REAL(KIND=JPRB) :: ZXTEMP, ZTMSK

! Critical Thickness of PBL
REAL(KIND=JPRB), PARAMETER :: PPBL=17500._JPRB

INTEGER(KIND=JPIM) :: JL, JI, JV, JLEVP, ILOCU, ILOCV, JCUF
INTEGER(KIND=JPIM) :: ILOC, ILPT, I, ICAPE, IPOL2, IMAXLEV, ITOP, IPPGFL, IEDR
INTEGER(KIND=JPIM),ALLOCATABLE :: IGFLCOD(:),IDIN_GFL(:)

TYPE(TYPE_IGFLFLD) :: YLPPIC

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "apache.intfb.h"
#include "ctstar.intfb.h"
#include "endpos_prepgfl.intfb.h"
#include "fpachmt.intfb.h"
#include "fpcica.intfb.h"
#include "fpps.intfb.h"
#include "gpept.intfb.h"
#include "gpgeo.intfb.h"
#include "gpgeo_sacc.intfb.h"
#include "gpiet.intfb.h"
#include "gphpre.intfb.h"
#include "gpprs0d.intfb.h"
#include "gprcp.intfb.h"
#include "gprh.intfb.h"
#include "gptet.intfb.h"






#include "ppleta.intfb.h"
#include "ppthpw.intfb.h"
#include "ppltw.intfb.h"
#include "fpsrh.intfb.h"
#include "fpstrmm.intfb.h"
#include "fcttrm.func.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ENDPOS',0,ZHOOK_HANDLE)
ASSOCIATE(TFP=>YDAFN%TFP, TFP_DYNDS=>YDAFN%TFP_DYNDS,YDCVER=>YDGEOMETRY%YRCVER, &
 & YDVAB=>YDGEOMETRY%YRVAB, YDPHY=>YDML_PHY_MF%YRPHY, YGFL=>YDML_GCONF%YGFL,YDDIMF=>YDML_GCONF%YRDIMF)
ASSOCIATE(LFPRH100=>YDNAMFPSCI%LFPRH100, LFPQ=>YDNAMFPSCI%LFPQ, NFPCAPE=>YDNAMFPSCI%NFPCAPE, &
 & FPRHMIN=>YDNAMFPSCI%FPRHMIN, FPRHMAX=>YDNAMFPSCI%FPRHMAX, LFPCAPEX=>YDNAMFPSCI%LFPCAPEX, &
 & LFPCLSTOGMV=>YDNAMFPSCI%LFPCLSTOGMV, RENTRA=>YDNAMFPSCI%RENTRA, &
 & LISOT_ABOVEG=>YDNAMFPSCI%LISOT_ABOVEG, &
 & NAERO=>YGFL%NAERO, &
 & NCHEM=>YGFL%NCHEM, NGFL_EXT=>YGFL%NGFL_EXT, NGHG=>YGFL%NGHG, &
 & NFLEVG=>YDGEOMETRY%YRDIMV%NFLEVG, &
 & STPRE=>YDGEOMETRY%YRSTA%STPRE, LMPHYS=>YDPHY%LMPHYS, LNEIGE=>YDPHY%LNEIGE)

!     ------------------------------------------------------------------

!*       1.    INITIALIZATIONS.
!              ---------------

!*    1.1  LOGICAL AND OTHER SCALAR INITIALIZATIONS

!     1.1.1  Define the "standard order" and pointers for post-processed GFL

! we compute here pointers YLPPIC%I[X] to locate GFL in LL_GFL and ZGFLPP,
! and the number of 2D fields (IPPGFL) we have to store in LL_GFL and ZGFLPP.
! Assumptions:
! - memory is conditionally allocated for multi-fields GFL.
! - for the other GFL variables allocation is always done for the time being,
!   even for not post-processed fields.
! These pointers will be re-used for the non-derivative fields stored in ZGFLT0.

IPPGFL=0
CALL SUPTRPPGFL_ENDPOS(1,IPPGFL,YLPPIC%IQ)
CALL SUPTRPPGFL_ENDPOS(1,IPPGFL,YLPPIC%IL)
CALL SUPTRPPGFL_ENDPOS(1,IPPGFL,YLPPIC%ILRAD)
CALL SUPTRPPGFL_ENDPOS(1,IPPGFL,YLPPIC%ILCONV)
CALL SUPTRPPGFL_ENDPOS(1,IPPGFL,YLPPIC%II)
CALL SUPTRPPGFL_ENDPOS(1,IPPGFL,YLPPIC%IIRAD)
CALL SUPTRPPGFL_ENDPOS(1,IPPGFL,YLPPIC%IICONV)
CALL SUPTRPPGFL_ENDPOS(1,IPPGFL,YLPPIC%IA)
CALL SUPTRPPGFL_ENDPOS(1,IPPGFL,YLPPIC%IO3)
CALL SUPTRPPGFL_ENDPOS(1,IPPGFL,YLPPIC%ICPF)
CALL SUPTRPPGFL_ENDPOS(1,IPPGFL,YLPPIC%ISPF)
CALL SUPTRPPGFL_ENDPOS(1,IPPGFL,YLPPIC%IS)
CALL SUPTRPPGFL_ENDPOS(1,IPPGFL,YLPPIC%ISCONV)
CALL SUPTRPPGFL_ENDPOS(1,IPPGFL,YLPPIC%IRR)
CALL SUPTRPPGFL_ENDPOS(1,IPPGFL,YLPPIC%IRCONV)
CALL SUPTRPPGFL_ENDPOS(1,IPPGFL,YLPPIC%IG)
CALL SUPTRPPGFL_ENDPOS(1,IPPGFL,YLPPIC%IH)
CALL SUPTRPPGFL_ENDPOS(1,IPPGFL,YLPPIC%ITKE)
CALL SUPTRPPGFL_ENDPOS(1,IPPGFL,YLPPIC%IUAL)
CALL SUPTRPPGFL_ENDPOS(1,IPPGFL,YLPPIC%IUOM)
CALL SUPTRPPGFL_ENDPOS(1,IPPGFL,YLPPIC%IDAL)
CALL SUPTRPPGFL_ENDPOS(1,IPPGFL,YLPPIC%IDOM)
CALL SUPTRPPGFL_ENDPOS(NGHG,IPPGFL,YLPPIC%IGHG)
CALL SUPTRPPGFL_ENDPOS(NAERO,IPPGFL,YLPPIC%IAERO)
CALL SUPTRPPGFL_ENDPOS(NCHEM,IPPGFL,YLPPIC%ICHEM)
CALL SUPTRPPGFL_ENDPOS(NGFL_EXT,IPPGFL,YLPPIC%IEXT)

! Extra 3D physical diagnostics
IEDR=IPPGFL+1
IPPGFL=IPPGFL+1

!     1.1.2  LLICV,LLPBL,LL1,LL2.

LLICV = YDQTYPE%LL(TFP%CAPE%ICOD) .OR. YDQTYPE%LL(TFP%CIEN%ICOD)&
 & .OR. YDQTYPE%LL(TFP%LCL%ICOD)  .OR. YDQTYPE%LL(TFP%FCL%ICOD)&
 & .OR. YDQTYPE%LL(TFP%EL%ICOD)   .OR. YDQTYPE%LL(TFP%TCVS%ICOD)

LLPBL = YDQTYPE%LL(TFP%UCLS%ICOD) .OR. YDQTYPE%LL(TFP%VCLS%ICOD)&
 & .OR. YDQTYPE%LL(TFP%TCLS%ICOD) .OR. YDQTYPE%LL(TFP%QCLS%ICOD)&
 & .OR. YDQTYPE%LL(TFP%RCLS%ICOD) .OR. YDQTYPE%LL(TFP%FCLS%ICOD)&
 & .OR. YDQTYPE%LL(TFP%TX%ICOD)   .OR. YDQTYPE%LL(TFP%TN%ICOD)&
 & .OR. YDQTYPE%LL(TFP%HUX%ICOD)  .OR. YDQTYPE%LL(TFP%HUN%ICOD)&
 & .OR. YDQTYPE%LL(TFP%UGST%ICOD) .OR. YDQTYPE%LL(TFP%VGST%ICOD)&
 & .OR. YDQTYPE%LL(TFP%FGST%ICOD) .OR. (LLICV .AND. (.NOT.LFPCAPEX))

LL_SRH= YDQTYPE%LL(TFP%SRH%ICOD)

LL1=YDQTYPE%LL(TFP%TH%ICOD).OR.YDQTYPE%LL(TFP%ETH%ICOD).OR.&
 & YDQTYPE%LL(TFP%THPW%ICOD).OR.YDQTYPE%LL(TFP%IET%ICOD).OR.&
 & YDQTYPE%LL(TFP%SRE%ICOD).OR.YDQTYPE%LL(TFP%SREDB%ICOD).OR.&
 & YDQTYPE%LL(TFP%RHO%ICOD).OR.YDQTYPE%LL(TFP%THV%ICOD).OR.&
 & YDQTYPE%LL(TFP%SREX%ICOD).OR.YDQTYPE%LL(TFP%SREDBX%ICOD)

LL2=YDQTYPE%LL(TFP%TH%ICOD).OR.YDQTYPE%LL(TFP%ETH%ICOD).OR.&
 & YDQTYPE%LL(TFP%SRE%ICOD).OR.YDQTYPE%LL(TFP%SREDB%ICOD).OR.&
 & YDQTYPE%LL(TFP%IET%ICOD).OR.YDQTYPE%LL(TFP%THV%ICOD).OR.&
 & YDQTYPE%LL(TFP%SREX%ICOD).OR.YDQTYPE%LL(TFP%SREDBX%ICOD)

!     1.1.3  LL_GFL,LLS1_GFL,LLS2_GFL,LLIN_GFL,IGFLCOD,IDIN_GFL

ALLOCATE(LL_GFL(IPPGFL))
ALLOCATE(LLS1_GFL(IPPGFL))
ALLOCATE(LLS2_GFL(IPPGFL))
ALLOCATE(LLIN_GFL(IPPGFL))
ALLOCATE(IGFLCOD(IPPGFL))
ALLOCATE(IDIN_GFL(IPPGFL))

CALL ENDPOS_PREPGFL(YDQTYPE,TFP,KLMOD,NFLEVG,YGFL,IPPGFL,YDIN_GFL,YLPPIC,LL1,LL2,LLPBL,LLICV,&
 & LL_GFL,LLS1_GFL,LLS2_GFL,LLIN_GFL,IGFLCOD,IDIN_GFL)

LL_GFL  (IEDR)=YDQTYPE%LL(TFP%EDR%ICOD)
LLS1_GFL(IEDR)=.FALSE.
LLS2_GFL(IEDR)=.FALSE.
LLIN_GFL(IEDR)=KLMOD(TFP%EDR%ICOD) >= 1
IGFLCOD (IEDR)=TFP%EDR%ICOD
IDIN_GFL(IEDR)=KEDR

!     1.1.4  LL_[X] for other variables

LL_ETH  = YDQTYPE%LL(TFP%ETH%ICOD)
LL_IET  = YDQTYPE%LL(TFP%IET%ICOD)
LL_SRE  = YDQTYPE%LL(TFP%SRE%ICOD)
LL_SREDB= YDQTYPE%LL(TFP%SREDB%ICOD)
LL_THV  = YDQTYPE%LL(TFP%THV%ICOD)
LL_PD   = YDQTYPE%LL(TFP%PD%ICOD)
LL_HU   = YDQTYPE%LL(TFP%HU%ICOD)
LL_TH   = YDQTYPE%LL(TFP%TH%ICOD)
LL_THPW = YDQTYPE%LL(TFP%THPW%ICOD)
LL_TPW  = YDQTYPE%LL(TFP%TPW%ICOD)
LL_VD   = YDQTYPE%LL(TFP%VD%ICOD)
LL_VW   = YDQTYPE%LL(TFP%VW%ICOD)
LL_WND  = YDQTYPE%LL(TFP%WND%ICOD)
LL_WWS  = YDQTYPE%LL(TFP%WWS%ICOD)
LL_Z    = YDQTYPE%LL(TFP%Z%ICOD)

LL_PTB = YDQTYPE%LL(TFP%PTB%ICOD)
LL_HTB = YDQTYPE%LL(TFP%HTB%ICOD)
LL_HTPW = YDQTYPE%LL(TFP%HTPW%ICOD)
LL_HTPW1 = YDQTYPE%LL(TFP%HTPW1%ICOD)
LL_HTPW2 = YDQTYPE%LL(TFP%HTPW2%ICOD)
LL_VJET = YDQTYPE%LL(TFP%VJET%ICOD)
LL_PJET = YDQTYPE%LL(TFP%PJET%ICOD)
LL_PCAO = YDQTYPE%LL(TFP%PCAO%ICOD)
LL_TCAO = YDQTYPE%LL(TFP%TCAO%ICOD)
LL_RHO  = YDQTYPE%LL(TFP%RHO%ICOD)
LL_SREX = YDQTYPE%LL(TFP%SREX%ICOD)
LL_SREDBX= YDQTYPE%LL(TFP%SREDBX%ICOD)

LLAERO = LL_UJET .OR. LL_VJET .OR. LL_PJET .OR. LL_PCAO .OR. LL_TCAO

!     1.1.5  Other variables.

Z1SG=1.0_JPRB/RG
ZDTDPHI=RDTDZ1/RG

IMAXLEV=MAX(NFLEVG,KOPLEV)

DO JI=KST,KND
  Z1SGM(JI)=1.0_JPRB/PGM(JI)
ENDDO

ZRHBNDS(1)=FPRHMIN
ZRHBNDS(2)=FPRHMAX

ZTMSK=288._JPRB


!*    1.2  ARRAYS INITIALIZATIONS

!*      1.2.1  GFL

ALLOCATE(ZGFLT0(KPROMA,NFLEVG,IPPGFL))
DO JV=1,IPPGFL
  IF (LLIN_GFL(JV)) THEN
    ZGFLT0(:,:,JV)=PABUF(:,IDIN_GFL(JV):IDIN_GFL(JV)+NFLEVG-1)
  ELSE
    ZGFLT0(:,:,JV)=0.0_JPRB
  ENDIF
ENDDO

DO JL=1,NFLEVG
  DO JI=KST,KND
    ZGRHLT0(JI,JL)=ZGFLT0(JI,JL,YLPPIC%IG)+ZGFLT0(JI,JL,YLPPIC%IH)
  ENDDO
ENDDO

!*      1.2.2  GMV and other non-GFL quantities

DO JL=1, NFLEVG
  DO JI=KST, KND
    ZUT0(JI,JL)=PABUF(JI,KUT0+JL-1)*PGM(JI)
    ZVT0(JI,JL)=PABUF(JI,KVT0+JL-1)*PGM(JI)
    ZTT0(JI,JL)=PABUF(JI,KTT0+JL-1)
    SELECT CASE (KDEP0)
    CASE (1:)
      ZDEP0(JI,JL)=PABUF(JI,KDEP0+JL-1)
    CASE DEFAULT
      ZDEP0(JI,JL)=0.0_JPRB
    END SELECT
  ENDDO
ENDDO

DO JI=KST, KND
  ZOROG(JI)=PABUF(JI,KOROG)
ENDDO

!*    1.3   MODEL LEVEL PRESSURES AND DYNAMICS

! We find there all calls to GP, GNH routines for departure system.

! Pressures:
ZI_PRESH(KST:KND,NFLEVG)=PSPT0(KST:KND)
CALL GPHPRE(KPROMA,NFLEVG,KST,KND,YDVAB,YDCVER,ZI_PRESH,PXYB=ZI_XYB,PRESF=ZI_PRESF)

DO JL=1,NFLEVG
  DO JI=KST,KND
    ZI_DELP(JI,JL)=ZI_PRESH(JI,JL-1)-ZI_PRESH(JI,JL)
  ENDDO
ENDDO

! R:
CALL GPRCP(KPROMA,KST,KND,NFLEVG,PQ=ZGFLT0(1,1,YLPPIC%IQ),&
 & PQI=ZGFLT0(1,1,YLPPIC%II),PQL=ZGFLT0(1,1,YLPPIC%IL),&
 & PQR=ZGFLT0(1,1,YLPPIC%IRR),PQS=ZGFLT0(1,1,YLPPIC%IS),PQG=ZGRHLT0(1,1),&
 & PCP=ZI_CP,PR=ZI_R,PKAP=ZI_KAPA)

! Compute (pre/prehyd)*ZI_R for input to GPGEO.
DO JL=1,NFLEVG
  DO JI=KST,KND
    ZI_RRED(JI,JL)=(1.0_JPRB+ZDEP0(JI,JL)/ZI_PRESF(JI,JL))*ZI_R(JI,JL)
  ENDDO
ENDDO

! Geopotential height:
ZI_GEOPH(KST:KND,NFLEVG)=ZOROG(KST:KND)
!! If LSACC=.T. the geop calc. is not strictly accurate, as should be calculated via GPGEO_SACC
!! but for dignostics purposes, makes sense to output the true hyd. GEOP
CALL GPGEO(KPROMA,KST,KND,NFLEVG,ZI_GEOPH,ZI_GEOPF,ZTT0(1,1),ZI_RRED,&
  & ZI_XYB(:,:,YYTXYB%M_LNPR),ZI_XYB(:,:,YYTXYB%M_ALPH),&
  & YDGEOMETRY%YRVERT_GEOM)

! Relative humidity RH:
CALL GPRH(.FALSE.,KPROMA,KST,KND,NFLEVG,ZRHBNDS(2),ZRHBNDS(1),&
 & ZGFLT0(1,1,YLPPIC%IQ),ZTT0(1,1),ZI_PRESF,ZUSL_ES,ZI_RHF(1,1))


!*    1.4   SET OUTPUT OROGRAPHY AND SURFACE TEMPERATURES

SELECT CASE (LMPHYS.AND.KITS0 > 0)
CASE (.TRUE.)
  DO JI=KST, KND
    ZTSI(JI)=PGAUX(JI,KITS0)
    ZTSO(JI)=PTS0(JI)
  ENDDO
CASE DEFAULT
  ! Compute T* :
  CALL CTSTAR(KPROMA,KST,KND,ZTT0(1,NLEXTRAP),ZI_PRESH(1,NFLEVG),&
   & ZI_PRESF(1,NLEXTRAP),ZOROG,ZTSI,ZTSEA)
  DO JI=KST, KND
    ZTSO(JI)=ZTSI(JI)+ZDTDPHI*(PCFIS(JI)-ZOROG(JI))
  ENDDO
END SELECT

DO JI=KST, KND
  ZTSA(JI)=ZTSI(JI)+ZDTDPHI*(PCFIS(JI)-ZOROG(JI))
ENDDO

!*    1.5   PREPARE FOR INTERPOLATIONS: COMPUTE POST-PROCESSING LEVELS

CALL FPPS(KPROMA,KST,KND,NFLEVG,1,ZI_GEOPH,PCFIS,ZTT0(1,1),&
 & ZI_R,ZI_PRESH,ZTSI,ZSPP)

DO JI=KST,KND
  ZLNSP(JI)=LOG(ZSPP(JI))
ENDDO

!*    1.6   ALLOCATIONS.

ALLOCATE(ZGFLPP(KPROMA,IMAXLEV,IPPGFL))

!     ------------------------------------------------------------------

!*       2.    POST-PROCESSING: CALL TO APACHE
!              -------------------------------

!*    2.1   SURFACE & PBL POST-PROCESSING

LLAPA(1)=LL_WWS.OR.LLPBL.OR.LLICV.OR.LL_SRH
LLAPA(2)=LL_PTB.OR.LL_HTB.OR.LL_HTPW.OR.LLAERO.OR.LL_HTPW1.OR.LL_HTPW2.OR.LL_TPW

IF (LLAPA(1).OR.LLAPA(2)) THEN

  ALLOCATE(ZWHPP2(KPROMA,NFLEVG,0:1))
  ALLOCATE(ZPREH2(KPROMA,NFLEVG,0:1))

  DO JL=1,NFLEVG
    ZXLEVG(JL)=REAL(JL,JPRB)
  ENDDO
  CALL PPLETA(YDGEOMETRY%YRVAB,YDGEOMETRY%YRCVER,KPROMA,KST,KND,NFLEVG,ZXLEVG,ZSPP,ZPREPP,&
   & ZPREH2(1,1,0),ZPREH2(1,1,1))

  LLSPCH=LFPQ
  LLWFULL=.FALSE.
  LLNHPD=.FALSE.

ENDIF

IF (LLAPA(1)) THEN

  LLWHALF=LL_WWS
  LLGEOP=.FALSE.

  CALL APACHE(LNEIGE,KPROMA,KST,KND,NFLEVG,NFLEVG,IPPGFL,JPNPPM,YLPPIC,&
   & LLNHPD,LLWHALF,LLWFULL,LLS1_GFL,LLGEOP,LLSPCH,ZRHBNDS,&
   & ZOROG,PCFIS,PPBL,ZTSI,ZTSO,ZSPP,ZPREPP,ZPREH2,YDGEOMETRY%YRVERT_GEOM,&
   & ZTT0,ZUT0,ZVT0,ZDEP0,PWVEL0,ZGFLT0,&
   & ZI_PRESH,ZI_PRESF,ZI_XYB(:,:,YYTXYB%M_LNPR),ZI_XYB(:,:,YYTXYB%M_ALPH),ZI_R,ZI_KAPA,ZI_RHF,&
   & ZTPP,ZUPP,ZVPP,ZPDEPP,ZWHPP2,ZWFPP,&
   & ZGFLPP(1:KPROMA,1:NFLEVG,1:IPPGFL),ZPHIPP,ZHUPP,ZPEORG)

  IF (LLS1_GFL(YLPPIC%IG)) THEN
    ZGRHLPP(KST:KND,1:NFLEVG)=ZGFLPP(KST:KND,1:NFLEVG,YLPPIC%IG)
  ELSE
    ZGRHLPP(KST:KND,1:NFLEVG)=0._JPRB
  ENDIF
  IF (LLS1_GFL(YLPPIC%IH)) THEN
    ZGRHLPP(KST:KND,1:NFLEVG)=&
     & ZGRHLPP(KST:KND,1:NFLEVG)+ZGFLPP(KST:KND,1:NFLEVG,YLPPIC%IH)
  ENDIF

  SELECT CASE (LL_WWS)
  CASE (.TRUE.)
    DO JI=KST, KND
      ZWSPP(JI)=ZWHPP2(JI,NFLEVG,1)
    ENDDO
  END SELECT

  IF (LLPBL) THEN

    ! We find there all calls to GP, GNH routines for NFLEVG-target system
    ! (using ._PP as input data).

    ! R:
    CALL GPRCP(KPROMA,KST,KND,NFLEVG,&
     & PQ=ZGFLPP(1,1,YLPPIC%IQ),PQI=ZGFLPP(1,1,YLPPIC%II),&
     & PQL=ZGFLPP(1,1,YLPPIC%IL),PQR=ZGFLPP(1,1,YLPPIC%IRR),&
     & PQS=ZGFLPP(1,1,YLPPIC%IS),PQG=ZGRHLPP,PCP=ZCPPP,PR=ZRPP,PKAP=ZKAPPP)

    ! Pressures:
    ZPRES_PPH(KST:KND,NFLEVG)=ZSPP(KST:KND)
    CALL GPHPRE(KPROMA,NFLEVG,KST,KND,YDVAB,YDCVER,ZPRES_PPH,PXYB=ZXYB_PP,PRESF=ZPRES_PPF)

    ! Compute (pre/prehyd)*ZRPP for input to GPGEO.
    DO JL=1,NFLEVG
      DO JI=KST,KND
        ZRRED_PP(JI,JL)=(1.0_JPRB+ZPDEPP(JI,JL)/ZPRES_PPF(JI,JL))*ZRPP(JI,JL)
      ENDDO
    ENDDO

    ! Geopotential height:
    ZGEO_PPH(KST:KND,NFLEVG)=PCFIS(KST:KND)
    IF (YDDYNA%LSACC) THEN
      CALL GPGEO_SACC(KPROMA,KST,KND,NFLEVG,&
       & ZGEO_PPH,ZGEO_PPF,ZTPP,ZUPP,ZRRED_PP,ZXYB_PP(:,:,YYTXYB%M_LNPR),ZXYB_PP(:,:,YYTXYB%M_ALPH),&
       & YDGEOMETRY%YRVERT_GEOM)
    ELSE
      CALL GPGEO(KPROMA,KST,KND,NFLEVG,&
       & ZGEO_PPH,ZGEO_PPF,ZTPP,ZRRED_PP,ZXYB_PP(:,:,YYTXYB%M_LNPR),ZXYB_PP(:,:,YYTXYB%M_ALPH),&
       & YDGEOMETRY%YRVERT_GEOM)
    ENDIF
    CALL FPACHMT(YDML_PHY_MF,LDZ0H,KPROMA,KST,KND,NFLEVG,&
    & ZTPP,ZGFLPP(1,1,YLPPIC%IQ),ZUPP,ZVPP,&
    & ZPRES_PPH,ZPRES_PPF,ZCPPP,ZRPP,ZGEO_PPH,ZGEO_PPF,&
    & PCLSM,ZTSO,PGAUX(:,KSR0),PGAUX(:,KTR0),PCVEG,PCARG,PCD2,PCIVE,PCSAB,&
    & PGAUX(:,KSD0),PGAUX(:,KSW0),PGAUX(:,KSF0),PGAUX(:,KHV0),ZQCLS,ZRHCLS,&
    & ZTCLS,ZUCLS,ZVCLS,ZNUCLS,ZNVCLS,ZFCLS,ZUGST,&
    & ZVGST,ZFGST,ZPCLS,ZQS)

  ENDIF

  ! * compute storm motion diagnostics and storm relative helicity
  IF (LL_SRH) THEN
    CALL FPSTRMM(KST,KND,KPROMA,NFLEVG,&
          & ZUT0,ZVT0,ZI_DELP,ZI_PRESF,ZOROG,ZI_GEOPF,ZSTRMMU,ZSTRMMV)
    CALL FPSRH(KST,KND,KPROMA,NFLEVG, &
         & ZUT0,ZVT0,ZI_DELP,ZI_PRESF,ZOROG,ZI_GEOPF,ZSTRMMU,ZSTRMMV,ZSRH)
  ENDIF

  IF (LLICV) THEN

    ! convert model levels to height levels
    ICAPE=NFPCAPE
    IF (LFPCAPEX) THEN
      DO JI=KST,KND
        ZPCLS(JI)=ZSPP(JI)
      ENDDO
      CALL FPCICA(YDML_PHY_MF,KST,KND,KPROMA,NFLEVG,ICAPE,RENTRA,PGAUX(:,KTCLS),ZPCLS,PGAUX(:,KRHCLS),ZTPP,&
       & ZPREPP,ZGFLPP(1,1,YLPPIC%IQ),ZI_KAPA,ZCAPE,ZCIN,ZTCVS,ILCL,IFCL,IEL)
    ELSE
      CALL FPCICA(YDML_PHY_MF,KST,KND,KPROMA,NFLEVG,ICAPE,RENTRA,ZTCLS,ZPCLS,ZRHCLS,ZTPP,&
       & ZPREPP,ZGFLPP(1,1,YLPPIC%IQ),ZI_KAPA,ZCAPE,ZCIN,ZTCVS,ILCL,IFCL,IEL)
    ENDIF

    IF (YDQTYPE%LL(TFP%LCL%ICOD)) THEN
      DO JI=KST,KND
        IF (ILCL(JI) > 0 .AND. ILCL(JI) <= NFLEVG) THEN
          ZLCL(JI)=ZI_GEOPF(JI,ILCL(JI))
        ELSE
          ZLCL(JI)=0._JPRB
        ENDIF
      ENDDO
    ENDIF

    IF (YDQTYPE%LL(TFP%FCL%ICOD)) THEN
      DO JI=KST,KND
        IF (IFCL(JI) > 0 .AND. IFCL(JI) <= NFLEVG) THEN
          ZFCL(JI)=ZI_GEOPF(JI,IFCL(JI))
        ELSE
          ZFCL(JI)=0._JPRB
        ENDIF
      ENDDO
    ENDIF

    IF (YDQTYPE%LL(TFP%EL%ICOD)) THEN
      DO JI=KST,KND
        IF (IEL(JI) > 0 .AND. IEL(JI) <= NFLEVG) THEN
          ZEL(JI)=ZI_GEOPF(JI,IEL(JI))
        ELSE
          ZEL(JI)=0._JPRB
        ENDIF
      ENDDO
    ENDIF

  ENDIF

ENDIF

IF (LLAPA(2)) THEN

  LLWHALF=.FALSE.
  LLGEOP=.TRUE.

  CALL APACHE(LNEIGE,KPROMA,KST,KND,NFLEVG,NFLEVG,IPPGFL,JPNPPM,YLPPIC,&
   & LLNHPD,LLWHALF,LLWFULL,LLS2_GFL,LLGEOP,LLSPCH,ZRHBNDS,&
   & ZOROG,PCFIS,PPBL,ZTSI,ZTSA,ZSPP,ZPREPP,ZPREH2,YDGEOMETRY%YRVERT_GEOM,&
   & ZTT0,ZUT0,ZVT0,ZDEP0,PWVEL0,ZGFLT0,&
   & ZI_PRESH,ZI_PRESF,ZI_XYB(:,:,YYTXYB%M_LNPR),ZI_XYB(:,:,YYTXYB%M_ALPH),ZI_R,ZI_KAPA,ZI_RHF,&
   & ZTPP,ZUPP,ZVPP,ZPDEPP,ZWHPP2,ZWFPP,&
   & ZGFLPP(1:KPROMA,1:NFLEVG,1:IPPGFL),ZPHIPP,ZHUPP,ZPEORG)

  IF (LL_HTB .OR. LL_PTB) THEN
    IF (CDCONF == 'K') THEN
      DO JLEVP=1,KOPLEV
        ZXTEMP=ABS(PXLEV(JLEVP))
        LLTOP = (PXLEV(JLEVP) < 0)
        CALL PPLTEMP(YDGEOMETRY%YRSTA,KPROMA,KST,KND,NFLEVG,ZPHIPP,ZTPP,ZXTEMP,LLTOP,ZHTB(1,JLEVP))
      ENDDO
      IF (LL_PTB) THEN
        CALL FPPS(KPROMA,KST,KND,NFLEVG,KOPLEV,ZI_GEOPH,ZHTB,ZTT0(1,1),&
         & ZI_R,ZI_PRESH,ZTSI,ZPTB)
      ENDIF
      IF (LL_HTB) THEN
        ! Convert from geopotentials to meters
        DO JLEVP=1,KOPLEV
          ZHTB(KST:KND,JLEVP)=ZHTB(KST:KND,JLEVP)*Z1SG
          IF ( LISOT_ABOVEG ) THEN
           DO JI=KST,KND
            IF ( ZHTB(JI,JLEVP) < ZOROG(JI)*Z1SG ) THEN
              ZHTB(JI,JLEVP)=-9999._JPRB
            ENDIF
           ENDDO
          ENDIF
        ENDDO
      ENDIF
    ELSE
      CALL ABOR1('ENDPOS: UNEXPECTED CDCONF FOR ISOTHERMIC LEVELS')
    ENDIF
  ENDIF

  IF ((LL_HTPW).OR.(LL_HTPW1).OR.(LL_HTPW2)) THEN
    ITOP=1
    DO JL=1,NFLEVG
      IF (STPRE(JL) >= PPROOF) THEN
        ITOP=JL
        EXIT
      ENDIF
    ENDDO
    DO JL=ITOP,NFLEVG
      CALL PPWETPOINT(YDPHY,KST,KND,KPROMA,ZPREPP(1,JL),ZTPP(1,JL),&
       & ZGFLPP(1,JL,YLPPIC%IQ),ZGFLPP(1,JL,YLPPIC%IL),ZGFLPP(1,JL,YLPPIC%II),&
       & ZTPWPP(1,JL))
    ENDDO
    DO JL=1,ITOP-1
      DO JI=KST,KND
        ZTPWPP(JI,JL)=ZTPWPP(JI,ITOP)
      ENDDO
    ENDDO
    LLTOP=.FALSE.
    IF (LL_HTPW) THEN
      ZXTEMP=RTT
      CALL PPLTW(YDGEOMETRY%YRSTA,KPROMA,KST,KND,NFLEVG,ZPHIPP,ZTPWPP,ZXTEMP,ZHTPW)
      ! Convert from geopotentials to meters
      DO JI=KST,KND
        ZHTPW(JI)=ZHTPW(JI)*Z1SG
      ENDDO
    ENDIF
    IF (LL_HTPW1) THEN
      ZXTEMP=RTT+1.0_JPRB
      CALL PPLTW(YDGEOMETRY%YRSTA,KPROMA,KST,KND,NFLEVG,ZPHIPP,ZTPWPP,ZXTEMP,ZHTPW1)
      ! Convert from geopotentials to meters
      DO JI=KST,KND
        ZHTPW1(JI)=ZHTPW1(JI)*Z1SG
      ENDDO
    ENDIF
    IF (LL_HTPW2) THEN
      ZXTEMP=RTT+1.5_JPRB
      CALL PPLTW(YDGEOMETRY%YRSTA,KPROMA,KST,KND,NFLEVG,ZPHIPP,ZTPWPP,ZXTEMP,ZHTPW2)
      ! Convert from geopotentials to meters
      DO JI=KST,KND
        ZHTPW2(JI)=ZHTPW2(JI)*Z1SG
      ENDDO
    ENDIF
  ENDIF

  IF (LLAERO) THEN
    IPOL2=4
    CALL POAERO(ZPHIPP,ZTPP,ZUPP,ZVPP,ZPICAO,ZTICAO,ZUJET,ZVJET,ZPJET,ZPREPP,&
     & KPROMA,KST,KND,NFLEVG,IPOL2)
  ENDIF

ENDIF

IF (ALLOCATED(ZWHPP2)) DEALLOCATE(ZWHPP2)
IF (ALLOCATED(ZPREH2)) DEALLOCATE(ZPREH2)

IF (YDQTYPE%LL(TFP%TX%ICOD)) THEN
  DO JI=KST,KND
    ZTMAX(JI)=ZTCLS(JI)+PGAUX(JI,KTXI)
  ENDDO
ENDIF

IF (YDQTYPE%LL(TFP%TN%ICOD)) THEN
  DO JI=KST,KND
    ZTMIN(JI)=ZTCLS(JI)+PGAUX(JI,KTNI)
  ENDDO
ENDIF

IF (YDQTYPE%LL(TFP%HUX%ICOD)) THEN
  DO JI=KST,KND
    ZHUMAX(JI)=MAX(0.0_JPRB,MIN(1.0_JPRB,ZRHCLS(JI)+PGAUX(JI,KHUXI)))
  ENDDO
ENDIF

IF (YDQTYPE%LL(TFP%HUN%ICOD)) THEN
  DO JI=KST,KND
    ZHUMIN(JI)=MAX(0.0_JPRB,MIN(1.0_JPRB,ZRHCLS(JI)+PGAUX(JI,KHUNI)))
  ENDDO
ENDIF

!*    2.2   UPPER AIR POST-PROCESSING

IF (KOPLEV > 0) THEN

  !      2.2.1   Compute pressure on post-processing levels
  !              ------------------------------------------

  ALLOCATE(ZWHPP2(KPROMA,KOPLEV,0:1))
  ALLOCATE(ZPREH2(KPROMA,KOPLEV,0:1))
  SELECT CASE (CDCONF)
  CASE ('S')
    CALL PPLETA(YDFPVAB,YDGEOMETRY%YRCVER,KPROMA,KST,KND,KOPLEV,PXLEV,ZSPP,ZPREPP,ZPREH2(1,1,0),ZPREH2(1,1,1))
  CASE ('H')
    DO JL=1,KOPLEV
      DO JI=KST,KND
        ZGEOP(JI,JL)=RG*PXLEV(JL)+PCFIS(JI)
      ENDDO
    ENDDO
    CALL FPPS(KPROMA,KST,KND,NFLEVG,KOPLEV,ZI_GEOPH,ZGEOP,ZTT0(1,1),&
     & ZI_R,ZI_PRESH,ZTSI,ZPREPP)
  CASE ('K')
    DO JL=1,KOPLEV
        ZXTEMP=ABS(PXLEV(JL))
        LLTOP = (PXLEV(JL) < 0)
        CALL PPLTEMP(YDGEOMETRY%YRSTA,KPROMA,KST,KND,NFLEVG,ZI_GEOPF,ZTT0(1,1),ZXTEMP,LLTOP,ZGEOP(1,JL))
    ENDDO
    CALL FPPS(KPROMA,KST,KND,NFLEVG,KOPLEV,ZI_GEOPH,ZGEOP,ZTT0(1,1),&
     & ZI_R,ZI_PRESH,ZTSI,ZPREPP)
  CASE ('F')
    DO JL=1,KOPLEV
      DO JI=KST,KND
        ZGEOP(JI,JL)=RG*PXLEV(JL)
      ENDDO
    ENDDO
    CALL FPPS(KPROMA,KST,KND,NFLEVG,KOPLEV,ZI_GEOPH,ZGEOP,ZTT0(1,1),&
     & ZI_R,ZI_PRESH,ZTSI,ZPREPP)
  CASE DEFAULT
    CALL ABOR1('ENDPOS : INTERNAL ERROR CDCONF')
  END SELECT

  !      2.2.2   Vertical interpolations on basic fields
  !              ---------------------------------------

  LLGEOP =LL_Z
  LLWHALF=LL_VD
  LLWFULL=LL_VW
  LLNHPD =LL_PD.OR.LL_RHO
  LLSPCH=LFPQ

  CALL APACHE(LNEIGE,KPROMA,KST,KND,NFLEVG,KOPLEV,IPPGFL,JPNPPM,YLPPIC,&
   & LLNHPD,LLWHALF,LLWFULL,LL_GFL,LLGEOP,LLSPCH,ZRHBNDS,&
   & ZOROG,PCFIS,PPBL,ZTSI,ZTSO,ZSPP,ZPREPP,ZPREH2,YDGEOMETRY%YRVERT_GEOM,&
   & ZTT0,ZUT0,ZVT0,ZDEP0,PWVEL0,ZGFLT0,&
   & ZI_PRESH,ZI_PRESF,ZI_XYB(:,:,YYTXYB%M_LNPR),ZI_XYB(:,:,YYTXYB%M_ALPH),ZI_R,ZI_KAPA,ZI_RHF,&
   & ZTPP,ZUPP,ZVPP,ZPDEPP,ZWHPP2,ZWFPP,&
   & ZGFLPP(1:KPROMA,1:KOPLEV,1:IPPGFL),ZPHIPP,ZHUPP,ZPEORG)

  IF (LL_GFL(YLPPIC%IG)) THEN
    ZGRHLPP(KST:KND,1:KOPLEV)=ZGFLPP(KST:KND,1:KOPLEV,YLPPIC%IG)
  ELSE
    ZGRHLPP(KST:KND,1:KOPLEV)=0._JPRB
  ENDIF
  IF (LL_GFL(YLPPIC%IH)) THEN
    ZGRHLPP(KST:KND,1:KOPLEV)=&
     & ZGRHLPP(KST:KND,1:KOPLEV)+ZGFLPP(KST:KND,1:KOPLEV,YLPPIC%IH)
  ENDIF

  ZGFLPP(KST:KND,1:KOPLEV,YLPPIC%IQ)=&
    & MAX(0.0_JPRB,ZGFLPP(KST:KND,1:KOPLEV,YLPPIC%IQ))

  ZHUPP(KST:KND,1:KOPLEV)=&
    & MAX(0.0_JPRB,MIN(1.0_JPRB,ZHUPP(KST:KND,1:KOPLEV)))


  !      2.2.3   Computation of elaborated fields from basic ones
  !              ------------------------------------------------

  !      2.2.3.1 Vertical divergence quantity -G*(Delta w).

  IF (LL_VD) THEN
    DO JL=1,KOPLEV
      DO JI=KST,KND
        ZMGDWPP(JI,JL)=RG*(ZWHPP2(JI,JL,1)-ZWHPP2(JI,JL,0))
      ENDDO
    ENDDO
  ENDIF

  !      2.2.3.2 Potential and equivalent potential temperatures
  !              and isobaric equivalent temperature.

  IF (LL_TH .OR. LL_ETH .OR. LL_IET .OR. LL_THV) THEN
    ! We find there all calls to GP, GNH routines for KOPLEV-target system
    ! (using ._PP as input data).
    CALL GPRCP(KPROMA,KST,KND,KOPLEV,PQ=ZGFLPP(1,1,YLPPIC%IQ),&
     & PQI=ZGFLPP(1,1,YLPPIC%II),PQL=ZGFLPP(1,1,YLPPIC%IL),&
     & PQR=ZGFLPP(1,1,YLPPIC%IRR),&
     & PQS=ZGFLPP(1,1,YLPPIC%IS),PQG=ZGRHLPP,PCP=ZCPPP,PR=ZRPP,PKAP=ZKAPPP)
    IF (LL_TH .OR. LL_ETH.OR. LL_THV) THEN
      CALL GPTET(KPROMA,KST,KND,KOPLEV,ZPREPP,ZTPP,ZKAPPP,ZTETAPP)
      IF (LL_ETH) THEN
        CALL GPEPT(YDPHY,KPROMA,KST,KND,KOPLEV,ZTETAPP,ZTPP,ZPREPP,ZEPTHPP)
      ENDIF
    ENDIF
    IF (LL_IET) THEN
      CALL GPIET(YDPHY,KPROMA,KST,KND,KOPLEV,&
       & ZGFLPP(1,1,YLPPIC%IQ),ZTPP,ZCPPP,ZPREPP,ZIBETPP)
    ENDIF
  ENDIF

  !      2.2.3.3 Relative moisture in %

  IF (LL_HU .AND. LFPRH100) THEN
    DO JL=1,KOPLEV
      DO JI=KST,KND
        ZHUPP(JI,JL)=ZHUPP(JI,JL)*100.0_JPRB
      ENDDO
    ENDDO
  ENDIF

  !      2.2.3.4 Wet bulb potential temperature

  IF (LL_THPW) THEN
    DO JL=1,KOPLEV
      CALL PPTHPW(YDPHY,KST,KND,KPROMA,ZPREPP(1,JL),ZTPP(1,JL),ZGFLPP(1,JL,YLPPIC%IQ),&
       & ZGFLPP(1,JL,YLPPIC%IL),ZGFLPP(1,JL,YLPPIC%II),ZTHPWPP(1,JL))
    ENDDO
  ENDIF

  ! And wet bulb temperature

  IF (LL_TPW) THEN
    DO JL=1, KOPLEV
      CALL PPWETPOINT(YDPHY,KST,KND,KPROMA,ZPREPP(1,JL),ZTPP(1,JL),&
      & ZGFLPP(1,JL,YLPPIC%IQ),ZGFLPP(1,JL,YLPPIC%IL),ZGFLPP(1,JL,YLPPIC%II),&
      & ZTTPWPP(1,JL))
    ENDDO
  ENDIF


  !      2.2.3.5 Wind velocity

  IF (LL_WND) THEN
    DO JL=1,KOPLEV
      DO JI=KST,KND
        ZWINDPP(JI,JL)=SQRT(ZUPP(JI,JL)*ZUPP(JI,JL)+ZVPP(JI,JL)*ZVPP(JI,JL))
      ENDDO
    ENDDO
  ENDIF

  !     2.2.3.6. Simulated Reflectivities

  LL_SREC=LL_SRE.OR.(LL_SREX.AND.(CDCONF=='H'))
  LL_SREDBC=LL_SREDB.OR.(LL_SREDBX.AND.(CDCONF=='H'))

  IF (LL_SREC.OR.LL_SREDBC) THEN
    LLUNIT=.TRUE.
    CALL GPPRS0D(KPROMA,KST,KND,KOPLEV,ZTPP,ZGFLPP(1,1,YLPPIC%IQ),&
     & ZGFLPP(1,1,YLPPIC%IL),ZGFLPP(1,1,YLPPIC%IRR),&
     & ZGFLPP(1,1,YLPPIC%II),ZGFLPP(1,1,YLPPIC%IS),ZGRHLPP,LLUNIT,LL_SREC,LL_SREDBC,ZSIMRPP,ZSIMRDBPP)
  ENDIF
  IF (LL_SREX.AND.(CDCONF=='H')) THEN
    DO JI=KST,KND
      ZSREMAX(JI)=MAXVAL(ZSIMRPP(JI,1:KOPLEV))
    ENDDO
  ENDIF
  IF (LL_SREDBX.AND.(CDCONF=='H')) THEN
    DO JI=KST,KND
      ZSREDBMAX(JI)=MAXVAL(ZSIMRDBPP(JI,1:KOPLEV))
    ENDDO
  ENDIF

  !     2.2.3.7. Virtual theat

  IF (LL_THV) THEN





    call abor1("OIFS - ppcvirt should never be called")

  ENDIF

  !      2.2.3.8 Velocity potential & stream function

  IF (YDQTYPE%LL(TFP%PSI%ICOD)) THEN
    DO JL=1,KOPLEV
      ZPSIPP(KST:KND,JL)=ZVPP(KST:KND,JL)*Z1SGM(KST:KND)
    ENDDO
  ENDIF
  IF (YDQTYPE%LL(TFP%KHI%ICOD)) THEN
    DO JL=1,KOPLEV
      ZKHIPP(KST:KND,JL)=ZUPP(KST:KND,JL)*Z1SGM(KST:KND)
    ENDDO
  ENDIF

  !      2.2.3.9 Vorticity & Divergence

  I=TFP%DIV%ICOD
  IF (YDQTYPE%LL(I)) THEN
    IF (YDQTYPE%ISKP(I) == 2) THEN
      IF (YDQTYPE%ISF(I) >= 1) THEN
        DO JL=1,YDQTYPE%ILEV(I)
          ILPT=YDQTYPE%ILVP(JL,I)
          ILOCU=YDQTYPE%IGT1(I,1)+JL-1
          ILOCV=YDQTYPE%IGT1(I,2)+JL-1
          PFP(KST:KND,ILOCU)=ZUPP(KST:KND,ILPT)*Z1SGM(KST:KND)
          PFP(KST:KND,ILOCV)=ZVPP(KST:KND,ILPT)*Z1SGM(KST:KND)
        ENDDO
      ELSE
        CALL ABOR1('ENDPOS : INTERNAL ERROR : DIV MUST BE SPECTRAL')
      ENDIF
    ELSE
      CALL ABOR1('ENDPOS : INTERNAL ERROR : DIV/ISKP')
    ENDIF
  ENDIF

  I=TFP%VOR%ICOD
  IF (YDQTYPE%LL(I)) THEN
    IF (YDQTYPE%ISKP(I) == 2) THEN
      IF (YDQTYPE%ISF(I) >= 1) THEN
        DO JL=1,YDQTYPE%ILEV(I)
          ILPT=YDQTYPE%ILVP(JL,I)
          ILOCU=YDQTYPE%IGT1(I,1)+JL-1
          ILOCV=YDQTYPE%IGT1(I,2)+JL-1
          PFP(KST:KND,ILOCU)=ZUPP(KST:KND,ILPT)*Z1SGM(KST:KND)
          PFP(KST:KND,ILOCV)=ZVPP(KST:KND,ILPT)*Z1SGM(KST:KND)
        ENDDO
      ELSE
        CALL ABOR1('ENDPOS : INTERNAL ERROR : VOR MUST BE SPECTRAL')
      ENDIF
    ELSE
      CALL ABOR1('ENDPOS : INTERNAL ERROR : VOR/ISKP')
    ENDIF
  ENDIF

  !     2.2.3.8. Humid density

  IF (LL_RHO) THEN
    ! R as function of Qi, QL, Qv only
    CALL GPRCP(KPROMA,KST,KND,KOPLEV,&
     & PQ=ZGFLPP(1,1,YLPPIC%IQ),PQI=ZGFLPP(1,1,YLPPIC%II),&
     & PQL=ZGFLPP(1,1,YLPPIC%IL),PR=ZR2PP)
    DO JL=1,KOPLEV
      DO JI= KST,KND
      ZRHOPP(JI,JL)=(ZPREPP(JI,JL) + ZPDEPP(JI,JL))/(&
     & ZR2PP(JI,JL)*ZTPP(JI,JL))
      ENDDO
    ENDDO
  ENDIF

  DEALLOCATE(ZWHPP2)
  DEALLOCATE(ZPREH2)

  !     2.3. Overwrite T at lowest level if option lfpclstogmv

! Copy CLS fields T and Q into last model level, used for varpack when
! the last model level correspond approximatively to 2m height.
  IF (LFPCLSTOGMV) THEN
!   Isolate non-vectorizing log
    DO JI=KST,KND
      ZEW(JI)=FOEW(PGAUX(JI,KTCLS),0._JPRB)
    ENDDO
!   Lowest output model level:
    ILOC=MAXLOC(PXLEV(1:KOPLEV),DIM=1)
    DO JI=KST,KND
      ZTPP(JI,ILOC)=PGAUX(JI,KTCLS)
      ZGFLPP(JI,ILOC,YLPPIC%IQ)=(PGAUX(JI,KRHCLS)*ZEW(JI)) /(&
       & ZSPP(JI)*(RETV+1.0_JPRB)-RETV*ZEW(JI)*PGAUX(JI,KRHCLS))
    ENDDO
  ENDIF

ENDIF ! (KOPLEV > 0)

!     ------------------------------------------------------------------

!*       3.     STORE POST-PROCESSED FIELDS.
!               ----------------------------

!*    3.1   GFL FIELDS

DO JV=1,IPPGFL
  ZPPIN(1:KPROMA,1:IMAXLEV)=ZGFLPP(1:KPROMA,1:IMAXLEV,JV)
  CALL STORE_DATA(IGFLCOD(JV),PDATA2=ZPPIN)
ENDDO

!*    3.2   OTHER FIELDS

! Store data according to the derivative order and the shape of array:

CALL STORE_DATA(TFP%FIS%ICOD ,PDATA1=PCFIS )
CALL STORE_DATA(TFP%SP%ICOD  ,PDATA1=ZSPP   )
CALL STORE_DATA(TFP%LNSP%ICOD,PDATA1=ZLNSP  )
CALL STORE_DATA(TFP%WWS%ICOD ,PDATA1=ZWSPP  )
CALL STORE_DATA(TFP%UCLS%ICOD,PDATA1=ZUCLS  )
CALL STORE_DATA(TFP%VCLS%ICOD,PDATA1=ZVCLS  )
CALL STORE_DATA(TFP%TCLS%ICOD,PDATA1=ZTCLS  )
CALL STORE_DATA(TFP%QCLS%ICOD,PDATA1=ZQCLS  )
CALL STORE_DATA(TFP%RCLS%ICOD,PDATA1=ZRHCLS )
CALL STORE_DATA(TFP%FCLS%ICOD,PDATA1=ZFCLS  )
CALL STORE_DATA(TFP%TN%ICOD  ,PDATA1=ZTMIN  )
CALL STORE_DATA(TFP%TX%ICOD  ,PDATA1=ZTMAX  )
CALL STORE_DATA(TFP%HUX%ICOD ,PDATA1=ZHUMAX )
CALL STORE_DATA(TFP%HUN%ICOD ,PDATA1=ZHUMIN )
CALL STORE_DATA(TFP%SREX%ICOD,PDATA1=ZSREMAX)
CALL STORE_DATA(TFP%SREDBX%ICOD,PDATA1=ZSREDBMAX)
CALL STORE_DATA(TFP%UGST%ICOD,PDATA1=ZUGST  )
CALL STORE_DATA(TFP%VGST%ICOD,PDATA1=ZVGST  )
CALL STORE_DATA(TFP%FGST%ICOD,PDATA1=ZFGST  )
CALL STORE_DATA(TFP%CAPE%ICOD,PDATA1=ZCAPE  )
CALL STORE_DATA(TFP%CIEN%ICOD,PDATA1=ZCIN   )
CALL STORE_DATA(TFP%TCVS%ICOD,PDATA1=ZTCVS  )
CALL STORE_DATA(TFP%HTPW%ICOD,PDATA1=ZHTPW  )
CALL STORE_DATA(TFP%HTPW1%ICOD,PDATA1=ZHTPW1)
CALL STORE_DATA(TFP%HTPW2%ICOD,PDATA1=ZHTPW2)
CALL STORE_DATA(TFP%UJET%ICOD,PDATA1=ZUJET  )
CALL STORE_DATA(TFP%VJET%ICOD,PDATA1=ZVJET  )
CALL STORE_DATA(TFP%PJET%ICOD,PDATA1=ZPJET  )
CALL STORE_DATA(TFP%PCAO%ICOD,PDATA1=ZPICAO )
CALL STORE_DATA(TFP%TCAO%ICOD,PDATA1=ZTICAO )
CALL STORE_DATA(TFP%STRMMU%ICOD,PDATA1=ZSTRMMU)
CALL STORE_DATA(TFP%STRMMV%ICOD,PDATA1=ZSTRMMV)
CALL STORE_DATA(TFP%SRH%ICOD,PDATA1=ZSRH)
DO JCUF=1,KFPCUFNR
  CALL STORE_DATA(KCUFCOD(JCUF),PDATA1=PABUF(:,KCUFT2(JCUF)))
ENDDO
CALL STORE_DATA(TFP%MSK%ICOD ,PDATA1=PMASK)

CALL STORE_DATA(TFP%PTB%ICOD ,PDATA2=ZPTB  )
CALL STORE_DATA(TFP%HTB%ICOD ,PDATA2=ZHTB  )
CALL STORE_DATA(TFP%U%ICOD   ,PDATA2=ZUPP   )
CALL STORE_DATA(TFP%V%ICOD   ,PDATA2=ZVPP   )
CALL STORE_DATA(TFP%KHI%ICOD ,PDATA2=ZKHIPP )
CALL STORE_DATA(TFP%PSI%ICOD ,PDATA2=ZPSIPP )
CALL STORE_DATA(TFP%T%ICOD   ,PDATA2=ZTPP   )
CALL STORE_DATA(TFP%Z%ICOD   ,PDATA2=ZPHIPP )
CALL STORE_DATA(TFP%HU%ICOD  ,PDATA2=ZHUPP  )
CALL STORE_DATA(TFP%TH%ICOD  ,PDATA2=ZTETAPP)
CALL STORE_DATA(TFP%THPW%ICOD,PDATA2=ZTHPWPP)
CALL STORE_DATA(TFP%TPW%ICOD ,PDATA2=ZTTPWPP)
CALL STORE_DATA(TFP%ETH%ICOD ,PDATA2=ZEPTHPP)
CALL STORE_DATA(TFP%IET%ICOD ,PDATA2=ZIBETPP)
CALL STORE_DATA(TFP%WND%ICOD ,PDATA2=ZWINDPP)
CALL STORE_DATA(TFP%P%ICOD   ,PDATA2=ZPREPP )
CALL STORE_DATA(TFP%PD%ICOD  ,PDATA2=ZPDEPP )
CALL STORE_DATA(TFP%VD%ICOD  ,PDATA2=ZMGDWPP)
CALL STORE_DATA(TFP%VW%ICOD  ,PDATA2=ZWFPP  )
CALL STORE_DATA(TFP%SRE%ICOD ,PDATA2=ZSIMRPP)
CALL STORE_DATA(TFP%SREDB%ICOD ,PDATA2=ZSIMRDBPP)
CALL STORE_DATA(TFP%THV%ICOD ,PDATA2=ZVIRTPP)
CALL STORE_DATA(TFP%RHO%ICOD ,PDATA2=ZRHOPP)
CALL STORE_DATA(TFP%LCL%ICOD ,PDATA1=ZLCL)
CALL STORE_DATA(TFP%FCL%ICOD ,PDATA1=ZFCL)
CALL STORE_DATA(TFP%EL%ICOD ,PDATA1=ZEL)

DO JV=1,JPOSVX2
  I=TFP%FUA(JV)%ICOD
  IF (YDQTYPE%LL(I)) THEN
    IF (YDQTYPE%ISF(I) >= 1) THEN
      DO JL=1,YDQTYPE%ILEV(I)
        ILOC=YDQTYPE%IGT1(I,1)+JL-1
        PFP(KST:KND,ILOC)=0.0_JPRB
      ENDDO
    ELSE
      DO JL=1,YDQTYPE%ILEV(I)
        ILOC=YDQTYPE%IGPX(I)+JL-1
        PGP(KST:KND,ILOC)=0.0_JPRB
      ENDDO
    ENDIF

  ENDIF
ENDDO

!     ------------------------------------------------------------------

!*       4.     DEALLOCATIONS.
!               --------------

IF (ALLOCATED(ZGFLPP)) DEALLOCATE(ZGFLPP)
IF (ALLOCATED(IGFLCOD)) DEALLOCATE(IGFLCOD)
IF (ALLOCATED(IDIN_GFL)) DEALLOCATE(IDIN_GFL)
IF (ALLOCATED(LL_GFL)) DEALLOCATE(LL_GFL)
IF (ALLOCATED(LLS1_GFL)) DEALLOCATE(LLS1_GFL)
IF (ALLOCATED(LLS2_GFL)) DEALLOCATE(LLS2_GFL)
IF (ALLOCATED(LLIN_GFL)) DEALLOCATE(LLIN_GFL)
IF (ALLOCATED(ZGFLT0)) DEALLOCATE(ZGFLT0)

!------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('ENDPOS',1,ZHOOK_HANDLE)
!------------------------------------------------------------------
CONTAINS

SUBROUTINE STORE_DATA(K,PDATA1,PDATA2)

!------------------------------------------------------------------
INTEGER(KIND=JPIM), INTENT(IN) :: K
REAL(KIND=JPRB)   , INTENT(IN), OPTIONAL :: PDATA1(:)
REAL(KIND=JPRB)   , INTENT(IN), OPTIONAL :: PDATA2(:,:)
!------------------------------------------------------------------
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ENDPOS:STORE_DATA',0,ZHOOK_HANDLE)
ASSOCIATE(TFP_DYNDS=>YDAFN%TFP_DYNDS)
!------------------------------------------------------------------
 IF (YDQTYPE%LL(K)) THEN
   IF (PRESENT(PDATA1)) THEN
     SELECT CASE (ABS(TFP_DYNDS(K)%IORDR))
     CASE (0)
       IF (YDQTYPE%ISF(K) >= 1) THEN
         ILOC=YDQTYPE%IGT1(K,1)
         PFP(KST:KND,ILOC) = PDATA1(KST:KND)
       ELSE
         ILOC=YDQTYPE%IGPX(K)
         PGP(KST:KND,ILOC) = PDATA1(KST:KND)
       ENDIF
     CASE (1)
       IF (YDQTYPE%ISF(K) >= 1) THEN
         ILOC=YDQTYPE%IGT1(K,1)
         PFP(KST:KND,ILOC) = PDATA1(KST:KND)*Z1SGM(KST:KND)
       ELSE
         ILOC=YDQTYPE%IGPX(K)
         PGP(KST:KND,ILOC) = PDATA1(KST:KND)*Z1SGM(KST:KND)
       ENDIF
     CASE (2)
       IF (YDQTYPE%ISF(K) >= 1) THEN
         ILOC=YDQTYPE%IGT1(K,1)
         PFP(KST:KND,ILOC) = PDATA1(KST:KND)*(Z1SGM(KST:KND)*Z1SGM(KST:KND))
       ELSE
         ILOC=YDQTYPE%IGPX(K)
         PGP(KST:KND,ILOC) = PDATA1(KST:KND)*(Z1SGM(KST:KND)*Z1SGM(KST:KND))
       ENDIF
     END SELECT
   ELSEIF (PRESENT(PDATA2)) THEN
     DO JL=1,YDQTYPE%ILEV(K)
       ILPT=YDQTYPE%ILVP(JL,K)
       SELECT CASE (ABS(TFP_DYNDS(K)%IORDR))
       CASE (0)
         IF (YDQTYPE%ISF(K) >= 1) THEN
           ILOC=YDQTYPE%IGT1(K,1)+JL-1
           PFP(KST:KND,ILOC) = PDATA2(KST:KND,ILPT)
         ELSE
           ILOC=YDQTYPE%IGPX(K)+JL-1
           PGP(KST:KND,ILOC) = PDATA2(KST:KND,ILPT)
         ENDIF
       CASE (1)
         IF (YDQTYPE%ISF(K) >= 1) THEN
           ILOC=YDQTYPE%IGT1(K,1)+JL-1
           PFP(KST:KND,ILOC) = PDATA2(KST:KND,ILPT)*Z1SGM(KST:KND)
         ELSE
           ILOC=YDQTYPE%IGPX(K)+JL-1
           PGP(KST:KND,ILOC) = PDATA2(KST:KND,ILPT)*Z1SGM(KST:KND)
         ENDIF
       CASE (2)
         IF (YDQTYPE%ISF(K) >= 1) THEN
           ILOC=YDQTYPE%IGT1(K,1)+JL-1
           PFP(KST:KND,ILOC) = PDATA2(KST:KND,ILPT)*(Z1SGM(KST:KND)*Z1SGM(KST:KND))
         ELSE
           ILOC=YDQTYPE%IGPX(K)+JL-1
           PGP(KST:KND,ILOC) = PDATA2(KST:KND,ILPT)*(Z1SGM(KST:KND)*Z1SGM(KST:KND))
         ENDIF
       END SELECT
     ENDDO
   ELSE
     CALL ABOR1('ENPOS.STORE_DATA : INTERFACE ERROR')
   ENDIF
ENDIF
!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('ENDPOS:STORE_DATA',1,ZHOOK_HANDLE)
END SUBROUTINE STORE_DATA

!     ------------------------------------------------------------------

SUBROUTINE SUPTRPPGFL_ENDPOS(KFL,KPPGFL,KMPP)
! computes pointers to store interpolated GFL variables.

INTEGER(KIND=JPIM), INTENT(IN)     :: KFL ! number of fields
INTEGER(KIND=JPIM), INTENT(INOUT)  :: KPPGFL ! cumulated number of fields
INTEGER(KIND=JPIM), INTENT(OUT)    :: KMPP ! pointer

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('ENDPOS:SUPTRPPGFL_ENDPOS',0,ZHOOK_HANDLE)

KMPP=KPPGFL+1
KPPGFL=KPPGFL+KFL

IF (LHOOK) CALL DR_HOOK('ENDPOS:SUPTRPPGFL_ENDPOS',1,ZHOOK_HANDLE)
END SUBROUTINE SUPTRPPGFL_ENDPOS

!     ------------------------------------------------------------------
END SUBROUTINE ENDPOS
