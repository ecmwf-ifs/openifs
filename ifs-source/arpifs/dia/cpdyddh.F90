! (C) Copyright 1989- Meteo-France.

!OCL  NOEVAL
SUBROUTINE CPDYDDH(YDVAB,YDCVER,YDGMV,&
 !     INPUT .
 & YDML_GCONF,YDDPHY,YDML_DIAG,YDPHY,YDDYNA,KPROMA,KPROMNH,KSTART,KPROF,KFLEV,&
 & YDGSGEOM,YDCSGEOM,PXYB,&
 & PRCP, PAPHI, PAPHIF, PAPHIFL, PAPHIFM,&
 & PKENE, PRTL, PRTM, PORL, PORM,&
 & PSPL, PSPM, PAPRSF, PAPRS,&
 & PNHPREF,PNHPREH,PQCHAL,PQCHAM,&
 & PCTY, POMEGA, PFPLCL, PFPLCN, PFPLSL, PFPLSN,&
 & PQVS, PTSFC, PRH,&
 & PNEB,PQL,PQI,PQR,PQS,PEXT,PCOVPTOT,PGFL,PGMV,&
 & PGMVTNDSI,PGMVTNDHD,PGFLTNDHD,PATND,&
 !     OUTPUT .
 & PUZ, PVM, PDHCV, PENTRA, PENTRV,YDDDH)
!---------------------------------------------------------------------------
!**** *CPDYDDH*  - CALCUL DES VARIABLES ET FLUX/TENDANCES DYNAMIQUES
!                  DIAGNOSTICS PAR DOMAINES HORIZONTAUX (DDH)

!                  COMPUTATION OF VARIABLES AND FLUX/DYNAMICAL TENDENCIES
!                  HORIZONTAL DOMAINS DIAGNOSTICS ("DDH" PACKAGE)
!---------------------------------------------------------------------------
!     BUT/ PURPOSE:
!     -------------

!     Fait les calculs dynamiques pour   |  Does dynamical calculations for
!     les diagnostics "DDH".             |  the DDH package.
!     Principalement, fait le calcul     |  Does principally the computation
!     des variables DDH atmospheriques   |  of DDH atmospheric variables at
!     a l'instant courant.               |  the current instant.

!     Pour une variable X, on part de    |  For a variable X, one starts from
!     la forme flux de l'equation        |  the flux form of its temporal
!     d'evolution temporelle de X, qui   |  evolution equation, which evaluates
!     consiste a evaluer le membre de    |  the RHS of the equation of
!     droite de l'equation de            |
!      " (1/g) d ( (d prehyd / d eta) X) / dt ".

!     Ce membre de droite fait           |  This RHS involves not only some
!     intervenir, outre des              |  physical contributions (not
!     contributions physiques (pas       |  computed here but in CPPHDDH),
!     traitees ici mais dans CPPHDDH),   |  but also the following dynamical
!     les contributions dynamiques       |  contributions:
!     suivantes:                         |
!     - les divergences des flux         |  - horizontal fluxes divergences.
!       horizontaux.                     |
!     - les divergences des flux         |  - vertical fluxes divergences
!       verticaux (ce qui est stocke ici |    (what is stored are actually
!       ce sont les flux verticaux).     |    vertical fluxes).
!     - les tendances adiabatiques       |  - adiabatic Lagrangian tendencies
!       lagrangiennes "DX/Dt"            |    "DX/Dt".
!     - les termes lies a "delta m=1"    |  - "delta m=1" additional terms
!       (la vitesse verticale "EVEL" est |    (vertical velocity "EVEL"
!       supposee avoir ete modifiee par  |    is assumed to have been
!       la routine CPMVPPS).             |    modified in routine CPMVVPS).
!     Il faut remarquer que ces          |  It has to be noticed that these
!     diagnostics sont pleinement        |  diagnostics are fully consistent
!     coherents avec une formulation     |  with an Eulerian discretisation
!     eulerienne des equations avec      |  of equations with an explicit
!     traitement explicite pour le terme |  treatment of Coriolis term.
!     de Coriolis. Les divergences de    |  Flux divergences related to
!     flux liees a la diffusion          |  horizontal diffusion and semi-
!     horizontale et a la correction de  |  implicit scheme corrections are
!     semi-implicite ne sont pas         |  not available here.
!     disponibles ici.                   |

!     Les variables X traitees ici:      |  Variables diagnosed are:
!      * X = 1 (cf. eq de continuite)    |  * X = 1 (cf. continuity eqn)
!      * X = q                           |  * X = q
!      * X = vec(V), traitee dans un     |  * X = vec(V), in a local repere
!        repere local lie aux lat-lon    |    defined by the true geographical
!        geographiques.                  |    lat-lon coordinates.
!      * X = 0.5*vec(V)*vec(V)           |  * X = 0.5*vec(V)*vec(V)
!        (energie cinetique)             |    (kinetic energy)
!      * X = cp*T                        |  * X = cp*T
!      * X = cp*T +gz +0.5*vec(V)*vec(V) |  * X = cp*T +gz +0.5*vec(V)*vec(V)
!        (enthalpie)                     |    (enthalpy)
!      * X = ql (eau liquide)            |  * X = ql (liquid water)
!      * X = qi (glace)                  |  * X = qi (ice)
!      * X = MT (moment cinetique)       |  * X = MT (kinetic momentum)
!               MT = vec(r) wedge ( vec(Omega) wedge vec(r) + vec(V) )
!        traitee dans un repere fixe     |    in a fixed equatorial repere.
!        equatorial.                     |
!        (vec(r) = r vec(k))             |    (vec(r) = r vec(k))
!      * X = S (entropie)                |  * X = S (entropy)
!     Les tendances lagrangiennes        |  Lagrangian adiabatic tendencies
!     adiabatiques "(DX/Dt)_adiab"       |  "(DX/Dt)_adiab" are assumed to be
!     sont supposees etre nulles pour    |  zero for the following variables:
!     les variables suivantes:           |
!      * X = 1 (cf. eq de continuite)    |  * X = 1 (cf. continuity eqn)
!      * X = q                           |  * X = q
!      * X = cp*T +gz +0.5*vec(V)*vec(V) |  * X = cp*T +gz +0.5*vec(V)*vec(V)
!        (enthalpie)                     |    (enthalpy)
!      * X = ql (eau liquide)            |  * X = ql (liquid water)
!      * X = qi (glace)                  |  * X = qi (ice)
!      * X = qr (pluie)                  |  * X = qr (rain)
!      * X = qs (neige)                  |  * X = qs (snow)
!      * X = S (entropie)                |  * X = S (entropy)
!     Pour l'energie cinetique,          |  For kinetic energy,
!     "(D (0.5*vec(V)*vec(V))/Dt)_adiab" |  "(D (0.5*vec(V)*vec(V))/Dt)_adiab"
!     est calcule par:                   |  is computed as:
!     "vec(V)*(D vec(V))/Dt)_adiab"      |  "vec(V)*(D vec(V))/Dt)_adiab"
!     Pour l'equation du vent et de      |  For wind and kinetic energy
!     l'energie cinetique la friction    |  equations, Rayleigh friction is
!     de Rayleigh est ici ignoree.       |  omitted here.

!     Surveiller le positionnement des   |  Check fields positions with the
!     champs au moyen de l'indicateur    |  index integer "IDHCV".
!     "IDHCV".

!     Remarques supplementaires:         |  Additional remarks:

!     Le codage des termes sous          |  Code for terms under "delta m=1"
!     "delta m=1" semble douteux, le RHS |  seem suspicious, the RHS of
!     de l'equation de membre de gauche  |  equation, the LHS of which is
!             " (1/g) d ( (d prehyd / d eta) X) / dt ".
!     fait apparaitre "- X d Fp / d eta".|  yields the term "- X d Fp / d eta".
!     La forme discretisee devrait donc  |  Its discretised form then sould
!     faire apparaitre 0 si NDPSFI=0 et  |  yield "0" if NDPSFI=0 and
!     "-X[l] (Fp[lbar]-Fp[lbar-1])"      |  "-X[l] (Fp[lbar]-Fp[lbar-1])"
!     si NDPSFI=1.                       |  if NDPSFI=1.
!     De plus ce code ne semble pas      |  Furthermore this code does not seem
!     convenablement remis a jour vis    |  to have been properly updated
!     a vis des modifications            |  relatively to the modifications
!     de la formulation "delta m=1"      |  of the "delta m=1" formulation
!     introduites dans cy29t3.           |  introduced in cy29t3.

!     Modele non hydrostatique:          |  Non hydrostatic model:
!     le code a ete mis a jour pour les  |  the code has been updated for
!     variables presentes dans un modele |  the variables present in a
!     hydrostatique, mais les            |  hydrostatic model, but diagnostics
!     diagnostics pour les variables     |  relative to pressure departure
!     liees a "pre-prehyd" et la         |  and vertical divergence variables
!     divergence verticale n'ont pas ete |  have not yet been introduced.
!     encore introduits.                 |

!     Notations:                         |  Notations:
!     * prehyd: pression hydrostatique.  |  * prehyd: hydrostatic pressure.
!     * prehyds: pression hydrostatique  |  * prehyds: surface hydrostatic
!       de surface.                      |    pressure.
!     * pre: pression totale.            |  * pre: total pressure.
!     * grad: operateur gradient         |  * grad: is the
!       horizontal                       |    "horizontal gradient" operator:
!                  (grad X = vnabla X = M vnabla' X).
!     * (theta,lambda) sont les          |  * (theta,lambda) are the
!       coordonnees geographiques        |    geographical lat-lon coordinates
!       lat-lon du point.                |    of the grid-point.

!     ARGUMENTS D ENTREE / INPUT
!     --------------------------

!      KPROMA         : horizontal dimension.
!      KPROMNH        : horizontal dimension for arrays used only in NH model.
!      KSTART         : start of work.
!      KPROF          : working length.
!      KFLEV          : number of levels.

!      YDGSGEOM       : structure for geographical sphere horizontal geometry.
!      YDCSGEOM       : structure for computational sphere horizontal geometry.
!      PXYB           : contains pressure depth, "delta", "alpha".
!      PRCP           : contains "cp", "R" and "Kap=R/Cp".
!      PAPHI          : geopotential height "gz" at half levels.
!      PAPHIF         : geopotential height "gz" at full levels.
!      PAPHIFL        : zonal component of " grad (gz) " at full levels.
!      PAPHIFM        : meridian component of " grad (gz) " at full levels.
!      PKENE          : kinetic energy at full levels.
!      PRTL, PRTM     : components of "grad RT" at full levels.
!      PORL, PORM     : components of "grad Phi_s". 
!      PSPL, PSPM     : components of "grad prehyds" of surface hydrostatic pressure.
!      PAPRS          : hydrostatic pressure at half levels.
!      PAPRSF         : hydrostatic pressure at full levels.
!      PNHPREF        : "pre" at full levels (time t).
!      PNHPREH        : "pre" at half levels (time t).
!      PQCHAL         : zonal comp grad(log(pre/prehyd)).
!      PQCHAM         : merid comp grad(log(pre/prehyd)).
!      PCTY           : contains vertical velocities, vertical integral of divergence.
!      POMEGA         : "omega" at full levels, including the "lrubc" and "delta m=1" effects.
!      PFPLCL         : convective liquid rainfall flux, at half levels.
!      PFPLCN         : convective solid rainfall flux, at half levels.
!      PFPLSL         : stratiform liquid rainfall flux, at half levels.
!      PFPLSN         : stratiform solid rainfall flux, at half levels.
!      PQVS           : surface vapour specific humidity "q[surf]".
!      PTSFC          : surface temperature "T[surf]".
!      PRH            : relative humidity "RH" at full levels.
!      PNEB           : cloud fraction "q_a" at full levels.
!      PQL            : liquid water "q_l" at full levels.
!      PQI            : ice water "q_i" at full levels.
!      PQR            : rain "q_r" at full levels.
!      PQS            : snow "q_s" at full levels.
!      PEXT           : extra GFL at full levels.
!      PCOVPTOT       : precip fraction 
!      PGFL           : GFL variables (note that PQL to PEXT may differ from PGFL content).
!      PGMV           : GMV variables.
!      PGMVTNDSI      : tendencies of the semi-implicit scheme for the 
!                       horizontal wind components and for temperature
!      PGMVTNDHD      : tendencies of the horizontal diffusion scheme for the 
!                       horizontal wind components and for temperature
!      PGFLTNDHD      : tendencies of the horizontal diffusion scheme for the 
!                       specific humidity spectrally treated
!      PATND          : adiabatic Lagrangian tendencies.


!     SORTIES / OUTPUT
!     ----------------
!       PUZ, PVM      : horizontal wind components in a geographical local
!                       repere (i.e. the unit vector vec(j) is directed towards
!                       the true North pole).
!       PDHCV         : result array containing flux divergences and fluxes.
!                       Entering this routine this array is assumed to be
!                       pre-initialised with zeros.
!       PENTRA        : specific entropy of dry air "Sd", at full levels.
!       PENTRV        : specific entropy of water vapour "Sv", at full levels.

!     ORIGINAL
!     --------
!     91-03-04, Alain Joly.

!     MODIFICATIONS
!     -------------
!   N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!   17-09-2008   O. Riviere New data flux for diag. if LFLEXDIA=.TRUE.
!   K. Yessad (Dec 2008): remove dummy CDLOCK
!   K. Yessad (March 2009): correct false comments for LRWSDLG=T
!   O.Riviere (Sept 09) : add Hail for DDH under LFLEXDIA
!   Y. Seity (Jan 11) bf for AROME runs without hail
!   K. Yessad (Nov 2009): DT/Dt now pre-computed in CPG_GP.
!   K. Yessad (Jan 2011): introduce INTDYN_MOD structures.
!   K. Yessad (Nov 2011): various modifications.
!   M. Ahlgrimm  31-Oct-2011 add rain, snow and PEXTRA to DDH output
!   M. Ahlgrimm Apr 2014: Add lake variables and precip fraction to DDH output
!   T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!   K. Yessad (July 2014): Move some variables.
!   R. El Khatib 03-Sep-2014 vectorizations
!   Y. Seity 25-Sep-2014 Bf if LAROME (missing fields)
!   Y. Seity 11-Oct-2016 Bf if LAROME (VNT already written under apl_arome)
!   K. Yessad (Feb 2018): remove deep-layer formulations.
!   M. Hrastinski Sep-2019: TKE and TTE terms for ALARO DDH (tendencies,
!                           horiz. and vert. fluxes, upper and bottom BC)
!     ------------------------------------------------------------------

USE MODEL_DIAGNOSTICS_MOD  , ONLY : MODEL_DIAGNOSTICS_TYPE
USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE YOMPHY                 , ONLY : TPHY, JPHYARO
USE YOMVERT                , ONLY : TVAB
USE YOMGMV                 , ONLY : TGMV
USE YOMDYNA                , ONLY : TDYNA
USE PARKIND1               , ONLY : JPIM, JPRB
USE YOMHOOK                , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST                 , ONLY : ROMEGA, RA, RG, RD, RV, RCPD, RCPV, RCW, RCS, RTT, RATM  
USE YOMRIP0                , ONLY : NSSSSS
USE YOMLUN                 , ONLY : NULOUT
USE YOMCSGEOM              , ONLY : TCSGEOM
USE YOMGSGEOM              , ONLY : TGSGEOM
USE YOMCT0                 , ONLY : LAROME, NUNDEFLD
USE YOMCVER           , ONLY : TCVER
USE YOMDPHY                , ONLY : TDPHY
USE INTDYN_MOD             , ONLY : YYTCTY0, YYTRCP0, YYTXYB0
USE DDH_MIX                , ONLY : ADD_FIELD_3D, NEW_ADD_FIELD_3D, TYP_DDH

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TVAB)                   , INTENT(IN)    :: YDVAB
TYPE(TCVER)                  , INTENT(IN)    :: YDCVER
TYPE(TGMV)                   , INTENT(INOUT) :: YDGMV
TYPE(TDPHY)                  , INTENT(IN)    :: YDDPHY
TYPE(MODEL_DIAGNOSTICS_TYPE) , INTENT(INOUT) :: YDML_DIAG
TYPE(MODEL_GENERAL_CONF_TYPE), INTENT(IN)    :: YDML_GCONF
TYPE(TPHY)                   , INTENT(IN)    :: YDPHY
TYPE(TDYNA)                  , INTENT(IN)    :: YDDYNA
INTEGER(KIND=JPIM)           , INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM)           , INTENT(IN)    :: KPROMNH
INTEGER(KIND=JPIM)           , INTENT(IN)    :: KFLEV 
INTEGER(KIND=JPIM)           , INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM)           , INTENT(IN)    :: KPROF 
TYPE(TGSGEOM)                , INTENT(IN)    :: YDGSGEOM
TYPE(TCSGEOM)                , INTENT(IN)    :: YDCSGEOM
REAL(KIND=JPRB)              , INTENT(IN)    :: PXYB(KPROMA,KFLEV,YYTXYB0%NDIM) 
REAL(KIND=JPRB)              , INTENT(IN)    :: PRCP(KPROMA,KFLEV,YYTRCP0%NDIM) 
REAL(KIND=JPRB)              , INTENT(IN)    :: PAPHI(KPROMA,0:KFLEV) 
REAL(KIND=JPRB)              , INTENT(IN)    :: PAPHIF(KPROMA,KFLEV) 
REAL(KIND=JPRB)              , INTENT(IN)    :: PAPHIFL(KPROMA,KFLEV) 
REAL(KIND=JPRB)              , INTENT(IN)    :: PAPHIFM(KPROMA,KFLEV) 
REAL(KIND=JPRB)              , INTENT(IN)    :: PKENE(KPROMA,KFLEV) 
REAL(KIND=JPRB)              , INTENT(IN)    :: PRTL(KPROMA,KFLEV) 
REAL(KIND=JPRB)              , INTENT(IN)    :: PRTM(KPROMA,KFLEV) 
REAL(KIND=JPRB)              , INTENT(IN)    :: PORL(KPROMA) 
REAL(KIND=JPRB)              , INTENT(IN)    :: PORM(KPROMA) 
REAL(KIND=JPRB)              , INTENT(IN)    :: PSPL(KPROMA) 
REAL(KIND=JPRB)              , INTENT(IN)    :: PSPM(KPROMA) 
REAL(KIND=JPRB)              , INTENT(IN)    :: PAPRSF(KPROMA,KFLEV) 
REAL(KIND=JPRB)              , INTENT(IN)    :: PAPRS(KPROMA,0:KFLEV) 
REAL(KIND=JPRB)              , INTENT(IN)    :: PNHPREF(KPROMNH,KFLEV)
REAL(KIND=JPRB)              , INTENT(IN)    :: PNHPREH(KPROMNH,0:KFLEV)
REAL(KIND=JPRB)              , INTENT(IN)    :: PQCHAL(KPROMNH,KFLEV)
REAL(KIND=JPRB)              , INTENT(IN)    :: PQCHAM(KPROMNH,KFLEV)
REAL(KIND=JPRB)              , INTENT(IN)    :: PCTY(KPROMA,0:KFLEV,YYTCTY0%NDIM)
REAL(KIND=JPRB)              , INTENT(IN)    :: POMEGA(KPROMA,KFLEV) 
REAL(KIND=JPRB)              , INTENT(IN)    :: PFPLCL(KPROMA,0:KFLEV) 
REAL(KIND=JPRB)              , INTENT(IN)    :: PFPLCN(KPROMA,0:KFLEV) 
REAL(KIND=JPRB)              , INTENT(IN)    :: PFPLSL(KPROMA,0:KFLEV) 
REAL(KIND=JPRB)              , INTENT(IN)    :: PFPLSN(KPROMA,0:KFLEV) 
REAL(KIND=JPRB)              , INTENT(IN)    :: PQVS(KPROMA) 
REAL(KIND=JPRB)              , INTENT(IN)    :: PTSFC(KPROMA) 
REAL(KIND=JPRB)              , INTENT(IN)    :: PRH(KPROMA,KFLEV) 
REAL(KIND=JPRB)              , INTENT(IN)    :: PNEB(KPROMA,KFLEV) 
REAL(KIND=JPRB)              , INTENT(IN)    :: PQL(KPROMA,KFLEV) 
REAL(KIND=JPRB)              , INTENT(IN)    :: PQI(KPROMA,KFLEV) 
REAL(KIND=JPRB)              , INTENT(IN)    :: PQR(KPROMA,KFLEV)
REAL(KIND=JPRB)              , INTENT(IN)    :: PQS(KPROMA,KFLEV)
REAL(KIND=JPRB)              , INTENT(IN)    :: PEXT(KPROMA,KFLEV,YDDPHY%NVEXTR)
REAL(KIND=JPRB)              , INTENT(IN)    :: PCOVPTOT(KPROMA,KFLEV,1) 
REAL(KIND=JPRB)              , INTENT(IN)    :: PGFL(KPROMA,KFLEV,YDML_GCONF%YGFL%NDIM)
REAL(KIND=JPRB)              , INTENT(IN)    :: PGMV(KPROMA,KFLEV,YDGMV%NDIMGMV)
REAL(KIND=JPRB)              , INTENT(IN)    :: PGMVTNDSI(KPROMA,KFLEV,2+YDML_GCONF%YRDIMF%NFTHER)
REAL(KIND=JPRB)              , INTENT(IN)    :: PGMVTNDHD(KPROMA,KFLEV,2+YDML_GCONF%YRDIMF%NFTHER)
REAL(KIND=JPRB)              , INTENT(IN)    :: PGFLTNDHD(KPROMA,KFLEV,1)
REAL(KIND=JPRB)              , INTENT(IN)    :: PATND(KPROMA,KFLEV,YDDYNA%YYTTND%NDIM)
REAL(KIND=JPRB)              , INTENT(INOUT) :: PUZ(KPROMA,KFLEV) 
REAL(KIND=JPRB)              , INTENT(INOUT) :: PVM(KPROMA,KFLEV) 
REAL(KIND=JPRB)              , INTENT(INOUT) :: PDHCV(KPROMA,0:KFLEV,YDML_DIAG%YRMDDH%NDHCVSU) 
REAL(KIND=JPRB)              , INTENT(INOUT) :: PENTRA(KPROMA,KFLEV) 
REAL(KIND=JPRB)              , INTENT(INOUT) :: PENTRV(KPROMA,KFLEV)
TYPE(TYP_DDH)                , INTENT(INOUT) :: YDDDH 
!     ------------------------------------------------------------------
REAL(KIND=JPRB) :: ZT_EVEL(KPROMA,0:KFLEV) 
REAL(KIND=JPRB) :: ZT_PSDIV(KPROMA,0:KFLEV) 
REAL(KIND=JPRB) :: ZT_DIVDP(KPROMA,KFLEV) 
REAL(KIND=JPRB) :: ZT_DELP(KPROMA,KFLEV) 
REAL(KIND=JPRB) :: ZT_RDELP(KPROMA,KFLEV) 

REAL(KIND=JPRB) :: ZCPT(KPROMA,KFLEV)
REAL(KIND=JPRB) :: ZGKEU(KPROMA,KFLEV), ZGKEV(KPROMA,KFLEV)
REAL(KIND=JPRB) :: ZSLON(KPROMA), ZCLON(KPROMA)
REAL(KIND=JPRB) :: ZM1(KPROMA,KFLEV), ZM2(KPROMA,KFLEV), ZM3(KPROMA,KFLEV)
REAL(KIND=JPRB) :: ZRTZON(KPROMA,KFLEV), ZRTMER(KPROMA,KFLEV)
REAL(KIND=JPRB) :: ZPHILZ(KPROMA), ZPHILM(KPROMA)
REAL(KIND=JPRB) :: ZFPSUM(KPROMA,0:KFLEV)
REAL(KIND=JPRB) :: ZFPTOT(KPROMA,0:KFLEV)
REAL(KIND=JPRB) :: ZDELLT(KPROMA,0:KFLEV)
REAL(KIND=JPRB) :: ZENTR(KPROMA,KFLEV)
REAL(KIND=JPRB) :: ZENTRL(KPROMA,0:KFLEV), ZENTRS(KPROMA,0:KFLEV)
REAL(KIND=JPRB) :: ZDELPREV(KPROMA,KFLEV)
REAL(KIND=JPRB) :: ZRHSK(KPROMA,KFLEV)
REAL(KIND=JPRB) :: ZRHSME(KPROMA)
REAL(KIND=JPRB) :: ZRHSU(KPROMA,KFLEV)
REAL(KIND=JPRB) :: ZRHSV(KPROMA,KFLEV)
REAL(KIND=JPRB) :: ZRHSZO(KPROMA)
REAL(KIND=JPRB) :: ZEVELH(KPROMA,0:KFLEV)
REAL(KIND=JPRB) :: ZAPRSF(KPROMA,KFLEV)
REAL(KIND=JPRB) :: ZAPRS(KPROMA,0:KFLEV)
REAL(KIND=JPRB) :: ZQCHAL(KPROMA,KFLEV)
REAL(KIND=JPRB) :: ZQCHAM(KPROMA,KFLEV)
REAL(KIND=JPRB) :: Z_CP_TIMES_CONV(KPROMA,KFLEV)
REAL(KIND=JPRB) :: ZTMPAF(KPROMA,KFLEV) 

INTEGER(KIND=JPIM) :: IDHCV, JLEV, JROF, JEXT
INTEGER(KIND=JPIM) :: IPTR_S, IPTR_SL, IPTR_SM, IPTR_R, IPTR_RL, IPTR_RM, IPTR_Q, IPTR_QL, IPTR_QM
INTEGER(KIND=JPIM) :: IPTR_G, IPTR_H, IPTR_LL, IPTR_LM, IPTR_IL, IPTR_IM
INTEGER(KIND=JPIM) :: IPTR_V, IPTR_U, IPTR_VL, IPTR_UL, IPTR_T, IPTR_TL, IPTR_TM
INTEGER(KIND=JPIM) :: IPTR_CP, IPTR_RD
INTEGER(KIND=JPIM) :: IPTR_TKE,IPTR_TTE,IPTR_TKEM,IPTR_TKEL,IPTR_TTEM,IPTR_TTEL

REAL(KIND=JPRB) :: ZADPHI, ZADVLNT, ZADVP,&
 & ZADVQL, ZADVQI, ZADVQR, ZADVQS, ZADVQV,&
 & ZADVT, ZCVMCD, &
 & ZDPREM, ZDPREZ, ZDPSFI, ZENTRA,&
 & ZENTRA0, ZENTRA00, ZENTRAS, ZENTRL0, ZENTRL00,&
 & ZENTRS0, ZENTRS00, ZENTRV, ZENTRV0, ZENTRV00,&
 & ZENTRVS, ZLTBOT, ZLTHL, ZLTTOP, ZORMER, ZORZON,&
 & ZPFIS2, ZPREV, ZPREVS, ZQTHRES,&
 & ZUSA, ZUS2RG, ZUSRG, ZADVTKE, ZADVTTE

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CPDYDDH',0,ZHOOK_HANDLE)
ASSOCIATE(NDIM=>YDML_GCONF%YGFL%NDIM, YG=>YDML_GCONF%YGFL%YG, YH=>YDML_GCONF%YGFL%YH, YI=>YDML_GCONF%YGFL%YI, &
 &  YL=>YDML_GCONF%YGFL%YL, &
 & YQ=>YDML_GCONF%YGFL%YQ, YR=>YDML_GCONF%YGFL%YR, YS=>YDML_GCONF%YGFL%YS, &
 & YTKE=>YDML_GCONF%YGFL%YTKE, YTTE=>YDML_GCONF%YGFL%YTTE, &
 & NFTHER=>YDML_GCONF%YRDIMF%NFTHER, &
 & NVEXTR=>YDDPHY%NVEXTR, NVECOUT=>YDDPHY%NVECOUT,&
 & NDIMGMV=>YDGMV%NDIMGMV, YT0=>YDGMV%YT0, &
 & LFLEXDIA=>YDML_DIAG%YRLDDH%LFLEXDIA, LHDENT=>YDML_DIAG%YRLDDH%LHDENT, LHDHKS=>YDML_DIAG%YRLDDH%LHDHKS, &
 & LHDLIST=>YDML_DIAG%YRLDDH%LHDLIST, LHDMCI=>YDML_DIAG%YRLDDH%LHDMCI, LRHDDDH=>YDML_DIAG%YRLDDH%LRHDDDH, &
 & LRSIDDH=>YDML_DIAG%YRLDDH%LRSIDDH, LRSLDDH=>YDML_DIAG%YRLDDH%LRSLDDH, LDDH_OMP=>YDML_DIAG%YRLDDH%LDDH_OMP, &
 & MHDDDH_NHX=>YDML_DIAG%YRMDDH%MHDDDH_NHX, MHDDDH_PD=>YDML_DIAG%YRMDDH%MHDDDH_PD, &
 & MHDDDH_Q=>YDML_DIAG%YRMDDH%MHDDDH_Q, MHDDDH_T=>YDML_DIAG%YRMDDH%MHDDDH_T, &
 & MHDDDH_U=>YDML_DIAG%YRMDDH%MHDDDH_U, MHDDDH_V=>YDML_DIAG%YRMDDH%MHDDDH_V, &
 & MHDDDH_VD=>YDML_DIAG%YRMDDH%MHDDDH_VD, MSIDDH_PD0=>YDML_DIAG%YRMDDH%MSIDDH_PD0, &
 & MSIDDH_PD1=>YDML_DIAG%YRMDDH%MSIDDH_PD1, MSIDDH_T0=>YDML_DIAG%YRMDDH%MSIDDH_T0, &
 & MSIDDH_T1=>YDML_DIAG%YRMDDH%MSIDDH_T1, MSIDDH_U0=>YDML_DIAG%YRMDDH%MSIDDH_U0, &
 & MSIDDH_U1=>YDML_DIAG%YRMDDH%MSIDDH_U1, MSIDDH_V0=>YDML_DIAG%YRMDDH%MSIDDH_V0, &
 & MSIDDH_V1=>YDML_DIAG%YRMDDH%MSIDDH_V1, MSIDDH_VD0=>YDML_DIAG%YRMDDH%MSIDDH_VD0, &
 & MSIDDH_VD1=>YDML_DIAG%YRMDDH%MSIDDH_VD1, NDHAHKD=>YDML_DIAG%YRMDDH%NDHAHKD, &
 & NDHAMCD=>YDML_DIAG%YRMDDH%NDHAMCD, NDHAVD=>YDML_DIAG%YRMDDH%NDHAVD, NDHBHKD=>YDML_DIAG%YRMDDH%NDHBHKD, &
 & NDHBMCD=>YDML_DIAG%YRMDDH%NDHBMCD, NDHCVSU=>YDML_DIAG%YRMDDH%NDHCVSU, NDHVHK=>YDML_DIAG%YRMDDH%NDHVHK, &
 & NDHVMC=>YDML_DIAG%YRMDDH%NDHVMC, NDHVV=>YDML_DIAG%YRMDDH%NDHVV, &
 & RSTATI=>YDML_GCONF%YRRIP%RSTATI, &
 & LHDPAS=>YDML_DIAG%YRSDDH%LHDPAS, LHDQLN=>YDML_DIAG%YRSDDH%LHDQLN, NHDPASVA=>YDML_DIAG%YRSDDH%NHDPASVA, &
 & NHDQLNVA=>YDML_DIAG%YRSDDH%NHDQLNVA, L3MT=>YDPHY%L3MT, LSTRAPRO=>YDPHY%LSTRAPRO, &
 & NDPSFI=>YDPHY%NDPSFI, NPHY=>YDPHY%NPHY, LPTKE=>YDPHY%LPTKE, &
 & LCOEFK_PTTE=>YDPHY%LCOEFK_PTTE)
!     ------------------------------------------------------------------

!     CHECK RELIABILITY OF INPUT ARGUMENTS.

! Avoid unneccessary complex addressing, which can break the vectorization :
IPTR_Q=YQ%MP
IPTR_QL=YQ%MPL
IPTR_QM=YQ%MPM
IPTR_S=YS%MP
IPTR_SL=YS%MPL
IPTR_SM=YS%MPM
IPTR_R=YR%MP
IPTR_RL=YR%MPL
IPTR_RM=YR%MPM
IPTR_G=YG%MP
IPTR_H=YH%MP
IPTR_TKE=YTKE%MP
IPTR_TKEL=YTKE%MPL
IPTR_TKEM=YTKE%MPM
IPTR_TTE=YTTE%MP
IPTR_TTEL=YTTE%MPL
IPTR_TTEM=YTTE%MPM
IPTR_LL=YL%MPL
IPTR_LM=YL%MPM
IPTR_IL=YI%MPL
IPTR_IM=YI%MPM
IPTR_V=YT0%MV
IPTR_U=YT0%MU
IPTR_VL=YT0%MVL
IPTR_UL=YT0%MUL
IPTR_T=YT0%MT
IPTR_TL=YT0%MTL
IPTR_TM=YT0%MTM
IPTR_CP=YYTRCP0%M_CP
IPTR_RD=YYTRCP0%M_R

!     ------------------------------------------------------------------
!*    0.- INITIALISATIONS DE SECURITE / SAFETY SET-UP.
!         AFIN DE CONTOURNER LE MESSAGE DE WARNING PRODUIT PAR
!         LE LOGICIEL DE CONTROLE DES INITIALISATIONS DE VARIABLES
!         EN CAS DE VARIABLES INITIALISEES A L'INTERIEUR DE
!         TESTS IF-ENDIF, ON PORTE UNE VALEUR QUELCONQUE
!         DANS LESDITES VARIABLES.

ZENTRA0=0.0_JPRB
ZENTRV0=0.0_JPRB
ZQTHRES=0.0_JPRB
IDHCV=0

ZT_EVEL(KSTART:KPROF,0:KFLEV)=PCTY(KSTART:KPROF,0:KFLEV,YYTCTY0%M_EVEL)
ZT_PSDIV(KSTART:KPROF,0:KFLEV)=PCTY(KSTART:KPROF,0:KFLEV,YYTCTY0%M_PSDIV)
ZT_DIVDP(KSTART:KPROF,1:KFLEV)=PCTY(KSTART:KPROF,1:KFLEV,YYTCTY0%M_DIVDP)

ZT_DELP(KSTART:KPROF,1:KFLEV)=PXYB(KSTART:KPROF,1:KFLEV,YYTXYB0%M_DELP)
ZT_RDELP(KSTART:KPROF,1:KFLEV)=PXYB(KSTART:KPROF,1:KFLEV,YYTXYB0%M_RDELP)

! * Calcul des epaisseurs de pression hydrostatique virtuelle sur les couches:
!   Computation of virtual hydrostatic pressure depths at full levels:
!   "Delta prehydv = Delta prehyd".

ZDELPREV(KSTART:KPROF,1:KFLEV)=ZT_DELP(KSTART:KPROF,1:KFLEV)

! * Variable intermediaire contenant / intermediate variable containing:
!    "(etapt d prehyd/d eta)_[lbar]".
IF (YDCVER%LVERTFE) THEN
  ! * Top and surface:
  ZEVELH(KSTART:KPROF,0)=0.0_JPRB     ! lrubc=t, ndpsfi=1 not yet coded
  ZEVELH(KSTART:KPROF,KFLEV)=0.0_JPRB ! lrubc=t, ndpsfi=1 not yet coded
  ! * half levels 1 to L-1:
  DO JLEV=1,KFLEV-1
    ZEVELH(KSTART:KPROF,JLEV)=&
     & 0.5_JPRB*(ZT_EVEL(KSTART:KPROF,JLEV)+ZT_EVEL(KSTART:KPROF,JLEV+1))  
  ENDDO
ELSE
  ZEVELH(KSTART:KPROF,0:KFLEV)=ZT_EVEL(KSTART:KPROF,0:KFLEV)
ENDIF

! * Total pressure or hydrostatic pressure according to LNHDYN:
IF (YDDYNA%LNHDYN) THEN
  ZAPRSF(KSTART:KPROF,1:KFLEV)=PNHPREF(KSTART:KPROF,1:KFLEV)
  ZAPRS(KSTART:KPROF,0:KFLEV)=PNHPREH(KSTART:KPROF,0:KFLEV)
  ZQCHAL(KSTART:KPROF,1:KFLEV)=PQCHAL(KSTART:KPROF,1:KFLEV)
  ZQCHAM(KSTART:KPROF,1:KFLEV)=PQCHAM(KSTART:KPROF,1:KFLEV)
ELSE
  ZAPRSF(KSTART:KPROF,1:KFLEV)=PAPRSF(KSTART:KPROF,1:KFLEV)
  ZAPRS(KSTART:KPROF,0:KFLEV)=PAPRS(KSTART:KPROF,0:KFLEV)
  ZQCHAL(KSTART:KPROF,1:KFLEV)=0.0_JPRB
  ZQCHAM(KSTART:KPROF,1:KFLEV)=0.0_JPRB
ENDIF

! * "cp*(conversion term)" (in the RHS of cp*T equation):
DO JLEV=1,KFLEV
  DO JROF=KSTART,KPROF
    Z_CP_TIMES_CONV(JROF,JLEV)=PRCP(JROF,JLEV,IPTR_CP)*PATND(JROF,JLEV,YDDYNA%YYTTND%M_TNDT)
  ENDDO
ENDDO

!     ------------------------------------------------------------------
!*    1.- CALCUL DES VARIABLES DYNAMIQUES / DYNAMICAL VARIABLES COMPUTATIONS.

ZUSA=1.0_JPRB/RA
ZUSRG = 1.0_JPRB/RG
ZCVMCD = RCPV - RCPD
ZDPSFI = REAL(NDPSFI,JPRB)

! * Calcul des composantes du vent geographique dans un repere local lie au
!   vrai pole nord.
!   Computation of geographical wind components in a local repere, the ordinate
!   of which is directed towards the true North pole.

IF ( LHDHKS .OR. LHDMCI .OR. LHDENT ) THEN
  DO JLEV = 1, KFLEV
    DO JROF=KSTART,KPROF
      PUZ(JROF,JLEV) =YDGSGEOM%GNORDM(JROF)*PGMV(JROF,JLEV,IPTR_U)-YDGSGEOM%GNORDL(JROF)*PGMV(JROF,JLEV,IPTR_V)
      PVM(JROF,JLEV) =YDGSGEOM%GNORDL(JROF)*PGMV(JROF,JLEV,IPTR_U)+YDGSGEOM%GNORDM(JROF)*PGMV(JROF,JLEV,IPTR_V)
    ENDDO
  ENDDO
ENDIF

! * Calcul du gradient horizontal de l'energie cinetique.
!   Computation of the horizontal gradient of kinetic energy.

IF ( LHDHKS ) THEN
  IF ( YDDYNA%LSLAG ) THEN
    DO JLEV = 1, KFLEV
      DO JROF=KSTART,KPROF
        ZGKEU(JROF,JLEV) = 0.0_JPRB
        ZGKEV(JROF,JLEV) = 0.0_JPRB
      ENDDO
    ENDDO
  ELSE
    DO JLEV = 1, KFLEV
      DO JROF=KSTART,KPROF
        ZGKEU(JROF,JLEV) = - PGMV(JROF,JLEV,IPTR_V)* PGMV(JROF,JLEV,YT0%MVOR)&
         & + PGMV(JROF,JLEV,IPTR_V)*PGMV(JROF,JLEV,IPTR_VL)&
         & + PGMV(JROF,JLEV,IPTR_U)*PGMV(JROF,JLEV,IPTR_UL)&
         & - PKENE(JROF,JLEV)*YDCSGEOM%RATATX(JROF)  
        ZGKEV(JROF,JLEV) = PGMV(JROF,JLEV,IPTR_V)* PGMV(JROF,JLEV,YT0%MDIV)&
         & - PGMV(JROF,JLEV,IPTR_V)*PGMV(JROF,JLEV,IPTR_UL)&
         & + PGMV(JROF,JLEV,IPTR_U)*PGMV(JROF,JLEV,IPTR_VL)&
         & + PKENE(JROF,JLEV)*YDCSGEOM%RATATH(JROF)  
      ENDDO
    ENDDO
  ENDIF
ENDIF

! * Calcul de "cp T" sur les couches.
!   Computation of "cp T" at full levels.

IF ( LHDHKS ) THEN
  DO JLEV = 1, KFLEV
    DO JROF=KSTART,KPROF
      ZCPT(JROF,JLEV) = PRCP(JROF,JLEV,IPTR_CP)*PGMV(JROF,JLEV,IPTR_T)
    ENDDO
  ENDDO
ENDIF

! * Calculs relatifs a l'entropie.
!    On calcule d'abord des entropies de reference (variables ZENTR..00
!    et ZENTR..0) puis on calcule:
!    - l'entropie de l'air sec (PENTRA).
!    - l'entropie de la vapeur d'eau (PENTRV).
!    - l'entropie de l'air humide (ZENTR).
!    - l'entropie de l'eau liquide (ZENTRL).
!    - l'entropie de la glace (ZENTRS).
!   Computations involving entropy.
!    One first computes some reference entropies ((variables ZENTR..00
!    and ZENTR..0), one then computes:
!    - dry air entropy (PENTRA).
!    - water vapour entropy (PENTRV).
!    - moist air entropy (ZENTR).
!    - liquid water entropy (ZENTRL).
!    - solid water entropy (ZENTRS).
!   Entropy = cp T + R log(pre), with "pre"=total pressure.

IF ( LHDENT ) THEN
  DO JLEV = 2, KFLEV-1
    DO JROF=KSTART,KPROF
      ZDELLT(JROF,JLEV) = ( PGMV(JROF,JLEV+1,IPTR_T)-PGMV(JROF,JLEV-1,IPTR_T) ) /&
       & (2.0_JPRB*PGMV(JROF,JLEV,IPTR_T))  
    ENDDO
  ENDDO
  DO JROF=KSTART,KPROF
    ZDELLT(JROF,1) = - 0.5_JPRB + PGMV(JROF,2,IPTR_T)/(2.0_JPRB*PGMV(JROF,1,IPTR_T))
    ZDELLT(JROF,KFLEV) = - (0.5_JPRB*PGMV(JROF,KFLEV-1,IPTR_T)-PTSFC(JROF)) /&
     & PGMV(JROF,KFLEV,IPTR_T) - 0.5_JPRB  
  ENDDO
! * THE REFERENCE SPECIFIC ENTROPIES
  ZENTRA00 = 6775._JPRB
  ZENTRV00 = 10320._JPRB
  ZENTRL00 = 3517._JPRB
  ZENTRS00 = 2296._JPRB
  ZENTRA0 = -RCPD*LOG(RTT)+RD*LOG(RATM)+ZENTRA00
  ZENTRV0 = -RCPV*LOG(RTT)+RV*LOG(RATM)+ZENTRV00
  ZENTRL0 = -RCW*LOG(RTT)+ZENTRL00
  ZENTRS0 = -RCS*LOG(RTT)+ZENTRS00
! * THE THRESHOLD VALUE AVOIDING NEGATIVE HUMIDITY
  ZQTHRES = 1.E-12_JPRB
! * CALCULATION OF THE SPECIFIC ENTROPIES
  DO JLEV = 1, KFLEV
    DO JROF=KSTART,KPROF
      ZPREV = MAX(ZQTHRES,PGFL(JROF,JLEV,IPTR_Q))*RV*ZAPRSF(JROF,JLEV)/PRCP(JROF,JLEV,IPTR_RD)
      ZENTRA = RCPD*LOG(PGMV(JROF,JLEV,IPTR_T))-&
       & RD*LOG(ZAPRSF(JROF,JLEV)-ZPREV)+ZENTRA0  
      ZENTRV = RCPV*LOG(PGMV(JROF,JLEV,IPTR_T))-RV*LOG(ZPREV)+ZENTRV0
      ZENTR(JROF,JLEV) = ZENTRA+(ZENTRV-ZENTRA)*PGFL(JROF,JLEV,IPTR_Q)
      PENTRA(JROF,JLEV) = ZENTRA
      PENTRV(JROF,JLEV) = ZENTRV
    ENDDO
  ENDDO
  DO JLEV=1,KFLEV-1
    DO JROF=KSTART,KPROF
      ZLTHL = LOG((PGMV(JROF,JLEV,IPTR_T)+PGMV(JROF,JLEV,IPTR_T))/2.0_JPRB)
      ZENTRL(JROF,JLEV) = RCW*ZLTHL + ZENTRL0
      ZENTRS(JROF,JLEV) = RCS*ZLTHL + ZENTRS0
    ENDDO
  ENDDO
  DO JROF=KSTART,KPROF
    ZLTTOP = LOG(PGMV(JROF,1,IPTR_T))
    ZLTBOT = LOG(PTSFC(JROF))
    ZENTRL(JROF,0) = RCW*ZLTTOP + ZENTRL0
    ZENTRL(JROF,KFLEV) = RCW*ZLTBOT + ZENTRL0
    ZENTRS(JROF,0) = RCS*ZLTBOT + ZENTRS0
    ZENTRS(JROF,KFLEV) = RCS*ZLTBOT + ZENTRS0
  ENDDO
! * TOTAL FLUX OF PRECIPITATIONS
  DO JLEV=0, KFLEV
    DO JROF=KSTART,KPROF
      ZFPTOT(JROF,JLEV) = PFPLCL(JROF,JLEV) + PFPLCN(JROF,JLEV) +&
       & PFPLSL(JROF,JLEV) + PFPLSN(JROF,JLEV)  
    ENDDO
  ENDDO
ENDIF

! * Calcul des trois composantes M1, M2 et M3 du moment angulaire,
!   puis des composantes de "grad (RT)" dans un repere local lie au
!   vrai pole nord.
!   Computation of the three components M1, M2 et M3 of angular moment,
!   then computation of the components of "grad (RT)" in a local repere,
!   the ordinate of which is directed towards the true North pole.

IF ( LHDMCI ) THEN
  DO JROF=KSTART,KPROF
    ZSLON(JROF) = SIN(YDGSGEOM%GELAM(JROF)+ROMEGA*(RSTATI-NSSSSS))
    ZCLON(JROF) = COS(YDGSGEOM%GELAM(JROF)+ROMEGA*(RSTATI-NSSSSS))
  ENDDO
  DO JLEV=1, KFLEV
    DO JROF=KSTART,KPROF
      ZM1(JROF,JLEV) = RA*( PVM(JROF,JLEV)*ZSLON(JROF)&
       & -(PUZ(JROF,JLEV)+RA*ROMEGA*YDGSGEOM%GSQM2(JROF))*&
       & YDGSGEOM%GEMU(JROF)*ZCLON(JROF) )  
      ZM2(JROF,JLEV) = RA*( -PVM(JROF,JLEV)*ZCLON(JROF)&
       & -(PUZ(JROF,JLEV)+RA*ROMEGA*YDGSGEOM%GSQM2(JROF))*&
       & YDGSGEOM%GEMU(JROF)*ZSLON(JROF) )  
      ZM3(JROF,JLEV) = RA*(PUZ(JROF,JLEV)&
       & +RA*ROMEGA*YDGSGEOM%GSQM2(JROF))*YDGSGEOM%GSQM2(JROF)  
      ZRTZON(JROF,JLEV) = YDGSGEOM%GNORDM(JROF)*PRTL(JROF,JLEV) -&
       & YDGSGEOM%GNORDL(JROF)*PRTM(JROF,JLEV)  
      ZRTMER(JROF,JLEV) = YDGSGEOM%GNORDL(JROF)*PRTL(JROF,JLEV) +&
       & YDGSGEOM%GNORDM(JROF)*PRTM(JROF,JLEV)  
    ENDDO
  ENDDO
ENDIF

! * Calcul de "(Delta prehyd) X" sur les couches, pour
!   X = 1; (-1/g); q; Ugeo; Vgeo; (U**2+V**2)/2; cp*T; g*z; RH; ql; qi; qa;
!   omega; M1; M2; M3; S.
!   Computation of "(Delta prehyd) X" at full levels,
!   for X = 1; (-1/g); q; Ugeo; Vgeo; (U**2+V**2)/2; cp*T; g*z; RH; ql; qi; qa;
!   omega; M1; M2; M3; S.

IF ( LHDHKS ) THEN
  IF(LFLEXDIA) THEN
    IF (LDDH_OMP) THEN
      CALL NEW_ADD_FIELD_3D(YDML_DIAG%YRMDDH,ZDELPREV,'VPP',YDDDH)
      ZTMPAF(KSTART:KPROF,1:KFLEV)=&
       & ZDELPREV(KSTART:KPROF,1:KFLEV)*PAPHIF(KSTART:KPROF,1:KFLEV)
      CALL NEW_ADD_FIELD_3D(YDML_DIAG%YRMDDH,ZTMPAF,'VEP',YDDDH)
      ZTMPAF(KSTART:KPROF,1:KFLEV)=&
       & ZDELPREV(KSTART:KPROF,1:KFLEV)*PRH(KSTART:KPROF,1:KFLEV)
      CALL NEW_ADD_FIELD_3D(YDML_DIAG%YRMDDH,ZTMPAF,'VHR',YDDDH)
      IF (LAROME) THEN
        ZTMPAF(KSTART:KPROF,1:KFLEV)=&
        & ZDELPREV(KSTART:KPROF,1:KFLEV)*PGFL(KSTART:KPROF,1:KFLEV,YQ%MP)
        CALL NEW_ADD_FIELD_3D(YDML_DIAG%YRMDDH,ZTMPAF,'VQV',YDDDH)
        ZTMPAF(KSTART:KPROF,1:KFLEV)=&
        & ZDELPREV(KSTART:KPROF,1:KFLEV)*PUZ(KSTART:KPROF,1:KFLEV)
        CALL NEW_ADD_FIELD_3D(YDML_DIAG%YRMDDH,ZTMPAF,'VUU',YDDDH)
        ZTMPAF(KSTART:KPROF,1:KFLEV)=&
        & ZDELPREV(KSTART:KPROF,1:KFLEV)*PVM(KSTART:KPROF,1:KFLEV)
        CALL NEW_ADD_FIELD_3D(YDML_DIAG%YRMDDH,ZTMPAF,'VVV',YDDDH)
        ZTMPAF(KSTART:KPROF,1:KFLEV)=&
        & ZDELPREV(KSTART:KPROF,1:KFLEV)*PKENE(KSTART:KPROF,1:KFLEV)
        CALL NEW_ADD_FIELD_3D(YDML_DIAG%YRMDDH,ZTMPAF,'VKK',YDDDH)
        ZTMPAF(KSTART:KPROF,1:KFLEV)=&
        & ZDELPREV(KSTART:KPROF,1:KFLEV)*ZCPT(KSTART:KPROF,1:KFLEV)
        CALL NEW_ADD_FIELD_3D(YDML_DIAG%YRMDDH,ZTMPAF,'VCT',YDDDH)
      ENDIF
    ELSE
     CALL ADD_FIELD_3D(YDML_DIAG%YRLDDH,ZDELPREV,'VPP','V','ARP',.TRUE.,.TRUE.)
     ZTMPAF(KSTART:KPROF,1:KFLEV)=&
      & ZDELPREV(KSTART:KPROF,1:KFLEV)*PAPHIF(KSTART:KPROF,1:KFLEV)
     CALL ADD_FIELD_3D(YDML_DIAG%YRLDDH,ZTMPAF,'VEP','V','ARP',.TRUE.,.TRUE.)
     ZTMPAF(KSTART:KPROF,1:KFLEV)=&
      & ZDELPREV(KSTART:KPROF,1:KFLEV)*PRH(KSTART:KPROF,1:KFLEV)
     CALL ADD_FIELD_3D(YDML_DIAG%YRLDDH,ZTMPAF,'VHR','V','ARP',.TRUE.,.TRUE.)
     IF (LAROME) THEN
       ZTMPAF(KSTART:KPROF,1:KFLEV)=&
        & ZDELPREV(KSTART:KPROF,1:KFLEV)*PGFL(KSTART:KPROF,1:KFLEV,IPTR_Q)
       CALL ADD_FIELD_3D(YDML_DIAG%YRLDDH,ZTMPAF,'VQV','V','ARP',.TRUE.,.TRUE.)
       ZTMPAF(KSTART:KPROF,1:KFLEV)=&
        & ZDELPREV(KSTART:KPROF,1:KFLEV)*PUZ(KSTART:KPROF,1:KFLEV)
       CALL ADD_FIELD_3D(YDML_DIAG%YRLDDH,ZTMPAF,'VUU','V','ARP',.TRUE.,.TRUE.)
       ZTMPAF(KSTART:KPROF,1:KFLEV)=&
        & ZDELPREV(KSTART:KPROF,1:KFLEV)*PVM(KSTART:KPROF,1:KFLEV)
       CALL ADD_FIELD_3D(YDML_DIAG%YRLDDH,ZTMPAF,'VVV','V','ARP',.TRUE.,.TRUE.)
       ZTMPAF(KSTART:KPROF,1:KFLEV)=&
        & ZDELPREV(KSTART:KPROF,1:KFLEV)*PKENE(KSTART:KPROF,1:KFLEV)
       CALL ADD_FIELD_3D(YDML_DIAG%YRLDDH,ZTMPAF,'VKK','V','ARP',.TRUE.,.TRUE.)
       ZTMPAF(KSTART:KPROF,1:KFLEV)=&
        & ZDELPREV(KSTART:KPROF,1:KFLEV)*ZCPT(KSTART:KPROF,1:KFLEV)
       CALL ADD_FIELD_3D(YDML_DIAG%YRLDDH,ZTMPAF,'VCT','V','ARP',.TRUE.,.TRUE.)
     ENDIF
    ENDIF
  ELSE
    IDHCV=0
    DO JLEV = 1, KFLEV
      DO JROF=KSTART,KPROF
!       VPP
        PDHCV(JROF,JLEV,IDHCV+1) = ZDELPREV(JROF,JLEV)
!       VQV
        PDHCV(JROF,JLEV,IDHCV+2) = PGFL(JROF,JLEV,IPTR_Q)*ZDELPREV(JROF,JLEV)
!       VUU
        PDHCV(JROF,JLEV,IDHCV+3) = PUZ(JROF,JLEV)*ZDELPREV(JROF,JLEV)
!       VVV
        PDHCV(JROF,JLEV,IDHCV+4) = PVM(JROF,JLEV)*ZDELPREV(JROF,JLEV)
!       VKK
        PDHCV(JROF,JLEV,IDHCV+5) = PKENE(JROF,JLEV)*ZDELPREV(JROF,JLEV)
!       VCT
        PDHCV(JROF,JLEV,IDHCV+6) = ZCPT(JROF,JLEV)*ZDELPREV(JROF,JLEV)
!       VEP
        PDHCV(JROF,JLEV,IDHCV+7) = PAPHIF(JROF,JLEV)*ZDELPREV(JROF,JLEV)
!       VHR
        PDHCV(JROF,JLEV,IDHCV+8) = PRH(JROF,JLEV)*ZDELPREV(JROF,JLEV)
      ENDDO
    ENDDO
  ENDIF
    
  SELECT CASE (NPHY)
    CASE(JPHYARO)
      IF (LFLEXDIA) THEN
        IF (LDDH_OMP) THEN
          ZTMPAF(KSTART:KPROF,1:KFLEV)=&
           & ZDELPREV(KSTART:KPROF,1:KFLEV)*PQL(KSTART:KPROF,1:KFLEV)
          CALL NEW_ADD_FIELD_3D(YDML_DIAG%YRMDDH,ZTMPAF,'VQL',YDDDH)
          ZTMPAF(KSTART:KPROF,1:KFLEV)=&
           & ZDELPREV(KSTART:KPROF,1:KFLEV)*PQI(KSTART:KPROF,1:KFLEV)
          CALL NEW_ADD_FIELD_3D(YDML_DIAG%YRMDDH,ZTMPAF,'VQI',YDDDH)
          ZTMPAF(KSTART:KPROF,1:KFLEV)=&
           & ZDELPREV(KSTART:KPROF,1:KFLEV)*PGFL(KSTART:KPROF,1:KFLEV,IPTR_R)
          CALL NEW_ADD_FIELD_3D(YDML_DIAG%YRMDDH,ZTMPAF,'VQR',YDDDH)
          ZTMPAF(KSTART:KPROF,1:KFLEV)=&
           & ZDELPREV(KSTART:KPROF,1:KFLEV)*PGFL(KSTART:KPROF,1:KFLEV,IPTR_S)
          CALL NEW_ADD_FIELD_3D(YDML_DIAG%YRMDDH,ZTMPAF,'VQS',YDDDH)
          ZTMPAF(KSTART:KPROF,1:KFLEV)=&
           & ZDELPREV(KSTART:KPROF,1:KFLEV)*PGFL(KSTART:KPROF,1:KFLEV,IPTR_G)
          CALL NEW_ADD_FIELD_3D(YDML_DIAG%YRMDDH,ZTMPAF,'VQG',YDDDH)
          IF (YH%MP /= NUNDEFLD ) THEN
            ZTMPAF(KSTART:KPROF,1:KFLEV)=&
            & ZDELPREV(KSTART:KPROF,1:KFLEV)*PGFL(KSTART:KPROF,1:KFLEV,IPTR_H)
            CALL NEW_ADD_FIELD_3D(YDML_DIAG%YRMDDH,ZTMPAF,'VQH',YDDDH)
          ENDIF
          !ZTMPAF(KSTART:KPROF,1:KFLEV)= &
          ! & ZDELPREV(KSTART:KPROF,1:KFLEV)*PNEB(KSTART:KPROF,1:KFLEV)
          !CALL NEW_ADD_FIELD_3D(ZTMPAF,'VNT',YDDDH)
          ZTMPAF(KSTART:KPROF,1:KFLEV)=&
           & ZDELPREV(KSTART:KPROF,1:KFLEV)*POMEGA(KSTART:KPROF,1:KFLEV)
          CALL NEW_ADD_FIELD_3D(YDML_DIAG%YRMDDH,ZTMPAF,'VOM',YDDDH)
          ZTMPAF(KSTART:KPROF,1:KFLEV)=&
           & ZDELPREV(KSTART:KPROF,1:KFLEV)*PGFL(KSTART:KPROF,1:KFLEV,IPTR_TKE)
          CALL NEW_ADD_FIELD_3D(YDML_DIAG%YRMDDH,ZTMPAF,'VKE',YDDDH)
        ELSE
          ZTMPAF(KSTART:KPROF,1:KFLEV)=&
           & ZDELPREV(KSTART:KPROF,1:KFLEV)*PQL(KSTART:KPROF,1:KFLEV)
          CALL ADD_FIELD_3D(YDML_DIAG%YRLDDH,ZTMPAF,'VQL','V','ARP',.TRUE.,.TRUE.)
          ZTMPAF(KSTART:KPROF,1:KFLEV)=&
           & ZDELPREV(KSTART:KPROF,1:KFLEV)*PQI(KSTART:KPROF,1:KFLEV)
          CALL ADD_FIELD_3D(YDML_DIAG%YRLDDH,ZTMPAF,'VQI','V','ARP',.TRUE.,.TRUE.)
          ZTMPAF(KSTART:KPROF,1:KFLEV)=&
           & ZDELPREV(KSTART:KPROF,1:KFLEV)*PGFL(KSTART:KPROF,1:KFLEV,IPTR_R)
          CALL ADD_FIELD_3D(YDML_DIAG%YRLDDH,ZTMPAF,'VQR','V','ARP',.TRUE.,.TRUE.)
          ZTMPAF(KSTART:KPROF,1:KFLEV)=&
           & ZDELPREV(KSTART:KPROF,1:KFLEV)*PGFL(KSTART:KPROF,1:KFLEV,IPTR_S)
          CALL ADD_FIELD_3D(YDML_DIAG%YRLDDH,ZTMPAF,'VQS','V','ARP',.TRUE.,.TRUE.)
          ZTMPAF(KSTART:KPROF,1:KFLEV)=&
           & ZDELPREV(KSTART:KPROF,1:KFLEV)*PGFL(KSTART:KPROF,1:KFLEV,IPTR_G)
          CALL ADD_FIELD_3D(YDML_DIAG%YRLDDH,ZTMPAF,'VQG','V','ARP',.TRUE.,.TRUE.)
          IF (IPTR_H /= NUNDEFLD ) THEN
             ZTMPAF(KSTART:KPROF,1:KFLEV)=&
             & ZDELPREV(KSTART:KPROF,1:KFLEV)*PGFL(KSTART:KPROF,1:KFLEV,IPTR_H)
             CALL ADD_FIELD_3D(YDML_DIAG%YRLDDH,ZTMPAF,'VQH','V','ARP',.TRUE.,.TRUE.)
          ENDIF 
          ZTMPAF(KSTART:KPROF,1:KFLEV)=&
           & ZDELPREV(KSTART:KPROF,1:KFLEV)*PNEB(KSTART:KPROF,1:KFLEV)
          CALL ADD_FIELD_3D(YDML_DIAG%YRLDDH,ZTMPAF,'VNT','V','ARP',.TRUE.,.TRUE.)
          ZTMPAF(KSTART:KPROF,1:KFLEV)=&
           & ZDELPREV(KSTART:KPROF,1:KFLEV)*POMEGA(KSTART:KPROF,1:KFLEV)
          CALL ADD_FIELD_3D(YDML_DIAG%YRLDDH,ZTMPAF,'VOM','V','ARP',.TRUE.,.TRUE.)
          ZTMPAF(KSTART:KPROF,1:KFLEV)=&
          & ZDELPREV(KSTART:KPROF,1:KFLEV)*PGFL(KSTART:KPROF,1:KFLEV,IPTR_TKE)
          CALL ADD_FIELD_3D(YDML_DIAG%YRLDDH,ZTMPAF,'VKE','V','ARP',.TRUE.,.TRUE.)
        ENDIF
      ELSE
         DO JLEV = 1, KFLEV
          DO JROF=KSTART,KPROF
!           VQL
            PDHCV(JROF,JLEV,IDHCV+9) = PQL(JROF,JLEV)*ZDELPREV(JROF,JLEV)
!           VQI 
            PDHCV(JROF,JLEV,IDHCV+10) = PQI(JROF,JLEV)*ZDELPREV(JROF,JLEV)
!           VQR
            PDHCV(JROF,JLEV,IDHCV+11) = PGFL(JROF,JLEV,IPTR_R)*ZDELPREV(JROF,JLEV)
!           VQS
            PDHCV(JROF,JLEV,IDHCV+12) = PGFL(JROF,JLEV,IPTR_S)*ZDELPREV(JROF,JLEV)
!           VQG
            PDHCV(JROF,JLEV,IDHCV+13) = PGFL(JROF,JLEV,IPTR_G)*ZDELPREV(JROF,JLEV)
!           VNT
            PDHCV(JROF,JLEV,IDHCV+9+NHDQLNVA+NHDPASVA)&
             & =PNEB(JROF,JLEV)*ZDELPREV(JROF,JLEV)  
!           VOM
!            PDHCV(JROF,JLEV,IDHCV+10+NHDQLNVA+NHDPASVA)&
!             & =POMEGA(JROF,JLEV)*ZDELPREV(JROF,JLEV)  
            PDHCV(JROF,JLEV,IDHCV+10+NHDQLNVA+NHDPASVA)=0.0_JPRB
!           VTE
            PDHCV(JROF,JLEV,IDHCV+11+NHDQLNVA+NHDPASVA)&
             & = PGFL(JROF,JLEV,IPTR_TKE)*ZDELPREV(JROF,JLEV)
          ENDDO
        ENDDO
      ENDIF
      
    CASE DEFAULT
      IF(LFLEXDIA) THEN
        IF(L3MT.OR.LSTRAPRO) THEN
          IF (LDDH_OMP) THEN
            !VQR
            ZTMPAF(KSTART:KPROF,1:KFLEV)=&
            &ZDELPREV(KSTART:KPROF,1:KFLEV)*PGFL(KSTART:KPROF,1:KFLEV,IPTR_R)
            CALL NEW_ADD_FIELD_3D(YDML_DIAG%YRMDDH,ZTMPAF,'VQR',YDDDH)
            !VQS
            ZTMPAF(KSTART:KPROF,1:KFLEV)=&
            &ZDELPREV(KSTART:KPROF,1:KFLEV)*PGFL(KSTART:KPROF,1:KFLEV,IPTR_S)
            CALL NEW_ADD_FIELD_3D(YDML_DIAG%YRMDDH,ZTMPAF,'VQS',YDDDH)
          ELSE 
            !VQR
            ZTMPAF(KSTART:KPROF,1:KFLEV)=&
             &ZDELPREV(KSTART:KPROF,1:KFLEV)*PGFL(KSTART:KPROF,1:KFLEV,IPTR_R)
            CALL ADD_FIELD_3D(YDML_DIAG%YRLDDH,ZTMPAF,'VQR','V','ARP',.TRUE.,.TRUE.)
            !VQS
            ZTMPAF(KSTART:KPROF,1:KFLEV)=&
            &ZDELPREV(KSTART:KPROF,1:KFLEV)*PGFL(KSTART:KPROF,1:KFLEV,IPTR_S)
            CALL ADD_FIELD_3D(YDML_DIAG%YRLDDH,ZTMPAF,'VQS','V','ARP',.TRUE.,.TRUE.)
          ENDIF  
        ENDIF
        IF (LDDH_OMP) THEN
          !VOM
          ZTMPAF(KSTART:KPROF,1:KFLEV)=&
          &ZDELPREV(KSTART:KPROF,1:KFLEV)*POMEGA(KSTART:KPROF,1:KFLEV)
          CALL NEW_ADD_FIELD_3D(YDML_DIAG%YRMDDH,ZTMPAF,'VOM',YDDDH)
        ELSE 
          !VOM
          ZTMPAF(KSTART:KPROF,1:KFLEV)=&
          &ZDELPREV(KSTART:KPROF,1:KFLEV)*POMEGA(KSTART:KPROF,1:KFLEV)
          CALL ADD_FIELD_3D(YDML_DIAG%YRLDDH,ZTMPAF,'VOM','V','ARP',.TRUE.,.TRUE.)
        ENDIF
      ELSE
      
        IF(LHDQLN) THEN
          DO JLEV = 1, KFLEV
            DO JROF=KSTART,KPROF
!             VQL
              PDHCV(JROF,JLEV,IDHCV+9)=PQL(JROF,JLEV)*ZDELPREV(JROF,JLEV)
!             VQI
              PDHCV(JROF,JLEV,IDHCV+10)=PQI(JROF,JLEV)*ZDELPREV(JROF,JLEV)
            ENDDO
          ENDDO
!         VQR and VQS
          IF ( L3MT.OR.LSTRAPRO )THEN
            DO JLEV = 1, KFLEV
              DO JROF=KSTART,KPROF
                PDHCV(JROF,JLEV,IDHCV+11)=PGFL(JROF,JLEV,IPTR_R)*ZDELPREV(JROF,JLEV)
                PDHCV(JROF,JLEV,IDHCV+12)=PGFL(JROF,JLEV,IPTR_S)*ZDELPREV(JROF,JLEV)
              ENDDO
            ENDDO
          ELSE
            DO JLEV = 1, KFLEV
              DO JROF=KSTART,KPROF
                PDHCV(JROF,JLEV,IDHCV+11)=PQR(JROF,JLEV)*ZDELPREV(JROF,JLEV)
                PDHCV(JROF,JLEV,IDHCV+12)=PQS(JROF,JLEV)*ZDELPREV(JROF,JLEV)
              ENDDO
            ENDDO
          ENDIF
        ENDIF
        IF ( LHDPAS ) THEN
!         VEXT
          IF(.NOT.(NUNDEFLD == 1 .AND. NVEXTR == 1))THEN
            DO JLEV = 1, KFLEV
              DO JROF=KSTART,KPROF
                IF (NVECOUT > 0) THEN
                  ! First place in extra/passive variables reserved for precip fraction
                  PDHCV(JROF,JLEV,IDHCV+8+NHDQLNVA+1)=PCOVPTOT(JROF,JLEV,1)*ZDELPREV(JROF,JLEV)
                ENDIF
                IF(.NOT.(NUNDEFLD == 1 .AND. NVEXTR == 1))THEN
                  DO JEXT=1,NVEXTR
                    PDHCV(JROF,JLEV,IDHCV+8+NHDQLNVA+NVECOUT+JEXT)=PEXT(JROF,JLEV,JEXT)*ZDELPREV(JROF,JLEV)
                  ENDDO
                ELSE
                  ! For bounds-checking case when there are no real extra variables - DS
                  DO JEXT=1,NVEXTR
                    PDHCV(JROF,JLEV,IDHCV+8+NHDQLNVA+NVECOUT+JEXT)=0.0_JPRB
                  ENDDO
                ENDIF
              ENDDO
            ENDDO
          ENDIF
        ENDIF
        DO JLEV = 1, KFLEV
          DO JROF=KSTART,KPROF
!           VNT
            PDHCV(JROF,JLEV,IDHCV+9+NHDQLNVA+NHDPASVA)&
             & =PNEB(JROF,JLEV)*ZDELPREV(JROF,JLEV)  
!           VOM
            PDHCV(JROF,JLEV,IDHCV+10+NHDQLNVA+NHDPASVA)&
             & =POMEGA(JROF,JLEV)*ZDELPREV(JROF,JLEV)  
          ENDDO
        ENDDO

        IF (LPTKE) THEN
          DO JLEV=1,KFLEV
            DO JROF=KSTART,KPROF
              ! TKE
              PDHCV(JROF,JLEV,IDHCV+10+NHDQLNVA+NHDPASVA+1)=PGFL(JROF,JLEV,IPTR_TKE)*ZDELPREV(JROF,JLEV)
            ENDDO
          ENDDO

          IF (LCOEFK_PTTE) THEN
            DO JLEV=1,KFLEV
              DO JROF=KSTART,KPROF
                ! TTE
                PDHCV(JROF,JLEV, IDHCV+10+NHDQLNVA+NHDPASVA+2) = PGFL(JROF,JLEV,IPTR_TTE)*ZDELPREV(JROF,JLEV)
              ENDDO
            ENDDO
          ENDIF

        ENDIF

      ENDIF
  END SELECT
ENDIF

IF ( LHDMCI ) THEN
  IDHCV=NDHVHK
  DO JLEV=1, KFLEV
    DO JROF=KSTART,KPROF
!     VA1
      PDHCV(JROF,JLEV,IDHCV+1) = ZM1(JROF,JLEV)*ZDELPREV(JROF,JLEV)
!     VA2
      PDHCV(JROF,JLEV,IDHCV+2) = ZM2(JROF,JLEV)*ZDELPREV(JROF,JLEV)
!     VA3
      PDHCV(JROF,JLEV,IDHCV+3) = ZM3(JROF,JLEV)*ZDELPREV(JROF,JLEV)
    ENDDO
  ENDDO
ENDIF

IF ( LHDENT ) THEN
  IDHCV=NDHVHK+NDHVMC
  DO JLEV = 1, KFLEV
    DO JROF=KSTART,KPROF
!     VSS
      PDHCV(JROF,JLEV,IDHCV+1) = ZENTR(JROF,JLEV)*ZDELPREV(JROF,JLEV)
    ENDDO
  ENDDO
ENDIF

!     ------------------------------------------------------------------
!*    2.- CALCUL DES TENDANCES DYNAMIQUES / DYNAMICAL TENDENCIES COMPUTATIONS

! 2.1/
! * Calcul du flux total de precipitations "Fp" sur les intercouches
!   pour traiter le cas "delta m=1".
! * Computation of the total rainfall flux "Fp" at half levels
!   to treat the case "delta m=1".

IF(NDPSFI == 1)THEN
  DO JLEV = 0, KFLEV
    DO JROF=KSTART,KPROF
      ZFPSUM(JROF,JLEV) =&
       & ( PFPLCL(JROF,JLEV) + PFPLCN(JROF,JLEV)&
       & + PFPLSL(JROF,JLEV) + PFPLSN(JROF,JLEV) )  
    ENDDO
  ENDDO
ELSE
  DO JLEV = 0, KFLEV
    DO JROF=KSTART,KPROF
      ZFPSUM(JROF,JLEV) = 0.0_JPRB
    ENDDO
  ENDDO
ENDIF

! 2.3/
! * Calcul sur les couches de/ computation at full levels of:
!    " - (1/g) grad { vec(V) X (Delta prehyd) } "
!    for X = 1; q; Ugeo; Vgeo; (U**2+V**2)/2; cp*T; g*z; ql; qi; M1; M2; M3; S.
!   Calcul sur les couches de/ computation at full levels of:
!    " - (1/g) (Delta prehyd) (DX/Dt)_adiab "
!    for X = Ugeo; Vgeo; 0.5*(U**2+V**2); cp*T; M1; M2; M3.

! 2.3.1/
! -- Variables "X" autres que "M1", "M2", "M3" et "S".
!    Variables "X" other than "M1", "M2", "M3" and "S".

IF ( LHDHKS ) THEN
  IDHCV=NDHVV

! 2.3.1.1/
! RHS de l'equation du vent horizontal.
! RHS for horizontal wind equation.

  DO JLEV=1,KFLEV
    DO JROF=KSTART,KPROF
      ZRHSU(JROF,JLEV)=PATND(JROF,JLEV,YDDYNA%YYTTND%M_TNDU)
      ZRHSV(JROF,JLEV)=PATND(JROF,JLEV,YDDYNA%YYTTND%M_TNDV)
      ZRHSZO(JROF)=YDGSGEOM%GNORDM(JROF)*ZRHSU(JROF,JLEV)-YDGSGEOM%GNORDL(JROF)*ZRHSV(JROF,JLEV)
      ZRHSME(JROF)=YDGSGEOM%GNORDL(JROF)*ZRHSU(JROF,JLEV)+YDGSGEOM%GNORDM(JROF)*ZRHSV(JROF,JLEV)
      ZRHSU(JROF,JLEV)=ZRHSZO(JROF)
      ZRHSV(JROF,JLEV)=ZRHSME(JROF)
    ENDDO
  ENDDO

! 2.3.1.2/
! RHS de l'equation de l'energie cinetique.
! RHS for kinetic energy equation.

  DO JLEV=1,KFLEV
    DO JROF=KSTART,KPROF
      ZRHSK(JROF,JLEV)=PGMV(JROF,JLEV,IPTR_U)*PATND(JROF,JLEV,YDDYNA%YYTTND%M_TNDU)&
       & +PGMV(JROF,JLEV,IPTR_V)*PATND(JROF,JLEV,YDDYNA%YYTTND%M_TNDV)
    ENDDO
  ENDDO

! 2.3.1.3/
! Le calcul et stockage des divergences de flux.
! Computation and storage of flux divergences.

! 2.3.1.3.a/ "X"=1.

  DO JLEV = 1, KFLEV
    DO JROF=KSTART,KPROF
!     TPPDIVFLUHOR
      PDHCV(JROF,JLEV,IDHCV+1) = -ZUSRG*ZT_DIVDP(JROF,JLEV)
    ENDDO
  ENDDO

! 2.3.1.3.b/ "X"=q.

  DO JLEV = 1, KFLEV
    DO JROF=KSTART,KPROF
!     TQVDIVFLUHOR
      IF(YQ%LCDERS) THEN
        ZADVQV =  PGMV(JROF,JLEV,IPTR_U)*PGFL(JROF,JLEV,IPTR_QL)+PGMV(JROF,JLEV,IPTR_V)*PGFL(JROF,JLEV,IPTR_QM)
      ELSE
        ZADVQV = 0.0_JPRB
      ENDIF
      PDHCV(JROF,JLEV,IDHCV+2) = ZUSRG*(PGFL(JROF,JLEV,IPTR_Q)*&
       & (-ZT_DIVDP(JROF,JLEV))-ZT_DELP(JROF,JLEV)*ZADVQV)  
    ENDDO
  ENDDO

! 2.3.1.3.c/ "X"=vec(V).

  DO JLEV = 1, KFLEV
    DO JROF=KSTART,KPROF
!     TUUDIVFLUHOR
      PDHCV(JROF,JLEV,IDHCV+3) = ZUSRG*(PUZ(JROF,JLEV)*&
       & (-ZT_DIVDP(JROF,JLEV))-ZT_DELP(JROF,JLEV)*&
       & (YDGSGEOM%GNORDM(JROF)*ZGKEU(JROF,JLEV)-YDGSGEOM%GNORDL(JROF)*ZGKEV(JROF,JLEV)))  
!     TVVDIVFLUHOR
      PDHCV(JROF,JLEV,IDHCV+4) = ZUSRG*(PVM(JROF,JLEV)*&
       & (-ZT_DIVDP(JROF,JLEV))-ZT_DELP(JROF,JLEV)*&
       & (YDGSGEOM%GNORDL(JROF)*ZGKEU(JROF,JLEV)+YDGSGEOM%GNORDM(JROF)*ZGKEV(JROF,JLEV)))  
!     TUUFFVGADPSG
      PDHCV(JROF,JLEV,IDHCV+5) = ZDELPREV(JROF,JLEV)*ZUSRG*ZRHSU(JROF,JLEV)
!     TUUFFVGADPSG
      PDHCV(JROF,JLEV,IDHCV+6) = ZDELPREV(JROF,JLEV)*ZUSRG*ZRHSV(JROF,JLEV)
    ENDDO
  ENDDO

! 2.3.1.3.d/ "X"= 0.5 (U**2+V**2).

  DO JLEV = 1, KFLEV
    DO JROF=KSTART,KPROF
!     TKKDIVFLUHOR
      PDHCV(JROF,JLEV,IDHCV+7) = ZUSRG*(PKENE(JROF,JLEV)*&
       & (-ZT_DIVDP(JROF,JLEV))-ZT_DELP(JROF,JLEV)*&
       & (ZGKEU(JROF,JLEV)*PGMV(JROF,JLEV,IPTR_U)+ZGKEV(JROF,JLEV)*PGMV(JROF,JLEV,IPTR_V)))  
!     TKKCONVERSI1
      PDHCV(JROF,JLEV,IDHCV+8) = ZDELPREV(JROF,JLEV)*ZUSRG*ZRHSK(JROF,JLEV)
    ENDDO
  ENDDO

! 2.3.1.3.e/ "X"= cp*T.

! TCTDIVFLUHOR
  IF(YQ%LCDERS) THEN
    DO JLEV = 1, KFLEV
      DO JROF=KSTART,KPROF
        ZADVQV =  PGMV(JROF,JLEV,IPTR_U)*PGFL(JROF,JLEV,IPTR_QL)+PGMV(JROF,JLEV,IPTR_V)*PGFL(JROF,JLEV,IPTR_QM)
        ZADVT =  PGMV(JROF,JLEV,IPTR_U)*PGMV(JROF,JLEV,IPTR_TL) +PGMV(JROF,JLEV,IPTR_V)*PGMV(JROF,JLEV,IPTR_TM)
        PDHCV(JROF,JLEV,IDHCV+9) = ZUSRG*(ZCPT(JROF,JLEV)*&
         & (-ZT_DIVDP(JROF,JLEV))-ZT_DELP(JROF,JLEV)*&
         & (PRCP(JROF,JLEV,IPTR_CP)*ZADVT+PGMV(JROF,JLEV,IPTR_T)*ZCVMCD*ZADVQV))  
!       TCTCONVERSI2
        PDHCV(JROF,JLEV,IDHCV+10) = ZDELPREV(JROF,JLEV)*ZUSRG*&
         & Z_CP_TIMES_CONV(JROF,JLEV)  
      ENDDO
    ENDDO
  ELSE
    DO JLEV = 1, KFLEV
      DO JROF=KSTART,KPROF
        ZADVQV = 0.0_JPRB
        ZADVT =  PGMV(JROF,JLEV,IPTR_U)*PGMV(JROF,JLEV,IPTR_TL) +PGMV(JROF,JLEV,IPTR_V)*PGMV(JROF,JLEV,IPTR_TM)
        PDHCV(JROF,JLEV,IDHCV+9) = ZUSRG*(ZCPT(JROF,JLEV)*&
         & (-ZT_DIVDP(JROF,JLEV))-ZT_DELP(JROF,JLEV)*&
         & (PRCP(JROF,JLEV,IPTR_CP)*ZADVT+PGMV(JROF,JLEV,IPTR_T)*ZCVMCD*ZADVQV))  
!       TCTCONVERSI2
        PDHCV(JROF,JLEV,IDHCV+10) = ZDELPREV(JROF,JLEV)*ZUSRG*&
         & Z_CP_TIMES_CONV(JROF,JLEV)  
      ENDDO
    ENDDO
  ENDIF
  IF(NDPSFI == 1) THEN
    ! K.Y.: this code does not seem to have been properly updated
    !       relatively to the modifications of the "delta m=1"
    !       formulation introduced in cy29t3.
    DO JLEV = 2, KFLEV-1
      DO JROF=KSTART,KPROF
!       TCTCONVERSI3 for layers 2 to L-1 if "ndpsfi=1".
        PDHCV(JROF,JLEV,IDHCV+11) =&
         & 0.5_JPRB*(ZFPSUM(JROF,JLEV)+ZFPSUM(JROF,JLEV-1))&
         & *(PAPHI(JROF,JLEV)-PAPHI(JROF,JLEV-1)&
         & +0.5_JPRB*(PKENE(JROF,JLEV+1)-PKENE(JROF,JLEV-1)))  
      ENDDO
    ENDDO
    DO JROF=KSTART,KPROF
!     TCTCONVERSI3 for layer l=1 if "ndpsfi=1".
      PDHCV(JROF,1,IDHCV+11) =0.5_JPRB*ZFPSUM(JROF,1)&
       & *(PAPHI(JROF,1)-PAPHI(JROF,0)&
       & +PKENE(JROF,2)-PKENE(JROF,1))  
!     TCTCONVERSI3 for layer l=L if "ndpsfi=1".
      PDHCV(JROF,KFLEV,IDHCV+11) =&
       & 0.5_JPRB*(ZFPSUM(JROF,KFLEV)+ZFPSUM(JROF,KFLEV-1))&
       & *(PAPHI(JROF,KFLEV)-PAPHI(JROF,KFLEV-1)&
       & -2.0_JPRB*PKENE(JROF,KFLEV))  
    ENDDO
  ELSE
    DO JLEV = 1, KFLEV
      DO JROF=KSTART,KPROF
!       TCTCONVERSI3 if "ndpsfi=0".
        PDHCV(JROF,JLEV,IDHCV+11) = 0.0_JPRB
      ENDDO
    ENDDO
  ENDIF

! 2.3.1.3.f/ "X"= gz.

  DO JLEV = 1, KFLEV
    DO JROF=KSTART,KPROF
!     TEPDIVFLUHOR
!     On suppose que c'est "cp*T + gz + 0.5*(U**2+V**2)" qui est la
!     quantite energetique conservee, donc on a besoin de calculer la
!     divergence des flux horizontaux et le gradient horizontal de "gz".
!     One assumes that the energetic quantity which is conserved is
!     "cp*T + gz + 0.5*(U**2+V**2)", so one has to compute the
!     horizontal flux divergence and the horizontal gradient of "gz".
      ZADPHI=PGMV(JROF,JLEV,IPTR_U)*PAPHIFL(JROF,JLEV)+PGMV(JROF,JLEV,IPTR_V)*PAPHIFM(JROF,JLEV)
      PDHCV(JROF,JLEV,IDHCV+12) = ZUSRG*(PAPHIF(JROF,JLEV)*&
       & (-ZT_DIVDP(JROF,JLEV))-ZT_DELP(JROF,JLEV)*ZADPHI)  
    ENDDO
  ENDDO

! 2.3.1.3.g/ "X"= 1 (addendum).

  IF (LFLEXDIA) THEN
    ZTMPAF(KSTART:KPROF,1:KFLEV)=ZDELPREV(KSTART:KPROF,1:KFLEV)*ZUSRG
    IF (LDDH_OMP) THEN
      CALL NEW_ADD_FIELD_3D(YDML_DIAG%YRMDDH,ZTMPAF,'PPP',YDDDH,CDTYPE='T')
    ELSE
      CALL ADD_FIELD_3D(YDML_DIAG%YRLDDH,ZTMPAF,'PPP','T','ARP',.TRUE.,.TRUE.)
    ENDIF
  ELSE
    DO JLEV = 1, KFLEV
      DO JROF=KSTART,KPROF
!       PPP
        PDHCV(JROF,JLEV,IDHCV+13) = ZDELPREV(JROF,JLEV)*ZUSRG
      ENDDO
    ENDDO
  ENDIF


! 2.3.1.3.h/ "X"= ql.

  IF(LHDQLN) THEN
    DO JLEV = 1, KFLEV
      DO JROF=KSTART,KPROF
!       TQLDIVFLUHOR
        IF(YL%LCDERS) THEN
          ZADVQL =  PGMV(JROF,JLEV,IPTR_U)*PGFL(JROF,JLEV,IPTR_LL)+&
           & PGMV(JROF,JLEV,IPTR_V)*PGFL(JROF,JLEV,IPTR_LM)
        ELSE
          ZADVQL=0.0_JPRB
        ENDIF
        PDHCV(JROF,JLEV,IDHCV+14) = ZUSRG*(PQL(JROF,JLEV)*&
         & (-ZT_DIVDP(JROF,JLEV))-ZT_DELP(JROF,JLEV)*ZADVQL)  
      ENDDO
    ENDDO
  ENDIF

! 2.3.1.3.i/ "X"= qi.

  IF(LHDQLN) THEN
    DO JLEV = 1, KFLEV
      DO JROF=KSTART,KPROF
!       TQIDIVFLUHOR
        IF(YI%LCDERS) THEN
          ZADVQI =  PGMV(JROF,JLEV,IPTR_U)*PGFL(JROF,JLEV,IPTR_IL)+&
           & PGMV(JROF,JLEV,IPTR_V)*PGFL(JROF,JLEV,IPTR_IM)
        ELSE
          ZADVQI=0.0_JPRB
        ENDIF
        PDHCV(JROF,JLEV,IDHCV+15) = ZUSRG*(PQI(JROF,JLEV)*&
         & (-ZT_DIVDP(JROF,JLEV))-ZT_DELP(JROF,JLEV)*ZADVQI)  
      ENDDO
    ENDDO
  ENDIF

! 2.3.1.3.j/ "X"= qr.

  IF(LHDQLN) THEN
!   TQRDIVFLUHOR
    IF(YR%LCDERS) THEN
      DO JLEV = 1, KFLEV
        DO JROF=KSTART,KPROF
          ZADVQR =  PGMV(JROF,JLEV,IPTR_U)*PGFL(JROF,JLEV,IPTR_RL)+&
           & PGMV(JROF,JLEV,IPTR_V)*PGFL(JROF,JLEV,IPTR_RM)
          PDHCV(JROF,JLEV,IDHCV+16) = ZUSRG*(PQR(JROF,JLEV)*&
           & (-ZT_DIVDP(JROF,JLEV))-ZT_DELP(JROF,JLEV)*ZADVQR)  
        ENDDO
      ENDDO
    ELSE
      DO JLEV = 1, KFLEV
        DO JROF=KSTART,KPROF
          ZADVQR=0.0_JPRB
          PDHCV(JROF,JLEV,IDHCV+16) = ZUSRG*(PQR(JROF,JLEV)*&
           & (-ZT_DIVDP(JROF,JLEV))-ZT_DELP(JROF,JLEV)*ZADVQR)  
        ENDDO
      ENDDO
    ENDIF
  ENDIF

! 2.3.1.3.k/ "X"= qs.

  IF(LHDQLN) THEN
    ! TQSDIVFLUHOR
    IF(YS%LCDERS) THEN
      DO JLEV = 1, KFLEV
        DO JROF=KSTART,KPROF
          ZADVQS =  PGMV(JROF,JLEV,IPTR_U)*PGFL(JROF,JLEV,IPTR_SL)+&
           & PGMV(JROF,JLEV,IPTR_V)*PGFL(JROF,JLEV,IPTR_SM)
          PDHCV(JROF,JLEV,IDHCV+17) = ZUSRG*(PQS(JROF,JLEV)*&
           & (-ZT_DIVDP(JROF,JLEV))-ZT_DELP(JROF,JLEV)*ZADVQS)  
        ENDDO
      ENDDO
    ELSE
      DO JLEV = 1, KFLEV
        DO JROF=KSTART,KPROF
          ZADVQS=0.0_JPRB
          PDHCV(JROF,JLEV,IDHCV+17) = ZUSRG*(PQS(JROF,JLEV)*&
           & (-ZT_DIVDP(JROF,JLEV))-ZT_DELP(JROF,JLEV)*ZADVQS)  
        ENDDO
      ENDDO
    ENDIF
  ENDIF

  ! TKE + TTE 
  IF (LPTKE) THEN

!   2.3.1.3.l/ "X"=tke

    DO JLEV = 1, KFLEV
      DO JROF=KSTART,KPROF
!       TTKDIVFLUHOR
        IF(YTKE%LCDERS) THEN
          ZADVTKE =  PGMV(JROF,JLEV,IPTR_U)*PGFL(JROF,JLEV,IPTR_TKEL)+&
           & PGMV(JROF,JLEV,IPTR_V)*PGFL(JROF,JLEV,IPTR_TKEM)
        ELSE
          ZADVTKE=0.0_JPRB
        ENDIF
        PDHCV(JROF,JLEV,IDHCV+18) = ZUSRG*(-PGFL(JROF,JLEV,IPTR_TKE)*&
         & ZT_DIVDP(JROF,JLEV)-ZT_DELP(JROF,JLEV)*ZADVTKE)
      ENDDO
    ENDDO

    IF (LCOEFK_PTTE) THEN 

!   2.3.1.3.m/ "X"=tte

      DO JLEV = 1, KFLEV
        DO JROF=KSTART,KPROF
!         TTTDIVFLUHOR
          IF(YTTE%LCDERS) THEN
            ZADVTTE =  PGMV(JROF,JLEV,IPTR_U)*PGFL(JROF,JLEV,IPTR_TTEL)+&
             & PGMV(JROF,JLEV,IPTR_V)*PGFL(JROF,JLEV,IPTR_TTEM)
          ELSE
            ZADVTTE=0.0_JPRB
          ENDIF
          PDHCV(JROF,JLEV,IDHCV+19) = ZUSRG*(-PGFL(JROF,JLEV,IPTR_TTE)*&
           & ZT_DIVDP(JROF,JLEV)-ZT_DELP(JROF,JLEV)*ZADVTTE)
        ENDDO
      ENDDO

    ENDIF ! LCOEFK_PTTE

  ENDIF ! LPTKE

ENDIF

! 2.3.2/
! -- Variables "X" = "M1", "M2", "M3".

IF ( LHDMCI ) THEN

! * Code not updated for NH model.

  IDHCV=NDHVV+NDHAHKD
  DO JLEV = 1, KFLEV
    DO JROF=KSTART,KPROF
!     DIVERGENCE OF ANGULAR MOMENTUM
!       Only the term M div(v delta p) is computed, the remaining term
!       delta p v.grad M can be computed as soon as the gradient of velocity
!       is available.

!     grad deltaP ON THE REAL SPHERE
      ZDPREZ = YDVAB%VDELB(JLEV)*(YDGSGEOM%GNORDM(JROF)*PSPL(JROF)-YDGSGEOM%GNORDL(JROF)*PSPM(JROF))
      ZDPREM = YDVAB%VDELB(JLEV)*(YDGSGEOM%GNORDL(JROF)*PSPL(JROF)+YDGSGEOM%GNORDM(JROF)*PSPM(JROF))

      PDHCV(JROF,JLEV,IDHCV+1) = - ZUSRG * ZM1(JROF,JLEV)*&
       & (ZT_DIVDP(JROF,JLEV)+ZDPREZ*RA*ROMEGA*YDGSGEOM%GSQM2(JROF))  
      PDHCV(JROF,JLEV,IDHCV+2) = - ZUSRG * ZM2(JROF,JLEV)*&
       & (ZT_DIVDP(JROF,JLEV)+ZDPREZ*RA*ROMEGA*YDGSGEOM%GSQM2(JROF))  
      PDHCV(JROF,JLEV,IDHCV+3) = -ZUSRG * ZM3(JROF,JLEV)*&
       & (ZT_DIVDP(JROF,JLEV)+ZDPREZ*RA*ROMEGA*YDGSGEOM%GSQM2(JROF))  
!     ANGULAR MOMENTUM TENDENCY DUE TO ADJUSTMENT
      PDHCV(JROF,JLEV,IDHCV+4) = ZUSRG * RA * (&
       & - ZSLON(JROF)*(ZT_DELP(JROF,JLEV)*ZRTMER(JROF,JLEV)+&
       & PRCP(JROF,JLEV,IPTR_RD)*PGMV(JROF,JLEV,IPTR_T)*ZDPREM)&
       & + ZCLON(JROF)*YDGSGEOM%GEMU(JROF)*(ZT_DELP(JROF,JLEV)*ZRTZON(JROF,&
       & JLEV)+&
       & PRCP(JROF,JLEV,IPTR_RD)*PGMV(JROF,JLEV,IPTR_T)*ZDPREZ) )  
      PDHCV(JROF,JLEV,IDHCV+5) = ZUSRG * RA * (&
       & ZCLON(JROF)*(ZT_DELP(JROF,JLEV)*ZRTMER(JROF,JLEV)+&
       & PRCP(JROF,JLEV,IPTR_RD)*PGMV(JROF,JLEV,IPTR_T)*ZDPREM)&
       & + ZSLON(JROF)*YDGSGEOM%GEMU(JROF)*(ZT_DELP(JROF,JLEV)*ZRTZON(JROF,&
       & JLEV)+&
       & PRCP(JROF,JLEV,IPTR_RD)*PGMV(JROF,JLEV,IPTR_T)*ZDPREZ) )  
      PDHCV(JROF,JLEV,IDHCV+6) = ZUSRG * RA *&
       & YDGSGEOM%GSQM2(JROF)*(-ZT_DELP(JROF,JLEV)*ZRTZON(JROF,JLEV)-&
       & PRCP(JROF,JLEV,IPTR_RD)*PGMV(JROF,JLEV,IPTR_T)*ZDPREZ)  
!     NONAXIAL MOMENTS TENDENCY
      PDHCV(JROF,JLEV,IDHCV+7) =   ZUSRG*ZT_DELP(JROF,JLEV)*RA**2 *&
       & ROMEGA**2 *YDGSGEOM%GEMU(JROF)*YDGSGEOM%GSQM2(JROF)*ZSLON(JROF)  
      PDHCV(JROF,JLEV,IDHCV+8) = - ZUSRG*ZT_DELP(JROF,JLEV)*RA**2 *&
       & ROMEGA**2 *YDGSGEOM%GEMU(JROF)*YDGSGEOM%GSQM2(JROF)*ZCLON(JROF)  
    ENDDO
  ENDDO
ENDIF

! 2.3.3/
! -- Variable "X" = "S".

IF ( LHDENT ) THEN
  IDHCV=NDHVV+NDHAHKD+NDHAMCD

  DO JLEV = 1, KFLEV
    DO JROF=KSTART,KPROF
!     HORIZONTAL DIVERGENCE OF ENTROPY
      IF(YQ%LCDERS) THEN
        ZADVQV =  PGMV(JROF,JLEV,IPTR_U)*PGFL(JROF,JLEV,IPTR_QL)+PGMV(JROF,JLEV,IPTR_V)*PGFL(JROF,JLEV,IPTR_QM)
      ELSE
        ZADVQV = 0.0_JPRB
      ENDIF
      ZADVLNT =  (PGMV(JROF,JLEV,IPTR_U)*PGMV(JROF,JLEV,IPTR_TL) +&
       & PGMV(JROF,JLEV,IPTR_V)*PGMV(JROF,JLEV,IPTR_TM))/PGMV(JROF,JLEV,IPTR_T)  
      ZADVP = PGMV(JROF,JLEV,IPTR_U)*(PXYB(JROF,JLEV,YYTXYB0%M_RTGR)*PSPL(JROF)+ZQCHAL(JROF,JLEV))&
       & +PGMV(JROF,JLEV,IPTR_V)*(PXYB(JROF,JLEV,YYTXYB0%M_RTGR)*PSPM(JROF)+ZQCHAM(JROF,JLEV))
      PDHCV(JROF,JLEV,IDHCV+1)= ZUSRG*(ZENTR(JROF,JLEV)*&
       & (-ZT_DIVDP(JROF,JLEV))&
       & -ZT_DELP(JROF,JLEV)*((PENTRV(JROF,JLEV)-PENTRA(JROF,JLEV))*&
       & ZADVQV+PRCP(JROF,JLEV,IPTR_CP)*ZADVLNT-PRCP(JROF,JLEV,IPTR_RD)*ZADVP))  
    ENDDO
  ENDDO

! TENDENCY OF ENTROPY DUE TO TRANSPORT OF ENERGY BY PRECIPITATION
  IF (NDPSFI == 1) THEN
    ! K.Y.: this code does not seem to have been properly updated
    !       relatively to the modifications of the "delta m=1"
    !       formulation introduced in cy29t3.
    DO JLEV = 1, KFLEV
      DO JROF=KSTART,KPROF
        PDHCV(JROF,JLEV,IDHCV+2) =&
         & -(ZFPTOT(JROF,JLEV)+ZFPTOT(JROF,JLEV-1)) *&
         & (PAPHI(JROF,JLEV)-PAPHI(JROF,JLEV-1))/(2.0_JPRB*PGMV(JROF,JLEV,IPTR_T))  
      ENDDO
    ENDDO
  ELSE
    DO JLEV = 1, KFLEV
      DO JROF=KSTART,KPROF
        PDHCV(JROF,JLEV,IDHCV+2) =&
         & ( PENTRA(JROF,JLEV)*(ZFPTOT(JROF,JLEV)-ZFPTOT(JROF,JLEV-1))&
         & + RCPD*0.5_JPRB*ZDELLT(JROF,JLEV)*&
         & (ZFPTOT(JROF,JLEV)+ZFPTOT(JROF,JLEV-1)) )  
      ENDDO
    ENDDO
  ENDIF

ENDIF

!     ------------------------------------------------------------------
!*    3.- CALCUL DES FLUX VERTICAUX DYNAMIQUES
!         DYNAMICAL VERTICAL FLUXES COMPUTATIONS

! * Calcul sur les intercouches de/ computation at half levels of:
!    " (1/g) X etapt (d prehyd / d eta) "
!   for X = 1; q; Ugeo; Vgeo; (U**2+V**2)/2; cp*T; ql; qi; M1; M2; M3; S.

! 3.1/
! -- Variables "X" autres que "M1", "M2", "M3" et "S".
!    Variables "X" other than "M1", "M2", "M3" and "S".

IF ( LHDHKS ) THEN
  IDHCV=NDHVV + NDHAVD

! * CONDITION A LA LIMITE INFERIEURE / BOTTOM BOUNDARY CONDITION

!   Hypotheses faites ici / Assumptions for bottom values:
!   (U,V)_[bottom]=(U,V)_[l=L]
!   T_[bottom]=T_[surf]; cp_[bottom]=cpd+(cpv-cpd)*q_[surf]
!   q_[bottom]=q_[surf]; ql_[bottom]=ql_[l=L]; qi_[bottom]=qi_[l=L]

  DO JROF=KSTART,KPROF
!   FPPFLUVERTDYN
    PDHCV(JROF,KFLEV,IDHCV+1) = ZEVELH(JROF,KFLEV)*ZUSRG
!   FQVFLUVERTDYN
    PDHCV(JROF,KFLEV,IDHCV+2) = ZEVELH(JROF,KFLEV)*PQVS(JROF)*ZUSRG
!   FUUFLUVERTDYN
    PDHCV(JROF,KFLEV,IDHCV+3) = 0.0_JPRB
!   FVVFLUVERTDYN
    PDHCV(JROF,KFLEV,IDHCV+4) = 0.0_JPRB
!   FUUFLUDUAPLUI
    PDHCV(JROF,KFLEV,IDHCV+5) = 0.0_JPRB
!   FVVFLUDUAPLUI
    PDHCV(JROF,KFLEV,IDHCV+6) = 0.0_JPRB
!   FKKFLUVERTDYN
    PDHCV(JROF,KFLEV,IDHCV+7) = 0.0_JPRB
!   FKKFLUDUAPLUI
    PDHCV(JROF,KFLEV,IDHCV+8) = 0.0_JPRB
!   FCTFLUVERTDYN
    PDHCV(JROF,KFLEV,IDHCV+9) = ZEVELH(JROF,KFLEV)*&
     & (RCPD + ZCVMCD*PQVS(JROF))*PTSFC(JROF)*ZUSRG  
!   FEPCONVERSIFL
    PDHCV(JROF,KFLEV,IDHCV+10) = - PAPHI(JROF,KFLEV)*ZT_PSDIV(JROF,KFLEV)*ZUSRG
    IF(LHDQLN) THEN
!     FQLFLUVERTDYN
      PDHCV(JROF,KFLEV,IDHCV+11) = 0.0_JPRB
!     FQIFLUVERTDYN
      PDHCV(JROF,KFLEV,IDHCV+12) = 0.0_JPRB
!     FQRFLUVERTDYN
      PDHCV(JROF,KFLEV,IDHCV+13) = 0.0_JPRB
!     FQSFLUVERTDYN
      PDHCV(JROF,KFLEV,IDHCV+14) = 0.0_JPRB
    ENDIF
    ! TKE + TTE
    IF (LPTKE) THEN
!     FTKFLUVERTDYN    
      PDHCV(JROF,KFLEV,IDHCV+15) = 0.0_JPRB
      IF (LCOEFK_PTTE) THEN
!       FTTFLUVERTDYN
        PDHCV(JROF,KFLEV,IDHCV+16) = 0.0_JPRB
      ENDIF
    ENDIF
  ENDDO

! * INTERNIVEAUX NORMAUX / HALF LEVELS

  ZUS2RG = 0.5_JPRB*ZUSRG
  ZPFIS2 = 0.5_JPRB*ZDPSFI
  ! K.Y.: lines containing "ZDPSFI":
  !       this code does not seem to have been properly updated
  !       relatively to the modifications of the "delta m=1"
  !       formulation introduced in cy29t3.
  DO JLEV = 1, KFLEV-1
    DO JROF=KSTART,KPROF
!     FPPFLUVERTDYN
      PDHCV(JROF,JLEV,IDHCV+1) = ZEVELH(JROF,JLEV)*ZUSRG
!     FQVFLUVERTDYN
      PDHCV(JROF,JLEV,IDHCV+2) = ZUS2RG*( PGFL(JROF,JLEV,IPTR_Q)+&
       & PGFL(JROF,JLEV+1,IPTR_Q) )*ZEVELH(JROF,JLEV)  
!     FUUFLUVERTDYN
      PDHCV(JROF,JLEV,IDHCV+3) = ZUS2RG*( PUZ(JROF,JLEV)+&
       & PUZ(JROF,JLEV+1) )*ZEVELH(JROF,JLEV)  
!     FVVFLUVERTDYN
      PDHCV(JROF,JLEV,IDHCV+4) = ZUS2RG*( PVM(JROF,JLEV)+&
       & PVM(JROF,JLEV+1) )*ZEVELH(JROF,JLEV)  
!     FUUFLUDUAPLUI
      PDHCV(JROF,JLEV,IDHCV+5) = ZPFIS2*( PUZ(JROF,JLEV)+&
       & PUZ(JROF,JLEV+1) )*ZFPSUM(JROF,JLEV)  
!     FVVFLUDUAPLUI
      PDHCV(JROF,JLEV,IDHCV+6) = ZPFIS2*( PVM(JROF,JLEV)+&
       & PVM(JROF,JLEV+1) )*ZFPSUM(JROF,JLEV)  
!     FKKFLUVERTDYN
      PDHCV(JROF,JLEV,IDHCV+7) = ZUS2RG*( PKENE(JROF,JLEV)+&
       & PKENE(JROF,JLEV+1) )*ZEVELH(JROF,JLEV)  
!     FKKFLUDUAPLUI
      PDHCV(JROF,JLEV,IDHCV+8) = ZPFIS2*( PKENE(JROF,JLEV)+&
       & PKENE(JROF,JLEV+1) )*ZFPSUM(JROF,JLEV)  
!     FCTFLUVERTDYN
      PDHCV(JROF,JLEV,IDHCV+9) = ZUS2RG*( ZCPT(JROF,JLEV)+&
       & ZCPT(JROF,JLEV+1) )*ZEVELH(JROF,JLEV)  
!     FEPCONVERSIFL
      PDHCV(JROF,JLEV,IDHCV+10) = -PAPHI(JROF,JLEV)*ZT_PSDIV(JROF,JLEV)*ZUSRG
    ENDDO
  ENDDO
  IF(LHDQLN) THEN
    DO JLEV = 1, KFLEV-1
      DO JROF=KSTART,KPROF
!       FQLFLUVERTDYN
        PDHCV(JROF,JLEV,IDHCV+11) = ZUS2RG*( PQL(JROF,JLEV)+&
         & PQL(JROF,JLEV+1) )*ZEVELH(JROF,JLEV)  
!       FQIFLUVERTDYN
        PDHCV(JROF,JLEV,IDHCV+12) = ZUS2RG*( PQI(JROF,JLEV)+&
         & PQI(JROF,JLEV+1) )*ZEVELH(JROF,JLEV)  
!       FQRFLUVERTDYN
        PDHCV(JROF,JLEV,IDHCV+13) = ZUS2RG*( PQR(JROF,JLEV)+&
         & PQR(JROF,JLEV+1) )*ZEVELH(JROF,JLEV)
!       FQSFLUVERTDYN
        PDHCV(JROF,JLEV,IDHCV+14) = ZUS2RG*( PQS(JROF,JLEV)+&
         & PQS(JROF,JLEV+1) )*ZEVELH(JROF,JLEV)
      ENDDO
    ENDDO
  ENDIF
  ! TKE + TTE
  IF (LPTKE) THEN
    DO JLEV=1,KFLEV-1
      DO JROF=KSTART,KPROF
!       FTKFLUVERTDYN
        PDHCV(JROF,JLEV,IDHCV+15) = ZUS2RG*( PGFL(JROF,JLEV,IPTR_TKE)+&
         & PGFL(JROF,JLEV+1,IPTR_TKE) )*ZEVELH(JROF,JLEV)
      ENDDO
    ENDDO
    IF (LCOEFK_PTTE) THEN
      DO JLEV=1,KFLEV-1
        DO JROF=KSTART,KPROF
!         FTTFLUVERTDYN
          PDHCV(JROF,JLEV,IDHCV+16) = ZUS2RG*( PGFL(JROF,JLEV,IPTR_TTE)+&
           & PGFL(JROF,JLEV+1,IPTR_TTE) )*ZEVELH(JROF,JLEV)
        ENDDO
      ENDDO
    ENDIF
  ENDIF

! * CONDITION A LA LIMITE SUPERIEURE (SEMI-IMPOSEE POUR GAGNER DU TEMPS)
!   TOP BOUNDARY CONDITION

!   Hypotheses faites ici / Assumptions for top values:
!   (U,V)_[top]=(U,V)_[l=1]
!   T_[top]=T_[l=1]; q_[top]=q_[l=1]; ql_[top]=ql_[l=1]; qi_[top]=qi_[l=1]

  DO JROF=KSTART,KPROF
!   FPPFLUVERTDYN
    PDHCV(JROF,0,IDHCV+1) = ZEVELH(JROF,0)*ZUSRG
!   AUTRES FLUX VERTICAUX AU SOMMET
!   FQVFLUVERTDYN
    PDHCV(JROF,0,IDHCV+2) = 0.0_JPRB
!   FUUFLUVERTDYN
    PDHCV(JROF,0,IDHCV+3) = 0.0_JPRB
!   FVVFLUVERTDYN
    PDHCV(JROF,0,IDHCV+4) = 0.0_JPRB
!   FUUFLUDUAPLUI
    PDHCV(JROF,0,IDHCV+5) = 0.0_JPRB
!   FVVFLUDUAPLUI
    PDHCV(JROF,0,IDHCV+6) = 0.0_JPRB
!   FKKFLUVERTDYN
    PDHCV(JROF,0,IDHCV+7) = 0.0_JPRB
!   FKKFLUDUAPLUI
    PDHCV(JROF,0,IDHCV+8) = 0.0_JPRB
!   FCTFLUVERTDYN
    PDHCV(JROF,0,IDHCV+9) = 0.0_JPRB
!   FEPCONVERSIFL
    PDHCV(JROF,0,IDHCV+10) = -PAPHI(JROF,0)*ZT_PSDIV(JROF,0)*ZUSRG
    IF(LHDQLN) THEN
!     FQLFLUVERTDYN
      PDHCV(JROF,0,IDHCV+11) = 0.0_JPRB
!     FQIFLUVERTDYN
      PDHCV(JROF,0,IDHCV+12) = 0.0_JPRB
!     FQRFLUVERTDYN
      PDHCV(JROF,0,IDHCV+13) = 0.0_JPRB
!     FQSFLUVERTDYN
      PDHCV(JROF,0,IDHCV+14) = 0.0_JPRB
    ENDIF
    ! TKE + TTE
    IF (LPTKE) THEN
!     FTKFLUVERTDYN    
      PDHCV(JROF,0,IDHCV+15) = 0.0_JPRB
      IF (LCOEFK_PTTE) THEN
!       FTTFLUVERTDYN
        PDHCV(JROF,0,IDHCV+16) = 0.0_JPRB
      ENDIF
    ENDIF
  ENDDO
ENDIF

! 3.2/
! -- Variables "X" = "M1", "M2", "M3".
!    Variables "X" = "M1", "M2", "M3".

IF ( LHDMCI ) THEN
  IDHCV=NDHVV+NDHAVD+NDHBHKD

! * Code not updated for NH model.

! * CONDITION A LA LIMITE INFERIEURE / BOTTOM BOUNDARY CONDITION

  DO JROF=KSTART,KPROF
!   FA1FLUVERTDYN
    PDHCV(JROF,KFLEV,IDHCV+1) = - ZUSRG*RA**2 *ROMEGA*&
     & YDGSGEOM%GEMU(JROF)*YDGSGEOM%GSQM2(JROF)*ZCLON(JROF)*ZT_EVEL(JROF,KFLEV)  
!   FA2FLUVERTDYN
    PDHCV(JROF,KFLEV,IDHCV+2) = - ZUSRG*RA**2 *ROMEGA*&
     & YDGSGEOM%GEMU(JROF)*YDGSGEOM%GSQM2(JROF)*ZSLON(JROF)*ZT_EVEL(JROF,KFLEV)  
!   FA3FLUVERTDYN
    PDHCV(JROF,KFLEV,IDHCV+3) =   ZUSRG*RA**2 *ROMEGA*&
     & YDGSGEOM%GSQM2(JROF)*YDGSGEOM%GSQM2(JROF)*ZT_EVEL(JROF,KFLEV)  
!   FA1GRAV
!   FA2GRAV
!   FA3GRAV
    ZORZON = YDGSGEOM%GNORDM(JROF)*PORL(JROF)-YDGSGEOM%GNORDL(JROF)*PORM(JROF)
    ZORMER = YDGSGEOM%GNORDL(JROF)*PORL(JROF)+YDGSGEOM%GNORDM(JROF)*PORM(JROF)
    ZPHILZ(JROF) = ZORZON
    ZPHILM(JROF) = ZORMER
    PDHCV(JROF,KFLEV,IDHCV+4) = ZUSRG*RA*PAPRS(JROF,KFLEV)*&
     & ( ZORMER*ZSLON(JROF)-ZORZON*ZCLON(JROF)*YDGSGEOM%GEMU(JROF))  
    PDHCV(JROF,KFLEV,IDHCV+5) = ZUSRG*RA*PAPRS(JROF,KFLEV)*&
     & (-ZORMER*ZCLON(JROF)-ZORZON*ZSLON(JROF)*YDGSGEOM%GEMU(JROF))  
    PDHCV(JROF,KFLEV,IDHCV+6) = ZUSRG*RA*PAPRS(JROF,KFLEV)*ZORZON*YDGSGEOM%GSQM2(JROF)
!   FA1FLUDUAPLUI
!   FA2FLUDUAPLUI
!   FA3FLUDUAPLUI
    PDHCV(JROF,KFLEV,IDHCV+7) = -RA**2 *ROMEGA*ZFPSUM(JROF,KFLEV)*&
     & YDGSGEOM%GEMU(JROF)*YDGSGEOM%GSQM2(JROF)*ZCLON(JROF)  
    PDHCV(JROF,KFLEV,IDHCV+8) = -RA**2 *ROMEGA*ZFPSUM(JROF,KFLEV)*&
     & YDGSGEOM%GEMU(JROF)*YDGSGEOM%GSQM2(JROF)*ZSLON(JROF)  
    PDHCV(JROF,KFLEV,IDHCV+9) =  RA**2 *ROMEGA*ZFPSUM(JROF,KFLEV)*&
     & YDGSGEOM%GSQM2(JROF)*YDGSGEOM%GSQM2(JROF)  
  ENDDO

! * INTERNIVEAUX NORMAUX / HALF LEVELS

  DO JLEV = KFLEV-1, 1, -1
    DO JROF=KSTART,KPROF
      ZUS2RG = 0.5_JPRB*ZUSRG
!     FA1FLUVERTDYN
!     FA2FLUVERTDYN
!     FA3FLUVERTDYN
      PDHCV(JROF,JLEV,IDHCV+1) = ZUS2RG*ZT_EVEL(JROF,JLEV)*&
       & (ZM1(JROF,JLEV+1)+ZM1(JROF,JLEV))  
      PDHCV(JROF,JLEV,IDHCV+2) = ZUS2RG*ZT_EVEL(JROF,JLEV)*&
       & (ZM2(JROF,JLEV+1)+ZM2(JROF,JLEV))  
      PDHCV(JROF,JLEV,IDHCV+3) = ZUS2RG*ZT_EVEL(JROF,JLEV)*&
       & (ZM3(JROF,JLEV+1)+ZM3(JROF,JLEV))  
!     FA1GRAV
!     FA2GRAV
!     FA3GRAV
      ZPHILZ(JROF) = ZPHILZ(JROF) +ZRTZON(JROF,JLEV+1)*PXYB(JROF,JLEV+1,YYTXYB0%M_LNPR)
      ZPHILM(JROF) = ZPHILM(JROF) +ZRTMER(JROF,JLEV+1)*PXYB(JROF,JLEV+1,YYTXYB0%M_LNPR)
      PDHCV(JROF,JLEV,IDHCV+4) = ZUSRG*RA*PAPRS(JROF,JLEV)*&
       & (ZPHILM(JROF)*ZSLON(JROF)-ZPHILZ(JROF)*ZCLON(JROF)*&
       & YDGSGEOM%GEMU(JROF))  
      PDHCV(JROF,JLEV,IDHCV+5) = ZUSRG*RA*PAPRS(JROF,JLEV)*&
       & (-ZPHILM(JROF)*ZCLON(JROF)-ZPHILZ(JROF)*ZSLON(JROF)*&
       & YDGSGEOM%GEMU(JROF))  
      PDHCV(JROF,JLEV,IDHCV+6) = ZUSRG*RA*PAPRS(JROF,JLEV)*&
       & ZPHILZ(JROF)*YDGSGEOM%GSQM2(JROF)  
!     FA1FLUDUAPLUI
!     FA2FLUDUAPLUI
!     FA3FLUDUAPLUI
      PDHCV(JROF,JLEV,IDHCV+7) = 0.5_JPRB*ZFPSUM(JROF,JLEV)*&
       & (ZM1(JROF,JLEV)+ZM1(JROF,JLEV+1))  
      PDHCV(JROF,JLEV,IDHCV+8) = 0.5_JPRB*ZFPSUM(JROF,JLEV)*&
       & (ZM2(JROF,JLEV)+ZM2(JROF,JLEV+1))  
      PDHCV(JROF,JLEV,IDHCV+9) = 0.5_JPRB*ZFPSUM(JROF,JLEV)*&
       & (ZM3(JROF,JLEV)+ZM3(JROF,JLEV+1))  
    ENDDO
  ENDDO

! * CONDITION A LA LIMITE SUPERIEURE (SEMI-IMPOSEE POUR GAGNER DU TEMPS)
!   TOP BOUNDARY CONDITION

!   Hypotheses faites ici / Assumptions for top values:
!   (M1,M2,M3)_[top]=(M1,M2,M3)_[l=1]

  DO JROF=KSTART,KPROF
    PDHCV(JROF,0,IDHCV+1) = 0.0_JPRB
    PDHCV(JROF,0,IDHCV+2) = 0.0_JPRB
    PDHCV(JROF,0,IDHCV+3) = 0.0_JPRB
    PDHCV(JROF,0,IDHCV+4) = 0.0_JPRB
    PDHCV(JROF,0,IDHCV+5) = 0.0_JPRB
    PDHCV(JROF,0,IDHCV+6) = 0.0_JPRB
    PDHCV(JROF,0,IDHCV+7) = 0.0_JPRB
    PDHCV(JROF,0,IDHCV+8) = 0.0_JPRB
    PDHCV(JROF,0,IDHCV+9) = 0.0_JPRB
  ENDDO
ENDIF

! 3.3/
! -- Variable "X" = "S".

IF ( LHDENT ) THEN
  IDHCV=NDHVV+NDHAVD+NDHBHKD+NDHBMCD

! * CONDITION A LA LIMITE INFERIEURE / BOTTOM BOUNDARY CONDITION

!   Hypotheses faites ici / Assumptions for bottom values:
!   S_[bottom]=S_[surf]=Sd_[surf]+(Sv_[surf]-Sd_[surf])*q_[surf]

  DO JROF=KSTART,KPROF
!   FSSFLUVERTDYN
    ZPREVS = MAX(ZQTHRES,PQVS(JROF))&
     & *RV*ZAPRS(JROF,KFLEV)/(RD+(RV-RD)*PQVS(JROF))  
    ZENTRAS = RCPD*LOG(PTSFC(JROF))-RD*LOG(ZAPRS(JROF,KFLEV)-ZPREVS)+ZENTRA0
    ZENTRVS = RCPV*LOG(PTSFC(JROF))-RV*LOG(ZPREVS)+ZENTRV0
    PDHCV(JROF,KFLEV,IDHCV+1) = ZUSRG*ZEVELH(JROF,KFLEV)*&
     & (ZENTRAS+(ZENTRVS-ZENTRAS)*PQVS(JROF))  
!   FSSPRECISS
    PDHCV(JROF,KFLEV,IDHCV+2) =&
     & ZENTRL(JROF,KFLEV)*(PFPLCL(JROF,KFLEV)+PFPLSL(JROF,KFLEV)) +&
     & ZENTRS(JROF,KFLEV)*(PFPLCN(JROF,KFLEV)+PFPLSN(JROF,KFLEV))  
  ENDDO

! * INTERNIVEAUX NORMAUX / HALF LEVELS

  ZUS2RG = 0.5_JPRB*ZUSRG
  DO JLEV = 1, KFLEV-1
    DO JROF=KSTART,KPROF
!     FSSFLUVERTDYN
      PDHCV(JROF,JLEV,IDHCV+1) = ZUS2RG*(ZENTR(JROF,JLEV)+&
       & ZENTR(JROF,JLEV+1))*ZEVELH(JROF,JLEV)  
!     FSSPRECISS
      PDHCV(JROF,JLEV,IDHCV+2) =&
       & ZENTRL(JROF,JLEV)*(PFPLCL(JROF,JLEV)+PFPLSL(JROF,JLEV)) +&
       & ZENTRS(JROF,JLEV)*(PFPLCN(JROF,JLEV)+PFPLSN(JROF,JLEV))  
    ENDDO
  ENDDO

! * CONDITION A LA LIMITE SUPERIEURE (SEMI-IMPOSEE POUR GAGNER DU TEMPS)
!   TOP BOUNDARY CONDITION

!   Hypotheses faites ici / Assumptions for top values:
!   S_[top]=S_[l=1]

  DO JROF=KSTART,KPROF
!   FSSFLUVERTDYN
    PDHCV(JROF,0,IDHCV+1) = 0.0_JPRB
!   FSSPRECISS
    PDHCV(JROF,0,IDHCV+2) = 0.0_JPRB
  ENDDO
ENDIF

!     ------------------------------------------------------------------
!*    4.- TENDANCES DYNAMIQUES DU SEMI-IMPLICITES/
!         SEMI-IMPLICIT DYNAMICAL TENDENCIES

IF (LRSIDDH) THEN
  IF (LFLEXDIA) THEN
    IF (LDDH_OMP) THEN
      ZTMPAF(KSTART:KPROF,1:KFLEV)=ZDELPREV(KSTART:KPROF,1:KFLEV)*(&
      & PGMVTNDSI(KSTART:KPROF,1:KFLEV,MSIDDH_U1)+PGMVTNDSI(KSTART:KPROF,1:KFLEV,MSIDDH_U0))
      CALL NEW_ADD_FIELD_3D(YDML_DIAG%YRMDDH,ZTMPAF,'TUUSI',YDDDH)
      ZTMPAF(KSTART:KPROF,1:KFLEV)=ZDELPREV(KSTART:KPROF,1:KFLEV)*(&
      & PGMVTNDSI(KSTART:KPROF,1:KFLEV,MSIDDH_V1)+PGMVTNDSI(KSTART:KPROF,1:KFLEV,MSIDDH_V0))
      CALL NEW_ADD_FIELD_3D(YDML_DIAG%YRMDDH,ZTMPAF,'TVVSI',YDDDH)
      IF (NFTHER > 0) THEN
        ZTMPAF(KSTART:KPROF,1:KFLEV)=ZDELPREV(KSTART:KPROF,1:KFLEV)*(&
        & PGMVTNDSI(KSTART:KPROF,1:KFLEV,MSIDDH_T1)+PGMVTNDSI(KSTART:KPROF,1:KFLEV,MSIDDH_T0))
        CALL NEW_ADD_FIELD_3D(YDML_DIAG%YRMDDH,ZTMPAF,'TTTSI',YDDDH)
      ENDIF
      IF (YDDYNA%LNHDYN) THEN
        ZTMPAF(KSTART:KPROF,1:KFLEV)=ZDELPREV(KSTART:KPROF,1:KFLEV)*(&
        &PGMVTNDSI(KSTART:KPROF,1:KFLEV,MSIDDH_PD1)+PGMVTNDSI(KSTART:KPROF,1:KFLEV,MSIDDH_PD0))
        CALL NEW_ADD_FIELD_3D(YDML_DIAG%YRMDDH,ZTMPAF,'TPDSI',YDDDH)
        ZTMPAF(KSTART:KPROF,1:KFLEV)=ZDELPREV(KSTART:KPROF,1:KFLEV)*(&
        & PGMVTNDSI(KSTART:KPROF,1:KFLEV,MSIDDH_VD1)+PGMVTNDSI(KSTART:KPROF,1:KFLEV,MSIDDH_VD0))
        CALL NEW_ADD_FIELD_3D(YDML_DIAG%YRMDDH,ZTMPAF,'TVDSI',YDDDH)
      ENDIF
    ELSE
      ZTMPAF(KSTART:KPROF,1:KFLEV)=ZDELPREV(KSTART:KPROF,1:KFLEV)*(&
      & PGMVTNDSI(KSTART:KPROF,1:KFLEV,MSIDDH_U1)+PGMVTNDSI(KSTART:KPROF,1:KFLEV,MSIDDH_U0))
      CALL ADD_FIELD_3D(YDML_DIAG%YRLDDH,ZTMPAF,'TUUSI','T','ARP',.TRUE.,.TRUE.)
      ZTMPAF(KSTART:KPROF,1:KFLEV)=ZDELPREV(KSTART:KPROF,1:KFLEV)*(&
      & PGMVTNDSI(KSTART:KPROF,1:KFLEV,MSIDDH_V1)+PGMVTNDSI(KSTART:KPROF,1:KFLEV,MSIDDH_V0))
      CALL ADD_FIELD_3D(YDML_DIAG%YRLDDH,ZTMPAF,'TVVSI','T','ARP',.TRUE.,.TRUE.)
      IF (NFTHER > 0) THEN
        ZTMPAF(KSTART:KPROF,1:KFLEV)=ZDELPREV(KSTART:KPROF,1:KFLEV)*(&
        & PGMVTNDSI(KSTART:KPROF,1:KFLEV,MSIDDH_T1)+PGMVTNDSI(KSTART:KPROF,1:KFLEV,MSIDDH_T0))
        CALL ADD_FIELD_3D(YDML_DIAG%YRLDDH,ZTMPAF,'TTTSI','T','ARP',.TRUE.,.TRUE.)
      ENDIF
      IF (YDDYNA%LNHDYN) THEN
        ZTMPAF(KSTART:KPROF,1:KFLEV)=ZDELPREV(KSTART:KPROF,1:KFLEV)*(&
        & PGMVTNDSI(KSTART:KPROF,1:KFLEV,MSIDDH_PD1)+PGMVTNDSI(KSTART:KPROF,1:KFLEV,MSIDDH_PD0))
        CALL ADD_FIELD_3D(YDML_DIAG%YRLDDH,ZTMPAF,'TPDSI','T','ARP',.TRUE.,.TRUE.)
        ZTMPAF(KSTART:KPROF,1:KFLEV)=ZDELPREV(KSTART:KPROF,1:KFLEV)*(&
        & PGMVTNDSI(KSTART:KPROF,1:KFLEV,MSIDDH_VD1)+PGMVTNDSI(KSTART:KPROF,1:KFLEV,MSIDDH_VD0))
        CALL ADD_FIELD_3D(YDML_DIAG%YRLDDH,ZTMPAF,'TVDSI','T','ARP',.TRUE.,.TRUE.)
      ENDIF
    ENDIF
  ENDIF
ENDIF


!     ------------------------------------------------------------------
!*    5.- TENDANCES DYNAMIQUES DE LA DIFFUSION HORIZONTALE/
!         DYNAMICAL TENDENCIES OF HORIZONTAL DIFFUSION

IF (LRHDDDH) THEN
  IF (LFLEXDIA) THEN
    IF (LDDH_OMP) THEN
      ZTMPAF(KSTART:KPROF,1:KFLEV)=ZDELPREV(KSTART:KPROF,1:KFLEV)*&
      & PGMVTNDHD(KSTART:KPROF,1:KFLEV,MHDDDH_U)
      CALL NEW_ADD_FIELD_3D(YDML_DIAG%YRMDDH,ZTMPAF,'TUUHD',YDDDH)
      ZTMPAF(KSTART:KPROF,1:KFLEV)=ZDELPREV(KSTART:KPROF,1:KFLEV)*&
      & PGMVTNDHD(KSTART:KPROF,1:KFLEV,MHDDDH_V)
      CALL NEW_ADD_FIELD_3D(YDML_DIAG%YRMDDH,ZTMPAF,'TVVHD',YDDDH)
      IF (NFTHER > 0) THEN
        ZTMPAF(KSTART:KPROF,1:KFLEV)=ZDELPREV(KSTART:KPROF,1:KFLEV)*&
        & PGMVTNDHD(KSTART:KPROF,1:KFLEV,MHDDDH_T)
        CALL NEW_ADD_FIELD_3D(YDML_DIAG%YRMDDH,ZTMPAF,'TTTHD',YDDDH)
      ENDIF
      IF (YDDYNA%LNHDYN) THEN
        ZTMPAF(KSTART:KPROF,1:KFLEV)=ZDELPREV(KSTART:KPROF,1:KFLEV)*&
        & PGMVTNDHD(KSTART:KPROF,1:KFLEV,MHDDDH_PD)
        CALL NEW_ADD_FIELD_3D(YDML_DIAG%YRMDDH,ZTMPAF,'TPDHD',YDDDH)
        ZTMPAF(KSTART:KPROF,1:KFLEV)=ZDELPREV(KSTART:KPROF,1:KFLEV)*&
        & PGMVTNDHD(KSTART:KPROF,1:KFLEV,MHDDDH_VD)
        CALL NEW_ADD_FIELD_3D(YDML_DIAG%YRMDDH,ZTMPAF,'TVDHD',YDDDH)
      ENDIF
      IF (YDDYNA%LNHX) THEN
        ZTMPAF(KSTART:KPROF,1:KFLEV)=ZDELPREV(KSTART:KPROF,1:KFLEV)*&
        & PGMVTNDHD(KSTART:KPROF,1:KFLEV,MHDDDH_NHX)
        CALL NEW_ADD_FIELD_3D(YDML_DIAG%YRMDDH,ZTMPAF,'TNHXHD',YDDDH)
      ENDIF
      ZTMPAF(KSTART:KPROF,1:KFLEV)=ZDELPREV(KSTART:KPROF,1:KFLEV)*&
      & PGFLTNDHD(KSTART:KPROF,1:KFLEV,MHDDDH_Q)
      CALL NEW_ADD_FIELD_3D(YDML_DIAG%YRMDDH,ZTMPAF,'TQVHD',YDDDH)
    ELSE
      ZTMPAF(KSTART:KPROF,1:KFLEV)=ZDELPREV(KSTART:KPROF,1:KFLEV)*&
      & PGMVTNDHD(KSTART:KPROF,1:KFLEV,MHDDDH_U)
      CALL ADD_FIELD_3D(YDML_DIAG%YRLDDH,ZTMPAF,'TUUHD','T','ARP',.TRUE.,.TRUE.)
      ZTMPAF(KSTART:KPROF,1:KFLEV)=ZDELPREV(KSTART:KPROF,1:KFLEV)*&
      & PGMVTNDHD(KSTART:KPROF,1:KFLEV,MHDDDH_V)
      CALL ADD_FIELD_3D(YDML_DIAG%YRLDDH,ZTMPAF,'TVVHD','T','ARP',.TRUE.,.TRUE.)
      IF (NFTHER > 0) THEN
        ZTMPAF(KSTART:KPROF,1:KFLEV)=ZDELPREV(KSTART:KPROF,1:KFLEV)*&
        & PGMVTNDHD(KSTART:KPROF,1:KFLEV,MHDDDH_T)
        CALL ADD_FIELD_3D(YDML_DIAG%YRLDDH,ZTMPAF,'TTTHD','T','ARP',.TRUE.,.TRUE.)
      ENDIF
      IF (YDDYNA%LNHDYN) THEN
        ZTMPAF(KSTART:KPROF,1:KFLEV)=ZDELPREV(KSTART:KPROF,1:KFLEV)*&
        & PGMVTNDHD(KSTART:KPROF,1:KFLEV,MHDDDH_PD)
        CALL ADD_FIELD_3D(YDML_DIAG%YRLDDH,ZTMPAF,'TPDHD','T','ARP',.TRUE.,.TRUE.)
        ZTMPAF(KSTART:KPROF,1:KFLEV)=ZDELPREV(KSTART:KPROF,1:KFLEV)*&
        & PGMVTNDHD(KSTART:KPROF,1:KFLEV,MHDDDH_VD)
        CALL ADD_FIELD_3D(YDML_DIAG%YRLDDH,ZTMPAF,'TVDHD','T','ARP',.TRUE.,.TRUE.)
      ENDIF
      IF (YDDYNA%LNHX) THEN
        ZTMPAF(KSTART:KPROF,1:KFLEV)=ZDELPREV(KSTART:KPROF,1:KFLEV)*&
        & PGMVTNDHD(KSTART:KPROF,1:KFLEV,MHDDDH_NHX)
        CALL ADD_FIELD_3D(YDML_DIAG%YRLDDH,ZTMPAF,'TNHXHD','T','ARP',.TRUE.,.TRUE.)
      ENDIF
      ZTMPAF(KSTART:KPROF,1:KFLEV)=ZDELPREV(KSTART:KPROF,1:KFLEV)*&
      & PGFLTNDHD(KSTART:KPROF,1:KFLEV,MHDDDH_Q)
      CALL ADD_FIELD_3D(YDML_DIAG%YRLDDH,ZTMPAF,'TQVHD','T','ARP',.TRUE.,.TRUE.)
    ENDIF
  ENDIF
ENDIF


!     ------------------------------------------------------------------

IF (LRSLDDH) THEN
  IF (.NOT. LDDH_OMP) THEN
    ZTMPAF(KSTART:KPROF,1:KFLEV)= 0.0_JPRB
    CALL ADD_FIELD_3D(YDML_DIAG%YRLDDH,ZTMPAF,'TSLUADV','T','ARP',.TRUE.,.TRUE.,LDSLFD=.TRUE.)
    CALL ADD_FIELD_3D(YDML_DIAG%YRLDDH,ZTMPAF,'TSLVADV','T','ARP',.TRUE.,.TRUE.,LDSLFD=.TRUE.)
    IF (NFTHER > 0) THEN
      CALL ADD_FIELD_3D(YDML_DIAG%YRLDDH,ZTMPAF,'TSLTADV','T','ARP',.TRUE.,.TRUE.,LDSLFD=.TRUE.)
    ENDIF
    IF (YDDYNA%LNHDYN) THEN
      CALL ADD_FIELD_3D(YDML_DIAG%YRLDDH,ZTMPAF,'TSLPDADV','T','ARP',.TRUE.,.TRUE.,LDSLFD=.TRUE.)
      CALL ADD_FIELD_3D(YDML_DIAG%YRLDDH,ZTMPAF,'TSLVDADV','T','ARP',.TRUE.,.TRUE.,LDSLFD=.TRUE.)
    ENDIF
  ENDIF
ENDIF

!     ------------------------------------------------------------------
!*    6.- VERIFICATION/MAINTENANCE EVENTUELLE
!         ADDITIONAL CHECK

IF ( LHDLIST ) THEN
  WRITE (NULOUT,*) ' DDH * CPDYDDH IDHCV=',IDHCV
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CPDYDDH',1,ZHOOK_HANDLE)
END SUBROUTINE CPDYDDH
