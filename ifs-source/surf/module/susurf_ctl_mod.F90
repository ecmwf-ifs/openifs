! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction
MODULE SUSURF_CTL_MOD
CONTAINS
SUBROUTINE SUSURF_CTL(KSW,KCSS,KCSNEC,KSIL,KCOM,KTILES,KTSW,KLWEMISS,&
 & LD_LEFLAKE, LD_LEOCML,LD_LOCMLTKE,&
 & LD_LWCOU, LD_LWCOU2W, LD_LWCOUHMF,&
 & LD_LLCCNL,LD_LLCCNO,LD_LEVGEN,LD_LESSRO,LD_LELAIV,&
 & LD_LECTESSEL,LD_LEAGS,LD_LEFARQUHAR, LD_LEAIRCO2COUP,LD_LESN09,LD_LESNML,&
 & LD_LEOCWA,LD_LEOCCO,LD_LEOCSA,LD_LEOCLA,KALBEDOSCHEME,KEMISSSCHEME,&
 & LD_LSCMEC,LD_LROUGH,PEXTZ0M,PEXTZ0H,&
 & PTHRFRTI,PTSTAND,PXP,PRCCNSEA,PRCCNLND,&
 & PRLAIINT,PRSUN,PRCORIOI,PRPLRG,PNSNMLWS,&
 & PRVR0VT,PRVCMAX25,PRHUMREL,PRA1,PRB1,PRG0,PRGM25,PRE_VCMAX,PRE_JMAX,&
 & YDDIM,YDEXC,YDCST,YDRAD,YDRDI,YDLW,YDSW,YDSOIL,&
 & YDVEG,YDAGS,YDAGF,YDMLM,YDFLAKE,YDOCEAN_ML,YDURB,PRALFMINPSN,PRCIMIN)

USE PARKIND1,     ONLY : JPIM, JPRB
USE YOMHOOK,      ONLY : LHOOK, DR_HOOK, JPHOOK

!ifndef INTERFACE

USE YOS_DIM,      ONLY : TDIM
USE YOS_EXC,      ONLY : TEXC
USE YOS_CST,      ONLY : TCST
USE YOS_RAD,      ONLY : TRAD
USE YOS_RDI,      ONLY : TRDI
USE YOS_LW,       ONLY : TLW
USE YOS_SW,       ONLY : TSW
USE YOS_SOIL,     ONLY : TSOIL
USE YOS_VEG,      ONLY : TVEG
USE YOS_AGS,      ONLY : TAGS
USE YOS_AGF,      ONLY : TAGF
USE YOS_MLM,      ONLY : TMLM
USE YOS_FLAKE,    ONLY : TFLAKE
USE YOS_OCEAN_ML, ONLY : TOCEAN_ML
USE YOS_URB,      ONLY : TURB

USE SUSCST_MOD
USE SUSTHF_MOD
USE SUSRAD_MOD
USE SUSSOIL_MOD
USE SUSVEG_MOD
USE SUCOTWO_MOD
USE SUFARQUHAR_MOD
USE SUVEXC_MOD
USE SUVEXCS_MOD
USE SUSFLAKE_MOD



!OIFS USE SUSOCEAN_ML_MOD

USE SUGRIDMLM_MOD
USE SUSURB_MOD

!**   *SUSURF* IS THE SET-UP ROUTINE FOR COMMON BLOCK *YOESOIL*

!     PURPOSE
!     -------
!          THIS ROUTINE INITIALIZES THE CONSTANTS IN COMMON BLOCK
!     *YOESOIL*

!     INTERFACE.
!     ----------
!     CALL *SUSURF* FROM *SUPHEC*

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!     SUSCST     Setup general constants
!     SUSTHF     Setup thermodynamic function constants
!     SUSRAD     Setup radiation constants
!     SUSSOIL    Setup soil constants
!     SUSVEG     Setup vegetation constants
!     SUVEXC     Setup surface exchange coefficients constants
!     SUVEXCS    Setup static stability constants

!     REFERENCE.
!     ----------

!     Original    A.C.M. BELJAARS         E.C.M.W.F.      89/11/02

!     MODIFICATIONS
!     -------------
!     J.-J. MORCRETTE         E.C.M.W.F.      91/07/14
!     P. VITERBO              E.C.M.W.F.       8/10/93
!     P. Viterbo     99-03-26    Tiling of the land surface
!     C. Fischer 00-12-20 Meteo-France recode initialization of rdat to avoid
!                         memory overflow on SUN workstation
!     J.F. Estrade *ECMWF* 03-10-01 move in surf vob
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!     P. Viterbo    24-05-2004      Change surface units
!     P. Viterbo   ECMWF   03-12-2004  Include user-defined RTHRFRTI
!     P. Viterbo   ECMWF   May 2005    Externalise surf
!     JJMorcrette 20060511 MODIS albedo
!     G. Balsamo   ECMWF   10-01-2006  Include Vangenucthen Hydro.
!     G. Balsamo   ECMWF   11-01-2006  Include sub-grid surface runoff
!     V. Stepanenko/G. Balsamo 01-05-2008 add lake tile
!     Y. Takaya    ECMWF   07-10-2008  Include setup for ocean mixed layer model
!     G. Balsamo   ECMWF   13-10-2008  Include switch for liquid water in snow
!     E. Dutra             12-11.2008  Include new snow parameterization
!     Y. Takaya    ECMWF   21-08-2009  Include Langmuir effect to skin layer model
!     S. Boussetta/G.Balsamo May 2009 Add lai
!     E. Dutra             16-11-2009  snow 2009 cleaning
!     S. Boussetta/G.Balsamo May 2010 Add CTESSEL for CO2
!     P. Bechtold          26-03-2012  Add small planet PRCORIOI PRPLRG
!     R. Hogan             14-01-2019  Replace LE4ALB with KALBEDOSCHEME; add KEMISSSCHEME, PALFMINPSN
!     A. Agusti-Panareda Nov 2020 couple atm CO2 tracer (LEAIRCO2COUP) with photosynthesis
!     V.Bastrikov,F.Maignan,P.Peylin,A.Agusti-Panareda/S.Boussetta Feb 2021 Add Farquhar photosynthesis model
!     A. Agusti-Panareda June 2021 Pass optimized photosynthesis parameters from namelist yoephy to sufarquhar
!
!  INTERFACE:

!    Integers (In):

!      KSW       : NUMBER OF SHORTWAVE SPECTRAL INTERVALS
!      KCSS      : Number of soil levels
!      KCSNEC    : Number of snow levels
!      KSIL      : NUMBER OF (infrared) SPECTRAL INTERVALS
!      KCOM      : Number of layers in mixed layer model
!      KTILES    : Number of surface tiles
!      KTSW      : Maximum possible number of sw spectral intervals
!      KLWEMISS  : Number of longwave emissivity spectral intervals
!      KALBEDOSCHEME : (0)ERBE,(1)4-comp MODIS,(2)6-comp MODIS,(3)2-comp MODIS (=4comp diffuse)
!      KEMISSSCHEME  : (0) 2-value emissivity scheme, (1) 6-value scheme

!    Logicals (In):

!      LD_LLCCNL : .T. IF CCN CONCENTRATION OVER LAND IS DIAGNOSED
!      LD_LLCCNO : .T. IF CCN CONCENTRATION OVER OCEAN IS DIAGNOSED
!      LD_LEVGEN : .T. IF VAN GENUCHTEN HYDRO IS ACTIVATED
!      LD_LESSRO : .T. IF SUB-GRID SURFACE RUNOFF IS ACTIVATED
!      LD_LELAIV : .T. IF LAI FROM CLIMATE FIELDS IS USED

!      LD_LECTESSEL : .T. IF CTESSEL schemme is used
!      LD_LEAGS     : .T. IF AGS canopy resistance is used else Jravis type
!      LD_LEFARQUHAR: .T. IF Farquhar model is used for photosynthesis
!      LD_LEAIRCO2COUP : .T. if variable CO2 is used in photosynthesis else a fixed global annual value is used
!      LD_LESN09 : .T. IF SNOW 2009 IS ACTIVATED
!      LD_LEOCWA : .T. if WARM OCEAN LAYER PARAMETRIZATION active
!      LD_LEOCCO : .T. if COOL OCEAN SKIN PARAMETRIZATION active
!      LD_LEOCSA : .T. if SALINTY EFFECT ON SATURATION AT OCEAN SURFACE active
!      LD_LEOCLA : .T. if LANGMUIR CIRCULATION EFFECT IN VOSKIN active
!      LD_LEFLAKE: .T. IF FLAKE PARAMETRIZATION USED FOR LAKES
!      LD_LEOCML : .T. IF OCEAN MIXED LAYER MODEL (KPP) ACTIVATE
!      LD_LOCMLTKE : .T. IF OCEAN MIXED LAYER MODEL WITH TKE SCHEME ACTIVATE

!      LD_LWCOU    : .T. IF COUPLED TO WAVE MODEL
!      LD_LWCOU2W  : .T. IF COUPLED TO WAVE MODEL WITH FEEDBACK TO ATMOSPHERE
!      LD_LWCOUHMF : .T. IF SEA STATE DEPENDENT HEAT AND MOISTURE FLUXES IF COUPLED TO WAVE MODEL

!    Reals (In):

!      PTHRFRTI  : ! MINIMUM THRESHOLD FOR TILE FRACTION
!      PTSTAND   : ! REFERENCE TEMPERATURE FOR TEMPERATURE DEPENDENCE
!      PXP       : ! POLYNOMIAL COEFFICIENTS OF PLANCK FUNCTION
!      PRCCNSEA  : ! NUMBER CONCENTRATION (CM-3) OF CCNs OVER SEA
!      PRCCNLND  : ! NUMBER CONCENTRATION (CM-3) OF CCNs OVER LAND
!      PRLAIINT  : ! INTERACTIVE LAI COEFFICIENT (1=Interactive ; 0=climatology)
!      PRSUN     : ! SOLAR FRACTION IN SPECTRAL INTERVALS
!      PRCORIOI  : ! REDUCTION FACTOR IN LENGTH OF DAY FOR SMALL PLANET
!      PRPLRG    : ! FACTOR FOR GRAVITY
!      PRALFMINPSN:! Albedo of permanent snow
!      PRCIMIN   : ! MINIMUM ICE FRACTION

!      * Farquhar photosynthesis model parameters optimized with observations *
!
!      PRVR0VT : reference ecosystem respiration [Kg CO2 m-2 s-1]
!      RVCMAX25 : MAX RATE of Rubisco activity-limited carboxylation at 25°C
!      RHUMREL : SCALING FACTOR for soil moisture stress (optimized per PFT with FLUXNET data)
!      RA1: Empirical factor in the calculation of fvpd (-)
!      RB1: Empirical factor in the calculation of fvpd (-)
!      RG0: Residual stomatal conductance when irradiance approaches zero (mol CO2 m−2 s−1 bar−1)
!      RGM25 Mesophyll diffusion conductance at 25oC (mol m-2 s-1 bar-1)
!      RE_VCMAX: Energy of activation for Vcmax (J mol-1)
!      RE_JMAX: Energy of activation for Jmax (J mol-1)
!     ------------------------------------------------------------------


IMPLICIT NONE

! Declaration of arguments

INTEGER(KIND=JPIM),INTENT(IN)    :: KSW
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTHRFRTI
INTEGER(KIND=JPIM),INTENT(IN)    :: KCSS
INTEGER(KIND=JPIM),INTENT(IN)    :: KCSNEC
INTEGER(KIND=JPIM),INTENT(IN)    :: KSIL
INTEGER(KIND=JPIM),INTENT(IN)    :: KCOM
INTEGER(KIND=JPIM),INTENT(IN)    :: KTILES
INTEGER(KIND=JPIM),INTENT(IN)    :: KTSW
INTEGER(KIND=JPIM),INTENT(IN)    :: KLWEMISS
INTEGER(KIND=JPIM),INTENT(IN)    :: KALBEDOSCHEME
INTEGER(KIND=JPIM),INTENT(IN)    :: KEMISSSCHEME
LOGICAL           ,INTENT(IN)    :: LD_LLCCNL
LOGICAL           ,INTENT(IN)    :: LD_LLCCNO
LOGICAL           ,INTENT(IN)    :: LD_LEVGEN
LOGICAL           ,INTENT(IN)    :: LD_LESSRO
LOGICAL           ,INTENT(IN)    :: LD_LELAIV
LOGICAL           ,INTENT(IN)    :: LD_LECTESSEL
LOGICAL           ,INTENT(IN)    :: LD_LEAGS
LOGICAL           ,INTENT(IN)    :: LD_LEFARQUHAR
LOGICAL           ,INTENT(IN)    :: LD_LEAIRCO2COUP
LOGICAL           ,INTENT(IN)    :: LD_LESN09
LOGICAL           ,INTENT(IN)    :: LD_LESNML
LOGICAL           ,INTENT(IN)    :: LD_LEOCWA
LOGICAL           ,INTENT(IN)    :: LD_LEOCCO
LOGICAL           ,INTENT(IN)    :: LD_LEOCSA
LOGICAL           ,INTENT(IN)    :: LD_LEOCLA
LOGICAL           ,INTENT(IN)    :: LD_LSCMEC
LOGICAL           ,INTENT(IN)    :: LD_LROUGH
LOGICAL           ,INTENT(IN)    :: LD_LEFLAKE
LOGICAL           ,INTENT(IN)    :: LD_LEOCML
LOGICAL           ,INTENT(IN)    :: LD_LOCMLTKE
LOGICAL           ,INTENT(IN)    :: LD_LWCOU
LOGICAL           ,INTENT(IN)    :: LD_LWCOU2W
LOGICAL           ,INTENT(IN)    :: LD_LWCOUHMF
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEXTZ0M
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEXTZ0H
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSTAND
REAL(KIND=JPRB)   ,INTENT(IN)    :: PXP(6,6)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRCCNSEA
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRCCNLND
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRLAIINT
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRSUN(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRCORIOI
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRPLRG
INTEGER(KIND=JPIM),INTENT(IN)    :: PNSNMLWS
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRVR0VT(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRVCMAX25(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRHUMREL(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRA1(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRB1(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRG0(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRGM25(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRE_VCMAX(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRE_JMAX(:,:)
TYPE(TDIM)        ,INTENT(INOUT) :: YDDIM
TYPE(TEXC)        ,INTENT(INOUT) :: YDEXC
TYPE(TCST)        ,INTENT(INOUT) :: YDCST
TYPE(TRAD)        ,INTENT(INOUT) :: YDRAD
TYPE(TRDI)        ,INTENT(INOUT) :: YDRDI
TYPE(TLW)         ,INTENT(INOUT) :: YDLW
TYPE(TSW)         ,INTENT(INOUT) :: YDSW
TYPE(TSOIL)       ,INTENT(INOUT) :: YDSOIL
TYPE(TVEG)        ,INTENT(INOUT) :: YDVEG
TYPE(TAGS)        ,INTENT(INOUT) :: YDAGS
TYPE(TAGF)        ,INTENT(INOUT) :: YDAGF
TYPE(TMLM)        ,INTENT(INOUT) :: YDMLM
TYPE(TFLAKE)      ,INTENT(INOUT) :: YDFLAKE
TYPE(TOCEAN_ML)   ,INTENT(INOUT) :: YDOCEAN_ML
TYPE(TURB)        ,INTENT(INOUT) :: YDURB
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRALFMINPSN
REAL(KIND=JPRB)   ,INTENT(IN), OPTIONAL :: PRCIMIN

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUSURF_CTL_MOD:SUSURF_CTL',0,ZHOOK_HANDLE)
ASSOCIATE(NCSS=>YDDIM%NCSS, NMONTH=>YDDIM%NMONTH, NTILES=>YDDIM%NTILES,NCSNEC=>YDDIM%NCSNEC)

! Initialise dimensions

NCSS=KCSS
NTILES=KTILES
NMONTH=12
NCSNEC=KCSNEC

! Initialise fundamental constants

CALL SUSCST(PRCORIOI,PRPLRG,YDCST)

! Initialise constants defined from fundamental constants

CALL SUSTHF(YDCST)

! Initialise radiation constants

CALL SUSRAD(KSW,KSIL,KTSW,KLWEMISS,&
          & LD_LLCCNL,LD_LLCCNO,KALBEDOSCHEME,KEMISSSCHEME,&
          & PTSTAND,PXP,PRCCNSEA,PRCCNLND,PRSUN,&
          & YDDIM,YDRAD,YDRDI,YDLW,YDSW)

! Initialise soil constants

IF(PRESENT(PRCIMIN)) THEN
  CALL SUSSOIL(PTHRFRTI,LD_LEVGEN,LD_LESSRO,LD_LESN09,LD_LESNML,PNSNMLWS,&
             & YDDIM,YDCST,YDSOIL,PRALFMINPSN,PRCIMIN=PRCIMIN)
ELSE
  CALL SUSSOIL(PTHRFRTI,LD_LEVGEN,LD_LESSRO,LD_LESN09,LD_LESNML,PNSNMLWS,&
             & YDDIM,YDCST,YDSOIL,PRALFMINPSN)
ENDIF

! Initialise vegetation constants

CALL SUSVEG(LD_LELAIV,LD_LECTESSEL,LD_LEAGS,LD_LEFARQUHAR,LD_LEAIRCO2COUP,PRLAIINT,&
          & YDDIM,YDCST,YDSOIL,YDVEG)

! Initialise CTESSEL constants
IF (LD_LECTESSEL) then
  CALL SUCOTWO(PRVR0VT,YDVEG,YDCST,YDAGS)
  IF (LD_LEFARQUHAR) THEN
     CALL SUFARQUHAR(PRVCMAX25,PRHUMREL,PRA1,PRB1,PRG0,PRGM25,PRE_VCMAX,&
            &  PRE_JMAX,YDVEG,YDAGS,YDAGF)
  ENDIF
ENDIF
! Initialise surface exchange coefficient constants

CALL SUVEXC(LD_LEOCWA,LD_LEOCCO,LD_LEOCSA,LD_LEOCLA,&
          & LD_LWCOU, LD_LWCOU2W, LD_LWCOUHMF,&
          & LD_LSCMEC,LD_LROUGH,PEXTZ0M,PEXTZ0H,PRPLRG,&
          & YDEXC)

! Initialise static stability functions constant

CALL SUVEXCS

CALL SUSFLAKE(LD_LEFLAKE,YDFLAKE)

! Initialize kpp ocean mixed layer model.




!OIFS_rm CALL SUSOCEAN_ML(LD_LEOCML,YDOCEAN_ML)


! Initialize tke based ocean mixed layer model.

IF (LD_LOCMLTKE) THEN
  CALL SUGRIDMLM(LD_LOCMLTKE,KCOM,YDMLM)
ELSE
  YDMLM%LOCMLTKE = .FALSE.
ENDIF

! Initialize urban constants. (Always initialises even if LEURBAN = FALSE)

CALL SUSURB(YDURB)

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUSURF_CTL_MOD:SUSURF_CTL',1,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

END SUBROUTINE SUSURF_CTL
END MODULE SUSURF_CTL_MOD
