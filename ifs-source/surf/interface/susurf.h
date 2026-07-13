INTERFACE
SUBROUTINE SUSURF(KSW,KCSS,KCSNEC,KSIL,KCOM,KTILES,KTSW,KLWEMISS,&
 & LD_LEFLAKE,LD_LEOCML,LD_LOCMLTKE,&       
 & LD_LWCOU, LD_LWCOU2W, LD_LWCOUHMF,&
 & LD_LLCCNL,LD_LLCCNO,LD_LEVGEN,LD_LESSRO,&
 & LD_LELAIV,LD_LECTESSEL,LD_LEAGS,LD_LEFARQUHAR,&
 & LD_LEAIRCO2COUP,LD_LESN09,LD_LESNML,&
 & LD_LEOCWA,LD_LEOCCO,LD_LEOCSA,LD_LEOCLA,KALBEDOSCHEME,KEMISSSCHEME,&
 & LD_LSCMEC,LD_LROUGH,PEXTZ0M,PEXTZ0H,&
 & PTHRFRTI,PTSTAND,PXP,PRCCNSEA,PRCCNLND,&
 & PRLAIINT,PRSUN,PRCORIOI,PRPLRG,PNSNMLWS,&
 & PRVR0VT,PRVCMAX25,PRHUMREL,PRA1,PRB1,PRG0,PRGM25,PRE_VCMAX,PRE_JMAX,&  
 & YDSURF,PRALFMINPSN,PRCIMIN)

USE PARKIND1, ONLY : JPIM, JPRB
USE, INTRINSIC :: ISO_C_BINDING

! (C) Copyright 2002- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!**   *SUSURF* IS THE SET-UP ROUTINE FOR surface modules contain constants

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
!     P. Viterbo   ECMWF   03-12-2004  Include user-defined RTHRFRTI
!     P. Viterbo   ECMWF   May 2005    Externalise surf
!     G. Balsamo   ECMWF   10-01-2006  Include Vangenucthen Hydro.
!     JJMorcrette 20060511 MODIS albedo
!     V. Stepanen/G. Balsamo May 2008  Lake tile
!     Y. Takaya    ECMWF   07-10-2008  Implement an ocean mixed layer model
!     G. Balsamo   ECMWF   13-10-2008  Include switch for liquid water in snow
!     G. Balsamo   ECMWF   25-02-2009  Include full set of snow switches
!     S. Boussetta/G.Balsamo May 2009  Add LAI monthly fields
!     E. Dutra             16-11-2009  snow 2009 cleaning
!     S. Boussetta/G.Balsamo May 2010  Include CTESSEL switch LECTESSEL
!     S. Boussetta/G.Balsamo June 2011  Include switch LEAGS (for modularity of evap&CO2)
!     P. Bechtold          26-03-2012  Include small planet PRCORIOI PRPLRG
!     R. Hogan             14-01-2019  Replace LE4ALB with KALBEDOSCHEME; add KEMISSSCHEME
!     A. Agusti-Panareda Nov 2020  Include switch LEAIRCO2COUP (to couple atm CO2 with photosynthesis)
!     A. Agusti-Panareda Jul 2021  Include switch LEFARQUHAR for Farquhar photosynthesis model


!  INTERFACE: 

!    Integers (In):

!      KSW       : NUMBER OF SHORTWAVE SPECTRAL INTERVALS
!      KCSS      : Number of soil levels
!      KCSS      : Number of snow levels
!      KSIL      : NUMBER OF (infrared) SPECTRAL INTERVALS
!      KCOM      : Number of layers in mixed layer model 
!      KTILES    : Number of surface tiles
!      KTSW      : Maximum possible number of shortwave spectral intervals
!      KLWEMISS  : Number of longwave emissivity spectral intervals
!      KALBEDOSCHEME : (0) ERBE, (1) 4-component MODIS, (2) 6-component MODIS
!      KEMISSSCHEME : (0) ERBE, (1) 4-component MODIS, (2) 6-component MODIS

!    Logicals (In):

!      LD_LLCCNL : .T. IF CCN CONCENTRATION OVER LAND IS DIAGNOSED
!      LD_LLCCNO : .T. IF CCN CONCENTRATION OVER OCEAN IS DIAGNOSED
!      LD_LEVGEN : .T. IF VAN GENUCHTEN HYDRO IS ACTIVATED
!      LD_LESSRO : .T. IF SUB-GRID SURFACE RUNOFF IS ACTIVATED
!      LD_LELAIV : .T. IF LAI IS READ FROM CLIMATE FIELDS

!      LD_LECTESSEL : .T. IF CTESSEL scheme is used for CO2 fluxes
!      LD_LEAGS  : .T. IF CTESSEL scheme is used for CO2 and Evap fluxes
!      LD_LEFARQUHAR: .T. IF CTESSEL scheme is used for CO2 and Evap fluxes
!      LD_LEAIRCO2COUP  : .T. IF variable atmospheric CO2 is used in photosynthesis
!      LD_LESN09 : .T. IF  SNOW 2009 IS ACTIVATED
!      LD_LEFLAKE: .T. IF FLAKE PARAMETRIZATION USED FOR LAKES 
!      LD_LEOCML : .T. IF OCEAN MIXED LAYER MODEL (KPP) IS ACTIVATED
!      LD_LOCMLTKE : .T. IF OCEAN MIXED LAYER MODEL (TKE) IS ACTIVATED

!      LD_LWCOU    : .T. COUPLING TO WAVE MODEL
!      LD_LWCOU2W  : .T. COUPLING TO WAVE MODEL WITH FEEDBACK TO ATMOSPHERE
!      LD_LWCOUHMF : .T. SEA STATE DEPENDENT HEAT AND MOISTURE FLUXES IF COUPLED TO WAVE MODEL 

!    Reals (In):

!      PTHRFRTI  : ! MINIMUM THRESHOLD FOR TILE FRACTION
!      PTSTAND   : ! REFERENCE TEMPERATURE FOR TEMPERATURE DEPENDENCE
!      PXP       : ! POLYNOMIAL COEFFICIENTS OF PLANCK FUNCTION
!      PRCCNSEA  : ! NUMBER CONCENTRATION (CM-3) OF CCNs OVER SEA
!      PRCCNLND  : ! NUMBER CONCENTRATION (CM-3) OF CCNs OVER LAND
!      PRLAIINT  : ! INTERACTIVE LAI COEFFICIENT (1=Interactive ; 0=climatology)
!      PRSUN     : ! SOLAR FRACTION IN SPECTRAL INTERVALS
!      PRCORIOI  : ! REDUCTION FACTOR FOR LEGNTH OF DAY FOR SMALL PLANET
!      PRPLRG    : ! GRAVITY FACTOR FOR FOR SMALL PLANET
!      PRALFMINPSN:! Albedo of permanent snow
!      PRCIMIN   : ! MINIMUM ICE FRACTION
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
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRVR0VT(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRVCMAX25(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRHUMREL(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRA1(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRB1(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRG0(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRGM25(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRE_VCMAX(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRE_JMAX(:,:)
TYPE(C_PTR)       ,INTENT(OUT)   :: YDSURF 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRALFMINPSN
REAL(KIND=JPRB)   ,INTENT(IN), OPTIONAL :: PRCIMIN
INTEGER(KIND=JPIM),INTENT(IN)    :: PNSNMLWS 




!     ------------------------------------------------------------------

END SUBROUTINE SUSURF
END INTERFACE
