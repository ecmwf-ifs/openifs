! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE AER_PHY2 &
 &( YDGEOMETRY, YDMODEL,KIDIA, KFDIA, KLON , KTDIA , KLEV , KFLDX, KLEVX, KTILES, KSTGLO, &
 &  KTRAC, KAERO, KSW, &
 &  PRS1 , PRSF1, PAPHI, PTP   , PVERVEL, PCEN , PGEOH, &
 &  PALB , PALBD, PALUVD,&
 &  PAERDEP,PAERLTS,PAERSCC    , PAERWS,PAERGUST, PAERUST,&
 &  PURBF, PRSURF, PSSURF, &
 &  PSO2DD, PDSF, PDSZ, PDMS, &
 &  PCI  , PFRTI, PLSM, PCLK, PQP, PSST, PSNS  , &
 &  PTL  ,PGELAM,PGELAT,PGEMU  , PHSDFOR, &
 &  PUP  , PVP  , PWSA1, PTSPHY, PZ0M, KCHEM, &
 &  KTVL, KTVH, PCVL, PCVH, PLAIL, PLAIH, &
!-- OUTPUTS
 &  PCFLX, PTENC, PDDVLC, &
 &  PDMSO, PLDAY, PLISS, PSO2 , PTDMS, PAERDDP, PAERSDM, PAERSRC, PAERMAP, PAERFLX, PAERLIF, &
 &  PDMSI, PODMS, PEXTRA )

!**** *AER_PHY2* - ROUTINE DEALING WITH AEROSOL SOURCES, DRY DEPOSITION 
!                  AND SEDIMENTATION


!      JJ.MORCRETTE , ECMWF


!**   INTERFACE.
!     ----------
!          *AER_PHY2* IS CALLED FROM *CALLPAR*.

! INPUTS:
! -------
! PRS1 (KLON,0:KLEV)  : HALF-LEVEL PRESSURE           (Pa)
! PRSF1(KLON,KLEV)    : FULL-LEVEL PRESSURE           (Pa)
! PAPHI(KLON,0:KLEV)  : GEOPOTENTIAL ON HALF-LEVELS
! PTP(KLON,KLEV)      : FULL-LEVEL TEMPERATURE (W. DYN.TEND.) (K)
! PVERVEL(KLON,KLEV)  : FULL-LEVEL VERTICAL VELOCITY  (Pa s-1)
! PCEN(KLON,KLEV,KTRAC): TRACER CONCENTRATION
! PGEOH(KLON,0:KLEV)  : GEOPOTENTIAL ON HALF-LEVELS
! PALB(KLON)          : SURFACE FORECAST ALBEDO
! PALUVD(KLON)        : MODIS ALBEDO UV-VIS. DIFFUSE
! PCI(KLON)           : FRACTION OF SEA-ICE
!-- field for biomass burning emission heights
! PFRTI(KLON,KTILES)  : FRACTION OF VARIOUS SURFACES (TILES)
! PLSM(KLON)          : LAND-SEA MASK (note: treats lakes as sea)
! PCLK(KLON)          : LAKE COVER
! PQP(KLON,KLEV)      : FULL-LEVEL HUMIDITY (W. DYN.TEND.) (kg kg-1)
! PSST(KLON)          : SST
! PSNS(KLON)          : SNOW DEPTH
! PTL(KLON)           : SURFACE TEMPERATURE
! PGELAM(KLON)        : LONGITUDE
! PGELAT(KLON)        : LATITUDE (RADIANS)
! PGEMU(KLON)         : SINE OF LATITUDE
! PUP(KLON)           : LOWEST MODEL U-COMPONENT OF WIND (W. DYN.TEND.) 
! PVP(KLON)           : LOWEST MODEL V-COMPONENT OF WIND (W. DYN.TEND.) 
! PWSA1(KLON)         : MOISTURE IN TOP LAYER OF SURFACE
! PTSPHY              : PHYSICS TIME-STEP
! PZ0M                : roughness length for momentum              (m)
!-- aerosol climatological fields
! PAERDEP, PAERLTS, PAERSCC, PDSF, PDMS
!-- dust aerosol predictors
! PAERWS              : Wind speed (average of horizontal wind speed)   (m/s)
! PAERGUST            : Wind gust (maximum 3 second gust in the hour)   (m/s)
! PAERUST             : Friction velocity 
! PHSDFOR             : Standard deviation of filtered sub-grid orography 
! PRSURF, PSSURF      : Surface values for rain and snow

! INPUTS/OUTPUTS:
! ---------------
! PTENC(KLON,KLEV,KTRAC)    : TENDENCY OF TRACER CONCENTRATION
! PDDVLCC(KLON,KLEV,KTRAC)  : DRY DEPOSITION VELOCITY
 

! OUTPUTS:
! --------
! PCFLX(KLON,KTRAC)      : SURFACE FLUX OF TRACERS              (xx m-2)
! PAERDDP(KLON,NACTAERO) : DRY DEPOSITION FLUX                  (xx m-2)
! PAERSDM(KLON,NACTAERO) : SEDIMENTATION FLUX                   (xx m-2)
! PAERSRC(KLON,NACTAERO) : SOURCE FLUX                          (xx m-2) 
! PAERMAP(KLON,5)        : DUST MASK-RELATED QUANTITIES
! PAERFLX(KLON,KBINDD)   : DIAGNOSTIC DUST SOURCE FLUXES
! PAERLIF(KLON,KBINDD)   : DIAGNOSTIC LIFTING THRESHOLD SPEED

!     EXTERNALS.
!     ----------
!          *AER_SRC*, *AER_DRYDEP*, *AER_SEDIMNT*

!     AUTHOR.
!     -------
!          Original JJ.Morcrette, 20060220
!

!     SWITCHES.
!     --------

!     MODEL PARAMETERS
!     ----------------

!     Modifications:
!     --------------
!      K. Yessad (July 2014): Move some variables.
!-----------------------------------------------------------------------

USE TYPE_MODEL   , ONLY : MODEL
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST   , ONLY : RD, RG, RPI
USE YOMCT3   , ONLY : NSTEP
USE YOMLUN   , ONLY : NULOUT
USE YOMCHEM  , ONLY : IEXTR_EM

!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY),    INTENT(IN)    :: YDGEOMETRY
TYPE(MODEL),       INTENT(INOUT) :: YDMODEL
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA, KFDIA, KLON, KFLDX, KLEVX
INTEGER(KIND=JPIM),INTENT(IN)    :: KTDIA, KLEV ,KSTGLO
INTEGER(KIND=JPIM),INTENT(IN)    :: KSW
INTEGER(KIND=JPIM),INTENT(IN)    :: KTILES 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTRAC
INTEGER(KIND=JPIM),INTENT(IN)    :: KAERO(YDMODEL%YRML_GCONF%YGFL%NAERO)
REAL(KIND=JPRB),   INTENT(IN)    :: PRS1(KLON,0:KLEV), PRSF1(KLON,KLEV), PAPHI(KLON,0:KLEV)
REAL(KIND=JPRB),   INTENT(IN)    :: PTP(KLON,KLEV)   , PVERVEL(KLON,KLEV)
REAL(KIND=JPRB),   INTENT(IN)    :: PGEOH(KLON,0:KLEV)
!!!!!!!! PCEN is IN or INOUT
REAL(KIND=JPRB),   INTENT(INOUT) :: PCEN(KLON,KLEV,KTRAC)

REAL(KIND=JPRB),   INTENT(IN)    :: PALB(KLON) , PALBD(KLON,KSW)
REAL(KIND=JPRB),   INTENT(IN)    :: PALUVD(KLON)
REAL(KIND=JPRB),   INTENT(IN)    :: PAERDEP(KLON), PAERLTS(KLON), PAERSCC(KLON)
REAL(KIND=JPRB),   INTENT(IN)    :: PAERWS(KLON) , PAERGUST(KLON), PAERUST(KLON)
REAL(KIND=JPRB),   INTENT(IN)    :: PURBF(KLON)
REAL(KIND=JPRB),   INTENT(IN)    :: PRSURF(KLON), PSSURF(KLON)
REAL(KIND=JPRB),   INTENT(IN)    :: PSO2DD(KLON)
REAL(KIND=JPRB),   INTENT(IN)    :: PDSF(KLON), PDSZ(KLON)
REAL(KIND=JPRB),   INTENT(IN)    :: PDMS(KLON)  
REAL(KIND=JPRB),   INTENT(IN)    :: PCI(KLON)  , PFRTI(KLON,KTILES), PLSM(KLON), PCLK(KLON), PSNS(KLON)
REAL(KIND=JPRB),   INTENT(IN)    :: PSST(KLON)
REAL(KIND=JPRB),   INTENT(IN)    :: PGELAM(KLON), PGELAT(KLON), PGEMU(KLON), PHSDFOR(KLON) 
REAL(KIND=JPRB),   INTENT(IN)    :: PTL(KLON)  , PQP(KLON,KLEV), PUP(KLON)  , PVP(KLON)  , PWSA1(KLON)
REAL(KIND=JPRB),   INTENT(IN)    :: PTSPHY
REAL(KIND=JPRB),   INTENT(IN)    :: PZ0M(KLON)
INTEGER(KIND=JPIM),INTENT(IN)    :: KCHEM(YDMODEL%YRML_GCONF%YGFL%NCHEM)
REAL(KIND=JPRB),   INTENT(IN)    :: PCVL(KLON) ! LOW VEGETATION COVER 
REAL(KIND=JPRB),   INTENT(IN)    :: PCVH(KLON) !  HIGH VEGETATIONCOVER 
REAL(KIND=JPRB),   INTENT(IN)    :: PLAIL(KLON) ! LOW VEGETATION LAI  
REAL(KIND=JPRB),   INTENT(IN)    :: PLAIH(KLON) ! HIGH VEGETATION LAI  
INTEGER(KIND=JPIM),INTENT(IN)    :: KTVL(KLON)  ! LOW VEGETATION TYPE   
INTEGER(KIND=JPIM),INTENT(IN)    :: KTVH(KLON)  ! HIGH VEGETATION TYPE 

REAL(KIND=JPRB),   INTENT(INOUT) :: PTENC(KLON,KLEV,KTRAC)
REAL(KIND=JPRB),   INTENT(OUT)   :: PCFLX(KLON,KTRAC)
REAL(KIND=JPRB),   INTENT(INOUT) :: PDDVLC(KLON,KTRAC)
REAL(KIND=JPRB),   INTENT(OUT)   :: PAERDDP(KLON,YDMODEL%YRML_GCONF%YGFL%NACTAERO)
REAL(KIND=JPRB),   INTENT(OUT)   :: PAERSDM(KLON,YDMODEL%YRML_GCONF%YGFL%NACTAERO)
REAL(KIND=JPRB),   INTENT(OUT)   :: PAERSRC(KLON,YDMODEL%YRML_GCONF%YGFL%NACTAERO)
REAL(KIND=JPRB),   INTENT(OUT)   :: PAERMAP(KLON,5)
REAL(KIND=JPRB),   INTENT(OUT)   :: PAERFLX(KLON,YDMODEL%YRML_PHY_RAD%YREAERATM%NTYPAER(2))
REAL(KIND=JPRB),   INTENT(OUT)   :: PAERLIF(KLON,YDMODEL%YRML_PHY_RAD%YREAERATM%NTYPAER(2))
REAL(KIND=JPRB),   INTENT(OUT)   :: PDMSO(KLON), PLDAY(KLON), PLISS(KLON), PSO2(KLON), PTDMS(KLON)
REAL(KIND=JPRB),   INTENT(OUT)   :: PDMSI(KLON), PODMS(KLON)
REAL(KIND=JPRB),   INTENT(INOUT) :: PEXTRA(KLON,KLEVX,KFLDX)

!-----------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JAER, JK, JL, JRES, ITYP, IBIN
INTEGER(KIND=JPIM) :: JTRESUS=2
INTEGER(KIND=JPIM) :: IINDRESUS(2)

INTEGER(KIND=JPIM) :: INBSU, INBVASH,INBVSO2, INBVSO4, JTAB, IRH(KLON,KLEV)

LOGICAL :: LLPRINT
REAL(KIND=JPRB) :: ZAER(KLON,KLEV) , ZAEROP(KLON,KLEV,YDMODEL%YRML_GCONF%YGFL%NACTAERO)
REAL(KIND=JPRB) :: ZTAERI(KLON,KLEV,YDMODEL%YRML_GCONF%YGFL%NACTAERO), ZTAERO(KLON,KLEV,YDMODEL%YRML_GCONF%YGFL%NACTAERO)
REAL(KIND=JPRB) :: ZTAERA(KLON,KLEV), ZTAERZ(KLON,KLEV)
REAL(KIND=JPRB) :: ZALT(KLON,0:KLEV), ZDP(KLON,KLEV), ZDZ(KLON,KLEV)
REAL(KIND=JPRB) :: ZFAER(KLON), ZFAERO(KLON,YDMODEL%YRML_GCONF%YGFL%NACTAERO), ZFDRYD(KLON,YDMODEL%YRML_GCONF%YGFL%NACTAERO)
REAL(KIND=JPRB) :: ZVDEP(KLON,YDMODEL%YRML_GCONF%YGFL%NACTAERO), ZVEG(KLON),ZVEG1(KLON),ZVEG2(KLON)
REAL(KIND=JPRB) :: ZRHO(KLON,KLEV), ZUSTRBARE, ZRATE_RESUS
REAL(KIND=JPRB) :: ZQSAT(KLON,KLEV)
REAL(KIND=JPRB) :: ZRHCL3D(KLON,KLEV), ZQRWP  , ZQSWP, ZAERORES, ZURBF
REAL(KIND=JPRB) :: ZRWPWP , ZRWSAT,ZSWETN, ZSEDIMV(KLON,KLEV)
REAL(KIND=JPRB) :: ZNG(KLON,KLEV),ZRHOP(KLON,KLEV),ZMFPA(KLON,KLEV),ZDVISC(KLON,KLEV)
REAL(KIND=JPRB) :: ZLNSQSG, ZWETD(KLON,KLEV,3)
REAL(KIND=JPRB) :: ZPREF(KLON,KLEV), ZVBA (KLON,KLEV)
REAL(KIND=JPRB) :: ZTERM1
REAL(KIND=JPRB) :: ZTERM2
REAL(KIND=JPRB) :: ZTERM3
REAL(KIND=JPRB) :: ZTERM4

REAL(KIND=JPRB), PARAMETER :: ZRVONKAR = 0.4_JPRB
REAL(KIND=JPRB), PARAMETER, DIMENSION(3)    :: ZD_SS=(/0.4E-6_JPRB,3.5E-6_JPRB,10.E-6_JPRB/)  ! SS
REAL(KIND=JPRB), PARAMETER :: ZMA=4.78E-26_JPRB
REAL(KIND=JPRB), PARAMETER :: ZBOLTZ = 1.38E-23_JPRB

REAL(KIND=JPRB) :: ZWND(KLON)
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------

#include "aer_src.intfb.h"
#include "aer_drydep.intfb.h"
#include "aer_sedimnt.intfb.h"
#include "chem_inext.intfb.h"
#include "surf_inq.h"
#include "aer_resuspension.intfb.h"
#include "satur.intfb.h"

!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('AER_PHY2',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM,YDMP=>YDGEOMETRY%YRMP, &
 & YGFL=>YDMODEL%YRML_GCONF%YGFL, YDML_GCONF=>YDMODEL%YRML_GCONF, YDML_CHEM=>YDMODEL%YRML_CHEM, &
 & YDCOMPO=>YDMODEL%YRML_CHEM%YRCOMPO,YDEAERATM=>YDMODEL%YRML_PHY_RAD%YREAERATM, &
 & YDEAERSNK=>YDMODEL%YRML_PHY_AER%YREAERSNK, &
 & YDEPHY=>YDMODEL%YRML_PHY_EC%YREPHY, &
 & YDEAERSRC=>YDMODEL%YRML_PHY_AER%YREAERSRC,YDRIP=>YDMODEL%YRML_GCONF%YRRIP)
ASSOCIATE(NACTAERO=>YGFL%NACTAERO, NAERO=>YGFL%NAERO, &
 & LAERDRYDP=>YDEAERATM%LAERDRYDP, LAERSEDIM=>YDEAERATM%LAERSEDIM, &
 & LAERSEDIMSS=>YDEAERATM%LAERSEDIMSS, &
 & LAERSURF=>YDEAERATM%LAERSURF, LAERVOL=>YDEAERATM%LAERVOL, &
 & NDRYDEP=>YDEAERSNK%NDRYDEP, &
 & RSSDENS_RHTAB=>YDEAERSNK%RSSDENS_RHTAB,                 &
 & RSSGROWTH_RHTAB=>YDEAERSNK%RSSGROWTH_RHTAB,             &
 & RRHTAB=>YDEAERSNK%RRHTAB,                               &
 & YR=>YGFL%YR, YS=>YGFL%YS,&
 & NTYPAER=>YDEAERATM%NTYPAER,&
 & LAERRESUSPENSION => YDCOMPO%LAERRESUSPENSION, &
  & YSURF=>YDEPHY%YSURF, &
 & LCOMPO_DDFLX_DIR=>YDCOMPO%LCOMPO_DDFLX_DIR, &
 & NSTART=>YDRIP%NSTART, LCHEM_DIA=>YDCOMPO%LCHEM_DIA, NCHEM=>YGFL%NCHEM, &
 & YAERO_DESC=>YDEAERATM%YAERO_DESC)
!     ------------------------------------------------------------------

!*         0.     INFORMATION ON (LMD-Z-TYPE) MASS-BASED AEROSOL SCHEME
!                 -----------------------------------------------------

! Useful references for various parts of the aerosol model:

! Boucher, O., M. Pham, and C. Venkataraman, 2002: Note Sci., IPSL 23, 27 pp
!  http://www.ipsl.jussieu.fr/poles/Modelisation/NotesSciences.htm

! Giorgi, F., and W.L. Chameides, 1986: J. Geophys. Res., 91, 14367-14376.

! Reddy, M.S., et al., 2005: J. Geophys. Res., 110, D10S16, doi:10.1029/2004JD004757.

! Huneeus, N., et al., 2009: Geosci. Model Dev., 2, 213-229.

! Morcrette, J.-J., et al., 2009: J. Geophys. Res., 114, D06206, doi:10.1029/2008JD011235.


!*         1.     PROGNOSTIC AEROSOLS - INITIAL COMPUTATIONS
!                 ------------------------------------------


ZAEROP(:,:,:) = 0.0_JPRB

IINDRESUS(1:JTRESUS)=(3,6)

LLPRINT=.FALSE.

ZAEROP(KIDIA:KFDIA,1:KLEV,1:NACTAERO) = PCEN(KIDIA:KFDIA,1:KLEV,KAERO(1):KAERO(NACTAERO))
ZTAERI(KIDIA:KFDIA,1:KLEV,1:NACTAERO) = PTENC(KIDIA:KFDIA,1:KLEV,KAERO(1):KAERO(NACTAERO))

CALL SATUR (KIDIA, KFDIA, KLON  , KTDIA, KLEV, YDMODEL%YRML_PHY_SLIN%YREPHLI%LPHYLIN, PRSF1, PTP, ZQSAT, 2)
ZWETD(:,:,1)=ZD_SS(1)
ZWETD(:,:,2)=ZD_SS(2)
ZWETD(:,:,3)=ZD_SS(3)
DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    ZRHO(JL,JK)=PRSF1(JL,JK)/(RD*PTP(JL,JK))
    ZDP(JL,JK)= PRS1(JL,JK) - PRS1(JL,JK-1)
    ZDZ(JL,JK)= ZDP(JL,JK) / (ZRHO(JL,JK)*RG)  
    ZRHCL3D(JL,JK)=PQP(JL,JK)/ZQSAT(JL,JK)
    IRH(JL,JK)=1
    DO JTAB=1,12
      IF (ZRHCL3D(JL,JK)*100._JPRB > RRHTAB(JTAB)) THEN
        IRH(JL,JK)=JTAB
      ENDIF
    ENDDO
    ZRHOP(JL,JK)=RSSDENS_RHTAB(IRH(JL,JK))
    ZWETD(JL,JK,1)=ZD_SS(1)*RSSGROWTH_RHTAB(IRH(JL,JK))
    ZWETD(JL,JK,2)=ZD_SS(2)*RSSGROWTH_RHTAB(IRH(JL,JK))
    ZWETD(JL,JK,3)=ZD_SS(3)*RSSGROWTH_RHTAB(IRH(JL,JK))
  ENDDO
ENDDO
DO JL=KIDIA,KFDIA
  ZALT(JL,KLEV)=PAPHI(JL,KLEV)/RG
ENDDO
DO JK=KLEV-1,0,-1
  DO JL=KIDIA,KFDIA
    ZALT(JL,JK)=ZALT(JL,JK+1)+ZDZ(JL,JK+1)
  ENDDO
ENDDO


! NB: The physical calculations for the prognostic aerosols are done in a 
!     sequential mode, with any given process adding its tendency to the 
!     initial tendency (giving ZTENC), the field itself (ZCEN) remaining 
!     the same.
 
! NB: for multi-bins, 1-3 is sea-salt, 4-6 is desert dust

! define a composite index for each bin of each different aerosol type to be used
! in source, sedimentation and deposition routines

! NB: now done in SU_AERW 


!-- set all fluxes relevant to surface emissions, dry deposition and sedimentation to zero
IF (.NOT.LAERSURF) THEN
  PAERMAP(:,:)  =0._JPRB
  PAERFLX(:,:)  =0._JPRB
  PAERLIF(:,:)  =0._JPRB
  PAERDDP(:,:)  =0._JPRB
  PAERSDM(:,:)  =0._JPRB
  PAERSRC(:,:)  =0._JPRB
  ZWND(:)       =0._JPRB
  DO JAER=1,NACTAERO
    DO JL=KIDIA,KFDIA
      ZFAERO(JL,JAER) =0._JPRB
      PCFLX(JL,KAERO(JAER))=0._JPRB
    ENDDO
  ENDDO 
ELSEIF (LAERSURF) THEN

!*         2.     SURFACE FLUXES OF AEROSOLS
!                 --------------------------

! N.B: COMPUTATIONS USE SURFACE PRECIPITATION AT T-DT

  DO JL=KIDIA,KFDIA
!    ZWND(JL)=SQRT(PUP(JL)*PUP(JL)+PVP(JL)*PVP(JL))
    ZWND(JL)=PAERWS(JL)
    PDMSO(JL)=0._JPRB
    PLDAY(JL)=0._JPRB
    PLISS(JL)=0._JPRB
    PSO2(JL) =0._JPRB
    PTDMS(JL)=0._JPRB
    PODMS(JL)=0._JPRB
    PDMSI(JL)=PDMS(JL)
  ENDDO  

!- surface (mainly) fluxes for different tropospheric aerosols


  CALL AER_SRC &
    &( YDGEOMETRY, YDMODEL, KIDIA  , KFDIA , KLON , KTDIA, KLEV , KTILES, NSTART, NSTEP , KSTGLO, &
    &  KSW    , KTRAC , KAERO, & 
    &  PALB   , PALBD , PALUVD, PAPHI , &
    &  PAERDEP,PAERLTS, PAERSCC, PAERGUST, ZALT, &
    &  PDSF   , PDSZ  , PDMS, &
    &  PRS1   , PRSF1 , PCI  , PGELAM, PGELAT, PGEMU, PFRTI, PHSDFOR, &
    &  PLSM   , PCLK  , PQP   , ZRHO , PSNS , PTP  , PTL  , PTSPHY, PZ0M, KCHEM, &
    &  ZWND   , PWSA1 , PSST, &
    &  PDMSO  , PLDAY , PLISS, PSO2 , PTDMS, PAERFLX, PAERLIF, PAERMAP,PCFLX, PCEN , PTENC, &
    &  PODMS)

!-- VERY IMPORTANT
! sea salt and desert dust fluxes in kg m-2 s-1, thus:
! dynamics, vert.diff, convection see mixing ratio in kg kg-1 (PCEN, ZAEROP) 
! and corresponding tendencies ZTENC are in kg kg-1 s-1 

  DO JAER=1,NACTAERO
    DO JL=KIDIA,KFDIA
      PAERSRC(JL,JAER)=PCFLX(JL,KAERO(JAER))
    ENDDO
  ENDDO
! collect distrubuted fluxes  
 IF (LCHEM_DIA) THEN
    CALL CHEM_INEXT( KIDIA , KFDIA  , KLON , KLEV , NACTAERO, NACTAERO ,  &
   &    ZDP, PTSPHY, PTENC(:,:,KAERO(1):KAERO(NACTAERO)) , ZTAERI, PEXTRA(:,NCHEM+1:NCHEM+NACTAERO,IEXTR_EM))
 ENDIF

!-- update tendencies if volcanic activity
IF (LAERVOL) THEN
  INBSU=NTYPAER(1)+NTYPAER(2)+NTYPAER(3)+NTYPAER(4)+NTYPAER(5)
  INBVASH=INBSU+1
  ZTAERI(KIDIA:KFDIA,1:KLEV,INBVASH) = PTENC(KIDIA:KFDIA,1:KLEV,KAERO(INBVASH))
  INBVSO4=INBSU+2
  ZTAERI(KIDIA:KFDIA,1:KLEV,INBVSO4) = PTENC(KIDIA:KFDIA,1:KLEV,KAERO(INBVSO4))
  INBVSO2=INBSU+3
  ZTAERI(KIDIA:KFDIA,1:KLEV,INBVSO2) = PTENC(KIDIA:KFDIA,1:KLEV,KAERO(INBVSO2))
ENDIF

!*         3.     DRY DEPOSITION INCLUDED AS MODIFICATION TO SURFACE FLUXES
!                 ---------------------------------------------------------

  IF (LAERDRYDP) THEN

!-- Dry deposition is applied to all aerosols: Reddy et al., 2005, $2.6, p.3 and Table 1

    PAERDDP(:,:) = 0.0_JPRB

    CALL AER_DRYDEP ( &
       & YDGEOMETRY, YDMODEL, &
       & KIDIA , KFDIA, KLON , KTDIA, KLEV , NSTEP , KSTGLO, KTILES,&
       & ZAEROP, PCI  , PCFLX(:,KAERO(1):KAERO(NACTAERO)), &
       & PRS1  , PRSF1, ZDP  , PGEOH, PLSM , ZRHO , ZTAERI, PTSPHY,&
       & PSO2DD, PGELAM, PGELAT, PAERUST, PZ0M, &
       & PTP, PQP(KIDIA:KFDIA,KLEV),   &
       & ZDZ(KIDIA:KFDIA,KLEV), &
       & PCVL, PCVH, KTVL, KTVH, PLAIL, PLAIH, PFRTI, PSNS, &
       & ZFAERO, ZTAERO, ZFDRYD, ZVDEP, ZVEG, ZVEG1, ZVEG2)

    IF (LCOMPO_DDFLX_DIR) THEN
      ! The dry deposition is passed through to the vertical diffusion as
      ! deposition velocity, with the surface flux untouched.  Diagnostics
      ! must still be output using either the deposition flux or rate
      ! calculated in AER_DRYDEP as per NDRYDEP, however.
      DO JAER=1,NACTAERO
        DO JL=KIDIA,KFDIA
          PDDVLC(JL,KAERO(JAER)) = ZVDEP(JL,JAER)
        ENDDO
      ENDDO

      IF (NDRYDEP == 1) THEN
        PAERDDP(KIDIA:KFDIA,1:NACTAERO)=ZFAERO(KIDIA:KFDIA,1:NACTAERO)-PCFLX(KIDIA:KFDIA,KAERO(1:NACTAERO))
      ELSE
        PAERDDP(KIDIA:KFDIA,1:NACTAERO)=ZFDRYD(KIDIA:KFDIA,1:NACTAERO)
      ENDIF
    ELSEIF (NDRYDEP == 1) THEN
!--   The dry deposition at the surface is the difference between the original 
!     emission flux and the flux from which the effect of the dry depoasition 
!     has been substracted. This latter flux is the one used afterwards in the 
!     vertical diffusion.

      DO JAER=1,NACTAERO
        DO JL=KIDIA,KFDIA
          PAERDDP(JL,JAER)=ZFAERO(JL,JAER)-PCFLX(JL,KAERO(JAER))
          PCFLX(JL,KAERO(JAER))=ZFAERO(JL,JAER)
        ENDDO
      ENDDO

    ELSE
!--   The dry deposition at the surface has been explicitly computed, and the 
!     tendency in the first layer above the surface is to be updated. The 
!     emission flux is kept untouched.

      DO JAER=1,NACTAERO
        DO JL=KIDIA,KFDIA
          PAERDDP(JL,JAER)=ZFDRYD(JL,JAER)
        ENDDO
      ENDDO
!-- update tendencies
      PTENC(KIDIA:KFDIA,1:KLEV,KAERO(1):KAERO(NACTAERO)) = ZTAERO(KIDIA:KFDIA,1:KLEV,1:NACTAERO)
    ENDIF

  ELSE
!-- If dry depo. is not called at all, the emission fluxes remain untouched, 
!    as do the lowest layer tendencies/mixing ratios 
    DO JAER=1,NACTAERO
      DO JL=KIDIA,KFDIA
        ZFAERO(JL,JAER) =0._JPRB
        PAERDDP(JL,JAER)=0._JPRB
      ENDDO
    ENDDO    
  ENDIF

ENDIF

!*         4.     SEDIMENTATION OF AEROSOLS
!                 -------------------------

IF (LAERSEDIM) THEN

  PAERSDM(:,:) = 0.0_JPRB

  DO JAER=1,NACTAERO
    ! Skip species if no sedimentation velocity configured
    ITYP=YAERO_DESC(JAER)%NTYP
    IBIN=YAERO_DESC(JAER)%NBIN
    IF (YAERO_DESC(JAER)%RSEDIMV == 0._JPRB .AND. (ITYP /=1)) CYCLE
    IF (YAERO_DESC(JAER)%RSEDIMV == 0._JPRB .AND. ITYP ==1 .AND. IBIN ==1) CYCLE

    DO JK=1,KLEV
      DO JL=KIDIA,KFDIA
        ZAER(JL,JK)  =PCEN(JL,JK,KAERO(JAER))
        ZTAERA(JL,JK)=PTENC(JL,JK,KAERO(JAER))
      ENDDO
    ENDDO
    ZSEDIMV(KIDIA:KFDIA,1:KLEV)=YAERO_DESC(JAER)%RSEDIMV
    ZDVISC(:,:)=0._JPRB
    IF (ITYP == 1 .AND. (IBIN==2 .OR. IBIN==3) .AND. LAERSEDIMSS) THEN
      ZLNSQSG=LOG(1.8_JPRB)*LOG(1.8_JPRB)
      ZTERM1=8.0_JPRB
      ZTERM2=3.5_JPRB
      ZTERM3=EXP(ZTERM1*ZLNSQSG)
      ZTERM4=1.246_JPRB*EXP(ZTERM2*ZLNSQSG)
      ZDVISC(KIDIA:KFDIA,1:KLEV)  = 1.83E-5_JPRB*(416.16_JPRB/(PTP(KIDIA:KFDIA,1:KLEV)+120.0_JPRB))* &
             & (SQRT(PTP(KIDIA:KFDIA,1:KLEV)/296.16_JPRB)**3.0_JPRB)
      ZVBA(KIDIA:KFDIA,1:KLEV)    = SQRT(8.0_JPRB*ZBOLTZ*PTP(KIDIA:KFDIA,1:KLEV)/(RPI*ZMA)) ! mean free speed of air molecules
      ZMFPA(KIDIA:KFDIA,1:KLEV)   = 2.0_JPRB*ZDVISC(KIDIA:KFDIA,1:KLEV)/ &
             & (ZRHO(KIDIA:KFDIA,1:KLEV)*ZVBA(KIDIA:KFDIA,1:KLEV)) ! mean free path of air molecules

      ZNG(KIDIA:KFDIA,1:KLEV)=2.0_JPRB*ZMFPA(KIDIA:KFDIA,1:KLEV)/ZWETD(KIDIA:KFDIA,1:KLEV,JAER)
      WHERE (ZDVISC(:,:) > 1.E-12_JPRB)
        ZPREF(:,:)=ZRHOP(:,:)*ZWETD(:,:,JAER)*ZWETD(:,:,JAER)*RG/(18.0_JPRB*ZDVISC(:,:))
        ZSEDIMV(:,:)=0.5_JPRB*ZPREF(:,:)*(ZTERM3+ZTERM4*ZNG(:,:))  
      ENDWHERE
    ENDIF

!-- sedimentation of aerosols follows closely sedimentation of ice particles by A. Tompkins

    CALL AER_SEDIMNT &
      &( KIDIA , KFDIA , KLON , KLEV, ZSEDIMV, &
      &  PCI   , PLSM  , &
      &  PTSPHY, PTP   , PRSF1, PRS1, PVERVEL, & 
      &  ZAER  , ZTAERA, &
      &  ZTAERZ, ZFAER ) 

    IF (LLPRINT .AND. NSTEP == NSTART) THEN
      DO JL=KIDIA,KFDIA
        WRITE(NULOUT,FMT='(" AER_PHY2 ",I4,3E12.5)') JAER,ZAER(JL,JK),ZTAERA(JL,JK),ZTAERZ(JL,JK)
      ENDDO
    ENDIF

    DO JK=1,KLEV
      DO JL=KIDIA,KFDIA
        PTENC(JL,JK,KAERO(JAER))=ZTAERZ(JL,JK)
      ENDDO
    ENDDO
    DO JL=KIDIA,KFDIA
      PAERSDM(JL,JAER)=ZFAER(JL)
    ENDDO

  ENDDO
ELSE
!-- If sedimentation is not called, tendencies after dry depo. stay the same.
  DO JL=KIDIA,KFDIA
    PAERSDM(JL,:)=0._JPRB
  ENDDO
ENDIF


!*            5     RESUSPENSION OF DRY DEPOSITED AEROSOLS
!                 ----------------------------------------------------



IF (LAERRESUSPENSION) THEN
  CALL SURF_INQ(YSURF,PRWPWP=ZRWPWP)
  CALL SURF_INQ(YSURF,PRWSAT=ZRWSAT)


  DO JL=KIDIA,KFDIA
    ZUSTRBARE = PAERGUST(JL)*ZRVONKAR/LOG(10._JPRB/0.0003_JPRB)
    ZSWETN = MIN(1._JPRB, MAX(0.001_JPRB, (PWSA1(JL)-ZRWPWP)/(ZRWSAT-ZRWPWP) ))
    ZAERORES=0._JPRB
    IF (PURBF(JL) >= 0.01_JPRB .AND. PSNS(JL) == 0.0_JPRB) THEN
      DO JRES=1,JTRESUS
        JAER=IINDRESUS(JRES)
        ! Increase urban fraction
        ZURBF=MIN(1.0_JPRB,PURBF(JL)*2.0_JPRB)
        ! add dry deposited and sedimented aerosol
        ZAERORES=(PAERDDP(JL,JAER)+PAERSDM(JL,JAER))*PTSPHY*ZURBF
        ! If saturated soil or rain, then set to 0
        ZQRWP   = MAX(0._JPRB, PRSURF(JL))
        ZQSWP   = MAX(0._JPRB, PSSURF(JL))
        IF (ZSWETN > 0.95_JPRB .OR. ZQRWP > 0._JPRB .OR. ZQSWP > 0._JPRB) THEN
          ZAERORES=0._JPRB
        ENDIF
        IF (ZAERORES > 0._JPRB) THEN
          ! Compute resuspension rate
          CALL AER_RESUSPENSION &
          &( YDEAERSNK, YDEAERSRC, YDEAERATM, JAER , PZ0M(JL), PTP(JL,KLEV), ZRHCL3D(JL,KLEV), ZRHO(JL,KLEV), &
          &  ZUSTRBARE, ZRATE_RESUS )
          ! check that resuspension doesn't go beyond reservoir capacity
          ZRATE_RESUS=MIN(1._JPRB/PTSPHY,ZRATE_RESUS)
          ZRATE_RESUS=MAX(0._JPRB,ZRATE_RESUS)
          ! Update flux and reservoir
          PCFLX(JL,KAERO(JAER))=PCFLX(JL,KAERO(JAER))-ZRATE_RESUS*ZAERORES
         ENDIF
      ENDDO
    ENDIF
  ENDDO
ENDIF



!-----------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('AER_PHY2',1,ZHOOK_HANDLE)
END SUBROUTINE AER_PHY2
