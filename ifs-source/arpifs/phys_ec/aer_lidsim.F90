! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE AER_LIDSIM &
  &( YDEAERATM,YDML_PHY_AER,YDCOMPO,YGFL,KIDIA , KFDIA, KLON , KLEV , KACTAERO, &
  &  PAERO , PAPH , PAP  , PCO2 , PCLDF   , &
  &  PDP   , PDZ  , PNO2 , PO3  , PRHCL   , PT , &
  &  PLISIM, PLISIA, PLISIS, PLISIT &
  & )

!**** *AER_LIDSIM* - SIMULATES THE LIDAR SIGNAL OUT OF THE ATMOSPHERE
!                    ASSUMED CLOUDLESS BUT INCLUDING THE GEMS-AEROSOLS.

!     PURPOSE.
!     --------

!          THIS ROUTINE COMPUTES THE TOTAL BACK-SCATTERING COEFFICIENT 
!          ACCOUNTING RAYLEIGH SCATTERING, AEROSOL SCATTERING AND GASEOUS 
!          ABSORPTION FOR THE THREE USUAL WAVELENGTHS (355, 532 AND 1064 NM)
!          FROM THE STANDARD "LIDAR EQUATION".

!**   INTERFACE.
!     ----------

!          *AER_LIDSIM* IS CALLED FROM *CALLPAR*

!        IMPLICIT ARGUMENTS :
!        --------------------

!     ==== INPUTS ===
! PAERO  (KLON,KLEV,KACTAERO): PROGNOSTIC AEROSOLS
! PAPH   (KLON,0:KLEV) : HALF-LEVEL PRESSURE (Pa)
! PAP    (KLON,KLEV)   : FULL-LEVEL PRESSURE (Pa)
! PCO2   (KLON,KLEV)   : CO2 CONCENTRATION (kg/kg)
! PCLDF  (KLON,KLEV)   : CLOUD FRACTION
! PDP    (KLON,KLEV)   : PRESSURE THICKNESS OF LAYER (Pa)
! PDZ    (KLON,KLEV)   : THICKNESS OF LAYER (m)
! PNO2   (KLON,KLEV)   :
! PO3    (KLON,KLEV)   : OZONE MIXING RATIO (kg/kg) -- prognostic --
! PRHCL  (KLON,KLEV)   : (CLEAR-SKY) RELATIVE HUMIDITY
! PT     (KLON,KLEV)   : TEMPERATURE (K)

!     ==== OUTPUTS ===
! PLISIM (KLON,NWLID,0:KLEV) : UNATTENUATED MOLECULAR BACK-SCATTERING COEFFICIENT
! PLISIA (KLON,NWLID,0:KLEV) : UNATTENUATED AEROSOL BACK-SCATTERING COEFFICIENT
! PLISIS (KLON,NWLID,0:KLEV) : BACK-SCATTERING SIGNAL FOR SURFACE LIDAR
! PLISIT (KLON,NWLID,0:KLEV) : BACK-SCATTERING SIGNAL FOR SATELLITE LIDAR

!     METHOD.
!     -------
!         from Ackerman, J., 1998, J.Atmos.Ocean.Technol., 15, 1043-1050
!         and Huneeus, N. and O. Boucher, 2007, J.Geophys.Res. 

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------

!     AUTHOR.
!     -------
!        JEAN-JACQUES MORCRETTE  *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 20090702
!        JJMorcrette 20130705 bug-fix molecular backscatter within layer
!                    following correction from A.Benedetti in analysis part (op_obs/hop.F90)
!        SRemy  20160916  Add nitrates
!     ------------------------------------------------------------------

USE MODEL_PHYSICS_AEROSOL_MOD , ONLY : MODEL_PHYSICS_AEROSOL_TYPE
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST   , ONLY : RG, RMD, RMCO2, RMNO2, RMO3, RNAVO, RPI
USE YOEAEROP , ONLY : ALF_SU, OMG_SU, ALF_OM, OMG_OM, &
  & ALF_DD, OMG_DD, ALF_SS, ALF_NI, ALF_AM, ALF_SOA, OMG_SS, ALF_BC, OMG_BC, &
  & OMG_NI, OMG_AM, OMG_SOA, RALI_BC,RALI_DD,RALI_OM,RALI_SU,RALI_SS, &
  & RALI_NI, RALI_AM, RALI_SOA
USE YOEAERATM, ONLY : TEAERATM
USE YOMCOMPO , ONLY : TCOMPO
USE YOM_YGFL , ONLY : TYPE_GFLD

IMPLICIT NONE

TYPE(TEAERATM)    ,INTENT(INOUT):: YDEAERATM
TYPE(MODEL_PHYSICS_AEROSOL_TYPE),INTENT(INOUT):: YDML_PHY_AER
TYPE(TCOMPO)      ,INTENT(INOUT):: YDCOMPO
TYPE(TYPE_GFLD)   ,INTENT(INOUT):: YGFL
INTEGER(KIND=JPIM),INTENT(IN) :: KLON 
INTEGER(KIND=JPIM),INTENT(IN) :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN) :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN) :: KFDIA
INTEGER(KIND=JPIM),INTENT(IN) :: KACTAERO

REAL(KIND=JPRB)   ,INTENT(IN) :: PAERO(KLON,KLEV,KACTAERO)
REAL(KIND=JPRB)   ,INTENT(IN) :: PAP(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN) :: PAPH(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(IN) :: PCO2(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN) :: PCLDF(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN) :: PDP(KLON,KLEV), PDZ(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN) :: PNO2(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN) :: PO3(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN) :: PRHCL(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN) :: PT(KLON,KLEV)

REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLISIM(KLON,YDML_PHY_AER%YREAERLID%NWLID,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLISIA(KLON,YDML_PHY_AER%YREAERLID%NWLID,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLISIS(KLON,YDML_PHY_AER%YREAERLID%NWLID,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLISIT(KLON,YDML_PHY_AER%YREAERLID%NWLID,0:KLEV) 
!     ------------------------------------------------------------------

!*       0.1   LOCAL ARRAYS
!              ------------

INTEGER(KIND=JPIM) :: IBIN, IEFRH, IIRH, IK, ITYP, IWAVL, IWL
INTEGER(KIND=JPIM) :: JAER, JL, JK, JTAB, JWL
INTEGER(KIND=JPIM) :: IRH(KLON,KLEV)

LOGICAL :: LLPRINT

REAL(KIND=JPRB) :: ZOTRAY(KLON,YDML_PHY_AER%YREAERLID%NWLID,KLEV), ZOTCO2(KLON,YDML_PHY_AER%YREAERLID%NWLID,KLEV), ZOTO2(KLON, &
 & YDML_PHY_AER%YREAERLID%NWLID,KLEV)
REAL(KIND=JPRB) :: ZOTNO2(KLON,YDML_PHY_AER%YREAERLID%NWLID,KLEV), ZOTO3(KLON,YDML_PHY_AER%YREAERLID%NWLID,KLEV) , ZPATHT(KLON, &
 & YDML_PHY_AER%YREAERLID%NWLID,0:KLEV)
REAL(KIND=JPRB) :: ZPATHS(KLON,YDML_PHY_AER%YREAERLID%NWLID,0:KLEV)

!REAL(KIND=JPRB) :: ZOSRAY(KLON,NWLID,KLEV), ZOSCO2(KLON,NWLID,KLEV), ZOSO2(KLON,NWLID,KLEV)
!REAL(KIND=JPRB) :: ZOSNO2(KLON,NWLID,KLEV), ZOSO3(KLON,NWLID,KLEV) , ZPATHS(KLON,NWLID,0:KLEV)

REAL(KIND=JPRB) :: ZAERMSS(KLON,KLEV,YGFL%NACTAERO), ZALF(YGFL%NACTAERO), ZOMG(YGFL%NACTAERO), ZLIR(YGFL%NACTAERO)
REAL(KIND=JPRB) :: ZAEROD(KLON,YDML_PHY_AER%YREAERLID%NWLID,KLEV)    , ZBSC(KLON,YDML_PHY_AER%YREAERLID%NWLID,KLEV)         , &
 &  ZEXT(KLON,YDML_PHY_AER%YREAERLID%NWLID,KLEV)
REAL(KIND=JPRB) :: ZBSCMOL(KLON,YDML_PHY_AER%YREAERLID%NWLID,KLEV)

REAL(KIND=JPRB) :: ZA1, ZA2, ZA3, ZA4, ZAA, ZAN, ZAN2, ZDP, ZDPG, ZEPSAER, ZFAC, ZFACT, ZLIDRAT, &
  & ZMD, ZNO2MOL, ZNO2TAU, ZNO2V, ZO3MOL, ZO3TAU, ZO3V, ZREF, ZSIGMA, &
  & ZTAUR, ZTE, ZTE2, &
  & ZTRSOT, ZTRTOT, ZUSOT, ZUTOT, ZWLUM, ZWLCM, ZWL2CM, ZWL4CM, ZWN, ZWN2

REAL(KIND=JPRB) :: ZGCO2, ZMCO2, ZPHIC, ZPSIC, ZUCO2, ZUPCO2, ZVCO2
REAL(KIND=JPRB) :: ZO2 , ZGO2 , ZMO2 , ZPHIO, ZPSIO, ZUO2 , ZUPO2, ZBOLTZ

REAL(KIND=JPRB) :: ZCAER(KLON), ZCRAY(KLON), ZCOZO(KLON), ZCOXY(KLON), ZCCO2(KLON), ZCNO2(KLON)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('AER_LIDSIM',0,ZHOOK_HANDLE)
ASSOCIATE(NACTAERO=>YGFL%NACTAERO, &
 & JWLID=>YDML_PHY_AER%YREAERLID%JWLID, NWLID=>YDML_PHY_AER%YREAERLID%NWLID, RLICLS=>YDML_PHY_AER%YREAERLID%RLICLS, &
 & RLICO2=>YDML_PHY_AER%YREAERLID%RLICO2, RLIDELT=>YDML_PHY_AER%YREAERLID%RLIDELT, &
 & RLINO2=>YDML_PHY_AER%YREAERLID%RLINO2, RLINS=>YDML_PHY_AER%YREAERLID%RLINS, RLIO2=>YDML_PHY_AER%YREAERLID%RLIO2, &
 & RLIO3=>YDML_PHY_AER%YREAERLID%RLIO3, RLIPREF=>YDML_PHY_AER%YREAERLID%RLIPREF, RLIT0=>YDML_PHY_AER%YREAERLID%RLIT0, &
 & RLITREF=>YDML_PHY_AER%YREAERLID%RLITREF, RWLID=>YDML_PHY_AER%YREAERLID%RWLID, &
 & RRHTAB=>YDML_PHY_AER%YREAERSNK%RRHTAB, &
 & RSS_RH80_MASSFAC=>YDEAERATM%RSS_RH80_MASSFAC, &
 & YAERO_DESC=>YDEAERATM%YAERO_DESC, &
 & NVOLOPTP=>YDML_PHY_AER%YREAERVOL%NVOLOPTP)

!*         0.     SPECTRAL INFORMATION OF RELEVANCE
!                 ---------------------------------

! Rayleigh |  355nm / 28169cm-1 |  532nm / 18797cm-1 | 1064nm / 9398.5cm-1 |
! Index    |          2         |          8         |         16          ! 
! CO2      |         no         |         no         |        yes          |
! O2       |         no         |         no         |        yes          |
! O3       |        yes         |        yes         |        yes          |
! NO2      |        yes         |        yes         |        yes          |

!*    --------------------------------------------------------------------------

!*         1.     INITIALIZATION AND PREPARATORY WORK
!                 -----------------------------------

!!!ZSIGMA=0.056032_JPRB
ZREF=RLITREF/RLIPREF
ZO2=RMO3*2._JPRB/3._JPRB
ZEPSAER=1.E-16_JPRB
!-- the effective relative hunidity is the low value (20%) assumed for hydrophobic component of OM
IEFRH=3

LLPRINT=.FALSE.

IRH(:,:)=1

!-- define RH index from "clear-sky" (not yet!) relative humidity
DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    DO JTAB=1,12
      IF (PRHCL(JL,JK)*100._JPRB > RRHTAB(JTAB)) THEN
        IRH(JL,JK)=JTAB
      ENDIF
    ENDDO
  ENDDO
ENDDO
!

ZAEROD(KIDIA:KFDIA,1:NWLID,1:KLEV)=0._JPRB
ZBSC  (KIDIA:KFDIA,1:NWLID,1:KLEV)=0._JPRB
ZEXT  (KIDIA:KFDIA,1:NWLID,1:KLEV)=0._JPRB
ZBSCMOL(KIDIA:KFDIA,1:NWLID,1:KLEV)=0._JPRB

!-- optical thickness Rayleigh, NO2, O3, O2, CO2 -------
ZOTRAY(KIDIA:KFDIA,1:NWLID,1:KLEV)=0._JPRB
ZOTNO2(KIDIA:KFDIA,1:NWLID,1:KLEV)=0._JPRB
ZOTO3 (KIDIA:KFDIA,1:NWLID,1:KLEV)=0._JPRB
ZOTO2 (KIDIA:KFDIA,1:NWLID,1:KLEV)=0._JPRB
ZOTCO2(KIDIA:KFDIA,1:NWLID,1:KLEV)=0._JPRB

ZPATHT(KIDIA:KFDIA,1:NWLID,0:KLEV)=0._JPRB
ZPATHS(KIDIA:KFDIA,1:NWLID,0:KLEV)=0._JPRB

!*    --------------------------------------------------------------------------

!*         2.      RAYLEIGH AND GASEOUS OPTICAL THICKNESSES LAYER-BY-LAYER
!                  -------------------------------------------------------

DO JWL=1,NWLID
!-- RWLID in m, ZWLUM in um, ZWLCM in cm!!!!
  ZWLUM=RWLID(JWL)*1.0E+06_JPRB
  ZWLCM=RWLID(JWL)*1.0E+02_JPRB
  ZWL2CM=ZWLCM*ZWLCM
  ZWL4CM=ZWL2CM*ZWL2CM

  ZWN=1._JPRB/ZWLUM
  ZWN2=ZWN*ZWN

!* Air Refractive Index: Edlen, 1966, Metrologia, 2, 71-80 w pw=0
  ZA1=130._JPRB-ZWN2
  ZA2=38.9_JPRB-ZWN2
  ZA3=2406030._JPRB/ZA1
  ZA4=15997._JPRB/ZA2
  ZAN=1._JPRB+(8342.13_JPRB+ZA3+ZA4)*1.E-08_JPRB
  ZAN2=ZAN*ZAN
  ZAA=(24._JPRB*RPI**3) * ((ZAN2-1._JPRB)/(ZAN2+2._JPRB))**2 * (6._JPRB+3._JPRB*RLIDELT)/(6._JPRB-7._JPRB*RLIDELT)

!-- According to Bohdaine et al. (1999), T and p dependence of the molecular density Ns and of the function of 
!   the refractive index (ns**2-1)/(ns**2+2) cancels out, and Rayleigh scattering coefficient ZSIGMA is independent 
!   of T and p; wavelength MUST be in cm to get ZSIGMA in cm2. 
!   ZTAUR for the entire atmosphere should be between 0.59 at 355 nm and < 0.0007 at 1064 nm.

  ZSIGMA=ZAA/(ZWL4CM*RLINS*RLINS)
  ZTAUR =ZSIGMA*101325._JPRB*RNAVO/(RMD*RG)/10._JPRB

  ZOTRAY(KIDIA:KFDIA,JWL,1:KLEV)=0._JPRB
  ZOTNO2(KIDIA:KFDIA,JWL,1:KLEV)=0._JPRB
  ZOTO3 (KIDIA:KFDIA,JWL,1:KLEV)=0._JPRB

  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      ZDP=PDP(JL,JK)

!-- ZFACT is the number of gas molecules in a layer of pressure ZDP at temperature T
      ZFACT=RNAVO*(PDP(JL,JK)*RLITREF)/(PT(JL,JK)*RLIPREF)

!-- Rayleigh: Bodhaine et al., 1999: JAOT, 16, 1854-1861.   tau= sigma*DelP*Avo / (mmw_air * g)
      ZVCO2=PCO2(JL,JK)*RMD/RMCO2
      ZMD = 28.9595_JPRB + 15.0556_JPRB * ZVCO2
      ZMD = RMD
      ZTAUR=ZSIGMA*PDP(JL,JK)*RNAVO/(ZMD*RG)/10._JPRB
      ZOTRAY(JL,JWL,JK) = ZTAUR

!-- O3 
      ZO3V=PO3(JL,JK)*RMD/RMO3                              ! volume mixing ratio in layer
      ZO3MOL=ZO3V*ZFACT                                     ! number of O3 molecules in layer (pressure unit: kg m-1 s-2)
      ZO3TAU=1.E-04_JPRB*RLIO3(3,JWL)*ZO3MOL/(RMO3*RG)      ! cross-section from cm2 molec-1 to m-2 molec-1
      ZOTO3(JL,JWL,JK) = ZO3TAU

!-- NO2
      ZNO2V=PNO2(JL,JK)*RMD/RMNO2                           ! volume mixing ratio in layer
      ZNO2MOL=ZNO2V*ZFACT                                   ! number of NO2 molecules in layer (pressure unit: kg m-1 s-2)
      ZNO2TAU=1.E-04_JPRB*RLINO2(3,JWL)*ZNO2MOL/(RMNO2*RG)  ! cross-section from cm2 molec-1 to m-2 molec-1
      ZOTNO2(JL,JWL,JK) = ZNO2TAU
    ENDDO
  ENDDO
ENDDO


!-- to start with (and to be replaced by more proper quantities once debugged!)
ZMCO2= 357.E-06_JPRB*RMCO2/RMD
ZMO2 = 0.20947_JPRB *ZO2  /RMD
   
IWL=3
ZOTCO2(KIDIA:KFDIA,1:NWLID,1:KLEV)=0._JPRB
ZOTO2 (KIDIA:KFDIA,1:NWLID,1:KLEV)=0._JPRB
DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    ZTE=PT(JL,JK)-RLIT0
    ZTE2=ZTE*ZTE
    ZDPG=PDP(JL,JK)/RG

!-- CO2 -------------------------------
    ZPHIC=EXP(RLICO2(3)*ZTE+RLICO2(4)*ZTE2)
    ZPSIC=EXP(RLICO2(5)*ZTE+RLICO2(6)*ZTE2)    
    ZUCO2 =ZPHIC*ZMCO2*ZDPG
    ZUPCO2=ZPSIC*ZMCO2*ZDPG*PAP(JL,JK)/RLIPREF

!-- Goody model for CO2
    ZGCO2=RLICO2(1)*ZUCO2/SQRT(1._JPRB+RLICO2(1)*ZUCO2*ZUCO2/(RLICO2(2)*ZUPCO2))
!-- Malkmus model for CO2
!    ZM1CO2=RLICO2(2)*ZUPCO2/(2._JPRB*ZUCO2)
!    ZM2CO2=1._JPRB+4._JPRB*RLICO2(1)*ZUCO2*ZUCO2/(RLICO2(2)*ZUPCO2)
!    ZMCO2=ZM1CO2*(SQRT(ZM2CO2)-1._JPRB)

    ZOTCO2(JL,IWL,JK) = ZGCO2

!-- O2 --------------------------------
    ZPHIO=EXP(RLIO2(3) *ZTE+RLIO2(4) *ZTE2)
    ZPSIO=EXP(RLIO2(5) *ZTE+RLIO2(6) *ZTE2)
    ZUO2 =ZPHIO*ZMO2*ZDPG
    ZUPO2=ZPSIO*ZMO2*ZDPG*PAP(JL,JK)/RLIPREF
!-- Goody model for O2
    ZGO2=RLIO2(1)*ZUO2/SQRT(1._JPRB+RLIO2(1)*ZUO2*ZUO2/(RLIO2(2)*ZUPO2))
!-- Malkmus model for O2
!    ZM1O2=RLIO2(2)*ZUPO2/(2._JPRB*ZUO2)
!    ZM2O2=1._JPRB+4._JPRB*RLIO2(1)*ZUO2*ZUO2/(RLIO2(2)*ZUPO2)
!    ZMO2=ZM1O2*(SQRT(ZM2O2)-1._JPRB)

    ZOTO2(JL,IWL,JK) = ZGO2
  ENDDO
ENDDO

!*    --------------------------------------------------------------------------

!*         3.      AEROSOL-RELATED QUANTITIES
!                  --------------------------

DO JAER=1,NACTAERO
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      ZAERMSS(JL,JK,JAER) = MAX(ZEPSAER,PAERO(JL,JK,JAER))*(PAPH(JL,JK)-PAPH(JL,JK-1))/RG 
    ENDDO
  ENDDO
ENDDO

!-- path from TOA
ZPATHS(KIDIA:KFDIA,1:NWLID,0:KLEV)=0._JPRB
ZPATHT(KIDIA:KFDIA,1:NWLID,0:KLEV)=0._JPRB

DO JWL=1,NWLID
  IWAVL=JWLID(JWL)

  DO JAER=1,NACTAERO
    ITYP=YAERO_DESC(JAER)%NTYP
    IBIN=YAERO_DESC(JAER)%NBIN

!-- ITYP is the aerosol type 1:SS,   2:DD,   3:OM    4:BC,   5:SU,   6:FA,   7:BS,   8=VS,
!   IBIN is the bin index: 1-3:SS, 1-3:DD,   2:OM,   2:BC,   2:SU,   1:FA,   1:SB    1=VS,
!   N.B.: extinction coefficients are in m2 g-1

    DO JK=1,KLEV
      DO JL=KIDIA,KFDIA
        IIRH=IRH(JL,JK)
        ZFAC = 1.0_JPRB
        IF (ITYP == 1) THEN
          ZALF(JAER)= ALF_SS(IIRH,IWAVL,IBIN)
          ZOMG(JAER)= OMG_SS(IIRH,IWAVL,IBIN)
          ZLIR(JAER)=RALI_SS(IIRH,IWAVL,IBIN)
          ZFAC = RSS_RH80_MASSFAC
        ELSEIF (ITYP == 2) THEN
          ZALF(JAER)= ALF_DD(IBIN,IWAVL)
          ZOMG(JAER)= OMG_DD(IBIN,IWAVL)
          ZLIR(JAER)=RALI_DD(IBIN,IWAVL)
        ELSEIF (ITYP == 3) THEN
          IF (IBIN == 2) IIRH=IEFRH
          ZALF(JAER)= ALF_OM(IIRH,IWAVL)
          ZOMG(JAER)= OMG_OM(IIRH,IWAVL)
          ZLIR(JAER)=RALI_OM(IIRH,IWAVL)
        ELSEIF (ITYP == 4) THEN
          ZALF(JAER)= ALF_BC(IWAVL)
          ZOMG(JAER)= OMG_BC(IWAVL)
          ZLIR(JAER)=RALI_BC(IWAVL)
        ELSEIF (ITYP == 5 .OR. ITYP == 10) THEN
!          ZALF(JAER)=ALF_SU(IIRH,IWAVL)/RMSO4
          ZALF(JAER)= ALF_SU(IIRH,IWAVL)
          ZOMG(JAER)= OMG_SU(IIRH,IWAVL)
          ZLIR(JAER)=RALI_SU(IIRH,IWAVL)
!-- SO2 does not contribute to optical depth, only SO4 does.
          IF (IBIN == 2) THEN
            ZALF(JAER)=0._JPRB
            ZOMG(JAER)=0._JPRB
            ZLIR(JAER)=0._JPRB
          ENDIF
        ELSEIF (ITYP == 6) THEN
          ZALF(JAER)= ALF_NI(IIRH,IWAVL,IBIN)
          ZOMG(JAER)= OMG_NI(IIRH,IWAVL,IBIN)
          ZLIR(JAER)= RALI_NI(IIRH,IWAVL,IBIN)
        ELSEIF (ITYP == 7) THEN
          ZALF(JAER)= ALF_AM(IIRH,IWAVL)
          ZOMG(JAER)= OMG_AM(IIRH,IWAVL)
          ZLIR(JAER)= RALI_AM(IIRH,IWAVL)
        ELSEIF (ITYP == 8) THEN
          IF (IBIN < 3) THEN
            ZALF(JAER)= ALF_SOA(IIRH,IWAVL,IBIN)
            ZOMG(JAER)= OMG_SOA(IIRH,IWAVL,IBIN)
            ZLIR(JAER)= RALI_SOA(IIRH,IWAVL,IBIN)
          ELSE
            ZALF(JAER)= 0._JPRB
            ZOMG(JAER)= 0._JPRB
            ZLIR(JAER)= 0._JPRB
          ENDIF
        ELSEIF (ITYP ==9) THEN
          IF (NVOLOPTP == 1) THEN
            IIRH=IEFRH
            ZALF(JAER)= ALF_SU(IIRH, IWAVL)
            ZOMG(JAER)= OMG_SU(IIRH, IWAVL)
            ZLIR(JAER)=RALI_SU(IIRH, IWAVL)
          ELSEIF (NVOLOPTP == 2) THEN
            ZALF(JAER)= ALF_BC(IWAVL)
            ZOMG(JAER)= OMG_BC(IWAVL)
            ZLIR(JAER)=RALI_BC(IWAVL)
          ELSEIF (NVOLOPTP == 3) THEN
            ZALF(JAER)= ALF_DD(3,IWAVL)
            ZOMG(JAER)= OMG_DD(3,IWAVL)
            ZLIR(JAER)=RALI_DD(3,IWAVL)
          ENDIF
        ENDIF

!        IF (NAERLISI == 1 .AND. ZLIR(JAER) /= 0._JPRB) THEN
        IF (ZLIR(JAER) /= 0._JPRB) THEN
          ZLIDRAT=1._JPRB/ZLIR(JAER)
        ELSE
          ZLIDRAT=ZLIR(JAER)
        ENDIF

!-- total aerosol optical depth within a layer (sum on aerosol contribution)
        ZAEROD(JL,JWL,JK) = ZAEROD(JL,JWL,JK) + ZAERMSS(JL,JK,JAER) * ZFAC * 1.E+03_JPRB * ZALF(JAER)

!-- back-scatter coefficient computed from:
!    1/ the extinction coefficient multiplied by the lidar ratio  (= backscatter / extinction)
!    2/ this extinction coefficient being a layer integral will have to be divided by the geometrical thickness
        ZBSC(JL,JWL,JK) = ZBSC(JL,JWL,JK) + ZAERMSS(JL,JK,JAER) * ZFAC * 1.E+03_JPRB * ZALF(JAER) * ZLIDRAT

      ENDDO
    ENDDO
  ENDDO

  ZBOLTZ = 1.38E-23_JPRB  !  (J K-1, Boltzmann constant)
!* calculation of molecular backscatter
  ZWLUM=RWLID(JWL)*1.E+6_JPRB      !   wavelength in um
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
! Molecular backscatter from Collis and Russel (1976)
      ZBSCMOL(JL,JWL,JK)=5.45E-32_JPRB * PAP(JL,JK) / (ZBOLTZ * PT(JL,JK)) * (ZWLUM/0.55_JPRB)**(-4.09_JPRB)
    ENDDO
  ENDDO
    
!ENDDO

!*    --------------------------------------------------------------------------

!*         4.      COMPUTATION OF LIDAR SIGNAL
!                  ---------------------------


!ZPATHS(KIDIA:KFDIA,1:NWLID,0:KLEV)=0._JPRB
!ZPATHT(KIDIA:KFDIA,1:NWLID,0:KLEV)=0._JPRB
!DO JWL=1,NWLID
  ZPATHS(KIDIA:KFDIA,JWL,KLEV)=0._JPRB
  ZPATHT(KIDIA:KFDIA,JWL,0   )=0._JPRB


!*         4.1     TOP-DOWN (SATELLITE VIEW)
!                  -------------------------

  ZCAER(:)=0._JPRB
  ZCRAY(:)=0._JPRB
  ZCOZO(:)=0._JPRB
  ZCNO2(:)=0._JPRB
  ZCOXY(:)=0._JPRB
  ZCCO2(:)=0._JPRB
  DO JK=1,KLEV
    IK=KLEV-JK
!DEC$ IVDEP
    DO JL=KIDIA,KFDIA
      ZCAER(JL)=ZCAER(JL)+ZAEROD(JL,JWL,JK)
      ZCRAY(JL)=ZCRAY(JL)+ZOTRAY(JL,JWL,JK)
      ZCOZO(JL)=ZCOZO(JL)+ZOTO3 (JL,JWL,JK)
      ZCNO2(JL)=ZCNO2(JL)+ZOTNO2(JL,JWL,JK)
      ZCOXY(JL)=ZCOXY(JL)+ZOTO2 (JL,JWL,JK)
      ZCCO2(JL)=ZCCO2(JL)+ZOTCO2(JL,JWL,JK)

!-- unattenuated lidar back-scattering coefficients
!--   Q: should RLICLS be included here or not? It's currently 1.0 anyway
!--      so that's a moot point.
      PLISIM(JL,JWL,JK) = ZBSCMOL(JL,JWL,JK)
      PLISIA(JL,JWL,JK) = ZBSC(JL,JWL,JK)/PDZ(JL,JK)

!-- total optical depth (aerosols + gases) from TOA to level JK
      ZPATHT(JL,JWL,JK) = ZPATHT(JL,JWL,JK-1) + 2._JPRB * ( ZAEROD(JL,JWL,JK)&
      & +ZOTRAY(JL,JWL,JK)+ZOTO3(JL,JWL,JK)+ZOTNO2(JL,JWL,JK)+ZOTO2(JL,JWL,JK)+ZOTCO2(JL,JWL,JK) )
      ZUTOT=MIN(200._JPRB,ZPATHT(JL,JWL,JK))
!-- attenuation seen from TOA to level JK
      ZTRTOT=EXP(-ZUTOT)

!-- lidar signal seen from TOA
      PLISIT(JL,JWL,JK)=RLICLS*( ZBSC(JL,JWL,JK)/PDZ(JL,JK)+ZBSCMOL(JL,JWL,JK) )*ZTRTOT


!*         4.2     BOTTOM-UP (GROUND STATION VIEW)
!                  -------------------------------

!-- lidar signal seen from the surface
      ZPATHS(JL,JWL,IK) = ZPATHS(JL,JWL,IK+1) + 2._JPRB * (ZAEROD(JL,JWL,IK+1)&
      & +ZOTRAY(JL,JWL,IK+1)+ZOTO3(JL,JWL,IK+1)+ZOTNO2(JL,JWL,IK+1)+ZOTO2(JL,JWL,IK+1)+ZOTCO2(JL,JWL,IK+1))
      ZUSOT=MIN(200._JPRB,ZPATHS(JL,JWL,IK))

!-- attenuation seen from surface to level IK+1
      ZTRSOT=EXP(-ZUSOT)

!-- lidar signal seen from the surface
      PLISIS(JL,JWL,IK)=RLICLS*( ZBSC(JL,JWL,IK+1)/PDZ(JL,IK+1)+ZBSCMOL(JL,JWL,IK+1) )*ZTRSOT

    ENDDO
  ENDDO
ENDDO

!CALL ABOR1("AER_LIDSIM: END OF ROUTINE")

!     ------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('AER_LIDSIM',1,ZHOOK_HANDLE)
END SUBROUTINE AER_LIDSIM
