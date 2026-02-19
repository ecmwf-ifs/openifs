! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE RADVIS &
  &( YDERAD,KIDIA  , KFDIA  , KLON  , KLEV , &
  &  PRSF1  , PAP    , PQL, PQI, PQR, PQS  , PTP  , PTSPHY, &
  &  PCLAERS, PPRAERS, &
  &  PVISICL, PVISIPR, PVISCAE, PVISPAE, PVISCLD &
  & )

!**** *RADVIS* - ROUTINE COMPUTING THE VISIBILITY

!      J.-J. MORCRETTE , ECMWF


!**   INTERFACE.
!     ----------
!          *RADVISD* IS CALLED FROM *CALLPAR*.

! INPUTS:
! -------

! OUTPUTS:
! --------

!     MODIFICATIONS.
!     -------------
!          Original: JJMorcrette, 20101125  

!-----------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM , JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK

USE YOMCST   , ONLY : RD  , RG
USE YOERAD   , ONLY : TERAD
USE YOESRTCOP, ONLY : RSASWA, RSASWB, RSFUA0, RSFUA1

IMPLICIT NONE

TYPE(TERAD)       ,INTENT(INOUT):: YDERAD
INTEGER(KIND=JPIM),INTENT(IN) :: KIDIA, KFDIA, KLON, KLEV

REAL(KIND=JPRB),INTENT(IN)    :: PRSF1(KLON,KLEV)
REAL(KIND=JPRB),INTENT(IN)    :: PAP(KLON,KLEV)  , PTP(KLON,KLEV)
REAL(KIND=JPRB),INTENT(IN)    :: PQL(KLON,KLEV)
REAL(KIND=JPRB),INTENT(IN)    :: PQI(KLON,KLEV)
REAL(KIND=JPRB),INTENT(IN)    :: PQR(KLON,KLEV)
REAL(KIND=JPRB),INTENT(IN)    :: PQS(KLON,KLEV)
REAL(KIND=JPRB),INTENT(IN)    :: PCLAERS(KLON)   , PPRAERS(KLON)
REAL(KIND=JPRB),INTENT(IN)    :: PTSPHY

REAL(KIND=JPRB),INTENT(OUT)   :: PVISICL(KLON), PVISIPR(KLON)    
REAL(KIND=JPRB),INTENT(OUT)   :: PVISCAE(KLON), PVISPAE(KLON), PVISCLD(KLON)

!-- Local variables

INTEGER(KIND=JPIM) :: JL, JSW
REAL(KIND=JPRB)    :: ZDENSVIS, ZGDT , ZRANGE, ZRELRA , ZDESIC , ZQIWP  , ZQLWP
REAL(KIND=JPRB)    :: Z1CLD   , ZQRWP, ZQSWP , ZSNOICE, ZCFSNIC, ZLIQRAI, ZCFLIRA
REAL(KIND=JPRB)    :: ZNS     , ZSIGAIR
REAL(KIND=JPRB)    :: ZCONDEN(KLON), ZRAYL(KLON), ZRHO(KLON)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('RADVIS',0,ZHOOK_HANDLE)
ASSOCIATE(RNS=>YDERAD%RNS, RSIGAIR=>YDERAD%RSIGAIR)
!-- if prognostic aerosols, take the scattering cross-section from AER_BDGTMSS
!-- if not, take the scattering cross-section from RADACA (climatological aerosols)


!-- VISIBILITY CALCULATIONS
!-- formulation from Hinkley (1976) Vm = (3.91/Sigma) * (0.55/Lambda)^1.3
!   is here used to get visibility at 0.55 um, where Sigma is the total extinction coefficient in km-1
!   0.012 km-1 is the approximate contribution from Rayleigh (could be improved)
!   ZSIGAER is the extinction coefficient in km-1 from the aerosols (computed in aer_bdgtmss.F90)
!   ZCONDEN is the extinction coefficient, sum of the effect of rain and snow 

ZGDT=RG*PTSPHY
ZRANGE = 3.912023_JPRB
ZRELRA=1000._JPRB
ZDESIC=1000._JPRB
! values for 0.55 um (in fact 0.6250 - 0.4415 um)
JSW=10                
DO JL=KIDIA,KFDIA
  ZRHO(JL)=PRSF1(JL,KLEV)/(RD*PTP(JL,KLEV))
  ZDENSVIS=ZRHO(JL)
! All ..WP in kg kg-1
  ZQIWP = 0._JPRB
  ZQLWP = 0._JPRB
  IF (PAP(JL,KLEV) > 0.001_JPRB) THEN
    Z1CLD   = 1._JPRB/PAP(JL,KLEV)
    ZQIWP   = MAX(0._JPRB, PQI(JL,KLEV) * Z1CLD)
    ZQLWP   = MAX(0._JPRB, PQL(JL,KLEV) * Z1CLD)
  ENDIF
  ZQRWP   = MAX(0._JPRB, PQR(JL,KLEV))
  ZQSWP   = MAX(0._JPRB, PQS(JL,KLEV))
! ice and snow together
  ZSNOICE = ZQIWP+ZQSWP
  ZCFSNIC = RSFUA0(JSW) + RSFUA1(JSW) / ZDESIC
! liquid and rain together
  ZLIQRAI = ZQLWP+ZQRWP
  ZCFLIRA = RSASWA(JSW) + RSASWB(JSW) / ZRELRA
! cloud+precipitation (1.E+6 is 1.E+3 to go from kg to g, and 1.E+3 from m-1 to km-1)
  ZCONDEN(JL) = (ZSNOICE * ZCFSNIC + ZLIQRAI * ZCFLIRA) * ZDENSVIS * 1.E+06_JPRB
! Molecular density in first layer
  ZNS = RNS*273.15_JPRB/PTP(JL,KLEV)
! Rayleigh scattering cross section (cm2 molec-1)
  ZSIGAIR = (RSIGAIR/ZNS)/ZNS ! not RSIGAIR/(ZNS*ZNS) because intermediate value busts single precision
! Rayleigh scattering coefficient (km-1)
  ZRAYL(JL) = 1.E+05_JPRB * ZSIGAIR * ZNS * PAP(JL,KLEV) / 101325._JPRB
!-- to follow the Met Office (Clark et al., 2008, QJ, 134, 1801; Claxton, 2008, QJ, 134, 1527) 
  ZRAYL(JL) = 0.03912023_JPRB
! "total" visibility, visi.Rayl+aerosols, visi.Rayl+cloud/precip
  PVISICL(JL)  = ZRANGE / (ZRAYL(JL) + PCLAERS(JL) + ZCONDEN(JL))
  PVISIPR(JL)  = ZRANGE / (ZRAYL(JL) + PPRAERS(JL) + ZCONDEN(JL))
  PVISCAE(JL)  = ZRANGE / (ZRAYL(JL) + PCLAERS(JL))
  PVISPAE(JL)  = ZRANGE / (ZRAYL(JL) + PPRAERS(JL))
  PVISCLD(JL)  = ZRANGE / (ZRAYL(JL) + ZCONDEN(JL))
ENDDO

!-----------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('RADVIS',1,ZHOOK_HANDLE)
END SUBROUTINE RADVIS
