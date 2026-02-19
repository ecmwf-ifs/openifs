#ifdef RS6K
@PROCESS HOT NOSTRICT
#endif
SUBROUTINE SRTM_SPCVRT&
 & (YDDIMV, YDERDI,KIDIA   , KFDIA    , KLEV    , KSW    , PONEMINUS,&
 & PALBD   , PALBP,&
 & PFRCL   , PTAUC   , PASYC  , POMGC  , PTAUA   , PASYA   , POMGA , PRMU0,&
 & KLAYTROP,&
 & PCOLCH4  , PCOLCO2 , PCOLH2O , PCOLMOL  , PCOLO2 , PCOLO3 ,&
 & PFORFAC , PFORFRAC , KINDFOR , PSELFFAC, PSELFFRAC, KINDSELF ,&
 & PFAC00  , PFAC01   , PFAC10  , PFAC11 ,&
 & KJP     , KJT      , KJT1 ,&
 !-- output arrays 
 &   PBBFD   , PBBFU    , PBBCD   , PBBCU  , PSUDU , PSUDUC, &
 &   PBBFDIR , PBBCDIR )  

!**** *SRTM_SPCVRT* - SPECTRAL LOOP TO COMPUTE THE SHORTWAVE RADIATION FLUXES.

!     PURPOSE.
!     --------

!          THIS ROUTINE COMPUTES THE TWO-STREAM METHOD OF BARKER

!**   INTERFACE.
!     ----------

!          *SRTM_SPCVRT* IS CALLED FROM *SRTM_SRTM_224GP*

!        IMPLICIT ARGUMENTS :
!        --------------------

!     ==== INPUTS ===
!     ==== OUTPUTS ===

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!          *SWVRTQDR*

!     REFERENCE.
!     ----------

!        SEE RADIATION'S PART OF THE ECMWF RESEARCH DEPARTMENT
!        DOCUMENTATION
!     AUTHOR.
!     -------
!        from Howard Barker
!        JEAN-JACQUES MORCRETTE  *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 03-02-27
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        JJMorcrette 20070614 bug-fix for solar duration
!        D.Salmond  31-Oct-2007 Vector version in the style of RRTM from Meteo France & NEC
!        JJMorcrette/MJIacono 20080724 Look-up table replacing exponential
!        Y. Bouteloup 21 October 2016 : Add new output fluxes : Solar Downward 
!     ------------------------------------------------------------------

USE YOMDIMV  , ONLY : TDIMV
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE PARSRTM  , ONLY : JPB1, JPB2
USE YOESRTM  , ONLY : JPGPT
USE YOESRTWN , ONLY : NGC
USE YOERDI   , ONLY : TERDI
USE YOESRTAB , ONLY : BPADE, TRANS, RODLOW, RTBLINT

IMPLICIT NONE

!     ------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

TYPE(TDIMV)       ,INTENT(IN)    :: YDDIMV
TYPE(TERDI)       ,INTENT(INOUT) :: YDERDI
INTEGER(KIND=JPIM),INTENT(IN)    :: KSW 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA, KFDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV
REAL(KIND=JPRB)   ,INTENT(IN)    :: PONEMINUS(KIDIA:KFDIA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PALBD(KIDIA:KFDIA,KSW) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PALBP(KIDIA:KFDIA,KSW) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFRCL(KIDIA:KFDIA,YDDIMV%NFLEVG)     ! bottom to top
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTAUC(KIDIA:KFDIA,YDDIMV%NFLEVG,KSW) ! bottom to top
REAL(KIND=JPRB)   ,INTENT(IN)    :: PASYC(KIDIA:KFDIA,YDDIMV%NFLEVG,KSW) ! bottom to top
REAL(KIND=JPRB)   ,INTENT(IN)    :: POMGC(KIDIA:KFDIA,YDDIMV%NFLEVG,KSW) ! bottom to top
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTAUA(KIDIA:KFDIA,YDDIMV%NFLEVG,KSW) ! bottom to top
REAL(KIND=JPRB)   ,INTENT(IN)    :: PASYA(KIDIA:KFDIA,YDDIMV%NFLEVG,KSW) ! bottom to top
REAL(KIND=JPRB)   ,INTENT(IN)    :: POMGA(KIDIA:KFDIA,YDDIMV%NFLEVG,KSW) ! bottom to top
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRMU0(KIDIA:KFDIA)
INTEGER(KIND=JPIM),INTENT(IN)    :: KLAYTROP(KIDIA:KFDIA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCOLCH4(KIDIA:KFDIA,YDDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCOLCO2(KIDIA:KFDIA,YDDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCOLH2O(KIDIA:KFDIA,YDDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCOLMOL(KIDIA:KFDIA,YDDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCOLO2(KIDIA:KFDIA,YDDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCOLO3(KIDIA:KFDIA,YDDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFORFAC(KIDIA:KFDIA,YDDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFORFRAC(KIDIA:KFDIA,YDDIMV%NFLEVG)
INTEGER(KIND=JPIM),INTENT(IN)    :: KINDFOR(KIDIA:KFDIA,YDDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSELFFAC(KIDIA:KFDIA,YDDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSELFFRAC(KIDIA:KFDIA,YDDIMV%NFLEVG)
INTEGER(KIND=JPIM),INTENT(IN)    :: KINDSELF(KIDIA:KFDIA,YDDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFAC00(KIDIA:KFDIA,YDDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFAC01(KIDIA:KFDIA,YDDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFAC10(KIDIA:KFDIA,YDDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFAC11(KIDIA:KFDIA,YDDIMV%NFLEVG)
INTEGER(KIND=JPIM),INTENT(IN)    :: KJP(KIDIA:KFDIA,YDDIMV%NFLEVG)
INTEGER(KIND=JPIM),INTENT(IN)    :: KJT(KIDIA:KFDIA,YDDIMV%NFLEVG)
INTEGER(KIND=JPIM),INTENT(IN)    :: KJT1(KIDIA:KFDIA,YDDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PBBFD(KIDIA:KFDIA,YDDIMV%NFLEVG+1) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PBBFU(KIDIA:KFDIA,YDDIMV%NFLEVG+1) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PBBCD(KIDIA:KFDIA,YDDIMV%NFLEVG+1) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PBBCU(KIDIA:KFDIA,YDDIMV%NFLEVG+1) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSUDU(KIDIA:KFDIA), PSUDUC(KIDIA:KFDIA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PBBFDIR(KIDIA:KFDIA,YDDIMV%NFLEVG+1)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PBBCDIR(KIDIA:KFDIA,YDDIMV%NFLEVG+1)
!     ------------------------------------------------------------------

!              ------------

LOGICAL :: LLRTCHK(KIDIA:KFDIA,YDDIMV%NFLEVG)

REAL(KIND=JPRB) ::  ZCLEAR  , ZCLOUD
REAL(KIND=JPRB) ::&
 &   ZDBT(KIDIA:KFDIA,YDDIMV%NFLEVG+1)&
 &   , ZGCC(KIDIA:KFDIA,YDDIMV%NFLEVG)   , ZGCO(KIDIA:KFDIA,YDDIMV%NFLEVG)&
 &   , ZOMCC(KIDIA:KFDIA,YDDIMV%NFLEVG)  , ZOMCO(KIDIA:KFDIA,YDDIMV%NFLEVG)&
 &   , ZRDND(KIDIA:KFDIA,YDDIMV%NFLEVG+1), ZRDNDC(KIDIA:KFDIA,YDDIMV%NFLEVG+1)&
 &   , ZREF(KIDIA:KFDIA,YDDIMV%NFLEVG+1) , ZREFC(KIDIA:KFDIA,YDDIMV%NFLEVG+1) , ZREFO(KIDIA:KFDIA,YDDIMV%NFLEVG+1)&
 &   , ZREFD(KIDIA:KFDIA,YDDIMV%NFLEVG+1), ZREFDC(KIDIA:KFDIA,YDDIMV%NFLEVG+1), ZREFDO(KIDIA:KFDIA,YDDIMV%NFLEVG+1)&
 &   , ZRUP(KIDIA:KFDIA,YDDIMV%NFLEVG+1) , ZRUPD(KIDIA:KFDIA,YDDIMV%NFLEVG+1)&
 &   , ZRUPC(KIDIA:KFDIA,YDDIMV%NFLEVG+1), ZRUPDC(KIDIA:KFDIA,YDDIMV%NFLEVG+1)&
 &   , ZTAUC(KIDIA:KFDIA,YDDIMV%NFLEVG)  , ZTAUO(KIDIA:KFDIA,YDDIMV%NFLEVG)&
 &   , ZTDBT(KIDIA:KFDIA,YDDIMV%NFLEVG+1)&
 &   , ZTRA(KIDIA:KFDIA,YDDIMV%NFLEVG+1) , ZTRAC(KIDIA:KFDIA,YDDIMV%NFLEVG+1) , ZTRAO(KIDIA:KFDIA,YDDIMV%NFLEVG+1)&
 & , ZTRAD(KIDIA:KFDIA,YDDIMV%NFLEVG+1), ZTRADC(KIDIA:KFDIA,YDDIMV%NFLEVG+1), ZTRADO(KIDIA:KFDIA,YDDIMV%NFLEVG+1)   
REAL(KIND=JPRB) ::&
 & ZDBTC(KIDIA:KFDIA,YDDIMV%NFLEVG+1), ZTDBTC(KIDIA:KFDIA,YDDIMV%NFLEVG+1), ZINCFLX(KIDIA:KFDIA,JPGPT)&
 & ,  ZINCF14(KIDIA:KFDIA,14)   , ZINCTOT(KIDIA:KFDIA)   

INTEGER(KIND=JPIM) :: IB1, IB2, IBM, IGT, IKL, IW(KIDIA:KFDIA), JB, JG, JK, I_KMODTS, JL

REAL(KIND=JPRB) :: ZDBTMC, ZDBTMO, ZF, ZINCFLUX(KIDIA:KFDIA), ZWF

!-- Output of SRTM_TAUMOLn routines

REAL(KIND=JPRB) :: ZTAUG(KIDIA:KFDIA,YDDIMV%NFLEVG,16), ZTAUR(KIDIA:KFDIA,YDDIMV%NFLEVG,16), ZSFLXZEN(KIDIA:KFDIA,16)

!-- Output of SRTM_VRTQDR routine
REAL(KIND=JPRB) ::&
 & ZCD(KIDIA:KFDIA,YDDIMV%NFLEVG+1,JPGPT), ZCU(KIDIA:KFDIA,YDDIMV%NFLEVG+1,JPGPT)&
 & ,  ZFD(KIDIA:KFDIA,YDDIMV%NFLEVG+1,JPGPT), ZFU(KIDIA:KFDIA,YDDIMV%NFLEVG+1,JPGPT)  

!--  Use of exponential look-up table
REAL(KIND=JPRB) :: ZE1, ZE2, ZTBLIND
INTEGER(KIND=JPIM) :: ITIND

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE


#include "srtm_taumol16.intfb.h"
#include "srtm_taumol17.intfb.h"
#include "srtm_taumol18.intfb.h"
#include "srtm_taumol19.intfb.h"
#include "srtm_taumol20.intfb.h"
#include "srtm_taumol21.intfb.h"
#include "srtm_taumol22.intfb.h"
#include "srtm_taumol23.intfb.h"
#include "srtm_taumol24.intfb.h"
#include "srtm_taumol25.intfb.h"
#include "srtm_taumol26.intfb.h"
#include "srtm_taumol27.intfb.h"
#include "srtm_taumol28.intfb.h"
#include "srtm_taumol29.intfb.h"
#include "srtm_reftra.intfb.h"
#include "srtm_vrtqdr.intfb.h"
!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SRTM_SPCVRT',0,ZHOOK_HANDLE)
ASSOCIATE(NFLEVG=>YDDIMV%NFLEVG, &
 & REPCLC=>YDERDI%REPCLC)
!-- Two-stream model 1: Eddington, 2: PIFM, Zdunkowski et al., 3: discret ordinates
! KMODTS is set in SWREFTRA
!NDBUG=4

IB1=JPB1
IB2=JPB2

DO JL = KIDIA, KFDIA
  IF (PRMU0(JL) > 0.0_JPRB) THEN
    IW(JL)=0
    ZINCFLUX(JL)=0.0_JPRB
    ZINCTOT(JL)=0.0_JPRB
  ENDIF
ENDDO

JB=IB1-1
DO JB = IB1, IB2
  DO JL = KIDIA, KFDIA
    IF (PRMU0(JL) > 0.0_JPRB) THEN
      IBM = JB-15
      IGT = NGC(IBM)
      ZINCF14(JL,IBM)=0.0_JPRB
    ENDIF
  ENDDO

  !-- for each band, computes the gaseous and Rayleigh optical thickness 
  !  for all g-points within the band

  IF (JB == 16) THEN
    CALL SRTM_TAUMOL16&
     & (YDDIMV, KIDIA   , KFDIA    , KLEV    ,&
     & PFAC00  , PFAC01   , PFAC10   , PFAC11   ,&
     & KJP     , KJT      , KJT1     , PONEMINUS,&
     & PCOLH2O , PCOLCH4  , PCOLMOL  ,&
     & KLAYTROP, PSELFFAC , PSELFFRAC, KINDSELF, PFORFAC  , PFORFRAC, KINDFOR ,&
     & ZSFLXZEN, ZTAUG    , ZTAUR    , PRMU0&
     & )  

  ELSEIF (JB == 17) THEN
    CALL SRTM_TAUMOL17&
     & (YDDIMV, KIDIA   , KFDIA   , KLEV    ,&
     & PFAC00  , PFAC01  , PFAC10   , PFAC11 ,&
     & KJP     , KJT     , KJT1     , PONEMINUS ,&
     & PCOLH2O , PCOLCO2 , PCOLMOL  ,&
     & KLAYTROP, PSELFFAC, PSELFFRAC, KINDSELF  , PFORFAC, PFORFRAC, KINDFOR ,&
     & ZSFLXZEN, ZTAUG   , ZTAUR    , PRMU0&
     & )  

  ELSEIF (JB == 18) THEN
    CALL SRTM_TAUMOL18&
     & (YDDIMV, KIDIA   , KFDIA   , KLEV    ,&
     & PFAC00  , PFAC01  , PFAC10   , PFAC11 ,&
     & KJP     , KJT     , KJT1     , PONEMINUS ,&
     & PCOLH2O , PCOLCH4 , PCOLMOL  ,&
     & KLAYTROP, PSELFFAC, PSELFFRAC, KINDSELF  , PFORFAC, PFORFRAC, KINDFOR ,&
     & ZSFLXZEN, ZTAUG   , ZTAUR    , PRMU0&
     & )  

  ELSEIF (JB == 19) THEN
    CALL SRTM_TAUMOL19&
     & (YDDIMV, KIDIA   , KFDIA   , KLEV    ,&
     & PFAC00  , PFAC01  , PFAC10   , PFAC11 ,&
     & KJP     , KJT     , KJT1     , PONEMINUS ,&
     & PCOLH2O , PCOLCO2 , PCOLMOL  ,&
     & KLAYTROP, PSELFFAC, PSELFFRAC, KINDSELF  , PFORFAC, PFORFRAC, KINDFOR ,&
     & ZSFLXZEN, ZTAUG   , ZTAUR    , PRMU0&
     & )  

  ELSEIF (JB == 20) THEN
    CALL SRTM_TAUMOL20&
     & (YDDIMV, KIDIA   , KFDIA   , KLEV    ,&
     & PFAC00  , PFAC01  , PFAC10   , PFAC11 ,&
     & KJP     , KJT     , KJT1     ,&
     & PCOLH2O , PCOLCH4 , PCOLMOL  ,&
     & KLAYTROP, PSELFFAC, PSELFFRAC, KINDSELF  , PFORFAC, PFORFRAC, KINDFOR ,&
     & ZSFLXZEN, ZTAUG   , ZTAUR    , PRMU0&
     & )  

  ELSEIF (JB == 21) THEN
    CALL SRTM_TAUMOL21&
     & (YDDIMV, KIDIA   , KFDIA   , KLEV    ,&
     & PFAC00  , PFAC01  , PFAC10   , PFAC11 ,&
     & KJP     , KJT     , KJT1     , PONEMINUS ,&
     & PCOLH2O , PCOLCO2 , PCOLMOL  ,&
     & KLAYTROP, PSELFFAC, PSELFFRAC, KINDSELF  , PFORFAC, PFORFRAC, KINDFOR ,&
     & ZSFLXZEN, ZTAUG   , ZTAUR    , PRMU0&
     & )  

  ELSEIF (JB == 22) THEN
    CALL SRTM_TAUMOL22&
     & (YDDIMV, KIDIA   , KFDIA   , KLEV    ,&
     & PFAC00  , PFAC01  , PFAC10   , PFAC11 ,&
     & KJP     , KJT     , KJT1     , PONEMINUS ,&
     & PCOLH2O , PCOLMOL , PCOLO2   ,&
     & KLAYTROP, PSELFFAC, PSELFFRAC, KINDSELF  , PFORFAC, PFORFRAC, KINDFOR ,&
     & ZSFLXZEN, ZTAUG   , ZTAUR    , PRMU0&
     & )  

  ELSEIF (JB == 23) THEN
    CALL SRTM_TAUMOL23&
     & (YDDIMV, KIDIA   , KFDIA   , KLEV    ,&
     & PFAC00  , PFAC01  , PFAC10   , PFAC11 ,&
     & KJP     , KJT     , KJT1     ,&
     & PCOLH2O , PCOLMOL ,&
     & KLAYTROP, PSELFFAC, PSELFFRAC, KINDSELF  , PFORFAC, PFORFRAC, KINDFOR ,&
     & ZSFLXZEN, ZTAUG   , ZTAUR    , PRMU0&
     & )  

  ELSEIF (JB == 24) THEN
    CALL SRTM_TAUMOL24&
     & (YDDIMV, KIDIA   , KFDIA   , KLEV    ,&
     & PFAC00  , PFAC01  , PFAC10   , PFAC11 ,&
     & KJP     , KJT     , KJT1     , PONEMINUS ,&
     & PCOLH2O , PCOLMOL , PCOLO2   , PCOLO3 ,&
     & KLAYTROP, PSELFFAC, PSELFFRAC, KINDSELF  , PFORFAC, PFORFRAC, KINDFOR ,&
     & ZSFLXZEN, ZTAUG   , ZTAUR    , PRMU0&
     & )  

  ELSEIF (JB == 25) THEN
    !--- visible 16000-22650 cm-1   0.4415 - 0.6250 um
    CALL SRTM_TAUMOL25&
     & (YDDIMV, KIDIA    , KFDIA   , KLEV     ,&
     & PFAC00   , PFAC01  , PFAC10 , PFAC11 ,&
     & KJP      , KJT     , KJT1   ,&
     & PCOLH2O  , PCOLMOL , PCOLO3 ,&
     & KLAYTROP ,&
     & ZSFLXZEN, ZTAUG   , ZTAUR   , PRMU0&
     & )  

  ELSEIF (JB == 26) THEN
    !--- UV-A 22650-29000 cm-1   0.3448 - 0.4415 um
    CALL SRTM_TAUMOL26&
     & (YDDIMV, KIDIA   , KFDIA   , KLEV    ,&
     & PCOLMOL  ,KLAYTROP,&
     & ZSFLXZEN, ZTAUG   , ZTAUR    , PRMU0&
     & )  

  ELSEIF (JB == 27) THEN
    !--- UV-B 29000-38000 cm-1   0.2632 - 0.3448 um
    CALL SRTM_TAUMOL27&
     & (YDDIMV, KIDIA   , KFDIA   , KLEV    ,&
     & PFAC00  , PFAC01  , PFAC10   , PFAC11 ,&
     & KJP     , KJT     , KJT1     ,&
     & PCOLMOL , PCOLO3 ,&
     & KLAYTROP ,&
     & ZSFLXZEN, ZTAUG   , ZTAUR    , PRMU0&
     & )  

  ELSEIF (JB == 28) THEN
    !--- UV-C 38000-50000 cm-1   0.2000 - 0.2632 um
    CALL SRTM_TAUMOL28&
     & (YDDIMV, KIDIA   , KFDIA   , KLEV   ,&
     & PFAC00  , PFAC01  , PFAC10 , PFAC11 ,&
     & KJP     , KJT     , KJT1   , PONEMINUS ,&
     & PCOLMOL , PCOLO2  , PCOLO3 ,&
     & KLAYTROP ,&
     & ZSFLXZEN, ZTAUG   , ZTAUR  , PRMU0&
     & )  

  ELSEIF (JB == 29) THEN
    CALL SRTM_TAUMOL29&
     & (YDDIMV, KIDIA    , KFDIA   , KLEV     ,&
     & PFAC00   , PFAC01  , PFAC10   , PFAC11 ,&
     & KJP      , KJT     , KJT1     ,&
     & PCOLH2O  , PCOLCO2 , PCOLMOL  ,&
     & KLAYTROP , PSELFFAC, PSELFFRAC, KINDSELF  , PFORFAC, PFORFRAC, KINDFOR ,&
     & ZSFLXZEN , ZTAUG   , ZTAUR    , PRMU0&
     & )  

  ENDIF

  DO JG=1,IGT
    DO JL = KIDIA, KFDIA
      IF (PRMU0(JL) > 0.0_JPRB) THEN
        IW(JL)=IW(JL)+1

        ZINCFLX(JL,IW(JL)) =ZSFLXZEN(JL,JG)*PRMU0(JL)
        ZINCFLUX(JL)    =ZINCFLUX(JL)+ZSFLXZEN(JL,JG)*PRMU0(JL)           
        ZINCTOT(JL)     =ZINCTOT(JL)+ZSFLXZEN(JL,JG)
        ZINCF14(JL,IBM)=ZINCF14(JL,IBM)+ZSFLXZEN(JL,JG)

        !-- CALL to compute layer reflectances and transmittances for direct 
        !  and diffuse sources, first clear then cloudy.
        !   Use direct/parallel albedo for direct radiation and diffuse albedo
        !   otherwise.

        ! ZREFC(JK)  direct albedo for clear
        ! ZREFO(JK)  direct albedo for cloud
        ! ZREFDC(JK) diffuse albedo for clear
        ! ZREFDO(JK) diffuse albedo for cloud
        ! ZTRAC(JK)  direct transmittance for clear
        ! ZTRAO(JK)  direct transmittance for cloudy
        ! ZTRADC(JK) diffuse transmittance for clear
        ! ZTRADO(JK) diffuse transmittance for cloudy

        ! ZREF(JK)   direct reflectance
        ! ZREFD(JK)  diffuse reflectance
        ! ZTRA(JK)   direct transmittance
        ! ZTRAD(JK)  diffuse transmittance

        ! ZDBTC(JK)  clear direct beam transmittance
        ! ZDBTO(JK)  cloudy direct beam transmittance
        ! ZDBT(JK)   layer mean direct beam transmittance
        ! ZTDBT(JK)  total direct beam transmittance at levels

        !-- clear-sky    
        !----- TOA direct beam    
        ZTDBTC(JL,1)=1._JPRB
        !----- surface values
        ZDBTC(JL,KLEV+1) =0.0_JPRB
        ZTRAC(JL,KLEV+1) =0.0_JPRB
        ZTRADC(JL,KLEV+1)=0.0_JPRB
        ZREFC(JL,KLEV+1) =PALBP(JL,IBM)
        ZREFDC(JL,KLEV+1)=PALBD(JL,IBM)
        ZRUPC(JL,KLEV+1) =PALBP(JL,IBM)
        ZRUPDC(JL,KLEV+1)=PALBD(JL,IBM)

        !-- total sky    
        !----- TOA direct beam    
        ZTDBT(JL,1)=1._JPRB
        !----- surface values
        ZDBT(JL,KLEV+1) =0.0_JPRB
        ZTRA(JL,KLEV+1) =0.0_JPRB
        ZTRAD(JL,KLEV+1)=0.0_JPRB
        ZREF(JL,KLEV+1) =PALBP(JL,IBM)
        ZREFD(JL,KLEV+1)=PALBD(JL,IBM)
        ZRUP(JL,KLEV+1) =PALBP(JL,IBM)
        ZRUPD(JL,KLEV+1)=PALBD(JL,IBM)
      ENDIF
    ENDDO

    DO JK=1,KLEV
      DO JL = KIDIA, KFDIA
        IF (PRMU0(JL) > 0.0_JPRB) THEN

          !-- NB: a two-stream calculations from top to bottom, but RRTM_SW quantities 
          !       are given bottom to top (argh!)
          !       Inputs for clouds and aerosols are bottom to top as inputs

          IKL=KLEV+1-JK

          !-- clear-sky optical parameters      
          LLRTCHK(JL,JK)=.TRUE.

          !-- clear-sky optical parameters including aerosols
          ZTAUC(JL,JK) = ZTAUR(JL,IKL,JG) + ZTAUG(JL,IKL,JG) + PTAUA(JL,IKL,IBM)
          ZOMCC(JL,JK) = ZTAUR(JL,IKL,JG)*1.0_JPRB + PTAUA(JL,IKL,IBM)*POMGA(JL,IKL,IBM)
          ZGCC(JL,JK) = PASYA(JL,IKL,IBM)*POMGA(JL,IKL,IBM)*PTAUA(JL,IKL,IBM) / ZOMCC(JL,JK)
          ZOMCC(JL,JK) = ZOMCC(JL,JK) / ZTAUC(JL,JK)

          !-- total sky optical parameters        
          ZTAUO(JL,JK) = ZTAUR(JL,IKL,JG) + ZTAUG(JL,IKL,JG) + PTAUA(JL,IKL,IBM) + PTAUC(JL,IKL,IBM)
          ZOMCO(JL,JK) = PTAUA(JL,IKL,IBM)*POMGA(JL,IKL,IBM) + PTAUC(JL,IKL,IBM)*POMGC(JL,IKL,IBM)&
           & + ZTAUR(JL,IKL,JG)*1.0_JPRB  
          ZGCO(JL,JK) = (PTAUC(JL,IKL,IBM)*POMGC(JL,IKL,IBM)*PASYC(JL,IKL,IBM)&
           & +  PTAUA(JL,IKL,IBM)*POMGA(JL,IKL,IBM)*PASYA(JL,IKL,IBM))&
           & /  ZOMCO(JL,JK)  
          ZOMCO(JL,JK) = ZOMCO(JL,JK) / ZTAUO(JL,JK)

        ENDIF
      ENDDO
    ENDDO

    !-- Delta scaling for clear-sky / aerosol optical quantities
    DO  JK=1,KLEV
      DO JL = KIDIA, KFDIA
        IF (PRMU0(JL) > 0.0_JPRB) THEN
          ZF=ZGCC(JL,JK)*ZGCC(JL,JK)
          ZWF=ZOMCC(JL,JK)*ZF
          ZTAUC(JL,JK)=(1._JPRB-ZWF)*ZTAUC(JL,JK)
          ZOMCC(JL,JK)=(ZOMCC(JL,JK)-ZWF)/(1.0_JPRB-ZWF)
          ZGCC(JL,JK)=(ZGCC(JL,JK)-ZF)/(1.0_JPRB-ZF)
        ENDIF
      ENDDO
    ENDDO

    CALL SRTM_REFTRA(YDDIMV, KIDIA, KFDIA, KLEV, I_KMODTS ,&
     &   LLRTCHK, ZGCC  , PRMU0, ZTAUC , ZOMCC ,&
     &   ZREFC  , ZREFDC, ZTRAC, ZTRADC )  

    !-- Delta scaling for cloudy quantities
    DO JK=1,KLEV
      DO JL = KIDIA, KFDIA
        IF (PRMU0(JL) > 0.0_JPRB) THEN
          IKL=KLEV+1-JK
          LLRTCHK(JL,JK)=.FALSE.
          ZF=ZGCO(JL,JK)*ZGCO(JL,JK)
          ZWF=ZOMCO(JL,JK)*ZF
          ZTAUO(JL,JK)=(1._JPRB-ZWF)*ZTAUO(JL,JK)
          ZOMCO(JL,JK)=(ZOMCO(JL,JK)-ZWF)/(1._JPRB-ZWF)
          ZGCO(JL,JK)=(ZGCO(JL,JK)-ZF)/(1._JPRB-ZF)
          LLRTCHK(JL,JK)=(PFRCL(JL,IKL) > REPCLC)
        ENDIF
      ENDDO
    ENDDO

    CALL SRTM_REFTRA(YDDIMV, KIDIA, KFDIA, KLEV, I_KMODTS ,&
     &   LLRTCHK, ZGCO  , PRMU0, ZTAUO , ZOMCO ,&
     &   ZREFO , ZREFDO, ZTRAO, ZTRADO )  

    DO JK=1,KLEV
      DO JL = KIDIA, KFDIA
        IF (PRMU0(JL) > 0.0_JPRB) THEN

          !-- combine clear and cloudy contributions for total sky

          IKL=KLEV+1-JK 
          ZCLEAR   = 1.0_JPRB - PFRCL(JL,IKL)
          ZCLOUD   = PFRCL(JL,IKL)

          ZREF(JL,JK) = ZCLEAR*ZREFC(JL,JK) + ZCLOUD*ZREFO(JL,JK)
          ZREFD(JL,JK)= ZCLEAR*ZREFDC(JL,JK)+ ZCLOUD*ZREFDO(JL,JK)
          ZTRA(JL,JK) = ZCLEAR*ZTRAC(JL,JK) + ZCLOUD*ZTRAO(JL,JK)
          ZTRAD(JL,JK)= ZCLEAR*ZTRADC(JL,JK)+ ZCLOUD*ZTRADO(JL,JK)

          !-- direct beam transmittance        

!          ZDBTMC     = EXP(-ZTAUC(JL,JK)/PRMU0(JL))
!          ZDBTMO     = EXP(-ZTAUO(JL,JK)/PRMU0(JL))

!-- Use exponential look-up table for transmittance, or expansion of exponential for 
!   low optical thickness
          ZE1 = ZTAUC(JL,JK)/PRMU0(JL)
          IF (ZE1 <= RODLOW) THEN
            ZDBTMC = 1._JPRB - ZE1 + 0.5_JPRB*ZE1*ZE1
          ELSE
            ZTBLIND = ZE1 / (BPADE+ZE1)
            ITIND = RTBLINT * ZTBLIND + 0.5_JPRB
            ZDBTMC = TRANS(ITIND)
          ENDIF

          ZE2 = ZTAUO(JL,JK)/PRMU0(JL)
          IF (ZE2 <= RODLOW) THEN
            ZDBTMO = 1._JPRB - ZE2 + 0.5_JPRB*ZE2*ZE2
          ELSE
            ZTBLIND = ZE2 / (BPADE+ZE2)
            ITIND = RTBLINT * ZTBLIND + 0.5_JPRB
            ZDBTMO = TRANS(ITIND)
          ENDIF
!---
          ZDBT(JL,JK)   = ZCLEAR*ZDBTMC+ZCLOUD*ZDBTMO
          ZTDBT(JL,JK+1)= ZDBT(JL,JK)*ZTDBT(JL,JK)

          !-- clear-sky        
          ZDBTC(JL,JK)   =ZDBTMC
          ZTDBTC(JL,JK+1)=ZDBTC(JL,JK)*ZTDBTC(JL,JK)

        ENDIF
      ENDDO
    ENDDO

    !-- vertical quadrature producing clear-sky fluxes

    CALL SRTM_VRTQDR(YDDIMV, KIDIA, KFDIA, KLEV, IW ,&
     &   ZREFC, ZREFDC, ZTRAC , ZTRADC ,&
     &   ZDBTC, ZRDNDC, ZRUPC , ZRUPDC, ZTDBTC ,&
     &   ZCD  , ZCU   , PRMU0 )  

    !-- vertical quadrature producing cloudy fluxes

    CALL SRTM_VRTQDR(YDDIMV, KIDIA, KFDIA, KLEV, IW ,&
     &   ZREF , ZREFD , ZTRA , ZTRAD ,&
     &   ZDBT , ZRDND , ZRUP , ZRUPD , ZTDBT ,&
     &   ZFD  , ZFU   , PRMU0)  

    !-- up and down-welling fluxes at levels
    DO JK=1,KLEV+1
      DO JL = KIDIA, KFDIA
        IF (PRMU0(JL) > 0.0_JPRB) THEN
          !-- accumulation of spectral fluxes          
          PBBFU(JL,JK) = PBBFU(JL,JK) + ZINCFLX(JL,IW(JL))*ZFU(JL,JK,IW(JL))
          PBBFD(JL,JK) = PBBFD(JL,JK) + ZINCFLX(JL,IW(JL))*ZFD(JL,JK,IW(JL))
          PBBCU(JL,JK) = PBBCU(JL,JK) + ZINCFLX(JL,IW(JL))*ZCU(JL,JK,IW(JL))
          PBBCD(JL,JK) = PBBCD(JL,JK) + ZINCFLX(JL,IW(JL))*ZCD(JL,JK,IW(JL))

          PBBFDIR(JL,JK)=PBBFDIR(JL,JK)+ZINCFLX(JL,IW(JL))*ZTDBT (JL,JK)
          PBBCDIR(JL,JK)=PBBCDIR(JL,JK)+ZINCFLX(JL,IW(JL))*ZTDBTC(JL,JK)
        ENDIF
      ENDDO
    ENDDO
    DO JL = KIDIA, KFDIA
      IF (PRMU0(JL) > 0.0_JPRB) THEN
        PSUDU(JL)=PSUDU(JL)     + ZINCFLX(JL,IW(JL))*ZTDBT (JL,KLEV+1)
        PSUDUC(JL)= PSUDUC(JL)  + ZINCFLX(JL,IW(JL))*ZTDBTC(JL,KLEV+1)
      ENDIF
    ENDDO

  ENDDO
  !-- end loop on JG

ENDDO
!-- end loop on JB

!     ------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SRTM_SPCVRT',1,ZHOOK_HANDLE)
END SUBROUTINE SRTM_SPCVRT
