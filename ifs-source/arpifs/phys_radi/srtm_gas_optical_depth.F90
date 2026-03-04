#ifdef RS6K
@PROCESS HOT(NOVECTOR) NOSTRICT
#endif
SUBROUTINE SRTM_GAS_OPTICAL_DEPTH &
 & ( KIDIA   , KFDIA   , KLEV    , PONEMINUS, &
 &   PRMU0, &
 &   KLAYTROP,&
 &   PCOLCH4  , PCOLCO2 , PCOLH2O , PCOLMOL  , PCOLO2 , PCOLO3 ,&
 &   PFORFAC , PFORFRAC , KINDFOR , PSELFFAC, PSELFFRAC, KINDSELF ,&
 &   PFAC00  , PFAC01   , PFAC10  , PFAC11 ,&
 &   KJP     , KJT      , KJT1 ,&
 !-- output arrays 
 &   POD, PSSA, PINCSOL)


!**** *SRTM_GAS_OPTICAL_DEPTH* - SPECTRAL LOOP TO COMPUTE THE SHORTWAVE RADIATION FLUXES.

!     PURPOSE.
!     --------

!          COMPUTE THE GAS OPTICAL DEPTH AT EACH SHORTWAVE G POINT

!**   INTERFACE.
!     ----------

!          *SRTM_GAS_OPTICAL_DEPTH* IS CALLED FROM THE NEW RADIATION SCHEME

!        IMPLICIT ARGUMENTS :
!        --------------------

!     ==== INPUTS ===
!     ==== OUTPUTS ===

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------

!        SEE RADIATION'S PART OF THE ECMWF RESEARCH DEPARTMENT
!        DOCUMENTATION
!     AUTHOR.
!     -------
!        ADAPTED FROM SRTM_SPCVRT_MCICA (BY JEAN-JACQUES MORCRETTE) BY
!        ROBIN HOGAN
!
!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 2015-07-16

!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE PARSRTM  , ONLY : JPB1, JPB2
USE YOESRTM  , ONLY : JPGPT
USE YOESRTWN , ONLY : NGC
USE YOMDIMV  , ONLY : TDIMV
!USE YOERAD   , ONLY : NSW

IMPLICIT NONE

!     ------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA, KFDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV
REAL(KIND=JPRB)   ,INTENT(IN)    :: PONEMINUS(KIDIA:KFDIA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRMU0(KIDIA:KFDIA)
INTEGER(KIND=JPIM),INTENT(IN)    :: KLAYTROP(KIDIA:KFDIA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCOLCH4(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCOLCO2(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCOLH2O(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCOLMOL(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCOLO2(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCOLO3(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFORFAC(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFORFRAC(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM),INTENT(IN)    :: KINDFOR(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSELFFAC(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSELFFRAC(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM),INTENT(IN)    :: KINDSELF(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFAC00(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFAC01(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFAC10(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFAC11(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM),INTENT(IN)    :: KJP(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM),INTENT(IN)    :: KJT(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM),INTENT(IN)    :: KJT1(KIDIA:KFDIA,KLEV)

REAL(KIND=JPRB)   ,INTENT(OUT)   :: POD(KIDIA:KFDIA,KLEV,JPGPT) ! Optical depth
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSSA(KIDIA:KFDIA,KLEV,JPGPT) ! Single scattering albedo
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PINCSOL(KIDIA:KFDIA,JPGPT) ! Incoming solar flux

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: IB1, IB2, IBM, IGT, IW(KIDIA:KFDIA), JB, JG, JK, JL, IC, ICOUNT

INTEGER(KIND=JPIM) :: INDEX(KFDIA-KIDIA+1)

!-- Output of SRTM_TAUMOLn routines
REAL(KIND=JPRB) :: ZTAUG(KIDIA:KFDIA,KLEV,16) ! Absorption optical depth
REAL(KIND=JPRB) :: ZTAUR(KIDIA:KFDIA,KLEV,16) ! Rayleigh optical depth
REAL(KIND=JPRB) :: ZSFLXZEN(KIDIA:KFDIA,16) ! Incoming solar flux

TYPE(TDIMV) :: YDDIMV

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

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SRTM_GAS_OPTICAL_DEPTH',0,ZHOOK_HANDLE)

IB1=JPB1
IB2=JPB2

IC=0
DO JL = KIDIA, KFDIA
  IF (PRMU0(JL) > 0.0_JPRB) THEN
    IC=IC+1
    INDEX(IC)=JL
    IW(JL)=0
  ENDIF
ENDDO
ICOUNT=IC
IF (ICOUNT /=0) THEN

YDDIMV%NFLEVG = KLEV

JB=IB1-1
DO JB = IB1, IB2
  DO IC=1,ICOUNT
    JL=INDEX(IC)
    IBM = JB-15
    IGT = NGC(IBM)
  ENDDO

  !-- for each band, computes the gaseous and Rayleigh optical thickness 
  !  for all g-points within the band

  IF (JB == 16) THEN
    CALL SRTM_TAUMOL16 &
     & ( YDDIMV  , KIDIA   , KFDIA    , KLEV    ,&
     &   PFAC00  , PFAC01   , PFAC10   , PFAC11   ,&
     &   KJP     , KJT      , KJT1     , PONEMINUS,&
     &   PCOLH2O , PCOLCH4  , PCOLMOL  ,&
     &   KLAYTROP, PSELFFAC , PSELFFRAC, KINDSELF, PFORFAC  , PFORFRAC, KINDFOR ,&
     &   ZSFLXZEN, ZTAUG    , ZTAUR    , PRMU0     &
     & )  

  ELSEIF (JB == 17) THEN
    CALL SRTM_TAUMOL17 &
     & ( YDDIMV  , KIDIA   , KFDIA   , KLEV    ,&
     &   PFAC00  , PFAC01  , PFAC10   , PFAC11 ,&
     &   KJP     , KJT     , KJT1     , PONEMINUS ,&
     &   PCOLH2O , PCOLCO2 , PCOLMOL  ,&
     &   KLAYTROP, PSELFFAC, PSELFFRAC, KINDSELF  , PFORFAC, PFORFRAC, KINDFOR ,&
     &   ZSFLXZEN, ZTAUG   , ZTAUR    , PRMU0     &
     & )  

  ELSEIF (JB == 18) THEN
    CALL SRTM_TAUMOL18 &
     & ( YDDIMV  , KIDIA   , KFDIA   , KLEV    ,&
     &   PFAC00  , PFAC01  , PFAC10   , PFAC11 ,&
     &   KJP     , KJT     , KJT1     , PONEMINUS ,&
     &   PCOLH2O , PCOLCH4 , PCOLMOL  ,&
     &   KLAYTROP, PSELFFAC, PSELFFRAC, KINDSELF  , PFORFAC, PFORFRAC, KINDFOR ,&
     &   ZSFLXZEN, ZTAUG   , ZTAUR    , PRMU0     &
     & )  

  ELSEIF (JB == 19) THEN
    CALL SRTM_TAUMOL19 &
     & ( YDDIMV  , KIDIA   , KFDIA   , KLEV    ,&
     &   PFAC00  , PFAC01  , PFAC10   , PFAC11 ,&
     &   KJP     , KJT     , KJT1     , PONEMINUS ,&
     &   PCOLH2O , PCOLCO2 , PCOLMOL  ,&
     &   KLAYTROP, PSELFFAC, PSELFFRAC, KINDSELF  , PFORFAC, PFORFRAC, KINDFOR ,&
     &   ZSFLXZEN, ZTAUG   , ZTAUR    , PRMU0     &
     & )  

  ELSEIF (JB == 20) THEN
    CALL SRTM_TAUMOL20 &
     & ( YDDIMV  , KIDIA   , KFDIA   , KLEV    ,&
     &   PFAC00  , PFAC01  , PFAC10   , PFAC11 ,&
     &   KJP     , KJT     , KJT1     ,&
     &   PCOLH2O , PCOLCH4 , PCOLMOL  ,&
     &   KLAYTROP, PSELFFAC, PSELFFRAC, KINDSELF  , PFORFAC, PFORFRAC, KINDFOR ,&
     &   ZSFLXZEN, ZTAUG   , ZTAUR    , PRMU0     &
     & )  

  ELSEIF (JB == 21) THEN
    CALL SRTM_TAUMOL21 &
     & ( YDDIMV  , KIDIA   , KFDIA   , KLEV    ,&
     &   PFAC00  , PFAC01  , PFAC10   , PFAC11 ,&
     &   KJP     , KJT     , KJT1     , PONEMINUS ,&
     &   PCOLH2O , PCOLCO2 , PCOLMOL  ,&
     &   KLAYTROP, PSELFFAC, PSELFFRAC, KINDSELF  , PFORFAC, PFORFRAC, KINDFOR ,&
     &   ZSFLXZEN, ZTAUG   , ZTAUR    , PRMU0     &
     & )  

  ELSEIF (JB == 22) THEN
    CALL SRTM_TAUMOL22 &
     & ( YDDIMV  , KIDIA   , KFDIA   , KLEV    ,&
     &   PFAC00  , PFAC01  , PFAC10   , PFAC11 ,&
     &   KJP     , KJT     , KJT1     , PONEMINUS ,&
     &   PCOLH2O , PCOLMOL , PCOLO2   ,&
     &   KLAYTROP, PSELFFAC, PSELFFRAC, KINDSELF  , PFORFAC, PFORFRAC, KINDFOR ,&
     &   ZSFLXZEN, ZTAUG   , ZTAUR    , PRMU0     &
     & )  

  ELSEIF (JB == 23) THEN
    CALL SRTM_TAUMOL23 &
     & ( YDDIMV  , KIDIA   , KFDIA   , KLEV    ,&
     &   PFAC00  , PFAC01  , PFAC10   , PFAC11 ,&
     &   KJP     , KJT     , KJT1     ,&
     &   PCOLH2O , PCOLMOL ,&
     &   KLAYTROP, PSELFFAC, PSELFFRAC, KINDSELF  , PFORFAC, PFORFRAC, KINDFOR ,&
     &   ZSFLXZEN, ZTAUG   , ZTAUR    , PRMU0     &
     & )  

  ELSEIF (JB == 24) THEN
    CALL SRTM_TAUMOL24 &
     & ( YDDIMV  , KIDIA   , KFDIA   , KLEV    ,&
     &   PFAC00  , PFAC01  , PFAC10   , PFAC11 ,&
     &   KJP     , KJT     , KJT1     , PONEMINUS ,&
     &   PCOLH2O , PCOLMOL , PCOLO2   , PCOLO3 ,&
     &   KLAYTROP, PSELFFAC, PSELFFRAC, KINDSELF  , PFORFAC, PFORFRAC, KINDFOR ,&
     &   ZSFLXZEN, ZTAUG   , ZTAUR    , PRMU0     &
     & )  

  ELSEIF (JB == 25) THEN
    !--- visible 16000-22650 cm-1   0.4415 - 0.6250 um
    CALL SRTM_TAUMOL25 &
     & ( YDDIMV  , KIDIA    , KFDIA   , KLEV     ,&
     &   PFAC00   , PFAC01  , PFAC10 , PFAC11 ,&
     &   KJP      , KJT     , KJT1   ,&
     &   PCOLH2O  , PCOLMOL , PCOLO3 ,&
     &   KLAYTROP ,&
     &   ZSFLXZEN, ZTAUG   , ZTAUR   , PRMU0     &
     & )  

  ELSEIF (JB == 26) THEN
    !--- UV-A 22650-29000 cm-1   0.3448 - 0.4415 um
    CALL SRTM_TAUMOL26 &
     & ( YDDIMV  , KIDIA   , KFDIA   , KLEV    ,&
     &   PCOLMOL ,KLAYTROP,&
     &   ZSFLXZEN, ZTAUG   , ZTAUR    , PRMU0     &
     & )  

  ELSEIF (JB == 27) THEN
    !--- UV-B 29000-38000 cm-1   0.2632 - 0.3448 um
    CALL SRTM_TAUMOL27 &
     & ( YDDIMV  , KIDIA   , KFDIA   , KLEV    ,&
     &   PFAC00  , PFAC01  , PFAC10   , PFAC11 ,&
     &   KJP     , KJT     , KJT1     ,&
     &   PCOLMOL , PCOLO3 ,&
     &   KLAYTROP ,&
     &   ZSFLXZEN, ZTAUG   , ZTAUR    , PRMU0     &
     & )  

  ELSEIF (JB == 28) THEN
    !--- UV-C 38000-50000 cm-1   0.2000 - 0.2632 um
    CALL SRTM_TAUMOL28 &
     & ( YDDIMV  , KIDIA   , KFDIA   , KLEV    ,&
     &   PFAC00  , PFAC01  , PFAC10 , PFAC11 ,&
     &   KJP     , KJT     , KJT1   , PONEMINUS ,&
     &   PCOLMOL , PCOLO2  , PCOLO3 ,&
     &   KLAYTROP ,&
     &   ZSFLXZEN, ZTAUG   , ZTAUR  , PRMU0     &
     & )  

  ELSEIF (JB == 29) THEN
    CALL SRTM_TAUMOL29 &
     & ( YDDIMV  , KIDIA    , KFDIA   , KLEV     ,&
     &   PFAC00   , PFAC01  , PFAC10   , PFAC11 ,&
     &   KJP      , KJT     , KJT1     ,&
     &   PCOLH2O  , PCOLCO2 , PCOLMOL  ,&
     &   KLAYTROP , PSELFFAC, PSELFFRAC, KINDSELF  , PFORFAC, PFORFRAC, KINDFOR ,&
     &   ZSFLXZEN , ZTAUG   , ZTAUR    , PRMU0     &
     & )  

  ENDIF
   
  DO JG=1,IGT
    DO IC=1,ICOUNT
      JL=INDEX(IC)
      IW(JL)=IW(JL)+1

      ! Incoming solar flux into plane perp to incoming radiation
      PINCSOL(JL,IW(JL)) = ZSFLXZEN(JL,JG)
    ENDDO

    DO JK=1,KLEV
      DO IC=1,ICOUNT
        JL=INDEX(IC)
        POD (JL,JK,IW(JL)) = ZTAUR(JL,JK,JG) + ZTAUG(JL,JK,JG)
        PSSA(JL,JK,IW(JL)) = ZTAUR(JL,JK,JG) / POD(JL,JK,IW(JL))
      ENDDO
    ENDDO

  ENDDO   !-- end loop on JG (g point)

ENDDO     !-- end loop on JB (band)

ENDIF
!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SRTM_GAS_OPTICAL_DEPTH',1,ZHOOK_HANDLE)
END SUBROUTINE SRTM_GAS_OPTICAL_DEPTH
