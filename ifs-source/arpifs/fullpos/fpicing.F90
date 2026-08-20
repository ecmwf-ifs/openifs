! (C) Copyright 1989- Meteo-France.

SUBROUTINE FPICING(KST,KEND,KPROMA,KLEV,PRESF,PT,PRHF,PICE)

! --------------------------------------------------------------
! **** *FPICING* COMPUTE ICING INDEX FOR AERONAUTIC
! --------------------------------------------------------------
! SUBJECT:
!    ROUTINE COMPUTING ICING INDEX FOR AERONAUTIC FROM T AND RELATIVE
!    HUMIDITY

! INTERFACE:
!    *CALL* *FPICING*

! --------------------------------------------------------------
! -   INPUT ARGUMENTS
!     ---------------

! - DIMENSIONING

! KST      : FIRST INDEX OF LOOPS
! KEND     : LAST INDEX OF LOOPS
! KPROMA   : DEPTH OF THE VECTORIZATION ARRAYS
! KLEV     : END OF VERTICAL LOOP AND VERTICAL DIMENSION

! PT       : TEMPERATURE (K)
!! PQV      : WATER VAPOUR SPECIFIC HUMIDITY (NO DIM) ????
! PRHF     : RELATIVE HUMIDITY

! --------------------------------------------------------------
! -   OUTPUT ARGUMENTS
!     ---------------
! - VARIABLES
! PICE   : ICING index

! --------------------------------------------------------------
! -   IMPLICITE ARGUMENTS
!     -------------------
! --------------------------------------------------------------
! EXTERNALS:

! METHOD:      
!  Compute index according to 3 classes of [T,HU]
!  Index is bounded in [0,10]      
      
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KEND     
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRESF(KPROMA,KLEV)   
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KPROMA,KLEV)   
!REAL(KIND=JPRB)   ,INTENT(IN)    :: PQV(KPROMA,KLEV)   
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRHF(KPROMA,KLEV)   
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PICE(KPROMA,KLEV)   

INTEGER(KIND=JPIM) :: JROF, JLEV
REAL(KIND=JPRB)    :: ZRHF(KPROMA,KLEV), ZES(KPROMA,KLEV) , ZT(KPROMA,KLEV) ,  ZRHFP(KPROMA,KLEV), ZICE(KPROMA,KLEV)!remove ZES
REAL(KIND=JPHOOK)  :: ZHOOK_HANDLE
REAL(KIND=JPRB)    :: ZEPS

! ==== PARAMETERS FOR ICING INDEX  ===================
! Definition of classes
REAL(KIND=JPRB), PARAMETER  :: PPT1 = 0.0_JPRB 
REAL(KIND=JPRB), PARAMETER  :: PPT2 = -5.0_JPRB
REAL(KIND=JPRB), PARAMETER  :: PPT3 = -15.0_JPRB
REAL(KIND=JPRB), PARAMETER  :: PPT4 = -21.0_JPRB
REAL(KIND=JPRB), PARAMETER  :: PPHUMIN = 80._JPRB
! Coefficients for combinaison 
REAL(KIND=JPRB), PARAMETER  :: PPALPHA1 = -2.0/5.0_JPRB
REAL(KIND=JPRB), PARAMETER  :: PPALPHA2 = 1.0/5.0_JPRB
REAL(KIND=JPRB), PARAMETER  :: PPALPHA3 = 4.0/3.0_JPRB
REAL(KIND=JPRB), PARAMETER  :: PPBETA1 = 2.0/5.0_JPRB
REAL(KIND=JPRB), PARAMETER  :: PPBETA2 = 2.0/5.0_JPRB
REAL(KIND=JPRB), PARAMETER  :: PPBETA3 = 4.0/5.0_JPRB
REAL(KIND=JPRB), PARAMETER  :: PPGAMMA1 = -32.0_JPRB
REAL(KIND=JPRB), PARAMETER  :: PPGAMMA2 = -29.0_JPRB
REAL(KIND=JPRB), PARAMETER  :: PPGAMMA3 = -52.0_JPRB
      

!     ------------------------------------------------------------------
!#include "gprh.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('FPICING',0,ZHOOK_HANDLE)


!     ------------------------------------------------------------------
ZEPS = 1.E-6_JPRB

! Icing index
DO JROF=KST,KEND
  DO JLEV=1,KLEV
    ZT(JROF,JLEV)= PT(JROF,JLEV) - 273.15_JPRB ! T in degres
    ZRHFP(JROF,JLEV)= PRHF(JROF,JLEV) * 100._JPRB ! Hu in percent
    ZICE(JROF,JLEV)=&
        & MAX(0._JPRB, SIGN(1._JPRB, PPT1-ZT(JROF,JLEV)))*MAX(0._JPRB,SIGN(1._JPRB, ZT(JROF,JLEV)-PPT2)) *&
        &  MAX(0._JPRB, SIGN(1._JPRB, ZRHFP(JROF,JLEV)-PPHUMIN)) * &
        &  (PPALPHA1*ZT(JROF,JLEV)+PPBETA1*ZRHFP(JROF,JLEV)+PPGAMMA1) +&
        & MAX(0._JPRB, SIGN(1._JPRB, PPT2-ZT(JROF,JLEV)-ZEPS))*MAX(0._JPRB, SIGN(1._JPRB, ZT(JROF,JLEV)-PPT3)) * &
        &  MAX(0._JPRB, SIGN(1._JPRB, ZRHFP(JROF,JLEV)-PPHUMIN)) * &
        &  (PPALPHA2*ZT(JROF,JLEV)+PPBETA2*ZRHFP(JROF,JLEV)+PPGAMMA2) +&
        & MAX(0._JPRB, SIGN(1._JPRB,PPT3-ZT(JROF,JLEV)-ZEPS))*MAX(0._JPRB, SIGN(1._JPRB,ZT(JROF,JLEV)-PPT4)) * &
        &  MAX(0._JPRB, SIGN(1._JPRB, ZRHFP(JROF,JLEV)-PPHUMIN)) * &
        &  (PPALPHA3*ZT(JROF,JLEV)+PPBETA3*ZRHFP(JROF,JLEV)+PPGAMMA3)

        PICE(JROF,JLEV)=MAX(0._JPRB, MIN(10._JPRB ,ZICE(JROF,JLEV))) * 10._JPRB 
  ENDDO
ENDDO     

IF (LHOOK) CALL DR_HOOK('FPICING',1,ZHOOK_HANDLE)

END SUBROUTINE FPICING
