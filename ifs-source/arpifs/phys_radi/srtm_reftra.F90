!!!#ifdef RS6K
!!!@PROCESS HOT NOSTRICT
!!!#endif
SUBROUTINE SRTM_REFTRA&
 & (YDDIMV, KIDIA , KFDIA, KLEV  , KMODTS, &
 & LDRTCHK,&
 & PGG   , PRMUZ, PTAU , PW,&
 & PREF  , PREFD, PTRA , PTRAD&
 & )  

!**** *SRTM_REFTRA* - REFLECTIVITY AND TRANSMISSIVITY

!     PURPOSE.
!     --------
!           COMPUTES THE REFLECTIVITY AND TRANSMISSIVITY OF A CLEAR OR 
!     CLOUDY LAYER USING A CHOICE OF VARIOUS APPROXIMATIONS.

!**   INTERFACE.
!     ----------
!          *SRTM_REFTRA* IS CALLED BY *SRTM_SPCVRT*

!        EXPLICIT ARGUMENTS :
!        --------------------
! INPUTS
! ------ 
!      KMODTS  = 1 EDDINGTON (JOSEPH ET AL., 1976)
!              = 2 PIFM (ZDUNKOWSKI ET AL., 1980)
!              = 3 DISCRETE ORDINATES (LIOU, 1973)
!      LDRTCHK = .T. IF CLOUDY
!              = .F. IF CLEAR-SKY
!      PGG     = ASSYMETRY FACTOR
!      PRMUZ   = COSINE SOLAR ZENITH ANGLE
!      PTAU    = OPTICAL THICKNESS
!      PW      = SINGLE SCATTERING ALBEDO

! OUTPUTS
! -------
!      PREF    : COLLIMATED BEAM REFLECTIVITY
!      PREFD   : DIFFUSE BEAM REFLECTIVITY 
!      PTRA    : COLLIMATED BEAM TRANSMISSIVITY
!      PTRAD   : DIFFUSE BEAM TRANSMISSIVITY

!     METHOD.
!     -------
!          STANDARD DELTA-EDDINGTON, P.I.F.M., OR D.O.M. LAYER CALCULATIONS.

!     EXTERNALS.
!     ----------
!          NONE

!     REFERENCE.
!     ----------

!     AUTHOR.
!     -------
!        JEAN-JACQUES MORCRETTE  *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 03-02-27
!        M.Hamrud   01-Oct-2003      CY28 Cleaning
!        Mike Iacono, AER, Mar 2004: bug fix 
!        D.Salmond  31-Oct-2007 Vector version in the style of RRTM from Meteo France & NEC
!        JJMorcrette/MJIacono 20080724 Look-up table replacing exponential
!        ABozzo 19Mar2014 fixed bug in the comptutation of ZDEND (removed the ZBETA variable)
!        F. Vana  05-Mar-2015  Support for single precision
!        R. Hogan 25-Apr-2015  Remove ZEP1 and ZEP2 to avoid problems in single precision
!        R. Hogan 07-Jun-2016  Increased ZWCRIT from 0.9995 to 0.9999995

!     ------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

USE YOMDIMV  , ONLY : TDIMV
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOERDU   , ONLY : REPLOG
USE YOESRTAB , ONLY : BPADE, TRANS, RODLOW, RTBLINT

IMPLICIT NONE

TYPE(TDIMV)       ,INTENT(IN)      :: YDDIMV
INTEGER(KIND=JPIM),INTENT(IN)      :: KIDIA, KFDIA
INTEGER(KIND=JPIM),INTENT(IN)      :: KLEV 
INTEGER(KIND=JPIM),INTENT(OUT)     :: KMODTS 
LOGICAL           ,INTENT(IN)      :: LDRTCHK(KIDIA:KFDIA,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)      :: PGG(KIDIA:KFDIA,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)      :: PRMUZ(KIDIA:KFDIA) 
REAL(KIND=JPRB)   ,INTENT(IN)      :: PTAU(KIDIA:KFDIA,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)      :: PW(KIDIA:KFDIA,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT)   :: PREF(KIDIA:KFDIA,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT)   :: PREFD(KIDIA:KFDIA,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT)   :: PTRA(KIDIA:KFDIA,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT)   :: PTRAD(KIDIA:KFDIA,YDDIMV%NFLEVG) 
!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JK, JL, IC, INDEX(KIDIA:KFDIA), ICOUNT, ITIND

REAL(KIND=JPRB) :: ZA, ZA1, ZA2
REAL(KIND=JPRB) :: ZDEND, ZDEN
REAL(KIND=JPRB) :: ZE1, ZE2, ZEM1, ZEM2, ZEMM
REAL(KIND=JPRB) :: ZG, ZG3, ZGAMMA1, ZGAMMA2, ZGAMMA3, ZGAMMA4, ZGT
REAL(KIND=JPRB) :: ZR1, ZR2, ZR3, ZR4, ZR5, ZRK, ZRK2, ZRKG, ZRM1, ZRP, ZRP1, ZRPP
REAL(KIND=JPRB) :: ZSR3, ZT1, ZT2, ZT3, ZTO1, ZTBLIND
REAL(KIND=JPRB) :: ZW, ZWCRIT,ZTEMP,ZSQRT_REPLOG,ZZZ
!!!REAL(KIND=JPRB) :: ZWO, ZEXP500,ZEXPM500,ZEXP500R

! The following has been reduced from 1.e-8 as it is applied to a
! quantity that has now been multiplied by ZEM1
! (FV: It doesn't seem to make any harm in single precision either.)
REAL(KIND=JPRB), PARAMETER :: ZEPS = 1.E-16_JPRB

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SRTM_REFTRA',0,ZHOOK_HANDLE)

!!!ZEXP500=EXP(500.0_JPRB)
!!!ZEXP500R=1.0_JPRB/ZEXP500
!!!ZEXPM500=EXP(-500.0_JPRB)

ZSR3=SQRT(3._JPRB)
ZWCRIT=0.9999995_JPRB
KMODTS=2

ZSQRT_REPLOG = SQRT(REPLOG)

IC=0
DO JL = KIDIA, KFDIA
  IF (PRMUZ(JL) > 0.0_JPRB) THEN
    IC=IC+1
    INDEX(IC)=JL
  ENDIF
ENDDO
ICOUNT=IC
IF (ICOUNT > 0) THEN

DO JK=1,KLEV
  DO IC=1,ICOUNT
    JL=INDEX(IC)
    IF (.NOT.LDRTCHK(JL,JK)) THEN
      PREF(JL,JK) =0.0_JPRB
      PTRA(JL,JK) =1.0_JPRB
      PREFD(JL,JK)=0.0_JPRB
      PTRAD(JL,JK)=1.0_JPRB
    ELSE
      ZTO1=PTAU(JL,JK)
      ZW  =PW(JL,JK)
      ZG  =PGG(JL,JK)  

      !-- GENERAL TWO-STREAM EXPRESSIONS

      ZG3= 3._JPRB * ZG
      IF (KMODTS == 1) THEN
        ZGAMMA1= (7._JPRB - ZW * (4._JPRB + ZG3)) * 0.25_JPRB
        ZGAMMA2=-(1._JPRB - ZW * (4._JPRB - ZG3)) * 0.25_JPRB
        ZGAMMA3= (2._JPRB - ZG3 * PRMUZ(JL) ) * 0.25_JPRB
      ELSEIF (KMODTS == 2) THEN  
        ZGAMMA1= (8._JPRB - ZW * (5._JPRB + ZG3)) * 0.25_JPRB
        ZGAMMA2=  3._JPRB *(ZW * (1._JPRB - ZG )) * 0.25_JPRB
        ZGAMMA3= (2._JPRB - ZG3 * PRMUZ(JL) ) * 0.25_JPRB
      ELSEIF (KMODTS == 3) THEN  
        ZGAMMA1= ZSR3 * (2._JPRB - ZW * (1._JPRB + ZG)) * 0.5_JPRB
        ZGAMMA2= ZSR3 * ZW * (1._JPRB - ZG ) * 0.5_JPRB
        ZGAMMA3= (1._JPRB - ZSR3 * ZG * PRMUZ(JL) ) * 0.5_JPRB
      ENDIF
      ZGAMMA4= 1._JPRB - ZGAMMA3

      !-- RECOMPUTE ORIGINAL S.S.A. TO TEST FOR CONSERVATIVE SOLUTION
      !   ZTEMP=(1._JPRB - ZG)**2
      !   ZWO= ZW*ZTEMP/ (ZTEMP - (1._JPRB - ZW)*(ZG **2))

!       ZWO= ZW / (1._JPRB - (1._JPRB - ZW) * (ZG / (1._JPRB - ZG))**2)
!       IF (ZWO >= ZWCRIT) THEN

      ZZZ=(1._JPRB - ZG)**2
      IF (ZW*ZZZ >= ZWCRIT*(ZZZ - (1._JPRB - ZW)*(ZG **2))) THEN
        !!-- conservative scattering

        ZA  = ZGAMMA1 * PRMUZ(JL) 
        ZA1 = ZA - ZGAMMA3
        ZGT = ZGAMMA1 * ZTO1

        !-- Homogeneous reflectance and transmittance

        ! collimated beam

!!         ZE1 = MIN ( ZTO1 / PRMUZ(JL) , 500._JPRB)
!!         ZE2 = EXP ( - ZE1 )
!        ZZZ= ZTO1/PRMUZ(JL)
!        IF(ZZZ <= 500._JPRB) THEN
!          ZE2 = EXP ( - ZZZ )
!        ELSE
!          ZE2 = ZEXPM500
!        ENDIF

!-- Use exponential look-up table for transmittance, or expansion of
!   exponential for low optical thickness
        ZE1 = MIN ( ZTO1 / PRMUZ(JL) , 500._JPRB)
        IF (ZE1 <= RODLOW) THEN
          ZE2 = 1._JPRB - ZE1 + 0.5_JPRB*ZE1*ZE1
        ELSE
          ZTBLIND = ZE1 / (BPADE + ZE1)
          ITIND = RTBLINT * ZTBLIND + 0.5_JPRB
          ZE2 = TRANS(ITIND)
        ENDIF
!---

        ZTEMP=1.0_JPRB/(1._JPRB + ZGT)
        PREF(JL,JK) = (ZGT - ZA1 * (1._JPRB - ZE2)) *ZTEMP
        PTRA(JL,JK) = 1._JPRB - PREF(JL,JK)

        ! isotropic incidence

        PREFD(JL,JK) = ZGT *ZTEMP
        PTRAD(JL,JK) = 1._JPRB - PREFD(JL,JK)        

      ELSE

        !-- non-conservative scattering

        ZA1 = ZGAMMA1 * ZGAMMA4 + ZGAMMA2 * ZGAMMA3
        ZA2 = ZGAMMA1 * ZGAMMA3 + ZGAMMA2 * ZGAMMA4
        !      ZRK = SQRT ( ZGAMMA1**2 - ZGAMMA2**2)
!         ZRK = SQRT ( MAX ( REPLOG, ZGAMMA1**2 - ZGAMMA2**2) )
        ZZZ = ZGAMMA1**2 - ZGAMMA2**2
        IF (ZZZ >= REPLOG ) THEN
          ZRK = SQRT ( ZZZ )
        ELSE
          ZRK = ZSQRT_REPLOG
        ENDIF

        ZRP = ZRK * PRMUZ(JL)               
        ZRP1 = 1._JPRB + ZRP
        ZRM1 = 1._JPRB - ZRP
        ZRK2 = 2._JPRB * ZRK
        ZRPP = 1._JPRB - ZRP*ZRP
        ZRKG = ZRK + ZGAMMA1
        ZR1  = ZRM1 * (ZA2 + ZRK * ZGAMMA3)
        ZR2  = ZRP1 * (ZA2 - ZRK * ZGAMMA3)
        ZR3  = ZRK2 * (ZGAMMA3 - ZA2 * PRMUZ(JL) )
        ZR4  = ZRPP * ZRKG
        ZR5  = ZRPP * (ZRK - ZGAMMA1)
        ZT1  = ZRP1 * (ZA1 + ZRK * ZGAMMA4)
        ZT2  = ZRM1 * (ZA1 - ZRK * ZGAMMA4)
        ZT3  = ZRK2 * (ZGAMMA4 + ZA1 * PRMUZ(JL) )

        !-- Homogeneous reflectance and transmittance

!!     ZE1 = MIN ( ZRK * ZTO1, 500._JPRB)
!!     ZE2 = MIN ( ZTO1 / PRMUZ , 500._JPRB)
!!
!!     ZEP1 = EXP( ZE1 )
!!     ZEM1 = EXP(-ZE1 )
!!     ZEM1=1.0_JPRB/ZEP1
!!
!!     ZEP2 = EXP( ZE2 )
!!     ZEM2 = EXP(-ZE2 )
!!     ZEM2=1.0_JPRB/ZEP2

!        IF(ZRK * ZTO1 > 500._JPRB)THEN
!          ZEP1=ZEXP500
!          ZEM1=ZEXP500R
!        ELSE
!          ZEP1=EXP(ZRK * ZTO1)
!          ZEM1=1.0_JPRB/ZEP1
!        ENDIF
!        IF(ZTO1 > 500._JPRB*PRMUZ(JL))THEN
!          ZEP2=ZEXP500
!          ZEM2=ZEXP500R
!        ELSE
!          ZEP2=EXP(ZTO1 / PRMUZ(JL))
!          ZEM2=1.0_JPRB/ZEP2
!        ENDIF

!-- Use exponential look-up table for transmittance, or expansion of
!   exponential for low optical thickness
        ZE1 = MIN ( ZRK * ZTO1 , 500._JPRB)
        ZE2 = MIN ( ZTO1 / PRMUZ(JL) , 500._JPRB)
        IF (ZE1 <= RODLOW) THEN
          ZEM1 = 1._JPRB - ZE1 + 0.5_JPRB*ZE1*ZE1
        ELSE
          ZTBLIND = ZE1 / (BPADE + ZE1)
          ITIND = RTBLINT * ZTBLIND + 0.5_JPRB
          ZEM1 = TRANS(ITIND)
        ENDIF

        IF (ZE2 <= RODLOW) THEN
          ZEM2 = 1._JPRB - ZE2 + 0.5_JPRB*ZE2*ZE2
        ELSE
          ZTBLIND = ZE2 / (BPADE + ZE2)
          ITIND = RTBLINT * ZTBLIND + 0.5_JPRB
          ZEM2 = TRANS(ITIND)
        ENDIF
!---

        ! collimated beam

        ZEMM = ZEM1*ZEM1
        ! Original code:
        !        ZDENR = ZR4*ZEP1 + ZR5*ZEM1
        !        ZDENT = ZT4*ZEP1 + ZT5*ZEM1
        ! Now multiply by ZEM1 to cancel ZEP1, but note also that both
        ! terms are the same so represent as a single denominator
        ! variable
        ZDEN = ZR4 + ZR5*ZEMM

        IF ((ZDEN > -ZEPS .AND. ZDEN < ZEPS)) THEN 
          PREF(JL,JK) = ZEPS
          PTRA(JL,JK) = ZEM2
        ELSE
          ! Original code:
          !  PREF(JL,JK) = ZW  * (ZR1*ZEP1 - ZR2*ZEM1 - ZR3*ZEM2) / ZDENR
          ! Now multiply through by ZEM1 so that ZEM1*ZEP1=1
          ! (denominator already multiplied by ZEM1)
          PREF(JL,JK) = ZW * (ZR1 - ZR2*ZEMM - ZR3*ZEM1*ZEM2) / ZDEN

          ! Original code:
          !  PTRA(JL,JK) = ZEM2 * (1._JPRB - ZW  * (ZT1*ZEP1 - ZT2*ZEM1 - ZT3*ZEP2) / &
          !  & SIGN(MAX(ABS(ZDENT),ZEPSILON),ZDENT))
          ! Now multiply through ZEM2 so that ZEM2*ZEP2=1 and multiply
          ! through by ZEM1 (denominator already multiplied by
          ! ZEM1). 
          PTRA(JL,JK) = ZEM2 - ZW * (ZT1*ZEM2 - ZT2*ZEMM*ZEM2 - ZT3*ZEM1) &
               & / ZDEN
        ENDIF

        ! Diffuse beam

        ZDEND =  1._JPRB / (ZRKG + (ZRK - ZGAMMA1)*ZEMM)
        !ZDEND = 1._JPRB / ( (1._JPRB - ZBETA*ZEMM ) * ZRKG) !ZBETA was -ZR5/ZR4. Potentially caused 0/0
        PREFD(JL,JK) =  ZGAMMA2 * (1._JPRB - ZEMM) * ZDEND
        PTRAD(JL,JK) =  ZRK2*ZEM1*ZDEND

      ENDIF

    ENDIF

  ENDDO
ENDDO

ENDIF

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SRTM_REFTRA',1,ZHOOK_HANDLE)
END SUBROUTINE SRTM_REFTRA
