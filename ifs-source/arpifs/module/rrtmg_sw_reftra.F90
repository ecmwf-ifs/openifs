!     path:      $Source: /storm/rc1/cvsroot/rc/rrtmg_sw/src/rrtmg_sw_reftra.f90,v $
!     author:    $Author: mike $
!     revision:  $Revision: 1.5 $
!     created:   $Date: 2009/05/22 22:22:22 $

      MODULE RRTMG_SW_REFTRA

!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2002-2009, Atmospheric & Environmental Research, Inc. (AER).  |
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
!  --------------------------------------------------------------------------

! ------- Modules -------

      USE PARKIND1  ,ONLY : JPIM     , JPRB
      USE YOESRTAB , ONLY : BPADE, TRANS, RODLOW, RTBLINT

      IMPLICIT NONE

      CONTAINS

! --------------------------------------------------------------------
      SUBROUTINE REFTRA_SW(YDERAD,KLAYERS, LDRTCHK, PGG, PRMUZ, PTAU, PW, &
                         &  PREF, PREFD, PTRA, PTRAD)
! --------------------------------------------------------------------
  
! Purpose: computes the reflectivity and transmissivity of a clear or 
!   cloudy layer using a choice of various approximations.
!
! Interface:  *rrtmg_sw_reftra* is called by *rrtmg_sw_spcvrt*
!
! Description:
! explicit arguments :
! --------------------
! inputs
! ------ 
!      ldrtchk = .t. for all layers in clear profile
!      ldrtchk = .t. for cloudy layers in cloud profile 
!              = .f. for clear layers in cloud profile
!      pgg     = assymetry factor
!      prmuz   = cosine solar zenith angle
!      ptau    = optical thickness
!      pw      = single scattering albedo
!
! outputs
! -------
!      pref    : collimated beam reflectivity
!      prefd   : diffuse beam reflectivity 
!      ptra    : collimated beam transmissivity
!      ptrad   : diffuse beam transmissivity
!
!
! Method:
! -------
!      standard delta-eddington, p.i.f.m., or d.o.m. layer calculations.
!      kmodts  = 1 eddington (joseph et al., 1976)
!              = 2 pifm (zdunkowski et al., 1980)
!              = 3 discrete ordinates (liou, 1973)
!
!
! Modifications:
! --------------
! Original: J-JMorcrette, ECMWF, Feb 2003
! Revised for F90 reformatting: MJIacono, AER, Jul 2006
! Revised to add exponential lookup table: MJIacono, AER, Aug 2007
!
! ABozzo Dec 2014: small changes (modules) to work 
! as monochromatic UV code in ECMWF IFS
! ABozzo Jun2018: modifiactions for single precision. Removed 
!                 ZEP1 and ZEP2 to avid problems
! ------------------------------------------------------------------

! ------- Declarations ------

! ------- Input -------
USE YOERAD   , ONLY : TERAD
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
  TYPE(TERAD)            ,INTENT(INOUT):: YDERAD
      INTEGER(KIND=JPIM), INTENT(IN) :: KLAYERS

      LOGICAL, INTENT(IN) :: LDRTCHK(:)                         ! Logical flag for reflectivity and
                                                               ! and transmissivity calculation; 
                                                               !   Dimensions: (klayers)

      REAL(KIND=JPRB), INTENT(IN) :: PGG(:)                      ! asymmetry parameter
                                                               !   Dimensions: (klayers)
      REAL(KIND=JPRB), INTENT(IN) :: PTAU(:)                     ! optical depth
                                                               !   Dimensions: (klayers)
      REAL(KIND=JPRB), INTENT(IN) :: PW(:)                       ! single scattering albedo 
                                                               !   Dimensions: (klayers)
      REAL(KIND=JPRB), INTENT(IN) :: PRMUZ                       ! cosine of solar zenith angle

! ------- Output -------

      REAL(KIND=JPRB), INTENT(INOUT) :: PREF(:)                  ! direct beam reflectivity
                                                               !   Dimensions: (klayers+1)
      REAL(KIND=JPRB), INTENT(INOUT) :: PREFD(:)                 ! diffuse beam reflectivity
                                                               !   Dimensions: (klayers+1)
      REAL(KIND=JPRB), INTENT(INOUT) :: PTRA(:)                  ! direct beam transmissivity
                                                               !   Dimensions: (klayers+1)
      REAL(KIND=JPRB), INTENT(INOUT) :: PTRAD(:)                 ! diffuse beam transmissivity
                                                               !   Dimensions: (klayers+1)

! ------- Local -------

      INTEGER(KIND=JPIM) :: JK!, kmodts
      INTEGER(KIND=JPIM) :: ITIND

      REAL(KIND=JPRB) :: ZTBLIND
      REAL(KIND=JPRB) :: ZA, ZA1, ZA2
      REAL(KIND=JPRB) :: ZBETA, ZDEND, ZDEN
      REAL(KIND=JPRB) :: ZE1, ZE2, ZEM1, ZEM2, ZEMM
      REAL(KIND=JPRB) :: ZG, ZG3, ZGAMMA1, ZGAMMA2, ZGAMMA3, ZGAMMA4, ZGT
      REAL(KIND=JPRB) :: ZR1, ZR2, ZR3, ZR4, ZR5
      REAL(KIND=JPRB) :: ZRK, ZRK2, ZRKG, ZRM1, ZRP, ZRP1, ZRPP
      REAL(KIND=JPRB) :: ZSR3, ZT1, ZT2, ZT3, ZT4, ZT5, ZTO1
      REAL(KIND=JPRB) :: ZW, ZWCRIT, ZWO

! The following has been reduced from 1.e-8 as it is applied to a
! quantity that has now been multiplied by ZEM1
! (FV: It doesn't seem to make any harm in single precision either.)
      REAL(KIND=JPRB), PARAMETER :: EPS = 1.E-016_JPRB

!for the exp lookup tables I use the arrays already available in IFS
      REAL(KIND=JPRB) :: ZTBLINT,ZOD_LO
      REAL(KIND=JPRB) , DIMENSION(0:10000) :: ZEXP_TBL
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!     ------------------------------------------------------------------

! Initialize
IF (LHOOK) CALL DR_HOOK('rrtmg_sw_reftra:reftra_sw',0,ZHOOK_HANDLE)
      ZOD_LO=RODLOW
      ZTBLINT=RTBLINT
      ZEXP_TBL=TRANS


      ZSR3=SQRT(3._JPRB)
      ZWCRIT=0.9999995_JPRB
      !kmodts=2 !AB passed through naerad

      DO JK=1, KLAYERS
         IF (.NOT.LDRTCHK(JK)) THEN
            PREF(JK) =0._JPRB
            PTRA(JK) =1._JPRB
            PREFD(JK)=0._JPRB
            PTRAD(JK)=1._JPRB
         ELSE
            ZTO1=PTAU(JK)
            ZW  =PW(JK)
            ZG  =PGG(JK)  

! General two-stream expressions

            ZG3= 3._JPRB * ZG
            IF (YDERAD%KMODTS == 1) THEN
               ZGAMMA1= (7._JPRB - ZW * (4._JPRB + ZG3)) * 0.25_JPRB
               ZGAMMA2=-(1._JPRB - ZW * (4._JPRB - ZG3)) * 0.25_JPRB
               ZGAMMA3= (2._JPRB - ZG3 * PRMUZ ) * 0.25_JPRB
            ELSEIF (YDERAD%KMODTS == 2) THEN  
               ZGAMMA1= (8._JPRB - ZW * (5._JPRB + ZG3)) * 0.25_JPRB
               ZGAMMA2=  3._JPRB *(ZW * (1._JPRB - ZG )) * 0.25_JPRB
               ZGAMMA3= (2._JPRB - ZG3 * PRMUZ ) * 0.25_JPRB
            ELSEIF (YDERAD%KMODTS == 3) THEN  
               ZGAMMA1= ZSR3 * (2._JPRB - ZW * (1._JPRB + ZG)) * 0.5_JPRB
               ZGAMMA2= ZSR3 * ZW * (1._JPRB - ZG ) * 0.5_JPRB
               ZGAMMA3= (1._JPRB - ZSR3 * ZG * PRMUZ ) * 0.5_JPRB
            ENDIF
            ZGAMMA4= 1._JPRB - ZGAMMA3
    
! Recompute original s.s.a. to test for conservative solution

            ZWO= ZW / (1._JPRB - (1._JPRB - ZW) * (ZG / (1._JPRB - ZG))**2)


            IF (ZWO >= ZWCRIT) THEN
! Conservative scattering

               ZA  = ZGAMMA1 * PRMUZ 
               ZA1 = ZA - ZGAMMA3
               ZGT = ZGAMMA1 * ZTO1
        
! Homogeneous reflectance and transmittance,
! collimated beam

               ZE1 = MIN ( ZTO1 / PRMUZ , 500._JPRB)
!               ze2 = exp( -ze1 )

! Use exponential lookup table for transmittance, or expansion of 
! exponential for low tau
               IF (ZE1 <= ZOD_LO) THEN 
                  ZE2 = 1._JPRB - ZE1 + 0.5_JPRB * ZE1 * ZE1
               ELSE
                  ZTBLIND = ZE1 / (BPADE + ZE1)
                  ITIND = ZTBLINT * ZTBLIND + 0.5_JPRB
                  ZE2 = ZEXP_TBL(ITIND)
               ENDIF
!

               PREF(JK) = (ZGT - ZA1 * (1._JPRB - ZE2)) / (1._JPRB + ZGT)
               PTRA(JK) = 1._JPRB - PREF(JK)

! isotropic incidence

               PREFD(JK) = ZGT / (1._JPRB + ZGT)
               PTRAD(JK) = 1._JPRB - PREFD(JK)        

! This is applied for consistency between total (delta-scaled) and direct (unscaled) 
! calculations at very low optical depths (tau < 1.e-4) when the exponential lookup
! table returns a transmittance of 1.0.
               IF (ZE2 == 1.0_JPRB) THEN 
                  PREF(JK) = 0.0_JPRB
                  PTRA(JK) = 1.0_JPRB
                  PREFD(JK) = 0.0_JPRB
                  PTRAD(JK) = 1.0_JPRB
               ENDIF

            ELSE
! Non-conservative scattering

               ZA1 = ZGAMMA1 * ZGAMMA4 + ZGAMMA2 * ZGAMMA3
               ZA2 = ZGAMMA1 * ZGAMMA3 + ZGAMMA2 * ZGAMMA4
               ZRK = SQRT ( ZGAMMA1**2 - ZGAMMA2**2)
               ZRP = ZRK * PRMUZ               
               ZRP1 = 1._JPRB + ZRP
               ZRM1 = 1._JPRB - ZRP
               ZRK2 = 2._JPRB * ZRK
               ZRPP = 1._JPRB - ZRP*ZRP
               ZRKG = ZRK + ZGAMMA1
               ZR1  = ZRM1 * (ZA2 + ZRK * ZGAMMA3)
               ZR2  = ZRP1 * (ZA2 - ZRK * ZGAMMA3)
               ZR3  = ZRK2 * (ZGAMMA3 - ZA2 * PRMUZ )
               ZR4  = ZRPP * ZRKG
               ZR5  = ZRPP * (ZRK - ZGAMMA1)
               ZT1  = ZRP1 * (ZA1 + ZRK * ZGAMMA4)
               ZT2  = ZRM1 * (ZA1 - ZRK * ZGAMMA4)
               ZT3  = ZRK2 * (ZGAMMA4 + ZA1 * PRMUZ )
               ZT4  = ZR4
               ZT5  = ZR5

! mji - reformulated code to avoid potential floating point exceptions
!               zbeta = - zr5 / zr4
               ZBETA = (ZGAMMA1 - ZRK) / ZRKG
!!
        
! Homogeneous reflectance and transmittance

               ZE1 = MIN ( ZRK * ZTO1, 500._JPRB)
               ZE2 = MIN ( ZTO1 / PRMUZ , 500._JPRB)
!
! Original
!              zep1 = exp( ze1 )
!              zem1 = exp(-ze1 )
!              zep2 = exp( ze2 )
!              zem2 = exp(-ze2 )
!
! Revised original, to reduce exponentials
!              zep1 = exp( ze1 )
!              zem1 = 1._JPRB / zep1
!              zep2 = exp( ze2 )
!              zem2 = 1._JPRB / zep2
!

! Use exponential lookup table for transmittance, or expansion of 
! exponential for low tau
               IF (ZE1 <= ZOD_LO) THEN 
                  ZEM1 = 1._JPRB - ZE1 + 0.5_JPRB * ZE1 * ZE1
!                  ZEP1 = 1._JPRB / ZEM1
               ELSE
                  ZTBLIND = ZE1 / (BPADE + ZE1)
                  ITIND = ZTBLINT * ZTBLIND + 0.5_JPRB
                  ZEM1 = ZEXP_TBL(ITIND)
!                  ZEP1 = 1._JPRB / ZEM1
               ENDIF



               IF (ZE2 <= ZOD_LO) THEN 
                  ZEM2 = 1._JPRB - ZE2 + 0.5_JPRB * ZE2 * ZE2
!                  ZEP2 = 1._JPRB / ZEM2
               ELSE
                  ZTBLIND = ZE2 / (BPADE + ZE2)
                  ITIND = ZTBLINT * ZTBLIND + 0.5_JPRB
                  ZEM2 = ZEXP_TBL(ITIND)
!                  ZEP2 = 1._JPRB / ZEM2
               ENDIF



! collimated beam

! mji - reformulated code to avoid potential floating point exceptions
!               zdenr = zr4*zep1 + zr5*zem1
!               pref(jk) = zw * (zr1*zep1 - zr2*zem1 - zr3*zem2) / zdenr
!               zdent = zt4*zep1 + zt5*zem1
!               ptra(jk) = zem2 - zem2 * zw * (zt1*zep1 - zt2*zem1 - zt3*zep2) / zdent

!Abozzo/RJHogan - reformulated to avoid problems with single precision. Original code:
!               ZDENR = ZR4*ZEP1 + ZR5*ZEM1
!               ZDENT = ZT4*ZEP1 + ZT5*ZEM1
! Now multiply by ZEM1 to cancel ZEP1, but note also that both
        ! terms are the same so represent as a single denominator
        ! variable
               ZEMM = ZEM1*ZEM1
               ZDEN = ZR4 + ZR5*ZEMM
               IF (ZDEN >= -EPS .AND. ZDEN <= EPS) THEN
                  PREF(JK) = EPS
                  PTRA(JK) = ZEM2
               ELSE 
                  ! Original code
                  !PREF(JK) = ZW * (ZR1*ZEP1 - ZR2*ZEM1 - ZR3*ZEM2) / ZDENR
                  ! Now multiply through by ZEM1 so that ZEM1*ZEP1=1
                  ! (denominator already multiplied by ZEM1)
                  PREF(JK) = ZW * (ZR1 - ZR2*ZEMM - ZR3*ZEM1*ZEM2) / ZDEN

                  !Original code
                  !PTRA(JK) = ZEM2 - ZEM2 * ZW * (ZT1*ZEP1 - ZT2*ZEM1 - ZT3*ZEP2) / ZDENT
                  ! Now multiply through ZEM2 so that ZEM2*ZEP2=1 and multiply
                  ! through by ZEM1 (denominator already multiplied by
                  ! ZEM1). 
                  PTRA(JK) = ZEM2 - ZW * (ZT1*ZEM2 - ZT2*ZEMM*ZEM2 - ZT3*ZEM1) &
                       & / ZDEN
               ENDIF
!!

! diffuse beam

               ZDEND = 1._JPRB / ( (1._JPRB - ZBETA*ZEMM ) * ZRKG)
               PREFD(JK) =  ZGAMMA2 * (1._JPRB - ZEMM) * ZDEND
               PTRAD(JK) =  ZRK2*ZEM1*ZDEND

            ENDIF

         ENDIF         



      ENDDO    

IF (LHOOK) CALL DR_HOOK('rrtmg_sw_reftra:reftra_sw',1,ZHOOK_HANDLE)

      END SUBROUTINE REFTRA_SW

      END MODULE RRTMG_SW_REFTRA

