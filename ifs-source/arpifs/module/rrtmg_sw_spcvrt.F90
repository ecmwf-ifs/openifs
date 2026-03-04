!     path:      $Source: /storm/rc1/cvsroot/rc/rrtmg_sw/src/rrtmg_sw_spcvrt.f90,v $
!     author:    $Author: mike $
!     revision:  $Revision: 1.5 $
!     created:   $Date: 2009/05/22 22:22:22 $

      MODULE RRTMG_SW_SPCVRT

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
      
      USE PARKIND1  ,ONLY : JPIM     ,JPRB
      USE YOESRTAB , ONLY : BPADE, TRANS, RODLOW, RTBLINT
      USE RRTMG_SW_REFTRA, ONLY: REFTRA_SW
      USE RRTMG_SW_VRTQDR, ONLY: VRTQDR_SW

      IMPLICIT NONE

      CONTAINS

! ---------------------------------------------------------------------------
      SUBROUTINE SPCVRT_SW &
             & (YDERAD,KLAYERS, KSTART, KEND, KCPR, KDELM, &
             & PAVEL, PTAVEL, PZ, PTZ, PTBOUND, PALBD, PALBP, &
             & PCLFR, PCLEAR, PTAUC, PASYC, POMGC, PTAUCORIG, &
             & PTAUR, PTAUA, PASYA, POMGA, PRMU0, PTAUG, PADJFLUX, &
             & PBBFD, PBBFU, PBBCD, PBBCU, &
             & PBBFDDIR, PBBCDDIR)
! ---------------------------------------------------------------------------
!
! Purpose: Contains spectral loop to compute the shortwave radiative fluxes, 
!          using the two-stream method of H. Barker. 
!
! Interface:  *spcvrt_sw* is called from *rrtmg_sw.F90* or rrtmg_sw.1col.F90*
!
! Method:
!    Adapted from two-stream model of H. Barker;
!    Two-stream model options (selected with kmodts in rrtmg_sw_reftra.F90):
!        1: Eddington, 2: PIFM, Zdunkowski et al., 3: discret ordinates
!
! Modifications:
!
! Original: H. Barker
! Revision: Merge with RRTMG_SW: J.-J.Morcrette, ECMWF, Feb 2003
! Revision: Add adjustment for Earth/Sun distance : MJIacono, AER, Oct 2003
! Revision: Bug fix for use of PALBP and PALBD: MJIacono, AER, Nov 2003
! Revision: Bug fix to apply delta scaling to clear sky: AER, Dec 2004
! Revision: Code modified so that delta scaling is not done in cloudy profiles
!           if routine cldprop is used; delta scaling can be applied by swithcing
!           code below if cldprop is not used to get cloud properties. 
!           AER, Jan 2005
! Revision: Uniform formatting for RRTMG: MJIacono, AER, Jul 2006 
! Revision: Use exponential lookup table for transmittance: MJIacono, AER, 
!           Aug 2007 

!!!!!!!!!!!!!!MODIFICATION FOR UV COMPUTATIONS ONLY!!!!!!!!!!!!!!!!!!!!!!!
! ABozzo dec 2014: modified to work at single wavelength in the UV range
!
! ------------------------------------------------------------------

! ------- Declarations ------

! -------- Input -------

USE YOERAD , ONLY : TERAD
  TYPE(TERAD)            ,INTENT(INOUT):: YDERAD
      INTEGER(KIND=JPIM), INTENT(IN) :: KLAYERS
      INTEGER(KIND=JPIM), INTENT(IN) :: KSTART
      INTEGER(KIND=JPIM), INTENT(IN) :: KEND
      INTEGER(KIND=JPIM), INTENT(IN) :: KCPR
      INTEGER(KIND=JPIM), INTENT(IN) :: KDELM   ! delta-m scaling flag
                                              ! [0 = direct and diffuse fluxes are unscaled]
                                              ! [1 = direct and diffuse fluxes are scaled]


      REAL(KIND=JPRB), INTENT(IN) :: PAVEL(:)                    ! layer pressure (hPa, mb) 
                                                               !   Dimensions: (klayers)
      REAL(KIND=JPRB), INTENT(IN) :: PTAVEL(:)                    ! layer temperature (K)
                                                               !   Dimensions: (klayers)
      REAL(KIND=JPRB), INTENT(IN) :: PZ(0:)                      ! level (interface) pressure (hPa, mb)
                                                               !   Dimensions: (0:klayers)
      REAL(KIND=JPRB), INTENT(IN) :: PTZ(0:)                      ! level temperatures (hPa, mb)
                                                               !   Dimensions: (0:klayers)
      REAL(KIND=JPRB), INTENT(IN) :: PTBOUND                      ! surface temperature (K)
                                                               !   Dimensions: (klayers)
      REAL(KIND=JPRB), INTENT(IN) :: PADJFLUX                  ! Earth/Sun distance adjusted solar flux
                                                                  

      REAL(KIND=JPRB), INTENT(IN) :: PALBD                    ! surface albedo (diffuse)
                                                               !   Dimensions: (nbndsw)
      REAL(KIND=JPRB), INTENT(IN) :: PALBP                    ! surface albedo (direct)
                                                               !   Dimensions: (nbndsw)
      REAL(KIND=JPRB), INTENT(IN) :: PRMU0                       ! cosine of solar zenith angle
      REAL(KIND=JPRB), INTENT(IN) :: PCLFR(:)                    ! cloud fraction
                                                               !   Dimensions: (klayers)
      REAL(KIND=JPRB), INTENT(IN) :: PCLEAR                    ! 1-total_cloud_cover
      REAL(KIND=JPRB), INTENT(IN) :: PTAUC(:)                  ! cloud optical depth
                                                               !   Dimensions: (klayers,nbndsw)
      REAL(KIND=JPRB), INTENT(IN) :: PASYC(:)                  ! cloud asymmetry parameter
                                                               !   Dimensions: (klayers,nbndsw)
      REAL(KIND=JPRB), INTENT(IN) :: POMGC(:)                  ! cloud single scattering albedo
                                                               !   Dimensions: (klayers,nbndsw)
      REAL(KIND=JPRB), INTENT(IN) :: PTAUCORIG(:)              ! cloud optical depth, non-delta scaled
                                                               !   Dimensions: (klayers,nbndsw)
      REAL(KIND=JPRB), INTENT(IN) :: PTAUA(:)                  ! aerosol optical depth
                                                               !   Dimensions: (klayers,nbndsw)
      REAL(KIND=JPRB), INTENT(IN) :: PASYA(:)                  ! aerosol asymmetry parameter
                                                               !   Dimensions: (klayers,nbndsw)
      REAL(KIND=JPRB), INTENT(IN) :: POMGA(:)                  ! aerosol single scattering albedo

      REAL(KIND=JPRB), INTENT(IN) :: PTAUR(:)                  !rayleigh optical depth
      REAL(KIND=JPRB), INTENT(IN) :: PTAUG(:)                  !gas optical depth. Here O3 only! 
 
! ------- Output -------
                                                               !   All Dimensions: (klayers+1)
      REAL(KIND=JPRB), INTENT(OUT) :: PBBCD(:)
      REAL(KIND=JPRB), INTENT(OUT) :: PBBCU(:)
      REAL(KIND=JPRB), INTENT(OUT) :: PBBFD(:)
      REAL(KIND=JPRB), INTENT(OUT) :: PBBFU(:)
      REAL(KIND=JPRB), INTENT(OUT) :: PBBFDDIR(:)
      REAL(KIND=JPRB), INTENT(OUT) :: PBBCDDIR(:)


!inactive output, reverted to local (got rid of the intent(out))
      REAL(KIND=JPRB) :: ZUVCD(KLAYERS+1)
      REAL(KIND=JPRB) :: ZUVFD(KLAYERS+1)
      REAL(KIND=JPRB) :: ZUVCDDIR(KLAYERS+1)
      REAL(KIND=JPRB) :: ZUVFDDIR(KLAYERS+1)

      REAL(KIND=JPRB) :: ZNICD(KLAYERS+1)
      REAL(KIND=JPRB) :: ZNIFD(KLAYERS+1)
      REAL(KIND=JPRB) :: ZNICDDIR(KLAYERS+1)
      REAL(KIND=JPRB) :: ZNIFDDIR(KLAYERS+1)

! Output - inactive                                            !   All Dimensions: (klayers+1)
!      real(kind=jprb), intent(out) :: puvcu(:)
!      real(kind=jprb), intent(out) :: puvfu(:)
!      real(kind=jprb), intent(out) :: pnicu(:)
!      real(kind=jprb), intent(out) :: pnifu(:)
!      real(kind=jprb), intent(out) :: pvscd(:)
!      real(kind=jprb), intent(out) :: pvscu(:)
!      real(kind=jprb), intent(out) :: pvsfd(:)
!      real(kind=jprb), intent(out) :: pvsfu(:)


! ------- Local -------

      LOGICAL :: LLRTCHKCLR(KLAYERS),LLRTCHKCLD(KLAYERS)

      INTEGER(KIND=JPIM)  :: ILEV
      INTEGER(KIND=JPIM) :: IB1, IB2, IGT, IKL
      INTEGER(KIND=JPIM) :: IW, JB, JG, JK
!      integer(kind=jpim), parameter :: nuv = ?? 
      INTEGER(KIND=JPIM), PARAMETER ::  IGPTSW=1 
      INTEGER(KIND=JPIM) :: ITIND

      REAL(KIND=JPRB) :: ZTBLIND, ZE1
      REAL(KIND=JPRB) :: ZCLEAR, ZCLOUD
      REAL(KIND=JPRB) :: ZDBT(KLAYERS+1), ZDBT_NODEL(KLAYERS+1)
      REAL(KIND=JPRB) :: ZGCC(KLAYERS), ZGCO(KLAYERS)
      REAL(KIND=JPRB) :: ZOMCC(KLAYERS), ZOMCO(KLAYERS)
      REAL(KIND=JPRB) :: ZRDND(KLAYERS+1), ZRDNDC(KLAYERS+1)
      REAL(KIND=JPRB) :: ZREF(KLAYERS+1), ZREFC(KLAYERS+1), ZREFO(KLAYERS+1)
      REAL(KIND=JPRB) :: ZREFD(KLAYERS+1), ZREFDC(KLAYERS+1), ZREFDO(KLAYERS+1)
      REAL(KIND=JPRB) :: ZRUP(KLAYERS+1), ZRUPD(KLAYERS+1)
      REAL(KIND=JPRB) :: ZRUPC(KLAYERS+1), ZRUPDC(KLAYERS+1)
      REAL(KIND=JPRB) :: ZTAUC(KLAYERS), ZTAUO(KLAYERS)
      REAL(KIND=JPRB) :: ZTDBT(KLAYERS+1)
      REAL(KIND=JPRB) :: ZTRA(KLAYERS+1), ZTRAC(KLAYERS+1), ZTRAO(KLAYERS+1)
      REAL(KIND=JPRB) :: ZTRAD(KLAYERS+1), ZTRADC(KLAYERS+1), ZTRADO(KLAYERS+1)
      REAL(KIND=JPRB) :: ZDBTC(KLAYERS+1), ZTDBTC(KLAYERS+1)
      REAL(KIND=JPRB) :: ZINCFLX, ZDBTC_NODEL(KLAYERS+1) 
      REAL(KIND=JPRB) :: ZTDBT_NODEL(KLAYERS+1), ZTDBTC_NODEL(KLAYERS+1)

      REAL(KIND=JPRB) :: ZDBTMC, ZDBTMO, ZF
      REAL(KIND=JPRB) :: ZWF, ZTAUORIG, ZREPCLC
!     real(kind=jprb) :: zincflux                                   ! inactive

! Arrays from rrtmg_sw_taumoln routines

!      real(kind=jprb) :: ptaug(klayers,16), ztaur(klayers,16)
!      real(kind=jprb) :: zsflxzen(16)



! Arrays from rrtmg_sw_vrtqdr routine

      REAL(KIND=JPRB) :: ZCD(KLAYERS+1,IGPTSW), ZCU(KLAYERS+1,IGPTSW)
      REAL(KIND=JPRB) :: ZFD(KLAYERS+1,IGPTSW), ZFU(KLAYERS+1,IGPTSW)

!for the exp lookup tables I use the arrays already available in IFS
!a bit dirty, it needs cleanup for final version
      REAL(KIND=JPRB) :: ZTBLINT,ZOD_LO
      REAL(KIND=JPRB) , DIMENSION(0:10000) :: ZEXP_TBL

! Inactive arrays
!     real(kind=jprb) :: zbbcd(klayers+1), zbbcu(klayers+1)
!     real(kind=jprb) :: zbbfd(klayers+1), zbbfu(klayers+1)
!     real(kind=jprb) :: zbbfddir(klayers+1), zbbcddir(klayers+1)

! ------------------------------------------------------------------

! Initializations

      ZOD_LO=RODLOW
      ZTBLINT=RTBLINT
      ZEXP_TBL=TRANS

      IB1 = KSTART
      IB2 = KEND
      ILEV = KLAYERS
      IW = 0
      ZREPCLC = 1.E-12_JPRB
!      zincflux = 0.0_jprb

      DO JK=1,ILEV+1
         PBBCD(JK)=0._JPRB
         PBBCU(JK)=0._JPRB
         PBBFD(JK)=0._JPRB
         PBBFU(JK)=0._JPRB
         PBBCDDIR(JK)=0._JPRB
         PBBFDDIR(JK)=0._JPRB
         ZUVCD(JK)=0._JPRB
         ZUVFD(JK)=0._JPRB
         ZUVCDDIR(JK)=0._JPRB
         ZUVFDDIR(JK)=0._JPRB
         ZNICD(JK)=0._JPRB
         ZNIFD(JK)=0._JPRB
         ZNICDDIR(JK)=0._JPRB
         ZNIFDDIR(JK)=0._JPRB
      ENDDO


! Calculate the optical depths for gaseous absorption and Rayleigh scattering

      
! Top of shortwave spectral band loop, jb = 16 -> 29; ibm = 1 -> 14
!UVmod (AB 2014): only one "band" between 250nm and 400nm. The g-points become here the
!monochromatic wavelengths. Their number depends on the resolution chosen in UVRAD.F90
!the call is done already in the spectral loop so no need here of any g-point loop

      DO JB = IB1, IB2

!        do jk=1,ilev+1
!           zbbcd(jk)=0.0_jprb
!           zbbcu(jk)=0.0_jprb
!           zbbfd(jk)=0.0_jprb
!           zbbfu(jk)=0.0_jprb
!        enddo

! Top of g-point interval loop within each band (iw is cumulative counter) 
         DO JG = 1,1 !igt
            IW = IW+1

! Apply adjustments for correct Earth/Sun distance and zenith angle to incoming solar flux
!            zincflx(iw) = padjflux(jb) * zsflxzen(iw) * prmu0
!already provided adjusted flux in input, computed in UVRAD
            ZINCFLX = PADJFLUX 
!             zincflux = zincflux + padjflux(jb) * zsflxzen(iw) * prmu0           ! inactive

! Compute layer reflectances and transmittances for direct and diffuse sources, 
! first clear then cloudy

! zrefc(jk)  direct albedo for clear
! zrefo(jk)  direct albedo for cloud
! zrefdc(jk) diffuse albedo for clear
! zrefdo(jk) diffuse albedo for cloud
! ztrac(jk)  direct transmittance for clear
! ztrao(jk)  direct transmittance for cloudy
! ztradc(jk) diffuse transmittance for clear
! ztrado(jk) diffuse transmittance for cloudy
!  
! zref(jk)   direct reflectance
! zrefd(jk)  diffuse reflectance
! ztra(jk)   direct transmittance
! ztrad(jk)  diffuse transmittance
!
! zdbtc(jk)  clear direct beam transmittance
! zdbto(jk)  cloudy direct beam transmittance
! zdbt(jk)   layer mean direct beam transmittance
! ztdbt(jk)  total direct beam transmittance at levels

! Clear-sky    
!   TOA direct beam    
            ZTDBTC(1)=1.0_JPRB
            ZTDBTC_NODEL(1)=1.0_JPRB
!   Surface values
            ZDBTC(ILEV+1) =0.0_JPRB
            ZTRAC(ILEV+1) =0.0_JPRB
            ZTRADC(ILEV+1)=0.0_JPRB
            ZREFC(ILEV+1) =PALBP
            ZREFDC(ILEV+1)=PALBD
            ZRUPC(ILEV+1) =PALBP
            ZRUPDC(ILEV+1)=PALBD
           
! Cloudy-sky    
!   Surface values
            ZTRAO(ILEV+1) =0.0_JPRB
            ZTRADO(ILEV+1)=0.0_JPRB
            ZREFO(ILEV+1) =PALBP
            ZREFDO(ILEV+1)=PALBD
           
! Total sky    
!   TOA direct beam    
            ZTDBT(1)=1.0_JPRB
            ZTDBT_NODEL(1)=1.0_JPRB
!   Surface values
            ZDBT(ILEV+1) =0.0_JPRB
            ZTRA(ILEV+1) =0.0_JPRB
            ZTRAD(ILEV+1)=0.0_JPRB
            ZREF(ILEV+1) =PALBP
            ZREFD(ILEV+1)=PALBD
            ZRUP(ILEV+1) =PALBP
            ZRUPD(ILEV+1)=PALBD
    

! Top of layer loop
            DO JK=1,ILEV

! Note: two-stream calculations proceed from top to bottom; 
!   RRTMG_SW quantities are given bottom to top and are reversed here
! UV PROCESSOR: NO NEED TO REVERSE, INPUT ALREADY TOP-DOWN
               IKL=JK !ilev+1-jk

! Set logical flag to do REFTRA calculation
!   Do REFTRA for all clear layers
               LLRTCHKCLR(JK)=.TRUE.

!   Do REFTRA only for cloudy layers in profile, since already done for clear layers
               LLRTCHKCLD(JK)=.FALSE.
               LLRTCHKCLD(JK)=(PCLFR(IKL) > ZREPCLC)

! Clear-sky optical parameters - this section inactive     
!   Original
!               ztauc(jk) = ztaur(ikl,iw) + ptaug(ikl,iw)
!               zomcc(jk) = ztaur(ikl,iw) / ztauc(jk)
!               zgcc(jk) = 0.0001_jprb
!   Total sky optical parameters        
!               ztauo(jk) = ztaur(ikl,iw) + ptaug(ikl,iw) + ptauc(ikl,ibm)
!               zomco(jk) = ptauc(ikl,ibm) * pomgc(ikl,ibm) + ztaur(ikl,iw)
!               zgco (jk) = (ptauc(ikl,ibm) * pomgc(ikl,ibm) * pasyc(ikl,ibm) + &
!                           ztaur(ikl,iw) * 0.0001_jprb) / zomco(jk)
!               zomco(jk) = zomco(jk) / ztauo(jk)

! Clear-sky optical parameters including aerosols
               ZTAUC(JK) = PTAUR(IKL) + PTAUG(IKL) + PTAUA(IKL)
               ZOMCC(JK) = PTAUR(IKL) * 1.0_JPRB + PTAUA(IKL) * POMGA(IKL)
               ZGCC(JK) = PASYA(IKL) * POMGA(IKL) * PTAUA(IKL) / ZOMCC(JK)
               ZOMCC(JK) = ZOMCC(JK) / ZTAUC(JK)

! Pre-delta-scaling clear and cloudy direct beam transmittance (must use 'orig', unscaled cloud OD)       
! \/\/ This block of code is only needed for unscaled direct beam calculation - NOT USED IN UV COMPUTATIONS
               IF (KDELM == 0) THEN
!     
                  ZCLEAR = 1.0_JPRB - PCLFR(IKL)
                  ZCLOUD = PCLFR(IKL)

! Clear
!                   zdbtmc = exp(-ztauc(jk) / prmu0)
 
! Use exponential lookup table for transmittance, or expansion of exponential for low tau
                  ZE1 = ZTAUC(JK) / PRMU0
                  IF (ZE1 <= ZOD_LO) THEN
                     ZDBTMC = 1._JPRB - ZE1 + 0.5_JPRB * ZE1 * ZE1
                  ELSE 
                     ZTBLIND = ZE1 / (BPADE + ZE1)
                     ITIND = ZTBLINT * ZTBLIND + 0.5_JPRB
                     ZDBTMC = ZEXP_TBL(ITIND)
                  ENDIF

                  ZDBTC_NODEL(JK) = ZDBTMC
                  ZTDBTC_NODEL(JK+1) = ZDBTC_NODEL(JK) * ZTDBTC_NODEL(JK)

! Clear + Cloud
                  ZTAUORIG = ZTAUC(JK) + PTAUCORIG(IKL)
!                   zdbtmo = exp(-ztauorig / prmu0)

! Use exponential lookup table for transmittance, or expansion of exponential for low tau
                  ZE1 = ZTAUORIG / PRMU0
                  IF (ZE1 <= ZOD_LO) THEN
                     ZDBTMO = 1._JPRB - ZE1 + 0.5_JPRB * ZE1 * ZE1
                  ELSE
                     ZTBLIND = ZE1 / (BPADE + ZE1)
                     ITIND = ZTBLINT * ZTBLIND + 0.5_JPRB
                     ZDBTMO = ZEXP_TBL(ITIND)
                  ENDIF

                  ZDBT_NODEL(JK) = ZCLEAR * ZDBTMC + ZCLOUD * ZDBTMO
                  ZTDBT_NODEL(JK+1) = ZDBT_NODEL(JK) * ZTDBT_NODEL(JK)

               ENDIF
! /\/\ Above code only needed for unscaled direct beam calculation - NOT USED IN UV COMPUTATIONS


! Delta scaling - clear   
               ZF = ZGCC(JK) * ZGCC(JK)
               ZWF = ZOMCC(JK) * ZF
               ZTAUC(JK) = (1.0_JPRB - ZWF) * ZTAUC(JK)
               ZOMCC(JK) = (ZOMCC(JK) - ZWF) / (1.0_JPRB - ZWF)
               ZGCC (JK) = (ZGCC(JK) - ZF) / (1.0_JPRB - ZF)

! Total sky optical parameters (cloud properties already delta-scaled)
!   Use this code if cloud properties are derived in rrtmg_sw_cldprop       
               IF (KCPR >= 1) THEN
                  ZTAUO(JK) = ZTAUC(JK) + PTAUC(IKL)
                  ZOMCO(JK) = ZTAUC(JK) * ZOMCC(JK) + PTAUC(IKL) * POMGC(IKL) 
                  ZGCO (JK) = (PTAUC(IKL) * POMGC(IKL) * PASYC(IKL) +&
                              & ZTAUC(JK) * ZOMCC(JK) * ZGCC(JK)) / ZOMCO(JK)
                  ZOMCO(JK) = ZOMCO(JK) / ZTAUO(JK)

! Total sky optical parameters (if cloud properties not delta scaled)
!   Use this code if cloud properties are not derived in rrtmg_sw_cldprop       
               ELSEIF (KCPR == 0) THEN
                  ZTAUO(JK) = PTAUR(IKL) + PTAUG(IKL) + PTAUA(IKL) + PTAUC(IKL)
                  ZOMCO(JK) = PTAUA(IKL) * POMGA(IKL) + PTAUC(IKL) * POMGC(IKL) +&
                              & PTAUR(IKL) * 1.0_JPRB
                  ZGCO (JK) = (PTAUC(IKL) * POMGC(IKL) * PASYC(IKL) +&
                              & PTAUA(IKL)*POMGA(IKL)*PASYA(IKL)) / ZOMCO(JK)
                  ZOMCO(JK) = ZOMCO(JK) / ZTAUO(JK)

! Delta scaling - clouds 
!   Use only if subroutine rrtmg_sw_cldprop is not used to get cloud properties and to apply delta scaling
                  ZF = ZGCO(JK) * ZGCO(JK)
                  ZWF = ZOMCO(JK) * ZF
                  ZTAUO(JK) = (1._JPRB - ZWF) * ZTAUO(JK)
                  ZOMCO(JK) = (ZOMCO(JK) - ZWF) / (1.0_JPRB - ZWF)
                  ZGCO (JK) = (ZGCO(JK) - ZF) / (1.0_JPRB - ZF)
               ENDIF 

! End of layer loop
            ENDDO    


! Clear sky reflectivities
            CALL REFTRA_SW(YDERAD,ILEV,LLRTCHKCLR,ZGCC,PRMU0,ZTAUC,ZOMCC,ZREFC,ZREFDC,ZTRAC,ZTRADC)

! Total sky reflectivities      
            CALL REFTRA_SW(YDERAD,ILEV,LLRTCHKCLD,ZGCO,PRMU0,ZTAUO,ZOMCO,ZREFO,ZREFDO,ZTRAO,ZTRADO)


            DO JK=1,ILEV

! Combine clear and cloudy contributions for total sky
               IKL = JK !ilev+1-jk again, no need to reverse cloud cover

! AB here clear and cloud fractions are recomputed taking into account
! the total cloud cover computed with one of the methods in the driver uvradi
               IF (PCLEAR < 1.0_JPRB) THEN
                  ZCLOUD = PCLFR(IKL)/(1._JPRB-PCLEAR)
                  ZCLEAR = 1._JPRB-PCLFR(IKL)/(1._JPRB-PCLEAR)
               ELSE
                  ZCLOUD = 0.0_JPRB
                  ZCLEAR = 1.0_JPRB
               ENDIF
               ZCLOUD = MAX(0.0_JPRB,MIN(1.0_JPRB,ZCLOUD))
               ZCLEAR = MAX(0.0_JPRB,MIN(1.0_JPRB,ZCLEAR))

               ZREF(JK) = ZCLEAR*ZREFC(JK) + ZCLOUD*ZREFO(JK)
               ZREFD(JK)= ZCLEAR*ZREFDC(JK) + ZCLOUD*ZREFDO(JK)
               ZTRA(JK) = ZCLEAR*ZTRAC(JK) + ZCLOUD*ZTRAO(JK)
               ZTRAD(JK)= ZCLEAR*ZTRADC(JK) + ZCLOUD*ZTRADO(JK)

! Direct beam transmittance        

! Clear
!                zdbtmc = exp(-ztauc(jk) / prmu0)

! Use exponential lookup table for transmittance, or expansion of 
! exponential for low tau
               ZE1 = ZTAUC(JK) / PRMU0
               IF (ZE1 <= ZOD_LO) THEN
                  ZDBTMC = 1._JPRB - ZE1 + 0.5_JPRB * ZE1 * ZE1
               ELSE
                  ZTBLIND = ZE1 / (BPADE + ZE1)
                  ITIND = ZTBLINT * ZTBLIND + 0.5_JPRB
                  ZDBTMC = ZEXP_TBL(ITIND)
               ENDIF

               ZDBTC(JK) = ZDBTMC
               ZTDBTC(JK+1) = ZDBTC(JK)*ZTDBTC(JK)

! Clear + Cloud
!                zdbtmo = exp(-ztauo(jk) / prmu0)

! Use exponential lookup table for transmittance, or expansion of 
! exponential for low tau
               ZE1 = ZTAUO(JK) / PRMU0
               IF (ZE1 <= ZOD_LO) THEN
                  ZDBTMO = 1._JPRB - ZE1 + 0.5_JPRB * ZE1 * ZE1
               ELSE
                  ZTBLIND = ZE1 / (BPADE + ZE1)
                  ITIND = ZTBLINT * ZTBLIND + 0.5_JPRB
                  ZDBTMO = ZEXP_TBL(ITIND)
               ENDIF

               ZDBT(JK) = ZCLEAR*ZDBTMC + ZCLOUD*ZDBTMO
               ZTDBT(JK+1) = ZDBT(JK)*ZTDBT(JK)
        
            ENDDO           
                 
! Vertical quadrature for clear-sky fluxes

            CALL VRTQDR_SW (ILEV, IW,&
                            & ZREFC, ZREFDC, ZTRAC, ZTRADC,&
                            & ZDBTC, ZRDNDC, ZRUPC, ZRUPDC, ZTDBTC,&
                            & ZCD, ZCU)
      
! Vertical quadrature for cloudy fluxes

            CALL VRTQDR_SW (ILEV, IW,&
                            & ZREF, ZREFD, ZTRA, ZTRAD,&
                            & ZDBT, ZRDND, ZRUP, ZRUPD, ZTDBT,&
                            & ZFD, ZFU)

! Upwelling and downwelling fluxes at levels
!   Two-stream calculations go from top to bottom; 
!   layer indexing is reversed to go bottom to top for output arrays

            DO JK=1,ILEV+1
               IKL=ILEV+2-JK 

! Accumulate spectral fluxes over whole spectrum  
               PBBFU(IKL) = PBBFU(IKL) + ZINCFLX*ZFU(JK,IW)
               PBBFD(IKL) = PBBFD(IKL) + ZINCFLX*ZFD(JK,IW)
               PBBCU(IKL) = PBBCU(IKL) + ZINCFLX*ZCU(JK,IW)
               PBBCD(IKL) = PBBCD(IKL) + ZINCFLX*ZCD(JK,IW)
               IF (KDELM == 0) THEN 
                  PBBFDDIR(IKL) = PBBFDDIR(IKL) + ZINCFLX*ZTDBT_NODEL(JK)
                  PBBCDDIR(IKL) = PBBCDDIR(IKL) + ZINCFLX*ZTDBTC_NODEL(JK)
               ELSEIF (KDELM == 1) THEN
                  PBBFDDIR(IKL) = PBBFDDIR(IKL) + ZINCFLX*ZTDBT(JK)
                  PBBCDDIR(IKL) = PBBCDDIR(IKL) + ZINCFLX*ZTDBTC(JK)
               ENDIF
               
               PBBFDDIR(IKL) = PBBFDDIR(IKL)*(1._JPRB-PCLEAR)&
                    & + PBBCDDIR(IKL)*PCLEAR
               PBBFU(IKL) = PBBFU(IKL)*(1._JPRB-PCLEAR) + PBBCU(IKL)*PCLEAR
               PBBFD(IKL) = PBBFD(IKL)*(1._JPRB-PCLEAR) + PBBCD(IKL)*PCLEAR
            ENDDO

! End loop on jg, g-point interval
         ENDDO             

! End loop on jb, spectral band
      ENDDO                    

      END SUBROUTINE SPCVRT_SW

      END MODULE RRTMG_SW_SPCVRT


