! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE BASCOE_HETCONST(YGFL,KTRACER,PTEMP, PRS, PDENS, LD_PSC_POSSIBLE, KTOP_PSC,KBOT_PSC, KLEV,  &
  & PCONC, PSA_SIZEDIST, PAER, PTSTEP, PRHET )


!**   DESCRIPTION
!     ----------
!    Update the heterogeneous chemical reaction rates prhet(:)
!    Here we also set the surface area density (SAD) of PSC particles 
!                                          (types NAT and ICE)
!
!**   AUTHORS
!     -------
!     pre-2014: BASCOE team *BIRA*
!     2014-02:  VINCENT HUIJNEN (VH) *KNMI* implemented into IFS
!               with old param "v4q30" for PSC SAD depending only on T
!
!     MODIFICATIONS.
!     --------------
!     2016   :  SIMON CHABRILLAT (SC) *BIRA, implemented by VH:
!               new param "cifts" for PSC SAD depending on saturation p of H2O and HNO3
!               as described in Huijnen et al. (GMD, 2016)
!     2021-12: SC and YVES CHRISTOPHE (YC) *BIRA* : 
!              - removed old PSC param "v4q30" (was unrealistic and not used since 2014)
!              - updated to BASCOE 09.00.01/sb15d as developed at BIRA
!                 including Shi et al. (JGR,2001) param to compute PRHET(1:2) 
!                 (i.e. rates for ClONO2+H2O(c) & ClONO2+HCl(c)) 
!
!-----------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM,    JPRB, JPRD
USE YOMLUN   , ONLY : NULERR
#ifdef WITH_COMPO_DR_HOOK
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
#endif
USE YOM_YGFL , ONLY : TYPE_GFLD
USE BASCOE_MODULE , ONLY :  NHET, NBINS, PTSIZE, NAER, IAER_SAD, &
  & RCHEM_SAD_NAT_PSC,RCHEM_SSR_NAT_PSC,RCHEM_SAD_ICE_PSC,RCHEM_SSR_ICE_PSC

IMPLICIT NONE

!-----------------------------------------------------------------------
!*       0.1  ARGUMENTS
!             ---------
TYPE(TYPE_GFLD)    ,INTENT(IN)  :: YGFL
INTEGER(KIND=JPIM), INTENT(IN)  ::KTRACER(9)
REAL(KIND=JPRB),INTENT(IN)      :: PTEMP, PRS, PDENS
LOGICAL,INTENT(IN)              :: LD_PSC_POSSIBLE
INTEGER(KIND=JPIM), INTENT(IN)  ::KLEV,KTOP_PSC,KBOT_PSC
REAL(KIND=JPRB),INTENT(IN)      :: PCONC(YGFL%NCHEM)
REAL(KIND=JPRB),INTENT(IN)      :: PSA_SIZEDIST(NBINS)
REAL(KIND=JPRB),INTENT(IN)      :: PAER(NAER)
REAL(KIND=JPRB),INTENT(IN)      :: PTSTEP
REAL(KIND=JPRB),INTENT(OUT)     :: PRHET(NHET)


! * LOCAL 
#ifdef WITH_COMPO_DR_HOOK
REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE
#endif
INTEGER(KIND=JPIM) :: IH2O,IN2O5,IHCL,IHOCL,ICLONO2,IHOBR,IHBR,IBRONO2,IHNO3
INTEGER(KIND=JPIM) :: IHET, IERR
REAL(KIND=JPRB)    :: ZLCLONO2, ZLN2O5, ZLHOCL, ZLBRONO2, ZLHOBR, ZLHCL, ZLHBR, ZRATIO
REAL(KIND=JPRB)    :: ZVMRWV, ZVMRHCL, ZVMRCLONO2
REAL(KIND=JPRD)    :: ZWT_H2SO4, ZMOLALITY_H2SO4   ! NOTE trying double-prec !!
REAL(KIND=JPRB), DIMENSION(NBINS) :: ZND
REAL(KIND=JPRB)    :: ZSLOW1,ZSLOW2,ZSLOW3
REAL(KIND=JPRB)    :: ZP_ICE,ZVMR_H2O,ZVMR_HNO3,ZPW,ZPN0,ZPN0T,ZMT,ZBT,ZHNO3EQ,ZTEMP 
REAL(KIND=JPRB),PARAMETER :: ZPREF = 101325._JPRB

REAL(KIND=JPRB)  :: ZS_STS,ZS_NAT,ZS_ICE,ZGH2O,ZGHCL,ZG2HCL,ZHCL_VMR,ZHBR_VMR
REAL(KIND=JPRB)  :: ZRH1, ZRH2, ZRH3, ZRH4, ZRH5, ZRH6, ZRH7, ZRH8, ZRH9

! Fudging parameters: for now keeping same values as in IFS(CBA) 2016-2021
REAL(KIND=JPRB), PARAMETER :: ZL_SLOW = 1.1_JPRB           ! factor to slow down destruction of HCl and HBr when already depleted

!  Switches to choose methods used in this routine:
!       IMODE_HETCONST        method to calc PRHET, especially w.r.t. uptake coeffs on liquid sulfate aerosols
!           1 = 'sb15b':  as used in IFS from 2014 to 2021 including Hanson & Ravishankara (1994: GLIQ) for PRHET(1:2)
!           2 = 'sb15d':  as in BASCOE 09.00.01/sb15d including Shi et al. (2001) param for PRHET(1:2)
!
INTEGER(KIND=JPIM),PARAMETER :: IMODE_HETCONST = 2_JPIM

!-------------------------------------------------------------------
#include "bascoe_gliq.intfb.h"
#include "bascoe_ga_shi.intfb.h"
!-------------------------------------------------------------------
#ifdef WITH_COMPO_DR_HOOK
IF (LHOOK) CALL DR_HOOK('BASCOE_HETCONST',0,ZHOOK_HANDLE )
#endif

PRHET(:) = 0._JPRB

IF( IMODE_HETCONST == 2_JPIM .AND. PRS < 500._JPRB ) THEN
#ifdef WITH_COMPO_DR_HOOK
  IF (LHOOK) CALL DR_HOOK('BASCOE_HETCONST',1,ZHOOK_HANDLE)
#endif
 RETURN  ! All Het chem skipped above 5 hPa ! 
ENDIF

  ZS_STS = 0._JPRB      ! Surface area dens (cm2/cm3) of liquid (sulfate) aerosols
  ZS_NAT = 0._JPRB      ! Surface area dens (cm2/cm3) of NAT PSC
  ZS_ICE = 0._JPRB      ! Surface area dens (cm2/cm3) of ICE PSC
  ZGH2O  = 0._JPRB      ! uptake coeff (gamma) for ClONO2 + H2O(c)
  ZGHCL  = 0._JPRB      ! uptake coeff (gamma) for ClONO2 + HCl(c)
  ZG2HCL = 0._JPRB      ! uptake coeff (gamma) for HOCl + HCl(c)

  IH2O   =KTRACER(1)
  IN2O5  =KTRACER(2)
  IHCL   =KTRACER(3)
  IHOCL  =KTRACER(4)
  ICLONO2=KTRACER(5)
  IHOBR  =KTRACER(6)
  IHBR   =KTRACER(7)
  IBRONO2=KTRACER(8)
  IHNO3  =KTRACER(9)
 
! ----------------------------------------------------------------------
!    Prepare surf area dens and nb dens of liquid sulfate aerosols 
! ----------------------------------------------------------------------
ZS_STS = MAX( 0.0_JPRB, PAER(IAER_SAD) )    ! cm2/cm3, read by BASCOE_GS_LIQ
ZND(:) = PSA_SIZEDIST(1:NBINS)              ! particles/kg air

! ----------------------------------------------------------------------
!    Calculate uptake coefficients on liquid sulfate aerosols: 
!     ZGH2O  for ClONO2 + H2O(c)
!     ZGHCL  for ClONO2 + HCl(c)
!     ZG2HCL for HOCl   + HCl(c)
! ----------------------------------------------------------------------

! Convert from PCONC [molec/cm3] to ZVMR [mole/mole]  through division by PDENS, enforcing strict range
ZVMRWV  = MAX( MIN( PCONC(IH2O)/PDENS, 6.E-6_JPRB ), 1.E-7_JPRB )
ZVMRHCL = MAX( MIN( PCONC(IHCL)/PDENS, 4.E-9_JPRB ), 1.E-10_JPRB )

IF( IMODE_HETCONST == 1_JPIM ) THEN          ! param 'gliq' i.e. Hanson & Ravishankara (1994) as in IFS(CBA) 2014-2021
  ZS_STS  = MIN( 5.E-7_JPRB, ZS_STS )
  CALL BASCOE_GLIQ( ZND(1:NBINS), PTSIZE(1:NBINS,1) ,&
                  & PTEMP, PRS, ZVMRWV, ZVMRHCL, ZGHCL, ZGH2O )
  ZG2HCL = 0.4_JPRB * ZGHCL

ELSEIF( IMODE_HETCONST == 2_JPIM ) THEN      ! param 'sb15d' i.e. Shi et al. (JGR,2001)
  IF( ALL(ZND(:) <= 1.E-2_JPRB )) ZS_STS  = 0._JPRB
  ZVMRCLONO2 = MAX( MIN( PCONC(ICLONO2)/PDENS, 4.E-6_JPRB ), 1.E-13_JPRB )
  CALL WT_M_H2SO4( PTEMP, PRS, ZVMRWV, ZWT_H2SO4, ZMOLALITY_H2SO4 )
  CALL BASCOE_GA_SHI( ZND(1:NBINS), PTSIZE(1:NBINS,1),         &
                    & PTEMP, PRS, ZVMRWV, ZVMRHCL, ZVMRCLONO2,  &
                    & ZWT_H2SO4, ZMOLALITY_H2SO4, ZGHCL, ZGH2O , ZG2HCL )
ELSE
  CALL ABOR1( 'BASCOE_HETCONST: invalid IMODE_HETCONST' )
ENDIF

! ----------------------------------------------------------------------
!    Set surface area density for PSC particles (types NAT and ICE)
! ----------------------------------------------------------------------
IF( .NOT. LD_PSC_POSSIBLE .OR.&
  &    KLEV < KTOP_PSC .OR. KLEV > KBOT_PSC ) THEN
  ZS_ICE = 0.0_JPRB 
  ZS_NAT = 0.0_JPRB
ELSE

  ZTEMP = PTEMP

  ! Koop and Murphy, QJRMS, 2005:
  ZP_ICE = EXP( 9.550426_JPRB - 5723.265_JPRB/ZTEMP + 3.53068_JPRB*LOG(ZTEMP) - 0.00728332_JPRB*ZTEMP ) 

  !* Compute H2O vmr from concentration in molec/cm3, through division by PDENS, i.e. air density in molec/cm3
  ZVMR_H2O=PCONC(IH2O)/PDENS
  IF( ZVMR_H2O*PRS > RCHEM_SSR_ICE_PSC*ZP_ICE ) THEN       ! All ICE PSC ICE, no NAT PSC
    ZS_ICE = RCHEM_SAD_ICE_PSC
    ZS_NAT = 0._JPRB

  ELSE                                   ! no ICE PSC ; any NAT PSC ?

    !* Compute HNO3 vmr from given concentration PCONC in units [molec/cm3]
    ! partial pressure of water vapor normalized by the standard pressure ZPREF ( / )
    !    2e-5 mb (=2e-3 Pa) < pw * ZPREF < 2e-3 mb (= 2e-1 Pa)
    ! pw * ZPREF > 2e-5 mb (= 2e-3 Pa), see line above
    ! "virtual" HNO3 partial pressure normalized by the standard pressure ZPREF ( / ).

    ZVMR_HNO3 = PCONC(IHNO3)/PDENS
    ZPW = ZVMR_H2O*PRS/ZPREF                              
    ZPW = MAX(MIN(ZPW,2.E-1_JPRB/ZPREF),2.E-3_JPRB/ZPREF) 
    ZPN0 = ZVMR_HNO3*PRS/ZPREF                            

    ! Conversion of normalized pressure from "/" to "Torr/Pa". Rem: 1 Pa = 7.5e-3 Torr
    ! Parameter m(T) in Eq.(1) of Hanson and Mauersberger (1988), p.857.                                                      
    ZTEMP = MIN(ZTEMP,273._JPRB)
    ZPN0T = ZPN0*ZPREF/100._JPRB*0.75_JPRB           
    ZMT  = -2.7836_JPRB - 0.00088_JPRB*ZTEMP         
    ! Parameter b(T) in Eq.(1) of Hanson and Mauersberger (1988), p.857.
    ZBT  = 38.9855_JPRB - 11397.0_JPRB/ZTEMP + 0.009179_JPRB*ZTEMP       
    ! HNO3 partial pressure (Torr/Pa !!!). Hanson and Mauersberger (1988), Eq.(1), p.857.
    ZHNO3EQ = 10.0_JPRB**(ZMT*LOG10(ZPW*ZPREF/100.*0.75_JPRB) + ZBT)     

    IF (ZPN0T > RCHEM_SSR_NAT_PSC*ZHNO3EQ) THEN               ! yes: NAT PSC
      ZS_ICE = 0._JPRB
      ZS_NAT = RCHEM_SAD_NAT_PSC
    ELSE                                    ! no PSC
      ZS_ICE = 0._JPRB
      ZS_NAT = 0._JPRB
    ENDIF
  ENDIF
ENDIF

! ----------------------------------------------------------------------
!    Update Heter. Reaction rates.
!    a) Compute pseudo-1st order heterogeneous reaction rates depending on chosen method
!    WARNING: RH9_SB15B is WRONG as it is for another reaction; kept only for back-compatibility!
! ----------------------------------------------------------------------
ZRH1 = RH1( PTEMP, ZGH2O,ZS_STS,ZS_NAT,ZS_ICE )
ZRH2 = RH2( PTEMP, ZGHCL,ZS_STS,ZS_NAT,ZS_ICE )
ZRH3 = RH3( PTEMP,       ZS_STS,ZS_NAT,ZS_ICE )
ZRH4 = RH4( PTEMP,       ZS_STS,ZS_NAT,ZS_ICE )
ZRH5 = RH5( PTEMP,ZG2HCL,ZS_STS,ZS_NAT,ZS_ICE )

IF( IMODE_HETCONST == 1_JPIM ) THEN          ! param 'SB15B' (as in IFS(CBA) 2014-2021)
    ZRH6 = RH6_SB15B( PTEMP,ZS_STS,ZS_NAT,ZS_ICE )
    ZRH7 = RH7_SB15B( PTEMP,ZS_STS,ZS_NAT,ZS_ICE )
    ZRH8 = 0._JPRB
    ZRH9 = RH9_SB15B( PTEMP,ZS_STS,ZS_NAT,ZS_ICE )  ! BUG kept only for back-compatibility!

ELSEIF( IMODE_HETCONST == 2_JPIM ) THEN      ! param 'sb15d' (i.e newer param in CTM 9.0.1)
    ZRH6 = RH6_SB15D( PTEMP,ZWT_H2SO4,ZS_STS,ZS_NAT,ZS_ICE )
    ZRH7 = RH7_SB15D( PTEMP,          ZS_STS,ZS_NAT,ZS_ICE )
    ZRH8 = RH8_SB15D( PTEMP,          ZS_STS,ZS_NAT,ZS_ICE )
    ZRH9 = RH9_SB15D( PTEMP,PRS,PDENS,PCONC(IHOCL),PCONC(IHBR),ZS_STS,ZS_NAT,ZS_ICE )
ENDIF

! -----------------------------------------------------------------------
!    Update Heter. Reaction rates.
!  b) Compute lhcl: number of molecules of HCl which would be lost 
!     in one timestep due to hetero chemistry, assuming its abundance is infinite
! -----------------------------------------------------------------------
ZLCLONO2 = PCONC(ICLONO2)*ZRH2 / (1.0_JPRB+ZRH2*PTSTEP)
ZLN2O5   = PCONC(IN2O5)*  ZRH4 / (1.0_JPRB+ZRH4*PTSTEP)
ZLHOCL   = PCONC(IHOCL)*  ZRH5 / (1.0_JPRB+ZRH5*PTSTEP)
ZLHOBR   = PCONC(IHOBR)*  ZRH7 / (1.0_JPRB+ZRH7*PTSTEP)
IF( IMODE_HETCONST == 1_JPIM ) THEN
  ZLBRONO2 = PCONC(IBRONO2)*ZRH9 / (1.0_JPRB+ZRH9*PTSTEP) ! BUG kept only for back-compatibility!
ELSEIF( IMODE_HETCONST == 2_JPIM ) THEN
  ZLBRONO2 = 0._JPRB    ! chem schemes sb15b and sb15d actually do not have BrONO2+HCl(c)
ENDIF

ZLHCL  = (ZLCLONO2+ZLN2O5+ZLHOCL+ZLBRONO2+ZLHOBR)*PTSTEP
      
! -----------------------------------------------------------------------
!    Finish and apply limitations of rates due to pseudo 1rst order reaction
! Two cases are possible: 
! a) There is enough HCl available wrt lhcl : divide by HCl number density 
!     number density since RH* are pseudo 1st order reaction rates 
! b) There is not enough HCl  -> slow down the reactions
!     to prevent them from leading to negative vmr for H2O, ClONO2, HOCl, HOBr.
!     i.e. the rates of these reactions are divided by lhcl ( instead of n(HCl)
!     Note the FUDGING FACTOR zl_slow (historically set to 1.1) 
!        to slow down the reactions further
!        ###### WARNING: This may need adjustment: 
!           model seems to underestimate HCl severely in PSC conditions
! -----------------------------------------------------------------------

!VH Set minimum to HCL /HBR concentrations - required for troposphere with low HCL due to wet dep.
ZHCL_VMR = MAX(PCONC(IHCL)/PDENS, 1E-10_JPRB)
ZHBR_VMR = MAX(PCONC(IHBR)/PDENS, 1E-15_JPRB)

ZRATIO   = ZLHCL/(ZHCL_VMR*PDENS)

IF(ZRATIO > 1._JPRB)THEN                              ! case (b)
  ZSLOW1=1.0_JPRB/(ZL_SLOW*ZLHCL)
ELSE                                                  ! case (a)
  ZSLOW1=1.0_JPRB/(ZHCL_VMR*PDENS)
ENDIF

! Strict inhibition of N2O5+HCl(c): removed from BASCOE CTM at sb15d
IF( IMODE_HETCONST == 1_JPIM .AND. PCONC(IN2O5) < 100._JPRB ) THEN   
  ZSLOW2 = 0.0_JPRB
ELSE
  ZSLOW2 = ZSLOW1
ENDIF

PRHET(1) = ZRH1
PRHET(2) = ZRH2*ZSLOW1
PRHET(3) = ZRH3
PRHET(4) = ZRH4*ZSLOW2
PRHET(5) = ZRH5*ZSLOW1
PRHET(6) = ZRH6
PRHET(7) = ZRH7*ZSLOW1

! Same for HBr (but simply inhibit reaction)
! None of this is necessary for H2O because it is not written as reactant
!  in the reaction list (its abundance assumed infinite)
IF( IMODE_HETCONST == 1_JPIM ) THEN    ! Old 'SB15B' approach does not account for HOCl+HBr(c)
  IF( (PCONC(IHOBR)) < 100._JPRB ) THEN
    ZSLOW3 = 0._JPRB
  ELSE 
    ZLHOBR=PTSTEP*PCONC(IHOBR)*ZRH8/(1.0_JPRB+ZRH8*PTSTEP)
    ZRATIO=ZLHOBR/PCONC(IHOBR)
    IF(ZRATIO > 1._JPRB) THEN
      ZSLOW3=1.0_JPRB/(ZL_SLOW*ZRATIO*ZHBR_VMR*PDENS)
    ELSE
      ZSLOW3=1.0_JPRB/(ZHBR_VMR*PDENS)
    ENDIF
  ENDIF 
  PRHET(8) = ZRH8*ZSLOW3
  PRHET(9) = ZRH9*ZSLOW1   ! BUG kept for back-compatibility
  
ELSEIF( IMODE_HETCONST == 2_JPIM ) THEN ! New 'SB15D' approach
  ZSLOW3 = 0.0_JPRB
  ZLHOBR = PCONC(IHOBR)*  ZRH8 / (1.0_JPRB+ZRH8*PTSTEP)
  ZLHOCL = PCONC(IHOCL)*  ZRH9 / (1.0_JPRB+ZRH9*PTSTEP)
  ZLHBR = (ZLHOCL+ZLHOBR)*PTSTEP
  IF( ZHBR_VMR > 1.E-20_JPRB .AND. ZLHBR > 0._JPRB &
    & .AND. PCONC(IHBR) > 10._JPRB * ZLHBR ) THEN
    ZSLOW3 = 1.0_JPRB / PCONC(IHBR)
  ENDIF

  PRHET(8) = ZRH8*ZSLOW3
  PRHET(9) = ZRH9*ZSLOW3

!  IF( LD_PSC_POSSIBLE .AND. KLEV>=KTOP_PSC .AND. KLEV<=KBOT_PSC ) THEN 
!    IERR = 0_JPIM
!    DO IHET = 1, NHET
!      IF( ISNAN(PRHET(IHET)) .OR. PRHET(IHET)<0._JPRB .OR. PRHET(IHET)>1.e4_JPRB ) THEN    
!        WRITE(NULERR,'(a,i3,a,es12.2)') 'BASCOE_HETCONST debugging - ERROR: PRHET( '  &
!           &   ,IHET,') reached *abnormal* value: ',PRHET(IHET)
!        IERR = 1_JPIM
!      ENDIF
!    ENDDO
!    IF(IERR>0_JPIM) CALL ABOR1( 'BASCOE_HETCONST: abnormal PRHET(1:9) at testing level (47)' )
!  ENDIF
  PRHET = MAX( 0._JPRB, PRHET )
  PRHET = MIN( 1._JPRB, PRHET ) ! ceiling added in April 2022
!  WHERE( ISNAN(PRHET(:)) .OR. PRHET(:) > 1._JPRB )  !  ### NOTE ARBITRARY CEILING ! (why not more than 1!?)
!    PRHET(:) = 0._JPRB
!  ENDWHERE

ENDIF


#ifdef WITH_COMPO_DR_HOOK
IF (LHOOK) CALL DR_HOOK('BASCOE_HETCONST',1,ZHOOK_HANDLE)
#endif

! ----------------------------------------------------------------------
CONTAINS
! ----------------------------------------------------------------------

SUBROUTINE WT_M_H2SO4(PTEMP,PPAIR,PVMRH2O,PWT,PM)
!     H2SO4 mass percentage (wt%) calculation for binary solutions (H2SO4/H2O)
!     where: - ZP0H2O is the saturation water vapor (mbar)
!            - PWT is the H2SO4 mass percentage (wt%)
!            - m is the  molality of H2SO4 (mol/kg)
!     see Tabazadeh et al. (1997)
! NOTE double-prec: parameters below require double-prec -> output may need double-prec
!                   -> all calculations here are done in double-prec!

#ifdef WITH_COMPO_DR_HOOK
  USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
#endif
  USE PARKIND1 , ONLY : JPRB, JPRD
  IMPLICIT NONE
  REAL(KIND=JPRB), INTENT(IN)  :: PTEMP, PPAIR, PVMRH2O
  REAL(KIND=JPRD), INTENT(OUT) :: PWT, PM
  INTEGER(KIND=JPIM) :: I
  REAL(KIND=JPRD) :: ZY1, ZY2, ZAW, ZPH2O, ZP0H2O
  REAL(KIND=JPRD), DIMENSION(3), PARAMETER :: &
       &   ZA1 = (/  12.37208932_JPRD,    11.820654354_JPRD ,   -180.06541028_JPRD /),  &
       &   ZA2 = (/  13.455394705_JPRD,   12.891938068_JPRD,    -176.95814097_JPRD /),  &
       &   ZB1 = (/  -0.16125516114_JPRD, -0.20786404244_JPRD,    -0.38601102596_JPRD /), &
       &   ZB2 = (/  -0.191312255_JPRD,   -0.23233847708_JPRD,    -0.36257048154_JPRD /), &
       &   ZC1 = (/ -30.490657554_JPRD,   -4.807306373_JPRD,     -93.317846778_JPRD /),  &
       &   ZC2 = (/ -34.285174607_JPRD,   -6.4261237757_JPRD,    -90.469744201_JPRD /),  &
       &   ZD1 = (/  -2.1133114241_JPRD,  -5.1727540348_JPRD,    273.88132245_JPRD /),   &
       &   ZD2 = (/  -1.7620073078_JPRD,  -4.900547139_JPRD,     267.45509988_JPRD /)
#ifdef WITH_COMPO_DR_HOOK
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

  IF (LHOOK) CALL DR_HOOK('BASCOE_HETCONST:WT_M_H2SO4',0,ZHOOK_HANDLE)
#endif

  ZPH2O=PVMRH2O*PPAIR/101325._JPRD                 
  ZP0H2O=EXP(18.452406985_JPRD-3505.1578807_JPRD/PTEMP  &
            & - 330918.55082_JPRD/(PTEMP**2)      &
            & + 12725068.262_JPRD/(PTEMP**3))

  ! WATER ACTIVITY
  ZAW = 1013.25_JPRD*ZPH2O/ZP0H2O                      
  ZAW = MIN(ZAW,1._JPRD)
  IF (ZAW<=0.05_JPRD) THEN
    I = 1_JPIM
  ELSEIF (0.05_JPRD<ZAW.AND.ZAW<0.85_JPRD)THEN
    I = 2_JPIM
  ELSE
    I = 3_JPIM
  ENDIF

  ZY1 = ZA1(I)*ZAW**ZB1(I) + ZC1(I)*ZAW + ZD1(I)
  ZY2 = ZA2(I)*ZAW**ZB2(I) + ZC2(I)*ZAW + ZD2(I)

  ! MOLALITY OF H2SO4 (MOL/KG)
  PM= ZY1 + (PTEMP-190.0_JPRD)*(ZY2-ZY1)/70.0_JPRD

  ! H2SO4 MASS PERCENTAGE (WT%)
  PWT = 9800.0_JPRD*PM/(98.0_JPRD*PM+1000.0_JPRD)

#ifdef WITH_COMPO_DR_HOOK
  IF (LHOOK) CALL DR_HOOK('BASCOE_HETCONST:WT_M_H2SO4',1,ZHOOK_HANDLE)
#endif

END SUBROUTINE WT_M_H2SO4


! ----------------------------------------------------------------------
!   Compute reaction rates for heterogeneous reactions
!     List of references:
!           JPL%Y:          i.e. JPL evaluation for year %Y
!           Shi2001:        Shi et al., JGR, 2001
!           HR94:           Hanson and Ravishankara, 1994
!           Crowley:        Crowley et al., ACP, 2010
!           Ammann:         Ammann et al., ACP, 2013
!           Dennis:         Dennison et al., GMD, 2019
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
!   RH1(TEMP) ; {ClONO2 + H2O(c) = HOCl + HNO3(c)}

FUNCTION RH1(PTEMP,PGH2O,PS_STS,PS_NAT,PS_ICE)

#ifdef WITH_COMPO_DR_HOOK      
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
#endif
      USE PARKIND1 , ONLY : JPRB
      IMPLICIT NONE
      REAL(KIND=JPRB) :: RH1
      REAL(KIND=JPRB),INTENT(IN) :: PTEMP,PGH2O,PS_STS,PS_NAT,PS_ICE
      REAL(KIND=JPRB) :: ZG_STS, ZG_NAT, ZG_ICE, ZGAMMA, ZTERM
      REAL(KIND=JPRB), PARAMETER :: ZMASS_CLONO2 = 97.4579    ! molar mass of ClONO2 (g/mol) 
#ifdef WITH_COMPO_DR_HOOK
      REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

      IF (LHOOK) CALL DR_HOOK('BASCOE_HETCONST:RH1',0,ZHOOK_HANDLE)
#endif

      !     0.25*mean thermal velocity of ClONO2 
      ZTERM=0.25*SQRT(8.E7*1.38*PTEMP/(3.1415*ZMASS_CLONO2/6.022))
      ZG_STS=PGH2O  ! factor surf area dens already in PGH2O
      ZG_NAT=0.004  ! JPL2000 (p.47); same in JPL2015, p.5-63
      ZG_ICE=0.3    ! JPL2015, p.5-63
      ZGAMMA=ZG_STS+ZG_NAT*PS_NAT+ZG_ICE*PS_ICE
      RH1=ZGAMMA*ZTERM
#ifdef WITH_COMPO_DR_HOOK
  IF (LHOOK) CALL DR_HOOK('BASCOE_HETCONST:RH1',1,ZHOOK_HANDLE)
#endif

END FUNCTION RH1

! ----------------------------------------------------------------------
!   RH2(TEMP) ; {ClONO2 + HCl(c) = Cl2 + HNO3(c)}
!
! *WARNING* : RH2 taken into account outside chem integ
FUNCTION RH2(PTEMP,PGHCL,PS_STS,PS_NAT,PS_ICE)


#ifdef WITH_COMPO_DR_HOOK
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
#endif
      USE PARKIND1 , ONLY : JPRB
      IMPLICIT NONE
      REAL(KIND=JPRB) :: RH2
      REAL(KIND=JPRB),INTENT(IN) :: PTEMP,PGHCL,PS_STS,PS_NAT,PS_ICE

      REAL(KIND=JPRB) :: ZG_STS, ZG_NAT, ZG_ICE, ZGAMMA, ZTERM
      REAL(KIND=JPRB), PARAMETER :: ZMASS_CLONO2 = 97.4579    ! molar mass of ClONO2 (g/mol) 
#ifdef WITH_COMPO_DR_HOOK
      REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

      IF (LHOOK) CALL DR_HOOK('BASCOE_HETCONST:RH2',0,ZHOOK_HANDLE)
#endif

      !     0.25*mean thermal velocity of ClONO2 
      ZTERM=0.25*SQRT(8.E7*1.38*PTEMP/(3.1415*ZMASS_CLONO2/6.022))
      ZG_STS=PGHCL  ! factor surf area dens already in PGHCL
      ZG_NAT=0.2    ! JPL2015, p.5-63
      ZG_ICE=0.3    ! JPL2000 (p.47); same in JPL2015, p.5-63
      ZGAMMA=ZG_STS+ZG_NAT*PS_NAT+ZG_ICE*PS_ICE
      RH2=ZTERM*ZGAMMA   

#ifdef WITH_COMPO_DR_HOOK
  IF (LHOOK) CALL DR_HOOK('BASCOE_HETCONST:RH2',1,ZHOOK_HANDLE)
#endif


END FUNCTION RH2

! ----------------------------------------------------------------------
!   RH3(TEMP) ; {N2O5 + H2O(c) = 2HNO3(c)}
!
FUNCTION RH3(PTEMP,PS_STS,PS_NAT,PS_ICE)

#ifdef WITH_COMPO_DR_HOOK
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
#endif
      USE PARKIND1 , ONLY : JPRB
      IMPLICIT NONE
      REAL(KIND=JPRB) :: RH3
      REAL(KIND=JPRB),INTENT(IN) :: PTEMP,PS_STS,PS_NAT,PS_ICE
      REAL(KIND=JPRB) :: ZG_STS, ZG_NAT, ZG_ICE, ZGAMMA, ZTERM
      REAL(KIND=JPRB), PARAMETER :: ZMASS_N2O5 = 108.0104 ! molar mass of N2O5 (g/mol) 
#ifdef WITH_COMPO_DR_HOOK
      REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

      IF (LHOOK) CALL DR_HOOK('BASCOE_HETCONST:RH3',0,ZHOOK_HANDLE)
#endif

      !     0.25*mean thermal velocity of N2O5
      ZTERM=0.25*SQRT(8.E7*1.38*PTEMP/(3.1415*ZMASS_N2O5/6.022))
      ZG_STS=0.100   ! JPL2011, p575; JPL2015, note 35, p5-86; Denni
      ZG_NAT=0.0004  ! JPL2000 (p.47) ; same in JPL2015, p.5-61
      ZG_ICE=0.02    ! JPL2015, p.5-61
      ZGAMMA=ZG_STS*PS_STS+ZG_NAT*PS_NAT+ZG_ICE*PS_ICE
      RH3=ZTERM*ZGAMMA

#ifdef WITH_COMPO_DR_HOOK
      IF (LHOOK) CALL DR_HOOK('BASCOE_HETCONST:RH3',1,ZHOOK_HANDLE)
#endif

      
END FUNCTION RH3

! ----------------------------------------------------------------------
!   RH4(TEMP) ; {N2O5 + HCl(c) = ClNO2 + HNO3(c)}
!
FUNCTION RH4(PTEMP,PS_STS,PS_NAT,PS_ICE)

#ifdef WITH_COMPO_DR_HOOK
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
#endif
      USE PARKIND1 , ONLY : JPRB
      IMPLICIT NONE
      REAL(KIND=JPRB) :: RH4
      REAL(KIND=JPRB),INTENT(IN) :: PTEMP,PS_STS,PS_NAT,PS_ICE
      REAL(KIND=JPRB) :: ZG_STS, ZG_NAT, ZG_ICE, ZGAMMA, ZTERM
      REAL(KIND=JPRB), PARAMETER :: ZMASS_N2O5 = 108.0104 ! molar mass of N2O5 (g/mol) 
#ifdef WITH_COMPO_DR_HOOK
      REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE
     
      IF (LHOOK) CALL DR_HOOK('BASCOE_HETCONST:RH4',0,ZHOOK_HANDLE)
#endif
     
      !     0.25*mean thermal velocity of N2O5
      ZTERM=0.25*SQRT(8.E7*1.38*PTEMP/(3.1415*ZMASS_N2O5/6.022))
      ZG_STS=0.00    ! JPL2000 (p.47)
      ZG_NAT=0.003   ! JPL2000 (p.47) ; same in JPL2015, p.5-61
      ZG_ICE=0.03    ! JPL2000 (p.47) ; same in JPL2015, p.5-61
      ZGAMMA=ZG_STS*PS_STS+ZG_NAT*PS_NAT+ZG_ICE*PS_ICE
      RH4=ZTERM*ZGAMMA

#ifdef WITH_COMPO_DR_HOOK
      IF (LHOOK) CALL DR_HOOK('BASCOE_HETCONST:RH4',1,ZHOOK_HANDLE)
#endif
      
END FUNCTION RH4

! ----------------------------------------------------------------------
!   RH5(TEMP) ; {HOCl + HCl(c) = Cl2 + H2O(c)}
!
FUNCTION RH5(PTEMP,PG2HCL,PS_STS,PS_NAT,PS_ICE)

#ifdef WITH_COMPO_DR_HOOK
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
#endif
      USE PARKIND1 , ONLY : JPRB
      IMPLICIT NONE
      REAL(KIND=JPRB) :: RH5
      REAL(KIND=JPRB),INTENT(IN) :: PTEMP,PG2HCL,PS_STS,PS_NAT,PS_ICE
      REAL(KIND=JPRB) :: ZG_STS, ZG_NAT, ZG_ICE, ZGAMMA, ZTERM
      REAL(KIND=JPRB), PARAMETER :: ZMASS_HOCL = 52.4603 ! molar mass of HOCl (g/mol) 
#ifdef WITH_COMPO_DR_HOOK
      REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

      IF (LHOOK) CALL DR_HOOK('BASCOE_HETCONST:RH5',0,ZHOOK_HANDLE)
#endif

      !     0.25*mean thermal velocity of HOCl
      ZTERM=0.25*SQRT(8.E7*1.38*PTEMP/(3.1415*ZMASS_HOCL/6.022))
      ZG_STS=PG2HCL  ! factor surf area dens already in PGHCL
      ZG_NAT=0.1     ! JPL2015, p.5-63
      ZG_ICE=0.2     ! JPL2015, p.5-63
      ZGAMMA=ZG_STS+ZG_NAT*PS_NAT+ZG_ICE*PS_ICE
      RH5=ZGAMMA*ZTERM

#ifdef WITH_COMPO_DR_HOOK
      IF (LHOOK) CALL DR_HOOK('BASCOE_HETCONST:RH5',1,ZHOOK_HANDLE)
#endif

      
END FUNCTION RH5

! ----------------------------------------------------------------------
!   RH6(TEMP) ; {BrONO2 + H2O(c) = HOBr+HNO3} -> Old (SB15B) or New (SB15D)
!
FUNCTION RH6_SB15B(PTEMP,PS_STS,PS_NAT,PS_ICE)

#ifdef WITH_COMPO_DR_HOOK
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
#endif
      USE PARKIND1 , ONLY : JPRB
      IMPLICIT NONE
      REAL(KIND=JPRB) :: RH6_SB15B
      REAL(KIND=JPRB),INTENT(IN) :: PTEMP,PS_STS,PS_NAT,PS_ICE
      REAL(KIND=JPRB) :: ZG_STS, ZG_NAT, ZG_ICE, ZGAMMA, ZTERM
      REAL(KIND=JPRB), PARAMETER :: ZMASS_BRONO2 = 141.9089 ! molar mass of BrONO2 (g/mol) 
#ifdef WITH_COMPO_DR_HOOK
      REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

      IF (LHOOK) CALL DR_HOOK('BASCOE_HETCONST:RH6',0,ZHOOK_HANDLE)
#endif

      !     0.25*mean thermal velocity of BrONO2
      ZTERM=0.25*SQRT(8.E7*1.38*PTEMP/(3.1415*ZMASS_BRONO2/6.022))
      ZG_STS=0.7  ! v3s13: checked JPL2000 OK (p.56)
      ZG_NAT=0.0
      ZG_ICE=0.3  ! v3s13: checked JPL2000 OK (p.47)
      ZGAMMA=ZG_STS*PS_STS+ZG_NAT*PS_NAT+ZG_ICE*PS_ICE
      RH6_SB15B=ZTERM*ZGAMMA

#ifdef WITH_COMPO_DR_HOOK
      IF (LHOOK) CALL DR_HOOK('BASCOE_HETCONST:RH6',1,ZHOOK_HANDLE)
#endif

END FUNCTION RH6_SB15B

FUNCTION RH6_SB15D(PTEMP,PWT_H2SO4,PS_STS,PS_NAT,PS_ICE)

      USE PARKIND1 , ONLY : JPRB, JPRD
      IMPLICIT NONE
      REAL(KIND=JPRB) :: RH6_SB15D
      REAL(KIND=JPRB),INTENT(IN) :: PTEMP,PS_STS,PS_NAT,PS_ICE
      REAL(KIND=JPRD),INTENT(IN) :: PWT_H2SO4
      REAL(KIND=JPRB) :: ZG_STS, ZG_NAT, ZG_ICE, ZGAMMA, ZTERM
      REAL(KIND=JPRB), PARAMETER :: ZMASS_BRONO2 = 141.9089_JPRB ! molar mass of BrONO2 (g/mol) 
     
      ! No DR_HOOK for performance

      !     0.25*mean thermal velocity of BrONO2
      ZTERM=0.25*SQRT(8.E7*1.38*PTEMP/(3.1415*ZMASS_BRONO2/6.022))
      ZG_STS = 0.11_JPRB + EXP(29.2_JPRB - 0.4_JPRB*PWT_H2SO4)   ! Ammann, p.8217
      ZG_STS = 1._JPRB/0.8_JPRB + 1._JPRB/ZG_STS         
      ZG_STS = 1._JPRB/ZG_STS
      ZG_NAT=0.0
      ZG_ICE = 5.3E-4_JPRB*EXP(1100._JPRB/PTEMP)                 ! Crowley, p. 9160
      ZGAMMA=ZG_STS*PS_STS+ZG_NAT*PS_NAT+ZG_ICE*PS_ICE
      RH6_SB15D=ZTERM*ZGAMMA

END FUNCTION RH6_SB15D

! ----------------------------------------------------------------------
!   RH7(TEMP) ; {HOBr + HCl(c) = BrCl + H2O(c)} 
!    -> Old (SB15B) or New (SB15D)   differ only wrt ZG_ICE (old 0.3, new 0.25)
!
FUNCTION RH7_SB15B(PTEMP,PS_STS,PS_NAT,PS_ICE)

#ifdef WITH_COMPO_DR_HOOK
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
#endif
      USE PARKIND1 , ONLY : JPRB
      IMPLICIT NONE
      REAL(KIND=JPRB) :: RH7_SB15B
      REAL(KIND=JPRB),INTENT(IN) :: PTEMP,PS_STS,PS_NAT,PS_ICE
      REAL(KIND=JPRB) :: ZG_STS, ZG_NAT, ZG_ICE, ZGAMMA, ZTERM
      REAL(KIND=JPRB), PARAMETER :: ZMASS_HOBR = 96.9113 ! molar mass of HOBr (g/mol) 
#ifdef WITH_COMPO_DR_HOOK
      REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE
     
      IF (LHOOK) CALL DR_HOOK('BASCOE_HETCONST:RH7',0,ZHOOK_HANDLE)
#endif

      !     0.25*mean thermal velocity of HOBr
      ZTERM=0.25*SQRT(8.E7*1.38*PTEMP/(3.1415*ZMASS_HOBR/6.022))
      ZG_STS=0.00 ! v3s13: checked JPL2000 OK (p.55)
      ZG_NAT=0.0 
      ZG_ICE=0.3  ! v3s13: checked JPL2000 OK (p.47)
      ZGAMMA=ZG_STS*PS_STS+ZG_NAT*PS_NAT+ZG_ICE*PS_ICE
      RH7_SB15B=ZTERM*ZGAMMA

#ifdef WITH_COMPO_DR_HOOK
      IF (LHOOK) CALL DR_HOOK('BASCOE_HETCONST:RH7',1,ZHOOK_HANDLE)
#endif
      
END FUNCTION RH7_SB15B

FUNCTION RH7_SB15D(PTEMP,PS_STS,PS_NAT,PS_ICE)

#ifdef WITH_COMPO_DR_HOOK
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
#endif
      USE PARKIND1 , ONLY : JPRB
      IMPLICIT NONE
      REAL(KIND=JPRB) :: RH7_SB15D
      REAL(KIND=JPRB),INTENT(IN) :: PTEMP,PS_STS,PS_NAT,PS_ICE
      REAL(KIND=JPRB) :: ZG_STS, ZG_NAT, ZG_ICE, ZGAMMA, ZTERM
      REAL(KIND=JPRB), PARAMETER :: ZMASS_HOBR = 96.9113 ! molar mass of HOBr (g/mol) 
#ifdef WITH_COMPO_DR_HOOK
      REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

      IF (LHOOK) CALL DR_HOOK('BASCOE_HETCONST:RH7_SB15D',0,ZHOOK_HANDLE)
#endif


      !     0.25*mean thermal velocity of HOBr
      ZTERM=0.25*SQRT(8.E7*1.38*PTEMP/(3.1415*ZMASS_HOBR/6.022))
      ZG_STS=0.00     ! v3s13: checked JPL2000 OK (p.55)
      ZG_NAT=0.0 
      ZG_ICE=0.25     ! Crowley, p9156, for T in [185-210]
      ZGAMMA=ZG_STS*PS_STS+ZG_NAT*PS_NAT+ZG_ICE*PS_ICE
      RH7_SB15D=ZTERM*ZGAMMA
      
#ifdef WITH_COMPO_DR_HOOK
      IF (LHOOK) CALL DR_HOOK('BASCOE_HETCONST:RH7_SB15D',1,ZHOOK_HANDLE)
#endif

END FUNCTION RH7_SB15D

! ----------------------------------------------------------------------
!   RH8(TEMP) ; {HOBr + HBr(c) = Br2 + H2O(c)} was set to ZERO in old (SB15B)
!
FUNCTION RH8_SB15D(PTEMP,PS_STS,PS_NAT,PS_ICE)

#ifdef WITH_COMPO_DR_HOOK
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
#endif
      USE PARKIND1 , ONLY : JPRB
      IMPLICIT NONE
      REAL(KIND=JPRB) :: RH8_SB15D
      REAL(KIND=JPRB),INTENT(IN) :: PTEMP,PS_STS,PS_NAT,PS_ICE
      REAL(KIND=JPRB) :: ZG_STS, ZG_NAT, ZG_ICE, ZGAMMA, ZTERM
      REAL(KIND=JPRB), PARAMETER :: ZMASS_HOBR = 96.9113 ! molar mass of HOBr (g/mol) 
#ifdef WITH_COMPO_DR_HOOK
      REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

      IF (LHOOK) CALL DR_HOOK('BASCOE_HETCONST:RH8_SB15D',0,ZHOOK_HANDLE)
#endif

      !     0.25*mean thermal velocity of HOBr
      ZTERM=0.25*SQRT(8.E7*1.38*PTEMP/(3.1415*ZMASS_HOBR/6.022))
      ZG_STS=0.25                                   ! JPL2015 p.5-1, note 113, gamma=0.25 (at pretty high T)
      ZG_NAT=0.0
      ZG_ICE = 4.8E-4_JPRB * EXP(1240._JPRB/PTEMP)  ! Crowley, p9159, for T in [180-230]
      ZGAMMA=ZG_STS*PS_STS+ZG_NAT*PS_NAT+ZG_ICE*PS_ICE
      RH8_SB15D=ZTERM*ZGAMMA

#ifdef WITH_COMPO_DR_HOOK
      IF (LHOOK) CALL DR_HOOK('BASCOE_HETCONST:RH8_SB15D',1,ZHOOK_HANDLE)
#endif
      
END FUNCTION RH8_SB15D

! ----------------------------------------------------------------------
!   RH9(TEMP)  in old (SB15B) method: should have been for 
!   {  HOCl + HBr(c) -> BrCl + H2O(c)} but was ERRONEOUSLY computed for
!   {BrONO2 + HCl(c) = BrCl + HNO3(c)} ; BUG kept only for back-compatibility!
!
FUNCTION RH9_SB15B(PTEMP,PS_STS,PS_NAT,PS_ICE)

#ifdef WITH_COMPO_DR_HOOK
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
#endif
      USE PARKIND1 , ONLY : JPRB
      IMPLICIT NONE
      REAL(KIND=JPRB) :: RH9_SB15B
      REAL(KIND=JPRB),INTENT(IN) :: PTEMP,PS_STS,PS_NAT,PS_ICE
      REAL(KIND=JPRB) :: ZG_STS, ZG_NAT, ZG_ICE, ZGAMMA, ZTERM
#ifdef WITH_COMPO_DR_HOOK
      REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE
     
      IF (LHOOK) CALL DR_HOOK('BASCOE_HETCONST:RH9_SB15B',0,ZHOOK_HANDLE)
#endif
     
      ZTERM=0.25*SQRT(8.E7*1.38*PTEMP/(3.1415*141.9089/6.022))
! since v3s13 (except v3s16), ZG_STS=ZG_ICE=0.5, since JPL97 p247 said ZG_STS & 
! ZG_ICE could be very high (not in JPL2000)
      ZG_STS=0.5
      ZG_NAT=0.0
      ZG_ICE=0.5
      ZGAMMA=ZG_STS*PS_STS+ZG_NAT*PS_NAT+ZG_ICE*PS_ICE
      RH9_SB15B=ZTERM*ZGAMMA

#ifdef WITH_COMPO_DR_HOOK
      IF (LHOOK) CALL DR_HOOK('BASCOE_HETCONST:RH8_SB15B',1,ZHOOK_HANDLE)
#endif

      
END FUNCTION RH9_SB15B

! ----------------------------------------------------------------------
!   RH9(TEMP) ; {HOCl + HBr(c) -> BrCl + H2O(c)} fixed implementation for SB15D
!     STS: g_STS=0.0 ** i.e. we neglect this reaction on STS because
!     n(HBr)<<n(HCl) -> insignificant w.r.t. RH5
!       JPL2015 says that a full reaction/solubility/liquid phase diffusion model is necessary 
!     NAT: g_NAT=0
!     ICE: g_ICE=0  ** Crowley, p9145
!
FUNCTION RH9_SB15D(PTEMP,PRS,PDENS,PN_HOCL,PN_HBR,PS_STS,PS_NAT,PS_ICE)

#ifdef WITH_COMPO_DR_HOOK
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
#endif
      USE PARKIND1 , ONLY : JPRB
      IMPLICIT NONE
      REAL(KIND=JPRB) :: RH9_SB15D
      REAL(KIND=JPRB),INTENT(IN) :: PTEMP,PRS,PDENS,PN_HOCL,PN_HBR,PS_STS,PS_NAT,PS_ICE
      REAL(KIND=JPRB) :: ZG_STS, ZG_NAT, ZG_ICE, ZGAMMA, ZTERM
      REAL(KIND=JPRB) :: ZK_LINC, ZK_LANGC, ZPHBR, ZHBR_S
      REAL(KIND=JPRB), PARAMETER :: ZKS    = 3.3E-15      ! cm2 molec-1
      REAL(KIND=JPRB), PARAMETER :: ZNMAX  = 3.3E14       ! molec cm-2
      REAL(KIND=JPRB), PARAMETER :: ZALPHA = 0.3          ! dimensionless
      REAL(KIND=JPRB), PARAMETER :: ZMASS_HOCL = 52.4603  ! molar mass of HOCl (g/mol) 
     
#ifdef WITH_COMPO_DR_HOOK
      REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE
     
      IF (LHOOK) CALL DR_HOOK('BASCOE_HETCONST:RH9_SB15D',0,ZHOOK_HANDLE)
#endif

      !     0.25*mean thermal velocity of HOCl
      ZTERM=0.25*SQRT(8.E7*1.38*PTEMP/(3.1415*ZMASS_HOCL/6.022))
      ZG_STS=0.0
      ZG_NAT=0.0      
      IF (PTEMP <= 220._JPRB) THEN
        ZK_LINC  = 3.6e-8*EXP(4760/PTEMP)   ! i.e. K_linC(HOCl), see Crowley, p9123
        ZK_LANGC = ZK_LINC/ZNMAX            ! i.e. K_LangC(HOCl), see Crowley, p9066
        ZPHBR = PRS* MAX(MIN( PN_HBR/PDENS, 4.E-9_JPRB),1.E-16_JPRB ) ! partial pressure of HBr, keeping its vmr in check
        IF (ZPHBR > 1.E-9_JPRB) THEN    
          ZHBR_S = 4.14E5_JPRB*PN_HBR**0.88
        ELSE
          ZHBR_S = 3.E14_JPRB
        ENDIF
        ZG_ICE = ZTERM*(1+ZK_LANGC*PN_HOCL)
        ZG_ICE = 4*ZKS*ZHBR_S*ZK_LANGC*ZNMAX / ZG_ICE
        ZG_ICE = 1./(1./ZALPHA+1./ZG_ICE)
      ELSE
        ZG_ICE = 0.0
      ENDIF      
      ZGAMMA=ZG_STS*PS_STS+ZG_NAT*PS_NAT+ZG_ICE*PS_ICE
      RH9_SB15D=ZTERM*ZGAMMA

#ifdef WITH_COMPO_DR_HOOK
      IF (LHOOK) CALL DR_HOOK('BASCOE_HETCONST:RH8_SB15D',1,ZHOOK_HANDLE)
#endif

      
END FUNCTION RH9_SB15D


END SUBROUTINE BASCOE_HETCONST
