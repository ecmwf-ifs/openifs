! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE BASCOE_PSC_PARAM( YGFL,KIDIA, KFDIA, KLON,  KTRACER, PTSTEP,KTROPOP,KL, PTP, &
  & PRSF1, PCONC)

!**   DESCRIPTION 
!   Simple param for PSC effect on chemistry. Same interface as PSC_FOR_CHEM
!   Replaces full PSC calc which is PSC advection + PSC_EVOLUTION + PSC_FOR_CHEM
!                                         simonc, v3s73, Jan 2005
! ----------------------------------------------------------------------
! A simple param of PSC's s.a.d. and the permanent loss of H2O and HNO3
! through sedimentation (quick for ICEPSC and slow for NAT PSC).
! Temporary losses through condensation are not modelled -> H2O and HNO3 will
! be overestimated -> using 186K & 194K i.o. TCICE & TCNAT
! ----------------------------------------------------------------------
USE PARKIND1  ,    ONLY : JPIM,   JPRB
USE YOMHOOK   ,    ONLY : LHOOK,  DR_HOOK, JPHOOK
! USE BASCOE_MODULE, ONLY : IHNO3,IH2O
USE YOM_YGFL , ONLY : TYPE_GFLD
IMPLICIT NONE

!-----------------------------------------------------------------------
!*       0.1  ARGUMENTS
!             ---------
TYPE(TYPE_GFLD)   , INTENT(INOUT) :: YGFL
INTEGER(KIND=JPIM), INTENT(IN)    :: KIDIA, KFDIA, KLON
INTEGER(KIND=JPIM), INTENT(IN)    :: KTRACER(2)
REAL(KIND=JPRB),    INTENT(IN)    :: PTSTEP
INTEGER(KIND=JPIM), INTENT(IN)    :: KTROPOP(KLON)
INTEGER(KIND=JPIM), INTENT(IN)    :: KL
REAL(KIND=JPRB),    INTENT(IN)    :: PTP(KLON),PRSF1(KLON)
REAL(KIND=JPRB),    INTENT(INOUT) :: PCONC(KLON,YGFL%NCHEM)


!-----------------------------------------------------------------------
!  Local variables
!-----------------------------------------------------------------------
REAL(KIND=JPRB),PARAMETER :: ZT_ICE = 186._JPRB, ZT_NAT = 194._JPRB ! Temperatures
REAL(KIND=JPRB),PARAMETER :: ZSEDIM_NAT = 1._JPRB/(100._JPRB*86400._JPRB) ! 1/s, char time: 100days ! used at 3s80b
REAL(KIND=JPRB),PARAMETER :: ZSEDIM_NAT_2 = 1._JPRB/(20._JPRB*86400._JPRB) ! 1/s, char time: 10days ! used with IMODE2
REAL(KIND=JPRB),PARAMETER :: ZSEDIM_ICE = 1._JPRB/(  9._JPRB*86400._JPRB) ! 1/s, char time: 9 days
REAL(KIND=JPRB),PARAMETER :: ZPREF = 101325._JPRB

REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE
INTEGER(KIND=JPIM) :: JL, IHNO3,IH2O
REAL(KIND=JPRB)    :: ZF_LOSS,ZTEMP
REAL(KIND=JPRB)    :: ZP_ICE,ZDENS,ZVMR_H2O,ZVMR_HNO3,ZPW,ZPN0,ZPN0T,ZMT,ZBT,ZHNO3EQ 
!- Switch to select the PSC parameterization. 
!-   1: very simple. 
!-   2: more realistic approach, but may need further tuning
INTEGER(KIND=JPIM) :: IMODE=2

IF (LHOOK) CALL DR_HOOK('BASCOE_PSC_PARAM',0,ZHOOK_HANDLE )

IH2O =KTRACER(1)
IHNO3=KTRACER(2)

IF ( IMODE == 1_JPIM ) THEN

DO JL=KIDIA,KFDIA
! ----------------------------------------------------------------------
!  Make sure we are only in stratosphere
! ----------------------------------------------------------------------
  IF ( KL < KTROPOP(JL) ) THEN
    ZTEMP = PTP(JL)
    IF( ZTEMP < ZT_ICE ) THEN
      ZF_LOSS = EXP( -ZSEDIM_ICE*PTSTEP )
      PCONC(JL,IH2O)  = ZF_LOSS * PCONC(JL,IH2O) ! commented out at 3s80c & 3s81c
      PCONC(JL,IHNO3) = ZF_LOSS * PCONC(JL,IHNO3) ! commented out at 3s80c & 3s81c
    ELSEIF ( ZTEMP < ZT_NAT ) THEN
      ZF_LOSS = EXP( -ZSEDIM_NAT*PTSTEP )
      ! PCONC(JL,iH2O)  = ZF_LOSS * PCONC(JL,iH2O) ! commented out for all v3s80 & v3s81 runs
      PCONC(JL,IHNO3) = ZF_LOSS * PCONC(JL,IHNO3) ! commented out at 3s80c & 3s81c
    ENDIF
  ENDIF
ENDDO  

ELSEIF (IMODE == 2_JPIM) THEN

DO JL=KIDIA,KFDIA
! ----------------------------------------------------------------------
!  Make sure we are only in stratosphere
! ----------------------------------------------------------------------
  IF ( KL < KTROPOP(JL) ) THEN
    ZTEMP = PTP(JL)
    ! Koop and Murphy, QJRMS, 2005
    ZP_ICE = EXP( 9.550426_JPRB - 5723.265_JPRB/ZTEMP + 3.53068_JPRB*LOG(ZTEMP) - 0.00728332_JPRB*ZTEMP ) 
    !*  Air density (molec/cm3) 
    ZDENS = 7.24291E16_JPRB*PRSF1(JL)/ZTEMP 
    !* Compute H2O vmr from concentration in molec/cm3
    ZVMR_H2O=PCONC(JL,IH2O)/ZDENS
    !VH IF( ZVMR_H2O*PRSF1(JL)/ZP_ICE > 1.0_JPRB ) THEN
    IF( ZVMR_H2O*PRSF1(JL) > ZP_ICE ) THEN
      ZF_LOSS = EXP( -ZSEDIM_ICE*PTSTEP )
      PCONC(JL,IH2O)  = ZF_LOSS * PCONC(JL,IH2O) ! commented out at 3s80c & 3s81c
      PCONC(JL,IHNO3) = ZF_LOSS * PCONC(JL,IHNO3) ! commented out at 3s80c & 3s81c
    ELSE
    
      !* Compute HNO3 vmr from given concentration PCONC in units [molec/cm3]
      ! partial pressure of water vapor normalized by the standard pressure ZPREF ( / )
      ! 2e-3Pa < pw * ZPREF < 2e-3 mb (= 2e-1 Pa)
      ! pw * ZPREF > 2e-5 mb (= 2e-3 Pa), see line above
      ! "virtual" HNO3 partial pressure normalized by the standard pressure ZPREF ( / ). 
      ! Note vmr(HNO3) should be for gas phase only!
      ZVMR_HNO3=PCONC(JL,IHNO3)/ZDENS
      ZPW = ZVMR_H2O*PRSF1(JL)/ZPREF                         
      ZPW = MAX(MIN(ZPW,2.E-1_JPRB/ZPREF),2.E-3_JPRB/ZPREF)  
      ! ZPW = max(pw,2.e-3/ZPREF)                            
      ZPN0 = ZVMR_HNO3*PRSF1(JL)/ZPREF                       
                                                             
      ZTEMP = MIN(ZTEMP,273._JPRB)
      ! Conversion of normalized pressure from "/" to "Torr/Pa". Rem: 1 Pa = 7.5e-3 Torr.
      ZPN0T = ZPN0*ZPREF/100._JPRB*0.75_JPRB           
      ! Parameter m(T) in Eq.(1) of Hanson and Mauersberger (1988), p.857.
      ZMT  = -2.7836_JPRB - 0.00088_JPRB*ZTEMP
      ! Parameter b(T) in Eq.(1) of Hanson and Mauersberger (1988), p.857.         
      ZBT  = 38.9855_JPRB - 11397.0_JPRB/ZTEMP + 0.009179_JPRB*ZTEMP        
      ! HNO3 partial pressure (Torr/Pa !!!). Hanson and Mauersberger (1988), Eq.(1), p.857.
      ZHNO3EQ = 10.0_JPRB**(ZMT*LOG10(ZPW*ZPREF/100.*0.75_JPRB) + ZBT)     

      !VH IF (ZPN0T/ZHNO3EQ > 1.0_JPRB) THEN
      IF (ZPN0T > ZHNO3EQ) THEN
         ZF_LOSS = EXP( -ZSEDIM_NAT_2*PTSTEP )
         ! PCONC(JL,iH2O)  = ZF_LOSS * PCONC(JL,iH2O) ! commented out for all v3s80 & v3s81 runs
         PCONC(JL,IHNO3) = ZF_LOSS * PCONC(JL,IHNO3) ! commented out at 3s80c & 3s81c
      ENDIF
    ENDIF
  ENDIF
ENDDO  

ENDIF


IF (LHOOK) CALL DR_HOOK('BASCOE_PSC_PARAM',1,ZHOOK_HANDLE )
END SUBROUTINE BASCOE_PSC_PARAM

