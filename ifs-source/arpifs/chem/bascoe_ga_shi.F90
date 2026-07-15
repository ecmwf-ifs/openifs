! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE BASCOE_GA_SHI(PNDENS, PRADIUS, PTEMP, PPAIR,         &
                        & PVMRH2O, PVMRHCL, PVMRCLONO2, PWT, PM, &
                        & PGA_CLONO2_H2O, PGA_CLONO2_HCL, PGA_HOCL_HCL)
!**   DESCRIPTION
!     ----------
!       Compute gamma factor (uptake coefficient) according to 
!           parameterization by Shi et al. 2001
!       extracted from bascoe > 08.
! 
!**   AUTHOR.
!     -------
!     2021-12-10: YVES CHRISTOPHE (YC) *BIRA*
!
!YC : **WARNING** coding conventions for IFS not yet fully implemented !
!
!-----------------------------------------------------------------------

  USE PARKIND1 , ONLY : JPIM,    JPRB, JPRD
  USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
  USE BASCOE_MODULE , ONLY :  NBINS
                    
  IMPLICIT NONE

!-----------------------------------------------------------------------
!*       0.1  ARGUMENTS
!             ---------
  REAL(KIND=JPRB), INTENT(IN)  :: PNDENS(NBINS), PRADIUS(NBINS), PTEMP, PPAIR, PVMRH2O, PVMRHCL, PVMRCLONO2
  REAL(KIND=JPRD), INTENT(IN)  :: PWT, PM   ! note double-prec !!
  REAL(KIND=JPRB), INTENT(OUT) :: PGA_CLONO2_H2O, PGA_CLONO2_HCL, PGA_HOCL_HCL

!*      LOCAL 
  REAL(KIND=JPHOOK)  :: ZHOOK_HANDLE
  INTEGER(KIND=JPIM) :: I
  REAL(KIND=JPRB)    :: ZCONV, ZGAMMA_HCL, ZGAMMA_H2O, ZGAMMA_HOCL, ZA
  REAL(KIND=JPRB), PARAMETER ::  ZPI=ACOS(-1._JPRB)

#ifdef WITH_COMPO_DR_HOOK
  IF (LHOOK) CALL DR_HOOK('BASCOE_GA_SHI',0,ZHOOK_HANDLE )
#endif

  ZCONV=(0.01_JPRB/1000._JPRB)*PPAIR*28.96_JPRB/(PTEMP*8.3_JPRB)    
  PGA_CLONO2_HCL = 0._JPRB
  PGA_CLONO2_H2O = 0._JPRB
  PGA_HOCL_HCL   = 0._JPRB
  DO I = 1, SIZE(PNDENS)
    ZA = 4._JPRB*ZPI * PRADIUS(i)**2 * PNDENS(I)
    CALL GAMMA_SHI( PRADIUS(I),PTEMP,PPAIR,PVMRH2O,PVMRHCL,PVMRCLONO2, &
                  & PWT,PM, &
                  & ZGAMMA_HCL,ZGAMMA_H2O,ZGAMMA_HOCL )

    PGA_CLONO2_HCL = PGA_CLONO2_HCL + ZGAMMA_HCL *ZA*ZCONV
    PGA_CLONO2_H2O = PGA_CLONO2_H2O + ZGAMMA_H2O *ZA*ZCONV
    PGA_HOCL_HCL   = PGA_HOCL_HCL   + ZGAMMA_HOCL*ZA*ZCONV
  ENDDO

#ifdef WITH_COMPO_DR_HOOK
  IF (LHOOK) CALL DR_HOOK('BASCOE_GA_SHI',1,ZHOOK_HANDLE )
#endif
! ----------------------------------------------------------------------

CONTAINS

!=======================================================================
  SUBROUTINE GAMMA_SHI(PRADIUS,PTEMP,PPAIR,PVMRH2O,PVMRHCL,PVMRCLONO2, &
                      & PWT,PM, &
                      & PGAMMA_CLONO2_HCL, PGAMMA_CLONO2_H2O, PGAMMA_HOCL)
  !**   DESCRIPTION
  !     ----------
  !**   AUTHOR
  !     ------
  !     2021-12-10: YVES CHRISTOPHE (YC) *BIRA*
  !-----------------------------------------------------------------------

    USE PARKIND1 , ONLY : JPIM,    JPRB,     JPRD
    USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

    IMPLICIT NONE

  !-----------------------------------------------------------------------
  !*       0.1  ARGUMENTS
  !             ---------
    REAL(KIND=JPRB), INTENT(IN)  :: PRADIUS,PTEMP,PPAIR,PVMRH2O ,PVMRCLONO2,PVMRHCL
    REAL(KIND=JPRD), INTENT(IN)  :: PWT, PM   ! note double-prec !!
    REAL(KIND=JPRB), INTENT(OUT) :: PGAMMA_CLONO2_HCL, PGAMMA_CLONO2_H2O, PGAMMA_HOCL

  !*       LOCAL: single precision
    REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
    REAL(KIND=JPRB), PARAMETER  :: ZCDIF1=5e-8_JPRB,   ZH01=1.6e-6_JPRB, &
                                 & ZB01=4710.0_JPRB, ZC01=0.306_JPRB , ZD01=24.0_JPRB  ! for ClONO2
    REAL(KIND=JPRB), PARAMETER  :: ZCDIF2=6.4e-8_JPRB, ZH02=1.91e-6_JPRB, &
                                 & ZB02=5862.4_JPRB, ZC02=0.0776_JPRB, ZD02=59.18_JPRB ! for HOCl
    REAL(KIND=JPRB), PARAMETER  :: R=0.082_JPRB

  !*       LOCAL: double precision
    REAL(KIND=JPRD) :: ZPH2O, ZPHCL, ZPCLONO2, ZP0H2O
    REAL(KIND=JPRD) :: ZAW, ZA, ZTEMP0, ZV, ZDIF, ZZ1, ZZ2, ZZ3, ZDE, ZMH2SO4, ZST, ZH, ZX
    REAL(KIND=JPRD) :: ZW, ZAH, ZKH, ZKH2O, ZKHYDR, ZGHYDR, ZHHCL, ZMHCL, ZKHCL, ZL_CLONO2, ZF_CLONO2 
    REAL(KIND=JPRD) :: ZGB1, ZGB2, ZGS1 ,ZF, ZGS2, ZGB3, ZGB, ZGAMMA_CLONO2
    REAL(KIND=JPRD) :: ZD_HOCL, ZW_HOCL, ZK_HOCL, ZST_HOCL, ZH_HOCL, ZL_HOCL, ZF_HOCL, ZG_HOCL


#ifdef WITH_COMPO_DR_HOOK
    IF (LHOOK) CALL DR_HOOK('BASCOE_GA_SHI:GAMMA_SHI',0,ZHOOK_HANDLE )
#endif

    ZPH2O    = PVMRH2O*PPAIR/101325_JPRB        ! in atm   
    ZPHCL    = PVMRHCL*PPAIR/101325_JPRB        ! in atm
    ZPCLONO2 = PVMRCLONO2*PPAIR/101325_JPRB     ! in atm
    ZP0H2O = EXP( 18.452406985_JPRB - 3505.1578807_JPRB/PTEMP &
                & - 330918.55082_JPRB/(PTEMP**2) + 12725068.262_JPRB/(PTEMP**3) )

    ZAW=1013.25_JPRB*ZPH2O/ZP0H2O
    ZAW = MIN( ZAW, 1._JPRB )
  ! ---------------------------------------------------------------------------
  !    Parameters for H2SO4 solutions 
  ! ---------------------------------------------------------------------------

    !   viscosity ZV(cP)
    !--------------------------------
    ZA = 169.5_JPRB+5.18_JPRB*PWT-0.0825_JPRB*PWT**2+3.27e-3_JPRB*PWT**3
    ZTEMP0 = 144.11_JPRB+0.166_JPRB*PWT-0.015_JPRB*PWT**2+2.18e-4_JPRB*PWT**3
    ZV = ZA*PTEMP**(-1.43_JPRB)*EXP(448.0_JPRB/(PTEMP-ZTEMP0))               


    !   molarity (mol/L)
    !--------------------------------
    ZZ1 =  0.12364_JPRB -   5.6e-7_JPRB*PTEMP**2
    ZZ2 = -0.02954_JPRB + 1.814e-7_JPRB*PTEMP**2
    ZZ3 = 2.343e-3_JPRB - 1.487e-6_JPRB*PTEMP - 1.324e-8_JPRB*PTEMP**2
    ZDE = 1._JPRB + ZZ1*PM + ZZ2*PM**1.5 + ZZ3*PM**2                    !H2SO4 density (g/cm3) 


    ZMH2SO4 = ZDE*PWT/9.8_JPRB                                        

    !   H2SO4 mole fraction
    !--------------------------------
    ZX = (PWT/(PWT+(100._JPRB-PWT)*98._JPRB/18._JPRB))

    !   acid activity "in molarity"
    !--------------------------------  
    ZAH = EXP( 60.51_JPRB -0.095_JPRB*PWT + 0.0077_JPRB*PWT**2 - 1.61e-5_JPRB*PWT**3 &
             & - (1.76_JPRB+2.52e-4_JPRB*PWT**2)*SQRT(PTEMP) &
             & + (-805.89_JPRB+253.05_JPRB*PWT**0.076_JPRB)/SQRT(PTEMP) )

  ! ---------------------------------------------------------------------------
  !    Parameters for the reactions of ClONO2 
  ! ---------------------------------------------------------------------------      

    !    Solubility of ClONO2(g) into the H2SO4 solution (mol/atm)
    !      where the Setchenow coefficient (ZST) is used to formulate
    !      the solubility as a fonction of the molarity
    ! ----------------------------------------------------------------     
    ZST = ZC01+ZD01/PTEMP                                         
    ZH = ZH01*EXP(ZB01/PTEMP)*EXP(-ZST*ZMH2SO4)                     

    !     mean molecular speed w=SQRT(8*R*PTEMP/(pi*m_ClONO2)) (cm/s)
    !-----------------------------------------------------------------            
    ZW = 1474.0_JPRB*SQRT(PTEMP)                                    


    !    Diffusion coefficient of ClONO2 (cm2/s)
    !      where ZCDIF1=5d-8 cm2*cp./(K*s)
    ! ----------------------------------------------------------------
    ZDIF = ZCDIF1*PTEMP/ZV

    !     Uptake coefficient of ClONO2 in the absence of HCL  
    !       Depends on the hydrolyse constant of ClONO2(ZKHYDR) 
    !       follwing two channels of reaction :
    !         -direct hydrolyse of ClONO2 (ZKH2O)
    !          ClONO2(aq) +H2O(l)  --> HNO3(aq) + HOCl(aq)
    !         -acid catalyse channel (ClONO2 protoned) (kh) 
    !          ClONO2(aq) + H+ = HClONO2+(aq) 
    !          HClONO2+(aq) +H2O(l)  --> HNO3(aq) + H2OCl+(aq)
    !          HClONO2+(aq) +HCl(aq)  --> HNO3(aq) + Cl2(aq,g)
    !------------------------------------------------------------------   
    ZKH2O  = 1.95e10_JPRB*EXP(-2800._JPRB/PTEMP)                            
    ZKH    = 1.22e12_JPRB*EXP(-6200._JPRB/PTEMP)                              
    ZKHYDR = ZKH2O*ZAW+ZKH*ZAH*ZAW

    ZGHYDR = 4._JPRB*ZH*R*PTEMP*SQRT(ZDIF*ZKHYDR)/ZW             

    !     Constant of the reaction ClONO2+HCl --> Cl2 + HNO3
    !---------------------------------------------------------   
    ZHHCL = (0.094_JPRB-0.61_JPRB*ZX+1.2_JPRB*ZX**2) &
          & * EXP(-8.68_JPRB+(8515.0_JPRB-10718.0_JPRB*ZX**0.7_JPRB)/PTEMP)
    ZMHCL = ZHHCL*ZPHCL

    ZKHCL = 7.9e11_JPRB*ZAH*ZDIF*ZMHCL

    !     Correction of experimental measurements through the parameter f 
    !     by comparing the reactio-diffusive length of ClONO2 (ZL_CLONO2) 
    !     to the aerosol radius (if l<<PRADIUS then f=1 and the reaction is 
    !     diffusion-limited near the surface).
    !--------------------------------------------------------------------------
    ZL_CLONO2 = SQRT(ZDIF/(ZKHYDR+ZKHCL))/100._JPRD
    ZF_CLONO2 = (1._JPRB/TANH(PRADIUS/ZL_CLONO2))-ZL_CLONO2/PRADIUS

    !     Reactive uptake coefficients for the hydrolysis of ClONO2 in the 
    !     presence of HCl (ZGB1), the reaction of ClONO2 + HCl (ZGB2), both for 
    !     bulk phase process, and the reaction of ClONO2 + HCl on the surface
    !     of the aerosl (ZGS1). 
    !-------------------------------------------------------------------------- 
    ZGB1 = ZF_CLONO2*ZGHYDR*SQRT(1._JPRB+ZKHCL/ZKHYDR)                   

    ZGB2 = ZGB1*ZKHCL/(ZKHCL+ZKHYDR)                               

    ZGS1 = 66.12_JPRB*EXP(-1374.0_JPRB/PTEMP)*ZH*ZMHCL                      

    !     HCl flow correction 
    !------------------------------------
    ZF = 1._JPRB / (1._JPRB+0.612_JPRB*(ZGS1+ZGB2)*ZPCLONO2/ZPHCL)                    

    ZGS2=ZF*ZGS1
    ZGB3=ZF*ZGB2                                               
    ZGB=ZGB3+ZGB1*ZKHYDR/(ZKHCL+ZKHYDR) 
 
    !     Total reactive uptake for ClONO2 (gamma_ClONO2) and reactive uptakes 
    !     for ClONO2 + HCl and for ClONO2 + H2O 
    !--------------------------------------------------------------------------
    IF( ZGS2+ZGB < 1.e-35_JPRD ) THEN
      PGAMMA_CLONO2_HCL = 0._JPRB
      PGAMMA_CLONO2_H2O = 0._JPRB
    ELSE
      ZGAMMA_CLONO2 = 1._JPRB/(1._JPRB+1._JPRB/(ZGS2+ZGB))        ! probability of reaction for ClONO2
      PGAMMA_CLONO2_HCL = ZGAMMA_CLONO2*(ZGS2+ZGB3)/(ZGS2+ZGB)
      PGAMMA_CLONO2_H2O = ZGAMMA_CLONO2-PGAMMA_CLONO2_HCL
    ENDIF

  ! ---------------------------------------------------------------------------
  !    Parameters for the reaction of HOCl 
  ! ---------------------------------------------------------------------------   

    !    Diffusion coefficient of HOCl (cm2/s)
    !      where ZCDIF2=6.4e-8 cm2*cp./(K*s)
    ! -----------------------------------------------------------------      
    ZD_HOCL = ZCDIF2*PTEMP/ZV

    !     mean molecular speed w=SQRT(8*R*PTEMP/(pi*m_HOCl)) (cm/s)
    !------------------------------------------------------------------      
    ZW_HOCL = 2009._JPRB*SQRT(PTEMP)

    !    Solubility of HOCl(g) into the H2SO4 solution (mol/atm)
    ! -----------------------------------------------------------------        
    ZST_HOCL = ZC02+ZD02/PTEMP
    ZH_HOCL  = ZH02*EXP(ZB02/PTEMP)*EXP(-ZST_HOCL*ZMH2SO4)

    !     Reactive uptake for the reaction HOCl + HCl (PGAMMA_HOCL) where :
    !           ZK_HOCL is diffusion limited
    !           ZF_HOCL is the correction to the experimental measurments
    !           F is the HCl flow correction
    !-------------------------------------------------------------------
    ZK_HOCL = 1.25e9_JPRB*ZAH*ZD_HOCL*ZMHCL    

    ZL_HOCL = SQRT(ZD_HOCL/ZK_HOCL)/100._JPRB
    ZF_HOCL = (1/TANH(PRADIUS/ZL_HOCL)) - ZL_HOCL/PRADIUS

    ZG_HOCL = ZF_HOCL*4._JPRB*ZH_HOCL*R*PTEMP*SQRT(ZD_HOCL*ZK_HOCL)/ZW_HOCL
    IF( ZF*ZG_HOCL < 1.E-35_JPRD ) THEN
      PGAMMA_HOCL = 0._JPRB
    ELSE
      PGAMMA_HOCL = 1._JPRB/(1._JPRB+1._JPRB/(ZF*ZG_HOCL))
    ENDIF

#ifdef WITH_COMPO_DR_HOOK
    IF (LHOOK) CALL DR_HOOK('BASCOE_GA_SHI:GAMMA_SHI',1,ZHOOK_HANDLE )
#endif

  END SUBROUTINE GAMMA_SHI

! ----------------------------------------------------------------------
END SUBROUTINE BASCOE_GA_SHI                   
