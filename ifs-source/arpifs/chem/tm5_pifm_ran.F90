! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE TM5_PIFM_RAN(KLEV,PANGLE,PALB,PCST_O3_COL,PDV2,PDV3, &
      &       PTAUA_CLD_COL,PTAUS_CLD_COL,PCLD_COL, &
      &       PTAUA_AER_COL,PTAUS_AER_COL,PAER_COL,PFACT,PFRAC)


!**   DESCRIPTION 
!     ----------
!
!   Part of TM5 routines for IFS chemistry: 
!   The evaluation of photolysis rates
!
!************************************************************************
!*            PRACTICAL IMPROVED FLUX METHOD (PIFM)                     * 
!*                 to calculate actinic fluxes                          *
!*     Zdunkowski,Welch,Korb: Beitr. Phys. Atmosph. vol. 53, p. 147 ff  *
!*                                                                      *
!*           Cloud effects added using the method of :                  *
!*     Geleyn and Hollingworth, Contribs. Atms.Phys,52(1),1979          *
!*                                                                      *
!************************************************************************      
!*     This version is not suitable for calculation for conserving      *
!*     scattering (ZW0=1). ZW0 is limited to ZW0 .le. 1. - 1.E-15.         *
!*     For ZW0 = 1, ZAL(4) and ZAL(5) has to be calculated differently.    *
!
!
!**   INTERFACE.
!     ----------
!          *TM5_PIFM_RAN* IS CALLED FROM  *TM5_PHOTO_FLUX* 

! INPUTS:
! -------
! PANGLE      : zenith angle          
! PALB        : albedo
! PFRAC       : cloud fraction per layer (0->1)
! PCST_O3_COL : O3 cross sections          
! pdv2, pdv3  : differential column info
! PTAUA_CLD_COL,taus_cld_col  : optical depth clouds       
! PTAUA_AER_COL, TAUS_AER_COL : optical depth aerosols      
! PFACT       : actinic flux 
! PAER_COL    : aerosol phase functions
! PCLD_COL    : cloud phase functions
!
!
! OUTPUTS:
! -------
!
! PFACT(KLEV,NBANDS_TROP)  : actinic flux 
!
! LOCAL:
! -------
! ZPRAY   : rayleigh phase function
! ZRW     : flux array for cloudy sky
! ZRF     : flux array for clear sky 
! ZTU1    : parallel solar flux   
! ZTU2    : matrix coefficient 
!
!
!     AUTHOR.
!     -------
!
!        Jason Williams     *KNMI*
!        VINCENT HUIJNEN    *KNMI*
!         
!        TM5-community    
!
!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 2012-07-23


USE PARKIND1 , ONLY : JPIM     ,JPRB, JPRD
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOMCST   , ONLY : RPI
USE TM5_PHOTOLYSIS , ONLY : MAXWAV,NBANDS_TROP,NGRID,LMID,LMID_GRIDA,  & 
      & SZA_LIMIT, CS_RAY, NMOM, FLUX

IMPLICIT NONE 

!-----------------------------------------------------------------------
!*       0.1  ARGUMENTS
!             ---------

INTEGER(KIND=JPIM),INTENT(IN) :: KLEV
REAL(KIND=JPRB),INTENT(IN)    :: PANGLE, PALB
REAL(KIND=JPRB),INTENT(IN)    :: PFRAC(KLEV)

REAL(KIND=JPRB),INTENT(IN)    :: PTAUA_CLD_COL(KLEV), PTAUS_CLD_COL(KLEV)  ! optical depth clouds     
REAL(KIND=JPRB),INTENT(IN)    :: PTAUA_AER_COL(KLEV,NBANDS_TROP,NGRID)     ! optical depth aerosols 
REAL(KIND=JPRB),INTENT(IN)    :: PTAUS_AER_COL(KLEV,NBANDS_TROP,NGRID)     ! optical depth aerosols 
REAL(KIND=JPRB),INTENT(IN)    :: PCLD_COL(KLEV)                            ! cloud phase functions
REAL(KIND=JPRB),INTENT(IN)    :: PAER_COL(KLEV,NBANDS_TROP,NGRID)          ! aerosol phase functions

REAL(KIND=JPRB),INTENT(IN)    :: PCST_O3_COL(KLEV,MAXWAV)
REAL(KIND=JPRB),INTENT(IN)    :: PDV2(KLEV),PDV3(KLEV)

REAL(KIND=JPRB),INTENT(OUT)   :: PFACT(KLEV,NBANDS_TROP)



! * LOCAL 
REAL(KIND=JPRB),PARAMETER     :: ZA_CNST = 0.50572 
REAL(KIND=JPRB),PARAMETER     :: ZB_CNST = 6.07995 
REAL(KIND=JPRB),PARAMETER     :: ZC_CNST = 1.63640
REAL(KIND=JPRB)               :: ZGAMMA, ZU0, ZALB

REAL(KIND=JPRB)              :: ZRF(3*(KLEV+1))          !flux array for clear sky 
REAL(KIND=JPRB)              :: ZTU1(KLEV),ZTU2(KLEV)    !matrix coefficient 

REAL(KIND=JPRB)           :: ZAL(KLEV,5)
REAL(KIND=JPRB)           :: ZSD(KLEV+1,NBANDS_TROP), ZFD(KLEV+1,NBANDS_TROP), ZFU(KLEV+1,NBANDS_TROP)
REAL(KIND=JPRB)           :: ZTS_PI_CLR(KLEV,NBANDS_TROP),ZTA_PI_CLR(KLEV,NBANDS_TROP)
REAL(KIND=JPRB)           :: ZG_PI_CLR(KLEV,NBANDS_TROP), ZTS_PI_CLD(KLEV,NBANDS_TROP)
REAL(KIND=JPRB)           :: ZTA_PI_CLD(KLEV), ZG_PI_CLD(KLEV)
REAL(KIND=JPRB)           :: ZTAUA_CLR_COL(KLEV,NBANDS_TROP), ZTAUS_CLR_COL(KLEV,NBANDS_TROP) ! optical depth clear sky

REAL(KIND=JPRB)      :: ZCS_O2 

INTEGER(KIND=JPIM)   :: IGRID,JMAXLEV,JMAXLEV3,JNL,J3 ! KMAXLAY2

REAL(KIND=JPRB),DIMENSION(0:NMOM)   :: ZPRAY  ! phase functions

REAL(KIND=JPRB) :: ZW0,ZP1,ZF,ZSMOOTH1,ZSMOOTH2,ZB0,ZBU0,ZALPH1,ZALPH2,ZALPH3,ZALPH4,ZTAUTOT,ZEPS,ZFACTOR,ZEMIN
REAL(KIND=JPRB) :: ZRM,ZGAM1,ZGAM2,ZE,ZTAUSCAT,ZTD1,ZTD2,ZTDS1,ZUEPS2,ZG,ZARG,ZFACT_BOT, ZFACT_TOP

! Help variables to prevent single-precision to fail
REAL(KIND=JPRB)  :: ZEPS_IN

!     diffusivity factor
REAL(KIND=JPRB),PARAMETER ::   ZU_CNST = 2._JPRB
REAL(KIND=JPRB),PARAMETER ::   ZDELU0=  1.E-3_JPRB 
REAL(KIND=JPRB),PARAMETER ::   ZRESONC= 1.E-6_JPRB

REAL(KIND=JPRB)    :: ZAL_CLR_1,ZAL_CLR_2,ZAL_CLR_3,ZAL_CLR_4,ZAL_CLR_5   !matrix coefficient for clear sky
REAL(KIND=JPRB)    :: ZAL_CLD_1,ZAL_CLD_2,ZAL_CLD_3,ZAL_CLD_4,ZAL_CLD_5   !matrix coefficient for cloudy sky


REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

! * counters
INTEGER(KIND=JPIM) :: JL, JK

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('TM5_PIFM_RAN',0,ZHOOK_HANDLE )

!-- Ensure realistic albedo values
ZALB = MIN(PALB,0.7_JPRB)

!-----correction of the air mass factor---------------------------------
!     F. Kasten and T. Young, Revised optical air mass tabels and       
!     approximation formula (1989) Aplied Optics Vol 28, No. 22 P. 4735 
!     and J. Lenoble, atmospheric radiative transfer (1993), P. 236     
!    define air mass factor in mu

ZGAMMA = 90._JPRB - PANGLE   
ZU0 = SIN(ZGAMMA*RPI/180._JPRB) + ZA_CNST *(ZGAMMA+ZB_CNST)**(-ZC_CNST)     
ZU0 = MIN(1._JPRB,ZU0)    


!---------------------------------------------------------------------     
!     internal integer variable 

JMAXLEV  = KLEV + 1
JMAXLEV3 = 3 * JMAXLEV

! initialize

PFACT = 0._JPRB ; ZSD = 0._JPRB ; ZFD = 0._JPRB ;  ZFU = 0._JPRB ; ZTD1 = 0._JPRB
   
ZRF = 0._JPRB 
!VH ZRW = 0. 
ZAL = 0._JPRB
ZTU1 = 0._JPRB
ZTU2 = 0._JPRB
   
! determine phase functions
   ZPRAY(0) = 1. ; ZPRAY(1) = 0. ; ZPRAY(2)=0.1 ; ZPRAY(3:NMOM) = 0. 

DO JL  = 1,NBANDS_TROP 

  IF (PANGLE <= 71.) THEN
    JNL = LMID(JL)
    IGRID = 1_JPIM
    ZCS_O2 = 7.28E-24_JPRB
  ELSEIF (PANGLE > 71. .AND. PANGLE<=SZA_LIMIT) THEN
    JNL = LMID_GRIDA(JL)
    IGRID = 2_JPIM
    ZCS_O2 = 6.6E-24_JPRB
  ENDIF
   
  ! fill array for absorption and scattering components before performing
  ! calculations.   
  DO JK = 1,KLEV

    ZTAUS_CLR_COL(JK,JL) = CS_RAY(JNL)*1._JPRB/0.2095*PDV2(JK)
  
    IF(JNL<17) THEN
      ZTAUA_CLR_COL(JK,JL) = ZCS_O2*PDV2(JK) + PCST_O3_COL(JK,JNL)*PDV3(JK)
    ELSE
      ZTAUA_CLR_COL(JK,JL) = PCST_O3_COL(JK,JNL)*PDV3(JK)
    ENDIF     

      !VH - set minimum value to OD values, otherwise computation of Actinic fluxes
      !     is not stable in single precision. Requires better solution eventually.
      ZTS_PI_CLR(JK,JL) = ZTAUS_CLR_COL(JK,JL)+PTAUS_AER_COL(JK,JL,IGRID)
      ZTA_PI_CLR(JK,JL) = MAX(1E-7_JPRB,ZTAUA_CLR_COL(JK,JL)+PTAUA_AER_COL(JK,JL,IGRID))
      ZTS_PI_CLD(JK,JL) = PTAUS_CLD_COL(JK)
      ZTA_PI_CLD(JK) = PTAUA_CLD_COL(JK)
     
    IF (ZTAUS_CLR_COL(JK,JL)+PTAUS_AER_COL(JK,JL,IGRID)>0.) THEN
      ZG_PI_CLR(JK,JL)  = (ZPRAY(1)*ZTAUS_CLR_COL(JK,JL) +   &
        &  PAER_COL(JK,JL,IGRID)*PTAUS_AER_COL(JK,JL,IGRID))/  &
        & (ZTAUS_CLR_COL(JK,JL)+PTAUS_AER_COL(JK,JL,IGRID))       
    ELSE
      ZG_PI_CLR(JK,JL)  = 0._JPRB
    ENDIF
    ZG_PI_CLD(JK)  = PCLD_COL(JK) 
  ENDDO

  DO JK = 1,KLEV       
         
!*****  first: clear sky *******************************************
 
    ZTAUTOT = ZTS_PI_CLR(JK,JL)+ZTA_PI_CLR(JK,JL)
 
    IF (ZTAUTOT /= 0.) THEN  
       ZW0=ZTS_PI_CLR(JK,JL) / ZTAUTOT
    ELSE
       ZW0=0.
    ENDIF

    ZW0 = MIN(ZW0,1.-1E-15_JPRB)

!   ZP1: first expansion coefficient of the phase function

    ZP1=3.*ZG_PI_CLR(JK,JL)                            

!   ZF: fraction of radiation contained in diffraction peak
 
    ZF=ZG_PI_CLR(JK,JL)**2_JPIM                        

!   ZB0:  fractional mean backward scattering coefficient of diffuse light
!   ZBU0: backward scattering coefficient of primary scattered parallel solar light

!   for small ZP1 ZSMOOTH1,ZSMOOTH2 manage the ZSMOOTH change of ZB0 and
!   ZBU0 to 0

    IF (ZP1<=0.1) THEN
       ZSMOOTH1=1.33333333-ZP1*3.3333333
       ZSMOOTH2=10.*ZP1
    ELSE
       ZSMOOTH1=1._JPRB
       ZSMOOTH2=1._JPRB
    ENDIF

    ZB0=(3._JPRD-ZP1)/8._JPRD  *ZSMOOTH1                                 
    ZBU0=0.5_JPRB-ZU0/4._JPRB*(ZP1-3._JPRB*ZF)/(1._JPRB-ZF)  *ZSMOOTH2    

!   alpha coefficient

    ZALPH1=ZU_CNST*(1._JPRD-(1._JPRD-ZB0)*ZW0)     
    ZALPH2=ZU_CNST*ZB0*ZW0               
    ZALPH3=ZW0*ZBU0*(1._JPRB-ZF)
    ZALPH4=ZW0*(1._JPRB-ZBU0)*(1._JPRB-ZF)

!   epsilon and gamma coefficient

    ZEPS=SQRT(ZALPH1**2_JPIM-ZALPH2**2_JPIM)          

    ZFACTOR=1._JPRB-ZW0*ZF

!   check for resonance condition in GAM1 and GAM2, if fulfil then
!   chance U0 and calculate ZUEPS2, ZBU0, ZALPH3, ZALPH4 again.

    ZUEPS2=(ZU0*ZEPS)**2_JPIM
    ZARG = ZUEPS2-ZFACTOR**2_JPIM
    IF (ZARG < 0._JPRB) ZARG = -1._JPRB * ZARG

    IF (ZARG < ZRESONC) THEN
       IF(ZUEPS2 < ZFACTOR**2_JPIM) THEN
          ZU0=ZU0-ZDELU0
       ELSE
          ZU0=ZU0+ZDELU0
       ENDIF
       ZUEPS2=(ZU0*ZEPS)**2_JPIM
       ZBU0=0.5_JPRB-ZU0/4._JPRB*(ZP1-3._JPRB*ZF)/(1._JPRB-ZF)  *ZSMOOTH2    
       ZALPH3=ZW0*ZBU0*(1._JPRB-ZF)
       ZALPH4=ZW0*(1.-ZBU0)*(1._JPRB-ZF)
    ENDIF

    ZGAM1=( ZFACTOR*ZALPH3-ZU0*(ZALPH1*ZALPH3+ZALPH2*ZALPH4) ) * &
      &  1/(ZFACTOR**2_JPIM-ZUEPS2)
    ZGAM2=(-ZFACTOR*ZALPH4-ZU0*(ZALPH1*ZALPH4+ZALPH2*ZALPH3) ) * &
      &  1/(ZFACTOR**2_JPIM-ZUEPS2)

    ZE=EXP(-ZEPS*ZTAUTOT)
    ZRM=ZALPH2/(ZALPH1+ZEPS)

    ZEMIN=MAX(1.E-18_JPRB, 1._JPRB-ZE**2_JPIM * ZRM**2_JPIM)
    ZAL_CLR_4= ZE*(1._JPRB-ZRM**2_JPIM)/ZEMIN
    ZAL_CLR_5= ZRM*(1._JPRB-ZE**2_JPIM)/ZEMIN
    ZAL_CLR_1= EXP(-ZFACTOR*ZTAUTOT/ZU0)
               
    ZAL_CLR_2=-ZAL_CLR_4*ZGAM2-ZAL_CLR_5*ZGAM1*     &
      &  ZAL_CLR_1+ZGAM2*ZAL_CLR_1  
    ZAL_CLR_3=-ZAL_CLR_5*ZGAM2-ZAL_CLR_4*ZGAM1*     &
      &  ZAL_CLR_1+ZGAM1   
       
!******************* second: cloudy sky *****************************
! For ECMWF input there is the possibility that cloud fraction occurs without a 
! corresponding value for lwc
 

    IF(  PFRAC(JK) > 0.02 .AND. ZTS_PI_CLD(JK,JL) > 1.E-5 ) THEN

       ZTAUSCAT = ZTS_PI_CLR(JK,JL) + PTAUS_CLD_COL(JK)
       ZTAUTOT = ZTA_PI_CLR(JK,JL) + PTAUA_CLD_COL(JK)+ZTAUSCAT
       ZG      = ZG_PI_CLD(JK)*PTAUS_CLD_COL(JK)/ZTAUSCAT 
   
       IF (ZTAUTOT > 0.) THEN
          ZW0=ZTAUSCAT/ZTAUTOT    
       ELSE
          ZW0=0.
       ENDIF

       ZW0 = MIN(ZW0,1.-1E-15_JPRB)
 
       ZP1=3.*ZG
       ZF=ZG**2_JPIM
       
       IF (ZP1<0.1) THEN
          ZSMOOTH1=1.33333333-ZP1*3.3333333
          ZSMOOTH2=10.*ZP1
       ELSE
          ZSMOOTH1=1._JPRB
          ZSMOOTH2=1._JPRB
       ENDIF

       ZB0=(3._JPRD-ZP1)/8._JPRD  *ZSMOOTH1                                 
       ZBU0=0.5_JPRB-ZU0/4._JPRB*(ZP1-3._JPRB*ZF)/(1._JPRB-ZF)  *ZSMOOTH2    
       
       ZALPH1=ZU_CNST*(1._JPRD-(1._JPRD-ZB0)*ZW0)     
       ZALPH2=ZU_CNST*ZB0*ZW0               
       ZALPH3=ZW0*ZBU0*(1._JPRB-ZF)
       ZALPH4=ZW0*(1._JPRB-ZBU0)*(1._JPRB-ZF)
       ! gets negative - precision issue ??
       !ZEPS=SQRT(ZALPH1**2-ZALPH2**2)  

       ! Add small number to prevent negative value.. (1e-30 seems not to work. Try 1e-20)
       ! Add another max-statement.
       ZEPS_IN=MAX(0._JPRB,ZALPH1**2_JPIM-ZALPH2**2_JPIM)
       ZEPS=SQRT(ZEPS_IN)
       !ZEPS=SQRT(ZALPH1*ZALPH1-ZALPH2*ZALPH2)  

       ZFACTOR=1._JPRB-ZW0*ZF
       
       ZUEPS2=(ZU0*ZEPS)**2_JPIM
       ZARG = ZUEPS2-ZFACTOR**2_JPIM
       IF (ZARG < 0.) ZARG = -1._JPRB * ZARG
    
       IF (ZARG < ZRESONC) THEN 
          IF(ZUEPS2 < ZFACTOR**2_JPIM) THEN
             ZU0=ZU0-ZDELU0
          ELSE
             ZU0=ZU0+ZDELU0
          ENDIF
          ZUEPS2=(ZU0*ZEPS)**2_JPIM
          ZBU0=0.5_JPRB-ZU0/4._JPRB*(ZP1-3._JPRB*ZF)/(1._JPRB-ZF)  *ZSMOOTH2    
          ZALPH3=ZW0*ZBU0*(1._JPRB-ZF)
          ZALPH4=ZW0*(1._JPRB-ZBU0)*(1._JPRB-ZF)
       ENDIF

       ZGAM1=( ZFACTOR*ZALPH3-ZU0*(ZALPH1*ZALPH3+ZALPH2*ZALPH4) ) * 1._JPRB/(ZFACTOR**2_JPIM-ZUEPS2)
       ZGAM2=(-ZFACTOR*ZALPH4-ZU0*(ZALPH1*ZALPH4+ZALPH2*ZALPH3) ) * 1._JPRB/(ZFACTOR**2_JPIM-ZUEPS2)

       ZE=EXP(-ZEPS*ZTAUTOT)
       ZRM=ZALPH2/(ZALPH1+ZEPS)

       ZEMIN=MAX(1.E-18_JPRB, 1._JPRB-ZE**2_JPIM * ZRM**2_JPIM)
       ZAL_CLD_4=ZE *(1._JPRB-ZRM**2_JPIM)/ZEMIN
       ZAL_CLD_5=ZRM*(1._JPRB-ZE**2_JPIM )/ZEMIN
       ZAL_CLD_1=EXP(-ZFACTOR*ZTAUTOT/ZU0) 
         
       ZAL_CLD_2=-ZAL_CLD_4*ZGAM2-ZAL_CLD_5*ZGAM1*ZAL_CLD_1+ZGAM2*ZAL_CLD_1  
       ZAL_CLD_3=-ZAL_CLD_5*ZGAM2-ZAL_CLD_4*ZGAM1*ZAL_CLD_1+ZGAM1   

       ZAL(JK,1) =(1._JPRB-PFRAC(JK))*ZAL_CLR_1 + PFRAC(JK)*ZAL_CLD_1
       ZAL(JK,2) =(1._JPRB-PFRAC(JK))*ZAL_CLR_2 + PFRAC(JK)*ZAL_CLD_2
       ZAL(JK,3) =(1._JPRB-PFRAC(JK))*ZAL_CLR_3 + PFRAC(JK)*ZAL_CLD_3
       ZAL(JK,4) =(1._JPRB-PFRAC(JK))*ZAL_CLR_4 + PFRAC(JK)*ZAL_CLD_4
       ZAL(JK,5) =(1._JPRB-PFRAC(JK))*ZAL_CLR_5 + PFRAC(JK)*ZAL_CLD_5
    ELSE
       
       ZAL(JK,1) =  ZAL_CLR_1
       ZAL(JK,2) =  ZAL_CLR_2
       ZAL(JK,3) =  ZAL_CLR_3
       ZAL(JK,4) =  ZAL_CLR_4
       ZAL(JK,5) =  ZAL_CLR_5
       
    ENDIF
  ENDDO !JK
                                               
!====================================================================
! matrix inversion
!====================================================================
      

!   direct solution of the first four equations       

  ZRF(1) = ZU0*FLUX(JNL)
  ZRF(2) = 0._JPRB      
      
!  5th to 10th equation: bring matrix elements on the left of the main
!  diagonal to the rhs:  save elements on the right of the main 
!      diagonal in array -ZTU(l,1)    

  ZRF(3) = ZAL(1,3) * ZRF(1)   
  ZRF(4) = ZAL(1,1) * ZRF(1)   
  ZRF(5) = ZAL(1,2) * ZRF(1)   

  ZTU1(1) = ZAL(1,4)   
  ZTU2(1) = ZAL(1,5)   

!  blocks of 6 equations: eliminate left matrix elements, save right 
!       matrix elements in array -ZTU(l,i), 
!       calculate rhs.
        
  DO  JK=2,KLEV      
    J3=3_JPIM*JK   
    ZTD1       = 1._JPRB/(1._JPRB-ZAL(JK,5_JPIM)*ZTU2(JK-1_JPIM))      
    ZTD2       = ZAL(JK,4_JPIM)*ZTU2(JK-1_JPIM)   
    ZTU1(JK)  = ZTD1*ZAL(JK,4_JPIM)
    ZTU2(JK)  = ZAL(JK,5_JPIM) + ZTD2*ZTU1(JK) 
    ZRF(J3)  = ZTD1 * (ZAL(JK,3_JPIM)*ZRF(J3-2_JPIM) + ZAL(JK,5_JPIM)*ZRF(J3-1_JPIM)) 
    ZRF(J3+1)= ZAL(JK,1)*ZRF(J3-2_JPIM)      
    ZRF(J3+2)= ZAL(JK,2)*ZRF(J3-2_JPIM) + ZAL(JK,4_JPIM)*ZRF(J3-1_JPIM)+ZTD2*ZRF(J3)
  ENDDO

!   last two equations: the same as before 

  ZTDS1   = 1._JPRB / (1._JPRB-ZALB*ZTU2(KLEV))  
  ZRF(JMAXLEV3) = ZTDS1 * (ZALB * ZRF(JMAXLEV3-2_JPIM)+ ZALB * ZRF(JMAXLEV3-1_JPIM))
  
!   now we have created an upper triangular matrix the elements of which
!   are -ZTU(l,i), 0, or 1 (in the main diagonal). the 0 and 1 elements 
!   are not stored in an array. let us solve the system now and store the
!   results in the arrays ZRF (fluxes clear sky) and ZRW (fluxes cloudy sky)

  DO  JK=KLEV,1,-1     
    J3=3_JPIM*JK      
    ZRF(J3+2_JPIM) = ZRF(J3+2_JPIM) + ZTU2(JK)*ZRF(J3+3_JPIM) 
    ZRF(J3)   = ZRF(J3) + ZTU1(JK)*ZRF(J3+3_JPIM)   
    ZSD(JK+1_JPIM,JL)  = ZRF(J3+1_JPIM)
    ZFD(JK+1_JPIM,JL)  = ZRF(J3+2_JPIM) 
    ZFU(JK+1_JPIM,JL)  = ZRF(J3+3_JPIM)
  ENDDO ! JK  

  ZSD(1,JL)  = ZRF(1) 
  ZFD(1,JL)  = ZRF(2) 
  ZFU(1,JL)  = ZRF(3)
       
! calculate the actinic flux                                     

  ZFACT_TOP = ZSD(1_JPIM,JL)/ZU0 + ZU_CNST * ZFD(1_JPIM,JL)   + ZU_CNST * ZFU(1_JPIM,JL)

  DO JK = 1,KLEV
    ZFACT_BOT = ZSD(JK+1,JL)/ZU0 + ZU_CNST * ZFD(JK+1,JL)   + ZU_CNST * ZFU(JK+1,JL)
    PFACT(JK,JL) = MAX(0._JPRB ,(ZFACT_TOP + ZFACT_BOT)/2._JPRB)
    ZFACT_TOP = ZFACT_BOT 
  ENDDO ! JK

ENDDO ! wavelength
     

IF (LHOOK) CALL DR_HOOK('TM5_PIFM_RAN',1,ZHOOK_HANDLE )
END SUBROUTINE TM5_PIFM_RAN

