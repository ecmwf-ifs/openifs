! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE TM5_PHOTORATES_TROPO(KLEV,PANGLE,PCST_O3_COL,PCST_NO2_COL,PCST_HNO3_COL, &
      &     PCST_H2O2_COL,PCST_CH2O_COL,PCST_N2O5_COL,PCST_PAN_COL,PCST_NO3_COL,&
      &     PQY_NO2_COL,PQY_O1D_COL,PQY_CH3COCHO_COL,PQY_CO_COL,PQY_C2O3_COL, &
      &     PFACT,PFDIR,PRJ,PT)


!**   DESCRIPTION 
!     ----------
!
!   Part of TM5 routines for IFS chemistry: 
!   The evaluation of photolysis rates
!
!=======================================================================
!  SUBROUTINE INDEX:                                                    
!     PHOTORATES:    calculation of photolysis and heating rates         
!=======================================================================
!     References:                                                        
!     J. Landgraf and P.J. Crutzen, 1998:                               
!     An Efficient Methode for Online Calculation of Photolysis and     
!     Heating Rates, J. Atmos. Sci., 55, 863-878
! 
!      Expanded for high angles and online:
!      J.E.Williams, J. Landgraf, B. Bregman and H. H. Walter, A modified band
!      approach for the accurate calculation of online photolysis rates in
!      stratospheric-tropospheric Chemistry Transport Models, Atms. Chem. Phys., 
!      6, 4137-4161, 2006.
!
!=======================================================================
!
!**   INTERFACE.
!     ----------
!          *TM5_PHOTORATES_TROPO* IS CALLED FROM  *TM5_PHOTO_FLUX* 

! INPUTS:
! -------
!
!
! OUTPUTS:
! -------
!
! RJ ! photolysis rates 
!
! LOCAL:
! -------
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


USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
USE TM5_PHOTOLYSIS , ONLY : MAXWAV,MAXW,NBANDS_TROP,                            & 
      & NWAV_NO2,  NWAV_CO,   NWAV_HNO3, NWAV_O3, NWAV_PAN, NWAV_H2O2,              &
      & NWAV_N2O5, NWAV_CH2O, NWAV_NO3, NWAV_CH3COCHO,                              &
      & SZA_LIMIT, NPHOTO,                                                          &
      & CS_H2O2,   CS_N2O5,   CS_CH2O,     CS_O2     , CS_HONO, CS_HNO4, CS_GLYAL,  &
      & CS_GLY,    CS_CH3OOH, CS_CH3COCHO, CS_ORGN   , CS_ALD2, CS_ISPD, CS_HALD,   &
      & CS_CH3COCH3_235,CS_CH3COCH3_254,CS_CH3COCH3_263,CS_CH3COCH3_280,CS_CH3COCH3_298, &
      & QYH_HCO,   QYH2_CO,   QYNO2_O, QYNO_O2,  QY_ISPD, QYPAN, QY_GLY, QY_GLY_TOT,&
      & JO3D,      JNO2,      JH2O2,     JHNO3,  JHNO4, JN2O5,  JACH2O,  JBCH2O,    &
      & JMEPE,     JANO3,     JBNO3,     JPANA,  JPANB, JORGN,  JMENO2,   JHONO,    &
      & J45,       J74,       JROOH,  JISOPOOH,    JO2, JISPD, JA_ACET, JB_ACET,    &
      & JGLYALD, JGLYA,       JGLYB,   JHPALDA,JHPALDB, JHYAC,                       &
      & LINI,LFIN,LMID, LINI_GRIDA, LFIN_GRIDA, LMID_GRIDA  
  
USE YOMLUN   , ONLY : NULOUT
  
IMPLICIT NONE 

!-----------------------------------------------------------------------
!*       0.1  ARGUMENTS
!             ---------

INTEGER(KIND=JPIM),INTENT(IN) :: KLEV
REAL(KIND=JPRB),INTENT(IN)    :: PANGLE
REAL(KIND=JPRB),DIMENSION(KLEV,NBANDS_TROP),INTENT(IN) :: PFACT  !actinic flux
REAL(KIND=JPRB),DIMENSION(KLEV,MAXW),INTENT(IN)        :: PFDIR  !flux o2-o3 abs.
 
REAL(KIND=JPRB),DIMENSION(KLEV),INTENT(IN)              :: PT
REAL(KIND=JPRB),DIMENSION(KLEV,MAXWAV),INTENT(IN)       :: PCST_O3_COL
REAL(KIND=JPRB),DIMENSION(KLEV,NWAV_NO2),INTENT(IN)     :: PCST_NO2_COL,PQY_NO2_COL
REAL(KIND=JPRB),DIMENSION(KLEV,NWAV_CO),INTENT(IN)      :: PQY_CO_COL,PQY_C2O3_COL
REAL(KIND=JPRB),DIMENSION(KLEV,NWAV_HNO3),INTENT(IN)    :: PCST_HNO3_COL
REAL(KIND=JPRB),DIMENSION(KLEV,NWAV_O3),INTENT(IN)      :: PQY_O1D_COL
REAL(KIND=JPRB),DIMENSION(KLEV,NWAV_N2O5),INTENT(IN)    :: PCST_N2O5_COL
REAL(KIND=JPRB),DIMENSION(KLEV,NWAV_H2O2),INTENT(IN)    :: PCST_H2O2_COL
REAL(KIND=JPRB),DIMENSION(KLEV,NWAV_CH2O),INTENT(IN)    :: PCST_CH2O_COL
REAL(KIND=JPRB),DIMENSION(KLEV,NWAV_PAN),INTENT(IN)     :: PCST_PAN_COL
REAL(KIND=JPRB),DIMENSION(KLEV,NWAV_NO3),INTENT(IN)     :: PCST_NO3_COL
REAL(KIND=JPRB),DIMENSION(KLEV,NWAV_CH3COCHO),INTENT(IN):: PQY_CH3COCHO_COL

REAL(KIND=JPRB),DIMENSION(KLEV,NPHOTO),INTENT(OUT):: PRJ  !photolysis rates 

! * LOCAL 
REAL(KIND=JPRB),DIMENSION(NBANDS_TROP)  :: & 
    & ZBJO1D,    ZBJNO2,  ZBJHONO, ZBJHNO3,    ZBJCOH2,    ZBJCHOH,  ZBJN2O5, ZBJHNO4, &
    & ZBJPANA,   ZBJPANB, ZBJNO2O, ZBJNOO2,    ZBJH2O2,    ZBJCH3OOH,ZBJALD2, ZBJCH3COCHO,&
    & ZBJCH3ONO2,ZBJO2,   ZBJISPD, ZBJ_A_ACET, ZBJ_B_ACET, ZBJGLYAL, ZBJHALD, ZBJGLYA, ZBJGLYB

REAL(KIND=JPRB), DIMENSION(77)  :: ZCS_CH3COCH3            
REAL(KIND=JPRB)                 :: ZDELTA1,ZDELTA2,ZDELTA3,ZDELTA4,ZDELTA5,ZDELTA6,ZDELTA7
REAL(KIND=JPRB)                 :: ZTSCALE,ZDELTA1_SPEC,ZDELTA3_SPEC
!VH REAL(KIND=JPRB)            :: esrm2
INTEGER(KIND=JPIM),DIMENSION(7)  :: IBAND_START,IBAND_MIDDLE,IBAND_END
 

REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

! * counters
INTEGER(KIND=JPIM) :: JL, JK, JC, JJ

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
#include "abor1.intfb.h"

IF (LHOOK) CALL DR_HOOK('TM5_PHOTORATES_TROPO',0,ZHOOK_HANDLE )
  
! ==============================================================        
! CALCULATION PHOTOLYSIS RATES FOR SCATTERING ATMOSPHERE            
! ============================================================== 

   
! the check for daylight is skipped as the caluclation of rj is only called within the photo_flux loop
! where a filter already applies
 
PRJ(:,:) = 0.     
!                   
! choose the grid settings dependent on the SZA for the column       
IF(PANGLE<=71.) THEN
  IBAND_START(:) =LINI(:)
  IBAND_MIDDLE(:)=LMID(:)
  IBAND_END(:)   =LFIN(:)
ELSEIF(PANGLE>71. .AND. PANGLE<=SZA_LIMIT) THEN
  IBAND_START(:) =LINI_GRIDA(:)
  IBAND_MIDDLE(:)=LMID_GRIDA(:)
  IBAND_END(:)   =LFIN_GRIDA(:)
ENDIF
  
DO JK = 1,KLEV  
  ! re-initialize ZBJ values               
  DO JL = 1,NBANDS_TROP    
    ZBJO1D(JL)   = 0.0 
    ZBJNO2(JL)   = 0.0 
    ZBJHNO3(JL)  = 0.0 
    ZBJHONO(JL)  = 0.0
    ZBJCOH2(JL)  = 0.0 
    ZBJCHOH(JL)  = 0.0
    ZBJN2O5(JL)  = 0.0 
    ZBJHNO4(JL)  = 0.0    
    ZBJNO2O(JL)  = 0.0 
    ZBJNOO2(JL)  = 0.0 
    ZBJH2O2(JL)  = 0.0 
    ZBJCH3OOH(JL)= 0.0  
    ZBJPANA(JL)  = 0.0
    ZBJPANB(JL)  = 0.0    
    ZBJALD2(JL)  = 0.0
    ZBJCH3COCHO(JL) = 0.0  
    ZBJCH3ONO2(JL)  = 0.0
    ZBJO2(JL)    = 0.0
    ZBJISPD(JL)  = 0.0
    ZBJ_A_ACET(JL) = 0.0
    ZBJ_B_ACET(JL) = 0.0
    ZBJGLYA(JL) = 0.0
    ZBJGLYB(JL) = 0.0        
    ZBJGLYAL(JL) = 0.0
    ZBJHALD(JL) = 0.0     
  ENDDO 
                       
  ! assign cs_ch3coch3  w/ a good value
  !if(PT(JK)<=245.) ZCS_CH3COCH3=cs_ch3coch3_235
  !if(PT(JK)>245. .and. PT(JK)<=255.) ZCS_CH3COCH3=cs_ch3coch3_254
  !if(PT(JK)>255. .and. PT(JK)<=265.) ZCS_CH3COCH3=cs_ch3coch3_263
  !if(PT(JK)>265. .and. PT(JK)<=285.) ZCS_CH3COCH3=cs_ch3coch3_280
  !if(PT(JK)>285.) ZCS_CH3COCH3=cs_ch3coch3_298                     

  IF(PT(JK)<=235.) THEN
    ZCS_CH3COCH3=CS_CH3COCH3_235
  ELSEIF(PT(JK)>235. .AND. PT(JK)<=254.) THEN
    ZTSCALE=(PT(JK)-235._JPRB)/19.0_JPRB
    ZCS_CH3COCH3=((1.0_JPRB-ZTSCALE)*CS_CH3COCH3_235)+(ZTSCALE*CS_CH3COCH3_254)
  ELSEIF(PT(JK)>254. .AND. PT(JK)<=263.) THEN
    ZTSCALE=(PT(JK)-254._JPRB)/9.0_JPRB
    ZCS_CH3COCH3=((1.0_JPRB-ZTSCALE)*CS_CH3COCH3_254)+(ZTSCALE*CS_CH3COCH3_263)
  ELSEIF(PT(JK)>263. .AND. PT(JK)<=280.) THEN
    ZTSCALE=(PT(JK)-264._JPRB)/16.0_JPRB
    ZCS_CH3COCH3=((1.0_JPRB-ZTSCALE)*CS_CH3COCH3_263)+(ZTSCALE*CS_CH3COCH3_280)
  ELSEIF(PT(JK)>280. .AND. PT(JL)<= 298.) THEN
    ZTSCALE=(PT(JK)-280._JPRB)/18.0_JPRB
    ZCS_CH3COCH3=((1.0_JPRB-ZTSCALE)*CS_CH3COCH3_280)+(ZTSCALE*CS_CH3COCH3_298)
  ELSE 
    ZCS_CH3COCH3=CS_CH3COCH3_298
  ENDIF

!==============================================================================================
! re-initialize the temporary arrays for the next layer
!=============================================================================================== 

  ZDELTA1 = 0. ; ZDELTA2 = 0. ;ZDELTA3 = 0. ; ZDELTA4 = 0. ; ZDELTA5 = 0. ; ZDELTA6 = 0. ; ZDELTA7 = 0.
  
!   ============================================================================================
!       Tropo band 1
!    ===========================================================================================
! temp indep species for this band: no2,h2o2,n2o5,hno4,ch2o,ch3ooh,ald2,orgntr                       
  DO JL = IBAND_START(1),IBAND_END(1) 
    ZBJO1D(1)    = ZBJO1D(1) + PCST_O3_COL(JK,JL)*PQY_O1D_COL(JK,JL)*PFDIR(JK,JL)
    ZBJNO2(1)    = ZBJNO2(1) + PCST_NO2_COL(JK,JL)*PQY_NO2_COL(JK,JL)*PFDIR(JK,JL)
    ZBJH2O2(1)   = ZBJH2O2(1) + CS_H2O2(JL)*PFDIR(JK,JL)
    ZBJHONO(1)   = ZBJHONO(1) + CS_HONO(JL)*PFDIR(JK,JL)    
    ZBJHNO3(1)   = ZBJHNO3(1) + PCST_HNO3_COL(JK,JL)*PFDIR(JK,JL)
    ZBJHNO4(1)   = ZBJHNO4(1) + CS_HNO4(JL)*PFDIR(JK,JL)        
    ZBJN2O5(1)   = ZBJN2O5(1) + 0.8*CS_N2O5(JL)*PFDIR(JK,JL)
    ZBJCHOH(1)   = ZBJCHOH(1) + CS_CH2O(JL)*QYH_HCO(JL)*PFDIR(JK,JL)    
    ZBJCH3OOH(1) = ZBJCH3OOH(1) + CS_CH3OOH(JL)*PFDIR(JK,JL)
    ZBJALD2(1)   = ZBJALD2(1) + CS_ALD2(JL)*PFDIR(JK,JL)
    ZBJPANA(1)   = ZBJPANA(1) + PCST_PAN_COL(JK,JL)*QYPAN(JL)*PFDIR(JK,JL)
    ZBJPANB(1)   = ZBJPANB(1) + PCST_PAN_COL(JK,JL)*(1._JPRB-QYPAN(JL))*PFDIR(JK,JL)      
    ZBJCH3ONO2(1)= ZBJCH3ONO2(1) + CS_ORGN(JL) * PFDIR(JK,JL)
    ZBJGLYAL(1)  = ZBJGLYAL(1) + CS_GLYAL(JL) * PFDIR(JK,JL)    
    ZBJO2(1)     = ZBJO2(1) + CS_O2(JL) * PFDIR(JK,JL)
!    ZBJ_a_acet(1)= ZBJ_a_acet(1) + ZCS_CH3COCH3(JL)*PQY_CO_COL(JK,JL)*PFDIR(JK,JL)
  ENDDO
                             
! ===================================================================================== 
!           Tropo band 2                                                      
! =====================================================================================
! temp indep species for this band: hno4,ch3ooh,ald2,orgntr                             
  DO JL = IBAND_START(2),IBAND_END(2) ! temp indep for : no2,hno4,ch2o,ch3ooh,orgntr,ispd
    ZBJO1D(2)    = ZBJO1D(2) + PCST_O3_COL(JK,JL)*PQY_O1D_COL(JK,JL)*PFDIR(JK,JL)
    ZBJNO2(2)    = ZBJNO2(2) + PCST_NO2_COL(JK,JL)*PQY_NO2_COL(JK,JL)*PFDIR(JK,JL)
    ZBJH2O2(2)   = ZBJH2O2(2) + PCST_H2O2_COL(JK,JL-17)*PFDIR(JK,JL)
    ZBJHONO(2)   = ZBJHONO(2) + CS_HONO(JL)*PFDIR(JK,JL)
    ZBJHNO3(2)   = ZBJHNO3(2) + PCST_HNO3_COL(JK,JL)*PFDIR(JK,JL)
    ZBJHNO4(2)   = ZBJHNO4(2) + CS_HNO4(JL)*PFDIR(JK,JL)
    ZBJN2O5(2)   = ZBJN2O5(2) + 0.8*PCST_N2O5_COL(JK,JL-17)*PFDIR(JK,JL)
    ZBJCOH2(2)   = ZBJCOH2(2) + PCST_CH2O_COL(JK,JL-17)*QYH2_CO(JL)*PFDIR(JK,JL)    
    ZBJCHOH(2)   = ZBJCHOH(2) + PCST_CH2O_COL(JK,JL-17)*QYH_HCO(JL)*PFDIR(JK,JL)   
    ZBJCH3OOH(2) = ZBJCH3OOH(2) + CS_CH3OOH(JL)*PFDIR(JK,JL)
    ZBJPANA(2)   = ZBJPANA(2) + PCST_PAN_COL(JK,JL)*QYPAN(JL)*PFDIR(JK,JL)
    ZBJPANB(2)   = ZBJPANB(2) + PCST_PAN_COL(JK,JL)*(1._JPRB-QYPAN(JL))*PFDIR(JK,JL)  
    ZBJALD2(2)   = ZBJALD2(2) + CS_ALD2(JL)* PFDIR(JK,JL)
    ZBJCH3COCHO(2)= ZBJCH3COCHO(2) + CS_CH3COCHO(JL)*PFDIR(JK,JL)
    ZBJCH3ONO2(2) = ZBJCH3ONO2(2) + CS_ORGN(JL) * PFDIR(JK,JL)
    ZBJISPD(2)    = ZBJISPD(2) + CS_ISPD(JL) * QY_ISPD(JL) * PFDIR(JK,JL)
    ZBJ_A_ACET(2) = ZBJ_A_ACET(2) + ZCS_CH3COCH3(JL)*PQY_CO_COL(JK,JL)*PFDIR(JK,JL)
    ZBJ_B_ACET(2) = ZBJ_B_ACET(2) + ZCS_CH3COCH3(JL)*PQY_C2O3_COL(JK,JL)*PFDIR(JK,JL)
    ZBJGLYA(2)    = ZBJGLYA(2) + CS_GLY(JL-17)*(QY_GLY_TOT(JL-17)-QY_GLY(JL-17))*PFDIR(JK,JL)
    ZBJGLYB(2)    = ZBJGLYB(2) + CS_GLY(JL-17)*QY_GLY(JL-17)*PFDIR(JK,JL)    
    ZBJGLYAL(2)   = ZBJGLYAL(2) + CS_GLYAL(JL) * PFDIR(JK,JL)
    ZBJHALD(2)    = ZBJHALD(2) + CS_HALD(JL) * PFDIR(JK,JL)            
  ENDDO 
    
! ======================================================================================= 
!           Tropo band 3                                                      
! ======================================================================================= 
! temp indep species for this band: hno4,ch3ooh,ald2,mgly,orgntr,ispd      
  DO JL = IBAND_START(3),IBAND_END(3)
    ZBJO1D(3)   = ZBJO1D(3) + PCST_O3_COL(JK,JL)*PQY_O1D_COL(JK,JL)*PFDIR(JK,JL)
    ZBJNO2(3)   = ZBJNO2(3) + PCST_NO2_COL(JK,JL)*PQY_NO2_COL(JK,JL)*PFDIR(JK,JL)
    ZBJH2O2(3)   = ZBJH2O2(3) + PCST_H2O2_COL(JK,JL-17)*PFDIR(JK,JL)
    ZBJHONO(3)   = ZBJHONO(3) + CS_HONO(JL)*PFDIR(JK,JL)    
    ZBJHNO3(3)   = ZBJHNO3(3) + PCST_HNO3_COL(JK,JL)*PFDIR(JK,JL)
    ZBJHNO4(3)   = ZBJHNO4(3) + CS_HNO4(JL)*PFDIR(JK,JL)
    ZBJN2O5(3)   = ZBJN2O5(3) + 0.8*PCST_N2O5_COL(JK,JL-17)*PFDIR(JK,JL)
    ZBJCOH2(3)   = ZBJCOH2(3) + PCST_CH2O_COL(JK,JL-17)*QYH2_CO(JL)*PFDIR(JK,JL)   
    ZBJCHOH(3)   = ZBJCHOH(3) + PCST_CH2O_COL(JK,JL-17)*QYH_HCO(JL)*PFDIR(JK,JL)
    ZBJCH3OOH(3) = ZBJCH3OOH(3) + CS_CH3OOH(JL)*PFDIR(JK,JL)
    ZBJPANA(3)   = ZBJPANA(3) + PCST_PAN_COL(JK,JL)*QYPAN(JL)*PFDIR(JK,JL)
    ZBJPANB(3)   = ZBJPANB(3) + PCST_PAN_COL(JK,JL)*(1._JPRB-QYPAN(JL))*PFDIR(JK,JL)  
    ZBJALD2(3)   = ZBJALD2(3) + CS_ALD2(JL)*PFDIR(JK,JL)  
    ZBJCH3COCHO(3) = ZBJCH3COCHO(3) + CS_CH3COCHO(JL)*PFDIR(JK,JL)
    ZBJCH3ONO2(3)  = ZBJCH3ONO2(3) + CS_ORGN(JL)*PFDIR(JK,JL)
    ZBJISPD(3)     = ZBJISPD(3) + CS_ISPD(JL) *QY_ISPD(JL) * PFDIR(JK,JL)
    ZBJ_A_ACET(3)  = ZBJ_A_ACET(3) + ZCS_CH3COCH3(JL)*PQY_CO_COL(JK,JL)*PFDIR(JK,JL)
    ZBJ_B_ACET(3)  = ZBJ_B_ACET(3) + ZCS_CH3COCH3(JL)*PQY_C2O3_COL(JK,JL)*PFDIR(JK,JL)
    ZBJGLYA(3)     = ZBJGLYA(3) + CS_GLY(JL-17)*(QY_GLY_TOT(JL-17)-QY_GLY(JL-17))*PFDIR(JK,JL)
    ZBJGLYB(3)     = ZBJGLYB(3) + CS_GLY(JL-17)*QY_GLY(JL-17)*PFDIR(JK,JL)       
    ZBJGLYAL(3)    = ZBJGLYAL(3) + CS_GLYAL(JL) * PFDIR(JK,JL)
    ZBJHALD(3)     = ZBJHALD(3) + CS_HALD(JL) * PFDIR(JK,JL)                  
  ENDDO
 
! ================================================================================ 
!           Tropo band 4                                                      
! ================================================================================
! temp indep species for this band: hno4,ch3ooh,ald2,orgntr,ispd  
  DO JL = IBAND_START(4),IBAND_END(4)
    ZBJO1D(4)  = ZBJO1D(4)  + PCST_O3_COL(JK,JL)*PQY_O1D_COL(JK,JL)*PFDIR(JK,JL)
    ZBJNO2(4)  = ZBJNO2(4)  + PCST_NO2_COL(JK,JL)*PQY_NO2_COL(JK,JL)*PFDIR(JK,JL)
    ZBJH2O2(4) = ZBJH2O2(4) + PCST_H2O2_COL(JK,JL-17)*PFDIR(JK,JL)
    ZBJHONO(4) = ZBJHONO(4) + CS_HONO(JL)*PFDIR(JK,JL)     
    ZBJHNO3(4) = ZBJHNO3(4) + PCST_HNO3_COL(JK,JL)*PFDIR(JK,JL)
    ZBJHNO4(4) = ZBJHNO4(4) + CS_HNO4(JL)*PFDIR(JK,JL)
    ZBJN2O5(4) = ZBJN2O5(4) + PCST_N2O5_COL(JK,JL-17)*PFDIR(JK,JL)
    ZBJCOH2(4) = ZBJCOH2(4) + PCST_CH2O_COL(JK,JL-17)*QYH2_CO(JL)*PFDIR(JK,JL)
    ZBJCHOH(4) = ZBJCHOH(4) + PCST_CH2O_COL(JK,JL-17)*QYH_HCO(JL)*PFDIR(JK,JL)   
    ZBJCH3OOH(4)= ZBJCH3OOH(4) + CS_CH3OOH(JL)*PFDIR(JK,JL)
    ZBJPANA(4)  = ZBJPANA(4) + PCST_PAN_COL(JK,JL)*QYPAN(JL)*PFDIR(JK,JL)
    ZBJPANB(4)  = ZBJPANB(4) + PCST_PAN_COL(JK,JL)*(1.0_JPRB-QYPAN(JL))*PFDIR(JK,JL) 
    ZBJALD2(4)  = ZBJALD2(4) + CS_ALD2(JL)*PFDIR(JK,JL)
    ZBJCH3COCHO(4) = ZBJCH3COCHO(4) + CS_CH3COCHO(JL)*PQY_CH3COCHO_COL(JK,JL-37)*PFDIR(JK,JL)
    ZBJCH3ONO2(4)  = ZBJCH3ONO2(4) + CS_ORGN(JL)*PFDIR(JK,JL)
    ZBJISPD(4)     = ZBJISPD(4) + CS_ISPD(JL) * QY_ISPD(JL) *PFDIR(JK,JL)
    ZBJ_A_ACET(4)  = ZBJ_A_ACET(4) + ZCS_CH3COCH3(JL)*PQY_CO_COL(JK,JL)*PFDIR(JK,JL)
    ZBJ_B_ACET(4)  = ZBJ_B_ACET(4) + ZCS_CH3COCH3(JL)*PQY_C2O3_COL(JK,JL)*PFDIR(JK,JL)
    ZBJGLYA(4)     = ZBJGLYA(4) + CS_GLY(JL-17)*(QY_GLY_TOT(JL-17)-QY_GLY(JL-17))*PFDIR(JK,JL)
    ZBJGLYB(4)     = ZBJGLYB(4) + CS_GLY(JL-17)*QY_GLY(JL-17)*PFDIR(JK,JL)       
    ZBJGLYAL(4)    = ZBJGLYAL(4) + CS_GLYAL(JL) * PFDIR(JK,JL)
    ZBJHALD(4)     = ZBJHALD(4) + CS_HALD(JL) * PFDIR(JK,JL)             
  ENDDO 
                            
! ======================================================================================== 
!           Tropo band 5                                                      
! ======================================================================================== 
! temp indep species for this band: hno4,ch3ooh,ald2,orgntr,ispd                             
  DO JL = IBAND_START(5),IBAND_END(5)
    ZBJO1D(5)      = ZBJO1D(5)+PCST_O3_COL(JK,JL)*PQY_O1D_COL(JK,JL)*PFDIR(JK,JL)
    ZBJNO2(5)      = ZBJNO2(5)+PCST_NO2_COL(JK,JL)*PQY_NO2_COL(JK,JL)*PFDIR(JK,JL)
    ZBJH2O2(5)     = ZBJH2O2(5)+PCST_H2O2_COL(JK,JL-17)*PFDIR(JK,JL)
    ZBJHONO(5)     = ZBJHONO(5)+CS_HONO(JL)*PFDIR(JK,JL)        
    ZBJHNO3(5)     = ZBJHNO3(5)+PCST_HNO3_COL(JK,JL)*PFDIR(JK,JL)
    ZBJHNO4(5)     = ZBJHNO4(5)+CS_HNO4(JL)*PFDIR(JK,JL)
    ZBJN2O5(5)     = ZBJN2O5(5)+PCST_N2O5_COL(JK,JL-17)*PFDIR(JK,JL)
    ZBJCOH2(5)     = ZBJCOH2(5)+PCST_CH2O_COL(JK,JL-17)*QYH2_CO(JL)*PFDIR(JK,JL)
    ZBJCHOH(5)     = ZBJCHOH(5)+PCST_CH2O_COL(JK,JL-17)*QYH_HCO(JL)*PFDIR(JK,JL)
    ZBJCH3OOH(5)   = ZBJCH3OOH(5)+CS_CH3OOH(JL)*PFDIR(JK,JL)
    ZBJPANA(5)     = ZBJPANA(5) + PCST_PAN_COL(JK,JL)*QYPAN(JL)*PFDIR(JK,JL)
    ZBJPANB(5)     = ZBJPANB(5) + PCST_PAN_COL(JK,JL)*(1.0_JPRB-QYPAN(JL))*PFDIR(JK,JL) 
    ZBJALD2(5)     = ZBJALD2(5)+CS_ALD2(JL)* PFDIR(JK,JL)
    ZBJCH3COCHO(5) = ZBJCH3COCHO(5)+CS_CH3COCHO(JL)*PQY_CH3COCHO_COL(JK,JL-37)*PFDIR(JK,JL)
    ZBJCH3ONO2(5)  = ZBJCH3ONO2(5)+CS_ORGN(JL)*PFDIR(JK,JL)
    ZBJISPD(5)     = ZBJISPD(5) + CS_ISPD(JL) *QY_ISPD(JL) * PFDIR(JK,JL)
    ZBJGLYA(5)     = ZBJGLYA(5) + CS_GLY(JL-17)*(QY_GLY_TOT(JL-17)-QY_GLY(JL-17))*PFDIR(JK,JL)
    ZBJGLYB(5)     = ZBJGLYB(5) + CS_GLY(JL-17)*QY_GLY(JL-17)*PFDIR(JK,JL)     
    ZBJGLYAL(5)    = ZBJGLYAL(5) + CS_GLYAL(JL) * PFDIR(JK,JL)
    ZBJHALD(5)     = ZBJHALD(5) + CS_HALD(JL) * PFDIR(JK,JL)              
  ENDDO 
!=========================================================================================
!           Tropo band 6
!=========================================================================================
 
  DO JL = IBAND_START(6),IBAND_END(6)    
    ZBJNO2(6)     = ZBJNO2(6) +PCST_NO2_COL(JK,JL)*PQY_NO2_COL(JK,JL)*PFDIR(JK,JL)
    ZBJCOH2(6)     = ZBJCOH2(6)+PCST_CH2O_COL(JK,JL-17)*QYH2_CO(JL)*PFDIR(JK,JL)
    ZBJCH3OOH(6)   = ZBJCH3OOH(6)+CS_CH3OOH(JL)*PFDIR(JK,JL)
    ZBJNO2O(6)     = ZBJNO2O(6)+PCST_NO3_COL(JK,JL-60)*QYNO2_O(JL-60)*PFDIR(JK,JL)
    ZBJCH3COCHO(6) = ZBJCH3COCHO(6)+CS_CH3COCHO(JL)*PQY_CH3COCHO_COL(JK,JL-37)*PFDIR(JK,JL)
    ZBJISPD(6)     = ZBJISPD(6) + CS_ISPD(JL) *QY_ISPD(JL) * PFDIR(JK,JL)
    IF(JL < 100) then
     ZBJGLYA(6)    = ZBJGLYA(6) + CS_GLY(JL-17)*(QY_GLY_TOT(JL-17)-QY_GLY(JL-17))*PFDIR(JK,JL)
     ZBJGLYB(6)    = ZBJGLYB(6) + CS_GLY(JL-17)*QY_GLY(JL-17)*PFDIR(JK,JL)    
    ENDIF
  ENDDO
     
! Add the contribution from band 6 separately for H2O2, N2O5 and HNO3 for high angle grid to allow reducion in cst arrays
! This has been tested using a box model to ensure that no errors are introduced compared to using main array

  IF(PANGLE>71.) THEN

    ZBJH2O2(6)     = ZBJH2O2(6)+PCST_H2O2_COL(JK,44)*PFDIR(JK,61)
    ZBJH2O2(6)     = ZBJH2O2(6)+PCST_H2O2_COL(JK,45)*PFDIR(JK,62)
  
    DO JC=63,68
      ZBJN2O5(6)   = ZBJN2O5(6)+PCST_N2O5_COL(JK,JC-17)*PFDIR(JK,JC)
    ENDDO
       
    ZBJHNO3(6)     = ZBJHNO3(6)+PCST_HNO3_COL(JK,63)*PFDIR(JK,63)
    ZBJPANA(6)      = 0._JPRB
    ZBJPANB(6)      = 0._JPRB    
  ELSEIF(PANGLE<=71.) THEN

    DO JC=61,63
      ZBJHNO3(6)     = ZBJHNO3(6)+PCST_HNO3_COL(JK,JC)*PFDIR(JK,JC)
    ENDDO

    DO JC=61,69
      ZBJN2O5(6)     = ZBJN2O5(6)+PCST_N2O5_COL(JK,JC-17)*PFDIR(JK,JC)
    ENDDO
       
    DO JC=61,62
      ZBJPANA(6)      = ZBJPANA(6)+PCST_PAN_COL(JK,JC)*QYPAN(JC)*PFDIR(JK,JC)
      ZBJPANB(6)      = ZBJPANB(6)+PCST_PAN_COL(JK,JC)*(1._JPRB-QYPAN(JC))*PFDIR(JK,JC)      
    ENDDO 

  ENDIF           
!======================================================================================== 
!           Tropo band 7                                                      
!========================================================================================= 
           
  DO JL = IBAND_START(7),IBAND_END(7)-2
    ZBJNO2O(7)     = ZBJNO2O(7)+PCST_NO3_COL(JK,JL-60)*QYNO2_O(JL-60)*PFDIR(JK,JL)
    ZBJNOO2(7)     = ZBJNOO2(7)+PCST_NO3_COL(JK,JL-60)*QYNO_O2(JL-60)*PFDIR(JK,JL)
  ENDDO  

 ! only set a scaling ratio for the first band once the sza > 71.
 ! only calculate limits to bands 2 and 4 for high sza         
             
  IF(PANGLE<=71.) THEN
    ZDELTA1 = PFACT(JK,1)/PFDIR(JK,IBAND_MIDDLE(1))
    ZDELTA2 = PFACT(JK,2)/PFDIR(JK,IBAND_MIDDLE(2))
    ZDELTA3 = PFACT(JK,3)/PFDIR(JK,IBAND_MIDDLE(3))
    ZDELTA4 = PFACT(JK,4)/PFDIR(JK,IBAND_MIDDLE(4))
    ZDELTA5 = PFACT(JK,5)/PFDIR(JK,IBAND_MIDDLE(5))
    ZDELTA6 = PFACT(JK,6)/PFDIR(JK,IBAND_MIDDLE(6))
    ZDELTA7 = PFACT(JK,7)/PFDIR(JK,IBAND_MIDDLE(7))
  ELSEIF(PANGLE>71. .AND. PANGLE<=SZA_LIMIT) THEN
    ZDELTA1 = PFACT(JK,1)/PFDIR(JK,IBAND_MIDDLE(1))
    ZDELTA2 = PFACT(JK,2)/PFDIR(JK,IBAND_MIDDLE(2))
    ZDELTA3 = PFACT(JK,3)/PFDIR(JK,IBAND_MIDDLE(3))
    ZDELTA4 = PFACT(JK,4)/PFDIR(JK,IBAND_MIDDLE(4))
    ZDELTA5 = PFACT(JK,5)/PFDIR(JK,IBAND_MIDDLE(5))
    ZDELTA6 = PFACT(JK,6)/PFDIR(JK,IBAND_MIDDLE(6))
    ZDELTA7 = PFACT(JK,7)/PFDIR(JK,IBAND_MIDDLE(7))
    IF(PFDIR(JK,IBAND_MIDDLE(1))>5.E9) ZDELTA1_SPEC = PFACT(JK,1)/PFDIR(JK,IBAND_MIDDLE(1))      
    IF(PFDIR(JK,IBAND_MIDDLE(3))>1.E9) ZDELTA3_SPEC = PFACT(JK,3)/PFDIR(JK,IBAND_MIDDLE(3))
  ENDIF  

  IF(PANGLE<=80._JPRB) THEN
                              
  !  ============================================================
  !  CALCULATION PHOTOLYSIS RATES FOR SCATTERING ATMOSPHERE    
  !  ============================================================

  ! O2 -> O3P
  PRJ(JK,JO2)   = ZBJO2(1)*ZDELTA1
   
  ! O3 -> 2OH 
  PRJ(JK,JO3D)  =  ZBJO1D(1)*ZDELTA1 + ZBJO1D(2)*ZDELTA2 &
    &  + ZBJO1D(3)*ZDELTA3 + ZBJO1D(4)*ZDELTA4 &
    &  + ZBJO1D(5)*ZDELTA5      

  ! NO2 -> O3
  PRJ(JK,JNO2)  = ZBJNO2(1)*ZDELTA1 + ZBJNO2(2)*ZDELTA2 + &
    &   ZBJNO2(3)*ZDELTA3 + ZBJNO2(4)*ZDELTA4 + &
    &   ZBJNO2(5)*ZDELTA5 + ZBJNO2(6)*ZDELTA6   
  
  ! HNO3 -> OH + NO2 
  PRJ(JK,JHNO3)  = ZBJHNO3(1)*ZDELTA1 + ZBJHNO3(2)*ZDELTA2 + &
      &    ZBJHNO3(3)*ZDELTA3 + ZBJHNO3(4)*ZDELTA4 + &
      &    ZBJHNO3(5)*ZDELTA5 + ZBJHNO3(6)*ZDELTA6

  ! HNO4 -> HO2 + NO2 
  PRJ(JK,JHNO4)  =  ZBJHNO4(1)*ZDELTA1 + ZBJHNO4(2)*ZDELTA2 + &
      &     ZBJHNO4(3)*ZDELTA3 + ZBJHNO4(4)*ZDELTA4 + &
      &     ZBJHNO4(5)*ZDELTA5 
      
  ! HONO -> OH + NO 
  PRJ(JK,JHONO)  =  ZBJHONO(1)*ZDELTA1 + ZBJHONO(2)*ZDELTA2 + &
      &     ZBJHONO(3)*ZDELTA3 + ZBJHONO(4)*ZDELTA4 + &
      &     ZBJHONO(5)*ZDELTA5       
      
  ! CH3O2NO2 -> 0.5*CH3O2 + 0.5*NO2 + 0.5*CH2O + 0.5*HO2 + 0.5*NO3 
  PRJ(JK,JMENO2) = PRJ(JK,JHNO4)
  
  ! N2O5 -> NO2 + NO2               
  PRJ(JK,JN2O5)  =  ZBJN2O5(1)*ZDELTA1 + ZBJN2O5(2)*ZDELTA2 +  &
      &     ZBJN2O5(3)*ZDELTA3 + ZBJN2O5(4)*ZDELTA4 +  &
      &     ZBJN2O5(5)*ZDELTA5 + ZBJN2O5(6)*ZDELTA6 

  ! PAN -> NO2 + C2O3
  PRJ(JK,JPANA) =   ZBJPANA(1)*ZDELTA1  + ZBJPANA(2)*ZDELTA2 +   &
      &    ZBJPANA(3)*ZDELTA3  + ZBJPANA(4)*ZDELTA4 +   &
      &    ZBJPANA(5)*ZDELTA5  + ZBJPANA(6)*ZDELTA6
      
  ! PAN -> NO3 + CH3O2
  PRJ(JK,JPANB) =   ZBJPANB(1)*ZDELTA1  + ZBJPANB(2)*ZDELTA2 +   &
      &    ZBJPANB(3)*ZDELTA3  + ZBJPANB(4)*ZDELTA4 +   &
      &    ZBJPANB(5)*ZDELTA5  + ZBJPANB(6)*ZDELTA6      

  ! CH2O -> CO
  PRJ(JK,JACH2O)  =  ZBJCOH2(2)*ZDELTA2  +  ZBJCOH2(3)*ZDELTA3 &
      &    + ZBJCOH2(4)*ZDELTA4 +  ZBJCOH2(5)*ZDELTA5 &
      &    + ZBJCOH2(6)*ZDELTA6   

  ! CH2O -> COH + H
  PRJ(JK,JBCH2O)  =  ZBJCHOH(1)*ZDELTA1 + ZBJCHOH(2)*ZDELTA2 &
      &    + ZBJCHOH(3)*ZDELTA3 + ZBJCHOH(4)*ZDELTA4 &
      &    + ZBJCHOH(5)*ZDELTA5 + ZBJCHOH(6)*ZDELTA6

  ! H2O2 -> 2OH
  PRJ(JK,JH2O2) = ZBJH2O2(1)*ZDELTA1 + ZBJH2O2(2)*ZDELTA2 &
      &    + ZBJH2O2(3)*ZDELTA3 + ZBJH2O2(4)*ZDELTA4 &
      &    + ZBJH2O2(5)*ZDELTA5 + ZBJH2O2(6)*ZDELTA6 

  ! NO3 -> NO2 + O
  PRJ(JK,JANO3)  = ZBJNO2O(6)*ZDELTA6 + ZBJNO2O(7)*ZDELTA7      

  ! NO3 -> NO
  PRJ(JK,JBNO3) = ZBJNOO2(7)*ZDELTA7

  ! CH3OOH -> HCHO + HO2 + OH     
  PRJ(JK,JMEPE) = ZBJCH3OOH(1)*ZDELTA1 + ZBJCH3OOH(2)*ZDELTA2 + &
      &    ZBJCH3OOH(3)*ZDELTA3 + ZBJCH3OOH(4)*ZDELTA4 + &
      &    ZBJCH3OOH(5)*ZDELTA5 + ZBJCH3OOH(6)*ZDELTA6 

  ! ALD2 -> HCHO + XO2 + CO + 2HO2 
  PRJ(JK,J45) = ZBJALD2(1)*ZDELTA1 + ZBJALD2(2)*ZDELTA2 &
      &      + ZBJALD2(3)*ZDELTA3 + ZBJALD2(4)*ZDELTA4 &
      &      + ZBJALD2(5)*ZDELTA5

  ! CH3ONO2 -> HO2 + NO2
  PRJ(JK,JORGN) = ZBJCH3ONO2(1)*ZDELTA1 + ZBJCH3ONO2(2)*ZDELTA2 + &
      &    ZBJCH3ONO2(3)*ZDELTA3 + ZBJCH3ONO2(4)*ZDELTA4 + &
      &    ZBJCH3ONO2(5)*ZDELTA5

  ! CH3COCHO -> C2O3 + HO2 + CO       
  PRJ(JK,J74) = ZBJCH3COCHO(2)*ZDELTA2 + ZBJCH3COCHO(3)*ZDELTA3 +  &
       &        ZBJCH3COCHO(4)*ZDELTA4 + ZBJCH3COCHO(5)*ZDELTA5 + &
       &        ZBJCH3COCHO(6)*ZDELTA6 + ZBJCH3COCHO(7)*ZDELTA7

  ! ispd -> 0.333CO + 0.067ALD2 + 0.9CH2O + 0.832PAR + 1.033HO2 + 0.7XO2 + 0.967C2O3
  PRJ(JK,JISPD) = ZBJISPD(2)*ZDELTA2 + ZBJISPD(3)*ZDELTA3 + &
       &          ZBJISPD(4)*ZDELTA4 + ZBJISPD(5)*ZDELTA5 + &
       &          ZBJISPD(6)*ZDELTA6

   !ch3coch3 -> 2CH3O2 + CO
  PRJ(JK,JA_ACET) =  ZBJ_A_ACET(2)*ZDELTA2 &
       &          +  ZBJ_A_ACET(3)*ZDELTA3 + ZBJ_A_ACET(4)*ZDELTA4
        
  PRJ(JK,JB_ACET) =  ZBJ_B_ACET(2)*ZDELTA2 & 
       &          +  ZBJ_B_ACET(3)*ZDELTA3 + ZBJ_B_ACET(4)*ZDELTA4 
       
   ! ch2ohcho -> 2HO2 + CO + HCHO
  PRJ(JK,JGLYALD) = ZBJGLYAL(1)*ZDELTA1 + ZBJGLYAL(2)*ZDELTA2 + &
       &            ZBJGLYAL(3)*ZDELTA3 + ZBJGLYAL(4)*ZDELTA4 + &
       &            ZBJGLYAL(5)*ZDELTA5 
          
   ! HPALD2 -> OH + HO2 + 0.5HYAC + 0.5MGLY + 0.5GLY + HCHO
  PRJ(JK,JHPALDA) = ZBJHALD(2)*ZDELTA2 + ZBJHALD(3)*ZDELTA3 + &
      &             ZBJHALD(4)*ZDELTA4 + ZBJHALD(5)*ZDELTA5
   ! GLY -> 2CO + 2HO2   
  PRJ(JK,JGLYA) = ZBJGLYA(2)*ZDELTA2 + ZBJGLYA(3)*ZDELTA3 + &
      &           ZBJGLYA(4)*ZDELTA4 + ZBJGLYA(5)*ZDELTA5 + &
      &           ZBJGLYA(6)*ZDELTA6
      
   ! GLY -> HCHO + CO 
  PRJ(JK,JGLYB) = ZBJGLYB(2)*ZDELTA2 + ZBJGLYB(3)*ZDELTA3 + &
      &           ZBJGLYB(4)*ZDELTA4 + ZBJGLYB(5)*ZDELTA5 + &
      &           ZBJGLYB(6)*ZDELTA6            
          
  DO JJ = 1,NPHOTO
     IF (PRJ(JK,JJ) < 0.0_JPRB) THEN
       WRITE(NULOUT,'(a30,2i5,e16.8)') 'NEG. PHOT. RATE - low angle',JK,JJ,PRJ(JK,JJ)
       WRITE(NULOUT,*)'ZDELTA1',ZDELTA1
       CALL ABOR1('tm5_photorates_tropo: error with photolysis')
     ENDIF
  ENDDO

  ELSEIF(PANGLE>80. .AND. PANGLE<=85.) THEN
   ! PRJ(JK,jo3d) = 0. ; PRJ(JK,jhno3) = 0. ; PRJ(JK,jhno4) = 0. ; PRJ(JK,jn2o5) = 0.
   ! PRJ(JK,jh2o2) = 0. ; PRJ(JK,jach2o) = 0. ; PRJ(JK,jbch2o) = 0.
   ! PRJ(JK,jpan) = 0. ; PRJ(JK,j45) = 0. ; PRJ(JK,jorgn) = 0. ; PRJ(JK,j74) = 0.
   ! PRJ(JK,jo2) = 0. ; PRJ(JK,jispd) = 0. ; PRJ(JK,ja_acet) = 0. ; PRJ(JK,JGLYAL) = 0.
                                
    ! O2 -> 2O3P
    PRJ(JK,JO2) = ZBJO2(1)*ZDELTA1_SPEC

    ! O3 -> O(1D) + O2            
    PRJ(JK,JO3D)  =  ZBJO1D(1)*ZDELTA1_SPEC + ZBJO1D(2)*ZDELTA2 +   &
       &            ZBJO1D(3)*ZDELTA3_SPEC + ZBJO1D(4)*ZDELTA4 +   &
       &            ZBJO1D(5)*ZDELTA5        
        
    ! HNO3 -> OH + NO2                
    PRJ(JK,JHNO3)  = ZBJHNO3(1)*ZDELTA1_SPEC + ZBJHNO3(2)*ZDELTA2 + &
       &            ZBJHNO3(3)*ZDELTA3_SPEC + ZBJHNO3(4)*ZDELTA4 + &
       &            ZBJHNO3(5)*ZDELTA5 + ZBJHNO3(6)*ZDELTA6      
    
    ! N2O5 -> NO2 + NO2      
    PRJ(JK,JN2O5)  =  ZBJN2O5(1)*ZDELTA1_SPEC + ZBJN2O5(2)*ZDELTA2 +  &
       &             ZBJN2O5(3)*ZDELTA3_SPEC + ZBJN2O5(4)*ZDELTA4 +  &
       &             ZBJN2O5(5)*ZDELTA5  + ZBJN2O5(6)*ZDELTA6 

    ! PAN -> C3O2 + NO2
    PRJ(JK,JPANA) = ZBJPANA(1)*ZDELTA1_SPEC + ZBJPANA(2)*ZDELTA2 +   &
       &          ZBJPANA(3)*ZDELTA3_SPEC + ZBJPANA(4)*ZDELTA4 +   &
       &          ZBJPANA(5)*ZDELTA5  + ZBJPANA(6)*ZDELTA6 
       
     ! PAN -> CH3O2 + NO3
    PRJ(JK,JPANB) = ZBJPANB(1)*ZDELTA1_SPEC + ZBJPANB(2)*ZDELTA2 +   &
       &          ZBJPANB(3)*ZDELTA3_SPEC + ZBJPANB(4)*ZDELTA4 +   &
       &          ZBJPANB(5)*ZDELTA5  + ZBJPANB(6)*ZDELTA6        
       
    ! ORGNTR -> HO2 + NO2
    PRJ(JK,JORGN) = ZBJCH3ONO2(1)*ZDELTA1_SPEC + ZBJCH3ONO2(2)*ZDELTA2 +   &
       &           ZBJCH3ONO2(3)*ZDELTA3_SPEC + ZBJCH3ONO2(4)*ZDELTA4 +   &
       &           ZBJCH3ONO2(5)*ZDELTA5 
       
    ! ALD2 -> HCHO + XO2 + CO + 2HO2
    PRJ(JK,J45) = ZBJALD2(1)*ZDELTA1_SPEC + ZBJALD2(2)*ZDELTA2 + &
       &         ZBJALD2(3)*ZDELTA3_SPEC + ZBJALD2(4)*ZDELTA4 + &
       &         ZBJALD2(5)*ZDELTA5
       
    ! NO2 -> O3
    PRJ(JK,JNO2)  =  ZBJNO2(1)*ZDELTA1_SPEC + ZBJNO2(2)*ZDELTA2 + &
       &            ZBJNO2(3)*ZDELTA3_SPEC + ZBJNO2(4)*ZDELTA4 + &
       &            ZBJNO2(5)*ZDELTA5 + ZBJNO2(6)*ZDELTA6
                                
    ! CH2O -> CO
    PRJ(JK,JACH2O)  =  ZBJCOH2(2)*ZDELTA2 + ZBJCOH2(3)*ZDELTA3_SPEC + &
        &             ZBJCOH2(4)*ZDELTA4 + ZBJCOH2(5)*ZDELTA5 + &
        &             ZBJCOH2(6)*ZDELTA6

    ! CH2O -> COH + H
    PRJ(JK,JBCH2O)  =  ZBJCHOH(1)*ZDELTA1_SPEC + ZBJCHOH(2)*ZDELTA2 + &
        &             ZBJCHOH(3)*ZDELTA3_SPEC + ZBJCHOH(4)*ZDELTA4 + &
        &             ZBJCHOH(5)*ZDELTA5 + ZBJCHOH(6)*ZDELTA6
  
    PRJ(JK,JHNO4)  =  ZBJHNO4(1)*ZDELTA1_SPEC + ZBJHNO4(2)*ZDELTA2 + &
        &            ZBJHNO4(3)*ZDELTA3_SPEC + ZBJHNO4(4)*ZDELTA4 + &
        &            ZBJHNO4(5)*ZDELTA5 
    ! HONO -> OH + NO    
    PRJ(JK,JHONO)  =  ZBJHONO(1)*ZDELTA1_SPEC + ZBJHONO(2)*ZDELTA2 + &
        &            ZBJHONO(3)*ZDELTA3_SPEC + ZBJHONO(4)*ZDELTA4 + &
        &            ZBJHONO(5)*ZDELTA5                        
      
  ! CH3O2NO2 -> 0.5*CH3O2 + 0.5*NO2 + 0.5*CH2O + 0.5*HO2 + 0.5*NO3 
    PRJ(JK,JMENO2) = PRJ(JK,JHNO4) 
       
    ! H2O2 -> 2OH
    PRJ(JK,JH2O2) = ZBJH2O2(1)*ZDELTA1_SPEC + ZBJH2O2(2)*ZDELTA2 + &
        &          ZBJH2O2(3)*ZDELTA3_SPEC + ZBJH2O2(4)*ZDELTA4 + &
        &          ZBJH2O2(5)*ZDELTA5 + ZBJH2O2(6)*ZDELTA6

    ! CH3OOH -> HCHO + HO2 + OH              
    PRJ(JK,JMEPE) = ZBJCH3OOH(1)*ZDELTA1_SPEC + ZBJCH3OOH(2)*ZDELTA2 + &
        &          ZBJCH3OOH(3)*ZDELTA3_SPEC + ZBJCH3OOH(4)*ZDELTA4 + &
        &          ZBJCH3OOH(5)*ZDELTA5 + ZBJCH3OOH(6)*ZDELTA6 

    ! CH3COCHO -> C2O3 + HO2 + CO        
    PRJ(JK,J74) = ZBJCH3COCHO(2)*ZDELTA2 + ZBJCH3COCHO(3)*ZDELTA3_SPEC +  &
        &        ZBJCH3COCHO(4)*ZDELTA4 + ZBJCH3COCHO(5)*ZDELTA5 +  &
        &        ZBJCH3COCHO(6)*ZDELTA6 + ZBJCH3COCHO(7)*ZDELTA7
    
    ! ispd -> 0.333CO + 0.067ALD2 + 0.9CH2O + 0.832PAR + 1.033HO2 + 0.7XO2 + 0.967C2O3
    PRJ(JK,JISPD) = ZBJISPD(2)*ZDELTA2 + ZBJISPD(3)*ZDELTA3_SPEC + &
        &          ZBJISPD(4)*ZDELTA4 + ZBJISPD(5)*ZDELTA5 + &
        &          ZBJISPD(6)*ZDELTA6
     
    PRJ(JK,JA_ACET) = ZBJ_A_ACET(2)*ZDELTA2 + &
        &            ZBJ_A_ACET(3)*ZDELTA3_SPEC + ZBJ_A_ACET(4)*ZDELTA4 
    
    PRJ(JK,JB_ACET) = ZBJ_B_ACET(2)*ZDELTA2 + &
        &            ZBJ_B_ACET(3)*ZDELTA3_SPEC + ZBJ_B_ACET(4)*ZDELTA4 
                       
   ! ch2ohcho -> 2HO2 + CO + HCHO
    PRJ(JK,JGLYALD) = ZBJGLYAL(1)*ZDELTA1_SPEC + ZBJGLYAL(2)*ZDELTA2 + &
        &             ZBJGLYAL(3)*ZDELTA3_SPEC + ZBJGLYAL(4)*ZDELTA4 + &
        &             ZBJGLYAL(5)*ZDELTA5 
        
   ! HPALD2 -> OH + HO2 + 0.5HYAC + 0.5MGLY + 0.5GLY + HCHO
    PRJ(JK,JHPALDA) = ZBJHALD(2)*ZDELTA2 + ZBJHALD(3)*ZDELTA3_SPEC + &
        &             ZBJHALD(4)*ZDELTA4 + ZBJHALD(5)*ZDELTA5 
        
   ! GLY -> 2CO + 2HO2   
     PRJ(JK,JGLYA) = ZBJGLYA(2)*ZDELTA2 + ZBJGLYA(3)*ZDELTA3_SPEC + &
      &             ZBJGLYA(4)*ZDELTA4 + ZBJGLYA(5)*ZDELTA5 + &
      &             ZBJGLYA(6)*ZDELTA6
      
   ! GLY -> HCHO + CO 
     PRJ(JK,JGLYB) = ZBJGLYB(2)*ZDELTA2 + ZBJGLYB(3)*ZDELTA3_SPEC + &
         &           ZBJGLYB(4)*ZDELTA4 + ZBJGLYB(5)*ZDELTA5 + &
         &           ZBJGLYB(6)*ZDELTA6                 
       
    DO JJ = 1,NPHOTO
     IF (PRJ(JK,JJ) < 0.0_JPRB) THEN
       WRITE(NULOUT,'(a30,2i5,e16.8)') 'NEG. PHOT. RATE - high angle',JK,JJ,PRJ(JK,JJ)
       WRITE(NULOUT,*)'ZDELTA1',ZDELTA1
       CALL ABOR1('tm5_photorates_tropo: error with photolysis')
     ENDIF
    ENDDO

          
  ENDIF ! sza > 80. 
  
  PRJ(JK,JROOH)   = PRJ(JK,JMEPE)
  PRJ(JK,JISOPOOH)= PRJ(JK,JMEPE)
  PRJ(JK,JHYAC)   = PRJ(JK,JA_ACET)+PRJ(JK,JB_ACET)
  !
  ! JEW (2020)
  ! Isoprene-derived hydroperoxy enals (HPALDs) photodissociate with a
  ! quantum yield of unity, i.e. more than two orders of magnitude higher than
  ! in the photolysis of monofunctional enones and enals such as
  ! methacrolein (MACR); Peeters and Muller et al., PCCP, 2011
  !
  PRJ(JK,JHPALDA) = 100.*PRJ(JK,JHPALDA)
  PRJ(JK,JHPALDB) = PRJ(JK,JHPALDA)

  ! VH test sensitivity to J-rate, considering that IFS(MOZ) rates are 20% lower than CB05
  ! PRJ(JK,JACH2O)=PRJ(JK,JACH2O)*0.8_JPRB
  ! PRJ(JK,JBCH2O)=PRJ(JK,JBCH2O)*0.8_JPRB
    
ENDDO ! layers

!
! correction photolysis for distance sun-earth:
!VH do this in CHEM_TM5.F90
!
!VH esrm2 = sundis(idate(2),idate(3))
!VH rj=rj*esrm2     
     

IF (LHOOK) CALL DR_HOOK('TM5_PHOTORATES_TROPO',1,ZHOOK_HANDLE )
END SUBROUTINE TM5_PHOTORATES_TROPO

