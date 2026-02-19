! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE TM5_DIRECTFLUX(KLEV, PANGLE,PT,PV2_COL,PV3_COL,PCST_O3_COL, &
       &   PFDIR,PDV2_COL,PDV3_COL,PTAUA_CLD, &
       &   PTAUA_AER,PGPH_DAT)


!**   DESCRIPTION 
!     ----------
!
!   Part of TM5 routines for IFS chemistry: 
!   The evaluation of photolysis rates
!
!         Schumann-Runge parameterization  
!         Koppers GAA, Murtagh DP, Ann. Geophysicae 14, 68-79          
!
! 
!**   INTERFACE.
!     ----------
!          *TM5_directflux* IS CALLED FROM *TM5_PHOTO_FLUX*.

! INPUTS:
! -------
! KLEV :  NMBER OF LEVELS         (INPUT)
! PANGLE : solar zenith angle
! PV2,PV3: oxygen and ozone column
!
!
! OUTPUTS:
! -------
!
! LOCAL:
! -------
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


USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOMCST   , ONLY : RPI
USE TM5_PHOTOLYSIS , ONLY : MAXWAV, MAXW, FLUX, LFIN, SZA_LIMIT, &
     & NBANDS_TROP,NGRID, CROSS_O2

IMPLICIT NONE

!-----------------------------------------------------------------------
!*       0.1  ARGUMENTS
!             ---------

INTEGER(KIND=JPIM),INTENT(IN) :: KLEV
REAL(KIND=JPRB),INTENT(IN)    :: PANGLE
REAL(KIND=JPRB),INTENT(IN)    :: PV2_COL(0:KLEV),PV3_COL(0:KLEV), PT(KLEV)
REAL(KIND=JPRB),INTENT(IN)    :: PCST_O3_COL(KLEV,MAXWAV)
REAL(KIND=JPRB),INTENT(OUT)   :: PFDIR(KLEV,MAXW)
REAL(KIND=JPRB),INTENT(OUT)   :: PDV2_COL(KLEV),PDV3_COL(KLEV)
REAL(KIND=JPRB),INTENT(IN)    :: PTAUA_CLD(KLEV)
REAL(KIND=JPRB),INTENT(IN)    :: PTAUA_AER(KLEV,NBANDS_TROP,NGRID)
REAL(KIND=JPRB),INTENT(IN)    :: PGPH_DAT(0:KLEV)


! * LOCAL 
REAL(KIND=JPRB),PARAMETER     :: ZA_CNST = 0.50572 
REAL(KIND=JPRB),PARAMETER     :: ZB_CNST = 6.07995 
REAL(KIND=JPRB),PARAMETER     :: ZC_CNST = 1.63640
REAL(KIND=JPRB),PARAMETER     :: Z_AE  = 6.371E6     ! earth radius [m] - take from ECMWF...
REAL(KIND=JPRB),PARAMETER     :: ZSM  = 1.0_JPRB
REAL(KIND=JPRB)               :: ZGAMMA, ZU0, ZRE
REAL(KIND=JPRB)               :: ZZE(0:KLEV)
REAL(KIND=JPRB)               :: ZRPSINZ,ZRJ,ZRJP1,ZDIFFJ,ZHEIGHT1,ZHEIGHT2
REAL(KIND=JPRB)               :: ZENRAD,ZDSJ
REAL(KIND=JPRB)               :: ZV2S(0:KLEV),ZV3S(0:KLEV),ZDV2S(KLEV),ZDV3S(KLEV),ZT(0:KLEV)
REAL(KIND=JPRB)               :: ZFDIR_TOP(MAXW),ZFDIR_BOT
REAL(KIND=JPRB)               :: ZTA_O2, ZTA_O3
REAL(KIND=JPRB)               :: ZDS(0:KLEV,KLEV)
INTEGER(KIND=JPIM)            :: JBANDNO


REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

! * counters
INTEGER(KIND=JPIM) :: JL, JK, JK1,JK2, JN

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('TM5_DIRECTFLUX',0,ZHOOK_HANDLE )

! initialise all output
PFDIR = 0._JPRB
PDV2_COL = 0._JPRB
PDV3_COL = 0._JPRB 
ZDS = 0._JPRB
!VH ds_tmp = 0._JPRB

 
! define temperature on grid levels. Note temperature on layers now has reversed vertical grid
ZT(0) = PT(1)
DO JK  = 1,KLEV-1
  ZT(JK) = (PT(JK) + PT(JK+1) ) * 0.5
ENDDO
ZT(KLEV) = PT(KLEV)

PDV2_COL(1)  = PV2_COL(1)      
PDV3_COL(1)  = PV3_COL(1) 
! Note: Original TM5 evaluation is performed upside down...
DO JK=2,KLEV
  PDV2_COL(JK)  = PV2_COL(JK)-PV2_COL(JK-1)
  PDV3_COL(JK)  = PV3_COL(JK)-PV3_COL(JK-1) 
ENDDO

!-----correction of the air mass factor---------------------------------
!     F. Kasten and T. Young, Revised optical air mass tabels and       
!     approximation formula (1989) Aplied Optics Vol 28, No. 22 P. 4735 
!     and J. Lenoble, atmospheric radiative transfer (1993), P. 236     
!    define air mass factor in mu

ZGAMMA = 90._JPRB - PANGLE   
ZU0 = SIN(ZGAMMA*RPI/180._JPRB) + ZA_CNST *(ZGAMMA+ZB_CNST)**(-ZC_CNST)     
ZU0 = MIN(1._JPRB,ZU0)    

IF(PANGLE <= 75._JPRB) THEN
  ZDS = 1._JPRB/ZU0
ELSEIF(PANGLE >75._JPRB .AND. PANGLE <=SZA_LIMIT) THEN
  !VH - no need for swap... PGPH_DAT(1:lm(region)+1)=PGPH_DAT(lm(region)+1:1:-1)
   
  ZRE = (Z_AE+PGPH_DAT(KLEV))*1.E-3
   
  DO JK=0,KLEV
      ZZE(JK) = (PGPH_DAT(JK)-PGPH_DAT(KLEV))*1.E-3
  ENDDO  
    
  ZENRAD=ACOS(ZU0)

  DO JK=0,KLEV
    ZRPSINZ = (ZRE+ZZE(JK))*SIN(ZENRAD)
    
    DO JN=1,JK
      ZRJ   = ZRE + ZZE(JN-1)
      ZRJP1 = ZRE + ZZE(JN)
    
      ZDIFFJ = ZZE(JN-1) - ZZE(JN)
      
      ZHEIGHT1 = ZRJ**2 - ZRPSINZ**2
      ZHEIGHT2 = ZRJP1**2 - ZRPSINZ**2

      ZHEIGHT1=MAX(0.0_JPRB,ZHEIGHT1)
      ZHEIGHT2=MAX(0.0_JPRB,ZHEIGHT2)
      
      !VH if(JK>k .and. n==JK) then
      !VH    dsj=sqrt(Zheight1)
      !VH else
        ZDSJ=SQRT(ZHEIGHT1)-ZSM*SQRT(ZHEIGHT2)
      !VH endif
        
      ZDS(JK,JN) = ZDSJ/ZDIFFJ
 
    ENDDO ! JN
  ENDDO !JK
  ! invert matrix to be compatible with lm(region)+1 being TOA
  ! VH - NOT NEEDED ds(1:lm(region)+1,1:lm(region))=ds_tmp(lm(region):0:-1,lm(region):1:-1)
ENDIF ! sza_rad  

  
! slant column : calculate in different way depending on PANGLE
IF(PANGLE <=75._JPRB) THEN 
  ZV2S(0)=ZDS(1,1)*PV2_COL(0)     
  ZV3S(0)=ZDS(1,1)*PV3_COL(0)    
  DO JK1=0,KLEV    
    ZV2S(JK1) = ZV2S(0) 
    ZV3S(JK1) = ZV3S(0) 
    DO JK2=1,JK1
      ZV2S(JK1) = ZV2S(JK1) + ZDS(JK1,JK2)*PDV2_COL(JK2)    
      ZV3S(JK1) = ZV3S(JK1) + ZDS(JK1,JK2)*PDV3_COL(JK2)    
    ENDDO          
  ENDDO
ELSEIF(PANGLE>75._JPRB .AND. PANGLE<=SZA_LIMIT) THEN
  ZV2S(0) = ZDS(1,1)*PV2_COL(0)            
  ZV3S(0) = ZDS(1,1)*PV3_COL(0)
  DO JK1=0,KLEV
    ZV2S(JK1) = ZV2S(0)
    ZV3S(JK1) = ZV3S(0)
    DO JK2=1,KLEV
      ZV2S(JK1) = ZV2S(JK1) + ZDS(JK1,JK2)*PDV2_COL(JK2)    
      ZV3S(JK1) = ZV3S(JK1) + ZDS(JK1,JK2)*PDV3_COL(JK2)
    ENDDO
  ENDDO    
ENDIF
 
!VH DO JK=0,KLEV
!VH   vo3s_col(JK)=V3S(JK)
!VH ENDDO 
  
! differential slant column
ZDV2S(1) = ZV2S(1)-ZV2S(0)
ZDV3S(1) = ZV3S(1)-ZV3S(0)    
DO JK = 2,KLEV
  ZDV2S(JK) = ZV2S(JK)-ZV2S(JK-1)     
  ZDV3S(JK) = ZV3S(JK)-ZV3S(JK-1)
  IF(ZDV2S(JK)<0.) ZDV2S(JK)=-1.0*ZDV2S(JK) 
  IF(ZDV3S(JK)<0.) ZDV3S(JK)=-1.0*ZDV3S(JK)   
ENDDO
 
! intialise direct flux
ZFDIR_TOP(1:MAXW) = FLUX(1:MAXW)     
! direct flux = actinic flux for a purely abs. atmosphere           
! here, layer averaged quantity 


JBANDNO=1

IF(PANGLE <= SZA_LIMIT) THEN
  DO JL = 1,MAXW
    DO JK = 1,KLEV          
      IF(JL<17_JPIM) THEN
        ZTA_O2 = CROSS_O2(JL)*ZDV2S(JK)
      ELSE
        ZTA_O2 = 0._JPRB
      ENDIF     
      ZTA_O3 = PCST_O3_COL(JK,JL)*ZDV3S(JK)
      ZFDIR_BOT  = ZFDIR_TOP(JL) * EXP(-ZTA_O2-ZTA_O3-PTAUA_CLD(JK)-PTAUA_AER(JK,JBANDNO,1_JPIM)) 
      !VH - limit flux to a non-zero positive value
      ZFDIR_BOT  = MAX(ZFDIR_BOT,1.E-10_JPRB)
      PFDIR(JK,JL) = (ZFDIR_TOP(JL) + ZFDIR_BOT)/2._JPRB
      ZFDIR_TOP(JL) = ZFDIR_BOT
      IF(JL>LFIN(JBANDNO)) JBANDNO=JBANDNO+1_JPIM
    ENDDO
  ENDDO
ENDIF 

IF (LHOOK) CALL DR_HOOK('TM5_DIRECTFLUX',1,ZHOOK_HANDLE )
END SUBROUTINE TM5_DIRECTFLUX

