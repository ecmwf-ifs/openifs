! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE TM5_AEROSOL_INFO ( KIDIA,KFDIA,KLON,KLEV, &
    &   PGEOH,PTP, PLSM,PRS1 ,PRSF1, PCC, PQP, PGELAT, &
    &   PTAUS_AER,PTAUA_AER, PMAER         )

 
!**   DESCRIPTION 
!     ----------
!
!   Part of TM5 routines for IFS chemistry: 
!   The evaluation of photolysis rates
!
!
! assignment of aerosol optical depths
! This is a crude method to provide some average values for the absorption and scattering
! terms by background aerosol choice of values defined according to the relative humidity.
! Absorption component can be set to zero throughout (M.van Weele, private comm.,2005).
! Scattering component chosen to be representative of background aerosol
! Moreover there is a choice as the whether the values defined by shettle and fenn are used
! This will require a look up table on 1x1,60 layers with indexes 1->4 with respect to
! location. At the moment background aerosol is chosen throughout
! JEW JUNE 2005
!
!------------------------------------------------------------------
!
!
!**   INTERFACE.
!     ----------
!          *TM5_AEROSOL_INFO* IS CALLED FROM *CHEM_tm5*.

! INPUTS:
! -------
! KIDIA :  Start of Array  
! KFDIA :  End  of Array 
! KLON  :  Length of Arrays 
! KLEV  :  NMBER OF LEVELS         (INPUT)
! PGEOH : Height info
! PTP   : Temp
! PLSM  : Land-sea mask (albedo??)
! PRS1  : Half-layer pressure
! PRSF1  : Full-layer pressure
! PCC     : Cloud fraction
! PQP     : Specific humidity
! PGELAT  : Latitude (radians)
!
! OUTPUTS:
! -------
!
! PMAER
! PTAUS_AER
! PTAUA_AER
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
!        ORIGINAL : 2012-07-25


USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOMCST   , ONLY : RPI, RG, RMD
USE TM5_PHOTOLYSIS , ONLY : PN_REF, NGRID,NBANDS_TROP, &
    &  LMID, LMID_GRIDA, SCA, ABS_EFF, GFAC
USE TM5_CHEM_MODULE, ONLY : XMH2O

IMPLICIT NONE

!-----------------------------------------------------------------------
!*       0.1  ARGUMENTS
!             ---------

INTEGER(KIND=JPIM),INTENT(IN) :: KIDIA, KFDIA, KLON, KLEV
REAL(KIND=JPRB), INTENT(IN)   :: PGEOH(KLON,0:KLEV)
REAL(KIND=JPRB), INTENT(IN)   :: PTP(KLON,KLEV)
REAL(KIND=JPRB), INTENT(IN)   :: PLSM(KLON)
REAL(KIND=JPRB), INTENT(IN)   :: PRS1(KLON,0:KLEV)
REAL(KIND=JPRB), INTENT(IN)   :: PRSF1(KLON,KLEV)
REAL(KIND=JPRB), INTENT(IN)   :: PCC(KLON,KLEV)
REAL(KIND=JPRB), INTENT(IN)   :: PQP(KLON,KLEV)
REAL(KIND=JPRB), INTENT(IN)   :: PGELAT(KLON)

REAL(KIND=JPRB), INTENT(OUT)   :: PTAUS_AER(KLON,KLEV,NBANDS_TROP,NGRID)
REAL(KIND=JPRB), INTENT(OUT)   :: PTAUA_AER(KLON,KLEV,NBANDS_TROP,NGRID)
REAL(KIND=JPRB), INTENT(OUT)   :: PMAER(KLON,KLEV,NBANDS_TROP,NGRID)



! * LOCAL 

LOGICAL                      :: LLAEROSOL,LLSHETTLE_AND_FENN

INTEGER(KIND=JPIM)           :: IAERO_CLIM(KLON)
INTEGER(KIND=JPIM)           :: I_TYPE, I_REF, I_REF2
REAL(KIND=JPRB)              :: ZPART(KLON,KLEV)
REAL(KIND=JPRB)              :: ZRHUM(KLON,KLEV)
REAL(KIND=JPRB)              :: ZLAY1,ZLAY2,ZSCALE_AERO, ZDZ
REAL(KIND=JPRB)              :: ZA_SC,ZB_SC,ZA_AB,ZB_AB,ZA_G,ZB_G
REAL(KIND=JPRB)              :: ZBSA,ZBAA,ZGA

REAL(KIND=JPRB)              :: ZRGI, ZTR, ZWV, ZDENS

REAL(KIND=JPRB),DIMENSION(8),PARAMETER  :: ZRH_REF = (/0., 50., 70., 80., 90., 95., 98., 99./)
! * Consider different aerosol types:                                          
!  1 = rural aerosol        
!  2 = maritime aerosol        
!  3 = urban aerosol        
!  4 = free troposphere aerosol        


REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

! * counters
INTEGER(KIND=JPIM) :: JL, JLEV, JK, JN, JB, JWAV

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('TM5_AEROSOL_INFO',0,ZHOOK_HANDLE )

LLAEROSOL = .FALSE.
LLSHETTLE_AND_FENN = .TRUE.

! initialize values

ZSCALE_AERO = 1.0

! the land-sea mask is used to ascribe either marine or rural aerosol for the bottom layers

ZRHUM = 0. 
  
ZRGI=1.0_JPRB/RG

         
IF (LLSHETTLE_AND_FENN) THEN
  
  ! * Definition of aerosol clim... 
  ! * Should idealy be taken from elsewhere
  DO JL=KIDIA,KFDIA
    IF(PGELAT(JL)/RPI * 180.<=-70._JPRB .OR. PGELAT(JL)/RPI * 180.>=70._JPRB) THEN
       ! Over NP and SP assume free tropospheric aerosol
       IAERO_CLIM(JL)=4
    ELSE    
       IF(PLSM(JL)<0.3)  IAERO_CLIM(JL)=2_JPIM
       IF(PLSM(JL)>=0.3) IAERO_CLIM(JL)=1_JPIM
       ! VH Find a way to select urban aerosol...
       ! VH IF( ??? ) IAERO_CLIM(JL)=3
    ENDIF 
  ENDDO
 
  DO JLEV=1,KLEV 
    DO JL=KIDIA,KFDIA 
      ZDENS = 7.2427E16 * PRSF1(JL,JLEV)/PTP(JL,JLEV)
      ZTR=1.-373.15/PTP(JL,JLEV)
      ZWV=EXP((((-.1299*ZTR-.6445)*ZTR-1.976)*ZTR+13.3185)*ZTR)  
      IF(PCC(JL,JLEV)<0.95_JPRB) ZRHUM(JL,JLEV)= PQP(JL,JLEV) *ZDENS*RMD/XMH2O*PTP(JL,JLEV)/(1013.25*ZWV*7.24E16)
      ! * assume saturation when cloud fraction is very high.
      IF(PCC(JL,JLEV)>=0.95_JPRB)  ZRHUM(JL,JLEV)=98.0
      ZRHUM(JL,JLEV)=MIN(98._JPRB,ZRHUM(JL,JLEV))
    ENDDO
  ENDDO
        
  !
  ! use land-sea mask to ascribe either marine or rural aerosol for different grid cells
  ! in the lower few KM
  !
  DO JLEV=1,KLEV
    DO JL=KIDIA,KFDIA
      ZDZ = (PGEOH(JL,JLEV-1) - PGEOH(JL,JLEV)) * ZRGI 
      IF(JLEV>KLEV*5/6) THEN
        ! ~ Near-surface...
        I_TYPE=IAERO_CLIM(JL)
      ELSE
        ! Free troposphere aerosol
        I_TYPE=4
      ENDIF 

      IF(I_TYPE<4) THEN
        ZPART(JL,JLEV) = ZSCALE_AERO*PN_REF(I_TYPE)*ZDZ*100.
      ELSE
        ZLAY1 = (PRS1(JL,JLEV)/PRS1(JL,KLEV))**3
        ZLAY2 = (PRS1(JL,JLEV-1)/PRS1(JL,KLEV))**3
        ZPART(JL,JLEV) = ZSCALE_AERO*PN_REF(I_TYPE)*0.5*(ZLAY1+ZLAY2)*ZDZ*100.
      ENDIF
      
      I_REF = 8
     
      DO JK = 1,8
        IF(ZRH_REF(JK) < ZRHUM(JL,JLEV)) I_REF = JK
      ENDDO
      
      DO JN=1,NGRID
        DO JB=1,NBANDS_TROP
        
          IF (JN==1) JWAV=LMID(JB)
          IF (JN==2) JWAV=LMID_GRIDA(JB) 
 
          I_REF2=MIN(I_REF+1,8_JPIM)
  
          ZA_SC = (SCA(JWAV,I_REF,I_TYPE)-SCA(JWAV,I_REF2,I_TYPE))/(ZRH_REF(I_REF)-ZRH_REF(I_REF2)+0.001)
    
          ZB_SC = SCA(JWAV,I_REF,I_TYPE)- ZA_SC*ZRH_REF(I_REF)
     
          ZA_AB = (ABS_EFF(JWAV,I_REF,I_TYPE)-ABS_EFF(JWAV,I_REF2,I_TYPE))/(ZRH_REF(I_REF)-ZRH_REF(I_REF2)+0.001)
    
          ZB_AB = ABS_EFF(JWAV,I_REF,I_TYPE)- ZA_AB*ZRH_REF(I_REF)
     
          ZA_G =(GFAC(JWAV,I_REF,I_TYPE)-GFAC(JWAV,I_REF2,I_TYPE))/(ZRH_REF(I_REF)-ZRH_REF(I_REF2)+0.001)
    
          ZB_G =GFAC(JWAV,I_REF,I_TYPE) - ZA_G*ZRH_REF(I_REF)    
     
          ZBSA = ZA_SC*ZRHUM(JL,JLEV) + ZB_SC
          ZBAA = ZA_AB*ZRHUM(JL,JLEV) + ZB_AB
          ZGA = ZA_G*ZRHUM(JL,JLEV) + ZB_G
       
          PTAUS_AER(JL,JLEV,JB,JN) = ZBSA*ZPART(JL,JLEV)
          PTAUA_AER(JL,JLEV,JB,JN) = ZBAA*ZPART(JL,JLEV)
        
          IF(PTAUS_AER(JL,JLEV,JN,1)>0.) THEN
            ! DO k = 1,nmom
            PMAER(JL,JLEV,JB,JN)=ZGA
            ! ENDDO
          ELSE
             PMAER(JL,JLEV,JB,JN) = 0._JPRB
          ENDIF 
            
        ENDDO !nbands_trop
      ENDDO !ngrid
          
    ENDDO ! JL
  ENDDO ! JLEV
   
   
! ENDIF ! shettle and fenn switch
! JEW switch for the aerosol absorption and scattering properties  
ELSEIF (LLAEROSOL) THEN

  DO JLEV=1,KLEV 
    DO JL=KIDIA,KFDIA 
      ZDENS = 7.2427E16_JPRB * PRSF1(JL,JLEV)/PTP(JL,JLEV)
      ZTR=1._JPRB-373.15/PTP(JL,JLEV)
      ZWV=EXP((((-.1299*ZTR-.6445)*ZTR-1.976)*ZTR+13.3185)*ZTR)  
      ZRHUM(JL,JLEV)=PQP(JL,JLEV)*ZDENS*RMD/XMH2O*PTP(JL,JLEV)/(1013.25*ZWV*7.24E16)
      IF(ZRHUM(JL,JLEV)>40. .AND. ZRHUM(JL,JLEV)<80.) PMAER(JL,JLEV,:,:) = 0.65
      IF(ZRHUM(JL,JLEV)>=80.)  PMAER(JL,JLEV,:,:) = 0.85
    ENDDO
  ENDDO    

ENDIF
   

IF (LHOOK) CALL DR_HOOK('TM5_AEROSOL_INFO',1,ZHOOK_HANDLE )
END SUBROUTINE TM5_AEROSOL_INFO

