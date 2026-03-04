! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE TM5_SLINGO(KIDIA,KFDIA,KLON,KLEV,PIWC,PLWC,PCC,PGEOH,PRS1,PT, &
   & PTAUA_CLD,PTAUS_CLD,PPMCLD,PCLOUD_REFF )


!**   DESCRIPTION 
!     ----------
!
!   Part of TM5 routines for IFS chemistry: 
!   The evaluation of photolysis rates
!
! A. Slingo's data for cloud particle radiative properties (from 'A GCM 
! Parameterization for the Shortwave Properties of Water Clouds' JAS    
! vol. 46 may 1989 pp 1419-1427)
!
!
!------------------------------------------------------------------
!
!
!**   INTERFACE.
!     ----------
!          *TM5_SLINGO_FLUX* IS CALLED FROM *CHEM_TM5*.

! INPUTS:
! -------
! KIDIA :  Start of Array  
! KFDIA :  End  of Array 
! KLON  :  Length of Arrays 
! KLEV  :  NMBER OF LEVELS         (INPUT)
!
! PIWC  : ICWC
! PLWC  : LCWC
! PCC   : Cloud Fraction 
! PGEOH : height (0:KLEV)
! PRS1  : Half level pressure(0:KLEV)
! PT    : Temperature
!
! OUTPUTS:
! -------
!
! PTAUA_CLD   ! total absorption optical depth for cloud layer (liquid_cirrus)
! PTAUS_CLD   ! total scattering optical depth for cloud layer (liquid+cirrus)
! PPMCLD      ! phase function (HG)
! PCLOUD_REFF ! Cloud effective radius
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
!     F. Vana  05-Mar-2015  Support for single precision


USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOMCST   , ONLY : RNAVO, RG, RMD
USE TM5_PHOTOLYSIS , ONLY : ABAR,BBAR,CBAR,DBAR,EBAR,FBAR,A0,A1

IMPLICIT NONE

!-----------------------------------------------------------------------
!*       0.1  ARGUMENTS
!             ---------

INTEGER(KIND=JPIM),INTENT(IN) :: KIDIA, KFDIA, KLON, KLEV
REAL(KIND=JPRB),INTENT(IN)    :: PIWC(KLON,KLEV),PLWC(KLON,KLEV)
REAL(KIND=JPRB),INTENT(IN)    :: PCC(KLON,KLEV)
REAL(KIND=JPRB),INTENT(IN)    :: PGEOH(KLON,0:KLEV)
REAL(KIND=JPRB),INTENT(IN)    :: PRS1(KLON,0:KLEV)
REAL(KIND=JPRB),INTENT(IN)    :: PT(KLON,KLEV)

REAL(KIND=JPRB),INTENT(OUT)   :: PTAUA_CLD(KLON,KLEV)
REAL(KIND=JPRB),INTENT(OUT)   :: PTAUS_CLD(KLON,KLEV)
REAL(KIND=JPRB),INTENT(OUT)   :: PPMCLD(KLON,KLEV)
REAL(KIND=JPRB),INTENT(OUT)   :: PCLOUD_REFF(KLON,KLEV)


! * LOCAL 
! ABCDEF coefficients for current spectral interval 
REAL(KIND=JPRB)    :: ZABARI, ZBBARI, ZCBARI, ZDBARI, ZEBARI, ZFBARI 
!ZGC: asymmetry factor
!ZWC: single scattering albedo
!ZTOT: total optical depth of cloud layer
REAL(KIND=JPRB)    :: ZGC, ZWC, ZTOT
REAL(KIND=JPRB)    :: ZTMP1 ,ZTMP2, ZTMP3
REAL(KIND=JPRB)    :: ZLWP !cloud liquid water path [g/m^2]
REAL(KIND=JPRB)    :: ZRGI

REAL(KIND=JPRB)    :: ZFACTOR 
REAL(KIND=JPRB)    :: ZEFF_RAD ! effective radius linked to LWP
REAL(KIND=JPRB)    :: ZRHODZ,ZCIWC,ZXSA,ZPRES,ZAIRN,ZD_GE, ZDZ


REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! * counters
INTEGER(KIND=JPIM) :: JL, JLEV,INDXSL

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('TM5_SLINGO',0,ZHOOK_HANDLE )


!-------------------------------------------------------------------------           
!     Set index for cloud particle properties based on the wavelength,  
!     according to A. Slingo (1989) equations 1-3:                      
!     Use index 1 (0.25 to 0.69 micrometers) for visible                
!     Use index 2 (0.69 - 1.19 micrometers) for near-infrared           
!     Use index 3 (1.19 to 2.38 micrometers) for near-infrared          
!     Use index 4 (2.38 to 4.00 micrometers) for near-infrared          
                                                                        
INDXSL = 1_JPIM
                                                                        
!     Set cloud extinction optical depth, single scatter albedo,        
!     asymmetry parameter, and forward scattered fraction:              

ZABARI = ABAR(INDXSL)
ZBBARI = BBAR(INDXSL)
ZCBARI = CBAR(INDXSL)
ZDBARI = DBAR(INDXSL)
ZEBARI = EBAR(INDXSL)
ZFBARI = FBAR(INDXSL)
!
!
!--------------------------------------------------------------------------
!     In the parameterization of Fu(1996) wavelength dependant extinction maybe
!     calculated using a set of pre-defined constants:
!     
!      nm               a0                a1
!    250-300         -.236447E-03      .253817E+01
!    300-330         -.266955E-03      .254719E+01
!    330-360         -.293599E-03      .254540E+01
!    360-400         -.258858E-03      .253815E+01
!    400-440         -.106451E-03      .252684E+01
!    440-480          .129121E-03      .250410E+01
!    480-520         -.954458E-04      .252061E+01
!    520-570         -.303108E-04      .251805E+01
!    570-640          .982244E-04      .250875E+01     
!
!    used in Beta=IWC(a0+a1/Dg)
!--------------------------------------------------------------------------     
    
! Initialize values

PTAUA_CLD(KIDIA:KFDIA,1:KLEV) = 0. 
PTAUS_CLD(KIDIA:KFDIA,1:KLEV) = 0. 
PPMCLD(KIDIA:KFDIA,1:KLEV) = 0. 
PCLOUD_REFF(KIDIA:KFDIA,1:KLEV) = 0.     
       
!JEW The ice and water particles are now treated seperately. For cloud the values are taken from slingo 
!JEW the refractive index for ICE is very low below 750nm therefore T~T(scatt).
!JEW avoid negative input if it occurs  
!VH Check not needed?
!VH lwc=max(0.,lwc)
!VH iwc=max(0.,iwc)  

! * read in new cloud data to feed into the slingo routine. 
! * The values of lwc are zero in top ~ 1/4 levels (Stratosphere) so limit layer loop
ZRGI=1.0_JPRB/RG
ZFACTOR = 7.24E16*1.E6*RMD*29.2605/RNAVO

DO JLEV=KLEV/4,KLEV
  DO JL=KIDIA,KFDIA 
    ZPRES =  (PRS1(JL,JLEV) + PRS1(JL,JLEV-1)) * 0.5_JPRB ! hPa
    ZDZ   = ( PGEOH(JL,JLEV-1) - PGEOH(JL,JLEV) ) * ZRGI  ! delta height in meter
    ZLWP = 0._JPRB
    IF(PLWC(JL,JLEV)>1.0E-10) THEN      
      ! * calculate total water path locally : convert to g/m(2) for slingo input          
      ! * following the conversion procedure for LWC from ECMWF input from old cloud subroutine.
      ZRHODZ = ZFACTOR* ZPRES*LOG(PRS1(JL,JLEV)/(1.0_JPRB+PRS1(JL,JLEV-1)))
      ZLWP = ZRHODZ*PLWC(JL,JLEV)
    ENDIF 
    
!========================================================================================================
! there is a potential problem in that cloud frac values may occur with no associated clp value
! on TM5 vertical resolutions, therefore a filter w.r.t. both fraction and cloud liquid path are included. 
!========================================================================================================
                                                                     
    IF (PCC(JL,JLEV)> 0.02_JPRB .AND. PLWC(JL,JLEV) > 1.E-10_JPRB ) THEN  

! calculate constants which are used in slingo

      ZEFF_RAD=(11.0_JPRB*ZLWP)+4.0_JPRB
      
! prevent the radius of non-precipitation clouds being too big      
      
      ZEFF_RAD=MIN(12._JPRB,ZEFF_RAD)
      ZEFF_RAD=MAX(4._JPRB,ZEFF_RAD)
      
      ZTMP1 = ZABARI + ZBBARI/ZEFF_RAD 
      ZTMP2 = 1._JPRB - ZCBARI - ZDBARI*ZEFF_RAD  
      ZTMP3 = ZFBARI*ZEFF_RAD
                                                                              
!     Do not let single scatter albedo be 1; delta-eddington solution
!     for non-conservative case:                                                                                              
      ZWC   = MIN(ZTMP2,.999999_JPRB)        
      ZTOT  = ZLWP*ZTMP1        
      ZGC   = ZEBARI + ZTMP3 
! * : no wavelength dependence for the absorption or scattering effects due to liquid clouds !!!!!           
      PTAUA_CLD(JL,JLEV) = (1.-ZWC)*ZTOT 
      PTAUS_CLD(JL,JLEV) = ZWC*ZTOT
! * avoid possible negatives due to input data            
      PTAUA_CLD(JL,JLEV)=MAX(0.0_JPRB,PTAUA_CLD(JL,JLEV))
      
      PCLOUD_REFF(JL,JLEV) = ZEFF_RAD

! JEW : for calculating the scattering component due to ice 
      IF(PIWC(JL,JLEV)>1.E-10) THEN
        ZAIRN=7.24E16*ZPRES/PT(JL,JLEV)
        ! ZCIWC in g/m3 from TM5 definition ; PIWC was in kg/kg.      
        ZCIWC=PIWC(JL,JLEV)*ZAIRN*1E6_JPRB/RNAVO 
        ZXSA=1.0E-4_JPRB*ZCIWC**0.9_JPRB
        !
        ! calculate D_ge using the relationship in Fu (1996) where Beta=extinction co-efficient (m-1)
        !
        !  D_ge = 2(3)**0.5/(3 Rho)*(IWC/ZXSA)
        !
        ZD_GE=(2_JPRB*1.73205_JPRB/(3_JPRB*0.917_JPRB))*(ZCIWC/ZXSA)
        ! convert to uM
        ZD_GE=ZD_GE/100_JPRB
        !
        ! Cirrus scattering has a wavelength dependancy
        !
        PTAUS_CLD(JL,JLEV)=PTAUS_CLD(JL,JLEV)+((ZCIWC*(A0(3_JPIM)+(A1(3_JPIM)/ZD_GE)))*ZDZ)  
      ENDIF
      
      IF (PTAUS_CLD(JL,JLEV) > 0.) THEN 
        PPMCLD(JL,JLEV)=ZGC 
      ELSE 
        PPMCLD(JL,JLEV)=0._JPRB
      ENDIF 
    ENDIF 
  ENDDO 
ENDDO
   

IF (LHOOK) CALL DR_HOOK('TM5_SLINGO',1,ZHOOK_HANDLE )
END SUBROUTINE TM5_SLINGO

