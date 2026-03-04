! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE TM5_PHOTO_FLUX(KIDIA, KFDIA, KLON, KLEV, PTP, PCSZA, PALB,   & 
   &          PV3, PRSF1, PRS1,  &
   &          PTAUA_CLD,PTAUS_CLD,PPMCLD,PTAUA_AER,PTAUS_AER,PPMAER,& 
   &          PCC,PGEOH, PRJ )


!**   DESCRIPTION 
!     ----------
!
!   Part of TM5 routines for IFS chemistry: 
!   The evaluation of photolysis rates
!
!********************************************************************** 
!                                                                     * 
!     contact:                                                        * 
!                                                                     *
!     Jason Williams                                                  * 
!     Koninklijk Nederlands Meteorologisch instituut (KNMI)           *
!     Postbus 201                                                     *
!     3730 AE De Bilt                                                 *
!     tel +31 (0)30 2206748                                           * 
!     e-mail : williams@knmi.nl                                       *
!                                                                     *
!********************************************************************** 
! purpose:
!
! The calculation of the radiation fields.
!
! method: 
!
! Optical properties (optical thickness,cross sections, quantum yields)
! are calculated for KLEV layers. Radiation fluxes are derived for (0:KLEV)
! levels, including the surface and top of atmosphere. Thus, level l overlays layer l.
! Photolysis rates in a layer are evaluated by taking the average of the
! photolysis rates evaluated at the upper and lower boundaries of the layer.
!
! First, spectral calculations are performed for an atmosphere with absorption
! only, including gaseous absorption by O2 and O3 and aerosol absorption.
! Second, a correction is applied for scattering and surface reflection, per
! scattering band and using radiative transfer calculations at representative
! wavelengths 
!
! This photolysis module is based on:
!
!     Landgraf and Crutzen, 1998, J. Atmos. Sci, 55, 863-878.
!     Williams et al., 2006, Atmos.Chem.Phys., 6, 4137-4161.   
!
!     J. E. Williams, A. Strunk, V. Huijnen, and M. van Weele
!     Geosci. Model Dev., 5, 15-35, doi:10.5194/gmd-5-15-2012, 2012.
!
!------------------------------------------------------------------
!
!
!**   INTERFACE.
!     ----------
!          *TM5_PHOTO_FLUX* IS CALLED FROM *CHEM_TM5*.

! INPUTS:
! -------
! KIDIA :  Start of Array  
! KFDIA :  End  of Array 
! KLON  :  Length of Arrays 
! KLEV :  NMBER OF LEVELS         (INPUT)
!
! 
! PTP     (KLON,KLEV)         : TEMPERATURE                   (K)
! PCSZA(KLON)                 : COS of Solar Zenit Angle
! PALB(KLON)                  : Surface albedo
! PV3                         : O3 overhead column (kg/m2)
! PRSF1(KLON,KLEV)            : FULL-LEVEL PRESSURE           (Pa)
! PRS1(KLON,0:KLEV)           : HALF-LEVEL PRESSURE           (Pa)
! PTAUA_CLD, PTAUS_CLD,PPMCLD : Cloud optical properties: absorption, scattering, SSA
! PTAUS_AER, PTAUA_AER, PPMAER: Aerosol optical properties: absorption, scattering, SSA
! PCC     (KLON,KLEV)         : CLOUD FRACTION                0..1  
! PGEOH(KLON,0:KLEV)           : GEOPOTENTIAL                 (m*m/s*s)
!
!
! OUTPUTS:
! -------
! PRJ  (KLON,KLEV,NPHOTO)     : resulting photolysis rates
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
USE YOMCST   , ONLY : RPI, RG, RMD, RKBOL
USE TM5_PHOTOLYSIS , ONLY : XS_O3_LOOK, XS_HNO3_LOOK,XS_PAN_LOOK,QY_O3_LOOK, &
      & XS_H2O2_LOOK, XS_N2O5_LOOK, XS_NO2_LOOK,QY_NO2_LOOK,QY_CO_LOOK, XS_NO3_LOOK, &
      & XS_CH2O_LOOK, NPHOTO, MAXW, MAXWAV,NBANDS_TROP,NGRID, NTEMP_MAX, &
      & NWAV_NO2,NWAV_CO,NWAV_HNO3,NWAV_O3,NWAV_PAN,NWAV_H2O2,NWAV_N2O5,NWAV_CH2O,&
      & NWAV_NO3,NWAV_CH3COCHO, &
      & SZA_LIMIT, SZA_WIDELIMIT, PHI0, RD, &
      & A1_ACET, A2_ACET, A3_ACET, A4_ACET

IMPLICIT NONE

!-----------------------------------------------------------------------
!*       0.1  ARGUMENTS
!             ---------

INTEGER(KIND=JPIM),INTENT(IN) :: KIDIA, KFDIA, KLON, KLEV
REAL(KIND=JPRB),INTENT(IN)    :: PTP(KLON,KLEV)
REAL(KIND=JPRB),INTENT(IN)    :: PCSZA(KLON),PALB(KLON)
REAL(KIND=JPRB),INTENT(IN)    :: PRSF1(KLON,KLEV)
REAL(KIND=JPRB),INTENT(IN)    :: PRS1(KLON,0:KLEV)
!VH REAL(KIND=JPRB),INTENT(IN)    :: PRS1_KLEV(KLON)
REAL(KIND=JPRB),INTENT(IN)    :: PV3(KLON,KLEV)

REAL(KIND=JPRB),DIMENSION(KLON,KLEV),INTENT(IN)     :: PTAUA_CLD,PTAUS_CLD,PPMCLD,PCC
REAL(KIND=JPRB),DIMENSION(KLON,0:KLEV),INTENT(IN)     :: PGEOH
REAL(KIND=JPRB),DIMENSION(KLON,KLEV,NBANDS_TROP,NGRID),INTENT(IN) :: PTAUA_AER,PTAUS_AER,PPMAER

REAL(KIND=JPRB),INTENT(OUT)   :: PRJ(KLON,KLEV,NPHOTO)    

! * LOCAL 
!VH logical       :: LLCLEAR, LLCLOUDY
REAL(KIND=JPRB)    :: ZV2_COL(0:KLEV),ZV3_COL(0:KLEV)
REAL(KIND=JPRB)    :: ZTEMP_IND(NTEMP_MAX)
REAL(KIND=JPRB)    :: ZSOL(KLON),ZSOL_IN(KLON)
REAL(KIND=JPRB)    :: ZDV2_COL(KLEV),ZDV3_COL(KLEV)
REAL(KIND=JPRB)    :: ZCST_O3_COL(KLEV,MAXWAV)
REAL(KIND=JPRB)    :: ZCST_NO2_COL(KLEV,NWAV_NO2),ZQY_NO2_COL(KLEV,NWAV_NO2)
REAL(KIND=JPRB)    :: ZCST_HNO3_COL(KLEV,NWAV_HNO3)
REAL(KIND=JPRB)    :: ZCST_H2O2_COL(KLEV,NWAV_H2O2)
REAL(KIND=JPRB)    :: ZCST_N2O5_COL(KLEV,NWAV_N2O5)
REAL(KIND=JPRB)    :: ZCST_CH2O_COL (KLEV,NWAV_CH2O)
REAL(KIND=JPRB)    :: ZCST_PAN_COL(KLEV,NWAV_PAN)
REAL(KIND=JPRB)    :: ZCST_NO3_COL(KLEV,NWAV_NO3)
REAL(KIND=JPRB)    :: ZQY_O1D_COL(KLEV,NWAV_O3)
REAL(KIND=JPRB)    :: ZQY_CH3COCHO_COL(KLEV,NWAV_CH3COCHO)
REAL(KIND=JPRB)    :: ZQY_CO_COL(KLEV,NWAV_CO),ZQY_C2O3_COL(KLEV,NWAV_CO) 
REAL(KIND=JPRB)    :: ZFDIR(KLEV,MAXW),ZFACT(KLEV,NBANDS_TROP)
REAL(KIND=JPRB)    :: ZRJ_COLUMN(KLEV,NPHOTO)
REAL(KIND=JPRB)    :: ZLEVELS(0:KLEV)
REAL(KIND=JPRB)    :: ZRGI, ZQY, ZA1, ZA2, ZA4, ZPTORR, ZY_AIR,ZRDIF_SZALIMS
INTEGER(KIND=JPIM) :: ITABLE_POS
!
REAL(KIND=JPRB), PARAMETER   :: ZTOMOL = 1.255E+21   ! O3 from kg/m2 --> mol/cm2
REAL(KIND=JPRB)              :: ZSP  ! O2 column

REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

! * counters
INTEGER(KIND=JPIM) :: JL, JLEV, JK

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
#include "tm5_directflux.intfb.h"
#include "tm5_pifm_ran.intfb.h"
#include "tm5_photorates_tropo.intfb.h"

IF (LHOOK) CALL DR_HOOK('TM5_PHOTO_FLUX',0,ZHOOK_HANDLE )

!VH LLCLEAR = .FALSE.
!VH LLCLOUDY = .TRUE.

ZRGI=1.0_JPRB/RG
ZSP = 6.022E23*1.E-4*0.2095/(RMD*1.E-3*RG)

! Inverse of SZA difference in limits,i.e. 1./(94-85)
ZRDIF_SZALIMS=1.0_JPRB/(SZA_WIDELIMIT-SZA_LIMIT)

!*   Convert solar elevation to SZA
!DO JL=KSTART,KPROF
!  IF (ZSOL(JL) == RMDI) ZSOL(JL)=45.0_JPRB  !Set missing SOE values to 45 deg
!  ZSOL(JL) = 90._JPRB-PSOL(JL)
!  ZSZA_H(JL)=COS(RPI/180.0_JPRB*(ZSOL(JL)))
!  ZSZA(JL) = MAX(0.01_JPRB,ZSZA_H(JL))
!ENDDO

!*   Convert cosine SZA to degrees zenith angle
DO JL=KIDIA,KFDIA
  ZSOL(JL)=ACOS(PCSZA(JL))*180.0_JPRB/RPI
  ZSOL_IN(JL)=MIN(SZA_LIMIT,ZSOL(JL))
ENDDO
   
! define temperature index array from the lookup table
DO JK = 1, NTEMP_MAX
   ZTEMP_IND(JK) = 182.5 + (JK-1) * 5.0
ENDDO 


DO JL=KIDIA,KFDIA

  ! Compute photorates up to SZA_WIDELIMIT, i.e. 94 deg SZA...
  ! But using ZSOL_IN, i.e. solar zenith angle up to 85 
  IF (ZSOL(JL) < SZA_WIDELIMIT) THEN
    !
    ! initialise
    !
    ZFDIR = 0._JPRB 
    ZDV2_COL = 0._JPRB
    ZDV3_COL = 0._JPRB
    ZFACT = 0._JPRB

    ! Calculate MGLY & ACETONE press/ temp dependent
    ! quantum yields here
    DO JLEV = 1,KLEV
      ZPTORR= PRSF1(JL,JLEV)*0.00750061683_JPRB
      !second absorption band , wave < 420 nm  , Koch and Moortgat  
      DO JK = 1,39
       ZQY_CH3COCHO_COL(JLEV,JK)  = MIN(1._JPRB,(PHI0(JK) * RD(JK))/(RD(JK) + ZPTORR*PHI0(JK)) )
      ENDDO
      ! * switch to the formulation by Chen et al, JPC, 104, 
      ! * 11126-11131, 2000 for wav > 420nm
      DO JK = 40,NWAV_CH3COCHO
        ZQY_CH3COCHO_COL(JLEV,JK)   = MIN(1.0_JPRB, PHI0(JK)/(1.+(PHI0(JK)*RD(JK)*ZPTORR)))
      ENDDO

      ! Air concentration (molec/cm3)
      !ZY_AIR=ZRHO * RNAVO / (RMD*1e3) * 1e-6
      ZY_AIR= PRSF1(JL,JLEV)/(RKBOL*PTP(JL,JLEV))*1E-6

      ITABLE_POS=INT( ( PTP(JL,JLEV)- (ZTEMP_IND(1)-0.5_JPRB*5._JPRB) )/5_JPRB )+1_JPIM      
      ITABLE_POS=MIN(NTEMP_MAX,MAX(1_JPIM,ITABLE_POS))     
      DO JK=1,25
        ZQY_C2O3_COL(JLEV,JK)=0.0_JPRB
      ENDDO
      DO JK=26,35
        ZQY=(1.0_JPRB-QY_CO_LOOK(JK,ITABLE_POS))/(1_JPRB+A1_ACET(JK,ITABLE_POS)*ZY_AIR)
        ZQY_C2O3_COL(JLEV,JK)=MAX(0.0_JPRB,ZQY)
        ZQY_C2O3_COL(JLEV,JK)=MIN(1.0_JPRB,ZQY_C2O3_COL(JLEV,JK))
      ENDDO
      DO JK=36,56
        ! determine constants for calculating qy_c2o3
        ZA1=1_JPRB+(ZY_AIR*A4_ACET(JK,ITABLE_POS))+A3_ACET(JK,ITABLE_POS)
        ZA2=1_JPRB+(ZY_AIR*A2_ACET(JK,ITABLE_POS))+A3_ACET(JK,ITABLE_POS)
        ZA4=1_JPRB+(ZY_AIR*A4_ACET(JK,ITABLE_POS))
        ZQY=(1.0_JPRB-QY_CO_LOOK(JK,ITABLE_POS))*(ZA1/(ZA2*ZA4))
        ZQY_C2O3_COL(JLEV,JK)=MAX(0.0_JPRB,ZQY)
        ZQY_C2O3_COL(JLEV,JK)=MIN(1.0_JPRB,ZQY_C2O3_COL(JLEV,JK))
      ENDDO
      DO JK=57,NWAV_CO
        ZQY_C2O3_COL(JLEV,JK)=0.0_JPRB
      ENDDO


    ENDDO


        
    ! determine oxygen columns (can directly be calculated on levels)
    ZV2_COL(1:KLEV) = PRS1(JL,1:KLEV)*ZSP
    ! upper BC
    ZV2_COL(0)=ZV2_COL(1)*0.5
    ! Set O3 columns to correct units...
    ZV3_COL(1:KLEV) = PV3(JL,1:KLEV) * ZTOMOL ! Conversion of kg O3/m2 -> molec O3/ cm2
    ! upper BC
    ZV3_COL(0)=ZV3_COL(1)*0.5
    ! Zqy_ch3cocho_col(:,:) =  pqy_ch3cocho(JL,:,:)
    ! Zqy_c2o3_col(:,:) = pqy_c2o3(JL,:,:)
    ! column_cloud(:) = PCC(JL,:)
    ZLEVELS=PGEOH(JL,:) * ZRGI

     
    !
    ! assign the cross-sections here to avoid the exporting of large arrays
    !
    DO JLEV=1,KLEV
      ITABLE_POS=INT((PTP(JL,JLEV) - (ZTEMP_IND(1)-0.5*5_JPRB)) / 5_JPRB) + 1_JPIM
      ! stops spurious temperatures affecting run         
      ITABLE_POS=MIN(NTEMP_MAX,MAX(1_JPIM,ITABLE_POS))      
      ZCST_O3_COL(JLEV,1:MAXWAV)      = XS_O3_LOOK(1:MAXWAV,ITABLE_POS)
      ZCST_NO2_COL(JLEV,1:NWAV_NO2)   = XS_NO2_LOOK(1:NWAV_NO2,ITABLE_POS)
      ZCST_HNO3_COL(JLEV,1:NWAV_HNO3) = XS_HNO3_LOOK(1:NWAV_HNO3,ITABLE_POS)
      ZCST_H2O2_COL(JLEV,1:NWAV_H2O2) = XS_H2O2_LOOK(1:NWAV_H2O2,ITABLE_POS)
      ZCST_CH2O_COL(JLEV,1:NWAV_CH2O) = XS_CH2O_LOOK(1:NWAV_CH2O,ITABLE_POS)
      ZCST_N2O5_COL(JLEV,1:NWAV_N2O5) = XS_N2O5_LOOK(1:NWAV_N2O5,ITABLE_POS)
      ZCST_PAN_COL(JLEV,1:NWAV_PAN)   = XS_PAN_LOOK(1:NWAV_PAN,ITABLE_POS)
      ZCST_NO3_COL(JLEV,1:NWAV_NO3)   = XS_NO3_LOOK(1:NWAV_NO3,ITABLE_POS)
      ZQY_CO_COL(JLEV,1:NWAV_CO)      = QY_CO_LOOK(1:NWAV_CO,ITABLE_POS)
      ZQY_O1D_COL(JLEV,1:NWAV_O3)     = QY_O3_LOOK(1:NWAV_O3,ITABLE_POS)
      ZQY_NO2_COL(JLEV,1:NWAV_NO2)    = QY_NO2_LOOK(1:NWAV_NO2,ITABLE_POS)
    ENDDO
    
    ! slant column, lyman-alpha and direct flux ; all are zenith angle dependent

    CALL TM5_DIRECTFLUX( KLEV, ZSOL_IN(JL),PTP(JL,1:KLEV),ZV2_COL,ZV3_COL,ZCST_O3_COL, &
       &   ZFDIR,ZDV2_COL,ZDV3_COL,PTAUA_CLD(JL,1:KLEV), &
       &   PTAUA_AER(JL,1:KLEV,1:NBANDS_TROP,1:NGRID),ZLEVELS)

    ! Now call the PIFM RT solver for the online calculation ; 
    ! both clear and cloudy conditions can be used

    !   IF (LLCLEAR) THEN
    !     call TM5_PIFM(ZSOL(JL),PALB(JL),cst_o3_col,ZDV2_COL,ZDV3_COL, &
    !   & Ptaua_aer(JL,:,:,:),Ptaus_aer(JL,:,:,:),ppmaer(JL,:,:,:),fact)
    !   ENDIF
    !  IF (LLCLOUDY) THEN
        CALL TM5_PIFM_RAN(KLEV, ZSOL_IN(JL),PALB(JL),ZCST_O3_COL,ZDV2_COL,ZDV3_COL, &
           & PTAUA_CLD(JL,1:KLEV),PTAUS_CLD(JL,1:KLEV),PPMCLD(JL,1:KLEV),PTAUA_AER(JL,1:KLEV,1:NBANDS_TROP,1:NGRID), &
           & PTAUS_AER(JL,1:KLEV,1:NBANDS_TROP,1:NGRID),PPMAER(JL,1:KLEV,1:NBANDS_TROP,1:NGRID),ZFACT,PCC(JL,1:KLEV))              
    !  ENDIF

    ZRJ_COLUMN=0._JPRB
               
    CALL TM5_PHOTORATES_TROPO(KLEV, ZSOL_IN(JL),ZCST_O3_COL,ZCST_NO2_COL,ZCST_HNO3_COL,ZCST_H2O2_COL, &
           &    ZCST_CH2O_COL,ZCST_N2O5_COL,ZCST_PAN_COL,ZCST_NO3_COL,ZQY_NO2_COL,ZQY_O1D_COL,        &
           &    ZQY_CH3COCHO_COL,ZQY_CO_COL,ZQY_C2O3_COL, &
           &    ZFACT,ZFDIR,ZRJ_COLUMN,PTP(JL,1:KLEV)) 
     
    IF ( ZSOL(JL) <= SZA_LIMIT) THEN

      ! If actual SZA is below 85 deg then use values as such...
      PRJ(JL,1:KLEV,1:NPHOTO) = ZRJ_COLUMN
      
    ELSE

      ! Otherwise, with SZA is between 85 deg and 95 then use interpolated values 
      PRJ(JL,1:KLEV,1:NPHOTO) = ZRJ_COLUMN*(SZA_WIDELIMIT - ZSOL(JL)) * ZRDIF_SZALIMS

    ENDIF

  ELSE

    PRJ(JL,1:KLEV,1:NPHOTO) = 0.0_JPRB

  ENDIF
  
ENDDO





IF (LHOOK) CALL DR_HOOK('TM5_PHOTO_FLUX',1,ZHOOK_HANDLE )
END SUBROUTINE TM5_PHOTO_FLUX

