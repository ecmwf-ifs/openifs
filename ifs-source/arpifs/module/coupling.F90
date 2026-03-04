! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction
! 
! (C) Copyright 1989- Meteo-France.
! 

MODULE COUPLING

   USE PARKIND1,ONLY : JPIM

   IMPLICIT NONE

   PRIVATE

   PUBLIC CPL_INIT
   PUBLIC CPL_CONFIG
   PUBLIC CPL_FINALIZE

   PUBLIC CPL_STAGE_OCE_SND
   PUBLIC CPL_STAGE_OCE_RCV
   PUBLIC CPL_STAGE_CHE_SND
   PUBLIC CPL_STAGE_CHE_RCV
   PUBLIC CPL_STAGE_VEG_SND
   PUBLIC CPL_STAGE_VEG_RCV

   PUBLIC CPL_NEMO_LIM
   PUBLIC LECEARTH

   LOGICAL, SAVE :: CPL_NEMO_LIM   = .FALSE.
   LOGICAL, SAVE :: LECEARTH       = .FALSE.

   NAMELIST /NAMCPLCFG/ CPL_NEMO_LIM, LECEARTH

   INTEGER(KIND=JPIM),PARAMETER :: CPL_STAGE_OCE_SND = 1
   INTEGER(KIND=JPIM),PARAMETER :: CPL_STAGE_OCE_RCV = 2
   INTEGER(KIND=JPIM),PARAMETER :: CPL_STAGE_CHE_SND = 3
   INTEGER(KIND=JPIM),PARAMETER :: CPL_STAGE_CHE_RCV = 4
   INTEGER(KIND=JPIM),PARAMETER :: CPL_STAGE_VEG_SND = 5
   INTEGER(KIND=JPIM),PARAMETER :: CPL_STAGE_VEG_RCV = 6

CONTAINS

! =============================================================================
! *** CPL_INIT
! =============================================================================
SUBROUTINE CPL_INIT(LDACTIVE)

   USE CPLNG, ONLY : CPLNG_INIT

   LOGICAL, INTENT(INOUT) :: LDACTIVE

   CHARACTER(LEN=255) :: CL_CPL_USE_CPLNG_ENV

   CALL GET_ENVIRONMENT_VARIABLE('CPL_USE_CPLNG',CL_CPL_USE_CPLNG_ENV)
    
   SELECT CASE (TRIM(CL_CPL_USE_CPLNG_ENV))

   CASE ('active','ACTIVE')

      LDACTIVE=.TRUE.

   CASE DEFAULT
       
      LDACTIVE=.FALSE.

   END SELECT
   CALL CPLNG_INIT(LDACTIVE=LDACTIVE)

END SUBROUTINE CPL_INIT

! =============================================================================
! *** CPL_CONFIG
! =============================================================================
SUBROUTINE CPL_CONFIG(YDGEOMETRY,YDMCC)

   USE GEOMETRY_MOD, ONLY : GEOMETRY

   USE YOMLUN, ONLY : NULOUT, NULNAM

   USE CPLNG, ONLY : CPLNG_IS_ACTIVE
   
   USE YOMMCC, ONLY : TMCC
   
   TYPE(GEOMETRY), INTENT(IN) :: YDGEOMETRY
   TYPE(TMCC), INTENT(IN) :: YDMCC

#include "posnam.intfb.h"
   
   ASSOCIATE(LNEMOCOUP=>YDMCC%LNEMOCOUP)

     ! -------------------------------------------------------------------------
     ! * (1) Read EC-Earth configuration namelist
     ! -------------------------------------------------------------------------

     IF (.NOT.CPLNG_IS_ACTIVE(YDMCC)) RETURN
     
     IF (.NOT.LNEMOCOUP) THEN
        CALL POSNAM(NULNAM,'NAMCPLCFG')
        READ(NULNAM,NAMCPLCFG)
     ENDIF
     WRITE(NULOUT,'(/,1X,''NAMELIST NAMCPLCFG'')')
     WRITE(NULOUT,'(5X,''CPL_NEMO_LIM   ='',L2)') CPL_NEMO_LIM
     WRITE(NULOUT,'(5X,''LECEARTH       ='',L2)') LECEARTH
     
     ! -------------------------------------------------------------------------
     ! * (2) Configure CPLNG
     ! -------------------------------------------------------------------------
     CALL CPL_CONFIG_COUPLING(YDGEOMETRY,YDMCC)

   END ASSOCIATE

END SUBROUTINE CPL_CONFIG

! =============================================================================
! *** CPL_CONFIG_COUPLING
! =============================================================================
SUBROUTINE CPL_CONFIG_COUPLING(YDGEOMETRY,YDMCC)
   
   USE GEOMETRY_MOD, ONLY : GEOMETRY
   
   USE YOMLUN, ONLY : NULOUT

   USE CPLNG, ONLY : CPLNG_IS_ACTIVE, CPLNG_ADD_FLD, CPLNG_ADD_FLD_COMPLETED, &
      & CPL_OUT, CPL_IN, CPL_OUTINST, CPLNG_FLD_TYPE_GRIDPOINT

   USE YOMMCC, ONLY : TMCC

   TYPE(GEOMETRY), INTENT(IN) :: YDGEOMETRY
   TYPE(TMCC) :: YDMCC

   ASSOCIATE(LNEMOCOUP=>YDMCC%LNEMOCOUP, LNEMO1WAY=>YDMCC%LNEMO1WAY, &
      &     LNEMOFLUXNC=>YDMCC%LNEMOFLUXNC, LMCCDYNSEAICE=>YDMCC%LMCCDYNSEAICE , &
      &     LNEMOLIMGET=>YDMCC%LNEMOLIMGET, LNEMOLIMPUT=>YDMCC%LNEMOLIMPUT, &
      &     LNEMOLIMALB=>YDMCC%LNEMOLIMALB, LNEMOLIMTEMP=>YDMCC%LNEMOLIMTEMP, &
      &     LNEMOLIMTHK=>YDMCC%LNEMOLIMTHK, LNEMOOCEICEMIX=>YDMCC%LNEMOOCEICEMIX, &
      &     LNEMOATMFLDS=>YDMCC%LNEMOATMFLDS,LNEMOLIMTLVL=>YDMCC%LNEMOLIMTLVL, &
      &     LNEMOACCUMFLUX=>YDMCC%LNEMOACCUMFLUX, L2DECV2NEMO=>YDMCC%L2DECV2NEMO )

     ! Coupling needs CPLNG. Check if it is active
     IF (.NOT.CPLNG_IS_ACTIVE(YDMCC)) THEN
        CALL ABOR1('CPL_CONFIG_COUPLING: CPLNG is needed for EC-Earth, but it is not active!')
     ENDIF
     
     ! -------------------------------------------------------------------------
     ! * (1) Configure coupling fields
     ! -------------------------------------------------------------------------
     
     ! 1.0 SST and sea-ice fraction: Needed in all setups
     CALL CPLNG_ADD_FLD(YDMCC,'A_SST',     CPLNG_FLD_TYPE_GRIDPOINT,CPL_IN,CPL_STAGE_OCE_RCV,&
        & KIDX=YDMCC%IP_A_SST)
     CALL CPLNG_ADD_FLD(YDMCC,'A_Ice_frac',CPLNG_FLD_TYPE_GRIDPOINT,CPL_IN,CPL_STAGE_OCE_RCV,&
        & KIDX=YDMCC%IP_A_ICE_FRAC)
     IF (.NOT.LECEARTH) THEN
        CALL CPLNG_ADD_FLD(YDMCC,'A_Curr_U',CPLNG_FLD_TYPE_GRIDPOINT,CPL_IN,CPL_STAGE_OCE_RCV,&
           & KIDX=YDMCC%IP_A_CURR_U)
        CALL CPLNG_ADD_FLD(YDMCC,'A_Curr_V',CPLNG_FLD_TYPE_GRIDPOINT,CPL_IN,CPL_STAGE_OCE_RCV,&
           & KIDX=YDMCC%IP_A_CURR_V)
     ENDIF

     ! 1.1 IFS-NEMO coupling
     IF (CPL_NEMO_LIM.OR.(LNEMOCOUP.AND.LMCCDYNSEAICE)) THEN
        
        IF (LECEARTH) THEN
           LNEMOLIMALB  = .TRUE.
           LNEMOLIMTEMP = .TRUE.
           LNEMOLIMTHK  = .TRUE.
           LNEMOACCUMFLUX = .FALSE.
           WRITE(UNIT=NULOUT,FMT='(" Resetting some YOMCC variables to configure EC-Earth coupling with NEMO")')
           WRITE(UNIT=NULOUT,FMT='(&
              &  " LNEMOCOUP = ",L2, &
              &  " LNEMO1WAY = ",L2," LNEMOFLUXNC = ",L2," LMCCDYNSEAICE = ",L2, &
              &  " LNEMOLIMGET = ",L2," LNEMOLIMPUT = ",L2, &
              &  " LNEMOLIMALB = ",L2," LNEMOLIMTEMP = ",L2," LNEMOLIMTHK = ",L2)')&
              &  LNEMOCOUP, LNEMO1WAY, LNEMOFLUXNC, LMCCDYNSEAICE, &
              &  LNEMOLIMGET, LNEMOLIMPUT,LNEMOLIMALB, LNEMOLIMTEMP,LNEMOLIMTHK
        ENDIF
       
        ! Fields sent from atmosphere to ocean
        CALL CPLNG_ADD_FLD(YDMCC,'A_TauX_oce',      CPLNG_FLD_TYPE_GRIDPOINT,CPL_OUT,CPL_STAGE_OCE_SND,&
           & KIDX=YDMCC%IP_A_TAUX_OCE)
        CALL CPLNG_ADD_FLD(YDMCC,'A_TauY_oce',      CPLNG_FLD_TYPE_GRIDPOINT,CPL_OUT,CPL_STAGE_OCE_SND,&
           & KIDX=YDMCC%IP_A_TAUY_OCE)
        CALL CPLNG_ADD_FLD(YDMCC,'A_TauX_ice',      CPLNG_FLD_TYPE_GRIDPOINT,CPL_OUT,CPL_STAGE_OCE_SND,&
           & KIDX=YDMCC%IP_A_TAUX_ICE)
        CALL CPLNG_ADD_FLD(YDMCC,'A_TauY_ice',      CPLNG_FLD_TYPE_GRIDPOINT,CPL_OUT,CPL_STAGE_OCE_SND,&
           & KIDX=YDMCC%IP_A_TAUY_ICE)
        IF (LNEMOOCEICEMIX) THEN
           CALL CPLNG_ADD_FLD(YDMCC,'A_Qs_mix',    CPLNG_FLD_TYPE_GRIDPOINT,CPL_OUT,CPL_STAGE_OCE_SND,&
              & KIDX=YDMCC%IP_A_QS_MIX)
           CALL CPLNG_ADD_FLD(YDMCC,'A_Qns_mix',   CPLNG_FLD_TYPE_GRIDPOINT,CPL_OUT,CPL_STAGE_OCE_SND,&
              & KIDX=YDMCC%IP_A_QNS_MIX)
        ELSE
           CALL CPLNG_ADD_FLD(YDMCC,'A_Qs_oce',    CPLNG_FLD_TYPE_GRIDPOINT,CPL_OUT,CPL_STAGE_OCE_SND,&
              & KIDX=YDMCC%IP_A_QS_OCE)
           CALL CPLNG_ADD_FLD(YDMCC,'A_Qns_oce',   CPLNG_FLD_TYPE_GRIDPOINT,CPL_OUT,CPL_STAGE_OCE_SND,&
              & KIDX=YDMCC%IP_A_QNS_OCE)
        ENDIF
        CALL CPLNG_ADD_FLD(YDMCC,'A_Qs_ice',        CPLNG_FLD_TYPE_GRIDPOINT,CPL_OUT,CPL_STAGE_OCE_SND,&
           & KIDX=YDMCC%IP_A_QS_ICE)
        CALL CPLNG_ADD_FLD(YDMCC,'A_Qns_ice',       CPLNG_FLD_TYPE_GRIDPOINT,CPL_OUT,CPL_STAGE_OCE_SND,&
           & KIDX=YDMCC%IP_A_QNS_ICE)
        CALL CPLNG_ADD_FLD(YDMCC,'A_Precip_liquid', CPLNG_FLD_TYPE_GRIDPOINT,CPL_OUT,CPL_STAGE_OCE_SND,&
           & KIDX=YDMCC%IP_A_PRECIP_LIQUID)
        CALL CPLNG_ADD_FLD(YDMCC,'A_Precip_solid',  CPLNG_FLD_TYPE_GRIDPOINT,CPL_OUT,CPL_STAGE_OCE_SND,&
           & KIDX=YDMCC%IP_A_PRECIP_SOLID)
        CALL CPLNG_ADD_FLD(YDMCC,'A_Evap_total',    CPLNG_FLD_TYPE_GRIDPOINT,CPL_OUT,CPL_STAGE_OCE_SND,&
           & KIDX=YDMCC%IP_A_EVAP_TOTAL)
        CALL CPLNG_ADD_FLD(YDMCC,'A_Evap_ice',      CPLNG_FLD_TYPE_GRIDPOINT,CPL_OUT,CPL_STAGE_OCE_SND,&
           & KIDX=YDMCC%IP_A_EVAP_ICE)
        CALL CPLNG_ADD_FLD(YDMCC,'A_dQns_dT',       CPLNG_FLD_TYPE_GRIDPOINT,CPL_OUT,CPL_STAGE_OCE_SND,&
           & KIDX=YDMCC%IP_A_DQNS_DT)
        CALL CPLNG_ADD_FLD(YDMCC,'A_Runoff',        CPLNG_FLD_TYPE_GRIDPOINT,CPL_OUT,CPL_STAGE_OCE_SND,&
           & KIDX=YDMCC%IP_A_RUNOFF)
        IF (.NOT.LECEARTH) THEN
           CALL CPLNG_ADD_FLD(YDMCC,'A_Total_CC',   CPLNG_FLD_TYPE_GRIDPOINT,CPL_OUT,CPL_STAGE_OCE_SND,&
              & KIDX=YDMCC%IP_A_TOTAL_CC)
           CALL CPLNG_ADD_FLD(YDMCC,'A_Low_CC',     CPLNG_FLD_TYPE_GRIDPOINT,CPL_OUT,CPL_STAGE_OCE_SND,&
              & KIDX=YDMCC%IP_A_LOW_CC)
           CALL CPLNG_ADD_FLD(YDMCC,'A_Ice_ATM',    CPLNG_FLD_TYPE_GRIDPOINT,CPL_OUT,CPL_STAGE_OCE_SND,&
              & KIDX=YDMCC%IP_A_IST_ATM)
        ENDIF
        ! Fields received by the atmosphere from the ocean                                            
        CALL CPLNG_ADD_FLD(YDMCC,'A_Ice_temp',      CPLNG_FLD_TYPE_GRIDPOINT,CPL_IN, CPL_STAGE_OCE_RCV,&
           & KIDX=YDMCC%IP_A_ICE_TEMP)
        CALL CPLNG_ADD_FLD(YDMCC,'A_Ice_albedo',    CPLNG_FLD_TYPE_GRIDPOINT,CPL_IN, CPL_STAGE_OCE_RCV,&
           & KIDX=YDMCC%IP_A_ICE_ALBEDO)
        CALL CPLNG_ADD_FLD(YDMCC,'A_Ice_thickness', CPLNG_FLD_TYPE_GRIDPOINT,CPL_IN, CPL_STAGE_OCE_RCV,&
           & KIDX=YDMCC%IP_A_ICE_THICKNESS)
        CALL CPLNG_ADD_FLD(YDMCC,'A_Snow_thickness',CPLNG_FLD_TYPE_GRIDPOINT,CPL_IN, CPL_STAGE_OCE_RCV,&
           & KIDX=YDMCC%IP_A_SNOW_THICKNESS)
        CALL CPLNG_ADD_FLD(YDMCC,'A_Ice_temp_lvls', CPLNG_FLD_TYPE_GRIDPOINT,CPL_IN, CPL_STAGE_OCE_RCV,&
           & KLVL=3, KIDX=YDMCC%IP_A_ICE_TEMP_LVLS)
        
     ELSEIF (LNEMOCOUP) THEN
        
        ! ECMWF no sea ice coupling.
        CALL CPLNG_ADD_FLD(YDMCC,'A_TauX',          CPLNG_FLD_TYPE_GRIDPOINT,CPL_OUT,CPL_STAGE_OCE_SND,&
           & KIDX=YDMCC%IP_A_TAUX)
        CALL CPLNG_ADD_FLD(YDMCC,'A_TauY',          CPLNG_FLD_TYPE_GRIDPOINT,CPL_OUT,CPL_STAGE_OCE_SND,&
           & KIDX=YDMCC%IP_A_TAUY)
        CALL CPLNG_ADD_FLD(YDMCC,'A_Qs',            CPLNG_FLD_TYPE_GRIDPOINT,CPL_OUT,CPL_STAGE_OCE_SND,&
           & KIDX=YDMCC%IP_A_QS)
        CALL CPLNG_ADD_FLD(YDMCC,'A_Qns',           CPLNG_FLD_TYPE_GRIDPOINT,CPL_OUT,CPL_STAGE_OCE_SND,&
           & KIDX=YDMCC%IP_A_QNS)
        CALL CPLNG_ADD_FLD(YDMCC,'A_Water',         CPLNG_FLD_TYPE_GRIDPOINT,CPL_OUT,CPL_STAGE_OCE_SND,&
           & KIDX=YDMCC%IP_A_WATER)
        
     ENDIF ! CPL_NEMO_LIM
    
     IF (LNEMOCOUP.AND.LNEMOATMFLDS) THEN
        ! Fields for diagnositics in NEMO
        CALL CPLNG_ADD_FLD(YDMCC,'A_SST_atm',  CPLNG_FLD_TYPE_GRIDPOINT,CPL_OUTINST,CPL_STAGE_OCE_SND,&
           & KIDX=YDMCC%IP_A_SST_ATM)
        CALL CPLNG_ADD_FLD(YDMCC,'A_TSK_atm',  CPLNG_FLD_TYPE_GRIDPOINT,CPL_OUTINST,CPL_STAGE_OCE_SND,&
           & KIDX=YDMCC%IP_A_TSK_ATM)
     ENDIF

     IF (L2DECV2NEMO) THEN
        ! Assimilation coupling variables
        CALL CPLNG_ADD_FLD(YDMCC,'A_2DECV_SKT', CPLNG_FLD_TYPE_GRIDPOINT,CPL_OUT,CPL_STAGE_OCE_SND,&
           & KIDX=YDMCC%IP_A_2DECV_SKT)
     ENDIF

     
     ! -------------------------------------------------------------------------
     ! * (2) Complete CPLNG configuration
     ! -------------------------------------------------------------------------
     CALL CPLNG_ADD_FLD_COMPLETED(YDMCC,YDGEOMETRY)
     
   END ASSOCIATE
  
END SUBROUTINE CPL_CONFIG_COUPLING

! =============================================================================
! *** CPL_FINALIZE
! =============================================================================
SUBROUTINE CPL_FINALIZE(LDACTIVE)

   USE CPLNG, ONLY : CPLNG_FINALIZE

   LOGICAL, INTENT(INOUT) :: LDACTIVE

   ! -------------------------------------------------------------------------
   ! * (1) Shut down CPLNG
   ! -------------------------------------------------------------------------
   CALL CPLNG_FINALIZE(LDACTIVE)

END SUBROUTINE CPL_FINALIZE

END MODULE COUPLING
