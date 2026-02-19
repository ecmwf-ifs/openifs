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

SUBROUTINE UPDCLIE_OASIS(YDGEOMETRY,YDSURF,YDMCC,YDRIP,YDERAD,YDDYNA,PTSTEP)

    USE PARKIND1, ONLY: JPRB, JPIM
    USE YOMHOOK,  ONLY: LHOOK, DR_HOOK, JPHOOK

    USE GEOMETRY_MOD , ONLY : GEOMETRY
    USE SURFACE_FIELDS_MIX , ONLY : TSURF

    USE YOMMCC   , ONLY : TMCC
    USE YOEPHY   , ONLY : YREPHY
    USE COUPLING , ONLY : CPL_NEMO_LIM, CPL_STAGE_OCE_RCV
    USE CPLNG    , ONLY : CPLNG_EXCHANGE, CPLNG_IDX
    USE YOMRIP   , ONLY : TRIP
    USE YOERAD   , ONLY : TERAD
    USE YOMDYNA  , ONLY : TDYNA

    IMPLICIT NONE

#include "surf_inq.h"

    ! Arguments
    TYPE(GEOMETRY) ,INTENT(IN)    :: YDGEOMETRY
    TYPE(TSURF)    ,INTENT(INOUT) :: YDSURF
    TYPE(TMCC)     ,INTENT(INOUT) :: YDMCC
    TYPE(TRIP)     ,INTENT(INOUT) :: YDRIP
    TYPE(TERAD)    ,INTENT(INOUT) :: YDERAD
    TYPE(TDYNA)    ,INTENT(INOUT) :: YDDYNA
    REAL(KIND=JPRB), INTENT(IN)   :: PTSTEP ! Time step

    ! Locals
    REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

    INTEGER(KIND=JPIM) :: IEND,IBL
    INTEGER(KIND=JPIM) :: JSTGLO,JROF

    REAL(KIND=JPRB) :: ZRTFREEZSICE,ZRCIMIN

    LOGICAL :: LL_SEA_POINTS(1:YDGEOMETRY%YRDIM%NPROMA)
    LOGICAL :: LL_ICE_POINTS(1:YDGEOMETRY%YRDIM%NPROMA)

    IF (LHOOK) CALL DR_HOOK('UPDCLIE_OASIS',0,ZHOOK_HANDLE)

    ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDGEM=>YDGEOMETRY%YRGEM)
    ASSOCIATE(NPROMA=>YDDIM%NPROMA,NGPTOT=>YDGEM%NGPTOT,YSURF=>YREPHY%YSURF, &
       & SD_VF=>YDSURF%SD_VF,SP_SB=>YDSURF%SP_SB,SP_RR=>YDSURF%SP_RR, &
       & YSD_VF=>YDSURF%YSD_VF,YSP_SB=>YDSURF%YSP_SB,YSP_RR=>YDSURF%YSP_RR)

    ! =========================================================================
    ! *** 0. Initialisation
    ! =========================================================================

    CALL SURF_INQ(YSURF,PRTFREEZSICE=ZRTFREEZSICE,PRCIMIN=ZRCIMIN)  

    ! =========================================================================
    ! *** 3. Update coupling fields (from CPLNG coupler)
    ! =========================================================================

    CALL CPLNG_EXCHANGE(YDMCC,YDRIP,YDERAD,YDDYNA,KSTAGE=CPL_STAGE_OCE_RCV)

    ! =========================================================================
    ! *** 4. Update IFS variables from coupling and CLIMR fields
    ! =========================================================================

    DO JSTGLO=1,NGPTOT,NPROMA

        IEND=MIN(NPROMA,NGPTOT-JSTGLO+1)
        IBL=(JSTGLO-1)/NPROMA+1

        ! Define sea and sea-ice points
        LL_SEA_POINTS(1:IEND) = SD_VF(1:IEND,YSD_VF%YLSM%MP,IBL) <= 0.5_JPRB
! alt.  LL_SEA_POINTS(1:IEND) = YDMCC%CPLNG_FLD(CPLNG_IDX(YDMCC,'A_SST'))%D(JSTGLO:JSTGLO+IEND-1,1,1) > 0.0_JPRB
        LL_ICE_POINTS(1:IEND) = LL_SEA_POINTS(1:IEND) .AND. &
           & YDMCC%CPLNG_FLD(CPLNG_IDX(YDMCC,'A_Ice_frac'))%D(JSTGLO:JSTGLO+IEND-1,1,1) > ZRCIMIN

        ! ---------------------------------------------------------------------
        ! *** 4.1 Set SST
        ! ---------------------------------------------------------------------
        WHERE (LL_SEA_POINTS(1:IEND))
            SD_VF(1:IEND,YSD_VF%YSST%MP,IBL) = &
               & YDMCC%CPLNG_FLD(CPLNG_IDX(YDMCC,'A_SST'))%D(JSTGLO:JSTGLO+IEND-1,1,1)
        ENDWHERE

        ! ---------------------------------------------------------------------
        ! *** 4.2 Set sea-ice fraction
        ! ---------------------------------------------------------------------
        WHERE (LL_ICE_POINTS(1:IEND))
            SD_VF(1:IEND,YSD_VF%YCI%MP,IBL) = MAX(0.0_JPRB,&
               & MIN(1.0_JPRB,YDMCC%CPLNG_FLD(CPLNG_IDX(YDMCC,'A_Ice_frac'))%D(JSTGLO:JSTGLO+IEND-1,1,1)))
        ELSEWHERE
            SD_VF(1:IEND,YSD_VF%YCI%MP,IBL) = 0.0_JPRB
        ENDWHERE

        ! ---------------------------------------------------------------------
        ! *** 4.3 Set sea-ice temperature
        ! ---------------------------------------------------------------------
        IF (CPL_NEMO_LIM) THEN
            WHERE (LL_ICE_POINTS(1:IEND))
                SP_SB(1:IEND,1,YSP_SB%YTL%MP,IBL) = &
                   & YDMCC%CPLNG_FLD(CPLNG_IDX(YDMCC,'A_Ice_temp'))%D(JSTGLO:JSTGLO+IEND-1,1,1)
            ELSEWHERE
                SP_SB(1:IEND,1,YSP_SB%YTL%MP,IBL) = ZRTFREEZSICE
                SP_SB(1:IEND,2,YSP_SB%YTL%MP,IBL) = ZRTFREEZSICE
                SP_SB(1:IEND,3,YSP_SB%YTL%MP,IBL) = ZRTFREEZSICE
                SP_SB(1:IEND,4,YSP_SB%YTL%MP,IBL) = ZRTFREEZSICE
            ENDWHERE
        ENDIF

        ! ---------------------------------------------------------------------
        ! *** 4.4 Set further derived variables
        ! ---------------------------------------------------------------------
        DO JROF=1,IEND

            ! Update surface and soil temp
            SP_SB(JROF,1,YSP_SB%YT%MP,IBL)=SD_VF(JROF,YSD_VF%YCI%MP,IBL)*SP_SB(JROF,1,YSP_SB%YTL%MP,IBL)+ &
             & (1.-SD_VF(JROF,YSD_VF%YCI%MP,IBL))*SD_VF(JROF,YSD_VF%YSST%MP,IBL)  
            SP_SB(JROF,2,YSP_SB%YT%MP,IBL)=SD_VF(JROF,YSD_VF%YCI%MP,IBL)*SP_SB(JROF,2,YSP_SB%YTL%MP,IBL)+ &
             & (1.-SD_VF(JROF,YSD_VF%YCI%MP,IBL))*SD_VF(JROF,YSD_VF%YSST%MP,IBL)  
            SP_SB(JROF,3,YSP_SB%YT%MP,IBL)=SD_VF(JROF,YSD_VF%YCI%MP,IBL)*SP_SB(JROF,3,YSP_SB%YTL%MP,IBL)+ &
             & (1.-SD_VF(JROF,YSD_VF%YCI%MP,IBL))*SD_VF(JROF,YSD_VF%YSST%MP,IBL)  
            SP_SB(JROF,4,YSP_SB%YT%MP,IBL)=SD_VF(JROF,YSD_VF%YCI%MP,IBL)*SP_SB(JROF,4,YSP_SB%YTL%MP,IBL)+ &
             & (1.-SD_VF(JROF,YSD_VF%YCI%MP,IBL))*SD_VF(JROF,YSD_VF%YSST%MP,IBL)  

            ! New skin temp
            SP_RR(JROF,YSP_RR%YT%MP,IBL) =                                                 &
            & (       SD_VF(JROF,YSD_VF%YCI%MP,IBL)  * SP_SB(JROF,1,YSP_SB%YTL%MP, IBL)**4 &
            &   + (1.-SD_VF(JROF,YSD_VF%YCI%MP,IBL)) * SD_VF(JROF,  YSD_VF%YSST%MP,IBL)**4 &
            & )**0.25_JPRB

        ENDDO ! JROF=1,IEND

    ENDDO ! JSTGLO=1,NGPTOT,NPROMA

    END ASSOCIATE
    END ASSOCIATE

    IF (LHOOK) CALL DR_HOOK('UPDCLIE_OASIS',1,ZHOOK_HANDLE)

END SUBROUTINE UPDCLIE_OASIS
