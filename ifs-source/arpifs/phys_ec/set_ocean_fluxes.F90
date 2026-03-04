! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SET_OCEAN_FLUXES(YDSURF,YDMCC,KDIM,SURFL,PSURF,FLUX)

   USE SURFACE_FIELDS_MIX , ONLY : TSURF
   USE PARKIND1, ONLY: JPRB, JPIM
   USE YOMHOOK,  ONLY: LHOOK, DR_HOOK, JPHOOK

   USE YOMPHYDER,ONLY: DIMENSION_TYPE, SURF_AND_MORE_LOCAL_TYPE, &
      &                SURF_AND_MORE_TYPE, FLUX_TYPE

   USE YOMCST,   ONLY: RLVTT, RLSTT, RSIGMA, RCPD

   USE YOMMCC,   ONLY: TMCC

   USE COUPLING, ONLY: LECEARTH





   IMPLICIT NONE

   ! Arguments
   TYPE(TSURF),                INTENT(INOUT) :: YDSURF
   TYPE(TMCC),                 INTENT(INOUT) :: YDMCC
   TYPE(DIMENSION_TYPE),          INTENT(IN) :: KDIM
   TYPE(SURF_AND_MORE_LOCAL_TYPE),INTENT(IN) :: SURFL
   TYPE(SURF_AND_MORE_TYPE),      INTENT(IN) :: PSURF
   TYPE(FLUX_TYPE),               INTENT(IN) :: FLUX

   ! Locals
   REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE
   INTEGER(KIND=JPIM) :: IL,IE,IG
   REAL(KIND=JPRB)    :: ZAHFLTI(KDIM%KLON,2)
   REAL(KIND=JPRB)    :: ZTS2(KDIM%KLON)
   REAL(KIND=JPRB)    :: ZTS3(KDIM%KLON)
   REAL(KIND=JPRB)    :: ZU10(KDIM%KLON)
   REAL(KIND=JPRB)    :: ZTMP(KDIM%KLON)

   IF (LHOOK) CALL DR_HOOK('SET_OCEAN_FLUXES',0,ZHOOK_HANDLE)
   ASSOCIATE(YSD_VD=>YDSURF%YSD_VD, LNEMOOCEICEMIX=>YDMCC%LNEMOOCEICEMIX, &
      &      LNEMOACCUMFLUX=>YDMCC%LNEMOACCUMFLUX, CPLNG_FLD=>YDMCC%CPLNG_FLD,&
      &      SP_SB=>YDSURF%SP_SB, YSP_SB=>YDSURF%YSP_SB, L2DECV2NEMO=>YDMCC%L2DECV2NEMO)

     ! =========================================================================
     ! *** Pre-compute indices
     ! =========================================================================

     IL = KDIM%KIDIA
     IE = KDIM%KFDIA - KDIM%KIDIA
     IG = KDIM%KSTGLO - 1 + KDIM%KIDIA

     ! =========================================================================
     ! *** Momentum fluxes (stresses)
     ! =========================================================================

     ZTMP(IL:IL+IE) = PSURF%PUSTRTI(IL:IL+IE,1)
     IF (LNEMOACCUMFLUX) THEN
        CPLNG_FLD(YDMCC%IP_A_TAUX_OCE)%D(IG:IG+IE,1,1) = &
           & CPLNG_FLD(YDMCC%IP_A_TAUX_OCE)%D(IG:IG+IE,1,1) + ZTMP(IL:IL+IE)
     ELSE
        CPLNG_FLD(YDMCC%IP_A_TAUX_OCE)%D(IG:IG+IE,1,1) = ZTMP(IL:IL+IE)
     ENDIF
     ZTMP(IL:IL+IE) = PSURF%PVSTRTI(IL:IL+IE,1)
     IF (LNEMOACCUMFLUX) THEN
        CPLNG_FLD(YDMCC%IP_A_TAUY_OCE)%D(IG:IG+IE,1,1) = &
           & CPLNG_FLD(YDMCC%IP_A_TAUY_OCE)%D(IG:IG+IE,1,1) + ZTMP(IL:IL+IE)
     ELSE
        CPLNG_FLD(YDMCC%IP_A_TAUY_OCE)%D(IG:IG+IE,1,1) = ZTMP(IL:IL+IE)
     ENDIF
     ZTMP(IL:IL+IE) = PSURF%PUSTRTI(IL:IL+IE,2)
     IF (LNEMOACCUMFLUX) THEN
        CPLNG_FLD(YDMCC%IP_A_TAUX_ICE)%D(IG:IG+IE,1,1) = &
           & CPLNG_FLD(YDMCC%IP_A_TAUX_ICE)%D(IG:IG+IE,1,1) + ZTMP(IL:IL+IE)
     ELSE
        CPLNG_FLD(YDMCC%IP_A_TAUX_ICE)%D(IG:IG+IE,1,1) = ZTMP(IL:IL+IE)
     ENDIF
     ZTMP(IL:IL+IE) = PSURF%PVSTRTI(IL:IL+IE,2)
     IF (LNEMOACCUMFLUX) THEN
        CPLNG_FLD(YDMCC%IP_A_TAUY_ICE)%D(IG:IG+IE,1,1) = &
           & CPLNG_FLD(YDMCC%IP_A_TAUY_ICE)%D(IG:IG+IE,1,1) + ZTMP(IL:IL+IE)
     ELSE
        CPLNG_FLD(YDMCC%IP_A_TAUY_ICE)%D(IG:IG+IE,1,1) = ZTMP(IL:IL+IE)
     ENDIF

     ! =========================================================================
     ! *** Radiative fluxes (solar, non-solar, dQ/dT)
     ! =========================================================================

     IF (LNEMOOCEICEMIX) THEN
        ZTMP(IL:IL+IE) = &
           &   SURFL%ZFRTI(IL:IL+IE,1)*SURFL%ZFRSOTI(IL:IL+IE,1) &
           & + SURFL%ZFRTI(IL:IL+IE,2)*SURFL%ZFRSOTI(IL:IL+IE,2)
        IF (LNEMOACCUMFLUX) THEN
           CPLNG_FLD(YDMCC%IP_A_QS_MIX)%D(IG:IG+IE,1,1) = &
              & CPLNG_FLD(YDMCC%IP_A_QS_MIX)%D(IG:IG+IE,1,1) + ZTMP(IL:IL+IE)
        ELSE
           CPLNG_FLD(YDMCC%IP_A_QS_MIX)%D(IG:IG+IE,1,1) = ZTMP(IL:IL+IE)
        ENDIF
     ELSE
        ZTMP(IL:IL+IE) = SURFL%ZFRSOTI(IL:IL+IE,1)
        IF (LNEMOACCUMFLUX) THEN
           CPLNG_FLD(YDMCC%IP_A_QS_OCE)%D(IG:IG+IE,1,1) = &
              & CPLNG_FLD(YDMCC%IP_A_QS_OCE)%D(IG:IG+IE,1,1) + ZTMP(IL:IL+IE)
        ELSE
           CPLNG_FLD(YDMCC%IP_A_QS_OCE)%D(IG:IG+IE,1,1) = ZTMP(IL:IL+IE)
        ENDIF
     ENDIF

     ZTMP(IL:IL+IE) = SURFL%ZFRSOTI(IL:IL+IE,2)
     IF (LNEMOACCUMFLUX) THEN
        CPLNG_FLD(YDMCC%IP_A_QS_ICE)%D(IG:IG+IE,1,1) = &
           & CPLNG_FLD(YDMCC%IP_A_QS_ICE)%D(IG:IG+IE,1,1) + ZTMP(IL:IL+IE)
     ELSE
        CPLNG_FLD(YDMCC%IP_A_QS_ICE)%D(IG:IG+IE,1,1) = ZTMP(IL:IL+IE)
     ENDIF

     ! Latent heat flux is computed from evaporation over water(1) and ice(2)
     ZAHFLTI(IL:IL+IE,1) = PSURF%PEVAPTI(IL:IL+IE,1) * RLVTT
     ZAHFLTI(IL:IL+IE,2) = PSURF%PEVAPTI(IL:IL+IE,2) * RLSTT

     ZTMP(IL:IL+IE) = PSURF%PAHFSTI(IL:IL+IE,2) + ZAHFLTI(IL:IL+IE,2) &
        &           + SURFL%ZAHFTRTI(IL:IL+IE,2)
     IF (LNEMOACCUMFLUX) THEN
        CPLNG_FLD(YDMCC%IP_A_QNS_ICE)%D(IG:IG+IE,1,1) = &
           & CPLNG_FLD(YDMCC%IP_A_QNS_ICE)%D(IG:IG+IE,1,1) + ZTMP(IL:IL+IE)
     ELSE
        CPLNG_FLD(YDMCC%IP_A_QNS_ICE)%D(IG:IG+IE,1,1) = ZTMP(IL:IL+IE)
     ENDIF

     IF (LNEMOOCEICEMIX) THEN
        ZTMP(IL:IL+IE) =  SURFL%ZFRTI(IL:IL+IE,2) * ZTMP(IL:IL+IE) &
           & + SURFL%ZFRTI(IL:IL+IE,1) * ( PSURF%PAHFSTI(IL:IL+IE,1)    &
           & + ZAHFLTI(IL:IL+IE,1) + SURFL%ZAHFTRTI(IL:IL+IE,1) )
        IF (LNEMOACCUMFLUX) THEN
           CPLNG_FLD(YDMCC%IP_A_QNS_MIX)%D(IG:IG+IE,1,1) = &
              & CPLNG_FLD(YDMCC%IP_A_QNS_MIX)%D(IG:IG+IE,1,1) + ZTMP(IL:IL+IE)
        ELSE
           CPLNG_FLD(YDMCC%IP_A_QNS_MIX)%D(IG:IG+IE,1,1) = ZTMP(IL:IL+IE)
        ENDIF
     ELSE
        ZTMP(IL:IL+IE) = ( PSURF%PAHFSTI(IL:IL+IE,1)    &
           & + ZAHFLTI(IL:IL+IE,1) + SURFL%ZAHFTRTI(IL:IL+IE,1) )
        IF (LNEMOACCUMFLUX) THEN
           CPLNG_FLD(YDMCC%IP_A_QNS_OCE)%D(IG:IG+IE,1,1) = &
              & CPLNG_FLD(YDMCC%IP_A_QNS_OCE)%D(IG:IG+IE,1,1) + ZTMP(IL:IL+IE)
        ELSE
           CPLNG_FLD(YDMCC%IP_A_QNS_OCE)%D(IG:IG+IE,1,1) = ZTMP(IL:IL+IE)
        ENDIF
     ENDIF

     ! Sensitivity of non-solar heat flux (only over ice)
     ZTS2(IL:IL+IE) = PSURF%PTSKTI(IL:IL+IE,2)**2
     ZTS3(IL:IL+IE) = PSURF%PTSKTI(IL:IL+IE,2)**3
     ZU10(IL:IL+IE) = SQRT(  PSURF%PSD_VD(IL:IL+IE,YSD_VD%Y10U%MP)**2 &
        &                     + PSURF%PSD_VD(IL:IL+IE,YSD_VD%Y10V%MP)**2 )

     ! From NEMO core bulk formulae
     ! Pay attention to the signs from the various contributions!
     ! NB 1.63e-3 is tuneable - ECEARTH have 1.4 e-3
     ZTMP(IL:IL+IE) = -4.00 * 0.95 * RSIGMA * ZTS3(IL:IL+IE)        &
        & -1.22 * RCPD * 1.63e-3 * ZU10(IL:IL+IE)       &
        & + RLSTT * 1.63e-3 * 11637800.                 &
        & * (-5897.8) * ZU10(IL:IL+IE)/ZTS2(IL:IL+IE) &
        & * EXP(-5897.8/PSURF%PTSKTI(IL:IL+IE,2))
     IF (LNEMOACCUMFLUX) THEN
        CPLNG_FLD(YDMCC%IP_A_DQNS_DT)%D(IG:IG+IE,1,1) = &
           & CPLNG_FLD(YDMCC%IP_A_DQNS_DT)%D(IG:IG+IE,1,1) + ZTMP(IL:IL+IE)
     ELSE
        CPLNG_FLD(YDMCC%IP_A_DQNS_DT)%D(IG:IG+IE,1,1) = ZTMP(IL:IL+IE)
     ENDIF

     ! =========================================================================
     ! *** Mass fluxes (runoff, precipitation, evaporation)
     ! =========================================================================

     ZTMP(IL:IL+IE) = FLUX%PFWRO1(IL:IL+IE) + FLUX%PFWROD(IL:IL+IE)
     IF (LNEMOACCUMFLUX) THEN
        CPLNG_FLD(YDMCC%IP_A_RUNOFF)%D(IG:IG+IE,1,1) = &
           & CPLNG_FLD(YDMCC%IP_A_RUNOFF)%D(IG:IG+IE,1,1) + ZTMP(IL:IL+IE)
     ELSE
        CPLNG_FLD(YDMCC%IP_A_RUNOFF)%D(IG:IG+IE,1,1) = ZTMP(IL:IL+IE)
     ENDIF

     ZTMP(IL:IL+IE) = FLUX%PFPLCL(IL:IL+IE,KDIM%KLEV) + FLUX%PFPLSL(IL:IL+IE,KDIM%KLEV)
     IF (LNEMOACCUMFLUX) THEN
        CPLNG_FLD(YDMCC%IP_A_PRECIP_LIQUID)%D(IG:IG+IE,1,1) = &
           & CPLNG_FLD(YDMCC%IP_A_PRECIP_LIQUID)%D(IG:IG+IE,1,1) + ZTMP(IL:IL+IE)
     ELSE
        CPLNG_FLD(YDMCC%IP_A_PRECIP_LIQUID)%D(IG:IG+IE,1,1) = ZTMP(IL:IL+IE)
     ENDIF

     ZTMP(IL:IL+IE) = FLUX%PFPLCN(IL:IL+IE,KDIM%KLEV) + FLUX%PFPLSN(IL:IL+IE,KDIM%KLEV)
     IF (LNEMOACCUMFLUX) THEN
        CPLNG_FLD(YDMCC%IP_A_PRECIP_SOLID)%D(IG:IG+IE,1,1) = &
           & CPLNG_FLD(YDMCC%IP_A_PRECIP_SOLID)%D(IG:IG+IE,1,1) + ZTMP(IL:IL+IE)
     ELSE
        CPLNG_FLD(YDMCC%IP_A_PRECIP_SOLID)%D(IG:IG+IE,1,1) = ZTMP(IL:IL+IE)
     ENDIF

     ZTMP(IL:IL+IE) =  - PSURF%PEVAPTI(IL:IL+IE,1) * SURFL%ZFRTI(IL:IL+IE,1) &
        &              - PSURF%PEVAPTI(IL:IL+IE,2) * SURFL%ZFRTI(IL:IL+IE,2)
     IF (LNEMOACCUMFLUX) THEN
        CPLNG_FLD(YDMCC%IP_A_EVAP_TOTAL)%D(IG:IG+IE,1,1) = &
           & CPLNG_FLD(YDMCC%IP_A_EVAP_TOTAL)%D(IG:IG+IE,1,1) + ZTMP(IL:IL+IE)
     ELSE
        CPLNG_FLD(YDMCC%IP_A_EVAP_TOTAL)%D(IG:IG+IE,1,1) = ZTMP(IL:IL+IE)
     ENDIF

     ZTMP(IL:IL+IE) = - PSURF%PEVAPTI(IL:IL+IE,2)
     IF (LNEMOACCUMFLUX) THEN
        CPLNG_FLD(YDMCC%IP_A_EVAP_ICE)%D(IG:IG+IE,1,1) = &
           & CPLNG_FLD(YDMCC%IP_A_EVAP_ICE)%D(IG:IG+IE,1,1) + ZTMP(IL:IL+IE)
     ELSE
        CPLNG_FLD(YDMCC%IP_A_EVAP_ICE)%D(IG:IG+IE,1,1) = ZTMP(IL:IL+IE)
     ENDIF

     IF (.NOT.LECEARTH) THEN

        ZTMP(IL:IL+IE) = PSURF%PSD_VD(KDIM%KIDIA:KDIM%KFDIA,YSD_VD%YTCC%MP)
        IF (LNEMOACCUMFLUX) THEN
           CPLNG_FLD(YDMCC%IP_A_TOTAL_CC)%D(IG:IG+IE,1,1) = &
              & CPLNG_FLD(YDMCC%IP_A_TOTAL_CC)%D(IG:IG+IE,1,1) + ZTMP(IL:IL+IE)
        ELSE
           CPLNG_FLD(YDMCC%IP_A_TOTAL_CC)%D(IG:IG+IE,1,1) = ZTMP(IL:IL+IE)
        ENDIF

        ZTMP(IL:IL+IE) = PSURF%PSD_VD(KDIM%KIDIA:KDIM%KFDIA,YSD_VD%YLCC%MP)
        IF (LNEMOACCUMFLUX) THEN
           CPLNG_FLD(YDMCC%IP_A_LOW_CC)%D(IG:IG+IE,1,1) = &
              & CPLNG_FLD(YDMCC%IP_A_LOW_CC)%D(IG:IG+IE,1,1) + ZTMP(IL:IL+IE)
        ELSE
           CPLNG_FLD(YDMCC%IP_A_LOW_CC)%D(IG:IG+IE,1,1) = ZTMP(IL:IL+IE)
        ENDIF

        ZTMP(IL:IL+IE) = SP_SB(KDIM%KIDIA:KDIM%KFDIA,1,YSP_SB%YTL%MP,KDIM%KBL)
        IF (LNEMOACCUMFLUX) THEN
           CPLNG_FLD(YDMCC%IP_A_IST_ATM)%D(IG:IG+IE,1,1) = &
              & CPLNG_FLD(YDMCC%IP_A_IST_ATM)%D(IG:IG+IE,1,1) + ZTMP(IL:IL+IE)
        ELSE
           CPLNG_FLD(YDMCC%IP_A_IST_ATM)%D(IG:IG+IE,1,1) = ZTMP(IL:IL+IE)
        ENDIF

     ENDIF
   END ASSOCIATE
   IF (LHOOK) CALL DR_HOOK('SET_OCEAN_FLUXES',1,ZHOOK_HANDLE)

END SUBROUTINE SET_OCEAN_FLUXES
