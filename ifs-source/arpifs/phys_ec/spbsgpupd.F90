! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SPBSGPUPD(YDDIM,YDML_PHY_STOCH,YDDYN,YDPHY2,KIDIA,KFDIA,KLON,KLEV,PVORTGRADX,PVORTGRADY,PVORT,PTENU,PTENV, &
  & PAPRSF,PT,PMFUDE_RATE,PMFU,PWMEAN,PTP,PQP,PRSF1,PSTREAM,PTEMP,PSTOPHCA,PTOTDISS,PTOTDISS_SMOOTH,&
  & PDISSGW,PDISSCU,PEXTRA,PEXTR2,KLEVX,KFLDX,KFLDX2,PGEMU,PU,PV,PAPRS,KTYPE)

!**** *SPBSGPUPD * - Update SPBS related grid point fields (dissipation, streamfunction forcing, ...)

!     PURPOSE.
!     --------
!       Update SPBS (Spectral stochastic backscatter) related grid point fields (dissipation, streamfunction forcing, ...)
!       Also inlcudes calculations for CA backscatter and vorticity confinement

!**   INTERFACE.
!     ----------
!        *CALL* *SPBSGPUPD(...)*

!     INPUT ARGUMENTS.

!     KIDIA            : start of horizontal loop
!     KFDIA            : start of horizontal loop
!     KLON             : horizontal dimension
!     KLEV             : end of vertical loop and vertical dimension
!     PVORTGRADX       : zonal vorticity gradient
!     PVORTGRADY       : meridional vorticity gradient
!     PVORT            : vorticity
!     PAPRSF           : PRESSURE ON FULL LEVELS.
!     PT               : TEMPERATURE
!     PMFUDE_RATE      : UD detrainmnet rate (KG/(M3*S))
!     PMFU             : Conv. mass flux up
!     PWMEAN           : vertically averaged conv. updraught speed (M/S)
!     PTP              :
!     PQP              :
!     PRSF1            : PROVISIONAL T+DT PRESSURE ON FULL LEVELS
!     PSTOPHCA         : CA pattern
!     PTOTDISS_SMOOTH  : smoothed total dissipation rate
!     PDISSGW          : gravity wave dissipation rate
!     PDISSCU          : cumulus convection dissipation rate diagnosed from convection scheme
!     PGEMU            : SINE OF LATITUDE
!     PU               : X-COMPONENT OF WIND.
!     PV               : Y-COMPONENT OF WIND.
!     PAPRS            : PRESSURE ON HALF-LEVELS.
!     KTYPE            : convection type

!
!     OUTPUT ARGUMENTS.

!     PTOTDISS         : total dissipation rate

!     IN+OUTPUT ARGUMENTS.

!     PTENU            : TENDENCY OF U-COMP. OF WIND
!     PTENV            : TENDENCY OF V-COMP. OF WIND
!     PSTREAM          : Streamfunction forcing/spectral pattern
!     PTEMP            : Temperature forcing/spectral pattern
!     PEXTRA           : EXTRA MULTI-LAYER FIELDS
!     PEXTR2           : EXTRA 2-D FIELDS

!
!     EXTERNALS.  NONE
!     ---------  

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     AUTHOR.
!     -------
!     09-Aug-2010  Martin Steinheimer

!     MODIFICATIONS.
!     --------------
!     13-Mar-2013 Martin Leutbecher     introduced LSPBS_DISSGW
!     K. Yessad (July 2014): Move some variables.
!     17-Oct-2014 M. Leutbecher: NSTOCHOPT=3, LSPBS_DISSCU, LSPBS_DISSNUM
!-----------------------------------------------------------------------

USE MODEL_PHYSICS_STOCHAST_MOD , ONLY : MODEL_PHYSICS_STOCHAST_TYPE
USE YOMDIM                     , ONLY : TDIM
USE PARKIND1                   , ONLY : JPIM     ,JPRB
!!USE YOMDYNA                  , ONLY : LRFRIC
USE YOMHOOK                    , ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOMCST                     , ONLY : RG, RCPD, RD, RETV, RA
USE YOMPHY2                    , ONLY : TPHY2
USE YOMCT3                     , ONLY : NSTEP
USE YOMDYN                     , ONLY : TDYN

!-----------------------------------------------------------------------

  IMPLICIT NONE

  TYPE(TDIM)        ,INTENT(IN)    :: YDDIM
 TYPE(TDYN)         ,INTENT(IN)    :: YDDYN
 TYPE(MODEL_PHYSICS_STOCHAST_TYPE),INTENT(IN)   :: YDML_PHY_STOCH
 TYPE(TPHY2)        ,INTENT(IN)    :: YDPHY2
  INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA
  INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA
  INTEGER(KIND=JPIM),INTENT(IN)    :: KLON
  INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV
  INTEGER(KIND=JPIM),INTENT(IN)    :: KLEVX 
  INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDX 
  INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDX2
  REAL(KIND=JPRB)   ,INTENT(IN)    :: PVORTGRADX(KLON,KLEV)
  REAL(KIND=JPRB)   ,INTENT(IN)    :: PVORTGRADY(KLON,KLEV)
  REAL(KIND=JPRB)   ,INTENT(IN)    :: PVORT(KLON,KLEV)
  REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTENU(KLON,KLEV)
  REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTENV(KLON,KLEV)
  REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPRSF(KLON,KLEV)
  REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KLON,KLEV)
  REAL(KIND=JPRB)   ,INTENT(IN)    :: PMFUDE_RATE(KLON,KLEV)     ! UD detrainmnet rate (KG/(M3*S))
  REAL(KIND=JPRB)   ,INTENT(IN)    :: PMFU(KLON,KLEV)
  REAL(KIND=JPRB)   ,INTENT(IN)    :: PWMEAN(KLON)
  REAL(KIND=JPRB)   ,INTENT(IN)    :: PTP(KLON,KLEV)
  REAL(KIND=JPRB)   ,INTENT(IN)    :: PQP(KLON,KLEV)
  REAL(KIND=JPRB)   ,INTENT(IN)    :: PRSF1(KLON,KLEV)
  REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTOTDISS(KLON,KLEV)
  REAL(KIND=JPRB)   ,INTENT(IN)    :: PTOTDISS_SMOOTH(KLON,KLEV)
  REAL(KIND=JPRB)   ,INTENT(IN)    :: PDISSGW(KLON,KLEV)
  REAL(KIND=JPRB)   ,INTENT(IN)    :: PDISSCU(KLON,KLEV)
  REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSTREAM(KLON,KLEV)
  REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTEMP(KLON,KLEV)
  REAL(KIND=JPRB)   ,INTENT(IN)    :: PSTOPHCA(KLON)
  REAL(KIND=JPRB)   ,INTENT(INOUT) :: PEXTRA(KLON,KLEVX,KFLDX) 
  REAL(KIND=JPRB)   ,INTENT(INOUT) :: PEXTR2(KLON,KFLDX2)
  REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEMU(KLON)
  REAL(KIND=JPRB)   ,INTENT(IN)    :: PU(KLON,KLEV) 
  REAL(KIND=JPRB)   ,INTENT(IN)    :: PV(KLON,KLEV) 
  REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPRS(KLON,0:KLEV)
  INTEGER(KIND=JPIM),INTENT(IN)    :: KTYPE(KLON)

!-----------------------------------------------------------------------

  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  INTEGER(KIND=JPIM) :: JK, JL

  REAL(KIND=JPRB) :: ZEGENVC(KLON,KLEV)         ! energy input rate due to vorticity confinement
  REAL(KIND=JPRB) :: ZDISSNUM(KLON,KLEV)        ! kinetic energy/unit mass dissipated in one t-step (numerics)
  REAL(KIND=JPRB) :: ZDISSCONV_DEEP(KLON,KLEV)  ! kinetic energy/unit mass detrained in one t-step (deep convection)
  REAL(KIND=JPRB) :: ZDISSCONV_SHAL(KLON,KLEV)  ! kinetic energy/unit mass detrained in one t-step (shallow convection)
  REAL(KIND=JPRB) :: ZDISS2D(KLON), ZDISSNUM2D(KLON), ZDISSGW2D(KLON), ZDISSCONV2D_DEEP(KLON)
  REAL(KIND=JPRB) :: ZDISSCONV2D_SHAL(KLON), ZFIELD2D(KLON)
  REAL(KIND=JPRB) :: ZEGEN2D(KLON)              ! depth-integrated KE production due to streamfunction forcing
  REAL(KIND=JPRB) :: ZEGENVC2D(KLON)            ! depth-integrated KE production due to vorticity con.
  REAL(KIND=JPRB) :: ZVCX(KLON,KLEV), ZVCY(KLON,KLEV) ! hold dudt and dvdt from vorticity confinement 
  REAL(KIND=JPRB) :: ZVC
  REAL(KIND=JPRB) :: ZVORTGRAD                  ! hold the modulus of the vorticity gradient 
  REAL(KIND=JPRB) :: ZRHO, ZDELP, ZSCALE
  REAL(KIND=JPRB) :: ZALPHAINVDEEP2,ZALPHAINVSHAL2,ZCONS_AMP_SPBS,ZCONS_AMP_CASBS,ZCONS_AMP_SPBS_T

  REAL(KIND=JPRB) :: ZRG, ZRCPD

!-------------------------------------------------------------

#include "abor1.intfb.h"

!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SPBSGPUPD',0,ZHOOK_HANDLE)
ASSOCIATE(NSMAX=>YDDIM%NSMAX, &
 & NMAXLEVRF=>YDDYN%NMAXLEVRF, &
 & ALPHA_DEEP_CONV=>YDML_PHY_STOCH%YRSTOPH%ALPHA_DEEP_CONV, &
 & ALPHA_SHAL_CONV=>YDML_PHY_STOCH%YRSTOPH%ALPHA_SHAL_CONV, &
 & AMAGSTOPH_CASBS=>YDML_PHY_STOCH%YRSTOPH%AMAGSTOPH_CASBS, BIHARM=>YDML_PHY_STOCH%YRSTOPH%BIHARM, &
 & LEXTRAFIELDS=>YDML_PHY_STOCH%YRSTOPH%LEXTRAFIELDS, LSPBS_DISSCU=>YDML_PHY_STOCH%YRSTOPH%LSPBS_DISSCU, &
 & LSPBS_DISSGW=>YDML_PHY_STOCH%YRSTOPH%LSPBS_DISSGW, LSPBS_DISSNUM=>YDML_PHY_STOCH%YRSTOPH%LSPBS_DISSNUM, &
 & LSTOPH_CASBS=>YDML_PHY_STOCH%YRSTOPH%LSTOPH_CASBS, LSTOPH_SPBS=>YDML_PHY_STOCH%YRSTOPH%LSTOPH_SPBS, &
 & LSTOPH_SPBS_T=>YDML_PHY_STOCH%YRSTOPH%LSTOPH_SPBS_T, LVORTCON=>YDML_PHY_STOCH%YRSTOPH%LVORTCON, &
 & NFRSTOPH_SPBS=>YDML_PHY_STOCH%YRSTOPH%NFRSTOPH_SPBS, NSTOCHOPT=>YDML_PHY_STOCH%YRSTOPH%NSTOCHOPT, &
 & RATIO_APE2KE=>YDML_PHY_STOCH%YRSTOPH%RATIO_APE2KE, RATIO_BACKSCAT=>YDML_PHY_STOCH%YRSTOPH%RATIO_BACKSCAT, &
 & RATIO_BACKSCAT_CON2NUM=>YDML_PHY_STOCH%YRSTOPH%RATIO_BACKSCAT_CON2NUM, &
 & RFLUX_DET_CLIP=>YDML_PHY_STOCH%YRSTOPH%RFLUX_DET_CLIP, VC_CON=>YDML_PHY_STOCH%YRSTOPH%VC_CON, &
 & TSPHY=>YDPHY2%TSPHY)
!-----------------------------------------------------------------------

  ZRG=1.0_JPRB/RG
  ZRCPD=1.0_JPRB/RCPD

!     -----------------------------------------------------------------
!     -----------------------------------------------------------------
!          1    stochastic backscatter calculation
!     -----------------------------------------------------------------
!     -----------------------------------------------------------------
  IF (LSTOPH_SPBS.OR.LSTOPH_CASBS) THEN
!     -----------------------------------------------------------------
!*        Initialize local fields
!     -----------------------------------------------------------------
    DO JL=KIDIA,KFDIA
      ZDISS2D(JL)    =   0._JPRB    ! initialize arrays to hold depth-means
      ZDISSNUM2D(JL) =   0._JPRB
      ZDISSGW2D(JL)  = 0._JPRB
      ZDISSCONV2D_DEEP(JL)= 0._JPRB
      ZDISSCONV2D_SHAL(JL)= 0._JPRB
      ZEGEN2D(JL)= 0._JPRB
    ENDDO
    DO JK=1,KLEV
      DO JL=KIDIA,KFDIA
        ZDISSNUM(JL,JK)  = 0.0_JPRB 
        ZDISSCONV_DEEP(JL,JK) = 0.0_JPRB 
        ZDISSCONV_SHAL(JL,JK) = 0.0_JPRB 
        ZFIELD2D(JL)=0.0_JPRB
      ENDDO
    ENDDO
!     -----------------------------------------------------------------
!      first set of extrafields (input fields)
!     -----------------------------------------------------------------
    IF (LEXTRAFIELDS) THEN
      DO JK=1,KLEV
        DO JL=KIDIA,KFDIA
          PEXTRA(JL,JK, 2)    = PTOTDISS(JL,JK)        !write before updating to be consistent with smoothed
          PEXTRA(JL,JK, 3)    = PTOTDISS_SMOOTH(JL,JK)
          IF (LSTOPH_SPBS) THEN
            PEXTRA(JL,JK, 4)    = PSTREAM(JL,JK)         !write Pattern
            IF (LSTOPH_SPBS_T) THEN
              PEXTRA(JL,JK, 5)    = PTEMP(JL,JK)
            ELSE
              PEXTRA(JL,JK, 5)    = 0
            ENDIF
          ENDIF
        ENDDO
      ENDDO
    ENDIF

!     -----------------------------------------------------------------
!      DISSIPATION RATE CALCULATION
!     -----------------------------------------------------------------
    IF (MOD(NSTEP+1,NFRSTOPH_SPBS) == 0) THEN
!     -----------------------------------------------------------------
!       estimate numerical dissipation rate using mod(vorticity gradient)**2
!     -----------------------------------------------------------------
      IF (LSPBS_DISSNUM) THEN
        DO JK=1,KLEV
          DO JL=KIDIA,KFDIA
            ZDISSNUM(JL,JK)=  BIHARM*(PVORTGRADX(JL,JK)**2 + PVORTGRADY(JL,JK)**2)  
          ENDDO
        ENDDO
      ENDIF

!     -----------------------------------------------------------------
!       estimate convective dissipation rate
!     -----------------------------------------------------------------
      IF (NSTOCHOPT==1) THEN
        ZALPHAINVDEEP2=1./ALPHA_DEEP_CONV**2 
        ZALPHAINVSHAL2=1./ALPHA_SHAL_CONV**2
        DO JK=1,KLEV
          DO JL=KIDIA,KFDIA
            ZRHO= PAPRSF(JL,JK)/(RD*PT(JL,JK))
            IF (KTYPE(JL) ==1) THEN  
              ZDISSCONV_DEEP(JL,JK)     = TSPHY*MIN(PMFUDE_RATE(JL,JK),RFLUX_DET_CLIP)/ &
                 &(ZRHO**3)*PMFU(JL,JK)**2*ZALPHAINVDEEP2
            ELSEIF (KTYPE(JL) ==2) THEN  !  correction *GJS* 2/1/08
              ZDISSCONV_SHAL(JL,JK) = TSPHY*MIN(PMFUDE_RATE(JL,JK),RFLUX_DET_CLIP)/  &
                 &(ZRHO**3)*PMFU(JL,JK)**2*ZALPHAINVSHAL2
            ENDIF
          ENDDO
        ENDDO
      ELSEIF (NSTOCHOPT==2) THEN
        DO JK=1,KLEV
          DO JL=KIDIA,KFDIA
            ZRHO= RD*PTP(JL,JK)*(1.0_JPRB+RETV*PQP(JL,JK))/PRSF1(JL,JK) ! 1/rho 
            IF (KTYPE(JL)==1) THEN
               ZDISSCONV_DEEP(JL,JK)=MIN(TSPHY*PMFUDE_RATE(JL,JK)*ZRHO,1.00_JPRB)*MIN(PWMEAN(JL),17._JPRB)**2
            ELSEIF (KTYPE(JL)==2) THEN
               ZDISSCONV_SHAL(JL,JK)=TSPHY*MIN(PMFUDE_RATE(JL,JK),RFLUX_DET_CLIP)*ZRHO*PWMEAN(JL)**2*ALPHA_SHAL_CONV
            ENDIF
          ENDDO
        ENDDO
      ELSEIF (NSTOCHOPT==3) THEN
        DO JK=1,KLEV
          DO JL=KIDIA,KFDIA
            IF (KTYPE(JL)==1) THEN
               ZDISSCONV_DEEP(JL,JK)=PDISSCU(JL,JK)
            ENDIF
          ENDDO
        ENDDO
      ELSE
        CALL ABOR1('SPBSGPUPD: NSTOCHOPT not in {1,2,3}.')
      ENDIF

!     -----------------------------------------------------------------
!       compute total dissipation rate
!     -----------------------------------------------------------------
      IF (LSPBS_DISSNUM) THEN
        ! contribution from horizontal diffusion/numerics included
        DO JK=1,KLEV
          DO JL=KIDIA,KFDIA
            PTOTDISS(JL,JK)   = ZDISSNUM(JL,JK)  
          ENDDO
        ENDDO
      ELSE
        PTOTDISS(:,:)=0._JPRB
      ENDIF
      IF (LSPBS_DISSGW) THEN
        ! orographic gravity wave drag included
        DO JK=1,KLEV
          DO JL=KIDIA,KFDIA
            PTOTDISS(JL,JK)   = PTOTDISS(JL,JK) + PDISSGW(JL,JK)
          ENDDO
        ENDDO
      ENDIF
      IF (LSPBS_DISSCU) THEN
        ! deep convection included
        DO JK=1,KLEV
          DO JL=KIDIA,KFDIA
            PTOTDISS(JL,JK)   = PTOTDISS(JL,JK) +    &
                 &   RATIO_BACKSCAT_CON2NUM*(PGEMU(JL))**2*ZDISSCONV_DEEP(JL,JK)
          ENDDO
        ENDDO
      ENDIF
      IF (YDDYN%LRFRIC) PTOTDISS(KIDIA:KFDIA,1:NMAXLEVRF)=0.0_JPRB
    ENDIF

!     ------------------------------------------------------------------
!      UPDATE STREAMFUNCTION (and TEMPERATURE) FORCING
!     -----------------------------------------------------------------

!     -----------------------------------------------------------------
!       CASBS
!     -----------------------------------------------------------------
    IF (LSTOPH_CASBS.AND.NSTEP >= 0) THEN
      ZCONS_AMP_CASBS=AMAGSTOPH_CASBS*RA/(REAL(NSMAX)*TSPHY)
      DO JK=1,KLEV
        DO JL=KIDIA,KFDIA
          IF (PTOTDISS_SMOOTH(JL,JK)>0._JPRB) THEN
            PSTREAM(JL,JK)=ZCONS_AMP_CASBS*(PSTOPHCA(JL))*SQRT(PTOTDISS_SMOOTH(JL,JK))
          ELSE
            PSTREAM(JL,JK)=0._JPRB     ! set -ve dissipation regions to zero dissipation
          ENDIF
        ENDDO
      ENDDO
!     -----------------------------------------------------------------
!      SPBS
!     -----------------------------------------------------------------
    ELSEIF (LSTOPH_SPBS.AND.NSTEP >= 0.AND.(MOD(NSTEP,NFRSTOPH_SPBS) == 0)) THEN
      ZCONS_AMP_SPBS=SQRT(RATIO_BACKSCAT/TSPHY) 
      DO JK=1,KLEV
        DO JL=KIDIA,KFDIA
          IF (PTOTDISS_SMOOTH(JL,JK) > 0._JPRB) THEN
            PSTREAM(JL,JK)= ZCONS_AMP_SPBS*SQRT(PTOTDISS_SMOOTH(JL,JK))*PSTREAM(JL,JK) 
          ELSE
            PSTREAM(JL,JK)=0._JPRB
          ENDIF
        ENDDO
      ENDDO
!     -----------------------------------------------------------------
!      T-backscatter
!     -----------------------------------------------------------------
      IF (LSTOPH_SPBS_T) THEN
        ZCONS_AMP_SPBS_T=  SQRT(RATIO_APE2KE*RATIO_BACKSCAT/TSPHY)
        DO JK=1,KLEV
          DO JL=KIDIA,KFDIA
            IF (PTOTDISS_SMOOTH(JL,JK) > 0._JPRB) THEN
              PTEMP(JL,JK)= ZCONS_AMP_SPBS_T*SQRT(PTOTDISS_SMOOTH(JL,JK))*PTEMP(JL,JK)
            ELSE
              PTEMP(JL,JK)=0._JPRB
            ENDIF
          ENDDO
        ENDDO
      ENDIF
    ENDIF

!     -----------------------------------------------------------------
!       second set of extrafields (output fields)
!     -----------------------------------------------------------------
    IF (LEXTRAFIELDS) THEN 
      ZSCALE=1.0_JPRB
      DO JK=1,KLEV
        DO JL=KIDIA,KFDIA
          ZDELP               =PAPRS(JL,JK)        -  PAPRS(JL,JK-1)
          ZDISS2D(JL)         =ZDISS2D(JL)         + PTOTDISS_SMOOTH(JL,JK)*ZDELP*ZRG*ZSCALE
          ZDISSNUM2D(JL)      =ZDISSNUM2D(JL)      + ZDISSNUM(JL,JK)       *ZDELP*ZRG*ZSCALE
          ZDISSGW2D(JL)       =ZDISSGW2D(JL)       + PDISSGW(JL,JK)        *ZDELP*ZRG*ZSCALE
          ZDISSCONV2D_DEEP(JL)=ZDISSCONV2D_DEEP(JL)+ ZDISSCONV_DEEP(JL,JK) *ZDELP*ZRG*ZSCALE
          ZDISSCONV2D_SHAL(JL)=ZDISSCONV2D_SHAL(JL)+ ZDISSCONV_SHAL(JL,JK) *ZDELP*ZRG*ZSCALE
          ZEGEN2D(JL)         =ZEGEN2D(JL)         + PSTREAM(JL,JK)*PVORT(JL,JK)*ZDELP*ZRG*ZSCALE

          PEXTRA(JL,JK, 1)    = PSTREAM(JL,JK)
          IF (LSTOPH_SPBS_T) THEN
            PEXTRA(JL,JK, 6)    = PTEMP(JL,JK)
          ELSE
            PEXTRA(JL,JK, 6)    = 0
          ENDIF
        ENDDO
      ENDDO
      ! 2D-fields
      DO JL=KIDIA,KFDIA
!        PEXTR2(JL,1)= PSTOPHCA(JL)
!        PEXTR2(JL,2)= PEXTR2(JL,2) + ZDISSNUM2D(JL)
!        PEXTR2(JL,3)=ZDISSGW2D(JL)/TSPHY 
!        PEXTR2(JL,4)=ZDISSCONV2D_DEEP(JL)/TSPHY 
!        PEXTR2(JL,5)= PEXTR2(JL,5) + ZEGEN2D(JL)*TSPHY
      ENDDO
    ENDIF
  ENDIF !SPBS or CASBS



!     -----------------------------------------------------------------
!     -----------------------------------------------------------------
!          2    vorticity confinement calculation
!     -----------------------------------------------------------------
!     -----------------------------------------------------------------
  IF (LVORTCON) THEN
!  zero diagnostic fields for vorticity confinement
    DO JL=KIDIA,KFDIA
      ZEGENVC2D(JL)= 0.0_JPRB
    ENDDO
    DO JK=1,KLEV
      DO JL=KIDIA,KFDIA
        ZEGENVC(JL,JK)= 0.0_JPRB    
        ZVCX(JL,JK)= 0.0_JPRB
        ZVCY(JL,JK)= 0.0_JPRB
      ENDDO
    ENDDO

  !-------------------------------------------------------------------------------------
  !  add vorticity confinement terms to the momentum equation
  !  note that (PVORTGRADY(),-PVORTGRADX())/ZVORTGRAD represents the unit vector pointing
  !  along vorticity isopleths with high vorticity to the left
  !------------------------------------------------------------------------------------
    DO JK=1,KLEV
      DO JL=KIDIA,KFDIA
        ZVORTGRAD= SQRT(PVORTGRADX(JL,JK)**2 + PVORTGRADY(JL,JK)**2) + 1.E-20_JPRB
        ZVCX(JL,JK)=   VC_CON*PVORTGRADY(JL,JK)*ABS(PVORT(JL,JK))/ZVORTGRAD
        ZVCY(JL,JK)=  -VC_CON*PVORTGRADX(JL,JK)*ABS(PVORT(JL,JK))/ZVORTGRAD
        PTENU(JL,JK)=  PTENU(JL,JK) +  ZVCX(JL,JK)
        PTENV(JL,JK)=  PTENV(JL,JK) +  ZVCY(JL,JK)
      ENDDO
    ENDDO

!     -----------------------------------------------------------------
!     extra fields from vorticity confinement
!                  -----------------------------------------------
    IF (LEXTRAFIELDS) THEN
      DO JK=1,KLEV
        DO JL=KIDIA,KFDIA
          ZEGENVC(JL,JK)= PU(JL,JK)*ZVCX(JL,JK) + PV(JL,JK)*ZVCY(JL,JK)
          ZDELP= PAPRS(JL,JK) - PAPRS(JL,JK-1)
          ZEGENVC2D(JL)= ZEGENVC2D(JL) + ZEGENVC(JL,JK)*ZDELP*ZRG
          ZVC= SQRT(ZVCX(JL,JK)**2 + ZVCY(JL,JK)**2)/8._JPRB + 1.E-20_JPRB  ! includes scaling factor for Metview
          PEXTRA(JL,JK, 3)=ZVCX(JL,JK)*1.E05_JPRB
          PEXTRA(JL,JK, 4)=ZVCY(JL,JK)*1.E05_JPRB
          PEXTRA(JL,JK, 5)= PVORT(JL,JK)
        ENDDO
      ENDDO
      DO  JL=KIDIA,KFDIA
        PEXTR2(JL,6)=  ZEGENVC2D(JL)
      ENDDO
    ENDIF
  ENDIF   !  test on LVORTCON and LEXTRAFIELDS

  END ASSOCIATE
!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SPBSGPUPD',1,ZHOOK_HANDLE)
END SUBROUTINE SPBSGPUPD
