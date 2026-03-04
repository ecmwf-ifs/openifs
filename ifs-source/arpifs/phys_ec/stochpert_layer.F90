! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE STOCHPERT_LAYER( &
 ! Input quantities
 &  YDMODEL,YGFL,YDPHY2,KDIM, PAUX, STATE, TENDENCY_DYN, &
 &  TENDENCY_PHY, PTENGFL, &
 ! Input/Output quantities
 &  PPERT, TENDENCY, GEMSL )

!**** *STOCHPERT_LAYER* - Layer routine calling stochatic perturbation to physics

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------

!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
! KDIM         : Derived variable for dimensions
! PAUX         : Derived variables for general auxiliary quantities
! STATE        : Derived variable for (updated) model state
! TENDENCY_DYN : D. V. for model tendencies from explicit dynamics
! TENDENCY_PHY : D. V. for model tendencies from physics parametrisations<
!     ==== Input/output ====
! PPERT        :  Perturbations, tendencies for SPPTGFIX
! TENDENCY     :  Cumulated tendencies
! GEMSL        :  Local GEMS arrays

!        --------------------

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     Note: Physics / dynamics tendencies come IN;
!           perturbed total tendencies come OUT.

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS

!     AUTHOR.
!     -------
!      Original : 2012-11-29  F. VANA (c) ECMWF

!     MODIFICATIONS.
!     --------------
!      Dec2015 S.Massart: perturbation of GFL tendencies
!      Dec2015 A.Inness: Add tracers
!      S. Lang, Jan 2016 - pass PPERT2 to SPPT for SPPTGFIX
!      SJ Lock: Jan-2016: Enabled multiple patterns for iSPPT
!      SJ Lock: Oct-2016: Option to leave clear-skies radiation UNperturbed
!-----------------------------------------------------------------------

USE PARKIND1  , ONLY : JPIM ,   JPRB
USE YOMHOOK   , ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOMPHYDER , ONLY : DIMENSION_TYPE, STATE_TYPE, AUX_TYPE, &
  &                    PERTURB_TYPE, GEMS_LOCAL_TYPE
USE YOMPHY2   , ONLY : TPHY2
USE YOM_YGFL  , ONLY : TYPE_GFLD
USE TYPE_MODEL, ONLY : MODEL

!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE(MODEL)                    , INTENT(INOUT) :: YDMODEL
TYPE(TPHY2)                    , INTENT(INOUT) :: YDPHY2
TYPE(TYPE_GFLD)                , INTENT(INOUT) :: YGFL
TYPE (DIMENSION_TYPE)          , INTENT (IN)   :: KDIM
TYPE (AUX_TYPE)                , INTENT (IN)   :: PAUX
TYPE (STATE_TYPE)              , INTENT (IN)   :: STATE
TYPE (STATE_TYPE)              , INTENT (IN)   :: TENDENCY_DYN
TYPE (STATE_TYPE)              , INTENT (IN)   :: TENDENCY_PHY(KDIM%K2DSDT)
REAL(KIND=JPRB)                , INTENT (IN)   :: PTENGFL(KDIM%KLON,KDIM%KLEV,YGFL%NDIM1)
TYPE (PERTURB_TYPE)            , INTENT (INOUT):: PPERT
TYPE (GEMS_LOCAL_TYPE)         , INTENT (INOUT):: GEMSL
TYPE (STATE_TYPE)              , INTENT (INOUT):: TENDENCY

!-----------------------------------------------------------------------

INTEGER(KIND=JPIM) :: J2D
REAL(KIND=JPRB)    :: ZMULNOISE(KDIM%KLON,KDIM%K2DSDT)
REAL(KIND=JPRB)    :: ZUNP_U(KDIM%KLON,KDIM%KLEV)        !For unperturbed U tendencies
REAL(KIND=JPRB)    :: ZUNP_V(KDIM%KLON,KDIM%KLEV)        ! "       "      V     "
REAL(KIND=JPRB)    :: ZUNP_T(KDIM%KLON,KDIM%KLEV)        ! "       "      T     "
REAL(KIND=JPRB)    :: ZUNP_Q(KDIM%KLON,KDIM%KLEV)        ! "       "      Q     "
REAL(KIND=JPRB)    :: ZPHY_U(KDIM%KLON,KDIM%KLEV,KDIM%K2DSDT) !For   perturbed U tendencies
REAL(KIND=JPRB)    :: ZPHY_V(KDIM%KLON,KDIM%KLEV,KDIM%K2DSDT) ! "       "      V     "
REAL(KIND=JPRB)    :: ZPHY_T(KDIM%KLON,KDIM%KLEV,KDIM%K2DSDT) ! "       "      T     "
REAL(KIND=JPRB)    :: ZPHY_Q(KDIM%KLON,KDIM%KLEV,KDIM%K2DSDT) ! "       "      Q     "
LOGICAL            :: LLSPPTMOD

REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

!-----------------------------------------------------------------------

#include "sppten.intfb.h"
#include "abor1.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('STOCHPERT_LAYER',0,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

! Prepare fields for SPPT:
IF (YDMODEL%YRML_GCONF%YRSPPT_CONFIG%LSPSDT) THEN

  !*         0.     INITIALISE FIELDS

  ZUNP_U(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV) = 0.0_JPRB
  ZUNP_V(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV) = 0.0_JPRB
  ZUNP_T(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV) = 0.0_JPRB
  ZUNP_Q(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV) = 0.0_JPRB

  LLSPPTMOD = .FALSE.  !Modified SPPT
  IF (.NOT.(YDMODEL%YRML_GCONF%YRSPPT_CONFIG%LRADCLR_SDT) .OR. &
    & .NOT.(YDMODEL%YRML_GCONF%YRSPPT_CONFIG%LSATADJ_SDT) .OR. KDIM%K2DSDT>1) LLSPPTMOD = .TRUE.

  !Modified SPPT: account for unperturbed physics tendencies
  IF (LLSPPTMOD) THEN
    !Physics NOT to be perturbed: total minus dynamics minus physics to be perturbed (more see below...)
    ZUNP_U(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV) = TENDENCY%U(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV) &
                                            & - TENDENCY_DYN%U(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV)
    ZUNP_V(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV) = TENDENCY%V(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV) &
                                            & - TENDENCY_DYN%V(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV)
    ZUNP_T(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV) = TENDENCY%T(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV) &
                                            & - TENDENCY_DYN%T(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV)
    ZUNP_Q(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV) = TENDENCY%Q(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV) &
                                            & - TENDENCY_DYN%Q(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV)
  ENDIF

  !*         1.     UNROLL THE DERIVED STRUCTURES AND CALL APPROPRIATE ROUTINE

  DO J2D=1,KDIM%K2DSDT   !Standard SPPT: K2DSDT=1
    !Random pattern field(s)
    ZMULNOISE(KDIM%KIDIA:KDIM%KFDIA,J2D)          = PPERT%PGP2DSDT(KDIM%KIDIA:KDIM%KFDIA,1,J2D)
    !Physics to be perturbed: from callpar
    ZPHY_U(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV,J2D) = TENDENCY_PHY(J2D)%U(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV)
    ZPHY_V(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV,J2D) = TENDENCY_PHY(J2D)%V(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV)
    ZPHY_T(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV,J2D) = TENDENCY_PHY(J2D)%T(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV)
    ZPHY_Q(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV,J2D) = TENDENCY_PHY(J2D)%Q(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV)

    IF (LLSPPTMOD) THEN !Modified SPPT (iSPPT or LRADCLR_SDT or LSATADJ_SDT): account for unperturbed physics tendencies
      !Physics tendencies NOT to be perturbed: total physics (from above) minus physics to be perturbed
      ZUNP_U(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV) = ZUNP_U(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV) &
                                              & - ZPHY_U(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV,J2D)
      ZUNP_V(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV) = ZUNP_V(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV) &
                                              & - ZPHY_V(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV,J2D)
      ZUNP_T(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV) = ZUNP_T(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV) &
                                              & - ZPHY_T(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV,J2D)
      ZUNP_Q(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV) = ZUNP_Q(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV) &
                                              & - ZPHY_Q(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV,J2D)
    ENDIF
  ENDDO !J2D

  CALL SPPTEN(YDMODEL,YGFL,KDIM%KIDIA,KDIM%KFDIA,KDIM%KLON,KDIM%KLEV,KDIM%K2DSDT,YDPHY2%TSPHY,   &
   &          PTSL=STATE%T,PQSL=STATE%Q,PA=STATE%A,PAP=PAUX%PRSF1,PAPH=PAUX%PRS1,                &
   &          PDYN_U=TENDENCY_DYN%U,PDYN_V=TENDENCY_DYN%V,PDYN_T=TENDENCY_DYN%T,PDYN_Q=TENDENCY_DYN%Q,&
   &          PUNP_U=ZUNP_U,PUNP_V=ZUNP_V,PUNP_T=ZUNP_T,PUNP_Q=ZUNP_Q,                           &
   &          PPHY_U=ZPHY_U,PPHY_V=ZPHY_V,PPHY_T=ZPHY_T,PPHY_Q=ZPHY_Q,                           &
   &          PMULNOISE=ZMULNOISE,   &
   &          PTENU=TENDENCY%U,PTENV=TENDENCY%V,PTENT=TENDENCY%T,PTENQ=TENDENCY%Q,               &
   &          PTENGFL=PTENGFL,KTRAC=GEMSL%ITRAC,PTENC=GEMSL%ZTENC,PPERT=PPERT)
ELSE
  CALL ABOR1('Should only get here if LSPSDT!')
ENDIF

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('STOCHPERT_LAYER',1,ZHOOK_HANDLE)
END SUBROUTINE STOCHPERT_LAYER
