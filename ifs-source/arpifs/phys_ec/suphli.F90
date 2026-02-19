! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

!OPTIONS XOPT(HSFUN)
SUBROUTINE SUPHLI(YDSTA,YDDIMV,YDMODEL)

!     ------------------------------------------------------------------

!**   *SUPHLI* IS THE SET-UP ROUTINE FOR COMMON BLOCK *YOEPHLI*

!     J.F. MAHFOUF         E.C.M.W.F.      96/06/23

!     PURPOSE
!     -------

!          THIS ROUTINE INITIALIZES THE CONSTANTS IN COMMON BLOCK
!     *YOEPHLI*

!     INTERFACE.
!     ----------

!     CALL *SUPHLI* FROM *SUPHEC*

!     METHOD.
!     -------

!         INITIALIZATION OF THE CONSTANTS USED IN THE LINEARIZED
!         PHYSICS

!     EXTERNALS.
!     ----------

!        NONE

!     REFERENCE.
!     ----------

!     MODIFICATIONS
!     -------------
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        M.Janiskova   23-Nov-2006 switch for transmission functions
!                                  for H2O + CO2
!        M.Janiskova   18-Mar-2008 switch for optimized TL/AD of SW radiation
!                      01-May-2008 frequency of calling linearized LW
!        H.Hersbach    04-Dec-2009 trigger LEPPCFLS on LVDFTRAJ
!        T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!----------------------------------------------------------------------

USE TYPE_MODEL , ONLY : MODEL
USE YOMSTA     , ONLY : TSTA
USE YOMDIMV    , ONLY : TDIMV
USE PARKIND1   , ONLY : JPIM, JPRB
USE YOMHOOK    , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOETHF     , ONLY : RTWAT, RTICE
USE YOMLUN     , ONLY : NULNAM
USE YOMOBS     , ONLY : LVDFTRAJ

IMPLICIT NONE


TYPE(TSTA)  ,INTENT(IN)          :: YDSTA
TYPE(TDIMV) ,INTENT(IN)          :: YDDIMV
TYPE(MODEL) ,INTENT(INOUT),TARGET:: YDMODEL
INTEGER(KIND=JPIM) :: JLEV
LOGICAL :: LLEVLWC

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "posnam.intfb.h"

LOGICAL, POINTER :: LPHYLIN,LPHYSFCLIN, LTLEVOL, LOPPTWINS

#include "naephli.nam.h"

!     ------------------------------------------------------------------

!*         1.     SET LOGICAL TO SWICH ON LINEARIZED PHYSICS
!                 ------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUPHLI',0,ZHOOK_HANDLE)
ASSOCIATE(YDEPHLI=>YDMODEL%YRML_PHY_SLIN%YREPHLI,YDPHNC=>YDMODEL%YRML_PHY_SLIN%YRPHNC, &
 & YDERAD=>YDMODEL%YRML_PHY_RAD%YRERAD,YDELWRAD=>YDMODEL%YRML_PHY_RAD%YRELWRAD,YDEPHY=>YDMODEL%YRML_PHY_EC%YREPHY)
ASSOCIATE(NFLEVG=>YDDIMV%NFLEVG, &
 & NLEVLWC=>YDELWRAD%NLEVLWC, &
 & LEPPCFLS=>YDEPHLI%LEPPCFLS, LRAISANEN=>YDEPHLI%LRAISANEN, &
 & RLPAL1=>YDEPHLI%RLPAL1, RLPAL2=>YDEPHLI%RLPAL2, RLPBB=>YDEPHLI%RLPBB, &
 & RLPBETA=>YDEPHLI%RLPBETA, RLPCC=>YDEPHLI%RLPCC, RLPDD=>YDEPHLI%RLPDD, &
 & RLPDRAG=>YDEPHLI%RLPDRAG, RLPEVAP=>YDEPHLI%RLPEVAP, RLPMIXL=>YDEPHLI%RLPMIXL, &
 & RLPP00=>YDEPHLI%RLPP00, RLPTRC=>YDEPHLI%RLPTRC, &
 & LEMWAVE=>YDEPHY%LEMWAVE, &
 & LRRTM=>YDERAD%LRRTM, NOVLP=>YDERAD%NOVLP, &
 & LERADLW2=>YDPHNC%LERADLW2, LH2OCO2=>YDPHNC%LH2OCO2, LWLOPT=>YDPHNC%LWLOPT, &
 & LWSOPT=>YDPHNC%LWSOPT, LERADI2=>YDPHNC%LERADI2, &
 & STZ=>YDSTA%STZ)

! Associate pointers for variables in namelist
LPHYLIN   => YDEPHLI%LPHYLIN
LPHYSFCLIN=> YDEPHLI%LPHYSFCLIN
LTLEVOL   => YDEPHLI%LTLEVOL
LOPPTWINS => YDEPHLI%LOPPTWINS

LPHYLIN   = .FALSE.
LPHYSFCLIN= .FALSE.
LTLEVOL   = .FALSE.
LOPPTWINS = .FALSE.

CALL POSNAM(NULNAM,'NAEPHLI')
READ(NULNAM,NAEPHLI)

!*         1.1 Processing of surface fields activated
!          ------------------------------------------

LEPPCFLS = .FALSE.
IF (LEMWAVE) LEPPCFLS = .TRUE.   ! activate 10m/2m-surface field calculation in 4D-Var rain assimilation
IF (LVDFTRAJ)LEPPCFLS = .TRUE.   ! activate 10m/2m-surface field calculation for scatt/conventional data

!*         2.     SET CONSTANTS RELATED TO WATER MIXED PHASE
!                 ------------------------------------------

RLPTRC=RTICE+(RTWAT-RTICE)/SQRT(2.0_JPRB)
RLPAL1=0.15_JPRB
RLPAL2=20._JPRB

!*         3.     SET CONSTANTS RELATED TO VERTICAL DIFFUSION
!                 -------------------------------------------

!     CONSTANTS OF THE LOUIS FORMULATION 

RLPBB=5._JPRB
RLPCC=5._JPRB
RLPDD=5._JPRB

!     PSEUDO DEPTH OF THE BOUNDARY LAYER

RLPMIXL=4000._JPRB

!     REDUCTION FACTOR OF THE ASYMPTOTIC MIXING LENGTH

RLPBETA=0.2_JPRB

!*         4.     SET CONSTANTS RELATED TO GRAVITY WAVE DRAG
!                 ------------------------------------------

!RLPDRAG=0.3_JPRB
!effectively switches off TL of gwdrag, should be on but bug somewhere, Nils
RLPDRAG=0._JPRB

!*         5.     SET CONSTANTS RELATED TO RAINFALL EVAPORATION 
!                 ---------------------------------------------

RLPEVAP=0.0_JPRB

!*         6.     SET CONSTANTS RELATED TO RADIATION
!                 ----------------------------------

!     Pressure level above which long-wave cooling is not applied

RLPP00=30000._JPRB

!     Using Raisanen overlap scheme

IF (LRRTM .AND. (NOVLP == 1)) THEN
  LRAISANEN = .TRUE.
ELSE
  LRAISANEN = .FALSE.
ENDIF

!*         7.     SET UP FOR OPTIMIZED TL/AD RADIATION
!                 ------------------------------------

!     Transmission functions computed only for H2O and CO2 absorbers
!     in LW radiation

IF(LERADI2) THEN
  LH2OCO2 = .TRUE.

!     Optimized version for TL/AD of LW radiation

  LWLOPT = .TRUE.

!     Optimized version for TL/AD of SW radiation

  LWSOPT = .TRUE.
ELSE
  LH2OCO2=.FALSE.
  LWLOPT=.FALSE.
  LWSOPT=.FALSE.
ENDIF

!     Level hight up to which cloud effects on LW radiation are accounted
!     for (based on the standard atmosphere height)

NLEVLWC = 5
LLEVLWC = .FALSE.

IF (LERADLW2) THEN
  DO JLEV=1,NFLEVG 
    IF (.NOT. LLEVLWC) THEN
      IF (STZ(JLEV) < 16600.0_JPRB) THEN
        NLEVLWC = JLEV
        LLEVLWC = .TRUE.
      ENDIF
    ENDIF
  ENDDO 
ENDIF

!     -------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUPHLI',1,ZHOOK_HANDLE)
END SUBROUTINE SUPHLI
