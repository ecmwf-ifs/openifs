! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE COND_LAYER( &
 ! Input quantities
 &  YDECLDP,YDECND,YDEPHLI,YDPHY2,KDIM, PAUX, STATE, TENDENCY,&
 ! Input/Output quantities
 &  FLUX, PDIAG, & 
 ! Output tendencies
 &  TENDENCY_LOC)

!**** *COND_LAYER* - Layer routine calling cond scheme

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------

!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
! KDIM     : Derived variable for dimensions
! KEYMASK  : Mask of GFL fields usage
! state    : Derived variable for model state
! tendency : D. V. for model tendencies (entering cloud) from processes before 
! PAUX     : Derived variables for general auxiliary quantities

!     ==== Input/output ====
! FLUX     : Derived variable for fluxes
! PDIAG    : Derived variable for diagnostics

!    ==== Output tendencies from convection ====
! tendency_loc :  Output process tendencies


!        --------------------

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS

!     AUTHOR.
!     -------
!      Original : 2012-11-28  F. VANA (c) ECMWF

!     MODIFICATIONS.
!     --------------
!      F. Vana  08-Apr-2016  Small fix

!-----------------------------------------------------------------------

USE YOECLDP  , ONLY : TECLDP, NCLDQR, NCLDQS, NCLDQI, NCLDQL
USE YOECND   , ONLY : TECND
USE YOEPHLI  , ONLY : TEPHLI
USE PARKIND1 , ONLY : JPRB
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMPHYDER, ONLY : DIMENSION_TYPE, STATE_TYPE, AUX_TYPE, &
  &                   FLUX_TYPE, MASK_GFL_TYPE, AUX_DIAG_TYPE
USE YOMPHY2  , ONLY : TPHY2

!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE(TECLDP)         , INTENT(INOUT) :: YDECLDP
TYPE(TECND)          , INTENT(INOUT) :: YDECND
TYPE(TEPHLI)         , INTENT(INOUT) :: YDEPHLI
TYPE(TPHY2)          , INTENT(INOUT) :: YDPHY2
TYPE (DIMENSION_TYPE), INTENT (IN)   :: KDIM
TYPE (AUX_TYPE)      , INTENT (IN)   :: PAUX
TYPE (STATE_TYPE)    , INTENT (IN)   :: STATE
TYPE (STATE_TYPE)    , INTENT (IN)   :: TENDENCY
TYPE (FLUX_TYPE)     , INTENT(INOUT) :: FLUX
TYPE (AUX_DIAG_TYPE) , INTENT(INOUT) :: PDIAG
TYPE (STATE_TYPE)    , INTENT(INOUT) :: TENDENCY_LOC

!-----------------------------------------------------------------------

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------

#include "cond.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('COND_LAYER',0,ZHOOK_HANDLE)
ASSOCIATE(TSPHY=>YDPHY2%TSPHY)
!     ------------------------------------------------------------------

!*         1.     UNROLL THE DERIVED STRUCTURES AND CALL COND

CALL COND &
 & ( YDEPHLI, YDECND, KDIM%KIDIA , KDIM%KFDIA , KDIM%KLON  , KDIM%KLEV, &
 & TSPHY,&
 & STATE%T    , STATE%Q    , PAUX%PRS1,   PAUX%PRSF1,&
 & TENDENCY_LOC%T, TENDENCY%T, TENDENCY_LOC%Q, TENDENCY%Q, &
 & FLUX%PFHPSL, FLUX%PFHPSN, FLUX%PFPLSL, FLUX%PFPLSN          )  

!  zero other tendencis of more sophisticated cloud schemes
IF (TENDENCY_LOC%KEYMASK%A ) TENDENCY_LOC%A(:,:)         =0.0_JPRB
IF (TENDENCY_LOC%KEYMASK%QL) TENDENCY_LOC%CLD(:,:,NCLDQL)=0.0_JPRB
IF (TENDENCY_LOC%KEYMASK%QI) TENDENCY_LOC%CLD(:,:,NCLDQI)=0.0_JPRB
IF (TENDENCY_LOC%KEYMASK%QR) TENDENCY_LOC%CLD(:,:,NCLDQR)=0.0_JPRB
IF (TENDENCY_LOC%KEYMASK%QS) TENDENCY_LOC%CLD(:,:,NCLDQS)=0.0_JPRB

PDIAG%PPRECTYPE(KDIM%KIDIA:KDIM%KFDIA) = 0._JPRB
PDIAG%PFZRA    (KDIM%KIDIA:KDIM%KFDIA) = 0._JPRB
!     ------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('COND_LAYER',1,ZHOOK_HANDLE)
END SUBROUTINE COND_LAYER
