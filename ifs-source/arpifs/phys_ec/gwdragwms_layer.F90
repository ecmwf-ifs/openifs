! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE GWDRAGWMS_LAYER(&
 ! Input quantities
 &  YDSTA, YDEGWD,YDEGWWMS,YGFL,KDIM, STATE, PAUX,&
 ! Input/Output quantities
 & FLUX, NOGW)

!**** *GWDRAGWMS_LAYER* - Routine computing the full non orographic gravity wave drag

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------

!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
! KDIM     : Derived variable for dimensions
! state    : Derived variable for updated model state
! PAUX     : Derived variables for general auxiliary quantities

!     ==== Input/output ====
! FLUX     : Derived variable for fluxes
! PGFL     : GFL atructure

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
!      Original : 2012-11-27  F. VANA (c) ECMWF

!     MODIFICATIONS.
!     --------------

!-----------------------------------------------------------------------

USE YOEGWD   , ONLY : TEGWD
USE YOMSTA   , ONLY : TSTA
USE PARKIND1 , ONLY : JPIM, JPRB
USE VARIABLE_MODULE , ONLY: VARIABLE_3D
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMPHYDER, ONLY : DIMENSION_TYPE, STATE_TYPE, AUX_TYPE, FLUX_TYPE
USE YOM_YGFL , ONLY : TYPE_GFLD
USE YOEGWWMS , ONLY : TEGWWMS

!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE(TSTA)            , INTENT (IN)   :: YDSTA
TYPE(TEGWD)            ,INTENT(INOUT) :: YDEGWD
TYPE(TEGWWMS)          ,INTENT(INOUT) :: YDEGWWMS
TYPE(TYPE_GFLD)        ,INTENT(INOUT) :: YGFL
TYPE (DIMENSION_TYPE) , INTENT (IN)   :: KDIM
TYPE (STATE_TYPE)     , INTENT (IN)   :: STATE
TYPE (AUX_TYPE)       , INTENT (IN)   :: PAUX
TYPE (FLUX_TYPE)      , INTENT(INOUT) :: FLUX
TYPE (VARIABLE_3D)    , INTENT(INOUT) :: NOGW(:)

!-----------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JK, JL, JLAUNCH
REAL(KIND=JPRB) :: ZRAINT(KDIM%KLON)
REAL(KIND=JPRB) :: ZFLXWGWU(KDIM%KLON,0:KDIM%KLEV), ZFLXWGWV(KDIM%KLON,0:KDIM%KLEV)
REAL(KIND=JPRB) :: ZNOGWU(KDIM%KLON,KDIM%KLEV), ZNOGWV(KDIM%KLON,KDIM%KLEV)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------

#include "gwdrag_wms.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GWDRAGWMS_LAYER',0,ZHOOK_HANDLE)
ASSOCIATE(NDIM=>YGFL%NDIM, YNOGW=>YGFL%YNOGW, &
 & GTPHYGWWMS=>YDEGWWMS%GTPHYGWWMS, NLAUNCHLEV=>YDEGWWMS%NLAUNCHLEV)
!     ------------------------------------------------------------------

!*         0.     INITIALIZATION

ZRAINT(KDIM%KIDIA:KDIM%KFDIA)=MAX(0.0_JPRB,&
  &  FLUX%PFPLCL(KDIM%KIDIA:KDIM%KFDIA,KDIM%KLEV)+FLUX%PFPLCN(KDIM%KIDIA:KDIM%KFDIA,KDIM%KLEV)+&
  &  FLUX%PFPLSL(KDIM%KIDIA:KDIM%KFDIA,KDIM%KLEV)+FLUX%PFPLSN(KDIM%KIDIA:KDIM%KFDIA,KDIM%KLEV))

DO JK=1,KDIM%KLEV
  DO JL=KDIM%KIDIA,KDIM%KFDIA
    NOGW(1)%P(JL,JK)=0.0_JPRB
    NOGW(2)%P(JL,JK)=0.0_JPRB
  ENDDO
ENDDO

!*         1.     UNROLL THE DERIVED STRUCTURES AND CALL GWDRAG_WMS

DO JLAUNCH=1,NLAUNCHLEV

  CALL GWDRAG_WMS(YDSTA, YDEGWD,YDEGWWMS, KDIM%KIDIA,  KDIM%KFDIA, KDIM%KLON  , KDIM%KLEV , JLAUNCH, GTPHYGWWMS,&
    & STATE%T, STATE%U, STATE%V, PAUX%PAPRSF, PAUX%PAPRS, PAUX%PAPHIF,&
    & PAUX%PGELAT, PAUX%PGAW, ZRAINT, ZNOGWU, ZNOGWV, ZFLXWGWU, ZFLXWGWV )

  DO JK=1,KDIM%KLEV
    DO JL=KDIM%KIDIA,KDIM%KFDIA
      NOGW(1)%P(JL,JK) = NOGW(1)%P(JL,JK) + ZNOGWU(JL,JK)
      NOGW(2)%P(JL,JK) = NOGW(2)%P(JL,JK) + ZNOGWV(JL,JK)
    ENDDO
  ENDDO

!*         2.     ADD CONTRIBUTION TO GWD DRAG 

  DO JK=0,KDIM%KLEV
    DO JL=KDIM%KIDIA,KDIM%KFDIA
      FLUX%PSTRDU(JL,JK)=FLUX%PSTRDU(JL,JK)-ZFLXWGWU(JL,JK)
      FLUX%PSTRDV(JL,JK)=FLUX%PSTRDV(JL,JK)-ZFLXWGWV(JL,JK)
    ENDDO
  ENDDO

ENDDO

!     ------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GWDRAGWMS_LAYER',1,ZHOOK_HANDLE)
END SUBROUTINE GWDRAGWMS_LAYER
