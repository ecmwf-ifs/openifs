! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE AER_GLOMAPDIAG_LAYER(KDIM, PTSPHY, PAUX, STATE, PDIAG,GEMSL)
use parkind1 , only:&
 & jprb
USE YOMPHYDER ,ONLY : DIMENSION_TYPE, STATE_TYPE, AUX_TYPE, &
 & AUX_DIAG_TYPE, &
 & GEMS_LOCAL_TYPE
TYPE (DIMENSION_TYPE) , INTENT (IN) :: KDIM
REAL(KIND=JPRB) ,INTENT(IN) :: PTSPHY
TYPE (AUX_TYPE) , INTENT (IN) :: PAUX
TYPE (STATE_TYPE) , INTENT (IN) :: STATE
TYPE (AUX_DIAG_TYPE) , INTENT (IN) :: PDIAG
TYPE (GEMS_LOCAL_TYPE) , INTENT(INOUT) :: GEMSL
call abor1("aer_glomapdiag_layer.F90 should never be called with OpenIFS - EXITING")
END SUBROUTINE AER_GLOMAPDIAG_LAYER
