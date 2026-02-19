! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE NEMOADDFLDS_LAYER(YDSURF, &
 ! Input quantities
  & YDMCC,KDIM, SURFL, &
 ! Input/Output quantities
  & PSURF)

!**** *NEMOADDFLDS_LAYER* - Layer routine for flouxes to be coupled with OPA/LIM

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------

!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
! KDIM     : Derived variable for dimensions
! SURFL    : Derived variable for local surface quantities 

!     ==== Input/output ====
! PSURF    : Derived variables for general surface quantities


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
!      Original : 10-Apr-2013 F. VANA (c) ECMWF

!     MODIFICATIONS.
!     --------------
!     P. Marguinaud : 04-10-2016 : Port to single precision

!-----------------------------------------------------------------------

USE YOMMCC             , ONLY : TMCC
USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE PARKIND1           , ONLY : JPIM,    JPRB
USE YOMHOOK            , ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMPHYDER          , ONLY : DIMENSION_TYPE, SURF_AND_MORE_TYPE, &
   &                            SURF_AND_MORE_LOCAL_TYPE

!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE(TSURF), INTENT(INOUT) :: YDSURF
TYPE(TMCC)  ,INTENT(INOUT) :: YDMCC
TYPE (DIMENSION_TYPE)          , INTENT (IN)   :: KDIM
TYPE (SURF_AND_MORE_LOCAL_TYPE), INTENT (IN)   :: SURFL
TYPE (SURF_AND_MORE_TYPE)      , INTENT(INOUT) :: PSURF
!-----------------------------------------------------------------------
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------

#include "nemoaddflds.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('NEMOADDFLDS_LAYER',0,ZHOOK_HANDLE)
ASSOCIATE(YSD_VF=>YDSURF%YSD_VF)

!*         1.  CALL NEMOADDFLDS
#ifndef PARKIND1_SINGLE
CALL NEMOADDFLDS(YDMCC,KDIM%KSTGLO,KDIM%KIDIA,KDIM%KFDIA, &
  & PSURF%PSD_VF(KDIM%KIDIA:KDIM%KFDIA,YSD_VF%YSST%MP),PSURF%PTSKTI(KDIM%KIDIA:KDIM%KFDIA,1))
#endif


END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('NEMOADDFLDS_LAYER',1,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

END SUBROUTINE NEMOADDFLDS_LAYER
