! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SURFWS_LAYER(YDSURF,YDEPHY,KDIM, PSURF,SURFL,PAUX)

!**** *SURFWS_LAYER* - Layer routine calling SURFWS scheme 

!     PURPOSE.
!     --------
!     Layer routine for parametrizing multi-layer snow fields
!     based on single-layer snow input.

!**   INTERFACE.
!     ----------
!       SURFWS_LAYER is called by *CALLPAR*

!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
! KDIM     : Derived variable for dimensions
! PAUX     : Derived variables for general auxiliary quantities

!     ==== Input/output ====
! PSURF    : Derived variable for general surface quantities
! SURFL    : Derived variable for local surface quantities


!        IMPLICIT ARGUMENTS :   NONE
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
!      Original : 2018 Jan  G.Arduini (c) ECMWF

!     MODIFICATIONS.
!     --------------

!-----------------------------------------------------------------------

USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE PARKIND1  , ONLY : JPRB,JPIM
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOEPHY    , ONLY : TEPHY
USE YOMPHYDER , ONLY : DIMENSION_TYPE, &
   & SURF_AND_MORE_TYPE, SURF_AND_MORE_LOCAL_TYPE,AUX_TYPE

!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE(TSURF), INTENT(INOUT) :: YDSURF
TYPE(TEPHY), INTENT(INOUT) :: YDEPHY
TYPE (DIMENSION_TYPE)          , INTENT (IN)   :: KDIM
TYPE (AUX_TYPE)                , INTENT (IN)   :: PAUX
TYPE (SURF_AND_MORE_TYPE)      , INTENT(INOUT) :: PSURF
TYPE (SURF_AND_MORE_LOCAL_TYPE), INTENT(INOUT) :: SURFL
!-----------------------------------------------------------------------
REAL(KIND=JPRB) :: ZSNM1M(KDIM%KLON, KDIM%KLEVSN), ZTSNM1M(KDIM%KLON, KDIM%KLEVSN), &
                 & ZRSNM1M(KDIM%KLON, KDIM%KLEVSN), ZWSNM1M(KDIM%KLON, KDIM%KLEVSN)

!REAL(KIND=JPRB) :: ZSNM1M(KDIM%KLON, 5), ZTSNM1M(KDIM%KLON, 5), ZRSNM1M(KDIM%KLON, 5), ZWSNM1M(KDIM%KLON, 5)

INTEGER(KIND=JPIM) :: JL

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------

#include "surfws.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SURFWS_LAYER',0,ZHOOK_HANDLE)
ASSOCIATE(YSURF=>YDEPHY%YSURF, &
 & YSD_VF=>YDSURF%YSD_VF, YSP_SL=>YDSURF%YSP_SL, &
 & YSP_SB=>YDSURF%YSP_SB, YSP_RR=>YDSURF%YSP_RR, &
 & YSP_SG=>YDSURF%YSP_SG)

!     ------------------------------------------------------------------

DO JL=KDIM%KIDIA,KDIM%KFDIA
  ZSNM1M(JL,  1:KDIM%KLEVSN) = SURFL%ZSNM1M(JL, 1:KDIM%KLEVSN)
  ZRSNM1M(JL, 1:KDIM%KLEVSN) = PSURF%PSP_SG(JL, 1, YSP_SG%YR%MP9)
  ZTSNM1M(JL, 1:KDIM%KLEVSN) = PSURF%PSP_SG(JL, 1, YSP_SG%YT%MP9)
  ZWSNM1M(JL, 1:KDIM%KLEVSN) = 0._JPRB

ENDDO

!*         1.     UNROLL THE DERIVED STRUCTURES AND CALL SURFWS

CALL SURFWS (YDSURF=YSURF,KIDIA=KDIM%KIDIA,KFDIA=KDIM%KFDIA,KLON=KDIM%KLON,KLEVS=KDIM%KLEVS,KLEVSN=KDIM%KLEVSN,&
    & PSDOR=PSURF%PHSTD, &
    & KTILES=KDIM%KTILES, PLSM=PSURF%PSD_VF(:,YSD_VF%YLSM%MP), PFRTI=SURFL%ZFRTI, PMU0=PAUX%PMU0, & 
    & PTSAM1M=PSURF%PSP_SB(:,:,YSP_SB%YT%MP9), PTSKIN=PSURF%PSP_RR(:,YSP_RR%YT%MP9), &
    & PALBSN=PSURF%PSP_SG(:,1,YSP_SG%YA%MP9), PTSNM1M=ZTSNM1M, PSNM1M=ZSNM1M, PRSNM1M=ZRSNM1M, PWSNM1M= ZWSNM1M)



! Update generalized variables PSP_SG at the t-1 step
DO JL=KDIM%KIDIA,KDIM%KFDIA
  SURFL%ZSNM1M(JL, 1:KDIM%KLEVSN)                = ZSNM1M(JL,  1:KDIM%KLEVSN) 
  PSURF%PSP_SG(JL, 1:KDIM%KLEVSN, YSP_SG%YF%MP9) = ZSNM1M(JL,  1:KDIM%KLEVSN)
  PSURF%PSP_SG(JL, 1:KDIM%KLEVSN, YSP_SG%YR%MP9) = ZRSNM1M(JL, 1:KDIM%KLEVSN) 
  PSURF%PSP_SG(JL, 1:KDIM%KLEVSN, YSP_SG%YT%MP9) = ZTSNM1M(JL, 1:KDIM%KLEVSN)
  PSURF%PSP_SG(JL, 1:KDIM%KLEVSN, YSP_SG%YW%MP9) = ZWSNM1M(JL, 1:KDIM%KLEVSN) 

! Update also the t+1 step 
  PSURF%PSP_SG(JL, 1:KDIM%KLEVSN, YSP_SG%YF%MP1) = ZSNM1M(JL,  1:KDIM%KLEVSN)
  PSURF%PSP_SG(JL, 1:KDIM%KLEVSN, YSP_SG%YR%MP1) = ZRSNM1M(JL, 1:KDIM%KLEVSN)
  PSURF%PSP_SG(JL, 1:KDIM%KLEVSN, YSP_SG%YT%MP1) = ZTSNM1M(JL, 1:KDIM%KLEVSN)
  PSURF%PSP_SG(JL, 1:KDIM%KLEVSN, YSP_SG%YW%MP1) = ZWSNM1M(JL, 1:KDIM%KLEVSN)
ENDDO


!     ------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SURFWS_LAYER',1,ZHOOK_HANDLE)
END SUBROUTINE SURFWS_LAYER
