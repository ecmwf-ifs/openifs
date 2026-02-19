! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SURFTSTP_S_LAYER(YDSURF,YDEPHY,YDPHY2,PAUX,KDIM,PSURF,SURFL,LLKEYS,FLUX)

!**** *SURFTSTP_S_LAYER* - Layer routine calling time stepping of simplified surface scheme

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------

!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
! KDIM     : Derived variable for dimensions

!     ==== Input/output ====
! PSURF    : Derived variables for general surface quantities
! SURFL    : Derived variables for local surface quantities
! LLKEYS   : Derived variable with keys
! FLUX     : Derived variable for fluxes
! PAUX     : Derived variables for general auxiliary quantities



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
!      Original : 2013-04-12 F. VANA (c) ECMWF

!     MODIFICATIONS.
!     --------------
!     E. Dutra/G.Arduini : Jan 2018: Temporary fix for multi-layer snow (should be revised)
!                           PSURF%PTSNE1(:,1)

!-----------------------------------------------------------------------

USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE PARKIND1 , ONLY : JPRB, JPIM
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOMPHYDER, ONLY : DIMENSION_TYPE,AUX_TYPE, SURF_AND_MORE_TYPE, SURF_AND_MORE_LOCAL_TYPE, &
   &                  KEYS_LOCAL_TYPE, FLUX_TYPE
USE YOMPHY2  , ONLY : TPHY2
USE YOEPHY   , ONLY : TEPHY


!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE(TSURF), INTENT(INOUT) :: YDSURF
TYPE(TEPHY) ,INTENT(INOUT) :: YDEPHY
TYPE(TPHY2) ,INTENT(INOUT) :: YDPHY2
TYPE (DIMENSION_TYPE)          , INTENT (IN)   :: KDIM
TYPE (AUX_TYPE)                , INTENT (IN)   :: PAUX
TYPE (SURF_AND_MORE_TYPE)      , INTENT(INOUT) :: PSURF
TYPE (SURF_AND_MORE_LOCAL_TYPE), INTENT(INOUT) :: SURFL
TYPE (KEYS_LOCAL_TYPE)         , INTENT(INOUT) :: LLKEYS
TYPE (FLUX_TYPE)               , INTENT(INOUT) :: FLUX

!-----------------------------------------------------------------------
! Local var
REAL(KIND=JPRB) :: ZTMP, ZTMP2
INTEGER(KIND=JPIM) :: JK, JL

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------

#include "surftstps.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SURFTSTP_S_LAYER',0,ZHOOK_HANDLE)

!     ------------------------------------------------------------------
!*         1.     UNROLL THE DERIVED STRUCTURES AND CALL SURFTSTPS

CALL SURFTSTPS(YDSURF=YDEPHY%YSURF, KIDIA=KDIM%KIDIA, KFDIA=KDIM%KFDIA, &
  & KLON=KDIM%KLON, KLEVS=KDIM%KLEVS, KLEVSN=YDSURF%YSP_SGD%NLEVS, KTILES=KDIM%KTILES,&
  & KSOTY=PSURF%ISOTY,&
  & PTSPHY=YDPHY2%TSPHY , PSDOR=PSURF%PHSTD, PFRTI=SURFL%ZFRTI,&
  & PAHFSTI=PSURF%PAHFSTI, PEVAPTI=PSURF%PEVAPTI, PSSRFLTI=SURFL%ZFRSOTI,&
  & LDLAND=LLKEYS%LLLAND, LDSICE=LLKEYS%LLSICE, LDSI=YDEPHY%LESICE, LDNH=LLKEYS%LLNH,&
  & PSNM1M=SURFL%ZSNM1M,PTSNM1M=PSURF%PSP_SG(:,:,YDSURF%YSP_SG%YT%MP9), &
  & PRSNM1M=PSURF%PSP_SG(:,:,YDSURF%YSP_SG%YR%MP9), &
  & PTSAM1M=PSURF%PSP_SB(:,:,YDSURF%YSP_SB%YT%MP9), PTIAM1M=PSURF%PSP_SB(:,:,YDSURF%YSP_SB%YTL%MP9),&
  & PWLM1M=SURFL%ZWLM1M, PWSAM1M=PSURF%PSP_SB(:,:,YDSURF%YSP_SB%YQ%MP9),&
  & PHLICEM1M=PSURF%PSP_SL(:,YDSURF%YSP_SL%YLICD%MP), PAPRS=PAUX%PAPRS(:,KDIM%KLEV),&
  & PRSFC=FLUX%PFPLCL(:,KDIM%KLEV),PRSFL=FLUX%PFPLSL(:,KDIM%KLEV),&
  & PSLRFL=FLUX%PFRTH(:,KDIM%KLEV), PSSFC=FLUX%PFPLCN(:,KDIM%KLEV), PSSFL=FLUX%PFPLSN(:,KDIM%KLEV),&
  & PCVL=PSURF%PCVL, PCVH=PSURF%PCVH, PCUR=PSURF%PCUR, PWLMX=SURFL%ZWLMX , PEVAPSNW=SURFL%ZEVAPSNW,&
!-TENDENCIES OUTPUT
  & PTSNE1=PSURF%PTSNE1, PTSAE1=PSURF%PTSAE1, PTIAE1=PSURF%PTIAE1)


!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SURFTSTP_S_LAYER',1,ZHOOK_HANDLE)
END SUBROUTINE SURFTSTP_S_LAYER
