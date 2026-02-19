! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SURFRAD_LAYER(YDSURF, YDMCC,YDERAD,YDEPHY,YDRIP,KDIM, PAUX, STATE, LLKEYS, &
  & AUXL,PSURF,SURFL)

!**** *SURFRAD_LAYER* - Layer routine calling SURFRAD scheme

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------

!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
! KDIM     : Derived variable for dimensions
! PAUX     : Derived variable for general auxiliary quantities
! state    : Derived variable for model state
! LLKEYS   : Derived variable for local switches

!     ==== Input/output ====
! AUXL     : Derived variable for local auxiliary quantities
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
!      Original : 2012-11-29  F. VANA (c) ECMWF

!     MODIFICATIONS.
!     --------------
!      K. Yessad (July 2014): Move some variables.
!      R. Hogan  (Feb  2019): Pass MODIS albedo coefficients as a block
!-----------------------------------------------------------------------

USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE PARKIND1 , ONLY : JPIM , JPRB, JPRD
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMPHYDER, ONLY : DIMENSION_TYPE, STATE_TYPE, AUX_TYPE, KEYS_LOCAL_TYPE, &
   &                  AUX_DIAG_LOCAL_TYPE, SURF_AND_MORE_TYPE, SURF_AND_MORE_LOCAL_TYPE
USE YOMRIP0  , ONLY : NINDAT   ,NSSSSS
USE YOMRIP   , ONLY : TRIP
USE YOMCST   , ONLY : RDAY
USE YOERAD   , ONLY : TERAD
USE YOEPHY   , ONLY : TEPHY
USE YOMMCC   , ONLY : TMCC

!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE(TSURF)                    , INTENT(INOUT) :: YDSURF
TYPE(TEPHY)                    , INTENT(INOUT) :: YDEPHY
TYPE(TERAD)                    , INTENT(INOUT) :: YDERAD
TYPE(TMCC)                     , INTENT(INOUT) :: YDMCC
TYPE(TRIP)                     , INTENT(INOUT) :: YDRIP
TYPE (DIMENSION_TYPE)          , INTENT (IN)   :: KDIM
TYPE (AUX_TYPE)                , INTENT (IN)   :: PAUX
TYPE (STATE_TYPE)              , INTENT (IN)   :: STATE
TYPE (KEYS_LOCAL_TYPE)         , INTENT (IN)   :: LLKEYS
TYPE (AUX_DIAG_LOCAL_TYPE)     , INTENT(INOUT) :: AUXL
TYPE (SURF_AND_MORE_TYPE)      , INTENT(INOUT) :: PSURF
TYPE (SURF_AND_MORE_LOCAL_TYPE), INTENT(INOUT) :: SURFL
!-----------------------------------------------------------------------
INTEGER(KIND=JPIM) :: ID0, IM0, IY0, IDINCR,ISEC,IDD,IMM,IYY,ILMON(12), JL

! Start and end of albedo ceofficient indexing
INTEGER(KIND=JPIM) :: IALSTART, IALEND

REAL(KIND=JPRB) :: ZALBICE(KDIM%KLON)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------

#include "updcal.intfb.h"
#include "surfrad.h"
#include "icestatenemo.intfb.h"
#include "fcttim.func.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SURFRAD_LAYER',0,ZHOOK_HANDLE)
ASSOCIATE(YSURF=>YDEPHY%YSURF, &
 & NSW=>YDERAD%NSW, NLWEMISS=>YDERAD%NLWEMISS, &
 & LNEMOLIMALB=>YDMCC%LNEMOLIMALB, &
 & RSTATI=>YDRIP%RSTATI, &
 & YSD_VD=>YDSURF%YSD_VD, YSD_VF=>YDSURF%YSD_VF, YSP_RR=>YDSURF%YSP_RR, &
 & YSP_SB=>YDSURF%YSP_SB, YSP_SG=>YDSURF%YSP_SG, YSP_SL=>YDSURF%YSP_SL)
!     ------------------------------------------------------------------

!*         1.     UNROLL THE DERIVED STRUCTURES AND CALL UPDCAL AND SURFRAD


DO JL=KDIM%KIDIA,KDIM%KFDIA
  AUXL%ZWND(JL)=SQRT( STATE%U(JL,KDIM%KLEV)*STATE%U(JL,KDIM%KLEV) &
    &                +STATE%V(JL,KDIM%KLEV)*STATE%V(JL,KDIM%KLEV) )
ENDDO  

IY0=NCCAA(NINDAT)
IM0=NMM(NINDAT)
ID0=NDD(NINDAT)
IDINCR=(NSSSSS+NINT(RSTATI))/NINT(RDAY)
ISEC=MOD(NSSSSS+NINT(RSTATI),NINT(RDAY))
CALL UPDCAL(ID0,IM0,IY0,IDINCR,IDD,IMM,IYY,ILMON,-1)
ZALBICE(KDIM%KIDIA:KDIM%KFDIA) = 0.0_JPRB
IF (LNEMOLIMALB) THEN
  CALL ICESTATENEMO(YDMCC,KDIM%KSTGLO,KDIM%KIDIA,KDIM%KFDIA,PALBICE=ZALBICE(KDIM%KIDIA:KDIM%KFDIA))
ENDIF

! Two longwave emissivity intervals are requested below (outside and
! inside atmospheric window).  The SURFL%PEMIR and SURFL%PEMIS arrays
! point to columns 1 and 2 of PSPECTRALEMISS, so will be populated
! when PSPECTRALEMISS is modified by this routine.

CALL SURFRAD(YDSURF=YSURF,KDD=IDD,KMM=IMM,KMON=ILMON,KSECO=ISEC,&
  & KIDIA=KDIM%KIDIA,KFDIA=KDIM%KFDIA,KLON=KDIM%KLON,KTILES=KDIM%KTILES,KSW=NSW,&
  & KLW=NLWEMISS,&
  & LDNH=LLKEYS%LLNH,&
  & PALBF=PSURF%PSD_VF(:,YSD_VF%YALBF%MP),PALBICEF=ZALBICE,PTVH=PSURF%PSD_VF(:,YSD_VF%YTVH%MP),&
  & PALCOEFF=PSURF%PSD_VF(:,YSD_VF%IALSTART:YSD_VF%IALEND),&
  & PCURC=PSURF%PSD_VF(:,YSD_VF%YCUR%MP),PCVH=PSURF%PSD_VF(:,YSD_VF%YCVH%MP),&
  & PASN=PSURF%PSP_SG(:,1,YSP_SG%YA%MP9),PMU0=PAUX%PMU0M,PTS=PSURF%PSP_RR(:,YSP_RR%YT%MP9),PWND=AUXL%ZWND,&
  & PWS1=PSURF%PSP_SB(:,1,YSP_SB%YQ%MP9),KSOTY=PSURF%ISOTY,PFRTI=SURFL%ZFRTI,&
  & PHLICE=PSURF%PSP_SL(:,YSP_SL%YLICD%MP9),PTLICE=PSURF%PSP_SL(:,YSP_SL%YLICT%MP9),&   !LAKE
  & PALBD=SURFL%ZALBD(:,1:NSW),PALBP=SURFL%ZALBP(:,1:NSW),PALB=PSURF%PSD_VD(:,YSD_VD%YALB%MP ),&
  & PSPECTRALEMISS=SURFL%ZSPECTRALEMISS(:,1:NLWEMISS),PEMIT=PSURF%PEMIS,&
  & PALBTI=SURFL%ZALBTI,PCCNL=AUXL%ZCCNL,PCCNO=AUXL%ZCCNO,LNEMOLIMALB=LNEMOLIMALB )  

!     ------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SURFRAD_LAYER',1,ZHOOK_HANDLE)
END SUBROUTINE SURFRAD_LAYER
