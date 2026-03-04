! (C) Copyright 1989- Meteo-France.

SUBROUTINE SPACONVERT(YDGEOMETRY,YDGFL,YDML_GCONF,YDDYNA,LDMODEL_TO_FILE,YDSPEC)

!**** *SPACONVERT*  - Spectral array conversion

!     Purpose.
!     --------
!        To convert model spectral fields into the fields actually written on file, 
!        and vice-versa.
!        Arpege/Ifs : Vor/Div in model is Psi/Khi in file
!        Aladin     : Vor/Div + mean wind in model is geographical wind in file
!        Both       : NH variables are scaled.
!        Both       : T-RT conversion

!**   Interface.
!     ----------
!        *CALL* *SPACONVERT(LDMODEL_TO_FILE)

!        Explicit arguments :
!        --------------------
!           LDMODEL_TO_FILE : .TRUE. to convert from model to file

!        Implicit arguments :
!        --------------------
!        None

!     Method.
!     -------
!        See documentation

!     Externals.    
!     ----------    

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        R. El Khatib  *METEO-FRANCE*
!        Original : 98-03-30

!     Modifications.
!     --------------
!   N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!   R. El Khatib   04-Sep-2008 T-RT conversion
!   O. Marsden     August 2016 Cleanups due to removal of SPA3
!   K. Yessad (Feb 2018): remove deep-layer formulations.
!     ------------------------------------------------------------------

USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMDYNA  , ONLY : TDYNA
USE YOMGFL   , ONLY : TGFL
USE PARKIND1 , ONLY : JPRB, JPIM
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCT0   , ONLY : LELAM
USE SPECTRAL_FIELDS_MOD, ONLY : SPECTRAL_FIELD, ASSIGNMENT(=)

IMPLICIT NONE

TYPE(GEOMETRY)      ,INTENT(IN)    :: YDGEOMETRY
TYPE(TGFL)          ,INTENT(INOUT) :: YDGFL
TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(IN):: YDML_GCONF
TYPE(TDYNA)         ,INTENT(IN)    :: YDDYNA
LOGICAL             ,INTENT(IN)    :: LDMODEL_TO_FILE 
TYPE(SPECTRAL_FIELD),INTENT(INOUT) :: YDSPEC

!     ------------------------------------------------------------------

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "espconvert.intfb.h"
#include "spconvert.intfb.h"
#include "spnh_conv_nhvar.intfb.h"
#include "especrt.intfb.h"
#include "specrt.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SPACONVERT',0,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

IF (LDMODEL_TO_FILE) THEN
! From model to file :
  IF (LELAM) THEN
!   Convert RT to T if relevent
    IF (YDDYNA%LSPRT) THEN
      CALL ESPECRT(YDGEOMETRY,YDGFL,YDML_GCONF%YGFL,1,YDSPEC)
    ENDIF
!   Scale NH variable AFTER RT->T conversion because it uses T not RT
    IF (YDDYNA%LNHDYN) THEN
      CALL SPNH_CONV_NHVAR(YDGEOMETRY,YDML_GCONF,YDDYNA,.TRUE.,YDSPEC)
    ENDIF
!   Convert (vor,div,mean-U,mean-V) to geographical (U,V) at the end
    CALL ESPCONVERT(YDGEOMETRY,.TRUE.,YDSPEC) 
  ELSE
!   Convert RT to T if relevent
    IF (YDDYNA%LSPRT) THEN
      CALL SPECRT(YDGEOMETRY,YDGFL,YDML_GCONF%YGFL,1,YDSPEC)
    ENDIF
!   Scale NH variable AFTER RT->T conversion because it uses T not RT
    IF (YDDYNA%LNHDYN) THEN
      CALL SPNH_CONV_NHVAR(YDGEOMETRY,YDML_GCONF,YDDYNA,.TRUE.,YDSPEC)
    ENDIF
!   Convert (vor,div) to (psi,khi) at the end
    CALL SPCONVERT(YDGEOMETRY,.TRUE., YDSPEC)
  ENDIF

ELSE
! From file to model :
  IF (LELAM) THEN
!   Convert geographical (U,V) to (vor,div,mean-U,mean-V) first
    CALL ESPCONVERT(YDGEOMETRY,.FALSE.,YDSPEC)
!   Unscale NH variable BEFORE T->RT conversion because it uses T not RT
    IF (YDDYNA%LNHDYN) THEN
      CALL SPNH_CONV_NHVAR(YDGEOMETRY,YDML_GCONF,YDDYNA,.FALSE.,YDSPEC)
    ENDIF
!   Convert T to RT if relevent
    IF (YDDYNA%LSPRT) THEN
      CALL ESPECRT(YDGEOMETRY,YDGFL,YDML_GCONF%YGFL,0,YDSPEC)
    ENDIF
  ELSE
!   Convert (psi,khi) to (vor,div) first
    CALL SPCONVERT(YDGEOMETRY,.FALSE., YDSPEC)
!   Unscale NH variable BEFORE T->RT conversion because it uses T not RT
    IF (YDDYNA%LNHDYN) THEN
      CALL SPNH_CONV_NHVAR(YDGEOMETRY,YDML_GCONF,YDDYNA,.FALSE.,YDSPEC)
    ENDIF
!   Convert T to RT if relevent
    IF (YDDYNA%LSPRT) THEN
      CALL SPECRT(YDGEOMETRY,YDGFL,YDML_GCONF%YGFL,0,YDSPEC)
    ENDIF
  ENDIF
ENDIF

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SPACONVERT',1,ZHOOK_HANDLE)
END SUBROUTINE SPACONVERT
