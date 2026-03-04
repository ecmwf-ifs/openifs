! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction
! 
! (C) Copyright 1989- Meteo-France.
! 

#ifdef RS6K
@PROCESS NOOPTIMIZE
#endif
SUBROUTINE GRIB_CODE_MESSAGE(PTSTEP,KGRIB_HANDLE,KGRIBCD,KLEV,CDREPR,CDLTYPE,&
 & LDGRAD,LDVALID,KTOP,KBOT,KSMAX)

! No heading title ??

!    Modifications:
!    --------------
!      R. El Khatib : 01-Mar-2012 LFPOS => NFPOS
!      K. Yessad (July 2014): Move some variables.
!      R. Forbes    : 01-Mar-2014 Added MXTP/MNTP
!      M. Chrust    : 03-Jan-2020 OOPS cleaning: Model error config type replaces globals
!--------------------------------------------------------------------------

USE PARKIND1          , ONLY : JPIM, JPRB
USE YOMCT0            , ONLY : NSTEPINI, NCONF, LECMWF, L_OOPS
USE YOM_GRIB_CODES    , ONLY : NGRBBTMP
USE YOMGRIB           , ONLY : JPNUMLPP, NSFLEVS,&
 &                             NSTEPLPP ,NLOCGRB  ,NSTREAM  ,&
 &                             NLEG     ,NREFERENCE,&
 &                             NWINOFF_4V, CTYPE
USE YOMANEB           , ONLY : NANERADS
USE YOMCT3            , ONLY : NSTEP
USE YOMVAR            , ONLY : MUPTRA, LTWGRA, LTWINC, LTWBGV, LTWCGL, NFGFCLEN, LBGM
USE ALGORITHM_STATE_MOD,ONLY : GET_NUPTRA, L_IN_WRITE_INCREMENT
USE YOMHOOK           , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLCZ            , ONLY : NRITZREAD, NRITZNUMB, NRITZGBH1, NRITZGBH2
USE YOMVAREPS         , ONLY : NFCHO_TRUNC_INI, NFCLENGTH_INI, LVAREPS
USE YOMMODERR         , ONLY : YGMODERRCONF
USE YOMDYNCORE        , ONLY : LPPSTEPS
USE GRIB_API_INTERFACE, ONLY : IGRIB_SET_VALUE
USE GRIB_UTILS_MOD    , ONLY : GRIB_SET_TIME

!--------------------------------------------------------------------------

IMPLICIT NONE

REAL(KIND=JPRB)   , INTENT(IN)          :: PTSTEP
INTEGER(KIND=JPIM), INTENT(IN)          :: KGRIB_HANDLE
INTEGER(KIND=JPIM), INTENT(IN)          :: KGRIBCD
INTEGER(KIND=JPIM), INTENT(IN)          :: KLEV
CHARACTER         , INTENT(IN)          :: CDREPR*(*)
CHARACTER         , INTENT(IN)          :: CDLTYPE*(*)
LOGICAL           , INTENT(OUT)         :: LDGRAD
LOGICAL           , INTENT(OUT)         :: LDVALID
INTEGER(KIND=JPIM), INTENT(OUT)         :: KTOP
INTEGER(KIND=JPIM), INTENT(OUT)         :: KBOT
INTEGER(KIND=JPIM), OPTIONAL,INTENT(IN) :: KSMAX

!--------------------------------------------------------------------------


INTEGER(KIND=JPIM), EXTERNAL :: ISRCHEQ
INTEGER(KIND=JPIM) :: ILST, ISTEPLPP, ISMAX, ILEN, ISTEP
INTEGER(KIND=JPIM) :: IFREQS(NANERADS)

LOGICAL :: LLBAFERROR
REAL(KIND=JPRB) :: ZTSTEP
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!--------------------------------------------------------------------------

#include "surf_inq.h"

#include "abor1.intfb.h"

!--------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('GRIB_CODE_MESSAGE',0,ZHOOK_HANDLE)
!--------------------------------------------------------------------------


IF (PRESENT(KSMAX)) THEN
  IF(CDREPR == 'SH') THEN
    ISMAX=KSMAX
    CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'pentagonalResolutionParameterJ',ISMAX)
    CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'pentagonalResolutionParameterK',ISMAX)
    CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'pentagonalResolutionParameterM',ISMAX)
  ENDIF
ENDIF

LLBAFERROR = CTYPE == 'ef' .OR. CTYPE == 'ea'
IF(CTYPE == 'eb') CALL ABOR1('GRIB_CODE_MESSAGE : Unknown type -  eb ')
! IF(LLBAFERROR) CALL ABOR1('GRIB_CODE_MESSAGE : type ef and ea no longer supported ')

!*        1.1    ECMWF LOCAL GRIB

IF(NLOCGRB == 1 .OR. NLOCGRB == 36) THEN
  CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'type',CTYPE)
ENDIF

IF(NLOCGRB == 9) THEN       ! Singular vectors
  CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'forecastOrSingularVectorNumber',NRITZNUMB)
  CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'numberOfSingularVectorsComputed',NRITZREAD)
  CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'NINT_LOG10_RITZ',NRITZGBH1)
  CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'NINT_RITZ_EXP',NRITZGBH2)
ENDIF

IF(LLBAFERROR) THEN
  IF(KGRIBCD == NGRBBTMP) THEN
    IF (NSTREAM == 1249 .OR. NSTREAM == 1247) THEN
      CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'localDefinitionNumber',37)
      CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'channelNumber',KLEV)
      CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'scalingFactorForFrequencies',1)
      CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'numberOfFrequencies',NANERADS)
      CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'anoffset',NWINOFF_4V)
    ELSE
      CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'localDefinitionNumber',14)
      CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'channelNumber',KLEV)
      CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'scalingFactorForFrequencies',1)
      CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'numberOfFrequencies',NANERADS)
    ENDIF
    IFREQS(:) = 0 !Dummy values
    CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'listOfScaledFrequencies',IFREQS)
  ENDIF
ENDIF

IF(CTYPE == 'me') THEN
  IF (NSTREAM == 1249 .OR. NSTREAM == 1247) THEN
!   long window
    CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'localDefinitionNumber',39)
    CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'type',35)
    CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'number',1)
    CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'numberOfComponents',YGMODERRCONF%NCOMP_MODERR)
    CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'modelErrorType',YGMODERRCONF%NTYPE_MODERR)
    CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'anoffset',NWINOFF_4V)
  ELSE
    CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'localDefinitionNumber',25)
    CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'type',35)
    CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'number',1)
    CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'numberOfComponents',YGMODERRCONF%NCOMP_MODERR)
    CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'modelErrorType',YGMODERRCONF%NTYPE_MODERR)
  ENDIF
ENDIF


IF ((NCONF==131 .AND. LTWINC) .OR. (L_OOPS .AND. L_IN_WRITE_INCREMENT())) THEN
  IF (NSTREAM == 1249 .OR. NSTREAM == 1247) THEN
    CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'localDefinitionNumber',38)
    CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'type',33)
    CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'iteration',GET_NUPTRA())
    CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'totalNumberOfIterations',MUPTRA)
    CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'anoffset',NWINOFF_4V)
  ELSE
    CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'localDefinitionNumber',20)
    CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'type',33)
    CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'iteration',GET_NUPTRA())
    CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'totalNumberOfIterations',MUPTRA)
  ENDIF
ENDIF
IF (.NOT.LECMWF) THEN
  ! remark KY+GD: this seems necessary for use at METEO-FRANCE.
  CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'generatingProcessIdentifier',1)
ENDIF
LDGRAD=LTWGRA .OR. LTWINC .OR. LTWBGV .OR. LTWCGL

!*       1.3   TIME INDICATORS

ILST = ISRCHEQ(JPNUMLPP,NSTEPLPP(:,1),1,KGRIBCD)
IF(ILST <= JPNUMLPP) THEN
  ISTEPLPP = NSTEPLPP(ILST,2)
ELSE
  ISTEPLPP = -1
ENDIF

IF(LLBAFERROR) THEN
  ISTEP=0
  IF(CTYPE == 'ef') THEN
    ISTEP=NFGFCLEN
  ELSEIF(CTYPE /= 'eb') THEN
    IF(LBGM) ISTEP=NFGFCLEN
  ENDIF
  ZTSTEP=3600
ELSE
  ISTEP=NSTEP
  ZTSTEP=PTSTEP
ENDIF

CALL GRIB_SET_TIME(KGRIB_HANDLE,LPPSTEPS,ISTEP,&
 & ZTSTEP,NSTEPINI,&
 & LVAREPS,NLEG,&
 & NFCHO_TRUNC_INI,NFCLENGTH_INI,&
 & NREFERENCE,NSTREAM,&
 & CTYPE,ISTEPLPP,KGRIBCD,LDVALID)


!*        1.4   LEVEL(S)
KTOP = -1
KBOT =- 1
IF(CDLTYPE(1:2) == 'SF') THEN
  ILEN = SIZE(NSFLEVS)/3
  ILST = ISRCHEQ(ILEN,NSFLEVS(:,1),1,KGRIBCD)
  IF(ILST <= ILEN) THEN
    KBOT = NSFLEVS(ILST,2)
    KTOP = NSFLEVS(ILST,3)
  ENDIF
ENDIF

!--------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GRIB_CODE_MESSAGE',1,ZHOOK_HANDLE)
END SUBROUTINE GRIB_CODE_MESSAGE

