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

MODULE GRIB_UTILS_MOD

!     Purpose.
!     --------

!      Set GRIB headers

!     Author.
!     -------
!     Mats Hamrud(ECMWF)

!     Modifications.
!     --------------
!        Original : 2015-01-01
USE PARKIND1  ,ONLY : JPIM, JPRB, JPRD
USE GRIB_API_INTERFACE, ONLY : IGRIB_SET_VALUE, IGRIB_GET_VALUE
USE YOM_GRIB_CODES, ONLY : NGRBMX2T3, NGRBMN2T3, NGRB10FG3, NGRBMXTPR3, &
 & NGRBMNTPR3, NGRBMX2T6, NGRBMN2T6, NGRB10FG6, NGRBMXTPR6, NGRBMNTPR6, &
 & NGRBZ, NGRBLNSP, NGRBSD, &
 & NGRBCLBT, NGRBCSBT, NGRBMXCAP6, NGRBMXCAPS6, &
 & NGRBLITOTI, NGRBLITOTA1, NGRBLITOTA3, NGRBLITOTA6, &
 & NGRBLICGI, NGRBLICGA1, NGRBLICGA3, NGRBLICGA6, &
 & NGRBPTYPEMODE1, NGRBPTYPEMODE3, NGRBPTYPEMODE6, &
 & NGRBPTYPESEVR1, NGRBPTYPESEVR3, NGRBPTYPESEVR6, &
 & NGRBVIWVE, NGRBVIWVN,&
 & NGRBMUCAPE, NGRBMLCAPE50, NGRBMLCAPE100, NGRBMLCIN50, NGRBMLCIN100, NGRBMUDEPL,&
 & NGRBTROPOTP, NGRBCH4WET, &
 & NGRBSDSL, NGRBRSN, NGRBTSN, NGRBWSN
USE FITSPECTRUM_MOD    , ONLY : FITSPECTRUM, LFORCEZ, LFORCELNSP
USE ALGORITHM_STATE_MOD, ONLY :GET_NSIM4D
USE YOMVAR             , ONLY : NSIM4DL
USE YOMGRIB            , ONLY : NBITSSHLNSP, NBITSEXPR
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
IMPLICIT NONE
CONTAINS
!================================================================
FUNCTION IGRIB_EDITION(KGRIBCD,KLEV)

! Purpose - Select edition for GRIB parameter


!     Modifications.
!     --------------
!        Original : 2017-05-30

IMPLICIT NONE
INTEGER(KIND=JPIM),INTENT(IN)    :: KGRIBCD      ! Grib paramid code
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV         ! level of KGRIBCD

INTEGER(KIND=JPIM) :: IGRIB_EDITION
INTEGER(KIND=JPIM) :: IEDITION

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!---------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('GRIB_UTILS_MOD:IGRIB_EDITION',0,ZHOOK_HANDLE)

SELECT CASE(KGRIBCD)
  CASE (256000:)
    IEDITION = 2
  CASE(NGRBMXCAP6, NGRBMXCAPS6, NGRBVIWVE, NGRBVIWVN, &
     & NGRBLITOTI, NGRBLITOTA1, NGRBLITOTA3, NGRBLITOTA6, &
     & NGRBLICGI, NGRBLICGA1, NGRBLICGA3, NGRBLICGA6,&
     & NGRBMUCAPE, NGRBMLCAPE50, NGRBMLCAPE100, NGRBMLCIN50, NGRBMLCIN100, NGRBMUDEPL,&
     & NGRBTROPOTP, NGRBCH4WET, &
     & NGRBSD, NGRBWSN)
    IEDITION = 2
  CASE(NGRBRSN, NGRBTSN)
    IF (KLEV /= 0) THEN
      IEDITION = 2
    ELSE
      IEDITION = 1
    ENDIF
  
  CASE(210170, 210252:210253, 215212:215226, 217057:217206, 217222:217232, &
     & 218057:218206, 218221:218232, 219057:219232, 221057:221206, 221221:221232, &
     & 222001:222256,223001:223256)
    ! These are newer atmospheric composition parameters that have
    ! only been defined as GRIB-2, under DGOV-83, DGOV-353 and later.
    IEDITION = 2
  CASE DEFAULT
    IEDITION = 1
END SELECT

IGRIB_EDITION = IEDITION

IF (LHOOK) CALL DR_HOOK('GRIB_UTILS_MOD:IGRIB_EDITION',1,ZHOOK_HANDLE)

END FUNCTION IGRIB_EDITION
!===========================================================================
SUBROUTINE GRIB_SET_PARAMETER(KGRIB_HANDLE,KGRIBCD, &
 & KIMG,KSERIES,KSATID,KINST,KCHAN,PCWN)

! Purpose - Set GRIB parameter


!     Modifications.
!     --------------
!        Original : 2015-01-01
IMPLICIT NONE
INTEGER(KIND=JPIM),INTENT(IN)    :: KGRIB_HANDLE ! GRIB_API handle
INTEGER(KIND=JPIM),INTENT(IN)    :: KGRIBCD      ! Grib paramid code
INTEGER(KIND=JPIM),INTENT(IN)    :: KIMG         ! Satellite image number
INTEGER(KIND=JPIM),INTENT(IN)    :: KSERIES(:)   ! WMO satellite series
INTEGER(KIND=JPIM),INTENT(IN)    :: KSATID(:)    ! WMO satellite ID
INTEGER(KIND=JPIM),INTENT(IN)    :: KINST(:)     ! WMO instrument ID
INTEGER(KIND=JPIM),INTENT(IN)    :: KCHAN(:)     ! WMO channel ID
REAL(KIND=JPRB),   INTENT(IN)    :: PCWN(:)      ! Central wave number (cm-1)

INTEGER(KIND=JPIM) :: ITABLE,IPARAM,ILOCGRB,IANOFFSET,IRCWN
CHARACTER(LEN=6)   :: CLTYPE
CHARACTER(LEN=128) :: CLERR
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!---------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('GRIB_UTILS_MOD:GRIB_SET_PARAMETER',0,ZHOOK_HANDLE)

ITABLE = 128
IPARAM = KGRIBCD
IF(IPARAM > 1000) THEN
  ITABLE = KGRIBCD/1000
  IPARAM = KGRIBCD - ITABLE*1000
ENDIF
CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'type',CLTYPE)
SELECT CASE (CLTYPE(1:2))
CASE('4i','me')
  IF(ITABLE == 128) IPARAM = 200000 + IPARAM
  IF(ITABLE == 210) IPARAM = 211000 + IPARAM
CASE('sg','sf')
  IF(ITABLE == 128) IPARAM = 129000 + IPARAM
CASE DEFAULT
  IF(ITABLE /= 128) IPARAM = ITABLE*1000 + IPARAM
END SELECT

IF(IPARAM /= 128) CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'paramId',IPARAM)

! Simulated satellite images only
IF (IPARAM == NGRBCLBT .OR. IPARAM == NGRBCSBT) THEN
  CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'localDefinitionNumber',ILOCGRB)
  SELECT CASE(ILOCGRB)
  CASE(1)
    CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'localDefinitionNumber',24)
  CASE(36)
    CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'anoffset',IANOFFSET)
    CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'localDefinitionNumber',192)
    CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'subLocalDefinitionNumber1',24)
    CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'subLocalDefinitionNumber2',36)
    CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'anoffset',IANOFFSET)
  CASE DEFAULT
    WRITE(CLERR,*) 'GRIB_UTILS_MOD:GRIB_SET_PARAMETERS - SIMULATED SATELLITE IMAGES UNSUPPORTED FOR LOCAL DEFINITION ',ILOCGRB
    CALL ABOR1(CLERR)
  END SELECT

  CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'type','ssd') ! Simulated satellite data
  CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'productDefinitionTemplateNumber',32)
  CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'numberOfContributingSpectralBands',1)
  CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'satelliteSeries',KSERIES(KIMG))
  CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'satelliteNumber',KSATID(KIMG))
  CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'instrumentType',KINST(KIMG))
  CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'channelNumber',KCHAN(KIMG))
  CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'scaleFactorOfCentralWaveNumber',0)
  IRCWN = 100*NINT(PCWN(KIMG))
  CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'scaledValueOfCentralWaveNumber',IRCWN)
ENDIF


IF (LHOOK) CALL DR_HOOK('GRIB_UTILS_MOD:GRIB_SET_PARAMETER',1,ZHOOK_HANDLE)

END SUBROUTINE GRIB_SET_PARAMETER
!===========================================================================
SUBROUTINE GRIB_SET_LEVELS(KGRIB_HANDLE,KGRIBCD,KLEV,KBOT,KTOP,CDLTYPE)

! Purpose - Set GRIB level info


!     Modifications.
!     --------------
!        Original : 2015-01-01
!        G. Arduini: added multi-layer surface field meta data
IMPLICIT NONE
INTEGER(KIND=JPIM),INTENT(IN)    :: KGRIB_HANDLE ! GRIB_API handle
INTEGER(KIND=JPIM),INTENT(IN)    :: KGRIBCD      ! Grib paramid code
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV         ! Level 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTOP         ! Top level (for soil fields)
INTEGER(KIND=JPIM),INTENT(IN)    :: KBOT         ! Bottom level (for soil fields)
CHARACTER         ,INTENT(IN)    :: CDLTYPE*(*)  ! Level type (model, pressure etc.)

CHARACTER(LEN=128) :: CLERR
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!--------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GRIB_UTILS_MOD:GRIB_SET_LEVELS',0,ZHOOK_HANDLE)

IF(CDLTYPE(1:2) == 'SF' .OR. CDLTYPE(1:1) == 's') THEN
  CALL IGRIB_SET_VALUE(KGRIB_HANDLE, 'typeOfLevel','surface')
  CALL IGRIB_SET_VALUE(KGRIB_HANDLE, 'level',0)
  IF (KLEV /= 0) THEN !sfc multi-layer
     IF   (KGRIBCD == NGRBSD .OR. KGRIBCD == NGRBTSN .OR.&
          &KGRIBCD == NGRBRSN .OR. KGRIBCD == NGRBWSN) THEN
       CALL IGRIB_SET_VALUE(KGRIB_HANDLE, 'typeOfLevel','snowLayer')
       CALL IGRIB_SET_VALUE(KGRIB_HANDLE, 'typeOfFirstFixedSurface',114)
       CALL IGRIB_SET_VALUE(KGRIB_HANDLE, 'typeOfSecondFixedSurface',114)
       CALL IGRIB_SET_VALUE(KGRIB_HANDLE, 'scaleFactorOfFirstFixedSurface',0)
       CALL IGRIB_SET_VALUE(KGRIB_HANDLE, 'scaledValueOfFirstFixedSurface',KLEV-1)
       CALL IGRIB_SET_VALUE(KGRIB_HANDLE, 'scaleFactorOfSecondFixedSurface',0)
       CALL IGRIB_SET_VALUE(KGRIB_HANDLE, 'scaledValueOfSecondFixedSurface',KLEV)
       CALL IGRIB_SET_VALUE(KGRIB_HANDLE, 'level',KLEV) ! level must be defined at the end
     ENDIF
  ENDIF
  IF(KTOP >= 0 .AND. KBOT >= 0) THEN
     CALL IGRIB_SET_VALUE(KGRIB_HANDLE, 'typeOfLevel','depthBelowLandLayer')
     CALL IGRIB_SET_VALUE(KGRIB_HANDLE, 'level',0)
     CALL IGRIB_SET_VALUE(KGRIB_HANDLE, 'topLevel',KTOP)
     CALL IGRIB_SET_VALUE(KGRIB_HANDLE, 'bottomLevel',KBOT)
  ENDIF
ELSEIF(CDLTYPE(1:2) == 'ML' .OR. CDLTYPE(1:1) == 'm') THEN
  CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'typeOfLevel','hybrid')
  CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'level',KLEV)
ELSEIF(CDLTYPE(1:2) == 'PL' .OR. CDLTYPE(1:1) == 'p') THEN
  IF(KLEV >= 100) THEN
    CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'typeOfLevel','isobaricInhPa')
    CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'level',KLEV/100)
  ELSE
    CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'typeOfLevel','isobaricInPa')
    CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'level',KLEV)
  ENDIF
ELSEIF(CDLTYPE(1:2) == 'PV' .OR. CDLTYPE(1:1) == 'v') THEN
  CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'typeOfLevel','potentialVorticity')
  CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'level',KLEV)
ELSEIF(CDLTYPE(1:2) == 'TH' .OR. CDLTYPE(1:1) == 't') THEN
  CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'typeOfLevel','theta')
  CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'level',KLEV)
ELSE
  WRITE(CLERR,'(2A)') &
    & 'GRIB_UTILS_MOD:GRIB_SET_LEVELS - LEVEL TYPE NOT SUPPORTED ',TRIM(CDLTYPE)
  CALL ABOR1(CLERR)
ENDIF

IF (LHOOK) CALL DR_HOOK('GRIB_UTILS_MOD:GRIB_SET_LEVELS',1,ZHOOK_HANDLE)
END SUBROUTINE GRIB_SET_LEVELS
!=================================================================================
SUBROUTINE GRIB_SET_PACKING(KGRIB_HANDLE,KGRIBCD,CDLTYPE,LDSPEC,LDGRAD,PDATA,KSMAX)
! Purpose - Set GRIB packing


!     Modifications.
!     --------------
!        Original : 2015-01-01
IMPLICIT NONE
INTEGER(KIND=JPIM),INTENT(IN)    :: KGRIB_HANDLE ! GRIB_API handle 
INTEGER(KIND=JPIM),INTENT(IN)    :: KGRIBCD      ! Grib paramid code
CHARACTER         ,INTENT(IN)    :: CDLTYPE*(*)  ! Level type (model, pressure etc.)
LOGICAL           ,INTENT(IN)    :: LDSPEC       ! TRUE if spectral field
LOGICAL           ,INTENT(IN)    :: LDGRAD       ! TRUE if gradient field
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDATA(:)     ! Field values
INTEGER(KIND=JPIM),INTENT(IN)    :: KSMAX        ! Spectral truncation

INTEGER(KIND=JPIM), EXTERNAL :: IMULTIO_ENCODE_BITSPERVALUE

INTEGER(KIND=JPIM) :: ISPEC2, IBITSSH, ISTRUNC, IBITS, IRET
CHARACTER(LEN=128) :: CLERR
REAL(KIND=JPRB)    :: ZP, ZBETA0, ZBETA1
REAL(KIND=JPRD)    :: ZMIN, ZMAX
REAL(KIND=JPRB), PARAMETER :: ZMDI=HUGE(ZMDI)
LOGICAL :: LLMISSING
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

!--------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('GRIB_UTILS_MOD:GRIB_SET_PACKING',0,ZHOOK_HANDLE)

IF(LDSPEC) THEN ! Spectral fields
  ISPEC2 = (KSMAX+1)*(KSMAX+2)
  IF (KGRIBCD == NGRBLNSP.AND. (CDLTYPE(1:2) == 'ML' .OR. CDLTYPE(1:1) == 'm') .AND. GET_NSIM4D() /= NSIM4DL) THEN
    IBITSSH = NBITSSHLNSP
  ELSE
    IBITSSH = 16
  ENDIF  
  IF(KSMAX >= 213) THEN
    ISTRUNC = 20
  ELSE
    ISTRUNC = MIN(10,KSMAX)
  ENDIF
  CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'numberOfBitsContainingEachPackedValue',IBITSSH)
  IF(LFORCEZ .AND. KGRIBCD == NGRBZ   .AND. (CDLTYPE(1:2) == 'ML' .OR. CDLTYPE(1:1) == 'm')) THEN
    CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'laplacianOperator',0.5_JPRB)
  ELSEIF (LFORCELNSP .AND. KGRIBCD == NGRBLNSP.AND. (CDLTYPE(1:2) == 'ML' .OR. CDLTYPE(1:1) == 'm') .AND. .NOT. LDGRAD) THEN
    CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'laplacianOperator',0.5_JPRB)
  ELSEIF(ISTRUNC < KSMAX ) THEN
    ! Fit line to spectrum for complex spectral packing
    CALL FITSPECTRUM(KSMAX,ISPEC2,PDATA,ISTRUNC,ZBETA0,ZBETA1)
    ZP        = -ZBETA1
    ZP        = MAX(-9.9_JPRB,MIN(9.9_JPRB,ZP))
    CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'laplacianOperator',ZP)
  ELSE
    CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'laplacianOperator',0.0_JPRB)
  ENDIF
  CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'subSetJ',ISTRUNC)
  CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'subSetK',ISTRUNC)
  CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'subSetM',ISTRUNC)

ELSE ! Grid point fields
  IF(ANY(PDATA(:) == ZMDI)) THEN
    CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'bitmapPresent',1)
    CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'missingValue',ZMDI)
  ENDIF

  IF ((NBITSEXPR>0).AND.(KGRIBCD>=80).AND.(KGRIBCD<=120)) THEN
    IBITS = NBITSEXPR
  ELSE
    ZMIN = REAL(MINVAL(PDATA(:)),JPRD)
    ZMAX = REAL(MAXVAL(PDATA(:)),JPRD)
    IRET = IMULTIO_ENCODE_BITSPERVALUE(IBITS,KGRIBCD,TRIM(CDLTYPE),ZMIN,ZMAX)
    IF (IRET  /= 0) THEN
      WRITE(CLERR,*) 'GRIB_UTILS_MOD:GRIB_SET_PACKING: IMULTIO_ENCODE_BITSPERVALUE failed with return code ',IRET
      CALL ABOR1(CLERR)
    ENDIF
  ENDIF

  CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'numberOfBitsContainingEachPackedValue',IBITS)
ENDIF

IF (LHOOK) CALL DR_HOOK('GRIB_UTILS_MOD:GRIB_SET_PACKING',1,ZHOOK_HANDLE)
END SUBROUTINE GRIB_SET_PACKING
!=================================================================================

SUBROUTINE GRIB_SET_TIME(KGRIB_HANDLE,LDPPSTEPS,KSTEP,PTSTEP,KSTEPINI,LDVAREPS,KLEG,&
 & KFCHO_TRUNC_INI,KFCLENGTH_INI,KREFERENCE,KSTREAM,CDTYPE,KPREVPP,KGRIBCD,LDVALID)

! Purpose - Set GRIB time (step) info


!     Modifications.
!     --------------
!        Original : 2015-01-01
!        R. Forbes: May-2022 Added precip-type "mode" and "most severe" for last 1/3/6 hours

IMPLICIT NONE
INTEGER(KIND=JPIM),INTENT(IN)    :: KGRIB_HANDLE ! GRIB_API handle
LOGICAL           ,INTENT(IN)    :: LDPPSTEPS    ! TRUE if fiddling model steps into hours
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTEP        ! Model step
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSTEP       ! Time step (s)
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTEPINI     ! Shift of step (s)
LOGICAL           ,INTENT(IN)    :: LDVAREPS     ! TRUE if VAREPS system with multiple legs
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEG         ! Leg number in VAREPS
INTEGER(KIND=JPIM),INTENT(IN)    :: KFCHO_TRUNC_INI ! ?? see code below
INTEGER(KIND=JPIM),INTENT(IN)    :: KFCLENGTH_INI   ! ?? see code below
INTEGER(KIND=JPIM),INTENT(IN)    :: KREFERENCE      ! ?? see code below
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTREAM         ! MARS stream
CHARACTER         ,INTENT(IN)    :: CDTYPE*(*)      ! MARS type (an, fc etc.)
INTEGER(KIND=JPIM),INTENT(IN)    :: KPREVPP         ! Time of last  post-processing step for parameter
INTEGER(KIND=JPIM),INTENT(IN)    :: KGRIBCD         ! GRIB paramId
LOGICAL           ,INTENT(OUT)   :: LDVALID         ! Valid time for this parameter

INTEGER(KIND=JPIM) :: ISEC,ISEC0,ISECOVERLAP,ISTEPLPP,ISECLPP,ISECSTART
LOGICAL            :: LLZEROINST    ! enable stepType='accum' for ISEC=0
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!--------------------------------------------------------------------------

#include "abor1.intfb.h"

IF (LHOOK) CALL DR_HOOK('GRIB_UTILS_MOD:GRIB_SET_TIME',0,ZHOOK_HANDLE)

LLZEROINST=.TRUE.

IF( LDPPSTEPS ) THEN
  ISEC  = KSTEP*3600._JPRB  
ELSE
  ISEC  = KSTEP*PTSTEP
ENDIF

IF(TRIM(CDTYPE) == 'fc') THEN
  ISEC = ISEC+KSTEPINI*3600
ENDIF

!* CHANGE GRIB HEADERS TO WRITE FC IF VAREPS

!   Consider a VAREPS system with two legs: ModelA and modelB

!           t=0              ISEC0  ISECOVERLAP
!   ModelA:  |-----------------|----x----|

!                             t=0  ISEC
!   ModelB:                    |----x-------------------|------------|
!
!   In ModelB, to define correctly the time-step in the grib-header,
!   ISEC is shifted by ISEC0. Then, if ISEC is less than ISECOVERLAP
!   data are written in the overlap stream efov (1034-1032), otherwise
!   in stream enfo (1035-1033).

IF(LDVAREPS.AND.(KLEG >= 2)) THEN
  ISEC0=KFCHO_TRUNC_INI*3600
  ISEC =ISEC0+ISEC
  ISECOVERLAP=KFCLENGTH_INI*3600
  IF(ISEC <= ISECOVERLAP .AND. KLEG == 2 ) THEN
    IF(KREFERENCE > 0 ) THEN
      CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'stream',1032)
    ELSE
      CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'stream',1034)
    ENDIF
  ELSE
    CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'stream',KSTREAM)
  ENDIF
ELSE
  ISEC0=0
ENDIF

LDVALID = .TRUE.

CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'timeRangeIndicator',0)

!   For these parameters apply to time since the last post-processing step
IF(KPREVPP >= 0) THEN
  ISTEPLPP = KPREVPP
  ISECLPP  = ISTEPLPP*PTSTEP+ISEC0
  CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'stepType','max')
  CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'stepUnits','s')
  CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'startStep',ISECLPP)
  CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'endStep',ISEC)
ELSE
  ISECSTART = -9999
  SELECT CASE (KGRIBCD)
  CASE(NGRBMX2T3, NGRBMN2T3, NGRB10FG3, NGRBMXTPR3, NGRBMNTPR3)
    !       these parameters are min or max over the last 3 hours in GRIB edition 1
    IF (MOD(ISEC,3600) /= 0) CALL ABOR1('GRIB_SET_TIME : 3 hour min/max only on whole hours ')
    ISECSTART  = MAX(ISEC-3*3600,ISEC0)
    CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'stepType','max')
  CASE(NGRBMX2T6, NGRBMN2T6, NGRB10FG6, NGRBMXTPR6, NGRBMNTPR6)
    !       these parameters are min or max over the last 6 hours in GRIB edition 1
    IF (MOD(ISEC,3600) /= 0) CALL ABOR1('GRIB_SET_TIME : 6 hour min/max only on whole hours ')
    ISECSTART  = MAX(ISEC-6*3600,ISEC0)
    CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'stepType','max')
  CASE(NGRBLITOTA1, NGRBLICGA1)
    !       these parameters are averages over the last 1 hour in GRIB edition 2
    IF (MOD(ISEC,3600) /= 0) CALL ABOR1('GRIB_SET_TIME : 1 hour avg only on whole hours ')
    CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'stepType','avg')
    ISECSTART  = ISEC-1*3600
    IF (ISECSTART < ISEC0) LDVALID = .FALSE.
  CASE(NGRBLITOTA3, NGRBLICGA3)
    !       these parameters are averages over the last 3 hours in GRIB edition 2
    IF (MOD(ISEC,3600) /= 0) CALL ABOR1('GRIB_SET_TIME : 3 hour avg only on whole hours ')
    CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'stepType','avg')
    ISECSTART  = ISEC-3*3600
    IF (ISECSTART < ISEC0) LDVALID = .FALSE.
  CASE(NGRBLITOTA6, NGRBLICGA6)
    !       these parameters are averages over the last 6 hours in GRIB edition 2
    IF (MOD(ISEC,3600) /= 0) CALL ABOR1('GRIB_SET_TIME : 6 hour avg only on whole hours ')
    CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'stepType','avg')
    ISECSTART  = ISEC-6*3600
    IF (ISECSTART < ISEC0) LDVALID = .FALSE.
  CASE(NGRBMXCAP6, NGRBMXCAPS6)
    !       these parameters are maximums over the last 6 hours in GRIB edition 2
    IF (MOD(ISEC,3600) /= 0) CALL ABOR1('GRIB_SET_TIME : 6 hour avg only on whole hours ')
    CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'stepType','max')
    ISECSTART  = ISEC-6*3600
    IF (ISECSTART < ISEC0) LDVALID = .FALSE.
  CASE(NGRBPTYPEMODE1)
    !       these parameters are most frequent over the last 1 hours in GRIB edition 2
    IF (MOD(ISEC,3600) /= 0) CALL ABOR1('GRIB_SET_TIME : 1 hour mode (most frequent) only on whole hours ')
    CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'stepType','mode')
    ISECSTART  = ISEC-3600
    IF (ISECSTART < ISEC0) LDVALID = .FALSE.
  CASE(NGRBPTYPEMODE3)
    !       these parameters are most frequent over the last 3 hours in GRIB edition 2
    IF (MOD(ISEC,3600) /= 0) CALL ABOR1('GRIB_SET_TIME : 3 hour mode (most frequent) only on whole hours ')
    CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'stepType','mode')
    ISECSTART  = ISEC-3*3600
    IF (ISECSTART < ISEC0) LDVALID = .FALSE.
  CASE(NGRBPTYPEMODE6)
    !       these parameters are most frequent over the last 6 hours in GRIB edition 2
    IF (MOD(ISEC,3600) /= 0) CALL ABOR1('GRIB_SET_TIME : 6 hour mode (most frequent) only on whole hours ')
    CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'stepType','mode')
    ISECSTART  = ISEC-6*3600
    IF (ISECSTART < ISEC0) LDVALID = .FALSE.
  CASE(NGRBPTYPESEVR1)
    !       these parameters are most severe over the last 1 hours in GRIB edition 2
    IF (MOD(ISEC,3600) /= 0) CALL ABOR1('GRIB_SET_TIME : 1 hour most severe only on whole hours ')
    CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'stepType','severity')
    ISECSTART  = ISEC-3600
    IF (ISECSTART < ISEC0) LDVALID = .FALSE.
  CASE(NGRBPTYPESEVR3)
    !       these parameters are most severe over the last 3 hours in GRIB edition 2
    IF (MOD(ISEC,3600) /= 0) CALL ABOR1('GRIB_SET_TIME : 3 hour most severe only on whole hours ')
    CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'stepType','severity')
    ISECSTART  = ISEC-3*3600
    IF (ISECSTART < ISEC0) LDVALID = .FALSE.
  CASE(NGRBPTYPESEVR6)
    !       these parameters are most severe over the last 6 hours in GRIB edition 2
    IF (MOD(ISEC,3600) /= 0) CALL ABOR1('GRIB_SET_TIME : 6 hour most severe only on whole hours ')
    CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'stepType','severity')
    ISECSTART  = ISEC-6*3600
    IF (ISECSTART < ISEC0) LDVALID = .FALSE.
  CASE (222001:222256, 223001:223256) ! DGOV-353, has to be 'accum' also for SEC=0
    CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'stepType','accum')
    LLZEROINST=.FALSE. ! to avoid setting stepType='inst' for ISEC=0 below  
    ISECSTART  = 0 ! is that required ?    
    !IF (ISECSTART == ISEC)  LDVALID = .FALSE. ! no output at time=0 
  CASE DEFAULT
    CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'stepType','instant')
  END SELECT
  CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'stepUnits','s')
  IF (LDVALID .AND. ISECSTART >= ISEC0)  THEN
    CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'startStep',ISECSTART)
  ENDIF
  CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'endStep',ISEC)
ENDIF

IF(ISEC<= 0 .AND. LDVALID .AND. LLZEROINST ) THEN
  CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'stepType','instant')
  IF(TRIM(CDTYPE) /= 'an') THEN
    CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'timeRangeIndicator',1)
  ENDIF
ENDIF

! Reset step-units to default(hours) for downstream inquiries
IF(ISEC > 0 .AND. MOD(ISEC,3600) == 0) THEN
  CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'stepUnits','h')
ENDIF

IF (LHOOK) CALL DR_HOOK('GRIB_UTILS_MOD:GRIB_SET_TIME',1,ZHOOK_HANDLE)
END SUBROUTINE GRIB_SET_TIME
!===========================================================================

END MODULE GRIB_UTILS_MOD
