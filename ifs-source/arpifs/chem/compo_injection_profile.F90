! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE COMPO_INJECTION_PROFILE &
  &( KIDIA, KFDIA, KLON, KLEV, CDTYPE, PPROFILE, &
  &  PBASEHEIGHT, PTOPHEIGHT, PTHRESHOLD, KBASELEV, KTOPLEV, &
  &  PPARAM, PBLH, PAPHIF, PDELP, PSFCFRAC )

!*** * COMPO_INJECTION_PROFILE* - COMPUTE INJECTION PROFILE FOR EMISSIONS

!**   INTERFACE.
!     ----------

!* INPUTS:
!  -------
! CDTYPE      - type of function
! PBASEHEIGHT - plume base height (m) above reference
! PTOPHEIGHT  - plume top height (m) above reference
! PTHRESHOLD  - height/altitude of auxiliary input field below which to treat as surface flux
! KBASELEV    - base level number or offset
! KTOPLEV     - top level number or offset
! PPARAM      - auxiliary parameter to use
! PBLH        - boundary layer height field
! PAPHIF      - geopotential on full levels 
! PSFCFRAC    - fraction of emission to leave as surface flux

!* OUTPUTS:
!  --------
! PPROFILE     Vertical profile scaling factor (level KLEV+1 remains as surface flux)

!     AUTHOR.
!     -------
!        Zak Kipling

USE PARKIND1  ,ONLY : JPIM, JPRB
USE YOMCST   , ONLY : RG
USE YOMHOOK   ,ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN    ,ONLY : NULOUT

!-----------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)  :: KLON, KIDIA, KFDIA, KLEV
CHARACTER(LEN=*)  ,INTENT(IN)  :: CDTYPE
REAL(KIND=JPRB)   ,INTENT(OUT) :: PPROFILE(KLON,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(IN), OPTIONAL  :: PBASEHEIGHT, PTOPHEIGHT, PTHRESHOLD
INTEGER(KIND=JPIM),INTENT(IN), OPTIONAL  :: KBASELEV, KTOPLEV
REAL(KIND=JPRB)   ,INTENT(IN), OPTIONAL  :: PPARAM(KLON)
REAL(KIND=JPRB)   ,INTENT(IN), OPTIONAL  :: PBLH(KLON)
REAL(KIND=JPRB)   ,INTENT(IN), OPTIONAL  :: PAPHIF(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN), OPTIONAL  :: PDELP(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN), OPTIONAL  :: PSFCFRAC
!-----------------------------------------------------------------------

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
REAL(KIND=JPRB) :: ZBASE(KLON), ZTOP(KLON)
REAL(KIND=JPRB) :: ZSFCFRAC, ZDELP
REAL(KIND=JPRB) :: ZTEST1(KLON,KLEV+1), ZTEST2(KLON), ZTEST3, ZTEST4
INTEGER(KIND=JPIM) :: IBASE(KLON),ITOP(KLON)
INTEGER(KIND=JPIM) :: JL, JK
LOGICAL :: LLTEST 

!-----------------------------------------------------------------------
#include "abor1.intfb.h"
!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('COMPO_INJECTION_PROFILE',0,ZHOOK_HANDLE)
!-----------------------------------------------------------------------
LLTEST=.FALSE.

! The integral of the function over the second index is supposed to be 1
PPROFILE(KIDIA:KFDIA,1:KLEV) = 0.0_JPRB 
PPROFILE(KIDIA:KFDIA,KLEV+1) = 1.0_JPRB 

IBASE(KIDIA:KFDIA) = KLEV+1
ITOP(KIDIA:KFDIA) = KLEV+1

IF (PRESENT(PSFCFRAC)) THEN
  ZSFCFRAC = PSFCFRAC
ELSE
  ZSFCFRAC = 0.0_JPRB
ENDIF

SELECT CASE (CDTYPE) 
  CASE ('Surface') ! All emissions at surface
    IBASE(KIDIA:KFDIA) = KLEV+1
    ITOP(KIDIA:KFDIA) = KLEV+1

  CASE ('HeightRange') ! Specified height range in metres
    IF (.NOT. ((PRESENT(PBASEHEIGHT) .OR. PRESENT(PTOPHEIGHT)) .AND. PRESENT(PAPHIF))) THEN
      CALL ABOR1("COMPO_INJECTION_PROFILE: CDTYPE='HeightRange' needs PBASEHEIGHT and/or PTOPHEIGHT, plus PAPHIF") 
    ENDIF
    IF (PRESENT(PBASEHEIGHT)) THEN
      ! FIXME: reference should probably be surface geopotential, not that on lowest full level
      IBASE(KIDIA:KFDIA)=MINLOC(ABS( (PAPHIF(KIDIA:KFDIA,1:KLEV) - SPREAD(PAPHIF(KIDIA:KFDIA,KLEV),2,KLEV))/RG &
                                   & - PBASEHEIGHT ), DIM=2)
    ELSE
      ! If PBASEHEIGHT missing, all the way to bottom level
      IBASE(KIDIA:KFDIA) = KLEV
    ENDIF
    IF (PRESENT(PTOPHEIGHT)) THEN
      ! FIXME: reference should probably be surface geopotential, not that on lowest full level
      ITOP(KIDIA:KFDIA)=MINLOC(ABS( (PAPHIF(KIDIA:KFDIA,1:KLEV) - SPREAD(PAPHIF(KIDIA:KFDIA,KLEV),2,KLEV))/RG &
                                  & - PTOPHEIGHT ), DIM=2)
    ELSE
      ! If PTOPHEIGHT missing, inject everything at PBASEHEIGHT
      ITOP(KIDIA:KFDIA) = IBASE(KIDIA:KFDIA)
    ENDIF

  CASE ('AltitudeRange') ! Specified altitude range in metres
    IF (.NOT. ((PRESENT(PBASEHEIGHT) .OR. PRESENT(PTOPHEIGHT)) .AND. PRESENT(PAPHIF))) THEN
      CALL ABOR1("COMPO_INJECTION_PROFILE: CDTYPE='AltitudeRange' needs PBASEHEIGHT and/or PTOPHEIGHT, plus PAPHIF") 
    ENDIF
    IF (PRESENT(PBASEHEIGHT)) THEN
      IBASE(KIDIA:KFDIA)=MINLOC(ABS((PAPHIF(KIDIA:KFDIA,1:KLEV))/RG - PBASEHEIGHT), DIM=2)
    ELSE
      ! If PBASEHEIGHT missing, all the way to bottom level
      IBASE(KIDIA:KFDIA) = KLEV
    ENDIF
    IF (PRESENT(PTOPHEIGHT)) THEN
      ITOP(KIDIA:KFDIA)=MINLOC(ABS((PAPHIF(KIDIA:KFDIA,1:KLEV))/RG - PTOPHEIGHT), DIM=2)
    ELSE
      ! If PTOPHEIGHT missing, inject everything at PBASEHEIGHT
      ITOP(KIDIA:KFDIA) = IBASE(KIDIA:KFDIA)
    ENDIF

  CASE ('LevelRange') ! Specified level range
    IF (.NOT. (PRESENT(KBASELEV) .AND. PRESENT(KTOPLEV))) THEN
      CALL ABOR1("COMPO_INJECTION_PROFILE: CDTYPE='LevelRange' needs KBASELEV and/or KTOPLEV") 
    ENDIF 
    IF (PRESENT(KBASELEV)) THEN
      IBASE(KIDIA:KFDIA) = KBASELEV
    ELSE
      ! If KBASELEV missing, all the way to bottom level
      IBASE(KIDIA:KFDIA) = KLEV
    ENDIF
    IF (PRESENT(KTOPLEV)) THEN
      ITOP(KIDIA:KFDIA) = KTOPLEV
    ELSE
      ! If KTOPLEV missing, inject everything at KBASELEV
      ITOP(KIDIA:KFDIA) = IBASE(KIDIA:KFDIA)
    ENDIF

  CASE ('HeightMap') ! Map of heights in metres, with optional height/level shifts
    IF (.NOT. (PRESENT(PPARAM) .AND. PRESENT(PAPHIF))) THEN
      CALL ABOR1("COMPO_INJECTION_PROFILE: CDTYPE='HeightMap' needs PPARAM and PAPHIF") 
    ENDIF
    IF (PRESENT(PBASEHEIGHT)) THEN
      ZBASE(KIDIA:KFDIA) = MAX(0.0_JPRB, PPARAM(KIDIA:KFDIA) + PBASEHEIGHT)
    ELSE
      ZBASE(KIDIA:KFDIA) = PPARAM(KIDIA:KFDIA)
    ENDIF
    IF (PRESENT(PTOPHEIGHT)) THEN
      ZTOP(KIDIA:KFDIA) = MAX(0.0_JPRB, PPARAM(KIDIA:KFDIA) + PTOPHEIGHT)
    ELSE
      ZTOP(KIDIA:KFDIA) = PPARAM(KIDIA:KFDIA)
    ENDIF
    IBASE(KIDIA:KFDIA)=MINLOC(ABS( (PAPHIF(KIDIA:KFDIA,1:KLEV) - SPREAD(PAPHIF(KIDIA:KFDIA,KLEV),2,KLEV))/RG &
                                 & - SPREAD(ZBASE(KIDIA:KFDIA),2,KLEV) ), DIM=2)
    ITOP(KIDIA:KFDIA)=MINLOC(ABS( (PAPHIF(KIDIA:KFDIA,1:KLEV) - SPREAD(PAPHIF(KIDIA:KFDIA,KLEV),2,KLEV))/RG &
                                & - SPREAD(ZTOP(KIDIA:KFDIA),2,KLEV) ), DIM=2)

    IF (PRESENT(KBASELEV)) IBASE(KIDIA:KFDIA)=MIN(KLEV, MAX(1, IBASE(KIDIA:KFDIA) + KBASELEV))
    IF (PRESENT(KTOPLEV)) ITOP(KIDIA:KFDIA)=MIN(KLEV, MAX(1, ITOP(KIDIA:KFDIA) + KTOPLEV))

    IF (PRESENT(PTHRESHOLD)) THEN
      WHERE (PPARAM(KIDIA:KFDIA) <= PTHRESHOLD)
        ZBASE(KIDIA:KFDIA) = 0.0_JPRB
        ZTOP(KIDIA:KFDIA) = 0.0_JPRB
        IBASE(KIDIA:KFDIA) = KLEV+1
        ITOP(KIDIA:KFDIA) = KLEV+1
      ENDWHERE
    ENDIF

  CASE ('AltitudeMap') ! Map of altitudes in metres, with optional height/level shifts
    IF (.NOT. (PRESENT(PPARAM) .AND. PRESENT(PAPHIF))) THEN
      CALL ABOR1("COMPO_INJECTION_PROFILE: CDTYPE='AltitudeMap' needs PPARAM and PAPHIF") 
    ENDIF
    IF (PRESENT(PBASEHEIGHT)) THEN
      ZBASE(KIDIA:KFDIA) = MAX(0.0_JPRB, PPARAM(KIDIA:KFDIA) + PBASEHEIGHT)
    ELSE
      ZBASE(KIDIA:KFDIA) = PPARAM(KIDIA:KFDIA)
    ENDIF
    IF (PRESENT(PTOPHEIGHT)) THEN
      ZTOP(KIDIA:KFDIA) = MAX(0.0_JPRB, PPARAM(KIDIA:KFDIA) + PTOPHEIGHT)
    ELSE
      ZTOP(KIDIA:KFDIA) = PPARAM(KIDIA:KFDIA)
    ENDIF
    IBASE(KIDIA:KFDIA)=MINLOC(ABS((PAPHIF(KIDIA:KFDIA,1:KLEV))/RG - SPREAD(ZBASE(KIDIA:KFDIA),2,KLEV)), DIM=2)
    ITOP(KIDIA:KFDIA)=MINLOC(ABS((PAPHIF(KIDIA:KFDIA,1:KLEV))/RG - SPREAD(ZTOP(KIDIA:KFDIA),2,KLEV)), DIM=2)

    IF (PRESENT(KBASELEV)) IBASE(KIDIA:KFDIA)=MIN(KLEV, MAX(1, IBASE(KIDIA:KFDIA) + KBASELEV))
    IF (PRESENT(KTOPLEV)) ITOP(KIDIA:KFDIA)=MIN(KLEV, MAX(1, ITOP(KIDIA:KFDIA) + KTOPLEV))

    IF (PRESENT(PTHRESHOLD)) THEN
      WHERE (PPARAM(KIDIA:KFDIA) <= PTHRESHOLD)
        ZBASE(KIDIA:KFDIA) = 0.0_JPRB
        ZTOP(KIDIA:KFDIA) = 0.0_JPRB
        IBASE(KIDIA:KFDIA) = KLEV+1
        ITOP(KIDIA:KFDIA) = KLEV+1
      ENDWHERE
    ENDIF

  CASE ('GFAS') ! Like HeightMap but with specific GFAS overrides for minimum height, PBL height and plume depth
    IF (.NOT. (PRESENT(PPARAM) .AND. PRESENT(PBLH) .AND. PRESENT(PAPHIF))) THEN
      CALL ABOR1("COMPO_INJECTION_PROFILE: CDTYPE='GFAS' needs PPARAM, PBLH and PAPHIF") 
    ENDIF 
    WHERE(PPARAM(KIDIA:KFDIA) > 200._JPRB .AND. PBLH(KIDIA:KFDIA) > 400._JPRB)
      ITOP(KIDIA:KFDIA)=MINLOC(ABS((PAPHIF(KIDIA:KFDIA,1:KLEV) - SPREAD(PAPHIF(KIDIA:KFDIA,KLEV),2,KLEV))/RG &
                                  & - SPREAD(PPARAM(KIDIA:KFDIA),2,KLEV) ), DIM=2)
      IBASE(KIDIA:KFDIA)=ITOP(KIDIA:KFDIA)
    ELSEWHERE
      ITOP(KIDIA:KFDIA)=MINLOC(ABS( (PAPHIF(KIDIA:KFDIA,1:KLEV) - SPREAD(PAPHIF(KIDIA:KFDIA,KLEV),2,KLEV))/RG &
                                  & - 50._JPRB ), DIM=2)
      IBASE(KIDIA:KFDIA)=KLEV
    ENDWHERE

  CASE DEFAULT 
    CALL ABOR1("COMPO_INJECTION_PROFILE: CDTYPE='"//TRIM(CDTYPE)//"' not recognised") 
END SELECT 

IF (ANY(IBASE(KIDIA:KFDIA) < KLEV+1 .AND. ITOP(KIDIA:KFDIA) <= IBASE(KIDIA:KFDIA))) THEN
  IF (.NOT. PRESENT(PDELP)) THEN
    CALL ABOR1("COMPO_INJECTION_PROFILE: injection above surface requires PDELP") 
  ENDIF
  DO JL=KIDIA,KFDIA
    IF (IBASE(JL) < KLEV+1 .AND. ITOP(JL) <= IBASE(JL)) THEN
      ! calculate total detltap over injected levels
      ZDELP=0.0_JPRB
      DO JK = ITOP(JL), IBASE(JL)
        ZDELP = ZDELP + PDELP(JL,JK)
      ENDDO
      DO JK = ITOP(JL), IBASE(JL)
        PPROFILE(JL,JK) = (1.0_JPRB-ZSFCFRAC) * RG / ZDELP
      ENDDO
      ! reduce surface flux  
      PPROFILE(JL,KLEV+1) = ZSFCFRAC
    ENDIF
  ENDDO 
ENDIF

IF (LLTEST) THEN 
  ZTEST1(KIDIA:KFDIA, 1:KLEV) = PPROFILE(KIDIA:KFDIA,1:KLEV) * PDELP(KIDIA:KFDIA,1:KLEV) / RG
  ZTEST1(KIDIA:KFDIA, KLEV+1) = PPROFILE(KIDIA:KFDIA,KLEV+1)
  ZTEST2(KIDIA:KFDIA) = SUM(ZTEST1(KIDIA:KFDIA,1:KLEV+1), DIM=2)
  ZTEST3 = MINVAL(ZTEST2(KIDIA:KFDIA))
  ZTEST4 = MAXVAL(ZTEST2(KIDIA:KFDIA))
  IF (ZTEST3 < 0.99_JPRB .OR. ZTEST4 > 1.01_JPRB) THEN
    WRITE(NULOUT, *) "COMPO_INJECTION_PROFILE: TYPE ",TRIM(CDTYPE)," INTEGRALS FROM ",ZTEST3," TO ",ZTEST4
    CALL ABOR1("COMPO_INJECTION_PROFILE: PROFILE DOES NOT INTEGRATE TO 1")
  ENDIF
ENDIF

!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('COMPO_INJECTION_PROFILE',1,ZHOOK_HANDLE)
END SUBROUTINE COMPO_INJECTION_PROFILE
