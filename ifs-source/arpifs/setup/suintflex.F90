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

SUBROUTINE SUINTFLEX(YGFL,YDARPHY,YDPHY)

! SUINTFLEX
!
!   Purpose: set up use of flexible interface
!
!   Should be done after physics and GFL setup
!

! modules
USE PARKIND1,  ONLY : JPRB
USE YOMLUN   , ONLY : NULNAM, NULOUT
USE YOM_YGFL , ONLY : TYPE_GFLD
USE YOMCST,    ONLY : RLSZER, RLVZER, RCS, RCW
USE YOMARPHY , ONLY : TARPHY
USE YOMPHY   , ONLY : TPHY
USE INTFLEX_MOD, ONLY : LINTFLEX, LENTHPREC, LRADFLEX
USE YOMHOOK,   ONLY : LHOOK, DR_HOOK, JPHOOK

IMPLICIT NONE

! interfaces
#include "posnam.intfb.h"

! namelist
#include "namintflex.nam.h"

! auxiliary variables
TYPE(TARPHY),INTENT(INOUT):: YDARPHY
TYPE(TPHY),INTENT(INOUT):: YDPHY
TYPE(TYPE_GFLD),INTENT(INOUT):: YGFL
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('SUINTFLEX',0,ZHOOK_HANDLE)
ASSOCIATE(LMPA=>YDARPHY%LMPA, &
 & YI=>YGFL%YI, YS=>YGFL%YS, YR=>YGFL%YR, YL=>YGFL%YL, YG=>YGFL%YG, &
 & LCONDWT=>YDPHY%LCONDWT, LRAY=>YDPHY%LRAY, NRAY=>YDPHY%NRAY)
! default values
LINTFLEX=.FALSE.
LENTHPREC=.FALSE.
LRADFLEX=.FALSE.

! read namelist to set nondefaults
CALL POSNAM(NULNAM,'NAMINTFLEX')
READ(NULNAM,NAMINTFLEX)

! communicate settings
WRITE (NULOUT,*) '=== Flexible interface setup ==='
WRITE (NULOUT,*) ' LINTFLEX       = ',LINTFLEX
WRITE (NULOUT,*) ' LENTHPREC      = ',LENTHPREC
WRITE (NULOUT,*) ' LRADFLEX       = ',LRADFLEX
WRITE (NULOUT,*) '================================'

! checks
IF (LRADFLEX.AND.(.NOT.LINTFLEX)) THEN
  WRITE (NULOUT,'(A)') 'LRADFLEX=.TRUE. requires LINTFLEX=.TRUE.!'
  CALL ABOR1('SUINTFLEX: ABOR1 CALLED')
ENDIF
IF ( LRAY.AND.NRAY == 2.AND.LMPA.AND.(.NOT.LRADFLEX) ) THEN
  WRITE(NULOUT,'(A)')&
   & 'ACRANEB2 radiation in AROME cannot run without LRADFLEX!'
  CALL ABOR1('SUINTFLEX: ABOR1 CALLED')
ENDIF

! change GFL settings if necessary: if water species are not prognostic (LCONDWT or LMPA), 
! the corresponding GFL variables should become inactive, but they still should have the
! correct thermodynamic properties
IF (LINTFLEX .AND. .NOT.(LCONDWT.OR.LMPA) ) THEN
  WRITE (NULOUT,*) 'SUINTFLEX: DESACTIVATING WATER SPECIES BECAUSE LCONDWT=.F. AND LMPA=.F.'
  
  ! cloud ice
  YI%LACTIVE=.FALSE.
  YI%LT1=.FALSE.
  YI%RCP=RCS
  YI%RLZER=RLSZER

  ! cloud water
  YL%LACTIVE=.FALSE.
  YL%LT1=.FALSE.
  YL%RCP=RCW
  YL%RLZER=RLVZER
  
  ! rain
  YR%LACTIVE=.FALSE.
  YR%LT1=.FALSE.
  YR%RCP=RCW
  YR%RLZER=RLVZER

  ! snow
  YS%LACTIVE=.FALSE.
  YS%LT1=.FALSE.
  YS%RCP=RCS
  YS%RLZER=RLSZER

  ! graupel
  YG%LACTIVE=.FALSE.
  YG%LT1=.FALSE.
  YG%RCP=RCS
  YG%RLZER=RLSZER
ENDIF

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUINTFLEX',1,ZHOOK_HANDLE)
END SUBROUTINE SUINTFLEX
