! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE WVCOUPLE_UPDATE_GRIB_HANDLES(LDWINIT, KSTPW, PTSTEP, KSTEP, LDPPSTEPS, KGRIB_HANDLE_FOR_WAM)
#ifdef WITH_WAVE

!----------------------------------------------------------------------

!**** *WVCOUPLE_UPDATE_GRIB_HANDLES* 

!     J. HAWKES    ECMWF MARCH 2018

!*    PURPOSE.
!     --------
!     	Initialize and update the GRIB handle (KGRIB_HANDLE_FOR_WAM) for WAM when coupled to IFS.
!		If LDWINIT == false, a new grib handle will be created.
!		If LDWINIT == true, the existing grib handle will be updated.


!**   INTERFACE.
!     ----------

!     SUBROUTINE WVCOUPLE_UPDATE_GRIB_HANDLES
!                INPUT: YDGEOMETRY
!                OUTPUT:

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------
!     WVCOUPLE_INIT_EARLY
!     WVXF2GB

!     REFERENCE.
!     ----------
!     EXTRACTED FROM WVXF2FB
!

!-------------------------------------------------------------------

USE YOMHOOK,            ONLY : LHOOK,   DR_HOOK, JPHOOK
USE PARKIND1,           ONLY : JPIM, JPRB
USE IOSTREAM_MIX,       ONLY : YGBH
USE GRIB_API_INTERFACE, ONLY : IGRIB_CLONE, IGRIB_SET_VALUE
USE GRIB_UTILS_MOD,     ONLY : GRIB_SET_PARAMETER

IMPLICIT NONE

LOGICAL,            INTENT(IN)    :: LDWINIT              ! Is wave model already initialized?
INTEGER(KIND=JPIM), INTENT(IN)    :: KSTPW                ! Frequency of call to the wave model (from YREWCOU)
REAL(KIND=JPRB),    INTENT(IN)    :: PTSTEP               ! IFS timestep (from YDRIP)
INTEGER(KIND=JPIM), INTENT(IN)    :: KSTEP                ! Current IFS step (from YOMCT3)
LOGICAL,            INTENT(IN)    :: LDPPSTEPS            ! (from YOMDYNCORE)
INTEGER(KIND=JPIM), INTENT(INOUT) :: KGRIB_HANDLE_FOR_WAM ! The grib handle to be created/updated

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

INTEGER(KIND=JPIM) :: ILEV, IGRIBCD, IDUMMY(1)
CHARACTER (LEN=2) :: CLREPR
CHARACTER (LEN=3) :: CLTYPE
INTEGER(KIND=JPIM) :: ITOP_DUM,IBOT_DUM
LOGICAL :: LLGRAD_DUM, LLDUMMY
REAL(KIND=JPRB) :: ZDUMMY(1)
INTEGER(KIND=JPIM) :: ISTEP_OUT

#include "grib_code_message.intfb.h"

!-------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('WVCOUPLE_UPDATE_GRIB_HANDLES',0,ZHOOK_HANDLE)

CALL GSTATS(1714,0)
IF (.NOT.LDWINIT) THEN
! INITIALISE GRIB HANDLE FOR WAM
IF(YGBH%NGRIB_HANDLE_GG < 0 ) THEN
CALL ABOR1("WVXF2GB: PROBLEM NGRIB_HANDLE_GG < 0 ")
ELSE
KGRIB_HANDLE_FOR_WAM=-99
CALL IGRIB_CLONE(YGBH%NGRIB_HANDLE_GG,KGRIB_HANDLE_FOR_WAM)
ENDIF
ENDIF

! UPDATE GRIB HANDLE FOR WAM
! ONLY THE FIRST TIME AND EVERY HOUR (3 or 6) for long runs AFTER THAT 
! The wave model does not output more often than every hour
! (it will be longer it TSEP is not a multiple of one hour)
! We could do it all the time but then there is a potential
! problem in *GRIB_CODE_MESSAGE* whereby grib_api is unable
! to find an appropriate time unit when coding the step (it is a GRIB-1 limitation)  
! There is another GRIB-1 limitation that limit the forecast to be less than 2**16=65536
! whatever the time unit is. This limits the hourly output to (2**16-1) hours
! Instead, one can hope to have 3 hourly output up to (2**16-1)*3 hours
! Followed by 6 hourly output up to 2**16-1)*6 hours and so on.
IF(MOD(3600,NINT(PTSTEP)) == 0 .AND. NINT(KSTEP*(PTSTEP/3600._JPRB)) <= 65535) THEN 
  ISTEP_OUT=MAX(NINT(3600._JPRB/PTSTEP),1)
ELSEIF(MOD(10800,NINT(PTSTEP)) == 0 .AND. NINT(KSTEP*(PTSTEP/10800._JPRB)) <= 65535 ) THEN 
  ISTEP_OUT=MAX(NINT(10800._JPRB/PTSTEP),1)
ELSEIF(MOD(21600,NINT(PTSTEP)) == 0 .AND. NINT(KSTEP*(PTSTEP/21600._JPRB)) <= 65535 ) THEN 
  ISTEP_OUT=MAX(NINT(21600._JPRB/PTSTEP),1)
ELSEIF(MOD(43200,NINT(PTSTEP)) == 0 .AND. NINT(KSTEP*(PTSTEP/43200._JPRB)) <= 65535 ) THEN 
  ISTEP_OUT=MAX(NINT(43200._JPRB/PTSTEP),1)
ELSEIF(MOD(86400,NINT(PTSTEP)) == 0 .AND. NINT(KSTEP*(PTSTEP/86400._JPRB)) <= 65535 ) THEN 
  ISTEP_OUT=MAX(NINT(86400._JPRB/PTSTEP),1)
ELSE
  ISTEP_OUT=1
ENDIF

IF ((KSTEP == KSTPW) .OR. MOD(KSTEP,ISTEP_OUT) == 0 .OR. LDPPSTEPS ) THEN
  IGRIBCD=165
  ILEV=0
  CLREPR='GG'
  CLTYPE='SFC'
  CALL GRIB_CODE_MESSAGE(PTSTEP,KGRIB_HANDLE_FOR_WAM,IGRIBCD,ILEV,CLREPR,CLTYPE,&
&   LLGRAD_DUM,LLDUMMY,ITOP_DUM,IBOT_DUM)
  CALL GRIB_SET_PARAMETER(KGRIB_HANDLE_FOR_WAM,IGRIBCD,ILEV,&
&   IDUMMY,IDUMMY,IDUMMY,IDUMMY,ZDUMMY)
  ZDUMMY=0.0_JPRB
  CALL IGRIB_SET_VALUE(KGRIB_HANDLE_FOR_WAM,'values',ZDUMMY)
ENDIF
CALL GSTATS(1714,1)

IF (LHOOK) CALL DR_HOOK('WVCOUPLE_UPDATE_GRIB_HANDLES',1,ZHOOK_HANDLE)

#endif
END SUBROUTINE WVCOUPLE_UPDATE_GRIB_HANDLES 
