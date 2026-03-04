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

SUBROUTINE SUINIMODERR(YDGEOMETRY,YDRIP,YDMODERRCONF)

!   Purpose.
!   --------
!     Initialise model error variables from file.

!   Author.
!   -------
!     Y. Tremolet

!   Modifications.
!   --------------
!     Original   09-Jul-2004
!     M. Chrust  (Jan 2020): OOPS cleaning
! ------------------------------------------------------------------

USE YOMRIP             , ONLY : TRIP
USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE PARKIND1           , ONLY : JPRB, JPIM
USE YOMHOOK            , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN             , ONLY : NULOUT
USE YOMCT0             , ONLY : LIFSMIN, LIFSTRAJ
USE ALGORITHM_STATE_MOD, ONLY : GET_NUPTRA
USE YOMMODERRCONF      , ONLY : TMODERR_CONF
USE YOMMODERR          , ONLY : SPMODERR, GPMODERR, LFGMODERR, SPFGMODERR, GPFGMODERR, &
  &                             SPCTLMODERR, SPGPMODERR
USE SPECTRAL_FIELDS_MOD, ONLY : ASSIGNMENT(=), SPECTRAL_FIELD
USE GRIDPOINT_FIELDS_MIX, ONLY : ASSIGNMENT(=)

IMPLICIT NONE
TYPE(GEOMETRY),     INTENT(INOUT) :: YDGEOMETRY
TYPE(TRIP),         INTENT(INOUT) :: YDRIP
TYPE(TMODERR_CONF), INTENT(IN)    :: YDMODERRCONF
INTEGER(KIND=JPIM) :: JJ
CHARACTER(LEN=14) :: CLFILE
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "rdmoderr.intfb.h"

! ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUINIMODERR',0,ZHOOK_HANDLE)

IF (LIFSTRAJ) THEN
  DO JJ=1,YDMODERRCONF%NDIM_MODERR
    CLFILE(1:14)='spmoderrfgXXYY'
    WRITE(CLFILE(11:12),'(I2.2)')GET_NUPTRA()
    WRITE(CLFILE(13:14),'(I2.2)')JJ
    WRITE(NULOUT,*)'SUINIMODERR: Reading model error from ',CLFILE
    CALL RDMODERR(YDMODERRCONF,YDGEOMETRY,YDRIP,CLFILE, &
      &           SPCTLMODERR,SPGPMODERR,SPMODERR(JJ),GPMODERR(JJ))
  ENDDO

ELSEIF (LIFSMIN) THEN
  IF (LFGMODERR) THEN
    WRITE(NULOUT,*)'SUINIMODERR: setting MODERR from FGMODERR'
    DO JJ=1,YDMODERRCONF%NDIM_MODERR
      SPMODERR(JJ)=SPFGMODERR(JJ)
      GPMODERR(JJ)=GPFGMODERR(JJ)
    ENDDO
  ELSE
    WRITE(NULOUT,*)'SUINIMODERR: setting MODERR to ZERO'
    DO JJ=1,YDMODERRCONF%NDIM_MODERR
      SPMODERR(JJ)=0.0_JPRB
      GPMODERR(JJ)=0.0_JPRB
    ENDDO
  ENDIF

ELSE  ! Forecast
  DO JJ=1,YDMODERRCONF%NDIM_MODERR
    CLFILE(1:14)='spmoderror_YYY'
    WRITE(CLFILE(12:14),'(I3.3)')JJ
    WRITE(NULOUT,*)'SUINIMODERR: Reading model error from ',CLFILE
    CALL RDMODERR(YDMODERRCONF,YDGEOMETRY,YDRIP,CLFILE, &
      &           SPCTLMODERR,SPGPMODERR,SPMODERR(JJ),GPMODERR(JJ))
  ENDDO
ENDIF

IF (LHOOK) CALL DR_HOOK('SUINIMODERR',1,ZHOOK_HANDLE)
END SUBROUTINE SUINIMODERR
