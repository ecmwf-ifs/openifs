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

SUBROUTINE WRRESF_TIME (YDRIP, CDTIME, CDFN, CD_CMYPROC)

USE PARKIND1, ONLY : JPRB, JPIM
USE YOMHOOK, ONLY : LHOOK, DR_HOOK, JPHOOK

USE YOMRES             , ONLY : LSDHM
USE YOMCT3             , ONLY : NSTEP
USE YOMMP0             , ONLY : MYPROC
USE YOMIOS             , ONLY : CIOSPRF
USE YOMOPH0            , ONLY : LINC
USE YOMRIP             , ONLY : TRIP

IMPLICIT NONE

TYPE (TRIP),       INTENT (IN)  :: YDRIP
CHARACTER (LEN=*), INTENT (OUT) :: CDTIME
CHARACTER (LEN=*), INTENT (OUT) :: CDFN
CHARACTER (LEN=*), INTENT (OUT) :: CD_CMYPROC

INTEGER(KIND=JPIM) :: IFCTIM, IDAY, IHOUR, IMIN, IINC

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('WRRESF_TIME',0,ZHOOK_HANDLE)

IF (LSDHM) THEN
  IFCTIM=(NSTEP*YDRIP%TSTEP+0.5_JPRB)/60.0_JPRB
  IDAY =IFCTIM/24/60
  IHOUR=IFCTIM/60-IDAY*24
  IMIN =IFCTIM-IDAY*24*60-IHOUR*60
  WRITE(CDTIME,'(I6.6,I2.2,I2.2)') IDAY,IHOUR,IMIN
ELSE
  IF (LINC) THEN
    IINC=NINT(REAL(NSTEP,JPRB)*YDRIP%TSTEP/3600._JPRB)
  ELSE
    IINC=NSTEP
  ENDIF
  WRITE(CDTIME,'(''+'',I4.4,''   '')') IINC
ENDIF

WRITE(CD_CMYPROC,'(A1,I4.4)')'.',MYPROC
CDFN=TRIM (CIOSPRF)//TRIM (CDTIME)//CD_CMYPROC

IF (LHOOK) CALL DR_HOOK('WRRESF_TIME',1,ZHOOK_HANDLE)

END SUBROUTINE WRRESF_TIME

