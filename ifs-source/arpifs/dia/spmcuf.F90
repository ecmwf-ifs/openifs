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

SUBROUTINE SPMCUF(YDGEOMETRY,YDSP,YDMCUF)

! Purpose:
! --------  
!   *SPMCUF*  - Filters ln P_s with a recursive high-pass filter
!               providing a diagnostic instantaneous field for 
!               monitoring the update frequency of the coupling files 
!               for ALADIN

! Interface:
! ----------
!   *CALL* *SPMCUF*

! Externals:
! ----------
!   None.

! Method:
! -------
!   See documentation.

! Reference:
! ----------
!   Mon. Wea. Rev., 132, 2130-2141

! Author:
! -------
!   01-Sep-2004 Piet Termonia                    *RMI Belgium*

! Modifications:
! --------------
!  O. Marsden  Sept 2016 : Removed use of SPA3
! End Modifications
!------------------------------------------------------------------

USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE PARKIND1           , ONLY : JPIM, JPRB
USE YOMCT3             , ONLY : NSTEP
USE YOMHOOK            , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMMP0             , ONLY : MYSETV
USE YOMMCUF            , ONLY : TMCUF
USE SPECTRAL_FIELDS_MOD, ONLY : SPECTRAL_FIELD,ASSIGNMENT(=)

IMPLICIT NONE

TYPE(GEOMETRY)      , INTENT(IN) :: YDGEOMETRY
TYPE(SPECTRAL_FIELD), INTENT(IN) :: YDSP
TYPE(TMCUF)         , INTENT(INOUT) :: YDMCUF
INTEGER(KIND=JPIM) :: JOR      ,JNR      ,JSP
REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

!------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SPMCUF',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NSPEC2=>YDDIM%NSPEC2, NBSETSP=>YDMP%NBSETSP)

IF (MYSETV /= NBSETSP) THEN
  IF (LHOOK) CALL DR_HOOK('SPMCUF',1,ZHOOK_HANDLE)
  RETURN
ENDIF

! 1. Filter
! ---------

DO JSP=1,NSPEC2
  YDMCUF%RMCUFSP(JSP,0)=YDSP%SP(JSP)
ENDDO
IF (NSTEP > YDMCUF%NCUFOR) THEN
  DO JNR=1,YDMCUF%NCUFNR
    DO JSP=1,NSPEC2
      YDMCUF%RMCUFFP(JSP,0,JNR)=YDMCUF%RMCUFA(0,JNR)*YDMCUF%RMCUFSP(JSP,0)
    ENDDO
    DO JOR=1,YDMCUF%NCUFOR
      DO JSP=1,NSPEC2
        YDMCUF%RMCUFFP(JSP,0,JNR)=YDMCUF%RMCUFFP(JSP,0,JNR)         &
             & + YDMCUF%RMCUFA(JOR,JNR)*YDMCUF%RMCUFSP(JSP,JOR)     &
             & - YDMCUF%RMCUFB(JOR,JNR)*YDMCUF%RMCUFFP(JSP,JOR,JNR)
      ENDDO
    ENDDO
  ENDDO
ENDIF

! 2. Time shifting
! ----------------

DO JOR=YDMCUF%NCUFOR,1,-1
  DO JSP=1,NSPEC2
    YDMCUF%RMCUFSP(JSP,JOR)=YDMCUF%RMCUFSP(JSP,JOR-1)
  ENDDO
ENDDO
DO JNR=1,YDMCUF%NCUFNR
  DO JOR=YDMCUF%NCUFOR,1,-1
    DO JSP=1,NSPEC2
      YDMCUF%RMCUFFP(JSP,JOR,JNR)=YDMCUF%RMCUFFP(JSP,JOR-1,JNR)
    ENDDO
  ENDDO
ENDDO

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SPMCUF',1,ZHOOK_HANDLE)

END SUBROUTINE SPMCUF
