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

MODULE YOMMCUF

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------

!     MAXIMA FOR THE FILTERS FOR THE MONITORING OF THE COUPLING UPDATES

!     JPMFNR : MAXIMUM NUMBER OF RECURSIVE FILTERS
!     JPMFOR : MAXIMUM ORDER OF RECURSIVE FILTERS

!     ------------------------------------------------------------------
INTEGER(KIND=JPIM), PARAMETER :: JPMFNR=5
INTEGER(KIND=JPIM), PARAMETER :: JPMFOR=8

!     ------------------------------------------------------------------

!     Specifications of the digital filter for monitoring 
!     coupling-update frequency

!     NCUFNR  : number of filters applied         (NCUFNR<=JPMFNR)
!     NCUFOR  : the order of the recursive filter (NCUFOR<=JPMFOR)
!     LMCUF   : LOGICAL to switch on the monitoring
!     LREACUF : LOGICAL to read Coupling Update Frequency fields. 
!     RMCUFI  : the coupling-update time interval(s) under interest
!     RMCUFSP : to store the previous surface pressures
!     RMCUFFP : to store the filtered surface pressures in Spectral Space.
!              - first index for the different coefficients
!              - second index: different filters
!     RMCUFA  : A coefficients of the recursive filter
!     RMCUFB  : B coefficients of the recursive filter
!     SPFSP   : the resulting filtered field

!     ------------------------------------------------------------------

TYPE :: TMCUF

REAL(KIND=JPRB), ALLOCATABLE :: RMCUFSP(:,:)
REAL(KIND=JPRB), ALLOCATABLE :: RMCUFFP(:,:,:)
REAL(KIND=JPRB)              :: RMCUFA(0:JPMFOR,1:JPMFNR)  , &
                              & RMCUFB(1:JPMFOR,1:JPMFNR)  , &
                              & RMCUFI(1:JPMFNR)
INTEGER(KIND=JPIM)           :: NCUFOR   ,NCUFNR
LOGICAL                      :: LMCUF    ,LREACUF


END TYPE TMCUF

CONTAINS

SUBROUTINE ZERO_MCUF(YDMCUF)

USE YOMHOOK,  ONLY : LHOOK, DR_HOOK, JPHOOK

IMPLICIT NONE

TYPE(TMCUF), INTENT(INOUT) :: YDMCUF

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('YOMMCUF:ZERO_MCUF',0,ZHOOK_HANDLE)

IF (YDMCUF%LMCUF) THEN
  YDMCUF%RMCUFSP(:,:) = 0._JPRB
  YDMCUF%RMCUFFP(:,:,:) = 0._JPRB
ENDIF

IF (LHOOK) CALL DR_HOOK('YOMMCUF:ZERO_MCUF',1,ZHOOK_HANDLE)

END SUBROUTINE ZERO_MCUF

SUBROUTINE COPY_MCUF(YDMCUFOUT,YDMCUFIN)

USE YOMHOOK,  ONLY : LHOOK, DR_HOOK, JPHOOK

IMPLICIT NONE

TYPE(TMCUF), INTENT(INOUT) :: YDMCUFOUT
TYPE(TMCUF), INTENT(IN)    :: YDMCUFIN

INTEGER(KIND=JPIM) :: JBL
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

IF (LHOOK) CALL DR_HOOK('YOMMCUF:COPY_MCUF',0,ZHOOK_HANDLE)

IF(YDMCUFOUT%LMCUF .AND. YDMCUFIN%LMCUF) THEN
!$OMP PARALLEL PRIVATE(JBL)
  IF (ALLOCATED(YDMCUFOUT%RMCUFFP) .AND. ALLOCATED(YDMCUFIN%RMCUFFP)) THEN
    IF (ALL(SHAPE(YDMCUFOUT%RMCUFFP) == SHAPE(YDMCUFIN%RMCUFFP))) THEN
!$OMP DO
      DO JBL = 1, SIZE(YDMCUFIN%RMCUFFP,DIM=3)
        YDMCUFOUT%RMCUFFP(:,:,JBL) = YDMCUFIN%RMCUFFP(:,:,JBL)
      ENDDO
!$OMP END DO NOWAIT
    ELSE
      CALL ABOR1('YOMMCUF:COPY_MCUF different shapes')
    ENDIF
  ELSE
    CALL ABOR1('YOMMCUF:COPY_MCUF trying to copy into unallocated array')
  ENDIF
  IF (ALLOCATED(YDMCUFOUT%RMCUFSP) .AND. ALLOCATED(YDMCUFIN%RMCUFSP)) THEN
    IF (ALL(SHAPE(YDMCUFOUT%RMCUFSP) == SHAPE(YDMCUFIN%RMCUFSP))) THEN
!$OMP DO
      DO JBL = 0, SIZE(YDMCUFIN%RMCUFSP,DIM=2)-1
        YDMCUFOUT%RMCUFSP(:,JBL) = YDMCUFIN%RMCUFSP(:,JBL)
      ENDDO
!$OMP END DO
    ELSE
      CALL ABOR1('YOMMCUF:COPY_MCUF different shapes')
    ENDIF
  ELSE
    CALL ABOR1('YOMMCUF:COPY_MCUF trying to copy into unallocated array')
  ENDIF
!$OMP END PARALLEL
ENDIF

IF (LHOOK) CALL DR_HOOK('YOMMCUF:COPY_MCUF',1,ZHOOK_HANDLE)

END SUBROUTINE COPY_MCUF
!     ------------------------------------------------------------------

END MODULE YOMMCUF
