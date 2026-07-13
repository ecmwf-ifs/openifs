MODULE YOMRIP
! (C) Copyright 2005- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

USE PARKIND1  ,ONLY : JPIM     ,JPRB, JPRD

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------

!*    Real time related variables (updated in UPDTIM)

!     NINDAT : run initial date in the form AAAAMMDD
!     NSSSSS : initial time in seconds (e.g. for 12h, 43200)
!     RTIMST : ABSOLUTE TIME OF THE MODEL AT START

!     NSTADD : NUMBER OF DAYS SINCE START OF THE MODEL
!     NSTASS : NUMBER OF SECONDS since start of model modulo(86400)
!     RSTATI : NUMBER OF SECONDS SINCE START OF THE MODEL
!     RTIMTR : ABSOLUTE TIME OF THE MODEL

!     RHGMT  : GMT TIME OF THE MODEL  (BETWEEN 0 AND 86400)
!     REQTIM : EQUATION OF TIME

!     ------------------------------------------------------------------

!*    Sun related variables (updated in UPDTIM)
!     RSOVR  : TRUE SOLAR TIME (GMT+EQUATION OF TIME)

!     RDEASO : DISTANCE EARTH-SUN
!     RDECLI : DECLINATION
!     RWSOVR : IN RADIANS, TRUE SOLAR TIME (GMT+EQUATION OF TIME)
!              THIS ANGLE IS ALSO PI - (LONGITUDE OF THE POINT
!              WHERE THE SUN IS AT ZENITH)
!     RIP0   : I0 WEIGHTED BY THE DISTANCE EARTH-SUN

!     RCODEC : COSINE OF THE DECLINATION
!     RSIDEC :   SINE OF THE DECLINATION

!     RCOVSR : COSINE OF TRUE SOLAR TIME
!     RSIVSR :   SINE OF TRUE SOLAR TIME

!     ------------------------------------------------------------------

!*    Moon related variables (updated in UPDTIM)

!     RDECLU : LUNAR DECLINATION
!     RTMOLT : IN RADIANS, TRUE LUNAR TIME (GMT+EQUATION OF TIME)
!              THIS ANGLE IS ALSO PI - (LONGITUDE OF THE POINT
!              WHERE THE MOON IS AT ZENITH)
!     RIP0LU : LUNAR I0 (DOWNWARD TOA LUNAR FLUX)

!     RCODECLU : COSINE OF THE LUNAR DECLINATION
!     RSIDECLU :   SINE OF THE LUNAR DECLINATION

!     RCOVSRLU : COSINE OF TRUE LUNAR TIME
!     RSIVSRLU :   SINE OF TRUE LUNAR TIME

!     ------------------------------------------------------------------

!*    Time step related variables

!     RDTSA  : TDT  /RA
!     RDTSA2 : RDTSA**2
!     RDTS62 : RDTSA**2/6
!     RDTS22 : RDTSA**2/2

!     RTDT   : TDT

INTEGER(KIND=JPIM) :: NINDAT
INTEGER(KIND=JPIM) :: NSSSSS
INTEGER(KIND=JPIM) :: NSTADD
INTEGER(KIND=JPIM) :: NSTASS
REAL(KIND=JPRD) :: RTIMST
REAL(KIND=JPRD) :: RSTATI
REAL(KIND=JPRD) :: RTIMTR
REAL(KIND=JPRB) :: RHGMT
REAL(KIND=JPRB) :: REQTIM
REAL(KIND=JPRB) :: RSOVR
REAL(KIND=JPRB) :: RDEASO
REAL(KIND=JPRB) :: RDECLI
REAL(KIND=JPRB) :: RWSOVR
REAL(KIND=JPRB) :: RIP0
REAL(KIND=JPRB) :: RCODEC
REAL(KIND=JPRB) :: RSIDEC
REAL(KIND=JPRB) :: RCOVSR
REAL(KIND=JPRB) :: RSIVSR
REAL(KIND=JPRB) :: RDTSA
REAL(KIND=JPRB) :: RDTSA2
REAL(KIND=JPRB) :: RDTS62
REAL(KIND=JPRB) :: RDTS22
REAL(KIND=JPRB) :: RTDT
REAL(KIND=JPRB) :: RDECLU
REAL(KIND=JPRB) :: RTMOLT
REAL(KIND=JPRB) :: RIP0LU
REAL(KIND=JPRB) :: RCODECLU
REAL(KIND=JPRB) :: RSIDECLU
REAL(KIND=JPRB) :: RCOVSRLU
REAL(KIND=JPRB) :: RSIVSRLU
!     ------------------------------------------------------------------
END MODULE YOMRIP
