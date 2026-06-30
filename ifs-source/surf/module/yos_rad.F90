MODULE YOS_RAD

USE PARKIND1  ,ONLY : JPIM     ,JPRB

! (C) Copyright 2005- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

IMPLICIT NONE

SAVE

!       ----------------------------------------------------------------
!*    ** *YOSRAD* - CONTROL OPTIONS FOR RADIATION CONFIGURATION
!       ----------------------------------------------------------------

TYPE :: TRAD
INTEGER(KIND=JPIM) :: NSW  ! Number of SW spectral intervals in which surface albedo is defined
INTEGER(KIND=JPIM) :: NTSW ! Max number of SW spectral intervals, used to dimension arrays
! Number of spectral intervals for the UV/Vis part of the spectrum,
! with the remainder for the Near-IR. Thus the UV/Vis intervals are
! indexed 1:NUVVIS and the Near-IR intervals are indexed NUVVIS+1:NSW.
INTEGER(KIND=JPIM) :: NUVVIS
INTEGER(KIND=JPIM) :: NLWEMISS ! Number of longwave emissivity spectral intervals
INTEGER(KIND=JPIM) :: NALBEDOSCHEME ! 0=ERBE, 1=4compMODIS, 2=6compMODIS, 3=2compMODIS(=4compdiffuse)
INTEGER(KIND=JPIM) :: NEMISSSCHEME  ! 0=2-interval,1=6-interval
LOGICAL :: LCCNL   ! .T. IF CCN CONCENTRATION OVER LAND IS DIAGNOSED
LOGICAL :: LCCNO   ! .T. IF CCN CONCENTRATION OVER OCEAN IS DIAGNOSED
REAL(KIND=JPRB) :: RCCNSEA ! NUMBER CONCENTRATION (CM-3) OF CCNs OVER SEA
REAL(KIND=JPRB) :: RCCNLND ! NUMBER CONCENTRATION (CM-3) OF CCNs OVER LAND
END TYPE TRAD

END MODULE YOS_RAD
