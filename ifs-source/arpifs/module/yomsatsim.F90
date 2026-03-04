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

MODULE YOMSATSIM

USE PARKIND1, ONLY : JPIM, JPRB
IMPLICIT NONE

SAVE

INTEGER(KIND=JPIM), PARAMETER :: JPSATSIM=20     ! Maximum number of requested images
INTEGER(KIND=JPIM) :: NSATSIM                    ! Actual number of requested images
INTEGER(KIND=JPIM) :: MSATSIM(JPSATSIM)          ! Image request in namelist


INTEGER(KIND=JPIM) :: NINST                      ! Number of instruments
INTEGER(KIND=JPIM), ALLOCATABLE :: NTOPLEVELS(:) ! Number of RTTOV levels above IFS top
INTEGER(KIND=JPIM), ALLOCATABLE :: NCHAN(:)      ! Number of channels per instrument

INTEGER(KIND=JPIM), ALLOCATABLE :: MSATID(:)     ! WMO satellite ID
INTEGER(KIND=JPIM), ALLOCATABLE :: MSERIES(:)    ! WMO satellite series
INTEGER(KIND=JPIM), ALLOCATABLE :: MINST(:)      ! WMO instrument ID
INTEGER(KIND=JPIM), ALLOCATABLE :: MCHAN(:)      ! WMO channel ID
INTEGER(KIND=JPIM), ALLOCATABLE :: MRTCHAN(:)    ! RTTOV channel ID
REAL(KIND=JPRB),    ALLOCATABLE :: RCWN(:)       ! Central wave number (cm-1)

END MODULE YOMSATSIM
