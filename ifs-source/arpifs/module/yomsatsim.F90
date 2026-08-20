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

INTEGER(KIND=JPIM), PARAMETER :: JPSATSIM=40        ! Maximum number of requested images
INTEGER(KIND=JPIM) :: NSATSIM=0                     ! Actual number of requested images
INTEGER(KIND=JPIM) :: MSATSIM(JPSATSIM)             ! Channels to be written listed in namelist
INTEGER(KIND=JPIM) :: MSATFIELDS(JPSATSIM)          ! Field to be written for each channel

INTEGER(KIND=JPIM) :: NINST                         ! Number of instruments
INTEGER(KIND=JPIM), ALLOCATABLE :: NTOPLEVELS(:)    ! Number of RTTOV levels above IFS top
INTEGER(KIND=JPIM), ALLOCATABLE :: NCHAN(:)         ! Number of channels per instrument

INTEGER(KIND=JPIM), ALLOCATABLE :: MSATID(:)        ! WMO satellite ID
INTEGER(KIND=JPIM), ALLOCATABLE :: MSERIES(:)       ! WMO satellite series
INTEGER(KIND=JPIM), ALLOCATABLE :: MINST(:)         ! WMO instrument ID
INTEGER(KIND=JPIM), ALLOCATABLE :: MCHAN(:)         ! WMO channel ID
INTEGER(KIND=JPIM), ALLOCATABLE :: MRTCHAN(:)       ! RTTOV channel ID
REAL(KIND=JPRB),    ALLOCATABLE :: RCWN(:)          ! Central wave number (cm-1)

INTEGER(KIND=JPIM), ALLOCATABLE :: MCLBT_FIELDS(:)  ! Channels to write cloudy brightness temperature
INTEGER(KIND=JPIM), ALLOCATABLE :: MCSBT_FIELDS(:)  ! Channels to write clear-sky brightness temperature
INTEGER(KIND=JPIM), ALLOCATABLE :: MCDRFL_FIELDS(:) ! Channels to write cloudy reflectance
INTEGER(KIND=JPIM), ALLOCATABLE :: MCRRFL_FIELDS(:) ! Channels to write clear-sky reflectance

! Control whether traj structures are reallocated for every call or once at initialisation
LOGICAL                         :: LREUSE_TRAJ

END MODULE YOMSATSIM
