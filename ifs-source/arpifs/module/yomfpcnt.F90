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

MODULE YOMFPCNT

USE PARKIND1  ,ONLY : JPIM

IMPLICIT NONE

SAVE

! High level control of post-processing management

! NFPCONF : configuration of the post-processing :
!           0 : vertical interpolation only (<CFPFMT='MODEL'>)
!           1 : gridpoint post-processing, possibly with spectral filters (<NFPOS=1>)
!           2 : gridpoint/spectral post-processing (spectral outputs possible) (<NFPOS=2>)

! LFPCNT ! control varying output variables according to time step
! LFPNAMELIST ! .TRUE. if concatenated pp namelist file exists.

! CFPNCF : control filename used to monitor output files achievement
! CNAM   : namelist filename attached to this object

! NFRFPOS    : frequency of post-processing events
! NFPOSTS    : array containing postprocessing events
! NFPOSTSMIN : array containing postprocessing events in minutes for sub-hour outputs

TYPE TFPCNT

INTEGER(KIND=JPIM) :: NFPCONF

LOGICAL :: LFPCNT
LOGICAL :: LFPNAMELIST

CHARACTER(LEN=180) :: CFPNCF='monitor'
CHARACTER(LEN=180) :: CNAM=' '

INTEGER(KIND=JPIM) :: NFRFPOS = 0
INTEGER(KIND=JPIM), ALLOCATABLE :: NFPOSTS(:)
INTEGER(KIND=JPIM), ALLOCATABLE :: NFPOSTSMIN(:)

END TYPE TFPCNT

!     ------------------------------------------------------------------
END MODULE YOMFPCNT

