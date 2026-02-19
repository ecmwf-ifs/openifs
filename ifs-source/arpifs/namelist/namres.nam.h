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

!     ------------------------------------------------------------------
!*    NAMELIST NAMRES




! NFRRES = RESTART FILE INTERVAL
! NRESTS = LIST OF RESTART TIMES
!       (1) = NR (positive),
!          THEN NRESTS(2..)*NR  DEFINES time-steps for RESTART FILE CREATION
!       (1) = NR (negative),
!          THEN NRESTS(2..)*NR  DEFINES hours      for RESTART FILE CREATION
! LSDHM = .TRUE. if time stamp is written as 'ddddhhmm' ; else time stamp
!         is controlled by LINC

!     ------------------------------------------------------------------
NAMELIST/NAMRES/NFRRES,NRESTS,LSDHM

!     ------------------------------------------------------------------


