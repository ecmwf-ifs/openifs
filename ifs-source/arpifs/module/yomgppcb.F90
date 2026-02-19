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

MODULE YOMGPPCB

USE PARKIND1, ONLY : JPRB

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------

!*    Buffer for gridpoint workfile needed for PC schemes 
!      This is an obsolescent buffer, which was used for NH dynamics
!      in the past, and which is now used also in the hydrostatic
!      dynamics: its purpose is to pass some quantities from predictor
!      to corrector step if one of the PC schemes is switched on.
!      For this reason it has been renamed in January 2008.
!     This buffer should later disappear, its content should be
!      moved in the buffers GMV, GMVS and GFL (new attributes to create).

! GPPCBUF  : buffer

REAL (KIND=JPRB), POINTER :: GPPCBUF (:,:,:) => NULL ()

!      containing pointers

!     ------------------------------------------------------------------
END MODULE YOMGPPCB

