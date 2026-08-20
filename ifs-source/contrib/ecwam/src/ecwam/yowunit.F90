! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      MODULE YOWUNIT

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      IMPLICIT NONE

!*    **  *UNITS* - INPUT / OUTPUT UNITS.


!!!! some are obsolete !!! will need to be cleaned

      INTEGER(KIND=JWIM), PARAMETER :: IREADG=1
      INTEGER(KIND=JWIM), PARAMETER :: NPROPAGS=2
      INTEGER(KIND=JWIM) :: IU02
      INTEGER(KIND=JWIM) :: IU07
      INTEGER(KIND=JWIM), DIMENSION(0:NPROPAGS) :: IU08
      INTEGER(KIND=JWIM) :: IU11
      INTEGER(KIND=JWIM) :: IU12
      INTEGER(KIND=JWIM) :: IU13
      INTEGER(KIND=JWIM) :: IU14
      INTEGER(KIND=JWIM) :: IU15
      INTEGER(KIND=JWIM), ALLOCATABLE :: IU19(:)
      INTEGER(KIND=JWIM) :: IU20

      LOGICAL :: LWVWAMINIT=.TRUE.

!      *IREADG*- PROCESSOR NUMBER DOING THE READING.
!      *IU02*  - LOGICAL UNIT FOR INPUT OF BOUDARY VALUES FROM A
!                PREVIOUS COARSE GRID IF THIS A FINE GRID RUN.
!      *IU08*  - LOGICAL UNITS FOR INPUT OF  UBUF
!                (OUTPUT OF PREPROC).
!      *IU11*  - LOGICAL UNIT FOR INPUT OF SPECTRA AT ALL GRID
!                POINTS. EACH PROPAGATION TIMESTEP THE FILES
!                CONNECTED TO IU11 AND IU12 ARE INTERCHANGED.
!      *IU12*  - LOGICAL UNIT FOR OUTPUT (SEE IU11).
!      *IU13*  - LOGICAL UNIT FOR INPUT OF SPECTRA ON LAST LATUTUDE
!                OF A BLOCK. SPECTRA ARE SAVED FROM THE SECOND
!                LATITUDE OF THE NEXT BLOCK.
!                EACH PROPAGATION TIMESTEP THE FILES CONNECTED TO
!                IU14 AND IU13 ARE INTERCHANGED.
!      *IU14*  - LOGICAL UNIT FOR OUTPUT (SEE IU13).
!      *IU15*  - LOGICAL UNIT FOR OUTPUT OF LAST WINDFIELDS FOR
!                RESTART.
!      *IU19*  - LOGICAL UNIT FOR OUTPUT OF OF BOUNDARY VALUES IF
!                THIS IS A FINE GRID RUN.
!      *IU20*  - LOGICAL UNIT FOR OUTPUT OF INTEGRATED PARAMETERS
!                OF THE TOTAL SPECTRUM (FIRST GUESS) in pure binary form (OBSOLETE option !)

!      *LWVWAMINIT* - TRUE IF LWVWAMINIT HAS TO BE CALLED.
! ----------------------------------------------------------------------

      END MODULE YOWUNIT
