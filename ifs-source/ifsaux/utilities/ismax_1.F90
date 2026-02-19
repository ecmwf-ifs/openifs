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

! Oct-2012 P. Marguinaud 64b LFI
      FUNCTION ISMAX_164  &
     &        ( KNELEM, PCHAMP) RESULT (IISMAX_1)
      USE LFI_PRECISION
      IMPLICIT NONE
 
!    Search for position of maximum value in array PCHAMP.
!    Simplified version (good for vector computers) assuming stride=1
!     Original  F. Vana - ONPP/CHMI - 23-Mar-2010

      INTEGER (KIND=JPLIKB) IISMAX_1
      INTEGER (KIND=JPLIKB) KNELEM
      REAL (KIND=JPDBLR) PCHAMP (*)

      IF (KNELEM <= 0) THEN
        IISMAX_1=0
      ELSE
        IISMAX_1=INT (MAXLOC(ARRAY=PCHAMP(1:KNELEM) ,DIM=1), JPLIKB)
      ENDIF
      RETURN
      END FUNCTION ISMAX_164
