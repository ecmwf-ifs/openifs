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

INTERFACE

SUBROUTINE FACILO_MT                                           &
&           (FA, KREP, KNUMER, CDPREF, KNIVAU, CDSUFF, PCHAMP, &
&           LDCOSP, LDUNDF, PUNDF)
USE FA_MOD, ONLY : FA_COM
USE LFI_PRECISION
IMPLICIT NONE
! Arguments
TYPE (FA_COM)          FA                                     ! INOUT
INTEGER (KIND=JPLIKM)  KREP                                   !   OUT
INTEGER (KIND=JPLIKM)  KNUMER                                 ! IN   
CHARACTER (LEN=*)      CDPREF                                 ! IN   
INTEGER (KIND=JPLIKM)  KNIVAU                                 ! IN   
CHARACTER (LEN=*)      CDSUFF                                 ! IN   
REAL (KIND=JPDBLR)     PCHAMP     (*)                         !   OUT
LOGICAL                LDCOSP                                 ! IN   
LOGICAL,               OPTIONAL :: LDUNDF                     !   OUT
REAL (KIND=JPDBLR),    OPTIONAL :: PUNDF                      !   OUT
END SUBROUTINE

END INTERFACE

