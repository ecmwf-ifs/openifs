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
SUBROUTINE HPOSSFX(KPROMA,KNUMB, KSTGLO, KPREP, PFPBUF1)
USE PARKIND1, ONLY : JPRB, JPIM
INTEGER (KIND=JPIM), INTENT (IN)    :: KPROMA
INTEGER (KIND=JPIM), INTENT (IN)    :: KNUMB
INTEGER (KIND=JPIM), INTENT (IN)    :: KSTGLO
INTEGER (KIND=JPIM), INTENT (IN)    :: KPREP
REAL (KIND=JPRB),    INTENT (INOUT) :: PFPBUF1(KPROMA)
END SUBROUTINE HPOSSFX
END INTERFACE
