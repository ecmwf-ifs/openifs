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

!     EFFICENT WAVE NUMBER FOR EACH WAVE NUMBER COUPLE (N,M)

REAL(KIND=JPRB) :: FKSTAR
INTEGER(KIND=JPIM) :: KN, KM, KNSMAX, KNMSMAX

FKSTAR(KN,KM,KNSMAX,KNMSMAX)=INT(KNSMAX*SQRT((REAL(KN,JPRB)&
                               &/REAL(KNSMAX,JPRB))**2 &
                               &+(REAL(KM,JPRB)&
                               &/REAL(KNMSMAX,JPRB))**2)&
                               &+ 0.49_JPRB)
!     ------------------------------------------------------------------



