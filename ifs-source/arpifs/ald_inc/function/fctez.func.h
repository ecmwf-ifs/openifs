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

! ---------------------------------------------------------------

! These functions are used to compute the relaxation coefficient:

! FEZBP and FEZBM are the profiles of the relaxation coefficient

! FEZE is the exponent e of z**e=x**e+y**e

REAL(KIND=JPRB) :: FEZBP, FEZBM, FEZE
REAL(KIND=JPRB) :: PZ, PEPA, PALP
INTEGER(KIND=JPIM) :: KNA, KMA

FEZBP(PZ,PEPA)=(PEPA+1.0_JPRB)*PZ**PEPA-PEPA*PZ**(PEPA+1.0_JPRB)

FEZBM(PZ,PEPA)=1.0_JPRB-(PEPA+1.0_JPRB)*(1.0_JPRB-PZ)**PEPA &
              &+PEPA*(1.0_JPRB-PZ)**(PEPA+1.0_JPRB)

FEZE(PZ,PALP,KNA,KMA)=1.0_JPRB/(PALP*(PZ**KNA)*((1.0_JPRB-PZ)**KMA))
! ------------------------------------------------------------------

