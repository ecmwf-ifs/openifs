! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

! FUNCTIONS ARISING FROM THE GRAVITY-CAPILLARY DISPERSION RELATION

REAL(KIND=JWRB) :: XWNB
REAL(KIND=JWRB) :: FOMEG_GC, FVG_GC, FC_GC

! DISPERSION RELATION: 
FOMEG_GC(XWNB)=SQRT(G*XWNB+SURFT*XWNB**3)

! GROUP SPEED:
FVG_GC(XWNB)=0.5_JWRB/FOMEG_GC(XWNB)*(G+3.0_JWRB*SURFT*XWNB**2)

! PHASE SPEED:
FC_GC(XWNB)=FOMEG_GC(XWNB)/XWNB

