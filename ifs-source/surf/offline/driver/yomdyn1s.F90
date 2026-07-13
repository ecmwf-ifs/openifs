MODULE YOMDYN1S
! (C) Copyright 2005- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE
SAVE

!     ------------------------------------------------------------------

!*    Control variables for the temporal evolution

!=========== TIME STEPPING ================

! NSTEP   : current timestep of model
! TSTEP   : length of the timestep in seconds
! TDT     : 2*TSTEP except at the first time step where it is TSTEP
! NACCTYPE: accumulation type for the forcing fluxes w.r.t. timestamp
!           0=centred ; 1=forward ; 2=backward
!           |--x--|     |x----|     |----x|
!           0 is default for 1d site (forcing Tstep is short)
!           1 is for GSWP2 runs ; 2 for IFS based runs (ERA-40, AMMA, ...)
! LSWINT   L : LOGICAL:  INTERPOLATE SOLAR DOWNWARD RADIATION FOLLOWING SOLAR ZENITH ANGLE, DEFAULT = FALSE 
! LPREINT  L : LOGICAL:  DISTRIBUTE PRECIPITATION, DEFAULT = FALSE
! LFLXINT  L : LOGICAL:  Linear Interpolation for solar rad and thermal rad (works for for NACCTYPE 0 and 2
! TCOUPFREQ  : Frequency of coupling with Cama-flood (in seconds)

INTEGER(KIND=JPIM) :: NSTEP
REAL(KIND=JPRB) :: TSTEP
REAL(KIND=JPRB) :: TDT
REAL(KIND=JPRB) :: TCOUPFREQ
INTEGER(KIND=JPIM) :: NACCTYPE
LOGICAL LPREINT 
LOGICAL LSWINT
LOGICAL LFLXINT


!     ------------------------------------------------------------------
END MODULE YOMDYN1S
