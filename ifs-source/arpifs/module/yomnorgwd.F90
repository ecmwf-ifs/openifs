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

MODULE YOMNORGWD

! YOMNORGWD -- Parameters/switches for the non-orographic GW parameterization (ACNORGWD)

USE PARKIND1, ONLY : JPIM, JPRB
   
IMPLICIT NONE
   
SAVE  

TYPE TNORGWD

CHARACTER(LEN=4)   :: NORGWD_SCHEME     ! 'RAND' : ORIGINAL PARAM.
                                        ! 'PREC' : LINKS TO CONVEC. SOURCES
                                        ! 'ALLS' : LINKS TO ALL SOURCES

REAL(KIND=JPRB)    :: NORGWD_PRMAX      ! maximum of rain for which our theory applies (in kg/m^2/s)
REAL(KIND=JPRB)    :: NORGWD_DZ         ! Characteristic depth of the source
REAL(KIND=JPRB)    :: NORGWD_PTROPO     ! Tropopause height for GW spectrum (Pa)
INTEGER(KIND=JPIM) :: NORGWD_NTROPO     ! Tropopause height for GW spectrum (level)
 
REAL(KIND=JPRB)    :: NORGWD_RUWMAX     ! Max EP-Flux at Launch altitude
REAL(KIND=JPRB)    :: NORGWD_SAT        ! Saturation parameter
REAL(KIND=JPRB)    :: NORGWD_RDISS      ! Dissipation coefficient
REAL(KIND=JPRB)    :: NORGWD_DELTAT     ! Time scale of the life cycle of the waves parameterized
REAL(KIND=JPRB)    :: NORGWD_KMIN       ! Min horizontal wavenumbers
REAL(KIND=JPRB)    :: NORGWD_KMAX       ! Max horizontal wavenumbers
REAL(KIND=JPRB)    :: NORGWD_CMIN       ! Min absolute ph. vel.
REAL(KIND=JPRB)    :: NORGWD_CMAX       ! Max absolute ph. vel.
REAL(KIND=JPRB)    :: NORGWD_PLAUNCH    ! Launch height of GW spectrum (Pa)
INTEGER(KIND=JPIM) :: NORGWD_NLAUNCH    ! Launch height of GW spectrum (level)
REAL(KIND=JPRB)    :: NORGWD_PNOVERDIF  ! Bottom height for no vertical diffusion (Pa)
INTEGER(KIND=JPIM) :: NORGWD_NNOVERDIF  ! Bottom height for no vertical diffusion (level)

REAL(KIND=JPRB)    :: NORGWD_DZFRON     ! Characteristic depth of the source (front and jets source)
REAL(KIND=JPRB)    :: NORGWD_GFRON      ! Parameter G_0 (~1) that controls the amplitude of the EP flux emitted by fronts and jets
REAL(KIND=JPRB)    :: NORGWD_GB         ! Parameter that controls the amplitude of the EP flux emitted by 'background' sources

END TYPE TNORGWD

END MODULE YOMNORGWD
