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

! This module define the Rosenbrock parameters given on input at the Integrator and
! AdIntegrator for the KPP chemistry solver.
MODULE CIFS_KPP_INTPARAM

  USE PARKIND1           , ONLY:   JPIM     ,JPRB

  IMPLICIT NONE

  SAVE

  
      REAL(KIND=JPRB), PARAMETER :: HMIN = 1.0E-1    &
     &                            , HMAX = 900.0     &     
     &                            , HSTART = 60.0
      REAL(KIND=JPRB), PARAMETER :: RTOLS_G= 1.0E-1
      INTEGER(KIND=JPIM), PARAMETER :: iautonom = 1  !0=non-autonom; 1=autonom
!     Choice of rosenbrock method: 0=default (Rodas3), 1=Ros2, 2=Ros3
!           3=Ros4, 4=Rodas3, 5=Rodas4
      INTEGER(KIND=JPIM), PARAMETER :: IROSMETH = 4
!     Type of adjoint algorithm: 0=default (discrete adjoint)
!           1=no adjoint, 2=discrete adjoint, 3=fully continuous 
!           4=simplified continuous
      INTEGER(KIND=JPIM), PARAMETER :: IADJALGO = 4
!     Maximum concentration before solver stops      
      REAL(KIND=JPRB), PARAMETER :: VMR_BAD_LARGE = 1E-2

END MODULE CIFS_KPP_INTPARAM
