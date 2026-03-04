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

MODULE YOMJBCHVAR

USE PARKIND1  ,ONLY : JPRB

IMPLICIT NONE

SAVE

PRIVATE
PUBLIC :: TYPE_JBCHVAR

!     ------------------------------------------------------------------
!*    Coefficients of the change of variable for the 
!     background constraint Jb. Calculated in sujbchvar.F90

!     Humidity change of variable:
!     RQS(0:L)            : TROPOPAUSE TRANSITION OF SIGMA REDUCTION FOR Q
!     RQT(0:K)            : TROPOPAUSE TRANSITION OF BALANCE COEFFICIENT
!                           BETWEEN Q AND T
!     RQTBAL (NFLEVG,0:N) : BALANCE COEFFICIENT BETWEEN Q AND T
!     RQLMIN(NFLEVG)      : MINIMUM Q/Q_s, SUBSATURATED   (RH^b<1.x)
!     RQLMAX(NFLEVG)      : MAXIMUM Q/Q_s, SUBSATURATED   (RH^b<1.x)
!     RQLSTD (NFLEVG,0:N) : NORMALIZATION STANDARD DEVIATION COEFFICIENTS (RH^b<1.x)
!     RQHMIN(NFLEVG)      : MINIMUM Q/Q_s, SUPERSATURATED (RH^b>1.x)
!     RQHMAX(NFLEVG)      : MAXIMUM Q/Q_s, SUPERSATURATED (RH^b>1.x)
!     RQHSTD (NFLEVG,0:N) : NORMALIZATION STANDARD DEVIATION COEFFICIENTS (RH^b<1.x)
!
! Modifications
! -------------
!    M. Fisher   7-March-2012 Use DEALLOCATE_IF_ASSOCIATED
!-------------------------------------------------------------------------------

TYPE TYPE_JBCHVAR
  REAL(KIND=JPRB),ALLOCATABLE :: RQS(:)
  REAL(KIND=JPRB),ALLOCATABLE :: RQT(:)
  REAL(KIND=JPRB),ALLOCATABLE :: RQTBAL(:,:)
  REAL(KIND=JPRB),ALLOCATABLE :: RQLMIN(:)
  REAL(KIND=JPRB),ALLOCATABLE :: RQLMAX(:)
  REAL(KIND=JPRB),ALLOCATABLE :: RQLSTD(:,:)
  REAL(KIND=JPRB),ALLOCATABLE :: RQHMIN(:)
  REAL(KIND=JPRB),ALLOCATABLE :: RQHMAX(:)
  REAL(KIND=JPRB),ALLOCATABLE :: RQHSTD(:,:)
END TYPE TYPE_JBCHVAR

!     ------------------------------------------------------------------

END MODULE YOMJBCHVAR



