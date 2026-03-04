! (C) Copyright 1988- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE YOETHF

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------
!*     *YOETHF* DERIVED CONSTANTS SPECIFIC TO ECMWF THERMODYNAMICS
!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: R2ES
REAL(KIND=JPRB) :: R3LES
REAL(KIND=JPRB) :: R3IES
REAL(KIND=JPRB) :: R4LES
REAL(KIND=JPRB) :: R4IES
REAL(KIND=JPRB) :: R5LES
REAL(KIND=JPRB) :: R5IES
REAL(KIND=JPRB) :: RVTMP2
REAL(KIND=JPRB) :: RHOH2O
REAL(KIND=JPRB) :: R5ALVCP
REAL(KIND=JPRB) :: R5ALSCP
REAL(KIND=JPRB) :: RALVDCP
REAL(KIND=JPRB) :: RALSDCP
REAL(KIND=JPRB) :: RALFDCP
REAL(KIND=JPRB) :: RTWAT
REAL(KIND=JPRB) :: RTBER
REAL(KIND=JPRB) :: RTBERCU
REAL(KIND=JPRB) :: RTICE
REAL(KIND=JPRB) :: RTICECU
REAL(KIND=JPRB) :: RTWAT_RTICE_R
REAL(KIND=JPRB) :: RTWAT_RTICECU_R
REAL(KIND=JPRB) :: RKOOP1
REAL(KIND=JPRB) :: RKOOP2

!     J.-J. MORCRETTE                   91/07/14  ADAPTED TO I.F.S.

!      NAME     TYPE      PURPOSE
!      ----     ----      -------

!     *R__ES*   REAL      *CONSTANTS USED FOR COMPUTATION OF SATURATION
!                         MIXING RATIO OVER LIQUID WATER(*R_LES*) OR
!                         ICE(*R_IES*).
!     *RVTMP2*  REAL      *RVTMP2=RCPV/RCPD-1.
!     *RHOH2O*  REAL      *DENSITY OF LIQUID WATER.   (RATM/100.)
!     *R5ALVCP* REAL      *R5LES*RLVTT/RCPD
!     *R5ALSCP* REAL      *R5IES*RLSTT/RCPD
!     *RALVDCP* REAL      *RLVTT/RCPD
!     *RALSDCP* REAL      *RLSTT/RCPD
!     *RALFDCP* REAL      *RLMLT/RCPD
!     *RTWAT*   REAL      *RTWAT=RTT
!     *RTBER*   REAL      *RTBER=RTT-0.05
!     *RTBERCU  REAL      *RTBERCU=RTT-5.0
!     *RTICE*   REAL      *RTICE=RTT-0.1
!     *RTICECU* REAL      *RTICECU=RTT-23.0
!     *RKOOP?   REAL      *CONSTANTS TO DESCRIBE KOOP FORM FOR NUCLEATION
!     *RTWAT_RTICE_R*   REAL      *RTWAT_RTICE_R=1./(RTWAT-RTICE)
!     *RTWAT_RTICECU_R* REAL      *RTWAT_RTICECU_R=1./(RTWAT-RTICECU)

!       ----------------------------------------------------------------
END MODULE YOETHF
