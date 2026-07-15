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

MODULE PARDIMO

USE PARKIND1  ,ONLY : JPIM

IMPLICIT NONE

SAVE

!****-------------------------------------------------------------------
!****  CD PARDIMO : PARAMETERS RELATED TO OBSERVATIONS
!****  ----------
!****  Auteurs    : J.PAILLEUX, D.VASILJEVIC, P.CAILLE    90/05-91/09
!****-------------------------------------------------------------------
!*  OBSERVATIONS ARRAYS DIMENSIONS
!*  ------------------------------
!*    JPNOTP  : MAXIMUM NUMBER OF OBSERVATIONS TYPES
!*    JPXTIM  : MAXIMUM NUMBER OF TIME SLOTS
!*
!*  OBSERVATIONS PROCESSING
!*  -----------------------
!*    JPXSOBT : MAXIMUM NUMBER OF SUB-OBSERVATIONS TYPES
!*    JPXVAR  : MAXIMUM NUMBER OF VARIABLES
!*    JPXAREA : MAXIMUM NUMBER OF AREAS
!*    JPHUBLEV: NUMBER OF DIFFERENT HUBER LAYERS IN THE VERTICAL
!*              (STRATOSPHERE=1,TROPOSPHERE=2,BOUNDARY LAYER=3) 
!*    JPXAMVPROD: MAXIMUM NUMBER OF AMV(SATAM) PRODUCERS
!*
!    N. Semane     24-Oct-2022 Add JPGNSSRO_SAT and JPGNSSRO_TYP dimensions
!****-------------------------------------------------------------------
INTEGER(KIND=JPIM), PARAMETER :: JPNOTP=19
INTEGER(KIND=JPIM), PARAMETER :: JPXTIM=99

INTEGER(KIND=JPIM), PARAMETER :: JPXSOBT=3
INTEGER(KIND=JPIM), PARAMETER :: JPXAREA=20

INTEGER(KIND=JPIM), PARAMETER :: JPHUBLEV=3
INTEGER(KIND=JPIM), PARAMETER :: JPXAMVPROD=10
INTEGER(KIND=JPIM), PARAMETER :: JPGNSSRO_SAT = 99 ! MAXIMUM NUMBER OF GNSSRO SATELLITES= Max_gnssro_reportype (7XXX) minus 7000
INTEGER(KIND=JPIM), PARAMETER :: JPGNSSRO_TYP = 9 ! MAXIMUM NUMBER OF GNSSRO TYPES = Max_gnssro_type (40X) minus 400

!****-------------------------------------------------------------------
END MODULE PARDIMO
