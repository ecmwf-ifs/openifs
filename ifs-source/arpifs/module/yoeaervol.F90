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

MODULE YOEAERVOL

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------
!*    ** *YOEAERVOL* - CONTROL PARAMETERS FOR VOLCANIC AEROSOLS IN THE ATMOSPHERE
!     ------------------------------------------------------------------

TYPE :: TEAERVOL
INTEGER(KIND=JPIM) :: NAERVOLC
INTEGER(KIND=JPIM) :: NAERVOLE
INTEGER(KIND=JPIM) :: NINTERPT
INTEGER(KIND=JPIM) :: NTVOLC
INTEGER(KIND=JPIM) :: NVOLERUP
INTEGER(KIND=JPIM) :: NVOLHOMO
INTEGER(KIND=JPIM) :: NVOLOPTP
INTEGER(KIND=JPIM) :: NVOLERUZ(100)
INTEGER(KIND=JPIM) :: NVOLDATS(100)
INTEGER(KIND=JPIM) :: NVOLDATE(100)

REAL(KIND=JPRB) :: RAERVOLC(100,8)
REAL(KIND=JPRB) :: RAERVOLE(100,8)
REAL(KIND=JPRB) :: RVOLERUZ(100)
REAL(KIND=JPRB) :: RVSEDVOL, RVDEPVOL, RWDEPVOL
!---------------------------------------------------------------------
CONTAINS
  PROCEDURE, PASS :: PRINT => PRINT_CONFIGURATION 
END TYPE TEAERVOL
!=====================================================================

TYPE(TEAERVOL), POINTER :: YREAERVOL => NULL()

!     ------------------------------------------------------------------
! AERVOLC  parameters for continuous volcanoes
! AERVOLE  parameters for explosive volcanoes
! NAERVOL.. number of .. volcanoes
! ...VOL.(.. , 1) latitude
! ...VOL.(.. , 2) longitude  
! ...VOL.(.. , 3) mass of particulates (to go into FlyAsh) (kg)
! ...VOL.(.. , 4) mass of gases (to go into SO2) (kg)
! ...VOL.(.. , 5) base of drifting plume  (m)
! ...VOL.(.. , 6) top of drifting plume   (m)
! NINTERPT = 0 takes latest profile, = 1 interpolates between profiles  
! NTVOLC    total number of volcanoes considered (continuous + explosive)
! NVOLERUP = 0 no erupting volcano
! NVOLERUP = 1 using AEROCOM-type data for location and ejected mass
! NVOLERUP = 2 erupting volcano has some input as f(time/height)
! NVOLERUZ = 0 for erupting volcano, homogeneous distribution of mass over the vertical
! NVOLERUZ = 1 interpolating the "observed/retrieved" mass over the vertical
! RVOLERUZ     emission weighting factor when emission profile is known (empirical)
! NVOLHOMO = 0 profile information is used, = 1 vertically homogeneous distribution 
! NVOLOPTP     index for choosing the optical properties assigned to the ash
! NVOLOPTP = 1 opt.prop. as SO4
! NVOLOPTP = 2 opt.prop. as BC 
! NVOLOPTP = 3 opt.prop. as DU3
! NVOLDATS     yyyymmddhh date and start of eruption
! RVSEDVOL     gravitational sedimentation coefficient for volcanic component
! RVDEPVOL     dry deposition coefficient for volcanic component
! RWDEPVOL     wet deposition coefficient for volcanic component
!     ------------------------------------------------------------------

CONTAINS

SUBROUTINE PRINT_CONFIGURATION(SELF, KDEPTH, KOUTNO)
IMPLICIT NONE
CLASS(TEAERVOL), INTENT(IN) :: SELF
INTEGER        , INTENT(IN) :: KDEPTH
INTEGER        , INTENT(IN) :: KOUTNO

INTEGER :: IDEPTHLOC

IDEPTHLOC = KDEPTH+2

WRITE(KOUTNO,*) REPEAT(' ',KDEPTH   ) // 'model%yrml_phy_aer%yreaervol : '
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NAERVOLC = ', SELF%NAERVOLC
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NAERVOLE = ', SELF%NAERVOLE
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NINTERPT = ', SELF%NINTERPT
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NTVOLC = ', SELF%NTVOLC
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NVOLERUP = ', SELF%NVOLERUP
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NVOLHOMO = ', SELF%NVOLHOMO
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NVOLOPTP = ', SELF%NVOLOPTP
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NVOLERUZ SUM = ', SUM(SELF%NVOLERUZ)
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NVOLDATS SUM = ', SUM(SELF%NVOLDATS)
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RAERVOLC SUM = ', SUM(SELF%RAERVOLC)
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RAERVOLE SUM = ', SUM(SELF%RAERVOLE)
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RVOLERUZ SUM = ', SUM(SELF%RVOLERUZ)
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RVSEDVOL = ', SELF%RVSEDVOL
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RVDEPVOL = ', SELF%RVDEPVOL
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RWDEPVOL = ', SELF%RWDEPVOL

END SUBROUTINE PRINT_CONFIGURATION

END MODULE YOEAERVOL
