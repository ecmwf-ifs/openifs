! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      MODULE YOWCURG

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      IMPLICIT NONE

!*    ** *CURGRD* - INPUT GRID CURRENT SPECFICATIONS.

      REAL(KIND=JWRB) :: DUCLO 
      REAL(KIND=JWRB) :: DUCLA 
      REAL(KIND=JWRB) :: UCLO 
      REAL(KIND=JWRB) :: UCLH 
      REAL(KIND=JWRB) :: UCLW 
      REAL(KIND=JWRB) :: UCLE 
      REAL(KIND=JWRB), ALLOCATABLE :: UCUR(:,:)
      REAL(KIND=JWRB), ALLOCATABLE :: VCUR(:,:)
      INTEGER(KIND=JWIM) :: KRCI 
      INTEGER(KIND=JWIM) :: KCCI 
      INTEGER(KIND=JWIM) :: IPCUR 

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *DUCLO*     REAL      STEPSIZE BETWEEN LONGITUDES IN DEG.
!      *DUCLA*     REAL      STEPSIZE BETWEEN LATITUDES  IN DEG.
!      *UCLO*      REAL      MOST SOUTHERN LATITUDE.
!      *UCLH*      REAL      MOST NORTHERN LATITUDE.
!      *UCLW*      REAL      LEFT MOST LONGITUDE.
!      *UCLE*      REAL      RIGHT MOST LONGITUDE.
!      *UCUR       REAL      U - GRIDDED COMPONENT OF CURRENT (M/S).
!      *VCUR       REAL      V - GRIDDED COMPONENT OF CURRENT (M/S).
!      *KCCI*      INTEGER   NUMBER OF COLUMNES IN CURRENT INPUT (USED).
!      *KRCI*      INTEGER   NUMBER OF ROWS     IN CURRENT INPUT (USED).
!      *IPCUR*     INTEGER   INDICATOR PERIODICAL(GLOBAL) GRID OR NOT
!                            0 = NON-PERIODICAL;  1 = PERIODICAL
! ----------------------------------------------------------------------
      END MODULE YOWCURG
