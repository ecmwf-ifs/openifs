! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      MODULE YOWTEMP

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      IMPLICIT NONE
!!!!
!!!! this module should now be obsolete !!!!!
!!!!

!*    ** *TEMP* TEMPERATURE GRID

      INTEGER(KIND=JWIM), PARAMETER :: NXT=720 
      INTEGER(KIND=JWIM), PARAMETER :: NYT=325 

      CHARACTER(LEN=14)  :: CTDATE 
      INTEGER(KIND=JWIM) :: NGXT 
      INTEGER(KIND=JWIM) :: NGYT 
      REAL(KIND=JWRB)    :: TMOWEP 
      REAL(KIND=JWRB)    :: TMOSOP 
      REAL(KIND=JWRB)    :: TMOEAP 
      REAL(KIND=JWRB)    :: TMONOP 
      REAL(KIND=JWRB)    :: DELLA 
      REAL(KIND=JWRB)    :: DELLO 

!-----------------------------------------------------------------------

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *CTDATE*    CHAR*14   DATE OF TEMPERATURE GRID
!      *NGXT*      INTEGER   NUMBER OF LONGITUDES IN T-GRID.
!      *NGYT*      INTEGER   NUMBER OF LATITUDES  IN T-GRID.
!      *TMOEAP*    REAL      MOST EASTERN LONGITUDE OF T-GRID (DEGREE).
!      *TMONOP*    REAL      MOST NORTHERN LATITUDE OF T-GRID (DEGREE).
!      *TMOSOP*    REAL      MOST SOUTHERN LATITUDE OF T-GRID (DEGREE).
!      *TMOWEP*    REAL      MOST WESTERN LONGITUDE OF T-GRID (DEGREE).
!      *DELLA*     REAL      INCREMENT FOR LATITUDE IN TEMP-GRID (DEG).
!      *DELLO*     REAL      INCREMENT FOR LONGITUDE IN TEMP-GRID (DEG).

!-----------------------------------------------------------------------

      END MODULE YOWTEMP
