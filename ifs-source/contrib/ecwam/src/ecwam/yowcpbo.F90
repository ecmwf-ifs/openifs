! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      MODULE YOWCPBO

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      IMPLICIT NONE

!*    ** *YOWCPBO* USED DURING THE COARSE GRID RUN FOR THE ORGANIZATION OF
!                  THE BOUNDARY POINTS (for the next fine grid(s)).

      INTEGER(KIND=JWIM), PARAMETER :: GBOUNC_MAX=20  ! should not exceed 99 (see USERIN)
      INTEGER(KIND=JWIM) :: IBOUNC
      INTEGER(KIND=JWIM) :: GBOUNC
      INTEGER(KIND=JWIM) :: NBOUNC
      INTEGER(KIND=JWIM), ALLOCATABLE :: IJARC(:)
      INTEGER(KIND=JWIM), ALLOCATABLE :: IGARC(:)  !!! obsolete
      INTEGER(KIND=JWIM), ALLOCATABLE :: IPOGBO(:)
      REAL(KIND=JWRB) :: DLAMAC
      REAL(KIND=JWRB) :: DPHIAC
      REAL(KIND=JWRB) :: AMOSOC(GBOUNC_MAX)
      REAL(KIND=JWRB) :: AMONOC(GBOUNC_MAX)
      REAL(KIND=JWRB) :: AMOEAC(GBOUNC_MAX)
      REAL(KIND=JWRB) :: AMOWEC(GBOUNC_MAX)
      REAL(KIND=JWRB), ALLOCATABLE :: BLATC(:)
      REAL(KIND=JWRB), ALLOCATABLE :: BLNGC(:)

      CHARACTER(LEN=3) :: CBCPREF(GBOUNC_MAX)


!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *GBOUNC_MAX*INTEGER   MAXIMUN NUMBER OF FINE GRIDS THAT COULD BE SPECIFIED.
!      *IBOUNC*    INTEGER   FLAG FOR THE COARSE GRID TO GENERATE BOUNDARY POINTS.
!                            = 1;  THE RUN INCLUDES THE PRODUCTION OF BOUNDARY POINTS.
!                            ELSE; NO BOUNDARY POINTS.
!      *GBOUNC*    INTEGER   ACTUAL NUMBER OF FINE GRIDS.
!      *NBOUNC*    INTEGER   TOTAL NUMBER OF BOUNDARY POINTS (i.e. including all fine grids).
!      *IGARC*     INTEGER   INDEX OF BLOCK FOR A BOUNDARY POINT.
!      *IJARC*     INTEGER   INDEX IN A BLOCK FOR A BOUNDARY POINT.
!      *IPOGBO*    INTEGER   LAST INDEX OF BOUNDARY POINTS FOR EACH FINE GRID.
!      *DLAMAC*    REAL      LONGITUDE INCREMENT OF COARSE GRID (DEG).
!      *DPHIAC*    REAL      LATITUDE INCREMENT OF COARSE GRID (DEG).
!      *AMOWEC*    REAL      MOST EASTERN LONGITUDE FOR THE FINE GRID.
!      *AMOSOC*    REAL      MOST SOUTHERN LONGITUDE FOR THE FINE GRID.
!      *AMOEAC*    REAL      MOST WESTERN LONGITUDE FOR THE FINE GRID.
!      *AMONOC*    REAL      MOST NORTHERN LONGITUDE FOR THE FINE GRID.
!      *BLATC*     REAL      LATITUDE OF COARSE GRID BOUNDARY POINTS.
!      *BLNGC*     REAL      LONGITUDE OF COARSE GRID BOUNDARY POINTS.
!      *CBCPREF*   CHAR*3    CHARACTER STRING INDICATING THE PREFIX
!                            OF BC FILES.

! ----------------------------------------------------------------------
      END MODULE YOWCPBO
