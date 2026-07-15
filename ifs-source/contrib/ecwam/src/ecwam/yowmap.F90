! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      MODULE YOWMAP

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWDRVTYPE  , ONLY : WVGRIDGLO, WVGRIDLOC

      IMPLICIT NONE

!*    ** *MAP*  LON/LAT INDEX OF EACH SEA POINT.

      INTEGER(KIND=JWIM)              :: NGX 
      INTEGER(KIND=JWIM)              :: NGY 
      INTEGER(KIND=JWIM)              :: NIBLO
      INTEGER(KIND=JWIM)              :: IPER 
      INTEGER(KIND=JWIM)              :: IRGG 
      INTEGER(KIND=JWIM)              :: IQGAUSS
      INTEGER(KIND=JWIM)              :: KMNOP
      INTEGER(KIND=JWIM)              :: KMSOP
      INTEGER(KIND=JWIM), PARAMETER   :: MIKOFST=4
      INTEGER(KIND=JWIM), ALLOCATABLE :: NLONRGG(:)
      INTEGER(KIND=JWIM), ALLOCATABLE :: KXLTMIN(:) 
      INTEGER(KIND=JWIM), ALLOCATABLE :: KXLTMAX(:) 
      TYPE(WVGRIDGLO) :: BLK2GLO
      TYPE(WVGRIDLOC) :: BLK2LOC

      REAL(KIND=JWRB)                 :: AMOWEP 
      REAL(KIND=JWRB)                 :: AMOSOP 
      REAL(KIND=JWRB)                 :: AMOEAP 
      REAL(KIND=JWRB)                 :: AMONOP 
      REAL(KIND=JWRB)                 :: XDELLA 
      REAL(KIND=JWRB)                 :: XDELLO 
      REAL(KIND=JWRB), ALLOCATABLE    :: ZDELLO(:) 

      CHARACTER(LEN=1)   :: CLDOMAIN

      LOGICAL              :: LLOBSTRCT
      LOGICAL              :: LAQUA

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *NGX*       INTEGER   NUMBER OF LONGITUDES IN GRID.
!      *NGY*       INTEGER   NUMBER OF LATITUDES  IN GRID.
!      *NIBLO*     INTEGER   NUMBER OF SEA POINTS.
!      *IPER*      INTEGER   = 1 IF GRID IS PERIODIC.
!      *IRGG*      INTEGER   GRID CODE: 0 = REGULAR, 1 = IRREGULAR.
!      *IQGAUSS*   INTEGER   =1 IF A QUASI GAUSSIAN GRID IS USED. 
!                            =0 OTHERWISE.
!      *KMNOP*     INTEGER   INDEX FOR THE MOST NORTHERN LATITUDE WITH
!                            MODEL SEA POINTS. IT SHOULD NORMALLY BE 1
!                            EXCEPT IF THE NORHTERN LATITUDES HAVES BEEN
!                            BLANKED (SEE QUASI GAUSSIAN GRID)
!      *KMSOP*     INTEGER   INDEX FOR THE MOST SOUTHERN LATITUDE WITH
!                            MODEL SEA POINTS. IT SHOULD NORMALLY BE NGY 
!      *MIKOFST*   INTEGER   PERIODIC OFFSET FOR ARRAY IJFROMIK USED IN CTUW.
!      *NLONRGG*   INTEGER   NUMBER OF GRID POINTS PER LATITUDES.
!      *KXLTMIN*   INTEGER   MINIMUM OF KXLT ON EACH PE.
!      *KXLTMAX*   INTEGER   MAXIMUM OF KXLT ON EACH PE.
!      *BLK2GLO*             TRANSFORMATION FROM BLOCK TO GLOBAL GRID POINTS (structured grid)
!      *BLK2LOC*             TRANSFORMATION FROM BLOCK TO LOCAL GRID POINTS (on same MPI task) (structured grid) 
!      *AMOWEP*    REAL      MOST WESTERN LONGITUDE IN GRID (DEGREE).
!      *AMOSOP*    REAL      MOST SOUTHERN LATITUDE IN GRID (DEGREE).
!      *AMOEAP*    REAL      MOST EASTERN LONGITUDE IN GRID (DEGREE).
!      *AMONOP*    REAL      MOST NORTHERN LATITUDE IN GRID (DEGREE).
!      *XDELLA*    REAL      GRID INCREMENT FOR LATITUDE (DEGREE).
!      *XDELLO*    REAL      CONSTANT GRID INCREMENT FOR LONGITUDE (DEG)
!      *ZDELLO*    REAL      VARIABLE GRID INCREMENT FOR LONGITUDE (DEG)
!      *CLDOMAIN*  CHARACTER DEFINES THE DOMAIN OF THE MODEL (for the
!                            FDB and for selection of some variables)
!      *LLOBSTRCT* LOGICAL   CONTROLS WHETHER THE NEW TYPE OF BATHYMETRY
!                            INPUT IS USED IN CONJUNCTION WITH THE
!                            OBSTRUCTION COEFFICIENT IN THE ADVECTION.
!      *LAQUA*     LOGICAL   TRUE IF MODEL IN AQUA PLANET CONFIGURATION
!                            (i.e. NO LAND AND DEEP WATER).

! ----------------------------------------------------------------------
      END MODULE YOWMAP
