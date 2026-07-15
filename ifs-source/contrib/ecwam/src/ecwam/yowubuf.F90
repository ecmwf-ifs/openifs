! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      MODULE YOWUBUF

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      IMPLICIT NONE

!*    ** *UBUF*  GRID POINT DEPENDENT CONSTANTS

      INTEGER(KIND=JWIM), ALLOCATABLE :: KLAT(:,:,:) 
      INTEGER(KIND=JWIM), ALLOCATABLE :: KLON(:,:) 
      INTEGER(KIND=JWIM), ALLOCATABLE :: KCOR(:,:,:) 
      INTEGER(KIND=JWIM), ALLOCATABLE :: KRLAT(:,:,:) 
      INTEGER(KIND=JWIM), ALLOCATABLE :: KRLON(:,:,:) 
      INTEGER(KIND=JWIM), ALLOCATABLE :: KPM(:,:) 
      INTEGER(KIND=JWIM), ALLOCATABLE :: MPM(:,:) 
      INTEGER(KIND=JWIM), ALLOCATABLE :: JXO(:,:) 
      INTEGER(KIND=JWIM), ALLOCATABLE :: JYO(:,:) 
      INTEGER(KIND=JWIM), ALLOCATABLE :: KCR(:,:) 

      REAL(KIND=JWRB), ALLOCATABLE :: WLAT(:,:) 
      REAL(KIND=JWRB), ALLOCATABLE :: WCOR(:,:) 
      REAL(KIND=JWRB), ALLOCATABLE :: SUMWN(:,:,:) 
      REAL(KIND=JWRB), ALLOCATABLE :: WLATN(:,:,:,:,:) 
      REAL(KIND=JWRB), ALLOCATABLE :: WLONN(:,:,:,:) 
      REAL(KIND=JWRB), ALLOCATABLE :: WCORN(:,:,:,:,:) 
      REAL(KIND=JWRB), ALLOCATABLE :: WKPMN(:,:,:,:) 
      REAL(KIND=JWRB), ALLOCATABLE :: WMPMN(:,:,:,:) 
      REAL(KIND=JWRB), ALLOCATABLE :: WRLAT(:,:) 
      REAL(KIND=JWRB), ALLOCATABLE :: WRLON(:,:) 
      REAL(KIND=JWRB), ALLOCATABLE :: OBSLAT(:,:,:) 
      REAL(KIND=JWRB), ALLOCATABLE :: OBSLON(:,:,:) 
      REAL(KIND=JWRB), ALLOCATABLE :: OBSCOR(:,:,:) 
      REAL(KIND=JWRB), ALLOCATABLE :: OBSRLAT(:,:,:) 
      REAL(KIND=JWRB), ALLOCATABLE :: OBSRLON(:,:,:) 

      LOGICAL :: LUPDTWGHT
      LOGICAL, ALLOCATABLE :: LLWLATN(:,:,:,:) 
      LOGICAL, ALLOCATABLE :: LLWLONN(:,:,:) 
      LOGICAL, ALLOCATABLE :: LLWCORN(:,:,:,:) 
      LOGICAL, ALLOCATABLE :: LLWKPMN(:,:,:) 
      LOGICAL, ALLOCATABLE :: LLWMPMN(:,:,:) 

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *KLAT*      INTEGER   KLAT(:,:,1) INDEX OF GRIDPOINT SOUTH AND NORTH
!                            LANDPOINTS ARE MARKED BY ZERO OR NINF+1.
!                            KLAT(:,:,2) INDEX OF 2nd closest GRIDPOINT
!                            SOUTH AND NORTH LANDPOINTS ARE MARKED BY ZERO
!                            OR NINF+1.
!                            !!!!!! THE SECOND INDEX IS 1 FOR SOUTH 
!                                                       2 FOR NORTH
!      *KLON*      INTEGER   INDEX OF GRIDPOINT WEST AND EAST
!                            LANDPOINTS ARE MARKED BY ZERO OR NINF+1.
!                            !!!! THE SECOND INDEX IS 1 FOR WEST 
!                                                     2 FOR EAST 
!      *KCOR*      INTEGER   INDEX OF THE GRID CORNER POINTS IN DIAGONAL
!                            DIRECTIONS. IT IS USED IN MPDECOMP TO DETERMINE
!                            THE MAXIMUM SIZE OF THE HALO WHEN IPROPAGS=2.
!                            BECAUSE OF THE CTUW SCHEME EACH CORNER HAS TO BE
!                            SPECIFIED BY ITS CLOSEST GRID POINTS AND THE POINT
!                            ON EACH SIDE OF IT 
!                            THE SECOND INDEX IS 1 FOR NORTH EAST CORNER
!                                                2 FOR SOUTH EAST CORNER
!                                                3 FOR SOUTH WEST CORNER
!                                                4 FOR NORTH WEST CORNER
!                            THE THIRD INDEX IS 1 FOR THE CLOSEST GRID POINT
!                                               2 FOR THE SECOND CLOSEST GRID POINT
!   *KRLAT(:,:,1)*INTEGER   INDEX OF GRIDPOINT SOUTH-EAST AND NORTH-WEST
!                            LANDPOINTS ARE MARKED BY ZERO OR NINF+1.
!   *KRLAT(:,:,2)*INTEGER   INDEX OF 2nd closest GRIDPOINT SOUTH-EAST
!                            NORTH-WEST. LANDPOINTS ARE MARKED BY ZERO
!                            OR NINF+1.
!   *KRLON(:,:,1)*INTEGER   INDEX OF GRIDPOINT SOUTH-WEST AND NORTH-EAST
!                            LANDPOINTS ARE MARKED BY ZERO OR NINF+1.
!   *KRLON(:,:,2)*INTEGER   INDEX OF 2nd closest GRIDPOINT SOUTH-WEST
!                            AND NORTH-EAST. LANDPOINTS ARE MARKED BY ZERO
!                            OR NINF+1.
!      *KPM*       INTEGER   INDEX FOR DIRECTION TERMS IN CTUW.
!      *MPM*       INTEGER   INDEX FOR FREQUENCY TERMS IN CTUW.
!      *JXO*       INTEGER   INDEX FOR EAST-WEST TERMS IN CTUW.
!      *JYO*       INTEGER   INDEX FOR NORTH-SOUTH TERMS IN CTUW.
!      *KCR*       INTEGER   INDEX FOR CORNER TERMS IN CTUW.
!      *WLAT*      REAL      WEIGHT IN ADVECTION SCHEME FOR CLOSEST
!                            AND THE SECOND CLOSEST GRIDPOINTS IN THE
!                            NORTH-SOUTH DIRECTION.
!      *WCOR*      REAL      WEIGHT IN ADVECTION SCHEME FOR CLOSEST
!                            AND THE SECOND CLOSEST GRIDPOINTS IN THE
!                            CORNER DIRECTIONS (SEE KCOR).
!      *SUMWN*     REAL      SUM OF WEIGHTS IN ADVECTION CTU SCHEME
!      *WLATN*     REAL      WEIGHTS IN ADVECTION CTU SCHEME
!                            FOR BOTH NORTH-SOUTH DIRECTIONS.
!                            WLATN(IJ,K,M,IC,ICL)
!                            IJ: GRID POINT
!                            K : DIRECTION
!                            M : FREQUENCY
!                            IC : NORTH SOUTH INDEX
!                            ICL : 1 FOR THE CLOSEST GRID POINT TO THE CORNER
!                                  2 FOR THE SECOND CLOSEST GRID POINT TO THE CORNER
!      *WLONN*     REAL      WEIGHTS IN ADVECTION CTU SCHEME
!                            FOR BOTH EAST-WEST DIRECTIONS.
!                            WLONN(IJ,K,M,IC)
!                            IJ: GRID POINT
!                            K : DIRECTION
!                            M : FREQUENCY
!                            IC : EAST WEST INDEX
!      *WCORN*     REAL      WEIGHTS IN ADVECTION CTU SCHEME
!                            FOR THE 4 CORNER DIRECTIONS.
!                            WCORN(IJ,K,M,ICR,ICL)
!                            IJ: GRID POINT
!                            K : DIRECTION
!                            M : FREQUENCY
!                            ICR : CORNER INDEX:
!                                1: CORNER CORRESPONDING TO UPWIND CORNER IF NO CURRENTS.
!                                2: CORNER CORRESPONDING TO THE CORNER ON THE SAME LATITUDE AS THE UPWIND CORNER IF NO CURRENTS.
!                                3: CORNER CORRESPONDING TO THE CORNER ON THE SAME LONGITUDE AS THE UPWIND CORNER IF NO CURRENTS.
!                                4: CORNER CORRESPONDING TO THE OPPOSITE CORNER TO THE UPWIND CORNER IF NO CURRENTS.
!                            ICL : 1 FOR THE CLOSEST GRID POINT TO THE CORNER
!                                  2 FOR THE SECOND CLOSEST GRID POINT TO THE CORNER
!
!      *WKPMN*     REAL      WEIGHTS WHEN USING ADVECTION CTU SCHEME
!                            FOR REFRACTION TERM K+-1.
!      *WMPMN*     REAL      WEIGHTS WHEN USING ADVECTION CTU SCHEME
!                            FOR FREQUENCY SHIFTING TERM M+-1.
!      *WRLAT*     REAL      WEIGHT IN ADVECTION SCHEME FOR CLOSEST
!                            AND THE SECOND CLOSEST GRIDPOINTS IN THE
!                            SOUTH-EAST AND NORTH-WEST DIRECTION.
!      *WRLON*     REAL      WEIGHT IN ADVECTION SCHEME FOR CLOSEST
!                            AND THE SECOND CLOSEST GRIDPOINTS IN THE
!                            SOUTH-WEST AND NORTH-EAST DIRECTION.
!      *OBSLAT*    REAL      TRANSMISSION COEFFICIENT DUE TO OBSTRUCTIONS
!                            IN THE MERIDIONAL DIRECTION. 
!      *OBSLON*    REAL      TRANSMISSION COEFFICIENT DUE TO OBSTRUCTIONS
!                            IN THE ZONAL DIRECTION. 
!      *OBSCOR*    REAL      TRANSMISSION COEFFICIENT DUE TO OBSTRUCTIONS
!                            IN THE GRID CORNER DIRECTIONS. 
!      *OBSRLAT*   REAL      TRANSMISSION COEFFICIENT DUE TO OBSTRUCTIONS
!                            IN THE SOUTH-EAST AND NORTH-WEST DIRECTION. 
!      *OBSRLON*   REAL      TRANSMISSION COEFFICIENT DUE TO OBSTRUCTIONS
!                            IN THE SOUTH-WEST AND NORTH-EAST DIRECTION. 
!      *LUPDTWGHT* LOGICAL   TRUE IF CTUW HAS TO BE CALLED
!      *LLWLATN*   LOGICAL ARRAY, TRUE IF WLATN > 0. AT ALL GRID POINTS. 
!      *LLWLONN*   LOGICAL ARRAY, TRUE IF WLONN > 0. AT ALL GRID POINTS.
!      *LLWCORN*   LOGICAL ARRAY, TRUE IF WCORN > 0. AT ALL GRID POINTS.
!      *LLWKPMN*   LOGICAL ARRAY, TRUE IF WKPNM > 0. AT ALL GRID POINTS.
!      *LLWMPMN*   LOGICAL ARRAY, TRUE IF WMPMN > 0. AT ALL GRID POINTS.

! ----------------------------------------------------------------------
      END MODULE YOWUBUF
