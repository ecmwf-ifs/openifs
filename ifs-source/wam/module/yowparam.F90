! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

      MODULE YOWPARAM

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      IMPLICIT NONE

      INTEGER(KIND=JWIM), PARAMETER :: IMDLGRDID=218
      INTEGER(KIND=JWIM), PARAMETER :: KWAMVER=8

!*    ** *PARAM*  GENERAL GRID SIZE INFORMATION.

      INTEGER(KIND=JWIM) :: NANG 
      INTEGER(KIND=JWIM) :: NFRE
      INTEGER(KIND=JWIM) :: NFRE_RED
      INTEGER(KIND=JWIM) :: NFRE_ODD
      INTEGER(KIND=JWIM) :: NGX
      INTEGER(KIND=JWIM) :: NGY
      INTEGER(KIND=JWIM) :: NIBLO
      INTEGER(KIND=JWIM) :: NOVER 
      INTEGER(KIND=JWIM) :: NIBL1 
      INTEGER(KIND=JWIM) :: NIBLD 
      INTEGER(KIND=JWIM) :: NIBLC 

      REAL(KIND=JWRB) :: SWAMPWIND 
      REAL(KIND=JWRB) :: SWAMPWIND2 
      REAL(KIND=JWRB) :: DTNEWWIND
      REAL(KIND=JWRB) :: SWAMPCIFR
      REAL(KIND=JWRB) :: SWAMPCITH

      CHARACTER(LEN=1)   :: CLDOMAIN

      LOGICAL :: LTURN90
      LOGICAL :: LL1D
      LOGICAL :: LWDINTS

!*     VARIABLE.   TYPE.     PURPOSE.
!       ---------   -------   --------
!      *IMDLGRDID* INTEGER   MODEL GRID IDENTIFICATION NUMBER.
!                            IT SHOULD BE UPDATED ONLY IF NEW gridglou
!                            AND ubufglou FILES ARE PRODUCED!!!!
!      *KWAMVER*   INTEGER   WAM SOFTWARE VERSION NUMBER.
!      *NANG*      INTEGER   NUMBER OF ANGLES.
!      *NFRE*      INTEGER   NUMBER OF FREQUENCIES FOR THE PHYSICS.
!      *NFRE_RED*  INTEGER   REDUCED NUMBER OF FREQUENCIES FOR THE PROPAGATION AND IO
!                            BY DEFAULT = NFRE
!      *NFRE_ODD*  INTEGER   NFRE-1 IF NFRE EVEN, BUT NFRE IS NFRE ODD 
!      *NGX*       INTEGER   NUMBER OF LONGITUDES IN GRID.
!      *NGY*       INTEGER   NUMBER OF LATITUDES IN GRID.
!      *NIBLO*     INTEGER   NUMBER OF SEA POINTS IN BLOCK.
!      *NOVER*     INTEGER   MAXIMUM NUMBER POINTS IN FIRST LATITUDE
!                            OF BLOCKS.
!      *NIBL1*     INTEGER   = NIBLO IF MULTI BLOCK VERSION.
!                            =     1 IF ONE BLOCK VERSION.
!      *NIBLD*     INTEGER   = NIBLO IF DEPTH OR CURRENT REFRACTION.
!                            = 1     ELSE.
!      *NIBLC*     INTEGER   = NIBLO IF CURRENT REFRACTION.
!                            = 1     ELSE.
!      *SWAMPWIND* REAL      CONSTANT WIND SPEED USED TO RUN SWAMP CASE.
!                            FIRST VALUE
!      *SWAMPWIND2*REAL      CONSTANT WIND SPEED USED TO RUN SWAMP CASE.
!                            SECOND VALUE PROVIDED IT'S NOT ZERO
!      *DTNEWWIND* REAL      TIME AFTER WHICH SWAMPWIND2 IS APPLIED
!                            IN HOURS
!      *SWAMPCIFR* REAL      SEA ICE COVER FOR THE NORTHERN HALF OF THE SWAMP DOMAIN 
!      *SWAMPCITH* REAL      SEA ICE THICKNESS FOR THE NORTHERN HALF OF THE SWAMP DOMAIN 
!      *LTURN90*   LOGICAL   IF TRUE THE NEW SWAMP WIND WILL BE TURNED
!                            BY 90 DEGREES.
!      *LWDINTS*   LOGICAL   IF TRUE A FILE CONTAINING WIND SPEED (m/s) AND DIRECTION
!                            (in degrees following the meteorological convention)
!                            TIME SERIES WILL BY USED TO SET THE WIND FORCING
!                            OVER THE ALL DOMAIN.
!                            THE PRESENCE OF THE FILE, windforcing_time_series
!                            WILL DETERMINE WHETHER OR NOT LWDINTS IS TRUE OR NOT
!      *CLDOMAIN*  CHARACTER DEFINES THE DOMAIN OF THE MODEL (for the
!                            FDB and for selection of some variables)
!      *LL1D*      LOGICAL   IF TRUE THEN THE DOMAIN DECOMPOSITION IS ONLY
!                            DONE IN LATITUNAL BANDS
!                            (like it used to be done).
! ----------------------------------------------------------------------

      END MODULE YOWPARAM
