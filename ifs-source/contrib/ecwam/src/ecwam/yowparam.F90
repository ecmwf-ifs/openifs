! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      MODULE YOWPARAM

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      IMPLICIT NONE

      INTEGER(KIND=JWIM), PARAMETER :: KWAMVER=9

!*    ** *PARAM*  GENERAL GRID SIZE INFORMATION.

      INTEGER(KIND=JWIM) :: NANG 
      INTEGER(KIND=JWIM) :: NFRE
      INTEGER(KIND=JWIM) :: NFRE_RED
      INTEGER(KIND=JWIM) :: NFRE_ODD

      REAL(KIND=JWRB) :: SWAMPWIND 
      REAL(KIND=JWRB) :: SWAMPWIND2 
      REAL(KIND=JWRB) :: DTNEWWIND
      REAL(KIND=JWRB) :: SWAMPCIFR
      REAL(KIND=JWRB) :: SWAMPCITH

      LOGICAL :: LTURN90
      LOGICAL :: LL1D
      LOGICAL :: LWDINTS

      ! Moved here from yowunpool
      LOGICAL :: LLUNSTR
      LOGICAL :: LLR8TOR4 = .FALSE. ! Is input in fact legacy REAL*8 & NKIND == 4 ?

!*     VARIABLE.   TYPE.     PURPOSE.
!       ---------   -------   --------
!      *KWAMVER*   INTEGER   WAM SOFTWARE VERSION NUMBER.
!      *NANG*      INTEGER   NUMBER OF ANGLES.
!      *NFRE*      INTEGER   NUMBER OF FREQUENCIES FOR THE PHYSICS.
!      *NFRE_RED*  INTEGER   REDUCED NUMBER OF FREQUENCIES FOR THE PROPAGATION AND IO
!                            BY DEFAULT = NFRE
!      *NFRE_ODD*  INTEGER   NFRE-1 IF NFRE EVEN, BUT NFRE IS NFRE ODD 
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
!      *LL1D*      LOGICAL   IF TRUE THEN THE DOMAIN DECOMPOSITION IS ONLY
!                            DONE IN LATITUNAL BANDS
!                            (like it used to be done).
! ----------------------------------------------------------------------

      END MODULE YOWPARAM
