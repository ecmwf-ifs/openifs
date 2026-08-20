! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      MODULE YOWSHAL

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWDRVTYPE  , ONLY : ENVIRONMENT, FREQUENCY

      IMPLICIT NONE

!*    ** *SHALLOW*   SHALLOW WATER TABLES and FEATURES.

      INTEGER(KIND=JWIM)              :: NDEPTH

      REAL(KIND=JWRB), PARAMETER   :: GAM_B_J = 0.8_JWRB
      REAL(KIND=JWRB), PARAMETER   :: BATHYMAX = 999.0_JWRB
      REAL(KIND=JWRB), PARAMETER   :: XKDMIN = 0.75_JWRB 

      REAL(KIND=JWRB)              :: DEPTHA
      REAL(KIND=JWRB)              :: DEPTHD
      REAL(KIND=JWRB)              :: TOOSHALLOW

      REAL(KIND=JWRB), ALLOCATABLE :: BATHY(:,:)   !!! should be removed as soon as possible
      LOGICAL, ALLOCATABLE         :: LLOCEANMASK(:,:) 

      REAL(KIND=JWRB), ALLOCATABLE :: TCGOND(:,:)
      REAL(KIND=JWRB), ALLOCATABLE :: TFAK(:,:)
      REAL(KIND=JWRB), ALLOCATABLE :: TSIHKD(:,:)
      REAL(KIND=JWRB), ALLOCATABLE :: TFAC_ST(:,:)

!*    ** GRID POINT FIELDS **
      TYPE(ENVIRONMENT) :: WVENVI

!*    ** GRID POINT AND FREQUENCY FIELDS **
      TYPE(FREQUENCY) :: WVPRPT

!*   ** FICTIOUS VALUE FOR LAND POINT (NSUP+1)

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *NDEPTH*    INTEGER   LENGTH OF SHALLOW WATER TABLES (see MTABS).
!      *GAM_B_J    REAL      FACTIO OF THE DEPTH THAT DETERMINES THE MAXIMUM SIGNIFICANT WAVE HEIGHT ALLOWED
!      *BATHYMAX*  REAL      MAXIMUM DEPTH, ASSUMED TO CORRESPOND TO DEEP WATER CONDITIONS
!      *BATHY*     READ      WATER DEPTH (> 0 FOR OCEAN POINTS
!      *LLOCEANMASK*         LOGICAL, TRUE IF BATHY > 0 (i.e. OCEAN POINTS)
!      *XKDMIN*    REAL      MINIMUM VALUE FOR K*DEPTH IN NON-LINEAR EFFECT FOR FREAK WAVE SOFTWARE. 
!      *DEPTHA*    REAL      MINIMUM DEPTH.
!      *TOOSHALLOW REAL      MINIMUM DEPTH THAT WILL BE ALLOWED
!                            USED TO POINT TO LAND POINTS
!      *DEPTHD*    REAL      DEPTH INCREMENT (METRES) IN TABLES.

!!     GRID POINT FIELDS:
!      -----------------
!      *WVENVI*    ENVIRONMENT IN WHICH WAVES EVOLVE

!!     GRID POINT AND FREQUENCY FIELDS:
!      --------------------------------
!      *WVPRPT*    REAL      WAVE PROPERTIES

! ----------------------------------------------------------------------
      END MODULE YOWSHAL
