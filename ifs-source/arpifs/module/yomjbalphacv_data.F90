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

MODULE YOMJBALPHACV_DATA

USE PARKIND1, ONLY: JPIM, JPRB

IMPLICIT NONE

SAVE

!*
!      /YOMJBALPHACV_DATA/ - ALPHA CONTROL VARIABLE RELEVANT PARAMETERS

!     SEBASTIEN MASSART     *ECMWF*

!     --------- definition of ACV through namelist NAMJBALPHACV --------------------

TYPE TJBALPHACV_CONFIG
  LOGICAL                         :: LALPHACV_TAPER_HUM
  LOGICAL                         :: LALPHACV_UNBAL
  LOGICAL                         :: LALPHACV_CVHUM
  LOGICAL                         :: LFILTER_ALPHA
  LOGICAL                         :: L_OZONE_IN_ALPHACV
  INTEGER(KIND=JPIM)              :: NFILTER_EXP
  INTEGER(KIND=JPIM)              :: NLOC_LAP_HUM 
  INTEGER(KIND=JPIM)              :: NLOC_LAP_LNSP
  INTEGER(KIND=JPIM)              :: NLOC_LAP_VOD
  INTEGER(KIND=JPIM)              :: NLOC_VSCALE
  INTEGER(KIND=JPIM)              :: NMEMBERS_ALPHACV
  INTEGER(KIND=JPIM)              :: NTIMES_ALPHACV
  INTEGER(KIND=JPIM)              :: NWEIGHTE_LAYER
  REAL(KIND=JPRB)                 :: REDNMC_ALPHACV
  REAL(KIND=JPRB)                 :: RFILTER_SCALE
  REAL(KIND=JPRB)                 :: R_STATIC_WEIGHT
END TYPE TJBALPHACV_CONFIG

TYPE :: TJBHYBINACV
  INTEGER(KIND=JPIM)              :: NINTERP_PERT = 4
  INTEGER(KIND=JPIM)              :: NPERT_ALPHACV
  INTEGER(KIND=JPIM)              :: NDIM2_ALPHACV
  INTEGER(KIND=JPIM), ALLOCATABLE :: MHYBRIDVAR(:)
END TYPE TJBHYBINACV


TYPE(TJBHYBINACV)               :: YRJBALPHACV
TYPE(TJBALPHACV_CONFIG)         :: YRJBALPHACV_CONF


END MODULE YOMJBALPHACV_DATA
