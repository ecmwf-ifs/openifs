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

MODULE YOMGBRAD
      
USE PARKIND1, ONLY: JPIM, JPRB

IMPLICIT NONE

SAVE

REAL (KIND=JPRB) :: RLIMIT_SD  ! obs-grid search distance limit  [m]
REAL (KIND=JPRB) :: REPSLOG    ! offset in rain rate log transform
REAL (KIND=JPRB) :: RR_MIN     ! threshold minimum RR [mm/h]

LOGICAL :: LEGBRAD_ACTIVE

!LLL INTEGER (KIND=JPIM), ALLOCATABLE :: OB_FROM    (:,:)
!LLL INTEGER (KIND=JPIM), ALLOCATABLE :: OB_JOBS    (:,:)
!LLL INTEGER (KIND=JPIM), ALLOCATABLE :: OB_JROF_GP (:,:)
!LLL INTEGER (KIND=JPIM), ALLOCATABLE :: OB_INDEX   (:)
!LLL INTEGER (KIND=JPIM), ALLOCATABLE :: GP_COUNT   (:)
!LLL INTEGER (KIND=JPIM), ALLOCATABLE :: OB_COUNT   (:)

!LLL TYPE GBRAD_PHYS_TYPE
!LLL   INTEGER (KIND=JPIM) :: NLEVELS             ! number of atmospheric levels
!LLL   INTEGER (KIND=JPIM) :: NPROFILES           ! number of profiles
!LLL   REAL    (KIND=JPRB), POINTER :: P   (:)    ! full-level model pressure (hPa)
!LLL   REAL    (KIND=JPRB), POINTER :: PH  (:)    ! half-level model pressure (hPa)
!LLL   REAL    (KIND=JPRB), POINTER :: T   (:)    ! temperature (K)
!LLL   REAL    (KIND=JPRB), POINTER :: Q   (:)    ! specific humidity (kg/kg)
!LLL   REAL    (KIND=JPRB) :: PRECS               ! accumulated surface precipitation (mm/s)
!LLL   REAL    (KIND=JPRB) :: T2M                 ! 2m temperature (K)
!LLL   REAL    (KIND=JPRB) :: LSM                 ! land-sea mask ()
!LLL   REAL    (KIND=JPRB) :: SDO                 ! orography standard deviation (m)
!LLL   REAL    (KIND=JPRB) :: SLAT                ! sin (latitude) ()
!LLL   REAL    (KIND=JPRB) :: SLON                ! sin (longitude) ()
!LLL   REAL    (KIND=JPRB) :: CLON                ! cos (longitude) ()
!LLL END TYPE GBRAD_PHYS_TYPE  

!LLL TYPE GBRAD_RR_TYPE
!LLL   INTEGER (KIND=JPIM) :: SCREEN  ! screening flag
!LLL   REAL    (KIND=JPRB) :: OBS     ! observed precipitation [log(mm/h +repslog)]
!LLL   REAL    (KIND=JPRB) :: FW      ! modelled precipitation [log(mm/h +repslog)], forward
!LLL   REAL    (KIND=JPRB) :: TL      ! modelled precipitation [log(mm/h +repslog)], perturbation 
!LLL   REAL    (KIND=JPRB) :: AD      ! modelled precipitation [log(mm/h +repslog)], gradient
!LLL   REAL    (KIND=JPRB) :: OERR    ! observation error [log(mm/h +repslog)]
!LLL   REAL    (KIND=JPRB) :: BIAS    ! observation bias  [log(mm/h +repslog)]
!LLL END TYPE GBRAD_RR_TYPE  

! Status code bitfield.  
INTEGER (KIND=JPIM), PARAMETER :: &
 & GBRAD_ACTIVE        = 0 , & ! Observation active if no other bits set
 & GBRAD_OBS_RR        = 1 , & ! Observed TB out of bounds
 & GBRAD_ANAPROP       = 2 , & ! Anomalous propagation in model
 & GBRAD_SEA           = 3 , & ! Sea or lake point in model
 & GBRAD_SNOW          = 4 , & ! Surface snowfall probable according to model 2m temp.
 & GBRAD_NO_OBS        = 5 , & ! No obs for this grid point (grid point space only)
 & GBRAD_FG_DEPAR      = 6 , & ! FG depar outside limit
 & GBRAD_CONTAMINATED  = 7 , & ! Obs rejected because neighbour on low-res grid rejected
 & GBRAD_RR_LOW        = 8 , & ! FG or OBS rain rate below rr_min
 & GBRAD_PASSIVE       = 9     ! Passive observation

CONTAINS
 
! ---------------------------------------------------------------------

FUNCTION GBRAD_OK (KSTATUS) 

! Test status flag to see if observation is OK to process.
! Function works for scalar
     
INTEGER (KIND=JPIM), INTENT (IN) :: KSTATUS
LOGICAL                          :: GBRAD_OK

GBRAD_OK = KSTATUS == IBSET(0,GBRAD_ACTIVE)  &
   &  .OR. KSTATUS == IBSET(IBSET(0,GBRAD_ACTIVE),GBRAD_SNOW)  &
   &  .OR. KSTATUS == IBSET(IBSET(0,GBRAD_ACTIVE),GBRAD_PASSIVE)  &
   &  .OR. KSTATUS == IBSET(IBSET(IBSET(0,GBRAD_ACTIVE),GBRAD_PASSIVE),GBRAD_SNOW) 

END FUNCTION GBRAD_OK 

! ---------------------------------------------------------------------

END MODULE YOMGBRAD
