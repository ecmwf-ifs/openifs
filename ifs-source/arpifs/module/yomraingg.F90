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

MODULE YOMRAINGG
      
USE PARKIND1, ONLY: JPIM, JPRB

IMPLICIT NONE

SAVE

REAL (KIND=JPRB) :: RLIMIT_SD  ! obs-grid search distance limit  [m]
REAL (KIND=JPRB) :: REPSLOG    ! offset in rain rate log transform
REAL (KIND=JPRB) :: RR_MIN     ! threshold minimum RR [mm/h]

LOGICAL :: LERAINGG_ACTIVE

!LLL INTEGER (KIND=JPIM), ALLOCATABLE :: OB_FROM    (:,:)
!LLL INTEGER (KIND=JPIM), ALLOCATABLE :: OB_JOBS    (:,:)
!LLL INTEGER (KIND=JPIM), ALLOCATABLE :: OB_JROF_GP (:,:)
!LLL INTEGER (KIND=JPIM), ALLOCATABLE :: OB_INDEX   (:)
!LLL INTEGER (KIND=JPIM), ALLOCATABLE :: GP_COUNT   (:)
!LLL INTEGER (KIND=JPIM), ALLOCATABLE :: OB_COUNT   (:)

! Accumulation periods definition
INTEGER (KIND=JPIM), PARAMETER   :: NACCLX = 3, NACCTOTX = 33
INTEGER (KIND=JPIM)              :: KACCL(NACCLX), NACCP(NACCLX), KREPTYP(NACCLX)
INTEGER (KIND=JPIM)              :: NACCTOT
INTEGER (KIND=JPIM), ALLOCATABLE :: NACCST(:), NACCEND(:), KACCLT(:)

!LLL TYPE RAINGG_PHYS_TYPE
!LLL   REAL    (KIND=JPRB) :: PRECS   (NACCLX)   ! accumulated surface precipitation (mm/s)
!LLL   REAL    (KIND=JPRB) :: PRECSAD (NACCTOTX) ! accumulated surface precipitation (mm/s) (in adjoint)
!LLL   REAL    (KIND=JPRB) :: T2M                ! 2m temperature (K)
!LLL   REAL    (KIND=JPRB) :: LSM                ! land-sea mask ()
!LLL   REAL    (KIND=JPRB) :: SDO                ! orography standard deviation (m)
!LLL   REAL    (KIND=JPRB) :: SLAT               ! sin (latitude) ()
!LLL   REAL    (KIND=JPRB) :: SLON               ! sin (longitude) ()
!LLL   REAL    (KIND=JPRB) :: CLON               ! cos (longitude) ()
!LLL END TYPE RAINGG_PHYS_TYPE      

!LLL TYPE RAINGG_RR_TYPE
!LLL   INTEGER (KIND=JPIM) :: SCREEN   (NACCLX)   ! screening flag
!LLL   INTEGER (KIND=JPIM) :: SCREENAD (NACCTOTX) ! screening flag (in adjoint)
!LLL   REAL    (KIND=JPRB) :: OBS      (NACCLX)   ! observed precipitation [log(mm/h +repslog)]
!LLL   REAL    (KIND=JPRB) :: FW       (NACCLX)   ! modelled precipitation [log(mm/h +repslog)], forward
!LLL   REAL    (KIND=JPRB) :: TL       (NACCLX)   ! modelled precipitation [log(mm/h +repslog)], perturbation 
!LLL   REAL    (KIND=JPRB) :: AD       (NACCTOTX) ! modelled precipitation [log(mm/h +repslog)], gradient
!LLL   REAL    (KIND=JPRB) :: OERR     (NACCLX)   ! observation error [log(mm/h +repslog)]
!LLL   REAL    (KIND=JPRB) :: BIAS     (NACCLX)   ! observation bias  [log(mm/h +repslog)]
!LLL END TYPE RAINGG_RR_TYPE  

! Status code bitfield.  
INTEGER (KIND=JPIM), PARAMETER :: &
 & RAINGG_ACTIVE        = 0 , & ! Observation active if no other bits set
 & RAINGG_OBS_RR        = 1 , & ! Observed rain amount out of bounds
 & RAINGG_ORO           = 2 , & ! Rugged terrain in model
 & RAINGG_SEA           = 3 , & ! Sea or lake point in model
 & RAINGG_SNOW          = 4 , & ! Surface snowfall probable according to model 2m temp.
 & RAINGG_NO_OBS        = 5 , & ! No obs for this grid point (grid point space only)
 & RAINGG_FG_DEPAR      = 6 , & ! FG depar outside limit
 & RAINGG_CONTAMINATED  = 7 , & ! Obs rejected because neighbour on low-res grid rejected
 & RAINGG_RR_LOW        = 8 , & ! FG or OBS rain rate below rr_min
 & RAINGG_PASSIVE       = 9     ! Passive observation

! Global array for precipitation accumulations in NL 
! and associated perturbations in TL & AD.
REAL (KIND=JPRB), ALLOCATABLE :: RACCPR3 (:,:,:)

CONTAINS
 
! ---------------------------------------------------------------------

FUNCTION RAINGG_OK (KSTATUS) 

! Test status flag to see if observation is OK to process.
! Function works for scalar
     
INTEGER (KIND=JPIM), INTENT (IN) :: KSTATUS
LOGICAL                          :: RAINGG_OK

RAINGG_OK = KSTATUS == IBSET(0,RAINGG_ACTIVE)  &
   &  .OR. KSTATUS == IBSET(IBSET(0,RAINGG_ACTIVE),RAINGG_PASSIVE) 

END FUNCTION RAINGG_OK 

! ---------------------------------------------------------------------

END MODULE YOMRAINGG
