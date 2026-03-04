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

MODULE YOMSMOS
      
USE PARKIND1,   ONLY : JPIM, JPRB
USE PARDIM,     ONLY : JPMXGL
!use parsmos,  only: npol_max, nang_max, nsmos_max
USE YOMCMEMPAR, ONLY : JPCHARLEN

IMPLICIT NONE

SAVE

INTEGER(KIND=JPIM)   :: NLENMAX      ! max number of obs for allocation 

LOGICAL :: LESMOS_ACTIVE, LESMAP_ACTIVE, LESMOS_SEKF

! Variables needed to assimilate SMOS data in the SEKF
INTEGER (KIND=JPIM) :: NPOL_MAX              ! number of pols monitored
INTEGER (KIND=JPIM) :: NANG_MAX              ! number of angles monitored
INTEGER (KIND=JPIM) :: NSMOS_MAX             ! maximum number of assimilated observations per gp
REAL    (KIND=JPRB), ALLOCATABLE :: NSMOS_ANG(:)  ! which angles monitored
REAL    (KIND=JPRB) :: DELTA_ANG             ! offset for incidence angle


TYPE SMOS_PHYS_TYPE
 INTEGER (KIND=JPIM) :: NLEVELSUR          ! number of surface levels
! prognostic/diagnostic meteo variables 
  REAL    (KIND=JPRB), POINTER :: SM (:)=>NULL()  ! soil moisture 
  REAL    (KIND=JPRB), POINTER :: ST (:)=>NULL()  ! soil temperature
  REAL    (KIND=JPRB)   :: SKT             ! skin temperature
  REAL    (KIND=JPRB)   :: T2M             ! 2m temperature
  REAL    (KIND=JPRB)   :: TLK             ! mixed-layer lake temperature
  REAL    (KIND=JPRB)   :: SNDT            ! snow density
  REAL    (KIND=JPRB)   :: SNDP            ! snow depth
  REAL    (KIND=JPRB)   :: CI              ! sea-ice fraction
! vegetation
  REAL    (KIND=JPRB)   :: LAIH            ! LAI high vegetation
  REAL    (KIND=JPRB)   :: LAIL            ! LAI low vegetation
  REAL    (KIND=JPRB)   :: FCH             ! fraction cover high vegetation
  REAL    (KIND=JPRB)   :: FCL             ! fraction cover low vegetation
 INTEGER  (KIND=JPIM)   :: VEGTH           ! high vegetation type
 INTEGER  (KIND=JPIM)   :: VEGTL           ! low vegetation type
! static variables
  REAL    (KIND=JPRB)   :: OR              ! geopotential at surface (orography)
  REAL    (KIND=JPRB)   :: LSM             ! land-sea mask
 INTEGER  (KIND=JPIM)   :: STP             ! soil type
END TYPE SMOS_PHYS_TYPE  

TYPE SMOS_TB_TYPE
  INTEGER (KIND=JPIM), POINTER :: SCREEN (:,:)=>NULL() ! screening flag 
  REAL    (KIND=JPRB), POINTER :: OBS    (:,:)=>NULL() ! observed brightness temperature
  REAL    (KIND=JPRB), POINTER :: FW     (:,:)=>NULL() ! modelled brightness temperature
  REAL    (KIND=JPRB), POINTER :: ANG    (:,:)=>NULL() ! incidence angle 
  REAL    (KIND=JPRB), POINTER :: FAR    (:,:)=>NULL() ! faraday rotational angle
  REAL    (KIND=JPRB), POINTER :: GEO    (:,:)=>NULL() ! geometrical rotational angle
  REAL    (KIND=JPRB), POINTER :: OERR        =>NULL()      ! observation error 
  REAL    (KIND=JPRB), POINTER :: BIAS        =>NULL()      ! observation bias  
END TYPE SMOS_TB_TYPE

TYPE SMOS_CMEM_TYPE
  CHARACTER(LEN=JPCHARLEN) :: CIDIEL
  CHARACTER(LEN=JPCHARLEN) :: CITEFF
  CHARACTER(LEN=JPCHARLEN) :: CISMR
  CHARACTER(LEN=JPCHARLEN) :: CIRGHR
  CHARACTER(LEN=JPCHARLEN) :: CIVEG
  CHARACTER(LEN=JPCHARLEN) :: CIATM
  CHARACTER(LEN=JPCHARLEN) :: CITVEG
  CHARACTER(LEN=JPCHARLEN) :: CIDVEG 
  CHARACTER(LEN=JPCHARLEN) :: CITDIEL 
END TYPE SMOS_CMEM_TYPE 

! Status code bitfield.  
INTEGER (KIND=JPIM), PARAMETER :: &
 & SMOS_USE            = 0 , & ! Observation active if no other bits set
 & smos_obs_tb         = 1 , & ! Observed TB out of bounds
 & smos_anaprop        = 2 , & ! Anomalous propagation in model
 & smos_sea            = 3 , & ! Sea or lake point in model
 & smos_snow           = 4 , & ! Surface snowfall probable according to model 2m temp.
 & smos_no_obs         = 5 , & ! No obs for this grid point (grid point space only)
 & smos_fg_depar       = 6 , & ! FG depar outside limit
 & smos_obs_pos        = 7 , & ! Position of obs in ILEN to active model
 & smos_coast          = 8     ! Coastal point defined by 0.01 < lsm < 0.95 

! Value of status code bitfield when observation is active
INTEGER (KIND=JPIM), PARAMETER :: SMOS_ACTIVE = IBSET(0,SMOS_USE)
INTEGER (KIND=JPIM), PARAMETER :: SMOS_ACTIVE_SEA = IBSET(0,SMOS_SEA)

!Variables for SMOS CDF matching 
!     ASMOS, BSMOS - parameters for linear scaling of TB 
!            derived from an offline analysis (stored in grib fields at T511 res.)
!     INROWS, ILONE, ZLAT0, ZLATI0 - define the geometry of the ASMOS, BSMOS grib fields
INTEGER(KIND=JPIM)             :: INROWS
INTEGER(KIND=JPIM)             :: ILONE(JPMXGL)
REAL  (KIND=JPRB)              :: ZLAT0, ZLATI0
REAL  (KIND=JPRB), ALLOCATABLE :: ASMOS(:,:,:), BSMOS(:,:,:)

END MODULE YOMSMOS
