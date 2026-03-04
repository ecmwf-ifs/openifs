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

MODULE YOMLIMB

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!*

!     /YOMLIMB/  N. Bormann  05/11/2004 

!     MODIFICATIONS:


!     VARIABLES DEFINING LOGISTICS FOR LIMB COMPUTATIONS.

!  Define sensor numbering
INTEGER(KIND=JPIM), PARAMETER :: MXSENSOR_LRAD   = 2  ! Maximum number of rtlimb sensors
INTEGER(KIND=JPIM), PARAMETER :: MSENSOR_MIPAS   = 32 ! RTLIMB MIPAS ID
INTEGER(KIND=JPIM), PARAMETER :: MSENSOR_MLS     = 33 ! RTLIMB MLS ID
INTEGER(KIND=JPIM), PARAMETER :: MSENSOR_MIPASSQ = 1  ! MIPAS sensor sequence number
INTEGER(KIND=JPIM), PARAMETER :: MXCHAN_LRAD      = 500 ! Maximum number of channels (used for ROERR_LRAD & ROBIAS_LRAD only)

! List linking BUFR-SATID, RT-PLATFORM/SERIES-NUMBER, RT-SATELLITE-NUMBER
INTEGER(KIND=JPIM), PARAMETER :: NRTLIMB_LIST=2
INTEGER(KIND=JPIM)            :: IRTLIMB_LIST(NRTLIMB_LIST, 3) 

! Thesholds used in cloud and top screening for MIPAS:
! PP_MIPAS_ZCLOUD_THRESH - Use cloud detect. for tan hgts below this level [m]
! PP_MIPAS_CI_THRESH     - Cloud index threshold (Kramer et al. at Envis.Symp.)
! PP_MIPAS_WR_THRESH   - Threshold on clearest MIPAS channel (960.7 cm-1)[r.u.]
! PP_MIPAS_ZTAN_MAX      - Maximum tangent height allowed [m]
! PP_MIPAS_ZTAN_MIN      - Minimum tangent height allowed [m]
REAL(KIND=JPRB), PARAMETER :: PP_MIPAS_ZCLOUD_THRESH = 39000._JPRB
REAL(KIND=JPRB), PARAMETER :: PP_MIPAS_CI_THRESH = 5._JPRB!
REAL(KIND=JPRB), PARAMETER :: PP_MIPAS_WR_THRESH = 100._JPRB ! (Dudhia: 125.)
REAL(KIND=JPRB), PARAMETER :: PP_MIPAS_ZTAN_MAX  = 55000._JPRB
REAL(KIND=JPRB), PARAMETER :: PP_MIPAS_ZTAN_MIN  = 12000._JPRB

! Parameters used in RTLIMB preparations
! PP_TOLERANCE - Allow use of RTLIMB coefficient outside limits by this fract.
! PP_MESO_P    - Pressure limit [Pa] above which mesosperic extrapolation of
!                temperature above model top is applied
! PP_MESO_DT_BY_DLNP - Fixed mesospheric temperature lapse rate used for 
!                extrapolation above model top
REAL(KIND=JPRB), PARAMETER :: PP_TOLERANCE       = 0.1_JPRB
REAL(KIND=JPRB), PARAMETER :: PP_MESO_P          = 60._JPRB
REAL(KIND=JPRB), PARAMETER :: PP_MESO_DT_BY_DLNP = 15._JPRB

!     ROERR_LRAD   R    OBSERVATION ERRORS FOR LIMB RADIANCES
!     ROBIAS_LRAD  R    OBSERVATION BIASES FOR LIMB RADIANCES
REAL(KIND=JPRB) :: ROERR_LRAD(0:MXSENSOR_LRAD,MXCHAN_LRAD)
REAL(KIND=JPRB) :: ROBIAS_LRAD(0:MXSENSOR_LRAD,MXCHAN_LRAD)

! Some dimensions. These are set via sulimb or rtl_setup, defined through
! the used data or through values defined in the RTLIMB code.

INTEGER(KIND=JPIM) :: NLIMBGRP ! Size of Y_LIMBGRP_TABLE
INTEGER(KIND=JPIM) :: NJPPF_RTL   ! Maximum number of profs that can be processed in one go
INTEGER(KIND=JPIM) :: NJPNSAT_RTL ! Maximum number of coefficient files
INTEGER(KIND=JPIM) :: NPNAV_RTL   ! Number of profile variables
INTEGER(KIND=JPIM) :: NPMXTAN  ! Maximum number of tangent heights

LOGICAL            :: LQTAN_PRES ! Use tangent pressures?


! Define table for administration of sat-ID, sensor and codetype...
TYPE Y_LIMBGRP_T
 INTEGER(KIND=JPIM)         :: BUFRSATID          ! BUFR Satellite Identifier
 INTEGER(KIND=JPIM)         :: CODETYPE           ! Obs codetype (as in ODB)
 INTEGER(KIND=JPIM)         :: SENSOR             ! Sensor used by RTLIMB
 INTEGER(KIND=JPIM)         :: SATGROUP           ! Satellite data group number
 INTEGER(KIND=JPIM)         :: RTID               ! SatID used by RTLIMB
 INTEGER(KIND=JPIM)         :: RTSERIES           ! Series/platform
 INTEGER(KIND=JPIM)         :: RTCOEF_POS         ! Pointer to RT-coefficients
 INTEGER(KIND=JPIM)         :: NLSAT              ! Number of levels used in RT model
 REAL(KIND=JPRB), POINTER   :: RADPRE(:)=>NULL()  ! Pressure levels used in RT model
 REAL(KIND=JPRB), POINTER   :: TMINRT(:)=>NULL()  ! Limit profiles for 
 REAL(KIND=JPRB), POINTER   :: TMAXRT(:)=>NULL()  ! validity of RT coefficients
 REAL(KIND=JPRB), POINTER   :: QMINRT(:)=>NULL()
 REAL(KIND=JPRB), POINTER   :: QMAXRT(:)=>NULL()
 REAL(KIND=JPRB), POINTER   :: O3MINRT(:)=>NULL()
 REAL(KIND=JPRB), POINTER   :: O3MAXRT(:)=>NULL()
 INTEGER(KIND=JPIM)         :: NMXCHAN            ! Number of channels in coef structure
 INTEGER(KIND=JPIM), POINTER:: IVALCHAN(:)=>NULL()! List of valid channel numbers
END TYPE Y_LIMBGRP_T

TYPE(Y_LIMBGRP_T), ALLOCATABLE, TARGET :: Y_LIMBGRP_TABLE(:)  ! Limb data group admin table

TYPE Y_LIMB_CHAN_SEL_T
 ! channel selection structure for limb sounding
 LOGICAL                  :: ACTIVE                 ! Use channel selection?
 CHARACTER(LEN=256)       :: FILE                   ! Name of channel selection file
 INTEGER(KIND=JPIM)       :: NCHAN                  ! Number of selected channels
 INTEGER(KIND=JPIM)       :: NHEIGHT                ! Number of heights used for mask
 REAL(KIND=JPRB), POINTER :: WNO(:)=>NULL()         ! Selected wavenumbers 
 REAL(KIND=JPRB), POINTER :: MSK_HEIGHTS(:)=>NULL() ! Height borders for mask
 LOGICAL, POINTER         :: MSK(:,:)=>NULL()       ! Mask
 REAL(KIND=JPRB), POINTER :: HEIGHT_LOW(:)=>NULL()  ! Lowest height
 REAL(KIND=JPRB), POINTER :: HEIGHT_UP(:)=>NULL()   ! Heighest height     
END TYPE Y_LIMB_CHAN_SEL_T

TYPE(Y_LIMB_CHAN_SEL_T), ALLOCATABLE, TARGET  :: Y_LIMB_CHAN_SEL(:) ! Channel selection 


! for setting up 2d plane


REAL(KIND=JPRB), PARAMETER :: DTHETA = 6.27844E-3_JPRB ! T511 in Radians
REAL(KIND=JPRB), PARAMETER :: DN_DX_MAX = 0.157_JPRB ! setting the maxiumin N gradients

! 2D bending angle computation SBH Sept 19, 2013

REAL(KIND=JPRB), PARAMETER :: ZED_2D_TRAJ = 5.0E5_JPRB
REAL(KIND=JPRB), PARAMETER :: ZED_2D_MIN  = 5.0E5_JPRB

INTEGER(KIND=JPIM), PARAMETER :: INUM_TRAJ   = 31
INTEGER(KIND=JPIM), PARAMETER :: ISPLIT_TRAJ = 2
INTEGER(KIND=JPIM), PARAMETER :: ISPLIT_MIN  = 2
INTEGER(KIND=JPIM), PARAMETER :: NGPS_GROUP  = 11
INTEGER(KIND=JPIM), PARAMETER :: ISPLINE=4

REAL(KIND=JPRB), PARAMETER :: RINT_HEIGHT = 1.0E4_JPRB

! radius of curvature values sbh 5-3-2013

REAL(KIND=JPRB), PARAMETER :: ROC_MAX = 6.6E6_JPRB
REAL(KIND=JPRB), PARAMETER :: ROC_MIN = 6.2E6_JPRB
REAL(KIND=JPRB), PARAMETER :: RUNDUL_MAX = 1.5E2_JPRB

! SBH 7-4-2006
! DR  10-7-2018
! GPSRO SAT IDs 

INTEGER(KIND=JPIM), PARAMETER :: MMETOP_1      = 3 
INTEGER(KIND=JPIM), PARAMETER :: MMETOP_2      = 4
INTEGER(KIND=JPIM), PARAMETER :: MMETOP_3      = 5
INTEGER(KIND=JPIM), PARAMETER :: MCHAMP        = 41
INTEGER(KIND=JPIM), PARAMETER :: MTERRASARX    = 42
INTEGER(KIND=JPIM), PARAMETER :: MTANDEMX      = 43
INTEGER(KIND=JPIM), PARAMETER :: MPAZ          = 44
INTEGER(KIND=JPIM), PARAMETER :: MSENT6A       = 66
INTEGER(KIND=JPIM), PARAMETER :: MSPIRE        = 269

INTEGER(KIND=JPIM), PARAMETER :: MOCEANSAT_2   = 421

INTEGER(KIND=JPIM), PARAMETER :: MFY_3C        = 522
INTEGER(KIND=JPIM), PARAMETER :: MFY_3D        = 523
INTEGER(KIND=JPIM), PARAMETER :: MFY_3E        = 524

INTEGER(KIND=JPIM), PARAMETER :: MGRACE_A      = 722
INTEGER(KIND=JPIM), PARAMETER :: MGRACE_B      = 723
INTEGER(KIND=JPIM), PARAMETER :: MCOSMIC_1     = 740
INTEGER(KIND=JPIM), PARAMETER :: MCOSMIC_2     = 741
INTEGER(KIND=JPIM), PARAMETER :: MCOSMIC_3     = 742
INTEGER(KIND=JPIM), PARAMETER :: MCOSMIC_4     = 743
INTEGER(KIND=JPIM), PARAMETER :: MCOSMIC_5     = 744
INTEGER(KIND=JPIM), PARAMETER :: MCOSMIC_6     = 745
INTEGER(KIND=JPIM), PARAMETER :: MCOSMIC2_E1   = 750
INTEGER(KIND=JPIM), PARAMETER :: MCOSMIC2_E2   = 751
INTEGER(KIND=JPIM), PARAMETER :: MCOSMIC2_E3   = 752
INTEGER(KIND=JPIM), PARAMETER :: MCOSMIC2_E4   = 753
INTEGER(KIND=JPIM), PARAMETER :: MCOSMIC2_E5   = 754
INTEGER(KIND=JPIM), PARAMETER :: MCOSMIC2_E6   = 755

INTEGER(KIND=JPIM), PARAMETER :: MCNOFS        = 786
INTEGER(KIND=JPIM), PARAMETER :: MGRACE_C      = 803
INTEGER(KIND=JPIM), PARAMETER :: MGRACE_D      = 804

INTEGER(KIND=JPIM), PARAMETER :: MCOMS_1       = 810
INTEGER(KIND=JPIM), PARAMETER :: MCOMS_2       = 811

INTEGER(KIND=JPIM), PARAMETER :: MSACC         = 820
INTEGER(KIND=JPIM), PARAMETER :: MKOMPSAT5     = 825

!     ------------------------------------------------------------------

END MODULE YOMLIMB
