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

MODULE YOMJBECV 

!     Purpose.
!     --------
!       Data and controls for extended control variable

!     Author.
!     -------
!       S. Massart

!     Modifications.
!     --------------
!       Original    April-2018
! ------------------------------------------------------------------


USE PARKIND1, ONLY: JPIM, JPRB

IMPLICIT NONE

SAVE

TYPE TECV_CONFIG
  !   ------------------------------------------------------------------
  !   NINTERPECV          Interpolation type
  !   CFNECVINC           Filename of the gridpoint ECV increment
  !   CFNECVHR        Filename of the gridpoint ECV trajectory
  !   ------------------------------------------------------------------
  !
  CHARACTER(LEN=1)   :: CFNECVINC = 'E'
  CHARACTER(LEN=5)   :: CFNECVFG = 'ICEGG'
  CHARACTER(LEN=16)  :: CFNECVHR = 'TRAJHR00/ecvgrid'
  CHARACTER(LEN=11)  :: CFNECV4V = 'ecv4v_ggml_'
END TYPE TECV_CONFIG

TYPE TECV_DATA
  !   ------------------------------------------------------------------
  !   CSETDESC          Description of each variable
  !   CECV_STD_STRINGS  Description for the errgrib files
  !   CECV_COR_STRINGS  Description of each variable in the wavelet file
  !   NECV              Total number of ECV fields
  !   NECV_2D           Total number of 2D ECV 
  !   NECV_3D           Total number of 3D ECV 
  !   NVIDS_ECV         Number of ECV field containers
  !   NVIDS_ECV_TOTFLDS Total Dimension of the ECV fields
  !   MGRIBS_ECV_FIDS   Grib numbers of the ECV fields
  !   NFLDS_VIDS        Dimesion of each ECV field containers
  !   VIDS_ECV          Id of the each ECV field containers
  !   LECV_1D           True if in 1D 
  !   LECV_2D           True if 2D
  !   LECV_3D           True if 3D
  !
  !   ------------------------------------------------------------------
  CHARACTER(LEN=20), ALLOCATABLE :: CSETDESC(:)
  CHARACTER(LEN=40), ALLOCATABLE :: CECV_STD_STRINGS(:)
  CHARACTER(LEN=40), ALLOCATABLE :: CECV_COR_STRINGS(:)
  INTEGER(KIND=JPIM)             :: NECV
  INTEGER(KIND=JPIM)             :: NECV_1D
  INTEGER(KIND=JPIM)             :: NECV_2D
  INTEGER(KIND=JPIM)             :: NECV_3D
  INTEGER(KIND=JPIM)             :: NVIDS_ECV=0
  INTEGER(KIND=JPIM)             :: NVIDS_ECV_TOTFLDS=0
  INTEGER(KIND=JPIM),ALLOCATABLE :: MGRIBS_ECV_FIDS(:)
  INTEGER(KIND=JPIM),ALLOCATABLE :: NFLDS_VIDS(:)
  INTEGER(KIND=JPIM),ALLOCATABLE :: VIDS_ECV(:)
  LOGICAL, ALLOCATABLE           :: LECV_1D(:)
  LOGICAL, ALLOCATABLE           :: LECV_2D(:)
  LOGICAL, ALLOCATABLE           :: LECV_3D(:)
END TYPE TECV_DATA

!     ------------------------------------------------------------------

!   LECPHYSPARECV       True if optimisation of parameters
!   LINVERACV           True if emission inversion (CAMS)
!   LJB_ALPHA_CV        True if hybrid B using Alpha Control Variable 
!   LSKTECV             True if skin temperature in ECV
!   LSSHECV             True if sea surface height in ECV
!   LTSLECV             True if surface temperature in ECV
!   NDIAECV             Level of output diagnostics
!
!   YRECVDATA           ECV data 
!   YRECVCONFIG         ECV configuration from namelist


!     ------------------------------------------------------------------

LOGICAL               :: LECPHYSPARECV
LOGICAL               :: LINVERACV
LOGICAL               :: LJB_ALPHA_CV
LOGICAL               :: LSKTECV
LOGICAL               :: LSSHECV
LOGICAL               :: LTSLECV
INTEGER(KIND=JPIM)    :: NDIAECV
INTEGER(KIND=JPIM)    :: NINTERPTRAJ_ECV
INTEGER(KIND=JPIM)    :: NINTERPINCR_ECV
REAL(KIND=JPRB)       :: WINDOW_LENGTH=-1.0

TYPE(TECV_DATA)       :: YRECVDATA
TYPE(TECV_CONFIG)     :: YRECVCONFIG

CHARACTER(LEN =20), PARAMETER  :: CKNOWNECV(5) = (/  &
  &   'SOLAR_CONSTANT      ', &
  &   'SKTECV              ', &
  &   'JB_HYBRID_ALPHACV   ', &
  &   'FLUX_INVER          ', &
  &   'EC_PHYS             '/)

!-----------------------------------------------------------------------

END MODULE YOMJBECV
