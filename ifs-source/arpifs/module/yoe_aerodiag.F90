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

MODULE YOE_AERODIAG

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------
!*    ** *YOE_AERODIAG* - FIXED INDICES IN AEROSOL DIAGNOSTICS ARRAYS
!     ------------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: NPAERODIAG = 8    ! per-aerosol-bin diagnostics (flux,od,mass)
INTEGER(KIND=JPIM), PARAMETER :: NPAERAOT = 3      ! per-layer aerosol optical thicknesses
INTEGER(KIND=JPIM), PARAMETER :: NPAERLISI_VAR = 6 ! lidar simulator variables
INTEGER(KIND=JPIM), PARAMETER :: NPAERLISI_WVL = 3 ! lidar simulator wavelengths
INTEGER(KIND=JPIM), PARAMETER :: NPAERLISI = NPAERLISI_VAR * NPAERLISI_WVL

! General per-aerosol-bin diagnostics
! WARNING - The numbering of these may affect the sequence of elements in
!   the PDHFSS array, as filled in postphy.F90
INTEGER(KIND=JPIM), PARAMETER :: JPAERODIAG_SRC = 1 ! Source/gain
INTEGER(KIND=JPIM), PARAMETER :: JPAERODIAG_DDP = 2 ! Dry deposition
INTEGER(KIND=JPIM), PARAMETER :: JPAERODIAG_SDM = 3 ! Sedimentation
INTEGER(KIND=JPIM), PARAMETER :: JPAERODIAG_WDL = 4 ! Wet deposition (large-scale)
INTEGER(KIND=JPIM), PARAMETER :: JPAERODIAG_WDC = 5 ! Wet deposition (convective)
INTEGER(KIND=JPIM), PARAMETER :: JPAERODIAG_NGT = 6 ! Negative fixer
INTEGER(KIND=JPIM), PARAMETER :: JPAERODIAG_OD  = 7 ! Optical depth (at 550nm)
INTEGER(KIND=JPIM), PARAMETER :: JPAERODIAG_MSS = 8 ! Total column mass

! Labels for the above for prefixing field names
CHARACTER(LEN=3), PARAMETER :: CPAERODIAG_LABEL(8) = &
  & (/ 'SRC', 'DDP', 'SDM', 'WDL', 'WDC', 'NGT', 'OD.', 'TC.' /)

! Optical thickness variables
INTEGER(KIND=JPIM), PARAMETER :: JPAERAOT_TOTAL   = 1 ! Total AOT
INTEGER(KIND=JPIM), PARAMETER :: JPAERAOT_NATURAL = 2 ! Natural AOT (SS+DU)
INTEGER(KIND=JPIM), PARAMETER :: JPAERAOT_ANTHRO  = 3 ! Anthropogenic AOT (OM+BC+SU)

! Lidar simulator: variables
INTEGER(KIND=JPIM), PARAMETER :: JPAERLISI_EXT    = 1 ! Extinction coefficient
INTEGER(KIND=JPIM), PARAMETER :: JPAERLISI_SSA    = 2 ! Single-scattering albedo
INTEGER(KIND=JPIM), PARAMETER :: JPAERLISI_BSTOA  = 3 ! Backscatter from TOA (satellite lidar)
INTEGER(KIND=JPIM), PARAMETER :: JPAERLISI_BSGND  = 4 ! Backscatter from ground-based lidar
INTEGER(KIND=JPIM), PARAMETER :: JPAERLISI_MOLBSC = 5 ! Unattenuated molecular backscattering coefficient
INTEGER(KIND=JPIM), PARAMETER :: JPAERLISI_AERBSC = 6 ! Unattenuated aerosol backscattering coefficient
!INTEGER(KIND=JPIM), PARAMETER :: JPAERLISI_EXTD  = 7 ! Dust-only extinction coefficient

! ... and wavelengths
INTEGER(KIND=JPIM), PARAMETER :: JPAERLISI_355  = 1 ! 355 nm
INTEGER(KIND=JPIM), PARAMETER :: JPAERLISI_532  = 2 ! 532 nm
INTEGER(KIND=JPIM), PARAMETER :: JPAERLISI_1064 = 3 ! 1064 nm


! -----------------------------------------------------------------------------
! And the structures for controlling per-wavelength aerosol optical diagnostics

INTEGER(KIND=JPIM), PARAMETER :: NPAERO_WVL_DIAG = 20 ! Maximum number of aerosol optical diagnostic wavelengths
INTEGER(KIND=JPIM), PARAMETER :: NPAERO_WVL_DIAG_TYPES = 5 ! Number of optical diagnostic types

TYPE TYPE_AERO_WVL_DIAG
SEQUENCE
INTEGER(KIND=JPIM) :: IWVL                             ! Wavelength in nm
INTEGER(KIND=JPIM) :: IGRIBDIAG(NPAERO_WVL_DIAG_TYPES) ! GRIB codes for diagnostics as numbered below
END TYPE TYPE_AERO_WVL_DIAG

! Optical depth MUST be number 1, since without interactive aerosol it's the only one used.
INTEGER(KIND=JPIM), PARAMETER :: JPAERO_WVL_AOD       = 1 ! Optical depth
INTEGER(KIND=JPIM), PARAMETER :: JPAERO_WVL_AODABS    = 2 ! Absorption optical depth
INTEGER(KIND=JPIM), PARAMETER :: JPAERO_WVL_AODFM     = 3 ! Fine-mode optical depth
INTEGER(KIND=JPIM), PARAMETER :: JPAERO_WVL_SSA       = 4 ! Single scattering albedo
INTEGER(KIND=JPIM), PARAMETER :: JPAERO_WVL_ASSIMETRY = 5 ! Asymmetry ("assimetry") parameter

! Labels for the above for prefixing field names
CHARACTER(LEN=8), PARAMETER :: CPAERO_WVL_DIAG_LABEL(5) = &
  & (/ 'OptDepth', 'AbOptDep', 'FMOptDep', 'SScatAlb', 'Asymmetr' /)

!------------------------------------------------------------------------------
END MODULE YOE_AERODIAG
