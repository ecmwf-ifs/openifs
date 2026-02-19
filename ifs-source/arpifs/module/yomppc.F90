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

MODULE YOMPPC

USE PARKIND1  ,ONLY : JPIM
USE PARFPOS   , ONLY : JPOSSGP

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------
!*    POST PROCESSING CONTROL 

! LRSUP  : .TRUE. if surface post-processing requested
! LRSACC : .TRUE. if accumulated fluxes are reset at post-processing time
! NO2DGG : number of 2-D gaussian grid fields requested
! M2DGGP : array containing field codes for 2-D gaussian grid fields
! NO3DGGM: number of 3-D gaussian grid fields requested
! M3DGGM : array containing field codes for 3-D gaussian grid fields
! NO2DSPE: number of 2-D ECV spectral fields requested
! M2DSPE : array containing field codes for 2-D ECV spectral fields
! NO3DSPE: number of 3-D ECV spectral fields requested
! M3DSPE : array containing field codes for 3-D ECV spectral fields
! NO3DSPM: number of spherical harmonic fields requested (model level pp.)
! M3DSPM : array containing field codes for 3-D and 2-D fields (model lev. pp.)
! NO2DGGL: number of gaussian grid fields requested (pp. in lagged mode)
! M2DGGPL: array cont. field codes for gaussian grid fields (pp. in lagged mode)
! NRSACCFRQ: if >0 accumulated fluxes are reset at this frequency
! NRSACCOFF: add to time for the NRSACCFRQ>0 option to e.g. allow reset at
!            at 0z for runs starting at 12z

! Backup from fullpos :
!     MFPPHY : Gribcodes of physical fields to processed
!     NFPPHY : useful dimension of CFPPHY

INTEGER(KIND=JPIM), PARAMETER :: JPOPLEV=200
INTEGER(KIND=JPIM), PARAMETER :: JPO2DGG=187
INTEGER(KIND=JPIM), PARAMETER :: JPO3DSPM=35
INTEGER(KIND=JPIM), PARAMETER :: JPO3DGG=38
INTEGER(KIND=JPIM), PARAMETER :: JPO2DSPA=30
INTEGER(KIND=JPIM), PARAMETER :: JPO3DSPA=10
INTEGER(KIND=JPIM), PARAMETER :: JPO3DGGA=10

INTEGER(KIND=JPIM) :: M2DGGP(JPO2DGG)
INTEGER(KIND=JPIM) :: M3DSPM(JPO3DSPM)
INTEGER(KIND=JPIM) :: M2DGGPL(JPO2DGG)
INTEGER(KIND=JPIM) :: M3DGGM(JPO3DGG)
INTEGER(KIND=JPIM) :: M2DSPE(JPO2DSPA)
INTEGER(KIND=JPIM) :: M3DSPE(JPO3DSPA)
INTEGER(KIND=JPIM) :: M3DGGE(JPO3DGGA)

INTEGER(KIND=JPIM) :: NO2DGG
INTEGER(KIND=JPIM) :: NO3DSPM
INTEGER(KIND=JPIM) :: NO2DGGL
INTEGER(KIND=JPIM) :: NO3DGGM
INTEGER(KIND=JPIM) :: NO2DSPE
INTEGER(KIND=JPIM) :: NO3DSPE
INTEGER(KIND=JPIM) :: NO3DGGE

LOGICAL :: LRSUP
LOGICAL :: LRSACC
INTEGER(KIND=JPIM) :: NRSACCFRQ, NRSACCOFF

INTEGER(KIND=JPIM) :: MFPPHY(JPOSSGP)
INTEGER(KIND=JPIM) :: NFPPHY
!     ------------------------------------------------------------------
END MODULE YOMPPC
