! Copyright 2006-2019 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities granted to it by
! virtue of its status as an intergovernmental organisation nor does it submit to any jurisdiction.
!=====================================================================================================

MODULE YOMCMEMPAR

! Purpose :
! -------
!   Universal constants of CMEM
!-----------------------------------------------------------------------------------------------------

USE PARKIND1, ONLY : JPIM, JPRB

IMPLICIT NONE

INTEGER(KIND=JPIM),       PARAMETER :: JPCHARLEN                = 100   ! Length of CMEM option 
INTEGER(KIND=JPIM),       PARAMETER :: JPNAMEIDLEN              = 100   ! Length of CNAMEID 
INTEGER(KIND=JPIM),       PARAMETER :: JPCPOLLEN                = 1     ! Length of polarization character
INTEGER(KIND=JPIM),       PARAMETER :: JPCMEMTILE               = 7     ! number of tiles in CMEM 
INTEGER(KIND=JPIM),       PARAMETER :: CMEM_ID_POL_H            = 1     ! H-Polaraization ID in CMEM
INTEGER(KIND=JPIM),       PARAMETER :: CMEM_ID_POL_V            = 2     ! V-Polaraization ID in CMEM
INTEGER(KIND=JPIM),       PARAMETER :: NLAY_SOIL_LS_DEFAULT     = 3     ! Default number of soil layer in input data 
INTEGER(KIND=JPIM),       PARAMETER :: CL_N                     = 0     ! Check Level: (No check)
INTEGER(KIND=JPIM),       PARAMETER :: CL_W                     = 1     ! Check Level: Warning
INTEGER(KIND=JPIM),       PARAMETER :: CL_R                     = 2     ! Check Level: Warning + Reset parameter
INTEGER(KIND=JPIM),       PARAMETER :: CL_A                     = 3     ! Check Level: Warning + Reset Parameter + Abort
CHARACTER(LEN=JPCHARLEN), PARAMETER :: FNAME_LSM_VERTICAL_RESOL = "LSM_VERTICAL_RESOL.asc" ! Input Filename of soil layers depths
REAL(KIND=JPRB),          PARAMETER :: Z_LSM_DEFAULT(NLAY_SOIL_LS_DEFAULT) = (/0.07_JPRB,    0.28_JPRB,   1.0_JPRB/)
                                                                        ! Default Soil layers depths (m) (HTESSEL)
                                                                        ! (/0.02151_JPRB, 0.1642_JPRB, 1.314_JPRB/) (ORCHIDEE)
 
END MODULE YOMCMEMPAR
