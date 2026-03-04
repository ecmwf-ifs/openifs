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

MODULE YOMRLX

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------

! LRLXG  : USE OF RELAXATION TERMS

! LRLXVO : relaxation switch for vorticity (altitude)
! LRLXDI : relaxation switch for divergence (altitude)
! LRLXTE : relaxation switch for temperature (altitude)
! LRLXQ : relaxation switch for specific humidity (altitude)
! LRLXQI : relaxation switch for ice water (altitude)
! LRLXQL : relaxation switch for liquid water (altitude)
! LRLXQC : relaxation switch for cloud fraction (altitude)
! LRLXLP : relaxation switch for log(surf pressure)
! LRLXO3 : relaxation switch for ozone (altitude)

! NFRLXG : number of analysis time steps in memory
! NFRLXU : interval between two references fields (in hours)
!          (this is equal to NFRHIS in nudging)
! NPRLX3D: number of 3D parameters being relaxed
! NPRLX2D: number of 2D parameters being relaxed
! NRLXAN : number of reference steps

! XRLXVO : relaxation coefficient for vorticity (altitude)
! XRLXDI : relaxation coefficient for divergence (altitude)
! XRLXTE : relaxation coefficient for temperature (altitude)
! XRLXQ : relaxation coefficient for temperature (altitude)
! XRLXLP : relaxation coefficient for log(surf pressure)
! XRLXO3 : relaxation coefficient for ozone (altitude)

! CRLXPATHGG : path to gridpoint relaxation file location
! CRLXPATHSH : path to spectral relaxation file location

! ALATRLX1 : northern limit
! ALATRLX2 : southern limit
! ALONRLX1 : western limit in [0,360]
! ALONRLX2 : eastern limit in [0,360] (E<W will span Greenwich)

! AXRLX : smoothing parameter in longitude
! AYRLX : smoothing parameter in latitude
! AZRLX : smoothing parameter in the vertical

! NRLXLMIN : uppermost model level
! NRLXLMAX : lowermost model level
! NRLXLMINU: uppermost model level for velocity variables
! NRLXLMAXU: lowermost model level for velocity variables

! NRLXSMAX : maximum spectral resolution that is relaxed

INTEGER(KIND=JPIM) :: NFRLXG, NFRLXU
INTEGER(KIND=JPIM) :: NPRLX3D, NPRLX2D
INTEGER(KIND=JPIM) :: NRLXAN

LOGICAL :: LRLXG

LOGICAL :: LRLXVO
LOGICAL :: LRLXDI
LOGICAL :: LRLXTE
LOGICAL :: LRLXLP
LOGICAL :: LRLXQ
LOGICAL :: LRLXQI
LOGICAL :: LRLXQL
LOGICAL :: LRLXQC
LOGICAL :: LRLXO3

CHARACTER(LEN=512) :: CRLXPATHGG
CHARACTER(LEN=512) :: CRLXPATHSH

REAL(KIND=JPRB) :: XRLXVO
REAL(KIND=JPRB) :: XRLXDI
REAL(KIND=JPRB) :: XRLXTE
REAL(KIND=JPRB) :: XRLXLP
REAL(KIND=JPRB) :: XRLXQ
REAL(KIND=JPRB) :: XRLXO3

REAL(KIND=JPRB) :: ALATRLX1
REAL(KIND=JPRB) :: ALATRLX2
REAL(KIND=JPRB) :: ALONRLX1
REAL(KIND=JPRB) :: ALONRLX2

REAL(KIND=JPRB) :: AXRLX
REAL(KIND=JPRB) :: AYRLX
REAL(KIND=JPRB) :: AZRLX

INTEGER(KIND=JPIM) :: NRLXLMIN
INTEGER(KIND=JPIM) :: NRLXLMAX
INTEGER(KIND=JPIM) :: NRLXLMINU
INTEGER(KIND=JPIM) :: NRLXLMAXU

INTEGER(KIND=JPIM) :: NRLXSMAX 
END MODULE YOMRLX
