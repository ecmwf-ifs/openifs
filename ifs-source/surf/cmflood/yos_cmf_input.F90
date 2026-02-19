MODULE YOS_CMF_INPUT
!==========================================================
!* PURPOSE: Shared variables for CaMa-Flood model configulation
! (C) D.Yamazaki & E. Dutra  (U-Tokyo/FCUL)  Aug 2019
!
! Licensed under the Apache License, Version 2.0 (the "License");
!   You may not use this file except in compliance with the License.
!   You may obtain a copy of the License at: http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software distributed under the License is 
!  distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
! See the License for the specific language governing permissions and limitations under the License.
!==========================================================
USE PARKIND1,                ONLY: JPIM, JPRB, JPRM
IMPLICIT NONE
SAVE 
!================================================
!*** CMF default files
LOGICAL                         :: LLOGOUT                 !! true: log output to file
INTEGER(KIND=JPIM)              :: LOGNAM                  !! default log    file FID
INTEGER(KIND=JPIM)              :: NSETFILE                !! input namelist file FID
INTEGER(KIND=JPIM)              :: TMPNAM                  !! temporal I/O   file FIG
CHARACTER(LEN=256)              :: CLOGOUT                 !! default log    file name
CHARACTER(LEN=256)              :: CSETFILE                !! input namelist file name

DATA LLOGOUT       /.TRUE./
DATA CLOGOUT       /'./log_CaMa.txt'/
DATA CSETFILE      /'input_cmf.nam'/

!================================================
!*** NAMELIST/NRUNVER/
LOGICAL                         :: LADPSTP                 !! true: use adaptive time step

LOGICAL                         :: LFPLAIN                 !! true: consider floodplain (false: only river channel)
LOGICAL                         :: LKINE                   !! true: use kinematic wave
LOGICAL                         :: LFLDOUT                 !! true: floodplain flow (high-water channel flow) active
LOGICAL                         :: LPTHOUT                 !! true: activate bifurcation scheme
LOGICAL                         :: LDAMOUT                 !! true: activate dam operation (under development)

LOGICAL                         :: LROSPLIT                !! true: input if surface (Qs) and sub-surface (Qsb) runoff
LOGICAL                         :: LWEVAP                  !! true: input water evaporation to extract from floodplain
LOGICAL                         :: LWEVAPFIX               !! true: water balance closure extracting water from evap when available
LOGICAL                         :: LWEXTRACTRIV            !! true: also extract water from rivers 

LOGICAL                         :: LGDWDLY                 !! true: Activate ground water reservoir and delay
LOGICAL                         :: LSLPMIX                 !! true: activate mixed kinematic and local inertia based on slope

LOGICAL                         :: LMEANSL                 !! true : boundary condition for mean sea level
LOGICAL                         :: LSEALEV                 !! true : boundary condition for variable sea level

LOGICAL                         :: LRESTART                !! true: initial condition from restart file
LOGICAL                         :: LSTOONLY                !! true: storage only restart (mainly for data assimilation)
LOGICAL                         :: LOUTPUT                 !! true: use standard output (to file)

LOGICAL                         :: LGRIDMAP                !! true: for standard XY gridded 2D map
LOGICAL                         :: LLEAPYR                 !! true: neglect leap year (Feb29 skipped)
LOGICAL                         :: LMAPEND                 !! true: for map data endian conversion
LOGICAL                         :: LBITSAFE                !! true: for Bit Identical simulation (avoid OSM ATOMIC)
LOGICAL                         :: LSLOPEMOUTH             !! true: prescribe water level slope == elevation slope on river month
                                                           

!================================================
!*** NAMELIST/NCONF/
CHARACTER(LEN=256)              :: CDIMINFO                !! Dimention Information

REAL(KIND=JPRB)                 :: DT                      !! Time Step Length [SEC] (should be multiple of 60)

INTEGER(KIND=JPIM)              :: IFRQ_OUT                !! [hour]: frequency to write output     e.g. (1,2,3,6,12,24) hour
INTEGER(KIND=JPIM)              :: IFRQ_INP                !! [hour]: frequency to update runoff    e.g. (1,2,3,6,12,24) hour
INTEGER(KIND=JPIM)              :: IFRQ_SL                 !! [min]:  frequency to update sea level e.g. (1,2,5,10,15,20,30,60) min

!*** set by CDIMINFO
INTEGER(KIND=JPIM)              :: NX                      !! NUMBER OF GRIDS IN HORIZONTAL
INTEGER(KIND=JPIM)              :: NY                      !! NUMBER OF GRIDS IN VERTICAL
INTEGER(KIND=JPIM)              :: NLFP                    !! NUMBER OF VERTICAL LEVELS DEFINING FLOODPLAIN

INTEGER(KIND=JPIM)              :: NXIN                    !! NUMBER OF GRIDS IN HORIZONTAL
INTEGER(KIND=JPIM)              :: NYIN                    !! NUMBER OF GRIDS IN VERTICAL
INTEGER(KIND=JPIM)              :: INPN                    !! MAX INPUT NUMBER

REAL(KIND=JPRB)                 :: WEST                    !! west, east, north, south edge of the domain [deg]
REAL(KIND=JPRB)                 :: EAST
REAL(KIND=JPRB)                 :: NORTH
REAL(KIND=JPRB)                 :: SOUTH

!*** calculated from IFRQ & DT
REAL(KIND=JPRB)                 :: DTIN                    !! SECOND IN INPUT TIME STEP [SEC]
REAL(KIND=JPRB)                 :: DTSL                    !! SECOND IN TIME STEP [SEC]

!================================================
!*** NAMELIST/PARAM/
!* parameters
REAL(KIND=JPRB)                 :: PMANRIV              !! manning roughness (river)
REAL(KIND=JPRB)                 :: PMANFLD              !! manning roughness (floodplain)
REAL(KIND=JPRB)                 :: PGRV                 !! gravity acceleration [m/s2]
REAL(KIND=JPRB)                 :: PDSTMTH              !! downstream distance at river mouth [m]
REAL(KIND=JPRB)                 :: PCADP                !! CFL coefficient
REAL(KIND=JPRB)                 :: PMINSLP              !! minimum topographic slope (kinematic wave) [m/m]
!* missing values
INTEGER(KIND=JPIM)              :: IMIS                !! integer undefined
REAL(KIND=JPRM)                 :: RMIS                !! real    undefined
REAL(KIND=JPRB)                 :: DMIS                !! double  undefined
!* file suffix
CHARACTER(LEN=256)              :: CSUFBIN            ! .bin suffix for binary (2D map)
CHARACTER(LEN=256)              :: CSUFVEC            ! .vec suffix for binary (1D vector)
CHARACTER(LEN=256)              :: CSUFPTH            ! .pth suffix for binary (1D bifurcation channel)
CHARACTER(LEN=256)              :: CSUFCDF            ! .nc  suffix for netCDF

!================================================
END MODULE YOS_CMF_INPUT
