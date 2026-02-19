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

MODULE YOMARG

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE PARDIM    ,ONLY : JPMXLE   ,JPMXGL

IMPLICIT NONE

SAVE

!     ----------------------------------------------------------------
!*    Former command line arguments - File frame information
!     ----------------------------------------------------------------

!    * NAMARG variables (these variables are in YOMCT0 too):

!    LELAM      : .TRUE. if LAM model; .FALSE. if global model.
!    NCONF      : Configuration for 3D integration.
!    LECMWF     : .T.: ECMWF configuration
!    CNMEXP     : Name of the experiment
!    CUSTOP     : Forecast range
!    UTSTEP     : Time step in second
!    NSUPERSEDE : =1 date and geometry are read in intial file,
!                 =0 date and geoemtry are read in namelists
!    NFPSERVER  : configuration of post-processing server.
!             = 0 : disabled unless configured by an external script
!                   (filename CSCRIPT_PPSERVER, defined in namelist NAMCT0)
!                   In this case the external script should handle the input,
!                   output, climatology and events files rotation.
!             = 1 : automatic configuration :
!                   input files and output files have the same suffix ("time
!                   stamp") but different prefix.
!                   The list of files is determined by the time stamps
!                   controlled by the model variables LINC, NTIMEFMT for the
!                   format ; and NFRPOS, NPOSTS for the number and values of
!                   iterations.
!                   The initial file (...INIT) is excluded.
!             =-1 : same as NFPSERVER=1, but the initial file (...INIT) is included.
!

LOGICAL:: LELAM=.FALSE.
INTEGER(KIND=JPIM) :: NCONF=-1
LOGICAL:: LECMWF=.TRUE.
CHARACTER (LEN = 4):: CNMEXP='-9'
!CHARACTER (LEN = 8)  :: CUSTOP='-9'
!REAL(KIND=JPRB) :: UTSTEP=-9._JPRB
INTEGER(KIND=JPIM) :: NSUPERSEDE
INTEGER (KIND=JPIM) :: NFPSERVER=0

!    * Complement to NAMARG variables:

!    NUSTOP     : Number of time step (can be given in hours)
!    NECMWF     : Linked with LECMWF

!INTEGER(KIND=JPIM) :: NUSTOP=-9
INTEGER(KIND=JPIM) :: NECMWF=-1

!    * Input file structure (0=FA ; 1=GRIB)

INTEGER(KIND=JPIM) :: NGRIBFILE

!    * Date:

!    NUDATE     : Initial date in the form AAAAMMDD
!    NUSSSS     : Initial time in seconds (e.g. for 12h, 43200)

INTEGER(KIND=JPIM) :: NUDATE
INTEGER(KIND=JPIM) :: NUSSSS

!    * Vertical geometry:

!    NUFLEV     : Number of vertical levels
!    UVALH      : Vertical function *A*
!    UVBH       : Vertical function *B*

INTEGER(KIND=JPIM) :: NUFLEV
REAL(KIND=JPRB) :: UVALH(0:JPMXLE)
REAL(KIND=JPRB) :: UVBH(0:JPMXLE)

!    * Horizontal geometry:

!    NUSMAX     : Nominal truncation
!    NUCMAX     : Upper truncation when using dilatation/contraction or filtering matrices
!    NUDGL      : Number of latitudes
!    NUDLON     : Number of longitudes
!    NUHTYP     : Type of grid (reduced or not)
!    NUSTTYP    : Type of transformation (rotated or not)
!    USTRET     : Stretching factor
!    UMUCEN     : Sine of latitude of the pole of stretching
!    ULOCEN     : Cosine of longitude of the pole of stretching
!    NULOEN     : number of points on each latitude
!    NUMEN      : number of Fourier wavenumbers on each latitude
!    NULIM      : Limited Area characteristics
!    UGEMU      : Horizontal geometry characteristics

INTEGER(KIND=JPIM) :: NUSMAX
INTEGER(KIND=JPIM) :: NUCMAX
INTEGER(KIND=JPIM) :: NUDGL
INTEGER(KIND=JPIM) :: NUDLON
INTEGER(KIND=JPIM) :: NUHTYP
INTEGER(KIND=JPIM) :: NUSTTYP
REAL(KIND=JPRB) :: USTRET
REAL(KIND=JPRB) :: UMUCEN
REAL(KIND=JPRB) :: ULOCEN
INTEGER(KIND=JPIM) :: NULOEN(JPMXGL/2)
INTEGER(KIND=JPIM) :: NUMEN(JPMXGL/2)
INTEGER(KIND=JPIM) :: NULIM(8)
REAL(KIND=JPRB) :: UGEMU(JPMXGL/2)

!     ------------------------------------------------------------------
END MODULE YOMARG
