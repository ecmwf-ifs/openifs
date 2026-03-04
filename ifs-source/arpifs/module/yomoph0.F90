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

MODULE YOMOPH0

USE PARKIND1  ,ONLY : JPIM

USE TYPE_FAOPH, ONLY : TFAOPH

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------

! CFNSH : pfn for trajectory spectral data file
! CFNPL : PFN FOR PRESSURE LEVEL SPECTRAL DATA FILE
! CFNVAREPS : pfn for accumulated fields for VAREPS
! CFNGG : pfn for trajectory grid point data file
! CFNUA : pfn for trajectory upper air grid point data file
! CFNRF : pfn for background=first guess/reference spectral data file
! CFNRFBDP : pfn for boundary pert in grid-point space
! CFNAN : pfn for analysis field during simulation spectral data file
! CFNGR : pfn for gradient field during simulation spectral data file
! CFNINSH : pfn for spectral part of increment field during simulation
! CFNINGG : pfn for grid point part of increment field during simulation
! CFNISH : pfn for initial spectral data file
! CFNIGG : pfn for initial surface grid point data file
! CFNIUA : pfn for initial upper air grid point data file
! CFNGSH : pfn for initial point of the minimization spectral data file
! CFNGSHBDP : pfn for boundary pert in spectral space
! CFNGGG : pfn for initial point of the minimization grid point data file
! CFNDDH : pfn for files produced by the DDH diagnostics.
! CFNPANSH : pfn for previous analysis spectral data file
! CFNPANGG : pfn for previous analysis grid-point data file
! CFNHWF : pfn for history witness file
! CFNBS  : pfn for bias spectral data file (dig. filtered guess - guess)
! CFANS  : pfn for files containing only surface fields assimilation.
! CFNBGV : pfn for random background vector file
! CFNCGL : pfn for (scaled) eigenvector of the Hessian (CGL algorithm)
! CFNFGI : pfn for background low resolution (incremental) SH
! CFNANI : pfn for analysis low resolution (incremental) SH
! CFNFGIG: pfn for background low resolution (incremental) GP
! CFNANIG: pfn for analysis low resolution (incremental) GP
! CFNTRAJHRGRID: pfn for files containing the high resolution trajectory:
!                GRIB grid-point fields.
! CFNTRAJHRSPEC: pfn for files containing the high resolution trajectory:
!                GRIB upper air spectral fields.
! CFNTRAJHRSURF: pfn for files containing the high resolution trajectory:
!                GRIB surface fields.
! CFNTRAJHR    : pfn for files containing the high resolution trajectory:
!                ARPEGE files.
! CFNTRAJBGGRID: pfn for files containing the background:
!                GRIB grid-point fields.
! CFNTRAJBGSPEC: pfn for files containing the background:
!                GRIB upper air spectral fields.
! CFNTRAJBG    : pfn for files containing the background:
!                ARPEGE files.
! CEFLS   : pfn for lateral boundary compacted files (special internal file)
! CEFNLSH : pfn for lateral boundary coupling files (ALADIN file)
! CETSTAMP : Time stamp for lateral boundary coupling files (ALADIN file)
! CFNBGHRSH    : pfn for the high resolution background at all times
!                where the increments are added: spectral part.
! CFNBGHRGG    : pfn for the high resolution background at all times
!                where the increments are added: gridpoint part.
! CFNCLIMIN  : absolute pfn for input climatology file(s)
! CFNCLIMOUT : absolute pfn for output climatology file(s)
! CNMCA    : name of standard frame of *FA*
! LINC     : incremental switch for files (true : hour, false : time step)
! CCLIMINC : suffix in front of month value for climatology files (no month value if empty string)
! CFPEXTSFX : extension used after the time stamp for the the Surfex files
! NTIMEFMT : time format in output files; 0: +HHHH, 1: +HHHH.mm
! LBCINC   : incremental switch for input LBC files (true : hour, false : number)
! CFPATH   : directory name (ending with a '/') of time-dependent output files
! NCADFORM:format of *FA* frames in new limited area files
!          0=old format, 1=new format
! LBCRESTART: True when reading LBC files on restart

CHARACTER(LEN=210) :: CFNSH
CHARACTER(LEN=210) :: CFNGG
CHARACTER(LEN=210) :: CFNUA
CHARACTER(LEN=16) :: CFNISH
CHARACTER(LEN=16) :: CFNIGG
CHARACTER(LEN=16) :: CFNIUA
CHARACTER(LEN=16) :: CFNPL
CHARACTER(LEN=16) :: CFNVAREPS
CHARACTER(LEN=16) :: CNMCA
CHARACTER(LEN=16) :: CFNRF
CHARACTER(LEN=16) :: CFNRFBDP
CHARACTER(LEN=16) :: CFANS
CHARACTER(LEN=16) :: CFNAN
CHARACTER(LEN=16) :: CFNGR
CHARACTER(LEN=16) :: CFNINSH
CHARACTER(LEN=16) :: CFNINGG
CHARACTER(LEN=256):: CFNHWF
CHARACTER(LEN=16) :: CFNGSH
CHARACTER(LEN=16) :: CFNGSHBDP
CHARACTER(LEN=16) :: CFNGGG
CHARACTER(LEN=16) :: CFNDDH
CHARACTER(LEN=16) :: CFNPANSH
CHARACTER(LEN=16) :: CFNPANGG
CHARACTER(LEN=200):: CFPATH
CHARACTER(LEN=16) :: CFNBS
CHARACTER(LEN=16) :: CFNBGV
CHARACTER(LEN=16) :: CFNCGL
CHARACTER(LEN=16) :: CFNFGI
CHARACTER(LEN=16) :: CFNANI
CHARACTER(LEN=16) :: CFNFGIG
CHARACTER(LEN=16) :: CFNANIG
CHARACTER(LEN=17) :: CFNTRAJHRGRID
CHARACTER(LEN=17) :: CFNTRAJHRSPEC
CHARACTER(LEN=20) :: CFNTRAJHRSURF
CHARACTER(LEN=16) :: CFNTRAJHR
CHARACTER(LEN=17) :: CFNTRAJBGGRID
CHARACTER(LEN=17) :: CFNTRAJBGSPEC
CHARACTER(LEN=16) :: CFNTRAJBG
CHARACTER(LEN=16) :: CEFLS
CHARACTER(LEN=200) :: CEFNLSH
CHARACTER(LEN=5) :: CETSTAMP
CHARACTER(LEN=16) :: CFNBGHRSH
CHARACTER(LEN=16) :: CFNBGHRGG
CHARACTER(LEN=256) :: CFNCLIMIN
CHARACTER(LEN=256) :: CFNCLIMOUT
CHARACTER(LEN=8) :: CCLIMINC
CHARACTER(LEN=8)   :: CFPEXTSFX

LOGICAL :: LBCRESTART = .FALSE.

LOGICAL :: LINC
INTEGER(KIND=JPIM) :: NTIMEFMT = 0


LOGICAL :: LBCINC

INTEGER(KIND=JPIM) :: NCADFORM

TYPE(TFAOPH) :: YMDLOPH(1)
!     ------------------------------------------------------------------
END MODULE YOMOPH0
