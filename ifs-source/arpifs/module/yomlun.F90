! (C) Copyright 1988- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE YOMLUN

USE PARKIND1  ,ONLY : JPIM
USE YOMLUN_IFSAUX, ONLY : NULOUT, NULERR

IMPLICIT NONE

SAVE

PRIVATE
PUBLIC :: NULERR, NULSTAT, NULNAM, NSCRTCH, NULFP10, NULCL1, NULFP11, &
        & NULCL2, NULFP09, NULFP08, NULFP07, NULASE, NULASS, &
        & NULFP06, NULFP05, NULOUT, NULRAD, NULRTL, NULFP04, &
        & NULFP12, NULFP13, NULFP14, NULFP15, NSCATAB, NSCASIG, &
        & NSCASPE, NEGASH, NULFP03, NULDILA, NULCONT, NULFP02, NPOSSH, &
        & NPODDH, NULFP01, NULCO, NTIDE, NTRJSH, NULTRAJHR, &
        & NEFLSS, NULFPOS, NEFLS, NINMSH, NINISH, &
        & NINIGG, NFGISH, NULTRAJBG, NBIAS, NPPPSH, NPDIRL, NULHWF, &
        & NULRCF, NULUSR1, NULUSR2, NULUSR3, NULUSR4, NULUSR5, &
        & NULBDESH, NULBDEGP, &
        & NPRT0,  NPRT1,  NPRT2,  NPRT3,  NPRT4,&
        & NUIO_SERV_LOG, RESERVE_LUN, FREE_LUN, FOPEN

!     ------------------------------------------------------------------

!*    Logical units used by code

! NB: Please use the RESERVE_LUN/FREE_LUN mechanism where possible.

!     NULNAM :   unit number for namelist
!     NULCL1 :   unit number for climatological fields (month before)
!     NULCL2 :   unit number for climatological fields (month after)
!     NTRJSH :   unit number for trajectory spectral data          WRTRA
!     NINMSH :   unit number for initial point of the minimization SUVAZX
!     NINISH :   unit number for initial spectral data             SUSPEC
!     NINIGG :   unit number for initial grid-point data           SUSPEC
!     NFGISH :   unit number for first-guess spectral data
!     NPOSSH :   output unit number (spectral fields)              CPREP1
!     NTIDE  :   unit number for the LFI file containing the total tendencies
!     NPDIRL :   unit number for post-processing directory listing
!     NPPPSH :   unit number for post-processed spherical harmonics WRPLPP

!     NULDILA:   unit number for dilatation matrix (SUDIL,DILAT,SPDILA)
!     NULCONT:   unit number for contraction matrix (SUDIL,DILAT,SPDILA)

!     NULCO  :   unit number for coupled fields (ICMCO)
!     NPODDH :   unit number for mask diagnostic files (DDH)
!     NULRCF :   unit number for restart control file
!     NULHWF :   unit number for history witness file
!     NBIAS  :   unit number for bias (dig. filt. guess - guess)

!     NEFLS  :   unit number for coupling ALADIN file
!     NEFLSS :   unit number for coupling ALADIN file (initialisation)

!     NULUSR1:   unit numbers for user defined files
!     NULSTAT:   unit number for status  file

!     NULASE :   unit number for CANARI statistics (forecast error s.d.)
!     NULASS :   unit number for CANARI statistics (analysis error s.d.)

!     NULUSR2
!     NULUSR3
!     NULUSR4
!     NULUSR5

!     NULFPxx    unit numbers for Full-POS output files
!     NSCRTCH:   unit number for Full-POS scratch file (for in-line post-proc.)
!     NULFPOS:   unit number for Full-POS control file (end of post-processing
!                in conf. 001 ; auxilary namelist file in conf. 927)
!     NULRAD :   unit number for writing radiation diagnostics
!     NULRTL :   unit number for reading RTLIMB coefficient files
!     NSCATAB    SCAT. SIGMA0 TABLE
!     NSCASIG    SCAT. SIGMA0 BIAS CORRECTION
!     NSCASPE    SCAT. SPEED BIAS CORRECTION
!     NEGASH     UNIT NUMBER FOR JK INPUT

!     NULTRAJHR: unit number for high resolution trajectory (option LTRAJHR)
!     NULTRAJBG: unit number for background (option LBACKGR)
!     NUIO_SERV_LOG: unit number of IO server log

! NB:
!   NULERR is initialised to zero in YOMLUN_IFSAUX
!   NULOUT is initialised to 6 in YOMLUN_IFSAUX, and is reset to 20 in SULUN.
!   NULDRHACK is initialised to 999 in YOMLUN_IFSAUX (see dr_hook_util.F90 for more infos)

! JP_RESERVE_FIRST = Lowest unit number that can be reserved by RESERVE_LUN
! JP_RESERVE_LAST  = Highest unit number that can be reserved by RESERVE_LUN

! NULBDESH  = offset for boundary pert spectral space
! NULBDEGP  = offset for boundary pert grid-point space


! Modifications:
!   2011-03-04 M. Fisher   : Unit numbers are now parameters.
!                            Added RESERVE_LUN and FREE_LUN.
!     ------------------------------------------------------------------

! --- Reserved for NULERR ----------------  0
INTEGER(KIND=JPIM), PARAMETER :: NULSTAT =  1
!                                           2
!                                           3
INTEGER(KIND=JPIM), PARAMETER :: NULNAM  =  4
!                                           5
! --- Reserved by Fortran                   6
INTEGER(KIND=JPIM), PARAMETER :: NULFP11 =  7
INTEGER(KIND=JPIM), PARAMETER :: NSCRTCH =  8
INTEGER(KIND=JPIM), PARAMETER :: NULFP10 =  9
INTEGER(KIND=JPIM), PARAMETER :: NULCL1  = 10
INTEGER(KIND=JPIM), PARAMETER :: NULCL2  = 11
INTEGER(KIND=JPIM), PARAMETER :: NULFP09 = 12
INTEGER(KIND=JPIM), PARAMETER :: NULFP08 = 13
INTEGER(KIND=JPIM), PARAMETER :: NULFP07 = 14
INTEGER(KIND=JPIM), PARAMETER :: NULASE  = 16
INTEGER(KIND=JPIM), PARAMETER :: NULASS  = 17
INTEGER(KIND=JPIM), PARAMETER :: NULFP06 = 18
INTEGER(KIND=JPIM), PARAMETER :: NULFP05 = 19
! --- Reserved for NULOUT ---------------- 20
INTEGER(KIND=JPIM), PARAMETER :: NULRAD  = 25
INTEGER(KIND=JPIM), PARAMETER :: NULRTL  = 26
INTEGER(KIND=JPIM), PARAMETER :: NULFP04 = 30
!                                          31
INTEGER(KIND=JPIM), PARAMETER :: NULFP12 = 32
INTEGER(KIND=JPIM), PARAMETER :: NULFP13 = 33
INTEGER(KIND=JPIM), PARAMETER :: NULFP14 = 34
INTEGER(KIND=JPIM), PARAMETER :: NULFP15 = 35
INTEGER(KIND=JPIM), PARAMETER :: NSCATAB = 36
INTEGER(KIND=JPIM), PARAMETER :: NSCASIG = 37
INTEGER(KIND=JPIM), PARAMETER :: NSCASPE = 38
INTEGER(KIND=JPIM), PARAMETER :: NEGASH  = 39
!                                          40
!                                          41
!                                          42
INTEGER(KIND=JPIM), PARAMETER :: NULFP03 = 43
!                                          44
!                                          45
INTEGER(KIND=JPIM), PARAMETER :: NULDILA = 46
INTEGER(KIND=JPIM), PARAMETER :: NULCONT = 47
!                                          48
!                                          49
INTEGER(KIND=JPIM), PARAMETER :: NULFP02 = 50
INTEGER(KIND=JPIM), PARAMETER :: NPOSSH  = 51
!                                          52
INTEGER(KIND=JPIM), PARAMETER :: NPODDH  = 53
INTEGER(KIND=JPIM), PARAMETER :: NULFP01 = 54
INTEGER(KIND=JPIM), PARAMETER :: NULCO   = 55
INTEGER(KIND=JPIM), PARAMETER :: NTIDE   = 56
! --- Reserved for RESERVE_LUN ----------- 57 
! --- Reserved for RESERVE_LUN ----------- 58
! --- Reserved for RESERVE_LUN ----------- 59
! --- Reserved for RESERVE_LUN ----------- 60
! --- Reserved for RESERVE_LUN ----------- 61
! --- Reserved for RESERVE_LUN ----------- 62
! --- Reserved for RESERVE_LUN ----------- 63
! --- Reserved for RESERVE_LUN ----------- 64
! --- Reserved for RESERVE_LUN ----------- 65
! --- Reserved for RESERVE_LUN ----------- 66
INTEGER(KIND=JPIM), PARAMETER :: NUIO_SERV_LOG = 67
INTEGER(KIND=JPIM), PARAMETER :: NTRJSH  = 71

! NB: NULTRAJHR is not a parameter because read_surfgrid_traj_fromfa resets it.
INTEGER(KIND=JPIM)            :: NULTRAJHR = 72

! NB: NEFLSS is not a parameter because ELSAC (in ALD) resets it.
INTEGER(KIND=JPIM)            :: NEFLSS  = 73

INTEGER(KIND=JPIM), PARAMETER :: NULFPOS = 74
!                                          75
!                                          76
!                                          77
INTEGER(KIND=JPIM), PARAMETER :: NEFLS   = 78
INTEGER(KIND=JPIM), PARAMETER :: NINMSH  = 79
!                                          80
INTEGER(KIND=JPIM), PARAMETER :: NINISH  = 81
INTEGER(KIND=JPIM), PARAMETER :: NINIGG  = 82
INTEGER(KIND=JPIM), PARAMETER :: NFGISH  = 83
INTEGER(KIND=JPIM), PARAMETER :: NULTRAJBG = NFGISH ! -- Dangerous!!
INTEGER(KIND=JPIM), PARAMETER :: NPRT0   = 85
INTEGER(KIND=JPIM), PARAMETER :: NPRT1   = 86
INTEGER(KIND=JPIM), PARAMETER :: NPRT2   = 87
INTEGER(KIND=JPIM), PARAMETER :: NPRT3   = 88
INTEGER(KIND=JPIM), PARAMETER :: NPRT4   = 89
INTEGER(KIND=JPIM), PARAMETER :: NBIAS   = 90
INTEGER(KIND=JPIM), PARAMETER :: NPPPSH  = 91
INTEGER(KIND=JPIM), PARAMETER :: NPDIRL  = 92
INTEGER(KIND=JPIM), PARAMETER :: NULHWF  = 93
INTEGER(KIND=JPIM), PARAMETER :: NULRCF  = 94
INTEGER(KIND=JPIM), PARAMETER :: NULUSR1 = 95
INTEGER(KIND=JPIM), PARAMETER :: NULUSR2 = 96
INTEGER(KIND=JPIM), PARAMETER :: NULUSR3 = 97
INTEGER(KIND=JPIM), PARAMETER :: NULUSR4 = 98
INTEGER(KIND=JPIM), PARAMETER :: NULUSR5 = 99
INTEGER(KIND=JPIM), PARAMETER :: NULBDEGP = 100
INTEGER(KIND=JPIM), PARAMETER :: NULBDESH = 200
! --- Reserved for NULDRHACK ----------- 999 defined in yomlun_ifsaux.F90
INTEGER(KIND=JPIM), PARAMETER :: JP_RESERVE_FIRST=57
INTEGER(KIND=JPIM), PARAMETER :: JP_RESERVE_LAST =70

INTEGER(KIND=JPIM) :: I
LOGICAL :: LRESERVED(JP_RESERVE_FIRST:JP_RESERVE_LAST) = &
        & (/ (.FALSE.,I=JP_RESERVE_FIRST,JP_RESERVE_LAST) /)
!     ------------------------------------------------------------------

CONTAINS

SUBROUTINE FOPEN(KUNIT,CDFILE,CDFORM,CDTEXT)

! Hans Hersbach 10-12-2010
! Open fortran file CDFILE and allocate unit KUNIT
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE PARKIND1  ,ONLY : JPRB
  IMPLICIT NONE

  INTEGER(KIND=JPIM)         ,INTENT(OUT) :: KUNIT
  CHARACTER(LEN=*)           ,INTENT(IN ) :: CDFILE
  CHARACTER(LEN=*)  ,OPTIONAL,INTENT(IN ) :: CDFORM
  CHARACTER(LEN=*)  ,OPTIONAL,INTENT(IN ) :: CDTEXT

! Local vars
  LOGICAL          :: LLOPENED
  INTEGER(KIND=JPIM), PARAMETER  :: FREE_FIRST = 301
  INTEGER(KIND=JPIM), PARAMETER  :: FREE_LAST  = 1000

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('YOMLUN:FOPEN',0,ZHOOK_HANDLE)
!- Assign a free unit number
   DO KUNIT=FREE_FIRST,FREE_LAST
      IF(MOD(KUNIT,100)==0)CYCLE   ! Avoid units that are a multiple of 100
      INQUIRE(UNIT=KUNIT,OPENED=LLOPENED)
      IF (.NOT.LLOPENED) EXIT
   ENDDO

!-Open file
   IF (PRESENT(CDFORM)) THEN
     OPEN(KUNIT,FILE=TRIM(CDFILE),FORM=TRIM(CDFORM))
   ELSE
     OPEN(KUNIT,FILE=TRIM(CDFILE))
   ENDIF

   IF (PRESENT(CDTEXT)) THEN
      WRITE(NULOUT,'(A,": OPEN UNIT, FILE: ",I4,", ",A)')TRIM(CDTEXT),KUNIT,TRIM(CDFILE)
   ENDIF
IF (LHOOK) CALL DR_HOOK('YOMLUN:FOPEN',1,ZHOOK_HANDLE)

END SUBROUTINE FOPEN

!     ------------------------------------------------------------------

    FUNCTION RESERVE_LUN()

!   Reserve a unit number that will not be used by anyone else.
!   The unit number remains reserved until freed by FREE_LUN.
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE PARKIND1  ,ONLY : JPRB
    INTEGER(KIND=JPIM) :: RESERVE_LUN
    INTEGER(KIND=JPIM) :: J
    LOGICAL :: LLOPEN
    CHARACTER*80 :: CLNAME
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('YOMLUN:RESERVE_LUN',0,ZHOOK_HANDLE)

    J = JP_RESERVE_FIRST
    DO WHILE (J <= JP_RESERVE_LAST)
      IF (.NOT.LRESERVED(J)) EXIT
      J=J+1
    ENDDO

    IF (J > JP_RESERVE_LAST) THEN
      CALL ABOR1 ('ALL AVAILABLE UNITS ARE RESERVED OR OPEN')
    ENDIF

    INQUIRE (UNIT=J,OPENED=LLOPEN)
    IF (LLOPEN) THEN
      WRITE (NULOUT, &
        &'(''ERROR: RESERVE_LUN THINKS UNIT '',I3,'' IS AVAILABLE'')') J
      WRITE (NULOUT,'(''BUT, THE UNIT IS OPENED FOR FILE:'')')
      INQUIRE (UNIT=J,NAME=CLNAME)
      WRITE (NULOUT,'(A80)') CLNAME
      CALL ABOR1 ('CANNOT RESERVE A UNIT THAT IS ALREADY OPEN')
    ENDIF

    LRESERVED(J) = .TRUE.
    RESERVE_LUN=J
    WRITE (NULOUT,'(''RESERVE_LUN reserves unit '',I3)') J
IF (LHOOK) CALL DR_HOOK('YOMLUN:RESERVE_LUN',1,ZHOOK_HANDLE)

  END FUNCTION RESERVE_LUN

!     ------------------------------------------------------------------

  SUBROUTINE FREE_LUN (KLUN)
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE PARKIND1  ,ONLY : JPRB
!   Return a unit number to the pool of reservable units.

    INTEGER(KIND=JPIM), INTENT(IN) :: KLUN
    LOGICAL :: LLOPEN
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('YOMLUN:FREE_LUN',0,ZHOOK_HANDLE)

    IF (KLUN < JP_RESERVE_FIRST .OR. KLUN > JP_RESERVE_LAST) THEN
      WRITE (NULOUT,'(''ERROR: UNIT NUMBER OUT OF RANGE: '',I3)') KLUN
      CALL ABOR1 ('UNIT NUMBER OUT OF RANGE')     
    ENDIF

    IF (.NOT.LRESERVED(KLUN)) THEN
      WRITE (NULOUT, &
        & '(''WARNING: FREE_LUN CALLED FOR ALREADY FREE UNIT: '',I3)') KLUN
    ENDIF

    LRESERVED(KLUN) = .FALSE.

    INQUIRE (UNIT=KLUN,OPENED=LLOPEN)
    IF (LLOPEN) THEN
      WRITE (NULOUT, &
        & '(''ERROR: FREE_LUN CALLED FOR STILL-OPEN UNIT: '',I3)') KLUN
      CALL ABOR1 ('FREE_LUN CALLED FOR A UNIT THAT IS STILL OPEN')
    ENDIF

    WRITE (NULOUT,'(''FREE_LUN returns unit '',I3,'' to the pool'')') KLUN

IF (LHOOK) CALL DR_HOOK('YOMLUN:FREE_LUN',1,ZHOOK_HANDLE)
  END SUBROUTINE FREE_LUN

!     ------------------------------------------------------------------
END MODULE YOMLUN

