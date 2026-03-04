! (C) Copyright 2000- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE SETUP_TRANS(KSMAX,KDGL,KDLON,KLOEN,LDSPLIT,PSTRET,&
&KTMAX,KRESOL,PWEIGHT,LDGRIDONLY,LDUSERPNM,LDKEEPRPNM,LDUSEFLT,&
&LDSPSETUPONLY,LDPNMONLY,LDUSEFFTW,&
&LDLL,LDSHIFTLL,CDIO_LEGPOL,CDLEGPOLFNAME,KLEGPOLPTR,KLEGPOLPTR_LEN)

!**** *SETUP_TRANS* - Setup transform package for specific resolution

!     Purpose.
!     --------
!     To setup for making spectral transforms. Each call to this routine
!     creates a new resolution up to a maximum of NMAX_RESOL set up in
!     SETUP_TRANS0. You need to call SETUP_TRANS0 before this routine can
!     be called.

!**   Interface.
!     ----------
!     CALL SETUP_TRANS(...)

!     Explicit arguments : KLOEN,LDSPLIT are optional arguments
!     --------------------
!     KSMAX - spectral truncation required
!     KDGL  - number of Gaussian latitudes
!     KDLON - number of points on each Gaussian latitude [2*KDGL]
!     KLOEN(:) - number of points on each Gaussian latitude [2*KDGL]
!     LDSPLIT - true if split latitudes in grid-point space [false]
!     KTMAX - truncation order for tendencies?
!     KRESOL - the resolution identifier
!     PWEIGHT - the weight per grid-point (for a weighted distribution)
!     LDGRIDONLY - true if only grid space is required

!     KSMAX,KDGL,KTMAX and KLOEN are GLOBAL variables desribing the resolution
!     in spectral and grid-point space

!     LDSPLIT describe the distribution among processors of grid-point data and
!     has no relevance if you are using a single processor

!     PSTRET     - stretching factor - for the case the Legendre polynomials are
!                  computed on the stretched sphere - works with LSOUTHPNM
!     LDUSEFLT   - use Fast Legandre Transform (Butterfly algorithm)
!     LDUSERPNM  - Use Belusov algorithm to compute legendre pol. (else new alg.)
!     LDKEEPRPNM - Keep Legendre Polynomials (only applicable when using
!                  FLT, otherwise always kept)
!     LDPNMONLY  - Compute the Legendre polynomials only, not the FFTs.
!     LDUSEFFTW    - Use FFTW for FFTs
!     LDLL                 - Setup second set of input/output latitudes
!                                 the number of input/output latitudes to transform is equal KDGL
!                                 or KDGL+2 in the case that includes poles + equator
!                                 the number of input/output longitudes to transform is 2*KDGL
!     LDSHIFTLL       - Shift output lon/lat data by 0.5*dx and 0.5*dy
!     CDIO_LEGPOL  - IO option on Legendre polinomials :  N.B. Only works for NPROC=1
!                    Options:
!                    'READF' -  read Leg.Pol. from file CDLEGPOLFNAME
!                    'WRITEF' - write Leg.Pol. to file CDLEGPOLFNAME
!                    'MEMBUF' - Leg. Pol provided in shared memory segment pointed to by KLEGPOLPTR of
!                               length KLEGPOLPTR_LEN
!     CDLEGPOLFNAME - file name for Leg.Pol. IO
!     KLEGPOLPTR    - pointer to Legendre polynomials memory segment
!     KLEGPOLPTR_LEN  - length of  Legendre polynomials memory segment

!     Method.
!     -------

!     Externals.  SET_RESOL   - set resolution
!     ----------  SETUP_DIMS  - setup distribution independent dimensions
!                 SUMP_TRANS_PRELEG - first part of setup of distr. environment
!                 SULEG - Compute Legandre polonomial and Gaussian
!                         Latitudes and Weights
!                 SUMP_TRANS - Second part of setup of distributed environment
!                 SUFFT - setup for FFT
!                 SHAREDMEM_CREATE - create memory buffer for Leg.pol.

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-03-03
!        Daan Degrauwe : Mar 2012 E'-zone dimensions
!        R. El Khatib 09-Aug-2012 %LAM in GEOM_TYPE
!        R. El Khatib 14-Jun-2013 PSTRET, LDPNMONLY, LENABLED
!        G. Mozdzynski : Oct 2014 Support f
!        N. Wedi       : Apr 2015 Support dual set of lat/lon
!        G. Mozdzynski : Jun 2015 Support alternative FFTs to FFTW
!        M.Hamrud/W.Deconinck : July 2015 IO options for Legenndre polynomials
!        R. El Khatib 07-Mar-2016 Better flexibility for Legendre polynomials computation in stretched mode
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB,  JPRD
USE, INTRINSIC :: ISO_C_BINDING, ONLY:  C_PTR, C_INT,C_ASSOCIATED,C_SIZE_T

!ifndef INTERFACE

USE TPM_GEN         ,ONLY : NOUT, MSETUP0, NCUR_RESOL, NDEF_RESOL, &
     &                      NMAX_RESOL, NPRINTLEV, LENABLED, NERR
USE TPM_DIM         ,ONLY : R, DIM_RESOL
USE TPM_DISTR       ,ONLY : D, DISTR_RESOL,NPROC
USE TPM_GEOMETRY    ,ONLY : G, GEOM_RESOL
USE TPM_FIELDS      ,ONLY : FIELDS_RESOL
USE TPM_FFT         ,ONLY : T, FFT_RESOL, TB, FFTB_RESOL
#ifdef WITH_FFTW
USE TPM_FFTW        ,ONLY : TW, FFTW_RESOL
#endif
USE TPM_FLT
USE TPM_CTL

USE SET_RESOL_MOD   ,ONLY : SET_RESOL
USE SETUP_DIMS_MOD  ,ONLY : SETUP_DIMS
USE SUMP_TRANS_MOD  ,ONLY : SUMP_TRANS
USE SUMP_TRANS_PRELEG_MOD ,ONLY : SUMP_TRANS_PRELEG
USE SULEG_MOD       ,ONLY : SULEG
USE PRE_SULEG_MOD   ,ONLY : PRE_SULEG
USE SUFFT_MOD       ,ONLY : SUFFT
USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS
USE SHAREDMEM_MOD    ,ONLY : SHAREDMEM_CREATE
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

!endif INTERFACE

IMPLICIT NONE

! Dummy arguments

INTEGER(KIND=JPIM) ,INTENT(IN) :: KSMAX,KDGL
INTEGER(KIND=JPIM) ,OPTIONAL,INTENT(IN) :: KDLON
INTEGER(KIND=JPIM) ,OPTIONAL,INTENT(IN) :: KLOEN(:)
LOGICAL   ,OPTIONAL,INTENT(IN) :: LDSPLIT
INTEGER(KIND=JPIM) ,OPTIONAL,INTENT(IN) :: KTMAX
INTEGER(KIND=JPIM) ,OPTIONAL,INTENT(OUT):: KRESOL
REAL(KIND=JPRB)    ,OPTIONAL,INTENT(IN) :: PWEIGHT(:)
REAL(KIND=JPRB)    ,OPTIONAL,INTENT(IN) :: PSTRET
LOGICAL   ,OPTIONAL,INTENT(IN):: LDGRIDONLY
LOGICAL   ,OPTIONAL,INTENT(IN):: LDUSEFLT
LOGICAL   ,OPTIONAL,INTENT(IN):: LDUSERPNM
LOGICAL   ,OPTIONAL,INTENT(IN):: LDKEEPRPNM
LOGICAL   ,OPTIONAL,INTENT(IN):: LDSPSETUPONLY
LOGICAL   ,OPTIONAL,INTENT(IN):: LDPNMONLY
LOGICAL   ,OPTIONAL,INTENT(IN):: LDUSEFFTW
LOGICAL   ,OPTIONAL,INTENT(IN):: LDLL
LOGICAL   ,OPTIONAL,INTENT(IN):: LDSHIFTLL
CHARACTER(LEN=*),OPTIONAL,INTENT(IN):: CDIO_LEGPOL
CHARACTER(LEN=*),OPTIONAL,INTENT(IN):: CDLEGPOLFNAME
TYPE(C_PTR) ,OPTIONAL,INTENT(IN) :: KLEGPOLPTR
INTEGER(C_SIZE_T) ,OPTIONAL,INTENT(IN) :: KLEGPOLPTR_LEN

!ifndef INTERFACE

! Local variables
INTEGER(KIND=JPIM) :: JGL,JRES,IDEF_RESOL

LOGICAL :: LLP1,LLP2, LLSPSETUPONLY
REAL(KIND=JPRD)    :: ZTIME0,ZTIME1,ZTIME2
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "user_clock.intfb.h"
!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SETUP_TRANS',0,ZHOOK_HANDLE)

IF(MSETUP0 == 0) THEN
  CALL ABORT_TRANS('SETUP_TRANS: SETUP_TRANS0 HAS TO BE CALLED BEFORE SETUP_TRANS')
ENDIF
LLP1 = NPRINTLEV>0
LLP2 = NPRINTLEV>1
IF(LLP1) WRITE(NOUT,*) '=== ENTER ROUTINE SETUP_TRANS ==='

! Allocate resolution dependent structures
IF(.NOT. ALLOCATED(DIM_RESOL)) THEN
  IDEF_RESOL = 1
  ALLOCATE(DIM_RESOL(NMAX_RESOL))
  ALLOCATE(FIELDS_RESOL(NMAX_RESOL))
  ALLOCATE(GEOM_RESOL(NMAX_RESOL))
  ALLOCATE(DISTR_RESOL(NMAX_RESOL))
  ALLOCATE(FFT_RESOL(NMAX_RESOL))
  ALLOCATE(FFTB_RESOL(NMAX_RESOL))
#ifdef WITH_FFTW
  ALLOCATE(FFTW_RESOL(NMAX_RESOL))
#endif
  ALLOCATE(FLT_RESOL(NMAX_RESOL))
  ALLOCATE(CTL_RESOL(NMAX_RESOL))
  GEOM_RESOL(:)%LAM=.FALSE.
  ALLOCATE(LENABLED(NMAX_RESOL))
  LENABLED(:)=.FALSE.
ELSE
  IDEF_RESOL = NMAX_RESOL+1
  DO JRES=1,NMAX_RESOL
    IF(.NOT.LENABLED(JRES)) THEN
      IDEF_RESOL = JRES
      EXIT
    ENDIF
  ENDDO
  IF(IDEF_RESOL > NMAX_RESOL) THEN
    CALL ABORT_TRANS('SETUP_TRANS:IDEF_RESOL > NMAX_RESOL')
  ENDIF
ENDIF

IF (PRESENT(KRESOL)) THEN
  KRESOL=IDEF_RESOL
ENDIF

! Point at structures due to be initialized
CALL SET_RESOL(IDEF_RESOL,LDSETUP=.TRUE.)

IF(LLP1) WRITE(NOUT,*) '=== DEFINING RESOLUTION ',NCUR_RESOL



! Defaults for optional arguments


G%LREDUCED_GRID = .FALSE.
G%RSTRET=1.0_JPRB
D%LGRIDONLY = .FALSE.
D%LSPLIT = .FALSE.
D%LCPNMONLY=.FALSE.
S%LUSE_BELUSOV=.TRUE. ! use Belusov algorithm to compute RPNM array instead of per m
S%LKEEPRPNM=.FALSE. ! Keep Legendre polonomials (RPNM)
S%LUSEFLT=.FALSE. ! Use fast legendre transforms
#ifdef WITH_FFTW
TW%LFFTW=.FALSE. ! Use FFTW interface for FFTs
#endif
LLSPSETUPONLY = .FALSE. ! Only create distributed spectral setup
S%LDLL = .FALSE. ! use mapping to/from second set of latitudes
S%LSHIFTLL = .FALSE. ! shift output lat-lon by 0.5dx, 0.5dy
C%LREAD_LEGPOL = .FALSE.
C%LWRITE_LEGPOL = .FALSE.


! NON-OPTIONAL ARGUMENTS
R%NSMAX = KSMAX
R%NDGL  = KDGL
! E'-defaults
R%NNOEXTZL=0
R%NNOEXTZG=0

! IMPLICIT argument :
G%LAM = .FALSE.

IF(PRESENT(KDLON)) THEN
  R%NDLON = KDLON
ELSE
  R%NDLON = 2*R%NDGL
ENDIF

IF(PRESENT(LDLL)) THEN
  S%LDLL=LDLL
  IF( LDLL ) THEN
    S%NDLON=R%NDLON
    ! account for pole + equator
    R%NDGL=R%NDGL+2
    IF(PRESENT(LDSHIFTLL)) THEN
      S%LSHIFTLL = LDSHIFTLL
      ! geophysical (shifted) lat-lon without pole and equator
      IF(S%LSHIFTLL) R%NDGL=R%NDGL-2
    ENDIF
    S%NDGL=R%NDGL
  ENDIF
ENDIF

IF (R%NDGL <= 0 .OR. MOD(R%NDGL,2) /= 0) THEN
  CALL ABORT_TRANS ('SETUP_TRANS: KDGL IS NOT A POSITIVE, EVEN NUMBER')
ENDIF

! Optional arguments

ALLOCATE(G%NLOEN(R%NDGL))
IF(LLP2)WRITE(NOUT,9) 'NLOEN   ',SIZE(G%NLOEN   ),SHAPE(G%NLOEN   )
IF(PRESENT(KLOEN)) THEN
  IF( MINVAL(KLOEN(:)) <= 0 )THEN
     CALL ABORT_TRANS ('SETUP_TRANS: KLOEN INVALID (ONE or MORE POINTS <= 0)')
  ENDIF
  R%NDLON=MAXVAL(KLOEN(:))
  DO JGL=1,R%NDGL
    IF(KLOEN(JGL) /= R%NDLON) THEN
      G%LREDUCED_GRID = .TRUE.
      EXIT
    ENDIF
  ENDDO
ENDIF

IF (G%LREDUCED_GRID) THEN
  G%NLOEN(:) = KLOEN(1:R%NDGL)
ELSE
  G%NLOEN(:) = R%NDLON
ENDIF

IF(PRESENT(LDSPLIT)) THEN
  D%LSPLIT = LDSPLIT
ENDIF

IF(PRESENT(KTMAX)) THEN
  R%NTMAX = KTMAX
ELSE
  R%NTMAX = R%NSMAX
ENDIF

IF(PRESENT(PWEIGHT)) THEN
  D%LWEIGHTED_DISTR = .TRUE.
  IF( D%LWEIGHTED_DISTR .AND. .NOT.D%LSPLIT )THEN
    CALL ABORT_TRANS('SETUP_TRANS: LWEIGHTED_DISTR=T AND LSPLIT=F NOT SUPPORTED')
  ENDIF
  IF(SIZE(PWEIGHT) /= SUM(G%NLOEN(:)) )THEN
    CALL ABORT_TRANS('SETUP_TRANS:SIZE(PWEIGHT) /= SUM(G%NLOEN(:))')
  ENDIF
  IF( MINVAL(PWEIGHT(:)) < 0.0_JPRB )THEN
    CALL ABORT_TRANS('SETUP_TRANS: INVALID WEIGHTS')
  ENDIF
  ALLOCATE(D%RWEIGHT(SIZE(PWEIGHT)))
  D%RWEIGHT(:)=PWEIGHT(:)
ELSE
  D%LWEIGHTED_DISTR = .FALSE.
ENDIF

IF(PRESENT(LDGRIDONLY)) THEN
  D%LGRIDONLY=LDGRIDONLY
ENDIF

IF(PRESENT(LDSPSETUPONLY)) THEN
  LLSPSETUPONLY=LDSPSETUPONLY
ENDIF

IF(PRESENT(LDPNMONLY)) THEN
  D%LCPNMONLY=LDPNMONLY
ENDIF


#ifdef WITH_FFTW
IF(PRESENT(LDUSEFFTW)) THEN
  TW%LFFTW=LDUSEFFTW
ENDIF
IF( LLSPSETUPONLY .OR. D%LGRIDONLY ) THEN
  TW%LFFTW = .FALSE.
ENDIF
#endif

S%LSOUTHPNM=.FALSE.
IF(PRESENT(PSTRET)) THEN
  IF (ABS(PSTRET-1.0_JPRB)>100._JPRB*EPSILON(1._JPRB)) THEN
    G%RSTRET=PSTRET
    S%LSOUTHPNM=.TRUE.
  ENDIF
ENDIF

IF(PRESENT(CDIO_LEGPOL)) THEN
  IF(NPROC > 1) CALL  ABORT_TRANS('SETUP_TRANS:CDIO_LEGPOL OPTIONS ONLY FOR NPROC=1 ')
  IF(R%NSMAX > 511 ) S%LUSEFLT = .TRUE. !To save IO and memory
  IF(TRIM(CDIO_LEGPOL) == 'readf' .OR. TRIM(CDIO_LEGPOL) == 'READF' ) THEN
    IF(.NOT.PRESENT(CDLEGPOLFNAME)) CALL  ABORT_TRANS('SETUP_TRANS: CDLEGPOLFNAME ARGUMENT MISSING')
    C%LREAD_LEGPOL = .TRUE.
    C%CLEGPOLFNAME = TRIM(CDLEGPOLFNAME)
    C%CIO_TYPE='file'
  ELSEIF(TRIM(CDIO_LEGPOL) == 'writef' .OR. TRIM(CDIO_LEGPOL) == 'WRITEF') THEN
    IF(.NOT.PRESENT(CDLEGPOLFNAME)) CALL  ABORT_TRANS('SETUP_TRANS: CDLEGPOLFNAME ARGUMENT MISSING')
    C%LWRITE_LEGPOL = .TRUE.
    C%CLEGPOLFNAME = TRIM(CDLEGPOLFNAME)
    C%CIO_TYPE='file'
  ELSEIF(TRIM(CDIO_LEGPOL) == 'membuf' .OR. TRIM(CDIO_LEGPOL) == 'MEMBUF') THEN
    IF(.NOT.PRESENT(KLEGPOLPTR)) CALL  ABORT_TRANS('SETUP_TRANS: KLEGPOLPTR  ARGUMENT MISSING')
    IF(.NOT.C_ASSOCIATED(KLEGPOLPTR))  CALL  ABORT_TRANS('SETUP_TRANS: KLEGPOLPTR NULL POINTER')
    IF(.NOT.PRESENT(KLEGPOLPTR_LEN)) CALL  ABORT_TRANS('SETUP_TRANS: KLEGPOLPTR_LEN ARGUMENT MISSING')
    C%LREAD_LEGPOL = .TRUE.
    C%CIO_TYPE='mbuf'
    CALL SHAREDMEM_CREATE( C%STORAGE,KLEGPOLPTR,KLEGPOLPTR_LEN)
  ELSE
    WRITE(NERR,*) 'CDIO_LEGPOL ', TRIM(CDIO_LEGPOL)
    CALL  ABORT_TRANS('SETUP_TRANS:CDIO_LEGPOL UNKNOWN METHOD ')
  ENDIF
ENDIF

IF(PRESENT(LDUSEFLT)) THEN
  S%LUSEFLT=LDUSEFLT
ENDIF
IF(PRESENT(LDUSERPNM)) THEN
  S%LUSE_BELUSOV=LDUSERPNM
ENDIF
IF(PRESENT(LDKEEPRPNM)) THEN
  IF(S%LUSEFLT) THEN
    IF(LDKEEPRPNM.AND..NOT.LDUSERPNM) THEN
      CALL ABORT_TRANS('SETUP_TRANS: LDKEEPRPNM=true with LDUSERPNM=false')
    ENDIF
  ENDIF
  S%LKEEPRPNM=LDKEEPRPNM
ENDIF
!     Setup resolution dependent structures
!     -------------------------------------

! Setup distribution independent dimensions
CALL SETUP_DIMS

! First part of setup of distributed environment
CALL SUMP_TRANS_PRELEG

IF( .NOT.LLSPSETUPONLY ) THEN

! Compute Legendre polonomial and Gaussian Latitudes and Weights
  CALL SULEG

! Second part of setup of distributed environment
  CALL SUMP_TRANS
  CALL GSTATS(1802,0)

! Initialize Fast Fourier Transform package
  IF (.NOT.D%LCPNMONLY) CALL SUFFT
  CALL GSTATS(1802,1)
ELSE
  CALL PRE_SULEG
ENDIF

! Signal the current resolution is active
LENABLED(IDEF_RESOL)=.TRUE.
NDEF_RESOL = COUNT(LENABLED)

IF (LHOOK) CALL DR_HOOK('SETUP_TRANS',1,ZHOOK_HANDLE)
!     ------------------------------------------------------------------
9 FORMAT(1X,'ARRAY ',A10,' ALLOCATED ',8I8)

!endif INTERFACE

END SUBROUTINE SETUP_TRANS
