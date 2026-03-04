! (C) Copyright 2000- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

INTERFACE
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
!     KDLON - number of points on each latitude [2*KDGL]
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

!     LDUSEFLT   - use Fast Legandre Transform (Butterfly algorithm)
!     LDUSERPNM  - Use Belusov to compute legendre pol. (else new alg.)
!     LDKEEPRPNM - Keep Legendre Polynomials (only applicable when using
!                  FLT, otherwise always kept)
!     LDPNMONLY  - Compute the Legendre polynomialsonly, not the FFTs.
!     LDUSEFFTW   - Use FFTW for FFTs
!     LDLL                 - Setup second set of input/output latitudes
!                                 the number of input/output latitudes to transform is equal KDGL 
!                                 or KDGL+2 in the case that includes poles + equator
!                                 the number of input/output longitudes to transform is 2*KDGL
!     LDSHIFTLL       - Shift output lon/lat data by 0.5*dx and 0.5*dy
 
!     Method.
!     -------

!     Externals.  SET_RESOL   - set resolution
!     ----------  SETUP_DIMS  - setup distribution independent dimensions
!                 SUMP_TRANS_PRELEG - first part of setup of distr. environment
!                 SULEG - Compute Legandre polonomial and Gaussian 
!                         Latitudes and Weights
!                 SUMP_TRANS - Second part of setup of distributed environment
!                 SUFFT - setup for FFT

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-03-03

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
    USE, INTRINSIC :: ISO_C_BINDING, ONLY:  C_PTR, C_INT,C_ASSOCIATED,C_SIZE_T


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
LOGICAL   ,OPTIONAL,INTENT(IN):: LDPNMONLY
LOGICAL   ,OPTIONAL,INTENT(IN):: LDSPSETUPONLY
LOGICAL   ,OPTIONAL,INTENT(IN):: LDUSEFFTW
LOGICAL   ,OPTIONAL,INTENT(IN):: LDLL
LOGICAL   ,OPTIONAL,INTENT(IN):: LDSHIFTLL
CHARACTER(LEN=*),OPTIONAL,INTENT(IN):: CDIO_LEGPOL
CHARACTER(LEN=*),OPTIONAL,INTENT(IN):: CDLEGPOLFNAME
TYPE(C_PTR) ,OPTIONAL,INTENT(IN) :: KLEGPOLPTR
INTEGER(C_SIZE_T) ,OPTIONAL,INTENT(IN) :: KLEGPOLPTR_LEN


END SUBROUTINE SETUP_TRANS


END INTERFACE
