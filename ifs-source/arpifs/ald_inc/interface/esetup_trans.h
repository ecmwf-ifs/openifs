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

INTERFACE
SUBROUTINE ESETUP_TRANS(KMSMAX,KSMAX,KDGL,KDGUX,KLOEN,LDSPLIT,&
 & KTMAX,KRESOL,PEXWN,PEYWN,PWEIGHT,LDGRIDONLY,KNOEXTZL,KNOEXTZG,&
 & LDUSEFFTW)
!**** *ESETUP_TRANS* - Setup transform package for specific resolution

!     Purpose.
!     --------
!     To setup for making spectral transforms. Each call to this routine
!     creates a new resolution up to a maximum of NMAX_RESOL set up in
!     SETUP_TRANS0. You need to call SETUP_TRANS0 before this routine can
!     be called.

!**   Interface.
!     ----------
!     CALL ESETUP_TRANS(...)

!     Explicit arguments : KLOEN,LDSPLIT are optional arguments
!     -------------------- 
!     KSMAX - spectral truncation required
!     KDGL  - number of Gaussian latitudes
!     KLOEN(:) - number of points on each Gaussian latitude [2*KDGL]
!     LDSPLIT - true if split latitudes in grid-point space [false]
!     KTMAX - truncation order for tendencies?
!     KRESOL - the resolution identifier
!     KSMAX,KDGL,KTMAX and KLOEN are GLOBAL variables desribing the resolution
!     in spectral and grid-point space
!     LDGRIDONLY - true if only grid space is required


!     LDSPLIT describe the distribution among processors of
!     grid-point data and has no relevance if you are using a single processor
 
!     LDUSEFFTW   - Use FFTW for FFTs

!     Method.
!     -------

!     Externals.  ESET_RESOL   - set resolution
!     ----------  ESETUP_DIMS  - setup distribution independent dimensions
!                 SUEMP_TRANS_PRELEG - first part of setup of distr. environment
!                 SULEG - Compute Legandre polonomial and Gaussian 
!                         Latitudes and Weights
!                 ESETUP_GEOM - Compute arrays related to grid-point geometry
!                 SUEMP_TRANS - Second part of setup of distributed environment
!                 SUEFFT - setup for FFT

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-03-03
!        02-04-11 A. Bogatchev: Passing of TCDIS
!        02-11-14 C. Fischer: soften test on KDGL
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        A.Nmiri       15-Nov-2007 Phasing with TFL 32R3
!        A.Bogatchev   16-Sep-2010 Phasing cy37
!        D. Degrauwe  (Feb 2012): Alternative extension zone (E')

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

! Dummy arguments
INTEGER(KIND=JPIM),INTENT(IN)    :: KMSMAX 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSMAX 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDGL 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDGUX 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLOEN(:) 
LOGICAL           ,OPTIONAL,INTENT(IN)    :: LDSPLIT 
LOGICAL           ,OPTIONAL,INTENT(IN)    :: LDGRIDONLY
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN)    :: KTMAX 
INTEGER(KIND=JPIM),OPTIONAL,INTENT(INOUT) :: KRESOL 
REAL(KIND=JPRB)   ,OPTIONAL,INTENT(IN)    :: PEXWN 
REAL(KIND=JPRB)   ,OPTIONAL,INTENT(IN)    :: PEYWN 
REAL(KIND=JPRB)   ,OPTIONAL,INTENT(IN)    :: PWEIGHT(:)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN)    :: KNOEXTZL
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN)    :: KNOEXTZG
LOGICAL   ,OPTIONAL,INTENT(IN)            :: LDUSEFFTW

END SUBROUTINE ESETUP_TRANS
END INTERFACE
