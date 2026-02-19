! (C) Copyright 2000- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE SPECNORM(PSPEC,KVSET,KMASTER,KRESOL,PMET,PNORM)

!**** *SPECNORM* - Compute global spectral norms

!     Purpose.
!     --------
!        Interface routine for computing spectral norms

!**   Interface.
!     ----------
!     CALL SPECNORM(...)

!     Explicit arguments : All arguments optional
!     --------------------
!     PSPEC(:,:)  - Spectral array
!     KVSET(:)    - "B-Set" for each field
!     KMASTER     - processor to recieve norms
!     KRESOL      - resolution tag  which is required ,default is the
!                   first defined resulution (input)
!     PMET(:)     - metric
!     PNORM(:)    - Norms (output for processor KMASTER)
!
!     Method.
!     -------

!     Externals.  SET_RESOL - set resolution
!     ----------  SPNORM_CTL - control routine

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-03-03

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

!ifndef INTERFACE

USE TPM_GEN         ,ONLY : NERR
!USE TPM_DIM
USE TPM_DISTR       ,ONLY : D, NPRTRV, MYSETV, MYPROC

USE SET_RESOL_MOD   ,ONLY : SET_RESOL
USE SPNORM_CTL_MOD  ,ONLY : SPNORM_CTL
USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS

!endif INTERFACE

IMPLICIT NONE

! Declaration of arguments


REAL(KIND=JPRB)    ,OPTIONAL, INTENT(IN)  :: PSPEC(:,:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KVSET(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KMASTER
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(IN)  :: PMET(:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PNORM(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KRESOL
!ifndef INTERFACE

INTEGER(KIND=JPIM) :: IMASTER,IFLD,IFLD_G,J

!     ------------------------------------------------------------------

! Set current resolution
CALL SET_RESOL(KRESOL)

! Set defaults
IMASTER = 1
IFLD    = 0


IF(PRESENT(KMASTER)) THEN
  IMASTER = KMASTER
ENDIF

IF(PRESENT(KVSET)) THEN
  IFLD_G = UBOUND(KVSET,1)
  DO J=1,IFLD_G
    IF(KVSET(J) > NPRTRV) THEN
      WRITE(NERR,*) 'SPECNORM:KVSET(J) > NPRTRV ',J,KVSET(J),NPRTRV
      CALL ABORT_TRANS('SPECNORM:KVSET TOO LONG OR CONTAINS VALUES OUTSIDE RANGE')
    ENDIF
    IF(KVSET(J) == MYSETV) THEN
      IFLD = IFLD+1
    ENDIF
  ENDDO
ELSE
  IF(PRESENT(PSPEC)) THEN
    IFLD = UBOUND(PSPEC,1)
  ENDIF
  IFLD_G = IFLD
ENDIF

IF(NPRTRV >1) THEN
  IF(IFLD > 0 .AND. .NOT. PRESENT(KVSET)) THEN
    WRITE(NERR,*)'NPRTRV >1 AND IFLD > 0 AND NOT PRESENT(KVSET)',&
                 &NPRTRV,IFLD
    CALL ABORT_TRANS('SPECNORM: SPECIFY VERTICAL SPECTRAL DISTRIBUTION!')
  ENDIF
ENDIF
IF(MYPROC == IMASTER) THEN
  IF(.NOT. PRESENT(PNORM)) THEN
    CALL ABORT_TRANS('SPECNORM: PNORM NOT PRESENT')
  ENDIF
  IF(UBOUND(PNORM,1) < IFLD_G) THEN
    CALL ABORT_TRANS('SPECNORM: PNORM TOO SMALL')
  ENDIF
ENDIF
IF(IFLD > 0 ) THEN
  IF(.NOT. PRESENT(PSPEC)) THEN
    CALL ABORT_TRANS('SPECNORM: PSPEC NOT PRESENT')
  ENDIF
  IF(UBOUND(PSPEC,1) < IFLD) THEN
    CALL ABORT_TRANS('SPECNORM: FIRST DIMENSION OF PSPEC TOO SMALL')
  ENDIF
  IF(UBOUND(PSPEC,2) < D%NSPEC2) THEN
    CALL ABORT_TRANS('SPECNORM: FIRST DIMENSION OF PSPEC TOO SMALL')
  ENDIF
ENDIF

CALL SPNORM_CTL(PSPEC,IFLD,IFLD_G,KVSET,IMASTER,PMET,PNORM)

!endif INTERFACE

!     ------------------------------------------------------------------

END SUBROUTINE SPECNORM

