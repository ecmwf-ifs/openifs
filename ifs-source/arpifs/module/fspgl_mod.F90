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

MODULE FSPGL_MOD

! Module which stores values necessary for the call to FSPGLH from Fourier space
! Values are filled in by call to FSPGL_FILL inside transinv_mdl

USE PARKIND1, ONLY : JPIM,JPRB

INTEGER(KIND=JPRB), ALLOCATABLE, SAVE :: NFSPGL_MYLEVS(:)
REAL(KIND=JPRB)   , ALLOCATABLE, SAVE :: RFSPGL_KRF(:)
REAL(KIND=JPRB)   ,              SAVE :: RFSPGL_TSTEP

CONTAINS

SUBROUTINE FSPGL_FILL(KMYLEVS,PKRF,PTSTEP)  

!
!        KMYLEVS - list of vertical levels handled by this processor/PE
!        PKRF    - vertical profile of Rayleigh friction values
!        PTSTEP  - time step

  
!     Externals.   None.
!     ----------
!     Author.
!     -------
!        Olivier Marsden *ECMWF*
!        Original : June   2017

!     Modifications.
!     --------------
!     ------------------------------------------------------------------

USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN) :: KMYLEVS(:)
REAL(KIND=JPRB)   , INTENT(IN) :: PKRF(:)
REAL(KIND=JPRB)   , INTENT(IN) :: PTSTEP

REAL(KIND=JPHOOK)  :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('FSPGL_FILL',0,ZHOOK_HANDLE)


IF (ALLOCATED(NFSPGL_MYLEVS)) DEALLOCATE(NFSPGL_MYLEVS)
IF (ALLOCATED(RFSPGL_KRF))    DEALLOCATE(RFSPGL_KRF)

ALLOCATE( NFSPGL_MYLEVS(LBOUND(KMYLEVS,1):UBOUND(KMYLEVS,1)) )
ALLOCATE( RFSPGL_KRF   (LBOUND(PKRF,1)   :UBOUND(PKRF,1)   ) )

NFSPGL_MYLEVS = KMYLEVS
RFSPGL_KRF    = PKRF
RFSPGL_TSTEP  = PTSTEP

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('FSPGL_FILL',1,ZHOOK_HANDLE)
END SUBROUTINE FSPGL_FILL

END MODULE FSPGL_MOD
