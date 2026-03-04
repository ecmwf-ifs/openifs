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
SUBROUTINE GPNORM_TRANS(PGP,KFIELDS,KPROMA,PAVE,PMIN,PMAX,LDAVE_ONLY,KRESOL)


!**** *GPNORM_TRANS* - calculate grid-point norms

!     Purpose.
!     --------
!        calculate grid-point norms using a 2 stage (NPRTRV,NPRTRW) communication rather
!        than an approach using a more expensive global gather collective communication

!**   Interface.
!     ----------
!     CALL GPNORM_TRANS(...)

!     Explicit arguments :
!     --------------------
!     PGP(:,:,:) - gridpoint fields (input)
!                  PGP is  dimensioned (NPROMA,KFIELDS,NGPBLKS) where
!                  NPROMA is the blocking factor, KFIELDS the total number
!                  of fields and NGPBLKS the number of NPROMA blocks.
!     KFIELDS     - number of fields (input)
!                   (these do not have to be just levels)
!     KPROMA      - required blocking factor (input)
!     PAVE        - average (output)
!     PMIN        - minimum (input/output)
!     PMAX        - maximum (input/output)
!     LDAVE_ONLY  - T : PMIN and PMAX already contain local MIN and MAX
!     KRESOL      -  resolution tag (optional)
!                    default assumes first defined resolution
!

!     Author.
!     -------
!        George Mozdzynski *ECMWF*

!     Modifications.
!     --------------
!        Original : 19th Sept 2008

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

! Declaration of arguments
  
REAL(KIND=JPRB),INTENT(IN)    :: PGP(:,:,:)
REAL(KIND=JPRB),INTENT(OUT)   :: PAVE(:)
REAL(KIND=JPRB),INTENT(INOUT) :: PMIN(:)
REAL(KIND=JPRB),INTENT(INOUT) :: PMAX(:)
INTEGER(KIND=JPIM),INTENT(IN) :: KFIELDS
INTEGER(KIND=JPIM),INTENT(IN) :: KPROMA
LOGICAL,INTENT(IN)            :: LDAVE_ONLY
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KRESOL

END SUBROUTINE GPNORM_TRANS
END INTERFACE
