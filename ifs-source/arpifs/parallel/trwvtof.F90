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

SUBROUTINE TRWVTOF(KRESOL,PSPECG,KFGATHG,KTO,KVSET,PSPEC)

!**** *TRWVTOF*  - TRANSPOSE A WAVE-DISTRIBUTED SET OF FIELDS 
!                  TO A FIELD-DISTRIBUTED SET OF SPECTRA

!     PURPOSE.
!     --------
!        To communicate a set of global spectral fields from W/V-set processors
!        to the processors in charge of global work prior to writing out fields 
!        (usually : packing).

!**   INTERFACE.
!     ----------
!       *CALL* *TRWVTOF*

!        EXPLICIT ARGUMENTS
!        ------------------
!     KRESOL      - resolution tag  which is required
!     PSPECG(:,:) - Global spectral array
!     KFGATHG     - Global number of fields to be gathered
!     KTO(:)      - Processor responsible for gathering each field
!     KVSET(:)    - "B-Set" for each field
!     PSPEC(:,:)  - Local spectral array

!        IMPLICIT ARGUMENTS
!        ------------------
!        See #include below.

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     AUTHOR.
!     -------
!      RYAD EL KHATIB *METEO-FRANCE*
!      ORIGINAL : 99-02-05

!     MODIFICATIONS.
!     --------------
!      R. El Khatib : 02-21-20 Fullpos B-level distribution
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      R. El Ouaraini & R. El Khatib: 03-Oct-06: The received global spectral array is ordered in according to JM=0,1,2,3,4,...
!      R. El Khatib : 20-Aug-2012 Rewrite
!      P.Marguinaud : Use LDIM1_IS_FLD=.FALSE.
!     --------------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK

IMPLICIT NONE

INTEGER(KIND=JPIM) ,INTENT(IN)  :: KRESOL
REAL(KIND=JPRB)    ,INTENT(IN)  :: PSPEC(:,:)
INTEGER(KIND=JPIM) ,INTENT(IN)  :: KFGATHG
INTEGER(KIND=JPIM) ,INTENT(IN)  :: KTO(:)
INTEGER(KIND=JPIM) ,INTENT(IN)  :: KVSET(:)
REAL(KIND=JPRB)    ,INTENT(OUT) :: PSPECG(:,:)

LOGICAL :: LLETRANS
REAL(KIND=JPRB), ALLOCATABLE :: ZSPEC (:,:)
INTEGER (KIND=JPIM) :: JFLD, JS
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "gath_spec.h"
#include "egath_spec.h"

#include "updtrans.intfb.h"
!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('TRWVTOF',0,ZHOOK_HANDLE)

ALLOCATE (ZSPEC (SIZE (PSPEC, 2), SIZE (PSPEC, 1)))

DO JFLD = 1, SIZE (PSPEC, 1)
  DO JS = 1, SIZE (PSPEC, 2)
    ZSPEC (JS, JFLD) = PSPEC (JFLD, JS)
  ENDDO
ENDDO

CALL UPDTRANS(KRESOL,LLETRANS)
IF (LLETRANS) THEN
  CALL EGATH_SPEC(KRESOL=KRESOL,PSPECG=PSPECG(:,:),KFGATHG=KFGATHG,KTO=KTO(:), &
   & KVSET=KVSET,PSPEC=ZSPEC(:,:),LDIM1_IS_FLD=.FALSE.)
ELSE
  CALL GATH_SPEC(KRESOL=KRESOL,PSPECG=PSPECG(:,:),KFGATHG=KFGATHG,KTO=KTO(:), &
   & KVSET=KVSET,PSPEC=ZSPEC(:,:),LDIM1_IS_FLD=.FALSE.)
ENDIF

DEALLOCATE (ZSPEC)

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('TRWVTOF',1,ZHOOK_HANDLE)
END SUBROUTINE TRWVTOF
