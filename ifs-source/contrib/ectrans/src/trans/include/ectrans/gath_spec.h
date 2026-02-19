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
SUBROUTINE GATH_SPEC(PSPECG,KFGATHG,KTO,KVSET,KRESOL,PSPEC,LDIM1_IS_FLD,KSMAX,LDZA0IP)

!**** *GATH_SPEC* - Gather global spectral array from processors

!     Purpose.
!     --------
!        Interface routine for gathering spectral array

!**   Interface.
!     ----------
!     CALL GATH_SPEC(...)

!     Explicit arguments : 
!     -------------------- 
!     PSPECG(:,:) - Global spectral array
!     KFGATHG     - Global number of fields to be gathered
!     KTO(:)      - Processor responsible for gathering each field
!     KVSET(:)    - "B-Set" for each field
!     KRESOL      - resolution tag  which is required ,default is the
!                   first defined resulution (input)
!     PSPEC(:,:)  - Local spectral array
!     LDIM1_IS_FLD - If TRUE first dimension of PSCPEC and PSPECG is the field dimension [.T.]
!
!     Method.
!     -------

!     Externals.  SET_RESOL   - set resolution
!     ----------  GATH_SPEC_CONTROL - control routine

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-03-03

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB


IMPLICIT NONE

! Declaration of arguments

REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PSPECG(:,:)
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KFGATHG
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KTO(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KVSET(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KRESOL
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(IN)  :: PSPEC(:,:)
LOGICAL            ,OPTIONAL, INTENT(IN)  :: LDIM1_IS_FLD
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KSMAX
LOGICAL            ,OPTIONAL, INTENT(IN)  :: LDZA0IP


!     ------------------------------------------------------------------

END SUBROUTINE GATH_SPEC

END INTERFACE
