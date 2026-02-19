! (C) Copyright 1998- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE PE2SET_MOD
CONTAINS
SUBROUTINE PE2SET(KPE,KPRGPNS,KPRGPEW,KPRTRW,KPRTRV)


!**** *PE2SET* - Convert from PE number to set numbers

!     Purpose.
!     --------
!        Convert from PE number to set numbers in both
!                  grid-point space and spectral space

!**   Interface.
!     ----------
!        *CALL* *PE2SET(KPE,KPRGPNS,KPRGPEW,KPRTRW,KPRTRV)

!        Explicit arguments :
!        --------------------
!                  input:   KPE     - integer processor number
!                                     in the range 1 .. NPROC
!                  output:  KPRGPNS - integer A set number in grid space
!                                     in the range 1 .. NPRGPNS
!                           KPRGPEW - integer B set number in grid space
!                                     in the range 1 .. NPRGPEW
!                           KPRTRW  - integer A set number in spectral space
!                                     in the range 1 .. NPRTRW
!                           KPRTRV  - integer B set number in spectral space
!                                     in the range 1 .. NPRTRV

!        Implicit arguments :  YOMMP parameters
!                              NPRGPNS,NPRGPEW,NPRTRW,NPRTRV,NPROC

!        --------------------
!     Method.
!     -------

!        PE allocation order is row oriented (e.g. NPRGPNS or NPRTRW = 4):

!                1  2  3  4
!                5  6  7  8
!                9 10 11 12
!               13 14 15 16
!                .  .  .  .

!     Externals.
!     ----------
!         NONE

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        David Dent *ECMWF*

!     Modifications.
!     --------------
!        Original : 98-08-19
!        Revision : 98-10-13 row ordering
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE TPM_DISTR       ,ONLY : LEQ_REGIONS, NPRGPEW, NPROC, NPRTRV
USE EQ_REGIONS_MOD  ,ONLY : N_REGIONS, N_REGIONS_NS
USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS
!

IMPLICIT NONE
INTEGER(KIND=JPIM),INTENT(IN)  :: KPE
INTEGER(KIND=JPIM),INTENT(OUT) :: KPRGPNS,KPRGPEW,KPRTRW,KPRTRV

INTEGER(KIND=JPIM) :: IPE,JA
!     ------------------------------------------------------------------

!*       1.    Check input argument for validity
!              ---------------------------------

IF(KPE <= 0.OR.KPE > NPROC) THEN
  WRITE(*,'(A,2I8)') ' PE2SET INVALID ARGUMENT ',KPE,NPROC
  CALL ABORT_TRANS(' PE2SET INVALID ARGUMENT ')

ELSE

!*       2.    Compute output parameters
!              -------------------------

  IF( LEQ_REGIONS )THEN
    KPRGPNS=1
    IPE=KPE
    DO JA=1,N_REGIONS_NS
      IF( IPE > N_REGIONS(JA) )THEN
        IPE=IPE-N_REGIONS(JA)
        KPRGPNS=KPRGPNS+1
        CYCLE
      ENDIF
      KPRGPEW=IPE
      EXIT
    ENDDO
  ELSE
    KPRGPEW=MOD(KPE-1,NPRGPEW)+1
    KPRGPNS=(KPE-1)/NPRGPEW+1
  ENDIF
  KPRTRV =MOD(KPE-1,NPRTRV)+1
  KPRTRW =(KPE-1)/NPRTRV+1

ENDIF

END SUBROUTINE PE2SET
END MODULE PE2SET_MOD
