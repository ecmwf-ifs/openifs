! (C) Copyright 1998- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE SET2PE_MOD
CONTAINS
SUBROUTINE SET2PE(KPE,KPRGPNS,KPRGPEW,KPRTRW,KPRTRV)


!**** *SET2PE* - Convert from set numbers to PE number

!     Purpose.
!     --------
!        Convert from set numbers in either grid-point space or spectral space
!        to PE number

!**   Interface.
!     ----------
!        *CALL* *SET2PE(KPRGPNS,KPRGPEW,KPRTRW,KPRTRV,KPE)

!        Explicit arguments :
!        --------------------

!                  input :  KPRGPNS - integer A set number in grid space
!                                     in the range 1 .. NPRGPNS
!                           KPRGPEW - integer B set number in grid space
!                                     in the range 1 .. NPRGPEW
!                           KPRTRW  - integer A set number in spectral space
!                                     in the range 1 .. NPRTRW
!                           KPRTRV  - integer B set number in spectral space
!                                     in the range 1 .. NPRTRV
!                  output:  KPE     - integer processor number
!                                     in the range 1 .. NPROC

!                  Normally, one pair of input set numbers will be set to zero
!                  SET2PE will compute KPE from the first pair if they are valid numbers.
!                  else from the other pair,

!        Implicit arguments :  YOMMP parameters
!                              NPRGPNS,NPRGPEW,NPRTRW,NPRTRV,NPROC

!        --------------------
!     Method.
!     -------

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
!     ------------------------------------------------------------------


USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE TPM_DISTR       ,ONLY : LEQ_REGIONS, NPRGPEW, NPRGPNS, NPRTRV, NPRTRW
USE EQ_REGIONS_MOD  ,ONLY : N_REGIONS, N_REGIONS_NS
USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS
!

IMPLICIT NONE
INTEGER(KIND=JPIM),INTENT(IN) :: KPRGPNS,KPRGPEW,KPRTRW,KPRTRV
INTEGER(KIND=JPIM),INTENT(OUT)  :: KPE

INTEGER(KIND=JPIM) :: IPE,JA
!     ------------------------------------------------------------------

!*       1.    Choose from input parameters
!              ----------------------------

IF(KPRGPNS > 0.AND.KPRGPEW > 0) THEN

  IF( LEQ_REGIONS )THEN
    IF( KPRGPNS > N_REGIONS_NS )THEN
      WRITE(*,'(A,2I8)') ' SET2PE INVALID ARGUMENT ',KPRGPNS,N_REGIONS_NS
      CALL ABOR1(' SET2PE INVALID ARGUMENT ')
    ENDIF
    IF( KPRGPEW > N_REGIONS(KPRGPNS) )THEN
      WRITE(*,'(A,2I8)') ' SET2PE INVALID ARGUMENT ',KPRGPEW,N_REGIONS(KPRGPNS)
      CALL ABOR1(' SET2PE INVALID ARGUMENT ')
    ENDIF
    KPE=0
    DO JA=1,KPRGPNS-1
      KPE=KPE+N_REGIONS(JA)
    ENDDO
    KPE=KPE+KPRGPEW
  ELSE
    IF(KPRGPNS <= NPRGPNS.AND.KPRGPEW <= NPRGPEW) THEN

!*       2.    Grid-space set values supplied
!              ------------------------------

      KPE=(KPRGPNS-1)*NPRGPEW + KPRGPEW
    ELSE
      WRITE(*,'(A,2I8)') ' SET2PE INVALID ARGUMENT ',KPRGPNS,KPRGPEW
      CALL ABORT_TRANS(' SET2PE INVALID ARGUMENT ')
    ENDIF
  ENDIF

ELSE

!*       3.    Spectral space set values supplied
!              ----------------------------------

  IF(KPRTRW <= NPRTRW.AND.KPRTRV <= NPRTRV) THEN
    KPE=(KPRTRW-1)*NPRTRV + KPRTRV
  ELSE
    WRITE(*,'(A,2I8)') ' SET2PE INVALID ARGUMENT ',KPRTRW,KPRTRV
    CALL ABORT_TRANS(' SET2PE INVALID ARGUMENT ')
  ENDIF

ENDIF

END SUBROUTINE SET2PE
END MODULE SET2PE_MOD
