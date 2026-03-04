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
!        Original : 98-08-19

!     Modifications.
!     --------------
!        Y.Tremolet: 02-03-13 use groups
!        R. El Khatib : 03-01-23 Case LMPOFF
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMMP0   , ONLY : NPRGPNS, NPRGPEW, NPRTRW, NPRTRV, LMPOFF, LEQ_REGIONS, N_REGIONS_NS, N_REGIONS
USE YOMLUN   , ONLY : NULERR
USE MPL_MODULE

IMPLICIT NONE
INTEGER(KIND=JPIM),INTENT(OUT)   :: KPE 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPRGPNS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPRGPEW 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPRTRW 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPRTRV 
INTEGER(KIND=JPIM) :: IPE,JA
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

!     ------------------------------------------------------------------

!*       1.    Choose from input parameters
!              ----------------------------

IF (LHOOK) CALL DR_HOOK('SET2PE',0,ZHOOK_HANDLE)

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
      CALL ABOR1(' SET2PE INVALID ARGUMENT ')
    ENDIF
  ENDIF

ELSE

!*       3.    Spectral space set values supplied
!              ----------------------------------

  IF(KPRTRW <= NPRTRW.AND.KPRTRV <= NPRTRV) THEN
    SELECT CASE (LMPOFF)
    CASE (.TRUE.)
      KPE=(KPRTRW-1)*NPRTRV + KPRTRV
    CASE (.FALSE.)
      KPE = MPL_CART_RANK(KPRTRW,KPRTRV)

!       Just checking for now...
      IPE=(KPRTRW-1)*NPRTRV + KPRTRV
      IF (IPE/=KPE) THEN
        WRITE(NULERR,*)'SET2PE kpe, ipe=',KPE,IPE
        CALL ABOR1('SET2PE wrong group values')
      ENDIF
    END SELECT

  ELSE
    WRITE(*,'(A,2I8)') ' SET2PE INVALID ARGUMENT ',KPRTRW,KPRTRV
    CALL ABOR1(' SET2PE INVALID ARGUMENT ')
  ENDIF

ENDIF

IF (LHOOK) CALL DR_HOOK('SET2PE',1,ZHOOK_HANDLE)
END SUBROUTINE SET2PE
