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

SUBROUTINE PE2SET(KPE,KPRGPNS,KPRGPEW,KPRTRW,KPRTRV)

!**** *PE2SET* - Convert from PE number to set numbers

!     Purpose.
!     --------
!        Convert from PE number to set numbers in both
!                  grid-point space and spectral space

!**   Interface.
!     ----------
!        *CALL* *PE2SET(...)

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
!      David Dent *ECMWF*
!      Original : 98-08-19

!     Modifications.
!     --------------
!      Y.Tremolet: 02-03-13 use groups
!      R. El Khatib : 03-01-23 Case LMPOFF
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      P.Marguinaud  07-Nov-2012 Make output arguments optional
!      T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMMP0   , ONLY : NPRGPEW, NPRTRV, NPROC, LMPOFF, LEQ_REGIONS, N_REGIONS_NS, N_REGIONS
USE YOMLUN   , ONLY : NULERR
USE MPL_MODULE

IMPLICIT NONE
INTEGER(KIND=JPIM),INTENT(IN)            :: KPE 
INTEGER(KIND=JPIM),INTENT(OUT), OPTIONAL :: KPRGPNS 
INTEGER(KIND=JPIM),INTENT(OUT), OPTIONAL :: KPRGPEW 
INTEGER(KIND=JPIM),INTENT(OUT), OPTIONAL :: KPRTRW 
INTEGER(KIND=JPIM),INTENT(OUT), OPTIONAL :: KPRTRV 

INTEGER(KIND=JPIM) :: IPRTRW,IPRTRV,IPE,JA

INTEGER(KIND=JPIM) :: JPRGPNS 
INTEGER(KIND=JPIM) :: JPRGPEW 
INTEGER(KIND=JPIM) :: JPRTRW 
INTEGER(KIND=JPIM) :: JPRTRV 

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

!     ------------------------------------------------------------------

!*       1.    Check input argument for validity 
!              ---------------------------------

IF (LHOOK) CALL DR_HOOK('PE2SET',0,ZHOOK_HANDLE)

IF(KPE <= 0.OR.KPE > NPROC) THEN

  WRITE(*,'(A,2I8)') ' PE2SET INVALID ARGUMENT ',KPE,NPROC
  CALL ABOR1(' PE2SET INVALID ARGUMENT ')

ELSE

!*       2.    Compute output parameters
!              -------------------------

  IF( LEQ_REGIONS )THEN
    JPRGPNS=1
    IPE=KPE
    DO JA=1,N_REGIONS_NS
      IF( IPE > N_REGIONS(JA) )THEN
        IPE=IPE-N_REGIONS(JA)
        JPRGPNS=JPRGPNS+1
        CYCLE
      ENDIF
      JPRGPEW=IPE
      EXIT
    ENDDO
  ELSE
    JPRGPEW=MOD(KPE-1,NPRGPEW)+1
    JPRGPNS=(KPE-1)/NPRGPEW+1
  ENDIF

  SELECT CASE (LMPOFF)
  CASE (.TRUE.)
    JPRTRW=(KPE-1)/NPRTRV+1
    JPRTRV=MOD(KPE-1,NPRTRV)+1
  CASE (.FALSE.)
    CALL MPL_CART_COORDS(KPE,JPRTRW,JPRTRV)

!     Just checking for now...
    IPRTRV =MOD(KPE-1,NPRTRV)+1
    IPRTRW =(KPE-1)/NPRTRV+1
    IF (IPRTRV/=JPRTRV .OR. IPRTRW/=JPRTRW) THEN
      WRITE(NULERR,*)'PE2SET IPRTRV,JPRTRV,IPRTRW,JPRTRW=', &
       & IPRTRV,JPRTRV,IPRTRW,JPRTRW  
      CALL ABOR1('PE2SET wrong group values')
    ENDIF
  END SELECT

ENDIF

IF (PRESENT (KPRGPNS)) KPRGPNS = JPRGPNS 
IF (PRESENT (KPRGPEW)) KPRGPEW = JPRGPEW 
IF (PRESENT (KPRTRW )) KPRTRW  = JPRTRW 
IF (PRESENT (KPRTRV )) KPRTRV  = JPRTRV 

IF (LHOOK) CALL DR_HOOK('PE2SET',1,ZHOOK_HANDLE)

END SUBROUTINE PE2SET
