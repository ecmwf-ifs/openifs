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

SUBROUTINE BALADSM(PNOR,KM,KSMAX,PVOR,PZ)

!**** *BALADSM*  - Selfadjoint linear balance equation

!     Purpose.
!     --------
!           vorticity to equivalent perturbation height conversion

!**   Interface.
!        *CALL* *BALADSM(...)

!     Explicit arguments :
!     --------------------
!      Input:
!        PNOR normalization factor of PZ (1 for Phi, 1/RG for height)
!        KM zonal wave number
!        KSMAX truncation n of spectral arrays
!        PVOR(2,KM:KSMAX) input vorticity fields
!        PZ(2,KM:KSMAX)  output geopotential

!     Implicit arguments :
!     --------------------
!        None

!     Method.
!     -------
!       Solve linear balance equation in spectral space
!       to convert vorticity into geopotential.
!       WARNING: the mean height is set to zero 
!          One should set it e.g. to 5560m for z500 field.
!          (In the adjoint, the mean vorticity gradient is always zero.)
!       Method. See Rochas,la Meteorologie, V-17-1971,p.32
!       solve lapin.div(f.grad.lapin.dzeta)=phi
!       using normalized spectral functions, expressions for
!       lapin,muYmn,dYmn/dmu,f in spectral coordinates,
!       spheric div-> 2romega.mu.dzeta+2romega(1-mu2)dpsi/dmu=lapphi
!       Nota Bene : the discretisation is such that the operator is
!        self adjoint and invertible (though it is ill-conditioned)

!     Externals.
!     ----------
!        None

!     Reference.
!     ----------
!        Rochas,la Meteorologie, V-17-1971,p.32 (in French).
!        METEO-FRANCE/GMAP internal documentation

!     Author.
!     -------
!      Original : 96-03-07 Philippe Courtier

!     Modifications.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMCST   , ONLY : ROMEGA   ,RA

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KM 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSMAX 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PNOR 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVOR(2,KM:KSMAX) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PZ(2,KM:KSMAX) 
INTEGER(KIND=JPIM) :: IN, IN1, IN2, JN

REAL(KIND=JPRB) :: ZA, ZB, ZCC
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     special case of km=ksmax: do nothing

IF (LHOOK) CALL DR_HOOK('BALADSM',0,ZHOOK_HANDLE)
PZ(:,:)=0.0_JPRB
IF(KM == KSMAX .AND. LHOOK) CALL DR_HOOK('BALADSM',1,ZHOOK_HANDLE)
IF(KM == KSMAX) RETURN

IN1=MAX(1,ABS(KM))
ZCC=-2.0_JPRB*ROMEGA*RA**2*PNOR

!     special case of m=0, imaginary part and mean set to 0         

IF(KM == 0)THEN
  DO JN=0,KSMAX
    PZ(2,JN)=0.0_JPRB
  ENDDO
  PZ(1,0)=0.0_JPRB
ENDIF

!        Last coefficient is not always computed (to ensure bijectivity)

IF(MOD(KSMAX-IN1+1,2) == 1)THEN
  IN2=KSMAX-1
  PZ(1,KSMAX)=0.0_JPRB
  PZ(1,KSMAX)=0.0_JPRB
ELSE
  IN2=KSMAX
ENDIF

!         first line of matrix

IN=IN1
ZB=SQRT(REAL((IN+1+KM)*(IN+1-KM),JPRB)/REAL((2*IN+3)*(2*IN+1),JPRB))*&
 & (1.0_JPRB-1.0_JPRB/REAL(IN+1,JPRB))  
PZ(1,IN)=ZCC/REAL(IN*(IN+1),JPRB)*ZB*PVOR(1,IN+1)
IF(KM /= 0)THEN
  PZ(2,IN)=ZCC/REAL(IN*(IN+1),JPRB)*ZB*PVOR(2,IN+1)
ENDIF

!         central lines of tridiagonal matrix

DO JN=IN1+1,IN2-1
  IN=JN
  ZA=SQRT(REAL((IN+KM)*(IN-KM),JPRB)/REAL((2*IN+1)*(2*IN-1),JPRB))*&
   & (1.0_JPRB+1.0_JPRB/REAL(IN,JPRB))  
  ZB=SQRT(REAL((IN+1+KM)*(IN+1-KM),JPRB)&
   & /REAL((2*IN+3)*(2*IN+1),JPRB))*(1.0_JPRB-1.0_JPRB/REAL(IN+1,JPRB))  
  PZ(1,IN)=ZCC/(IN*(IN+1))*(ZA*PVOR(1,IN-1)+ZB*PVOR(1,IN+1))
  IF(KM /= 0)THEN
    PZ(2,IN)=ZCC/(IN*(IN+1))*(ZA*PVOR(2,IN-1)+ZB*PVOR(2,IN+1))
  ENDIF
ENDDO

!         last line of matrix (drop high order term)

IN=IN2
ZA=SQRT(REAL((IN+KM)*(IN-KM),JPRB)/REAL((2*IN+1)*(2*IN-1),JPRB))*&
 & (1.0_JPRB+1.0_JPRB/REAL(IN,JPRB))  
PZ(1,IN)=ZCC/(IN*(IN+1))*ZA*PVOR(1,IN-1)
IF(KM == 0)THEN
  PZ(2,IN)=ZCC/(IN*(IN+1))*ZA*PVOR(2,IN-1)
ENDIF
IF (LHOOK) CALL DR_HOOK('BALADSM',1,ZHOOK_HANDLE)
END SUBROUTINE BALADSM
