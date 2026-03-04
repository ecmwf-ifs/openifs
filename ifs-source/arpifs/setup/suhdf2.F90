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

#ifdef RS6K
@PROCESS NOOPTIMIZE
#endif
!OPTION! -O nochg
!OCL  NOUNROLL,NOPREEX,NOEVAL
SUBROUTINE SUHDF2(YDDIM,YDDYN)

!**** *SUHDF2*   - Initialize horizontal diffusion for 2D models (shallow-water, vorticity models)

!     Purpose.
!     --------

!         COMPUTES HORIZONTAL DIFFUSION COEFFICIENTS BY USE OF SHAPE FUNCTION

!**   Interface.
!     ----------
!        *CALL* *SUHDF2

!        Explicit arguments :
!        --------------------

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ARPEGE/ALADIN DOCUMENTATION

!     Author.
!     -------
!      K. Yessad (CNRM/GMAP) from SUHDF
!      Original : Jan 2012

!     Modifications.
!     --------------
!     ------------------------------------------------------------------

USE YOMDIM   , ONLY : TDIM
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMLUN   , ONLY : NULOUT
USE YOMDYN   , ONLY : TDYN

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TDIM) , INTENT(IN) :: YDDIM
TYPE(TDYN)  ,INTENT(INOUT):: YDDYN

INTEGER(KIND=JPIM) :: JN
LOGICAL :: LLGRID
REAL(KIND=JPRB) :: ZEFOLD, ZEFOLV

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

! * FUNCTIONS AND FUNCTION ARGUMENTS:
REAL(KIND=JPRB) :: PDISPE
REAL(KIND=JPRB) :: PDISPVOR
INTEGER(KIND=JPIM) :: KMAX, KN
REAL(KIND=JPRB) :: PEXPDH

! * HORIZONTAL DIFFUSION SHAPE FUNCTIONS:
PDISPE(KN,KMAX,PEXPDH)=( MAX ( 0.0_JPRB ,&
 & (SQRT(REAL(KN*(KN+1),JPRB)/REAL(KMAX*(KMAX+1),JPRB))-YDDYN%FRANDH)&
 & /(1.0_JPRB-YDDYN%FRANDH) ) )**PEXPDH  
PDISPVOR(KN,KMAX,PEXPDH)=( MAX ( 0.0_JPRB ,&
 & (SQRT(MAX(0.0_JPRB,REAL(KN*(KN+1)-MIN(REAL(KMAX,JPRB),2.0_JPRB),JPRB) &
 & /REAL(KMAX*(KMAX+1)-MIN(REAL(KMAX,JPRB),2.0_JPRB),JPRB))) &
 & -YDDYN%FRANDH)&
 & /(1.0_JPRB-YDDYN%FRANDH) ) )**PEXPDH  

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUHDF2',0,ZHOOK_HANDLE)
ASSOCIATE(NDLON=>YDDIM%NDLON, NSMAX=>YDDIM%NSMAX, &
 & FRANDH=>YDDYN%FRANDH, HDIRDIV=>YDDYN%HDIRDIV, HDIRSP=>YDDYN%HDIRSP, &
 & HDIRVOR=>YDDYN%HDIRVOR, LRDISPE_EC=>YDDYN%LRDISPE_EC, RDIDIV=>YDDYN%RDIDIV, &
 & RDISP=>YDDYN%RDISP, RDIVOR=>YDDYN%RDIVOR, REXPDH=>YDDYN%REXPDH)
!     ------------------------------------------------------------------

!*       1.    MAIN HORIZONTAL DIFFUSION.
!              --------------------------

LLGRID=NSMAX > (NDLON+3)/3
WRITE(NULOUT,*) ' SUHDF2, LLGRID=',LLGRID,' NSMAX=',NSMAX,' NDLON=',NDLON

!*       1.1   Computation of RDI[X] for ECMWF type set-up:

IF (LRDISPE_EC) THEN

  WAVE_LOOP_EC : DO JN=0,NSMAX
    IF(HDIRVOR > 1.0_JPRB) THEN
      RDIVOR(1,JN)=PDISPE(JN,NSMAX,REXPDH)/HDIRVOR
    ELSE
      RDIVOR(1,JN)=0.0_JPRB
    ENDIF
    IF(HDIRDIV > 1.0_JPRB) THEN
      RDIDIV(1,JN)=PDISPE(JN,NSMAX,REXPDH)/HDIRDIV
    ELSE
      RDIDIV(1,JN)=0.0_JPRB
    ENDIF
    IF (HDIRSP > 1.0_JPRB) THEN
      RDISP(JN)=PDISPE(JN,NSMAX,REXPDH)/HDIRSP
    ELSE
      RDISP(JN)=0.0_JPRB
    ENDIF
  ENDDO WAVE_LOOP_EC

ENDIF ! LRDISPE_EC

!*       1.2   Computation of RDI[X] for METEO-FRANCE type set-up:

IF (.NOT.LRDISPE_EC) THEN

  WAVE_LOOP_MF : DO JN=0,NSMAX
    IF(HDIRVOR > 1.0_JPRB) THEN
      RDIVOR(1,JN)=PDISPVOR(JN,NSMAX,REXPDH)/HDIRVOR
    ELSE
      RDIVOR(1,JN)=0.0_JPRB
    ENDIF
    IF(HDIRDIV > 1.0_JPRB) THEN
      RDIDIV(1,JN)=PDISPE(JN,NSMAX,REXPDH)/HDIRDIV
    ELSE
      RDIDIV(1,JN)=0.0_JPRB
    ENDIF
    IF (HDIRSP > 1.0_JPRB) THEN
      RDISP(JN)=PDISPE(JN,NSMAX,REXPDH)/HDIRSP
    ELSE
      RDISP(JN)=0.0_JPRB
    ENDIF
  ENDDO WAVE_LOOP_MF

ENDIF ! .NOT.LRDISPE_EC

!*       1.3   Diagnostic of e-folding time:

IF(RDIVOR(1,NSMAX) > 1E-7_JPRB) THEN
  ZEFOLV=1.0_JPRB/RDIVOR(1,NSMAX)/3600._JPRB
ELSE
  ZEFOLV=0._JPRB
ENDIF
IF(RDIDIV(1,NSMAX) > 1E-7_JPRB) THEN
  ZEFOLD=1.0_JPRB/RDIDIV(1,NSMAX)/3600._JPRB
ELSE
  ZEFOLD=0._JPRB
ENDIF
WRITE(NULOUT,*) '  E-FOLDING TIME FOR VORTICITY=',ZEFOLV,' AND FOR DIVERGENCE=',ZEFOLD,' HOURS '

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUHDF2',1,ZHOOK_HANDLE)
END SUBROUTINE SUHDF2
