! (C) Copyright 2000- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE LEINV_MOD
CONTAINS
SUBROUTINE LEINV(KM,KMLOC,KFC,KIFC,KF_OUT_LT,KSL,KDGLU,PIA,PAOA1,PSOA1)

!**** *LEINV* - Inverse Legendre transform.

!     Purpose.
!     --------
!        Inverse Legendre tranform of all variables(kernel).

!**   Interface.
!     ----------
!        CALL LEINV(...)

!        Explicit arguments :  KM - zonal wavenumber (input-c)
!        --------------------  KFC - number of fields to tranform (input-c)
!                              PIA - spectral fields
!                              for zonal wavenumber KM (input)
!                              PAOA1 - antisymmetric part of Fourier
!                              fields for zonal wavenumber KM (output)
!                              PSOA1 - symmetric part of Fourier
!                              fields for zonal wavenumber KM (output)

!        Implicit arguments :  None.
!        --------------------

!     Method.    use butterfly or dgemm
!     -------

!     Externals.   
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      Nils Wedi + Mats Hamrud + George Modzynski
!
!     Modifications.
!     --------------
!        J.Hague : Oct 2012 DR_HOOK round calls to DGEMM:
!      F. Vana  05-Mar-2015  Support for single precision
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPRD, JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE TPM_DIM         ,ONLY : R
USE TPM_FLT
USE TPM_GEN ! Fpr nout
USE BUTTERFLY_ALG_MOD

use, intrinsic :: ieee_exceptions

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN)  :: KM
INTEGER(KIND=JPIM), INTENT(IN)  :: KMLOC
INTEGER(KIND=JPIM), INTENT(IN)  :: KFC
INTEGER(KIND=JPIM), INTENT(IN)  :: KIFC
INTEGER(KIND=JPIM), INTENT(IN)  :: KDGLU
INTEGER(KIND=JPIM), INTENT(IN)  :: KSL
INTEGER(KIND=JPIM), INTENT(IN)  :: KF_OUT_LT
REAL(KIND=JPRB),    INTENT(IN)  :: PIA(:,:)
REAL(KIND=JPRB),    INTENT(OUT) :: PSOA1(:,:)
REAL(KIND=JPRB),    INTENT(OUT) :: PAOA1(:,:)

!     LOCAL 
INTEGER(KIND=JPIM) :: IA, ILA, ILS, IS, ISKIP, ISL, J1, IF, JGL,JK, J,JI, IEND
INTEGER(KIND=JPIM) :: ITHRESHOLD
REAL(KIND=JPRB)    :: ZBA((R%NSMAX-KM+2)/2,KIFC), ZBS((R%NSMAX-KM+3)/2,KIFC), ZC(KDGLU,KIFC)
LOGICAL :: LL_HALT_INVALID
#ifdef WITH_IEEE_HALT
LOGICAL, PARAMETER :: LL_IEEE_HALT = .TRUE.
#else
LOGICAL, PARAMETER :: LL_IEEE_HALT = .FALSE.
#endif
LOGICAL, PARAMETER :: LLDOUBLE = (JPRB == JPRD)
CHARACTER(LEN=1) :: CLX
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

!*       1.       PERFORM LEGENDRE TRANFORM.
!                 --------------------------

!*       1.1      PREPARATIONS.

CLX = 'S'
IF (LLDOUBLE) CLX = 'D'

!ISL = MAX(R%NDGNH-G%NDGLU(KM)+1,1)
ISL = KSL
IEND = KSL + KDGLU - 1

IA  = 1+MOD(R%NSMAX-KM+2,2)
IS  = 1+MOD(R%NSMAX-KM+1,2)
ILA = (R%NSMAX-KM+2)/2
ILS = (R%NSMAX-KM+3)/2

IF(KM == 0)THEN
  ISKIP = 2
  DO J1=2,KFC,2
    DO JGL=ISL,IEND
      PSOA1(J1,JGL) = 0.0_JPRB
      PAOA1(J1,JGL) = 0.0_JPRB
    ENDDO
  ENDDO
ELSE
  ISKIP = 1
ENDIF

IF( KDGLU > 0 ) THEN

  ITHRESHOLD=S%ITHRESHOLD

  ! 1. +++++++++++++ anti-symmetric

  IF=0
  DO JK=1,KFC,ISKIP
    IF=IF+1
    DO J=1,ILA
      ZBA(J,IF)=PIA(IA+1+(J-1)*2,JK)
    ENDDO
  ENDDO
  
  IF(ILA <= ITHRESHOLD .OR. .NOT.S%LUSEFLT) THEN

    IF (LHOOK) CALL DR_HOOK('LEINV_'//CLX//'GEMM_1',0,ZHOOK_HANDLE)
    IF (LLDOUBLE) THEN
      CALL DGEMM('N','N',KDGLU,KIFC,ILA,1.0_JPRB,S%FA(KMLOC)%RPNMA,KDGLU,&
       &ZBA,ILA,0._JPRB,ZC,KDGLU)
    ELSE
       IF (LL_IEEE_HALT) THEN
          call ieee_get_halting_mode(ieee_invalid,LL_HALT_INVALID)
          if (LL_HALT_INVALID) call ieee_set_halting_mode(ieee_invalid,.false.)
       ENDIF
       CALL SGEMM('N','N',KDGLU,KIFC,ILA,1.0_JPRB,S%FA(KMLOC)%RPNMA,KDGLU,&
            &ZBA,ILA,0._JPRB,ZC,KDGLU)
       if (LL_IEEE_HALT .and. LL_HALT_INVALID) call ieee_set_halting_mode(ieee_invalid,.true.)
    ENDIF
    IF (LHOOK) CALL DR_HOOK('LEINV_'//CLX//'GEMM_1',1,ZHOOK_HANDLE)
  
  ELSE

    IF (LHOOK) CALL DR_HOOK('LEINV_'//CLX//'BUTM_1',0,ZHOOK_HANDLE)
    CALL MULT_BUTM('N',S%FA(KMLOC)%YBUT_STRUCT_A,KIFC,ZBA,ZC)
    IF (LHOOK) CALL DR_HOOK('LEINV_'//CLX//'BUTM_1',1,ZHOOK_HANDLE)
    
  ENDIF

  ! we need the transpose of C
  IF=0
  DO JK=1,KFC,ISKIP
    IF=IF+1
    DO JI=1,KDGLU
      PAOA1(JK,ISL+JI-1) = ZC(JI,IF)
    ENDDO
  ENDDO

  ! 2. +++++++++++++ symmetric

  IF=0
  DO JK=1,KFC,ISKIP
    IF=IF+1
    DO J=1,ILS
      ZBS(J,IF)=PIA(IS+1+(J-1)*2,JK)
    ENDDO
  ENDDO
  
  IF(ILS <= ITHRESHOLD .OR. .NOT.S%LUSEFLT ) THEN

    IF (LHOOK) CALL DR_HOOK('LEINV_'//CLX//'GEMM_2',0,ZHOOK_HANDLE)
    IF (LLDOUBLE) THEN
       CALL DGEMM('N','N',KDGLU,KIFC,ILS,1.0_JPRB,S%FA(KMLOC)%RPNMS,KDGLU,&
            &ZBS,ILS,0._JPRB,ZC,KDGLU)
    ELSE
       IF (LL_IEEE_HALT) THEN
          call ieee_get_halting_mode(ieee_invalid,LL_HALT_INVALID)
          if (LL_HALT_INVALID) call ieee_set_halting_mode(ieee_invalid,.false.)
       ENDIF
       CALL SGEMM('N','N',KDGLU,KIFC,ILS,1.0_JPRB,S%FA(KMLOC)%RPNMS,KDGLU,&
            &ZBS,ILS,0._JPRB,ZC,KDGLU)
       if (LL_IEEE_HALT .and. LL_HALT_INVALID) call ieee_set_halting_mode(ieee_invalid,.true.)
    ENDIF
    IF (LHOOK) CALL DR_HOOK('LEINV_'//CLX//'GEMM_2',1,ZHOOK_HANDLE)
    
  ELSE

    IF (LHOOK) CALL DR_HOOK('LEINV_'//CLX//'BUTM_2',0,ZHOOK_HANDLE)
    CALL MULT_BUTM('N',S%FA(KMLOC)%YBUT_STRUCT_S,KIFC,ZBS,ZC)
    IF (LHOOK) CALL DR_HOOK('LEINV_'//CLX//'BUTM_2',1,ZHOOK_HANDLE)

  ENDIF

  ! we need the transpose of C 
  IF=0
  DO JK=1,KFC,ISKIP
    IF=IF+1
    DO JI=1,KDGLU
      PSOA1(JK,ISL+JI-1) = ZC(JI,IF)
    ENDDO
  ENDDO
  
ENDIF
!     ------------------------------------------------------------------

END SUBROUTINE LEINV
END MODULE LEINV_MOD
