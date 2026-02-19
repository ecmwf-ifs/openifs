! (C) Copyright 2001- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE LEINVAD_MOD
CONTAINS
SUBROUTINE LEINVAD(KM,KMLOC,KFC,KIFC,KF_OUT_LT,KDGLU,PIA,PAOA1,PSOA1)

!**** *LEINVAD* - Inverse Legendre transform.

!     Purpose.
!     --------
!        Inverse Legendre tranform of all variables(kernel).

!**   Interface.
!     ----------
!        CALL LEINVAD(...)

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

!     Method.
!     -------

!     Externals.   MXMAOP - calls SGEMVX (matrix multiply)
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-02-01 From LEINVAD in IFS CY22R1
!        Modified ! 16/10/12 J.Hague : DR_HOOK round calls to DGEMM:
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB     ,JPRD
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE TPM_DIM         ,ONLY : R
USE TPM_GEOMETRY    ,ONLY : G
USE TPM_FIELDS      ,ONLY : F
!USE TPM_TRANS
USE TPM_DISTR       ,ONLY : D
!
USE TPM_FLT
USE BUTTERFLY_ALG_MOD

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN)    :: KM
INTEGER(KIND=JPIM), INTENT(IN)    :: KMLOC
INTEGER(KIND=JPIM), INTENT(IN)    :: KFC
INTEGER(KIND=JPIM), INTENT(IN)    :: KIFC
INTEGER(KIND=JPIM), INTENT(IN)    :: KDGLU
INTEGER(KIND=JPIM), INTENT(IN)    :: KF_OUT_LT
REAL(KIND=JPRB),    INTENT(OUT)   :: PIA(:,:)
REAL(KIND=JPRB),    INTENT(INOUT) :: PSOA1(:,:)
REAL(KIND=JPRB),    INTENT(INOUT) :: PAOA1(:,:)

!     LOCAL VARIABLES
INTEGER(KIND=JPIM) :: IA, ILA, ILS, IS, ISKIP, ISL, IOAD1, JK,JI
INTEGER(KIND=JPIM) :: IF,ITHRESHOLD
REAL(KIND=JPRB)    :: ZBA((R%NSMAX-KM+2)/2,KIFC), ZBS((R%NSMAX-KM+3)/2,KIFC), ZC(KDGLU,KIFC)
LOGICAL, PARAMETER :: LLDOUBLE = (JPRD == JPRB)
CHARACTER(LEN=1) :: CLX
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

!*       1.       PERFORM LEGENDRE TRANFORM.
!                 --------------------------

!*       1.1      PREPARATIONS.

CLX = 'S'
IF (LLDOUBLE) CLX = 'D'

IA  = 1+MOD(R%NSMAX-KM+2,2)
IS  = 1+MOD(R%NSMAX-KM+1,2)
ILA = (R%NSMAX-KM+2)/2
ILS = (R%NSMAX-KM+3)/2
ISL = MAX(R%NDGNH-G%NDGLU(KM)+1,1)
IOAD1 = 2*KF_OUT_LT

IF(KM == 0)THEN
  ISKIP = 2
ELSE
  ISKIP = 1
ENDIF

IF( KDGLU > 0 ) THEN

  ITHRESHOLD=S%ITHRESHOLD


! 1. +++++++++++++ anti-symmetric

 ! we need the transpose of C
  IF=0
  DO JK=1,KFC,ISKIP
    IF=IF+1
    DO JI=1,KDGLU
      ZC(JI,IF) = PAOA1(JK,ISL+JI-1)
    ENDDO
  ENDDO

  IF(ILA <= ITHRESHOLD .OR. .NOT.S%LUSEFLT) THEN
     IF (LHOOK) CALL DR_HOOK('LE_'//CLX//'GEMM_1',0,ZHOOK_HANDLE)
     IF(LLDOUBLE)THEN
        CALL DGEMM('T','N',ILA,KIFC,KDGLU,1.0_JPRB,S%FA(KMLOC)%RPNMA,KDGLU,&
             &ZC,KDGLU,0._JPRB,ZBA,ILA)
     ELSE
        CALL SGEMM('T','N',ILA,KIFC,KDGLU,1.0_JPRB,S%FA(KMLOC)%RPNMA,KDGLU,&
             &ZC,KDGLU,0._JPRB,ZBA,ILA)
     END IF
     IF (LHOOK) CALL DR_HOOK('LE_'//CLX//'GEMM_1',1,ZHOOK_HANDLE)

  ELSE
    
    CALL MULT_BUTM('T',S%FA(KMLOC)%YBUT_STRUCT_A,KIFC,ZC,ZBA)

  ENDIF

  IF=0
  DO JK=1,KFC,ISKIP
    IF=IF+1
    DO JI=1,ILA
      PIA(IA+1+(JI-1)*2,JK) = ZBA(JI,IF)
    ENDDO
  ENDDO

! 2. +++++++++++++ symmetric

 ! we need the transpose of C
  IF=0
  DO JK=1,KFC,ISKIP
    IF=IF+1
    DO JI=1,KDGLU
      ZC(JI,IF) = PSOA1(JK,ISL+JI-1)
    ENDDO
  ENDDO

  IF(ILS <= ITHRESHOLD .OR. .NOT.S%LUSEFLT ) THEN

    IF (LHOOK) CALL DR_HOOK('LE_'//CLX//'GEMM_2',0,ZHOOK_HANDLE)
    IF(LLDOUBLE)THEN
       CALL DGEMM('T','N',ILS,KIFC,KDGLU,1.0_JPRB,S%FA(KMLOC)%RPNMS,KDGLU,&
            &ZC,KDGLU,0._JPRB,ZBS,ILS)
    ELSE
       CALL SGEMM('T','N',ILS,KIFC,KDGLU,1.0_JPRB,S%FA(KMLOC)%RPNMS,KDGLU,&
            &ZC,KDGLU,0._JPRB,ZBS,ILS)
    END IF
    IF (LHOOK) CALL DR_HOOK('LE_'//CLX//'GEMM_2',1,ZHOOK_HANDLE)

  ELSE

    CALL MULT_BUTM('T',S%FA(KMLOC)%YBUT_STRUCT_S,KIFC,ZC,ZBS)

  ENDIF

  IF=0
  DO JK=1,KFC,ISKIP
    IF=IF+1
    DO JI=1,ILS
      PIA(IS+1+(JI-1)*2,JK) = ZBS(JI,IF)
    ENDDO
  ENDDO


ENDIF
!
!     ------------------------------------------------------------------


END SUBROUTINE LEINVAD
END MODULE LEINVAD_MOD
