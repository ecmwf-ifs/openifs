! (C) Copyright 1988- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE LEDIRAD_MOD
CONTAINS
SUBROUTINE LEDIRAD(KM,KMLOC,KFC,KIFC,KDGLU,KLED2,PAIA,PSIA,POA1)

!**** *LEDIRAD* - Direct Legendre transform.

!     Purpose.
!     --------
!        Direct Legendre tranform of state variables.

!**   Interface.
!     ----------
!        CALL LEDIRAD(...)

!        Explicit arguments :  KM - zonal wavenumber
!        --------------------  KFC - number of field to transform
!                              PAIA - antisymmetric part of Fourier
!                              fields for zonal wavenumber KM
!                              PSIA - symmetric part of Fourier
!                              fields for zonal wavenumber KM
!                              POA1 -  spectral
!                              fields for zonal wavenumber KM
!                              PLEPO - Legendre polonomials

!        Implicit arguments :  None.
!        --------------------

!     Method.
!     -------

!     Externals.   MXMAOP - matrix multiply
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 88-01-28
!        Modified : 91-07-01 Philippe Courtier/Mats Hamrud - Rewrite
!                            for uv formulation
!        Modified : 93-03-19 D. Giard - NTMAX instead of NSMAX
!        Modified : 04/06/99 D.Salmond : change order of AIA and SIA
!        Modified ! 16/10/12 J.Hague : DR_HOOK round calls to DGEMM:
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB     ,JPRD
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE TPM_DIM         ,ONLY : R
USE TPM_GEOMETRY    ,ONLY : G
!USE TPM_TRANS
!
USE TPM_FLT
USE TPM_FIELDS
USE TPM_DISTR
USE BUTTERFLY_ALG_MOD

IMPLICIT NONE


!     DUMMY ARGUMENTS
INTEGER(KIND=JPIM), INTENT(IN)  :: KM
INTEGER(KIND=JPIM), INTENT(IN)  :: KMLOC
INTEGER(KIND=JPIM), INTENT(IN)  :: KFC
INTEGER(KIND=JPIM), INTENT(IN)  :: KIFC
INTEGER(KIND=JPIM), INTENT(IN)  :: KDGLU
INTEGER(KIND=JPIM), INTENT(IN)  :: KLED2

REAL(KIND=JPRB),    INTENT(OUT)  :: PSIA(:,:),   PAIA(:,:)
REAL(KIND=JPRB),    INTENT(IN)   :: POA1(:,:)

INTEGER(KIND=JPIM) :: IA, ILA, ILS, IS, ISKIP, ISL, J, JK,JGL,J1
INTEGER(KIND=JPIM) :: IF,ITHRESHOLD
REAL(KIND=JPRB)    :: ZB(KDGLU,KIFC), ZCA((R%NTMAX-KM+2)/2,KIFC), ZCS((R%NTMAX-KM+3)/2,KIFC)
LOGICAL, PARAMETER :: LLDOUBLE = (JPRD == JPRB)
CHARACTER(LEN=1) :: CLX
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

!*       1.       PERFORM LEGENDRE TRANFORM.
!                 --------------------------

!*       1.1      PREPARATIONS.

CLX = 'S'
IF (LLDOUBLE) CLX = 'D'

IA  = 1+MOD(R%NTMAX-KM+2,2)
IS  = 1+MOD(R%NTMAX-KM+1,2)
ILA = (R%NTMAX-KM+2)/2
ILS = (R%NTMAX-KM+3)/2
ISL = MAX(R%NDGNH-G%NDGLU(KM)+1,1)

IF(KM == 0)THEN
  ISKIP = 2
  DO JGL=ISL,R%NDGNH
    DO J1=2,KFC,2
      PSIA(J1,JGL)=0.0_JPRB
      PAIA(J1,JGL)=0.0_JPRB
    ENDDO
  ENDDO
ELSE
  ISKIP = 1
ENDIF


IF (KIFC > 0 .AND. KDGLU > 0 ) THEN

  ITHRESHOLD=S%ITHRESHOLD
 
!*       1. ANTISYMMETRIC PART.

  IF=0
  DO JK=1,KFC,ISKIP
    IF=IF+1
    DO J=1,ILA
      ZCA(J,IF) = POA1(IA+(J-1)*2,JK) 
    ENDDO
  ENDDO
  
  IF(ILA <= ITHRESHOLD .OR. .NOT.S%LUSEFLT) THEN
     IF (LHOOK) CALL DR_HOOK('LE_'//CLX//'GEMM_1',0,ZHOOK_HANDLE)
     IF(LLDOUBLE)THEN
        CALL DGEMM('N','N',KDGLU,KIFC,ILA,1.0_JPRB,S%FA(KMLOC)%RPNMA,KDGLU,&
             &ZCA,ILA,0._JPRB,ZB,KDGLU)
     ELSE
        CALL SGEMM('N','N',KDGLU,KIFC,ILA,1.0_JPRB,S%FA(KMLOC)%RPNMA,KDGLU,&
             &ZCA,ILA,0._JPRB,ZB,KDGLU)
     END IF
     IF (LHOOK) CALL DR_HOOK('LE_'//CLX//'GEMM_1',1,ZHOOK_HANDLE)

  ELSE

    CALL MULT_BUTM('N',S%FA(KMLOC)%YBUT_STRUCT_A,KIFC,ZCA,ZB)

  ENDIF

  IF=0
  DO JK=1,KFC,ISKIP
    IF=IF+1
    DO J=1,KDGLU
      PAIA(JK,ISL+J-1) = ZB(J,IF)*F%RW(ISL+J-1)
    ENDDO
  ENDDO

  
!*       1.3      SYMMETRIC PART.

  IF=0
  DO JK=1,KFC,ISKIP
    IF=IF+1
    DO J=1,ILS
      ZCS(J,IF) = POA1(IS+(J-1)*2,JK)
    ENDDO
  ENDDO
  
  
  IF(ILS <= ITHRESHOLD .OR. .NOT.S%LUSEFLT) THEN

    IF (LHOOK) CALL DR_HOOK('LE_'//CLX//'GEMM_2',0,ZHOOK_HANDLE)
    IF(LLDOUBLE)THEN
       CALL DGEMM('N','N',KDGLU,KIFC,ILS,1.0_JPRB,S%FA(KMLOC)%RPNMS,KDGLU,&
            &ZCS,ILS,0._JPRB,ZB,KDGLU)
    ELSE
       CALL SGEMM('N','N',KDGLU,KIFC,ILS,1.0_JPRB,S%FA(KMLOC)%RPNMS,KDGLU,&
            &ZCS,ILS,0._JPRB,ZB,KDGLU)

    END IF
    IF (LHOOK) CALL DR_HOOK('LE_'//CLX//'GEMM_2',1,ZHOOK_HANDLE)
    
  ELSE

    CALL MULT_BUTM('N',S%FA(KMLOC)%YBUT_STRUCT_S,KIFC,ZCS,ZB)
    
  ENDIF

  IF=0
  DO JK=1,KFC,ISKIP
    IF=IF+1
    DO J=1,KDGLU
      PSIA(JK,ISL+J-1) = ZB(J,IF)*F%RW(ISL+J-1)
    ENDDO
  ENDDO
  
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE LEDIRAD
END MODULE LEDIRAD_MOD
