! (C) Copyright 2000- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE TRANS_PNM(KRESOL,KM,PRPNM,LDTRANSPOSE,LDCHEAP)

!**** *TRANS_PNM* - Compute Legendre polynomials for a given wavenember

!     Purpose.
!     --------
!     Interface routine for computing Legendre polynomials for a given wavenember

!**   Interface.
!     ----------
!     CALL TRANS_PNM(...)

!     Explicit arguments : All arguments are optional.
!     --------------------
!     KRESOL   - resolution tag for which info is required ,default is the
!                first defined resulution (input)
!     KM       - wave number
!     PRPNM    - Legendre polynomials
!     LDTRANSPOSE - Legendre polynomials array is transposed
!     LDCHEAP   - cheapest but less accurate computation

!     Method.
!     -------

!     Externals.  SET_RESOL - set resolution
!     ----------

!     Author.
!     -------
!        R. El Khatib *METEO-FRANCE*

!     Modifications.
!     --------------
!        Original : 22-Jan-2016 from G. Mozdzynski's getpnm

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPRD, JPIM, JPRB

!ifndef INTERFACE

USE TPM_DIM         ,ONLY : R
USE TPM_DISTR       ,ONLY : D
USE TPM_GEOMETRY    ,ONLY : G
USE TPM_FIELDS      ,ONLY : F
USE TPM_FLT         ,ONLY : S

USE SET_RESOL_MOD   ,ONLY : SET_RESOL
USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS
USE TPM_POL
USE SUPOLF_MOD

!endif INTERFACE

IMPLICIT NONE

INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KRESOL
INTEGER(KIND=JPIM) ,INTENT(IN)  :: KM
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PRPNM(:,:)
LOGICAL, OPTIONAL, INTENT(IN) :: LDTRANSPOSE
LOGICAL, OPTIONAL, INTENT(IN) :: LDCHEAP

!ifndef INTERFACE

INTEGER(KIND=JPIM) :: IU1, IU2, IMAXN, INMAX, ICHEAP_SYM, ICHEAP_ANTISYM
INTEGER(KIND=JPIM) :: IC, JN, JMLOC, JGL, JI
INTEGER(KIND=JPIM) :: IA, IS, IDGLU, ILA, ILS, ISL
REAL(KIND=JPRD), ALLOCATABLE :: ZLPOL(:)
LOGICAL :: LLTRANSPOSE, LLCHEAP
!     ------------------------------------------------------------------

! Set current resolution
IF (PRESENT(KRESOL)) THEN
  CALL SET_RESOL(KRESOL)
ENDIF

IF (PRESENT(LDTRANSPOSE)) THEN
  LLTRANSPOSE=LDTRANSPOSE
ELSE
  LLTRANSPOSE=.FALSE.
ENDIF

IF (PRESENT(LDCHEAP)) THEN
  LLCHEAP=LDCHEAP
ELSE
  LLCHEAP=.FALSE.
ENDIF
IF (LLCHEAP) THEN
  ICHEAP_SYM    =2
  ICHEAP_ANTISYM=3
ELSE
  ICHEAP_SYM    =1
  ICHEAP_ANTISYM=1
ENDIF

IF (PRESENT(PRPNM)) THEN
  IF(D%LGRIDONLY) THEN
    CALL ABORT_TRANS('TRANS_PNM: PRPNM REQUIRED BUT LGRIDONLY=T')
  ENDIF
ENDIF

IU1 = UBOUND(PRPNM,1)
IU2 = UBOUND(PRPNM,2)

IF (LLTRANSPOSE) THEN

  IF(IU2 < R%NLEI3) THEN
    CALL ABORT_TRANS('TRANS_PNM : FIRST DIM. OF PRPNM TOO SMALL')
  ENDIF
  IF(IU1 < R%NTMAX-KM+3) THEN
    CALL ABORT_TRANS('TRANS_PNM : SECOND DIM. OF PRPNM TOO SMALL')
  ENDIF

  IF (IU2 >= R%NLEI3) THEN
    PRPNM(:,R%NLEI3) = 0.0_JPRB
  ENDIF

ELSE

  IF(IU1 < R%NLEI3) THEN
    CALL ABORT_TRANS('TRANS_PNM : FIRST DIM. OF PRPNM TOO SMALL')
  ENDIF
  IF(IU2 < R%NTMAX-KM+3) THEN
    CALL ABORT_TRANS('TRANS_PNM : SECOND DIM. OF PRPNM TOO SMALL')
  ENDIF

  IF (IU1 >= R%NLEI3) THEN
    PRPNM(R%NLEI3,:) = 0.0_JPRB
  ENDIF

ENDIF

ILA = (R%NTMAX-KM+2)/2
ILS = (R%NTMAX-KM+3)/2

CALL INI_POL(R%NTMAX+2,LDFAST=.TRUE.)

IMAXN=R%NTMAX+1

IA  = 1+MOD(R%NTMAX-KM+2,2)
IS  = 1+MOD(R%NTMAX-KM+1,2)

ISL = MAX(R%NDGNH-G%NDGLU(KM)+1,1)
IF (S%LSOUTHPNM) THEN
  IDGLU = 2*MIN(R%NDGNH,G%NDGLU(KM))
ELSE
  IDGLU = MIN(R%NDGNH,G%NDGLU(KM))
ENDIF

IF(MOD(IMAXN-KM,2) == 0) THEN
  INMAX=IMAXN+1
ELSE
  INMAX=IMAXN
ENDIF

ALLOCATE(ZLPOL(0:R%NTMAX+2))

!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JGL,ZLPOL,JI,JN)
DO JGL=1,IDGLU
  CALL SUPOLF(KM,INMAX,REAL (F%RMU(ISL+JGL-1), JPRD),ZLPOL(0:INMAX),KCHEAP=ICHEAP_ANTISYM)
  IF (LLTRANSPOSE) THEN
    DO JI=1,ILA
      PRPNM(IA+(JI-1)*2,ISL+JGL-1) = ZLPOL(KM+2*(ILA-JI)+1)
    ENDDO
  ELSE
    DO JI=1,ILA
      PRPNM(ISL+JGL-1,IA+(JI-1)*2) = ZLPOL(KM+2*(ILA-JI)+1)
    ENDDO
  ENDIF
  CALL SUPOLF(KM,INMAX,REAL (F%RMU(ISL+JGL-1), JPRD),ZLPOL(0:INMAX),KCHEAP=ICHEAP_SYM)
  IF (LLTRANSPOSE) THEN
    DO JI=1,ILS
      PRPNM(IS+(JI-1)*2,ISL+JGL-1) = ZLPOL(KM+2*(ILS-JI))
    ENDDO
  ELSE
    DO JI=1,ILS
      PRPNM(ISL+JGL-1,IS+(JI-1)*2) = ZLPOL(KM+2*(ILS-JI))
    ENDDO
  ENDIF
ENDDO
!$OMP END PARALLEL DO

CALL END_POL

!     ------------------------------------------------------------------

!endif INTERFACE

END SUBROUTINE TRANS_PNM
