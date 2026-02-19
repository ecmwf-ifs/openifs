! (C) Copyright 2000- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE SUFFT_MOD
CONTAINS
SUBROUTINE SUFFT

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE TPM_DIM         ,ONLY : R
USE TPM_GEN         ,ONLY : NOUT, NPRINTLEV
USE TPM_DISTR       ,ONLY : D, MYSETW
USE TPM_GEOMETRY    ,ONLY : G
USE TPM_FFT         ,ONLY : T, TB
#ifdef WITH_FFTW
USE TPM_FFTW        ,ONLY : TW, INIT_PLANS_FFTW
#endif
USE BLUESTEIN_MOD   ,ONLY : BLUESTEIN_INIT, FFTB_TYPE
!

IMPLICIT NONE

INTEGER(KIND=JPIM) :: JGL,IGLG, ILATS
LOGICAL :: LLP1,LLP2

!     ------------------------------------------------------------------

IF(.NOT.D%LGRIDONLY) THEN

  LLP1 = NPRINTLEV>0
  LLP2 = NPRINTLEV>1
  IF(LLP1) WRITE(NOUT,*) '=== ENTER ROUTINE SUFFT ==='

#ifdef WITH_FFTW
  IF(TW%LFFTW)THEN

    CALL INIT_PLANS_FFTW(R%NDLON)

  ELSE

    NULLIFY(TW%FFTW_PLANS)
#endif

    ALLOCATE(T%TRIGS(R%NDLON,D%NDGL_FS))
    IF(LLP2)WRITE(NOUT,9) 'T%TRIGS    ',SIZE(T%TRIGS),SHAPE(T%TRIGS)
    ALLOCATE(T%NFAX(19,D%NDGL_FS))
    IF(LLP2)WRITE(NOUT,9) 'T%NFAX     ',SIZE(T%NFAX),SHAPE(T%NFAX)
    ALLOCATE(T%LUSEFFT992(D%NDGL_FS))
    IF(LLP2)WRITE(NOUT,9) 'T%LUSEFFT992',SIZE(T%LUSEFFT992),SHAPE(T%LUSEFFT992)

    !
    ! create TRIGS and NFAX for latitude lengths supported by FFT992,
    ! that is just with factors 2, 3 or 5
    !

    T%LBLUESTEIN=.FALSE.
    ILATS=0
    DO JGL=1,D%NDGL_FS
      IGLG = D%NPTRLS(MYSETW)+JGL-1
      CALL SET99B(T%TRIGS(1,JGL),T%NFAX(1,JGL),G%NLOEN(IGLG),T%LUSEFFT992(JGL))
      IF( .NOT.T%LUSEFFT992(JGL) )THEN
        ILATS=ILATS+1
        T%LBLUESTEIN=.TRUE.
      ENDIF
    ENDDO
    
    !
    ! we only initialise for bluestein if there are latitude lengths 
    ! not supported by FFT992
    !

    IF( T%LBLUESTEIN )THEN
      TB%NDLON=R%NDLON
      TB%NLAT_COUNT=ILATS
      ILATS=0
      ALLOCATE(TB%NLATS(TB%NLAT_COUNT))
      DO JGL=1,D%NDGL_FS
        IF( .NOT.T%LUSEFFT992(JGL) )THEN
          ILATS=ILATS+1
          IGLG = D%NPTRLS(MYSETW)+JGL-1
          TB%NLATS(ILATS)=G%NLOEN(IGLG)
        ENDIF
      ENDDO
      CALL BLUESTEIN_INIT(TB)
    ENDIF
  
#ifdef WITH_FFTW

  ENDIF
#endif

ENDIF

!     ------------------------------------------------------------------

9 FORMAT(1X,'ARRAY ',A10,' ALLOCATED ',8I8)

END SUBROUTINE SUFFT
END MODULE SUFFT_MOD
