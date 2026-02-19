! (C) Copyright 2000- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE FTINVAD_MOD
CONTAINS
SUBROUTINE FTINVAD(PREEL,KFIELDS,KGL)


!**** *FTINVAD - Inverse Fourier transform - adjoint

!     Purpose. Routine for Fourier to Grid-point transform
!     --------

!**   Interface.
!     ----------
!        CALL FTINVAD(..)

!        Explicit arguments :  PREEL   - Fourier/grid-point array
!        --------------------  KFIELDS - number of fields

!     Method.
!     -------

!     Externals.  FFT992 - FFT routine
!     ----------
!

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-03-03
!        D. Degrauwe  (Feb 2012): Alternative extension zone (E')
!        G. Mozdzynski (Oct 2014): support for FFTW transforms
!        G. Mozdzynski (Jun 2015): Support alternative FFTs to FFTW

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM, JPIB, JPRB

USE TPM_DISTR       ,ONLY : D, MYSETW
USE TPM_DIM         ,ONLY : R
USE TPM_GEOMETRY    ,ONLY : G
USE TPM_FFT         ,ONLY : T, TB
USE BLUESTEIN_MOD   ,ONLY : BLUESTEIN_FFT
#ifdef WITH_FFTW
USE TPM_FFTW        ,ONLY : TW, EXEC_FFTW
#endif
USE TPM_DIM         ,ONLY : R
!

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KFIELDS,KGL
REAL(KIND=JPRB), INTENT(OUT)  :: PREEL(:,:)

INTEGER(KIND=JPIM) :: IGLG,IST,ILEN,IJUMP,JJ,JF,ILOEN
INTEGER(KIND=JPIM) :: IOFF,IRLEN,ICLEN,ITYPE
LOGICAL :: LL_ALL=.FALSE. ! T=do kfields ffts in one batch, F=do kfields ffts one at a time

!     ------------------------------------------------------------------

ITYPE =-1
IJUMP = 1
IGLG  = D%NPTRLS(MYSETW)+KGL-1
ILOEN = G%NLOEN(IGLG)+R%NNOEXTZL
IST   = 2*(G%NMEN(IGLG)+1)+1
ILEN  = ILOEN+3-IST
IOFF  = D%NSTAGTF(KGL)+1
IRLEN = ILOEN
ICLEN = (IRLEN/2+1)*2

  ! Change of metric (not in forward routine)
DO JJ=1,ILOEN
  DO JF=1,KFIELDS
    PREEL(JF,IOFF-1+JJ) = PREEL(JF,IOFF-1+JJ)*ILOEN
  ENDDO
ENDDO
  
#ifdef WITH_FFTW
IF( .NOT. TW%LFFTW )THEN
#endif

  IF( T%LUSEFFT992(KGL) )THEN

    CALL FFT992(PREEL(1,IOFF),T%TRIGS(1,KGL),&
     &T%NFAX(1,KGL),KFIELDS,IJUMP,ILOEN,KFIELDS,ITYPE)

  ELSE

    CALL BLUESTEIN_FFT(TB,ILOEN,ITYPE,KFIELDS,PREEL(1:KFIELDS,IOFF:IOFF+ICLEN-1))
    DO JJ=1,ICLEN
      DO JF=1,KFIELDS
        PREEL(JF,IOFF-1+JJ)=PREEL(JF,IOFF-1+JJ)/REAL(ILOEN,JPRB)
      ENDDO
    ENDDO

  ENDIF

#ifdef WITH_FFTW
ELSE

  CALL EXEC_FFTW(ITYPE,IRLEN,ICLEN,IOFF,KFIELDS,LL_ALL,PREEL)

ENDIF
#endif

DO JJ=1,ILEN
  DO JF=1,KFIELDS
    PREEL(JF,IST+IOFF-1+JJ-1) = 0.0_JPRB
  ENDDO
ENDDO

!     ------------------------------------------------------------------

END SUBROUTINE FTINVAD
END MODULE FTINVAD_MOD
