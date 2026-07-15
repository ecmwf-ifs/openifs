! (C) Copyright 2000- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE FTDIR_MOD
CONTAINS
SUBROUTINE FTDIR(PREEL,KFIELDS,KGL)


!**** *FTDIR - Direct Fourier transform

!     Purpose. Routine for Grid-point to Fourier transform
!     --------

!**   Interface.
!     ----------
!        CALL FTDIR(..)

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
!        G. Radnoti 01-04-24 2D model (NLOEN=1)
!        D. Degrauwe  (Feb 2012): Alternative extension zone (E')
!        G. Mozdzynski (Oct 2014): support for FFTW transforms
!        G. Mozdzynski (Jun 2015): Support alternative FFTs to FFTW

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM, JPIB, JPRB

USE TPM_DISTR       ,ONLY : D, MYSETW
USE TPM_GEOMETRY    ,ONLY : G
USE TPM_FFT         ,ONLY : T, TB
USE BLUESTEIN_MOD   ,ONLY : BLUESTEIN_FFT
#ifdef WITH_FFTW
USE TPM_FFTW        ,ONLY : TW, EXEC_FFTW
#endif
USE TPM_DIM         ,ONLY : R
!

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)  :: KFIELDS,KGL
REAL(KIND=JPRB), POINTER, CONTIGUOUS, INTENT(INOUT) :: PREEL(:,:)

INTEGER(KIND=JPIM) :: IGLG,IST,ILEN,IJUMP,JJ,JF,IST1
INTEGER(KIND=JPIM) :: IOFF,IRLEN,ICLEN, ITYPE
LOGICAL :: LL_ALL=.FALSE. ! T=do kfields ffts in one batch, F=do kfields ffts one at a time

!     ------------------------------------------------------------------

ITYPE=-1
IJUMP= 1
IGLG = D%NPTRLS(MYSETW)+KGL-1
IST  = 2*(G%NMEN(IGLG)+1)+1
ILEN = G%NLOEN(IGLG)+R%NNOEXTZL+3-IST

IF (G%NLOEN(IGLG)>1) THEN 
  IOFF=D%NSTAGTF(KGL)+1
  IRLEN=G%NLOEN(IGLG)+R%NNOEXTZL
  ICLEN=(IRLEN/2+1)*2

#ifdef WITH_FFTW
  IF( .NOT. TW%LFFTW )THEN
#endif

    IF( T%LUSEFFT992(KGL) )THEN

      CALL FFT992(PREEL(1,IOFF),T%TRIGS(1,KGL),&
       &T%NFAX(1,KGL),KFIELDS,IJUMP,IRLEN,KFIELDS,ITYPE)

    ELSE

      CALL BLUESTEIN_FFT(TB,IRLEN,ITYPE,KFIELDS,PREEL(1:KFIELDS,IOFF:IOFF+ICLEN-1))
      DO JJ=1,ICLEN
        DO JF=1,KFIELDS
          PREEL(JF,IOFF+JJ-1)=PREEL(JF,IOFF+JJ-1)/REAL(IRLEN,JPRB)
        ENDDO
      ENDDO

    ENDIF
  
#ifdef WITH_FFTW
  ELSE

    CALL EXEC_FFTW(ITYPE,IRLEN,ICLEN,IOFF,KFIELDS,LL_ALL,PREEL)
  
  ENDIF
#endif

ENDIF

IST1=1
IF (G%NLOEN(IGLG)==1) IST1=0
DO JJ=IST1,ILEN
  PREEL(1:KFIELDS,IST+D%NSTAGTF(KGL)+JJ-1) = 0.0_JPRB
ENDDO

!     ------------------------------------------------------------------

END SUBROUTINE FTDIR
END MODULE FTDIR_MOD
