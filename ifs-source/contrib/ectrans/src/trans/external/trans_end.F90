! (C) Copyright 2000- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE TRANS_END(CDMODE)

!**** *TRANS_END* - Terminate transform package

!     Purpose.
!     --------
!     Terminate transform package. Release all allocated arrays.

!**   Interface.
!     ----------
!     CALL TRANS_END

!     Explicit arguments : None
!     --------------------

!     Method.
!     -------

!     Externals.  None
!     ----------

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-03-03
!          G. Radnoti: 19-03-2009: intermediate end of transf to allow to switch to mono-task transforms
!        R. El Khatib 09-Jul-2013 LENABLED

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

!ifndef INTERFACE

USE TPM_GEN         ,ONLY : MSETUP0, NCUR_RESOL, NMAX_RESOL, LENABLED,NDEF_RESOL
USE TPM_DIM         ,ONLY : R, DIM_RESOL
USE TPM_DISTR       ,ONLY : D, DISTR_RESOL, NPRCIDS
USE TPM_GEOMETRY    ,ONLY : G, GEOM_RESOL
USE TPM_FIELDS      ,ONLY : F, FIELDS_RESOL
USE TPM_FFT         ,ONLY : T, FFT_RESOL, TB, FFTB_RESOL
USE TPM_CTL         ,ONLY : C, CTL_RESOL
#ifdef WITH_FFTW
USE TPM_FFTW        ,ONLY : TW, FFTW_RESOL
#endif
USE TPM_FLT
USE TPM_TRANS       ,ONLY : FOUBUF, FOUBUF_IN

USE EQ_REGIONS_MOD  ,ONLY : N_REGIONS
USE SET_RESOL_MOD   ,ONLY : SET_RESOL
USE DEALLOC_RESOL_MOD   ,ONLY : DEALLOC_RESOL
!

IMPLICIT NONE
CHARACTER*5, OPTIONAL,  INTENT(IN) :: CDMODE
! Local variables
INTEGER(KIND=JPIM) :: JRES
CHARACTER*5 :: CLMODE
!     ------------------------------------------------------------------
CLMODE='FINAL'
IF (PRESENT(CDMODE)) CLMODE=CDMODE
IF (CLMODE == 'FINAL') THEN

  IF( ALLOCATED( LENABLED ) ) THEN
    DO JRES=1,NMAX_RESOL
      IF(LENABLED(JRES)) THEN
        CALL DEALLOC_RESOL(JRES)
      ENDIF
    ENDDO
    DEALLOCATE(LENABLED)
  ENDIF

  NULLIFY(R)
  IF( ALLOCATED(DIM_RESOL) ) DEALLOCATE(DIM_RESOL)

  NULLIFY(D)
  IF( ALLOCATED(DISTR_RESOL) ) DEALLOCATE(DISTR_RESOL)

  !TPM_FFT
  NULLIFY(T)
  IF( ALLOCATED(FFT_RESOL) ) DEALLOCATE(FFT_RESOL)
  NULLIFY(TB)
  IF( ALLOCATED(FFTB_RESOL) ) DEALLOCATE(FFTB_RESOL)

#ifdef WITH_FFTW
  !TPM_FFTW
  NULLIFY(TW)
  IF( ALLOCATED(FFTW_RESOL) ) DEALLOCATE(FFTW_RESOL)
#endif

  !TPM_FLT
  NULLIFY(S)
  IF( ALLOCATED(FLT_RESOL) ) DEALLOCATE(FLT_RESOL)

  !TPM_CTL
  NULLIFY(C)
  IF( ALLOCATED(CTL_RESOL) ) DEALLOCATE(CTL_RESOL)

  !TPM_FIELDS
  NULLIFY(F)
  IF( ALLOCATED(FIELDS_RESOL) ) DEALLOCATE(FIELDS_RESOL)


  !TPM_GEOMETRY
  NULLIFY(G)
  IF( ALLOCATED(GEOM_RESOL) ) DEALLOCATE(GEOM_RESOL)

  !TPM_TRANS
  IF(ALLOCATED(FOUBUF_IN)) DEALLOCATE(FOUBUF_IN)
  IF(ALLOCATED(FOUBUF)) DEALLOCATE(FOUBUF)

  MSETUP0 = 0
  NMAX_RESOL = 0
  NCUR_RESOL = 0
  NDEF_RESOL = 0
ENDIF
IF (CLMODE == 'FINAL' .OR. CLMODE == 'INTER') THEN
  !EQ_REGIONS
  IF( ASSOCIATED(N_REGIONS) ) DEALLOCATE(N_REGIONS)
  !TPM_DISTR
  IF( ALLOCATED(NPRCIDS) ) DEALLOCATE(NPRCIDS)
ENDIF

!     ------------------------------------------------------------------

!endif INTERFACE

END SUBROUTINE TRANS_END
