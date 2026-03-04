! (C) Copyright 1995- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE SUMPLATF_MOD
CONTAINS
SUBROUTINE SUMPLATF(KDGL,KPROCA,KMYSETA,&
                   &KULTPP,KPROCL,KPTRLS)

!**** *SUMPLATF * - Initialize fourier space distibution in N-S direction

!     Purpose.
!     --------


!**   Interface.
!     ----------
!        *CALL* *SUMPLATF *

!     Explicit arguments - input :
!     --------------------
!                          KDGL       -last  latitude
!                          KPROCA     -number of processors in A direction
!                          KMYSETA    -process number in A direction

!     Explicit arguments - output:
!     --------------------

!                          KULTPP     -number of latitudes in process
!                                      (in Fourier space)
!                          KPROCL     -process responsible for latitude
!                                      (in Fourier space)
!                          KPTRLS     -pointer to first global latitude
!                                      of process (in Fourier space)

!        Implicit arguments :
!        --------------------


!     Method.
!     -------
!        See documentation

!     Externals.   SUMPLATB and SUEMPLATB.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        MPP Group *ECMWF*

!     Modifications.
!     --------------
!        Original : 95-10-01
!        David Dent:97-06-02 parameters KFRSTLAT etc added
!        JF. Estrade:97-11-13 Adaptation to ALADIN case
!        J.Boutahar: 98-07-06  phasing with CY19
!        Modified 98-08-10 by K. YESSAD: removal of LRPOLE option + cleanings
!         (correct computation of extrapolar latitudes for KPROCL).
!        Modified 98-12-07 by K. YESSAD and C. FISCHER: cleaning.
!         - merge old sumplat.F and suemplat.F
!         - gather 'lelam' code and 'not lelam' code.
!         - clean (useless duplication of variables, non doctor features).
!         - remodularise according to lelam/not lelam
!           -> lelam features in new routine suemplatb.F,
!              not lelam features in new routine sumplatb.F
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE TPM_GEOMETRY    ,ONLY : G

USE SUMPLATB_MOD    ,ONLY : SUMPLATB
!

IMPLICIT NONE

!     * DUMMY:
INTEGER(KIND=JPIM),INTENT(IN)  :: KDGL
INTEGER(KIND=JPIM),INTENT(IN)  :: KPROCA
INTEGER(KIND=JPIM),INTENT(IN)  :: KMYSETA
INTEGER(KIND=JPIM),INTENT(OUT) :: KULTPP(:)
INTEGER(KIND=JPIM),INTENT(OUT) :: KPROCL(:)
INTEGER(KIND=JPIM),INTENT(OUT) :: KPTRLS(:)

!     * LOCAL:
INTEGER(KIND=JPIM) :: INDIC(KPROCA),ILAST(KPROCA)

!     LOCAL INTEGER SCALARS
INTEGER(KIND=JPIM) :: IA, ILAT, ISTART, IMEDIAP,IRESTM, JA,  JLTLOC

LOGICAL :: LLSPLIT,LLFOURIER

!      -----------------------------------------------------------------

!*       1.    CODE DEPENDING ON 'LELAM': COMPUTATION OF
!              KMEDIAP, KRESTM, INDIC, ILAST.
!              -----------------------------------------

LLSPLIT = .FALSE.
LLFOURIER = .TRUE.

CALL SUMPLATB(1,KDGL,KPROCA,G%NLOEN,LLSPLIT,LLFOURIER,&
 &IMEDIAP,IRESTM,INDIC,ILAST)

!      -----------------------------------------------------------------

!*       2.    CODE NOT DEPENDING ON 'LELAM':
!              ------------------------------



!     * Definitions related to distribution of latitudes along sets
!       ------------ in fourier-space -----------------------------
ISTART = 0
KULTPP(1) = ILAST(1)
DO JA=1,KPROCA
  IF(JA > 1) THEN
    IF(ILAST(JA) /= 0) THEN
      KULTPP(JA) = ILAST(JA)-ILAST(JA-1)
    ELSE
      KULTPP(JA) = 0
    ENDIF
  ENDIF
  DO JLTLOC=1,KULTPP(JA)
    ILAT = ISTART + JLTLOC
    KPROCL(ILAT) = JA
  ENDDO
  ISTART = ISTART + KULTPP(JA)
ENDDO

!     * Computes KPTRLS.

IA = KPROCL(1)
KPTRLS(IA) = 1
DO JA=IA+1,KPROCA
  KPTRLS(JA) = KPTRLS(JA-1) + KULTPP(JA-1)
ENDDO

END SUBROUTINE SUMPLATF
END MODULE SUMPLATF_MOD
