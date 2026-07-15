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

#ifdef NECSX
!option! -O extendreorder
#endif
SUBROUTINE LAITLI(KPROMA,KPROMB,KSTART,KPROF,KFLEV,KFLDN,KFLDX,PDLAT,PDLO,KL0,PDVER,&
 & PXSL,PXF)

!**** *LAITLI  -  semi-LAgrangian scheme:
!                 Trilinear interpolations for one variable.

!     Purpose.
!     --------
!       Performs trilinear interpolations for one variable.

!**   Interface.
!     ----------
!        *CALL* *LAITLI(KPROMA,KPROMB,KSTART,KPROF,KFLEV
!                      ,KFLDN,KFLDX
!                      ,PDLAT,PDLO,KL0,PDVER
!                      ,PXSL,PXF)

!        Explicit arguments :
!        --------------------

!        INPUT:
!          KPROMA  - horizontal dimension for grid-point quantities.
!          KPROMB  - horizontal dimension for interpolation point
!                    quantities.
!          KSTART  - first element of arrays where
!                    computations are performed.
!          KPROF   - depth of work.
!          KFLEV   - vertical dimension.
!          KFLDN   - number of the first field.
!          KFLDX   - number of the last field.
!          PDLAT   - weight (distance) for horizontal linear interpolation
!                    on a same latitude.
!          PDLO    - weights (distances) for horizontal linear interpolation
!                    on a same longitude.
!          KL0     - indices of the four western points
!                    of the 16 points interpolation grid.
!          PDVER   - weights (distances) for vertical linear interpolation
!                    on a same vertical.
!          PXSL    - semi-lagrangian variable.

!        OUTPUT:
!          PXF     - interpolated variable.

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!        No external.
!        Called by LARCIN.

!     Reference.
!     ----------

!     Author.
!     -------
!        K. YESSAD, after the subroutine LAGINL3
!        written by Maurice IMBARD, Alain CRAPLET and Michel ROCHAS
!        METEO-FRANCE, CNRM/GMAP.

!     Modifications.
!     --------------
!        Original : FEBRUARY 1992.
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        F. Vana       26-Aug-2008 optimization for NEC
!        F. Courteille 16-Sep-2009 optimization for NEC SX-9
!        H Petithomme (Dec 2020): optimisation, no SIMD anymore
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB     ,JPIA
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMB
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDN
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDX
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDLAT(KPROMB,KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDLO(KPROMB,KFLEV,1:2)
INTEGER(KIND=JPIM),INTENT(IN)    :: KL0(KPROMB,KFLEV,1:2)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDVER(KPROMB,KFLEV)
REAL(KIND=JPRB),TARGET,INTENT(IN):: PXSL(KPROMA*(KFLDX-KFLDN+1))
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PXF(KPROMB,KFLEV)
!     ------------------------------------------------------------------

! optim: prefetch next 6 levels (best), use JPIA (64-bit) and PX1/2 for pxsl base address
#ifdef __INTEL_COMPILER
INTEGER, PARAMETER :: JPREFETCH=6
#else
INTEGER, PARAMETER :: JPREFETCH=0
#endif
INTEGER(KIND=JPIM) :: JLEV,JROF,KLEVP
INTEGER(KIND=JPIA) :: KP

REAL(KIND=JPRB) :: ZINFLO1,ZINFLO2,ZSUPLO1,ZSUPLO2,ZSUP,ZINF
REAL(KIND=JPRB),CONTIGUOUS,POINTER :: PX1(:),PX2(:)
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('LAITLI',0,ZHOOK_HANDLE)

! optim: initial conversion for 64-bit indexing (important)
KP = KPROMA

! optim: use KP+KP, not 2*KP (relatively important)
PX1 => PXSL(KP+1:)
PX2 => PXSL(KP+KP+1:)

KLEVP = MIN(JPREFETCH,KFLEV)

#ifdef __INTEL_COMPILER
! optim: prefetch the first 2 levels
DO JROF=KSTART,KPROF,2
  CALL MM_PREFETCH(PX1(KL0(JROF,KLEVP,1)+1),0)
  CALL MM_PREFETCH(PX1(KL0(JROF,KLEVP,2)+1),0)
  CALL MM_PREFETCH(PX2(KL0(JROF,KLEVP,1)+1),0)
  CALL MM_PREFETCH(PX2(KL0(JROF,KLEVP,2)+1),0)
ENDDO
#endif

#ifdef _CRAYFTN
  !DIR$ PREFERVECTOR
#endif
DO JLEV=1,KFLEV-JPREFETCH
  ! optim: prefetch the next level if needed
    !DIR$ IVDEP
    !CDIR NODEP
    DO JROF=KSTART,KPROF
#ifdef __INTEL_COMPILER
      ! Prefetches for the first and second loop
      CALL MM_PREFETCH(PX1(KL0(JROF,JLEV+JPREFETCH,1)+1),2)
      CALL MM_PREFETCH(PX1(KL0(JROF,JLEV+JPREFETCH,2)+1),2)
      CALL MM_PREFETCH(PX2(KL0(JROF,JLEV+JPREFETCH,1)+1),2)
      CALL MM_PREFETCH(PX2(KL0(JROF,JLEV+JPREFETCH,2)+1),2)
#endif

      ZSUPLO1 = DWEIGHT(PDLO(JROF,JLEV,1),PX1(KL0(JROF,JLEV,1)+1),PX1(KL0(JROF,JLEV,1)+2))
      ZINFLO1 = DWEIGHT(PDLO(JROF,JLEV,1),PX2(KL0(JROF,JLEV,1)+1),PX2(KL0(JROF,JLEV,1)+2))

      ZSUPLO2 = DWEIGHT(PDLO(JROF,JLEV,2),PX1(KL0(JROF,JLEV,2)+1),PX1(KL0(JROF,JLEV,2)+2))
      ZINFLO2 = DWEIGHT(PDLO(JROF,JLEV,2),PX2(KL0(JROF,JLEV,2)+1),PX2(KL0(JROF,JLEV,2)+2))

      ZSUP = DWEIGHT(PDLAT(JROF,JLEV),ZSUPLO1,ZSUPLO2)
      ZINF = DWEIGHT(PDLAT(JROF,JLEV),ZINFLO1,ZINFLO2)

      PXF(JROF,JLEV) = DWEIGHT(PDVER(JROF,JLEV),ZSUP,ZINF)
    ENDDO
ENDDO

DO JLEV=KFLEV-KLEVP+1,KFLEV
#ifdef _CRAYFTN
    !DIR$ NEXTSCALAR
#endif
    !DIR$ IVDEP
    !CDIR NODEP
    DO JROF=KSTART,KPROF
      ZSUPLO1 = DWEIGHT(PDLO(JROF,JLEV,1),PX1(KL0(JROF,JLEV,1)+1),PX1(KL0(JROF,JLEV,1)+2))
      ZINFLO1 = DWEIGHT(PDLO(JROF,JLEV,1),PX2(KL0(JROF,JLEV,1)+1),PX2(KL0(JROF,JLEV,1)+2))

      ZSUPLO2 = DWEIGHT(PDLO(JROF,JLEV,2),PX1(KL0(JROF,JLEV,2)+1),PX1(KL0(JROF,JLEV,2)+2))
      ZINFLO2 = DWEIGHT(PDLO(JROF,JLEV,2),PX2(KL0(JROF,JLEV,2)+1),PX2(KL0(JROF,JLEV,2)+2))

      ZSUP = DWEIGHT(PDLAT(JROF,JLEV),ZSUPLO1,ZSUPLO2)
      ZINF = DWEIGHT(PDLAT(JROF,JLEV),ZINFLO1,ZINFLO2)

      PXF(JROF,JLEV) = DWEIGHT(PDVER(JROF,JLEV),ZSUP,ZINF)
    ENDDO
ENDDO

IF (LHOOK) CALL DR_HOOK('LAITLI',1,ZHOOK_HANDLE)
CONTAINS
  ! basic function for linear interpolation, efficiently inlined
  ! optim: elemental seems important
  ELEMENTAL REAL(KIND=JPRB) FUNCTION DWEIGHT(PD,Z1,Z2)
    REAL(KIND=JPRB),INTENT(IN) :: PD,Z1,Z2

    DWEIGHT = Z1+PD*(Z2-Z1)
  END FUNCTION
END SUBROUTINE
