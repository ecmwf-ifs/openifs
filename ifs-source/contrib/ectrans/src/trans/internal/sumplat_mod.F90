! (C) Copyright 1995- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE SUMPLAT_MOD
CONTAINS
SUBROUTINE SUMPLAT(KDGL,KPROC,KPROCA,KMYSETA,LDSPLIT,LDEQ_REGIONS,&
                   &KFRSTLAT,KLSTLAT,KFRSTLOFF,KPTRLAT,&
                   &KPTRFRSTLAT,KPTRLSTLAT,KPTRFLOFF,&
                   &PWEIGHT,LDWEIGHTED_DISTR,PMEDIAP,KPROCAGP,&
                   &KMEDIAP,KRESTM,LDSPLITLAT,KMYPROC,KLOEN)

!**** *SUMPLAT * - Initialize gridpoint distrbution in N-S direction

!     Purpose.
!     --------


!**   Interface.
!     ----------
!        *CALL* *SUMPLAT *

!     Explicit arguments - input :
!     --------------------
!                          KDGL       -last  latitude
!                          KPROC      -total number of processors
!                          KPROCA     -number of processors in A direction
!                          KMYSETA    -process number in A direction
!                          LDSPLIT    -true for latitudes shared between sets
!                          LDEQ_REGIONS -true if eq_regions partitioning
!                          PWEIGHT    -weight per grid-point if weighted distribution
!                          LDWEIGHTED_DISTR -true if weighted distribution

!     Explicit arguments - output:
!     --------------------
!                          PMEDIAP    -mean weight per PE if weighted distribution
!                          KMEDIAP    -mean number of grid points per PE
!                          KPROCAGP   -number of grid points per A set
!                          KRESTM     -number of PEs with one extra point
!                          KFRSTLAT   -first latitude row on processor
!                          KLSTLAT    -last  latitude row on processor
!                          KFRSTLOFF  -offset for first latitude in set
!                          KPTRLAT    -pointer to start of latitude
!                          KPTRFRSTLAT-pointer to first latitude
!                          KPTRLSTLAT -pointer to last  latitude
!                          KPTRFLOFF  -offset for pointer to first latitude
!                          LDSPLITLAT -true for latitudes which are split

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
USE TPM_DISTR       ,ONLY : MYPROC

USE SUMPLATB_MOD    ,ONLY : SUMPLATB
USE SUMPLATBEQ_MOD  ,ONLY : SUMPLATBEQ
USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS
!

IMPLICIT NONE


!     * DUMMY:
REAL(KIND=JPRB),INTENT(OUT)    :: PMEDIAP
INTEGER(KIND=JPIM),INTENT(OUT) :: KMEDIAP
INTEGER(KIND=JPIM),INTENT(OUT) :: KRESTM
INTEGER(KIND=JPIM),INTENT(IN)  :: KDGL
INTEGER(KIND=JPIM),INTENT(IN)  :: KPROC
INTEGER(KIND=JPIM),INTENT(IN)  :: KPROCA
INTEGER(KIND=JPIM),INTENT(IN)  :: KMYSETA
REAL(KIND=JPRB),INTENT(IN)     :: PWEIGHT(:)
LOGICAL,INTENT(INOUT)          :: LDWEIGHTED_DISTR
INTEGER(KIND=JPIM),INTENT(OUT) :: KFRSTLAT(:)
INTEGER(KIND=JPIM),INTENT(OUT) :: KLSTLAT(:)
INTEGER(KIND=JPIM),INTENT(OUT) :: KFRSTLOFF
INTEGER(KIND=JPIM),INTENT(OUT) :: KPTRLAT(:)
INTEGER(KIND=JPIM),INTENT(OUT) :: KPTRFRSTLAT(:)
INTEGER(KIND=JPIM),INTENT(OUT) :: KPTRLSTLAT(:)
INTEGER(KIND=JPIM),INTENT(OUT) :: KPTRFLOFF
INTEGER(KIND=JPIM),INTENT(OUT) :: KPROCAGP(KPROCA)
LOGICAL,INTENT(IN)  :: LDSPLIT
LOGICAL,INTENT(IN)  :: LDEQ_REGIONS
LOGICAL,INTENT(OUT) :: LDSPLITLAT(:)
INTEGER(KIND=JPIM),INTENT(IN) :: KMYPROC
INTEGER(KIND=JPIM),INTENT(IN) :: KLOEN(KDGL)

!     * LOCAL:
! === END OF INTERFACE BLOCK ===
INTEGER(KIND=JPIM) :: INDIC(KPROCA),ILAST(KPROCA)

!     LOCAL INTEGER SCALARS
INTEGER(KIND=JPIM) :: IPTRLATITUDE,  JA, JGL

LOGICAL :: LLFOURIER
LOGICAL :: LLDEBUG=.FALSE.

!      -----------------------------------------------------------------

!*       1.    CODE DEPENDING ON 'LELAM': COMPUTATION OF
!              KMEDIAP, KRESTM, INDIC, ILAST.
!              -----------------------------------------
INDIC(:)=0
ILAST(:)=0

IF(LDWEIGHTED_DISTR.AND..NOT.LDEQ_REGIONS)THEN
  CALL ABORT_TRANS ('SUMPLAT: LDWEIGHTED_DISTR=T AND  LDEQ_REGIONS=F NOT SUPPORTED')
ENDIF

IF( LDEQ_REGIONS )THEN
  CALL SUMPLATBEQ(1,KDGL,KPROC,KPROCA,KLOEN,LDSPLIT,LDEQ_REGIONS,&
   &PWEIGHT,LDWEIGHTED_DISTR,PMEDIAP,KPROCAGP,&
   &KMEDIAP,KRESTM,INDIC,ILAST)
ELSE
  LLFOURIER=.FALSE.
  CALL SUMPLATB(1,KDGL,KPROCA,KLOEN,LDSPLIT,LLFOURIER,&
   &KMEDIAP,KRESTM,INDIC,ILAST)
ENDIF

!      -----------------------------------------------------------------

!*       2.    CODE NOT DEPENDING ON 'LELAM': COMPUTATION OF
!              KFRSTLAT TO LDSPLITLAT.
!              ---------------------------------------------


!     * Computation of first and last latitude of processor sets
!       -----------  in grid-point-space -----------------------

IF(KMYPROC==1.AND.LLDEBUG)THEN
  WRITE(0,'("")')
  WRITE(0,'("SUMPLAT_MOD:LDWEIGHTED_DISTR=",L1)')LDWEIGHTED_DISTR
  WRITE(0,'("")')
  DO JA=1,KPROCA
    WRITE(0,'("SUMPLAT_MOD: JA=",I5," ILAST=",I5," INDIC=",I5)')&
    &JA,ILAST(JA),INDIC(JA)
  ENDDO
  WRITE(0,'("")')
  IF( LDEQ_REGIONS .AND. LDSPLIT )THEN
    DO JA=1,KPROCA
      WRITE(0,'("SUMPLAT_MOD: JA=",I5," KPROCAGP=",I12)')&
      &JA,KPROCAGP(JA)
    ENDDO
    WRITE(0,'("")')
  ENDIF
ENDIF

KFRSTLAT(1) = 1
KLSTLAT(KPROCA) = KDGL
DO JA=1,KPROCA-1
  IF ((.NOT. LDSPLIT) .OR. INDIC(JA) == 0) THEN
    KFRSTLAT(JA+1) = ILAST(JA) + 1
    KLSTLAT(JA) = ILAST(JA)
  ELSE
    KFRSTLAT(JA+1) = INDIC(JA)
    KLSTLAT(JA) = INDIC(JA)
  ENDIF
ENDDO
KFRSTLOFF=KFRSTLAT(KMYSETA)-1

!     * Initialise following data structures:-
!       NPTRLAT     (pointer to the start of each latitude)
!       LSPLITLAT   (TRUE if latitude is split over two A sets)
!       NPTRFRSTLAT (pointer to the first latitude of each A set)
!       NPTRLSTLAT  (pointer to the last  latitude of each A set)

DO JGL=1,KDGL
  KPTRLAT  (JGL)=-999
  LDSPLITLAT(JGL)=.FALSE.
ENDDO
IPTRLATITUDE=0
DO JA=1,KPROCA
  DO JGL=KFRSTLAT(JA),KLSTLAT(JA)
    IPTRLATITUDE=IPTRLATITUDE+1
    LDSPLITLAT(JGL)=.TRUE.
    IF( KPTRLAT(JGL) == -999 )THEN
      KPTRLAT(JGL)=IPTRLATITUDE
      LDSPLITLAT(JGL)=.FALSE.
    ENDIF
  ENDDO
ENDDO
DO JA=1,KPROCA
  IF( LDSPLITLAT(KFRSTLAT(JA)) .AND. JA /= 1)THEN
    KPTRFRSTLAT(JA)=KPTRLAT(KFRSTLAT(JA))+1
  ELSE
    KPTRFRSTLAT(JA)=KPTRLAT(KFRSTLAT(JA))
  ENDIF
  IF( LDSPLITLAT(KLSTLAT(JA)) .AND. JA == KPROCA)THEN
    KPTRLSTLAT(JA)=KPTRLAT(KLSTLAT(JA))+1
  ELSE
    KPTRLSTLAT(JA)=KPTRLAT(KLSTLAT(JA))
  ENDIF
ENDDO
KPTRFLOFF=KPTRFRSTLAT(KMYSETA)-1

IF(KMYPROC==1.AND.LLDEBUG)THEN
  DO JGL=1,KDGL
    WRITE(0,'("SUMPLAT_MOD: JGL=",I5," KPTRLAT=",I5," LDSPLITLAT=",L4)')&
    & JGL,KPTRLAT(JGL),LDSPLITLAT(JGL)
  ENDDO
  DO JA=1,KPROCA
    WRITE(0,'("SUMPLAT_MOD: JA=",I5," KFRSTLAT=",I5," KLSTLAT=",I5,&
    & " KPTRFRSTLAT=",I5," KPTRLSTLAT=",I5," KLSTLAT-KFRSTLAT=",I5,&
    & " SUM(G%NLOEN(KFRSTLAT:KLSTLAT))=",I10)')&
    & JA,KFRSTLAT(JA),KLSTLAT(JA),KPTRFRSTLAT(JA),KPTRLSTLAT(JA),&
    & KLSTLAT(JA)-KFRSTLAT(JA),SUM(G%NLOEN(KFRSTLAT(JA):KLSTLAT(JA)))
  ENDDO
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE SUMPLAT
END MODULE SUMPLAT_MOD



