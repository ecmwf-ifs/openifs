! (C) Copyright 1989- Meteo-France.

SUBROUTINE SUFPCIP(YDCLIMO,YDGFP,KFPXFLD,YDFPGEOMETRY,YDFPWSTD,YDFPSTRUCT,PLSMIN,KFPMASK,PLSMOUT,KWIC)

! Purpose :
! -------
!   *SUFPCIP* Setup localization of created isolated points (lake or island)
!             in the horizontal interpolation step.

! Interface :
! ---------

!      KFPXFLD : maximum number of fields to extract at a time
!      PLSMIN    : source land/sea mask, including halo
!      PLSMOUT   : target land/sea mask.
!      KWIC      : indicator of missing lake/island

! Externals :
! ---------
!   MPL_ALLREDUCE

! Method :
! ------

! Distributed computation.
! Always consider the neighbouring 12 rather than 4 points.
! Rely upon land-sea mask possible values : 1 = land ; 0 = sea
! Action depends on the vegetation index (from climatology)
! Numbering of the points (I is the interpolation point):
!                     13       5       6      14

!                      7       1       2       8
!                                 (I)
!                      9       3       4      10

!                     15      11      12      16

! Reference :
! ---------

! Author :
! ------
!   03-Aug-2005 R. El Khatib  *METEO-FRANCE*

! Modifications :
! -------------
!   K. Yessad: 28-Feb-2007: optimisation of distributed memory in FULL-POS.
!   G.Mozdzynski (Jan 2011): OOPS cleaning, use of derived type SL_STRUCT
!   R. El Khatib  24-Jul-2012 NFPDISTRIB replaced by LFPDISTRIB
!   R. El Khatib 13-Dec-2012 Fullpos buffers reshaping
!   T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!   R. El Khatib 27-Sep-2013 Differentiation between input grid and output grid
! End Modifications
!-----------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMMP0   , ONLY : NPROC
USE YOMLUN   , ONLY : NULOUT
USE YOMFPGEOMETRY, ONLY : TFPGEOMETRY, LFPDISTRIB
USE TYPE_FPOSBUF, ONLY : FPOSBUF
USE YOMAFN, ONLY : ALL_FPDSPHY_TYPES
USE YOMCLI   , ONLY : NTPLAC   ,NTPMER
USE MPL_MODULE, ONLY : MPL_ALLREDUCE
USE EINT_MOD, ONLY : SL_STRUCT
USE YOMWFPB  , ONLY : TFPWSTD
USE YOMFP4L, ONLY : IFPSEARCH

!-----------------------------------------------------------------------------

IMPLICIT NONE

TYPE (FPOSBUF),  INTENT(IN) :: YDCLIMO
TYPE(ALL_FPDSPHY_TYPES), INTENT(IN) :: YDGFP
INTEGER(KIND=JPIM),INTENT(IN)  :: KFPXFLD
TYPE (TFPGEOMETRY),INTENT(IN)  :: YDFPGEOMETRY
TYPE (TFPWSTD),    INTENT(IN)  :: YDFPWSTD
TYPE(SL_STRUCT),   INTENT(IN)  :: YDFPSTRUCT
REAL(KIND=JPRB)   ,INTENT(IN)  :: PLSMIN(YDFPSTRUCT%NASLB1)
INTEGER(KIND=JPIM),INTENT(IN)  :: KFPMASK
REAL(KIND=JPRB)   ,INTENT(IN)  :: PLSMOUT(YDFPGEOMETRY%YFPGEO_DEP%NFPROMA,KFPMASK,YDFPGEOMETRY%YFPGEO_DEP%NFPBLOCS)
INTEGER(KIND=JPIM),INTENT(OUT) :: KWIC(YDFPGEOMETRY%YFPGEO%NFPROMA,YDFPGEOMETRY%YFPGEO%NFPBLOCS)

!-----------------------------------------------------------------------------

INTEGER(KIND=JPIM) :: IT, J, JBLOC, IST, IEND, ILAKES, ISLANDS, IFLDS, IOFF, IPTR
INTEGER(KIND=JPIM) :: IN(12), IMAX(1), IMIN(1)
INTEGER(KIND=JPIM) :: IVEG_DEP(YDFPGEOMETRY%YFPGEO_DEP%NFPROMA,YDFPGEOMETRY%YFPGEO_DEP%NFPBLOCS)
INTEGER(KIND=JPIM) :: IWIC_DEP(YDFPGEOMETRY%YFPGEO_DEP%NFPROMA,YDFPGEOMETRY%YFPGEO_DEP%NFPBLOCS)
INTEGER(KIND=JPIM) :: ISNGLL(4) ! local counter of singular points
INTEGER(KIND=JPIM) :: ISNGLG(4) ! global counter of singular points

REAL(KIND=JPRB) :: ZVEG(YDFPGEOMETRY%YFPGEO%NFPROMA,1,YDFPGEOMETRY%YFPGEO%NFPBLOCS)
REAL(KIND=JPRB) :: ZVEG_DEP(YDFPGEOMETRY%YFPGEO_DEP%NFPROMA,1,YDFPGEOMETRY%YFPGEO_DEP%NFPBLOCS)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------------

#include "fptratod.intfb.h"
#include "fptrdtoa.intfb.h"

#include "abor1.intfb.h"

!-----------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUFPCIP',0,ZHOOK_HANDLE)
ASSOCIATE(YDFPGEO_DEP=>YDFPGEOMETRY%YFPGEO_DEP, YDFPGEO=>YDFPGEOMETRY%YFPGEO, &
 & YDFPGIND=>YDFPGEOMETRY%YFPGIND)
ASSOCIATE(RFPBUF=>YDCLIMO%FPBUF, YDRQCLI=>YDCLIMO%YRQPHY, ML0=>YDFPWSTD%ML0, &
 & NFPRGPL_DEP=>YDFPGEO_DEP%NFPRGPL, NFPBLOCS_DEP=>YDFPGEO_DEP%NFPBLOCS, &
 & NFPROMA_DEP=>YDFPGEO_DEP%NFPROMA, NFPEND_DEP=>YDFPGEO_DEP%NFPEND, &
 & NFPBLOCS=>YDFPGEO%NFPBLOCS, NFPROMA=>YDFPGEO%NFPROMA, NFPEND=>YDFPGEO%NFPEND)

!-----------------------------------------------------------------------------

IFLDS=1
IPTR=IFPSEARCH(YDRQCLI,YDGFP%IVEG%ICOD)
ZVEG(:,1,:)=RFPBUF(:,IPTR,:)

! Transposition between arrival geometry and departure geometry DM-partitions:
IF (LFPDISTRIB(YDFPGEOMETRY)) THEN
! provisional padding with zeros :
  ZVEG_DEP(:,:,:)=0._JPRB
  CALL FPTRATOD(KFPXFLD,YDFPGEO_DEP,YDFPGIND,YDFPGEO,IFLDS,ZVEG,ZVEG_DEP)
  DO JBLOC=1,NFPBLOCS_DEP
    DO J=1,NFPEND_DEP(JBLOC)
      IVEG_DEP(J,JBLOC)=NINT(ZVEG_DEP(J,1,JBLOC))
    ENDDO
  ENDDO
ELSE
  DO JBLOC=1,NFPBLOCS
    DO J=1,NFPEND(JBLOC)
      IVEG_DEP(J,JBLOC)=NINT(ZVEG(J,1,JBLOC))
    ENDDO
  ENDDO
ENDIF

ISNGLL(:)=0

IF (SIZE(ML0,DIM=2) < 4) CALL ABOR1('SUFPCIP : ML0 TOO SMALL')

DO JBLOC=1,NFPBLOCS_DEP
  IST=1
  IEND=NFPEND_DEP(JBLOC)
  IOFF=NFPROMA_DEP*(JBLOC-1)
  DO J=IST,IEND
    IT     = NINT(PLSMOUT(J,1,JBLOC))
    IN( 1) = NINT(PLSMIN(ML0(J,2,JBLOC)+1))
    IN( 2) = NINT(PLSMIN(ML0(J,2,JBLOC)+2))
    IN( 3) = NINT(PLSMIN(ML0(J,3,JBLOC)+1))
    IN( 4) = NINT(PLSMIN(ML0(J,3,JBLOC)+2))
    IN( 5) = NINT(PLSMIN(ML0(J,1,JBLOC)+1))
    IN( 6) = NINT(PLSMIN(ML0(J,1,JBLOC)+2))
    IN( 7) = NINT(PLSMIN(ML0(J,2,JBLOC)))
    IN( 8) = NINT(PLSMIN(ML0(J,2,JBLOC)+3))
    IN( 9) = NINT(PLSMIN(ML0(J,3,JBLOC)))
    IN(10) = NINT(PLSMIN(ML0(J,3,JBLOC)+3))
    IN(11) = NINT(PLSMIN(ML0(J,4,JBLOC)+1))
    IN(12) = NINT(PLSMIN(ML0(J,4,JBLOC)+2))
    IMAX(1)=MAXVAL(IN)
    IMIN(1)=MINVAL(IN)
    IF (IT == 0 .AND. IMIN(1) == 1 .AND. IVEG_DEP(J,JBLOC) == NTPLAC) THEN
      ! An identified lake has been created
      IWIC_DEP(J,JBLOC)=1
      ISNGLL(1)=ISNGLL(1)+1
    ELSEIF (IT == 0 .AND. IMIN(1) == 1 .AND. IVEG_DEP(J,JBLOC) == NTPMER) THEN
      ! An unidentified lake has been created
      IWIC_DEP(J,JBLOC)=1
      ISNGLL(2)=ISNGLL(2)+1
    ELSEIF (IT == 1 .AND. IMAX(1) == 0) THEN
      ! An island has been created
      IWIC_DEP(J,JBLOC)=2
      ISNGLL(3)=ISNGLL(3)+1
    ELSEIF (IVEG_DEP(J,JBLOC) == NTPLAC) THEN
      ! An identified lake has NOT been created => force climatology
      IWIC_DEP(J,JBLOC)=1
      ISNGLL(4)=ISNGLL(4)+1
    ELSE
      ! At least 1 of the neighbouring points is same kind as target point
      IWIC_DEP(J,JBLOC)=0
    ENDIF
  ENDDO
ENDDO

! Transposition between departure geometry and arrival geometry DM-partitions:
IF (LFPDISTRIB(YDFPGEOMETRY)) THEN
  DO JBLOC=1,NFPBLOCS_DEP
    DO J=1,NFPEND_DEP(JBLOC)
      ZVEG_DEP(J,1,JBLOC)=REAL(IWIC_DEP(J,JBLOC))
    ENDDO
  ENDDO
  CALL FPTRDTOA(KFPXFLD,YDFPGEO_DEP,YDFPGIND,YDFPGEO,IFLDS,ZVEG_DEP,ZVEG)
  DO JBLOC=1,NFPBLOCS
    DO J=1,NFPEND(JBLOC)
      KWIC(J,JBLOC)=NINT(ZVEG(J,1,JBLOC))
    ENDDO
  ENDDO
ELSE
  KWIC(:,:)=IWIC_DEP(:,:)
ENDIF

ISNGLG(:) = ISNGLL(:)
IF (NPROC > 1) CALL MPL_ALLREDUCE(ISNGLG,'SUM',CDSTRING='SUFPCIP:')

ILAKES=ISNGLG(1)+ISNGLG(2)+ISNGLG(4)
ISLANDS=ISNGLG(3)

WRITE(UNIT=NULOUT,FMT='('' SINGULAR POINTS FOR HORIZONTAL INTERPOLATIONS : '',&
 & I4, '' ISOLATED LAKE(S) + '',I4,'' ISOLATED ISLAND(S)'')') ILAKES,ISLANDS
WRITE(UNIT=NULOUT,FMT='('' IDENTIFIED LAKES   : '',I4)') ISNGLG(1)
WRITE(UNIT=NULOUT,FMT='('' UNIDENTIFIED LAKES : '',I4)') ISNGLG(2)
WRITE(UNIT=NULOUT,FMT='('' MISSED LAKES       : '',I4)') ISNGLG(4)

!     ---------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUFPCIP',1,ZHOOK_HANDLE)
END SUBROUTINE SUFPCIP
