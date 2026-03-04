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

SUBROUTINE PKSURFA (YDGEOMETRY,YDSURF,CDNAME,PGP,LD2)

!     *PKSURFA*  - PACK SURFACE DATA THROUGH ARPEGE FILE

!     PURPOSE.
!     --------
!        To pack model surface fields by writing to, then reading from a file
!        ARPEGE.
!        Fields are distributed between the processors ; each processor
!        performs its own loop on fields to write & read back on a local file.

!     INTERFACE.
!     ----------
!       *CALL* *PKSURFA*

!        EXPLICIT ARGUMENTS
!        --------------------
!        PGP    - data array

!        IMPLICIT ARGUMENTS
!        ------------------

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     AUTHOR.
!     -------
!        RYAD EL KHATIB *METEO-FRANCE*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 16-Oct-2006
!     Apr 2008  K. Yessad: use DISGRID instead of DISGRID_C + cleanings
!        R. El Khatib 30-Jul-2009 Support for GRIB 1 needs names
!    R. El Khatib : 01-Apr-2010 Overhead reduction
!    P.Marguinaud : 28-05-2010 Change SUMPIOH interface
!    R. El Khatib : 14-Mar-2012 Allow blank-named fields
!    P.Marguinaud : Initialize INBARI (avoid valgrind warning)
!    T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!    R. El khatib 16-May-2014 Optimization of in-line/off-line post-processing reproducibility
!    P.Marguinaud 04-Oct-2016 Port to single precision
!     ------------------------------------------------------------------

USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN   , ONLY : NSCRTCH, NULOUT
USE YOMOPH0  , ONLY : CNMCA
USE YOMCT0   , ONLY : LELAM
USE YOMMP0   , ONLY : NPROC, MYPROC
USE FA_MOD   , ONLY : JPPRCM

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TSURF), INTENT(INOUT) :: YDSURF
CHARACTER(LEN=*)  , INTENT(IN) :: CDNAME(:)
REAL(KIND=JPRB)   , INTENT(INOUT) :: PGP(:,:,:)
LOGICAL,INTENT(IN), OPTIONAL :: LD2(:)
!     IASPIO(j): Processor number which packs fields #j
!     IFIELDS  : Number of fields packed by MYPROC
!     ZREALG   : global array of fields

REAL(KIND=JPRB),ALLOCATABLE :: ZREALG(:,:), ZVALCO(:)
INTEGER(KIND=JPIM) :: INFD(NPROC), IFLDOFF(NPROC)
INTEGER(KIND=JPIM) :: IEND, IFIELDS, IST, IFLDS,JROC
INTEGER(KIND=JPIM) :: JF, INBARI, INOMA, ILONGA, IG, ILEV
INTEGER(KIND=JPIM) :: INGRIB,INBPDG,INBCSP,ISTRON,IPUILA,IDMOPL
INTEGER(KIND=JPIM), ALLOCATABLE :: ITO(:)
INTEGER(KIND=JPIM) :: IREP, IPFAOVSZ


CHARACTER(LEN=16) :: CLSCRTCH, CLNOMA
CHARACTER(LEN=12) :: CLSUFFIX
CHARACTER(LEN=4)  :: CLPREFIX

LOGICAL :: LLPRINT, LLPACKING

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "sumpioh.intfb.h"
#include "abor1.intfb.h"

#include "gath_grid.h"
#include "egath_grid.h"
#include "dist_grid.h"
#include "edist_grid.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('PKSURFA',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDGEM=>YDGEOMETRY%YRGEM)
ASSOCIATE(NPROMA=>YDDIM%NPROMA,NRESOL=>YDDIM%NRESOL,NGPTOTG=>YDGEM%NGPTOTG)

!        1.  CONTROL ARGUMENTS
!            -----------------

IFLDS=SIZE(PGP,DIM=2)
IF (SIZE(CDNAME) /= IFLDS) CALL ABOR1('PKSURFA: SIZE(CDNAME) /= SIZE(PGP,DIM=2)')
IF (SIZE(LD2) /= IFLDS) CALL ABOR1('PKSURFA: SIZE(LD2) /= SIZE(PGP,DIM=2)')

!        2.  SETUP DISTIBUTION OF FIELDS
!            ---------------------------

CALL SUMPIOH(NPROC,NPROC,IFLDS,INFD,IFLDOFF)
IFIELDS=INFD(MYPROC)
ALLOCATE(ITO(IFLDS))
DO JROC=1,NPROC
  IST=IFLDOFF(JROC)+1
  IEND=IFLDOFF(JROC)+INFD(JROC)
  ITO(IST:IEND)=JROC
ENDDO

!        3.  GATHER FIELDS
!            -------------

IF (LELAM) THEN
  ALLOCATE(ZREALG(NGPTOTG,IFIELDS))
  CALL EGATH_GRID(PGPG=ZREALG,KPROMA=NPROMA,KRESOL=NRESOL,KFGATHG=IFLDS,KTO=ITO,PGP=PGP)
ELSE
  ALLOCATE(ZREALG(NGPTOTG,IFIELDS))
  CALL GATH_GRID(PGPG=ZREALG,KPROMA=NPROMA,KRESOL=NRESOL,KFGATHG=IFLDS,KTO=ITO,PGP=PGP)
ENDIF

!        4.  PACK THEN UNPACK DATA
!            ---------------------

IF (IFIELDS > 0) THEN
  INBARI=0
  CALL FANOUV(IREP,NSCRTCH,.FALSE.,CLSCRTCH,'NEW',.TRUE.,.TRUE.,1,1,INBARI,&
   & CNMCA)  
  CALL FAVEUR(IREP,NSCRTCH,INGRIB,INBPDG,INBCSP,ISTRON,IPUILA,IDMOPL)
! This prit was for debugging  LLPRINT=(INGRIB==-1 .OR. INGRIB==3)
  LLPRINT=.FALSE.
  IF (LLPRINT) WRITE(NULOUT,*) 'BEFORE/AFTER PACKING : MIN AVE MAX'
  ILEV=0
  CALL FASGRA (IREP, CNMCA, IPFAOVSZ)
  ALLOCATE(ZVALCO(NGPTOTG+JPPRCM*IPFAOVSZ))
ENDIF
DO JF=1,IFIELDS
  IG=IFLDOFF(MYPROC)+JF
  LLPACKING=.TRUE.
  IF (LD2(IG)) THEN
    LLPACKING=(CDNAME(IG)(1:16)/=YDSURF%YSD_VF%YZ0F%CNAME(1:16) .AND. &
             & CDNAME(IG)(1:16)/=YDSURF%YSD_VF%YZ0RLF%CNAME(1:16) .AND. &
             & CDNAME(IG)(1:16)/=YDSURF%YSD_VV%YZ0H%CNAME(1:16))
    IF (CDNAME(IG)(1:3)=='CLS' .OR. CDNAME(IG)(1:3)=='CLP') THEN
      CLPREFIX=CDNAME(IG)(1:3)
      CLSUFFIX=CDNAME(IG)(4:MIN(16,LEN(TRIM(CDNAME(IG)))))
    ELSE
      CLPREFIX=CDNAME(IG)(1:4)
      CLSUFFIX=CDNAME(IG)(5:MIN(16,LEN(TRIM(CDNAME(IG)))))
    ENDIF
    IF (CLPREFIX == ' ') CLPREFIX='SURF'
    IF (CLSUFFIX == ' ') CLSUFFIX='FIELD       '
  ELSE
    IF (CDNAME(IG)(1:4)=='SURF') THEN
      CLPREFIX=CDNAME(IG)(1:4)
      CLSUFFIX=CDNAME(IG)(5:MIN(16,LEN(TRIM(CDNAME(IG)))))
    ELSEIF (CDNAME(IG)(1:4)=='PROF') THEN
      CLPREFIX=CDNAME(IG)(1:4)
      CLSUFFIX=CDNAME(IG)(5:MIN(16,LEN(TRIM(CDNAME(IG)))))
    ELSEIF (CDNAME(IG)(1:4)=='SOMM') THEN
      CLPREFIX=CDNAME(IG)(1:4)
      CLSUFFIX=CDNAME(IG)(5:MIN(16,LEN(TRIM(CDNAME(IG)))))
    ELSEIF (CDNAME(IG)(1:3)=='CLS') THEN
      CLPREFIX=CDNAME(IG)(1:3)
      CLSUFFIX=CDNAME(IG)(4:MIN(16,LEN(TRIM(CDNAME(IG)))))
    ELSE
      CLPREFIX='S   '
      CLSUFFIX=CDNAME(IG)(1:MIN(12,LEN(TRIM(CDNAME(IG)))))
    ENDIF
    IF (CLPREFIX == ' ') CLPREFIX='SURF'
    IF (CLSUFFIX == ' ') CLSUFFIX='FIELD       '
  ENDIF
  IF (LLPACKING) THEN
    IF (LLPRINT) THEN
      IF (LD2(IG)) THEN
        WRITE(NULOUT,*) TRIM(CLPREFIX),TRIM(CLSUFFIX), &
         & MINVAL(ZREALG(1:NGPTOTG,JF)), &
         & SUM(ZREALG(1:NGPTOTG,JF))/NGPTOTG,MAXVAL(ZREALG(1:NGPTOTG,JF))
      ELSE
        WRITE(NULOUT,*) TRIM(CLPREFIX),TRIM(CLSUFFIX), &
         & MINVAL(ZREALG(1:NGPTOTG,JF)), &
         & SUM(ZREALG(1:NGPTOTG,JF))/NGPTOTG,MAXVAL(ZREALG(1:NGPTOTG,JF))
      ENDIF
    ENDIF
    CALL FACOND(IREP,NSCRTCH,CLPREFIX,ILEV,CLSUFFIX,ZREALG(1,JF),.FALSE., &
     & CLNOMA,INOMA,ZVALCO,ILONGA)
    CALL FADECO(IREP,NSCRTCH,CLPREFIX,ILEV,CLSUFFIX,.FALSE.,CLNOMA,INOMA,&
     & ZVALCO,ILONGA,ZREALG(1,JF))  
    IF (LLPRINT) THEN
      IF (LD2(IG)) THEN
        WRITE(NULOUT,*) TRIM(CLPREFIX),TRIM(CLSUFFIX), &
         & MINVAL(ZREALG(1:NGPTOTG,JF)) &
         & ,SUM(ZREALG(1:NGPTOTG,JF))/NGPTOTG,MAXVAL(ZREALG(1:NGPTOTG,JF))
      ELSE
        WRITE(NULOUT,*) TRIM(CLPREFIX),TRIM(CLSUFFIX), &
         & MINVAL(ZREALG(1:NGPTOTG,JF)), &
         & SUM(ZREALG(1:NGPTOTG,JF))/NGPTOTG,MAXVAL(ZREALG(1:NGPTOTG,JF))
      ENDIF
    ENDIF
  ENDIF
ENDDO
IF (IFIELDS > 0) THEN
  DEALLOCATE(ZVALCO)
  CALL FAIRNO(IREP,NSCRTCH,'DELETE')
ENDIF

!        5.  SCATTER FIELDS
!            --------------

IF (LELAM) THEN
  CALL EDIST_GRID(PGPG=ZREALG,KPROMA=NPROMA,KRESOL=NRESOL,KFDISTG=IFLDS,KFROM=ITO,PGP=PGP)
ELSE
  CALL DIST_GRID(PGPG=ZREALG,KPROMA=NPROMA,KRESOL=NRESOL,KFDISTG=IFLDS,KFROM=ITO,PGP=PGP)
ENDIF
DEALLOCATE(ZREALG)

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('PKSURFA',1,ZHOOK_HANDLE)
END SUBROUTINE PKSURFA

