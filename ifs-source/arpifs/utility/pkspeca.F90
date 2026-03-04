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

SUBROUTINE PKSPECA(YDGEOMETRY,PSPEC,KFIELDG)

!**** *PKSPECA*  - PACK SPECTRAL DATA THROUGH ARPEGE FILE

!     PURPOSE.
!     --------
!        To pack model spectral fields by writing to, then reading from a
!         file ARPEGE
!        In SM version : loop on all fields to write & read back
!        In DM version : Fields are distributed between the processors ;
!        Each processor performs its own loop on fields to write & read back
!        on a local file.

!**   INTERFACE.
!     ----------
!       *CALL* *PKSPECA*

!        EXPLICIT ARGUMENTS
!        --------------------
!        PSPEC   : (Local) spectral fields array
!        KFIELDG : actual number of fields in spectral array

!        IMPLICIT ARGUMENTS
!        --------------------
!        None

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
!      RYAD EL KHATIB *METEO-FRANCE*
!      ORIGINAL : 94-04-08

!     MODIFICATIONS.
!     --------------
!      R. El Khatib : 01-03-16 No I/O in packing ; (LELAM and NGRSP) moved to PRESPFPOS
!      Modified : 01-03-20 S. Ivatek-Sahdan call DIWRSPE -> call DIWRSPE0
!                                         & call DISSPEC -> call DISSPEC0
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      D. Paradis & R. El Khatib : 04-07-22 GRIBEX
!      R. El Khatib & O.Spaniel : B_level + fullpos inline, correction
!      R. El Khatib, thanks to Alex Deckmin 12-Jan-2009 Bugfix for large
!                    number of mpi tasks
!      R. El Khatib 07-Aug-2009 Bugfix for GRIB 1 encoding
!      R. El Khatib 25-Mar-2010 Reduce communications overhead
!      P.Marguinaud 28-May-2010 Change SUMPIOH interface
!      P.Marguinaud 11-Sep-2012 Use DIWRSPEC_MOD instead of DIWRSPE0
!      T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!      R. El khatib 16-May-2014 Optimization of in-line/off-line post-processing reproducibility
!      P.Marguinaud 04-Oct-2016 Port to single precision
!      H Petithomme February 2020: initialize ilonga (bug fix)
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN   , ONLY : NSCRTCH
USE YOMOPH0  , ONLY : CNMCA
USE YOMMP0   , ONLY : NPRTRW, MYPROC, MYSETV
USE YOMCT0   , ONLY : LELAM
USE FA_MOD   , ONLY : JPPRCM

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
INTEGER(KIND=JPIM),INTENT(IN)    :: KFIELDG
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPEC(:,:)

REAL(KIND=JPRB), ALLOCATABLE :: ZSPECBUF(:,:), ZVALCO(:)

INTEGER(KIND=JPIM) :: IFIELDS, JF, INBARI, INOMA, ILONGA
INTEGER(KIND=JPIM) :: INGRIB, INBPDG, INBCSP, ISTRON, IPUILA, IDMOPL
INTEGER(KIND=JPIM) :: IFLDS, ISETW
INTEGER(KIND=JPIM), ALLOCATABLE :: ITO(:), IVSET(:)
INTEGER(KIND=JPIM) :: IREP, IPFAOVSZ


CHARACTER(LEN=16) :: CLSCRTCH, CLNOMA
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "set2pe.intfb.h"

#include "gath_spec.h"
#include "egath_spec.h"
#include "dist_spec.h"
#include "edist_spec.h"

#include "fadoco.h"
#include "facono.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('PKSPECA',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NSEFRE=>YDDIM%NSEFRE, NSPEC2=>YDDIM%NSPEC2, NSPEC2G=>YDDIM%NSPEC2G, NRESOL=>YDDIM%NRESOL, &
 & NPOSSP=>YDMP%NPOSSP)
!     ------------------------------------------------------------------

!*       1. SETUP DISTRIBUTION
!           ------------------

IFLDS=KFIELDG

ALLOCATE(ITO(IFLDS))
ALLOCATE(IVSET(IFLDS))
DO JF=1,IFLDS
  IVSET(JF)=MYSETV
  ISETW=MOD(JF-1,NPRTRW)+1
  CALL SET2PE(ITO(JF),0,0,ISETW,MYSETV)
ENDDO
IFIELDS=COUNT(ITO==MYPROC)


!*       2. GATHER FIELDS
!           -------------

ALLOCATE(ZSPECBUF(NSPEC2G,IFIELDS))
IF (LELAM) THEN
  CALL EGATH_SPEC(PSPECG=ZSPECBUF,KFGATHG=IFLDS,KTO=ITO,KVSET=IVSET,PSPEC=PSPEC,LDIM1_IS_FLD=.FALSE.,KRESOL=NRESOL)
ELSE
  CALL GATH_SPEC(PSPECG=ZSPECBUF,KFGATHG=IFLDS,KTO=ITO,KVSET=IVSET,PSPEC=PSPEC,LDIM1_IS_FLD=.FALSE.,KRESOL=NRESOL)
ENDIF

!*       3.PACK THEN UNPACK DATA
!          ---------------------

IF (IFIELDS > 0) THEN
  CALL FASGRA (IREP, CNMCA, IPFAOVSZ)
  ALLOCATE(ZVALCO(NSPEC2G+JPPRCM*IPFAOVSZ))
  INBARI=0
  CALL FANOUV(IREP,NSCRTCH,.FALSE.,CLSCRTCH,'NEW',.TRUE.,.TRUE.,1,1,INBARI,CNMCA)
  DO JF=1,IFIELDS
    ILONGA = NSPEC2G
    CALL FACONO(IREP,NSCRTCH,'S',0,'FIELD',ZSPECBUF(:,JF),.TRUE.,CLNOMA,INOMA,ZVALCO,ILONGA)
    CALL FADOCO(IREP,NSCRTCH,'S',0,'FIELD',.TRUE.,CLNOMA,INOMA,ZVALCO,ILONGA,ZSPECBUF(:,JF))
  ENDDO

  CALL FAIRNO(IREP,NSCRTCH,'DELETE')
  DEALLOCATE(ZVALCO)
ENDIF

!*       2. SCATTER FIELDS
!           --------------

IF (LELAM) THEN
  CALL EDIST_SPEC(PSPECG=ZSPECBUF,KFDISTG=IFLDS,KFROM=ITO,KVSET=IVSET,PSPEC=PSPEC,LDIM1_IS_FLD=.FALSE.,KRESOL=NRESOL)
ELSE
  CALL DIST_SPEC(PSPECG=ZSPECBUF,KFDISTG=IFLDS,KFROM=ITO,KVSET=IVSET,PSPEC=PSPEC,LDIM1_IS_FLD=.FALSE.,KRESOL=NRESOL)
ENDIF

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('PKSPECA',1,ZHOOK_HANDLE)

END SUBROUTINE PKSPECA
