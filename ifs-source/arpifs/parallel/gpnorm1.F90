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

SUBROUTINE GPNORM1(YDGEOMETRY,PGP,KFIELDS,LDLEVELS,CDLABEL,PNORMS)

!**** *gpnorm1 * - Routine to calculate reproducible grid point norm

!     Purpose.
!     --------
! Routine to calculate reproducible grid point norm.

!**   Interface.
!     ----------
!        *CALL* *GPNORM1(...) *

!        Explicit arguments :
!        --------------------
!        PGP     : if nproc is > 1 then PGP is the distributed field
!        KFIELDS  : number of levels
!        LDLEVELS : =F Only print average for all levels. Min and max are global

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      L. Isaksen *ECMWF*
!      Original : 96-07-07

!     Modifications.
!     --------------
!      Modified : 02-05-10 by S.Checker : changed MIMIMUM to MINIMUM
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      L.Isaksen : Support printing of average for all levels
!      G.Mozdzynski : 25-Oct-2005 use all available tasks and threads for norms
!      Y.Tremolet    03-Aug-2005 Optional output argument
!      N.Wedi    : Proper area-weighted mean
!      G.Mozdzynski  19-Sept-2008 gpnorm optimisations
!      A.Bogatchev: 12-Jun-2009 call to egpnorm_trans
!      T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!      K. Yessad (July 2014): Move some variables.
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB, JPRD
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMMP0       , ONLY : NPRINTLEV, MYPROC, NPROC
USE YOMCT0       , ONLY : LALLOPR, LELAM
USE YOMLUN       , ONLY : NULOUT
USE MPL_MODULE   , ONLY : MPL_BROADCAST

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
INTEGER(KIND=JPIM),INTENT(IN)    :: KFIELDS 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGP(YDGEOMETRY%YRDIM%NPROMA,KFIELDS,YDGEOMETRY%YRDIM%NGPBLKS) 
LOGICAL           ,INTENT(IN)    :: LDLEVELS
CHARACTER(LEN=*)  ,INTENT(IN),OPTIONAL :: CDLABEL
REAL(KIND=JPRB), INTENT(OUT), OPTIONAL :: PNORMS(3,KFIELDS) 

!     ------------------------------------------------------------------

LOGICAL :: LLP
REAL(KIND=JPRB),ALLOCATABLE :: ZGP(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE :: ZAVE(:), ZMAX(:), ZMIN(:)
REAL(KIND=JPRB) :: ZMAXBLK(YDGEOMETRY%YRDIM%NGPBLKS), ZMINBLK(YDGEOMETRY%YRDIM%NGPBLKS)
REAL(KIND=JPRD) :: ZAVEALL, ZMAXALL, ZMINALL
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

INTEGER(KIND=JPIM) :: IGPNMASTER,IFIELDS
INTEGER(KIND=JPIM) :: JF,JL,ITAG,JKGLO,IST,IEND,IBL

CHARACTER(LEN=16) :: CLABEL
LOGICAL :: LLAVE_ONLY

!     ------------------------------------------------------------------

#include "egpnorm_trans.h"
#include "gpnorm_trans.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('GPNORM1',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NGPBLKS=>YDDIM%NGPBLKS, NPROMA=>YDDIM%NPROMA, NRESOL=>YDDIM%NRESOL, &
 & NGPTOT=>YDGEM%NGPTOT)
!     ------------------------------------------------------------------

LLP = NPRINTLEV >= 1.OR. LALLOPR

CLABEL='OUTPUT'
IF(PRESENT(CDLABEL)) CLABEL = CDLABEL
LLAVE_ONLY = .NOT. (PRESENT(PNORMS) .OR. LDLEVELS)
IGPNMASTER=1

IF(LLAVE_ONLY) THEN
  ALLOCATE(ZGP(NPROMA,1,NGPBLKS))
  ALLOCATE(ZAVE(1))
  ALLOCATE(ZMIN(1))
  ALLOCATE(ZMAX(1))
  ZMAXBLK(:)=PGP(1,1,1)
  ZMINBLK(:)=PGP(1,1,1)
  CALL GSTATS(1427,0)
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JKGLO,IST,IEND,IBL,JF,JL)
  DO JKGLO=1,NGPTOT,NPROMA
    IST=1
    IEND=MIN(NPROMA,NGPTOT-JKGLO+1)
    IBL=(JKGLO-1)/NPROMA+1
    DO JL=IST,IEND
      ZGP(JL,1,IBL) = 0.0_JPRB
    ENDDO
    DO JF=1,KFIELDS
      DO JL=IST,IEND
        ZGP(JL,1,IBL) = ZGP(JL,1,IBL)+PGP(JL,JF,IBL)
        ZMINBLK(IBL) = MIN(ZMINBLK(IBL),PGP(JL,JF,IBL))
        ZMAXBLK(IBL) = MAX(ZMAXBLK(IBL),PGP(JL,JF,IBL))
      ENDDO
    ENDDO
  ENDDO
!$OMP END PARALLEL DO
  CALL GSTATS(1427,1)
  ZMIN(1)=MINVAL(ZMINBLK(1:NGPBLKS))
  ZMAX(1)=MAXVAL(ZMAXBLK(1:NGPBLKS))
  IFIELDS=1
  IF( .NOT. LELAM ) THEN
   CALL GPNORM_TRANS(ZGP,IFIELDS,NPROMA,ZAVE,ZMIN,ZMAX,LLAVE_ONLY,NRESOL)
  ELSE
   CALL EGPNORM_TRANS(ZGP,IFIELDS,NPROMA,ZAVE,ZMIN,ZMAX,LLAVE_ONLY,NRESOL)
  ENDIF
ELSE
  ALLOCATE(ZAVE(KFIELDS))
  ALLOCATE(ZMIN(KFIELDS))
  ALLOCATE(ZMAX(KFIELDS))
  IFIELDS=KFIELDS
  IF( .NOT. LELAM ) THEN
   CALL GPNORM_TRANS(PGP,IFIELDS,NPROMA,ZAVE,ZMIN,ZMAX,LLAVE_ONLY,NRESOL)
  ELSE
   CALL EGPNORM_TRANS(PGP,IFIELDS,NPROMA,ZAVE,ZMIN,ZMAX,LLAVE_ONLY,NRESOL)
  ENDIF
ENDIF

IF(ALLOCATED(ZGP)) DEALLOCATE(ZGP)

IF (PRESENT(PNORMS)) THEN
  IF (NPROC > 1) THEN
    ITAG=12345
    CALL MPL_BROADCAST (ZAVE,ITAG,IGPNMASTER,CDSTRING='GPNORM1-AVE:')
    ITAG=ITAG+100
    CALL MPL_BROADCAST (ZMIN,ITAG,IGPNMASTER,CDSTRING='GPNORM1-MIN:')
    ITAG=ITAG+100
    CALL MPL_BROADCAST (ZMAX,ITAG,IGPNMASTER,CDSTRING='GPNORM1-MAX:')
  ENDIF
  PNORMS(1,:)=ZAVE(:)
  PNORMS(2,:)=ZMIN(:)
  PNORMS(3,:)=ZMAX(:)
ELSE
  IF (MYPROC == IGPNMASTER.OR.NPROC == 1) THEN
    WRITE(NULOUT,'(5(1X,A))') 'GPNORM',CLABEL,'    AVERAGE         ',&
     & '     MINIMUM        ','      MAXIMUM        '
    IF (KFIELDS>0) THEN
      ZAVEALL=ZAVE(1)
      ZMINALL=ZMIN(1)
      ZMAXALL=ZMAX(1)
      DO JF=2,IFIELDS
        ZAVEALL=ZAVE(JF)+ZAVEALL
        ZMINALL=MIN(ZMINALL,ZMIN(JF))
        ZMAXALL=MAX(ZMAXALL,ZMAX(JF))
      ENDDO
      ZAVEALL=ZAVEALL/KFIELDS
      WRITE(NULOUT,'(9X,A5,3(1X,E21.15))') 'AVE  ',ZAVEALL,ZMINALL,ZMAXALL
      IF (LDLEVELS) THEN
        DO JF=1,KFIELDS
          WRITE(NULOUT,'(9X,I5,3(1X,E21.15))') JF,ZAVE(JF),ZMIN(JF),ZMAX(JF)
        ENDDO
      ENDIF
    ENDIF
  ENDIF
ENDIF
IF(ALLOCATED(ZAVE)) DEALLOCATE(ZAVE)
IF(ALLOCATED(ZMIN)) DEALLOCATE(ZMIN)
IF(ALLOCATED(ZMAX)) DEALLOCATE(ZMAX)

!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GPNORM1',1,ZHOOK_HANDLE)
END SUBROUTINE GPNORM1
