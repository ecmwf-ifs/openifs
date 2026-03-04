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

SUBROUTINE SUSLAD1(YDGEOMETRY,YDSLREP,YDSL,YDAD)

!**** *SUSLAD1*  Initialise data structures for SL adjoint

!     Purpose.
!     -------
!           Initialise NGPTOTAD and NADCORE required for SL adjoint.
!           NGPTOTAD represents the number of grid points spanning 
!           a standard width halo around this processors core points.
!           NADCORE are the offsets of the above points in a double
!           standard width semi-langrangian buffer.
!           Note the great care is taken that these data structures 
!           contain no duplicate points (i.e. mirror latitudes or 
!           latitude extension points).
!           It is also important for bit reproducibility that the order 
!           of points in NADCORE is identical to that when executed
!           on a single processor.

!     Method.
!     -------
!        See documentation

!     Author.
!     ------
!        George Mozdzynski  * ECMWF *
!        Original  : 00-08-01

!     Modifications.
!     -------------
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        G.Mozdzynski (Jan 2011): OOPS cleaning, use of derived type SL_STRUCT
!        T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!        K. Yessad (July 2014): Move some variables.
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMLUN   , ONLY : NULOUT
USE YOMCT0   , ONLY : LALLOPR
USE YOMMP0   , ONLY : NPROC, NPRINTLEV
USE YOMSLREP , ONLY : TSLREP
USE EINT_MOD , ONLY : SL_STRUCT

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY) ,INTENT(INOUT) :: YDGEOMETRY
TYPE(TSLREP)   ,INTENT(INOUT) :: YDSLREP
TYPE(SL_STRUCT),INTENT(INOUT) :: YDSL
TYPE(SL_STRUCT),INTENT(INOUT) :: YDAD
INTEGER(KIND=JPIM) :: JSL, JAD, JOFF, IOFF, ILAT, IU

LOGICAL :: LLERROR, LLFOUND, LLP, LLDEBUG

INTEGER(KIND=JPIM), ALLOCATABLE :: IADCORE(:)
INTEGER(KIND=JPIM), ALLOCATABLE :: IADLAT(:)
INTEGER(KIND=JPIM), ALLOCATABLE :: ISLBUF(:,:)

INTEGER(KIND=JPIM) :: IMAP(YDGEOMETRY%YRDIM%NDLON)

LOGICAL :: LLUSED(YDGEOMETRY%YRDIM%NDLON)
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
#ifdef __INTEL_COMPILER
INTRINSIC :: SHAPE ! Fails with Intel compiler as it thinks we use ATLAS shape function
#endif


!     ------------------------------------------------------------------

#include "abor1.intfb.h"

!    -------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUSLAD1',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NDLON=>YDDIM%NDLON, &
 & NGPTOT=>YDGEM%NGPTOT, NLOENG=>YDGEM%NLOENG, &
 & NGPTOTAD=>YDSLREP%NGPTOTAD)
!     ------------------------------------------------------------------

LLDEBUG=.FALSE.

LLP = NPRINTLEV >= 1.OR. LALLOPR
IU = NULOUT

ALLOCATE(IADCORE(YDAD%NASLB1))
ALLOCATE(IADLAT(YDAD%NASLB1))
IADLAT(:)=99999999
IF( NPROC>1 )THEN
  NGPTOTAD=0
  DO JSL=1,SIZE(YDSL%NSLMAP,DIM=2)
    LLFOUND=.FALSE.
    DO JAD=1,SIZE(YDAD%NSLMAP,DIM=2)
      IF( YDSL%NSLMAP(1,JSL) == YDAD%NSLMAP(1,JAD) )THEN
        LLFOUND=.TRUE.
        ILAT=YDSL%NSLMAP(1,JSL)
        IMAP(:)=99999999
        DO JOFF=YDAD%NSLMAP(2,JAD),YDAD%NSLMAP(2,JAD)+YDAD%NSLMAP(3,JAD)-1
          IOFF=MOD(JOFF+NLOENG(ILAT)-1,NLOENG(ILAT))+1
          IMAP(IOFF)=JOFF-YDAD%NSLMAP(2,JAD)+1
        ENDDO
        LLUSED(:)=.FALSE.
        DO JOFF=YDSL%NSLMAP(2,JSL),YDSL%NSLMAP(2,JSL)+YDSL%NSLMAP(3,JSL)-1
          IOFF=MOD(JOFF+NLOENG(ILAT)-1,NLOENG(ILAT))+1
          LLUSED(IOFF)=.TRUE.
        ENDDO
        DO JOFF=1,NLOENG(ILAT)
          IF( LLUSED(JOFF) )THEN
            NGPTOTAD=NGPTOTAD+1
            IF( NGPTOTAD > YDAD%NASLB1 )THEN
              CALL ABOR1('SUSLAD1: INTERNAL ERROR NGPTOTAD > YDAD%NASLB1 ')
            ENDIF
            IF( IMAP(JOFF) > NLOENG(ILAT) )THEN
              CALL ABOR1('SUSLAD1: INTERNAL ERRROR IN BUILDING IADCORE)')
            ENDIF
            IADCORE(NGPTOTAD)=YDAD%NSLMAP(4,JAD)+IMAP(JOFF)
            IADLAT(NGPTOTAD)=ILAT
          ENDIF
        ENDDO
      ENDIF
    ENDDO
    IF( .NOT.LLFOUND )THEN
      CALL ABOR1('SUSLAD1: INTERNAL ERROR DURING SEARCH FOR NADCORE POINTS')
    ENDIF
  ENDDO
ELSE
  NGPTOTAD=NGPTOT
ENDIF

WRITE(NULOUT,'("SUSLAD1: YDAD%NASLB1=",I10," NGPTOTAD=",I10," NGPTOT=",I10)')&
 & YDAD%NASLB1,NGPTOTAD,NGPTOT  

ALLOCATE(YDSLREP%NADCORE(NGPTOTAD))
IF(LLP)WRITE(IU,9) 'NADCORE ',SIZE(YDSLREP%NADCORE),SHAPE(YDSLREP%NADCORE)
IF( NPROC>1 )THEN
  YDSLREP%NADCORE(1:NGPTOTAD)=IADCORE(1:NGPTOTAD)
ELSE
  YDSLREP%NADCORE(1:NGPTOTAD)=YDAD%NSLCORE(1:NGPTOTAD)
ENDIF
IF(ALLOCATED(IADCORE)) DEALLOCATE(IADCORE)

IF( LLDEBUG )THEN
  ! confirm that all points in yrad%nslcore are in nadcore
  DO JAD=1,NGPTOTAD
    WRITE(NULOUT,'("SUSLAD1: JAD=",I10," NADCORE(JAD)=",I10)')&
     & JAD,YDSLREP%NADCORE(JAD)  
  ENDDO
  DO JSL=1,NGPTOT
    WRITE(NULOUT,'("SUSLAD1: JSL=",I10," YDAD%NSLCORE(JSL)=",I10)')&
     & JSL,YDAD%NSLCORE(JSL)  
  ENDDO
  LLERROR=.FALSE.
  DO JSL=1,NGPTOT
    LLFOUND=.FALSE.
    DO JAD=1,NGPTOTAD
      IF( YDAD%NSLCORE(JSL) == YDSLREP%NADCORE(JAD) )THEN
        IF( .NOT.LLFOUND )THEN
          LLFOUND=.TRUE.
        ELSE
          WRITE(NULOUT,'("SUSLAD1: JSL=",I10," YDAD%NSLCORE(JSL)=",I10,&
           & " OCCURS MORE THAN ONCE IN NADCORE")')JSL,YDAD%NSLCORE(JSL)  
          LLERROR=.TRUE.
        ENDIF
      ENDIF
    ENDDO
    IF( .NOT.LLFOUND )THEN
      WRITE(NULOUT,'("SUSLAD1: JSL=",I10," YDAD%NSLCORE(JSL)=",I10,&
       & " DOES NOT OCCUR IN NADCORE")')JSL,YDAD%NSLCORE(JSL)  
      LLERROR=.TRUE.
    ENDIF
  ENDDO
  IF( LLERROR )THEN
    ALLOCATE (ISLBUF(YDAD%NASLB1,3))
    ISLBUF(:,:)=99999999
    DO JAD=1,NGPTOTAD
      ISLBUF(YDSLREP%NADCORE(JAD),1)=JAD
      ISLBUF(YDSLREP%NADCORE(JAD),3)=IADLAT(JAD)
    ENDDO
    DO JSL=1,NGPTOT
      ISLBUF(YDAD%NSLCORE(JSL),2)=JSL
    ENDDO
    DO JOFF=1,YDAD%NASLB1
      WRITE(NULOUT,'("SUSLAD1: LAT=",I4," JOFF=",I10," NADCORE POINT ",I10,&
       & "   YDAD%NSLCORE POINT ",I10)')&
       & ISLBUF(JOFF,3),JOFF,ISLBUF(JOFF,1),ISLBUF(JOFF,2)  
    ENDDO
    IF (ALLOCATED(ISLBUF)) DEALLOCATE(ISLBUF)
    CALL FLUSH(NULOUT)
    CALL ABOR1('SUSLAD1: YDAD%NSLCORE NOT A SUBSET OF NADCORE')
  ENDIF

ENDIF

IF(ALLOCATED(IADLAT)) DEALLOCATE(IADLAT)

IF (NPROC > 1) THEN
  ALLOCATE(YDAD%MASK_SLTOT(YDAD%NASLB1))
  IF(LLP)WRITE(IU,9) 'YDAD%MASK_SLTOT',SIZE(YDAD%MASK_SLTOT),SHAPE(YDAD%MASK_SLTOT)
  YDAD%MASK_SLTOT(:)=0
  ALLOCATE(YDSLREP%LADCORE(YDAD%NASLB1))
  IF(LLP)WRITE(IU,9) 'LADCORE',SIZE(YDSLREP%LADCORE),SHAPE(YDSLREP%LADCORE)
  YDSLREP%LADCORE(:)=.FALSE.
  DO JSL=1,NGPTOTAD
    YDSLREP%LADCORE(YDSLREP%NADCORE(JSL))=.TRUE.
  ENDDO
ENDIF

9 FORMAT(1X,'ARRAY ',A10,' ALLOCATED ',8I8)
!    -------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUSLAD1',1,ZHOOK_HANDLE)
END SUBROUTINE SUSLAD1
