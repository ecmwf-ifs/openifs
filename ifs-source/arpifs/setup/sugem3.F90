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

SUBROUTINE SUGEM3(YDGEOMETRY)

!**** *SUGEM3*  Set-up of geometry: part 3

!     Purpose.
!     -------
!        Computes NDGSAH and NDGENH.
!        Reallocates NLOEN, NMEN, RLATI and fills them.

!     Interface.
!     ---------
!        *CALL*  *SUGEM3*

!        Explicit arguments : None
!        ------------------
!        Implicit arguments :
!        ------------------
!        see above "USE MODULE"

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Author.
!     ------
!        K. Yessad (Dec 2013) after SUSC2B.

!     Modifications.
!     -------------
!      K. Yessad (July 2014): Move some variables.
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK

USE PARDIM   , ONLY : JPSLWIDE
USE YOMLUN   , ONLY : NULOUT
USE YOMMP0   , ONLY : NPRINTLEV, MYSETW
USE YOMCT0   , ONLY : LALLOPR

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY), INTENT(INOUT) :: YDGEOMETRY
INTEGER(KIND=JPIM) :: IU, JGLLOC, IGLGLO
LOGICAL :: LLP

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
#ifdef __INTEL_COMPILER
INTRINSIC :: SHAPE ! Fails with Intel compiler as it thinks we use ATLAS shape function
#endif


!    -------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUGEM3',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM, &
  & YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, YDCSGLEG=>YDGEOMETRY%YRCSGLEG)
ASSOCIATE(NDGENG=>YDDIM%NDGENG, NDGENH=>YDDIM%NDGENH, NDGENL=>YDDIM%NDGENL, &
 & NDGLL=>YDDIM%NDGLL, NDGSAG=>YDDIM%NDGSAG, NDGSAH=>YDDIM%NDGSAH, &
 & NDGSAL=>YDDIM%NDGSAL, &
 & NGPTOT=>YDGEM%NGPTOT, NLOENG=>YDGEM%NLOENG, NMENG=>YDGEM%NMENG, &
 & NFRSTLOFF=>YDMP%NFRSTLOFF, NPTRLS=>YDMP%NPTRLS)
!    -------------------------------------------------------------------

!*        1. SETUP.
!            ------

LLP = NPRINTLEV >= 1.OR. LALLOPR
IU = NULOUT

NDGSAH=MAX(NDGSAG,NDGSAL+NFRSTLOFF-JPSLWIDE)-NFRSTLOFF
NDGENH=MIN(NDGENG,NDGENL+NFRSTLOFF+JPSLWIDE)-NFRSTLOFF

IF (ALLOCATED(YDGEM%NLOEN)) DEALLOCATE(YDGEM%NLOEN)
IF (ALLOCATED(YDGEM%NMEN)) DEALLOCATE(YDGEM%NMEN)
IF (ALLOCATED(YDCSGLEG%RLATI)) DEALLOCATE(YDCSGLEG%RLATI)

WRITE(NULOUT,'('' Allocations in SUGEM3:'')')
ALLOCATE(YDGEM%NLOEN  (NDGSAL-JPSLWIDE:NDGLL+JPSLWIDE))
IF(LLP)WRITE(IU,9) 'NLOEN    ',SIZE(YDGEM%NLOEN    ),SHAPE(YDGEM%NLOEN    )
ALLOCATE(YDGEM%NMEN   (NDGSAL-JPSLWIDE:NDGLL+JPSLWIDE))
IF(LLP)WRITE(IU,9) 'NMEN     ',SIZE(YDGEM%NMEN     ),SHAPE(YDGEM%NMEN     )
ALLOCATE(YDCSGLEG%RLATI(NDGSAH:NDGENH))
IF(LLP)WRITE(IU,9) 'RLATI    ',SIZE(YDCSGLEG%RLATI),SHAPE(YDCSGLEG%RLATI)

DO JGLLOC=NDGSAL-JPSLWIDE,NDGLL+JPSLWIDE
  IGLGLO=NPTRLS(MYSETW)+JGLLOC-1
  IF( IGLGLO >= NDGSAG.AND.IGLGLO <= NDGENG )THEN
    YDGEM%NLOEN(JGLLOC)=NLOENG(IGLGLO)
    YDGEM%NMEN (JGLLOC)=NMENG(IGLGLO)
  ENDIF
ENDDO
DO JGLLOC=NDGSAH,NDGENH
  IGLGLO=JGLLOC+NFRSTLOFF
  YDCSGLEG%RLATI(JGLLOC)=YDCSGLEG%RLATIG(IGLGLO)
ENDDO

!    -------------------------------------------------------------------

!*        2. PRINTINGS.
!            ----------

WRITE(NULOUT,'('' Printings in SUGEM3:'')')
WRITE(NULOUT,'(2X,A,I5,A,I5)') ' NDGSAH = ',NDGSAH,' NDGENH = ',NDGENH
WRITE(NULOUT,'('' NGPTOT ='',I9)') NGPTOT

9 FORMAT(1X,'ARRAY ',A10,' ALLOCATED ',8I8)

!    -------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUGEM3',1,ZHOOK_HANDLE)
END SUBROUTINE SUGEM3
