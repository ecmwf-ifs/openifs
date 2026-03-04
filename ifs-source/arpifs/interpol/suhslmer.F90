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

SUBROUTINE SUHSLMER(YDGEOMETRY,KGLN,KGLX,LDEQMER,PSLHDEPSH,YDHSLMER)

!**** *SUHSLMER*  Sets-up quantities involved in interpolator:
!                 intermediate quantities to compute horizontal meridian weights.
!                 (currently RIPI,RSLD,RSLDW,R3DTW)

!     Purpose.
!     -------

!     Interface.
!     ---------
!        *CALL*  *SUHSLMER*

!        Explicit arguments :
!        ------------------
!         KGLN    : lower latitude   (in)
!         KGLX    : upper latitude   (in)
!         LDEQMER : T if regular meridian spacing (like for plane geometry)  (in)

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
!        Mats Hamrud  * ECMWF * (91-01-15) in SUHSLMER.
!        K. Yessad (Jun 2012) from SUHSLMER code.

!     Modifications.
!     -------------
!      K. Yessad (July 2014): Move some variables.
!      F. Vana  30-Oct-2019: More precision to LSLHDQUAD
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMHSLMER,     ONLY : THSLMER
USE PARKIND1 , ONLY : JPIM, JPRB, JPRD
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK

USE YOMLUN_IFSAUX,ONLY: NULOUT

! arp/ifs dependencies to be solved later
USE YOMMP0   , ONLY : NPRINTLEV
USE YOMCT0   , ONLY : LALLOPR

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
INTEGER(KIND=JPIM),INTENT(IN)    :: KGLN
INTEGER(KIND=JPIM),INTENT(IN)    :: KGLX
LOGICAL           ,INTENT(IN)    :: LDEQMER
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSLHDEPSH
TYPE(THSLMER)     ,INTENT(INOUT) :: YDHSLMER
!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: IGLGLO, IU, JGL
LOGICAL :: LLP
REAL(KIND=JPRD) :: ZD1, ZD2, ZD3, ZD4, ZD5, ZD6

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
#ifdef __INTEL_COMPILER
INTRINSIC :: SHAPE ! Fails with Intel compiler as it thinks we use ATLAS shape function
#endif


!    -------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUHSLMER',0,ZHOOK_HANDLE)
!    -------------------------------------------------------------------

! Allocations:
ALLOCATE(YDHSLMER%RIPI(KGLN:KGLX,3))
ALLOCATE(YDHSLMER%RSLD(KGLN:KGLX,3))
ALLOCATE(YDHSLMER%RSLDW(3,3,KGLN:KGLX))
ALLOCATE(YDHSLMER%R3DTW(3,3,KGLN:KGLX))

ASSOCIATE(YDVAB=>YDGEOMETRY%YRVAB,YDVETA=>YDGEOMETRY%YRVETA,YDVFE=>YDGEOMETRY%YRVFE,YDSTA=>YDGEOMETRY%YRSTA, &
 & YDLAP=>YDGEOMETRY%YRLAP,YDCSGLEG=>YDGEOMETRY%YRCSGLEG,&
 & YDCSGEOM=>YDGEOMETRY%YRCSGEOM,YDCSGEOM_NB=>YDGEOMETRY%YRCSGEOM_NB,YDGSGEOM=>YDGEOMETRY%YRGSGEOM, &
 & YDGSGEOM_NB=>YDGEOMETRY%YRGSGEOM_NB, YDSPGEOM=>YDGEOMETRY%YSPGEOM)
ASSOCIATE(NDGLG=>YDGEOMETRY%YRDIM%NDGLG, &
 & R3DTW=>YDHSLMER%R3DTW, RIPI=>YDHSLMER%RIPI, RSLD=>YDHSLMER%RSLD, &
 & RSLDW=>YDHSLMER%RSLDW, &
 & NFRSTLOFF=>YDGEOMETRY%YRMP%NFRSTLOFF)

LLP = NPRINTLEV >= 1.OR. LALLOPR
IU = NULOUT

IF (LLP) THEN
  WRITE(IU,9) 'RIPI     ',SIZE(RIPI),SHAPE(RIPI)
  WRITE(IU,9) 'RSLD     ',SIZE(RSLD),SHAPE(RSLD)
  WRITE(IU,9) 'RSLDW    ',SIZE(RSLDW),SHAPE(RSLDW)
  WRITE(IU,9) 'R3DTW    ',SIZE(R3DTW),SHAPE(R3DTW)
ENDIF

! set-up:
IF(LDEQMER) THEN

  ! Meridian spacing is regular.

  DO JGL= KGLN,KGLX
    IGLGLO=JGL+NFRSTLOFF

    ! * RIPI arrays
    RIPI(JGL,1)=-1.0_JPRD
    RIPI(JGL,2)= 1.0_JPRD
    RIPI(JGL,3)=-1.0_JPRD

    ! * RSLD arrays
    RSLD(JGL,1)= 0.25_JPRB
    RSLD(JGL,2)= 0.25_JPRB
    RSLD(JGL,3)=-0.25_JPRB

    RSLDW(1,1,JGL)=1.0_JPRB-2.0_JPRB*PSLHDEPSH
    RSLDW(1,2,JGL)=PSLHDEPSH
    RSLDW(1,3,JGL)=0.0_JPRB
    RSLDW(2,1,JGL)=PSLHDEPSH
    RSLDW(2,2,JGL)=1.0_JPRB-2.0_JPRB*PSLHDEPSH
    RSLDW(2,3,JGL)=0.0_JPRB
    RSLDW(3,1,JGL)=0.0_JPRB
    RSLDW(3,2,JGL)=PSLHDEPSH
    RSLDW(3,3,JGL)=1.0_JPRB

    R3DTW(1,1,JGL)=0.0_JPRB
    R3DTW(1,2,JGL)=0.5_JPRB
    R3DTW(1,3,JGL)=0.0_JPRB
    R3DTW(2,1,JGL)=0.5_JPRB
    R3DTW(2,2,JGL)=0.0_JPRB
    R3DTW(2,3,JGL)=0.0_JPRB
    R3DTW(3,1,JGL)=0.0_JPRB
    R3DTW(3,2,JGL)=0.5_JPRB
    R3DTW(3,3,JGL)=1.0_JPRB
  ENDDO ! JGL

ELSE

  ! Meridian spacing is irregular and meridian coordinate is computational space latitude.
  DO JGL= KGLN,KGLX
    IGLGLO=JGL+NFRSTLOFF

    IF(IGLGLO >= 0.AND.IGLGLO <= NDGLG) THEN
      ! assumes that ndgsag <= -1 and ndgeng >= ndglg+2

      ! * RIPI arrays
      ZD1=REAL(YDCSGLEG%RLATIG(IGLGLO-1),JPRD)-REAL(YDCSGLEG%RLATIG(IGLGLO),JPRD)
      ZD2=REAL(YDCSGLEG%RLATIG(IGLGLO-1),JPRD)-REAL(YDCSGLEG%RLATIG(IGLGLO+1),JPRD)
      ZD3=REAL(YDCSGLEG%RLATIG(IGLGLO-1),JPRD)-REAL(YDCSGLEG%RLATIG(IGLGLO+2),JPRD)
      ZD4=REAL(YDCSGLEG%RLATIG(IGLGLO  ),JPRD)-REAL(YDCSGLEG%RLATIG(IGLGLO+1),JPRD)
      ZD5=REAL(YDCSGLEG%RLATIG(IGLGLO  ),JPRD)-REAL(YDCSGLEG%RLATIG(IGLGLO+2),JPRD)
      ZD6=REAL(YDCSGLEG%RLATIG(IGLGLO+1),JPRD)-REAL(YDCSGLEG%RLATIG(IGLGLO+2),JPRD)

      RIPI(JGL,1)=-1.0_JPRD/(ZD1*ZD4*ZD5)
      RIPI(JGL,2)= 1.0_JPRD/(ZD2*ZD4*ZD6)
      RIPI(JGL,3)=-1.0_JPRD/(ZD3*ZD5*ZD6)

      ! * RSLD arrays
      RSLD(JGL,1)=(ZD2-ZD6)*ZD4/(ZD1*ZD2+ZD5*ZD6)
      RSLD(JGL,2)=(ZD5-ZD1)*ZD4/(ZD1*ZD2+ZD5*ZD6)
      RSLD(JGL,3)=    -ZD4 *ZD4/(ZD1*ZD2+ZD5*ZD6)

      RSLDW(1,1,JGL)=1.0_JPRD-2.0_JPRD*PSLHDEPSH*ZD4/ZD1
      RSLDW(1,2,JGL)=2.0_JPRD*PSLHDEPSH*ZD4/ZD5
      RSLDW(1,3,JGL)=0.0_JPRD
      RSLDW(2,1,JGL)=2.0_JPRD*PSLHDEPSH*ZD4/ZD2
      RSLDW(2,2,JGL)=1.0_JPRD-2.0_JPRD*PSLHDEPSH*ZD4/ZD6
      RSLDW(2,3,JGL)=0.0_JPRD
      RSLDW(3,1,JGL)=0.0_JPRD
      RSLDW(3,2,JGL)=2.0_JPRD*PSLHDEPSH*(ZD4*ZD4)/(ZD5*ZD6)
      RSLDW(3,3,JGL)=1.0_JPRD

      R3DTW(1,1,JGL)=1.0_JPRD-1.0_JPRD*ZD4/ZD1
      R3DTW(1,2,JGL)=1.0_JPRD*ZD4/ZD5
      R3DTW(1,3,JGL)=0.0_JPRD
      R3DTW(2,1,JGL)=1.0_JPRD*ZD4/ZD2
      R3DTW(2,2,JGL)=1.0_JPRD-1.0_JPRD*ZD4/ZD6
      R3DTW(2,3,JGL)=0.0_JPRD
      R3DTW(3,1,JGL)=0.0_JPRD
      R3DTW(3,2,JGL)=1.0_JPRD*(ZD4*ZD4)/(ZD5*ZD6)
      R3DTW(3,3,JGL)=1.0_JPRD
    ENDIF ! IGLGLO
  ENDDO ! JGL

ENDIF ! LDEQMER

9 FORMAT(1X,'ARRAY ',A10,' ALLOCATED ',8I8)

!    -------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUHSLMER',1,ZHOOK_HANDLE)
END SUBROUTINE SUHSLMER
