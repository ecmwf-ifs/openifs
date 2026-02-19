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

SUBROUTINE SUGEM1B(YDGEOMETRY)

!**** *SUGEM1B*  - Initialize geometry parameters - first part - b

!     Purpose.
!     --------
!           Initialize geometry
!           Fills some quantities in YOMGEM

!        Initialize some quantities of YOMGEM:
!         * part 1: R4JP,RC2P1,RC2M1,RCOR0,RCOR1
!         * part 2: NLOEN,NMEN
!         * part 3: NSTAGP,NTSTAGP
!         * Part 4: RDELXN.

!**   Interface.
!     ----------
!        *CALL* *SUGEM1B

!        Explicit arguments :
!        --------------------
!        None

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!        None

!     Reference
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      MPP Group *ECMWF*
!      Original : 95-10-01

!     Modifications.
!     --------------
!      K. Yessad (Aug 2009): prune conf 912; externalise conf 911.
!      G. Mozdzynski (Feb 2011): OOPS cleaning, use of derived type TCSGEOM and TCSGLEG
!      K. Yessad (jul 2011): move YRGSGEOM and YRCSGEOM calculations in SUGEM2.
!      K. Yessad (jul 2011): reorder calculations for better consistency LAM/global model.
!      K. Yessad (dec 2011): ensure consistency with LAM model (comments); rm useless RCOR2.
!      T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!      K. Yessad (July 2014): Move some variables.
!      K. Yessad (Dec 2016): Prune obsolete options.
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST   , ONLY : RPI, RA
USE YOMCT0   , ONLY : LALLOPR
USE YOMLUN   , ONLY : NULOUT
USE YOMMP0   , ONLY : NPRINTLEV, LOUTPUT, MYSETW, MY_REGION_EW

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY), INTENT(INOUT) :: YDGEOMETRY
INTEGER(KIND=JPIM) :: IENDL, IENDLAT, IGLGLO,&
 & IROF, ISTAL, ISTALAT, IU, JGL
LOGICAL :: LLP
REAL(KIND=JPRB) :: ZC, ZX
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
#ifdef __INTEL_COMPILER
INTRINSIC :: SHAPE ! Fails with Intel compiler as it thinks we use ATLAS shape function
#endif


!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUGEM1B',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NDGENL=>YDDIM%NDGENL, NDGLG=>YDDIM%NDGLG, NDGLL=>YDDIM%NDGLL, &
 & NDGSAL=>YDDIM%NDGSAL, NDLON=>YDDIM%NDLON, &
 & NLOEN=>YDGEM%NLOEN, NLOENG=>YDGEM%NLOENG, NMEN=>YDGEM%NMEN, &
 & NMENG=>YDGEM%NMENG, &
 & R4JP=>YDGEM%R4JP, RC2M1=>YDGEM%RC2M1, RC2P1=>YDGEM%RC2P1, RCOR0=>YDGEM%RCOR0, &
 & RCOR1=>YDGEM%RCOR1, RDELXN=>YDGEM%RDELXN, RMUCEN=>YDGEM%RMUCEN, &
 & RSTRET=>YDGEM%RSTRET, &
 & NONL=>YDMP%NONL, NPTRFLOFF=>YDMP%NPTRFLOFF, NPTRLS=>YDMP%NPTRLS)
!     ------------------------------------------------------------------

!*       1.     Initialize constants
!                --------------------

! ky: LAM counterpart of these calculations is in part 1 of SUEGEM1B.

R4JP=(4._JPRB*REAL(NDGLG,JPRB)+2)/(4._JPRB*RPI)
ZC=RSTRET
RC2P1=ZC*ZC+1.0_JPRB
RC2M1=ZC*ZC-1.0_JPRB
! if c > 1. :
! RCOR0/sin(latP)= a/b + (1-(a/b)**2)*ln(sqrt((a+b)/(a-b)))
! RCOR1/sin(latP)= sqrt(3)* (1-(a/b)**2) * (1-(a/b)*ln(sqrt((a+b)/(a-b))) )
! if c = 1. :
! RCOR0= 0
! RCOR1= sin(latP)/sqrt(3)
! a=(c+1/c)/2 ; b=(c-1/c)/2 ; latP = latitude of the pole on the real sphere
IF (ABS(RSTRET-1.0_JPRB) < 1.E-14_JPRB) THEN
  RCOR0=0._JPRB
  RCOR1=RMUCEN/SQRT(3._JPRB)
ELSE
  ZX=RC2P1/RC2M1
  RCOR0=RMUCEN*(ZX+(1.0_JPRB-ZX*ZX)*LOG(ZC))
  RCOR1=RMUCEN*SQRT(3._JPRB)*(1.0_JPRB-ZX*ZX)*(1.0_JPRB-ZX*LOG(ZC))
ENDIF
WRITE(UNIT=NULOUT,FMT='('' CORIOLIS EXPANSION : '',&
 & ''RCOR0 = '',E13.7,''  RCOR1 = '',E13.7)') RCOR0,RCOR1

!     ------------------------------------------------------------------

!*       2.    Initialize NLOEN, NMEN.
!              -----------------------

! ky: LAM counterpart of these calculations is in part 2 of SUEGEM1B.

ISTAL=1
IENDL=NDGLL

! * Compute NMEN and NLOEN:
DO JGL=ISTAL,IENDL
  IGLGLO = NPTRLS(MYSETW)+JGL-1
  NMEN(JGL)=NMENG(IGLGLO)
  NLOEN(JGL)=NLOENG(IGLGLO)
ENDDO

! * Print NLOEN, NMEN.
IF (LOUTPUT.AND.NPRINTLEV >= 1) THEN
  WRITE(UNIT=NULOUT,FMT='('' (JGL,NLOEN,NMEN) '')')
  WRITE(UNIT=NULOUT,FMT='(8(1X,''('',I4,I4,I4,'')''))')&
   & (JGL,NLOEN(JGL),NMEN(JGL),JGL=ISTAL,IENDL)  
ENDIF

!     ------------------------------------------------------------------

!*       3.     Compute NSTAGP
!               --------------

! ky: LAM counterpart of these calculations is in part 3 of SUEGEM1B.

ISTALAT=NDGSAL
IENDLAT=NDGENL

LLP = NPRINTLEV >= 1.OR. LALLOPR
IU = NULOUT
ALLOCATE(YDGEM%NSTAGP(ISTALAT:IENDLAT))
IF(LLP)WRITE(IU,'(1X,"ARRAY ",A10," ALLOCATED ",8I8)') 'NSTAGP   ',&
 & SIZE(YDGEM%NSTAGP ),SHAPE(YDGEM%NSTAGP )

IROF=1
DO JGL=ISTALAT,IENDLAT
  IF (NONL(NPTRFLOFF+JGL,MY_REGION_EW) > 0) THEN
    YDGEM%NSTAGP(JGL)=IROF
    IROF=IROF+NONL(NPTRFLOFF+JGL,MY_REGION_EW)
  ELSE
    YDGEM%NSTAGP(JGL)=0
  ENDIF
ENDDO

IF (LOUTPUT.AND.NPRINTLEV >= 1) THEN
  WRITE(UNIT=NULOUT,FMT='('' NSTAGP'')')
  WRITE(UNIT=NULOUT,FMT='(16(1X,I6))')(YDGEM%NSTAGP(JGL),JGL=ISTALAT,IENDLAT)
ENDIF

!     ------------------------------------------------------------------

!*       4.     Compute RDELXN
!               --------------

RDELXN=(2.0_JPRB*RPI*RA)/REAL(NDLON,JPRB)
WRITE(UNIT=NULOUT,FMT='(A,E16.7)') ' RDELXN= ',RDELXN

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUGEM1B',1,ZHOOK_HANDLE)
END SUBROUTINE SUGEM1B
