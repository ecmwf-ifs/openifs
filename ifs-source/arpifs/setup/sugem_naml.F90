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

SUBROUTINE SUGEM_NAML(YDGEOMETRY,KSUPERSEDE)

!**** *SUGEM_NAML*  - Initialize geometry parameters: namelist variables

!     Purpose.
!     --------
!           Initialize geometry: namelist variables

!**   Interface.
!     ----------
!        *CALL* *SUGEM_NAML

!        Explicit arguments :
!        --------------------

!        Implicit arguments :
!        --------------------
!        See #include below

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!                   SURGRI -  read NLOENG on namelist
!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      Mats Hamrud and Philippe Courtier  *ECMWF*
!      Original : 87-10-15

!     Modifications.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      R. El Khatib 05-02-25 Bugfix after the dimension of NUMEN
!      F. Bouyssel 05-03-01 REFLRHC
!      K. Yessad (Sept 2008): Gaussian lat and weights => dummy arg
!      K. Yessad (Sept 2008): Prune conf 951.
!      K. Yessad (Nov 2008): rename arp/SUGAW into arp/SUGAWA.
!      K. Yessad (Aug 2009): prune conf 912, externalise conf 911.
!      K. Yessad (May 2012): split; first part in SUGEM_NAML.
!      T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!      K. Yessad (July 2014): Move some variables.
!      R. El Khatib 04-Aug-2014 Pruning of the conf. 927/928
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMARG   , ONLY : NULOEN, NSUPERSEDE, NUHTYP, NUSTTYP, USTRET, UMUCEN, ULOCEN  
USE YOMCT0   , ONLY : LECMWF, LALLOPR
USE YOMMP0   , ONLY : NPRINTLEV, LOUTPUT
USE YOMLUN   , ONLY : NULOUT, NULNAM
USE YOMMP0   , ONLY : MYPROC

!     ------------------------------------------------------------------

IMPLICIT NONE

!     ------------------------------------------------------------------

TYPE(GEOMETRY), INTENT(INOUT),TARGET :: YDGEOMETRY
INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: KSUPERSEDE
REAL(KIND=JPRB), POINTER :: RMUCEN, RLOCEN, RSTRET, RNLGINC
INTEGER(KIND=JPIM), POINTER :: NSTTYP, NHTYP
LOGICAL, POINTER :: LNONHYD_GEOM

#include "namgem.nam.h"

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZNLGINC
INTEGER(KIND=JPIM) ::  JGL, IDNLGINC
LOGICAL :: LLGRID, LLP, LLSUPERSEDE

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "posnam.intfb.h"
#include "surgri.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUGEM_NAML',0,ZHOOK_HANDLE)

!     ------------------------------------------------------------------
!*       1.    Allocations (moved from SUALLO)
!              -------------------------------

LLP = NPRINTLEV >= 1.OR. LALLOPR

ALLOCATE(YDGEOMETRY%YRGEM%NLOEN  (YDGEOMETRY%YRDIM%NDGSAG:YDGEOMETRY%YRDIM%NDGENG))
IF(LLP)WRITE(NULOUT,9) 'NLOEN    ',SIZE(YDGEOMETRY%YRGEM%NLOEN),SHAPE(YDGEOMETRY%YRGEM%NLOEN)
ALLOCATE(YDGEOMETRY%YRGEM%NLOENG (YDGEOMETRY%YRDIM%NDGSAG:YDGEOMETRY%YRDIM%NDGENG))
IF(LLP)WRITE(NULOUT,9) 'NLOENG   ',SIZE(YDGEOMETRY%YRGEM%NLOENG),SHAPE(YDGEOMETRY%YRGEM%NLOENG)

ASSOCIATE(NCMAX=>YDGEOMETRY%YRDIM%NCMAX, NDGENG=>YDGEOMETRY%YRDIM%NDGENG, NDGLG=>YDGEOMETRY%YRDIM%NDGLG, &
 & NDGSAG=>YDGEOMETRY%YRDIM%NDGSAG, NDGSUR=>YDGEOMETRY%YRDIM%NDGSUR, NDLON=>YDGEOMETRY%YRDIM%NDLON, &
 & NSMAX=>YDGEOMETRY%YRDIM%NSMAX, &
 & NLOEN=>YDGEOMETRY%YRGEM%NLOEN, NLOENG=>YDGEOMETRY%YRGEM%NLOENG)

! Associate pointers for variables in namelist
LNONHYD_GEOM => YDGEOMETRY%LNONHYD_GEOM
RMUCEN  => YDGEOMETRY%YRGEM%RMUCEN
RLOCEN  => YDGEOMETRY%YRGEM%RLOCEN
RSTRET  => YDGEOMETRY%YRGEM%RSTRET
RNLGINC => YDGEOMETRY%YRGEM%RNLGINC
NSTTYP  => YDGEOMETRY%YRGEM%NSTTYP
NHTYP   => YDGEOMETRY%YRGEM%NHTYP 

!     ------------------------------------------------------------------
!*       2.    Initialize YOMGEM.
!              ------------------

!*    2.1   SET DEFAULT VALUES.

LNONHYD_GEOM = .FALSE.
IF (PRESENT(KSUPERSEDE)) THEN
  LLSUPERSEDE=(KSUPERSEDE == 1)
ELSE
  LLSUPERSEDE=(NSUPERSEDE == 1)
ENDIF
!        2.1.1 Set implicit default values

IF (LLSUPERSEDE) THEN
  RMUCEN =UMUCEN
  RLOCEN =ULOCEN
  RSTRET =USTRET
  NSTTYP =NUSTTYP
  NHTYP  =NUHTYP
  ! We try to guess the value of RNLGINC : 0., 0.1, 0.2, ...,0.5, ... or 2. : 
  IF (ABS(RSTRET-1.0_JPRB)>100.0_JPRB*TINY(1.0_JPRB)) THEN
    ZNLGINC=REAL(NDLON-3-2*NSMAX,JPRB)/REAL(NSMAX,JPRB)
  ELSE
    ZNLGINC=REAL(NDLON-1-2*NSMAX,JPRB)/REAL(NSMAX,JPRB)
  ENDIF
  IDNLGINC=NINT(ZNLGINC*10._JPRB)
  RNLGINC=REAL(IDNLGINC,JPRB)/10._JPRB
ELSE
  RMUCEN =1._JPRB
  RLOCEN =0._JPRB
  RSTRET =1._JPRB
  NSTTYP =1
  NHTYP  =0
  RNLGINC=1._JPRB
ENDIF

!        2.1.2 Modify default values according to LECMWF

IF (LECMWF) THEN
  LLGRID=NSMAX > (NDLON+3)/3
  IF (LLGRID) THEN
    RNLGINC=0._JPRB
  ELSE
    RNLGINC=1._JPRB
  ENDIF
ENDIF

!*    2.2 MODIFY VALUES

!        2.2.1  READ NAMELIST.

CALL POSNAM(NULNAM,'NAMGEM')
READ(NULNAM,NAMGEM)

!        2.2.2  OVERWRITE NAMELIST WITH FILE FRAME ARGUMENTS
IF (LLSUPERSEDE) THEN
  IF (RSTRET /= USTRET) THEN
    RSTRET=USTRET
    WRITE(NULOUT,*) 'RSTRET  OVERWRITTEN BY FILE FRAME'
  ENDIF
  IF (RMUCEN /= UMUCEN) THEN
    RMUCEN=UMUCEN
    WRITE(NULOUT,*) 'RMUCEN  OVERWRITTEN BY FILE FRAME'
  ENDIF
  IF (RLOCEN /= ULOCEN) THEN
    RLOCEN=ULOCEN
    WRITE(NULOUT,*) 'RLOCEN  OVERWRITTEN BY FILE FRAME'
  ENDIF
  IF (NHTYP /= NUHTYP) THEN
    NHTYP=NUHTYP
    WRITE(NULOUT,*) 'NHTYP  OVERWRITTEN BY FILE FRAME'
  ENDIF
  IF (NSTTYP /= NUSTTYP) THEN
    NSTTYP=NUSTTYP
    WRITE(NULOUT,*) 'NSTTYP  OVERWRITTEN BY FILE FRAME'
  ENDIF
ENDIF

!        2.2.3  RESET VARIABLES AND TEST

IF (NSTTYP == 1) THEN
  RMUCEN=1._JPRB
  RLOCEN=0._JPRB
ENDIF

RNLGINC=MIN(2.0_JPRB,MAX(0.0_JPRB,RNLGINC))

! RSTRET < 1. not allowed
IF ((RSTRET-1.0_JPRB) < -1.E-14_JPRB) THEN
  CALL ABOR1( ' SUGEM_NAML: RSTRET < 1. NOT ALLOWED, CHANGE THE POLE !')  
ENDIF

! Checkings on NCMAX.
IF (ABS(RSTRET-1.0_JPRB) < 1.E-14_JPRB) THEN
  ! No stretching
  IF (NCMAX /= NSMAX) CALL ABOR1(' SUGEM_NAML: IF NO STRETCHING NCMAX=NSMAX !')
ELSE
  ! Stretching
  IF (NCMAX<NSMAX) CALL ABOR1(' SUGEM_NAML: CHECK NCMAX !')
ENDIF

!        2.2.4  PRINTINGS

WRITE(NULOUT,*) ' '
WRITE(NULOUT,*) ' --- Printings in SUGEM_NAML:'
WRITE(NULOUT,FMT='(''LNONHYD_GEOM = '',L2)') LNONHYD_GEOM
WRITE(NULOUT,FMT='('' NSTTYP = '',I6,&
 & '' RMUCEN = '',E13.7,'' RLOCEN = '',E13.7,&
 & '' RSTRET = '',E13.7)')NSTTYP,RMUCEN,RLOCEN,RSTRET  
WRITE(NULOUT,FMT='('' NHTYP  = '',I6)') NHTYP
WRITE(NULOUT,FMT='('' RNLGINC  = '',E13.7)') RNLGINC

!     ------------------------------------------------------------------
!*       3.    Initialize NLOENG.
!              ------------------

IF (NHTYP == 0) THEN
  DO JGL = NDGSAG, NDGENG
    NLOENG(JGL)=NDLON
  ENDDO
ELSEIF (NHTYP == 2) THEN
  IF (LLSUPERSEDE) THEN
    ! read from file frame.
    DO JGL = 1, (NDGLG+1)/2
      NLOENG(JGL)=NULOEN(JGL)
    ENDDO
    DO JGL = (NDGLG+1)/2+1, NDGLG
      NLOENG(JGL)=NULOEN(NDGLG-JGL+1)
    ENDDO
    IF (NDGSAG <= 0) THEN
      NLOENG(0)=NLOENG(1)
      NLOENG(NDGLG+1)=NLOENG(NDGLG)
    ENDIF
  ELSE
    ! read from namelist.
    CALL SURGRI(YDGEOMETRY,NULNAM,NULOUT,NDGLG,NDGSAG,NDGENG,NLOENG)
  ENDIF

  ! extra-polar latitudes extension.
  IF (NDGSUR > 0) THEN
!DIR$ IVDEP
!OCL NOVREC
    DO JGL = 1, NDGSUR
      NLOENG(1-JGL)=NLOENG(JGL)
      NLOENG(NDGLG+JGL)=NLOENG(NDGLG+1-JGL)
    ENDDO
  ENDIF
ELSE
  WRITE(NULOUT,FMT='('' ERROR NHTYP '')')
  CALL ABOR1(' SUGEM_NAML ')
ENDIF

!*    Prints NLOENG

IF (LOUTPUT.AND.NPRINTLEV >= 1) THEN
  WRITE(NULOUT,FMT='('' (JGL,NLOENG) '')')
  WRITE(NULOUT,FMT='(8(1X,''('',I4,I5,'')''))')&
   & (JGL,NLOENG(JGL),JGL=NDGSAG,NDGENG)  
ENDIF

!     ------------------------------------------------------------------
9 FORMAT(1X,'ARRAY ',A10,' ALLOCATED ',8I8)

!     ------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUGEM_NAML',1,ZHOOK_HANDLE)
END SUBROUTINE SUGEM_NAML
