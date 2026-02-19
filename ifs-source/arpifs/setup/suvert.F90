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

#ifdef RS6K
! compiler bug on RS6000
@PROCESS NOOPTIMIZE
#endif
SUBROUTINE SUVERT(YDGEOMETRY)

!**** *SUVERT*  - Routine to initialize vertical coordinate

!     Purpose.
!     --------
!           Initialize the hybrid-cordinate system of the model.

!**   Interface.
!     ----------

!     *CALL* SUVERT
!        Explicit arguments :
!        --------------------

!        Implicit arguments :
!        --------------------
!        see the modules used above.

!     Method.
!     -------

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      Mats Hamrud and Philippe Courtier  *ECMWF*
!      Original : 87-10-15

!     Modifications.
!     --------------
!      30-Jun-2008 J. Masek   Auxiliary quantities for SLHD interpolators.
!      15-Sep-2008 K. Yessad  LREGETA -> LREGETA+LVFE_REGETA.
!      31-Mar-2011 M.Hamrud Intruduce NAMVV0 to have the possibility to
!       force reading NAMVV1. This is because of issues with representing
!       the A and B in GRIB1 and GRIB2 (not bit-identicle)
!      G.Mozdzynski (Feb 2011): OOPS cleaning, use of derived type TVAB, TVETA, TVSLETA
!      R. El Khatib 10-Aug-2011 NIOLEVG management
!      K. Yessad (Mar 2012): code reorganisation
!      R. El Khatib 26-Jul-2012 part 1 moved in SUVV1
!      P. Smolikova and J. Vivoda (Oct 2013): new options for VFE-NH
!      T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!      K. Yessad (July 2014): Move some variables.
!      J. Vivoda and P. Smolikova (Sep 2017): new options for VFE-NH
! End Modifications
!-------------------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN   , ONLY : NULOUT
USE YOMMP0   , ONLY : LOUTPUT, NPRINTLEV
USE SUVFE_HLP, ONLY : RTMIN, RTMAX, FX2T
USE YOMVERT  , ONLY : VP00
USE YOMCST   , ONLY : RPI, RD

!-------------------------------------------------------------------------------

IMPLICIT NONE

!-------------------------------------------------------------------------------

TYPE(GEOMETRY), INTENT(INOUT) :: YDGEOMETRY
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
REAL(KIND=JPRB) :: ZS(0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZP1,ZP2,ZETA_REG,ZETA_CHEB,ZWGH
REAL(KIND=JPRB) :: ZA, ZB, ZC, ZG, ZGD, ZETA, ZTOP, ZSUR
INTEGER(KIND=JPIM) :: JLEV, ITER

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "suvertfe.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUVERT',0,ZHOOK_HANDLE)
ASSOCIATE(YDVAB=>YDGEOMETRY%YRVAB,YDVETA=>YDGEOMETRY%YRVETA,YDVFE=>YDGEOMETRY%YRVFE,YDSTA=>YDGEOMETRY%YRSTA, &
 & YDLAP=>YDGEOMETRY%YRLAP,YDCSGLEG=>YDGEOMETRY%YRCSGLEG,YDCVER=>YDGEOMETRY%YRCVER, &
 & YDCSGEOM=>YDGEOMETRY%YRCSGEOM,YDCSGEOM_NB=>YDGEOMETRY%YRCSGEOM_NB,YDGSGEOM=>YDGEOMETRY%YRGSGEOM, &
 & YDGSGEOM_NB=>YDGEOMETRY%YRGSGEOM_NB, YDSPGEOM=>YDGEOMETRY%YSPGEOM)
ASSOCIATE(NFLEVG=>YDGEOMETRY%YRDIMV%NFLEVG, &
  & NIOLEVG=>YDGEOMETRY%YRDIMV%NIOLEVG,YDDIMV=>YDGEOMETRY%YRDIMV)

!*            1. SET UP VERTICAL SYSTEM: YRVETA and VFE OPERATORS.
!             ----------------------------------------------------

IF (NFLEVG > 1) THEN

  DO JLEV=0,NFLEVG
    IF (YDCVER%LREGETA) THEN
      YDVETA%VETAH(JLEV)=REAL(JLEV,JPRB)/REAL(NFLEVG,JPRB)
    ELSE
      YDVETA%VETAH(JLEV)=YDVAB%VALH(JLEV)+YDVAB%VBH(JLEV)
    ENDIF
    IF(JLEV > 0) YDVETA%VETAF(JLEV)=(YDVETA%VETAH(JLEV)+YDVETA%VETAH(JLEV-1))*0.5_JPRB
  ENDDO
  YDVETA%VETAF(0)=YDVETA%VETAH(0) ! equal to zero if the top pressure is zero.
  YDVETA%VETAF(NFLEVG+1)=YDVETA%VETAH(NFLEVG) ! equal to one.

  IF(YDCVER%LVERTFE) THEN

    RTMIN = FX2T(0.0_JPRB, 1.0_JPRB, 0.0_JPRB)
    RTMAX = FX2T(0.0_JPRB, 1.0_JPRB, 1.0_JPRB)

    WRITE(NULOUT,'("RTMIN = ",F6.3," RTMAX = ",F6.3)') RTMIN, RTMAX

    IF(NIOLEVG /= NFLEVG) THEN
      CALL ABOR1('SUVERT: HANDLING NIOLEVG /= NFLEVG WITH LVERTFE IS NOT YET CODED')
    ENDIF 

    ! half level eta definition for VFE
    IF (YDCVER%LVFE_CENTRI) THEN

       ! centripetal definition
       ZS(0) = 0.0_JPRB
 
      WRITE(NULOUT,*) "VFE eta definition = centripetal - regeta"
 
      DO JLEV = 1,NFLEVG
        ZP1 = YDVAB%VAH(JLEV-1) + YDVAB%VBH(JLEV-1) * VP00
        ZP2 = YDVAB%VAH(JLEV  ) + YDVAB%VBH(JLEV  ) * VP00
        ZS(JLEV) = ZS(JLEV-1) + (ZP2-ZP1)**YDCVER%RVFE_ALPHA
      ENDDO
 
      DO JLEV = 0,NFLEVG
        YDVETA%VFE_ETAH(JLEV) = ZS(JLEV) / ZS(NFLEVG)
      ENDDO
 
      WRITE(NULOUT,*) "VFE eta definition = chebyshev - regeta"
 
      DO JLEV = 0,NFLEVG
        ZETA_REG   = YDVETA%VFE_ETAH(JLEV)
        ZETA_CHEB  = 0.5_JPRB - 0.5_JPRB * COS(RPI * ZETA_REG)
        YDVETA%VFE_ETAH(JLEV) = (1.0_JPRB - YDCVER%RVFE_BETA) * ZETA_REG + YDCVER%RVFE_BETA * ZETA_CHEB
      ENDDO

      ! to avoid rounding arrors we prefer exact definition of boundary values
      YDVETA%VFE_ETAH(0) = 0.0_JPRB
      YDVETA%VFE_ETAH(NFLEVG) = 1.0_JPRB

    ELSEIF( YDCVER%LVFE_CHEB )THEN

      WRITE(NULOUT,*) "VFE eta definition = dpi/deta = a + b eta + c eta^2"

      ! RVFE_ALPHA - value of dpi/deta_top
      ! RVFE_BETA  - value of dpi/deta_surf

      ! RVFE_ALPHA and RVFE_BETA represents here relative
      ! density of layers in eta space.
      ! RVFE_ALPHA = 1, regular distribution
      ! RVFE_ALPHA > 1, denser close top/bottom BC
      ! RVFE_ALPHA < 1, denser towards inner domain 
      ZTOP = YDCVER%RVFE_ALPHA 
      ZSUR = YDCVER%RVFE_BETA  

      ZA = ZTOP
      ZB = -2.0_JPRB * (- 3.0_JPRB + ZSUR + 2.0_JPRB * ZTOP)
      ZC =  3.0_JPRB * (-2.0_JPRB + ZSUR + ZTOP)

      YDVETA%VFE_ETAH(0) = 0.0_JPRB

      DO JLEV = 1,NFLEVG
        ! ZP2  = YDVAB%VAH(JLEV  )/VP00 + YDVAB%VBH(JLEV  )
        ZP2  = REAL(JLEV,JPRB) / REAL(NFLEVG,JPRB)
        ZETA = YDVETA%VFE_ETAH(JLEV - 1)
        DO  ITER = 1, 10
          ZG   = ZETA * (ZA + ZETA * (ZB / 2.0_JPRB + ZETA * ZC / 3.0_JPRB)) - ZP2
          ZGD  = ZA + ZETA * (ZB + ZETA * ZC)
          ZETA = ZETA - ZG / ZGD
          WRITE(NULOUT,*) "DBG :", JLEV, ITER, ZETA, - ZG / ZGD
        ENDDO
        YDVETA%VFE_ETAH(JLEV) = ZETA
      ENDDO

      ! to avoid rounding arrors we prefer exact definition of boundary values
      YDVETA%VFE_ETAH(NFLEVG) = 1.0_JPRB

    ELSEIF( YDCVER%LVFE_REGETA )THEN

      WRITE(NULOUT,*) "VFE eta definition = i/L"

      ! regular definition
      DO JLEV=0,NFLEVG
        YDVETA%VFE_ETAH(JLEV)=REAL(JLEV,JPRB)/REAL(NFLEVG,JPRB)
      ENDDO

    ELSE

      WRITE(NULOUT,*) "VFE definition of eta = A/p_r + B (dpi/deta =~ 1)"

      ! standard definition (so called chordal);
      ! 1/p_r * dpi/deta == 1 because with this definition pi = eta * p_r
      DO JLEV=0,NFLEVG
        YDVETA%VFE_ETAH(JLEV)=YDVAB%VALH(JLEV)+YDVAB%VBH(JLEV)
      ENDDO

    ENDIF

    ! full level eta definition for VFE
    DO JLEV=1,NFLEVG
      YDVETA%VFE_ETAF(JLEV)=(YDVETA%VFE_ETAH(JLEV)+YDVETA%VFE_ETAH(JLEV-1))*0.5_JPRB
    ENDDO
    YDVETA%VFE_ETAF(0)=YDVETA%VFE_ETAH(0) ! equal to zero if the top pressure is zero.
    YDVETA%VFE_ETAF(NFLEVG+1)=YDVETA%VFE_ETAH(NFLEVG) ! equal to one.
    ! inverse of layers depth for VFE
    DO JLEV=1,NFLEVG
      YDVETA%VFE_RDETAH(JLEV)=1.0_JPRB/(YDVETA%VFE_ETAH(JLEV)-YDVETA%VFE_ETAH(JLEV-1))
    ENDDO

    IF (.NOT.YDCVER%LVFE_ECMWF) THEN
    ! initialize internal knots (first guess)
      IF( MOD(YDCVER%NVFE_ORDER,2) == 0 )THEN
        DO JLEV = 1, YDCVER%NVFE_INTERNALS
          YDVFE%VFE_KNOT(JLEV) = YDVETA%VFE_ETAF(JLEV+1)
        ENDDO
      ELSE
        DO JLEV = 1, YDCVER%NVFE_INTERNALS
          YDVFE%VFE_KNOT(JLEV) = YDVETA%VFE_ETAH(JLEV+1)
        ENDDO
      ENDIF
    ENDIF

    ! check consistency of half and full level layers
    IF (LOUTPUT.AND.NPRINTLEV >= 1) THEN
      WRITE(NULOUT,*) "HALF AND FULL LEVELS CONSISTENCY"
      DO JLEV=1,NFLEVG
        IF (YDVETA%VFE_ETAH(JLEV-1)<YDVETA%VFE_ETAF(JLEV) .AND. &
         & YDVETA%VFE_ETAF(JLEV)<YDVETA%VFE_ETAH(JLEV)) THEN
          WRITE(NULOUT,*) "LEVEL :", JLEV, " OK"
        ELSE
          WRITE(NULOUT,*) "LEVEL :", JLEV, " INCONSISTENCY; CORRECTION OF RELEVANT FULL LEVEL."
          ! YDVETA%VFE_ETAF(JLEV) = 0.5_JPRB*(YDVETA%VFE_ETAH(JLEV) + YDVETA%VFE_ETAH(JLEV-1))
       ENDIF
      ENDDO
    ENDIF

    ! * compute VFE operators (RDERI, RINTE, ...)
    CALL SUVERTFE(YDGEOMETRY)
  ENDIF

  IF (LOUTPUT.AND.NPRINTLEV >= 1) THEN
    WRITE(NULOUT,*) ''
    WRITE(NULOUT,*) ' * ETA at half levels:'
    WRITE(UNIT=NULOUT,FMT='(1X,''JLEV'',1X,5X,''     ETA            '')')
    DO JLEV=0,NFLEVG
      WRITE(UNIT=NULOUT,FMT='(1X,I3,2X,(5X,F20.10))')JLEV,YDVETA%VETAH(JLEV)
    ENDDO
    IF (YDCVER%LVERTFE) THEN
      WRITE(NULOUT,*) ''
      WRITE(NULOUT,*) ' * VFE_ETA at half levels:'
      WRITE(UNIT=NULOUT,FMT='(1X,''JLEV'',1X,5X,''     VFE_ETAH       '')')
      DO JLEV=0,NFLEVG
        WRITE(UNIT=NULOUT,FMT='(1X,I3,2X,(5X,F20.10))')JLEV,YDVETA%VFE_ETAH(JLEV)
      ENDDO
      WRITE(NULOUT,*) ''
      WRITE(NULOUT,*) ' * VFE_ETA at full levels:'
      WRITE(UNIT=NULOUT,FMT='(1X,''JLEV'',1X,5X,''     VFE_ETAF       '')')
      DO JLEV=1,NFLEVG
        WRITE(UNIT=NULOUT,FMT='(1X,I3,2X,(5X,F20.10))')JLEV,YDVETA%VFE_ETAF(JLEV)
      ENDDO
    ENDIF
  ENDIF

ELSE

  !*  Set-up vertical system for 2D models (coherence with SUOPH)
  YDVETA%VETAF(0)=0._JPRB
  YDVETA%VETAF(1)=0.5_JPRB
  YDVETA%VETAF(2)=1.0_JPRB
  DO JLEV=0,1
    YDVETA%VETAH(JLEV)=YDVAB%VALH(JLEV)+YDVAB%VBH(JLEV)
  ENDDO

ENDIF

!     ------------------------------------------------------------------

!*            2. SET UP VERTICAL SYSTEM: YRVAB.
!             ---------------------------------

IF(NFLEVG > 1 .AND. YDCVER%LVERTFE) THEN

  CALL SUVFE_ADJUST_AB(YDGEOMETRY)

ELSEIF (NFLEVG > 1 .AND. .NOT.YDCVER%LVERTFE) THEN

  ! attributes VDELA, VAF, VBF are not used in this case.
  DO JLEV=1,NFLEVG
    YDVAB%VC(JLEV)=YDVAB%VAH(JLEV)*YDVAB%VBH(JLEV-1)-YDVAB%VAH(JLEV-1)*YDVAB%VBH(JLEV)
    YDVAB%VDELB(JLEV)=YDVAB%VBH(JLEV)-YDVAB%VBH(JLEV-1)
  ENDDO

ELSE

  !*  Set-up vertical system for 2D models (coherence with SUOPH)
  ! attributes VC, VDELA, VDELB, VAF, VBF are not used in this case.
  YDVAB%VAH(0)=0._JPRB
  YDVAB%VAH(1)=0._JPRB
  YDVAB%VALH(0)=0._JPRB
  YDVAB%VALH(1)=0._JPRB
  YDVAB%VBH (0)=0._JPRB
  YDVAB%VBH (1)=1._JPRB

ENDIF

IF (LOUTPUT) THEN

  WRITE(NULOUT,*) ''
  WRITE(NULOUT,*) ' * A and B at half levels:'
  WRITE(UNIT=NULOUT,FMT='(1X,''JLEV'',1X,&
   & 5X,''     ALPHA          '',&
   & 5X,''     B              '',&
   & 5X,''     A              '')')  

  DO JLEV=0,NFLEVG
    WRITE(UNIT=NULOUT,FMT='(1X,I3,2X,3(5X,F20.10))')&
     & JLEV,YDVAB%VALH(JLEV),YDVAB%VBH(JLEV),YDVAB%VAH(JLEV)  
  ENDDO

ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUVERT',1,ZHOOK_HANDLE)
END SUBROUTINE SUVERT

!-------------------------------------------------------
! externalized VFE computation of A and B on full levels
! and VDELB and VDELA
!-------------------------------------------------------

SUBROUTINE SUVFE_ADJUST_AB(YDGEOMETRY)

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN   , ONLY : NULOUT
USE YOMMP0   , ONLY : LOUTPUT
USE YOMVERT  , ONLY : VP00

IMPLICIT NONE

TYPE(GEOMETRY), INTENT(INOUT) :: YDGEOMETRY
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
REAL(KIND=JPRB) :: ZVALH(0:YDGEOMETRY%YRDIMV%NIOLEVG)
REAL(KIND=JPRB) :: ZVBH (0:YDGEOMETRY%YRDIMV%NIOLEVG)
REAL(KIND=JPRB) :: ZIN  (0:YDGEOMETRY%YRDIMV%NFLEVG+1)
REAL(KIND=JPRB) :: ZOUT (YDGEOMETRY%YRDIMV%NFLEVG+1)

INTEGER(KIND=JPIM) :: JLEV, IREF, JIT, JITM
REAL(KIND=JPRB)    :: ZFAC, ZDETA, ZONE(YDGEOMETRY%YRDIMV%NFLEVG), ZM, ZPRES
LOGICAL            :: LLITER, LLA

#include "verdisint.intfb.h"

IF (LHOOK) CALL DR_HOOK('SUVFE_ADJUST_AB',0,ZHOOK_HANDLE)
ASSOCIATE(YDVAB=>YDGEOMETRY%YRVAB,YDVETA=>YDGEOMETRY%YRVETA,YDVFE=>YDGEOMETRY%YRVFE,&
 & YDCVER=>YDGEOMETRY%YRCVER)
ASSOCIATE(NFLEVG=>YDGEOMETRY%YRDIMV%NFLEVG, NIOLEVG=>YDGEOMETRY%YRDIMV%NIOLEVG)

DO JLEV=1,NFLEVG
  YDVAB%VDELB(JLEV)=YDVAB%VBH(JLEV)-YDVAB%VBH(JLEV-1)
  ZVBH(JLEV)=YDVAB%VDELB(JLEV)*YDVETA%VFE_RDETAH(JLEV)
  YDVAB%VDELA(JLEV)=YDVAB%VAH(JLEV)-YDVAB%VAH(JLEV-1)
ENDDO

ZIN(0) = 0.0_JPRB
ZIN(1:NFLEVG) = ZVBH(1:NFLEVG)
ZIN(NFLEVG+1) = 0.0_JPRB
CALL VERDISINT(YDVFE,YDCVER,'ITOP','11',1,1,1,NFLEVG,ZIN,ZOUT)

! NORMALIZE VDELB
DO JLEV=1,NFLEVG
  YDVAB%VDELB(JLEV)=YDVAB%VDELB(JLEV)/ZOUT(NFLEVG+1)
  ZVBH(JLEV)=YDVAB%VDELB(JLEV)*YDVETA%VFE_RDETAH(JLEV)
ENDDO

! COMPUTE FULL LEVEL VBF INTEGRATING dB/dETA (in ZVBH)
ZIN(0)=0.0_JPRB
ZIN(1:NFLEVG)=ZVBH(1:NFLEVG)
ZIN(NFLEVG+1)=0.0_JPRB
CALL VERDISINT(YDVFE,YDCVER,'ITOP','11',1,1,1,NFLEVG,ZIN,ZOUT)
DO JLEV=1,NFLEVG
  YDVAB%VBF(JLEV)=ZOUT(JLEV)
ENDDO
WRITE(NULOUT,*) ''
WRITE(NULOUT,*) ' INTEGRAL OF VDELB=',ZOUT(NFLEVG+1)

!------------------------------------
! COMPUTE VDELA
!------------------------------------
IF (YDCVER%LVFE_ECMWF.OR.(YDCVER%NVFE_INTBC<=1)) THEN
  JITM=3
  LLITER=.TRUE.
  LLA = .FALSE.
ELSE
  JITM=3
  LLITER=.FALSE.
  LLA = .TRUE.
ENDIF

IF(LLITER)THEN

  IREF = 0
  ! ITERATIVE SEARCH FOR CONSISTENT PROFILE OF A
  ! Integral(A,{eta,0,1}) must be 0
  DO JIT=1,JITM
    DO JLEV=1,NFLEVG
      ZVBH(JLEV)=YDVAB%VDELA(JLEV)*YDVETA%VFE_RDETAH(JLEV)
      IF(YDVAB%VDELA(JLEV) < 0.0_JPRB) THEN
        IF(JLEV==1) THEN
          IREF=1
        ELSEIF(YDVAB%VDELA(JLEV-1)>=0._JPRB.AND.IREF==0) THEN
          IREF=JLEV
         ENDIF
      ENDIF
    ENDDO
    IF( IREF > 0 )THEN
      ZIN(0) = 0.0_JPRB
      ZIN(1:NFLEVG) = ZVBH(1:NFLEVG)
      ZIN(NFLEVG+1) = 0.0_JPRB
      CALL VERDISINT(YDVFE,YDCVER,'ITOP','11',1,1,1,NFLEVG,ZIN,ZOUT)
      WRITE(NULOUT,*) ' ITERATION ',JIT,', INTEGRAL OF VDELA=',ZOUT(NFLEVG+1)
      ZFAC=ZOUT(IREF)/(ZOUT(IREF)-ZOUT(NFLEVG+1))
      DO JLEV=1,IREF
        YDVAB%VDELA(JLEV)=YDVAB%VDELA(JLEV)/ZFAC
        ZVBH(JLEV)=YDVAB%VDELA(JLEV)*YDVETA%VFE_RDETAH(JLEV)
      ENDDO
    ENDIF
  ENDDO

  ZIN(0) = 0.0_JPRB
  ZIN(1:NFLEVG) = ZVBH(1:NFLEVG)
  ZIN(NFLEVG+1) = 0.0_JPRB
  CALL VERDISINT(YDVFE,YDCVER,'ITOP','11',1,1,1,NFLEVG,ZIN,ZOUT)
  WRITE(NULOUT,*) ' AFTER ADJUSTING, INTEGRAL OF VDELA=',ZOUT(NFLEVG+1)
  DO JLEV=1,NFLEVG
    YDVAB%VAF(JLEV)=ZOUT(JLEV)
  ENDDO

ELSE

  !------------------------
  ! ZONE is constant for which
  ! Integral_0^1 ZONE deta = 1.0
  !------------------------

  ZIN(0)=0.0_JPRB
  IF(LLA)THEN
    ZIN = 0.0_JPRB
    DO JLEV=1,NFLEVG
      ZIN(JLEV) = YDVAB%VDELB(JLEV) * YDVETA%VFE_RDETAH(JLEV)
    ENDDO
  ELSE
    DO JLEV=1,NFLEVG
      ZIN(JLEV)=1.0_JPRB 
    ENDDO
  ENDIF
  ZIN(NFLEVG+1)=0.0_JPRB

  CALL VERDISINT(YDVFE,YDCVER,'ITOP','11',1,1,1,NFLEVG,ZIN,ZOUT)

  DO JLEV=1,NFLEVG
    ZONE(JLEV)=ZIN(JLEV) / ZOUT(NFLEVG + 1)
  ENDDO

  WRITE(NULOUT,*) ' INTEGRAL OF 1.0 = ',ZOUT(NFLEVG + 1)

  ZIN(0)=0.0_JPRB
  DO JLEV=1,NFLEVG
    ZIN(JLEV)=ZONE(JLEV)
  ENDDO
  ZIN(NFLEVG+1)=0.0_JPRB
  CALL VERDISINT(YDVFE,YDCVER,'ITOP','11',1,1,1,NFLEVG,ZIN,ZOUT)

  WRITE(NULOUT,*) ' INTEGRAL OF NORMALIZED ONE FUNCTION (shall be 1.0) = ',ZOUT(NFLEVG + 1)

  ! define quantity with vertical integral equal to 1.0
  DO JLEV=1,NFLEVG
    ZVBH(JLEV)=YDVAB%VDELA(JLEV)*YDVETA%VFE_RDETAH(JLEV)/VP00 + ZONE(JLEV)
  ENDDO

  ! integrate 
  ZIN(0) = 0.0_JPRB
  ZIN(1:NFLEVG) = ZVBH(1:NFLEVG)
  ZIN(NFLEVG+1) = 0.0_JPRB
  CALL VERDISINT(YDVFE,YDCVER,'ITOP','11',1,1,1,NFLEVG,ZIN,ZOUT)

  ! normalize arbitrary quantity
  DO JLEV=1,NFLEVG
    ZVBH(JLEV)= ZVBH(JLEV) /  ZOUT(NFLEVG + 1)  -  ZONE(JLEV)
  ENDDO

  ! integrate arbitrary quantity
  ZIN(0) = 0.0_JPRB
  ZIN(1:NFLEVG) = ZVBH(1:NFLEVG)
  ZIN(NFLEVG+1) = 0.0_JPRB
  CALL VERDISINT(YDVFE,YDCVER,'ITOP','11',1,1,1,NFLEVG,ZIN,ZOUT)

  WRITE(NULOUT,*) ' INTEGRAL OF FUNCTION OF VDELA (shall be 0.0) =',ZOUT(NFLEVG+1)

  ! redefine VDELA
  DO JLEV=1,NFLEVG
    YDVAB%VDELA(JLEV)= ZVBH(JLEV) * VP00 / YDVETA%VFE_RDETAH(JLEV)
    ZVBH(JLEV)       = YDVAB%VDELA(JLEV) * YDVETA%VFE_RDETAH(JLEV)
  ENDDO

  ! integrate to get full level values od A
  ZIN(0) = 0.0_JPRB
  ZIN(1:NFLEVG) = ZVBH(1:NFLEVG)
  ZIN(NFLEVG+1) = 0.0_JPRB
  CALL VERDISINT(YDVFE,YDCVER,'ITOP','11',1,1,1,NFLEVG,ZIN,ZOUT)
   
  WRITE(NULOUT,*) ' INTEGRAL OF VDELA AFTER NORMALISATION=',ZOUT(NFLEVG+1)

  ! integral of 1 shall be 1

  DO JLEV=1,NFLEVG
    YDVAB%VAF(JLEV)=ZOUT(JLEV)
  ENDDO
ENDIF

DO JLEV=1,NFLEVG
  YDVAB%VC(JLEV)=YDVAB%VAH(JLEV)*YDVAB%VBH(JLEV-1)-YDVAB%VAH(JLEV-1)*YDVAB%VBH(JLEV)
ENDDO

IF (LOUTPUT) THEN
  WRITE(NULOUT,*) ''
  WRITE(NULOUT,*) ' * A and B and m at full levels:'
  WRITE(UNIT=NULOUT,FMT='(1X,''JLEV'',1X,&
   & 2X,''     A              '',&
   & 2X,''     B              '',&
   & 2X,''     eta            '',&
   & 2X,''     deta           '',&
   & 2X,''     dA/deta        '',&
   & 2X,''     dB/deta        '',&
   & 2X,''     m/p_r = 1/p_r * dpi/deta'')')

  DO JLEV=1,NFLEVG
    ZPRES = VP00
    ZM = (YDVAB%VDELA(JLEV)/ZPRES + YDVAB%VDELB(JLEV)) * YDVETA%VFE_RDETAH(JLEV)
    WRITE(UNIT=NULOUT,FMT='(1X,I3,2X,7(2X,F20.10))')&
     & JLEV,YDVAB%VAF(JLEV),YDVAB%VBF(JLEV),YDVETA%VFE_ETAF(JLEV), &
     & 1.0_JPRB / YDVETA%VFE_RDETAH(JLEV), &
     & YDVAB%VDELA(JLEV)/ZPRES * YDVETA%VFE_RDETAH(JLEV), &
     & YDVAB%VDELB(JLEV)/ZPRES * YDVETA%VFE_RDETAH(JLEV), ZM
  ENDDO

  ! DO JLEV=1,NFLEVG
  !   ZPRES = 110000.0_JPRB
  !   ZM = (YDVAB%VDELA(JLEV)/ZPRES + YDVAB%VDELB(JLEV)) * YDVETA%VFE_RDETAH(JLEV)
  !   WRITE(UNIT=NULOUT,FMT='(1X,I3,2X,4(5X,F20.10))')&
  !    & JLEV,YDVAB%VAF(JLEV),YDVAB%VBF(JLEV),YDVETA%VFE_ETAF(JLEV), ZM
  ! ENDDO

ENDIF

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUVFE_ADJUST_AB',1,ZHOOK_HANDLE)
END SUBROUTINE SUVFE_ADJUST_AB

