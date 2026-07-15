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
!      K. Yessad (Aug 2018): add VRATH and VRATF in structure TVAB
!      J. Vivoda and P. Smolikova (Sep 2020): suvfe_adjust_ab rewrite for VFE.
! End Modifications
!-------------------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM, JPRB
USE SUVFE_HLP, ONLY : RTMIN, RTMAX, FX2T
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN   , ONLY : NULOUT
USE YOMCST   , ONLY : RPI
USE YOMMP0   , ONLY : LOUTPUT
USE YOMVERT  , ONLY : VP00

!-------------------------------------------------------------------------------

IMPLICIT NONE

!-------------------------------------------------------------------------------

TYPE(GEOMETRY), INTENT(INOUT) :: YDGEOMETRY
INTEGER(KIND=JPIM) :: JLEV,ITER
REAL(KIND=JPRB) :: ZS(0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZP1,ZP2,ZETA_REG,ZETA_COS,ZETA
REAL(KIND=JPRB) :: ZA,ZB,ZC,ZG,ZGD
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "suvertfe.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUVERT',0,ZHOOK_HANDLE)
ASSOCIATE(YDVAB=>YDGEOMETRY%YRVAB,YDVETA=>YDGEOMETRY%YRVETA,YDVFE=>YDGEOMETRY%YRVFE,&
 & YDCVER=>YDGEOMETRY%YRCVER)
ASSOCIATE(NFLEVG=>YDGEOMETRY%YRDIMV%NFLEVG, &
  & NIOLEVG=>YDGEOMETRY%YRDIMV%NIOLEVG,YDDIMV=>YDGEOMETRY%YRDIMV,&
  & REXP_VRAT=>YDCVER%REXP_VRAT)

!*            1. SET UP VERTICAL SYSTEM: YRVETA and VFE OPERATORS.
!             ----------------------------------------------------

IF (NFLEVG > 1) THEN

  IF (YDCVER%LREGETA) THEN
    DO JLEV=0,NFLEVG
      YDVETA%VETAH(JLEV)=REAL(JLEV,JPRB)/REAL(NFLEVG,JPRB)
    ENDDO
  ELSE
    DO JLEV=0,NFLEVG
      YDVETA%VETAH(JLEV)=YDVAB%VALH(JLEV)+YDVAB%VBH(JLEV)
    ENDDO
  ENDIF

  DO JLEV=1,NFLEVG
    YDVETA%VETAF(JLEV)=(YDVETA%VETAH(JLEV)+YDVETA%VETAH(JLEV-1))*0.5_JPRB
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

    ! Half level eta definition for VFE
    IF (YDCVER%CVFE_ETAH=='REGETA') THEN

      ! regular definition
      WRITE(NULOUT,*) "VFE eta definition = i/L"
      DO JLEV=0,NFLEVG
        YDVETA%VFE_ETAH(JLEV)=REAL(JLEV,JPRB)/REAL(NFLEVG,JPRB)
      ENDDO

    ELSEIF (YDCVER%CVFE_ETAH=='MCONST') THEN

      ! centripetal definition
      WRITE(NULOUT,*) "VFE eta definition = centripetal - regeta"
      ZS(0) = 0.0_JPRB
      DO JLEV = 1,NFLEVG
        ZP1 = YDVAB%VALH(JLEV-1) + YDVAB%VBH(JLEV-1)
        ZP2 = YDVAB%VALH(JLEV  ) + YDVAB%VBH(JLEV  )
        ZS(JLEV) = ZS(JLEV-1) + (ZP2-ZP1)**YDCVER%RVFE_ALPHA
      ENDDO

      WRITE(NULOUT,*) "VFE eta definition: density increased on boundaries &
       & acc.to RVFE_BETA"
      DO JLEV = 0,NFLEVG
        ZETA_REG   = ZS(JLEV)/ZS(NFLEVG)
        ZETA_COS  = 0.5_JPRB - 0.5_JPRB * COS(RPI * ZETA_REG)
        YDVETA%VFE_ETAH(JLEV) = (1.0_JPRB - YDCVER%RVFE_BETA) * ZETA_REG &
         & + YDCVER%RVFE_BETA * ZETA_COS
      ENDDO
      ! to avoid rounding arrors we prefer exact definition of boundary values
      YDVETA%VFE_ETAH(0) = 0.0_JPRB
      YDVETA%VFE_ETAH(NFLEVG) = 1.0_JPRB

    ELSEIF (YDCVER%CVFE_ETAH=='MCNTRL') THEN

      WRITE(NULOUT,*) "VFE eta definition = dpi/deta = a + b eta + c eta^2"
      ! RVFE_ALPHA (dpi/deta_top) and RVFE_BETA (dpi/deta_surf) represent here 
      !   relative density of layers in eta space:
      !   RVFE_ALPHA/BETA = 0, regular distribution
      !   RVFE_ALPHA/BETA > 0, denser close to top/bottom boundary
      !   RVFE_ALPHA/BETA < 0, denser towards inner domain 
      ZA = YDCVER%RVFE_ALPHA + 1.0_JPRB
      ZB = -2.0_JPRB * (2.0_JPRB * YDCVER%RVFE_ALPHA + YDCVER%RVFE_BETA)
      ZC =  3.0_JPRB * (YDCVER%RVFE_ALPHA + YDCVER%RVFE_BETA)

      YDVETA%VFE_ETAH(0) = 0.0_JPRB
      DO JLEV = 1,NFLEVG
        ZP2  = REAL(JLEV,JPRB)/REAL(NFLEVG,JPRB)
        ZETA = YDVETA%VFE_ETAH(JLEV-1)
        DO  ITER = 1, 10
          ZG   = ZETA * (ZA + ZETA * (ZB / 2.0_JPRB + ZETA * ZC / 3.0_JPRB)) - ZP2
          ZGD  = ZA + ZETA * (ZB + ZETA * ZC)
          ZETA = ZETA - ZG / ZGD
        ENDDO
        YDVETA%VFE_ETAH(JLEV) = ZETA
      ENDDO
      ! to avoid rounding arrors we prefer exact definition of boundary values
      YDVETA%VFE_ETAH(NFLEVG) = 1.0_JPRB

    ELSE

      ! Standard definition (so called chordal);
      ! 1/pref * dpi/deta == 1 by definition since A/pref + B = pi/pref = eta
      WRITE(NULOUT,*) "VFE eta definition = A/p_r + B (dpi/deta =~ 1)"
      DO JLEV=0,NFLEVG
        YDVETA%VFE_ETAH(JLEV)=YDVAB%VALH(JLEV)+YDVAB%VBH(JLEV)
      ENDDO
    ENDIF

    ! full level eta definition and inverse of layers depth for VFE
    DO JLEV=1,NFLEVG
      YDVETA%VFE_ETAF(JLEV)=(YDVETA%VFE_ETAH(JLEV)+YDVETA%VFE_ETAH(JLEV-1))*0.5_JPRB
      YDVETA%VFE_RDETAH(JLEV)=1.0_JPRB/(YDVETA%VFE_ETAH(JLEV)-YDVETA%VFE_ETAH(JLEV-1))
    ENDDO
    YDVETA%VFE_ETAF(0)=YDVETA%VFE_ETAH(0) ! equal to zero if the top pressure is zero.
    YDVETA%VFE_ETAF(NFLEVG+1)=YDVETA%VFE_ETAH(NFLEVG) ! equal to one.

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

    IF (YDCVER%LVFE_VERBOSE) THEN
      WRITE(NULOUT,*) ''
      WRITE(NULOUT,*) ' * VFE_ETA at half and full levels:'
      WRITE(UNIT=NULOUT,FMT='(1X,''JLEV'',1X,5X,''     VFE_ETAH '', &
       &                                  1X,5X,''     VFE_ETAF '')')
      WRITE(UNIT=NULOUT,FMT='(1X,I3,2X,1X,F18.10)') 0,YDVETA%VFE_ETAH(0)
      DO JLEV=1,NFLEVG
        WRITE(UNIT=NULOUT,FMT='(1X,I3,2X,2(1X,F18.10))') &
         & JLEV,YDVETA%VFE_ETAH(JLEV),YDVETA%VFE_ETAF(JLEV)
      ENDDO
    ENDIF

    ! * compute VFE operators (RDERI, RINTE, ...)
    CALL SUVERTFE(YDGEOMETRY)

  ENDIF

ELSE

  !*  Set-up vertical system for 2D models (coherence with SUOPH)
  YDVETA%VETAF(0)=0._JPRB
  YDVETA%VETAF(1)=0.5_JPRB
  YDVETA%VETAF(2)=1.0_JPRB

  YDVETA%VETAH(0)=YDVAB%VALH(0)+YDVAB%VBH(0)
  YDVETA%VETAH(1)=YDVAB%VALH(1)+YDVAB%VBH(1)

ENDIF

!     ------------------------------------------------------------------

!*            2. SET UP VERTICAL SYSTEM: YRVAB.
!             ---------------------------------

IF (NFLEVG > 1) THEN 
  IF (YDCVER%LVERTFE) THEN
    CALL SUVFE_ADJUST_AB(YDGEOMETRY)

    ! YDVAB%VRATH and YDVAB%VRATF:
    YDVAB%VRATH(0)=0.0_JPRB
    DO JLEV=1,NFLEVG
      YDVAB%VRATH(JLEV)=(YDVAB%VBH(JLEV)*VP00/(YDVAB%VAH(JLEV)+YDVAB%VBH(JLEV)*VP00))**REXP_VRAT
      YDVAB%VRATF(JLEV)=(YDVAB%VBF(JLEV)*VP00/(YDVAB%VAF(JLEV)+YDVAB%VBF(JLEV)*VP00))**REXP_VRAT
    ENDDO
    YDVAB%VRATH(1)=MIN(0.5_JPRB,YDVAB%VRATH(1))
    YDVAB%VRATF(1)=MIN(0.25_JPRB,YDVAB%VRATF(1))
    IF (NFLEVG > 2) YDVAB%VRATF(2)=MIN(0.75_JPRB,YDVAB%VRATF(2))
  ELSE
    ! attributes VDELA, VAF, VBF are not used in this case.
    DO JLEV=1,NFLEVG
      YDVAB%VC(JLEV)=YDVAB%VAH(JLEV)*YDVAB%VBH(JLEV-1)-YDVAB%VAH(JLEV-1)*YDVAB%VBH(JLEV)
      YDVAB%VDELB(JLEV)=YDVAB%VBH(JLEV)-YDVAB%VBH(JLEV-1)
    ENDDO

    YDVAB%VRATH(0)=0.0_JPRB
    DO JLEV=1,NFLEVG
      YDVAB%VRATH(JLEV)=(YDVAB%VBH(JLEV)*VP00/(YDVAB%VAH(JLEV)+YDVAB%VBH(JLEV)*VP00))**REXP_VRAT
      ZB = YDVAB%VBH(JLEV)+YDVAB%VBH(JLEV-1)
      YDVAB%VRATF(JLEV)=(ZB*VP00/(YDVAB%VAH(JLEV)+YDVAB%VAH(JLEV-1)+ZB*VP00))**REXP_VRAT
    ENDDO
    YDVAB%VRATH(1)=MIN(0.5_JPRB,YDVAB%VRATH(1))
    YDVAB%VRATF(1)=MIN(0.25_JPRB,YDVAB%VRATF(1))
    IF (NFLEVG>2) YDVAB%VRATF(2)=MIN(0.75_JPRB,YDVAB%VRATF(2))
  ENDIF
ELSE
  !*  Set-up vertical system for 2D models (coherence with SUOPH)
  ! attributes VC, VDELA, VDELB, VAF, VBF are not used in this case.
  YDVAB%VAH(0)=0._JPRB
  YDVAB%VAH(1)=0._JPRB
  YDVAB%VALH(0)=0._JPRB
  YDVAB%VALH(1)=0._JPRB
  YDVAB%VBH (0)=0._JPRB
  YDVAB%VBH (1)=1._JPRB
  YDVAB%VRATH(0)=0.0_JPRB
  YDVAB%VRATH(1)=1.0_JPRB
  YDVAB%VRATF(1)=0.5_JPRB
ENDIF

IF (LOUTPUT) THEN
  WRITE(NULOUT,*) ''
  WRITE(NULOUT,*) ' * A and B at half levels:'
  WRITE(UNIT=NULOUT,FMT="('JLEV',5(X,A15))") "A","B","ETA","ALPHA","S"

  DO JLEV=0,NFLEVG
    WRITE(UNIT=NULOUT,FMT='(I4,X,F15.6,4(X,F15.12))') JLEV,YDVAB%VAH(JLEV),YDVAB%VBH(JLEV),&
     & YDVETA%VETAH(JLEV),YDVAB%VALH(JLEV),YDVAB%VRATH(JLEV)  
  ENDDO
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUVERT',1,ZHOOK_HANDLE)
CONTAINS

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
REAL(KIND=JPRB) :: ZIN(0:YDGEOMETRY%YRDIMV%NFLEVG+1)
REAL(KIND=JPRB) :: ZOUT(YDGEOMETRY%YRDIMV%NFLEVG+1)
REAL(KIND=JPRB) :: ZONE(YDGEOMETRY%YRDIMV%NFLEVG)

INTEGER(KIND=JPIM) :: JLEV, JIT, JITM
REAL(KIND=JPRB)    :: ZFAC, ZM

#include "verdisint.intfb.h"

IF (LHOOK) CALL DR_HOOK('SUVERT:SUVFE_ADJUST_AB',0,ZHOOK_HANDLE)
ASSOCIATE(YDVAB=>YDGEOMETRY%YRVAB,YDVETA=>YDGEOMETRY%YRVETA,YDVFE=>YDGEOMETRY%YRVFE,YDCVER=>YDGEOMETRY%YRCVER)
ASSOCIATE(NFLEVG=>YDGEOMETRY%YRDIMV%NFLEVG)

DO JLEV=1,NFLEVG
  YDVAB%VDELB(JLEV)=YDVAB%VBH(JLEV)-YDVAB%VBH(JLEV-1)
  ZIN(JLEV)=YDVAB%VDELB(JLEV)*YDVETA%VFE_RDETAH(JLEV)
  YDVAB%VDELA(JLEV)=YDVAB%VAH(JLEV)-YDVAB%VAH(JLEV-1)
ENDDO
ZIN(0) = 0.0_JPRB
ZIN(NFLEVG+1) = 0.0_JPRB
CALL VERDISINT(YDVFE,YDCVER,'ITOP','11',1,1,1,NFLEVG,ZIN,ZOUT)

! NORMALIZE VDELB
DO JLEV=1,NFLEVG
  YDVAB%VDELB(JLEV)=YDVAB%VDELB(JLEV)/ZOUT(NFLEVG+1)
  ZIN(JLEV)=YDVAB%VDELB(JLEV)*YDVETA%VFE_RDETAH(JLEV)
ENDDO

! COMPUTE FULL LEVEL VBF INTEGRATING dB/dETA (in ZIN)
ZIN(0)=0.0_JPRB
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
IF (YDCVER%LVFE_ECMWF) THEN
  JITM=3

  ! ITERATIVE SEARCH FOR CONSISTENT PROFILE OF A
  ! Integral(A,{eta,0,1}) must be 0

  ZIN(1:NFLEVG)=YDVAB%VDELA(1:NFLEVG)*YDVETA%VFE_RDETAH(1:NFLEVG)
  ZIN(0) = 0.0_JPRB
  ZIN(NFLEVG+1) = 0.0_JPRB
  DO JIT=1,JITM
    DO JLEV=1,NFLEVG
      IF (YDVAB%VDELA(JLEV) < 0.0_JPRB) EXIT
    ENDDO
    IF (JLEV > NFLEVG) EXIT

    CALL VERDISINT(YDVFE,YDCVER,'ITOP','11',1,1,1,NFLEVG,ZIN,ZOUT)
    WRITE(NULOUT,*) " ITER/JLEV/INTEGRAL OF VDELA:",JIT,JLEV,ZOUT(NFLEVG+1)

    ZFAC=ZOUT(JLEV)/(ZOUT(JLEV)-ZOUT(NFLEVG+1))
    YDVAB%VDELA(1:JLEV)=YDVAB%VDELA(1:JLEV)/ZFAC
    ZIN(1:JLEV)=YDVAB%VDELA(1:JLEV)*YDVETA%VFE_RDETAH(1:JLEV)
  ENDDO

  ZIN(0) = 0.0_JPRB
  ZIN(NFLEVG+1) = 0.0_JPRB
  CALL VERDISINT(YDVFE,YDCVER,'ITOP','11',1,1,1,NFLEVG,ZIN,ZOUT)
  WRITE(NULOUT,*) ' AFTER ADJUSTING, INTEGRAL OF VDELA=',ZOUT(NFLEVG+1)

  DO JLEV=1,NFLEVG
    YDVAB%VAF(JLEV)=ZOUT(JLEV)
    YDVAB%VC(JLEV)=YDVAB%VAH(JLEV)*YDVAB%VBH(JLEV-1)-YDVAB%VAH(JLEV-1)*YDVAB%VBH(JLEV)
  ENDDO

ELSE

  !------------------------
  ! ZONE is constant for which
  ! Integral_0^1 ZONE deta = 1.0
  !------------------------

  DO JLEV=1,NFLEVG
    ZIN(JLEV) = YDVAB%VDELB(JLEV) * YDVETA%VFE_RDETAH(JLEV)
  ENDDO
  ZIN(0)=0.0_JPRB
  ZIN(NFLEVG+1)=0.0_JPRB
  CALL VERDISINT(YDVFE,YDCVER,'ITOP','11',1,1,1,NFLEVG,ZIN,ZOUT)
  WRITE(NULOUT,*) ' INTEGRAL OF 1.0 = ',ZOUT(NFLEVG+1)

  DO JLEV=1,NFLEVG
    ZONE(JLEV)=ZIN(JLEV)/ZOUT(NFLEVG+1)
    ZIN(JLEV)=ZONE(JLEV)
  ENDDO
  ZIN(0)=0.0_JPRB
  ZIN(NFLEVG+1)=0.0_JPRB
  CALL VERDISINT(YDVFE,YDCVER,'ITOP','11',1,1,1,NFLEVG,ZIN,ZOUT)
  WRITE(NULOUT,*) ' INTEGRAL OF NORMALIZED ONE FUNCTION (shall be 1.0) = ',ZOUT(NFLEVG+1)

  ! define quantity with vertical integral equal to 1.0
  DO JLEV=1,NFLEVG
    ZIN(JLEV)=YDVAB%VDELA(JLEV)*YDVETA%VFE_RDETAH(JLEV)/VP00 + ZONE(JLEV)
  ENDDO

  ! integrate 
  ZIN(0) = 0.0_JPRB
  ZIN(NFLEVG+1) = 0.0_JPRB
  CALL VERDISINT(YDVFE,YDCVER,'ITOP','11',1,1,1,NFLEVG,ZIN,ZOUT)

  ! normalize arbitrary quantity
  DO JLEV=1,NFLEVG
    ZIN(JLEV) = ZIN(JLEV)/ZOUT(NFLEVG+1) - ZONE(JLEV)
  ENDDO

  ! integrate arbitrary quantity
  ZIN(0) = 0.0_JPRB
  ZIN(NFLEVG+1) = 0.0_JPRB
  CALL VERDISINT(YDVFE,YDCVER,'ITOP','11',1,1,1,NFLEVG,ZIN,ZOUT)
  WRITE(NULOUT,*) ' INTEGRAL OF FUNCTION OF VDELA (shall be 0.0) =',ZOUT(NFLEVG+1)

  ! redefine VDELA
  DO JLEV=1,NFLEVG
    YDVAB%VDELA(JLEV) = ZIN(JLEV) * VP00 / YDVETA%VFE_RDETAH(JLEV)
    ZIN(JLEV)         = YDVAB%VDELA(JLEV) * YDVETA%VFE_RDETAH(JLEV)
  ENDDO

  ! integrate to get full level values of A
  ZIN(0) = 0.0_JPRB
  ZIN(NFLEVG+1) = 0.0_JPRB
  CALL VERDISINT(YDVFE,YDCVER,'ITOP','11',1,1,1,NFLEVG,ZIN,ZOUT)
  WRITE(NULOUT,*) ' AFTER ADJUSTING, INTEGRAL OF VDELA=',ZOUT(NFLEVG+1)

  DO JLEV=1,NFLEVG
    YDVAB%VAF(JLEV)=ZOUT(JLEV)
    YDVAB%VC(JLEV)=YDVAB%VAH(JLEV)*YDVAB%VBH(JLEV-1)-YDVAB%VAH(JLEV-1)*YDVAB%VBH(JLEV)
  ENDDO
ENDIF

IF (LOUTPUT) THEN
  IF (YDCVER%LVFE_VERBOSE) THEN
    WRITE(NULOUT,*) ''
    WRITE(NULOUT,*) ' * A, B and some more things at full levels:'
    WRITE(UNIT=NULOUT,FMT="('JLEV',7(X,A18))") "A","B","eta","deta","dA/deta",&
     & "dB/deta","m/p_r = 1/p_r * dpi/deta"

    DO JLEV=1,NFLEVG
      ZM = (YDVAB%VDELA(JLEV)/VP00 + YDVAB%VDELB(JLEV)) * YDVETA%VFE_RDETAH(JLEV)
      WRITE(UNIT=NULOUT,FMT='(1X,I3,7(X,F18.10))') &
       & JLEV,YDVAB%VAF(JLEV),YDVAB%VBF(JLEV),YDVETA%VFE_ETAF(JLEV), &
       & 1.0_JPRB / YDVETA%VFE_RDETAH(JLEV), &
       & YDVAB%VDELA(JLEV)/VP00 * YDVETA%VFE_RDETAH(JLEV), &
       & YDVAB%VDELB(JLEV)/VP00 * YDVETA%VFE_RDETAH(JLEV), ZM
    ENDDO
  ELSE
    WRITE(NULOUT,*) ''
    WRITE(NULOUT,*) ' * A and B at full levels:'
    WRITE(UNIT=NULOUT,FMT="('JLEV',2(X,A18))") "A","B"

    DO JLEV=1,NFLEVG
      WRITE(UNIT=NULOUT,FMT='(1X,I3,2(X,F18.10))') JLEV,YDVAB%VAF(JLEV),YDVAB%VBF(JLEV)
    ENDDO
  ENDIF
ENDIF

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUVERT:SUVFE_ADJUST_AB',1,ZHOOK_HANDLE)
END SUBROUTINE SUVFE_ADJUST_AB

END SUBROUTINE SUVERT
