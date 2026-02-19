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

SUBROUTINE LAINOR2(YDGEOMETRY,YDDYN,YDDYNA,YDSL,YDSLINT,KPROMA,KSTART,KPROF,KFLDN,KFLDX,&
 & KSTABUF,&
 & KVLAG,KWLAG,LDPLANE,&
 & LDQMHW,LDQMHP,&
 & PUL9,PVL9,PCL9,PUL0,PVL0,PCL0,&
 & KIBL,&
 & KSTTYP,PDSTRET,PC2M1,PC2P1,P4JP,PI,PDEPI,PIS2,&
 & PLOCEN,PMUCEN,PLSDEPI,PLATI,&
 & PCOSCO,PSINCO,PSINLA,PCOPHI,&
 & PLON,PLAT,&
 & PO,PQ,&
 & PU9,PV9,PC9,&
 & PZUZ9,PZVZ9)  

!**** *LAINOR2- semi-LAgrangian scheme:
!               INterpolations corresponding to the ORigin.
!               2D model.

!     Purpose.
!     --------

!      Once the medium point M is known, this subroutine determines
!      the origin O in computing a bow on a great circle and calls
!      a subroutine that performs interpolations at the origin point
!      of the SL trajectory.

!**   Interface.
!     ----------
!        *CALL* *LAINOR2(...)

!        Explicit arguments :
!        --------------------

!        INPUT:
!          YDSL    - SL_STRUCT definition.
!          KPROMA  - horizontal dimension of point handled by this call.
!          KSTART  - first element of arrays where
!                    computations are performed.
!          KPROF   - depth of work.
!          KFLDN   - number of the first field.
!          KFLDX   - number of the last field.
!          KSTABUF - for a latitude IGL, KSTABUF(IGL) is the
!                    address of the element corresponding to
!                    (ILON=1,IGL) in the KPROMA arrays.
!          KVLAG   - switch for discretization of continuity equation.
!          KWLAG   - switch for discretization of momentum equation.
!          LDPLANE - switch: .T. = plane geometry; .F. = spherical geometry.
!          LDQMH[X]- Use quasi-monotone interpolation in the horizontal
!                    for field X (Wind, eq. height).
!          PUL9    - Quantity to interpolate at O in the U-wind equation.
!          PVL9    - Quantity to interpolate at O in the V-wind equation.
!          PCL9    - Quantity to interpolate at O in the continuity equation.
!          PUL0    - Second quantity to interpolate at O
!                    (if NWLAG=3) in the U-component of the wind equation.
!          PVL0    - Second quantity to interpolate at O
!                    (if NWLAG=3) in the V-component of the wind equation.
!          PCL0    - Second quantity to interpolate at O
!                    (if NVLAG=3) in the continuity equation.
!          KIBL    - index into YRGSGEOM/YRCSGEOM in YDGEOMETRY 
!          KSTTYP  - 1: Not tilted pole;  2: Tilted pole.
!          PDSTRET - 2*c (where c is the stretching factor).
!          PC2M1   - c*c-1.
!          PC2P1   - c*c+1.
!          P4JP    - Approximative inverse of the differences of latitudes.
!          PI      - number PI
!          PDEPI   - 2 * PI
!          PIS2    - PI / 2
!          PLOCEN  - geographical longitude of the stretching pole.
!          PMUCEN  - sine of the geographical latitude of the pole.
!          PLSDEPI - (Number of points by latitude) / (2 * PI) .
!          PLATI   - latitude on the computational sphere.

!        INPUT/OUTPUT:
!          PCOSCO  - cos(Longitude-Longitude(grid-point))*cos(Latitude)
!                    of the interpolation point (geographical longitude
!                    and latitude).
!          PSINCO  - sin(Longitude-Longitude(grid-point))*cos(Latitude)
!                    of the interpolation point (geographical longitude
!                    and latitude).
!          PSINLA  - sine of the interpolation point geographical latitude.
!          PCOPHI  - cosine of the geographical angle between the
!                    interpolation point and the grid-point.

!        OUTPUT:
!          PLON    - computational sphere longitude of interpolation point.
!          PLAT    - computational sphere latitude of interpolation point.
!          PO      - first element of the wind displacement matrix.
!          PQ      - second element of the wind displacement matrix.
!          PU9     - Interpolated quantity at O in the U-wind equation.
!          PV9     - Interpolated quantity at O in the V-wind equation.
!          PC9     - Interpolated quantity at O in the continuity equation.
!          PZUZ9   - Interpolated quantity at O in the U-wind equation
!                    in the SL2TL scheme for case L2TLFF=.T.
!          PZVZ9   - Interpolated quantity at O in the V-wind equation
!                    in the SL2TL scheme for case L2TLFF=.T.

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation about semi-Lagrangian scheme.

!     Externals.
!     ----------
!        Calls LARCIN2.
!        Is called by LADINE (2D shallow water model).

!     Reference.
!     ----------

!     Author.
!     -------
!      C. Temperton (ECMWF), by simplification of LAINOR
!      Original : 98/01/19

!     Modifications.
!     --------------
!      Modified 01-08-30 by K. YESSAD: pruning and some other cleanings.
!      01-Oct-2003 M. Hamrud   CY28 Cleaning
!      30-Jul-2008 J. Masek    Dataflow for new SLHD interpolators.
!      K. Yessad Nov 2008: rationalisation of dummy argument interfaces
!      G.Mozdzynski (Jan 2011): OOPS cleaning, use of derived type SL_STRUCT
!      G.Mozdzynski (Feb 2011): OOPS cleaning, use of derived types TGSGEOM and TCSGEOM
!      G. Mozdzynski (May 2012): further cleaning
!      T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
! End Modifications
!-------------------------------------------------------------------------------

USE YOMDYN       , ONLY : TDYN
USE YOMDYNA      , ONLY : TDYNA
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMSLINT     , ONLY : TSLINT
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE EINT_MOD     , ONLY : SL_STRUCT

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TDYN)        ,INTENT(IN)    :: YDDYN
TYPE(TDYNA)       ,INTENT(IN)    :: YDDYNA
TYPE(SL_STRUCT)   ,INTENT(INOUT) :: YDSL
TYPE(TSLINT)      ,INTENT(IN)    :: YDSLINT
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDN 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDX 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTABUF(YDSL%NDGSAH:YDSL%NDGENH) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KVLAG 
INTEGER(KIND=JPIM),INTENT(IN)    :: KWLAG 
LOGICAL           ,INTENT(IN)    :: LDPLANE 
LOGICAL           ,INTENT(IN)    :: LDQMHW 
LOGICAL           ,INTENT(IN)    :: LDQMHP 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUL9(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVL9(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCL9(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUL0(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVL0(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCL0(YDSL%NASLB1,KFLDN:KFLDX) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIBL
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTTYP 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDSTRET
REAL(KIND=JPRB)   ,INTENT(IN)    :: PC2M1 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PC2P1
REAL(KIND=JPRB)   ,INTENT(IN)    :: P4JP 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PI 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDEPI 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PIS2 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLOCEN
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMUCEN
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLSDEPI(YDSL%NDGSAH:YDSL%NDGENH) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLATI(YDSL%NDGSAH:YDSL%NDGENH) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCOSCO(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSINCO(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSINLA(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCOPHI(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PLON(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLAT(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PO(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQ(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PU9(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PV9(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PC9(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PZUZ9(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PZVZ9(KPROMA) 

!     -----------------------------------------------------------------

INTEGER(KIND=JPIM) :: IWISA
INTEGER(KIND=JPIM) :: IXLAG(6)

INTEGER(KIND=JPIM) :: IROT, JROF

REAL(KIND=JPRB) :: ZUF0(1), ZURL0(1), ZVF0(1), ZVRL0(1)
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "larcin2.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LAINOR2',0,ZHOOK_HANDLE)
ASSOCIATE(YDGSGEOM=>YDGEOMETRY%YRGSGEOM(KIBL))
!     ------------------------------------------------------------------

!*       1.    COMPUTATION OF THE ORIGIN POINT "O" COORDINATES.
!              ------------------------------------------------

IXLAG(1)=KVLAG
IXLAG(4)=KWLAG

IF(LDPLANE) THEN
  ! not coded
ELSE

  DO JROF=KSTART,KPROF

!   * computations on horizontal plans.

    PSINLA(JROF)=2.0_JPRB*PCOPHI(JROF)*PSINLA(JROF)-YDGSGEOM%GEMU (JROF)
    PSINCO(JROF)=2.0_JPRB*PCOPHI(JROF)*PSINCO(JROF)
    PCOSCO(JROF)=2.0_JPRB*PCOPHI(JROF)*PCOSCO(JROF)-YDGSGEOM%GSQM2(JROF)

!   * Computation of the angle <0,G> before calling LARCIN2:

    PCOPHI(JROF)=2.0_JPRB*PCOPHI(JROF)*PCOPHI(JROF)-1.0_JPRB

  ENDDO

ENDIF

!     ------------------------------------------------------------------

!*       2.    CALL OF LARCIN2 FOR INTERPOLATIONS.
!              -----------------------------------

IWISA=3
IROT=1

!     Warning: ZURL0,ZVRL0,ZUF0,ZVF0 are not
!              used in this call to LARCIN2.

CALL LARCIN2(YDGEOMETRY,YDDYN,YDDYNA,YDSL,YDSLINT,KPROMA,KSTART,KPROF,KFLDN,KFLDX,&
 & KSTABUF,IWISA,IXLAG,IROT,&
 & LDPLANE,LDQMHW,LDQMHP,&
 & KSTTYP,PDSTRET,PC2M1,PC2P1,P4JP,PI,PDEPI,PIS2,&
 & PLOCEN,PMUCEN,PLSDEPI,PLATI,&
 & KIBL,&
 & PCOSCO,PSINCO,PSINLA,PCOPHI,&
 & ZURL0,ZVRL0,&
 & PUL9,PVL9,PCL9,PUL0,PVL0,PCL0,&
 & PLON,PLAT,&
 & PO,PQ,&
 & ZUF0,ZVF0,&
 & PU9,PV9,PC9,&
 & PZUZ9,PZVZ9)

!     ------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('LAINOR2',1,ZHOOK_HANDLE)
END SUBROUTINE LAINOR2

