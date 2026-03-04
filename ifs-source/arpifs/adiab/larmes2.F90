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

SUBROUTINE LARMES2(YDGEOMETRY,YDDYN,YDDYNA,YDSL,YDSLINT,KPROMA,KSTART,KPROF,KFLDN,KFLDX,&
 & KULOUT,KSTABUF,&
 & KVLAG,KWLAG,&
 & KITER,PVMAX1,PVMAX2,LDPLANE,LDSETTLS,LDELTRA,&
 & PURL0,PVRL0,PURL,PVRL,&
 & KSTTYP,PDSTRET,PC2M1,PC2P1,P4JP,PI,PDEPI,PIS2,&
 & PDTSA,PDTS62,PDTS22,&
 & PLOCEN,PMUCEN,PLSDEPI,PLATI,&
 & KIBL,&
 & PLON,PLAT,&
 & PCOSCO,PSINCO,PSINLA,PCOPHI,&
 & PR,PS)  

!**** *LARMES2 - semi-LAgrangian scheme:
!                Research of the MEdium point on the Sphere (2-d model)

!     Purpose.
!     --------

!     The computation of the location of the medium point I of
!     the lagrangian trajectory is performed by an iterative
!     method described by Robert and adapted to the sphere by M. Rochas.
!     Trajectories are great circles.

!**   Interface.
!     ----------
!        *CALL* *LARMES2(...)

!        Explicit arguments :
!        --------------------

!        INPUT:
!          YDSL    - SL_STRUCT definition.
!          KPROMA  - horizontal dimension.
!          KSTART  - first element of arrays where
!                    computations are performed.
!          KPROF   - depth of work.
!          KFLDN   - number of the first field.
!          KFLDX   - number of the last field.
!          KULOUT  - output logical unit.
!          KSTABUF - for a latitude IGL, KSTABUF(IGL) is the
!                    address of the element corresponding to
!                    (ILON=1,IGL) in the KPROMA arrays.
!          KVLAG   - Switch for discretization of continuity equation.
!          KWLAG   - Switch for discretization of momentum equation.
!          KITER   - Number of iterations for computing the medium point.
!          PVMAX1  - For a wind above this value, print of a warning.
!          PVMAX2  - For a wind above this value, CALL ABOR1.
!          LDPLANE - switch: .T. = plane geometry; .F. = spherical geometry.
!          LDSETTLS- .T./.F.: Stable/Conventional algorithm of trajectory
!                    research in the 2TLSL scheme.
!          LDELTRA - .T./.F.: "Elegant"/Conventional algorithm of trajectory
!                    research in the 2TLSL scheme.
!          PURL0   - U-component of the wind, to be interpolated.
!          PVRL0   - V-component of the wind, to be interpolated.
!          PURL    - U-component of the wind.
!          PVRL    - V-component of the wind.
!          KSTTYP  - 1: Not tilted pole;  2: Tilted pole.
!          PDSTRET - 2*c (where c is the stretching factor).
!          PC2M1   - c*c-1.
!          PC2P1   - c*c+1.
!          P4JP    - Approximative inverse of the differences of latitudes.
!          PI      - number PI
!          PDEPI   - 2 * PI
!          PIS2    - PI / 2
!          PDTSA   - (0.5 PDT /a)
!          PDTS62  - (0.5 PDT /a)*(0.5 PDT /a)/6
!          PDTS22  - (0.5 PDT /a)*(0.5 PDT /a)/2
!          PLOCEN  - geographical longitude of the stretching pole.
!          PMUCEN  - sine of the geographical latitude of the HR pole.
!          PLSDEPI - (Number of points by latitude) / (2 * PI) .
!          PLATI   - latitude on the computational sphere.
!          KIBL    - index into YRGSGEOM/YRCSGEOM instances in YDGEOMETRY

!        OUTPUT:
!          PLON    - computational sphere longitude of interpolation point.
!          PLAT    - computational sphere latitude of interpolation point.
!          PCOSCO  - cos(Longitude-Longitude(grid-point))*cos(Latitude)
!                    of the interpolation point (geographical longitude
!                    and latitude).
!          PSINCO  - sin(Longitude-Longitude(grid-point))*cos(Latitude)
!                    of the interpolation point (geographical longitude
!                    and latitude).
!          PSINLA  - sine of the interpolation point geographical latitude.
!          PCOPHI  - cosine of the geographical angle between the
!                    interpolation point and the grid-point.
!          PR      - first element of the wind displacement matrix.
!          PS      - second element of the wind displacement matrix.

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!        Calls  LARCIN2.
!        Is called by LADINE (2D shallow water model).

!     Reference.
!     ----------

!     Author.
!     -------
!      C. Temperton (ECMWF), by simplification of LARMES
!      Original : 98/01/19

!     Modifications.
!     --------------
!      Modified 01-08-30 by K. YESSAD: pruning and some other cleanings.
!      R. El Khatib : 01-08-07 Pruning options
!      01-Oct-2003 M. Hamrud   CY28 Cleaning
!      30-Jun-2008 J. Masek    Dataflow for new SLHD interpolators.
!      K. Yessad Nov 2008: rationalisation of dummy argument interfaces
!      G. Mozdzynski (Jan 2011): OOPS cleaning, use of derived type SL_STRUCT
!      G. Mozdzynski (Feb 2011): OOPS cleaning, use of derived types TGSGEOM and TCSGEOM
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

!-------------------------------------------------------------------------------

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
INTEGER(KIND=JPIM),INTENT(IN)    :: KULOUT 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTABUF(YDSL%NDGSAH:YDSL%NDGENH) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KVLAG 
INTEGER(KIND=JPIM),INTENT(IN)    :: KWLAG 
INTEGER(KIND=JPIM),INTENT(IN)    :: KITER 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVMAX1 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVMAX2 
LOGICAL           ,INTENT(IN)    :: LDPLANE 
LOGICAL           ,INTENT(IN)    :: LDSETTLS 
LOGICAL           ,INTENT(IN)    :: LDELTRA 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PURL0(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVRL0(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PURL(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVRL(KPROMA) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTTYP 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDSTRET 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PC2M1 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PC2P1 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P4JP 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PI 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDEPI 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PIS2 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDTSA 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDTS62 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDTS22 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLOCEN 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMUCEN 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLSDEPI(YDSL%NDGSAH:YDSL%NDGENH) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLATI(YDSL%NDGSAH:YDSL%NDGENH) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIBL
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PLON(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLAT(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCOSCO(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSINCO(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSINLA(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCOPHI(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PR(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PS(KPROMA) 

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZUF(KPROMA)
REAL(KIND=JPRB) :: ZVF(KPROMA)

INTEGER(KIND=JPIM) :: IWISA
INTEGER(KIND=JPIM) :: IXLAG(6)

INTEGER(KIND=JPIM) :: IROFEXP, IROT, JITER, JROF

LOGICAL :: LLQMHP, LLQMHW

REAL(KIND=JPRB) :: ZC0(1), ZCL0(1), ZDTS22, ZDTS62, ZDTSA, ZEW,&
 & ZEWX, ZINT, ZNOR, ZNOR2, ZNORX, ZPU, ZPV,&
 & ZSPHSV, ZU0(1), ZUL0(1), ZV0(1),&
 & ZVL0(1), ZVMAX1, ZVMAX2, ZZUZ0(1), ZZVZ0(1)  
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "larcin2.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LARMES2',0,ZHOOK_HANDLE)

!     ------------------------------------------------------------------
ASSOCIATE(YDCSGEOM=>YDGEOMETRY%YRCSGEOM(KIBL), YDGSGEOM=>YDGEOMETRY%YRGSGEOM(KIBL)) 
!*       1.    PRELIMINARY INITIALISATIONS AND TESTS.
!              --------------------------------------

ZNOR=0.0_JPRB
ZEW=0.0_JPRB
ZVMAX1=PVMAX1*PVMAX1
ZVMAX2=PVMAX2*PVMAX2

DO JROF=KSTART,KPROF
  ZNOR=MAX(PVRL(JROF)*PVRL(JROF),ZNOR)
  ZEW =MAX(PURL(JROF)*PURL(JROF),ZEW )
ENDDO

IF (ZNOR > ZVMAX1) THEN
  WRITE(KULOUT,*) ' MAX V WIND=',SQRT(ZNOR)
ENDIF
IF (ZEW > ZVMAX1) THEN
  WRITE(KULOUT,*) ' MAX U WIND=',SQRT(ZEW)
ENDIF
IF (ZNOR > ZVMAX2) THEN
  ZNOR=0.0_JPRB
  DO JROF=KSTART,KPROF
    ZNORX=PVRL(JROF)*PVRL(JROF)
    ZNOR=MAX(ZNORX,ZNOR)
    IF(ZNOR == ZNORX) THEN
      IROFEXP=JROF
    ENDIF
  ENDDO
  WRITE(KULOUT,*) ' V WIND =',SQRT(ZNOR),' IS TOO STRONG, EXPLOSION.'
  WRITE(KULOUT,*) ' POINT= ',IROFEXP
  WRITE(KULOUT,*) ' YDCSGEOM%RCOLON= ',YDCSGEOM%RCOLON(IROFEXP)
  WRITE(KULOUT,*) ' YDGSGEOM%GEMU = ',YDGSGEOM%GEMU(IROFEXP)
  CALL ABOR1(' !V WIND TOO STRONG, EXPLOSION!!!')
ENDIF
IF (ZEW > ZVMAX2) THEN
  ZEW=0.0_JPRB
  DO JROF=KSTART,KPROF
    ZEWX=PURL(JROF)*PURL(JROF)
    ZEW=MAX(ZEWX,ZEW)
    IF(ZEW == ZEWX) THEN
      IROFEXP=JROF
    ENDIF
  ENDDO
  WRITE(KULOUT,*) ' U WIND =',SQRT(ZEW),' IS TOO STRONG, EXPLOSION.'
  WRITE(KULOUT,*) ' POINT= ',IROFEXP
  WRITE(KULOUT,*) ' YDCSGEOM%RCOLON= ',YDCSGEOM%RCOLON(IROFEXP)
  WRITE(KULOUT,*) ' YDGSGEOM%GEMU = ',YDGSGEOM%GEMU(IROFEXP)
  CALL ABOR1(' !U WIND TOO STRONG, EXPLOSION!!!')
ENDIF

IXLAG(1)=KVLAG
IXLAG(4)=KWLAG

IF((LDSETTLS.OR.LDELTRA).AND.(KITER > 1)) THEN
  ZDTS22=PDTS22*4
  ZDTSA=PDTSA*2
  ZDTS62=PDTS62*4
ELSE
  ZDTS22=PDTS22
  ZDTSA=PDTSA
  ZDTS62=PDTS62
ENDIF

!     ------------------------------------------------------------------

!*       2.    ITERATIONS.
!              -----------

DO JITER=1,KITER

!*       2.1   DETERMINATION OF THE MEDIUM POINT "M".

!       Computation of the norm of the real wind vector.
!       Computation of the angle (PHI=DT . NORM OF V ON RADIUS)**2
!       then computation of the coordinates of the medium point "M".
!       If LDSETTLS=T and JITER < KITER the origin point "O" is computed
!       instead of the medium point "M", and "M" is computed only at
!       the last iteration (JITER = KITER).

  IF (JITER == 1) THEN

    DO JROF=KSTART,KPROF

!             computations on horizontal plans.

      ZNOR2= ( PURL(JROF)*PURL(JROF)+PVRL(JROF)*PVRL(JROF) )
      PCOPHI(JROF)=1.0_JPRB-ZDTS22*ZNOR2
      ZSPHSV           =ZDTSA*(1.0_JPRB-ZDTS62*ZNOR2)
      ZINT= ( PURL(JROF)*YDGSGEOM%GNORDL(JROF)+PVRL(JROF)*YDGSGEOM%GNORDM(JROF) )*ZSPHSV
      PSINLA(JROF)= YDGSGEOM%GEMU(JROF)*PCOPHI(JROF)-ZINT*YDGSGEOM%GSQM2(JROF)
      PCOSCO(JROF)= YDGSGEOM%GSQM2(JROF)*PCOPHI(JROF)+ZINT* YDGSGEOM%GEMU(JROF)
      PSINCO(JROF)=-( PURL(JROF)*YDGSGEOM%GNORDM(JROF)&
       & -PVRL(JROF)*YDGSGEOM%GNORDL(JROF) )*ZSPHSV  

    ENDDO

  ELSE

    DO JROF=KSTART,KPROF

!           ZPU,ZPV are the coordinates of VM in the local repere
!           related to G.

      IF(LDSETTLS) THEN
        ZPU= 0.5_JPRB*(PR(JROF)*ZUF(JROF)&
         & +PS(JROF)*ZVF(JROF)&
         & +PURL(JROF)*YDGSGEOM%GNORDM(JROF)&
         & -PVRL(JROF)*YDGSGEOM%GNORDL(JROF))  
        ZPV=0.5_JPRB*(-PS(JROF)*ZUF(JROF)&
         & +PR(JROF)*ZVF(JROF)&
         & +PURL(JROF)*YDGSGEOM%GNORDL(JROF)&
         & +PVRL(JROF)*YDGSGEOM%GNORDM(JROF))  
      ELSE
        ZPU= PR(JROF)*ZUF(JROF)+PS(JROF)*ZVF(JROF)
        ZPV=-PS(JROF)*ZUF(JROF)+PR(JROF)*ZVF(JROF)
      ENDIF
      ZNOR2             =(ZPU*ZPU+ZPV*ZPV)
      IF (JITER == KITER) THEN
        PCOPHI(JROF) =1.0_JPRB-PDTS22*ZNOR2
        ZSPHSV       =PDTSA*(1.0_JPRB-PDTS62*ZNOR2)
      ELSE
        PCOPHI(JROF) =1.0_JPRB-ZDTS22*ZNOR2
        ZSPHSV       =ZDTSA*(1.0_JPRB-ZDTS62*ZNOR2)
      ENDIF
      ZINT              =ZPV*ZSPHSV
      PSINLA(JROF) = YDGSGEOM%GEMU(JROF)*PCOPHI(JROF)-ZINT*YDGSGEOM%GSQM2(JROF)
      PCOSCO(JROF) = YDGSGEOM%GSQM2(JROF)*PCOPHI(JROF)+ZINT* YDGSGEOM%GEMU(JROF)
      PSINCO(JROF) =-ZPU*ZSPHSV

    ENDDO

  ENDIF

!*       2.2   DETERMINATION OF THE WIND AT THE MEDIUM POINT "I".

  IF (JITER /= KITER) THEN

!         If LDSETTLS=T and JITER < KITER the wind is interpolated at
!         the origin point "O" instead of at the medium point "M".
!         In the other cases the wind is interpolated at "M".

    IWISA=1
    IROT=1

!          Warning: ZUL0,ZVL0,ZCL0,
!                   ZU0,ZV0,ZC0,
!                   ZZUZ0,ZZVZ0,LLQMHW,LLQMHP
!                   are not used in this call of LARCIN2.

    CALL LARCIN2(YDGEOMETRY,YDDYN,YDDYNA,YDSL,YDSLINT,KPROMA,KSTART,KPROF,KFLDN,KFLDX,&
     & KSTABUF,&
     & IWISA,IXLAG,IROT,LDPLANE,&
     & LLQMHW,LLQMHP,&
     & KSTTYP,PDSTRET,PC2M1,PC2P1,P4JP,PI,PDEPI,PIS2,&
     & PLOCEN,PMUCEN,PLSDEPI,PLATI,&
     & KIBL,&
     & PCOSCO,PSINCO,PSINLA,PCOPHI,&
     & PURL0,PVRL0,&
     & ZUL0,ZVL0,ZCL0,ZUL0,ZVL0,ZCL0,&
     & PLON,PLAT,&
     & PR,PS,&
     & ZUF,ZVF,&
     & ZU0,ZV0,ZC0,&
     & ZZUZ0,ZZVZ0)

  ENDIF

ENDDO

!     ------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('LARMES2',1,ZHOOK_HANDLE)
END SUBROUTINE LARMES2
