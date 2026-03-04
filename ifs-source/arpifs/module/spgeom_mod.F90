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

MODULE SPGEOM_MOD

! Module for spectral geometry structures.

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YEMGEO   , ONLY : TEGEO
USE YOMDIM   , ONLY : TDIM
USE YOMLUN   , ONLY : NULOUT
USE YOMCST   , ONLY : RPI, RA
USE YOMMP0   , ONLY : NPRINTLEV
USE YOMCT0   , ONLY : LELAM, LALLOPR

!     ------------------------------------------------------------------

IMPLICIT NONE
SAVE

!!! GMR     : coefficients for spectral multiplication by GM.
!!! SCGMAP  : coefficients for multiplication by (GM**2) in spectral space (global model).
!!! ESCGMAP : coefficients for multiplication by (GM**2) in spectral space (LAM model).

CONTAINS

SUBROUTINE SUSPGEOM(YDGEOMETRY)

!--------------------------------------------------------------------------
! Sets-up and allocates spectral geometry structures
!--------------------------------------------------------------------------

USE TYPE_GEOMETRY , ONLY : GEOMETRY
IMPLICIT NONE
TYPE(GEOMETRY), INTENT(INOUT) :: YDGEOMETRY

LOGICAL :: LLP
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SPGEOM_MOD:SUSPGEOM',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM, &
 & YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, &
 & YDMP=>YDGEOMETRY%YRMP, YDLAP=>YDGEOMETRY%YRLAP, &
 &  YDSPGEOM=>YDGEOMETRY%YSPGEOM)
ASSOCIATE(NSPEC=>YDDIM%NSPEC)
! ----------------------------------------------------------------------

!*       1.    ALLOCATION OF ARRAYS.
!              ---------------------


LLP = NPRINTLEV >= 1.OR. LALLOPR
ALLOCATE(YDSPGEOM%GMR(NSPEC,2))
IF (LLP) WRITE(NULOUT,"(1X,'ARRAY ',A10,' ALLOCATED ',8I8)") 'GMR      ',SIZE(YDSPGEOM%GMR),SHAPE(YDSPGEOM%GMR)
ALLOCATE(YDSPGEOM%SCGMAP(NSPEC,3))
IF (LLP) WRITE(NULOUT,"(1X,'ARRAY ',A10,' DECLARED  ',8I8)") 'SCGMAP   ',SIZE(YDSPGEOM%SCGMAP),SHAPE(YDSPGEOM%SCGMAP)

!*       2.    INITIALIZATION OF ARRAYS.
!              -------------------------

IF (.NOT. LELAM) THEN
  CALL SUSPGM(YDGEOMETRY,NULOUT,YDSPGEOM%GMR(1,1),YDSPGEOM%GMR(1,2))
  CALL SUSMAP(YDGEOMETRY,NULOUT,YDSPGEOM%SCGMAP(1,1),YDSPGEOM%SCGMAP(1,2),YDSPGEOM%SCGMAP(1,3))
ELSE
  CALL SUESMAP(YDGEOMETRY%YREGEO,NULOUT,YDSPGEOM%ESCGMAP(1),YDSPGEOM%ESCGMAP(2),YDSPGEOM%ESCGMAP(3))
ENDIF

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SPGEOM_MOD:SUSPGEOM',1,ZHOOK_HANDLE)
END SUBROUTINE SUSPGEOM

! ----------------------------------------------------------------------
! Inner subroutines
! ----------------------------------------------------------------------

SUBROUTINE SUSMAP(YDGEOMETRY,KULOUT,PD,PE,PF)

!**** *SUSMAP*   - Initialisation of coefficients allowing spectral
!                  multiplication by GM**2 (GM is the map factor).
!                  Caution! Only valid with a triangular truncation.

!    ((GM**2)*SPX)(m,n)=PD(m,n  )*SPX(m,n  )
!     +PE(m,n-1)*SPX(m,n-1)+PE(m,n  )*SPX(m,n+1)
!     +PF(m,n-2)*SPX(m,n-2)+PF(m,n  )*SPX(m,n+2)
!    where:
!         a=0.5*(c+1/c); b=0.5*(c-1/c); c=RSTRET: stretching coefficient.
!         e(m,n)=sqrt((n**2-m**2)/(4n**2-1))
!         PD(m,n)=(a**2)+(b**2)[(e(m,n)**2)+(e(m,n+1)**2)]
!         PE(m,n)=(2ab)[e(m,n+1)]
!         PF(m,n)=(b**2)[e(m,n+1)][e(m,n+2)]

!    Explicit arguments :
!    --------------------
!     INPUT:
!      KULOUT:    - Output logical unit.

!     OUTPUT:
!      PD,PE,PF:  - Coefficients to be computed (see above).

!    Author.
!    -------
!      K. YESSAD: OCTOBER 1993.
!     ------------------------------------------------------------------

USE TYPE_GEOMETRY , ONLY : GEOMETRY
IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)  :: YDGEOMETRY
INTEGER(KIND=JPIM),INTENT(IN)  :: KULOUT 
REAL(KIND=JPRB)   ,INTENT(OUT) :: PD(YDGEOMETRY%YRDIM%NSPEC) 
REAL(KIND=JPRB)   ,INTENT(OUT) :: PE(YDGEOMETRY%YRDIM%NSPEC) 
REAL(KIND=JPRB)   ,INTENT(OUT) :: PF(YDGEOMETRY%YRDIM%NSPEC) 

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: Z2AB, ZA, ZA2, ZB, ZB2
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SPGEOM_MOD:SUSMAP',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM, &
  & YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, YDLAP=>YDGEOMETRY%YRLAP)
ASSOCIATE(RSTRET=>YDGEM%RSTRET)
!     ------------------------------------------------------------------

!*       1.    COMPUTATION OF PD,PE,PF.
!              ------------------------

WRITE (KULOUT,*)
WRITE (KULOUT,'(A,A)') ' SUSMAP: CALL TO SUGMRE TO ',&
 & 'INITIALIZE COEFFICIENTS FOR SPECTRAL MULT. BY GM**2.'  

ZA=0.5_JPRB*(RSTRET+1.0_JPRB/RSTRET)
ZB=0.5_JPRB*(RSTRET-1.0_JPRB/RSTRET)
ZA2=ZA*ZA
ZB2=ZB*ZB
Z2AB=2.0_JPRB*ZA*ZB
CALL SUGMRE(YDLAP,YDGEOMETRY%YRDIM,KULOUT,ZA2,Z2AB,ZB2,PD,PE,PF)

!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SPGEOM_MOD:SUSMAP',1,ZHOOK_HANDLE)
END SUBROUTINE SUSMAP

SUBROUTINE SUGMRE(YDLAP,YDDIM,KULOUT,PGMRA0,PGMRA1,PGMRA2,PD,PE,PF)

!**** *SUGMRE*   - Initialisation of coefficients allowing spectral
!                  multiplication by a second-order polynomial of "mu"
!                  (for example (map factor)**2).

!    ((PGMR)*SPX)(m,n)=PD(m,n  )*SPX(m,n  )
!     +PE(m,n-1)*SPX(m,n-1)+PE(m,n  )*SPX(m,n+1)
!     +PF(m,n-2)*SPX(m,n-2)+PF(m,n  )*SPX(m,n+2)

!    where:
!         pgmr = a2 mu**2 + a1 mu + a0
!         (mu = sinus of latitude)

!    Explicit arguments :
!    --------------------
!     INPUT:
!      KULOUT:    - Output logical unit.
!      PGMRA0:    - Coefficient a0 of second order polynomial.
!      PGMRA1:    - Coefficient a1 of second order polynomial.
!      PGMRA2:    - Coefficient a2 of second order polynomial.

!     OUTPUT:
!      PD,PE,PF:  - Coefficients to be computed (see above).

!    Author.
!    -------
!      K. YESSAD: APRIL 1994.
!     ------------------------------------------------------------------

USE YOMLAP , ONLY : TLAP
USE YOMDIM , ONLY : TDIM
IMPLICIT NONE

TYPE(TLAP)        , INTENT(IN)   :: YDLAP
TYPE(TDIM)        , INTENT(IN)   :: YDDIM
INTEGER(KIND=JPIM), INTENT(IN)   :: KULOUT 
REAL(KIND=JPRB)   , INTENT(IN)   :: PGMRA0 
REAL(KIND=JPRB)   , INTENT(IN)   :: PGMRA1 
REAL(KIND=JPRB)   , INTENT(IN)   :: PGMRA2 
REAL(KIND=JPRB)   , INTENT(OUT)  :: PD(YDDIM%NSPEC) 
REAL(KIND=JPRB)   , INTENT(OUT)  :: PE(YDDIM%NSPEC) 
REAL(KIND=JPRB)   , INTENT(OUT)  :: PF(YDDIM%NSPEC) 

!      ----------------------------------------------------------------

REAL(KIND=JPRB) :: ZEE(0:YDDIM%NSMAX+2)
INTEGER(KIND=JPIM) :: ISE, ISE0, JMLOC, IM, JN
REAL(KIND=JPRB) :: ZM2, ZN2
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!      ----------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SPGEOM_MOD:SUGMRE',0,ZHOOK_HANDLE)
ASSOCIATE(NUMP=>YDDIM%NUMP, NSMAX=>YDDIM%NSMAX, &
 & MYMS=>YDLAP%MYMS, NSE0L=>YDLAP%NSE0L)
!      ----------------------------------------------------------------

!*       1.    COMPUTATION OF PD,PE,PF.
!              ------------------------

DO JMLOC=1,NUMP
  IM=MYMS(JMLOC)
  ZM2=REAL(IM*IM,JPRB)
  DO JN=IM,NSMAX+2
    ZN2=REAL(JN*JN,JPRB)
    ZEE(JN)=SQRT(MAX(0.0_JPRB,(ZN2-ZM2)/(4._JPRB*ZN2-1.0_JPRB)))
  ENDDO
  ISE0=NSE0L(JMLOC)
  DO JN=IM,NSMAX
    ISE=ISE0+JN+1-IM
    PD(ISE)=PGMRA0+PGMRA2*(ZEE(JN)**2+ZEE(JN+1)**2)
  ENDDO
  DO JN=IM,NSMAX-1
    ISE=ISE0+JN+1-IM
    PE(ISE)=PGMRA1*ZEE(JN+1)
  ENDDO
  ISE=ISE0+NSMAX+1-IM
  PE(ISE)=0.0_JPRB
  DO JN=IM,NSMAX-2
    ISE=ISE0+JN+1-IM
    PF(ISE)=PGMRA2*ZEE(JN+1)*ZEE(JN+2)
  ENDDO
  DO JN=NSMAX-1,NSMAX
    ISE=ISE0+JN+1-IM
    PF(ISE)=0.0_JPRB
  ENDDO
ENDDO

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SPGEOM_MOD:SUGMRE',1,ZHOOK_HANDLE)
END SUBROUTINE SUGMRE

SUBROUTINE SUSPGM(YDGEOMETRY,KULOUT,PD,PE)

!**** *SUSPGM*   - Initialisation of coefficients allowing spectral
!                  multiplication by the mapping factor M.

!     (M*SPX)(m,n)=PD(m,n  )*SPX(m,n  )
!      +PE(m,n-1)*SPX(m,n-1)+PE(m,n  )*SPX(m,n+1)

!     M = 0.5(c+1/c) + 0.5 (c-1/c) mu
!     (mu = sinus of latitude)

!    Explicit arguments :
!    --------------------
!     INPUT:
!      KULOUT:    - Output logical unit.

!     OUTPUT:
!      PD,PE:     - Coefficients to be computed (see above).

!    Author.
!    -------
!        K. YESSAD: DEC 2003.
!     ------------------------------------------------------------------

USE TYPE_GEOMETRY , ONLY : GEOMETRY
IMPLICIT NONE

TYPE(GEOMETRY)    , INTENT(IN)  :: YDGEOMETRY
INTEGER(KIND=JPIM), INTENT(IN)  :: KULOUT 
REAL(KIND=JPRB)   , INTENT(OUT) :: PD(YDGEOMETRY%YRDIM%NSPEC) 
REAL(KIND=JPRB)   , INTENT(OUT) :: PE(YDGEOMETRY%YRDIM%NSPEC) 

!      ----------------------------------------------------------------

REAL(KIND=JPRB) :: ZEE(0:YDGEOMETRY%YRDIM%NSMAX+2)
INTEGER(KIND=JPIM) :: ISE, ISE0, JMLOC, IM, JN
REAL(KIND=JPRB) :: ZM2, ZN2
REAL(KIND=JPRB) :: ZGMRA0 
REAL(KIND=JPRB) :: ZGMRA1 
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!      ----------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SPGEOM_MOD:SUSPGM',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM, &
  & YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, YDLAP=>YDGEOMETRY%YRLAP)
ASSOCIATE(NUMP=>YDDIM%NUMP, NSMAX=>YDDIM%NSMAX, &
 & RSTRET=>YDGEM%RSTRET, &
 & MYMS=>YDLAP%MYMS, NSE0L=>YDLAP%NSE0L)
!      ----------------------------------------------------------------

!*       1.    COMPUTATION OF PD,PE.
!              ---------------------

WRITE (KULOUT,*)
WRITE (KULOUT,'(A)') ' ROUTINE SUSPGM:'

ZGMRA0=0.5_JPRB*(RSTRET + 1.0_JPRB/RSTRET)
ZGMRA1=0.5_JPRB*(RSTRET - 1.0_JPRB/RSTRET)

DO JMLOC=1,NUMP
  IM=MYMS(JMLOC)
  ZM2=REAL(IM*IM,JPRB)
  DO JN=IM,NSMAX+2
    ZN2=REAL(JN*JN,JPRB)
    ZEE(JN)=SQRT(MAX(0.0_JPRB,(ZN2-ZM2)/(4._JPRB*ZN2-1.0_JPRB)))
  ENDDO
  ISE0=NSE0L(JMLOC)
  DO JN=IM,NSMAX
    ISE=ISE0+JN+1-IM
    PD(ISE)=ZGMRA0
  ENDDO
  DO JN=IM,NSMAX-1
    ISE=ISE0+JN+1-IM
    PE(ISE)=ZGMRA1*ZEE(JN+1)
  ENDDO
  ISE=ISE0+NSMAX+1-IM
  PE(ISE)=0.0_JPRB
ENDDO

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SPGEOM_MOD:SUSPGM',1,ZHOOK_HANDLE)
END SUBROUTINE SUSPGM

SUBROUTINE SUESMAP(YDEGEO,KULOUT,PD,PE,PF)

!**** *SUESMAP*   - Initialisation of coefficients PD,PE,PF allowing spectral
!                   multiplication by GM**2 (GM is the map factor), used for
!                   the RTM variable map factor.

!     ((GM**2)*SPX)(m,n)=PD*SPX(m,n  )
!      +PE*SPX(m,n-1)+PE*SPX(m,n+1)
!      +PF*SPX(m,n-2)+PF*SPX(m,n+2)

!     where:
!      PD=(1/2f)(exp(f)-exp(-f))+1
!      PE=(-f/2)*(exp(f)-exp(-f))/(pi**2+f**2)
!      PF=(f/2)*(exp(f)-exp(-f))/((2*pi)**2+f**2)
!     are the Fourier coefficients for the RTM map factor and:
!      f=factor=Ly/a
!      Ly=size of the area
!      a=radius of the earth
 
!    Explicit arguments :
!    --------------------
!     INPUT:
!      KULOUT:    - Output logical unit.
 
!     OUTPUT:
!      PD,PE,PF:  - Coefficients to be computed (see above).
 
!    Author.
!    -------
!        I. Santos and I. Martinez (MARCH 2010)
!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TEGEO)       , INTENT(IN)   :: YDEGEO
INTEGER(KIND=JPIM), INTENT(IN)   :: KULOUT
REAL(KIND=JPRB)   , INTENT(OUT)  :: PD
REAL(KIND=JPRB)   , INTENT(OUT)  :: PE
REAL(KIND=JPRB)   , INTENT(OUT)  :: PF

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZFACTOR
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SPGEOM_MOD:SUESMAP',0,ZHOOK_HANDLE)
!      ----------------------------------------------------------------

!*             COMPUTATION OF PD,PE,PF.
!              ------------------------

WRITE (KULOUT,*) ' SUESMAP: INITIALIZE COEFFICIENTS FOR SPECTRAL MULT. BY GM**2.'

ZFACTOR = YDEGEO%ELY/RA
PD=(EXP(ZFACTOR)-EXP(-ZFACTOR))/(2.0_JPRB*ZFACTOR)+1.0_JPRB
PE=-ZFACTOR/2.0_JPRB*(EXP(ZFACTOR)-EXP(-ZFACTOR))/(RPI**2+ZFACTOR**2)
PF=ZFACTOR/2.0_JPRB*(EXP(ZFACTOR)-EXP(-ZFACTOR))/((RPI*2.0_JPRB)**2+ZFACTOR**2)

PD = PD/2.0_JPRB
PE = PE/2.0_JPRB
PF = PF/2.0_JPRB

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SPGEOM_MOD:SUESMAP',1,ZHOOK_HANDLE)
END SUBROUTINE SUESMAP

END MODULE SPGEOM_MOD
