! (C) Copyright 1987- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE SUGAW_MOD
CONTAINS
SUBROUTINE SUGAW(KDGL,KM,KN,PL,PW,PANM,PFN)

USE PARKIND1  ,ONLY : JPRD, JPIM
USE PARKIND2  ,ONLY : JPRH

USE TPM_CONSTANTS   ,ONLY : RA

USE TPM_GEN         ,ONLY : NOUT
USE GAWL_MOD        ,ONLY : GAWL
USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS
USE SUPOLF_MOD
USE TPM_POL

!**** *SUGAW * - Routine to initialize the Gaussian
!                 abcissa and the associated weights

!     Purpose.
!     --------
!           Initialize arrays PL, and PW (quadrature abscissas and weights)
!**   Interface.
!     ----------
!        *CALL* *SUGAW(KN,PFN,PL,PW) *

!        Explicit arguments :
!        --------------------
!           INPUT:
!              KDGL     :  Number of Gauss  abscissas 
!              KM       :  Polynomial order m
!              KN       :  Polynomial degree n
!             PFN       :  Fourier coefficients of series expansion for
!                          the ordinary Legendre polynomials
!           OUTPUT:
!              PL (KN)  :  abscissas of Gauss
!              PW (KN)  :  Weights of the Gaussian integration

!     PL (i) is the abscissa i starting from the northern pole, it is
! the cosine of the colatitude of the corresponding row of the collocation
! grid.

!        Implicit arguments :
!        --------------------
!       None

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------

!     S.L. Belousov, Tables of normalized associated Legendre Polynomials, Pergamon Press (1962)
!     P.N. Swarztrauber, On computing the points and weights for Gauss-Legendre quadrature,
!     SIAM J. Sci. Comput. Vol. 24 (3) pp. 945-954 (2002)

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original          : 87-10-15
!        Michel Rochas     : 90-08-30
!        Philippe Courtier : 92-12-19 Multitasking
!        Ryad El Khatib    : 94-04-20 Remove unused comdecks pardim and yomdim
!        Mats Hamrud       : 94-08-12 Printing level
!        K. Yessad (Sep 2008): cleaning, improve comments.
!        Nils Wedi + Mats Hamrud, 2009-02-05 revised following Swarztrauber, 2002
!      F. Vana  05-Mar-2015  Support for single precision
!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KDGL
INTEGER(KIND=JPIM),INTENT(IN) :: KM
INTEGER(KIND=JPIM),INTENT(IN) :: KN

REAL(KIND=JPRD)   ,INTENT(IN)  :: PANM

REAL(KIND=JPRD),INTENT(OUT) :: PW(KDGL)
REAL(KIND=JPRD),INTENT(OUT) :: PL(KDGL)

REAL(KIND=JPRD)   ,OPTIONAL, INTENT(IN)  :: PFN(0:KDGL,0:KDGL)

!     ------------------------------------------------------------------

REAL(KIND=JPRD) :: ZLI(KDGL),ZT(KDGL),ZFN(0:KDGL/2),ZL(KDGL)
REAL(KIND=JPRD) :: ZREG(KDGL),ZMOD(KDGL),ZM(KDGL),ZRR(KDGL)
INTEGER(KIND=JPIM) :: ITER(KDGL)

INTEGER(KIND=JPIM) :: IALLOW, INS2, ISYM, JGL, IK, IODD, I, IMAX

REAL(KIND=JPRD) :: Z, ZEPS, Z0, ZPI

! computations in extended precision for alternative root finding
! which also works for associated polynomials (m>0)
REAL(KIND=JPRH) :: ZLK, ZLK1, ZLLDN, ZANM
REAL(KIND=JPRH) :: ZTHETA, ZTHETA0, ZX, ZX0, ZDX0, ZH, ZPIH, ZS0
REAL(KIND=JPRH) :: ZK1, ZK2, ZK3, ZK4
REAL(KIND=JPRH) :: ZF1, ZF2, ZF3
REAL(KIND=JPRH) :: FP, FQ, FP1, FQ1
REAL(KIND=JPRH) :: X, ZXOLD, ZBIG, ZEPSH

INTEGER(KIND=JPIM) :: ISTEPMAX

LOGICAL :: LLP2, LLREF, LLOLD

REAL(KIND=JPRD) :: ZDDPOL(0:KN)

INTEGER(KIND=JPIM), PARAMETER :: JPKD=KIND(ZLK)

FP(X) = 1._JPRH-X**2
FQ(X) = REAL(KN*(KN+1),JPRH)-REAL(KM**2,JPRH)/(1._JPRH-X**2)
FP1(X) = -2._JPRH*X
FQ1(X) = -2._JPRH*X*REAL(KM**2,JPRH)/SQRT(1._JPRH-X**2)

!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!*       1. Initialization + root + weight computation
!           ------------------------------------------

LLP2 = .FALSE.
INS2 = KDGL/2

LLOLD=( KM == 0 .AND. KN == KDGL ).AND.PRESENT(PFN)


CALL GSTATS(1650,0)

ZEPS  = EPSILON(Z)
ZEPSH = EPSILON(X)

ZBIG  = SQRT(HUGE(X))

!*       1.1 Find the roots of the ordinary
!           Legendre polynomial of degree KN using an analytical first guess
!           and then refine to machine precision via Newton's method
!           in double precision following Swarztrauber (2002)

!        Nils Comment: in principle the else case could also be used for this but
!                      this is slightly more accurate and consistent with the past

IF( LLOLD ) THEN

  ZPI  = 2.0_JPRD*ASIN(1.0_JPRD)
  IODD=MOD(KDGL,2)
  IK=IODD
  DO JGL=IODD,KDGL,2
    ZFN(IK)=PFN(KDGL,JGL)
    IK=IK+1
  ENDDO

  DO JGL=1,INS2
    Z = REAL(4*JGL-1,JPRD)*ZPI/REAL(4*KN+2,JPRD)
    ! analytic initial guess for cos(theta) (same quality as RK below)
    ! ZX = 1._JPRD-REAL(KN-1,JPRD)/REAL(8*KN*KN*KN,JPRD)-(1._JPRD/REAL(384*KN*KN*KN*KN))*(39._JPRD-28._JPRD/SIN(Z)**2)
    ! PL(JGL) = ACOS(ZX*COS(Z))
    ZL(JGL) = Z+1.0_JPRD/(TAN(Z)*REAL(8*KN**2,JPRD))
    ZREG(JGL) = COS(Z)
    ZLI(JGL) = COS(ZL(JGL))
  ENDDO

  ! refine PL here via Newton's method

  !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JGL)
  DO JGL=INS2,1,-1
    CALL GAWL(ZFN,ZL(JGL),PW(JGL),ZEPS,KN,ITER(JGL),ZMOD(JGL))
  ENDDO
  !$OMP END PARALLEL DO

  ! convert to physical latitude space PMU
  !DIR$ IVDEP
  !OCL NOVREC
  DO JGL=1,INS2
    PL(JGL) = COS(ZL(JGL))
  ENDDO

ELSE

!*       1.2 Find the roots of the associated
!            Legendre polynomial of degree KN and the associated Gaussian weights 
!            using a Runge-Kutta 4 integration of the Pruefer transformed Sturm-Liouville problem 
!            (Tygert (J. Comput. Phys. 2008) and Glaser et al., SIAM J. SCI. COMPUT. Vol. 29 (4) 1420-1438)
! 

  ISTEPMAX=10

  ZANM = REAL(PANM, JPKD)
  ZPIH = 2.0_JPRH*ASIN(1.0_JPRH)

  ZX0 = 0._JPRH
  Z0  = 0._JPRD

  ! first guess starting point
  IF( MOD(KN-KM,2) == 0 ) THEN
    ! even, extremum at X == 0 
    ZTHETA0 = 0._JPRH
    ZH = -0.5_JPRH*ZPIH/REAL(ISTEPMAX,JPRH)
  ELSE
    ! odd, root at X == 0
    ZTHETA0 = 0.5_JPRH*ZPIH
    ZX0 = 0._JPRH
    ZH = -ZPIH/REAL(ISTEPMAX,JPRH)
  ENDIF
  
  ZX = ZX0
  ZTHETA = ZTHETA0

  ZF1 = SQRT(FQ(ZX)/FP(ZX))
  ZF2 = FQ1(ZX)/FQ(ZX)
  ZF3 = FP1(ZX)/FP(ZX)

  ! Formula (81) in Tygert
  ZDX0=-1._JPRH/(ZF1 + 0.25_JPRH*(ZF2 + ZF3)*SIN(2._JPRH*ZTHETA))

  ! loop over all roots
  LLREF=.TRUE.
  DO JGL=INS2,1,-1
    
    ! runge-kutta
    DGL:DO IK=1,ISTEPMAX

      ZK1 = ZDX0
      ZTHETA  = ZTHETA + 0.5_JPRH*ZH

      ZX = ZX0 + 0.5_JPRH*ZH*ZK1
      
      ZF1 = SQRT(FQ(ZX)/FP(ZX))
      ZF2 = FQ1(ZX)/FQ(ZX)
      ZF3 = FP1(ZX)/FP(ZX)
      
      ZK2 = -1._JPRH/(ZF1 + 0.25_JPRH*(ZF2 + ZF3)*SIN(2._JPRH*ZTHETA))
      ZX = ZX0 + 0.5_JPRH*ZH*ZK2

      ZF1 = SQRT(FQ(ZX)/FP(ZX))
      ZF2 = FQ1(ZX)/FQ(ZX)
      ZF3 = FP1(ZX)/FP(ZX)

      ZK3 = -1._JPRH/(ZF1 + 0.25_JPRH*(ZF2 + ZF3)*SIN(2._JPRH*ZTHETA))
      ZTHETA  = ZTHETA + 0.5_JPRH*ZH
      ZX = ZX0 + ZH*ZK3

      ZF1 = SQRT(FQ(ZX)/FP(ZX))
      ZF2 = FQ1(ZX)/FQ(ZX)
      ZF3 = FP1(ZX)/FP(ZX)

      ZK4 = -1._JPRH/(ZF1 + 0.25_JPRH*(ZF2 + ZF3)*SIN(2._JPRH*ZTHETA))
      ZX = ZX0 + (1._JPRH/6._JPRH)*ZH*(ZK1+2._JPRH*ZK2+2._JPRH*ZK3+ZK4)
      ZXOLD = ZX0

      ZX0 = ZX

      IF( .NOT.ZX==ZX ) THEN
        WRITE(NOUT,*) 'invoke overflow ...ZX ',KM, KN, JGL
        ZX  = ZXOLD
        ZX0 = ZXOLD
        EXIT DGL
      ENDIF
      
      ZF1 = SQRT(FQ(ZX)/FP(ZX))
      ZF2 = FQ1(ZX)/FQ(ZX)
      ZF3 = FP1(ZX)/FP(ZX)  

      ZDX0 = -1._JPRH/(ZF1 + 0.25_JPRH*(ZF2 + ZF3)*SIN(2._JPRH*ZTHETA))

    ENDDO DGL

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Everything from here until <<END>> is to refine the 
! root and compute the starting point for the next root search
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    ! should not happen, but does if loss of accuracy in supolf occurs (useful for debugging)
    IF( JGL < INS2 ) LLREF = PW(JGL+1).GT.ZEPSH
      
    IF( LLREF ) THEN

      ! chosen for speed/accuracy compromise
      IMAX=3
      LOOP: DO I=1,IMAX
        ! supol fast
        ZS0 = ACOS(ZX0)
        CALL SUPOLF(KM,KN,REAL(ZX0,JPRD),ZDDPOL)
        ZLK=REAL(ZDDPOL(KN),JPKD)
        ZLK1= REAL(ZDDPOL(KN-1),JPKD)
        ZLLDN= -(ZANM*ZLK1-DDI(KN)*COS(ZS0)*ZLK)/SIN(ZS0)
        
        IF( ABS(ZLLDN) > ZEPSH ) THEN
          ! single Newton refinement in theta
          ZS0 = ZS0 - ZLK/ZLLDN
          ZX  = COS(ZS0)
        ELSE
          ! do nothing
          ZX = ZX0
        ENDIF
        
        IF( ABS(ZX-ZX0) > 1000._JPRD*ZEPS ) THEN
          ZX0 = ZX
        ELSE
          EXIT LOOP
        ENDIF
      ENDDO LOOP
      
      ! recompute for accuracy weights
      CALL SUPOLF(KM,KN,REAL(ZX,JPRD),ZDDPOL)
      ! option f in Schwarztrauber to compute the weights
      ZS0 = ACOS(ZX)
      ZLK=REAL(ZDDPOL(KN),JPKD)
      ZLK1= REAL(ZDDPOL(KN-1),JPKD)
      ZLLDN= -(ZANM*ZLK1-DDI(KN)*COS(ZS0)*ZLK)/SIN(ZS0)
      
      PW(JGL) = REAL(REAL(2*KN+1,JPRH)/ZLLDN**2,JPRD)
      
      ! catch overflow, should never happen
      IF( .NOT.(PW(JGL)==PW(JGL)) ) THEN
        WRITE(NOUT,*) 'invoke overflow ...PW ',KM, KN, JGL
        PW(JGL) = 0.0_JPRD
      ENDIF
      
    ELSE
      ! should never happen ...
      WRITE(NOUT,*) 'Refinement not possible ... PW set to 0',KM, KN, JGL
      PW(JGL) = 0.0_JPRD
    ENDIF
    
    ZX0 = ZX
    PL(JGL) = REAL(ZX0,JPRD)
    
    ! catch overflow, should never happen
    IF( .NOT.(PW(JGL)==PW(JGL)) ) THEN
      WRITE(NOUT,*) 'invoke overflow ...PW ',KM, KN, JGL
      PW(JGL) = 0.0_JPRD
    ENDIF
    
! ++++++++++++++++++++++++++++++++++++++++++++++++
! <<<< END REFINEMENT >>>> 
! ++++++++++++++++++++++++++++++++++++++++++++++++

    ZF1 = SQRT(FQ(ZX0)/FP(ZX0))
    ZF2 = FQ1(ZX0)/FQ(ZX0)
    ZF3 = FP1(ZX0)/FP(ZX0)  

    ! continue to next root with refined ZX,ZR as initial condition
    ZH = -ZPIH/REAL(ISTEPMAX,JPRH)
    ZTHETA = 0.5_JPRH*ZPIH
    ZDX0 = -1._JPRH/(ZF1 + 0.25_JPRH*(ZF2 + ZF3)*SIN(2._JPRH*ZTHETA))
  ENDDO

ENDIF

CALL GSTATS(1650,1)
!     ------------------------------------------------------------------

!DIR$ IVDEP
!OCL NOVREC
DO JGL=1,KDGL/2
  ISYM = KDGL-JGL+1
  PL(ISYM) = -PL(JGL)
  PW(ISYM) = PW(JGL)
ENDDO

!     ------------------------------------------------------------------

!*      3. Diagnostics.
!          ------------

IF( LLOLD ) THEN

  IF(LLP2)THEN
    DO JGL=1,INS2
      ZM(JGL) = (ACOS(PL(JGL))-ACOS(ZLI(JGL)))*RA
      ZRR(JGL) = (ACOS(PL(JGL))-ACOS(ZREG(JGL)))*RA
      ZT(JGL) = ACOS(PL(JGL))*180._JPRD/ZPI
    ENDDO
  ENDIF
  
  IALLOW = 20
  DO JGL=1,INS2
  
    IF(LLP2)THEN
      WRITE(UNIT=NOUT,FMT=&
       &'('' M ='',I4,'' ROW ='',I4,'' ITERATIONS='',I4,'' ROOT='',F30.20,&
       &'' WEIGHT='',F30.20,'' MODIF :'',E8.2)')KM,JGL,ITER(JGL),PL(JGL)&
       &,PW(JGL),PL(JGL)-ZLI(JGL)
      WRITE(UNIT=NOUT,FMT=&
       &'(10X,'' LAST INC. : '',E8.2,'' MODIF IN M : '',F10.3,&
       &'' FROM THE REGULAR GRID : '',F10.3,'' COLAT '',F10.3)')&
       &ZMOD(JGL),ZM(JGL),ZRR(JGL),ZT(JGL)
    ENDIF

    IF(ITER(JGL) > IALLOW)THEN
      WRITE(UNIT=NOUT,FMT='('' CONVERGENCE FAILED IN SUGAW '')')
      WRITE(UNIT=NOUT,FMT='('' ALLOWED : '',I4,''&
       &NECESSARY : '',&
       &I4)')IALLOW,ITER(JGL)
      CALL ABORT_TRANS(' FAILURE IN SUGAW ')
    ENDIF

  ENDDO

ELSE

  IF(LLP2)THEN
    DO JGL=1,INS2
      WRITE(UNIT=NOUT,FMT=&
       &'('' M ='',I4,'' ROW ='',I4,'' ITERATIONS='',I4,'' ROOT='',F30.20,&
       &'' WEIGHT='',F30.20,'' COLAT '',F10.3)')KM,JGL,0,PL(JGL),PW(JGL),&
       & ACOS(PL(JGL))*180._JPRD/ZPIH
    ENDDO
  ENDIF

ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE SUGAW
END MODULE SUGAW_MOD
