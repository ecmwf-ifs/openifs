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

SUBROUTINE SUVERTFEB(&
  & YDVETA,YDDIMV,YDVFE,YDCVER,&
  & CDTAG,LDCOMP,LDINT_FROM_SURF,LDFIX_ORDER,&
  & KTYPE,KFLEV_IN,PETA_IN,KTYPE_IN,&
  & KFLEV_OUT,PETA_OUT,KTYPE_OUT,&
  & KTBC_IN,KBBC_IN,KTBC_OUT,KBBC_OUT,&
  & POPER)

!**** *SUVERTFEB*  - define VFE vertical operator (B-splines)

!**   Interface.
!     ----------

!     *CALL* SUVERTFEB

!     Explicit arguments :
!     --------------------
!      * INPUT:
!        CDTAG                     : name of the operator
!        LDCOMP                    : back compatibility
!        LDINT_FROM_SURF           : integral calculated from surface/top
!        LDFIX_ORDER               : fix order / fixed knots
!        KTYPE                     : TYPE OF OPERATOR:
!                                    KTYPE = -1 INTEGRAL OPERATOR
!                                    KTYPE =  1 FIRTS  ORDER DERIVATIVE OPERATOR
!                                    KTYPE =  2 SECOND ORDER DERIVATIVE OPERATOR
!        KFLEV_IN                  : NUMBER OF INPUT LEVELS OF OPERATOR
!        PETA_IN                   : VALUES OF INPUT ETA LEVELS
!        KTYPE_IN                  : 0 = value , 1=derivative, 2=second derivative
!        KFLEV_OUT                 : NUMBER OF OUTPUT LEVELS OF OPERATOR
!        PETA_OUT                  : VALUES OF OUTPUT ETA LEVELS
!        KTYPE_OUT                 : 0 = value , 1=derivative, 2=second derivative
!        KTBC_IN                   : TOP BOUNDARY TYPE CONDITION OF INPUT VECTOR
!        KBBC_IN                   : BOTTOM BOUNDARY TYPE CONDITION OF INPUT VECTOR
!        KTBC_OUT                  : TOP BOUNDARY TYPE CONDITION OF OUTPUT VECTOR
!        KBBC_OUT                  : BOTTOM BOUNDARY TYPE CONDITION OF OUTPUT VECTOR

!        Additional info about KTBC_IN to KBBC_OUT:
!        K[TBC|BBC]_[IN|OUT] =  1  -  value of function at relevant point is zero
!                            =  2  -  value of derivative at relevant point is zero
!                            =  0  -  no boundary condition assumed

!      * OUTPUT:
!        POPER                     : VFE OPERATOR (ACCORDING KTYPE)

!     Method.
!     -------

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ALADIN-NH documentation.

!     Author.
!     -------
!        Jozef Vivoda, SHMU/LACE 
!        Original : 2010-09

!     Modifications.
!     --------------
!      K. Yessad (July 2014): Move some variables.
!      J. Vivoda and P. Smolikova (Sep 2017): new options for VFE-NH
!     ------------------------------------------------------------------

USE YOMVERT  , ONLY : TVETA, TVFE
USE YOMDIMV  , ONLY : TDIMV
USE PARKIND1  ,ONLY : JPIM     ,JPRB    ,JPRD
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOMCVER  , ONLY : TCVER
USE YOMLUN   , ONLY : NULOUT
USE YOMMP0   , ONLY : LOUTPUT, NPRINTLEV
USE SUVFE_HLP, ONLY : SETGP2VFE, SETTM, RTMIN, RTMAX

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TVETA), INTENT(INOUT)     :: YDVETA
TYPE(TDIMV), INTENT(INOUT)     :: YDDIMV
TYPE(TVFE) , INTENT(INOUT)     :: YDVFE
TYPE(TCVER),INTENT(IN)         :: YDCVER
CHARACTER(LEN=15),INTENT(IN)   :: CDTAG
LOGICAL, INTENT(IN)            :: LDCOMP
LOGICAL, INTENT(IN)            :: LDINT_FROM_SURF
LOGICAL, INTENT(IN)            :: LDFIX_ORDER
INTEGER(KIND=JPIM), INTENT(IN) :: KTYPE
INTEGER(KIND=JPIM), INTENT(IN) :: KFLEV_IN
REAL   (KIND=JPRB), INTENT(IN) :: PETA_IN(KFLEV_IN)
INTEGER(KIND=JPIM), INTENT(IN) :: KTYPE_IN(KFLEV_IN)
INTEGER(KIND=JPIM), INTENT(IN) :: KFLEV_OUT
REAL   (KIND=JPRB), INTENT(IN) :: PETA_OUT(KFLEV_OUT)
INTEGER(KIND=JPIM), INTENT(IN) :: KTYPE_OUT(KFLEV_OUT)
INTEGER(KIND=JPIM), INTENT(IN) :: KTBC_IN(2)
INTEGER(KIND=JPIM), INTENT(IN) :: KBBC_IN(2)
INTEGER(KIND=JPIM), INTENT(IN) :: KTBC_OUT(2)
INTEGER(KIND=JPIM), INTENT(IN) :: KBBC_OUT(2)
REAL   (KIND=JPRD), INTENT(OUT):: POPER(KFLEV_OUT,KFLEV_IN)

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM)           :: JLEV, IJ
INTEGER(KIND=JPIM)           :: IORDER_IN, IORDER_OUT, IORDER_W
INTEGER(KIND=JPIM)           :: IBASIS_IN, IBASIS_OUT, IBASIS_W
INTEGER(KIND=JPIM)           :: IKNOT_IN , IKNOT_OUT , IKNOT_W
INTEGER(KIND=JPIM)           :: IOFF_IN  , IOFF_OUT  , IOFF_W
INTEGER(KIND=JPIM)           :: INTERNALS_IN, INTERNALS_OUT, INTERNALS_W
INTEGER(KIND=JPIM)           :: ITBC_W(2), IBBC_W(2)

!-------
! Z_IN(I,J) - I - i-th spline, J - j-th segment, k-th coefficient
!    - linear spline (i-th spline on j-th segment)
!    Spline_(i,j) = Z_IN(I,J,1) + Z_IN(I,J,2)*t 
!    - cubic spline 
!    Spline_(i,j) = Z_IN(I,J,1) + Z_IN(I,J,2)*t + Z_IN(I,J,3)*t*t + Z_IN(I,J,4)*t*t*t
!-------

! knots sequences
REAL(KIND=JPRB), ALLOCATABLE :: ZKNOT_IN  (:)
REAL(KIND=JPRB), ALLOCATABLE :: ZKNOT_OUT (:)
REAL(KIND=JPRB), ALLOCATABLE :: ZKNOT_W   (:)

! position of maximas of splines
REAL(KIND=JPRB), ALLOCATABLE :: Z_IN_ETAMAX  (:)
REAL(KIND=JPRB), ALLOCATABLE :: Z_OUT_ETAMAX (:)
REAL(KIND=JPRB), ALLOCATABLE :: Z_W_ETAMAX   (:)

! spline coefficients
REAL(KIND=JPRB), ALLOCATABLE :: Z_IN (:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: Z_OUT(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: Z_W  (:,:,:)

! local interval coordinate transformation matrices
REAL(KIND=JPRB), ALLOCATABLE :: ZTM_IN (:,:)
REAL(KIND=JPRB), ALLOCATABLE :: ZTM_OUT(:,:)
REAL(KIND=JPRB), ALLOCATABLE :: ZTM_W  (:,:)

! aux.
REAL(KIND=JPRB) :: ZDIFF

! VFE matrices
REAL(KIND=JPRB) :: ZA(KFLEV_IN,KFLEV_IN)
REAL(KIND=JPRB) :: ZB(KFLEV_IN,KFLEV_IN)
REAL(KIND=JPRB) :: ZAO(KFLEV_IN,KFLEV_IN)

! transformation matrix GP -> VFE
REAL(KIND=JPRB) :: ZGP2VFE(KFLEV_IN ,KFLEV_IN)

! transformation matrix VFE -> GP
REAL(KIND=JPRB) :: ZVFE2GP(KFLEV_OUT,KFLEV_IN)

! working array used for single/double precision conversion
REAL(KIND=JPRB) :: ZOPER(KFLEV_OUT,KFLEV_IN)

CHARACTER(LEN=240) :: CTAG

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "minv_caller.h"
#include "suvfe_matrix.intfb.h"
#include "suvfe_cpsplines.intfb.h"
#include "suvfe_knot.intfb.h"
#include "suvfe_implicitbc.intfb.h"
#include "suvfe_oper_setup.intfb.h"
#include "suvfe_testoper.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUVERTFEB',0,ZHOOK_HANDLE)
ASSOCIATE(NFLEVG=>YDDIMV%NFLEVG)
!     ------------------------------------------------------------------

WRITE(NULOUT,*) "----- SUVERTFEB ----------------"
IF( KTYPE == -1 )THEN
  WRITE(NULOUT,*) " VFE integral operator"
ELSEIF( KTYPE == 1 )THEN
  WRITE(NULOUT,*) " VFE first derivative operator"
ELSEIF( KTYPE == 2 )THEN
  WRITE(NULOUT,*) " VFE second derivative operator"
ELSE
  CALL ABOR1("ERROR IN SUVERTFEB: unknown operator type required ...")
ENDIF

! ------------------------------------------
! input basis functions setup
! ------------------------------------------
CALL SUVFE_OPER_SETUP(YDCVER,LDFIX_ORDER, KTBC_IN, KBBC_IN, KFLEV_IN, &
 & IORDER_IN, IBASIS_IN, INTERNALS_IN, IOFF_IN, IKNOT_IN)

! ------------------------------------------
! output basis functions 
! ------------------------------------------
CALL SUVFE_OPER_SETUP(YDCVER,LDFIX_ORDER, KTBC_OUT, KBBC_OUT, KFLEV_IN, &
 & IORDER_OUT, IBASIS_OUT, INTERNALS_OUT, IOFF_OUT, IKNOT_OUT)

! ------------------------------------------
! weighting functions (this is important for control of laplacian eigenvalues)
! ------------------------------------------
IF( LDFIX_ORDER .AND. .NOT.LDCOMP )THEN

  ! zero TBC for weighting function
  ! ITBC_W(1) = 0
  ! ITBC_W(2) = 0

  ! zero BBC for weighting function
  ! IBBC_W(1) = 0
  ! IBBC_W(2) = 0

  ! CALL SUVFE_OPER_SETUP(YDCVER,LDFIX_ORDER, ITBC_W, IBBC_W, KFLEV_IN, &
  !  & IORDER_W, IBASIS_W, INTERNALS_W, IOFF_W, IKNOT_W)

  ! pure galerkin 

  ITBC_W        = KTBC_OUT
  IBBC_W        = KBBC_OUT
  IORDER_W      = IORDER_OUT
  IBASIS_W      = IBASIS_OUT
  INTERNALS_W   = INTERNALS_OUT
  IKNOT_W       = IKNOT_OUT
  IOFF_W        = IOFF_OUT

  ! ITBC_W        = KTBC_IN
  ! IBBC_W        = KBBC_IN
  ! IORDER_W      = IORDER_IN
  ! IBASIS_W      = IBASIS_IN
  ! INTERNALS_W   = INTERNALS_IN
  ! IKNOT_W       = IKNOT_IN
  ! IOFF_W        = IOFF_IN

ELSE

  IORDER_W      = IORDER_OUT
  IBASIS_W      = IBASIS_OUT
  INTERNALS_W   = INTERNALS_OUT
  IKNOT_W       = IKNOT_OUT
  IBBC_W        = KBBC_OUT
  ITBC_W        = KTBC_OUT
  IOFF_W        = IOFF_OUT

ENDIF

!-------------------------
! transformation matrices (for coordinate transformations):
! this is done to minimize rounding errors
! during polynomial multiplications and other operations.
! If we perform operation on any polynomial above given 
! interval we always transform polynomial to be expressed
! in coordinate t=<0,1> on given interval.
!-------------------------
ALLOCATE( ZTM_IN(IORDER_IN,IORDER_IN) )
CALL SETTM(RTMIN, RTMAX, IORDER_IN,ZTM_IN)

ALLOCATE( ZTM_OUT(IORDER_OUT,IORDER_OUT) )
CALL SETTM(RTMIN, RTMAX, IORDER_OUT,ZTM_OUT)

ALLOCATE( ZTM_W(IORDER_W,IORDER_W) )
CALL SETTM(RTMIN, RTMAX, IORDER_W,ZTM_W)

!------------------------------------
!  compute knots
!------------------------------------
ALLOCATE( ZKNOT_IN (IKNOT_IN ) )
ALLOCATE( ZKNOT_OUT(IKNOT_OUT) )
ALLOCATE( ZKNOT_W  (IKNOT_W  ) )

CALL SUVFE_KNOT(YDVFE,YDCVER,LDFIX_ORDER, KTBC_IN, KBBC_IN, IBASIS_IN, IORDER_IN, &
 & KFLEV_IN, PETA_IN, ZKNOT_IN)

CALL SUVFE_KNOT(YDVFE,YDCVER,LDFIX_ORDER, KTBC_OUT, KBBC_OUT, IBASIS_OUT, IORDER_OUT, &
 & KFLEV_IN, PETA_IN, ZKNOT_OUT)
 
CALL SUVFE_KNOT(YDVFE,YDCVER,LDFIX_ORDER, ITBC_W, IBBC_W, IBASIS_W, IORDER_W, &
 & KFLEV_IN, PETA_IN, ZKNOT_W  )

! CALL SUVFE_KNOT(YDVFE,YDCVER,LDFIX_ORDER, KTBC_OUT, KBBC_OUT, IBASIS_OUT, IORDER_OUT, &
!  & KFLEV_OUT, PETA_OUT, ZKNOT_OUT)
!  
! CALL SUVFE_KNOT(YDVFE,YDCVER,LDFIX_ORDER, ITBC_W, IBBC_W, IBASIS_W, IORDER_W, &
!  & KFLEV_OUT, PETA_OUT, ZKNOT_W  )

!-------------------------
! compute basis functions
!-------------------------
ALLOCATE( Z_IN    (IBASIS_IN ,IORDER_IN ,IORDER_IN ))
ALLOCATE( Z_OUT   (IBASIS_OUT,IORDER_OUT,IORDER_OUT))
ALLOCATE( Z_W     (IBASIS_W  ,IORDER_W  ,IORDER_W  ))

IF (YDCVER%LVFE_VERBOSE) THEN
  WRITE(NULOUT,*) "SHAPE Z_W  :",SHAPE(Z_W  )
  WRITE(NULOUT,*) "SHAPE Z_IN :",SHAPE(Z_IN )
  WRITE(NULOUT,*) "SHAPE Z_OUT:",SHAPE(Z_OUT)
ENDIF

ALLOCATE( Z_IN_ETAMAX (IBASIS_IN ))
ALLOCATE( Z_OUT_ETAMAX(IBASIS_OUT))
ALLOCATE( Z_W_ETAMAX  (IBASIS_W)  )

CALL SUVFE_CPSPLINES(YDCVER,IORDER_IN ,IBASIS_IN ,ZTM_IN ,ZKNOT_IN ,Z_IN ,Z_IN_ETAMAX)
CALL SUVFE_CPSPLINES(YDCVER,IORDER_OUT,IBASIS_OUT,ZTM_OUT,ZKNOT_OUT,Z_OUT,Z_OUT_ETAMAX)
CALL SUVFE_CPSPLINES(YDCVER,IORDER_W  ,IBASIS_W  ,ZTM_W  ,ZKNOT_W  ,Z_W  ,Z_W_ETAMAX)

!-------------------------
! boundary conditions involvement
! this contains linear compbination of basis function
! to fullfill BCs of input and output function
!-------------------------
CALL SUVFE_IMPLICITBC(KTBC_IN, KBBC_IN, IORDER_IN, IBASIS_IN, Z_IN)
CALL SUVFE_IMPLICITBC(KTBC_OUT,KBBC_OUT,IORDER_OUT,IBASIS_OUT,Z_OUT)
CALL SUVFE_IMPLICITBC(ITBC_W,  IBBC_W,  IORDER_W,  IBASIS_W,  Z_W)

! -----
! loop over knot intervals that cover the domain
! (we skip zero size intervals)
! -----
ZA    = 0.0_JPRB
ZB    = 0.0_JPRB

!-----------------
! derivate operator
!-----------------
! mass matrix 
CALL SUVFE_MATRIX(YDCVER,LDINT_FROM_SURF,0, &
 & IKNOT_OUT,ZKNOT_OUT,IORDER_OUT,ZTM_OUT,IBASIS_OUT,IOFF_OUT,Z_OUT,  &
 & IKNOT_W  ,ZKNOT_W  ,IORDER_W  ,ZTM_W  ,IBASIS_W  ,IOFF_W  ,Z_W  ,  &
 & KFLEV_IN,ZA)

! stiffness matrix
CALL SUVFE_MATRIX(YDCVER,LDINT_FROM_SURF,KTYPE, &
 & IKNOT_IN,ZKNOT_IN,IORDER_IN,ZTM_IN,IBASIS_IN,IOFF_IN,Z_IN,  &
 & IKNOT_W ,ZKNOT_W ,IORDER_W ,ZTM_W ,IBASIS_W ,IOFF_W ,Z_W ,  &
 & KFLEV_IN,ZB)

! GP -> VFE (der)
CALL SETGP2VFE(YDCVER%LVFE_APPROX,YDCVER%RMINDETA,KFLEV_IN,KFLEV_IN,PETA_IN,KTYPE_IN,IORDER_IN, &
 & IBASIS_IN,Z_IN,IOFF_IN,IKNOT_IN,ZKNOT_IN,ZGP2VFE)

! VFE -> GP (der) (on output all values are derivatives or integrals)
CALL SETGP2VFE(.FALSE.,YDCVER%RMINDETA,KFLEV_IN,KFLEV_OUT,PETA_OUT,KTYPE_OUT,IORDER_OUT, &
 & IBASIS_OUT,Z_OUT,IOFF_OUT,IKNOT_OUT,ZKNOT_OUT,ZVFE2GP)

! check diagonal dominance of transformation matrix
! WRITE(NULOUT,*) "ZGP2VFE:: diagonal dominance ", TRIM(CDTAG)
! DO JLEV=1,KFLEV_IN
!   DO IJ=1,KFLEV_IN
!     IF(JLEV /= IJ)THEN
!       IF(ABS(ZGP2VFE(JLEV,JLEV)) < ABS(ZGP2VFE(JLEV,IJ)))THEN
!, FT2X         WRITE(NULOUT,*) "PROBLEM WITH DIAGONAL DOMINANCE ", TRIM(CDTAG)
!         ZDIFF = ABS(ZGP2VFE(JLEV,JLEV)) - ABS(ZGP2VFE(JLEV,IJ))
!         WRITE(NULOUT,'(2I7,F15.9,E15.5)') JLEV, IJ, ABS(ZGP2VFE(JLEV,JLEV)), ZDIFF
!       ENDIF
!     ENDIF
!   ENDDO
! ENDDO

IF (YDCVER%LVFE_VERBOSE) THEN
  WRITE(NULOUT,*) TRIM(CDTAG), " PETA_IN "
  DO JLEV=1,KFLEV_IN
    WRITE(NULOUT,'(A,1X,I5,1X,300F12.7)') TRIM(CDTAG),JLEV,PETA_IN(JLEV)
  ENDDO

  WRITE(NULOUT,*) TRIM(CDTAG), " ZKNOT_IN "
  DO JLEV=1,IKNOT_IN
    WRITE(NULOUT,'(A,1X,I5,1X,300F12.7)') TRIM(CDTAG),JLEV,ZKNOT_IN(JLEV)
  ENDDO

  WRITE(CTAG,*) TRIM(CDTAG), " ZGP2VFE:: GP -> VFE (der)"

  CALL PRINT_MATRIX(CTAG, KFLEV_IN, KFLEV_IN, ZGP2VFE)

  WRITE(CTAG,*) TRIM(CDTAG), " ZA :: MASS MATRIX to be inverted"
  CALL PRINT_MATRIX(CTAG, KFLEV_IN, KFLEV_IN, ZA)

  WRITE(CTAG,*) TRIM(CDTAG), " ZB:: STIFF matrix"
  CALL PRINT_MATRIX(CTAG, KFLEV_IN, KFLEV_IN, ZB)

  WRITE(CTAG,*) TRIM(CDTAG), " ZVFE2GP:: VFE -> GP (der)"
  CALL PRINT_MATRIX(CTAG, KFLEV_OUT, KFLEV_IN, ZVFE2GP)
ENDIF

WRITE(NULOUT,*) "VFE - Inversion of mass matrix MA"

! ensure that ZA is symmetric
! DO JLEV = 1, KFLEV_IN
!   DO IJ = JLEV +1, KFLEV_IN
!     ZDIFF = 0.5_JPRB * (ZA(JLEV, IJ) + ZA(IJ, JLEV))
!     ZA(JLEV, IJ  ) = ZDIFF
!     ZA(IJ  , JLEV) = ZDIFF
!   ENDDO
! ENDDO

CALL MINV_CALLER(.TRUE.,KFLEV_IN,ZA,ZAO)
ZA(:,:)=ZAO(:,:)

WRITE(NULOUT,*) "VFE - Inversion of transform matrix ZGP2VFE"

CALL MINV_CALLER(.TRUE.,KFLEV_IN,ZGP2VFE,ZAO)
ZGP2VFE(:,:)=ZAO(:,:)

IF (YDCVER%LVFE_VERBOSE) THEN
  WRITE(CTAG,*) TRIM(CDTAG), " ZGP2VFE^(-1):: GP -> VFE (der)"
  CALL PRINT_MATRIX(CTAG, KFLEV_IN, KFLEV_IN, ZGP2VFE)

  WRITE(CTAG,*) TRIM(CDTAG), " ZA:: MASS matrix (inverted)"
  CALL PRINT_MATRIX(CTAG, KFLEV_IN, KFLEV_IN, ZA)
ENDIF

! operator
POPER = MATMUL(ZVFE2GP,MATMUL(ZA,MATMUL(ZB,ZGP2VFE)))

IF (YDCVER%LVFE_VERBOSE) THEN
  IF( KTYPE == -1 )THEN
    WRITE(CTAG,*) TRIM(CDTAG), " POPER::VFE integral operator"
  ELSEIF( KTYPE == 1 )THEN
    WRITE(CTAG,*) TRIM(CDTAG), " POPER::VFE first derivative operator"
  ELSEIF( KTYPE == 2 )THEN
    WRITE(CTAG,*) TRIM(CDTAG), " POPER::VFE second derivative operator"
  ENDIF
  ZOPER(:,:)=REAL(POPER(:,:),JPRB)
  CALL PRINT_MATRIX(CTAG, KFLEV_OUT, KFLEV_IN, ZOPER)
ENDIF

! ------------------------
! Test of operator quality
! ------------------------

IF (YDCVER%LVFE_VERBOSE) THEN
  CALL SUVFE_TESTOPER(CDTAG,'CST',LDINT_FROM_SURF,KTYPE,KFLEV_IN, &
   & PETA_IN,KTYPE_IN,KFLEV_OUT,PETA_OUT,POPER)
  CALL SUVFE_TESTOPER(CDTAG,'LIN',LDINT_FROM_SURF,KTYPE,KFLEV_IN, &
   & PETA_IN,KTYPE_IN,KFLEV_OUT,PETA_OUT,POPER)
  CALL SUVFE_TESTOPER(CDTAG,'QUA',LDINT_FROM_SURF,KTYPE,KFLEV_IN, &
   & PETA_IN,KTYPE_IN,KFLEV_OUT,PETA_OUT,POPER)
  CALL SUVFE_TESTOPER(CDTAG,'SIN',LDINT_FROM_SURF,KTYPE,KFLEV_IN, &
   & PETA_IN,KTYPE_IN,KFLEV_OUT,PETA_OUT,POPER)
ENDIF

!-------------------------
! Deallocate arrays
!-------------------------

IF (ALLOCATED(ZKNOT_IN))    DEALLOCATE(ZKNOT_IN)
IF (ALLOCATED(ZKNOT_OUT))   DEALLOCATE(ZKNOT_OUT)
IF (ALLOCATED(ZKNOT_W))     DEALLOCATE(ZKNOT_W)
IF (ALLOCATED(Z_IN))        DEALLOCATE(Z_IN)
IF (ALLOCATED(Z_OUT))       DEALLOCATE(Z_OUT)
IF (ALLOCATED(Z_W))         DEALLOCATE(Z_W)
IF (ALLOCATED(ZTM_IN))      DEALLOCATE(ZTM_IN)
IF (ALLOCATED(ZTM_OUT))     DEALLOCATE(ZTM_OUT)
IF (ALLOCATED(ZTM_W))       DEALLOCATE(ZTM_W)
IF (ALLOCATED(Z_IN_ETAMAX)) DEALLOCATE(Z_IN_ETAMAX)
IF (ALLOCATED(Z_OUT_ETAMAX))DEALLOCATE(Z_OUT_ETAMAX)
IF (ALLOCATED(Z_W_ETAMAX))  DEALLOCATE(Z_W_ETAMAX)

!     ------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUVERTFEB',1,ZHOOK_HANDLE)
!     ------------------------------------------------------------------
END SUBROUTINE SUVERTFEB

!---------------------

SUBROUTINE PRINT_MATRIX(CDTAG, NX, NY, PM)

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMLUN   , ONLY : NULOUT

IMPLICIT NONE

CHARACTER(LEN=240)  ,INTENT(IN) :: CDTAG
INTEGER(KIND=JPIM), INTENT(IN) :: NX
INTEGER(KIND=JPIM), INTENT(IN) :: NY
REAL   (KIND=JPRB), INTENT(IN) :: PM(NX, NY)

INTEGER(KIND=JPIM) :: I, J, K

! K - how much elements to print
K = NX

WRITE(NULOUT,*) TRIM(CDTAG), " LEFT UPPER CORNER"
CALL FLUSH(NULOUT)
DO I=1, MIN(K,NX)
   WRITE(NULOUT,'(I5,100(1X,G12.5))') I,(PM(I,J),J=1, MIN(K,NY))
   CALL FLUSH(NULOUT)
ENDDO

WRITE(NULOUT,*) TRIM(CDTAG), " RIGHT BOTTOM CORNER"
CALL FLUSH(NULOUT)
DO I=MAX(1,NX - K + 1), NX
   WRITE(NULOUT,'(I5,100(1X,G12.5))') I,(PM(I,J),J=MAX(1,NY - K + 1),NY)
   CALL FLUSH(NULOUT)
ENDDO

END SUBROUTINE PRINT_MATRIX

