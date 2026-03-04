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

MODULE ORDER_INDEPENDENT_SUMMATION_MOD

!**** ORDER_INDEPENDENT_SUMMATION_MOD

!     Purpose.
!     --------
!     Functions to perform global (over all processors) and 
!     local (per-processor) order-independent accurate summation
!     and order-independent inner products.

!**   Interface.
!     ----------

!        result = ORDER_INDEP_GLOBAL_SUM (P1)
!        result = ORDER_INDEP_LOCAL_SUM (P1)
!        result = ORDER_INDEP_DOT_PRODUCT (P1,P2,PW)

!        Input required arguments :
!        -------------------------
!           P1       -  A 1d array of KIND=JPRB reals
!           P2       -  A 1d array of KIND=JPRB reals. Same length as P1.

!        Input optional arguments :
!        -------------------------
!           PW       -  A 1d array of KIND=JPRB reals. Same length as P1.
!                       If specified. ORDER_INDEP_DOT_PRODUCT returns
!                       a weighted inner product. PW defines the weights.
!                       If not specified, ORDER_INDEP_DOT_PRODUCT returns
!                       an unweighted dot product of P and P2.
!
!           KNG (global sum only) - Global length of input array.
!
!           LD_ABORT_IFNOT_REPROD - Abort if results are not guaranteed
!                                   to be bit-reproducible. (Default=T)
!
!           LD_OPENMP             - .TRUE. => use openMP (Default is
!                                   .TRUE. for global sum, .FALSE. for
!                                   local sum.)

!        Output required arguments :
!        -------------------------
!           none

!     Author.
!     -------
!        Mike Fisher  ECMWF

!     Modifications.
!     --------------
!        Original: 2006-20-22

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE COMPENSATED_SUMMATION_MOD, ONLY : COMPENSATED_SUM, COMPENSATED_SUM_OMP, &
                                & COMPENSATED_DOT_PRODUCT, &
                                & COMPENSATED_DOT_PRODUCT_OMP
USE MPL_MODULE, ONLY : MPL_ALLREDUCE, MPL_ALLGATHERV, MPL_MYRANK, MPL_NPROC, &
                     & MPL_MESSAGE, MPL_SEND, MPL_RECV, MPL_WAIT, &
                     & JP_NON_BLOCKING_STANDARD
USE YOMHOOK   , ONLY : LHOOK,   DR_HOOK, JPHOOK

SAVE
PRIVATE
PUBLIC ORDER_INDEP_LOCAL_SUM, &
     & ORDER_INDEP_GLOBAL_SUM, &
     & ORDER_INDEP_GLOBAL_SUM2, &
     & ORDER_INDEP_ALLREDUCE, &
     & ORDER_INDEP_DOT_PRODUCT

INTERFACE ORDER_INDEP_LOCAL_SUM
  MODULE PROCEDURE ORDER_INDEP_LOCAL_SUM
END INTERFACE

INTERFACE ORDER_INDEP_GLOBAL_SUM
  MODULE PROCEDURE ORDER_INDEP_GLOBAL_SUM
END INTERFACE

INTERFACE ORDER_INDEP_GLOBAL_SUM2
  MODULE PROCEDURE ORDER_INDEP_GLOBAL_SUM2
END INTERFACE

INTERFACE ORDER_INDEP_ALLREDUCE
  MODULE PROCEDURE ORDER_INDEP_ALLREDUCE
END INTERFACE

INTERFACE ORDER_INDEP_DOT_PRODUCT
  MODULE PROCEDURE ORDER_INDEP_DOT_PRODUCT
END INTERFACE

CONTAINS

FUNCTION ORDER_INDEP_LOCAL_SUM (PIN,LD_ABORT_IFNOT_REPROD,LD_OPENMP)

!-----------------------------------------------------------------
!  Returns an accurate local (i.e. on a single processor) sum of
!  the elements of PIN. The sum is bit-reproducible for any
!  ordering of the elements of PIN.
!
!  NB: PIN is unmodified on return
!
! Algorithm:
! ----------
!
!  The algorithm is based on Ogita et al. (2005) SIAM J. Sci. Computing,
!  Vol.26, No.6, pp1955-1988. This is based in turn on an algorithm
!  by Knuth (1969, seminumerical algorithms).
!
!  This version iterates the compensated sum algorithm until the
!  result is guaranteed to be within 4*eps of the true sum. It
!  then rounds the result to the nearest floating-point number
!  whose last three bits are zero, thereby guaranteeing an
!  order-independent result.
!
!  Author: Mike Fisher ECMWF 2006/02/08
!
!-----------------------------------------------------------------

IMPLICIT NONE

REAL(KIND=JPRB) :: ORDER_INDEP_LOCAL_SUM
REAL(KIND=JPRB), INTENT(IN) :: PIN(:)
LOGICAL,OPTIONAL,INTENT(IN) :: LD_ABORT_IFNOT_REPROD, LD_OPENMP

INTEGER(KIND=JPIM) :: IN
REAL(KIND=JPRB) :: ZCORR,ZERR,ZOLDERR,ZBETA,ZRES
REAL(KIND=JPRB), ALLOCATABLE :: ZP(:)
LOGICAL :: LLABORT, LL_OPENMP
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

INTEGER(KIND=JPIM), SAVE :: INMSG=0

IF (LHOOK) CALL DR_HOOK ('ORDER_INDEPENDENT_SUMMATION_MOD:ORDER_INDEP_LOCAL_SUM', &
                      &  0,ZHOOK_HANDLE)

IF (PRESENT(LD_ABORT_IFNOT_REPROD)) THEN
  LLABORT = LD_ABORT_IFNOT_REPROD
ELSE
  LLABORT = .TRUE.
ENDIF

IF (PRESENT(LD_OPENMP)) THEN
  LL_OPENMP = LD_OPENMP
ELSE
  LL_OPENMP = .FALSE.
ENDIF

IN = SIZE(PIN)

IF (REAL(2*IN,JPRB)*EPSILON(ZRES) >= 1.0) THEN
  CALL MPL_MESSAGE (CDMESSAGE='n is too large to guarantee error bounds', &
                    & CDSTRING='ORDER_INDEP_LOCAL_SUM',LDABORT=.TRUE.)
ENDIF

TEST_ARRAY_LENGTH: IF (IN>0) THEN
  ZOLDERR = HUGE(ZERR)

!--- Copy the input array. This avoids some tricky indexing, at the
!--- expense of some inefficency.

  ALLOCATE (ZP(IN))
  ZP(:) = PIN(:)

  K_LOOP: DO

!--- transform local arrays

    IF (LL_OPENMP) THEN
      CALL COMPENSATED_SUM_OMP (ZP,IN,ZCORR,ZERR)
    ELSE
      CALL COMPENSATED_SUM (ZP,IN,ZCORR,ZERR)
    ENDIF

!--- Calculate final result

    ZRES = ZP(IN) + ZCORR

!--- Calculate error bound. This is corollary 4.7 from Ogita et al. (2005)

    ZBETA = ZERR*(REAL(2*IN,JPRB)*EPSILON(ZRES)) &
         & /(1.0_JPRB - REAL(2*IN,JPRB)*EPSILON(ZRES))

    ZERR = EPSILON(ZRES)*ABS(ZRES) &
       & +(ZBETA + ( 2.0_JPRB*EPSILON(ZRES)*EPSILON(ZRES)*ABS(ZRES) &
       &            +3.0_JPRB*TINY(ZRES)))

!--- exit if the error is small enough

    IF (ZERR<4.0_JPRB*SPACING(ZRES)) EXIT K_LOOP

!--- Take appropriate action if ZRES cannot be sufficiently refined.

    IF (ZERR >= ZOLDERR) THEN
      INMSG=INMSG+1

      IF (INMSG<=100) THEN
        CALL MPL_MESSAGE ( &
            & CDMESSAGE= 'ORDER_INDEP_LOCAL_SUM: FALIED TO REFINE SUM', &
            & CDSTRING='ORDER_INDEP_LOCAL_SUM')
        CALL MPL_MESSAGE ( &
            & CDMESSAGE='WARNING: POSSIBLITY OF NON-REPRODUCIBLE RESULTS',&
            & CDSTRING='ORDER_INDEP_LOCAL_SUM')
      ENDIF

      IF (INMSG==100) THEN
        CALL MPL_MESSAGE ( &
            & CDMESSAGE='ORDER_INDEP_LOCAL_SUM: INMSG>100. OUTPUT SUPPRESSED',&
            & CDSTRING='ORDER_INDEP_LOCAL_SUM')
      ENDIF

      IF (LLABORT) THEN
        CALL MPL_MESSAGE (CDMESSAGE= &
                        & 'ABORT BECAUSE LD_ABORT_IFNOT_REPROD WAS SET', &
                        & CDSTRING='ORDER_INDEP_LOCAL_SUM',LDABORT=.TRUE.)
      ENDIF
    ENDIF

    ZOLDERR = ZERR

  ENDDO K_LOOP

!--- At this stage, we have guaranteed that ZRES is less than 4*EPS
!--- away from the exact sum. There are only eight floating point
!--- numbers in this range. So, if we find the nearest number that
!--- has its last three bits zero, then we have a reproducible result.

  ORDER_INDEP_LOCAL_SUM = ROUND (ZRES)

  DEALLOCATE (ZP)
ELSE TEST_ARRAY_LENGTH

  ORDER_INDEP_LOCAL_SUM = 0.0_JPRB

ENDIF TEST_ARRAY_LENGTH

IF (LHOOK) CALL DR_HOOK ('ORDER_INDEPENDENT_SUMMATION_MOD:ORDER_INDEP_LOCAL_SUM', &
                      &  1,ZHOOK_HANDLE)

END FUNCTION ORDER_INDEP_LOCAL_SUM

FUNCTION ORDER_INDEP_GLOBAL_SUM (PIN,KNG,LD_ABORT_IFNOT_REPROD, LD_OPENMP)

!-----------------------------------------------------------------
!
!  Returns an accurate global sum of the elements of PIN. The
!  sum is bit-reproducible for any distribution of PIN over
!  threads and tasks, and is independent of the ordering of the
!  elements of PIN.
!
!  NB: PIN is unmodified on return
!
! Arguments:
! ----------
!
! Required:
!
!   PIN                  - INTENT(IN) - The array to be summed.
!
!
! Optional:
!
!   KNG                  - INTENT(IN) - Global length of array.
!                 
!   LD_ABORT_IFNOT_REPROD - INTENT(IN) - Defines behaviour in case
!                                        a reproducible result cannot
!                                        be guaranteed.
!
!   LD_OPENMP            - INTENT(IN) - Use OpenMP parallelization.
!              
!
! Algorithm:
! ----------
!
!  The algorithm is based on Ogita et al. (2005) SIAM J. Sci. Computing,
!  Vol.26, No.6, pp1955-1988. This is based in turn on an algorithm
!  by Knuth (1969, seminumerical algorithms).
!
!  This version adds a second layer of parallelism on top of that
!  provided by COMPENSATED_SUM_OMP. It iterates the compensated 
!  summation until the result is guaranteed to be within 4*eps
!  of the true sum. It then rounds the result to the nearest
!  floating-point number whose last three bits are zero, thereby
!  guaranteeing an order-independent result.
!
!  Author: Mike Fisher ECMWF 2006/02/08
!
!-----------------------------------------------------------------

IMPLICIT NONE

REAL(KIND=JPRB) :: ORDER_INDEP_GLOBAL_SUM

REAL(KIND=JPRB),             INTENT(IN) :: PIN(:)
INTEGER(KIND=JPIM),OPTIONAL, INTENT(IN) :: KNG
LOGICAL,           OPTIONAL, INTENT(IN) :: LD_ABORT_IFNOT_REPROD, LD_OPENMP

INTEGER(KIND=JPIM) :: J,IN,ING,INPROC
REAL(KIND=JPRB) :: ZCORR,ZERR,ZOLDERR,ZBUFFL(3),ZBETA,ZRES
REAL(KIND=JPRB), ALLOCATABLE :: ZPSUMS(:),ZPERRS(:),ZPCORS(:), &
                              & ZBUFFG(:),ZP(:)
INTEGER(KIND=JPIM), ALLOCATABLE :: IRECVCOUNTS(:)
LOGICAL :: LLABORT, LL_OPENMP
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

INTEGER(KIND=JPIM), SAVE :: INMSG=0

IF (LHOOK) CALL DR_HOOK ('ORDER_INDEPENDENT_SUMMATION_MOD:ORDER_INDEP_GLOBAL_SUM', &
                      &  0,ZHOOK_HANDLE)

INPROC = MPL_NPROC()

IF (PRESENT(LD_ABORT_IFNOT_REPROD)) THEN
  LLABORT = LD_ABORT_IFNOT_REPROD
ELSE
  LLABORT = .TRUE.
ENDIF

IF (PRESENT(LD_OPENMP)) THEN
  LL_OPENMP = LD_OPENMP
ELSE
  LL_OPENMP = .TRUE.
ENDIF

IN = SIZE(PIN)

!--- global length of vector (needed for error bound calculation)

IF (.NOT.PRESENT(KNG)) THEN
  ING = IN
  IF (INPROC>1) THEN
    CALL MPL_ALLREDUCE (ING,'SUM',CDSTRING='ORDER_INDEP_GLOBAL_SUM')
  ENDIF
ELSE
  ING = KNG
  IF (KNG<IN) THEN
    CALL MPL_MESSAGE (CDMESSAGE='Specified KNG < SIZE(PIN)', &
                    & CDSTRING='ORDER_INDEP_GLOBAL_SUM',LDABORT=.TRUE.)
  ENDIF
ENDIF

IF (REAL(2*ING,JPRB)*EPSILON(ZRES) >= 1.0) THEN
  CALL MPL_MESSAGE (CDMESSAGE='n is too large to guarantee error bounds', &
                  & CDSTRING='ORDER_INDEP_GLOBAL_SUM',LDABORT=.TRUE.)
ENDIF

ALLOCATE (ZP(MAX(IN,1_JPIM)))
ALLOCATE (ZBUFFG(INPROC*SIZE(ZBUFFL)))
ALLOCATE (ZPSUMS(INPROC))
ALLOCATE (ZPERRS(INPROC))
ALLOCATE (ZPCORS(INPROC))
ALLOCATE (IRECVCOUNTS(INPROC))

ZOLDERR = HUGE(ZERR)

!--- Copy the input array. This avoids some tricky indexing, at the
!--- expense of some inefficency.

IF (IN>0) THEN
  ZP(:) = PIN(:)
ELSE
  ZP(1) = 0.0_JPRB
ENDIF

K_LOOP: DO

!--- transform local arrays

  IF (IN>0) THEN
    IF (LL_OPENMP) THEN
      CALL COMPENSATED_SUM_OMP (ZP,IN,ZCORR,ZERR)
    ELSE
      CALL COMPENSATED_SUM (ZP,IN,ZCORR,ZERR)
    ENDIF
  ENDIF

!--- gather partial sums and error bounds to all processors

  ZBUFFL(1) = ZP(MAX(IN,1_JPIM))

  IF (IN>0) THEN
    ZBUFFL(2) = ZERR
    ZBUFFL(3) = ZCORR
  ELSE
    ZBUFFL(2) = 0.0_JPRB
    ZBUFFL(3) = 0.0_JPRB
  ENDIF

  IF (INPROC>1) THEN
!-- could use MPL_ALLGATHER here, if it existed!

    IRECVCOUNTS(:) = SIZE(ZBUFFL)

    CALL MPL_ALLGATHERV (ZBUFFL,ZBUFFG,IRECVCOUNTS, &
                       & CDSTRING='ORDER_INDEP_GLOBAL_SUM')

    DO J=1,INPROC
      ZPSUMS(J) = ZBUFFG(1+(J-1)*SIZE(ZBUFFL))
      ZPERRS(J) = ZBUFFG(2+(J-1)*SIZE(ZBUFFL))
      ZPCORS(J) = ZBUFFG(3+(J-1)*SIZE(ZBUFFL))
    ENDDO
  ELSE
    ZPSUMS(1) = ZBUFFL(1)
    ZPERRS(1) = ZBUFFL(2)
    ZPCORS(1) = ZBUFFL(3)
  ENDIF

!--- transform partial sums

  CALL COMPENSATED_SUM (ZPSUMS,INPROC,ZCORR,ZERR)
  ZERR  = ZERR  + SUM(ZPERRS)
  ZCORR = ZCORR + SUM(ZPCORS)

!--- Calculate final result

  ZRES = ZPSUMS(INPROC) + ZCORR

!--- Calculate error bound. This is corollary 4.7 from Ogita et al. (2005)

  ZBETA = ZERR*(REAL(2*ING,JPRB)*EPSILON(ZRES)) &
       & /(1.0_JPRB - REAL(2*ING,JPRB)*EPSILON(ZRES))

  ZERR = EPSILON(ZRES)*ABS(ZRES) &
     & +(ZBETA + ( 2.0_JPRB*EPSILON(ZRES)*EPSILON(ZRES)*ABS(ZRES) &
     &            +3.0_JPRB*TINY(ZRES)))

!--- update the last element of the local array

  ZP(MAX(IN,1_JPIM)) = ZPSUMS(MPL_MYRANK())

!--- exit if the global error is small enough

  IF (ZERR<4.0_JPRB*SPACING(ZRES)) EXIT K_LOOP

!--- Take appropriate action if ZRES cannot be sufficiently refined.

  IF (ZERR >= ZOLDERR) THEN
    INMSG=INMSG+1

    IF (INMSG<=100) THEN
      CALL MPL_MESSAGE ( &
          & CDMESSAGE= 'ORDER_INDEP_GLOBAL_SUM: FALIED TO REFINE SUM', &
          & CDSTRING='ORDER_INDEP_GLOBAL_SUM')
      CALL MPL_MESSAGE ( &
          & CDMESSAGE='WARNING: POSSIBLITY OF NON-REPRODUCIBLE RESULTS',&
          & CDSTRING='ORDER_INDEP_GLOBAL_SUM')
    ENDIF

    IF (INMSG==100) THEN
      CALL MPL_MESSAGE ( &
          & CDMESSAGE='ORDER_INDEP_GLOBAL_SUM: INMSG>100. OUTPUT SUPPRESSED',&
          & CDSTRING='ORDER_INDEP_GLOBAL_SUM')
    ENDIF

    IF (LLABORT) THEN
      CALL MPL_MESSAGE (CDMESSAGE= &
                      & 'ABORT BECAUSE LD_ABORT_IFNOT_REPROD WAS SET', &
                      & CDSTRING='ORDER_INDEP_GLOBAL_SUM',LDABORT=.TRUE.)
    ENDIF
  ENDIF

  ZOLDERR = ZERR

ENDDO K_LOOP

!--- At this stage, we have guaranteed that ZRES less than 4*EPS
!--- away from the exact sum. There are only four floating point
!--- numbers in this range. So, if we find the nearest number that
!--- has its last three bits zero, then we have a reproducible result.

ORDER_INDEP_GLOBAL_SUM = ROUND (ZRES)

DEALLOCATE (IRECVCOUNTS)
DEALLOCATE (ZPCORS)
DEALLOCATE (ZPERRS)
DEALLOCATE (ZPSUMS)
DEALLOCATE (ZBUFFG)
DEALLOCATE (ZP)

IF (LHOOK) CALL DR_HOOK ('ORDER_INDEPENDENT_SUMMATION_MOD:ORDER_INDEP_GLOBAL_SUM', &
                      &  1,ZHOOK_HANDLE)

END FUNCTION ORDER_INDEP_GLOBAL_SUM

SUBROUTINE ORDER_INDEP_GLOBAL_SUM2 (PIN,POUT,KNVEC,KDIM,KNL,LD_ABORT_IFNOT_REPROD,LD_OPENMP)

!-----------------------------------------------------------------
!
!  This is a vector version of ORDER_INDEP_GLOBAL_SUM, which 
!  returns a vector of accurate global sums of the elements of matrix
!  PIN along dimension KDIM. The  sum is bit-reproducible for any
!  distribution of PIN over  threads and tasks, and is independent of
!  the ordering of the elements of PIN.
!
!  NB: PIN is unmodified on return
!
! Arguments:
! ----------
!
! Required:
!
!   PIN                  - INTENT(IN)  - The array to be summed.
!
!   POUT                 - INTENT(OUT) - The vector with sums
!
!   KDIM                 - INTENT(IN)  - The dimension to sum along.
!
!   KNL                  - INTENT(IN)  - Local lengths of vectors.
!                 
!
! Optional:
!
!   LD_ABORT_IFNOT_REPROD - INTENT(IN) - Defines behaviour in case
!                                        a reproducible result cannot
!                                        be guaranteed.
!
!   LD_OPENMP            - INTENT(IN)  - Use OpenMP parallelization.
!              
!
! Algorithm:
! ----------
!
!  The algorithm is based on Ogita et al. (2005) SIAM J. Sci. Computing,
!  Vol.26, No.6, pp1955-1988. This is based in turn on an algorithm
!  by Knuth (1969, seminumerical algorithms).
!
!  This version adds a second layer of parallelism on top of that
!  provided by COMPENSATED_SUM_OMP. It iterates the compensated 
!  summation until the result is guaranteed to be within 4*eps
!  of the true sum. It then rounds the result to the nearest
!  floating-point number whose last three bits are zero, thereby
!  guaranteeing an order-independent result.
!
!  Author: Tomas Wilhelmsson ECMWF 2010/03/30
!
!-----------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),           INTENT(IN)  :: KNVEC
REAL(KIND=JPRB),              INTENT(IN)  :: PIN(:,:)
REAL(KIND=JPRB),              INTENT(OUT) :: POUT(KNVEC)
INTEGER(KIND=JPIM),           INTENT(IN)  :: KDIM
INTEGER(KIND=JPIM),           INTENT(IN)  :: KNL(KNVEC)
LOGICAL,            OPTIONAL, INTENT(IN)  :: LD_ABORT_IFNOT_REPROD, LD_OPENMP

INTEGER(KIND=JPIM) :: J,JL,JP,IBUFLEN,INVEC,INPROC,ING(KNVEC)
REAL(KIND=JPRB), DIMENSION(KNVEC) :: ZCORR,ZERR,ZOLDERR,ZBETA,ZRES
REAL(KIND=JPRB), ALLOCATABLE :: ZPSUMS(:,:),ZPERRS(:,:),ZPCORS(:,:), &
                              & ZBUFFL(:),ZBUFFG(:),ZP(:,:)
INTEGER(KIND=JPIM), ALLOCATABLE :: IRECVCOUNTS(:)
LOGICAL :: LLABORT, LL_OPENMP, LLDONE(KNVEC)
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

INTEGER(KIND=JPIM), SAVE :: INMSG=0

IF (LHOOK) CALL DR_HOOK ('ORDER_INDEPENDENT_SUMMATION_MOD:ORDER_INDEP_GLOBAL_SUM2', &
                      &  0,ZHOOK_HANDLE)

IF (KDIM<1 .OR. KDIM>2) THEN
  CALL MPL_MESSAGE (CDMESSAGE='Invalid KDIM value', &
                  & CDSTRING='ORDER_INDEP_GLOBAL_SUM2',LDABORT=.TRUE.)
ENDIF

INPROC = MPL_NPROC()

IF (PRESENT(LD_ABORT_IFNOT_REPROD)) THEN
  LLABORT = LD_ABORT_IFNOT_REPROD
ELSE
  LLABORT = .TRUE.
ENDIF

IF (PRESENT(LD_OPENMP)) THEN
  LL_OPENMP = LD_OPENMP
ELSE
  LL_OPENMP = .TRUE.
ENDIF

!--- global lengths of vectors (needed for error bound calculation)

ING(:) = KNL(:)
IF (INPROC>1) THEN
  CALL MPL_ALLREDUCE (ING,'SUM',CDSTRING='ORDER_INDEP_GLOBAL_SUM2')
ENDIF

IF (ANY(REAL(2*ING(:),JPRB)*EPSILON(ZRES) >= 1.0)) THEN
  CALL MPL_MESSAGE (CDMESSAGE='n is too large to guarantee error bounds', &
                  & CDSTRING='ORDER_INDEP_GLOBAL_SUM2',LDABORT=.TRUE.)
ENDIF

IBUFLEN=3

ALLOCATE (ZP(MAX(MAXVAL(KNL),1_JPIM),KNVEC))
ALLOCATE (ZBUFFL(IBUFLEN*KNVEC))
ALLOCATE (ZBUFFG(INPROC*IBUFLEN*KNVEC))
ALLOCATE (ZPSUMS(INPROC,KNVEC))
ALLOCATE (ZPERRS(INPROC,KNVEC))
ALLOCATE (ZPCORS(INPROC,KNVEC))
ALLOCATE (IRECVCOUNTS(INPROC))

ZOLDERR(:) = HUGE(ZERR)
LLDONE(:)  = .FALSE. 

!--- Copy the input array. This avoids some tricky indexing, at the
!--- expense of some inefficency.


DO J=1,KNVEC
  IF (KNL(J)>0) THEN
    IF (KDIM==1) ZP(1:KNL(J),J) = PIN(1:KNL(J),J) 
    IF (KDIM==2) ZP(1:KNL(J),J) = PIN(J,1:KNL(J)) 
  ELSE
    ZP(1,J) = 0.0_JPRB
  ENDIF
ENDDO
  
K_LOOP: DO

!--- transform local arrays

  JL=0
  DO J=1,KNVEC
    IF (KNL(J)>0 .AND. .NOT. LLDONE(J)) THEN
      IF (LL_OPENMP) THEN
        CALL COMPENSATED_SUM_OMP (ZP(:,J),KNL(J),ZCORR(J),ZERR(J))
      ELSE
        CALL COMPENSATED_SUM     (ZP(:,J),KNL(J),ZCORR(J),ZERR(J))
      ENDIF
    ENDIF

!--- gather partial sums and error bounds to all processors

    ZBUFFL(JL+1) = ZP(MAX(KNL(J),1_JPIM),J)

    IF (KNL(J)>0) THEN
      ZBUFFL(JL+2) = ZERR(J)
      ZBUFFL(JL+3) = ZCORR(J)
    ELSE
      ZBUFFL(JL+2) = 0.0_JPRB
      ZBUFFL(JL+3) = 0.0_JPRB
    ENDIF
    JL = JL + IBUFLEN
  ENDDO

  IF (INPROC>1) THEN
!-- could use MPL_ALLGATHER here, if it existed!

    IRECVCOUNTS(:) = SIZE(ZBUFFL)

    CALL MPL_ALLGATHERV (ZBUFFL,ZBUFFG,IRECVCOUNTS, &
                       & CDSTRING='ORDER_INDEP_GLOBAL_SUM2')

    DO JP=1,INPROC
      JL = 0
      DO J = 1,KNVEC
        ZPSUMS(JP,J) = ZBUFFG(JL+1+(JP-1)*SIZE(ZBUFFL))
        ZPERRS(JP,J) = ZBUFFG(JL+2+(JP-1)*SIZE(ZBUFFL))
        ZPCORS(JP,J) = ZBUFFG(JL+3+(JP-1)*SIZE(ZBUFFL))
        JL = JL + IBUFLEN
      ENDDO
    ENDDO
  ELSE
    JL = 0
    DO J = 1,KNVEC
      ZPSUMS(1,J) = ZBUFFL(JL+1)
      ZPERRS(1,J) = ZBUFFL(JL+2)
      ZPCORS(1,J) = ZBUFFL(JL+3)
    ENDDO
  ENDIF

!--- transform partial sums

  DO J = 1,KNVEC
    IF (LLDONE(J)) CYCLE

    CALL COMPENSATED_SUM (ZPSUMS(:,J),INPROC,ZCORR(J),ZERR(J))
    ZERR(J)  = ZERR(J)  + SUM(ZPERRS(:,J))
    ZCORR(J) = ZCORR(J) + SUM(ZPCORS(:,J))
    
!--- Calculate final result

    ZRES(J) = ZPSUMS(INPROC,J) + ZCORR(J)

!--- Calculate error bound. This is corollary 4.7 from Ogita et al. (2005)

    ZBETA(J) = ZERR(J)*(REAL(2*ING(J),JPRB)*EPSILON(ZRES(J))) &
       & /(1.0_JPRB - REAL(2*ING(J),JPRB)*EPSILON(ZRES(J)))

    ZERR(J) = EPSILON(ZRES(J))*ABS(ZRES(J)) &
     & +(ZBETA(J) + ( 2.0_JPRB*EPSILON(ZRES(J))*EPSILON(ZRES(J))*ABS(ZRES(J)) &
     &            +3.0_JPRB*TINY(ZRES(J))))

!--- update the last element of the local array

    ZP(MAX(KNL(J),1_JPIM),J) = ZPSUMS(MPL_MYRANK(),J)
  ENDDO

!--- exit if the global error is small enough

  LLDONE(:) = (ZERR(:)<4.0_JPRB*SPACING(ZRES(:))) .OR. LLDONE(:)

  IF (ALL(LLDONE(:))) EXIT K_LOOP

!--- Take appropriate action if ZRES cannot be sufficiently refined.

  DO J = 1,KNVEC
    IF (ZERR(J) >= ZOLDERR(J) .AND. .NOT. LLDONE(J)) THEN
      INMSG=INMSG+1

      IF (INMSG<=100) THEN
        CALL MPL_MESSAGE ( &
          & CDMESSAGE= 'ORDER_INDEP_GLOBAL_SUM2: FALIED TO REFINE SUM', &
          & CDSTRING='ORDER_INDEP_GLOBAL_SUM2')
        CALL MPL_MESSAGE ( &
          & CDMESSAGE='WARNING: POSSIBLITY OF NON-REPRODUCIBLE RESULTS',&
          & CDSTRING='ORDER_INDEP_GLOBAL_SUM2')
      ENDIF
      
      IF (INMSG==100) THEN
        CALL MPL_MESSAGE ( &
          & CDMESSAGE='ORDER_INDEP_GLOBAL_SUM2: INMSG>100. OUTPUT SUPPRESSED',&
          & CDSTRING='ORDER_INDEP_GLOBAL_SUM2')
      ENDIF

      IF (LLABORT) THEN
        CALL MPL_MESSAGE (CDMESSAGE= &
          & 'ABORT BECAUSE LD_ABORT_IFNOT_REPROD WAS SET', &
          & CDSTRING='ORDER_INDEP_GLOBAL_SUM2',LDABORT=.TRUE.)
      ENDIF
    ENDIF

    ZOLDERR(J) = ZERR(J)
  ENDDO

ENDDO K_LOOP

!--- At this stage, we have guaranteed that ZRES less than 4*EPS
!--- away from the exact sum. There are only four floating point
!--- numbers in this range. So, if we find the nearest number that
!--- has its last three bits zero, then we have a reproducible result.

DO J=1,KNVEC
  POUT(J) = ROUND (ZRES(J))
ENDDO

DEALLOCATE (IRECVCOUNTS)
DEALLOCATE (ZPCORS)
DEALLOCATE (ZPERRS)
DEALLOCATE (ZPSUMS)
DEALLOCATE (ZBUFFG)
DEALLOCATE (ZBUFFL)
DEALLOCATE (ZP)

IF (LHOOK) CALL DR_HOOK ('ORDER_INDEPENDENT_SUMMATION_MOD:ORDER_INDEP_GLOBAL_SUM2', &
                      &  1,ZHOOK_HANDLE)

END SUBROUTINE ORDER_INDEP_GLOBAL_SUM2

SUBROUTINE ORDER_INDEP_ALLREDUCE (PIN,POUT,LD_ABORT_IFNOT_REPROD,LD_OPENMP)

!-----------------------------------------------------------------
!
!  Returns in POUT an accurate global sum of the elements of PIN across tasks. 
!  This has a similar functionality to MPL_allreduce where an array is supplied
!  and we want the individual elements of the array to be summed across tasks.
!  This is different to ORDER_INDEP_GLOBAL_SUM where PIN is considered part of
!  a global array.
!
! Arguments:
! ----------
!
! Required:
!
!   PIN                  - INTENT(IN) - input array to be summed.
!   POUT                 - INTENT(OUT) - output array of same size
!
!
! Optional:
!
!   LD_ABORT_IFNOT_REPROD - INTENT(IN) - Defines behaviour in case
!                                        a reproducible result cannot
!                                        be guaranteed.
!
!   LD_OPENMP            - INTENT(IN) - Use OpenMP parallelization.
!              
!
!  Author: George Mozdzynski ECMWF June 2009
!
!-----------------------------------------------------------------

IMPLICIT NONE
  
REAL(KIND=JPRB) :: ORDER_INDEP_GLOBAL_SUM
  
REAL(KIND=JPRB),             INTENT(IN) :: PIN(:)
REAL(KIND=JPRB),             INTENT(OUT):: POUT(:)
LOGICAL,           OPTIONAL, INTENT(IN) :: LD_ABORT_IFNOT_REPROD, LD_OPENMP

INTEGER(KIND=JPIM) :: INPROC,MYPROC,IN,ITAG,I,J,IR
INTEGER(KIND=JPIM), ALLOCATABLE :: ICOUNT(:),IND(:),IREQ(:)
REAL(KIND=JPRB), ALLOCATABLE :: ZBUFF(:),ZIN(:),ZOUT(:)

LOGICAL :: LLABORT, LL_OPENMP
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK ('ORDER_INDEPENDENT_SUMMATION_MOD:ORDER_INDEP_ALLREDUCE', &
                      &  0,ZHOOK_HANDLE) 

INPROC = MPL_NPROC()

IF (PRESENT(LD_ABORT_IFNOT_REPROD)) THEN
  LLABORT = LD_ABORT_IFNOT_REPROD
ELSE
  LLABORT = .TRUE.
ENDIF

IF (PRESENT(LD_OPENMP)) THEN
  LL_OPENMP = LD_OPENMP
ELSE
  LL_OPENMP = .TRUE.
ENDIF

IF( SIZE(PIN) /= SIZE(POUT) )THEN
  CALL MPL_MESSAGE (CDMESSAGE='SIZE(PIN) /= SIZE(POUT)', &
                  & CDSTRING='ORDER_INDEP_ALLREDUCE',LDABORT=.TRUE.)
ENDIF

IN = SIZE(PIN)

IF (INPROC==1) THEN
  POUT(:)=PIN(:)
ELSE
  ITAG=1234
  ALLOCATE(ICOUNT(INPROC))
  ICOUNT(:) = 0

! Determine distribution of input array over tasks

  DO J=1,IN
    I=MOD(J-1,INPROC)+1
    ICOUNT(I)=ICOUNT(I)+1
  ENDDO
  ALLOCATE(IND(INPROC))
  IND(:)=0
  IND(1)=1
  DO J=2,INPROC
    IND(J)=IND(J-1)+ICOUNT(J-1)
  ENDDO
  MYPROC = MPL_MYRANK()
  ALLOCATE(ZBUFF(ICOUNT(MYPROC)*INPROC))
  ALLOCATE(IREQ(2*INPROC))
  ALLOCATE(ZIN(INPROC))

! Distribute input array over tasks

  IR=0
  IF(ICOUNT(MYPROC) /= 0)THEN
    DO J=1,INPROC
      IR=IR+1
      CALL MPL_RECV (ZBUFF((J-1)*ICOUNT(MYPROC)+1:J*ICOUNT(MYPROC)),&
                    &KSOURCE=J,&
                    &KTAG=ITAG,&
                    &KMP_TYPE=JP_NON_BLOCKING_STANDARD,&
                    &KREQUEST=IREQ(IR),&
                    &CDSTRING='ORDER_INDEP_ALLREDUCE')
    ENDDO
  ENDIF
  DO J=1,INPROC
    IF(ICOUNT(J) /= 0)THEN
      IR=IR+1
      CALL MPL_SEND(PIN(IND(J):IND(J)+ICOUNT(J)-1),&
                   &KDEST=J,&
                   &KTAG=ITAG,&
                   &KMP_TYPE=JP_NON_BLOCKING_STANDARD,&
                   &KREQUEST=IREQ(IR),&
                   &CDSTRING='ORDER_INDEP_ALLREDUCE')
    ENDIF
  ENDDO
  IF(IR > 0)THEN
    CALL MPL_WAIT(KREQUEST=IREQ(1:IR),&
                 &CDSTRING='ORDER_INDEP_ALLREDUCE')
  ENDIF

! Perform local order independent sums for myproc's part of input array

  ALLOCATE(ZOUT(IN))
  DO J=1,ICOUNT(MYPROC)
    DO I=1,INPROC
      ZIN(I)=ZBUFF((I-1)*ICOUNT(MYPROC)+J)
    ENDDO
    ZOUT(J)=ORDER_INDEP_LOCAL_SUM(ZIN,LLABORT,LL_OPENMP)
  ENDDO

! Gather results of order independent sums over tasks

  CALL MPL_ALLGATHERV (ZOUT(1:ICOUNT(MYPROC)),POUT,ICOUNT, &
                     & CDSTRING='ORDER_INDEP_ALLREDUCE')

  DEALLOCATE (ICOUNT)
  DEALLOCATE (IND)
  DEALLOCATE (IREQ)
  DEALLOCATE (ZBUFF)
  DEALLOCATE (ZIN)
  DEALLOCATE (ZOUT)

ENDIF

IF (LHOOK) CALL DR_HOOK ('ORDER_INDEPENDENT_SUMMATION_MOD:ORDER_INDEP_ALLREDUCE', &
                      &  1,ZHOOK_HANDLE)

END SUBROUTINE ORDER_INDEP_ALLREDUCE

FUNCTION ORDER_INDEP_DOT_PRODUCT (P1,P2,PW,KNG,LD_ABORT_IFNOT_REPROD, &
                                & LD_OPENMP)

!-----------------------------------------------------------------
!
!  Returns an accurate global sum of the elements of P1*P2, or
!  P1*P2*PW. The result is identical to the result that would be
!  obtained by the following:
!
!  IF (PRESENT(PW)) THEN
!    PTEMP(:) = P1(:)*P2(:)*PW(:)
!  ELSE
!    PTEMP(:) = P1(:)*P2(:)
!  ENDIF
!  CALL ORDER_INDEP_GLOBAL_SUM (PTEMP,KNG,LD_ABORT_IFNOT_REPROD, &
!                              & LD_OPENMP)
!
!   This routine is provided only because the above is not very
!   cache-friendly.
!
!  NB: P1, P2 and PW are unmodified on return
!
! Arguments:
! ----------
!
! Required:
!
!   P1,P2                - INTENT(IN) - Arrays whose inner product is
!                                       to be calculated
!
!
! Optional:
!
!   PW                   - INTENT(IN) - Weight array defining the
!                                       metric for the inner product.
!
!   KNG                  - INTENT(IN) - Global length of array.
!                 
!   LD_ABORT_IFNOT_REPRO - INTENT(IN) - Defines behaviour in case
!                                       a reproducible result cannot
!                                       be guaranteed.
!
!   LD_OPENMP            - INTENT(IN) - Use OpenMP parallelization.
!              
!
! Algorithm:
! ----------
!
!  The algorithm is based on Ogita et al. (2005) SIAM J. Sci. Computing,
!  Vol.26, No.6, pp1955-1988. This is based in turn on an algorithm
!  by Knuth (1969, seminumerical algorithms).
!
!  This version adds a second layer of parallelism on top of that
!  provided by COMPENSATED_SUM_OMP. It iterates the compensated 
!  summation until the result is guaranteed to be within 4*eps
!  of the true sum. It then rounds the result to the nearest
!  floating-point number whose last three bits are zero, thereby
!  guaranteeing an order-independent result.
!
!  Author: Mike Fisher ECMWF 2006/02/08
!
!-----------------------------------------------------------------

IMPLICIT NONE

REAL(KIND=JPRB) :: ORDER_INDEP_DOT_PRODUCT

REAL(KIND=JPRB),             INTENT(IN) :: P1(:), P2(:)
REAL(KIND=JPRB),   OPTIONAL, INTENT(IN) :: PW(:)
INTEGER(KIND=JPIM),OPTIONAL, INTENT(IN) :: KNG
LOGICAL,           OPTIONAL, INTENT(IN) :: LD_ABORT_IFNOT_REPROD, LD_OPENMP

INTEGER(KIND=JPIM) :: J,IN,ING,INPROC
REAL(KIND=JPRB) :: ZCORR,ZERR,ZOLDERR,ZBUFFL(3),ZBETA,ZRES
REAL(KIND=JPRB), ALLOCATABLE :: ZPSUMS(:),ZPERRS(:),ZPCORS(:), &
                              & ZBUFFG(:),ZP(:)
INTEGER(KIND=JPIM), ALLOCATABLE :: IRECVCOUNTS(:)
LOGICAL :: LLABORT, LL_OPENMP, LL_FIRST_ITER
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

INTEGER(KIND=JPIM), SAVE :: INMSG=0

IF (LHOOK) CALL DR_HOOK ( &
              &'ORDER_INDEPENDENT_SUMMATION_MOD:ORDER_INDEP_DOT_PRODUCT', &
              &  0,ZHOOK_HANDLE)

INPROC = MPL_NPROC()

IF (PRESENT(LD_ABORT_IFNOT_REPROD)) THEN
  LLABORT = LD_ABORT_IFNOT_REPROD
ELSE
  LLABORT = .TRUE.
ENDIF

IF (PRESENT(LD_OPENMP)) THEN
  LL_OPENMP = LD_OPENMP
ELSE
  LL_OPENMP = .TRUE.
ENDIF

IN = SIZE(P1)

IF (SIZE(P2)/=IN) THEN
  CALL MPL_MESSAGE (CDMESSAGE='SIZE(P2)/=SIZE(P1)', &
                  & CDSTRING='ORDER_INDEP_DOT_PRODUCT',LDABORT=.TRUE.)
ENDIF

IF (PRESENT(PW)) THEN
  IF (SIZE(PW)/=IN) THEN
    CALL MPL_MESSAGE (CDMESSAGE='SIZE(PW)/=SIZE(P1)', &
                    & CDSTRING='ORDER_INDEP_DOT_PRODUCT',LDABORT=.TRUE.)
  ENDIF
ENDIF

!--- global length of vector (needed for error bound calculation)

IF (.NOT.PRESENT(KNG)) THEN
  ING = IN
  IF (INPROC>1) THEN
    CALL MPL_ALLREDUCE (ING,'SUM',CDSTRING='ORDER_INDEP_DOT_PRODUCT')
  ENDIF
ELSE
  ING = KNG
  IF (KNG<IN) THEN
    CALL MPL_MESSAGE (CDMESSAGE='Specified KNG < SIZE(PIN)', &
                    & CDSTRING='ORDER_INDEP_DOT_PRODUCT',LDABORT=.TRUE.)
  ENDIF
ENDIF

IF (REAL(2*ING,JPRB)*EPSILON(ZRES) >= 1.0) THEN
  CALL MPL_MESSAGE (CDMESSAGE='n is too large to guarantee error bounds', &
                  & CDSTRING='ORDER_INDEP_DOT_PRODUCT',LDABORT=.TRUE.)
ENDIF

ALLOCATE (ZP(MAX(IN,1_JPIM)))
ALLOCATE (ZBUFFG(INPROC*SIZE(ZBUFFL)))
ALLOCATE (ZPSUMS(INPROC))
ALLOCATE (ZPERRS(INPROC))
ALLOCATE (ZPCORS(INPROC))
ALLOCATE (IRECVCOUNTS(INPROC))

ZOLDERR = HUGE(ZERR)

!--- Copy the input array. This avoids some tricky indexing, at the
!--- expense of some inefficency.

IF (IN==0) THEN
  ZP(1) = 0.0_JPRB
ENDIF

LL_FIRST_ITER = .TRUE.

K_LOOP: DO

!--- transform local arrays

  IF (IN>0) THEN
    IF (LL_FIRST_ITER) THEN
      IF (PRESENT(PW)) THEN
        IF (LL_OPENMP) THEN
          CALL COMPENSATED_DOT_PRODUCT_OMP (P1,P2,PW,ZP,IN,ZCORR,ZERR)
        ELSE
          CALL COMPENSATED_DOT_PRODUCT     (P1,P2,PW,ZP,IN,ZCORR,ZERR)
        ENDIF
      ELSE
        IF (LL_OPENMP) THEN
          CALL COMPENSATED_DOT_PRODUCT_OMP (P1=P1,P2=P2,POUT=ZP, &
                                          & KN=IN,PCORR=ZCORR,PERR=ZERR)
        ELSE
          CALL COMPENSATED_DOT_PRODUCT     (P1=P1,P2=P2,POUT=ZP, &
                                          & KN=IN,PCORR=ZCORR,PERR=ZERR)
        ENDIF
      ENDIF
    ELSE
      IF (LL_OPENMP) THEN
        CALL COMPENSATED_SUM_OMP (ZP,IN,ZCORR,ZERR)
      ELSE
        CALL COMPENSATED_SUM     (ZP,IN,ZCORR,ZERR)
      ENDIF
    ENDIF
  ENDIF

!--- gather partial sums and error bounds to all processors

  ZBUFFL(1) = ZP(MAX(IN,1_JPIM))

  IF (IN>0) THEN
    ZBUFFL(2) = ZERR
    ZBUFFL(3) = ZCORR
  ELSE
    ZBUFFL(2) = 0.0_JPRB
    ZBUFFL(3) = 0.0_JPRB
  ENDIF

  IF (INPROC>1) THEN
!-- could use MPL_ALLGATHER here, if it existed!

    IRECVCOUNTS(:) = SIZE(ZBUFFL)

    CALL MPL_ALLGATHERV (ZBUFFL,ZBUFFG,IRECVCOUNTS, &
                       & CDSTRING='ORDER_INDEP_DOT_PRODUCT')

    DO J=1,INPROC
      ZPSUMS(J) = ZBUFFG(1+(J-1)*SIZE(ZBUFFL))
      ZPERRS(J) = ZBUFFG(2+(J-1)*SIZE(ZBUFFL))
      ZPCORS(J) = ZBUFFG(3+(J-1)*SIZE(ZBUFFL))
    ENDDO
  ELSE
    ZPSUMS(1) = ZBUFFL(1)
    ZPERRS(1) = ZBUFFL(2)
    ZPCORS(1) = ZBUFFL(3)
  ENDIF

!--- transform partial sums

  CALL COMPENSATED_SUM (ZPSUMS,INPROC,ZCORR,ZERR)
  ZERR  = ZERR  + SUM(ZPERRS)
  ZCORR = ZCORR + SUM(ZPCORS)

!--- Calculate final result

  ZRES = ZPSUMS(INPROC) + ZCORR

!--- Calculate error bound. This is corollary 4.7 from Ogita et al. (2005)

  ZBETA = ZERR*(REAL(2*ING,JPRB)*EPSILON(ZRES)) &
       & /(1.0_JPRB - REAL(2*ING,JPRB)*EPSILON(ZRES))

  ZERR = EPSILON(ZRES)*ABS(ZRES) &
     & +(ZBETA + ( 2.0_JPRB*EPSILON(ZRES)*EPSILON(ZRES)*ABS(ZRES) &
     &            +3.0_JPRB*TINY(ZRES)))

!--- update the last element of the local array

  ZP(MAX(IN,1_JPIM)) = ZPSUMS(MPL_MYRANK())

!--- exit if the global error is small enough

  IF (ZERR<4.0_JPRB*SPACING(ZRES)) EXIT K_LOOP

!--- Take appropriate action if ZRES cannot be sufficiently refined.

  IF (ZERR >= ZOLDERR) THEN
    INMSG=INMSG+1

    IF (INMSG<=100) THEN
      CALL MPL_MESSAGE ( &
          & CDMESSAGE= 'ORDER_INDEP_DOT_PRODUCT: FALIED TO REFINE SUM', &
          & CDSTRING='ORDER_INDEP_DOT_PRODUCT')
      CALL MPL_MESSAGE ( &
          & CDMESSAGE='WARNING: POSSIBLITY OF NON-REPRODUCIBLE RESULTS',&
          & CDSTRING='ORDER_INDEP_DOT_PRODUCT')
    ENDIF

    IF (INMSG==100) THEN
      CALL MPL_MESSAGE ( &
       & CDMESSAGE='ORDER_INDEP_DOT_PRODUCT: INMSG>100. OUTPUT SUPPRESSED',&
       & CDSTRING='ORDER_INDEP_DOT_PRODUCT')
    ENDIF

    IF (LLABORT) THEN
      CALL MPL_MESSAGE (CDMESSAGE= &
                      & 'ABORT BECAUSE LD_ABORT_IFNOT_REPROD WAS SET', &
                      & CDSTRING='ORDER_INDEP_DOT_PRODUCT',LDABORT=.TRUE.)
    ENDIF
  ENDIF

  ZOLDERR = ZERR

  LL_FIRST_ITER = .FALSE.
ENDDO K_LOOP

!--- At this stage, we have guaranteed that ZRES less than 4*EPS
!--- away from the exact sum. There are only four floating point
!--- numbers in this range. So, if we find the nearest number that
!--- has its last three bits zero, then we have a reproducible result.

ORDER_INDEP_DOT_PRODUCT = ROUND (ZRES)

DEALLOCATE (IRECVCOUNTS)
DEALLOCATE (ZPCORS)
DEALLOCATE (ZPERRS)
DEALLOCATE (ZPSUMS)
DEALLOCATE (ZBUFFG)
DEALLOCATE (ZP)

IF (LHOOK) CALL DR_HOOK ( &
         &'ORDER_INDEPENDENT_SUMMATION_MOD:ORDER_INDEP_DOT_PRODUCT', &
         &  1,ZHOOK_HANDLE)

END FUNCTION ORDER_INDEP_DOT_PRODUCT

FUNCTION ROUND (PRES)

!-----------------------------------------------------------------
!
!  Returns the value of PRES rounded to the nearest floating-point
!  number that has its last three bits zero

!  The code to do this in Fortran is not nice, because Fortran
!  does not proved access to the binary representation for REALs.
!  Perhaps we should code it in c?

!  This works on big-endian and little-endian machines.
 
!  Author: Mike Fisher ECMWF 2006/02/08
!
!-----------------------------------------------------------------

IMPLICIT NONE

REAL(KIND=JPRB), INTENT(IN) :: PRES
REAL(KIND=JPRB) :: ROUND

INTEGER(KIND=JPIM) :: II(2),IEQUIV(8),INTS_PER_REAL,J,I_LOW_WORD
REAL(KIND=JPRB)    :: ZZ(2),ZUP,ZDOWN


II(:)=1
ZZ(:)=1.0_JPRB
INTS_PER_REAL=STORAGE_SIZE(ZZ)/STORAGE_SIZE(II)

IF (INTS_PER_REAL>SIZE(IEQUIV)) THEN
  CALL MPL_MESSAGE (CDMESSAGE='INTS_PER_REAL>SIZE(IEQUIV)', &
                    & CDSTRING='ORDER_INDEP_GLOBAL_SUM',LDABORT=.TRUE.)
ENDIF

!--- Test whether big-endian or little-endian

ZUP = -1.0_JPRB
IEQUIV(1:INTS_PER_REAL) = TRANSFER(ZUP,IEQUIV(1:INTS_PER_REAL))

IF (IEQUIV(1)==0) THEN
  I_LOW_WORD = 1                ! Little-endian
ELSE
  I_LOW_WORD = INTS_PER_REAL    ! Big-endian
ENDIF

!--- Find the nearest number with all 3 lowest-order bits zeroed

IEQUIV(1:INTS_PER_REAL) = TRANSFER(PRES,IEQUIV(1:INTS_PER_REAL))
ZUP    = PRES
ZDOWN  = PRES

IF (IBITS(IEQUIV(I_LOW_WORD),0,3)/=0) THEN
  DO J=1,4
    ZUP=NEAREST(ZUP,1.0_JPRB)
    IEQUIV(1:INTS_PER_REAL) = TRANSFER(ZUP,IEQUIV(1:INTS_PER_REAL))
    IF (IBITS(IEQUIV(I_LOW_WORD),0,3)==0) EXIT

    ZDOWN=NEAREST(ZDOWN,-1.0_JPRB)
    IEQUIV(1:INTS_PER_REAL) = TRANSFER(ZDOWN,IEQUIV(1:INTS_PER_REAL))
    IF (IBITS(IEQUIV(I_LOW_WORD),0,3)==0) EXIT
  ENDDO

  IF (IBITS(IEQUIV(I_LOW_WORD),0,3)/=0) THEN
    CALL MPL_MESSAGE (CDMESSAGE='THIS IS NOT POSSIBLE', &
                    & CDSTRING='ORDER_INDEP_GLOBAL_SUM',LDABORT=.TRUE.)
  ENDIF
ENDIF

ROUND = TRANSFER(IEQUIV(1:INTS_PER_REAL),PRES)

END FUNCTION ROUND

END MODULE ORDER_INDEPENDENT_SUMMATION_MOD
