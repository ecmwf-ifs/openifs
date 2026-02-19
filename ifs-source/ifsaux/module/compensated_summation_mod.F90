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

MODULE COMPENSATED_SUMMATION_MOD
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK

!**** COMPENSATED_SUMMATION_MOD
!     Purpose.
!     --------
!     Functions to perform compensated (i.e. accurate) summation.
!     These functions are used by ORDER_INDEPENDENT_SUMMATION_MOD

!**   Interface.
!     ----------

!      CALL COMPENSATED_SUM     (P,KN,PCORR,PERR)
!      CALL COMPENSATED_SUM_OMP (P,KN,PCORR,PERR)
! 
!        These routines transform the elements of the array P, such that:
!
!        1)  p(kn)         contains sum(p)
!
!        2)  p(1)...p(kn-1) contain the rounding errors that were made
!                           in calculating sum(p).
!        3)  The exact sum of the elements of p is unmodified.
!
!        On return, pcorr contains the sum of the rounding errors, perr
!        contains the sum of their absolute values.
!
!        After calling this routine, an accurate sum of the elements of p
!        can be calculated as res=p(n)+pcorr.
!
!        CALL COMPENSATED_SUM_OMP is an OpenMP-parallelized version of
!        CALL COMPENSATED_SUM.

     
!      CALL COMPENSATED_DOT_PRODUCT     (P1,P2,   POUT,PCORR,PERR)
!      CALL COMPENSATED_DOT_PRODUCT     (P1,P2,PW,POUT,PCORR,PERR)
!      CALL COMPENSATED_DOT_PRODUCT_OMP (P1,P2,   POUT,PCORR,PERR)
!      CALL COMPENSATED_DOT_PRODUCT_OMP (P1,P2,PW,POUT,PCORR,PERR)
!
!        These routines are variations on COMPENSATED_SUM that are
!        provided simply to reduce memory-access time.
 
!      CALL COMPENSATED_DOT_PRODUCT (P1,P2,POUT,PCORR,PERR)
!
!      is functionally equivalent to the following:
!
!          POUT(:) = P1(:)*P2(:)
!          CALL CALL COMPENSATED_SUM (POUT,KN,PCORR,PERR)

!      CALL COMPENSATED_DOT_PRODUCT_OMP (P1,P2,POUT,PCORR,PERR)
!
!      is functionally equivalent to the following:
!
!          POUT(:) = P1(:)*P2(:)
!          CALL CALL COMPENSATED_SUM_OMP (POUT,KN,PCORR,PERR)
 
!      CALL COMPENSATED_DOT_PRODUCT (P1,P2,PW,POUT,PCORR,PERR)
!
!      is functionlly equivalent to the following:
!
!          POUT(:) = P1(:)*P2(:)*PW(:)
!          CALL CALL COMPENSATED_SUM (POUT,KN,PCORR,PERR)

!      CALL COMPENSATED_DOT_PRODUCT_OMP (P1,P2,PW,POUT,PCORR,PERR)
!
!      is functionlly equivalent to the following:
!
!          POUT(:) = P1(:)*P2(:)*PW(:)
!          CALL CALL COMPENSATED_SUM_OMP (POUT,KN,PCORR,PERR)

!**   Algorithm
!     ---------
 
!  The algorithm is based on Ogita et al. (2005) SIAM J. Sci. Computing,
!  Vol.26, No.6, pp1955-1988. This is based in turn on an algorithm
!  by Knuth (1969, seminumerical algorithms).
!
!  The basic idea is that we can transform a pair of floating-point
!  numbers "a" and "b" into a new pair "x" and "y", such that:
!
!      x = add(a,b)  and x+y = a+b
!
!  where "add" denotes floating-point addition, and "+" denotes
!  exact, mathematical (i.e. infinite-precision) addition.
!
!  Applying this to an array, p, we can transform the array such that:
!
!      p(kn) := add(p(1),p(2),...,p(kn))
!  and p(1) + p(2) + ... + p(kn) is unchanged.


!     Author.
!     -------
!        Mike Fisher  ECMWF

!     Modifications.
!     --------------
!        Original: 2006-20-22

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

SAVE
PRIVATE
PUBLIC COMPENSATED_SUM, &
     & COMPENSATED_SUM_OMP, &
     & COMPENSATED_DOT_PRODUCT, &
     & COMPENSATED_DOT_PRODUCT_OMP

INTERFACE COMPENSATED_SUM
  MODULE PROCEDURE COMPENSATED_SUM
END INTERFACE COMPENSATED_SUM

INTERFACE COMPENSATED_SUM_OMP
  MODULE PROCEDURE COMPENSATED_SUM_OMP
END INTERFACE COMPENSATED_SUM_OMP

INTERFACE COMPENSATED_DOT_PRODUCT
  MODULE PROCEDURE COMPENSATED_DOT_PRODUCT
END INTERFACE COMPENSATED_DOT_PRODUCT

INTERFACE COMPENSATED_DOT_PRODUCT_OMP
  MODULE PROCEDURE COMPENSATED_DOT_PRODUCT_OMP
END INTERFACE COMPENSATED_DOT_PRODUCT_OMP

CONTAINS

SUBROUTINE COMPENSATED_SUM (P,KN,PCORR,PERR)
  IMPLICIT NONE

  INTEGER(KIND=JPIM), INTENT(IN)    :: KN
  REAL(KIND=JPRB),    INTENT(INOUT) :: P(KN)
  REAL(KIND=JPRB),    INTENT(OUT)   :: PCORR, PERR

  REAL(KIND=JPRB) :: ZX,ZZ,ZPSUM
  INTEGER(KIND=JPIM) :: J

  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  IF (LHOOK) CALL DR_HOOK('COMPENSATED_SUMMATION_MOD:COMPENSATED_SUM',0,ZHOOK_HANDLE)
  PCORR = 0.0
  PERR  = 0.0

  ZPSUM = P(1)
  DO J=2,KN
!--- It is vital that these 4 lines are not optimized in any way that
!--- changes the results.
    ZX     = P(J) + ZPSUM
    ZZ     = ZX   - P(J)
    P(J-1) = (P(J)-(ZX-ZZ)) + (ZPSUM-ZZ)
    ZPSUM  = ZX
!--- accumulate the correction and the error
    PCORR = PCORR + P(J-1)
    PERR  = PERR  + ABS(P(J-1))
  ENDDO
  P(KN) = ZPSUM
 
!-----------------------------------------------------------------
!  Vectorization
!  -------------
!
!  NB: As coded, the above loop may not run very well on a vector
!      computer. However, any loop ordering that preserves the exact
!      sum, and ends up with the floating-point sum in the last
!      element and rounding errors in the rest of the array,
!      should be OK. For example:
!
!      ILEN=KN
!      DO
!        IHALF=ILEN/2
!      !--- no vector dependency
!        DO J=1,IHALF
!--- It is vital that these 4 lines are not optimized in any way that
!--- changes the results.
!          ZX     = P(KN-IHALF+J) + P(KN-ILEN+J)
!          ZZ     = ZX   - P(KN-IHALF+J)
!          P(KN-ILEN+J)  = (P(KN-IHALF+J)-(ZX-ZZ)) + (P(KN-ILEN+J)-ZZ)
!          P(KN-IHALF+J) = ZX
!--- accumulate the correction and the error
!          PCORR = PCORR + P(KN-ILEN+J)
!          PERR  = PERR  + ABS(P(KN-ILEN+J))
!        ENDDO
!        ILEN=ILEN-IHALF
!        IF (ILEN<=1) EXIT
!      ENDDO
!-----------------------------------------------------------------
 
IF (LHOOK) CALL DR_HOOK('COMPENSATED_SUMMATION_MOD:COMPENSATED_SUM',1,ZHOOK_HANDLE)
END SUBROUTINE COMPENSATED_SUM

SUBROUTINE COMPENSATED_SUM_OMP (P,KN,PCORR,PERR)
  USE OML_MOD

  IMPLICIT NONE

  INTEGER(KIND=JPIM), INTENT(IN)    :: KN
  REAL(KIND=JPRB),    INTENT(INOUT) :: P(KN)
  REAL(KIND=JPRB),    INTENT(OUT)   :: PCORR, PERR

  REAL(KIND=JPRB), ALLOCATABLE :: ZERRS(:),ZCORS(:)
  REAL(KIND=JPRB) :: ZX,ZZ
  INTEGER(KIND=JPIM) :: J,JCHUNK,ILEN,INCHUNKS,IMINLEN,ILENCHUNK, &
                    & INTHREADS,I,ISTART,IEND

!--- IMINLEN is a tunable parameter. It represents the vector length
!--- below which there is too little work to make it worth spawning threads

  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  IF (LHOOK) CALL DR_HOOK('COMPENSATED_SUMMATION_MOD:COMPENSATED_SUM_OMP',0,ZHOOK_HANDLE)
  IMINLEN=1000  !-- this value is pure guesswork (Mike Fisher)

  INTHREADS = OML_MAX_THREADS()

  ILENCHUNK = MAX(IMINLEN,(KN+INTHREADS-1)/INTHREADS)
  INCHUNKS=1+(KN-1)/ILENCHUNK

  ALLOCATE(ZERRS(INCHUNKS))
  ALLOCATE(ZCORS(INCHUNKS))

!--- First, we split the array into chunks, and apply compensated_sum
!--- to each chunk independently.

!$OMP PARALLEL DO PRIVATE(ISTART,IEND), SCHEDULE(STATIC), IF(INCHUNKS>1)
  DO JCHUNK=1,INCHUNKS
    ISTART = 1+(JCHUNK-1)*ILENCHUNK
    IEND   = MIN(JCHUNK*ILENCHUNK,KN)
    CALL COMPENSATED_SUM (P(ISTART:IEND),1+IEND-ISTART, &
                      & ZCORS(JCHUNK), ZERRS(JCHUNK))
  ENDDO
!$OMP END PARALLEL DO

  PCORR = SUM(ZCORS)
  PERR  = SUM(ZERRS)

!--- The final element of each chunk contains a partial sum. We apply
!--- compensated summation to the vector of the final elements.

  DO JCHUNK=2,INCHUNKS
    I = MIN(JCHUNK*ILENCHUNK,KN)
    ILEN = I - (JCHUNK-1)*ILENCHUNK
!--- It is vital that these 4 lines are not optimized
    ZX        = P(I) + P(I-ILEN)
    ZZ        = ZX   - P(I)
    P(I-ILEN) = (P(I)-(ZX-ZZ)) + (P(I-ILEN)-ZZ)
    P(I)      = ZX
!--- accumulate the result and its error bound
    PCORR     = PCORR + P(I-ILEN)
    PERR      = PERR  + ABS(P(I-ILEN))
  ENDDO

  DEALLOCATE(ZERRS)
  DEALLOCATE(ZCORS)
IF (LHOOK) CALL DR_HOOK('COMPENSATED_SUMMATION_MOD:COMPENSATED_SUM_OMP',1,ZHOOK_HANDLE)
END SUBROUTINE COMPENSATED_SUM_OMP

SUBROUTINE COMPENSATED_DOT_PRODUCT (P1,P2,PW,POUT,KN,PCORR,PERR)
  IMPLICIT NONE

  INTEGER(KIND=JPIM), INTENT(IN)           :: KN
  REAL(KIND=JPRB),    INTENT(IN)           :: P1(KN), P2(KN)
  REAL(KIND=JPRB),    INTENT(IN), OPTIONAL :: PW(KN)
  REAL(KIND=JPRB),    INTENT(OUT)          :: POUT(KN)
  REAL(KIND=JPRB),    INTENT(OUT)          :: PCORR, PERR

  REAL(KIND=JPRB) :: ZX,ZZ,ZPJ,ZPSUM
  INTEGER(KIND=JPIM) :: J
!==================================================================
!INTEGER(KIND=JPIM) :: ILEN, IHALF
!==================================================================

  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  IF (LHOOK) CALL DR_HOOK('COMPENSATED_SUMMATION_MOD:COMPENSATED_DOT_PRODUCT',0,ZHOOK_HANDLE)
  PCORR = 0.0
  PERR  = 0.0

  IF (PRESENT(PW)) THEN
    ZPSUM = P1(1)*P2(1)*PW(1)
  ELSE
    ZPSUM = P1(1)*P2(1)
  ENDIF

  DO J=2,KN
    IF (PRESENT(PW)) THEN
      ZPJ = P1(J)*P2(J)*PW(J)
    ELSE
      ZPJ = P1(J)*P2(J)
    ENDIF
!--- It is vital that these 4 lines are not optimized in any way that
!--- changes the results.
    ZX     = ZPJ + ZPSUM
    ZZ     = ZX   - ZPJ
    POUT(J-1) = (ZPJ-(ZX-ZZ)) + (ZPSUM-ZZ)
    ZPSUM     = ZX
!--- accumulate the correction and the error
    PCORR = PCORR + POUT(J-1)
    PERR  = PERR  + ABS(POUT(J-1))
  ENDDO
  POUT(KN) = ZPSUM

!==================================================================
!== vectorized version
!      DO J=1,KN
!        IF (PRESENT(PW)) THEN
!          POUT(J) = P1(J)*P2(J)*PW(J)
!        ELSE
!          POUT(J) = P1(J)*P2(J)
!        ENDIF
!      ENDDO
!
!      ILEN=KN
!      DO
!        IHALF=ILEN/2
!      !--- no vector dependency
!        DO J=1,IHALF
!!--- It is vital that these 4 lines are not optimized in any way that
!!--- changes the results.
!          ZX     = POUT(KN-IHALF+J) + POUT(KN-ILEN+J)
!          ZZ     = ZX   - POUT(KN-IHALF+J)
!          POUT(KN-ILEN+J)  = (POUT(KN-IHALF+J)-(ZX-ZZ)) + (POUT(KN-ILEN+J)-ZZ)
!          POUT(KN-IHALF+J) = ZX
!!--- accumulate the correction and the error
!          PCORR = PCORR + POUT(KN-ILEN+J)
!          PERR  = PERR  + ABS(POUT(KN-ILEN+J))
!        ENDDO
!        ILEN=ILEN-IHALF
!        IF (ILEN<=1) EXIT
!      ENDDO
!==================================================================
IF (LHOOK) CALL DR_HOOK('COMPENSATED_SUMMATION_MOD:COMPENSATED_DOT_PRODUCT',1,ZHOOK_HANDLE)
END SUBROUTINE COMPENSATED_DOT_PRODUCT

SUBROUTINE COMPENSATED_DOT_PRODUCT_OMP (P1,P2,PW,POUT,KN,PCORR,PERR)
  USE OML_MOD

  IMPLICIT NONE

  INTEGER(KIND=JPIM), INTENT(IN)           :: KN
  REAL(KIND=JPRB),    INTENT(IN)           :: P1(KN), P2(KN)
  REAL(KIND=JPRB),    INTENT(IN), OPTIONAL :: PW(KN)
  REAL(KIND=JPRB),    INTENT(OUT)          :: POUT(KN)
  REAL(KIND=JPRB),    INTENT(OUT)          :: PCORR, PERR

  REAL(KIND=JPRB), ALLOCATABLE :: ZERRS(:),ZCORS(:)
  REAL(KIND=JPRB) :: ZX,ZZ
  INTEGER(KIND=JPIM) :: J,JCHUNK,ILEN,INCHUNKS,IMINLEN,ILENCHUNK, &
                    & INTHREADS,I,ISTART,IEND

!--- IMINLEN is a tunable parameter. It represents the vector length
!--- below which there is too little work to make it worth spawning threads

  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  IF (LHOOK) CALL DR_HOOK('COMPENSATED_SUMMATION_MOD:COMPENSATED_DOT_PRODUCT_OMP',0,ZHOOK_HANDLE)
  IMINLEN=1000  !-- this value is pure guesswork (Mike Fisher)

  INTHREADS = OML_MAX_THREADS()

  ILENCHUNK = MAX(IMINLEN,(KN+INTHREADS-1)/INTHREADS)
  INCHUNKS=1+(KN-1)/ILENCHUNK

  ALLOCATE(ZERRS(INCHUNKS))
  ALLOCATE(ZCORS(INCHUNKS))

!--- First, we split the array into chunks, and apply compensated_sum
!--- to each chunk independently.

  IF (PRESENT(PW)) THEN
!$OMP PARALLEL DO PRIVATE(ISTART,IEND), SCHEDULE(STATIC), IF(INCHUNKS>1)
    DO JCHUNK=1,INCHUNKS
      ISTART = 1+(JCHUNK-1)*ILENCHUNK
      IEND   = MIN(JCHUNK*ILENCHUNK,KN)
      CALL COMPENSATED_DOT_PRODUCT (P1(ISTART:IEND),P2(ISTART:IEND),   &
                        &           PW(ISTART:IEND),POUT(ISTART:IEND), &
                        &           1+IEND-ISTART, &
                        &           ZCORS(JCHUNK), ZERRS(JCHUNK))
    ENDDO
!$OMP END PARALLEL DO
  ELSE
!$OMP PARALLEL DO PRIVATE(ISTART,IEND), SCHEDULE(STATIC), IF(INCHUNKS>1)
    DO JCHUNK=1,INCHUNKS
      ISTART = 1+(JCHUNK-1)*ILENCHUNK
      IEND   = MIN(JCHUNK*ILENCHUNK,KN)
      CALL COMPENSATED_DOT_PRODUCT (P1=P1(ISTART:IEND),P2=P2(ISTART:IEND),  &
                        &           POUT=POUT(ISTART:IEND),           &
                        &           KN=1+IEND-ISTART, &
                        &           PCORR=ZCORS(JCHUNK), PERR=ZERRS(JCHUNK))
    ENDDO
!$OMP END PARALLEL DO
  ENDIF

  PCORR = SUM(ZCORS)
  PERR  = SUM(ZERRS)

!--- The final element of each chunk contains a partial sum. We apply
!--- compensated summation to the vector of the final elements.

  DO JCHUNK=2,INCHUNKS
    I = MIN(JCHUNK*ILENCHUNK,KN)
    ILEN = I - (JCHUNK-1)*ILENCHUNK
!--- It is vital that these 4 lines are not optimized
    ZX           = POUT(I) + POUT(I-ILEN)
    ZZ           = ZX   - POUT(I)
    POUT(I-ILEN) = (POUT(I)-(ZX-ZZ)) + (POUT(I-ILEN)-ZZ)
    POUT(I)      = ZX
!--- accumulate the result and its error bound
    PCORR     = PCORR + POUT(I-ILEN)
    PERR      = PERR  + ABS(POUT(I-ILEN))
  ENDDO

  DEALLOCATE(ZERRS)
  DEALLOCATE(ZCORS)
IF (LHOOK) CALL DR_HOOK('COMPENSATED_SUMMATION_MOD:COMPENSATED_DOT_PRODUCT_OMP',1,ZHOOK_HANDLE)
END SUBROUTINE COMPENSATED_DOT_PRODUCT_OMP

END MODULE COMPENSATED_SUMMATION_MOD
