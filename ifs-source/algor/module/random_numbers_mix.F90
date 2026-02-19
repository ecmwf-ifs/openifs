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
@PROCESS HOT(NOVECTOR) NOSTRICT
#endif
MODULE RANDOM_NUMBERS_MIX
USE PARKIND1    , ONLY : JPIM, JPRB
USE YOMHOOK     , ONLY : LHOOK, DR_HOOK, JPHOOK
USE BYTES_IO_MOD, ONLY : BYTES_IO_OPEN, BYTES_IO_CLOSE, BYTES_IO_READ, BYTES_IO_WRITE

!**** *RANDOM_NUMBERS_MIX*  - Portable Random Number Generator

!     Purpose.
!     --------
!           Generate machine-independent pseudo-random numbers

!**   Interface.
!     ----------
!        CALL initialize_random_numbers (kseed, yd_stream)
!        CALL uniform_distribution      (px   , yd_stream)
!        CALL gaussian_distribution     (px   , yd_stream)
!        CALL random_number_restartfile (fname, action)
!        CALL wr_rangen_state           (nunit)

!        Explicit arguments :
!        --------------------
!        kseed  (input)    : integer seed in the range [0,HUGE(kseed)]
!        yd_stream (optional) : the state of the random number generator
!        px     (output)   : array to receive random numbers in the range

!        In the case of uniform_distribution, px has values in the range [0.0,1.0)

!        Implicit arguments :
!        --------------------
!        None

!     Method.
!     -------
!        Based loosly on ZUFALL (Petersen, 1994).

!        The main difference between this generator and ZUFALL is that integer arithmetic
!        is used. This ensures portability to vector machines that implement different
!        real arithmetic. In particular, vector machines often implement non-IEEE
!        arithmetic for their vector pipes. This routine will give identical results for
!        any integer type with at least 32 bits.

!        The generator is a lagged-Fibonacci generator: x(i) = x(i-p) + x(i-q) mod 2**m.
!        Lagged-Fibonacci generators have very long repeat periods: (2**q -1) * 2**(m-1)
!        (i.e about 2.85E191 for q=607, m=30). They pass most tests for randomness.

!        p and q must be chosen carefully. Values from the following table are OK.
!        Larger values give better random numbers, but smaller values are more
!        cache-friendly.

!          q         p
!        9689      4187
!        4423      2098
!        2281      1029
!        1279       418
!         607       273
!         521       168
!         250       103
!         127        63
!          97        33
!          55        24
!          43        22
!          31        13
!          24        10

!        The initial q values of x are set using the binary shirt register method of
!        Burns and Pryor 1999.

!        Mascagni et al (1995) show how to choose different sets of initial values that
!        are guaranteed to be drawn from different maximal-length cycles. This requires
!        the initial values of x(1)...x(q) to be in "canonical form". Specifically,
!        x(1) must be zero and all but a particular one or two values of x must be
!        even. For q=607 and p=273, only one element (jpq-jps) must be odd.

!     Externals.
!     ----------
!        None

!     Reference.
!     ----------
!        Burns P.J. and Pryor D.V. 1999,
!                             Surface Radiative Transport at Large Scale via Monte Carlo.
!                             Annual Review of Heat Transfer, Vol 9.
!
!        Petersen W.P., 1994, Lagged Fibonacci Series Random Number Generator
!                             for the NEC SX-3. International Journal of High Speed Computing
!                             Vol. 6, No. 3, pp387-398.
!
!        Mascagni M., Cuccaro S.A., Pryor D.V., Robinson M.L., 1995,
!                             A Fast, High Quality and Reproducible Parallel Lagged-Fibonacci
!                             Pseudorandom Number Generator. Journal of Computational Physics
!                             Vol 119. pp211-219.

!     Author.
!     -------
!        Mike Fisher *ECMWF*

!     Modifications.
!     --------------
!        Original : 2002-09-25
!        Made parallel friendly: 2003-08-11 Robert Pincus
!        M Leutbecher: 2004-05-10 restart capability
!        M Fisher:     2005-03-30 replaced LCG initialization with shift register
!        F. Vana:      2015-12-17 single precision version 
!     ------------------------------------------------------------------

IMPLICIT NONE

SAVE

PRIVATE
PUBLIC RANDOMNUMBERSTREAM, &
     & INITIALIZE_RANDOM_NUMBERS, UNIFORM_DISTRIBUTION, GAUSSIAN_DISTRIBUTION ,&
     & RANDOM_NUMBER_RESTARTFILE, WR_RANGEN_STATE

INTEGER(KIND=JPIM), PARAMETER      :: JPP=273, JPQ=607, JPS=105
INTEGER(KIND=JPIM), PARAMETER      :: JPMM=30
INTEGER(KIND=JPIM), PARAMETER      :: JPM=2**JPMM
INTEGER(KIND=JPIM), PARAMETER      :: JPNUMSPLIT=(JPQ-2)/(JPP-1)
INTEGER(KIND=JPIM), PARAMETER      :: JPLENSPLIT=(JPQ-JPP+JPNUMSPLIT-1)/JPNUMSPLIT
INTEGER(KIND=JPIM), PARAMETER      :: INITVALUE = 12345678

TYPE RANDOMNUMBERSTREAM
  PRIVATE
  INTEGER(KIND=JPIM)                 :: IUSED
  INTEGER(KIND=JPIM)                 :: INITTEST ! Should initialize to zero, but can't in F90
  INTEGER(KIND=JPIM), DIMENSION(JPQ) :: IX 
  REAL(KIND=JPRB)                    :: ZRM
END TYPE RANDOMNUMBERSTREAM

CONTAINS
!-------------------------------------------------------------------------------
SUBROUTINE INITIALIZE_RANDOM_NUMBERS (KSEED, YD_STREAM) 
  !-------------------------------------------------------------------------------
  ! Initialize fibgen
  !-------------------------------------------------------------------------------
  INTEGER(KIND=JPIM)      , INTENT(IN   ) :: KSEED
  TYPE(RANDOMNUMBERSTREAM), INTENT(INOUT) :: YD_STREAM
  
  INTEGER, PARAMETER :: JPMASK=123459876
  INTEGER(KIND=JPIM), PARAMETER     :: JPWARMUP_SHFT=64, JPWARMUP_LFG=999
  INTEGER(KIND=JPIM)                :: IDUM,JK,JJ,JBIT
  REAL(KIND=JPRB), DIMENSION(JPWARMUP_LFG)   :: ZWARMUP

  !-------------------------------------------------------------------------------
  ! Initialize the buffer using a binary shift register (Burns and Pryor, 1999).
  ! The Galois representation is used for the shift register as it is more
  ! efficient than the Fibonacci representation. The magic numbers 31 and 87
  ! define the shift register primitive polynomial=(32,7,5,3,2,1,0).
  !
  ! To ensure that different seeds produce distinct initial buffer states in
  ! canonical form, bits 0...jpmm-2 of the initial seed (after XORing with jpmask
  ! and spinning up using the linear congruential generator) are used to construct
  ! x(2), and the remaining bits are used to construct x(jpq).
  !-------------------------------------------------------------------------------
  
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  IF (LHOOK) CALL DR_HOOK('RANDOM_NUMBERS_MIX:INITIALIZE_RANDOM_NUMBERS',0,ZHOOK_HANDLE)
  IDUM = ABS(IEOR(KSEED,JPMASK))
  IF (IDUM==0) IDUM=JPMASK

  DO JJ=1,JPWARMUP_SHFT
    IF (BTEST(IDUM,31)) THEN
      IDUM=IBSET(ISHFT(IEOR(IDUM,87),1),0)
    ELSE
      IDUM=IBCLR(ISHFT(IDUM,1),0)
    ENDIF
  ENDDO

  YD_STREAM%IX(1:JPQ-1)= 0
  YD_STREAM%IX(2)      = ISHFT(IBITS(IDUM,0,JPMM-1),1)
  YD_STREAM%IX(JPQ)    = IBITS(IDUM,JPMM-1,BIT_SIZE(IDUM)+1-JPMM)

  DO JBIT=1,JPMM-1
    DO JJ=3,JPQ-1
      IF (BTEST(IDUM,31)) THEN
        IDUM=IBSET(ISHFT(IEOR(IDUM,87),1),0)
        YD_STREAM%IX(JJ)=IBSET(YD_STREAM%IX(JJ),JBIT)
      ELSE
        IDUM=IBCLR(ISHFT(IDUM,1),0)
      ENDIF
    ENDDO
  ENDDO

  YD_STREAM%IX(JPQ-JPS) = IBSET(YD_STREAM%IX(JPQ-JPS),0)
  
  !-------------------------------------------------------------------------------
  ! Initialize some constants
  !-------------------------------------------------------------------------------
  
  YD_STREAM%IUSED=JPQ
  YD_STREAM%ZRM=1.0_JPRB/REAL(JPM,JPRB)
  
  !-------------------------------------------------------------------------------
  ! Check the calculation of jpnumsplit and jplensplit.
  !-------------------------------------------------------------------------------
  
  IF (JPP+JPNUMSPLIT*JPLENSPLIT < JPQ) THEN
    CALL ABOR1 ('initialize_random_numbers: upper limit of last loop < jpq')
  ENDIF
  
  IF (JPLENSPLIT >=JPP) THEN
    CALL ABOR1 ('initialize_random_numbers: loop length > jpp')
  ENDIF
  
  IF (JPNUMSPLIT>1) THEN
    IF ((JPQ-JPP+JPNUMSPLIT-2)/(JPNUMSPLIT-1) < JPP) THEN
      CALL ABOR1 ('initialize_random_numbers: jpnumsplit is bigger than necessary')
    ENDIF
  ENDIF

  !-------------------------------------------------------------------------------
  ! Set initTest to show that the stream is initialized.
  !-------------------------------------------------------------------------------

  YD_STREAM%INITTEST = INITVALUE
  
  !-------------------------------------------------------------------------------
  ! Warm up the generator.
  !-------------------------------------------------------------------------------

  CALL UNIFORM_DISTRIBUTION (ZWARMUP, YD_STREAM)

IF (LHOOK) CALL DR_HOOK('RANDOM_NUMBERS_MIX:INITIALIZE_RANDOM_NUMBERS',1,ZHOOK_HANDLE)
END SUBROUTINE INITIALIZE_RANDOM_NUMBERS

!@PROCESS HOT NOSTRICT
SUBROUTINE UNIFORM_DISTRIBUTION (PX,YD_STREAM)
  !--------------------------------------------------------------------------------
  ! Generate uniformly distributed random numbers in the range 0.0<= px < 1.0
  !--------------------------------------------------------------------------------
  INTEGER(KIND=JPIM), PARAMETER :: IVAR=Z"3FFFFFFF"
  TYPE(RANDOMNUMBERSTREAM)     , INTENT(INOUT) :: YD_STREAM
  REAL(KIND=JPRB), DIMENSION(:), INTENT(  OUT) :: PX

  INTEGER(KIND=JPIM)                :: JJ, JK, IN, IFILLED
  
  ! This test is a little dirty but Fortran 90 doesn't allow for the initialization
  !   of components of derived types. 
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
! DR_HOOK removed to reduce overhead
! IF (LHOOK) CALL DR_HOOK('RANDOM_NUMBERS_MIX:UNIFORM_DISTRIBUTION',0,ZHOOK_HANDLE)
  IF(YD_STREAM%INITTEST /= INITVALUE) &
    & CALL ABOR1 ('uniform_distribution called before initialize_random_numbers')
  
  !--------------------------------------------------------------------------------
  ! Copy numbers that were generated during the last call, but not used.
  !--------------------------------------------------------------------------------
  
  IN=SIZE(PX)
  IFILLED=0
  
  DO JJ=YD_STREAM%IUSED+1,MIN(JPQ,IN+YD_STREAM%IUSED)
    PX(JJ-YD_STREAM%IUSED) = YD_STREAM%IX(JJ)*YD_STREAM%ZRM
    IFILLED=IFILLED+1
  ENDDO
  
  YD_STREAM%IUSED=YD_STREAM%IUSED+IFILLED
  
  IF (IFILLED==IN)  THEN 
!   IF (LHOOK) CALL DR_HOOK('RANDOM_NUMBERS_MIX:UNIFORM_DISTRIBUTION',1,ZHOOK_HANDLE)
! DR_HOOK removed to reduce overhead
    RETURN
  ENDIF
  
  !--------------------------------------------------------------------------------
  ! Generate batches of jpq numbers until px has been filled
  !--------------------------------------------------------------------------------
  
  DO WHILE (IFILLED<IN)
  
  !--------------------------------------------------------------------------------
  ! Generate jpq numbers in vectorizable loops. The first loop is length jpp. The
  ! remaining jpq-jpp elements are calculated in loops of length shorter than jpp.
  !--------------------------------------------------------------------------------
  
  !OCL NOVREC
    DO JJ=1,JPP
!     yd_stream%ix(jj) = yd_stream%ix(jj) + yd_stream%ix(jj-jpp+jpq)
!     if (yd_stream%ix(jj)>=jpm) yd_stream%ix(jj) = yd_stream%ix(jj)-jpm
      YD_STREAM%IX(JJ) = IAND(IVAR,YD_STREAM%IX(JJ) + YD_STREAM%IX(JJ-JPP+JPQ))
    ENDDO
  
    DO JK=1,JPNUMSPLIT
  !OCL NOVREC
      DO JJ=1+JPP+(JK-1)*JPLENSPLIT,MIN(JPQ,JPP+JK*JPLENSPLIT)
!       yd_stream%ix(jj) = yd_stream%ix(jj) + yd_stream%ix(jj-jpp)
!       if (yd_stream%ix(jj)>=jpm) yd_stream%ix(jj) = yd_stream%ix(jj)-jpm
        YD_STREAM%IX(JJ) = IAND(IVAR,YD_STREAM%IX(JJ) + YD_STREAM%IX(JJ-JPP))
      ENDDO
    ENDDO
  
    YD_STREAM%IUSED = MIN(JPQ,IN-IFILLED)
    PX(IFILLED+1:IFILLED+YD_STREAM%IUSED) = YD_STREAM%IX(1:YD_STREAM%IUSED)*YD_STREAM%ZRM
    IFILLED = IFILLED+YD_STREAM%IUSED
  ENDDO
  
!IF (LHOOK) CALL DR_HOOK('RANDOM_NUMBERS_MIX:UNIFORM_DISTRIBUTION',1,ZHOOK_HANDLE)
! DR_HOOK removed to reduce overhead
END SUBROUTINE UNIFORM_DISTRIBUTION
!-------------------------------------------------------------------------------
SUBROUTINE GAUSSIAN_DISTRIBUTION (PX, YD_STREAM)
  TYPE(RANDOMNUMBERSTREAM), INTENT(INOUT) :: YD_STREAM
  REAL(KIND=JPRB)         , INTENT(  OUT) :: PX(:)
  !--------------------------------------------------------------------------------
  ! Generate normally-distributed random numbers using the Box-Muller method.
  !
  ! NB: this routine does not use buffering. This means that the following calls:
  !     call gaussian_distribution (zx(1:k))
  !     call gaussian_distribution (zx(k+1:n))
  ! will produce different numbers for elements k+1 onwards than the single call:
  !     call gaussian_distribution (zx(1:n))
  !--------------------------------------------------------------------------------
  
  INTEGER(KIND=JPIM) :: ILEN, J
  REAL(KIND=JPRB) :: ZFAC, ZTWOPI
  REAL(KIND=JPRB) :: ZX(SIZE(PX)+1)
  REAL(KIND=JPRB) :: ZEPSILON
  
  !--------------------------------------------------------------------------------
  ! Generate uniform random points in the range [0,1)
  !--------------------------------------------------------------------------------

    REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
    IF (LHOOK) CALL DR_HOOK('RANDOM_NUMBERS_MIX:GAUSSIAN_DISTRIBUTION',0,ZHOOK_HANDLE)
    
    ZEPSILON=100._JPRB*EPSILON(ZEPSILON)
    CALL UNIFORM_DISTRIBUTION (ZX, YD_STREAM)

  !--------------------------------------------------------------------------------
  ! Generate gaussian deviates using Box-Muller method
  !--------------------------------------------------------------------------------
  
  ZTWOPI = 8.0_JPRB*ATAN(1.0_JPRB)
  ILEN=SIZE(PX)
  
  DO J=1,ILEN-1,2
    ZFAC = SQRT(-2.0_JPRB*LOG(MAX(1.0_JPRB-ZX(J),ZEPSILON)))
    PX(J  ) = ZFAC*COS(ZTWOPI*ZX(J+1))
    PX(J+1) = ZFAC*SIN(ZTWOPI*ZX(J+1))
  ENDDO
  
  !--------------------------------------------------------------------------------
  ! Generate the last point if ilen is odd
  !--------------------------------------------------------------------------------
  
  IF (MOD(ILEN,2) /= 0) THEN
    ZFAC = SQRT(-2.0_JPRB*LOG(MAX(1.0_JPRB-ZX(ILEN),ZEPSILON)))
    PX(ILEN) = ZFAC*COS(ZTWOPI*ZX(ILEN+1))
  ENDIF
  
IF (LHOOK) CALL DR_HOOK('RANDOM_NUMBERS_MIX:GAUSSIAN_DISTRIBUTION',1,ZHOOK_HANDLE)
END SUBROUTINE GAUSSIAN_DISTRIBUTION
!-------------------------------------------------------------------------------
SUBROUTINE RANDOM_NUMBER_RESTARTFILE( CDFNAME, CDACTION,YD_STREAM )
!--------------------------------------------------------------------------------
!
! read (cdaction='r') or write (cdaction='w') restart file
! for random number generator
!
!--------------------------------------------------------------------------------
CHARACTER (LEN=*)       , INTENT(IN)    :: CDFNAME
CHARACTER (LEN=1  )     , INTENT(IN)    :: CDACTION
TYPE(RANDOMNUMBERSTREAM), INTENT(INOUT) :: YD_STREAM
  
INTEGER(KIND=JPIM) :: IUNIT, IRET, IBYTES_IN_JPIM

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('RANDOM_NUMBERS_MIX:RANDOM_NUMBER_RESTARTFILE',0,ZHOOK_HANDLE)
IBYTES_IN_JPIM= CEILING(REAL(BIT_SIZE(YD_STREAM%IUSED))/8.0_JPRB - TINY(1.0_JPRB))

IF (IBYTES_IN_JPIM /= 4) THEN
  CALL ABOR1('random_number_restartfile: number of bytes for JPIM is not 4 ')        
ENDIF

CALL BYTES_IO_OPEN(IUNIT, CDFNAME, CDACTION, IRET)
IF (IRET /= 0) THEN
  CALL ABOR1('random_number_restartfile: BYTES_IO_OPEN FAILED opening '//CDFNAME)    
ENDIF


IF (CDACTION=='r' .OR. CDACTION=='R') THEN
  CALL BYTES_IO_READ(IUNIT, YD_STREAM%IX,    IBYTES_IN_JPIM*JPQ, IRET)
  IF (IRET < 0) THEN
    CALL ABOR1('random_number_restartfile: BYTES_IO_READ could not read ix from '//CDFNAME)    
  ENDIF
  CALL BYTES_IO_READ(IUNIT, YD_STREAM%IUSED, IBYTES_IN_JPIM    , IRET)
  IF (IRET < 0) THEN
    CALL ABOR1('random_number_restartfile: BYTES_IO_READ could not read iused from '//CDFNAME)    
  ENDIF

!  l_initialized = .TRUE.

  YD_STREAM%INITTEST = INITVALUE
  YD_STREAM%ZRM=1.0_JPRB/REAL(JPM,JPRB)
ELSEIF(CDACTION=='w' .OR. CDACTION=='W') THEN
  CALL BYTES_IO_WRITE(IUNIT, YD_STREAM%IX, IBYTES_IN_JPIM*JPQ, IRET)
  IF (IRET < 0) THEN
    CALL ABOR1('random_number_restartfile: BYTES_IO_WRITE could not write ix on '//CDFNAME)    
  ENDIF
  CALL BYTES_IO_WRITE(IUNIT, YD_STREAM%IUSED, IBYTES_IN_JPIM , IRET)
  IF (IRET < 0) THEN
    CALL ABOR1('random_number_restartfile: BYTES_IO_WRITE could not write iused on '//CDFNAME)    
  ENDIF

ELSE
  CALL ABOR1 ('random_number_restartfile: cdaction = '//CDACTION//' is undefined.')
ENDIF

CALL BYTES_IO_CLOSE(IUNIT, IRET)
IF (IRET /= 0) THEN
  CALL ABOR1('random_number_restartfile: BYTES_IO_CLOSE FAILED closing '//CDFNAME)    
ENDIF

IF (LHOOK) CALL DR_HOOK('RANDOM_NUMBERS_MIX:RANDOM_NUMBER_RESTARTFILE',1,ZHOOK_HANDLE)
END SUBROUTINE RANDOM_NUMBER_RESTARTFILE


SUBROUTINE WR_RANGEN_STATE( KUNIT, YD_STREAM )
!--------------------------------------------------------------------------------
! write state of random number generator to unit kunit
!--------------------------------------------------------------------------------
INTEGER(KIND=JPIM)      , INTENT(IN) :: KUNIT
TYPE(RANDOMNUMBERSTREAM), INTENT(IN) :: YD_STREAM

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('RANDOM_NUMBERS_MIX:WR_RANGEN_STATE',0,ZHOOK_HANDLE)
WRITE( KUNIT, * ) 'module random_numbers_mix, generator state is'
WRITE( KUNIT, '(8I10)') YD_STREAM%IX
WRITE( KUNIT, '(I10)')  YD_STREAM%IUSED

IF (LHOOK) CALL DR_HOOK('RANDOM_NUMBERS_MIX:WR_RANGEN_STATE',1,ZHOOK_HANDLE)
END SUBROUTINE WR_RANGEN_STATE

END MODULE RANDOM_NUMBERS_MIX
