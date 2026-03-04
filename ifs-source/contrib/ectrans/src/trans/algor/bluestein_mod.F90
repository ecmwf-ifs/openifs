! (C) Copyright 2015- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE BLUESTEIN_MOD

! Implementation of the Bluestein FFT algorithm as described in a paper titled
! "Bluestein's FFT for Arbitrary N on the Hypercube", Paul N. Swarztrauber et al.,
! Parallel Computing, 17 (1991), pp. 607-617.
!
! George Mozdzynski and Nils Wedi, June 2015
!
! The naming convention follows the algorithm description in the above paper.
!
USE PARKIND1, ONLY : JPIM, JPRB

IMPLICIT NONE

PRIVATE
PUBLIC BLUESTEIN_FFT, BLUESTEIN_INIT, BLUESTEIN_TERM, FFTB_TYPE

TYPE FFTB_PLAN
  INTEGER(KIND=JPIM)             :: NSIZE          ! latitude length security check
  REAL(KIND=JPRB),ALLOCATABLE :: HS(:,:,:)
  REAL(KIND=JPRB),ALLOCATABLE :: H2xT(:,:,:)
END TYPE FFTB_PLAN

TYPE FFTB_TYPE
  INTEGER(KIND=JPIM) :: NDLON ! maximum number of points on a latitude
  REAL(KIND=JPRB)   ,ALLOCATABLE :: TRIGS(:,:) ! list of trigonometric function values (PO2)
  INTEGER(KIND=JPIM),ALLOCATABLE :: NFAX(:,:)  ! list of factors of truncation (PO2)
  INTEGER(KIND=JPIM)             :: NLAT_COUNT ! number of lats requiring bluestein FFT
  INTEGER(KIND=JPIM),ALLOCATABLE :: NLATS(:)   ! the latitude lengths of these latitudes
  TYPE(FFTB_PLAN),ALLOCATABLE :: FFTB(:)
END TYPE FFTB_TYPE

CONTAINS
!-----------------------------------------------------------------------------
SUBROUTINE BLUESTEIN_FFT(TB,N,KSIGN,KLOT,PDAT)
! N   : FFT LENGTH  
! KSIGN   : FFT DIRECTION
!           -1  DIRECT  (R2C)
!            1  INVERSE (C2R)
IMPLICIT NONE
TYPE(FFTB_TYPE),INTENT(INOUT) :: TB
INTEGER,INTENT(IN) :: N,KSIGN,KLOT
REAL(KIND=JPRB),INTENT(INOUT)  :: PDAT (:,:)
REAL(KIND=JPRB),ALLOCATABLE :: ZDATAR(:,:), ZDATAI(:,:),ZY(:,:)
REAL(KIND=JPRB) :: ZR(KLOT),ZI(KLOT),ZX0(KLOT)
REAL(KIND=JPRB) :: ZWR,ZWI
INTEGER(KIND=JPIM) :: I,K,M,JLOT,NN,II,IR,IPO2
INTEGER(KIND=JPIM) :: IJUMP,ILOT,IINC,ISIGN,IFFTSIGN

!WRITE(*,'("BLUESTEIN_FFT: N=",I6," KSIGN=",I2," KLOT=",I4)')&
!  & N,KSIGN,KLOT

IF( KSIGN/=-1 .AND. KSIGN/=1 )THEN
  CALL ABOR1('BLUESTEIN_FFT: INVALID KSIGN')
ENDIF

NN=N/2+1

IF( TB%FFTB(N)%NSIZE /= N )THEN
  WRITE(0,'("BLUESTEIN_FFT: UNEXPECTED PLAN LATITUDE LENGTH, N=",I6," TB%FFTB(N)%NSIZE=",I6)')&
   & N,TB%FFTB(N)%NSIZE
  CALL ABOR1('BLUESTEIN_FFT: UNEXPECTED PLAN LATITUDE LENGTH')
ENDIF

IF( KSIGN==-1 )THEN
  ISIGN=1
ELSE
  ISIGN=2
ENDIF

! input data preparation

ALLOCATE(ZDATAR(KLOT,0:2*NN-1))
ALLOCATE(ZDATAI(KLOT,0:2*NN-1))
ZDATAR(:,:)=0.0D0
ZDATAI(:,:)=0.0D0
IF( KSIGN==-1 )THEN
  DO K=0,N-1
    DO JLOT=1,KLOT
      ZDATAR(JLOT,K)=PDAT(JLOT,K+1)
    ENDDO
  ENDDO
ELSEIF( KSIGN==1 )THEN
  DO JLOT=1,KLOT
    DO K=0,NN-1
      ZDATAR(JLOT,K)=PDAT(JLOT,K*2+1)
      ZDATAR(JLOT,N-K)=PDAT(JLOT,K*2+1)
      ZDATAI(JLOT,K)=PDAT(JLOT,K*2+2)
      ZDATAI(JLOT,N-K) = -PDAT(JLOT,K*2+2)
    ENDDO
    ZDATAI(JLOT,0)=0._JPRB
  ENDDO

ENDIF

!
! Compute M as the smallest power of two that is greater than or equal to 2N-2
! and compute the M vector H2 from equations (2.16)
!

M=1
IPO2=0
DO WHILE( M<=2*N-2)
  M=M*2
  IPO2=IPO2+1
ENDDO

ALLOCATE(ZY(2*KLOT,0:(M/2+1)*2))

! create Y by mult with bluestein n**2

ZX0(1:KLOT) = ZDATAR(1:KLOT,0)
DO I=0,N-1

  ZR=ZDATAR(1:KLOT,I)
  ZI=ZDATAI(1:KLOT,I)

  ZWR=TB%FFTB(N)%HS(1,I,ISIGN)
  ZWI=TB%FFTB(N)%HS(2,I,ISIGN)

  DO K=1,KLOT
    ZY((K-1)*2+1,I) = ZR(K)*ZWR + ZI(K)*ZWI
    ZY((K-1)*2+2,I) = ZI(K)*ZWR - ZR(K)*ZWI
  ENDDO

ENDDO

! zero padding of Y

DO I=N,(M/2+1)*2
  ZY(:,I) = 0._JPRB
ENDDO

! FFT of Y

ILOT=2*KLOT
IINC=ILOT
IJUMP=1
IFFTSIGN=-1  ! R->C

CALL FFT992(ZY,TB%TRIGS(1,IPO2),TB%NFAX(1,IPO2),IINC,IJUMP,M,ILOT,IFFTSIGN)
CALL FFT992_CC(ZY, IINC, IJUMP, M, ILOT, IFFTSIGN)

! convert real FFT output, pointwise multiplication with h_hat(n-k) and real/imag
! swap in preparation for inverse FFT

DO I=0,M-1
  DO K=1,KLOT
    ZR(K)=ZY((K-1)*2+1,I)
    ZI(K)=ZY((K-1)*2+2,I)
  ENDDO

  ZWR=TB%FFTB(N)%H2xT(1,I,ISIGN)
  ZWI=TB%FFTB(N)%H2xT(2,I,ISIGN)

! swap
  DO K=1,KLOT
    ZY((K-1)*2+1,I) = ZI(K)*ZWR + ZR(K)*ZWI
    ZY((K-1)*2+2,I) = ZR(K)*ZWR - ZI(K)*ZWI
  ENDDO
ENDDO

! IFFT as a FFT with swapped input and swapped output

CALL FFT992(ZY,TB%TRIGS(1,IPO2),TB%NFAX(1,IPO2),IINC,IJUMP,M,ILOT,IFFTSIGN)
CALL FFT992_CC (ZY, IINC, IJUMP, M, ILOT, IFFTSIGN)

! create final by mult with another bluestein n**2 and swap output of prev FFT
! postprocessing

IF( KSIGN==-1) THEN

  DO I=0,N/2
    DO K=1,KLOT
      ZI(K)=ZY((K-1)*2+1,I)
      ZR(K)=ZY((K-1)*2+2,I)
    ENDDO

    ZWR=TB%FFTB(N)%HS(1,I,ISIGN)
    ZWI=TB%FFTB(N)%HS(2,I,ISIGN)
    IR=I*2+1
    II=I*2+2
    DO K=1,KLOT
      PDAT(K,IR) = ZR(K)*ZWR + ZI(K)*ZWI
      PDAT(K,II) = ZI(K)*ZWR - ZR(K)*ZWI
    ENDDO
  ENDDO

ELSE

  DO I=0,N-1

    DO K=1,KLOT
      ZI(K)=ZY((K-1)*2+1,I)
      ZR(K)=ZY((K-1)*2+2,I)
    ENDDO

    ZWR=TB%FFTB(N)%HS(1,I,ISIGN)
    ZWI=TB%FFTB(N)%HS(2,I,ISIGN)

    DO K=1,KLOT
      PDAT(K,I+1) = ZR(K)*ZWR + ZI(K)*ZWI
    ENDDO
  ENDDO
  DO K=1,KLOT
    PDAT(K,N) =PDAT(K,N) + ZX0(K)
  ENDDO

ENDIF

DEALLOCATE(ZY)
DEALLOCATE(ZDATAR)
DEALLOCATE(ZDATAI)

RETURN
END SUBROUTINE BLUESTEIN_FFT


!=============================================================================
SUBROUTINE BLUESTEIN_INIT(TB)
!
! Initialize data structures required by Bluestein FFT
!
!
TYPE(FFTB_TYPE),INTENT(INOUT) :: TB
INTEGER(KIND=JPIM) :: N,M,IPO2,JLAT,J,K,ISIGN,KSIGN
INTEGER(KIND=JPIM) :: ICURR,IPREV
INTEGER(KIND=JPIM) :: IJUMP,ILOT,IINC,IFFTSIGN

LOGICAL :: LLUSEFFT992
REAL(KIND=JPRB) :: DEL,ANGLE,ZSIGN

! determine number of PO2 FFT sizes needed by Bluestein FFTs
M=1
IPO2=0
DO WHILE( M<=2*TB%NDLON-2)
  M=M*2
  IPO2=IPO2+1
ENDDO

!WRITE(*,'("BLUESTEIN_INIT: 2*KLON-2=",I5," M=",I5," IPO2=",I2)')2*TB%NDLON-2,M,IPO2

! now go and generate the trigs for the above number of PO2 FFT sizes
ALLOCATE(TB%TRIGS(M,IPO2))
ALLOCATE(TB%NFAX(19,IPO2))
TB%TRIGS(:,:)=0.0D0
TB%NFAX (:,:)=0.0D0

M=1
IPO2=0
DO WHILE( M<=2*TB%NDLON-2)
  M=M*2
  IPO2=IPO2+1
  CALL SET99B(TB%TRIGS(1,IPO2),TB%NFAX(1,IPO2),M,LLUSEFFT992)
  IF( .NOT.LLUSEFFT992 )THEN
    CALL ABOR1("BLUESTEIN_INIT: UNEXPECTED LLUSEFFT992=F")
  ENDIF
ENDDO

ALLOCATE(TB%FFTB(TB%NDLON))
DO J=1,TB%NDLON
  TB%FFTB(J)%NSIZE=-1
ENDDO

DO JLAT=1,TB%NLAT_COUNT

  N=TB%NLATS(JLAT)

  IF( TB%FFTB(N)%NSIZE==N )THEN
  ! we have already initialised this latitude length
  ! WRITE(0,'("BLUESTEIN_INIT: WARNING - LATITUDE LENGTH ",I6," ALREADY INITIALIZED")')N
    CYCLE
  ENDIF

  IF( N > TB%NDLON )THEN
    CALL ABOR1("BLUESTEIN_INIT: N > TB%NDLON UNEXPECTED")
  ENDIF

  ! now set the specific PO2 (i.e. M and IPO2) for the N length of 
  ! this latitude being initialized
  M=1
  IPO2=0
  DO WHILE( M<=2*N-2)
    M=M*2
    IPO2=IPO2+1
  ENDDO

  TB%FFTB(N)%NSIZE=N


  DEL=2.0D0*ASIN(1.0D0)/REAL(N,JPRB)

  ALLOCATE(TB%FFTB(N)%HS(2,0:N-1,2))
  ALLOCATE(TB%FFTB(N)%H2xT(2,0:(M/2+1)*2,2))

  DO ISIGN=1,2

    IF( ISIGN==1 )THEN
      KSIGN=-1
    ELSE
      KSIGN= 1
    ENDIF

    ZSIGN=-KSIGN
    
    ! conjugate bluestein sequence

    DO K=0,N-1
      ANGLE=REAL(K*K,JPRB)*DEL
      TB%FFTB(N)%HS(1,K,ISIGN)=COS(ANGLE)
      TB%FFTB(N)%HS(2,K,ISIGN)=ZSIGN*SIN(ANGLE)
    ENDDO
  
    DO K=0,(M/2+1)*2
      TB%FFTB(N)%H2xT(1,K,ISIGN) = 0._JPRB
      TB%FFTB(N)%H2xT(2,K,ISIGN) = 0._JPRB
    ENDDO
    TB%FFTB(N)%H2xT(1,0,ISIGN) = TB%FFTB(N)%HS(1,0,ISIGN)
    TB%FFTB(N)%H2xT(2,0,ISIGN) = TB%FFTB(N)%HS(2,0,ISIGN)
  
    DO K=1,N-1
      TB%FFTB(N)%H2xT(1,K,ISIGN)   = TB%FFTB(N)%HS(1,K,ISIGN)
      TB%FFTB(N)%H2xT(1,M-K,ISIGN) = TB%FFTB(N)%HS(1,K,ISIGN)
      TB%FFTB(N)%H2xT(2,K,ISIGN)   = TB%FFTB(N)%HS(2,K,ISIGN)
      TB%FFTB(N)%H2xT(2,M-K,ISIGN) = TB%FFTB(N)%HS(2,K,ISIGN)
    ENDDO
    IF( M > 2*N-2 ) THEN
      DO K=N,M-N+1
        TB%FFTB(N)%H2xT(1,K,ISIGN)   = 0._JPRB
        TB%FFTB(N)%H2xT(2,K,ISIGN)   = 0._JPRB
      ENDDO
    ENDIF


    !
    ! Compute an unnormalized discrete Fourier transform of H2 -> H_hat
    !
    ILOT=2
    IINC=ILOT
    IJUMP=1
    IFFTSIGN=1  ! C->R
    CALL FFT992_CC(TB%FFTB(N)%H2xT(:,:,ISIGN),IINC,IJUMP,M,ILOT,IFFTSIGN)
    CALL FFT992(TB%FFTB(N)%H2xT(:,:,ISIGN),TB%TRIGS(1,IPO2),TB%NFAX(1,IPO2),IINC,IJUMP,M,ILOT,IFFTSIGN)

  ENDDO ! ISIGN

ENDDO ! JLAT

RETURN
END SUBROUTINE BLUESTEIN_INIT


!=============================================================================
SUBROUTINE BLUESTEIN_TERM(TB)
!
! Remove data structures used by Bluestein FFT
!
!
TYPE(FFTB_TYPE),INTENT(INOUT) :: TB
INTEGER(KIND=JPIM) :: N,JLAT

IF( ALLOCATED(TB%TRIGS) ) DEALLOCATE(TB%TRIGS)
IF( ALLOCATED(TB%NFAX)  ) DEALLOCATE(TB%NFAX)
DO JLAT=1,TB%NLAT_COUNT
  N=TB%NLATS(JLAT)
  IF( ALLOCATED(TB%FFTB(N)%HS) )   DEALLOCATE(TB%FFTB(N)%HS)
  IF( ALLOCATED(TB%FFTB(N)%H2xT) ) DEALLOCATE(TB%FFTB(N)%H2xT)
ENDDO
IF( ALLOCATED(TB%NLATS) ) DEALLOCATE(TB%NLATS)
IF( ALLOCATED(TB%FFTB)  ) DEALLOCATE(TB%FFTB)

RETURN
END SUBROUTINE BLUESTEIN_TERM


!=============================================================================

END MODULE BLUESTEIN_MOD
