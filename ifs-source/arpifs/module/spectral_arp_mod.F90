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

MODULE SPECTRAL_ARP_MOD

!   Author.
!   -------
!     Martin Leutbecher
!
!   Modifications.
!   -------------- 
!     Original      02-Jan-2009
!     M.Leutbecher  09-Jan-2009 parallelisation of rs2d/rs3d by zonal wavenumber
!     F.Bouttier    31-Jul-2012 added LAM support
!     F.Vana        17-Dec-2015 support for single precision
! ------------------------------------------------------------------

USE GEOMETRY_MOD        , ONLY : GEOMETRY
USE PARKIND1            , ONLY : JPIM, JPRB, JPRD
USE YOMHOOK             , ONLY : LHOOK, DR_HOOK, JPHOOK
USE SPECTRAL_FIELDS_MOD
USE SPECTRAL_FIELDS_DATA, ONLY : LELAM
USE RANDOM_NUMBERS_MIX

IMPLICIT NONE
PRIVATE
PUBLIC SPECTRAL_ARP, ALLOCATE_ARP, DEALLOCATE_ARP, SET_ARP, SET_ARP2D, EVOLVE_ARP, SUM_ARPS, SET_SEED_ARP

TYPE SPECTRAL_ARP
  !
  !   spectral field sf containing the current state of the auto-regressive process
  !
  TYPE(SPECTRAL_FIELD) :: SF
  !
  !   spectral field sf containing the current state of sum of the 
  !   auto-regressive processes and switch if this should be activated
  !
  LOGICAL              :: LSUM
  TYPE(SPECTRAL_FIELD) :: SFSUM
  !
  !   master seed that generates all other seeds
  !
  INTEGER(KIND=JPIM)   :: ISEED0
  !
  !   initial master seed that generates all other seeds (set by allocate_arp only)
  !
  INTEGER(KIND=JPIM)   :: ISEED0T0
  !
  !   seeds (global variable) for local randomNumberStreams
  !     (depending on [level] spherical harmonics index and field
  !
  INTEGER(KIND=JPIM),       DIMENSION(:,:),   POINTER :: ISEED2D=>NULL()
  INTEGER(KIND=JPIM),       DIMENSION(:,:,:), POINTER :: ISEED3D=>NULL() 
  !
  !   (local) randomNumberStreams for the noise used in the AR(1) processes
  !
  TYPE(RANDOMNUMBERSTREAM), DIMENSION(:,:),   POINTER :: RS2D=>NULL() 
  TYPE(RANDOMNUMBERSTREAM), DIMENSION(:,:,:), POINTER :: RS3D=>NULL() 
  !
  !   standard deviation as function of [level,] total wavenumber and field
  !
  REAL(KIND=JPRB),          DIMENSION(:,:),   POINTER :: SDEV2D=>NULL() 
  REAL(KIND=JPRB),          DIMENSION(:,:,:), POINTER :: SDEV3D=>NULL() 
  !   
  !   switch for clipping and
  !   amplitude for clipping the noise (aclip) and for clipping the AR1 (aclipar) as function of [level,] total wavenumber and field
  !
  LOGICAL :: LCLIP
  REAL(KIND=JPRB),          DIMENSION(:,:),   POINTER :: ACLIP2D=>NULL() 
  REAL(KIND=JPRB),          DIMENSION(:,:,:), POINTER :: ACLIP3D=>NULL() 
  REAL(KIND=JPRB),          DIMENSION(:,:),   POINTER :: ACLIPAR2D=>NULL() 
  REAL(KIND=JPRB),          DIMENSION(:,:,:), POINTER :: ACLIPAR3D=>NULL() 
  !
  !   correlation between successive steps as function of [level,] total wavenumber and field
  !
  REAL(KIND=JPRB),          DIMENSION(:,:),   POINTER :: PHI2D=>NULL() 
  REAL(KIND=JPRB),          DIMENSION(:,:,:), POINTER :: PHI3D=>NULL() 
END TYPE SPECTRAL_ARP


! ------------------------------------------------------------------
CONTAINS
! ------------------------------------------------------------------

SUBROUTINE ALLOCATE_ARP(YDGEOMETRY,YDARP,KS3D,KS2D,KS1D,KGRIB,KSEED,LDCLIP,LDSUM)

TYPE (GEOMETRY), INTENT(IN) :: YDGEOMETRY
TYPE (SPECTRAL_ARP), INTENT(OUT) :: YDARP
INTEGER(KIND=JPIM), INTENT(IN) :: KS3D, KS2D, KS1D, KGRIB(KS3D+KS2D), KSEED
LOGICAL, OPTIONAL, INTENT(IN) :: LDCLIP, LDSUM

INTEGER(KIND=JPIM) :: INFLEVG,INFLEVL,ISPEC2,ISPEC2G
INTEGER(KIND=JPIM) :: JS, JFLD, JL, IL
INTEGER(KIND=JPIM) :: JM, IM, IJSG, IJS, JNX, ILENX
INTEGER(KIND=JPIM) :: JJJ
LOGICAL            :: LLCLIP, LLSUM

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('SPECTRAL_ARP_MOD:ALLOCATE_ARP',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM, YDDIMV=>YDGEOMETRY%YRDIMV, &
 & YDLAP=>YDGEOMETRY%YRLAP, YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NMSMAX=>YDDIM%NMSMAX, NSMAX=>YDDIM%NSMAX, NUMP=>YDDIM%NUMP, &
 & NFLEVG=>YDDIMV%NFLEVG, NFLEVL=>YDDIMV%NFLEVL, &
 & MYMS=>YDLAP%MYMS, &
 & NALLMS=>YDMP%NALLMS, NPSURF=>YDMP%NPSURF, NPTRLL=>YDMP%NPTRLL, &
 & NPTRMS=>YDMP%NPTRMS, NUMLL=>YDMP%NUMLL)

!
!    initial checks
!
IF (PRESENT(LDCLIP)) THEN
  LLCLIP=LDCLIP
ELSE
  LLCLIP=.FALSE.
ENDIF
YDARP%LCLIP=LLCLIP

IF (PRESENT(LDSUM)) THEN
  LLSUM=LDSUM
ELSE
  LLSUM=.FALSE.
ENDIF
YDARP%LSUM=LLSUM

YDARP%ISEED0=KSEED
YDARP%ISEED0T0=KSEED ! first seed, not changed by set_seed_arp
!
!   allocate the spectral field component and
!        initialize AR(1) with 0.
!
CALL ALLOCATE_SPEC(YDARP%SF,NFLEVL,NFLEVG,NUMP,MYMS,NSMAX,NMSMAX, &
 & NALLMS,NPTRMS,NUMLL,NPTRLL,NPSURF,KS3D,KS2D,KGRIB)
YDARP%SF = 0.0_JPRB
IF (LLSUM) THEN
  CALL ALLOCATE_SPEC(YDARP%SFSUM,NFLEVL,NFLEVG,NUMP,MYMS,NSMAX,NMSMAX, &
 & NALLMS,NPTRMS,NUMLL,NPTRLL,NPSURF,MIN(1,KS3D),MIN(1,KS2D),KGRIB)
  YDARP%SFSUM = 0.0_JPRB
ENDIF
!
!   some shortcuts
!
INFLEVG=YDARP%SF%NFLEVG
INFLEVL=YDARP%SF%NFLEVL
ISPEC2 =YDARP%SF%NSPEC2
ISPEC2G=YDARP%SF%NSPEC2G
!
!   allocate the seeds (version with one random number generator for each m-wavenumber)
!
ALLOCATE(YDARP%ISEED2D(          0:NSMAX, KS2D))
ALLOCATE(YDARP%ISEED3D( INFLEVG, 0:NSMAX, KS3D))
!   allocate the random streams
!     (local to each task)
!
IF (KS2D>0) ALLOCATE(YDARP%RS2D(         NUMP, KS2D) )
IF (KS3D>0) ALLOCATE(YDARP%RS3D(INFLEVL, NUMP, KS3D) )
!
!   initialize random streams with master seed
!
CALL SET_SEED_ARP( YDARP, KSEED)
!
!   allocate arrays for standard deviation and correlation time scale
!
IF (KS2D>0) THEN
  ALLOCATE( YDARP%SDEV2D( 0:NSMAX, KS2D ))
  ALLOCATE( YDARP%PHI2D(  0:NSMAX, KS2D ))
  IF (LLCLIP) THEN
    ALLOCATE(YDARP%ACLIP2D(0:NSMAX, KS2D ))
    ALLOCATE(YDARP%ACLIPAR2D(0:NSMAX, KS2D ))
  ENDIF
ENDIF
IF (KS3D>0) THEN
  ALLOCATE( YDARP%SDEV3D(INFLEVG,0:NSMAX, KS3D ))
  ALLOCATE( YDARP%PHI3D( INFLEVG,0:NSMAX, KS3D ))
  IF (LLCLIP) THEN 
    ALLOCATE( YDARP%ACLIP3D(INFLEVG,0:NSMAX, KS3D ))
    ALLOCATE( YDARP%ACLIPAR3D(INFLEVG,0:NSMAX, KS3D ))
  ENDIF
ENDIF

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SPECTRAL_ARP_MOD:ALLOCATE_ARP',1,ZHOOK_HANDLE)
END SUBROUTINE ALLOCATE_ARP
!
!
!
SUBROUTINE SET_SEED_ARP(YDARP,KSEED)

TYPE (SPECTRAL_ARP), INTENT(INOUT) :: YDARP
INTEGER(KIND=JPIM), INTENT(IN) ::  KSEED

INTEGER(KIND=JPIM) :: INFLEVG,INFLEVL,ISPEC2,ISPEC2G
INTEGER(KIND=JPIM) :: JS, JFLD, JL, IL, ISMAX,IS2D,IS3D
INTEGER(KIND=JPIM) :: JM, IM, IJSG, IJS, JNX, ILENX
INTEGER(KIND=JPIM) :: JJJ
LOGICAL            :: LLCLIP

TYPE(RANDOMNUMBERSTREAM) :: YL_SEEDS
REAL(KIND=JPRB), DIMENSION(:), ALLOCATABLE :: ZSEED
REAL(KIND=JPRD), PARAMETER :: ZHUGE = REAL(HUGE(1_JPIM)-1,JPRD)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('SPECTRAL_ARP_MOD:SET_SEED_ARP',0,ZHOOK_HANDLE)

!
!    initial checks
!
YDARP%ISEED0=KSEED
!
!   some shortcuts
!
INFLEVG=YDARP%SF%NFLEVG
INFLEVL=YDARP%SF%NFLEVL
ISPEC2 =YDARP%SF%NSPEC2
ISPEC2G=YDARP%SF%NSPEC2G
ISMAX  =YDARP%SF%NSMAX
IS2D   =YDARP%SF%NS2D
IS3D   =YDARP%SF%NS3D
!
!   generate seeds (global)
!
CALL INITIALIZE_RANDOM_NUMBERS( KSEED, YL_SEEDS )
IF (IS2D>0) THEN
  ALLOCATE(ZSEED( (ISMAX+1)*IS2D ))
  CALL UNIFORM_DISTRIBUTION (ZSEED, YL_SEEDS)
  YDARP%ISEED2D= FLOOR(RESHAPE(ZSEED, (/(ISMAX+1), IS2D /) ) * ZHUGE)
  DEALLOCATE(ZSEED)
ENDIF
IF (IS3D>0) THEN
  ALLOCATE(ZSEED( INFLEVG*(ISMAX+1)*IS3D ))
  CALL UNIFORM_DISTRIBUTION (ZSEED, YL_SEEDS)
  YDARP%ISEED3D= FLOOR(RESHAPE(ZSEED, (/INFLEVG, ISMAX+1, IS3D /) ) * ZHUGE)
  DEALLOCATE(ZSEED)
ENDIF
!
!   (re-)initialise the random streams for 2D fields
!
IF (IS2D>0) THEN
  DO JFLD=1, IS2D
    DO JM=1,YDARP%SF%NUMP
      IM=YDARP%SF%MYMS(JM)
      CALL INITIALIZE_RANDOM_NUMBERS( YDARP%ISEED2D(IM,JFLD), YDARP%RS2D(JM, JFLD))
    ENDDO
  ENDDO
ENDIF
!
!   initialise the random streams for 3D fields
!
IF (IS3D>0) THEN
  DO JFLD=1, IS3D
    DO JM=1,YDARP%SF%NUMP
      IM=YDARP%SF%MYMS(JM)
      DO JL=1,YDARP%SF%NFLEVL
        IL=YDARP%SF%MYLEVS(JL)
        CALL INITIALIZE_RANDOM_NUMBERS( YDARP%ISEED3D(IL,IM,JFLD), YDARP%RS3D(JL,JM, JFLD))
      ENDDO
    ENDDO
  ENDDO
ENDIF

IF (LHOOK) CALL DR_HOOK('SPECTRAL_ARP_MOD:SET_SEED_ARP',1,ZHOOK_HANDLE)

END SUBROUTINE SET_SEED_ARP
  
!
!
!
SUBROUTINE DEALLOCATE_ARP(YDARP)
! ----------------------------------------------------------------------
! deallocate components of ydarp
! ----------------------------------------------------------------------
TYPE (SPECTRAL_ARP), INTENT(INOUT) :: YDARP

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('SPECTRAL_ARP_MOD:DEALLOCATE_ARP',0,ZHOOK_HANDLE)
CALL DEALLOCATE_SPEC( YDARP%SF )
IF (YDARP%LSUM) THEN
  CALL DEALLOCATE_SPEC( YDARP%SFSUM )
ENDIF

IF (ASSOCIATED( YDARP%ISEED2D )) DEALLOCATE(YDARP%ISEED2D)
IF (ASSOCIATED( YDARP%ISEED3D )) DEALLOCATE(YDARP%ISEED3D)

IF (ASSOCIATED( YDARP%RS2D    )) DEALLOCATE(YDARP%RS2D)
IF (ASSOCIATED( YDARP%RS3D    )) DEALLOCATE(YDARP%RS3D)

IF (ASSOCIATED( YDARP%SDEV2D  )) DEALLOCATE(YDARP%SDEV2D)
IF (ASSOCIATED( YDARP%SDEV3D  )) DEALLOCATE(YDARP%SDEV3D)

IF (ASSOCIATED( YDARP%ACLIP2D  )) DEALLOCATE(YDARP%ACLIP2D)
IF (ASSOCIATED( YDARP%ACLIP3D  )) DEALLOCATE(YDARP%ACLIP3D)

IF (ASSOCIATED( YDARP%ACLIPAR2D  )) DEALLOCATE(YDARP%ACLIPAR2D)
IF (ASSOCIATED( YDARP%ACLIPAR3D  )) DEALLOCATE(YDARP%ACLIPAR3D)

IF (ASSOCIATED( YDARP%PHI2D   )) DEALLOCATE(YDARP%PHI2D)
IF (ASSOCIATED( YDARP%PHI3D   )) DEALLOCATE(YDARP%PHI3D)

IF (LHOOK) CALL DR_HOOK('SPECTRAL_ARP_MOD:DEALLOCATE_ARP',1,ZHOOK_HANDLE)
END SUBROUTINE DEALLOCATE_ARP
!
!
!
SUBROUTINE SET_ARP2D(YDARP,SDEV2D,PHI2D,ACLIP2D)
! ----------------------------------------------------------------------
! set standard deviation and correlation 
! ----------------------------------------------------------------------

TYPE (SPECTRAL_ARP), INTENT(INOUT) :: YDARP
REAL(KIND=JPRB), DIMENSION(:,:),   INTENT(IN) :: SDEV2D, PHI2D
REAL(KIND=JPRB), DIMENSION(:,:),   OPTIONAL, INTENT(IN) :: ACLIP2D

INTEGER(KIND=JPIM) :: JN , JS
REAL(KIND=JPRB)    :: ZPHI_FAC0, ZPHI, ZEPSILON

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('SPECTRAL_ARP_MOD:SET_ARP2D',0,ZHOOK_HANDLE)

ZEPSILON=100._JPRB*TINY(1._JPRB)
YDARP%SDEV2D=SDEV2D
YDARP%PHI2D =PHI2D
IF (PRESENT(ACLIP2D)) THEN
  YDARP%ACLIP2D =ACLIP2D
ELSE
  IF (YDARP%LCLIP) THEN
    YDARP%ACLIP2D=10._JPRB * YDARP%SDEV2D
  ENDIF
ENDIF


IF (YDARP%LCLIP) THEN
  DO JS=LBOUND(YDARP%SDEV2D,DIM=2), UBOUND(YDARP%SDEV2D,DIM=2)
    DO JN=LBOUND(YDARP%SDEV2D,DIM=1), UBOUND(YDARP%SDEV2D,DIM=1)
      ZPHI=YDARP%PHI2D( JN, JS)
      ZPHI_FAC0 = 1._JPRB / MAX(SQRT( 1._JPRB - ZPHI*ZPHI),ZEPSILON)
      YDARP%ACLIPAR2D(JN,JS) = ZPHI_FAC0 * YDARP%ACLIP2D(JN,JS)
    ENDDO
  ENDDO
ENDIF

IF (LHOOK) CALL DR_HOOK('SPECTRAL_ARP_MOD:SET_ARP2D',1,ZHOOK_HANDLE)

END SUBROUTINE SET_ARP2D
!
!
!
SUBROUTINE SET_ARP(YDARP,SDEV2D,SDEV3D,PHI2D,PHI3D,ACLIP2D,ACLIP3D)
! ----------------------------------------------------------------------
! set standard deviation and correlation 
! ----------------------------------------------------------------------

TYPE (SPECTRAL_ARP), INTENT(INOUT) :: YDARP
REAL(KIND=JPRB), DIMENSION(:,:),   INTENT(IN) :: SDEV2D, PHI2D
REAL(KIND=JPRB), DIMENSION(:,:,:), INTENT(IN) :: SDEV3D, PHI3D
REAL(KIND=JPRB), DIMENSION(:,:),    OPTIONAL,INTENT(IN) :: ACLIP2D
REAL(KIND=JPRB), DIMENSION(:,:,:),  OPTIONAL,INTENT(IN) :: ACLIP3D

INTEGER(KIND=JPIM) :: JN , JS, JLEV
REAL(KIND=JPRB)    :: ZPHI_FAC0, ZPHI

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('SPECTRAL_ARP_MOD:SET_ARP',0,ZHOOK_HANDLE)

YDARP%SDEV2D=SDEV2D; YDARP%SDEV3D=SDEV3D
YDARP%PHI2D =PHI2D;  YDARP%PHI3D =PHI3D
IF (PRESENT(ACLIP2D)) THEN
  YDARP%ACLIP2D=ACLIP2D
ELSE
  IF (YDARP%LCLIP) THEN
    YDARP%ACLIP2D=10._JPRB * YDARP%SDEV2D
  ENDIF
ENDIF

IF (YDARP%LCLIP) THEN
  DO JS=LBOUND(YDARP%SDEV2D,DIM=2), UBOUND(YDARP%SDEV2D,DIM=2)
    DO JN=LBOUND(YDARP%SDEV2D,DIM=1), UBOUND(YDARP%SDEV2D,DIM=1)
      ZPHI=YDARP%PHI2D( JN, JS)
      ZPHI_FAC0 = 1._JPRB / SQRT( 1._JPRB - ZPHI*ZPHI)
      YDARP%ACLIPAR2D(JN,JS) = ZPHI_FAC0 * YDARP%ACLIP2D(JN,JS)
    ENDDO
  ENDDO
ENDIF

IF (PRESENT(ACLIP3D)) THEN
  YDARP%ACLIP3D=ACLIP3D
ELSE
  IF (YDARP%LCLIP) THEN
    YDARP%ACLIP3D=10._JPRB * YDARP%SDEV3D
  ENDIF
ENDIF

IF (YDARP%LCLIP) THEN
  DO JS=LBOUND(YDARP%SDEV3D,DIM=3), UBOUND(YDARP%SDEV3D,DIM=3)
    DO JN=LBOUND(YDARP%SDEV3D,DIM=2), UBOUND(YDARP%SDEV3D,DIM=2)
      DO JLEV=LBOUND(YDARP%SDEV3D,DIM=1), UBOUND(YDARP%SDEV3D,DIM=1)
        ZPHI=YDARP%PHI3D(JLEV, JN, JS)
        ZPHI_FAC0 = 1._JPRB / SQRT( 1._JPRB - ZPHI*ZPHI)
        YDARP%ACLIPAR3D(JLEV,JN,JS) = ZPHI_FAC0 * YDARP%ACLIP3D(JLEV,JN,JS)
      ENDDO
    ENDDO
  ENDDO
ENDIF

IF (LHOOK) CALL DR_HOOK('SPECTRAL_ARP_MOD:SET_ARP',1,ZHOOK_HANDLE)

END SUBROUTINE SET_ARP
!
!
!
SUBROUTINE EVOLVE_ARP(YDARP,KARP,LDINIT)
! ----------------------------------------------------------------------
! evolve the autoregressive process by karp steps
! ----------------------------------------------------------------------
TYPE (SPECTRAL_ARP),          INTENT(INOUT) :: YDARP
INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN)    :: KARP
LOGICAL, OPTIONAL,            INTENT(IN)    :: LDINIT
!
!   local variables
!
LOGICAL :: LLINIT
INTEGER(KIND=JPIM) :: ISTEPS,JSTEP
INTEGER(KIND=JPIM) :: JMLOC, IM,JIR,JN, INM, ILEN1, ILEN2, IDXNAR, IRMAX, INMAX, INN
INTEGER(KIND=JPIM) :: JLEV,JFLD,ILEV
REAL(KIND=JPRB)    :: ZNOISE, ZPHI, ZPHI_FAC0, ZLIMIT, ZEPSILON
REAL(KIND=JPRB), DIMENSION(:,:),   ALLOCATABLE :: ZNAR2D
REAL(KIND=JPRB), DIMENSION(:,:,:), ALLOCATABLE :: ZNAR3D
!
#include "fkstar.func.h"
!
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('SPECTRAL_ARP_MOD:EVOLVE_ARP',0,ZHOOK_HANDLE)

ZEPSILON=100._JPRB*TINY(1.0_JPRB)
IF (PRESENT(LDINIT)) THEN
  LLINIT=LDINIT
ELSE
  LLINIT=.FALSE.
ENDIF
IF (PRESENT(KARP)) THEN
  ISTEPS=KARP
ELSE
  ISTEPS=1
ENDIF

IF (LELAM) THEN
  ALLOCATE(ZNAR2D(                 YDARP%SF%NS2D,4*(1+YDARP%SF%NSMAX)))
  ALLOCATE(ZNAR3D(YDARP%SF%NFLEVL, YDARP%SF%NS3D,4*(1+YDARP%SF%NSMAX)))
ELSE
  ALLOCATE(ZNAR2D(                 YDARP%SF%NS2D,2*(1+YDARP%SF%NSMAX)))
  ALLOCATE(ZNAR3D(YDARP%SF%NFLEVL, YDARP%SF%NS3D,2*(1+YDARP%SF%NSMAX)))
ENDIF

DO JSTEP=1, ISTEPS
  !
  !   generate new noise and add to spectral field
  !
  DO JMLOC=1,YDARP%SF%NUMP
    IM=YDARP%SF%MYMS(JMLOC)
    !
    !   generate all random numbers for M=IM and N=M, M+1, ..., NSMAX
    !    real and imaginary part (factor min(1,im))
    !
    !   To avoid code duplication, the innermost part of the loop is shareed by global and lam models,
    !   only the spectral indexing changes: [FBt]
    !                       in global                    in lam
    !      jir:          scans 1 or 2 coeffs       scans 4 coeffs: sinx,cosx,siny,cosy
    !      n-index:      goes from m to nsmax      goes from 0 to the spectral m-truncation f(n) i.e. the ellipse width
    !      totalwavenum: is n                      is fkstar(n,m) as in Jb
    !
    IF (LELAM) THEN
      ILEN2=YDARP%SF%NCPL4M(IM) ! number of real and imaginary coefficients for M=im
    ELSE
      ILEN1=YDARP%SF%NSMAX-IM+1 ! number of N wavenumbers for M=im
      ILEN2=(1+MIN(1,IM))*ILEN1 ! number of real and imaginary coefficients for M=im
    ENDIF
    DO JFLD=1,YDARP%SF%NS2D
      CALL GAUSSIAN_DISTRIBUTION(  ZNAR2D(JFLD,     1:ILEN2), YDARP%RS2D(JMLOC,     JFLD) )
    ENDDO
    DO JFLD=1,YDARP%SF%NS3D
      DO JLEV=1,YDARP%SF%NFLEVL
        CALL GAUSSIAN_DISTRIBUTION(ZNAR3D(JLEV,JFLD,1:ILEN2), YDARP%RS3D(JLEV,JMLOC,JFLD) )
      ENDDO
    ENDDO
    IF(LELAM)THEN
      IRMAX=3                  ! lam : 4 spectral coeffs per (n,m)
    ELSE
      IRMAX=MIN(1,IM)          ! global: 1 or 2 spectral coeffs per (n,m)
    ENDIF
    DO JIR=0,IRMAX
      IF(LELAM)THEN
        INMAX=YDARP%SF%NCPL4M(IM)/4 - 1  ! lam : ncpl4(m)/4 complex coeffs for each m
      ELSE
        INMAX=YDARP%SF%NSMAX-IM          ! global: nsmax-m+1 complex coeffs for each m
      ENDIF
      DO JN=0,INMAX   
        IF(LELAM) THEN
          INN=FKSTAR(JN,IM,YDARP%SF%NSMAX,YDARP%SF%NMSMAX)   ! lam : total wavenumber=f(n,m)
          INM=YDARP%SF%NESM0(IM)+JN*4+JIR
          IDXNAR=1+JN*4+JIR
        ELSE
          INN=JN+IM                                          ! global: total wavenumber
          INM=YDARP%SF%NASM0(IM)+JN*2+JIR
          IDXNAR=1+JN*(1+MIN(1,IM))+JIR
        ENDIF
        !
        !   2D fields
        !
        DO JFLD=1,YDARP%SF%NS2D

          ZNOISE = YDARP%SDEV2D( INN, JFLD ) * ZNAR2D(JFLD, IDXNAR)
          !
          !   clipping of forcing noise
          !
          IF (YDARP%LCLIP) THEN
            ZLIMIT=YDARP%ACLIP2D( INN, JFLD )
            ZNOISE=MAX(-ZLIMIT,MIN(ZLIMIT, ZNOISE))
          ENDIF
          ZPHI   = YDARP%PHI2D(  INN, JFLD )
          IF (LLINIT) THEN
            ZPHI_FAC0 = 1._JPRB / MAX(SQRT( 1._JPRB - ZPHI*ZPHI),ZEPSILON) ! stdev ratio
            !                                            between AR(1) and forcing noise
            ZNOISE = ZNOISE * ZPHI_FAC0 
          ENDIF
          YDARP%SF%SP2D( INM, JFLD ) = ZPHI * YDARP%SF%SP2D( INM, JFLD ) + ZNOISE
          !
          !   clipping of AR(1)
          !
          IF (YDARP%LCLIP) THEN
            ZLIMIT=YDARP%ACLIPAR2D( INN, JFLD )
            YDARP%SF%SP2D( INM, JFLD ) = MAX(-ZLIMIT, MIN(ZLIMIT,YDARP%SF%SP2D( INM, JFLD )))
          ENDIF
        ENDDO
        !
        !   3D fields
        !
        DO JFLD=1,YDARP%SF%NS3D
          DO JLEV=1,YDARP%SF%NFLEVL
            ILEV=YDARP%SF%MYLEVS(JLEV)
            ZNOISE = YDARP%SDEV3D( ILEV, INN, JFLD ) * ZNAR3D(JLEV,JFLD,IDXNAR)
            !
            !   clipping of forcing noise
            !
            IF (YDARP%LCLIP) THEN
              ZLIMIT=YDARP%ACLIP3D( ILEV, INN, JFLD )
              ZNOISE=MAX(-ZLIMIT,MIN(ZLIMIT, ZNOISE))
            ENDIF
            ZPHI   = YDARP%PHI3D(  ILEV, INN, JFLD )
            IF (LLINIT) THEN
              ZPHI_FAC0 = 1._JPRB / MAX(SQRT( 1._JPRB - ZPHI*ZPHI),ZEPSILON) ! stdev ratio
              !                                            between AR(1) and forcing noise
              ZNOISE = ZNOISE * ZPHI_FAC0 
            ENDIF
            YDARP%SF%SP3D( JLEV, INM, JFLD ) = ZPHI * YDARP%SF%SP3D( JLEV, INM, JFLD ) + ZNOISE
            !
            !   clipping of AR(1)
            !
            IF (YDARP%LCLIP) THEN
              ZLIMIT=YDARP%ACLIPAR3D( ILEV, INN, JFLD )
              YDARP%SF%SP3D( JLEV, INM, JFLD ) = MAX(-ZLIMIT, MIN(ZLIMIT,YDARP%SF%SP3D( JLEV, INM, JFLD )))
            ENDIF
          ENDDO ! JLEV
        ENDDO ! JFLD
      ENDDO ! JN
    ENDDO ! JIR
  ENDDO ! JMLOC
ENDDO ! JSTEP








DEALLOCATE(ZNAR2D)
DEALLOCATE(ZNAR3D)

IF (LHOOK) CALL DR_HOOK('SPECTRAL_ARP_MOD:EVOLVE_ARP',1,ZHOOK_HANDLE)
END SUBROUTINE EVOLVE_ARP
!
!
!
SUBROUTINE SUM_ARPS(YDARP)
! ----------------------------------------------------------------------
! sum autoregressive processes
! ----------------------------------------------------------------------
TYPE (SPECTRAL_ARP),          INTENT(INOUT) :: YDARP
!
!   local variables
!
INTEGER(KIND=JPIM) :: JMLOC, IM,JIR,JN, INM
INTEGER(KIND=JPIM) :: JS,JL,JF,ISPEC2,IFLEVL

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('SPECTRAL_ARP_MOD:SUM_ARPS',0,ZHOOK_HANDLE)

IF (.NOT.YDARP%LSUM) THEN
  CALL ABOR1('ydarp%lsum has to be true in order to use sum_arps.')
ENDIF

ISPEC2=YDARP%SF%NSPEC2
IFLEVL=YDARP%SF%NFLEVL

IF (YDARP%SF%NS2D > 0) THEN
  YDARP%SFSUM%SP2D(1:ISPEC2,1)=YDARP%SF%SP2D(1:ISPEC2,1)
ENDIF

!yt!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jf,js)
DO JF=2,YDARP%SF%NS2D
  DO JS=1,ISPEC2
    YDARP%SFSUM%SP2D(JS,1) = YDARP%SFSUM%SP2D(JS,1) + YDARP%SF%SP2D(JS,JF)
  ENDDO
ENDDO
!yt!$OMP END PARALLEL DO

IF (YDARP%SF%NS3D > 0) THEN
  YDARP%SFSUM%SP3D(1:IFLEVL,1:ISPEC2,1)=YDARP%SF%SP3D(1:IFLEVL,1:ISPEC2,1)
ENDIF
!yt!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jf,js,jl)
DO JF=2,YDARP%SF%NS3D
  DO JS=1,ISPEC2
    DO JL=1,IFLEVL
      YDARP%SFSUM%SP3D(JL,JS,1) = YDARP%SFSUM%SP3D(JL,JS,1) + YDARP%SF%SP3D(JL,JS,JF)
    ENDDO
  ENDDO
ENDDO
!yt!$OMP END PARALLEL DO

IF (LHOOK) CALL DR_HOOK('SPECTRAL_ARP_MOD:SUM_ARPS',1,ZHOOK_HANDLE)

END SUBROUTINE SUM_ARPS

! ------------------------------------------------------------------

END MODULE SPECTRAL_ARP_MOD
