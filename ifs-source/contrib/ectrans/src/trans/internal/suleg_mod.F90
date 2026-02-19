! (C) Copyright 1987- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE SULEG_MOD
#ifdef __NEC__
#define SIZEOF(x) STORAGE_SIZE(x)/KIND(x)
#endif
CONTAINS
SUBROUTINE SULEG
!DEC$ OPTIMIZE:1

USE PARKIND1  ,ONLY : JPRD, JPIM, JPRB
USE PARKIND2  ,ONLY : JPRH
USE MPL_MODULE

USE TPM_GEN
USE TPM_DIM
USE TPM_CONSTANTS
USE TPM_DISTR
USE TPM_FIELDS
USE TPM_FLT
USE TPM_GEOMETRY
USE TPM_CTL
USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS

USE PRE_SULEG_MOD
USE SUGAW_MOD
USE SUPOL_MOD
USE SUPOLF_MOD
USE TPM_POL
USE SUTRLE_MOD
USE SETUP_GEOM_MOD
USE BUTTERFLY_ALG_MOD
USE SEEFMM_MIX
USE SET2PE_MOD
USE ABORT_TRANS_MOD
USE PREPSNM_MOD     ,ONLY : PREPSNM
USE WRITE_LEGPOL_MOD
USE READ_LEGPOL_MOD

!**** *SULEG * - initialize the Legendre polynomials

!     Purpose.
!     --------
!           Initialize COMMON YOMLEG

!**   Interface.
!     ----------
!        *CALL* *SULEG*

!        Explicit arguments :
!        --------------------

!        Implicit arguments :
!        --------------------
!              COMMON YOMLEG

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!        SUGAW  (Gaussian latitudes)
!        SUPOLM (polynomials)
!        LFI routines for external IO's
!        Called by SUGEM.

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS
!     
!     S.L. Belousov, Tables of normalized associated Legendre Polynomials, Pergamon Press (1962)
!     P.N. Swarztrauber, On computing the points and weights for Gauss-Legendre quadrature,
!     SIAM J. Sci. Comput. Vol. 24 (3) pp. 945-954 (2002)

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 87-10-15
!        MODIFICATION : 91-04  J.M. Piriou:
!                       - Read gaussian latitudes and PNM on LFI
!                       - If file missing, computes
!                       91-04 M.Hamrud:
!                       - IO Scheme introduced
!        MODIFICATION : 91-07-03 P.Courtier suppress derivatives
!        MODIFICATION : 91-07-03 P.Courtier computes RATATH and RACTHE
!        MODIFICATION : 91-07-03 P.Courtier change upper limit (NSMAX+1)
!        MODIFICATION : 91-07-03 P.Courtier change ordering
!     Order of the PNM in the file, as in the model :
!         - increasing wave numbers m
!         - for a given m, from n=NSMAX+1 to m
!        MODIFICATION : 92-07-02 R. Bubnova: shift RATATH calculation
!                                            to SUGEM1
!        MODIFICATION : 92-12-17 P.Courtier multitask computations
!        Modified by R. EL Khatib : 93-04-02 Set-up defaults controled by LECMWF
!        MODIFICATION : 93-03-19 D.Giard : n <= NTMAX
!        K. YESSAD    : 93-05-11 : DLMU --> global array DRMU(NDGSA:NDGEN).
!                       (not stored currently on LFI files).
!        MODIFICATION : 94-02-03 R. El Khatib : subroutine SULEG2 to write out
!                       the Leg. polynomials on workfile or LFI file
!        Modification : 94-08-31 M. Tolstykh: Setup for CUD interpolation
!        Modified by K. YESSAD (MARCH 1995): Extra-latitudes computations
!                 according to value of NDGSUR and LRPOLE only.
!                 + change fancy loop numbering.
!        Modified 98-08-10 by K. YESSAD: removal of LRPOLE option.
!          - removal of LRPOLE in YOMCT0.
!          - removal of code under LRPOLE.
!        R. El Khatib 11-Apr-2007 Emulation of vectorized quadruple precision
!                                 on NEC
!        Nils Wedi + Mats Hamrud, 2009-02-05 revised following Swarztrauber, 2002
!        G.Mozdzynski: March 2011 Support 2D (RW,RV) initialisation of legendre coeffs
!        G.Mozdzynski: July 2012  distribute FLT initialisation over NPRTRV
!        R. El Khatib 14-Jun-2013 optional computation on the stretched latitudes
!        F. Vana  05-Mar-2015  Support for single precision
!        Nils Wedi, 20-Apr-2015 Support dual latitude/longitude set
!        T. Wilhelmsson, 22-Sep-2016 Support single precision for dual too
!     ------------------------------------------------------------------

IMPLICIT NONE

!     LOCAL 
!     ------------------------------------------------------------------
REAL(KIND=JPRD),ALLOCATABLE :: ZPNMG(:)
REAL(KIND=JPRD),ALLOCATABLE :: ZFN(:,:)
REAL(KIND=JPRD),ALLOCATABLE :: ZLRMUZ2(:)
REAL(KIND=JPRB) :: ZEPSNM(0:R%NTMAX+2)
REAL(KIND=JPRD) :: ZLRMUZ(R%NDGL)
REAL(KIND=JPRD) :: ZW(R%NDGL)

REAL(KIND=JPRD) :: ZANM
REAL(KIND=JPRD) :: ZFNN
REAL(KIND=JPRD) :: ZPI, ZINC, ZOFF, ZTEMP, ZORIG, ZTHETA, ZCOS

REAL(KIND=JPRD), ALLOCATABLE :: ZSNDBUFV(:),ZRCVBUFV(:,:)
REAL(KIND=JPRD), ALLOCATABLE :: ZPNMCDO(:,:),ZPNMCDD(:,:)
REAL(KIND=JPRB), ALLOCATABLE :: ZRCVBUTFV(:,:)
REAL(KIND=JPRB) :: ZDUM(2)
REAL(KIND=KIND(ZRCVBUTFV)) :: ZBYTES
INTEGER(KIND=JPIM) :: IBYTES
INTEGER(KIND=JPIM) :: ISENDREQ(NPRTRV)
INTEGER(KIND=JPIM) :: IRECVREQ(NPRTRV)
INTEGER(KIND=JPIM) :: IKOUNT(NPRTRV)
INTEGER(KIND=JPIM) :: IRECVLENMAXV(NPRTRV)

INTEGER(KIND=JPIM) :: INM, IM, IRECV, ISEND, ISREQ, IRREQ, &
             &JGL,  JM, JMLOC, IMLOC, JN, JNM, IODD, INN, INMAX, JI, IMAXN, ITAG, ITAG1, &
             &INX, ISL, ISTART, ITHRESHOLD, INSMAX, IMAXCOLS,ILATSMAX,JW,JV,J, &
             &IDGLU, ILA, ILS, IA, IS, I, ILATS, ILOOP, IPRTRV, JSETV, JH, &
             &IMAXRECVA, IMAXRECVS, IRECVLENMAX, ICLONELEN, IHEMIS, INNH, IGL, IGL1, IGL2, &
             &IDGLU2, ISYM, INZ

REAL(KIND=JPRD) :: ZEPS_INT_DEC
REAL(KIND=JPRD) :: ZEPS
REAL(KIND=JPRD),ALLOCATABLE :: ZLFPOL(:,:)
REAL(KIND=JPRD),ALLOCATABLE :: ZLPOL(:)

TYPE(CLONE),ALLOCATABLE :: ZCLONEA(:),ZCLONES(:)

LOGICAL :: LLP1,LLP2

! For latitudes on the stretched geometry
REAL(KIND=JPRH) :: ZTAN
REAL(KIND=JPRH) :: ZSTRETMU(R%NDGL)

!     ------------------------------------------------------------------

!*       0.    Some initializations.
!              ---------------------

IBYTES = MPL_BYTES(ZBYTES)

ZEPS = 1000._JPRD*EPSILON(ZEPS)
!ZEPS_INT_DEC = EPSILON(ZEPS)
ZEPS_INT_DEC = 1.0E-7_JPRD
!ZEPS_INT_DEC = 1.0E-5_JPRD

IHEMIS=1
IF (S%LSOUTHPNM) IHEMIS=2
LLP1 = NPRINTLEV>0
LLP2 = NPRINTLEV>1
IF(LLP1) WRITE(NOUT,*) '=== ENTER ROUTINE SULEG ==='

IF( NPROC > 1 )THEN
  CALL GSTATS(798,0)
  CALL MPL_BARRIER(CDSTRING='SULEG:')
  CALL GSTATS(798,1)
ENDIF

CALL GSTATS(140,0)
CALL GSTATS(1801,0)

IF(.NOT.D%LGRIDONLY) THEN
  CALL PRE_SULEG
ENDIF

ALLOCATE(F%RMU(R%NDGL))
IF (LLP2) WRITE(NOUT,9) 'F%RMU     ',SIZE(F%RMU ),SHAPE(F%RMU ) 
ALLOCATE(F%RW(R%NDGL))
IF (LLP2) WRITE(NOUT,9) 'F%RW      ',SIZE(F%RW  ),SHAPE(F%RW  ) 


!*       1.0 Initialize Fourier coefficients for ordinary Legendre polynomials
!     ------------------------------------------------------------------------

ALLOCATE(ZFN(0:R%NDGL,0:R%NDGL))
IF (LLP2) WRITE(NOUT,9) 'ZFN       ',SIZE(ZFN   ),SHAPE(ZFN   )



! determines the number of stripes in butterfly NSMAX/IMAXCOLS 
! IMAXCOLS = (R%NSMAX - 1)/4 + 1
! IMAXCOLS=64 (min flops)
IMAXCOLS=64

! the threshold of efficiency
IF(NPROC == 1 .OR. R%NDGNH <= 2560) THEN
   ITHRESHOLD = R%NDGNH/4
   DO
      IF(ITHRESHOLD >= IMAXCOLS*4) EXIT
      IMAXCOLS = IMAXCOLS/2
   ENDDO
ELSE
   ITHRESHOLD = 900
ENDIF

ITHRESHOLD = MAX(ITHRESHOLD,IMAXCOLS+1)
S%ITHRESHOLD = ITHRESHOLD

!*       3.1   Gaussian latitudes and weights
!     ---------------------------------------

!IF (S%LUSE_BELUSOV.OR.S%LUSEFLT) 
CALL INI_POL(R%NTMAX+3)

IF(.NOT.D%LGRIDONLY) THEN
  ISTART=1
ELSE
  ISTART=R%NDGL
ENDIF

INMAX=R%NDGL
! Belousov, Swarztrauber use ZFN(0,0)=SQRT(2._JPRD)
! IFS normalisation chosen to be 0.5*Integral(Pnm**2) = 1
ZFN(0,0)=2._JPRD
DO JN=ISTART,R%NDGL
  ZFNN=ZFN(0,0)
  DO JGL=1,JN
    ZFNN=ZFNN*SQRT(1._JPRD-0.25_JPRD/REAL(JGL**2,JPRD))
  ENDDO

  IODD=MOD(JN,2)
  ZFN(JN,JN)=ZFNN
  DO JGL=2,JN-IODD,2
    ZFN(JN,JN-JGL)=ZFN(JN,JN-JGL+2)*REAL((JGL-1)*(2*JN-JGL+2),JPRD)/REAL(JGL*(2*JN-JGL+1),JPRD)
  ENDDO
ENDDO

! compute latitudes and weights for original Gaussian latitudes
ZANM=SQRT(REAL(2*INMAX+1,JPRD)*REAL(INMAX**2,JPRD)/REAL(2*INMAX-1,JPRD))
INN=R%NDGL
CALL GSTATS(1801,2)
CALL SUGAW(INN,0,INMAX,ZLRMUZ(1:INN),ZW(1:INN),ZANM,ZFN)
CALL GSTATS(1801,3)

IF (ABS(G%RSTRET-1.0_JPRD)>100._JPRD*EPSILON(1._JPRD)) THEN
  WRITE(NOUT,*) '=== SULEG: Change Gaussian latitudes to the transformed sphere ==='
  INNH=(INN+1)/2
  ZTAN=(1.0_JPRD-G%RSTRET**2)/(1.0_JPRD+G%RSTRET**2)
! North hemisphere
  DO JGL=1,INNH
    ZSTRETMU(JGL)=(ZTAN+REAL(ZLRMUZ(JGL),JPRH))/(1.0_JPRD+ZTAN*REAL(ZLRMUZ(JGL),JPRH))
  ENDDO
! South hemisphere
  DO JGL=1,INNH
    IGL=2*INNH-JGL+1
    ZSTRETMU(IGL)=(ZTAN-REAL(ZLRMUZ(JGL),JPRH))/(1.0_JPRD-ZTAN*REAL(ZLRMUZ(JGL),JPRH))
  ENDDO
  DO JGL=1,INN
    ZLRMUZ(JGL)=REAL(ZSTRETMU(JGL),JPRD)
  ENDDO
ENDIF

DO JGL=1,R%NDGL
  F%RW(JGL)     = ZW(JGL)
  F%RMU(JGL)    = ZLRMUZ(JGL)
ENDDO

IF (LLP1) WRITE(NOUT,*) '=== SULEG: Finished Gaussian latitudes ==='

!*       3.1.1  specify a dual set of output (inv_trans) or input (dir_trans) latitudes / longitudes

IF( S%LDLL ) THEN

  INMAX = S%NDGL
  INN= S%NDGL

  S%NDGNHD = (INMAX+1)/2
  ALLOCATE(ZLRMUZ2(INN))

  ! here we want to use the positions of the specified dual grid
  ! accuracy requirement is ZLRMUZ2(JGL) < F%RMU(1)
  ! so we use approximations for the remaining latitudes outside this range
  ! we approximate the vicinity to the pole/equator

  ZPI = 2.0_JPRD*ASIN(1.0_JPRD)
  
  ZORIG = ASIN(F%RMU(1))
  IF( S%LSHIFTLL )  THEN
    ZINC = ZPI/REAL(INN,JPRD)
    ZOFF = 0.5_JPRD*ZINC
    ZTEMP = ZOFF + ZINC*REAL(S%NDGNHD-1,JPRD)
    ZLRMUZ2(1) =  SIN(MIN(ZTEMP,ZORIG) - 0.5_JPRD*MAX(0._JPRD,ZTEMP - ZORIG))
    ZLRMUZ2(S%NDGNHD) = SIN(ZOFF)
  ELSE
    ZINC = ZPI/REAL(INN-2,JPRD)
    ZOFF=-0.5_JPRD*ZINC
    ZTEMP = ZOFF + ZINC*REAL(S%NDGNHD-1,JPRD)
    ZLRMUZ2(1) =  SIN(MIN(ZTEMP,ZORIG) - 0.5_JPRD*MAX(0._JPRD,ZTEMP - ZORIG))
    ZOFF=0.01_JPRD*ZINC
    ZLRMUZ2(S%NDGNHD) = SIN(ZOFF)
    ZOFF=0._JPRD
  ENDIF
  DO JGL=2, S%NDGNHD-1
    ZLRMUZ2(JGL) = SIN(ZOFF + ZINC*REAL(S%NDGNHD-JGL,JPRD))
  ENDDO
  DO JGL=1, S%NDGNHD
    ISYM = INN-JGL+1
    ZLRMUZ2(ISYM) = -ZLRMUZ2(JGL)
  ENDDO

  IF( LLP2 ) THEN
    WRITE(NOUT,*) 'dual latitudes'
    DO JGL= 1, INN
      WRITE(NOUT,*)  'dual  JGL=',JGL,(180._JPRD/ZPI)*ZINC,(180._JPRD/ZPI)*ASIN(ZLRMUZ2(JGL)),(180._JPRD/ZPI)*ASIN(F%RMU(JGL))
    ENDDO
  ENDIF
  
  ALLOCATE(F%RMU2(INMAX))
  IF (LLP2) WRITE(NOUT,9) 'F%RMU2     ',SIZE(F%RMU2 ),SHAPE(F%RMU2 )
  ALLOCATE(F%RACTHE2(INMAX))
  IF (LLP2) WRITE(NOUT,9) 'F%RACTHE2  ',SIZE(F%RACTHE2),SHAPE(F%RACTHE2 )
  DO JGL=1,INN
    F%RMU2(JGL)    = ZLRMUZ2(JGL)
    F%RACTHE2(JGL) = 1.0_JPRD/(SQRT(1.0_JPRD-ZLRMUZ2(JGL)*ZLRMUZ2(JGL))+ZEPS)/REAL(RA,JPRD)
  ENDDO

  IF (LLP1) WRITE(NOUT,*)  '=== SULEG: Finished dual Gaussian latitudes ==='

  ! inverse + direct map for FMM
  INX=2*R%NDGNH
  INZ=2*S%NDGNHD
  ALLOCATE(S%FMM_INTI)
  CALL SETUP_SEEFMM(INX,F%RMU,INZ,F%RMU2,S%FMM_INTI)

ENDIF

!*       3.2   Computes related arrays

IF(.NOT.D%LGRIDONLY) THEN

  ALLOCATE(S%FA(D%NUMP))

  ALLOCATE(F%R1MU2(R%NDGL))
  IF (LLP2) WRITE(NOUT,9) 'F%R1MU2   ',SIZE(F%R1MU2),SHAPE(F%R1MU2 ) 
  ALLOCATE(F%RACTHE(R%NDGL))
  IF (LLP2) WRITE(NOUT,9) 'F%RACTHE  ',SIZE(F%RACTHE),SHAPE(F%RACTHE ) 

  IF( S%LUSE_BELUSOV) THEN
    ALLOCATE(F%RPNM(R%NLEI3,D%NSPOLEGL))
    IF (LLP2) WRITE(NOUT,9) 'F%RPNM    ',SIZE(F%RPNM),SHAPE(F%RPNM) 
    DO JNM=1,D%NSPOLEGL
      F%RPNM(R%NLEI3,JNM) = 0.0_JPRD
    ENDDO
  ENDIF

!*       3.2   Computes related arrays

  DO JGL=1,R%NDGL
! test cosine differently
    ZTHETA = ASIN(ZLRMUZ(JGL))
    ZCOS   = COS(ZTHETA)
    F%R1MU2(JGL)  = ZCOS**2
    F%RACTHE(JGL) = 1.0_JPRD/ZCOS/REAL(RA,JPRD)
!    F%R1MU2(JGL)  = 1.0_JPRD-ZLRMUZ(JGL)*ZLRMUZ(JGL)
!    F%RACTHE(JGL) = 1.0_JPRD/SQRT(1.0_JPRD-ZLRMUZ(JGL)*ZLRMUZ(JGL))/REAL(RA,JPRD)
  ENDDO

!*       3.3   Working arrays

! compute the Legendre polynomials as a function of the z_k (Gaussian Latitudes)
! this may be faster than calling supolf for each m but uses extra communication
! and the parallelism is more limited ? Nils

  IF( S%LUSE_BELUSOV .AND. .NOT. C%LREAD_LEGPOL  ) THEN

    INSMAX = R%NTMAX+1

    IF( INSMAX /= R%NDGL) THEN
      DEALLOCATE(ZFN)
      ALLOCATE(ZFN(0:INSMAX,0:INSMAX))
      ! Belousov, Swarztrauber use ZFN(0,0)=SQRT(2._JPRD)
      ! IFS normalisation chosen to be 0.5*Integral(Pnm**2) = 1
      ZFN(0,0)=2._JPRD
      DO JN=1,INSMAX
        ZFNN=ZFN(0,0)
        DO JGL=1,JN
          ZFNN=ZFNN*SQRT(1._JPRD-0.25_JPRD/REAL(JGL**2,JPRD))
        ENDDO

        IODD=MOD(JN,2)
        ZFN(JN,JN)=ZFNN
        DO JGL=2,JN-IODD,2
          ZFN(JN,JN-JGL)=ZFN(JN,JN-JGL+2)*REAL((JGL-1)*(2*JN-JGL+2),JPRD)/REAL(JGL*(2*JN-JGL+1),JPRD)
        ENDDO
      ENDDO
    ENDIF

    ALLOCATE(ZLFPOL(0:INSMAX,0:INSMAX))
    ALLOCATE(ZPNMG(R%NSPOLEG))

    DO JH=1,IHEMIS

      IF (JH==1) THEN
        IGL1=D%NLATLS(MYSETW,MYSETV)
        IGL2=D%NLATLE(MYSETW,MYSETV)
      ELSE
        IGL1=R%NDGL-D%NLATLE(MYSETW,MYSETV)+1
        IGL2=R%NDGL-D%NLATLS(MYSETW,MYSETV)+1
      ENDIF

      ILOOP=0
      DO JGL=IGL1,IGL2

        INM = 0
        CALL SUPOL(INSMAX,ZLRMUZ(JGL),ZFN,ZLFPOL)
        DO JM=0,R%NSMAX
          DO JN=INSMAX,JM,-1
            INM = INM+1
            ZPNMG(INM) = ZLFPOL(JM,JN)
          ENDDO
        ENDDO

        CALL GSTATS(1801,2)
        ILOOP = JGL-IGL1+1
        CALL SUTRLE(ZPNMG,JGL,ILOOP)
        CALL GSTATS(1801,3)

      ENDDO

      ILATSMAX=0
      DO JW=1,NPRTRW
        DO JV=1,NPRTRV
          ILATSMAX=MAX(ILATSMAX,D%NLATLE(JW,JV)-D%NLATLS(JW,JV)+1)
        ENDDO
      ENDDO

      ILATS=IGL2-IGL1+1
      IF (S%LSOUTHPNM .AND. IHEMIS==1 .AND. ILATSMAX-1 >= ILATS) THEN
        ! I don't know what to do for south pole. But isn't this piece of code
        ! a dead stuff for poles rows ? 
        CALL ABORT_TRANS('SULEG: WILL BE BROKEN FOR SOUTH HEMISPHERE')
      ENDIF

      DO J=ILATS,ILATSMAX-1
        ILOOP=ILOOP+1
        CALL GSTATS(1801,2)
        CALL SUTRLE(ZPNMG,-1,ILOOP)
        CALL GSTATS(1801,3)
      ENDDO

    ENDDO

    DEALLOCATE(ZLFPOL)
    IF( ALLOCATED(ZFN) ) DEALLOCATE(ZFN)

    DEALLOCATE(ZPNMG)

    IF(LLP1) WRITE(NOUT,*) '=== SULEG: Finished RPNM ==='

  ENDIF

  CALL SETUP_GEOM

  IMAXN=R%NTMAX+1

  ITAG=MTAGLETR
  ITAG1=MTAGLETR+1

  IMAXRECVA=0
  IMAXRECVS=0
  DO JMLOC=1,D%NUMP
    IM = D%MYMS(JMLOC)
    ILA = (R%NSMAX-IM+2)/2
    ILS = (R%NSMAX-IM+3)/2
    IDGLU = MIN(R%NDGNH,G%NDGLU(IM))
    ISL = MAX(R%NDGNH-G%NDGLU(IM)+1,1)
    IMAXRECVA = MAX(IDGLU*ILA,IMAXRECVA)
    IMAXRECVS = MAX(IDGLU*ILS,IMAXRECVS)

    !find nearest starting latitude of the dual set
    IF( S%LDLL ) THEN

      INMAX=MIN(R%NTMAX+1,S%NDGL)
      IDGLU2=S%NDGNHD
      S%FA(JMLOC)%ISLD = 1
      LLA:DO JGL=1,S%NDGNHD-1
        IF( (ZLRMUZ2(JGL) < ZLRMUZ(ISL)) ) THEN
          S%FA(JMLOC)%ISLD = JGL
          IDGLU2 = S%NDGNHD-S%FA(JMLOC)%ISLD+1
          EXIT LLA
        ENDIF
      ENDDO LLA

      IF( .NOT. C%LREAD_LEGPOL  ) THEN
       ! compute auxiliary quantities for the dual mapping

       ! output data latitudes
        ALLOCATE(ZPNMCDO(2*IDGLU2,2))
      !$OMP PARALLEL PRIVATE(JGL,ZLPOL)
        IF (.NOT.ALLOCATED(ZLPOL)) ALLOCATE(ZLPOL(0:INMAX))
      !$OMP DO SCHEDULE(DYNAMIC,1)
        DO JGL=1,2*IDGLU2
          CALL SUPOLF(IM,INMAX,ZLRMUZ2(S%FA(JMLOC)%ISLD+JGL-1),ZLPOL(0:INMAX))
          ZPNMCDO(JGL,1)=ZLPOL(INMAX-1)
          ZPNMCDO(JGL,2)=ZLPOL(INMAX)
        ENDDO
        !$OMP END DO
        IF (ALLOCATED(ZLPOL)) DEALLOCATE(ZLPOL)
        !$OMP END PARALLEL

        ! internal (gg-roots) latitudes
        ALLOCATE(ZPNMCDD(2*IDGLU,2))
        !$OMP PARALLEL PRIVATE(JGL,ZLPOL,JI,JN)
        IF (.NOT.ALLOCATED(ZLPOL)) ALLOCATE(ZLPOL(0:INMAX))
        !$OMP DO SCHEDULE(DYNAMIC,1)
        DO JGL=1,2*IDGLU
          CALL SUPOLF(IM,INMAX,ZLRMUZ(ISL+JGL-1),ZLPOL(0:INMAX))
          ZPNMCDD(JGL,1)=ZLPOL(INMAX-1)
          ZPNMCDD(JGL,2)=ZLPOL(INMAX)
        ENDDO
        !$OMP END DO
        IF (ALLOCATED(ZLPOL)) DEALLOCATE(ZLPOL)
        !$OMP END PARALLEL
        
        CALL PREPSNM(IM,JMLOC,ZEPSNM)
        ALLOCATE(S%FA(JMLOC)%RPNMWI(2*IDGLU,1:2))
        DO JGL=1,2*IDGLU
          ! inverse trafo
          S%FA(JMLOC)%RPNMWI(JGL,1) = F%RW(ISL+JGL-1)*ZPNMCDD(JGL,1)
          S%FA(JMLOC)%RPNMWI(JGL,2) = F%RW(ISL+JGL-1)*ZPNMCDD(JGL,2)
          ! direct trafo needed if mapping to another set of gg roots
          !S%FA(JMLOC)%RPNMWI(JGL,3) = -ZEPSNM(IMAXN)*ZPNMCDD(JGL,2)
          !S%FA(JMLOC)%RPNMWI(JGL,4) = -ZEPSNM(IMAXN)*ZPNMCDD(JGL,1)
        ENDDO
        DEALLOCATE(ZPNMCDD)
        ALLOCATE(S%FA(JMLOC)%RPNMWO(2*IDGLU2,1:2))
        DO JGL=1,2*IDGLU2
          ! inverse trafo
          S%FA(JMLOC)%RPNMWO(JGL,1) = -ZEPSNM(IMAXN)*ZPNMCDO(JGL,2)
          S%FA(JMLOC)%RPNMWO(JGL,2) = -ZEPSNM(IMAXN)*ZPNMCDO(JGL,1)
          ! only needed in direct trafo, need if mapping to another set of roots
          !S%FA(JMLOC)%RPNMWO(JGL,3) = F%RW2(S%FA(JMLOC)%ISLD+JGL-1)*ZPNMCDO(JGL,1)
          !S%FA(JMLOC)%RPNMWO(JGL,4) = F%RW2(S%FA(JMLOC)%ISLD+JGL-1)*ZPNMCDO(JGL,2)
        ENDDO
        DEALLOCATE(ZPNMCDO)
      ENDIF ! LREAD_LEGPOL
    ENDIF ! LDLL

  ENDDO

  IF( S%LDLL ) THEN
    DEALLOCATE(ZLRMUZ2)
  ENDIF

!!$  IF( S%LUSEFLT.AND.LMPOFF ) THEN
!!$    CALL ABORT_TRANS('SULEG: LUSEFLT=T and LMPOFF=T not supported')
!!$  ENDIF
  CALL GSTATS(1801,2)

  IF(.NOT.C%LREAD_LEGPOL) THEN
  IF( S%LUSEFLT )THEN
    ALLOCATE(ZCLONEA(D%NUMP))
    ALLOCATE(ZCLONES(D%NUMP))
  ENDIF

! not correct logic
!  IF (S%LUSE_BELUSOV.OR.S%LUSEFLT) THEN

  DO JMLOC=1,D%NUMP,NPRTRV  ! +++++++++++++++++++++ JMLOC LOOP +++++++++++++++++++++++
!     IF( S%LUSEFLT )THEN
!       ZCLONES(JMLOC)%COMMSBUF => NULL()
!       ZCLONEA(JMLOC)%COMMSBUF => NULL()
!     ENDIF

    IPRTRV=MIN(NPRTRV,D%NUMP-JMLOC+1)

    ! --------------------anti-symmetric-----------------------
    ! --------------------anti-symmetric-----------------------
    ! --------------------anti-symmetric-----------------------

    DO JSETV=1,IPRTRV
      IMLOC=JMLOC+JSETV-1
      IM = D%MYMS(IMLOC)
      ILA = (R%NSMAX-IM+2)/2
      IDGLU = MIN(R%NDGNH,G%NDGLU(IM))
      ALLOCATE(S%FA(IMLOC)%RPNMA(IDGLU,ILA))  
    ENDDO

    IF( .NOT. S%LUSE_BELUSOV ) THEN 

      ISREQ = 0
      IRREQ = 0

      ALLOCATE (ZRCVBUFV(IMAXRECVA,IPRTRV))
      CALL GSTATS(851,0)
      DO JSETV=1,IPRTRV
        CALL SET2PE(IRECV,0,0,MYSETW,JSETV)
        IF( .NOT.LMPOFF )THEN
          IRREQ = IRREQ+1
          CALL MPL_RECV(ZRCVBUFV(:,JSETV),KSOURCE=NPRCIDS(IRECV), &
           & KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=IRECVREQ(IRREQ),&
           & KTAG=ITAG,CDSTRING='SULEG:')
        ENDIF
      ENDDO
      CALL GSTATS(851,1)

      IF( JMLOC+MYSETV-1 <= D%NUMP )THEN

        IMLOC=JMLOC+MYSETV-1
        IM = D%MYMS(IMLOC)
        ISL = MAX(R%NDGNH-G%NDGLU(IM)+1,1)
        IA  = 1+MOD(R%NSMAX-IM+2,2)
        ILA = (R%NSMAX-IM+2)/2
        IDGLU = MIN(R%NDGNH,G%NDGLU(IM))

        ALLOCATE(ZSNDBUFV(IDGLU*ILA))
      
        IF(MOD(IMAXN-IM,2) == 0) THEN
          INMAX=IMAXN+1
        ELSE
          INMAX=IMAXN
        ENDIF

        CALL GSTATS(1251,0)
        IF (.NOT.ALLOCATED(ZLPOL)) ALLOCATE(ZLPOL(0:INMAX))
        !$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JGL,ZLPOL,JI,JN)
        DO JGL=1,IDGLU
          CALL SUPOLF(IM,INMAX,ZLRMUZ(ISL+JGL-1),ZLPOL(0:INMAX),KCHEAP=3)
          DO JI=1,ILA
            JN=IM+2*(JI-1)+1
            ZSNDBUFV((JGL-1)*ILA+JI)=ZLPOL(JN)
          ENDDO
        ENDDO
        !$OMP END PARALLEL DO
        IF (ALLOCATED(ZLPOL)) DEALLOCATE(ZLPOL)
        CALL GSTATS(1251,1)

        CALL GSTATS(851,0)
        DO JSETV=1,NPRTRV
          CALL SET2PE(ISEND,0,0,MYSETW,JSETV)
          IF( .NOT.LMPOFF )THEN
            ISREQ = ISREQ+1
            CALL MPL_SEND(ZSNDBUFV(:),KDEST=NPRCIDS(ISEND), &
             & KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=ISENDREQ(ISREQ),&
            & KTAG=ITAG,CDSTRING='SULEG:')
          ENDIF
        ENDDO
        CALL GSTATS(851,1)

      ENDIF

      CALL GSTATS(851,0)
      IF(IRREQ > 0) THEN
        CALL MPL_WAIT(KREQUEST=IRECVREQ(1:IRREQ), &
         & CDSTRING='SUTRLE: SULEG')
      ENDIF

      IF(ISREQ > 0) THEN
        CALL MPL_WAIT(KREQUEST=ISENDREQ(1:ISREQ), &
         & CDSTRING='SUTRLE: SULEG')
      ENDIF

      IF( NPROC==1.AND.LMPOFF )THEN
        ZRCVBUFV(1:SIZE(ZSNDBUFV(:)),1)=ZSNDBUFV(:)
      ENDIF
      CALL GSTATS(851,1)

      CALL GSTATS(1251,0)
      !$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JSETV,IMLOC,IM,ISL,IA,ILA,IDGLU,JGL,JI)
      DO JSETV=1,IPRTRV
        IMLOC=JMLOC+JSETV-1
        IM = D%MYMS(IMLOC)
        ISL = MAX(R%NDGNH-G%NDGLU(IM)+1,1)
        IA  = 1+MOD(R%NSMAX-IM+2,2)
        ILA = (R%NSMAX-IM+2)/2
        IDGLU = MIN(R%NDGNH,G%NDGLU(IM))
        IF( S%LUSEFLT .AND. ILA > ITHRESHOLD ) THEN
          IF( .NOT. S%LKEEPRPNM ) DEALLOCATE(S%FA(IMLOC)%RPNMA)
          ALLOCATE(S%FA(IMLOC)%RPNMDA(IDGLU,ILA))
          DO JGL=1,IDGLU
            DO JI=1,ILA
              S%FA(IMLOC)%RPNMDA(JGL,ILA-JI+1)=ZRCVBUFV((JGL-1)*ILA+JI,JSETV)
            ENDDO
          ENDDO
          IF( S%LKEEPRPNM ) THEN
             DO JGL=1,IDGLU
                DO JI=1,ILA
                   S%FA(IMLOC)%RPNMA(JGL,ILA-JI+1)=ZRCVBUFV((JGL-1)*ILA+JI,JSETV)
                ENDDO
             ENDDO
          ENDIF
       ELSE
          DO JGL=1,IDGLU
            DO JI=1,ILA
              S%FA(IMLOC)%RPNMA(JGL,ILA-JI+1)=ZRCVBUFV((JGL-1)*ILA+JI,JSETV)
            ENDDO
          ENDDO
        ENDIF
      ENDDO
     !$OMP END PARALLEL DO
      CALL GSTATS(1251,1)
        
      IF( ALLOCATED(ZSNDBUFV) ) DEALLOCATE(ZSNDBUFV)
      IF( ALLOCATED(ZRCVBUFV) ) DEALLOCATE(ZRCVBUFV)

    ELSE    

       CALL GSTATS(1251,0)
       !$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JSETV,IMLOC,IM,ISL,IA,ILA,IDGLU,JGL,JI)
       DO JSETV=1,IPRTRV
          IMLOC=JMLOC+JSETV-1
          IM = D%MYMS(IMLOC)
          ISL = MAX(R%NDGNH-G%NDGLU(IM)+1,1)
          IA  = 1+MOD(R%NSMAX-IM+2,2)
          ILA = (R%NSMAX-IM+2)/2
          IDGLU = MIN(R%NDGNH,G%NDGLU(IM))
          IF( S%LUSEFLT .AND. ILA > ITHRESHOLD ) THEN
             IF( .NOT. S%LKEEPRPNM ) DEALLOCATE(S%FA(IMLOC)%RPNMA)
             ALLOCATE(S%FA(IMLOC)%RPNMDA(IDGLU,ILA))
             DO JI=1,ILA
                DO JGL=1,IDGLU
                   S%FA(IMLOC)%RPNMDA(JGL,JI) = F%RPNM(ISL+JGL-1,D%NPMS(IM)+IA+(JI-1)*2)
                ENDDO
             ENDDO
             IF( S%LKEEPRPNM ) THEN
               DO JI=1,ILA
                 DO JGL=1,IDGLU
                   S%FA(IMLOC)%RPNMA(JGL,JI) = F%RPNM(ISL+JGL-1,D%NPMS(IM)+IA+(JI-1)*2)
                 ENDDO
               ENDDO
            ENDIF
          ELSE
             DO JI=1,ILA
                DO JGL=1,IDGLU
                   S%FA(IMLOC)%RPNMA(JGL,JI) = F%RPNM(ISL+JGL-1,D%NPMS(IM)+IA+(JI-1)*2)
                ENDDO
             ENDDO
          END IF
       ENDDO
       !$OMP END PARALLEL DO
       CALL GSTATS(1251,1)

    ENDIF

    ! -------------------- anti-symmetric FLT iniitialisation -----------------------

    IF( S%LUSEFLT) THEN
      IRECVLENMAX=0
      ISREQ = 0
      IRREQ = 0
      IF( JMLOC+MYSETV-1 <= D%NUMP )THEN
        IMLOC=JMLOC+MYSETV-1
        IM = D%MYMS(IMLOC)
        ILA = (R%NSMAX-IM+2)/2
        IDGLU = MIN(R%NDGNH,G%NDGLU(IM))
        IF( ILA > ITHRESHOLD ) THEN
          S%LSYM = .FALSE.
          INX = IDGLU
          CALL CONSTRUCT_BUTTERFLY(ZEPS_INT_DEC,IMAXCOLS,INX,ILA,S%FA(IMLOC)%RPNMDA,&
           & S%FA(IMLOC)%YBUT_STRUCT_A)
          CALL PACK_BUTTERFLY_STRUCT(S%FA(IMLOC)%YBUT_STRUCT_A,ZCLONEA(IMLOC))
          IRECVLENMAX=SIZE(ZCLONEA(IMLOC)%COMMSBUF)
          CALL GSTATS(852,0)
          DO JSETV=1,NPRTRV
            CALL SET2PE(ISEND,0,0,MYSETW,JSETV)
            IF(.NOT.LMPOFF) THEN
              ISREQ = ISREQ+1
              CALL MPL_SEND(ZCLONEA(IMLOC)%COMMSBUF(:),KDEST=NPRCIDS(ISEND),&
               & KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=ISENDREQ(ISREQ),&
               & KTAG=ITAG,CDSTRING='SULEG:')
            ENDIF
          ENDDO
          CALL GSTATS(852,1)
        ELSE
          IRECVLENMAX=2
          ZDUM(:)=0.0_JPRB
          CALL GSTATS(852,0)
          DO JSETV=1,NPRTRV
            CALL SET2PE(ISEND,0,0,MYSETW,JSETV)
            IF(.NOT.LMPOFF) THEN
              ISREQ = ISREQ+1
              CALL MPL_SEND(ZDUM(:),KDEST=NPRCIDS(ISEND), &
               & KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=ISENDREQ(ISREQ),&
               & KTAG=ITAG,CDSTRING='SULEG:')
            ENDIF
          ENDDO
          CALL GSTATS(852,1)
        ENDIF
      ENDIF
      IF(.NOT.LMPOFF) THEN
        CALL GSTATS(852,0)
        DO JSETV=1,NPRTRV
          CALL SET2PE(ISEND,0,0,MYSETW,JSETV)
          CALL MPL_SEND(IRECVLENMAX,KDEST=NPRCIDS(ISEND),KTAG=ITAG1,CDSTRING='SULEG:')
        ENDDO
        IRECVLENMAX=0
        DO JSETV=1,NPRTRV ! ? should this be IPRTRV ?
          CALL SET2PE(IRECV,0,0,MYSETW,JSETV)
          CALL MPL_RECV(IRECVLENMAXV(JSETV),KSOURCE=NPRCIDS(IRECV),KTAG=ITAG1,CDSTRING='SULEG:')
          IRECVLENMAX=MAX(IRECVLENMAX,IRECVLENMAXV(JSETV))
        ENDDO
        IF( MYPROC == 1 )THEN
          IF( IRECVLENMAX > 2 )THEN
            WRITE(NOUT,'("SULEG: ANTI-SYM MAX BUTTERFLY CLONE LEN=",I8)')IRECVLENMAX
          ENDIF
        ENDIF
        IF( IRECVLENMAX == 0 )THEN
          WRITE(NOUT,'("SULEG: ANTI-SYM WARNING CLONE LEN=",I8,I8)') MYPROC, IRECVLENMAX
        ENDIF
        IF( IRECVLENMAX > 0 )THEN
           ALLOCATE (ZRCVBUTFV(IRECVLENMAX,IPRTRV))
           DO JSETV=1,IPRTRV
              IRREQ = IRREQ+1
              CALL SET2PE(IRECV,0,0,MYSETW,JSETV)
              CALL MPL_RECV(ZRCVBUTFV(:,JSETV),KSOURCE=NPRCIDS(IRECV), &
                   &KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=IRECVREQ(IRREQ),&
                   & KTAG=ITAG,CDSTRING='SULEG:')
           ENDDO
        END IF
        IF(ISREQ > 0) THEN
          CALL MPL_WAIT(KREQUEST=ISENDREQ(1:ISREQ), &
           & CDSTRING='SUTRLE: SULEG')
        ENDIF
        IF(IRREQ > 0) THEN
          CALL MPL_WAIT(KBYTES=IBYTES,KOUNT=IKOUNT(1:IRREQ),KREQUEST=IRECVREQ(1:IRREQ), &
           & CDSTRING='SUTRLE: SULEG')
        ENDIF
        CALL GSTATS(852,1)
        IF( IRECVLENMAX > 0 )THEN
          CALL GSTATS(1252,0)
          !$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JSETV,IMLOC,IM,ILA,IDGLU,INX,ICLONELEN)
          DO JSETV=1,IPRTRV
            IMLOC=JMLOC+JSETV-1
            IM = D%MYMS(IMLOC)
            ILA = (R%NSMAX-IM+2)/2
            IF( ILA > ITHRESHOLD ) THEN
              IDGLU = MIN(R%NDGNH,G%NDGLU(IM))
              INX=IDGLU
              IF( .NOT.ALLOCATED(ZCLONEA(IMLOC)%COMMSBUF) )THEN
                ICLONELEN=IKOUNT(JSETV)
                ALLOCATE(ZCLONEA(IMLOC)%COMMSBUF(ICLONELEN))
                ZCLONEA(IMLOC)%COMMSBUF(1:ICLONELEN) = ZRCVBUTFV(1:ICLONELEN,JSETV)
                CALL UNPACK_BUTTERFLY_STRUCT(S%FA(IMLOC)%YBUT_STRUCT_A,ZCLONEA(IMLOC))
              ENDIF
              IF(ALLOCATED(ZCLONEA(IMLOC)%COMMSBUF) ) THEN
                IF( SIZEOF(ZCLONEA(IMLOC)%COMMSBUF) > 0 ) DEALLOCATE(ZCLONEA(IMLOC)%COMMSBUF)
                !              ZCLONEA(IMLOC)%COMMSBUF=>NULL()
              ENDIF
              IF( ASSOCIATED(S%FA(IMLOC)%RPNMA) .AND. .NOT. S%LKEEPRPNM ) DEALLOCATE(S%FA(IMLOC)%RPNMA)
              IF( ASSOCIATED(S%FA(IMLOC)%RPNMDA) ) DEALLOCATE(S%FA(IMLOC)%RPNMDA)
            ENDIF
          ENDDO
          !$OMP END PARALLEL DO
          CALL GSTATS(1252,1)
          DEALLOCATE(ZRCVBUTFV)
        ENDIF
      ENDIF
   ENDIF

    ! --------------------symmetric-----------------------
    ! --------------------symmetric-----------------------
    ! --------------------symmetric-----------------------

    DO JSETV=1,IPRTRV
      IMLOC=JMLOC+JSETV-1
      IM = D%MYMS(IMLOC)
      ILS = (R%NSMAX-IM+3)/2
      IDGLU = MIN(R%NDGNH,G%NDGLU(IM))
      ALLOCATE(S%FA(IMLOC)%RPNMS(IDGLU,ILS))  
    ENDDO

    IF( .NOT. S%LUSE_BELUSOV ) THEN 

      ISREQ = 0
      IRREQ = 0

      ALLOCATE (ZRCVBUFV(IMAXRECVS,IPRTRV))
      CALL GSTATS(851,0)
      DO JSETV=1,IPRTRV
        CALL SET2PE(IRECV,0,0,MYSETW,JSETV)
        IF( .NOT.LMPOFF )THEN
          IRREQ = IRREQ+1
          CALL MPL_RECV(ZRCVBUFV(:,JSETV),KSOURCE=NPRCIDS(IRECV), &
           & KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=IRECVREQ(IRREQ),&
           & KTAG=ITAG,CDSTRING='SULEG:')
        ENDIF
      ENDDO
      CALL GSTATS(851,1)

      IF( JMLOC+MYSETV-1 <= D%NUMP )THEN

        IMLOC=JMLOC+MYSETV-1
        IM = D%MYMS(IMLOC)
        ISL = MAX(R%NDGNH-G%NDGLU(IM)+1,1)
        IS  = 1+MOD(R%NSMAX-IM+1,2)
        ILS = (R%NSMAX-IM+3)/2
        IDGLU = MIN(R%NDGNH,G%NDGLU(IM))

        ALLOCATE(ZSNDBUFV(IDGLU*ILS))
      
        IF(MOD(IMAXN-IM,2) == 0) THEN
          INMAX=IMAXN
        ELSE
          INMAX=IMAXN+1
        ENDIF

        IF (.NOT.ALLOCATED(ZLPOL)) ALLOCATE(ZLPOL(0:INMAX))
        CALL GSTATS(1251,0)
        !$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JGL,ZLPOL,JI,JN)
        DO JGL=1,IDGLU
          CALL SUPOLF(IM,INMAX,ZLRMUZ(ISL+JGL-1),ZLPOL(0:INMAX),KCHEAP=2)
          DO JI=1,ILS
            JN=IM+2*(JI-1)
            ZSNDBUFV((JGL-1)*ILS+JI)=ZLPOL(JN)
          ENDDO
        ENDDO
        !$OMP END PARALLEL DO
        CALL GSTATS(1251,1)
        IF (ALLOCATED(ZLPOL)) DEALLOCATE(ZLPOL)

        CALL GSTATS(851,0)
        DO JSETV=1,NPRTRV
          CALL SET2PE(ISEND,0,0,MYSETW,JSETV)
          IF( .NOT.LMPOFF )THEN
            ISREQ = ISREQ+1
            CALL MPL_SEND(ZSNDBUFV(:),KDEST=NPRCIDS(ISEND), &
             & KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=ISENDREQ(ISREQ),&
             & KTAG=ITAG,CDSTRING='SULEG:')
          ENDIF
        ENDDO
        CALL GSTATS(851,1)

      ENDIF

      CALL GSTATS(851,0)
      IF(IRREQ > 0) THEN
        CALL MPL_WAIT(KREQUEST=IRECVREQ(1:IRREQ), &
         & CDSTRING='SUTRLE: SULEG')
      ENDIF

      IF(ISREQ > 0) THEN
        CALL MPL_WAIT(KREQUEST=ISENDREQ(1:ISREQ), &
         & CDSTRING='SUTRLE: SULEG')
      ENDIF
      IF( NPROC==1.AND.LMPOFF )THEN
        ZRCVBUFV(1:SIZE(ZSNDBUFV(:)),1)=ZSNDBUFV(:)
      ENDIF
      CALL GSTATS(851,1)

      CALL GSTATS(1251,0)
      !$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JSETV,IMLOC,IM,ISL,IS,ILS,IDGLU,JGL,JI)
      DO JSETV=1,IPRTRV
        IMLOC=JMLOC+JSETV-1
        IM = D%MYMS(IMLOC)
        ISL = MAX(R%NDGNH-G%NDGLU(IM)+1,1)
        IS  = 1+MOD(R%NSMAX-IM+1,2)
        ILS = (R%NSMAX-IM+3)/2
        IDGLU = MIN(R%NDGNH,G%NDGLU(IM))
        IF( S%LUSEFLT .AND. ILS > ITHRESHOLD ) THEN
          IF( .NOT. S%LKEEPRPNM ) DEALLOCATE(S%FA(IMLOC)%RPNMS)
          ALLOCATE(S%FA(IMLOC)%RPNMDS(IDGLU,ILS))
          DO JGL=1,IDGLU
            DO JI=1,ILS
              S%FA(IMLOC)%RPNMDS(JGL,ILS-JI+1)=ZRCVBUFV((JGL-1)*ILS+JI,JSETV)
            ENDDO
          ENDDO
          IF( S%LKEEPRPNM ) THEN
            DO JGL=1,IDGLU
              DO JI=1,ILS
                S%FA(IMLOC)%RPNMS(JGL,ILS-JI+1)=ZRCVBUFV((JGL-1)*ILS+JI,JSETV)
              ENDDO
            ENDDO
          ENDIF
        ELSE
          DO JGL=1,IDGLU
            DO JI=1,ILS
              S%FA(IMLOC)%RPNMS(JGL,ILS-JI+1)=ZRCVBUFV((JGL-1)*ILS+JI,JSETV)
            ENDDO
          ENDDO
       ENDIF
    ENDDO
    !$OMP END PARALLEL DO
      CALL GSTATS(1251,1)
        
      IF( ALLOCATED(ZSNDBUFV) ) DEALLOCATE(ZSNDBUFV)
      IF( ALLOCATED(ZRCVBUFV) ) DEALLOCATE(ZRCVBUFV)

    ELSE    
      CALL GSTATS(1251,0)
      !$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JSETV,IMLOC,IM,ISL,IS,ILS,IDGLU,JGL,JI)
      DO JSETV=1,IPRTRV
        IMLOC=JMLOC+JSETV-1
        IM = D%MYMS(IMLOC)
        ISL = MAX(R%NDGNH-G%NDGLU(IM)+1,1)
        IS  = 1+MOD(R%NSMAX-IM+1,2)
        ILS = (R%NSMAX-IM+3)/2
        IDGLU = MIN(R%NDGNH,G%NDGLU(IM))
        IF( S%LUSEFLT .AND. ILS > ITHRESHOLD ) THEN
           IF( .NOT. S%LKEEPRPNM ) DEALLOCATE(S%FA(IMLOC)%RPNMS)
           ALLOCATE(S%FA(IMLOC)%RPNMDS(IDGLU,ILS))
           DO JI=1,ILS
              DO JGL=1,IDGLU
                 S%FA(IMLOC)%RPNMDS(JGL,JI) = F%RPNM(ISL+JGL-1,D%NPMS(IM)+IS+(JI-1)*2)
              ENDDO
           ENDDO
           IF( S%LKEEPRPNM ) THEN
             DO JI=1,ILS
               DO JGL=1,IDGLU
                 S%FA(IMLOC)%RPNMS(JGL,JI) = F%RPNM(ISL+JGL-1,D%NPMS(IM)+IS+(JI-1)*2)
               ENDDO
             ENDDO
           ENDIF
         ELSE
           DO JI=1,ILS
              DO JGL=1,IDGLU
                 S%FA(IMLOC)%RPNMS(JGL,JI) = F%RPNM(ISL+JGL-1,D%NPMS(IM)+IS+(JI-1)*2)
              ENDDO
           ENDDO
        END IF
     END DO
     !$OMP END PARALLEL DO
     CALL GSTATS(1251,1)
        
  ENDIF

    ! -------------------- symmetric FLT iniitialisation -----------------------

    IF( S%LUSEFLT) THEN
      IRECVLENMAX=0
      ISREQ = 0
      IRREQ = 0
      IF( JMLOC+MYSETV-1 <= D%NUMP )THEN
         IMLOC=JMLOC+MYSETV-1
         IM = D%MYMS(IMLOC)
         ILS = (R%NSMAX-IM+3)/2
         IDGLU = MIN(R%NDGNH,G%NDGLU(IM))
         IF( ILS > ITHRESHOLD ) THEN
            S%LSYM = .TRUE.
            INX = IDGLU
            CALL CONSTRUCT_BUTTERFLY(ZEPS_INT_DEC,IMAXCOLS,INX,ILS,S%FA(IMLOC)%RPNMDS,&
                 & S%FA(IMLOC)%YBUT_STRUCT_S)
          CALL PACK_BUTTERFLY_STRUCT(S%FA(IMLOC)%YBUT_STRUCT_S,ZCLONES(IMLOC))
          IRECVLENMAX=SIZE(ZCLONES(IMLOC)%COMMSBUF)
          CALL GSTATS(852,0)
          DO JSETV=1,NPRTRV
            CALL SET2PE(ISEND,0,0,MYSETW,JSETV)
            IF(.NOT.LMPOFF) THEN
              ISREQ = ISREQ+1
              CALL MPL_SEND(ZCLONES(IMLOC)%COMMSBUF(:),KDEST=NPRCIDS(ISEND),&
               &KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=ISENDREQ(ISREQ),&
               & KTAG=ITAG,CDSTRING='SULEG:')
            ENDIF
          ENDDO
          CALL GSTATS(852,1)
        ELSE
          IRECVLENMAX=2
          ZDUM(:)=0.0_JPRB
          CALL GSTATS(852,0)
          DO JSETV=1,NPRTRV
            CALL SET2PE(ISEND,0,0,MYSETW,JSETV)
            IF(.NOT.LMPOFF) THEN
              ISREQ = ISREQ+1
              CALL MPL_SEND(ZDUM(:),KDEST=NPRCIDS(ISEND), &
               &KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=ISENDREQ(ISREQ),&
               & KTAG=ITAG,CDSTRING='SULEG:')
            ENDIF
          ENDDO
          CALL GSTATS(852,1)
        ENDIF
      ENDIF
      IF(.NOT. LMPOFF) THEN
        CALL GSTATS(852,0)
        DO JSETV=1,NPRTRV
          CALL SET2PE(ISEND,0,0,MYSETW,JSETV)
          CALL MPL_SEND(IRECVLENMAX,KDEST=NPRCIDS(ISEND),KTAG=ITAG1,CDSTRING='SULEG:')
        ENDDO
        IRECVLENMAX=0
        DO JSETV=1,NPRTRV ! ? should this be IPRTRV ?
          CALL SET2PE(IRECV,0,0,MYSETW,JSETV)
          CALL MPL_RECV(IRECVLENMAXV(JSETV),KSOURCE=NPRCIDS(IRECV),KTAG=ITAG1,CDSTRING='SULEG:')
          IRECVLENMAX=MAX(IRECVLENMAX,IRECVLENMAXV(JSETV))
        ENDDO
        IF( MYPROC == 1 )THEN
          IF( IRECVLENMAX > 2 )THEN
            WRITE(NOUT,'("SULEG: SYM      MAX BUTTERFLY CLONE LEN=",I8)')IRECVLENMAX
          ENDIF
        ENDIF
        IF( IRECVLENMAX == 0 )THEN
          WRITE(NOUT,'("SULEG: SYM WARNING CLONE LEN=",I8,I8)')MYPROC, IRECVLENMAX
        ENDIF
        IF( IRECVLENMAX > 0 )THEN
          ALLOCATE (ZRCVBUTFV(IRECVLENMAX,IPRTRV))
          DO JSETV=1,IPRTRV
            IRREQ = IRREQ+1
            CALL SET2PE(IRECV,0,0,MYSETW,JSETV)
            CALL MPL_RECV(ZRCVBUTFV(:,JSETV),KSOURCE=NPRCIDS(IRECV), &
             &KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=IRECVREQ(IRREQ),&
             & KTAG=ITAG,CDSTRING='SULEG:')
          ENDDO
        ENDIF
        IF(ISREQ > 0) THEN
          CALL MPL_WAIT(KREQUEST=ISENDREQ(1:ISREQ), &
           & CDSTRING='SUTRLE: SULEG')
        ENDIF
        IF(IRREQ > 0) THEN
          CALL MPL_WAIT(KBYTES=IBYTES,KOUNT=IKOUNT(1:IRREQ),KREQUEST=IRECVREQ(1:IRREQ), &
           & CDSTRING='SUTRLE: SULEG')
        ENDIF
        CALL GSTATS(852,1)


        IF( IRECVLENMAX > 0 )THEN
          CALL GSTATS(1252,0)
          !$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JSETV,IMLOC,IM,ILS,IDGLU,INX,ICLONELEN)
          DO JSETV=1,IPRTRV
            IMLOC=JMLOC+JSETV-1
            IM = D%MYMS(IMLOC)
            ILS = (R%NSMAX-IM+3)/2
            IF( ILS > ITHRESHOLD ) THEN
              IDGLU = MIN(R%NDGNH,G%NDGLU(IM))
              INX=IDGLU
              IF( .NOT.ALLOCATED(ZCLONES(IMLOC)%COMMSBUF) )THEN
                ICLONELEN=IKOUNT(JSETV)
                ALLOCATE(ZCLONES(IMLOC)%COMMSBUF(ICLONELEN))
                ZCLONES(IMLOC)%COMMSBUF(1:ICLONELEN) = ZRCVBUTFV(1:ICLONELEN,JSETV)
                CALL UNPACK_BUTTERFLY_STRUCT( S%FA(IMLOC)%YBUT_STRUCT_S,ZCLONES(IMLOC))
              ENDIF
              IF( ALLOCATED(ZCLONES(IMLOC)%COMMSBUF) ) THEN
                IF( SIZEOF(ZCLONES(IMLOC)%COMMSBUF) > 0 ) DEALLOCATE(ZCLONES(IMLOC)%COMMSBUF)
                !              ZCLONES(IMLOC)%COMMSBUF=>NULL()
              ENDIF
              IF( ASSOCIATED(S%FA(IMLOC)%RPNMS) .AND. .NOT. S%LKEEPRPNM ) DEALLOCATE(S%FA(IMLOC)%RPNMS)
              IF( ASSOCIATED(S%FA(IMLOC)%RPNMDS) ) DEALLOCATE(S%FA(IMLOC)%RPNMDS)
            ENDIF
          ENDDO
          !$OMP END PARALLEL DO
          CALL GSTATS(1252,1)
          DEALLOCATE(ZRCVBUTFV)
        ENDIF
      ENDIF
    ENDIF

  ENDDO                     ! +++++++++++++++++++++ END JMLOC LOOP +++++++++++++++++++++++

!  ENDIF ! (S%LUSE_BELUSOV.OR.S%LUSEFLT)
    
  IF( S%LUSEFLT )THEN
    DEALLOCATE(ZCLONEA)
    DEALLOCATE(ZCLONES)
  ENDIF

  IF( LLP1 .AND. S%LUSEFLT ) THEN    
    WRITE(NOUT,*) '=== SULEG: Finished SETUP_BUTTERFLY ==='
  ENDIF
  ENDIF

  CALL GSTATS(1801,3)
  IF(S%LUSE_BELUSOV) DEALLOCATE(F%RPNM)

  IF(C%LWRITE_LEGPOL) CALL WRITE_LEGPOL
  IF(C%LREAD_LEGPOL)  CALL READ_LEGPOL


ENDIF
CALL GSTATS(1801,1)
CALL GSTATS(140,1)

!     ------------------------------------------------------------------
9 FORMAT(1X,'ARRAY ',A10,' ALLOCATED ',8I8)

!IF (S%LUSE_BELUSOV.OR.S%LUSEFLT) 
CALL END_POL

END SUBROUTINE SULEG
END MODULE SULEG_MOD
