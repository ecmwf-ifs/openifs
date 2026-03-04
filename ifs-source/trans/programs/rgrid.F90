! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

PROGRAM RGRID

! Purpose :
! -------
!   A calling program to make a gaussian reduced grid from several parameters

! Interface :
! ---------
!   None

! Externals :
! ---------
!   REDUCED_GRID  - to compute the reduced grid

! Method :
! ------

! Reference :
! ---------

! Author :
! ------
!   12-Feb-2007 R. El Khatib  *METEO-FRANCE*

! Modifications :
! -------------
!   08-Nov 2007 R. El Khatib : Fix arguments to get_opt

! End Modifications
!-----------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB, JPRD
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE MPL_MODULE


IMPLICIT NONE


INTEGER(KIND=JPIM) :: IULOUT, IULERR ! stdout, stderr unit numbers
INTEGER(KIND=JPIM) :: IULNAM ! unit number of output file containing NAMRGRI
INTEGER(KIND=JPIM) :: IDGLG ! number of latidudes
INTEGER(KIND=JPIM) :: IDLON ! number of longitudes
INTEGER(KIND=JPIM) :: ISMAXG ! nominal truncation
INTEGER(KIND=JPIM) :: IXMAX  ! truncation used in the spectral transforms,quad grid
INTEGER(KIND=JPIM) :: INMAX  ! truncation used in the spectral transforms, lin grid
INTEGER(KIND=JPIM) :: IALIAS ! allowed aliasing, as Log10() value
INTEGER(KIND=JPIM) :: IORTHO ! orthogonality precision, as Log10() value
INTEGER(KIND=JPIM) :: IODD ! Odd numbers allowed (1) or not (0)
INTEGER(KIND=JPIM) :: IVERBOSE ! Verbosity level (0 or 1)
INTEGER(KIND=JPIM) :: IPROC ! Number of MPI tasks
INTEGER(KIND=JPIM) :: IMYPROC ! current MPI task

CHARACTER(LEN=64) :: CLARG
CHARACTER(LEN=32) :: CLENV

INTEGER(KIND=JPIM) :: ICOMMAND
LOGICAL :: LLMPI=.TRUE. ! MPI is initialized
! LLREAD and LLWRITE can be used for very high resolutions (>T10K)
! to save on expensive recomputation of ZGR when simply experimenting with a different 
! ALIAS value
LOGICAL :: LLREAD=.FALSE. ! ZGR not read from file per task
LOGICAL :: LLWRITE=.FALSE.! ZGR not written to file per task

REAL(KIND=JPRB) :: ZALIAS, ZORTHO

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------------

! Initialize message passing
CALL MPL_INIT(LDINFO=.FALSE.) ! Do not produce any output

! First call to DR_HOOK must be after MPI_INIT
IF (LHOOK) CALL DR_HOOK('RGRID',0,ZHOOK_HANDLE)

IPROC   = MPL_NPROC()
IMYPROC = MPL_MYRANK()

! Default parameters

IULOUT=6
IULNAM=6
IULERR=0
IDLON=64
IDGLG=HUGE(IDLON)
ISMAXG=21
IALIAS=4
!IORTHO=12
IORTHO=HUGE(IALIAS)
IODD=1
IVERBOSE=0

! Crack options

ICOMMAND=1
DO
  CALL GET_COMMAND_ARGUMENT(ICOMMAND,CLARG)
  IF (LEN_TRIM(CLARG) == 0) EXIT
  IF (CLARG(1:2) == '-o') THEN
!   -o bsolute value of orthogonality exponent threshold
    IF (LEN_TRIM(CLARG) == 2) THEN
      ICOMMAND=ICOMMAND+1
      CALL GET_COMMAND_ARGUMENT(ICOMMAND,CLARG)
      READ(UNIT=CLARG,FMT='(I2)') IORTHO
    ELSE
      READ(UNIT=CLARG(3:),FMT='(I2)') IORTHO
    ENDIF
  ELSEIF (CLARG(1:2) == '-a') THEN
!   -a absolute value of aliasing exponent
    IF (LEN_TRIM(CLARG) == 2) THEN
      ICOMMAND=ICOMMAND+1
      CALL GET_COMMAND_ARGUMENT(ICOMMAND,CLARG)
      READ(UNIT=CLARG,FMT='(I2)') IALIAS
    ELSE
      READ(UNIT=CLARG(3:),FMT='(I2)') IALIAS
    ENDIF
  ELSEIF (CLARG(1:2) == '-n') THEN
!   - numbers allowed (odd or even)
    IF (LEN_TRIM(CLARG) == 2) THEN
      ICOMMAND=ICOMMAND+1
      CALL GET_COMMAND_ARGUMENT(ICOMMAND,CLARG)
      READ(UNIT=CLARG,FMT='(I2)') IODD 
    ELSE
      READ(UNIT=CLARG(3:),FMT='(I2)') IODD
    ENDIF
    IODD=MAX(0,MIN(1,IODD))
  ELSEIF (CLARG(1:2) == '-t') THEN
!   -t nominal truncation
    IF (LEN_TRIM(CLARG) == 2) THEN
      ICOMMAND=ICOMMAND+1
      CALL GET_COMMAND_ARGUMENT(ICOMMAND,CLARG)
      READ(UNIT=CLARG,FMT='(I5)') ISMAXG 
    ELSE
      READ(UNIT=CLARG(3:),FMT='(I5)') ISMAXG
    ENDIF
  ELSEIF (CLARG(1:2) == '-l') THEN
!   -l number of longitudes
    IF (LEN_TRIM(CLARG) == 2) THEN
      ICOMMAND=ICOMMAND+1
      CALL GET_COMMAND_ARGUMENT(ICOMMAND,CLARG)
      READ(UNIT=CLARG,FMT='(I5)') IDLON 
    ELSE
      READ(UNIT=CLARG(3:),FMT='(I5)') IDLON
    ENDIF
  ELSEIF (CLARG(1:2) == '-g') THEN
!   -g number of Gaussian latitudes
    IF (LEN_TRIM(CLARG) == 2) THEN
      ICOMMAND=ICOMMAND+1
      CALL GET_COMMAND_ARGUMENT(ICOMMAND,CLARG)
      READ(UNIT=CLARG,FMT='(I5)') IDGLG
    ELSE
      READ(UNIT=CLARG(3:),FMT='(I5)') IDGLG
    ENDIF
  ELSEIF (CLARG(1:2) == '-v') THEN
!   -v verbosity
    IF (LEN_TRIM(CLARG) == 2) THEN
      ICOMMAND=ICOMMAND+1
      CALL GET_COMMAND_ARGUMENT(ICOMMAND,CLARG)
      READ(UNIT=CLARG,FMT='(I4)') IVERBOSE
    ELSE
      READ(UNIT=CLARG(3:),FMT='(I2)') IVERBOSE
    ENDIF
    IVERBOSE=MAX(0,MIN(1,IVERBOSE))
  ELSEIF (CLARG(1:2) == '-f') THEN
!   -f unit number of namelist file NAMRGRI
    IF (LEN_TRIM(CLARG) == 2) THEN
      ICOMMAND=ICOMMAND+1
      CALL GET_COMMAND_ARGUMENT(ICOMMAND,CLARG)
      READ(UNIT=CLARG,FMT='(I4)') IULNAM
    ELSE
      READ(UNIT=CLARG(3:),FMT='(I4)') IULNAM
    ENDIF
    IULNAM=MAX(0,IULNAM)
  ELSEIF (CLARG(1:2) == '-r') THEN
!   -r read ZGR from file (i.e. do not compute ZGR)
    LLREAD=.TRUE.
  ELSEIF (CLARG(1:2) == '-w') THEN
!   -w write ZGR to file (i.e. so that we can read it on a subsequent run)
    LLWRITE=.TRUE.
  ELSE
    IF (IMYPROC == 1) THEN
      PRINT*, ' USAGE:'
      PRINT*, ' -t nominal truncation                        [',ISMAXG,']'
      PRINT*, ' -l number of longitudes                      [',IDLON,']'
      PRINT*, ' -g number of Gaussian latitudes              (self-determined if not specified)'
      PRINT*, ' -o orthogonality precision, as Log10() value (self-determined if not specified)'
      PRINT*, ' -a allowed aliasing, as a Log10() value      [',IALIAS,']'
      PRINT*, ' -n odd numbers allowed (1) or not (0)        [',IODD,']'
      PRINT*, ' -v verbosity (0 or 1)                        [',IVERBOSE,']'
      PRINT*, ' -f unit number of namelist file NAMRGRI      [',IULNAM,']'
      PRINT*, ' -r read ZGR from file per task               [',LLREAD,']'
      PRINT*, ' -w write ZGR to file per task                [',LLWRITE,']'
      PRINT*, ' -h displays options and arguments'
    ENDIF
    IF (CLARG(1:2) == '-h') THEN
      CALL MPL_END
      IF (LHOOK) CALL DR_HOOK('RGRID',1,ZHOOK_HANDLE)
      STOP
    ELSE
      CALL ABOR1('RGRID : ERROR IN COMMAND LINE ARGUMENTS')
    ENDIF
  ENDIF
  ICOMMAND=ICOMMAND+1
ENDDO

CALL MPL_BARRIER(CDSTRING='RGRID:')

IF (LLREAD.AND.LLWRITE) CALL ABOR1('RGRID : LLREAD AND LLWRITE ARE EXCLUSIVE')

ZALIAS=10._JPRB**(-REAL(IALIAS,KIND=JPRB))
IF (IORTHO==HUGE(IALIAS)) THEN
  ZORTHO=1000._JPRB*EPSILON(ZALIAS)
ELSE
  ZORTHO=10._JPRB**(-REAL(IORTHO,KIND=JPRB))
ENDIF

IF (IDGLG==HUGE(IDLON)) THEN
  IDGLG=(IDLON+1)/2
ENDIF

IF (ISMAXG > (IDLON+3)/3) THEN
! Linear grid : reproduce ECMWF code :
  IXMAX=ISMAXG
  INMAX=MIN((IDLON-1)/3,ISMAXG)
ELSE
! Quadratic grid : reproduce Meteo-France code :
  IXMAX=2*ISMAXG
  INMAX=MIN((IDLON-1)/3,IXMAX)
ENDIF

! Computation

CALL REDUCED_GRID(IULOUT,IULERR,IDGLG,IDLON,ISMAXG,IXMAX,INMAX,ZORTHO,ZALIAS, &
 & IODD,IPROC,IMYPROC,IVERBOSE,IULNAM,LLREAD,LLWRITE)

! Finalize

CALL MPL_END()

IF (LHOOK) CALL DR_HOOK('RGRID',1,ZHOOK_HANDLE)

STOP

CONTAINS

SUBROUTINE REDUCED_GRID(KULOUT,KULERR,KDGLG,KDLON,KSMAXG,KXMAX,KNMAX,PORTHO,PALIAS, &
 & KODD,KPROC,KMYPROC,KVERBOSE,KULNAM,LDREAD,LDWRITE)

! Purpose :
! -------
!   *REDUCED_GRID* Compute a global reduced grid, from a given truncation and
!   maximum aliasing rate.

! Interface :
! ---------
!   KULOUT : standard output unit number
!   KULERR : standard error unit number
!   KDGLG  : total number of latitude
!   KDLON  : maximum number of longitude
!   KSMAXG : spectral truncation
!   PORTHO : Orthogonality threshold (abs. exponent)
!   PALIAS : maximum allowed aliasing rate (exponent)
!   KPROC  : number of MPI tasks
!   KMYPROC : current MPI task
!   KVERBOSE : verbosity option
!   KULNAM : unit number for writing rtable file
!   LDREAD : F - (default) compute ZGR, 
!          : T - do not compute ZGR, read from save file per task
!   LDWRITE: F - (default) do not write ZGR, 
!          : T - write ZGR to a separate file per task

! Externals :
! ---------
!   SETUP_TRANS0 - basic initialization
!   SETUP_TRANS  - resolution dependent initialization

! Method :
! ------

! Use the spectral transform package to compute the Legendre polynomials at
! twice the given truncation.
! Check orthogonality of polynomials up to IOMAX with those within KNMAX
! IOMAX=KXMAX ==> quadratic grid if KXMAX=KNMAX*1.5 (oversized KXMAX)
! IOMAX=KNMAX ==> linear grid
! No need for quadruple precision from cycle 37 onwards, thanks to the Swartztrauber algorithm

! Recommended values :
!   For orthogonality : 
!     1.E-12 
!     but the best (and default) value is now 1000*EPSILON(Z) from cycle 37 onwards.
!   For Aliasing : 
!     1.E-4 at Meteo-France
!     1.E-2 at ECMWF, except at T7999 (1.E-1)

! Notice : Fast Legendre Transforms not working yet.

! Reference :
! ---------
!   Courtier and Naughton (1994)  *ECMWF/METEO-FRANCE*

! Author :
! ------
!   12-Feb-2007 R. El Khatib  *METEO-FRANCE*

! Modifications :
! -------------
!   R. El Khatib 30-Nov-2012 Optimization, OPEN-MP parallelization, cleanings.
!   R. El Khatib 29-Apr-2013 More cleanings and distribution, merge with ECMWF
!   G. Mozdzynski Jan 2014   Optimization for reducing memory requirement
!   G. Mozdzynski Jan 2014   Options -w and -r to tune the orthogonality
!                            threshold without recomputing ZGR
!   R. El Khatib 15-Feb-2016 Merge and finalization of G. Mozdzynski's mods

! End Modifications
!-----------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN)  :: KULOUT, KULERR, KDGLG, KSMAXG, KXMAX, KNMAX, KDLON
REAL(KIND=JPRB),    INTENT(IN)  :: PORTHO, PALIAS
INTEGER(KIND=JPIM), INTENT(IN)  :: KODD, KPROC, KMYPROC, KVERBOSE, KULNAM
LOGICAL,            INTENT(IN)  :: LDREAD, LDWRITE

INTEGER(KIND=JPIM) :: IOMAX, IDGNH, INUMP, IM, IPRTRV
INTEGER(KIND=JPIM) :: JOR, JM, JN1, JN2, JGL, J, IGPTOT, IGPTOTS
INTEGER(KIND=JPIM) :: IPRGPNS, IPRGPEW, IPRTRW, IPROC, ITAG1, ITAG2, JROC, IMSGLEN

REAL(KIND=JPRB) :: ZMAX(KDGLG)
REAL(KIND=JPRB) :: ZGW(KDGLG) ! Gaussian weights & sines of latitudes
REAL(KIND=JPRD) :: ZMU(KDGLG) 
REAL(KIND=JPRB) :: Z2GW(KDGLG) ! Gaussian weights *2
REAL(KIND=JPRB) :: ZUNITY(0:MAX(KXMAX,KNMAX)+MOD(MAX(KXMAX,KNMAX)+1,2))
REAL(KIND=JPRB) :: ZORTH(0:MAX(KXMAX,KNMAX)+MOD(MAX(KXMAX,KNMAX)+1,2))

REAL(KIND=JPRB), ALLOCATABLE :: ZPNMT(:,:)  ! transposed Legendre polynomials
REAL(KIND=JPRB), ALLOCATABLE :: ZGR(:,:,:)

INTEGER(KIND=JPIM) :: IFFTP0(KDLON), IPRCIDS(KPROC)
INTEGER(KIND=JPIM), ALLOCATABLE :: ILOENG(:) ! number of longitudes for each latitude (global)
INTEGER(KIND=JPIM), ALLOCATABLE :: ILOENS(:)
INTEGER(KIND=JPIM), ALLOCATABLE :: IILOEN(:) ! number of longitudes for each latitude (local contrib)
INTEGER(KIND=JPIM), ALLOCATABLE :: IMYMS(:) ! wave numbers on local task

INTEGER(KIND=JPIM) :: ILEI3    ! First dimension of Legendre polynomials
INTEGER(KIND=JPIM) :: ISPOLEG  ! Second dimension of Legendre polynomials
INTEGER(KIND=JPIM) :: IULZGR   ! Logical unit number for ZGR
 
CHARACTER(LEN=5) :: CLMODE
CHARACTER(LEN=13) :: CFILE

LOGICAL :: LLODD, LLINEAR

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE0
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE1
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE2


#include "setup_trans0.h"
#include "setup_trans.h"
#include "trans_inq.h"
#include "trans_pnm.h"

#include "sufftp.h"

!-----------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('RGRID:REDUCED_GRID',0,ZHOOK_HANDLE0)

! 1. Setup
!    -----

IDGNH=(KDGLG+1)/2
LLINEAR=(KXMAX==KSMAXG)

! Provisional "A distribution" only :
IPROC=KPROC
IPRTRV=1
IPRGPNS=IPROC/IPRTRV
IPRGPEW=IPRTRV
IPRTRW=IPROC/IPRTRV

ITAG1=2000
ITAG2=2001
DO JROC=1,IPRTRW
  IPRCIDS(JROC)=JROC
ENDDO

IULZGR=44
IF (LDREAD .OR. LDWRITE) WRITE(CFILE,'("ZGR.TASK",I5.5)') KMYPROC
IF (LDWRITE) OPEN(UNIT=IULZGR,FILE=CFILE,FORM='UNFORMATTED',STATUS='NEW')
IF (LDREAD) OPEN(IULZGR,FILE=CFILE,FORM='UNFORMATTED',STATUS='OLD')

! 2. Compute Legendre polynomials
!    ----------------------------

CALL SETUP_TRANS0(KOUT=KULOUT,KERR=KULERR,LDMPOFF=.FALSE., &
 & KPRGPNS=IPRGPNS,KPRGPEW=IPRTRV,KPRTRW=IPRTRW)
CALL SETUP_TRANS(KSMAX=KXMAX,KDGL=KDGLG,LDUSEFLT=.FALSE.,LDUSERPNM=.FALSE., &
 & LDKEEPRPNM=.FALSE.,LDPNMONLY=.TRUE.)

CALL TRANS_INQ(KLEI3=ILEI3,KSPOLEGL=ISPOLEG,KNUMP=INUMP)

IF (KMYPROC == 1) THEN
  WRITE(KULOUT,*) 'dimensions: IPROC, ILEI3, ISPOLEG, INUMP, KXMAX, KNMAX, IDGNH '
  WRITE(KULOUT,*) IPROC, ILEI3, ISPOLEG, INUMP, KXMAX, KNMAX, IDGNH
ENDIF

ALLOCATE(IMYMS(INUMP)) 
CALL TRANS_INQ(KMYMS=IMYMS)
ALLOCATE(ZGR(INUMP,IDGNH,2))

IF( .NOT. LDREAD ) THEN

  CALL TRANS_INQ(PGW=ZGW,PMU=ZMU)
  DO JGL=1,IDGNH
    Z2GW(JGL)=2._JPRB*ZGW(JGL)
  ENDDO
  IF (KVERBOSE == 1 .AND. KMYPROC == 1) THEN
    WRITE(KULERR,'(A,E10.3)') ' Error in the sum of weights:',0.5_JPRB-SUM(ZGW(1:IDGNH))
    DO JGL=IDGNH,1,-1
      WRITE(KULERR,*) ' Gaussian weight: ', JGL, ZGW(JGL)
      WRITE(KULERR,*) ' Gaussian latitude (degrees): ', JGL, (180._JPRB/3.14159_JPRB)*ASIN(ZMU(JGL))
    ENDDO
    WRITE(KULOUT,*) ' STARTING COMPUTATIONS '
  ENDIF

! 3. Control the polynomials orthogonality
!    -------------------------------------

  ALLOCATE(ZPNMT(KXMAX+3,ILEI3))
  IF (LHOOK) CALL DR_HOOK('RGRID:REDUCED_GRID>LOOP1',0,ZHOOK_HANDLE1)
  DO JM=1,INUMP
    IM=IMYMS(JM)
    CALL TRANS_PNM(KM=IM,PRPNM=ZPNMT,LDTRANSPOSE=.TRUE.,LDCHEAP=.TRUE.)
    IF (.NOT.LLINEAR .OR. IM <= KNMAX) THEN
      DO JOR=1,2
        IF (JOR == 1) THEN
          IOMAX=KXMAX
        ELSE
          IOMAX=KNMAX
        ENDIF
        DO JGL=1,IDGNH
          ZMAX(JGL)=-1._JPRB
        ENDDO
! FIXME:am:deactivated: not reproducible ! $OMP PARALLEL PRIVATE(JN1,JN2,JGL,ZORTH,ZUNITY,ZHOOK_HANDLE2)
        IF (LHOOK) CALL DR_HOOK('RGRID:REDUCED_GRID>LOOP2',0,ZHOOK_HANDLE2)
!!$OMP DO SCHEDULE(DYNAMIC,1)
        DO JN1=IM,KNMAX
          ZORTH(IM:IOMAX)=0._JPRB
          ZUNITY(IM:IOMAX)=0._JPRB
          DO JN2=IM,IOMAX
            IF (JN2==JN1) THEN
              ZUNITY(JN2)=1._JPRB
            ENDIF
          ENDDO
          DO JGL=IDGNH,1,-1
            IF (MOD(JN1-IM,2) == 0) THEN
              DO JN2=IM,IOMAX,2
                ZORTH(JN2)=ZORTH(JN2)+ZPNMT(KXMAX+2-JN1,JGL)*ZPNMT(KXMAX+2-JN2,JGL)*Z2GW(JGL)
                ZMAX(JGL)=MAX(ZMAX(JGL),ABS(ZORTH(JN2)-ZUNITY(JN2)))
              ENDDO
            ELSE
              DO JN2=IM+1,IOMAX,2
                ZORTH(JN2)=ZORTH(JN2)+ZPNMT(KXMAX+2-JN1,JGL)*ZPNMT(KXMAX+2-JN2,JGL)*Z2GW(JGL)
                ZMAX(JGL)=MAX(ZMAX(JGL),ABS(ZORTH(JN2)-ZUNITY(JN2)))
              ENDDO
            ENDIF
          ENDDO
        ENDDO
!!$OMP END DO
        IF (LHOOK) CALL DR_HOOK('RGRID:REDUCED_GRID>LOOP2',1,ZHOOK_HANDLE2)
! am:deactivated: not reproducible ! $OMP END PARALLEL
        DO JGL=1,IDGNH
          ZGR(JM,JGL,JOR)=ZMAX(JGL)
        ENDDO
        IF (KVERBOSE == 1) THEN
          IF(JOR == 2)THEN
            WRITE(KULERR,*) 'Filling ZGR JGL==1 IM=', IM, JM, ZGR(JM,1,JOR)
            WRITE(KULERR,*) 'Filling ZGR JGL==IDGNH/2 IM=', IM, JM, ZGR(JM,IDGNH/2,JOR)
          ENDIF
        ENDIF
      ENDDO
    ENDIF
  ENDDO
  IF (LHOOK) CALL DR_HOOK('RGRID:REDUCED_GRID>LOOP1',1,ZHOOK_HANDLE1)
  DEALLOCATE(ZPNMT)
  IF (LDWRITE) WRITE(IULZGR) ZGR(:,:,:)

ENDIF ! .NOT. LDREAD

CLMODE='FINAL'
CALL TRANS_END(CLMODE)


! 4. Compute the reduced grid
!    ------------------------

IF (KVERBOSE == 1) THEN
  IF (KMYPROC == 1) THEN
    WRITE(KULERR,*) ' start checking ffts '
  ENDIF
ENDIF

IF( LDREAD ) READ(IULZGR) ZGR(:,:,:)

LLODD=(KODD==1)
CALL SUFFTP(KDLON,IFFTP0,LLODD)
ALLOCATE(ILOENG(KDGLG))
ALLOCATE(ILOENS(KDGLG))
ILOENG(:)=0
ILOENS(:)=0

DO JGL=1,IDGNH
  DO JM=INUMP,1,-1
    IM=IMYMS(JM)
    IF (.NOT.LLINEAR .OR. IM <= KNMAX) THEN
      IF ((ZGR(JM,MIN(IDGNH,JGL+1),2)) >= PORTHO) THEN
        ILOENG(JGL)=2*MIN(IM+1,KNMAX)+1
        IF (KVERBOSE == 1) THEN
          WRITE(KULERR,*) 'decision1 ILOEN JGL IM', ILOENG(JGL), JGL, IM, ZGR(JM,MIN(IDGNH,JGL+1),2)
        ENDIF
        EXIT
      ENDIF
    ENDIF
  ENDDO
! fft compatibility :
  DO J=ILOENG(JGL),KDLON
    IF (IFFTP0(J) == 0)THEN
      ILOENS(JGL)=J
      EXIT
    ENDIF
  ENDDO
  DO JM=INUMP,1,-1
    IM=IMYMS(JM)
    IF (.NOT.LLINEAR .OR. IM <= KNMAX) THEN
      IF ((ZGR(JM,MIN(IDGNH,JGL+1),1)) >= PALIAS) THEN
        ILOENG(JGL)=MAX(ILOENG(JGL),3*MIN(IM+1,KNMAX)+1)
        IF (KVERBOSE == 1) THEN
          WRITE(KULERR,*) 'decision2 ILOEN JGL IM', ILOENG(JGL), JGL, IM, ZGR(JM,MIN(IDGNH,JGL+1),1)
        ENDIF
        EXIT
      ENDIF
    ENDIF
  ENDDO
! fft compatibility :
  DO J=ILOENG(JGL),KDLON
    IF (IFFTP0(J) == 0)THEN
      ILOENG(JGL)=J
      EXIT
    ENDIF
  ENDDO
ENDDO

DEALLOCATE(ZGR)
DEALLOCATE(IMYMS)

! 5. Gather contributions to task #1 and find the maximum per latitude
!    -----------------------------------------------------------------

IF (KMYPROC /= 1) THEN

  CALL MPL_SEND(ILOENG(1:IDGNH),KDEST=IPRCIDS(1),KTAG=ITAG1,CDSTRING='RGRID:')
  CALL MPL_SEND(ILOENS(1:IDGNH),KDEST=IPRCIDS(1),KTAG=ITAG2,CDSTRING='RGRID:')

ELSE

  ALLOCATE(IILOEN(IDGNH))
  DO JROC=1,IPRTRW-1
    CALL MPL_RECV(IILOEN,KTAG=ITAG1,KOUNT=IMSGLEN,CDSTRING='RGRID:')
    IF (IMSGLEN /= IDGNH) THEN
      CALL ABOR1('RGRID : RECEIVED MESSAGE LENGTH OF WRONG SIZE')
    ENDIF
    DO J=1,IDGNH
      ILOENG(J)=MAX(IILOEN(J),ILOENG(J))
    ENDDO
    CALL MPL_RECV(IILOEN,KTAG=ITAG2,KOUNT=IMSGLEN,CDSTRING='RGRID:')
    IF (IMSGLEN /= IDGNH) THEN
      CALL ABOR1('RGRID : RECEIVED MESSAGE LENGTH OF WRONG SIZE')
    ENDIF
    DO J=1,IDGNH
      ILOENS(J)=MAX(IILOEN(J),ILOENS(J))
    ENDDO
  ENDDO
  DEALLOCATE(IILOEN)

  DO J=1,KDGLG/2
    ILOENG(KDGLG+1-J)=ILOENG(J)
    ILOENS(KDGLG+1-J)=ILOENS(J)
  ENDDO
  IGPTOT=SUM(ILOENG)
  IGPTOTS=SUM(ILOENS)
  IF (KVERBOSE == 1) THEN
    DO JGL=1,IDGNH
      WRITE(KULOUT,FMT='(''   ILOEN('',I5,'')='',2I5)') JGL,ILOENG(JGL),ILOENS(JGL)
    ENDDO
  ENDIF
  WRITE(KULOUT,FMT='('' REDUCED GRID FOR NSMAX = '', I5, '' NDLON = '', I5, &
   & '' NDGL = '',I5,  '' IODD = '',I1)') KSMAXG,KDLON,KDGLG,KODD
  WRITE(KULOUT,FMT='('' ORTHOGONALITY = '',E16.9, '' ALIASING = '',E16.9)') PORTHO, PALIAS
  WRITE(KULOUT,FMT='('' IGPTOT = '',I8,'' IGPTOTS = '',I8)') IGPTOT, IGPTOTS
  WRITE(KULOUT,FMT='('' REDUCTION from a full grid = '',F4.1,''%'')') (1-REAL(IGPTOT,JPRB)/(KDGLG*KDLON))*100
  WRITE(KULOUT,FMT='('' KULNAM = '',I2)') KULNAM

  IF (KULNAM /= KULOUT .AND. KULNAM /= KULERR) OPEN(KULNAM,FORM="FORMATTED")
  WRITE(KULNAM,'('' &NAMRGRI'')')
  IF (KDLON.LE.99) THEN
    DO JGL=1,KDGLG
      WRITE(KULNAM,FMT='(''   NRGRI('',I2,'')='',I2,'','')') JGL,ILOENG(JGL)
    ENDDO
  ELSEIF (KDLON.LE.999) THEN
    DO JGL=1,KDGLG
      WRITE(KULNAM,FMT='(''   NRGRI('',I3,'')='',I3,'','')') JGL,ILOENG(JGL)
    ENDDO
  ELSEIF (KDLON.LE.9999) THEN 
    DO JGL=1,KDGLG
      WRITE(KULNAM,FMT='(''   NRGRI('',I4,'')='',I4,'','')') JGL,ILOENG(JGL)
    ENDDO
  ELSE
    DO JGL=1,KDGLG
      WRITE(KULNAM,FMT='(''   NRGRI('',I5,'')='',I5,'','')') JGL,ILOENG(JGL)
    ENDDO
  ENDIF
  WRITE(KULNAM,'('' /'')')
  IF (KULNAM /= KULOUT .AND. KULNAM /= KULERR) CLOSE(KULNAM)
ENDIF

DEALLOCATE(ILOENG)
DEALLOCATE(ILOENS)

IF (LDREAD .OR. LDWRITE) CLOSE(IULZGR)

!-----------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('RGRID:REDUCED_GRID',1,ZHOOK_HANDLE0)
END SUBROUTINE REDUCED_GRID

END PROGRAM RGRID
