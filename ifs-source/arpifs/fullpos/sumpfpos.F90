! (C) Copyright 1989- Meteo-France.

SUBROUTINE SUMPFPOS(CD,KFPRGPG,KNUMPROCFP,KPROMA_IN,KFPRGPNUM,KFPRGPL,KFPRGPLX,KFPRGPIND,KFPROMA,KFPBLOCS,KFPEND)

!**** *SUMPFPOS*  - SET UP MESSAGE PASSING FOR FULL-POS

!     PURPOSE.
!     --------
!        To initialize parameters of the distributed post-processing buffers

!**   INTERFACE.
!     ----------
!       *CALL* *SUMPFPOS*

!        EXPLICIT ARGUMENTS
!        ------------------
!         CD             : label for print
!         KFPRGPG        : global number of gridpoints
!         KNUMPROCFP     : MPI task rank for each point on the global grid
!         KPROMA_IN      : Suggested blocking factor
!         KFPRGPNUM      : Number of gridpoints on each task
!         KFPRGPL        : local number of gridpoints
!         KFPRGPLX       : maximum local number of gridpoints among all tasks
!         KFPRGPIND      : global address of each local gridpoint for each task
!         KFPROMA        : cache-blocking factor
!         KFPBLOCS       : Number of blocks
!         KFPEND         : actual number of gridpoint in each block
!         

!        IMPLICIT ARGUMENTS
!        ------------------
!        See modules above

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------
!       SUPROCFP

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS
!        Documentation about FULL-POS.

!     AUTHOR.
!     -------
!      RYAD EL KHATIB *METEO-FRANCE*
!      Original : 98-10-05 from sufpg

!     MODIFICATIONS.
!     --------------
!      R. El Khatib : 02-21-20 Fullpos B-level distribution + remove IO scheme
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      K. Yessad: 27-Feb-2007 Move old code in SUMPFPOS_DEP, adapt to DM-arrival geometry.
!      R. El Khatib  24-Jul-2012 LFPDISTRIB replaced by NFPDISTRIB
!      R. El Khatib & Tayfun Dalkilic 13-Sep-2012 NFPDISTRIB=2
!      R. El Khatib 27-Sep-2013 Boyd periodization in Fullpos-2 (differentiation
!      between interpolation grid and output grid)
!      T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!      R. El Khatib 27-Jul-2016 bugfix for the case NFPBW* > 0 + recode LWIDER_DOM + optimize for Boyd
!-----------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN   , ONLY : NULOUT

USE YOMMP0    ,ONLY : NPROC, MYPROC, NPRINTLEV
USE YOMFPC   , ONLY : LTRACEFP ,LALLOFP

IMPLICIT NONE

CHARACTER(LEN=*),   INTENT(IN)  :: CD
INTEGER(KIND=JPIM), INTENT(IN)  :: KFPRGPG
INTEGER(KIND=JPIM), INTENT(IN)  :: KNUMPROCFP(KFPRGPG)
INTEGER(KIND=JPIM), INTENT(IN)  :: KPROMA_IN
INTEGER(KIND=JPIM), ALLOCATABLE, INTENT(OUT) :: KFPRGPNUM(:)
INTEGER(KIND=JPIM), INTENT(OUT) :: KFPRGPL
INTEGER(KIND=JPIM), INTENT(OUT) :: KFPRGPLX
INTEGER(KIND=JPIM), ALLOCATABLE, INTENT(OUT) :: KFPRGPIND(:,:)
INTEGER(KIND=JPIM), INTENT(OUT) :: KFPROMA
INTEGER(KIND=JPIM), INTENT(OUT) :: KFPBLOCS
INTEGER(KIND=JPIM), ALLOCATABLE, INTENT(OUT) :: KFPEND(:)

INTEGER(KIND=JPIM) :: J, IROC, IPLP(NPROC)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUMPFPOS',0,ZHOOK_HANDLE)

! Number of output point each processor treats : 
ALLOCATE(KFPRGPNUM(NPROC))
IF (LALLOFP) THEN
  WRITE(NULOUT,'(''ARRAY NFPRGPNUM'',A,'' ALLOCATED '',8I8)') CD,SIZE(KFPRGPNUM),SHAPE(KFPRGPNUM)
ENDIF
DO J=1,NPROC
  KFPRGPNUM(J)=COUNT(KNUMPROCFP(:)==J)
ENDDO
KFPRGPL=KFPRGPNUM(MYPROC)
KFPRGPLX=MAXVAL(KFPRGPNUM)
IF (LTRACEFP) THEN
  WRITE(UNIT=NULOUT,FMT='('' NFPRGPL'',A,'' = '',I6,'' NFPRGPLX'',A,'' = '',I6)') CD, KFPRGPL, CD, KFPRGPLX
ENDIF
! Global index of each output point in each processor
ALLOCATE(KFPRGPIND(KFPRGPLX,NPROC))
IF (LALLOFP) THEN
  WRITE(NULOUT,'(''ARRAY NFPRGPIND'',A,'' ALLOCATED '',8I8)') CD,SIZE(KFPRGPIND),SHAPE(KFPRGPIND)
ENDIF
IF (NPROC == 1) THEN
  DO J=1,KFPRGPG
    KFPRGPIND(J,1)=J
  ENDDO
ELSE
  IF (KFPRGPLX > 0) THEN
    ! perhaps it is better to save numprocfp than nfprgpind ... REK
    KFPRGPIND(:,:)=HUGE(1_JPIM)
    IPLP(:)=0
    DO J=1,KFPRGPG
      IROC=KNUMPROCFP(J)
      IPLP(IROC)=IPLP(IROC)+1
      KFPRGPIND(IPLP(IROC),IROC)=J
    ENDDO
 ENDIF
ENDIF
IF (KPROMA_IN < 0) THEN
  ! Keep the absolute value
  KFPROMA=-KPROMA_IN
ELSEIF (KPROMA_IN > 0) THEN
  ! Ensure that the value is odd. We could make it more clever later (a divider of the number of blocks).
  KFPROMA=KPROMA_IN+1-MOD(KPROMA_IN,2)
ENDIF
KFPROMA=MAX(1,MIN(KFPROMA,KFPRGPL))
IF (KFPRGPL > 0) THEN
  KFPBLOCS=(KFPRGPL-1)/KFPROMA +1
ELSE
  KFPBLOCS=0
ENDIF
ALLOCATE(KFPEND(KFPBLOCS))
DO J=1,KFPBLOCS-1
  KFPEND(J)=KFPROMA
ENDDO
IF (KFPBLOCS > 0) THEN
  KFPEND(KFPBLOCS)=KFPRGPL-KFPROMA*(KFPBLOCS-1)
ENDIF
IF (LALLOFP) THEN
  WRITE(NULOUT,'(''ARRAY NFPEND'',A,'' ALLOCATED '',8I8)') CD,SIZE(KFPEND),SHAPE(KFPEND)
ENDIF
IF (LTRACEFP) THEN
  WRITE(NULOUT,'('' NFPROMA'',A,'' = '',I6,'' NFPBLOCS'',A,'' = '',I6)') CD, KFPROMA, CD, KFPBLOCS
  IF (NPRINTLEV > 1) THEN
    WRITE(NULOUT,'('' (IBLOC NFPEND'',A,'') '')') CD
    WRITE(NULOUT,'(8(1X,''('',I4,I6,'')''))') (J,KFPEND(J), J=1, KFPBLOCS)  
  ENDIF
ENDIF

IF (LHOOK) CALL DR_HOOK('SUMPFPOS',1,ZHOOK_HANDLE)
END SUBROUTINE SUMPFPOS
